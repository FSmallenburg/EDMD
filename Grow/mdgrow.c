#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <stdint.h>

#define MAXNEIGH 24

#include "mt19937ar.c"
#include "mdgrow.h"


//Size of array allocated for the event tree (currently overkill)
#define MAXEVENTS (N+12)

//Maximum number of cells in each direction
#define CEL 75

//Pi (if not already defined)
#ifndef M_PI
#define M_PI 3.1415926535897932
#endif

//Code for creating densely packed initial configurations. Particles start small in a random configuration, then grow to the desired size.
//The code is set up for binary mixtures. The parameters below control the number of particles, size ratio, composition, target packing fraction, etc.

//Number of particles:
#define N 100000


double targetpackfrac = 0.5; //Target packing fraction (if too high, simulation will not finish or crash)  
double composition = 1.0;     //Fraction of large particles
double sizeratio = 0.85;      //small diameter / large diameter (must be <= 1)
double growthspeed = 0.1;     //Factor determining growth speed (slower growth means higher packing fractions can be reached)
double thermostatinterval = 0.001;  //Time between applications of thermostat, which gets rid of excess heat generated while growing

int makesnapshots = 0;        //Whether to make snapshots during the run (yes = 1, no = 0)
double writeinterval = 1;     //Time between output to screen / data file
double snapshotinterval = 1;  //Time between snapshots (should be a multiple of writeinterval)


//Variables related to the event queueing system. These can affect efficiency.
//The system schedules only events in the current block of time with length "eventlisttime" into a sorted binary search tree. 
//The rest are scheduled in unordered linked lists associated with the "numeventlists" next blocks.
//"numeventlists" is roughly equal to maxscheduletime / eventlisttime
//Any events occurring even later are put into an overflow list
//After every time block with length "eventlisttime", the set of events in the next linear list is moved into the binary search try.
//All events in the overflow list are also rescheduled.

//After every "writeinterval", the code will output two listsizes to screen. 
//The first is the average number of events in the first that gets moved into the event tree after each block.
//The second is the average length of the overflow list.
//Ideally, we set maxscheduletime large enough that the average overflow list size is negligible (i.e. <10 events)
//Also, there is some optimum value for the number of events per block (scales approximately linearly with "eventlisttime").
//I seem to get good results with an eventlisttime chosen such that there are a few hundred events per block, and dependence is pretty weak (similar performance in the range of e.g. 5 to 500 events per block...)
#define MAXNUMEVENTLISTS (5*N)
double maxscheduletime = 5;
int numeventlists;
double eventlisttime = 2 / (float)N;

double maxdt = 1;

//Neighbor lists
double shellsize = 1.5;


//Internal variables
double time = 0;
double maxt = -1;

double reftime = 0;
int currentlist = 0;
const double never = 9999999999999999999999.9;


int listcounter1 = 0, listcounter2 = 0, mergecounter = 0;
int nbig;

event* eventlists[MAXNUMEVENTLISTS + 1]; //Last one is overflow list

particle particles[N];
particle* celllist[CEL][CEL][CEL];
event eventlist[MAXEVENTS];
event* root;
event* eventpool[MAXEVENTS];
int nempty = 0;
double density;
double xsize, ysize, zsize; //Box size
double hx, hy, hz; //Half box size
double cxsize, cysize, czsize; //Cell size
int    cx, cy, cz;  //Number of cells
double dvtot = 0;   //Momentum transfer (for calculating pressure)
unsigned int colcounter = 0; //Collision counter (will probably overflow in a long run...)
int stop = 0;


int usethermostat = 1; //Whether to use a thermostat



int main()
{
    init();
    if (stop) return 1;				//Stop is set to 1 whenever the simulation should stop
    printf("Starting\n");
    while (!stop)
    {
        step();
    }
    printstuff();
    return 0;
}

/**************************************************
**                 PRINTSTUFF
** Some data at the end of the simulation
**************************************************/
void printstuff()
{
    int i;
    particle* p;
    double v2tot = 0;
    double vfilled = 0;

    for (i = 0; i < N; i++)
    {
        p = particles + i;
        v2tot += p->mass * (p->vx * p->vx + p->vy * p->vy + p->vz * p->vz);
        vfilled += p->r * p->r * p->r * 8;
    }
    vfilled *= M_PI / 6.0;
    printf("Average kinetic energy: %lf\n", 0.5 * v2tot / N);
    double volume = xsize * ysize * zsize;
    double dens = N / volume;
    double press = -dvtot / (3.0 * volume * time);
    double pressid = dens;
    double presstot = press + pressid;
    printf("Total time simulated  : %lf\n", time);
    //  printf ("Density               : %lf\n", (double) N / volume);
    printf("Packing fraction      : %lf\n", vfilled / volume);
    printf("Measured pressure     : %lf + %lf = %lf\n", press, pressid, presstot);

}


/**************************************************
**                    INIT
**************************************************/
void init()
{
    int i;
    unsigned long seed = 1;
    //   FILE *fp=fopen("/dev/urandom","r");
    //   int tmp = fread(&seed,1,sizeof(unsigned long),fp);
    //   if (tmp != sizeof(unsigned long)) printf ("error with seed\n");
    //   fclose(fp);
    printf("Seed: %u\n", (int)seed);
    init_genrand(seed);
    initeventpool();

    for (i = 0; i < N; i++)
    {
        particle* p = particles + i;
        p->number = i;
        p->boxestraveledx = 0;
        p->boxestraveledy = 0;
        p->boxestraveledz = 0;
        p->nneigh = 0;
        p->counter = 0;
        p->t = 0;
    }


    randomparticles();
    randommovement();
    hx = 0.5 * xsize; hy = 0.5 * ysize; hz = 0.5 * zsize;	//Values used for periodic boundary conditions



    for (i = 0; i < N; i++)
    {
        particle* p = particles + i;
        p->xn = p->x;
        p->yn = p->y;
        p->zn = p->z;
    }

    for (i = 0; i < N; i++)
    {
        makeneighborlist(particles + i, 1);
    }
    printf("Done adding collisions: %d events\n", MAXEVENTS - nempty);


    if (usethermostat)
    {
        thermostat(NULL);
        printf("Started thermostat: %d events\n", MAXEVENTS - nempty);
    }

}


/******************************************************
**               MYGETLINE
** Reads a single line, skipping over lines
** commented out with #
******************************************************/
int mygetline(char* str, FILE* f)
{
    int comment = 1;
    while (comment)
    {
        if (!fgets(str, 255, f)) return -1;
        if (str[0] != '#') comment = 0;
    }
    return 0;
}




/**************************************************
**                    RANDOMPARTICLES
** Positions particles randomly in the box
**************************************************/
void randomparticles()
{
    double x = composition;
    double alpha = sizeratio;
    double vol = N * M_PI / 6 * (x + (1 - x) * alpha * alpha * alpha) / targetpackfrac;

    printf("Volume: %lf\n", vol);

    xsize = cbrt(vol);
    ysize = xsize;
    zsize = ysize;
    initcelllist();
    int i;
    particle* p;
    for (i = 0; i < N; i++)				//First put particles at zero
    {
        particles[i].x = 0; particles[i].y = 0; particles[i].z = 0;
    }
    for (i = 0; i < N; i++)
    {
        p = &(particles[i]);
        p->rtarget = 1;
        if (i >= x * N - 0.00000001) p->rtarget = alpha;
        p->rtarget *= 0.5;
        p->r = 0.5 * p->rtarget;
        p->mass = 1;
        p->type = 0;
        if (p->rtarget < 0.5) p->type = 1;
        p->number = i;
        p->vr = growthspeed * p->r;
        double mt = (p->rtarget - p->r) / p->vr;
        if (mt < 0) printf("Particles should not exceed size 1!\n");
        if (maxt < 0 || maxt > mt) maxt = mt;
        do
        {
            p->x = genrand_real2() * xsize;			//Random location and speed
            p->y = genrand_real2() * ysize;
            p->z = genrand_real2() * zsize;
            p->cellx = p->x / cxsize;				//Find particle's cell
            p->celly = p->y / cysize;
            p->cellz = p->z / czsize;
        } while (overlaplist(p, 0));
        double sqm = 1.0 / sqrt(p->mass);
        p->xn = p->x; p->yn = p->y; p->zn = p->z;
        p->vx = (genrand_real2() - 0.5) / sqm;
        p->vy = (genrand_real2() - 0.5) / sqm;
        p->vz = (genrand_real2() - 0.5) / sqm;
        p->t = 0;						//r and v known at t=0
        p->next = celllist[p->cellx][p->celly][p->cellz];	//Add particle to celllist
        if (p->next) p->next->prev = p;			//Link up list
        celllist[p->cellx][p->celly][p->cellz] = p;
        p->prev = NULL;
    }
}


/**************************************************
**                RANDOMMOVEMENT
**************************************************/
void randommovement()
{
    particle* p;
    double v2tot = 0, vxtot = 0, vytot = 0, vztot = 0;
    double mtot = 0;
    int i;

    for (i = 0; i < N; i++)
    {
        p = particles + i;
        double imsq = 1.0 / sqrt(p->mass);

        p->vx = imsq * random_gaussian();
        p->vy = imsq * random_gaussian();
        p->vz = imsq * random_gaussian();
        vxtot += p->mass * p->vx;					//Keep track of total v
        vytot += p->mass * p->vy;
        vztot += p->mass * p->vz;
        mtot += p->mass;
    }


    vxtot /= mtot; vytot /= mtot; vztot /= mtot;
    for (i = 0; i < N; i++)
    {
        p = &(particles[i]);
        p->vx -= vxtot;					//Make sure v_cm = 0
        p->vy -= vytot;
        p->vz -= vztot;
        v2tot += p->mass * (p->vx * p->vx + p->vy * p->vy + p->vz * p->vz);
    }
    double fac = sqrt(3.0 / (v2tot / N));
    v2tot = 0;
    vxtot = vytot = vztot = 0;
    for (i = 0; i < N; i++)
    {
        p = &(particles[i]);
        p->vx *= fac;					//Fix energy
        p->vy *= fac;
        p->vz *= fac;
        v2tot += p->mass * (p->vx * p->vx + p->vy * p->vy + p->vz * p->vz);
        vxtot += p->mass * p->vx;					//Keep track of total v
        vytot += p->mass * p->vy;
        vztot += p->mass * p->vz;

    }
    printf("average v2: %lf (%lf, %lf, %lf)\n", v2tot / N, vxtot / N, vytot / N, vztot / N);
}

/**************************************************
**                UPDATE
**************************************************/
void update(particle* p1)
{
    double dt = time - p1->t;
    p1->t = time;
    p1->x += dt * p1->vx;
    p1->y += dt * p1->vy;
    p1->z += dt * p1->vz;
    p1->r += dt * p1->vr;

}



/**************************************************
**                 INITCELLLIST
**************************************************/
void initcelllist()
{
    int i, j, k;
    cx = (int)(xsize - 0.0001) / shellsize;				//Set number of cells
    cy = (int)(ysize - 0.0001) / shellsize;
    cz = (int)(zsize - 0.0001) / shellsize;
    printf("Cells: %d, %d, %d\n", cx, cy, cz);
    if (cx >= CEL || cy >= CEL || cz >= CEL)
    {
        printf("Too many cells!\n");
        stop = 1; return;
    }
    cxsize = xsize / cx;						//Set cell size
    cysize = ysize / cy;
    czsize = zsize / cz;
    for (i = 0; i < CEL; i++)					//Clear celllist
        for (j = 0; j < CEL; j++)
            for (k = 0; k < CEL; k++)
            {
                celllist[i][j][k] = NULL;
            }

}


/**************************************************
**               REMOVEFROMCELLLIST
**************************************************/
void removefromcelllist(particle* p1)
{
    if (p1->prev) p1->prev->next = p1->next;    //Remove particle from celllist
    else          celllist[p1->cellx][p1->celly][p1->cellz] = p1->next;
    if (p1->next) p1->next->prev = p1->prev;
}

/**************************************************
**                    ADDTOCELLLIST
**************************************************/
void addtocelllist(particle* p)
{
    p->cellx = p->x / cxsize;				//Find particle's cell
    p->celly = p->y / cysize;
    p->cellz = p->z / czsize;
    p->next = celllist[p->cellx][p->celly][p->cellz];	//Add particle to celllist
    if (p->next) p->next->prev = p;			//Link up list
    celllist[p->cellx][p->celly][p->cellz] = p;
    p->prev = NULL;
    p->edge = (p->cellx == 0 || p->celly == 0 || p->cellz == 0 || p->cellx == cx - 1 || p->celly == cy - 1 || p->cellz == cz - 1);

}

/**************************************************
**                     STEP
**************************************************/
void step()
{
    event* ev;
    ev = root->child2;
    if (ev == NULL)
    {
        addnexteventlist();
        ev = root->child2;
    }

    while (ev->child1) ev = ev->child1;		//Find first event


    //if (ev->type == 0)
    //{
      //printf("Time: %lf, ev: %d, part: %d, %d\n", ev->time, ev->type, ev->p1 ? ev->p1->number : -1, ev->p2 ? ev->p2->number : -1);
    //}

   //if ((ev->p1 && ev->p1->number == 591) && (ev->p2 && ev->p2->number == 701) || 
   //     (ev->p1 && ev->p1->number == -1) || (ev->p2 && ev->p2->number == -1))
   // {
   //      printf ("Time: %.10lf, ev: %d, part: %d, %d (%lf, %d)\n", ev->time, ev->type, ev->p1?ev->p1->number:-1, ev->p2?ev->p2->number:-1, particles[621].z, particles[374].cellz);
   // }

   // if (ev->time < time)
   // {
   //     printf("Time: %lf, ev: %d, part: %d, %d\n", ev->time, ev->type, ev->p1 ? ev->p1->number : -1, ev->p2 ? ev->p2->number : -1);
   //     exit(3);
   // }

   // if (ev->time < time)
   // {
   //     printf(" time error: %lf, %lf \n", time, ev->time);
   //     printf("Time: %.10lf, ev: %d, part: %d, %d (%lf, %d)\n", ev->time, ev->type, ev->p1 ? ev->p1->number : -1, ev->p2 ? ev->p2->number : -1, particles[621].z, particles[374].cellz);

   //     exit(3);
   // }

    if (ev->time > maxt)
    {
        time = maxt;
        write();
        writelast();
        printf("Time is up!\n");
        stop = 1;
    }

    if (ev->type == 100)
    {
        time = ev->time;
        removeevent(ev);
        write();
    }
    else if (ev->type == 200)
    {
        thermostat(ev);
    }
    else if (ev->type == 8)
    {
        time = ev->time;
        removeevent(ev);
        makeneighborlist(ev->p1, 0);
    }
    else
    {
        collision(ev);
    }
}



/**************************************************
**                MAKENEIGHBORLIST
**************************************************/
void makeneighborlist(particle* p1, int firsttime)
{

    int cdx, cdy, cdz, cellx, celly, cellz;
    particle* p2;
    double dx, dy, dz, r2, rm;

    update(p1);

    if (p1->x >= xsize) { p1->x -= xsize; p1->boxestraveledx++; }
    else if (p1->x < 0) { p1->x += xsize; p1->boxestraveledx--; }
    if (p1->y >= ysize) { p1->y -= ysize; p1->boxestraveledy++; }
    else if (p1->y < 0) { p1->y += ysize; p1->boxestraveledy--; }
    if (p1->z >= zsize) { p1->z -= zsize; p1->boxestraveledz++; }
    else if (p1->z < 0) { p1->z += zsize; p1->boxestraveledz--; }
    p1->xn = p1->x;
    p1->yn = p1->y;
    p1->zn = p1->z;



    removefromcelllist(p1);
    addtocelllist(p1);


    int i, j;
    for (i = 0; i < p1->nneigh; i++)
    {
        p2 = p1->neighbors[i];
        for (j = 0; j < p2->nneigh; j++)
        {
            if (p2->neighbors[j] == p1)
            {
                p2->nneigh--;
                p2->neighbors[j] = p2->neighbors[p2->nneigh];
                break;
            }
        }
    }


    cellx = p1->cellx + cx;
    celly = p1->celly + cy;
    cellz = p1->cellz + cz;

    p1->nneigh = 0;

    for (cdx = cellx - 1; cdx < cellx + 2; cdx++)
        for (cdy = celly - 1; cdy < celly + 2; cdy++)
            for (cdz = cellz - 1; cdz < cellz + 2; cdz++)
            {
                p2 = celllist[cdx % cx][cdy % cy][cdz % cz];
                while (p2)
                {
                    if (p2 != p1)
                    {
                        update(p2);
                        dx = p1->xn - p2->xn;
                        dy = p1->yn - p2->yn;
                        dz = p1->zn - p2->zn;
                        if (p1->edge)
                        {
                            if (dx > hx) dx -= xsize; else if (dx < -hx) dx += xsize;  //periodic boundaries
                            if (dy > hy) dy -= ysize; else if (dy < -hy) dy += ysize;
                            if (dz > hz) dz -= zsize; else if (dz < -hz) dz += zsize;
                        }
                        r2 = dx * dx + dy * dy + dz * dz;
                        rm = (p1->r + p1->vr*maxdt + p2->rtarget + p2->vr*maxdt) * shellsize;
                        if (r2 < rm * rm)
                        {
                            p1->neighbors[p1->nneigh++] = p2;
                            p2->neighbors[p2->nneigh++] = p1;
                            if (p1->nneigh >= MAXNEIGH || p2->nneigh >= MAXNEIGH)
                            {
                                printf ("Too many neighbors; increase MAXNEIGH\n");
                                exit(3);
                            }
			}
                    }
                    p2 = p2->next;
                }
            }

    //printf("Neighbors: %d, %d\n", p1->number, p1->nneigh);
    //if (p1->nneigh > MAXNEIGH) exit(3);

    findcollisions(p1);


}


/**************************************************
**                FINDNEIGHBORLISTUPDATE
** Assumes p1 is up to date
** Note that the particle is always in the same
** box as its neighborlist position (p->xn)
**************************************************/
double findneighborlistupdate(particle* p1)
{
    double dx = p1->x - p1->xn;
    double dy = p1->y - p1->yn;
    double dz = p1->z - p1->zn;

    double dvx = p1->vx, dvy = p1->vy, dvz = p1->vz;

    double b = dx * dvx + dy * dvy + dz * dvz;                  //dr.dv

    double dv2 = dvx * dvx + dvy * dvy + dvz * dvz;
    double dr2 = dx * dx + dy * dy + dz * dz;
    double md = (shellsize - 1) * p1->rtarget;

    double disc = b * b - dv2 * (dr2 - md * md);
    double t = (-b + sqrt(disc)) / dv2;
    //printf("Predicting nlistupdate %d, %lf\n", p1->number, t);
    return t;
}

/**************************************************
**                FINDCOLLISION
** Detect the next collision for two particles
** Note that p1 is always up to date in
** findcollision
**************************************************/
double findcollision(particle* p1, particle* p2, double tmin)
{
    double dt2 = time - p2->t;
    double dx = p1->x - p2->x - dt2 * p2->vx;    //relative distance at current time
    double dy = p1->y - p2->y - dt2 * p2->vy;
    double dz = p1->z - p2->z - dt2 * p2->vz;
    if (dx > hx) dx -= xsize; else if (dx < -hx) dx += xsize;  //periodic boundaries
    if (dy > hy) dy -= ysize; else if (dy < -hy) dy += ysize;
    if (dz > hz) dz -= zsize; else if (dz < -hz) dz += zsize;
    double dvx = p1->vx - p2->vx;                               //relative velocity
    double dvy = p1->vy - p2->vy;
    double dvz = p1->vz - p2->vz;
    double dvr = p1->vr + p2->vr;
    double md = p1->r + p2->r + dt2 * p2->vr;


    double b = dx * dvx + dy * dvy + dz * dvz - dvr * md;     //dr.dv
    if (b > 0) return never;
    double dv2 = dvx * dvx + dvy * dvy + dvz * dvz;
    double dr2 = dx * dx + dy * dy + dz * dz;

    double twoa = dv2 - dvr * dvr;

    double disc = b * b - twoa * (dr2 - md * md);
    if (disc < 0) return never;
    double t = (-b - sqrt(disc)) / twoa;
    return t;

}




/**************************************************
**                FINDCOLLISIONS
** Find all collisions for particle p1.
** The particle 'not' isn't checked.
**************************************************/
void findcollisions(particle* p1)    //All collisions of particle p1
{
    int i;
    double t;
    double tmin = findneighborlistupdate(p1);
    if (tmin > maxdt) tmin = maxdt;
    int type = 8;
    particle* partner = p1;
    particle* p2;
    for (i = 0; i < p1->nneigh; i++)
    {
        p2 = p1->neighbors[i];
        t = findcollision(p1, p2, tmin);
        if (t < tmin)
        {
            tmin = t;
            partner = p2;
            type = 0;
        }
    }
    event* ev = createevent(tmin + time, p1, partner, type);
    p1->firstcollision = ev;
    ev->counter2 = partner->counter;
}



/**************************************************
**                FINDALLCOLLISION
** All collisions of all particle pairs
**************************************************/
void findallcollisions()       //All collisions of all particle pairs
{
    int i, j;


    for (i = 0; i < N; i++)
    {
        particle* p1 = particles + i;
        particle* partner = p1;
        double tmin = findneighborlistupdate(p1);
        int type = 8;
        for (j = 0; j < p1->nneigh; j++)
        {
            particle* p2 = p1->neighbors[j];
            if (p2 > p1)
            {
                double t = findcollision(p1, p2, tmin);
                if (t < tmin)
                {
                    tmin = t;
                    partner = p2;
                    type = 0;
                }
            }
        }
        if (partner)
        {
            event* ev = createevent(tmin, p1, partner, type);
            p1->firstcollision = ev;
            ev->counter2 = partner->counter;
        }

    }
}







/**************************************************
**                  COLLISION
** Process a single collision event
**************************************************/
void collision(event* ev)
{
    time = ev->time;
    particle* p1 = ev->p1;
    particle* p2 = ev->p2;
    update(p1);
    removeevent(ev);
    if (ev->counter2 != p2->counter)
    {
        findcollisions(p1);
        return;
    }

    //     if (ev->counter1 != p1->counter)
    //     {
    //         printf("Huh?\n");
    //     }
    update(p2);
    p1->counter++;
    p2->counter++;



    double m1 = p1->mass, r1 = p1->r;
    double m2 = p2->mass, r2 = p2->r;

    double r = r1 + r2;
    double rinv = 1.0 / r;
    double dx = (p1->x - p2->x);			//Normalized distance vector
    double dy = (p1->y - p2->y);
    double dz = (p1->z - p2->z);
    if (p1->edge)
    {
        if (dx > hx) dx -= xsize; else if (dx < -hx) dx += xsize;  //periodic boundaries
        if (dy > hy) dy -= ysize; else if (dy < -hy) dy += ysize;
        if (dz > hz) dz -= zsize; else if (dz < -hz) dz += zsize;
    }
    dx *= rinv;  dy *= rinv;  dz *= rinv;

    double dvx = p1->vx - p2->vx;                               //relative velocity
    double dvy = p1->vy - p2->vy;
    double dvz = p1->vz - p2->vz;
    double dvr = p1->vr + p2->vr;

    double b = dx * dvx + dy * dvy + dz * dvz - dvr;                  //dr.dv
    b *= 1.0 / (m1 + m2);
    double dv1 = 2 * b * m2;
    double dv2 = 2 * b * m1;
    dvtot += 2 * m1 * m2 * b;

    p1->vx -= dv1 * dx;         //Change velocities after collision
    p1->vy -= dv1 * dy;         //delta v = (-) dx2.dv2
    p1->vz -= dv1 * dz;
    p2->vx += dv2 * dx;
    p2->vy += dv2 * dy;
    p2->vz += dv2 * dz;

    colcounter++;

    if (p2->firstcollision && p2->firstcollision != ev)
    {
        removeevent(p2->firstcollision);
    }

    findcollisions(p1);
    findcollisions(p2);





}






/**************************************************
**                 INITEVENTPOOL
** Creates two first events, and sets up
**************************************************/
void initeventpool()
{
    numeventlists = ceil(maxscheduletime / eventlisttime);
    maxscheduletime = numeventlists * eventlisttime;
    if (numeventlists > MAXNUMEVENTLISTS)
    {
        printf("Number of event lists too large: increase MAXNUMEVENTLISTS to at least %d\n", numeventlists);
        exit(3);
    }
    printf("number of lists: %d\n", numeventlists);


    int i;
    event* e;
    for (i = 0; i < MAXEVENTS; i++)			//Clear eventpool
    {
        e = &(eventlist[MAXEVENTS - i - 1]);			//Fill in backwards, so the first few events are 'special'
        eventpool[i] = e;					//This includes root, the write events, in that order
        eventpool[i]->child1 = NULL;			//...  Not really used for now, but it might be useful at some point
        eventpool[i]->child2 = NULL;			//Clear children
        nempty++;						//All events empty so far
    }
    root = eventpool[--nempty];				//Create root event
    root->time = -99999999999.99;				//Root event is empty, but makes sure every other event has a parent
    root->type = 200;					//This makes sure we don't have to keep checking this when adding/removing events
    root->parent = NULL;
    event* writeevent = eventpool[--nempty];		//Pick first unused event
    writeevent->time = 0;
    writeevent->type = 100;
    root->child2 = writeevent;
    writeevent->parent = root;
    printf("Event tree initialized: %d events\n", MAXEVENTS - nempty);
}

/**************************************************
**                  ADDEVENTTOTREE
**************************************************/
void addeventtotree(event* newevent)
{
    double time = newevent->time;
    event* loc = root;
    int busy = 1;
    while (busy)						//Find location to add event into tree (loc)
    {
        if (time < loc->time)				//Go left
        {
            if (loc->child1) loc = loc->child1;
            else
            {
                loc->child1 = newevent;
                busy = 0;
            }
        }
        else						//Go right
        {
            if (loc->child2) loc = loc->child2;
            else
            {
                loc->child2 = newevent;
                busy = 0;
            }
        }
    }
    newevent->parent = loc;

}

/**************************************************
**                  ADDEVENT
**************************************************/
void addevent(event* newevent)
{
    double dt = newevent->time - reftime;

    if (dt < eventlisttime) //Put it in the tree
    {
        newevent->queue = currentlist;
        addeventtotree(newevent);
    }
    else
    {
        int list_id;
        if (dt >= numeventlists * eventlisttime) list_id = numeventlists;    //This also handles int overflow when calculating list_id
        else
        {
            list_id = currentlist + dt / eventlisttime;
            if (list_id >= numeventlists)
            {
                list_id -= numeventlists;
            }
        }

        newevent->queue = list_id;
        newevent->nextq = eventlists[list_id]; //Add to linear list
        newevent->prevq = NULL;
        if (newevent->nextq) newevent->nextq->prevq = newevent;
        eventlists[list_id] = newevent;
    }
}
/**************************************************
**                  CREATEEVENT
**************************************************/
event* createevent(double time, particle* p1, particle* p2, int type)
{
    event* newevent = eventpool[--nempty];		//Pick first unused event
    newevent->time = time;
    newevent->p1 = p1;
    newevent->type = type;
    newevent->p2 = p2;
    addevent(newevent);
    return newevent;
}

/**************************************************
**                     ADDNEXTEVENTLIST
**************************************************/
void addnexteventlist()
{
    do
    {
        currentlist++;
        if (currentlist == numeventlists) currentlist = 0;
        reftime += eventlisttime;
    } while (eventlists[currentlist] == NULL);

    //   printf("Currentlist is now %d (%lf)\n", currentlist, reftime);

    event* ev = eventlists[currentlist];
    event* nextev;
    while (ev)
    {
        nextev = ev->nextq;
        //         if (ev->type != 0 || ev->counter2 == ev->p2->counter) 
        addeventtotree(ev);
        ev = nextev;
        listcounter1++;
    }
    eventlists[currentlist] = NULL;
    ev = eventlists[numeventlists];//Overflow queue
    eventlists[numeventlists] = NULL;
    while (ev)
    {
        nextev = ev->nextq;
        addevent(ev);
        ev = nextev;
        listcounter2++;
    }
    mergecounter++;
}

/**************************************************
**                  REMOVEEVENT
**************************************************/
void removeevent(event* oldevent)
{

    //event* ev = oldevent;
    //if (ev->type != 8) printf("Removing event: %lf, ev: %d, part: %d, %d\n", ev->time, ev->type, ev->p1 ? ev->p1->number : -1, ev->p2 ? ev->p2->number : -1);

    if (oldevent->queue != currentlist)
    {
        if (oldevent->nextq) oldevent->nextq->prevq = oldevent->prevq;
        if (oldevent->prevq) oldevent->prevq->nextq = oldevent->nextq;
        else
        {
            eventlists[oldevent->queue] = oldevent->nextq;
        }
        eventpool[nempty++] = oldevent;     //Put the removed event back in the event pool.
        return;
    }

    event* parent = oldevent->parent;
    event* node;					//This node will be attached to parent in the end


    if (oldevent->child1 == NULL)			//Only one child: easy to delete
    {
        node = oldevent->child2;			//Child2 is attached to parent
        if (node)
        {
            node->parent = parent;
            oldevent->child2 = NULL;			//Clear child, so createevent doesn't have to do it
        }
    }
    else if (oldevent->child2 == NULL)		//Only one child again
    {
        node = oldevent->child1;			//Child1 is attached to parent
        node->parent = parent;
        oldevent->child1 = NULL;
    }
    else		  //Node to delete has 2 children
    {               //In this case: a) Find first node after oldevent     (This node will have no child1)
                    //              b) Remove this node from the tree     (Attach node->child2 to node->parent)
                    //              c) Put this node in place of oldevent (Oldevent's children are adopted by node)
        node = oldevent->child2;
        while (node->child1) node = node->child1;	//Find first node of right tree of descendants of oldevent
        event* pnode = node->parent;
        if (pnode != oldevent)			//node is not a child of oldevent
        {						//Both of oldevent's children should be connected to node
            pnode->child1 = node->child2;		//Remove node from right tree
            if (node->child2) node->child2->parent = pnode;
            oldevent->child1->parent = node;
            node->child1 = oldevent->child1;
            oldevent->child2->parent = node;
            node->child2 = oldevent->child2;
        }
        else					//This means node == oldevent->child2
        {						//Only child1 has to be attached to node
            oldevent->child1->parent = node;
            node->child1 = oldevent->child1;
        }
        node->parent = parent;
        oldevent->child1 = NULL;
        oldevent->child2 = NULL;
    }
    if (parent->child1 == oldevent) parent->child1 = node;
    else                            parent->child2 = node;
    eventpool[nempty++] = oldevent;     //Put the removed event back in the event pool.
}

/**************************************************
**                  SHOWTREE
** Gives a rough view of the event tree.
** Not so useful except for very small trees
**************************************************/
void showtree()
{
    shownode(root);
}

void shownode(event* ev)
{
    int c1 = 0, c2 = 0, p = 0;
    if (ev->child1) c1 = ev->child1->type;
    if (ev->child2) c2 = ev->child2->type;
    if (ev->parent) p = ev->parent->type;
    printf("%3d => %3d => %d (p: %d, %d)\n", p, ev->type, c1, ev->p1->number, ev->p2->number);
    printf("           => %d \n", c2);

    if (ev->child1) shownode(ev->child1);
    if (ev->child2) shownode(ev->child2);
}

/**************************************************
**                  CHECKTREE
** Checks the tree for possible errors.
**  1 ) node's parent doesn't point to node
**  1b) node->t > parent->t
**  1c) node->t < parent->t
**  2 ) A non-root node lacks a parent
**  3 ) node's child1 doesn't point to parent
**  3b) node->t < child1->t
**  4 ) node's child2 doesn't point to parent
**  4b) node->t > child2->t
** Also checks if all events are in the tree
**************************************************/
void checktree()
{
    return;
    int t = checknode(root);
    if (t != MAXEVENTS - nempty) printf("Error: %d, %d\n", t, MAXEVENTS - nempty);
}

int checknode(event* node)
{
    static int count = 0;
    if (node->parent)
    {
        if (node->parent->child1 != node && node->parent->child2 != node) printf("Error 1\n");
        if (node->parent->child1 == node && node->time > node->parent->time) printf("Error 1b\n");
        if (node->parent->child2 == node && node->time < node->parent->time) printf("Error 1c\n");
    }
    else
    {
        if (root != node) printf("Error 2\n");
        count = 0;
    }
    if (node->child1)
    {
        checknode(node->child1);
        if (node->child1->parent != node) printf("Error 3\n");
        if (node->child1->time > node->time) printf("Error 3b\n");
    }
    if (node->child2)
    {
        checknode(node->child2);
        if (node->child2->parent != node) printf("Error 4\n");
        if (node->child2->time < node->time) printf("Error 4b\n");
    }
    count++;
    return count;
}

/**************************************************
**                    OVERLAP
** Checks for overlaps
** Should write one that uses the cell list...
**************************************************/
int overlap(particle* part)
{
    particle* p;
    double dx, dy, dz, r2, rm;
    int i;
    double dl = pow(10.0, -10);

    for (i = 0; i < N; i++)
    {
        if (i == part->number) continue;
        p = &(particles[i]);
        dx = part->x - p->x;
        dy = part->y - p->y;
        dz = part->z - p->z;
        if (dx > 0.5 * xsize) dx -= xsize; else if (dx < -0.5 * xsize) dx += xsize;  //periodic boundaries
        if (dy > 0.5 * ysize) dy -= ysize; else if (dy < -0.5 * ysize) dy += ysize;
        if (dz > 0.5 * zsize) dz -= zsize; else if (dz < -0.5 * zsize) dz += zsize;
        r2 = dx * dx + dy * dy + dz * dz;
        rm = p->r + part->r;
        if (r2 < rm * rm - dl)
        {
            printf("Overlap: %lf, %d, %d\n", r2, part->number, p->number);
            return 1;
        }
    }
    return 0;
}

/**************************************************
**                    OVERLAPLIST
** Checks for overlaps
** if error is one, allow a small margin of error
**
**************************************************/
int overlaplist(particle* part, int error)
{
    int cdx, cdy, cdz, cellx, celly, cellz, num;
    particle* p;
    double dx, dy, dz, r2, rm;
    double dl = error * pow(10.0, -10);

    cellx = part->cellx + cx;
    celly = part->celly + cy;
    cellz = part->cellz + cz;
    num = part->number;

    for (cdx = cellx - 1; cdx < cellx + 2; cdx++)
        for (cdy = celly - 1; cdy < celly + 2; cdy++)
            for (cdz = cellz - 1; cdz < cellz + 2; cdz++)
            {
                p = celllist[cdx % cx][cdy % cy][cdz % cz];
                while (p)
                {
                    if (p->number != num)
                    {
                        dx = part->x - p->x;
                        dy = part->y - p->y;
                        dz = part->z - p->z;
                        if (dx > 0.5 * xsize) dx -= xsize; else if (dx < -0.5 * xsize) dx += xsize;  //periodic boundaries
                        if (dy > 0.5 * ysize) dy -= ysize; else if (dy < -0.5 * ysize) dy += ysize;
                        if (dz > 0.5 * zsize) dz -= zsize; else if (dz < -0.5 * zsize) dz += zsize;
                        r2 = dx * dx + dy * dy + dz * dz;
                        rm = p->r + part->r;
                        if (r2 < rm * rm - dl)
                        {
                            //            printf ("Overlap: %lf, %d, %d\n", r2, part->number, p->number);
                            return 1;
                        }
                    }
                    p = p->next;
                }
            }
    return 0;
}




/**************************************************
**                    WRITE
** Writes a movie
**************************************************/
void write()
{
    static int counter = 0;
    static int first = 1;
    static double lastsnapshottime = -999999999.9;
    int i;
    particle* p;
    FILE* file;

    double en = 0;
    for (i = 0; i < N; i++)
    {
        p = particles + i;
        update(p);
        en += p->mass * (p->vx * p->vx + p->vy * p->vy + p->vz * p->vz);
    }
    double temperature = 0.5 * en / (float)N / 1.5;


    //     checktree();
    //     checkcells();
    double volume = xsize * ysize * zsize;
    double dens = N / volume;
    double timeint = time;
    if (timeint < 0) timeint = time;
    double press = -dvtot / (3.0 * volume * timeint);
    double pressid = dens;
    double presstot = colcounter ? press + pressid : 0;

    double listsize1 = (double)listcounter1 / mergecounter;
    double listsize2 = (double)listcounter2 / mergecounter;
    if (mergecounter == 0) listsize1 = listsize2 = 0.0;
    listcounter1 = listcounter2 = mergecounter = 0;

    printf("Simtime: %lf, Collisions: %u, Press: %lf (%lf), T: %lf, Listsizes: (%lf, %lf)\n", time, colcounter, presstot, presstot / dens, temperature, listsize1, listsize2);
    char filename[200];
    if (makesnapshots && time - lastsnapshottime > snapshotinterval - 0.001)
    {
        sprintf(filename, "mov.n%d.v%.4lf.sph", N, xsize * ysize * zsize);
        if (first) { first = 0; file = fopen(filename, "w"); }
        else                     file = fopen(filename, "a");
        fprintf(file, "%d\n%.12lf %.12lf %.12lf\n", (int)N, xsize, ysize, zsize);
        for (i = 0; i < N; i++)
        {
            p = &(particles[i]);

            fprintf(file, "%c %.12lf  %.12lf  %.12lf  %lf\n", 'a' + p->type, p->x + xsize * p->boxestraveledx, p->y + ysize * p->boxestraveledy, p->z + zsize * p->boxestraveledz, p->r);
        }
        fclose(file);
        lastsnapshottime = time;
    }
    if (temperature > 1.5) thermostatinterval *= 0.5;

    counter++;

    createevent(time + writeinterval, NULL, NULL, 100);     //Add next write interval
}


/**************************************************
**                    CHECKCELLS
** Checks the cell list for possible errors
**************************************************/
void checkcells()
{
    int x, y, z, count = 0;
    double dl = 0.0000001;
    particle* p;
    for (x = 0; x < cx; x++)
        for (y = 0; y < cy; y++)
            for (z = 0; z < cz; z++)
            {
                p = celllist[x][y][z];
                if (p)
                {
                    if (p->prev) printf("First part has a prev (%d, %d, %d) %d, %d\n",
                        x, y, z, p->number, p->prev->number);
                    while (p)
                    {
                        count++;
                        if (p->cellx != x || p->celly != y || p->cellz != z)
                        {
                            printf("Cell error: %d, %d, %d / %d, %d, %d\n", x, y, z, p->cellx, p->celly, p->cellz);
                            exit(3);
                        }
                        if (p->x < cxsize * x - dl || p->x > cxsize * (x + 1) + dl)
                        {
                            printf("wrong cell x: %lf, %d, %d, %lf\n", p->x, x, p->number, p->vx);
                            exit(3);
                        }
                        if (p->y < cysize * y - dl || p->y > cysize * (y + 1) + dl)
                        {
                            printf("wrong cell y: %lf, %d, %d, %lf\n", p->y, y, p->number, p->vy);
                            exit(3);
                        }
                        if (p->z < czsize * z - dl || p->z > czsize * (z + 1) + dl)
                        {
                            printf("wrong cell z: %lf, %d, %d, %lf\n", p->z, z, p->number, p->vz);
                            exit(3);
                        }
                        if (p->next)
                        {
                            if (p->next->prev != p) printf("link error: %d, %d, %d\n",
                                p->number, p->next->number, p->next->prev->number);
                        }
                        p = p->next;
                    }
                }
            }
    if (count != N) printf("error in number of particles (%d)\n", count);
}

/**************************************************
**                    BACKINBOX
** Just for initialization
**************************************************/
void backinbox(particle* p)
{
    p->x -= xsize * floor(p->x / xsize);
    p->y -= ysize * floor(p->y / ysize);
    p->z -= zsize * floor(p->z / zsize);
}


/**************************************************
**                    THERMOSTAT
**************************************************/
void thermostat(event* ev)
{
    if (ev)
    {
        int i, num;
        particle* p;
        time = ev->time;
        int freq = N / 100;
        if (freq == 0) freq = 1;
        for (i = 0; i < freq; i++)
        {
            num = genrand_real2() * N;			//Random particle
            p = particles + num;
            double imsq = 1.0 / sqrt(p->mass);
            update(p);
            p->vx = random_gaussian() * imsq;			//Kick it
            p->vy = random_gaussian() * imsq;
            p->vz = random_gaussian() * imsq;
            p->counter++;
            removeevent(p->firstcollision);
            findcollisions(p);
        }
        removeevent(ev);
    }
    createevent(time + thermostatinterval, NULL, NULL, 200);     //Add next write interval
}

/**************************************************
**                    WRITELAST
**************************************************/
void writelast()
{
    int i;
    particle* p;
    FILE* file;
    char filename[200];
    sprintf(filename, "last.sph");
    file = fopen(filename, "w");
    fprintf(file, "%d\n%lf %lf %lf\n", (int)N, xsize, ysize, zsize);
    for (i = 0; i < N; i++)
    {
        p = &(particles[i]);
        fprintf(file, "%c %.12lf  %.12lf  %.12lf  %.12lf\n", 'a' + p->type, p->x, p->y, p->z, p->r);
    }
    fclose(file);
}

double random_gaussian()
{
    static int have_deviate = 0;
    static double u1, u2;
    double  x1, x2, w;

    if (have_deviate)
    {
        have_deviate = 0;
        return u2;
    }
    else
    {
        do
        {
            x1 = 2.0 * genrand_real2() - 1.0;
            x2 = 2.0 * genrand_real2() - 1.0;
            w = x1 * x1 + x2 * x2;
        } while (w >= 1.0);
        w = sqrt((-2.0 * log(w)) / w);
        u1 = x1 * w;
        u2 = x2 * w;
        have_deviate = 1;
        return u1;
    }
}


