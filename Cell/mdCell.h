

//Event structure
typedef struct sevent
{
    double eventtime;
    struct sevent* left;		//Left child in tree or previous event in event list
    struct sevent* right;		//Right child in tree or next event in event list
    struct sevent* parent;		//Parent in tree
    struct sparticle* p1;		//Particles involved in the event
    struct sparticle* p2;
    struct sevent* prevp1, *nextp1;	//Circular linked list for all events involving p1 as the first particle
    struct sevent* prevp2, *nextp2; //Circular linked list for all events involving p2 as the second particle
    int eventtype;
    int queue;					//Index of the event queue this event is in
} event;

//Particle structure
typedef struct sparticle
{
    double x, y, z;				//Position
    double vx, vy, vz;			//Velocity
  	double t;					//Last update time
    double radius;
    double mass;
    uint8_t nearboxedge;		//Is this particle in a cell near the box edge?
    int cellx, celly, cellz; 					//Current cell
  	int boxestraveledx, boxestraveledy, boxestraveledz;	//Keeps track of dynamics across periodic boundaries
  	event* cellcrossing;		//Cell crossing event
	  struct sparticle* prev, *next;	//Doubly linked cell list
	  uint8_t type;				//Particle type
} particle;

int main();
void printstuff();
void init();

void randomparticles();
void initeventpool();
void fcc();
void loadparticles();
void randommovement();
void initcelllist();
void addtocelllist(particle* p);
void step();
void findcollision(particle*, particle*);
void findallcollisions();
void findcollisions(particle*, particle*);
void findcollisioncell(particle*, int);
void findcrossing(particle* part);
void collision(event*);
void fakecollision(event*);
void addfakecollision(particle*);
void findlistupdate(particle*);
void cellcross(event*);
void addevent (event*);
void createevent (double time, particle* p1, particle* p2, int type);
void addnexteventlist();

void removeevent (event*);
void outputsnapshot();
void write();
void thermostat();
double random_gaussian();
void backinbox(particle* p);
