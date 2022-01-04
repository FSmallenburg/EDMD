

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
	uint8_t eventtype;
	int queue;					//Index of the event queue this event is in
} event;

//Particle structure
typedef struct sparticle
{
	double x, y, z;				//Position
	double vx, vy, vz;			//Velocity
	double xn, yn, zn;			//Neighbor list center
	struct sparticle* neighbors[MAXNEIGH];
	uint8_t nneigh;				//Number of neighbors
	double t;					//Last update time
	double radius;
	double mass;
    uint8_t nearboxedge;		//Is this particle in a cell near the box edge?
	int cell;					//Current cell
	int boxestraveledx, boxestraveledy, boxestraveledz;	//Keeps track of dynamics across periodic boundaries
	event* cellcrossing;		//Cell crossing event
	struct sparticle* prev, *next;	//Doubly linked cell list
	uint8_t type;				//Particle type
} particle;

int main();
void printstuff();
void init();

void initeventpool();
void fcc();
void loadparticles();
void randommovement();
void initcelllist();
void addtocelllist(particle* p, int cellx, int celly, int cellz);
void step();
void findcollision(particle*, particle*);
void findallcollisions();
void findcollisions(particle*, particle*);
void collision(event*);

void findlistupdate(particle*);
void addevent(event*);
void removeevent(event*);
void createevent(double time, particle* p1, particle* p2, int type);
void addnexteventlist();
double findneighborlistupdate(particle* p1);
void makeneighborlist(particle* p1, int firsttime);


void outputsnapshot();
void write();
void thermostat();
double random_gaussian();
void backinbox(particle* p);
