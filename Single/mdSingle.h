
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
	unsigned int counter;		//Number of collision events experienced
	unsigned int counter2;		//Number of collision events of collision partner at time of scheduling
	struct sparticle* prev, *next;	//Doubly linked cell list
	uint8_t type;				//Particle type
	double eventtime;
	struct sparticle* left;		//Left child in tree or previous event in event list
	struct sparticle* right;	//Right child in tree or next event in event list
	struct sparticle* parent;	//Parent in tree
	struct sparticle* p2;		//Collision partner
	int queue;					//Index of the event queue the event is in
	unsigned char eventtype;
} particle;

int main();
void printstuff();
void init();


void initevents();
void fcc();
void loadparticles();
void randommovement();
void initcelllist();
void addtocelllist(particle* p, int cellx, int celly, int cellz);
int celloffset(int a, int b, int c);

void step();
int findcollision(particle*, particle*, double*);
void findallcollisions();
void findcollisions(particle*);
void collision(particle*);

void addevent(particle*);
void removeevent(particle*);
void createevent(double time, particle* p1, particle* p2, int type);
void addnexteventlist();
double findneighborlistupdate(particle* p1);
void makeneighborlist(particle* p1);

void outputsnapshot();
void write(particle* ev);
void thermostat(particle* ev);
double random_gaussian();
void backinbox(particle* p);
