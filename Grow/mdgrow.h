// Event driven MD, headers.

typedef struct sevent
{
	double time;
	struct sevent* child1;
	struct sevent* child2;
	struct sevent* parent;
	struct sparticle* p1;
	struct sparticle* p2;
	struct sevent* prevq, * nextq;
	uint8_t type;
	int queue;
	unsigned int counter2;
} event;

typedef struct sparticle
{
	double x, y, z;
	double vx, vy, vz;
	double xn, yn, zn;
	double vr;
	struct sparticle* neighbors[MAXNEIGH];
	uint8_t nneigh;
	double t;
	double r;
	double rtarget;
	double mass;
	uint8_t edge;
	uint8_t cellx, celly, cellz;
	int boxestraveledx, boxestraveledy, boxestraveledz;
	event* firstcollision;
	unsigned int counter;
	struct sparticle* prev, * next;
	int number;
	uint8_t type;
} particle;

int main();
void printstuff();
void init();
void cubicfcc();
void cubicfcc2();
void cubicfcc3();

void initeventpool();
void fcc();
void randomparticles();
void randommovement();
void initcelllist();
void addtocelllist(particle*);

void step();
double findcollision(particle*, particle*, double);
void findallcollisions();
void findcollisions(particle*);
void collision(event*);

void addevent(event*);
void removeevent(event*);
event* createevent(double time, particle* p1, particle* p2, int type);
void addnexteventlist();
double findneighborlistupdate(particle* p1);
void makeneighborlist(particle* p1, int firsttime);



void showtree();
void shownode(event*);
void checktree();
int checknode(event*);
int overlap(particle*);
int overlaplist(particle* part, int error);
void outputsnapshot();
void write();
void writelast();
void checkcells();
void thermostat(event* ev);
double random_gaussian();
void backinbox(particle* p);
void writelast();
