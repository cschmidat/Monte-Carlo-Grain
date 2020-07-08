
#define unlikely(expr) __builtin_expect(!!(expr),0)
#define likely(expr) __builtin_expect(!!(expr),1)
#define POSORNE(a,b) ((pts[a].nb[0]==b)?0:((pts[a].nb[1]==b)?1:2))
#define DOUBCON(a) (pts[a].nb[0]==pts[a].nb[1]||pts[a].nb[0]==pts[a].nb[2]||pts[a].nb[1]==pts[a].nb[2])

/*
Simulation: 
Events: T1...Swap
T2...Disappearance of small grain
Jump length is chosen s.t. 30% of jumps pass without considering Boltzmann factor
Temperatur is chosen s.t. points fluctuate an average of 1% of the typical initial length.
*/


int read(char *, char *, struct pt *);
int jump(int i,struct pt *points);
int write(char *suffix, struct pt *pts);
