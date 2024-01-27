/* constants */

/* used only for static allocation */
#define MAX_CONSTRAINT 50
#define MAX_VARIABLE 50

/* other constants */
#define INF 1.0e14
#define TRUE  1
#define FALSE 0

/* basic variables for seqpen */

 int nreal; /* number of variables */
 int ncon;  /* number of constraints */

/* max and min definition */
#define max(x,y) (((x) > (y)) ? (x) : (y))
#define min(x,y) (((x) < (y)) ? (x) : (y))

/* math library */
double pow(double,double);

/************** last record of basic.h *************/
