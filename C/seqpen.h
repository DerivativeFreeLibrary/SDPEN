/* variables, arrays, and routines for seqpen */

  /********************* variables ***********************/
  int nf_max, num_fal, num_iter, num_funct;
  int iprint, istop; 
  int i_corr, i_corr_fall;
  double alfa, alfa_max, alfa_stop, obj, violvalue, zobj;

  /********************* arrays **************************/
 /* for dynamic allocation:
  * (1) uncomment the left hand side definitions, and comment
  *     out the right hand side definitions
  * (2) in seqpmain.c: uncomment 'malloc' statements at beginning
  *     of file and 'free' statements at end of file
  *
  * caution: all arrays are used starting with index = 1
  * hence all dimensions must be increased by 1
  *
  *  dynamic alloc     min dimension      static allocation
  *-------------------------------------------------------*/
  /*double *alfa_d;*/  /* nreal+1 */ double alfa_d[MAX_VARIABLE+1];
  /*double *dir;*/     /* nreal+1 */ double dir[MAX_VARIABLE+1];
  /*double *fstop;*/   /* nreal+2 */ double fstop[MAX_VARIABLE+2];

  /*double *lb;*/      /* nreal+1 */ double lb[MAX_VARIABLE+1];
  /*double *ub;*/      /* nreal+1 */ double ub[MAX_VARIABLE+1];
  /*double *xreal;*/   /* nreal+1 */ double xreal[MAX_VARIABLE+1];
  /*double *zvec;*/    /* nreal+1 */ double zvec[MAX_VARIABLE+1];

  /*double *eps;*/     /* ncon+1  */ double eps[MAX_CONSTRAINT+1];
  /*double *epsiniz;*/ /* ncon+1  */ double espiniz[MAX_CONSTRAINT+1];
  /*double *constr;*/  /* ncon+1  */ double constr[MAX_CONSTRAINT+1];
/*eject*/
/************************ routines ************************/

/* seqpenmain.c */
/* main program */

/* seqpen.c */
void seqpen();
void displaysolution(char *label);
double funct(double *x);
void linesearchbox();
void teststop();

/* problem.c */
void fconstr(double *x);
double fobj(double *x);
void setbounds();
void setdim();
void startp(double *x);

/************* last record of seqpen.h **********/
