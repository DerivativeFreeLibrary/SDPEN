/* routines defining problem */

# include <stdio.h>
# include <stdlib.h>
# include <math.h>

# include "basic.h"
# include "seqpen.h"

/***************************************************************
 *
 * subroutines in this file:
 *       void fconstr(double *x)
 *       void fobj(double *x)
 *       void setbounds()
 *       void setdim()
 *       void startp(double *x)
 **************************************************************/
/*eject*/
/**************************************************************
 *   void fconstr(double *x): defines constr values for given 
 *     problem
 *     all constraints are assumed to be <= 0
 *   uses or modifies: nreal, ncon, constr          
 **************************************************************/
void fconstr(double *x) {

  constr[1] = -x[1] - pow(x[2],2.0);
  constr[2] = -x[2] - pow(x[1],2.0);
  constr[3] = 1.0 - pow(x[1],2.0) - pow(x[2],2.0);

  return;

}
/*eject*/
/**************************************************************
 *   double fobj(double *x): returns obj value for given vector x
 *   uses or modifies: nreal       
 **************************************************************/
double fobj(double *x){

  num_funct++;
  return (100.0 * pow((x[2] - pow(x[1],2.0)),2.0) + 
          pow((1.0-x[1]),2.0));

}
/*eject*/
/**************************************************************
 *   void setbounds(): defines lb and ub bounds for given problem
 *   uses or modifies: nreal,lb,ub
 **************************************************************/
void setbounds() {

  int j;

  for (j=1; j<=nreal; j++) {
    lb[j] = -1.e+6;
    ub[j] = 1.e+6;
  }
  lb[1] = -0.5;
  ub[1] = 0.5;

  return;

}
/*eject*/
/**************************************************************
 *   void setdim(): defines number of variables nreal
 *                          number of constraints ncon
 *   uses or modifies: nreal, ncon
 **************************************************************/
void setdim() {

  nreal = 2;
  ncon = 3;

  return;

}
/*eject*/
/**************************************************************
 *   void startp(double *x): defines initial vector x
 *   uses or modifies: nreal
 **************************************************************/
void startp(double *x) {

  x[1] =-0.5;
  x[2] = 1.0;

  return;

}
/***************** last record of problem.c ****************/

