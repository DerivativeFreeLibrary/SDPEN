/* main program */
/***************************************************************
 *   Seqpen: Sequential Penalty Derivative-free Method for Nonlinear
 *   Constrained Optimization
 *
 *   This is an implementation in C of the algorithm described in
 *       G. Liuzzi, S. Lucidi, M. Sciandrone. Sequential Penalty 
 *       Derivative-free Methods for Nonlinear Constrained 
 *       Optimization, SIAM Journal on Optimization,
 *       20 (2010) 2614-2635.
 *   Copyright (C) K. Truemper
 *
 *   A Fortran90 implementation by G.Liuzzi, S.Lucidi, M.Sciandrone
 *   called SDPEN is available on the website 
 *   http://www.dis.uniroma1.it/~lucidi/DFL
 *
 *   This program is free software: You can redistribute it and/or 
 *   modify it under the terms of the GNU General Public License as 
 *   published by the Free Software Foundation, either version 3 of 
 *   the License, or (at your option) any later version.
 *
 *   For details of the GNU General Public License, see
 *   <http://www.gnu.org/licenses/>.
 *
 *   We do not make any representation or warranty, 
 *   expressed or implied, as to the condition,
 *   merchantability, title, design, operation, or fitness 
 *   of seqpen for a particular purpose.
 *
 *   We do not assume responsibility for any errors, 
 *   including mechanics or logic, in the operation 
 *   or use of seqpen, and have no liability whatsoever 
 *   for any loss or damages suffered by any user as a result of 
 *   seqpen.
 * 
 *   In particular, in no event shall we be liable
 *   for special, incidental, consequential, or
 *   tort damages, even if we have been advised of 
 *   the possibility of such damages.
 *
 **************************************************************/
# include <stdio.h>
# include <string.h>
# include <stdlib.h>
# include <math.h>

# include "basic.h"
# include "seqpen.h"
/*eject*/
int main () {
	
  int i;

  double finiz, violiniz; 

  /* initialize arrays */

  /* define nreal and ncon */
  setdim();

  /* allocate arrays; for activation, see seqpen.h */
  /* caution: all arrays are used starting with index = 1 */
  /* hence all dimensions are increased by 1 */
  /* alfa_d =  (double *)malloc((nreal+1)*sizeof(double));
  dir =     (double *)malloc((nreal+1)*sizeof(double));
  fstop =   (double *)malloc((nreal+2)*sizeof(double));

  lb =      (double *)malloc((nreal+1)*sizeof(double));
  ub =      (double *)malloc((nreal+1)*sizeof(double));
  xreal =   (double *)malloc((nreal+1)*sizeof(double));
  zvec =    (double *)malloc((nreal+1)*sizeof(double));

  eps =     (double *)malloc((ncon+1)*sizeof(double));
  epsiniz = (double *)malloc((ncon+1)*sizeof(double));
  constr =  (double *)malloc((ncon+1)*sizeof(double)); */ 

  /* define bounds lb and ub */
  setbounds();

  /* define starting xreal vector */
  startp(xreal);


/*eject*/
  /* evaluate initial solution */

  /* obj for xreal */
  obj = fobj(xreal);

  /* constr for xreal */
  fconstr(xreal);

  /* initialize eps using constr */
  for (i=1; i<=ncon; i++) {
    if (constr[i] < 1.0) {
      eps[i] = 1.e-3;
    } else {
      eps[i] = 1.e-1;
    }
  }

  /* compute violvalue */
  violvalue = 0.0;
  for (i=1; i<=ncon; i++) {
    violvalue = max(violvalue,constr[i]);
  }  

  /* display initial solution */
  displaysolution("initial");
/*eject*/
  /* initialize variables */
  num_funct = 0; 
  num_iter = 0;

  finiz = obj;
  violiniz = violvalue;

  alfa_stop = 1.e-6;
  nf_max = 5000;
  iprint = 2; /* = 0 if no iteration output */
              /* = 1 or 2 if iteration output */

  /* solve problem */
  printf("\nStart the optimizer");
  seqpen();
  printf("\nTermination condition: ");
  if (istop == 1) {
    printf("convergence\n");
  } else if (istop == 2) {
    printf("max number of evaluations\n");
  } else {
    printf(" seqpen: error, unknown termination code = %d",istop);
    exit(1);
  }
/*eject*/
  /* evaluate final solution */

  /* obj, constr, violvalue for xreal */
  obj = fobj(xreal);
  fconstr(xreal);
  violvalue = 0.0;
  for (i=1; i<=ncon; i++) {
    violvalue = max(violvalue,constr[i]);
  }

  /* display final solution */
  displaysolution("final");

  /* write LaTeX output line */
  printf("\n & %2d & %2d & %7d & %10.3e &",
         nreal,ncon,num_funct,finiz);
  printf("\n %8.1e & %10.3e & %8.1e \\\\ \\hline",
         violiniz,obj,violvalue);
  printf("\n------------------------\n");

  /* free arrays; for activation, see seqpen.h */
  /*free(alfa_d);
  free(dir);
  free(fstop);

  free(lb);
  free(ub);
  free(xreal);
  free(zvec);

  free(eps);
  free(epsiniz);
  free(constr);*/

  exit(0);

}
/************ last record of seqpenmain.c **************/
