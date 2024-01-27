/*============================================================================================
!    SDPEN - A Sequential Penalty Derivative-free Method for
!    Nonlinear Constrained Optimization
!    Copyright (C) 2011  G.Liuzzi, S.Lucidi, M.Sciandrone
!
!    This program is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    This program is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
!    The C version of the SDPEN has been kindly provided by Prof. 
!    Klaus Truemper from University of Texas at Dallas
!
!    G. Liuzzi, S. Lucidi, M. Sciandrone. Sequential Penalty Derivative-free Methods for
!    Nonlinear Constrained Optimization, SIAM J. on Optimization, 20(5): 2614-2635 (2010)
!
!!============================================================================================
*/
/* routines solving problem */
# include <stdio.h>
# include <stdlib.h>
# include <math.h>

# include "basic.h"
# include "seqpen.h"
/**************************************************************
 *
 * subroutines in this file:
 *       void seqpen()
 *       void displaysolution(char *label)
 *       void funct()
 *       void linesearchbox()
 *       void stop()
 **************************************************************/
/*eject*/
/**************************************************************
 *   void seqpen(): main solution subroutine
 *   uses or modifies: nreal,xreal,obj,lb,ub,alfa_stop,nf_max,
 *                     num_iter,num_funct,iprint,eps          
 **************************************************************/
void seqpen() {

  int cambio_eps, i;

  double maxeps; 

  num_fal = 0;
  istop = 0;

  for (i=1; i<=nreal+1; i++) {
    fstop[i] = 0.0;
  }

  /* initialize alfa_d and alfa_max */
  for (i=1; i<=nreal; i++) {
    alfa_d[i] = max(1.e-3,min(1.0,fabs(xreal[i])));    
    if (iprint >= 2) {
      printf("\nalfainiz[%d] = %g",i, alfa_d[i]);
      fflush(stdout);
    }
  }
  alfa_max = alfa_d[1];
  for (i=1; i<=nreal; i++) {
    alfa_max = max(alfa_max,alfa_d[i]);
  }
  /* initialize dir */
  for (i=1; i<=nreal; i++) {      
    dir[i] = 1.0; 
  }  
  /* obj and constr for xreal */   
  obj = funct(xreal);
  i_corr = 1;
  fstop[i_corr] = obj;

  for (i=1; i<=nreal; i++) {
    zvec[i]=xreal[i];
  }

  if (iprint >= 2) {
    printf("\n ----------------------------------");
    printf("\n finiz = %g",obj);
    for (i=1; i<=nreal; i++) {
      printf("\n xiniz[%d] = %g",i,xreal[i]);
      fflush(stdout);
    }
  }

  /* begin main loop */

  while (1) { /* while #1 */

    if (iprint >= 1) {
      printf("\n------------------------------");
      printf("\n num_iter = %d num_funct = %d obj = %g alfamax = %g",
             num_iter,num_funct,obj,alfa_max);
      fflush(stdout);
    }
    if (iprint >=2 ) {
      for (i=1; i<= nreal; i++) {
        printf("\n xreal[%d] = %g",i,xreal[i]);
        fflush(stdout);
      }
    }

    /* search along axis i_corr */
    linesearchbox();

    /* update xreal and obj */    
    if(fabs(alfa) >= 1.e-12) {
      xreal[i_corr] += alfa*dir[i_corr];
      obj = zobj;
      fstop[i_corr] = obj;
      num_fal = 0;
      num_iter++;     
    } else {
      if (i_corr_fall <= 2) { 
        fstop[i_corr] = zobj;         
        num_fal++;
        num_iter++;
      }
    }

    /* update zvec */
    zvec[i_corr] = xreal[i_corr];

    /* update i_corr */
    if (i_corr < nreal) {
      i_corr++;
    } else {
      i_corr = 1;
    }

    /* test for stop condition */
    teststop();
    /* termination; istop >= 1 is condition code */
    if (istop >= 1) {
      return;
    }

    /* update eps */       
    cambio_eps = FALSE;
    maxeps = -INF;
    for (i=1; i<=ncon; i++) {
      maxeps = max(maxeps,eps[i]);
    }
    for (i=1; i<=ncon; i++) {
      if(eps[i] == maxeps) {
        if(eps[i] > 1.0e-2*sqrt(alfa_max)) {
          printf("\n**************************************");
          printf("\n************ update eps **************");
          eps[i] = min(1.e-2*eps[i], 1.0e-1*sqrt(alfa_max));
          cambio_eps = TRUE;
          obj = funct(xreal);
        }
      }
    }
    if (cambio_eps == TRUE) { 
      for (i=1; i<=nreal; i++) {
        alfa_d[i] = max(1.e-3,min(1.0,fabs(xreal[i]))); 
      }
    }

  } /* end while #1 */

  return;

}
/*eject*/
/**************************************************************
 *   void displaysolution(char *label): display
 *   obj, xreal, constr, eps, violvalue
 *   under the heading given by label
 **************************************************************/
void displaysolution(char *label) {

  int i;

  printf("\n------------------------");
  printf("\n---- %s values ----",label);
  printf("\n------------------------");
  /* obj */
  printf("\nobj = %g\n",obj);
  /* xreal */
  for (i=1; i<=nreal; i++) {
    printf("\nxreal[%d] = %g",i,xreal[i]);
  }
  printf("\n");
  /* constr and eps */
  for (i=1; i<=ncon; i++) {
    printf("\nconstr[%d] = %g\teps[%d] = %g",
           i,constr[i],i,eps[i]);
  }
  printf("\n");
  /* violvalue */
  printf("\nmax violation value = %g",violvalue);
  if (violvalue == 0.0) {
    printf("\nsolution is feasible");
  } else {
    printf("\nsolution is infeasible");
  } 
  printf("\n------------------------");
  printf("\n");
  fflush(stdout);

  return;

}
/*eject*/
/**************************************************************
 *   double funct(double *x): evaluate obj function and 
 *          constraints and return combined measure fmax
 *   uses or modifies: nreal, constr          
 **************************************************************/
double funct(double *x) {

  int i;
  double fmax;
  /* initialize fmax as obj value of x */
  fmax = fobj(x);
  /* compute constr for x */
  fconstr(x);
  /* update fmax using constr and eps */
  for (i=1; i<=ncon; i++) {
      if (eps[i] <= 0.0) {
        printf("\n funct: error, eps[%d] = %g <= 0.0",i,eps[i]);
        exit(1);
      }
      fmax += pow(max(0.0,constr[i]),1.1)/eps[i];
  }

  return fmax;

}
/*eject*/
/**************************************************************
 *   void linesearchbox(): carry out line search
 *   uses or modifies: nreal,xreal,obj,dir,alfa,alfa_d,
 *                     zvec,zobj,i_corr,
 *                     num_fal,alfa_max,i_corr_fall,iprint,
 *                     lb,ub,num_iter,num_funct          
 **************************************************************/
void linesearchbox() {

  int i,j;
  int ifront,ielle;
  double alfaex,gamma;
  double delta,delta1,fpar,fzdelta;

  gamma =1.e-6;
  delta = 0.5;
  delta1 = 0.5;

  i_corr_fall = 0;
  ifront = 0;

  /* index of current direction */
  j = i_corr;

  if (iprint >= 1) {
    printf("\n j = %d    dir[j] = %g",j,dir[j]);
    fflush(stdout);
  }

  if(fabs(alfa_d[j]) <= 1.e-3*min(1.0,alfa_max)) {
    alfa=0.0;
    if (iprint >= 1) {
      printf("\n alfa small");
      printf("\n alfa_d[j] = %g    alfamax = %g",
             alfa_d[j], alfa_max);
      fflush(stdout);
    }
      return;
  }
     
  for (ielle=1; ielle<=2; ielle++) {

    /* update alfa, ifront */
    if (dir[j] > 0.0) {
      if((alfa_d[j]-(ub[j]-xreal[j])) < (-1.e-6)) {
        alfa = max(1.e-24,alfa_d[j]);
      } else {
        alfa = ub[j]-xreal[j];
        ifront = 1;
      }
    } else {
      if((alfa_d[j]-(xreal[j]-lb[j])) <(-1.e-6)) {
        alfa = max(1.e-24,alfa_d[j]);
      } else {
        alfa = xreal[j]-lb[j];
        ifront = 1;
      }
    }

    if (fabs(alfa) <= 1.e-3*min(1.0,alfa_max)) {
      dir[j] = -dir[j];
      i_corr_fall++;
      alfa = 0.0;
      ifront = 0;
       if(iprint >= 1) {
         printf("\nopposite direction to alfa small");
         printf("\n j = %d    dir[j] = %g    ",j,dir[j]);
         printf("alfa = %g    alfamax = %g",alfa,alfa_max);
         fflush(stdout);
       }           
      continue; 
    }
    alfaex = alfa;

    /* update zvec[j] */
    zvec[j] = xreal[j] + alfa*dir[j];
    zobj = funct(zvec);
    if (iprint >= 1) {
      printf("\nzobj = %g   alfa = %g",zobj, alfa);
      fflush(stdout);
    }
    if (iprint >= 2) {
      for (i=1; i<=nreal; i++) {
        printf("\nzvec[%d] = %g",i,zvec[i]);
        fflush(stdout);
      }
    }

    fpar = obj-gamma*alfa*alfa;

    if (zobj < fpar) {
      /* expansion */
      while (1) { /* while # 2 */
        if ((ifront == 1) || (num_fal > nreal-1)) {
          alfa_d[j] = delta * alfa;
          return;
        }

        if (dir[j] > 0.0) {
          if((alfa/delta1-(ub[j]-xreal[j])) < (-1.e-6)) {
            alfaex = alfa/delta1;
          } else {
            alfaex = ub[j]-xreal[j];
            ifront = 1;
            if(iprint >= 1) {
              printf("\npoint expansion on the front");
              fflush(stdout);
            }
          }
        } else {
          if ((alfa/delta1-(xreal[j]-lb[j])) < (-1.e-6)) {
            alfaex=alfa/delta1;
          } else {
            alfaex = xreal[j]-lb[j];
            ifront = 1;
            if (iprint >= 1) {
              printf("\npoint expansion on the front");
              fflush(stdout);
            }
          }
        }
		             
        zvec[j] = xreal[j] + alfaex*dir[j];
        fzdelta = funct(zvec);
        if (iprint >= 1) {
          printf("\nfzex =%g  alfaex = %g",fzdelta,alfaex);
          fflush(stdout); 
        }
        if(iprint >= 2) {
          for (i=1; i<=nreal; i++) {
            printf("\nzvec[%d] = %g",i,zvec[i]);
            fflush(stdout);
          }
        }

        fpar = obj - gamma*alfaex*alfaex;
        if (fzdelta < fpar) {
          zobj = fzdelta;
          alfa = alfaex;
        } else {               
          alfa_d[j] = delta*alfa;
          return;
        }
      } /* end while #2 */

    } else { 
      dir[j] = -dir[j];
      ifront = 0;
      if (iprint >= 1) {
        printf("\nopposite direction ");
	printf("\nj = %d    dir[j] = %g",j,dir[j]);
        fflush(stdout);
      }
    }
  } /* end for ielle=1 */

  if (i_corr_fall == 2) {
    alfa_d[j] = alfa_d[j];
  } else {
    alfa_d[j] = delta*alfa_d[j];
  }

   alfa = 0.0;

   return;

}
/*eject*/
/**************************************************************
 *   void teststop(): test termination conditions
 *   uses or modifies:  nreal,alfa_d,istop,alfa_max,num_funct,
 *                      num_iter,fstop,obj,alfa_stop,nf_max          
 **************************************************************/
void teststop() {

  int i;

  double ffstop,ffm;

  istop = 0;

  alfa_max = alfa_d[1];
  for (i=1; i<=nreal; i++) {
    alfa_max = max(alfa_max,alfa_d[i]);
  }
     
  if (num_iter >= (nreal+1)) {
    ffm = obj;
    for (i=1; i<=nreal; i++) {
       ffm = ffm + fstop[i];
    }
    ffm = ffm/((double)(nreal+1));
    ffstop = (obj-ffm) * (obj-ffm);
    for (i=1; i<=nreal; i++) {
      ffstop += (fstop[i]-ffm)*(fstop[i]-ffm);
    }
    ffstop=sqrt(ffstop/((double)(nreal+1)));
    /* if(ffstop <= alfa_stop) {
      istop = 1;
    } */
  }

  if(alfa_max <= alfa_stop) {
    istop = 1;
  }
  if(num_funct > nf_max) {
    istop = 2;
  }

  return;

}
/************** last record of seqpen.c *************/

