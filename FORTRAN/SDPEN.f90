!============================================================================================
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
!    G. Liuzzi, S. Lucidi, M. Sciandrone. Sequential Penalty Derivative-free Methods for
!    Nonlinear Constrained Optimization, SIAM J. on Optimization, 20(5): 2614-2635 (2010)
!
!============================================================================================
subroutine sdpen(n,x,f,bl,bu,alfa_stop,nf_max,ni,nf,iprint,qq,eps)
	implicit none
	logical :: cambio_eps
	integer :: n,i,j,i_corr,nf,ni,nf_max,qq
	integer :: num_fal,istop
	integer :: iprint,i_corr_fall

	real*8 :: x(n),z(n),d(n)
	real*8 :: alfa_d(n),alfa,alfa_max
	real*8 :: f,fz 
	real*8 :: bl(n),bu(n),eps(qq-1),alfa_stop,maxeps 
	real*8 :: fstop(n+1)

	num_fal=0

	istop = 0

	fstop=0.d0

!     ---- scelta iniziale dei passi lungo le direzioni --------

      do i=1,n

        alfa_d(i)=dmax1(1.d-3,dmin1(1.d0,dabs(x(i))))
      
        if(iprint.ge.2) then
          write(*,*) ' alfainiz(',i,')=',alfa_d(i)
        endif

      end do
!     -----------------------------------------------------------


!     ---- scelta iniziale delle direzioni ----------------------

      do i=1,n      
        d(i)=1.d0 
      end do
!     -----------------------------------------------------------  
     
      call funct(n,x,f)
	  nf=nf+1

	  i_corr=1

      fstop(i_corr)=f

      do i=1,n
	    z(i)=x(i)
      end do

      if(iprint.ge.2) then
        write(*,*) ' ----------------------------------'
        write(*,*) ' finiz =',f
        do i=1,n
          write(*,*) ' xiniz(',i,')=',x(i)
        enddo
      endif

!---------------------------   
!     ciclo principale
!---------------------------

      do 

         if(iprint.ge.1) then
           write(*,*) '----------------------------------------------'
           write(*,100) ni,nf,f,alfa_max
100        format(' ni=',i4,'  nf=',i5,'   f=',d12.5,'   alfamax=',d12.5)
!           pause
         endif
         if(iprint.ge.2) then
	       do i=1,n
                write(*,*) ' x(',i,')=',x(i)
            enddo
         endif
!-------------------------------------
!    campionamento lungo asse i_corr
!-------------------------------------
 
         call linesearchbox(n,x,f,d,alfa,alfa_d,z,fz,i_corr,num_fal,&
                       alfa_max,i_corr_fall,iprint,bl,bu,ni,nf)
              
         if(dabs(alfa).ge.1.d-12) then
		               
            x(i_corr) = x(i_corr)+alfa*d(i_corr)
            f=fz
 	        fstop(i_corr)=f
			     
            num_fal=0
            ni=ni+1
      
         else
      
	        if(i_corr_fall.lt.2) then 

		      fstop(i_corr)=fz         

              num_fal=num_fal+1
              ni=ni+1

	        endif

	     end if

		 z(i_corr) = x(i_corr)

         if(i_corr.lt.n) then
            i_corr=i_corr+1
         else
            i_corr=1
         end if 

         call stop(n,alfa_d,istop,alfa_max,nf,ni,fstop,f,alfa_stop,nf_max)

         if (istop.ge.1) exit

         !------------------------------------------------
         ! Aggiornamento parametro di smoothing eps
         !------------------------------------------------
         cambio_eps=.false.
  	     maxeps = maxval(eps)
	     do i = 1,qq-1
		    if(eps(i) == maxeps) then
		       if(eps(i) > 1.0d-2*sqrt(alfa_max)) then
			  	 write(*,*) '**************************************'
				 write(*,*) '*********** aggiorno eps *************'
				 eps(i) =min(1.d-2*eps(i), 1.0d-1*sqrt(alfa_max))
				 cambio_eps=.true.
			     	 call funct(n,x,f)
			   endif
		    endif
	     enddo
         if(cambio_eps) then 
		   do i=1,n 
		      alfa_d(i)=dmax1(1.d-3,dmin1(1.d0,dabs(x(i)))) 
		   enddo
         endif
      enddo
      return
    


      end
        


!     #######################################################

      subroutine stop(n,alfa_d,istop,alfa_max,nf,ni,fstop,f,alfa_stop,nf_max)

      implicit none
      
      integer :: n,istop,i,nf,ni,nf_max
      real*8 :: alfa_d(n),alfa_max,fstop(n+1),ffstop,ffm,f,alfa_stop

      istop=0

      alfa_max=alfa_d(1)
      do i=1,n
        if(alfa_d(i).gt.alfa_max) then
          alfa_max=alfa_d(i)
        end if
      end do
     
      if(ni.ge.(n+1)) then
        ffm=f
        do i=1,n
          ffm=ffm+fstop(i)
        enddo
        ffm=ffm/dfloat((n+1))

        ffstop=(f-ffm)*(f-ffm)
        do i=1,n
           ffstop=ffstop+(fstop(i)-ffm)*(fstop(i)-ffm)
        enddo
 
        ffstop=dsqrt(ffstop/dfloat(n+1))

!        if(ffstop.le.alfa_stop) then
!         istop = 1
!       end if

	  endif

      if(alfa_max.le.alfa_stop) then
        istop = 1
      end if


      if(nf.gt.nf_max) then
        istop = 1
      end if

      return

      end




!     *********************************************************
!     *         
!     *         linesearch
!     *
!     ********************************************************
           
 
      subroutine linesearchbox(n,x,f,d,alfa,alfa_d,z,fz,i_corr,num_fal,&
                                 alfa_max,i_corr_fall,iprint,bl,bu,ni,nf)
      

      implicit none

      integer :: n,i_corr,nf
      integer :: i,j
      integer :: ni,num_fal
      integer :: iprint,i_corr_fall
	  integer :: ifront,ielle
      real*8 :: x(n),d(n),alfa_d(n),z(n),bl(n),bu(n)
      real*8 :: f,alfa,alfa_max,alfaex, fz,gamma
      real*8 :: delta,delta1,fpar,fzdelta

      gamma=1.d-6

      delta =0.5d0
      delta1 =0.5d0

      i_corr_fall=0

	  ifront=0

!     indice della direzione corrente

      j=i_corr

      if(iprint.ge.1) then
         write(*,*) ' j =',j,'    d(j) =',d(j)
      endif

      if(dabs(alfa_d(j)).le.1.d-3*dmin1(1.d0,alfa_max)) then
           alfa=0.d0
		   if(iprint.ge.1) then
               write(*,*) '  alfa piccolo'
               write(*,*) ' alfa_d(j)=',alfa_d(j),'    alfamax=',alfa_max
            endif
           return
      endif
      
	  do ielle=1,2

         if(d(j).gt.0.d0) then

            if((alfa_d(j)-(bu(j)-x(j))).lt.(-1.d-6)) then                 
               alfa=dmax1(1.d-24,alfa_d(j))
            else
               alfa=bu(j)-x(j)
               ifront=1
            endif

	     else

            if((alfa_d(j)-(x(j)-bl(j))).lt.(-1.d-6)) then
               alfa=dmax1(1.d-24,alfa_d(j))
            else
               alfa=x(j)-bl(j)
               ifront=1
            endif

         endif

         if(dabs(alfa).le.1.d-3*dmin1(1.d0,alfa_max)) then
  
            d(j)=-d(j)
            i_corr_fall=i_corr_fall+1
			alfa=0.d0
            ifront=0

            if(iprint.ge.1) then
               write(*,*) ' direzione opposta per alfa piccolo'
			   write(*,*) ' j =',j,'    d(j) =',d(j)
               write(*,*) ' alfa=',alfa,'    alfamax=',alfa_max
            endif
            
            cycle

         endif

         alfaex=alfa

         z(j) = x(j)+alfa*d(j)
       
         call funct(n,z,fz)
	     nf=nf+1

         if(iprint.ge.1) then
            write(*,*) ' fz =',fz,'   alfa =',alfa
         endif
         if(iprint.ge.2) then
            do i=1,n
              write(*,*) ' z(',i,')=',z(i)
            enddo
         endif

         fpar= f-gamma*alfa*alfa

         if(fz.lt.fpar) then

!           espansione

            do

               if((ifront.eq.1).or.(num_fal.gt.n-1)) then
!			   if((ifront.eq.1)) then

	              alfa_d(j)=delta*alfa

                  return

               end if

               if(d(j).gt.0.d0) then
		                
                  if((alfa/delta1-(bu(j)-x(j))).lt.(-1.d-6)) then
                     alfaex=alfa/delta1
                  else
                     alfaex=bu(j)-x(j)
                     ifront=1
                     if(iprint.ge.1) then
                        write(*,*) ' punto espan. sulla front.'
                     endif
                  end if

               else

                  if((alfa/delta1-(x(j)-bl(j))).lt.(-1.d-6)) then
                     alfaex=alfa/delta1
                  else
                     alfaex=x(j)-bl(j)
                     ifront=1
                     if(iprint.ge.1) then
                        write(*,*) ' punto espan. sulla front.'
                     endif
                  end if

               endif
		             
               z(j) = x(j)+alfaex*d(j)    
		    
               call funct(n,z,fzdelta)
	           nf=nf+1

               if(iprint.ge.1) then
                  write(*,*) ' fzex=',fzdelta,'  alfaex=',alfaex  
               endif
               if(iprint.ge.2) then
                  do i=1,n
                     write(*,*) ' z(',i,')=',z(i)
                  enddo
               endif

               fpar= f-gamma*alfaex*alfaex

               if(fzdelta.lt.fpar) then

                  fz=fzdelta
                  alfa=alfaex

               else               
                  alfa_d(j)=delta*alfa
                  return
               end if

            enddo

         else 

	        d(j)=-d(j)
            ifront=0

            if(iprint.ge.1) then
               write(*,*) ' direzione opposta'
	       write(*,*) ' j =',j,'    d(j) =',d(j)
            endif

         endif
		  
      enddo

      if(i_corr_fall.eq.2) then
         alfa_d(j)=alfa_d(j)
      else
         alfa_d(j)=delta*alfa_d(j)
      end if

      alfa=0.d0
      return         

      end



