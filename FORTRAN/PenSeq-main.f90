module vincoli
	integer			:: m
	real*8, allocatable	:: eps(:)
end module vincoli

program penseq
	use vincoli	
	implicit none
	
	integer 		:: num_funct, num_iter, i, n, qq
	integer			:: nf_max, iprint, icheck

	real*8 			:: f, alfamax, fob
	real*8			:: violiniz, finiz, alfa_stop

	real*8, allocatable 	:: constr(:),lb(:),ub(:)
	real*8, allocatable	:: x(:), epsiniz(:)

	character*15 nomefun

	common /num/f
	common /calfamax/alfamax

	!---------------------------------------
	!  bounds, punto iniziale, valore ottimo
	!---------------------------------------

	call setdim(n,m)
	qq = m+1

	allocate(constr(qq-1),lb(n),ub(n),x(n), eps(qq-1), epsiniz(qq-1))

	call setbounds(n,lb,ub)
	call startp(n,x)

	call funob(n,x,fob)
	call fconstr(n,m,x,constr)

	write(*,*) '------------------------'
	write(*,*) '---- Initial values ----'
	write(*,*) '------------------------'
	write(*,*) 'fob = ',fob

	write(*,*)
	if(n <= 10) then
		do i=1,n
			write(*,*) 'x(',i,')=',x(i)
		enddo
		write(*,*)
	endif

	do i = 1,qq-1
		if(max(0.d0,constr(i)) < 1.d-0) then
			eps(i) = 1.d-3
		else
			eps(i) = 1.d-1
		endif
	enddo

	if(m <=10) then
		do i=1,m
			write(*,*) 'con(',i,')=',constr(i),' eps(',i,')=',eps(i)
		enddo
	endif
	write(*,*) '------------------------'

	num_funct   = 0 
	num_iter    = 0

	epsiniz     = eps
	finiz       = fob
	violiniz    = max(0.d0,maxval(constr))

1	alfa_stop=1.d-6
	nf_max=5000
	iprint=0

	write(*,*)
	write(*,*) 'Start the optimizer:'

	call sdpen(n,x,f,lb,ub,alfa_stop,nf_max,num_iter,num_funct,iprint,qq,eps)

	write(*,*) 'Done'
	write(*,*)


	call funob(n,x,fob)
	call fconstr(n,m,x,constr)
		
	write(*,*) '------------------------'
	write(*,*) '----   Final values ----'
	write(*,*) '------------------------'
	write(*,*)
	write(*,*) 'fob  = ',fob

	write(*,*)
	if(n <= 10) then
		do i=1,n
			write(*,*) 'x(',i,')=',x(i)
		enddo
		write(*,*)
	endif

	if(m <=10) then
		do i=1,m
			write(*,*) 'con(',i,')=',constr(i),' eps(',i,')=',eps(i)
		enddo
	endif
	write(*,*) '------------------------'

	write(*,1010) n,m,num_funct,	&
		      finiz,violiniz,fob,max(0.d0,maxval(constr))

	write(*,*) '------------------------'

1010 format(1x,2(' & ', i2),1(' & ',i7), ' & ',es10.3,' & ',es8.1,' & ',es10.3,' & ',es8.1,' \\\\ \\hline'  )
END program penseq

!===================================================================
!===================================================================

subroutine funct(n,x,f)
	use vincoli
	implicit none
	integer n,i
	real*8  x(n),f,fob,fmax
	real*8	constr(m)
	
	call funob(n,x,fob)
	call fconstr(n,m,x,constr)

	fmax = 0.d0

	do i = 1,m
		fmax = fmax + (max(0.d0,constr(i))**1.1d0)/eps(i)
	enddo

	f = fob + fmax

	return
end subroutine funct

