subroutine setdim(n,m)
	implicit none
	integer		:: n, m

	n = 2
	m = 3

	return
end subroutine setdim

subroutine funob(n, x, f)
	implicit none
	integer		:: n
	real*8		:: x(n), f

	f = 100.d0*(x(2) - x(1)**2.d0)**2.d0 + (1.d0-x(1))**2.d0

	return
end subroutine funob

subroutine fconstr(n, m, x, constr)
	implicit none
	integer		:: n, m
	real*8		:: x(n), constr(m)

	integer		:: i

	constr(1) = -x(1) -x(2)**2.d0
	constr(2) = -x(2) -x(1)**2.d0
	constr(3) = 1.d0 -x(1)**2.d0 -x(2)**2.d0

	return
end subroutine fconstr

subroutine startp(n,x)
	implicit none
	integer		n
	real*8		x(n)

	X(1) =-0.5D0
	X(2) = 1.D0

	return
end subroutine startp

subroutine setbounds(n,lb,ub)
	implicit none
	integer		n
	real*8		lb(n), ub(n)
	lb=-1.d+6
	ub=1.d+6
	lb(1)=-0.5d0
	ub(1)=0.5d0
	return
end subroutine setbounds

