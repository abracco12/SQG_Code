c	
c	Checking whether the numerical stability is obtained, 
c	using the Courant number (CFL criteria). This number should be
c	less than 0.5.
c
	subroutine courant(u, v)

	implicit none

	include 'paran.h'
	include 'commfis.h'
	include 'commsim.h'

c
c	input:	u, v	velocity of field
c
c	output: cfl	the maximum courant number
c
	real*8	u(n,n)
	real*8  v(n,n)

	real*8	maxu, maxv, dx, dy
	real*8	tempu, tempv
	integer	i,j

	maxu = 0.0d0
	maxv = 0.0d0
	dx   = xlx/n
	dy   = xly/n
	do j=1,n
	  do i=1,n
	    tempu = dabs(u(i,j))
	    if (tempu.gt.maxu) maxu = tempu
	  enddo
	enddo
	maxu = maxu*dt/dx
	do j=1,n
	  do i=1,n
	    tempv = dabs(v(i,j))
	    if (tempv.gt.maxv) maxv = tempv
	  enddo
	enddo
	maxv = maxv*dt/dy
	if (maxu.gt.maxv) then
	  cfl = maxu
	else
	  cfl = maxv
	endif

	return
	end

