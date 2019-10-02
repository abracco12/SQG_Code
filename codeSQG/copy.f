c	Copy.f
c
c	To copy the first vector in the second.
c
	subroutine copy(a,b,nd)
	implicit none

	integer nd
	real*8 a(nd),b(nd)
	integer i

	do i=1,nd
	 b(i)=a(i)
	end do

	return
	end
c --------------------------------------------------------------------
