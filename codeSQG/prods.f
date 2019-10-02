c	Prods.f
c
c       Product d=a*b, e=a*c
c
c       Dimension: (ND)
c
	subroutine prods(a,b,c,nd)
	implicit none

	integer nd
	real*8 a(nd),b(nd),c(nd)
	integer i

	do i=1,nd
	 b(i)=a(i)*b(i)
	 c(i)=a(i)*c(i)
	end do

	end
c --------------------------------------------------------------
