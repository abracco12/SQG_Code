c	Lire.f
c
c	Read n real*8 from the file ku  and set vector u.
c
c	(file ASCII)
c
	subroutine lire(ku,u,nd)
	implicit none

	integer ku,nd,i
	real*8 u(nd)

	do i=1,nd
           read(ku,*) u(i)
	enddo

	return
	end
c ---------------------------------------------------------------
