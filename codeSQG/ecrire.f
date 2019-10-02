c	Ecrire.f
c
c	Write an unformatted file of dimension nd.
c	
c
	subroutine ecrire(ku,u,nd)
	implicit none

	integer nd,i
	real*8 u(nd)
	character*12 ku

        open (unit=50,status='unknown',file=ku,
     -       form='unformatted',access='append')
	      write(50) (real(u(i)), i=1,nd)
	close(50)

	return
	end
c ---------------------------------------------------------------
