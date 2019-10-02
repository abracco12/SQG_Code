c	Forcing
c
c	Forza ad un dato numero d onda imponendo che
c       il modulo della trasformata di Fourier sia pari 
c       a un valore dato
c
	subroutine forcing(a)
	implicit none
        include "paran.h"
        include "parafor.h"

	real*8 a(n,n), fase
	integer j

c        do j=1,n2p1
c          print *,j-1,a(1,j),a(2,j)
c        enddo

        fase=a(2,kkk+1)/a(1,kkk+1)
        a(1,kkk+1)=dsqrt(fkcost/(1+fase**2))
        a(2,kkk+1)=fase*a(1,kkk+1)
	
        return
	end
c --------------------------------------------------------------------
