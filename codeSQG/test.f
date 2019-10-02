c	***********************************************
c
c	Write energy and enstrophy
c	Z in spectral space
c
	subroutine test(z,t)
	implicit none

	include "paran.h"
	include "commk.h"
	include "commfis.h"
	include "commsim.h"
	include "paraio.h"
	include "paraqg.h"

	real*8 z(N,N)
	real*8 t,k2
	real*8 ee,zz,uu
	integer i,j

	ee= 0.0
	zz= z(1,1)**2 + z(2,1)**2

	do j = 1, N
	   do i = 3, N, 2
	      k2 = fks(i/2+1) + fks(j) + F
	      uu = z(i,j)**2 + z(i+1,j)**2
	      zz = zz + uu
	      ee = ee + uu/k2
	   enddo
	enddo

	do j = 2, N2P1
	   uu = z(1,j)**2 + z(2,j)**2
	   zz = zz + uu
	   ee = ee + uu/(fks(j) + F)
	enddo

	do j = N2P2, N
	   uu = z(1,j)**2 + z(2,j)**2
	   zz = zz + uu
	   ee = ee + uu/(fks(N2P1) + fks(j) + F)
	enddo

C the sum in just over 1/2 of the spectrum (is simmetric respect to 0)
C because of the factor 1/2
	ee=ee/NS/NS
	zz=zz/NS/NS
	write(KUPHYS,*) real(t),ee,zz,real(cfl)

	return
	end
