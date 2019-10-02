c	Nlt.f
c
c	Non-linear term : spectral part
c
c	In input: zp1 = z*d1psi = z*v
c		  zp2 = z*d2psi =  z*um (=-z*u)
c
c	in output: znl = - d(psi,z)/d(x,y) = ikx * zp2 - iky * zp1 
c                        
c
c	Dimension: (n,n) - Fourier space (codified)
c
c	Version NX=NY=N
c
c
	subroutine nlt(zp1,zp2,znl)
	implicit none

	include "paran.h"
	include "commk.h"
        include "commfis.h"

	real*8 zp1(N,N),zp2(N,N),znl(N,N)
        real*8 fk2,scra
	integer i,j,ic
 

	do j = 1, N
	   do i = 3, N, 2
	      ic = i/2 + 1
              znl(i,j)   = -fk(ic)*zp2(i+1,j) + fk(j)*zp1(i+1,j) 
              znl(i+1,j) =  fk(ic)*zp2(i,j)   - fk(j)*zp1(i,j)
     
	   end do
	end do

c       fkx=0.0
	do j = 2, N2P1
           znl(1,j) =  fk(j)*zp1(2,j)
           znl(2,j) = -fk(j)*zp1(1,j) 
	end do

c        fkx = fk(N2P1)
	do j = N2P2, N
           znl(1,j) = -fk(N2P1)*zp2(2,j) + fk(j)*zp1(2,j)
           znl(2,j) =  fk(N2P1)*zp2(1,j) - fk(j)*zp1(1,j)
	end do

	return
	end
c ------------------------------------------------------------------
