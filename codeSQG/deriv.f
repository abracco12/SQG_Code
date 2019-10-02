c	Deriv.f
c
c	To evalue (u,-v) from  z in the Fourier space. 
c
c	 d1psi=ikx*z/sqrt(-k**2)=v
c	 d2psi=iky*z/sqrt(-k**2)=-u
c
c	Dimension: (n,n) -  Fourier space
c
c	Version NX=NY=N
c
	subroutine deriv(z,d1psi,d2psi)
	implicit none

	include "paran.h"
	include "commk.h"
        include "paraqg.h"

        
	real*8 z(N,N),d1psi(N,N),d2psi(N,N)
	integer i, j, ic
	real*8 k2, kx, ky

c       Put the mean = 0.0
        d1psi(1,1)=0.0d0
        d1psi(2,1)=0.0d0
        d2psi(1,1)=0.0d0
        d2psi(2,1)=0.0d0

	do j = 1, N
	   do i = 3, N, 2
	      ic = i/2 + 1
	      k2 = sqrt(fks(ic) + fks(j) + F)
	      kx = fk(ic)/k2
	      ky = fk(j)/k2
	      d2psi(i,j)   =   ky*z(i+1,j)                 !Re(-u)
	      d2psi(i+1,j) =  -ky*z(i,j)                   !Im(-u)
	      d1psi(i,j)   =   kx*z(i+1,j)                 !Re(v)
	      d1psi(i+1,j) =  -kx*z(i,j)                   !Im(v)
	   enddo
	enddo

	do j = 2, N2P1
	   ky = fk(j)/sqrt(fks(j) + F)
	   d2psi(1,j) =   ky*z(2,j)
	   d2psi(2,j) =  -ky*z(1,j)
	   d1psi(1,j) =   0.d0
	   d1psi(2,j) =   0.d0
	enddo

	do j = N2P2, N
	   k2 = sqrt(fks(N2P1) + fks(j) + F)
	   kx = fk(N2P1)/k2
	   ky = fk(j)/k2
	   d2psi(1,j) =   ky*z(2,j)
	   d2psi(2,j) =  -ky*z(1,j)
	   d1psi(1,j) =   kx*z(2,j)
	   d1psi(2,j) =  -kx*z(1,j)
	enddo

	end
c ----------------------------------------------------------------
