c	fft_inv.f
c
c	Invers FFt for a real vector of dimension (n,n)
c
c	Version NX=NY=N
c
c	Version codification/decodification real vector 
c
c       in input  z(n,n): real vector with positive components
c                         codified in the Fourier space
c       in output z(n,n): real vector (real field)
c
c       the two components (K_x= +/- N/2; K_y= 0) and
c       (K_x= +/- N/2; K_y= +/- N/2) are = 0.0
c
	subroutine fft_inv(z)
	implicit none

	include "paran.h"
	include "commfft.h"

	real*8 z(N,N)
	complex*16 wa(N,N)
	complex*16 IM
	real*8 C
	parameter (IM=(0.d0,1.d0))
	integer i,j, ia, ic, jc
	parameter (C=1.d0/(N*N))

C set corners
	j = N2P1
	wa(1,1) = z(1,1) + IM*z(2,1)
	wa(1,j) = z(1,j) + IM*z(2,j)
	wa(j,1) = (0.d0,0.d0)
	wa(j,j) = (0.d0,0.d0)

c	i = 1
	do j = 2, N2
	   jc = N - j + 2
	   wa(1,j)   =   z(1,j) + IM*z(2,j)
	   wa(1,jc)  =   z(1,j) - IM*z(2,j)
	end do

c	j = 1
	do i = 3, N, 2
	   ic = i/2 + 1
	   wa(ic,1) =   z(i,1) + IM*z(i+1,1)
	   ic = N - ic + 2
	   wa(ic,1) =   z(i,1) - IM*z(i+1,1)
	end do

c       i = 1
	ic = N2P1
	do j = 2, N2
	   jc = N - j + 2
	   wa(ic,j)  = z(1,jc) - IM*z(2,jc)
	   wa(ic,jc) = z(1,jc) + IM*z(2,jc)
	end do

	j = N2P1
	do i = 3, N, 2
	   ic = i/2 + 1
	   wa(ic,j) =   z(i,j) + IM*z(i+1,j)
	   ic = N - ic + 2
	   wa(ic,j) =   z(i,j) - IM*z(i+1,j)
	end do

        do j = 2, N2
	   jc = N - j + 2
	   do i = 3, N-1, 2
	      ic = i/2 + 1
	      wa(ic,j)  = z(i,j) + IM*z(i+1,j)
	      ic = N - ic + 2
	      wa(ic,jc) = z(i,j) - IM*z(i+1,j)
	   end do
	end do

	do j = N2P2, N
	   jc = N - j + 2
	   do i = 3, N-1, 2
	      ic = i/2 + 1
	      wa(ic,j)  = z(i,j) + IM*z(i+1,j)
	      ic = N - ic + 2
	      wa(ic,jc) = z(i,j) - IM*z(i+1,j)
	   end do
	end do

c   Inverse FFT	: w = exp(-isign*sqrt(-1)*2pi)
c	w = exp(-sqrt(-1)*2pi)

	call four2d(wa,1)
	do j=1,N
	   do i=1,N
	      z(i,j)=C*wa(i,j)
	   end do
	end do
	
	return
	end
c --------------------------------------------------------------

