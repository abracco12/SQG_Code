c	fft_dir.f
c
c	direct fft for a real vector of dimension (n,n)
c
c	Version NX=NY=N
c
c	codification/decodification of a real vector
c
c       in input  z(N,N): real vector
c       in output z(N,N): real vector with positive component
c                         (almost) in the Fourier space
c                         (odd real, even: imaginary)
c     
	subroutine fft_dir(z)
	implicit none

	include "paran.h"
	include "commfft.h"

	integer i,j, ia
	real*8 z(N,N)
	real*8 wa(2,N,N)

	do j = 1, N
	   do i = 1, N
	      wa(1,i,j) = z(i,j)
	      wa(2,i,j) = 0.d0 
	   end do
	end do

c       Direct FFT
c       w = exp(-isign*sqrt(-1)*2pi) = exp(+sqrt(-1)*2pi)

	call four2d(wa,-1)

	do j = 1, N
	   do i = 1, N-1, 2
	      ia = i/2 + 1
	      z(i,j)   = wa(1,ia,j)
	      z(i+1,j) = wa(2,ia,j)
	   end do
	end do

C Move the components K_x= +/- N/2 e K_y = -1, ..., -(N/2-1)

	i = N2P1
	do j = N2P2, N
	   z(1,j) = wa(1,i,j)
	   z(2,j) = wa(2,i,j)
	end do

	return
	end
c --------------------------------------------------------------
