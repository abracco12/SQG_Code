c	Trunc.f
c
c       Everything is outside of the circunference with radius = fksmax
c	is truncated.
c
c	Dimension: (N,N) (real)
c
c	Version NX=NY=N
c
c
	subroutine trunc(z)
	implicit none

	include "paran.h"
	include "commk.h"
	include "commtrunc.h"
	include "commfis.h"

	real*8 z(N,N)
	real*8 k2
	integer i,j

	do j = 1, N
	   do i = 3, N, 2
	      k2 = fks(i/2+1) + fks(j)
	      if(k2 .gt. fksmax) then
		 z(i,j)   = 0.d0
		 z(i+1,j) = 0.d0
	      end if
	   end do
	end do

	do j = 2, N2P1
	   if(fks(j) .gt. fksmax) then
	      z(1,j) = 0.d0
	      z(2,j) = 0.d0
	   end if
	end do

	do j = N2P2, N
	   k2 = fks(N2P1)+fks(j)
	   if(k2 .gt. fksmax) then
	      z(1,j) = 0.d0
	      z(2,j) = 0.d0
	   end if
	end do

	return
	end
c --------------------------------------------------------------	
