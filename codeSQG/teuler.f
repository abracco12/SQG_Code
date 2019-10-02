c	Teuler.f
c
c	Integrazione secondo Euler sulla parte non lineare ed
c	esatta su quella lineare.
c
c	Z=z(t)
c	Znl=Nlt(t)
c	Znew=z(t+1)
c
c	Dimensioni: (n,n) (codificato)
c
c	Versione NX=NY=N
c
c	Settembre 1994
c
	subroutine teuler(z,znl,znew,dt)
	implicit none

	include "paran.h"
	include "commk.h"

	real*8 z(N,N),znl(N,N),znew(N,N)
	real*8 dt
	integer i,j, ic

	do j = 1, N
	   do i = 1, N-1, 2
	      ic = i/2+1
	      znew(i,j)   = (z(i,j)+dt*znl(i,j))*lre(ic,j) -
     &                      (z(i+1,j)+ dt*znl(i+1,j))*lim(ic,j)
	      znew(i+1,j) = (z(i+1,j)+dt*znl(i+1,j))*lre(ic,j) +
     &                      (z(i,j)+dt*znl(i,j))*lim(ic,j)
	   end do
	end do
        
c        call forcing(znew)
     
	return
	end
c -----------------------------------------------------------

