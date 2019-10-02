c	Tadba.f
c
c	Integrazione Adams-Bashforth sulla parte non lineare ed
c	esatta su quella lineare.
c
c	Schema di Legras
c
c	Z=z(t)
c	Znl=Nlt(t)
c	Znlold=Nlt(t-1)
c       Znlold2=Nlt(t-2)
c	Znew=z(t+1)
c
c	Dimensioni: (n,n) (codificato)
c
c	Versione NX=NY=N
c
c	Settembre 1994
c
c ------------------------------------------------------------
	subroutine tadba2 (z,znl,znlold,znew,dt)
	implicit none

	include "paran.h"
	include "commk.h"

	real*8 z(N,N),znl(N,N),znlold(N,N),znew(N,N)
	integer i,j, ic
	real*8 dt2,dt

        dt2=dt*0.5

        do j = 1, N
	   do i = 1, N-1, 2
	      ic = i/2+1
	      znew(i,j)   = z(i,j)*lre(ic,j) - z(i+1,j)*lim(ic,j)
     &                    + dt2*(3.0*znl(i,j)*lre(ic,j) - 
     &                      znlold(i,j)*(lre(ic,j)*lre(ic,j) -
     &                      lim(ic,j)*lim(ic,j)))
     &                    - dt2*(3.0*znl(i+1,j)*lim(ic,j) - 
     &                      znlold(i+1,j)*2.0*lre(ic,j)*lim(ic,j))
	      znew(i+1,j) = z(i+1,j)*lre(ic,j) + z(i,j)*lim(ic,j)
     &                    + dt2*(3.0*znl(i+1,j)*lre(ic,j) - 
     &                      znlold(i+1,j)*(lre(ic,j)*lre(ic,j) -
     &                      lim(ic,j)*lim(ic,j)))
     &                    + dt2*(3.0*znl(i,j)*lim(ic,j) - 
     &                      znlold(i,j)*2.0*lre(ic,j)*lim(ic,j))
	   end do
        end do
c        call forcing(znew)
	return
	end
c ------------------------------------------------------------
c
c       Implemented by MRS
c
c ------------------------------------------------------------
c
	subroutine tadba3 (z,znl,znlold,znlold2,znew,dt)
	implicit none

	include "paran.h"
	include "commk.h"

	real*8  z(N,N)
	real*8  znl(N,N),znlold(N,N),znlold2(N,N)
	real*8  znew(N,N)
	integer i,j,ic
	real*8  c, dt

        c = dt/12.0

        do j = 1, N
	   do i = 1, N-1, 2
	      ic = i/2+1
	      znew(i,j) = z(i,j)*lre(ic,j) - z(i+1,j)*lim(ic,j)
     &                  + c*(23.0*znl(i,j)*lre(ic,j) - 
     &       16.0*znlold(i,j)*(lre(ic,j)*lre(ic,j)-lim(ic,j)*lim(ic,j))+
     &	      5.0*znlold2(i,j)*lre(ic,j)*
     &                    (lre(ic,j)*lre(ic,j)-3.0*lim(ic,j)*lim(ic,j)))
     &                  - c*(23.0*znl(i+1,j)*lim(ic,j) - 
     &        16.0*znlold(i+1,j)*2.0*lre(ic,j)*lim(ic,j) +
     &	      5.0*znlold2(i+1,j)*lim(ic,j)*
     &                    (3.0*lre(ic,j)*lre(ic,j)-lim(ic,j)*lim(ic,j)))
	      znew(i+1,j) = z(i+1,j)*lre(ic,j) + z(i,j)*lim(ic,j)
     &                  + c*(23.0*znl(i+1,j)*lre(ic,j) - 
     &     16.0*znlold(i+1,j)*(lre(ic,j)*lre(ic,j)-lim(ic,j)*lim(ic,j))+
     &	      5.0*znlold2(i+1,j)*lre(ic,j)*
     &                    (lre(ic,j)*lre(ic,j)-3.0*lim(ic,j)*lim(ic,j)))
     &                  + c*(23.0*znl(i,j)*lim(ic,j) - 
     &        16.0*znlold(i,j)*2.0*lre(ic,j)*lim(ic,j) +
     &	      5.0*znlold2(i,j)*lim(ic,j)*
     &                    (3.0*lre(ic,j)*lre(ic,j)-lim(ic,j)*lim(ic,j)))
	   end do
        end do

c        call forcing(znew)

	return
	end
c ------------------------------------------------------------


