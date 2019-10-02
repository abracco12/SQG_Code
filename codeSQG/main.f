c	Main4.f
c
c	Principal programme for integrating  SQG eqs. (year 1999)
c       Modified from a QG code
c
c	Versione NX=NY=N
c
c
	implicit none

	include "paran.h"
	include "paraio.h"
	include "commsim.h"
        include "commtrunc.h"
	include "commfis.h"
	include "commk.h"
        include "paraqg.h"                            
	include "beta.h"			
	include "commfrozen.h" 

c....	Large data structures.
	real*8 z(NS)	 ! The spectral vorticity field
	real*8 z1(NS)    ! The spectral vorticity field at t1=dt/2
	real*8 psi(NS)	 ! The streamfunction field (diagnostic)
	real*8 zcur(NS)  ! The physical vorticity field
	real*8 znl0(NS)  ! The nonl. part of the time derivative
	real*8 znl1(NS)  ! The nonl. part of the previous time deriv.
	real*8 znl2(NS)  ! The nonl. part of the preprevious time deriv.
	real*8 um(NS)     ! The -x component of velocity field
	real*8 v(NS)     ! The v component of velocity field

c....	Time variables.
	real*8 time, dt1
        integer ipas,i     ! Step counter

c....   Get the initial field (spectral)
	call inieuler(z)
c....   Computes and writes the streamfunction of the initial field
c	call calcpsi(z,psi)
c	call fft_inv(psi)
c        call ecrire(KUPSI,psi,NS)


c....	Write information about the parameters in the calculation.
	call writeinfo

c....	Initialize the time variables.
	time=0.d0
	ipas=0

c....	Find the invariants.
	call test(z,time)
c	First steps <<<<<<<<<<<<<<<<<<<<<<<<<<<

c	0: no evolution
	if (firststep.eq.0) then
	   stop
	end if

c 	Spectral: z(t)
	dt1=dt/2.
c	1: leap-frog step
	if (firststep.ge.1) then
c       Do half a step by forward Euler
	   call copy(z,zcur,NS)
	   call jacob(zcur,znl1,ipas,v,um)
	   if(nofrozen) then 
	    if (truncate) call trunc(znl1)
            call copy(znl1,znl2,NS)
            call teuler(z,znl1,zcur,dt1)
	    call copy(zcur,z1,NS)
	   end if

c       Do one step by Leap-frog without forcing
	   if(nofrozen) then   
	    call jacob(zcur,znl0,ipas,v,um)
	    if (truncate) call trunc(znl0)
	    call teuler(z,znl0,zcur,dt)
	   end if

	   ipas=ipas+1
           time=time+dt
        end if

c 	Spectral: z(dt),zcur(dt),znl0(dt/2),znl1(0),znl2(0)
c	2: 2nd. order Adam-Bashforth
	if (firststep.ge.1) then

c	Spectral: z(t+1),zcur(t+1),znl(t),znl0(t)
	   if (nofrozen) then
	    call copy(zcur,z,NS)
	    call jacob(zcur,znl0,ipas,v,um)
	    if (truncate) call trunc(znl0)
            call tadba2 (z,znl0,znl1,zcur,dt)
            call copy(znl0,znl1,NS)
	   end if
	   ipas=ipas+1
	   time=time+dt
	end if

c 	Spectral: z(2dt),zcur(2dt),znl0(dt),znl1(dt),znl2(0)

c	3: 3rd. order Adam-Bashforth
	if (firststep.eq.2) then
	   stop
	end if

c        call test(z,time)

c	>>>>>>>>>>  Time evolution <<<<<<<<<<<<<<<<<<<<<<<<

c 	Spectral: z(t),zcur(t),znl0(t-1),znl1(t-1),znl2(t-2)

 1000	if(nofrozen) then
       	  call copy(zcur,z,NS)
	  call jacob(zcur,znl0,ipas,v,um)
 	  if (truncate) call trunc(znl0)
	end if

c 	Spectral: z(t),znl0(t),znl1(t-1),znl2(t-2)
c	Real    : zcur(t)


c	Write the output (in real space) 
	if (mod(ipas,iout).eq.0) then
	 if (nofrozen) call ecrire(KUPLOT,zcur,NS)
	end if


	if (mod(ipas,check).eq.0) then
	  if (nofrozen) call test(z,time)
	end if

c	Time integration part

	if (nofrozen) then
	 call tadba3 (z,znl0,znl1,znl2,zcur,dt)
         call copy(znl1,znl2,NS)
         call copy(znl0,znl1,NS)
	end if
	
c	Spectral: z(t),zcur(t),znl(t-1),znl(t-1),znl(t-2) 

	if (mod(ipas,100).eq.0) then
            write(6,*)' Step just concluded: ',ipas
        endif

	ipas=ipas+1
	time=time+dt

c	>>>>>>>>>>  End time evolution <<<<<<<<<<<<<<<<<<<<<<<<
c	  write (*,*) 'Step',ipas,' passed'
	if(ipas.le.npas) goto 1000

c ----------------------------------------------------------------------
  61	format(A40,I4)
  62	format(A40,L4)
  63    format(A40,E16.8)
  64    format(A21,I4,A1,I4)
  65    format(A1)
  66	format(A70)
c  67	format(A40,C32.8)
  68	format(A32)
  69	format(A6)

	end
c ----------------------------------------------------------------------
