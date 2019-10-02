c	Inieuler.f
c
c	Initialization of z, zcur, zold e znl.
c	Dimension: z(ns)
c	
c	Version NX=NY=N
c
c
	subroutine inieuler(zold)
	implicit none

	include "paran.h"
	include "paraio.h"
	include "commsim.h"
	include "commtrunc.h"
	include "commfis.h"

        real*8 zold(NS)


c       Initialization according with the file: fort.9
        read(KUPAR,*)xlx,xly,xnu,alpu,xni,alpi,dt
        read(KUPAR,*)npas,iout,ioutl
        read(KUPAR,*)rtrunc,truncate
        read(KUPAR,*)firststep

	dk=pi2/xlx

c	Initialization of wave-numbers
	call defkbeta

c	Read the vorticity field real -> spectral
	call lire(KUIN,zold,NS)
	call fft_dir(zold)

c	Circular truncation (Nyquist included)
	if (truncate) then
	 call trunc(zold)
	end if

c	Subtract the mean
	zold(1)=0.0d0
	zold(2)=0.0d0
        print*, zold(4), zold(5)
        	return
	end
c ----------------------------------------------------------------	


