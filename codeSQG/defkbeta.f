c	Defkbeta.f
c
c	Version nx=ny=n
c
c	Initialization of spectral vectors (codified)
c
c       The beta effect is included in this calculation.
c
c
	subroutine defkbeta
	implicit none

	include "paran.h"
	include "commtrunc.h"
	include "commfis.h"
	include "paraqg.h"
	include "commk.h"
	include "beta.h"

	real*8 fk2
	real*8 fkx
	real*8 scra
	real*8 aniso
	integer i,j

c
c	Fk
	do i=1,N2P1
	 fk(i)=dble(i-1)*dk
	enddo
	do i=N2P2,N
	 fk(i)=dble(i-N-1)*dk
	enddo

c	Fks = fk**2
	do i=1,N
	 fks(i)=fk(i)*fk(i)
	enddo

c
c       Radius for the truncation (fksmax).
c       if  rtrunc=1.0 the Nyquist is included in the truncation.
        rtrunc=min(1.d0,rtrunc)
        if (.not.truncate) rtrunc=1.0
c       New definition for the truncation radius
        fksmax=(fk(n2p1)*rtrunc)**2

c
c       Integration term
c
c       xnu: ultraviolet viscosity
c       xni: infrared viscosity
c       beta: strength of coriolis force
c
c       Re{Linear operator} =
c       exp[-dt*(xnu*(k**2)**alpu + xni*(k**2)**alpi) * (k**2)/(k**2+F)] *
c       cos(dt*beta*kx/(k**2+F))
c
c       Im{Linear operator} =
c       exp[-dt*(xnu*(k**2)**alpu + xni*(k**2)**alpi) * (k**2)/(k**2+F)] *
c       sin(dt*beta*kx/(k**2+F))
c
        lre(1,1) = 0.0d0
        lim(1,1) = 0.0d0

        do j = 1, N
           do i = 3, N, 2
              fk2 = fks(i/2+1) + fks(j)
	      fkx = fk(i/2+1)
              scra = (xnu*fk2**alpu + xni*fk2**alpi)*fk2/sqrt(fk2+F)
	      aniso = beta*fkx/sqrt(fk2+F)
              lre(i/2+1,j) = dexp(-dt*scra)*dcos(dt*aniso)
              lim(i/2+1,j) = dexp(-dt*scra)*dsin(dt*aniso)
           enddo
        enddo

        fkx = 0.0d0
	do j = 2, N2P1
          scra=(xnu*fks(j)**alpu+xni*fks(j)**alpi)*fks(j)/sqrt(fks(j)+F)
	   aniso = beta*fkx/sqrt(fks(j)+F)
           lre(1,j) = dexp(-dt*scra)*dcos(dt*aniso)
           lim(1,j) = dexp(-dt*scra)*dsin(dt*aniso)
        enddo

	fkx = fk(N2P1)
        do j = N2P2, N
           fk2 = fks(N2P1) + fks(j)
           scra = (xnu*fk2**alpu + xni*fk2**alpi)*fk2/sqrt(fk2+F)
	   aniso = beta*fkx/sqrt(fk2)
           lre(1,j) = dexp(-dt*scra)*dcos(dt*aniso)
           lim(1,j) = dexp(-dt*scra)*dsin(dt*aniso)
        enddo
	return
	end
