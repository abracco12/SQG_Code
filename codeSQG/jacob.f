c	Jacob.f
c
c
c       Non-linear term  J=D(psi,z)/D(x,y)
c
c       In input: z
c
c       in output: znl = - d(psi,z)/d(x,y)
c                  v (=v) e um(=-u)
c       Dimension: (n,n) - Fourier space (codified)
c
c       Version NX=NY=N
c
	subroutine jacob(z,znl,ipas,v,um)
	implicit none

	include "paran.h"
        include "commsim.h"                                     
	include "paraio.h"                                     

        integer ipas,j,i
	real*8 z(N,N),znl(N,N),v(NS),um(NS)
	real*8 d1psi(NS),d2psi(NS)

c	 d1psi=ikx*z/|(-k)|=v
c	 d2psi=iky*z/|(-k)|=um
	call deriv(z,d1psi,d2psi)

	   
c       From spectral to real space
        call fft_inv(z)
        call fft_inv(d1psi)
        call fft_inv(d2psi)


	call copy (d1psi,v,NS)
	call copy (d2psi,um,NS)

c	Check whether stability is obtained.
	if (mod(ipas,check).eq.0) call courant(d2psi,d1psi)

c       Products d1psi*z e d2psi*z
        call prods(z,d1psi,d2psi,NS)

c       From real to spectral space
        call fft_dir(d1psi)
        call fft_dir(d2psi)

c       Znl= - Jacobiano
	call nlt(d1psi,d2psi,znl)

	return
	end
