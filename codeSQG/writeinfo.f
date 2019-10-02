c	writeinfo.f
c
c	Writes the parameters for a integration in 2D-QG turbulence
c
c       Programmed by MRS - Gennaio 1997
c
	subroutine writeinfo

	implicit none

	include "paran.h"
	include "paraio.h"
	include "commsim.h"
        include "commtrunc.h"
	include "commfis.h"
	include "commk.h"
        include "paraqg.h"                            
	include "beta.h"			

	real*8 B, A, C
	integer log	 ! Logical no. for info. file.
	character	progid*32, runid*6, comment*70


c....	Write information about the parameters in the calculation.
	log = 21
	B = beta
	progid = 'SQG TURBULENCE on the beta-plane '
	
	runid  = 'sqgbeta01'
	comment = 'SQG run # 1 
     &             '
	open(unit=log,file='parameters.info',
     &        form='formatted',status='unknown')
	write (log,68) progid
	write (log,65) ' ' 
	write (log,66) comment
	write (log,65) ' '
	write (log,69) 'Run ID: ',runid
	write (log,65) ' '
        write (log,64) ' Resolution of field ',n,'x',n
        write (log,65) ' '
	write (log,61) 'Integration parameters: '
	write (log,63) 'Domain size x                             = ',xlx
	write (log,63) 'Domain size y                             = ',xly
	write (log,63) 'Viscosity (UV) coefficient                = ',xnu
	write (log,63) 'Viscosity (UV) order                      = ',alpu
	write (log,63) 'Viscosity (IR) coefficient                = ',xni
	write (log,63) 'Viscosity (IR) order                      = ',alpi
	write (log,63) 'Time step                                 = ',dt
	write (log,61) 'Maximum time steps                        = ',npas
	write (log,61) 'Steps betweens shapshots of field         = ',iout
	write (log,62) 'Truncation used ?                      ',truncate
	write (log,63) 'Truncation ratio                       = ',rtrunc
        write (log,63) 'Rossby parameter (gamma^2 = F)            = ',F
        write (log,63) 'Coriolis parameter (beta)                 = ',B
	close(log)

c ----------------------------------------------------------------------

  61	format(A44,I8)
  62	format(A44,L8)
  63    format(A44,E16.8)
  64    format(A21,I4,A1,I4)
  65    format(A1)
  66	format(A70)
  68	format(A32)
  69	format(A8,A6)

	return
	end
c ----------------------------------------------------------------------
