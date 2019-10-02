c	Paraio.h
c
c	Input/Output files.
c
	integer KUPAR,KUIN,KUINPASS,KUINPASS2,KUINLAGR
	character*12 KUPLOT,KUPSI,KUPPASS,KUPPASS2,KUPLAGR,KUPLAGR2
	parameter(KUPAR=9,KUIN=10,KUINPASS=12,
     >            KUINLAGR=16)
	parameter(KUPLOT='vort.unf')
	parameter(KUPSI= 'psi.unf')
	parameter(KUPPASS='scalar1.unf')
	parameter(KUPLAGR='lagrange.unf')
	parameter(KUPLAGR2='psilagr.unf')
	integer KUPHYS,KUPSCA,KUPSCA2
	parameter(KUPHYS=20,KUPSCA=22)
c ----------------------------------------------------------------
