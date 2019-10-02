c Parametri per il modello N,P,Z per l'ecosistema (ref: capitolo Proceedings
c International Symposium Shallow Flows 2003)

	real*8 s0,Vmax,kn,c0,g,gamma,eps,mup,muz,mud,Amp
        parameter (s0=0.00125d0,Vmax=0.66d0,kn=0.5d0)
        parameter (c0=8.d0,g=2.d0,gamma=0.75d0,eps=1.d0)
	parameter (mup=0.03,muz=0.2,mud=0,Amp=100.)

c Choose the type of simulation: either fixed or relaxation flux,
c either Center (one active blob that covers 12% of the area)
c or Patches (several patches that cover the same area)
c or Vortex (Patches correlated with vorticity field: Choose the 
c proper values for OkuboWeiss thresholds OWmin and OWmax). 

        character*30 file0
        integer*4 NF
        real*8 OWmin,OWmax
        parameter (file0="pp_lang_fixed_NF1.dat")
        parameter (OWmin=-32.,OWmax=32.) !Only for VORTEX (K=10)
c        parameter (OWmin=-2.7d5,OWmax=6.d5) !Only for VORTEX (K=40)
        parameter (NF=4) !radius of patches is Rblob/NF, Rblob=100Km

c#define FIXED_FLUX
#define RELAX_FLUX
#define PATCHES
c#define VORTEX   !Remember to choose the Okubo-Weiss thresholds

c Choose the advecting field 

#define TURBO
c#define LANGEVIN

