cc
c  commspli.h
c
c   KSPLI SPLINE ORDER ( polynomium of order KSPLI-1)
c   LLOT number of points treated at the same time
c
cc

      integer kspli,kspli1,llot,nwg
      complex*16 p(0:N+1)
      parameter (kspli=4,kspli1=kspli-1)
      parameter (llot=4096)
      parameter (nwg=(n+2+Kspli1)*(n+2+kspli1))
      common /spli/ p
*******************************************************************az
