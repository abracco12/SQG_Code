c     subroutine four2d
c
c     the cooley-tukey fast fourier transform in usasi basic fortran
c
      subroutine four2d(data,isign)
      implicit none

      include "paran.h"
      include "commfft.h"

      real*8 data(2*NS)
      integer isign

      integer ntot, idim
      integer n1, np1, np1tw, np2, np2hf
      integer j, i, i1, i1max, i2, i3, j3
      integer m, mmax, ipar
      integer k1, k2, k3, k4, lmax, l
      integer kmin, kdif, kstep

      real*8 twopi
      real*8 tempr, tempi, theta
      real*8 wr, wi, wstpr, wstpi, w2r, w2i, w3r, w3i
      real*8 u1r, u1i, u2r, u2i, u3r, u3i, u4r, u4i
      real*8 t2r, t2i, t3r, t3i, t4r, t4i

      twopi=pi2
      ntot=2
      do 2 idim=1,NDIM
      if(NN(idim))700,700,2
 2    ntot=ntot*NN(idim)
c
c     main loop for each dimension
c
      np1=2
      do 600 idim=1,NDIM
      n1=NN(idim)
      np2=np1*n1
      if(n1-1)700,600,100
c
c     shuffle data by bit reversal, since n=2**k. as the shuffling
c     can be done by simple interchange, no working array is needed
c
100   np2hf=np2/2
      j=1
      do 160 i2=1,np2,np1
      if(j-i2)110,130,130
110   i1max=i2+np1-2
      do 120 i1=i2,i1max,2
      do 120 i3=i1,ntot,np2
      j3=j+i3-i2
      tempr=data(i3)
      tempi=data(i3+1)
      data(i3)=data(j3)
      data(i3+1)=data(j3+1)
      data(j3)=tempr
120   data(j3+1)=tempi
130   m=np2hf
140   if(j-m)160,160,150
150   j=j-m
      m=m/2
      if(m-np1)160,140,140
160   j=j+m
c
c     main loop. perform fourier trasforms of lenght four, with one
c     of lenght two if needed. the twiddle factor w=exp(-isign*2*pi*
c     sqrt(-1)) and repeat for w=exp(-isign*sqrt(-1))*conjugate(w).
c
      np1tw=np1+np1
      ipar=n1
310   if(ipar-2)350,330,320
320   ipar=ipar/4
      goto 310
330   do 340 i1=1,np1,2
      do 340 k1=i1,ntot,np1tw
      k2=k1+np1
      tempr=data(k2)
      tempi=data(k2+1)
      data(k2)=data(k1)-tempr
      data(k2+1)=data(k1+1)-tempi
      data(k1)=data(k1)+tempr
340   data(k1+1)=data(k1+1)+tempi
350   mmax=np1
360   if(mmax-np2hf)370,600,600
370   lmax=max0(np1tw,mmax/2)
      if(mmax-np1)405,405,380
380   theta=-twopi*dble(np1)/dble(4*mmax)
      if(isign)400,390,390
390   theta=-theta
400   wr=cos( theta )
      wi=sin( theta )
      wstpr=-2.d0*wi*wi
      wstpi=2.d0*wr*wi
405   do 570 l=np1,lmax,np1tw
      m=l
      if(mmax-np1)420,420,410
410   w2r=wr*wr-wi*wi
      w2i=2.d0*wr*wi
      w3r=w2r*wr-w2i*wi
      w3i=w2r*wi+w2i*wr
420   do 530 i1=1,np1,2
      kmin=ipar*m+i1
      if(mmax-np1)430,430,440
430   kmin=i1
440   kdif=ipar*mmax
450   kstep=4*kdif
      do 520 k1=kmin,ntot,kstep
      k2=k1+kdif
      k3=k2+kdif
      k4=k3+kdif
      if(mmax-np1)460,460,480
460   u1r=data(k1)+data(k2)
      u1i=data(k1+1)+data(k2+1)
      u2r=data(k3)+data(k4)
      u2i=data(k3+1)+data(k4+1)
      u3r=data(k1)-data(k2)
      u3i=data(k1+1)-data(k2+1)
      if(isign)470,475,475
470   u4r=data(k3+1)-data(k4+1)
      u4i=data(k4)-data(k3)
      goto 510
475   u4r=data(k4+1)-data(k3+1)
      u4i=data(k3)-data(k4)
      goto 510
480   t2r=w2r*data(k2)-w2i*data(k2+1)
      t2i=w2r*data(k2+1)+w2i*data(k2)
      t3r=wr*data(k3)-wi*data(k3+1)
      t3i=wr*data(k3+1)+wi*data(k3)
      t4r=w3r*data(k4)-w3i*data(k4+1)
      t4i=w3r*data(k4+1)+w3i*data(k4)
      u1r=data(k1)+t2r
      u1i=data(k1+1)+t2i
      u2r=t3r+t4r
      u2i=t3i+t4i
      u3r=data(k1)-t2r
      u3i=data(k1+1)-t2i
      if(isign)490,500,500
490   u4r=t3i-t4i
      u4i=t4r-t3r
      goto 510
500   u4r=t4i-t3i
      u4i=t3r-t4r
510   data(k1)=u1r+u2r
      data(k1+1)=u1i+u2i
      data(k2)=u3r+u4r
      data(k2+1)=u3i+u4i
      data(k3)=u1r-u2r
      data(k3+1)=u1i-u2i
      data(k4)=u3r-u4r
520   data(k4+1)=u3i-u4i
      kmin=4*(kmin-i1)+i1
      kdif=kstep
      if(kdif-np2)450,530,530
530   continue
      m=mmax-m
      if(isign)540,550,550
540   tempr=wr
      wr=-wi
      wi=-tempr
      goto 560
550   tempr=wr
      wr=wi
      wi=tempr
560   if(m-lmax)565,565,410
565   tempr=wr
      wr=wr*wstpr-wi*wstpi+wr
570   wi=wi*wstpr+tempr*wstpi+wi
      ipar=3-ipar
      mmax=mmax+mmax
      goto 360
c
c     end of loop over each dimension
c
600   np1=np2
700   return
      end
