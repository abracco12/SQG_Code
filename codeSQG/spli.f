c    SPLINE INTERPOLATION
c    (this version can be used on a cray) 
c 
*******************************************************************az
       subroutine spdiag(p,kspli)
       implicit none

       include "paran.h"

       complex*16 p(0:n+1),data(n)
       real*8 blk(20)
       integer k,j,kspli,i
       real*8 h,s,t
  
       k=kspli
       h=mod(k,2)*.5d0
       j=1
       blk(1)=1.d0
 20    if (j.lt.k) then
         s=0.d0
         do i=1,j
            t=blk(i)/dble(j)
            blk(i)=s+t*(dble(i)-h)
            s=t*(dble(j-i)+h)
         enddO
         j=j+1
         blk(j)=s
         goto 20
      endif

      do i=1,n
         data(i)=dcmplx(0.d0,0.d0)
      enddo
      
      do i=1,k
         j=i-k
         if(j.lt.0) j=j+n
         data(j+1)=blk(i)
      enddo
      call dfour1(data,n,-1)

cc
c     coefficients have to be divided by p because of the FFT used in the
c     program
c
cc

      do i=0,n-1
         p(i)=1.d0/data(i+1)/dble(n)
      enddo
      
      p(n)=0.d0
      p(n+1)=0.d0

      return
      end

*******************************************************************az
        subroutine spcoef(a,b)
        implicit none
cc
c     
c     
c     The spectral field is read and the data multiplied for the  
c     spectral coefficients.
c
c     1:   values on the vertical axis (kx=0; b(0,j) with j=0,n-1)
c     2:   Nyquist vertical values (kx=n/2 ; b(n,j) with j=0,n-1)
c     3:   values with positive x and y; and with negative x and 
c          positive y
c     4:   limit nyquist points.
c          The codification of the spectral field set to 0.0
c          (kx,ky)=(n/2,0) e (kx,ky)=(n/2,n/2), but they are in the region
c	   of truncation. 
c
cc

        include "paran.h"
        include "commspli.h"

        complex*16 p0,pn2,b(0:n2,0:n+1)
        real*8 a(0:n-1,0:n-1)
        integer i,j

        p0=p(0)
        pn2=p(n2)

        do j=0,n2                                                 !1
            b(0,j)=dcmplx(a(0,j),-a(1,j))*p(j)*p0
        enddo
        do j=n2+1,n-1
            b(0,j)=dcmplx(a(0,n-j),a(1,n-j))*p(j)*p0
        enddo

        do j=1,n2m1                                               !2
            b(n2,j)=dcmplx(a(0,n-j),a(1,n-j))*pn2*p(j)
        enddo
        do j=n2p1,n-1
           b(n2,j)=dcmplx(a(0,j),-a(1,j))*pn2*p(j)
        enddo

        do i=1,n2m1                                               !3
           do j=0,n-1
              b(i,j)=dcmplx(a(2*i,j),-a(2*i+1,j))*p(i)*p(j)
           enddo
        enddo

                                                                  !4
        b(n2,0) =dcmplx(a(1,0),0.d0)*pn2*p0		!test senza trunc
        b(n2,n2)=dcmplx(a(1,n2),0.d0)*pn2*pn2		!test senza trunc
        

        call fftcr(b)

        return
        end

********************************************************************az

      subroutine spli2l(a,xp,yp,r)
      implicit none

      include "paran.h"
      include "commspli.h"
     
      real*8 a(*),xp(np),yp(np),r(np)
      real*8 b(nwg),xspli(np),yspli(np)
      integer i,j,npt,nt,imod

      imod=mod(kspli,2)

      if (imod.eq.0) then
         do i=1,np
            xspli(i)=xp(i)-1.d0
            yspli(i)=yp(i)-1.d0
         enddo
      elseif (imod.ne.0) then     
         do i=1,np 
            xspli(i)=xp(i)-.5d0
            yspli(i)=yp(i)-.5d0
         enddo
      endif
      
      j=1
      call copyp(b,a)
      npt=np
1     nt=min(llot,npt)
      if(nt.le.0) return 
      call spli2c(b,xspli(j),yspli(j),r(j),nt)
      j=j+nt
      npt=npt-nt
      goto 1
      
      end
      
*******************************************************************az

      subroutine spli2c(a,xp,yp,r,npt1)
      implicit none
 
      include "paran.h"
      include "commspli.h"

      integer npt1
      real*8 a(-kspli1:n+1,-kspli1:n+1)
      integer ix(llot),iy(llot)
      real*8 xp(npt1),yp(npt1),r(npt1)
      real*8 xx,yy,xj,xj1,xl,xlj,yl,ylj
      integer lx,ly,i,j
      complex*16 w(1:llot,-kspli1:0,-kspli1:0)  

      do i=1,npt1
         xx=xp(i)/dble(n)
         xx=(xx-dint(xx))*dble(n)
         yy=yp(i)/dble(n)
         yy=(yy-dint(yy))*dble(n) 
         ix(i)=int(xx)
         xp(i)=xx-dble(ix(i))
         iy(i)=int(yy)
         yp(i)=yy-dble(iy(i))
      enddo

      do lx=-kspli1,0
         do ly=-kspli1,0
            do i=1,npt1
               w(i,lx,ly)=a(ix(i)+lx,iy(i)+ly)
            enddo
         enddo
      enddo

      do j=kspli1,1,-1
         xj=dble(j)
         xj1=1.d0/(xj*xj)
         do lx=0,-j+1,-1
            xl=dble(lx)
            xlj=xl+xj
            do ly=0,-j+1,-1
               yl=dble(ly)
               ylj=yl+xj
               do i=1,npt1
                  w(i,lx,ly)=xj1*((yp(i)-yl)*((xp(i)-xl)*w(i,lx,ly)
     &             +(xlj-xp(i))*w(i,lx-1,ly))+(ylj-yp(i))*((xp(i)-xl)
     &             *w(i,lx,ly-1)+(xlj-xp(i))*w(i,lx-1,ly-1)))
               enddo
            enddo
         enddo
      enddo

      do i=1,npt1
         r(i) = w(i,0,0)
      enddo
      
      return
      end

*******************************************************************az

      subroutine inispli
      implicit none

      include "paran.h"
      include "commspli.h"

      call spdiag(p,kspli)

      return
      end

*******************************************************************az

      subroutine interpo(a,afl,xp,yp)

      include "paran.h"

      real*8 a(n,n)
      real*8 b(0:n+1,0:n+1)
      real*8 afl(np),xp(np),yp(np) 

      call spcoef(a,b)
      call spli2l(b,xp,yp,afl)

      return
      end

*******************************************************************az
      subroutine copyp(a,b)
      implicit none

      include "paran.h"
      include "commspli.h"

      real*8 a(-kspli1:n+1,-kspli1:n+1)
      real*8 b(0:n+1,0:n+1)
      integer I,J

      do i=0,n+1
         do j=0,n+1
            a(i,j)=b(i,j)
         enddo
      enddo
      
      do i=-kspli1,-1
         do j=0,n+1
            a(i,j)=b(i+n,j)
         enddo
      enddo
      
      do j=-kspli1,-1
         do i=-kspli1,n+1
            a(i,j)=a(i,j+n)
         enddo
      enddo
      
      return
      end

*******************************************************************az





