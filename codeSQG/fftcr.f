*******************************************************************
*
*	Routine: Fast Fourier Transform
*
*	It has benn used in the spline interpolation 
*

      subroutine fftcr(z)

      include "paran.h"
      include "commfft.h"

      complex*16 wa(n,n)
      complex*16 z(n2+1,n+2)
      integer i,j
    
      do i=1,n2+1
         do j=1,n
            wa(i,j)=z(i,j)
         enddo
      enddo
  
      do i=n2+2,n
         do j=2,n
            wa(i,j)=conjg(wa(n-i+2,n-j+2))
         enddo
      enddo

      do i=n2+2,n
         wa(i,1)=conjg(wa(n-i+2,1))
      enddo

      do i=1,n
         do j=1,n
            wa(i,j)=conjg(wa(i,j))
         enddo
      enddo

      call four2d(wa,1)
      call transf(wa,z)
 

      return
      end

*******************************************************************az    
      subroutine transf(wa,z)

      include "paran.h"

      integer i,j
      real*8 z(n+2,n+2),wa(2,n,n)

      do i=1,n
         do j=1,n
            z(i,j)=wa(1,i,j) 
         enddo
      enddo

      end

*******************************************************************az

