c To evaluate the spectral streamfunction from the spectral vorticity.
c
      subroutine calcpsi(z,psi)
      implicit none

      include "paran.h"
      include "commk.h"
      include "paraqg.h"

      real*8 z(N,N), psi(N,N)
      real*8 k2
      integer i, j, ic

        do j = 1, N
           do i = 3, N, 2
              ic = i/2 + 1
              k2 = fks(ic) + fks(j) + F
              psi(i,j)  = -z(i,j)/k2
              psi(i+1,j)= -z(i+1,j)/k2
           enddo
        enddo

        do j = 2, N2P1
           psi(1,j) = -z(1,j)/(fks(j) + F)
           psi(2,j) = -z(2,j)/(fks(j) + F)
        enddo

        do j = N2P2, N
           k2 = fks(N2P1) + fks(j) + F
           psi(1,j)= -z(1,j)/k2
           psi(2,j)= -z(2,j)/k2
        enddo

      end
