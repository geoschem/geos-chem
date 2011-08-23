! $Id: EFOLD.f,v 1.1 2009/09/16 14:06:48 bmy Exp $
      subroutine EFOLD (F0, F1, N, F)
!-----------------------------------------------------------------------
!---  calculate the e-fold between two boundaries, given the value
!---     at both boundaries F0(x=0) = top, F1(x=1) = bottom.
!---  presume that F(x) proportional to exp[-A*x] for x=0 to x=1
!---          d2F/dx2 = A*A*F  and thus expect F1 = F0 * exp[-A]
!---           alternatively, could define A = ln[F0/F1]
!---  let X = A*x, d2F/dX2 = F
!---  assume equal spacing (not necessary, but makes this easier)
!---      with N-1 intermediate points (and N layers of thickness dX = A/N)
!---
!---  2nd-order finite difference:  (F(i-1) - 2F(i) + F(i+1)) / dX*dX = F(i)
!---      let D = 1 / dX*dX:
!
!  1  |   1        0        0        0        0        0   |    | F0 |
!     |                                                    |    | 0  |
!  2  |  -D      2D+1      -D        0        0        0   |    | 0  |
!     |                                                    |    | 0  |
!  3  |   0       -D      2D+1      -D        0        0   |    | 0  |
!     |                                                    |    | 0  |
!     |   0        0       -D      2D+1      -D        0   |    | 0  |
!     |                                                    |    | 0  |
!  N  |   0        0        0       -D      2D+1      -D   |    | 0  |
!     |                                                    |    | 0  |
! N+1 |   0        0        0        0        0        1   |    | F1 |
!      
!-----------------------------------------------------------------------
!  Advantage of scheme over simple attenuation factor: conserves total
!  number of photons - very useful when using scheme for heating rates.
!  Disadvantage: although reproduces e-folds very well for small flux
!  differences, starts to drift off when many orders of magnitude are
!  involved.
!-----------------------------------------------------------------------
      implicit none
      real*8 F0,F1,F(250)  !F(N+1)
      integer N
      integer I
      real*8 A,DX,D,DSQ,DDP1, B(101),R(101)
!
      if(F0.eq.0.d0) then
        do I=1,N
          F(I)=0.d0
        enddo
        return
      elseif(F1.eq.0.d0) then
        A = DLOG(F0/1.d-250)
      else
        A = DLOG(F0/F1)
      endif
!
      DX = float(N)/A
      D = DX*DX
      DSQ = D*D
      DDP1 = D+D+1.d0
!
      B(2) = DDP1
      R(2) = +D*F0
      do I=3,N
        B(I) = DDP1 - DSQ/B(I-1)
        R(I) = +D*R(I-1)/B(I-1)
      enddo
      F(N+1) = F1
      do I=N,2,-1
        F(I) = (R(I) + D*F(I+1))/B(I)
      enddo
      F(1) = F0
      return
      end
