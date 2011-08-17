C $Id: EFOLD.f,v 1.1 2009/09/16 14:06:48 bmy Exp $
      subroutine EFOLD (F0, F1, N, F)
C-----------------------------------------------------------------------
C---  calculate the e-fold between two boundaries, given the value
C---     at both boundaries F0(x=0) = top, F1(x=1) = bottom.
C---  presume that F(x) proportional to exp[-A*x] for x=0 to x=1
C---          d2F/dx2 = A*A*F  and thus expect F1 = F0 * exp[-A]
C---           alternatively, could define A = ln[F0/F1]
C---  let X = A*x, d2F/dX2 = F
C---  assume equal spacing (not necessary, but makes this easier)
C---      with N-1 intermediate points (and N layers of thickness dX = A/N)
C---
C---  2nd-order finite difference:  (F(i-1) - 2F(i) + F(i+1)) / dX*dX = F(i)
C---      let D = 1 / dX*dX:
C
C  1  |   1        0        0        0        0        0   |    | F0 |
C     |                                                    |    | 0  |
C  2  |  -D      2D+1      -D        0        0        0   |    | 0  |
C     |                                                    |    | 0  |
C  3  |   0       -D      2D+1      -D        0        0   |    | 0  |
C     |                                                    |    | 0  |
C     |   0        0       -D      2D+1      -D        0   |    | 0  |
C     |                                                    |    | 0  |
C  N  |   0        0        0       -D      2D+1      -D   |    | 0  |
C     |                                                    |    | 0  |
C N+1 |   0        0        0        0        0        1   |    | F1 |
C      
C-----------------------------------------------------------------------
C  Advantage of scheme over simple attenuation factor: conserves total
C  number of photons - very useful when using scheme for heating rates.
C  Disadvantage: although reproduces e-folds very well for small flux
C  differences, starts to drift off when many orders of magnitude are
C  involved.
C-----------------------------------------------------------------------
      implicit none
      real*8 F0,F1,F(250)  !F(N+1)
      integer N
      integer I
      real*8 A,DX,D,DSQ,DDP1, B(101),R(101)
C
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
C
      DX = float(N)/A
      D = DX*DX
      DSQ = D*D
      DDP1 = D+D+1.d0
C
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
