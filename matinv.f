C $Id: matinv.f,v 1.1 2003/06/30 20:26:06 bmy Exp $
      SUBROUTINE MATINV(B)

      IMPLICIT NONE

C**********************************************************************
C                                                                     *
C  HARVARD TROPOSPHERIC CHEMISTRY MODULE FOR 3-D APPLICATIONS         *
C  by Larry Horowitz, Jinyou Liang, Gerry Gardner, Prof. Daniel Jacob *
C  of HARVARD UNIVERSITY    (Release V2.0)                            *
C**********************************************************************
C*                  
C*    Matrix inversion routine for 3*3 matrix; use LU decomposition
C*                        
C*            
      REAL*8 B(3,3)
C*  
C*    Set up L and U
      B(2,1) = B(2,1)/B(1,1)
      B(2,2) = B(2,2) - B(2,1)*B(1,2)
      B(2,3) = B(2,3) - B(2,1)*B(1,3)
      B(3,1) = B(3,1)/B(1,1)
      B(3,2) = (B(3,2) - B(3,1)*B(1,2))/B(2,2)
      B(3,3) = B(3,3) - B(3,1)*B(1,3) - B(3,2)*B(2,3)
C*   
C*    Invert L
      B(3,2) = -B(3,2)
      B(3,1) = -B(3,1) - B(3,2)*B(2,1)
      B(2,1) = -B(2,1)
C*   
C*    Invert U
      B(3,3) = 1.0D0/B(3,3)
      B(2,3) = -B(2,3)*B(3,3)/B(2,2)
      B(2,2) = 1.0D0/B(2,2)
      B(1,3) = -(B(1,2)*B(2,3)+B(1,3)*B(3,3))/B(1,1)
      B(1,2) = -B(1,2)*B(2,2)/B(1,1)
      B(1,1) = 1.0D0/B(1,1)
C*   
C*    Multiply (U-Inverse)*(L-Inverse)
      B(1,1) = B(1,1)+B(1,2)*B(2,1)+B(1,3)*B(3,1)
      B(1,2) = B(1,2)+B(1,3)*B(3,2)
      B(2,1) = B(2,2)*B(2,1)+B(2,3)*B(3,1)
      B(2,2) = B(2,2)+B(2,3)*B(3,2)
      B(3,1) = B(3,3)*B(3,1)
      B(3,2) = B(3,3)*B(3,2)
C*   
C*    Return to SCATTR
      RETURN
      END
