C $Id: AerosolOH.f,v 1.1 2003/06/30 20:26:03 bmy Exp $
C
C**********************************************************************
      SUBROUTINE AerosolOH(BOH,I,J,L,K)
C**********************************************************************
      IMPLICIT NONE
#     include "CMN_OH"
      INTEGER I,J,L,K
C****************************************************
C Created by Bryan Duncan.
C****************************************************
C This SR calculates the percent reduction of OH from the
C presence of absorbing/reflecting mineral dust.  The percentages
C are from the work of Randall Martin.
C
C Randall Martin has not completed his runs for each month
C so PercentRed = 1 for those months.
C
C****************************************************
C List of Variables & Arrays
C
C BOH = array holding tropospheric parameterized OH.
C
      REAL*8 BOH(MLONBX,MLATBX,MVRTBX)
C
C PercentRed = array holding the percent reduction in OH
C     from an atmosphere with and without mineral dust 
C
C****************************************************
C
      BOH(I,J,L) = BOH(I,J,L) * PercentRed(I,J,L,K) 
C
C***************************************************************
C
      RETURN
      END
