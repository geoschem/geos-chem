C $Id: read_percentages.f,v 1.1 2003/06/30 20:26:04 bmy Exp $
C
C**********************************************************************
      SUBROUTINE read_percentages
C**********************************************************************
      IMPLICIT NONE
#     include "CMN_OH"
      INTEGER I,J,L,K
C****************************************************
C Created by Bryan Duncan.
C****************************************************
C This SR reads in information needed for SR AerosolOH.
C****************************************************
C List of Variables & Arrays
C
C PercentRed = array holding the percent reduction in OH
C     from an atmosphere with and without mineral dust 
C
C Randall Martin has not completed his runs for each month
C so PercentRed = 1 for those months.
C
C****************************************************
      print*,'Reading in info for OH reduction by mineral dust'
      print*,'in SR read_percentages.'
C
cbnd      OPEN(NGENERICFILE,FILE='PerRedOH',STATUS='OLD') 
C
cbnd      READ(NGENERICFILE,12) PercentRed
cbnd 12   FORMAT(12F5.2)
C
cbnd      CLOSE(NGENERICFILE)
       
      print*,'**********************'
      print*,'All mineral dust reductions set to 1.!!!!!' 
      print*,'**********************'
          PercentRed(:,:,:,:)=1.
C
C***************************************************************
C
      RETURN
      END
