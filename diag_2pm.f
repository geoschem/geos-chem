! $Id: diag_2pm.f,v 1.1 2003/06/30 20:26:03 bmy Exp $
      SUBROUTINE DIAG_2PM
!
!*****************************************************************************
!  Subroutine DIAG_2PM (bmy, 3/26/99, 2/27/02) constructs the diagnostic 
!  flag arrays : 
!      LTJV  : J-values           (ND22)
!      LTOH  : OH concentrations  (ND43)
!      LTNO  : NO concentrations  (ND43)
!      LTNO2 : NO2 concentrations (ND43)
!      LTHO2 : HO2 concentrations (ND43)
!      LTOTH : used for tracers   (ND45) 
!
!  These arrays are either 1 (if it is within a certain time interval)   
!  or 0 (if it is not within a certain time interval).  The limits of
!  the time intervals for CTOTH and CTJV are now defined in input.geos
!  The arrays CTOTH, CTOH, CTNO, CTJV count the number of times the 
!  diagnostics are accumulated for each grid box (i.e LTOTH is 1)
!
!  Arguments as input:
!  ===========================================================================
!  (1) NMIN : Number of minutes elapsed since the run began 
!
!  NOTES:
!  (1 ) Now use F90 syntax (bmy, 3/26/99)
!  (2 ) Now reference LTNO2, CTNO2, LTHO2, CTHO2 arrays from "diag_mod.f".  
!        Updated comments, cosmetic changes.  (rvm, bmy, 2/27/02)
!  (3 ) Now removed NMIN from the arg list.  Now use functions GET_LOCALTIME,
!        ITS_TIME_FOR_CHEM, ITS_TIME_FOR_DYN from "time_mod.f" (bmy, 2/11/03)
!*****************************************************************************
!     
      ! References to F90 modules
      USE DIAG_MOD, ONLY : LTJV,  CTJV,  LTNO,  CTNO,
     &                     LTOH,  CTOH,  LTOTH, CTOTH, LTNO2,
     &                     CTNO2, LTHO2, CTHO2, LTNO3, CTNO3
      USE TIME_MOD, ONLY : GET_LOCALTIME, 
     &                     ITS_TIME_FOR_DYN, ITS_TIME_FOR_CHEM

      IMPLICIT NONE

#     include "CMN_SIZE"
#     include "CMN_DIAG"
#     include "CMN_SETUP"

      ! Local variables
      INTEGER             :: I, J, NMIN
      REAL*8              :: XLOCTM(IIPAR,JJPAR) !, TMP(IIPAR) 

      !=================================================================
      ! DIAG_2PM begins here!
      !=================================================================

      ! Get local time
      DO J = 1, JJPAR
      DO I = 1, IIPAR
         XLOCTM(I,J) = GET_LOCALTIME( I ) 
      ENDDO
      ENDDO
         
      !=================================================================
      ! Diagnostics that are done every dynamic timestep
      !=================================================================
      IF ( ITS_TIME_FOR_DYN() ) THEN

         !-----------------------------
         ! ND45 -- mixing ratios
         !-----------------------------
         IF ( ND45 > 0 ) THEN

            ! LTOTH denotes where LT is between HR1_OTH and HR2_OTH
            ! CTOTH counts the times when LT was between HR1_OTH and HR2_OTH
            WHERE( XLOCTM >= HR1_OTH .and. XLOCTM <= HR2_OTH ) 
               LTOTH = 1
               CTOTH = CTOTH + 1
            ELSEWHERE
               LTOTH = 0
            ENDWHERE
         ENDIF
      ENDIF

      !=================================================================
      ! Diagnostics that are done every chemistry timestep
      !=================================================================
      IF ( ITS_TIME_FOR_CHEM() ) THEN 
        
         !-----------------------------
         ! ND22 -- J-Value diagnostic
         !-----------------------------
         IF ( ND22 > 0 ) THEN

            ! LTJV denotes where LT is between HR1_JV and HR2_JV
            ! CTJV counts the times when LT was between HR1_JV and HR2_JV
            WHERE( XLOCTM >= HR1_JV .and. XLOCTM <= HR2_JV )
               LTJV = 1
               CTJV = CTJV + 1               
            ELSEWHERE
               LTJV = 0 
            ENDWHERE
         ENDIF
            
         !-----------------------------
         ! ND43 -- OH, NO, NO2, HO2
         !-----------------------------
         IF ( ND43 > 0 ) THEN

            ! LTNO denotes where LT is between HR1_NO and HR2_NO
            ! CTNO counts the times when LT was between HR1_NO and HR2_NO
            ! Now set LTNO2, CTNO2 based on the NO times (rvm, bmy, 2/27/02)
            WHERE( XLOCTM >= HR1_NO .and. XLOCTM <= HR2_NO ) 
               LTNO  = 1
               CTNO  = CTNO + 1
               LTNO2 = 1
               CTNO2 = CTNO2 + 1
            ELSEWHERE
               LTNO  = 0
               LTNO2 = 0
            ENDWHERE

            ! LTNO denotes where LT is between HR1_OH and HR2_OH
            ! CTNO counts the times when LT was between HR1_OH and HR2_OH
            ! Now set LTHO2, CTHO2 based on the OH times (rvm, bmy, 2/27/02)
            WHERE( XLOCTM >= HR1_OH .and. XLOCTM <= HR2_OH ) 
               LTOH  = 1
               CTOH  = CTOH + 1
               LTHO2 = 1
               CTHO2 = CTHO2 + 1
               LTNO3 = 1
               CTNO3 = CTNO3 + 1
            ELSEWHERE
               LTOH  = 0
               LTHO2 = 0
               LTNO3 = 0
            ENDWHERE
         ENDIF
      ENDIF

      ! Return to calling program
      END SUBROUTINE DIAG_2PM
