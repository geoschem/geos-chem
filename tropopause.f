! $Id: tropopause.f,v 1.8 2006/11/07 19:02:08 bmy Exp $
      SUBROUTINE TROPOPAUSE
!
!******************************************************************************
!  Subroutine TROPOPAUSE defines the tropopause layer in terms of temperature 
!  lapse rates. (hyl, bmy, 11/30/99, 10/17/06)
!
!  NOTES:
!  (1 ) Make sure the DO-loops go in the order L-J-I, wherever possible.
!  (2 ) Now archive ND55 diagnostic here rather than in DIAG1.F.  Also,
!        use an allocatable array (AD55) to archive tropopause heights.
!  (3 ) HTPAUSE is now a local variable, since it is only used here.
!  (4 ) Make LTPAUSE a local variable, since LPAUSE is used to store
!        the annual mean tropopause. (bmy, 4/17/00)
!  (5 ) Replace PW(I,J) with P(I,J).  Also updated comments. (bmy, 10/3/01)
!  (6 ) Removed obsolete code from 9/01 and 10/01 (bmy, 10/24/01)
!  (7 ) Added polar tropopause for GEOS-3 in #if defined( GEOS_3 ) block 
!        (bmy, 5/20/02) 
!  (8 ) Replaced all instances of IM with IIPAR and JM with JJPAR, in order
!        to prevent namespace confusion for the new TPCORE (bmy, 6/25/02)
!  (9 ) Now use GET_PCENTER from "pressure_mod.f" to compute the pressure
!        at the midpoint of box (I,J,L).  Also deleted obsolete, commented-out
!        code. (dsa, bdf, bmy, 8/21/02)
!  (10) Now reference BXHEIGHT and T from "dao_mod.f".  Also reference routine
!        ERROR_STOP from "error_mod.f" (bmy, 10/15/02)
!  (11) Now uses routine GET_YMID from "grid_mod.f" to compute grid box 
!        latitude. (bmy, 2/3/03)
!  (12) Add proper polar tropopause level for GEOS-4 (bmy, 6/18/03)
!  (13) Remove support for GEOS-1 and GEOS-STRAT met fields (bmy, 8/4/06)
!  (14) Get tropopause level from TROPOPAUSE_MOD.F routines (phs, 10/17/06)
!******************************************************************************
!
      ! References to F90 modules.
      USE DAO_MOD,        ONLY : BXHEIGHT  !, T
      USE DIAG_MOD,       ONLY : AD55
      USE LOGICAL_MOD,    ONLY : LVARTROP
      USE PRESSURE_MOD,   ONLY : GET_PCENTER
      USE TROPOPAUSE_MOD, ONLY : GET_TPAUSE_LEVEL

      IMPLICIT NONE

#     include "CMN_SIZE"  ! Size parameters
#     include "CMN"       ! LPAUSE
#     include "CMN_DIAG"  ! Diagnostic switches
      
      ! Local variables
      INTEGER :: I, J, L 
      REAL*8  :: H(IIPAR,JJPAR,LLPAR)

      !=================================================================
      ! TROPOPAUSE begins here! 
      !
      ! H (in m) is the height of the midpoint of layer L (hyl, 03/28/99) 
      !=================================================================

      ! Find height of the midpoint of the first level
      DO J = 1, JJPAR
      DO I = 1, IIPAR
         H(I,J,1) = BXHEIGHT(I,J,1) / 2.d0
      ENDDO
      ENDDO

      ! Add to H 1/2 of the sum of the two adjacent boxheights
      DO L = 1, LLPAR-1
      DO J = 1, JJPAR
      DO I = 1, IIPAR
         H(I,J,L+1) = H(I,J,L) + 
     &               ( BXHEIGHT(I,J,L) + BXHEIGHT(I,J,L+1) ) / 2.d0
      ENDDO
      ENDDO
      ENDDO

      !=================================================================
      ! ND55: Tropopause level, height [ km ], and pressure [ mb ]
      !       Recall that PW(I,J) = PS(I,J) - PTOP
      !=================================================================
      IF ( ND55 > 0 ) THEN
         DO J = 1, JJPAR
         DO I = 1, IIPAR
            L           = GET_TPAUSE_LEVEL( I, J )
            IF ( LVARTROP ) L = L+1
            AD55(I,J,1) = AD55(I,J,1) + L
            AD55(I,J,2) = AD55(I,J,2) + H(I,J,L) / 1.0d3 ! m --> km
            AD55(I,J,3) = AD55(I,J,3) + GET_PCENTER(I,J,L)
         ENDDO
         ENDDO
      ENDIF

      ! Return to calling program
      END SUBROUTINE TROPOPAUSE
