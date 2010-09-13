!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: tropopause
!
! !DESCRIPTION: Subroutine TROPOPAUSE archives the ND55 tropopause diagnostic.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE TROPOPAUSE
!
! !USES:
!
      USE DAO_MOD,        ONLY : BXHEIGHT
      USE DAO_MOD,        ONLY : TROPP
      USE DIAG_MOD,       ONLY : AD55
      USE LOGICAL_MOD,    ONLY : LVARTROP
      USE PRESSURE_MOD,   ONLY : GET_PCENTER
      USE PRESSURE_MOD,   ONLY : GET_PEDGE
      USE TROPOPAUSE_MOD, ONLY : GET_TPAUSE_LEVEL

      IMPLICIT NONE

#     include "CMN_SIZE"  ! Size parameters
#     include "CMN_DIAG"  ! Diagnostic switches
!
! !REMARKS:
!  For GEOS-4, GEOS-5, 'MERRA', we use the tropopause pressure from the met 
!  field archive to determine if we are in the tropopause or not.  Therefore, 
!  the 3rd slot of AD55 should be archived with the tropopause pressure from 
!  the met fields.
!                                                                             .
!  For other met fields, we have to estimate the tropopause pressure from the
!  tropopause level.  Archive the pressure at the midpoint of the level in 
!  which the tropopause occurs.  NOTE: this may result in lower minimum 
!  tropopause pressure than reality. 
!
! !REVISION HISTORY:
!  30 Nov 1999 - H. Liu, R. Yantosca - Initial version
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
!  10 Sep 2010 - R. Yantosca - Added ProTeX headers
!  10 Sep 2010 - R. Yantosca - For GEOS-4, GEOS-5, MERRA met fields, take the
!                              the tropopause pressure directly from the
!                              met fields rather than computing it here.
!  10 Sep 2010 - R. Yantosca - Remove reference to LPAUSE, it's obsolete
!  10 Sep 2010 - R. Yantosca - Reorganize #if blocks for clarity
!EOP
!------------------------------------------------------------------------------
!BOC

#if   defined( GEOS_4 ) || defined( GEOS_5 ) || defined( MERRA )
!
! !LOCAL VARIABLES:
! 
      ! Scalars
      INTEGER :: I, J,    L,  L_TP
      REAL*8  :: H, FRAC, Pb, Pt

      !=================================================================
      ! %%%%% GEOS-4, GEOS-5, MERRA met fields %%%%%
      !
      ! We get tropopause pressure directly from the met field archive
      ! Compute tropopause height to be consistent w/ the pressure
      !=================================================================
      IF ( ND55 > 0 ) THEN

         ! Loop over surface grid boxes
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L_TP, H, Pb, Pt, FRAC )
         DO J = 1, JJPAR
         DO I = 1, IIPAR

            !---------------------------
            ! Compute quantities
            !---------------------------
    
            ! For this (I,J) column, get the level where the t'pause occurs
            L_TP = GET_TPAUSE_LEVEL( I, J )

            ! Get height (from surface to top edge) of all boxes that lie
            ! totally w/in the troposphere.  NOTE: Grid box (I,J,L_TP-1)
            ! is the highest purely tropospheric grid box in the column.
            H    = SUM( BXHEIGHT( I, J, 1:L_TP-1 ) )

            ! Get the pressures [hPa] at the bottom and top edges
            ! of the grid box in which the tropopause occurs
            Pb   = GET_PEDGE( I, J, L_TP   )  
            Pt   = GET_PEDGE( I, J, L_TP+1 )

            ! FRAC is the fraction of the grid box (I,J,L_TP) 
            ! that lies totally within the troposphere
            FRAC = ( Pb - TROPP(I,J) ) / ( Pb - Pt ) 

            ! Add to H the height [m] of the purely tropospheric 
            ! fraction of grid box (I,J,L_TP)
            H    = H + ( FRAC * BXHEIGHT(I,J,L_TP) )

            !---------------------------
            ! Archive into ND55 array
            !---------------------------
            AD55(I,J,1) = AD55(I,J,1) + L_TP        ! T'pause level
            AD55(I,J,2) = AD55(I,J,2) + H/1.0d3     ! T'pause height [km]
            AD55(I,J,3) = AD55(I,J,3) + TROPP(I,J)  ! T'pause pressure [hPa]

         ENDDO
         ENDDO
!$OMP END PARALLEL DO

      ENDIF

#else

!
! !LOCAL VARIABLES:
! 
      ! Scalars
      INTEGER :: I, J, L

      ! Arrays
      REAL*8  :: H(IIPAR,JJPAR,LLPAR)

      !=================================================================
      ! %%%%% ALL OTHER MET FIELDS %%%%%
      !
      ! We compute tropopause pressure from the tropopause level (which 
      ! is taken from the thermally-derived annual mean tropopause data 
      ! read from disk).
      !
      ! NOTE: Keep the existing algorithm for backwards compatibility.
      !=================================================================

      ! Find height of the midpoint of the first level
      ! H (in m) is the height of the midpoint of layer L (hyl, 03/28/99)
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
      !=================================================================
      IF ( ND55 > 0 ) THEN
         DO J = 1, JJPAR
         DO I = 1, IIPAR

            ! Get the tropopause level
            L           = GET_TPAUSE_LEVEL( I, J )

            ! If we are using the variable tropopause, then (I,J,L) is the
            ! highest purely tropospheric grid box.  The grid box in which
            ! the tropopause actually occurs is then (I,J,L+1).
            IF ( LVARTROP ) L = L + 1

            ! Archive level at which tropopause occurs
            AD55(I,J,1) = AD55(I,J,1) + L

            ! Archive tropopause height [km]
            AD55(I,J,2) = AD55(I,J,2) + H(I,J,L) / 1.0d3 ! m --> km

            ! We have to estimate the tropopause pressure from the 
            ! tropopause level.  Archive the pressure at the midpoint
            ! of the level in which the tropopause occurs.  NOTE: this may
            ! result in lower minimum tropopause pressure than reality.
            AD55(I,J,3) = AD55(I,J,3) + GET_PCENTER(I,J,L)

         ENDDO
         ENDDO
      ENDIF

#endif

      END SUBROUTINE TROPOPAUSE
!EOC
