! $Id: diag41.f,v 1.4 2003/11/06 21:07:18 bmy Exp $
      SUBROUTINE DIAG41 
!
!******************************************************************************
!  Subroutine DIAG41 produces monthly mean boundary layer height in meters 
!  between 1200-1600 local time for the U.S. geographical domain. 
!  (amf, swu, bmy, 11/18/99, 11/6/03)
!
!  Input via "CMN" header file:
!  ===========================================================================
!  (1 ) XTRA2 : Height of PBL in boxes
!
!  NOTES:
!  (1 ) DIAG41 is written in Fixed-Format F90. 
!  (2 ) XTRA2 must be computed by turning TURBDAY on first.  Also,
!        XTRA2 is a global-size array, so use window offsets IREF, JREF
!        to index it correctly. (bmy, 11/18/99)
!  (3 ) Do a little rewriting so that the DO-loops get executed
!        in the correct order (J first, then I). (bmy, 11/18/99)
!  (4 ) AD41 is now declared allocatable in "diag_mod.f". (bmy, 12/6/99)
!  (5 ) AFTTOT is now declared allocatable in "diag_mod.f". (bmy, 3/17/00)
!  (6 ) Remove NYMD from the argument list -- it wasn't used (bmy, 6/22/00) 
!  (7 ) XTRA2(IREF,JREF,5) is now XTRA2(I,J).  Also updated comments. 
!        (bmy, 9/25/01)
!  (8 ) Removed obsolete code from 9/01 (bmy, 10/23/01)
!  (9 ) Now reference BXHEIGHT from "dao_mod.f".  Also removed obsolete
!        code. (bmy, 9/18/02)
!  (10) Now use function GET_LOCALTIME from "dao_mod.f" (bmy, 2/11/03)
!  (11) Bug fix in DO-loop for calculating local time (bmy, 7/9/03)
!  (12) For GEOS-4, PBL depth is already in meters, so we only have to
!        multiply that by the GOOD array.  Also now references PBL array
!        from "dao_mod.f".  Bug fix: now use barometric law to compute PBL 
!        height in meters for GEOS-1, GEOS-STRAT, GEOS-3.  This eliminates an 
!        overprediction of the PBL height. (swu, bmy, 11/6/03)
!******************************************************************************
!
      ! References to F90 modules
      USE DAO_MOD,      ONLY : PBL
      USE DIAG_MOD,     ONLY : AD41, AFTTOT
      USE PRESSURE_MOD, ONLY : GET_PEDGE
      USE TIME_MOD,     ONLY : GET_LOCALTIME

      IMPLICIT NONE

#     include "CMN_SIZE"   ! IIPAR, DISIZE
#     include "CMN_GCTM"   ! SCALE_HEIGHT
#     include "CMN_DIAG"   ! ND41, AFTTOT
#     include "CMN"        ! XTRA2

      ! Local variables
      INTEGER             :: GOOD(IIPAR)
      !-----------------------------------------------------------------
      ! Prior to 11/6/03:
      !INTEGER             :: I, IREF, J, JREF, L, M, PBLINT
      !-----------------------------------------------------------------
      INTEGER             :: I, J, L, M

      REAL*8              :: XLOCTM(IIPAR)
      !-----------------------------------------------------------------
      ! Prior to 11/6/03:
      !REAL*8              :: TEMPBL, PBLDEC
      !-----------------------------------------------------------------
      REAL*8              :: TEMPBL, P_SURF, P_BLTOP

      !=================================================================
      ! DIAG41 begins here!
      !=================================================================

      ! XLOCTM: Local time at each longitude 
      DO I = 1, IIPAR
         XLOCTM(I) = GET_LOCALTIME( I ) 
      ENDDO

      ! GOOD=1 denotes longitudes where it is between 1200-1600 LT
      GOOD = 0
      WHERE ( XLOCTM >= 12d0 .and. XLOCTM <= 16d0 ) GOOD = 1

      ! Loop over surface grid boxes
      DO J = 1, JJPAR
      DO I = 1, IIPAR

#if   defined( GEOS_4 ) 
      
         !==============================================================
         ! GEOS-4: PBL is already in [m]
         !==============================================================
         TEMPBL = PBL(I,J) * GOOD(I)

#else

!----------------------------------------------------------------------------
! Prior to 11/6/03:
! Old code caused an overestimate of PBL height (swu, bmy, 11/6/03)
!         TEMPBL = 0d0
!
!         ! Integer and fractional parts of XTRA2
!         PBLINT = FLOOR( XTRA2(I,J) )
!         PBLDEC = XTRA2(I,J) - PBLINT
!
!         ! "Full" grid boxes
!         IF ( PBLINT > 0 ) THEN
!            DO M = 1, PBLINT 
!               TEMPBL = TEMPBL + ( BXHEIGHT(I,J,M) * GOOD(I) )
!            ENDDO
!         ENDIF
!
!         ! "Fraction" of the highest grid box
!         TEMPBL = TEMPBL + 
!     &            ( PBLDEC * BXHEIGHT(I,J,PBLINT+1) * GOOD(I) )
!----------------------------------------------------------------------------

         !==============================================================
         ! GEOS-1, GEOS-STRAT, GEOS-3; PBL is in [hPa]
         ! so we must use the barometric law to convert to [m]
         !==============================================================

         ! Surface pressure [hPa]
         P_SURF  = GET_PEDGE(I,J,1)   

         ! Pressure at BL top [hPa]
         P_BLTOP = P_SURF - PBL(I,J)    

         ! Convert PBL top pressure to [m]
         TEMPBL  = LOG( P_SURF / P_BLTOP ) * SCALE_HEIGHT * GOOD(I)

#endif

         !==============================================================
         ! Store to diagnostic arrays 
         !==============================================================
         AD41(I,J,1) = AD41(I,J,1) + TEMPBL
         AD41(I,J,2) = AD41(I,J,2) + ( XTRA2(I,J) * GOOD(I) )

         ! Increment counter of afternoon boxes
         AFTTOT(I,J) = AFTTOT(I,J) + GOOD(I)

      ENDDO    
      ENDDO    

      ! Return to calling program
      END SUBROUTINE DIAG41
