! $Id: fast_j.f,v 1.2 2003/07/21 15:09:26 bmy Exp $
      SUBROUTINE FAST_J( SUNCOS, OD, ALBD )  
!
!******************************************************************************
!  Subroutine FAST_J loops over longitude and latitude, and calls PHOTOJ 
!  to compute J-Values for each column at every chemistry time-step.  
!  (ppm, 4/98; bmy, rvm, 9/99, 7/17/03)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) SUNCOS (REAL*8) : Cosine of solar zenith angle [unitless]
!  (2 ) OD     (REAL*8) : Cloud optical depth          [unitless]
!  (3 ) ALBD   (REAL*8) : UV albedo                    [unitless]
!  
!  NOTES:
!  (1 ) Call this routine EACH chemistry time-step, before solver.
!  (2 ) This routine must know IMAX, JMAX, LMAX. 
!  (3 ) Now use new !$OMP compiler directives for parallelization (bmy, 5/2/00)
!  (4 ) Now reference "cmn_fj.h" and "jv_cmn.h" for the aerosol
!        optical depths (bmy, 10/2/00)
!  (5 ) Add OPTDUST as a local variable -- make OPTDUST private for
!        the parallel DO-loop, since it stores 1 column of aerosol optical
!        depth for each dust type (bmy, rvm, 10/2/00)
!  (6 ) For now, LPAR in "cmn_fj.h" = LGLOB in "CMN_SIZE".  Therefore we 
!        assume that we are always doing global runs. (bmy, 10/2/00)
!  (7 ) Removed obsolete code from 10/2/00 (bmy, 12/21/00)
!  (8 ) Replace {IJL}GLOB w/ IIPAR,JJPAR,LLPAR everywhere.  Also YLMID(NLAT)
!        needs to be referenced by YLMID(NLAT+J0). (bmy, 9/26/01)
!  (9 ) Remove obsolete code from 9/01.  Updated comments. (bmy, 10/24/01)
!  (10) Add OPTAER as a local variable, make it private for the parallel
!        DO loop, since it stores 1 column of aerosol optical depths for each
!        aerosol type.  Pass OPTAER to PHOTOJ via the argument list.  Declare
!        OPTAER as PRIVATE for the parallel DO-loop. (rvm, bmy, 2/27/02)
!  (11) Now reference GET_PEDGE from "pressure_mod.f", which returns the
!        correct "floating" pressure. (dsa, bdf, bmy, 8/20/02)
!  (12) Now reference T from "dao_mod.f" (bmy, 9/23/02)
!  (13) Now uses routine GET_YMID from "grid_mod.f" to compute grid box 
!        latitude.  Now make IDAY, MONTH local variables.  Now use function 
!        GET_DAY_OF_YEAR from "time_mod.f".  Bug fix: now IDAY (as passed to
!        photoj.f) is day of year rather than cumulative days since Jan 1, 
!        1985. (bmy, 2/11/03)
!  (14) Now reference routine GET_YEAR from "time_mod.f".  Added LASTMONTH
!        as a SAVEd variable.  Now call READ_TOMSO3 from "toms_mod.f" at the
!        beginning of a new month (or the first timestep) to read TOMS O3
!        columns which will be used by "set_prof.f".  Now also reference
!        routine GET_DAY from "time_mod.f".  Rename IDAY to DAY_OF_YR. Pass 
!        day of month to PHOTOJ.  Updated comments, cosmetic changes.
!        (bmy, 7/17/03)
!******************************************************************************
!
      ! References to F90 modules
      USE DAO_MOD,      ONLY : T
      USE GRID_MOD,     ONLY : GET_YMID
      USE PRESSURE_MOD, ONLY : GET_PEDGE
      USE TIME_MOD,     ONLY : GET_MONTH, GET_DAY, GET_DAY_OF_YEAR, 
     &                         GET_TAU,   GET_YEAR
      USE TOMS_MOD,     ONLY : READ_TOMS
      
      IMPLICIT NONE

#     include "cmn_fj.h" ! IPAR, JPAR, LPAR, CMN_SIZE
#     include "CMN"      ! P, T, YLMID            
#     include "jv_cmn.h" ! ODMDUST 

      ! Arguments
      REAL*8, INTENT(IN) :: SUNCOS(MAXIJ)
      REAL*8, INTENT(IN) :: OD(LLPAR,IIPAR,JJPAR) 
      REAL*8, INTENT(IN) :: ALBD(IIPAR,JJPAR)     
      
      ! Local variables
      INTEGER, SAVE      :: LASTMONTH = -1
      INTEGER            :: NLON, NLAT, DAY,  MONTH, DAY_OF_YR
      REAL*8             :: CSZA, PRES, SFCA, YLAT
      REAL*8             :: TEMP(LLPAR), OPTD(LLPAR)
      REAL*8             :: OPTDUST(LLPAR,NDUST)
      REAL*8             :: OPTAER(LLPAR,NAER*NRH)

      !=================================================================
      ! FAST_J begins here!
      !=================================================================

      ! Get day of year (0-365 or 0-366)
      DAY_OF_YR = GET_DAY_OF_YEAR()

      ! Get current month
      MONTH     = GET_MONTH()

      ! Get day of month
      DAY       = GET_DAY()

      ! Read TOMS O3 columns if it's a new month
      IF ( MONTH /= LASTMONTH ) THEN
         CALL READ_TOMS( MONTH, GET_YEAR() )
         LASTMONTH = MONTH
      ENDIF

      !=================================================================
      ! For each (NLON,NLAT) location, call subroutine PHOTOJ (in a 
      ! parallel loop to compute J-values for the entire column.  
      ! J-values will be stored in the common-block variable ZPJ, and 
      ! will be later accessed via function FJFUNC. 
      !=================================================================
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( NLON, NLAT, YLAT, CSZA, OPTAER  )
!$OMP+PRIVATE( PRES, TEMP, OPTD, SFCA, OPTDUST )
!$OMP+SCHEDULE( DYNAMIC )
      DO NLAT = 1, JJPAR

         ! Grid box latitude [degrees]
         YLAT = GET_YMID( NLAT )

         DO NLON = 1, IIPAR

            ! Cosine of solar zenith angle [unitless] at (NLON,NLAT)
            CSZA         = SUNCOS( (NLAT-1)*IIPAR + NLON ) 

            ! Surface pressure - PTOP [hPa] at (NLON,NLAT)
            PRES         = GET_PEDGE(NLON,NLAT,1) - PTOP

            ! Temperature profile [K] at (NLON,NLAT)
            TEMP         = T(NLON,NLAT,1:LLPAR)

            ! Surface albedo [unitless] at (NLON,NLAT)
            SFCA         = ALBD(NLON,NLAT)

            ! Aerosol OD profile [unitless] at (NLON,NLAT)
            OPTAER(:,:)  = ODAER(NLON,NLAT,:,:)

            ! Mineral dust OD profile [unitless] at (NLON,NLAT)
            OPTDUST(:,:) = ODMDUST(NLON,NLAT,:,:)

            ! Cloud OD profile [unitless] at (NLON,NLAT)
            OPTD         = OD(1:LLPAR,NLON,NLAT)

            !-----------------------------------------------------------
            !### If you want to exclude aerosol OD, mineral dust OD,
            !### or cloud OD, then uncomment the following lines:
            !OPTAER  = 0d0
            !OPTDUST = 0d0
            !OPTD    = 0d0
            !-----------------------------------------------------------

            ! Call FAST-J routines to compute J-values
!-----------------------------------------------------------------------
! Prior to 7/16/03:
! Now pass DAY_OF_YR and DAY to PHOTOJ.  Also swap order of
! CZSA and SFCA arguments in the call. (bmy, 7/16/03)
!            CALL PHOTOJ( NLON, NLAT, YLAT, IDAY, MONTH,   CSZA, 
!     &                   PRES, TEMP, OPTD, SFCA, OPTDUST, OPTAER )
!-----------------------------------------------------------------------
            CALL PHOTOJ( NLON, NLAT, YLAT, DAY_OF_YR,  MONTH,   DAY,
     &                   CSZA, PRES, TEMP, SFCA, OPTD, OPTDUST, OPTAER )

         ENDDO
      ENDDO
!$OMP END PARALLEL DO
 
      ! Return to calling program
      END SUBROUTINE FAST_J

 
