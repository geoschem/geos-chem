! $Id: fast_j.f,v 1.1 2003/06/30 20:26:03 bmy Exp $
      SUBROUTINE FAST_J( SUNCOS, OD, ALBD )  
!
!******************************************************************************
!  Subroutine FAST_J loops over longitude and latitude, and calls PHOTOJ 
!  to compute J-Values for each column at every chemistry time-step.  
!  (ppm, 4/98; bmy, rvm, 9/99, 2/11/03)
!
!  Inputs to PHOTOJ:
!  (take from argument list or common blocks, or a combination of both)
!  ===========================================================================
!  Variable  Type    Dimension  Units   Description
!  --------  ----    ---------  -----   -----------
!  nlon      int        -         -     Longitude index
!  nlat      int        -         -     Lattitude index
!  month     int        -         -     Month of year (1-12)
!  csza      dble       -         -     Cosine of solar zenith angle 
!                                        at nlon, nlat
!  pres      dble       -        [mb]   Column pressure at nlon, nlat
!  temp      dble    [LMAX]      [K]    Layer temperatures at nlon, nlat 
!  optd      dble    [LMAX]       -     Layer optical depths at nlon, nlat
!  sfca      dble       -         -     Surface albedo at nlon, nlat  
!  optdust   dble    [LMAX,NDUST] -     Dust optical depths 
!                                        (for NDUST dust types)
!  optaer    dble [LMAX,NAER*NRH]  -     Aerosol optical depths
!                                        (for NAER aerosol types)
!  
!
!  NOTES:
!  ======
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
!******************************************************************************
!
      ! References to F90 modules
       USE DAO_MOD,      ONLY : T
      USE GRID_MOD,     ONLY : GET_YMID
      USE PRESSURE_MOD, ONLY : GET_PEDGE
      USE TIME_MOD,     ONLY : GET_MONTH, GET_DAY_OF_YEAR, GET_TAU

      IMPLICIT NONE

#     include "cmn_fj.h" ! IPAR, JPAR, LPAR, CMN_SIZE
#     include "CMN"      ! P, T, YLMID            1
#     include "jv_cmn.h" ! ODMDUST 

      ! Arguments
      REAL*8, INTENT(IN) :: SUNCOS(MAXIJ)
      REAL*8, INTENT(IN) :: OD(LLPAR,IIPAR,JJPAR) 
      REAL*8, INTENT(IN) :: ALBD(IIPAR,JJPAR)     
      
      ! Local variables
      INTEGER            :: NLON, NLAT, IDAY, MONTH
      REAL*8             :: CSZA, PRES, SFCA, YLAT
      REAL*8             :: TEMP(LLPAR), OPTD(LLPAR)
      REAL*8             :: OPTDUST(LLPAR,NDUST)
      REAL*8             :: OPTAER(LLPAR,NAER*NRH)

      !=================================================================
      ! FAST_J begins here!
      !=================================================================

      ! Get day of year (0-365 or 0-366)
      IDAY  = GET_DAY_OF_YEAR()

      ! Get current month
      MONTH = GET_MONTH()

      !=================================================================
      ! For each (I,J) location, call subroutine PHOTOJ (in a parallel 
      ! loop to compute J-values for the entire column.  J-values will 
      ! be stored in the common-block variable ZPJ, and will be later 
      ! accessed via function FJFUNC. 
      !
      ! The parallel loop only takes effect if you compile with the 
      ! f90 "-mp" switch.  Otherwise the compiler will interpret the 
      ! parallel-processing directives as comments. (bmy, 5/2/00)
      !
      ! Add OPTDUST, OPTAER to PHOTOJ as arguments -- declare them
      ! PRIVATE for the parallel DO loop (rvm, bmy, 10/2/00, 2/27/02)
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
            CSZA = SUNCOS( (NLAT-1)*IIPAR + NLON ) 
            PRES = GET_PEDGE(NLON,NLAT,1) - PTOP
            TEMP = T(NLON,NLAT,1:LLPAR)
            !-----------------------------------------------------------
            ! To include aerosols, uncomment the following line:
            OPTAER(:,:) = ODAER(NLON,NLAT,:,:)
            !-----------------------------------------------------------
            ! To exclude aerosols, uncomment the following line:
            !OPTAER = 0d0
            !-----------------------------------------------------------
            ! To include mineral dust, uncomment the following line:
            OPTDUST(:,:) = ODMDUST(NLON,NLAT,:,:)
            !-----------------------------------------------------------
            ! To exclude mineral dust, uncomment the following line:
            !OPTDUST = 0d0
            !-----------------------------------------------------------
            ! To include clouds, uncomment the following line:              
            OPTD = OD(1:LLPAR,NLON,NLAT)
            !-----------------------------------------------------------
            ! For clear sky, uncomment the following line:                  
            !OPTD = 0d0
            !-----------------------------------------------------------
            SFCA = ALBD(NLON,NLAT)

            ! Add OPTDUST to the argument list (bmy, rvm, 10/2/00)
            ! Add OPTAER to the argument list (rvm, bmy, 2/27/02)
            CALL PHOTOJ( NLON, NLAT, YLAT, IDAY, MONTH,   CSZA, 
     &                   PRES, TEMP, OPTD, SFCA, OPTDUST, OPTAER )

         ENDDO
      ENDDO
!$OMP END PARALLEL DO
 
      ! Return to calling program
      END SUBROUTINE FAST_J

 
