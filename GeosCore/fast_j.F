!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: fast_j.f
!
! !DESCRIPTION: Subroutine FAST\_J loops over longitude and latitude, and
!  calls PHOTOJ to compute J-Values for each column at every chemistry 
!  time-step.  
!\\
!\\
!  References:
!
! \begin{enumerate}
! \item H. Liu, J.H. Crawford, R.B. Pierce, P. Norris, S.E. Platnick, G. Chen,
!       J.A. Logan, R.M. Yantosca, M.J. Evans, C. Kittaka, Y. Feng, and 
!       X. Tie, \emph{Radiative effect of clouds on tropospheric chemistry 
!       in a global three-dimensional chemical transport model}, 
!       \underline{J. Geophys. Res.}, \textbf{111}, D20303, 
!       doi:10.1029/2005JD006403, 2006. 
!       http://research.nianet.org/~hyl/publications/liu2006\_cloud1.abs.html
! \end{enumerate}
!
! !INTERFACE:
!
      SUBROUTINE FAST_J( SUNCOS, OD, ALBD )  
!
! !USES:
!
      USE DAO_MOD,      ONLY : T, CLDF
      USE ERROR_MOD,    ONLY : ERROR_STOP
      USE GRID_MOD,     ONLY : GET_YMID
      USE PRESSURE_MOD, ONLY : GET_PEDGE
      USE TIME_MOD,     ONLY : GET_MONTH, GET_DAY, GET_DAY_OF_YEAR
      USE TIME_MOD,     ONLY : GET_TAU,   GET_YEAR
      USE TOMS_MOD,     ONLY : READ_TOMS

#     include "define.h"

      USE CMN_SIZE_MOD, ONLY : IIPAR, JJPAR, LLPAR
      USE CMN_SIZE_MOD, ONLY : NDUST, MAXIJ, NAER, NRH
      USE CMN_FJ_MOD, ONLY : IPAR, JPAR, LPAR, JPMAX, JPPJ     ! IPAR, JPAR, LPAR, CMN_SIZE
      USE JV_CMN_MOD, ONLY : ODMDUST, PJ, NB, ODAER

      IMPLICIT NONE

!
! !INPUT PARAMETERS: 
!
      ! Cosine of solar zenith angle [unitless]
      REAL*8, INTENT(IN)    :: SUNCOS(MAXIJ)   

      ! Cloud optical depth [unitless]
      REAL*8, INTENT(IN)    :: OD(LLPAR,IIPAR,JJPAR) 

      ! UV albedo [unitless]
      REAL*8, INTENT(IN)    :: ALBD(IIPAR,JJPAR)    
!
! !REMARKS:
!
!  Parameter to choose cloud overlap algorithm:
!  ============================================================================
!  (1 ) OVERLAP (INTEGER) : 1 - Linear Approximation (used up to v7-04-12)
!                           2 - Approximate Random Overlap (default)
!                           3 - Maximum Random Overlap (computation intensive)
!
! !REVISION HISTORY:
!  01 Apr 1998 - P. Murti, R. Martin, R. Yantosca - Initial version
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
!  (15) Bug fix: PRES needs to be the true surface pressure for GEOS-4, but
!        PS-PTOP for all prior GEOS models.  (bmy, 2/6/04)
!  (16) Now account for cloud overlap (Maximum-Random Overlap and Random 
!        Overlap) in each column (hyl, phs, bmy, 9/18/07)
!  (17) Now initialize the PJ array here, instead of two layers below in
!        "set_prof.f".  Now no longer pass PRES to "photoj.f". (bmy, 11/29/07)
!  (18) Now switch to approx. random overlap option (hyl, phs, bmy, 10/7/08)
!  (19) Now can handle GEOS-5 reprocessed met data with OPTDEPTH being
!        in-cloud optical depths. (bmy, hyl, 10/24/08)
!  (10) Remove references to IN_CLOUD_OD (bmy, 10/15/09)
!  13 Aug 2010 - R. Yantosca - Added ProTeX headers
!  13 Aug 2010 - R. Yantosca - Treat MERRA in the same way as GEOS-5
!  08 Feb 2012 - R. Yantosca - Treat GEOS-5.7.x in the same way as MERRA
!  01 Mar 2012 - R. Yantosca - Now use GET_YMID(I,J,L) from grid_mod.F90
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER, SAVE      :: LASTMONTH = -1
      INTEGER            :: NLON, NLAT, DAY,  MONTH, DAY_OF_YR, L
      REAL*8             :: CSZA, PRES, SFCA, YLAT
      REAL*8             :: TEMP(LLPAR), OPTD(LLPAR)
      REAL*8             :: OPTDUST(LLPAR,NDUST)
      REAL*8             :: OPTAER(LLPAR,NAER*NRH)

      ! Local variables for cloud overlap (hyl, phs)
      INTEGER           :: NUMB, KK, I
      INTEGER            :: INDIC(LLPAR+1)
      INTEGER            :: INDGEN(LLPAR+1) = (/ (i,i=1,LLPAR+1) /)
      INTEGER            :: KBOT(LLPAR)
      INTEGER            :: KTOP(LLPAR)
      INTEGER            :: INDICATOR(LLPAR+2)
      REAL*8             :: FMAX(LLPAR)  ! maximum cloud fraction 
                                            !  in a block, size can be to 
                                            !  FIX(LLPAR)+1
      REAL*8             :: CLDF1D(LLPAR)
      REAL*8             :: ODNEW(LLPAR)

      ! NOTE: Switch from linear approximation (OVERLAP=1) to approximate
      ! random overlap (OVERLAP=2) because we have re-processed the GEOS-5
      ! met data such that OPTDEPTH, TAUCLI, and TAUCLW are now the in-cloud
      ! optical depths rather than the grid-box optical depths. 
      ! (hyl, phs, bmy, 10/7/08)
      INTEGER, PARAMETER :: OVERLAP = 2
      
      LOGICAL, SAVE      :: FIRST = .true.

      !=================================================================
      ! FAST_J begins here!
      !=================================================================

      ! Get day of year (0-365 or 0-366)
      DAY_OF_YR = GET_DAY_OF_YEAR()

      ! Get current month
      MONTH = GET_MONTH()

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
!$OMP+PRIVATE( NLON, NLAT,   YLAT,  CSZA,      OPTAER  )
!$OMP+PRIVATE( PRES, TEMP,   OPTD,  SFCA,      OPTDUST )
!$OMP+PRIVATE( FMAX, CLDF1D, KK,    NUMB,      L       )
!$OMP+PRIVATE( KBOT, KTOP,   ODNEW, INDICATOR, INDIC   )
!$OMP+SCHEDULE( DYNAMIC )

      ! Loop over latitudes
      DO NLAT = 1, JJPAR

         ! Loop over longitudes
         DO NLON = 1, IIPAR

            ! Grid box latitude [degrees]
            YLAT         = GET_YMID( NLON, NLAT, 1 )

            ! Cosine of solar zenith angle [unitless] at (NLON,NLAT)
            CSZA         = SUNCOS( (NLAT-1)*IIPAR + NLON ) 

            ! Define the PJ array here (bmy, 11/16/07)
            DO L = 1, NB
               PJ(L)     = GET_PEDGE( NLON, NLAT, L )
            ENDDO

            ! Top edge of PJ is top of atmosphere (bmy, 2/13/07)
            PJ(NB+1)     = 0d0

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

            !===========================================================
            ! CLOUD OVERLAP : LINEAR ASSUMPTION 
            ! Directly use OPTDEPTH = TAUCLD * CLDTOT
            ! 
            ! NOTE: Use this option if you want to compare to results
            !       from GEOS-Chem v7-04-12 and prior versions.
            !===========================================================
            IF ( OVERLAP == 1 ) then

#if   defined( GEOS_5 ) || defined( MERRA ) || defined( GEOS_57 )

               ! Column cloud fraction (not less than zero)
               CLDF1D = CLDF(1:LLPAR,NLON,NLAT)
               WHERE ( CLDF1D < 0d0 ) CLDF1D = 0d0
               
               ! NOTE: for the reprocessed GEOS-5 met fields (i.e. with
               ! optical depth & cloud fractions regridded with RegridTau)
               ! OPTD is the in-cloud optical depth.  At this point it has
               ! NOT been multiplied by cloud fraction yet.  Therefore,
               ! we can just apply the linear overlap formula as written 
               ! above (i.e. multiply by cloud fraction). (hyl, bmy, 10/24/08)
               OPTD = OPTD * CLDF1D
#endif

               ! Call FAST-J routines to compute J-values
               CALL PHOTOJ( NLON,  NLAT, YLAT,    DAY_OF_YR,  
     &                      MONTH, DAY,  CSZA,    TEMP,    
     &                      SFCA,  OPTD, OPTDUST, OPTAER )

            !===========================================================
            ! CLOUD OVERLAP : APPROXIMATE RANDOM OVERLAP
            ! Use OPTDEPTH = TAUCLD * CLDTOT**1.5
            !===========================================================
            ELSE IF ( OVERLAP == 2 ) THEN

               ! Column cloud fraction (not less than zero)
               CLDF1D = CLDF(1:LLPAR,NLON,NLAT)
               WHERE ( CLDF1D < 0d0 ) CLDF1D = 0d0
               
#if   defined( GEOS_5 ) || defined( MERRA ) || defined( GEOS_57 )

               ! NOTE: for the reprocessed GEOS-5 met fields (i.e. with
               ! optical depth & cloud fractions regridded with RegridTau)
               ! OPTD is the in-cloud optical depth.  At this point it has
               ! NOT been multiplied by cloud fraction yet.  Therefore,
               ! we can just apply the approximate random overlap formula
               ! as written above (i.e. multiply by cloud fraction^1.5).
               ! (hyl, bmy, 10/24/08)
               OPTD = OPTD * ( CLDF1D )**1.5d0
               
#else
               ! Otherwise, OPTD is the grid-box optical depth and has 
               ! already been multiplied by  the cloud fraction.  Therefore 
               ! we only need to multiply by the square root of the cloud 
               ! fraction here for the approximate random overlap option. 
               ! (hyl, bmy, 10/24/08)
               OPTD = OPTD * SQRT( CLDF1D )

#endif

               ! Call FAST-J routines to compute J-values
               CALL PHOTOJ( NLON,  NLAT, YLAT,    DAY_OF_YR,  
     &                      MONTH, DAY,  CSZA,    TEMP,  
     &                      SFCA,  OPTD, OPTDUST, OPTAER )

            !===========================================================
            ! CLOUD OVERLAP : MAXIMUM RANDOM OVERLAP
            !
            ! The Maximum-Random Overlap (MRAN) scheme assumes that 
            ! clouds in adjacent layers are maximally overlapped to 
            ! form a cloud block and that blocks of clouds separated by 
            ! clear layers are randomly overlapped.  A vertical profile 
            ! of fractional cloudiness is converted into a series of 
            ! column configurations with corresponding fractions 
            ! (see Liu et al., JGR 2006; hyl,3/3/04). 
            !
            ! For more details about cloud overlap assumptions and 
            ! their effect on photolysis frequencies and key oxidants 
            ! in the troposphere, refer to the following articles:
            ! 
            ! (1) Liu, H., et al., Radiative effect of clouds on 
            !      tropospheric chemistry in a global three-dimensional 
            !      chemical transport model, J. Geophys. Res., vol.111, 
            !      D20303, doi:10.1029/2005JD006403, 2006.
            ! (2) Tie, X., et al., Effect of clouds on photolysis and 
            !      oxidants in the troposphere, J. Geophys. Res., 
            !      108(D20), 4642, doi:10.1029/2003JD003659, 2003.
            ! (3) Feng, Y., et al., Effects of cloud overlap in 
            !      photochemical models, J. Geophys. Res., 109, 
            !      D04310, doi:10.1029/2003JD004040, 2004.
            ! (4) Stubenrauch, C.J., et al., Implementation of subgrid 
            !      cloud vertical structure inside a GCM and its effect 
            !      on the radiation budget, J. Clim., 10, 273-287, 1997.
            !-----------------------------------------------------------
            ! MMRAN needs IN-CLOUD optical depth (ODNEW) as input 
            ! Use cloud fraction, instead of OPTD, to form cloud blocks
            ! (hyl,06/19/04)
            !===========================================================
            ELSE IF ( OVERLAP == 3 ) THEN

               ! Initialize
               FMAX(:)   = 0d0  ! max cloud fraction in each cloud block
               ODNEW(:)  = 0d0  ! in-cloud optical depth
               CLDF1D    = CLDF(1:LLPAR,NLON,NLAT)
               INDICATOR = 0

               ! set small negative CLDF or OPTD to zero. 
               ! Set indicator vector.
               WHERE ( CLDF1D <= 0d0 ) 
                  CLDF1D               = 0d0
                  OPTD                 = 0D0
               ELSEWHERE
                  INDICATOR(2:LLPAR+1) = 1
               ENDWHERE

               ! Prevent negative opt depth
               WHERE ( OPTD < 0D0 ) OPTD   = 0D0

               !--------------------------------------------------------
               ! Generate cloud blocks & get their Bottom and Top levels
               !--------------------------------------------------------
               INDICATOR = CSHIFT(INDICATOR, 1) - INDICATOR
               INDIC     = INDICATOR(1:LLPAR+1)

               ! Number of cloud block
               NUMB      = COUNT( INDIC == 1 ) 
               
               ! Bottom layer of each block
               KBOT(1:NUMB) = PACK(INDGEN, (INDIC == 1 ) ) 

               ! Top layer of each block
               KTOP(1:NUMB) = PACK(INDGEN, (INDIC == -1) ) - 1 
             
               !--------------------------------------------------------
               ! For each cloud block, get Max Cloud Fractions, and 
               ! in-cloud optical depth vertical distribution.
               !--------------------------------------------------------
               DO KK = 1, NUMB

                  ! Max cloud fraction
                  FMAX(KK) = MAXVAL( CLDF1D(KBOT(KK):KTOP(KK)) )

#if   defined( GEOS_5 ) || defined( MERRA ) || defined( GEOS_57 )

                  ! NOTE: for the reprocessed GEOS-5 met fields (i.e. with
                  ! optical depth & cloud fractions regridded with RegridTau)
                  ! OPTD is the in-cloud optical depth.  At this point it has
                  ! NOT been multiplied by cloud fraction yet.  Therefore,
                  ! we can just set ODNEW = OPTD. (bmy, hyl, 10/24/08)

                  ! ODNEW is adjusted in-cloud OD vertical distrib.
                  ODNEW(KBOT(KK):KTOP(KK)) = OPTD(KBOT(KK):KTOP(KK))

#else

                  ! Otherwise, OPTD is the grid-box optical depth.  
                  ! Therefore, we must divide out by the cloud fraction
                  ! and thus set ODNEW = OPTD / FMAX. (bmy, hyl, 10/24/08)

                  ! ODNEW is adjusted in-cloud OD vertical distrib.
                  ODNEW(KBOT(KK):KTOP(KK)) = OPTD(KBOT(KK):KTOP(KK)) / 
     &                                       FMAX(KK)

#endif
               ENDDO
            
               !--------------------------------------------------------
               ! Apply Max RANdom if 1-6 clouds blocks, else use linear 
               !--------------------------------------------------------
               SELECT CASE( NUMB ) 
           
                  CASE( 0,7: )
                     CALL PHOTOJ( NLON,  NLAT, YLAT,    DAY_OF_YR, 
     &                            MONTH, DAY,  CSZA,    TEMP,  
     &                            SFCA,  OPTD, OPTDUST, OPTAER )

                  CASE( 1:6 ) 
                     CALL MMRAN_16( NUMB,  NLON,  NLAT,      YLAT,   
     &                              DAY,   MONTH, DAY_OF_YR, CSZA,    
     &                              TEMP,  SFCA,  OPTDUST,   OPTAER, 
     &                              LLPAR, FMAX,  ODNEW,     KBOT,   
     &                              KTOP )

               END SELECT
            ENDIF 
         ENDDO
      ENDDO
!$OMP END PARALLEL DO

      !-----------------------------------------------------------
      ! END OF SUBROUTINE FAST-J
      !-----------------------------------------------------------
      END SUBROUTINE FAST_J
!EOC
