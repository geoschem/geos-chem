!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: fast_jx_mod
!
! !DESCRIPTION: Module FAST\_JX\_MOD contains routines and variables for 
!  calculating photolysis rates using the Fast-JX scheme (Prather et al).
!  Current implementation is version 6.2.
!\\
!\\
! !INTERFACE: 
!
      MODULE FAST_JX_MOD
!
! !USES:
!
      USE CMN_FJX_MOD

      IMPLICIT NONE
#     include "define.h"

      PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
      PUBLIC  :: EXITC
      PUBLIC  :: SOLAR_JX
      PUBLIC  :: ACLIM_FJX
      PUBLIC  :: JP_ATMO
      PUBLIC  :: PHOTO_JX
      PUBLIC  :: INIT_FJX
      PUBLIC  :: FAST_JX

!      PUBLIC  :: FAST_J
!      PUBLIC  :: FJFUNC
!      PUBLIC  :: INPHOT
!      PUBLIC  :: JRATET
!      PUBLIC  :: OPMIE
!      PUBLIC  :: PHOTOJ
!
! !PRIVATE MEMBER FUNCTIONS:
!
      PRIVATE :: OPMIE
      PRIVATE :: MIESCT 
      PRIVATE :: LEGND0
      PRIVATE :: BLKSLV 
      PRIVATE :: GEN_ID 
      PRIVATE :: OPTICL
      PRIVATE :: OPTICA
      PRIVATE :: OPTICM 
      PRIVATE :: JRATET
      PRIVATE :: X_INTERP 
      PRIVATE :: JP_ATM
      PRIVATE :: SPHERE2 
      PRIVATE :: EXTRAL 
      PRIVATE :: RD_XXX
      PRIVATE :: RD_CLD
      PRIVATE :: RD_MIE
      PRIVATE :: RD_UM
      PRIVATE :: RD_PROF
      PRIVATE :: RD_JS_JX
!
! !REVISION HISTORY:
!  27 Mar 2013 - S. D. Eastham - Initial version (based on original GEOS-
!  Chem files and instructions from J. Mao)
!  28 Mar 2013 - S. D. Eastham - Upgraded to Fast-JX v7.0
!  
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !PRIVATE TYPES:
!
      CONTAINS      
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: fast_jx
!
! !DESCRIPTION: Subroutine FAST\_JX loops over longitude and latitude, and
!  calls PHOTO\_JX to compute J-Values for each column at every chemistry 
!  time-step.  
!\\
!\\
! !INTERFACE:
!
#if defined( DEVEL ) || defined( EXTERNAL_GRID ) || defined( EXTERNAL_FORCING )
      SUBROUTINE FAST_JX( SUNCOS, OD, ALBD, PEDGE, am_I_Root )  
#else
      SUBROUTINE FAST_JX( SUNCOS, OD, ALBD, am_I_Root )  
#endif
!
! !USES:
!
#     include "define.h"

      USE DAO_MOD,      ONLY : T, CLDF
      USE ERROR_MOD,    ONLY : ERROR_STOP, ALLOC_ERR
      USE GRID_MOD,     ONLY : GET_YMID
      USE PRESSURE_MOD, ONLY : GET_PEDGE
      USE TIME_MOD,     ONLY : GET_MONTH, GET_DAY, GET_DAY_OF_YEAR
      USE TIME_MOD,     ONLY : GET_TAU,   GET_YEAR
      USE TOMS_MOD,     ONLY : GET_OVERHEAD_O3

      USE CMN_SIZE_MOD, ONLY : IIPAR,   JJPAR, LLPAR
      USE CMN_SIZE_MOD, ONLY : NDUST,   MAXIJ, NAER,  NRH

      IMPLICIT NONE
!
! !INPUT PARAMETERS: 
!
      ! Is this the root CPU?
      LOGICAL, INTENT(IN) :: am_I_Root

      ! Cosine of solar zenith angle [unitless]
      REAL*8,  INTENT(IN) :: SUNCOS(MAXIJ)   

      ! Cloud optical depth [unitless]
      REAL*8,  INTENT(IN) :: OD(LLPAR,IIPAR,JJPAR) 

      ! UV albedo [unitless]
      REAL*8,  INTENT(IN) :: ALBD(IIPAR,JJPAR)    

#if defined( DEVEL ) || defined( EXTERNAL_GRID ) || defined( EXTERNAL_FORCING )
      ! Pressure at gridbox interfaces
      REAL*8,  INTENT(IN) :: PEDGE(IIPAR, JJPAR, LLPAR+1)
#endif
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
!  06 Mar 2012 - R. Yantosca - Now call GET_OVERHEAD_O3 to get the total
!                              overhead O3 column for FAST-J
!  30 Jul 2012 - R. Yantosca - Now accept am_I_Root as an argument when
!                              running with the traditional driver main.F
!  10 Aug 2012 - R. Yantosca - Replace IPAR, JPAR, LPAR w/ IIPAR, JJPAR, LLPAR
!  27 Mar 2013 - S. D. Eastham - Rolled into FAST_JX_MOD
!  29 Mar 2013 - S. D. Eastham - Updated to Fast-JX from Fast-J
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER, SAVE      :: LASTMONTH = -1
      INTEGER            :: NLON, NLAT, DAY,  MONTH, DAY_OF_YR, L
      REAL*8             :: CSZA, PRES, SFCA, YLAT,  O3COL
      REAL*8             :: TEMP(LLPAR), OPTD(LLPAR)
      REAL*8             :: OPTDUST(LLPAR,NDUST)
      REAL*8             :: OPTAER(LLPAR,NAER*NRH)

      ! Local variables for cloud overlap (hyl, phs)
      INTEGER           :: NUMB, KK, I
      INTEGER            :: INDIC(LLPAR+1)
!#if defined( DEVEL )
!      INTEGER, ALLOCATABLE :: INDGEN(:)
!      INTEGER              :: AS
!#else
      INTEGER            :: INDGEN(LLPAR+1)! = (/ (i,i=1,LLPAR+1) /)
!#endif
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
      ! FAST_JX begins here!
      !=================================================================

#if defined( DEVEL ) || defined( EXTERNAL_GRID ) || defined( EXTERNAL_FORCING )
      IF ( FIRST ) THEN

!         ALLOCATE ( INDGEN(LLPAR+1), STAT=AS)
!         IF ( AS /=0 ) CALL ALLOC_ERR( 'INDGEN' )

         DO i = 1,LLPAR+1
            INDGEN(i) = i       !(/(i,i=1,LLPAR+1)/)
         ENDDO
         FIRST = .FALSE.

      ENDIF
#endif

      ! Get day of year (0-365 or 0-366)
      DAY_OF_YR = GET_DAY_OF_YEAR()

      ! Get current month
      MONTH = GET_MONTH()

      ! Get day of month
      DAY       = GET_DAY()

      !=================================================================
      ! For each (NLON,NLAT) location, call subroutine PHOTO_JX (in a 
      ! parallel loop to compute J-values for the entire column.  
      ! J-values will be stored in the common-block variable ZPJ, and 
      ! will be later accessed via function FJXFUNC. 
      !=================================================================
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( NLON,  NLAT,   YLAT,  CSZA,      OPTAER  )
!$OMP+PRIVATE( PRES,  TEMP,   OPTD,  SFCA,      OPTDUST )
!$OMP+PRIVATE( FMAX,  CLDF1D, KK,    NUMB,      L       )
!$OMP+PRIVATE( KBOT,  KTOP,   ODNEW, INDICATOR, INDIC   )
!$OMP+PRIVATE( O3COL                                    )
!$OMP+SCHEDULE( DYNAMIC )

      ! Loop over latitudes
      DO NLAT = 1, JJPAR

         ! Loop over longitudes
         DO NLON = 1, IIPAR

            ! Grid box latitude [degrees]
            YLAT         = GET_YMID( NLON, NLAT, 1 )

            ! Cosine of solar zenith angle [unitless] at (NLON,NLAT)
            CSZA         = SUNCOS( (NLAT-1)*IIPAR + NLON ) 

#if defined( EXTERNAL_GRID ) || defined( EXTERNAL_FORCING )
            ! Define the PJ array here, using GCM pressures
            DO L = 1, LLPAR+1
               PJ(L)     = PEDGE( NLON, NLAT, L )
            ENDDO
#else
            ! Define the PJ array here, using Ap and Bp from GEOS-Chem
            DO L = 1, NB
               PJ(L)     = GET_PEDGE( NLON, NLAT, L )
            ENDDO
#endif

            ! Top edge of PJ is top of atmosphere (bmy, 2/13/07)
            PJ(NB+1)     = 0d0

            ! Temperature profile [K] at (NLON,NLAT)
            TEMP         = T(NLON,NLAT,1:LLPAR)

            ! Surface albedo [unitless] at (NLON,NLAT)
            SFCA         = ALBD(NLON,NLAT)

            ! Overhead ozone column [DU] at (NLON, NLAT)
            O3COL        = GET_OVERHEAD_O3( NLON, NLAT )

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

               ! Call FAST-JX routines to compute J-values
               CALL PHOTO_JX(NLON,NLAT,SUNCOS,DAY_OF_YR,
     &                       DAY,SFCA,.FALSE.,PJ,--ZZZ--,TEMP,
     &                   --DDD--,--RRR--,  O3COL,  --CLDWP--,
     &                       AERSP,NDXCLD,   NDXAER,   L1U,
     &                       ANU,  VALJXX,   NJXU,     am_I_Root)
!               CALL PHOTOJ( NLON,  NLAT,     YLAT,    DAY_OF_YR,  
!     &                      MONTH, DAY,      CSZA,    TEMP,    
!     &                      SFCA,  OPTD,     OPTDUST, OPTAER, 
!     &                      O3COL, am_I_Root                     )

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
               ! Call FAST-JX routines to compute J-values
               CALL PHOTO_JX(U0,   SZA,      REFLB,    SOLF,  
     &                       LPRTJ,PPP,      ZZZ,      TTT,
     &                       DDD,  RRR,      OOO,      CLDWP,
     &                       AERSP,NDXCLD,   NDXAER,   L1U,
     &                       ANU,  VALJXX,   NJXU,     am_I_Root)
!               CALL PHOTOJ( NLON,  NLAT,     YLAT,    DAY_OF_YR,  
!     &                      MONTH, DAY,      CSZA,    TEMP,  
!     &                      SFCA,  OPTD,     OPTDUST, OPTAER,
!     &                      O3COL, am_I_Root                    )

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
               ! Call FAST-JX routines to compute J-values
                     CALL PHOTO_JX(U0,   SZA,      REFLB,    SOLF,  
     &                             LPRTJ,PPP,      ZZZ,      TTT,
     &                             DDD,  RRR,      OOO,      CLDWP,
     &                            AERSP,NDXCLD,   NDXAER,   L1U,
     &                            ANU,  VALJXX,   NJXU,     am_I_Root)
!                     CALL PHOTOJ( NLON,  NLAT,     YLAT,    DAY_OF_YR, 
!     &                            MONTH, DAY,      CSZA,    TEMP,  
!     &                            SFCA,  OPTD,     OPTDUST, OPTAER,
!     &                            O3COL, am_I_Root                    )

                  CASE( 1:6 ) 
                     CALL MMRAN_16( NUMB,  NLON,  NLAT,      YLAT,   
     &                              DAY,   MONTH, DAY_OF_YR, CSZA,    
     &                              TEMP,  SFCA,  OPTDUST,   OPTAER, 
     &                              LLPAR, FMAX,  ODNEW,     KBOT,   
     &                              KTOP,  O3COL, am_I_Root )

               END SELECT
            ENDIF 
         ENDDO
      ENDDO
!$OMP END PARALLEL DO

      !-----------------------------------------------------------
      ! END OF SUBROUTINE FAST-J
      !-----------------------------------------------------------
      END SUBROUTINE FAST_JX
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: blkslv
!
! !DESCRIPTION: Subroutine BLKSLV solves the block tri-diagonal system
!  $A(I)\timesX(I-1) + B(I)\timesX(I) + C(I)\timesX(I+1) = H(I)$.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE BLKSLV
     &     (FJ,POMEGA,FZ,ZTAU,ZFLUX,RFL,PM,PM0,FJTOP,FJBOT,ND)
!
! !USES:
!
!
! !INPUT PARAMETERS:
!
      INTEGER, INTENT(IN) ::  ND
      REAL*8, INTENT(IN)  ::  POMEGA(M2_,N_,W_),FZ(N_,W_),ZTAU(N_,W_)
     &                       ,PM(M_,M2_),PM0(M2_)
     &                       ,RFL(W_),ZFLUX(W_)
!
! !OUTPUT PARAMETERS:
!
      REAL*8, INTENT(OUT) ::  FJ(N_,W_),FJTOP(W_),FJBOT(W_)
!
! !REMARKS:
! 
! !REVISION HISTORY: 
!  27 Mar 2013 - S. D. Eastham - Copied from GEOS-Chem v9-01-03
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      REAL*8, DIMENSION(M_,N_,W_)    ::  A,C,H,   RR

      REAL*8, DIMENSION(M_,M_,N_,W_) ::  B,AA,CC,  DD
      REAL*8, DIMENSION(M_,M_) ::  E
      REAL*8  SUMB,SUMBX,SUMT
      INTEGER I, J, K, L

      !=================================================================
      ! BLKSLV begins here!
      !=================================================================

      do K = 1,W_
       call GEN_ID (POMEGA(1,1,K),FZ(1,K),ZTAU(1,K),ZFLUX(K),RFL(K),
     &     PM,PM0, B(1,1,1,K),CC(1,1,1,K),AA(1,1,1,K),
     &             A(1,1,K),H(1,1,K),C(1,1,K), ND)
      enddo

      do K = 1,W_
C-----------UPPER BOUNDARY L=1
       L = 1
        do J = 1,M_
         do I = 1,M_
          E(I,J) = B(I,J,1,K)
         enddo
        enddo

c---setup L & U matrices
         E(2,1) = E(2,1)/E(1,1)
         E(2,2) = E(2,2)-E(2,1)*E(1,2)
         E(2,3) = E(2,3)-E(2,1)*E(1,3)
         E(2,4) = E(2,4)-E(2,1)*E(1,4)
         E(3,1) = E(3,1)/E(1,1)
         E(3,2) = (E(3,2)-E(3,1)*E(1,2))/E(2,2)
         E(3,3) = E(3,3)-E(3,1)*E(1,3)-E(3,2)*E(2,3)
         E(3,4) = E(3,4)-E(3,1)*E(1,4)-E(3,2)*E(2,4)
         E(4,1) = E(4,1)/E(1,1)
         E(4,2) = (E(4,2)-E(4,1)*E(1,2))/E(2,2)
         E(4,3) = (E(4,3)-E(4,1)*E(1,3)-E(4,2)*E(2,3))/E(3,3)
         E(4,4) = E(4,4)-E(4,1)*E(1,4)-E(4,2)*E(2,4)-E(4,3)*E(3,4)
c---invert L
         E(4,3) = -E(4,3)
         E(4,2) = -E(4,2)-E(4,3)*E(3,2)
         E(4,1) = -E(4,1)-E(4,2)*E(2,1)-E(4,3)*E(3,1)
         E(3,2) = -E(3,2)
         E(3,1) = -E(3,1)-E(3,2)*E(2,1)
         E(2,1) = -E(2,1)
c---invert U
         E(4,4) = 1.d0/E(4,4)
         E(3,4) = -E(3,4)*E(4,4)/E(3,3)
         E(3,3) = 1.d0/E(3,3)
         E(2,4) = -(E(2,3)*E(3,4)+E(2,4)*E(4,4))/E(2,2)
         E(2,3) = -E(2,3)*E(3,3)/E(2,2)
         E(2,2) = 1.d0/E(2,2)
         E(1,4) = -(E(1,2)*E(2,4)+E(1,3)*E(3,4)+E(1,4)*E(4,4))/E(1,1)
         E(1,3) = -(E(1,2)*E(2,3)+E(1,3)*E(3,3))/E(1,1)
         E(1,2) = -E(1,2)*E(2,2)/E(1,1)
         E(1,1) = 1.d0/E(1,1)
c---multiply U-invers * L-inverse
         E(1,1) = E(1,1)+E(1,2)*E(2,1)+E(1,3)*E(3,1)+E(1,4)*E(4,1)
         E(1,2) = E(1,2)+E(1,3)*E(3,2)+E(1,4)*E(4,2)
         E(1,3) = E(1,3)+E(1,4)*E(4,3)
         E(2,1) = E(2,2)*E(2,1)+E(2,3)*E(3,1)+E(2,4)*E(4,1)
         E(2,2) = E(2,2)+E(2,3)*E(3,2)+E(2,4)*E(4,2)
         E(2,3) = E(2,3)+E(2,4)*E(4,3)
         E(3,1) = E(3,3)*E(3,1)+E(3,4)*E(4,1)
         E(3,2) = E(3,3)*E(3,2)+E(3,4)*E(4,2)
         E(3,3) = E(3,3)+E(3,4)*E(4,3)
         E(4,1) = E(4,4)*E(4,1)
         E(4,2) = E(4,4)*E(4,2)
         E(4,3) = E(4,4)*E(4,3)

        do J = 1,M_
         do I = 1,M_
          DD(I,J,1,K) = -E(I,1)*CC(1,J,1,K)-E(I,2)*CC(2,J,1,K)
     &                  -E(I,3)*CC(3,J,1,K)-E(I,4)*CC(4,J,1,K)
         enddo
          RR(J,1,K) = E(J,1)*H(1,1,K)+E(J,2)*H(2,1,K)
     &              +E(J,3)*H(3,1,K)+E(J,4)*H(4,1,K)
        enddo

C----------CONTINUE THROUGH ALL DEPTH POINTS ID=2 TO ID=ND-1
       do L = 2,ND-1

        do J = 1,M_
         do I = 1,M_
          B(I,J,L,K) = B(I,J,L,K) + A(I,L,K)*DD(I,J,L-1,K)
         enddo
          H(J,L,K) = H(J,L,K) - A(J,L,K)*RR(J,L-1,K)
        enddo

        do J = 1,M_
         do I = 1,M_
          E(I,J) = B(I,J,L,K)
         enddo
        enddo

c---setup L & U matrices
         E(2,1) = E(2,1)/E(1,1)
         E(2,2) = E(2,2)-E(2,1)*E(1,2)
         E(2,3) = E(2,3)-E(2,1)*E(1,3)
         E(2,4) = E(2,4)-E(2,1)*E(1,4)
         E(3,1) = E(3,1)/E(1,1)
         E(3,2) = (E(3,2)-E(3,1)*E(1,2))/E(2,2)
         E(3,3) = E(3,3)-E(3,1)*E(1,3)-E(3,2)*E(2,3)
         E(3,4) = E(3,4)-E(3,1)*E(1,4)-E(3,2)*E(2,4)
         E(4,1) = E(4,1)/E(1,1)
         E(4,2) = (E(4,2)-E(4,1)*E(1,2))/E(2,2)
         E(4,3) = (E(4,3)-E(4,1)*E(1,3)-E(4,2)*E(2,3))/E(3,3)
         E(4,4) = E(4,4)-E(4,1)*E(1,4)-E(4,2)*E(2,4)-E(4,3)*E(3,4)
c---invert L
         E(4,3) = -E(4,3)
         E(4,2) = -E(4,2)-E(4,3)*E(3,2)
         E(4,1) = -E(4,1)-E(4,2)*E(2,1)-E(4,3)*E(3,1)
         E(3,2) = -E(3,2)
         E(3,1) = -E(3,1)-E(3,2)*E(2,1)
         E(2,1) = -E(2,1)
c---invert U
         E(4,4) = 1.d0/E(4,4)
         E(3,4) = -E(3,4)*E(4,4)/E(3,3)
         E(3,3) = 1.d0/E(3,3)
         E(2,4) = -(E(2,3)*E(3,4)+E(2,4)*E(4,4))/E(2,2)
         E(2,3) = -E(2,3)*E(3,3)/E(2,2)
         E(2,2) = 1.d0/E(2,2)
         E(1,4) = -(E(1,2)*E(2,4)+E(1,3)*E(3,4)+E(1,4)*E(4,4))/E(1,1)
         E(1,3) = -(E(1,2)*E(2,3)+E(1,3)*E(3,3))/E(1,1)
         E(1,2) = -E(1,2)*E(2,2)/E(1,1)
         E(1,1) = 1.d0/E(1,1)
c---multiply U-invers * L-inverse
         E(1,1) = E(1,1)+E(1,2)*E(2,1)+E(1,3)*E(3,1)+E(1,4)*E(4,1)
         E(1,2) = E(1,2)+E(1,3)*E(3,2)+E(1,4)*E(4,2)
         E(1,3) = E(1,3)+E(1,4)*E(4,3)
         E(2,1) = E(2,2)*E(2,1)+E(2,3)*E(3,1)+E(2,4)*E(4,1)
         E(2,2) = E(2,2)+E(2,3)*E(3,2)+E(2,4)*E(4,2)
         E(2,3) = E(2,3)+E(2,4)*E(4,3)
         E(3,1) = E(3,3)*E(3,1)+E(3,4)*E(4,1)
         E(3,2) = E(3,3)*E(3,2)+E(3,4)*E(4,2)
         E(3,3) = E(3,3)+E(3,4)*E(4,3)
         E(4,1) = E(4,4)*E(4,1)
         E(4,2) = E(4,4)*E(4,2)
         E(4,3) = E(4,4)*E(4,3)

        do J = 1,M_
         do I = 1,M_
          DD(I,J,L,K) = - E(I,J)*C(J,L,K)
         enddo
          RR(J,L,K) = E(J,1)*H(1,L,K)+E(J,2)*H(2,L,K)
     &              + E(J,3)*H(3,L,K)+E(J,4)*H(4,L,K)
        enddo

       enddo

C---------FINAL DEPTH POINT: L=ND
       L = ND
        do J = 1,M_
         do I = 1,M_
          B(I,J,L,K) = B(I,J,L,K)
     &     + AA(I,1,L,K)*DD(1,J,L-1,K) + AA(I,2,L,K)*DD(2,J,L-1,K)
     &     + AA(I,3,L,K)*DD(3,J,L-1,K) + AA(I,4,L,K)*DD(4,J,L-1,K)
         enddo
          H(J,L,K) = H(J,L,K)
     &     - AA(J,1,L,K)*RR(1,L-1,K) - AA(J,2,L,K)*RR(2,L-1,K)
     &     - AA(J,3,L,K)*RR(3,L-1,K) - AA(J,4,L,K)*RR(4,L-1,K)
        enddo

        do J = 1,M_
         do I = 1,M_
          E(I,J) = B(I,J,L,K)
         enddo
        enddo

c---setup L & U matrices
         E(2,1) = E(2,1)/E(1,1)
         E(2,2) = E(2,2)-E(2,1)*E(1,2)
         E(2,3) = E(2,3)-E(2,1)*E(1,3)
         E(2,4) = E(2,4)-E(2,1)*E(1,4)
         E(3,1) = E(3,1)/E(1,1)
         E(3,2) = (E(3,2)-E(3,1)*E(1,2))/E(2,2)
         E(3,3) = E(3,3)-E(3,1)*E(1,3)-E(3,2)*E(2,3)
         E(3,4) = E(3,4)-E(3,1)*E(1,4)-E(3,2)*E(2,4)
         E(4,1) = E(4,1)/E(1,1)
         E(4,2) = (E(4,2)-E(4,1)*E(1,2))/E(2,2)
         E(4,3) = (E(4,3)-E(4,1)*E(1,3)-E(4,2)*E(2,3))/E(3,3)
         E(4,4) = E(4,4)-E(4,1)*E(1,4)-E(4,2)*E(2,4)-E(4,3)*E(3,4)
c---invert L
         E(4,3) = -E(4,3)
         E(4,2) = -E(4,2)-E(4,3)*E(3,2)
         E(4,1) = -E(4,1)-E(4,2)*E(2,1)-E(4,3)*E(3,1)
         E(3,2) = -E(3,2)
         E(3,1) = -E(3,1)-E(3,2)*E(2,1)
         E(2,1) = -E(2,1)
c---invert U
         E(4,4) = 1.d0/E(4,4)
         E(3,4) = -E(3,4)*E(4,4)/E(3,3)
         E(3,3) = 1.d0/E(3,3)
         E(2,4) = -(E(2,3)*E(3,4)+E(2,4)*E(4,4))/E(2,2)
         E(2,3) = -E(2,3)*E(3,3)/E(2,2)
         E(2,2) = 1.d0/E(2,2)
         E(1,4) = -(E(1,2)*E(2,4)+E(1,3)*E(3,4)+E(1,4)*E(4,4))/E(1,1)
         E(1,3) = -(E(1,2)*E(2,3)+E(1,3)*E(3,3))/E(1,1)
         E(1,2) = -E(1,2)*E(2,2)/E(1,1)
         E(1,1) = 1.d0/E(1,1)
c---multiply U-invers * L-inverse
         E(1,1) = E(1,1)+E(1,2)*E(2,1)+E(1,3)*E(3,1)+E(1,4)*E(4,1)
         E(1,2) = E(1,2)+E(1,3)*E(3,2)+E(1,4)*E(4,2)
         E(1,3) = E(1,3)+E(1,4)*E(4,3)
         E(2,1) = E(2,2)*E(2,1)+E(2,3)*E(3,1)+E(2,4)*E(4,1)
         E(2,2) = E(2,2)+E(2,3)*E(3,2)+E(2,4)*E(4,2)
         E(2,3) = E(2,3)+E(2,4)*E(4,3)
         E(3,1) = E(3,3)*E(3,1)+E(3,4)*E(4,1)
         E(3,2) = E(3,3)*E(3,2)+E(3,4)*E(4,2)
         E(3,3) = E(3,3)+E(3,4)*E(4,3)
         E(4,1) = E(4,4)*E(4,1)
         E(4,2) = E(4,4)*E(4,2)
         E(4,3) = E(4,4)*E(4,3)

        do J = 1,M_
         RR(J,L,K) = E(J,1)*H(1,L,K)+E(J,2)*H(2,L,K)
     &              +E(J,3)*H(3,L,K)+E(J,4)*H(4,L,K)
        enddo

C-----------BACK SOLUTION
       do L = ND-1,1,-1
        do J = 1,M_
         RR(J,L,K) = RR(J,L,K)
     &    + DD(J,1,L,K)*RR(1,L+1,K) + DD(J,2,L,K)*RR(2,L+1,K)
     &    + DD(J,3,L,K)*RR(3,L+1,K) + DD(J,4,L,K)*RR(4,L+1,K)
        enddo
       enddo

C----------mean J & H
       do L = 1,ND,2
        FJ(L,K) = RR(1,L,K)*WT(1) + RR(2,L,K)*WT(2)
     &          + RR(3,L,K)*WT(3) + RR(4,L,K)*WT(4)
       enddo
       do L = 2,ND,2
        FJ(L,K) = RR(1,L,K)*WT(1)*EMU(1) + RR(2,L,K)*WT(2)*EMU(2)
     &          + RR(3,L,K)*WT(3)*EMU(3) + RR(4,L,K)*WT(4)*EMU(4)
       enddo

c---FJTOP = scaled diffuse flux out top-of-atmosphere (limit = mu0)
c---FJBOT = scaled diffuse flux onto surface:
c---ZFLUX = reflect/(1 + reflect) * mu0 * Fsolar(lower boundary)
c---SUMBX = flux from Lambert reflected I+
       SUMT = RR(1, 1,K)*WT(1)*EMU(1) + RR(2, 1,K)*WT(2)*EMU(2)
     &      + RR(3, 1,K)*WT(3)*EMU(3) + RR(4, 1,K)*WT(4)*EMU(4)
       SUMB = RR(1,ND,K)*WT(1)*EMU(1) + RR(2,ND,K)*WT(2)*EMU(2)
     &      + RR(3,ND,K)*WT(3)*EMU(3) + RR(4,ND,K)*WT(4)*EMU(4)
       SUMBX = 4.d0*SUMB*RFL(K)/(1.0d0 + RFL(K)) + ZFLUX(K)

       FJTOP(K) = 4.d0*SUMT
       FJBOT(K) = 4.d0*SUMB - SUMBX

      enddo

      END SUBROUTINE BLKSLV
!
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: gen_id
!
! !DESCRIPTION: Subroutine GEN generates coefficient matrices for the block 
!  tri-diagonal system described in BLKSLV.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE GEN_ID(POMEGA,FZ,ZTAU,ZFLUX,RFL,PM,PM0
     &              ,B,CC,AA,A,H,C,  ND)
!
! !INPUT PARAMETERS:
!
      INTEGER, INTENT(IN) ::  ND
      REAL*8, INTENT(IN)  ::  POMEGA(M2_,N_),PM(M_,M2_),PM0(M2_)
      REAL*8, INTENT(IN)  ::  ZFLUX,RFL
      REAL*8, INTENT(IN),DIMENSION(N_) :: FZ,ZTAU
!
! !OUTPUT PARAMETERS:
!
      REAL*8, INTENT(OUT),DIMENSION(M_,M_,N_) ::  B,AA,CC
      REAL*8, INTENT(OUT),DIMENSION(M_,N_) ::  A,C,H

!
! !REVISION HISTORY: 
!  28 Mar 2013 - S. D. Eastham - Copied from Fast-JX v7.0
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER I, J, K, L1,L2,LL
      REAL*8  SUM0, SUM1, SUM2, SUM3
      REAL*8  DELTAU, D1, D2, SURFAC

      REAL*8, DIMENSION(M_,M_) :: S,T,U,V,W
 
      !=================================================================
      ! GEN_ID begins here!
      !=================================================================

C---------upper boundary:  2nd-order terms
       L1 = 1
       L2 = 2
       do I = 1,M_
        SUM0 =
     &   POMEGA(1,L1)*PM(I,1)*PM0(1) + POMEGA(3,L1)*PM(I,3)*PM0(3)
     & + POMEGA(5,L1)*PM(I,5)*PM0(5) + POMEGA(7,L1)*PM(I,7)*PM0(7)
        SUM2 =
     &   POMEGA(1,L2)*PM(I,1)*PM0(1) + POMEGA(3,L2)*PM(I,3)*PM0(3)
     & + POMEGA(5,L2)*PM(I,5)*PM0(5) + POMEGA(7,L2)*PM(I,7)*PM0(7)
        SUM1 =
     &   POMEGA(2,L1)*PM(I,2)*PM0(2) + POMEGA(4,L1)*PM(I,4)*PM0(4)
     & + POMEGA(6,L1)*PM(I,6)*PM0(6) + POMEGA(8,L1)*PM(I,8)*PM0(8)
        SUM3 =
     &   POMEGA(2,L2)*PM(I,2)*PM0(2) + POMEGA(4,L2)*PM(I,4)*PM0(4)
     & + POMEGA(6,L2)*PM(I,6)*PM0(6) + POMEGA(8,L2)*PM(I,8)*PM0(8)
         H(I,L1) = 0.5d0*(SUM0*FZ(L1) + SUM2*FZ(L2))
         A(I,L1) = 0.5d0*(SUM1*FZ(L1) + SUM3*FZ(L2))
       enddo

       do I = 1,M_
        do J = 1,I
         SUM0 =
     &   POMEGA(1,L1)*PM(I,1)*PM(J,1) + POMEGA(3,L1)*PM(I,3)*PM(J,3)
     & + POMEGA(5,L1)*PM(I,5)*PM(J,5) + POMEGA(7,L1)*PM(I,7)*PM(J,7)
         SUM2 =
     &   POMEGA(1,L2)*PM(I,1)*PM(J,1) + POMEGA(3,L2)*PM(I,3)*PM(J,3)
     & + POMEGA(5,L2)*PM(I,5)*PM(J,5) + POMEGA(7,L2)*PM(I,7)*PM(J,7)
         SUM1 =
     &   POMEGA(2,L1)*PM(I,2)*PM(J,2) + POMEGA(4,L1)*PM(I,4)*PM(J,4)
     & + POMEGA(6,L1)*PM(I,6)*PM(J,6) + POMEGA(8,L1)*PM(I,8)*PM(J,8)
         SUM3 =
     &   POMEGA(2,L2)*PM(I,2)*PM(J,2) + POMEGA(4,L2)*PM(I,4)*PM(J,4)
     & + POMEGA(6,L2)*PM(I,6)*PM(J,6) + POMEGA(8,L2)*PM(I,8)*PM(J,8)
         S(I,J) = - SUM2*WT(J)
         S(J,I) = - SUM2*WT(I)
         T(I,J) = - SUM1*WT(J)
         T(J,I) = - SUM1*WT(I)
         V(I,J) = - SUM3*WT(J)
         V(J,I) = - SUM3*WT(I)
         B(I,J,L1) = - 0.5d0*(SUM0 + SUM2)*WT(J)
         B(J,I,L1) = - 0.5d0*(SUM0 + SUM2)*WT(I)
        enddo
       enddo

       do I = 1,M_
         S(I,I)   = S(I,I)   + 1.0d0
         T(I,I)   = T(I,I)   + 1.0d0
         V(I,I)   = V(I,I)   + 1.0d0
         B(I,I,L1)= B(I,I,L1) + 1.0d0

         C(I,L1)= S(I,1)*A(1,L1)/EMU(1) + S(I,2)*A(2,L1)/EMU(2)
     &          + S(I,3)*A(3,L1)/EMU(3) + S(I,4)*A(4,L1)/EMU(4)
       enddo

       do I = 1,M_
        do J = 1,M_
         W(J,I) = S(J,1)*T(1,I)/EMU(1) + S(J,2)*T(2,I)/EMU(2)
     &          + S(J,3)*T(3,I)/EMU(3) + S(J,4)*T(4,I)/EMU(4)
         U(J,I) = S(J,1)*V(1,I)/EMU(1) + S(J,2)*V(2,I)/EMU(2)
     &          + S(J,3)*V(3,I)/EMU(3) + S(J,4)*V(4,I)/EMU(4)
        enddo
       enddo
C-------------upper boundary, 2nd-order, C-matrix is full (CC)
         DELTAU = ZTAU(L2) - ZTAU(L1)
         D2 = 0.25d0*DELTAU
       do I = 1,M_
        do J = 1,M_
         B(I,J,L1) = B(I,J,L1) + D2*W(I,J)
         CC(I,J,L1) = D2*U(I,J)
        enddo
         H(I,L1) = H(I,L1) + 2.0d0*D2*C(I,L1)
         A(I,L1) = 0.0d0
       enddo
       do I = 1,M_
        D1 = EMU(I)/DELTAU
        B(I,I,L1)  = B(I,I,L1) + D1
        CC(I,I,L1) = CC(I,I,L1) - D1
       enddo

c------------intermediate points:  can be even or odd, A & C diagonal
c---mid-layer h-points, Legendre terms 2,4,6,8
       do LL=2,ND-1,2
        DELTAU = ZTAU(LL+1) - ZTAU(LL-1)
        do I = 1,M_
          A(I,LL) = EMU(I)/DELTAU
          C(I,LL) = -A(I,LL)
          H(I,LL) = FZ(LL)*(
     &     POMEGA(2,LL)*PM(I,2)*PM0(2) + POMEGA(4,LL)*PM(I,4)*PM0(4)
     &   + POMEGA(6,LL)*PM(I,6)*PM0(6) + POMEGA(8,LL)*PM(I,8)*PM0(8))
        enddo
        do I = 1,M_
         do J=1,I
          SUM0 =
     &     POMEGA(2,LL)*PM(I,2)*PM(J,2) + POMEGA(4,LL)*PM(I,4)*PM(J,4)
     &    +POMEGA(6,LL)*PM(I,6)*PM(J,6) + POMEGA(8,LL)*PM(I,8)*PM(J,8)
          B(I,J,LL) =  - SUM0*WT(J)
          B(J,I,LL) =  - SUM0*WT(I)
         enddo
        enddo
        do I = 1,M_
          B(I,I,LL) = B(I,I,LL) + 1.0d0
        enddo
       enddo

c---odd-layer j-points, Legendre terms 1,3,5,7
       do LL=3,ND-2,2
        DELTAU = ZTAU(LL+1) - ZTAU(LL-1)
        do I = 1,M_
          A(I,LL) = EMU(I)/DELTAU
          C(I,LL) = -A(I,LL)
          H(I,LL) = FZ(LL)*(
     &     POMEGA(1,LL)*PM(I,1)*PM0(1) + POMEGA(3,LL)*PM(I,3)*PM0(3)
     &   + POMEGA(5,LL)*PM(I,5)*PM0(5) + POMEGA(7,LL)*PM(I,7)*PM0(7))
        enddo
        do I = 1,M_
         do J=1,I
          SUM0 =
     &     POMEGA(1,LL)*PM(I,1)*PM(J,1) + POMEGA(3,LL)*PM(I,3)*PM(J,3)
     &    +POMEGA(5,LL)*PM(I,5)*PM(J,5) + POMEGA(7,LL)*PM(I,7)*PM(J,7)
          B(I,J,LL) =  - SUM0*WT(J)
          B(J,I,LL) =  - SUM0*WT(I)
         enddo
        enddo
        do I = 1,M_
          B(I,I,LL) = B(I,I,LL) + 1.0d0
        enddo
       enddo

C---------lower boundary:  2nd-order terms
       L1 = ND
       L2 = ND-1
       do I = 1,M_
        SUM0 =
     &   POMEGA(1,L1)*PM(I,1)*PM0(1) + POMEGA(3,L1)*PM(I,3)*PM0(3)
     & + POMEGA(5,L1)*PM(I,5)*PM0(5) + POMEGA(7,L1)*PM(I,7)*PM0(7)
        SUM2 =
     &   POMEGA(1,L2)*PM(I,1)*PM0(1) + POMEGA(3,L2)*PM(I,3)*PM0(3)
     & + POMEGA(5,L2)*PM(I,5)*PM0(5) + POMEGA(7,L2)*PM(I,7)*PM0(7)
        SUM1 =
     &   POMEGA(2,L1)*PM(I,2)*PM0(2) + POMEGA(4,L1)*PM(I,4)*PM0(4)
     & + POMEGA(6,L1)*PM(I,6)*PM0(6) + POMEGA(8,L1)*PM(I,8)*PM0(8)
        SUM3 =
     &   POMEGA(2,L2)*PM(I,2)*PM0(2) + POMEGA(4,L2)*PM(I,4)*PM0(4)
     & + POMEGA(6,L2)*PM(I,6)*PM0(6) + POMEGA(8,L2)*PM(I,8)*PM0(8)
         H(I,L1) = 0.5d0*(SUM0*FZ(L1) + SUM2*FZ(L2))
         A(I,L1) = 0.5d0*(SUM1*FZ(L1) + SUM3*FZ(L2))
       enddo

       do I = 1,M_
        do J = 1,I
         SUM0 =
     &    POMEGA(1,L1)*PM(I,1)*PM(J,1) + POMEGA(3,L1)*PM(I,3)*PM(J,3)
     &  + POMEGA(5,L1)*PM(I,5)*PM(J,5) + POMEGA(7,L1)*PM(I,7)*PM(J,7)
         SUM2 =
     &    POMEGA(1,L2)*PM(I,1)*PM(J,1) + POMEGA(3,L2)*PM(I,3)*PM(J,3)
     &  + POMEGA(5,L2)*PM(I,5)*PM(J,5) + POMEGA(7,L2)*PM(I,7)*PM(J,7)
         SUM1 =
     &    POMEGA(2,L1)*PM(I,2)*PM(J,2) + POMEGA(4,L1)*PM(I,4)*PM(J,4)
     &  + POMEGA(6,L1)*PM(I,6)*PM(J,6) + POMEGA(8,L1)*PM(I,8)*PM(J,8)
         SUM3 =
     &    POMEGA(2,L2)*PM(I,2)*PM(J,2) + POMEGA(4,L2)*PM(I,4)*PM(J,4)
     &  + POMEGA(6,L2)*PM(I,6)*PM(J,6) + POMEGA(8,L2)*PM(I,8)*PM(J,8)
         S(I,J) = - SUM2*WT(J)
         S(J,I) = - SUM2*WT(I)
         T(I,J) = - SUM1*WT(J)
         T(J,I) = - SUM1*WT(I)
         V(I,J) = - SUM3*WT(J)
         V(J,I) = - SUM3*WT(I)
         B(I,J,L1) = - 0.5d0*(SUM0 + SUM2)*WT(J)
         B(J,I,L1) = - 0.5d0*(SUM0 + SUM2)*WT(I)
        enddo
       enddo

       do I = 1,M_
         S(I,I)   = S(I,I)   + 1.0d0
         T(I,I)   = T(I,I)   + 1.0d0
         V(I,I)   = V(I,I)   + 1.0d0
         B(I,I,L1)= B(I,I,L1) + 1.0d0

         C(I,L1)= S(I,1)*A(1,L1)/EMU(1) + S(I,2)*A(2,L1)/EMU(2)
     &          + S(I,3)*A(3,L1)/EMU(3) + S(I,4)*A(4,L1)/EMU(4)
       enddo

       do I = 1,M_
        do J = 1,M_
         W(J,I) = S(J,1)*T(1,I)/EMU(1) + S(J,2)*T(2,I)/EMU(2)
     &          + S(J,3)*T(3,I)/EMU(3) + S(J,4)*T(4,I)/EMU(4)
         U(J,I) = S(J,1)*V(1,I)/EMU(1) + S(J,2)*V(2,I)/EMU(2)
     &          + S(J,3)*V(3,I)/EMU(3) + S(J,4)*V(4,I)/EMU(4)
        enddo
       enddo

C------------lower boundary, 2nd-order, A-matrix is full (AA)
         DELTAU = ZTAU(L1) - ZTAU(L2)
         D2 = 0.25d0*DELTAU
         SURFAC = 4.0d0*RFL/(1.0d0 + RFL)
       do I = 1,M_
          D1 = EMU(I)/DELTAU
          SUM0 = D1 + D2*(W(I,1)+W(I,2)+W(I,3)+W(I,4))
          SUM1 = SURFAC*SUM0
        do J = 1,M_
         AA(I,J,L1) = - D2*U(I,J)
         B(I,J,L1) = B(I,J,L1) + D2*W(I,J) - SUM1*EMU(J)*WT(J)
        enddo
         H(I,L1) = H(I,L1) - 2.0d0*D2*C(I,L1) + SUM0*ZFLUX
       enddo

       do I = 1,M_
          D1 = EMU(I)/DELTAU
        AA(I,I,L1) = AA(I,I,L1) + D1
        B(I,I,L1)  = B(I,I,L1) + D1
        C(I,L1) = 0.0d0
       enddo

      END SUBROUTINE GEN_ID
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: opticl
!
! !DESCRIPTION: Subroutine OPTICL calculated cloud properties at 5 wavelengths.
! 01 W_C02   (C1/Deir)GAMMA:r-m=2.0/alf=6 n=1.335   reff=3.000___G=19.55_rho=1.000
! 02 W_C04   (C1/Deir)GAMMA:r-m=4.0/alf=6 n=1.335   reff=6.000___G=78.19_rho=1.000
! 03 W_C08   (C1/Deir)GAMMA:r-m=8.0/alf=2 n=1.335   reff=12.00___G=301.1_rho=1.000
! 04 W_C13   (C1/Deir)GAMMA:r-m=13./alf=2 n=1.335   reff=20.00___G=472.9_rho=1.000
! 05 W_L06   (W/Lacis)GAMMA:r-m=5.5/alf=11/3        reff=10.00___G=183.9_rho=1.000
! 06 Ice-Hexagonal (Mishchencko)                    reff=50.00___G=999.9_rho=0.917
! 07 Ice-Irregular (Mishchencko)                    reff=50.00___G=999.9_rho=0.917
! 08 S-Bkg   LOGN:r=.090 s=.600 n=1.514/.../1.435   reff=0.221___G=.0523_rho=1.630
! 09 S-Vol   LOGN:r=.080 s=.800 n=1.514/.../1.435   reff=0.386___G=.0721_rho=1.630
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE OPTICL (OPTD,SSALB,SLEG, ODCLD,NDCLD)
!
! !USES:
!
! !INPUT PARAMETERS:
!
      real*8, intent(in)::     ODCLD      ! optical depth of cloud layer @600 nm
      integer,intent(inout)::  NDCLD      ! index of cloud layer:  4:13
!
! !OUTPUT VARIABLES:
!
      real*8, intent(out)::    OPTD(5)    ! optical depth of layer
      real*8, intent(out)::    SSALB(5)   ! single-scattering albedo
      real*8, intent(out)::    SLEG(8,5)  ! scatt phase fn (Leg coeffs)
!
! !REMARKS:
! 
! !REVISION HISTORY: 
!  28 Mar 2013 - S. D. Eastham - Copied from Fast-JX v7.0
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      integer I,J
      real*8  XTINCT, REFF,RHO

      !=================================================================
      ! OPTICL begins here!
      !=================================================================

c---default cloud type C1, Reff = 12 microns
      if (NDCLD .lt. 1 .or. NDCLD .gt. 9) then
         NDCLD = 3
      endif

c--rescale OD by Qext at 600 nm (J=4)
      do J=1,5
         OPTD(J) = ODCLD * QCC(J,NDCLD)/QCC(4,NDCLD)
         SSALB(J) = SCC(J,NDCLD)
        do I=1,8
         SLEG(I,J) =  PCC(I,J,NDCLD)
        enddo
      enddo

      END SUBROUTINE OPTICL
!
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: optica
!
! !DESCRIPTION: Subroutine OPTICA calculates aerosol optical properties.       
!---UCI aersols optical data  v-7.0+
! 01 S-Bkg   LOGN:r=.090 s=.600 n=1.514/.../1.435   reff=0.221___G=.0523_rho=1.630
! 02 S-Vol   LOGN:r=.080 s=.800 n=1.514/.../1.435   reff=0.386___G=.0721_rho=1.630
! 03 UT-sulfate LOGN:r=0.05 s=.693 n=1.44           reff=0.166___G=.0205_rho=1.769
! 04 UT-sulfate LOGN:r=0.05 s=.693 n=1.46           reff=0.166___G=.0205_rho=1.769
! 05 UT-sulfatM LOGN:r=.050 s=.642 n=1.53           reff=0.140___G=.0179_rho=1.769
! 06 UM-BC1     LOGN:r=.050 s=.642 n=1.80+0.50i     reff=0.140___G=.0179_rho=1.500
! 07 UM-BC2     LOGN:r=.080 s=.501 n=1.80+0.50i     reff=0.150___G=.0332_rho=1.500
! 08 UM-BB08 (%BC)LOGN:r=.080 s=.500 n=1.552+0.04i  reff=0.149___G=.0331_rho=1.230
! 09 UM-FF04 (%BC) LOGN:r=.050 s=.642 n=1.541+0.02i reff=0.140___G=.0179_rho=1.212
! 10 UM-FF10 (%BC)LOGN:r=.050 s=.642 n=1.557+0.05i  reff=0.140___G=.0179_rho=1.230
! 11 MDust.15  (R.V. Martin generated phase fns)    reff=0.150___G=1.000_rho=2.600
! 12 MDust.25  (R.V. Martin generated phase fns)    reff=0.250___G=1.000_rho=2.600
! 13 MDust.40  (R.V. Martin generated phase fns)    reff=0.400___G=1.000_rho=2.600
! 14 MDust.80  (R.V. Martin generated phase fns)    reff=0.800___G=1.000_rho=2.6001
! 15 MDust1.5  (R.V. Martin generated phase fns)    reff=1.500___G=1.000_rho=2.600
! 16 MDust2.5  (R.V. Martin generated phase fns)    reff=2.500___G=1.000_rho=2.600
! 17 MDust4.0  (R.V. Martin generated phase fns)    reff=4.000___G=1.000_rho=2.600
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE OPTICA (OPTD,SSALB,SLEG,PATH,RELH,K)
!
! !USES:
!
! !INPUT PARAMETERS:
!
      real*8, intent(in)::     PATH       ! path (g/m2) of aerosol/cloud
      real*8, intent(in)::     RELH       ! relative humidity (0.00->1.00+)
      integer,intent(inout)::     K       ! index of cloud/aerosols
!
! !OUTPUT VARIABLES:
!
      real*8, intent(out)::    OPTD(5)    ! optical depth of layer
      real*8, intent(out)::    SSALB(5)   ! single-scattering albedo
      real*8, intent(out)::    SLEG(8,5)  ! scatt phase fn (Leg coeffs)
!
! !REMARKS:
! 
! !REVISION HISTORY: 
!  28 Mar 2013 - S. D. Eastham - Copied from Fast-JX v7.0
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      integer I,J
      real*8  XTINCT, REFF,RHO

      !=================================================================
      ! OPTICA begins here!
      !=================================================================

      if (K .gt. NAA .or. K .lt. 1) then
         write(6,*) ' aerosol index out-of-range: K/NAA',K,NAA
         K = 18
      endif

         REFF = RAA(K)
         RHO = DAA(K)
      do J=1,5
c---extinction K(m2/g) = Q(wvl) / [4/3 * Reff(micron) * aerosol-density(g/cm3)]
         XTINCT = 0.75d0*QAA(J,K)/(REFF*RHO)
         OPTD(J) = PATH*XTINCT
         SSALB(J) = SAA(J,K)
       do I=1,8
         SLEG(I,J) =  PAA(I,J,K)
       enddo
      enddo

      END SUBROUTINE OPTICA
!
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: opticm
!
! !DESCRIPTION: Subroutine OPTICM calculates aerosol optical properties.       
!  For the U. Michigan aerosol data sets, this generate Fast-JX data formats.
!  NB: Approximates the Legendre expansion(L) of the scattering phase fn
!      as (2*L+1)*g**L
!   UMAER(I,J,K,L):
!   I=1:3 = [SSAbldeo, g, k-ext(m2/g)]
!   J=1:5 = [200, 300, 400, (550,) 600 , 1000 nm]
!   K=1:21= [0, 5, 10, 15, ..., 90, 95, 99 %RelHum]
!   L=1:33= UM aerosol types [SULF, SS-1,-2,-3,-4, DD-1,-2,-3,-4, FF00(0%BC),
!                      FF02, ...FF14(14%BC), BB00, BB02, ...BB30(30%BC)]
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE OPTICM (OPTD,SSALB,SLEG,PATH,RELH,LL)
!
! !USES:
!
! !INPUT PARAMETERS:
!
      real*8, intent(in)       :: PATH       ! path (g/m2) of aerosol/cloud
      real*8, intent(in)       :: RELH       ! relative humidity (0.00->1.00+)
      integer,intent(in)       :: LL         ! index of cloud/aerosols
!
! !OUTPUT VARIABLES:
!
      real*8, intent(out)      :: OPTD(5)    ! optical depth of layer
      real*8, intent(out)      :: SSALB(5)   ! single-scattering albedo
      real*8, intent(out)      :: SLEG(8,5)  ! scatt phase fn (Leg coeffs)
!
! !REMARKS:
! 
! !REVISION HISTORY: 
!  28 Mar 2013 - S. D. Eastham - Copied from Fast-JX v7.0
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      integer KR,I,J
      real*8  XTINCT,R,FRH,GCOS

      !=================================================================
      ! OPTICM begins here!
      !=================================================================

c---calculate fast-JX properties at the std 5 wavelengths:200-300-400-600-999nm
c---extrapolate phase fn from first term (g)
      L = LL
      if (L .lt. 1 .or. L .gt. 33) then
cccc         write(6,*) ' UM aer index too large: L',L
         L = 1
      endif

c---pick nearest Relative Humidity
      KR =  20.d0*RELH  + 1.5d0
      KR = max(1, min(21, KR))

      do J=1,5
       SSALB(J) = UMAER(1,J,KR,L)
         XTINCT = UMAER(3,J,KR,L)
       OPTD(J) = PATH*XTINCT
         GCOS   = UMAER(2,J,KR,L)
       SLEG(1,J) =  1.d0
       SLEG(2,J) =  3.d0*GCOS
       SLEG(3,J) =  5.d0*GCOS**2
       SLEG(4,J) =  7.d0*GCOS**3
       SLEG(5,J) =  9.d0*GCOS**4
       SLEG(6,J) = 11.d0*GCOS**5
       SLEG(7,J) = 13.d0*GCOS**6
       SLEG(8,J) = 15.d0*GCOS**7
      enddo

      END SUBROUTINE OPTICM
!
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: jratet
!
! !DESCRIPTION: Subroutine JRATET calculates temperature-dependent J-rates.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE JRATET(PPJ,TTJ,FFF,VALJL,LU,NJXU)
!
! !USES:
!
!
! !INPUT PARAMETERS:
!
      integer, intent(in)  :: LU,NJXU
      real*8, intent(in)  ::  PPJ(LU+1),TTJ(LU+1)
      real*8, intent(inout)  ::  FFF(W_,LU)
!
! !OUTPUT VARIABLES:
!
      real*8, intent(out), dimension(LU,NJXU) ::  VALJL
!
! !REMARKS:
! 
! !REVISION HISTORY: 
!  28 Mar 2013 - S. D. Eastham - Copied from Fast-JX v7.0
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      real*8  VALJ(X_)
      real*8  QO2TOT, QO3TOT, QO31DY, QO31D, QQQT, TFACT
      real*8  TT,PP,DD,TT200,TFACA,TFAC0,TFAC1,TFAC2,QQQA,QQ2,QQ1A,QQ1B
      integer J,K,L, IV

      !=================================================================
      ! JRATET begins here!
      !=================================================================

      if (NJXU .lt. NJX) then
        call EXITC(' JRATET:  CTM has not enough J-values dimensioned')
      endif
      do L = 1,LU
c---need temperature, pressure, and density at mid-layer (for some quantum yields):
          TT   = TTJ(L)
         if (L .eq. 1) then
          PP = PPJ(1)
         else
          PP  = (PPJ(L)+PPJ(L+1))*0.5d0
         endif
          DD = 7.24e18*PP/TT

c---if W_=18/12, must zero bin-11/5 below 100 hPa, since O2 e-fold is too weak
c        and does not represent the decay of 215.5-221.5 nm sunlight.
        if (PP .gt. 100.d0) then
          if (W_ .eq. 18) then
            FFF(11,L) = 0.d0
          elseif (W_ .eq. 12) then
            FFF(5,L) = 0.d0
          endif
        endif

        do J = 1,NJX
          VALJ(J) = 0.d0
        enddo

         do K = 1,W_
          call X_interp (TT,QO2TOT, TQQ(1,1),QO2(K,1),
     &            TQQ(2,1),QO2(K,2), TQQ(3,1),QO2(K,3), LQQ(1))
          call X_interp (TT,QO3TOT, TQQ(1,2),QO3(K,1),
     &            TQQ(2,2),QO3(K,2), TQQ(3,2),QO3(K,3), LQQ(2))
          call X_interp (TT,QO31DY, TQQ(1,3),Q1D(K,1),
     &            TQQ(2,3),Q1D(K,2), TQQ(3,3),Q1D(K,3), LQQ(3))
            QO31D  = QO31DY*QO3TOT
           VALJ(1) = VALJ(1) + QO2TOT*FFF(K,L)
           VALJ(2) = VALJ(2) + QO3TOT*FFF(K,L)
           VALJ(3) = VALJ(3) + QO31D*FFF(K,L)
         enddo

        do J = 4,NJX
         do K = 1,W_
c---also need to allow for Pressure interpolation if SQQ(J) = 'p'
          if (SQQ(J) .eq.'p') then
           call X_interp (PP,QQQT, TQQ(1,J),QQQ(K,1,J),
     &        TQQ(2,J),QQQ(K,2,J), TQQ(3,J),QQQ(K,3,J), LQQ(J))
          else
           call X_interp (TT,QQQT, TQQ(1,J),QQQ(K,1,J),
     &        TQQ(2,J),QQQ(K,2,J), TQQ(3,J),QQQ(K,3,J), LQQ(J))
          endif
            VALJ(J) = VALJ(J) + QQQT*FFF(K,L)
         enddo
        enddo

        do J=1,NJX
          VALJL(L,J) = VALJ(J)
        enddo

      enddo

      END SUBROUTINE JRATET
!
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: x_interp
!
! !DESCRIPTION: Subroutine X\_INTERP is an up-to-three-point linear interp.
!  function for cross-sections.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE X_INTERP (TINT,XINT,T1,X1,T2,X2,T3,X3,L123)
!
! !USES:
!
!
! !INPUT PARAMETERS:
!
      REAL*8, INTENT(IN)::  TINT,T1,T2,T3, X1,X2,X3
      INTEGER,INTENT(IN)::  L123
!
! !OUTPUT VARIABLES:
!
      REAL*8,INTENT(OUT)::  XINT
!
! !REMARKS:
! 
! !REVISION HISTORY: 
!  28 Mar 2013 - S. D. Eastham - Copied from Fast-JX v7.0
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      REAL*8  TFACT

      !=================================================================
      ! X_INTERP begins here!
      !=================================================================

      if (L123 .le. 1) then
           XINT = X1
      elseif (L123 .eq. 2) then
             TFACT = max(0.d0,min(1.d0,(TINT-T1)/(T2-T1) ))
           XINT = X1 + TFACT*(X2 - X1)
      else
        if (TINT.le. T2) then
             TFACT = max(0.d0,min(1.d0,(TINT-T1)/(T2-T1) ))
           XINT = X1 + TFACT*(X2 - X1)
        else

             TFACT = max(0.d0,min(1.d0,(TINT-T2)/(T3-T2) ))
           XINT = X2 + TFACT*(X3 - X2)
        endif
      endif

      END SUBROUTINE X_INTERP
!
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: jp_atm
!
! !DESCRIPTION: Subroutine JP\_ATM sets up the atmosphere for photolysis.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE JP_ATM(PPJ,TTJ,DDJ,OOJ,ZZJ,DTAU6,POMEG6,JXTRA,LU)
!
! !USES:
!
!
! !INPUT PARAMETERS:
!
      INTEGER,INTENT(IN)                  :: LU
      REAL*8, INTENT(IN), DIMENSION(LU+2) :: PPJ,ZZJ
      REAL*8, INTENT(IN), DIMENSION(LU+1) :: TTJ,DDJ,OOJ,DTAU6
      REAL*8, INTENT(IN), DIMENSION(8,LU+1) :: POMEG6
      INTEGER,INTENT(IN), DIMENSION(LU+LU+3) :: JXTRA
!
! !OUTPUT VARIABLES:
!
!
! !REMARKS:
! 
! !REVISION HISTORY: 
!  28 Mar 2013 - S. D. Eastham - Copied from Fast-JX v7.0
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER  I,J,K,L
      REAL*8   COLO2,COLO3,ZKM,DELZ,ZTOP

      !=================================================================
      ! JP_ATM begins here!
      !=================================================================

      write(6,'(4a)') '   L z(km)     p      T   ',
     & '    d(air)   d(O3)','  col(O2)  col(O3)     d-TAU   SS-alb',
     & '  g(cos) CTM lyr=>'

      L = LU+2
      write(6,'(1x,i3,0p,f6.2,f10.3,f7.2,1p,4e9.2,0p,f10.4,2f8.5,2i3)')
     &      L,ZZJ(L)*1.d-5,PPJ(L)

          COLO2 = 0.d0
          COLO3 = 0.d0
          ZTOP = ZZJ(LU+2)

        do L = LU+1,1,-1
          COLO2 = COLO2 + DDJ(L)*0.20948d0
          COLO3 = COLO3 + OOJ(L)
          DELZ = ZTOP-ZZJ(L)
          ZTOP = ZZJ(L)
          ZKM = ZZJ(L)*1.d-5

      write(6,'(1x,i3,0p,f6.2,f10.3,f7.2,1p,4e9.2,0p,f10.4,2f8.5,2i3)')
     &      L,ZKM,PPJ(L),TTJ(L),DDJ(L)/DELZ,OOJ(L)/DELZ,
     &      COLO2,COLO3,DTAU6(L),POMEG6(1,L),POMEG6(2,L)/3.d0,
     &      JXTRA(L+L),JXTRA(L+L-1)

        enddo

      END SUBROUTINE JP_ATM
!
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: jp_atm0
!
! !DESCRIPTION: Subroutine JP\_ATM0 sets up the atmosphere for photolysis.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE JP_ATM0(PPJ,TTJ,DDJ,OOJ,ZZJ,LU)
!
! !USES:
!
!
! !INPUT PARAMETERS:
!
      INTEGER,INTENT(IN)                  :: LU
      REAL*8, INTENT(IN), DIMENSION(LU+2) :: PPJ,ZZJ
      REAL*8, INTENT(IN), DIMENSION(LU+1) :: TTJ,DDJ,OOJ
!
! !OUTPUT VARIABLES:
!
!
! !REMARKS:
! 
! !REVISION HISTORY: 
!  28 Mar 2013 - S. D. Eastham - Copied from Fast-JX v7.0
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER  I,J,K,L
      REAL*8   COLO2,COLO3,ZKM,DELZ,ZTOP

      !=================================================================
      ! JP_ATM0 begins here!
      !=================================================================

      write(6,'(4a)') '   L z(km)     p      T   ',
     & '    d(air)   d(O3)','  col(O2)  col(O3)     d-TAU   SS-alb',
     & '  g(cos) CTM lyr=>'
      L = LU+2
      write(6,'(1x,i3,0p,f6.2,f10.3,f7.2,1p,4e9.2,0p,f10.4,2f8.5,2i3)')
     &      L,ZZJ(L)*1.d-5,PPJ(L)
          COLO2 = 0.d0
          COLO3 = 0.d0
          ZTOP = ZZJ(LU+2)
        do L = LU+1,1,-1
          COLO2 = COLO2 + DDJ(L)*0.20948d0
          COLO3 = COLO3 + OOJ(L)
          DELZ = ZTOP-ZZJ(L)
          ZTOP = ZZJ(L)
          ZKM = ZZJ(L)*1.d-5
      write(6,'(1x,i3,0p,f6.2,f10.3,f7.2,1p,4e9.2,0p,f10.4,2f8.5,2i3)')
     &      L,ZKM,PPJ(L),TTJ(L),DDJ(L)/DELZ,OOJ(L)/DELZ,
     &      COLO2,COLO3
        enddo

      END SUBROUTINE JP_ATM0
!
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: sphere2
!
! !DESCRIPTION: Subroutine SPHERE2 is an AMF2. Quoting from the original:
!  New v6.2: does AirMassFactors for mid-layer, needed for SZA ~ 90
!  This new AMF2 does each of the half-layers of the CTM separately,
!     whereas the original, based on the pratmo code did the whole layers
!     and thus calculated the ray-path to the CTM layre edges, NOT the middle.
!  Since fast-JX is meant to calculate the intensity at the mid-layer, the
!     solar beam at low sun (interpolated between layer edges) was incorrect.
!  This new model does make some approximations of the geometry of the layers:
!     the CTM layer is split evenly in mass (good) and in height (approx).
!
!  Calculation of spherical geometry; derive tangent heights, slant path
!  lengths and air mass factor for each layer. Not called when
!  SZA > 98 degrees.  Beyond 90 degrees, include treatment of emergent
!  beam (where tangent height is below altitude J-value desired at).
!  ---------------------------------------------------------------------
!  Inputs:
!     U0      cos(solar zenith angle)
!     RAD     radius of Earth mean sea level (cm)
!     ZHL(L)  height (cm) of the bottome edge of CTM level L
!     L1U     dimension of CTM = levels +1 (L+1 = above-CTM level)
!  Outputs:
!     AMF2(I,J) = air mass factor for CTM level I for sunlight reaching J
!         ( these are calculated for both layer middle and layer edge)
!  ---------------------------------------------------------------------
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE SPHERE2 (U0,RAD,ZHL,AMF2,L1U,LJX1U)
!
! !USES:
!
!
! !INPUT PARAMETERS:
!
      INTEGER, INTENT(IN) ::   L1U, LJX1U
      REAL*8, INTENT(IN)  ::   U0,RAD,ZHL(L1U+1)
!
! !OUTPUT VARIABLES:
!
      REAL*8, INTENT(OUT) ::   AMF2(2*LJX1U+1,2*LJX1U+1)
!
! !REMARKS:
! 
! !REVISION HISTORY: 
!  28 Mar 2013 - S. D. Eastham - Copied from Fast-JX v7.0
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER, PARAMETER  ::  LSPH_ = 100

!     RZ      Distance from centre of Earth to each point (cm)
!     RQ      Square of radius ratios
!     SHADHT  Shadow height for the current SZA
!     XL      Slant path between points

      INTEGER  I, J, K, II, L2
      REAL*8   XMU1,XMU2,XL,DIFF,SHADHT,RZ(LSPH_+1)
      REAL*8   RZ2(2*LSPH_+1),RQ2(2*LSPH_+1)

      !=================================================================
      ! SPHERE2 begins here!
      !=================================================================

c--- must have top-of-atmos (NOT top-of-CTM) defined

      if (L1U .gt. LSPH_) then
        call EXITC(' SPHERE2: temp arrays not large enough')
      endif

        RZ(1) = RAD + ZHL(1)
      do II = 2,L1U+1
        RZ(II)   = RAD + ZHL(II)
      enddo

c---calculate heights for edges of split CTM-layers
      L2 = 2*L1U
      do II = 2,L2,2
        I = II/2
        RZ2(II-1) = RZ(I)
        RZ2(II) = 0.5d0*(RZ(I)+RZ(I+1))
      enddo
        RZ2(L2+1) = RZ(L1U+1)
      do II = 1,L2
        RQ2(II) = (RZ2(II)/RZ2(II+1))**2
      enddo

c---shadow height for SZA > 90
      if (U0 .lt. 0.0d0)  then
        SHADHT = RZ2(1)/dsqrt(1.0d0 - U0**2)
      else
        SHADHT = 0.d0
      endif

c---up from the surface calculating the slant paths between each level
c---  and the level above, and deriving the appropriate Air Mass Factor
         AMF2(:,:) = 0.d0

      do 16 J = 1,2*L1U+1

c  Air Mass Factors all zero if below the tangent height
        if (RZ2(J) .lt. SHADHT) goto 16
c  Ascend from layer J calculating AMF2s
        XMU1 = abs(U0)
        do I = J,2*L1U
          XMU2     = dsqrt(1.0d0 - RQ2(I)*(1.0d0-XMU1**2))
          XL       = RZ2(I+1)*XMU2 - RZ2(I)*XMU1
          AMF2(I,J) = XL / (RZ2(I+1)-RZ2(I))
          XMU1     = XMU2
        enddo
c--fix above top-of-atmos (L=L1U+1), must set DTAU(L1U+1)=0
          AMF2(2*L1U+1,J) = 1.d0
c
c  Twilight case - Emergent Beam, calc air mass factors below layer
        if (U0 .ge. 0.0d0) goto 16

c  Descend from layer J
          XMU1       = abs(U0)
         do II = J-1,1,-1
          DIFF        = RZ2(II+1)*sqrt(1.0d0-XMU1**2)-RZ2(II)
          if (II.eq.1)  DIFF = max(DIFF,0.d0)   ! filter
c  Tangent height below current level - beam passes through twice
          if (DIFF .lt. 0.0d0)  then
            XMU2      = sqrt(1.0d0 - (1.0d0-XMU1**2)/RQ2(II))
            XL        = abs(RZ2(II+1)*XMU1-RZ2(II)*XMU2)
            AMF2(II,J) = 2.d0*XL/(RZ2(II+1)-RZ2(II))
            XMU1      = XMU2
c  Lowest level intersected by emergent beam
          else
            XL        = RZ2(II+1)*XMU1*2.0d0
            AMF2(II,J) = XL/(RZ2(II+1)-RZ2(II))
            goto 16
          endif
         enddo

   16 continue

      END SUBROUTINE SPHERE2
!
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: extral
!
! !DESCRIPTION: Subroutine EXTRAL adds sub-layers to thick cloud/aerosol layers
!  using log-spacing for sub-layers of increasing thickness ATAU.
!
!     DTAUX(L=1:L1X) = Optical Depth in layer L (generally 600 nm OD)
!        This can be just cloud or cloud+aerosol, it is used only to set
!        the number in levels to insert in each layer L
!        Set for log-spacing of tau levels, increasing top-down.
!
!     N.B. the TTAU, etc calculated here are NOT used elsewhere
!
!   The log-spacing parameters have been tested for convergence and chosen
!     to be within 0.5% for ranges OD=1-500, rflect=0-100%, mu0=0.1-1.0
!     use of ATAU = 1.18 and min = 0.01, gives at most +135 pts for OD=100
!     ATAU = 1.12 now recommended for more -accurate heating rates (not J's)
!   --------------------------------------------------------------------
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE EXTRAL (DTAUX,L1X,L2X,NX,JTAUMX,ATAU,ATAU0,JXTRA)
!
! !USES:
!
!
! !INPUT PARAMETERS:
!
      INTEGER, INTENT(IN) ::  JTAUMX,L1X,L2X  !index of cloud/aerosol
      integer, intent(in) ::  NX              !Mie scattering array size
      real*8,  intent(in) ::  DTAUX(L1X)      !cloud+3aerosol OD in each layer
      real*8,  intent(in) ::  ATAU,ATAU0
!
! !OUTPUT VARIABLES:
!
      integer, intent(out)::  JXTRA(L2X+1)    !number of sub-layers to be added
!
! !REMARKS:
! 
! !REVISION HISTORY: 
!  28 Mar 2013 - S. D. Eastham - Copied from Fast-JX v7.0
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER JTOTL,I,L,L2
      REAL*8  TTAU(L2X+1),DTAUJ, ATAU1,ATAULN,ATAUM,ATAUN1

      !=================================================================
      ! EXTRAL begins here!
      !=================================================================

C---Reinitialize arrays
      TTAU(:)  = 0.d0
      JXTRA(:) = 0
c
c---combine these edge- and mid-layer points into grid of size:
c---              L2X+1 = 2*L1X+1 = 2*L_+3
c---calculate column optical depths above each level, TTAU(1:L2X+1)
c---      note that TTAU(L2X+1)=0 and TTAU(1)=total OD
c
c---Divide thick layers to achieve better accuracy in the scattering code
c---In the original fast-J, equal sub-layers were chosen, this is wasteful
c---and this new code (ver 5.3) uses log-scale:
c---        Each succesive layer (down) increase thickness by ATAU > 1
c---        e.g., if ATAU = 2, a layer with OD = 15 could be divided into
c---        4 sub-layers with ODs = 1 - 2 - 4 - 8
c---The key parameters are:
c---        ATAU = factor increase from one layer to the next
c---        ATAUMN = the smallest OD layer desired
c---        JTAUMX = maximum number of divisions (i.e., may not get to ATAUMN)
c---These are hardwired below, can be changed, but have been tested/optimized

      ATAU1  = ATAU - 1.d0
      ATAULN = log(ATAU)
        TTAU(L2X+1)  = 0.0d0
      do L2 = L2X,1,-1
        L         = (L2+1)/2
        DTAUJ     = 0.5d0 * DTAUX(L)
        TTAU(L2)  = TTAU(L2+1) + DTAUJ
c---Now compute the number of log-spaced sub-layers to be added in
c---   the interval TTAU(L2) > TTAU(L2+1)
c---The objective is to have successive TAU-layers increasing by factor ATAU >1
c---the number of sub-layers + 1
        if (TTAU(L2) .lt. ATAU0) then
          JXTRA(L2) = 0
        else
          ATAUM    = max(ATAU0, TTAU(L2+1))
          ATAUN1 = log(TTAU(L2)/ATAUM) / ATAULN
          JXTRA(L2) = min(JTAUMX, max(0, int(ATAUN1 - 0.5d0)))
        endif
      enddo

c---check on overflow of arrays, cut off JXTRA at lower L if too many levels
      JTOTL    = L2X + 2
      do L2 = L2X,1,-1
        JTOTL  = JTOTL + JXTRA(L2)
        if (JTOTL .gt. NX/2)  then
          write(6,'(A,2I5,F9.2)') 'N_/L2_/L2-cutoff JXTRA:',NX,L2X,L2
          do L = L2,1,-1
            JXTRA(L) = 0
          enddo
          go to 10
        endif
      enddo
  10  continue

      END SUBROUTINE EXTRA
!
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: exitc
!
! !DESCRIPTION: Subroutine EXITC simply writes T\_EXIT to stdout.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE EXITC (T_EXIT)
!
! !USES:
!
!
! !INPUT PARAMETERS:
!
      CHARACTER(LEN=*), INTENT(IN) ::  T_EXIT
!
! !OUTPUT VARIABLES:
!
!
! !REMARKS:
! 
! !REVISION HISTORY: 
!  28 Mar 2013 - S. D. Eastham - Copied from Fast-JX v7.0
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      !=================================================================
      ! EXITC begins here!
      !=================================================================

      WRITE(6, '(a)') T_EXIT

      END SUBROUTINE EXITC
!
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: int_fjx
!
! !DESCRIPTION: Subroutine INIT\_FJX initializes Fast-JX variables.
!\\
!\\
! !INTERFACE:
!
!      SUBROUTINE INIT_FJX(TITLEJXX,NJXU,NJXX,am_I_Root)
      SUBROUTINE INIT_FJX(NCHEM,NJXU,am_I_Root)
!
! !USES:
!
      USE inquireMod,   ONLY : findFreeLUN
!
! !INPUT PARAMETERS:
!
      INTEGER, INTENT(IN) :: NCHEM       ! Number of chemistry layers
      INTEGER, INTENT(IN) :: NJXU
      LOGICAL, INTENT(IN) :: am_I_Root   ! Is this the root CPU?
!
! !OUTPUT VARIABLES:
!
!      INTEGER, INTENT(OUT) ::NJXX
!      CHARACTER*6, INTENT(OUT), DIMENSION(NJXU) :: TITLEJXX
!
! !REMARKS:
! 
! !REVISION HISTORY: 
!  28 Mar 2013 - S. D. Eastham - Copied from Fast-JX v7.0
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER                     :: JXUNIT,J
      CHARACTER(LEN=255)          :: FILENAME
      INTEGER                     :: NJXX     ! Actual no. of rxns FJX will run
      CHARACTER*6,DIMENSION(NJXU) :: TITLEJXX

      !=================================================================
      ! INIT_FJX begins here!
      !=================================================================

      IF ( am_I_Root) THEN
         write(6,*) ' Initializing Fast-JX v7.0 standalone CTM code.'

         if (W_.ne.8 .and. W_.ne.12 .and. W_.ne.18) then
           call EXITC(' INIT_FJX: invalid no. wavelengths')
         endif
      ENDIF

      ! Set shared variables according to input arguments
      JVLACTIVE = NCHEM

      ! Get a free LUN
      JXUNIT = findFreeLUN()

      ! Read in fast-J X-sections (spectral data)
      FILENAME = 'FJX_spec.dat'

#if defined( DEVEL ) || defined( EXTERNAL_GRID ) || defined( EXTERNAL_FORCING )
      ! For the ESMF interface to GEOS-5 GCM, we need to append *___.rc
      FILENAME = TRIM( FILENAME ) // '___.rc'
#endif

      ! Echo info
      IF ( am_I_Root ) THEN
         WRITE( 6, 100 ) TRIM( FILENAME ), LUN
 100     FORMAT( 'GEOSCHEM::INIT_FJX: Reading ', a ', on unit ', i4 )
      ENDIF

      call RD_XXX(JXUNIT,TRIM( FILENAME ))

      ! Read in cloud scattering data
      FILENAME = 'FJX_scat-cld.dat'

#if defined( DEVEL ) || defined( EXTERNAL_GRID ) || defined( EXTERNAL_FORCING )
      ! For the ESMF interface to GEOS-5 GCM, we need to append *___.rc
      FILENAME = TRIM( FILENAME ) // '___.rc'
#endif

      ! Echo info
      IF ( am_I_Root ) THEN
         WRITE( 6, 100 ) TRIM( FILENAME ), LUN
      ENDIF

      call RD_CLD(JXUNIT,TRIM( FILENAME ))

      ! Read in aerosols scattering data
      FILENAME = 'FJX_scat-aer.dat'

#if defined( DEVEL ) || defined( EXTERNAL_GRID ) || defined( EXTERNAL_FORCING )
      ! For the ESMF interface to GEOS-5 GCM, we need to append *___.rc
      FILENAME = TRIM( FILENAME ) // '___.rc'
#endif

      ! Echo info
      IF ( am_I_Root ) THEN
         WRITE( 6, 100 ) TRIM( FILENAME ), LUN
      ENDIF

      call RD_MIE(JXUNIT,TRIM( FILENAME ))

      ! Read in UMich aerosol scattering data
      FILENAME = 'FJX_scat-UMa.dat'

#if defined( DEVEL ) || defined( EXTERNAL_GRID ) || defined( EXTERNAL_FORCING )
      ! For the ESMF interface to GEOS-5 GCM, we need to append *___.rc
      FILENAME = TRIM( FILENAME ) // '___.rc'
#endif

      ! Echo info
      IF ( am_I_Root ) THEN
         WRITE( 6, 100 ) TRIM( FILENAME ), LUN
      ENDIF

      call RD_UM (JXUNIT,TRIM( FILENAME ))

      ! Read in T & O3 climatology used to fill e.g. upper layers or if O3 not calc.
      FILENAME = 'atmos_std.dat'

#if defined( DEVEL ) || defined( EXTERNAL_GRID ) || defined( EXTERNAL_FORCING )
      ! For the ESMF interface to GEOS-5 GCM, we need to append *___.rc
      FILENAME = TRIM( FILENAME ) // '___.rc'
#endif

      ! Echo info
      IF ( am_I_Root ) THEN
         WRITE( 6, 100 ) TRIM( FILENAME ), LUN
      ENDIF

      call RD_PROF(JXUNIT,TRIM( FILENAME ))

        NJXX = NJX
      do J = 1,NJX
        TITLEJXX(J) = TITLEJX(J)
      enddo

      ! Read in photolysis rates used in chemistry code and mapping onto FJX J's
      ! CTM call:  read in J-values names and link to fast-JX names
      FILENAME = 'FJX_j2j.dat'

#if defined( DEVEL ) || defined( EXTERNAL_GRID ) || defined( EXTERNAL_FORCING )
      ! For the ESMF interface to GEOS-5 GCM, we need to append *___.rc
      FILENAME = TRIM( FILENAME ) // '___.rc'
#endif

      ! Echo info
      IF ( am_I_Root ) THEN
         WRITE( 6, 100 ) TRIM( FILENAME ), LUN
      ENDIF

      call RD_JS_JX(JXUNIT,TRIM( FILENAME ), TITLEJXX,NJXX)

      END SUBROUTINE INIT_FJX
!
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: fjxfunc
!
! !DESCRIPTION: Function FJXFUNC supplies J-values to the chemistry solver.  
!\\
!\\
! !INTERFACE:
!
      REAL*8 FUNCTION FJXFUNC(I,J,L,NREAC,BRCH,NAME)
!
! !USES:
!
      USE ERROR_MOD, ONLY : ERROR_STOP
!
! !INPUT PARAMETERS: 
!
      INTEGER,           INTENT(IN) :: I,J,L,NREAC,BRCH
      CHARACTER (LEN=4), INTENT(IN) :: NAME
!
! !REVISION HISTORY: 
!  (1  ) "cmn_fj.h" also includes "CMN_SIZE" and "define.h".
!  (2  ) J-values are stored in array "ZPJ" from "cmn_fjx_mod.F".
!  (3  ) Now references ERROR_STOP from "error_mod.F".  Updated comments,
!         and made some cosmetic changes. (bmy, 10/15/02)
!  27 Mar 2013 - S. D. Eastham - Rolled into FAST_JX_MOD
!  29 Mar 2013 - S. D. Eastham - Converted fro fjfunc into fjxfunc
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER                       :: N
      
      !=================================================================
      ! FJXFUNC begins here!
      !
      ! If your compiler has subscript-range checking (-C or 
      ! -check_bounds) then it is recommended to use this option to 
      ! test for the validity of (I,J,L), since repeated IF statements 
      ! are computationally expensive.
      !
      ! If your compiler does not have subscript-range checking, then 
      ! uncomment the following lines to do a manual test for the 
      ! validity of (I,J,L).  
      !=================================================================

!      IF ( I > IIPAR .OR. J > JJPAR .OR. L > JPNL ) THEN 
!         CALL ERROR_STOP('Invalid grid-box # in call to FJXFUNC',
!     & 'FAST_JX_MOD.F')
!      ENDIF

      !=================================================================
      ! RINDEX converts the J-value index as read from "globchem.dat" to 
      ! the J-value index as read from "fjx_j2j.dat". 
      ! (bmy, 10/5/98, SDE 03/29/13)
      !
      ! Make sure that we have taken the proper reaction! 
      !=================================================================
      N = RINDEX(NREAC)

      IF ( N > JVN_ ) THEN
         WRITE(6,*) 'RXN for ',name,', branch ',brch,' not found!'
         CALL ERROR_STOP( 'Check FJXFUNC', 'FAST_JX_MOD.F' )
      ENDIF

      !=================================================================
      ! Return the appropriate J-value as the value of the function 
      !=================================================================
      FJXFUNC = ZPJ(L,N,I,J)

      ! Return to calling program
      END FUNCTION FJXFUNC
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: rd_xxx
!
! !DESCRIPTION: Subroutine RD\_XXX reads in wavelength bins, solar fluxes,
!  Rayleigh and temperature-dependent cross-sections.
!    NEW v-6.8  now allow 1 to 3 sets of X-sects for T or P
!           LQQ = 1, 2, or 3 to determine interpolation with T or P
!           IF the temperatures TQQQ are <0, then use as pressure interp (hPa)
!           NB - the temperatures and pressures must be increasing
!    NEW v-6.4  changed to collapse wavelengths & x-sections to Trop-only:
!           WX_ = 18 should match the JX_spec.dat wavelengths
!           W_ = 12 (Trop-only) or 18 (std) is set in (cmn_FXJ.f).
!       if W_=12 then drop strat wavels, and drop x-sects (e.g. N2O, ...)
!           W_ = 8, reverts to quick fix:  fast-J (12-18) plus bin (5) scaled
!
!   --------------------------------------------------------------------
!     NAMFIL   Name of spectral data file (JX_spec.dat) >> j2 for fast-J2
!     NUN      Channel number for reading data file
!
!     NJX    Number of species to calculate J-values for
!     NWWW     Number of wavelength bins, from 1:NWWW
!     WBIN     Boundaries of wavelength bins
!     WL       Centres of wavelength bins - 'effective wavelength'
!     FL       Solar flux incident on top of atmosphere (cm-2.s-1)
!     QRAYL    Rayleigh parameters (effective cross-section) (cm2)
!     QO2      O2 cross-sections
!     QO3      O3 cross-sections
!     Q1D      O3 => O(1D) quantum yield
!     TQQ      Temperature for supplied cross sections
!     QQQ      Supplied cross sections in each wavelength bin (cm2)
!   --------------------------------------------------------------------
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE RD_XXX (NUN,NAMFIL)
!
! !USES:
!
!
! !INPUT PARAMETERS:
!
      INTEGER, INTENT(IN) :: NUN
      CHARACTER(*), INTENT(IN) ::  NAMFIL
!
! !OUTPUT VARIABLES:
!
!
! !REMARKS:
! 
! !REVISION HISTORY: 
!  28 Mar 2013 - S. D. Eastham - Copied from Fast-JX v7.0
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER  I, J, JJ, K, IW, NQRD, NWWW,   LQ
      REAL*8  QQ2(199), TQQ2

      character*78 TITLE0
      character*6  TITLEJ2,TITLEJ3
      character*1  TSTRAT

      !=================================================================
      ! RD_XXX begins here!
      !=================================================================

      TQQ(:,:) = 0.d0

C----------spectral data----set for new format data------------------
c         note that X_ = max # Xsects read in
c                   NJX = # fast-JX J-values derived from this (.le. X_)

c >>>> W_ = 12 <<<< means trop-only, discard WL #1-4 and #9-10, some X-sects

      open (NUN,FILE=NAMFIL,status='old',form='formatted')
      read (NUN,100) TITLE0

c----note that NQRD is not used any more, a read until 'endofJ' is performed
      read (NUN,101) NQRD,NWWW
         NW1 = 1
         NW2 = NWWW
      write(6,'(1x,a)') TITLE0
      write(6,'(i8)') NWWW
C----J-values:  1=O2, 2=O3P,3=O3D 4=readin Xsects
      read (NUN,102) (WL(IW),IW=1,NWWW)
      read (NUN,102) (FL(IW),IW=1,NWWW)
      read (NUN,102) (QRAYL(IW),IW=1,NWWW)

C---Read O2 X-sects, O3 X-sects, O3=>O(1D) quant yields (each at 3 temps)
C---NB the O3 and q-O3-O1D are at different temperatures and cannot be combined
      read (NUN,103) TITLEJX(1),TQQ(1,1), (QO2(IW,1),IW=1,NWWW)
      read (NUN,103) TITLEJ2,  TQQ(2,1), (QO2(IW,2),IW=1,NWWW)
      read (NUN,103) TITLEJ3,  TQQ(3,1), (QO2(IW,3),IW=1,NWWW)

      read (NUN,103) TITLEJX(2),TQQ(1,2), (QO3(IW,1),IW=1,NWWW)
      read (NUN,103) TITLEJ2,  TQQ(2,2), (QO3(IW,2),IW=1,NWWW)
      read (NUN,103) TITLEJ3,  TQQ(3,2), (QO3(IW,3),IW=1,NWWW)

      read (NUN,103) TITLEJX(3),TQQ(1,3), (Q1D(IW,1),IW=1,NWWW)
      read (NUN,103) TITLEJ2,  TQQ(2,3), (Q1D(IW,2),IW=1,NWWW)
      read (NUN,103) TITLEJ3,  TQQ(3,3), (Q1D(IW,3),IW=1,NWWW)

      SQQ(1) = ' '
      SQQ(2) = ' '
      SQQ(3) = ' '

      LQQ(1) = 3
      LQQ(2) = 3
      LQQ(3) = 3

C---Read remaining species:  X-sections at 1-2-3 T_s
        JJ = 3
      do I=1,9999

c--try to read in 3 X-sects per J-value (JJ)
        read (NUN,104) TITLEJ2,TSTRAT,TQQ2,(QQ2(IW),IW=1,NWWW)
        if (TITLEJ2 .eq. 'endofJ') goto 1
c---skip stratosphere only J's (denoted by 'x')if W_<18 => trop-only J's
        if (W_.eq.18 .or. TSTRAT.ne.'x') then
          if (TITLEJ2 .ne. TITLEJX(JJ)) then
               JJ = JJ+1

         if (JJ .gt. X_) then
         call EXITC(' RD_XXX:  X_ not large enough for Xsects read in')
         endif

               TITLEJX(JJ) = TITLEJ2
               LQQ(JJ) = 1
               SQQ(JJ) = TSTRAT
                 LQ = LQQ(JJ)
               TQQ(LQ,JJ) = TQQ2
             do IW = 1,NWWW
               QQQ(IW,LQ,JJ) = QQ2(IW)
             enddo
           else
               LQQ(JJ) = LQQ(JJ)+1
            if (LQQ(JJ) .le. 3) then
               LQ = LQQ(JJ)
               TQQ(LQ,JJ) = TQQ2
             do IW = 1,NWWW
               QQQ(IW,LQ,JJ) = QQ2(IW)
             enddo
            endif
           endif
        endif
      enddo
    1 continue
        NJX = JJ

      do J = 1,NJX
        write(6,200) J,TITLEJX(J),SQQ(J),LQQ(J),(TQQ(I,J),I=1,LQQ(J))
c---need to check that TQQ is monotonically increasing:
        if (LQQ(J) .eq. 3) then
          if (TQQ(2,J) .ge. TQQ(3,J)) then
            call EXITC ('TQQ out of order')
          endif
          if (TQQ(1,J) .ge. TQQ(2,J)) then
            call EXITC ('TQQ out of order')
          endif
        endif
        if (LQQ(J) .eq. 2) then
          if (TQQ(1,J) .ge. TQQ(2,J)) then
            call EXITC ('TQQ out of order')
          endif
        endif
      enddo

c---check on doingpressure interp
c---check on consolidating Qo2 and others into
c---wrte a newFJX_J2J.dat for mapping on fjx Xsects


c---truncate number of wavelengths to do troposphere-only
      if (W_ .ne. WX_) then
c---TROP-ONLY
       if (W_ .eq. 12) then
        write(6,'(a)')
     &   ' >>>TROP-ONLY reduce wavelengths to 12, drop strat X-sects'
        NW2 = 12
        do IW = 1,4
          WL(IW) = WL(IW+4)
          FL(IW) = FL(IW+4)
          QRAYL(IW) = QRAYL(IW+4)
         do K = 1,3
          QO2(IW,K) = QO2(IW+4,K)
          QO3(IW,K) = QO3(IW+4,K)
          Q1D(IW,K) = Q1D(IW+4,K)
         enddo
         do J = 4,NJX
          QQQ(IW,1,J) = QQQ(IW+4,1,J)
          QQQ(IW,2,J) = QQQ(IW+4,2,J)
         enddo
        enddo
        do IW = 5,12
          WL(IW) = WL(IW+6)
          FL(IW) = FL(IW+6)
          QRAYL(IW) = QRAYL(IW+6)
         do K = 1,3
          QO2(IW,K) = QO2(IW+6,K)
          QO3(IW,K) = QO3(IW+6,K)
          Q1D(IW,K) = Q1D(IW+6,K)
         enddo
         do J = 4,NJX
          QQQ(IW,1,J) = QQQ(IW+6,1,J)
          QQQ(IW,2,J) = QQQ(IW+6,2,J)
         enddo
        enddo
c---TROP-QUICK  (must scale solar flux for W=5)
       elseif (W_ .eq. 8) then
        write(6,'(a)')
     &   ' >>>TROP-QUICK reduce wavelengths to 8, drop strat X-sects'
        NW2 = 8
        do IW = 1,1
          WL(IW) = WL(IW+4)
          FL(IW) = FL(IW+4)  * 2.d0
          QRAYL(IW) = QRAYL(IW+4)
         do K = 1,3
          QO2(IW,K) = QO2(IW+4,K)
          QO3(IW,K) = QO3(IW+4,K)
          Q1D(IW,K) = Q1D(IW+4,K)
         enddo
         do J = 4,NJX
          QQQ(IW,1,J) = QQQ(IW+4,1,J)
          QQQ(IW,2,J) = QQQ(IW+4,2,J)
         enddo
        enddo
        do IW = 2,8
          WL(IW) = WL(IW+10)
          FL(IW) = FL(IW+10)
          QRAYL(IW) = QRAYL(IW+10)
         do K = 1,3
          QO2(IW,K) = QO2(IW+10,K)
          QO3(IW,K) = QO3(IW+10,K)
          Q1D(IW,K) = Q1D(IW+10,K)
         enddo
         do J = 4,NJX
          QQQ(IW,1,J) = QQQ(IW+10,1,J)
          QQQ(IW,2,J) = QQQ(IW+10,2,J)
         enddo
        enddo

       else
         call EXITC(' no. wavelengths wrong: W_ .ne. 8,12,18')
       endif
      endif

      close(NUN)

  100 format(a)
  101 format(10x,5i5)
  102 format(10x,    6e10.3/(10x,6e10.3)/(10x,6e10.3))
  103 format(a6,1x,f3.0,6e10.3/(10x,6e10.3)/(10x,6e10.3))
  104 format(a6,a1,f3.0,6e10.3/(10x,6e10.3)/(10x,6e10.3))
  200 format(1x,' x-sect:',i3,a10,a4,i5,3(3x,f6.2))
  201 format(' Number of x-sections supplied to Fast-J2: ',i3,/,
     &       ' Maximum number allowed (X_) only set to: ',i3,
     &       ' - increase in cmn_FJX.f')

      END SUBROUTINE RD_XXX
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: rd_cld
!
! !DESCRIPTION: Subroutine RD\_CLD reads in aerosol/cloud scattering data.
!   --------------------------------------------------------------------
!     NAMFIL   Name of scattering data file (e.g., FJX_scat.dat)
!     NUN      Channel number for reading data file
!     NCC      Number of categories for cloud scattering phase functions
!     QCC      Cloud scattering phase functions
!     WCC      5 Wavelengths for supplied phase functions
!     PCC      Phase function: first 8 terms of expansion
!     RCC      Effective radius associated with cloud type
!     SCC      Single scattering albedo
!     DCC      density (g/cm^3)
!   --------------------------------------------------------------------
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE RD_CLD (NUN,NAMFIL)
!
! !USES:
!
!
! !INPUT PARAMETERS:
!
      INTEGER, INTENT(IN) :: NUN
      CHARACTER(*), INTENT(IN) ::  NAMFIL
!
! !OUTPUT VARIABLES:
!
!
! !REMARKS:
! 
! !REVISION HISTORY: 
!  28 Mar 2013 - S. D. Eastham - Copied from Fast-JX v7.0
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER  I, J, K
      CHARACTER*78 TITLE0
      CHARACTER*20 TITLAA(A_)   ! TITLAA: Title for scatering data

      !=================================================================
      ! RD_CLD begins here!
      !=================================================================

      open (NUN,FILE=NAMFIL,status='old',form='formatted')

      read (NUN,'(i2,a78)') NCC,TITLE0
        if (NCC .gt. C_) then
          call EXITC(' too many cld-data sets: NCC > C_')
        endif

      write(6,'(a,2f9.5,i5)') ' ATAU/ATAU0/JMX',ATAU,ATAU0,JTAUMX

      read (NUN,*)
      read (NUN,*)
      do J = 1,NCC
          read (NUN,'(3x,a8,1x,2f6.3)') TITLAA(J),RCC(J),DCC(J)
        do K = 1,5
          read (NUN,'(f4.0,f7.4,f7.4,7f6.3)')
     &  WCC(K,J),QCC(K,J),SCC(K,J),(PCC(I,K,J),I=2,8)
          PCC(1,K,J) = 1.d0
        enddo
      enddo

      close(NUN)

      write(6,'(a,9f8.1)') ' Aerosol optical: r-eff/rho/Q(@wavel):'
     &             ,(WCC(K,1),K=1,5)
      write(6,*) TITLE0
      do J=1,NCC
      write(6,'(i3,1x,a8,7f8.3)')
     &              J,TITLAA(J),RCC(J),DCC(J),(QCC(K,J),K=1,5)
      enddo

      END SUBROUTINE RD_CLD
!
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: rd_mie
!
! !DESCRIPTION: Subroutine RD\_MIE retrieves aerosol scattering data for FJX.
!   --------------------------------------------------------------------
!     NAMFIL   Name of scattering data file (e.g., FJX_scat.dat)
!     NUN      Channel number for reading data file
!     NAA      Number of categories for scattering phase functions
!     QAA      Aerosol scattering phase functions
!     WAA      5 Wavelengths for the supplied phase functions
!     PAA      Phase function: first 8 terms of expansion
!     RAA      Effective radius associated with aerosol type
!     SAA      Single scattering albedo
!     DAA      density (g/cm^3)
!   --------------------------------------------------------------------
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE RD_MIE(NUN,NAMFIL)
!
! !USES:
!
!
! !INPUT PARAMETERS:
!
      INTEGER, INTENT(IN) :: NUN
      CHARACTER(*), INTENT(IN) ::  NAMFIL
!
! !OUTPUT VARIABLES:
!
!
! !REMARKS:
! 
! !REVISION HISTORY: 
!  28 Mar 2013 - S. D. Eastham - Copied from Fast-JX v7.0
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER  I, J, K
      CHARACTER*78 TITLE0
      CHARACTER*20 TITLAA(A_)   ! TITLAA: Title for scatering data

      !=================================================================
      ! RD_MIE begins here!
      !=================================================================

      open (NUN,FILE=NAMFIL,status='old',form='formatted')

      read (NUN,'(i2,a78)') NAA,TITLE0
        if (NAA .gt. A_) then
          call EXITC(' too many aerosol-data sets: NAA > A_')
        endif

      write(6,'(a,2f9.5,i5)') ' ATAU/ATAU0/JMX',ATAU,ATAU0,JTAUMX

      read (NUN,*)
      read (NUN,*)
      do J = 1,NAA
          read (NUN,'(3x,a8,1x,2f6.3)') TITLAA(J),RAA(J),DAA(J)
        do K = 1,5
          read (NUN,'(f4.0,f7.4,f7.4,7f6.3)')
     &  WAA(K,J),QAA(K,J),SAA(K,J),(PAA(I,K,J),I=2,8)
          PAA(1,K,J) = 1.d0
        enddo
      enddo

      close(NUN)

      write(6,'(a,9f8.1)') ' Aerosol optical: r-eff/rho/Q(@wavel):'
     &             ,(WAA(K,1),K=1,5)
      write(6,*) TITLE0
      do J=1,NAA
      write(6,'(i3,1x,a8,7f8.3)')
     &              J,TITLAA(J),RAA(J),DAA(J),(QAA(K,J),K=1,5)
      enddo

      END SUBROUTINE RD_MIE
!
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: rd_um
!
! !DESCRIPTION: Subroutine RD\_UM retrieves U. Michigan optical data.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE RD_UM(NUN,NAMFIL)
!
! !USES:
!
!
! !INPUT PARAMETERS:
!
      INTEGER, INTENT(IN) :: NUN
      CHARACTER(*), INTENT(IN) ::  NAMFIL
!
! !OUTPUT VARIABLES:
!
!
! !REMARKS:
! 
! !REVISION HISTORY: 
!  28 Mar 2013 - S. D. Eastham - Copied from Fast-JX v7.0
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER  I, J, K, L
      CHARACTER*78 TITLE0
      CHARACTER*20 TITLUM(33)   ! TITLUM: Title for U Michigan aerosol data set

      !=================================================================
      ! RD_UM begins here!
      !=================================================================

      open (NUN,FILE=NAMFIL,status='old',form='formatted')

      read (NUN,'(a78)') TITLE0
        write(6,*) 'UMichigan Aerosols', TITLE0
      read(NUN,'(5x,10f5.0)') WMM
        write(6,'(a,10f7.1)') ' UMIchigan aerosol wavelengths:',WMM

c---33 Different UM Aerosol Types:  SULF, SS-1,-2,-3,-4, DD-1,-2,-3,-4,
c---      FF00(0%BC), FF02, ...FF14(14%BC),  BB00, BB02, ...BB30(30%BC)
      do L=1,33
          read(NUN,'(a4)') TITLUM(L)
c---21 Rel Hum:    K=1=0%, =2=5%, ... =20=95%, =21=99%
        do K=1,21
c---6 wavelengths: J=1=200nm, 2=300nm, 3=400nm, (4'=550nm) 5=600nm, 6=1000nm
c---3 optic vars:  I=1=SSAlbedo,  =2=g,  =3=k-ext
          read(NUN,'(18f9.5)')  ((UMAER(I,J,K,L),I=1,3),J=1,6)
        enddo
      enddo

      close(NUN)

        write(6,'(a)') 'collapse UM wavelengths, drop 550 nm'
          WMM(4) = WMM(5)
          WMM(5) = WMM(6)
       do L=1,33
       do K=1,21
       do I=1,3
          UMAER(I,4,K,L) = UMAER(I,5,K,L)
          UMAER(I,5,K,L) = UMAER(I,6,K,L)
       enddo
       enddo
       enddo

        write(6,'(7(i5,1x,a4))') (L,TITLUM(L), L=1,33)

      END SUBROUTINE RD_UM

!
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: rd_prof
!
! !DESCRIPTION: Subroutine RD\_PROF inputs temperature and O3 reference 
!  profiles from 'atmos_std.dat'.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE RD_PROF(NJ2,NAMFIL)
!
! !USES:
!
!
! !INPUT PARAMETERS:
!
      INTEGER, INTENT(IN) ::  NJ2
      CHARACTER(*), INTENT(IN) ::  NAMFIL
!
! !OUTPUT VARIABLES:
!
!
! !REMARKS:
! 
! !REVISION HISTORY: 
!  28 Mar 2013 - S. D. Eastham - Copied from Fast-JX v7.0
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      integer IA, I, M, L, LAT, MON, NTLATS, NTMONS, N216
      REAL*8  OFAC, OFAK
      CHARACTER*78 TITLE0

      !=================================================================
      ! RD_PROF begins here!
      !=================================================================

      open (NJ2,file=NAMFIL,status='old',form='formatted')
      read (NJ2,'(A)') TITLE0
      read (NJ2,'(2I5)') NTLATS,NTMONS
c      write(6,'(1X,A)') TITLE0
      write(6,1000) NTLATS,NTMONS
      N216  = min(216, NTLATS*NTMONS)
      do IA = 1,N216
        read (NJ2,'(1X,I3,3X,I2)') LAT, MON
        M = min(12, max(1, MON))
        L = min(18, max(1, (LAT+95)/10))
        read (NJ2,'(3X,11F7.1)') (TREF(I,L,M), I=1,41)
        read (NJ2,'(3X,11F7.4)') (OREF(I,L,M), I=1,31)
      enddo
      close (NJ2)

c  Extend climatology to 100 km
      OFAC = exp(-2.d5/5.d5)
      do I = 32,51
        OFAK = OFAC**(I-31)
        do M = 1,NTMONS
        do L = 1,NTLATS
          OREF(I,L,M) = OREF(31,L,M)*OFAK
        enddo
        enddo
      enddo
      do L = 1,NTLATS
      do M = 1,NTMONS
      do I = 42,51
        TREF(I,L,M) = TREF(41,L,M)
      enddo
      enddo
      enddo

 1000 format(1x,'std atmos profiles: ',i3,' lat x ',i2,' mon')

      END SUBROUTINE RD_PROF
!
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: rd_js_jx
!
! !DESCRIPTION: Subroutine RD\_JS\_JX reads in 'FJX_j2j.dat', which defines 
!  the mapping of Fast-JX J's (TITLEJX(1:NJX)) onto the CTM reactions.
!  Reaction number JJ, named T_REACT, uses Fast-JX's T_FJX (including scaling
!  factor F_FJX).
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE RD_JS_JX(NUNIT,NAMFIL,TITLEJX,NJX)
!
! !USES:
!
!
! !INPUT PARAMETERS:
!
      INTEGER, INTENT(IN)                    ::  NUNIT, NJX
      CHARACTER(*), INTENT(IN)               ::  NAMFIL
      CHARACTER*6, INTENT(IN),DIMENSION(NJX) :: TITLEJX
!
! !OUTPUT VARIABLES:
!
!
! !REMARKS:
! 
! !REVISION HISTORY: 
!  28 Mar 2013 - S. D. Eastham - Copied from Fast-JX v7.0
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER   J,JJ,K
      CHARACTER*120 CLINE
      CHARACTER*50 T_REACT
      CHARACTER*6  JMAP(JVN_), T_FJX
      REAL*8 F_FJX

      !=================================================================
      ! RD_JS_JX begins here!
      !=================================================================

c Read the FJX_j2j.dat file to map model specific J's onto fast-JX J's
c The chemistry code title describes fully the reaction (a50)
c Blank (unfilled) chemistry J's are unmapped
c The number NRATJ is the last JJ readin that is .le. JVN
c   include fractional quantum yield for the fast-JX J's

      JLABEL(:) = '------'
      JMAP(:)   = '------'
      JFACTA(:) = 0.d0

      open (NUNIT,file=NAMFIL,status='old',form='formatted')

       read (NUNIT,'(a)') CLINE
         write(6,'(a)') CLINE
      do J = 1,JVN_
       read (NUNIT,'(i4,1x,a50,4x,f5.3,2x,a6)') JJ,T_REACT,F_FJX,T_FJX
       if (JJ.gt.JVN_) goto 20
        JLABEL(JJ) = T_REACT
        JFACTA(JJ) = F_FJX
        JMAP(JJ) = T_FJX
        NRATJ = JJ
      enddo

 20   close(NUNIT)

c---Zero / Set index arrays that map Jvalue(j) onto rates
      do K = 1,NRATJ
         JIND(K) = 0
       do J = 1,NJX
         T_FJX = TITLEJX(J)
        if (JMAP(K) .eq. TITLEJX(J)) then
         JIND(K)=J
        endif
       enddo
      enddo

      write(6,'(a,i4,a)')' Photochemistry Scheme with',NRATJ,' J-values'
      do K=1,NRATJ
       if (JMAP(K) .ne. '------' ) then
        J = JIND(K)
        if (J.eq.0) then
         write(6,'(i5,a50,f6.3,a,1x,a6)') K,JLABEL(K),JFACTA(K),
     &         ' no mapping onto fast-JX',JMAP(K)
        else
         write(6,'(i5,a50,f6.3,a,i4,1x,a6)') K,JLABEL(K),JFACTA(K),
     &         ' mapped to FJX:',J,TITLEJX(J)
        endif
       endif
      enddo

      END SUBROUTINE RD_JS_JX
!
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: solar_jx
!
! !DESCRIPTION: Subroutine SOLAR\_JX handles solar zenith angles.
!  ---------------------------------------------------------------------
!     GMTIME = UT for when J-values are wanted
!           (for implicit solver this is at the end of the time step)
!     NDAY   = integer day of the year (used for solar lat and declin)
!     YGRDJ  = laitude (radians) for grid (I,J)
!     XGDRI  = longitude (radians) for grid (I,J)
!
!     SZA = solar zenith angle in degrees
!     COSSZA = U0 = cos(SZA)
!  ---------------------------------------------------------------------
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE SOLAR_JX(NDAY,COSSZA,SZA,SOLFX)
!
! !USES:
!
!
! !INPUT PARAMETERS:
!
      REAL*8, INTENT(IN) ::   COSSZA
      INTEGER, INTENT(IN) ::  NDAY
!
! !OUTPUT VARIABLES:
!
      REAL*8, INTENT(OUT) ::  SZA,SOLFX
!
! !REMARKS:
! 
! !REVISION HISTORY: 
!  28 Mar 2013 - S. D. Eastham - Adapted from Fast-JX v7.0
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      REAL*8  PI, PI180!, LOCT
!      REAL*8  SINDEC, SOLDEK, COSDEC, SINLAT, SOLLAT, COSLAT, COSZ

      !=================================================================
      ! SOLAR_JX begins here!
      !=================================================================

      PI     = 3.141592653589793d0
      PI180  = PI/180.d0
!      SINDEC = 0.3978d0*sin(0.9863d0*(dble(NDAY)-80.d0)*PI180)
!      SOLDEK = asin(SINDEC)
!      COSDEC = cos(SOLDEK)
!      SINLAT = sin(YGRDJ)
!      SOLLAT = asin(SINLAT)
!      COSLAT = cos(SOLLAT)
c
!      LOCT   = (((GMTIME)*15.d0)-180.d0)*PI180 + XGRDI
!      COSSZA = COSDEC*COSLAT*cos(LOCT) + SINDEC*SINLAT
      SZA    = acos(COSSZA)/PI180
c
c      write(6,*) ' XGRDI,YGRDJ',XGRDI,YGRDJ
c      write(6,*) ' LOCT (rad)',LOCT
c      write(6,*) ' SINDEC,COSDEC', SINDEC,COSDEC
c      write(6,*) ' SINLAT,COSLAT', SINLAT,COSLAT
c      write(6,*) ' COS, SZA',COSSZA,SZA
c
      ! Offset used for GEOS-Chem slightly different
      !SOLFX  = 1.d0-(0.034d0*cos(dble(NDAY-186)*2.d0*PI/365.d0))
      SOLFX  = 1.d0-(0.034d0*cos(dble(NDAY-172)*2.d0*PI/365.d0))

      END SUBROUTINE SOLAR_JX
!
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: aclim_fjx
!
! !DESCRIPTION: Subroutine ACLIM\_FJX loads Fast-JX climatology for a given 
!  latitude and month (for a specific pressure grid)
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE ACLIM_FJX (YLATD, MONTH, PPP,TTT,ZZZ,DDD,OOO, L1U)
!
! !USES:
!
!
! !INPUT PARAMETERS:
!
      REAL*8, INTENT(IN)  :: YLATD
      INTEGER, INTENT(IN) :: MONTH, L1U
      REAL*8, INTENT(IN),  DIMENSION(L1U+1) :: PPP
!
! !OUTPUT VARIABLES:
!
      REAL*8, INTENT(OUT), DIMENSION(L1U+1) :: ZZZ
      REAL*8, INTENT(OUT), DIMENSION(L1U) :: TTT,DDD,OOO
!
! !REMARKS:
! 
! !REVISION HISTORY: 
!  28 Mar 2013 - S. D. Eastham - Copied from Fast-JX v7.0
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      REAL*8, DIMENSION(51) :: OREF2,TREF2
      REAL*8, DIMENSION(52) :: PSTD
      INTEGER  K, L, M, N
      REAL*8   DLOGP,F0,T0,PB,PC,XC,MASFAC,SCALEH

      !=================================================================
      ! ACLIM_FJX begins here!
      !=================================================================

c  Select appropriate month
      M = max(1,min(12,MONTH))
c  Select appropriate latitudinal profiles
      N = max(1, min(18, (int(YLATD+99)/10 )))
      do K = 1,51
        OREF2(K) = OREF(K,N,M)
        TREF2(K) = TREF(K,N,M)
      enddo

c  Apportion O3 and T on supplied climatology z levels onto CTM levels +1
c  with mass (pressure) weighting, assuming constant mixing ratio and
c  temperature half a layer on either side of the point supplied.
c   PPP(L=1:L1_)=edge-pressure of CTM layer, PPP(L1_+1)=0 (top-of-atmos)
c  Mass factor - delta-Pressure (mbars) to delta-Column (molecules.cm-2)
      MASFAC = 100.d0*6.022d+23/(28.97d0*9.8d0*10.d0)

c  Set up pressure levels for O3/T climatology - assume that value
c  given for each 2 km z* level applies from 1 km below to 1 km above,
c  so select pressures at these boundaries. Surface level values at
c  1000 mb are assumed to extend down to the actual PSURF (if > 1000)
           PSTD(1) = max(PPP(1),1000.d0)
           PSTD(2) = 1000.d0 * 10.d0**(-1.d0/16.d0)
           DLOGP   = 10.d0**(-2.d0/16.d0)
      do K = 3,51
        PSTD(K) = PSTD(K-1)*DLOGP
      enddo
        PSTD(52)  = 0.d0
      do L = 1,L1U
        F0 = 0.d0
        T0 = 0.d0
        do K = 1,51
          PC   = min(PPP(L),PSTD(K))
          PB   = max(PPP(L+1),PSTD(K+1))
          if (PC .gt. PB) then
            XC = (PC-PB)/(PPP(L)-PPP(L+1))
            F0 = F0 + OREF2(K)*XC
            T0 = T0 + TREF2(K)*XC
          endif
        enddo
        TTT(L)  = T0
        DDD(L)  = (PPP(L)-PPP(L+1))*MASFAC
        OOO(L) = F0*1.d-6*DDD(L)
      enddo

c  Calculate effective altitudes using scale height at each level
        ZZZ(1) = 0.d0
      do L = 1,L1U-1
        SCALEH      = 1.3806d-19*MASFAC*TTT(L)
        ZZZ(L+1) = ZZZ(L) -( LOG(PPP(L+1)/PPP(L)) * SCALEH )
      enddo
        ZZZ(L1U+1) = ZZZ(L1U) + ZZHT

      END SUBROUTINE ACLIM_FJX
!
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: photo_jx
!
! !DESCRIPTION: Subroutine PHOTO\_JX is the core subroutine of Fast-JX.
!    calc J's for a single column atmosphere (aka Indep Colm Atmos or ICA)
!    needs P, T, O3, clds, aersls; adds top-of-atmos layer from climatology
!    needs day-of-year for solar functions.
!    Aerosol data is taken from GEOS-Chem's existing arrays, rather than
!    calculating using the U. Michigan techniques (as is the case in
!    original Fast-JX)
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE PHOTO_JX(NLON,NLAT,U0,DAY_YR,DAY_M,REFLB,LPRTJ,
     & PPP,ZZZ,TTT,DDD,RRR,CLDWP,AERSP,NDXCLD,NDXAER,L1U,ANU,
     &  VALJXX,NJXU,am_I_Root)
!
! !USES:
!
!
! !INPUT PARAMETERS:
!
      INTEGER, INTENT(IN)                   :: L1U,ANU,NJXU
      INTEGER, INTENT(IN)                   :: DAY_YR,DAY_M
      INTEGER, INTENT(IN)                   :: NLON, NLAT
      REAL*8,  INTENT(IN)                   :: U0,REFLB
      LOGICAL, INTENT(IN)                   :: LPRTJ
      REAL*8,  INTENT(IN),DIMENSION(L1U+1)  :: PPP
      REAL*8,  INTENT(IN),DIMENSION(L1U  )  :: TTT,DDD,RRR,OOO,CLDWP
      REAL*8,  INTENT(IN),DIMENSION(L1U,ANU):: AERSP
      INTEGER, INTENT(IN),DIMENSION(L1U  )  :: NDXCLD
      INTEGER, INTENT(IN),DIMENSION(L1U,ANU):: NDXAER
      LOGICAL, INTENT(IN)                   :: am_I_Root
!
! !OUTPUT VARIABLES:
!
      REAL*8, INTENT(OUT), DIMENSION(L1U-1,NJXU)::  VALJXX
!
! !REMARKS:
! 
! !REVISION HISTORY: 
!  28 Mar 2013 - S. D. Eastham - Copied from Fast-JX v7.0
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      REAL*8                     :: SZA     ! Solar zenith angle (deg)
      REAL*8                     :: SOLF    ! Solar function
c-----------------------------------------------------------------------
c---key LOCAL atmospheric data needed to solve plane-parallel J----
c-----these are dimensioned JXL_, and must have JXL_ .ge. L_
      REAL*8, DIMENSION(JXL1_+1)    :: TTJ,DDJ,OOJ,ZZJ
      REAL*8, DIMENSION(JXL1_+1)    :: PPJ,RELH
      INTEGER,DIMENSION(JXL2_+1)    :: JXTRA
      REAL*8, DIMENSION(W_)         :: FJTOP,FJBOT,FSBOT,FLXD0,RFL
      REAL*8, DIMENSION(JXL_, W_)   ::  AVGF, FJFLX
      REAL*8, DIMENSION(JXL1_,W_)   ::  DTAUX, FLXD
      REAL*8, DIMENSION(8,JXL1_,W_) ::  POMEGAX
      REAL*8, DIMENSION(JXL1_)      ::  DTAU600
      REAL*8, DIMENSION(8,JXL1_)    ::  POMG600
      REAL*8, DIMENSION(W_,JXL1_)   ::  FFX
      REAL*8, DIMENSION(W_,8)       ::  FFXNET
      REAL*8, DIMENSION(L1U+1)      ::  ZZZ
c---flux/heating arrays (along with FJFLX,FLXD,FLXD0)
      real*8  FLXJ(JXL1_),FFX0,FXBOT,FABOT
      real*8             :: ODABS,ODRAY
      real*8             :: RFLECT,FREFS,FREFL,FREFI
      real*8             :: AMF2(2*JXL1_+1,2*JXL1_+1)
c------------key SCATTERING arrays for clouds+aerosols------------------
      real*8             :: OPTX(5),SSAX(5),SLEGX(8,5)
      real*8             :: OD(5,JXL1_),SSA(5,JXL1_),SLEG(8,5,JXL1_)
      real*8             :: OD600(JXL1_)
      real*8             :: PATH,DENSWP,DENS,RH,RHO,REFF,XTINCT,ODCLD
c------------key arrays AFTER solving for J's---------------------------
      real*8  FFF(W_,JXL_),VALJ(X_)
      real*8  FLXUP(W_),FLXDN(W_),DIRUP(W_),DIRDN(W_)
      real*8  VALJL(JXL_,X_) !2-D array of J_s returned by JRATET

      integer  LU,L2U, I,J,K,L,M,KMIE,KW,NAER,NDCLD, RATIO(W_)
      real*8   XQO3,XQO2,DTAUC,WAVE, TTTX
C-----------------------------------------------------------------------

      !=================================================================
      ! PHOTO_JX begins here!
      !=================================================================

      ! Input conversion (SDE 03/29/13)
      ! Calculate solar zenith angle (degrees)
      CALL SOLAR_JX(DAY_YR,U0,SZA,SOLF)

      if (L1U .gt. JXL1_) then
        call EXITC(' PHOTO_JX: not enough levels in JX')
      endif

      LU = L1U - 1
      L2U = LU + LU + 2

      FFF(:,:) = 0.d0
      FREFI = 0.d0
      FREFL = 0.d0
      FREFS = 0.d0

c---check for dark conditions SZA > 98.0 deg => tan ht = 63 km
c                        or         99.0                 80 km
      if (SZA .gt. 98.d0) goto 99

      !-----------------------------------------------------------------
      !### If you want to set SZA = 0 degrees for testing,
      !### then uncomment the following lines (bmy, 9/13/99) 
      !U0  = 1.0d0
      !SZA = 0.0d0
      !-----------------------------------------------------------------

      ! Use SET_PROF to acquire atmospheric column data
      CALL SET_PROF( NLON, NLAT, YLAT, MONTH, DAY,
     &          PPP ,TTT, REFLB, ODCOL, OPTDUST, OPTAER, 
     &          TTJ, DO3, DBC, DM, ZZZ, AER )

      DO L = 1, L1U
         ZZJ(L) = ZZZ(L)
      ENDDO
      ZZJ(L1U+1) = ZZJ(L1U) + ZZHT

        RFLECT = REFLB

C-----------------------------------------------------------------------
      call SPHERE2 (U0,RAD,ZZJ,AMF2,L1U,JXL1_)
C-----------------------------------------------------------------------

      DO L=1,LU
         OD600(L) = ODCOL(L)
      ENDDO

      OD600(LU+1) = 0.d0

C-----------------------------------------------------------------------
      call EXTRAL(OD600,L1U,L2U,N_,JTAUMX,ATAU,ATAU0, JXTRA)
C-----------------------------------------------------------------------

c---set surface reflectance
        RFL(:) = max(0.d0,min(1.d0,RFLECT))

c-----------------------------------------------------------------------
C---Loop over all wavelength bins to calc mean actinic flux AVGF(L)
c-----------------------------------------------------------------------

      do K = 1,W_

        WAVE = WL(K)
C---Pick nearest Mie wavelength to get scattering properites------------
                               KMIE=1  ! use 200 nm prop for <255 nm
        if( WAVE .gt. 255.d0 ) KMIE=2  ! use 300 nm prop for 255-355 nm
        if( WAVE .gt. 355.d0 ) KMIE=3  ! use 400 nm prop for 355-500 nm
        if( WAVE .gt. 500.d0 ) KMIE=4  ! use 600 nm prop for 500-800 nm
        ! Truncate here (KMIE=5 not appropriate in our applicaation)

C---For Mie code scale extinction at 1000 nm to wavelength WAVEL (QXMIE)
      do I=1,MX
        QXMIE(I) = QAA(KMIE,MIEDX(I))/QAA(IND999,MIEDX(I))
        SSALB(I) = SSA(KMIE,MIEDX(I))
      enddo

c---Combine: Rayleigh scatters & O2 & O3 absorbers to get optical properties
c---values at L1_=L_+1 are a pseudo/climatol layer above the top CTM layer (L_)
        do L = 1,L1U
         TTTX     = TTJ(L)
         call X_interp (TTTX,XQO2, TQQ(1,1),QO2(K,1), TQQ(2,1),QO2(K,2),
     &         TQQ(3,1),QO2(K,3), LQQ(1))
         call X_interp (TTTX,XQO3, TQQ(1,2),QO3(K,1), TQQ(2,2),QO3(K,2),
     &         TQQ(3,2),QO3(K,3), LQQ(2))
        enddo
        XLO3=DO3(L)*XQO3(L)
        XLO2=DM(L)*XQO2(L)*0.20948d0
        XLRAY=DM(L)*QRAYL(K)
        do I=1,MX
          XLAER(I)=AER(I,L)*QXMIE(I)
        enddo
        DTAUX(L)=XLO3+XLO2+XLRAY
        do I=1,MX
          DTAUX(L)=DTAUX(L)+XLAER(I)
        enddo
c  Fractional extinction for Rayleigh scattering and each aerosol type
        PIRAY(L)=XLRAY/DTAUX(L)
        do I=1,MX
          PIAER(I,L)=SSALB(I)*XLAER(I)/DTAUX(L)
        enddo
       enddo !This is for the L=1,L1_
        MFIT = 2*M_
	  J1 = NLBATM
        do J=J1,NB
          do i=1,MFIT
              POMEGAX(I,J) = PIRAY(J)*PAA(I,KMIE,1)
           do k0=1,MX
         POMEGAX(I,J) = POMEGAX(I,J) + PIAER(K0,J)*PAA(I,KMIE,MIEDX(K0))
           enddo
          enddo
      enddo

C-----------------------------------------------------------------------

      call OPMIE (DTAUX,POMEGAX,U0,RFL,AMF2,JXTRA,
     &        AVGF,FJTOP,FJBOT,FSBOT,FJFLX,FLXD,FLXD0, LU)

C-----------------------------------------------------------------------

      do K = 1,W_

c----direct(DIR) and diffuse(FLX) fluxes at top(UP) (solar = negative by convention)
c----     also at bottom (DN), does not include diffuse reflected flux.

        do L = 1,LU
          FFF(K,L) = FFF(K,L) + SOLF*FL(K)*AVGF(L,K)
        enddo
!-----------------------------------------------------------------------
!     BEGINNING OF HEATING RATES CODE (Not relevant to GEOS-Chem)      !
!-----------------------------------------------------------------------
!        FLXUP(K) =  FJTOP(K)
!        DIRUP(K) = -FLXD0(K)
!        FLXDN(K) = -FJBOT(K)
!        DIRDN(K) = -FSBOT(K)
!
!          FREFI = FREFI + SOLF*FL(K)*FLXD0(K)/WL(K)
!          FREFL = FREFL + SOLF*FL(K)*FJTOP(K)/WL(K)
!          FREFS = FREFS + SOLF*FL(K)/WL(K)
!
!c---for each wavelength calculate the flux budget/heating rates:
!c  FLXD(L) = direct flux deposited in layer L  [approx = MU0*(F(L+1) -F(L))]
!c            but for spherical atmosphere!
!c  FJFLX(L) = diffuse flux across top of layer L
!
!c---calculate divergence of diffuse flux in each CTM layer (& t-o-a)
!c---     need special fix at top and bottom:
!c---FABOT = total abs at L.B. &  FXBOT = net diffusive flux at L.B.
!!         FABOT = (1.d0-RFL(K))*(FJBOT(K)+FSBOT(K))
!         FXBOT = -FJBOT(K) + RFL(K)*(FJBOT(K)+FSBOT(K))
!         FLXJ(1) = FJFLX(1,K) - FXBOT
!       do L=2,LU
!         FLXJ(L) = FJFLX(L,K) - FJFLX(L-1,K)
!       enddo
!         FLXJ(LU+1) = FJTOP(K) - FJFLX(LU,K)
!c---calculate net flux deposited in each CTM layer (direct & diffuse):
!         FFX0 = 0.d0
!       do L=1,L1U
!         FFX(K,L) = FLXD(L,K) - FLXJ(L)
!         FFX0 = FFX0 + FFX(K,L)
!      enddo
!
!c  NB: the radiation level ABOVE the top CTM level is included in these budgets
!c      these are the flux budget/heating terms for the column:
!c  FFXNET(K,1) = FLXD0        direct(solar) flux dep into atmos (spherical)
!c  FFXNET(K,2) = FSBOT        direct(solar) flux dep onto LB (surface)
!c  FFXNET(K,3) = FLXD0+FSBOT  TOTAL solar into atmopshere+surface
!c  FFXNET(K,4) = FJTOP        diffuse flux leaving top-of-atmos
!c  FFXNET(K,5) = FFX0         diffuse flux absorbed in atmos
!c  FFXNET(K,6) = FABOT        total (dir+dif) absorbed at LB (surface)
!c       these are surface fluxes to compare direct vs. diffuse:
!c  FFXNET(K,7) = FSBOT        direct flux dep onto LB (surface) - for srf diags
!c  FFXNET(K,8) = FJBOT        diffuse flux dep onto LB (surface)
!
!      FFXNET(K,1) = FLXD0(K)
!      FFXNET(K,2) = FSBOT(K)
!      FFXNET(K,3) = FLXD0(K) + FSBOT(K)
!      FFXNET(K,4) = FJTOP(K)
!      FFXNET(K,5) = FFX0
!      FFXNET(K,6) = FABOT
!      FFXNET(K,7) = FSBOT(K)
!      FFXNET(K,8) = FJBOT(K)
!
!-----------------------------------------------------------------------
!           END OF HEATING RATES CODE (Not relevant to GEOS-Chem)      !
!-----------------------------------------------------------------------
c-----------------------------------------------------------------------
      enddo       ! end loop over wavelength K
c-----------------------------------------------------------------------
!          FREFL = FREFL/FREFS      !calculate reflected flux (energy weighted)
!          FREFI = FREFI/FREFS

c---NB UVB = 280-320 = bins 12:15, UVA = 320-400 = bins 16:17, VIS = bin 18 (++)

C-----------------------------------------------------------------------
      call JRATET(PPJ,TTJ,FFF, VALJXX, LU,NJXU)
C-----------------------------------------------------------------------

c---mapping J-values from fast-JX onto CTM chemistry is done in main

C-----------------------------------------------------------------------

   99 continue

      END SUBROUTINE PHOTO_JX
!
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: opmie
!
! !DESCRIPTION: Subroutine OPMIE is a core Fast-JX scattering subroutine, 
!  specifically for Mie scattering.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE OPMIE (DTAUX,POMEGAX,U0,RFL,AMF2,JXTRA,
     &        FJACT,FJTOP,FJBOT,FSBOT,FJFLX,FLXD,FLXD0, LU)
!
! !USES:
!
!
! !INPUT PARAMETERS:
!
      REAL*8, INTENT(IN) ::   DTAUX(JXL1_,W_),POMEGAX(8,JXL1_,W_)
      REAL*8, INTENT(IN) ::   AMF2(2*JXL1_+1,2*JXL1_+1)
      REAL*8, INTENT(IN) ::   U0,RFL(W_)
      INTEGER, INTENT(IN) ::  JXTRA(JXL2_+1), LU
!
! !OUTPUT VARIABLES:
!
      REAL*8, INTENT(OUT) ::FJACT(JXL_,W_),FJTOP(W_),FJBOT(W_),FSBOT(W_)
      REAL*8, INTENT(OUT) ::  FJFLX(JXL_,W_),FLXD(JXL1_,W_),FLXD0(W_)
!
! !REMARKS:
! 
! !REVISION HISTORY: 
!  28 Mar 2013 - S. D. Eastham - Copied from Fast-JX v7.0
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER JNDLEV(JXL_),JNELEV(JXL1_)
      INTEGER JADDLV(JXL2_+1),JADDTO(JXL2_+1),L2LEV(JXL2_+1)
      INTEGER JTOTL,I,II,J,K,L,LL,IX,JK,   L2,L2L,L22,LZ,LZZ,ND
      INTEGER L1U,L2U,   LZ0,LZ1,LZMID
      REAL*8   SUMT,SUMJ

      REAL*8  DTAU(JXL1_+1,W_),POMEGAJ(M2_,JXL2_+1,W_),TTAU(JXL2_+1,W_)
      REAL*8  FTAU2(JXL2_+1,W_),POMEGAB(M2_,W_)
      REAL*8  ATAUA,ATAUZ,XLTAU,TAUDN,TAUUP,DTAUJ,FJFLX0
      REAL*8, DIMENSION(W_) :: TAUBTM,TAUTOP,FBTM,FTOP,ZFLUX
c--- variables used in mie code-----------------------------------------
      REAL*8, DIMENSION(W_)         :: FJT,FJB
      REAL*8, DIMENSION(N_,W_)      :: FJ,FZ,ZTAU
      REAL*8, DIMENSION(M2_,N_,W_)  :: POMEGA
      REAL*8, DIMENSION(2*JXL1_,W_)  :: FLXD2

      !=================================================================
      ! OPMIE begins here!
      !=================================================================

c---there is a parallel correspondence:
c  dimension of JX arrays JXL_ .ge. dimension that CTM is using = L_
c  but calculation is done for L_=LU, L1_=L1U, L2_=L2U lengths of CTM
c
C  fast-J Mie code for J_s, only uses 8-term expansion, 4-Gauss pts
c
c in:
c     DTAUX(1:L1_,1:W_) = optical depth of each layer
c     POMEGAX(1:8,1:L1_,1:W_) = scattering phase fn (multiplied by s-s albedo)
c     U0  = cos (SZA)
c     RFL(1:W_) = Lambertian albedo of surface
c     AMF2(1:2*L1_+1,1:2*L1_+1) = air mass factor (I,L)=wt of layer-I to layer-L
c        AMF2 now does both edges and middle of CTM layers
c     JXTRA(1:L1_) = number 0:J = no. of additional levels to be inserted
c out:
c     FJACT(1:L_,1:W_) = mean actinic flux(diff+direct) at std CTM levels(mid-lyr)
c  (new ver 5.7 diagnostics for fluxes, deposition)  fluxes 'down' are <0
c     FJTOP(1:W_) = diffuse flux out top-of-atmosphere (TAU=0 above top model lyr)
c     FJBOT(1:W_) = diffuse flux onto surface (<0 by definition)
c     FSBOT(1:W_) = direct/solar flux onto surface  (<0 by definition)
c     FJFLX(1:L_,1:W_) = diffuse flux across top of model layer L
C        this connects with FJBOT = FJFLX(0) & FJTOP = FJFLX(L_+1) (not dim!!)
c     FLXD(1:L_+1,1:W_) = solar flux deposited in layer L (includes lyr above CTM)
c        this should take into account sphericity, and is not just = mu0
c     FLXD0(1:W_) = sum of solar flux deposited in atmos
c        does NOT include flux on lower surface, does NOT mean absorbed!
C-----------------------------------------------------------------------
c
c     DTAU     Local optical depth of each CTM level
c     TTAU     Optical depth of air vertically above each point (to top of atm)
c     FTAU2     Attenuation of solar beam
c     POMEGAJ  Scattering phase function
c
c---new ver 5.3 code adds sub-layers (# = JXTRA(L2)) using ATAU as the
c   factor increase from sub-layer to sub-layer
c
C---------------------SET UP FOR MIE CODE-------------------------------
c
c-----------------wavelength independent--------------------------------
c
c  Transpose the ascending TTAU grid to a descending ZTAU grid.
c  Double the resolution - TTAU points become the odd points on the
c  ZTAU grid, even points needed for asymm phase fn soln, contain 'h'.
c  Odd point added at top of grid for unattenuated beam   (Z='inf')
c
c  The following mapping holds for JADDLV=0
c        Surface:   TTAU(1)    ==> ZTAU(2*L2_+1)
c        Top:       TTAU(L2_)  ==> ZTAU(3)
c        Infinity:     0.0     ==> ZTAU(1)
c        index: 2*(L2_+1-L2)+1 ==> LZ
c
c  Mie scattering code only used from surface to level L2_
C------------------------------------------------------------------------
c
C------------------------------------------------------------------------
c  Insert new levels, working downwards from the top of the atmosphere
c  to the surface (down in 'LZ', up in 'L2'). This allows ztau and pomega
c  to be incremented linearly, and the flux fz to be attenuated top-down
c    (avoiding problems where lower level fluxes are zero).
C------------------------------------------------------------------------
c
c  Ascend through atmosphere transposing grid and adding extra points
c  remember L2=1 is surface of CTM, but last layer (LZ) in scattering code.
c  there are twice the number of layers in the LZ arrays (2*L2_ + 2*JADDTO + 1)
c    because we need to insert the intermediate layers (even LZ) for the
c    asymmetric scattering code.
c
c  Transfer the L2=1:L2_+1 values (TTAU,FTAU2,POMEGAJ) onto the reverse
c    order, expanded, doubled-level scatter grid.
c    Note that we need to deal with the expansion by JADD levels (L2L).
c      These JADDLV levels are skipped and need to be interpolated later.
c    Note that only odd LZ levels are filled,
c
C----------------------re-grid data---------------------------------------------
c  Calculate cumulative total and define levels we want J-values at.
c  Sum upwards for levels, and then downwards for Mie code readjustments.
c
c     JXTRA(L2)  Number of new levels to add between (L2) and (L2+1)
c           ***JXTRA(1:L2_+1) is calculated based on the aerosol+cloud OD_s
c     JADDLV(L2)  Number of new levels actually added at each wavelength
c            where JADDLV = 0 when there is effectively no FTAU2
c     JADDTO(L2)   Total number of new levels to add to and above level (L2)
c     JNDLEV(L) = L2 index that maps on CTM mid-layer L
c
c---JADDLV(L2=1:L2_) = number of levels to add between TTAU2(L2) and TTAU(L2+1)
c---    JADDLV is taken from JXTRA, which is based on visible OD.
c---    JADDTO(L2=1:L2_+1) is the cumulative number of levels to be added
c---these should be fixed for all wavelengths to lock-in the array sizes

      if (LU .gt. JXL_) then
        call EXITC (' OPMIE:  JXL_ .lt. L_')
      endif

      L1U = LU + 1
      L2U = 2*LU + 2


      do L2 = 1,L2U,1
        JADDLV(L2) = JXTRA(L2)
      enddo
        JADDTO(L2U+1) = 0
      do L2 = L2U,1,-1
        JADDTO(L2) = JADDTO(L2+1) + JADDLV(L2)
      enddo

c---expanded grid now included CTM edge and mid layers plus expanded
c---    grid to allow for finer delta-tau at tops of clouds.
c---    DIM of new grid = L2U + JADDTO(1) + 1

c---L2LEV(L2) = L2-index for old level L2 in expanded J-grid (w/JADDLV)
c     in absence of JADDLV, L2LEV(L2) = L2
        L2LEV(1)  = 1
      do L2 = 2,L2U+1
        L2LEV(L2) = L2LEV(L2-1) + 1 + JADDLV(L2-1)
      enddo

c---JNDLEV(L=1:L_) = L2-index in expanded grid for CTM mid-layer L
c---JNELEV(L=1:L_) = L2-index for top of layer L
      do L = 1,LU
        JNDLEV(L) = L2LEV(2*L)
        JNELEV(L) = L2LEV(2*L+1)
      enddo
        JNELEV(LU+1) = 0  !need to set this to top-of-atmosphere

      ND = 2*L2U + 2*JADDTO(1) + 1

      if(ND .gt. N_) then
        call EXITC (' overflow of scatter arrays: ND > N_')
      endif

c----------------begin wavelength dependent set up------------------------------

C---Reinitialize arrays
      ZTAU(:,:)     = 0.d0
      FZ(:,:)       = 0.d0
      POMEGA(:,:,:) = 0.d0

      do K=1,W_

C---Set up optical depth DTAU(L)
       do L = 1,L1U
        DTAU(L,K) = DTAUX(L,K)
       enddo
        DTAU(L1U+1,K) = 0.d0

c---Define the total scattering phase fn for each CTM layer L=1:L_+1
c---   from a DTAU-wt_d mix of aerosols, cloud & Rayleigh
C---No. of quadrature pts fixed at 4(M_), expansion of phase fn @ 8
       do L = 1,L1U
        do I = 1,M2_
          POMEGAJ(I,L,K) = POMEGAX(I,L,K)
        enddo
       enddo

C---Calculate attenuated incident beam exp(-TTAU/U0 = DTAU * AirMassFactor)
c---      at the middle & edges of the CTM layers L=1:2*L1_+1
c---  L1_ is top-edge of CTM (ie, L=38 = 2 hPa) which has TAU > 0
c---  note that DTAU(L1_) is optical depth in the FULL CTM layer just above
        FTAU2(:,:) = 0.d0
        FTAU2(L2U+1,:) = 1.0d0
       do LL = 1,2*L1U+1
         L = (LL+1)/2
        if (AMF2(LL,LL) .gt. 0.0d0) then
           XLTAU = 0.0d0
         do II = 1,2*L1U+1
           I = (II+1)/2
           XLTAU = XLTAU + 0.5d0*DTAU(I,K)*AMF2(II,LL)
         enddo
         if (XLTAU .lt. 76.d0) then   ! zero out flux at 1e-33
          FTAU2(LL,K) = exp(-XLTAU)
         endif
        endif
       enddo

c---calculate direct solar flux deposited in each CTM half-layer: L=1:L2_
c---     use FSBOT for surface flux, cannot do layer above CTM (L_+1)
          FLXD2(:,:) = 0.d0
       do LL = 1,2*L1U
        if (AMF2(LL,LL) .gt. 0.d0) then
          FLXD2(LL,K) = (FTAU2(LL+1,K) - FTAU2(LL,K))/AMF2(LL,LL)
        endif
       enddo
        if (AMF2(1,1) .gt. 0.d0) then
          FSBOT(K) = FTAU2(1,K)/AMF2(1,1)
        else
          FSBOT(K) = 0.d0
        endif

       do LL = 2,2*L1U,2
         L=LL/2
         FLXD(L,K) = FLXD2(LL,K)+FLXD2(LL-1,K)
       enddo

c---integrate solar flux depositied in CTM layers L=1:L_, cannot do top layer
c---  note FLXD0 .ne. (1.d0 - FTAU(L_+1))/AMF(L_+1,L_+1) with spherical atmos
        FLXD0(K) = 0.d0
       if (AMF2(2*L1U,2*L1U) .gt. 0.d0) then
        do L=1,L1U
         FLXD0(K) = FLXD0(K) + FLXD(L,K)
        enddo
       endif

C------------------------------------------------------------------------
c  Take optical properties on CTM layers and convert to a photolysis
c  level grid corresponding to layer centres and boundaries. This is
c  required so that J-values can be calculated for the centre of CTM
c  layers; the index of these layers is kept in the JNDLEV array.
C------------------------------------------------------------------------
c---Now combine the CTM layer edges (1:L_+2) with the CTM mid-layer
c---    points (1:L_) plus 1 for the mid point of added top layer.
c---combine these edge- and mid-layer points into grid of size:
c---              L2_+1 = 2*L1_+1 = 2*L_+3
c---calculate column optical depths above each level, TTAU(1:L2_+1)
c---      note that TTAU(L2_+1)=0 and TTAU(1)=total OD

        TTAU(L2U+1,K) = 0.0d0
       do L2 = L2U,1,-1
        L          = (L2+1)/2
        DTAUJ      = 0.5d0 * DTAU(L,K)
        TTAU(L2,K)   = TTAU(L2+1,K) + DTAUJ
       enddo

c----solar flux incident on lower boundary & Lambertian reflect factor:
       if (FSBOT(K) .gt. 0.d0) then
        ZFLUX(K) = FSBOT(K)*RFL(K)/(1.d0+RFL(K))
       else
        ZFLUX(K) = 0.d0
       endif

c  Calculate scattering properties, level centres then level boundaries
c>>>>>be careful of order, we are overwriting/shifting the 'POMEGAJ' upward in index
       do L2 = L2U,2,-2
        L   = L2/2
        do I = 1,M2_
          POMEGAJ(I,L2,K) = POMEGAJ(I,L,K)
        enddo
       enddo
c---lower boundary value is set (POMEGAJ(I,1)), but set upper:
       do I = 1,M2_
         POMEGAJ(I,L2U+1,K) = POMEGAJ(I,L2U,K)
       enddo
c---now have POMEGAJ filled at even points from L2=3:L2_-1
c---use inverse interpolation for correct tau-weighted values at edges
       do L2 = 3,L2U-1,2
        TAUDN = TTAU(L2-1,K)-TTAU(L2,K)
        TAUUP = TTAU(L2,K)-TTAU(L2+1,K)
        do I = 1,M2_
          POMEGAJ(I,L2,K) = (POMEGAJ(I,L2-1,K)*TAUDN +
     &           POMEGAJ(I,L2+1,K)*TAUUP) / (TAUDN+TAUUP)
        enddo
       enddo

C---at this point FTAU2(1:L2_+1) and POMEAGJ(1:8, 1:L2_+1)
c---    where FTAU2(L2_+1) = 1.0 = top-of-atmos, FTAU2(1) = surface

       do L2 = 1,L2U+1          ! L2 = index of CTM edge- and mid-layers
        L2L = L2LEV(L2)        ! L2L = index for L2 in expanded scale(JADD)
        LZ  = ND + 2 - 2*L2L  ! LZ = index for L2 in scatt arrays
          ZTAU(LZ,K) = TTAU(L2,K)
          FZ(LZ,K)   = FTAU2(L2,K)
        do I=1,M2_
          POMEGA(I,LZ,K) = POMEGAJ(I,L2,K)
        enddo
       enddo

c   Now go thru the pairs of L2 levels to see if we need JADD levels
       do L2 = 1,L2U             ! L2 = index of CTM edge- and mid-layers
         L2L = L2LEV(L2)         ! L2L = index for L2 in expanded scale(JADD)
         LZ  = ND + 2 - 2*L2L   ! LZ = index for L2 in scatt arrays
         L22 = L2LEV(L2+1) - L2LEV(L2) - 1   ! L22 = 0 if no added levels

        if (L22 .gt. 0) then
          TAUBTM(K) = TTAU(L2,K)
          TAUTOP(K) = TTAU(L2+1,K)
          FBTM(K)   = FTAU2(L2,K)
          FTOP(K)   = FTAU2(L2+1,K)
         do I = 1,M2_
          POMEGAB(I,K) = POMEGAJ(I,L2,K)
         enddo

c---to fit L22 new layers between TAUBOT > TAUTOP, calculate new 1/ATAU factor
c---  such that TAU(just above TAU-btm) = ATUAZ * TAUBTM < TAUBTM
         ATAUZ = exp(-log(TAUBTM(K)/max(TAUTOP(K),ATAU0))/float(L22+1))
         do L = 1,L22           ! add odd levels between L2LEV(L2) & L2LEV(L2+1)
          LZZ = LZ - 2*L       ! LZZ = index(odd) of added level in scatt arrays
          ZTAU(LZZ,K) = TAUBTM(K) * ATAUZ

c---fraction from TAUBTM=>TAUTOP
          ATAUA=(TAUBTM(K)-ZTAU(LZZ,K))/(TAUBTM(K)-TAUTOP(K))
c---solar flux at interp-levels: use exp(TAU/U0) if U0>0.02 (89 deg),
c---else scale by TAU
          if (U0 .gt. 0.02d0) then
            FZ(LZZ,K) = FTOP(K) * exp((TAUTOP(K)-ZTAU(LZZ,K))/U0)
          else
            if (FBTM(K) .lt. 1.d-32) then
              FZ(LZZ,K) = 0.d0
            else
              FZ(LZZ,K) = FBTM(K) * (FTOP(K)/FBTM(K))**ATAUA
            endif
          endif
          do I = 1,M2_
            POMEGA(I,LZZ,K) = POMEGAB(I,K) +
     &               ATAUA*(POMEGAJ(I,L2+1,K)-POMEGAB(I,K))
          enddo
            TAUBTM(K)    = ZTAU(LZZ,K)
            FBTM(K)      = FZ(LZZ,K)
          do I = 1,M2_
            POMEGAB(I,K) = POMEGA(I,LZZ,K)
          enddo
         enddo
        endif
       enddo

c   Now fill in the even points with simple interpolation in scatter arrays:
       do LZ = 2,ND-1,2
         ZTAU(LZ,K) = 0.5d0*(ZTAU(LZ-1,K)+ZTAU(LZ+1,K))
         FZ(LZ,K)   = sqrt(FZ(LZ-1,K)*FZ(LZ+1,K))
        do I=1,M2_
         POMEGA(I,LZ,K) = 0.5d0*(POMEGA(I,LZ-1,K)+POMEGA(I,LZ+1,K))
        enddo
       enddo

      enddo  ! wavelength loop!

C-----------------------------------------------------------------------
       call MIESCT(FJ,FJT,FJB,POMEGA,FZ,ZTAU,ZFLUX,RFL,U0,ND)
C-----------------------------------------------------------------------

c---Move mean intensity from scatter array FJ(LZ=1:ND)
c---              to CTM mid-level array FJACT(L=1:L_)

      do K=1,W_

c---mean intensity at mid-layer:  4*<I> + solar
c       do L = 1,LU
c        L2L = JNDLEV(L)
c        LZ  = ND+2 - 2*L2L
c        FJACT(L,K) = 4.d0*FJ(LZ,K) + FZ(LZ,K)
c       enddo

c---mean intensity averaged throughout layer:
       do L = 1,LU
         LZ0 = ND+2 - 2*JNELEV(L)
        if (L .gt. 1) then
         LZ1 = ND+2 - 2*JNELEV(L-1)
        else
         LZ1 = ND
        endif
         SUMJ = (4.d0*FJ(LZ0,K)+FZ(LZ0,K))*(ZTAU(LZ0+2,K)-ZTAU(LZ0,K))
     &         +(4.d0*FJ(LZ1,K)+FZ(LZ1,K))*(ZTAU(LZ1,K)-ZTAU(LZ1-2,K))
         SUMT = ZTAU(LZ0+2,K)-ZTAU(LZ0,K) + ZTAU(LZ1,K)-ZTAU(LZ1-2,K)

        do LZ = LZ0+2,LZ1-2,2
         SUMJ =SUMJ+(4.d0*FJ(LZ,K)+FZ(LZ,K))*(ZTAU(LZ+2,K)-ZTAU(LZ-2,K))
         SUMT =SUMT + ZTAU(LZ+2,K)-ZTAU(LZ-2,K)
        enddo
        FJACT(L,K) = SUMJ/SUMT
       enddo

c---mean diffuse flux:  4<I*mu> (not solar) at top of layer L
c---      average (tau-wtd) the h's just above and below the L-edge
       do L = 1,LU
        L2L = JNELEV(L)
        LZ  = ND+2 - 2*L2L
        FJFLX0 = (ZTAU(LZ+1,K)-ZTAU(LZ,K))/(ZTAU(LZ+1,K)-ZTAU(LZ-1,K))
        FJFLX(L,K)=4.d0*(FJ(LZ-1,K)*FJFLX0 + FJ(LZ+1,K)*(1.d0-FJFLX0))
       enddo

c---diffuse fluxes reflected at top, incident at bottom
         FJTOP(K) = FJT(K)
         FJBOT(K) = FJB(K)

      enddo  ! wavelength loop!

      END SUBROUTINE OPMIE
!
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: miesct
!
! !DESCRIPTION: Subroutine MIESCT is an adaptation of the Prather radiative
!  transfer code (mjp, 10/95).
!     Prather, 1974, Astrophys. J. 192, 787-792.
!         Solution of inhomogeneous Rayleigh scattering atmosphere.
!         (original Rayleigh w/ polarization)
!     Cochran and Trafton, 1978, Ap.J., 219, 756-762.
!         Raman scattering in the atmospheres of the major planets.
!         (first use of anisotropic code)
!     Jacob, Gottlieb and Prather, 1989, J.Geophys.Res., 94, 12975-13002.
!         Chemistry of a polluted cloudy boundary layer,
!         (documentation of extension to anisotropic scattering)
!
!    takes atmospheric structure and source terms from std J-code
!    ALSO limited to 4 Gauss points, only calculates mean field! (M=1)
!   --------------------------------------------------------------------
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE MIESCT(FJ,FJT,FJB, POMEGA,FZ,ZTAU,ZFLUX,RFL,U0,ND)
!
! !USES:
!
!
! !INPUT PARAMETERS:
!
      INTEGER, INTENT(IN) ::  ND
      REAL*8, INTENT(IN)  ::  POMEGA(M2_,N_,W_),FZ(N_,W_),ZTAU(N_,W_)
     &                       ,RFL(W_),U0,ZFLUX(W_)
!
! !OUTPUT VARIABLES:
!
      REAL*8, INTENT(OUT) ::  FJ(N_,W_),FJT(W_),FJB(W_)
!
! !REMARKS:
! 
! !REVISION HISTORY: 
!  28 Mar 2013 - S. D. Eastham - Copied from Fast-JX v7.0
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      REAL*8  PM(M_,M2_),PM0(M2_)
      INTEGER I, IM  ,K

      !=================================================================
      ! MIESCT begins here!
      !=================================================================

      do I = 1,M_
       call LEGND0 (EMU(I),PM0,M2_)
       do IM = 1,M2_
         PM(I,IM) = PM0(IM)
       enddo
      enddo

       call LEGND0 (-U0,PM0,M2_)
       do IM=1,M2_
         PM0(IM) = 0.25d0*PM0(IM)
       enddo

c---BLKSLV now called with all the wavelength arrays (K=1:W_)

      call BLKSLV(FJ,POMEGA,FZ,ZTAU,ZFLUX,RFL,PM,PM0,FJT,FJB, ND)

      END SUBROUTINE MIESCT
!
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: legnd0
!
! !DESCRIPTION: Subroutine LEGND0 calculates ordinary Legendre functions 
!  of X (real) from $P[0] = PL(1) = 1, P[1] = X, \dots, P[N-1] = PL(N)$
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE LEGND0 (X,PL,N)
!
! !USES:
!
!
! !INPUT PARAMETERS:
!
      INTEGER, INTENT(IN) :: N
      REAL*8, INTENT(IN)  :: X
!
! !OUTPUT VARIABLES:
!
      REAL*8, INTENT(OUT) :: PL(N)
!
! !REMARKS:
! 
! !REVISION HISTORY: 
!  28 Mar 2013 - S. D. Eastham - Copied from Fast-JX v7.0
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER I
      REAL*8  DEN

      !=================================================================
      ! LEGND0 begins here!
      !=================================================================

C---Always does PL(2) = P[1]
        PL(1) = 1.d0
        PL(2) = X
        do I = 3,N
         DEN = (I-1)
         PL(I) = PL(I-1)*X*(2.d0-1.0/DEN) - PL(I-2)*(1.d0-1.d0/DEN)
        enddo

      END SUBROUTINE LEGND0
!
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: set_prof
!
! !DESCRIPTION: Subroutine SET\_PROF sets up the atmospheric column profile
!  for Fast-JX.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE SET_PROF( NLON, NLAT, YLAT,  MONTH,   DAY, 
     &              PCOL, TCOL,SFCALB,ODCOL, OPTDUST, OPTAER,
     &              TINTERP, DO3, DBC, DM, Z, AER )
!
! !USES:
!
      USE TOMS_MOD, ONLY : TOMS, DTOMS1, DTOMS2
      USE ERROR_MOD
!
! !INPUT PARAMETERS:
!
      INTEGER, INTENT(IN)    :: DAY, MONTH, NLAT, NLON
      REAL*8,  INTENT(IN)    :: SFCALB, YLAT
      REAL*8,  INTENT(IN)    :: PCOL(LLPAR+1),TCOL(LLPAR)
      REAL*8,  INTENT(INOUT) :: ODCOL(LLPAR)
      REAL*8,  INTENT(IN)    :: OPTDUST(LLPAR,NDUST)  
      REAL*8,  INTENT(IN)    :: OPTAER(LLPAR,NAER*NRH)  
!
! !OUTPUT VARIABLES:
!
      REAL*8,  INTENT(OUT)   :: TINTERP(LLPAR)
      REAL*8,  INTENT(OUT)   :: DO3(LLPAR)
      REAL*8,  INTENT(OUT)   :: DBC(LLPAR)
      REAL*8,  INTENT(OUT)   :: DM(LLPAR)
      REAL*8,  INTENT(OUT)   :: Z(LLPAR+1)
      REAL*8, INTENT(INOUT)  :: AER(A_,LLPAR)
!
! !REMARKS:
!  (1 ) Since we parallelize over columns, T, ODCOL, OPTDUST, and OPTAER
!        are 1-D vectors. In the original code from Oliver Wild, these were 
!        3-D arrays.  Also P and SA are just scalars since we just pass one 
!        surface location at a time w/in the parallel loop. (bmy, 9/13/99)
!  (2 ) Mineral dust profiles are also constructed (rvm, 06/04/00)
!  (3 ) Other aerosol profiles are also constructed (rvm, bmy, 2/27/02)
!  (4 ) Added NLON, NLAT, DAY to the arg list.  Now weight the O3 column by 
!        the observed monthly mean EP-TOMS data.  Also updated comments and 
!        added standard GEOS-CHEM documentation header. (mje, bmy, 7/13/03)
! 
! !REVISION HISTORY: 
!  28 Mar 2013 - S. D. Eastham - Addapted from J. Mao's implementation of
!                                Fast-JX v6.2 in GEOS-Chem
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER                :: I, K, L, M, N
      REAL*8                 :: DLOGP,F0,T0,B0,PB,PC,XC,MASFAC,SCALEH
      REAL*8                 :: PSTD(52),OREF2(51),TREF2(51),BREF2(51)
      REAL*8                 :: O3COL_DU, DAYTOMS

      !=================================================================
      ! SET_PROF begins here!
      !=================================================================

      ! Set up cloud and surface properties
      CALL CLDSRF( ODCOL, SFCALB )

      ! Zero out aerosol column
      DO K=1,A_
        DO I=1,LLPAR
          AER(K,I) = 0.D0
        ENDDO
      ENDDO

      !=================================================================      
      ! Set up pressure levels for O3/T climatology - assume that value
      ! given for each 2 km z* level applies from 1 km below to 1 km 
      ! above, so select pressures at these boundaries. Surface level 
      ! values at 1000 mb are assumed to extend down to the actual 
      ! P(nslon,nslat).
      !=================================================================      
      PSTD(1) = MAX( PCOL(1), 1000.D0 )
      PSTD(2) = 1000.D0 * 10.D0**( -1.D0/16.D0 )
      DLOGP   = 10.D0**( -2.D0/16.D0 )
      DO I = 3, 51
         PSTD(I) = PSTD(I-1) * DLOGP
      ENDDO
      PSTD(52) = 0.D0

      ! Mass factor - delta-Pressure [hPa] to delta-Column [molec/cm2]
      MASFAC = 100.D0 * 6.022D+23 / ( 28.97D0 * 9.8D0 * 10.D0 )

      ! Select appropriate monthly and latitudinal profiles
      ! Now use YLAT instead of Oliver's YDGRD(NSLAT) (bmy, 9/13/99) 
      M = MAX( 1, MIN( 12, MONTH                   ) )
      L = MAX( 1, MIN( 18, ( INT(YLAT) + 99 ) / 10 ) )

      ! Temporary arrays for climatology data
      DO I = 1, 51
          OREF2(I) = OREF(I,L,M)
          TREF2(I) = TREF(I,L,M)
          BREF2(I) = BREF(I)
      ENDDO

      ! Apportion O3 and T on supplied climatology z* levels onto CTM levels 
      ! with mass (pressure) weighting, assuming constant mixing ratio and
      ! temperature half a layer on either side of the point supplied.

      DO I = 1, NB
         F0 = 0.D0
         T0 = 0.D0
         B0 = 0.D0
         DO K = 1, 51
            PC = MIN( PCOL(I),   PSTD(K)   )
            PB = MAX( PCOL(I+1), PSTD(K+1) )
            IF ( PC .GT. PB ) THEN
               XC = ( PC - PB ) / ( PCOL(I) - PCOL(I+1) )
               F0 = F0 + OREF2(K)*XC
               T0 = T0 + TREF2(K)*XC
               B0 = B0 + BREF2(K)*XC
            ENDIF
         ENDDO
         TINTERP(I)  = T0
         DO3(I) = F0 * 1.D-6
         DBC(I) = B0
      ENDDO

      !=================================================================
      ! Calculate effective altitudes using scale height at each level
      !=================================================================
      Z(1) = 0.D0
      DO I = 1, LLPAR
         SCALEH = 1.3806D-19 * MASFAC * TINTERP(I)
         Z(I+1) = Z(I) - ( LOG( PCOL(I+1) / PCOL(I) ) * SCALEH )
      ENDDO

      !=================================================================
      ! Add Aerosol Column - include aerosol types here. Currently use 
      ! soot water and ice; assume black carbon x-section of 10 m2/g, 
      ! independent of wavelength; assume limiting temperature for 
      ! ice of -40 deg C.
      !=================================================================
      DO I = 1, LLPAR
         ! Turn off uniform black carbon profile (rvm, bmy, 2/27/02)
         !AER(1,I) = DBC(I) * 10.D0 * ( Z(I+1) - Z(I) )
         AER(1,I) = 0D0

         IF ( TCOL(I) .GT. 233.D0 ) THEN
            AER(2,I) = ODCOL(I)
            AER(3,I) = 0.D0
         ELSE
            AER(2,I) = 0.D0
            AER(3,I) = ODCOL(I)
         ENDIF   

         ! Also add in aerosol optical depth columns (rvm, bmy, 9/30/00)
         DO N = 1, NDUST
            AER(3+N,I) = OPTDUST(I,N)
         ENDDO
        
         ! Also add in other aerosol optical depth columns (rvm, bmy, 2/27/02)
         DO N = 1, NAER*NRH
            AER(3+N+NDUST,I) = OPTAER(I,N)
         ENDDO

      ENDDO

      DO K = 1, A_
         AER(K,LLPAR+1) = 0.D0
      ENDDO

      !=================================================================
      ! Calculate column quantities for Fast-JX
      !=================================================================

      O3COL_DU = 0d0

      DO I = 1, NB

         ! Monthly mean air Column [molec/cm2]
         DM(I)  = ( PCOL(I) - PCOL(I+1) ) * MASFAC

         ! Monthly mean O3 column [molec/cm2]
         DO3(I) = DO3(I) * DM(I)
    
         ! Monthly mean O3 column [DU] 
         O3COL_DU = O3COL_DU + ( DO3(I) / 2.69d16 )
      ENDDO
      
      !=================================================================
      ! Now weight the O3 column by the observed monthly mean TOMS.
      ! Missing data is denoted by the flag -999. (mje, bmy, 7/15/03)
      ! 
      ! TOMS/SBUV MERGED TOTAL OZONE DATA, Version 8, Revision 3.
      ! Resolution:  5 x 10 deg.
      !
      ! Methodology (bmy, 2/12/07)
      ! ----------------------------------------------------------------
      ! FAST-J comes with its own default O3 column climatology (from 
      ! McPeters 1992 & Nagatani 1991), which is stored in the input 
      ! file "jv_atms.dat".  These "FAST-J default" O3 columns are used 
      ! in the computation of the actinic flux and other optical 
      ! quantities for the FAST-J photolysis.  
      !
      ! The TOMS/SBUV O3 columns and 1/2-monthly O3 trends (contained 
      ! in the TOMS_200701 directory) are read into GEOS-Chem by routine 
      ! READ_TOMS in "toms_mod.f".  Missing values (i.e. locations where 
      ! there are no data) in the TOMS/SBUV O3 columns are defined by 
      ! the flag -999.  
      ! 
      ! After being read from disk in routine READ_TOMS, the TOMS/SBUV 
      ! O3 data are then passed to the FAST-J routine "set_prof.f".  In 
      ! "set_prof.f", a test is done to make sure that the TOMS/SBUV O3 
      ! columns and 1/2-monthly trends do not have any missing values 
      ! for (lat,lon) location for the given month.  If so, then the 
      ! TOMS/SBUV O3 column data is interpolated to the current day and 
      ! is used to weight the "FAST-J default" O3 column.  This 
      ! essentially "forces" the "FAST-J default" O3 column values to 
      ! better match the observations, as defined by TOMS/SBUV.
      !
      ! If there are no TOMS/SBUV O3 columns (and 1/2-monthly trends) 
      ! at a (lat,lon) location for given month, then FAST-J will revert 
      ! to its own "default" climatology for that location and month.  
      ! Therefore, the TOMS O3 can be thought of as an  "overlay" data 
      ! -- it is only used if it exists.
      !
      ! Note that there are no TOMS/SBUV O3 columns at the higher 
      ! latitudes.  At these latitudes, the code will revert to using 
      ! the "FAST-J default" O3 columns.
      !
      ! As of February 2007, we have TOMS/SBUV data for 1979 thru 2005.  
      ! 2006 TOMS/SBUV data is incomplete as of this writing.  For years
      ! 2006 and onward, we use 2005 TOMS O3 columns.
      !
      ! This methodology was originally adopted by Mat Evans.  Symeon 
      ! Koumoutsaris was responsible for creating the downloading and 
      ! processing the TOMS O3 data files from 1979 thru 2005 in the 
      ! TOMS_200701 directory.
      !=================================================================
      DAYTOMS = 0d0

      IF ( DAY <= 15 ) THEN 

         ! Interpolate O3 to current day (w/in first half of month)
         IF ( TOMS(NLON,NLAT)   > -999d0  .AND.
     &        DTOMS1(NLON,NLAT) > -999d0 ) THEN  
            DAYTOMS = TOMS(NLON,NLAT) + DTOMS1(NLON,NLAT) * ( DAY - 15 )
         ENDIF

      ELSE

         ! Interpolate O3 to current day (w/in 2nd half of month)
         IF ( TOMS(NLON,NLAT)   > -999d0  .AND.
     &        DTOMS2(NLON,NLAT) > -999d0 ) THEN  
            DAYTOMS = TOMS(NLON,NLAT) + DTOMS2(NLON,NLAT) * ( DAY - 15 )
         ENDIF

      ENDIF
      
      ! Scale monthly O3 profile to the daily O3 profile (if available)
      IF ( DAYTOMS > 0d0 ) THEN 
         DO I = 1, NB
            DO3(I) = DO3(I) * ( DAYTOMS / O3COL_DU )
         ENDDO
      ENDIF

      END SUBROUTINE SET_PROF
!
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: cldsrf
!
! !DESCRIPTION: Subroutine CLDSRF calculates cloud surface properties.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE CLDSRF( ODCOL, SFCALB )
!
! !USES:
!
!
! !INPUT PARAMETERS:
!
      REAL*8, INTENT(INOUT) :: ODCOL(LLPAR)
      REAL*8, INTENT(IN)    :: SFCALB
!
! !OUTPUT VARIABLES:
!
!
! !REMARKS:
! 
! !REVISION HISTORY: 
!  28 Mar 2013 - S. D. Eastham - Adapted from Fast-JX v6.2/GC implementation
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER I
      REAL*8  ODSUM, ODMAX, ODTOT

      !=================================================================
      ! CLDSRF begins here!
      !=================================================================
c
c Default lower photolysis boundary as bottom of level 1
!      nlbatm = 1
c
c Scale optical depths as appropriate - limit column to 'odmax'
      odmax = 200.d0
      odsum =   0.d0
      do i=1,llpar
        odcol(i) = dble(odcol(i))
        odsum    = odsum + odcol(i)
      enddo
      if(odsum.gt.odmax) then
        odsum = odmax/odsum
        do i=1,llpar
          odcol(i) = odcol(i)*odsum
        enddo
        odsum = odmax
      endif

      return
      END SUBROUTINE CLDSRF
!
!EOC
!!------------------------------------------------------------------------------
!!          Harvard University Atmospheric Chemistry Modeling Group            !
!!------------------------------------------------------------------------------
!!BOP
!!
!! !IROUTINE: srtn
!!
!! !DESCRIPTION: Subroutine SRTN is a reference subroutine, with more than one
!!  line of description.
!!\\
!!\\
!! !INTERFACE:
!!
!      SUBROUTINE SRTN (INARG1,IOARG2,OARG3)
!!
!! !USES:
!!
!      USE THIS_MOD, ONLY : THIS_FN
!!
!! !INPUT PARAMETERS:
!!
!      REAL*8, INTENT(IN)          :: INARG1     ! Description of INARG1
!      REAL*8, INTENT(INOUT)       :: IOARG2     ! Description of IOARG2
!!
!! !OUTPUT VARIABLES:
!!
!      REAL*8, INTENT(OUT)         :: OARG3      ! Description of OARG3
!!
!! !REMARKS:
!! (1) A remark
!! 
!! !REVISION HISTORY: 
!!  28 Mar 2013 - S. D. Eastham - Copied from Fast-JX v7.0
!!EOP
!!------------------------------------------------------------------------------
!!BOC
!!
!! !LOCAL VARIABLES:
!!
!      INTEGER, PARAMETER    :: IPARAM1=2
!
!      !=================================================================
!      ! SRTN begins here!
!      !=================================================================
!
!      OARG3 = INARG1 + IOARG2
!      IOARG2 = OARG3/2.d0*IPARAM1
!
!      END SUBROUTINE SRTN
!!
!!EOC
!!------------------------------------------------------------------------------
!!          Harvard University Atmospheric Chemistry Modeling Group            !
!!------------------------------------------------------------------------------
!!BOP
!!
!! !IROUTINE: srtn
!!
!! !DESCRIPTION: Subroutine SRTN is a reference subroutine, with more than one
!!  line of description.
!!\\
!!\\
!! !INTERFACE:
!!
!      SUBROUTINE SRTN (INARG1,IOARG2,OARG3)
!!
!! !USES:
!!
!      USE THIS_MOD, ONLY : THIS_FN
!!
!! !INPUT PARAMETERS:
!!
!      REAL*8, INTENT(IN)          :: INARG1     ! Description of INARG1
!      REAL*8, INTENT(INOUT)       :: IOARG2     ! Description of IOARG2
!!
!! !OUTPUT VARIABLES:
!!
!      REAL*8, INTENT(OUT)         :: OARG3      ! Description of OARG3
!!
!! !REMARKS:
!! (1) A remark
!! 
!! !REVISION HISTORY: 
!!  28 Mar 2013 - S. D. Eastham - Copied from Fast-JX v7.0
!!EOP
!!------------------------------------------------------------------------------
!!BOC
!!
!! !LOCAL VARIABLES:
!!
!      INTEGER, PARAMETER    :: IPARAM1=2
!
!      !=================================================================
!      ! SRTN begins here!
!      !=================================================================
!
!      OARG3 = INARG1 + IOARG2
!      IOARG2 = OARG3/2.d0*IPARAM1
!
!      END SUBROUTINE SRTN
!!
!!EOC
!!------------------------------------------------------------------------------
!!          Harvard University Atmospheric Chemistry Modeling Group            !
!!------------------------------------------------------------------------------
!!BOP
!!
!! !IROUTINE: srtn
!!
!! !DESCRIPTION: Subroutine SRTN is a reference subroutine, with more than one
!!  line of description.
!!\\
!!\\
!! !INTERFACE:
!!
!      SUBROUTINE SRTN (INARG1,IOARG2,OARG3)
!!
!! !USES:
!!
!      USE THIS_MOD, ONLY : THIS_FN
!!
!! !INPUT PARAMETERS:
!!
!      REAL*8, INTENT(IN)          :: INARG1     ! Description of INARG1
!      REAL*8, INTENT(INOUT)       :: IOARG2     ! Description of IOARG2
!!
!! !OUTPUT VARIABLES:
!!
!      REAL*8, INTENT(OUT)         :: OARG3      ! Description of OARG3
!!
!! !REMARKS:
!! (1) A remark
!! 
!! !REVISION HISTORY: 
!!  28 Mar 2013 - S. D. Eastham - Copied from Fast-JX v7.0
!!EOP
!!------------------------------------------------------------------------------
!!BOC
!!
!! !LOCAL VARIABLES:
!!
!      INTEGER, PARAMETER    :: IPARAM1=2
!
!      !=================================================================
!      ! SRTN begins here!
!      !=================================================================
!
!      OARG3 = INARG1 + IOARG2
!      IOARG2 = OARG3/2.d0*IPARAM1
!
!      END SUBROUTINE SRTN
!!
!!EOC
!!------------------------------------------------------------------------------
!!          Harvard University Atmospheric Chemistry Modeling Group            !
!!------------------------------------------------------------------------------
!!BOP
!!
!! !IROUTINE: srtn
!!
!! !DESCRIPTION: Subroutine SRTN is a reference subroutine, with more than one
!!  line of description.
!!\\
!!\\
!! !INTERFACE:
!!
!      SUBROUTINE SRTN (INARG1,IOARG2,OARG3)
!!
!! !USES:
!!
!      USE THIS_MOD, ONLY : THIS_FN
!!
!! !INPUT PARAMETERS:
!!
!      REAL*8, INTENT(IN)          :: INARG1     ! Description of INARG1
!      REAL*8, INTENT(INOUT)       :: IOARG2     ! Description of IOARG2
!!
!! !OUTPUT VARIABLES:
!!
!      REAL*8, INTENT(OUT)         :: OARG3      ! Description of OARG3
!!
!! !REMARKS:
!! (1) A remark
!! 
!! !REVISION HISTORY: 
!!  28 Mar 2013 - S. D. Eastham - Copied from Fast-JX v7.0
!!EOP
!!------------------------------------------------------------------------------
!!BOC
!!
!! !LOCAL VARIABLES:
!!
!      INTEGER, PARAMETER    :: IPARAM1=2
!
!      !=================================================================
!      ! SRTN begins here!
!      !=================================================================
!
!      OARG3 = INARG1 + IOARG2
!      IOARG2 = OARG3/2.d0*IPARAM1
!
!      END SUBROUTINE SRTN
!!
!!EOC
      END MODULE FAST_JX_MOD
