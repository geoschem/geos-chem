!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: fast_jx_mod
!
! !DESCRIPTION: Module FAST\_JX\_MOD contains routines and variables for 
!  calculating photolysis rates using the Fast-JX scheme (Prather et al).
!  Current implementation is version 7.0a.
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

      PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
      PUBLIC  :: EXITC
      PUBLIC  :: PHOTO_JX
      PUBLIC  :: INIT_FJX
      PUBLIC  :: FAST_JX
      PUBLIC  :: FJXFUNC

!
! !PRIVATE MEMBER FUNCTIONS:
!
      PRIVATE :: SOLAR_JX
      PRIVATE :: OPMIE
      PRIVATE :: MIESCT 
      PRIVATE :: LEGND0
      PRIVATE :: BLKSLV 
      PRIVATE :: GEN_ID 
      PRIVATE :: JRATET
      PRIVATE :: X_INTERP 
      PRIVATE :: SPHERE2 
      PRIVATE :: EXTRAL
      PRIVATE :: RD_PROF_NC
      PRIVATE :: RD_XXX
      PRIVATE :: RD_AOD
      PRIVATE :: RD_MIE
      PRIVATE :: RD_JS_JX
      PRIVATE :: SET_AER
!
! !REVISION HISTORY:
!  27 Mar 2013 - S. D. Eastham - Initial version (based on original GEOS-Chem
!                                files and instructions from J. Mao)
!  28 Mar 2013 - S. D. Eastham - Upgraded to Fast-JX v7.0
!  20 Feb 2014 - M. Sulprizio  - Removed "define.h", this is now obsolete
!EOP
!------------------------------------------------------------------------------
!BOC
      CONTAINS      
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
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
      SUBROUTINE FAST_JX( WLAOD, am_I_Root, Input_Opt, State_Met,
     &                    State_Chm, RC )
!
! !USES:
!
      USE CHEMGRID_MOD,       ONLY : GET_CHEMGRID_LEVEL
      USE CMN_SIZE_MOD,       ONLY : IIPAR,   JJPAR, LLPAR
      USE CMN_SIZE_MOD,       ONLY : NDUST,   MAXIJ
      USE COMODE_MOD,         ONLY : JLOP, CSPEC, AIRDENS, ABSHUM
      USE COMODE_LOOP_MOD,    ONLY : CONSVAP,AVG,RSTARG
      USE ERROR_MOD,          ONLY : ERROR_STOP, ALLOC_ERR
      USE ERROR_MOD,          ONLY : DEBUG_MSG
      USE GIGC_ErrCode_Mod
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
      USE GIGC_State_Chm_Mod, ONLY : ChmState
      USE GIGC_State_Met_Mod, ONLY : MetState
      USE GRID_MOD,           ONLY : GET_YMID
      USE PRESSURE_MOD,       ONLY : GET_PEDGE
      USE TIME_MOD,           ONLY : GET_MONTH, GET_DAY, GET_DAY_OF_YEAR
      USE TIME_MOD,           ONLY : GET_TAU,   GET_YEAR
      USE TOMS_MOD,           ONLY : GET_OVERHEAD_O3
      USE TRACERID_MOD,       ONLY : IDTO3

      IMPLICIT NONE

!==============================================================================
! Uncomment the appropriate #define statement to denote which of the
! available cloud overlap options that you wish to use.

!! Linear overlap
!#define USE_LINEAR_OVERLAP 1

! Approximate random overlap (balance between accuracy & speed)
#define USE_APPROX_RANDOM_OVERLAP 1

!! Maximum random cloud overlap (most computationally intensive)
!#define USE_MAXIMUM_RANDOM_OVERLAP 1
!==============================================================================
!
! !INPUT PARAMETERS: 
!
      ! How was AOD calculated? (1: 550 nm, 0: 999 nm)
      INTEGER,        INTENT(IN)  :: WLAOD

      ! Is this the root CPU?
      LOGICAL,        INTENT(IN)  :: am_I_Root

      ! Input options
      TYPE(OptInput), INTENT(IN)  :: Input_Opt

      ! Meteorology State object
      TYPE(MetState), INTENT(IN)  :: State_Met

      ! Chemistry State object
      TYPE(ChmState), INTENT(IN)  :: State_Chm
!
! !OUTPUT PARAMETERS:
!
      INTEGER,        INTENT(OUT) :: RC
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
!  27 Mar 2013 - S.D. Eastham- Rolled in FAST_JX_MOD
!  21 Feb 2014 - M. Sulprizio- Replaced all met field arrays with State_Met
!                              derived type object
!  28 Feb 2014 - M. Sulprizio- Now obtain O3_CTM directly from State_Chm%Tracers
!  05 Mar 2014 - R. Yantosca - Cosmetic changes
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER, SAVE :: LASTMONTH = -1
      INTEGER       :: NLON, NLAT, DAY,  MONTH, DAY_OF_YR, L
      REAL*8        :: CSZA, PRES, SFCA, YLAT,  O3_TOMS
      REAL*8        :: O3_CTM(LLPAR+1)
      REAL*8        :: T_CTM(LLPAR+1), OPTD(LLPAR)
      REAL*8        :: OPTDUST(LLPAR,NDUST)
      REAL*8        :: OPTAER(LLPAR,NALLAER)

      ! Local variables for cloud overlap (hyl, phs)
      INTEGER       :: NUMB, KK, I
      INTEGER       :: INDIC(LLPAR+1)
      INTEGER       :: INDGEN(LLPAR+1)! = (/ (i,i=1,LLPAR+1) /)
      INTEGER       :: KBOT(LLPAR)
      INTEGER       :: KTOP(LLPAR)
      INTEGER       :: INDICATOR(LLPAR+2)
      REAL*8        :: FMAX(LLPAR)     ! maximum cloud fraction 
                                       !  in a block, size can be to 
                                       !  FIX(LLPAR)+1
      REAL*8        :: CLDF1D(LLPAR)
      REAL*8        :: ODNEW(LLPAR)
      REAL*8        :: P_CTM(LLPAR+2)

      LOGICAL       :: AOD999

      LOGICAL, SAVE :: FIRST = .true.

      ! Local variables for quantities from Input_Opt
      LOGICAL       :: LPRT

      !=================================================================
      ! FAST_JX begins here!
      !=================================================================

      ! Assume success
      RC = GIGC_SUCCESS

      ! Copy fields from INPUT_OPT
      LPRT = ( am_I_Root .and. Input_Opt%LPRT )

#if defined( USE_MAXIMUM_RANDOM_OVERLAP )
      IF ( FIRST ) THEN

         DO i = 1,LLPAR+1
            INDGEN(i) = i       !(/(i,i=1,LLPAR+1)/)
         ENDDO
         FIRST = .FALSE.

      ENDIF
#endif

      ! Get day of year (0-365 or 0-366)
      DAY_OF_YR = GET_DAY_OF_YEAR()

      ! Get current month
      MONTH     = GET_MONTH()

      ! Get day of month
      DAY       = GET_DAY()

      ! Was AOD calculated at 999 nm or reference?
      AOD999    = ( WLAOD == 0 )

      !=================================================================
      ! For each (NLON,NLAT) location, call subroutine PHOTO_JX (in a 
      ! parallel loop to compute J-values for the entire column.  
      ! J-values will be stored in the common-block variable ZPJ, and 
      ! will be later accessed via function FJXFUNC. 
      !=================================================================
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( NLAT,    NLON,   YLAT,      CSZA,   P_CTM   )
!$OMP+PRIVATE( T_CTM,   SFCA,   O3_TOMS,   O3_CTM, OPTAER  )
!$OMP+PRIVATE( OPTDUST, OPTD,   CLDF1D                     )
#if defined( USE_MAXIMUM_RANDOM_OVERLAP )
!$OMP+PRIVATE( FMAX,    KK,     NUMB,      L       KBOT    )
!$OMP+PRIVATE( KTOP     ODNEW,  INDICATOR, INDIC           )
!$OMP+PRIVATE( O3_TOMS, P_CTM,  JLOOP                      ) 
#endif
!$OMP+SCHEDULE( DYNAMIC )

      ! Loop over latitudes
      DO NLAT = 1, JJPAR

         ! Loop over longitudes
         DO NLON = 1, IIPAR

            ! Grid box latitude [degrees]
            YLAT            = GET_YMID( NLON, NLAT, 1 )

            ! Cosine of solar zenith angle [unitless] at (NLON,NLAT)
            CSZA            = State_Met%SUNCOSmid(NLON,NLAT)

#if defined( EXTERNAL_GRID ) || defined( EXTERNAL_FORCING )
            ! Define the P array here, using GCM pressures
            DO L = 1, LLPAR+1
               P_CTM(L)     = State_Met%PEDGE( NLON, NLAT, L )
            ENDDO
#else
            ! Define the P array here, using Ap and Bp from GEOS-Chem
            DO L = 1, LLPAR+1
               P_CTM(L)     = GET_PEDGE( NLON, NLAT, L )
            ENDDO
#endif

            ! Top edge of P_CTM is top of atmosphere (bmy, 2/13/07)
            P_CTM(LLPAR+2)  = 0d0

            ! Temperature profile [K] at (NLON,NLAT)
            T_CTM(1:LLPAR)  = State_Met%T(NLON,NLAT,1:LLPAR)

            ! Top of atmosphere
            T_CTM(LLPAR+1)  = T_CTM(LLPAR)

            ! Surface albedo [unitless] at (NLON,NLAT)
            SFCA            = State_Met%UVALBEDO(NLON,NLAT)

            ! Overhead ozone column [DU] at (NLON, NLAT)
            O3_TOMS         = GET_OVERHEAD_O3( NLON, NLAT )
            
            ! CTM ozone densities (molec/cm3) at (NLON, NLAT)
            O3_CTM          = 0d0
            L               = GET_CHEMGRID_LEVEL(NLON,NLAT,State_Met)
            O3_CTM(1:L)     = State_Chm%Tracers(NLON,NLAT,1:L,IDTO3)

            ! Aerosol OD profile [unitless] at (NLON,NLAT)
            OPTAER          = ODAER(NLON,NLAT,:,:)

            ! Mineral dust OD profile [unitless] at (NLON,NLAT)
            OPTDUST         = ODMDUST(NLON,NLAT,:,:)

            ! Cloud OD profile [unitless] at (NLON,NLAT)
            OPTD            = State_Met%OPTD(1:LLPAR,NLON,NLAT)

            !-----------------------------------------------------------
            !### If you want to exclude aerosol OD, mineral dust OD,
            !### or cloud OD, then uncomment the following lines:
            !OPTAER  = 0d0
            !OPTDUST = 0d0
            !OPTD    = 0d0
            !-----------------------------------------------------------

#if defined( USE_LINEAR_OVERLAP )
            !===========================================================
            ! %%%% CLOUD OVERLAP: LINEAR ASSUMPTION %%%%
            !
            ! Directly use OPTDEPTH = TAUCLD * CLDTOT
            !===========================================================

#if defined( GEOS_5 ) || defined( MERRA ) || defined( GEOS_FP )

            ! Column cloud fraction (not less than zero)
            CLDF1D = State_Met%CLDF(1:LLPAR,NLON,NLAT)
            WHERE ( CLDF1D < 0d0 ) CLDF1D = 0d0
               
            ! NOTE: For GEOS-5, GEOS-FP, and MERRA met fields, the optical
            ! depth is the in-cloud optical depth.  At this point it has 
            ! NOT been multiplied by cloud fraction yet.  Therefore, we can
            ! just apply the ! we can just apply the linear overlap formula 
            ! as written above (i.e. multiply by cloud fraction).
            OPTD = OPTD * CLDF1D

#endif

            ! Call FAST-JX routines to compute J-values
            CALL PHOTO_JX( CSZA,      SFCA,      P_CTM,  T_CTM,
     &                     O3_CTM,    O3_TOMS,   AOD999, OPTAER,
     &                     OPTDUST,   OPTD,      NLON,   NLAT,
     &                     YLAT,      DAY_OF_YR, MONTH,  DAY,
     &                     am_I_Root, Input_Opt  )


#elif defined( USE_APPROX_RANDOM_OVERLAP )
            !===========================================================
            ! %%%% CLOUD OVERLAP: APPROX RANDOM OVERLAP ASSUMPTION %%%%
            !
            ! Use OPTDEPTH = TAUCLD * CLDTOT**1.5
            !===========================================================

            ! Column cloud fraction (not less than zero)
            CLDF1D = State_Met%CLDF(1:LLPAR,NLON,NLAT)
            WHERE ( CLDF1D < 0d0 ) CLDF1D = 0d0
               
#if   defined( GEOS_5 ) || defined( MERRA ) || defined( GEOS_FP )

            ! NOTE: For GEOS-5, GEOS-FP, and MERRA met fields, the optical
            ! depth is the in-cloud optical depth.  At this point it has 
            ! NOT been multiplied by cloud fraction yet.  Therefore, we can
            ! just apply the approximate random overlap formula as written 
            ! above (i.e. multiply by cloud fraction^1.5).
            OPTD = OPTD * ( CLDF1D )**1.5d0
            
#else
            ! Otherwise, OPTD is the grid-box optical depth and has 
            ! already been multiplied by  the cloud fraction.  Therefore 
            ! we only need to multiply by the square root of the cloud 
            ! fraction here for the approximate random overlap option. 
            OPTD = OPTD * SQRT( CLDF1D )

#endif

            ! Call FAST-JX routines to compute J-values
            CALL PHOTO_JX( CSZA,      SFCA,      P_CTM,  T_CTM,
     &                     O3_CTM,    O3_TOMS,   AOD999, OPTAER,
     &                     OPTDUST,   OPTD,      NLON,   NLAT,
     &                     YLAT,      DAY_OF_YR, MONTH,  DAY,
     &                     am_I_Root, Input_Opt  )

#elif defined( USE_MAXIMUM_RANDOM_OVERLAP )
            !===========================================================
            ! %%%% CLOUD OVERLAP: MAXIMUM RANDOM OVERLAP %%%%
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

               ! Sort this out later
               CALL ERROR_STOP('MMRAN_16 not yet FJX compatible.',
     &                         'fast_jx_mod.F')
!               ! Initialize
!               FMAX(:)   = 0d0  ! max cloud fraction in each cloud block
!               ODNEW(:)  = 0d0  ! in-cloud optical depth
!               CLDF1D    = State_Met%CLDF(1:LLPAR,NLON,NLAT)
!               INDICATOR = 0
!
!               ! set small negative CLDF or OPTD to zero. 
!               ! Set indicator vector.
!               WHERE ( CLDF1D <= 0d0 ) 
!                  CLDF1D               = 0d0
!                  OPTD                 = 0D0
!               ELSEWHERE
!                  INDICATOR(2:LLPAR+1) = 1
!               ENDWHERE
!
!               ! Prevent negative opt depth
!               WHERE ( OPTD < 0D0 ) OPTD   = 0D0
!
!               !--------------------------------------------------------
!               ! Generate cloud blocks & get their Bottom and Top levels
!               !--------------------------------------------------------
!               INDICATOR = CSHIFT(INDICATOR, 1) - INDICATOR
!               INDIC     = INDICATOR(1:LLPAR+1)
!
!               ! Number of cloud block
!               NUMB      = COUNT( INDIC == 1 ) 
!               
!               ! Bottom layer of each block
!               KBOT(1:NUMB) = PACK(INDGEN, (INDIC == 1 ) ) 
!
!               ! Top layer of each block
!               KTOP(1:NUMB) = PACK(INDGEN, (INDIC == -1) ) - 1 
!             
!               !--------------------------------------------------------
!               ! For each cloud block, get Max Cloud Fractions, and 
!               ! in-cloud optical depth vertical distribution.
!               !--------------------------------------------------------
!               DO KK = 1, NUMB
!
!                  ! Max cloud fraction
!                  FMAX(KK) = MAXVAL( CLDF1D(KBOT(KK):KTOP(KK)) )
!
!#if   defined( GEOS_5 ) || defined( MERRA ) || defined( GEOS_57 )
!
!                  ! NOTE: for the reprocessed GEOS-5 met fields (i.e. with
!                  ! optical depth & cloud fractions regridded with RegridTau)
!                  ! OPTD is the in-cloud optical depth.  At this point it has
!                  ! NOT been multiplied by cloud fraction yet.  Therefore,
!                  ! we can just set ODNEW = OPTD. (bmy, hyl, 10/24/08)
!
!                  ! ODNEW is adjusted in-cloud OD vertical distrib.
!                  ODNEW(KBOT(KK):KTOP(KK)) = OPTD(KBOT(KK):KTOP(KK))
!
!#else
!
!                  ! Otherwise, OPTD is the grid-box optical depth.  
!                  ! Therefore, we must divide out by the cloud fraction
!                  ! and thus set ODNEW = OPTD / FMAX. (bmy, hyl, 10/24/08)
!
!                  ! ODNEW is adjusted in-cloud OD vertical distrib.
!                  ODNEW(KBOT(KK):KTOP(KK)) = OPTD(KBOT(KK):KTOP(KK)) / 
!     &                                       FMAX(KK)
!
!#endif
!               ENDDO
!            
!               !--------------------------------------------------------
!               ! Apply Max RANdom if 1-6 clouds blocks, else use linear 
!               !--------------------------------------------------------
!               SELECT CASE( NUMB ) 
!           
!                  CASE( 0,7: )
!                     CALL PHOTOJ( NLON,  NLAT,     YLAT,    DAY_OF_YR, 
!     &                            MONTH, DAY,      CSZA,    TEMP,  
!     &                            SFCA,  OPTD,     OPTDUST, OPTAER,
!     &                            O3COL, am_I_Root                    )
!
!                  CASE( 1:6 ) 
!                     CALL MMRAN_16( NUMB,  NLON,  NLAT,      YLAT,   
!     &                              DAY,   MONTH, DAY_OF_YR, CSZA,    
!     &                              TEMP,  SFCA,  OPTDUST,   OPTAER, 
!     &                              LLPAR, FMAX,  ODNEW,     KBOT,   
!     &                              KTOP,  O3COL, am_I_Root )
!
!               END SELECT
!           ENDIF
#endif
         ENDDO
      ENDDO
!$OMP END PARALLEL DO

      FIRST=.FALSE.

      !-----------------------------------------------------------
      ! END OF SUBROUTINE FAST-J
      !-----------------------------------------------------------
      END SUBROUTINE FAST_JX
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: blkslv
!
! !DESCRIPTION: Subroutine BLKSLV solves the block tri-diagonal system
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
! The block tri-diagonal system:
!       A(I)*X(I-1) + B(I)*X(I) + C(I)*X(I+1) = H(I)
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
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
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
!                  GEOS-Chem Global Chemical Transport Model                  !
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
      SUBROUTINE JRATET(PPJ,TTJ,FFF,VALJL,LCTM,LCHEM,NJXU)
!
! !USES:
!
!
! !INPUT PARAMETERS:
!
      integer, intent(in)   :: LCTM,LCHEM,NJXU
      real*8, intent(in)    ::  PPJ(JXL1_+1),TTJ(JXL1_+1)
      real*8, intent(inout) ::  FFF(W_,LCTM)
!
! !OUTPUT VARIABLES:
!
      real*8, intent(out), dimension(LCTM,NJXU) ::  VALJL
!
! !REMARKS:
! 
! !REVISION HISTORY: 
!  28 Mar 2013 - S. D. Eastham - Copied from Fast-JX v7.0
!  02 Apr 2013 - S. D. Eastham - Now only assign J-rates on
!                                levels within the chemistry grid.
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
      do L = 1,LCTM
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

        do J = 1,NJXU
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

        do J = 4,NJXU
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

        do J=1,NJXU
          VALJL(L,J) = VALJ(J)
        enddo

      enddo

      ! Zero non-chemistry layers
      if (LCHEM.lt.LCTM) then
         do L=(LCTM+1),LCHEM
            do J=1,NJXU
               VALJL(L,J) = 0.d0
            enddo
         enddo
      endif

      END SUBROUTINE JRATET
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
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
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: sphere2
!
! !DESCRIPTION: Subroutine SPHERE2 is an AMF2.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE SPHERE2 (U0,ZHL,AMF2,L1U,LJX1U)
!
! !USES:
!
!
! !INPUT PARAMETERS:
!
      INTEGER, INTENT(IN) ::   L1U, LJX1U
      REAL*8, INTENT(IN)  ::   U0,ZHL(L1U+1)
!
! !OUTPUT VARIABLES:
!
      REAL*8, INTENT(OUT) ::   AMF2(2*LJX1U+1,2*LJX1U+1)
!
! !REMARKS:
! Quoting from the original:
!  New v6.2: does AirMassFactors for mid-layer, needed for SZA ~ 90
!  This new AMF2 does each of the half-layers of the CTM separately,
!     whereas the original, based on the pratmo code did the whole layers
!     and thus calculated the ray-path to the CTM layre edges, NOT the middle.
!  Since fast-JX is meant to calculate the intensity at the mid-layer, the
!     solar beam at low sun (interpolated between layer edges) was incorrect.
!  This new model does make some approximations of the geometry of the layers:
!     the CTM layer is split evenly in mass (good) and in height (approx).
!                                                                             .
!  Calculation of spherical geometry; derive tangent heights, slant path
!  lengths and air mass factor for each layer. Not called when
!  SZA > 98 degrees.  Beyond 90 degrees, include treatment of emergent
!  beam (where tangent height is below altitude J-value desired at).
!                                                                             .
!  ---------------------------------------------------------------------
!  Inputs:
!     U0      cos(solar zenith angle)
!     RAD  radius of Earth mean sea level (cm)
!     ZHL(L)  height (cm) of the bottom edge of CTM level L
!     ZZHT    scale height (cm) used above top of CTM (ZHL(L_+1))
!     L1U     dimension of CTM = levels +1 (L+1 = above-CTM level)
!  Outputs:
!     AMF2(I,J) = air mass factor for CTM level I for sunlight reaching J
!         ( these are calculated for both layer middle and layer edge)
!  ---------------------------------------------------------------------
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
c      ZHL(L1U+1) = ZHL(L1U) + ZZHT

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
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: extral
!
! !DESCRIPTION: Subroutine EXTRAL adds sub-layers to thick cloud/aerosol layers
!  using log-spacing for sub-layers of increasing thickness ATAU.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE EXTRAL (DTAUX,L1X,L2X,NX,JXTRA)
!
! !USES:
!
!
! !INPUT PARAMETERS:
!
      INTEGER, INTENT(IN) ::  L1X,L2X         !index of cloud/aerosol
      integer, intent(in) ::  NX              !Mie scattering array size
      real*8,  intent(in) ::  DTAUX(L1X)      !cloud+3aerosol OD in each layer
!
! !OUTPUT VARIABLES:
!
      integer, intent(out)::  JXTRA(L2X+1)    !number of sub-layers to be added
!
! !REMARKS:
!     DTAUX(L=1:L1X) = Optical Depth in layer L (generally 600 nm OD)
!        This can be just cloud or cloud+aerosol, it is used only to set
!        the number in levels to insert in each layer L
!        Set for log-spacing of tau levels, increasing top-down.
!                                                                             .
!     N.B. the TTAU, etc calculated here are NOT used elsewhere
!                                                                             .
!   The log-spacing parameters have been tested for convergence and chosen
!     to be within 0.5% for ranges OD=1-500, rflect=0-100%, mu0=0.1-1.0
!     use of ATAU = 1.18 and min = 0.01, gives at most +135 pts for OD=100
!     ATAU = 1.12 now recommended for more -accurate heating rates (not J's)
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
c---These are set in CMN_FJX_MOD, and have been tested/optimized

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

      END SUBROUTINE EXTRAL
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: exitc
!
! !DESCRIPTION: Subroutine EXITC forces an error in GEOS-Chem and quits.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE EXITC (T_EXIT)
!
! !USES:
!
      USE ERROR_MOD, ONLY : ERROR_STOP
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

      CALL ERROR_STOP( T_EXIT, 'fast_jx_mod.F' )

      END SUBROUTINE EXITC
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
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
      SUBROUTINE INIT_FJX(am_I_Root)
!
! !USES:
!
      USE inquireMod,   ONLY : findFreeLUN
      USE CMN_DIAG_MOD, ONLY : ND64
!
! !INPUT PARAMETERS:
!
      LOGICAL, INTENT(IN) :: am_I_Root   ! Is this the root CPU?
!
! !OUTPUT VARIABLES:
!
!
! !REMARKS:
! 
! !REVISION HISTORY: 
!  28 Mar 2013 - S. D. Eastham - Copied from Fast-JX v7.0
!  05 Mar 2014 - R. Yantosca - Use correct filenames for use w/ ESMF
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER                     :: JXUNIT,J
      CHARACTER(LEN=255)          :: FILENAME
      INTEGER                     :: NJXX
      CHARACTER*6,DIMENSION(JVN_)  :: TITLEJXX
      REAL*8                      :: ND64MULT

      !=================================================================
      ! INIT_FJX begins here!
      !=================================================================

      IF ( am_I_Root) THEN
         write(6,*) ' Initializing Fast-JX v7.0 standalone CTM code.'

         if (W_.ne.8 .and. W_.ne.12 .and. W_.ne.18) then
           call EXITC(' INIT_FJX: invalid no. wavelengths')
         endif
      ENDIF

      ! Set shared variable according to input argument
      

      ! Get a free LUN
      JXUNIT = findFreeLUN()

      ! Read in fast-J X-sections (spectral data)
      FILENAME = 'FJX_spec.dat'

#if defined( ESMF_ )
      ! For the ESMF interface to GEOS-5 GCM, we need to append *___.rc
      FILENAME = TRIM( FILENAME ) // '.rc'
#endif

      ! Echo info
      IF ( am_I_Root ) THEN
         WRITE( 6, 100 ) TRIM( FILENAME ), JXUNIT
 100     FORMAT( 'GEOSCHEM::INIT_FJX: Reading ', a ', on unit ', i4 )
      ENDIF

      call RD_XXX(JXUNIT,TRIM( FILENAME ), am_I_Root)

      ! If we are using ND64, calculate exposure factors
      IF ( ND64 > 0 ) THEN
         UVXFACTOR(:) = 0d0
         ND64MULT = UVXPLANCK*UVXCCONST*1.0d13
         DO  J=1,W_
            UVXFACTOR(J) = ND64MULT/WL(J)
         ENDDO
      ENDIF

      ! Read in 5-wavelength scattering data
      FILENAME = 'jv_spec_mie.dat'

#if defined( ESMF_ )
      ! For the ESMF interface to GEOS-5 GCM, we need to append *___.rc
      FILENAME = TRIM( FILENAME ) // '.rc'
#endif

      ! Echo info
      IF ( am_I_Root ) THEN
         WRITE( 6, 100 ) TRIM( FILENAME ), JXUNIT
      ENDIF

      call RD_MIE(JXUNIT,TRIM( FILENAME ), am_I_Root)

      ! Read in AOD data
      FILENAME = 'jv_spec_aod.dat'

#if defined( ESMF_ )
      ! For the ESMF interface to GEOS-5 GCM, we need to append *___.rc
      FILENAME = TRIM( FILENAME ) // '.rc'
#endif

      ! Echo info
      IF ( am_I_Root ) THEN
         WRITE( 6, 100 ) TRIM( FILENAME ), JXUNIT
      ENDIF

      call RD_AOD(JXUNIT,TRIM( FILENAME ), am_I_Root)

      ! Set up MIEDX array to interpret between GC and FJX aerosol indexing
      call SET_AER(am_I_Root)

      ! Read in T & O3 climatology used to fill e.g. upper layers or if O3 not calc.
      call RD_PROF_NC( am_I_Root )

      NJXX = NJX
      do J = 1,NJX
        TITLEJXX(J) = TITLEJX(J)
      enddo

      ! Read in photolysis rates used in chemistry code and mapping onto FJX J's
      ! CTM call:  read in J-values names and link to fast-JX names
      FILENAME = 'FJX_j2j.dat'

#if defined( ESMF_ )
      ! For the ESMF interface to GEOS-5 GCM, we need to append *___.rc
      FILENAME = TRIM( FILENAME ) // '.rc'
#endif

      ! Echo info
      IF ( am_I_Root ) THEN
         WRITE( 6, 100 ) TRIM( FILENAME ), JXUNIT
      ENDIF

      call RD_JS_JX(JXUNIT,TRIM( FILENAME ), TITLEJXX,NJXX,am_I_Root)

      ! Build cross-referencing indices
      call JV_INDEX(am_I_Root)

      END SUBROUTINE INIT_FJX
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: rd_xxx
!
! !DESCRIPTION: Subroutine RD\_XXX reads in wavelength bins, solar fluxes,
!  Rayleigh and temperature-dependent cross-sections.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE RD_XXX (NUN,NAMFIL,am_I_Root)
!
! !USES:
!
!
! !INPUT PARAMETERS:
!
      INTEGER, INTENT(IN) :: NUN
      CHARACTER(*), INTENT(IN) ::  NAMFIL
      LOGICAL,      INTENT(IN) :: am_I_Root
!
! !OUTPUT VARIABLES:
!
!
! !REMARKS:
!    NEW v-6.8  now allow 1 to 3 sets of X-sects for T or P
!           LQQ = 1, 2, or 3 to determine interpolation with T or P
!           IF the temperatures TQQQ are <0, then use as pressure interp (hPa)
!           NB - the temperatures and pressures must be increasing
!    NEW v-6.4  changed to collapse wavelengths & x-sections to Trop-only:
!           WX_ = 18 should match the JX_spec.dat wavelengths
!           W_ = 12 (Trop-only) or 18 (std) is set in (CMN_FJX.F).
!       if W_=12 then drop strat wavels, and drop x-sects (e.g. N2O, ...)
!           W_ = 8, reverts to quick fix:  fast-J (12-18) plus bin (5) scaled
!                                                                             .
!   --------------------------------------------------------------------
!     NAMFIL   Name of spectral data file (FJX_spec.dat) >> j2 for fast-J2
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
! 
! !REVISION HISTORY: 
!  28 Mar 2013 - S. D. Eastham - Copied from Fast-JX v7.0
!  09 May 2014 - S. D. Eastham - Fixed bug when using <18 wavelengths
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
        if (am_I_Root) then
           write(6,200) J,TITLEJX(J),SQQ(J),LQQ(J),(TQQ(I,J),I=1,LQQ(J))
        endif
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
        if (am_I_Root) then
           write(6,'(a)')
     &      ' >>>TROP-ONLY reduce wavelengths to 12, drop strat X-sects'
        endif
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
          do LQ=1,LQQ(J) 
           QQQ(IW,LQ,J) = QQQ(IW+4,LQ,J)
          enddo
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
          do LQ=1,LQQ(J) 
           QQQ(IW,LQ,J) = QQQ(IW+6,LQ,J)
          enddo
         enddo
        enddo
c---TROP-QUICK  (must scale solar flux for W=5)
       elseif (W_ .eq. 8) then
        if (am_I_Root) then
        write(6,'(a)')
     &   ' >>>TROP-QUICK reduce wavelengths to 8, drop strat X-sects'
        endif
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
          do LQ=1,LQQ(J) 
           QQQ(IW,LQ,J) = QQQ(IW+4,LQ,J)
          enddo
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
          do LQ=1,LQQ(J) 
           QQQ(IW,LQ,J) = QQQ(IW+10,LQ,J)
          enddo
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
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: rd_mie
!
! !DESCRIPTION: Subroutine RD\_MIE retrieves aerosol scattering data for FJX.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE RD_MIE(NUN,NAMFIL,am_I_Root)
!
! !USES:
!
      USE FILE_MOD,   ONLY : IOERROR
!
! !INPUT PARAMETERS:
!
      INTEGER,      INTENT(IN) :: NUN
      CHARACTER(*), INTENT(IN) :: NAMFIL
      LOGICAL,      INTENT(IN) :: am_I_Root
!
! !OUTPUT VARIABLES:
!
!
! !REMARKS:
!   --------------------------------------------------------------------
!     NAMFIL   Name of scattering data file (e.g., FJX_scat.dat)
!     NUN      Channel number for reading data file
!     NAA      Number of categories for scattering phase functions
!     QAA      Aerosol scattering phase functions
!     WAA      5 Wavelengths for the supplied phase functions
!     PAA      Phase function: first 8 terms of expansion
!     RAA      Effective radius associated with aerosol type
!     SAA      Single scattering albedo
!   --------------------------------------------------------------------
!
! !REVISION HISTORY: 
!  28 Mar 2013 - S. D. Eastham - Adapted from GEOS-Chem v9-1-3
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER  I, J, K, NK
      CHARACTER*78 TITLE0

      !=================================================================
      ! RD_MIE begins here!
      !=================================================================

      open (NUN,FILE=NAMFIL,status='old',form='formatted')

      ! Read header lines
      READ( NUN,'(A)' ) TITLE0
      IF  ( am_I_Root ) WRITE( 6, '(1X,A)' ) TITLE0
      READ( NUN,'(A)' ) TITLE0

!---Read aerosol phase functions:
      read(NUN,'(A10,I5,/)') TITLE0,NAA
      NK=5
      do j=1,NAA
        read(NUN,110) TITLEAA(j)
        do k=1,NK
          read(NUN,*) WAA(k,j),QAA(k,j),RAA(k,j),SAA(k,j),
     &                                             (PAA(i,k,j),i=1,8)
        enddo
      enddo

      close(NUN)

      IF ( am_I_Root ) THEN 
         write(6,'(a,9f8.1)') ' Aerosol optical: r-eff/rho/Q(@wavel):'
     &             ,(WAA(K,1),K=1,5)
         do J=1,NAA
            write(6,'(1x,A8,I2,A,9F8.1)')
     &           TITLEAA(J),J,'  wavel=',(WAA(K,J),K=1,NK)
            write(6,'(9x,I2,A,9F8.4)') J,'  Qext =',(QAA(K,J),K=1,NK)
         enddo
      ENDIF

  110 format(3x,a20)

      END SUBROUTINE RD_MIE
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: rd_aod
!
! !DESCRIPTION: Subroutine RD\_AOD reads aerosol phase functions that are
!  used to scale diagnostic output to an arbitrary wavelengh.  This
!  facilitates comparing with satellite observations.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE RD_AOD( NJ1, NAMFIL, am_I_Root )
!
! !USES:
!
      USE FILE_MOD,   ONLY : IOERROR

      IMPLICIT NONE

!
! !INPUT PARAMETERS:
!     
      INTEGER,          INTENT(IN) :: NJ1         ! Unit # of file to open
      CHARACTER(LEN=*), INTENT(IN) :: NAMFIL      ! Name of file to open
      LOGICAL,          INTENT(IN) :: am_I_Root   ! Is this the root CPU?
!
! !INPUT/OUTPUT PARAMETERS:
!
!
! !REMARKS:
!  The jv_spec_aod.dat file contains the optical properties for aerosols 
!  at a single wavelength to be used in the online calculation of the aerosol 
!  optical depth diagnostics.  The default properties are provided at 550 nm. 
!  These properties have been calculated using the same size and optical 
!  properties as the jv_spec.dat file used for the FAST-J photolysis 
!  calculations.  The user can exchange this set of properties with those at 
!  another wavelength.  We recommend that the wavelength used be included 
!  in the first line of the header for traceability (this line is output to 
!  the GEOS-Chem log file during run time). A complete set of optical 
!  properties from 250-2000 nm for aerosols is available at:  
!  ftp://ftp.as.harvard.edu/geos-chem/data/aerosol_optics/hi_spectral_res
!                                                                             .
!     -- Colette L. Heald, 05/10/10)
!
!  Important variables:
!                                                                             .
!     NAMFIL       Name of spectral data file (jv_spec_aod.dat)
!     NJ1          Channel number for reading data file
!     NAA2         Number of categories for scattering phase functions
!     QAA_AOD      Aerosol scattering phase functions
!     WAA_AOD      Wavelengths for the NK supplied phase functions
!     PAA_AOD      Phase function: first 8 terms of expansion
!     RAA_AOD      Effective radius associated with aerosol type
!     SAA_AOD      Single scattering albedo
!
! !REVISION HISTORY: 
!  10 May 2010 - C. Heald    - Initial version
!  06 Aug 2010 - C. Carouge  - Add an error check when opening the file
!  01 Aug 2012 - R. Yantosca - Now restore NJ1 to INTENT(IN) status
!  31 Mar 2013 - S. D. Eastham - Added to Fast-JX v7.0 implementation
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES
!
      INTEGER :: I, J, K, NAA2
      INTEGER :: IOS
      CHARACTER*78 TITLE0
      CHARACTER*20 TITLEA(NAA)   ! TITLEAA: Title for scattering data

      !================================================================
      ! RD_AOD begins here!
      !================================================================

      ! open file
      OPEN( NJ1, FILE=TRIM( NAMFIL ), STATUS='OLD', IOSTAT=IOS )

      ! Error check
      IF ( IOS /= 0 ) THEN
         IF ( am_I_Root ) THEN
            WRITE(6,100) trim(NAMFIL)
 100        FORMAT('Error opening filename=', a )
            CALL FLUSH(6)
         ENDIF
         CALL IOERROR( IOS, NJ1, 'RD_AOD:1')
      ENDIF

      
      ! Read header lines
      READ( NJ1,'(A)' ) TITLE0
      IF  ( am_I_Root ) WRITE( 6, '(1X,A)' ) TITLE0
      READ( NJ1,'(A)' ) TITLE0

      ! Read aerosol phase functions (one wavelength only):
      READ( NJ1,'(A10,I5,/)' ) TITLE0,NAA2
      DO J = 1, NAA
         READ(NJ1,110) TITLEA(j)
 110     FORMAT( 3x, a20 )
         IF ( am_I_Root ) WRITE(6,*) TITLEA(j)
         READ(NJ1,*) WAA_AOD(j),QAA_AOD(j),RAA_AOD(j),SAA_AOD(j),
     &               (PAA_AOD(i,j),i=1,8)
      ENDDO

      ! Echo info to stdout
      IF ( am_I_Root ) THEN
         WRITE( 6, '(a)' ) 'Aerosol Qext for AOD calculations'
         DO J=1,NAA
            WRITE( 6, * ) TITLEA(J),J,'  Qext =',(QAA_AOD(J))
         ENDDO
      ENDIF

      ! Close file
      CLOSE( NJ1 )

      END SUBROUTINE RD_AOD
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: rd_js_jx
!
! !DESCRIPTION: Subroutine RD\_JS\_JX reads in 'FJX\_j2j.dat', which defines 
!  the mapping of Fast-JX J's (TITLEJX(1:NJX)) onto the CTM reactions.
!  Reaction number JJ, named T\_REACT, uses Fast-JX's T\_FJX (including scaling
!  factor F\_FJX).
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE RD_JS_JX(NUNIT,NAMFIL,TITLEJX,NJXX,am_I_Root)
!
! !USES:
!
!
! !INPUT PARAMETERS:
!
      INTEGER, INTENT(IN)                    ::  NUNIT, NJXX
      CHARACTER(*), INTENT(IN)               ::  NAMFIL
      CHARACTER*6, INTENT(IN),DIMENSION(NJXX):: TITLEJX
      LOGICAL, INTENT(IN)                    :: am_I_Root
!
! !OUTPUT VARIABLES:
!
!
! !REMARKS:
! 
! !REVISION HISTORY: 
!  28 Mar 2013 - S. D. Eastham - Copied from Fast-JX v7.0
!  31 Mar 2013 - S. D. Eastham - Modified to track number of branches
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

! Read the FJX_j2j.dat file to map model specific J's onto fast-JX J's
! The chemistry code title describes fully the reaction (a50)
! Blank (unfilled) chemistry J's are unmapped
! The number NRATJ is the last JJ readin that is .le. JVN
!   include fractional quantum yield for the fast-JX J's

      JLABEL(:) = '------'
      JMAP(:)   = '------'
      JFACTA(:) = 0.d0

      open (NUNIT,file=NAMFIL,status='old',form='formatted')

       read (NUNIT,'(a)') CLINE
      IF ( am_I_Root ) THEN 
         write(6,'(a)') CLINE
      ENDIF
      do J = 1,JVN_
       read (NUNIT,'(i4,1x,a50,4x,f5.3,2x,a6)') JJ,T_REACT,F_FJX,T_FJX
       if (JJ.gt.JVN_) goto 20
        JLABEL(JJ) = T_REACT
        JFACTA(JJ) = F_FJX
        JMAP(JJ) = T_FJX
        NRATJ = JJ
        ! SDE 03/31/13: Check number of branches
        ! Note that the order of the branches in
        ! globchem.dat must match the order in 
        ! FJX_j2j.dat
        READ (T_REACT(1:10),"(a10)") RNAMES(JJ)
        RNAMES(JJ) = TRIM(RNAMES(JJ))
        BRANCH(JJ) = 1
        DO K=1,(JJ-1)
           IF (RNAMES(JJ) == RNAMES(K)) THEN
              BRANCH(JJ) = BRANCH(K) + 1
           ENDIF
        ENDDO
      enddo

 20   close(NUNIT)

c---Zero / Set index arrays that map Jvalue(j) onto rates
      do K = 1,NRATJ
         JIND(K) = 0
       do J = 1,NJXX
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
            IF ( am_I_Root ) THEN 
               if (J.eq.0) then
         write(6,'(i5,a50,f6.3,a,1x,a6)') K,JLABEL(K),JFACTA(K),
     &         ' no mapping onto fast-JX',JMAP(K)
               else
         write(6,'(i5,a50,f6.3,a,i4,1x,a6)') K,JLABEL(K),JFACTA(K),
     &         ' mapped to FJX:',J,TITLEJX(J)
               endif
            ENDIF
         endif
      enddo

      END SUBROUTINE RD_JS_JX
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: solar_jx
!
! !DESCRIPTION: Subroutine SOLAR\_JX handles solar zenith angles.
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
!  ---------------------------------------------------------------------
!     NDAY   = integer day of the year (used for solar lat and declin)
!     SZA = solar zenith angle in degrees
!     COSSZA = U0 = cos(SZA)
!     SOLFX = Solar function
!  ---------------------------------------------------------------------
! 
! !REVISION HISTORY: 
!  28 Mar 2013 - S. D. Eastham - Adapted from Fast-JX v7.0
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      REAL*8  PI, PI180

      !=================================================================
      ! SOLAR_JX begins here!
      !=================================================================

      PI     = 3.141592653589793d0
      PI180  = PI/180.d0
      SZA    = acos(COSSZA)/PI180

      ! Offset used for GEOS-Chem slightly different
      !SOLFX  = 1.d0-(0.034d0*cos(dble(NDAY-186)*2.d0*PI/365.d0))
      SOLFX  = 1.d0-(0.034d0*cos(dble(NDAY-172)*2.d0*PI/365.d0))

      END SUBROUTINE SOLAR_JX
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: photo_jx
!
! !DESCRIPTION: Subroutine PHOTO\_JX is the core subroutine of Fast-JX.
!    calc J's for a single column atmosphere (aka Indep Colm Atmos or ICA)
!    needs P, T, O3, clds, aersls; adds top-of-atmos layer from climatology
!    needs day-fo-year for sun distance, SZA (not lat or long)
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE PHOTO_JX( U0,             REFLB,     P_COL,
     &                     T_COL,          O3_COL,    O3_TOMS,
     &                     AOD999,         ODAER_COL, ODMDUST_COL,
     &                     ODCLOUD_COL_IN, ILON,      ILAT,
     &                     YLAT,           DAY_OF_YR, MONTH,
     &                     DAY,            am_I_Root, Input_Opt )
!
! !USES:
!
      USE CMN_DIAG_MOD,       ONLY : LD64, ND64
      USE DIAG_MOD,           ONLY : AD64
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
!
! !INPUT PARAMETERS:
!
      INTEGER, INTENT(IN)                           :: ILON,ILAT
      INTEGER, INTENT(IN)                           :: DAY_OF_YR
      INTEGER, INTENT(IN)                           :: MONTH,DAY
      REAL*8,  INTENT(IN)                           :: YLAT
      REAL*8,  INTENT(IN)                           :: U0,REFLB
      REAL*8,  INTENT(IN),DIMENSION(L1_+1)          :: P_COL
      REAL*8,  INTENT(IN),DIMENSION(L1_)            :: T_COL
      REAL*8,  INTENT(IN),DIMENSION(L1_)            :: O3_COL
      REAL*8,  INTENT(IN)                           :: O3_TOMS
      REAL*8,  INTENT(IN),DIMENSION(L_,NALLAER)     :: ODAER_COL
      REAL*8,  INTENT(IN),DIMENSION(L_,NDUST)       :: ODMDUST_COL
      REAL*8,  INTENT(IN),DIMENSION(L_)             :: ODCLOUD_COL_IN
      LOGICAL, INTENT(IN)                           :: AOD999
      LOGICAL, INTENT(IN)                           :: am_I_Root
      TYPE(OptInput), INTENT(IN)                    :: Input_Opt
!
! !REMARKS:
! 
! !REVISION HISTORY: 
!  28 Mar 2013 - S. D. Eastham - Copied from Fast-JX v7.0
!  05 Mar 2014 - R. Yantosca - Now only save ND64 with LD64 specified levels
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
c-----------------------------------------------------------------------
c---key LOCAL atmospheric data needed to solve plane-parallel J----
c-----these are dimensioned JXL_, and must have JXL_ .ge. LLPAR
      real*8, dimension(JXL1_)      :: DDJ,OOJ
      real*8, dimension(JXL1_+1)    :: PPJ,ZZJ
      integer,dimension(JXL2_+1)    :: JXTRA
      real*8, dimension(W_)         :: FJTOP,FJBOT,FSBOT,FLXD0,RFL
      real*8, dimension(JXL_, W_)   :: AVGF, FJFLX
      real*8, dimension(JXL1_,W_)   :: DTAUX, FLXD
      real*8, dimension(8,JXL1_,W_) :: POMEGAX

c---flux/heating arrays (along with FJFLX,FLXD,FLXD0)
      real*8 :: ODABS,ODRAY
      real*8 :: RFLECT
      real*8 :: AMF2(2*JXL1_+1,2*JXL1_+1)

c------------key SCATTERING arrays for clouds+aerosols------------------
      real*8 :: OD(5,JXL1_),SSA(5,JXL1_),SLEG(8,5,JXL1_)
      real*8 :: OD600(JXL1_)

c------------key arrays AFTER solving for J's---------------------------
      real*8   FFF(W_,JXL_),VALJ(X_)
      real*8   VALJL(JXL_,X_) !2-D array of J_s returned by JRATET

      integer  L2EDGE, I,J,K,L,M,KMIE,IDXAER
      real*8   XQO3,XQO2,WAVE, TTTX
C-----------------------------------------------------------------------
      ! For compatibility with GEOS-Chem (SDE 03/30/13)
      REAL*8   QSCALING,LOCALOD,LOCALSSA,SZA,SOLF
      REAL*8   AERX_COL(A_,L1_) ! Accumulated aerosol array
      REAL*8   ODCLOUD_COL(L_)
      REAL*8   VALJXX(L_,JVN_)

      ! T_CLIM: Climatology data temperature (K)
      REAL*8   T_CLIM(JXL1_)

      ! T_INPUT: Input temperature (K) with extra layer for compatibility
      REAL*8   T_INPUT(JXL1_+1)

      ! AD64 diagnostics
      REAL*8, DIMENSION(JXL1_) :: FDIRECT, FDIFFUSE

      !=================================================================
      ! PHOTO_JX begins here!
      !=================================================================

      if (LLPAR+1 .gt. JXL1_) then
        call EXITC(' PHOTO_JX: not enough levels in JX')
      endif

      ! Copy cloud OD data to a variable array
      DO L=1,L_
         ODCLOUD_COL(L) = ODCLOUD_COL_IN(L)
      ENDDO

      ! Input conversion (SDE 03/29/13)
      ! Calculate solar zenith angle (degrees)
      CALL SOLAR_JX(DAY_OF_YR,U0,SZA,SOLF)

      L2EDGE = L_ + L_ + 2
      FFF(:,:) = 0.d0

c---check for dark conditions SZA > 98.0 deg => tan ht = 63 km
c                        or         99.0                 80 km
      if (SZA .gt. 98.d0) goto 99

      ! Use GEOS-Chem methodology to set vertical profiles of:
      ! Pressure      (PPJ)    [hPa]
      ! Temperature   (T_CLIm) [K]
      ! Path density  (DDJ)    [# molec/cm2]
      ! New methodology for:
      ! Ozone density (OOJ)    [# O3 molec/cm2]
      CALL SET_PROF (YLAT,        MONTH,     DAY,
     &               T_COL,       P_COL,     ODCLOUD_COL,
     &               ODMDUST_COL, ODAER_COL, O3_COL,
     &               O3_TOMS,     AERX_COL,  T_CLIM,
     &               OOJ,         ZZJ,       DDJ,
     &               Input_Opt )

      ! Fill out PPJ and TTJ with CTM data to replace fixed climatology
      DO L=1,L1_
         PPJ(L) = P_COL(L)
         T_INPUT(L) = T_COL(L)
      ENDDO
      
      ! Ensure TOA pressure is zero
      PPJ(L1_+1) = 0.d0
      T_INPUT(L1_+1) = T_CLIM(L1_)

c---calculate spherical weighting functions (AMF: Air Mass Factor)
        RFLECT = REFLB

C-----------------------------------------------------------------------
      call SPHERE2 (U0,ZZJ,AMF2,L1_,JXL1_)
C-----------------------------------------------------------------------

!-----------------------------------------------------------------------
! Modification for GEOS-Chem: Optical depths are calculated at a single
! wavelength for GEOS-Chem, so we perform scaling in this routine.
! elsewhere.
! (SDE 03/31/13)
!-----------------------------------------------------------------------
c---calculate the optical properties (opt-depth, single-scat-alb, phase-fn(1:8))
c---  at the 5 std wavelengths 200-300-400-600-999 nm for cloud+aerosols
      do L = 1,L1_
       OD600(L) = 0.d0
       ! ODs stored with fine-grain (DTAUX) and coarse (OD)
       do K=1,W_
         DTAUX(L,K) = 0.d0
       enddo
       do K=1,5
         OD(K,L)  = 0.d0
         SSA(K,L) = 0.d0
        do I=1,8
         SLEG(I,K,L) = 0.d0
        enddo
       enddo
      enddo

       ! Don't bother on extraneous level - leave as zero
       do L = 1,L_
          DO KMIE=1,5
             ! Clouds, aerosols, mineral dust
             DO M=1,AN_
                IF (AERX_COL(M,L).gt.0d0) THEN
                   IDXAER=MIEDX(M)
                   IF (IDXAER.lt.15) THEN
                      ! Cloud (600 nm scaling)
                      QSCALING = QAA(KMIE,IDXAER)/QAA(4,IDXAER)
                   ELSEIF (AOD999) THEN
                      ! Aerosol/dust (999 nm scaling)
                      QSCALING = QAA(KMIE,IDXAER)/QAA(IND999,IDXAER)
                   ELSE
                      ! Aerosol/dust (550 nm scaling)
                      QSCALING = QAA(KMIE,IDXAER)/QAA_AOD(IDXAER)
                   ENDIF
                   LOCALOD = QSCALING*AERX_COL(M,L)
                   LOCALSSA = SAA(KMIE,IDXAER)*LOCALOD
                   OD(KMIE,L) = OD(KMIE,L) + LOCALOD
                   SSA(KMIE,L)= SSA(KMIE,L) + LOCALSSA
                   DO I=1,8
                      SLEG(I,KMIE,L) = SLEG(I,KMIE,L) + 
     &                            (PAA(I,KMIE,IDXAER)*LOCALSSA)
                   ENDDO ! I (Phase function)
                ENDIF
             ENDDO ! M (Aerosol)
          ENDDO ! KMIE (Mie scattering wavelength bin)

          ! Normalize
          DO KMIE=1,5
             IF (OD(KMIE,L).gt.0.d0) THEN
                SSA(KMIE,L) = SSA(KMIE,L)/OD(KMIE,L)
                DO I=1,8
                   SLEG(I,KMIE,L) = SLEG(I,KMIE,L)/OD(KMIE,L)
                ENDDO
             ENDIF
          ENDDO
          ! Retrieve 600 nm OD to determine added layers
          OD600(L) = OD(4,L)
       ENDDO ! L (Layer)

c---when combining with Rayleigh and O2-O3 abs, remember the SSA and
c---  phase fn SLEG are weighted by OD and OD*SSA, respectively.
c---Given the aerosol+cloud OD/layer in visible (600 nm) calculate how to add
C       additonal levels at top of clouds (now uses log spacing)

C-----------------------------------------------------------------------
      call EXTRAL(OD600,L1_,L2EDGE,N_,JXTRA)
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
        if( WAVE .gt. 500.d0 ) KMIE=4
        if( WAVE .gt. 800.d0 ) KMIE=5

c---Combine: Rayleigh scatters & O2 & O3 absorbers to get optical properties
c---values at L1_=L_+1 are a pseudo/climatol layer above the top CTM layer (L_)
        do L = 1,L1_
         TTTX     = T_CLIM(L) ! Following GEOS-Chem v9-1-3
         call X_interp (TTTX,XQO2, TQQ(1,1),QO2(K,1), TQQ(2,1),QO2(K,2),
     &         TQQ(3,1),QO2(K,3), LQQ(1))
         call X_interp (TTTX,XQO3, TQQ(1,2),QO3(K,1), TQQ(2,2),QO3(K,2),
     &         TQQ(3,2),QO3(K,3), LQQ(2))
         ODABS = XQO3*OOJ(L) + XQO2*DDJ(L)*0.20948d0
         ODRAY = DDJ(L)*QRAYL(K)

         DTAUX(L,K) = OD(KMIE,L) + ODABS + ODRAY

         ! Aerosols + clouds + O2 + O3
         do I=1,8
           POMEGAX(I,L,K) = SLEG(I,KMIE,L)*OD(KMIE,L)
         enddo
         ! Add Rayleigh scattering effects
         ! Only non-zero for 1st and 3rd phase functions
         POMEGAX(1,L,K) = POMEGAX(1,L,K) + 1.0d0*ODRAY
         POMEGAX(3,L,K) = POMEGAX(3,L,K) + 0.5d0*ODRAY
         ! Normalize
         do I=1,8
           POMEGAX(I,L,K) = POMEGAX(I,L,K)/DTAUX(L,K)
         enddo
        enddo
      enddo
C-----------------------------------------------------------------------
      call OPMIE (DTAUX,POMEGAX,U0,RFL,AMF2,JXTRA,
     &        AVGF,FJTOP,FJBOT,FSBOT,FJFLX,FLXD,FLXD0,LLPAR)

!C-----------------------------------------------------------------------

      do K = 1,W_

        do L = 1,L_
           FFF(K,L) = FFF(K,L) + SOLF*FL(K)*AVGF(L,K)
        enddo

      enddo

      ! Calculate photolysis rates
      call JRATET(PPJ,T_INPUT,FFF, VALJXX,L_,LLCHEM,NJX)

      ! Fill out common-block array of J-rates
      DO L=1,LLCHEM
         DO J=1,NRATJ
            IF (JIND(J).gt.0) THEN
               ZPJ(L,J,ILON,ILAT) = VALJXX(L,JIND(J))*JFACTA(J)
            ELSE
               ZPJ(L,J,ILON,ILAT) = 0.d0
            ENDIF
         ENDDO
      ENDDO
      
      ! Set J-rates outside the chemgrid to zero
      IF (LLCHEM.lt.L_) THEN
         DO L=LLCHEM+1,L_
            DO J=1,NRATJ
               ZPJ(L,J,ILON,ILAT) = 0.d0
            ENDDO
         ENDDO
      ENDIF

      !---------------------------------------------------------------
      ! New diagnostic (SDE 2013-10-25) 
      ! Use it to calculate fluxes for output if necessary
      !---------------------------------------------------------------
      IF ( ND64 > 0 .and. L <= LD64 ) THEN
         ! Get net direct and net diffuse fluxes separately
         ! Order:
         !    1 - Net flux
         !    2 - Direct flux
         !    3 - Diffuse flux
         ! Convention: -ve is downwards
         DO K=1,W_
            ! Direct flux at the bottom
            FDIRECT(:) = 0.0
            FDIFFUSE(:) = 0.0
            FDIRECT(1) = FSBOT(K)
            FDIFFUSE(1) = FJBOT(K)
            DO L=2,LD64
               FDIRECT(L) = FDIRECT(L-1) + FLXD(L-1,K)
               FDIFFUSE(L) = FJFLX(L-1,K)
            ENDDO
            DO L=1,L1_
               AD64(ILON,ILAT,L,K,1) = AD64(ILON,ILAT,L,K,1) +
     &            ((FDIRECT(L)+FDIFFUSE(L))*SOLF*FL(K)*UVXFACTOR(K))
               AD64(ILON,ILAT,L,K,2) = AD64(ILON,ILAT,L,K,2) +
     &            (FDIRECT(L)*SOLF*FL(K)*UVXFACTOR(K))
               AD64(ILON,ILAT,L,K,3) = AD64(ILON,ILAT,L,K,3) +
     &            (FDIFFUSE(L)*SOLF*FL(K)*UVXFACTOR(K))
            ENDDO
         ENDDO
      ENDIF

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%% COMMENT OUT FAST-J DIAGNOSTICS FOR NOW (bmy, 3/5/14)
!%%% LEAVE CODE HERE SO THAT IT CAN BE RESTORED 
!      !---------------------------------------------------------------
!      ! Majority of heating code ignored by GEOS-Chem
!      ! Everything from here until statement number 99 is ignored
!      ! (SDE 03/31/13)
!      !---------------------------------------------------------------
!      IF (LFJXDIAG) THEN
!         DO K=1,W_
!c----direct(DIR) and diffuse(FLX) fluxes at top(UP) (solar = negative by convention)
!c----     also at bottom (DN), does not include diffuse reflected flux.
!            FLXUP(K) =  FJTOP(K)
!            DIRUP(K) = -FLXD0(K)
!            FLXDN(K) = -FJBOT(K)
!            DIRDN(K) = -FSBOT(K)
!
!            FREFI = FREFI + SOLF*FL(K)*FLXD0(K)/WL(K)
!            FREFL = FREFL + SOLF*FL(K)*FJTOP(K)/WL(K)
!            FREFS = FREFS + SOLF*FL(K)/WL(K)
!
!c---for each wavelength calculate the flux budget/heating rates:
!c  FLXD(L) = direct flux deposited in layer L  [approx = MU0*(F(L+1) -F(L))]
!c            but for spherical atmosphere!
!c  FJFLX(L) = diffuse flux across top of layer L
!
!c---calculate divergence of diffuse flux in each CTM layer (& t-o-a)
!c---     need special fix at top and bottom:
!c---FABOT = total abs at L.B. &  FXBOT = net diffusive flux at L.B.
!            FABOT = (1.d0-RFL(K))*(FJBOT(K)+FSBOT(K))
!            FXBOT = -FJBOT(K) + RFL(K)*(FJBOT(K)+FSBOT(K))
!            FLXJ(1) = FJFLX(1,K) - FXBOT
!            do L=2,LU
!               FLXJ(L) = FJFLX(L,K) - FJFLX(L-1,K)
!            enddo
!            FLXJ(LU+1) = FJTOP(K) - FJFLX(LU,K)
!c---calculate net flux deposited in each CTM layer (direct & diffuse):
!            FFX0 = 0.d0
!            do L=1,L1U
!               FFX(K,L) = FLXD(L,K) - FLXJ(L)
!               FFX0 = FFX0 + FFX(K,L)
!            enddo
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
!            FFXNET(K,1) = FLXD0(K)
!            FFXNET(K,2) = FSBOT(K)
!            FFXNET(K,3) = FLXD0(K) + FSBOT(K)
!            FFXNET(K,4) = FJTOP(K)
!            FFXNET(K,5) = FFX0
!            FFXNET(K,6) = FABOT
!            FFXNET(K,7) = FSBOT(K)
!            FFXNET(K,8) = FJBOT(K)
!
!c-----------------------------------------------------------------------
!         enddo       ! end loop over wavelength K
!c-----------------------------------------------------------------------
!         FREFL = FREFL/FREFS      !calculate reflected flux (energy weighted)
!         FREFI = FREFI/FREFS
!
!c---NB UVB = 280-320 = bins 12:15, UVA = 320-400 = bins 16:17, VIS = bin 18 (++)
!
!
!c---mapping J-values from fast-JX onto CTM chemistry is done in main
!
!C-----------------------------------------------------------------------
!         if ((am_I_Root).and.(LPRTJ)) then
!c---diagnostics below are NOT returned to the CTM code
!            write(6,*)'fast-JX-(7.0)---PHOTO_JX internal print:',
!     &               ' Atmosphere---'
!c---used last called values of DTAUX and POMEGAX, should be 600 nm
!         do L=1,L1U
!            DTAU600(L) = DTAUX(L,W_)
!            do I=1,8
!               POMG600(I,L) = POMEGAX(I,L,W_)
!            enddo
!         enddo
!
!      !   call JP_ATM(PPJ,TTJ,DDJ,OOJ,ZZJ,DTAU600,POMG600,JXTRA, LU)
!
!C---PRINT SUMMARY of mean intensity, flux, heating rates:
!         if (am_I_Root) then
!            write(6,*)
!            write(6,*)'fast-JX(7.0)---PHOTO_JX internal print:',
!     &               ' Mean Intens---'
!            write(6,'(a,5f10.4)')
!     &       ' SUMMARY fast-JX: albedo/SZA/u0/F-incd/F-refl/',
!     &        RFLECT,SZA,U0,FREFI,FREFL
!
!            write(6,'(a5,18i8)')   ' bin:',(K, K=NW2,NW1,-1)
!            write(6,'(a5,18f8.1)') ' wvl:',(WL(K), K=NW2,NW1,-1)
!            write(6,'(a,a)') ' ----  100000=Fsolar   ',
!     &                  'MEAN INTENSITY per wvl bin'
!         endif
!         do L = LU,1,-1
!            do K=NW1,NW2
!               RATIO(K) = (1.d5*FFF(K,L)/FL(K))
!            enddo
!            if (am_I_Root) then
!               write(6,'(i3,2x,18i8)') L,(RATIO(K),K=NW2,NW1,-1)
!            endif
!         enddo
!
!         if (am_I_Root) then
!            write(6,*)
!            write(6,*)'fast-JX(7.0)---PHOTO_JX internal print:',
!                              ' Net Fluxes---'
!            write(6,'(a11,18i8)')   ' bin:',(K, K=NW2,NW1,-1)
!            write(6,'(a11,18f8.1)') ' wvl:',(WL(K), K=NW2,NW1,-1)
!c            write(6,'(a11,18f8.4)') ' sol in atm',(FFXNET(K,1), 
!c     &                        K=NW2,NW1,-1)
!c            write(6,'(a11,18f8.4)') ' sol at srf',(FFXNET(K,2), 
!c     &                        K=NW2,NW1,-1)
!            write(6,*) ' ---NET FLUXES--- '
!            write(6,'(a11,18f8.4)') ' sol TOTAL ',(FFXNET(K,3), 
!     &                        K=NW2,NW1,-1)
!            write(6,'(a11,18f8.4)') ' dif outtop',(FFXNET(K,4), 
!     &                        K=NW2,NW1,-1)
!            write(6,'(a11,18f8.4)') ' abs in atm',(FFXNET(K,5), 
!     &                        K=NW2,NW1,-1)
!            write(6,'(a11,18f8.4)') ' abs at srf',(FFXNET(K,6), 
!     &                        K=NW2,NW1,-1)
!            write(6,*) ' ---SRF FLUXES--- '
!            write(6,'(a11,18f8.4)') ' srf direct',(FFXNET(K,7), 
!     &                        K=NW2,NW1,-1)
!            write(6,'(a11,18f8.4)') ' srf diffus',(FFXNET(K,8), 
!     &                        K=NW2,NW1,-1)
!            write(6,'(4a)') '  ---NET ABS per layer:',
!     &         '       10000=Fsolar',
!     &         '  [NB: values <0 = numerical error w/clouds',
!     &         ' or SZA>90, colm OK]'
!         endif
!         do L = LU,1,-1
!            do K=NW1,NW2
!               RATIO(K) = 1.d5*FFX(K,L)
!            enddo
!            if (am_I_Root) then
!               write(6,'(i9,2x,18i8)') L,(RATIO(K),K=NW2,NW1,-1)
!            endif
!         enddo
!         if (am_I_Root) then
!            write(6,'(a)')
!            write(6,'(a)') ' fast-JX (7.0)----J-values----'
!            write(6,'(1x,a,72(a6,3x))') 'L=  ',(TITLEJX(K), K=1,NJX)
!            do L = LU,1,-1
!              write(6,'(i3,1p, 72e9.2)') L,(VALJXX(L,K),K=1,NJX)
!            enddo
!         endif
!
!      ENDIF
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


   99 continue

      END SUBROUTINE PHOTO_JX
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
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
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: miesct
!
! !DESCRIPTION: Subroutine MIESCT is an adaptation of the Prather radiative
!  transfer code (mjp, 10/95).
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
!     Prather, 1974, Astrophys. J. 192, 787-792.
!         Solution of inhomogeneous Rayleigh scattering atmosphere.
!         (original Rayleigh w/ polarization)
!     Cochran and Trafton, 1978, Ap.J., 219, 756-762.
!         Raman scattering in the atmospheres of the major planets.
!         (first use of anisotropic code)
!     Jacob, Gottlieb and Prather, 1989, J.Geophys.Res., 94, 12975-13002.
!         Chemistry of a polluted cloudy boundary layer,
!         (documentation of extension to anisotropic scattering)
!                                                                             .
!    takes atmospheric structure and source terms from std J-code
!    ALSO limited to 4 Gauss points, only calculates mean field! (M=1)
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
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
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
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: set_prof
!
! !DESCRIPTION: Subroutine SET\_PROF sets vertical profiles for a given
!  latitude and longitude. 
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE SET_PROF (YLAT,     MONTH,  DAY,     T_CTM,  P_CTM,
     &                     CLDOD,    DSTOD,  AEROD,   O3_CTM, O3_TOMS,
     &                     AERCOL,   T_CLIM, O3_CLIM, Z_CLIM, AIR_CLIM,
     &                     Input_Opt )
!
! !USES:
!
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
!
! !INPUT PARAMETERS:
!
      REAL*8, INTENT(IN)         :: YLAT              ! Latitude (degrees)
      INTEGER,INTENT(IN)         :: MONTH             ! Month
      INTEGER,INTENT(IN)         :: DAY               ! Day *of month*
      REAL*8, INTENT(IN)         :: T_CTM(L1_)        ! CTM temperatures (K)
      REAL*8, INTENT(IN)         :: O3_TOMS           ! O3 column (DU)
      REAL*8, INTENT(IN)         :: P_CTM(L1_+1)      ! CTM edge pressures (hPa)
      REAL*8, INTENT(INOUT)      :: CLDOD(L_)         ! Cloud optical depth
      REAL*8, INTENT(IN)         :: DSTOD(L_,NDUST)   ! Mineral dust OD
      REAL*8, INTENT(IN)         :: AEROD(L_,NALLAER) ! Aerosol OD
      REAL*8, INTENT(IN)         :: O3_CTM(L1_)       ! CTM ozone (molec/cm3)
      TYPE(OptInput), INTENT(IN) :: Input_Opt         ! Input options
!
! !OUTPUT VARIABLES:
!
      REAL*8, INTENT(OUT)        :: AERCOL(A_,L1_)    ! Aerosol column
      REAL*8, INTENT(OUT)        :: T_CLIM(L1_)       ! Clim. temperatures (K)
      REAL*8, INTENT(OUT)        :: Z_CLIM(L1_+1)     ! Edge altitudes (cm)
      REAL*8, INTENT(OUT)        :: O3_CLIM(L1_)      ! O3 column depth (#/cm2)
      REAL*8, INTENT(OUT)        :: AIR_CLIM(L1_)     ! O3 column depth (#/cm2)
!
! !REMARKS:
! 
! !REVISION HISTORY: 
!  30 Mar 2013 - S. D. Eastham - Adapted from J. Mao code
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER                :: I, K, L, M, N, LCTM
      REAL*8                 :: DLOGP,F0,T0,B0,PB,PC,XC,MASFAC,SCALEH
      REAL*8                 :: PSTD(52),OREF2(51),TREF2(51)
      REAL*8                 :: PROFCOL, ODSUM
      REAL*8, PARAMETER      :: ODMAX = 200.0d0

      ! Local variables for quantities from Input_Opt
      LOGICAL :: LO3FJX

      !=================================================================
      ! SET_PROF begins here!
      !=================================================================

      ! Copy fields from INPUT_OPT
      LO3FJX = Input_Opt%LO3FJX

      ! Zero aerosol column
      DO K=1,A_
         DO I=1,L1_   
            AERCOL(K,I) = 0.d0
         ENDDO
      ENDDO

      ! Scale optical depths to stay within limits
      ODSUM = 0.d0
      DO I=1,L_   
         CLDOD(I) = DBLE(CLDOD(I))
         ODSUM = ODSUM + CLDOD(I) 
      ENDDO
      IF (ODSUM.gt.ODMAX) THEN
         ODSUM = ODMAX/ODSUM ! Temporary
         DO I=1,L_    
            CLDOD(I) = CLDOD(I)*ODSUM
         ENDDO
         ODSUM = ODMAX
      ENDIF

      !=================================================================      
      ! Set up pressure levels for O3/T climatology - assume that value
      ! given for each 2 km z* level applies from 1 km below to 1 km 
      ! above, so select pressures at these boundaries. Surface level 
      ! values at 1000 mb are assumed to extend down to the actual 
      ! surface pressure for this lat/lon.
      !=================================================================      
      PSTD(1)  = MAX(P_CTM(1),1000.d0)
      PSTD(2)  = 1000.d0 * 10.d0 ** (-1.d0/16.d0)
      DLOGP    = 10.d0**(-2.d0/16.d0)
      DO I=3,51
         PSTD(I) = PSTD(I-1) * DLOGP
      ENDDO
      PSTD(52) = 0.d0

      ! Mass factor - delta-Pressure [hPa] to delta-Column [molec/cm2]
      MASFAC = 100.d0 * 6.022d+23 / ( 28.97d0 * 9.8d0 * 10.d0 )

      ! Select appropriate monthly and latitudinal profiles
      ! Now use YLAT instead of Oliver's YDGRD(NSLAT) (bmy, 9/13/99) 
      M = MAX( 1, MIN( 12, MONTH                   ) )
      L = MAX( 1, MIN( 18, ( INT(YLAT) + 99 ) / 10 ) )

      ! Temporary arrays for climatology data
      DO I = 1, 51
          OREF2(I) = OREF(I,L,M)
          TREF2(I) = TREF(I,L,M)
      ENDDO
 
      ! Apportion O3 and T on supplied climatology z* levels onto CTM levels 
      ! with mass (pressure) weighting, assuming constant mixing ratio and
      ! temperature half a layer on either side of the point supplied.

      DO I = 1, L1_  
         F0 = 0.D0
         T0 = 0.D0
         DO K = 1, 51
            PC = MIN( P_CTM(I),   PSTD(K)   )
            PB = MAX( P_CTM(I+1), PSTD(K+1) )
            IF ( PC .GT. PB ) THEN
               XC = ( PC - PB ) / ( P_CTM(I) - P_CTM(I+1) )
               F0 = F0 + OREF2(K)*XC
               T0 = T0 + TREF2(K)*XC
            ENDIF
         ENDDO
         T_CLIM(I)  = T0
         O3_CLIM(I) = F0 * 1.D-6
      ENDDO

      !=================================================================
      ! Calculate effective altitudes using scale height at each level
      !=================================================================
      Z_CLIM(1) = 0.D0
      DO I = 1, L_   
         SCALEH = 1.3806D-19 * MASFAC * T_CLIM(I)
         Z_CLIM(I+1) = Z_CLIM(I) - 
     &        ( LOG( P_CTM(I+1) / P_CTM(I) ) * SCALEH )
      ENDDO
      Z_CLIM(L1_+1)=Z_CLIM(L1_) + ZZHT

      !=================================================================
      ! Add Aerosol Column - include aerosol types here. Currently use 
      ! soot water and ice; assume black carbon x-section of 10 m2/g, 
      ! independent of wavelength; assume limiting temperature for 
      ! ice of -40 deg C.
      !=================================================================
      DO I = 1, L_   
         ! Turn off uniform black carbon profile (rvm, bmy, 2/27/02)
         AERCOL(1,I) = 0D0

         IF ( T_CTM(I) .GT. 233.D0 ) THEN
            AERCOL(2,I) = CLDOD(I)
            AERCOL(3,I) = 0.D0
         ELSE
            AERCOL(2,I) = 0.D0
            AERCOL(3,I) = CLDOD(I)
         ENDIF

         ! Also add in aerosol optical depth columns (rvm, bmy, 9/30/00)
         DO N = 1, NDUST
            AERCOL(3+N,I) = DSTOD(I,N)
         ENDDO
        
         ! Also add in other aerosol optical depth columns (rvm, bmy, 2/27/02)
         DO N = 1, NALLAER
            AERCOL(3+N+NDUST,I) = AEROD(I,N)
         ENDDO

      ENDDO

      DO K = 1,(3+NDUST+(NALLAER))
         AERCOL(K,L1_    ) = 0.D0
      ENDDO

      !=================================================================
      ! Calculate column quantities for FAST-JX
      !=================================================================
      PROFCOL = 0d0

      DO I = 1, L1_   

         ! Monthly mean air Column [molec/cm2]
         AIR_CLIM(I)  = ( P_CTM(I) - P_CTM(I+1) ) * MASFAC

         ! Monthly mean O3 column [molec/cm2]
         O3_CLIM(I) = O3_CLIM(I) * AIR_CLIM(I)
    
         ! Monthly mean O3 column [DU] 
         PROFCOL = PROFCOL + ( O3_CLIM(I) / 2.69d16 )
      ENDDO

      !! Top values are special (do not exist in CTM data)
      !AIR_CLIM(L1_)     = P_CTM(L1_) * MASFAC
      !O3_CLIM(L1_) = O3_CLIM(L1_) * AIR_CLIM(L1_)
      
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

      ! Updated with UCX
      ! Since we now have stratospheric ozone calculated online, use
      ! this instead of archived profiles for all chemistry-grid cells
      ! The variable O3_CTM is obtained from State_Met%Tracers, and will be 0
      ! outside the chemgrid (in which case we use climatology)

      ! Scale monthly O3 profile to the daily O3 profile (if available)
      DO I = 1, L1_
         IF ((LO3FJX) .and. (I.le.LLCHEM) .and. (O3_CTM(I) > 0d0)) THEN
            ! Convert from molec/cm3 to molec/cm2
            O3_CLIM(I) = O3_CTM(I) * (Z_CLIM(I+1)-Z_CLIM(I))
         ELSEIF (O3_TOMS > 0d0) THEN
            O3_CLIM(I) = O3_CLIM(I) * ( O3_TOMS / PROFCOL )
         ENDIF
      ENDDO

      END SUBROUTINE SET_PROF
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: set_aer
!
! !DESCRIPTION: Subroutine SET\_AER fills out the array MIEDX.
!  Each entry connects a GEOS-Chem aerosol to its Fast-JX counterpart:
!  MIEDX(Fast-JX index) = (GC index)
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE SET_AER(am_I_Root)
!
! !USES:
!
!
! !INPUT PARAMETERS:
!
      LOGICAL, INTENT(IN) :: am_I_Root
!
! !OUTPUT VARIABLES:
!
!
! !REMARKS:
! 
! !REVISION HISTORY: 
!  31 Mar 2013 - S. D. Eastham - Adapted from J. Mao FJX v6.2 implementation
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER               :: I, J, K
      INTEGER               :: IND(NRHAER)

      !=================================================================
      ! SER_AER begins here!
      !=================================================================

      ! Taken from aerosol_mod.F
      IND = (/22,29,36,43,50/)

      DO I=1,AN_
         MIEDX(I) = 0
      ENDDO

!     Select Aerosol/Cloud types to be used - define types here
!     Each of these types must be listed in the order used by OPMIE.F

      ! Clouds
      MIEDX(1)  =  3   !  Black carbon absorber
      MIEDX(2)  = 10   !  Water Cloud (Deirmenjian 8 micron)
      MIEDX(3)  = 14   !  Irregular Ice Cloud (Mishchenko)

      ! Dust
      MIEDX(4)  = 15   !  Mineral Dust  .15 micron    (rvm, 9/30/00)
      MIEDX(5)  = 16   !  Mineral Dust  .25 micron    (rvm, 9/30/00)
      MIEDX(6)  = 17   !  Mineral Dust  .4  micron    (rvm, 9/30/00)
      MIEDX(7)  = 18   !  Mineral Dust  .8  micron    (rvm, 9/30/00)
      MIEDX(8)  = 19   !  Mineral Dust 1.5  micron    (rvm, 9/30/00)
      MIEDX(9)  = 20   !  Mineral Dust 2.5  micron    (rvm, 9/30/00)
      MIEDX(10) = 21   !  Mineral Dust 4.0  micron    (rvm, 9/30/00)

      ! Aerosols
      DO I=1,NRHAER
         DO J=1,NRH
            MIEDX(10+((I-1)*NRH)+J)=IND(I)+J-1
         ENDDO
      ENDDO

#if defined( UCX )
      ! Stratospheric aerosols - SSA/STS and solid PSCs
      MIEDX(10+(NRHAER*NRH)+1) = 4  ! SSA/LBS/STS
      MIEDX(10+(NRHAER*NRH)+2) = 14 ! NAT/ice PSCs
#endif
c
c  Ensure all 'AN_' types are valid selections
      do i=1,AN_
        IF (am_I_Root) write(6,1000) MIEDX(i),TITLEAA(MIEDX(i))
        if(MIEDX(i).gt.NAA.or.MIEDX(i).le.0) then
          if (am_I_Root) then
             write(6,1200) MIEDX(i),NAA
          endif
          CALL EXITC('Bad MIEDX value.')
        endif
      enddo

 1000 format('Using Aerosol type: ',i3,1x,a)
 1200 format('Aerosol type ',i3,' unsuitable; supplied values must be ',
     $       'between 1 and ',i3)

      END SUBROUTINE SET_AER
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
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
!  29 Mar 2013 - S. D. Eastham - Converted from fjfunc into fjxfunc
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

!      IF ( I > IIPAR .OR. J > JJPAR .OR. L > JVL_ ) THEN 
!         CALL EXITC('Invalid grid-box # in call to FJXFUNC')
!      ENDIF

      !=================================================================
      ! RINDEX converts the J-value index as read from "globchem.dat" to 
      ! the J-value index as read from "fjx_j2j.dat". 
      ! (bmy, 10/5/98, SDE 03/29/13)
      !
      ! Make sure that we have taken the proper reaction! 
      !=================================================================
      N = RINDEX(NREAC)

      IF ( (N > NRATJ ) .or. (N .eq. 0) ) THEN
         WRITE(6,*) 'RXN for ',name,', branch ',brch,' not found!'
         CALL EXITC( 'Check FJXFUNC' )
      ENDIF

      !=================================================================
      ! Return the appropriate J-value as the value of the function 
      !=================================================================
      FJXFUNC = ZPJ(L,N,I,J)

!      ! Uncomment following lines for more rigorous checking
!      IF (ISNAN(FJXFUNC)) THEN
!         WRITE(6,*) 'Photolysis error for branch ', brch, ' of ', name
!         WRITE(6,*) 'NaN photorate in grid I J L: ',I, J, L
!         CALL EXITC('Bad photorate')
!      ENDIF

      END FUNCTION FJXFUNC
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: jv_index
!
! !DESCRIPTION: Subroutine JV\_INDEX computes the mapping between the CTM 
!  indices (from "globchem.dat") for J-values to the FAST-JX indices 
!  (from "FJX\_j2j.dat") for J-values.  (bmy, 10/5/98, 4/27/10)
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE JV_INDEX( am_I_Root )
!
! !USES:
!
      USE CMN_SIZE_MOD
      USE COMODE_LOOP_MOD  ! SMVGEAR II arrays

! !INPUT PARAMETERS: 
!
      LOGICAL, INTENT(IN) :: am_I_Root   ! Is this the root CPU?
! 
! !REVISION HISTORY:
!  05 Oct 1998 - R. Yantosca - Initial version
!  (1 ) Assumes the ordering of a species with several branches in 
!        "ratj.d" is the same as in "chem.dat".
!  (2 ) Updated comments, cosmetic changes (bmy, 11/15/01)
!  (3 ) NAMESPEC is now NAMEGAS for SMVGEAR II.   We don't need to reference 
!        CMN anymore. Now loop from NCS = 1..NCSGAS (bdf, bmy, 4/8/03)
!  (4 ) Now reset NCS to NCSURBAN after loop (dbm, bmy, 10/16/06)
!  (5 ) Increase species name length (fp, 6/09)
!  (6 ) Change "Harvard #" to "GEOS-Chem #" in output (bmy, 4/27/10)
!  30 Jul 2012 - R. Yantosca - Now accept am_I_Root as an argument when
!                              running with the traditional driver main.F
!  30 Jul 2012 - R. Yantosca - Added ProTeX headers
!  31 Mar 2013 - S. D. Eastham - Adapted for Fast-JX v7.0
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      ! Local variables
      INTEGER            :: I, IFNC, IBRCH, N, NK

      ! Needed longer name limits (SDE 03/31/13)
      CHARACTER (LEN=14)  :: FULLNAME
      CHARACTER (LEN=10)  :: SPECNAME

      !=================================================================
      ! JV_INDEX begins here!
      !=================================================================

      ! Zero the RINDEX array
      RINDEX(:) = 0

      ! Loop over photolysis rxns (urban chemistry only)
      DO NCS = 1, NCSGAS
      DO I   = 1, NPHOT

         !==============================================================
         ! I        = Index of photo rxns    from "globchem.dat"
         ! NK       = Absolute rxn number (adds offset to I)
         ! SPECNAME = Name of species I,     from "globchem.dat"
         ! IBRCH    = Branch # of species I, from "globchem.dat"
         !==============================================================
         NK       = NRATES(NCS) + I
         FULLNAME = TRIM(NAMEGAS(IRM(1,NK,NCS)))
         ! Get first 10 characters only
         READ (FULLNAME(1:10),"(a10)") SPECNAME
         IFNC     = DEFPRAT(NK,NCS) + 0.01d0
         IBRCH    = 10d0*( DEFPRAT(NK,NCS) - IFNC ) + 0.5d0

         !==============================================================
         ! N      = Index of photolysis reaction in "FJX_j2j.dat"
         ! RNAMES = Name of species N as listed in "FJX_j2j.dat" 
         ! BRANCH = Relative count of reaction in "FJX_j2j.dat"
         !  
         ! If the species names and branch numbers from both "globchem" 
         ! and "FJX_j2j" match, then store N (the FJX index) in the 
         ! Ith element of RINDEX.
         !  
         ! Thus, when looping over I (the "globchem" indices), as is 
         ! done in FJXFUNC, RINDEX(I) will access the correct J-value 
         ! according to the ordering in "FJX_j2j.dat".
         !==============================================================
         DO N = 1, NRATJ
            IF ( SPECNAME == RNAMES(N) .and. IBRCH == BRANCH(N) ) THEN
               RINDEX(I) = N

               ! change formatting for longer name (FP 6/2009)
               IF ( am_I_Root ) THEN
                  WRITE ( 6, 100 ) I,         FULLNAME,  IBRCH, 
     &                             RINDEX(I), RNAMES(N), BRANCH(N)
               ENDIF
 100           FORMAT('GEOS-Chem #: ', i3, 1x, a14, ' Branch: ', i2, 
     &               ' --->  Fast-JX #: ', i3, 1x, a10, ' Branch: ',i2 )
               EXIT
            ENDIF
         ENDDO
      ENDDO  
      ENDDO

      ! Reset NCS to NCSURBAN for safety's sake (bmy, 10/16/06)
      NCS = NCSURBAN

      ! Return to calling program      
      END SUBROUTINE JV_INDEX
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: rd_prof_nc
!
! !DESCRIPTION: Subroutine RAD\_PROF\_NC reads in the reference climatology  
!  from a NetCDF file rather than an ASCII .dat.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE RD_PROF_NC( am_I_Root )
!
! !USES:
!
      ! Modules for netCDF read
      USE m_netcdf_io_open
      USE m_netcdf_io_read
      USE m_netcdf_io_readattr
      USE m_netcdf_io_close

      ! GEOS-Chem modules
      USE DIRECTORY_MOD, ONLY : DATA_DIR_1x1    ! Data directory 
!
! !INPUT PARAMETERS:
!
      LOGICAL, INTENT(IN)         :: am_I_Root  ! Is this the root CPU?
!
! !OUTPUT VARIABLES:
!
!
! !REMARKS:
!  This file was automatically generated by the Perl scripts in the
!  NcdfUtilities package (which ships w/ GEOS-Chem) and was subsequently
!  hand-edited.
!
! !REVISION HISTORY:
!  19 Apr 2012 - R. Yantosca - Initial version
!  30 Jul 2012 - R. Yantosca - Now accept am_I_Root as an argument when
!                              running with the traditional driver main.F
!  10 Apr 2013 - S. D. Eastham - Adapted for Fast-JX v7.0
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      ! netCDF file handle
      INTEGER            :: fId                 ! netCDF file ID
    
      ! Character strings
      CHARACTER(LEN=255) :: nc_dir              ! netCDF directory name
      CHARACTER(LEN=255) :: nc_file             ! netCDF file name
      CHARACTER(LEN=255) :: nc_path             ! netCDF path name
      CHARACTER(LEN=255) :: v_name              ! netCDF variable name 
      CHARACTER(LEN=255) :: a_name              ! netCDF attribute name
      CHARACTER(LEN=255) :: a_val               ! netCDF attribute value
    
      ! Arrays for netCDF start and count values
      INTEGER            :: st3d(3), ct3d(3)    ! For 3D arrays 
   
      !=================================================================
      ! RD_PROF_NC begins here!
      !=================================================================

      !=========================================================================
      ! Open and read data from the netCDF file
      !=========================================================================
    
      ! Directory and file names
      nc_dir  = TRIM( DATA_DIR_1x1 ) // 'FastJ_201204/'      
      nc_file = 'fastj.jv_atms_dat.nc'  
      nc_path = TRIM( nc_dir ) // TRIM( nc_file )
    
      ! Open netCDF file
      CALL Ncop_Rd( fId, TRIM(nc_path) )
      
      ! Echo info to stdout
      IF ( am_I_Root ) THEN
         WRITE( 6, 100 ) REPEAT( '%', 79 )
         WRITE( 6, 110 ) TRIM(nc_file)
         WRITE( 6, 120 ) TRIM(nc_dir)
      ENDIF
    
      !----------------------------------------
      ! VARIABLE: T
      !----------------------------------------
      
      ! Variable name
      v_name = "T"
      
      ! Read T from file
      st3d   = (/  1,  1,  1 /)
      ct3d   = (/ 51, 18, 12 /)
      CALL NcRd( TREF, fId, TRIM(v_name), st3d, ct3d )
      
      ! Read the T:units attribute
      a_name = "units"
      CALL NcGet_Var_Attributes( fId,TRIM(v_name),TRIM(a_name),a_val )
      
      ! Echo info to stdout
      IF ( am_I_Root ) THEN 
         WRITE( 6, 130 ) TRIM(v_name), TRIM(a_val)
      ENDIF
    
      !----------------------------------------
      ! VARIABLE: O3
      !----------------------------------------
      
      ! Variable name
      v_name = "O3"
      
      ! Read O3 from file
      st3d   = (/  1,  1,  1 /)
      ct3d   = (/ 51, 18, 12 /)
      CALL NcRd( OREF, fId, TRIM(v_name), st3d, ct3d )
      
      ! Read the O3:units attribute
      a_name = "units"
      CALL NcGet_Var_Attributes( fId,TRIM(v_name),TRIM(a_name),a_val )
          
      ! Echo info to stdout
      IF ( am_I_Root ) THEN 
         WRITE( 6, 130 ) TRIM(v_name), TRIM(a_val)
      ENDIF
    
      !=================================================================
      ! Cleanup and quit
      !=================================================================
    
      ! Close netCDF file
      CALL NcCl( fId )
    
      ! Echo info to stdout
      IF ( am_I_Root ) THEN 
         WRITE( 6, 140 )
         WRITE( 6, 100 ) REPEAT( '%', 79 )
      ENDIF
    
      ! FORMAT statements
 100  FORMAT( a                                              )
 110  FORMAT( '%% Opening file  : ',         a               )
 120  FORMAT( '%%  in directory : ',         a, / , '%%'     )
 130  FORMAT( '%% Successfully read ',       a, ' [', a, ']' )
 140  FORMAT( '%% Successfully closed file!'                 )

      END SUBROUTINE RD_PROF_NC
!EOC
      END MODULE FAST_JX_MOD
