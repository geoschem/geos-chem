!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: hcox_dust_dead_mod.F
!
! !DESCRIPTION: Module hcox\_dust\_dead\_mod.F contains routines and
!  variables from Charlie Zender's DEAD dust mobilization model.
!  Most routines are from Charlie Zender, but have been modified and/or
!  cleaned up for inclusion into GEOS-Chem.
!\\
!\\
! This is a HEMCO extension module that uses many of the HEMCO core
! utilities.
!\\
!\\
! NOTE: The current (dust) code was validated at 2 x 2.5 resolution.
!  We have found that running at 4x5 we get much lower (~50%) dust
!  emissions than at 2x2.5.  Recommend we either find a way to scale
!  the U* computed in the dust module, or run a 1x1 and store the the
!  dust emissions, with which to drive lower resolution runs.
!  -- Duncan Fairlie, 1/25/07
!\\
!\\
!  (We'll) implement the [dust] code in the standard [GEOS-Chem]
!  model and put a warning about expected low bias when the simulation
!  is run at 4x5.  Whoever is interested in running dust at 4x5 in the
!  future can deal with making the fix.
!  -- Daniel Jacob, 1/25/07
!\\
!\\
! !REFERENCES:
!
! \begin{itemize}
! \item Zender, C. S., Bian, H., and Newman, D.: Mineral Dust Entrainment and
!   Deposition (DEAD) model: Description and 1990s dust climatology,
!   Journal of Geophysical Research: Atmospheres, 108, 2003.
! \end{itemize}
!
! !INTERFACE:
!
      MODULE HCOX_DUSTDEAD_MOD
!
! !USES:
!
      USE HCO_ERROR_MOD
      USE HCO_DIAGN_MOD
      USE HCOX_State_MOD,    ONLY : Ext_State
      USE HCO_STATE_MOD,     ONLY : HCO_State

      IMPLICIT NONE
      PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
      PUBLIC :: HCOX_DustDead_Run
      PUBLIC :: HCOX_DustDead_Init
      PUBLIC :: HCOX_DustDead_Final
      PUBLIC :: HCOX_DustDead_GetFluxTun
!
! !REVISION HISTORY:
!  (1 ) Added parallel DO loop in GET_ORO (bmy, 4/14/04)
!  (2 ) Now references "directory_mod.f" (bmy, 7/20/04)
!  (3 ) Fixed typo in ORO_IS_LND for PGI compiler (bmy, 3/1/05)
!  (4 ) Modified for GEOS-5 and GCAP met fields (swu, bmy, 8/16/05)
!  (5 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (6 ) Now uses GOCART source function (tdf, bmy, 1/25/07)
!  (7 ) Modifications for 0.5 x 0.667 grid (yxw, dan, bmy, 11/6/08)
!  (8 ) Updates for nested grids (amv, bmy, 12/18/09)
!  01 Mar 2012 - R. Yantosca - Now reference new grid_mod.F90
!  25 Nov 2013 - C. Keller   - Now a HEMCO extension
!  06 Oct 2014 - C. Keller   - Allow mass flux tuning factor be set in
!                              configuration file.
!  08 Jul 2015 - M. Sulprizio- Now include dust alkalinity source (tdf 04/10/08)
!  07 Jan 2016 - E. Lundgren - Change dry air gas constant and molec wt to
!                              match GC values and update acc due to gravity
!                              and universal gas constant to NIST 2014 values
!  14 Oct 2016 - C. Keller   - Now use HCO_EvalFld instead of HCO_GetPtr.
!  24 Aug 2017 - M. Sulprizio- Remove support for GEOS-4, GEOS-5, and MERRA
!  17 Oct 2017 - C. Keller   - Now estimate default scale factor from
!                              model resolution.
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !MODULE VARIABLES:
!
      ! Now pack all local variables into customized instance
      TYPE :: MyInst

       ! Fields required by module
       INTEGER                   :: Instance
       INTEGER                   :: ExtNr           ! Extension num for DustDead
       INTEGER                   :: ExtNrAlk        ! Extension num for DustAlk
       INTEGER, ALLOCATABLE      :: HcoIDs(:)       ! tracer IDs for DustDead
       INTEGER, ALLOCATABLE      :: HcoIDsAlk(:)    ! tracer IDs for DustAlk
       REAL*8                    :: FLX_MSS_FDG_FCT

       !---------------------------------------
       ! 2-D pointers pointing to netCDF arrays
       !---------------------------------------

       ! Time-invariant fields
       REAL(hp), POINTER         :: ERD_FCT_GEO  (:,:) => NULL()
!      REAL,    POINTER          :: ERD_FCT_HYDRO(:,:,:,:)
!      REAL,    POINTER          :: ERD_FCT_TOPO (:,:,:,:)
!      REAL,    POINTER          :: ERD_FCT_UNITY(:,:,:,:)
!      REAL,    POINTER          :: MBL_BSN_FCT  (:,:,:,:)

       ! GOCART source function (tdf, bmy, 1/25/07)
       REAL(hp), POINTER          :: SRCE_FUNC(:,:) => NULL()

       ! Land surface that is not lake or wetland (by area)
       REAL(hp), POINTER          :: LND_FRC_DRY  (:,:) => NULL()
       REAL(hp), POINTER          :: MSS_FRC_CACO3(:,:) => NULL()
       REAL(hp), POINTER          :: MSS_FRC_CLY  (:,:) => NULL()
       REAL(hp), POINTER          :: MSS_FRC_SND  (:,:) => NULL()
       REAL(hp), POINTER          :: SFC_TYP      (:,:) => NULL()

       REAL(hp), POINTER          :: VAI_DST(:,:) => NULL()

       ! Time-varying surface info from CTM
!       REAL*8,  ALLOCATABLE :: FLX_LW_DWN_SFC(:,:)
!       REAL*8,  ALLOCATABLE :: FLX_SW_ABS_SFC(:,:)
!       REAL*8,  ALLOCATABLE :: TPT_GND(:,:)
!       REAL*8,  ALLOCATABLE :: TPT_SOI(:,:)
!       REAL*8,  ALLOCATABLE :: VWC_SFC(:,:)

       ! Variables initialized in dst_tvbds_ntp() and dst_tvbds_ini()
!       REAL*8,  ALLOCATABLE :: SRC_STR(:,:)

       ! LSM plant type, 28 land surface types plus 0 for ocean
       ! Also account for 3 different land types in each grid box
       ! NN_SFCTYP denotes the highest possible surface type number.
       ! (ckeller, 07/24/2014)
       INTEGER, ALLOCATABLE :: PLN_TYP(:,:)
       REAL*8,  ALLOCATABLE :: PLN_FRC(:,:)
       REAL*8,  ALLOCATABLE :: TAI(:,:)

       ! Other fields
       REAL*8,  ALLOCATABLE :: DMT_VWR(:)
!       REAL*8,  ALLOCATABLE :: DNS_AER(:)
       REAL*8,  ALLOCATABLE :: OVR_SRC_SNK_FRC(:,:)
       REAL*8,  ALLOCATABLE :: OVR_SRC_SNK_MSS(:,:)
!       INTEGER, ALLOCATABLE :: OROGRAPHY(:,:)
       REAL*8,  ALLOCATABLE :: DMT_MIN(:)
       REAL*8,  ALLOCATABLE :: DMT_MAX(:)
       REAL*8,  ALLOCATABLE :: DMT_VMA_SRC(:)
       REAL*8,  ALLOCATABLE :: GSD_ANL_SRC(:)
       REAL*8,  ALLOCATABLE :: MSS_FRC_SRC(:)
       TYPE(MyInst), POINTER :: NextInst => NULL()
      END TYPE MyInst

      ! Pointer to instances
      TYPE(MyInst), POINTER :: AllInst => NULL()

      !---------------------------------------
      ! MODULE PARAMETER
      !---------------------------------------
      INTEGER, PARAMETER   :: NBINS = 4       ! # of dust bins
      INTEGER, PARAMETER   :: NN_SFCTYP = 28

      ! Fundamental physical constants
      REAL*8,  PARAMETER   :: GAS_CST_UNV      = 8.3144598d0
      REAL*8,  PARAMETER   :: MMW_H2O          = 1.8015259d-02
      REAL*8,  PARAMETER   :: MMW_DRY_AIR      = 28.97d-3
      REAL*8,  PARAMETER   :: CST_VON_KRM      = 0.4d0
      REAL*8,  PARAMETER   :: GRV_SFC          = 9.80665d0
      REAL*8,  PARAMETER   :: GAS_CST_DRY_AIR  = 287.0d0
      REAL*8,  PARAMETER   :: RDS_EARTH        = 6.37122d+6
      REAL*8,  PARAMETER   :: GAS_CST_H2O      = 461.65D0
      REAL*8,  PARAMETER   :: SPC_HEAT_DRY_AIR = 1005.0d0
      REAL*8,  PARAMETER   :: TPT_FRZ_PNT      = 273.15d0

      ! Derived quantities
      REAL*8,  PARAMETER   :: GRV_SFC_RCP      = 1.0d0   / GRV_SFC
      REAL*8,  PARAMETER   :: CST_VON_KRM_RCP  = 1.0d0   / CST_VON_KRM
      REAL*8,  PARAMETER   :: EPS_H2O          = MMW_H2O / MMW_DRY_AIR
      REAL*8,  PARAMETER   :: EPS_H2O_RCP_M1   = -1.0d0  + MMW_DRY_AIR
     &                                                   / MMW_H2O
      REAL*8,  PARAMETER   :: KAPPA_DRY_AIR    = GAS_CST_DRY_AIR
     &                                         / SPC_HEAT_DRY_AIR

      ! Fixed-size grid information
      INTEGER, PARAMETER   :: DST_SRC_NBR      = 3
      INTEGER, PARAMETER   :: MVT              = 14

      CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCOX_DustDead_Run
!
! !DESCRIPTION: Subroutine HcoX\_DustDead\_Run is the driver routine
! for the HEMCO DEAD dust extension.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE HCOX_DustDead_Run( ExtState, HcoState, RC )
!
! !USES:
!
      USE HCO_CALC_MOD,      ONLY : HCO_EvalFld, HCO_CalcEmis
      USE HCO_FLUXARR_MOD,   ONLY : HCO_EmisAdd
      USE HCO_CLOCK_MOD,     ONLY : HcoClock_Get
      USE HCO_CLOCK_MOD,     ONLY : HcoClock_First
!
! !INPUT PARAMETERS:
!
      TYPE(Ext_State), POINTER        :: ExtState    ! Module options
      TYPE(HCO_State), POINTER        :: HcoState    ! Hemco state
!
! !INPUT/OUTPUT PARAMETERS:
!
      INTEGER,         INTENT(INOUT)  :: RC

! !REVISION HISTORY:
!  08 Apr 2004 - T. D. Fairlie - Initial version
!  (1 ) Added OpenMP parallelization, added comments (bmy, 4/8/04)
!  (2 ) Bug fix: DSRC needs to be held PRIVATE (bmy, 4/14/04)
!  (3 ) Now references DATA_DIR from "directory_mod.f" (bmy, 7/20/04)
!  (4 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (5 ) Bug fix: It should be SNOW/1d3 not SNOW*1d3 (tdf, bmy, 11/18/05)
!  (6 ) Updated output statement (bmy, 1/23/07)
!  (7 ) Use SNOMAS (m H2O) for GEOS-5 (bmy, 1/24/07)
!  25 Aug 2010 - R. Yantosca - Treat MERRA in the same way as for GEOS-5
!  25 Aug 2010 - R. Yantosca - Added ProTeX headers
!  03 Sep 2010 - R. Yantosca - Bug fix, SNOMAS was mislabled in GEOS-5
!                              and has units of mm H2O instead of m H2O
!                              so we need to convert to m H2O.
!  08 Feb 2012 - R. Yantosca - Treat GEOS-5.7.x in the same way as MERRA
!  01 Mar 2012 - R. Yantosca - Now use GET_YMID_R(I,J,L) from grid_mod.F90
!  09 Nov 2012 - M. Payer    - Replaced all met field arrays with State_Met
!                              derived type object
!  25 Nov 2013 - C. Keller   - Now a HEMCO extension
!  06 Oct 2014 - C. Keller   - Now calculate pressure center from edges.
!  26 Jun 2015 - E. Lundgren - Add L. Zhang new dust size distribution scheme
!  08 Jul 2015 - M. Sulprizio- Now include dust alkalinity source (tdf 04/10/08)
!  26 Oct 2016 - R. Yantosca - Don't nullify local ptrs in declaration stmts
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      ! Scalars
      LOGICAL                :: ERR
      INTEGER                :: I,      J,      L,       N
      INTEGER                :: M,      IOS,    INC,     LAT_IDX
      INTEGER                :: NDB,    NSTEP,  intDOY
      REAL*8                 :: W10M,   DEN,    DIAM,    U_TS0
      REAL*8                 :: U_TS,   SRCE_P, Reynol,  YMID_R
      REAL*8                 :: ALPHA,  BETA,   GAMMA,   CW
      REAL*8                 :: XTAU,   P1,      P2
      REAL*8                 :: AREA_M2, DTSRCE,  DOY

      ! Arrays
      INTEGER                :: OROGRAPHY(HcoState%NX,HcoState%NY)
      REAL*8                 :: PSLON(HcoState%NX)         ! surface pressure
      REAL*8                 :: PTHICK(HcoState%NX)        ! delta P (L=1)
      REAL*8                 :: PMID(HcoState%NX)          ! mid layer P (L=1)
      REAL*8                 :: TLON(HcoState%NX)          ! temperature (L=1)
      REAL*8                 :: THLON(HcoState%NX)         ! pot. temp. (L=1)
      REAL*8                 :: ULON(HcoState%NX)          ! U-wind (L=1)
      REAL*8                 :: VLON(HcoState%NX)          ! V-wind (L=1)
      REAL*8                 :: BHT2(HcoState%NX)          ! half box height (L=1)
      REAL*8                 :: Q_H2O(HcoState%NX)         ! specific humidity (L=1)
      REAL*8                 :: ORO(HcoState%NX)           ! "orography"
      REAL*8                 :: SNW_HGT_LQD(HcoState%NX)   ! equivalent snow ht.
      REAL*8                 :: DSRC(HcoState%NX,NBINS)    ! dust mixing ratio incr.
      REAL*8                 :: DUST_EMI_TOTAL(HcoState%NX) ! total dust emiss

      ! Flux array [kg/m2/s]
      REAL(hp), TARGET       :: FLUX(HcoState%NX,
     &                               HcoState%NY,
     &                               NBINS)

      ! Flux array for dust alkalinity [kg/m2/s]
      REAL(hp), TARGET       :: FLUX_ALK(HcoState%NX,
     &                                   HcoState%NY,
     &                                   NBINS)

      ! Pointers
      TYPE(MyInst), POINTER  :: Inst

      ! Strings
      CHARACTER(LEN=255)     :: MSG
!
! !DEFINED PARAMETERS:
!
!      REAL*8, PARAMETER      :: Ch_dust = 9.375d-10
!      REAL*8, PARAMETER      :: g0      = 9.80665d0
!      REAL*8, PARAMETER      :: G       = g0 * 1.D2
!      REAL*8, PARAMETER      :: RHOA    = 1.25D-3
      REAL*8, PARAMETER      :: CP      = 1004.16d0
      REAL*8, PARAMETER      :: RGAS    = 8314.3d0 / 28.97d0
      REAL*8, PARAMETER      :: AKAP    = RGAS     / CP
      REAL*8, PARAMETER      :: P1000   = 1000d0

      !=================================================================
      ! HCOX_DUSTDEAD_RUN begins here!
      !=================================================================

      ! Return if extension disabled
      IF ( ExtState%DustDead <= 0 ) RETURN

      ! Enter
      CALL HCO_ENTER( HcoState%Config%Err,
     &               'HCOX_DustDead_Run (hcox_dustdead_mod.F90)', RC )
      IF ( RC /= HCO_SUCCESS ) RETURN

      ! Get instance
      Inst => NULL()
      CALL InstGet ( ExtState%DustDead, Inst, RC )
      IF ( RC /= HCO_SUCCESS ) THEN
       WRITE(MSG,*) 'Cannot find DEAD instance Nr. ', ExtState%DustDead
       CALL HCO_ERROR(HcoState%Config%Err,MSG,RC)
       RETURN
      ENDIF

      !=================================================================
      ! Get pointers to gridded data imported through config. file
      !=================================================================
      !
      ! The following time-invariant fields are read in
      ! ERD_FCT_GEO    ; geomorphic erodibility:       HcoState%NX HcoState%NY
      ! ERD_FCT_HYDRO  ; hydrologic erodibility:       HcoState%NX HcoState%NY
      ! ERD_FCT_TOPO   ; topog. erodibility (Ginoux):  HcoState%NX HcoState%NY
      ! ERD_FCT_UNITY  ; uniform erodibility:          HcoState%NX HcoState%NY
      ! MBL_BSN_FCT    ; overall erodibility factor :  HcoState%NX HcoState%NY
      !
      ! Erodibility field should be copied onto mbl_bsn_fct
      ! which is the one used by the DEAD code   Duncan 8/1/2003
      !
      ! LND_FRC_DRY    ; dry land fraction:            HcoState%NX HcoState%NY
      ! MSS_FRC_CACO3  ; mass fraction of soil CaCO3:  HcoState%NX HcoState%NY
      ! MSS_FRC_CLY    ; mass fraction of clay:        HcoState%NX HcoState%NY
      ! MSS_FRC_SND    ; mass fraction of sand:        HcoState%NX HcoState%NY
      ! SFC_TYP        ; surface type:                 HcoState%NX HcoState%NY
      !=================================================================
      !IF ( HcoClock_First(HcoState%Clock,.TRUE.) ) THEN
         CALL HCO_EvalFld( HcoState, 'DEAD_EF_GEO',
     &                     Inst%ERD_FCT_GEO, RC)
         IF ( RC /= HCO_SUCCESS ) RETURN

         CALL HCO_EvalFld( HcoState, 'DEAD_LF_DRY',
     &                     Inst%LND_FRC_DRY, RC)
         IF ( RC /= HCO_SUCCESS ) RETURN

         CALL HCO_EvalFld( HcoState, 'DEAD_MF_CACO3',
     &                     Inst%MSS_FRC_CACO3,  RC )
         IF ( RC /= HCO_SUCCESS ) RETURN

         CALL HCO_EvalFld( HcoState, 'DEAD_MF_CLY',
     &                     Inst%MSS_FRC_CLY, RC)
         IF ( RC /= HCO_SUCCESS ) RETURN

         CALL HCO_EvalFld( HcoState, 'DEAD_MF_SND',
     &                     Inst%MSS_FRC_SND, RC)
         IF ( RC /= HCO_SUCCESS ) RETURN

         CALL HCO_EvalFld( HcoState, 'DEAD_SFC_TYP',
     &                     Inst%SFC_TYP, RC )
         IF ( RC /= HCO_SUCCESS ) RETURN

         CALL HCO_EvalFld( HcoState, 'DEAD_GOC_SRC',
     &                     Inst%SRCE_FUNC, RC )
         IF ( RC /= HCO_SUCCESS ) RETURN

         CALL HCO_EvalFld( HcoState, 'DEAD_VAI',
     &                     Inst%VAI_DST, RC )
         IF ( RC /= HCO_SUCCESS ) RETURN

!         FIRST = .FALSE.
      !ENDIF

      !=================================================================
      ! CALL DUST MOBILIZATION SCHEME
      !=================================================================

      ! Make OROGRAPHY array from GEOS-CHEM LWI
      CALL GET_ORO( HcoState, ExtState, OROGRAPHY, RC )
      IF ( RC /= HCO_SUCCESS ) RETURN

      ! Get emissions time step
      DTSRCE = HcoState%TS_EMIS

      ! Get day of year, convert to real!!
      CALL HcoClock_Get( HcoState%Clock, cDOY = intDOY, RC=RC )
      IF ( RC /= HCO_SUCCESS ) RETURN
      DOY = intDOY

      ! Init
      FLUX(:,:,:)     = 0.0_hp
      FLUX_ALK(:,:,:) = 0.0_hp

      ! Error check
      ERR = .FALSE.

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I,     J,      P1,    P2,   PTHICK,  PMID, TLON        )
!$OMP+PRIVATE( THLON, ULON,   VLON,  BHT2, Q_H2O,   ORO,  SNW_HGT_LQD )
!$OMP+PRIVATE( N,     YMID_R, DSRC,  RC,   AREA_M2, DUST_EMI_TOTAL    )

      ! Loop over latitudes
      DO J = 1, HcoState%NY

         ! Don't do calculations if there has been an error
         IF ( ERR ) CYCLE

         ! Loop over longitudes
         DO I = 1, HcoState%NX

            ! Pressure [Pa] at bottom and top edge of level 1
            P1             = HcoState%Grid%PEDGE%Val(I,J,1)
            P2             = HcoState%Grid%PEDGE%Val(I,J,2)

            ! Pressure thickness of 1st layer [Pa]
            PTHICK(I)      = ( P1 - P2 )

            ! Pressure at midpt of surface layer [Pa]
            PMID(I)        = ( P1 + P2 ) / 2.0_hp

            ! Temperature [K] at midpoint of surface layer
            TLON(I)        = ExtState%TK%Arr%Val(I,J,1)

            ! Potential temperature [K] at midpoint
            THLON(I)       = TLON(I) * ( P1000 / PMID(I) )**AKAP

            ! U and V winds at surface [m/s]
            ! --> These variables won't be used at all...
            ULON(I)        = ExtState%U10M%Arr%Val(I,J)
            VLON(I)        = ExtState%V10M%Arr%Val(I,J)

            ! Half box height at surface [m]
            BHT2(I)        = HcoState%Grid%BXHEIGHT_M%Val(I,J,1) / 2.d0

            ! Specific humidity at midpoint of surface layer [kg H2O/kg air]
            Q_H2O(I)       = ExtState%SPHU%Arr%Val(I,J,1)

            ! Orography at surface
            ! Ocean is 0; land is 1; ice is 2
            ORO(I)         = REAL(OROGRAPHY(I,J),KIND=8)

            ! Snow [m H2O]. SNOWHGT is in kg H2O/m2, which is equivalent to
            ! mm H2O. Convert to m H2O here.
            SNW_HGT_LQD(I) = ExtState%SNOWHGT%Arr%Val(I,J) / 1000.d0

            ! Dust tracer and increments
            DSRC(I,:) = 0.0d0
         ENDDO !I

         !==============================================================
         ! CALL DUST MOBILIZATION DRIVER (DST_MBL) FOR LATITUDE J
         !==============================================================

         ! Latitude in RADIANS
         YMID_R = HcoState%Grid%YMID%Val(1,J) * HcoState%Phys%PI /180.d0

         ! Call DEAD dust mobilization
         CALL DST_MBL( HcoState, ExtState,  Inst, DOY,
     &                 BHT2,   J,      YMID_R, ORO,
     &                 PTHICK, PMID,  Q_H2O,  DSRC,   SNW_HGT_LQD,
     &                 DTSRCE, TLON,  THLON,  VLON,   ULON,
     &                 J,      RC )

         ! Error check
         IF ( RC /= HCO_SUCCESS ) THEN
            ERR = .TRUE.
            CYCLE
         ENDIF

         ! Redistribute dust emissions using new dust size distribution
         ! scheme (L. Zhang, 6/26/15)
         DUST_EMI_TOTAL = 0.0d0
         DO N = 1, NBINS
            DUST_EMI_TOTAL(:) = DUST_EMI_TOTAL(:) + DSRC(:,N)
         ENDDO
         DSRC(:,1) = DUST_EMI_TOTAL(:) * 0.0766d0
         DSRC(:,2) = DUST_EMI_TOTAL(:) * 0.1924d0
         DSRC(:,3) = DUST_EMI_TOTAL(:) * 0.3491d0
         DSRC(:,4) = DUST_EMI_TOTAL(:) * 0.3819d0

         ! Write to emissions array
         DO I = 1, HcoState%NX

            ! Loop over dust tracers
            ! Write into flux array: kg/box --> kg/m2/s
            AREA_M2 = HcoState%Grid%AREA_M2%Val( I, J )
            DO N = 1, NBINS

               IF ( Inst%HcoIDs(N) > 0 ) THEN
                  FLUX(I,J,N) = ( DSRC(I,N) / AREA_M2 / DTSRCE )
               ENDIF

               ! Include DUST Alkalinity SOURCE, assuming an alkalinity
               ! of 4% by weight [kg].                  !tdf 05/10/08
               !tdf with 3% Ca, there's also 1% equ. Mg, makes 4%
               IF ( Inst%ExtNrAlk > 0 ) THEN
                  FLUX_ALK(I,J,N) = 0.04 * ( DSRC(I,N) / AREA_M2 /
     &                              DTSRCE )
               ENDIF

            ENDDO !N
         ENDDO !I
      ENDDO !J
!$OMP END PARALLEL DO

      ! Error check
      IF ( ERR ) THEN
         RC = HCO_FAIL
         RETURN
      ENDIF

      !=================================================================
      ! PASS TO HEMCO STATE AND UPDATE DIAGNOSTICS
      !=================================================================
      DO N = 1, NBINS

         IF ( Inst%HcoIDs(N) > 0 ) THEN

            ! Add to emissions array
            CALL HCO_EmisAdd( HcoState, FLUX(:,:,N),
     &                        Inst%HcoIDs(N), RC,  ExtNr=Inst%ExtNr )
            IF ( RC /= HCO_SUCCESS ) THEN
               WRITE(MSG,*) 'HCO_EmisAdd error: dust bin ', N
               CALL HCO_ERROR(HcoState%Config%Err,MSG, RC )
               RETURN
            ENDIF

         ENDIF

         IF ( Inst%ExtNrAlk > 0 ) THEN
            IF ( Inst%HcoIDsAlk(N) > 0 ) THEN

               ! Add to dust alkalinity emissions array
               CALL HCO_EmisAdd( HcoState, FLUX_Alk(:,:,N),
     &                           Inst%HcoIDsAlk(N), RC, 
     &                           ExtNr=Inst%ExtNrAlk )
               IF ( RC /= HCO_SUCCESS ) THEN
                  WRITE(MSG,*) 'HCO_EmisAdd error: dust alk bin ', N
                  CALL HCO_ERROR(HcoState%Config%Err,MSG, RC )
                  RETURN
               ENDIF

            ENDIF
         ENDIF

      ENDDO !N

      ! Return w/ success
      Inst => NULL()
      CALL HCO_LEAVE( HcoState%Config%Err, RC )

      END SUBROUTINE HCOX_DustDead_Run
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCOX_DustDead_Init
!
! !DESCRIPTION: Subroutine HcoX\_DustDead\_Init initializes the HEMCO
! DUST\_DEAD extension.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE HCOX_DustDead_Init ( HcoState, ExtName,
     &                                ExtState,  RC                )
!
! !USES:
!
      USE HCO_ExtList_Mod,    ONLY : GetExtNr, GetExtOpt
      USE HCO_STATE_MOD,      ONLY : HCO_GetExtHcoID
!
! !INPUT PARAMETERS:
!
      TYPE(HCO_State),  POINTER        :: HcoState   ! Hemco state
      CHARACTER(LEN=*), INTENT(IN   )  :: ExtName    ! Extension name
      TYPE(Ext_State),  POINTER        :: ExtState     ! Module options
!
! !INPUT/OUTPUT PARAMETERS:
!
      INTEGER,          INTENT(INOUT)  :: RC

! !REVISION HISTORY:
!  25 Nov 2013 - C. Keller   - Now a HEMCO extension
!  14 Aug 2014 - R. Yantosca - Now always print out extension info
!  26 Oct 2016 - R. Yantosca - Don't nullify local ptrs in declaration stmts
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      CHARACTER(LEN=255)             :: MSG
      INTEGER                        :: I, J, N, AS
      INTEGER                        :: ExtNr, nSpc
      INTEGER                        :: nSpcAlk
      CHARACTER(LEN=31), ALLOCATABLE :: SpcNames(:)
      CHARACTER(LEN=31), ALLOCATABLE :: SpcNamesAlk(:)
      REAL(dp)                       :: TmpScal
      LOGICAL                        :: FOUND
      TYPE(MyInst), POINTER          :: Inst

      !=================================================================
      ! HCOX_DUST_DEAD_INIT begins here!
      !=================================================================

      ! Extension Nr.
      ExtNr = GetExtNr( HcoState%Config%ExtList, TRIM(ExtName) )
      IF ( ExtNr <= 0 ) RETURN

      ! Enter
      CALL HCO_ENTER ( HcoState%Config%Err,
     &                'HCOX_DustDead_Init (hcox_dustdead_mod.F90)', RC )
      IF ( RC /= HCO_SUCCESS ) RETURN

      ! Create AeroCom instance for this simulation
      Inst => NULL()
      CALL InstCreate ( ExtNr, ExtState%DustDead, Inst, RC )
      IF ( RC /= HCO_SUCCESS ) THEN
       CALL HCO_ERROR ( HcoState%Config%Err,
     &                 'Cannot create DEAD instance', RC )
       RETURN
      ENDIF

      ! Check for dust alkalinity option
      Inst%ExtNrAlk = GetExtNr( HcoState%Config%ExtList, 'DustAlk')

      ! Get horizontal dimensions
      I = HcoState%NX
      J = HcoState%NY

      !-----------------------------------------------------------------
      ! Get species IDs
      !-----------------------------------------------------------------

      CALL HCO_GetExtHcoID( HcoState, ExtNr, Inst%HcoIDs,
     &                      SpcNames, nSpc, RC)
      IF ( RC /= HCO_SUCCESS ) RETURN

      ! Get the dust alkalinity species defined for DustAlk option
      IF ( Inst%ExtNrAlk > 0 ) THEN
         CALL HCO_GetExtHcoID( HcoState,    Inst%ExtNrAlk,
     &                         Inst%HcoIDsAlk,
     &                         SpcNamesAlk, nSpcAlk,  RC)
         IF ( RC /= HCO_SUCCESS ) RETURN
      ENDIF

      ! Sanity check
      IF ( nSpc /= NBINS ) THEN
         MSG = 'Dust DEAD model does not have four species!'
         CALL HCO_ERROR(HcoState%Config%Err,MSG, RC )
         RETURN
      ENDIF

      ! Set scale factor: first try to read from configuration file. If
      ! not specified, call wrapper function which sets teh scale factor
      ! based upon compiler switches.
      CALL GetExtOpt( HcoState%Config, ExtNr,
     &                 'Mass tuning factor',
     &                 OptValDp=TmpScal, Found=FOUND, RC=RC )
      IF ( RC /= HCO_SUCCESS ) RETURN

      ! Set parameter FLX_MSS_FDG_FCT to specified tuning factor. Get from
      ! wrapper routine if not defined in configuration file
      IF ( FOUND ) THEN
         Inst%FLX_MSS_FDG_FCT = TmpScal
      ELSE
         Inst%FLX_MSS_FDG_FCT = HCOX_DustDead_GetFluxTun( HcoState )
         IF ( Inst%FLX_MSS_FDG_FCT < 0.0_dp ) THEN
            MSG = 'Error in getting dust DEAD mass tuning factor'
            CALL HCO_ERROR(HcoState%Config%Err,MSG,RC)
            RETURN
         ENDIF
      ENDIF

      ! Verbose mode
      IF ( HcoState%amIRoot ) THEN
         MSG = 'Use DEAD dust emissions (extension module)'
         CALL HCO_MSG(HcoState%Config%Err,MSG, SEP1='-' )

         IF ( Inst%ExtNrAlk > 0 ) THEN
            MSG = 'Use dust alkalinity option'
            CALL HCO_MSG(HcoState%Config%Err,MSG, SEP1='-' )
         ENDIF

         MSG = 'Use the following species (Name: HcoID):'
         CALL HCO_MSG(HcoState%Config%Err,MSG)
         DO N = 1, nSpc
            WRITE(MSG,*) TRIM(SpcNames(N)), ':', Inst%HcoIDs(N)
            CALL HCO_MSG(HcoState%Config%Err,MSG)
         ENDDO
         IF ( Inst%ExtNrAlk > 0 ) THEN
            DO N = 1, nSpcAlk
               WRITE(MSG,*) TRIM(SpcNamesAlk(N)), ':', Inst%HcoIDsAlk(N)
               CALL HCO_MSG(HcoState%Config%Err,MSG)
            ENDDO
         ENDIF

         WRITE(MSG,*) 'Global mass flux tuning factor: ',
     &                 Inst%FLX_MSS_FDG_FCT
         CALL HCO_MSG(HcoState%Config%Err,MSG,SEP2='-')

      ENDIF

      !-----------------------------------------------------------------
      ! Init module arrays
      !-----------------------------------------------------------------

      ALLOCATE ( Inst%ERD_FCT_GEO   ( HcoState%NX, HcoState%NY),
     &           Inst%SRCE_FUNC     ( HcoState%NX, HcoState%NY),
     &           Inst%LND_FRC_DRY   ( HcoState%NX, HcoState%NY),
     &           Inst%MSS_FRC_CACO3 ( HcoState%NX, HcoState%NY),
     &           Inst%MSS_FRC_CLY   ( HcoState%NX, HcoState%NY),
     &           Inst%MSS_FRC_SND   ( HcoState%NX, HcoState%NY),
     &           Inst%SFC_TYP       ( HcoState%NX, HcoState%NY),
     &           Inst%VAI_DST       ( HcoState%NX, HcoState%NY),
     &           STAT=AS )
      IF ( AS /= 0 ) THEN
        CALL HCO_ERROR ( HcoState%Config%Err, 'Allocation error', RC )
        RETURN
      ENDIF
      Inst%ERD_FCT_GEO    = 0.0_hp
      Inst%SRCE_FUNC      = 0.0_hp
      Inst%LND_FRC_DRY    = 0.0_hp
      Inst%MSS_FRC_CACO3  = 0.0_hp
      Inst%MSS_FRC_CLY    = 0.0_hp
      Inst%MSS_FRC_SND    = 0.0_hp
      Inst%SFC_TYP        = 0.0_hp
      Inst%VAI_DST        = 0.0_hp

!      ! Allocate arrays
!      ALLOCATE( Inst%FLX_LW_DWN_SFC( I, J ), STAT=AS )
!      IF ( AS /= 0 ) THEN
!        CALL HCO_ERROR ( HcoState%Config%Err, 'FLX_LW_DWN_SFC', RC )
!        RETURN
!      ENDIF
!      Inst%FLX_LW_DWN_SFC = 0d0

!      ALLOCATE( Inst%FLX_SW_ABS_SFC( I, J ), STAT=AS )
!      IF ( AS /= 0 ) THEN
!        CALL HCO_ERROR ( HcoState%Config%Err, 'FLX_SW_ABS_SFC', RC )
!        RETURN
!      ENDIF
!      Inst%FLX_SW_ABS_SFC = 0d0

!      ALLOCATE( Inst%TPT_GND( I, J ), STAT=AS )
!      IF ( AS /= 0 ) THEN
!        CALL HCO_ERROR ( HcoState%Config%Err, 'TPT_GND', RC )
!        RETURN
!      ENDIF
!      Inst%TPT_GND = 0d0

!      ALLOCATE( Inst%TPT_SOI( I, J ), STAT=AS )
!      IF ( AS /= 0 ) THEN
!        CALL HCO_ERROR ( HcoState%Config%Err, 'TPT_SOI', RC )
!        RETURN
!      ENDIF
!      Inst%TPT_SOI = 0d0

!      ALLOCATE( Inst%VWC_SFC( I, J ), STAT=AS )
!      IF ( AS /= 0 ) THEN
!        CALL HCO_ERROR ( HcoState%Config%Err, 'VWC_SFC', RC )
!        RETURN
!      ENDIF
!      Inst%VWC_SFC = 0d0

!      ALLOCATE( Inst%SRC_STR( I, J ), STAT=AS )
!      IF ( AS /= 0 ) THEN
!        CALL HCO_ERROR ( HcoState%Config%Err, 'SRC_STR', RC )
!        RETURN
!      ENDIF
!      Inst%SRC_STR = 0d0

      ALLOCATE( Inst%PLN_TYP( 0:28, 3 ), STAT=AS )
      IF ( AS /= 0 ) THEN
        CALL HCO_ERROR ( HcoState%Config%Err, 'PLN_TYP', RC )
        RETURN
      ENDIF
      Inst%PLN_TYP = 0

      ALLOCATE( Inst%PLN_FRC( 0:28, 3 ), STAT=AS )
      IF ( AS /= 0 ) THEN
        CALL HCO_ERROR ( HcoState%Config%Err, 'PLN_FRC', RC )
        RETURN
      ENDIF
      Inst%PLN_FRC = 0d0

      ALLOCATE( Inst%TAI( MVT, 12 ), STAT=AS )
      IF ( AS /= 0 ) THEN
        CALL HCO_ERROR ( HcoState%Config%Err, 'TAI', RC )
        RETURN
      ENDIF
      Inst%TAI = 0d0

      ALLOCATE( Inst%DMT_VWR( NBINS ), STAT=AS )
      IF ( AS /= 0 ) THEN
        CALL HCO_ERROR ( HcoState%Config%Err, 'DMT_VWR', RC )
        RETURN
      ENDIF
      Inst%DMT_VWR = 0d0

!      ALLOCATE( Inst%DNS_AER( NBINS ), STAT=AS )
!      IF ( AS /= 0 ) THEN
!        CALL HCO_ERROR ( HcoState%Config%Err, 'DNS_AER', RC )
!        RETURN
!      ENDIF
!      Inst%DNS_AER = 0d0

      ALLOCATE( Inst%OVR_SRC_SNK_FRC( DST_SRC_NBR, NBINS ), STAT=AS )
      IF ( AS /= 0 ) THEN
        CALL HCO_ERROR ( HcoState%Config%Err, 'OVR_SRC_SNK_FRC', RC )
        RETURN
      ENDIF
      Inst%OVR_SRC_SNK_FRC = 0d0

      ALLOCATE( Inst%OVR_SRC_SNK_MSS( DST_SRC_NBR, NBINS ), STAT=AS )
      IF ( AS /= 0 ) THEN
        CALL HCO_ERROR ( HcoState%Config%Err, 'OVR_SRC_SNK_MSS', RC )
        RETURN
      ENDIF
      Inst%OVR_SRC_SNK_MSS = 0d0

!      ALLOCATE( Inst%OROGRAPHY( I, J ), STAT=AS )
!      IF ( AS /= 0 ) THEN
!        CALL HCO_ERROR ( HcoState%Config%Err, 'OROGRAPHY', RC )
!        RETURN
!      ENDIF
!      Inst%OROGRAPHY = 0

      ! Bin size min diameter [m]
      ALLOCATE( Inst%DMT_MIN( NBINS ), STAT=AS )
      IF ( AS /= 0 ) THEN
        CALL HCO_ERROR ( HcoState%Config%Err, 'DMT_MIN', RC )
        RETURN
      ENDIF
      Inst%DMT_MIN(1) = 0.2d-6
      Inst%DMT_MIN(2) = 2.0d-6
      Inst%DMT_MIN(3) = 3.6d-6
      Inst%DMT_MIN(4) = 6.0d-6

      ! Bin size max diameter [m]
      ALLOCATE( Inst%DMT_MAX( NBINS ), STAT=AS )
      IF ( AS /= 0 ) THEN
        CALL HCO_ERROR ( HcoState%Config%Err, 'DMT_MAX', RC )
        RETURN
      ENDIF
      Inst%DMT_MAX(1) = 2.0d-6
      Inst%DMT_MAX(2) = 3.6d-6
      Inst%DMT_MAX(3) = 6.0d-6
      Inst%DMT_MAX(4) = 1.2d-5

      ! DMT_VMA_SRC: D'Almeida's (1987) "Background" modes
      ! as default [m]  (Zender et al. p.5 Table 1)
      ! These modes also summarized in BSM96 p. 73 Table 2
      ! Mass median diameter BSM96 p. 73 Table 2
      ALLOCATE( Inst%DMT_VMA_SRC( DST_SRC_NBR ), STAT=AS )
      IF ( AS /= 0 ) THEN
        CALL HCO_ERROR ( HcoState%Config%Err, 'DMT_VMA_SRC', RC )
        RETURN
      ENDIF
      Inst%DMT_VMA_SRC(1) = 0.832d-6
      Inst%DMT_VMA_SRC(2) = 4.82d-6
      Inst%DMT_VMA_SRC(3) = 19.38d-6

      ! GSD_ANL_SRC: Geometric standard deviation [fraction]
      ! BSM96 p. 73 Table 2
      ALLOCATE( Inst%GSD_ANL_SRC( DST_SRC_NBR ), STAT=AS )
      IF ( AS /= 0 ) THEN
        CALL HCO_ERROR ( HcoState%Config%Err, 'GSD_ANL_SRC', RC )
        RETURN
      ENDIF
      Inst%GSD_ANL_SRC(1) = 2.10d0
      Inst%GSD_ANL_SRC(2) = 1.90d0
      Inst%GSD_ANL_SRC(3) = 1.60d0

      ! MSS_FRC_SRC:  Mass fraction BSM96 p. 73 Table 2
      ALLOCATE( Inst%MSS_FRC_SRC( DST_SRC_NBR ), STAT=AS )
      IF ( AS /= 0 ) THEN
        CALL HCO_ERROR ( HcoState%Config%Err, 'MSS_FRC_SRC', RC )
        RETURN
      ENDIF
      Inst%MSS_FRC_SRC(1) = 0.036d0
      Inst%MSS_FRC_SRC(2) = 0.957d0
      Inst%MSS_FRC_SRC(3) = 0.007d0

      !=================================================================
      ! Compute mass overlaps, Mij, between "source" PDFs
      ! and size bins (Zender et al., 2K3, Equ. 12, and Table 1)
      !=================================================================
      CALL OVR_SRC_SNK_FRC_GET( HcoState,
     &                          DST_SRC_NBR,   Inst%DMT_VMA_SRC,
     &                          Inst%GSD_ANL_SRC,   NBINS,
     &                          Inst%DMT_MIN,  Inst%DMT_MAX,
     &                          Inst%OVR_SRC_SNK_FRC, RC )
      IF ( RC /= HCO_SUCCESS ) RETURN

      !=================================================================
      ! Compute OVR_SRC_SNK_MSS, the fraction of dust transported, given
      ! the mass overlap, OVR_SRC_SNK_FRC, and the mass fraction
      ! MSS_FRC_SRC.  OVR_SRC_SNK_MSS is used in routine
      ! FLX_MSS_VRT_DST_PRT which partitions the total vertical
      ! dust flux into transport
      !==============================================================
      CALL DST_PSD_MSS( Inst%OVR_SRC_SNK_FRC, Inst%MSS_FRC_SRC,
     &                  Inst%OVR_SRC_SNK_MSS, NBINS, DST_SRC_NBR )

      !=================================================================
      ! Get plant type, cover, and Leaf area index from land sfc model
      !=================================================================
      CALL PLN_TYP_GET( Inst%PLN_TYP, Inst%PLN_FRC, Inst%TAI )

      ! Activate met fields used by this extension
      ExtState%SPHU%DoUse    = .TRUE.
      ExtState%TK%DoUse      = .TRUE.
      ExtState%U10M%DoUse    = .TRUE.
      ExtState%V10M%DoUse    = .TRUE.
      ExtState%T2M%DoUse     = .TRUE.
      ExtState%GWETTOP%DoUse = .TRUE.
      ExtState%SNOWHGT%DoUse = .TRUE.
      ExtState%USTAR%DoUse   = .TRUE.
      ExtState%Z0%DoUse      = .TRUE.
      ExtState%ALBD%DoUse    = .TRUE.
      ExtState%WLI%DoUse     = .TRUE.

      ! Leave w/ success
      Inst => NULL()
      IF ( ALLOCATED(SpcNames   ) ) DEALLOCATE(SpcNames   )
      IF ( ALLOCATED(SpcNamesAlk) ) DEALLOCATE(SpcNamesAlk)
      CALL HCO_LEAVE( HcoState%Config%Err, RC )

      END SUBROUTINE HCOX_DustDead_Init
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCOX_DustDead_Final
!
! !DESCRIPTION: Subroutine HcoX\_DustDead\_Final finalizes the HEMCO
! DUST\_DEAD extension.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE HCOX_DustDead_Final ( ExtState )
!
! !INPUT PARAMETERS:
!
      TYPE(Ext_State),  POINTER       :: ExtState   ! Module options
!
! !REVISION HISTORY:
!  25 Nov 2013 - C. Keller - Now a HEMCO extension
!
! !NOTES:
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      CALL InstRemove ( ExtState%DustDead )

      END SUBROUTINE HCOX_DustDead_Final
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!
!******************************************************************************
! ORIGINAL ROUTINES FOLLOW BELOW
!******************************************************************************

      SUBROUTINE DST_MBL( HcoState,    ExtState, Inst,
     &                    DOY,         HGT_MDP,     LAT_IDX,
     &                    LAT_RDN,     ORO,         PRS_DLT,
     &                    PRS_MDP,     Q_H2O_VPR,   DSRC,
     &                    SNW_HGT_LQD, TM_ADJ,      TPT_MDP,
     &                    TPT_PTN_MDP, WND_MRD_MDP, WND_ZNL_MDP,
     &                    NSTEP,       RC )
!
!******************************************************************************
!  Subroutine DST_MBL is the driver for aerosol mobilization (DEAD model).
!  It is designed to require only single layer surface fields, allowing for
!  easier implementation.  DST_MBL is called once per latitude.  Modified
!  for GEOS-CHEM by Duncan Fairlie and Bob Yantosca.
!  (tdf, bmy, 1/25/07, 12/18/09)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) DOY         (REAL*8 ) : Day of year [1.0..366.0)            [unitless]
!  (2 ) HGT_MDP     (REAL*8 ) : Midpoint height above surface       [m       ]
!  (3 ) LAT_IDX     (INTEGER) : Model latitude index                [unitless]
!  (4 ) LAT_RDN     (REAL*8 ) : Model latitude                      [radians ]
!  (5 ) ORO         (REAL*8 ) : Orography                           [fraction]
!  (6 ) PRS_DLT     (REAL*8 ) : Pressure thickness of grid box      [Pa      ]
!  (7 ) PRS_MDP     (REAL*8 ) : Pressure @ midpoint of grid box     [Pa      ]
!  (8 ) Q_H2O_VPR,  (REAL*8 ) : Water vapor mixing ratio            [kg/kg   ]
!  (9 ) SNW_HGT_LQD (REAL*8 ) : Equivalent liquid water snow depth  [m       ]
!  (10) TM_ADJ,     (REAL*8 ) : Adjustment timestep                 [s       ]
!  (11) TPT_MDP,    (REAL*8 ) : Temperature                         [K       ]
!  (12) TPT_PTN_MDP (REAL*8 ) : Midlayer local potential temp.      [K       ]
!  (13) WND_MRD_MDP (REAL*8 ) : Meridional wind component (V-wind)  [m/s     ]
!  (14) WND_ZNL_MDP (REAL*8 ) : Zonal wind component (U-wind)       [m/s     ]
!  (15) FIRST,      (LOGICAL) : Logical used ot open output dataset [unitless]
!  (16) NSTEP       (INTEGER) : Iteration counter                   [unitless]
!
!  Arguments as Output:
!  ============================================================================
!  (10) DSRC                ! O [kg kg-1] Dust mixing ratio increment
!
!  NOTES:
!  (1 ) Cleaned up and added comments.  Also force double precision with
!        "D" exponents. (bmy, 3/30/04)
!  (2 ) Now get GOCART source function. (tdf, bmy, 1/25/07)
!  (3 ) Tune nested-domain emissions dust to the same as 2x2.5 simulation
!        Also tune GEOS-3 1x1 N. America nested-grid dust emissions to
!        the 4x5 totals from the GEOS-5 4x5 v8-01-01-Run0 benchmark.
!        (yxw, bmy, dan, 11/6/08)
!  (4 ) New scale parameter for 2x2.5 GEOS-5 (tdf, jaf, phs, 10/30/09)
!  (5 ) Defined FLX_MSS_FDG_FCT for GEOS_4 2x2.5, GEOS_5 2x2.5, NESTED_NA and
!        NESTED_EU.  Redefined FLX_MSS_FDG_FCT for NESTED_CH, based upon above
!        changes. (amv, bmy, 12/18/09)
!  (6 ) For now treat MERRA like GEOS-5 (bmy, 8/13/10)
!  29 Oct 2010 - T. D. Fairlie, R. Yantosca - Retune dust for MERRA 4x5
!  08 Feb 2012 - R. Yantosca - For now, use same FLX_MSS_FDG_FCT for
!                              GEOS-5.7.x as for MERRA
!  01 Mar 2012 - R. Yantosca - Now use GET_AREA_M2(I,J,L) from grid_mod.F90
!  09 Nov 2012 - M. Payer    - Replaced all met field arrays with State_Met
!                              derived type object
!   5 Jun 2013 - K. Yu       - Use 0.5 x 0.666 NA scale factor for the
!                              0.25 x 0.3125 NA nested simulation
!******************************************************************************
!
      ! Arguments
      TYPE(HCO_State), POINTER      :: HcoState   ! Hemco state
      TYPE(Ext_State), POINTER      :: ExtState    ! Module options
      TYPE(MyInst),    POINTER      :: Inst
      INTEGER,        INTENT(IN)    :: LAT_IDX
      REAL*8,         INTENT(IN)    :: DOY
      REAL*8,         INTENT(IN)    :: HGT_MDP(HcoState%NX)
      REAL*8,         INTENT(IN)    :: LAT_RDN
      REAL*8,         INTENT(IN)    :: ORO(HcoState%NX)
      REAL*8,         INTENT(IN)    :: PRS_DLT(HcoState%NX)
      REAL*8,         INTENT(IN)    :: PRS_MDP(HcoState%NX)
      REAL*8,         INTENT(IN)    :: Q_H2O_VPR(HcoState%NX)
      REAL*8,         INTENT(IN)    :: SNW_HGT_LQD(HcoState%NX)
      REAL*8,         INTENT(IN)    :: TM_ADJ
      REAL*8,         INTENT(IN)    :: TPT_MDP(HcoState%NX)
      REAL*8,         INTENT(IN)    :: TPT_PTN_MDP(HcoState%NX)
      REAL*8,         INTENT(IN)    :: WND_MRD_MDP(HcoState%NX)
      REAL*8,         INTENT(IN)    :: WND_ZNL_MDP(HcoState%NX)
      INTEGER,        INTENT(IN)    :: NSTEP
      REAL*8,         INTENT(INOUT) :: DSRC(HcoState%NX, NBINS)
      INTEGER,        INTENT(INOUT)  :: RC

      !--------------
      ! Parameters
      !--------------

      ! Reference height for mobilization processes [m]
      REAL*8,  PARAMETER     :: HGT_RFR         = 10.0d0

      ! Zero plane displacement for erodible surfaces [m]
      REAL*8,  PARAMETER     :: HGT_ZPD_MBL     = 0.0d0

      ! Set roughness length momentum for erodible surfaces, S&P, p. 858. [m]
      REAL*8,  PARAMETER     :: RGH_MMN_MBL     = 1.0d-3

      ! rgh_mmn_smt set to 33.3e-6 um, MaB95 p. 16426 recommend 10.0e-6
      ! Smooth roughness length MaB95 p. 16426, MaB97 p. 4392, GMB98 p. 6207
      ! [m]  Z0,m,s
      REAL*8,  PARAMETER     :: RGH_MMN_SMT     = 33.3d-6

      ! Minimum windspeed used for mobilization [m/s]
      REAL*8,  PARAMETER     :: WND_MIN_MBL     = 1.0d0

      !--------------
      ! Local Output
      !--------------
      REAL*8 DST_SLT_FLX_RAT_TTL(HcoState%NX) ! [m-1] Ratio of vertical dust flux to
                                        !       streamwise mass flux
      REAL*8 FLX_MSS_HRZ_SLT_TTL(HcoState%NX) ! [kg/m/s] Vertically integrated
                                        !              streamwise mass flux
      REAL*8 FLX_MSS_VRT_DST_TTL(HcoState%NX) ! [kg/m2/s] Total vertical mass
                                        !           flux of dust
      REAL*8 FRC_THR_NCR_DRG(HcoState%NX)     ! [frc] Threshold friction velocity
                                        !       increase from roughness
      REAL*8 FRC_THR_NCR_WTR(HcoState%NX)     ! [frc] Threshold friction velocity
                                        !       increase from moisture
      REAL*8 FLX_MSS_VRT_DST(HcoState%NX,NBINS) ! [kg/m2/s] Vertical mass flux
                                            !           of dust
      REAL*8 HGT_ZPD(HcoState%NX)             ! [m] Zero plane displacement
      REAL*8 LND_FRC_MBL_SLICE(HcoState%NX)   ! [frc] Bare ground fraction
      REAL*8 MNO_LNG(HcoState%NX)             ! [m] Monin-Obukhov length
      REAL*8 WND_FRC(HcoState%NX)             ! [m/s] Friction velocity
      REAL*8 WND_FRC_GEOS(HcoState%NX)        ! [m/s] Friction velocity
      REAL*8 Z0_GEOS(HcoState%NX)             ! [m] roughness height
      REAL*8 SNW_FRC(HcoState%NX)             ! [frc] Fraction of surface covered
                                        !       by snow
      REAL*8 TRN_FSH_VPR_SOI_ATM(HcoState%NX) ! [frc] Transfer efficiency of vapor
                                        !       from soil to atmosphere
      REAL*8 wnd_frc_slt(HcoState%NX)      ! [m/s] Saltating friction velocity
      REAL*8 WND_FRC_THR_SLT(HcoState%NX)  ! [m/s] Threshold friction velocity
                                     !       for saltation
      REAL*8 WND_MDP(HcoState%NX)          ! [m/s] Surface layer mean wind speed
      REAL*8 WND_RFR(HcoState%NX)          ! [m/s] Wind speed at reference height
      REAL*8 WND_RFR_THR_SLT(HcoState%NX)  ! [m/s] Threshold 10 m wind speed for
                                     !       saltation

      LOGICAL FLG_CACO3            ! [FLG] Activate CaCO3 tracer
      LOGICAL FLG_MBL_SLICE(HcoState%NX) ! [flg] Mobilization candidates
      CHARACTER(80) FL_OUT         ! [sng] Name of netCDF output file
      INTEGER I                    ! [idx] Counting index
      INTEGER M                    ! [idx] Counting index
      INTEGER MBL_NBR              ! [nbr] Number of mobilization candidates
      INTEGER SFC_TYP_SLICE(HcoState%NX) ! [idx] LSM surface type lat slice (0..28)
      REAL*8 CND_TRM_SOI(HcoState%NX)          ! [W/m/K] Soil thermal conductivity
      REAL*8 DNS_MDP(HcoState%NX)              ! [kg/m3] Midlayer density
      REAL*8 FLX_LW_DWN_SFC_SLICE(HcoState%NX) ! [W/m2] Longwave downwelling flux
                                         !        at surface
      REAL*8 FLX_SW_ABS_SFC_SLICE(HcoState%NX) ! [W/m2] Solar flux absorbed by ground

      REAL*8 LND_FRC_DRY_SLICE(HcoState%NX)   ! [frc] Dry land fraction
      REAL*8 MBL_BSN_FCT_SLICE(HcoState%NX)   ! [frc] Erodibility factor
      REAL*8 MSS_FRC_CACO3_SLICE(HcoState%NX) ! [frc] Mass fraction of CaCO3
      REAL*8 MSS_FRC_CLY_SLICE(HcoState%NX)   ! [frc] Mass fraction of clay
      REAL*8 MSS_FRC_SND_SLICE(HcoState%NX)   ! [frc] Mass fraction of sand

      ! GOCART source function (tdf, bmy, 1/25/07)
      REAL*8 SRCE_FUNC_SLICE(HcoState%NX)     ! GOCART source function

      REAL*8 LVL_DLT(HcoState%NX) ! [m] Soil layer thickness
      REAL*8 MPL_AIR(HcoState%NX) ! [kg/m2] Air mass path in layer

      REAL*8 TM_DLT                ! [s] Mobilization timestep
      REAL*8 TPT_GND_SLICE(HcoState%NX)  ! [K] Ground temperature
      REAL*8 TPT_SOI_SLICE(HcoState%NX)  ! [K] Soil temperature
      REAL*8 TPT_SOI_FRZ           ! [K] Temperature of frozen soil
      REAL*8 TPT_VRT_MDP           ! [K] Midlayer virtual temperature
      REAL*8 VAI_DST_SLICE(HcoState%NX)  ! [m2/m2] Vegetation area index,
                                   !         one-sided
      REAL*8 VWC_DRY(HcoState%NX)        ! [m3/s] Dry volumetric water content
                                   !        (no E-T)
      REAL*8 VWC_OPT(HcoState%NX)        ! [m3/m3] E-T optimal volumetric water
                                   !         content
      REAL*8 VWC_SAT(HcoState%NX)        ! [m3/m3] Saturated volumetric water
                                   !         content (sand-dependent)
      REAL*8 VWC_SFC_SLICE(HcoState%NX)  ! [m3/m3] Volumetric water content
      REAL*8 GWC_SFC(HcoState%NX)        ! [kg/kg] Gravimetric water content
      REAL*8 RGH_MMN(HcoState%NX)        ! [m] Roughness length momentum
      REAL*8 W10M

      ! GCM diagnostics
      ! Dust tendency due to gravitational settling [kg/kg/s]
      REAL*8 Q_DST_TND_MBL(HcoState%NX,NBINS)

      ! Total dust tendency due to gravitational settling [kg/kg/s]
      REAL*8 Q_DST_TND_MBL_TTL(HcoState%NX)

      ! Temperature
      REAL(dp) :: TMP

      ! For error handling
      CHARACTER(LEN=255)   :: MSG

      !=================================================================
      ! DST_MBL begins here!
      !=================================================================

      ! Start
      RC = HCO_SUCCESS

      ! Time step [s]
      TM_DLT                 = TM_ADJ

      ! Freezing pt of soil [K] -- assume it's 0C
      TPT_SOI_FRZ            = TPT_FRZ_PNT

      ! Initialize output fluxes and tendencies
      Q_DST_TND_MBL(:,:)     = 0.0D0       ! [kg kg-1 s-1]
      Q_DST_TND_MBL_TTL(:)   = 0.0D0       ! [kg kg-1 s-1]
      FLX_MSS_VRT_DST(:,:)   = 0.0D0       ! [kg m-2 s-1]
      FLX_MSS_VRT_DST_TTL(:) = 0.0D0       ! [kg m-2 s-1]
      FRC_THR_NCR_WTR(:)     = 0.0D0       ! [frc]
      WND_RFR(:)             = 0.0D0       ! [m s-1]
      WND_FRC(:)             = 0.0D0       ! [m s-1]
      WND_FRC_SLT(:)         = 0.0D0       ! [m s-1]
      WND_FRC_THR_SLT(:)     = 0.0D0       ! [m s-1]
      WND_RFR_THR_SLT(:)     = 0.0D0       ! [m s-1]
      HGT_ZPD(:)             = HGT_ZPD_MBL ! [m]

      DSRC(:,:)              = 0.0D0

      !=================================================================
      ! Compute necessary derived fields
      !=================================================================
      DO I = 1, HcoState%NX

         ! Stop occasional haywire model runs
!         IF ( TPT_MDP(I) > 350.0d0 ) THEN
!            MSG = 'TPT_MDP(i) > 350.0'
!            CALL HCO_ERROR(HcoState%Config%Err,MSG, RC, THISLOC='DST_MBL' )
!            RETURN
!         ENDIF
         ! Now simply restrict to 350K, rather than crashing
         IF ( TPT_MDP(I) > 350.0d0 ) THEN
            TMP = 350.0d0
         ELSE
            TMP = TPT_MDP(I)
         ENDIF

         ! Midlayer virtual temperature [K]
         TPT_VRT_MDP = TMP
     &               * (1.0d0 + EPS_H2O_RCP_M1 * Q_H2O_VPR(I))

         ! Density at center of gridbox [kg/m3]
         DNS_MDP(I) = PRS_MDP(I)
     &              / (TPT_VRT_MDP * GAS_CST_DRY_AIR)

         ! Commented out
         !cApproximate surface virtual temperature (uses midlayer moisture)
         !c tpt_vrt_sfc=tpt_sfc(i)*(1.0+eps_H2O_rcp_m1*q_H2O_vpr(i)) ! [K]
         !c
         !c Surface density
         !c dns_sfc(i)=prs_sfc(i)/(tpt_vrt_sfc*gas_cst_dry_air) ! [kg m-3]

         ! Mass of air currently in gridbox [kg/m2]
         MPL_AIR(I) = PRS_DLT(I) * GRV_SFC_RCP

         ! Mean surface layer horizontal wind speed
         WND_MDP(I) = SQRT( WND_ZNL_MDP(I)*WND_ZNL_MDP(I)
     &              +       WND_MRD_MDP(I)*WND_MRD_MDP(I) )

      ENDDO

      !=================================================================
      ! Gather input variables from GEOS-CHEM modules etc.
      !=================================================================

      ! Get LSM Surface type (0..28)
      CALL SFC_TYP_GET( HcoState, ExtState, Inst,
     &                  LAT_IDX,  SFC_TYP_SLICE, RC )
      IF ( RC /= HCO_SUCCESS ) RETURN

      ! Get erodability and mass fractions
      CALL SOI_TXT_GET(
     &    HcoState,            ! Hemco state object
     &    ExtState, Inst,       ! Extension options
     &    LAT_IDX,             ! I [idx] Latitude index
     &    LND_FRC_DRY_SLICE,   ! O [frc] Dry land fraction
     &    MBL_BSN_FCT_SLICE,   ! O [frc] Erodibility factor
     &    MSS_FRC_CACO3_SLICE, ! O [frc] Mass fraction of CaCO3
     &    MSS_FRC_CLY_SLICE,   ! O [frc] Mass fraction of clay
     &    MSS_FRC_SND_SLICE )  ! O [frc] Mass fraction of sand

      ! Get GOCART source function (tdf, bmy, 1/25/07)
      CALL SRCE_FUNC_GET(      ! GOCART source function
     &    HcoState, Inst,      ! Hemco state object
     &    LAT_IDX,             ! I [idx] Latitude index
     &    SRCE_FUNC_SLICE )    ! O [frc] GOCART source function

      ! Get volumetric water content from GWET
      CALL VWC_SFC_GET(
     &    HcoState,            ! Hemco state object
     &    LAT_IDX,             ! I [idx] Latitude index
     &    ExtState%GWETTOP%Arr%Val, ! I [unitless] Top soil moisture
     &    VWC_SFC_SLICE )      ! O [m3 m-3] Volumetric water content

      ! Get surface and soil temperature
      CALL TPT_GND_SOI_GET(
     &    HcoState,              ! Hemco state object
     &    LAT_IDX,               ! I [idx] Latitude index!
     &    ExtState%T2M%Arr%Val, ! I [K] Sfc temperature at 2m
     &    TPT_GND_SLICE,         ! O [K] Ground temperature
     &    TPT_SOI_SLICE )        ! O [K] Soil temperature

      ! Get time-varying vegetation area index
      CALL DST_TVBDS_GET(
     &    Inst,                ! # of lons
     &    HcoState%NX,         ! # of lons
     &    LAT_IDX,             ! I [idx] Latitude index
     &    VAI_DST_SLICE)       ! O [m2 m-2] Vegetation area index, one-sided

      ! Get fraction of surface covered by snow
      CALL SNW_FRC_GET(
     &    HcoState,            ! Hemco state object
     &    SNW_HGT_LQD,         ! I [m] Equivalent liquid water snow depth
     &    SNW_FRC )            ! O [frc] Fraction of surface covered by snow

      !=================================================================
      ! Use the variables retrieved above to compute the fraction
      ! of each gridcell suitable for dust mobilization
      !=================================================================
      CALL LND_FRC_MBL_GET(
     %    HcoState,
     &    DOY,                 ! I [day] Day of year [1.0..366.0)
     &    FLG_MBL_SLICE,       ! O [flg] Mobilization candidate flag
     &    LAT_RDN,             ! I [rdn] Latitude
     &    LND_FRC_DRY_SLICE,   ! I [frc] Dry land fraction
     &    LND_FRC_MBL_SLICE,   ! O [frc] Bare ground fraction
     &    MBL_NBR,             ! O [flg] Number of mobilization candidates
     &    ORO,                 ! I [frc] Orography
     &    SFC_TYP_SLICE,       ! I [idx] LSM surface type (0..28)
     &    SNW_FRC,             ! I [frc] Fraction of surface covered by snow
     &    TPT_SOI_SLICE,       ! I [K] Soil temperature
     &    TPT_SOI_FRZ,         ! I [K] Temperature of frozen soil
     &    VAI_DST_SLICE,       ! I [m2 m-2] Vegetation area index, one-sided
     &    RC )
      IF ( RC /= HCO_SUCCESS ) RETURN

      ! Much ado about nothing
      if (mbl_nbr == 0) then
        goto 737
      endif

      !=================================================================
      ! Compute time-invariant hydrologic properties
      ! NB flg_mbl IS time-dependent, so keep this in time loop.
      !=================================================================
      CALL HYD_PRP_GET(        ! NB: These properties are time-invariant
     &    HcoState,
     &    FLG_MBL_SLICE,       ! I [flg] Mobilization candidate flag
     &    MSS_FRC_CLY_SLICE,   ! I [frc] Mass fraction clay
     &    MSS_FRC_SND_SLICE,   ! I [frc] Mass fraction sand
     &    VWC_DRY,             ! O [m3/m3] Dry vol'mtric water content (no E-T)
     &    VWC_OPT,             ! O [m3/m3] E-T optimal volumetric water content
     &    VWC_SAT)             ! O [m3/m3] Saturated volumetric water content

      CND_TRM_SOI(:) = 0.0D0
      LVL_DLT(:)     = 0.0D0

      !=================================================================
      ! Get reference wind at 10m
      !=================================================================
      DO I = 1, HcoState%NX
         W10M = ExtState%U10M%Arr%Val(I,LAT_IDX)**2 +
     &          ExtState%V10M%Arr%Val(I,LAT_IDX)**2
         W10M = SQRT(W10M)

         ! add mobilisation criterion flag
         IF ( FLG_MBL_SLICE(I) ) THEN
            WND_RFR(I) = W10M
         ENDIF
      ENDDO

      !=================================================================
      ! Compute standard roughness length.   This call is probably
      ! unnecessary, because we are only concerned with mobilisation
      ! candidates, for which roughness length is imposed in blm_mbl
      !=================================================================
      CALL RGH_MMN_GET(      ! Set roughness length w/o zero plane displacement
     &       HcoState, Inst,
     &       ORO,            ! I [frc] Orography
     &       RGH_MMN,        ! O [m] Roughness length momentum
     &       SFC_TYP_SLICE,  ! I [idx] LSM surface type (0..28)
     &       SNW_FRC,        ! I [frc] Fraction of surface covered by snow
     &       WND_RFR,
     &       RC )       ! I [m s-1] 10 m wind speed
      IF ( RC /= HCO_SUCCESS ) RETURN

      !=================================================================
      ! Introduce Ustar and Z0 from GEOS data
      !=================================================================
      DO I = 1, HcoState%NX

         ! Just assign for flag mobilisation candidates
         IF ( FLG_MBL_SLICE(I) ) THEN
            WND_FRC_GEOS(I) = ExtState%USTAR%Arr%Val(I,LAT_IDX)
            Z0_GEOS(I)      = ExtState%Z0%Arr%Val(I,LAT_IDX)
         ELSE
            WND_FRC_GEOS(I) = 0.0D0
            Z0_GEOS(I)      = 0.0D0
         ENDIF
      ENDDO

      !=================================================================
      ! Surface exchange properties over erodible surfaces
      ! DO NEED THIS: Compute Monin-Obukhov and Friction velocities
      ! appropriate for dust producing regions.
      !
      ! Now calling Stripped down (adiabatic) version     tdf 10/27/2K3
      ! rgh_mmn_mbl parameter included directly in blm_mbl
      !=================================================================
      CALL BLM_MBL(
     &    HcoState,
     &    FLG_MBL_SLICE,       ! I [flg] Mobilization candidate flag
     &    RGH_MMN,             ! I [m] Roughness length momentum, Z0,m
     &    WND_RFR,             ! I [m s-1] 10 m wind speed
     &    MNO_LNG,             ! O [m] Monin-Obukhov length
     &    WND_FRC,
     &    RC )                 ! O [m s-1] Surface friction velocity, U*
      IF ( RC /= HCO_SUCCESS ) RETURN

      !=================================================================
      ! Factor by which surface roughness increases threshold friction
      ! velocity.  The sink of atrmospheric momentum into non-erodible
      ! roughness elements Zender et al., expression (3)
      !=================================================================
!-----------------------------------------------------------------------------
! Prior to 1/25/07:
! For now, instead of calling this routine to get FRC_THR_NCR_DRG, we will
! just set it to 1 (tdf, bmy, 1/25/07)
!
! %%%%% DO NOT DELETE -- LEAVE THIS CODE COMMENTED OUT %%%%%
!
!      CALL FRC_THR_NCR_DRG_GET(
!     &    HcoState,
!     &    FRC_THR_NCR_DRG,     ! O [frc] Factor increases thresh. fric. veloc.
!     &    FLG_MBL_SLICE,       ! I [flg] Mobilization candidate flag
!     &    RGH_MMN_MBL,         ! I [m] Rgh length momentum for erodible sfcs
!     &    RGH_MMN_SMT,         ! I [m] Smooth roughness length, Z0,m,s
!     &    RC )
!-----------------------------------------------------------------------------

      ! Now set roughness factor to 1.0 (tdf, bmy, 1/25/07)
      FRC_THR_NCR_DRG(:) = 1.0d0

      !=================================================================
      ! Convert volumetric water content to gravimetric water content
      ! NB: Owen effect included in wnd_frc_slt_get
      !=================================================================
      CALL VWC2GWC(
     &    HcoState,
     &    FLG_MBL_SLICE,       ! I [flg] Mobilization candidate flag
     &    GWC_SFC,             ! O [kg kg-1] Gravimetric water content
     &    VWC_SAT,             ! I [m3 m-3] Saturated VWC (sand-dependent)
     &    VWC_SFC_SLICE )      ! I [m3 m-3] Volumetric water content

      !=================================================================
      ! Factor by which soil moisture increases threshold friction
      ! velocity -- i.e. the inhibition of saltation by soil mositure,
      ! Zender et al., exp(5).
      !=================================================================
      CALL FRC_THR_NCR_WTR_GET(
     &    HcoState,
     &    FLG_MBL_SLICE,       ! I [flg] Mobilization candidate flag
     &    FRC_THR_NCR_WTR,     ! O [frc] Factor by which moisture increases
                               !         threshold friction velocity
     &    MSS_FRC_CLY_SLICE,   ! I [frc] Mass fraction of clay
     &    GWC_SFC)             ! I [kg kg-1] Gravimetric water content

      !=================================================================
      ! Now, compute basic threshold friction velocity for saltation
      ! over dry, bare, smooth ground.  fxm: Use surface density not
      ! midlayer density
      !=================================================================
      CALL WND_FRC_THR_SLT_GET(
     &    HcoState,
     &    FLG_MBL_SLICE,       ! I mobilisation flag
     &    DNS_MDP,             ! I [kg m-3] Midlayer density
     &    WND_FRC_THR_SLT,     ! O [m s-1] Threshold friction velocity
     &    RC )
      IF ( RC /= HCO_SUCCESS ) RETURN

      ! Adjust threshold friction velocity to account
      ! for moisture and roughness
      DO I = 1, HcoState%NX
         WND_FRC_THR_SLT(I) =      ! [m s-1] Threshold friction velocity
                                   !         for saltation
     &        WND_FRC_THR_SLT(i)   ! [m s-1] Threshold for dry, flat ground
     &        * FRC_THR_NCR_WTR(i) ! [frc] Adjustment for moisture
     &        * FRC_THR_NCR_DRG(i) ! [frc] Adjustment for roughness
      ENDDO

      ! Threshold saltation wind speed at reference height, 10m
      DO I = 1, HcoState%NX
         IF ( FLG_MBL_SLICE(I) ) THEN
           WND_RFR_THR_SLT(I) =  ! [m s-1] Threshold 10 m wind speed
                                 !         for saltation
     &     WND_RFR(I) * WND_FRC_THR_SLT(I) / WND_FRC(i)
         ENDIF
      ENDDO

      !=================================================================
      ! Saltation increases friction speed by roughening surface
      ! i.e. Owen effect, Zender et al., expression (4)
      !
      ! Compute the wind friction velocity due to saltation, U*,s
      ! accounting for the Owen effect.
      !=================================================================
      CALL WND_FRC_SLT_GET(
     &    HcoState,
     &    FLG_MBL_SLICE,     ! I [flg] Mobilization candidate flag
     &    WND_FRC,           ! I [m s-1] Surface friction velocity
     &    WND_FRC_SLT,       ! O [m s-1] Saltating friction velocity
     &    WND_RFR,           ! I [m s-1] Wind speed at reference height
     &    WND_RFR_THR_SLT )  ! I [m s-1] Thresh. 10 m wind speed for saltation

      !=================================================================
      ! Compute horizontal streamwise mass flux, Zender et al., expr. (10)
      !=================================================================
      CALL FLX_MSS_HRZ_SLT_TTL_WHI79_GET(
     &    HcoState,
     &    DNS_MDP,             ! I [kg m-3] Midlayer density
     &    FLG_MBL_SLICE,       ! I [flg] Mobilization candidate flag
     &    FLX_MSS_HRZ_SLT_TTL, ! O [kg m-1 s-1] Vertically integrated
                               !                streamwise mass flux
     &    WND_FRC_SLT,         ! I [m s-1] Saltating friction velocity
     &    WND_FRC_THR_SLT )    ! I [m s-1] Threshold friction vel for saltation

!-----------------------------------------------------------------------------
! Prior to 1/25/07:
! We now multiply by the GOCART source function, and we will ignore
! the MBL_BSN_FCT_SLICE.  (tdf, bmy, 1/25/07)
!
! %%%%% DO NOT DELETE -- LEAVE THIS CODE COMMENTED OUT %%%%%
!
!ctdf...prior to Apr/05/06
!      ! Apply land surface and vegetation limitations
!      ! and global tuning factor
!      DO I = 1, HcoState%NX
!         FLX_MSS_HRZ_SLT_TTL(I) = FLX_MSS_HRZ_SLT_TTL(I) ! [kg m-2 s-1]
!     &       * LND_FRC_MBL_SLICE(i)   ! [frc] Bare ground fraction
!     &       * MBL_BSN_FCT_SLICE(i)   ! [frc] Erodibility factor
!     &       * FLX_MSS_FDG_FCT        ! [frc] Global mass flux tuning
!                                      !       factor (empirical)
!      ENDDO
!-----------------------------------------------------------------------------

      ! Now simply multiply by the GOCART source function.
      ! The vegetation effect has been eliminated in LND_FRC_MBL_GET
      ! and we also ignore MBL_BSN_FCT. (tdf, bmy, 1/25/07)
      DO I = 1, HcoState%NX
         FLX_MSS_HRZ_SLT_TTL(I) = FLX_MSS_HRZ_SLT_TTL(I) ! [kg m-2 s-1]
     &       * LND_FRC_MBL_SLICE(i)   ! [frc] Bare ground fraction
     &       * Inst%FLX_MSS_FDG_FCT   ! [frc] Global mass flux tuning
     &       * SRCE_FUNC_SLICE(I)     ! GOCART source function
      ENDDO

      !=================================================================
      ! Compute vertical dust mass flux, see Zender et al., expr. (11).
      !=================================================================
      CALL FLX_MSS_VRT_DST_TTL_MAB95_GET(
     &    HcoState,
     &    DST_SLT_FLX_RAT_TTL, ! O [m-1] Ratio of vertical dust flux to
                               !         streamwise mass flux
     &    FLG_MBL_SLICE,       ! I [flg] Mobilization candidate flag
     &    FLX_MSS_HRZ_SLT_TTL, ! I [kg/m/s] Vertically integrated
                               !            streamwise mass flux
     &    FLX_MSS_VRT_DST_TTL, ! O [kg/m2/s] Total vertical mass flux of dust
     &    MSS_FRC_CLY_SLICE )  ! I [frc] Mass fraction clay

      !=================================================================
      ! Now, partition vertical dust mass flux into transport bins
      !
      ! OVR_SRC_SNK_MSS needed in FLX_MSS_VRT_DST_PRT
      ! computed in DST_PSD_MSS, called from "dust_mod.f" (tdf, 3/30/04)
      !=================================================================
      CALL FLX_MSS_VRT_DST_PRT( Inst,
     &    HcoState%NX,
     &    FLG_MBL_SLICE,       ! I [flg] Mobilization candidate flag
     &    FLX_MSS_VRT_DST,     ! O [kg m-2 s-1] Vertical mass flux of dust
     &    FLX_MSS_VRT_DST_TTL) ! I [kg m-2 s-1] Total vertical mass flux of dus

      !=================================================================
      ! Mask dust mass flux by tracer mass fraction at source
      !=================================================================
      FLG_CACO3 = .FALSE.          ! [flg] Activate CaCO3 tracer
      IF ( FLG_CACO3 ) THEN
         CALL FLX_MSS_CACO3_MSK(
     &        HcoState,
     &        ExtState,
     &        Inst%DMT_VWR,        ! I [m] Mass weighted diameter resolved
     &        FLG_MBL_SLICE,       ! I [flg] Mobilization candidate flag
     &        FLX_MSS_VRT_DST,     ! I/O [kg m-2 s-1] Vert. mass flux of dust
     &        MSS_FRC_CACO3_SLICE, ! I [frc] Mass fraction of CaCO3
     &        MSS_FRC_CLY_SLICE,   ! I [frc] Mass fraction of clay
     &        MSS_FRC_SND_SLICE,   ! I [frc] Mass fraction of sand
     &        RC )
         IF ( RC /= HCO_SUCCESS ) RETURN
      endif

      ! Now, flx_mss_vrt_dst has units of kg/m2/sec

      ! Fluxes are known, so adjust mixing ratios
      DO  I=1, HcoState%NX            ! NB: Inefficient loop order
         IF (FLG_MBL_SLICE(I)) THEN

            ! Loop over dust bins
            DO M = 1, NBINS

               !========================================================
               ! Compute dust mobilisation tendency.  Recognise that
               ! what GEOS-CHEM wants is an increment in kg...So,
               ! multiply by DXYP [m2] and tm_adj [sec]
               !========================================================

               ! [kg/sec]
               Q_DST_TND_MBL(I,M) = FLX_MSS_VRT_DST(I,M)
     &            *HcoState%Grid%AREA_M2%Val(I,LAT_IDX)

               ! Introduce DSRC: dust mixing ratio increment   12/9/2K3
               ! [kg]
               DSRC(I,M) = TM_ADJ * Q_DST_TND_MBL(I,M)

           ENDDO
         ENDIF
      ENDDO

      ! Jump to here when no points are mobilization candidates
  737 CONTINUE

      RC = HCO_SUCCESS

      ! Return to calling program
      END SUBROUTINE DST_MBL

!------------------------------------------------------------------------------

      SUBROUTINE SRCE_FUNC_GET( HcoState, Inst, LAT_IDX, SRCE_FUNC_OUT )
!
!******************************************************************************
!  Subroutine SRCE_FUNC_GET returns a latitude slice of the GOCART source
!  function.  This routine is called by DST_MBL. (tdf, bmy, 1/25/07)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) LAT_IDX       (INTEGER) : GEOS-Chem latitude index
!
!  Arguments as Output:
!  ============================================================================
!  (1 ) SRCE_FUNC_OUT (REAL*8 ) : GOCART source function [fraction]
!
!  NOTES:
!******************************************************************************
!
      ! Arguments
      TYPE(HCO_State), POINTER :: HcoState   ! Hemco state
      TYPE(MyInst),    POINTER :: Inst
      INTEGER, INTENT(IN)      :: LAT_IDX
      REAL*8,  INTENT(OUT)     :: SRCE_FUNC_OUT(HcoState%NX)

      ! Local variables
      INTEGER              :: LON_IDX

      !=================================================================
      ! SRCE_FUNC_GET begins here!
      !=================================================================

      ! Loop over longitudes
      DO LON_IDX = 1, HcoState%NX

         ! Save latitude slice in SRCE_FUNC_OUT
         SRCE_FUNC_OUT(LON_IDX) = Inst%SRCE_FUNC(LON_IDX,LAT_IDX)

      ENDDO

      ! Return to calling program
      END SUBROUTINE SRCE_FUNC_GET

!------------------------------------------------------------------------------

      SUBROUTINE SOI_TXT_GET( HcoState, ExtState, Inst, J,
     &                        LND_FRC_DRY_OUT,
     &                        MBL_BSN_FCT_OUT, MSS_FRC_CACO3_OUT,
     &                        MSS_FRC_CLY_OUT, MSS_FRC_SND_OUT )
!
!******************************************************************************
!  Subroutine SOI_GET_TXT returns a latitude slice of soil texture to the
!  calling program DST_MBL.  (tdf, bmy, 3/30/04)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) J                 (INTEGER) : Grid box latitude index
!
!  Arguments as Output:
!  ============================================================================
!  (2 ) lnd_frc_dry_out   (REAL*8 ) : Dry land fraction      [fraction]
!  (3 ) mbl_bsn_fct_out   (REAL*8 ) : Erodibility factor     [fraction]
!  (4 ) mss_frc_CaCO3_out (REAL*8 ) : Mass fraction of CaCO3 [fraction]
!  (5 ) mss_frc_cly_out   (REAL*8 ) : Mass fraction of clay  [fraction]
!  (6 ) mss_frc_snd_out   (REAL*8 ) : Mass fraction of sand  [fraction]
!
!  NOTES:
!  (1 ) Updated comments, cosmetic changes (bmy, 3/30/04)
!******************************************************************************
!

      ! Arguments
      TYPE(HCO_State), POINTER :: HcoState   ! Hemco state
      TYPE(Ext_State), POINTER :: ExtState    ! Module options
      TYPE(MyInst),    POINTER :: Inst
      INTEGER, INTENT(IN)  :: J
      REAL*8,  INTENT(OUT) :: LND_FRC_DRY_OUT(HcoState%NX)
      REAL*8,  INTENT(OUT) :: MBL_BSN_FCT_OUT(HcoState%NX)
      REAL*8,  INTENT(OUT) :: MSS_FRC_CACO3_OUT(HcoState%NX)
      REAL*8,  INTENT(OUT) :: MSS_FRC_CLY_OUT(HcoState%NX)
      REAL*8,  INTENT(OUT) :: MSS_FRC_SND_OUT(HcoState%NX)

      ! Local variables
      INTEGER              :: I

      ! Ad hoc globally uniform clay mass fraction [kg/kg]
      REAL*8,  PARAMETER   :: MSS_FRC_CLY_GLB = 0.20d0

      !=================================================================
      ! SOI_GET_TXT begins here!
      !=================================================================
      DO I = 1, HcoState%NX

         ! Save dry land fraction slice
         LND_FRC_DRY_OUT(I) = Inst%LND_FRC_DRY(I,J)

         ! Change surface source distribution to "geomorphic"  tdf 12/12/2K3
         MBL_BSN_FCT_OUT(I) = Inst%ERD_FCT_GEO(I,J)

         !fxm: CaCO3 currently has missing value of
         !     1.0e36 which causes problems
         IF ( Inst%MSS_FRC_CACO3(I,J) <= 1.0D0 ) THEN
            MSS_FRC_CACO3_OUT(I) = Inst%MSS_FRC_CACO3(I,J)
         ELSE
            MSS_FRC_CACO3_OUT(I) = 0.0D0
         ENDIF

         ! fxm Temporarily set mss_frc_cly used in mobilization to globally
         !     uniform SGS value of 0.20, and put excess mass fraction
         !     into sand
         MSS_FRC_CLY_OUT(I) = MSS_FRC_CLY_GLB
         MSS_FRC_SND_OUT(I) = Inst%MSS_FRC_SND(I,J) +
     &                        Inst%MSS_FRC_CLY(I,J) -
     &                        MSS_FRC_CLY_GLB

      ENDDO

      ! Return to calling program
      END SUBROUTINE SOI_TXT_GET

!------------------------------------------------------------------------------

      SUBROUTINE SFC_TYP_GET( HcoState, ExtState,
     &                        Inst, J, SFC_TYP_OUT, RC )
!
!******************************************************************************
!  Subroutine SFC_TYP_GET returns a latitude slice of LSM surface type
!  to the calling programs DST_MBL & DST_DPS_DRY. (tdf, bmy, 3/30/04)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) J           (INTEGER) : Grid box latitude index
!
!  Arguments as Output:
!  ============================================================================
!  (1 ) sfc_typ_out (REAL*8 ) : LSM surface type (0..28) [unitless]
!
!  NOTES
!  (1 ) Updated comments & cosmetic changes (bmy, 3/30/04)
!  (2 ) Added error trap (ckeller, 7/24/2014)
!******************************************************************************
!

      ! Arguments
      TYPE(HCO_State), POINTER :: HcoState   ! Hemco state
      TYPE(Ext_State), POINTER :: ExtState
      TYPE(MyInst),    POINTER :: Inst
      INTEGER, INTENT(IN)      :: J
      INTEGER, INTENT(OUT)     :: SFC_TYP_OUT(HcoState%NX)
      INTEGER, INTENT(INOUT)   :: RC

      ! Local variables
      INTEGER                  :: I, TMP
      CHARACTER(LEN=255)       :: MSG

      !=================================================================
      ! SFC_TYP_GET begins here!
      !=================================================================
      DO I = 1, HcoState%NX
         TMP = INT(Inst%SFC_TYP(I,J))

         ! Make sure value is within valid range (1 - NN_SFCTYP).
         SFC_TYP_OUT(I) = MAX( MIN(TMP,NN_SFCTYP), 0 )
      ENDDO

      ! Return with success
      RC = HCO_SUCCESS

      ! Return to calling program
      END SUBROUTINE SFC_TYP_GET                       ! end sfc_typ_get()

!------------------------------------------------------------------------------

      SUBROUTINE TPT_GND_SOI_GET( HcoState, J, TS,
     &                            TPT_GND_OUT, TPT_SOI_OUT )
!
!******************************************************************************
!  Subroutine TPT_GND_SOI_GET returns a latitude slice of soil temperature and
!  ground temperature to the calling program DST_MBL. (tdf, bmy, 3/30/04)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) J           (INTEGER) : Grid box latitude index
!  (2 ) TS          (REAL*8)  : Surface temperature at 2m [K]
!
!  Arguments as Output:
!  ============================================================================
!  (2 ) TPT_GND_OUT (REAL*8 ) : Ground temperature array slice [K]
!  (3 ) tpt_soi_out (REAL*8 ) : Soil temperature array slice   [K]
!
!  NOTES
!  (1 ) Updated comments & cosmetic changes (bmy, 3/30/04)
!******************************************************************************
!

      ! Arguments
      INTEGER, INTENT(IN)  :: J
      TYPE(HCO_State), POINTER :: HcoState   ! Hemco state
      REAL(hp),INTENT(IN)  :: TS(HcoState%NX,HcoState%NY)
      REAL*8,  INTENT(OUT) :: TPT_GND_OUT(HcoState%NX)
      REAL*8,  INTENT(OUT) :: TPT_SOI_OUT(HcoState%NX)

      ! Local variables
      INTEGER              :: I

      !=================================================================
      ! TPT_GND_SOI_GET begins here!
      !=================================================================

      ! Use TS from GEOS-CHEM (tdf, 3/30/04)
      DO I = 1, HcoState%NX
         TPT_GND_OUT(I) = TS(I,J)
         TPT_SOI_OUT(I) = TS(I,J)
      ENDDO

      ! Return to calling program
      END SUBROUTINE TPT_GND_SOI_GET

!------------------------------------------------------------------------------

      SUBROUTINE VWC_SFC_GET( HcoState, J, GWETTOP, VWC_SFC_OUT )
!
!******************************************************************************
!  Subroutine TPT_GND_SOI_GET returns a latitude slice of volumetric water
!  content to the calling program DST_MBL. (tdf, bmy, 3/30/04)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) J       (INTEGER) : Grid box latitude index
!  (2 ) GWETTOP (REAL*8)  : Top soil moisture [unitless]
!
!  Arguments as Output:
!  ============================================================================
!  VWC_SFC_OUT  (REAL*8 ) : Volumetric water content [m3/m3]
!
!  NOTES
!  (1 ) Updated comments & cosmetic changes (bmy, 3/30/04)
!******************************************************************************
!

      ! Arguments
      INTEGER, INTENT(IN)  :: J
      TYPE(HCO_State), POINTER        :: HcoState   ! Hemco state
      REAL(hp), INTENT(IN)  :: GWETTOP(HcoState%NX,HcoState%NY)
      REAL*8,  INTENT(OUT) :: VWC_SFC_OUT(HcoState%NX)

      ! Local variables
      INTEGER              :: I

      !=================================================================
      ! VWC_SFC_GET begins here!
      !=================================================================
      DO I = 1, HcoState%NX
         VWC_SFC_OUT(I) = GWETTOP(I,J)
      ENDDO

      ! Return to calling program
      END SUBROUTINE VWC_SFC_GET

!------------------------------------------------------------------------------

      REAL*8 FUNCTION DSVPDT_H2O_LQD_PRK78_FST_SCL( TPT_CLS )
!
!******************************************************************************
!  Function DSVPDT_H2O_LQD_PRK78_FST_SCL returns the derivative of saturation
!  vapor pressure [Pa] over planar liquid water (tdf, bmy, 3/30/04)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) TPT_CLS (REAL*8) : Temperature in Celsius [C]
!
!  NOTES:
!  (1 ) Updated comments, cosmetic changes.  Also now force double-precision
!        with "D" exponents. (bmy, 3/30/04)
!******************************************************************************
!
      ! Arguments
      REAL*8, INTENT(IN) :: TPT_CLS

      ! Local variables
      REAL*8, PARAMETER  :: C0 = 4.438099984d-01
      REAL*8, PARAMETER  :: C1 = 2.857002636d-02
      REAL*8, PARAMETER  :: C2 = 7.938054040d-04
      REAL*8, PARAMETER  :: C3 = 1.215215065d-05
      REAL*8, PARAMETER  :: C4 = 1.036561403d-07
      REAL*8, PARAMETER  :: C5 = 3.532421810d-10
      REAL*8, PARAMETER  :: C6 =-7.090244804d-13

      !=================================================================
      ! DSVPDT_H2O_LQD_PRK78_FST_SCL begins here!
      !=================================================================

      ! Return deriv. of saturation vapor pressure [Pa]
      DSVPDT_H2O_LQD_PRK78_FST_SCL = 100.0d0 * ( C0+TPT_CLS *
     &                                         ( C1+TPT_CLS *
     &                                         ( C2+TPT_CLS *
     &                                         ( C3+TPT_CLS *
     &                                         ( C4+TPT_CLS *
     &                                         ( C5+TPT_CLS * C6 ))))))

      ! Return to calling program
      END FUNCTION DSVPDT_H2O_LQD_PRK78_FST_SCL

!------------------------------------------------------------------------------

      REAL*8 FUNCTION DSVPDT_H2O_ICE_PRK78_FST_SCL( TPT_CLS )
!
!******************************************************************************
!  Function DSVPDT_H2O_ICE_PRK78_FST_SCL returns the derivative of saturation
!  vapor pressure [Pa] over planar ice water (tdf, bmy, 3/30/04)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) TPT_CLS (REAL*8) : Temperature in Celsius [C]
!
!  NOTES:
!  (1 ) Updated comments, cosmetic changes.  Also now force double-precision
!        with "D" exponents. (bmy, 3/30/04)
!******************************************************************************
!
      ! Arguments
      REAL*8, INTENT(IN) :: TPT_CLS

      ! Local variables
      REAL*8, PARAMETER  :: D0 = 5.030305237d-01
      REAL*8, PARAMETER  :: D1 = 3.773255020d-02
      REAL*8, PARAMETER  :: D2 = 1.267995369d-03
      REAL*8, PARAMETER  :: D3 = 2.477563108d-05
      REAL*8, PARAMETER  :: D4 = 3.005693132d-07
      REAL*8, PARAMETER  :: D5 = 2.158542548d-09
      REAL*8, PARAMETER  :: D6 = 7.131097725d-12

      !=================================================================
      ! DSVPDT_H2O_ICE_PRK78_FST_SCL begins here!
      !=================================================================

      ! Return deriv. of sat vapor pressure [Pa]
      DSVPDT_H2O_ICE_PRK78_FST_SCL = 100.0D0 * ( D0+TPT_CLS *
     &                                         ( D1+TPT_CLS *
     &                                         ( D2+TPT_CLS *
     &                                         ( D3+TPT_CLS *
     &                                         ( D4+TPT_CLS *
     &                                         ( D5+TPT_CLS * D6 ))))))

      ! Return to calling program
      END FUNCTION DSVPDT_H2O_ICE_PRK78_FST_SCL

!------------------------------------------------------------------------------

      REAL*8 FUNCTION SVP_H2O_LQD_PRK78_FST_SCL( TPT_CLS )
!
!******************************************************************************
!  Function SVP_H2O_LQD_PRK78_FST_SCL returns the saturation vapor pressure
!  over planer liquid water [Pa]  See Lowe and Ficke (1974) as reported in
!  PrK78 p. 625. Range of validity is -50 C < T < 50 C. (tdf, bmy, 3/30/04)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) TPT_CLS (REAL*8) : Temperature in Celsius [C]
!
!  NOTES:
!  (1 ) Updated comments, cosmetic changes.  Also now force double-precision
!        with "D" exponents. (bmy, 3/30/04)
!******************************************************************************
!
      ! Arguments
      REAL*8, INTENT(IN) :: TPT_CLS

      ! Local variables
      REAL*8, PARAMETER  :: A0 = 6.107799961d0
      REAL*8, PARAMETER  :: A1 = 4.436518521d-01
      REAL*8, PARAMETER  :: A2 = 1.428945805d-02
      REAL*8, PARAMETER  :: A3 = 2.650648471d-04
      REAL*8, PARAMETER  :: A4 = 3.031240396d-06
      REAL*8, PARAMETER  :: A5 = 2.034080948d-08
      REAL*8, PARAMETER  :: A6 = 6.136820929d-11

      !=================================================================
      ! SVP_H2O_LQD_PRK78_FST_SCL begins here!
      !=================================================================

      ! Return saturation vapor pressure over liquid water [Pa]
      SVP_H2O_LQD_PRK78_FST_SCL = 100.0D0 * ( A0+TPT_CLS *
     &                                      ( A1+TPT_CLS *
     &                                      ( A2+TPT_CLS *
     &                                      ( A3+TPT_CLS *
     &                                      ( A4+TPT_CLS *
     &                                      ( A5+TPT_CLS * A6 ))))))

      ! Return to calling program
      END FUNCTION SVP_H2O_LQD_PRK78_FST_SCL

!------------------------------------------------------------------------------

      REAL*8 FUNCTION SVP_H2O_ICE_PRK78_FST_SCL( TPT_CLS )
!
!******************************************************************************
!  Function SVP_H2O_ICE_PRK78_FST_SCL returns the saturation vapor pressure
!  [Pa] over planar ice water (tdf, bmy, 3/30/04)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) TPT_CLS (REAL*8) : Temperature in Celsius [C]
!
!  NOTES:
!  (1 ) Updated comments, cosmetic changes.  Also now force double-precision
!        with "D" exponents. (bmy, 3/30/04)
!******************************************************************************
!

      ! Arguments
      REAL*8, INTENT(IN) :: TPT_CLS

      ! Local variables
      REAL*8, PARAMETER  :: B0 = 6.109177956d0
      REAL*8, PARAMETER  :: B1 = 5.034698970d-01
      REAL*8, PARAMETER  :: B2 = 1.886013408d-02
      REAL*8, PARAMETER  :: B3 = 4.176223716d-04
      REAL*8, PARAMETER  :: B4 = 5.824720280d-06
      REAL*8, PARAMETER  :: B5 = 4.838803174d-08
      REAL*8, PARAMETER  :: B6 = 1.838826904d-10

      !=================================================================
      ! SVP_H2O_ICE_PRK78_FST_SCL begins here!
      !=================================================================

      ! Return saturation vapor pressure over ice [Pa]
      SVP_H2O_ICE_PRK78_FST_SCL = 100.0D0 * ( B0+TPT_CLS *
     &                                      ( B1+TPT_CLS *
     &                                      ( B2+TPT_CLS *
     &                                      ( B3+TPT_CLS *
     &                                      ( B4+TPT_CLS *
     &                                      ( B5+TPT_CLS * B6 ))))))

      ! Return to calling program
      END FUNCTION SVP_H2O_ICE_PRK78_FST_SCL

!------------------------------------------------------------------------------

      REAL*8 FUNCTION TPT_BND_CLS_GET( TPT )
!
!******************************************************************************
!  Function TPT_BND_CLS_GET returns the bounded temperature in [C],
!  (i.e., -50 < T [C] < 50 C), given the temperature in [K].
!  (tdf, bmy, 3/30/04)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) TPT (REAL*8) : Temperature in Kelvin [K]
!
!  NOTES:
!******************************************************************************
!
      ! Arguments
      REAL*8, INTENT(IN) :: TPT

      ! Local variables
      REAL*8, PARAMETER  :: TPT_FRZ_PNT=273.15

      !=================================================================
      ! TPT_BND_CLS_GET begins here!
      !=================================================================
      TPT_BND_CLS_GET = MIN( 50.0D0, MAX( -50.0D0, ( TPT-TPT_FRZ_PNT)) )

      ! Return to calling program
      END FUNCTION TPT_BND_CLS_GET

!------------------------------------------------------------------------------

      SUBROUTINE GET_ORO( HcoState, ExtState, OROGRAPHY, RC )
!
! !USES:
!
      USE HCO_GEOTOOLS_MOD, ONLY : HCO_LANDTYPE
!
!******************************************************************************
!  Subroutine GET_ORO creates a 2D orography array, OROGRAPHY, from the
!  GMAO LWI fields.  Ocean= 0; Land=1; ice=2. (tdf, bmy, 3/30/04, 8/17/05)
!
!  Arguments as Output:
!  ============================================================================
!  (1 ) OROGRAPHY (INTEGER) : Array for orography flags
!
!  NOTES:
!  (1 ) Added parallel DO-loop (bmy, 4/14/04)
!  (2 ) Now modified for GCAP and GEOS-5 met fields (swu, bmy, 6/9/05)
!  (3 ) Now use IS_LAND, IS_WATER, IS_ICE functions from "dao_mod.f"
!        (bmy, 8/17/05)
!  09 Nov 2012 - M. Payer    - Replaced all met field arrays with State_Met
!                              derived type object
!******************************************************************************
!

      ! Arguments
      TYPE(HCO_State), POINTER       :: HcoState
      Type(Ext_State), POINTER       :: ExtState
      INTEGER,         INTENT(OUT  ) :: OROGRAPHY(HcoState%NX,
     &                                            HcoState%NY)
      INTEGER,         INTENT(INOUT) :: RC

      ! Local variables
      INTEGER             :: I, J

      !=================================================================
      ! GET_ORO begins here!
      !=================================================================

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J )
      DO J = 1, HcoState%NY
      DO I = 1, HcoState%NX

         ! Set orography to land type
         OROGRAPHY (I,J) = HCO_LANDTYPE( ExtState%WLI%Arr%Val(I,J),
     &                                   ExtState%ALBD%Arr%Val(I,J) )

      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      ! Return w/ success
      RC = HCO_SUCCESS

      END SUBROUTINE GET_ORO

!------------------------------------------------------------------------------

      SUBROUTINE HYD_PRP_GET( HcoState,     FLG_MBL, MSS_FRC_CLY_SLC,
     &                        MSS_FRC_SND_SLC, VWC_DRY, VWC_OPT,
     &                        VWC_SAT                            )
!
!******************************************************************************
!  Subroutine HYD_PRP_GET determines hydrologic properties from soil texture.
!  (tdf, bmy, 3/30/04)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) FLG_MBL     (LOGICAL) : Mobilization candidate flag [unitless]
!  (2 ) MSS_FRC_CLY (REAL*8 ) : Mass fraction clay          [fraction]
!  (3 ) MSS_FRC_SND (REAL*8 ) : Mass fraction sand          [fraction]
!
!  Arguments as Output:
!  ============================================================================
!  (4 ) VWC_DRY     (REAL*8 ) : Dry volumetric water content (no E-T) [m3/m3]
!  (5 ) VWC_OPT     (REAL*8 ) : E-T optimal volumetric water content  [m3/m3]
!  (6 ) VWC_SAT     (REAL*8 ) : Saturated volumetric water content    [m3/m3]
!
!  NOTES:
!  (1 ) All I/O for this routine is time-invariant, thus, the hydrologic
!        properties could be computed once at initialization.  However,
!        FLG_MBL is time-dependent, so we should keep this as-is.
!        (tdf, 10/27/03)
!******************************************************************************
!

      ! Arguments
      TYPE(HCO_State), POINTER :: HcoState
      LOGICAL, INTENT(IN)  :: FLG_MBL(HcoState%NX)
      REAL*8,  INTENT(IN)  :: MSS_FRC_CLY_SLC(HcoState%NX)
      REAL*8,  INTENT(IN)  :: MSS_FRC_SND_SLC(HcoState%NX)
      REAL*8,  INTENT(OUT) :: VWC_DRY(HcoState%NX)
      REAL*8,  INTENT(OUT) :: VWC_OPT(HcoState%NX)
      REAL*8,  INTENT(OUT) :: VWC_SAT(HcoState%NX)

      ! Local variables
      INTEGER              :: LON_IDX

      ! [frc] Exponent "b" for smp (clay-dependent)
      REAL*8               :: SMP_XPN_B(HcoState%NX)

      ! [mm H2O] Saturated soil matric potential (sand-dependent)
      REAL*8               :: SMP_SAT(HcoState%NX)

      !=================================================================
      ! HYD_PRP_GET begins here
      !=================================================================

      ! Initialize output values
      VWC_DRY(:) = 0.0D0
      VWC_OPT(:) = 0.0D0
      VWC_SAT(:) = 0.0D0

      ! Time-invariant soil hydraulic properties
      ! See Bon96 p. 98, implemented in CCM:lsm/lsmtci()
      DO LON_IDX = 1, HcoState%NX

         IF ( FLG_MBL(LON_IDX) ) THEN

           ! Exponent "b" for smp (clay-dependent) [fraction]
           SMP_XPN_B(LON_IDX) =
     &         2.91D0 +0.159D0 * MSS_FRC_CLY_SLC(LON_IDX) * 100.0D0

           ! NB: Adopt convention that matric potential is positive definite
           ! Saturated soil matric potential (sand-dependent) [mm H2O]
           SMP_SAT(LON_IDX) =
     &         10.0D0 * (10.0D0**(1.88D0-0.0131D0
     &                          * MSS_FRC_SND_SLC(LON_IDX)*100.0D0))

           ! Saturated volumetric water content (sand-dependent) ! [m3 m-3]
           VWC_SAT(LON_IDX)=
     &         0.489D0 - 0.00126D0 * MSS_FRC_SND_SLC(LON_IDX)*100.0D0

           ! [m3 m-3]
           VWC_DRY(LON_IDX) =

                ! Dry volumetric water content (no E-T)
     &          VWC_SAT(LON_IDX)*(316230.0D0/SMP_SAT(LON_IDX))
     &                       **(-1.0D0/SMP_XPN_B(LON_IDX))

           ! E-T optimal volumetric water content! [m3 m-3]
           VWC_OPT(LON_IDX) =
     &         VWC_SAT(LON_IDX)*(158490.0D0/SMP_SAT(LON_IDX))
     &                        **(-1.0D0/SMP_XPN_B(LON_IDX))
         ENDIF
      ENDDO

      ! Return to calling program
      END SUBROUTINE HYD_PRP_GET

!------------------------------------------------------------------------------

      SUBROUTINE CND_TRM_SOI_GET( HcoState,CND_TRM_SOI,FLG_MBL,LVL_DLT,
     &                            MSS_FRC_CLY_SLC, MSS_FRC_SND_SLC,
     &                            TPT_SOI,
     &                            VWC_DRY,     VWC_OPT,     VWC_SAT,
     &                            VWC_SFC )

!
!******************************************************************************
!  Subroutine CND_TRM_SOI_GET gets thermal properties of soil.  Currently this
!  routine is optimized for ground without snow-cover.  Although snow
!  thickness is read in, it is not currently used. (tdf, bmy, 3/30/04)
!
!  Arguments as Input:
!  ============================================================================
!  (3 ) lvl_dlt     (REAL*8 ) : Soil layer thickness                  [m    ]
!  (4 ) mss_frc_cly (REAL*8 ) : Mass fraction clay                    [frac.]
!  (5 ) mss_frc_snd (REAL*8 ) : Mass fraction sand                    [frac.]
!  (6 ) tpt_soi     (REAL*8 ) : Soil temperature                      [K    ]
!  (7 ) vwc_dry     (REAL*8 ) : Dry volumetric water content (no E-T) [m3/m3]
!  (8 ) vwc_opt     (REAL*8 ) : E-T optimal volumetric water content  [m3/m3]
!  (9 ) vwc_sat     (REAL*8 ) : Saturated volumetric water content    [m3/m3]
!  (10) vwc_sfc     (REAL*8 ) : Volumetric water content              [m3/m3]
!
!  Arguments as Output:
!  ============================================================================
!  (1 ) CND_TRM_SOI (REAL*8 ) : Soil thermal conductivity             [W/m/K]
!  (2 ) FLG_MBL     (LOGICAL) : Mobilization candidate flag           [flag ]
!
!  NOTES:
!******************************************************************************
!

      ! Arguments
      TYPE(HCO_State), POINTER :: HcoState
      LOGICAL, INTENT(IN)  :: FLG_MBL(HcoState%NX)
      REAL*8,  INTENT(IN)  :: MSS_FRC_CLY_SLC(HcoState%NX)
      REAL*8,  INTENT(IN)  :: MSS_FRC_SND_SLC(HcoState%NX)
      REAL*8,  INTENT(IN)  :: TPT_SOI(HcoState%NX)
      REAL*8,  INTENT(IN)  :: VWC_DRY(HcoState%NX)
      REAL*8,  INTENT(IN)  :: VWC_OPT(HcoState%NX)
      REAL*8,  INTENT(IN)  :: VWC_SAT(HcoState%NX)
      REAL*8,  INTENT(IN)  :: VWC_SFC(HcoState%NX)
      REAL*8,  INTENT(OUT) :: CND_TRM_SOI(HcoState%NX)
      REAL*8,  INTENT(OUT) :: LVL_DLT(HcoState%NX)

      !------------
      ! Parameters
      !------------

      ! Thermal conductivity of ice water [W m-1 K-1]
      REAL*8, PARAMETER    :: CND_TRM_H2O_ICE      = 2.2d0

      ! Thermal conductivity of liquid water [W m-1 K-1]
      REAL*8, PARAMETER    :: CND_TRM_H2O_LQD      = 0.6d0

      ! Thermal conductivity of snow Bon96 p. 77 [W m-1 K-1]
      REAL*8, PARAMETER    :: CND_TRM_SNW          = 0.34d0

      ! Soil layer thickness, top layer! [m]
      REAL*8, PARAMETER    :: LVL_DLT_SFC          = 0.1d0

      ! Temperature range of mixed phase soil [K]
      REAL*8, PARAMETER    :: TPT_DLT              = 0.5d0

      ! Latent heat of fusion of H2O at 0 C, standard [J kg-1]
      REAL*8, PARAMETER    :: LTN_HEAT_FSN_H2O_STD = 0.3336d06

      ! Liquid water density [kg/m3]
      REAL*8, PARAMETER    :: DNS_H2O_LQD_STD      = 1000.0d0

      ! Kelvin--Celsius scale offset Bol80 [K]
      REAL*8, PARAMETER    :: TPT_FRZ_PNT          = 273.15d0

      !-----------------
      ! Local variables
      !-----------------

      ! Longitude index
      INTEGER              :: LON_IDX

      ! Thermal conductivity of dry soil [W m-1 K-1]
      REAL*8               :: CND_TRM_SOI_DRY(HcoState%NX)

      ! Soil thermal conductivity, frozen [W m-1 K-1]
      REAL*8               :: CND_TRM_SOI_FRZ(HcoState%NX)

      ! Thermal conductivity of soil solids [W m-1 K-1]
      REAL*8               :: CND_TRM_SOI_SLD(HcoState%NX)

      ! Soil thermal conductivity, unfrozen [W m-1 K-1]
      REAL*8               :: CND_TRM_SOI_WRM(HcoState%NX)

      ! Volumetric latent heat of fusion [J m-3]
      REAL*8               :: LTN_HEAT_FSN_VLM(HcoState%NX)

      ! Bounded geometric bulk thickness of snow [m]
      REAL*8               :: SNW_HGT_BND

      !=================================================================
      ! CND_TRM_SOI_GET begins here!
      !=================================================================

      ! [m] Soil layer thickness
      LVL_DLT(:) = LVL_DLT_SFC

      ! [W m-1 K-1] Soil thermal conductivity
      CND_TRM_SOI(:) = 0.0D0

      ! Loop over longitude
      DO LON_IDX = 1, HcoState%NX
         IF ( FLG_MBL(LON_IDX) ) THEN

           ! Volumetric latent heat of fusion [J m-3]
           LTN_HEAT_FSN_VLM(LON_IDX) = VWC_SFC(LON_IDX)
     &         * LTN_HEAT_FSN_H2O_STD * DNS_H2O_LQD_STD

           !Thermal conductivity of soil solids Bon96 p. 77 [W/m/K]
           CND_TRM_SOI_SLD(LON_IDX) =
     &         ( 8.80D0 *MSS_FRC_SND_SLC(LON_IDX)
     &         + 2.92D0 *MSS_FRC_CLY_SLC(LON_IDX) )
     &         / (MSS_FRC_SND_SLC(LON_IDX)
     &         + MSS_FRC_CLY_SLC(LON_IDX))

           ! Thermal conductivity of dry soil Bon96 p. 77 [W/m/K]
           cnd_trm_soi_dry(lon_idx) = 0.15D0

           ! Soil thermal conductivity, unfrozen [W/m/K]
           CND_TRM_SOI_WRM(LON_IDX) =
     &          CND_TRM_SOI_DRY(LON_IDX)
     &         + ( CND_TRM_SOI_SLD(LON_IDX)
     &         ** (1.0D0-VWC_SAT(LON_IDX))
     &         * (CND_TRM_H2O_LQD ** VWC_SFC(LON_IDX) )
     &         - CND_TRM_SOI_DRY(LON_IDX) )
     &         * VWC_SFC(LON_IDX) / VWC_SAT(lon_idx)

           ! Soil thermal conductivity, frozen [W/m/K]
           CND_TRM_SOI_FRZ(LON_IDX) =
     &          CND_TRM_SOI_DRY(LON_IDX)
     &         + ( CND_TRM_SOI_SLD(LON_IDX)
     &         ** (1.0D0-VWC_SAT(LON_IDX))
     &         * (CND_TRM_H2O_ICE ** VWC_SFC(LON_IDX) )
     &         - CND_TRM_SOI_DRY(LON_IDX) )
     &         * VWC_SFC(LON_IDX) / VWC_SAT(LON_IDX)

           IF (TPT_SOI(LON_IDX) < TPT_FRZ_PNT-TPT_DLT) THEN
               ! Soil thermal conductivity [W/m/K]
               CND_TRM_SOI(LON_IDX) = CND_TRM_SOI_FRZ(LON_IDX)
           ENDIF

           IF ( (TPT_SOI(LON_IDX) >= TPT_FRZ_PNT-TPT_DLT)
     &          .AND. (TPT_SOI(LON_IDX) <= TPT_FRZ_PNT+TPT_DLT) )
     &     THEN

              ! Soil thermal conductivity [W/m/K]
              CND_TRM_SOI(LON_IDX) =
     &            CND_TRM_SOI_FRZ(LON_IDX)
     &            + (CND_TRM_SOI_FRZ(LON_IDX)
     &            - CND_TRM_SOI_WRM(LON_IDX) )
     &            * (TPT_SOI(LON_IDX)
     &              -TPT_FRZ_PNT+TPT_DLT)
     &            / (2.0D0*TPT_DLT)
           ENDIF

           IF (TPT_SOI(LON_IDX) > TPT_FRZ_PNT+TPT_DLT) THEN
              ! Soil thermal conductivity[W/m/K]
              CND_TRM_SOI(LON_IDX)=CND_TRM_SOI_WRM(LON_IDX)
           ENDIF

! Implement this later(??)
!cZ Blend snow into first soil layer
!cZ Snow is not allowed to cover dust mobilization regions
!cZ snw_hgt_bnd=min(snw_hgt(lon_idx),1.0D0) ! [m] Bounded geometric bulk thickness of snow
!cZ lvl_dlt_snw(lon_idx)=lvl_dlt(lon_idx)+snw_hgt_bnd ! O [m] Soil layer thickness
!cZ including snow Bon96 p. 77
!
!cZ cnd_trm_soi(lon_idx)= & ! [W m-1 K-1] Soil thermal conductivity Bon96 p. 77
!cZ cnd_trm_snw*cnd_trm_soi(lon_idx)*lvl_dlt_snw(lon_idx) &
!cZ       /(cnd_trm_snw*lvl_dlt(lon_idx)+cnd_trm_soi(lon_idx)*snw_hgt_bnd)

         ENDIF
      ENDDO

      END SUBROUTINE CND_TRM_SOI_GET

!------------------------------------------------------------------------------

      SUBROUTINE TRN_FSH_VPR_SOI_ATM_GET( HcoState, FLG_MBL,
     &                                    TPT_SOI,
     &                                    TPT_SOI_FRZ,
     &                                    TRN_FSH_VPR_SOI_ATM,
     &                                    VWC_DRY,
     &                                    VWC_OPT,
     &                                    VWC_SFC )
!
!******************************************************************************
!  Subroutine TRN_FSH_VPR_SOI_ATM_GET computes the factor describing effects
!  of soil texture and moisture on vapor transfer between soil and atmosphere.
!  Taken from Bon96 p. 59, CCM:lsm/surphys. (tdf, bmy, 3/30/04)
!
!  The TRN_FSH_VPR_SOI_ATM efficiency factor attempts to tie soil texture and
!  moisture properties to the vapor conductance of the soil-atmosphere system.
!  When the soil temperature is sub-freezing, the conductance describes the
!  resistance to vapor sublimation (or deposition) and transport through the
!  open soil pores to the atmosphere.
!
!  For warm soils, vapor transfer is most efficient at the optimal VWC for E-T
!  Thus when vwc_sfc = vwc_opt, soil vapor transfer is perfectly efficient
!  (trn_fsh_vpr_soi_atm = 1.0) so the soil does not contribute any resistance
!  to the surface vapor transfer.
!
!  When vwc_sfc > vwc_opt, the soil has an excess of moisture and, again,
!  vapor transfer is not limited by soil characteristics.
!  In fact, according to Bon96 p. 98, vwc_dry is only slightly smaller than
!  vwc_opt, so trn_fsh_vpr_soi_atm is usually either 0 or 1 and intermediate
!  efficiencies occur over only a relatively small range of VWC.
!
!  When vwc_sfc < vwc_dry, the soil matrix is subsaturated and acts as a
!  one-way sink for vapor through osmotic and capillary potentials.
!  In this case trn_fsh_vpr_soi_atm = 0, which would cause the surface
!  resistance rss_vpr_sfc to blow up, but this is guarded against and
!  rss_sfc_vpr is set to ~1.0e6*rss_aer_vpr instead.
!
!  Note that this formulation does not seem to allow vapor transfer from
!  the atmosphere to the soil when vwc_sfc < vwc_dry, even when
!  e_atm > esat(Tg).
!
!  Air at the apparent sink for moisture is has vapor pressure e_sfc
!  e_atm = Vapor pressure of ambient air at z = hgt_mdp
!  e_sfc = Vapor pressure at apparent sink for moisture at z = zpd + rgh_vpr
!  e_gnd = Vapor pressure at air/ground interface temperature
!  Air at the soil interface is assumed saturated, i.e., e_gnd = esat(Tg)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) FLG_MBL             (LOGICAL) : Mobilization candidate flag [unitless]
!  (2 ) TPT_SOI             (REAL*8 ) : Soil temperature            [K       ]
!  (3 ) TPT_SOI_FRZ         (REAL*8 ) : Temperature of frozen soil  [K       ]
!  (5 ) VWC_DRY             (REAL*8 ) : Dry volumetric WC (no E-T)  [m3/m3   ]
!  (6 ) VWC_OPT             (REAL*8 ) : E-T optimal volumetric WC   [m3/m3   ]
!  (7 ) VWC_SFC             (REAL*8 ) : Volumetric water content    [m3/m3   ]
!
!  Arguments as Output:
!  ============================================================================
!  (4 ) TRN_FSH_VPR_SOI_ATM (REAL*8 ) : Transfer efficiency of vapor from
!                                       soil to atmosphere [fraction]
!
!  NOTES:
!  (1 ) Updated comments, cosmetic changes.  Also force double-precision
!        with "D" exponents. (tdf, bmy, 3/30/04)
!******************************************************************************
!

      !----------------
      ! Arguments
      !----------------
      TYPE(HCO_State), POINTER :: HCoState
      LOGICAL, INTENT(IN)  :: FLG_MBL(HcoState%NX)
      REAL*8,  INTENT(IN)  :: TPT_SOI(HcoState%NX)
      REAL*8,  INTENT(IN)  :: TPT_SOI_FRZ
      REAL*8,  INTENT(IN)  :: VWC_DRY(HcoState%NX)
      REAL*8,  INTENT(IN)  :: VWC_OPT(HcoState%NX)
      REAL*8,  INTENT(IN)  :: VWC_SFC(HcoState%NX)
      REAL*8,  INTENT(OUT) :: TRN_FSH_VPR_SOI_ATM(HcoState%NX)

      !----------------
      ! Parameters
      !----------------

      ! Transfer efficiency of vapor from frozen soil to
      ! atmosphere CCM:lsm/surphy()  [fraction]
      REAL*8, PARAMETER    :: TRN_FSH_VPR_SOI_ATM_FRZ = 0.01D0

      !-----------------
      ! Local variables
      !-----------------
      INTEGER              :: LON_IDX

      !=================================================================
      ! TRN_FSH_VPR_SOI_ATM_GET
      !=================================================================
      TRN_FSH_VPR_SOI_ATM(:) = 0.0D0

      ! Loop over longitudes
      DO LON_IDX = 1, HcoState%NX

         ! If this is a mobilization candidate ...
         IF ( FLG_MBL(LON_IDX) ) THEN

           ! ... and if the soil is above freezing ...
           IF ( TPT_SOI(LON_IDX) > TPT_SOI_FRZ ) THEN

              ! Transfer efficiency of cvapor from soil to atmosphere [frac]
              ! CCM:lsm/surphys Bon96 p. 59
              TRN_FSH_VPR_SOI_ATM(LON_IDX) =
     &             MIN ( MAX(VWC_SFC(LON_IDX)-VWC_DRY(LON_IDX), 0.0D0)
     &             /(VWC_OPT(LON_IDX)-VWC_DRY(LON_IDX)), 1.0D0)

           ELSE

              ! [frc] Bon96 p. 59
              TRN_FSH_VPR_SOI_ATM(LON_IDX) = TRN_FSH_VPR_SOI_ATM_FRZ

           ENDIF
         ENDIF
      ENDDO

      ! Return to calling program
      END SUBROUTINE TRN_FSH_VPR_SOI_ATM_GET

!------------------------------------------------------------------------------

      SUBROUTINE BLM_MBL( HcoState, FLG_MBL, RGH_MMN,
     &                    WND_10M,  MNO_LNG, WND_FRC, RC )
!
!******************************************************************************
!  Subroutine BLM_MBL computes the boundary-layer exchange properties, given
!  the meteorology at the GEOS-CHEM layer midpoint.  This routine is optimized
!  for dust source regions: dry, bare, uncovered land.  Theory and algorithms:
!  Bonan (1996) CCM:lsm/surtem().  Stripped down version, based on adiabatic
!  approximation to U*.  (tdf, bmy, 3/30/04)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) FLG_MBL (LOGICAL) : Mobilization candidate flag  [unitless]
!  (2 ) RGH_MMN (REAL*8 ) : Roughness length momentum    [m       ]
!  (3 ) WND_10M (REAL*8 ) : 10 m wind speed              [m/s     ]
!
!  Arguments as Output:
!  ============================================================================
!  (4 ) MNO_LNG (REAL*8 ) : Monin-Obukhov length         [m       ]
!  (5 ) WND_FRC (REAL*8 ) : Surface friction velocity    [m/s     ]
!
!  NOTES:
!  (1 ) Updated comments, cosmetic changes.  Also force double-precision with
!        "D" exponents. (tdf, bmy, 3/30/04)
!******************************************************************************
!
      !-----------------
      ! Arguments
      !-----------------
      TYPE(HCO_State), POINTER :: HcoState
      LOGICAL, INTENT(IN)  :: FLG_MBL(HcoState%NX)
      REAL*8,  INTENT(IN)  :: RGH_MMN(HcoState%NX)
      REAL*8,  INTENT(IN)  :: WND_10M(HcoState%NX)
      REAL*8,  INTENT(OUT) :: MNO_LNG(HcoState%NX)
      REAL*8,  INTENT(OUT) :: WND_FRC(HcoState%NX)
      INTEGER, INTENT(INOUT) :: RC

      !-----------------
      ! Parameters
      !-----------------

      ! Prevents division by zero [unitless]
      REAL*8,  PARAMETER  :: EPS_DBZ     = 1.0d-6

      ! Minimum windspeed used for mobilization [m/s]
      REAL*8,  PARAMETER  :: WND_MIN_MBL = 1.0d0

      ! Roughness length momentum for erodible surfaces [m]
      ! MaB95 p. 16420, GMB98 p. 6205
      REAL*8,  PARAMETER  :: RGH_MMN_MBL  = 100.0d-6

      ! Reference height for mobilization processes [m]
      REAL*8, PARAMETER   :: HGT_RFR       = 10.0d0

      !-----------------
      ! Local variables
      !-----------------

      ! Counting index for lon
      INTEGER             :: LON_IDX

      ! Denominator of Monin-Obukhov length Bon96 p. 49
      REAL*8              :: MNO_DNM

      ! Surface layer mean wind speed [m/s]
      REAL*8              :: WND_MDP_BND(HcoState%NX)

      ! denominator for wind friction velocity
      REAL*8              :: WND_FRC_DENOM

      ! For error handling
      CHARACTER(LEN=255)  :: MSG

      !=================================================================
      ! BLM_MBL begins here!
      !=================================================================

      ! Initialize
      MNO_LNG(:) = 0.0D0
      WND_FRC(:) = 0.0D0

      ! Loop over longitudes
      DO LON_IDX = 1, HcoState%NX

         ! Surface layer mean wind speed bounded [m/s]
         WND_MDP_BND(LON_IDX) =
     &        MAX(WND_10M(LON_IDX), WND_MIN_MBL)

         ! Friction velocity (adiabatic approximation  S&P equ. 16.57,
         ! tdf 10/27/2K3 -- Sanity check
         IF ( RGH_MMN(LON_IDX) <= 0.0 ) THEN
            MSG = 'RGH_MMN <= 0.0'
            CALL HCO_ERROR(HcoState%Config%Err,MSG,RC,THISLOC='BLM_MBL')
            RETURN
         ENDIF

         ! Distinguish between mobilisation candidates and noncandidates
         IF ( FLG_MBL(LON_IDX) ) THEN
            WND_FRC_DENOM = HGT_RFR / RGH_MMN_MBL      ! z = 10 m
         ELSE
            WND_FRC_DENOM = HGT_RFR / RGH_MMN(LON_IDX) ! z = 10 m
         ENDIF

         ! Sanity check
         IF ( WND_FRC_DENOM <= 0.0 ) THEN
            MSG = 'WND_FRC_DENOM <= 0.0'
            CALL HCO_ERROR(HcoState%Config%Err,MSG,RC,THISLOC='BLM_MBL')
            RETURN
         ENDIF

         ! Take natural LOG of WND_FRC_DENOM
         WND_FRC_DENOM    = LOG(WND_FRC_DENOM)

         ! Convert to [m/s]
         WND_FRC(LON_IDX) = WND_MDP_BND(LON_IDX) * CST_VON_KRM
     &                    / WND_FRC_DENOM

         ! Denominator of Monin-Obukhov length Bon96 p. 49
         ! Set denominator of Monin-Obukhov length to minimum value
         MNO_DNM = EPS_DBZ

         ! Monin-Obukhov length Bon96 p. 49 [m]
         MNO_LNG(LON_IDX) = -1.0D0 * (WND_FRC(LON_IDX)**3.0D0)
     &                       /MNO_DNM

         ! Override for non mobilisation candidates
         IF ( .NOT. FLG_MBL(LON_IDX) ) THEN
            WND_FRC(LON_IDX) = 0.0D0
         ENDIF
      ENDDO

      ! Return w/ success
      RC = HCO_SUCCESS

      END SUBROUTINE BLM_MBL

!------------------------------------------------------------------------------

      LOGICAL FUNCTION ORO_IS_OCN( ORO_VAL )
!
!******************************************************************************
!  Function ORO_IS_OCN returns TRUE if a grid box contains more than 50%
!  ocean. (tdf, bmy, 3/30/04)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) ORO_VAL (REAL*8) : Orography at a grid box (0=ocean; 1=land; 2=ice)
!
!  NOTES:
!******************************************************************************
!
      ! Arguments
      REAL*8, INTENT(IN) :: ORO_VAL

      !=================================================================
      ! ORO_IS_OCN begins here!
      !=================================================================
      ORO_IS_OCN = ( NINT( ORO_VAL ) == 0 )

      ! Return to calling program
      END FUNCTION ORO_IS_OCN

!------------------------------------------------------------------------------

      LOGICAL FUNCTION ORO_IS_LND( ORO_VAL )
!
!******************************************************************************
!  Function ORO_IS_LND returns TRUE if a grid box contains more than 50%
!  land. (tdf, bmy, 3/30/04, 3/1/05)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) ORO_VAL (REAL*8) : Orography at a grid box (0=ocean; 1=land; 2=ice)
!
!  NOTES:
!  (1 ) Bug fix: Replaced ": :" with "::" in order to prevent compile error
!        on Linux w/ PGI compiler.  (bmy, 3/1/05)
!******************************************************************************
!
      ! Arguments
      REAL*8, INTENT(IN) :: ORO_VAL

      !=================================================================
      ! ORO_IS_OCN begins here!
      !=================================================================
      ORO_IS_LND = ( NINT( ORO_VAL ) == 1 )

      ! Return to calling program
      END FUNCTION ORO_IS_LND

!------------------------------------------------------------------------------

      LOGICAL FUNCTION ORO_IS_ICE( ORO_VAL )
!
!******************************************************************************
!  Function ORO_IS_LND returns TRUE if a grid box contains more than 50%
!  ice. (tdf, bmy, 3/30/04)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) ORO_VAL (REAL*8) : Orography at a grid box (0=ocean; 1=land; 2=ice)
!
!  NOTES:
!******************************************************************************
!
      ! Arguments
      REAL*8, INTENT(IN) :: ORO_VAL

      !=================================================================
      ! ORO_IS_ICE begins here!
      !=================================================================
      ORO_IS_ICE = ( NINT( ORO_VAL ) == 2 )

      ! Return to calling program
      END FUNCTION ORO_IS_ICE

!------------------------------------------------------------------------------

      REAL*8 FUNCTION MNO_STB_CRC_HEAT_UNS_GET( SML_FNC_MMN_UNS_RCP )
!
!******************************************************************************
!  Function MNO_STB_CRC_HEAT_UNS_GET returns the stability correction factor
!  for heat (usually called PSI), given the reciprocal of the Monin-Obukhov
!  similarity function  (usually called PHI) for momentum in an unstable
!  atmosphere. (tdf, bmy, 3/30/04)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) sml_fnc_mmn_uns_rcp (REAL*8) : 1/(M-O similarity function) [fraction]
!
!  References:
!  ============================================================================
!  References are Ary88 p. 167, Bru82 p. 71, SeP97 p. 869,
!  Bon96 p. 52, BKL97 p. F1, LaP81 p. 325, LaP82 p. 466
!  Currently this function is BFB with CCM:dom/flxoce()
!
!  NOTES:
!  (1 ) Updated comments, cosmetic changes (bmy, 3/30/04)
!******************************************************************************
!
      ! Arguments
      REAL*8, INTENT(IN) :: SML_FNC_MMN_UNS_RCP

      !=================================================================
      ! MNO_STB_CRC_HEAT_UNS_GET
      !=================================================================
      MNO_STB_CRC_HEAT_UNS_GET = 2.0D0 *
     & LOG( ( 1.0D0+SML_FNC_MMN_UNS_RCP * SML_FNC_MMN_UNS_RCP) / 2.0D0 )

      ! Return to calling program
      END FUNCTION MNO_STB_CRC_HEAT_UNS_GET

!------------------------------------------------------------------------------

      REAL*8 FUNCTION MNO_STB_CRC_MMN_UNS_GET( SML_FNC_MMN_UNS_RCP )
!
!******************************************************************************
!  Function MNO_STB_CRC_MMN_UNS_GET returns the  stability correction factor
!  for momentum (usually called PSI), given the reciprocal of the
!  Monin-Obukhov similarity function (usually called PHI), for momentum in
!  an unstable atmosphere. (tdf, bmy, 3/30/04)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) SML_FNC_MMN_UNS_RCP (REAL*8) : 1/(M-O similarity function) [fraction]
!
!  References:
!  ============================================================================
!  References are Ary88 p. 167, Bru82 p. 71, SeP97 p. 869,
!  Bon96 p. 52, BKL97 p. F1, LaP81 p. 325, LaP82 p. 466
!  Currently this function is BFB with CCM:dom/flxoce()
!
!  NOTES:
!  (1 ) Updated comments, cosmetic changes (bmy, 3/30/04)
!******************************************************************************
!
      ! Arguments
      REAL*8, INTENT(IN) :: SML_FNC_MMN_UNS_RCP

      !=================================================================
      ! MNO_STB_CRC_MMN_UNS_GET begins here!
      !=================================================================
      MNO_STB_CRC_MMN_UNS_GET =
     &    LOG((1.0D0+SML_FNC_MMN_UNS_RCP*(2.0D0+SML_FNC_MMN_UNS_RCP))
     &       *(1.0D0+SML_FNC_MMN_UNS_RCP*SML_FNC_MMN_UNS_RCP)/8.0D0)
     &       -2.0D0*ATAN(SML_FNC_MMN_UNS_RCP)+1.571D0

      ! Return to calling program
      END FUNCTION MNO_STB_CRC_MMN_UNS_GET

!------------------------------------------------------------------------------

      REAL*8 FUNCTION XCH_CFF_MMN_OCN_NTR_GET( WND_10M_NTR )
!
!******************************************************************************
!  Function XCH_CFF_MMN_OCN_NTR_GET returns the Neutral 10m drag coefficient
!  over oceans. (tdf, bmy, 3/30/04)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) WIND_10M_NTR (REAL*8) : Wind speed @ 10 m[m/s]
!
!  References:
!  ============================================================================
!  LaP82 CCM:dom/flxoce(), NOS97 p. I2
!
!  NOTES:
!  (1 ) Updated comments, cosmetic changes (bmy, 3/30/04)
!******************************************************************************
!
      ! Arguments
      REAL*8, INTENT(IN) :: WND_10M_NTR

      !=================================================================
      ! XCH_CFF_MMN_OCN_NTR_GET begins here!
      !=================================================================
      XCH_CFF_MMN_OCN_NTR_GET = 0.0027D0    / WND_10M_NTR + 0.000142D0
     &                        + 0.0000764D0 * WND_10M_NTR

      ! REturn to calling program
      END FUNCTION XCH_CFF_MMN_OCN_NTR_GET

!------------------------------------------------------------------------------

      SUBROUTINE RGH_MMN_GET( HcoState,Inst,ORO, RGH_MMN,
     &                        SFC_TYP_SLC, SNW_FRC, WND_10M, RC )
!
!******************************************************************************
!  Subroutine RGH_MMN_GET sets the roughness length. (tdf, bmy, 3/30/04)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) ORO     (INTEGER) : Orography (0=ocean; 1=land; 2=ice)    [unitless]
!  (3 ) SFC_TYP (REAL*8 ) : LSM surface type (0..28)              [unitless]
!  (4 ) SNW_FRC (REAL*8 ) : Fraction of surface covered by snow   [fraction]
!  (5 ) WND_10M (REAL*8 ) : 10 m wind speed                       [m/s     ]
!
!  Arguments as Output:
!  ============================================================================
!  (2 ) RGH_MMN (REAL*8 ) : Roughness length momentu              [m       ]
!
!  NOTES:
!  (1 ) Updated comments, cosmetic changes.  Also now force double-precision
!        with "D" exponents (bmy, 3/30/04)
!******************************************************************************
!

      !-----------------
      ! Arguments
      !-----------------
      TYPE(HCO_State), POINTER  :: HcoState
      TYPE(MyInst),    POINTER  :: Inst
      INTEGER, INTENT(IN)  :: SFC_TYP_SLC(HcoState%NX)
      REAL*8,  INTENT(IN)  :: ORO(HcoState%NX)
      REAL*8,  INTENT(IN)  :: SNW_FRC(HcoState%NX)
      REAL*8,  INTENT(IN)  :: WND_10M(HcoState%NX)
      REAL*8,  INTENT(OUT) :: RGH_MMN(HcoState%NX)
      INTEGER, INTENT(INOUT) :: RC

      !-----------------
      ! Parameters
      !-----------------

      ! Roughness length over frozen lakes Bon96 p. 59 [m]
      REAL*8,  PARAMETER   :: RGH_MMN_ICE_LAK = 0.04d0

      ! Roughness length over ice, bare ground, wetlands Bon96 p. 59 [m]
      REAL*8,  PARAMETER   :: RGH_MMN_ICE_LND = 0.05d0

      ! Roughness length over sea ice BKL97 p. F-3 [m]
      REAL*8,  PARAMETER   :: RGH_MMN_ICE_OCN = 0.0005d0

      ! Roughness length over unfrozen lakes Bon96 p. 59 [m]
      REAL*8,  PARAMETER   :: RGH_MMN_LAK_WRM = 0.001d0

      ! Roughness length over snow Bon96 p. 59 CCM:lsm/snoconi.F ! [m]
      REAL*8,  PARAMETER   :: RGH_MMN_SNW     = 0.04d0

      ! Minimum windspeed for momentum exchange
      REAL*8,  PARAMETER   :: WND_MIN_DPS     = 1.0d0

      !-----------------
      ! Local variables
      !-----------------

      ! [idx] Longitude index array (sea ice)
      INTEGER              :: ICE_IDX(HcoState%NX)

      ! [nbr] Number of sea ice points
      INTEGER              :: ICE_NBR

      ! [Idx] Counting index
      INTEGER              :: IDX_IDX

      ! [idx] Longitude index array (land)
      INTEGER              :: LND_IDX(HcoState%NX)

      ! [nbr] Number of land points
      INTEGER              :: LND_NBR

      ! [idx] Counting index
      INTEGER              :: LON_IDX

      ! [idx] Longitude index array (ocean)
      INTEGER              :: OCN_IDX(HcoState%NX)

      ! [nbr] Number of ocean points
      INTEGER              :: OCN_NBR

      ! [idx] Plant type index
      INTEGER              :: PLN_TYP_IDX

      ! [idx] Surface type index
      INTEGER              :: SFC_TYP_IDX

      ! [idx] Surface sub-gridscale index
      INTEGER              :: SGS_IDX

      ! [m] Roughness length of current sub-gridscale
      REAL*8               :: RLM_CRR

      ! [m s-1] Bounded wind speed at 10 m
      REAL*8               :: WND_10M_BND

      ! [frc] Neutral 10 m drag coefficient over ocean
      REAL*8               :: XCH_CFF_MMN_OCN_NTR

      ! Momentum roughness length [m]
      REAL*8 :: Z0MVT(MVT) = (/ 0.94d0, 0.77d0, 2.62d0, 1.10d0, 0.99d0,
     &                          0.06d0, 0.06d0, 0.06d0, 0.06d0, 0.06d0,
     &                          0.06d0, 0.06d0, 0.06d0, 0.00d0 /)

      ! Displacement height (fn of plant type)
      REAL*8 :: ZPDVT(MVT)  = (/ 11.39d0, 9.38d0, 23.45d0, 13.40d0,
     &                           12.06d0, 0.34d0,  0.34d0,  0.34d0,
     &                            0.34d0, 0.34d0,  0.34d0,  0.34d0,
     &                            0.34d0, 0.00d0 /)

      !=================================================================
      ! RGH_MMN_GET begins here
      !=================================================================
      RGH_MMN(:) = 0.0D0

      ! Count ocean grid boxes
      OCN_NBR = 0
      DO LON_IDX = 1, HcoState%NX
         IF ( ORO_IS_OCN( ORO(LON_IDX) ) ) THEN
            OCN_NBR          = OCN_NBR + 1
            OCN_IDX(OCN_NBR) = LON_IDX
         ENDIF
      ENDDO

      ! Count ice grid boxes
      ICE_NBR = 0
      DO LON_IDX = 1, HcoState%NX
         IF ( ORO_IS_ICE( ORO(LON_IDX) ) ) THEN
            ICE_NBR          = ICE_NBR+1
            ICE_IDX(ICE_NBR) = LON_IDX
         ENDIF
      ENDDO

      ! Count land grid boxes
      LND_NBR = 0
      DO LON_IDX = 1, HcoState%NX
         IF ( ORO_IS_LND( ORO(LON_IDX) ) ) THEN
            LND_NBR          = LND_NBR + 1
            LND_IDX(LND_NBR) = LON_IDX
         ENDIF
      ENDDO

      !=================================================================
      ! Ocean points
      !=================================================================
      DO IDX_IDX = 1, OCN_NBR

         ! Longitude index of the ocean point
         LON_IDX = OCN_IDX(IDX_IDX)

         ! Convert wind speed to roughness length over ocean [m/s]
         WND_10M_BND = MAX( WND_MIN_DPS, WND_10M(LON_IDX) )

         !Approximation: neutral 10 m wind speed unavailable,
         ! use 10 m wind speed [fraction]
         XCH_CFF_MMN_OCN_NTR = XCH_CFF_MMN_OCN_NTR_GET(WND_10M_BND)

         ! BKL97 p. F-4, LaP81 p. 327 (14)  Ocean Points [m]
         RGH_MMN(LON_IDX)=10.0D0
     &       * EXP(-CST_VON_KRM / SQRT(XCH_CFF_MMN_OCN_NTR))
      ENDDO

      !=================================================================
      ! Sea ice points
      !=================================================================
      DO IDX_IDX = 1, ICE_NBR
         LON_IDX = ICE_IDX(IDX_IDX)
         RGH_MMN(LON_IDX) = SNW_FRC(LON_IDX) * RGH_MMN_SNW
     &      +(1.0D0-SNW_FRC(LON_IDX)) * RGH_MMN_ICE_OCN ! [m] Bon96 p. 59
      ENDDO

      !=================================================================
      ! Land points
      !=================================================================
      DO IDX_IDX = 1, LND_NBR

         ! Longitude
         LON_IDX = LND_IDX(IDX_IDX)

         ! Store surface blend for current gridpoint, sfc_typ(lon_idx)
         SFC_TYP_IDX = SFC_TYP_SLC(LON_IDX)

         ! Inland lakes
         IF ( SFC_TYP_IDX == 0 ) THEN

            !fxm: Add temperature input and so ability to discriminate warm
            !     from frozen lakes here [m] Bon96 p. 59
            RGH_MMN(LON_IDX) = RGH_MMN_LAK_WRM

         ! Land ice
         ELSE IF ( SFC_TYP_IDX == 1 ) THEN

           ! [m] Bon96 p. 59
           RGH_MMN(LON_IDX) = SNW_FRC(LON_IDX)*RGH_MMN_SNW
     &        + (1.0D0-SNW_FRC(LON_IDX))*RGH_MMN_ICE_LND


         ! Normal land
         ELSE
           DO SGS_IDX = 1, 3

              ! Bare ground is pln_typ=14, ocean is pln_typ=0
              PLN_TYP_IDX = Inst%PLN_TYP(SFC_TYP_IDX,SGS_IDX)

              ! Bare ground
              IF ( PLN_TYP_IDX == 14 ) THEN

                 ! Bon96 p. 59 (glacial ice is same as bare ground)
                 RLM_CRR = SNW_FRC(LON_IDX) * RGH_MMN_SNW
     &           + (1.0D0-SNW_FRC(LON_IDX)) * RGH_MMN_ICE_LND ! [m]

              ! Regular plant type
              ELSE IF ( PLN_TYP_IDX > 0 ) THEN
                 RLM_CRR = SNW_FRC(LON_IDX) * RGH_MMN_SNW
     &           + (1.0D0-SNW_FRC(LON_IDX)) * Z0MVT(PLN_TYP_IDX)
                                                      ! [m] Bon96 p. 59

              ! Presumably ocean snuck through
              ELSE
                 CALL HCO_ERROR( HcoState%Config%Err,
     &                          'pln_typ_idx == 0', RC,
     &                           THISLOC='RGH_MMN_GET' )
                 RETURN
              ENDIF            ! endif

              ! Roughness length for normal land
              RGH_MMN(LON_IDX) = RGH_MMN(LON_IDX)       ! [m]
     &              + Inst%PLN_FRC(SFC_TYP_IDX,SGS_IDX) ! [frc]
     &              * RLM_CRR                           ! [m]

           ENDDO
         ENDIF
      ENDDO

      ! Return w/ success
      RC = HCO_SUCCESS

      ! Return to calling program
      END SUBROUTINE RGH_MMN_GET

!------------------------------------------------------------------------------

      SUBROUTINE SNW_FRC_GET( HcoState, SNW_HGT_LQD, SNW_FRC )
!
!******************************************************************************
!  Subroutine SNW_FRC_GET converts equivalent liquid water snow depth to
!  fractional snow cover.  Uses the snow thickness -> fraction algorithm of
!  Bon96.  (tdf bmy, 3/30/04)
!
!  Arguments as Input:
!  ===========================================================================
!  (1 ) snw_hgt_lqd (REAL*8) : Equivalent liquid water snow depth [m]
!
!  Arguments as Output:
!  ===========================================================================
!  (2 ) snw_frc     (REAL*8 ) : Fraction of surface covered by snow
!
!  NOTES:
!  (1 ) Updated comments, cosmetic changes.  Also now force double-precision
!        with "D" exponents. (bmy, 3/30/04)
!******************************************************************************
!

      !----------------
      ! Arguments
      !----------------
      TYPE(HCO_State), POINTER :: HcoState
      REAL*8, INTENT(IN)  :: SNW_HGT_LQD(HcoState%NX)
      REAL*8, INTENT(OUT) :: SNW_FRC(HcoState%NX)

      !----------------
      ! Parameters
      !----------------

      ! Note disparity in bulk snow density between CCM and LSM
      ! WiW80 p. 2724, 2725 has some discussion of bulk snow density
      !
      ! Bulk density of snow [kg m-3]
      REAL*8,  PARAMETER  :: DNS_H2O_SNW_GND_LSM = 250.0D0

      ! Standard bulk density of snow on ground [kg m-3]
      REAL*8,  PARAMETER  :: DNS_H2O_SNW_GND_STD = 100.0D0

      ! Geometric snow thickness for 100% coverage ! [m]
      REAL*8,  PARAMETER  :: SNW_HGT_THR         = 0.05D0

      ! Liquid water density! [kg/m3]
      REAL*8,  PARAMETER  :: DNS_H2O_LQD_STD     = 1000.0D0

      !-----------------
      ! Local variables
      !-----------------

      ! [idx] Counting index for lon
      INTEGER             :: LON_IDX

      ! [m] Geometric bulk thickness of snow
      REAL*8              :: SNW_HGT(HcoState%NX)

      ! Conversion factor from liquid water depth
      ! to geometric snow thickness [fraction]
      REAL*8              :: HGT_LQD_SNW_CNV

      !=================================================================
      ! SNW_FRC_GET begins here!
      !=================================================================

      ! Conversion factor from liquid water depth to
      ! geometric snow thickness [fraction]
      HGT_LQD_SNW_CNV = DNS_H2O_LQD_STD
     &                / DNS_H2O_SNW_GND_STD

      ! Fractional snow cover
      DO LON_IDX = 1, HcoState%NX

         ! Snow height [m]
         SNW_HGT(LON_IDX) = SNW_HGT_LQD(LON_IDX)
     &                    * HGT_LQD_SNW_CNV

         ! Snow fraction
         ! NB: CCM and LSM seem to disagree on this
         SNW_FRC(LON_IDX) = MIN(SNW_HGT(LON_IDX)/SNW_HGT_THR, 1.0D0)
      ENDDO

      ! Return to calling program
      END SUBROUTINE SNW_FRC_GET

!------------------------------------------------------------------------------

      SUBROUTINE WND_RFR_GET( HcoState, FLG_ORO, HGT_MDP, HGT_RFR,
     &                        HGT_ZPD,  MNO_LNG, WND_FRC, WND_MDP,
     &                        WND_MIN,  WND_RFR )
!
!******************************************************************************
!  Subroutine WND_RFR_GET interpolates wind speed at given height to wind
!  speed at reference height. (tdf, bmy, 3/30/04)
!
!  Arguments as Input:
!  ===========================================================================
!  (1 ) FLG_ORO (LOGICAL)  : Orography flag (mobilization flag)       [flag]
!  (2 ) HGT_MDP (REAL*8 )  : Midpoint height above surface            [m   ]
!  (3 ) HGT_RFR (REAL*8 )  : Reference height                         [m   ]
!  (4 ) HGT_ZPD (REAL*8 )  : Zero plane displacement                  [m   ]
!  (5 ) MNO_LNG (REAL*8 )  : Monin-Obukhov length                     [m   ]
!  (6 ) WND_FRC (REAL*8 )  : Surface friction velocity                [m/s ]
!  (7 ) WND_MDP (REAL*8 )  : Surface layer mean wind speed            [m/s ]
!  (8 ) WND_MIN (REAL*8 )  : Minimum windspeed                        [m/s ]
!
!  Arguments as Output:
!  ===========================================================================
!  (9 ) WND_RFR (REAL*8 )  : Wind speed at reference height           [m/s ]
!
!  NOTES:
!  (1 ) Updated comments, cosmetic changes.  Also now force double-precision
!        with "D" exponents. (bmy, 3/30/04)
!******************************************************************************
!

      !------------------
      ! Arguments
      !------------------
      TYPE(HCO_State), POINTER :: HcoState
      LOGICAL, INTENT(IN)  :: FLG_ORO(HcoState%NX)
      REAL*8,  INTENT(IN)  :: HGT_MDP(HcoState%NX)
      REAL*8,  INTENT(IN)  :: HGT_RFR
      REAL*8,  INTENT(IN)  :: HGT_ZPD(HcoState%NX)
      REAL*8,  INTENT(IN)  :: MNO_LNG(HcoState%NX)
      REAL*8,  INTENT(IN)  :: WND_FRC(HcoState%NX)
      REAL*8,  INTENT(IN)  :: WND_MDP(HcoState%NX)
      REAL*8,  INTENT(IN)  :: WND_MIN
      REAL*8,  INTENT(OUT) :: WND_RFR(HcoState%NX)

      !------------------
      ! Parameters
      !------------------

      ! Named index for lower (target) hght
      INTEGER, PARAMETER   :: RFR_HGT_IDX=1

      ! Named index for upper (known) hght
      INTEGER, PARAMETER   :: GCM_HGT_IDX=2

      !------------------
      ! Local variables
      !------------------

      ! [idx] Counting index
      INTEGER              :: IDX_IDX

      ! [idx] Counting index for lon
      INTEGER              :: LON_IDX

      ! Stability computation loop index
      INTEGER              :: LVL_IDX

      ! Valid indices
      INTEGER              :: VLD_IDX(HcoState%NX)

      ! [nbr] Number of valid indices
      INTEGER              :: VLD_NBR

      ! [frc] Monin-Obukhov stability correction momentum
      REAL*8               :: MNO_STB_CRC_MMN(HcoState%NX,2)

      ! [frc] Monin-Obukhov stability parameter
      REAL*8               :: MNO_STB_PRM(HcoState%NX,2)

      ! [frc] Reciprocal of similarity function
      !       for momentum, unstable atmosphere
      REAL*8               :: SML_FNC_MMN_UNS_RCP

      ! Term in stability correction computation
      REAL*8               :: TMP2

      ! Term in stability correction computation
      REAL*8               :: TMP3

      ! Term in stability correction computation
      REAL*8               :: TMP4

      ! [frc] Wind correction factor
      REAL*8               :: WND_CRC_FCT(HcoState%NX)

      ! [m-1] Reciprocal of reference height
      REAL*8               :: HGT_RFR_RCP

      !=================================================================
      ! WND_RFR_GET begins here!
      !=================================================================

      HGT_RFR_RCP = 1.0D0 / HGT_RFR ! [m-1]
      WND_RFR = WND_MIN             ! [m s-1]

      ! Compute horizontal wind speed at reference height
      DO LON_IDX = 1, HcoState%NX
         IF (FLG_ORO(LON_IDX) .AND. HGT_ZPD(LON_IDX) < HGT_RFR) THEN

            ! Code uses notation of Bon96 p. 50, where lvl_idx=1
            ! is 10 m ref. hgt, lvl_idx=2 is atm. hgt
            MNO_STB_PRM(LON_IDX,RFR_HGT_IDX) =
     &           MIN((HGT_RFR-HGT_ZPD(LON_IDX))
     &           /MNO_LNG(LON_IDX),1.0D0) ! [frc]

            MNO_STB_PRM(LON_IDX,GCM_HGT_IDX) =
     &           MIN((HGT_MDP(LON_IDX)-HGT_ZPD(LON_IDX))
     &           /MNO_LNG(LON_IDX),1.0D0) ! [frc]

            DO LVL_IDX = 1, 2
               IF (MNO_STB_PRM(LON_IDX,LVL_IDX) < 0.0D0) THEN
                  SML_FNC_MMN_UNS_RCP = (1.0D0 - 16.0D0
     &                 * MNO_STB_PRM(LON_IDX,LVL_IDX))**0.25D0
                  TMP2 = LOG((1.0D0 + SML_FNC_MMN_UNS_RCP
     &                 * SML_FNC_MMN_UNS_RCP)/2.0D0)
                  TMP3 = LOG((1.0D0 + SML_FNC_MMN_UNS_RCP)/2.0D0)
                  MNO_STB_CRC_MMN(LON_IDX,LVL_IDX) =
     &                 2.0D0 * TMP3 + TMP2 - 2.0D0
     &                 * ATAN(SML_FNC_MMN_UNS_RCP) + 1.5707963
               ELSE             ! not stable
                  MNO_STB_CRC_MMN(LON_IDX,LVL_IDX) = -5.0D0
     &                 * MNO_STB_PRM(LON_IDX,LVL_IDX)
               ENDIF            ! stable
            ENDDO              ! end loop over lvl_idx

           TMP4 = LOG( (HGT_MDP(LON_IDX)-HGT_ZPD(LON_IDX))
     &          / (HGT_RFR-HGT_ZPD(LON_IDX)) )

           ! Correct neutral stability assumption
           WND_CRC_FCT(LON_IDX) = TMP4
     &             - MNO_STB_CRC_MMN(LON_IDX,GCM_HGT_IDX)
     &             + MNO_STB_CRC_MMN(LON_IDX,RFR_HGT_IDX) ! [frc]
           WND_RFR(LON_IDX) = WND_MDP(LON_IDX)-WND_FRC(LON_IDX)
     &             * CST_VON_KRM_RCP * WND_CRC_FCT(LON_IDX) ! [m s-1]
           WND_RFR(LON_IDX) = MAX(WND_RFR(LON_IDX),WND_MIN) ! [m s-1]
         ENDIF
      ENDDO

      ! Return to calling program
      END SUBROUTINE WND_RFR_GET

!------------------------------------------------------------------------------

      SUBROUTINE WND_FRC_THR_SLT_GET( HcoState, FLG_MBL,
     &                                DNS_MDP, WND_FRC_THR_SLT, RC )
!
!******************************************************************************
!  Subroutine WND_FRC_THR_SLT_GET ccmputes the dry threshold friction velocity
!  for saltation -- See Zender et al. expression (1) (tdf, bmy, 3/30/04)
!
!  Arguments as Input:
!  ===========================================================================
!  (1 ) FLG_MBL         (LOGICAL) : mobilisation flag
!  (2 ) DNS_MDP         (REAL*8 ) : Midlayer density [kg/m3]
!
!  Arguments as Output:
!  ===========================================================================
!  (3 ) WND_FRC_THR_SLT (REAL*8 ) : Threshold friction velocity
!                                    for saltation [m/s]
!
!  NOTES:
!  (1 ) Updated comments, cosmetic changes.  Also now force double-precision
!        with "D" exponents. (bmy, 3/30/04)
!******************************************************************************
!

      !----------------
      ! Arguments
      !----------------
      TYPE(HCO_State), POINTER :: HcoState
      LOGICAL, INTENT(IN)  :: FLG_MBL(HcoState%NX)
      REAL*8,  INTENT(IN)  :: DNS_MDP(HcoState%NX)
      REAL*8,  INTENT(OUT) :: WND_FRC_THR_SLT(HcoState%NX)
      INTEGER, INTENT(INOUT) :: RC

      !-----------------
      ! Parameters
      !-----------------

      ! [m] Optimal diameter for saltation,
      ! IvW82 p. 117 Fgr. 8, Pye87 p. 31, MBA97 p. 4388, SRL96 (2)
      REAL*8,  PARAMETER   :: DMT_SLT_OPT = 75.0d-6

      ! [kg m-3] Density of optimal saltation particles,
      ! MBA97 p. 4388 friction velocity for saltation
      REAL*8,  PARAMETER   :: DNS_SLT     = 2650.0d0

      !-----------------
      ! Local variables
      !-----------------

      ! [idx] Longitude Counting Index
      INTEGER              :: LON_IDX

      ! Threshold friction Reynolds number
      ! approximation for optimal size [frc]
      REAL*8               :: RYN_NBR

      !  Density ratio factor for saltation calculation
      REAL*8               :: DNS_FCT

      ! Interparticle cohesive forces factor for saltation calculation
      REAL*8               :: ALPHA, BETA, GAMMA, TMP1


      !=================================================================
      ! WND_FRC_THR_SLT_GET begins here!
      !=================================================================

      ! Initialize some variables
      ! MaB95 pzn. for Re*t(D_opt) circumvents iterative solution
      ! [frc] "B" MaB95 p. 16417 (5)

      ! [m/s] Threshold velocity
      WND_FRC_THR_SLT(:) = 0.0D0

      ! Threshold friction Reynolds number approximation for optimal size
      RYN_NBR = 0.38D0 + 1331.0D0
     &        * (100.0D0*DMT_SLT_OPT)**1.56D0

      ! tdf NB conversion of Dp to [cm]
      ! Given Re*t(D_opt), compute time independent factors contributing
      ! to u*t. IvW82 p. 115 (6) MaB95 p. 16417 (4) Interparticle cohesive
      ! forces. see Zender et al., Equ. (1).

      ! tdf introduced beta [fraction]
      BETA = 1.0D0+6.0D-07 / (DNS_SLT*GRV_SFC*(DMT_SLT_OPT**2.5D0))

      ! IvW82 p. 115 (6) MaB95 p. 16417 (4)
      DNS_FCT = DNS_SLT * GRV_SFC * DMT_SLT_OPT

      ! Error check
      IF ( RYN_NBR < 0.03D0 ) THEN
         CALL HCO_ERROR ( HcoState%Config%Err, 'RYN_NBR < 0.03', RC,
     &      THISLOC='WND_FRC_THR_SLT_GET' )
         RETURN

      ELSE IF ( RYN_NBR < 10.0D0 ) THEN

        ! IvW82 p. 114 (3), MaB95 p. 16417 (6)
        ! tdf introduced gamma [fraction]
        GAMMA = -1.0D0 + 1.928D0 * (RYN_NBR**0.0922D0)
        TMP1 = 0.129D0*0.129D0 * BETA / GAMMA

      ELSE

        ! ryn_nbr > 10.0D0
        ! IvW82 p. 114 (3), MaB95 p. 16417 (7)
        ! tdf introduced gamma [fraction]
        GAMMA = 1.0D0-0.0858D0 * EXP(-0.0617D0*(RYN_NBR-10.0D0))
        TMP1 = 0.12D0*0.12D0 * BETA * GAMMA * GAMMA

      ENDIF

      DO LON_IDX = 1, HcoState%NX

         ! Threshold friction velocity for saltation dry ground
         ! tdf introduced alpha
         ALPHA = DNS_FCT / DNS_MDP(LON_IDX)

         ! Added mobilisation constraint
         IF ( FLG_MBL(LON_IDX) ) THEN
            WND_FRC_THR_SLT(LON_IDX) =  SQRT(TMP1) * SQRT(ALPHA) ! [m s-1]
         ENDIF
      ENDDO

      ! Return w/ success
      RC = HCO_SUCCESS

      END SUBROUTINE WND_FRC_THR_SLT_GET

!------------------------------------------------------------------------------

      SUBROUTINE WND_RFR_THR_SLT_GET( HcoState, WND_FRC,
     &                                WND_FRC_THR_SLT, WND_MDP, WND_RFR,
     &                                WND_RFR_THR_SLT )
!
!******************************************************************************
!  Subroutine WND_RFR_THR_SLT_GET computes the threshold horizontal wind
!  speed at reference height for saltation. (tdf, bmy, 3/30/04)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) wnd_frc         (REAL*8) : Surface friction velocity              [m/s]
!  (2 ) wnd_frc_thr_slt (REAL*8) : Threshold friction vel. for saltation  [m/s]
!  (3 ) wnd_mdp         (REAL*8) : Surface layer mean wind speed          [m/s]
!  (4 ) wnd_rfr         (REAL*8) : Wind speed at reference height         [m/s]
!
!  Arguments as Output:
!  ============================================================================
!  (5 ) wnd_rfr_thr_slt (REAL*8) : Threshold 10m wind speed for saltation [m/s]
!
!  NOTES:
!  (1 ) Updated comments, cosmetic changes.
!******************************************************************************
!
      ! Arguments
      TYPE(HCO_State), POINTER :: HcoState
      REAL*8, INTENT(IN)  :: WND_FRC(HcoState%NX)
      REAL*8, INTENT(IN)  :: WND_FRC_THR_SLT(HcoState%NX)
      REAL*8, INTENT(IN)  :: WND_MDP(HcoState%NX)
      REAL*8, INTENT(IN)  :: WND_RFR(HcoState%NX)
      REAL*8, INTENT(OUT) :: WND_RFR_THR_SLT(HcoState%NX)

      ! Local variables
      INTEGER             :: I

      !=================================================================
      ! WND_RFR_THR_SLT_GET begins here
      !=================================================================
      DO I = 1, HcoState%NX

         ! A more complicated procedure would recompute mno_lng for
         ! wnd_frc_thr, and then integrate vertically from rgh_mmn+hgt_zpd
         ! to hgt_rfr.
         !
         ! wnd_crc_fct is (1/k)*[ln(z-D)/z0 - psi(zeta2) + psi(zeta1)]
         WND_RFR_THR_SLT(I) = WND_FRC_THR_SLT(I)
     &                      * WND_RFR(I) / WND_FRC(I)

      ENDDO

      ! Return to calling program
      END SUBROUTINE WND_RFR_THR_SLT_GET

!------------------------------------------------------------------------------

      SUBROUTINE VWC2GWC( HcoState, FLG_MBL, GWC_SFC, VWC_SAT, VWC_SFC )
!
!******************************************************************************
!  Subroutine VWC2GWC converts volumetric water content to gravimetric water
!  content -- assigned only for mobilisation candidates. (tdf, bmy, 3/30/04)
!
!  Arguments as Input:
!  ===========================================================================
!  (1 ) FLG_MBL (LOGICAL) : Mobilization candidate flag     [flag]
!  (3 ) VWC_SAT (REAL*8 ) : Saturated VWC (sand-dependent)  [m3/m3]
!  (4 ) VWC_SFC (REAL*8 ) : Volumetric water content!       [m3/m3
!
!  Arguments as Output:
!  ===========================================================================
!  (2 ) gwc_sfc (REAL*8 ) : Gravimetric water content       [kg/kg]
!
!  NOTES:
!  (1 ) Updated comments, cosmetic changes.  Also now forces double-precision
!        with "D" exponents. (tdf, bmy, 3/30/04)
!******************************************************************************
!

      !----------------
      ! Arguments
      !----------------
      TYPE(HCO_State), POINTER :: HcoState
      LOGICAL, INTENT(IN)  :: FLG_MBL(HcoState%NX)
      REAL*8,  INTENT(IN)  :: VWC_SAT(HcoState%NX)
      REAL*8,  INTENT(IN)  :: VWC_SFC(HcoState%NX)
      REAL*8,  INTENT(OUT) :: GWC_SFC(HcoState%NX)

      !----------------
      ! Parameters
      !----------------

      ! Dry density of soil ! particles (excluding pores) [kg/m3]
      REAL*8,  PARAMETER   :: DNS_PRT_SFC     = 2650.0d0

      ! liq. H2O density [kg/m3]
      REAL*8,  PARAMETER   :: DNS_H2O_LQD_STD = 1000.0d0

      !-----------------
      ! Local variables
      !-----------------

      ! Longitude index
      INTEGER              :: LON_IDX

      ! [kg m-3] Bulk density of dry surface soil
      REAL*8               :: DNS_BLK_DRY(HcoState%NX)

      !=================================================================
      ! VWC2GWC begins here!
      !=================================================================
      GWC_SFC(:)     = 0.0D0
      DNS_BLK_DRY(:) = 0.0D0

      ! Loop over longitudes
      DO LON_IDX = 1, HcoState%NX

         ! If this is a mobilization candidate then...
         IF ( FLG_MBL(LON_IDX) ) THEN

            ! Assume volume of air pores when dry equals saturated VWC
            ! This implies air pores are completely filled by water in
            ! saturated soil

            ! Bulk density of dry surface soil  [kg m-3]
            DNS_BLK_DRY(LON_IDX) = DNS_PRT_SFC
     &                           * ( 1.0d0 - VWC_SAT(LON_IDX) )

            ! Gravimetric water content [ kg kg-1]
            GWC_SFC(LON_IDX) = VWC_SFC(LON_IDX)
     &                       * DNS_H2O_LQD_STD
     &                       / DNS_BLK_DRY(LON_IDX)

         ENDIF
      ENDDO

      ! Return to calling program
      END SUBROUTINE VWC2GWC

!------------------------------------------------------------------------------

      SUBROUTINE FRC_THR_NCR_WTR_GET( HcoState,    FLG_MBL,
     &               FRC_THR_NCR_WTR, MSS_FRC_CLY_SLC, GWC_SFC )
!
!******************************************************************************
!  Subroutine FRC_THR_NCR_WTR_GET computes the factor by which soil moisture
!  increases threshold friction velocity. This parameterization is based on
!  FMB99. Zender et al., exp. (5). (tdf, bmy, 4/5/04)
!
!  Arguments as Input:
!  ===========================================================================
!  (1 ) FLG_MBL         (LOGICAL) : Mobilization candidate flag  [flags   ]
!  (3 ) MSS_FRC_CLY     (REAL*8 ) : Mass fraction of clay        [fraction]
!  (4 ) GWC_SFC         (REAL*8 ) : Gravimetric water content    [kg/kg   ]
!
!  Arguments as Output:
!  ===========================================================================
!  (2 ) FRC_THR_NCR_WTR (REAL*8 ) : Factor by which moisture increases
!                                    threshold friction velocity [fraction]
!
!  NOTES:
!  (1 ) Updated comments, cosmetic changes.  Also now forces double-precision
!        with "D" exponents. (tdf, bmy, 4/5/04)
!******************************************************************************
!

      ! Arguments
      TYPE(HCO_State), POINTER :: HcoState
      LOGICAL, INTENT(IN)  :: FLG_MBL(HcoState%NX)
      REAL*8,  INTENT(IN)  :: MSS_FRC_CLY_SLC(HcoState%NX)
      REAL*8,  INTENT(IN)  :: GWC_SFC(HcoState%NX)
      REAL*8,  INTENT(OUT) :: FRC_THR_NCR_WTR(HcoState%NX)

      ! local variables
      INTEGER              :: LON_IDX        ! [idx] Counting index
      REAL*8               :: GWC_THR(HcoState%NX) ! [kg/kg] Threshold GWC

      !=================================================================
      ! FRC_THR_NCR_WTR_GET begins here!
      !=================================================================

      ! Initialize
      frc_thr_ncr_wtr(:) = 1.0D0
      gwc_thr(:)         = 0.0D0

      ! Loop over longitudes
      DO LON_IDX = 1, HcoState%NX

         ! If this is a candidate for mobilization...
         IF ( FLG_MBL(LON_IDX) ) THEN

            !===========================================================
            ! Adjust threshold velocity for inhibition by moisture
            ! frc_thr_ncr_wtr(lon_idx)=exp(22.7D0*vwc_sfc(lon_idx))
            ! [frc] SRL96
            !
            ! Compute threshold soil moisture based on clay content
            ! GWC_THR=MSS_FRC_CLY*(0.17D0+0.14D0*MSS_FRC_CLY) [m3/m3]
            ! FMB99 p. 155 (14)
            !
            ! 19991105 remove factor of mss_frc_cly from gwc_thr to
            ! improve large scale behavior.
            !===========================================================

            ! [m3 m-3]
            GWC_THR(LON_IDX) = 0.17D0 + 0.14D0* MSS_FRC_CLY_SLC(LON_IDX)

            IF ( GWC_SFC(LON_IDX) > GWC_THR(LON_IDX) )
     &           FRC_THR_NCR_WTR(LON_IDX) = SQRT(1.0D0+1.21D0
     &           * (100.0D0 * (GWC_SFC(LON_IDX)-GWC_THR(LON_IDX)))
     &           ** 0.68D0)     ! [frc] FMB99 p. 155 (15)
         ENDIF
      ENDDO

      ! Return to calling program
      END SUBROUTINE FRC_THR_NCR_WTR_GET

!------------------------------------------------------------------------------

      SUBROUTINE FRC_THR_NCR_DRG_GET( HcoState, FRC_THR_NCR_DRG,
     &                                FLG_MBL,  Z0M, ZS0M, RC )
!
!******************************************************************************
!  Subroutine FRC_THR_NCR_DRG_GET computes factor by which surface roughness
!  increases threshold friction velocity. Zender et al., expression (3)
!  This parameterization is based on MaB95 and GMB98. (tdf, bmy, 4/5/04)
!
!  Arguments as Input:
!  ===========================================================================
!  (2 ) FLG_MBL         (LOGICAL) : Mobilization candidate flag
!  (3 ) Z0M             (REAL*8 ) : Roughness length momentum
!                                 :  for erodible surfaces [m]
!  (4 ) ZS0M            (REAL*8 ) : Smooth roughness length [m]
!
!  Arguments as Output:
!  ===========================================================================
!  (1 ) FRC_THR_NCR_DRG (REAL*8 ) : Factor by which surface roughness
!                                    increases threshold fric. velocity [frac]
!
!  NOTES:
!  (1 ) Updated comments, cosmetic changes.  Also now forces double-precision
!        with "D" exponents. (tdf, bmy, 4/5/04)
!******************************************************************************
!

      !-----------------
      ! Arguments
      !-----------------
      TYPE(HCO_State), POINTER :: HcoState
      LOGICAL, INTENT(IN)  :: FLG_MBL(HcoState%NX)
      REAL*8,  INTENT(IN)  :: Z0M
      REAL*8,  INTENT(IN)  :: ZS0M
      REAL*8,  INTENT(OUT) :: FRC_THR_NCR_DRG(HcoState%NX)
      INTEGER, INTENT(INOUT) :: RC

      !-----------------
      ! Local variables
      !-----------------

      ! [idx] Counting index
      integer lon_idx

      ! [frc] Efficient fraction of wind friction
      real*8 Feff

      ! [frc] Reciprocal of Feff
      real*8 Feff_rcp

      ! for error handling
      CHARACTER(LEN=255) :: MSG

      !=================================================================
      ! FRC_THR_NCR_DRG_GET begins here!
      !=================================================================

      FRC_THR_NCR_DRG(:) = 1.0D0

      ! Adjust threshold velocity for inhibition by roughness elements
      ! Zender et al. Equ. (3), fd.

      ! [frc] MaB95 p. 16420, GMB98 p. 6207
      FEFF = 1.0D0  - LOG( Z0M /ZS0M )
     &              / LOG( 0.35D0*( (0.1D0/ZS0M)**0.8D0) )

      ! Error check
      if ( FEFF <= 0.0D0 .OR. FEFF > 1.0D0 ) THEN
         MSG = 'Feff out of range!'
         CALL HCO_ERROR(HcoState%Config%Err,MSG, RC,
     &      THISLOC='FRC_THR_NC_DRG_GET' )
         RETURN
      ENDIF

      ! Reciprocal of FEFF [fraction]
      FEFF_RCP = 1.0D0 / FEFF

      ! Loop over longitudes
      DO LON_IDX = 1, HcoState%NX

         ! If this is a mobilization candidate...
         IF ( FLG_MBL(LON_IDX) ) THEN

            ! Save into FRC_THR_NCR_DRG
            FRC_THR_NCR_DRG(LON_IDX) = FEFF_RCP

            ! fxm: 19991012
            ! Set frc_thr_ncr_drg=1.0, equivalent to assuming mobilization
            ! takes place at smooth roughness length
            FRC_THR_NCR_DRG(LON_IDX) = 1.0D0

         ENDIF
      ENDDO

      ! Return w/ success
      RC = HCO_SUCCESS

      END SUBROUTINE FRC_THR_NCR_DRG_GET

!------------------------------------------------------------------------------

      SUBROUTINE WND_FRC_SLT_GET( HcoState, FLG_MBL, WND_FRC,
     &                            WND_FRC_SLT, WND_RFR, WND_RFR_THR_SLT)
!
!******************************************************************************
!  Subroutine WND_FRC_SLT_GET computes the saltating friction velocity.
!  Saltation increases friction speed by roughening surface, AKA "Owen's
!  effect".  This acts as a positive feedback to the friction speed.  GMB98
!  parameterized this feedback in terms of 10 m windspeeds, Zender et al.
!  equ. (4).  (tdf, bmy, 4/5/04, 1/25/07)
!
!  Arguments as Input:
!  ===========================================================================
!  (1 ) FLG_MBL         (LOGICAL) : Mobilization candidate flag
!  (2 ) WND_FRC         (REAL*8 ) : Surface friction velocity            [m/s]
!  (4 ) WND_RFR         (REAL*8 ) : Wind speed at reference height       [m/s]
!  (5 ) WND_RFR_THR_SLT (REAL*8 ) : Thresh. 10m wind speed for saltation [m/s]
!
!  Arguments as Output:
!  ===========================================================================
!  (3 ) WND_FRC_SLT     (REAL*8 ) : Saltating friction velocity          [m/s]
!
!  NOTES:
!  (1 ) Updated comments, cosmetic changes.  Also now forces double-precision
!        with "D" exponents. (tdf, bmy, 4/5/04)
!  (2 ) Now eliminate Owen effect (tdf, bmy, 1/25/07)
!******************************************************************************
!

      !-------------------
      ! Arguments
      !-------------------
      TYPE(HCO_State), POINTER :: HcoState
      LOGICAL, INTENT(IN)  :: FLG_MBL(HcoState%NX)
      REAL*8,  INTENT(IN)  :: WND_FRC(HcoState%NX)
      REAL*8,  INTENT(IN)  :: WND_RFR(HcoState%NX)
      REAL*8,  INTENT(IN)  :: WND_RFR_THR_SLT(HcoState%NX)
      REAL*8,  INTENT(OUT) :: WND_FRC_SLT(HcoState%NX)

      !-------------------
      ! Local variables
      !-------------------

      ! [idx] Counting index
      INTEGER              :: LON_IDX

      !---------------------------------------------------------------------
      ! Prior to 1/25/07:
      ! Eliminate Owen effect, so comment out this code (tdf, bmy, 1/25/07)
      !
      ! %%%%% DO NOT DELETE -- LEAVE THIS CODE COMMENTED OUT %%%%%
      !
      !! [m/s] Reference windspeed excess over threshold
      !REAL*8               :: WND_RFR_DLT
      !
      !! [m/s] Friction velocity increase from saltation
      !REAL*8               :: WND_FRC_SLT_DLT
      !---------------------------------------------------------------------

      !=================================================================
      ! WND_FRC_SLT_GET begins here!
      !=================================================================

      ! [m/s] Saltating friction velocity
      WND_FRC_SLT(:) = WND_FRC(:)

!------------------------------------------------------------------------------
! Prior to 1/25/07:
! Eliminate the Owen effect.  Note that the more computationally
! efficient way to do this is to just comment out the entire IF block.
! (tdf, bmy, 1/25/07)
!
! %%%%% DO NOT DELETE -- LEAVE THIS CODE COMMENTED OUT %%%%%
!
!      ! Loop over longitudes
!      DO LON_IDX = 1, HcoState%NX
!
!         ! If this is a mobilization candidate, then only
!         ! only apply Owen effect only when Uref > Ureft (tdf 4/5/04)
!         IF ( FLG_MBL(LON_IDX) .AND.
!     &        WND_RFR(LON_IDX) >= WND_RFR_THR_SLT(LON_IDX) ) THEN
!
!            !==================================================================
!            ! Saltation roughens the boundary layer, AKA "Owen's effect"
!            ! GMB98 p. 6206 Fig. 1 shows observed/computed u* dependence
!            ! on observed U(1 m).  GMB98 p. 6209 (12) has u* in cm s-1 and
!            ! U, Ut in m s-1, personal communication, D. Gillette, 19990529
!            ! With everything in MKS, the 0.3 coefficient in GMB98 (12)
!            ! becomes 0.003.  Increase in friction velocity due to saltation
!            ! varies as square of difference between reference wind speed
!            ! and reference threshold speed.
!            !==================================================================
!            WND_RFR_DLT = WND_RFR(LON_IDX) - WND_RFR_THR_SLT(LON_IDX)
!
!            ! Friction velocity increase from saltation GMB98 p. 6209 [m/s]
!            wnd_frc_slt_dlt = 0.003D0 * wnd_rfr_dlt * wnd_rfr_dlt
!
!            ! Saltation friction velocity, U*,s, Zender et al. Equ. (4).
!            WND_FRC_SLT(LON_IDX) = WND_FRC(LON_IDX)
!     &                           + WND_FRC_SLT_DLT ! [m s-1]
!
!            !
!ctdf Eliminate Owen effect                        tdf 01/13/2K5
!            wnd_frc_slt(:) = wnd_frc(:)
!
!         ENDIF
!      ENDDO
!------------------------------------------------------------------------------

      ! Return to calling program
      END SUBROUTINE WND_FRC_SLT_GET

!------------------------------------------------------------------------------

      SUBROUTINE FLX_MSS_CACO3_MSK( HcoState, ExtState,
     &                              DMT_VWR,
     &                              FLG_MBL,
     &                              FLX_MSS_VRT_DST_CACO3,
     &                              MSS_FRC_CACO3_SLC,
     &                              MSS_FRC_CLY_SLC,
     &                              MSS_FRC_SND_SLC, RC )
!
!******************************************************************************
!  Subroutine FLX_MSS_CACO3_MSK masks dust mass flux by CaCO3 mass fraction at
!  source.  Theory: Uses soil CaCO3 mass fraction from Global Soil Data Task,
!  1999 (Sch99).  Uses size dependent apportionment of CaCO3 from Claquin et
!  al, 1999 (CSB99). (tdf, bmy, 4/5/04)
!
!  Arguments as Input:
!  ===========================================================================
!  (1 ) DMT_VWR               (REAL*8 ) : Mass weighted diameter resolved [m]
!  (2 ) FLG_MBL               (LOGICAL) : Mobilization candidate flag
!  (3 ) FLX_MSS_VRT_DST_CACO3 (REAL*8 ) : Vert. mass flux of dust [kg/m2/s ]
!  (4 ) MSS_FRC_CACO3         (REAL*8 ) : Mass fraction of CaCO3  [fraction]
!  (5 ) MSS_FRC_CLY           (REAL*8 ) : Mass fraction of clay   [fraction]
!  (6 ) MSS_FRC_SND           (REAL*8 ) : Mass fraction of sand   [fraction]
!
!  Arguments as Output:
!  ===========================================================================
!  (3 ) FLX_MSS_VRT_DST_CACO3 (REAL*8 ) : Vertical mass flux of CaCO3 [kg/m2/s]
!
!  NOTES:
!  (1 ) Updated comments, cosmetic changes.  Also now forces double-precision
!        with "D" exponents. (tdf, bmy, 4/5/04)
!******************************************************************************
!

      !------------------
      ! Arguments
      !------------------
      TYPE(HCO_State), POINTER :: HcoState
      TYPE(Ext_State), POINTER :: ExtState
      LOGICAL, INTENT(IN)    :: FLG_MBL(HcoState%NX)
      REAL*8,  INTENT(IN)    :: DMT_VWR(NBINS)
      REAL*8,  INTENT(IN)    :: MSS_FRC_CACO3_SLC(HcoState%NX)
      REAL*8,  INTENT(IN)    :: MSS_FRC_CLY_SLC(HcoState%NX)
      REAL*8,  INTENT(IN)    :: MSS_FRC_SND_SLC(HcoState%NX)
      REAL*8,  INTENT(INOUT) :: FLX_MSS_VRT_DST_CACO3(HcoState%NX,NBINS)
      INTEGER, INTENT(INOUT) :: RC

      !------------------
      ! Parameters
      !------------------

      ! Maximum diameter of Clay soil texture CSB99 p. 22250 [m]
      REAL*8, PARAMETER      :: DMT_CLY_MAX = 2.0d-6

      ! Maximum diameter of Silt soil texture CSB99 p. 22250 [m]
      REAL*8, PARAMETER      :: DMT_SLT_MAX = 50.0d-6

      ! Density of CaCO3 http://www.ssc.on.ca/mandm/calcit.htm [kg/m3]
      REAL*8, PARAMETER      :: DNS_CACO3 = 2950.0d0

      !------------------
      ! Local variables
      !------------------

      ! [idx] Counting index
      INTEGER                :: M

      ! [idx] Counting index for lon
      INTEGER                :: LON_IDX

      ! [frc] Mass fraction of silt
      REAL*8                 :: MSS_FRC_SLT_SLC(HcoState%NX)

      ! [frc] Fraction of soil CaCO3 in size bin
      REAL*8                 :: MSS_FRC_CACO3_SZ_CRR

      ! [frc] Fraction of CaCO3 in clay
      REAL*8                 :: MSS_FRC_CACO3_CLY

      ! [frc] Fraction of CaCO3 in silt
      REAL*8                 :: MSS_FRC_CACO3_SLT

      ! [frc] Fraction of CaCO3 in sand
      REAL*8                 :: MSS_FRC_CACO3_SND

      ! Error handling
      CHARACTER(LEN=255)     :: MSG

      !=================================================================
      ! FLX_MSS_CACO3_MSK
      !=================================================================

      ! INITIALIZE
      MSS_FRC_SLT_SLC(:) = 0.0D0

      ! Loop over dust bins
      DO M = 1, NBINS

         ! Loop over longitudes
         DO LON_IDX = 1, HcoState%NX

            !===========================================================
            ! Simple technique is to mask dust mass by tracer mass
            ! fraction.  The model transports (hence conserves) CaCO3
            ! rather than total dust itself.  The method assumes source,
            ! transport, and removal processes are linear with tracer
            ! mass
            !===========================================================

            ! If this is a mobilization candidate, then...
            IF ( FLG_MBL(LON_IDX) ) THEN

               ! 20000320: Currently this is only process in
               ! dust model requiring mss_frc_slt

               ! [frc] Mass fraction of silt
               MSS_FRC_SLT_SLC(LON_IDX) =
     &              MAX(0.0D0, 1.0D0 -MSS_FRC_CLY_SLC(LON_IDX)
     &                               -MSS_FRC_SND_SLC(LON_IDX))

               ! CSB99 showed that CaCO3 is not uniformly distributed
               ! across sizes.  There is more CaCO3 per unit mass of
               ! silt than per unit mass of clay.

               ! Fraction of CaCO3 in clay CSB99 p. 22249 Figure 1b
               MSS_FRC_CACO3_CLY = MAX(0.0D0,-0.045D0+0.5D0
     &                           * MIN(0.5D0,MSS_FRC_CLY_SLC(LON_IDX)))

               ! Fraction of CaCO3 in silt CSB99 p. 22249 Figure 1a
               MSS_FRC_CACO3_SLT = MAX(0.0D0,-0.175D0+1.4D0
     &                           * MIN(0.5D0,MSS_FRC_SLT_SLC(LON_IDX)))

               ! Fraction of CaCO3 in sand CSB99 p. 22249 Figure 1a
               MSS_FRC_CACO3_SND = 1.0D0 - MSS_FRC_CACO3_CLY
     &                           - MSS_FRC_CACO3_SND

               ! Set CaCO3 fraction of total CaCO3 for each transport bin
               IF ( DMT_VWR(M) < DMT_CLY_MAX ) THEN

                  ! Transport bin carries Clay
                  ! Fraction of soil CaCO3 in size bin
                  MSS_FRC_CACO3_SZ_CRR = MSS_FRC_CACO3_CLY

               ELSE IF ( DMT_VWR(M) < DMT_SLT_MAX ) THEN

                  ! Transport bin carries Silt
                  ! Fraction of soil CaCO3 in size bin
                  MSS_FRC_CACO3_SZ_CRR = MSS_FRC_CACO3_SLT

               ELSE

                  ! Transport bin carries Sand
                  ! Fraction of soil CaCO3 in size bin
                  MSS_FRC_CACO3_SZ_CRR = MSS_FRC_CACO3_SND

               ENDIF

               ! Error checks
               IF ( MSS_FRC_CACO3_SZ_CRR < 0.0D0  .OR.
     &              MSS_FRC_CACO3_SZ_CRR > 1.0D0 ) THEN
                  MSG = 'mss_frc_CaC_s < 0.0.or.mss_frc_CaC_s > 1.0!'
                  CALL HCO_ERROR(HcoState%Config%Err,MSG, RC,
     &               THISLOC='FLX_MSS_CACO3_MSK' )
                  RETURN
               ENDIF

               IF ( MSS_FRC_CACO3_SLC(LON_IDX) < 0.0D0  .OR.
     &              MSS_FRC_CACO3_SLC(LON_IDX) > 1.0D0 ) THEN
                  MSG = 'mss_frc_CaCO3_s < 0.0.or.mss_frc_CaCO3 > 1.0!'
                  CALL HCO_ERROR(HcoState%Config%Err,MSG, RC,
     &               THISLOC='FLX_MSS_CACO3_MSK' )
                  RETURN
               ENDIF

               ! Convert dust flux to CaCO3 flux
               FLX_MSS_VRT_DST_CACO3(LON_IDX,M) =
     &              FLX_MSS_VRT_DST_CACO3(LON_IDX,M) ! [KG m-2 s-1]
     &              * MSS_FRC_CACO3_SLC(LON_IDX) ! [frc] Mass fraction of
                                            !       CaCO3 (at this location)
                    ! 20020925 fxm: Remove size dependence of CaCO3
     &              * 1.0D0

            ENDIF
         ENDDO
      ENDDO

      ! Return w/ success
      RC = HCO_SUCCESS

      END SUBROUTINE FLX_MSS_CACO3_MSK

!------------------------------------------------------------------------------

      SUBROUTINE FLX_MSS_HRZ_SLT_TTL_WHI79_GET( HcoState, DNS_MDP,
     &                                    FLG_MBL, QS_TTL,  U_S,  U_ST )
!
!******************************************************************************
!  Subroutine FLX_MSS_HRZ_SLT_TTL_WHI79_GET computes vertically integrated
!  streamwise mass flux of particles.  Theory: Uses method proposed by White
!  (1979). See Zender et al., expr (10).  fxm: use surface air density not
!  midlayer density (tdf, bmy, 4/5/04)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) DNS_MDP (REAL*8 ) : Midlayer density                           [g/m3  ]
!  (2 ) FLG_MBL (LOGICAL) : Mobilization candidate flag                [flag  ]
!  (4 ) U_S     (REAL*8 ) : Surface friction velocity                  [m/s   ]
!  (5 ) U_ST    (REAL*8 ) : Threshold friction spd for saltation       [m/s   ]
!
!  Arguments as Output:
!  ============================================================================
!  (3 ) QS_TTL  (REAL*8 ) : Vertically integrated streamwise mass flux [kg/m/s]
!
!  NOTES:
!  (1 ) Updated comments, cosmetic changes.  Also now forces double-precision
!        with "D" exponents. (tdf, bmy, 4/5/04)
!******************************************************************************
!

      !------------------
      ! Arguments
      !------------------
      TYPE(HCO_State), POINTER :: HcoState
      LOGICAL, INTENT(IN)  :: FLG_MBL(HcoState%NX)
      REAL*8,  INTENT(IN)  :: DNS_MDP(HcoState%NX)
      REAL*8,  INTENT(IN)  :: U_S(HcoState%NX)
      REAL*8,  INTENT(IN)  :: U_ST(HcoState%NX)
      REAL*8,  INTENT(OUT) :: QS_TTL(HcoState%NX)

      !------------------
      ! Parameters
      !------------------

      ! [frc] Saltation constant Whi79 p. 4648, MaB97 p. 16422
      REAL*8,  PARAMETER   :: CST_SLT = 2.61d0

      !------------------
      ! Local variables
      !------------------

      ! [frc] Ratio of wind friction threshold to wind friction
      real*8               :: U_S_rat

      ! [idx] Counting index for lon
      integer              :: lon_idx

      !=================================================================
      ! FLX_MSS_HRZ_SLT_TTL_WHI79_GET begins here!
      !=================================================================

      ! Initialize
      QS_TTL(:) = 0.0D0

      ! Loop over longitudes
      DO LON_IDX = 1, HcoState%NX

         ! If this is a mobilization candidate and the friction
         ! velocity is above the threshold for saltation...
         IF ( FLG_MBL(LON_IDX) .AND.
     &        U_S(LON_IDX) > U_ST(LON_IDX) ) THEN

            ! Ratio of wind friction threshold to wind friction
            U_S_RAT = U_ST(LON_IDX) / U_S(LON_IDX)

            ! Whi79 p. 4648 (19), MaB97 p. 16422 (28)
            QS_TTL(LON_IDX) =   ! [kg m-1 s-1]
     &           CST_SLT * DNS_MDP(LON_IDX) * (U_S(LON_IDX)**3.0D0)
     &           * (1.0D0-U_S_RAT) * (1.0D0+U_S_RAT)
     &            * (1.0D0+U_S_RAT) / GRV_SFC

         ENDIF
      ENDDO

      ! Return to calling program
      END SUBROUTINE FLX_MSS_HRZ_SLT_TTL_WHI79_GET

!------------------------------------------------------------------------------

      SUBROUTINE FLX_MSS_VRT_DST_TTL_MAB95_GET( HcoState,
     &                                          DST_SLT_FLX_RAT_TTL,
     &                                          FLG_MBL,
     &                                          FLX_MSS_HRZ_SLT_TTL,
     &                                          FLX_MSS_VRT_DST_TTL,
     &                                          MSS_FRC_CLY_SLC )
!
!******************************************************************************
!  Subroutine FLX_MSS_VRT_DST_TTL_MAB95_GET diagnoses total vertical mass flux
!  of dust from vertically integrated streamwise mass flux, Zender et al.,
!  expr. (11). (tdf, bmy, 4/5/04)
!
!  Theory: Uses clay-based method proposed by Marticorena & Bergametti (1995)
!  Their parameterization is based only on data for mss_frc_cly < 0.20
!  For clayier soils, dst_slt_flx_rat_ttl may behave dramatically differently
!  Whether this behavior changes when mss_frc_cly > 0.20 is unknown
!  Anecdotal evidence suggests vertical flux decreases for mss_frc_cly > 0.20
!  Thus we use min[mss_frc_cly,0.20] in MaB95 parameterization
!
!  Arguments as Input:
!  ============================================================================
!  (2 ) FLG_MBL             (LOGICAL) : Mobilization candidate flag
!  (3 ) FLX_MSS_HRZ_SLT_TTL (REAL*8 ) : Vertically integrated streamwise
!                                        mass flux [kg/m/s]
!  (5 ) MSS_FRC_CLY         (REAL*8 ) : Mass fraction clay [fraction]
!
!  Arguments as Output:
!  ============================================================================
!  (1 ) DST_SLT_FLX_RAT_TTL (REAL*8 ) : Ratio of vertical dust flux t
!                                       to streamwise mass flux [1/m]
!  (4 ) FX_MSS_VRT_DST_TTL  (REAL*8 ) : Total vert. mass flux of dust [kg/m2/s]
!
!  NOTES:
!  (1 ) Updated comments, cosmetic changes.  Also now forces double-precision
!        with "D" exponents. (tdf, bmy, 4/5/04)
!******************************************************************************
!

      !-----------------
      ! Arguments
      !-----------------
      TYPE(HCO_State), POINTER :: HcoState
      LOGICAL, INTENT(IN)  :: FLG_MBL(HcoState%NX)
      REAL*8,  INTENT(IN)  :: FLX_MSS_HRZ_SLT_TTL(HcoState%NX)
      REAL*8,  INTENT(IN)  :: MSS_FRC_CLY_SLC(HcoState%NX)
      REAL*8,  INTENT(OUT) :: DST_SLT_FLX_RAT_TTL(HcoState%NX)
      REAL*8,  INTENT(OUT) :: FLX_MSS_VRT_DST_TTL(HcoState%NX)

      !-----------------
      ! Local variables
      !-----------------

      ! [idx] Counting index for lon
      INTEGER              :: LON_IDX

      ! [frc] Mass fraction clay limited to 0.20
      REAL*8               :: MSS_FRC_CLY_VLD

      ! [frc] Natural log of 10
      REAL*8               :: LN10

      !=================================================================
      ! FLX_MSS_VRT_DST_TTL_MAB95_GET
      !=================================================================

      ! Initialize
      LN10                   = LOG(10.0D0)
      DST_SLT_FLX_RAT_TTL(:) = 0.0D0
      FLX_MSS_VRT_DST_TTL(:) = 0.0D0

      ! Loop over longitudes
      DO LON_IDX = 1, HcoState%NX

         ! If this is a mobilization candidate...
         IF ( FLG_MBL(LON_IDX) ) then

            ! 19990603: fxm: Dust production is EXTREMELY sensitive to
            ! this parameter, which changes flux by 3 orders of magnitude
            ! in 0.0 < mss_frc_cly < 0.20
            MSS_FRC_CLY_VLD = MIN(MSS_FRC_CLY_SLC(LON_IDX),0.2D0)  ! [frc]

            DST_SLT_FLX_RAT_TTL(LON_IDX) =           ! [m-1]
     &         100.0D0 * EXP(LN10*(13.4D0*MSS_FRC_CLY_VLD-6.0D0))
                                                     ! MaB95 p. 16423 (47)

            FLX_MSS_VRT_DST_TTL(LON_IDX) =           ! [kg M-1 s-1]
     &           FLX_MSS_HRZ_SLT_TTL(LON_IDX)
     &         * DST_SLT_FLX_RAT_TTL(LON_IDX)

         ENDIF
      ENDDO

      ! Return to calling program
      END SUBROUTINE FLX_MSS_VRT_DST_TTL_MAB95_GET

!------------------------------------------------------------------------------

      SUBROUTINE DST_PSD_MSS( OVR_SRC_SNK_FRC, MSS_FRC_SRC,
     &                        OVR_SRC_SNK_MSS, NBINS, DST_SRC_NBR )
!
!******************************************************************************
!  Subroutine DST_PSD_MSS computes OVR_SRC_SNK_MSS from OVR_SRC_SNK_FRC
!  and MSS_FRC_SRC. (tdf, bmy, 4/5/04)
!
!  Multiply ovr_src_snk_frc(src_idx,*) by mss_frc(src_idx) to obtain
!  absolute mass fraction mapping from source dists. to sink bins
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) OVR_SRC_SNK_FRC (REAL*8 ) : Mass overlap, Mij, Zender p. 5, Equ. 12
!  (2 ) MSS_FRC_SRC     (REAL*8 ) : Mass fraction in each mode (Table 1, M)
!  (4 ) NBINS         (INTEGER) : Number of GEOS_CHEM dust bins
!  (5 ) DST_SRC_NBR     (INTEGER) : Number of source modes
!
!  Arguments as Output:
!  ============================================================================
!  (3 ) OVR_SRC_SNK_MSS (REAL*8 ) : Mass of stuff ???
!
!  NOTES:
!  (1 ) Updated comments, cosmetic changes.  Also now forces double-precision
!        with "D" exponents. (tdf, bmy, 4/5/04)
!******************************************************************************
!
      !-----------------
      ! Arguments
      !-----------------
      INTEGER, INTENT(IN)  :: DST_SRC_NBR, NBINS
      REAL*8,  INTENT(IN)  :: OVR_SRC_SNK_FRC(DST_SRC_NBR,NBINS)
      REAL*8,  INTENT(IN)  :: MSS_FRC_SRC(DST_SRC_NBR)
      REAL*8,  INTENT(OUT) :: OVR_SRC_SNK_MSS(DST_SRC_NBR,NBINS)

      !-----------------
      ! Local variables
      !-----------------
      INTEGER              :: SRC_IDX, SNK_IDX
      REAL*8               :: MSS_FRC_TRN_DST_SRC(NBINS)
      REAL*8               :: OVR_SRC_SNK_MSS_TTL

      !=================================================================
      ! DST_PSD_MSS begins here!
      !=================================================================

      ! Fraction of vertical dust flux which is transported
      OVR_SRC_SNK_MSS_TTL = 0.0D0

      ! Fraction of transported dust mass at source
      DO SNK_IDX = 1, NBINS
         MSS_FRC_TRN_DST_SRC(SNK_IDX) = 0.0D0
      ENDDO

      DO SNK_IDX = 1, NBINS
      DO SRC_IDX = 1, DST_SRC_NBR
         OVR_SRC_SNK_MSS (SRC_IDX,SNK_IDX) = ! [frc]
     &        OVR_SRC_SNK_FRC (SRC_IDX,SNK_IDX)
     &        * MSS_FRC_SRC (SRC_IDX) ! [frc]
      ENDDO
      ENDDO

      ! Split double do loop into 2 parts      tdf 10/22/2K3
      DO SNK_IDX = 1, NBINS
      DO SRC_IDX = 1, DST_SRC_NBR

         ! [frc] Fraction of transported dust mass at source
         MSS_FRC_TRN_DST_SRC(SNK_IDX) =
     &        MSS_FRC_TRN_DST_SRC(SNK_IDX)
     &        + OVR_SRC_SNK_MSS(SRC_IDX,SNK_IDX)

         ! [frc] Compute total transported mass fraction of dust flux
         OVR_SRC_SNK_MSS_TTL = OVR_SRC_SNK_MSS_TTL
     &                       + OVR_SRC_SNK_MSS (SRC_IDX,snk_idx)
      ENDDO
      ENDDO

      ! Convert fraction of mobilized mass to fraction of transported mass
      DO SNK_IDX = 1, NBINS
         MSS_FRC_TRN_DST_SRC (SNK_IDX) =
     &        MSS_FRC_TRN_DST_SRC (SNK_IDX) / OVR_SRC_SNK_MSS_TTL
      ENDDO

      ! Return to calling program
      END SUBROUTINE DST_PSD_MSS

!------------------------------------------------------------------------------

      SUBROUTINE FLX_MSS_VRT_DST_PRT( Inst, NX, FLG_MBL,
     &                                FLX_MSS_VRT_DST,
     &                                FLX_MSS_VRT_DST_TTL )
!
!******************************************************************************
!  Subroutine FLX_MSS_VRT_DST_PRT partitions total vertical mass flux of dust
!  into transport bins.  Assumes a trimodal lognormal probability density
!  function (see Zender et al., p. 5). (tdf, bmy, 4/5/04)
!
!  DST_SRC_NBR  = 3 - trimodal size distribution in source c regions (p. 5)
!  OVR_SRC_SNK_MSS  [frc] computed in dst_psd_mss, called from dust_mod.f
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) FLG_MBL             (LOGICAL) : Mobilization candidate flag
!  (3 ) FLX_MSS_VRT_DST_TTL (REAL*8 ) : Total vert. mass flux of dust [kg/m2/s]
!
!  Arguments as Output:
!  ============================================================================
!  (2 ) FLX_MSS_VRT_DST     (REAL*8 ) : Vertical mass flux of dust [kg/m2/s]
!
!  NOTES:
!  (1 ) Updated comments, cosmetic changes.  Also now forces double-precision
!        with "D" exponents. (tdf, bmy, 4/5/04)
!******************************************************************************
!

      ! Arguments
      TYPE(MyInst), POINTER :: Inst
      INTEGER, INTENT(IN)   :: NX
      LOGICAL, INTENT(IN)   :: FLG_MBL(NX)
      REAL*8,  INTENT(IN)   :: FLX_MSS_VRT_DST_TTL(NX)
      REAL*8,  INTENT(OUT)  :: FLX_MSS_VRT_DST(NX,NBINS)

      ! Local variables
      INTEGER               :: LON_IDX   ! [idx] Counting index for lon
      INTEGER               :: SRC_IDX   ! [idx] Counting index for src
      INTEGER               :: SNK_IDX   ! [idx] Counting index for snk
      INTEGER               :: SNK_NBR   ! [nbr] Dimension size

      !=================================================================
      ! FLX_MSS_VRT_DST_PRT begins here!
      !=================================================================

      ! Initialize
      FLX_MSS_VRT_DST(:,:) = 0.0D0    ! [frc]

      ! Loop over longitudes (NB: Inefficient loop order)
      DO LON_IDX = 1, NX

         ! If this is a mobilization candidate...
         IF ( FLG_MBL(LON_IDX) ) THEN

            ! Loop over source & sink indices
            DO SNK_IDX = 1, NBINS
            DO SRC_IDX = 1, DST_SRC_NBR
               FLX_MSS_VRT_DST(LON_IDX,SNK_IDX) = ! [kg m-2 s-1]
     &              FLX_MSS_VRT_DST(LON_IDX,SNK_IDX)
     &              + Inst%OVR_SRC_SNK_MSS(SRC_IDX,SNK_IDX)
     &              * FLX_MSS_VRT_DST_TTL(LON_IDX)
            ENDDO
            ENDDO
         ENDIF
      ENDDO

      ! Return to calling program
      END SUBROUTINE FLX_MSS_VRT_DST_PRT

!------------------------------------------------------------------------------

      SUBROUTINE TM_2_IDX_WGT()

      ! routine eliminated: see original code
      END SUBROUTINE TM_2_IDX_WGT

!------------------------------------------------------------------------------

      SUBROUTINE LND_FRC_MBL_GET( HcoState,    DOY,
     &                            FLG_MBL,     LAT_RDN,
     &                            LND_FRC_DRY_SLC, LND_FRC_MBL, MBL_NBR,
     &                            ORO,         SFC_TYP_SLC,     SNW_FRC,
     &                            TPT_SOI,     TPT_SOI_FRZ, VAI_DST_SLC,
     &                            RC)
!
!******************************************************************************
!  Subroutine LND_FRC_MBL_GET returns the fraction of each GEOS-CHEM grid
!  box which is suitable for dust mobilization.  This routine is called
!  by DST_MBL. (tdf, bmy, 4/5/04, 1/13/10)
!
!  The DATE is used to obtain the time-varying vegetation cover.
!  Routine currently uses latitude slice of VAI from time-dependent surface
!  boundary dataset (tdf, 10/27/03).  LAI/VAI algorithm is from CCM:lsm/phenol
!  () Bon96.  The LSM data are mid-month values, i.e., valid on the 15th of !
!  the month.!
!
!  Criterion for mobilisation candidate (tdf, 4/5/04):
!  (1) first, must be a land point, not ocean, not ice
!  (2) second, it cannot be an inland lake, wetland or ice
!  (3) modulated by vegetation type
!  (4) modulated by subgridscale wetness
!  (5) cannot be snow covered
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) DOY         (REAL*8 ) : Day of year                         [1.0-366.0]
!  (3 ) LAT_RDN     (REAL*8 ) : Latitude                            [radians  ]
!  (4 ) LND_FRC_DRY (REAL*8 ) : Dry land fraction                   [fraction ]
!  (7 ) ORO         (REAL*8 ) : Orography: land/ocean/ice           [flags    ]
!  (8 ) SFC_TYP     (INTEGER) : LSM surface type (0..28)            [unitless ]
!  (9 ) SNW_FRC     (REAL*8 ) : Fraction of surface covered by snow [fraction ]
!  (10) TPT_SOI     (REAL*8 ) : Soil temperature                    [K        ]
!  (11) TPT_SOI_FRZ (REAL*8 ) : Temperature of frozen soil          [K        ]
!  (12) VAI_DST     (REAL*8 ) : Vegetation area index, one-sided    [m2/m2    ]
!
!  Arguments as Output:
!  ============================================================================
!  (2 ) FLG_MBL     (LOGICAL) : Mobilization candidate flag         [flag     ]
!  (5 ) LND_FRC_MBL (REAL*8 ) : Bare ground fraction                [fraction ]
!  (6 ) MBL_NBR     (INTEGER) : Number of mobilization candidates   [unitless ]
!
!  NOTES:
!  (1 ) Updated comments, cosmetic changes.  Also now forces double-precision
!        with "D" exponents. (tdf, bmy, 4/5/04)
!  (2 ) For the GOCART source function, we don't use VAI, so set FLG_VAI_TVBDS
!         = .FALSE. and disable calls to ERROR_STOP (tdf, bmy, 1/25/07)
!  (3 ) Modification for GEOS-4 1 x 1.25 grids (lok, bmy, 1/13/10)
!******************************************************************************
!

      !------------------
      ! Arguments
      !------------------
      TYPE(HCO_State), POINTER  :: HcoState
      INTEGER, INTENT(IN)  :: SFC_TYP_SLC(HcoState%NX)
      REAL*8,  INTENT(IN)  :: DOY
      REAL*8,  INTENT(IN)  :: LAT_RDN
      REAL*8,  INTENT(IN)  :: LND_FRC_DRY_SLC(HcoState%NX)
      REAL*8,  INTENT(IN)  :: ORO(HcoState%NX)
      REAL*8,  INTENT(IN)  :: SNW_FRC(HcoState%NX)
      REAL*8,  INTENT(IN)  :: TPT_SOI(HcoState%NX)
      REAL*8,  INTENT(IN)  :: TPT_SOI_FRZ
      REAL*8,  INTENT(IN)  :: VAI_DST_SLC(HcoState%NX)
      INTEGER, INTENT(OUT) :: MBL_NBR
      LOGICAL, INTENT(OUT) :: FLG_MBL(HcoState%NX)
      REAL*8,  INTENT(OUT) :: LND_FRC_MBL(HcoState%NX)
      INTEGER, INTENT(INOUT) :: RC

      !------------------
      ! Parameters
      !------------------

      ! VAI threshold quench [m2/m2]
      REAL*8,  PARAMETER   :: VAI_MBL_THR = 0.30D0

      !------------------
      ! Local variables
      !------------------

      ! [idx] Counting index
      INTEGER              :: IDX_IDX

      ! [idx] Interpolation month, future
      INTEGER              :: IDX_MTH_GLB

      ! [idx] Interpolation month, past
      INTEGER              :: IDX_MTH_LUB

      ! [idx] Longitude index array (land)
      INTEGER              :: LND_IDX(HcoState%NX)

      ! [nbr] Number of land points
      INTEGER              :: LND_NBR

      ! [idx] Counting index for longitude
      INTEGER              :: LON_IDX

      ! [idx] Surface type index
      INTEGER              :: SFC_TYP_IDX

      ! [idx] Surface sub-gridscale index
      INTEGER              :: SGS_IDX

      !-------------------------------------------------------------------
      ! Prior to 1/25/07:
      ! For GOCART source function, we don't use VAI (tdf, bmy, 1/25/07)
      !
      ! %%%%% DO NOT DELETE -- LEAVE THIS CODE COMMENTED OUT %%%%%
      !
      !! [flg] Use VAI data from time-varying boundary dataset
      ! LOGICAL              :: FLG_VAI_TVBDS = .TRUE.
      !-------------------------------------------------------------------

      ! For GOCART source function, we do not use VAI (tdf, bmy, 1/25/07)
      LOGICAL              :: FLG_VAI_TVBDS = .FALSE.

      ! [flg] Add 182 days in southern hemisphere
      LOGICAL              :: FLG_SH_ADJ = .TRUE.

      ! [dgr] Latitude
      REAL*8               :: LAT_DGR

      ! [m2 m-2] Leaf + stem area index, one-sided
      REAL*8               :: VAI_SGS

      ! Error handling
      CHARACTER(LEN=255)   :: MSG

      !=================================================================
      ! LND_FRC_MBL_GET begins here!
      !=================================================================

      ! Error check
      IF ( VAI_MBL_THR <= 0.0d0 ) THEN
         MSG = 'VAI_MBL_THR <= 0.0'
         CALL HCO_ERROR(HcoState%Config%Err,MSG, RC,
     &        THISLOC='LND_FRC_MBL_GET' )
         RETURN
      ENDIF

      ! Latitude (degrees)
      LAT_DGR = 180.0D0 * LAT_RDN/HcoState%Phys%PI

      ! Initialize outputs
      MBL_NBR = 0

      DO LON_IDX = 1, HcoState%NX
         FLG_MBL(LON_IDX) = .FALSE.
      ENDDO

      LND_FRC_MBL(:) = 0.0D0

      !=================================================================
      ! For dust mobilisation, we need to have land!  tdf 10/27/2K3
      ! Set up lnd_idx to hold the longitude indices for land
      ! Land ahoy!
      !=================================================================
      LND_NBR = 0
      DO LON_IDX = 1, HcoState%NX
         IF ( ORO_IS_LND( ORO(LON_IDX)) ) THEN
            LND_NBR          = LND_NBR + 1
            LND_IDX(LND_NBR) = LON_IDX
         ENDIF
      ENDDO

      ! Much ado about nothing (no land points)
      IF ( LND_NBR == 0 ) RETURN

!-----------------------------------------------------------------------------
! Prior to 1/25/07:
! When GOCART source function is used, VAI flag is NOT used, so
! we need to disable the ERROR_STOP call (tdf, bmy, 1/25/07)
!
! %%%%% DO NOT DELETE -- LEAVE THIS CODE COMMENTED OUT %%%%%
!
!      ! Introduce error message for flg_vai_tvbds=F (VAI not used!)
!      IF ( .not. FLG_VAI_TVBDS ) THEN
!c         print *,' FLG_VAI_TVBDS is false: GOCART source function used'
!         CALL ERROR_STOP( 'FLG_VAI_TVBDS=F',
!     &                    'LND_FRC_MBL_GET ("dust_dead_mod.f")' )
!      ENDIF
!-----------------------------------------------------------------------------

      !=================================================================
      ! Only land points are possible candidates for dust mobilization
      !=================================================================

      ! Loop over land points
      DO IDX_IDX = 1, LND_NBR
         LON_IDX = LND_IDX(IDX_IDX)

         ! Store surface blend of current gridpoint
         SFC_TYP_IDX = SFC_TYP_SLC(LON_IDX)

         ! Check for wet or frozen conditions - no mobilisation allowed
         ! Surface type 1  = inland lakes & land ice
         ! Surface type 27 = wetlands
         IF ( SFC_TYP_IDX <= 1  .OR. SFC_TYP_IDX >= 27 .OR.
     &        TPT_SOI(LON_IDX) < TPT_SOI_FRZ )          THEN

              ! SET bare ground fraction to zero
              LND_FRC_MBL(LON_IDX) = 0.0D0

         ELSE

           !-------------------------
           ! If we are using VAI...
           !-------------------------
           IF ( FLG_VAI_TVBDS ) THEN

              ! "bare ground" fraction of current gridcell decreases
              ! linearly from 1.0 to 0.0 as VAI increases from 0.0 to
              ! vai_mbl_thr.  NOTE: vai_mbl_thr set to 0.3  (tdf, 4/5/04)
              LND_FRC_MBL(LON_IDX) =
     &            1.0D0 - MIN(1.0D0, MIN(VAI_DST_SLC(LON_IDX),
     &                       VAI_MBL_THR) / VAI_MBL_THR)

           !---------------------------
           ! If we're not using VAI...
           !---------------------------
           ELSE

!-----------------------------------------------------------------------------
! Prior to 1/25/07:
! When GOCART source function is used, VAI flag is NOT used, so
! we need to disable the ERROR_STOP call. (tdf, bmy, 1/25/07)
!
! %%%%% DO NOT DELETE -- LEAVE THIS CODE COMMENTED OUT %%%%%
!
!              CALL ERROR_STOP( 'FLG_VAI_TVBDS=F',
!     &                         'LND_FRC_MBL_GET ("dust_dead_mod.f")' )
!-----------------------------------------------------------------------------

              ! For GOCART source function, set the bare
              ! ground fraction to 1 (tdf, bmy, 1/25/07)
              LND_FRC_MBL(LON_IDX) = 1.0D0

           ENDIF

         ENDIF                 ! endif normal land

         !==============================================================
         ! We have now filled "lnd_frc_mbl" the land fraction suitable
         ! for mobilisation.  Adjust for factors which constrain entire
         ! gridcell  LND_FRC_MBL modulated by LND_FRC_DRY and SNW_FRC.
         ! (tdf, 4/5/04)
         !==============================================================

         ! Take the bare ground fraction, multiply by the fraction
         ! that is dry and that is NOT covered by snow
         LND_FRC_MBL(LON_IDX) = LND_FRC_MBL(LON_IDX)
     &                        * LND_FRC_DRY_SLC(LON_IDX)
     &                        * ( 1.0D0 - SNW_FRC(LON_IDX) )

         ! Temporary fix for 1 x 1.25 grids -- Lok Lamsal 1/13/10
         IF ( LND_FRC_MBL(LON_IDX) .GT. 1.0D0 ) THEN
            LND_FRC_MBL(LON_IDX) = 0.99D0
         ENDIF

         ! Error check
         IF ( LND_FRC_MBL(lon_idx) > 1.0D0 ) THEN
            MSG = 'LND_FRC_MBL > 1'
            CALL HCO_ERROR(HcoState%Config%Err,MSG, RC,
     &         THISLOC='LND_FRC_MBL_GET' )
            RETURN
         ENDIF

         IF ( LND_FRC_MBL(LON_IDX) < 0.0D0 )   then
            MSG = 'LND_FRC_MBL < 0'
            CALL HCO_ERROR(HcoState%Config%Err,MSG, RC,
     &         THISLOC='LND_FRC_MBL_GET' )
            RETURN
         ENDIF

         ! If there is dry land in this longitude
         if ( LND_FRC_MBL(LON_IDX) > 0.0D0 ) then

            ! Set flag, we have a candidate!
            FLG_MBL(LON_IDX) = .TRUE.

            ! Increment # of candidates
            MBL_NBR          = MBL_NBR + 1
         ENDIF

      ENDDO

      ! Return w/ success
      RC = HCO_SUCCESS

      ! Return to calling program
      END SUBROUTINE LND_FRC_MBL_GET

!------------------------------------------------------------------------------

      SUBROUTINE DST_ADD_LON( NX, NBINS, Q, Q_TTL )
!
!******************************************************************************
!  Subroutine DST_ADD_LON dst_add_lon() computes and returns the total
!  property (e.g., mixing ratio, flux), obtained by simply adding along the
!  (dust) constituent dimension, when given an 3-D array of an additive
!  property (e.g., mixing ratio, flux). (tdf, bmy, 4/5/04)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) q     (REAL*8) : Total property
!
!  Arguments as Output:
!  ============================================================================
!  (2 ) q_ttl (REAL*8) : Property for each size class
!
!  NOTES:
!  (1 ) Updated comments, cosmetic changes.  Also now forces double-precision
!        with "D" exponents. (tdf, bmy, 4/5/04)
!******************************************************************************
!

      ! Arguments
      INTEGER, INTENT(IN) :: NX, NBINS
      REAL*8, INTENT(IN)  :: Q(NX,NBINS)
      REAL*8, INTENT(OUT) :: Q_TTL(NX)

      ! Local variables
      INTEGER             :: I, M

      !=================================================================
      ! DST_ADD_LON begins here!
      !=================================================================

      ! Initialize
      Q_TTL = 0d0

      ! Loop over dust bins
      DO M = 1, NBINS

         ! Loop over longitudes
         DO I = 1, NX

            ! Integrate!
            Q_TTL(I) = Q_TTL(I) + Q(I,M)

         ENDDO
      ENDDO

      ! Return to calling program
      END SUBROUTINE DST_ADD_LON

!------------------------------------------------------------------------------

      SUBROUTINE DST_TVBDS_GET( Inst, NX, LAT_IDX, VAI_DST_OUT )
!
!******************************************************************************
!  Subroutine DST_TVBDS_GET returns a specifed latitude slice of VAI data.
!  (tdf, bmy, 4/5/04)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) LAT_IDX     (INTEGER) : Latitude index
!
!  Arguments as Output:
!  ============================================================================
!  (2 ) VAI_DST_OUT (REAL*8 ) : Vegetation area index, 1-sided, current [m2/m2]
!
!  NOTES:
!  (1 ) Updated comments, cosmetic changes.  Also now forces double-precision
!        with "D" exponents. (tdf, bmy, 4/5/04)
!******************************************************************************
!

      ! Arguments
      TYPE(MyInst), POINTER  :: Inst
      INTEGER, INTENT(IN)    :: NX
      INTEGER, INTENT(IN)    :: LAT_IDX
      REAL*8,  INTENT(OUT)   :: VAI_DST_OUT(:)

      ! Local variables
      INTEGER              :: LON_IDX

      !=================================================================
      ! DST_TVBDS_GET begins here!
      !=================================================================

      ! Return lat slice of VAI [m2/m2]
      DO LON_IDX = 1, NX
         VAI_DST_OUT(LON_IDX) = Inst%VAI_DST(LON_IDX,LAT_IDX)
      ENDDO

      ! Return to calling program
      END SUBROUTINE DST_TVBDS_GET

!------------------------------------------------------------------------------

      SUBROUTINE OVR_SRC_SNK_FRC_GET( HcoState,
     &                                SRC_NBR,         MDN_SRC,
     &                                GSD_SRC,         SNK_NBR,
     &                                DMT_MIN_SNK,     DMT_MAX_SNK,
     &                                OVR_SRC_SNK_FRC, RC )

      USE HCO_CLOCK_MOD, ONLY : HcoClock_First
!
!******************************************************************************
!  Subroutine OVR_SRC_SNK_FRC_GET, given one set (the "source") of lognormal
!  distributions, and one set of bin boundaries (the "sink"), computes and
!  returns the overlap factors between the source distributions and the sink
!  bins.  (tdf, bmy, 4/5/04)
!
!  The output is a matrix, Mij, OVR_SRC_SNK_FRC(SRC_NBR,SNK_NBR)
!  Element ovr_src_snk_frc(i,j) is the fraction of size distribution i
!  in group src that overlaps sink bin j
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) SRC_NBR        (INTEGER)  : Dimension size                [unitless]
!  (2 ) MDN_SRC        (REAL*8 )  : Mass median particle size     [m       ]
!  (3 ) GSD_SRC        (REAL*8 )  : Geometric standard deviation  [fraction]
!  (4 ) SNK_NBR        (INTEGER)  : Dimension size                [unitless]
!  (5 ) DMT_MIN_SNK    (REAL*8 )  : Minimum diameter in bin       [m       ]
!  (6 ) DMT_MAX_SNK    (REAL*8 )  : Maximum diameter in bin       [m       ]
!
!  Arguments as Output:
!  ============================================================================
!  (7 ) OVR_SRC_SNK_FRC (REAL*8 ) : Fractional overlap of src with snk, Mij.
!
!  NOTES
!  (1 ) Updated comments, cosmetic changes.  Also now forces double-precision
!        with "D" exponents. (tdf, bmy, 4/5/04)
!******************************************************************************
!

      ! Arguments
      TYPE(HCO_State), POINTER :: HcoState
      INTEGER, INTENT(IN)      :: SRC_NBR
      REAL*8,  INTENT(IN)      :: MDN_SRC(SRC_NBR)
      REAL*8,  INTENT(IN)      :: GSD_SRC(SRC_NBR)
      INTEGER, INTENT(IN)      :: SNK_NBR
      REAL*8,  INTENT(IN)      :: DMT_MIN_SNK(SNK_NBR)
      REAL*8,  INTENT(IN)      :: DMT_MAX_SNK(SNK_NBR)
      REAL*8,  INTENT(OUT)     :: OVR_SRC_SNK_FRC(SRC_NBR,SNK_NBR)
      INTEGER, INTENT(INOUT)   :: RC

      ! Local
!      LOGICAL              :: FIRST = .TRUE.
      INTEGER              :: SRC_IDX         ! [idx] Counting index for src
      INTEGER              :: SNK_IDX         ! [idx] Counting index for snk
      REAL*8               :: LN_GSD          ! [frc] ln(gsd)
      REAL*8               :: SQRT2LNGSDI     ! [frc] Factor in erf() argument
      REAL*8               :: LNDMAXJOVRDMDNI ! [frc] Factor in erf() argument
      REAL*8               :: LNDMINJOVRDMDNI ! [frc] Factor in erf() argument
      CHARACTER(LEN=255)   :: MSG

      !=================================================================
      ! OVR_SRC_SNK_FRC_GET begins here
      !=================================================================

      IF ( HcoClock_First(HcoState%Clock,.TRUE.) ) THEN

         ! Test if ERF is implemented OK on this platform
         ! 19990913: erf() in SGI /usr/lib64/mips4/libftn.so is bogus
         IF ( ABS( 0.8427d0 - ERF(1.0d0) ) / 0.8427d0 > 0.001d0 ) THEN
            MSG = 'ERF error 1 in OVR_SRC_SNK_FRC_GET!'
            CALL HCO_ERROR(HcoState%Config%Err,MSG, RC,
     &         THISLOC='OVR_SRC_SNK_FRC_GET' )
            RETURN
         ENDIF

         ! Another ERF check
         IF ( ERF( 0.0D0 ) /= 0.0D0 ) THEN
            MSG = 'ERF error 2 in OVR_SRC_SNK_FRC_GET!'
            CALL HCO_ERROR(HcoState%Config%Err,MSG, RC,
     &         THISLOC='OVR_SRC_SNK_FRC_GET' )
            RETURN
         ENDIF

         ! Reset first-time flag
         !FIRST = .FALSE.
      ENDIF


      ! Loop over source index (cf Zender et al eq 12)
      DO SRC_IDX = 1, SRC_NBR

         ! Fraction
         SQRT2LNGSDI = SQRT(2.0D0) * LOG( GSD_SRC(SRC_IDX) )

         ! Loop over sink index
         DO SNK_IDX = 1, SNK_NBR

            ! [fraction]
            LNDMAXJOVRDMDNI = LOG(DMT_MAX_SNK(SNK_IDX)/MDN_SRC(SRC_IDX))

            ! [fraction]
            LNDMINJOVRDMDNI = LOG(DMT_MIN_SNK(SNK_IDX)/MDN_SRC(SRC_IDX))

            ! [fraction]
            OVR_SRC_SNK_FRC (SRC_IDX,SNK_IDX)=  ! [frc]
     &            0.5D0 * (ERF(LNDMAXJOVRDMDNI/SQRT2LNGSDI)
     &                   - ERF(LNDMINJOVRDMDNI/SQRT2LNGSDI) )
         ENDDO
      ENDDO

      ! Return w/ success
      RC = HCO_SUCCESS

      END SUBROUTINE OVR_SRC_SNK_FRC_GET

!------------------------------------------------------------------------------

       FUNCTION ERF( X ) RESULT( ERF_VAL )
!
!******************************************************************************
!  Function ERF returns the error function erf(x).  See comments heading
!  routine CALERF below.  Author/Date: W. J. Cody, January 8, 1985
!  (tdf, bmy, 4/5/04)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) X (REAL*8) : Argument to erf(x)
!
!  NOTES:
!  (1 ) Updated comments (bmy, 4/5/04)
!******************************************************************************
!
       IMPLICIT NONE

       ! Arguments
       REAL*8, INTENT(IN) :: X

       ! Local variables
       INTEGER            :: JINT
       REAL*8             :: RESULT, ERF_VAL

       !================================================================
       ! ERF begins here!
       !================================================================
       JINT = 0
       CALL CALERF( X, RESULT, JINT )
       ERF_VAL = RESULT

       ! Return to calling program
       END FUNCTION ERF

!------------------------------------------------------------------------------

       SUBROUTINE CALERF( ARG, RESULT, JINT )
!
!******************************************************************************
!  This packet evaluates erf(x), erfc(x), and exp(x*x)*erfc(x)
!  for a real argument  x.  It contains three function type
!  subprograms: erf, erfc, and erfcx (or derf, derfc, and derfcx),
!  and one subroutine type subprogram, calerf.  The calling
!  statements for the primary entries are:
!
!  y=erf(x)     (or   y=derf(x)),
!  y=erfc(x)    (or   y=derfc(x)),
!  and
!  y=erfcx(x)   (or   y=derfcx(x)).
!
!  The routine  calerf  is intended for internal packet use only,
!  all computations within the packet being concentrated in this
!  routine.  The function subprograms invoke  calerf  with the
!  statement
!  call calerf(arg,result,jint)
!  where the parameter usage is as follows
!
!  Function                     Parameters for calerf
!  Call              Arg                  Result          Jint
!
!  erf(arg)      any real argument         erf(arg)          0
!  erfc(arg)     abs(arg)  <  xbig        erfc(arg)          1
!  erfcx(arg)    xneg  <  arg  <  xmax   erfcx(arg)          2
!
!  The main computation evaluates near-minimax approximations:
!  from "Rational Chebyshev Approximations for the Error Function"
!  by W. J. Cody, Math. Comp., 1969, pp. 631-638.  This
!  transportable program uses rational functions that theoretically
!  approximate  erf(x)  and  erfc(x)  to at least 18 significant
!  decimal digits.  The accuracy achieved depends on the arithmetic
!  system, the compiler, the intrinsic functions, and proper
!  selection of the machine-dependent constants.
!
!  Explanation of machine-dependent constants:
!  xmin   = The smallest positive floating-point number.
!  xinf   = The largest positive finite floating-point number.
!  xneg   = The largest negative argument acceptable to erfcx;
!  the negative of the solution to the equation
!  2*exp(x*x) = xinf.
!  xsmall = Argument below which erf(x) may be represented by
!  2*x/sqrt(pi)  and above which  x*x  will not underflow.
!  A conservative value is the largest machine number x
!  such that   1.0 + x = 1.0   to machine precision.
!  xbig   = Largest argument acceptable to erfc;  solution to
!  the equation:  w(x)* (1-0.5/x**2) = xmin,  where
!  w(x) = exp(-x*x)/[x*sqrt(pi)].
!  xhuge  = Argument above which  1.0 - 1/(2*x*x) = 1.0  to
!  machine precision.  a conservative value is
!  1/[2*sqrt(xsmall)]
!  xmax   = Largest acceptable argument to erfcx; the minimum
!  of xinf and 1/[sqrt(pi)*xmin].
!
!  Approximate values for some important machines are:
!  xmin       xinf        xneg     xsmall
!  CDC 7600      (s.p.)  3.13e-294   1.26e+322   -27.220  7.11e-15
!  Cray-1        (s.p.)  4.58e-2467  5.45e+2465  -75.345  7.11e-15
!  IEEE (IBM/XT,
!  Sun, etc.)  (s.p.)  1.18e-38    3.40e+38     -9.382  5.96e-8
!  IEEE (IBM/XT,
!  Sun, etc.)  (d.p.)  2.23d-308   1.79d+308   -26.628  1.11d-16
!  IBM 195       (d.p.)  5.40d-79    7.23e+75    -13.190  1.39d-17
!  Univac 1108   (d.p.)  2.78d-309   8.98d+307   -26.615  1.73d-18
!  Vax d-format  (d.p.)  2.94d-39    1.70d+38     -9.345  1.39d-17
!  Vax g-format  (d.p.)  5.56d-309   8.98d+307   -26.615  1.11d-16
!
!  xbig       xhuge       xmax
!  CDC 7600      (s.p.)  25.922      8.39e+6     1.80x+293
!  Cray-1        (s.p.)  75.326      8.39e+6     5.45e+2465
!  IEEE (IBM/XT,
!  Sun, etc.)  (s.p.)   9.194      2.90e+3     4.79e+37
!  IEEE (IBM/XT,
!  Sun, etc.)  (d.p.)  26.543      6.71d+7     2.53d+307
!  IBM 195       (d.p.)  13.306      1.90d+8     7.23e+75
!  Univac 1108   (d.p.)  26.582      5.37d+8     8.98d+307
!  Vax d-format  (d.p.)   9.269      1.90d+8     1.70d+38
!  Vax g-format  (d.p.)  26.569      6.71d+7     8.98d+307
!
!  Error returns:
!  The program returns  erfc = 0      for  arg  >=  xbig;
!  erfcx = xinf  for  arg  <  xneg;
!  and
!  erfcx = 0     for  arg  >=  xmax.
!
!  Intrinsic functions required are:
!  abs, aint, exp
!
!  Author: W. J. Cody
!  Mathematics And Computer Science Division
!  Argonne National Laboratory
!  Argonne, IL 60439
!  Latest modification: March 19, 1990
!
!  NOTES:
!  (1 ) Now force double-precision w/ "D" exponents (bmy, 4/5/04)
!******************************************************************************
!
      IMPLICIT NONE
      INTEGER I,JINT
      REAL*8  A,ARG,B,C,D,DEL,FOUR,HALF,P,ONE,Q,RESULT,SIXTEN,SQRPI,
     &   TWO,THRESH,X,XBIG,XDEN,XHUGE,XINF,XMAX,XNEG,XNUM,XSMALL,
     &   Y,YSQ,ZERO
      DIMENSION A(5),B(4),C(9),D(8),P(6),Q(5)

      ! Mathematical constants
      data four,one,half,two,zero/4.0d0,1.0d0,0.5d0,2.0d0,0.0d0/,
     &         sqrpi/5.6418958354775628695d-1/,thresh/0.46875d0/,
     &         sixten/16.0d0/

      ! Machine-dependent constants
      data xinf,xneg,xsmall/3.40d+38,-9.382d0,5.96d-8/,
     &      xbig,xhuge,xmax/9.194d0,2.90d3,4.79d37/

      ! Coefficients for approximation to  erf  in first interval
      data a /3.16112374387056560d00,1.13864154151050156d02,
     &     3.77485237685302021d02,3.20937758913846947d03,
     &     1.85777706184603153d-1/

      data b /2.36012909523441209d01,2.44024637934444173d02,
     &     1.28261652607737228d03,2.84423683343917062d03/

      ! Coefficients for approximation to  erfc  in second interval
      data c /5.64188496988670089d-1,8.88314979438837594d0,
     &     6.61191906371416295d01,2.98635138197400131d02,
     &     8.81952221241769090d02,1.71204761263407058d03,
     &     2.05107837782607147d03,1.23033935479799725d03,
     &     2.15311535474403846d-8/

      data d /1.57449261107098347d01,1.17693950891312499d02,
     &     5.37181101862009858d02,1.62138957456669019d03,
     &     3.29079923573345963d03,4.36261909014324716d03,
     &     3.43936767414372164d03,1.23033935480374942d03/

      ! Coefficients for approximation to  erfc  in third interval
      data p /3.05326634961232344d-1,3.60344899949804439d-1,
     &     1.25781726111229246d-1,1.60837851487422766d-2,
     &     6.58749161529837803d-4,1.63153871373020978d-2/

      data q /2.56852019228982242d00,1.87295284992346047d00,
     &     5.27905102951428412d-1,6.05183413124413191d-2,
     &     2.33520497626869185d-3/

c Main Code
      x=arg
      y=abs(x)
      if (y <= thresh) then
c Evaluate  erf  for  |x| <= 0.46875
        ysq=zero
        if (y > xsmall) ysq=y*y
        xnum=a(5)*ysq
        xden=ysq
        do i=1,3
          xnum=(xnum+a(i))*ysq
          xden=(xden+b(i))*ysq
        end do
        result=x*(xnum+a(4))/(xden+b(4))
        if (jint /= 0) result=one-result
        if (jint == 2) result=exp(ysq)*result
        go to 800

c Evaluate  erfc  for 0.46875 <= |x| <= 4.0
      else if (y <= four) then
        xnum=c(9)*y
        xden=y
        do i=1,7
          xnum=(xnum+c(i))*y
          xden=(xden+d(i))*y
        end do
        result=(xnum+c(8))/(xden+d(8))
        if (jint /= 2) then
          ysq=aint(y*sixten)/sixten
          del=(y-ysq)*(y+ysq)
          result=exp(-ysq*ysq)*exp(-del)*result
        end if

c Evaluate  erfc  for |x| > 4.0
      else
        result=zero
        if (y >= xbig) then
          if ((jint /= 2).or.(y >= xmax)) go to 300
          if (y >= xhuge) then
             result=sqrpi/y
             go to 300
          end if
        end if
        ysq=one/(y*y)
        xnum=p(6)*ysq
        xden=ysq
        do i=1,4
          xnum=(xnum+p(i))*ysq
          xden=(xden+q(i))*ysq
        end do
        result=ysq*(xnum+p(5))/(xden+q(5))
        result=(sqrpi-result)/y
        if (jint /= 2) then
          ysq=aint(y*sixten)/sixten
          del=(y-ysq)*(y+ysq)
          result=exp(-ysq*ysq)*exp(-del)*result
        end if
      end if

c Fix up for negative argument, erf, etc.
  300 if (jint == 0) then
        result=(half-result)+half
        if (x < zero) result=-result
      else if (jint == 1) then
        if (x < zero) result=two-result
      else
        if (x < zero) then
          if (x < xneg) then
             result=xinf
          else
             ysq=aint(x*sixten)/sixten
             del=(x-ysq)*(x+ysq)
             y=exp(ysq*ysq)*exp(del)
             result=(y+y)-result
          end if
        end if
      end if
  800 return

      ! Return to calling program
      END SUBROUTINE CALERF

!------------------------------------------------------------------------------

      SUBROUTINE PLN_TYP_GET( PLN_TYP, PLN_FRC, TAI )

!
!******************************************************************************
!  Subroutine PLN_TYPE_GET returns LSM information needed by the DEAD
!  dust parameterization. (tdf, bmy, 4/5/04)
!
!  Arguments as Output:
!  ============================================================================
!  (1 ) PLN_TYP (INTEGER) : LSM plant type index (1..14)
!  (2 ) PLN_TYP (REAL*8 ) : Weight of corresponding plant type (sums to 1.0)
!  (3 ) TAI     (REAL*8 ) : Leaf-area index (one sided) [index]
!
!  NOTES:
!  (1 ) Updated comments.  Now force double-precision w/ "D" exponents.
!        (bmy, 4/5/04)
!******************************************************************************
!
      ! Arguments
      INTEGER, INTENT(OUT) :: PLN_TYP(0:28,3)
      REAL*8,  INTENT(OUT) :: PLN_FRC(0:28,3)
      REAL*8,  INTENT(OUT) :: TAI(14,12)

      ! Local variables
      INTEGER              :: I, J

      !=================================================================
      ! There are 29 land surface types: 0 = ocean, 1 to 28 = land.
      ! Each land point has up to three vegetation types, ranging in
      ! value from 1 to 14.  PLN_TYPE contains the vegetation type of
      ! the 3 subgrid points for each surface type.  PLN_FRC contains
      ! the fractional area of the 3 subgrid points for each surface
      ! type.
      !=================================================================
      PLN_TYP(0:28,1) = (/   0,
     &                      14,  14,   1,   2,   4,   1  , 1,
     &                       4,   1,   3,   5,  13,   1,   2,
     &                      11,  11,   6,  13,   9,   7,   8,
     &                       8,  12,  11,  12,  11,   3,  14/)

      PLN_FRC(0:28,1) = (/ 0.00d0,
     &                     1.00d0, 1.00d0, 0.75d0, 0.50d0,
     &                     0.75d0, 0.37d0, 0.75d0,
     &                     0.75d0, 0.37d0, 0.95d0, 0.75d0,
     &                     0.70d0, 0.25d0, 0.25d0,
     &                     0.40d0, 0.40d0, 0.60d0, 0.60d0,
     &                     0.30d0, 0.80d0, 0.80d0,
     &                     0.10d0, 0.85d0, 0.85d0, 0.85d0,
     &                     0.85d0, 0.80d0, 1.00d0/)


      PLN_TYP(0:28,2) = (/   0,
     &                      14,  14,  14,  14,  14,   4  ,14,
     &                      14,   4,  14,  14,   5,  10,  10,
     &                       4,   4,  13,   6,  10,  14,  14,
     &                      14,  14,  14,  14,  14,  14,  14/)

      PLN_FRC(0:28,2) = (/ 0.00d0,
     &                     0.00d0, 0.00d0, 0.25d0, 0.50d0,
     &                     0.25d0, 0.37d0, 0.25d0,
     &                     0.25d0, 0.37d0, 0.05d0, 0.25d0,
     &                     0.30d0, 0.25d0, 0.25d0,
     &                     0.30d0, 0.30d0, 0.20d0, 0.20d0,
     &                     0.30d0, 0.20d0, 0.20d0,
     &                     0.90d0, 0.15d0, 0.15d0, 0.15d0,
     &                     0.15d0, 0.20d0, 0.00d0/)

      PLN_TYP(0:28,3) = (/   0,
     &                      14,  14,  14,  14,  14,  14,  14,
     &                      14,  14,  14,  14,  14,  14,  14,
     &                       1,   1,  14,  14,  14,  14,  14,
     &                      14,  14,  14,  14,  14,  14,  14/)

      PLN_FRC(0:28,3) = (/ 0.00d0,
     &                     0.00d0, 0.00d0, 0.00d0, 0.00d0,
     &                     0.00d0, 0.26d0, 0.00d0,
     &                     0.00d0, 0.26d0, 0.00d0, 0.00d0,
     &                     0.00d0, 0.50d0, 0.50d0,
     &                     0.30d0, 0.30d0, 0.20d0, 0.20d0,
     &                     0.40d0, 0.00d0, 0.00d0,
     &                     0.00d0, 0.00d0, 0.00d0, 0.00d0,
     &                     0.00d0, 0.00d0, 0.00d0/)

      !=================================================================
      ! ----------------------------------------------------------------
      ! description of the 29 surface types
      ! ----------------------------------------------------------------
      !
      ! no vegetation
      ! -------------
      !  0 ocean
      !  1 land ice (glacier)
      !  2 desert
      !
      ! forest vegetation
      ! -----------------
      !  3 cool needleleaf evergreen tree
      !  4 cool needleleaf deciduous tree
      !  5 cool broadleaf  deciduous tree
      !  6 cool mixed needleleaf evergreen and broadleaf deciduous tree
      !  7 warm needleleaf evergreen tree
      !  8 warm broadleaf  deciduous tree
      !  9 warm mixed needleleaf evergreen and broadleaf deciduous tree
      ! 10 tropical broadleaf evergreen tree
      ! 11 tropical seasonal deciduous tree
      !
      ! interrupted woods
      ! ----------------
      ! 12 savanna
      ! 13 evergreen forest tundra
      ! 14 deciduous forest tundra
      ! 15 cool forest crop
      ! 16 warm forest crop
      !
      ! non-woods
      ! ---------
      ! 17 cool grassland
      ! 18 warm grassland
      ! 19 tundra
      ! 20 evergreen shrub
      ! 21 deciduous shrub
      ! 22 semi-desert
      ! 23 cool irrigated crop
      ! 24 cool non-irrigated crop
      ! 25 warm irrigated crop
      ! 26 warm non-irrigated crop
      !
      ! wetlands
      ! --------
      ! 27 forest (mangrove)
      ! 28 non-forest
      !
      ! ----------------------------------------------------------------
      ! description of the 14 plant types. see vegconi.F for
      ! parameters that depend on vegetation type
      ! ----------------------------------------------------------------
      !
      !  1 = needleleaf evergreen tree
      !  2 = needleleaf deciduous tree
      !  3 = broadleaf evergreen tree
      !  4 = broadleaf deciduous tree
      !  5 = tropical seasonal tree
      !  6 = cool grass (c3)
      !  7 = evergreen shrub
      !  8 = deciduous shrub
      !  9 = arctic deciduous shrub
      ! 10 = arctic grass
      ! 11 = crop
      ! 12 = irrigated crop
      ! 13 = warm grass (c4)
      ! 14 = not vegetated
      !=================================================================

      ! TAI = monthly leaf area index + stem area index, one-sided
      TAI(1,1:12) =  (/ 4.5d0, 4.7d0, 5.0d0, 5.1d0, 5.3d0, 5.5d0,
     &                  5.3d0, 5.3d0, 5.2d0, 4.9d0, 4.6d0, 4.5d0 /)

      TAI(2,1:12) =  (/ 0.3d0, 0.3d0, 0.3d0, 1.0d0, 1.6d0, 2.4d0,
     &                  4.3d0, 2.9d0, 2.0d0, 1.3d0, 0.8d0, 0.5d0 /)

      TAI(3,1:12) =  (/ 5.0d0, 5.0d0, 5.0d0, 5.0d0, 5.0d0, 5.0d0,
     &                  5.0d0, 5.0d0, 5.0d0, 5.0d0, 5.0d0, 5.0d0 /)

      TAI(4,1:12) =  (/ 0.4d0, 0.4d0, 0.7d0, 1.6d0, 3.5d0, 5.1d0,
     &                  5.4d0, 4.8d0, 3.8d0, 1.7d0, 0.6d0, 0.4d0 /)

      TAI(5,1:12) =  (/ 1.2d0, 1.0d0, 0.9d0, 0.8d0, 0.8d0, 1.0d0,
     &                  2.0d0, 3.7d0, 3.2d0, 2.7d0, 1.9d0, 1.2d0 /)

      TAI(6,1:12) =  (/ 0.7d0, 0.8d0, 0.9d0, 1.0d0, 1.5d0, 3.4d0,
     &                  4.3d0, 3.8d0, 1.8d0, 1.0d0, 0.9d0, 0.8d0 /)

      TAI(7,1:12) =  (/ 1.3d0, 1.3d0, 1.3d0, 1.3d0, 1.3d0, 1.3d0,
     &                  1.3d0, 1.3d0, 1.3d0, 1.3d0, 1.3d0, 1.3d0 /)

      TAI(8,1:12) =  (/ 1.0d0, 1.0d0, 0.8d0, 0.3d0, 0.6d0, 0.0d0,
     &                  0.1d0, 0.3d0, 0.5d0, 0.6d0, 0.7d0, 0.9d0 /)

      TAI(9,1:12) =  (/ 0.1d0, 0.1d0, 0.1d0, 0.1d0, 0.1d0, 0.3d0,
     &                  1.5d0, 1.7d0, 1.4d0, 0.1d0, 0.1d0, 0.1d0 /)

      TAI(10,1:12) = (/ 0.7d0, 0.8d0, 0.9d0, 1.0d0, 1.5d0, 3.4d0,
     &                  4.3d0, 3.8d0, 1.8d0, 1.0d0, 0.9d0, 0.8d0 /)

      TAI(11,1:12) = (/ 0.0d0, 0.0d0, 0.0d0, 0.0d0, 1.0d0, 2.0d0,
     &                  3.0d0, 3.0d0, 1.5d0, 0.0d0, 0.0d0, 0.0d0 /)

      TAI(12,1:12) = (/ 0.0d0, 0.0d0, 0.0d0, 0.0d0, 1.0d0, 2.0d0,
     &                  3.0d0, 3.0d0, 1.5d0, 0.0d0, 0.0d0, 0.0d0 /)

      TAI(13,1:12) = (/ 0.7d0, 0.8d0, 0.9d0, 1.0d0, 1.5d0, 3.4d0,
     &                  4.3d0, 3.8d0, 1.8d0, 1.0d0, 0.9d0, 0.8d0 /)

      TAI(14,1:12) = (/ 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0,
     &                  0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0 /)

      ! Return to calling program
      END SUBROUTINE PLN_TYP_GET

!******************************************************************************

      FUNCTION HCOX_DustDead_GetFluxTun( HcoState )
     & RESULT( FLX_MSS_FDG_FCT )
!
!******************************************************************************
! Returns global mass flux tuning factor (a posteriori).
!
!******************************************************************************
!
      ! Input argument
      TYPE(HCO_State), POINTER :: HcoState    ! Hemco state

      ! Output argument
      REAL*8  :: FLX_MSS_FDG_FCT
      REAL*8  :: DY

      ! Local variables
      CHARACTER(LEN=255)  :: MSG
      INTEGER             :: RC

      !=================================================================
      ! HCOX_DustDead_GetFluxTun begins here!
      !=================================================================

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Prior to 4/13/19:
! Need to set tuning factor in HEMCO_Config.rc for FlexGrid. Keep factors for
! 2x2.5 and 4x5 for now. (mps, 4/13/19)
!      ! Global mass flux tuning factor (a posteriori) [frc]
!#if   defined( GEOS_FP ) && defined( GRID025x03125 )
!
!#if defined(NESTED_CH)
!
!      !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!      !%%%       NOTE: NEED TO SCALE THIS FOR GEOS-FP nested CH!         %%%
!      !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!      ! Copy the 05x0666 scale factors for now (skim, 5/20/12)
!      FLX_MSS_FDG_FCT = 3.23d-4
!
!#elif defined(NESTED_NA)
!
!      !----------------------------------------------------------------
!      ! Based on results from GEOS-FP 0.25x0.3125 NA for 08/2012-08/2013
!      ! compared to 2006 GEOS-5 0.5x0.666 NA
!      !
!      !   (GEOS-FP - GEOS-5)/GEOS-5 * 100 is 80.04% in each size bin.
!      !
!      ! We need to scale to the parameter FLX_MSS_FDG_FCT to make the
!      ! dust emissions consistent. Consequently, to bring GEOS-FP nested
!      ! NA dust emissions down to GEOS-5 nested NA levels, we need use
!      ! a FLX_MSS_FDG_FCT of:
!      !
!      !             1 / (1 + 0.8004) = 0.5554
!      !
!      !             2.16d-4 * 0.5554 = 1.20d-4   ( 2.16d-4 is factor for
!      !                                            GEOS-5 nested NA )
!      !
!      !    -- M. Sulprizio, 15 Jan 2014
!      !----------------------------------------------------------------
!      FLX_MSS_FDG_FCT = 1.20d-4
!
!#endif
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      IF ( TRIM(HcoState%Config%MetField) == 'GEOSFP'   .and.
     &     TRIM(HcoState%Config%GridRes)  == '2.0x2.5'  ) THEN

      !----------------------------------------------------------------
      ! Based on results from GEOS-FP 2x25 for 08/2012-08/2013
      !
      !   (GEOS-FP - GEOS-5)/GEOS-5 * 100 is -2.810% in each size bin.
      !
      ! We need to scale to the parameter FLX_MSS_FDG_FCT to make the
      ! dust emissions consistent.  Consequently, to bring GEOS-FP 2x25
      ! dust emissions up to 2005 GEOS-5 levels, we need to MULTIPLY the
      ! FLX_MSS_FDG_FCT used for GEOS-5 by:
      !
      !             1 / (1. - 0.02810) = 1 / 0.9719 = 1.0289
      !
      !    -- M. Sulprizio, 15 Jan 2014
      !----------------------------------------------------------------
      FLX_MSS_FDG_FCT = 4.9d-4 * 1.0289d0

      ELSE IF ( TRIM(HcoState%Config%MetField) == 'GEOSFP'   .and.
     &          TRIM(HcoState%Config%GridRes)  == '4.0x5.0'  ) THEN

      !----------------------------------------------------------------
      ! Based on results from GEOS-FP 4x5 for 08/2012-08/2013
      !
      !   (GEOS-FP - GEOS-5)/GEOS-5 * 100 is -15.95% in each size bin.
      !
      ! We need to scale to the parameter FLX_MSS_FDG_FCT to make the
      ! dust emissions consistent.  Consequently, to bring GEOS-FP 4x5
      ! dust emissions up to 2005 GEOS-5 levels, we need to MULTIPLY the
      ! FLX_MSS_FDG_FCT used for GEOS-5 by:
      !
      !             1 / (1. - 0.1595) = 1 / 0.8405 = 1.1898
      !
      !    -- R. Yantosca, M. Sulprizio, 06 Jan 2014
      !----------------------------------------------------------------
      FLX_MSS_FDG_FCT = 7.0d-4 * 1.1898d0

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Prior to 4/13/19:
! Need to set tuning factor in HEMCO_Config.rc for FlexGrid. Keep factors for
! 2x2.5 and 4x5 for now. (mps, 4/13/19)
!#elif defined( MERRA2 ) && defined( GRID05x0625 )
!
!#if defined(NESTED_AS)
!
!      !----------------------------------------------------------------
!      ! Based on results from MERRA2 0.5x0.625 CH for 2011 compared to
!      ! GEOS-5 0.5x0.666 CH for 2011
!      !
!      !   (MERRA2 - GEOS-5)/GEOS-5 * 100 is 5.788% in each size bin.
!      !
!      ! We need to scale to the parameter FLX_MSS_FDG_FCT to make the
!      ! dust emissions consistent. Consequently, to bring MERRA2 nested
!      ! CH dust emissions down to GEOS-5 nested CH levels, we need use
!      ! a FLX_MSS_FDG_FCT of:
!      !
!      !             1 / (1 + 0.05788) = 0.9453
!      !
!      !             3.23d-4 * 0.9453 = 3.05d-4   ( 3.23d-4 is factor for
!      !                                            GEOS-5 nested CH )
!      !
!      ! NOTE: The MERRA2 0.5x0.625 CH grid extends further west (to 60E)
!      ! than the GEOS-5 0.5x0.666 CH grid (to 70E). For this calculation,
!      ! any dust emissions in the 60E-70E band in the MERRA2 CH output
!      ! were ignored.
!      !
!      !    -- M. Sulprizio, 22 Feb 2016
!      !----------------------------------------------------------------
!      FLX_MSS_FDG_FCT = 3.05d-4
!
!#elif defined(NESTED_NA)
!
!      !----------------------------------------------------------------
!      ! Based on results from MERRA2 0.5x0.625 NA for 2012 compared to
!      ! GEOS-5 0.5x0.666 CH for 2012
!      !
!      !   (MERRA2 - GEOS-5)/GEOS-5 * 100 is 25.02% in each size bin.
!      !
!      ! We need to scale to the parameter FLX_MSS_FDG_FCT to make the
!      ! dust emissions consistent. Consequently, to bring MERRA2 nested
!      ! NA dust emissions down to GEOS-5 nested NA levels, we need use
!      ! a FLX_MSS_FDG_FCT of:
!      !
!      !             1 / (1 + 0.2502) = 0.7999
!      !
!      !             2.16d-4 * 0.7999 = 1.73d-4   ( 2.16d-4 is factor for
!      !                                            GEOS-5 nested NA )
!      !
!      !    -- M. Sulprizio, 22 Feb 2016
!      !----------------------------------------------------------------
!      FLX_MSS_FDG_FCT = 1.73d-4
!
!#endif
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      ELSE IF ( TRIM(HcoState%Config%MetField) == 'MERRA2'   .and.
     &          TRIM(HcoState%Config%GridRes)  == '2.0x2.5'  ) THEN

      !----------------------------------------------------------------
      ! Based on results from MERRA2 2x2.5 for 2012 (Chi Li, 7/8/16)
      !----------------------------------------------------------------
      FLX_MSS_FDG_FCT = 4.9d-4 * 0.971141

      ELSE IF ( TRIM(HcoState%Config%MetField) == 'MERRA2'   .and.
     &          TRIM(HcoState%Config%GridRes)  == '4.0x5.0'    ) THEN

      !----------------------------------------------------------------
      ! Based on results from MERRA2 4x5 for 2012
      !
      !   (MERRA2 - GEOS-5)/GEOS-5 * 100 is -10.86% in each size bin.
      !
      ! We need to scale to the parameter FLX_MSS_FDG_FCT to make the
      ! dust emissions consistent.  Consequently, to bring 2012 MERRA2
      ! dust emissions up to 2012 GEOS-5 levels, we need to MULTIPLY the
      ! FLX_MSS_FDG_FCT used for GEOS-5 by:
      !
      !             1 / (1. - 0.1086) = 1 / 0.8913 = 1.1219
      !
      !    -- M. Sulprizio, 12 Feb 2016
      !----------------------------------------------------------------
      FLX_MSS_FDG_FCT = 7.0d-4 * 1.1219d0

      ELSE

      ! Default value
      ! ==> Don't allow default value anymore. The tuning factor can
      ! now also be set in the HEMCO configuration file. This should be
      ! done whenever the compiler combinations above don't apply!
!      FLX_MSS_FDG_FCT = 7.0d-4
      FLX_MSS_FDG_FCT = -999.0d0

      ENDIF

      ! If not defined yet, try to estimate it from grid spacing
      ! ckeller, 6/6/2017
      IF ( FLX_MSS_FDG_FCT == -999.0d0 ) THEN
         ! average latitude spacing
         DY = ABS(MAXVAL(HcoState%Grid%YMID%Val)
     &          - MINVAL(HcoState%Grid%YMID%Val)) / ( HcoState%NY - 1 )

         ! 0.125 degrees / C720
         IF ( DY < 0.175d0 ) THEN
            FLX_MSS_FDG_FCT = 2.3d-4
         ! 0.25 degrees / C360
         ELSEIF ( DY < 0.375d0 ) THEN
            FLX_MSS_FDG_FCT = 2.35d-4
         ! 0.5 degrees / C180
         ELSEIF ( DY < 0.75d0 ) THEN
            FLX_MSS_FDG_FCT = 3.23d-4
         ! 1.0 degrees / C90
         ELSEIF ( DY < 1.5d0 ) THEN
            FLX_MSS_FDG_FCT = 4.5d-4
         ! 2.0 degrees / C48
         ELSEIF ( DY < 3.0d0 ) THEN
            FLX_MSS_FDG_FCT = 5.0416d-4
         ENDIF
      ENDIF

      ! Error
      IF ( FLX_MSS_FDG_FCT == -999.0d0 ) THEN
         MSG = 'Cannot assign mass flux tuning factor - ' //
     &         'please explicitly set it by adding ' //
     &         '` --> Mass tuning factor: XX.X` to the extension ' //
     &         'options in your config. file.'
            CALL HCO_ERROR(HcoState%Config%Err, MSG,
     &                     RC, THISLOC='HCOX_DustDead_GetFluxTun')
         RETURN
      ENDIF

      ! Return to calling program
      END FUNCTION HCOX_DustDead_GetFluxTun

!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: InstGet
!
! !DESCRIPTION: Subroutine InstGet returns a pointer to the desired instance.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE InstGet ( Instance, Inst, RC, PrevInst )
!
! !INPUT PARAMETERS:
!
      INTEGER                             :: Instance
      TYPE(MyInst),     POINTER           :: Inst
      INTEGER                             :: RC
      TYPE(MyInst),     POINTER, OPTIONAL :: PrevInst
!
! !REVISION HISTORY:
!  18 Feb 2016 - C. Keller   - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
      TYPE(MyInst),     POINTER    :: PrvInst

      !=================================================================
      ! InstGet begins here!
      !=================================================================

      ! Get instance. Also archive previous instance.
      PrvInst => NULL()
      Inst    => AllInst
      DO WHILE ( ASSOCIATED(Inst) )
         IF ( Inst%Instance == Instance ) EXIT
         PrvInst => Inst
         Inst    => Inst%NextInst
      END DO
      IF ( .NOT. ASSOCIATED( Inst ) ) THEN
         RC = HCO_FAIL
         RETURN
      ENDIF

      ! Pass output arguments
      IF ( PRESENT(PrevInst) ) PrevInst => PrvInst

      ! Cleanup & Return
      PrvInst => NULL()
      RC = HCO_SUCCESS

      END SUBROUTINE InstGet
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: InstCreate
!
! !DESCRIPTION: Subroutine InstCreate creates a new instance.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE InstCreate ( ExtNr, Instance, Inst, RC )
!
! !INPUT PARAMETERS:
!
      INTEGER,       INTENT(IN)       :: ExtNr
!
! !OUTPUT PARAMETERS:
!
      INTEGER,       INTENT(  OUT)    :: Instance
      TYPE(MyInst),  POINTER          :: Inst
!
! !INPUT/OUTPUT PARAMETERS:
!
      INTEGER,       INTENT(INOUT)    :: RC
!
! !REVISION HISTORY:
!  18 Feb 2016 - C. Keller   - Initial version
!  26 Oct 2016 - R. Yantosca - Don't nullify local ptrs in declaration stmts
!EOP
!------------------------------------------------------------------------------
!BOC
      TYPE(MyInst), POINTER          :: TmpInst
      INTEGER                        :: nnInst

      !=================================================================
      ! InstCreate begins here!
      !=================================================================

      ! ----------------------------------------------------------------
      ! Generic instance initialization
      ! ----------------------------------------------------------------

      ! Initialize
      Inst => NULL()

      ! Get number of already existing instances
      TmpInst => AllInst
      nnInst = 0
      DO WHILE ( ASSOCIATED(TmpInst) )
         nnInst  =  nnInst + 1
         TmpInst => TmpInst%NextInst
      END DO

      ! Create new instance
      ALLOCATE(Inst)
      Inst%Instance = nnInst + 1
      Inst%ExtNr    = ExtNr

      ! Attach to instance list
      Inst%NextInst => AllInst
      AllInst       => Inst

      ! Update output instance
      Instance = Inst%Instance

      ! ----------------------------------------------------------------
      ! Type specific initialization statements follow below
      ! ----------------------------------------------------------------
      Inst%ERD_FCT_GEO     => NULL()
      Inst%SRCE_FUNC       => NULL()
      Inst%LND_FRC_DRY     => NULL()
      Inst%MSS_FRC_CACO3   => NULL()
      Inst%MSS_FRC_SND     => NULL()
      Inst%SFC_TYP         => NULL()
      Inst%VAI_DST         => NULL()

      ! Return w/ success
      RC = HCO_SUCCESS

      END SUBROUTINE InstCreate
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: InstRemove
!
! !DESCRIPTION: Subroutine InstRemove creates a new instance.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE InstRemove ( Instance )
!
! !INPUT PARAMETERS:
!
      INTEGER                         :: Instance
!
! !REVISION HISTORY:
!  18 Feb 2016 - C. Keller   - Initial version
!  26 Oct 2016 - R. Yantosca - Don't nullify local ptrs in declaration stmts
!EOP
!------------------------------------------------------------------------------
!BOC
      INTEGER                     :: RC
      TYPE(MyInst), POINTER       :: PrevInst
      TYPE(MyInst), POINTER       :: Inst

      !=================================================================
      ! InstRemove begins here!
      !=================================================================

      ! Get instance. Also archive previous instance.
      PrevInst => NULL()
      Inst     => NULL()
      CALL InstGet ( Instance, Inst, RC, PrevInst=PrevInst )

      ! Instance-specific deallocation
      IF ( ASSOCIATED(Inst) ) THEN

      !-----------------------------------------------------------------
      ! Nullify pointers
      !-----------------------------------------------------------------
      IF(ASSOCIATED( Inst%ERD_FCT_GEO  )) DEALLOCATE(Inst%ERD_FCT_GEO  )
      IF(ASSOCIATED( Inst%SRCE_FUNC    )) DEALLOCATE(Inst%SRCE_FUNC    )
      IF(ASSOCIATED( Inst%LND_FRC_DRY  )) DEALLOCATE(Inst%LND_FRC_DRY  )
      IF(ASSOCIATED( Inst%MSS_FRC_CACO3)) DEALLOCATE(Inst%MSS_FRC_CACO3)
      IF(ASSOCIATED( Inst%MSS_FRC_SND  )) DEALLOCATE(Inst%MSS_FRC_SND  )
      IF(ASSOCIATED( Inst%SFC_TYP      )) DEALLOCATE(Inst%SFC_TYP      )
      IF(ASSOCIATED( Inst%VAI_DST      )) DEALLOCATE(Inst%VAI_DST      )

      !-----------------------------------------------------------------
      ! Cleanup module arrays
      !-----------------------------------------------------------------

!      IF ( ALLOCATED( Inst%FLX_LW_DWN_SFC  ) )
!     &    DEALLOCATE( Inst%FLX_LW_DWN_SFC  )
!      IF ( ALLOCATED( Inst%FLX_SW_ABS_SFC  ) )
!     &    DEALLOCATE( Inst%FLX_SW_ABS_SFC  )
!      IF ( ALLOCATED( Inst%TPT_GND         ) )
!     &    DEALLOCATE( Inst%TPT_GND         )
!      IF ( ALLOCATED( Inst%TPT_SOI         ) )
!     &    DEALLOCATE( Inst%TPT_SOI         )
!      IF ( ALLOCATED( Inst%VWC_SFC         ) )
!     &    DEALLOCATE( Inst%VWC_SFC         )
!      IF ( ALLOCATED( Inst%SRC_STR         ) )
!     &    DEALLOCATE( Inst%SRC_STR         )
      IF ( ALLOCATED( Inst%PLN_TYP         ) )
     &    DEALLOCATE( Inst%PLN_TYP         )
      IF ( ALLOCATED( Inst%PLN_FRC         ) )
     &    DEALLOCATE( Inst%PLN_FRC         )
      IF ( ALLOCATED( Inst%TAI             ) )
     &    DEALLOCATE( Inst%TAI             )
      IF ( ALLOCATED( Inst%DMT_VWR         ) )
     &    DEALLOCATE( Inst%DMT_VWR         )
!      IF ( ALLOCATED( Inst%DNS_AER         ) )
!     &    DEALLOCATE( Inst%DNS_AER         )
      IF ( ALLOCATED( Inst%OVR_SRC_SNK_FRC ) )
     &    DEALLOCATE( Inst%OVR_SRC_SNK_FRC )
      IF ( ALLOCATED( Inst%OVR_SRC_SNK_MSS ) )
     &    DEALLOCATE( Inst%OVR_SRC_SNK_MSS )
!      IF ( ALLOCATED( Inst%OROGRAPHY       ) )
!     &    DEALLOCATE( Inst%OROGRAPHY       )
      IF ( ALLOCATED( Inst%DMT_MIN         ) )
     &    DEALLOCATE( Inst%DMT_MIN         )
      IF ( ALLOCATED( Inst%DMT_MAX         ) )
     &    DEALLOCATE( Inst%DMT_MAX         )
      IF ( ALLOCATED( Inst%DMT_VMA_SRC     ) )
     &    DEALLOCATE( Inst%DMT_VMA_SRC     )
      IF ( ALLOCATED( Inst%GSD_ANL_SRC     ) )
     &    DEALLOCATE( Inst%GSD_ANL_SRC     )
      IF ( ALLOCATED( Inst%MSS_FRC_SRC     ) )
     &    DEALLOCATE( Inst%MSS_FRC_SRC     )
      IF ( ALLOCATED( Inst%HcoIDs          ) )
     &    DEALLOCATE( Inst%HcoIDs          )
      IF ( ALLOCATED( Inst%HcoIDsALK       ) )
     &    DEALLOCATE( Inst%HcoIDsALK       )

       ! Pop off instance from list
       IF ( ASSOCIATED(PrevInst) ) THEN
          PrevInst%NextInst => Inst%NextInst
       ELSE
          AllInst => Inst%NextInst
       ENDIF
       DEALLOCATE(Inst)
       Inst => NULL()
      ENDIF

      END SUBROUTINE InstRemove
!EOC
      END MODULE HCOX_DUSTDEAD_MOD
!EOM
