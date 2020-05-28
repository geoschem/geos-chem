!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: hcox_finn_mod.F90
!
! !DESCRIPTION: Module HCOX\_FINN\_MOD contains routines and variables to
! calculate FINN biomass burning emissions in HEMCO.
!
! !INTERFACE:
!
MODULE HcoX_FINN_Mod
!
! !USES:
!
  USE HCO_Error_Mod
  USE HCO_Diagn_Mod
  USE HCOX_TOOLS_MOD
  USE HCO_State_Mod,  ONLY : HCO_State
  USE HCOX_State_Mod, ONLY : Ext_State

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC :: HCOX_FINN_Init
  PUBLIC :: HCOX_FINN_Run
  PUBLIC :: HCOX_FINN_Final
!
! !REMARKS:
!  Emissions of biomass burning species are read at monthly or daily
!  resolution. Note: no emission factors are used here - emissions of
!  individual species are given in input files. Emissions on the FINN 0.5x0.5
!  degree grid are regridded to the current model grid.
!                                                                             .
!  FINN biomass burning emissions are computed for the following gas-phase
!  and aerosol-phase species:
!                                                                             .
!     (1 ) NOx  [ kg/m2/s]     (13) BC   [kgC/m2/s]
!     (2 ) CO   [ kg/m2/s]     (14) OC   [kgC/m2/s]
!     (3 ) ALK4 [kgC/m2/s]     (15) MGLY [ kg/m2/s]
!     (4 ) ACET [kgC/m2/s]     (16) BENZ [kgC/m2/s]
!     (5 ) MEK  [kgC/m2/s]     (17) TOLU [kgC/m2/s]
!     (6 ) ALD2 [kgC/m2/s]     (18) C2H4 [kgC/m2/s]
!     (7 ) PRPE [kgC/m2/s]     (19) C2H2 [kgC/m2/s]
!     (8 ) C3H8 [kgC/m2/s]     (20) GLYC [ kg/m2/s]
!     (9 ) CH2O [ kg/m2/s]     (21) HAC  [ kg/m2/s]
!     (10) C2H6 [kgC/m2/s]     (22) CO2  [ kg/m2/s]
!     (11) SO2  [ kg/m2/s]     (23) CH4  [ kg/m2/s]
!     (12) NH3  [ kg/m2/s]     (24)
!
!
! All species to be used must be listed in the settings section of the HEMCO
! configuration file. For every listed species, individual scale factors as
! well as masks can be defined. For example, to scale FINN CO emissions by a
! factor of 1.05 and restrict them to North America, as well as to scale NO
! emissions by a factor of 1.5:
!
!114     FINN              : on    NO/CO/ALK4/ACET/MEK/ALD2/PRPE/C3H8/CH2O/C2H6/SO2/NH3/BCPI/BCPO/OCPI/OCPO/GLYC/HAC
!    --> FINN_daily        :       false
!    --> hydrophilic BC    :       0.2
!    --> hydrophilic OC    :       0.5
!    --> Mask_CO           :       NAMASK
!    --> Scaling_CO        :       1.05
!    --> Scaling_NO        :       1.5
!
! Field NAMASK must be defined in section mask of the HEMCO configuration file.
!
!  References:
!  ============================================================================
!  (1 ) Original FINN database from Christine Wiedinmyer
!        http://bai.acd.ucar.edu/Data/fire/
!  (2 ) Wiedinmyer, C., Akagi, S.K., Yokelson, R.J., Emmons, L.K.,
!       Al-Saadi, J.A., Orlando, J.J., and Soja, A.J.: The Fire
!       INventory from NCAR (FINN): a high resolution global model to
!       estimate the emissions from open burning, Geoscientific Model
!       Development, 4, 625-641, doi:10.5194/gmd-4-625-2011, 2011.
!
! !REVISION HISTORY:
!  02 Jan 2013 - J. Mao & J.A. Fisher - Initial version, based on GFED3
!  01 Oct 2013 - J.A. Fisher - Update to only use one input file
!  05 May 2014 - J.A. Fisher - Replace NOx emissions with NO emissions as part
!                              of removal of NOx-Ox partitioning
!  18 Jun 2014 - C. Keller   - Now a HEMCO extension.
!  03 Jul 2014 - C. Keller   - Added 13 new FINN species
!  11 Aug 2014 - R. Yantosca - Now get emission factors and NMOC ratios from
!                              hard-coded statements in hcox_finn_include.H
!  11 Aug 2014 - R. Yantosca - Now use F90 free-form indentation
!  11 Aug 2014 - R. Yantosca - Cosmetic changes to ProTeX subroutine headers
!  11 Jun 2015 - C. Keller   - Update to include individual scale factors and
!                              masks.
!  14 Oct 2016 - C. Keller   - Now use HCO_EvalFld instead of HCO_GetPtr.
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
  !=================================================================
  ! MODULE PARAMETERS
  !
  ! nSpcMax : Maximum number of emitted species
  ! N_EMFAC : Number of emission factors per species
  ! N_SPEC  : Number of FINN species
  ! MW_CO2  : Molecular weight of CO2  (g/mol)
  ! MW_NMOC : Molecular weight of NMOC (g/mol). Assumed MW for NMOC
  !           is 68 g/mol.
  !=================================================================
  INTEGER,           PARAMETER   :: nSpcMax = 100
  INTEGER,           PARAMETER   :: N_EMFAC = 6
  INTEGER,           PARAMETER   :: N_SPEC  = 58
  REAL(dp),          PARAMETER   :: MW_CO2  = 44.01_dp
  REAL(dp),          PARAMETER   :: MW_NMOC = 68.00_dp
!
! !PRIVATE TYPES:
!
  TYPE :: MyInst
   !=================================================================
   ! HEMCO VARIABLES
   !
   ! ExtNr   : Extension number
   ! UseDay  : True if daily data is used
   !=================================================================
   INTEGER                        :: Instance
   INTEGER                        :: ExtNr
   LOGICAL                        :: UseDay

   !=================================================================
   ! SPECIES VARIABLES
   !
   ! nSpc           : Number of used species (specified in config. file)
   ! SpcNames       : Names of all used species
   ! HcoIDs         : HEMCO species IDs of all used species
   ! FinnIDs        : Index of used species within FINN
   ! FINN_SPEC_NAME : Names of all FINN species
   !=================================================================
   INTEGER                        :: nSpc
   CHARACTER(LEN=31), POINTER     :: SpcNames(:)
   CHARACTER(LEN=61), POINTER     :: SpcScalFldNme(:)
   INTEGER,           POINTER     :: HcoIDs(:)
   INTEGER,           POINTER     :: FinnIDs(:)
   CHARACTER(LEN=6),  POINTER     :: FINN_SPEC_NAME(:)

   !=================================================================
   ! SCALE FACTORS
   !
   ! FINN_EMFAC: emission scale factors for each species and
   !             emission factor type. The filename of the emissions
   !             emissions factor table is specified in the HEMCO
   !             configuration file. The scale factors are converted
   !             to kg species/kg CO2 when reading them from disk.
   ! OCPIfrac  : Fraction of OC that converts into hydrophilic OC.
   !             Can be set in HEMCO configuration file (default=0.5)
   ! BCPIfrac  : Fraction of BC that converts into hydrophilic BC.
   !             Can be set in HEMCO configuration file (default=0.2)
   ! SpcScal  : Additional scaling factors assigned to species through
   !            the HEMCO configuration file (e.g. Scaling_CO).
   !=================================================================
   REAL(dp),          POINTER     :: FINN_EMFAC(:,:)
   REAL(sp),          POINTER     :: SpcScal(:)
   REAL(sp)                       :: OCPIfrac
   REAL(sp)                       :: BCPIfrac

   !=================================================================
   ! DATA ARRAY POINTERS
   !
   ! These are the pointers to the 6 vegetation type data arrays
   ! specified in the configuration file
   !=================================================================
   REAL(hp),          POINTER     :: VEGTYP1(:,:) => NULL()
   REAL(hp),          POINTER     :: VEGTYP2(:,:) => NULL()
   REAL(hp),          POINTER     :: VEGTYP3(:,:) => NULL()
   REAL(hp),          POINTER     :: VEGTYP4(:,:) => NULL()
   REAL(hp),          POINTER     :: VEGTYP5(:,:) => NULL()
   REAL(hp),          POINTER     :: VEGTYP9(:,:) => NULL()

   TYPE(MyInst), POINTER           :: NextInst => NULL()
  END TYPE MyInst

  ! Pointer to instances
  TYPE(MyInst), POINTER            :: AllInst => NULL()

CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCOX_FINN_Run
!
! !DESCRIPTION: Subroutine HCOX\_FINN\_Run computes the FINN biomass
!  burning emissions for the current date.
!\\
!\\
! !INTERFACE:
  !
  SUBROUTINE HCOX_FINN_Run( ExtState, HcoState, RC )
!
! !USES:
!
    USE HCO_EmisList_mod,  ONLY : HCO_GetPtr
    USE HCO_Calc_Mod,      ONLY : HCO_EvalFld
    USE HCO_FluxArr_mod,   ONLY : HCO_EmisAdd
    USE HCO_State_mod,     ONLY : HCO_GetHcoID
    USE HCO_Clock_mod,     ONLY : HcoClock_Get
    USE HCO_Clock_mod,     ONLY : HcoClock_First
    USE HCO_Clock_mod,     ONLY : HcoClock_NewMonth, HcoClock_NewDay
!
! !INPUT PARAMETERS:
!
    TYPE(Ext_State), POINTER        :: ExtState   ! Module options
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(HCO_State), POINTER        :: HcoState   ! Output obj
    INTEGER,         INTENT(INOUT)  :: RC         ! Success or failure?
!
! !REVISION HISTORY:
!  02 Jan 2012 - J. Mao & J. Fisher - Initial version, based on GFED3
!  18 Jun 2014 - C. Keller          - Now a HEMCO extension.
!  26 Oct 2016 - R. Yantosca - Don't nullify local ptrs in declaration stmts
!  10 Mar 2017 - M. Sulprizio- Add SpcArr3D for emitting 65% of biomass
!                              burning emissions into the PBL and 35% into the
!                              free troposphere, following code from E.Fischer
!  24 Apr 2017 - M. Sulprizio- Comment out vertical distribution of biomass
!                              burning emissions for now.
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER             :: N, M, NF
    INTEGER             :: FinnID, HcoID
!    LOGICAL, SAVE       :: FIRST = .TRUE.
    LOGICAL             :: DoRepeat
    INTEGER             :: Cnt
    CHARACTER(LEN=31)   :: PREFIX, FLDNME
    INTEGER             :: NDAYS, cYYYY, cMM, cDD
    REAL(dp)            :: TOTAL
    CHARACTER(LEN=255)  :: MSG

    ! Arrays
    REAL(hp), TARGET    :: SpcArr(HcoState%NX,HcoState%NY)
    REAL(hp), TARGET    :: TypArr(HcoState%NX,HcoState%NY)

!==============================================================================
! This code is required for the vertical distribution of biomass burning emiss.
! We will keep it here for a future implementation. (mps, 4/24/17)
!    INTEGER             :: I, J, L, N, M
!    INTEGER             :: PBL_MAX
!    REAL(hp)            :: PBL_FRAC, F_OF_PBL, F_OF_FT
!    REAL(hp)            :: DELTPRES, TOTPRESFT
!    REAL(hp), TARGET    :: SpcArr3D(HcoState%NX,HcoState%NY,HcoState%NZ)
!==============================================================================

    ! Pointers
    REAL(hp), POINTER   :: THISTYP(:,:)

    ! Local instance
    TYPE(MyInst), POINTER :: Inst

    !=======================================================================
    ! HCOX_FINN_Run begins here!
    !=======================================================================

    ! Return if extension disabled
    IF ( ExtState%FINN <= 0 ) RETURN

    ! Enter
    CALL HCO_ENTER( HcoState%Config%Err, 'HCOX_FINN_RUN (hcox_finn_mod.F90)', RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Init
    THISTYP => NULL()

    ! Get instance
    Inst => NULL()
    CALL InstGet ( ExtState%FINN, Inst, RC )
    IF ( RC /= HCO_SUCCESS ) THEN
       WRITE(MSG,*) 'Cannot find FINN instance Nr. ', ExtState%FINN
       CALL HCO_ERROR(HcoState%Config%Err,MSG,RC)
       RETURN
    ENDIF

!==============================================================================
! This code is required for the vertical distribution of biomass burning emiss.
! We will keep it here for a future implementation. (mps, 4/24/17)
!    ! Add only 65% biomass burning source to boundary layer, the
!    ! rest is emitted into the free troposphere (mps from evf+tjb, 3/10/17)
!    PBL_FRAC = 0.65_hp
!==============================================================================

    !-----------------------------------------------------------------------
    ! Get pointers to data arrays
    !-----------------------------------------------------------------------
    !IF ( HcoClock_First(HcoState%Clock,.TRUE.) ) THEN
       IF ( Inst%UseDay ) THEN
          PREFIX = 'FINN_DAILY_'
       ELSE
          PREFIX = 'FINN_'
       ENDIF

       FLDNME = TRIM(PREFIX) // 'VEGTYP1'
       CALL HCO_EvalFld( HcoState, TRIM(FLDNME), Inst%VEGTYP1, RC )
       IF ( RC /= HCO_SUCCESS ) RETURN

       FLDNME = TRIM(PREFIX) // 'VEGTYP2'
       CALL HCO_EvalFld( HcoState, TRIM(FLDNME), Inst%VEGTYP2, RC )
       IF ( RC /= HCO_SUCCESS ) RETURN

       FLDNME = TRIM(PREFIX) // 'VEGTYP3'
       CALL HCO_EvalFld( HcoState, TRIM(FLDNME), Inst%VEGTYP3, RC )
       IF ( RC /= HCO_SUCCESS ) RETURN

       FLDNME = TRIM(PREFIX) // 'VEGTYP4'
       CALL HCO_EvalFld( HcoState, TRIM(FLDNME), Inst%VEGTYP4, RC )
       IF ( RC /= HCO_SUCCESS ) RETURN

       FLDNME = TRIM(PREFIX) // 'VEGTYP5'
       CALL HCO_EvalFld( HcoState, TRIM(FLDNME), Inst%VEGTYP5, RC )
       IF ( RC /= HCO_SUCCESS ) RETURN

       FLDNME = TRIM(PREFIX) // 'VEGTYP9'
       CALL HCO_EvalFld( HcoState, TRIM(FLDNME), Inst%VEGTYP9, RC )
       IF ( RC /= HCO_SUCCESS ) RETURN

!       FIRST = .FALSE.
    !ENDIF

    ! For logfile
    IF ( HcoState%amIRoot ) THEN
       IF ( Inst%UseDay ) THEN
          IF ( HcoClock_NewDay( HcoState%Clock, .TRUE. ) ) THEN
             CALL HcoClock_Get( HcoState%Clock, &
                                cYYYY=cYYYY, cMM=cMM, cDD=cDD, RC=RC )
             IF ( RC/=HCO_SUCCESS ) RETURN
             WRITE(MSG, 100) cYYYY, cMM, cDD
             CALL HCO_MSG(HcoState%Config%Err,MSG)
100          FORMAT( 'FINN daily emissions for year, month, day: ', &
                      i4, '/', i2.2, '/', i2.2 )
          ENDIF
       ELSE
          IF ( HcoClock_NewMonth( HcoState%Clock, .TRUE. ) ) THEN
             CALL HcoClock_Get( HcoState%Clock, &
                                cYYYY=cYYYY, cMM=cMM, LMD=NDAYS, RC=RC)
             IF ( RC/=HCO_SUCCESS ) RETURN
             WRITE(MSG, 110) cYYYY, cMM
             CALL HCO_MSG(HcoState%Config%Err,MSG)
110          FORMAT( 'FINN monthly emissions for year, month: ', &
                      i4, '/', i2.2 )
          ENDIF
       ENDIF
    ENDIF

    !-----------------------------------------------------------------------
    ! Calculate emissions for all selected species
    !-----------------------------------------------------------------------

    ! Loop over all emitted species
    DO N = 1, Inst%nSpc

       ! ID is the FINN species index of this species
       FinnID = Inst%FinnIDs(N)
       IF ( FinnID <= 0 ) CYCLE

       ! HcoID is the HEMCO species index of this species
       HcoID = Inst%HcoIDs(N)
       IF ( HcoID < 0 ) CYCLE

       ! Species with no emission factor have FINN_EMFAC=0
       IF ( MAXVAL(Inst%FINN_EMFAC(FinnID,:)) <= 0.0_hp ) CYCLE

       ! SpcArr are the total biomass burning emissions for this
       ! species. TypArr are the emissions from a given vegetation type.
       SpcArr   = 0.0_hp
!==============================================================================
! This code is required for the vertical distribution of biomass burning emiss.
! We will keep it here for a future implementation. (mps, 4/24/17)
!       SpcArr3D = 0.0_hp
!==============================================================================

       ! Calculate emissions for all source types
       DO NF = 1, N_EMFAC

          ! Select emission factor array
          IF ( NF == 1 ) THEN
             THISTYP => Inst%VEGTYP1
          ELSEIF ( NF == 2 ) THEN
             THISTYP => Inst%VEGTYP2
          ELSEIF ( NF == 3 ) THEN
             THISTYP => Inst%VEGTYP3
          ELSEIF ( NF == 4 ) THEN
             THISTYP => Inst%VEGTYP4
          ELSEIF ( NF == 5 ) THEN
             THISTYP => Inst%VEGTYP5
          ELSEIF ( NF == 6 ) THEN
             THISTYP => Inst%VEGTYP9
          ELSE
             CALL HCO_ERROR ( HcoState%Config%Err, 'Undefined emission factor', RC )
             RETURN
          ENDIF

          ! Multiply CO2 emissions by appropriate ratio for each land
          ! type and sum to get total emissions for the species on the
          ! native grid - emissions are in [kg CO2/m2/s[. FINN_EMFAC is
          ! in [kg X]/[kg CO2].
          TypArr(:,:) = THISTYP(:,:) * Inst%FINN_EMFAC(FinnID,NF)

          ! TODO: Add to diagnostics here

          ! Add to species array
          SpcArr = SpcArr + TypArr
       ENDDO !NF

       ! Apply species specific scale factors
       SpcArr = SpcArr * Inst%SpcScal(N)

       ! Check for masking
       CALL HCOX_SCALE( HcoState, SpcArr, TRIM(Inst%SpcScalFldNme(N)), RC )
       IF ( RC /= HCO_SUCCESS ) RETURN

       SELECT CASE ( Inst%SpcNames(N) )
          CASE ( 'OCPI' )
             SpcArr = SpcArr * Inst%OCPIfrac
          CASE ( 'OCPO' )
             SpcArr = SpcArr * (1.0_sp - Inst%OCPIfrac)
          CASE ( 'BCPI' )
             SpcArr = SpcArr * Inst%BCPIfrac
          CASE ( 'BCPO' )
             SpcArr = SpcArr * (1.0_sp - Inst%BCPIfrac)
       END SELECT

!==============================================================================
! This code is required for the vertical distribution of biomass burning emiss.
! We will keep it here for a future implementation. (mps, 4/24/17)
!
!       !--------------------------------------------------------------------
!       ! For grid boxes with emissions, distribute 65% to PBL and 35% to FT
!       !--------------------------------------------------------------------
!       DO J = 1, HcoState%Ny
!       DO I = 1, HcoState%Nx
!
!          IF ( SpcArr(I,J) > 0e+0_hp ) THEN
!
!             ! Initialize
!             PBL_MAX  = 1
!             F_OF_PBL = 0e+0_hp
!             F_OF_FT  = 0e+0_hp
!             DELTPRES = 0e+0_hp
!
!             ! Determine PBL height
!             DO L = HcoState%NZ, 1, -1
!                IF ( ExtState%FRAC_OF_PBL%Arr%Val(I,J,L) > 0.0_hp ) THEN
!                   PBL_MAX = L
!                   EXIT
!                ENDIF
!             ENDDO
!
!             ! Loop over the boundary layer
!             DO L = 1, PBL_MAX
!
!                ! Fraction of PBL that box (I,J,L) makes up [unitless]
!                F_OF_PBL = ExtState%FRAC_OF_PBL%Arr%Val(I,J,L)
!
!                ! Add only 65% biomass burning source to PBL
!                ! Distribute emissions thru the entire boundary layer
!                ! (mps from evf+tjb, 3/10/17)
!                SpcArr3D(I,J,L) = SpcArr(I,J) * PBL_FRAC * F_OF_PBL
!
!             ENDDO
!
!             ! Total thickness of the free troposphere [hPa]
!             ! (considered here to be 10 levels above the PBL)
!             TOTPRESFT = HcoState%Grid%PEDGE%Val(I,J,PBL_MAX+1) - &
!                         HcoState%Grid%PEDGE%Val(I,J,PBL_MAX+11)
!
!             ! Loop over the free troposphere
!             DO L = PBL_MAX+1, PBL_MAX+10
!
!                ! Thickness of level L [hPa]
!                DELTPRES= HcoState%Grid%PEDGE%Val(I,J,L) - &
!                          HcoState%Grid%PEDGE%Val(I,J,L+1)
!
!                ! Fraction of FT that box (I,J,L) makes up [unitless]
!                F_OF_FT = DELTPRES / TOTPRESFT
!
!                ! Add 35% of biomass burning source to free troposphere
!                ! Distribute emissions thru 10 model levels above the BL
!                ! (mps from evf+tjb, 3/10/17)
!                SpcArr3D(I,J,L) = SpcArr(I,J) * (1.0-PBL_FRAC) * F_OF_FT
!
!             ENDDO
!
!          ENDIF
!
!       ENDDO
!       ENDDO
!
!       ! Add flux to HEMCO emission array
!       ! Now 3D flux (mps, 3/10/17)
!       CALL HCO_EmisAdd( HcoState, SpcArr3D, HcoID, &
!                         RC,       ExtNr=ExtNr, Cat=-1,   Hier=-1 )
!==============================================================================

       ! Add flux to HEMCO emission array
       CALL HCO_EmisAdd( HcoState, SpcArr, HcoID, &
                         RC,       ExtNr=Inst%ExtNr, Cat=-1, Hier=-1 )
       IF ( RC /= HCO_SUCCESS ) THEN
          MSG = 'HCO_EmisAdd error: ' // TRIM(HcoState%Spc(HcoID)%SpcName)
          CALL HCO_ERROR(HcoState%Config%Err,MSG, RC )
          RETURN
       ENDIF

       ! Write out total (daily or monthly) emissions to log-file
       IF ( HcoState%amIRoot ) THEN
          IF ( Inst%UseDay ) THEN
             IF ( HcoClock_NewDay( HcoState%Clock, .TRUE. ) ) THEN
                TOTAL = SUM(SpcArr(:,:)*HcoState%Grid%AREA_M2%Val(:,:))
                TOTAL = TOTAL * 86400.0_hp * 1e-9_hp
                WRITE(MSG, 120) HcoState%Spc(HcoID)%SpcName, TOTAL
                CALL HCO_MSG(HcoState%Config%Err,MSG)
120             FORMAT( 'SUM biomass ', a4,1x,': ', f11.4,1x,'[Tg]' )
             ENDIF
          ELSE
             IF ( HcoClock_NewMonth( HcoState%Clock, .TRUE. ) ) THEN
                TOTAL = SUM(SpcArr(:,:)*HcoState%Grid%AREA_M2%Val(:,:))
                TOTAL = TOTAL * NDAYS * 86400.0_hp * 1e-9_hp
                WRITE(MSG, 130) HcoState%Spc(HcoID)%SpcName, TOTAL
                CALL HCO_MSG(HcoState%Config%Err,MSG)
130             FORMAT( 'SUM biomass ', a4,1x,': ', f11.4,1x,'[Tg]' )
             ENDIF
          ENDIF
       ENDIF

    ENDDO !N

    ! Nullify pointers
    THISTYP   => NULL()
    Inst      => NULL()

    ! Leave w/ success
    CALL HCO_LEAVE( HcoState%Config%Err,RC )

  END SUBROUTINE HCOX_FINN_Run
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCOX_FINN_Init
!
! !DESCRIPTION: Subroutine HCOX\_FINN\_INIT initializes all module
! arrays and variables.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCOX_FINN_Init( HcoState, ExtName, ExtState, RC )
!
! !USES:
!
    USE HCO_State_Mod,   ONLY : HCO_GetHcoID
    USE HCO_State_Mod,   ONLY : HCO_GetExtHcoID
    USE HCO_ExtList_Mod, ONLY : GetExtNr, GetExtOpt
    USE HCO_ExtList_Mod, ONLY : GetExtSpcVal
!
! !INPUT PARAMETERS:
!
    TYPE(HCO_State),  POINTER        :: HcoState    ! HEMCO state object
    CHARACTER(LEN=*), INTENT(IN   )  :: ExtName     ! Extension name
    TYPE(Ext_State),  POINTER        :: ExtState    ! Extensions object
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(INOUT)  :: RC          ! Return status
!
! !REVISION HISTORY:
!  02 Jan 2013 - J. Mao & J. Fisher - Initial version, based on GFED3
!  05 May 2014 - J.A. Fisher - Replace NOx emissions with NO emissions as part
!                              of removal of NOx-Ox partitioning
!  18 Jun 2014 - C. Keller   - Now a HEMCO extension.
!  11 Aug 2014 - R. Yantosca - Now get FINN emission factors and species names
!                              from include file hcox_finn_include.H.
!  11 Nov 2014 - C. Keller   - Now get hydrophilic fractions through config file
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES
!
    ! Scalars
    INTEGER               :: ExtNr
    INTEGER               :: N_SPEC_EMFAC  ! # of CO2 file emission species
    INTEGER               :: N_NMOC        ! # of VOC file NMOC ratios
    INTEGER               :: IU_FILE, L, N_LUMPED, tmpNr
    INTEGER               :: AS, IOS, M, N, NDUM
    INTEGER               :: N_SPECSTRS, N_NMOCSTRS
    INTEGER               :: NCHAR
    LOGICAL               :: IS_NMOC, Matched, Missing, FOUND
    CHARACTER(LEN=1023)   :: ADUM
    CHARACTER(LEN=255)    :: SDUM(255)
    CHARACTER(LEN=255)    :: IN_SPEC_NAME(255)
    CHARACTER(LEN=255)    :: IN_NMOC_NAME(255)
    CHARACTER(LEN=255)    :: TMPNAME
    CHARACTER(LEN=  6)    :: SPCNAME
    REAL*8                :: C_MOLEC
    REAL(dp)              :: AdjFact
    REAL(sp)              :: ValSp
    CHARACTER(LEN=255)    :: MSG, EF_CO2_FILE, VOC_SPEC_FILE

    ! Temporary variables. These values will be passed to module
    ! array nSpc, SpcNames, etc.
    INTEGER                        :: tnSpc
    CHARACTER(LEN=31), ALLOCATABLE :: tSpcNames(:)
    CHARACTER(LEN=61), ALLOCATABLE :: tSpcScalFldNme(:)
    REAL(sp),          ALLOCATABLE :: tSpcScal(:)
    INTEGER,           ALLOCATABLE :: tHcoIDs(:)

    ! Arrays
    REAL(dp), ALLOCATABLE :: EMFAC_IN(:,:)
    REAL(dp), ALLOCATABLE :: NMOC_RATIO_IN(:,:)
    REAL*8                :: NMOC_EMFAC(N_EMFAC), NMOC_RATIO(N_EMFAC)

    ! Local instance
    TYPE(MyInst), POINTER :: Inst

    !=======================================================================
    ! HCOX_FINN_INIT begins here!
    !=======================================================================

    ! Extension Nr.
    ExtNr = GetExtNr( HcoState%Config%ExtList, TRIM(ExtName) )
    IF ( ExtNr <= 0 ) RETURN

    ! Enter
    CALL HCO_ENTER( HcoState%Config%Err, 'HCOX_FINN_INIT (hcox_finn_mod.F90)', RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Create local instance for this simulation
    Inst => NULL()
    CALL InstCreate ( ExtNr, ExtState%FINN, Inst, RC )
    IF ( RC /= HCO_SUCCESS ) THEN
       CALL HCO_ERROR ( HcoState%Config%Err, 'Cannot create FINN instance', RC )
       RETURN
    ENDIF

    ! Check if this is GFED4
    !-----------------------------------------------------------------------
    ! Get settings
    ! The CO scale factor (to account for oxidation from VOCs) as well as
    ! the speciation of carbon aerosols into hydrophilic and hydrophobic
    ! fractions can be specified in the configuration file, e.g.:
    ! 100     GFED3           : on    NO/CO/OCPI/OCPO/BCPI/BCPO
    !     --> hydrophilic BC  :       0.2
    !     --> hydrophilic OC  :       0.5
    !
    ! Setting these values is optional and default values are applied if
    ! they are not specified. The values only take effect if the
    ! corresponding species (CO, BCPI/BCPO, OCPI/OCPO) are listed as species
    ! to be used.
    !-----------------------------------------------------------------------

    ! Try to read hydrophilic fractions of BC. Defaults to 0.2.
    CALL GetExtOpt( HcoState%Config, ExtNr, 'hydrophilic BC', &
                     OptValSp=ValSp, FOUND=FOUND, RC=RC )
    IF ( RC /= HCO_SUCCESS ) RETURN
    IF ( .NOT. FOUND ) THEN
       Inst%BCPIfrac = 0.2
    ELSE
       Inst%BCPIfrac = ValSp
    ENDIF

    ! Try to read hydrophilic fractions of OC. Defaults to 0.5.
    CALL GetExtOpt( HcoState%Config, ExtNr, 'hydrophilic OC', &
                     OptValSp=ValSp, FOUND=FOUND, RC=RC )
    IF ( RC /= HCO_SUCCESS ) RETURN
    IF ( .NOT. FOUND ) THEN
       Inst%OCPIfrac = 0.5
    ELSE
       Inst%OCPIfrac = ValSp
    ENDIF

    ! Error check: OCPIfrac and BCPI frac must be between 0 and 1
    IF ( Inst%OCPIfrac < 0.0_sp .OR. Inst%OCPIfrac > 1.0_sp .OR. &
         Inst%BCPIfrac < 0.0_sp .OR. Inst%BCPIfrac > 1.0_sp     ) THEN
       WRITE(MSG,*) 'hydrophilic fractions must be between 0-1: ', &
          Inst%OCPIfrac, Inst%BCPIfrac
       CALL HCO_ERROR(HcoState%Config%Err,MSG, RC )
       RETURN
    ENDIF

    ! Use daily data?
    CALL GetExtOpt( HcoState%Config, ExtNr, 'FINN_daily', &
                     OptValBool=Inst%UseDay, FOUND=FOUND, RC=RC )
    IF ( RC /= HCO_SUCCESS ) RETURN
    IF ( .NOT. FOUND ) THEN
       Inst%UseDay = .FALSE.
    ENDIF

    !-----------------------------------------------------------------------
    ! Allocate arrays
    !-----------------------------------------------------------------------

    ! FINN species names
    ALLOCATE ( Inst%FINN_SPEC_NAME ( N_SPEC ), STAT=AS )
    IF ( AS/=0 ) THEN
       CALL HCO_ERROR( HcoState%Config%Err, 'Cannot allocate FINN_SPEC_NAME', RC )
       RETURN
    ENDIF
    Inst%FINN_SPEC_NAME = ''

    ! Allocate scale factors table: FINN_EMFAC holds the species/CO2
    ! scale factors for all FINN species.
    ALLOCATE ( Inst%FINN_EMFAC ( N_SPEC, N_EMFAC ), STAT=AS )
    IF ( AS/=0 ) THEN
       CALL HCO_ERROR( HcoState%Config%Err, 'Cannot allocate FINN_EMFAC', RC )
       RETURN
    ENDIF
    Inst%FINN_EMFAC = 0.0_dp

    ! Allocate and initialize vectors holding species information for
    ! all species to be emitted
    ALLOCATE ( Inst%FinnIDs(nSpcMax), Inst%HcoIDs(nSpcMax), Inst%SpcNames(nSpcMax), &
               Inst%SpcScal(nSpcMax), Inst%SpcScalFldNme(nSpcMax), STAT=AS )

    IF ( AS/=0 ) THEN
       CALL HCO_ERROR( HcoState%Config%Err, 'Cannot allocate FinnIDs', RC )
       RETURN
    ENDIF
    Inst%nSpc             = 0
    Inst%FinnIDs(:)       = -1
    Inst%HcoIDs(:)        = -1
    Inst%SpcScal          = 1.0_sp
    Inst%SpcNames(:)      = ''
    Inst%SpcScalFldNme(:) = HCOX_NOSCALE

    ALLOCATE ( Inst%VEGTYP1(HcoState%NX,HcoState%NY), &
               Inst%VEGTYP2(HcoState%NX,HcoState%NY), &
               Inst%VEGTYP3(HcoState%NX,HcoState%NY), &
               Inst%VEGTYP4(HcoState%NX,HcoState%NY), &
               Inst%VEGTYP5(HcoState%NX,HcoState%NY), &
               Inst%VEGTYP9(HcoState%NX,HcoState%NY), STAT=AS )
    IF ( AS/=0 ) THEN
       CALL HCO_ERROR( HcoState%Config%Err, 'Cannot allocate VEGTYP', RC )
       RETURN
    ENDIF
    Inst%VEGTYP1 = 0.0_hp
    Inst%VEGTYP2 = 0.0_hp
    Inst%VEGTYP3 = 0.0_hp
    Inst%VEGTYP4 = 0.0_hp
    Inst%VEGTYP5 = 0.0_hp
    Inst%VEGTYP9 = 0.0_hp

    !-----------------------------------------------------------------------
    ! Define FINN species names
    !-----------------------------------------------------------------------

    ! Species listed in emission factor ratios (CO2/X) table (except NMOC,
    ! which is speciated as specified in the VOC speciation table).
    Inst%FINN_SPEC_NAME(1)  = 'CO2'
    Inst%FINN_SPEC_NAME(2)  = 'CO'
    Inst%FINN_SPEC_NAME(3)  = 'CH4'
    Inst%FINN_SPEC_NAME(4)  = 'NOx'
    Inst%FINN_SPEC_NAME(5)  = 'SO2'
    Inst%FINN_SPEC_NAME(6)  = 'OC'
    Inst%FINN_SPEC_NAME(7)  = 'BC'
    Inst%FINN_SPEC_NAME(8)  = 'NH3'
    Inst%FINN_SPEC_NAME(9)  = 'NO'    ! Currently not used
    Inst%FINN_SPEC_NAME(10) = 'NO2'   ! Currently not used

    ! Species listed in VOC speciation table
    Inst%FINN_SPEC_NAME(11) = 'ACET'
    Inst%FINN_SPEC_NAME(12) = 'ACTA'   ! Not currently emitted by BB in GC
    Inst%FINN_SPEC_NAME(13) = 'ALD2'
    Inst%FINN_SPEC_NAME(14) = 'ALK4'
    Inst%FINN_SPEC_NAME(15) = 'APINE'  ! Currently lumped into MTPA
    Inst%FINN_SPEC_NAME(16) = 'AROM'   ! Currently not used
    Inst%FINN_SPEC_NAME(17) = 'BENZ'
    Inst%FINN_SPEC_NAME(18) = 'BPINE'  ! Currently lumped into MTPA
    Inst%FINN_SPEC_NAME(19) = 'C2H2'
    Inst%FINN_SPEC_NAME(20) = 'C2H4'
    Inst%FINN_SPEC_NAME(21) = 'C2H6'
    Inst%FINN_SPEC_NAME(22) = 'C3H8'
    Inst%FINN_SPEC_NAME(23) = 'CARENE' ! Currently lumped into MTPA
    Inst%FINN_SPEC_NAME(24) = 'CH2Br2'
    Inst%FINN_SPEC_NAME(25) = 'CH2O'
    Inst%FINN_SPEC_NAME(26) = 'CH3Br'
    Inst%FINN_SPEC_NAME(27) = 'CH3CN'
    Inst%FINN_SPEC_NAME(28) = 'CH3I'
    Inst%FINN_SPEC_NAME(29) = 'DMS'
    Inst%FINN_SPEC_NAME(30) = 'EOH'    ! Not currently emitted in GC
    Inst%FINN_SPEC_NAME(31) = 'ETBENZ' ! Currently lumped with TOLU
    Inst%FINN_SPEC_NAME(32) = 'FUR'    ! Currently not used
    Inst%FINN_SPEC_NAME(33) = 'GLYC'
    Inst%FINN_SPEC_NAME(34) = 'GLYX'
    Inst%FINN_SPEC_NAME(35) = 'HAC'
    Inst%FINN_SPEC_NAME(36) = 'HCN'    ! Not currently emitted in GC
    Inst%FINN_SPEC_NAME(37) = 'HCOOH'  ! Not currently emitted by BB in GC
    Inst%FINN_SPEC_NAME(38) = 'HNO2'   ! Not currently emitted in GC
    Inst%FINN_SPEC_NAME(39) = 'ISOP'   ! Not currently emitted by BB in GC
    Inst%FINN_SPEC_NAME(40) = 'LIMO'
    Inst%FINN_SPEC_NAME(41) = 'MACR'   ! Not currently emitted in GC
    Inst%FINN_SPEC_NAME(42) = 'MEK'
    Inst%FINN_SPEC_NAME(43) = 'MGLY'
    Inst%FINN_SPEC_NAME(44) = 'MNO3'
    Inst%FINN_SPEC_NAME(45) = 'MOH'    ! Not currently emitted in GC
    Inst%FINN_SPEC_NAME(46) = 'MTPO'   ! Not currently emitted in GC
    Inst%FINN_SPEC_NAME(47) = 'MVK'    ! Not currently emitted in GC
    Inst%FINN_SPEC_NAME(48) = 'PRPE'
    Inst%FINN_SPEC_NAME(49) = 'R4N2'   ! Not currently emitted in GC
    Inst%FINN_SPEC_NAME(50) = 'RCHO'   ! Not currently emitted by BB in GC
    Inst%FINN_SPEC_NAME(51) = 'RCOOH'  ! Currently not used
    Inst%FINN_SPEC_NAME(52) = 'ROH'    ! Currently not used
    Inst%FINN_SPEC_NAME(53) = 'SESQ'   ! Currently not used
    Inst%FINN_SPEC_NAME(54) = 'STYR'   ! Currently lumped with TOLU
    Inst%FINN_SPEC_NAME(55) = 'TMB'    ! Currently lumped with XYLE
    Inst%FINN_SPEC_NAME(56) = 'TOLU'
    Inst%FINN_SPEC_NAME(57) = 'XYLE'
    Inst%FINN_SPEC_NAME(58) = 'H2'     ! Currently not used

    !=======================================================================
    ! We now get the following input information from hard-coded F90
    ! assignment statements in the include file "hcox_finn_include.H":
    !
    ! Quantities formerly defined in the "FINN_EFratios_CO2.csv" file:
    ! ----------------------------------------------------------------------
    ! (1 ) N_SPEC_EMFAC  : # of species in the FINN_EFratios_CO2.csv file
    ! (2 ) N_SPECSTRS    : Synonym for N_SPEC_EMFAC
    ! (3 ) IN_SPEC_NAME  : Name of emissions species
    ! (4 ) EMFAC_IN      : Emission ratios for each species
    !
    ! Quantities formerly defined in the "FINN_VOC_speciation.csv" file:
    ! ----------------------------------------------------------------------
    ! (5 ) N_NMOC_       : # of species in the FINN_EFratios_CO2.csv file
    ! (6 ) N_NMOCSTRS    : Synonym for N_NMOC
    ! (7 ) IN_NMOC_NAME  : Name of NMOC ratios
    ! (8 ) NMOC_RATIO_IN : NMOC ratios for each species
    !
    ! Furthermore, the F90 statements to allocate the arrays IN_SPEC_NAME
    ! and IN_NMOC_NAME are included in "hcox_finn_include.H".
    !
    ! NOTE: If new FINN emisison factors and NMOC ratios are issued in the
    ! future, you can regenerate the include file "hcox_finn_include.H"
    ! with the Perl script HEMCO/Extensions/Preprocess/finn.pl.
    !=======================================================================
#include "hcox_finn_include.H"

    !-----------------------------------------------------------------------
    ! Match specified species with FINN species. The species to be used are
    ! specified in the HEMCO configuration file.
    !-----------------------------------------------------------------------

    ! Write to log file
    IF ( HcoState%amIRoot ) THEN
       MSG = 'Use FINN extension'
       CALL HCO_MSG(HcoState%Config%Err,MSG, SEP1='-' )
       WRITE(MSG,*) '   - Use daily data          : ', Inst%UseDay
       CALL HCO_MSG(HcoState%Config%Err,MSG )
    ENDIF

    ! Get HEMCO species IDs of all species specified in configuration file
    CALL HCO_GetExtHcoID( HcoState, ExtNr, tHcoIDs, tSpcNames, tnSpc, RC)
    IF ( RC /= HCO_SUCCESS ) RETURN
    IF ( tnSpc == 0 ) THEN
       MSG = 'No FINN species specified'
       CALL HCO_ERROR(HcoState%Config%Err,MSG, RC )
       RETURN
    ENDIF

    ! Get species scale factors
    CALL GetExtSpcVal( HcoState%Config, ExtNr, tnSpc, &
                       tSpcNames, 'Scaling', 1.0_sp, tSpcScal, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Get species mask fields
    CALL GetExtSpcVal( HcoState%Config, ExtNr, tnSpc, &
                       tSpcNames, 'ScaleField', HCOX_NOSCALE, tSpcScalFldNme, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Error trap: in previous versions, CO, POA and NAP scale factor were given as
    ! 'CO scale factor', etc. Make sure those attributes do not exist any more!
    CALL GetExtOpt( HcoState%Config, ExtNr, 'CO scale factor', &
                     OptValSp=ValSp, FOUND=FOUND, RC=RC )
    IF ( .NOT. FOUND ) THEN
       CALL GetExtOpt( HcoState%Config, ExtNr, 'POA scale factor', &
                        OptValSp=ValSp, FOUND=FOUND, RC=RC )
    ENDIF
    IF ( .NOT. FOUND ) THEN
       CALL GetExtOpt( HcoState%Config, ExtNr, 'NAP scale factor', &
                        OptValSp=ValSp, FOUND=FOUND, RC=RC )
    ENDIF
    IF ( FOUND ) THEN
       MSG = 'Found old definition of CO, POA and/or NAP scale factor! '  // &
             'This version of HEMCO expects species scale factors to be ' // &
             'set as `Scaling_XX` instead of `XX scale factor`. '         // &
             'Please update the FINN settings section accordingly.'
       CALL HCO_ERROR(HcoState%Config%Err,MSG, RC )
       RETURN
    ENDIF

    ! Find matching FINN index for each specified species.
    ! Also get appropriate emission ratios to CO2 (jaf, 10/2/13).
    ! Do this only for species selected for emission calculation. For
    ! all others, keep default values in FINN_EMFAC.
    DO L = 1, tnSpc
       IF ( tHcoIDs(L) < 0 ) CYCLE
       SpcName  = tSpcNames(L)
       N_LUMPED = 0
       Matched  = .FALSE.
       Missing  = .TRUE.

       ! Reduce species if needed
       NCHAR   = LEN(SpcName)
       IF ( NCHAR > 3 ) THEN
          IF ( SpcName(1:3) == 'CO2' ) THEN
             SpcName = 'CO2'
          ELSEIF ( SpcName(1:3) == 'CH4' ) THEN
             SpcName = 'CH4'
          ELSEIF ( SpcName(1:3) == 'CO_' ) THEN
             SpcName = 'CO'
          ELSEIF ( SpcName(1:2) == 'BC' ) THEN
             SpcName = 'BC'
          ELSEIF ( SpcName(1:2) == 'OC' ) THEN
             SpcName = 'OC'
          ENDIF
       ENDIF
       ! For model species NO, the emission factors are taken from FINN
       ! species NOx.  For model species MTPA, the emission factors are
       ! taken from FINN species APINE (BPINE and CARENE will be lumped
       ! into it as well).
       IF ( TRIM(SpcName) == 'POA1' ) SpcName = 'OC'
       IF ( TRIM(SpcName) == 'NAP'  ) SpcName = 'CO'
       IF ( TRIM(SpcName) == 'NO'   ) SpcName = 'NOx'
       IF ( TRIM(SpcName) == 'MTPA' ) SpcName = 'APINE'
       IF ( TRIM(SpcName) == 'Hg0'  ) SpcName = 'CO'
       IF ( TRIM(SpcName) == 'SOAP' ) SpcName = 'CO'

       ! For lumped species, we have to repeat the lookup multiple times,
       ! so use a while loop here.  For example, for species TOLU this will
       ! make sure that FINN species 'TOLU', 'ETBENZ', and 'STYR' are
       ! associated with HEMCO species TOLU. Variable nSpc keeps track of
       ! the total number of species emitted by FINN. All species vectors
       ! (FinnIDs, HcoIDs, SpcNames, SpcScal, etc.) contain nSpc valid
       ! elements.
       DO WHILE ( Missing )

          ! Search for SpcName in FINN
          DO N = 1, N_SPEC
             IF ( TRIM(SpcName) == TRIM(Inst%FINN_SPEC_NAME(N)) ) THEN

                ! Update number of species to be emitted via FINN and
                ! archive all related information in vectors FinnIDs,
                ! HcoIDs, SpcNames, etc.

                ! nSpc is the total number of emitted FINN species. Must
                ! not exceed nSpcMax.
                Inst%nSpc = Inst%nSpc + 1
                IF ( Inst%nSpc > nSpcMax ) THEN
                   MSG = 'nSpc greater than nSpcMax, please increase ' // &
                         'parameter `nSpcMax` in hcox_finn_mod.F90'
                   CALL HCO_ERROR ( HcoState%Config%Err, MSG, RC )
                   RETURN
                ENDIF

                ! Archive corresponding FINN species ID, HEMCO species ID,
                ! scale factor, etc.
                Matched             = .TRUE.
                Inst%FinnIDs(Inst%nSpc)       = N
                Inst%HcoIDs (Inst%nSpc)       = tHcoIDs(L)
                Inst%SpcNames(Inst%nSpc)      = tSpcNames(L)
                Inst%SpcScalFldNme(Inst%nSpc) = tSpcScalFldNme(L)
                Inst%SpcScal(Inst%nSpc)       = tSpcScal(L)

                ! Verbose
                IF ( HcoState%amIRoot ) THEN
                   MSG = '   - FINN species ' // TRIM(Inst%FINN_SPEC_NAME(N)) // &
                         '     will be emitted as ' // TRIM(Inst%SpcNames(Inst%nSpc))
                   CALL HCO_MSG(HcoState%Config%Err,MSG )
                WRITE(MSG,*) '     --> Uniform scale factor : ', Inst%SpcScal(Inst%nSpc)
                CALL HCO_MSG(HcoState%Config%Err,MSG )
                WRITE(MSG,*) '     --> Scale field          : ', TRIM(Inst%SpcScalFldNme(Inst%nSpc))
                CALL HCO_MSG(HcoState%Config%Err,MSG )
                ENDIF

                ! Reset variables
                IS_NMOC    = .FALSE.
                C_MOLEC    = 1d0
                NMOC_RATIO = 0d0

                ! Get emission factor in [kg X]/[kg CO2].
                DO M = 1, N_SPECSTRS
                   TMPNAME = IN_SPEC_NAME(M)
                   IF ( TRIM(Inst%FINN_SPEC_NAME(N)) == TRIM(TMPNAME(5:8)) ) THEN
                      ! First two entries are not species. Also, EMFAC
                      ! is stored as [mole CO2]/[mole X], but we want the
                      ! inverse.  This gives us [mole X]/[mole CO2].
                      ! To convert this to  [kg X]/[kg CO2], we also need
                      ! to adjust for the molecular weights of species X
                      ! and CO2.  The EF ratios of OC and BC are in
                      ! [mole CO2]/[g X], so the adjustment factor is
                      ! calculated slightly differently for those two
                      ! species!
                      IF ( TRIM(Inst%FINN_SPEC_NAME(N)) == 'OC' .OR. &
                           TRIM(Inst%FINN_SPEC_NAME(N)) == 'BC'       ) THEN
                         AdjFact = 1.0_dp / MW_CO2

                      ! Make sure that adjustment factor for CO is always
                      ! computed using the MW of CO. CO might be used as
                      ! proxy for other species (e.g. Hg0), in which case
                      ! we still want to normalize by the MW of CO.
                      ELSEIF ( TRIM(Inst%FINN_SPEC_NAME(N)) == 'CO' ) THEN
                         AdjFact = 28.01_dp / MW_CO2

                      ! Normalize by species' molecular weight.
                      ELSE
                         AdjFact = 1.0_dp / MW_CO2 * &
                                   HcoState%Spc(Inst%HcoIDs(Inst%nSpc))%EmMW_g
                      ENDIF
                      Inst%FINN_EMFAC(N,:) = AdjFact / EMFAC_IN(M,:)
                      IF ( HcoState%amIRoot ) THEN
                         WRITE( MSG, 200 ) TRIM( Inst%FINN_SPEC_NAME(N))
                         CALL HCO_MSG(HcoState%Config%Err,MSG )
                      ENDIF
                      EXIT

                   ! NMOC_EMFAC is converted to [kg NMOC]/[kg CO2].
                   ! Input unit is [mole CO2]/[mole NMOC].
                   ELSE IF ( TRIM(TMPNAME(5:8)) == 'NMOC' ) THEN
                      AdjFact = MW_NMOC / MW_CO2
                      NMOC_EMFAC = AdjFact / EMFAC_IN(M,:)

                   ENDIF
                ENDDO
200             FORMAT( 'Found FINN emission ratio for species ',a5 )

                DO M = 1, N_NMOCSTRS
                   TMPNAME = IN_NMOC_NAME(M)
                   IF ( TRIM(Inst%FINN_SPEC_NAME(N)) == TRIM(TMPNAME) ) THEN
                      ! First two entries are not species
                      NMOC_RATIO = NMOC_RATIO_IN(M,:)
                      IS_NMOC = .TRUE.
                      IF ( HcoState%amIRoot ) THEN
                         WRITE( MSG, 201 ) TRIM( Inst%FINN_SPEC_NAME(N) )
                         CALL HCO_MSG(HcoState%Config%Err,MSG )
                      ENDIF
                      EXIT
                   ENDIF
                ENDDO
201             FORMAT( 'Found FINN NMOC factor for species ',a5 )

                ! Create emission factor for NMOC species
                ! NMOC_EMFAC is [kg NMOC] / [kg CO2]
                ! NMOC_RATIO is [mole X] / [kg NMOC]
                ! To convert NMOC_RATIO to [kg X] / [kg NMOC], we need to
                ! multiply by the MW of X (kg/mol this time).
                ! Most (not all) of these species are carried as atoms C,
                ! so we also multiply here by the number of carbon atoms/molec.
                IF ( IS_NMOC ) THEN
                   DO M = 1, N_EMFAC
                      C_MOLEC              = HcoState%Spc(Inst%HcoIDs(Inst%nSpc))%MolecRatio
                      AdjFact              = HcoState%Spc(Inst%HcoIDs(Inst%nSpc))%EmMW_g
                      Inst%FINN_EMFAC(N,M) = NMOC_EMFAC(M)             * &
                                           ( NMOC_RATIO(M) * C_MOLEC ) * &
                                           ( AdjFact       * 1e-3_hp )
                   ENDDO
                ENDIF
             ENDIF

          ENDDO !N

          ! Update variable Missing. Missing has to be False to exit the
          ! while loop.
          Missing = .FALSE.

          ! For lumped species, we have to repeat the lookup for all
          ! lumped species. For lumped species, we just assign the same
          ! HEMCO species ID to multiple FINN species, so that all of
          ! them will be added to the same model species.

          ! --> TMB is lumped into XYLE
          IF ( Inst%SpcNames(Inst%nSpc) == 'XYLE' ) THEN
             IF ( N_LUMPED == 0 ) THEN
                SpcName  = 'TMB'
                Missing  = .TRUE.
                N_LUMPED = N_LUMPED + 1
             ENDIF
          ENDIF

          ! --> ETBENZ and STYR are lumped into TOLU
          IF ( Inst%SpcNames(Inst%nSpc) == 'TOLU' ) THEN
             IF ( N_LUMPED == 0 ) THEN
                SpcName  = 'ETBENZ'
                Missing  = .TRUE.
                N_LUMPED = N_LUMPED + 1
             ELSEIF ( N_LUMPED == 1 ) THEN
                SpcName  = 'STYR'
                Missing  = .TRUE.
                N_LUMPED = N_LUMPED + 1
             ENDIF
          ENDIF

          ! --> BPINE and CARENE are lumped into MTPA
          IF ( Inst%SpcNames(Inst%nSpc) == 'MTPA' ) THEN
             IF ( N_LUMPED == 0 ) THEN
                SpcName  = 'BPINE'
                Missing  = .TRUE.
                N_LUMPED = N_LUMPED + 1
             ELSEIF ( N_LUMPED == 1 ) THEN
                SpcName  = 'CARENE'
                Missing  = .TRUE.
                N_LUMPED = N_LUMPED + 1
             ENDIF
          ENDIF

       ENDDO !While missing

       ! Error check: we must not specify a species that is not defined
       ! in FINN.
       IF ( .NOT. Matched ) THEN
          MSG = 'Species '// TRIM(SpcName) //' not found in FINN'
          CALL HCO_ERROR(HcoState%Config%Err,MSG, RC )
          RETURN
       ENDIF
    ENDDO !L

    !=======================================================================
    ! Activate this module and the fields of ExtState that it uses
    !=======================================================================

!==============================================================================
! This code is required for the vertical distribution of biomass burning emiss.
! We will keep it here for a future implementation. (mps, 4/24/17)
!    ! Activate met fields required by this extension
!    ExtState%FRAC_OF_PBL%DoUse = .TRUE.
!==============================================================================

    ! Cleanup
    IF ( ALLOCATED(EMFAC_IN        )) DEALLOCATE( EMFAC_IN       )
    IF ( ALLOCATED(NMOC_RATIO_IN   )) DEALLOCATE( NMOC_RATIO_IN  )
    IF ( ALLOCATED(tHcoIDs         )) DEALLOCATE( tHcoIDs        )
    IF ( ALLOCATED(tSpcNames       )) DEALLOCATE( tSpcNames      )
    IF ( ALLOCATED(tSpcScalFldNme  )) DEALLOCATE( tSpcScalFldNme )
    IF ( ALLOCATED(tSpcScal        )) DEALLOCATE( tSpcScal       )

    ! Return w/ success
    Inst => NULL()
    CALL HCO_LEAVE( HcoState%Config%Err,RC )

  END SUBROUTINE HCOX_FINN_Init
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCOX_FINN_Final
!
! !DESCRIPTION: Subroutine HCOX\_FINN\_FINAL deallocates all module
!  arrays.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCOX_FINN_FINAL( ExtState )
!
! !INPUT PARAMETERS:
!
    TYPE(Ext_State),  POINTER       :: ExtState   ! Module options
!
! !REVISION HISTORY:
!  02 Jan 2013 - J. Mao & J. Fisher - Initial version, based on GFED3
!  18 Jun 2014 - C. Keller          - Now a HEMCO extension.
!EOP
!------------------------------------------------------------------------------
!BOC
!
    !=================================================================
    ! HCOX_FINN_FINAL begins here!
    !=================================================================

    CALL InstRemove ( ExtState%FINN )

  END SUBROUTINE HCOX_FINN_Final
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: InstGet
!
! !DESCRIPTION: Subroutine InstGet returns a poiner to the desired instance.
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

    ! Init
    PrevInst => NULL()
    Inst     => NULL()

    ! Get instance. Also archive previous instance.
    CALL InstGet ( Instance, Inst, RC, PrevInst=PrevInst )

    ! Instance-specific deallocation
    IF ( ASSOCIATED(Inst) ) THEN

       ! Pop off instance from list
       IF ( ASSOCIATED(PrevInst) ) THEN

          ! Free pointers
          IF ( ASSOCIATED( Inst%VEGTYP1 ) ) DEALLOCATE( Inst%VEGTYP1 )
          IF ( ASSOCIATED( Inst%VEGTYP2 ) ) DEALLOCATE( Inst%VEGTYP2 )
          IF ( ASSOCIATED( Inst%VEGTYP3 ) ) DEALLOCATE( Inst%VEGTYP3 )
          IF ( ASSOCIATED( Inst%VEGTYP4 ) ) DEALLOCATE( Inst%VEGTYP4 )
          IF ( ASSOCIATED( Inst%VEGTYP5 ) ) DEALLOCATE( Inst%VEGTYP5 )
          IF ( ASSOCIATED( Inst%VEGTYP9 ) ) DEALLOCATE( Inst%VEGTYP9 )

          ! Cleanup module arrays
          IF ( ASSOCIATED( Inst%FINN_EMFAC     )) DEALLOCATE( Inst%FINN_EMFAC     )
          IF ( ASSOCIATED( Inst%FinnIDs        )) DEALLOCATE( Inst%FinnIDs        )
          IF ( ASSOCIATED( Inst%HcoIDs         )) DEALLOCATE( Inst%HcoIDs         )
          IF ( ASSOCIATED( Inst%SpcNames       )) DEALLOCATE( Inst%SpcNames       )
          IF ( ASSOCIATED( Inst%SpcScalFldNme  )) DEALLOCATE( Inst%SpcScalFldNme  )
          IF ( ASSOCIATED( Inst%SpcScal        )) DEALLOCATE( Inst%SpcScal        )
          IF ( ASSOCIATED( Inst%FINN_SPEC_NAME )) DEALLOCATE( Inst%FINN_SPEC_NAME )

          PrevInst%NextInst => Inst%NextInst
       ELSE
          AllInst => Inst%NextInst
       ENDIF
       DEALLOCATE(Inst)
       Inst => NULL()
    ENDIF

   END SUBROUTINE InstRemove
!EOC
END MODULE HCOX_FINN_Mod
