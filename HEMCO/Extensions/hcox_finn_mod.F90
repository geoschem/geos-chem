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
!                                                                             .
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
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
  !=================================================================
  ! MODULE PARAMETERS
  !
  ! N_EMFAC : Number of emission factors per species
  ! N_SPEC  : Number of FINN species
  ! MW_CO2  : Molecular weight of CO2  (g/mol)
  ! MW_NMOC : Molecular weight of NMOC (g/mol). Assumed MW for NMOC
  !           is 68 g/mol.
  !=================================================================
  INTEGER,           PARAMETER   :: N_EMFAC = 6
  INTEGER,           PARAMETER   :: N_SPEC  = 58
  REAL(dp),          PARAMETER   :: MW_CO2  = 44.01_dp
  REAL(dp),          PARAMETER   :: MW_NMOC = 68.00_dp
!
! !PRIVATE TYPES:
!
  !=================================================================
  ! HEMCO VARIABLES 
  !
  ! ExtNr   : Extension number 
  ! UseDay  : True if daily data is used
  !=================================================================
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
  CHARACTER(LEN=31), ALLOCATABLE :: SpcNames(:)
  INTEGER,           ALLOCATABLE :: HcoIDs(:)
  INTEGER,           ALLOCATABLE :: FinnIDs(:)
  CHARACTER(LEN=6),  ALLOCATABLE :: FINN_SPEC_NAME(:)
  
  !=================================================================
  ! SCALE FACTORS 
  !
  ! FINN_EMFAC: emission scale factors for each species and 
  !             emission factor type. The filename of the emissions
  !             emissions factor table is specified in the HEMCO
  !             configuration file. The scale factors are converted
  !             to kg species/kg CO2 when reading them from disk.
  ! COScale   : CO scale factor to account for production from 
  !             VOCs. Read from HEMCO configuration file.
  ! OCPIfrac  : Fraction of OC that converts into hydrophilic OC.
  !             Can be set in HEMCO configuration file (default=0.5)
  ! BCPIfrac  : Fraction of BC that converts into hydrophilic BC.
  !             Can be set in HEMCO configuration file (default=0.2)
  ! POASCALE  : Scale factor for POA. If tracer POA1 is specified, 
  !             emissions are calculated from OC, multiplied by a
  !             POA scale factor that must be specified in the HEMCO
  !             configuration file (POA scale).
  !=================================================================
  REAL(dp),          ALLOCATABLE :: FINN_EMFAC(:,:)
  REAL(sp)                       :: COScale
  REAL(sp)                       :: OCPIfrac
  REAL(sp)                       :: BCPIfrac
  REAL(sp)                       :: POASCALE 

  !=================================================================
  ! DATA ARRAY POINTERS 
  !
  ! These are the pointers to the 6 vegetation type data arrays
  ! specified in the configuration file
  !=================================================================
  REAL(sp),          POINTER     :: VEGTYP1(:,:) => NULL()
  REAL(sp),          POINTER     :: VEGTYP2(:,:) => NULL()
  REAL(sp),          POINTER     :: VEGTYP3(:,:) => NULL()
  REAL(sp),          POINTER     :: VEGTYP4(:,:) => NULL()
  REAL(sp),          POINTER     :: VEGTYP5(:,:) => NULL()
  REAL(sp),          POINTER     :: VEGTYP9(:,:) => NULL()

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
  SUBROUTINE HCOX_FINN_Run( am_I_Root, ExtState, HcoState, RC )
!
! !USES:
!
    USE HCO_EmisList_mod,  ONLY : HCO_GetPtr
    USE HCO_FluxArr_mod,   ONLY : HCO_EmisAdd
    USE HCO_State_mod,     ONLY : HCO_GetHcoID
    USE HCO_Clock_mod,     ONLY : HcoClock_Get
    USE HCO_Clock_mod,     ONLY : HcoClock_NewMonth, HcoClock_NewDay
!
! !INPUT PARAMETERS:
!
    LOGICAL,         INTENT(IN   )  :: am_I_Root  ! root CPU?
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
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER             :: N, NF, ID, HcoID
    LOGICAL, SAVE       :: FIRST = .TRUE.
    LOGICAL             :: DoRepeat
    INTEGER             :: Cnt
    CHARACTER(LEN=31)   :: PREFIX, FLDNME
    INTEGER             :: NDAYS, cYYYY, cMM, cDD
    REAL(dp)            :: TOTAL
    CHARACTER(LEN=255)  :: MSG

    ! Arrays
    REAL(hp), TARGET    :: SpcArr(HcoState%NX,HcoState%NY)
    REAL(hp), TARGET    :: TypArr(HcoState%NX,HcoState%NY)

    ! Pointers
    REAL(sp), POINTER   :: THISTYP(:,:) => NULL()

    !=======================================================================
    ! HCOX_FINN_Run begins here!
    !=======================================================================

    ! Return if extension disabled 
    IF ( .NOT. ExtState%FINN ) RETURN

    ! Enter 
    CALL HCO_ENTER( 'HCOX_FINN_RUN (hcox_finn_mod.F90)', RC ) 
    IF ( RC /= HCO_SUCCESS ) RETURN

    !-----------------------------------------------------------------------
    ! Get pointers to data arrays 
    !-----------------------------------------------------------------------
    IF ( FIRST ) THEN
       IF ( UseDay ) THEN
          PREFIX = 'FINN_DAILY_'
       ELSE
          PREFIX = 'FINN_'
       ENDIF
   
       FLDNME = TRIM(PREFIX) // 'VEGTYP1'
       CALL HCO_GetPtr( am_I_Root, TRIM(FLDNME), VEGTYP1, RC )
       IF ( RC /= HCO_SUCCESS ) RETURN

       FLDNME = TRIM(PREFIX) // 'VEGTYP2'
       CALL HCO_GetPtr( am_I_Root, TRIM(FLDNME), VEGTYP2, RC )
       IF ( RC /= HCO_SUCCESS ) RETURN
   
       FLDNME = TRIM(PREFIX) // 'VEGTYP3'
       CALL HCO_GetPtr( am_I_Root, TRIM(FLDNME), VEGTYP3, RC )
       IF ( RC /= HCO_SUCCESS ) RETURN
   
       FLDNME = TRIM(PREFIX) // 'VEGTYP4'
       CALL HCO_GetPtr( am_I_Root, TRIM(FLDNME), VEGTYP4, RC )
       IF ( RC /= HCO_SUCCESS ) RETURN

       FLDNME = TRIM(PREFIX) // 'VEGTYP5'
       CALL HCO_GetPtr( am_I_Root, TRIM(FLDNME), VEGTYP5, RC )
       IF ( RC /= HCO_SUCCESS ) RETURN
   
       FLDNME = TRIM(PREFIX) // 'VEGTYP9'
       CALL HCO_GetPtr( am_I_Root, TRIM(FLDNME), VEGTYP9, RC )
       IF ( RC /= HCO_SUCCESS ) RETURN

       FIRST = .FALSE.
    ENDIF

    ! For logfile
    IF ( am_I_Root ) THEN
       IF ( UseDay ) THEN
          IF ( HcoClock_NewDay( .TRUE. ) ) THEN
             CALL HcoClock_Get( cYYYY = cYYYY, cMM=cMM, cDD=cDD, RC=RC )
             IF ( RC/=HCO_SUCCESS ) RETURN
             WRITE(MSG, 100) cYYYY, cMM, cDD
             CALL HCO_MSG(MSG)
100          FORMAT( 'FINN daily emissions for year, month, day: ', &
                      i4, '/', i2.2, '/', i2.2 )
          ENDIF
       ELSE 
          IF ( HcoClock_NewMonth( .TRUE. ) ) THEN
             CALL HcoClock_Get( cYYYY = cYYYY, cMM=cMM, LMD=NDAYS, RC=RC)
             IF ( RC/=HCO_SUCCESS ) RETURN
             WRITE(MSG, 110) cYYYY, cMM
             CALL HCO_MSG(MSG)
110          FORMAT( 'FINN monthly emissions for year, month: ', &
                      i4, '/', i2.2 )
          ENDIF
       ENDIF
    ENDIF

    !-----------------------------------------------------------------------
    ! Calculate emissions for all selected species
    !-----------------------------------------------------------------------

    ! Loop over FINN species
    DO N = 1, N_SPEC

       ! ID is the index of the suite of defined species. 
       ID = FinnIDs(N)
       IF ( ID <= 0 ) CYCLE

       ! HcoID is the species index in HEMCO
       HcoID = HcoIDs(ID)
       IF ( HcoID < 0 ) CYCLE

       ! Species with no emission factor have FINN_EMFAC=0
       IF ( MAXVAL(FINN_EMFAC(N,:)) <= 0.0_hp ) CYCLE

       ! SpcArr are the total biomass burning emissions for this
       ! species. TypArr are the emissions from a given vegetation type. 
       SpcArr = 0.0_hp

       ! Calculate emissions for all source types
       DO NF = 1, N_EMFAC
         
          ! Select emission factor array
          IF ( NF == 1 ) THEN
             THISTYP => VEGTYP1 
          ELSEIF ( NF == 2 ) THEN
             THISTYP => VEGTYP2 
          ELSEIF ( NF == 3 ) THEN
             THISTYP => VEGTYP3
          ELSEIF ( NF == 4 ) THEN
             THISTYP => VEGTYP4
          ELSEIF ( NF == 5 ) THEN
             THISTYP => VEGTYP5
          ELSEIF ( NF == 6 ) THEN
             THISTYP => VEGTYP9
          ELSE
             CALL HCO_ERROR ( 'Undefined emission factor', RC )
             RETURN
          ENDIF

          ! Multiply CO2 emissions by appropriate ratio for each land
          ! type and sum to get total emissions for the species on the
          ! native grid - emissions are in [kg CO2/m2/s[. FINN_EMFAC is
          ! in [kg X]/[kg CO2].
          TypArr(:,:) = THISTYP(:,:) * FINN_EMFAC(N,NF)

          ! TODO: Add to diagnostics here

          ! Add to species array
          SpcArr = SpcArr + TypArr
       ENDDO !NF

       ! Apply species specific scale factors
       SELECT CASE ( SpcNames(N) )
          CASE ( 'CO' )
             SpcArr = SpcArr * COScale
          CASE ( 'OCPI' )
             SpcArr = SpcArr * OCPIfrac
          CASE ( 'OCPO' )
             SpcArr = SpcArr * (1.0_sp - OCPIfrac)
          CASE ( 'BCPI' )
             SpcArr = SpcArr * BCPIfrac
          CASE ( 'BCPO' )
             SpcArr = SpcArr * (1.0_sp - BCPIfrac)
          CASE ( 'POA1' )
             SpcArr = POASCALE * SpcArr
       END SELECT

       ! Add flux to HEMCO emission array
       CALL HCO_EmisAdd( am_I_Root, HcoState,    SpcArr, HcoID, &
                         RC,        ExtNr=ExtNr, Cat=-1, Hier=-1 ) 
       IF ( RC /= HCO_SUCCESS ) THEN
          MSG = 'HCO_EmisAdd error: ' // TRIM(HcoState%Spc(HcoID)%SpcName)
          CALL HCO_ERROR( MSG, RC )
          RETURN 
       ENDIF
 
       ! Write out total (daily or monthly) emissions to log-file
       IF ( am_I_Root ) THEN
          IF ( UseDay ) THEN
             IF ( HcoClock_NewDay( .TRUE. ) ) THEN
                TOTAL = SUM(SpcArr(:,:)*HcoState%Grid%AREA_M2%Val(:,:))
                TOTAL = TOTAL * 86400.0_hp * 1e-9_hp
                WRITE(MSG, 120) HcoState%Spc(HcoID)%SpcName, TOTAL
                CALL HCO_MSG(MSG)
120             FORMAT( 'SUM biomass ', a4,1x,': ', f11.4,1x,'[Tg]' )
             ENDIF
          ELSE
             IF ( HcoClock_NewMonth( .TRUE. ) ) THEN
                TOTAL = SUM(SpcArr(:,:)*HcoState%Grid%AREA_M2%Val(:,:))
                TOTAL = TOTAL * NDAYS * 86400.0_hp * 1e-9_hp
                WRITE(MSG, 130) HcoState%Spc(HcoID)%SpcName, TOTAL
                CALL HCO_MSG(MSG)
130             FORMAT( 'SUM biomass ', a4,1x,': ', f11.4,1x,'[Tg]' )
             ENDIF
          ENDIF
       ENDIF
  
    ENDDO !N

    ! Nullify pointers
    THISTYP   => NULL()

    ! Leave w/ success
    CALL HCO_LEAVE ( RC )

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
  SUBROUTINE HCOX_FINN_Init( am_I_Root, HcoState, ExtName, ExtState, RC ) 
!
! !USES:
!
    USE HCO_State_Mod,   ONLY : HCO_GetHcoID
    USE HCO_State_Mod,   ONLY : HCO_GetExtHcoID
    USE HCO_ExtList_Mod, ONLY : GetExtNr, GetExtOpt
!
! !INPUT PARAMETERS:
!
    LOGICAL,          INTENT(IN   )  :: am_I_Root   ! root CPU?
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
    INTEGER               :: N_SPEC_EMFAC  ! # of CO2 file emission species
    INTEGER               :: N_NMOC        ! # of VOC file NMOC ratios
    INTEGER               :: IU_FILE, L, N_LUMPED, tmpNr
    INTEGER               :: AS, IOS, M, N, NDUM
    INTEGER               :: N_SPECSTRS, N_NMOCSTRS
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

    ! Arrays
    REAL(dp), ALLOCATABLE :: EMFAC_IN(:,:)
    REAL(dp), ALLOCATABLE :: NMOC_RATIO_IN(:,:)
    REAL*8                :: NMOC_EMFAC(N_EMFAC), NMOC_RATIO(N_EMFAC)


    !=======================================================================
    ! HCOX_FINN_INIT begins here!
    !=======================================================================

    ! Extension Nr.
    ExtNr = GetExtNr( TRIM(ExtName) )
    IF ( ExtNr <= 0 ) RETURN
 
    ! Enter 
    CALL HCO_ENTER( 'HCOX_FINN_INIT (hcox_finn_mod.F90)', RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    !-----------------------------------------------------------------------
    ! Get settings 
    ! The CO scale factor (to account for oxidation from VOCs) as well as 
    ! the speciation of carbon aerosols into hydrophilic and hydrophobic
    ! fractions can be specified in the configuration file, e.g.:
    ! 100     GFED3           : on    NO/CO/OCPI/OCPO/BCPI/BCPO
    !     --> CO scale factor :       1.05
    !     --> hydrophilic BC  :       0.2
    !     --> hydrophilic OC  :       0.5
    !
    ! Setting these values is optional and default values are applied if 
    ! they are not specified. The values only take effect if the
    ! corresponding species (CO, BCPI/BCPO, OCPI/OCPO) are listed as species
    ! to be used.
    !----------------------------------------------------------------------- 

    ! Try to read CO scale factor. Defaults to 1.0
    CALL GetExtOpt ( ExtNr, 'CO scale factor', &
                     OptValSp=ValSp, FOUND=FOUND, RC=RC )
    IF ( RC /= HCO_SUCCESS ) RETURN
    IF ( .NOT. FOUND ) THEN
       COScale = 1.0
    ELSE
       COScale = ValSp
    ENDIF

    ! Try to read hydrophilic fractions of BC. Defaults to 0.2.
    CALL GetExtOpt ( ExtNr, 'hydrophilic BC', &
                     OptValSp=ValSp, FOUND=FOUND, RC=RC )
    IF ( RC /= HCO_SUCCESS ) RETURN
    IF ( .NOT. FOUND ) THEN
       BCPIfrac = 0.2
    ELSE
       BCPIfrac = ValSp
    ENDIF

    ! Try to read hydrophilic fractions of OC. Defaults to 0.5.
    CALL GetExtOpt ( ExtNr, 'hydrophilic OC', &
                     OptValSp=ValSp, FOUND=FOUND, RC=RC )
    IF ( RC /= HCO_SUCCESS ) RETURN
    IF ( .NOT. FOUND ) THEN
       OCPIfrac = 0.5
    ELSE
       OCPIfrac = ValSp
    ENDIF

    ! Error check: OCPIfrac and BCPI frac must be between 0 and 1
    IF ( OCPIfrac < 0.0_sp .OR. OCPIfrac > 1.0_sp .OR. &
         BCPIfrac < 0.0_sp .OR. BCPIfrac > 1.0_sp     ) THEN
       WRITE(MSG,*) 'hydrophilic fractions must be between 0-1: ', &
          OCPIfrac, BCPIfrac
       CALL HCO_ERROR( MSG, RC )
       RETURN
    ENDIF
 
    ! Use daily data?
    CALL GetExtOpt ( ExtNr, 'FINN_daily', OptValBool=UseDay, FOUND=FOUND, RC=RC ) 
    IF ( RC /= HCO_SUCCESS ) RETURN
    IF ( .NOT. FOUND ) THEN 
       UseDay = .FALSE.
    ENDIF

    !----------------------------------------------------------------------- 
    ! Allocate arrays
    !----------------------------------------------------------------------- 

    ! FINN species names
    ALLOCATE ( FINN_SPEC_NAME ( N_SPEC ), STAT=AS )
    IF ( AS/=0 ) THEN
       CALL HCO_ERROR( 'Cannot allocate FINN_SPEC_NAME', RC )
       RETURN
    ENDIF
    FINN_SPEC_NAME = ''

    ! Allocate scale factors table: FINN_EMFAC holds the species/CO2
    ! scale factors for all FINN species.
    ALLOCATE ( FINN_EMFAC ( N_SPEC, N_EMFAC ), STAT=AS )
    IF ( AS/=0 ) THEN
       CALL HCO_ERROR( 'Cannot allocate FINN_EMFAC', RC )
       RETURN
    ENDIF
    FINN_EMFAC = 0.0_dp

    ! FinnIDs maps the FINN species onto the suite of specified species.
    ALLOCATE ( FinnIDs ( N_SPEC ), STAT=AS )
    IF ( AS/=0 ) THEN
       CALL HCO_ERROR( 'Cannot allocate FinnIDs', RC )
       RETURN
    ENDIF
    FinnIDs = -1

    !----------------------------------------------------------------------- 
    ! Define FINN species names
    !----------------------------------------------------------------------- 
      
    ! Species listed in emission factor ratios (CO2/X) table (except NMOC,
    ! which is speciated as specified in the VOC speciation table).
    FINN_SPEC_NAME(1)  = 'CO2'
    FINN_SPEC_NAME(2)  = 'CO'
    FINN_SPEC_NAME(3)  = 'CH4'
    FINN_SPEC_NAME(4)  = 'NOx'
    FINN_SPEC_NAME(5)  = 'SO2'
    FINN_SPEC_NAME(6)  = 'OC'
    FINN_SPEC_NAME(7)  = 'BC'
    FINN_SPEC_NAME(8)  = 'NH3'
    FINN_SPEC_NAME(9)  = 'NO'    ! Currently not used
    FINN_SPEC_NAME(10) = 'NO2'   ! Currently not used

    ! Species listed in VOC speciation table
    FINN_SPEC_NAME(11) = 'ACET'
    FINN_SPEC_NAME(12) = 'ACTA'   ! Not currently emitted by BB in GC
    FINN_SPEC_NAME(13) = 'ALD2'
    FINN_SPEC_NAME(14) = 'ALK4'
    FINN_SPEC_NAME(15) = 'APINE'  ! Currently lumped into MTPA
    FINN_SPEC_NAME(16) = 'AROM'   ! Currently not used
    FINN_SPEC_NAME(17) = 'BENZ'
    FINN_SPEC_NAME(18) = 'BPINE'  ! Currently lumped into MTPA
    FINN_SPEC_NAME(19) = 'C2H2'
    FINN_SPEC_NAME(20) = 'C2H4'
    FINN_SPEC_NAME(21) = 'C2H6'
    FINN_SPEC_NAME(22) = 'C3H8'
    FINN_SPEC_NAME(23) = 'CARENE' ! Currently lumped into MTPA
    FINN_SPEC_NAME(24) = 'CH2Br2'
    FINN_SPEC_NAME(25) = 'CH2O'
    FINN_SPEC_NAME(26) = 'CH3Br'
    FINN_SPEC_NAME(27) = 'CH3CN'
    FINN_SPEC_NAME(28) = 'CH3I'
    FINN_SPEC_NAME(29) = 'DMS'
    FINN_SPEC_NAME(30) = 'EOH'    ! Not currently emitted in GC
    FINN_SPEC_NAME(31) = 'ETBENZ' ! Currently lumped with TOLU
    FINN_SPEC_NAME(32) = 'FUR'    ! Currently not used
    FINN_SPEC_NAME(33) = 'GLYC'
    FINN_SPEC_NAME(34) = 'GLYX'
    FINN_SPEC_NAME(35) = 'HAC'
    FINN_SPEC_NAME(36) = 'HCN'    ! Not currently emitted in GC
    FINN_SPEC_NAME(37) = 'HCOOH'  ! Not currently emitted by BB in GC
    FINN_SPEC_NAME(38) = 'HNO2'   ! Not currently emitted in GC
    FINN_SPEC_NAME(39) = 'ISOP'   ! Not currently emitted by BB in GC
    FINN_SPEC_NAME(40) = 'LIMO'
    FINN_SPEC_NAME(41) = 'MACR'   ! Not currently emitted in GC
    FINN_SPEC_NAME(42) = 'MEK'
    FINN_SPEC_NAME(43) = 'MGLY'
    FINN_SPEC_NAME(44) = 'MNO3'
    FINN_SPEC_NAME(45) = 'MOH'    ! Not currently emitted in GC
    FINN_SPEC_NAME(46) = 'MTPO'   ! Not currently emitted in GC
    FINN_SPEC_NAME(47) = 'MVK'    ! Not currently emitted in GC
    FINN_SPEC_NAME(48) = 'PRPE'
    FINN_SPEC_NAME(49) = 'R4N2'   ! Not currently emitted in GC
    FINN_SPEC_NAME(50) = 'RCHO'   ! Not currently emitted by BB in GC
    FINN_SPEC_NAME(51) = 'RCOOH'  ! Currently not used
    FINN_SPEC_NAME(52) = 'ROH'    ! Currently not used
    FINN_SPEC_NAME(53) = 'SESQ'   ! Currently not used
    FINN_SPEC_NAME(54) = 'STYR'   ! Currently lumped with TOLU
    FINN_SPEC_NAME(55) = 'TMB'    ! Currently lumped with XYLE
    FINN_SPEC_NAME(56) = 'TOLU'
    FINN_SPEC_NAME(57) = 'XYLE'
    FINN_SPEC_NAME(58) = 'H2'     ! Currently not used

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
    IF ( am_I_Root ) THEN
       MSG = 'Use FINN extension'
       CALL HCO_MSG( MSG, SEP1='-' )
       WRITE(MSG,*) '   - Use daily data          : ', UseDay
       CALL HCO_MSG( MSG )
    ENDIF

    ! Get HEMCO species IDs of all species specified in configuration file
    CALL HCO_GetExtHcoID( HcoState, ExtNr, HcoIDs, SpcNames, nSpc, RC)
    IF ( RC /= HCO_SUCCESS ) RETURN
    IF ( nSpc == 0 ) THEN
       MSG = 'No FINN species specified'
       CALL HCO_ERROR ( MSG, RC ) 
       RETURN
    ENDIF

    ! Find matching FINN index for each specified species. 
    ! Also get appropriate emission ratios to CO2 (jaf, 10/2/13).
    ! Do this only for species selected for emission calculation. For
    ! all others, keep default values in FINN_EMFAC.
    DO L = 1, nSpc
       IF ( HcoIDs(L) < 0 ) CYCLE
       SpcName  = SpcNames(L)
       N_LUMPED = 0
       Matched  = .FALSE.
       Missing  = .TRUE.

       ! For model species NO, the emission factors are taken from FINN  
       ! species NOx.  For model species MTPA, the emission factors are 
       ! taken from FINN species APINE (BPINE and CARENE will be lumped 
       ! into it as well).
       ! Also reduce 'CO2bb' to 'CO2'.
       SELECT CASE ( TRIM(SpcName) )
          CASE ( 'NO' )
            SpcName = 'NOx'
          CASE ('MTPA' )
            SpcName = 'APINE'
          CASE ('CO2bb' )
            SpcName = 'CO2'
          CASE ('CH4_bb', 'CH4_tot' )
            SpcName = 'CH4'
          CASE ( 'BC', 'BCPI', 'BCPO' )
             SpcName = 'BC'
          CASE ( 'OC', 'OCPI', 'OCPO', 'POA1' )
             SpcName = 'OC'
       END SELECT

       ! For lumped species, we have to repeat the lookup multiple times,
       ! so use a while loop here.  For example, for species TOLU this will 
       ! make sure that FINN species 'TOLU', 'ETBENZ', and 'STYR' are 
       ! associated with HEMCO species TOLU.
       DO WHILE ( Missing )

          ! Search for SpcName in FINN
          DO N = 1, N_SPEC 
             IF ( TRIM(SpcName) == TRIM(FINN_SPEC_NAME(N)) ) THEN
                FinnIDs(N) = L
                Matched    = .TRUE.
  
                IF ( am_I_Root ) THEN 
                   MSG = '   - FINN species ' // TRIM(FINN_SPEC_NAME(N)) // &
                         '     will be emitted as ' // TRIM(SpcNames(L))
                   CALL HCO_MSG( MSG )
                ENDIF   

                ! Reset variables
                IS_NMOC    = .FALSE.
                C_MOLEC    = 1d0
                NMOC_RATIO = 0d0

                ! Get emission factor in [kg X]/[kg CO2]. 
                DO M = 1, N_SPECSTRS
                   TMPNAME = IN_SPEC_NAME(M)
                   IF ( TRIM(FINN_SPEC_NAME(N)) == TRIM(TMPNAME(5:8)) ) THEN
                      ! First two entries are not species. Also, EMFAC 
                      ! is stored as [mole CO2]/[mole X], but we want the 
                      ! inverse.  This gives us [mole X]/[mole CO2].  
                      ! To convert this to  [kg X]/[kg CO2], we also need 
                      ! to adjust for the molecular weights of species X 
                      ! and CO2.  The EF ratios of OC and BC are in 
                      ! [mole CO2]/[g X], so the adjustment factor is 
                      ! calculated slightly differently for those two 
                      ! species!
                      IF ( TRIM(FINN_SPEC_NAME(N)) == 'OC' .OR. &
                           TRIM(FINN_SPEC_NAME(N)) == 'BC'       ) THEN
                         AdjFact = 1.0_dp / MW_CO2
                      ELSE
                         AdjFact = 1.0_dp / MW_CO2 * &
                                   HcoState%Spc(HcoIDs(L))%EmMW_g
                      ENDIF
                      FINN_EMFAC(N,:) = AdjFact / EMFAC_IN(M,:)
                      IF ( am_I_Root ) THEN
                         WRITE( MSG, 200 ) TRIM( FINN_SPEC_NAME(N)) 
                         CALL HCO_MSG( MSG )
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
                   IF ( TRIM(FINN_SPEC_NAME(N)) == TRIM(TMPNAME) ) THEN
                      ! First two entries are not species
                      NMOC_RATIO = NMOC_RATIO_IN(M,:)
                      IS_NMOC = .TRUE.
                      IF ( am_I_Root ) THEN
                         WRITE( MSG, 201 ) TRIM( FINN_SPEC_NAME(N) )
                         CALL HCO_MSG( MSG )
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
                      C_MOLEC         = HcoState%Spc(HcoIDs(L))%MolecRatio
                      AdjFact         = HcoState%Spc(HcoIDs(L))%EmMW_g
                      FINN_EMFAC(N,M) = NMOC_EMFAC(M)             * &
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
          IF ( SpcNames(L) == 'XYLE' ) THEN
             IF ( N_LUMPED == 0 ) THEN
                SpcName  = 'TMB'
                Missing  = .TRUE.
                N_LUMPED = N_LUMPED + 1
             ENDIF
          ENDIF

          ! --> ETBENZ and STYR are lumped into TOLU
          IF ( SpcNames(L) == 'TOLU' ) THEN
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
          IF ( SpcNames(L) == 'MTPA' ) THEN
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

       ! For tracer POA1, the POA scale factor must be defined in the HEMCO 
       ! configuration file
       IF ( TRIM(SpcNames(L)) == 'POA1' ) THEN
          CALL GetExtOpt ( ExtNr, 'POA scale', &
                           OptValSp=ValSp, FOUND=FOUND, RC=RC )
          IF ( RC /= HCO_SUCCESS ) RETURN
          IF ( .NOT. FOUND ) THEN
             MSG = 'You must specify a POA scale factor for species POA1'
             CALL HCO_ERROR( MSG, RC )
             RETURN
          ENDIF
          POASCALE = ValSp
       ENDIF

       ! Error check: we must not specify a species that is not defined
       ! in FINN.
       IF ( .NOT. Matched ) THEN
          MSG = 'Species '// TRIM(SpcName) //' not found in FINN'
          CALL HCO_ERROR( MSG, RC )
          RETURN
       ENDIF
    ENDDO !L

    ! Enable module
    ExtState%FINN = .TRUE.

    ! Cleanup
    IF ( ALLOCATED(EMFAC_IN     ) ) DEALLOCATE( EMFAC_IN      )
    IF ( ALLOCATED(NMOC_RATIO_IN) ) DEALLOCATE( NMOC_RATIO_IN )
 
    ! Return w/ success
    CALL HCO_LEAVE( RC ) 
 
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
  SUBROUTINE HCOX_FINN_FINAL()
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

    ! Free pointers
    VEGTYP1   => NULL()
    VEGTYP2   => NULL()
    VEGTYP3   => NULL()
    VEGTYP4   => NULL()
    VEGTYP5   => NULL()
    VEGTYP9   => NULL()

    ! Cleanup module arrays
    IF ( ALLOCATED( FINN_EMFAC     )) DEALLOCATE( FINN_EMFAC     )
    IF ( ALLOCATED( FinnIDs        )) DEALLOCATE( FinnIDs        )
    IF ( ALLOCATED( HcoIDs         )) DEALLOCATE( HcoIDs         )
    IF ( ALLOCATED( SpcNames       )) DEALLOCATE( SpcNames       )
    IF ( ALLOCATED( FINN_SPEC_NAME )) DEALLOCATE( FINN_SPEC_NAME )

  END SUBROUTINE HCOX_FINN_Final
!EOC
END MODULE HCOX_FINN_Mod
