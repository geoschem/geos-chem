!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: hcox_finn_mod
!
! !DESCRIPTION: Module HCOX\_FINN\_MOD contains routines and variables to 
! calculate FINN biomass burning emissions in HEMCO.
!
! !INTERFACE: 
!
      MODULE HCOX_FINN_MOD
!
! !USES:
! 
      USE HCO_ERROR_MOD
      USE HCO_DIAGN_MOD
      USE HCO_STATE_MOD,        ONLY : HCO_State
      USE HCOX_State_MOD,       ONLY : Ext_State

      IMPLICIT NONE
      PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
      PUBLIC :: HCOX_FINN_INIT
      PUBLIC :: HCOX_FINN_RUN
      PUBLIC :: HCOX_FINN_FINAL
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
      INTEGER,          PARAMETER   :: N_EMFAC = 6
      INTEGER,          PARAMETER   :: N_SPEC  = 58
      REAL(dp),         PARAMETER   :: MW_CO2  = 44.01_dp
      REAL(dp),         PARAMETER   :: MW_NMOC = 68.00_dp
!
! !PRIVATE TYPES:
!
      !=================================================================
      ! HEMCO VARIABLES 
      !
      ! ExtNr   : Extension number 
      ! UseDay  : True if daily data is used
      !=================================================================
      INTEGER                     :: ExtNr
      LOGICAL                     :: UseDay

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
      !=================================================================
      REAL(dp),          ALLOCATABLE :: FINN_EMFAC(:,:)

      !=================================================================
      ! DATA ARRAY POINTERS 
      !
      ! These are the pointers to the 6 vegetation type data arrays
      ! specified in the configuration file
      !=================================================================
      REAL(hp), POINTER   :: VEGTYP1(:,:) => NULL()
      REAL(hp), POINTER   :: VEGTYP2(:,:) => NULL()
      REAL(hp), POINTER   :: VEGTYP3(:,:) => NULL()
      REAL(hp), POINTER   :: VEGTYP4(:,:) => NULL()
      REAL(hp), POINTER   :: VEGTYP5(:,:) => NULL()
      REAL(hp), POINTER   :: VEGTYP9(:,:) => NULL()

      CONTAINS
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HcoX_FINN_RUN
!
! !DESCRIPTION: Subroutine HCOX\_FINN\_RUN computes the FINN biomass
! burning emissions for the current date.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE HCOX_FINN_RUN ( am_I_Root, ExtState, HcoState, RC )
!
! !USES:
!
      USE HCO_EMISLIST_MOD,  ONLY : EmisList_GetDataArr
      USE HCO_FLUXARR_MOD,   ONLY : HCO_EmisAdd
      USE HCO_STATE_MOD,     ONLY : HCO_GetHcoID
      USE HCO_CLOCK_MOD,     ONLY : HcoClock_Get
      USE HCO_CLOCK_MOD,     ONLY : HcoClock_NewMonth, HcoClock_NewDay
!
! !INPUT/OUTPUT PARAMETERS:
!
      LOGICAL,         INTENT(IN   )  :: am_I_Root  ! root CPU?
      TYPE(HCO_State), POINTER        :: HcoState   ! Output obj
      TYPE(Ext_State), POINTER        :: ExtState   ! Module options  
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
      INTEGER             :: N, NF, ID, HcoID
      REAL(hp), POINTER   :: THISTYP(:,:) => NULL()
      REAL(hp), POINTER   :: Arr2D  (:,:) => NULL()
      
      REAL(hp), TARGET    :: SpcArr(HcoState%NX,HcoState%NY)
      REAL(hp), TARGET    :: TypArr(HcoState%NX,HcoState%NY)

      LOGICAL, SAVE       :: FIRST = .TRUE.
 
      ! For OC/BC splitting:
      LOGICAL             :: DoRepeat
      INTEGER             :: Cnt

      ! Get field names
      CHARACTER(LEN=31)   :: PREFIX, FLDNME

      ! Write totals to log file 
      INTEGER             :: NDAYS, cYYYY, cMM, cDD
      REAL(dp)            :: TOTAL
      CHARACTER(LEN=255)  :: MSG
 
      !=================================================================
      ! HCOX_FINN_RUN begins here!
      !=================================================================

      ! Return if extension disabled 
      IF ( .NOT. ExtState%FINN ) RETURN

      ! Enter 
      CALL HCO_ENTER ( 'HCOX_FINN_RUN (HcoX_FINN_Mod.F90)', RC ) 
      IF ( RC /= HCO_SUCCESS ) RETURN

      !-----------------------------------------------------------------
      ! Get pointers to data arrays 
      !-----------------------------------------------------------------
      IF ( FIRST ) THEN
         IF ( UseDay ) THEN
            PREFIX = 'FINN_DAILY_'
         ELSE
            PREFIX = 'FINN_'
         ENDIF
   
         FLDNME = TRIM(PREFIX) // 'VEGTYP1'
         CALL EmisList_GetDataArr( am_I_Root, TRIM(FLDNME), VEGTYP1, RC )
         IF ( RC /= HCO_SUCCESS ) RETURN
   
         FLDNME = TRIM(PREFIX) // 'VEGTYP2'
         CALL EmisList_GetDataArr( am_I_Root, TRIM(FLDNME), VEGTYP2, RC )
         IF ( RC /= HCO_SUCCESS ) RETURN
   
         FLDNME = TRIM(PREFIX) // 'VEGTYP3'
         CALL EmisList_GetDataArr( am_I_Root, TRIM(FLDNME), VEGTYP3, RC )
         IF ( RC /= HCO_SUCCESS ) RETURN
   
         FLDNME = TRIM(PREFIX) // 'VEGTYP4'
         CALL EmisList_GetDataArr( am_I_Root, TRIM(FLDNME), VEGTYP4, RC )
         IF ( RC /= HCO_SUCCESS ) RETURN

         FLDNME = TRIM(PREFIX) // 'VEGTYP5'
         CALL EmisList_GetDataArr( am_I_Root, TRIM(FLDNME), VEGTYP5, RC )
         IF ( RC /= HCO_SUCCESS ) RETURN
   
         FLDNME = TRIM(PREFIX) // 'VEGTYP9'
         CALL EmisList_GetDataArr( am_I_Root, TRIM(FLDNME), VEGTYP9, RC )
         IF ( RC /= HCO_SUCCESS ) RETURN

         FIRST = .FALSE.
      ENDIF

      ! For logfile
      IF ( UseDay ) THEN
         IF ( HcoClock_NewDay() ) THEN
            CALL HcoClock_Get( cYYYY = cYYYY, cMM=cMM, cDD=cDD, RC=RC )
            IF ( RC/=HCO_SUCCESS ) RETURN
            WRITE(MSG, 100) cYYYY, cMM, cDD
            CALL HCO_MSG(MSG)
 100        FORMAT( 'FINN daily emissions for year, month, day: ', &
                     i4, '/', i2.2, '/', i2.2 )
         ENDIF
      ELSE 
         IF ( HcoClock_NewMonth() ) THEN
            CALL HcoClock_Get( cYYYY = cYYYY, cMM=cMM, LMD=NDAYS, RC=RC)
            IF ( RC/=HCO_SUCCESS ) RETURN
            WRITE(MSG, 110) cYYYY, cMM
            CALL HCO_MSG(MSG)
 110        FORMAT( 'FINN monthly emissions for year, month: ', &
                     i4, '/', i2.2 )
         ENDIF
      ENDIF  

      !-----------------------------------------------------------------
      ! Calculate emissions for all selected species
      !-----------------------------------------------------------------

      ! Loop over FINN species
      DO N = 1, N_SPEC

         ! ID is the index of the suite of defined species. 
         ID = FinnIDs(N)
         IF ( ID <= 0 ) CYCLE

         ! HcoID is the species index in the atm. model
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

         ! Add flux to HEMCO emission array
         CALL HCO_EmisAdd( HcoState, SpcArr, HcoID, RC ) 
         IF ( RC /= HCO_SUCCESS ) RETURN
   
         ! Write out total (daily or monthly) emissions to log-file
         IF ( UseDay ) THEN
            IF ( HcoClock_NewDay() ) THEN
               TOTAL = SUM(SpcArr(:,:)*HcoState%Grid%Area_M2(:,:))
               TOTAL = TOTAL * 86400.0_hp * 1e-9_hp
               WRITE(MSG, 120) HcoState%Spc(HcoID)%SpcName, TOTAL
               CALL HCO_MSG(MSG)
 120           FORMAT( 'SUM biomass ', a4,1x,': ', f11.4,1x,'[Tg]' )
            ENDIF
         ELSE
            IF ( HcoClock_NewMonth() ) THEN
               TOTAL = SUM(SpcArr(:,:)*HcoState%Grid%Area_M2(:,:))
               TOTAL = TOTAL * NDAYS * 86400.0_hp * 1e-9_hp
               WRITE(MSG, 130) HcoState%Spc(HcoID)%SpcName, TOTAL
               CALL HCO_MSG(MSG)
 130           FORMAT( 'SUM biomass ', a4,1x,': ', f11.4,1x,'[Tg]' )
            ENDIF
         ENDIF 

         ! Eventually update diagnostics
         IF ( Diagn_AutoFillLevelDefined(2) ) THEN
            Arr2D => SpcArr
            CALL Diagn_Update(am_I_Root, HcoState,     ExtNr=ExtNr,&
                              Cat=-1,    Hier=-1,      HcoID=HcoID,&
                              AutoFill=1,Array2D=Arr2D,RC=RC        )
            IF ( RC /= HCO_SUCCESS ) RETURN
            Arr2D => NULL()
         ENDIF
  
      ENDDO !N

      ! Nullify pointers
      THISTYP   => NULL()
      Arr2D     => NULL()

      ! Leave w/ success
      CALL HCO_LEAVE ( RC )

      END SUBROUTINE HCOX_FINN_RUN
!EOC
!-----------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group     !
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCOX_FINN_INIT
!
! !DESCRIPTION: Subroutine HCOX\_FINN\_INIT initializes all module 
! arrays and variables.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE HCOX_FINN_INIT ( am_I_Root, HcoState, ExtName, &
                                  ExtState,  RC                  ) 
!
! !USE:
!
!----------------------------------------------------------------------------
! Prior to 8/11/14:
! We now get FINN emission factors from an include file (bmy, 8/11/14)
!      USE CHARPAK_MOD,            ONLY : STRSPLIT
!      USE inquireMod,             ONLY : findfreeLUN
!----------------------------------------------------------------------------
      USE HCO_STATE_MOD,          ONLY : HCO_GetHcoID
      USE HCO_STATE_MOD,          ONLY : HCO_GetExtHcoID
      USE HCO_ExtList_Mod,        ONLY : GetExtNr, GetExtOpt
!
! !INPUT/OUTPUT PARAMETERS:
!
      LOGICAL,          INTENT(IN   )  :: am_I_Root         ! root CPU?
      TYPE(HCO_State),  POINTER        :: HcoState          ! HEMCO state object 
      CHARACTER(LEN=*), INTENT(IN   )  :: ExtName           ! Extension name
      TYPE(Ext_State),  POINTER        :: ExtState          ! Extensions object
      INTEGER,          INTENT(INOUT)  :: RC                ! Return status
!
! !REVISION HISTORY:
!  02 Jan 2013 - J. Mao & J. Fisher - Initial version, based on GFED3
!  05 May 2014 - J.A. Fisher - Replace NOx emissions with NO emissions as part
!                              of removal of NOx-Ox partitioning
!  18 Jun 2014 - C. Keller   - Now a HEMCO extension.
!  11 Aug 2014 - R. Yantosca - Now get FINN emission factors and species names
!                              from include file hcox_finn_include.H.
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES
!
      INTEGER               :: N_SPEC_EMFAC       ! # of cols in EF_CO2_FILE   (w/o first two)
      INTEGER               :: N_NMOC             ! # of cols in VOC_SPEC_FILE (w/o first two)
      INTEGER               :: IU_FILE, L, N_LUMPED, tmpNr
      INTEGER               :: AS, IOS, M, N, NDUM
      INTEGER               :: N_SPECSTRS, N_NMOCSTRS
      LOGICAL               :: IS_NMOC, Matched, Missing
      CHARACTER(LEN=1023)   :: ADUM
      CHARACTER(LEN=255)    :: SDUM(255)
      CHARACTER(LEN=255)    :: IN_SPEC_NAME(255)
      CHARACTER(LEN=255)    :: IN_NMOC_NAME(255)
      CHARACTER(LEN=255)    :: TMPNAME
      CHARACTER(LEN=  6)    :: SPCNAME
      REAL*8                :: C_MOLEC
      REAL(dp), ALLOCATABLE :: EMFAC_IN(:,:)
      REAL(dp), ALLOCATABLE :: NMOC_RATIO_IN(:,:)
      REAL*8                :: NMOC_EMFAC(N_EMFAC), NMOC_RATIO(N_EMFAC)
      REAL(dp)              :: AdjFact
      CHARACTER(LEN=255)    :: MSG, EF_CO2_FILE, VOC_SPEC_FILE

      !=================================================================
      ! HCOX_FINN_INIT begins here!
      !=================================================================

      ! Extension Nr.
      ExtNr = GetExtNr( TRIM(ExtName) )
      IF ( ExtNr <= 0 ) RETURN
 
      ! Enter 
      CALL HCO_ENTER ( 'HCOX_FINN_INIT (HcoX_FINN_Mod.F90)', RC )
      IF ( RC /= HCO_SUCCESS ) RETURN

      ! ---------------------------------------------------------------------- 
      ! Get settings 
      ! ---------------------------------------------------------------------- 
 
      ! Get file with CO2 emission factor ratios, as set in configuration file,
      ! and read corresponding number of columns
      CALL GetExtOpt ( ExtNr, 'EF ratios CO2', &
                       OptValChar=EF_CO2_File, RC=RC )
      IF ( RC /= HCO_SUCCESS ) RETURN
      CALL GetExtOpt ( ExtNr, 'EF rat columns', &
                       OptValInt=N_SPEC_EMFAC, RC=RC )
      IF ( RC /= HCO_SUCCESS ) RETURN
      N_SPEC_EMFAC = N_SPEC_EMFAC - 2 ! First two columns are not used

      ! Get file with VOC speciations, as set in configuration file 
      CALL GetExtOpt ( ExtNr, 'VOC speciation', &
                       OptValChar=VOC_SPEC_FILE, RC=RC )
      IF ( RC /= HCO_SUCCESS ) RETURN
      CALL GetExtOpt ( ExtNr, 'VOC spec columns', &
                       OptValInt=N_NMOC, RC=RC )
      IF ( RC /= HCO_SUCCESS ) RETURN
      N_NMOC = N_NMOC - 2 ! First two columns are not used

      ! Use daily data?
      tmpName = TRIM(ExtName) // "_daily"
      tmpNr   = GetExtNr( TRIM(tmpName) )
      IF ( tmpNr > 0 ) THEN
         UseDay = .TRUE.
      ELSE
         UseDay = .FALSE.
      ENDIF

      ! ---------------------------------------------------------------------- 
      ! Allocate arrays
      ! ---------------------------------------------------------------------- 

      ! temporary arrays
      ALLOCATE ( EMFAC_IN(N_SPEC_EMFAC, N_EMFAC), STAT=AS )
      IF ( AS/=0 ) THEN
         CALL HCO_ERROR( 'Cannot allocate EMFAC_IN', RC )
         RETURN
      ENDIF
      ALLOCATE ( NMOC_RATIO_IN(N_NMOC, N_EMFAC), STAT=AS )
      IF ( AS/=0 ) THEN
         CALL HCO_ERROR( 'Cannot allocate NMOC_RATIO_IN', RC )
         RETURN
      ENDIF
      EMFAC_IN     = 0.0_dp
      N_SPEC_EMFAC = 0.0_dp

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

!-----------------------------------------------------------------------------
! Prior to 8/11/14
! Now initialize this data with preprocessed data file hcox_finn_include.H
!      ! ---------------------------------------------------------------------- 
!      ! Define FINN species names
!      ! ---------------------------------------------------------------------- 
!      
!      ! Species listed in emission factor ratios (CO2/X) table (except NMOC,
!      ! which is speciated as specified in the VOC speciation table).
!      FINN_SPEC_NAME(1)  = 'CO2'
!      FINN_SPEC_NAME(2)  = 'CO'
!      FINN_SPEC_NAME(3)  = 'CH4'
!      FINN_SPEC_NAME(4)  = 'NOx'
!      FINN_SPEC_NAME(5)  = 'SO2'
!      FINN_SPEC_NAME(6)  = 'OC'
!      FINN_SPEC_NAME(7)  = 'BC'
!      FINN_SPEC_NAME(8)  = 'NH3'
!      FINN_SPEC_NAME(9)  = 'NO'    ! Currently not used
!      FINN_SPEC_NAME(10) = 'NO2'   ! Currently not used
!
!      ! Species listed in VOC speciation table
!      FINN_SPEC_NAME(11) = 'ACET'
!      FINN_SPEC_NAME(12) = 'ACTA'   ! Not currently emitted by BB in GC
!      FINN_SPEC_NAME(13) = 'ALD2'
!      FINN_SPEC_NAME(14) = 'ALK4'
!      FINN_SPEC_NAME(15) = 'APINE'  ! Currently lumped into MTPA
!      FINN_SPEC_NAME(16) = 'AROM'   ! Currently not used
!      FINN_SPEC_NAME(17) = 'BENZ'
!      FINN_SPEC_NAME(18) = 'BPINE'  ! Currently lumped into MTPA
!      FINN_SPEC_NAME(19) = 'C2H2'
!      FINN_SPEC_NAME(20) = 'C2H4'
!      FINN_SPEC_NAME(21) = 'C2H6'
!      FINN_SPEC_NAME(22) = 'C3H8'
!      FINN_SPEC_NAME(23) = 'CARENE' ! Currently lumped into MTPA
!      FINN_SPEC_NAME(24) = 'CH2Br2'
!      FINN_SPEC_NAME(25) = 'CH2O'
!      FINN_SPEC_NAME(26) = 'CH3Br'
!      FINN_SPEC_NAME(27) = 'CH3CN'
!      FINN_SPEC_NAME(28) = 'CH3I'
!      FINN_SPEC_NAME(29) = 'DMS'
!      FINN_SPEC_NAME(30) = 'EOH'    ! Not currently emitted in GC
!      FINN_SPEC_NAME(31) = 'ETBENZ' ! Currently lumped with TOLU
!      FINN_SPEC_NAME(32) = 'FUR'    ! Currently not used
!      FINN_SPEC_NAME(33) = 'GLYC'
!      FINN_SPEC_NAME(34) = 'GLYX'
!      FINN_SPEC_NAME(35) = 'HAC'
!      FINN_SPEC_NAME(36) = 'HCN'    ! Not currently emitted in GC
!      FINN_SPEC_NAME(37) = 'HCOOH'  ! Not currently emitted by BB in GC
!      FINN_SPEC_NAME(38) = 'HNO2'   ! Not currently emitted in GC
!      FINN_SPEC_NAME(39) = 'ISOP'   ! Not currently emitted by BB in GC
!      FINN_SPEC_NAME(40) = 'LIMO'
!      FINN_SPEC_NAME(41) = 'MACR'   ! Not currently emitted in GC
!      FINN_SPEC_NAME(42) = 'MEK'
!      FINN_SPEC_NAME(43) = 'MGLY'
!      FINN_SPEC_NAME(44) = 'MNO3'
!      FINN_SPEC_NAME(45) = 'MOH'    ! Not currently emitted in GC
!      FINN_SPEC_NAME(46) = 'MTPO'   ! Not currently emitted in GC
!      FINN_SPEC_NAME(47) = 'MVK'    ! Not currently emitted in GC
!      FINN_SPEC_NAME(48) = 'PRPE'
!      FINN_SPEC_NAME(49) = 'R4N2'   ! Not currently emitted in GC
!      FINN_SPEC_NAME(50) = 'RCHO'   ! Not currently emitted by BB in GC
!      FINN_SPEC_NAME(51) = 'RCOOH'  ! Currently not used
!      FINN_SPEC_NAME(52) = 'ROH'    ! Currently not used
!      FINN_SPEC_NAME(53) = 'SESQ'   ! Currently not used
!      FINN_SPEC_NAME(54) = 'STYR'   ! Currently lumped with TOLU
!      FINN_SPEC_NAME(55) = 'TMB'    ! Currently lumped with XYLE
!      FINN_SPEC_NAME(56) = 'TOLU'
!      FINN_SPEC_NAME(57) = 'XYLE'
!      FINN_SPEC_NAME(58) = 'H2'     ! Currently not used
!
!      ! ---------------------------------------------------------------------- 
!      ! Read emission factors ([mole CO2]/[mole X])
!      ! ---------------------------------------------------------------------- 
!
!      ! Find a free file LUN
!      IU_FILE = findFreeLUN()
!
!      ! Open emission factor file (ASCII format)
!      OPEN( IU_FILE, FILE=TRIM(EF_CO2_File), STATUS='OLD', IOSTAT=IOS )
!      IF ( IOS /= 0 ) THEN
!         MSG = 'Error 1 reading ' // TRIM(EF_CO2_FILE)
!         CALL HCO_ERROR( MSG, RC )
!         RETURN
!      ENDIF 
!
!      ! Skip unnecessary header lines
!      DO N = 1, 2
!         READ( IU_FILE, *, IOSTAT=IOS )
!         IF ( IOS /= 0 ) THEN
!            MSG = 'Error 2 reading ' // TRIM(EF_CO2_FILE)
!            CALL HCO_ERROR( MSG, RC )
!            RETURN
!         ENDIF 
!      ENDDO
!
!      ! Read species names for emission ratio file 
!      READ( IU_FILE, '(A)', IOSTAT=IOS ) ADUM
!      CALL STRSPLIT(ADUM,',',IN_SPEC_NAME,N_SPECSTRS)
!
!      ! Read emission factors for each species and land type
!      DO N = 1, N_EMFAC
!!         READ( IU_FILE, *, IOSTAT=IOS ) NDUM, ADUM, EMFAC_IN(:,N)
!         READ( IU_FILE, '(A)', IOSTAT=IOS ) ADUM
!         IF ( IOS /= 0 ) THEN
!            MSG = 'Error 3 reading ' // TRIM(EF_CO2_FILE)
!            CALL HCO_ERROR( MSG, RC )
!            RETURN
!         ENDIF
!         CALL STRSPLIT(ADUM,',',SDUM,NDUM)
!         ! PASS TO EMFAC_IN
!
!         DO M = 1, (NDUM-2)
!            READ( SDUM(M+2), * ) EMFAC_IN(M,N)
!         ENDDO
!      ENDDO
!
!      ! Close file
!      CLOSE( IU_FILE )      
!
!      !----------------------------------------------------------------------- 
!      ! Read NMOC factors ([mole X]/[kg NMOC])
!      !----------------------------------------------------------------------- 
!
!      ! Find a free file LUN
!      IU_FILE = findFreeLUN()
!
!      ! Open emission factor file (ASCII format)
!      OPEN( IU_FILE, FILE=TRIM(VOC_SPEC_File), STATUS='OLD', IOSTAT=IOS )
!      IF ( IOS /= 0 ) THEN
!         MSG = 'Error 4 reading ' // TRIM(VOC_SPEC_FILE)
!         CALL HCO_ERROR( MSG, RC )
!         RETURN
!      ENDIF 
!
!      ! Skip unnecessary header lines
!      DO N = 1, 2
!         READ( IU_FILE, *, IOSTAT=IOS )
!         IF ( IOS /= 0 ) THEN
!            MSG = 'Error 5 reading ' // TRIM(VOC_SPEC_FILE)
!            CALL HCO_ERROR( MSG, RC )
!            RETURN
!         ENDIF 
!      ENDDO
!
!      ! Read species names for emission ratio file 
!      READ( IU_FILE, '(A)', IOSTAT=IOS ) ADUM
!      CALL STRSPLIT(ADUM,',',IN_NMOC_NAME,N_NMOCSTRS)
!
!      ! Read emission factors for each species and land type
!      DO N = 1, N_EMFAC
!         READ( IU_FILE, '(A)', IOSTAT=IOS ) ADUM
!         IF ( IOS /= 0 ) THEN
!            MSG = 'Error 6 reading ' // TRIM(VOC_SPEC_FILE)
!            CALL HCO_ERROR( MSG, RC )
!            RETURN
!         ENDIF 
!         CALL STRSPLIT(ADUM,',',SDUM,NDUM)
!         ! PASS TO NMOC_RATIO_IN
!         DO M = 1, (NDUM-2)
!            READ( SDUM(M+2), * ) NMOC_RATIO_IN(M,N)
!         ENDDO
!      ENDDO
!
!      ! Close file
!      CLOSE( IU_FILE )      
!-----------------------------------------------------------------------------

      !-----------------------------------------------------------------------
      ! FINN now initializes the species names and emissions factors with
      ! hard-coded F90 assignment statements in the include file 
      ! "hcox_finn_include.H".  This include file takes the species names
      ! and emission factors from the input files "FINN_EFratios.csv" and 
      ! "FINN_VOC_speciation.csv".  If new emission factors are issued, then
      ! you can regenerate this include file simply by running the Perl
      ! script HEMCO/Extensions/Preprocess/finn.pl. (bmy, 8/11/14)
      !-----------------------------------------------------------------------
#include "hcox_finn_include.H"

      ! ---------------------------------------------------------------------- 
      ! Match specified species with FINN species. The species to be used are 
      ! specified in the HEMCO configuration file.
      ! ---------------------------------------------------------------------- 

      ! Write to log file
      MSG = 'Use FINN extension'
      CALL HCO_MSG( MSG, SEP1='-' )
      WRITE(MSG,*) '   - CO2 EF scale factors    : ', TRIM(EF_CO2_FILE) 
      CALL HCO_MSG( MSG )
      WRITE(MSG,*) '   - VOC speciations         : ', TRIM(VOC_SPEC_FILE)
      CALL HCO_MSG( MSG )
      WRITE(MSG,*) '   - Use daily data          : ', UseDay
      CALL HCO_MSG( MSG )

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

         ! For model species NO, the emission factors are taken from FINN species NOx
         ! For model species MTPA, the emission factors are taken from FINN species APINE (BPINE and CARENE will be lumped into it as well) 
         SELECT CASE ( TRIM(SpcName) )
            CASE ( 'NO' )
               SpcName = 'NOx'
            CASE ('MTPA' )
               SpcName = 'APINE'
         END SELECT

         ! For lumped species, we have to repeat the lookup multiple times,
         ! so use a while loop here.
         ! For example, for species TOLU this will make sure that FINN species
         ! 'TOLU', 'ETBENZ', and 'STYR' are associated with HEMCO species TOLU.
         DO WHILE ( Missing )

            ! Search for SpcName in FINN
            DO N = 1, N_SPEC 
               IF ( TRIM(SpcName) == TRIM(FINN_SPEC_NAME(N)) ) THEN
                  FinnIDs(N) = L
                  Matched    = .TRUE.
   
                  MSG = '   - FINN species ' // TRIM(FINN_SPEC_NAME(N)) // &
                        '     will be emitted as ' // TRIM(SpcNames(L))
                  CALL HCO_MSG( MSG )
   
                  ! Reset variables
                  IS_NMOC    = .FALSE.
                  C_MOLEC    = 1d0
                  NMOC_RATIO = 0d0

                  ! Get emission factor in [kg X]/[kg CO2]. 
                  DO M = 1, N_SPECSTRS
                     TMPNAME = IN_SPEC_NAME(M)
                     IF ( TRIM(FINN_SPEC_NAME(N)) == TRIM(TMPNAME(5:8)) ) THEN
                        ! First two entries are not species. Also, EMFAC is
                        ! stored as [mole CO2]/[mole X], but we want the inverse.
                        ! This gives us [mole X]/[mole CO2]. To convert this to
                        ! [kg X]/[kg CO2], we also need to adjust for the molecular
                        ! weights of species X and CO2.
                        ! The EF ratios of OC and BC are in [mole CO2]/[g X], so 
                        ! the adjustment factor is calculated slightly differently
                        ! for those two species!
                        IF ( TRIM(FINN_SPEC_NAME(N)) == 'OC' .OR. &
                             TRIM(FINN_SPEC_NAME(N)) == 'BC'       ) THEN
                           AdjFact = 1.0_dp / MW_CO2
                        ELSE
                           AdjFact = 1.0_dp / MW_CO2 * &
                                     HcoState%Spc(HcoIDs(L))%MW_g
                        ENDIF
                        FINN_EMFAC(N,:) = AdjFact / EMFAC_IN(M-2,:)
                        WRITE( MSG, 200 ) TRIM( FINN_SPEC_NAME(N)) 
                        CALL HCO_MSG( MSG )
                        EXIT
                     ! NMOC_EMFAC is converted to [kg NMOC]/[kg CO2].
                     ! Input unit is [mole CO2]/[mole NMOC].
                     ELSE IF ( TRIM(TMPNAME(5:8)) == 'NMOC' ) THEN
                        AdjFact = MW_NMOC / MW_CO2
                        NMOC_EMFAC = AdjFact / EMFAC_IN(M-2,:)
                     ENDIF
                  ENDDO
 200              FORMAT( 'Found FINN emission ratio for species ',a5 )
   
                  DO M = 1, N_NMOCSTRS
                     TMPNAME = IN_NMOC_NAME(M)
                     IF ( TRIM(FINN_SPEC_NAME(N)) == TRIM(TMPNAME) ) THEN
                        ! First two entries are not species
                        NMOC_RATIO = NMOC_RATIO_IN(M-2,:)
                        IS_NMOC = .TRUE.
                        WRITE( MSG, 201 ) TRIM( FINN_SPEC_NAME(N) )
                        CALL HCO_MSG( MSG )
                        EXIT
                     ENDIF
                  ENDDO
 201              FORMAT( 'Found FINN NMOC factor for species ',a5 )
   
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
                       AdjFact         = HcoState%Spc(HcoIDs(L))%MW_g
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
      CALL HCO_LEAVE ( RC ) 
 
      END SUBROUTINE HCOX_FINN_INIT
!EOC
!-----------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group     !
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCOX_FINN_FINAL
!
! !DESCRIPTION: Subroutine HCOX\_FINN\_FINAL deallocates all module
! arrays.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE HCOX_FINN_FINAL
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

      END SUBROUTINE HCOX_FINN_FINAL
!EOC
      END MODULE HCOX_FINN_MOD
!EOM
