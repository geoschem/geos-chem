!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Model                            !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: GEOS_Analysis 
!
! !DESCRIPTION: Module to apply analysis fields to GEOS-Chem tracers. 
!\\
!\\
! !INTERFACE:
!
MODULE GEOS_Analysis 
!
! !USES:
!
  ! GEOS-Chem
  USE Error_Mod
  USE ErrCode_Mod
  USE Precision_Mod
  USE HCO_Error_Mod

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC   :: GEOS_AnaInit 
  PUBLIC   :: GEOS_AnaRun
  PUBLIC   :: GEOS_AnaFinal
!
! !PRIVATE MEMBER FUNCTIONS:
!
  PRIVATE  :: ReadSettings_
  PRIVATE  :: DoAnalysis_
  PRIVATE  :: GetAnaTime_
  PRIVATE  :: ReplaceChar_ 
!
! !PRIVATE TYPES:
!
  ! Number of species with analysis on
  INTEGER                         :: nAnaSpec 

  ! Options for dependent species 
  TYPE Spec2Opt
     CHARACTER(LEN=63)            :: Spec2Name
     LOGICAL                      :: Spec2Strat 
     LOGICAL                      :: Spec2Trop 
     REAL                         :: Spec2MinRatio
     REAL                         :: Spec2MaxRatio
  END TYPE

  ! Analysis options object. A separate object will be created for each analysed species/family
  TYPE AnaOptions
     CHARACTER(LEN=63)            :: SpecName
     LOGICAL                      :: Active
     LOGICAL                      :: DryRun
     INTEGER                      :: AnalysisFreq
     INTEGER                      :: AnalysisHour
     INTEGER                      :: AnalysisMinute
     LOGICAL                      :: ForwardLooking
     LOGICAL                      :: SkipPredictor
     CHARACTER(LEN=127)           :: FldNameHco
     CHARACTER(LEN=63)            :: FileVarUnit
     INTEGER                      :: DryFlag 
     LOGICAL                      :: IsIncrement 
     INTEGER                      :: IAU
     INTEGER                      :: AnalysisWindow
     LOGICAL                      :: HasMask
     CHARACTER(LEN=127)           :: MskNameHco
     REAL                         :: MaskThreshold 
     LOGICAL                      :: InStrat
     LOGICAL                      :: InTrop
     INTEGER                      :: AnaL1 
     INTEGER                      :: AnaL2 
     INTEGER                      :: AnaL3 
     INTEGER                      :: AnaL4
     REAL                         :: AnaFraction
     INTEGER                      :: StratSponge 
     REAL                         :: MaxChangeStrat
     REAL                         :: MaxChangeTrop
     REAL                         :: MaxRatioStrat
     REAL                         :: MaxRatioTrop
     REAL                         :: MinRatioStrat
     REAL                         :: MinRatioTrop
     REAL                         :: MinConc
     INTEGER                      :: nSpec2
     TYPE(Spec2Opt), POINTER      :: Spec2(:) => NULL()
     INTEGER                      :: ErrorMode
     INTEGER                      :: Verbose
  END TYPE AnaOptions

  ! List holding all analysis information
  TYPE(AnaOptions), POINTER  :: AnaConfig(:) => NULL()

  ! Main configuration file for analysis options
  CHARACTER(LEN=127), PARAMETER :: AnaConfigFile = './geoschem_analysis.yml'

!
! !REVISION HISTORY:
!  25 May 2022 - C. Keller - initial version (refactored Chem_GridCompMod)
!  See https://github.com/geoschem/geos-chem for history
!EOP
!------------------------------------------------------------------------------
!BOC
CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Model                            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: GEOS_AnaInit 
!
! !DESCRIPTION: Initialize the GEOS analysis module by reading all analysis 
!  settings from configuration files.
!
! !INTERFACE:
!
  SUBROUTINE GEOS_AnaInit( am_I_Root, AnaPhase, RC )
!
! !USES:
!
  USE QfYaml_Mod
!
! !INPUT PARAMETERS:
!
    LOGICAL,              INTENT(IN)    :: am_I_Root ! Root PET?
!
! !OUTPUT PARAMETERS:
!
    INTEGER, INTENT(OUT)                :: AnaPhase  ! Do analysis after run phase 1 or 2?
    INTEGER, INTENT(INOUT)              :: RC        ! Success or failure
!
! !REMARKS:
!
! !REVISION HISTORY:
!  25 May 2022 - C. Keller - initial version (refactored Chem_GridCompMod)
!  See https://github.com/geoschem/geos-chem for history
!EOP
!------------------------------------------------------------------------------
!BOC
! 
! !LOCAL VARIABLES:
!
    CHARACTER(LEN=255)            :: Iam
    CHARACTER(LEN=511)            :: errMsg
    LOGICAL                       :: FileExists
    TYPE(QFYAML_t)                :: Config, ConfigAnchored
    CHARACTER(LEN=QFYAML_NamLen)  :: key
    INTEGER                       :: N, v_int, ix
    CHARACTER(LEN=63)             :: iSpecName

    !=======================================================================
    ! GEOS_AnaInit begins here 
    !=======================================================================

    ! callback name
    Iam = 'GEOS-Chem::GEOS_AnaInit'

    ! Initialize
    RC       = GC_SUCCESS
    ANAPHASE = 2
    nAnaSpec = 0

    ! Check if file exists
    INQUIRE( FILE=AnaConfigFile, EXIST=FileExists )
    IF ( .NOT. FileExists ) THEN
       IF ( am_I_Root ) THEN
          WRITE( 6, * ) TRIM(Iam)//": analysis configuration file not there: "//TRIM(AnaConfigFile)
          WRITE( 6, * ) " --> no GEOS-Chem species will be nudged"
       ENDIF
       RETURN
    ENDIF

    ! Read YAML file into Config object
    CALL QFYAML_Init( AnaConfigFile, Config, ConfigAnchored, RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error reading analysis configuration file: ' // TRIM( AnaConfigFile )
       CALL GC_Error( errMsg, RC, Iam )
       RETURN
    ENDIF

    ! Run phase after which to apply analysis
    key = "general%runphase"
    v_int = 2
    CALL QFYAML_Add_Get( Config, TRIM( key ), v_int, "", RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error parsing ' // TRIM( key ) // '!'
       CALL GC_Error( errMsg, RC, Iam ) 
       RETURN
    ENDIF
    ANAPHASE = v_int

    ! Get list of species 
    key = "general%nspecies"
    v_int = 0 
    CALL QFYAML_Add_Get( Config, TRIM( key ), v_int, "", RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error parsing ' // TRIM( key ) // '!'
       CALL GC_Error( errMsg, RC, Iam ) 
       RETURN
    ENDIF
    nAnaSpec = v_int 

    ! Verbose
    IF ( am_I_Root ) THEN
       WRITE(6,*) 'Number of analysis species: ',nAnaSpec
       WRITE(6,*) 'Analysis phase set to ',ANAPHASE
    ENDIF

    ! Read settings for all species from config file
    IF ( nAnaSpec > 0 ) THEN
       ALLOCATE( AnaConfig(nAnaSpec), STAT=RC )
       IF ( RC /= 0 ) THEN
          errMsg = 'AnaConfig could not be allocated!'
          CALL GC_Error( errMsg, RC, Iam ) 
          RETURN
       ENDIF

       ! Now read settings for each species
       DO N=1,nAnaSpec
          CALL ReadSettings_( am_I_Root, Config, N, RC=RC )
          IF ( RC /= GC_SUCCESS ) EXIT 
       ENDDO
    ENDIF

    ! Cleanup
    CALL QFYAML_CleanUp( Config         )
    CALL QFYAML_CleanUp( ConfigAnchored )

  END SUBROUTINE GEOS_AnaInit 
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Model                            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: GEOS_AnaRun  
!
! !DESCRIPTION: Driver routine to run the GEOS analysis module. 
!
! !INTERFACE:
!
  SUBROUTINE GEOS_AnaRun( Input_Opt, State_Met, State_Chm, State_Grid, State_Diag, RC )
!
! !USES:
!
  USE Input_Opt_Mod,         ONLY : OptInput
  USE State_Chm_Mod,         ONLY : ChmState         ! Chemistry State obj
  USE State_Met_Mod,         ONLY : MetState         ! Meteorology State obj
  USE State_Grid_Mod,        ONLY : GrdState         ! Meteorology State obj
  USE State_Diag_Mod,        ONLY : DgnState         ! Meteorology State obj
  USE UnitConv_Mod,          ONLY : Convert_Spc_Units
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(OptInput)                             :: Input_Opt
    TYPE(MetState)                             :: State_Met
    TYPE(ChmState)                             :: State_Chm
    TYPE(GrdState)                             :: State_Grid
    TYPE(DgnState)                             :: State_Diag
!
! !OUTPUT PARAMETERS:
!
    INTEGER, INTENT(INOUT)                     :: RC        ! Success or failure
!
! !REMARKS:
!
! !REVISION HISTORY:
!  25 May 2022 - C. Keller - initial version (refactored Chem_GridCompMod)
!  See https://github.com/geoschem/geos-chem for history
!EOP
!------------------------------------------------------------------------------
!BOC
! 
! !LOCAL VARIABLES:
!
    INTEGER                       :: ispec
    CHARACTER(LEN=63)             :: OrigUnit
    CHARACTER(LEN=255)            :: ErrMsg, ThisLoc 

    !=======================================================================
    ! GEOS_AnaRun begins here 
    !=======================================================================

    ! callback name
    RC      = GC_SUCCESS
    ThisLoc = ' -> GEOS-Chem AnaRun'

    ! Convert to total mixing ratio
    CALL Convert_Spc_Units ( Input_Opt, State_Chm, State_Grid, State_Met, &
                             'kg/kg total', RC, OrigUnit=OrigUnit )
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Unit conversion error (start)!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN              
    ENDIF  

    ! Reset diagnostics
    IF ( State_Diag%Archive_AnaInc     ) State_Diag%AnaInc(:,:,:,:)     = 0.0
    IF ( State_Diag%Archive_AnaIncFrac ) State_Diag%AnaIncFrac(:,:,:,:) = 0.0
    IF ( State_Diag%Archive_AnaMask    ) State_Diag%AnaMask(:,:,:,:)    = 0.0
    IF ( State_Diag%Archive_AnaMaskSum ) State_Diag%AnaMaskSum(:,:,:)   = 0.0

    ! Do analysis for all analysis species
    IF ( nAnaSpec > 0 .AND. ASSOCIATED(AnaConfig) ) THEN 
       DO ispec=1,nAnaSpec
          IF ( AnaConfig(ispec)%Active ) THEN
             CALL DoAnalysis_( ispec, Input_Opt, State_Met, State_Chm, &
                               State_Grid, State_Diag, RC )
             IF ( RC /= GC_SUCCESS ) THEN
                ErrMsg = 'DoAnalysis_ error!'
                CALL GC_Error( ErrMsg, RC, ThisLoc )
                RETURN              
             ENDIF  
          ENDIF
       ENDDO
    ENDIF

    ! Convert species back to original unit 
    CALL Convert_Spc_Units ( Input_Opt, State_Chm, State_Grid, State_Met, &
                             OrigUnit, RC )
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Unit conversion error (end)!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN              
    ENDIF  

    ! Successful return
    RC = GC_SUCCESS 

  END SUBROUTINE GEOS_AnaRun  
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Model                            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: GEOS_AnaFinal
!
! !DESCRIPTION: Finalize the GEOS analysis module
!
! !INTERFACE:
!
  SUBROUTINE GEOS_AnaFinal( RC )
!
! !OUTPUT PARAMETERS:
!
    INTEGER, INTENT(OUT)                :: RC        ! Success or failure
!
! !REMARKS:
!
! !REVISION HISTORY:
!  25 May 2022 - C. Keller - initial version (refactored Chem_GridCompMod)
!  See https://github.com/geoschem/geos-chem for history
!EOP
!------------------------------------------------------------------------------
!BOC
    CHARACTER(LEN=255)            :: Iam
    INTEGER                       :: N

    !=======================================================================
    ! GEOS_AnaInit begins here 
    !=======================================================================

    ! callback name
    Iam = 'GEOS_AnaFinal'
    RC = GC_SUCCESS

    ! Clean up 
    IF ( ASSOCIATED(AnaConfig) ) THEN
       DO N=1,nAnaSpec
          IF ( ASSOCIATED(AnaConfig(N)%Spec2) ) DEALLOCATE(AnaConfig(N)%Spec2)
       ENDDO
       DEALLOCATE( AnaConfig )
    ENDIF
    AnaConfig => NULL()
    nAnaSpec = 0

  END SUBROUTINE GEOS_AnaFinal
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Model                            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: DoAnalysis_ 
!
! !DESCRIPTION: Routine to perform the analysis for a given species. 
!
! !INTERFACE:
!
  SUBROUTINE DoAnalysis_( ispec, Input_Opt, State_Met, State_Chm, State_Grid, State_Diag, RC )
!
! !USES:
!
  USE Input_Opt_Mod,         ONLY : OptInput
  USE State_Chm_Mod,         ONLY : ChmState         ! Chemistry State obj
  USE State_Met_Mod,         ONLY : MetState         ! Meteorology State obj
  USE State_Grid_Mod,        ONLY : GrdState         ! Grid State obj
  USE State_Diag_Mod,        ONLY : DgnState, DgnMap ! Diagnostics State obj
  USE HCO_Utilities_GC_Mod,  ONLY : HCO_GC_EvalFld
  USE TIME_MOD,              ONLY : GET_TS_CHEM
  Use PhysConstants,         ONLY : AIRMW
  USE HCOIO_Util_Mod,        ONLY : HCOIO_IsValid 
  USE HCO_State_GC_Mod,      ONLY : HcoState
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,             INTENT(IN)            :: ispec     ! analysis species index
    TYPE(OptInput)                             :: Input_Opt
    TYPE(MetState)                             :: State_Met
    TYPE(ChmState)                             :: State_Chm
    TYPE(GrdState)                             :: State_Grid
    TYPE(DgnState)                             :: State_Diag
!
! !OUTPUT PARAMETERS:
!
    INTEGER, INTENT(OUT)                       :: RC        ! Success or failure
!
! !REMARKS:
!
! !REVISION HISTORY:
!  25 May 2022 - C. Keller - initial version (refactored Chem_GridCompMod)
!  See https://github.com/geoschem/geos-chem for history
!EOP
!------------------------------------------------------------------------------
!BOC
! 
! !LOCAL VARIABLES:
!
    LOGICAL                    :: am_I_Root
    TYPE(AnaOptions), POINTER  :: iopt => NULL() 
    LOGICAL                    :: TimeForAna, HasField
    INTEGER                    :: yy, mm, dd, h, m, s
    INTEGER                    :: VarID
    INTEGER                    :: StratCount
    REAL                       :: ThisHour
    CHARACTER(LEN=255)         :: SpecName, Spec2Name, FldName
    REAL, ALLOCATABLE          :: AnaMask2d(:,:),  AnaMask3d(:,:,:)
    REAL, ALLOCATABLE          :: SpcBkg(:,:,:), SpcAsm(:,:,:)
    REAL, ALLOCATABLE          :: Spc2Bkg(:,:,:,:), Spc2Asm(:,:,:,:)
    INTEGER                    :: I, J, L, N, IM, JM, LM, indSpc
    INTEGER, ALLOCATABLE       :: indSpc2(:)
    INTEGER                    :: UnitFlag, DryFlag, NNEG
    REAL                       :: OldRatio, NewRatio 
    REAL                       :: wgt, tropwgt, stratwgt
    REAL                       :: DilFact, tsChem
    REAL                       :: frac, diff, maxChange, maxRatio, minRatio
    REAL                       :: mwSpc
    REAL, ALLOCATABLE          :: mwSpc2(:)
    REAL                       :: SpcAna, SpcNew
    REAL                       :: MinConc
    LOGICAL                    :: UpdateSpec2
    CHARACTER(LEN=255)         :: Iam
    CHARACTER(LEN=255)         :: ErrMsg, ThisLoc 

    ! Diagnostics
    INTEGER                    :: DgnId

    ! temporary arrays
    REAL, ALLOCATABLE          :: Qtmp(:,:,:)

    REAL(hp), POINTER          :: AnaPtr(:,:,:), MskPtr(:,:,:)

    !=======================================================================
    ! DoAnalysis_ begins here 
    !=======================================================================

    ! Callback
    Iam     = 'GEOS_Analysis::DoAnalysis_'
    ThisLoc = 'GEOS_Analysis::DoAnalysis_'
    RC      = GC_SUCCESS

    ! Root CPU?
    am_I_Root = Input_Opt%amIRoot 

    ! Get settings
    iopt => AnaConfig(ispec)
    SpecName = iopt%SpecName
    MinConc  = iopt%MinConc

    ! Check if it's time to do the analysis
    TimeForAna = .FALSE. 
    CALL GetAnaTime_( iopt%ForwardLooking, yy, mm, dd, h, m, s, RC )
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Error in GetAnaTime_'
       CALL GC_Error( ErrMsg, RC, ThisLoc)
       RETURN
    ENDIF

    ThisHour = real(h)
    DilFact  = 1.0

    ! Always do analysis if spreading increment evenly
    IF ( iopt%IAU ) THEN
       TimeForAna = .TRUE.
       ! Calculate dilution factor, to be applied to analysis/increment weight
       tsChem  = GET_TS_CHEM()
       DilFact = real(iopt%AnalysisWindow)*(3600./tsChem)
    ! Otherwise, use specified analysis frequency and hour/minute offsets
    ELSE 
       IF ( m==iopt%AnalysisMinute .AND. MOD(h,iopt%AnalysisFreq)==iopt%AnalysisHour ) TimeForAna = .TRUE. 
    ENDIF

    ! Eventually skip during predictor step
    IF ( iopt%SkipPredictor .AND. State_Grid%PredictorIsActive ) THEN
       TimeForAna = .FALSE.
    ENDIF

    ! Check if file exists (only if it's time to do the analysis)
    HasField = .FALSE.
    AnaPtr => NULL()
    IF ( TimeForAna ) THEN

       ! HEMCO field name 
       FldName = iopt%FldNameHco

       ! Check first if the field exists / is valid (in ESMF, fields can become invalid)
       CALL HCOIO_IsValid( HcoState, FldName, HasField, RC )
       IF ( RC /= GC_SUCCESS ) HasField = .FALSE. 

       ! Evaluate field if it is valid to do so
       IF ( HasField ) THEN

          ALLOCATE(AnaPtr(State_Grid%NX,State_Grid%NY,State_Grid%NZ))
          CALL HCO_GC_EvalFld ( Input_Opt, State_Grid, FldName, AnaPtr, RC )
          IF ( RC /= GC_SUCCESS ) HasField = .FALSE. 

          ! Eventually read mask field
          IF ( iopt%HasMask ) THEN
             FldName = iopt%MskNameHco
             ALLOCATE(MskPtr(State_Grid%NX,State_Grid%NY,State_Grid%NZ))
             CALL HCO_GC_EvalFld ( Input_Opt, State_Grid, FldName, MskPtr, RC )
             IF ( RC /= GC_SUCCESS ) THEN
                ErrMsg = 'Error getting mask: '//TRIM(FldName)
                CALL GC_Error( ErrMsg, RC, ThisLoc)
                RETURN
             ENDIF
          ENDIF

       ! Handle cases where field is not found. If error mode is set to 0, stop with error.
       ! Print a warning otherwise.
       ELSE
          IF ( iopt%ErrorMode == 0 ) THEN 
             ErrMsg = 'Field not found: '//TRIM(iopt%FldNameHco)//'; to get past this error, set error mode to > 0'
             CALL GC_Error( ErrMsg, RC, ThisLoc)
             RETURN
          ELSE
             IF ( am_I_Root .AND. (iopt%Verbose>1) ) THEN
                WRITE(*,*) 'GEOS-Chem analysis: field not found: '//TRIM(iopt%FldNameHco)//'; will skip analysis'
                !CALL GC_Warning( ErrMsg, RC, ThisLoc=ThisLoc)
             ENDIF
          ENDIF
       ENDIF
    ENDIF

    ! Apply increments if it's time to do so and if file exists
    ! ---------------------------------------------------------
    IF ( HasField ) THEN 

       ! Verbose
       IF ( am_I_Root .AND. (iopt%Verbose>0) ) THEN
          IF ( .NOT. iopt%DryRun ) THEN
             WRITE(*,100) SpecName,yy,mm,dd,h,m
100          FORMAT( "GEOS-Chem: apply analysis for species ",A5," for ",I4.4,"-",I2.2,"-",I2.2," ",I2.2,":",I2.2)
          ELSE
             WRITE(*,110) SpecName,yy,mm,dd,h,m
110          FORMAT( "GEOS-Chem: do analysis dry-run for species ",A5," for ",I4.4,"-",I2.2,"-",I2.2," ",I2.2,":",I2.2)
          ENDIF
       ENDIF

       ! Select GEOS-Chem index and molecular weight for analysis species. Also get the same for 2nd species (if used) 
       indSpc   = -1
       IF ( iopt%nSpec2 > 0 ) THEN
          ALLOCATE(indSpc2(iopt%nSpec2))
          ALLOCATE(mwSpc2(iopt%nSpec2))
       ENDIF 
       DO I = 1, State_Chm%nSpecies
          IF ( TRIM(State_Chm%SpcData(I)%Info%Name) == TRIM(SpecName) ) THEN
             indSpc = I
             mwSpc  = State_Chm%SpcData(I)%Info%MW_g
          ENDIF
          IF ( iopt%nSpec2 > 0 ) THEN
             DO N=1,iopt%nSpec2
                IF ( TRIM(State_Chm%SpcData(I)%Info%Name) == TRIM(iopt%Spec2(N)%Spec2Name) ) THEN
                   indSpc2(N) = I
                   mwSpc2(N)  = State_Chm%SpcData(I)%Info%MW_g
                ENDIF
             ENDDO
          ENDIF
       ENDDO
       ! Make sure species IDs and MWs are valid
       IF ( indSpc <= 0 .OR. mwSpc <= 0.0 ) THEN
          ErrMsg = 'Invalid species index/MW: '//TRIM(SpecName)
          CALL GC_Error( ErrMsg, RC, ThisLoc)
          RETURN
       ENDIF
       IF ( iopt%nSpec2 > 0 ) THEN
          DO N=1,iopt%nSpec2
             IF ( indSpc2(N) <= 0 .OR. mwSpc2(N) <= 0.0 ) THEN
                ErrMsg = 'Invalid species index/MW: '//TRIM(iopt%Spec2(N)%Spec2Name)
                CALL GC_Error( ErrMsg, RC, ThisLoc)
                RETURN
             ENDIF
          ENDDO
       ENDIF

       ! array dimensions 
       IM = SIZE(AnaPtr,1)
       JM = SIZE(AnaPtr,2)
       LM = SIZE(AnaPtr,3)

       ! Set dry flag
       DryFlag = iopt%DryFlag

       ! Set unit flag. This is to prevent parsing the unit string within the loop below
       SELECT CASE ( TRIM(iopt%FileVarUnit) )
          CASE ( 'kg/kg' )
             UnitFlag = 1
             IF ( DryFlag < 0 ) DryFlag = 0 ! assume kg/kg is total, not dry
          CASE ( 'mol/mol', 'v/v' )
             UnitFlag = 2
             IF ( DryFlag < 0 ) DryFlag = 1 ! assume is dry 
          CASE ( 'ppmv', 'ppm', 'PPMV', 'PPM' )
             UnitFlag = 3
             IF ( DryFlag < 0 ) DryFlag = 1 ! assume dry 
          CASE ( 'ppbv', 'ppb', 'PPBV', 'PPB' )
             UnitFlag = 4
             IF ( DryFlag < 0 ) DryFlag = 1 ! assume dry 
          CASE DEFAULT
             UnitFlag = 1
             IF ( DryFlag < 0 ) DryFlag = 0 ! assume kg/kg is total, not dry
       END SELECT

       ! State_Chm%Species are in kg/kg total. Make local copy in v/v dry before applying increments.
       ALLOCATE(SpcBkg(IM,JM,LM),SpcAsm(IM,JM,LM))
       SpcBkg(:,:,:) = State_Chm%Species(indSpc)%Conc / (1.-State_Met%SPHU/1000.) * AIRMW / mwSpc
       SpcAsm(:,:,:) = SpcBkg(:,:,:)
       IF ( iopt%nSpec2 > 0 ) THEN
          ALLOCATE(Spc2Bkg(IM,JM,LM,iopt%nSpec2),Spc2Asm(IM,JM,LM,iopt%nSpec2))
          DO N=1,iopt%nSpec2
             Spc2Bkg(:,:,:,N) = State_Chm%Species(indSpc2(N))%Conc / (1.-State_Met%SPHU/1000.) * AIRMW / mwSpc2(N)
             Spc2Asm(:,:,:,N) = Spc2Bkg(:,:,:,N)
          ENDDO
       ENDIF

       ! Eventually initialize diagnostic fields for analysis masks
       IF ( State_Diag%Archive_AnaMask ) THEN
          ALLOCATE(AnaMask3d(IM,JM,LM)) 
          AnaMask3d = 0.0
       ENDIF
       IF ( State_Diag%Archive_AnaMaskSum ) THEN
          ALLOCATE(AnaMask2d(IM,JM)) 
          AnaMask2d = 0.0
       ENDIF
 
       ! Number of negative cells
       NNEG = 0

       DO J=1,JM
       DO I=1,IM
          ! Loop over vertical
          StratCount = 0
          DO L=1,LM
             ! Fraction of cell in troposphere / stratosphere 
             tropwgt  = MAX(0.0,                                           &
                        MIN(1.0,                                           &
                         (State_Met%PEDGE(I,J,L)-State_Met%TROPP(I,J)) /   &
                         (State_Met%PEDGE(I,J,L)-State_Met%PEDGE(I,J,L+1)) &
                        )                                                  &
                        )
             stratwgt = 1.0 - tropwgt    

             ! Count number of cells since vertical loop start that have been (at least partly) in stratosphere
             IF ( stratwgt > 0.1 ) StratCount = StratCount + 1 

             ! Skip cell if concentration change is too small
             IF ( iopt%IsIncrement .AND. ABS(AnaPtr(I,J,L)) < MinConc ) CYCLE 

             ! Skip cell if masked out 
             IF ( iopt%HasMask ) THEN
                IF ( MskPtr(I,J,L) <= iopt%MaskThreshold ) CYCLE
             ENDIF

             ! Default weight to be given to analysis.
             wgt = iopt%AnaFraction

             ! Adjust weight based on stratosphere / troposphere flag 
             IF ( .NOT. iopt%InStrat ) wgt = wgt * tropwgt
             IF ( .NOT. iopt%InTrop  ) wgt = wgt * stratwgt
            
             ! Adjust weight based on the specified analysis levels, with gradual transition from L1 to L2 and L3 to L4
             IF ( L < iopt%AnaL2 ) wgt = wgt * ( (L-iopt%AnaL1) / (iopt%AnaL2-iopt%AnaL1) )
             IF ( L > iopt%AnaL3 ) wgt = wgt * ( (iopt%AnaL4-L) / (iopt%AnaL4-iopt%AnaL4) )

             ! Check for tropopause sponge layer when applying increments in strat
             IF ( iopt%InStrat .AND. .NOT. iopt%InTrop .AND. iopt%StratSponge > 0 ) THEN
                IF ( stratwgt > 0.0 .AND. StratCount <= iopt%StratSponge ) wgt = 0.0 
             ENDIF

             ! Adjust weight by # of time steps if spreading evenly using IAU.
             IF ( iopt%IAU ) wgt = wgt / DilFact

             ! Fraction must be between 0 and 1
             wgt = max(0.0,min(1.0,wgt))
             IF ( wgt == 0.0 ) CYCLE

             ! Get target concentration in v/v dry
             SpcAna = AnaPtr(I,J,L)
             IF ( UnitFlag == 1 ) SpcAna = SpcAna * ( AIRMW / mwSpc )
             IF ( UnitFlag == 3 ) SpcAna = SpcAna * 1.0e-6
             IF ( UnitFlag == 4 ) SpcAna = SpcAna * 1.0e-9
             IF ( DryFlag  == 0 ) SpcAna = SpcAna / ( 1. - State_Met%SPHU(I,J,L)/1000.0 )

             ! Update field
             SpcNew = SpcBkg(I,J,L)
             IF ( iopt%IsIncrement ) THEN
                SpcNew = SpcBkg(I,J,L) + wgt*SpcAna
                IF ( SpcNew <= MinConc ) THEN
                   SpcNew = MinConc
                   NNEG   = NNEG + 1
                ENDIF
             ELSE
                IF ( SpcAna >= MinConc ) THEN
                   SpcNew = wgt*SpcAna + (1.0-wgt)*SpcBkg(I,J,L)
                ELSE
                   SpcNew = MinConc
                   NNEG   = NNEG + 1
                ENDIF 
             ENDIF 

             ! Check for absolute change limit 
             IF ( stratwgt >= 0.5 ) THEN
                maxChange = iopt%MaxChangeStrat
                maxRatio  = iopt%MaxRatioStrat
                minRatio  = iopt%MinRatioStrat
             ELSE 
                maxChange = iopt%MaxChangeTrop
                maxRatio  = iopt%MaxRatioTrop 
                minRatio  = iopt%MinRatioTrop 
             ENDIF
             IF ( maxChange >= 0.0 ) THEN
                diff = SpcNew - SpcBkg(I,J,L)
                IF ( ABS(diff) > maxChange ) THEN
                   IF ( diff > 0.0 ) SpcNew = SpcBkg(I,J,L) + maxChange
                   IF ( diff < 0.0 ) SpcNew = SpcBkg(I,J,L) - maxChange
                ENDIF
             ENDIF

             ! Check for relative change limit
             IF ( maxRatio > 0.0 .AND. minRatio > 0.0 ) THEN 
                frac = SpcNew / MAX(SpcBkg(I,J,L),MinConc)
                ! If change is greater than maximum allowed fraction, restrict to max. fraction
                IF ( frac > maxRatio ) THEN 
                   SpcNew = MAX(SpcBkg(I,J,L),MinConc) * maxRatio
                ! If change is smaller than maximum allowed fraction, restrict to min. fraction
                ELSEIF ( frac < minRatio ) THEN
                   SpcNew = MAX(SpcBkg(I,J,L),MinConc) * minRatio 
                ENDIF
             ENDIF

             ! Update assimilated field
             SpcAsm(I,J,L) = MAX(SpcNew,MinConc)

             ! Update diagnostics
             IF ( ALLOCATED(AnaMask2d ) ) AnaMask2d(I,J)   = AnaMask2d(I,J) + 1.0
             IF ( ALLOCATED(AnaMask3d ) ) AnaMask3d(I,J,L) = 1.0 

             ! Eventually update dependent species to maintain concentration ratio of species 2 / species 1
             IF ( iopt%nSpec2>0 ) THEN
                DO N=1,iopt%nSpec2
                   UpdateSpec2 = .FALSE.
                   ! Default is to use background field
                   Spc2Asm(I,J,L,N) = Spc2Bkg(I,J,L,N)
                   ! Use background field if in stratosphere and no adjustment to be done in stratosphere
                   IF     ( stratwgt >= 0.5 .AND. .NOT. iopt%Spec2(N)%Spec2Strat ) THEN
                      CYCLE 
                   ! Use background field if in troposphere and no adjustment to be done in troposphere
                   ELSEIF ( tropwgt  >= 0.5 .AND. .NOT. iopt%Spec2(N)%Spec2Trop  ) THEN
                      CYCLE
                   ! Calculate Spc2/Spc1 ratio before update and maintain that ratio in assimilation field
                   ELSE    
                      OldRatio = Spc2Bkg(I,J,L,N) / MAX(SpcBkg(I,J,L),MinConc)
                      NewRatio = Spc2Bkg(I,J,L,N) / SpcAsm(I,J,L)
                      ! Update species only if the ratio is within the specified limits. Otherwise, we assume
                      ! that species 2 is so abundant or missing that updating it is not meaningful.
                      IF ( ( OldRatio<iopt%Spec2(N)%Spec2MaxRatio .AND. OldRatio>iopt%Spec2(N)%Spec2MinRatio ) .OR. &
                           ( NewRatio<iopt%Spec2(N)%Spec2MaxRatio .AND. NewRatio>iopt%Spec2(N)%Spec2MinRatio )       ) THEN
                         Spc2Asm(I,J,L,N) = SpcAsm(I,J,L) * OldRatio 
                         UpdateSpec2 = .TRUE.
                      ENDIF
                   ENDIF
                ENDDO
             ENDIF
          ENDDO
       ENDDO
       ENDDO

       ! Print warning if at least one negative cell
       IF ( NNEG > 0 .AND. iopt%Verbose > 2 ) THEN
          WRITE(6,*) '*** DoAnalysis_ warning: encountered concentration below threshold, set to minimum: ',TRIM(SpecName),NNEG,MinConc,' ***'
       ENDIF

       ! Pass back to State_Chm%Species array: convert v/v dry to kg/kg total
       ! -------------------------------------------------------------------------------------------
       IF ( .NOT. iopt%DryRun ) THEN
          State_Chm%Species(indSpc)%Conc(:,:,:) = SpcAsm(:,:,:) * (1.-State_Met%SPHU/1000.) / AIRMW * mwSpc
          IF ( iopt%nSpec2 > 0 ) THEN
             DO N=1,iopt%nSpec2 
                State_Chm%Species(indSpc2(N))%Conc(:,:,:) = Spc2Asm(:,:,:,N) * (1.-State_Met%SPHU/1000.) / AIRMW * mwSpc2(N)
             ENDDO
          ENDIF
       ENDIF

       ! Write out diagnostics as needed
       ! -------------------------------------------------------------------------------------------

       ! AnaInc
       IF ( State_Diag%Archive_AnaInc ) THEN
          DgnID = GetDiagnID ( State_Diag%Map_AnaInc, indSpc )
          IF ( DgnID > 0 ) State_Diag%AnaInc(:,:,:,DgnID) = SpcAsm - SpcBkg
          IF ( iopt%nSpec2 > 0 ) THEN
             DO N=1,iopt%nSpec2
                DgnID = GetDiagnID ( State_Diag%Map_AnaInc, indSpc2(N) )
                IF ( DgnID > 0 ) State_Diag%AnaInc(:,:,:,DgnID) = Spc2Asm(:,:,:,N) - Spc2Bkg(:,:,:,N)
             ENDDO
          ENDIF
       ENDIF

       ! AnaIncFrac
       IF ( State_Diag%Archive_AnaIncFrac ) THEN
          DgnID = GetDiagnID ( State_Diag%Map_AnaIncFrac, indSpc )
          IF ( DgnID > 0 ) State_Diag%AnaIncFrac(:,:,:,DgnID) = SpcAsm / MAX(SpcBkg,MinConc)
          IF ( iopt%nSpec2 > 0 ) THEN
             DO N=1,iopt%nSpec2
                DgnID = GetDiagnID ( State_Diag%Map_AnaIncFrac, indSpc2(N) )
                IF ( DgnID > 0 ) State_Diag%AnaIncFrac(:,:,:,DgnID) = Spc2Asm(:,:,:,N) / MAX(Spc2Bkg(:,:,:,N),MinConc)
             ENDDO
          ENDIF
       ENDIF

       ! AnaMask
       IF ( State_Diag%Archive_AnaMask ) THEN
          DgnID = GetDiagnID ( State_Diag%Map_AnaMask, indSpc )
          IF ( DgnID > 0 ) State_Diag%AnaMask(:,:,:,DgnID) = AnaMask3d(:,:,:)
          IF ( iopt%nSpec2 > 0 ) THEN
             DO N=1,iopt%nSpec2
                DgnID = GetDiagnID ( State_Diag%Map_AnaMask, indSpc2(N) )
                IF ( DgnID > 0 ) State_Diag%AnaMask(:,:,:,DgnID) = AnaMask3d(:,:,:)
             ENDDO
          ENDIF
       ENDIF

       ! AnaMaskSum
       IF ( State_Diag%Archive_AnaMaskSum ) THEN
          DgnID = GetDiagnID ( State_Diag%Map_AnaMaskSum, indSpc )
          IF ( DgnID > 0 ) State_Diag%AnaMaskSum(:,:,DgnID) = AnaMask2d(:,:)
          IF ( iopt%nSpec2 > 0 ) THEN
             DO N=1,iopt%nSpec2
                DgnID = GetDiagnID ( State_Diag%Map_AnaMaskSum, indSpc2(N) )
                IF ( DgnID > 0 ) State_Diag%AnaMaskSum(:,:,DgnID) = AnaMask2d(:,:)
             ENDDO
          ENDIF
       ENDIF

    ENDIF ! HasField

    ! Cleanup
    ! -------
    IF ( ALLOCATED(SpcBkg    ) ) DEALLOCATE(SpcBkg)
    IF ( ALLOCATED(SpcAsm    ) ) DEALLOCATE(SpcAsm)
    IF ( ALLOCATED(Spc2Bkg   ) ) DEALLOCATE(Spc2Bkg)
    IF ( ALLOCATED(Spc2Asm   ) ) DEALLOCATE(Spc2Asm)
    IF ( ALLOCATED(indSpc2   ) ) DEALLOCATE(indSpc2)
    IF ( ALLOCATED(mwSpc2    ) ) DEALLOCATE(mwSpc2)
    IF ( ASSOCIATED(AnaPtr   ) ) DEALLOCATE(AnaPtr)
    IF ( ASSOCIATED(MskPtr   ) ) DEALLOCATE(MskPtr)
    IF ( ALLOCATED(AnaMask2d ) ) DEALLOCATE(AnaMask2d)
    IF ( ALLOCATED(AnaMask3d ) ) DEALLOCATE(AnaMask3d) 

    iopt => NULL()

    ! Successful return
    RC = GC_SUCCESS

  END SUBROUTINE DoAnalysis_ 
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Model                            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: GetAnaTime_ 
!
! !DESCRIPTION: Get analysis time. This is either the current date/time or one
! GEOS-Chem time step ahead, depending on the Fwd input argument.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE GetAnaTime_( Fwd, yy, mm, dd, h, m, s, RC )
!
! !USES:
!
    USE TIME_MOD,  ONLY : GET_TS_CHEM
    USE TIME_MOD,  ONLY : GET_TIME_AHEAD
    USE TIME_MOD,  ONLY : GET_NYMD, GET_NHMS
    USE TIME_MOD,  ONLY : YMD_EXTRACT 
!
! !INPUT/OUTPUT PARAMETERS:
!
    LOGICAL,             INTENT(IN)            :: Fwd        ! Adjust time one time step forward?
    INTEGER,             INTENT(OUT)           :: yy, mm, dd ! year, month, day
    INTEGER,             INTENT(OUT)           :: h,  m,  s  ! hour, minute, second
    INTEGER,             INTENT(OUT)           :: RC         ! Success or failure?
!
! !REMARKS:
!
! !REVISION HISTORY:
!  01 Mar 2022 - C. Keller   - Initial version
!  See https://github.com/geoschem/geos-chem for history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! LOCAL VARIABLES:
!
    REAL                       :: tsChem
    INTEGER                    :: cNYMD, cNHMS
    INTEGER                    :: fDate(2)

    ! Get time, eventually adjust forward
    IF ( Fwd ) THEN
       tsChem = GET_TS_CHEM()
       fDate  = GET_TIME_AHEAD(INT(tsChem))
       cNYMD  = fDate(1)
       cNHMS  = fDate(2)
    ELSE
       cNYMD  = GET_NYMD()
       cNHMS  = GET_NHMS()
    ENDIF
    CALL YMD_EXTRACT( cNYMD, yy, mm, dd ) 
    CALL YMD_EXTRACT( cNHMS, h, m, s ) 

    ! All done
    RC = GC_SUCCESS 

    END SUBROUTINE GetAnaTime_
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Model                            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ReadSettings_ 
!
! !DESCRIPTION: Reads the analysis settings from a given configuration file
!
!
! !INTERFACE:
!
  SUBROUTINE ReadSettings_( am_I_Root, Config, ispec, RC )
!
! !USES:
!
  USE QfYaml_Mod
  USE RoundOff_Mod,  ONLY : Cast_and_RoundOff
!
! !INPUT PARAMETERS:
!
    LOGICAL,              INTENT(IN)    :: am_I_Root  ! Root PET?
    INTEGER,              INTENT(IN)    :: ispec      ! species count 
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(QFYAML_t),       INTENT(INOUT) :: Config     ! yaml config file 
!
! !OUTPUT PARAMETERS:
!
    INTEGER, INTENT(INOUT)              :: RC         ! Success or failure
!
! !REMARKS:
!
! !REVISION HISTORY:
!  25 May 2022 - C. Keller - initial version (refactored Chem_GridCompMod)
!  See https://github.com/geoschem/geos-chem for history
!EOP
!------------------------------------------------------------------------------
!BOC
! 
! !LOCAL VARIABLES:
!
    CHARACTER(LEN=255)            :: Iam
    CHARACTER(LEN=63)             :: pkey, SpecName
    CHARACTER(LEN=3)              :: intStr
    INTEGER                       :: C, ix

    CHARACTER(LEN=QFYAML_NamLen)  :: key
    CHARACTER(LEN=QFYAML_StrLen)  :: v_str
    LOGICAL                       :: found, v_bool
    INTEGER                       :: v_int

    INTEGER                       :: N, nSpec2
    CHARACTER(LEN=255)            :: ThisStr, Spec2Name, Spec2Strat, Spec2Trop 
    CHARACTER(LEN=255)            :: Spec2MinRatio, Spec2MaxRatio

    !=======================================================================
    ! ReadSettings_ begins here 
    !=======================================================================
    Iam = 'GEOS_Analysis:: ReadSettings_'

    ! Prefix for this entry
    WRITE( intStr, '(I3.3)' ) ispec
    pkey = "species%Spc"//TRIM(intStr)

    ! Get species name. Eventually remove single quotes from species
    key = TRIM(pkey)//"%SpeciesName"
    CALL GetKey_( Config, key, RC, vstr=v_str, vstr_default='N/A' )
    IF ( RC /= GC_SUCCESS ) RETURN
    SpecName = TRIM(v_str)
    C = INDEX( SpecName, "'" )
    IF ( C > 0 ) THEN
       SpecName = SpecName(C+1:)
       C = INDEX( SpecName, "'" )
       IF ( C > 0 ) SpecName = Specname(1:C-1)
    ENDIF
    AnaConfig(ispec)%SpecName = SpecName
    IF ( am_I_Root ) write(6,*) 'Reading analysis settings for species: '//TRIM(AnaConfig(ispec)%SpecName)

    ! Active?
    key = TRIM(pkey)//"%Active"
    CALL GetKey_( Config, key, RC, vbool=v_bool, vbool_default=.FALSE. )
    IF ( RC /= GC_SUCCESS ) RETURN
    AnaConfig(ispec)%Active = v_bool

    key = TRIM(pkey)//"%DryRun"
    CALL GetKey_( Config, key, RC, vbool=v_bool, vbool_default=.FALSE. )
    IF ( RC /= GC_SUCCESS ) RETURN
    AnaConfig(ispec)%DryRun = v_bool

    ! Analysis frequency
    key = TRIM(pkey)//"%AnalysisFreq"
    CALL GetKey_( Config, key, RC, vint=v_int, vint_default=6 )
    AnaConfig(ispec)%AnalysisFreq = v_int

    ! Analysis hour 
    key = TRIM(pkey)//"%AnalysisHour"
    CALL GetKey_( Config, key, RC, vint=v_int, vint_default=0 )
    AnaConfig(ispec)%AnalysisHour = v_int

    ! Analysis minute 
    key = TRIM(pkey)//"%AnalysisMinute"
    CALL GetKey_( Config, key, RC, vint=v_int, vint_default=0 )
    IF ( RC /= GC_SUCCESS ) RETURN
    AnaConfig(ispec)%AnalysisMinute = v_int

    ! Forward looking? 
    key = TRIM(pkey)//"%ForwardLooking"
    CALL GetKey_( Config, key, RC, vbool=v_bool, vbool_default=.TRUE. )
    IF ( RC /= GC_SUCCESS ) RETURN
    AnaConfig(ispec)%ForwardLooking = v_bool

    ! Skip predictor? 
    key = TRIM(pkey)//"%SkipPredictor"
    CALL GetKey_( Config, key, RC, vbool=v_bool, vbool_default=.TRUE. )
    IF ( RC /= GC_SUCCESS ) RETURN
    AnaConfig(ispec)%SkipPredictor = v_bool

    ! HEMCO field name 
    key = TRIM(pkey)//"%FldNameHco"
    CALL GetKey_( Config, key, RC, vstr=v_str, vstr_default='N/A' )
    IF ( RC /= GC_SUCCESS ) RETURN
    AnaConfig(ispec)%FldNameHco = v_str

    ! Field units
    key = TRIM(pkey)//"%FileVarUnit"
    CALL GetKey_( Config, key, RC, vstr=v_str, vstr_default='v/v' )
    IF ( RC /= GC_SUCCESS ) RETURN
    AnaConfig(ispec)%FileVarUnit = v_str

    ! Dry air flag 
    key = TRIM(pkey)//"%DryFlag"
    CALL GetKey_( Config, key, RC, vint=v_int, vint_default=1 )
    IF ( RC /= GC_SUCCESS ) RETURN
    AnaConfig(ispec)%DryFlag = v_int

    ! Is increment?
    key = TRIM(pkey)//"%IsIncrement"
    CALL GetKey_( Config, key, RC, vbool=v_bool, vbool_default=.FALSE. )
    IF ( RC /= GC_SUCCESS ) RETURN
    AnaConfig(ispec)%IsIncrement = v_bool

    ! IAU? 
    key = TRIM(pkey)//"%IAU"
    CALL GetKey_( Config, key, RC, vbool=v_bool, vbool_default=.FALSE. )
    IF ( RC /= GC_SUCCESS ) RETURN
    AnaConfig(ispec)%IAU = v_bool

    ! Analysis window length (hours) 
    key = TRIM(pkey)//"%AnalysisWindow"
    CALL GetKey_( Config, key, RC, vint=v_int, vint_default=6 )
    IF ( RC /= GC_SUCCESS ) RETURN
    AnaConfig(ispec)%AnalysisWindow = v_int

    ! Apply in stratosphere? 
    key = TRIM(pkey)//"%InStrat"
    CALL GetKey_( Config, key, RC, vbool=v_bool, vbool_default=.TRUE. )
    IF ( RC /= GC_SUCCESS ) RETURN
    AnaConfig(ispec)%InStrat = v_bool

    ! Apply in troposphere? 
    key = TRIM(pkey)//"%InTrop"
    CALL GetKey_( Config, key, RC, vbool=v_bool, vbool_default=.TRUE. )
    IF ( RC /= GC_SUCCESS ) RETURN
    AnaConfig(ispec)%InTrop = v_bool

    ! Has mask field 
    key = TRIM(pkey)//"%HasMask"
    CALL GetKey_( Config, key, RC, vbool=v_bool, vbool_default=.FALSE. )
    IF ( RC /= GC_SUCCESS ) RETURN
    AnaConfig(ispec)%HasMask = v_bool

    ! Mask name (from HEMCO) 
    key = TRIM(pkey)//"%MskNameHco"
    CALL GetKey_( Config, key, RC, vstr=v_str, vstr_default='N/A' )
    IF ( RC /= GC_SUCCESS ) RETURN
    AnaConfig(ispec)%MskNameHco = v_str

    ! Mask threshold 
    key = TRIM(pkey)//"%MaskThreshold"
    CALL GetKey_( Config, key, RC, vstr=v_str, vstr_default='0.1' )
    IF ( RC /= GC_SUCCESS ) RETURN
    AnaConfig(ispec)%MaskThreshold = Cast_and_RoundOff( v_str, places=2 )

    ! Analysis level 1 
    key = TRIM(pkey)//"%AnaL1"
    CALL GetKey_( Config, key, RC, vint=v_int, vint_default=1 )
    IF ( RC /= GC_SUCCESS ) RETURN
    AnaConfig(ispec)%AnaL1 = v_int

    ! Analysis level 2 
    key = TRIM(pkey)//"%AnaL2"
    CALL GetKey_( Config, key, RC, vint=v_int, vint_default=1 )
    IF ( RC /= GC_SUCCESS ) RETURN
    AnaConfig(ispec)%AnaL2 = v_int

    ! Analysis level 3 
    key = TRIM(pkey)//"%AnaL3"
    CALL GetKey_( Config, key, RC, vint=v_int, vint_default=72 )
    IF ( RC /= GC_SUCCESS ) RETURN
    AnaConfig(ispec)%AnaL3 = v_int

    ! Analysis level 4 
    key = TRIM(pkey)//"%AnaL4"
    CALL GetKey_( Config, key, RC, vint=v_int, vint_default=72 )
    IF ( RC /= GC_SUCCESS ) RETURN
    AnaConfig(ispec)%AnaL4 = v_int

    ! Analysis fraction 
    key = TRIM(pkey)//"%AnaFraction"
    CALL GetKey_( Config, key, RC, vstr=v_str, vstr_default='1.0' )
    IF ( RC /= GC_SUCCESS ) RETURN
    AnaConfig(ispec)%AnaFraction = Cast_and_RoundOff( v_str, places=2 )

    ! Stratosphere sponge layer 
    key = TRIM(pkey)//"%StratSponge"
    CALL GetKey_( Config, key, RC, vint=v_int, vint_default=0 )
    IF ( RC /= GC_SUCCESS ) RETURN
    AnaConfig(ispec)%StratSponge = v_int

    ! Maximum change in stratosphere 
    key = TRIM(pkey)//"%MaxChangeStrat"
    CALL GetKey_( Config, key, RC, vstr=v_str, vstr_default='-1.0' )
    IF ( RC /= GC_SUCCESS ) RETURN
    AnaConfig(ispec)%MaxChangeStrat = Cast_and_RoundOff( v_str, places=2 )

    ! Maximum change in troposphere 
    key = TRIM(pkey)//"%MaxChangeTrop"
    CALL GetKey_( Config, key, RC, vstr=v_str, vstr_default='-1.0' )
    IF ( RC /= GC_SUCCESS ) RETURN
    AnaConfig(ispec)%MaxChangeTrop = Cast_and_RoundOff( v_str, places=2 )

    ! Maximum ratio change in stratosphere 
    key = TRIM(pkey)//"%MaxRatioStrat"
    CALL GetKey_( Config, key, RC, vstr=v_str, vstr_default='-1.0' )
    IF ( RC /= GC_SUCCESS ) RETURN
    AnaConfig(ispec)%MaxRatioStrat = Cast_and_RoundOff( v_str, places=2 )

    ! Maximum ratio change in troposphere 
    key = TRIM(pkey)//"%MaxRatioTrop"
    CALL GetKey_( Config, key, RC, vstr=v_str, vstr_default='-1.0' )
    IF ( RC /= GC_SUCCESS ) RETURN
    AnaConfig(ispec)%MaxRatioTrop = Cast_and_RoundOff( v_str, places=2 )

    ! Minimum ratio change in stratosphere 
    key = TRIM(pkey)//"%MinRatioStrat"
    CALL GetKey_( Config, key, RC, vstr=v_str, vstr_default='-1.0' )
    IF ( RC /= GC_SUCCESS ) RETURN
    AnaConfig(ispec)%MinRatioStrat = Cast_and_RoundOff( v_str, places=2 )

    ! Minimum ratio change in troposphere 
    key = TRIM(pkey)//"%MinRatioTrop"
    CALL GetKey_( Config, key, RC, vstr=v_str, vstr_default='-1.0' )
    IF ( RC /= GC_SUCCESS ) RETURN
    AnaConfig(ispec)%MinRatioTrop = Cast_and_RoundOff( v_str, places=2 )

    ! Minimum concentration 
    key = TRIM(pkey)//"%MinConc"
    CALL GetKey_( Config, key, RC, vstr=v_str, vstr_default='1.0e-20' )
    IF ( RC /= GC_SUCCESS ) RETURN
    READ( v_str, * ) AnaConfig(ispec)%MinConc 
    !AnaConfig(ispec)%MinConc = Cast_and_RoundOff( v_str, places=2 )

    ! Error mode 
    key = TRIM(pkey)//"%ErrorMode"
    CALL GetKey_( Config, key, RC, vint=v_int, vint_default=1 )
    IF ( RC /= GC_SUCCESS ) RETURN
    AnaConfig(ispec)%ErrorMode = v_int

    ! Verbose flag
    key = TRIM(pkey)//"%Verbose"
    CALL GetKey_( Config, key, RC, vint=v_int, vint_default=2 )
    IF ( RC /= GC_SUCCESS ) RETURN
    AnaConfig(ispec)%Verbose = v_int

    ! Check for "dependent" species
    key = TRIM(pkey)//"%HasSpec2"
    CALL GetKey_( Config, key, RC, vint=v_int, vint_default=0 )
    IF ( RC /= GC_SUCCESS ) RETURN
    nSpec2 = v_int
    AnaConfig(ispec)%nSpec2 = nSpec2 

    IF ( nSpec2 > 0 ) THEN
       ALLOCATE(AnaConfig(ispec)%Spec2(nSpec2)) 
       ! Read parameter as string (can be different for multiple dependent species, separated by comma')

       ! species names 
       key = TRIM(pkey)//"%Spec2Name"
       CALL GetKey_( Config, key, RC, vstr=v_str, vstr_default='unknown' )
       IF ( RC /= GC_SUCCESS ) RETURN
       Spec2Name = v_str 
       key = TRIM(pkey)//"%Spec2Strat"
       CALL GetKey_( Config, key, RC, vstr=v_str, vstr_default='1' )
       IF ( RC /= GC_SUCCESS ) RETURN
       Spec2Strat = v_str 
       key = TRIM(pkey)//"%Spec2Trop"
       CALL GetKey_( Config, key, RC, vstr=v_str, vstr_default='1' )
       IF ( RC /= GC_SUCCESS ) RETURN
       Spec2Trop = v_str 
       key = TRIM(pkey)//"%Spec2MinRatio"
       CALL GetKey_( Config, key, RC, vstr=v_str, vstr_default='-1.0' )
       IF ( RC /= GC_SUCCESS ) RETURN
       Spec2MinRatio = v_str 
       key = TRIM(pkey)//"%Spec2MaxRatio"
       CALL GetKey_( Config, key, RC, vstr=v_str, vstr_default='-1.0' )
       IF ( RC /= GC_SUCCESS ) RETURN
       Spec2MaxRatio = v_str 
       ! Assign parameter to various slots
       DO N = 1, nSpec2
          ! Species name
          CALL Spec2Parse_( Spec2Name, N, ThisStr )
          AnaConfig(ispec)%Spec2(N)%Spec2Name = ThisStr 
          ! Use in troposphere/stratosphere?
          CALL Spec2Parse_( Spec2Strat, N, ThisStr )
          AnaConfig(ispec)%Spec2(N)%Spec2Strat = ( TRIM(ThisStr)=='1' )
          CALL Spec2Parse_( Spec2Trop, N, ThisStr )
          AnaConfig(ispec)%Spec2(N)%Spec2Trop  = ( TRIM(ThisStr)=='1' )
          ! Minimum / maximum ratios 
          CALL Spec2Parse_( Spec2MinRatio, N, ThisStr )
          READ(ThisStr,*) AnaConfig(ispec)%Spec2(N)%Spec2MinRatio
          CALL Spec2Parse_( Spec2MaxRatio, N, ThisStr )
          READ(ThisStr,*) AnaConfig(ispec)%Spec2(N)%Spec2MaxRatio
       ENDDO
    ELSE
       AnaConfig(ispec)%Spec2 => NULL()
    ENDIF

    ! Verbose
    IF ( am_I_Root ) THEN
       WRITE(6,*) '----------------------------------------'
       WRITE(6,*) 'Analysis settings for GEOS-Chem species ',TRIM(AnaConfig(ispec)%SpecName),':'
       WRITE(6,*) 'Active: ',AnaConfig(ispec)%Active
       IF ( AnaConfig(ispec)%Active ) THEN
          WRITE(6,*) '- Dry run mode                  : ',AnaConfig(ispec)%DryRun
          WRITE(6,*) '- Analysis frequency            : ',AnaConfig(ispec)%AnalysisFreq
          WRITE(6,*) '- Analysis hour                 : ',AnaConfig(ispec)%AnalysisHour
          WRITE(6,*) '- Analysis minute               : ',AnaConfig(ispec)%AnalysisMinute
          WRITE(6,*) '- Forward looking file read     : ', AnaConfig(ispec)%ForwardLooking
          WRITE(6,*) '- Ignore during predictor step  : ', AnaConfig(ispec)%SkipPredictor
          WRITE(6,*) '- HEMCO field name              : ', TRIM(AnaConfig(ispec)%FldNameHco)   
          WRITE(6,*) '- Input variable unit           : ', TRIM(AnaConfig(ispec)%FileVarUnit)
          WRITE(6,*) '- Dry air flag (0=dry, 1=total) : ', AnaConfig(ispec)%DryFlag
          WRITE(6,*) '- Field are increments          : ', AnaConfig(ispec)%IsIncrement
          WRITE(6,*) '- Spread ana/inc (IAU)          : ', AnaConfig(ispec)%IAU
          WRITE(6,*) '- Analysis window length [h]    : ', AnaConfig(ispec)%AnalysisWindow
          WRITE(6,*) '- Restrict analysis to mask area: ', AnaConfig(ispec)%HasMask
          IF ( AnaConfig(ispec)%HasMask ) THEN
             WRITE(6,*) '- Mask field name: ', TRIM(AnaConfig(ispec)%MskNameHco)
             WRITE(6,*) '- Mask threshold : ', AnaConfig(ispec)%MaskThreshold
          ENDIF
          WRITE(6,*) '- Apply analysis in stratosphere: ', AnaConfig(ispec)%InStrat
          WRITE(6,*) '- Apply analysis in troposphere : ', AnaConfig(ispec)%InTrop
          WRITE(6,*) '- Tropopause sponge layer       : ', AnaConfig(ispec)%StratSponge
          WRITE(6,*) '- Analysis level 1              : ', AnaConfig(ispec)%AnaL1  
          WRITE(6,*) '- Analysis level 2              : ', AnaConfig(ispec)%AnaL2
          WRITE(6,*) '- Analysis level 3              : ', AnaConfig(ispec)%AnaL3
          WRITE(6,*) '- Analysis level 4              : ', AnaConfig(ispec)%AnaL4
          WRITE(6,*) '- Analysis fraction             : ', AnaConfig(ispec)%AnaFraction
          WRITE(6,*) '- Max. absolute change in strat : ', AnaConfig(ispec)%MaxChangeStrat
          WRITE(6,*) '- Max. absolute change in trop  : ', AnaConfig(ispec)%MaxChangeTrop  
          WRITE(6,*) '- Max. relative change in strat : ', AnaConfig(ispec)%MaxRatioStrat
          WRITE(6,*) '- Max. relative change in trop  : ', AnaConfig(ispec)%MaxRatioTrop  
          WRITE(6,*) '- Min. relative change in strat : ', AnaConfig(ispec)%MinRatioStrat
          WRITE(6,*) '- Min. relative change in trop  : ', AnaConfig(ispec)%MinRatioTrop  
          WRITE(6,*) '- Min. concentration (for ratio): ', AnaConfig(ispec)%MinConc
          WRITE(6,*) '- # of dependent species        : ', AnaConfig(ispec)%nSpec2
          IF ( AnaConfig(ispec)%nSpec2 > 0 ) THEN
             DO N = 1, AnaConfig(ispec)%nSpec2
                WRITE(6,*) '- Name of dependent species     : ', TRIM(AnaConfig(ispec)%Spec2(N)%Spec2Name)
                WRITE(6,*) '- Update in stratosphere        : ', AnaConfig(ispec)%Spec2(N)%Spec2Strat
                WRITE(6,*) '- Update in troposphere         : ', AnaConfig(ispec)%Spec2(N)%Spec2Trop 
                WRITE(6,*) '- Minimum ratio (x/parent)      : ', AnaConfig(ispec)%Spec2(N)%Spec2MinRatio
                WRITE(6,*) '- Maximum ratio (x/parent)      : ', AnaConfig(ispec)%Spec2(N)%Spec2MaxRatio
             ENDDO
          ENDIF
          WRITE(6,*) '- Error mode                    : ', AnaConfig(ispec)%ErrorMode
          WRITE(6,*) '- Verbose flag                  : ', AnaConfig(ispec)%Verbose
       ENDIF ! Active 
    ENDIF

  END SUBROUTINE ReadSettings_ 
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Model                            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Spec2Parse_ 
!
! !DESCRIPTION: Helper routine to get the Nth index of a string (separated
!  by comma).  
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Spec2Parse_( instr, N, outstr )
!
! !INPUT PARAMETERS:
!
    CHARACTER(LEN=*), INTENT(IN)    :: instr 
    INTEGER,          INTENT(IN)    :: N 
!                                                             
! !INPUT/OUTPUT PARAMETERS:                                         
!              
    CHARACTER(LEN=*), INTENT(INOUT) :: outstr 
!
! !REVISION HISTORY:
!  07 Jul 2022 - C. Keller   - Initial version
!  See https://github.com/geoschem/geos-chem for history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! LOCAL VARIABLES:
!
    INTEGER                    :: lenstr
    INTEGER                    :: I, LPOS, RPOS
    CHARACTER(LEN=255)         :: tmpstr

    ! Local copy of original string
    tmpstr = TRIM(instr)
    lenstr = LEN(instr)

    ! Check for first separator in string
    LPOS = 0
    RPOS = INDEX(TRIM(tmpstr),',')

    ! Use full string if no separator in string 
    IF ( RPOS>0 .AND. N>1 ) THEN
       DO I = 2, N
          IF(RPOS<=0) CYCLE
          LPOS = RPOS
          tmpstr(LPOS:LPOS) = '.'
          RPOS = INDEX(TRIM(tmpstr),',')
       ENDDO
    ENDIF
    IF ( RPOS <= 0 ) RPOS = lenstr+1 
    outstr = TRIM(instr(LPOS+1:RPOS-1))

  END SUBROUTINE Spec2Parse_ 
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Model                            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ReplaceChar_ 
!
! !DESCRIPTION: Replaces all characters in a string. 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE ReplaceChar_ ( str, pattern, replace )
!
! !INPUT PARAMETERS:
!
    CHARACTER(LEN=*), INTENT(IN)    :: pattern 
    CHARACTER(LEN=*), INTENT(IN)    :: replace
!                                                             
! !INPUT/OUTPUT PARAMETERS:                                         
!              
    CHARACTER(LEN=*), INTENT(INOUT) :: str 
!
! !REVISION HISTORY:
!  20 Dec 2018 - C. Keller   - Initial version
!  See https://github.com/geoschem/geos-chem for history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! LOCAL VARIABLES:
!
    INTEGER  :: I, LP, LS

    LP = LEN(TRIM(pattern))
    LS = LEN(TRIM(str))+1
    I  = INDEX(TRIM(str),TRIM(pattern))
    DO WHILE ( I > 0 ) 
       str = TRIM(str(1:I-1))//TRIM(replace)//TRIM(str(I+LP:LS))
       I = INDEX(TRIM(str),TRIM(pattern))
    ENDDO

  END SUBROUTINE ReplaceChar_
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: GetDiagnID 
!
! !DESCRIPTION: Helper function to get the diagnostics ID 
!\\
!\\
! !INTERFACE:
!
  FUNCTION GetDiagnID ( mapDgn, SpcID ) RESULT ( SlotNumber ) 
!
! USES:
! 
  USE State_Diag_Mod,        ONLY : DgnMap 
!
! INPUTS:
!
    TYPE(DgnMap), POINTER      :: mapDgn
    INTEGER,      INTENT(IN)   :: SpcID
!
! !RETURN VALUE:
!
    INTEGER                    :: SlotNumber
!
! !REVISION HISTORY:
!  14 Dec 2023 - C. Keller - Initial version.
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
! 
    INTEGER :: N

    SlotNumber = -1
    DO N = 1,mapDgn%nSlots
       IF ( mapDgn%slot2id(N) == SpcID ) THEN
          SlotNumber = N
          EXIT
       ENDIF
    ENDDO

  END FUNCTION GetDiagnID
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Model                            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: GetKey_ 
!
! !DESCRIPTION: Helper routine to get key 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE GetKey_( Config, key, RC, vint, vint_default, vbool, vbool_default, vstr, vstr_default ) 
!
! !USES:
!
  USE QfYaml_Mod
!
! !INPUT PARAMETERS:
!
    TYPE(QFYAML_t),    INTENT(INOUT) :: Config     ! yaml config file 
    CHARACTER(LEN=*),  INTENT(IN)    :: key 
!                                                             
! !INPUT/OUTPUT PARAMETERS:                                         
!              
    INTEGER,           INTENT(INOUT)                      :: RC 
    INTEGER, OPTIONAL, INTENT(OUT)                        :: vint
    INTEGER, OPTIONAL, INTENT(IN)                         :: vint_default
    LOGICAL, OPTIONAL, INTENT(OUT)                        :: vbool
    LOGICAL, OPTIONAL, INTENT(IN)                         :: vbool_default
    CHARACTER(LEN=QFYAML_StrLen), OPTIONAL, INTENT(OUT)   :: vstr
    CHARACTER(LEN=*), OPTIONAL, INTENT(IN)                :: vstr_default
!
! !REVISION HISTORY:
!  15 Dec 2023 - C. Keller   - Initial version
!  See https://github.com/geoschem/geos-chem for history
!EOP
!------------------------------------------------------------------------------
!BOC
!   
! LOCAL VARIABLES:
!   
    INTEGER :: v_int
    LOGICAL :: v_bool
    CHARACTER(LEN=QFYAML_StrLen) :: v_str
    CHARACTER(LEN=255) :: errMsg
    CHARACTER(LEN=255), PARAMETER :: Iam = 'GetKey_'

    ! Starts here
    RC = GC_SUCCESS

    ! integer
    IF ( PRESENT(vint) ) THEN
       v_int = vint_default
       CALL QFYAML_Add_Get( Config, TRIM( key ), v_int, "", RC )
       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = 'Error parsing ' // TRIM( key ) // '!'
          CALL GC_Error( errMsg, RC, Iam )
          RETURN
       ENDIF
       vint = v_int
    ENDIF
    
    ! bool 
    IF ( PRESENT(vbool) ) THEN
       v_bool = vbool_default
       CALL QFYAML_Add_Get( Config, TRIM( key ), v_bool, "", RC )
       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = 'Error parsing ' // TRIM( key ) // '!'
          CALL GC_Error( errMsg, RC, Iam )
          RETURN
       ENDIF
       vbool = v_bool
    ENDIF
    
    ! character 
    IF ( PRESENT(vstr) ) THEN
       v_str = TRIM(vstr_default)
       CALL QFYAML_Add_Get( Config, TRIM( key ), v_str, "", RC )
       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = 'Error parsing ' // TRIM( key ) // '!'
          CALL GC_Error( errMsg, RC, Iam )
          RETURN
       ENDIF
       vstr = v_str
    ENDIF
       
  END SUBROUTINE GetKey_ 
!EOC
END MODULE GEOS_Analysis 
