#include "MAPL_Generic.h"

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
  USE ESMF     
  USE MAPL_Mod 

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
  PRIVATE  :: GetAnaBundle_
  PRIVATE  :: ReplaceChar_ 
!
! !PRIVATE TYPES:
!
  ! Number of species with analysis on
  INTEGER                         :: nAnaSpec 

  ! Analysis options object. A separate object will be created for each analysed species/family
  TYPE AnaOptions
     CHARACTER(LEN=ESMF_MAXSTR)   :: SpecName
     LOGICAL                      :: Active
     INTEGER                      :: AnalysisFreq
     INTEGER                      :: AnalysisHour
     INTEGER                      :: AnalysisMinute
     LOGICAL                      :: ForwardLooking
     LOGICAL                      :: ReadAnaTime 
     CHARACTER(LEN=ESMF_MAXSTR)   :: FileTemplate 
     CHARACTER(LEN=ESMF_MAXSTR)   :: FileVarName
     CHARACTER(LEN=ESMF_MAXSTR)   :: FileVarUnit
     LOGICAL                      :: ApplyIncrement
     LOGICAL                      :: NonZeroIncOnly
     CHARACTER(LEN=ESMF_MAXSTR)   :: FileVarNameInc
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
     LOGICAL                      :: UseObsHour
     CHARACTER(LEN=ESMF_MAXSTR)   :: ObsHourName 
     LOGICAL                      :: HasSpec2
     CHARACTER(LEN=ESMF_MAXSTR)   :: Spec2Name
     LOGICAL                      :: Spec2Strat 
     LOGICAL                      :: Spec2Trop 
     REAL                         :: Spec2MinRatio
     REAL                         :: Spec2MaxRatio
     INTEGER                      :: ErrorMode
  END TYPE AnaOptions

  ! List holding all analysis information
  TYPE(AnaOptions), POINTER  :: AnaConfig(:) => NULL()
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
  SUBROUTINE GEOS_AnaInit( am_I_Root, GC, GEOSCF, AnaPhase, RC )
!
! !USES:
!
!
! !INPUT PARAMETERS:
!
    LOGICAL,              INTENT(IN)    :: am_I_Root ! Root PET?
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ESMF_GridComp),  INTENT(INOUT) :: GC        ! GridComp 
    TYPE(ESMF_Config),    INTENT(INOUT) :: GEOSCF    ! GEOSCHEMchem_GridComp.rc
!
! !OUTPUT PARAMETERS:
!
    INTEGER, INTENT(OUT)                :: AnaPhase  ! Do analysis after run phase 1 or 2?
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
! 
! !LOCAL VARIABLES:
!
    TYPE(ESMF_Config)             :: AnaSpecCF
    CHARACTER(LEN=ESMF_MAXSTR)    :: compName
    CHARACTER(LEN=ESMF_MAXSTR)    :: Iam
    CHARACTER(LEN=ESMF_MAXSTR)    :: ConfigName
    CHARACTER(LEN=ESMF_MAXSTR)    :: SpecName
    INTEGER                       :: I, N, NDIAG, ThisInt
    INTEGER                       :: STATUS 

    !=======================================================================
    ! GEOS_AnaInit begins here 
    !=======================================================================

    ! Get configuration
    CALL ESMF_GridCompGet( GC, name=compName, __RC__ )

    ! callback name
    Iam = TRIM(compName)//'::GEOS_AnaInit'

    ! Run phase after which to apply analysis
    CALL ESMF_ConfigGetAttribute( GEOSCF, ThisInt, Label="ANAPHASE:", Default=2, __RC__ )
    ANAPHASE = ThisInt

    ! Get number of analysis species
    CALL ESMF_ConfigGetAttribute( GEOSCF, nAnaSpec, Label="Analysis_nSpecies:", Default=0, __RC__ )
    IF ( am_I_Root ) THEN
       WRITE(*,*) 'Number of analysis species: ',nAnaSpec
       WRITE(*,*) 'Analysis phase set to ',ANAPHASE
    ENDIF

    ! Read settings for all species from config file
    IF ( nAnaSpec > 0 ) THEN
       ALLOCATE( AnaConfig(nAnaSpec), STAT=STATUS )
       _ASSERT( STATUS==0, 'AnaConfig could not be allocated' )

       DO N=1,nAnaSpec
          CALL ReadSettings_( am_I_Root, GEOSCF, N, __RC__ )
       ENDDO
    ENDIF

    ! Initialize diagnostics
    IF ( nAnaSpec > 0 ) THEN
       DO N=1,nAnaSpec
          NDIAG = 1
          IF ( AnaConfig(N)%HasSpec2 ) NDIAG = 2 
          DO I=1,NDIAG
             IF ( I==1 ) THEN
                SpecName = AnaConfig(N)%SpecName
             ELSE
                SpecName = AnaConfig(N)%Spec2Name
             ENDIF
             CALL MAPL_AddExportSpec(GC,                                                                   &
                   SHORT_NAME         = 'GCC_ANA_INC_'//TRIM(SpecName),                                    &
                   LONG_NAME          = TRIM(SpecName)//'_analysis_increment_volume_mixing_ratio_dry_air', &
                   UNITS              = 'mol mol-1',                                                       &
                   DIMS               = MAPL_DimsHorzVert,                                                 &
                   VLOCATION          = MAPL_VLocationCenter,                                              &
                                                             __RC__ )
             CALL MAPL_AddExportSpec(GC,                                                                   &
                   SHORT_NAME         = 'GCC_ANA_INC_FRAC_'//TRIM(SpecName),                               &
                   LONG_NAME          = TRIM(SpecName)//'_analysis_increment_ratio_volume_mixing_ratio_dry_air', &
                   UNITS              = '1',                                                               &
                   DIMS               = MAPL_DimsHorzVert,                                                 &
                   VLOCATION          = MAPL_VLocationCenter,                                              &
                                                             __RC__ )
             IF ( I==1 ) THEN 
                CALL MAPL_AddExportSpec(GC,                                        &
                      SHORT_NAME         = 'GCC_ANA_MASK_VSUM_'//TRIM(SpecName),   &
                      LONG_NAME          = TRIM(SpecName)//'_analysis_counts',     &
                      UNITS              = '1',                                    &
                      DIMS               = MAPL_DimsHorzOnly,                      &
                      VLOCATION          = MAPL_VLocationNone,                     &
                                                                __RC__ )
                CALL MAPL_AddExportSpec(GC,                                   &
                      SHORT_NAME         = 'GCC_ANA_MASK_'//TRIM(SpecName),   &
                      LONG_NAME          = TRIM(SpecName)//'_analysis_mask',  &
                      UNITS              = '1',                               &
                      DIMS               = MAPL_DimsHorzVert,                 &
                      VLOCATION          = MAPL_VLocationCenter,              &
                                                                __RC__ )
             ENDIF
             IF ( I==2 ) THEN
                CALL MAPL_AddExportSpec(GC,                                   &
                      SHORT_NAME         = 'GCC_ANA_RATIO_'//TRIM(SpecName)//'_TO_'//TRIM(AnaConfig(N)%SpecName),   &
                      LONG_NAME          = TRIM(SpecName)//'_to_'//TRIM(AnaConfig(N)%SpecName)//'_species_ratio_after_analysis',  &
                      UNITS              = '1',                               &
                      DIMS               = MAPL_DimsHorzVert,                 &
                      VLOCATION          = MAPL_VLocationCenter,              &
                                                                __RC__ )
             ENDIF
          ENDDO
       ENDDO
    ENDIF

    ! Successful return
    RETURN_(ESMF_SUCCESS) 

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
  SUBROUTINE GEOS_AnaRun( GC, Import, Internal, Export, Clock, &
                          Input_Opt,  State_Met, State_Chm, Q, PLE, TROPP, RC )
!
! !USES:
!
  USE Input_Opt_Mod,         ONLY : OptInput
  USE State_Chm_Mod,         ONLY : ChmState         ! Chemistry State obj
  USE State_Met_Mod,         ONLY : MetState         ! Meteorology State obj
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ESMF_GridComp), INTENT(INOUT)         :: GC       ! Ref. to this GridComp
    TYPE(ESMF_State),    INTENT(INOUT)         :: Import   ! Import State
    TYPE(ESMF_STATE),    INTENT(INOUT)         :: Internal ! Internal state
    TYPE(ESMF_State),    INTENT(INOUT)         :: Export   ! Export State
    TYPE(ESMF_Clock),    INTENT(INOUT)         :: Clock    ! ESMF Clock object
    TYPE(OptInput)                             :: Input_Opt
    TYPE(MetState)                             :: State_Met
    TYPE(ChmState)                             :: State_Chm
    REAL,                INTENT(INOUT)         :: Q(:,:,:)
    REAL,                POINTER               :: PLE(:,:,:)
    REAL,                POINTER               :: TROPP(:,:)
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
! 
! !LOCAL VARIABLES:
!
    INTEGER                       :: ispec
    CHARACTER(LEN=ESMF_MAXSTR)    :: compName
    CHARACTER(LEN=ESMF_MAXSTR)    :: Iam
    INTEGER                       :: STATUS 

    !=======================================================================
    ! GEOS_AnaRun begins here 
    !=======================================================================

    ! Do nothing if

    ! Get configuration
    CALL ESMF_GridCompGet( GC, name=compName, __RC__ )

    ! callback name
    Iam = TRIM(compName)//'::GEOS_AnaRun'

    ! Do analysis for all analysis species
    IF ( nAnaSpec > 0 .AND. ASSOCIATED(AnaConfig) ) THEN 
       DO ispec=1,nAnaSpec
          IF ( AnaConfig(ispec)%Active ) THEN
             CALL DoAnalysis_( GC, Import, Internal, Export, Clock, ispec, &
                               Input_Opt,  State_Met, State_Chm, Q, PLE, TROPP, __RC__ )
          ENDIF
       ENDDO
    ENDIF

    ! Successful return
    RETURN_(ESMF_SUCCESS) 

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
    CHARACTER(LEN=ESMF_MAXSTR)    :: Iam
    INTEGER                       :: STATUS 

    !=======================================================================
    ! GEOS_AnaInit begins here 
    !=======================================================================

    ! callback name
    Iam = 'GEOS_AnaFinal'

    ! Clean up 
    IF ( ASSOCIATED(AnaConfig) ) DEALLOCATE( AnaConfig )
    AnaConfig => NULL()
    nAnaSpec = 0

    ! Successful return
    RETURN_(ESMF_SUCCESS) 

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
  SUBROUTINE DoAnalysis_( GC, Import, Internal, Export, Clock, ispec, &
                         Input_Opt,  State_Met, State_Chm, Q, PLE, TROPP, RC )
!
! !USES:
!
  USE Input_Opt_Mod,         ONLY : OptInput
  USE State_Chm_Mod,         ONLY : ChmState         ! Chemistry State obj
  USE State_Met_Mod,         ONLY : MetState         ! Meteorology State obj
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ESMF_GridComp), INTENT(INOUT)         :: GC        ! Ref. to this GridComp
    TYPE(ESMF_State),    INTENT(INOUT)         :: Import    ! Import State
    TYPE(ESMF_STATE),    INTENT(INOUT)         :: Internal  ! Internal state
    TYPE(ESMF_State),    INTENT(INOUT)         :: Export    ! Export State
    TYPE(ESMF_Clock),    INTENT(INOUT)         :: Clock     ! ESMF Clock object
    INTEGER,             INTENT(IN)            :: ispec     ! analysis species index
    TYPE(OptInput)                             :: Input_Opt
    TYPE(MetState)                             :: State_Met
    TYPE(ChmState)                             :: State_Chm
    REAL,                INTENT(INOUT)         :: Q(:,:,:)
    REAL,                POINTER               :: PLE(:,:,:)
    REAL,                POINTER               :: TROPP(:,:)
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
    CHARACTER(LEN=ESMF_MAXSTR) :: compName
    TYPE(ESMF_Grid)            :: grid
    LOGICAL                    :: am_I_Root
    TYPE(AnaOptions), POINTER  :: iopt => NULL() 
    LOGICAL                    :: TimeForAna, HasBundle
    INTEGER                    :: yy, mm, dd, h, m, s
    INTEGER                    :: VarID
    INTEGER                    :: StratCount
    REAL                       :: ThisHour
    CHARACTER(LEN=ESMF_MAXSTR) :: SpecName, Spec2Name, FldName
    REAL, POINTER              :: DiagInc(:,:,:),  DiagIncFrac(:,:,:)
    REAL, POINTER              :: DiagInc2(:,:,:), DiagIncFrac2(:,:,:), DiagSpcRatio(:,:,:)
    REAL, POINTER              :: DiagMsk2d(:,:),  DiagMsk3d(:,:,:)
    REAL, POINTER              :: AnaPtr(:,:,:), IncPtr(:,:,:), ObsHour(:,:)
    REAL, ALLOCATABLE          :: SpcBkg(:,:,:), SpcAsm(:,:,:)
    REAL, ALLOCATABLE          :: Spc2Bkg(:,:,:), Spc2Asm(:,:,:)
    TYPE(MAPL_SimpleBundle)    :: VarBundle, VarBundleH
    CHARACTER(LEN=ESMF_MAXSTR) :: ifile, only_vars
    TYPE(ESMF_TIME)            :: fileTime
    INTEGER                    :: I, J, L, IM, JM, LM, LB, indSpc, indSpc2
    INTEGER                    :: UnitFlag, NNEG
    REAL                       :: Spc2toSpc, wgt, tropwgt, stratwgt
    REAL                       :: frac, diff, maxChange, maxRatio, minRatio
    REAL                       :: mwSpc, mwSpc2
    REAL                       :: SpcAna, SpcNew
    REAL                       :: MinConc
    LOGICAL                    :: UpdateSpec2

    CHARACTER(LEN=ESMF_MAXSTR) :: Iam
    INTEGER                    :: STATUS 

    !=======================================================================
    ! DoAnalysis_ begins here 
    !=======================================================================

    Iam = 'GEOS_Analysis::DoAnalysis_'

    ! Get configuration
    CALL ESMF_GridCompGet( GC, name=compName, grid=grid, __RC__ )

    ! Root CPU?
    am_I_Root = MAPL_am_I_Root()

    ! Get settings
    iopt => AnaConfig(ispec)
    SpecName = iopt%SpecName
    IF ( iopt%HasSpec2 ) THEN
       Spec2Name = iopt%Spec2Name 
    ENDIF
    MinConc   = iopt%MinConc

    ! Check if it's time to do the analysis
    TimeForAna = .FALSE. 
    CALL GetAnaTime_( Clock, iopt%ForwardLooking, yy, mm, dd, h, m, s, __RC__ )
    ThisHour = real(h)

    ! If using observation hours, apply analysis every (full) hour
    IF ( iopt%UseObsHour ) THEN
       IF ( m==0 ) TimeForAna = .TRUE.
    ! Otherwise, use specified analysis frequency and hour/minute offsets
    ELSE 
       IF ( m==iopt%AnalysisMinute .AND. MOD(h,iopt%AnalysisFreq)==iopt%AnalysisHour ) TimeForAna = .TRUE. 
    ENDIF

    ! Initialize diagnostics
    ! ----------------------
    FldName = 'GCC_ANA_INC_'//TRIM(SpecName)
    CALL MAPL_GetPointer ( Export, DiagInc, TRIM(FldName), NotFoundOk=.TRUE., __RC__ )
    IF ( ASSOCIATED(DiagInc) ) DiagInc = 0.0 
    FldName = 'GCC_ANA_INC_FRAC_'//TRIM(SpecName)
    CALL MAPL_GetPointer ( Export, DiagIncFrac, TRIM(FldName), NotFoundOk=.TRUE., __RC__ )
    IF ( ASSOCIATED(DiagIncFrac) ) DiagIncFrac = 1.0 
    IF ( iopt%HasSpec2 ) THEN
       FldName = 'GCC_ANA_INC_'//TRIM(Spec2Name)
       CALL MAPL_GetPointer ( Export, DiagInc2, TRIM(FldName), NotFoundOk=.TRUE., __RC__ )
       IF ( ASSOCIATED(DiagInc2) ) DiagInc2 = 0.0 
       FldName = 'GCC_ANA_INC_FRAC_'//TRIM(Spec2Name)
       CALL MAPL_GetPointer ( Export, DiagIncFrac2, TRIM(FldName), NotFoundOk=.TRUE., __RC__ )
       IF ( ASSOCIATED(DiagIncFrac2) ) DiagIncFrac2 = 1.0 
       FldName = 'GCC_ANA_RATIO_'//TRIM(Spec2Name)//'_TO_'//TRIM(SpecName)
       CALL MAPL_GetPointer ( Export, DiagSpcRatio, TRIM(FldName), NotFoundOk=.TRUE., __RC__ )
       IF ( ASSOCIATED(DiagIncFrac2) ) DiagSpcRatio = -999.0 
    ENDIF
    FldName = 'GCC_ANA_MASK_VSUM_'//TRIM(SpecName)
    CALL MAPL_GetPointer ( Export, DiagMsk2d, TRIM(FldName), NotFoundOk=.TRUE., __RC__ )
    IF ( ASSOCIATED(DiagMsk2d) ) DiagMsk2d = 0.0 
    FldName = 'GCC_ANA_MASK_'//TRIM(SpecName)
    CALL MAPL_GetPointer ( Export, DiagMsk3d, TRIM(FldName), NotFoundOk=.TRUE., __RC__ )
    IF ( ASSOCIATED(DiagMsk3d) ) DiagMsk3d = 0.0 

    ! Check if file exists (only if it's time to do the analysis
    HasBundle = .FALSE.
    IF ( TimeForAna ) THEN
       only_vars = TRIM(iopt%FileVarName)
       IF ( iopt%NonZeroIncOnly ) only_vars = TRIM(only_vars)//','//TRIM(iopt%FileVarNameInc)
       CALL GetAnaBundle_( am_I_Root, iopt%FileTemplate, 'AnaFld', yy, mm, dd, h, m, grid, &
                           VarBundle, HasBundle, ifile=ifile, fileTime=fileTime,    &
                           only_vars=only_vars, err_mode=iopt%ErrorMode, anatime=iopt%ReadAnaTime, __RC__ )

       ! Read obs time using voting regridding method 
       IF ( HasBundle .AND. iopt%UseObsHour ) THEN
          VarBundleH =  MAPL_SimpleBundleRead ( TRIM(ifile), 'AnaHour', grid, fileTime, &
                                                ONLY_VARS=TRIM(iopt%ObsHourName), voting=.TRUE., __RC__ )
!                                               ONLY_VARS=TRIM(iopt%ObsHourName), regrid_method=REGRID_METHOD_VOTE, __RC__ )
       ENDIF
    ENDIF

    ! Apply increments if it's time to do so and if file exists
    ! ---------------------------------------------------------
    IF ( HasBundle ) THEN 

       ! Verbose
       IF ( am_I_Root ) THEN
          WRITE(*,100) SpecName,yy,mm,dd,h,m
100       FORMAT( "GEOS-Chem: apply analysis for species ",A5," for ",I4.4,"-",I2.2,"-",I2.2," ",I2.2,":",I2.2)
       ENDIF

       ! Get analysis field 
       VarID   = MAPL_SimpleBundleGetIndex ( VarBundle, TRIM(iopt%FileVarName), 3, RC=STATUS, QUIET=.TRUE. )
       ASSERT_(RC==ESMF_SUCCESS .AND. VarID > 0)
       AnaPtr => VarBundle%r3(VarID)%q
       IF ( iopt%nonZeroIncOnly ) THEN
          VarID   = MAPL_SimpleBundleGetIndex ( VarBundle, TRIM(iopt%FileVarNameInc), 3, RC=STATUS, QUIET=.TRUE. )
          ASSERT_(RC==ESMF_SUCCESS .AND. VarID > 0)
          IncPtr => VarBundle%r3(VarID)%q
       ENDIF
       ! Observation hour
       ObsHour => NULL()
       IF ( iopt%UseObsHour ) THEN
          VarID   =  MAPL_SimpleBundleGetIndex ( VarBundleH, TRIM(iopt%ObsHourName), 2, RC=STATUS, QUIET=.TRUE. )
          ASSERT_(RC==ESMF_SUCCESS .AND. VarID > 0)
          ObsHour => VarBundleH%r2(VarID)%q
       ENDIF

       ! Select GEOS-Chem index and molecular weight for analysis species. Also get the same for 2nd species (if used) 
       indSpc   = -1
       indSpc2  = -1
       DO I = 1, State_Chm%nSpecies
          IF ( TRIM(State_Chm%SpcData(I)%Info%Name) == TRIM(SpecName) ) THEN
             indSpc = I
             mwSpc  = State_Chm%SpcData(I)%Info%MW_g
          ENDIF
          IF ( iopt%HasSpec2 ) THEN
             IF ( TRIM(State_Chm%SpcData(I)%Info%Name) == TRIM(Spec2Name) ) THEN
                indSpc2 = I
                mwSpc2  = State_Chm%SpcData(I)%Info%MW_g
             ENDIF
          ENDIF
       ENDDO
       ASSERT_(indSpc > 0  )
       ASSERT_( mwSpc > 0.0)
       IF ( iopt%HasSpec2 ) THEN
          ASSERT_(indSpc2 > 0  )
          ASSERT_( mwSpc2 > 0.0)
       ENDIF

       ! array dimensions 
       IM = SIZE(AnaPtr,1)
       JM = SIZE(AnaPtr,2)
       LM = SIZE(AnaPtr,3)

       ! Get lower bound of PLE array
       LB = LBOUND(PLE,3)

       ! Set unit flag. This is to prevent parsing the unit string within the loop below
       SELECT CASE ( TRIM(iopt%FileVarUnit) )
          CASE ( 'kg/kg' )
             UnitFlag = 1
          CASE ( 'mol/mol', 'v/v' )
             UnitFlag = 2
          CASE ( 'ppmv', 'ppm', 'PPMV', 'PPM' )
             UnitFlag = 3
          CASE ( 'ppbv', 'ppb', 'PPBV', 'PPB' )
             UnitFlag = 4
          CASE DEFAULT
             UnitFlag = 1
       END SELECT

       ! State_Chm%Species are in kg/kg total. Make local copy in v/v dry before applying increments.
       ! Also flip vertical axis to be consistent with GEOS
       ALLOCATE(SpcBkg(IM,JM,LM),SpcAsm(IM,JM,LM))
       SpcBkg(:,:,:) = State_Chm%Species(indSpc)%Conc(:,:,LM:1:-1) / (1.-Q) * MAPL_AIRMW / mwSpc
       ! testing only
       !SpcBkg(:,:,:) = State_Chm%Species(indSpc)%Conc(:,:,LM:1:-1)
       SpcAsm(:,:,:) = SpcBkg(:,:,:)
       IF ( iopt%HasSpec2 ) THEN
          ALLOCATE(Spc2Bkg(IM,JM,LM),Spc2Asm(IM,JM,LM))
          Spc2Bkg(:,:,:) = State_Chm%Species(indSpc2)%Conc(:,:,LM:1:-1) / (1.-Q) * MAPL_AIRMW / mwSpc2
          Spc2Asm(:,:,:) = Spc2Bkg(:,:,:)
          IF ( ASSOCIATED(DiagSpcRatio) ) THEN
             WHERE ( SpcBkg > MinConc ) 
                DiagSpcRatio = Spc2Bkg / SpcBkg
             ELSEWHERE
                DiagSpcRatio = Spc2Bkg / MinConc 
             ENDWHERE
          ENDIF
       ENDIF
 
       ! Number of negative cells
       NNEG = 0

       DO J=1,JM
       DO I=1,IM
          ! Move to next grid box if there was no observation in this cell for the given hour and the obshour flag is on 
          IF ( iopt%UseObsHour ) THEN
             IF ( ObsHour(I,J) /= ThisHour ) CYCLE
          ENDIF

          ! Loop over vertical
          StratCount = 0
          DO L=LM,1,-1

             ! Fraction of cell in troposphere / stratosphere 
             tropwgt  = MAX(0.0,MIN(1.0,(PLE(I,J,L+LB)-TROPP(I,J))/(PLE(I,J,L+LB)-PLE(I,J,L+LB-1))))
             stratwgt = 1.0 - tropwgt    

             ! Count number of cells since vertical loop start that have been (at least partly) in stratosphere
             IF ( stratwgt > 0.1 ) StratCount = StratCount + 1 

             ! Skip cell if concentration change is too small
             IF ( iopt%ApplyIncrement .AND. ABS(AnaPtr(I,J,L)) < MinConc ) CYCLE 
             IF ( iopt%NonZeroIncOnly .AND. ABS(IncPtr(I,J,L)) < MinConc ) CYCLE

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

             ! Fraction must be between 0 and 1
             wgt = max(0.0,min(1.0,wgt))
             IF ( wgt == 0.0 ) CYCLE

             ! Get target concentration in v/v dry
             SpcAna = AnaPtr(I,J,L)
             IF ( UnitFlag == 1 ) SpcAna = SpcAna * ( MAPL_AIRMW / mwSpc ) / ( 1. - Q(I,J,L) )
             IF ( UnitFlag == 3 ) SpcAna = SpcAna * 1.0e-6
             IF ( UnitFlag == 4 ) SpcAna = SpcAna * 1.0e-9

             ! Update field
             SpcNew = SpcBkg(I,J,L)
             IF ( iopt%ApplyIncrement ) THEN
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

             ! Eventually update second species to maintain concentration ratio of species 2 / species 1
             UpdateSpec2 = .FALSE.
             IF ( iopt%HasSpec2 ) THEN
                ! Use background field if in stratosphere and no adjustment to be done in stratosphere
                IF     ( stratwgt >= 0.5 .AND. .NOT. iopt%Spec2Strat ) THEN
                   Spc2Asm(I,J,L) = Spc2Bkg(I,J,L)
                ! Use background field if in troposphere and no adjustment to be done in troposphere
                ELSEIF ( tropwgt  >= 0.5 .AND. .NOT. iopt%Spec2Trop  ) THEN
                   Spc2Asm(I,J,L) = Spc2Bkg(I,J,L)
                ! Calculate Spc2/Spc1 ratio before update and maintain that ratio in assimilation field
                ELSE    
                   UpdateSpec2 = .TRUE.
                   Spc2toSpc = Spc2Bkg(I,J,L) / MAX(SpcBkg(I,J,L),MinConc)
                   ! Restrict ratio to limits specified in settings
                   Spc2toSpc = MIN(iopt%Spec2MaxRatio,MAX(Spc2toSpc,iopt%Spec2MinRatio)) 
                   Spc2Asm(I,J,L) = SpcAsm(I,J,L) * Spc2toSpc
                ENDIF
             ENDIF
 
             ! Update diagnostics
             IF ( ASSOCIATED(DiagInc        ) ) DiagInc(I,J,L)      = SpcAsm(I,J,L) - SpcBkg(I,J,L)
             IF ( ASSOCIATED(DiagIncFrac    ) ) DiagIncFrac(I,J,L)  = SpcAsm(I,J,L) / MAX(SpcBkg(I,J,L),MinConc)
             IF( UpdateSpec2 ) THEN 
                IF ( ASSOCIATED(DiagInc2    ) ) DiagInc2(I,J,L)     = Spc2Asm(I,J,L) - Spc2Bkg(I,J,L)
                IF ( ASSOCIATED(DiagIncFrac2) ) DiagIncFrac2(I,J,L) = Spc2Asm(I,J,L) / MAX(Spc2Bkg(I,J,L),MinConc)
                IF ( ASSOCIATED(DiagSpcRatio) ) DiagSpcRatio(I,J,L) = Spc2toSpc
             ENDIF 
             IF ( ASSOCIATED(DiagMsk2d      ) ) DiagMsk2d(I,J)      = DiagMsk2d(I,J) + 1.0
             IF ( ASSOCIATED(DiagMsk3d      ) ) DiagMsk3d(I,J,L)    = 1.0 
          ENDDO
       ENDDO
       ENDDO

       ! Print warning if at least one negative cell
       IF ( NNEG > 0 ) THEN
          WRITE(*,*) '*** DoAnalysis_ warning: encountered concentration below threshold, set to minimum: ',TRIM(SpecName),NNEG,MinConc,' ***'
       ENDIF

       ! Pass back to State_Chm%Species array: flip vertical axis and convert v/v dry to kg/kg total
       ! -------------------------------------------------------------------------------------------
       State_Chm%Species(indSpc)%Conc(:,:,LM:1:-1)     = SpcAsm(:,:,:)  * (1.-Q) / MAPL_AIRMW * mwSpc
       ! testing only
       !State_Chm%Species(indSpc)%Conc(:,:,LM:1:-1) = SpcAsm(:,:,:)
       IF ( iopt%HasSpec2 ) THEN
          State_Chm%Species(indSpc2)%Conc(:,:,LM:1:-1) = Spc2Asm(:,:,:) * (1.-Q) / MAPL_AIRMW * mwSpc2
       ENDIF

    ENDIF

    ! Cleanup
    ! -------
    IF ( ALLOCATED(SpcBkg ) ) DEALLOCATE(SpcBkg)
    IF ( ALLOCATED(SpcAsm ) ) DEALLOCATE(SpcAsm)
    IF ( ALLOCATED(Spc2Bkg) ) DEALLOCATE(Spc2Bkg)
    IF ( ALLOCATED(Spc2Asm) ) DEALLOCATE(Spc2Asm)

    IF ( ASSOCIATED(DiagInc     ) ) DiagInc      => NULL()
    IF ( ASSOCIATED(DiagIncFrac ) ) DiagIncFrac  => NULL()
    IF ( ASSOCIATED(DiagInc2    ) ) DiagInc2     => NULL()
    IF ( ASSOCIATED(DiagIncFrac2) ) DiagIncFrac2 => NULL()
    IF ( ASSOCIATED(DiagSpcRatio) ) DiagSpcRatio => NULL()
    IF ( ASSOCIATED(DiagMsk2d   ) ) DiagMsk2d    => NULL()
    IF ( ASSOCIATED(DiagMsk3d   ) ) DiagMsk3d    => NULL()

    IF ( HasBundle ) THEN
       CALL MAPL_SimpleBundleDestroy ( VarBundle, __RC__ )
       IF ( iopt%UseObsHour ) CALL MAPL_SimpleBundleDestroy ( VarBundleH, __RC__ )
    ENDIF
    iopt => NULL()

    ! Successful return
    RETURN_(ESMF_SUCCESS) 

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
  SUBROUTINE GetAnaTime_( Clock, Fwd, yy, mm, dd, h, m, s, RC )
!
! !USES:
!
    USE TIME_MOD,  ONLY : GET_TS_CHEM
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ESMF_Clock),    INTENT(INOUT)         :: Clock      ! ESMF Clock object
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
    TYPE(ESMF_TIME)            :: currTime
    TYPE(ESMF_TimeInterval)    :: tsChemInt
    REAL                       :: tsChem

    ! Begins here
    __Iam__('GetAnaTime_')

    ! Get current time 
    CALL ESMF_ClockGet( Clock, currTime = currTime, __RC__ )

    ! Eventually adjust time
    IF ( Fwd ) THEN
       tsChem = GET_TS_CHEM()
    ELSE
       tsChem = 0.0
    ENDIF
    CALL ESMF_TimeIntervalSet(tsChemInt, s_r8=real(tsChem,8), __RC__ )
    CALL ESMF_TimeGet( currTime+tsChemInt, yy=yy, mm=mm, dd=dd, h=h, m=m, s=s, __RC__ )

    ! All done
    RETURN_(ESMF_SUCCESS)

    END SUBROUTINE GetAnaTime_
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Model                            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: GetAnaBundle_ 
!
! !DESCRIPTION: Get analysis data bundle. 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE GetAnaBundle_( am_I_Root, FileTmpl,  bName, yy, mm, dd, h, m, grid, &
                            VarBundle, HasBundle, ifile, fileTime, only_vars,   &
                            err_mode,  anatime,   RC )
!
! !USES:
!
!
! !INPUT/OUTPUT PARAMETERS:
!
    LOGICAL,             INTENT(IN)            :: am_I_Root  ! Root CPU?
    CHARACTER(LEN=*),    INTENT(IN)            :: FileTmpl   ! file template
    CHARACTER(LEN=*),    INTENT(IN)            :: bName      ! bundle name 
    INTEGER,             INTENT(IN)            :: yy, mm, dd ! year, month, day
    INTEGER,             INTENT(IN)            :: h,  m      ! hour, minute, second
    TYPE(ESMF_Grid),     INTENT(INOUT)         :: grid       ! output grid
    TYPE(MAPL_SimpleBundle)                    :: VarBundle  ! Bundle
    LOGICAL,             INTENT(INOUT)         :: HasBundle  ! Was bundle found?
    CHARACTER(LEN=*),    INTENT(OUT), OPTIONAL :: ifile      ! file name 
    TYPE(ESMF_TIME),     INTENT(OUT), OPTIONAL :: fileTime   ! file time
    CHARACTER(LEN=*),    INTENT(IN),  OPTIONAL :: only_vars  ! variables to read
    INTEGER,             INTENT(IN),  OPTIONAL :: err_mode   ! error mode
    LOGICAL,             INTENT(IN),  OPTIONAL :: anatime    ! round time to analysis time? 
    INTEGER,             INTENT(OUT), OPTIONAL :: RC         ! Success or failure?
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
    CHARACTER(LEN=ESMF_MAXSTR) :: ifile_
    CHARACTER(LEN=4)           :: syy
    CHARACTER(LEN=2)           :: smm, sdd, sh, sm
    INTEGER                    :: nymd, nhms, incSecs
    INTEGER                    :: yy_, mm_, dd_, h_, m_, s_, fid
    TYPE(ESMF_TIME)            :: currTime, fileTime_
    TYPE(ESMF_TimeInterval)    :: tsInt
    LOGICAL                    :: HasFile,  anatime_
    INTEGER                    :: errmode_

    ! Begins here
    __Iam__('GetAnaBundle_')

    ! Initialize
    HasBundle = .FALSE.
    errmode_ = 2
    anatime_ = .FALSE.
    if ( present(err_mode) ) errmode_ = err_mode
    if ( present(anatime ) ) anatime_ = anatime

    ! Get date & time of file. These are the passed values by default 
    yy_ = yy
    mm_ = mm
    dd_ = dd
    h_  = h
    m_  = m
    ! If anatime is true, set time to closest analysis hour (0z, 6z, 12z, 18z)
    IF ( anatime_ ) THEN
       IF (     h < 3  ) THEN
          h_ = 0
       ELSEIF ( h < 9  ) THEN
          h_ = 6
       ELSEIF ( h < 15 ) THEN
          h_ = 12
       ELSEIF ( h < 21 ) THEN
          h_ = 18
       ! If 21z, get next day (but keep minutes)
       ELSE
          call ESMF_TimeSet(currTime, yy=yy_, mm=mm_, dd=dd_, h=23, m=m_, s=0)
          call ESMF_TimeIntervalSet(tsInt, s_r8=real(7200.0,8), __RC__ )
          call ESMF_TimeGet( currTime+tsInt, yy=yy_, mm=mm_, dd=dd_)
          h_ = 0
       ENDIF
    ENDIF

    ! Parse file name
    ifile_ = FileTmpl
    write(syy,'(I4.4)') yy_
    CALL ReplaceChar_ ( ifile_, '%y4', syy )
    write(smm,'(I2.2)') mm_
    CALL ReplaceChar_ ( ifile_, '%m2', smm )
    write(sdd,'(I2.2)') dd_
    CALL ReplaceChar_ ( ifile_, '%d2', sdd )
    write(sh,'(I2.2)') h_
    CALL ReplaceChar_ ( ifile_, '%h2', sh  )
    write(sm,'(I2.2)') m_
    CALL ReplaceChar_ ( ifile_, '%n2', sm  )

    ! set default file time 
    s_ = 0
    call ESMF_TimeSet(fileTime_, yy=yy_, mm=mm_, dd=dd_, h=h_, m=m_, s=s_)

    ! Check if file exists
    INQUIRE( FILE=TRIM(ifile_), EXIST=HasFile )
    IF ( HasFile ) THEN
       IF ( am_I_Root ) WRITE(*,*) 'GCC GetAnaBundle_: Reading '//TRIM(ifile_)
       ! Try reading current time stamp on file 
       VarBundle =  MAPL_SimpleBundleRead ( TRIM(ifile_), TRIM(bname), grid, fileTime_, ONLY_VARS=only_vars, RC=STATUS )
       IF ( STATUS == ESMF_SUCCESS ) HasBundle = .TRUE.
       ! If current time stamp not found in file, just read the first entry (dangerous!)
       IF ( .NOT. HasBundle ) THEN
          ! If error mode is 0 or 1, stop with error
          IF ( errmode_ <= 1 ) THEN
             IF ( am_I_Root ) THEN
                WRITE(*,*) 'Error: current time not found in file: ',TRIM(ifile_),yy_,mm_,dd_,h_,m_
                WRITE(*,*) 'You can get past this error by setting the error mode to > 1'
             ENDIF
             ASSERT_(.FALSE.)
          ELSE
             IF ( am_I_Root ) THEN
                WRITE(*,*) 'Warning: current time not found in file - will read first time slice on file!! ',yy_,mm_,dd_,h_,m_
             ENDIF
             ! Get time stamp on file
             call GFIO_Open( ifile_, 1, fid, STATUS )
             ASSERT_(STATUS==0)
             call GetBegDateTime ( fid, nymd, nhms, incSecs, STATUS )
             ASSERT_(STATUS==0)
             caLL GFIO_Close( fid, STATUS )
             ASSERT_(STATUS==0)
             yy_ = nymd/10000
             mm_ = (nymd-yy_*10000) / 100
             dd_ = nymd - (10000*yy_ + mm_*100)
             h_  = nhms/10000
             m_  = (nhms- h_*10000) / 100
             s_  = nhms - (10000*h_  + m_*100)
             call ESMF_TimeSet(fileTime_, yy=yy_, mm=mm_, dd=dd_, h=h_, m=m_, s=s_)
             VarBundle =  MAPL_SimpleBundleRead ( TRIM(ifile_), TRIM(bname), grid, fileTime_, ONLY_VARS=only_vars, RC=STATUS )
             IF ( STATUS==ESMF_SUCCESS ) HasBundle = .TRUE.
          ENDIF
       ENDIF
    ! error handling if file not found 
    ELSE
       ! If file not found and error mode is zero, stop with error
       IF ( errmode_ == 0 ) THEN
          IF ( am_I_Root ) THEN
             WRITE(*,*) 'ERROR: file not found: '//TRIM(ifile_)
             WRITE(*,*) 'You can get past this error setting the error mode to > 0'
          ENDIF
          ASSERT_(.FALSE.)
       ! If file not found and error mode is not zero, just skip nudging 
       ELSE
          IF ( am_I_Root ) WRITE(*,*) '*** GCC warning in GetAnaBundle_, file not found: '//TRIM(ifile_)
       ENDIF
    ENDIF

    ! Return
    IF ( present(ifile   ) ) ifile    = ifile_
    IF ( present(fileTime) ) fileTime = fileTime_
    RETURN_(ESMF_SUCCESS)

  END SUBROUTINE GetAnaBundle_
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
  SUBROUTINE ReadSettings_( am_I_Root, GEOSCF, ispec, RC )
!
! !USES:
!
!
! !INPUT PARAMETERS:
!
    LOGICAL,              INTENT(IN)    :: am_I_Root  ! Root PET?
    INTEGER,              INTENT(IN)    :: ispec      ! species number 
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ESMF_CONFIG),    INTENT(INOUT) :: GEOSCF     ! GCC RC file 
!
! !OUTPUT PARAMETERS:
!
    INTEGER, INTENT(OUT)                :: RC         ! Success or failure
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
    TYPE(ESMF_Config)             :: CF
    CHARACTER(LEN=ESMF_MAXSTR)    :: ConfigNameLabel, ConfigName
    CHARACTER(LEN=3)              :: intStr
    INTEGER                       :: ThisInt 
    CHARACTER(LEN=ESMF_MAXSTR)    :: Iam
    INTEGER                       :: STATUS 

    !=======================================================================
    ! ReadSettings_ begins here 
    !=======================================================================
    Iam = 'ReadSettings_'

    ! Get name of configuration file with settings
    WRITE( intStr, '(I3.3)' ) ispec
    ConfigNameLabel = 'Analysis_Settings_Spec'//TRIM(intStr)//':'
    CALL ESMF_ConfigGetAttribute( GEOSCF, ConfigName, Label=TRIM(ConfigNameLabel), __RC__ )

    ! Load configuration file
    CF = ESMF_ConfigCreate (__RC__)
    IF ( am_I_Root ) write(*,*) 'Reading analysis settings from file '//TRIM(ConfigName)
    call ESMF_ConfigLoadFile (CF, TRIM(ConfigName), __RC__ ) 

    ! Read settings and write to configuration list
    CALL ESMF_ConfigGetAttribute( CF, AnaConfig(ispec)%SpecName,       Label='SpeciesName:'   ,                __RC__ ) 
    CALL ESMF_ConfigGetAttribute( CF, ThisInt,                         Label='Active:'        , Default=1,     __RC__ )
    AnaConfig(ispec)%Active = ( ThisInt == 1 )
    CALL ESMF_ConfigGetAttribute( CF, AnaConfig(ispec)%AnalysisFreq,   Label='AnalysisFreq:'  , Default=6,     __RC__ )
    CALL ESMF_ConfigGetAttribute( CF, AnaConfig(ispec)%AnalysisHour,   Label='AnalysisHour:'  , Default=0,     __RC__ )
    CALL ESMF_ConfigGetAttribute( CF, AnaConfig(ispec)%AnalysisMinute, Label='AnalysisMinute:', Default=0,     __RC__ )
    CALL ESMF_ConfigGetAttribute( CF, ThisInt,                         Label='ForwardLooking:', Default=1,     __RC__ )
    AnaConfig(ispec)%ForwardLooking = ( ThisInt == 1 )
    CALL ESMF_ConfigGetAttribute( CF, ThisInt,                         Label='ReadAnaTime:'   , Default=0,     __RC__ )
    AnaConfig(ispec)%ReadAnaTime = ( ThisInt == 1 )
    CALL ESMF_ConfigGetAttribute( CF, AnaConfig(ispec)%FileTemplate,   Label='FileTemplate:'  ,                __RC__ )
    CALL ESMF_ConfigGetAttribute( CF, AnaConfig(ispec)%FileVarName,    Label='FileVarName:'   ,                __RC__ )
    CALL ESMF_ConfigGetAttribute( CF, AnaConfig(ispec)%FileVarUnit,    Label='FileVarUnit:'   , Default='v/v', __RC__ )
    CALL ESMF_ConfigGetAttribute( CF, ThisInt,                         Label='ApplyIncrement:', Default=0,     __RC__ )
    AnaConfig(ispec)%ApplyIncrement = ( ThisInt == 1 )
    CALL ESMF_ConfigGetAttribute( CF, ThisInt,                         Label='InStrat:'       , Default=1,     __RC__ )
    AnaConfig(ispec)%InStrat = ( ThisInt == 1 )
    CALL ESMF_ConfigGetAttribute( CF, ThisInt,                         Label='InTrop:'        , Default=1,     __RC__ )
    AnaConfig(ispec)%InTrop = ( ThisInt == 1 )
    CALL ESMF_ConfigGetAttribute( CF, ThisInt,                         Label='NonZeroIncOnly:', Default=1,     __RC__ )
    AnaConfig(ispec)%NonZeroIncOnly = ( ThisInt == 1 )
    CALL ESMF_ConfigGetAttribute( CF, AnaConfig(ispec)%FileVarNameInc, Label='FileVarNameInc:',                __RC__ )
    CALL ESMF_ConfigGetAttribute( CF, AnaConfig(ispec)%AnaL1,          Label='AnaL1:'         , Default=1,     __RC__ )
    CALL ESMF_ConfigGetAttribute( CF, AnaConfig(ispec)%AnaL2,          Label='AnaL2:'         , Default=1,     __RC__ )
    CALL ESMF_ConfigGetAttribute( CF, AnaConfig(ispec)%AnaL3,          Label='AnaL3:'         , Default=72,    __RC__ )
    CALL ESMF_ConfigGetAttribute( CF, AnaConfig(ispec)%AnaL4,          Label='AnaL4:'         , Default=72,    __RC__ )
    CALL ESMF_ConfigGetAttribute( CF, AnaConfig(ispec)%AnaFraction,    Label='AnaFraction:'   , Default=1.0,   __RC__ )
    CALL ESMF_ConfigGetAttribute( CF, AnaConfig(ispec)%StratSponge,    Label='StratSponge:'   , Default=0,     __RC__ )
    CALL ESMF_ConfigGetAttribute( CF, AnaConfig(ispec)%MaxChangeStrat, Label='MaxChangeStrat:', Default=-1.0,  __RC__ )
    CALL ESMF_ConfigGetAttribute( CF, AnaConfig(ispec)%MaxChangeTrop , Label='MaxChangeTrop:' , Default=-1.0,  __RC__ )
    CALL ESMF_ConfigGetAttribute( CF, AnaConfig(ispec)%MaxRatioStrat , Label='MaxRatioStrat:' , Default=-1.0,  __RC__ )
    CALL ESMF_ConfigGetAttribute( CF, AnaConfig(ispec)%MaxRatioTrop  , Label='MaxRatioTrop:'  , Default=-1.0,  __RC__ )
    CALL ESMF_ConfigGetAttribute( CF, AnaConfig(ispec)%MinRatioStrat , Label='MinRatioStrat:' , Default=-1.0,  __RC__ )
    CALL ESMF_ConfigGetAttribute( CF, AnaConfig(ispec)%MinRatioTrop  , Label='MinRatioTrop:'  , Default=-1.0,  __RC__ )
    CALL ESMF_ConfigGetAttribute( CF, ThisInt,                         Label='UseObsHour:'    , Default=0,     __RC__ )
    AnaConfig(ispec)%UseObsHour = ( ThisInt == 1 )
    CALL ESMF_ConfigGetAttribute( CF, AnaConfig(ispec)%ObsHourName,    Label='ObsHourName:'   , Default='ana_hour', __RC__ )
    CALL ESMF_ConfigGetAttribute( CF, ThisInt,                         Label='HasSpec2:'      , Default=0,     __RC__ )
    AnaConfig(ispec)%HasSpec2 = ( ThisInt == 1 )
    CALL ESMF_ConfigGetAttribute( CF, AnaConfig(ispec)%Spec2Name,      Label='Spec2Name:'     , Default='N/A', __RC__ )
    CALL ESMF_ConfigGetAttribute( CF, ThisInt,                         Label='Spec2Strat:'    , Default=1,     __RC__ )
    AnaConfig(ispec)%Spec2Strat = ( ThisInt == 1 )
    CALL ESMF_ConfigGetAttribute( CF, ThisInt,                         Label='Spec2Trop:'     , Default=1,     __RC__ )
    AnaConfig(ispec)%Spec2Trop = ( ThisInt == 1 )
    CALL ESMF_ConfigGetAttribute( CF, AnaConfig(ispec)%Spec2MinRatio,  Label='Spec2MinRatio:' , Default=-1.0,  __RC__ )
    CALL ESMF_ConfigGetAttribute( CF, AnaConfig(ispec)%Spec2MaxRatio,  Label='Spec2MaxRatio:' , Default=-1.0,  __RC__ )
    CALL ESMF_ConfigGetAttribute( CF, AnaConfig(ispec)%MinConc,        Label='MinConc:'       , Default=1.0e-20, __RC__ )
    CALL ESMF_ConfigGetAttribute( CF, AnaConfig(ispec)%ErrorMode,      Label='ErrorMode:'     , Default=1   ,  __RC__ )

    ! Some logical checks
    IF ( AnaConfig(ispec)%ApplyIncrement ) AnaConfig(ispec)%NonZeroIncOnly = .FALSE.

    ! Verbose
    IF ( am_I_Root ) THEN
       WRITE(*,*) '----------------------------------------'
       WRITE(*,*) 'Analysis settings for GEOS-Chem species ',TRIM(AnaConfig(ispec)%SpecName),':'
       WRITE(*,*) 'Active: ',AnaConfig(ispec)%Active
       IF ( AnaConfig(ispec)%Active ) THEN
          WRITE(*,*) '- Analysis frequency            : ',AnaConfig(ispec)%AnalysisFreq
          WRITE(*,*) '- Analysis hour                 : ',AnaConfig(ispec)%AnalysisHour
          WRITE(*,*) '- Analysis minute               : ',AnaConfig(ispec)%AnalysisMinute
          WRITE(*,*) '- Forward looking file read     : ', AnaConfig(ispec)%ForwardLooking
          WRITE(*,*) '- Read file analysis time stamp : ', AnaConfig(ispec)%ReadAnaTime
          WRITE(*,*) '- Use observation hour          : ', AnaConfig(ispec)%UseObsHour
          WRITE(*,*) '- File template                 : ', TRIM(AnaConfig(ispec)%FileTemplate)
          WRITE(*,*) '- Variable name on file         : ', TRIM(AnaConfig(ispec)%FileVarName)
          WRITE(*,*) '- Variable unit on file         : ', TRIM(AnaConfig(ispec)%FileVarUnit)
          IF ( AnaConfig(ispec)%UseObsHour ) THEN
             WRITE(*,*) '- Observation hour name on file : ', TRIM(AnaConfig(ispec)%ObsHourName)
          ENDIF
          WRITE(*,*) '- Apply increments              : ', AnaConfig(ispec)%ApplyIncrement
          WRITE(*,*) '- Analysis where inc is not zero: ', AnaConfig(ispec)%NonZeroIncOnly
          IF ( AnaConfig(ispec)%NonZeroIncOnly ) THEN
             WRITE(*,*) '- Analysis inc variable name    : ', TRIM(AnaConfig(ispec)%FileVarNameInc)
          ENDIF
          WRITE(*,*) '- Apply analysis in stratosphere: ', AnaConfig(ispec)%InStrat
          WRITE(*,*) '- Apply analysis in troposphere : ', AnaConfig(ispec)%InTrop
          WRITE(*,*) '- Tropopause sponge layer       : ', AnaConfig(ispec)%StratSponge
          WRITE(*,*) '- Analysis level 1              : ', AnaConfig(ispec)%AnaL1  
          WRITE(*,*) '- Analysis level 2              : ', AnaConfig(ispec)%AnaL2
          WRITE(*,*) '- Analysis level 3              : ', AnaConfig(ispec)%AnaL3
          WRITE(*,*) '- Analysis level 4              : ', AnaConfig(ispec)%AnaL4
          WRITE(*,*) '- Analysis fraction             : ', AnaConfig(ispec)%AnaFraction
          WRITE(*,*) '- Max. absolute change in strat : ', AnaConfig(ispec)%MaxChangeStrat
          WRITE(*,*) '- Max. absolute change in trop  : ', AnaConfig(ispec)%MaxChangeTrop  
          WRITE(*,*) '- Max. relative change in strat : ', AnaConfig(ispec)%MaxRatioStrat
          WRITE(*,*) '- Max. relative change in trop  : ', AnaConfig(ispec)%MaxRatioTrop  
          WRITE(*,*) '- Min. relative change in strat : ', AnaConfig(ispec)%MinRatioStrat
          WRITE(*,*) '- Min. relative change in trop  : ', AnaConfig(ispec)%MinRatioTrop  
          WRITE(*,*) '- Min. concentration (for ratio): ', AnaConfig(ispec)%MinConc
          WRITE(*,*) '- Also update 2nd species       : ', AnaConfig(ispec)%HasSpec2
          IF ( AnaConfig(ispec)%HasSpec2 ) THEN
             WRITE(*,*) '- Name of 2nd species           : ', TRIM(AnaConfig(ispec)%Spec2Name)
             WRITE(*,*) '- Do 2nd species in stratosphere: ', AnaConfig(ispec)%Spec2Strat
             WRITE(*,*) '- Do 2nd species in troposphere : ', AnaConfig(ispec)%Spec2Trop 
             WRITE(*,*) '- Minimum 2nd/1st species ratio : ', AnaConfig(ispec)%Spec2MinRatio
             WRITE(*,*) '- Maximum 2nd/1st species ratio : ', AnaConfig(ispec)%Spec2MaxRatio
          ENDIF
          WRITE(*,*) '- Error mode                    : ', AnaConfig(ispec)%ErrorMode
       ENDIF ! Active 
    ENDIF

    ! Successful return
    RETURN_(ESMF_SUCCESS) 

  END SUBROUTINE ReadSettings_ 
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Model                            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ReplaceChar_ 
!
! !DESCRIPTION: Replaces all characters in a string 
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
END MODULE GEOS_Analysis 
