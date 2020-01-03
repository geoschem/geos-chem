!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: hcox_volcano_mod.F90
!
! !DESCRIPTION: Module HCOX\_Volcano\_Mod.F90 is a HEMCO extension to use
!  volcano emissions (such as AeroCom or OMI) from ascii tables. This module
!  reads the daily data tables and emits the emissions according to the
!  information in this file.
!\\
!\\
! !INTERFACE:
!
MODULE HCOX_Volcano_Mod
!
! !USES:
!
  USE HCO_Error_MOD
  USE HCO_Diagn_MOD
  USE HCOX_TOOLS_MOD
  USE HCOX_State_MOD, ONLY : Ext_State
  USE HCO_State_MOD,  ONLY : HCO_State

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC :: HCOX_Volcano_Run
  PUBLIC :: HCOX_Volcano_Init
  PUBLIC :: HCOX_Volcano_Final
!
! !PRIVATE MEMBER FUNCTIONS:
!
  PRIVATE :: ReadVolcTable
  PRIVATE :: EmitVolc
!
! !REMARKS:
! Each volcano table (e.g. AeroCom or OMI) is expected to list the volcano
! location, sulfur emissions (in kg S/s), and the volcano elevation as well
! as the volcano plume column height. These entries need be separated by space
! characters. For example:
!                                                                             .
! ###  LAT (-90,90), LON (-180,180), SULFUR [kg S/s], ELEVATION [m], CLOUD_COLUMN_HEIGHT [m]
! ### If elevation=cloud_column_height, emit in layer of elevation
! ### else, emit in top 1/3 of cloud_column_height
! volcano::
! 50.170 6.850 3.587963e-03 600. 600.
! ::
!                                                                             .
! The sulfur read from table is emitted as the species defined in the
! Volcano settings section. More than one species can be provided. Mass
! sulfur is automatically converted to mass of emitted species (using the
! emitted molecular weight and molecular ratio of the corresponding HEMCO
! species). Additional scale factors can be defined in the settings section
! by using the (optional) setting 'Scaling_<SpecName>'.
! For example, to emit SO2 and BrO from volcanoes, with an additional scale
! factor of 1e-4 kg BrO / kgS for BrO, use the following setting:
!115     Volcano           : on    SO2/BrO
!    --> Scaling_BrO       :       1.0e-4
!    --> Volcano_Source    :       OMI
!    --> Volcano_Table     :       $ROOT/VOLCANO/v2018-03/$YYYY/so2_volcanic_emissions_Carns.$YYYY$MM$DD.rc
!                                                                             .
! This extension was originally added for usage within GEOS-5 and AeroCom
! volcanic emissions, but has been modified to work with OMI-based volcanic
! emissions from Ge et al. (2016).
!                                                                             .
! When using this extension, you should turn off any other volcano emission
! inventories!
!                                                                             .
!  References:
!  ============================================================================
!  (1 ) Ge, C., J. Wang, S. Carn, K. Yang, P. Ginoux, and N. Krotkov,
!       Satellite-based global volcanic SO2 emissions and sulfate direct
!       radiative forcing during 2005-2012, J. Geophys. Res. Atmos., 121(7),
!       3446-3464, doi:10.1002/2015JD023134, 2016.
!
! !REVISION HISTORY:
!  04 Jun 2015 - C. Keller   - Initial version
!  28 Mar 2018 - M. Sulprizio- Update to allow for OMI-based volcanic emissions;
!                              Added Volcano_Source option to specify AeroCom
!                              or OMI emissions because the latter are only
!                              currently available for 2005-2012
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !MODULE VARIABLES:
!
  TYPE :: MyInst
   INTEGER                         :: Instance
   INTEGER                         :: ExtNr     = -1   ! Extension number
   INTEGER                         :: CatErupt  = -1   ! Category of eruptive emissions
   INTEGER                         :: CatDegas  = -1   ! Category of degassing emissions
   INTEGER                         :: nSpc      =  0   ! # of species
   INTEGER                         :: nVolc     =  0   ! # of volcanoes in buffer
   INTEGER,  ALLOCATABLE           :: SpcIDs(:)        ! HEMCO species IDs
   REAL(sp), ALLOCATABLE           :: SpcScl(:)        ! Species scale factors
   REAL(sp), ALLOCATABLE           :: VolcSlf(:)       ! Sulface emissions [kg S/s]
   REAL(sp), ALLOCATABLE           :: VolcElv(:)       ! Elevation [m]
   REAL(sp), ALLOCATABLE           :: VolcCld(:)       ! Cloud column height [m]
   INTEGER,  ALLOCATABLE           :: VolcIdx(:)       ! Lon grid index
   INTEGER,  ALLOCATABLE           :: VolcJdx(:)       ! Lat grid index
   INTEGER,  ALLOCATABLE           :: VolcBeg(:)       ! Begin time (optional)
   INTEGER,  ALLOCATABLE           :: VolcEnd(:)       ! End time   (optional)
   CHARACTER(LEN=255)              :: FileName         ! Volcano file name
   CHARACTER(LEN=255)              :: VolcSource       ! Volcano data source
   CHARACTER(LEN=61), ALLOCATABLE  :: SpcScalFldNme(:) ! Names of scale factor fields
   TYPE(MyInst), POINTER           :: NextInst => NULL()
  END TYPE MyInst

  ! Pointer to instances
  TYPE(MyInst), POINTER            :: AllInst => NULL()

  ! Volcano data is in kgS. Will be converted to kg emitted species.
  ! MW_S is the molecular weight of sulfur
  REAL(hp), PARAMETER             :: MW_S = 32.0_hp

CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCOX_Volcano_Run
!
! !DESCRIPTION: Subroutine HCOX\_Volcano\_Run is the driver routine
! for the customizable HEMCO extension.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCOX_Volcano_Run( ExtState, HcoState, RC )
!
! !USES:
!
    USE HCO_FluxArr_Mod,  ONLY : HCO_EmisAdd
!
! !INPUT PARAMETERS:
!
    TYPE(Ext_State), POINTER       :: ExtState    ! Module options
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(HCO_State), POINTER       :: HcoState    ! Hemco state
    INTEGER,         INTENT(INOUT) :: RC          ! Success or failure
!
! !REVISION HISTORY:
!  04 Jun 2015 - C. Keller   - Initial version
!  26 Oct 2016 - R. Yantosca - Don't nullify local ptrs in declaration stmts
!  20 Mar 2019 - R. Yantosca - Remove MinDiagLev=2 from calls to HCO_EmisAdd
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER               :: N
    LOGICAL               :: ERR

    ! Strings
    CHARACTER(LEN=255)    :: ErrMsg, ThisLoc

    ! Arrays
    REAL(sp)              :: SO2degas(HcoState%NX,HcoState%NY,HcoState%NZ)
    REAL(sp)              :: SO2erupt(HcoState%NX,HcoState%NY,HcoState%NZ)
    REAL(sp)              :: iFlx    (HcoState%NX,HcoState%NY,HcoState%NZ)

    ! Pointers
    TYPE(MyInst), POINTER :: Inst

    !=================================================================
    ! HCOX_VOLCANO_RUN begins here!
    !=================================================================

    ! Assume success
    RC = HCO_SUCCESS

    ! Sanity check: return if extension not turned on
    IF ( ExtState%Volcano <= 0 ) RETURN

    ! Enter
    CALL HCO_Enter( HcoState%Config%Err,                                     &
                    'HCOX_Volcano_Run (hcox_volcano_mod.F90)', RC           )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Define strings for error messgaes
    ErrMsg = ''
    ThisLoc =  &
    ' -> in HCOX_Volcano_Run (in module HEMCO/Extensions/hcox_volcano_mod.F90)'

    ! Get instance
    Inst => NULL()
    CALL InstGet( ExtState%Volcano, Inst, RC )
    IF ( RC /= HCO_SUCCESS ) THEN
       WRITE( ErrMsg, * ) 'Cannot find Volcano instance Nr. ', ExtState%Volcano
       CALL HCO_Error( HcoState%Config%Err, ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    !----------------------------------------------
    ! Read/update the volcano data
    ! (will be done only if this is a new day)
    !----------------------------------------------
    CALL ReadVolcTable( HcoState, ExtState, Inst, RC )
    IF ( RC /= HCO_SUCCESS ) THEN
       ErrMsg = 'Error encountered in "ReadVolcTable"!'
       CALL HCO_Error( HcoState%Config%Err, ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    !=======================================================================
    ! Compute volcano emissions for non dry-run simulations
    ! (Skip for GEOS-Chem dry-run or HEMCO-standalone dry-run)
    !=======================================================================
    IF ( .not. HcoState%Options%IsDryRun ) THEN

       ! Emit volcanos into SO2degas and SO2erupt arrays [kg S/m2/s]
       CALL EmitVolc( HcoState, ExtState, Inst, &
                      SO2degas, SO2erupt, RC    )
       IF ( RC /= HCO_SUCCESS ) THEN
          ErrMsg = 'Error encountered in "EmitVolc"!'
          CALL HCO_Error( HcoState%Config%Err, ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

       ! Add eruptive and degassing emissions to emission arrays & diagnostics
       DO N = 1, Inst%nSpc

          !-------------------------------------------
          ! Add degassing emissions
          !-------------------------------------------

          ! Convert from [kg S/m2/s] to [kg species/m2/s]
          iFlx = SO2degas * Inst%SpcScl(N)

          ! Apply user-defined scaling (if any) for this species
          CALL HCOX_Scale( HcoState, iFlx, &
                           TRIM(Inst%SpcScalFldNme(N)), RC )
          IF ( RC /= HCO_SUCCESS ) THEN
             ErrMsg = 'Error encountered in "HCOX_Scale (degassing)"!'
             CALL HCO_Error( HcoState%Config%Err, ErrMsg, RC, ThisLoc )
             RETURN
          ENDIF

          ! Add degassing emissions into the HEMCO state
          CALL HCO_EmisAdd( HcoState, iFlx, Inst%SpcIDs(N), &
                            RC, ExtNr=Inst%ExtNr, Cat=Inst%CatDegas )
          IF ( RC /= HCO_SUCCESS ) THEN
             ErrMsg = 'Error encountered in "HCO_EmisAdd" (degassing)!'
             CALL HCO_Error( HcoState%Config%Err, ErrMsg, RC, ThisLoc )
             RETURN
          ENDIF

          !-------------------------------------------
          ! Add eruptive emissions
          !-------------------------------------------

          ! Convert from [kg S/m2/s] to [kg species/m2/s]
          iFlx = SO2erupt * Inst%SpcScl(N)

          ! Apply user-defined scaling (if any) for this species
          CALL HCOX_Scale( HcoState, iFlx, &
                           TRIM(Inst%SpcScalFldNme(N)), RC )
          IF ( RC /= HCO_SUCCESS ) THEN
             ErrMsg = 'Error encountered in "HCOX_Scale" (eruptive"!'
             CALL HCO_Error( HcoState%Config%Err, ErrMsg, RC, ThisLoc )
             RETURN
          ENDIF

          ! Add eruptive emissions to the HEMCO state
          CALL HCO_EmisAdd( HcoState, iFlx, Inst%SpcIDs(N), &
                            RC, ExtNr=Inst%ExtNr, Cat=Inst%CatErupt )
          IF ( RC /= HCO_SUCCESS ) THEN
             ErrMsg = 'Error encountered in "HCO_EmisAdd" (eruptive)!'
             CALL HCO_Error( HcoState%Config%Err, ErrMsg, RC, ThisLoc )
             RETURN
          ENDIF

       ENDDO !N
    ENDIF

    !=======================================================================
    ! Exit
    !=======================================================================

    ! Cleanup
    Inst => NULL()

    ! Return w/ success
    CALL HCO_LEAVE( HcoState%Config%Err, RC )

  END SUBROUTINE HCOX_Volcano_Run
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCOX_Volcano_Init
!
! !DESCRIPTION: Subroutine HCOX\_Volcano\_Init initializes the HEMCO
! CUSTOM extension.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCOX_Volcano_Init( HcoState, ExtName,ExtState, RC )
!
! !USES:
!
    USE HCO_ExtList_Mod,    ONLY : GetExtNr
    USE HCO_ExtList_Mod,    ONLY : GetExtOpt
    USE HCO_ExtList_Mod,    ONLY : GetExtSpcVal
    USE HCO_STATE_MOD,      ONLY : HCO_GetExtHcoID
!
! !INPUT PARAMETERS:
!
    CHARACTER(LEN=*), INTENT(IN   ) :: ExtName    ! Extension name
    TYPE(Ext_State),  POINTER       :: ExtState   ! Module options
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(HCO_State),  POINTER       :: HcoState   ! Hemco state
    INTEGER,          INTENT(INOUT) :: RC

! !REVISION HISTORY:
!  04 Jun 2015 - C. Keller   - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    TYPE(MyInst), POINTER          :: Inst
    REAL(sp)                       :: ValSp
    INTEGER                        :: ExtNr, N, Dum
    LOGICAL                        :: FOUND
    CHARACTER(LEN=31), ALLOCATABLE :: SpcNames(:)
    CHARACTER(LEN=255)             :: MSG, Str

    !=================================================================
    ! HCOX_VOLCANO_INIT begins here!
    !=================================================================

    ! Extension Nr.
    ExtNr = GetExtNr( HcoState%Config%ExtList, TRIM(ExtName) )
    IF ( ExtNr <= 0 ) THEN
       MSG = 'The Volcano extension is turned off.'
       CALL HCO_MSG( HcoState%Config%Err,  MSG )
       RETURN
    ENDIF

    ! Enter
    CALL HCO_Enter( HcoState%Config%Err,                                     &
                    'HCOX_Volcano_Init (hcox_volcano_mod.F90)', RC          )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Create Volcano instance for this simulation
    Inst => NULL()
    CALL InstCreate( ExtNr, ExtState%Volcano, Inst, RC )
    IF ( RC /= HCO_SUCCESS ) THEN
       CALL HCO_Error( HcoState%Config%Err,                                  &
                      'Cannot create Volcano instance', RC                  )
       RETURN
    ENDIF

    ! Get species IDs.
    CALL HCO_GetExtHcoID( HcoState, ExtNr,     Inst%SpcIDs,                  &
                          SpcNames, Inst%nSpc, RC                           )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! There must be at least one species
    IF ( Inst%nSpc == 0 ) THEN
       CALL HCO_Error( HcoState%Config%Err,                                  &
                      'No Volcano species specified', RC                    )
       RETURN
    ENDIF

    ! Determine scale factor to be applied to each species. This is 1.00
    ! by default, but can be set in the HEMCO configuration file via setting
    ! Scaling_<SpcName>.
    CALL GetExtSpcVal( HcoState%Config, ExtNr,  Inst%nSpc,   SpcNames,       &
                       'Scaling',       1.0_sp, Inst%SpcScl, RC             )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Get species mask fields
    CALL GetExtSpcVal( HcoState%Config,    ExtNr,        Inst%nSpc,          &
                       SpcNames,           'ScaleField', HCOX_NOSCALE,       &
                       Inst%SpcScalFldNme, RC                               )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Add conversion factor from kg S to kg of emitted species
    DO N = 1, Inst%nSpc
       Inst%SpcScl(N) = Inst%SpcScl(N) * HcoState%Spc(Inst%SpcIDs(N))%EmMW_g &
                      * HcoState%Spc(Inst%SpcIDs(N))%MolecRatio / MW_S
    ENDDO

    ! Get location of volcano table. This must be provided.
    CALL GetExtOpt( HcoState%Config, ExtNr, 'Volcano_Table',                 &
                    OptValChar=Inst%FileName, FOUND=FOUND, RC=RC            )

    IF ( RC /= HCO_SUCCESS .OR. .NOT. FOUND ) THEN
       MSG = 'Cannot read Volcano table file name. Please provide '       // &
             'the Volcano table as a setting to the Volcano extension. '  // &
             'The name of this setting must be `Volcano_Table`.'
       CALL HCO_Error( HcoState%Config%Err, MSG, RC )
       RETURN
    ENDIF

    ! See if emissions data source is given
    ! As of v11-02f, options are AeroCom or OMI
    Inst%VolcSource = 'AeroCom'
    CALL GetExtOpt( HcoState%Config, ExtNr,       'Volcano_Source',          &
                    OptValChar=Str,  FOUND=FOUND, RC=RC                     )
    IF ( RC /= HCO_SUCCESS ) RETURN
    IF ( FOUND ) Inst%VolcSource = Str

    ! See if eruptive and degassing hierarchies are given
    Inst%CatErupt = 51
    Inst%CatDegas = 52
    CALL GetExtOpt( HcoState%Config, ExtNr,       'Cat_Degassing',           &
                    OptValInt=Dum,   FOUND=FOUND, RC=RC                     )
    IF ( RC /= HCO_SUCCESS ) RETURN
    IF ( FOUND ) Inst%CatDegas = Dum
    CALL GetExtOpt( HcoState%Config, ExtNr,       'Cat_Eruptive',            &
                    OptValInt=Dum,   FOUND=FOUND, RC=RC                    )
    IF ( RC /= HCO_SUCCESS ) RETURN
    IF ( FOUND ) Inst%CatErupt = Dum

    ! Verbose mode
    IF ( HcoState%amIRoot ) THEN
       MSG = 'Use emissions extension `Volcano`:'
       CALL HCO_MSG( HcoState%Config%Err,  MSG )

       MSG = ' - use the following species (Name, HcoID, Scaling relative to kgS):'
       CALL HCO_MSG( HcoState%Config%Err, MSG)
       DO N = 1, Inst%nSpc
          WRITE(MSG,*) TRIM(SpcNames(N)), ', ', Inst%SpcIDs(N), ', ', Inst%SpcScl(N)
          CALL HCO_MSG( HcoState%Config%Err, MSG)
          WRITE(MSG,*) 'Apply scale field: ', TRIM(Inst%SpcScalFldNme(N))
          CALL HCO_MSG( HcoState%Config%Err, MSG)
       ENDDO
       WRITE(MSG,*) ' - Emissions data source is ', TRIM(Inst%VolcSource)
       CALL HCO_MSG( HcoState%Config%Err,  MSG )
       WRITE(MSG,*) ' - Emit eruptive emissions as category ', Inst%CatErupt
       CALL HCO_MSG( HcoState%Config%Err,  MSG )
       WRITE(MSG,*) ' - Emit degassing emissions as category ', Inst%CatDegas
       CALL HCO_MSG( HcoState%Config%Err,  MSG )
    ENDIF

    ! Cleanup
    Inst => NULL()
    IF ( ALLOCATED(SpcNames) ) DEALLOCATE(SpcNames)

    CALL HCO_Leave( HcoState%Config%Err, RC )

  END SUBROUTINE HCOX_Volcano_Init
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCOX_Volcano_Final
!
! !DESCRIPTION: Subroutine HCOX\_AeroCom\_Final finalizes the HEMCO
!  AeroCom extension.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCOX_Volcano_Final( ExtState )
!
! !INPUT PARAMETERS:
!
    TYPE(Ext_State), POINTER :: ExtState   ! Module options
!
! !REVISION HISTORY:
!  04 Jun 2015 - C. Keller   - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
    !=================================================================
    ! HCOX_VOLCANO_FINAL begins here!
    !=================================================================
    CALL InstRemove( ExtState%Volcano )

  END SUBROUTINE HCOX_Volcano_Final
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ReadVolcTable
!
! !DESCRIPTION: Subroutine ReadVolcTable reads the AeroCom volcano table of the
!  current day.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE ReadVolcTable( HcoState, ExtState, Inst, RC )
!
! !USES:
!
    USE HCO_CharTools_Mod
    USE inquireMod,         ONLY : findfreeLun
    USE HCO_CLOCK_MOD,      ONLY : HcoClock_NewDay
    USE HCO_CLOCK_MOD,      ONLY : HcoClock_Get
    USE HCO_GeoTools_MOD,   ONLY : HCO_GetHorzIJIndex
    USE HCO_EXTLIST_MOD,    ONLY : HCO_GetOpt
!
! !INPUT PARAMETERS:
!
    TYPE(Ext_State),  POINTER       :: ExtState   ! Module options
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(HCO_State),  POINTER       :: HcoState   ! Hemco state
    TYPE(MyInst),     POINTER       :: Inst
    INTEGER,          INTENT(INOUT) :: RC

! !REVISION HISTORY:
!  04 Jun 2015 - C. Keller   - Initial version
!  28 Mar 2018 - M. Sulprizio- Add check for OMI-based emissions and use closest
!                              available year if simulation year is outside
!                              2005-2012
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER               :: YYYY, MM, DD
    INTEGER               :: N, LUN, IOS, AS
    INTEGER               :: nVolc, nCol
    REAL(sp)              :: Dum(10)
    REAL(hp), ALLOCATABLE :: VolcLon(:)      ! Volcano longitude [deg E]
    REAL(hp), ALLOCATABLE :: VolcLat(:)      ! Volcano latitude  [deg N]
    LOGICAL               :: FileExists, EOF
    CHARACTER(LEN=255)    :: ThisFile, ThisLine
    CHARACTER(LEN=255)    :: MSG,      FileMsg
    CHARACTER(LEN=255)    :: LOC = 'ReadVolcTable (hcox_volcano_mod.F90)'

    !=================================================================
    ! ReadVolcTable begins here!
    !=================================================================

    ! Do only if it's a new day...
    IF ( HcoClock_NewDay( HcoState%Clock, EmisTime=.TRUE. ) ) THEN

       ! Get current year, month, day
       CALL HcoClock_Get ( HcoState%Clock, cYYYY=YYYY, cMM=MM, cDD=DD, RC=RC )
       IF ( RC /= HCO_SUCCESS ) RETURN

#if defined( MODEL_GEOS )
       ! Error trap: skip leap days
       IF ( MM == 2 .AND. DD > 28 ) DD = 28
#endif

       ! Get file name
       ThisFile = Inst%FileName
       CALL HCO_CharParse( HcoState%Config, ThisFile, YYYY, MM, DD, 0, 0, RC )
       IF ( RC /= HCO_SUCCESS ) RETURN

       !--------------------------------------------------------------------
       ! In dry-run mode, print file path to dryrun log and exit.
       ! Otherwise, print file path to the HEMCO log file and continue.
       !--------------------------------------------------------------------

       ! Test if the file exists
       INQUIRE( FILE=TRIM( ThisFile ), EXIST=FileExists )

       ! Create a display string based on whether or not the file is found
       IF ( FileExists ) THEN
          FileMsg = 'HEMCO (VOLCANO): Opening'
       ELSE
          FileMsg = 'HEMCO (VOLCANO): REQUIRED FILE NOT FOUND'
       ENDIF

       ! Write file status to stdout and the HEMCO log
       IF ( Hcostate%amIRoot ) THEN
          WRITE( 6,   300 ) TRIM( FileMsg ), TRIM( ThisFile )
          WRITE( MSG, 300 ) TRIM( FileMsg ), TRIM( ThisFile )
          CALL HCO_MSG( HcoState%Config%Err, MSG )
 300      FORMAT( a, ' ', a )
       ENDIF

       ! For dry-run simulations, return to calling program.
       ! For regular simulations, throw an error if we can't find the file.
       IF ( HcoState%Options%IsDryRun ) THEN
          RETURN
       ELSE
          IF ( .not. FileExists ) THEN
             WRITE( MSG, 300 ) TRIM( FileMsg ), TRIM( ThisFile )
             CALL HCO_ERROR( HcoState%Config%Err, MSG, RC )
             RETURN
          ENDIF
       ENDIF

       !--------------------------------------------------------------------
       ! Read data from files
       !--------------------------------------------------------------------

       ! Open file
       LUN = findFreeLun()
       OPEN ( LUN, FILE=TRIM(ThisFile), STATUS='OLD', IOSTAT=IOS )
       IF ( IOS /= 0 ) THEN
          MSG = 'Error reading ' // TRIM(ThisFile)
          CALL HCO_ERROR( HcoState%Config%Err,  MSG, RC, THISLOC=LOC )
          RETURN
       ENDIF

       ! Get number of volcano records
       nVolc = 0
       DO
          CALL GetNextLine( LUN, ThisLine, EOF, RC )
          IF ( RC /= HCO_SUCCESS ) RETURN
          IF ( EOF ) EXIT

          ! Skip any entries that contain '::'
          IF ( INDEX( TRIM(ThisLine), '::') > 0 ) CYCLE

          ! If we make it to here, this is a valid entry
          nVolc = nVolc + 1
       ENDDO

       ! Verbose
       IF ( HCO_IsVerb(HcoState%Config%Err,2) ) THEN
          WRITE(MSG,*) 'Number of volcanoes: ', nVolc
          CALL HCO_MSG( HcoState%Config%Err, MSG)
       ENDIF

       ! Allocate arrays
       IF ( nVolc > 0 ) THEN
          ! Eventually deallocate previously allocated data
          IF ( ALLOCATED(Inst%VolcSlf) ) DEALLOCATE(Inst%VolcSlf)
          IF ( ALLOCATED(Inst%VolcElv) ) DEALLOCATE(Inst%VolcElv)
          IF ( ALLOCATED(Inst%VolcCld) ) DEALLOCATE(Inst%VolcCld)
          IF ( ALLOCATED(Inst%VolcIdx) ) DEALLOCATE(Inst%VolcIdx)
          IF ( ALLOCATED(Inst%VolcJdx) ) DEALLOCATE(Inst%VolcJdx)
          IF ( ALLOCATED(Inst%VolcBeg) ) DEALLOCATE(Inst%VolcBeg)
          IF ( ALLOCATED(Inst%VolcEnd) ) DEALLOCATE(Inst%VolcEnd)

          ALLOCATE(     VolcLon(nVolc), &
                        VolcLat(nVolc), &
                   Inst%VolcSlf(nVolc), &
                   Inst%VolcElv(nVolc), &
                   Inst%VolcCld(nVolc), &
                   Inst%VolcIdx(nVolc), &
                   Inst%VolcJdx(nVolc), &
                   Inst%VolcBeg(nVolc), &
                   Inst%VolcEnd(nVolc), &
                   STAT=AS )
          IF ( AS /= 0 ) THEN
             CALL HCO_ERROR ( HcoState%Config%Err, &
                              'Volc allocation error', RC, THISLOC=LOC )
             RETURN
          ENDIF
               VolcLon = 0.0_hp
               VolcLat = 0.0_hp
          Inst%VolcSlf = 0.0_sp
          Inst%VolcElv = 0.0_sp
          Inst%VolcCld = 0.0_sp
          Inst%VolcBeg = 0
          Inst%VolcEnd = 0

       ELSE
          WRITE(MSG,*) 'No volcano data found for year/mm/dd: ', YYYY, MM, DD
          CALL HCO_WARNING(HcoState%Config%Err,MSG,RC,WARNLEV=1,THISLOC=LOC)
       ENDIF

       ! Now read records
       IF ( nVolc > 0 ) THEN
          REWIND( LUN )

          N = 0
          DO
             CALL GetNextLine( LUN, ThisLine, EOF, RC )
             IF ( RC /= HCO_SUCCESS ) RETURN
             IF ( EOF ) EXIT

             ! Skip any entries that contain '::'
             IF ( INDEX( TRIM(ThisLine), '::') > 0 ) CYCLE

             ! Write this data into the following vector element
             N = N + 1
             IF ( N > nVolc ) THEN
                WRITE(MSG,*) 'N exceeds nVolc: ', N, nVolc, &
                             ' - This error occurred when reading ', &
                             TRIM(ThisFile), '. This line: ', TRIM(ThisLine)
                CALL HCO_ERROR ( HcoState%Config%Err, MSG, RC, THISLOC = LOC )
                RETURN
             ENDIF

             CALL HCO_CharSplit( TRIM(ThisLine), ' ', &
                                 HCO_GetOpt(HcoState%Config%ExtList,'Wildcard'), &
                                 Dum, nCol, RC )
             IF ( RC /= HCO_SUCCESS ) RETURN

             ! Allow for 5 or 7 values
             IF ( nCol /= 5 .and. nCol /= 7 ) THEN
                WRITE(MSG,*) 'Cannot parse line ', TRIM(ThisLine), &
                             'Expected five or seven entries, separated by ', &
                             'space character, instead found ', nCol
                CALL HCO_ERROR ( HcoState%Config%Err, MSG, RC, THISLOC = LOC )
                RETURN
             ENDIF

             ! Now pass to vectors
                  VolcLat(N) = Dum(1)
                  VolcLon(N) = Dum(2)
             Inst%VolcSlf(N) = Dum(3)
             Inst%VolcElv(N) = Dum(4)
             Inst%VolcCld(N) = Dum(5)

             ! Some lines also include start and end time
             IF ( nCol == 7 ) THEN
                Inst%VolcBeg(N) = Dum(6)
                Inst%VolcEnd(N) = DUM(7)
             ENDIF
          ENDDO

          ! At this point, we should have read exactly nVolc entries!
          IF ( N /= nVolc ) THEN
             WRITE(MSG,*) 'N /= nVolc: ', N, nVolc, &
                          ' - This error occurred when reading ', TRIM(ThisFile)
             CALL HCO_ERROR ( HcoState%Config%Err, MSG, RC, THISLOC = LOC )
             RETURN
          ENDIF

       ENDIF

       ! All done
       CLOSE ( LUN )

       ! Get grid box indeces for each location
       IF ( nVolc > 0 ) THEN
          CALL HCO_GetHorzIJIndex( HcoState, nVolc, VolcLon, &
                                   VolcLat, Inst%VolcIdx, Inst%VolcJdx, RC )
          IF ( RC /= HCO_SUCCESS ) RETURN
       ENDIF

       ! Save # of volcanoes in archive
       Inst%nVolc = nVolc

    ENDIF ! new day

    ! Cleanup
    IF ( ALLOCATED(VolcLon) ) DEALLOCATE(VolcLon)
    IF ( ALLOCATED(VolcLat) ) DEALLOCATE(VolcLat)

    ! Return w/ success
    RC = HCO_SUCCESS

  END SUBROUTINE ReadVolcTable
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: EmitVolc
!
! !DESCRIPTION: Subroutine EmitVolc reads the AeroCom volcano table of the
!  current day.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE EmitVolc( HcoState, ExtState, Inst, SO2d, SO2e, RC )
!
! !USES:
!
    USE HCO_CLOCK_MOD,      ONLY : HcoClock_Get
!
! !INPUT PARAMETERS:
!
    TYPE(Ext_State),  POINTER       :: ExtState   ! Module options
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(HCO_State),  POINTER       :: HcoState   ! Hemco state
    TYPE(MyInst),     POINTER       :: Inst
    INTEGER,          INTENT(INOUT) :: RC
!
! !OUTPUT PARAMETERS:
!
    REAL(sp),         INTENT(  OUT) :: SO2e(HcoState%NX,HcoState%NY,HcoState%NZ)
    REAL(sp),         INTENT(  OUT) :: SO2d(HcoState%NX,HcoState%NY,HcoState%NZ)
!
! !REVISION HISTORY:
!  04 Jun 2015 - C. Keller   - Initial version
!  20 Jun 2019 - M. Sulprizio- Update to accommodate beginning and start time
!                              entries for some volcanoes
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER            :: I, J, L, N, HH, MN, hhmmss
    LOGICAL            :: Erupt
    REAL(sp)           :: nSO2, zTop, zBot, PlumeHgt
    REAL(sp)           :: z1,   z2
    REAL(sp)           :: tmp1, tmp2, Frac
    REAL(sp)           :: totE, totD, volcE, volcD
    CHARACTER(LEN=255) :: MSG
    CHARACTER(LEN=255) :: LOC = 'EmitVolc (hcox_volcano_mod.F90)'

    !=================================================================
    ! EmitVolc begins here!
    !=================================================================

    ! Init
    SO2e = 0.0_sp
    SO2d = 0.0_sp
    totE = 0.0_sp
    totD = 0.0_sp

    ! Make sure all required grid quantities are defined
    IF ( .NOT. ASSOCIATED(HcoState%Grid%AREA_M2%Val) ) THEN
       CALL HCO_ERROR ( HcoState%Config%Err, &
                       'Grid box areas not defined', RC, THISLOC=LOC )
       RETURN
    ENDIF
    IF ( .NOT. ASSOCIATED(HcoState%Grid%ZSFC%Val) ) THEN
       CALL HCO_ERROR ( HcoState%Config%Err, &
                       'Surface heights not defined', RC, THISLOC=LOC )
       RETURN
    ENDIF
    IF ( .NOT. ASSOCIATED(HcoState%Grid%BXHEIGHT_M%Val) ) THEN
       CALL HCO_ERROR ( HcoState%Config%Err, &
                       'Grid box heights not defined', RC, THISLOC=LOC )
       RETURN
    ENDIF

    ! Get current hour, minute and save as hhmmss
    CALL HcoClock_Get ( HcoState%Clock, cH=HH, cM=MN, RC=RC )
    IF ( RC /= HCO_SUCCESS ) RETURN
    hhmmss = HH*10000 + MN*100

    ! Do for every volcano
    IF ( Inst%nVolc > 0 ) THEN
       DO N = 1, Inst%nVolc

          ! Grid box index for this volcano
          I = Inst%VolcIdx(N)
          J = Inst%VolcJdx(N)

          ! Skip if outside of domain
          IF( I < 1 .OR. J < 1 ) CYCLE

          ! Check if beginning and end time are specified
          ! Do not include emissions for this volcano outside of
          ! start and end times (mps, 6/20/19)
          IF ( Inst%VolcBeg(N) > 0 .or. Inst%VolcEnd(N) > 0 ) THEN
             IF ( hhmmss <  Inst%VolcBeg(N) ) CYCLE
             IF ( hhmmss >= Inst%VolcEnd(N) ) CYCLE
          ENDIF

          ! total emissions of this volcano
          volcE = 0.0_sp
          volcD = 0.0_sp

          z1 = HcoState%Grid%ZSFC%Val(I,J)

          ! Get total emitted kgS/m2/s. Data in table is in kgS/s.
          nSo2 = Inst%VolcSlf(N) / HcoState%Grid%AREA_M2%Val(I,J)

          ! Elevation of volcano base and volcano cloud top height [m]
          ! Make sure that the bottom / top are at least at surface level
          zBot = MAX(Inst%VolcElv(N),z1)
          zTop = MAX(Inst%VolcCld(N),z1)

          ! If volcano is eruptive, zBot /= zTop. In this case, evenly
          ! distribute emissions in top 1/3 of the plume
          IF ( zBot /= zTop ) THEN
             zBot  = zTop - ( ( zTop - zBot ) / 3.0_sp )
             Erupt = .TRUE.
          ELSE
             Erupt = .FALSE.
          ENDIF

          ! Volcano plume height
          PlumeHgt = zTop - zBot

          ! Distribute emissions into emission arrays. The volcano plume
          ! ranges from zBot to zTop.
          DO L = 1, HcoState%NZ

             ! Get top height of this box
             z2 = z1 + HcoState%Grid%BXHEIGHT_M%Val(I,J,L)

             ! Skip if the plume bottom is above this grid box top
             IF ( zBot >= z2 ) THEN
                z1 = z2
                CYCLE
             ENDIF

             ! If the plume top is below this grid box bottom, we can exit
             ! since there will be no more emissions to distribute.
             IF ( zTop < z1 ) EXIT

             ! If we make it to here, the volcano plume is at least partly
             ! within this level. Determine the fraction of the plume that
             ! is within heights z1 to z2.

             ! Get the bottom and top height of the plume within this layer.
             tmp1 = MAX(z1,zBot)  ! this layer's plume bottom
             tmp2 = MIN(z2,zTop)  ! this layer's plume top

             ! Special case that zTop is heigher than the highest level: make
             ! sure that all emissions are going to be used.
             IF ( ( L == HcoState%NZ ) .AND. ( zTop > z2 ) ) THEN
                tmp2 = zTop
             ENDIF

             ! Fraction of total plume that is within this layer
             IF ( PlumeHgt == 0.0_sp ) THEN
                Frac = 1.0_sp
             ELSE
                Frac = (tmp2-tmp1) / PlumeHgt
             ENDIF

             ! Distribute emissions
             IF ( Erupt ) THEN
                SO2e(I,J,L) = SO2e(I,J,L) + ( Frac * nSo2 )
                volcE       = volcE &
                            + ( Frac * nSo2 * HcoState%Grid%AREA_M2%Val(I,J) )
             ELSE
                SO2d(I,J,L) = SO2d(I,J,L) + ( Frac * nSo2 )
                volcD       = volcD &
                            + ( Frac * nSo2 * HcoState%Grid%AREA_M2%Val(I,J) )
             ENDIF

             ! The top height is the new bottom
             z1 = z2
          ENDDO

          ! testing
          !IF ( HCO_IsVerb(HcoState%Config%Err,3) ) THEN
          !   WRITE(MSG,*) 'Total eruptive  emissions of volcano ', N, ' [kgS/s]: ', volcE
          !   CALL HCO_MSG(HcoState%Config%Err,MSG)
          !   WRITE(MSG,*) 'Total degassing emissions of volcano ', N, ' [kgS/s]: ', volcD
          !   CALL HCO_MSG(HcoState%Config%Err,MSG)
          !ENDIF

          ! total
          totE = totE + volcE
          totD = totD + volcD

       ENDDO
    ENDIF

    ! verbose
    IF ( HCO_IsVerb(HcoState%Config%Err,3) ) THEN
       WRITE(MSG,*) 'Total eruptive  emissions [kgS/s]: ', totE
       CALL HCO_MSG(HcoState%Config%Err,MSG)
       WRITE(MSG,*) 'Total degassing emissions [kgS/s]: ', totD
       CALL HCO_MSG(HcoState%Config%Err,MSG)
    ENDIF

    ! Return w/ success
    RC = HCO_SUCCESS

  END SUBROUTINE EmitVolc
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
  SUBROUTINE InstGet( Instance, Inst, RC, PrevInst )
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
  SUBROUTINE InstCreate( ExtNr, Instance, Inst, RC )
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
  SUBROUTINE InstRemove( Instance )
!
! !INPUT PARAMETERS:
!
    INTEGER :: Instance
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
       IF ( ALLOCATED(Inst%VolcSlf      ) ) DEALLOCATE ( Inst%VolcSlf       )
       IF ( ALLOCATED(Inst%VolcElv      ) ) DEALLOCATE ( Inst%VolcElv       )
       IF ( ALLOCATED(Inst%VolcCld      ) ) DEALLOCATE ( Inst%VolcCld       )
       IF ( ALLOCATED(Inst%VolcIdx      ) ) DEALLOCATE ( Inst%VolcIdx       )
       IF ( ALLOCATED(Inst%VolcJdx      ) ) DEALLOCATE ( Inst%VolcJdx       )
       IF ( ALLOCATED(Inst%SpcIDs       ) ) DEALLOCATE ( Inst%SpcIDs        )
       IF ( ALLOCATED(Inst%SpcScl       ) ) DEALLOCATE ( Inst%SpcScl        )
       IF ( ALLOCATED(Inst%SpcScalFldNme) ) DEALLOCATE ( Inst%SpcScalFldNme )

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
END MODULE HCOX_Volcano_Mod
