!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: hco_state_mod.F90
!
! !DESCRIPTION: Module HCO\_State\_Mod contains definitions and sub-
! routines for the HEMCO state derived type. The HEMCO state object
! (HcoState) contains all information related to the HEMCO run, such
! as the HEMCO clock, information on the emission grid and the data
! fields to be read, details on all used species, various physical
! constants, etc.
! It also contains the final assembled 3D flux and 2D deposition
! arrays (to be passed to the overlaying model) and a pointer to the
! HEMCO configuration object (Config). The latter contains error and
! traceback information and holds the data fields (in the data list
! ConfigList).
!\\
!\\
! The HEMCO state object (typically called HcoState) for a given HEMCO
! run must be defined on the HEMCO-model interface level (subroutine
! HcoState\_Init).
!\\
!\\
! !INTERFACE:
!
MODULE HCO_State_Mod
!
! USES:
!
  USE HCO_Types_Mod
  USE HCO_Error_Mod
  USE HCO_Arr_Mod
  USE HCO_VertGrid_Mod

#if defined(ESMF_)
  USE ESMF
#endif

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC :: HcoState_Init
  PUBLIC :: HcoState_Final
  PUBLIC :: HCO_GetModSpcID
  PUBLIC :: HCO_GetHcoID
  PUBLIC :: HCO_GetExtHcoID

  !=========================================================================
  ! HCO_State: Main HEMCO State derived type
  !=========================================================================
  TYPE, PUBLIC :: HCO_State

     !%%%%% Species information %%%%%
     LOGICAL                     :: amIRoot    ! Is this the root CPU?

     !%%%%% Species information %%%%%
     INTEGER                     :: nSpc       ! # of species
     TYPE(HcoSpc),       POINTER :: Spc(:)     ! list of species

     !%%%%% Emission grid information %%%%%
     INTEGER                     :: NX         ! # of x-pts (lons) on this CPU
     INTEGER                     :: NY         ! # of y-pts (lats) on this CPU
     INTEGER                     :: NZ         ! # of z-pts (levs) on this CPU
     TYPE(HcoGrid),      POINTER :: Grid       ! HEMCO grid information
     TYPE(HcoClock),     POINTER :: Clock      ! HEMCO clock

     ! Data array
     TYPE(Arr3D_HP),     POINTER :: Buffer3D   ! Placeholder to store temporary
                                               ! 3D array.  Emissions will be
                                               ! written into this array if
                                               ! option FillBuffer = .TRUE.

     !%%%%% Constants and timesteps %%%%%
     TYPE(HcoPhys),      POINTER :: Phys       ! Physical constants
     REAL(sp)                    :: TS_EMIS    ! Emission timestep [s]
     REAL(sp)                    :: TS_CHEM    ! Chemical timestep [s]
     REAL(sp)                    :: TS_DYN     ! Dynamic  timestep [s]

     !%%%%% Aerosol quantities %%%%%
     INTEGER                     :: nDust      ! # of dust species
     LOGICAL                     :: MarinePOA  ! MUse marine organic aerosols?
     TYPE(HcoMicroPhys), POINTER :: MicroPhys  ! Microphysics settings

     !%%%%%  Run time options %%%%%
     TYPE(HcoOpt),       POINTER :: Options    ! HEMCO run options

     !%%%%% ReadLists %%%%%
     TYPE(RdList),      POINTER  :: ReadLists
     LOGICAL                     :: SetReadListCalled

     !%%%%% Emissions linked list %%%%%%
     TYPE(ListCont), POINTER     :: EmisList
     INTEGER                     :: nnEmisCont =  0 ! # of container in EmisList

     !%%%%% Data container indeces %%%%%
     ! Element i of cIDList will point to data-container with container
     ! ID i (e.g. cIDList(3) points to data-container with cID = 3).
     TYPE(cIDListPnt), POINTER   :: cIDList(:) => NULL()

     ! # of defined data containers. Will be automatically increased
     ! by one when creating a new data container (DataCont_Init)
     INTEGER                     :: nnDataCont = 0

     ! Define object based on TimeIdxCollection derived type
     TYPE(TimeIdxCollection), POINTER :: AlltIDx  => NULL()

     ! HEMCO configuration object
     TYPE(ConfigObj), POINTER    :: Config => NULL()

     ! Pointer to beginning of collections linked list
     TYPE(DiagnBundle),  POINTER :: Diagn  => NULL()

     !%%%%%  ESMF objects
#if defined(ESMF_)
     TYPE(ESMF_GridComp), POINTER :: GridComp
     TYPE(ESMF_State),    POINTER :: IMPORT
     TYPE(ESMF_State),    POINTER :: EXPORT
#endif
  END TYPE HCO_State
!
! !REVISION HISTORY:
!  20 Aug 2013 - C. Keller   - Initial version, adapted from
!                              gigc_state_chm_mod.F90
!  07 Jul 2014 - R. Yantosca - Cosmetic changes
!  30 Sep 2014 - R. Yantosca - Add HcoMicroPhys derived type to HcoState
!  08 Apr 2015 - C. Keller   - Added MaskFractions to HcoState options.
!  13 Jul 2015 - C. Keller   - Added option 'Field2Diagn'.
!  15 Feb 2016 - C. Keller   - Update to v2.0
!  02 Nov 2019 - H.P. Lin    - Add a HEMCO isDryRun option which is intended to flag
!                              that all "meaningful" IO is skipped and files should
!                              only be checked. If file does not exist DO NOT STOP
!                              THE RUN. (This is for GC Classic for now)
!EOP
!------------------------------------------------------------------------------
!BOC
CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HcoState_Init
!
! !DESCRIPTION: Routine HcoState\_Init initializes the HEMCO state object.
! This initializes (nullifies) all pointers and sets all HEMCO settings
! and options to default values.
! The here defined pointers are defined/connected at the HEMCO-model
! interface level.
! The passed HEMCO configuration object (HcoConfig) must be defined,
! e.g. this subroutine must be called after having read (at least
! stage 1 of) the HEMCO configuration file (Config\_ReadFile in
! hco\_config\_mod.F90).
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HcoState_Init( HcoState, HcoConfig, nSpecies, RC )
!
! !USES:
!
    USE HCO_EXTLIST_MOD,    ONLY : GetExtOpt, CoreNr
    USE HCO_UNIT_MOD,       ONLY : HCO_UnitTolerance
!
! !INPUT PARAMETERS:
!
    INTEGER,          INTENT(IN)    :: nSpecies  ! # HEMCO species
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(HCO_State),  POINTER       :: HcoState  ! HEMCO State object
    TYPE(ConfigObj),  POINTER       :: HcoConfig ! HEMCO Config object
    INTEGER,          INTENT(INOUT) :: RC        ! Return code
!
! !REVISION HISTORY:
!  20 Aug 2013 - C. Keller - Adapted from gigc_state_chm_mod.F90
!  07 Jan 2016 - E. Lundgren - Add physical constant RSTARG and updated
!                              Avgdr and g0 to NIST 2014 values
!  15 Feb 2016 - C. Keller - Now pass HcoConfig object
!  01 Nov 2016 - C. Keller - Now nullify all pointers
!  12 May 2017 - C. Keller - Added option ScaleEmis
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER            :: I, AS
    INTEGER            :: UnitTolerance
    LOGICAL            :: FOUND
    CHARACTER(LEN=255) :: MSG

    !=====================================================================
    ! HcoState_Init begins here!
    !=====================================================================

    ! For error handling
    CALL HCO_ENTER (HcoConfig%Err,'Init_HCO_State (hco_state_mod.F90)', RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    !=====================================================================
    ! Allocate emission field vectors
    !=====================================================================

    ! Check if already allocated
    IF ( ASSOCIATED(HcoState)) THEN
       CALL HCO_WARNING( HcoConfig%Err,'HcoState already allocated!', RC )
       RETURN
    ENDIF
    ALLOCATE ( HcoState )

    ! Is this the Root CPU?
    HcoState%amIRoot = HcoConfig%amIRoot

    ! Initialize vector w/ species information
    HcoState%nSpc = nSpecies
    IF ( nSpecies > 0 ) THEN
       ALLOCATE ( HcoState%Spc (nSpecies ), STAT=AS )
       IF ( AS /= 0 ) THEN
          CALL HCO_ERROR( HcoConfig%Err, 'Species', RC )
          RETURN
       ENDIF
    ENDIF

    ! Initalize species information. The effective values for species
    ! names, model IDs, etc. are set in the HEMCO-model interface
    ! routine.
    DO I = 1, nSpecies
       HcoState%Spc(I)%HcoID      = I
       HcoState%Spc(I)%ModID      = -1
       HcoState%Spc(I)%SpcName    = ''
       HcoState%Spc(I)%MW_g       = 0.0_dp
       HcoState%Spc(I)%EmMW_g     = 0.0_dp
       HcoState%Spc(I)%MolecRatio = 1.0_dp
       HcoState%Spc(I)%HenryK0    = 0.0_dp
       HcoState%Spc(I)%HenryCR    = 0.0_dp
       HcoState%Spc(I)%HenryPKA   = 0.0_dp

       ! Initialize data arrays. Pass dimension zero, which
       ! will just create a pointer to the data array (XX%Val).
       ! Will specify the arrays in HEMCO-model interface routine
       ! or when writing to them for the first time.
       CALL Hco_ArrInit( HcoState%Spc(I)%Emis, 0, 0, 0, RC )
       IF ( RC /= HCO_SUCCESS ) RETURN

       CALL Hco_ArrInit( HcoState%Spc(I)%Conc, 0, 0, 0, RC )
       IF ( RC /= HCO_SUCCESS ) RETURN

       CALL Hco_ArrInit( HcoState%Spc(I)%Depv, 0, 0, RC )
       IF ( RC /= HCO_SUCCESS ) RETURN
    ENDDO !I

    !=====================================================================
    ! Initialize grid
    !=====================================================================

    ! Initialize grid dimensions.
    HcoState%NX   = 0
    HcoState%NY   = 0
    HcoState%NZ   = 0
    ALLOCATE ( HcoState%Grid, STAT = AS )
    IF ( AS /= 0 ) THEN
       CALL HCO_ERROR( HcoConfig%Err, 'HEMCO grid', RC )
       RETURN
    ENDIF

    ! Initialize grid arrays.
    CALL HCO_ArrInit ( HcoState%Grid%XMID,       0, 0,    RC )
    IF ( RC /= HCO_SUCCESS ) RETURN
    CALL HCO_ArrInit ( HcoState%Grid%YMID,       0, 0,    RC )
    IF ( RC /= HCO_SUCCESS ) RETURN
    CALL HCO_ArrInit ( HcoState%Grid%XEDGE,      0, 0,    RC )
    IF ( RC /= HCO_SUCCESS ) RETURN
    CALL HCO_ArrInit ( HcoState%Grid%YEDGE,      0, 0,    RC )
    IF ( RC /= HCO_SUCCESS ) RETURN
    CALL HCO_ArrInit ( HcoState%Grid%PEDGE,      0, 0, 0, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN
    CALL HCO_ArrInit ( HcoState%Grid%YSIN,       0, 0,    RC )
    IF ( RC /= HCO_SUCCESS ) RETURN
    CALL HCO_ArrInit ( HcoState%Grid%AREA_M2,    0, 0,    RC )
    IF ( RC /= HCO_SUCCESS ) RETURN
    CALL HCO_ArrInit ( HcoState%Grid%PBLHEIGHT,  0, 0,    RC )
    IF ( RC /= HCO_SUCCESS ) RETURN
    CALL HCO_ArrInit ( HcoState%Grid%BXHEIGHT_M, 0, 0, 0, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN
    CALL HCO_ArrInit ( HcoState%Grid%ZSFC,       0, 0,    RC )
    IF ( RC /= HCO_SUCCESS ) RETURN
    CALL HCO_ArrInit ( HcoState%Grid%PSFC,       0, 0,    RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Initialize vertical grid
    HcoState%Grid%ZGRID => NULL()
    CALL HCO_VertGrid_Init( HcoState%Grid%ZGRID, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    !=====================================================================
    ! Set misc. parameter
    !=====================================================================

    ! Physical constants
    ALLOCATE ( HcoState%Phys, STAT = AS )
    IF ( AS /= 0 ) THEN
       CALL HCO_ERROR( HcoConfig%Err, 'HEMCO physical constants', RC )
       RETURN
    ENDIF
    HcoState%Phys%Avgdr  = 6.022140857e23_dp
    HcoState%Phys%PI     = 3.14159265358979323_dp
    HcoState%Phys%PI_180 = HcoState%Phys%PI / 180.0_dp
    HcoState%Phys%Re     = 6.375e6_dp
    HcoState%Phys%AIRMW  = 28.97_dp
    HcoState%Phys%g0     = 9.80665_dp
    HcoState%Phys%Rd     = 287.0_dp
    HcoState%Phys%Rdg0   = HcoState%Phys%Rd / HcoState%Phys%g0
    HcoState%Phys%RSTARG = 8.31450_dp

    ! Timesteps
    HcoState%TS_EMIS = 0.0_sp
    HcoState%TS_CHEM = 0.0_sp
    HcoState%TS_DYN  = 0.0_sp

    ! Nullify temporary array. This array may be used as temporary
    ! place to write emissions into.
    HcoState%Buffer3D => NULL()
    CALL HCO_ArrInit( HcoState%Buffer3D, 0, 0, 0, RC )
    IF ( RC /= 0 ) RETURN

    ! Dust bins (set default to 4)
    HcoState%nDust = 4

    ! Turn off marine POA by default
    HcoState%MarinePOA = .FALSE.

    ! Aerosol options
    ALLOCATE ( HcoState%MicroPhys, STAT = AS )
    IF ( AS /= 0 ) THEN
       CALL HCO_ERROR( HcoConfig%Err, 'HEMCO aerosol microphysics options', RC )
       RETURN
    ENDIF
    HcoState%MicroPhys%nBins           = 0
    HcoState%MicroPhys%nActiveModeBins = 0
    NULLIFY( HcoState%MicroPhys%BinBound )

    ! Default HEMCO options
    ! ==> execute HEMCO core; use all species, categories; not ESMF; not dryrun
    ALLOCATE( HcoState%Options )
    HcoState%Options%ExtNr          =  0
    HcoState%Options%SpcMin         =  1
    HcoState%Options%SpcMax         = -1
    HcoState%Options%CatMin         =  1
    HcoState%Options%CatMax         = -1
    HcoState%Options%AutoFillDiagn  = .TRUE.
    HcoState%Options%HcoWritesDiagn = .FALSE.
    HcoState%Options%FillBuffer     = .FALSE.
    HcoState%Options%isESMF         = .FALSE.
    HcoState%Options%isDryRun       = .FALSE.

    ! SetReadList has not been called yet
    HcoState%SetReadListCalled      = .FALSE.

    ! Get negative flag value from configuration file. If not found, set to 0.
    CALL GetExtOpt ( HcoConfig, CoreNr, 'Negative values', &
                     OptValInt=HcoState%Options%NegFlag, Found=Found, RC=RC )
    IF ( RC /= HCO_SUCCESS ) RETURN
    IF ( .NOT. Found ) HcoState%Options%NegFlag = 0

    ! Get PBL_DRYDEP flag from configuration file. If not found, set to default
    ! value of false.
    CALL GetExtOpt ( HcoConfig, CoreNr, 'PBL dry deposition', &
                     OptValBool=HcoState%Options%PBL_DRYDEP, Found=Found, RC=RC )
    IF ( RC /= HCO_SUCCESS ) RETURN
    IF ( .NOT. Found ) HcoState%Options%PBL_DRYDEP = .FALSE.

    ! Apply uniform scale factors specified in HEMCO_Config.rc?
    CALL GetExtOpt ( HcoConfig, CoreNr, 'Scale emissions', &
                     OptValBool=HcoState%Options%ScaleEmis, Found=Found, RC=RC )
    IF ( RC /= HCO_SUCCESS ) RETURN
    IF ( .NOT. Found ) HcoState%Options%ScaleEmis = .TRUE.

    ! Only shift hh/mm when applying time shift?
    CALL GetExtOpt ( HcoConfig, CoreNr, 'Cap time shift', &
                     OptValBool=HcoState%Options%TimeShiftCap, &
                     Found=Found, RC=RC )
    IF ( RC /= HCO_SUCCESS ) RETURN
    IF ( .NOT. Found ) HcoState%Options%TimeShiftCap = .FALSE.

    ! Get MaxDepExp from configuration file. If not found, set to default
    ! value of 20.
    CALL GetExtOpt ( HcoConfig, CoreNr, 'Maximum dep x ts', &
                     OptValHp=HcoState%Options%MaxDepExp, Found=Found, RC=RC )
    IF ( RC /= HCO_SUCCESS ) RETURN
    IF ( .NOT. Found ) HcoState%Options%MaxDepExp = 20.0_hp

    ! Get binary mask flag from configuration file. If not found, set to default
    ! value of TRUE.
    CALL GetExtOpt ( HcoConfig, CoreNr, 'Mask fractions', &
                     OptValBool=HcoState%Options%MaskFractions, Found=Found, RC=RC )
    IF ( RC /= HCO_SUCCESS ) RETURN
    IF ( .NOT. Found ) HcoState%Options%MaskFractions = .FALSE.

    CALL GetExtOpt ( HcoConfig, CoreNr, 'ConfigField to diagnostics', &
                     OptValBool=HcoState%Options%Field2Diagn, Found=Found, RC=RC )
    IF ( RC /= HCO_SUCCESS ) RETURN
    IF ( .NOT. Found ) HcoState%Options%Field2Diagn = .FALSE.

    CALL GetExtOpt ( HcoConfig, CoreNr, 'Vertical weights', &
                     OptValBool=HcoState%Options%VertWeight, Found=Found, RC=RC )
    IF ( RC /= HCO_SUCCESS ) RETURN
    IF ( .NOT. Found ) HcoState%Options%VertWeight = .TRUE.

    ! Make sure ESMF pointers are not dangling
#if defined(ESMF_)
    HcoState%GridComp => NULL()
    HcoState%IMPORT   => NULL()
    HcoState%EXPORT   => NULL()
#endif

    ! Read unit tolerance
    UnitTolerance = HCO_UnitTolerance( HcoConfig )

    ! Connect to config object
    HcoState%Config => HcoConfig

    ! Make sure pointers are not dangling
    HcoState%Diagn     => NULL()
    HcoState%EmisList  => NULL()
    HcoState%ReadLists => NULL()
    HcoState%Clock     => NULL()
    HcoState%cIDList   => NULL()
    HcoState%AlltIDx   => NULL()

    ! Verbose mode
    IF ( HCO_IsVerb(HcoConfig%Err,1) ) THEN
       WRITE(MSG,'(A68)') 'Initialized HEMCO state. Will use the following settings:'
       CALL HCO_MSG(HcoConfig%Err,MSG)
       WRITE(MSG,'(A33,I2)') 'Unit tolerance                 : ', UnitTolerance
       CALL HCO_MSG(HcoConfig%Err,MSG)
       WRITE(MSG,'(A33,I2)') 'Negative values                : ', HcoState%Options%NegFlag
       CALL HCO_MSG(HcoConfig%Err,MSG)
       WRITE(MSG,'(A33,L2)') 'Mask fractions                 : ', HcoState%Options%MaskFractions
       CALL HCO_MSG(HcoConfig%Err,MSG)
       WRITE(MSG,'(A33,L2)') 'Do drydep over entire PBL      : ', HcoState%Options%PBL_DRYDEP
       CALL HCO_MSG(HcoConfig%Err,MSG)
       WRITE(MSG,'(A33,F6.2)') 'Upper limit for deposition x ts: ', HcoState%Options%MaxDepExp
       CALL HCO_MSG(HcoConfig%Err,MSG)
       WRITE(MSG,'(A33,L2)') 'Scale emissions                : ', HcoState%Options%ScaleEmis
       CALL HCO_MSG(HcoConfig%Err,MSG)
       WRITE(MSG,'(A33,L2)') 'Cap time shift                 : ', HcoState%Options%TimeShiftCap
       CALL HCO_MSG(HcoConfig%Err,MSG)
    ENDIF

    ! Leave w/ success
    CALL HCO_LEAVE ( HcoConfig%Err, RC )

  END SUBROUTINE HcoState_Init
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HcoState_Final
!
! !DESCRIPTION: Routine HcoState\_CLEANUP cleans up HcoState.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HcoState_Final( HcoState )
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(HCO_State), POINTER  :: HcoState    ! HEMCO State object
!
! !REVISION HISTORY:
!  20 Aug 2013 - C. Keller   - Adapted from gigc_state_chm_mod.F90
!  24 Sep 2014 - R. Yantosca - Add an extra safety check when deallocating
!                              the pointer field HcoState%Spc
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER  :: I

    !=====================================================================
    ! HcoState_Final begins here!
    !=====================================================================

    ! Deallocate buffer array
    CALL HCO_ArrCleanup ( HcoState%Buffer3D )

    ! Deallocate all species arrays
    IF ( ASSOCIATED ( HcoState%Spc ) .and. HcoState%nSpc > 0 ) THEN
       DO I = 1, HcoState%nSpc
          CALL HCO_ArrCleanup( HcoState%Spc(I)%Emis )
          CALL HCO_ArrCleanup( HcoState%Spc(I)%Conc )
          CALL HCO_ArrCleanup( HcoState%Spc(I)%Depv )
       ENDDO
       DEALLOCATE( HcoState%Spc )
    ENDIF

    ! Deallocate grid information
    IF ( ASSOCIATED ( HcoState%Grid) ) THEN
       CALL HCO_VertGrid_Cleanup( HcoState%Grid%ZGRID )
       CALL HCO_ArrCleanup( HcoState%Grid%XMID        )
       CALL HCO_ArrCleanup( HcoState%Grid%YMID        )
       CALL HCO_ArrCleanup( HcoState%Grid%XEDGE       )
       CALL HCO_ArrCleanup( HcoState%Grid%YEDGE       )
       CALL HCO_ArrCleanup( HcoState%Grid%PEDGE       )
       CALL HCO_ArrCleanup( HcoState%Grid%YSIN        )
       CALL HCO_ArrCleanup( HcoState%Grid%AREA_M2     )
       CALL HCO_ArrCleanup( HcoState%Grid%PBLHEIGHT   )
       CALL HCO_ArrCleanup( HcoState%Grid%BXHEIGHT_M  )
       CALL HCO_ArrCleanup( HcoState%Grid%ZSFC        )
       CALL HCO_ArrCleanup( HcoState%Grid%PSFC        )
       DEALLOCATE(HcoState%Grid)
    ENDIF

    ! Deallocate microphysics information
    IF ( ASSOCIATED( HcoState%MicroPhys ) ) THEN
       IF ( HcoState%MicroPhys%nBins > 0 ) THEN
          IF ( ASSOCIATED( HcoState%MicroPhys%BinBound ) ) THEN
             NULLIFY( HcoState%MicroPhys%BinBound )
          ENDIF
          DEALLOCATE( HcoState%MicroPhys )
       ENDIF
    ENDIf

    ! Cleanup various types
    IF ( ASSOCIATED ( HcoState%Options ) ) DEALLOCATE ( HcoState%Options )
    IF ( ASSOCIATED ( HcoState%Phys    ) ) DEALLOCATE ( HcoState%Phys    )

#if defined(ESMF_)
    HcoState%GridComp => NULL()
    HcoState%IMPORT   => NULL()
    HcoState%EXPORT   => NULL()
#endif

  END SUBROUTINE HcoState_Final
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCO_GetModSpcId
!
! !DESCRIPTION: Function HCO\_GetModSpcId returns the model species index
! of a species by name. Returns -1 if given species is not found, 0 if
! name corresponds to the HEMCO wildcard character.
!\\
!\\
! !INTERFACE:
!
  FUNCTION HCO_GetModSpcID( name, HcoState ) RESULT( Indx )
!
! !USES:
!
      USE HCO_EXTLIST_MOD,     ONLY : HCO_GetOpt
!
! !INPUT PARAMETERS:
!
    CHARACTER(LEN=*), INTENT(IN)    :: name      ! Species name
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(HCO_State),  INTENT(INOUT) :: HcoState  ! HEMCO State
!
! !RETURN VALUE:
!
    INTEGER                         :: Indx      ! Index of this species
!
! !REVISION HISTORY:
!  20 Aug 2013 - C. Keller - Adapted from gigc_state_chm_mod.F90
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER :: N

    ! Default
    Indx = -1

    ! Return 0 if wildcard character
    IF ( TRIM(name) == TRIM(HCO_GetOpt(HcoState%Config%ExtList,'Wildcard')) ) THEN
       Indx = 0
       RETURN
    ENDIF

    ! Loop over all species names
    DO N = 1, HcoState%nSpc

       ! Return the index of the sought-for species
       IF( TRIM( name ) == TRIM( HcoState%Spc(N)%SpcName ) ) THEN
          Indx = HcoState%Spc(N)%ModID
          EXIT
       ENDIF

    ENDDO

  END FUNCTION HCO_GetModSpcID
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCO_GetHcoId
!
! !DESCRIPTION: Function HCO\_GetHcoIdHCO returns the HEMCO species index
! of a species by name. Returns -1 if given species is not found, 0 if
! name corresponds to the HEMCO wildcard character.
!\\
!\\
! !INTERFACE:
!
  FUNCTION HCO_GetHcoID( name, HcoState ) RESULT( Indx )
!
! !USES:
!
      USE HCO_EXTLIST_MOD,   ONLY : HCO_GetOpt
!
! !INPUT PARAMETERS:
!
    CHARACTER(LEN=*), INTENT(IN)   :: name         ! Species name
    TYPE(HCO_State), INTENT(INOUT) :: HcoState     ! HEMCO State
!
! !RETURN VALUE:
!
    INTEGER                      :: Indx         ! Index of this species
!
! !REVISION HISTORY:
!  20 Aug 2013 - C. Keller - Adapted from gigc_state_chm_mod.F90
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER :: N

    ! Default
    Indx = -1

    ! Return 0 if wildcard character
    IF ( TRIM(name) == TRIM(HCO_GetOpt(HcoState%Config%ExtList,'Wildcard')) ) THEN
       Indx = 0
       RETURN
    ENDIF

    ! Loop over all species names
    DO N = 1, HcoState%nSpc

       ! Return the index of the sought-for species
       IF( TRIM( name ) == TRIM( HcoState%Spc(N)%SpcName ) ) THEN
          Indx = N
          EXIT
       ENDIF
    ENDDO

  END FUNCTION HCO_GetHcoID
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: HCO_GetExtHcoID
!
! !DESCRIPTION: Subroutine HCO\_GetExtHcoID returns the HEMCO species IDs
! and names for all species assigned to the given extension (identified by
! its extension number).
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCO_GetExtHcoID( HcoState, ExtNr, HcoIDs, &
                              SpcNames, nSpc,  RC       )
!
! !USES:
!
    USE CHARPAK_MOD,         ONLY : STRSPLIT
    USE HCO_EXTLIST_MOD,     ONLY : GetExtSpcStr
    USE HCO_EXTLIST_MOD,     ONLY : HCO_GetOpt
!
! !INPUT PARAMETERS:
!
    TYPE(HCO_State),               POINTER       :: HcoState
    INTEGER,                       INTENT(IN   ) :: ExtNr       ! Extension #
!
! !OUTPUT PARAMETERS:
!
    INTEGER,          ALLOCATABLE, INTENT(  OUT) :: HcoIDs(:)   ! Species IDs
!
! !INPUT/OUTPUT PARAMETERS:
!
    CHARACTER(LEN=*), ALLOCATABLE, INTENT(INOUT) :: SpcNames(:) ! Species names
    INTEGER,                       INTENT(INOUT) :: nSpc        ! # of species
    INTEGER,                       INTENT(INOUT) :: RC          ! Success/fail
!
! !REVISION HISTORY:
!  10 Jan 2014 - C. Keller: Initialization (update)
!  29 Sep 2014 - C. Keller: Now allows species lists up to 2047 instead of 255
!                           characters.
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER             :: I,      AS
    CHARACTER(LEN=255)  :: MSG,    LOC
    CHARACTER(LEN=2047) :: SpcStr, SUBSTR(255)
    CHARACTER(LEN=2047) :: TmpStr

    !======================================================================
    ! HCO_GetExtHcoID begins here
    !======================================================================

    ! Enter
    LOC = 'HCO_GetExtHcoID (hco_state_mod.F90)'

    ! Get all species names belonging to extension Nr. ExtNr
    CALL GetExtSpcStr( HcoState%Config, ExtNr, SpcStr, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Split character into species string.
    CALL STRSPLIT( SpcStr, HCO_GetOpt(HcoState%Config%ExtList,'Separator'), SUBSTR, nSpc )

    ! nothing to do if there are no species
    IF ( nSpc == 0 ) RETURN

    ! Allocate arrays
    IF ( ALLOCATED(HcoIDs  ) ) DEALLOCATE(HcoIDs  )
    IF ( ALLOCATED(SpcNames) ) DEALLOCATE(SpcNames)
    ALLOCATE(HcoIDs(nSpc), SpcNames(nSpc), STAT=AS)
#if defined( MODEL_GEOS )
    SpcNames(:) = ''
    HcoIDs(:)   = -1
#endif
    IF ( AS/=0 ) THEN
       CALL HCO_ERROR(HcoState%Config%Err,'HcoIDs allocation error', RC, THISLOC=LOC)
       RETURN
    ENDIF

    ! Extract species information
    DO I = 1, nSpc
       !---------------------------------------------------------------------
       ! Prior to 6/26/18:
       ! This code can cause issues with certain compiler versions,
       ! so let's rewrite it slightly (bmy, 6/26/18)
       !SpcNames(I) = SUBSTR(I)
       !HcoIDs(I)   = HCO_GetHcoID( TRIM(SpcNames(I)), HcoState )
       !---------------------------------------------------------------------

       ! Rewrite this code to be a little more friendly to compilers with
       ! strict string-parsing syntax, such as ifort 17. ALSO NOTE: We don't
       ! necessarily have to do the TRIM in the call to HCO_GetHcoID, because
       ! the species name will be TRIMmed internally.  We have noticed that
       ! some compilers don't like taking the TRIM of an array element as
       ! an argument to a function call. (bmy, 6/26/18)
       TmpStr      = SubStr(I)
       SpcNames(I) = TRIM( TmpStr )
       HcoIDs(I)   = HCO_GetHcoID( TmpStr, HcoState )
    ENDDO

    ! Return w/ success
    RC = HCO_SUCCESS

  END SUBROUTINE HCO_GetExtHcoID
!EOC
END MODULE HCO_STATE_MOD
