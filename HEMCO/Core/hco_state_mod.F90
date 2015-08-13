!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: hco_state_mod.F90
!
! !DESCRIPTION: Module HCO\_State\_Mod contains definitions and sub- 
! routines for the HEMCO state derived type. The HEMCO state HcoState 
! contains information about the emissions grid, all used species, 
! various physical constants, etc. It also contains the final assembled 
! 3D flux and 2D deposition arrays (to be passed to the overlaying 
! model). 
! HcoState is defined on the HEMCO-model interface level.
!\\
!\\
! !INTERFACE: 
!
MODULE HCO_State_Mod
!
! USES:
!
  USE HCO_Error_Mod
  USE HCO_Arr_Mod

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
     INTEGER                     :: nSpc       ! # of species
     TYPE(HcoSpc),       POINTER :: Spc(:)     ! list of species

     !%%%%% Emission grid information %%%%%
     INTEGER                     :: NX         ! # of x-pts (lons) on this CPU
     INTEGER                     :: NY         ! # of y-pts (lats) on this CPU
     INTEGER                     :: NZ         ! # of z-pts (levs) on this CPU
     TYPE(HcoGrid),      POINTER :: Grid       ! HEMCO grid information
  
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
     TYPE(HcoMicroPhys), POINTER :: MicroPhys  ! Microphysics settings

     !%%%%%  Run time options %%%%%
     CHARACTER(LEN=255)          :: ConfigFile ! Full path to HEMCO Config file
     LOGICAL                     :: isESMF     ! Are we using ESMF?
     TYPE(HcoOpt),       POINTER :: Options    ! HEMCO run options

     !%%%%%  ESMF objects
#if defined(ESMF_)
     TYPE(ESMF_GridComp), POINTER :: GridComp 
     TYPE(ESMF_State),    POINTER :: IMPORT
     TYPE(ESMF_State),    POINTER :: EXPORT
#endif
  END TYPE HCO_State
!
! !PRIVATE TYPES:
!
  !=========================================================================
  ! HcoSpc: Derived type for HEMCO species
  !
  ! Notes:
  ! **1 The emission molecular weight is the molecular weight of the 
  !     emitted compound. This value is only different to MW_g if the 
  !     emitted compound does not correspond to the transported species, 
  !     e.g. if emissions are in kg C4H10 but the corresponding species 
  !     is transported as mass Carbon. 
  ! **2 MolecRatio is the ratio between # of species molecules per emitted 
  !       molecule, e.g. 4 if emissions are kg C4H10 but model species 
  !       are kg C.
  !=========================================================================
  TYPE :: HcoSpc
     INTEGER                 :: HcoID      ! HEMCO species ID
     INTEGER                 :: ModID      ! Model species ID
     CHARACTER(LEN= 31)      :: SpcName    ! species names
     REAL(hp)                :: MW_g       ! species molecular wt.     (g/mol)
     REAL(hp)                :: EmMW_g     ! emission molecular wt.**1 (g/mol)
     REAL(hp)                :: MolecRatio ! molecule emission ratio**2 (-)
     REAL(hp)                :: HenryK0    ! liq. over gas Henry const [M/atm]
     REAL(hp)                :: HenryCR    ! K0 temp. dependency [K] 
     REAL(hp)                :: HenryPKA   ! pKa for Henry const. correction
     TYPE(Arr2D_HP), POINTER :: Depv       ! Deposition velocity [1/s]
     TYPE(Arr3D_HP), POINTER :: Emis       ! Emission flux [kg/m2/s]
     TYPE(Arr3D_HP), POINTER :: Conc       ! Concentration [v/v]
  END TYPE HcoSpc

  !=========================================================================
  ! HcoOpt: Derived type for HEMCO run options
  !=========================================================================
  TYPE :: HcoOpt
     INTEGER  :: ExtNr          ! ExtNr to be used 
     INTEGER  :: SpcMin         ! Smallest HEMCO species ID to be considered 
     INTEGER  :: SpcMax         ! Highest HEMCO species ID to be considered
     INTEGER  :: CatMin         ! Smallest category to be considered
     INTEGER  :: CatMax         ! Highest category to be considered
     LOGICAL  :: HcoWritesDiagn ! If set to .TRUE., HEMCO will schedule the
                                ! output of the default HEMCO diagnostics 
                                ! (in hco_driver_mod.F90).
     LOGICAL  :: AutoFillDiagn  ! Write into AutoFill diagnostics?
     LOGICAL  :: FillBuffer     ! Write calculated emissions into buffer
                                ! instead of emission array? 
     INTEGER  :: NegFlag        ! Negative value flag (from configfile):
                                ! 2 = allow negative values
                                ! 1 = set neg. values to zero and prompt warning 
                                ! 0 = return w/ error if neg. value
     LOGICAL  :: PBL_DRYDEP     ! If true, dry deposition frequencies will
                                ! be calculated over the full PBL. If false, 
                                ! they are calculated over the first layer only.
     REAL(hp) :: MaxDepExp      ! Maximum value of deposition freq. x time step.
     LOGICAL  :: MaskFractions  ! If TRUE, masks are treated as binary, e.g.  
                                ! grid boxes are 100% inside or outside of a
                                ! mask. 
  END TYPE HcoOpt

  !=========================================================================
  ! HcoGrid: Derived type for HEMCO grid. The grid edges are used for data
  ! interpolation. The pressure midpoints are not needed by HEMCO core but
  ! can be specified for the extensions through ExtList.
  !
  ! NOTES:
  ! *  Not used in ESMF environment
  ! ** Only used by some extensions
  !=========================================================================
  TYPE :: HcoGrid
     TYPE(Arr2D_Hp), POINTER :: XMID       ! mid-points in x-direction (lon)
     TYPE(Arr2D_Hp), POINTER :: YMID       ! mid-points in y-direction (lat)
     TYPE(Arr2D_Hp), POINTER :: XEDGE      ! grid edges in x-direction (lon)*
     TYPE(Arr2D_Hp), POINTER :: YEDGE      ! grid edges in y-direction (lat)*
     TYPE(Arr3D_Hp), POINTER :: PEDGE      ! pressure edges (Pa) 
     TYPE(Arr2D_Hp), POINTER :: YSIN       ! sin of y-direction grid edges*
     TYPE(Arr2D_Hp), POINTER :: AREA_M2    ! grid box areas (m2)
     TYPE(Arr2D_Hp), POINTER :: ZSFC       ! surface geopotential height (m)**
     TYPE(Arr3D_Hp), POINTER :: BXHEIGHT_M ! grid box heights (m)** 
  END TYPE HcoGrid

  !=========================================================================
  ! HcoPhys: Derived type for HEMCO physical constants
  !=========================================================================
  TYPE :: HcoPhys
     REAL(dp) :: Avgdr   ! Avogadro number (mol-1)
     REAL(dp) :: PI      ! Pi
     REAL(dp) :: PI_180  ! Pi / 180
     REAL(dp) :: Re      ! Earth radius [m] 
     REAL(dp) :: AIRMW   ! Molecular weight of air (g/mol)
     REAL(dp) :: g0      ! Gravity at surface of earth (m/s2)
     REAL(dp) :: Rd      ! Gas Constant (R) in dry air (J/K/kg)
     REAL(dp) :: Rdg0    ! Rd/g0
  END TYPE HcoPhys 

  !=========================================================================
  ! HcoMicroPhys: Derived type for aerosol microphysics settings
  !=========================================================================
  TYPE :: HcoMicroPhys
     INTEGER           :: nBins              ! # of size-resolved bins
     INTEGER           :: nActiveModebins    ! # of active mode bins
     REAL(dp), POINTER :: BinBound(:)        ! Size bin boundaries
  END TYPE HcoMicroPhys
!                                                                             
! !REVISION HISTORY:
!  20 Aug 2013 - C. Keller   - Initial version, adapted from 
!                              gigc_state_chm_mod.F90
!  07 Jul 2014 - R. Yantosca - Cosmetic changes
!  30 Sep 2014 - R. Yantosca - Add HcoMicroPhys derived type to HcoState
!  08 Apr 2015 - C. Keller   - Added MaskFractions to HcoState options.
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
! and options to default values. The here defined pointers are connected
! at the HEMCO-model interface level. 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HcoState_Init( am_I_Root, HcoState, nSpecies, RC ) 
!
! !USES:
!
    USE HCO_EXTLIST_MOD,    ONLY : GetExtOpt, CoreNr
    USE HCO_UNIT_MOD,       ONLY : HCO_UnitTolerance
!
! !INPUT PARAMETERS:
! 
    LOGICAL,          INTENT(IN)    :: am_I_Root ! root CPU?
    INTEGER,          INTENT(IN)    :: nSpecies  ! # HEMCO species 
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(HCO_State),  POINTER       :: HcoState  ! HEMCO State object
    INTEGER,          INTENT(INOUT) :: RC        ! Return code
! 
! !REVISION HISTORY: 
!  20 Aug 2013 - C. Keller - Adapted from gigc_state_chm_mod.F90
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
    CALL HCO_ENTER ('Init_HCO_State (hco_state_mod.F90)', RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    !=====================================================================
    ! Allocate emission field vectors
    !=====================================================================

    ! Check if already allocated
    IF ( ASSOCIATED(HcoState)) THEN
       CALL HCO_WARNING( 'HcoState already allocated!', RC ) 
       RETURN
    ENDIF
    ALLOCATE ( HcoState )

    ! Initialize vector w/ species information
    HcoState%nSpc = nSpecies
    IF ( nSpecies > 0 ) THEN
       ALLOCATE ( HcoState%Spc (nSpecies ), STAT=AS )
       IF ( AS /= 0 ) THEN
          CALL HCO_ERROR( 'Species', RC )
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
       CALL HCO_ERROR( 'HEMCO grid', RC )
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
    CALL HCO_ArrInit ( HcoState%Grid%ZSFC,       0, 0,    RC )
    IF ( RC /= HCO_SUCCESS ) RETURN
    CALL HCO_ArrInit ( HcoState%Grid%BXHEIGHT_M, 0, 0, 0, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    !=====================================================================
    ! Set misc. parameter
    !=====================================================================

    HcoState%ConfigFile = '' 
    HcoState%isESMF     = .FALSE.

    ! Physical constants
    ALLOCATE ( HcoState%Phys, STAT = AS )
    IF ( AS /= 0 ) THEN
       CALL HCO_ERROR( 'HEMCO physical constants', RC )
       RETURN
    ENDIF
    HcoState%Phys%Avgdr  = 6.022e23_dp
    HcoState%Phys%PI     = 3.14159265358979323_dp
    HcoState%Phys%PI_180 = HcoState%Phys%PI / 180.0_dp 
    HcoState%Phys%Re     = 6.375e6_dp
    HcoState%Phys%AIRMW  = 28.97_dp
    HcoState%Phys%g0     = 9.8_dp
    HcoState%Phys%Rd     = 287.0_dp
    HcoState%Phys%Rdg0   = HcoState%Phys%Rd / HcoState%Phys%g0

    ! Timesteps
    HcoState%TS_EMIS = 0.0_sp 
    HcoState%TS_CHEM = 0.0_sp
    HcoState%TS_DYN  = 0.0_sp

    ! Nullify temporary array. This array may be used as temporary
    ! place to write emissions into.
    HcoState%Buffer3D => NULL()
    CALL HCO_ArrInit( HcoState%Buffer3D, 0, 0, 0, RC )
    IF ( RC /= 0 ) RETURN

    ! Aerosol options
    HcoState%nDust = 0
    ALLOCATE ( HcoState%MicroPhys, STAT = AS )
    IF ( AS /= 0 ) THEN
       CALL HCO_ERROR( 'HEMCO aerosol microphysics options', RC )
       RETURN
    ENDIF
    HcoState%MicroPhys%nBins           = 0
    HcoState%MicroPhys%nActiveModeBins = 0
    NULLIFY( HcoState%MicroPhys%BinBound )

    ! Default HEMCO options
    ! ==> execute HEMCO core; use all species and categories
    ALLOCATE( HcoState%Options )
    HcoState%Options%ExtNr          =  0
    HcoState%Options%SpcMin         =  1
    HcoState%Options%SpcMax         = -1
    HcoState%Options%CatMin         =  1
    HcoState%Options%CatMax         = -1
    HcoState%Options%AutoFillDiagn  = .TRUE.
    HcoState%Options%HcoWritesDiagn = .FALSE.
    HcoState%Options%FillBuffer     = .FALSE.

    ! Get negative flag value from configuration file. If not found, set to 0. 
    CALL GetExtOpt ( CoreNr, 'Negative values', &
                     OptValInt=HcoState%Options%NegFlag, Found=Found, RC=RC )
    IF ( RC /= HCO_SUCCESS ) RETURN
    IF ( .NOT. Found ) HcoState%Options%NegFlag = 0

    ! Get PBL_DRYDEP flag from configuration file. If not found, set to default
    ! value of false. 
    CALL GetExtOpt ( CoreNr, 'PBL dry deposition', &
                     OptValBool=HcoState%Options%PBL_DRYDEP, Found=Found, RC=RC )
    IF ( RC /= HCO_SUCCESS ) RETURN
    IF ( .NOT. Found ) HcoState%Options%PBL_DRYDEP = .FALSE. 

    ! Get MaxDepExp from configuration file. If not found, set to default
    ! value of 20. 
    CALL GetExtOpt ( CoreNr, 'Maximum dep x ts', &
                     OptValHp=HcoState%Options%MaxDepExp, Found=Found, RC=RC )
    IF ( RC /= HCO_SUCCESS ) RETURN
    IF ( .NOT. Found ) HcoState%Options%MaxDepExp = 20.0_hp 

    ! Get binary mask flag from configuration file. If not found, set to default
    ! value of TRUE. 
    CALL GetExtOpt ( CoreNr, 'Mask fractions', &
                     OptValBool=HcoState%Options%MaskFractions, Found=Found, RC=RC )
    IF ( RC /= HCO_SUCCESS ) RETURN
    IF ( .NOT. Found ) HcoState%Options%MaskFractions = .FALSE.

    ! Make sure ESMF pointers are not dangling 
#if defined(ESMF_)
    HcoState%GridComp => NULL()
    HcoState%IMPORT   => NULL()
    HcoState%EXPORT   => NULL()
#endif

    ! Read unit tolerance
    UnitTolerance = HCO_UnitTolerance()

    ! Verbose mode 
    IF ( HCO_IsVerb(1) ) THEN
       WRITE(MSG,'(A68)') 'Initialized HEMCO state. Will use the following settings:'
       CALL HCO_MSG(MSG)
       WRITE(MSG,'(A33,I2)') 'Unit tolerance                 : ', UnitTolerance 
       CALL HCO_MSG(MSG)
       WRITE(MSG,'(A33,I2)') 'Negative values                : ', HcoState%Options%NegFlag 
       CALL HCO_MSG(MSG)
       WRITE(MSG,'(A33,L2)') 'Mask fractions                 : ', HcoState%Options%MaskFractions
       CALL HCO_MSG(MSG)
       WRITE(MSG,'(A33,L2)') 'Do drydep over entire PBL      : ', HcoState%Options%PBL_DRYDEP
       CALL HCO_MSG(MSG)
       WRITE(MSG,'(A33,F6.2)') 'Upper limit for deposition x ts: ', HcoState%Options%MaxDepExp
       CALL HCO_MSG(MSG,SEP2='-')
    ENDIF

    ! Leave w/ success
    CALL HCO_LEAVE ( RC ) 

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
       CALL HCO_ArrCleanup( HcoState%Grid%XMID       )
       CALL HCO_ArrCleanup( HcoState%Grid%YMID       )
       CALL HCO_ArrCleanup( HcoState%Grid%XEDGE      )
       CALL HCO_ArrCleanup( HcoState%Grid%YEDGE      )
       CALL HCO_ArrCleanup( HcoState%Grid%PEDGE      )
       CALL HCO_ArrCleanup( HcoState%Grid%YSIN       )
       CALL HCO_ArrCleanup( HcoState%Grid%AREA_M2    )
       CALL HCO_ArrCleanup( HcoState%Grid%ZSFC       )
       CALL HCO_ArrCleanup( HcoState%Grid%BXHEIGHT_M )
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
      USE HCO_CHARTOOLS_MOD,   ONLY : HCO_WCD
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
    IF ( TRIM(name) == HCO_WCD() ) THEN
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
      USE HCO_CHARTOOLS_MOD,   ONLY : HCO_WCD
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
    IF ( TRIM(name) == HCO_WCD() ) THEN
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
    USE HCO_CHARTOOLS_MOD,   ONLY : HCO_SEP
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

    !======================================================================
    ! HCO_GetExtHcoID begins here
    !======================================================================

    ! Enter
    LOC = 'HCO_GetExtHcoID (hco_state_mod.F90)'

    ! Get all species names belonging to extension Nr. ExtNr
    CALL GetExtSpcStr( ExtNr, SpcStr, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Split character into species string. 
    CALL STRSPLIT( SpcStr, HCO_SEP(), SUBSTR, nSpc )

    ! nothing to do if there are no species
    IF ( nSpc == 0 ) RETURN 

    ! Allocate arrays 
    IF ( ALLOCATED(HcoIDs  ) ) DEALLOCATE(HcoIDs  ) 
    IF ( ALLOCATED(SpcNames) ) DEALLOCATE(SpcNames) 
    ALLOCATE(HcoIDs(nSpc), SpcNames(nSpc), STAT=AS)
    IF ( AS/=0 ) THEN
       CALL HCO_ERROR('HcoIDs allocation error', RC, THISLOC=LOC)
       RETURN
    ENDIF

    ! Extract species information
    DO I = 1, nSpc
       SpcNames(I) = SUBSTR(I)
       HcoIDs(I)   = HCO_GetHcoID( TRIM(SpcNames(I)), HcoState )
    ENDDO

    ! Return w/ success
    RC = HCO_SUCCESS 

  END SUBROUTINE HCO_GetExtHcoID 
!EOC
END MODULE HCO_STATE_MOD
