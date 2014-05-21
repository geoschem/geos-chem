!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: hco_state_mod
!
! !DESCRIPTION: Module HCO\_STATE\_MOD contains definitions and sub- 
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
MODULE HCO_STATE_MOD
!
! USES:
!
  USE HCO_ERROR_MOD
  USE HCO_ARR_MOD

#if defined(ESMF_)
  USE ESMF_Mod
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

  !=========================================================================
  ! HEMCO State derived type
  !=========================================================================
  TYPE, PUBLIC :: HCO_State

     ! Species information
     INTEGER                     :: nSpc         ! # of species
     TYPE(HcoSpc),       POINTER :: Spc(:)       ! list of species

     !
     ! Emission grid information 
     INTEGER                     :: NX     ! # of x-points (lons) on this CPU
     INTEGER                     :: NY     ! # of y-points (lats) on this CPU
     INTEGER                     :: NZ     ! # of z-points (levs) on this CPU
     TYPE(HcoGrid),     POINTER  :: Grid   ! HEMCO grid information

     ! Placeholder to store temporary 3D array
     ! Emissions will be written into this array only if option FillBuffer is enabled
     TYPE(Arr3D_HP),    POINTER  :: Buffer3D 

     ! Physical constants
     TYPE(HcoPhys),     POINTER  :: Phys 

     ! Emission & dynamic time step (seconds)
     REAL(sp)                   :: TS_EMIS
     REAL(sp)                   :: TS_CHEM
     REAL(sp)                   :: TS_DYN

     ! Settings
     CHARACTER(LEN=255)          :: ConfigFile    ! Path + Filename of configuration file 
     LOGICAL                     :: isESMF        ! ESMF application? 

     ! Run options
     TYPE(HcoOpt),     POINTER   :: Options

     ! If run w/in ESMF, also need to point to IMPORT state
     ! so that data arrays can be imported.
#if defined(ESMF_)
     TYPE(ESMF_State), POINTER   :: IMPORT
#endif

  END TYPE HCO_State
!
! !PRIVATE MODULE TYPES:
!

  ! HCO species
  TYPE :: HcoSpc
     INTEGER                     :: HcoID      ! HEMCO species ID
     INTEGER                     :: ModID      ! Model species ID
     CHARACTER(LEN= 31)          :: SpcName    ! species names
     REAL(hp)                    :: MW_g       ! species molecular weight (g/mol)
     REAL(hp)                    :: EmMW_g     ! emission molecular weight**1 (g/mol)
     REAL(hp)                    :: MolecRatio ! molecule emission ratio **2 (-)
     REAL(hp)                    :: HenryK0    ! liq. over gas Henry const [M/atm]
     REAL(hp)                    :: HenryCR    ! K0 temp. dependency [K] 
     REAL(hp)                    :: HenryPKA   ! pKa for Henry const. correction
     TYPE(Arr2D_HP),     POINTER :: Depv       ! Deposition velocity [m/s]
     TYPE(Arr3D_HP),     POINTER :: Emis       ! Emission flux [kg/m2/s]
  END TYPE HcoSpc
  ! Notes:
  ! **1 The emission molecular weight is the molecular weight of the emitted 
  !     compound. This value is only different to MW_g if the emitted compound
  !     does not correspond to the transported species, e.g. if emissions are 
  !     in kg C4H10 but the corresponding species is transported as mass Carbon. 
  ! **2 MolecRatio is the ratio between # of species molecules per emitted molecule, 
  !     e.g. 4 if emissions are kg C4H10 but model species are kg C.

  ! HEMCO run options
  TYPE :: HcoOpt
     INTEGER                    :: ExtNr         ! ExtNr to be used 
     INTEGER                    :: SpcMin        ! Smallest HEMCO species ID to be considered 
     INTEGER                    :: SpcMax        ! Highest HEMCO species ID to be considered
     INTEGER                    :: CatMin        ! Smallest category to be considered
     INTEGER                    :: CatMax        ! Highest category to be considered
     LOGICAL                    :: AutoFillDiagn ! Write into AutoFill diagnostics?
     LOGICAL                    :: FillBuffer    ! Write calculated emissions into buffer
                                                 ! instead of emission array? 
  END TYPE HcoOpt

  ! HEMCO grid
  TYPE :: HcoGrid
     REAL(df),           POINTER :: XMID       (:,:)   ! mid-points in x-direction (lon)
     REAL(df),           POINTER :: YMID       (:,:)   ! mid-points in y-direction (lat)
     REAL(df),           POINTER :: XEDGE      (:,:)   ! grid edges in x-direction (lon)*
     REAL(df),           POINTER :: YEDGE      (:,:)   ! grid edges in y-direction (lat)*
     REAL(df),           POINTER :: YSIN       (:,:)   ! sin of grid edges in y-direction (lat)*
     REAL(df),           POINTER :: AREA_M2    (:,:)   ! grid box areas (m2)
     REAL(df),           POINTER :: BXHEIGHT_M (:,:,:) ! grid box heights (m)**
  END TYPE HcoGrid
  ! *  Not used in ESMF environment
  ! ** Only used by some extensions

! Physical constants 
  TYPE :: HcoPhys
     REAL(dp)                    :: Avgdr   ! Avogadro number (mol-1)
     REAL(dp)                    :: PI      ! Pi
     REAL(dp)                    :: Re      ! Earth radius [m] 
     REAL(dp)                    :: AIRMW   ! Molecular weight of air (g/mol)
     REAL(dp)                    :: g0      ! Gravity at surface of earth (m/s2)
     REAL(dp)                    :: Rd      ! Gas Constant (R) in dry air (J/K/kg)
     REAL(dp)                    :: Rdg0    ! Rd/g0
  END TYPE HcoPhys 

!                                                                             
! !REVISION HISTORY:
!  20 Aug 2013 - C. Keller - Initial version, adapted from gigc_state_chm_mod.F90
!EOP
!------------------------------------------------------------------------------
!BOC
CONTAINS
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
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
  SUBROUTINE HcoState_Init ( am_I_Root, HcoState, nSpecies, RC ) 
!
! !USES:
!
!
! !PARAMETERS:
! 
    LOGICAL,          INTENT(IN)    :: am_I_Root ! root CPU?
    INTEGER,          INTENT(IN)    :: nSpecies  ! # HEMCO species 
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
    INTEGER                :: I, AS

    !=====================================================================
    ! HcoState_Init begins here!
    !=====================================================================

    ! For error handling
    CALL HCO_ENTER ('Init_HCO_State (HCO_STATE_MOD.F90)', RC )
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
    ALLOCATE ( HcoState%Spc (nSpecies ), STAT=AS )
    IF ( AS /= 0 ) THEN
       CALL HCO_ERROR( 'Species', RC )
       RETURN
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

       ! Initialize emission arrays. Pass dimension zero, which
       ! will just create a pointer to the emission/depostion
       ! array (Emis%Val, Depv%Val). Will specify the arrays in 
       ! HEMCO-model interface routine.
       HcoState%Spc(I)%Emis => NULL()
       CALL Hco_ArrInit( HcoState%Spc(I)%Emis, 0, 0, 0, RC )
       IF ( RC /= HCO_SUCCESS ) RETURN

       HcoState%Spc(I)%Depv => NULL()
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

    ! Nullify grid arrays. 
    NULLIFY ( HcoState%Grid%XMID       )
    NULLIFY ( HcoState%Grid%YMID       )
    NULLIFY ( HcoState%Grid%XEDGE      )
    NULLIFY ( HcoState%Grid%YEDGE      )
    NULLIFY ( HcoState%Grid%YSIN       )
    NULLIFY ( HcoState%Grid%AREA_M2    )
    NULLIFY ( HcoState%Grid%BXHEIGHT_M )

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
    HcoState%Phys%Avgdr = 6.022e23_dp
    HcoState%Phys%PI    = 3.14159265358979323_dp
    HcoState%Phys%Re    = 6.375e6_dp
    HcoState%Phys%AIRMW = 28.97_dp
    HcoState%Phys%g0    = 9.8_dp
    HcoState%Phys%Rd    = 287.0_dp
    HcoState%Phys%Rdg0  = HcoState%Phys%Rd / HcoState%Phys%g0

    ! Timesteps
    HcoState%TS_EMIS = 0.0_sp 
    HcoState%TS_CHEM = 0.0_sp
    HcoState%TS_DYN  = 0.0_sp

    ! Nullify temporary array. This array may be used as temporary
    ! place to write emissions into.
    HcoState%Buffer3D => NULL()
    CALL HCO_ArrInit( HcoState%Buffer3D, 0, 0, 0, RC )
    IF ( RC /= 0 ) RETURN

    ! Default HEMCO options
    ! ==> execute HEMCO core; use all species and categories
    ALLOCATE( HcoState%Options )
    HcoState%Options%ExtNr         =  0
    HcoState%Options%SpcMin        =  1
    HcoState%Options%SpcMax        = -1
    HcoState%Options%CatMin        =  1
    HcoState%Options%CatMax        = -1
    HcoState%Options%AutoFillDiagn = .TRUE.
    HcoState%Options%FillBuffer    = .FALSE.

    ! Leave w/ success
    CALL HCO_LEAVE ( RC ) 

  END SUBROUTINE HcoState_Init
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
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
  SUBROUTINE HcoState_Final ( HcoState )
!
! !USES:
!
!
! !INPUT PARAMETERS:
! 
    TYPE(HCO_State), POINTER  :: HcoState    ! HEMCO State object
!
! !REVISION HISTORY: 
!  20 Aug 2013 - C. Keller - Adapted from gigc_state_chm_mod.F90 
!EOP
!------------------------------------------------------------------------------
!BOC

    INTEGER                :: I

    !=====================================================================
    ! HcoState_Final begins here!
    !=====================================================================

    ! Deallocate buffer array
    CALL HCO_ArrCleanup ( HcoState%Buffer3D )

    ! Deallocate all species arrays 
    IF ( ASSOCIATED ( HcoState%Spc ) ) THEN
       DO I = 1, HcoState%nSpc
          CALL HCO_ArrCleanup( HcoState%Spc(I)%Emis )
          CALL HCO_ArrCleanup( HcoState%Spc(I)%Depv )
       ENDDO
       DEALLOCATE ( HcoState%Spc )
    ENDIF

    ! Deallocate grid information
    IF ( ASSOCIATED ( HcoState%Grid) ) THEN
     IF ( ASSOCIATED ( HcoState%Grid%XMID      ) ) DEALLOCATE ( HcoState%Grid%XMID      ) 
     IF ( ASSOCIATED ( HcoState%Grid%YMID      ) ) DEALLOCATE ( HcoState%Grid%YMID      )
     IF ( ASSOCIATED ( HcoState%Grid%XEDGE     ) ) DEALLOCATE ( HcoState%Grid%XEDGE     )
     IF ( ASSOCIATED ( HcoState%Grid%YEDGE     ) ) DEALLOCATE ( HcoState%Grid%YEDGE     )
     IF ( ASSOCIATED ( HcoState%Grid%YSIN      ) ) DEALLOCATE ( HcoState%Grid%YSIN      )
     IF ( ASSOCIATED ( HcoState%Grid%AREA_M2   ) ) DEALLOCATE ( HcoState%Grid%AREA_M2   )
     IF ( ASSOCIATED ( HcoState%Grid%BXHEIGHT_M) ) DEALLOCATE ( HcoState%Grid%BXHEIGHT_M)
     DEALLOCATE(HcoState%Grid)
    ENDIF

    ! Cleanup various types
    IF ( ASSOCIATED ( HcoState%Options ) ) DEALLOCATE ( HcoState%Options )
    IF ( ASSOCIATED ( HcoState%Phys    ) ) DEALLOCATE ( HcoState%Phys    )

#if defined(ESMF_)
    HcoState%IMPORT => NULL()
#endif

  END SUBROUTINE HcoState_Final
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: hco_getmodspcid 
!
! !DESCRIPTION: Function HCO\_GETMODSPCID returns the model species index 
! of a species by name. Returns -1 if given species is not found, 0 if 
! name corresponds to the HEMCO wildcard character.
!\\
!\\
! !INTERFACE:
!
      FUNCTION HCO_GetModSpcID( name, HcoState ) RESULT( Indx )
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
      IF ( TRIM(name) == HCO_WILDCARD() ) THEN
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
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: hco_gethcoid 
!
! !DESCRIPTION: Function HCO\_GETHCOID returns the HEMCO species index 
! of a species by name. Returns -1 if given species is not found, 0 if 
! name corresponds to the HEMCO wildcard character.
!\\
!\\
! !INTERFACE:
!
      FUNCTION HCO_GetHcoID( name, HcoState ) RESULT( Indx )
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
      IF ( TRIM(name) == HCO_WILDCARD() ) THEN
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
END MODULE HCO_STATE_MOD
