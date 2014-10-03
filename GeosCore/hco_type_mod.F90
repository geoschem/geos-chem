!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: hco_type_mod
!
! !DESCRIPTION: Module HCO\_TYPE\_MOD contains type definitions used by 
!  the Harvard Emissions Component (HEMCO). 
!\\
!\\
!\\
! !INTERFACE: 
!
MODULE HCO_TYPE_MOD
!
! USES:
!
  USE HCO_ERROR_MOD

#if defined(ESMF_)
  USE ESMF
#endif

  IMPLICIT NONE

  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC :: Init_HCO_State
  PUBLIC :: Cleanup_HCO_State
  PUBLIC :: AllocCheck 
  PUBLIC :: ResetArrays
  PUBLIC :: HCO_GetSpecID 
!
! !PUBLIC DATA MEMBERS:
!
! 3D container. 
! Contains the 3D flux array (x,y,z) of a given compound.
  TYPE :: Arr3D
     REAL*8,            POINTER :: Arr3d(:,:,:)
  END TYPE Arr3D
!
! 2D container.
! Contains the 2D flux array (x,y) of a given compound.
  TYPE :: Arr2D
     REAL*8,            POINTER :: Arr2d(:,:)
  END TYPE Arr2D

  !=========================================================================
  ! Derived type for HEMCO State
  !=========================================================================
  TYPE, PUBLIC :: HCO_State

     ! Vector containing the 3D emission flux arrays (kg/m2/s) 
     TYPE(Arr3D),        POINTER :: Emsr3D(:) 

     ! Vector containing the 2D deposition velocity arrays (m/s)
     TYPE(Arr2D),        POINTER :: Depv2D(:) 

     ! Tracer information
     INTEGER,            POINTER :: nSpecies
     INTEGER,            POINTER :: SpecIDs   (:)
     CHARACTER(LEN=255), POINTER :: SpecNames (:)
     REAL*8,             POINTER :: SpecMW    (:)
     REAL*8,             POINTER :: EmSpecMW  (:)
     REAL*8,             POINTER :: MolecRatio(:)

     ! HEMCO grid specifications
     LOGICAL                     :: OnSimGrid   ! Emission grid = Simulation grid?
     INTEGER                     :: ISIZE
     INTEGER                     :: JSIZE
     INTEGER                     :: LSIZE
     REAL*8, POINTER             :: XMID       (:,:,:)
     REAL*8, POINTER             :: YMID       (:,:,:)
     REAL*8, POINTER             :: XEDGE      (:,:,:)
     REAL*8, POINTER             :: YSIN       (:,:,:)
     REAL*8, POINTER             :: AREA_M2    (:,:,:)
     REAL*8, POINTER             :: BXHEIGHT_M (:,:,:)

     ! Pointer to temporary array
     REAL*8, POINTER             :: Temp3D(:,:,:) 

     ! Physical constants
     REAL*8                      :: Avgdr   ! Avogadro number (mol-1)
     REAL*8                      :: PI      ! pi
     REAL*8                      :: AIRMW   ! molecular weight of air (g/mol)
     REAL*8                      :: g0      ! gravity at surface of earth (m/s2)

     ! HEMCO parameter
     CHARACTER(LEN=1)            :: WILDCARD

     ! Emission time step in seconds
     REAL*8                     :: TS_EMIS
     REAL*8                     :: TS_DYN

     ! HEMCO options
     INTEGER                    :: ExtNr
     INTEGER                    :: TrcMin
     INTEGER                    :: TrcMax
     INTEGER                    :: CatMin
     INTEGER                    :: CatMax
     INTEGER                    :: sYear
     INTEGER                    :: sMonth
     INTEGER                    :: sDay
     INTEGER                    :: sHour
     INTEGER                    :: sMin
     INTEGER                    :: sSec
     INTEGER                    :: sDayOfYear
     INTEGER                    :: sWeekDay
     LOGICAL                    :: FillTemp3D

     ! HEMCO settings
     CHARACTER(LEN=255)          :: ConfigFile 
     LOGICAL                     :: verbose

     ! If run w/in ESMF, also need to point to IMPORT state 
#if defined(ESMF_)
     TYPE(ESMF_State), POINTER   :: IMPORT
#endif

  END TYPE HCO_State
!
! Container for reading data. This container is mainly used in 
! HCO_READLIST_MOD.F90.
  TYPE, PUBLIC :: RdCont
     CHARACTER(LEN=255)    :: ncFile
     CHARACTER(LEN= 31)    :: ncPara
     CHARACTER(LEN= 31)    :: cName
     INTEGER               :: DataType
     INTEGER               :: ExtNr
     INTEGER               :: cID
     INTEGER               :: ncYrs(2)
     INTEGER               :: ncMts(2)
     INTEGER               :: ncDys(2)
     INTEGER               :: ncHrs(2)
     INTEGER               :: TrcID
     INTEGER               :: ScalID
     INTEGER               :: Cat
     INTEGER               :: Hier
     INTEGER               :: Oper
     INTEGER, POINTER      :: UseScalIDs(:)
     REAL*8,  POINTER      :: Array(:,:,:,:)
     LOGICAL               :: Add
     LOGICAL               :: ncRead
     CHARACTER(LEN= 31)    :: ESMF_Dim
     CHARACTER(LEN= 31)    :: ESMF_Unit
     TYPE(RdCont), POINTER :: NextRdCont
  ENDTYPE RdCont
!
! !REMARKS:
!                                                                             
! !REVISION HISTORY:
!  20 Aug 2013 - C. Keller - Initial version, moved from gigc_state_chm_mod.F90
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
! !IROUTINE: init_hco_state
!
! !DESCRIPTION: Routine INIT\_HCO\_STATE allocates and initializes the 
!  pointer fields of the HEMCO state object.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Init_HCO_State( am_I_Root, nSpecies, SpecIDs,    SpecNames,  &
                             SpecMW,    EmSpecMW, MolecRatio, ConfigFile, &
#if defined(ESMF_)
                             IMPORT,                                      &
#endif
                             TS_EMIS,   TS_DYN, HcoState, RC )
!
! !USES:
!
!
! !INPUT PARAMETERS:
! 
    LOGICAL,          INTENT(IN)          :: am_I_Root            ! root CPU?
    INTEGER,          INTENT(IN)          :: nSpecies             ! # advected tracers
    INTEGER,          INTENT(IN)          :: SpecIDs(nSpecies)    ! Species IDs
    CHARACTER(LEN=*), INTENT(IN)          :: SpecNames(nSpecies)  ! Species names
    REAL*8,           INTENT(IN)          :: SpecMW(nSpecies)     ! Species MW (g/mol) 
    REAL*8,           INTENT(IN)          :: EmSpecMW(nSpecies)   ! Emitted MW (g/mol) 
    REAL*8,           INTENT(IN)          :: MolecRatio(nSpecies) ! emitted mol /
                                                                  ! species mol
    CHARACTER(LEN=*), INTENT(IN)          :: ConfigFile           ! Full path name 
#if defined(ESMF_)
    TYPE(ESMF_State), INTENT(IN), TARGET  :: IMPORT               ! ESMF IMPORT object
#endif
    REAL*8,           INTENT(IN)          :: TS_EMIS              ! Emission timestep [s] 
    REAL*8,           INTENT(IN)          :: TS_DYN               ! Dynamics timestep [s] 
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(HCO_State),  POINTER             :: HcoState             ! HEMCO State object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(INOUT)        :: RC                  ! Return code
!
! !REMARKS:
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
    CHARACTER(LEN=255) :: LOC

    !=====================================================================
    ! Init_HCO_State begins here!
    !=====================================================================

    ! For error handling
    LOC = 'AllocCheck (HCO_TYPE_MOD.F90)'

    !=====================================================================
    ! Allocate emission field vectors
    !=====================================================================

    ! Check if already allocated
    IF ( ASSOCIATED(HcoState)) THEN
       CALL HCO_WARNING( 'HcoState already allocated!', LOC, RC ) 
       RETURN
    ENDIF

    ! Allocate pointer first
    ALLOCATE ( HcoState )

    ! Allocate vectors of flux arrays and deposition velocities
    ALLOCATE( HcoState%Emsr3d ( nSpecies ), STAT=AS )
    IF ( AS /= 0 ) THEN
       CALL HCO_ERROR( 'Emsr3d', LOC, RC ); RETURN
    ENDIF
    ALLOCATE( HcoState%Depv2d ( nSpecies ), STAT=AS )
    IF ( AS /= 0 ) THEN
       CALL HCO_ERROR( 'Depv2d', LOC, RC ); RETURN
    ENDIF

    ! Nullify all flux arrays
    DO I = 1, nSpecies
       NULLIFY ( HcoState%Emsr3d(I)%Arr3d )
       NULLIFY ( HcoState%Depv2d(I)%Arr2d )
    ENDDO !I

    ! Pass HEMCO path name
    HcoState%ConfigFile = TRIM(ConfigFile)

    ! Initalize HEMCO settings
    HcoState%verbose   = .FALSE.

    ! Set tracer settings
    ALLOCATE ( HcoState%nSpecies )

    ALLOCATE ( HcoState%SpecNames (nSpecies ), STAT=AS )
    IF ( AS /= 0 ) THEN
       CALL HCO_ERROR( 'SpecNames', LOC, RC ); RETURN
    ENDIF

    ALLOCATE ( HcoState%SpecIDs   (nSpecies ), STAT=AS )
    IF ( AS /= 0 ) THEN
       CALL HCO_ERROR( 'SpecIDs', LOC, RC ); RETURN
    ENDIF

    ALLOCATE ( HcoState%SpecMW    (nSpecies ), STAT=AS )
    IF ( AS /= 0 ) THEN
       CALL HCO_ERROR( 'SpecMW', LOC, RC ); RETURN
    ENDIF

    ALLOCATE ( HcoState%EmSpecMW  (nSpecies ), STAT=AS )
    IF ( AS /= 0 ) THEN
       CALL HCO_ERROR( 'EmSpecMW', LOC, RC ); RETURN
    ENDIF

    ALLOCATE ( HcoState%MolecRatio(nSpecies ), STAT=AS )
    IF ( AS /= 0 ) THEN
       CALL HCO_ERROR( 'nSpecies', LOC, RC ); RETURN
    ENDIF

    ! Set values
    HcoState%nSpecies      = nSpecies
    HcoState%SpecNames(:)  = SpecNames(:)
    HcoState%SpecIDs(:)    = SpecIDs(:)
    HcoState%SpecMW(:)     = SpecMW(:)
    HcoState%EmSpecMW(:)   = EmSpecMW(:)
    HcoState%MolecRatio(:) = MolecRatio(:)

    ! Initialize GRID settings
    HcoState%OnSimGrid = .TRUE.

    ! Initialize grid dimensions. They will be set in hco_grid_mod
    HcoState%ISIZE   = 0
    HcoState%JSIZE   = 0
    HcoState%LSIZE   = 0

    ! Define physical constants
    HcoState%Avgdr = 6.022d23
    HcoState%PI    = 3.14159265358979323d0 
    HcoState%AIRMW = 28.97d0
    HcoState%g0    = 9.8d0

    ! Define HEMCO parameter
    HcoState%WILDCARD = '*'

    ! Set emission timestep
    HcoState%TS_EMIS = TS_EMIS 
    HcoState%TS_DYN  = TS_DYN

    ! Nullify grid arrays. They will be initialized in hco_grid_mod
    NULLIFY ( HcoState%XMID       )
    NULLIFY ( HcoState%YMID       )
    NULLIFY ( HcoState%XEDGE      )
    NULLIFY ( HcoState%YSIN       )
    NULLIFY ( HcoState%AREA_M2    )
    NULLIFY ( HcoState%BXHEIGHT_M )

    ! Nullify temporary array. This array may be used as temporary
    ! place to write emissions into. 
    NULLIFY ( HcoState%Temp3D )

    ! Set default HEMCO options
    HcoState%ExtNr      =  0
    HcoState%TrcMin     =  1
    HcoState%TrcMax     = -1
    HcoState%CatMin     =  1
    HcoState%CatMax     = -1
    HcoState%sYear      = -1
    HcoState%sMonth     = -1
    HcoState%sDay       = -1
    HcoState%sHour      = -1
    HcoState%sMin       = -1
    HcoState%sSec       = -1
    HcoState%SDayOfYear = -1
    HcoState%sWeekDay   = -1
    HcoState%FillTemp3D = .FALSE.

    ! Point to ESMF IMPORT object 
#if defined(ESMF_)
    HcoState%IMPORT => IMPORT
#endif

    ! Leave w/ success
    RC = HCO_SUCCESS 

  END SUBROUTINE Init_HCO_State
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: cleanup_hco_state
!
! !DESCRIPTION: Routine CLEANUP\_HCO\_STATE deallocates the fields 
!  of the chemistry state object.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Cleanup_HCO_State ( HcoState )
!
! !INPUT PARAMETERS:
! 
    TYPE(HCO_State), POINTER  :: HcoState    ! HEMCO State object
!
! !REMARKS:
!
! !REVISION HISTORY: 
!  20 Aug 2013 - C. Keller - Adapted from gigc_state_chm_mod.F90 
!EOP
!------------------------------------------------------------------------------
!BOC

    INTEGER  :: I

    !=====================================================================
    ! Cleanup_HCO_State begins here!
    !=====================================================================

    ! Deallocate Emsr3D vector
    IF ( ASSOCIATED ( HcoState%Emsr3D ) ) THEN
       DO I = 1, SIZE(HcoState%Emsr3D) 
          IF ( ASSOCIATED(HcoState%Emsr3D(I)%Arr3D) ) THEN
             DEALLOCATE ( HcoState%Emsr3D(I)%Arr3D )
          ENDIF
       ENDDO

       ! Deallocate vector
       DEALLOCATE( HcoState%Emsr3D )
    ENDIF

    ! Deallocate Depv2D vector
    IF ( ASSOCIATED ( HcoState%Depv2D ) ) THEN
       DO I = 1, SIZE(HcoState%Depv2D) 
          IF ( ASSOCIATED(HcoState%Depv2D(I)%Arr2D) ) THEN
             DEALLOCATE ( HcoState%Depv2D(I)%Arr2D )
          ENDIF
       ENDDO

       ! Deallocate vector
       DEALLOCATE( HcoState%Depv2D )
    ENDIF

    ! Cleanup species information
    IF ( ASSOCIATED ( HcoState%nSpecies  ) ) DEALLOCATE ( HcoState%nSpecies  )
    IF ( ASSOCIATED ( HcoState%SpecIDs   ) ) DEALLOCATE ( HcoState%SpecIDs   )
    IF ( ASSOCIATED ( HcoState%SpecNames ) ) DEALLOCATE ( HcoState%SpecNames )
    IF ( ASSOCIATED ( HcoState%SpecMW    ) ) DEALLOCATE ( HcoState%SpecMW    )
    IF ( ASSOCIATED ( HcoState%EmSpecMW  ) ) DEALLOCATE ( HcoState%EmSpecMW  )
    IF ( ASSOCIATED ( HcoState%MolecRatio) ) DEALLOCATE ( HcoState%MolecRatio)

    ! Nullify grid box information 
    HcoState%XMID       => NULL()
    HcoState%YMID       => NULL()
    HcoState%XEDGE      => NULL()
    HcoState%YSIN       => NULL()
    HcoState%AREA_M2    => NULL()
    HcoState%BXHEIGHT_M => NULL()

    ! Nullify temporary array
    HcoState%Temp3D     => NULL()

    ! Eventually remove pointer to ESMF IMPORT object
#if defined(ESMF_)
    HcoState%IMPORT => NULL()
#endif

  END SUBROUTINE Cleanup_HCO_State
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: AllocCheck
!
! !DESCRIPTION: Routine AllocCheck checks if the emission array of a given
! tracer is allocated or not and allocates it if necessary. 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE AllocCheck ( HcoState, TrcID, TPE, RC )
!
! !INPUT/OUTPUT PARAMETERS:
! 
    TYPE(HCO_State), INTENT(INOUT)   :: HcoState
    INTEGER,         INTENT(IN   )   :: TrcID
    INTEGER,         INTENT(IN   )   :: TPE
    INTEGER,         INTENT(INOUT)   :: RC 
!
! !REMARKS:
!
! !REVISION HISTORY: 
!  01 May 2013 - C. Keller - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
    INTEGER            :: I, J, L, AS
    CHARACTER(LEN=255) :: LOC

    !=====================================================================
    ! AllocCheck begins here!
    !=====================================================================

    ! For error handling
    LOC = 'AllocCheck (HCO_TYPE_MOD.F90)'

    ! Check flux array
    IF ( TPE == 1 ) THEN
       IF ( .NOT. ASSOCIATED ( HcoState%Emsr3D(TrcID)%Arr3D ) ) THEN
          I = HcoState%ISIZE
          J = HcoState%JSIZE
          L = HcoState%LSIZE
          ALLOCATE ( HcoState%Emsr3D(TrcID)%Arr3D(I,J,L), STAT=AS )
          IF ( AS /= 0 ) THEN
             CALL HCO_ERROR ( 'Emsr3D allocation', LOC, RC )
             RETURN
          ENDIF

          HcoState%Emsr3D(TrcID)%Arr3D = 0d0
       ENDIF

    ! Check drydep array 
    ELSEIF ( TPE == 2 ) THEN
       IF ( .NOT. ASSOCIATED ( HcoState%Depv2D(TrcID)%Arr2D ) ) THEN
          I = HcoState%ISIZE
          J = HcoState%JSIZE
          ALLOCATE ( HcoState%Depv2D(TrcID)%Arr2D(I,J), STAT=AS )
          IF ( AS /= 0 ) THEN
             CALL HCO_ERROR ( 'Depv2D allocation', LOC, RC )
             RETURN
          ENDIF
          HcoState%Depv2D(TrcID)%Arr2D = 0d0
       ENDIF

    ENDIF

    ! Return w/ success
    RC = HCO_SUCCESS

  END SUBROUTINE AllocCheck
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ResetArrays 
!
! !DESCRIPTION: Routine ResetArrays sets all defined 3D flux arrays (Emsr3D)
! and 2D deposition arrays (Depv2D) of the passed HEMCO state object to zero. 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE ResetArrays ( HcoState, RC )
!
! !INPUT/OUTPUT PARAMETERS:
! 
    TYPE(HCO_State), INTENT(INOUT)   :: HcoState
    INTEGER,         INTENT(INOUT)   :: RC 
!
! !REMARKS:
!
! !REVISION HISTORY: 
!  01 May 2013 - C. Keller - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
    INTEGER            :: N 

    !=====================================================================
    ! ResetArrays begins here!
    !=====================================================================

    ! Loop over all arrays. 
    DO N = 1, HcoState%nSpecies

       ! 3D flux rates array
       IF ( ASSOCIATED(HcoState%Emsr3D(N)%Arr3D) ) THEN
          HcoState%Emsr3D(N)%Arr3D = 0d0
       ENDIF

       ! 2D deposition velocity array
       IF ( ASSOCIATED(HcoState%Depv2D(N)%Arr2D) ) THEN
          HcoState%Depv2D(N)%Arr2D = 0d0
       ENDIF
    ENDDO

    ! Return w/ success
    RC = HCO_SUCCESS

  END SUBROUTINE ResetArrays 
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: hco_getspecid 
!
! !DESCRIPTION: Function HCO\_GETSPECID returns the index of an species 
!  contained in HEMCO state object by name.
!\\
!\\
! !INTERFACE:
!
      FUNCTION HCO_GetSpecID( name, HcoState ) RESULT( Indx )
!
! !INPUT PARAMETERS:
!
      CHARACTER(LEN=*), INTENT(IN)   :: name         ! Species or tracer name
      TYPE(HCO_State), INTENT(INOUT) :: HcoState     ! HEMCO State
!
! !RETURN VALUE:
!
      INTEGER                      :: Indx         ! Index of this species 
!
! !REVISION HISTORY: 
!  09 Oct 2012 - M. Long     - Initial version, based on gc_esmf_utils_mod.F90
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER :: M

      ! Initialize
      Indx = -1

      ! Return 0 if wildcard character
      IF ( TRIM(name) == TRIM(HcoState%WILDCARD) ) THEN
         Indx = 0
         RETURN
      ENDIF

      ! Loop over all species names
      DO M = 1, SIZE( HcoState%SpecNames )

         ! Return the index of the sought-for species
         IF( TRIM( name ) == TRIM( HcoState%SpecNames(M) ) ) THEN
            Indx = HcoState%SpecIDs(M)
            EXIT
         ENDIF

      ENDDO

      END FUNCTION HCO_GetSpecID
!EOC
END MODULE HCO_TYPE_MOD
