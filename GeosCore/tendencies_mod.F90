#if defined( USE_TEND )
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: tendencies_mod.F90
!
! !DESCRIPTION: Module tendencies\_mod.F90 is a module to define, archive and
! write species tendencies. This module is still under development and only
! active if the DEVEL compiler flag is enabled. Also, the species to be
! diagnosed as well as the processes for which tendencies shall be calculated
! are currently hardcoded below. If enabled, this module will calculate 
! concentration tendencies for all defined species and processes and write 
! these into netCDF diagnostics. All tendencies are given as v/v/s.
!\\
!\\
! Species tendencies can be archived for as many processes are desired. For 
! each tendency process, an own instance of a tendency class has to be 
! defined (subroutine Tend\_CreateClass) and all active species for this class
! need be defined via subroutine Tend\_Add. Once a tendencies class is defined,
! the entry and exit 'checkpoint' of this tendency need be manually set in the
! code (subroutines Tend\_Stage1 and Tend\_Stage2).
!\\
!\\
! For example, suppose there is routine PROCESS in module example_mod.F90 and 
! we are interested in the species tendencies of O3 and CO by this process. We
! can then define a new tendency class (named 'PROCESS') during initialization
! of the tendencies (i.e. in tend\_init):
!
!    ! Create new class
!    CALL Tend_CreateClass( am_I_Root, Input_Opt, State_Chm, 'PROCESS', RC )
!    IF ( RC /= GC_SUCCESS ) RETURN
!
! The second step is to assign the species of interest to this tendency class:
!
!    ! Add species to classes
!    CALL Tend_Add ( am_I_Root, Input_Opt, 'PROCESS', id_O3, RC )
!    IF ( RC /= GC_SUCCESS ) RETURN
!    CALL Tend_Add ( am_I_Root, Input_Opt, 'PROCESS', id_CO, RC )
!    IF ( RC /= GC_SUCCESS ) RETURN
!
! The last step then involves the definition of the entry and exit points of
! the tendencies, e.g. the interval in between the tendencies shall be 
! calculated. To do so, we can wrap the Tend\_Stage routines around the 
! process of interest, e.g. in module example_mod.F90:
!
! CALL Tend\_Stage1 ( ... TendName='PROCESS', ... )
! CALL PROCESS ( ... )
! CALL Tend\_Stage2 ( ..., TendName='PROCESS', ... ) 
!\\
!\\
! The following six tendency classes are implemented by default: ADV (transport),
! CONV (convection), CHEM (chemistry), WETD (wet deposition), PBLMIX (PBL mixing, 
! includes emissions and dry deposition below PBL if non-local PBL is enabled),
! FLUX (emissions and dry depositions not coverd in PBLMIX).
! Subroutine Tend\_Init contains some example tendencies that are calculated 
! if flag 'DoTend' (subroutine Tend\_Init) is enabled. 
!
! !INTERFACE:
!
MODULE Tendencies_Mod 
!
! !USES:
!
  USE ErrCode_Mod
  USE HCO_Error_Mod
  USE HCO_Diagn_Mod
  USE Input_Opt_Mod,      ONLY : OptInput
  USE Precision_Mod
  USE State_Chm_Mod,      ONLY : ChmState
  USE State_Met_Mod,      ONLY : MetState

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC :: Tend_Init
  PUBLIC :: Tend_CreateClass
  PUBLIC :: Tend_Add
  PUBLIC :: Tend_Stage1
  PUBLIC :: Tend_Stage2
  PUBLIC :: Tend_Get
  PUBLIC :: Tend_Cleanup
!
! !PRIVATE MEMBER FUNCTIONS:
!
  PRIVATE :: Tend_FindClass
!
! !PUBLIC DATA MEMBERS:
!
  ! maximum string length of tendency name
  INTEGER, PARAMETER   :: MAXSTR = 31 

  ! Number of GEOS-Chem species
  INTEGER              :: nSpc = 0

  ! vector of tendency arrays
  TYPE :: TendArr
     REAL(f4),       POINTER   :: Arr(:,:,:) => NULL()
  END TYPE TendArr

  ! type holding tendencies of one type (class). Will
  ! be linked together via linked list.
  TYPE :: TendClass
     CHARACTER(LEN=MAXSTR)     :: TendName
     INTEGER                   :: Stage
     INTEGER,         POINTER  :: SpcInUse(:) => NULL()
     TYPE(TendArr),   POINTER  :: Tendency(:) => NULL()
     TYPE(TendClass), POINTER  :: NextTend    => NULL()
  END TYPE TendClass

  ! Tendency class linked list
  TYPE(TendClass),   POINTER   :: TendList    => NULL()
!
! !REVISION HISTORY:
!  14 Jul 2015 - C. Keller   - Initial version. 
!  26 Oct 2015 - C. Keller   - Now organize in linked list for more flexibility.
!  19 Jul 2016 - R. Yantosca - Now block out this routine with #ifdef USE_TEND
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
! !IROUTINE: Tend_Init
!
! !DESCRIPTION: Subroutine Tend\_Init is the wrapper routine to initialize the
! tendencies. At the moment, all tendencies are hardcoded and tendencies will 
! only be written if the manual flag `DoTend` is enabled. 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Tend_Init ( am_I_Root, Input_Opt, State_Met, State_Chm, RC ) 
!
! !USES:
!
     USE State_Chm_Mod,      ONLY : Ind_
!
! !INPUT PARAMETERS:
!
    LOGICAL,          INTENT(IN   ) :: am_I_Root  ! Are we on the root CPU?
    TYPE(OptInput),   INTENT(IN   ) :: Input_Opt  ! Input opts
    TYPE(MetState),   INTENT(IN   ) :: State_Met  ! met. state 
    TYPE(ChmState),   INTENT(IN   ) :: State_Chm  ! chm. state 
!
! !OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(OUT)   :: RC         ! Failure or success
!
! !REMARKS:
!  This subroutine gets called from Diagnostics_Init (in diagnostics_mod.F90).
!  It is only executed once.
!
! !REVISION HISTORY: 
!  26 Oct 2015 - C. Keller   - Initial version 
!  16 Jun 2016 - J. Kaiser   - Move tracer IDS to variable names
!  20 Jun 2016 - R. Yantosca - Renamed IDTCO, IDTO3 to id_CO and id_O3
!  19 Jul 2016 - R. Yantosca - Activate more tendency classes
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Species ID flags
    INTEGER            :: id_CO
    INTEGER            :: id_O3
    INTEGER            :: id_NO
    INTEGER            :: id_NO2
    INTEGER            :: id_HNO3
    INTEGER            :: id_tmp

    ! Strings
    CHARACTER(LEN=255) :: MSG
    CHARACTER(LEN=255) :: LOC = 'Tend_Init (tendencies_mod.F)' 
!
! !DEFINED PARAMETERS:
!
    ! Set this to .TRUE. to enable some test diagnostics
    LOGICAL, PARAMETER :: DoTend = .TRUE.
 
    !=======================================================================
    ! Tend_Init begins here!
    !=======================================================================

    ! Assume successful return
    RC = GC_SUCCESS

    ! Define species ID flags
    id_CO   = Ind_('CO')
    id_O3   = Ind_('O3') 
    id_NO   = Ind_('NO') 
    id_NO2  = Ind_('NO2') 
    id_HNO3 = Ind_('HNO3') 
       
    ! Execute only if DoTend is enabled
    IF ( DoTend ) THEN

       !--------------------------------------------------------------------
       ! Define tendency classes (add more as you wish)
       !--------------------------------------------------------------------

       CALL Tend_CreateClass( am_I_Root, Input_Opt, State_Chm, 'ADV' ,   RC )
       IF ( RC /= GC_SUCCESS ) RETURN

       CALL Tend_CreateClass( am_I_Root, Input_Opt, State_Chm, 'CHEM',   RC )
       IF ( RC /= GC_SUCCESS ) RETURN

       CALL Tend_CreateClass( am_I_Root, Input_Opt, State_Chm, 'CONV',   RC )
       IF ( RC /= GC_SUCCESS ) RETURN

       CALL Tend_CreateClass( am_I_Root, Input_Opt, State_Chm, 'FLUX',   RC )
       IF ( RC /= GC_SUCCESS ) RETURN

       CALL Tend_CreateClass( am_I_Root, Input_Opt, State_Chm, 'PBLMIX', RC )
       IF ( RC /= GC_SUCCESS ) RETURN

       CALL Tend_CreateClass( am_I_Root, Input_Opt, State_Chm, 'WETD',   RC )
       IF ( RC /= GC_SUCCESS ) RETURN

       !--------------------------------------------------------------------
       ! Activate tendency computations for O3 (add more as you wish)
       !--------------------------------------------------------------------

       CALL Tend_Add ( am_I_Root, Input_Opt, State_Chm, 'ADV',    id_O3, RC )
       IF ( RC /= GC_SUCCESS ) RETURN

       CALL Tend_Add ( am_I_Root, Input_Opt, State_Chm, 'CHEM',   id_O3, RC )
       IF ( RC /= GC_SUCCESS ) RETURN

       CALL Tend_Add ( am_I_Root, Input_Opt, State_Chm, 'CONV',   id_O3, RC )
       IF ( RC /= GC_SUCCESS ) RETURN

       CALL Tend_Add ( am_I_Root, Input_Opt, State_Chm, 'FLUX',   id_O3, RC )
       IF ( RC /= GC_SUCCESS ) RETURN

       CALL Tend_Add ( am_I_Root, Input_Opt, State_Chm, 'WETD',   id_O3, RC )
       IF ( RC /= GC_SUCCESS ) RETURN

       CALL Tend_Add ( am_I_Root, Input_Opt, State_Chm, 'PBLMIX', id_O3, RC )
       IF ( RC /= GC_SUCCESS ) RETURN

       !--------------------------------------------------------------------
       ! Activate tendency computations for CO (add more as you wish)
       !--------------------------------------------------------------------

       CALL Tend_Add ( am_I_Root, Input_Opt, State_Chm, 'ADV',    id_CO, RC )
       IF ( RC /= GC_SUCCESS ) RETURN

       CALL Tend_Add ( am_I_Root, Input_Opt, State_Chm, 'CHEM',   id_CO, RC )
       IF ( RC /= GC_SUCCESS ) RETURN

       CALL Tend_Add ( am_I_Root, Input_Opt, State_Chm, 'CONV',   id_CO, RC )
       IF ( RC /= GC_SUCCESS ) RETURN

       CALL Tend_Add ( am_I_Root, Input_Opt, State_Chm, 'PBLMIX', id_CO, RC )
       IF ( RC /= GC_SUCCESS ) RETURN

       CALL Tend_Add ( am_I_Root, Input_Opt, State_Chm, 'WETD', id_CO, RC )
       IF ( RC /= GC_SUCCESS ) RETURN

       CALL Tend_Add ( am_I_Root, Input_Opt, State_Chm, 'FLUX', id_CO, RC )
       IF ( RC /= GC_SUCCESS ) RETURN

       !--------------------------------------------------------------------
       ! Activate tendency computations for NO (add more as you wish)
       !--------------------------------------------------------------------

       CALL Tend_Add ( am_I_Root, Input_Opt, State_Chm, 'ADV',    id_NO, RC )
       IF ( RC /= GC_SUCCESS ) RETURN

       CALL Tend_Add ( am_I_Root, Input_Opt, State_Chm, 'CHEM',   id_NO, RC )
       IF ( RC /= GC_SUCCESS ) RETURN

       CALL Tend_Add ( am_I_Root, Input_Opt, State_Chm, 'CONV',   id_NO, RC )
       IF ( RC /= GC_SUCCESS ) RETURN

       CALL Tend_Add ( am_I_Root, Input_Opt, State_Chm, 'PBLMIX', id_NO, RC )
       IF ( RC /= GC_SUCCESS ) RETURN

       CALL Tend_Add ( am_I_Root, Input_Opt, State_Chm, 'WETD', id_NO, RC )
       IF ( RC /= GC_SUCCESS ) RETURN

       CALL Tend_Add ( am_I_Root, Input_Opt, State_Chm, 'FLUX', id_NO, RC )
       IF ( RC /= GC_SUCCESS ) RETURN

       !--------------------------------------------------------------------
       ! Activate tendency computations for NO2 (add more as you wish)
       !--------------------------------------------------------------------

       CALL Tend_Add ( am_I_Root, Input_Opt, State_Chm, 'ADV',    id_NO2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN

       CALL Tend_Add ( am_I_Root, Input_Opt, State_Chm, 'CHEM',   id_NO2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN

       CALL Tend_Add ( am_I_Root, Input_Opt, State_Chm, 'CONV',   id_NO2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN

       CALL Tend_Add ( am_I_Root, Input_Opt, State_Chm, 'PBLMIX', id_NO2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN

       CALL Tend_Add ( am_I_Root, Input_Opt, State_Chm, 'WETD', id_NO2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN

       CALL Tend_Add ( am_I_Root, Input_Opt, State_Chm, 'FLUX', id_NO2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN

       !--------------------------------------------------------------------
       ! Activate tendency computations for HNO3 (add more as you wish)
       !--------------------------------------------------------------------

       CALL Tend_Add ( am_I_Root, Input_Opt, State_Chm, 'ADV',    id_HNO3, RC )
       IF ( RC /= GC_SUCCESS ) RETURN

       CALL Tend_Add ( am_I_Root, Input_Opt, State_Chm, 'CHEM',   id_HNO3, RC )
       IF ( RC /= GC_SUCCESS ) RETURN

       CALL Tend_Add ( am_I_Root, Input_Opt, State_Chm, 'CONV',   id_HNO3, RC )
       IF ( RC /= GC_SUCCESS ) RETURN

       CALL Tend_Add ( am_I_Root, Input_Opt, State_Chm, 'PBLMIX', id_HNO3, RC )
       IF ( RC /= GC_SUCCESS ) RETURN

       CALL Tend_Add ( am_I_Root, Input_Opt, State_Chm, 'WETD', id_HNO3, RC )
       IF ( RC /= GC_SUCCESS ) RETURN

       CALL Tend_Add ( am_I_Root, Input_Opt, State_Chm, 'FLUX', id_HNO3, RC )
       IF ( RC /= GC_SUCCESS ) RETURN

       id_tmp = Ind_('OH') 
       CALL Tend_Add ( am_I_Root, Input_Opt, State_Chm, 'CHEM', id_tmp, RC )
       IF ( RC /= GC_SUCCESS ) RETURN

    ENDIF ! test toggle

  END SUBROUTINE Tend_Init
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Tend_CreateClass 
!
! !DESCRIPTION: Subroutine Tend\_CreateClass creates a new tendency class. 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Tend_CreateClass ( am_I_Root, Input_Opt, State_Chm, TendName, RC ) 
!
! !USES:
!
!
! !INPUT PARAMETERS:
!
    LOGICAL,          INTENT(IN   ) :: am_I_Root  ! Are we on the root CPU?
    TYPE(OptInput),   INTENT(IN   ) :: Input_Opt  ! Input Options object
    TYPE(ChmState),   INTENT(IN   ) :: State_Chm  ! Chemistry State object
    CHARACTER(LEN=*), INTENT(IN   ) :: TendName   ! Tendency class name
!
! !OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(OUT)   :: RC         ! Failure or success
!
! !REVISION HISTORY: 
!  26 Oct 2015 - C. Keller   - Initial version 
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    TYPE(TendClass), POINTER :: NewTend => NULL()
    LOGICAL                  :: FOUND
    CHARACTER(LEN=255)       :: MSG
    CHARACTER(LEN=255)       :: LOC = 'Tend_CreateClass (tendencies_mod.F)' 
    
    !=======================================================================
    ! Tend_CreateClass begins here!
    !=======================================================================

    ! Assume successful return
    RC = GC_SUCCESS
    
    ! Check if class already exists
    CALL Tend_FindClass( am_I_Root, TendName, FOUND, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    IF ( .NOT. FOUND ) THEN 

       ! Eventually set local # of species
       IF ( nSpc <= 0 ) THEN
          !nSpc = State_Chm%nAdvect
          nSpc = State_Chm%nSpecies
       ENDIF

       ! Initialize class
       ALLOCATE(NewTend)

       ! Set tendency class name
       NewTend%TendName = TRIM(TendName)

       ! Initialize stage
       newTend%Stage    = -1

       ! Initialize vector with species flags
       ALLOCATE(NewTend%SpcInUse(nSpc))
       NewTend%SpcInUse(:) = 0

       ! Initialize tendency arrays (only allocated when needed)
       ALLOCATE(NewTend%Tendency(nSpc))

       ! Add tendency class to linked list
       NewTend%NextTend => TendList
       TendList         => NewTend 
    ENDIF

  END SUBROUTINE Tend_CreateClass
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Tend_FindClass
!
! !DESCRIPTION: Subroutine Tend\_FindClass searches for a tendency class. 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Tend_FindClass ( am_I_Root, TendName, FOUND, RC, ThisTend ) 
!
! !USES:
!
!
! !INPUT PARAMETERS:
!
    LOGICAL,          INTENT(IN   )          :: am_I_Root  ! Are we on the root CPU?
    CHARACTER(LEN=*), INTENT(IN   )          :: TendName   ! tendency class name
!
! !OUTPUT PARAMETERS:
!
    LOGICAL,          INTENT(  OUT)          :: FOUND      ! class found
    INTEGER,          INTENT(  OUT)          :: RC         ! Failure or success
    TYPE(TendClass),  POINTER,      OPTIONAL :: ThisTend   ! Pointer to this class
!
! !REVISION HISTORY: 
!  26 Oct 2015 - C. Keller   - Initial version
!  19 Jul 2016 - R. Yantosca - Don't nullify local pointers in declarations
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    TYPE(TendClass), POINTER :: TmpTend
    
    !=======================================================================
    ! Tend_FindClass begins here!
    !=======================================================================

    ! Assume successful return
    RC = GC_SUCCESS

    ! Init
    FOUND   = .FALSE.
    TmpTend => NULL()

    ! Loop through linked list and search for class with same name
    TmpTend => TendList
    DO WHILE ( ASSOCIATED(TmpTend) ) 
    
       ! Is this the tendency of interest?  
       IF ( TRIM(TmpTend%TendName) == TRIM(TendName) ) THEN
          FOUND = .TRUE.
          EXIT
       ENDIF
 
       ! Advance in list
       TmpTend => TmpTend%NextTend
    END DO

    ! Eventually 
    IF ( PRESENT(ThisTend) ) ThisTend => TmpTend

    ! Cleanup
    TmpTend => NULL()

  END SUBROUTINE Tend_FindClass
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Tend_Cleanup
!
! !DESCRIPTION: Subroutine Tend\_Cleanup cleans up the tendencies linked list. 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Tend_Cleanup ( ) 
!
! !USES:
!
!
! !INPUT PARAMETERS:
!
!
! !REVISION HISTORY: 
!  26 Oct 2015 - C. Keller   - Initial version
!  19 Jul 2016 - R. Yantosca - Don't nullify local pointers in declarations
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER                  :: I
    TYPE(TendClass), POINTER :: ThisTend
    TYPE(TendClass), POINTER :: NextTend
    
    !=======================================================================
    ! Tend_Cleanup begins here!
    !=======================================================================

    ! Initialize
    ThisTend => NULL()
    NextTend => NULL()

    ! Loop through linked list and search for class with same name
    ThisTend => TendList
    DO WHILE ( ASSOCIATED(ThisTend) ) 
   
       ! Get pointer to next tendency
       NextTend => ThisTend%NextTend

       ! Cleanup every array
       DO I = 1, nSpc
          IF ( ASSOCIATED(ThisTend%Tendency(I)%Arr) ) THEN
             DEALLOCATE( ThisTend%Tendency(I)%Arr)
          ENDIF
       ENDDO
       DEALLOCATE(ThisTend%Tendency)
       DEALLOCATE(ThisTend%SpcInUse)

       ! Cleanup
       ThisTend%NextTend => NULL()
       NULLIFY(ThisTend)
 
       ! Advance in list
       ThisTend => NextTend
    END DO

    ! Cleanup
    ThisTend => NULL()
    NextTend => NULL()
    nSpc = 0

  END SUBROUTINE Tend_Cleanup
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Tend_Add
!
! !DESCRIPTION: Subroutine Tend\_Add adds a species to a tendency class. 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Tend_Add ( am_I_Root, Input_Opt, State_Chm, TendName, &
                        SpcID, RC, CreateClass )
!
! !USES:
!
    USE HCO_INTERFACE_MOD,  ONLY : HcoState
    USE CMN_SIZE_MOD,       ONLY : IIPAR, JJPAR, LLPAR
!
! !INPUT PARAMETERS:
!
    LOGICAL,          INTENT(IN   )           :: am_I_Root   ! Are we on the root CPU?
    TYPE(OptInput),   INTENT(IN   )           :: Input_Opt   ! Input opts
    TYPE(ChmState),   INTENT(IN   )           :: State_Chm   ! Chemistry State object
    CHARACTER(LEN=*), INTENT(IN   )           :: TendName    ! Tendency class name 
    INTEGER,          INTENT(IN   )           :: SpcID       ! Species ID 
    LOGICAL,          INTENT(IN   ), OPTIONAL :: CreateClass ! Create class if missing?
!
! !OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(  OUT)           :: RC          ! Failure or success
!
! !REVISION HISTORY: 
!  14 Jul 2015 - C. Keller   - Initial version 
!  26 Oct 2015 - C. Keller   - Update for linked list 
!  19 Jul 2016 - R. Yantosca - Don't nullify local pointers in declarations
!  29 Dec 2017 - C. Keller   - Don't use HEMCO diagnostics in ESMF env.
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER                  :: AS
    INTEGER                  :: Collection
    LOGICAL                  :: FOUND 

    ! Pointers
    TYPE(TendClass), POINTER :: ThisTend

    ! Strings
    CHARACTER(LEN=63)        :: DiagnName
    CHARACTER(LEN=255)       :: ErrMsg
    CHARACTER(LEN=255)       :: ThisLoc 
    
    !=======================================================================
    ! Tend_Add begins here!
    !=======================================================================

    ! Assume successful return
    RC      = GC_SUCCESS
    ErrMsg  = ''
    ThisLoc = '- > at Tend_Add (in module GeosCore/tendencies_mod.F90)' 

    ! Initialize
    ThisTend => NULL()

    ! Ignore invalid species IDs
    IF ( SpcID <= 0 ) RETURN

    ! Search for diagnostics class
    CALL Tend_FindClass( am_I_Root, TendName, FOUND, RC, ThisTend=ThisTend ) 

    ! Eventually create this class if it does not exist yet
    IF ( .NOT. FOUND .AND. PRESENT( CreateClass ) ) THEN
       IF ( CreateClass ) THEN

          ! Create class
          CALL Tend_CreateClass( am_I_Root, Input_Opt, State_Chm, TendName, RC )
          IF ( RC /= GC_SUCCESS ) RETURN

          ! Get pointer to class object 
          CALL Tend_FindClass( am_I_Root, TendName, FOUND, RC, &
                               ThisTend=ThisTend ) 

       ENDIF
    ENDIF

    ! Return here if class not found
    IF ( .NOT. FOUND ) RETURN

    ! Species ID must not exceed # of tendency species
    IF ( SpcID > nSpc ) THEN
       WRITE(ErrMsg,*) 'Species ID exceeds number of tendency species: ', SpcID, ' > ', nSpc
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF 
 
    ! Name of diagnostics 
    DiagnName = 'TEND_' // TRIM(TendName) // '_' //   &
                TRIM( State_Chm%SpcData(SpcID)%Info%Name )

    ! Mark species as being used
    ThisTend%SpcInUse(SpcID) = 1

    ! Make sure array is allocated
    IF ( .NOT. ASSOCIATED(ThisTend%Tendency(SpcID)%Arr) ) THEN
       ALLOCATE(ThisTend%Tendency(SpcID)%Arr(IIPAR,JJPAR,LLPAR),STAT=RC)
       IF ( RC /= 0 ) THEN
          ErrMsg = 'Tendency allocation error: ' // TRIM(DiagnName)
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
    ENDIF

#if !defined(ESMF_) && defined( MODEL_GEOS )
    ! Get diagnostic parameters from the Input_Opt object
    Collection = Input_Opt%DIAG_COLLECTION

    ! Create container for tendency
    CALL Diagn_Create( am_I_Root, &
                       HcoState  = HcoState,          & 
                       Col       = Collection,        & 
!                       cID       = cID,               &
                       cName     = TRIM( DiagnName ), &
                       AutoFill  = 0,                 &
                       ExtNr     = -1,                &
                       Cat       = -1,                &
                       Hier      = -1,                &
                       HcoID     = -1,                &
                       SpaceDim  =  3,                &
                       OutUnit   = 'v/v/s',           &
                       OutOper   = 'Mean',            &
                       OkIfExist = .TRUE.,            &
                       RC        = RC )
   
    IF ( RC /= HCO_SUCCESS ) THEN
       ErrMsg = 'Cannot create diagnostics: ' // TRIM(DiagnName)
       CALL GC_Error( ErrMsg, RC, ThisLoc ) 
       RETURN
    ENDIF
#endif

    ! Initialize
    ThisTend => NULL()

  END SUBROUTINE Tend_Add
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Tend_Stage1
!
! !DESCRIPTION: Subroutine Tend\_Stage1 archives the current species 
! concentrations into the local tendency arrays.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Tend_Stage1( am_I_Root, Input_Opt, State_Met, &
                          State_Chm, TendName,  RC ) 
!
! !USES:
!
    USE UNITCONV_MOD
    USE PHYSCONSTANTS,      ONLY : AIRMW
!
! !INPUT PARAMETERS:
!
    LOGICAL,          INTENT(IN   ) :: am_I_Root  ! Are we on the root CPU?
    TYPE(OptInput),   INTENT(IN   ) :: Input_Opt  ! Input opts
    TYPE(MetState),   INTENT(IN   ) :: State_Met  ! Met state
    CHARACTER(LEN=*), INTENT(IN   ) :: TendName   ! tendency name 
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ChmState),   INTENT(INOUT) :: State_Chm  ! Chemistry state 
!
! !OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(OUT)   :: RC         ! Failure or success
!
! !REVISION HISTORY: 
!  14 Jul 2015 - C. Keller   - Initial version 
!  26 Oct 2015 - C. Keller   - Update to include tendency classes
!  22 Jun 2016 - M. Yannetti - Replace TCVV with species db MW and phys constant
!  19 Jul 2016 - R. Yantosca - Don't nullify local pointers in declarations
!  19 Jul 2016 - R. Yantosca - Now use State_Chm%Species
!  17 Oct 2017 - C. Keller   - Stage2 now updates internal array to tendency
!                              array (instead of resetting it to zero).
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER                  :: I
    LOGICAL                  :: FOUND

    ! Pointers
    REAL(f4),        POINTER :: Ptr3D(:,:,:)
    TYPE(TendClass), POINTER :: ThisTend

    ! Strings
    CHARACTER(LEN=63)        :: OrigUnit 
    CHARACTER(LEN=255)       :: MSG
    CHARACTER(LEN=255)       :: LOC = 'TEND_STAGE1 (tendencies_mod.F)' 
   
    !=======================================================================
    ! TEND_STAGE1 begins here!
    !=======================================================================

    ! Assume successful return
    RC = GC_SUCCESS

    ! Initialize
    Ptr3D    => NULL()
    ThisTend => NULL()

    ! Find tendency class
    CALL Tend_FindClass( am_I_Root, TendName, FOUND, RC, ThisTend=ThisTend )
    IF ( .NOT. FOUND .OR. .NOT. ASSOCIATED(ThisTend) ) RETURN

    ! Convert tracers to kg/kg dry
    CALL Convert_Spc_Units( am_I_Root, Input_Opt, State_Met, &
                            State_Chm, 'kg/kg dry', RC, OrigUnit=OrigUnit )
    IF ( RC/= HCO_SUCCESS ) RETURN

    ! Loop over # of tendencies species
    DO I = 1, nSpc 

       ! Skip if species is not in use    
       IF ( ThisTend%SpcInUse(I) <= 0 ) CYCLE

       ! Get pointer to 3D array to be filled 
       Ptr3D => ThisTend%Tendency(I)%Arr

       ! Fill 3D array with current values. 
       ! Convert from kg/kg dry to v/v dry
       Ptr3D = 0.0_f4
       Ptr3D = State_Chm%Species(:,:,:,I) &
             * ( AIRMW / State_Chm%SpcData(I)%Info%emMW_g )
!!!             / State_Met%AD(:,:,:)

       ! Cleanup
       Ptr3D => NULL()
    ENDDO !I

    ! Update stage 
    ThisTend%Stage = 1

    ! Convert tracers back to original unit 
    CALL Convert_Spc_Units( am_I_Root, Input_Opt, State_Met, &
                            State_Chm, OrigUnit,  RC )
    IF ( RC/= HCO_SUCCESS ) RETURN

    ! Cleanup
    ThisTend => NULL()

  END SUBROUTINE Tend_Stage1
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Tend_Stage2
!
! !DESCRIPTION: Subroutine Tend_Stage2 calculates the tendencies and 
! writes them into the diagnostics arrays. 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Tend_Stage2( am_I_Root, Input_Opt, State_Met, &
                          State_Chm, TendName,  DT, RC ) 
!
! !USES:
!
    USE UNITCONV_MOD
    USE PHYSCONSTANTS,      ONLY : AIRMW
    USE HCO_INTERFACE_MOD,  ONLY : HcoState
    USE CMN_SIZE_MOD,       ONLY : IIPAR, JJPAR, LLPAR
!
! !INPUT PARAMETERS:
!
    LOGICAL,          INTENT(IN   ) :: am_I_Root  ! Are we on the root CPU?
    TYPE(OptInput),   INTENT(IN   ) :: Input_Opt  ! Input opts
    TYPE(MetState),   INTENT(IN   ) :: State_Met  ! Met state
    CHARACTER(LEN=*), INTENT(IN   ) :: TendName   ! tendency name 
    REAL(fp),         INTENT(IN   ) :: DT         ! delta time, in seconds 
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ChmState),   INTENT(INOUT) :: State_Chm  ! Chemistry state 
!
! !OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(OUT)   :: RC         ! Failure or success
!
! !REVISION HISTORY: 
!  14 Jul 2015 - C. Keller   - Initial version 
!  26 Oct 2015 - C. Keller   - Update to include tendency classes
!  05 Jan 2015 - C. Keller   - Small updates to remove spurious tendencies
!                              caused by floating point errors.
!  22 Jun 2016 - M. Yannetti - Replace TCVV with species db MW and phys constant
!  19 Jul 2016 - R. Yantosca - Don't nullify local pointers in declarations
!  19 Jul 2016 - R. Yantosca - Now use State_Chm%Species
!  29 Dec 2017 - C. Keller   - Don't use HEMCO diagnostics in ESMF env.
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    LOGICAL                  :: ZeroTend
    LOGICAL                  :: FOUND
    INTEGER                  :: cID, I

    ! Arrays
    REAL(f4)                 :: TEND(IIPAR,JJPAR,LLPAR)
    REAL(f4)                 :: TMP (IIPAR,JJPAR,LLPAR)

    ! Pointers
    REAL(f4),        POINTER :: Ptr3D(:,:,:)
    TYPE(TendClass), POINTER :: ThisTend

    ! Strings
    CHARACTER(LEN=63)        :: DiagnName
    CHARACTER(LEN=63)        :: OrigUnit 
    CHARACTER(LEN=255)       :: ErrMsg
    CHARACTER(LEN=255)       :: ThisLoc 
  
    !=======================================================================
    ! TEND_STAGE2 begins here!
    !=======================================================================

    ! Assume successful return
    RC      = GC_SUCCESS
    ErrMsg  = ''
    ThisLoc = ' -> at TEND_STAGE2 (in module GeosCore/tendencies_mod.F90)' 

    ! Initialize
    Ptr3d    => NULL()
    ThisTend => NULL()

    ! Find tendency class
    CALL Tend_FindClass( am_I_Root, TendName, FOUND, RC, ThisTend=ThisTend )
    IF ( .NOT. FOUND .OR. .NOT. ASSOCIATED(ThisTend) ) RETURN

    ! Error check: stage 2 must be called after stage 1
    ZeroTend = .FALSE.
    IF ( ThisTend%Stage /= 1 ) THEN
       IF ( am_I_Root ) THEN
          WRITE(*,*) 'Warning: cannot do tendency stage 2 - stage 1 not yet called: ', TRIM(TendName) 
       ENDIF
       ZeroTend = .TRUE.
    ENDIF

    ! Error check: DT must not be 0
    IF ( DT == 0.0_fp ) THEN
       IF ( am_I_Root ) THEN
          WRITE(*,*) 'Warning: cannot calculate tendency for DT = 0: ', TRIM(TendName)
       ENDIF
       ZeroTend = .TRUE.
    ENDIF

    ! Convert tracers to kg/kg dry
    CALL Convert_Spc_Units( am_I_Root, Input_Opt, State_Met, &
                            State_Chm, 'kg/kg dry', RC, OrigUnit=OrigUnit )
    IF ( RC/= HCO_SUCCESS ) RETURN

    ! Loop over # of tendencies species
    DO I = 1, nSpc 

       ! Skip if not used    
       IF ( ThisTend%SpcInUse(I) <= 0 ) CYCLE

       ! Name of diagnostics 
       DiagnName = 'TEND_' // TRIM(ThisTend%TendName) // '_' //   &
                   TRIM( State_Chm%SpcData(I)%Info%Name )

       ! Get pointer to 3D array, define time interval
       Ptr3D => ThisTend%Tendency(I)%Arr 

       ! Calculate tendency in v/v/s
       IF ( ZeroTend ) THEN
          Tend = 0.0_f4
       ELSE

          ! TMP is the current concentration in v/v
          ! Convert from kg/kg dry to v/v dry
          TMP = State_Chm%Species(:,:,:,I)                   &
              * ( AIRMW / State_Chm%SpcData(I)%Info%emMW_g )
!!!              / State_Met%AD(:,:,:) ) &

          ! Calculate tendency 
          Tend = ( TMP - Ptr3D ) / REAL(DT,f4)
       ENDIF

#if !defined( ESMF_ )
       ! Update diagnostics array
       CALL Diagn_Update( am_I_Root, HcoState, cName=DiagnName, &
               Array3D=Tend, COL=Input_Opt%DIAG_COLLECTION, RC=RC )
                          
       IF ( RC /= HCO_SUCCESS ) THEN 
          WRITE(ErrMsg,*) 'Error in updating diagnostics with ID ', cID
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
#endif

       !! Reset values 
       !Ptr3D = 0.0_fp
       ! Update values in the internal array to current tendency
       Ptr3D =  Tend
       Ptr3D => NULL()

    ENDDO !I

    ! Update stage 
    ThisTend%Stage = 2

    ! Convert tracers back to original unit 
    CALL Convert_Spc_Units( am_I_Root, Input_Opt, State_Met, &
                            State_Chm, OrigUnit,  RC )
    IF ( RC/= HCO_SUCCESS ) RETURN

    ! Cleanup
    ThisTend => NULL()

  END SUBROUTINE Tend_Stage2
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Tend_Get
!
! !DESCRIPTION: Subroutine Tend_Get returns the current tendency for the 
! given species and tendency class. 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Tend_Get( am_I_Root, Input_Opt, TendName, SpcID, Stage, Tend, RC ) 
!
! !USES:
!
    USE CMN_SIZE_MOD,      ONLY : IIPAR, JJPAR, LLPAR
!
! !INPUT PARAMETERS:
!
    LOGICAL,          INTENT(IN   ) :: am_I_Root   ! Are we on the root CPU?
    TYPE(OptInput),   INTENT(IN   ) :: Input_Opt   ! Input opts
    CHARACTER(LEN=*), INTENT(IN   ) :: TendName    ! tendency name 
    INTEGER,          INTENT(IN   ) :: SpcID       ! Species ID 
!
! !OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(  OUT) :: Stage       ! Stage of tendency:
                                                   ! 0=does not exist; 1=stage 1; 2=stage 2 
    REAL(f4),         POINTER       :: Tend(:,:,:) ! Tendency array
    INTEGER,          INTENT(OUT)   :: RC          ! Failure or success
!
! !REVISION HISTORY: 
!  14 Jul 2015 - C. Keller   - Initial version 
!  26 Oct 2015 - C. Keller   - Update to include tendency classes
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    LOGICAL                  :: FOUND
    INTEGER                  :: cID, I

    ! Pointers
    TYPE(TendClass), POINTER :: ThisTend

    ! Strings
    CHARACTER(LEN=63)        :: DiagnName
    CHARACTER(LEN=255)       :: MSG
    CHARACTER(LEN=255)       :: LOC = 'TEND_GET (tendencies_mod.F)' 
   
    !=======================================================================
    ! TEND_GET begins here!
    !=======================================================================

    ! Assume successful return
    RC = GC_SUCCESS

    ! Init
    ThisTend => NULL()
    Tend     => NULL()
    Stage    =  0

    ! Find tendency class
    CALL Tend_FindClass( am_I_Root, TendName, FOUND, RC, ThisTend=ThisTend )
    IF ( .NOT. FOUND .OR. .NOT. ASSOCIATED(ThisTend) ) RETURN

    ! Skip if not used    
    IF ( ThisTend%SpcInUse(SpcID) <= 0 ) RETURN 

    ! Get pointer to tendency
    Tend  => ThisTend%Tendency(SpcID)%Arr
    Stage =  ThisTend%Stage

    ! Cleanup
    ThisTend => NULL()

  END SUBROUTINE Tend_Get
!EOC
END MODULE Tendencies_Mod
#endif
