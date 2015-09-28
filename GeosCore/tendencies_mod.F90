!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: tendencies_mod.F90
!
! !DESCRIPTION: Module tendencies\_mod.F90 is a module to define, archive and
! write tracer tendencies. 
!
! !INTERFACE:
!
MODULE Tendencies_Mod 
!
! !USES:
!
  USE Precision_Mod
  USE HCO_Error_Mod
  USE HCO_Diagn_Mod
  USE GIGC_ErrCode_Mod
  USE Error_Mod,          ONLY : Error_Stop
  USE GIGC_Input_Opt_Mod, ONLY : OptInput
  USE GIGC_State_Met_Mod, ONLY : MetState
  USE GIGC_State_Chm_Mod, ONLY : ChmState

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!

  ! All of the following is currently only active in development mode:
#if defined( DEVEL )
  PUBLIC :: Tendencies_Init   
  PUBLIC :: Tendencies_Stage1
  PUBLIC :: Tendencies_Stage2
  PUBLIC :: Tendencies_Cleanup
!
! !PRIVATE MEMBER FUNCTIONS:
!
  REAL(f4), ALLOCATABLE, TARGET :: TendAdv(:,:,:,:)
  REAL(f4), ALLOCATABLE, TARGET :: TendCnv(:,:,:,:)
  REAL(f4), ALLOCATABLE, TARGET :: TendChm(:,:,:,:)
  REAL(f4), ALLOCATABLE, TARGET :: TendTrb(:,:,:,:)
  REAL(f4), ALLOCATABLE, TARGET :: TendSfc(:,:,:,:)
  REAL(f4), ALLOCATABLE, TARGET :: TendWet(:,:,:,:)
  INTEGER,  ALLOCATABLE         :: TendIds(:)
!
! !DEFINED PARAMETERS:
!
  ! Tendencies species - hardcoded for now
  INTEGER, PARAMETER            :: NTEND          = 1
  CHARACTER(LEN=15), PARAMETER  :: TENDSPC(NTEND) = (/ 'O3' /) 

  ! Switches - hardcoded for now
  LOGICAL, PARAMETER            :: DO_ADV = .TRUE.
  LOGICAL, PARAMETER            :: DO_CNV = .TRUE.
  LOGICAL, PARAMETER            :: DO_CHM = .TRUE.
  LOGICAL, PARAMETER            :: DO_TRB = .TRUE.
  LOGICAL, PARAMETER            :: DO_SFC = .TRUE.
  LOGICAL, PARAMETER            :: DO_WET = .TRUE.
!
! !REVISION HISTORY:
!  14 Jul 2015 - C. Keller   - Initial version. 
!EOP
!------------------------------------------------------------------------------
!BOC
CONTAINS
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Tendencies_Init 
!
! !DESCRIPTION: Subroutine Tendencies\_Init initializes the tendencies 
! diagnostics. 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Tendencies_Init( am_I_Root, Input_Opt, State_Met, State_Chm, RC ) 
!
! !USES:
!
    USE CMN_SIZE_MOD,       ONLY : IIPAR, JJPAR, LLPAR
    USE GIGC_State_Chm_Mod, ONLY : Get_Indx
!
! !INPUT PARAMETERS:
!
    LOGICAL,          INTENT(IN   ) :: am_I_Root  ! Are we on the root CPU?
    TYPE(OptInput),   INTENT(IN   ) :: Input_Opt  ! Input opts
    TYPE(MetState),   INTENT(IN   ) :: State_Met  ! Met state
    TYPE(ChmState),   INTENT(IN   ) :: State_Chm  ! Chemistry state 
!
! !OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(OUT)   :: RC         ! Failure or success
!
! !REVISION HISTORY: 
!  14 Jul 2015 - C. Keller   - Initial version 
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER            :: AS
    INTEGER            :: cID, Collection, I, N
    LOGICAL            :: SKIP
    CHARACTER(LEN=60)  :: DiagnName
    CHARACTER(LEN=255) :: MSG
    CHARACTER(LEN=255) :: LOC = 'TENDENCIES_INIT (tendencies_mod.F)' 
    
    !=======================================================================
    ! TENDENCIES_INIT begins here!
    !=======================================================================

    ! Assume successful return
    RC = GIGC_SUCCESS

    ! Get diagnostic parameters from the Input_Opt object
    Collection = Input_Opt%DIAG_COLLECTION
    
    ! Initialize arrays
    ALLOCATE( TendIds(NTEND), STAT = AS )
    IF ( AS==0 .AND. DO_ADV ) THEN
       ALLOCATE( TendAdv(IIPAR,JJPAR,LLPAR,NTEND), STAT = AS )
    ENDIF
    IF ( AS==0 .AND. DO_CNV ) THEN
       ALLOCATE( TendCnv(IIPAR,JJPAR,LLPAR,NTEND), STAT = AS )
    ENDIF
    IF ( AS==0 .AND. DO_CHM ) THEN
       ALLOCATE( TendChm(IIPAR,JJPAR,LLPAR,NTEND), STAT = AS )
    ENDIF
    IF ( AS==0 .AND. DO_TRB ) THEN
       ALLOCATE( TendTrb(IIPAR,JJPAR,LLPAR,NTEND), STAT = AS )
    ENDIF
    IF ( AS==0 .AND. DO_SFC ) THEN
       ALLOCATE( TendSfc(IIPAR,JJPAR,LLPAR,NTEND), STAT = AS )
    ENDIF
    IF ( AS==0 .AND. DO_WET ) THEN
       ALLOCATE( TendWet(IIPAR,JJPAR,LLPAR,NTEND), STAT = AS )
    ENDIF
    IF ( AS /= 0 ) THEN
       CALL ERROR_STOP ( 'Tendencies allocation error', LOC )
       RC = GIGC_FAILURE
       RETURN
    ENDIF

    IF( ALLOCATED(TendAdv) ) TendAdv = -999.0_f4 
    IF( ALLOCATED(TendCnv) ) TendCnv = -999.0_f4 
    IF( ALLOCATED(TendChm) ) TendChm = -999.0_f4 
    IF( ALLOCATED(TendTrb) ) TendTrb = -999.0_f4 
    IF( ALLOCATED(TendSfc) ) TendSfc = -999.0_f4 
    IF( ALLOCATED(TendWet) ) TendWet = -999.0_f4 

    ! Loop over # of tendencies species
    DO I = 1, NTEND 
    
       ! Get tendencies species index
       TendIds(I) = Get_Indx( TRIM(TendSpc(I)), Input_Opt%ID_TRACER, Input_Opt%TRACER_NAME )
       IF ( TendIds(I) <= 0 ) THEN
          IF ( am_I_Root ) THEN
             WRITE(*,*) ''
             WRITE(*,*) 'WARNING: cannot write tendencies - this is not a tracer: ', TRIM(TendSpc(I))
             WRITE(*,*) ''
          ENDIF
          CYCLE
       ENDIF 

       !----------------------------------------------------------------
       ! Create containers for tendencies 
       !----------------------------------------------------------------
       DO N = 1, 6

          ! Diagnostic name and unique ID
          SELECT CASE ( N ) 
             CASE ( 1 ) 
                DiagnName = 'TEND_ADV_' // TRIM( Input_Opt%TRACER_NAME( TendIds(I) ) )
                cID       = 731000 + TendIDs(I)
                SKIP      = .NOT. DO_ADV
             CASE ( 2 ) 
                DiagnName = 'TEND_CONV_' // TRIM( Input_Opt%TRACER_NAME( TendIds(I) ) )
                cID       = 732000 + TendIDs(I)
                SKIP      = .NOT. DO_CNV
             CASE ( 3 ) 
                DiagnName = 'TEND_CHEM_' // TRIM( Input_Opt%TRACER_NAME( TendIds(I) ) )
                cID       = 733000 + TendIDs(I)
                SKIP      = .NOT. DO_CHM
             CASE ( 4 ) 
                DiagnName = 'TEND_TURB_' // TRIM( Input_Opt%TRACER_NAME( TendIds(I) ) )
                cID       = 734000 + TendIDs(I)
                SKIP      = .NOT. DO_TRB
             CASE ( 5 ) 
                DiagnName = 'TEND_SURF_' // TRIM( Input_Opt%TRACER_NAME( TendIds(I) ) )
                cID       = 735000 + TendIDs(I)
                SKIP      = .NOT. DO_SFC
             CASE ( 6 ) 
                DiagnName = 'TEND_WETD_' // TRIM( Input_Opt%TRACER_NAME( TendIds(I) ) )
                cID       = 736000 + TendIDs(I)
                SKIP      = .NOT. DO_WET
          END SELECT

          ! Skip if not used
          IF ( SKIP ) CYCLE

          ! Create container
          CALL Diagn_Create( am_I_Root,                     &
                             Col       = Collection,        & 
                             cID       = cID,               &
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
             MSG = 'Cannot create diagnostics: ' // TRIM(DiagnName)
             CALL ERROR_STOP( MSG, LOC ) 
          ENDIF
 
       ENDDO !N
    ENDDO !I

  END SUBROUTINE Tendencies_Init 
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Tendencies_Stage1
!
! !DESCRIPTION: Subroutine Tendencies\_Stage1 archives the current tracer 
! concentrations into the local tendency arrays.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Tendencies_Stage1( am_I_Root, Input_Opt, State_Met, &
                                State_Chm, TendType,  IsInvv,    RC ) 
!
! !USES:
!
!
! !INPUT PARAMETERS:
!
    LOGICAL,          INTENT(IN   ) :: am_I_Root  ! Are we on the root CPU?
    TYPE(OptInput),   INTENT(IN   ) :: Input_Opt  ! Input opts
    TYPE(MetState),   INTENT(IN   ) :: State_Met  ! Met state
    TYPE(ChmState),   INTENT(IN   ) :: State_Chm  ! Chemistry state 
    INTEGER,          INTENT(IN   ) :: TendType   ! 1=Adv; 2=Conv; 3=Chem
    LOGICAL,          INTENT(IN   ) :: IsInvv     ! Is tracer in v/v? 
!
! !OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(OUT)   :: RC         ! Failure or success
!
! !REVISION HISTORY: 
!  14 Jul 2015 - C. Keller   - Initial version 
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER            :: I
    REAL(f4), POINTER  :: Ptr3D(:,:,:) => NULL()
    CHARACTER(LEN=255) :: MSG
    CHARACTER(LEN=255) :: LOC = 'TENDENCIES_STAGE1 (tendencies_mod.F)' 
    
    !=======================================================================
    ! TENDENCIES_STAGE1 begins here!
    !=======================================================================

    ! Assume successful return
    RC = GIGC_SUCCESS

    ! Check if this diagnostics type is indeed defined
    IF ( TendType == 1 .AND. .NOT. DO_ADV ) RETURN
    IF ( TendType == 2 .AND. .NOT. DO_CNV ) RETURN
    IF ( TendType == 3 .AND. .NOT. DO_CHM ) RETURN
    IF ( TendType == 4 .AND. .NOT. DO_TRB ) RETURN
    IF ( TendType == 5 .AND. .NOT. DO_SFC ) RETURN
    IF ( TendType == 6 .AND. .NOT. DO_WET ) RETURN

    ! Loop over # of tendencies species
    DO I = 1, NTEND 
    
       IF ( TendIds(I) <= 0 ) THEN
          CYCLE
       ENDIF 

       ! Get pointer to 3D array to be filled 
       IF ( TendType == 1 ) THEN
          Ptr3D => TendAdv(:,:,:,I) 
       ELSEIF ( TendType == 2 ) THEN
          Ptr3D => TendCnv(:,:,:,I) 
       ELSEIF ( TendType == 3 ) THEN
          Ptr3D => TendChm(:,:,:,I) 
       ELSEIF ( TendType == 4 ) THEN
          Ptr3D => TendTrb(:,:,:,I) 
       ELSEIF ( TendType == 5 ) THEN
          Ptr3D => TendSfc(:,:,:,I) 
       ELSEIF ( TendType == 6 ) THEN
          Ptr3D => TendWet(:,:,:,I) 
       ENDIF

       ! Fill 3D array with current values. Make sure it's in v/v
       IF ( IsInvv ) THEN
          Ptr3D = State_Chm%Tracers(:,:,:,TendIds(I))
       ELSE
          Ptr3D = State_Chm%Tracers(:,:,:,TendIds(I)) &
                * Input_Opt%TCVV(TendIds(I)) / State_Met%AD(:,:,:)
       ENDIF

       ! Cleanup
       Ptr3D => NULL()

    ENDDO !I

  END SUBROUTINE Tendencies_Stage1
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Tendencies_Stage2
!
! !DESCRIPTION: Subroutine Tendencies\_Stage2 calculates the tendencies and 
! writes them into the diagnostics arrays. 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Tendencies_Stage2( am_I_Root, Input_Opt, State_Met, &
                                State_Chm, TendType,  IsInvv, DT, RC ) 
!
! !USES:
!
    USE CMN_SIZE_MOD,      ONLY : IIPAR, JJPAR, LLPAR
!
! !INPUT PARAMETERS:
!
    LOGICAL,          INTENT(IN   ) :: am_I_Root  ! Are we on the root CPU?
    TYPE(OptInput),   INTENT(IN   ) :: Input_Opt  ! Input opts
    TYPE(MetState),   INTENT(IN   ) :: State_Met  ! Met state
    TYPE(ChmState),   INTENT(IN   ) :: State_Chm  ! Chemistry state 
    INTEGER,          INTENT(IN   ) :: TendType   ! 1=Adv; 2=Conv; 3=Chem
    LOGICAL,          INTENT(IN   ) :: IsInvv     ! Is tracer in v/v? 
    REAL(fp),         INTENT(IN   ) :: DT         ! delta time, in seconds 
!
! !OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(OUT)   :: RC         ! Failure or success
!
! !REVISION HISTORY: 
!  14 Jul 2015 - C. Keller   - Initial version 
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER            :: cID, I
    REAL(f4), POINTER  :: Ptr3D(:,:,:) => NULL()
    REAL(f4)           :: TEND(IIPAR,JJPAR,LLPAR)
    CHARACTER(LEN=255) :: MSG
    CHARACTER(LEN=255) :: LOC = 'TENDENCIES_STAGE2 (tendencies_mod.F)' 
   
    !=======================================================================
    ! TENDENCIES_STAGE2 begins here!
    !=======================================================================

    ! Assume successful return
    RC = GIGC_SUCCESS

    ! Check if this diagnostics type is indeed defined
    IF ( TendType == 1 .AND. .NOT. DO_ADV ) RETURN
    IF ( TendType == 2 .AND. .NOT. DO_CNV ) RETURN
    IF ( TendType == 3 .AND. .NOT. DO_CHM ) RETURN
    IF ( TendType == 4 .AND. .NOT. DO_TRB ) RETURN
    IF ( TendType == 5 .AND. .NOT. DO_SFC ) RETURN
    IF ( TendType == 6 .AND. .NOT. DO_WET ) RETURN

    ! Loop over # of tendencies species
    DO I = 1, NTEND 
    
       IF ( TendIds(I) <= 0 ) THEN
          CYCLE
       ENDIF 

       ! Get pointer to 3D array, define time interval, and get diagnostics ID 
       IF ( TendType == 1 ) THEN
          Ptr3D => TendAdv(:,:,:,I)
          cID   = 731000 + TendIds(I)
       ELSEIF ( TendType == 2 ) THEN
          Ptr3D => TendCnv(:,:,:,I) 
          cID   = 732000 + TendIds(I)
       ELSEIF ( TendType == 3 ) THEN
          Ptr3D => TendChm(:,:,:,I) 
          cID   = 733000 + TendIds(I)
       ELSEIF ( TendType == 4 ) THEN
          Ptr3D => TendTrb(:,:,:,I) 
          cID   = 734000 + TendIds(I)
       ELSEIF ( TendType == 5 ) THEN
          Ptr3D => TendSfc(:,:,:,I) 
          cID   = 735000 + TendIds(I)
       ELSEIF ( TendType == 6 ) THEN
          Ptr3D => TendWet(:,:,:,I) 
          cID   = 736000 + TendIds(I)
       ENDIF

       ! Error check: stage 2 must be called after stage 1
       IF ( ALL ( Ptr3D == -999.0_f4 ) ) THEN
          IF ( am_I_Root ) THEN
             WRITE(*,*) 'Warning: cannot do tendency stage 2 - stage 1 not yet called: ', TendType
          ENDIF
          CYCLE
       ENDIF

       ! Calculate tendency in v/v/s
       IF ( IsInvv ) THEN
          Tend = ( State_Chm%Tracers(:,:,:,TendIds(I)) - Ptr3D(:,:,:) ) / DT
       ELSE
          Tend = ( ( State_Chm%Tracers(:,:,:,TendIds(I))                &
                   * Input_Opt%TCVV(TendIds(I)) / State_Met%AD(:,:,:) ) &
                   - Ptr3D(:,:,:) ) / DT
       ENDIF

       ! Update diagnostics array
       CALL Diagn_Update( am_I_Root, cID=cID, Array3D=Tend, &
                          COL=Input_Opt%DIAG_COLLECTION, RC=RC )
       IF ( RC /= HCO_SUCCESS ) THEN 
          WRITE(MSG,*) 'Error in updating diagnostics with ID ', cID
          CALL ERROR_STOP ( MSG, LOC )
          RC = GIGC_FAILURE
          RETURN
       ENDIF

       ! Reset and cleanup
       Ptr3D = -999.0_f4
       Ptr3D => NULL()

    ENDDO !I

  END SUBROUTINE Tendencies_Stage2
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Tendencies_Cleanup
!
! !DESCRIPTION: Subroutine Tendencies\_Cleanup cleans up the tendencies module.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Tendencies_Cleanup
!
! !REVISION HISTORY: 
!  14 Jul 2015 - C. Keller   - Initial version 
!EOP
!------------------------------------------------------------------------------
!BOC

    !=======================================================================
    ! TENDENCIES_CLEANUP begins here!
    !=======================================================================

    IF ( ALLOCATED( TendAdv ) ) DEALLOCATE( TendAdv )
    IF ( ALLOCATED( TendCnv ) ) DEALLOCATE( TendCnv )
    IF ( ALLOCATED( TendChm ) ) DEALLOCATE( TendChm )
    IF ( ALLOCATED( TendTrb ) ) DEALLOCATE( TendTrb )
    IF ( ALLOCATED( TendSfc ) ) DEALLOCATE( TendSfc )
    IF ( ALLOCATED( TendWet ) ) DEALLOCATE( TendWet )
    IF ( ALLOCATED( TendIds ) ) DEALLOCATE( TendIds )

  END SUBROUTINE Tendencies_Cleanup
!EOC
#endif
END MODULE Tendencies_Mod
