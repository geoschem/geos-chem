!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: sfcvmr_mod.F90
!
! !DESCRIPTION: Module sfcvmr\_mod.F90 is a simple module which forces
!  surface concentrations of relevant species to values read from an external
!  file (via HEMCO). The names of the corresponding HEMCO configuration file
!  entries need to be composed of the below defined prefix and the species
!  name, e.g.:
!
!  * SfcVMR_CH3Cl  /discover/nobackup/cakelle2/data/SfcVMR/sfcvmr.nc CH3Cl  2000/1/1/0 C xy ppbv * - 1 1
!  * SfcVMR_CH2Cl2 /discover/nobackup/cakelle2/data/SfcVMR/sfcvmr.nc CH2Cl2 2000/1/1/0 C xy ppbv * - 1 1
!  * SfcVMR_CHCl3  /discover/nobackup/cakelle2/data/SfcVMR/sfcvmr.nc CHCl3  2000/1/1/0 C xy ppbv * - 1 1
!
!  The concentrations in the file are expected to be in units of ppbv.
!  It is also possible to apply scale factors to these fields, e.g. (to scale surface concentrations by 2):
!  * SfcVMR_CH3Cl  /discover/nobackup/cakelle2/data/SfcVMR/sfcvmr.nc CH3Cl  2000/1/1/0 C xy ppbv * 500 1 1
!  ...
!  500 SfcVMR_ScaleFactor 2.0 - - - xy unitless 1
!
!\\
!\\
! !INTERFACE:
!
MODULE SFCVMR_MOD
!
! !USES:
!
    USE PhysConstants       ! Physical constants
    USE inquireMod,    ONLY : findFreeLUN
    USE PRECISION_MOD       ! For GEOS-Chem Precision (fp)
    USE HCO_ERROR_MOD

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
    PRIVATE :: fixSfcVMR_Init
    PUBLIC  :: fixSfcVMR_Run
    PUBLIC  :: fixSfcVMR_Final
!
! !REVISION HISTORY:
!  24 Dec 2016 - S. D. Eastham - Initial version.
!  16 Aug 2019 - C. Keller  - Now read 2D fields with boundary concentrations
!                             via HEMCO. Applicable to any species.
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !PRIVATE TYPES:
!
  TYPE :: SfcMrObj
     CHARACTER(LEN=63)       :: FldName        ! Field name
     INTEGER                 :: SpcID          ! ID in species database
     TYPE(SfcMrObj), POINTER :: Next => NULL() ! Next element in list
  END TYPE SfcMrObj

  ! Heat of linked list with SfcMrObj objects
  TYPE(SfcMrObj), POINTER :: SfcMrHead => NULL()

  ! Field prefix
  CHARACTER(LEN=63), PARAMETER :: Prefix = 'SfcVMR_'

CONTAINS
!EOC
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
! !IROUTINE: FIXSFCVMR_Init
!
! !DESCRIPTION: Subroutine FIXSFCVMR_Init initializes the SfcMr objects.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE fixSfcVMR_Init( am_I_Root, Input_Opt, State_Met, State_Chm, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE CMN_SIZE_MOD
    USE PRECISION_MOD       ! For GEOS-Chem Precision (fp)
    USE ERROR_MOD,          ONLY : ERROR_STOP
    USE Input_Opt_Mod,      ONLY : OptInput
    USE State_Met_Mod,      ONLY : MetState
    USE State_Chm_Mod,      ONLY : ChmState
    USE Species_Mod,        ONLY : Species
    USE HCO_Interface_Mod,  ONLY : HcoState
    USE HCO_Calc_Mod,       ONLY : HCO_EvalFld
!
! !INPUT PARAMETERS:
!
    LOGICAL,          INTENT(IN   )  :: am_I_Root  ! root CPU?
    TYPE(MetState),   INTENT(IN   )  :: State_Met  ! Met state
    TYPE(ChmState),   INTENT(IN   )  :: State_Chm  ! Chemistry state
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(OptInput),   INTENT(INOUT)  :: Input_Opt  ! Input opts
    INTEGER,          INTENT(INOUT)  :: RC         ! Failure or success
!
! !REVISION HISTORY:
!  16 Aug 2019 - C. Keller   - Updated version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Loop indices
    Integer                 :: N
    TYPE(Species), POINTER  :: SpcInfo
    REAL(hp)                :: Arr2D(IIPAR,JJPAR)
    CHARACTER(LEN=63)       :: FldName
    LOGICAL                 :: FOUND
    TYPE(SfcMrObj), POINTER :: iSfcMrObj => NULL()

    !=================================================================
    ! FIXSFCVMR_Init begins here!
    !=================================================================

    ! Verbose
    IF ( am_I_Root ) THEN
       WRITE(6,*) '--- Initialize surface boundary conditions from input file ---'
    ENDIF

    ! Assume success
    RC = GC_SUCCESS

    ! Head of linked list
    SfcMrHead => NULL()

    ! Loop over all species
    DO N = 1, State_Chm%nSpecies
       ! Species information
       SpcInfo => State_Chm%SpcData(N)%Info

       ! Check if file exists
       FldName = TRIM(Prefix)//TRIM(SpcInfo%Name)
       CALL HCO_EvalFld( am_I_Root, HcoState, Trim(FldName), Arr2D, RC, FOUND=FOUND )
       IF ( RC /= HCO_SUCCESS ) THEN
          RC = GC_FAILURE
          EXIT
       ENDIF

       ! Add to linked list if necessary
       IF ( FOUND ) THEN
           ! Must have positive, non-zero MW
           IF ( SpcInfo%emMW_g <= 0.0_fp ) THEN
               WRITE(6,*) 'Cannot use surface boundary condition for species '//TRIM(SpcInfo%Name)//' due to invalid MW'
               RC = GC_FAILURE
               EXIT
           ENDIF
           ! Create new object, add to list
           ALLOCATE(iSfcMrObj)
           iSfcMrObj%SpcID   =  N
           iSfcMrObj%FldName =  TRIM(Prefix)//TRIM(SpcInfo%Name)
           iSfcMrObj%Next    => SfcMrHead
           SfcMrHead         => iSfcMrObj
           IF ( am_I_Root ) THEN
               WRITE(6,*) '--> '//TRIM(SpcInfo%Name)//': Will use prescribed surface boundary conditions from field '//TRIM(iSfcMrObj%FldName)
           ENDIF
           iSfcMrObj => NULL()
       ENDIF
       RC = GC_SUCCESS
    ENDDO

    ! Error check
    IF ( RC == GC_FAILURE ) RETURN

    ! Verbose
    IF ( am_I_Root ) THEN
       WRITE(6,*) '--- Finished initializing surface boundary conditions ---'
    ENDIF

    ! Return with success
    RC = GC_SUCCESS

  END SUBROUTINE fixSfcVMR_Init
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: FIXSFCVMR_Run
!
! !DESCRIPTION: Subroutine FIXSFCVMR_Run fixes the VMR of selected species
! throughout the PBL to observed values.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE fixSfcVMR_Run( am_I_Root, Input_Opt, State_Met, State_Chm, RC )
!
! !USES:
!
    USE CMN_SIZE_MOD
    USE ErrCode_Mod
    USE ERROR_MOD,          ONLY : ERROR_STOP
    USE Input_Opt_Mod,      ONLY : OptInput
    USE State_Met_Mod,      ONLY : MetState
    USE State_Chm_Mod,      ONLY : ChmState
    USE State_Chm_Mod,      ONLY : Ind_
    USE Species_Mod,        ONLY : Species
    USE HCO_Interface_Mod,  ONLY : HcoState
    USE HCO_Calc_Mod,       ONLY : Hco_EvalFld
    USE HCO_Error_Mod,      ONLY : HCO_SUCCESS
    USE TIME_MOD,           ONLY : Get_Month
    USE HCO_INTERFACE_MOD,  ONLY : HcoState
    USE HCO_STATE_MOD,      ONLY : HCO_GetHcoID

    ! Needed for the new CHxCly boundary condition
    USE PBL_MIX_MOD,        ONLY : GET_FRAC_UNDER_PBLTOP
    Use PhysConstants,      Only : AirMW
!
! !INPUT PARAMETERS:
!
    LOGICAL,          INTENT(IN   )  :: am_I_Root  ! root CPU?
    TYPE(MetState),   INTENT(IN   )  :: State_Met  ! Met state
    TYPE(ChmState),   INTENT(INOUT)  :: State_Chm  ! Chemistry State object
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(OptInput),   INTENT(INOUT)  :: Input_Opt  ! Input opts
    INTEGER,          INTENT(INOUT)  :: RC         ! Failure or success
!
! !REVISION HISTORY:
!  27 Aug 2014 - C. Keller   - Initial version
!  16 Jun 2016 - J. Sheng    - Added tracer index retriever
!  20 Jun 2016 - R. Yantosca - Now define species IDs only in the INIT phase
!  16 Aug 2019 - C. Keller   - Updated version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Loop indices
    Integer                     :: I, J, L, MONTH
    ! Species index
    Integer                     :: id_Spc
    ! Saved (main) scalars
    LOGICAL,          SAVE      :: FIRST         = .TRUE.
    ! Strings
    CHARACTER(LEN=255)          :: LOC = 'fixSfcVMR_Run (sfcvmr_mod.F90)'

    ! Pointer to the species array
    Real(fp),    Pointer       :: Spc(:,:,:,:)
    ! Species information
    TYPE(Species), POINTER     :: SpcInfo
    ! Species information
    REAL(hp)                   :: Arr2D(IIPAR,JJPAR)
    ! Linked list
    TYPE(SfcMrObj), POINTER    :: iObj => NULL()

    !=================================================================
    ! FIXSFCVMR_Run begins here!
    !=================================================================

    ! Assume success
    RC        = GC_SUCCESS

    ! Initialize object if needed
    IF ( FIRST ) THEN
       CALL fixSfcVMR_Init( am_I_Root, Input_Opt, State_Met, State_Chm, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       FIRST = .FALSE.
    ENDIF

    ! Get a pointer to the species array
    Spc => State_Chm%Species

    ! Loop over all objects
    iObj => SfcMrHead
    DO WHILE( ASSOCIATED(iObj) )

       ! Get concentration for this species
       CALL HCO_EvalFld( am_I_Root, HcoState, Trim(iObj%FldName), Arr2D, RC )
       IF ( RC /= HCO_SUCCESS ) RETURN

       ! Set mixing ratio in PBL
       SpcInfo => State_Chm%SpcData(iObj%SpcID)%Info
       id_Spc = SpcInfo%ModelID
       IF ( id_Spc > 0 ) THEN

          DO L = 1, LLPAR
          DO J = 1, JJPAR
          DO I = 1, IIPAR
             IF (GET_FRAC_UNDER_PBLTOP(I,J,L)>0e+0_fp) THEN
                Spc(I,J,L,id_Spc) = Arr2d(I,J)*1.0e-9 / ( AIRMW / SpcInfo%emMW_g )
             ENDIF  ! end selection of PBL boxes
          ENDDO
          ENDDO
          ENDDO

       ENDIF

       ! Point to next element in list
       iObj => iObj%Next
    ENDDO

    ! Free pointer
    Spc => NULL()

    ! Return w/ success
    RC = GC_SUCCESS

  END SUBROUTINE fixSfcVMR_Run
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: FIXSFCVMR_Final
!
! !DESCRIPTION: Subroutine FIXSFCVMR_Final cleans up the FixSfcMR linked list.
!\\
!\\
! !INTERFACE:
!
    SUBROUTINE fixSfcVMR_Final( am_I_Root, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE ERROR_MOD,          ONLY : ERROR_STOP
!
! !INPUT PARAMETERS:
!
    LOGICAL,          INTENT(IN   )  :: am_I_Root  ! root CPU?
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(INOUT)  :: RC         ! Failure or success
!
! !REVISION HISTORY:
!  16 Aug 2019 - C. Keller   - Updated version
!EOP
!------------------------------------------------------------------------------
!BOC

    ! Linked list
    TYPE(SfcMrObj), POINTER :: iObj     => NULL()
    TYPE(SfcMrObj), POINTER :: iObjNext => NULL()

    ! Loop over all objects
    iObj => SfcMrHead
    DO WHILE( ASSOCIATED(iObj) )
        iObjNext => iObj%Next
        iObj%Next => NULL()
        DEALLOCATE(iObj)
        iObj => iObjNext
    ENDDO

    ! Return w/ success
    RC = GC_SUCCESS

  END SUBROUTINE fixSfcVMR_Final
!EOC
END MODULE SFCVMR_MOD
