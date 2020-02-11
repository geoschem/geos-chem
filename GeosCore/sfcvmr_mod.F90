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
! * SfcVMR_CH3Cl  $ROOT/CMIP6/v2019-09//LIVE/CMIP6_GHG_surface_VMR_1750_2014_for_CH3Cl.nc CH3Cl   1750-2014/1-12/1/0 C xy ppbv * 801 1 1
! * SfcVMR_CH2Cl2 $ROOT/CMIP6/v2019-09//LIVE/CMIP6_GHG_surface_VMR_1750_2014_for_CH2Cl2.nc CH2Cl2 1750-2014/1-12/1/0 C xy ppbv * 801 1 1
! * SfcVMR_CHCl3  $ROOT/CMIP6/v2019-09//LIVE/CMIP6_GHG_surface_VMR_1750_2014_for_CHCl3.nc CHCl3   1750-2014/1-12/1/0 C xy ppbv * 801 1 1
! * SfcVMR_CH3Br  /LIVE/CMIP6_GHG_surface_VMR_1750_2014_for_CH3Br.nc CH3Br   1750-2014/1-12/1/0 C xy ppbv * 801 1 1
!
!  The concentrations in the file are expected to be in units of ppbv.
!  It is also possible to apply scale factors to these fields, e.g. (to scale surface concentrations by 2):
!  * SfcVMR_CH3Cl  $ROOT/CMIP6/v2019-09//LIVE/CMIP6_GHG_surface_VMR_1750_2014_for_CH3Cl.nc CH3Cl   1750-2014/1-12/1/0 C xy ppbv * 801 1 1

!  ...
!  # Scale the CMIP6 values in pptv to ppbv
!  802 SfcVMR_ScaleFactor 0.001 - - - xy unitless 1
!
!\\
!\\
! !INTERFACE:
!
MODULE SfcVmr_Mod
!
! !USES:
!
  USE PhysConstants       ! Physical constants
  USE Precision_Mod       ! For GEOS-Chem Precision (fp)
  USE HCO_Error_Mod       ! HEMCO error handling variables & functions

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC  :: FixSfcVmr_Run
  PUBLIC  :: FixSfcVmr_Final
!
! !PRIVATE MEMBER FUNCTIONS:
!
  PRIVATE :: FixSfcVmr_Init
!
! !REVISION HISTORY:
!  24 Dec 2016 - S. D. Eastham - Initial version.
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !PRIVATE TYPES:
!
  ! Linked list type
  TYPE :: SfcMrObj
     CHARACTER(LEN=63)         :: FldName        ! Field name
     INTEGER                   :: SpcID          ! ID in species database
     TYPE(SfcMrObj), POINTER   :: Next           ! Next element in list
  END TYPE SfcMrObj

  ! Heat of linked list with SfcMrObj objects
  TYPE(SfcMrObj),    POINTER   :: SfcMrHead => NULL()

  ! Field prefix
  CHARACTER(LEN=63), PARAMETER :: Prefix = 'SfcVMR_'

CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: FixSfcVmr_Init
!
! !DESCRIPTION: Subroutine FixSfcVmr_Init initializes the SfcMr objects.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE FixSfcVmr_Init( Input_Opt, State_Chm, State_Grid, State_Met, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE Input_Opt_Mod,      ONLY : OptInput
    USE State_Met_Mod,      ONLY : MetState
    USE State_Chm_Mod,      ONLY : ChmState
    USE State_Grid_Mod,     ONLY : GrdState
    USE Species_Mod,        ONLY : Species
    USE HCO_Interface_Mod,  ONLY : HcoState
    USE HCO_Calc_Mod,       ONLY : HCO_EvalFld
!
! !INPUT PARAMETERS:
!
    TYPE(MetState),   INTENT(IN)     :: State_Met   ! Met state
    TYPE(GrdState),   INTENT(IN)     :: State_Grid  ! Grid State object
    TYPE(ChmState),   INTENT(IN)     :: State_Chm   ! Chemistry state
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(OptInput),   INTENT(INOUT)  :: Input_Opt   ! Input opts
!
! !OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(OUT)    :: RC          ! Failure or success
!
! !REVISION HISTORY:
!  16 Aug 2019 - C. Keller   - Updated version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    LOGICAL                 :: FOUND
    INTEGER                 :: N

    ! Strings
    CHARACTER(LEN=63)       :: FldName
    CHARACTER(LEN=255)      :: ErrMsg
    CHARACTER(LEN=255)      :: ThisLoc

    ! Arrays
    REAL(hp)                :: Arr2D(State_Grid%NX,State_Grid%NY)

    ! Pointers
    TYPE(Species),  POINTER :: SpcInfo
    TYPE(SfcMrObj), POINTER :: iSfcMrObj

    !=================================================================
    ! FIXSFCVMR_Init begins here!
    !=================================================================

    ! Initialize
    RC        = HCO_SUCCESS
    ErrMsg    = ''
    ThisLoc   = ' --> at fixSfcVMR_Init (in module GeosCore/sfcvmr_mod.F90)'
    iSfcMrObj => NULL()
    SpcInfo   => NULL()

    ! Verbose
    IF ( Input_Opt%amIRoot ) THEN
       WRITE( 6, 100 ) 
 100   FORMAT('--- Initialize surface boundary conditions from input file ---')
    ENDIF

    ! Head of linked list
    SfcMrHead => NULL()

    ! Loop over all species
    DO N = 1, State_Chm%nSpecies
       ! Species information
       SpcInfo => State_Chm%SpcData(N)%Info

       ! Check if file exists
       FldName = TRIM( Prefix ) // TRIM( SpcInfo%Name )
       CALL HCO_EvalFld( HcoState, TRIM(FldName), Arr2D, RC, FOUND=FOUND )
       IF ( RC /= HCO_SUCCESS ) THEN
          ErrMsg = 'Could not find field : ' // TRIM( FldName )
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

       ! Add to linked list if necessary
       IF ( FOUND ) THEN

           ! Must have positive, non-zero MW
           IF ( SpcInfo%emMW_g <= 0.0_fp ) THEN
              ErrMsg = 'Cannot use surface boundary condition for species '  &
                     // TRIM(SpcInfo%Name) // ' due to invalid MW!'
              CALL GC_Error( ErrMsg, RC, ThisLoc )
              RETURN
           ENDIF

           ! Create new object, add to list
           ALLOCATE( iSfcMrObj, STAT=RC )
           CALL GC_CheckVar( 'sfcvmr_mod.F90:iSfcMrObj', 0, RC )
           IF ( RC /= HCO_SUCCESS ) RETURN

           iSfcMrObj%SpcID   =  N
           iSfcMrObj%FldName =  TRIM(Prefix)//TRIM(SpcInfo%Name)
           iSfcMrObj%Next    => SfcMrHead
           SfcMrHead         => iSfcMrObj
           IF ( Input_Opt%amIRoot ) THEN
              WRITE( 6, 110 ) TRIM( SpcInfo%Name ), TRIM( iSfcMrObj%FldName )
 110          FORMAT( '--> ', a, ' will use prescribed surface boundary ',   &
                      'conditions from field ', a )
           ENDIF

           ! Free the pointer
           iSfcMrObj => NULL()
       ENDIF

       ! Indicate success
       RC = HCO_SUCCESS
    ENDDO

    ! Exit if unsuccessful
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! If successful, print message
    IF ( Input_Opt%amIRoot ) THEN
       WRITE( 6, 120 )
 120   FORMAT( '--- Finished initializing surface boundary conditions ---' )
    ENDIF

  END SUBROUTINE fixSfcVMR_Init
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: FixSfcVmr_Run
!
! !DESCRIPTION: Subroutine FIXSFCVMR_Run fixes the VMR of selected species
! throughout the PBL to observed values.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE FixSfcVmr_Run( Input_Opt, State_Chm, State_Grid, State_Met, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE Input_Opt_Mod,      ONLY : OptInput
    USE State_Met_Mod,      ONLY : MetState
    USE State_Chm_Mod,      ONLY : ChmState
    USE State_Grid_Mod,     ONLY : GrdState
    USE State_Chm_Mod,      ONLY : Ind_
    USE Species_Mod,        ONLY : Species
    USE HCO_Interface_Mod,  ONLY : HcoState
    USE HCO_Calc_Mod,       ONLY : Hco_EvalFld
    USE HCO_Error_Mod,      ONLY : HCO_SUCCESS
    USE TIME_MOD,           ONLY : Get_Month
    USE HCO_INTERFACE_MOD,  ONLY : HcoState
    USE HCO_STATE_MOD,      ONLY : HCO_GetHcoID

    ! Needed for the new CHxCly boundary condition
    Use PhysConstants,      Only : AirMW
!
! !INPUT PARAMETERS:
!
    TYPE(GrdState),   INTENT(IN)     :: State_Grid  ! Grid State object
    TYPE(MetState),   INTENT(IN)     :: State_Met   ! Met State object
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ChmState),   INTENT(INOUT)  :: State_Chm   ! Chemistry State object
    TYPE(OptInput),   INTENT(INOUT)  :: Input_Opt   ! Input Options object
    INTEGER,          INTENT(INOUT)  :: RC          ! Failure or success
!
! !REVISION HISTORY:
!  27 Aug 2014 - C. Keller   - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! SAVEd scalars
    LOGICAL,        SAVE     :: FIRST = .TRUE.

    ! Scalars
    INTEGER                  :: I, J, L, MONTH
    INTEGER                  :: id_Spc

    ! Strings
    CHARACTER(LEN=255)       :: ErrMsg
    CHARACTER(LEN=255)       :: ThisLoc

    ! Arrays
    REAL(hp)                 :: Arr2D(State_Grid%NX,State_Grid%NY)

    ! Linked list
    Real(fp),       POINTER  :: Spc(:,:,:,:)   ! Ptr to species array
    TYPE(Species),  POINTER  :: SpcInfo        ! Ptr to species database
    TYPE(SfcMrObj), POINTER  :: iObj           ! Linked list

    !=======================================================================
    ! FIXSFCVMR_Run begins here!
    !=======================================================================

    ! Assume success
    RC        = HCO_SUCCESS
    ErrMsg    = ''
    ThisLoc   = ' -> at FixSfcVmrRun (in module GeosCore/sfcvmr_mod.F90)'

    ! Initialize object if needed
    IF ( FIRST ) THEN
       CALL FixSfcVMR_Init( Input_Opt, State_Chm, State_Grid, State_Met, RC )
       IF ( RC /= HCO_SUCCESS ) THEN
          ErrMsg = 'Error encountered in routine "FixSfcVmrInit"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
       
       ! Reset first-time flag
       FIRST = .FALSE.
    ENDIF

    ! Get a pointer to the species array
    Spc => State_Chm%Species

    ! Loop over all objects
    iObj => SfcMrHead
    DO WHILE( ASSOCIATED( iObj ) )

       ! Get concentration for this species
       CALL HCO_EvalFld( HcoState, Trim(iObj%FldName), Arr2D, RC )
       IF ( RC /= HCO_SUCCESS ) THEN
          ErrMsg = 'Could not get surface VMR for species: '//               &
                   TRIM( iObj%FldName ) // '!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

       ! Set mixing ratio in PBL
       SpcInfo => State_Chm%SpcData(iObj%SpcID)%Info
       id_Spc = SpcInfo%ModelID
       IF ( id_Spc > 0 ) THEN

          DO L = 1, State_Grid%NZ
          DO J = 1, State_Grid%NY
          DO I = 1, State_Grid%NX
             IF ( State_Met%F_UNDER_PBLTOP(I,J,L) > 0.0_fp ) THEN
                Spc(I,J,L,id_Spc) = ( Arr2d(I,J) * 1.0e-9_fp      )          &
                                  / ( AIRMW      / SpcInfo%emMW_g )
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

  END SUBROUTINE FixSfcVmr_Run
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: FixSfcVmr_Final
!
! !DESCRIPTION: Subroutine FIXSFCVMR_Final cleans up the FixSfcMR linked list.
!\\
!\\
! !INTERFACE:
!
    SUBROUTINE FixSfcVmr_Final( RC )
!
! !USES:
!
    USE ErrCode_Mod
!
! !OUTPUT PARAMETERS:
!
    INTEGER, INTENT(OUT) :: RC          ! Failure or success
!
! !REVISION HISTORY:
!  16 Aug 2019 - C. Keller   - Updated version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Pointers
    TYPE(SfcMrObj), POINTER :: iObj
    TYPE(SfcMrObj), POINTER :: iObjNext

    ! Initialize
    RC       = GC_SUCCESS
    iObj     => NULL()
    iObjNext => NULL()

    ! Loop over all objects and deallocate
    iObj => SfcMrHead
    DO WHILE( ASSOCIATED( iObj ) )
        iObjNext  => iObj%Next
        iObj%Next => NULL()
        IF ( ASSOCIATED( iObj ) ) DEALLOCATE(iObj)
        iObj => iObjNext
    ENDDO

  END SUBROUTINE FixSfcVmr_Final
!EOC
END MODULE SfcVmr_Mod
