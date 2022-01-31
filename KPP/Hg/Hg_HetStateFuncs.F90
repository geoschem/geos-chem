!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: Hg_HetStateFuncs.F90
!
! !DESCRIPTION: Module for initializing the HetState object, which passes
!  arguments from GEOS-Chem to the heterogeneous chemistry routines.
!\\
!\\
! !INTERFACE:

MODULE Hg_HetStateFuncs
!
! !USES:
!
  USE GcKpp_Precision

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC :: Hg_SetStateHet

! !REVISION HISTORY:
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Hg_SetStateHet
!
! !DESCRIPTION: Initializes the State_Het object with gridbox values passed
!  from Hg_mod.  These values are used in the heterogenous chemistry
!  reaction rate computations.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Hg_SetStateHet( I,          J,         L,                       &
                             Input_Opt,  State_Chm, State_Met,               &
                             fracOrgAer, H,         RC                      )
!
! !USES:
!
    USE ErrCode_Mod
    USE GcKpp_Global
    USE Input_Opt_Mod,    ONLY : OptInput
    USE rateLawUtilFuncs, ONLY : Cld_Params
    USE State_Chm_Mod,    ONLY : ChmState
    USE State_Met_Mod,    ONLY : MetState
!
! !INPUT PARAMETERS:
!
    INTEGER,        INTENT(IN)    :: I           ! Lon index
    INTEGER,        INTENT(IN)    :: J           ! Lat index
    INTEGER,        INTENT(IN)    :: L           ! Level index
    TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
    TYPE(ChmState), INTENT(IN)    :: State_Chm   ! Chemistry State object
    TYPE(MetState), INTENT(IN)    :: State_Met   ! Meteorology State object
    REAL(sp),       INTENT(IN)    :: fracOrgAer  ! Frac forming organic HgIIP
!
! INPUT/OUTPUT PARAMETERS:
!
    TYPE(HetState), INTENT(INOUT) :: H           ! Hetchem State object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
!EOP
!------------------------------------------------------------------------------
!BOC

    ! Initialization
    RC = GC_SUCCESS

    !========================================================================
    ! Populate fields of the HetState object in gckpp_Global
    !========================================================================

    ! Identify a box for debug printout within rate-law functions
    H%debugBox     = .FALSE.

    ! Identify if this grid box is cloudy
    H%cloudBox     = ( ( State_Met%InTroposphere(I,J,L)              )  .or. &
                       ( State_Met%T(I,J,L)             >= 258.0_dp  )  .or. &
                       ( State_Met%CLDF(I,J,L)          >= 1.0e-3_dp )      )

    ! Meteorology-related quantities
    H%CldFr        = MIN( MAX( State_Met%CLDF(I,J,L), 0.0_dp ), 1.0_dp )
    H%ClearFr      = 1.0_dp - State_Het%CldFr
    H%QICE         = State_Met%QI(I,J,L)
    H%QLIQ         = State_Met%QL(I,J,L)
    H%vAir         = State_Met%AIRVOL(I,J,L) * 1.0e6_dp

    ! Fraction of species forming organic or inorganic HgII aerosol
    H%fracOrgAer   = fracOrgAer
    H%fracInorgAer = 1.0_dp - fracOrgAer

    ! Cloud fields
    CALL Cld_Params(                                                         &
         AD        = State_Met%AD(I,J,L),                                    &
         CLDF      = State_Met%CLDF(I,J,L),                                  &
         FRLAND    = State_Met%FRLAND(I,J),                                  &
         FROCEAN   = State_Met%FROCEAN(I,J),                                 &
         QI        = State_Met%QI(I,J,L),                                    &
         QL        = State_Met%QL(I,J,L),                                    &
         T         = State_Met%T(I,J,L),                                     &
         H         = H                                                      )

  END SUBROUTINE Hg_SetStateHet
!EOC
END MODULE Hg_HetStateFuncs
