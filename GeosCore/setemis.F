!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: setemis
!
! !DESCRIPTION: Subroutine SETEMIS places emissions from the 3-D 
! emission fields in State\_Chm%Trac\_Tend into the REMIS array (used
! for SMVGEAR II chemistry). Emissions are converted from kg/m2/s to
! molec/cm3/s.
!\\
!\\
! If the full-mixing PBL scheme is used, the entire emission arrays
! of State\_Chm%Trac\_Tend are passed to REMIS. For the non-local PBL
! scheme, only the emissions above the PBL are passed since all
! emissions w/in the PBL are already considered in the non-local PBL
! mixing routine (vdiff\_mod.F90). 
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE SETEMIS( am_I_Root, Input_Opt, 
     &                    State_Met, State_Chm, RC )
!
! !USES:
!
      USE GIGC_ErrCode_Mod
      USE COMODE_MOD,         ONLY : REMIS
      USE GIGC_State_Met_Mod, ONLY : MetState
      USE GIGC_State_Chm_Mod, ONLY : ChmState
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
      USE PRECISION_MOD
!
! !INPUT PARAMETERS:
!
      LOGICAL,        INTENT(IN)  :: am_I_Root   ! Are we on the root CPU
      TYPE(OptInput), INTENT(IN)  :: Input_Opt   ! Input Options object
      TYPE(MetState), INTENT(IN)  :: State_Met   ! Meteorology State object
      TYPE(ChmState), INTENT(IN)  :: State_Chm   ! Chemistry State object
!
! !OUTPUT PARAMETERS:
!
      INTEGER,        INTENT(OUT) :: RC          ! Success or failure

! !REVISION HISTORY:
!  13 Apr 2012 - C. Keller   - Initial Version
!  31 Jan 2013 - C. Keller   - Now pass chemistry and met state as argument
!  22 Aug 2013 - C. Keller   - Now pass Trac_Tend2 to REMIS
!  21 Oct 2014 - C. Keller   - Now use GetHcoVal instead of Trac_Tend
!  20 Nov 2014 - M. Yannetti - Added PRECISION_MOD
!  22 May 2015 - R. Yantosca - Remove variables made obsolete by HEMCO
!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!

      !=================================================================
      ! SETEMIS begins here!
      !=================================================================

      ! Assume success
      RC        = GIGC_SUCCESS

      ! Reset 
      REMIS     = 0e+0_fp

      ! Emissions are now completely treated outside of the chemistry
      ! solver. Therefore return here
      RETURN

      END SUBROUTINE SETEMIS
!EOC
