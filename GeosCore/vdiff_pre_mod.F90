!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: vdiff_pre_mod.F90
!
! !DESCRIPTION: Module VDIFF\_PRE\_MOD contains variables used in VDIFF\_MOD.
!\\
!\\
! !INTERFACE: 
!
MODULE VDIFF_PRE_MOD
! 
! !USES:
!
  IMPLICIT NONE
  
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
! 
  PUBLIC :: Init_VDIFF_PRE
  PUBLIC :: Cleanup_VDIFF_PRE
  PUBLIC :: Set_VDIFF_VALUES
!
! !PUBLIC DATA MEMBERS:
!
  LOGICAL, PUBLIC      :: LPRT               ! Passes LPRT to vdiff_mod
  LOGICAL, PUBLIC      :: LTURB              ! Passes LTURB to vdiff_mod
  INTEGER, PUBLIC      :: PCNST              ! Passes N_TRACERS to vdiff_mod
!
! !REVISION HISTORY:
!  01 Jun 2009 - C. Carouge & J. Lin - Initial version  
!  07 Oct 2009 - R. Yantosca         - Added CVS Id tag  
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
!  24 Jun 2014 - R. Yantosca - Now add PCNST as a module variable
!  24 Jun 2014 - R. Yantosca - Add routine SET_VDIFF_VALUES, since we need to
!                              pass N_TRACERS, LPRT, LTURB to vdiff_mod.F90
!                              now that logical_mod.F, tracer_mod.F are gone.
!  24 Jun 2014 - R. Yantosca - Renamed to vdiff_pre_mod.F90
!  24 Nov 2014 - M. Yannetti - Added PRECISION_MOD
!  30 Jun 2016 - M. Sulprizio- Remove NCS, NDRYDEP, ND15, ND44, IIPAR, JJPAR
!                              since they are not used anymore
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
! !IROUTINE: init_vdiff_pre
!
! !DESCRIPTION: Subroutine INIT\_VDIFF\_PRE allocates all module arrays.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Init_VDIFF_PRE( am_I_Root, RC )
!
! !USES:
!
    USE ErrCode_Mod
!
! !INPUT PARAMETERS: 
!
    LOGICAL, INTENT(IN)  :: am_I_Root   ! Are we on the root CPU?
!
! !OUTPUT PARAMETERS:
!
    INTEGER, INTENT(OUT) :: RC          ! Success or failure?
!
! !REMARKS:
!  Need to add error-checking on the allocation statements, so that we
!  exit the code upon error.
! 
! !REVISION HISTORY: 
!  19 Nov 2012 - R. Yantosca - Added ProTeX headers
!  24 Jun 2014 - R. Yantosca - Now accept Input_Opt via the arg list
!  24 Jun 2014 - R. Yantosca - Now allocate EMIS_SAVE to the # of tracers
!                              in the simulation (i.e. INIT_OPT)
!  22 May 2015 - R. Yantosca - Remove variables made obsolete by HEMCO
!EOP
!------------------------------------------------------------------------------
!BOC
    ! Assume success
    RC = GC_SUCCESS
      
  END SUBROUTINE Init_VDIFF_PRE
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: cleanup_vdiff_pre
!
! !DESCRIPTION: Subroutine CLEANUP\_VDIFF\_PRE deallocates all module arrays.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Cleanup_VDIFF_PRE( am_I_Root, RC )
!
! !USES:
!
    USE ErrCode_Mod

    IMPLICIT NONE
!
! !INPUT PARAMETERS: 
!
    LOGICAL, INTENT(IN)  :: am_I_Root   ! Are we on the root CPU?
!
! !OUTPUT PARAMETERS:
!
    INTEGER, INTENT(OUT) :: RC          ! Success or failure?
!
! !REMARKS:
!
! !REVISION HISTORY: 
!  19 Nov 2012 - R. Yantosca - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
    ! Assume success
    RC = GC_SUCCESS

  END SUBROUTINE Cleanup_VDIFF_PRE
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: set_vdiff_values
!
! !DESCRIPTION: Subroutine SET\_VDIFF\_VALUES initializes the PCNST value, which
!  is the number of advected species.  This is needed in vdiff\_mod.F90.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Set_VDIFF_VALUES( am_I_Root, Input_Opt, State_Chm, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE Input_Opt_Mod,      ONLY : OptInput
    USE State_Chm_Mod,      ONLY : ChmState
!
! !INPUT PARAMETERS:
!
    LOGICAL,        INTENT(IN)  :: am_I_Root   ! Are we on the root CPU?
    TYPE(OptInput), INTENT(IN)  :: Input_Opt   ! Input Options object
    TYPE(ChmState), INTENT(IN)  :: State_Chm   ! Chemistry State object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT) :: RC          ! Success or failure?
!
! !REMARKS:
!  This routine has to be called after routine READ_INPUT_FILE.
!
! !REVISION HISTORY: 
!  24 Jun 2014 - R. Yantosca - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
    ! Assume success
    RC    = GC_SUCCESS

    !=====================================================================
    ! The following fields of Input_Opt have to be set in module
    ! variables.  This will allow us to pass these to vdiff_mod.F90,
    ! now that logical_mod.F and tracer_mod.F have been retired.
    !=====================================================================

    ! Number of advected species
    PCNST = State_Chm%nAdvect

    ! Debug print?
    LPRT  = ( Input_Opt%LPRT .and. am_I_Root )
      
    ! Use PBL mixing?
    LTURB = Input_Opt%LTURB

  END SUBROUTINE SET_VDIFF_VALUES
!EOC
END MODULE VDIFF_PRE_MOD

