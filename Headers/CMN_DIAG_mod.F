#ifdef BPCH_DIAG
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: CMN_DIAG_mod.F
!
! !DESCRIPTION: Module CMN\_DIAG\_mod contains size parameters and global
!  variables for the GEOS-Chem diagnostic arrays.  This is mostly historical
!  baggage.
!\\
!\\
! !INTERFACE:
!
      MODULE CMN_DIAG_MOD
!
! !USES:
!
      USE CMN_SIZE_MOD

      IMPLICIT NONE
      PUBLIC
!
! !DEFINED PARAMETERS:
!
      INTEGER            :: ND40
      INTEGER, PARAMETER :: MAX_TRACER = 600      ! Large placeholder value

#if defined(TOMAS)
      !=================================================================
      ! Settings for TOMAS aerosol microphysics (win, bmy, 1/22/10)
      !=================================================================
      INTEGER, PARAMETER :: PD59=TOMASBIN*TOMASSPEC
      INTEGER, PARAMETER :: PD60=TOMASBIN*TOMASSPEC
      INTEGER, PARAMETER :: PD61=2
      INTEGER            :: PD65
      INTEGER, PARAMETER :: MAX_DIAG   = 80

#elif defined(RRTMG)
      !=================================================================
      ! Settings for RRTMG radiative transfer model
      !=================================================================

      !number of rad flux and optics output types
      !8 flux and 3*3=9 optics
      INTEGER, PARAMETER :: PD72R=17

      !total number of possible rad outputs (types*specs)
      !there are 11 possible flux output 'species' but
      !only 8 possible optics output 'species'
      !for simplicity we take the largest and put up with
      !some redundancy (should be 88+72=160)
      INTEGER, PARAMETER :: PD72=187 ! Radiation (Ridley 10/2012)
      INTEGER, PARAMETER :: MAX_DIAG   = 187

#else
      !=================================================================
      ! For non-TOMAS and non-RRTMG with BPCH_DIAG=y
      !=================================================================
      INTEGER, PARAMETER :: MAX_DIAG   = 80

#endif
!
! !PUBLIC DATA MEMBERS:
!
      !=================================================================
      ! NDxx diagnostic flags
      !=================================================================
#ifdef TOMAS
      INTEGER :: LD06, LD44, LD65, LD59, LD60, LD61
      INTEGER :: ND06, ND44, ND65, ND59, ND60, ND61
#endif

#ifdef RRTMG
      INTEGER :: ND72
#endif

      !=================================================================
      ! Variables for printing out selected tracers in diagnostic output
      !=================================================================
      INTEGER :: TINDEX(MAX_DIAG,MAX_TRACER)
      INTEGER :: TCOUNT(MAX_DIAG)
      INTEGER :: TMAX(MAX_DIAG)
!
! !REMARKS:
!  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!  %%%  NOTE: THIS MODULE WILL BE A STUB UNLESS GEOS-Chem IS COMPILED    %%%
!  %%%  WITH THE BPCH_DIAG=y OPTION. (bmy, 10/4/19)                      %%%
!  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! !REVISION HISTORY:
!  See the Git history with the gitk browser!
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
! !IROUTINE: Init_Cmn_Diag
!
! !DESCRIPTION: Subroutine INIT\_CMN\_DIAG initializes quantities based on
!  the grid-independent size parameters.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE Init_CMN_DIAG( am_I_Root, State_Grid, RC )
!
! !USES:
!
      USE ErrCode_Mod
      USE State_Grid_Mod, ONLY : GrdState
!
! !INPUT PARAMETERS:
!
      LOGICAL,        INTENT(IN)  :: am_I_Root   ! Are we on the root CPU?
      TYPE(GrdState), INTENT(IN)  :: State_Grid  ! Grid State object
!
! !OUTPUT PARAMETERS:
!
      INTEGER,        INTENT(OUT) :: RC          ! Success or failure?
!
! !REVISION HISTORY:
!  19 Nov 2012 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
      ! Assume success
      RC   = GC_SUCCESS

#ifdef TOMAS
      PD65 = State_Grid%NZ * MAXFAM
#endif

      END SUBROUTINE Init_CMN_DIAG
!EOC
      END MODULE CMN_DIAG_MOD
#endif
