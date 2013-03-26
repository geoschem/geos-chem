! $Id: CO_strat_pl.f,v 1.1 2009/09/16 14:06:49 bmy Exp $
#if defined( DEVEL )
      SUBROUTINE CO_STRAT_PL( COPROD, COLOSS, State_Met, State_Chm )
#else
      SUBROUTINE CO_STRAT_PL( COPROD, COLOSS, State_Met )
#endif
!
!******************************************************************************
!  Subroutine CO_STRAT_PL computes net production of CO above the 
!  annual mean tropopause using archived rates for P(CO) and L(CO).
!  (bnd, qli, bmy, 12/9/99, 10/25/05)
!
!  Arguments as Input:
!  ===========================================================================
!  (1 ) COPROD : (REAL*4) Zonally averaged P(CO) in [v/v/s]
!  (2 ) COLOSS : (REAL*4) Zonally averaged L(CO) in [1/s]
!
!  NOTES:
!  (1 ) P(CO) and L(CO) rates were provided by Dylan Jones.  Bob Yantosca
!        has regridded these rates to the GEOS-1 and GEOS-STRAT vertical 
!        levels, for both 2 x 2.5 and 4 x 5 degree resolution.
!  (2 ) CO_STRAT_PL was specially adapted from Bryan Duncan's routine 
!        "CO_strat.f" for use with the simple chemistry module "schem.f".
!  (3 ) Now reference AD from "dao_mod.f".  Now reference IDTCO and IDTCH2O
!        from "tracerid_mod.f". (bmy, 11/6/02)
!  (4 ) Now use function GET_TS_CHEM from "time_mod.f".  Updated comments.
!        (bmy, 2/11/03)
!  (5 ) Now references STT from "tracer_mod.f" (bmy, 7/20/04)
!  (6 ) Now use functions from "tropopause_mod.f" to diagnose whether a box
!        is in the stratosphere or not.  Remove reference to CMN, it's 
!        obsolete. (bmy, 8/22/05)
!  (7 ) Now reference XNUMOL and XNUMOLAIR from "tracer_mod.f" (bmy, 10/25/05)
!  09 Nov 2012 - M. Payer    - Replaced all met field arrays with State_Met
!                              derived type object
!******************************************************************************
!
      ! References to F90 modules 
      USE GIGC_State_Met_Mod, ONLY : MetState
      USE TIME_MOD,           ONLY : GET_TS_CHEM
#if defined( DEVEL )
      USE GIGC_State_Chm_Mod, ONLY : ChmState
#else
      USE TRACER_MOD,         ONLY : STT
#endif
      USE TRACER_MOD,         ONLY : XNUMOL, XNUMOLAIR
      USE TRACERID_MOD,       ONLY : IDTCO,  IDTCH2O
      USE TROPOPAUSE_MOD,     ONLY : ITS_IN_THE_TROP
      USE TROPOPAUSE_MOD,     ONLY : GET_MIN_TPAUSE_LEVEL

      USE CMN_SIZE_MOD             ! Size parameters

      IMPLICIT NONE
#     include "define.h"

      ! Arguments
      REAL*4,         INTENT(IN) :: COPROD(JJPAR,LLPAR)
      REAL*4,         INTENT(IN) :: COLOSS(JJPAR,LLPAR)
      TYPE(MetState), INTENT(IN) :: State_Met   ! Meteorology State object
#if defined( DEVEL )
      TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
#endif
      ! Local variables
      INTEGER             :: I, J, L, M, N, LMIN

      REAL*8              :: BAIRDENS, DT, GCO, STTTOGCO

      ! Pointers
#if defined( DEVEL )
      ! We need to define local arrays to hold corresponding values 
      ! from the Chemistry State (State_Chm) object. (mpayer, 12/6/12)
      REAL*8, POINTER :: STT(:,:,:,:)
#endif

      ! External functions
      REAL*8, EXTERNAL    :: BOXVL

      !=================================================================
      ! CO_STRAT_PL begins here!
      !=================================================================

#if defined( DEVEL )
      ! Initialize GEOS-Chem tracer array [kg] from Chemistry State object
      ! (mpayer, 12/6/12)
      STT => State_Chm%Tracers
#endif

      ! Chemistry timestep [s]
      DT = GET_TS_CHEM() * 60d0

      !=================================================================
      ! Loop over all stratospheric grid boxes ( L >= LPAUSE(I,J) ). 
      !
      ! Compute the net CO from the P(CO) and L(CO) rates that are 
      ! stored in the COPROD and COLOSS arrays.
      !
      ! Unit conversion to/from [kg/box] and [molec/cm3] is required.
      ! The conversion factor is STTTOGCO, which is given below.
      !
      !   kg CO       box     |   mole CO   | 6.022e23 molec CO       
      !  ------- * -----------+-------------+-------------------  
      !    box      BOXVL cm3 | 28e-3 kg CO |     mole CO             
      !
      !  =  molec CO
      !     --------
      !       cm3
      !=================================================================

      ! Get the minimum extent of the tropopause
      LMIN = GET_MIN_TPAUSE_LEVEL()

      DO L = LMIN, LLPAR
      DO J = 1,    JJPAR
      DO I = 1,    IIPAR

         ! Skip tropospheric grid boxes
         IF ( ITS_IN_THE_TROP( I, J, L, State_Met ) ) CYCLE

         ! conversion factor from [kg/box] to [molec/cm3]
         STTTOGCO = 6.022d23 / ( 28d-3 * BOXVL(I,J,L,State_Met) )

         ! Convert STT from [kg/box] to [molec/cm3]
         GCO = STT(I,J,L,IDTCO) * STTTOGCO

         ! Air density in molec/cm3
         BAIRDENS = State_Met%AD(I,J,L) * XNUMOLAIR /
     &              BOXVL(I,J,L,State_Met)

         ! Apply P(CO) and L(CO) rates to GCO
         GCO = GCO * ( 1d0 - COLOSS(J,L) * DT ) +
     &               ( COPROD(J,L) * DT * BAIRDENS )

         ! Compute production of CH2O (qli, 12/9/99)
         STT(I,J,L,IDTCH2O) = STT(I,J,L,IDTCH2O) + 
     &                        XNUMOL(IDTCO) / XNUMOL(IDTCH2O) * 
     &                        COPROD(J,L)   * BAIRDENS        / 
     &                        STTTOGCO

         ! Convert STT from [molec/cm3] to [kg/box]
         STT(I,J,L,IDTCO) = GCO / STTTOGCO
      ENDDO
      ENDDO
      ENDDO

#if defined( DEVEL )
      ! Free pointer
      NULLIFY( STT )
#endif

      ! Return to calling program
      END SUBROUTINE CO_STRAT_PL
