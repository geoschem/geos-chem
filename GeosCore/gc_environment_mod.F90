#if defined( DEVEL )
! $Id: gc_interface_mod.F,v 1.1 2009/12/11 20:57:28 dasilva Exp $
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: gc_environment_mod
!
! !DESCRIPTION: Module GC\_ENVIRONMENT\_MOD establishes the runtime environment
! for the GEOS-Chem model. It is designed to receive model parameter and geophys.
! environment information and allocate memory based upon it.
!
! It provides routines to do the following:
! (1) Allocate geo-spatial arrays
! (2) Initialize met. field derived type.
! (3) Initialize CHEM, PHYS, and EMISSIONS states
! (4) ...
!\\
!\\
!  NOTE: This is mostly for testing the grid-independent code in the current 
!  GEOS-Chem.  Many of these inputs will come from the GEOS-5 interface. It will
!  remain in DEVEL state for some time.
!\\
!\\
! !INTERFACE: 
!
      MODULE GC_ENVIRONMENT_MOD
!
! !USES
!        
        USE GC_TYPE_MOD                  ! Various derived type definitions

        IMPLICIT NONE
#include "define.h"

        PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
        PUBLIC :: ALLOCATE_ALL
        PUBLIC :: INIT_ALL
!
! !PRIVATE MEMBER FUNCTIONS:
!
        
!
! !REVISION HISTORY:
!  26 Jan 2012 - M. Long - Created module file
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !PRIVATE TYPES:
!
      !=================================================================
      ! MODULE ROUTINES -- follow below the "CONTAINS" statement
      !=================================================================
      CONTAINS
        
!------------------------------------------------------------------------------
        SUBROUTINE ALLOCATE_ALL
!
!******************************************************************************
! Subroutine ALLOCATE_ALL allocates all LAT/LON ALLOCATABLE arrays for global
! use by the GEOS-Chem either as a standalone program or module.
!
! NOTES:
!******************************************************************************

          USE GC_TYPE2_MOD,      ONLY : CHEM_STATE, INIT_CHEMSTATE
          USE CMN_DEP_MOD,       ONLY : SET_CMN_DEP_MOD
          USE CMN_ISOP_MOD,      ONLY : SET_CMN_ISOP_MOD
          USE CMN_MONOT_MOD,     ONLY : SET_CMN_MONOT_MOD
          USE CMN_NOX_MOD,       ONLY : SET_CMN_NOX_MOD
          USE CMN_O3_MOD,        ONLY : SET_CMN_O3_MOD
          USE CMN_VEL_MOD,       ONLY : SET_CMN_VEL_MOD
          USE CMN_MOD,           ONLY : SET_CMN_MOD
          USE CMN_FJ_MOD,        ONLY : SET_CMN_FJ_MOD
          USE JV_CMN_MOD,        ONLY : SET_JV_CMN_MOD
          USE COMMSOIL_MOD,      ONLY : SET_COMMSOIL_MOD
          USE VDIFF_PRE_MOD,     ONLY : SET_VDIFF_PRE_MOD
          USE CMN_SIZE_MOD,      ONLY : SET_CMN_SIZE_MOD
          USE CMN_DIAG_MOD,      ONLY : SET_CMN_DIAG_MOD
          USE COMODE_LOOP_MOD,   ONLY : SET_COMODE_LOOP_MOD
          
          IMPLICIT NONE
          
          CALL SET_CMN_SIZE_MOD
          CALL SET_CMN_DEP_MOD
          CALL SET_CMN_DIAG_MOD
          CALL SET_CMN_ISOP_MOD
          CALL SET_CMN_MONOT_MOD
          CALL SET_CMN_NOX_MOD
          CALL SET_CMN_O3_MOD
          CALL SET_CMN_VEL_MOD
          CALL SET_CMN_MOD
          CALL SET_CMN_FJ_MOD
          CALL SET_COMMSOIL_MOD
          CALL SET_COMODE_LOOP_MOD
          CALL SET_JV_CMN_MOD
          
          CALL SET_VDIFF_PRE_MOD
          
        END SUBROUTINE ALLOCATE_ALL
!------------------------------------------------------------------------------

        SUBROUTINE INIT_ALL(LOCAL_MET, CHEM_STATE) 
!
!******************************************************************************
! Subroutine INIT_ALL initializes the top-level data structures that are either
! passed to/from GC or between GC components (emis->transport->chem->etc)
!
! NOTES:
!******************************************************************************
          USE GC_TYPE_MOD
          USE GC_TYPE2_MOD, ONLY : INIT_CHEMSTATE, CHEMSTATE

          IMPLICIT NONE

          TYPE(GC_MET_LOCAL), INTENT(OUT) :: LOCAL_MET
          TYPE(CHEMSTATE   ), INTENT(OUT) :: CHEM_STATE

          CALL INIT_LOCAL_MET(LOCAL_MET)    ! Initializes the Met derived type bundle
          CALL INIT_CHEMSTATE(CHEM_STATE)   ! Initializes the Met derived type bundle

        END SUBROUTINE INIT_ALL

        SUBROUTINE INIT_LOCAL_MET(LOCAL_MET)
          USE GC_TYPE_MOD
          USE CMN_SIZE_MOD, ONLY : IIPAR, JJPAR, LLPAR
          USE ERROR_MOD,    ONLY : ALLOC_ERR
          USE PRESSURE_MOD, ONLY : GET_PEDGE, GET_PCENTER

          IMPLICIT NONE

          TYPE(GC_MET_LOCAL), INTENT(OUT) :: LOCAL_MET
          INTEGER :: I,J,L,AS

          ALLOCATE( &
               LOCAL_MET%AD(IIPAR,JJPAR,LLPAR), &        ! Air mass [kg]
               LOCAL_MET%AIRDENS(LLPAR,IIPAR,JJPAR), &   ! Air density [kg/m3]
               LOCAL_MET%AIRVOL(IIPAR,JJPAR,LLPAR), &    ! Grid box volume [m3]
               LOCAL_MET%BXHEIGHT(IIPAR,JJPAR,LLPAR), &  ! Grid box height [m]
               LOCAL_MET%CLDF(LLPAR,IIPAR,JJPAR), &      ! 3-D cloud fraction [unitless]
               LOCAL_MET%CMFMC(IIPAR,JJPAR,LLPAR), &     ! Cloud mass flux [kg/m2/s]
               LOCAL_MET%DQIDTMST(IIPAR,JJPAR,LLPAR), &  ! Ice tendency, mst proc [kg/kg/s]
               LOCAL_MET%DQLDTMST(IIPAR,JJPAR,LLPAR), &  ! H2O tendency, mst proc [kg/kg/s]
               LOCAL_MET%DQVDTMST(IIPAR,JJPAR,LLPAR), &  ! Vapor tendency, mst proc [kg/kg/s]
               LOCAL_MET%DTRAIN(IIPAR,JJPAR,LLPAR), &    ! Detrainment flux [kg/m2/s]
               LOCAL_MET%MOISTQ(LLPAR,IIPAR,JJPAR), &    ! Tendency in sp. humidity [kg/kg/s]
               LOCAL_MET%OPTD(LLPAR,IIPAR,JJPAR), &      ! Visible optical depth [unitless]
               LOCAL_MET%PEDGE(IIPAR,JJPAR,LLPAR), &     ! Pressure @ level edges [Pa]
               LOCAL_MET%PMID(IIPAR,JJPAR,LLPAR), &      ! Pressure @ level centers [Pa]
               LOCAL_MET%DELP(LLPAR,IIPAR,JJPAR), &      !
               LOCAL_MET%RH(IIPAR,JJPAR,LLPAR), &        ! Relative humidity [unitless]
               LOCAL_MET%SPHU(IIPAR,JJPAR,LLPAR), &      ! Specific humidity [kg/kg]
               LOCAL_MET%T(IIPAR,JJPAR,LLPAR), &         ! Temperature [K]
               LOCAL_MET%TAUCLI(IIPAR,JJPAR,LLPAR), &    ! Opt depth of ice clouds [unitless]
               LOCAL_MET%TAUCLW(IIPAR,JJPAR,LLPAR), &    ! Opt depth of H2O clouds [unitless]
               STAT=AS)
          
          IF (AS /= 0) CALL ALLOC_ERR('LOCAL_MET')

          DO I = 1,IIPAR
             DO J = 1, JJPAR
                DO L = 1, LLPAR
                   LOCAL_MET%PEDGE(I,J,L) = GET_PEDGE(I,J,L)
                   LOCAL_MET%PMID (I,J,L) = GET_PCENTER(I,J,L)
                ENDDO
             ENDDO
          ENDDO

        END SUBROUTINE INIT_LOCAL_MET
           
      END MODULE GC_ENVIRONMENT_MOD
!EOC
#endif
!-----------------------------------------------------------------------
! TO BE DONE:
!
!  ARRAY     !  COMPLETE?  ! NOTE:
!------------!-------------!--------------------------------------------
!  AD        !   YES       ! SET IN DAO_MOD.AIRQNT IN STD GC
!  AIRDENS   !   YES       ! SET IN DAO_MOD.AIRQNT IN STD GC
!  AIRVOL    !   YES       ! SET IN DAO_MOD.AIRQNT IN STD GC
!  BXHEIGHT  !   YES       ! SET IN DAO_MOD.AIRQNT IN STD GC
!  CLDF      !   YES       ! SET IN A6_READMOD
!  CMFMC     !   YES       ! SET IN A6_READMOD
!  DQIDTMST  !   YES       ! SET IN A6_READMOD
!  DQLDTMST  !   YES       ! SET IN A6_READMOD
!  DQVDTMST  !   YES       ! SET IN A6_READMOD
!  DTRAIN    !   YES       ! SET IN A6_READMOD
!  MOISTQ    !   YES       ! SET IN A6_READMOD
!  OPTD      !   YES       ! SET IN A6_READMOD
!  PEDGE     !   YES       ! SET LOCALLY (ABOVE)
!  PMID      !   YES       ! SET LOCALLY (ABOVE)
!  DELP      !   YES       ! SET IN DAO_MOD.AIRQNT IN STD GC
!  RH        !   YES       ! SET IN A6_READMOD
!  SPHU      !   YES       ! SET IN A6_READMOD
!  T         !   YES       ! SET IN A6_READMOD
!  TAUCLI    !   YES       ! SET IN A6_READMOD
!  TAUCLW    !   YES       ! SET IN A6_READMOD
!-----------------------------------------------------------------------
