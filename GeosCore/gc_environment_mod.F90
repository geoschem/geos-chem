#if defined( DEVEL )
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: gc_environment_mod
!
! !DESCRIPTION: Module GC\_ENVIRONMENT\_MOD establishes the runtime 
!  environment for the GEOS-Chem model. It is designed to receive model 
!  parameter and geophysical environment information and allocate memory 
!  based upon it.
!
!  It provides routines to do the following:
!  (1) Allocate geo-spatial arrays
!  (2) Initialize met. field derived type.
!  (3) Initialize CHEM, PHYS, and EMISSIONS states
!  (4) ...
!\\
!\\
!  NOTE: This is mostly for testing the grid-independent code in the current 
!  GEOS-Chem.  Many of these inputs will come from the GEOS-5 interface. It will
!  remain in DEVEL state for some time.
!\\
!\\
! !INTERFACE: 
!
MODULE GC_Environment_Mod
!
! !USES
!        
  USE GC_TYPE_MOD                  ! Various derived type definitions

  IMPLICIT NONE
# include "define.h"

  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC :: ALLOCATE_ALL
  PUBLIC :: INIT_ALL
!
! !PRIVATE MEMBER FUNCTIONS:
!
  PRIVATE :: INIT_LOCAL_MET
!
! !REMARKS:
!  For consistency, we should probably move the met state initialization
!  to the same module where the met state derived type is contained.
!
! !REVISION HISTORY:
!  26 Jan 2012 - M. Long     - Created module file
!  13 Aug 2012 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
CONTAINS
!EOC        
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: allocate_all
!
! !DESCRIPTION: Subroutine ALLOCATE\_ALL allocates all LAT/LON ALLOCATABLE 
!  arrays for global use by the GEOS-Chem either as a standalone program or
!  module.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE ALLOCATE_ALL
!
! !USES:
!
    USE CMN_DEP_MOD,       ONLY : SET_CMN_DEP_MOD
    USE CMN_NOX_MOD,       ONLY : SET_CMN_NOX_MOD
    USE CMN_O3_MOD,        ONLY : SET_CMN_O3_MOD
    USE CMN_MOD,           ONLY : SET_CMN_MOD
    USE CMN_FJ_MOD,        ONLY : SET_CMN_FJ_MOD
    USE JV_CMN_MOD,        ONLY : SET_JV_CMN_MOD
    USE COMMSOIL_MOD,      ONLY : SET_COMMSOIL_MOD
    USE VDIFF_PRE_MOD,     ONLY : SET_VDIFF_PRE_MOD
    USE CMN_SIZE_MOD,      ONLY : SET_CMN_SIZE_MOD
    USE CMN_DIAG_MOD,      ONLY : SET_CMN_DIAG_MOD
    USE COMODE_LOOP_MOD,   ONLY : SET_COMODE_LOOP_MOD
          
    IMPLICIT NONE
! 
! !REVISION HISTORY: 
!  26 Jan 2012 - M. Long     - Initial version
!  13 Aug 2012 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
    CALL SET_CMN_SIZE_MOD
    CALL SET_CMN_DEP_MOD
    CALL SET_CMN_DIAG_MOD
    CALL SET_CMN_NOX_MOD
    CALL SET_CMN_O3_MOD
    CALL SET_CMN_MOD
    CALL SET_CMN_FJ_MOD
    CALL SET_COMMSOIL_MOD
    CALL SET_COMODE_LOOP_MOD
    CALL SET_JV_CMN_MOD
    
    CALL SET_VDIFF_PRE_MOD
          
  END SUBROUTINE ALLOCATE_ALL
!EOC
 !------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_all
!
! !DESCRIPTION: Subroutine INIT\_ALL initializes the top-level data structures
!  that are either passed to/from GC or between GC components 
!  (emis->transport->chem->etc)
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Init_All( State_Met, State_Chm, am_I_Root, RC ) 

!
! !USES:
!
    USE GC_TYPE_MOD
    USE GC_TYPE2_MOD, ONLY : ChemState
    USE GC_TYPE2_MOD, ONLY : Init_Chemistry_State

    IMPLICIT NONE
!
! !INPUT PARAMETERS:
!
    LOGICAL,            INTENT(IN)    :: am_I_Root   ! Are we on the root CPU?
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(GC_MET_LOCAL), INTENT(INOUT) :: State_Met   ! Meteorology state
    TYPE(CHEMSTATE   ), INTENT(INOUT) :: State_Chm   ! Chemistry state
!
! !OUTPUT PARAMETERS:
!
    INTEGER,            INTENT(OUT)   :: RC          ! Success or failure
!
! !REVISION HISTORY: 
!  26 Jan 2012 - M. Long     - Initial version
!  13 Aug 2012 - R. Yantosca - Added ProTeX headers
!  16 Oct 2012 - R. Yantosca - Renamed LOCAL_MET argument to State_Met
!  16 Oct 2012 - R. Yantosca - Renamed GC_STATE  argument to State_Chm
!  16 Oct 2012 - R. Yantosca - Call Init_Chemistry_State (in gc_type2_mod.F90,
!                              which was renamed from INIT_CHEMSTATE)
!EOP
!------------------------------------------------------------------------------
!BOC

    ! Initialize object for met fields
    CALL INIT_LOCAL_MET( State_Met )

    ! Initialize object for chemical state
    CALL Init_Chemistry_State( State_Chm, am_I_Root, RC )
    
  END SUBROUTINE Init_All
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_local_met
!
! !DESCRIPTION: Subroutine INIT\_LOCAL\_MET allocates all fields of an
!  object based on derived type GC_MET_LOCAL (from gc_type_mod.F90).
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE INIT_LOCAL_MET( State_Met )
!
! !USES:
!
    USE GC_TYPE_MOD
    USE CMN_SIZE_MOD, ONLY : IIPAR, JJPAR, LLPAR
    USE ERROR_MOD,    ONLY : ALLOC_ERR
    USE PRESSURE_MOD, ONLY : GET_PEDGE, GET_PCENTER
    USE TRACER_MOD,   ONLY : N_TRACERS
    
    IMPLICIT NONE
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(GC_MET_LOCAL), INTENT(INOUT) :: State_Met
!
! !REMARKS:
!  For consistency, maybe this should be moved to a different module.
!
! !REVISION HISTORY: 
!  01 Oct 1995 - R. Yantosca - Initial version
!  08 Dec 2009 - R. Yantosca - Added ProTeX headers
!  16 Oct 2012 - R. Yantosca - Renamed LOCAL_MET argument to State_Met
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER :: I,J,L,AS

    !=======================================================================
    ! Allocate 2-D Arrays
    !=======================================================================
    ALLOCATE( &
       State_Met%ALBD    (IIPAR,JJPAR), & ! Visible surface albedo [unitless]
       State_Met%CLDFRC  (IIPAR,JJPAR), & ! Column cloud fraction [unitless]
       State_Met%FRCLND  (IIPAR,JJPAR), & ! Olson land fraction [unitless]
       State_Met%GWETTOP (IIPAR,JJPAR), & ! Top soil moisture [unitless]
       State_Met%HFLUX   (IIPAR,JJPAR), & ! Sensible heat flux [W/m2]
       State_Met%LWI     (IIPAR,JJPAR), & ! Land/water indices [unitless]
       State_Met%PARDR   (IIPAR,JJPAR), & ! Direct  photsyn active rad [W/m2]
       State_Met%PARDF   (IIPAR,JJPAR), & ! Diffuse photsyn active rad [W/m2]
       State_Met%PBLH    (IIPAR,JJPAR), & ! PBL height [m]
       State_Met%PRECCON (IIPAR,JJPAR), & ! Conv  precip @ ground [kg/m2/s]
       State_Met%PRECTOT (IIPAR,JJPAR), & ! Total precip @ ground [kg/m2/s]
       State_Met%RADSWG  (IIPAR,JJPAR), & ! Solar radiation @ ground [W/m2]
       State_Met%SST     (IIPAR,JJPAR), & ! Sea surface temperature [K]
       State_Met%SUNCOS  (IIPAR*JJPAR), & ! Cosine of solar zenith angle
       State_Met%TO3     (IIPAR,JJPAR), & ! Total overhead O3 column [DU]
       State_Met%TROPP   (IIPAR,JJPAR), & ! Tropopause pressure [hPa]
       State_Met%TS      (IIPAR,JJPAR), & ! Surface temperature [K]
       State_Met%U10M    (IIPAR,JJPAR), & ! E/W wind speed @ 10m height [m/s]
       State_Met%USTAR   (IIPAR,JJPAR), & ! Friction velocity [m/s]
       State_Met%UVALBEDO(IIPAR,JJPAR), & ! UV surface albedo [unitless]
       State_Met%V10M    (IIPAR,JJPAR), & ! N/S wind speed @ 10m height [m/s]
       State_Met%Z0      (IIPAR,JJPAR), &
       STAT = AS )
    IF (AS /= 0) CALL ALLOC_ERR('INIT_LOCAL_MET 2D')

    !=======================================================================
    ! Allocate 3-D Arrays
    !=======================================================================
    ALLOCATE( &
       State_Met%AD      (IIPAR,JJPAR,LLPAR  ), & ! Air mass [kg]
       State_Met%AIRDENS (LLPAR,IIPAR,JJPAR  ), & ! Air density [kg/m3]
       State_Met%AIRVOL  (IIPAR,JJPAR,LLPAR  ), & ! Grid box volume [m3]
       State_Met%AREA_M2 (IIPAR,JJPAR,LLPAR  ), & ! Grid box surface area [cm2]
       State_Met%BXHEIGHT(IIPAR,JJPAR,LLPAR  ), & ! Grid box height [m]
       State_Met%CLDF    (LLPAR,IIPAR,JJPAR  ), & ! 3-D cloud fraction [1]
       State_Met%CMFMC   (IIPAR,JJPAR,LLPAR+1), & ! Cloud mass flux [kg/m2/s]
       State_Met%DQIDTMST(IIPAR,JJPAR,LLPAR  ), & ! Ice tndcy, mst [kg/kg/s]
       State_Met%DQLDTMST(IIPAR,JJPAR,LLPAR  ), & ! H2O tndcy, mst [kg/kg/s]
       State_Met%DQVDTMST(IIPAR,JJPAR,LLPAR  ), & ! Vapor tndcy, mst [kg/kg/s]
       State_Met%DTRAIN  (IIPAR,JJPAR,LLPAR  ), & ! Detrainment flux [kg/m2/s]
       State_Met%MOISTQ  (LLPAR,IIPAR,JJPAR  ), & ! Tndcy in spc hum [kg/kg/s]
       State_Met%OPTD    (LLPAR,IIPAR,JJPAR  ), & ! Visible opt depth [1]
       State_Met%PEDGE   (IIPAR,JJPAR,LLPAR+1), & ! Pressure @ edges [Pa]
       State_Met%PMID    (IIPAR,JJPAR,LLPAR  ), & ! Pressure @ centers [Pa]
       State_Met%DELP    (LLPAR,IIPAR,JJPAR  ), & ! Pressure thickness [Pa]
       State_Met%RH      (IIPAR,JJPAR,LLPAR  ), & ! Relative humidity [1]
       State_Met%SPHU    (IIPAR,JJPAR,LLPAR  ), & ! Specific humidity [kg/kg]
       State_Met%T       (IIPAR,JJPAR,LLPAR  ), & ! Temperature [K]
       State_Met%TAUCLI  (IIPAR,JJPAR,LLPAR  ), & ! Opt depth of ice clouds [1]
       State_Met%TAUCLW  (IIPAR,JJPAR,LLPAR  ), & ! Opt depth of H2O clouds [1]
       STAT=AS)
    
    IF (AS /= 0) CALL ALLOC_ERR('INIT_LOCAL_MET 3D')
    
! Comment out DO loop
!    DO I = 1,IIPAR
!       DO J = 1, JJPAR
!          DO L = 1, LLPAR
!!             State_Met%PEDGE(I,J,L) = GET_PEDGE(I,J,L)
!!                   State_Met%PMID (I,J,L) = GET_PCENTER(I,J,L)
!                ENDDO
!             ENDDO
!          ENDDO

!          ALLOCATE( TRACER_INDEX(N_TRACERS), STAT = AS )
!          write(*,*) 'TRACER_INDEX: ', shape(tracer_index)

  END SUBROUTINE INIT_LOCAL_MET
!EOC           
END MODULE GC_Environment_Mod
#endif
