#if defined( DEVEL )
! $Id: gc_environment_mod.F
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
        PUBLIC :: TRACER_INDEX
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

        INTEGER, ALLOCATABLE :: TRACER_INDEX(:)

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
          USE TRACER_MOD,   ONLY : N_TRACERS

          IMPLICIT NONE

          TYPE(GC_MET_LOCAL), INTENT(OUT) :: LOCAL_MET
          INTEGER :: I,J,L,AS

          ! Allocate 2-D Arrays
          ALLOCATE( &
               LOCAL_MET%ALBD    (IIPAR, JJPAR), & ! Visible surface albedo [unitless]
               LOCAL_MET%AREA_M2 (IIPAR, JJPAR), & ! Grid box surface area [cm2]
               LOCAL_MET%CLDFRC  (IIPAR, JJPAR), & ! Column cloud fraction [unitless]
               LOCAL_MET%CLDTOPS (IIPAR, JJPAR), & ! Max cloud top height [levels]
               LOCAL_MET%EFLUX   (IIPAR, JJPAR), & ! Latent heat flux [W/m2]
               LOCAL_MET%EVAP    (IIPAR, JJPAR), & ! Surface evap [kg/m2/s]
               LOCAL_MET%FRCLND  (IIPAR, JJPAR), & ! Olson land fraction [unitless]
               LOCAL_MET%FRLAKE  (IIPAR, JJPAR), & ! Fraction of lake [unitless]
               LOCAL_MET%FRLAND  (IIPAR, JJPAR), & ! Fraction of land [unitless]
               LOCAL_MET%FRLANDIC(IIPAR, JJPAR), & ! Fraction of land ice [unitless]
               LOCAL_MET%FROCEAN (IIPAR, JJPAR), & ! Fraction of ocean [unitless]
               LOCAL_MET%FRSEAICE(IIPAR, JJPAR), & ! Sfc sea ice fraction
               LOCAL_MET%FRSNO   (IIPAR, JJPAR), & ! Sfc snow fraction
               LOCAL_MET%GRN     (IIPAR, JJPAR), & ! Greenness fraction
               LOCAL_MET%GWETROOT(IIPAR, JJPAR), & ! Root soil wetness [unitless]
               LOCAL_MET%GWETTOP (IIPAR, JJPAR), & ! Top soil moisture [unitless]
               LOCAL_MET%HFLUX   (IIPAR, JJPAR), & ! Sensible heat flux [W/m2]
               LOCAL_MET%LAI     (IIPAR, JJPAR), & ! Leaf area index [m2/m2]
               LOCAL_MET%LWI     (IIPAR, JJPAR), & ! Land/water indices [unitless]
               LOCAL_MET%LWI_GISS(IIPAR, JJPAR), & ! Land fraction [unitless]
               LOCAL_MET%MOLENGTH(IIPAR, JJPAR), & ! Monin-Obhukov length [m]
               LOCAL_MET%OICE    (IIPAR, JJPAR), & ! Fraction of ocean ice [unitless]
               LOCAL_MET%PARDR   (IIPAR, JJPAR), & ! Direct  photsyn active rad [W/m2]
               LOCAL_MET%PARDF   (IIPAR, JJPAR), & ! Diffuse photsyn active rad [W/m2]
               LOCAL_MET%PBLH    (IIPAR, JJPAR), & ! PBL height [m]
               LOCAL_MET%PHIS    (IIPAR, JJPAR), & ! Surface geopotential height [m2/s2]
               LOCAL_MET%PRECANV (IIPAR, JJPAR), & ! Anvil precip @ ground [kg/m2/s]
               LOCAL_MET%PRECCON (IIPAR, JJPAR), & ! Conv  precip @ ground [kg/m2/s]
               LOCAL_MET%PRECTOT (IIPAR, JJPAR), & ! Total precip @ ground [kg/m2/s]
               LOCAL_MET%PRECLSC (IIPAR, JJPAR), & ! Large-scale precip @ ground [kg/m2/s]
               LOCAL_MET%PRECSNO (IIPAR, JJPAR), & ! Snow precip [kg/m2/s]
               LOCAL_MET%PS1     (IIPAR, JJPAR), & ! Surface pressure at start of timestep [hPa]
               LOCAL_MET%PS2     (IIPAR, JJPAR), & ! Surface pressure at end of timestep [hPa]
               LOCAL_MET%PSC2    (IIPAR, JJPAR), & ! Interpolated sfc pressure [hPa]
               LOCAL_MET%RADLWG  (IIPAR, JJPAR), & ! Net LW radiation @ ground [W/m2]
               LOCAL_MET%RADSWG  (IIPAR, JJPAR), & ! Solar radiation @ ground [W/m2]
               LOCAL_MET%SEAICE00(IIPAR, JJPAR), & ! Sea ice coverage 00-10%
               LOCAL_MET%SEAICE10(IIPAR, JJPAR), & ! Sea ice coverage 10-20%
               LOCAL_MET%SEAICE20(IIPAR, JJPAR), & ! Sea ice coverage 20-30%
               LOCAL_MET%SEAICE30(IIPAR, JJPAR), & ! Sea ice coverage 30-40%
               LOCAL_MET%SEAICE40(IIPAR, JJPAR), & ! Sea ice coverage 40-50%
               LOCAL_MET%SEAICE50(IIPAR, JJPAR), & ! Sea ice coverage 50-60%
               LOCAL_MET%SEAICE60(IIPAR, JJPAR), & ! Sea ice coverage 60-70%
               LOCAL_MET%SEAICE70(IIPAR, JJPAR), & ! Sea ice coverage 70-80%
               LOCAL_MET%SEAICE80(IIPAR, JJPAR), & ! Sea ice coverage 80-90%
               LOCAL_MET%SEAICE90(IIPAR, JJPAR), & ! Sea ice coverage 90-100%
               LOCAL_MET%SLP     (IIPAR, JJPAR), & ! Sea level pressure [hPa]
               LOCAL_MET%SNICE   (IIPAR, JJPAR), & ! Fraction of snow/ice [unitless]
               LOCAL_MET%SNODP   (IIPAR, JJPAR), & ! Snow depth [m]
               LOCAL_MET%SNOMAS  (IIPAR, JJPAR), & ! Snow mass [kg/m2]
               LOCAL_MET%SNOW    (IIPAR, JJPAR), & ! Snow depth (H2O equivalent) [mm H2O]
               LOCAL_MET%SST     (IIPAR, JJPAR), & ! Sea surface temperature [K]
               LOCAL_MET%SUNCOS  (IIPAR, JJPAR), & ! Cosine of solar zenith angle
               LOCAL_MET%TO3     (IIPAR, JJPAR), & ! Total overhead O3 column [DU]
               LOCAL_MET%TO31    (IIPAR, JJPAR), & ! Total overhead O3 at start of timestep [DU]
               LOCAL_MET%TO32    (IIPAR, JJPAR), & ! Total overhead O3 at end of timestep [DU]
               LOCAL_MET%TROPP   (IIPAR, JJPAR), & ! Tropopause pressure [hPa]
               LOCAL_MET%TROPP1  (IIPAR, JJPAR), & ! Tropopause P at start of timestep [hPa]
               LOCAL_MET%TROPP2  (IIPAR, JJPAR), & ! Tropopause P at end of timestep [hPa]
               LOCAL_MET%TS      (IIPAR, JJPAR), & ! Surface temperature at 2m [K]
               LOCAL_MET%TSKIN   (IIPAR, JJPAR), & ! Surface skin temperature [K]
               LOCAL_MET%TTO3    (IIPAR, JJPAR), & ! Tropospheric ozone column [DU]
               LOCAL_MET%U10M    (IIPAR, JJPAR), & ! E/W wind speed @ 10m height [m/s]
               LOCAL_MET%USTAR   (IIPAR, JJPAR), & ! Friction velocity [m/s]
               LOCAL_MET%UVALBEDO(IIPAR, JJPAR), & ! UV surface albedo [unitless]
               LOCAL_MET%V10M    (IIPAR, JJPAR), & ! N/S wind speed @ 10m height [m/s]
               LOCAL_MET%Z0      (IIPAR, JJPAR), & ! Surface roughness height [m]
               STAT = AS )
          IF (AS /= 0) CALL ALLOC_ERR('LOCAL_MET 2D')

          ! Allocate 3-D Arrays
          ALLOCATE( &
               ! Fields dimensioned as (I,J,L)
               LOCAL_MET%AD      (IIPAR, JJPAR, LLPAR), & ! Air mass [kg]
               LOCAL_MET%AIRVOL  (IIPAR, JJPAR, LLPAR), & ! Grid box volume [m3]
               LOCAL_MET%AVGW    (IIPAR, JJPAR, LLPAR), & ! Mixing ratio of water vapor
               LOCAL_MET%BXHEIGHT(IIPAR, JJPAR, LLPAR), & ! Grid box height [m]
               LOCAL_MET%CMFMC   (IIPAR, JJPAR, LLPAR), & ! Cloud mass flux [kg/m2/s]
               LOCAL_MET%DETRAINE(IIPAR, JJPAR, LLPAR), & ! GCAP detrainment (entraining plume) [Pa/s]
               LOCAL_MET%DETRAINN(IIPAR, JJPAR, LLPAR), & ! GCAP detrainment (non-entr'n plume) [Pa/s]
               LOCAL_MET%DNDE    (IIPAR, JJPAR, LLPAR), & ! GCAP downdraft   (entraining plume) [Pa/s]
               LOCAL_MET%DNDN    (IIPAR, JJPAR, LLPAR), & ! GCAP downdraft   (non-entr'n plume) [Pa/s]
               LOCAL_MET%DQRCU   (IIPAR, JJPAR, LLPAR), & ! Convective precip production rate [kg/kg/s]
               LOCAL_MET%DQRLSAN (IIPAR, JJPAR, LLPAR), & ! Large-scale precip production rate [kg/kg/s]
               LOCAL_MET%DQIDTMST(IIPAR, JJPAR, LLPAR), & ! Ice tendency, mst proc [kg/kg/s]
               LOCAL_MET%DQLDTMST(IIPAR, JJPAR, LLPAR), & ! H2O tendency, mst proc [kg/kg/s]
               LOCAL_MET%DQVDTMST(IIPAR, JJPAR, LLPAR), & ! Vapor tendency, mst proc [kg/kg/s]
               LOCAL_MET%DTRAIN  (IIPAR, JJPAR, LLPAR), & ! Detrainment flux [kg/m2/s]
               LOCAL_MET%ENTRAIN (IIPAR, JJPAR, LLPAR), & ! GCAP entrainment [Pa/s]
               LOCAL_MET%HKBETA  (IIPAR, JJPAR, LLPAR), & ! Hack overshoot parameter [unitless]
               LOCAL_MET%HKETA   (IIPAR, JJPAR, LLPAR), & ! Hack convective mass flux [kg/m2/s]
               LOCAL_MET%PEDGE   (IIPAR, JJPAR, LLPAR), & ! Pressure @ level edges [Pa]
               LOCAL_MET%PMID    (IIPAR, JJPAR, LLPAR), & ! Pressure @ level centers [Pa]
               LOCAL_MET%PFICU   (IIPAR, JJPAR, LLPAR), & ! Downward flux of ice precip: conv   [kg/m2/s]
               LOCAL_MET%PFILSAN (IIPAR, JJPAR, LLPAR), & ! Downward flux of ice precip: LS+anv [kg/m2/s]
               LOCAL_MET%PFLCU   (IIPAR, JJPAR, LLPAR), & ! Downward flux of liq precip: conv   [kg/m2/s]
               LOCAL_MET%PFLLSAN (IIPAR, JJPAR, LLPAR), & ! Downward flux of ice precip: LS+anv [kg/m2/s]
               LOCAL_MET%PV      (IIPAR, JJPAR, LLPAR), & ! Potential vorticity [kg*m2/kg/s]
               LOCAL_MET%QI      (IIPAR, JJPAR, LLPAR), & ! Cloud ice mixing ratio [kg/kg]
               LOCAL_MET%QL      (IIPAR, JJPAR, LLPAR), & ! Cloud water mixing ratio [kg/kg]
               LOCAL_MET%REEVAPCN(IIPAR, JJPAR, LLPAR), & ! Evap of precip conv condensate [kg/kg/s]
               LOCAL_MET%REEVAPLS(IIPAR, JJPAR, LLPAR), & ! Evap of precip LS+anvil condenstate [kg/kg/s]
               LOCAL_MET%RH      (IIPAR, JJPAR, LLPAR), & ! Relative humidity [unitless]
               LOCAL_MET%RH1     (IIPAR, JJPAR, LLPAR), & ! RH at start of timestep [unitless]
               LOCAL_MET%RH2     (IIPAR, JJPAR, LLPAR), & ! RH at end of timestep [unitless]
               LOCAL_MET%SPHU    (IIPAR, JJPAR, LLPAR), & ! Specific humidity [kg/kg]
               LOCAL_MET%SPHU1   (IIPAR, JJPAR, LLPAR), & ! Specific humidity at start of timestep[kg/kg]
               LOCAL_MET%SPHU2   (IIPAR, JJPAR, LLPAR), & ! Specific humidity at end of timestep [kg/kg]
               LOCAL_MET%T       (IIPAR, JJPAR, LLPAR), & ! Temperature [K]
               LOCAL_MET%TAUCLI  (IIPAR, JJPAR, LLPAR), & ! Optical depth of ice clouds [unitless]
               LOCAL_MET%TAUCLW  (IIPAR, JJPAR, LLPAR), & ! Optical depth of H2O clouds [unitless]
               LOCAL_MET%TMPU1   (IIPAR, JJPAR, LLPAR), & ! Temperature at start of timestep [K]
               LOCAL_MET%TMPU2   (IIPAR, JJPAR, LLPAR), & ! Temperature at end of timestep [K]
               LOCAL_MET%U       (IIPAR, JJPAR, LLPAR), & ! E/W component of wind [m s-1]
               LOCAL_MET%UPDE    (IIPAR, JJPAR, LLPAR), & ! GCAP updraft (entraining plume) [Pa/s]
               LOCAL_MET%UPDN    (IIPAR, JJPAR, LLPAR), & ! GCAP updraft (non-entr'n plume) [Pa/s]
               LOCAL_MET%V       (IIPAR, JJPAR, LLPAR), & ! N/S component of wind [m s-1]
               LOCAL_MET%ZMEU    (IIPAR, JJPAR, LLPAR), & ! Zhang/McFarlane updraft entrainment [Pa/s]
               LOCAL_MET%ZMMD    (IIPAR, JJPAR, LLPAR), & ! Zhang/McFarlane downdraft mass flux [Pa/s]
               LOCAL_MET%ZMMU    (IIPAR, JJPAR, LLPAR), & ! Zhang/McFarlane updraft   mass flux [Pa/s]
               ! Fields dimensioned as (L,I,J)
               LOCAL_MET%AIRDENS (LLPAR, IIPAR, JJPAR), & ! Air density [kg/m3]
               LOCAL_MET%CLDF    (LLPAR, IIPAR, JJPAR), & ! 3-D cloud fraction [unitless]
               LOCAL_MET%DELP    (LLPAR, IIPAR, JJPAR), & ! Pressure thickness for layer [Pa]
               LOCAL_MET%MOISTQ  (LLPAR, IIPAR, JJPAR), & ! Tendency in sp. humidity [kg/kg/s]
               LOCAL_MET%OPTD    (LLPAR, IIPAR, JJPAR), & ! Visible optical depth [unitless]
               STAT=AS)
          
          IF (AS /= 0) CALL ALLOC_ERR('LOCAL_MET 3D')

          DO I = 1,IIPAR
             DO J = 1, JJPAR
                DO L = 1, LLPAR
!                   LOCAL_MET%PEDGE(I,J,L) = GET_PEDGE(I,J,L)
!                   LOCAL_MET%PMID (I,J,L) = GET_PCENTER(I,J,L)
                ENDDO
             ENDDO
          ENDDO

!          ALLOCATE( TRACER_INDEX(N_TRACERS), STAT = AS )

        END SUBROUTINE INIT_LOCAL_MET
           
      END MODULE GC_ENVIRONMENT_MOD
!EOC
#endif
