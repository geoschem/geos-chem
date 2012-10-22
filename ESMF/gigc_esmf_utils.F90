#if defined( ESMF_)
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: gigc_esmf_utils
!
! !DESCRIPTION: Module GIGC\_ESMF\_UTILS is the module that ...
!\\
!\\
! !INTERFACE: 
!      
MODULE GIGC_ESMF_Utils
!
! !USES:
!     
  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC  :: ESMF2GeosChem
!
! !REMARKS:
!  This module may just be for the BCC interface. 
!
! !REVISION HISTORY:
!  22 Jun 2009 - R. Yantosca & P. Le Sager - Chunkized & cleaned up.
!  16 Oct 2012 - R. Yantosca - Renamed GC_MET argument to State_Met
!  16 Oct 2012 - R. Yantosca - Renamed GC_STATE argument to State_Chm
!  22 Oct 2012 - R. Yantosca - Renamed to gigc_esmf_utils.F90

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
! !IROUTINE: esmf2geoschem 
!
! !DESCRIPTION: Routine ESMF2GeosChem initializes ...
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE ESMF2GeosChem( State_Met, State_Chm )
!
! !USES:
!
    USE GIGC_State_Chm_Mod, ONLY : ChmState
    USE GIGC_State_Met_Mod, ONLY : MetState
!
! !INPUT/OUTPUT PARAMETERS: 
!
    TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
    TYPE(MetState), INTENT(INOUT) :: State_Met   ! Meteorology State object
!
! !REMARKS:
!  This routine may just be for the BCC interface. 
! 
! !REVISION HISTORY: 
!  16 Oct 2012 - M. Long     - Initial version
!  16 Oct 2012 - R. Yantosca - Added ProTeX headers
!  19 Oct 2012 - R. Yantosca - Now reference gigc_state_chm_mod.F90
!  19 Oct 2012 - R. Yantosca - Now reference gigc_state_met_mod.F90
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    !----------------------------------------
    ! LINK DATA FIELDS FROM BCC TO GEOS-CHEM
    ! THESE ARE CALCULATED ONLINE
    !----------------------------------------
    
    ! TRACER MASS MIXING-RATIO ARRAY
    State_Chm%TRACERS = 1.
    
    ! ALBEDO
    State_Met%ALBD = 1.
    
    ! COLUMN-INTEGRATED CLOUD FRACTION
    State_Met%CLDFRC = 1.
    
    ! SPECIFIC HUMIDITY
    State_Met%SPHU = 1.
    
    ! LAND FRACTION
    State_Met%FRCLND = 1.
    
    ! LAND-WATER-ICE FLAGS
    
    ! PBL HEIGHT (M)
    State_Met%PBLH = 1.

    ! CONVECTIVE PRECIP
    State_Met%PRECCON = 1.

    ! LS PRECIP  
    State_Met%PRECTOT = 1.

    ! AN PRECIP
    !!!!State_Met%PTROP(:) ! IS THIS CORRECT?

    ! TROPOPAUSE PRESSURE
    State_Met%TROPP = 1.

    ! SURFACE TEMPERATURE (2-METER AIR TEMPERATURE??)
    State_Met%TS = 1.

    ! RADIATION @ GROUND (NET, LONG OR SHORTWAVE?)
    State_Met%RADSWG = 1.

    ! SEA-SURFACE TEMPERATURE
    State_Met%SST = 1.

    ! 10-METER MERIDIONAL WIND-SPEED
    State_Met%U10M = 1.

    ! 10-METER ZONAL WIND-SPEED
    State_Met%V10M = 1.

    ! FRICTION VELOCITY (U*)
    State_Met%USTAR = 1.

    ! ROUGHNESS HEIGHT
    State_Met%Z0 = 1.

    ! CLOUD FRACTION
    State_Met%CLDF = 1.

    ! CONVECTIVE (CLOUD) MASS FLUX
    State_Met%CMFMC = 1.

    ! DETRAINMENT FLUX (CONVECTIVE: SHALLOW & DEEP?)
    State_Met%DTRAIN = 1.

    ! VISIBLE (ICE) OPTICAL DEPTH
    State_Met%TAUCLI = 1.

    ! VISIBLE (LIQUID WATER) OPTICAL DEPTH
    State_Met%TAUCLW = 1.

    ! RELATIVE HUMIDITY
    State_Met%RH = 1.

    ! SENSIBLE HEAT FLUX
    State_Met%HFLUX = 1.

    ! AIR TEMPERATURE
    State_Met%T = 1.

    ! AIR DENSITY
    State_Met%AIRDENS = 1.

    ! COSINE OF SOLAR ZENITH-ANGLE
    State_Met%SUNCOS = 1.

    ! MOIST (Q) TENDENCY
    State_Met%MOISTQ = 1.

    ! MODEL LAYER PRESSURE-THICKNESS (PDEL)
    State_Met%DELP = 1.

    ! AIR PRESSURE (GRIDBOX EDGES)
    State_Met%PEDGE = 1.

    ! MID-LEVEL PERSSURES
    State_Met%PMID = 1.

    ! TOTAL (COLUMN INTEGRATED) OVERHEAD OZONE COLUMN
    !State_Met%

    ! GRID-BOX AREA
    !State_Met%AREA_M2(:,:)

    ! TOTAL OPTICAL DEPTH
    State_Met%OPTD = 1.

    ! DIRECT PAR (THIS MAY NOT BE THE CORRECT VALUE)
    State_Met%PARDR = 1.
    
    ! DIFFUSE PAR
    State_Met%PARDF = 1.

    ! UV ALBEDO
    State_Met%UVALBEDO = 1.

    State_Met%BXHEIGHT = 1.

    State_Met%DQIDTMST = 1.
    State_Met%DQLDTMST = 1.
    
    State_Met%DQVDTMST = 1.
    
    State_Met%TS = 1.
    
    !State_Met%LAT = 1.

    !State_Met%LON = 1.

    State_Met%GWETTOP = 1.

    !State_Met%ZMID = 1.

  END SUBROUTINE ESMF2GeosChem
!EOC
END MODULE GIGC_Esmf_Utils
#endif
