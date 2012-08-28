#if defined( ESMF_)
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: gc_chunk_mod
!
! !DESCRIPTION: Module GC\_CHUNK\_MOD is the module that contains 
!  the GEOS-Chem chunk code init, run and finalize methods.
!\\
!\\
! !INTERFACE: 
!      
  MODULE GC_ESMF_UTILS
!
! !USES:
!      
    USE GC_TYPE_MOD          
    
    IMPLICIT NONE
    PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
    PUBLIC  :: ESMF2GEOSCHEM
!
! !PRIVATE MEMBER FUNCTIONS:
!

!
! !REVISION HISTORY:
!  22 Jun 2009 - R. Yantosca & P. Le Sager - Chunkized & cleaned up.
!
!EOP
!------------------------------------------------------------------------------
!BOC
  CONTAINS

    SUBROUTINE ESMF2GEOSCHEM(GC_MET, GC_STATE)

      USE GC_TYPE2_MOD
      TYPE(CHEMSTATE),    INTENT(OUT) :: GC_STATE
      TYPE(GC_MET_LOCAL), INTENT(OUT) :: GC_MET

      !----------------------------------------
      ! LINK DATA FIELDS FROM BCC TO GEOS-CHEM
      ! THESE ARE CALCULATED ONLINE
      !----------------------------------------
      ! TRACER MASS MIXING-RATIO ARRAY
      GC_STATE%TRACERS = 1.
      ! ALBEDO
      GC_MET%ALBD = 1.
      ! COLUMN-INTEGRATED CLOUD FRACTION
      GC_MET%CLDFRC = 1.
      ! SPECIFIC HUMIDITY
      GC_MET%SPHU = 1.
      ! LAND FRACTION
      GC_MET%FRCLND = 1.
      ! LAND-WATER-ICE FLAGS
      
      ! PBL HEIGHT (M)
      GC_MET%PBLH = 1.
      ! CONVECTIVE PRECIP
      GC_MET%PRECCON = 1.
      ! LS PRECIP
      GC_MET%PRECTOT = 1.
      ! AN PRECIP
!!!!      GC_MET%PTROP(:) ! IS THIS CORRECT?
      ! TROPOPAUSE PRESSURE
      GC_MET%TROPP = 1.
      ! SURFACE TEMPERATURE (2-METER AIR TEMPERATURE??)
      GC_MET%TS = 1.
      ! RADIATION @ GROUND (NET, LONG OR SHORTWAVE?)
      GC_MET%RADSWG = 1.
      ! SEA-SURFACE TEMPERATURE
      GC_MET%SST = 1.
      ! 10-METER MERIDIONAL WIND-SPEED
      GC_MET%U10M = 1.
      ! 10-METER ZONAL WIND-SPEED
      GC_MET%V10M = 1.
      ! FRICTION VELOCITY (U*)
      GC_MET%USTAR = 1.
      ! ROUGHNESS HEIGHT
      GC_MET%Z0 = 1.
      ! CLOUD FRACTION
      GC_MET%CLDF = 1.
      ! CONVECTIVE (CLOUD) MASS FLUX
      GC_MET%CMFMC = 1.
      ! DETRAINMENT FLUX (CONVECTIVE: SHALLOW & DEEP?)
      GC_MET%DTRAIN = 1.
      ! VISIBLE (ICE) OPTICAL DEPTH
      GC_MET%TAUCLI = 1.
      ! VISIBLE (LIQUID WATER) OPTICAL DEPTH
      GC_MET%TAUCLW = 1.
      ! RELATIVE HUMIDITY
      GC_MET%RH = 1.
      ! SENSIBLE HEAT FLUX
      GC_MET%HFLUX = 1.
      ! AIR TEMPERATURE
      GC_MET%T = 1.
      ! AIR DENSITY
      GC_MET%AIRDENS = 1.
      ! COSINE OF SOLAR ZENITH-ANGLE
      GC_MET%SUNCOS = 1.
      ! MOIST (Q) TENDENCY
      GC_MET%MOISTQ = 1.
      ! MODEL LAYER PRESSURE-THICKNESS (PDEL)
      GC_MET%DELP = 1.
      ! AIR PRESSURE (GRIDBOX EDGES)
      GC_MET%PEDGE = 1.
      ! MID-LEVEL PERSSURES
      GC_MET%PMID = 1.
      ! TOTAL (COLUMN INTEGRATED) OVERHEAD OZONE COLUMN
!      GC_MET%
      ! GRID-BOX AREA
!      GC_MET%AREA_M2(:,:)
      ! TOTAL OPTICAL DEPTH
      GC_MET%OPTD = 1.
      ! DIRECT PAR (THIS MAY NOT BE THE CORRECT VALUE)
      GC_MET%PARDR = 1.
      ! DIFFUSE PAR
      GC_MET%PARDF = 1.
      ! UV ALBEDO
      GC_MET%UVALBEDO = 1.
      GC_MET%BXHEIGHT = 1.
      GC_MET%DQIDTMST = 1.
      GC_MET%DQLDTMST = 1.
      GC_MET%DQVDTMST = 1.
      GC_MET%TS = 1.
!      GC_MET%LAT = 1.
!      GC_MET%LON = 1.
      GC_MET%GWETTOP = 1.
!      GC_MET%ZMID = 1.


    END SUBROUTINE ESMF2GEOSCHEM

  END MODULE GC_ESMF_UTILS
#endif
