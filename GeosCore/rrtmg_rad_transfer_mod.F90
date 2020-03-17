#ifdef RRTMG
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: rrtmg_rad_transfer_mod.F90
!
! !DESCRIPTION: Module RRTMG\_RAD\_TRANSFER\_MOD contains arrays and routines
!  for performing online radiative transfer in GEOS-Chem.
!\\
!\\
! !INTERFACE:
!
MODULE RRTMG_RAD_TRANSFER_MOD
!
! !USES:
!
  USE CMN_FJX_MOD,  ONLY : RTODAER, RTSSAER, RTASYMAER, WVAA, SPECMASK
  USE CMN_SIZE_MOD, ONLY : NDUST, NAER
#ifdef BPCH_DIAG
  USE DIAG_MOD,     ONLY : AD72 !RAD OUTPUT DIAGNOSTIC ARRAY
#endif
  USE OMP_LIB
  USE PARRRTM,      ONLY : NBNDLW
  USE PARRRSW,      ONLY : NBNDSW

  IMPLICIT NONE
!
! !PUBLIC MEMBER FUNCTIONS
!
  PUBLIC  :: Do_RRTMG_Rad_Transfer
  PUBLIC  :: Init_RRTMG_Rad_Transfer
  PUBLIC  :: Cleanup_RRTMG_Rad_Transfer
!
! !PRIVATE:
!
  PRIVATE :: Init_MCICA_Clouds
  PRIVATE :: Init_Surface_Rad
  PRIVATE :: AttachPointersFromHemco
!
! !PUBLIC DATA MEMBERS:
!
  ! NOTE: Changed to pointers to get inputs from HEMCO (bmy, 10/30/18)
  ! NOTE: These should eventually go into fields of State_Chm
  REAL*4,  POINTER,     PUBLIC         :: MODIS_ALBDFNIR(:,:  )
  REAL*4,  POINTER,     PUBLIC         :: MODIS_ALBDFVIS(:,:  )
  REAL*4,  POINTER,     PUBLIC         :: MODIS_ALBDRNIR(:,:  )
  REAL*4,  POINTER,     PUBLIC         :: MODIS_ALBDRVIS(:,:  )
  REAL*4,  POINTER,     PUBLIC         :: MODIS_EMISS_01(:,:  )
  REAL*4,  POINTER,     PUBLIC         :: MODIS_EMISS_02(:,:  )
  REAL*4,  POINTER,     PUBLIC         :: MODIS_EMISS_03(:,:  )
  REAL*4,  POINTER,     PUBLIC         :: MODIS_EMISS_04(:,:  )
  REAL*4,  POINTER,     PUBLIC         :: MODIS_EMISS_05(:,:  )
  REAL*4,  POINTER,     PUBLIC         :: MODIS_EMISS_06(:,:  )
  REAL*4,  POINTER,     PUBLIC         :: MODIS_EMISS_07(:,:  )
  REAL*4,  POINTER,     PUBLIC         :: MODIS_EMISS_08(:,:  )
  REAL*4,  POINTER,     PUBLIC         :: MODIS_EMISS_09(:,:  )
  REAL*4,  POINTER,     PUBLIC         :: MODIS_EMISS_10(:,:  )
  REAL*4,  POINTER,     PUBLIC         :: MODIS_EMISS_11(:,:  )
  REAL*4,  POINTER,     PUBLIC         :: MODIS_EMISS_12(:,:  )
  REAL*4,  POINTER,     PUBLIC         :: MODIS_EMISS_13(:,:  )
  REAL*4,  POINTER,     PUBLIC         :: MODIS_EMISS_14(:,:  )
  REAL*4,  POINTER,     PUBLIC         :: MODIS_EMISS_15(:,:  )
  REAL*4,  POINTER,     PUBLIC         :: MODIS_EMISS_16(:,:  )
  REAL*4,  POINTER,     PUBLIC         :: CH4CLIM       (:,:,:)
  REAL*4,  POINTER,     PUBLIC         :: N2OCLIM       (:,:,:)
  REAL*4,  POINTER,     PUBLIC         :: CFC11CLIM     (:,:,:)
  REAL*4,  POINTER,     PUBLIC         :: CFC12CLIM     (:,:,:)
  REAL*4,  POINTER,     PUBLIC         :: CCL4CLIM      (:,:,:)
  REAL*4,  POINTER,     PUBLIC         :: CFC22CLIM     (:,:,:)

  !MCICA cloud variables now stored for reuse
  !NOTE: These should eventually go into fields of State_Chm
  REAL*8,  ALLOCATABLE, PUBLIC, TARGET :: CLDFMCL_LW(:,:,:,:)
  REAL*8,  ALLOCATABLE, PUBLIC, TARGET :: CIWPMCL_LW(:,:,:,:)
  REAL*8,  ALLOCATABLE, PUBLIC, TARGET :: CLWPMCL_LW(:,:,:,:)
  REAL*8,  ALLOCATABLE, PUBLIC, TARGET :: TAUCMCL_LW(:,:,:,:)
  REAL*8,  ALLOCATABLE, PUBLIC, TARGET :: CLDFMCL_SW(:,:,:,:)
  REAL*8,  ALLOCATABLE, PUBLIC, TARGET :: CIWPMCL_SW(:,:,:,:)
  REAL*8,  ALLOCATABLE, PUBLIC, TARGET :: CLWPMCL_SW(:,:,:,:)
  REAL*8,  ALLOCATABLE, PUBLIC, TARGET :: TAUCMCL_SW(:,:,:,:)
  REAL*8,  ALLOCATABLE, PUBLIC, TARGET :: SSACMCL   (:,:,:,:)
  REAL*8,  ALLOCATABLE, PUBLIC, TARGET :: ASMCMCL   (:,:,:,:)
  REAL*8,  ALLOCATABLE, PUBLIC, TARGET :: FSFCMCL   (:,:,:,:)
  REAL*8,  ALLOCATABLE, PUBLIC, TARGET :: REICMCL   (:,:,:  )
  REAL*8,  ALLOCATABLE, PUBLIC, TARGET :: RELQMCL   (:,:,:  )
!
! !REVISION HISTORY:
!  18 Jun 2013 - D.A. Ridley - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
  REAL*8,  ALLOCATABLE  :: LW_UFLUX (:,:,:)
  REAL*8,  ALLOCATABLE  :: LW_DFLUX (:,:,:)
  REAL*8,  ALLOCATABLE  :: SW_UFLUX (:,:,:)
  REAL*8,  ALLOCATABLE  :: SW_DFLUX (:,:,:)
  REAL*8,  ALLOCATABLE  :: LW_UFLUXC(:,:,:)
  REAL*8,  ALLOCATABLE  :: LW_DFLUXC(:,:,:)
  REAL*8,  ALLOCATABLE  :: SW_UFLUXC(:,:,:)
  REAL*8,  ALLOCATABLE  :: SW_DFLUXC(:,:,:)

  REAL*8  :: RRTMG_LMB(NBNDLW+NBNDSW)

  INTEGER :: ID_AER_LMB0 (NBNDLW+NBNDSW)
  INTEGER :: ID_AER_LMB1 (NBNDLW+NBNDSW)

CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: do_rrtmg_rad_transfer
!
! !DESCRIPTION: Interface between GEOS-Chem and the RRTMG radiative
!  transfer model.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE DO_RRTMG_RAD_TRANSFER( ThisDay,    ThisMonth,  &
                                    iCld,       iSpecMenu,  &
                                    iNcDiag,    iSeed,      &
                                    Input_Opt,  State_Chm,  &
                                    State_Diag, State_Grid, &
                                    State_Met,  RC )
!
! !USES:
!
    !-----------------------------------------------------------------
    ! Modules from GeosRad
    !-----------------------------------------------------------------
    USE MCICA_SUBCOL_GEN_LW, ONLY : MCICA_SUBCOL_LW
    USE MCICA_SUBCOL_GEN_SW, ONLY : MCICA_SUBCOL_SW
    USE PARKIND,             ONLY : IM=>KIND_IM, RB=>KIND_RB
    USE RRLW_CON,            ONLY : GASCON, AVOGAD
    USE PARRRTM,             ONLY : NBNDLW, NGPTLW
    USE PARRRSW,             ONLY : NBNDSW, NGPTSW,NAEREC
    USE RRTMG_LW_RAD,        ONLY : RRTMG_LW
    USE RRTMG_SW_RAD,        ONLY : RRTMG_SW

    !-----------------------------------------------------------------
    ! GEOS-Chem modules
    !-----------------------------------------------------------------
    USE CMN_FJX_MOD,         ONLY : NSPECRAD  ! NUMBER OF SPECIES FOR RT
    USE CMN_FJX_MOD,         ONLY : NASPECRAD ! NUMBER OF AEROSOL SPECIES
    USE CMN_FJX_MOD,         ONLY : SPECMASK,   IRTWVSELECT
    USE CMN_FJX_MOD,         ONLY : ACOEF_RTWV, BCOEF_RTWV, CCOEF_RTWV
    USE CMN_FJX_MOD,         ONLY : WVAA,       NWVAA
    USE CMN_FJX_MOD,         ONLY : NWVAA0
    USE ErrCode_Mod
    USE ERROR_MOD
    USE Input_Opt_Mod,       ONLY : OptInput
    USE PhysConstants,       ONLY : AIRMW, PI, AVO
    USE PRESSURE_MOD,        ONLY : GET_PCENTER,      GET_PEDGE
    USE State_Chm_Mod,       ONLY : ChmState
    USE State_Chm_Mod,       ONLY : Ind_
    USE State_Diag_Mod,      ONLY : DgnState
    USE State_Grid_Mod,      ONLY : GrdState
    USE State_Met_Mod,       ONLY : MetState
    USE TIME_MOD,            ONLY : GET_DAY_OF_YEAR, GET_HOUR
    USE TOMS_MOD,            ONLY : GET_OVERHEAD_O3
    USE UnitConv_Mod,        ONLY : Convert_Spc_Units
#ifdef BPCH_DIAG
    USE DIAG_MOD,            ONLY : AD72
#endif
!
! !INPUT PARAMETERS:
!
    INTEGER,        INTENT(IN)    :: ThisDay    ! CURRENT DAY
    INTEGER,        INTENT(IN)    :: ThisMonth  ! CURRENT MONTH
    INTEGER,        INTENT(IN)    :: iSpecMenu  ! THE SPECIES BEING INCLUDED
                                                ! NEEDED FOR OUTPUT PURPOSES
    INTEGER,        INTENT(IN)    :: iNcDiag    ! Index for netCDF diag arrays
    INTEGER,        INTENT(IN)    :: iSeed      ! Seed value
    TYPE(OptInput), INTENT(IN)    :: Input_Opt  ! Input Options object
    TYPE(GrdState), INTENT(IN)    :: State_Grid ! Grid State object
    TYPE(MetState), INTENT(IN)    :: State_Met  ! Meteorology State object
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ChmState), INTENT(INOUT) :: State_Chm  ! Chemistry State object
    TYPE(DgnState), INTENT(INOUT) :: State_Diag ! Diagnostics State object
    INTEGER,        INTENT(INOUT) :: iCld       ! CLOUD FLAG FOR RRTMG
                                                  ! 0-NOCLOUD, 1-GREY CLOUD
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT) :: RC          ! Failure or success
!
! !REMARKS:
!  THIS ROUTINE PASSES INPUTS TO THE RRTMG DRIVER ROUTINE "RAD_DRIVER"
!  VIA THE ARGUMENT LIST.  THIS PREVENTS CIRCULAR REFERENCES.
!
! !REVISION HISTORY:
!  17 AUG 2012 - R. YANTOSCA - INITIAL VERSION
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS
!
    INTEGER, PARAMETER :: NWV=37
!
! !LOCAL VARIABLES:
!
    LOGICAL            :: DOAERAD      ! INCLUDE AEROSOL DETERMINED FROM
                                           ! SPECMASK
    LOGICAL            :: LOUTPUTAERO  ! OUTPUT AEROSOL DIAGNOSTICS?
    INTEGER            :: NAD72        ! NUMBER OF OUTPUTS PER FIELD
    INTEGER            :: ITIMEVALS(8)
    INTEGER            :: IDIAGOUT     ! INDEX OF SPC OPTICS FOR OUTPUT
    REAL*8             :: OLDSECS, NEWSECS
    CHARACTER(LEN=63)  :: OrigUnit

    ! SAVEd scalars
    LOGICAL, SAVE      :: FIRST = .TRUE.
    INTEGER, SAVE      :: id_O3,    id_CH4,  id_N2O, id_CFC11
    INTEGER, SAVE      :: id_CFC12, id_CCL4, id_HCFC22

    !-----------------------------------------------------------------
    ! TEMPORARY AEROSOL VARIABLES
    !-----------------------------------------------------------------
    REAL*8             :: LAMBDA(NWV)

    !-----------------------------------------------------------------
    ! SCALARS
    !-----------------------------------------------------------------
    INTEGER            :: I, J ,L, LL, N, W
    INTEGER            :: IL, OMPID, LCHEM
    INTEGER            :: OUTIDX,IOUTWV
    INTEGER            :: IB,IBX,IB_SW,IS,NBNDS,NSPEC
    INTEGER            :: IS_ON,NASPECRAD_ON
    INTEGER            :: IASPECRAD_ON(NASPECRAD)
    REAL*8             :: RHOICE=0.9167, RHOLIQ=1.    ! G/CM3

    !-----------------------------------------------------------------
    ! REL AND REI FROM PERSONAL COMMUNICATION FROM LAZAROS OREOPOULOS
    ! (GSFC) 12/12/12
    !-----------------------------------------------------------------
    REAL*8             :: REL_DEF = 14.2, REI_DEF=24.8    ! MICRONS
    INTEGER            :: DOY

    INTEGER            :: IHR
    CHARACTER(LEN=2)   :: CHR

    !-----------------------------------------------------------------
    ! ARRAYS FROM GC
    !-----------------------------------------------------------------
    REAL(KIND=RB)      :: PCENTER(State_Grid%NX,State_Grid%NY,State_Grid%NZ)
    REAL(KIND=RB)      :: PEDGE  (State_Grid%NX,State_Grid%NY,State_Grid%NZ)

    !-----------------------------------------------------------------
    !ARRAYS FOR RRTMG
    !-----------------------------------------------------------------
    REAL(KIND=RB)      :: O3VMR (State_Grid%NX,State_Grid%NY,State_Grid%NZ)
    REAL(KIND=RB)      :: CH4VMR(State_Grid%NX,State_Grid%NY,State_Grid%NZ)
    REAL(KIND=RB)      :: N2OVMR(State_Grid%NX,State_Grid%NY,State_Grid%NZ)
    REAL(KIND=RB)      :: CFC11VMR(State_Grid%NX,State_Grid%NY,State_Grid%NZ)
    REAL(KIND=RB)      :: CFC12VMR(State_Grid%NX,State_Grid%NY,State_Grid%NZ)
    REAL(KIND=RB)      :: CFC22VMR(State_Grid%NX,State_Grid%NY,State_Grid%NZ)
    REAL(KIND=RB)      :: CCL4VMR(State_Grid%NX,State_Grid%NY,State_Grid%NZ)

    REAL(KIND=RB)      :: TAUCLD(State_Grid%NX,State_Grid%NY,State_Grid%NZ)
    REAL(KIND=RB)      :: CLDFR(State_Grid%NX,State_Grid%NY,State_Grid%NZ)
    REAL(KIND=RB)      :: RELIQ(State_Grid%NX,State_Grid%NY,State_Grid%NZ)
    REAL(KIND=RB)      :: REICE(State_Grid%NX,State_Grid%NY,State_Grid%NZ)
    REAL(KIND=RB)      :: CLIQWP(State_Grid%NX,State_Grid%NY,State_Grid%NZ)
    REAL(KIND=RB)      :: CICEWP(State_Grid%NX,State_Grid%NY,State_Grid%NZ)
    REAL(KIND=RB)      :: ALBVIS(State_Grid%NX,State_Grid%NY)

    REAL(KIND=RB)      :: TAUAER_LW(State_Grid%NX,State_Grid%NY,State_Grid%NZ,&
                                    NBNDLW)
    REAL(KIND=RB)      :: TAUAER_SW(State_Grid%NX,State_Grid%NY,State_Grid%NZ,&
                                    NBNDSW)
    REAL(KIND=RB)      :: SSAAER(State_Grid%NX,State_Grid%NY,State_Grid%NZ,&
                                 NBNDSW)
    REAL(KIND=RB)      :: ASMAER(State_Grid%NX,State_Grid%NY,State_Grid%NZ,&
                                 NBNDSW)

    !-----------------------------------------------------------------
    !TO STORE THE OPTICS FOR THE AEROSOLS WE ARE INTERESTED IN
    !-----------------------------------------------------------------
    REAL*8             :: TAUAERDIAG(State_Grid%NX,State_Grid%NY,State_Grid%NZ,&
                                     NBNDSW)
    REAL*8             :: SSAAERDIAG(State_Grid%NX,State_Grid%NY,State_Grid%NZ,&
                                     NBNDSW)
    REAL*8             :: ASMAERDIAG(State_Grid%NX,State_Grid%NY,State_Grid%NZ,&
                                     NBNDSW)

    REAL(KIND=RB)      :: H2OVMR   (State_Grid%NX,State_Grid%NY,State_Grid%NZ)
    REAL(KIND=RB)      :: TLAY     (State_Grid%NX,State_Grid%NY,State_Grid%NZ)
    REAL(KIND=RB)      :: PLAY     (State_Grid%NX,State_Grid%NY,State_Grid%NZ)
    REAL(KIND=RB)      :: SUNCOS   (State_Grid%NX,State_Grid%NY,State_Grid%NZ)
    REAL(KIND=RB)      :: TSFC     (State_Grid%NX,State_Grid%NY)

    !-----------------------------------------------------------------
    !VARIABLES FOR OBTAINING STRATOSPHERIC VARIABLES
    !-----------------------------------------------------------------
    REAL*8             :: O3COL, YLAT, AIR_TMP

    !-----------------------------------------------------------------
    !SURFACE
    !-----------------------------------------------------------------
    REAL(KIND=RB)      :: ALBDIRVIS(State_Grid%NX,State_Grid%NY)
    REAL(KIND=RB)      :: ALBDIFVIS(State_Grid%NX,State_Grid%NY)
    REAL(KIND=RB)      :: ALBDIRNIR(State_Grid%NX,State_Grid%NY)
    REAL(KIND=RB)      :: ALBDIFNIR(State_Grid%NX,State_Grid%NY)
    REAL(KIND=RB)      :: RTEMISS  (State_Grid%NX,State_Grid%NY, NBNDLW)

    REAL*8             :: NUMER,DENOM
    REAL*4             :: AODTMP,AODOUT,SSATMP,SSAOUT
    REAL*4             :: ASYMTMP,ASYMOUT
    INTEGER            :: FLG_FIRST_STRAT(State_Grid%NX,State_Grid%NY)
    INTEGER            :: ONECOL
    REAL*4             :: CH4SCL(State_Grid%NX,State_Grid%NY)

    !-----------------------------------------------------------------
    ! FROM RAD_DRIVER... TO BE MERGED
    ! FLAGS AND DIMENSIONS
    !-----------------------------------------------------------------
    INTEGER (KIND=IM)  :: IDRV
    INTEGER (KIND=IM)  :: INFLGLW, ICEFLGLW,LIQFLGLW
    INTEGER (KIND=IM)  :: INFLGSW, ICEFLGSW,LIQFLGSW

    !-----------------------------------------------------------------
    ! PROFILE VARIABLES
    !-----------------------------------------------------------------
    REAL (KIND=RB)     :: PLEV(State_Grid%NX,State_Grid%NY,State_Grid%NZ+1)
    REAL (KIND=RB)     :: TLEV(State_Grid%NX,State_Grid%NY,State_Grid%NZ+1)
    REAL (KIND=RB)     :: CO2VMR(State_Grid%NX,State_Grid%NY,State_Grid%NZ)
    REAL (KIND=RB)     :: O2VMR(State_Grid%NX,State_Grid%NY,State_Grid%NZ)
    REAL (KIND=RB)     :: T_CTM(State_Grid%NZ+1)
    REAL (KIND=RB)     :: P_CTM(State_Grid%NZ+2)
    REAL (KIND=RB)     :: O3_CTM(State_Grid%NZ+1)
    REAL (KIND=RB)     :: T_CLIM(State_Grid%NZ+1)
    REAL (KIND=RB)     :: O3_CLIM(State_Grid%NZ+1)
    REAL (KIND=RB)     :: Z_CLIM(State_Grid%NZ+2)
    REAL (KIND=RB)     :: AIR_CLIM(State_Grid%NZ+1)

    !-----------------------------------------------------------------
    ! SW SOLAR VARIABLES
    !-----------------------------------------------------------------
    REAL(KIND=RB)      :: ADJES=1.0     ! FLUX ADJUSTMENT FOR EARTH/SUN DIST
    REAL(KIND=RB)      :: SCON=1368.22  ! SOLAR CONSTANT (W/M2)

    !-----------------------------------------------------------------
    ! SW CLOUD VARIABLES
    !-----------------------------------------------------------------
    REAL(KIND=RB)      :: TAUCLD_SW(NBNDSW,State_Grid%NX,State_Grid%NY,State_Grid%NZ) ! IN-CLOUD OPTICAL DEPTH
    REAL(KIND=RB)      :: TAUCLD_LW(NBNDLW,State_Grid%NX,State_Grid%NY,State_Grid%NZ) ! NOT USED BUT PASSED TO MCICA_LW
    REAL(KIND=RB)      :: SSACLD(NBNDSW,State_Grid%NX,State_Grid%NY,State_Grid%NZ)    ! IN-CLOUD SINGLE SCATTERING ALBEDO
    REAL(KIND=RB)      :: ASMCLD(NBNDSW,State_Grid%NX,State_Grid%NY,State_Grid%NZ)    ! IN-CLOUD ASYMMETRY PARAMETER
    REAL(KIND=RB)      :: FSFCLD(NBNDSW,State_Grid%NX,State_Grid%NY,State_Grid%NZ)    ! IN-CLOUD FORWARD SCATTERING FRACTION
    REAL(KIND=RB)      :: ECAER(1,State_Grid%NZ,NAEREC)               ! AEROSOL OPTICAL DEPTH AT 0.55UM (IAER=6 ONLY)

    !-----------------------------------------------------------------
    ! LONGWAVE FLUX VARIABLES
    !-----------------------------------------------------------------
    REAL(KIND=RB)      :: UFLX(1,State_Grid%NZ+1)  ! TOTAL SKY LONGWAVE UPWARD FLUX (W/M2)
    REAL(KIND=RB)      :: DFLX(1,State_Grid%NZ+1)  ! TOTAL SKY LONGWAVE DOWNWARD FLUX (W/M2)
    REAL(KIND=RB)      :: HR(1,State_Grid%NZ)      ! TOTAL SKY LONGWAVE RADIATIVE HEATING RATE (K/D)
    REAL(KIND=RB)      :: UFLXC(1,State_Grid%NZ+1) ! CLEAR SKY LONGWAVE UPWARD FLUX (W/M2)
    REAL(KIND=RB)      :: DFLXC(1,State_Grid%NZ+1) ! CLEAR SKY LONGWAVE DOWNWARD FLUX (W/M2)
    REAL(KIND=RB)      :: HRC(1,State_Grid%NZ)     ! CLEAR SKY LONGWAVE RADIATIVE HEATING RATE (K/D)

    !-----------------------------------------------------------------
    !- OPTIONAL OUTPUT
    !-----------------------------------------------------------------
    REAL(KIND=RB)      :: DUFLX_DT(1,State_Grid%NZ) ! CHANGE IN UPWARD LONGWAVE FLUX (W/M2/K)
    REAL(KIND=RB)      :: DUFLXC_DT(1,State_Grid%NZ)! CHANGE IN CLEAR SKY UPWARD LONGWAVE FLUX (W/M2/K)

    !-----------------------------------------------------------------
    ! SHORTWAVE FLUX VARIABLES
    !-----------------------------------------------------------------
    ! ----- OUTPUT -----
    REAL(KIND=RB)      :: SWUFLX(1,State_Grid%NZ+1)  ! TOTAL SKY SHORTWAVE UPWARD FLUX (W/M2)
    REAL(KIND=RB)      :: SWDFLX(1,State_Grid%NZ+1)  ! TOTAL SKY SHORTWAVE DOWNWARD FLUX (W/M2)
    REAL(KIND=RB)      :: SWHR(1,State_Grid%NZ)      ! TOTAL SKY SHORTWAVE RADIATIVE HEATING RATE (K/D)
    REAL(KIND=RB)      :: SWUFLXC(1,State_Grid%NZ+1) ! CLEAR SKY SHORTWAVE UPWARD FLUX (W/M2)
    REAL(KIND=RB)      :: SWDFLXC(1,State_Grid%NZ+1) ! CLEAR SKY SHORTWAVE DOWNWARD FLUX (W/M2)
    REAL(KIND=RB)      :: SWHRC(1,State_Grid%NZ)     ! CLEAR SKY SHORTWAVE RADIATIVE HEATING RATE (K/D)

    !-----------------------------------------------------------------
    ! LOCAL VARIABLES
    !-----------------------------------------------------------------
    REAL*8             :: GCAIR
    REAL*8             :: RHOA, RHOB, RHOSUM
    REAL*8             :: HR_TEMP

    !-----------------------------------------------------------------
    ! MCICA VARIABLES
    !-----------------------------------------------------------------
    INTEGER(KIND=IM)   :: SEEDSW, SEEDLW
    INTEGER(KIND=IM)   :: IRNG=1  ! MERSENNE TWISTER RANDOM NUMBER GENERATOR
    INTEGER(KIND=IM)   :: ICLDMCL
    REAL(KIND=RB)      :: RELQMCL0(1,State_Grid%NZ)
    REAL(KIND=RB)      :: REICMCL0(1,State_Grid%NZ)

    !-----------------------------------------------------------------
    ! MCICA LW SPECIFIC
    !-----------------------------------------------------------------
    REAL(KIND=RB)      :: CLDFMCL_LW0(NGPTLW,1,State_Grid%NZ)
    REAL(KIND=RB)      :: CIWPMCL_LW0(NGPTLW,1,State_Grid%NZ)
    REAL(KIND=RB)      :: CLWPMCL_LW0(NGPTLW,1,State_Grid%NZ)
    REAL(KIND=RB)      :: TAUCMCL_LW0(NGPTLW,1,State_Grid%NZ)

    !-----------------------------------------------------------------
    ! MCICA SW SPECIFIC
    !-----------------------------------------------------------------
    REAL(KIND=RB)      :: CLDFMCL_SW0(NGPTSW,1,State_Grid%NZ)
    REAL(KIND=RB)      :: CIWPMCL_SW0(NGPTSW,1,State_Grid%NZ)
    REAL(KIND=RB)      :: CLWPMCL_SW0(NGPTSW,1,State_Grid%NZ)
    REAL(KIND=RB)      :: TAUCMCL_SW0(NGPTSW,1,State_Grid%NZ)
    REAL(KIND=RB)      :: SSACMCL0(NGPTSW,1,State_Grid%NZ)
    REAL(KIND=RB)      :: ASMCMCL0(NGPTSW,1,State_Grid%NZ)
    REAL(KIND=RB)      :: FSFCMCL0(NGPTSW,1,State_Grid%NZ)
    REAL(KIND=RB)      :: PCENTER0(1,State_Grid%NZ)
    REAL(KIND=RB)      :: CLDFR0(1,State_Grid%NZ)
    REAL(KIND=RB)      :: CICEWP0(1,State_Grid%NZ)
    REAL(KIND=RB)      :: CLIQWP0(1,State_Grid%NZ)
    REAL(KIND=RB)      :: REICE0(1,State_Grid%NZ)
    REAL(KIND=RB)      :: RELIQ0(1,State_Grid%NZ)
    REAL(KIND=RB)      :: TAUCLD_SW0(NBNDSW,1,State_Grid%NZ)
    REAL(KIND=RB)      :: SSACLD0(NBNDSW,1,State_Grid%NZ)
    REAL(KIND=RB)      :: ASMCLD0(NBNDSW,1,State_Grid%NZ)
    REAL(KIND=RB)      :: FSFCLD0(NBNDSW,1,State_Grid%NZ)

    !-----------------------------------------------------------------
    ! Variables used to avoid array temporaries (bmy, 6/3/15)
    !
    ! NOTE: Use temporary arrays instead of pointers.  For unknown
    ! reasons the pointer references incur segfaults. (bmy, 6/3/15)
    !-----------------------------------------------------------------

    ! For MCICA_SUBCOL_LW and MCICA_SUBCOL_LW
    REAL(KIND=RB)      :: p_PCENTER   (         State_Grid%NZ  )
    REAL(KIND=RB)      :: p_CLDFR     (         State_Grid%NZ  )
    REAL(KIND=RB)      :: p_CICEWP    (         State_Grid%NZ  )
    REAL(KIND=RB)      :: p_CLIQWP    (         State_Grid%NZ  )
    REAL(KIND=RB)      :: p_REICE     (         State_Grid%NZ  )
    REAL(KIND=RB)      :: p_RELIQ     (         State_Grid%NZ  )
    REAL(KIND=RB)      :: p_TAUCLD_LW ( NBNDLW, State_Grid%NZ  )
    REAL(KIND=RB)      :: p_TAUCLD_SW ( NBNDSW, State_Grid%NZ  )
    REAL(KIND=RB)      :: p_SSACLD    ( NBNDSW, State_Grid%NZ  )
    REAL(KIND=RB)      :: p_ASMCLD    ( NBNDSW, State_Grid%NZ  )
    REAL(KIND=RB)      :: p_FSFCLD    ( NBNDSW, State_Grid%NZ  )

    ! For RRTMG_LW and RRTMG_SW
    REAL(KIND=RB)      :: p_PLEV      (         State_Grid%NZ+1)
    REAL(KIND=RB)      :: p_TLAY      (         State_Grid%NZ  )
    REAL(KIND=RB)      :: p_TLEV      (         State_Grid%NZ+1)
    REAL(KIND=RB)      :: p_H2OVMR    (         State_Grid%NZ  )
    REAL(KIND=RB)      :: p_O3VMR     (         State_Grid%NZ  )
    REAL(KIND=RB)      :: p_CO2VMR    (         State_Grid%NZ  )
    REAL(KIND=RB)      :: p_CH4VMR    (         State_Grid%NZ  )
    REAL(KIND=RB)      :: p_N2OVMR    (         State_Grid%NZ  )
    REAL(KIND=RB)      :: p_O2VMR     (         State_Grid%NZ  )
    REAL(KIND=RB)      :: p_CFC11VMR  (         State_Grid%NZ  )
    REAL(KIND=RB)      :: p_CFC12VMR  (         State_Grid%NZ  )
    REAL(KIND=RB)      :: p_CFC22VMR  (         State_Grid%NZ  )
    REAL(KIND=RB)      :: p_CCL4VMR   (         State_Grid%NZ  )
    REAL(KIND=RB)      :: p_REICMCL   (         State_Grid%NZ  )
    REAL(KIND=RB)      :: p_RELQMCL   (         State_Grid%NZ  )
    REAL(KIND=RB)      :: p_SUNCOS    (         State_Grid%NZ  )
    REAL(KIND=RB)      :: p_CLDFMCL_LW( NGPTLW, State_Grid%NZ  )
    REAL(KIND=RB)      :: p_TAUCMCL_LW( NGPTLW, State_Grid%NZ  )
    REAL(KIND=RB)      :: p_CIWPMCL_LW( NGPTLW, State_Grid%NZ  )
    REAL(KIND=RB)      :: p_CLWPMCL_LW( NGPTLW, State_Grid%NZ  )
    REAL(KIND=RB)      :: p_CLDFMCL_SW( NGPTSW, State_Grid%NZ  )
    REAL(KIND=RB)      :: p_TAUCMCL_SW( NGPTSW, State_Grid%NZ  )
    REAL(KIND=RB)      :: p_SSACMCL   ( NGPTSW, State_Grid%NZ  )
    REAL(KIND=RB)      :: p_ASMCMCL   ( NGPTSW, State_Grid%NZ  )
    REAL(KIND=RB)      :: p_FSFCMCL   ( NGPTSW, State_Grid%NZ  )
    REAL(KIND=RB)      :: p_CIWPMCL_SW( NGPTSW, State_Grid%NZ  )
    REAL(KIND=RB)      :: p_CLWPMCL_SW( NGPTSW, State_Grid%NZ  )
    REAL(KIND=RB)      :: p_RTEMISS   (         NBNDLW )
    REAL(KIND=RB)      :: p_TAUAER_LW ( State_Grid%NZ,  NBNDLW )
    REAL(KIND=RB)      :: p_TAUAER_SW ( State_Grid%NZ,  NBNDSW )
    REAL(KIND=RB)      :: p_SSAAER    ( State_Grid%NZ,  NBNDSW )
    REAL(KIND=RB)      :: p_ASMAER    ( State_Grid%NZ,  NBNDSW )

    ! Strings
    CHARACTER(LEN=255) :: ErrMsg, ThisLoc

    !=================================================================
    ! DO_RRTMG_RAD_TRANSFER begins here!
    !=================================================================

    ! Initialize
    RC      = GC_SUCCESS
    ErrMsg  = ''
    ThisLoc = ' -> at DO_RRTMG_RAD_TRANSFER (in rrtmg_rad_transfer_mod.F90)'

    ! Convert species units to kg/kg dry for RRTMG
    CALL Convert_Spc_Units( Input_Opt, State_Chm, State_Grid, State_Met, &
                            'kg/kg dry', RC, OrigUnit=OrigUnit )

    ! Trap potential errors
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Unit conversion error in DO_RRTMG_RAD_TRANSFER!"'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    ! Also make sure that the ncDiag arguement is valid,
    ! as this is the index for the netCDF diagnostic arrays.
    IF ( iNcDiag <= 0 ) THEN
       ErrMsg = 'The iNcDiag argument is <= 0!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    !=================================================================
    ! PREPARE INPUTS FOR RAD_DRIVER
    !=================================================================

    !CALL GET_SPECIES( .TRUE., THISMONTH, 'O3',  GMI_O3  )
    !CALL GET_SPECIES( .TRUE., THISMONTH, 'CH4', GMI_CH4 )

    !=================================================================
    ! PREPARE INPUTS FOR RRTMG
    !=================================================================

    ! INITIALIZE
    NSPEC                = NAER+NDUST+4
    FLG_FIRST_STRAT(:,:) = 0 !FLAG TO DETERMINE IF THE FIRST STRATOSPHERIC
                             ! LEVEL HAS BEEN REACHED

    !DETERMINE IF WE ARE RUNNING WITH AEROSOL
    !CREATE INDEX FOR AEROSOLS REQUIRED
    LOUTPUTAERO     = .TRUE. !SET AEROSOL DIAG OUTPUT TO TRUE INITIALLY
    DOAERAD         = .FALSE.
    NASPECRAD_ON    = 0
    IASPECRAD_ON(:) = 0

    DO N=1,NASPECRAD
       IF (SPECMASK(N).GT.0) THEN
          DOAERAD = .TRUE.
          NASPECRAD_ON = NASPECRAD_ON +1
          !create list of species required and tag with index
          IASPECRAD_ON(NASPECRAD_ON) = N
          IDIAGOUT = MAX(IDIAGOUT,SPECMASK(N))
       ENDIF
    ENDDO
    !write(6,*) 'SPECMASK:',SPECMASK

    ! Initialize arrays
    TAUCLD(:,:,:)       = 0.0
    CLDFR(:,:,:)        = 0.0
    RELIQ(:,:,:)        = 0.0
    REICE(:,:,:)        = 0.0
    CLIQWP(:,:,:)       = 0.0
    CICEWP(:,:,:)       = 0.0
    TAUAER_LW(:,:,:,:)  = 0.0
    TAUAER_SW(:,:,:,:)  = 0.0
    SSAAER(:,:,:,:)     = 0.0
    ASMAER(:,:,:,:)     = 0.0
    TAUAERDIAG(:,:,:,:) = 0.0D0
    SSAAERDIAG(:,:,:,:) = 0.0D0
    ASMAERDIAG(:,:,:,:) = 0.0D0
    UFLX(:,:)           = 0.0
    DFLX(:,:)           = 0.0
    HR(:,:)             = 0.0
    UFLXC(:,:)          = 0.0
    DFLXC(:,:)          = 0.0
    HRC(:,:)            = 0.0
    DUFLX_DT(:,:)       = 0.0
    DUFLXC_DT(:,:)      = 0.0
    SWUFLX(:,:)         = 0.0
    SWDFLX(:,:)         = 0.0
    SWHR(:,:)           = 0.0
    SWUFLXC(:,:)        = 0.0
    SWDFLXC(:,:)        = 0.0
    SWHRC(:,:)          = 0.0
    O3VMR(:,:,:)        = 0.0
    CH4VMR(:,:,:)       = 0.0
    NBNDS               = NBNDLW+NBNDSW

    !=================================================================
    ! First-time setup
    !=================================================================
    IF ( FIRST ) THEN

       ! Define species ID flags
       id_O3     = Ind_('O3')
       id_CH4    = Ind_('CH4')
       id_N2O    = Ind_('N2O')
       id_CFC11  = Ind_('CFC11')
       id_CFC12  = Ind_('CFC12')
       id_CCL4   = Ind_('CCL4')
       id_HCFC22 = Ind_('HCFC22')

       ! Get pointers to data fields that are read by HEMCO
       ! NOTE: This has to be done here and not in initialization
       ! because we have to wait for HEMCO to read the data from disk,
       ! which is done after initialization. (bmy, 10/30/18)
       CALL AttachPointersFromHemco( RC )
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Error encountered in "AttachPointersFromHemco"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

       ! Reset first-time flag
       FIRST = .FALSE.
    ENDIF

    !$OMP PARALLEL DO          &
    !$OMP DEFAULT( SHARED )    &
    !$OMP PRIVATE( I, J, IB  ) &
    !$OMP SCHEDULE( DYNAMIC )
    DO J = 1, State_Grid%NY
    DO I = 1, State_Grid%NX

       ! BROADBAND ALBEDO
       ! Attache the HEMCO pointers to RRTMG albedo variables (bmy, 11/1/18)
       ALBDIRVIS(I,J) = MODIS_ALBDRVIS(I,J)
       ALBDIFVIS(I,J) = MODIS_ALBDFVIS(I,J)
       ALBDIRNIR(I,J) = MODIS_ALBDRNIR(I,J)
       ALBDIFNIR(I,J) = MODIS_ALBDFNIR(I,J)

       IF ((ALBDIRVIS(I,J).GT.0.999)  .OR. &
           (ALBDIRVIS(I,J).LT.0.001)) THEN
          WRITE(6,*) 'ALBEDO DRVIS OUT OF RANGE',I,J,ALBDIRVIS(I,J)
       ENDIF
       IF ((ALBDIFVIS(I,J).GT.0.999)  .OR. &
           (ALBDIFVIS(I,J).LT.0.001)) THEN
          WRITE(6,*) 'ALBEDO DFVIS OUT OF RANGE',I,J,ALBDIFVIS(I,J)
       ENDIF
       IF ((ALBDIRNIR(I,J).GT.0.999)  .OR. &
           (ALBDIRNIR(I,J).LT.0.001)) THEN
          WRITE(6,*) 'ALBEDO DRNIR OUT OF RANGE',I,J,ALBDIRNIR(I,J)
       ENDIF
       IF ((ALBDIFNIR(I,J).GT.0.999)  .OR. &
           (ALBDIFNIR(I,J).LT.0.001)) THEN
          WRITE(6,*) 'ALBEDO DFNIR OUT OF RANGE',I,J,ALBDIFNIR(I,J)
       ENDIF

       ! Assign the MODIS emissivity pointers from HEMCO to the
       ! different slots of the RTEMISS array (for each spectral band)
       ! Hardcode the assignments, which is much faster (bmy, 11/2/18)
       RTEMISS(I,J,1 ) = MODIS_EMISS_01(I,J)
       RTEMISS(I,J,2 ) = MODIS_EMISS_02(I,J)
       RTEMISS(I,J,3 ) = MODIS_EMISS_03(I,J)
       RTEMISS(I,J,4 ) = MODIS_EMISS_04(I,J)
       RTEMISS(I,J,5 ) = MODIS_EMISS_05(I,J)
       RTEMISS(I,J,6 ) = MODIS_EMISS_06(I,J)
       RTEMISS(I,J,7 ) = MODIS_EMISS_07(I,J)
       RTEMISS(I,J,8 ) = MODIS_EMISS_08(I,J)
       RTEMISS(I,J,9 ) = MODIS_EMISS_09(I,J)
       RTEMISS(I,J,10) = MODIS_EMISS_10(I,J)
       RTEMISS(I,J,11) = MODIS_EMISS_11(I,J)
       RTEMISS(I,J,12) = MODIS_EMISS_12(I,J)
       RTEMISS(I,J,13) = MODIS_EMISS_13(I,J)
       RTEMISS(I,J,14) = MODIS_EMISS_14(I,J)
       RTEMISS(I,J,15) = MODIS_EMISS_15(I,J)
       RTEMISS(I,J,16) = MODIS_EMISS_16(I,J)

    ENDDO
    ENDDO
    !$OMP END PARALLEL DO

    !GET PCENTER, PEDGE AND DETERMINE IF IN TROP
    !%%% NOTE: LOOPS ARE GOING IN WRONG ORDER (bmy, 1/8/18)
    DO I = 1, State_Grid%NX
    DO J = 1, State_Grid%NY
       DO L = 1, State_Grid%NZ
          PCENTER(I,J,L) = GET_PCENTER( I, J, L )
          PEDGE  (I,J,L) = GET_PEDGE  ( I, J, L )
          H2OVMR (I,J,L) = State_Met%AVGW(I,J,L)
          TLAY   (I,J,L) = State_Met%T(I,J,L)
          SUNCOS (I,J,L) = State_Met%SUNCOS(I,J)
       ENDDO
       TSFC  (I,J)   = State_Met%TSKIN(I,J)
    ENDDO
    ENDDO

    !$OMP PARALLEL DO       &
    !$OMP DEFAULT( SHARED ) &
    !$OMP PRIVATE( I,       J,      L                         ) &
    !$OMP PRIVATE( AIR_TMP, YLAT,   O3COL,  O3_CTM,  T_CTM    ) &
    !$OMP PRIVATE( P_CTM,   T_CLIM, Z_CLIM, O3_CLIM, AIR_CLIM ) &
    !$OMP SCHEDULE( DYNAMIC )
    DO J = 1, State_Grid%NY
    DO I = 1, State_Grid%NX

       ALBVIS(I,J) = State_Met%ALBD(I,J)

       ! Grid box latitude [degrees]
       YLAT = State_Grid%YMid(I,J)

#if defined( EXTERNAL_GRID ) || defined( EXTERNAL_FORCING )
       ! Define the P array here, using GCM pressures
       DO L = 1, State_Grid%NZ+1
          P_CTM(L) = State_Met%PEDGE(I,J,L)
       ENDDO
#else
       ! Define the P array here, using Ap and Bp from GEOS-Chem
       DO L = 1, State_Grid%NZ+1
          P_CTM(L) = GET_PEDGE(I,J,L)
       ENDDO
#endif

       ! Top edge of P_CTM is top of atmosphere
       P_CTM(State_Grid%NZ+2)  = 0d0

       ! Temperature profile [K]
       T_CTM(1:State_Grid%NZ)  = State_Met%T(I,J,1:State_Grid%NZ)

       ! Top of atmosphere
       T_CTM(State_Grid%NZ+1)  = T_CTM(State_Grid%NZ)

       ! Overhead ozone column [DU]
       ! These values are either from the met fields or TOMS/SBUV,
       ! depending on the settings in input.geos
       O3COL = GET_OVERHEAD_O3(I,J)

       ! CTM ozone densities (molec/cm3)
       O3_CTM          = 0d0
       LCHEM           = State_Met%ChemGridLev(I,J)
       DO L = 1, LCHEM
          O3_CTM(L)    = State_Chm%Species(I,J,L,id_O3)
       ENDDO

       DO L = 1, State_Grid%NZ

          !-----------------------------
          ! GET CLOUD PROPERTIES BY SETTING REASONABLE VALUES FOR REL
          ! AND REI IN MICRONS AND CALCULATING LWP AND IWP FROM
          ! VISIBLE OPTICAL DEPTH  (IN G/M2)
          !-----------------------------

          IF (ICLD.NE.0) THEN
             ! LIQUID
             CLIQWP(I,J,L) = 0.667*State_Met%TAUCLW(I,J,L)*RHOLIQ*REL_DEF
             RELIQ(I,J,L)  = REL_DEF
             ! ICE
             CICEWP(I,J,L) = 0.667*State_Met%TAUCLI(I,J,L)*RHOICE*REI_DEF
             REICE(I,J,L)  = REI_DEF
             !TAUCLD DERIVED IN MCICA SUB, NOT NEEDED
             CLDFR(I,J,L)  = State_Met%CLDF(I,J,L)
          ENDIF !CLOUDS

          IF ( State_Met%InTroposphere(I,J,L) ) THEN
             !-----------------------------
             ! WE ARE IN THE TROPOSPHERE
             !-----------------------------

             ! SET O3, CH4, N2O AND CFC PROFILES
             ! G-C CHEMISTRY IS ONLY DONE IN THE TROP
             ! THEREFORE State_Chm%Species WILL ONLY BE DEFINED IN THE TROP

             !IF O3 REQUESTED THEN SPECMASK WILL BE SET TO ZERO
             !SO THAT O3 WILL BE REMOVED RELATIVE TO THE BASELINE CASE
             !(WHEN SPECMASK DEFAULTS TO 1)
             !I.E. WE WANT TO RUN WITHOUT THE GAS IF IT HAS BEEN
             !REQUESTED SO THAT WE CAN DIFFERENCE WITH THE BASELINE RUN

             IF ( Input_Opt%LUCX ) THEN

                !--------------------------------------------------------
                !          %%%%%%% UCX-based mechanisms %%%%%%%
                !--------------------------------------------------------

                IF (SPECMASK(NASPECRAD+1).EQ.1) THEN
                   O3VMR(I,J,L)  = State_Chm%Species(I,J,L,id_O3) * AIRMW / &
                                   State_Chm%SpcData(id_O3)%Info%emMW_g

                ENDIF

                IF (SPECMASK(NASPECRAD+2).EQ.1) THEN
                   CH4VMR(I,J,L) = State_Chm%Species(I,J,L,id_CH4) * AIRMW / &
                                   State_Chm%SpcData(id_CH4)%Info%emMW_g

                ENDIF

                N2OVMR(I,J,L) = State_Chm%Species(I,J,L,id_N2O) * AIRMW / &
                                State_Chm%SpcData(id_N2O)%Info%emMW_g

                CFC11VMR(I,J,L) = State_Chm%Species(I,J,L,id_CFC11) * AIRMW / &
                                  State_Chm%SpcData(id_CFC11)%Info%emMW_g

                CFC12VMR(I,J,L) = State_Chm%Species(I,J,L,id_CFC12) * AIRMW / &
                                  State_Chm%SpcData(id_CFC12)%Info%emMW_g

                CCL4VMR(I,J,L)  = State_Chm%Species(I,J,L,id_CCL4) * AIRMW / &
                                  State_Chm%SpcData(id_CCL4)%Info%emMW_g

                CFC22VMR(I,J,L) = State_Chm%Species(I,J,L,id_HCFC22) * AIRMW / &
                                  State_Chm%SpcData(id_HCFC22)%Info%emMW_g

             ELSE

                !--------------------------------------------------------
                !        %%%%%%% Tropchem only mechanisms %%%%%%%
                !--------------------------------------------------------

                IF (SPECMASK(NASPECRAD+1).EQ.1) THEN
                   O3VMR(I,J,L)  = State_Chm%Species(I,J,L,id_O3) * AIRMW / &
                                   State_Chm%SpcData(id_O3)%Info%emMW_g
                ENDIF
                IF (SPECMASK(NASPECRAD+2).EQ.1) THEN
                   CH4VMR(I,J,L) = State_Chm%Species(I,J,L,id_CH4) * AIRMW / &
                                   State_Chm%SpcData(id_CH4)%Info%emMW_g
                ENDIF
                N2OVMR(I,J,L) = N2OCLIM(I,J,L)/1E9

                !CFC CLIMATOLOGY FROM UARS AND MIPAS
                CFC11VMR(I,J,L) = CFC11CLIM(I,J,L)/1E9
                CFC12VMR(I,J,L) = CFC12CLIM(I,J,L)/1E9
                CCL4VMR(I,J,L)  = CCL4CLIM(I,J,L)/1E9
                CFC22VMR(I,J,L) = CFC22CLIM(I,J,L)/1E9
             ENDIF

          ELSE
             !-----------------------------
             ! WE ARE IN THE STRATOSPHERE
             !-----------------------------

             IF ( Input_Opt%LUCX ) THEN

                !--------------------------------------------------------
                !          %%%%%%% UCX-based mechanisms %%%%%%%
                !--------------------------------------------------------

                !N.B. STRAT CH4 NOT CURRENTLY INCLUDED IN THE DRE OF CH4
                !N.B. STRAT O3  NOT CURRENTLY INCLUDED IN THE DRE OF O3
                O3VMR(I,J,L)  = State_Chm%Species(I,J,L,id_O3) * AIRMW / &
                                State_Chm%SpcData(id_O3)%Info%emMW_g

                CH4VMR(I,J,L) = State_Chm%Species(I,J,L,id_CH4) * AIRMW / &
                                State_Chm%SpcData(id_CH4)%Info%emMW_g

                N2OVMR(I,J,L) = State_Chm%Species(I,J,L,id_N2O) * AIRMW / &
                                State_Chm%SpcData(id_N2O)%Info%emMW_g

                CFC11VMR(I,J,L) =State_Chm%Species(I,J,L,id_CFC11) * AIRMW / &
                                 State_Chm%SpcData(id_CFC11)%Info%emMW_g

                CFC12VMR(I,J,L) =State_Chm%Species(I,J,L,id_CFC12) * AIRMW / &
                                 State_Chm%SpcData(id_CFC12)%Info%emMW_g

                CCL4VMR(I,J,L)  =State_Chm%Species(I,J,L,id_CCL4) * AIRMW / &
                                 State_Chm%SpcData(id_CCL4)%Info%emMW_g

                CFC22VMR(I,J,L) =State_Chm%Species(I,J,L,id_HCFC22) * AIRMW/ &
                                 State_Chm%SpcData(id_HCFC22)%Info%emMW_g


                ! TEST IMPACT OF STRAT CHEM
                !O3VMR(I,J,L)  = 0.0d0
                !CH4VMR(I,J,L) = 0.0d0
                !N2OVMR(I,J,L) = 0.0d0
                !CFC11VMR(I,J,L) = 0.0d0
                !CFC12VMR(I,J,L) = 0.0d0
                !CCL4VMR(I,J,L)  = 0.0d0
                !CFC22VMR(I,J,L) = 0.0d0

             ELSE

                !--------------------------------------------------------
                !        %%%%%%% Tropchem only mechanisms %%%%%%%
                !--------------------------------------------------------

                !N.B. STRAT CH4 NOT CURRENTLY INCLUDED IN THE DRE OF CH4
                !N.B. STRAT O3  NOT CURRENTLY INCLUDED IN THE DRE OF O3

                !! DENSITY OF AIR IN G/CM2
                AIR_TMP = State_Met%AIRDEN(I,J,L)* &
                          State_Met%BXHEIGHT(I,J,L)*1.0E-1
                !! DENSITY OF AIR IN MOLEC/CM2
                AIR_TMP = AVO*AIR_TMP/AIRMW

                CALL SET_PROF_O3 (YLAT,THISMONTH,THISDAY,      &
                                  T_CTM,  P_CTM,               &
                                  O3_CTM, O3COL, T_CLIM,       &
                                  O3_CLIM,  Z_CLIM,  AIR_CLIM, &
                                  Input_Opt, State_Grid )
                O3VMR(I,J,L) = O3_CLIM(L)/AIR_TMP

                !GET SCALINGS IF THIS IS THE FIRST LEVEL IN THE STRAT
                IF (FLG_FIRST_STRAT(I,J).EQ.0) THEN
                   FLG_FIRST_STRAT(I,J) = 1
                   CH4SCL(I,J) = State_Chm%Species(I,J,L,id_CH4) * AIRMW / &
                                 State_Chm%SpcData(id_CH4)%Info%emMW_g /   &
                                 (CH4CLIM(I,J,L)/1E9)
                ENDIF

                !TES PROFILES INTERPOLATED TO GC GRID WHEN SAVED
                !SO WE JUST NEED TO SCALE TO CURRENT CONC AT TOP OF TROP

                CH4VMR(I,J,L) = CH4SCL(I,J)*CH4CLIM(I,J,L)/1E9
                N2OVMR(I,J,L) = N2OCLIM(I,J,L)/1E9

                !CFC CLIMATOLOGY FROM UARS AND MIPAS
                CFC11VMR(I,J,L) = CFC11CLIM(I,J,L)/1E9
                CFC12VMR(I,J,L) = CFC12CLIM(I,J,L)/1E9
                CCL4VMR(I,J,L)  = CCL4CLIM(I,J,L)/1E9
                CFC22VMR(I,J,L) = CFC22CLIM(I,J,L)/1E9
             ENDIF
          ENDIF
       ENDDO
    ENDDO
    ENDDO
    !$OMP END PARALLEL DO

    IF (DOAERAD) THEN
       DO IB = 1,NBNDS
          !RRTMG WAVEBANDS START AFTER WVAA0 STANDARD WAVELNGTHS IN GC ARRAYS
          !BASED ON LUT ORDER. JUST APPLY OFFSET
          IBX=IB+NWVAA0
          IB_SW = IB-NBNDLW
          DO IS = 1,NASPECRAD
             !THE AEROSOL SPECIES WE ARE CURRENTLY CALCULATING FOR WILL BE
             !SET TO THE LSPECRADMENU VALUE FOR THAT SPECIES.
             !THIS MEANS THAT RRTMG REQUIRES *ALL OTHER* SPECIES SO THAT THE
             !FLUX IN ABSENCE OF THE SPECIES CAN BE CALCULATED (THE
             !DIFFERENCE OF THIS WITH THE BASELINE GIVES THE FLUX CHANGE FOR
             !THAT SPECIES).
             !
             !THEREFORE WE COMPILE TWO SETS OF AEROSOL PROPERTIES:
             !(1) ALL BUT THE CURRENT SPECIES TO SEND TO RRTMG
             !(2) THE CURRENT SPECIES FOR OUTPUT TO THE RT DIAGNOSTICS
             ! ALSO, WE MUST MERGE AEROSOL PROPERTIES FOR THE SPECIES TO BE OUTPUT
             ! (I.E. COMBINE HYDROPHILIC/PHOBIC AND MULTIPLE SIZES)

             !$OMP PARALLEL DO        &
             !$OMP DEFAULT( SHARED )  &
             !$OMP PRIVATE( I, J, L ) &
             !$OMP SCHEDULE( DYNAMIC )
             DO L = 1, State_Grid%NZ
             DO J = 1, State_Grid%NY
             DO I = 1, State_Grid%NX

                !if UCX on, we need to go above the tropopause to get
                !the strat AOD, but only for IS=8 and IS=9
                IF ( State_Met%InTroposphere(I,J,L) .OR. &
                     (Input_Opt%LUCX .and. ((IS.EQ.8).OR.(IS.EQ.9)))) THEN

                   !MAKE SURE WE HAVE SENSIBLE DATA
                   !DONT WASTE TIME IF VIRTUALLY NO AEROSOL
                   IF (RTODAER(I,J,L,IBX,IS).GT.1e-10) THEN
                      IF (IB.LE.16) THEN !LW
                         IF (SPECMASK(IS).EQ.1) THEN
                            TAUAER_LW(I,J,L,IB) = TAUAER_LW(I,J,L,IB) + &
                                 RTODAER(I,J,L,IBX,IS)
                         ENDIF
                      ELSE !SW
                         !IF SPECMASK(IS)=1 THEN WE AGGREGATE THAT SPECIES FOR RRTMG
                         !IF SPECMASK(IS)>1 THEN WE SAVE THAT SPECIES FOR DIAG OUTPUT
                         IF (SPECMASK(IS).EQ.1) THEN
                            TAUAER_SW(I,J,L,IB_SW)=TAUAER_SW(I,J,L,IB_SW)+ &
                                 RTODAER(I,J,L,IBX,IS)
                            SSAAER(I,J,L,IB_SW) =  SSAAER(I,J,L,IB_SW) + &
                                 RTSSAER(I,J,L,IBX,IS)*RTODAER(I,J,L,IBX,IS)
                            ASMAER(I,J,L,IB_SW) = ASMAER(I,J,L,IB_SW) + &
                                 RTASYMAER(I,J,L,IBX,IS) * &
                                 RTODAER(I,J,L,IBX,IS)*RTSSAER(I,J,L,IBX,IS)
                         ENDIF
                         IF (SPECMASK(IS).GT.1) THEN
                            TAUAERDIAG(I,J,L,IB_SW)=TAUAERDIAG(I,J,L,IB_SW)+ &
                                 RTODAER(I,J,L,IBX,IS)
                            SSAAERDIAG(I,J,L,IB_SW) = SSAAERDIAG(I,J,L,IB_SW) +&
                                 RTSSAER(I,J,L,IBX,IS)*RTODAER(I,J,L,IBX,IS)
                            ASMAERDIAG(I,J,L,IB_SW) = ASMAERDIAG(I,J,L,IB_SW) +&
                                 RTASYMAER(I,J,L,IBX,IS) * &
                                 RTODAER(I,J,L,IBX,IS)*RTSSAER(I,J,L,IBX,IS)
                            !IF ((IS.EQ.9).AND.(L.GT.30).AND.(IB_SW.EQ.10).AND.
                            !   (RTODAER(I,J,L,IBX,IS).GT.0.0d0)) THEN
                            ! write(6,*) 'STS',I,J,L,IBX,IS,RTODAER(I,J,L,IBX,IS), &
                            !            RTSSAER(I,J,L,IBX,IS)
                            !ENDIF

                         ENDIF
                      ENDIF
                   ENDIF
                ENDIF
             ENDDO
             ENDDO
             ENDDO
             !$OMP END PARALLEL DO
          ENDDO !SPECIES

          !NOW AEROSOL HAVE BEEN SUMMED AND WEIGHTED BY AOD AND SSA
          !DIVIDE THROUGH BY TOTAL AOD (FOR SSA) AND AOD*SSA (FOR ASYM)
          IF (IB.GT.16) THEN !SW

             !$OMP PARALLEL DO        &
             !$OMP DEFAULT( SHARED )  &
             !$OMP PRIVATE( I, J, L ) &
             !$OMP SCHEDULE( DYNAMIC )
             DO L = 1, State_Grid%NZ
             DO J = 1, State_Grid%NY
             DO I = 1, State_Grid%NX

                !if UCX on, we need to go above the tropopause to get
                !the strat AOD, but only for IS=8 and IS=9
                IF ( State_Met%InTroposphere(I,J,L) .OR. &
                   (Input_Opt%LUCX .and. ((IS.EQ.8).OR.(IS.EQ.9)))) THEN

                   IF ((TAUAER_SW(I,J,L,IB_SW).GT.0).AND. &
                      (    SSAAER(I,J,L,IB_SW).GT.0)) THEN
                      !DIVIDE SUM(ASYM*SSA*OD) BY SUM(SSA*OD) TO GET
                      !OD*SSA WEIGHTED ASYM
                      ASMAER(I,J,L,IB_SW) = ASMAER(I,J,L,IB_SW) / &
                                            SSAAER(I,J,L,IB_SW)
                      !DIVIDE SUM(SSA*OD) BY SUM(OD) TO GET OD WEIGHTED SSA
                      SSAAER(I,J,L,IB_SW) = SSAAER(I,J,L,IB_SW) / &
                                            TAUAER_SW(I,J,L,IB_SW)
                   ENDIF
                   !AND DO THE SAME FOR THE SPECIES WE'RE INTERESTED IN
                   IF ((TAUAERDIAG(I,J,L,IB_SW).GT.0).AND. &
                      ( SSAAERDIAG(I,J,L,IB_SW).GT.0)) THEN
                      !DIVIDE SUM(ASYM*SSA*OD) BY SUM(SSA*OD) TO GET
                      !OD*SSA WEIGHTED ASYM
                      ASMAERDIAG(I,J,L,IB_SW) = ASMAERDIAG(I,J,L,IB_SW) / &
                                                SSAAERDIAG(I,J,L,IB_SW)
                      !DIVIDE SUM(SSA*OD) BY SUM(OD) TO GET OD WEIGHTED SSA
                      SSAAERDIAG(I,J,L,IB_SW) = SSAAERDIAG(I,J,L,IB_SW) / &
                                                TAUAERDIAG(I,J,L,IB_SW)
                   ENDIF
                ENDIF
             ENDDO
             ENDDO
             ENDDO
             !$OMP END PARALLEL DO

          ENDIF
       ENDDO !BAND
    ELSE

       !NO AEROSOL, SET ALL TO SAFE VALUES
       !$OMP PARALLEL DO       &
       !$OMP DEFAULT( SHARED ) &
       !$OMP PRIVATE( I, J, L, IB, IB_SW ) &
       !$OMP SCHEDULE( DYNAMIC )
       DO IB= 1, NBNDS
          IB_SW = IB-NBNDLW
          DO L = 1, State_Grid%NZ
          DO J = 1, State_Grid%NY
          DO I = 1, State_Grid%NX

             !if UCX on, we need to go above the tropopause to get
             !the strat AOD, but only for IS=8 and IS=9
             IF ( State_Met%InTroposphere(I,J,L) .OR. &
                (Input_Opt%LUCX .and. ((IS.EQ.8).OR.(IS.EQ.9)))) THEN

                IF (IB.LE.16) THEN
                   TAUAER_LW(I,J,L,IB)    = 0.0
                ELSE
                   TAUAER_SW(I,J,L,IB_SW) = 0.0D0
                   SSAAER(I,J,L,IB_SW)    = 0.99D0
                   ASMAER(I,J,L,IB_SW)    = 0.2D0
                   TAUAERDIAG(I,J,L,IB_SW) = 0.0D0
                   SSAAERDIAG(I,J,L,IB_SW) = 0.99D0
                   ASMAERDIAG(I,J,L,IB_SW) = 0.2D0
                ENDIF
             ENDIF
          ENDDO
          ENDDO
          ENDDO
       ENDDO
       !$OMP END PARALLEL DO

    ENDIF

    ! checking values
    DO IB= NBNDLW+1, NBNDS
       IB_SW = IB-NBNDLW
       DO L = 1, State_Grid%NZ
       DO J = 1, State_Grid%NY
       DO I = 1, State_Grid%NX
          IF (ASMAER(I,J,L,IB_SW).GT.0.999d0) THEN
             ASMAER(I,J,L,IB_SW) = 0.999d0
          ENDIF
          !IF (ASMAER(I,J,L,IB_SW).LT.0.001d0) THEN
          ! ASMAER(I,J,L,IB_SW) = 0.001
          !ENDIF
          IF ((SSAAER(I,J,L,IB_SW).LT.0.001d0).OR. &
              (SSAAER(I,J,L,IB_SW).GT.1.0d0)) THEN
             SSAAER(I,J,L,IB_SW) = 0.99
          ENDIF
          IF (TAUAER_SW(I,J,L,IB_SW).GT.1.0) THEN
             TAUAER_SW(I,J,L,IB_SW) = 1.0
          ENDIF
          IF (ASMAERDIAG(I,J,L,IB_SW).GT.0.999d0) THEN
             ASMAERDIAG(I,J,L,IB_SW) = 0.999d0
          ENDIF
          !IF (ASMAERDIAG(I,J,L,IB_SW).LT.0.001d0) THEN
          ! ASMAERDIAG(I,J,L,IB_SW) = 0.001
          !ENDIF
          IF ((SSAAERDIAG(I,J,L,IB_SW).LT.0.001d0).OR. &
              (SSAAERDIAG(I,J,L,IB_SW).GT.1.0d0)) THEN
             SSAAERDIAG(I,J,L,IB_SW) = 0.99
          ENDIF
          IF (TAUAERDIAG(I,J,L,IB_SW).GT.1.0) THEN
             TAUAERDIAG(I,J,L,IB_SW) = 1.0
          ENDIF

       ENDDO
       ENDDO
       ENDDO
    ENDDO

    DOY = GET_DAY_OF_YEAR()
    ONECOL = 1

    ! GET LEVEL VALUES
    GCAIR = 1.0E-3*GASCON/AVOGAD
    DO J=1,State_Grid%NY
    DO I=1,State_Grid%NX
       PLEV(I,J,1) = PEDGE(I,J,1) ! SET LOWEST LEVEL TO SURFACE PRESSURE
       TLEV(I,J,1) = TLAY(I,J,1)  ! SET LOWEST LEVEL TO LAYER TEMPERATURE  (KLUDGE)
       PLEV(I,J,State_Grid%NZ+1) = PCENTER(I,J,State_Grid%NZ)
       TLEV(I,J,State_Grid%NZ+1) = TLAY(I,J,State_Grid%NZ)
       DO L=2,State_Grid%NZ
          RHOA = PCENTER(I,J,L-1)/(GCAIR*TLAY(I,J,L-1))
          RHOB = PCENTER(I,J,L)/(GCAIR*TLAY(I,J,L))
          RHOSUM = RHOA+RHOB
          PLEV(I,J,L) = (RHOA*PCENTER(I,J,L-1)+RHOB*PCENTER(I,J,L))/RHOSUM
          TLEV(I,J,L) = (RHOA*TLAY(I,J,L-1)+RHOB*TLAY(I,J,L))/RHOSUM
       END DO
    END DO
    END DO

    ! FILL CO2, N2O AND O2 ARRAYS WITH REASONABLE ATMOSPHERIC VALUES
    CO2VMR(:,:,:) = 3.90E-4
    O2VMR(:,:,:)  = 0.209

    SELECT CASE (ICLD)
    ! CLOUD SETUP FOR CLEAR
    CASE (0)
       IDRV = 0
       ICLDMCL = 0
       INFLGLW = 0
       INFLGSW = 0
       TAUCMCL_LW(:,:,:,:) = 0.0
       TAUCMCL_SW(:,:,:,:) = 0.0
       ICEFLGLW = 0
       LIQFLGLW = 0
       ICEFLGSW = 0
       LIQFLGSW = 0
       !PRINT *,'CLEAR'

    !  CLOUD SETUP FOR MCICA CLOUD (ONLY OPTION NOW)
    CASE (1)
       IDRV = 0
       ICLDMCL = 2                  !MAXIMUM RANDOM OVERLAP
       INFLGLW = 2
       INFLGSW = 2
       TAUCLD_LW(:,:,:,:)  = 0.0    ! TAUCLD NOT USED
       TAUCLD_SW(:,:,:,:)  = 0.0
       TAUCMCL_LW(:,:,:,:) = 0.0   ! USED ONLY AS A CHECK
       TAUCMCL_SW(:,:,:,:) = 0.0
       SSACLD(:,:,:,:) = 0.0
       ASMCLD(:,:,:,:) = 0.0
       FSFCLD(:,:,:,:) = 0.0
       ICEFLGLW = 2       !STREAMER
       LIQFLGLW = 1       !HU AND STAMNES
       ICEFLGSW = 2       !STREAMER
       LIQFLGSW = 1       !HU AND STAMNES
       RELIQ(:,:,:) = REL_DEF
       REICE(:,:,:) = REI_DEF
    END SELECT

    ! WE ONLY NEED TO CALC CLOUDS ONCE PER RT TIMESTEP
    ! DO THIS ON BASELINE CALL IF ALL-SKY IS REQUESTED
    IF (ISPECMENU.EQ.0) THEN
       SEEDLW=ISEED+NGPTSW+1
       SEEDSW=SEEDLW+NGPTLW+1

       !$OMP PARALLEL DO       &
       !$OMP DEFAULT( SHARED ) &
       !$OMP PRIVATE( I,           J,           PCENTER0,    CLDFR0      ) &
       !$OMP PRIVATE( CLIQWP0,     CICEWP0,     REICE0,      RELIQ0      ) &
       !$OMP PRIVATE( TAUCLD_SW0,  SSACLD0,     ASMCLD0,     FSFCLD0     ) &
       !$OMP PRIVATE( CLDFMCL_LW0, CIWPMCL_LW0, CLWPMCL_LW0, REICMCL0    ) &
       !$OMP PRIVATE( RELQMCL0,    TAUCMCL_LW0, CLDFMCL_SW0, CIWPMCL_SW0 ) &
       !$OMP PRIVATE( CLWPMCL_SW0, TAUCMCL_SW0, SSACMCL0,    ASMCMCL0    ) &
       !$OMP PRIVATE( FSFCMCL0,    p_PCENTER,   p_CLDFR,     p_CICEWP    ) &
       !$OMP PRIVATE( p_CLIQWP,    p_REICE,     p_RELIQ,     p_TAUCLD_LW ) &
       !$OMP PRIVATE( p_TAUCLD_SW, p_SSACLD,    p_ASMCLD,    p_FSFCLD    ) &
       !$OMP SCHEDULE( DYNAMIC )
       DO J=1, State_Grid%NY
       DO I=1, State_Grid%NX

          ! Avoid array temporaries in the subroutine call (bmy, 6/3/15)
          ! These arrays are used by both MCICA_SUBCOL_LW and MCICA_SUBCOL_SW
          p_PCENTER = PCENTER(I,J,:)
          p_CLDFR   = CLDFR  (I,J,:)
          p_CICEWP  = CICEWP (I,J,:)
          p_CLIQWP  = CLIQWP (I,J,:)
          p_REICE   = REICE  (I,J,:)
          p_RELIQ   = RELIQ  (I,J,:)

          !-------------------------------------------------------------
          ! Long-wave radiation
          !-------------------------------------------------------------
          IF (Input_Opt%LLWRAD) THEN

             ! Avoid array temporaries in the subroutine call (bmy, 6/3/15)
             ! These arrays are only used in MCICA_SUBCOL_LW
             p_TAUCLD_LW = TAUCLD_LW(:,I,J,:)

             ! Call MCICA longwave
             CALL MCICA_SUBCOL_LW( &
                  !-------------------------------------
                  ! Inputs
                  ONECOL,        &
                  State_Grid%NZ, &
                  ICLDMCL,       &
                  SEEDLW,        &
                  IRNG,          &
                  p_PCENTER,     &
                  p_CLDFR,       &
                  p_CICEWP,      &
                  p_CLIQWP,      &
                  p_REICE,       &
                  p_RELIQ,       &
                  p_TAUCLD_LW,   &
                  !-------------------------------------
                  ! Outputs
                  CLDFMCL_LW0,  &
                  CIWPMCL_LW0,  &
                  CLWPMCL_LW0,  &
                  REICMCL0,     &
                  RELQMCL0,     &
                  TAUCMCL_LW0 )

             ! Copy back into 3-D arrays
             CLDFMCL_LW(:,I,J,:) = CLDFMCL_LW0(:,1,:)
             CIWPMCL_LW(:,I,J,:) = CIWPMCL_LW0(:,1,:)
             CLWPMCL_LW(:,I,J,:) = CLWPMCL_LW0(:,1,:)
             TAUCMCL_LW(:,I,J,:) = TAUCMCL_LW0(:,1,:)

          ENDIF

          !-------------------------------------------------------------
          ! Short-wave radiation
          !-------------------------------------------------------------
          IF (Input_Opt%LSWRAD) THEN

             ! Avoid array temporaries in the subroutine call (bmy, 6/3/15)
             ! These arrays are only used in MCICA_SUBCOL_SW
             p_TAUCLD_SW = TAUCLD_SW(:,I,J,:)
             p_SSACLD    = SSACLD   (:,I,J,:)
             p_ASMCLD    = ASMCLD   (:,I,J,:)
             p_FSFCLD    = FSFCLD   (:,I,J,:)

             ! Call MCICA shortwave
             CALL MCICA_SUBCOL_SW( &
                  !-------------------------------------
                  ! Inputs
                  ONECOL,        &
                  State_Grid%NZ, &
                  ICLDMCL,       &
                  SEEDSW,        &
                  IRNG,          &
                  p_PCENTER,     &
                  p_CLDFR,       &
                  p_CICEWP,      &
                  p_CLIQWP,      &
                  p_REICE,       &
                  p_RELIQ,       &
                  p_TAUCLD_SW,   &
                  p_SSACLD,      &
                  p_ASMCLD,      &
                  p_FSFCLD,      &
                  !-------------------------------------
                  ! Outputs
                  CLDFMCL_SW0,   &
                  CIWPMCL_SW0,   &
                  CLWPMCL_SW0,   &
                  REICMCL0,      &
                  RELQMCL0,      &
                  TAUCMCL_SW0,   &
                  SSACMCL0,      &
                  ASMCMCL0,      &
                  FSFCMCL0  )

             ! Copy back into 3-D arrays
             CLDFMCL_SW(:,I,J,:) = CLDFMCL_SW0(:,1,:)
             CIWPMCL_SW(:,I,J,:) = CIWPMCL_SW0(:,1,:)
             CLWPMCL_SW(:,I,J,:) = CLWPMCL_SW0(:,1,:)
             TAUCMCL_SW(:,I,J,:) = TAUCMCL_SW0(:,1,:)
             SSACMCL   (:,I,J,:) = SSACMCL0   (:,1,:)
             ASMCMCL   (:,I,J,:) = ASMCMCL0   (:,1,:)
             FSFCMCL   (:,I,J,:) = FSFCMCL0   (:,1,:)

          ENDIF

          ! these should be independent of LW and SW
          ! simply rearranged by the MCICA routine
          REICMCL(I,J,:) = REICMCL0(1,:)
          RELQMCL(I,J,:) = RELQMCL0(1,:)
       ENDDO
       ENDDO
       !$OMP END PARALLEL DO

    ENDIF !DO MCICA CLOUDS

    !$OMP PARALLEL DO       &
    !$OMP DEFAULT( SHARED ) &
    !$OMP PRIVATE( I,   J,       UFLX,         DFLX,         HR           ) &
    !$OMP PRIVATE( UFLXC,        DFLXC,        HRC,          DUFLX_DT     ) &
    !$OMP PRIVATE( DUFLXC_DT,    ECAER,        SWUFLX,       SWDFLX       ) &
    !$OMP PRIVATE( SWHR,         SWUFLXC,      SWDFLXC,      SWHRC        ) &
    !$OMP PRIVATE( p_PCENTER,    p_PLEV,       p_TLAY,       p_TLEV       ) &
    !$OMP PRIVATE( p_H2OVMR,     p_O3VMR,      p_CO2VMR,     p_CH4VMR     ) &
    !$OMP PRIVATE( p_N2OVMR,     p_O2VMR,      p_CFC11VMR,   p_CFC12VMR   ) &
    !$OMP PRIVATE( p_CFC22VMR,   p_CCL4VMR,    p_RTEMISS,    p_REICMCL    ) &
    !$OMP PRIVATE( p_RELQMCL,    p_CLDFMCL_LW, p_TAUCMCL_LW, p_CIWPMCL_LW ) &
    !$OMP PRIVATE( p_CLWPMCL_LW, p_TAUAER_LW,  p_CLDFMCL_SW, p_TAUCMCL_SW ) &
    !$OMP PRIVATE( p_SSACMCL,    p_ASMCMCL,    p_FSFCMCL,    p_CIWPMCL_SW ) &
    !$OMP PRIVATE( p_CLWPMCL_SW, p_TAUAER_SW,  p_SSAAER,     p_ASMAER     ) &
    !$OMP PRIVATE( p_SUNCOS                                               ) &
    !$OMP SCHEDULE( DYNAMIC )
    DO J=1, State_Grid%NY
    DO I=1, State_Grid%NX

       ! Avoid arrray temporaries in subroutines (bmy, 6/3/15)
       ! These variables are used in both RRTMG_LW and RRTMG_SW
       p_PCENTER  = PCENTER (I,J,:)
       p_PLEV     = PLEV    (I,J,:)
       p_TLAY     = TLAY    (I,J,:)
       p_TLEV     = TLEV    (I,J,:)
       p_H2OVMR   = H2OVMR  (I,J,:)
       p_O3VMR    = O3VMR   (I,J,:)
       p_CO2VMR   = CO2VMR  (I,J,:)
       p_CH4VMR   = CH4VMR  (I,J,:)
       p_N2OVMR   = N2OVMR  (I,J,:)
       p_O2VMR    = O2VMR   (I,J,:)
       p_CFC11VMR = CFC11VMR(I,J,:)
       p_CFC12VMR = CFC12VMR(I,J,:)
       p_CFC22VMR = CFC22VMR(I,J,:)
       p_CCL4VMR  = CCL4VMR (I,J,:)
       p_RTEMISS  = RTEMISS (I,J,:)
       p_REICMCL  = REICMCL (I,J,:)
       p_RELQMCL  = RELQMCL (I,J,:)

       !--------------------------------------------------------------
       ! RRTMG - Longwave radiation
       !--------------------------------------------------------------
       IF (Input_Opt%LLWRAD) THEN

          ! Avoid array temporaries in subroutine calls (bmy, 6/3/15)
          ! These arrays are only used in RRTMG_LW
          p_CLDFMCL_LW = CLDFMCL_LW(:,I,J,:  )
          p_TAUCMCL_LW = TAUCMCL_LW(:,I,J,:  )
          p_CIWPMCL_LW = CIWPMCL_LW(:,I,J,:  )
          p_CLWPMCL_LW = CLWPMCL_LW(:,I,J,:  )
          p_TAUAER_LW  = TAUAER_LW (  I,J,:,:)

          ! Call RRTMG for longwave radiation
          CALL RRTMG_LW( &
               !-------------------------------------
               ! Inputs
               ONECOL,        &
               State_Grid%NZ, &
               ICLDMCL,       &
               IDRV,          &
               p_PCENTER,     &
               p_PLEV,        &
               p_TLAY,        &
               p_TLEV,        &
               TSFC(I,J),     &
               p_H2OVMR,      &
               p_O3VMR,       &
               p_CO2VMR,      &
               p_CH4VMR,      &
               p_N2OVMR,      &
               p_O2VMR,       &
               p_CFC11VMR,    &
               p_CFC12VMR,    &
               p_CFC22VMR,    &
               p_CCL4VMR,     &
               p_RTEMISS,     &
               INFLGLW,       &
               ICEFLGLW,      &
               LIQFLGLW,      &
               p_CLDFMCL_LW,  &
               p_TAUCMCL_LW,  &
               p_CIWPMCL_LW,  &
               p_CLWPMCL_LW,  &
               p_REICMCL,     &
               p_RELQMCL,     &
               p_TAUAER_LW,   &
               !-------------------------------------
               ! Outputs
               UFLX,          &
               DFLX,          &
               HR,            &
               UFLXC,         &
               DFLXC,         &
               HRC,           &
               DUFLX_DT,      &
               DUFLXC_DT )

          ! Copy back into 3-D arrays
          LW_UFLUX (I,J,:) = UFLX (1,:)
          LW_DFLUX (I,J,:) = DFLX (1,:)
          LW_UFLUXC(I,J,:) = UFLXC(1,:)
          LW_DFLUXC(I,J,:) = DFLXC(1,:)

       ENDIF                  !LW

       !--------------------------------------------------------------
       ! RRTMG - Shortwave radiation
       !--------------------------------------------------------------
       IF (Input_Opt%LSWRAD) THEN

          !### Debug output
          !write(6,*) 'SWSHIZ',TAUAER_SW(I,J,1,5),SWUFLX(1,1),SWDFLX(1,1), &
          !           SWUFLXC(1,1),SWDFLXC(1,1)

          ! Avoid array temporaries in subroutine calls (bmy, 6/3/15)
          ! These arrays are only used in RRTMG_SW
          p_SUNCOS     = SUNCOS    (  I,J,:  )
          p_CLDFMCL_SW = CLDFMCL_SW(:,I,J,:  )
          p_TAUCMCL_SW = TAUCMCL_SW(:,I,J,:  )
          p_SSACMCL    = SSACMCL   (:,I,J,:  )
          p_ASMCMCL    = ASMCMCL   (:,I,J,:  )
          p_FSFCMCL    = FSFCMCL   (:,I,J,:  )
          p_CIWPMCL_SW = CIWPMCL_SW(:,I,J,:  )
          p_CLWPMCL_SW = CLWPMCL_SW(:,I,J,:  )
          p_TAUAER_SW  = TAUAER_SW (  I,J,:,:)
          p_SSAAER     = SSAAER    (  I,J,:,:)
          p_ASMAER     = ASMAER    (  I,J,:,:)

          ! Call RRTMG for shortwave radiation
          CALL RRTMG_SW( &
               !-------------------------------------
               ! Inputs
               ONECOL,         &
               State_Grid%NZ,  &
               ICLDMCL,        &
               p_PCENTER,      &
               p_PLEV,         &
               p_TLAY,         &
               p_TLEV,         &
               TSFC(I,J),      &
               p_H2OVMR,       &
               p_O3VMR,        &
               p_CO2VMR,       &
               p_CH4VMR,       &
               p_N2OVMR,       &
               p_O2VMR,        &
               ALBDIRVIS(I,J), &
               ALBDIFVIS(I,J), &
               ALBDIRNIR(I,J), &
               ALBDIFNIR(I,J), &
               p_SUNCOS,       &
               ADJES,          &
               DOY,            &
               SCON,           &
               INFLGSW,        &
               ICEFLGSW,       &
               LIQFLGSW,       &
               p_CLDFMCL_SW,   &
               p_TAUCMCL_SW,   &
               p_SSACMCL,      &
               p_ASMCMCL,      &
               p_FSFCMCL,      &
               p_CIWPMCL_SW,   &
               p_CLWPMCL_SW,   &
               p_REICMCL,      &
               p_RELQMCL,      &
               p_TAUAER_SW,    &
               p_SSAAER,       &
               p_ASMAER,       &
               !-------------------------------------
               ! Outputs
               ECAER,          &
               SWUFLX,         &
               SWDFLX,         &
               SWHR,           &
               SWUFLXC,        &
               SWDFLXC,        &
               SWHRC    )

          ! Copy back into 3-D arrays
          SW_UFLUX (I,J,:) = SWUFLX (1,:)
          SW_DFLUX (I,J,:) = SWDFLX (1,:)
          SW_UFLUXC(I,J,:) = SWUFLXC(1,:)
          SW_DFLUXC(I,J,:) = SWDFLXC(1,:)

       ENDIF !SW
    ENDDO !State_Grid%NX
    ENDDO !State_Grid%NY
    !$OMP END PARALLEL DO

    ! OUTPUT RADIATION VARIABLES TO DIAGNOSTIC
    ! IF CALC WITH AEROSOLS AND GASES COMPLETED
    ! USE ISPECMENU (REFERENCES THE INPUT.GEOS.RAD LIST)
    ! TO DETERMINE WHICH FLUX HAS BEEN CALCULATED
    ! OUTPUT DIAGNOSTIC INDEX IS ISPECMENU+1 (ISPECMENU=0 FOR BASELINE)
    OUTIDX = ISPECMENU + 1

    !THE NUMBER OF ND72 OUTPUTS PER FIELD
    IF ( Input_Opt%LUCX ) THEN
       NAD72 = Input_Opt%NSPECRADMENU + 1
    ELSE
       NAD72 = Input_Opt%NSPECRADMENU
    ENDIF

    !FIRST CHECK IF WE HAVE ALREADY OUTPUT AEROSOL DIAGNOSTICS
    !(I.E. IF BOTH ALL-SKY AND CLEAR-SKY ARE SWITCHED ON)
    IF ((Input_Opt%LSKYRAD(1)).AND.(Input_Opt%LSKYRAD(2))) THEN
       !WE ONLY NEED TO OUTPUT DURING ONE OF THESE, SO DONT WHEN ICLD=0
       IF (ICLD.EQ.0) THEN
          LOUTPUTAERO=.FALSE.
       ENDIF
    ENDIF

    !$OMP PARALLEL DO       &
    !$OMP DEFAULT( SHARED ) &
    !$OMP PRIVATE( I, J, LL, W )            &
    !$OMP PRIVATE( AODTMP, SSATMP, ASYMTMP) &
    !$OMP PRIVATE( AODOUT, SSAOUT, ASYMOUT) &
    !$OMP SCHEDULE( DYNAMIC )
    DO J=1,State_Grid%NY
    DO I=1,State_Grid%NX
#ifdef BPCH_DIAG
       !================================================================
       ! %%%%% ND72 (bpch) DIAGNOSTIC %%%%%
       !
       ! Save clear-sky and all-sky fluxes from RRTMG [W/m2]
       !================================================================

       IF (ICLD.GT.0) THEN
          !-------------------------------------------------------
          !ALL-SKY (WE GET CLEAR-SKY WITH THIS TOO)
          !N.B. UPWELLING SHOULD BE NEGATIVE AS DOWN IS +VE
          !-------------------------------------------------------

          ! All-sky SW flux @ TOA [W/m2]
          AD72(I,J,OUTIDX) = AD72(I,J,OUTIDX) - &
                             SNGL(SW_UFLUX(I,J,State_Grid%NZ+1))

          ! All-sky SW flux @ surface [W/m2]
          AD72(I,J,OUTIDX+NAD72) = AD72(I,J,OUTIDX+NAD72) + &
                                   SNGL(SW_DFLUX(I,J,1))

          ! All-sky LW flux @ TOA [W/m2]
          AD72(I,J,OUTIDX+2*NAD72) = AD72(I,J,OUTIDX+2*NAD72) - &
                                     SNGL(LW_UFLUX(I,J,State_Grid%NZ+1))

          ! All-sky LW flux @ surface [W/m2]
          AD72(I,J,OUTIDX+3*NAD72) = AD72(I,J,OUTIDX+3*NAD72) + &
                                     SNGL(LW_DFLUX(I,J,1))

          ! Clear-sky SW flux @ TOA [W/m2]
          AD72(I,J,OUTIDX+4*NAD72) = AD72(I,J,OUTIDX+4*NAD72) - &
                                     SNGL(SW_UFLUXC(I,J,State_Grid%NZ+1))

          ! Clear-sky SW flux @ surface [W/m2]
          AD72(I,J,OUTIDX+5*NAD72) = AD72(I,J,OUTIDX+5*NAD72) + &
                                     SNGL(SW_DFLUXC(I,J,1))

          ! Clear-sky LW flux @ TOA [W/m2]
          AD72(I,J,OUTIDX+6*NAD72) = AD72(I,J,OUTIDX+6*NAD72) - &
                                     SNGL(LW_UFLUXC(I,J,State_Grid%NZ+1))

          ! Clear-sky LW flux @ surface [W/m2]
          AD72(I,J,OUTIDX+7*NAD72) = AD72(I,J,OUTIDX+7*NAD72) + &
                                     SNGL(LW_DFLUXC(I,J,1))
       ELSE
          !-------------------------------------------------------
          ! CLEAR-SKY (RUNNING WITH CLOUDS OFF)
          !-------------------------------------------------------

          ! Clear-sky SW flux @ TOA [w/m2]
          AD72(I,J,OUTIDX+4*NAD72) = AD72(I,J,OUTIDX+4*NAD72) - &
                                     SNGL(SW_UFLUX(I,J,State_Grid%NZ+1))

          ! Clear-sky SW flux @ surface [W/m2]
          AD72(I,J,OUTIDX+5*NAD72) = AD72(I,J,OUTIDX+5*NAD72) + &
                                     SNGL(SW_DFLUX(I,J,1))

          ! Clear-sky LW flux @ TOA [W/m2]
          AD72(I,J,OUTIDX+6*NAD72) = AD72(I,J,OUTIDX+6*NAD72) - &
                                     SNGL(LW_UFLUX(I,J,State_Grid%NZ+1))

          ! Clear-sky LW flux @ surface [W/m2]
          AD72(I,J,OUTIDX+7*NAD72) = AD72(I,J,OUTIDX+7*NAD72) + &
                                     SNGL(LW_DFLUX(I,J,1))
       ENDIF
#endif

       !================================================================
       ! %%%%% HISTORY (aka netCDF diagnostics) %%%%%
       !
       ! Save clear-sky and all-sky fluxes from RRTMG [W/m2]
       !================================================================
       IF ( ICLD > 0 ) THEN

          !-------------------------------------------------------
          !ALL-SKY (WE GET CLEAR-SKY WITH THIS TOO)
          !N.B. UPWELLING SHOULD BE NEGATIVE AS DOWN IS +VE
          !-------------------------------------------------------

          ! All-sky SW flux @ TOA [W/m2]
          IF ( State_Diag%Archive_RadAllSkySWTOA ) THEN
             State_Diag%RadAllSkySWTOA(I,J,iNcDiag) = &
                  -SW_UFLUX(I,J,State_Grid%NZ+1)
          ENDIF

          ! All-sky SW flux @ surface [W/m2]
          IF ( State_Diag%Archive_RadAllSkySWSurf ) THEN
             State_Diag%RadAllSkySWSurf(I,J,iNcDiag) = &
                  SW_DFLUX(I,J,1)
          ENDIF

          ! All-sky LW flux @ TOA [W/m2]
          IF ( State_Diag%Archive_RadAllSkyLWTOA ) THEN
             State_Diag%RadAllSkyLWTOA(I,J,iNcDiag) = &
                  -LW_UFLUX(I,J,State_Grid%NZ+1)
          ENDIF

          ! All-sky LW flux @ surface [W/m2]
          IF ( State_Diag%Archive_RadAllSkyLWSurf ) THEN
             State_Diag%RadAllSkyLWSurf(I,J,iNcDiag) = &
                  LW_DFLUX(I,J,1)
          ENDIF

          ! Clear-sky SW flux @ TOA [W/m2]
          IF ( State_Diag%Archive_RadClrSkySWTOA ) THEN
             State_Diag%RadClrSkySWTOA(I,J,iNcDiag) = &
                  -SW_UFLUXC(I,J,State_Grid%NZ+1)
          ENDIF

          ! Clear-sky SW flux @ surface [W/m2]
          IF ( State_Diag%Archive_RadClrSkySWSurf ) THEN
             State_Diag%RadClrSkySWSurf(I,J,iNcDiag) = &
                  SW_DFLUXC(I,J,1)
          ENDIF

          ! Clear-sky LW flux @ TOA [W/m2]
          IF ( State_Diag%Archive_RadClrSkyLWTOA ) THEN
             State_Diag%RadClrSkyLWTOA(I,J,iNcDiag) = &
                  -LW_UFLUXC(I,J,State_Grid%NZ+1)
          ENDIF

          ! Clear-sky LW flux @ surface [W/m2]
          IF ( State_Diag%Archive_RadClrSkyLWSurf ) THEN
             State_Diag%RadClrSkyLWSurf(I,J,iNcDiag) = &
                  LW_DFLUXC(I,J,1)
          ENDIF

       ELSE

          !-------------------------------------------------------
          ! CLEAR-SKY (RUNNING WITH CLOUDS OFF)
          !-------------------------------------------------------

          ! Clear-sky SW flux @ TOA [W/m2]
          IF ( State_Diag%Archive_RadClrSkySWTOA ) THEN
             State_Diag%RadClrSkySWTOA(I,J,iNcDiag) = &
                  -SW_UFLUXC(I,J,State_Grid%NZ+1)
          ENDIF

          ! Clear-sky SW flux @ surface [W/m2]
          IF ( State_Diag%Archive_RadClrSkySWSurf ) THEN
             State_Diag%RadClrSkySWSurf(I,J,iNcDiag) = &
                  SW_DFLUXC(I,J,1)
          ENDIF

          ! Clear-sky LW flux @ TOA [W/m2]
          IF ( State_Diag%Archive_RadClrSkyLWTOA ) THEN
             State_Diag%RadClrSkyLWTOA(I,J,iNcDiag) = &
                  -LW_UFLUXC(I,J,State_Grid%NZ+1)
          ENDIF

          ! Clear-sky LW flux @ surface [W/m2]
          IF ( State_Diag%Archive_RadClrSkyLWSurf ) THEN
             State_Diag%RadClrSkyLWSurf(I,J,iNcDiag) = &
                  LW_DFLUXC(I,J,1)
          ENDIF

       ENDIF

       !OUTPUT OPTICS FOR EACH AEROSOL...
       !CHECK THAT WE HAVE SOME AEROSOL TO OUTPUT
       !SKIP OUTIDX=1,2,3 (BASELINE, OZONE AND CH4)
       IF ((OUTIDX.GE.4).AND.(LOUTPUTAERO)) THEN
          !INTERPOLATE TO THE REQUESTED WAVELENGTH
          DO W=1,Input_Opt%NWVSELECT
             AODTMP  = 0.0D0
             SSATMP  = 0.0D0
             ASYMTMP = 0.0D0
             AODOUT  = 0.0D0
             SSAOUT  = 0.0D0
             ASYMOUT = 0.0D0
             DO LL=1,State_Grid%NZ
                !CHECK AOD IS NON-ZERO BEFORE LOG...
                IF((TAUAERDIAG(I,J,LL,IRTWVSELECT(2,W)).GT.0).AND. &
                   (TAUAERDIAG(I,J,LL,IRTWVSELECT(1,W)).GT.0)) THEN
                   AODTMP=SNGL(TAUAERDIAG(I,J,LL,IRTWVSELECT(2,W))* &
                        ACOEF_RTWV(W)**(BCOEF_RTWV(W)*              &
                        LOG(TAUAERDIAG(I,J,LL,IRTWVSELECT(1,W))/    &
                        TAUAERDIAG(I,J,LL,IRTWVSELECT(2,W)))))
                   SSATMP=SNGL( CCOEF_RTWV(W)*                      &
                        SSAAERDIAG(I,J,LL,IRTWVSELECT(2,W))+        &
                        (1.0D0-CCOEF_RTWV(W))*                      &
                        SSAAERDIAG(I,J,LL,IRTWVSELECT(1,W)))
                   ASYMTMP=SNGL( CCOEF_RTWV(W)*                     &
                        ASMAERDIAG(I,J,LL,IRTWVSELECT(2,W))+        &
                        (1.0D0-CCOEF_RTWV(W))*                      &
                        ASMAERDIAG(I,J,LL,IRTWVSELECT(1,W)))
                   AODOUT=AODOUT+AODTMP
                   SSAOUT=SSAOUT+SSATMP*AODTMP
                   ASYMOUT=ASYMOUT+ASYMTMP*SSATMP*AODTMP
                ENDIF
             ENDDO !State_Grid%NZ
             !WE ARE SAVING COLUMN AVERAGED VALUES FOR EACH SPECIES
             !DIVIDE THROUGH BY AOD*SSA (AOD-SSA WEIGHTING ACCOUNTS FOR
             !GRIDBOXES)
             ASYMOUT=ASYMOUT/SSAOUT
             !DIVIDE THROUGH BY AOD
             SSAOUT=SSAOUT/AODOUT
             !offsetting output depending on wavelength
#ifdef BPCH_DIAG
             AD72(I,J,OUTIDX+(8+3*(W-1))*NAD72) = &
                  AD72(I,J,OUTIDX+(8+3*(W-1))*NAD72) + AODOUT
             AD72(I,J,OUTIDX+(9+3*(W-1))*NAD72) = &
                  AD72(I,J,OUTIDX+(9+3*(W-1))*NAD72) + SSAOUT
             AD72(I,J,OUTIDX+(10+3*(W-1))*NAD72)= &
                  AD72(I,J,OUTIDX+(10+3*(W-1))*NAD72) + ASYMOUT
#endif
          ENDDO !NWVSELECT
       ENDIF
    ENDDO
    ENDDO
    !$OMP END PARALLEL DO

    ! Convert species units back to original unit
    CALL Convert_Spc_Units( Input_Opt, State_Chm, State_Grid, State_Met, &
                            OrigUnit,  RC )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL GC_Error('Unit conversion error', RC, 'DO_RRTMG_RAD_TRANSFER')
       RETURN
    ENDIF

  END SUBROUTINE DO_RRTMG_RAD_TRANSFER
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: set_specmask
!
! !DESCRIPTION: Subroutine SET\_SPECMASK converts the species switches in the
!  input.mod radiation section into the list of species that should be passed
!  through to RRTMG. This must be done in a subtractive way, e.g. If we require
!  the DRE of sulfate then the baseline will contain all species and the
!  sulfate run will contain everything but sulfate, this way the contribution
!  of sulfate can be inferred. Therefore, all species are initially set to 1
!  and their inclusion results in SPECMASK for the particular species being
!  set to zero. (dar 10/2013)
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Set_SpecMask( iSpecRadMenu )
!
! !USES:
!
    USE CMN_FJX_MOD, ONLY : SPECMASK, NSPECRAD, NASPECRAD
!
! !INPUT PARAMETERS:
!
    INTEGER, INTENT(IN) :: iSpecRadMenu  ! Index of RRTMG flux output
!
! !REVISION HISTORY:
!  18 Jun 2013 - D.A. Ridley - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER :: N0,N,I,II,NXTRA

    !=================================================================
    ! SET_SPECMASK begins here!
    !=================================================================

    !ISPECRADMENU IS THE INDEX OF THE SPECIES WITHIN THE INPUT MENU
    !THE INDEX OF SPECMASK INDICATES THE POSITION OF THE SPECIES IN
    !THE RT OPTICS ARRAY FOR THE AEROSOL
    !E.G. SO4 = 1 I.E. RTODAER(*,*,*,*,1)
    !
    !===FUNCTIONALITY FOR ADDING NEW SPECIES===
    !EXTRA SPECIES ARE ADDED AFTER NAER (BEFORE NDUST AND GASES)
    !SO WE NEED TO BUMP ALL THE SPECIES AFTER THAT BY NXTRA
    !WHERE NXTRA=NUMBER OF NEW SPECIES ADDED ABOVE THE STANDARD CODE
    !E.G. IF UCX=YES THEN NASPECRAD=18 AND STS AND NAT ARE INCLUDED
    !IN RTODAER INDEX 8 AND 9, BEFORE DUST
    NXTRA=NSPECRAD-16

    !CONVERT THE CURRENT SPECIES SELECTION FROM THE INPUT MENU INTO
    !THE REQUIRED SPECIES TO BE INCLUDED IN THE RRTMG CALCULATION
    SPECMASK(:)=1

    !IF ISPECRADMENU IS ZERO, WE JUST WANTED BASELINE, I.E. SPECMASK(:)=1
    IF ( ISpecRadMenu .GT. 0 ) THEN

       SELECT CASE( ISpecRadMenu )

       ! O3 = Ozone
       CASE( 1 )
          SPECMASK(15+NXTRA)=0

       ! ME = Methane
       CASE( 2 )
          SPECMASK(16+NXTRA)=0

       ! SU = Sulfate
       CASE( 3 )
          SPECMASK(1)=3

       ! NI = Nitrate
       CASE( 4 )
          SPECMASK(2)=4

       ! AM = Ammonium
       CASE( 5 )
          SPECMASK(3)=5

       ! BC = Black carbon (Hydrophilic+phobic)
       CASE( 6 )
          SPECMASK(4)=6

       ! OA = Organic aerosol (!Hydrophilic+phobic)
       CASE( 7 )
          SPECMASK(5)=7

       ! SS = Sea salt
       CASE( 8 )
          SPECMASK(6)=8
          SPECMASK(7)=8

       ! DU = Mineral dust
       CASE( 9 )
          SPECMASK(8+NXTRA)=9
          SPECMASK(9+NXTRA)=9
          SPECMASK(10+NXTRA)=9
          SPECMASK(11+NXTRA)=9
          SPECMASK(12+NXTRA)=9
          SPECMASK(13+NXTRA)=9
          SPECMASK(14+NXTRA)=9

       ! PM = All particulate matter
       ! add all aerosols but not gases here
       CASE( 10 )
          DO II = 1, NASPECRAD
             SPECMASK(II)=10
          ENDDO

       ! ST = STRAT AEROSOL (UCX only)
       CASE( 11 )

          !LSA
          SPECMASK(8) = 11

          !NAT
          SPECMASK(9) = 11

       END SELECT
    ENDIF

  END SUBROUTINE Set_SpecMask
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_rrtmg_rad_transfer
!
! !DESCRIPTION: Initializes all RRTMG module varaiables.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Init_RRTMG_Rad_Transfer( Input_Opt, State_Diag, State_Grid, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE Input_Opt_Mod,  ONLY : OptInput
    USE State_Diag_Mod, ONLY : DgnState
    USE State_Grid_Mod, ONLY : GrdState
!
! !INPUT PARAMETERS:
!
    TYPE(GrdState), INTENT(IN)    :: State_Grid  ! Grid State object
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(INOUT) :: Input_Opt   ! Input Options object
    TYPE(DgnState), INTENT(INOUT) :: State_Diag  ! Diagnostics State object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure
!
! !REMARKS:
!  The index fields Input_Opt%RadFluxCt, Input_Opt%RadFluxName, and
!  Input_Opt%RadFluxInd are populated from information obtained in
!  Headers/diaglist_mod.F90.  But the input.geos file is read before
!  the diaglist is constructed.  Therefore, we have to delay population
!  of these fields until after the call to Init_DiagList.
!
! !REVISION HISTORY:
!  09 Nov 2018 - R. Yantosca - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER :: N

    ! Strings
    CHARACTER(LEN=255) :: ErrMsg, ThisLoc, arrayId

    !=================================================================
    ! Init_RRTMG_Inputs begins here
    !=================================================================

    ! Initialize
    RC      = GC_SUCCESS
    ErrMsg  = ''
    ThisLoc = ' -> at Init_RRTMG_Rad_Transfer (in rrtmg_rad_transfer_mod.F90)'

    !=================================================================
    ! For backwards compatibility with existing RRTMG code, we need
    ! to populate the Input_Opt%LSpecRadMenu based on the flux
    ! outputs requested in the HISTORY.rc file.  Loop over all
    ! possible types of RRTMG flux outputs (excluding baseline,
    ! which is type 0).
    !
    ! Save the name of each flux output in Input_Opt%RadFluxName
    ! and its expected index value in Input_Opt%RadFluxInd.
    ! Expected index values for flux ouptuts are:
    !
    ! Optional outputs (requested via HISTORY.rc)
    !   1=O3  2=ME  3=SU   4=NI  5=AM  6=BC
    !   7=OA  8=SS  9=DU  10=PM  11=ST (UCX only)
    !
    ! This is a bit convoluted but we need to do this in order to
    ! keep track of the slot of the netCDF diagnostic arrays in
    ! State_Diag in which to archive the various flux outputs.
    ! Also this allows us to keep backwards compatibility with the
    ! existing code to the greatest extent.
    !
    ! NOTE: We can get rid of Input_Opt%LSPECRADMENU once all of
    ! the bpch code is removed from GEOS-Chem.  This array is still
    ! used in diag3.F90, so we need to keep it for the time being.
    ! (bmy, 11/9/18)
    !=================================================================

    ! Loop over all of the flux outputs requested in HISTORY.rc
    DO N = 1, State_Diag%nRadFlux

       SELECT CASE( State_Diag%RadFluxName(N) )
       CASE( 'O3' )
          Input_Opt%LSpecRadMenu(1)  = 1
       CASE( 'ME' )
          Input_Opt%LSpecRadMenu(2)  = 1
       CASE( 'SU' )
          Input_Opt%LSpecRadMenu(3)  = 1
       CASE( 'NI' )
          Input_Opt%LSpecRadMenu(4)  = 1
       CASE( 'AM' )
          Input_Opt%LSpecRadMenu(5)  = 1
       CASE( 'BC' )
          Input_Opt%LSpecRadMenu(6)  = 1
       CASE( 'OA' )
          Input_Opt%LSpecRadMenu(7)  = 1
       CASE( 'SS' )
          Input_Opt%LSpecRadMenu(8)  = 1
       CASE( 'DU' )
          Input_Opt%LSpecRadMenu(9)  = 1
       CASE( 'PM' )
          Input_Opt%LSpecRadMenu(10) = 1
       CASE( 'ST' )
          Input_Opt%LSpecRadMenu(11) = 1
       CASE DEFAULT
          ! Nothing
       END SELECT
    ENDDO

    !=================================================================
    ! Allocate arrays
    !=================================================================
    CALL Init_Surface_Rad ( State_Grid )
    CALL Init_MCICA_Clouds( State_Grid )

  END SUBROUTINE Init_RRTMG_Rad_Transfer
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_surface_rad
!
! !DESCRIPTION: Subroutine INIT\_SURFACE\_RAD initializes all allocatable
!  module arrays.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Init_Surface_Rad( State_Grid )
!
! !USES:
!
    USE CMN_FJX_MOD
    USE ERROR_MOD,      ONLY : ALLOC_ERR
    USE State_Grid_Mod, ONLY : GrdState
!
! !INPUT PARAMETERS:
!
    TYPE(GrdState), INTENT(IN) :: State_Grid  ! Grid State objecy
!
! !REVISION HISTORY:
!  18 Jun 2013 - D.A. Ridley - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER :: AS

    !=================================================================
    ! INIT_SURFACE_RAD begins here!
    !=================================================================

    ALLOCATE( LW_UFLUX(State_Grid%NX,State_Grid%NY,State_Grid%NZ+1), &
              STAT=AS )
    IF ( AS /= 0 ) CALL ALLOC_ERR( 'LW_UFLUX' )
    LW_UFLUX = 0D0

    ALLOCATE( LW_DFLUX(State_Grid%NX,State_Grid%NY,State_Grid%NZ+1), &
              STAT=AS )
    IF ( AS /= 0 ) CALL ALLOC_ERR( 'LW_DFLUX' )
    LW_DFLUX = 0D0

    ALLOCATE( SW_UFLUX(State_Grid%NX,State_Grid%NY,State_Grid%NZ+1), &
              STAT=AS )
    IF ( AS /= 0 ) CALL ALLOC_ERR( 'SW_UFLUX' )
    SW_UFLUX = 0D0

    ALLOCATE( SW_DFLUX(State_Grid%NX,State_Grid%NY,State_Grid%NZ+1), &
              STAT=AS )
    IF ( AS /= 0 ) CALL ALLOC_ERR( 'SW_DFLUX' )
    SW_DFLUX = 0D0

    ALLOCATE( LW_UFLUXC(State_Grid%NX,State_Grid%NY,State_Grid%NZ+1), &
              STAT=AS )
    IF ( AS /= 0 ) CALL ALLOC_ERR( 'LW_UFLUXC' )
    LW_UFLUXC = 0D0

    ALLOCATE( LW_DFLUXC(State_Grid%NX,State_Grid%NY,State_Grid%NZ+1), &
              STAT=AS )
    IF ( AS /= 0 ) CALL ALLOC_ERR( 'LW_DFLUXC' )
    LW_DFLUXC = 0D0

    ALLOCATE( SW_UFLUXC(State_Grid%NX,State_Grid%NY,State_Grid%NZ+1), &
              STAT=AS )
    IF ( AS /= 0 ) CALL ALLOC_ERR( 'SW_UFLUXC' )
    SW_UFLUXC = 0D0

    ALLOCATE( SW_DFLUXC(State_Grid%NX,State_Grid%NY,State_Grid%NZ+1), &
              STAT=AS )
    IF ( AS /= 0 ) CALL ALLOC_ERR( 'SW_DFLUXC' )
    SW_DFLUXC = 0D0

  END SUBROUTINE INIT_SURFACE_RAD
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: AttachPointersFromHemco
!
! !DESCRIPTION: Subroutine AttachPointersFromHemco attaches pointers for
!  various RRTMG input data that is read in via HEMCO to module variables.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE AttachPointersFromHemco( RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE HCO_Interface_Mod, ONLY : HcoState
    USE HCO_EmisList_Mod,  ONLY : HCO_GetPtr
!
! !OUTPUT PARAMETERS:
!
    INTEGER, INTENT(OUT) :: RC          ! Success or failure
!
! !REVISION HISTORY:
!  18 Jun 2013 - D.A. Ridley - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    CHARACTER(LEN=20)  :: FieldName
    CHARACTER(LEN=255) :: ErrMsg, ThisLoc

    !=================================================================
    ! AttachPointersFromHemco begins here!
    !=================================================================

    ! Initialize
    RC      = GC_SUCCESS
    ErrMsg  = ''
    ThisLoc = ' -> at AttachPointersFromHemco (in rrtmg_rad_transfer_mod.F90)'

    !-------------------------------
    ! CCl4 [ppb]
    !-------------------------------
    FieldName = 'TES_CLIM_CCL4'
    CALL HCO_GetPtr( HcoState, FieldName, CCL4CLIM,  RC )
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Cannot get pointer to ' // TRIM( FieldName )
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    !-------------------------------
    ! CFC11 [ppb]
    !-------------------------------
    FieldName = 'TES_CLIM_CFC11'
    CALL HCO_GetPtr( HcoState, FieldName, CFC11CLIM, RC )
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Cannot get pointer to ' // TRIM( FieldName )
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    !-------------------------------
    ! CFC12 [ppb]
    !-------------------------------
    FieldName = 'TES_CLIM_CFC12'
    CALL HCO_GetPtr( HcoState, FieldName, CFC12CLIM, RC  )
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Cannot get pointer to ' // TRIM( FieldName )
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    !-------------------------------
    ! CFC22 [ppb]
    !-------------------------------
    FieldName = 'TES_CLIM_CFC22'
    CALL HCO_GetPtr( HcoState, FieldName, CFC22CLIM, RC  )

    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Cannot get pointer to ' // TRIM( FieldName )
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    !-------------------------------
    ! CH4 [ppb]
    !-------------------------------
    FieldName = 'TES_CLIM_CH4'
    CALL HCO_GetPtr( HcoState, FieldName, CH4CLIM,   RC  )
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Cannot get pointer to ' // TRIM( FieldName )
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    !-------------------------------
    ! N2O [ppb]
    !-------------------------------
    FieldName = 'TES_CLIM_N2O'
    CALL HCO_GetPtr( HcoState, FieldName, N2OCLIM,   RC  )
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Cannot get pointer to ' // TRIM( FieldName )
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    !-------------------------------
    ! Diffuse Near-IR albedo [1]
    !-------------------------------
    FieldName = 'MODIS_ALBDFNIR'
    CALL HCO_GetPtr( HcoState, FieldName, MODIS_ALBDFNIR, RC  )
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Cannot get pointer to ' // TRIM( FieldName )
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    !-------------------------------
    ! Diffuse visible albedo [1]
    !-------------------------------
    FieldName = 'MODIS_ALBDFVIS'
    CALL HCO_GetPtr( HcoState, FieldName, MODIS_ALBDFVIS, RC  )
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Cannot get pointer to ' // TRIM( FieldName )
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    !-------------------------------
    ! Direct Near-IR albedo [1]
    !-------------------------------
    FieldName = 'MODIS_ALBDRNIR'
    CALL HCO_GetPtr( HcoState, FieldName, MODIS_ALBDRNIR, RC  )
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Cannot get pointer to ' // TRIM( FieldName )
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    !-------------------------------
    ! Direct visible albedo [1]
    !-------------------------------
    FieldName = 'MODIS_ALBDRVIS'
    CALL HCO_GetPtr( HcoState, FieldName, MODIS_ALBDRVIS, RC  )
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Cannot get pointer to ' // TRIM( FieldName )
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    !-------------------------------
    ! MODIS emissivity, band 1 [1]
    !-------------------------------
    FieldName = 'MODIS_EMISSIVITY_01'
    CALL HCO_GetPtr( HcoState, FieldName, MODIS_EMISS_01, RC  )
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Cannot get pointer to ' // TRIM( FieldName )
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    !-------------------------------
    ! MODIS emissivity, band 2 [1]
    !-------------------------------
    FieldName = 'MODIS_EMISSIVITY_02'
    CALL HCO_GetPtr( HcoState, FieldName, MODIS_EMISS_02, RC  )
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Cannot get pointer to ' // TRIM( FieldName )
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    !-------------------------------
    ! MODIS emissivity, band 3 [1]
    !-------------------------------
    FieldName = 'MODIS_EMISSIVITY_03'
    CALL HCO_GetPtr( HcoState, FieldName, MODIS_EMISS_03, RC  )
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Cannot get pointer to ' // TRIM( FieldName )
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    !-------------------------------
    ! MODIS emissivity, band 4 [1]
    !-------------------------------
    FieldName = 'MODIS_EMISSIVITY_04'
    CALL HCO_GetPtr( HcoState, FieldName, MODIS_EMISS_04, RC  )
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Cannot get pointer to ' // TRIM( FieldName )
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    !-------------------------------
    ! MODIS emissivity, band 5 [1]
    !-------------------------------
    FieldName = 'MODIS_EMISSIVITY_05'
    CALL HCO_GetPtr( HcoState, FieldName, MODIS_EMISS_05, RC  )
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Cannot get pointer to ' // TRIM( FieldName )
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    !-------------------------------
    ! MODIS emissivity, band 6 [1]
    !-------------------------------
    FieldName = 'MODIS_EMISSIVITY_06'
    CALL HCO_GetPtr( HcoState, FieldName, MODIS_EMISS_06, RC  )
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Cannot get pointer to ' // TRIM( FieldName )
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    !-------------------------------
    ! MODIS emissivity, band 7 [1]
    !-------------------------------
    FieldName = 'MODIS_EMISSIVITY_07'
    CALL HCO_GetPtr( HcoState, FieldName, MODIS_EMISS_07, RC  )
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Cannot get pointer to ' // TRIM( FieldName )
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    !-------------------------------
    ! MODIS emissivity, band 8 [1]
    !-------------------------------
    FieldName = 'MODIS_EMISSIVITY_08'
    CALL HCO_GetPtr( HcoState, FieldName, MODIS_EMISS_08, RC  )
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Cannot get pointer to ' // TRIM( FieldName )
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    !-------------------------------
    ! MODIS emissivity, band 9 [1]
    !-------------------------------
    FieldName = 'MODIS_EMISSIVITY_09'
    CALL HCO_GetPtr( HcoState, FieldName, MODIS_EMISS_09, RC  )
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Cannot get pointer to ' // TRIM( FieldName )
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    !-------------------------------
    ! MODIS emissivity, band 10 [1]
    !-------------------------------
    FieldName = 'MODIS_EMISSIVITY_10'
    CALL HCO_GetPtr( HcoState, FieldName, MODIS_EMISS_10, RC  )
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Cannot get pointer to ' // TRIM( FieldName )
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    !-------------------------------
    ! MODIS emissivity, band 11 [1]
    !-------------------------------
    FieldName = 'MODIS_EMISSIVITY_11'
    CALL HCO_GetPtr( HcoState, FieldName, MODIS_EMISS_11, RC  )
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Cannot get pointer to ' // TRIM( FieldName )
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    !-------------------------------
    ! MODIS emissivity, band 12 [1]
    !-------------------------------
    FieldName = 'MODIS_EMISSIVITY_12'
    CALL HCO_GetPtr( HcoState, FieldName, MODIS_EMISS_12, RC  )
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Cannot get pointer to ' // TRIM( FieldName )
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    !-------------------------------
    ! MODIS emissivity, band 13 [1]
    !-------------------------------
    FieldName = 'MODIS_EMISSIVITY_13'
    CALL HCO_GetPtr( HcoState, FieldName, MODIS_EMISS_13, RC  )
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Cannot get pointer to ' // TRIM( FieldName )
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    !-------------------------------
    ! MODIS emissivity, band 14 [1]
    !-------------------------------
    FieldName = 'MODIS_EMISSIVITY_14'
    CALL HCO_GetPtr( HcoState, FieldName, MODIS_EMISS_14, RC  )
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Cannot get pointer to ' // TRIM( FieldName )
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    !-------------------------------
    ! MODIS emissivity, band 15 [1]
    !-------------------------------
    FieldName = 'MODIS_EMISSIVITY_15'
    CALL HCO_GetPtr( HcoState, FieldName, MODIS_EMISS_15, RC  )
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Cannot get pointer to ' // TRIM( FieldName )
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    !-------------------------------
    ! MODIS emissivity, band 16 [1]
    !-------------------------------
    FieldName = 'MODIS_EMISSIVITY_16'
    CALL HCO_GetPtr( HcoState, FieldName, MODIS_EMISS_16, RC  )
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Cannot get pointer to ' // TRIM( FieldName )
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

  END SUBROUTINE AttachPointersFromHemco
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_mcica_clouds
!
! !DESCRIPTION: Subroutine INIT\_MCICA\_CLOUDS initializes all allocatable
!  module arrays.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE INIT_MCICA_CLOUDS( State_Grid )
!
! !USES:
!
    USE CMN_FJX_MOD
    USE ERROR_MOD,      ONLY : ALLOC_ERR
    USE PARRRTM,        ONLY : NGPTLW
    USE PARRRSW,        ONLY : NGPTSW
    USE State_Grid_Mod, ONLY : GrdState
!
! !INPUT PARAMETERS:
!
    TYPE(GrdState), INTENT(IN) :: State_Grid  ! Grid State object
!
! !REVISION HISTORY:
!  18 Jun 2013 - D.A. Ridley - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER :: AS

    !=================================================================
    ! INIT_MCICA_CLOUDS begins here!
    !=================================================================

    ALLOCATE( CLDFMCL_LW(NGPTLW,State_Grid%NX,State_Grid%NY,State_Grid%NZ ), &
              STAT=AS )
    IF ( AS /= 0 ) CALL ALLOC_ERR( 'CLDFMCL_LW' )
    CLDFMCL_LW = 0D0

    ALLOCATE( CIWPMCL_LW(NGPTLW, State_Grid%NX, State_Grid%NY, State_Grid%NZ ),&
              STAT=AS )
    IF ( AS /= 0 ) CALL ALLOC_ERR( 'CIWPMCL_LW' )
    CIWPMCL_LW = 0D0

    ALLOCATE( CLWPMCL_LW( NGPTLW, State_Grid%NX, State_Grid%NY, State_Grid%NZ ),&
              STAT=AS )
    IF ( AS /= 0 ) CALL ALLOC_ERR( 'CLWPMCL_LW' )
    CLWPMCL_LW = 0D0

    ALLOCATE( TAUCMCL_LW( NGPTLW, State_Grid%NX, State_Grid%NY,State_Grid%NZ ),&
              STAT=AS )
    IF ( AS /= 0 ) CALL ALLOC_ERR( 'TAUCMCL_LW' )
    TAUCMCL_LW = 0D0

    ALLOCATE( CLDFMCL_SW( NGPTSW, State_Grid%NX, State_Grid%NY, State_Grid%NZ ),&
              STAT=AS )
    IF ( AS /= 0 ) CALL ALLOC_ERR( 'CLDFMCL_SW' )
    CLDFMCL_SW = 0D0

    ALLOCATE( CIWPMCL_SW( NGPTSW, State_Grid%NX, State_Grid%NY, State_Grid%NZ ),&
              STAT=AS )
    IF ( AS /= 0 ) CALL ALLOC_ERR( 'CIWPMCL_SW' )
    CIWPMCL_SW = 0D0

    ALLOCATE( CLWPMCL_SW( NGPTSW, State_Grid%NX, State_Grid%NY, State_Grid%NZ ),&
              STAT=AS )
    IF ( AS /= 0 ) CALL ALLOC_ERR( 'CLWPMCL_SW' )
    CLWPMCL_SW = 0D0

    ALLOCATE( TAUCMCL_SW( NGPTSW, State_Grid%NX, State_Grid%NY, State_Grid%NZ ),&
              STAT=AS )
    IF ( AS /= 0 ) CALL ALLOC_ERR( 'TAUCMCL_SW' )
    TAUCMCL_SW = 0D0

    ALLOCATE( SSACMCL( NGPTSW, State_Grid%NX, State_Grid%NY, State_Grid%NZ ), &
              STAT=AS )
    IF ( AS /= 0 ) CALL ALLOC_ERR( 'SSACMCL' )
    SSACMCL = 0D0

    ALLOCATE( ASMCMCL( NGPTSW, State_Grid%NX, State_Grid%NY, State_Grid%NZ ), &
              STAT=AS )
    IF ( AS /= 0 ) CALL ALLOC_ERR( 'ASMCMCL' )
    ASMCMCL = 0D0

    ALLOCATE( FSFCMCL( NGPTSW, State_Grid%NX, State_Grid%NY, State_Grid%NZ ), &
              STAT=AS )
    IF ( AS /= 0 ) CALL ALLOC_ERR( 'FSFCMCL' )
    FSFCMCL = 0D0

    ALLOCATE( RELQMCL( State_Grid%NX, State_Grid%NY, State_Grid%NZ ), &
              STAT=AS )
    IF ( AS /= 0 ) CALL ALLOC_ERR( 'RELQMCL' )
    RELQMCL = 0D0

    ALLOCATE( REICMCL( State_Grid%NX, State_Grid%NY, State_Grid%NZ ), &
              STAT=AS )
    IF ( AS /= 0 ) CALL ALLOC_ERR( 'REICMCL' )
    REICMCL = 0D0

  END SUBROUTINE INIT_MCICA_CLOUDS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: cleanup_rrtmg_rad_transfer
!
! !DESCRIPTION: Deallocates all RRTMG module variables.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Cleanup_RRTMG_Rad_Transfer( RC )
!
! !USES:
!
    USE ErrCode_Mod
!
! !OUTPUT PARAMETERS:
!
    INTEGER, INTENT(OUT) :: RC          ! Success or failure
!
! !REVISION HISTORY:
!  09 Nov 2018 - R. Yantosca - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
    ! Assume success
    RC = GC_SUCCESS

    !=================================================================
    ! Nullify pointers to HEMCO fields
    !=================================================================

    ! Climatology
    CH4CLIM        => NULL()
    N2OCLIM        => NULL()
    CFC11CLIM      => NULL()
    CFC12CLIM      => NULL()
    CCL4CLIM       => NULL()
    CFC22CLIM      => NULL()

    ! Albedoes
    MODIS_ALBDFNIR => NULL()
    MODIS_ALBDFVIS => NULL()
    MODIS_ALBDRNIR => NULL()
    MODIS_ALBDRVIS => NULL()

    ! Emissivity
    MODIS_EMISS_01 => NULL()
    MODIS_EMISS_02 => NULL()
    MODIS_EMISS_03 => NULL()
    MODIS_EMISS_04 => NULL()
    MODIS_EMISS_05 => NULL()
    MODIS_EMISS_06 => NULL()
    MODIS_EMISS_07 => NULL()
    MODIS_EMISS_08 => NULL()
    MODIS_EMISS_09 => NULL()
    MODIS_EMISS_10 => NULL()
    MODIS_EMISS_11 => NULL()
    MODIS_EMISS_12 => NULL()
    MODIS_EMISS_13 => NULL()
    MODIS_EMISS_14 => NULL()
    MODIS_EMISS_15 => NULL()
    MODIS_EMISS_16 => NULL()

    !=================================================================
    ! Deallocate surface radiation arrays
    !=================================================================
    IF ( ALLOCATED( LW_UFLUX        ) ) DEALLOCATE( LW_UFLUX        )
    IF ( ALLOCATED( LW_DFLUX        ) ) DEALLOCATE( LW_DFLUX        )
    IF ( ALLOCATED( SW_UFLUX        ) ) DEALLOCATE( SW_UFLUX        )
    IF ( ALLOCATED( SW_DFLUX        ) ) DEALLOCATE( SW_DFLUX        )
    IF ( ALLOCATED( LW_UFLUXC       ) ) DEALLOCATE( LW_UFLUXC       )
    IF ( ALLOCATED( LW_DFLUXC       ) ) DEALLOCATE( LW_DFLUXC       )
    IF ( ALLOCATED( SW_UFLUXC       ) ) DEALLOCATE( SW_UFLUXC       )
    IF ( ALLOCATED( SW_DFLUXC       ) ) DEALLOCATE( SW_DFLUXC       )

    !=================================================================
    ! Deallocate MCICA cloud arrays
    !=================================================================
    IF ( ALLOCATED( CLDFMCL_LW     ) ) DEALLOCATE( CLDFMCL_LW       )
    IF ( ALLOCATED( CIWPMCL_LW     ) ) DEALLOCATE( CIWPMCL_LW       )
    IF ( ALLOCATED( CLWPMCL_LW     ) ) DEALLOCATE( CLWPMCL_LW       )
    IF ( ALLOCATED( TAUCMCL_LW     ) ) DEALLOCATE( TAUCMCL_LW       )
    IF ( ALLOCATED( CLDFMCL_SW     ) ) DEALLOCATE( CLDFMCL_SW       )
    IF ( ALLOCATED( CIWPMCL_SW     ) ) DEALLOCATE( CIWPMCL_SW       )
    IF ( ALLOCATED( CLWPMCL_SW     ) ) DEALLOCATE( CLWPMCL_SW       )
    IF ( ALLOCATED( TAUCMCL_SW     ) ) DEALLOCATE( TAUCMCL_SW       )
    IF ( ALLOCATED( SSACMCL        ) ) DEALLOCATE( SSACMCL          )
    IF ( ALLOCATED( ASMCMCL        ) ) DEALLOCATE( ASMCMCL          )
    IF ( ALLOCATED( FSFCMCL        ) ) DEALLOCATE( FSFCMCL          )
    IF ( ALLOCATED( REICMCL        ) ) DEALLOCATE( REICMCL          )
    IF ( ALLOCATED( RELQMCL        ) ) DEALLOCATE( RELQMCL          )

  END SUBROUTINE Cleanup_RRTMG_Rad_Transfer
!EOC
END MODULE RRTMG_RAD_TRANSFER_MOD
#endif
