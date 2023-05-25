!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: phot_container_mod.F90
!
! !DESCRIPTION: Module PHOT\_CONTAINER\_MOD contains the derived type used
!  to store photolysis and optics data in GEOS-Chem.
!\\
!\\
! !INTERFACE:
!
MODULE Phot_Container_Mod
!
! USES:
!
  USE ErrCode_Mod
  USE Precision_Mod

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC :: Init_Phot_Container
  PUBLIC :: Cleanup_Phot_Container
!
! !PUBLIC DATA MEMBERS:
!
  ! Parameters used for allocation

  !=========================================================================
  ! Derived type for photolysis and optics data container
  !=========================================================================
  TYPE, PUBLIC :: PhotContainer

     ! Scalars set during object initialization
     INTEGER :: IND999    ! Index in RAA & QAA of 999 nm
     INTEGER :: NWVAA     ! # LUT wavelengths   (RRTMG)
     INTEGER :: NSPAA     ! # LUT species       (RRTMG)
     INTEGER :: NRAA      ! # LUT aerosol sizes (RRTMG)   
     INTEGER :: NWVAA0    ! # non-RRTMG wavelengths
     INTEGER :: NALBD     ! ??                     
     INTEGER :: NEMISS    ! ??                     
     INTEGER :: NASPECRAD ! # RRTMG aerosol species
     INTEGER :: NSPECRAD  ! # RRTMG aerosol+gas species

     ! Scalars
     INTEGER  :: JTAUMX      ! max # divisions

     ! For RRTMG?
     INTEGER :: NWVAART      ! # RRTMG wavelengths
     INTEGER :: NWVREQUIRED  ! # WLs needed for interpolation
     INTEGER :: NRTWVREQUIRED! # WLs needed for RT interpolation
     INTEGER :: IWV1000      ! WL index for 1000 nm

     ! Renamed from Fast-JX module variables for clarity
     INTEGER  :: nWLbins      ! # WL bins (W_)
     INTEGER  :: nPhotRxns    ! # photolysis reactions in CTM chemistry (NRATJ)
     INTEGER  :: nMaxPhotRxns ! Maximum # of photolysis reactions (JVN_?)

     ! Photo-reaction flags for reactions adjusted in PhotRate_Adj
     INTEGER  :: RXN_O2      ! O2  + jv --> O   + O
     INTEGER  :: RXN_O3_1    ! O3  + hv --> O2  + O
     INTEGER  :: RXN_O3_2    ! O3  + hv --> O2  + O(1D)
     INTEGER  :: RXN_H2SO4   ! SO4 + hv --> SO2 + 2OH
     INTEGER  :: RXN_NO2     ! NO2 + hv --> NO  + O
     INTEGER  :: RXN_JHNO3   ! HNO3 + hv --> OH + NO2
     INTEGER  :: RXN_JNITSa  ! NITs  + hv --> HNO2
     INTEGER  :: RXN_JNITSb  ! NITs  + hv --> NO2
     INTEGER  :: RXN_JNITa   ! NIT + hv --> HNO2
     INTEGER  :: RXN_JNITb   ! NIT + hv --> NO2
     INTEGER  :: RXN_NO      ! For ucx_mod
     INTEGER  :: RXN_NO3     ! For ucx_mod
     INTEGER  :: RXN_N2O     ! For ucx_mod
     INTEGER  :: RXN_BrO     ! For Hg chem
     INTEGER  :: RXN_ClO     ! For Hg chem

     ! Arrays
     INTEGER,  ALLOCATABLE :: RINDEX     (:) ! GC to UCI spc name index mapping
     INTEGER,  ALLOCATABLE :: GC_Photo_Id(:) ! GC id per photolysis species
     INTEGER,  ALLOCATABLE :: MIEDX      (:) ! Interface indices for GC/FJX spc

     REAL(fp), ALLOCATABLE :: UVXFACTOR(:) ! Photons/cm2s -> W/m2 conv factors
     REAL(fp), ALLOCATABLE :: QAA_AOD  (:) ! Single scattering albedo        
     REAL(fp), ALLOCATABLE :: WAA_AOD  (:) ! Aerosol scattering phase fnctns 
     REAL(fp), ALLOCATABLE :: PAA_AOD  (:) ! WLs for supplied phase functions
     REAL(fp), ALLOCATABLE :: RAA_AOD  (:) ! Phase fnctn (first 8 terms)     
     REAL(fp), ALLOCATABLE :: SAA_AOD  (:) ! Aerosol type effective radius

     REAL(fp), ALLOCATABLE :: TREF     (:,:,:)     ! Temp reference profile
     REAL(fp), ALLOCATABLE :: OREF     (:,:,:)     ! Ozone reference profile
     REAL(fp), ALLOCATABLE :: ZPJ      (:,:,:,:)   ! J-values

     ! RRTMG allocatable arrays
     INTEGER, ALLOCATABLE :: SPECMASK     (:)     ! binary switches for spc flux
     INTEGER, ALLOCATABLE :: IWVREQUIRED  (:)     ! WL indexes for interpolation
     INTEGER, ALLOCATABLE :: IRTWVREQUIRED(:)     ! WL indexes for RT interp
     INTEGER, ALLOCATABLE :: IWVSELECT    (:,:)   ! Indexes of requested WLs
     INTEGER, ALLOCATABLE :: IRTWVSELECT  (:,:)   ! Indexes of requested RT WLs
     INTEGER, ALLOCATABLE :: IRHARR       (:,:,:) ! Relative humidity indices

     REAL*8,  ALLOCATABLE :: ACOEF_WV  (:)   ! Coeffs for WL interpolation
     REAL*8,  ALLOCATABLE :: BCOEF_WV  (:)   ! Coeffs for WL interpolation
     REAL*8,  ALLOCATABLE :: CCOEF_WV  (:)   ! Coeffs for WL interpolation
     REAL*8,  ALLOCATABLE :: ACOEF_RTWV(:)   ! Coeffs for RT WL interpolation
     REAL*8,  ALLOCATABLE :: BCOEF_RTWV(:)   ! Coeffs for RT WL interpolation
     REAL*8,  ALLOCATABLE :: CCOEF_RTWV(:)   ! Coeffs for RT WL interpolation
     REAL*8,  ALLOCATABLE :: WVAA      (:,:)     ! ??
     REAL*8,  ALLOCATABLE :: RHAA      (:,:)     ! ??
     REAL*8,  ALLOCATABLE :: RDAA      (:,:)     ! ??
     REAL*8,  ALLOCATABLE :: RWAA      (:,:)     ! ??
     REAL*8,  ALLOCATABLE :: SGAA      (:,:)     ! ??
     REAL*8,  ALLOCATABLE :: REAA      (:,:)     ! ??
     REAL*8,  ALLOCATABLE :: NRLAA     (:,:,:)   ! ??
     REAL*8,  ALLOCATABLE :: NCMAA     (:,:,:)   ! ??
     REAL*8,  ALLOCATABLE :: QQAA      (:,:,:)   ! ??
     REAL*8,  ALLOCATABLE :: ALPHAA    (:,:,:)   ! ??
     REAL*8,  ALLOCATABLE :: SSAA      (:,:,:)   ! ??
     REAL*8,  ALLOCATABLE :: ASYMAA    (:,:,:)   ! ??
     REAL*8,  ALLOCATABLE :: PHAA      (:,:,:,:) ! ??

     ! For optical depth diagnostics
     REAL(fp), ALLOCATABLE :: ISOPOD   (:,:,:,:)   ! Isoprene optical depth
     REAL(fp), ALLOCATABLE :: ODMDUST  (:,:,:,:,:) ! Dust optical depth
     REAL(fp), ALLOCATABLE :: ODAER    (:,:,:,:,:) ! Aerosol optical depth

#ifdef RRTMG
     REAL*8,  ALLOCATABLE :: RTODAER   (:,:,:,:,:) ! Optical dust
     REAL*8,  ALLOCATABLE :: RTSSAER   (:,:,:,:,:) ! ??
     REAL*8,  ALLOCATABLE :: RTASYMAER (:,:,:,:,:) ! ??
#endif

  END TYPE PhotContainer
!
! !REMARKS:
! 
! !REVISION HISTORY:
!  28 Nov 2022 - E. Lundgren- Initial version, based on state_grid_mod.F90
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Init_Phot_Container
!
! !DESCRIPTION: Subroutine INIT\_PHOT\_Container allocates and initializes
! the Phot container object.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Init_Phot_Container( Input_Opt, State_Grid, Phot, RC )
!
! !USES:
!
    USE CMN_FJX_Mod,    ONLY : A_, AN_, W_, WX_, JVN_, N_, JXL_
    USE CMN_Size_Mod,   ONLY : NDUST, NAER
    USE ErrCode_Mod
    USE Input_Opt_Mod,  ONLY : OptInput
    USE State_Grid_Mod, ONLY : GrdState
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput),      INTENT(IN)  :: Input_Opt  ! Input Options object
    TYPE(GrdState),      INTENT(IN)  :: State_Grid ! Grid object
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(PhotContainer), POINTER     :: Phot       ! Phot data container
!
! !OUTPUT PARAMETERS:
!
    INTEGER,             INTENT(OUT) :: RC         ! Success or failure?
!
! !REVISION HISTORY:
!  28 Nov 2022 - E. Lundgren- Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    CHARACTER(LEN=255)   :: errMsg, thisLoc

    !======================================================================
    ! Allocate and initialize module variables
    !======================================================================

    ! Assume success
    RC      = GC_SUCCESS
    thisLoc = ' -> at Init_Phot_Container (in module Headers/phot_container_mod.F90)'

    ! Constants
    Phot%IND999 = 5     ! Index in RAA & QAA of 999 nm
    Phot%NWVAA  = 41    ! # LUT wavelengths   (RRTMG)
    Phot%NSPAA  = 8     ! # LUT species       (RRTMG)
    Phot%NRAA   = 7     ! # LUT aerosol sizes (RRTMG)   
    Phot%NWVAA0 = 11    ! # non-RRTMG wavelengths
    Phot%NALBD  = 2     ! ??                     
    Phot%NEMISS = 16    ! ??                     
    Phot%NASPECRAD = 16 ! # RRTMG aerosol species
    Phot%NSPECRAD  = 18 ! # RRTMG aerosol+gas species

    ! Store certain values from Fast-JX with more intuitive name
    Phot%nWLbins      = W_
    Phot%nPhotRxns    = 0 ! Set during photolysis initialization
    Phot%nMaxPhotRxns = JVN_! Maximum # of photolysis reactions (JVN_?)

    ! Integer scalars
    Phot%JTAUMX = (N_-4*JXL_)/2

    Phot%RXN_O2     = -1
    Phot%RXN_O3_1   = -1
    Phot%RXN_O3_2   = -1
    Phot%RXN_H2SO4  = -1
    Phot%RXN_NO2    = -1
    Phot%RXN_JHNO3  = -1
    Phot%RXN_JNITSa = -1
    Phot%RXN_JNITSb = -1
    Phot%RXN_JNITa  = -1
    Phot%RXN_JNITb  = -1
    Phot%RXN_NO     = -1
    Phot%RXN_NO3    = -1
    Phot%RXN_N2O    = -1
    Phot%RXN_BrO    = -1
    Phot%RXN_ClO    = -1

    ! Allocate arrays
    IF ( .not. Input_Opt%DryRun ) THEN

       ! Integer arrays

       ! Phot%RINDEX      (:)
       ALLOCATE( Phot%RINDEX( Phot%nMaxPhotRxns ), STAT=RC )
       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = 'Error allocating array RINDEX!'
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF
       Phot%RINDEX = 0

       ! Phot%GC_Photo_Id (:)
       ALLOCATE( Phot%GC_Photo_Id( Phot%nMaxPhotRxns ), STAT=RC )
       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = 'Error allocating array GC_Photo_Id!'
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF
       Phot%GC_Photo_Id = 0

       ! Phot%MIEDX       (:)
       ALLOCATE( Phot%MIEDX( AN_ ), STAT=RC )
       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = 'Error allocating array MIEDX!'
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF
       Phot%MIEDX = 0

       ! Real(fp) arrays

       ! Phot%UVXFACTOR(:)
       ALLOCATE( Phot%UVXFACTOR( WX_ ), STAT=RC )
       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = 'Error allocating array UVXFACTOR!'
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF
       Phot%UVXFACTOR = 0e+0_fp

       ! Phot%QAA_AOD  (:)
       ALLOCATE( Phot%QAA_AOD( A_ ), STAT=RC )
       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = 'Error allocating array QAA_AOD!'
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF
       Phot%QAA_AOD = 0e+0_fp

       ! Phot%WAA_AOD  (:)
       ALLOCATE( Phot%WAA_AOD( A_ ), STAT=RC )
       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = 'Error allocating array WAA_AOD!'
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF
       Phot%WAA_AOD = 0e+0_fp

       ! Phot%PAA_AOD  (:)
       ALLOCATE( Phot%PAA_AOD( A_ ), STAT=RC )
       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = 'Error allocating array PAA_AOD!'
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF
       Phot%PAA_AOD = 0e+0_fp

       ! Phot%RAA_AOD  (:)
       ALLOCATE( Phot%RAA_AOD( A_ ), STAT=RC )
       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = 'Error allocating array RAA_AOD!'
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF
       Phot%RAA_AOD = 0e+0_fp

       ! Phot%SAA_AOD  (:)
       ALLOCATE( Phot%SAA_AOD( A_ ), STAT=RC )
       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = 'Error allocating array SAA_AOD!'
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF
       Phot%SAA_AOD = 0e+0_fp

       ! Phot%TREF     (:,:,:)
       ALLOCATE( Phot%TREF( 51, 18, 12 ), STAT=RC )
       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = 'Error allocating array TREF!'
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF
       Phot%TREF = 0e+0_fp

       ! Phot%OREF     (:,:,:)
       ALLOCATE( Phot%OREF( 51, 18, 12 ), STAT=RC )
       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = 'Error allocating array OREF!'
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF
       Phot%OREF = 0e+0_fp

       ! Phot%ZPJ      (:,:,:,:)
       ALLOCATE( Phot%ZPJ( State_Grid%NZ, Phot%nMaxPhotRxns, State_Grid%NX, &
                 State_Grid%NY ), STAT=RC )
       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = 'Error allocating array ZPJ!'
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF
       Phot%ZPJ = 0e+0_fp

    ENDIF

    !--------------------------------------------------
    ! Fields for RRTMG and optical depth diagnostics
    !--------------------------------------------------
    
    ! Other scalars
    Phot%NWVAART = Phot%NWVAA-Phot%NWVAA0 ! # RRTMG wavelengths    

    ! Scalars set in subroutine CALC_AOD
    Phot%NWVREQUIRED = 0
    Phot%NRTWVREQUIRED = 0

    ! Scalars set in subroutine RD_AOD
    Phot%IWV1000 = 0

    ! Allocate arrays
    IF ( .not. Input_Opt%DryRun ) THEN
    
       ! RRTMG integer arrays

       ! Phot%SPECMASK     (:)
       ALLOCATE( Phot%SPECMASK( AN_ ), STAT=RC )
       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = 'Error allocating array SPECMASK!'
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF
       Phot%SPECMASK = 0
    
       ! Phot%IWVREQUIRED  (:)
       ALLOCATE( Phot%IWVREQUIRED( AN_ ), STAT=RC )
       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = 'Error allocating array IWVREQUIRED!'
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF
       Phot%IWVREQUIRED = 0
    
       ! Phot%IRTWVREQUIRED(:)
       ALLOCATE( Phot%IRTWVREQUIRED( AN_ ), STAT=RC )
       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = 'Error allocating array IRTWVREQUIRED!'
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF 
       Phot%IRTWVREQUIRED = 0
    
       ! Phot%IWVSELECT    (:,:)
       ALLOCATE( Phot%IWVSELECT( 2, 3 ), STAT=RC )
       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = 'Error allocating array IWVSELECT!'
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF
       Phot%IWVSELECT = 0
    
       ! Phot%IRTWVSELECT  (:,:)
       ALLOCATE( Phot%IRTWVSELECT( 2, 3 ), STAT=RC )
       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = 'Error allocating array IRTWVSELECT!'
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF
       Phot%IRTWVSELECT = 0
    
       ! Phot%IRHARR    (:,:,:)
       ALLOCATE( Phot%IRHARR( State_Grid%NX, State_Grid%NY, &
                              State_Grid%NZ ), STAT=RC )
       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = 'Error allocating array IRHARR!'
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF
       Phot%IRHARR = 0d0
    
       ! RRTMG real*8 arrays
    
       ! Phot%ACOEF_WV  (:)
       ALLOCATE( Phot%ACOEF_WV( AN_ ), STAT=RC )
       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = 'Error allocating array ACOEF_WV!'
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF
       Phot%ACOEF_WV = 0d0
    
       ! Phot%BCOEF_WV  (:)
       ALLOCATE( Phot%BCOEF_WV( AN_ ), STAT=RC )
       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = 'Error allocating array BCOEF_WV!'
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF
       Phot%BCOEF_WV = 0d0
    
       ! Phot%CCOEF_WV  (:)
       ALLOCATE( Phot%CCOEF_WV( AN_ ), STAT=RC )
       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = 'Error allocating array CCOEF_WV!'
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF
       Phot%CCOEF_WV = 0d0
    
       ! Phot%ACOEF_RTWV(:)
       ALLOCATE( Phot%ACOEF_RTWV( AN_ ), STAT=RC )
       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = 'Error allocating array ACOEF_RTWV!'
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF
       Phot%ACOEF_RTWV = 0d0
    
       ! Phot%BCOEF_RTWV(:)
       ALLOCATE( Phot%BCOEF_RTWV( AN_ ), STAT=RC )
       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = 'Error allocating array BCOEF_RTWV!'
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF
       Phot%BCOEF_RTWV = 0d0
    
       ! Phot%CCOEF_RTWV(:)
       ALLOCATE( Phot%CCOEF_RTWV( AN_ ), STAT=RC )
       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = 'Error allocating array CCOEF_RTWV!'
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF
       Phot%CCOEF_RTWV = 0d0
    
       ! Phot%WVAA      (:,:)
       ALLOCATE( Phot%WVAA( Phot%NWVAA, Phot%NSPAA ), STAT=RC )
       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = 'Error allocating array WVAA!'
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF
       Phot%WVAA = 0d0
    
       ! Phot%RHAA      (:,:)
       ALLOCATE( Phot%RHAA( Phot%NRAA, Phot%NSPAA ), STAT=RC )
       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = 'Error allocating array RHAA!'
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF
       Phot%RHAA = 0d0
    
       ! Phot%RDAA      (:,:)!
       ALLOCATE( Phot%RDAA( Phot%NRAA, Phot%NSPAA ), STAT=RC )
       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = 'Error allocating array RDAA!'
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF
       Phot%RDAA = 0d0
    
       ! Phot%RWAA      (:,:)
       ALLOCATE( Phot%RWAA( Phot%NRAA, Phot%NSPAA ), STAT=RC )
       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = 'Error allocating array RWAA!'
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF
       Phot%RWAA = 0d0
    
       ! Phot%SGAA      (:,:)
       ALLOCATE( Phot%SGAA( Phot%NRAA, Phot%NSPAA ), STAT=RC )
       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = 'Error allocating array SGAA!'
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF
       Phot%SGAA = 0d0
    
       ! Phot%REAA      (:,:)
       ALLOCATE( Phot%REAA( Phot%NRAA, Phot%NSPAA ), STAT=RC )
       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = 'Error allocating array REAA!'
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF
       Phot%REAA = 0d0
    
       ! Phot%NRLAA     (:,:,:)
       ALLOCATE( Phot%NRLAA( Phot%NWVAA, Phot%NRAA, Phot%NSPAA ), STAT=RC )
       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = 'Error allocating array NRLAA!'
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF
       Phot%NRLAA = 0d0
    
       ! Phot%NCMAA     (:,:,:)
       ALLOCATE( Phot%NCMAA( Phot%NWVAA, Phot%NRAA, Phot%NSPAA ), STAT=RC )
       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = 'Error allocating array NCMAA!'
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF
       Phot%NCMAA = 0d0
    
       ! Phot%QQAA      (:,:,:)
       ALLOCATE( Phot%QQAA( Phot%NWVAA, Phot%NRAA, Phot%NSPAA ), STAT=RC )
       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = 'Error allocating array QQAA!'
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF
       Phot%QQAA = 0d0
    
       ! Phot%ALPHAA    (:,:,:)
       ALLOCATE( Phot%ALPHAA( Phot%NWVAA, Phot%NRAA, Phot%NSPAA ), &
                 STAT=RC )
       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = 'Error allocating array ALPHAA!'
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF
       Phot%ALPHAA = 0d0
    
       ! Phot%SSAA      (:,:,:)
       ALLOCATE( Phot%SSAA( Phot%NWVAA, Phot%NRAA, Phot%NSPAA ), STAT=RC )
       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = 'Error allocating array SSAA!'
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF
       Phot%SSAA = 0d0
    
       ! Phot%ASYMAA    (:,:,:)
       ALLOCATE( Phot%ASYMAA( Phot%NWVAA, Phot%NRAA, Phot%NSPAA ), &
                 STAT=RC )
       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = 'Error allocating array ASYMAA!'
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF
       Phot%ASYMAA = 0d0
    
       ! Phot%PHAA      (:,:,:,:)
       ALLOCATE( Phot%PHAA( Phot%NWVAA, Phot%NRAA, Phot%NSPAA, 8 ), &
                 STAT=RC )
       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = 'Error allocating array PHAA!'
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF
       Phot%PHAA = 0d0

       ! Phot%ISOPOD   (:,:,:,:)
       ALLOCATE( Phot%ISOPOD( State_Grid%NX, State_Grid%NY, State_Grid%NZ, &
                              Phot%NWVAA ), STAT=RC )
       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = 'Error allocating array ISOPOD!'
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF
       Phot%ISOPOD = 0e+0_fp
       
       ! Phot%ODMDUST  (:,:,:,:,:)
       ALLOCATE( Phot%ODMDUST( State_Grid%NX, State_Grid%NY, State_Grid%NZ, &
                               Phot%NWVAA, NDUST ), STAT=RC )
       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = 'Error allocating array ODMDUST!'
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF
       Phot%ODMDUST = 0e+0_fp
       
       ! Phot%ODAER    (:,:,:,:,:)
       ALLOCATE( Phot%ODAER( State_Grid%NX, State_Grid%NY, State_Grid%NZ, &
                             Phot%NWVAA, NAER ), STAT=RC )
       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = 'Error allocating array ODAER!'
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF
       Phot%ODAER = 0e+0_fp

#ifdef RRTMG
       ! Phot%RTODAER   (:,:,:,:,:)
       ! +2 to split SNA into SU, NI and AM
       ALLOCATE( Phot%RTODAER( State_Grid%NX, State_Grid%NY, State_Grid%NZ, &
                               Phot%NWVAA, NAER+2+NDUST ), STAT=RC )
       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = 'Error allocating array RTODAER!'
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF
       Phot%RTODAER = 0d0

       ! Phot%RTSSAER   (:,:,:,:,:)
       ALLOCATE( Phot%RTSSAER( State_Grid%NX, State_Grid%NY, State_Grid%NZ, &
                               Phot%NWVAA, NAER+2+NDUST ), STAT=RC )
       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = 'Error allocating array RTSSAER!'
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF
       Phot%RTSSAER = 0d0
       
       ! Phot%RTASYMAER (:,:,:,:,:)
       ALLOCATE( Phot%RTASYMAER( State_Grid%NX, State_Grid%NY, State_Grid%NZ, &
                                 Phot%NWVAA, NAER+2+NDUST ), STAT=RC )
       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = 'Error allocating array RTASYMAER!'
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF
       Phot%RTASYMAER = 0d0
#endif
     
    ENDIF

  END SUBROUTINE Init_Phot_Container
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Cleanup_Phot_Container
!
! !DESCRIPTION: Subroutine CLEANUP\_PHOT\_CONTAINER deallocates all fields
!  of the phot container object.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Cleanup_Phot_Container( Phot, RC )
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(PhotContainer), POINTER   :: Phot  ! Phot data container
!
! !OUTPUT PARAMETERS:
!
    INTEGER,           INTENT(OUT) :: RC    ! Return code
!
! !REVISION HISTORY:
!  28 Nov 2022 - E. Lundgren- Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC

    ! Assume success
    RC      = GC_SUCCESS

    !=======================================================================
    ! Deallocate arrays
    !=======================================================================
    ! Will need to change this to just do arrays etc
    IF ( ASSOCIATED( Phot ) ) THEN
       IF (ALLOCATED(Phot%RINDEX        )) DEALLOCATE(Phot%RINDEX        )
       IF (ALLOCATED(Phot%GC_Photo_Id   )) DEALLOCATE(Phot%GC_Photo_Id   )
       IF (ALLOCATED(Phot%MIEDX         )) DEALLOCATE(Phot%MIEDX         )
       IF (ALLOCATED(Phot%UVXFACTOR     )) DEALLOCATE(Phot%UVXFACTOR     )
       IF (ALLOCATED(Phot%QAA_AOD       )) DEALLOCATE(Phot%QAA_AOD       )
       IF (ALLOCATED(Phot%WAA_AOD       )) DEALLOCATE(Phot%WAA_AOD       )
       IF (ALLOCATED(Phot%PAA_AOD       )) DEALLOCATE(Phot%PAA_AOD       )
       IF (ALLOCATED(Phot%RAA_AOD       )) DEALLOCATE(Phot%RAA_AOD       )
       IF (ALLOCATED(Phot%SAA_AOD       )) DEALLOCATE(Phot%SAA_AOD       )
       IF (ALLOCATED(Phot%TREF          )) DEALLOCATE(Phot%TREF          )
       IF (ALLOCATED(Phot%OREF          )) DEALLOCATE(Phot%OREF          )
       IF (ALLOCATED(Phot%ISOPOD        )) DEALLOCATE(Phot%ISOPOD        )
       IF (ALLOCATED(Phot%ZPJ           )) DEALLOCATE(Phot%ZPJ           )
       IF (ALLOCATED(Phot%ODMDUST       )) DEALLOCATE(Phot%ODMDUST       )
       IF (ALLOCATED(Phot%ODAER         )) DEALLOCATE(Phot%ODAER         )
       IF (ALLOCATED(Phot%SPECMASK      )) DEALLOCATE(Phot%SPECMASK      )
       IF (ALLOCATED(Phot%IWVREQUIRED   )) DEALLOCATE(Phot%IWVREQUIRED   )
       IF (ALLOCATED(Phot%IRTWVREQUIRED )) DEALLOCATE(Phot%IRTWVREQUIRED )
       IF (ALLOCATED(Phot%IWVSELECT     )) DEALLOCATE(Phot%IWVSELECT     )
       IF (ALLOCATED(Phot%IRTWVSELECT   )) DEALLOCATE(Phot%IRTWVSELECT   )
       IF (ALLOCATED(Phot%IRHARR        )) DEALLOCATE(Phot%IRHARR        )
       IF (ALLOCATED(Phot%ACOEF_WV      )) DEALLOCATE(Phot%ACOEF_WV      )
       IF (ALLOCATED(Phot%BCOEF_WV      )) DEALLOCATE(Phot%BCOEF_WV      )
       IF (ALLOCATED(Phot%CCOEF_WV      )) DEALLOCATE(Phot%CCOEF_WV      )
       IF (ALLOCATED(Phot%ACOEF_RTWV    )) DEALLOCATE(Phot%ACOEF_RTWV    )
       IF (ALLOCATED(Phot%BCOEF_RTWV    )) DEALLOCATE(Phot%BCOEF_RTWV    )
       IF (ALLOCATED(Phot%CCOEF_RTWV    )) DEALLOCATE(Phot%CCOEF_RTWV    )
       IF (ALLOCATED(Phot%WVAA          )) DEALLOCATE(Phot%WVAA          )
       IF (ALLOCATED(Phot%RHAA          )) DEALLOCATE(Phot%RHAA          )
       IF (ALLOCATED(Phot%RDAA          )) DEALLOCATE(Phot%RDAA          )
       IF (ALLOCATED(Phot%RWAA          )) DEALLOCATE(Phot%RWAA          )
       IF (ALLOCATED(Phot%SGAA          )) DEALLOCATE(Phot%SGAA          )
       IF (ALLOCATED(Phot%REAA          )) DEALLOCATE(Phot%REAA          )
       IF (ALLOCATED(Phot%NRLAA         )) DEALLOCATE(Phot%NRLAA         )
       IF (ALLOCATED(Phot%NCMAA         )) DEALLOCATE(Phot%NCMAA         )
       IF (ALLOCATED(Phot%QQAA          )) DEALLOCATE(Phot%QQAA          )
       IF (ALLOCATED(Phot%ALPHAA        )) DEALLOCATE(Phot%ALPHAA        )
       IF (ALLOCATED(Phot%SSAA          )) DEALLOCATE(Phot%SSAA          )
       IF (ALLOCATED(Phot%ASYMAA        )) DEALLOCATE(Phot%ASYMAA        )
       IF (ALLOCATED(Phot%PHAA          )) DEALLOCATE(Phot%PHAA          )
#ifdef RRTMG 
       IF (ALLOCATED(Phot%RTODAER       )) DEALLOCATE(Phot%RTODAER   )
       IF (ALLOCATED(Phot%RTSSAER       )) DEALLOCATE(Phot%RTSSAER   )
       IF (ALLOCATED(Phot%RTASYMAER     )) DEALLOCATE(Phot%RTASYMAER )
#endif

       Phot => NULL()
    ENDIF

  END SUBROUTINE Cleanup_Phot_Container
!EOC

END MODULE Phot_Container_Mod
