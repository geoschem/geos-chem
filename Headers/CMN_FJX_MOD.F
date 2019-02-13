!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: CMN_FJX_mod.F
!
! !DESCRIPTION: Module CMN\_FJX\_MOD contains parameters and global variables
!  used to interface between Harvard chemistry and UC-Irvine photolysis 
!  programs (Fast-J/Fast-JX), along with all Fast-J(X) global variables
!  and some physical constants for the GEOS-Chem chemistry code.
!\\
!\\
! !INTERFACE: 
!
      MODULE CMN_FJX_MOD
!
! !USES:
!
      USE CMN_SIZE_MOD,      ONLY : IIPAR, JJPAR, LLPAR, LLCHEM
      USE CMN_SIZE_MOD,      ONLY : NDUST, NAER, NRHAER, NSTRATAER, NRH

      USE PRECISION_MOD    ! For GEOS-Chem Precision (fp)

      IMPLICIT NONE
      PUBLIC
!
! !DEFINED PARAMETERS:
!
      ! New (SDE 03/28/13)
      ! Index in RAA & QAA of 999 nm wavelength
      INTEGER, PARAMETER :: IND999  = 5

      ! Required size of aerosol arrays
      INTEGER            :: L_             ! Number of CTM layers

      INTEGER            :: L1_            ! Number of CTM layer edges

      INTEGER            :: L2_            ! Number of levels in FJX grid that
                                           ! inc. both edges and mid-points

      INTEGER            :: JVL_           ! Vertical levels for J-values

      INTEGER, PARAMETER :: JVN_ = 130     ! Max number of J-values

      INTEGER            :: AN_            ! # of separate aerosols per layer
                                           ! Now set in Init_CMN_FJX below

      ! Variables used to interface GEOS-Chem and Fast-JX at runtime
      ! Branches for photolysis species
      INTEGER            :: BRANCH(JVN_)

      ! Names of photolysis species
      CHARACTER (LEN=10) :: RNAMES(JVN_)

      ! Mapping array from Harvard species names to UCI species names
      INTEGER            :: RINDEX(JVN_)                

      ! GEOS-Chem "ModelId" corresponding to each photolysis species
      INTEGER            :: GC_Photo_Id(JVN_)    

      ! Output J values
      REAL(fp), ALLOCATABLE :: ZPJ(:,:,:,:)

      !-----------------------------------------------------------------------
      ! variables used to map fast-JX J's onto CTM J's
      !-----------------------------------------------------------------------

      ! Multiplication factor for fast-JX calculated J
      REAL(fp)             :: JFACTA(JVN_)

      ! Index arrays that map Jvalue(j) onto rates
      INTEGER            :: JIND(JVN_)

      ! Mumber of Photolysis reactions in CTM chemistry, derived here NRATJ
      ! must be .le. JVN_
      INTEGER            :: NRATJ

      ! Label of J-value used in the main chem model
      CHARACTER*50       :: JLABEL(JVN_)

      ! JXL_: vertical(levels) dim for J-values computed within fast-JX
      INTEGER            ::  JXL_
      INTEGER            ::  JXL1_

      ! JXL2_: 2*JXL_ + 2 = mx no. levels in the basic Fast-JX grid (mid-level)
      INTEGER            ::  JXL2_

      ! WX_  = dim = no. of wavelengths in input file
      INTEGER, PARAMETER ::  WX_   = 18

      ! X_   = dim = max no. of X-section data sets (input data)
      ! SDE 2016-11-04: Increased from 75 to 123 for new halogen
      ! chemistry (iodine cross-sections)
      INTEGER, PARAMETER ::  X_    = 123

      ! A_   = dim = no. of Aerosol/cloud Mie sets (input data)
      INTEGER, PARAMETER ::  A_    = 56

      ! W_   = dim = no. of Wavelength bins:  =18 std, =12 trop only
      ! Now set below in Init_CMN_FJX (mps, 1/3/18)
      INTEGER            ::  W_

      ! N_  = no. of levels in Mie scattering arrays
      !     = 2*NC+1 = 4*(L_+1) + 1 + 2*sum(JADDLV)
#if defined( MODEL_GEOS )
      !!!INTEGER, PARAMETER ::  N_    = 601
      !!!INTEGER, PARAMETER ::  N_    = 1201
      INTEGER            :: N_
#else
      INTEGER, PARAMETER ::  N_    = 601
#endif

      ! M_  = no. of Gauss points used, must = 4 in fast_JX (no option)
      INTEGER, PARAMETER ::  M_    = 4

      ! M2_ = 2*M_ = 8, replaces MFIT
      INTEGER, PARAMETER ::  M2_   = 2*M_

      !-----------------------------------------------------------------------
      ! 4 Gauss pts = 8-stream
      REAL(fp), DIMENSION(M_), PARAMETER  ::
     &             EMU = [.06943184420297e+0_fp, .33000947820757e+0_fp,
     &                    .66999052179243e+0_fp, .93056815579703e+0_fp]
      REAL(fp), DIMENSION(M_), PARAMETER  ::
     &             WT  = [.17392742256873e+0_fp, .32607257743127e+0_fp,
     &                    .32607257743127e+0_fp, .17392742256873e+0_fp]
      !-----------------------------------------------------------------------

      ! ZZHT: scale height (cm)
      REAL(fp), PARAMETER  :: ZZHT   = 5.e+5_fp

      ! RAD: Radius of Earth (cm)
      REAL(fp), PARAMETER  :: RAD    = 6375.e+5_fp

      ! ATAU: heating rate (factor increase from one layer to the next)
      REAL(fp), PARAMETER  :: ATAU   = 1.120e+0_fp
#if defined( MODEL_GEOS )
      !REAL(fp), PARAMETER  :: ATAU   = 1.180e+0_fp
#endif

      ! ATAU0: minimum heating rate
      REAL(fp), PARAMETER  :: ATAU0  = 0.010e+0_fp

      ! JTAUMX = maximum number of divisions (i.e., may not get to ATAUMN)
      INTEGER              :: JTAUMX

      ! Physical constants
      REAL(fp), PARAMETER  :: UVXPLANCK = 6.62606957e-34
      REAL(fp), PARAMETER  :: UVXCCONST = 2.99792458e8

      ! Conversion factors from photons/cm2s to W/m2
      REAL(fp), DIMENSION(WX_) :: UVXFACTOR

      !-----------------------------------------------------------------------
      ! Variables in file 'FJX_spec.dat' (RD_XXX)
      !-----------------------------------------------------------------------

      ! WBIN: Boundaries of wavelength bins
      REAL(fp)             :: WBIN(WX_+1)

      ! WL: Centres of wavelength bins - 'effective wavelength'
      REAL(fp)             :: WL(WX_)

      ! FL: Solar flux incident on top of atmosphere (cm-2.s-1)
      REAL(fp)             :: FL(WX_)

      REAL(fp)             :: QO2(WX_,3)   ! QO2: O2 cross-sections
      REAL(fp)             :: QO3(WX_,3)   ! QO3: O3 cross-sections
      REAL(fp)             :: Q1D(WX_,3)   ! Q1D: O3 => O(1D) quantum yield

      ! QQQ: Supplied cross sections in each wavelength bin (cm2)
      REAL(fp)             :: QQQ(WX_,3,X_)

      ! QRAYL: Rayleigh parameters (effective cross-section) (cm2)
      REAL(fp)             :: QRAYL(WX_+1)

      ! TQQ: Temperature for supplied cross sections
      REAL(fp)             :: TQQ(3,X_)

      ! LQQ = 1, 2, or 3 to determine interpolation with T or P
      INTEGER              :: LQQ(X_)

      ! TITLEJX: Title for supplied cross sections, from 'FJX_spec.dat'
      CHARACTER*6          :: TITLEJX(X_)

      ! SQQ: Flag for supplied cross sections, from 'FJX_spec.dat'
      CHARACTER*1          :: SQQ(X_)

      !-----------------------------------------------------------------------
      ! Variables in file 'jv_spec_mie.dat' (RD_MIE)
      !-----------------------------------------------------------------------

      ! QAA: Aerosol scattering phase functions
      REAL(fp)             :: QAA(5,A_)

      ! WAA: 5 Wavelengths for the supplied phase functions
      REAL(fp)             :: WAA(5,A_)

      ! PAA: Phase function: first 8 terms of expansion
      REAL(fp)             :: PAA(8,5,A_)

      ! RAA: Effective radius associated with aerosol type
      REAL(fp)             :: RAA(5,A_)

      ! SAA: Single scattering albedo
      REAL(fp)             :: SAA(5,A_)

      ! NAA: Number of categories for scattering phase functions
      INTEGER              :: NAA

      !-----------------------------------------------------------------------
      ! Variables in file 'jv_spec_aod.dat' (RD_AOD)
      !-----------------------------------------------------------------------

      ! QAA_AOD: Aerosol scattering phase functions
      REAL(fp)             :: QAA_AOD(A_)

      ! WAA: 5 Wavelengths for the supplied phase functions
      REAL(fp)             :: WAA_AOD(A_)

      ! PAA: Phase function: first 8 terms of expansion
      REAL(fp)             :: PAA_AOD(8,A_)

      ! RAA: Effective radius associated with aerosol type
      REAL(fp)             :: RAA_AOD(A_)

      ! SAA: Single scattering albedo
      REAL(fp)             :: SAA_AOD(A_)

      !-----------------------------------------------------------------------
      ! Variables in file 'atmos_std.dat' (RD_PROF)
      !-----------------------------------------------------------------------

      ! T and O3 reference profiles
      REAL(fp), DIMENSION(51,18,12) :: TREF, OREF

      ! Interfacing indices for GC and FJX aerosols
      INTEGER, ALLOCATABLE  :: MIEDX(:)

      ! TITLEAA: Title for scattering data
      CHARACTER*80, DIMENSION(A_) :: TITLEAA

      ! Dust and aerosol optical depths
      REAL(fp), ALLOCATABLE :: ODMDUST(:,:,:,:,:)
      REAL(fp), ALLOCATABLE :: ODAER(:,:,:,:,:)
      REAL(fp), ALLOCATABLE :: ISOPOD(:,:,:,:)   ! eam, 2014

      INTEGER NJX,NW1,NW2

      !-----------------------------------------------------------------------
      !  Variables added for RRTMG (dar, mps, 12/5/14)
      !-----------------------------------------------------------------------

      INTEGER, PARAMETER :: NWVAA   = 41     !number of wavelengths in LUT
      INTEGER, PARAMETER :: NWVAA0  = 11     !number of non-RRTMG wavelengths
      INTEGER, PARAMETER :: NWVAART = NWVAA-NWVAA0 !number of RRTMG wvs
      INTEGER, PARAMETER :: NRAA    = 7      !number of aer sizes in LUT
      INTEGER, PARAMETER :: NALBD   = 2
      INTEGER, PARAMETER :: NEMISS  = 16
      INTEGER, PARAMETER :: NSPECRADMENU = 11!number of choices for RT flux

      ! Now set the following in Init_CMN_FJX below (mps, 1/3/18)
      INTEGER            :: NSPAA            !number of species in LUT
      INTEGER            :: NASPECRAD        !aerosol species in RT
      INTEGER            :: NSPECRAD         !aerosol+gas species in RT

      ! New optical arrays
      REAL*8, ALLOCATABLE :: WVAA(:,:)
      REAL*8, ALLOCATABLE :: RHAA(:,:)
      REAL*8, ALLOCATABLE :: NRLAA(:,:,:)
      REAL*8, ALLOCATABLE :: NCMAA(:,:,:)
      REAL*8, ALLOCATABLE :: RDAA(:,:)
      REAL*8, ALLOCATABLE :: RWAA(:,:)
      REAL*8, ALLOCATABLE :: SGAA(:,:)
      REAL*8, ALLOCATABLE :: QQAA(:,:,:)
      REAL*8, ALLOCATABLE :: ALPHAA(:,:,:)
      REAL*8, ALLOCATABLE :: REAA(:,:)
      REAL*8, ALLOCATABLE :: SSAA(:,:,:)
      REAL*8, ALLOCATABLE :: ASYMAA(:,:,:)
      REAL*8, ALLOCATABLE :: PHAA(:,:,:,:)
      INTEGER :: IWVSELECT(2,3) !index of requested wavelengths
      INTEGER :: IRTWVSELECT(2,3) !index of requested RT wavelengths

      ! max of 3 but need 2 per wavelength if interpolating
      INTEGER :: NWVREQUIRED !number of wvs required for interpolation
      INTEGER :: IWVREQUIRED(6) !index of wavelengths for interpo.
      INTEGER :: NRTWVREQUIRED !number of wvs required for RT interpolation
      INTEGER :: IRTWVREQUIRED(6) !index of wavelengths for RT interpo.
      ! list of required wavelengths, up to max of 3 x 2

      INTEGER :: IWV1000 !Store the wavelength index for 1000nm for Fast-J

      !coefficients for interpolation of wavelength (and for RT too)
      REAL*8      :: ACOEF_WV(3),BCOEF_WV(3),CCOEF_WV(3)
      REAL*8      :: ACOEF_RTWV(3),BCOEF_RTWV(3),CCOEF_RTWV(3)
      INTEGER     :: LSPECRADMENU(NSPECRADMENU) !indices for number of species
                                                ! on input menu 
      INTEGER, ALLOCATABLE :: SPECMASK(:) !list of binary switches for different
                                          !species flux output

      ! RH indices 
      INTEGER, ALLOCATABLE :: IRHARR(:,:,:)

#if defined( RRTMG )
      !to pass to RT code
      !one for each hydrophilic/hydrophobic aerosol and optical dust bin
      !and also sulfate, nitrate and ammonia are separate too
      REAL*8, ALLOCATABLE :: RTODAER(:,:,:,:,:)
      REAL*8, ALLOCATABLE :: RTSSAER(:,:,:,:,:)
      REAL*8, ALLOCATABLE :: RTASYMAER(:,:,:,:,:)
#endif

      CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Init_Cmn_Fjx
!
! !DESCRIPTION: Routine INIT\_CMN\_FJX initializes quantities based on
!  the grid-independent size parameters.
!\\
!\\
! !INTERFACE:

      SUBROUTINE Init_CMN_FJX( am_I_Root, Input_Opt, RC )
!
! !USES:
!
      USE ErrCode_Mod
      USE Input_Opt_Mod, ONLY : OptInput
!
! !INPUT PARAMETERS: 
! 
      LOGICAL,        INTENT(IN)  :: am_I_Root   ! Are we on the root CPU?
      TYPE(OptInput), INTENT(IN)  :: Input_Opt   ! Input Options object
!
! !OUTPUT PARAMETERS:
!
      INTEGER,        INTENT(OUT) :: RC          ! Success or failure?
! 
! !REVISION HISTORY: 
!  21 Feb 2014 - M. Sulprizio- Renamed from Set_CMN_FJX_MOD to Init_CMN_FJX.
!                              We now set several variables here that depend on
!                              LLPAR, since LLPAR is no longer a parameter and
!                              can't be used to define a parameter.
!  29 Dec 2017 - C. Keller   - Initialize RNAMES to empty strings.
!  03 Jan 2018 - M. Sulprizio- Replace UCX CPP switch with Input_Opt%LUCX
!  23 Mar 2018 - E. Lundgren - Move NWVSELECT, WVSELECT, and STRWVSELECT to 
!                              Input_Opt since read in from input.geos
!EOP
!------------------------------------------------------------------------------
!BOC

      !=================================================================
      ! INIT_CMN_FJX begins here!
      !=================================================================

      L_     = LLPAR    ! Number of CTM layers
      L1_    = L_+1     ! Number of CTM layer edges
      L2_    = L1_*2    ! Number of levels in FJX grid that
                        ! inc. both edges and mid-points
      JVL_   = LLPAR    ! Vertical levels for J-values

      JXL_   = LLPAR    ! Vertical levels for J-values computed within Fast-JX
      JXL1_  = JXL_+1   ! Vertical levels edges for J-values
      JXL2_  = 2*JXL_+2 ! Max # levels in the basic Fast-JX grid (mid-level)

      
#if defined( MODEL_GEOS )
      ! N_  = no. of levels in Mie scattering arrays
      IF ( Input_Opt%LLFASTJX > 0 ) THEN
         N_ = Input_Opt%LLFASTJX
      ELSE
         N_ = 601
      ENDIF
#endif

      JTAUMX = ( N_ - 4*JXL_ ) / 2  ! Maximum number of divisions ( i.e., may
                                    ! not get to ATAUMN)

      !-----------------------------------------------------------------------
      !  Variables that differ for UCX-based mechanisms (mps, 1/3/18)
      !-----------------------------------------------------------------------

      IF ( Input_Opt%LUCX ) THEN
         AN_          = 37  ! # of separate aerosols per layer; Including PSCs
         W_           = 18  ! # of wavelength bins 

         ! For RRTMG:
         NSPAA        = 8   ! number of species in LUT
         NASPECRAD    = 16  ! aerosol species in RT
         NSPECRAD     = 18  ! aerosol+gas species in RT
      ELSE
         AN_          = 35  ! # of separate aerosols per layer
         W_           = 12  ! # of wavelength bins 

         ! For RRTMG:
         NSPAA        = 6   ! number of species in LUT
         NASPECRAD    = 14  ! aerosol species in RT
         NSPECRAD     = 16  ! aerosol+gas species in RT
      ENDIF

      !-----------------------------------------------------------------------
      !  Allocate arrays
      !-----------------------------------------------------------------------

      ALLOCATE( ZPJ( LLPAR, JVN_, IIPAR, JJPAR ), STAT=RC )
      IF ( RC /= GC_SUCCESS ) RETURN
      ZPJ = 0e+0_fp

      ALLOCATE( ODMDUST( IIPAR, JJPAR, LLPAR, NWVAA, NDUST), STAT=RC )
      IF ( RC /= GC_SUCCESS ) RETURN
      ODMDUST = 0e+0_fp

      ALLOCATE( ODAER( IIPAR, JJPAR, LLPAR, NWVAA, NAER), STAT=RC )
      IF ( RC /= GC_SUCCESS ) RETURN
      ODAER = 0e+0_fp

      ALLOCATE( MIEDX(AN_), STAT=RC )
      IF ( RC /= GC_SUCCESS ) RETURN
      MIEDX = 0

      ! Allocate array for isoprene SOA AOD (eam, 2014):
      ALLOCATE( ISOPOD( IIPAR, JJPAR, LLPAR, NWVAA), STAT=RC )
      IF ( RC /= GC_SUCCESS ) RETURN
      ISOPOD = 0e+0_fp

      !-----------------------------------------------------------------------
      !  Variables added for RRTMG (dar, mps, 12/5/14)
      !-----------------------------------------------------------------------

      ALLOCATE( IRHARR(IIPAR,JJPAR,LLPAR), STAT=RC )
      IF ( RC /= GC_SUCCESS ) RETURN
      IRHARR = 0d0

      ALLOCATE( WVAA(NWVAA,NSPAA), STAT=RC )
      IF ( RC /= GC_SUCCESS ) RETURN
      WVAA = 0d0

      ALLOCATE( RHAA(NRAA,NSPAA), STAT=RC )
      IF ( RC /= GC_SUCCESS ) RETURN
      RHAA = 0d0

      ALLOCATE( NRLAA(NWVAA,NRAA,NSPAA), STAT=RC )
      IF ( RC /= GC_SUCCESS ) RETURN
      NRLAA = 0d0

      ALLOCATE( NCMAA(NWVAA,NRAA,NSPAA), STAT=RC )
      IF ( RC /= GC_SUCCESS ) RETURN
      NCMAA = 0d0

      ALLOCATE( RDAA(NRAA,NSPAA), STAT=RC )
      IF ( RC /= GC_SUCCESS ) RETURN
      RDAA = 0d0

      ALLOCATE( RWAA(NRAA,NSPAA), STAT=RC )
      IF ( RC /= GC_SUCCESS ) RETURN
      RWAA = 0d0

      ALLOCATE( SGAA(NRAA,NSPAA), STAT=RC )
      IF ( RC /= GC_SUCCESS ) RETURN
      SGAA = 0d0

      ALLOCATE( QQAA(NWVAA,NRAA,NSPAA), STAT=RC )
      IF ( RC /= GC_SUCCESS ) RETURN
      QQAA = 0d0

      ALLOCATE( ALPHAA(NWVAA,NRAA,NSPAA), STAT=RC )
      IF ( RC /= GC_SUCCESS ) RETURN
      ALPHAA = 0d0

      ALLOCATE( REAA(NRAA,NSPAA), STAT=RC )
      IF ( RC /= GC_SUCCESS ) RETURN
      REAA = 0d0

      ALLOCATE( SSAA(NWVAA,NRAA,NSPAA), STAT=RC )
      IF ( RC /= GC_SUCCESS ) RETURN
      SSAA = 0d0

      ALLOCATE( ASYMAA(NWVAA,NRAA,NSPAA), STAT=RC )
      IF ( RC /= GC_SUCCESS ) RETURN
      ASYMAA = 0d0

      ALLOCATE( PHAA(NWVAA,NRAA,NSPAA,8), STAT=RC )
      IF ( RC /= GC_SUCCESS ) RETURN
      PHAA = 0d0

      ALLOCATE( SPECMASK(NSPECRAD), STAT=RC )
      IF ( RC /= GC_SUCCESS ) RETURN
      SPECMASK = 0

#if defined( RRTMG )
      ! +2 to split SNA into SU, NI and AM
      ALLOCATE( RTODAER(IIPAR,JJPAR,LLPAR,NWVAA,NAER+2+NDUST), STAT=RC )
      IF ( RC /= GC_SUCCESS ) RETURN
      RTODAER = 0d0

      ALLOCATE( RTSSAER(IIPAR,JJPAR,LLPAR,NWVAA,NAER+2+NDUST), STAT=RC )
      IF ( RC /= GC_SUCCESS ) RETURN
      RTSSAER = 0d0

      ALLOCATE( RTASYMAER(IIPAR,JJPAR,LLPAR,NWVAA,NAER+2+NDUST),STAT=RC)
      IF ( RC /= GC_SUCCESS ) RETURN
      RTASYMAER = 0d0
#endif

      ! Initialize RNAMES to empty string (ckeller,12/29/17)
      RNAMES(:) = ""

      ! Return w/ success
      RC = GC_SUCCESS

      END SUBROUTINE Init_CMN_FJX
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Cleanup_Cmn_Fjx
!
! !DESCRIPTION: Subroutine CLEANUP\_CMN\_FJX deallocates all module arrays.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE Cleanup_CMN_FJX( am_I_Root, RC )
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
! !REVISION HISTORY: 
!  21 Feb 2014 - M. Sulprizio- Initial version
!  10 Aug 2018 - H.P. Lin    - Add missing allocations for RRTMG
!EOP
!------------------------------------------------------------------------------
!BOC
      !=================================================================
      ! CLEANUP_CMN_FJX begins here!
      !=================================================================
      IF ( ALLOCATED( ZPJ       ) ) DEALLOCATE( ZPJ       )
      IF ( ALLOCATED( ODMDUST   ) ) DEALLOCATE( ODMDUST   )
      IF ( ALLOCATED( ODAER     ) ) DEALLOCATE( ODAER     )
      IF ( ALLOCATED( MIEDX     ) ) DEALLOCATE( MIEDX     )
      IF ( ALLOCATED( ISOPOD    ) ) DEALLOCATE( ISOPOD    )
      IF ( ALLOCATED( IRHARR    ) ) DEALLOCATE( IRHARR    )
      IF ( ALLOCATED( WVAA      ) ) DEALLOCATE( WVAA      )
      IF ( ALLOCATED( RHAA      ) ) DEALLOCATE( RHAA      )
      IF ( ALLOCATED( NRLAA     ) ) DEALLOCATE( NRLAA     )
      IF ( ALLOCATED( NCMAA     ) ) DEALLOCATE( NCMAA     )
      IF ( ALLOCATED( RDAA      ) ) DEALLOCATE( RDAA      )
      IF ( ALLOCATED( RWAA      ) ) DEALLOCATE( RWAA      )
      IF ( ALLOCATED( SGAA      ) ) DEALLOCATE( SGAA      )
      IF ( ALLOCATED( QQAA      ) ) DEALLOCATE( QQAA      )
      IF ( ALLOCATED( ALPHAA    ) ) DEALLOCATE( ALPHAA    )
      IF ( ALLOCATED( REAA      ) ) DEALLOCATE( REAA      )
      IF ( ALLOCATED( SSAA      ) ) DEALLOCATE( SSAA      )
      IF ( ALLOCATED( ASYMAA    ) ) DEALLOCATE( ASYMAA    )
      IF ( ALLOCATED( PHAA      ) ) DEALLOCATE( PHAA      )
      IF ( ALLOCATED( SPECMASK  ) ) DEALLOCATE( SPECMASK  )


#if defined( RRTMG )
      IF ( ALLOCATED( RTODAER   ) ) DEALLOCATE( RTODAER   )
      IF ( ALLOCATED( RTSSAER   ) ) DEALLOCATE( RTSSAER   )
      IF ( ALLOCATED( RTASYMAER ) ) DEALLOCATE( RTASYMAER )
#endif

      ! Return successfully
      RC = GC_SUCCESS

      END SUBROUTINE Cleanup_CMN_FJX
!EOC
      END MODULE CMN_FJX_MOD
