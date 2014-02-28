!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: cmn_fjx_mod
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

      IMPLICIT NONE
      PUBLIC
!
! !DEFINED PARAMETERS:
!
      ! New (SDE 03/28/13)
      ! Index in RAA & QAA of 999 nm wavelength
      INTEGER, PARAMETER :: IND999  = 5

      ! Required size of aerosol arrays
      INTEGER, PARAMETER :: NALLAER = (NRHAER*NRH)+NSTRATAER

!------------------------------------------------------------------------------
! Prior to 2/21/14:
! These are no longer parameters, but are set in Init_CMN_FJ below. LLPAR is
! no longer a parameter in CMN_SIZE_mod, so it can't be used to define a
! parameter (mps, 2/21/14)
!      INTEGER, PARAMETER :: L_=LLPAR   ! Number of CTM layers
!      INTEGER, PARAMETER :: L1_=L_+1   ! Number of CTM layer edges
!      INTEGER, PARAMETER :: L2_=L1_*2  ! Number of levels in FJX grid that
!                                       ! inc. both edges and mid-points
!      INTEGER, PARAMETER :: JVL_=LLPAR ! Vertical levels for J-values
!------------------------------------------------------------------------------

      INTEGER            :: L_             ! Number of CTM layers

      INTEGER            :: L1_            ! Number of CTM layer edges

      INTEGER            :: L2_            ! Number of levels in FJX grid that
                                           ! inc. both edges and mid-points

      INTEGER            :: JVL_           ! Vertical levels for J-values

      INTEGER, PARAMETER :: JVN_ = 121     ! Max number of J-values

#if defined( UCX )
      INTEGER, PARAMETER :: AN_  = 37      ! Including PSCs
#else
      INTEGER, PARAMETER :: AN_  = 35      ! # of separate aerosols per layer
#endif

      ! Variables used to interface GEOS-Chem and Fast-JX at runtime
      ! Branches for photolysis species
      INTEGER            :: BRANCH(JVN_)

      ! Names of photolysis species
      CHARACTER (LEN=10) :: RNAMES(JVN_)

      ! Mapping array from Harvard species names to UCI species names
      INTEGER            :: RINDEX(JVN_)                

      ! Output J values
      REAL*8, ALLOCATABLE :: ZPJ(:,:,:,:)

      !-----------------------------------------------------------------------
      ! variables used to map fast-JX J's onto CTM J's
      !-----------------------------------------------------------------------

      ! Multiplication factor for fast-JX calculated J
      REAL*8             :: JFACTA(JVN_)

      ! Index arrays that map Jvalue(j) onto rates
      INTEGER            :: JIND(JVN_)

      ! Mumber of Photolysis reactions in CTM chemistry, derived here NRATJ
      ! must be .le. JVN_
      INTEGER            :: NRATJ

      ! Label of J-value used in the main chem model
      CHARACTER*50       :: JLABEL(JVN_)

!------------------------------------------------------------------------------
! Prior to 2/21/14:
! These are no longer parameters, but are set in Init_CMN_FJ below. LLPAR is
! no longer a parameter in CMN_SIZE_mod, so it can't be used to define a
! parameter (mps, 2/21/14)
!      ! JXL_: vertical(levels) dim for J-values computed within fast-JX
!!      INTEGER, PARAMETER ::  JXL_=100, JXL1_=JXL_+1
!     INTEGER, PARAMETER ::  JXL_=LLPAR, JXL1_=JXL_+1
!     ! JXL2_: 2*JXL_ + 2 = mx no. levels in the basic Fast-JX grid (mid-level)
!     INTEGER, PARAMETER ::  JXL2_=2*JXL_+2
!------------------------------------------------------------------------------
      ! JXL_: vertical(levels) dim for J-values computed within fast-JX
      INTEGER            ::  JXL_
      INTEGER            ::  JXL1_

      ! JXL2_: 2*JXL_ + 2 = mx no. levels in the basic Fast-JX grid (mid-level)
      INTEGER            ::  JXL2_

      ! WX_  = dim = no. of wavelengths in input file
      INTEGER, PARAMETER ::  WX_   = 18

      ! X_   = dim = max no. of X-section data sets (input data)
      INTEGER, PARAMETER ::  X_    = 72

      ! A_   = dim = no. of Aerosol/cloud Mie sets (input data)
      INTEGER, PARAMETER ::  A_    = 56

      ! W_   = dim = no. of Wavelength bins:  =18 std, =12 trop only
#if defined( UCX )
      INTEGER, PARAMETER ::  W_    = 18
#else
      INTEGER, PARAMETER ::  W_    = 12
#endif

      ! N_  = no. of levels in Mie scattering arrays
      !     = 2*NC+1 = 4*(L_+1) + 1 + 2*sum(JADDLV)
      INTEGER, PARAMETER ::  N_    = 601

      ! M_  = no. of Gauss points used, must = 4 in fast_JX (no option)
      INTEGER, PARAMETER ::  M_    = 4

      ! M2_ = 2*M_ = 8, replaces MFIT
      INTEGER, PARAMETER ::  M2_   = 2*M_

      !-----------------------------------------------------------------------
      ! 4 Gauss pts = 8-stream
      REAL*8, DIMENSION(M_), PARAMETER  ::
     &                   EMU = [.06943184420297d0, .33000947820757d0,
     &                          .66999052179243d0, .93056815579703d0]
      REAL*8, DIMENSION(M_), PARAMETER  ::
     &                   WT  = [.17392742256873d0, .32607257743127d0,
     &                          .32607257743127d0, .17392742256873d0]
      !-----------------------------------------------------------------------

      ! ZZHT: scale height (cm)
      REAL*8, PARAMETER  :: ZZHT   = 5.d5

      ! RAD: Radius of Earth (cm)
      REAL*8, PARAMETER  :: RAD    = 6375.d5

      ! ATAU: heating rate (factor increase from one layer to the next)
      REAL*8, PARAMETER  :: ATAU   = 1.120d0

      ! ATAU0: minimum heating rate
      REAL*8, PARAMETER  :: ATAU0  = 0.010d0

!------------------------------------------------------------------------------
! Prior to 2/21/14:
! These are no longer parameters, but are set in Init_CMN_FJ below. JXL is no
! longer a parameter, so it can't be used to define a parameter (mps, 2/21/14)
!      ! JTAUMX = maximum number of divisions (i.e., may not get to ATAUMN)
!      INTEGER, PARAMETER  :: JTAUMX = (N_ - 4*JXL_)/2
!------------------------------------------------------------------------------
      ! JTAUMX = maximum number of divisions (i.e., may not get to ATAUMN)
      INTEGER            :: JTAUMX

      ! Physical constants
      REAL*8, PARAMETER  :: UVXPLANCK = 6.62606957e-34
      REAL*8, PARAMETER  :: UVXCCONST = 2.99792458e8

      ! Conversion factors from photons/cm2s to W/m2
      REAL*8, DIMENSION(WX_) :: UVXFACTOR

      !-----------------------------------------------------------------------
      ! Variables in file 'FJX_spec.dat' (RD_XXX)
      !-----------------------------------------------------------------------

      ! WBIN: Boundaries of wavelength bins
      REAL*8             :: WBIN(WX_+1)

      ! WL: Centres of wavelength bins - 'effective wavelength'
      REAL*8             :: WL(WX_)

      ! FL: Solar flux incident on top of atmosphere (cm-2.s-1)
      REAL*8             :: FL(WX_)

      REAL*8             :: QO2(WX_,3)   ! QO2: O2 cross-sections
      REAL*8             :: QO3(WX_,3)   ! QO3: O3 cross-sections
      REAL*8             :: Q1D(WX_,3)   ! Q1D: O3 => O(1D) quantum yield

      ! QQQ: Supplied cross sections in each wavelength bin (cm2)
      REAL*8             :: QQQ(WX_,3,X_)

      ! QRAYL: Rayleigh parameters (effective cross-section) (cm2)
      REAL*8             :: QRAYL(WX_+1)

      ! TQQ: Temperature for supplied cross sections
      REAL*8             :: TQQ(3,X_)

      ! LQQ = 1, 2, or 3 to determine interpolation with T or P
      INTEGER            :: LQQ(X_)

      ! TITLEJX: Title for supplied cross sections, from 'FJX_spec.dat'
      CHARACTER*6        :: TITLEJX(X_)

      ! SQQ: Flag for supplied cross sections, from 'FJX_spec.dat'
      CHARACTER*1        :: SQQ(X_)

      !-----------------------------------------------------------------------
      ! Variables in file 'jv_spec_mie.dat' (RD_MIE)
      !-----------------------------------------------------------------------

      ! QAA: Aerosol scattering phase functions
      REAL*8             :: QAA(5,A_)

      ! WAA: 5 Wavelengths for the supplied phase functions
      REAL*8             :: WAA(5,A_)

      ! PAA: Phase function: first 8 terms of expansion
      REAL*8             :: PAA(8,5,A_)

      ! RAA: Effective radius associated with aerosol type
      REAL*8             :: RAA(5,A_)

      ! SAA: Single scattering albedo
      REAL*8             :: SAA(5,A_)

      ! NAA: Number of categories for scattering phase functions
      INTEGER            :: NAA

      !-----------------------------------------------------------------------
      ! Variables in file 'jv_spec_aod.dat' (RD_AOD)
      !-----------------------------------------------------------------------

      ! QAA_AOD: Aerosol scattering phase functions
      REAL*8             :: QAA_AOD(A_)

      ! WAA: 5 Wavelengths for the supplied phase functions
      REAL*8             :: WAA_AOD(A_)

      ! PAA: Phase function: first 8 terms of expansion
      REAL*8             :: PAA_AOD(8,A_)

      ! RAA: Effective radius associated with aerosol type
      REAL*8             :: RAA_AOD(A_)

      ! SAA: Single scattering albedo
      REAL*8             :: SAA_AOD(A_)

      !-----------------------------------------------------------------------
      ! Variables in file 'atmos_std.dat' (RD_PROF)
      !-----------------------------------------------------------------------

      ! T and O3 reference profiles
      REAL*8, DIMENSION(51,18,12) :: TREF, OREF

      ! Interfacing indices for GC and FJX aerosols
      INTEGER, DIMENSION(AN_)     :: MIEDX

      ! TITLEAA: Title for scattering data
      CHARACTER*20, DIMENSION(A_) :: TITLEAA

      ! Dust and aerosol optical depths
      REAL*8, ALLOCATABLE :: ODMDUST(:,:,:,:)
      REAL*8, ALLOCATABLE :: ODAER(:,:,:,:)

      INTEGER NJX,NW1,NW2

      CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_cmn_fjx
!
! !DESCRIPTION: Routine INIT\_CMN\_FJX initializes quantities based on
!  the grid-independent size parameters.
!\\
!\\
! !INTERFACE:

      SUBROUTINE Init_CMN_FJX( am_I_Root, RC )
!
! !USES:
!
      USE GIGC_ErrCode_Mod
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
!  21 Feb 2014 - M. Sulprizio- Renamed from Set_CMN_FJX_MOD to Init_CMN_FJX.
!                              We now set several variables here that depend on
!                              LLPAR, since LLPAR is no longer a parameter and
!                              can't be used to define a parameter.
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

      JTAUMX = ( N_ - 4*JXL_ ) / 2  ! Maximum number of divisions ( i.e., may
                                    ! not get to ATAUMN)

      ALLOCATE( ZPJ( LLPAR, JVN_, IIPAR, JJPAR ), STAT=RC )
      IF ( RC /= GIGC_SUCCESS ) RETURN
      ZPJ = 0d0

      ALLOCATE( ODMDUST( IIPAR, JJPAR, LLPAR, NDUST), STAT=RC )
      IF ( RC /= GIGC_SUCCESS ) RETURN
      ODMDUST = 0d0

      ALLOCATE( ODAER( IIPAR, JJPAR, LLPAR, NALLAER), STAT=RC )
      IF ( RC /= GIGC_SUCCESS ) RETURN
      ODAER = 0d0

      ! Return w/ success
      RC = GIGC_SUCCESS

      END SUBROUTINE Init_CMN_FJX
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: cleanup_cmn_fjx
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
      USE GIGC_ErrCode_Mod
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
!EOP
!------------------------------------------------------------------------------
!BOC
      !=================================================================
      ! CLEANUP_CMN_FJX begins here!
      !=================================================================
      IF ( ALLOCATED( ZPJ    ) ) DEALLOCATE( ZPJ    )

      ! Return successfully
      RC = GIGC_SUCCESS

      END SUBROUTINE Cleanup_CMN_FJX
!EOC
      END MODULE CMN_FJX_MOD
