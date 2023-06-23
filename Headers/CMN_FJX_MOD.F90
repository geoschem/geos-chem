!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: CMN_FJX_mod.F90
!
! !DESCRIPTION: Module CMN\_FJX\_MOD contains parameters and global variables
!  used to interface between Harvard chemistry and UC-Irvine photolysis
!  programs (Fast-J/Fast-JX). Contents of this module previously were in
!  Headers/CMN_FJX_MOD.F90. That module was split into this file and
!  Headers/CMN_Phot_mod.F90 for the development of Cloud-J to replace
!  Fast-JX. This module will be used only in Fast-JX and not in Cloud-J.
!\\
!\\
! !INTERFACE:
!
MODULE CMN_FJX_MOD
!
! !USES:
!
  USE PRECISION_MOD      ! For GEOS-Chem Precision (fp)

  IMPLICIT NONE
  PUBLIC
!
! !DEFINED PARAMETERS:
!

  !-----------------------------------------------------------------------
  ! Parameters
  !-----------------------------------------------------------------------

  INTEGER, PARAMETER :: JVN_ = 166     ! Max number of J-values

  INTEGER, PARAMETER :: WX_   = 18     ! # wavelengths in file FJX_spec.dat

  INTEGER, PARAMETER :: X_    = 123    ! Max # X-section data in FJX_spec.dat

  INTEGER, PARAMETER :: A_    = 56     ! # aerosol/cloud Mie sets

#ifdef MODEL_GEOS
!!!INTEGER, PARAMETER :: N_ = 601
!!!INTEGER, PARAMETER :: N_ = 1201
  INTEGER            :: N_
#else
  INTEGER, PARAMETER :: N_    = 601    ! # levels in Mie scattering arrays
                                       ! = 2*NC+1 = 4*(L_+1) + 1 + 2*sum(JADDLV)
#endif

  INTEGER, PARAMETER :: M_    = 4      ! # Gauss points used (must be 4 in FJX)

  INTEGER, PARAMETER :: M2_   = 2*M_   ! M2_ = 2*M_ = 8, replaces MFIT

  ! 4 Gauss pts = 8-stream
  REAL(fp), DIMENSION(M_), PARAMETER  :: &
       EMU = [.06943184420297e+0_fp, .33000947820757e+0_fp, &
              .66999052179243e+0_fp, .93056815579703e+0_fp]
  REAL(fp), DIMENSION(M_), PARAMETER  :: &
       WT  = [.17392742256873e+0_fp, .32607257743127e+0_fp, &
              .32607257743127e+0_fp, .17392742256873e+0_fp]

  ! Physical constants
  REAL(fp), PARAMETER  :: ZZHT   = 5.e+5_fp      ! ZZHT: scale height (cm)

  REAL(fp), PARAMETER  :: RAD    = 6375.e+5_fp   ! RAD: Radius of Earth (cm)

  REAL(fp), PARAMETER  :: ATAU   = 1.120e+0_fp   ! ATAU: heating rate (factor
                                                 ! increase layer to layer
#ifdef MODEL_GEOS
  !REAL(fp), PARAMETER  :: ATAU   = 1.180e+0_fp
#endif

  REAL(fp), PARAMETER  :: ATAU0  = 0.010e+0_fp   ! ATAU0: minimum heating rate

  !-----------------------------------------------------------------------
  ! Array dimension sizes set in Init_CMN_FJX at run-time
  !-----------------------------------------------------------------------

  INTEGER :: L_             ! Number of CTM layers

  INTEGER :: L1_            ! Number of CTM layer edges

  INTEGER :: L2_            ! Number of levels in FJX grid that
                            ! inc. both edges and mid-points

  INTEGER :: JVL_           ! Vertical levels for J-values

  INTEGER :: AN_            ! # of separate aerosols per layer

  INTEGER :: JXL_           ! Vertical level for Jvals (mid)

  INTEGER :: JXL1_          ! Vertical level for Jvals (edge)

  INTEGER :: JXL2_          ! Max # levels in Fast-JX grid (mid)

  INTEGER :: W_             ! # wavelength bins

  INTEGER :: JTAUMX         ! Max # divisions, i.e. < ATAUMN

  !-----------------------------------------------------------------------
  ! RD_XXX variables (file 'FJX_spec.dat')
  !-----------------------------------------------------------------------

  INTEGER              :: NW1
  INTEGER              :: NW2

  ! WL: Centres of wavelength bins - 'effective wavelength'
  REAL(fp)             :: WL(WX_)

  ! FL: Solar flux incident on top of atmosphere (cm-2.s-1)
  REAL(fp)             :: FL(WX_)

  ! QRAYL: Rayleigh parameters (effective cross-section) (cm2)
  REAL(fp)             :: QRAYL(WX_+1)

  ! TITLEJX: Title for supplied cross sections, from 'FJX_spec.dat'
  CHARACTER*6          :: TITLEJX(X_)

  ! LQQQ = 1, 2, or 3 to determine interpolation with T or P
  INTEGER              :: LQQ(X_)

  ! SQQ: Flag for supplied cross sections, from 'FJX_spec.dat'
  CHARACTER*1          :: SQQ(X_)

  ! TQQ: Temperature for supplied cross sections
  REAL(fp)             :: TQQ(3,X_)

  ! QQQ: Supplied cross sections in each wavelength bin (cm2)
  REAL(fp)             :: QQQ(WX_,3,X_)

  ! NJX: Number of species to calculate J-values for
  INTEGER              :: NJX

  ! O2 and O3 cross-sections, and O3 => O(1D) quantum yield
  REAL(fp)             :: QO2(WX_,3)
  REAL(fp)             :: QO3(WX_,3)
  REAL(fp)             :: Q1D(WX_,3)

  ! WBIN: Boundaries of wavelength bins
  REAL(fp)             :: WBIN(WX_+1)

  !-----------------------------------------------------------------------
  ! RD_MIE variables (file 'jv_spec_mie.dat')
  !-----------------------------------------------------------------------

  ! NAA: Number of categories for scattering phase functions
  INTEGER              :: NAA

  ! TITLAA: Title per scattering data set, e.g. "01 RAYLE  = Rayleigh phase"
  CHARACTER*80, DIMENSION(A_) :: TITLAA

  ! WAA: 5 Wavelengths for the supplied phase functions. 1st col in file.
  REAL(fp)             :: WAA(5,A_)

  ! QAA: Aerosol scattering phase functions. 2nd col in file.
  REAL(fp)             :: QAA(5,A_)

  ! RAA: Effective radius associated with aerosol type. 3rd col in file.
  REAL(fp)             :: RAA(5,A_)

  ! SAA: Single scattering albedo. 4th col in file.
  REAL(fp)             :: SAA(5,A_)

  ! PAA: Phase function: first 8 terms of expansion. Cols 5:12 in file.
  REAL(fp)             :: PAA(8,5,A_)

  !-----------------------------------------------------------------------
  ! RD_JS_JX variables (file 'FJX_j2j.dat')
  !-----------------------------------------------------------------------

  ! Label of J-value used in the main chem model
  CHARACTER*50       :: JLABEL(JVN_)

  ! Multiplication factor for fast-JX calculated J
  REAL(fp)           :: JFACTA(JVN_)

  ! Mumber of Photolysis reactions in CTM chemistry, derived here NRATJ
  ! must be .le. JVN_
  INTEGER            :: NRATJ

  ! Names of photolysis species
  CHARACTER (LEN=10) :: RNAMES(JVN_)

  ! Branches for photolysis species
  INTEGER            :: BRANCH(JVN_)

  ! Index arrays that map Jvalue(j) onto rates
  INTEGER            :: JIND(JVN_)

CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Init_Cmn_FJX
!
! !DESCRIPTION: Routine INIT\_CMN\_FJX initializes quantities based on
!  the grid-independent size parameters.
!\\
!\\
! !INTERFACE:

  SUBROUTINE Init_CMN_FJX( Input_Opt, State_Grid, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE Input_Opt_Mod,  ONLY : OptInput
    USE State_Grid_Mod, ONLY : GrdState
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN)  :: Input_Opt   ! Input Options object
    TYPE(GrdState), INTENT(IN)  :: State_Grid  ! Grid State object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT) :: RC          ! Success or failure?
!
! !REVISION HISTORY:
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC

    !=================================================================
    ! INIT_CMN_FJX begins here!
    !=================================================================

    L_     = State_Grid%NZ ! Number of CTM layers
    L1_    = L_+1          ! Number of CTM layer edges
    L2_    = L1_*2         ! Number of levels in FJX grid that
                           ! inc. both edges and mid-points
    JVL_   = State_Grid%NZ ! Vertical levs for J-values

    JXL_   = State_Grid%NZ ! Vertical levs for J-values computed w/in Fast-JX
    JXL1_  = JXL_+1        ! Vertical levs edges for J-values
    JXL2_  = 2*JXL_+2      ! Max # levs in the basic Fast-JX grid (mid-level)

#ifdef MODEL_GEOS
    ! N_  = no. of levels in Mie scattering arrays
    IF ( Input_Opt%LLFASTJX > 0 ) THEN
       N_ = Input_Opt%LLFASTJX
    ELSE
       N_ = 601
    ENDIF
#endif

    JTAUMX = ( N_ - 4*JXL_ ) / 2  ! Maximum number of divisions ( i.e., may
                                  ! not get to ATAUMN)

    AN_       = 37  ! # of separate aerosols per layer; Including PSCs
    W_        = 18  ! # of wavelength bins

    ! Initialize RNAMES to empty string (ckeller,12/29/17)
    RNAMES(:) = ""

    ! Return w/ success
    RC = GC_SUCCESS

  END SUBROUTINE Init_CMN_FJX
!EOC

END MODULE CMN_FJX_MOD
