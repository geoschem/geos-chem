!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: fjx_mod.F90
!
! !DESCRIPTION: Module FJX\_MOD contains routines for the Fast-JX scheme
!  (Prather et al). Current implementation is version 7.0a. Content in this
!  file used to be in fast_jx_mod.F90 but was moved here for development of
!  Cloud-J which will replace Fast-JX in GEOS-Chem. This module interfaces
!  with fjx_interface_mod.
!\\
!\\
! !INTERFACE:
!
MODULE FJX_MOD
!
! !USES:
!
  USE PRECISION_MOD    ! For GEOS-Chem Precision (fp)
#if defined( MODEL_CESM )
  USE cam_abortutils, ONLY : endrun
  USE spmd_utils,     ONLY : mpicom, masterprocid, mpi_success
  USE spmd_utils,     ONLY : mpi_character, mpi_integer, mpi_real8
#endif

  IMPLICIT NONE

  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC :: PHOTO_JX  ! Computes J-values  
  PUBLIC :: SOLAR_JX  ! Computes solar zenith angle and solar function
  PUBLIC :: RD_MIE    ! Called in init_fjx
  PUBLIC :: RD_XXX    ! Called in init_fjx
  PUBLIC :: RD_JS_JX  ! Called in init_fjx
!
! !PRIVATE MEMBER FUNCTIONS:
!
  PRIVATE :: OPMIE     ! Called in Photo_JX
  PRIVATE :: SPHERE2   ! Called in Photo_JX
  PRIVATE :: EXTRAL    ! Called in Photo_JX
  PRIVATE :: JRATET    ! Called in Photo_JX
  PRIVATE :: X_INTERP  ! Called in both JRATET and PHOTO_JX
  PRIVATE :: MIESCT    ! Called in OPMIE
  PRIVATE :: BLKSLV    ! Called in MIESCT
  PRIVATE :: GEN_ID    ! Called in BLKSLV
  PRIVATE :: LEGND0    ! Called in MIESCT
  PRIVATE :: EXITC
!
! !REVISION HISTORY:
!  14 Dec 2022 - E. Lundgren - Adapted from fast_jx_mod
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC

CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: photo_jx
!
! !DESCRIPTION: Subroutine PHOTO\_JX is the core subroutine of Fast-JX.
!    calc J's for a single column atmosphere (aka Indep Colm Atmos or ICA)
!    needs P, T, O3, clds, aersls; adds top-of-atmos layer from climatology
!    needs day-fo-year for sun distance, SZA (not lat or long)
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE PHOTO_JX( amIRoot,   dryrun,                             &
                       U0,        REFLB,      SZA,        SOLF,       &
                       P_COL,     T_COL,      AOD999,     ILON,       &
                       ILAT,      AERX_COL,   T_CLIM,     OOJ,        &
                       ZZJ,       DDJ,        maxChemLev, State_Chm,  &
                       VALJXX,    FSBOT,      FJBOT,      FLXD,       &
                       FJFLX,     Input_Opt,  State_Diag )
!
! !USES:
!
    USE CMN_FJX_Mod,    ONLY : L_, L1_, A_, N_, W_, X_, AN_
    USE CMN_FJX_Mod,    ONLY : JXL_, JXL1_, JXL2_, JVN_
    USE CMN_FJX_Mod,    ONLY : QO2, QO3, NJX, FL, WL, QRAYL
    USE CMN_FJX_Mod,    ONLY : LQQ, TQQ, QAA, PAA, SAA
    USE CMN_SIZE_Mod,   ONLY : NRH, NRHAER
    USE Input_Opt_Mod,  ONLY : OptInput
    USE State_Chm_Mod,  ONLY : ChmState ! For NSPAA, QQAA, SSAA, PHAA
    USE State_Diag_Mod, ONLY : DgnState

    IMPLICIT NONE
!
! !INPUT PARAMETERS:
!
    LOGICAL,  INTENT(IN)                       :: amIRoot
    LOGICAL,  INTENT(IN)                       :: dryrun
    REAL(fp), INTENT(IN)                       :: U0,REFLB
    REAL(fp), INTENT(IN)                       :: SZA,SOLF
    REAL(fp), INTENT(IN), DIMENSION(L1_+1   )  :: P_COL
    REAL(fp), INTENT(IN), DIMENSION(L1_     )  :: T_COL
    LOGICAL,  INTENT(IN)                       :: AOD999
    INTEGER,  INTENT(IN)                       :: ILON, ILAT
    REAL(fp), INTENT(IN), DIMENSION(A_,L1_)    :: AERX_COL ! Aerosol column
    REAL(fp), INTENT(IN), DIMENSION(L1_   )    :: T_CLIM   ! Clim. temps (K)
    REAL(fp), INTENT(IN), DIMENSION(L1_   )    :: OOJ      ! O3 col depth (#/cm2)
    REAL(fp), INTENT(IN), DIMENSION(L1_+1 )    :: ZZJ      ! Edge alts (cm)
    REAL(fp), INTENT(IN), DIMENSION(L1_   )    :: DDJ
    INTEGER,  INTENT(IN)                       :: maxChemLev
    TYPE(ChmState), INTENT(IN)                 :: State_Chm
    TYPE(OptInput), INTENT(IN)                 :: Input_Opt
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(DgnState), INTENT(INOUT)              :: State_Diag
!
! OUTPUT PARAMETERS:
!
    REAL(fp), INTENT(OUT), DIMENSION(L_,JVN_ ) :: VALJXX
    REAL(fp), INTENT(OUT), DIMENSION(W_      ) :: FSBOT
    REAL(fp), INTENT(OUT), DIMENSION(W_      ) :: FJBOT
    REAL(fp), INTENT(OUT), DIMENSION(JXL1_,W_) :: FLXD
    REAL(fp), INTENT(OUT), DIMENSION(JXL_,W_ ) :: FJFLX
!
! !REMARKS:
!
! !REVISION HISTORY:
!  08 Dec 2022 - E. Lundgren - Adapted from S.D.Eastham's adapted Fast-JX v7.0
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! --------------------------------------------------------------------
    ! key LOCAL atmospheric data needed to solve plane-parallel J----
    ! --these are dimensioned JXL_, and must have JXL_ .ge. State_Grid%NZ
    real(fp), dimension(JXL1_+1)    :: PPJ
    integer,dimension(JXL2_+1)      :: JXTRA
    real(fp), dimension(W_)         :: FJTOP,FLXD0,RFL
    real(fp), dimension(JXL_, W_)   :: AVGF
    real(fp), dimension(JXL1_,W_)   :: DTAUX
    real(fp), dimension(8,JXL1_,W_) :: POMEGAX

    ! flux/heating arrays (along with FJFLX,FLXD,FLXD0)
    real(fp) :: ODABS,ODRAY
    real(fp) :: RFLECT
    real(fp) :: AMF2(2*JXL1_+1,2*JXL1_+1)

    ! Pointers
    REAL*8, POINTER :: QQAA(:,:,:)
    REAL*8, POINTER :: SSAA(:,:,:)
    REAL*8, POINTER :: PHAA(:,:,:,:)

    ! ---------key SCATTERING arrays for clouds+aerosols------------------
    real(fp) :: OD(5,JXL1_),SSA(5,JXL1_),SLEG(8,5,JXL1_)
    real(fp) :: OD600(JXL1_)

    ! ---------key arrays AFTER solving for J's---------------------------
    real(fp) :: FFF(W_,JXL_),VALJ(X_)
    real(fp) :: VALJL(JXL_,X_) !2-D array of J_s returned by JRATET

    integer  :: L2EDGE, I,J,K,L,M,KMIE,IDXAER,IM,LU
    INTEGER  :: KMIE2, IR
    real(fp) :: XQO3,XQO2,WAVE, TTTX
    ! --------------------------------------------------------------------
    ! For compatibility with GEOS-Chem (SDE 03/30/13)
    REAL(fp) :: QSCALING,LOCALOD,LOCALSSA

    ! T_INPUT: Input temperature (K) with extra layer for compatibility
    REAL(fp) :: T_INPUT(JXL1_+1)

    !Maps the new LUT optics wavelengths on to
    !the 5 jv_spec_mie.dat wavelengths
    ! N.B. currently 200nm and 300nm data is the same in
    ! jv_spec_mie.dat, so we copy from new LUT 300nm for both
    INTEGER :: LUTIDX(5)
    LUTIDX = (/1,1,2,6,10/)

    !=================================================================
    ! PHOTO_JX begins here!
    !=================================================================

    ! Initialize
    L2EDGE = L_ + L_ + 2
    FFF(:,:) = 0.e+0_fp

    ! Set pointers
    QQAA  => State_Chm%Phot%QQAA
    SSAA  => State_Chm%Phot%SSAA
    PHAA  => State_Chm%Phot%PHAA

    ! Fill out PPJ and TTJ with CTM data to replace fixed climatology
    DO L=1,L1_
       PPJ(L) = P_COL(L)
       T_INPUT(L) = T_COL(L)
    ENDDO

    ! Ensure TOA pressure is zero
    PPJ(L1_+1) = 0.e+0_fp
    T_INPUT(L1_+1) = T_CLIM(L1_)

    ! calculate spherical weighting functions (AMF: Air Mass Factor)
    RFLECT = REFLB

    ! --------------------------------------------------------------------
    call SPHERE2 (U0,ZZJ,AMF2,L1_,JXL1_)
    ! --------------------------------------------------------------------

    !-----------------------------------------------------------------------
    ! Modification for GEOS-Chem: Optical depths are calculated at a single
    ! wavelength for GEOS-Chem, so we perform scaling in this routine.
    ! elsewhere.
    ! (SDE 03/31/13)
    !-----------------------------------------------------------------------
    ! calculate the optical properties (opt-depth, single-scat-alb,
    ! phase-fn(1:8)) at the 5 std wavelengths 200-300-400-600-999 nm
    ! for cloud+aerosols
    do L = 1,L1_
       OD600(L) = 0.e+0_fp
       ! ODs stored with fine-grain (DTAUX) and coarse (OD)
       do K=1,W_
          DTAUX(L,K) = 0.e+0_fp
       enddo
       do K=1,5
          OD(K,L)  = 0.e+0_fp
          SSA(K,L) = 0.e+0_fp
          do I=1,8
             SLEG(I,K,L) = 0.e+0_fp
          enddo
       enddo
    enddo

    ! Clunky fix to accomodate RRTMG and UCX (DAR 01/2015)
    ! Using a combination of old optics LUT format and
    ! new optics format (greater spectral resolution and range
    ! for RRTMG). Clouds, non-species aerosols and strat aerosols
    ! are not incorporated into the new LUT so must still use the
    ! old LUT for these.

    ! Don't bother on extraneous level - leave as zero
    do L = 1,L_
       DO KMIE=1,5
          ! Clouds and non-species aerosols
          DO M=1,3
             IF (AERX_COL(M,L).gt.0e+0_fp) THEN
                IDXAER=State_Chm%Phot%MIEDX(M)
                ! Cloud (600 nm scaling)
                QSCALING = QAA(KMIE,IDXAER)/QAA(4,IDXAER)
                LOCALOD = QSCALING*AERX_COL(M,L)
                LOCALSSA = SAA(KMIE,IDXAER)*LOCALOD
                OD(KMIE,L) = OD(KMIE,L) + LOCALOD
                SSA(KMIE,L)= SSA(KMIE,L) + LOCALSSA
                DO I=1,8
                   SLEG(I,KMIE,L) = SLEG(I,KMIE,L) + &
                                    (PAA(I,KMIE,IDXAER)*LOCALSSA)
                ENDDO ! I (Phase function)
             ENDIF
          ENDDO ! M (Aerosol)
          !transpose wavelength indices from the mie LUT
          !to the new speciated LUT
          KMIE2=LUTIDX(KMIE)

          ! Stratospheric aerosols
          IM=10+(NRHAER*NRH)+1
          DO M=IM,IM+1
             IDXAER=M-IM+6 !6-STS, 7-NAT

             IF (AERX_COL(M,L).gt.0d0) THEN
                IF (AOD999) THEN
                   ! Aerosol/dust (999 nm scaling)
                   ! Fixed to dry radius
                   QSCALING = QQAA(KMIE2,1,IDXAER)/QQAA(10,1,IDXAER)
                ELSE
                   ! Aerosol/dust (550 nm scaling)
                   QSCALING = QQAA(KMIE2,1,IDXAER)/QQAA(5,1,IDXAER)
                ENDIF
                LOCALOD    = QSCALING*AERX_COL(M,L)
                LOCALSSA   = SSAA(KMIE2,1,IDXAER)*LOCALOD
                OD(KMIE,L) = OD(KMIE,L) + LOCALOD
                SSA(KMIE,L)= SSA(KMIE,L) + LOCALSSA
                DO I=1,8
                   SLEG(I,KMIE,L) = SLEG(I,KMIE,L) + &
                                    (PAA(I,KMIE,IDXAER)*LOCALSSA)
                ENDDO     ! I (Phase function)
             ENDIF

          ENDDO           ! M (Aerosol)

          ! Mineral dust (from new optics LUT)
          DO M=4,10
             IF (AERX_COL(M,L).gt.0d0) THEN
                IDXAER=State_Chm%Phot%NSPAA !dust is last in LUT
                IR=M-3
                IF (AOD999) THEN
                   QSCALING = QQAA(KMIE2,IR,IDXAER)/ &
                              QQAA(10,IR,IDXAER) !1000nm in new .dat
                ELSE
                   ! Aerosol/dust (550 nm scaling)
                   QSCALING = QQAA(KMIE2,IR,IDXAER)/ &
                              QQAA(5,IR,IDXAER)  !550nm in new .dat
                ENDIF
                LOCALOD = QSCALING*AERX_COL(M,L)
                LOCALSSA = SSAA(KMIE2,IR,IDXAER)*LOCALOD
                OD(KMIE,L) = OD(KMIE,L) + LOCALOD
                SSA(KMIE,L)= SSA(KMIE,L) + LOCALSSA
                DO I=1,8
                   SLEG(I,KMIE,L) = SLEG(I,KMIE,L) + &
                                    (PHAA(KMIE2,IR,IDXAER,I)*LOCALSSA)
                ENDDO ! I (Phase function)
             ENDIF
          ENDDO ! M (Aerosol)

          ! Other aerosol (from new optics LUT)
          DO M=1,5
             DO IR=1,5
                IDXAER=10+(M-1)*NRH+IR
                IF (AERX_COL(IDXAER,L).gt.0d0) THEN
                   IF (AOD999) THEN
                      QSCALING = QQAA(KMIE2,IR,M)/ &
                                 QQAA(10,IR,M) !1000nm in new .dat
                   ELSE
                      ! Aerosol/dust (550 nm scaling)
                      QSCALING = QQAA(KMIE2,IR,M)/ &
                                 QQAA(5,IR,M)  !550nm in new .dat
                   ENDIF
                   LOCALOD = QSCALING*AERX_COL(IDXAER,L)
                   LOCALSSA = SSAA(KMIE2,IR,M)*LOCALOD
                   OD(KMIE,L) = OD(KMIE,L) + LOCALOD
                   SSA(KMIE,L)= SSA(KMIE,L) + LOCALSSA
                   DO I=1,8
                      SLEG(I,KMIE,L) = SLEG(I,KMIE,L) + &
                                       (PHAA(KMIE2,IR,M,I)*LOCALSSA)
                   ENDDO ! I (Phase function)
                ENDIF
             ENDDO    ! IR (RH bins)
          ENDDO ! M (Aerosol)
       ENDDO ! KMIE (Mie scattering wavelength bin)

       ! Normalize
       DO KMIE=1,5
          IF (OD(KMIE,L).gt.0.e+0_fp) THEN
             SSA(KMIE,L) = SSA(KMIE,L)/OD(KMIE,L)
             DO I=1,8
                SLEG(I,KMIE,L) = SLEG(I,KMIE,L)/OD(KMIE,L)
             ENDDO
          ENDIF
       ENDDO
       ! Retrieve 600 nm OD to determine added layers
       OD600(L) = OD(4,L)
    ENDDO ! L (Layer)

    ! when combining with Rayleigh and O2-O3 abs, remember the SSA and
    !  phase fn SLEG are weighted by OD and OD*SSA, respectively.
    ! Given the aerosol+cloud OD/layer in visible (600 nm) calculate how to add
    !  additonal levels at top of clouds (now uses log spacing)
    ! --------------------------------------------------------------------
    call EXTRAL(Input_Opt,State_Diag,OD600,L1_,L2EDGE,N_,JXTRA,ILON,ILAT)
    ! --------------------------------------------------------------------

    ! set surface reflectance
    RFL(:) = max(0.e+0_fp,min(1.e+0_fp,RFLECT))

    ! --------------------------------------------------------------------
    ! Loop over all wavelength bins to calc mean actinic flux AVGF(L)
    ! --------------------------------------------------------------------
    do K = 1,W_

       WAVE = WL(K)
       ! Pick nearest Mie wavelength to get scattering properites------------
       KMIE=1                             ! use 200 nm prop for <255 nm
       if( WAVE .gt. 255.e+0_fp ) KMIE=2  ! use 300 nm prop for 255-355 nm
       if( WAVE .gt. 355.e+0_fp ) KMIE=3  ! use 400 nm prop for 355-500 nm
       if( WAVE .gt. 500.e+0_fp ) KMIE=4
       if( WAVE .gt. 800.e+0_fp ) KMIE=5

       ! Combine: Rayleigh scatters & O2 & O3 absorbers to get optical
       ! properties values at L1_=L_+1 are a pseudo/climatol layer above
       ! the top CTM layer (L_)
       do L = 1,L1_
          TTTX     = T_CLIM(L) ! Following GEOS-Chem v9-1-3
          call X_interp (TTTX,XQO2, TQQ(1,1),QO2(K,1), TQQ(2,1),QO2(K,2), &
                         TQQ(3,1),QO2(K,3), LQQ(1))
          call X_interp (TTTX,XQO3, TQQ(1,2),QO3(K,1), TQQ(2,2),QO3(K,2), &
                         TQQ(3,2),QO3(K,3), LQQ(2))
          ODABS = XQO3*OOJ(L) + XQO2*DDJ(L)*0.20948e+0_fp
          ODRAY = DDJ(L)*QRAYL(K)

          DTAUX(L,K) = OD(KMIE,L) + ODABS + ODRAY

          ! Aerosols + clouds + O2 + O3
          do I=1,8
             POMEGAX(I,L,K) = SLEG(I,KMIE,L)*OD(KMIE,L)
          enddo
          ! Add Rayleigh scattering effects
          ! Only non-zero for 1st and 3rd phase functions
          POMEGAX(1,L,K) = POMEGAX(1,L,K) + 1.0e+0_fp*ODRAY
          POMEGAX(3,L,K) = POMEGAX(3,L,K) + 0.5e+0_fp*ODRAY
          ! Normalize
          do I=1,8
             POMEGAX(I,L,K) = POMEGAX(I,L,K)/DTAUX(L,K)
          enddo
       enddo
    enddo
    ! --------------------------------------------------------------------

    LU = L_
    call OPMIE(DTAUX,POMEGAX,U0,RFL,AMF2,JXTRA, &
               AVGF,FJTOP,FJBOT,FSBOT,FJFLX,FLXD,FLXD0,LU)

    !! --------------------------------------------------------------------

    do K = 1,W_
       do L = 1,L_
          FFF(K,L) = FFF(K,L) + SOLF*FL(K)*AVGF(L,K)
       enddo
    enddo

    ! Calculate photolysis rates
    call JRATET(PPJ,T_INPUT,FFF, VALJXX,L_,maxChemLev,NJX)

    ! Free pointers
    QQAA  => NULL()
    SSAA  => NULL()
    PHAA  => NULL()

  END SUBROUTINE PHOTO_JX
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: solar_jx
!
! !DESCRIPTION: Subroutine SOLAR\_JX handles solar zenith angles.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE SOLAR_JX(NDAY,COSSZA,SZA,SOLFX)
!
! !USES:
!
    USE PhysConstants, ONLY : PI
!
! !INPUT PARAMETERS:
!
    REAL(fp), INTENT(IN)  ::  COSSZA
    INTEGER,  INTENT(IN)  ::  NDAY
!
! !OUTPUT VARIABLES:
!
    REAL(fp), INTENT(OUT) ::  SZA,SOLFX
!
! !REMARKS:
!  ---------------------------------------------------------------------
!     NDAY   = integer day of the year (used for solar lat and declin)
!     SZA = solar zenith angle in degrees
!     COSSZA = U0 = cos(SZA)
!     SOLFX = Solar function
!  ---------------------------------------------------------------------
!
! !REVISION HISTORY:
!  28 Mar 2013 - S. D. Eastham - Adapted from Fast-JX v7.0
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    REAL(fp)  PI180

    !=================================================================
    ! SOLAR_JX begins here!
    !=================================================================

    PI180  = PI/180.e+0_fp
    SZA    = acos(MIN(MAX(COSSZA,-1._fp),1._fp))/PI180

    ! Offset used for GEOS-Chem slightly different
    !SOLFX  = 1.e+0_fp-(0.034e+0_fp*cos(dble(NDAY-186)*2.e+0_fp*PI/365.e+0_fp))
    SOLFX  = 1.e+0_fp-(0.034e+0_fp*cos(dble(NDAY-172) &
            *2.e+0_fp*PI/365.e+0_fp))

  END SUBROUTINE SOLAR_JX
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: rd_mie
!
! !DESCRIPTION: Subroutine RD\_MIE retrieves aerosol scattering data for FJX.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE RD_MIE( amIRoot, dryrun, LBRC, NUN, NAMFIL, RC )
!
! !USES:
!
    USE CMN_FJX_Mod, ONLY : TITLAA, NAA, PAA, QAA, RAA, SAA, WAA
    USE ErrCode_Mod
!
! !INPUT PARAMETERS:
!
    LOGICAL,        INTENT(IN)  :: amIRoot     ! On root thread?
    LOGICAL,        INTENT(IN)  :: dryrun      ! Dry run to print inputs?
    LOGICAL,        INTENT(IN)  :: LBRC        ! Brown carbon?
    INTEGER,        INTENT(IN)  :: NUN         ! Logical unit #
    CHARACTER(*),   INTENT(IN)  :: NAMFIL      ! File name
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT) :: RC          ! Success or failure?
!
! !REMARKS:
!   --------------------------------------------------------------------
!     NAMFIL   Name of scattering data file (e.g., FJX_scat.dat)
!     NUN      Channel number for reading data file
!     NAA      Number of categories for scattering phase functions
!     QAA      Aerosol scattering phase functions
!     WAA      5 Wavelengths for the supplied phase functions
!     PAA      Phase function: first 8 terms of expansion
!     RAA      Effective radius associated with aerosol type
!     SAA      Single scattering albedo
!   --------------------------------------------------------------------
!
! !REVISION HISTORY:
!  28 Mar 2013 - S. D. Eastham - Adapted from GEOS-Chem v9-1-3
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    LOGICAL            :: FileExists
    INTEGER            :: I, J, K, NK
#if defined( MODEL_CESM )
    INTEGER            :: ierr
#endif

    ! Strings
    CHARACTER(LEN=78 ) :: TITLE0
    CHARACTER(LEN=255) :: FileMsg, ErrMsg, ThisLoc
#if defined( MODEL_CESM )
    CHARACTER(LEN=*), PARAMETER :: subname = 'rd_mie'
#endif

    !=================================================================
    ! In dry-run mode, print file path to dryrun log and exit.
    ! Otherwise, print file path to stdout and continue.
    !=================================================================

    ! Assume success
    RC = GC_SUCCESS
    ErrMsg = ''
    ThisLoc = ' -> at RD_MIE (in module GeosCore/fast_jx_mod.F90)'

    ! Test if the file exists
    INQUIRE( FILE=TRIM( NamFil ), EXIST=FileExists )

    ! Test if the file exists and define an output string
    IF ( FileExists ) THEN
       FileMsg = 'FAST-JX (RD_MIE): Opening'
    ELSE
       FileMsg = 'FAST-JX (RD_MIE): REQUIRED FILE NOT FOUND'
    ENDIF

    ! Write to stdout for both regular and dry-run simulations
    IF ( amIRoot ) THEN
       WRITE( 6, 300 ) TRIM( FileMsg ), TRIM( NamFil )
300    FORMAT( a, ' ', a )
    ENDIF

    ! For dry-run simulations, return to calling program
    ! For regular simulations, throw an error if we can't find the file.
    IF ( dryrun ) THEN
       RETURN
    ELSE
       IF ( .not. FileExists ) THEN
          WRITE( ErrMsg, 300 ) TRIM( FileMsg ), TRIM( NamFil )
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
    ENDIF

    !=================================================================
    ! RD_MIE begins here -- read data from file
    !=================================================================

#if defined( MODEL_CESM )
    ! Only read file on root thread if using CESM
    IF ( amIRoot ) THEN
#endif

    ! Open file
    open (NUN,FILE=NAMFIL,status='old',form='formatted')

    ! Read header lines
    READ( NUN,'(A)' ) TITLE0
    IF  ( amIRoot ) WRITE( 6, '(1X,A)' ) TITLE0
    READ( NUN,'(A)' ) TITLE0

    !---Read aerosol phase functions:
    read(NUN,'(A10,I5,/)') TITLE0,NAA

    NK=5
    do j=1,NAA
       read(NUN,110) TITLAA(j)
       do k=1,NK
          read(NUN,*) WAA(k,j),QAA(k,j),RAA(k,j),SAA(k,j), &
                      (PAA(i,k,j),i=1,8)
       enddo
    enddo

    ! Brown carbon option
    IF (LBRC) THEN

       ! Overwrite OC entries (36-42 in jv_spec_mie.dat)
       ! with BR entries at end of file (labeled 57-63)
       do j= 36, 42
          read(NUN,110) TITLAA(j)
          do k=1,NK
             read(NUN,*) WAA(k,j),QAA(k,j),RAA(k,j),SAA(k,j), &
                         (PAA(i,k,j),i=1,8)
          enddo
       enddo

    ENDIF

    close(NUN)

#if defined( MODEL_CESM )
    ENDIF

    CALL MPI_BCAST( QAA,    Size(QAA),  mpi_real8,     masterprocid, mpicom, ierr )
    IF ( ierr /= mpi_success ) CALL endrun(subname//': MPI_BCAST ERROR: QAA')
    CALL MPI_BCAST( WAA,    Size(WAA),  mpi_real8,     masterprocid, mpicom, ierr )
    IF ( ierr /= mpi_success ) CALL endrun(subname//': MPI_BCAST ERROR: WAA')
    CALL MPI_BCAST( PAA,    Size(PAA),  mpi_real8,     masterprocid, mpicom, ierr )
    IF ( ierr /= mpi_success ) CALL endrun(subname//': MPI_BCAST ERROR: PAA')
    CALL MPI_BCAST( RAA,    Size(RAA),  mpi_real8,     masterprocid, mpicom, ierr )
    IF ( ierr /= mpi_success ) CALL endrun(subname//': MPI_BCAST ERROR: RAA')
    CALL MPI_BCAST( SAA,    Size(SAA),  mpi_real8,     masterprocid, mpicom, ierr )
    IF ( ierr /= mpi_success ) CALL endrun(subname//': MPI_BCAST ERROR: SAA')
    CALL MPI_BCAST( NAA,    1,          mpi_integer,   masterprocid, mpicom, ierr )
    IF ( ierr /= mpi_success ) CALL endrun(subname//': MPI_BCAST ERROR: NAA')
    CALL MPI_BCAST( TITLAA, 80*A_,      mpi_character, masterprocid, mpicom, ierr )
    IF ( ierr /= mpi_success ) CALL endrun(subname//': MPI_BCAST ERROR: TITLAA')
#endif

    IF ( amIRoot ) THEN
       write(6,'(a,9f8.1)') ' Aerosol optical: r-eff/rho/Q(@wavel):', &
                            (WAA(K,1),K=1,5)
       do J=1,NAA
          write(6,'(1x,A)') TRIM(TITLAA(J))
          write(6,'(3x,I2,A,9F8.1)') J,'  wavel=',(WAA(K,J),K=1,NK)
          write(6,'(3x,I2,A,9F8.4)') J,'  Qext =',(QAA(K,J),K=1,NK)
       enddo
    ENDIF

110 format(3x,a80)

  END SUBROUTINE RD_MIE
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: rd_xxx
!
! !DESCRIPTION: Subroutine RD\_XXX reads in wavelength bins, solar fluxes,
!  Rayleigh and temperature-dependent cross-sections.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE RD_XXX ( amIRoot, dryrun, NUN, NAMFIL, RC )
!
! !USES:
!
    USE CMN_FJX_Mod, ONLY : W_, WX_, X_, NJX, NW1, NW2
    USE CMN_FJX_Mod, ONLY : TITLEJX, WL, FL, QRAYL, QO2, QO3, Q1D
    USE CMN_FJX_Mod, ONLY : LQQ, QQQ, SQQ, TQQ
    USE ErrCode_Mod
!
! !INPUT PARAMETERS:
!
    LOGICAL,        INTENT(IN)  :: amIRoot
    LOGICAL,        INTENT(IN)  :: dryrun
    INTEGER,        INTENT(IN)  :: NUN
    CHARACTER(*),   INTENT(IN)  :: NAMFIL
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT) :: RC
!
! !REMARKS:
!    NEW v-6.8  now allow 1 to 3 sets of X-sects for T or P
!           LQQ = 1, 2, or 3 to determine interpolation with T or P
!           IF the temperatures TQQQ are <0, then use as pressure interp (hPa)
!           NB - the temperatures and pressures must be increasing
!    NEW v-6.4  changed to collapse wavelengths & x-sections to Trop-only:
!           WX_ = 18 should match the JX_spec.dat wavelengths
!           W_ = 12 (Trop-only) or 18 (std) is set in (CMN_FJX.F).
!       if W_=12 then drop strat wavels, and drop x-sects (e.g. N2O, ...)
!           W_ = 8, reverts to quick fix:  fast-J (12-18) plus bin (5) scaled
!                                                                             .
!   --------------------------------------------------------------------
!     NAMFIL   Name of spectral data file (FJX_spec.dat) >> j2 for fast-J2
!     NUN      Channel number for reading data file
!
!     NJX    Number of species to calculate J-values for
!     NWWW     Number of wavelength bins, from 1:NWWW
!     WBIN     Boundaries of wavelength bins
!     WL       Centres of wavelength bins - 'effective wavelength'
!     FL       Solar flux incident on top of atmosphere (cm-2.s-1)
!     QRAYL    Rayleigh parameters (effective cross-section) (cm2)
!     QO2      O2 cross-sections
!     QO3      O3 cross-sections
!     Q1D      O3 => O(1D) quantum yield
!     TQQ      Temperature for supplied cross sections
!     QQQ      Supplied cross sections in each wavelength bin (cm2)
!   --------------------------------------------------------------------
!
! !REVISION HISTORY:
!  28 Mar 2013 - S. D. Eastham - Copied from Fast-JX v7.0
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    LOGICAL            :: FileExists
    INTEGER            :: I, J, JJ, K, IW, NQRD, NWWW, LQ
#if defined( MODEL_CESM )
    INTEGER            :: ierr
#endif
    REAL(fp)           :: TQQ2

    ! Arrays
    REAL(fp)           :: QQ2(199)

    ! Strings
    CHARACTER(LEN=255) :: FileMsg, FileStatus
    CHARACTER(LEN=255) :: ErrMsg,  ThisLoc
    CHARACTER(LEN=78)  :: TITLE0
    CHARACTER(LEN=6 )  :: TITLEJ2, TITLEJ3
    CHARACTER(LEN=1 )  :: TSTRAT
#if defined( MODEL_CESM )
    CHARACTER(LEN=*), PARAMETER :: subname = 'rd_xxx'
#endif

    !=================================================================
    ! In dry-run mode, print file path to dryrun log and exit.
    ! Otherwise, print file path to stdout and continue.
    !=================================================================

    ! Initialize
    RC      = GC_SUCCESS
    ErrMsg  = ''
    ThisLoc = '-> at RD_XXX (in module GeosCore/fast_jx_mod.F)'

    ! Test if the file exists
    INQUIRE( FILE=TRIM( NamFil ), EXIST=FileExists )

    ! Test if the file exists and define an output string
    IF ( FileExists ) THEN
       FileMsg    = 'FAST-JX (RD_XXX): Opening'
    ELSE
       FileMsg    = 'FAST-JX (RD_XXX): REQUIRED FILE NOT FOUND'
    ENDIF

    ! Write to stdout for both regular and dry-run simulations
    IF ( amIRoot ) THEN
       WRITE( 6, 300 ) TRIM( FileMsg ), TRIM( NamFil )
300    FORMAT( a, ' ', a )
    ENDIF

    ! For dry-run simulations, return to calling program.
    ! For regular simulations, throw an error if we can't find the file.
    IF ( dryrun ) THEN
       RETURN
    ELSE
       IF ( .not. FileExists ) THEN
          WRITE( ErrMsg, 300 ) TRIM( FileMsg ), TRIM( NamFil )
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
    ENDIF

    !=================================================================
    ! RD_XXX begins here!
    !=================================================================

    ! Initialization
    TQQ(:,:) = 0.e+0_fp

    ! -------spectral data----set for new format data------------------
    !         note that X_ = max # Xsects read in
    !                   NJX = # fast-JX J-values derived from this (.le. X_)

    ! >>>> W_ = 12 <<<< means trop-only, discard WL #1-4 and #9-10, some X-sects

#if defined( MODEL_CESM )
    ! Only read file on root thread if using CESM
    IF ( amIRoot ) THEN
#endif

    ! Open file
    open (NUN,FILE=NAMFIL,status='old',form='formatted')
    read (NUN,100) TITLE0

    ! -note that NQRD is not used any more, a read until 'endofJ' is performed
    read (NUN,101) NQRD,NWWW
    NW1 = 1
    NW2 = NWWW
    IF ( amIRoot ) THEN
       write(6,'(1x,a)') TITLE0
       write(6,'(i8)') NWWW
    ENDIF
    ! -J-values:  1=O2, 2=O3P,3=O3D 4=readin Xsects
    read (NUN,102) (WL(IW),IW=1,NWWW)
    read (NUN,102) (FL(IW),IW=1,NWWW)
    read (NUN,102) (QRAYL(IW),IW=1,NWWW)

    ! Read O2 X-sects, O3 X-sects, O3=>O(1D) quant yields (each at 3 temps)
    ! NB the O3 and q-O3-O1D are at different temperatures and cannot be combined
    read (NUN,103) TITLEJX(1),TQQ(1,1), (QO2(IW,1),IW=1,NWWW)
    read (NUN,103) TITLEJ2,   TQQ(2,1), (QO2(IW,2),IW=1,NWWW)
    read (NUN,103) TITLEJ3,   TQQ(3,1), (QO2(IW,3),IW=1,NWWW)

    read (NUN,103) TITLEJX(2),TQQ(1,2), (QO3(IW,1),IW=1,NWWW)
    read (NUN,103) TITLEJ2,   TQQ(2,2), (QO3(IW,2),IW=1,NWWW)
    read (NUN,103) TITLEJ3,   TQQ(3,2), (QO3(IW,3),IW=1,NWWW)

    read (NUN,103) TITLEJX(3),TQQ(1,3), (Q1D(IW,1),IW=1,NWWW)
    read (NUN,103) TITLEJ2,   TQQ(2,3), (Q1D(IW,2),IW=1,NWWW)
    read (NUN,103) TITLEJ3,   TQQ(3,3), (Q1D(IW,3),IW=1,NWWW)

    SQQ(1) = ' '
    SQQ(2) = ' '
    SQQ(3) = ' '

    LQQ(1) = 3
    LQQ(2) = 3
    LQQ(3) = 3

    ! Read remaining species:  X-sections at 1-2-3 T_s
    JJ = 3
    do I=1,9999

       ! try to read in 3 X-sects per J-value (JJ)
       read (NUN,104) TITLEJ2,TSTRAT,TQQ2,(QQ2(IW),IW=1,NWWW)
       if (TITLEJ2 .eq. 'endofJ') goto 1

       ! skip stratosphere only J's (denoted by 'x')if W_<18 => trop-only J's
       if (W_.eq.18 .or. TSTRAT.ne.'x') then
          if (TITLEJ2 .ne. TITLEJX(JJ)) then
             JJ = JJ+1

             if (JJ .gt. X_) then
                call EXITC(' RD_XXX:  X_ not large enough for Xsects read in')
             endif

             TITLEJX(JJ) = TITLEJ2
             LQQ(JJ) = 1
             SQQ(JJ) = TSTRAT
             LQ = LQQ(JJ)
             TQQ(LQ,JJ) = TQQ2
             do IW = 1,NWWW
                QQQ(IW,LQ,JJ) = QQ2(IW)
             enddo
          else
             LQQ(JJ) = LQQ(JJ)+1
             if (LQQ(JJ) .le. 3) then
                LQ = LQQ(JJ)
                TQQ(LQ,JJ) = TQQ2
                do IW = 1,NWWW
                   QQQ(IW,LQ,JJ) = QQ2(IW)
                enddo
             endif
          endif
       endif
    enddo
1   continue
    NJX = JJ

    do J = 1,NJX
       if ( amIRoot ) then
          write(6,200) J,TITLEJX(J),SQQ(J),LQQ(J),(TQQ(I,J),I=1,LQQ(J))
       endif
       ! need to check that TQQ is monotonically increasing:
       if (LQQ(J) .eq. 3) then
          if (TQQ(2,J) .ge. TQQ(3,J)) then
             call EXITC ('TQQ out of order')
          endif
          if (TQQ(1,J) .ge. TQQ(2,J)) then
             call EXITC ('TQQ out of order')
          endif
       endif
       if (LQQ(J) .eq. 2) then
          if (TQQ(1,J) .ge. TQQ(2,J)) then
             call EXITC ('TQQ out of order')
          endif
       endif
    enddo

    ! check on doingpressure interp
    ! check on consolidating Qo2 and others into
    ! wrte a newFJX_J2J.dat for mapping on fjx Xsects

    ! truncate number of wavelengths to do troposphere-only
    if (W_ .ne. WX_) then
       ! TROP-ONLY
       if (W_ .eq. 12) then
          if ( amIRoot ) then
             write(6,'(a)') &
                  ' >>>TROP-ONLY reduce wavelengths to 12, drop strat X-sects'
          endif
          NW2 = 12
          do IW = 1,4
             WL(IW) = WL(IW+4)
             FL(IW) = FL(IW+4)
             QRAYL(IW) = QRAYL(IW+4)
             do K = 1,3
                QO2(IW,K) = QO2(IW+4,K)
                QO3(IW,K) = QO3(IW+4,K)
                Q1D(IW,K) = Q1D(IW+4,K)
             enddo
             do J = 4,NJX
                do LQ=1,LQQ(J)
                   QQQ(IW,LQ,J) = QQQ(IW+4,LQ,J)
                enddo
             enddo
          enddo
          do IW = 5,12
             WL(IW) = WL(IW+6)
             FL(IW) = FL(IW+6)
             QRAYL(IW) = QRAYL(IW+6)
             do K = 1,3
                QO2(IW,K) = QO2(IW+6,K)
                QO3(IW,K) = QO3(IW+6,K)
                Q1D(IW,K) = Q1D(IW+6,K)
             enddo
             do J = 4,NJX
                do LQ=1,LQQ(J)
                   QQQ(IW,LQ,J) = QQQ(IW+6,LQ,J)
                enddo
             enddo
          enddo
          ! TROP-QUICK  (must scale solar flux for W=5)
       elseif (W_ .eq. 8) then
          if ( amIRoot ) then
             write(6,'(a)') &
                  ' >>>TROP-QUICK reduce wavelengths to 8, drop strat X-sects'
          endif
          NW2 = 8
          do IW = 1,1
             WL(IW) = WL(IW+4)
             FL(IW) = FL(IW+4)  * 2.e+0_fp
             QRAYL(IW) = QRAYL(IW+4)
             do K = 1,3
                QO2(IW,K) = QO2(IW+4,K)
                QO3(IW,K) = QO3(IW+4,K)
                Q1D(IW,K) = Q1D(IW+4,K)
             enddo
             do J = 4,NJX
                do LQ=1,LQQ(J)
                   QQQ(IW,LQ,J) = QQQ(IW+4,LQ,J)
                enddo
             enddo
          enddo
          do IW = 2,8
             WL(IW) = WL(IW+10)
             FL(IW) = FL(IW+10)
             QRAYL(IW) = QRAYL(IW+10)
             do K = 1,3
                QO2(IW,K) = QO2(IW+10,K)
                QO3(IW,K) = QO3(IW+10,K)
                Q1D(IW,K) = Q1D(IW+10,K)
             enddo
             do J = 4,NJX
                do LQ=1,LQQ(J)
                   QQQ(IW,LQ,J) = QQQ(IW+10,LQ,J)
                enddo
             enddo
          enddo

       else
          call EXITC(' no. wavelengths wrong: W_ .ne. 8,12,18')
       endif
    endif

    close(NUN)

#if defined( MODEL_CESM )
    ENDIF

    CALL MPI_BCAST( NJX,       1,            mpi_real8,     masterprocid, mpirun, ierr )
    IF ( ierr /= mpi_success ) CALL endrun(subname//': MPI_BCAST ERROR: NJX')
    CALL MPI_BCAST( NW1,       1,            mpi_real8,     masterprocid, mpirun, ierr )
    IF ( ierr /= mpi_success ) CALL endrun(subname//': MPI_BCAST ERROR: NW1')
    CALL MPI_BCAST( NW2,       1,            mpi_real8,     masterprocid, mpirun, ierr )
    IF ( ierr /= mpi_success ) CALL endrun(subname//': MPI_BCAST ERROR: NW2')
    CALL MPI_BCAST( WBIN,      Size(WBIN),   mpi_real8,     masterprocid, mpirun, ierr )
    IF ( ierr /= mpi_success ) CALL endrun(subname//': MPI_BCAST ERROR: WBIN')
    CALL MPI_BCAST( WL,        Size(WL),     mpi_real8,     masterprocid, mpirun, ierr )
    IF ( ierr /= mpi_success ) CALL endrun(subname//': MPI_BCAST ERROR: WL')
    CALL MPI_BCAST( FL,        Size(FL),     mpi_real8,     masterprocid, mpirun, ierr )
    IF ( ierr /= mpi_success ) CALL endrun(subname//': MPI_BCAST ERROR: FL')
    CALL MPI_BCAST( QO2,       Size(QO2),    mpi_real8,     masterprocid, mpirun, ierr )
    IF ( ierr /= mpi_success ) CALL endrun(subname//': MPI_BCAST ERROR: QO2')
    CALL MPI_BCAST( QO3,       Size(QO3),    mpi_real8,     masterprocid, mpirun, ierr )
    IF ( ierr /= mpi_success ) CALL endrun(subname//': MPI_BCAST ERROR: QO3')
    CALL MPI_BCAST( Q1D,       Size(Q1D),    mpi_real8,     masterprocid, mpirun, ierr )
    IF ( ierr /= mpi_success ) CALL endrun(subname//': MPI_BCAST ERROR: Q1D')
    CALL MPI_BCAST( QQQ,       Size(QQQ),    mpi_real8,     masterprocid, mpirun, ierr )
    IF ( ierr /= mpi_success ) CALL endrun(subname//': MPI_BCAST ERROR: QQQ')
    CALL MPI_BCAST( QRAYL,     Size(QRAYL),  mpi_real8,     masterprocid, mpirun, ierr )
    IF ( ierr /= mpi_success ) CALL endrun(subname//': MPI_BCAST ERROR: QRAYL')
    CALL MPI_BCAST( TQQ,       Size(TQQ),    mpi_real8,     masterprocid, mpirun, ierr )
    IF ( ierr /= mpi_success ) CALL endrun(subname//': MPI_BCAST ERROR: TQQ')
    CALL MPI_BCAST( LQQ,       Size(LQQ),    mpi_integer,   masterprocid, mpirun, ierr )
    IF ( ierr /= mpi_success ) CALL endrun(subname//': MPI_BCAST ERROR: LQQ')
    CALL MPI_BCAST( TITLEJX,   X_*6,         mpi_character, masterprocid, mpirun, ierr )
    IF ( ierr /= mpi_success ) CALL endrun(subname//': MPI_BCAST ERROR: TITLEJX')
    CALL MPI_BCAST( SQQ,       X_*1,         mpi_character, masterprocid, mpirun, ierr )
    IF ( ierr /= mpi_success ) CALL endrun(subname//': MPI_BCAST ERROR: SQQ')
#endif

100 format(a)
101 format(10x,5i5)
102 format(10x,    6e10.3/(10x,6e10.3)/(10x,6e10.3))
103 format(a6,1x,f3.0,6e10.3/(10x,6e10.3)/(10x,6e10.3))
104 format(a6,a1,f3.0,6e10.3/(10x,6e10.3)/(10x,6e10.3))
200 format(1x,' x-sect:',i3,a10,a4,i5,3(3x,f6.2))
201 format(' Number of x-sections supplied to Fast-J2: ',i3,/, &
           ' Maximum number allowed (X_) only set to: ',i3,    &
           ' - increase in cmn_FJX.f')

  END SUBROUTINE RD_XXX
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: rd_js_jx
!
! !DESCRIPTION: Subroutine RD\_JS\_JX reads in 'FJX\_j2j.dat', which defines
!  the mapping of Fast-JX J's (TITLEJX(1:NJX)) onto the CTM reactions.
!  Reaction number JJ, named T\_REACT, uses Fast-JX's T\_FJX (including scaling
!  factor F\_FJX).
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE RD_JS_JX( amIRoot, dryrun, NUNIT, NAMFIL, TITLEJX, NJXX, RC )
!
! !USES:
!
    USE CMN_FJX_Mod, ONLY : M2_, JVN_, JIND
    USE CMN_FJX_Mod, ONLY : JLABEL, JFACTA, NRATJ, BRANCH, RNAMES
    USE Charpak_Mod,    ONLY : CStrip
    USE ErrCode_Mod
!
! !INPUT PARAMETERS:
!
    LOGICAL,          INTENT(IN)                  :: amIRoot
    LOGICAL,          INTENT(IN)                  :: dryrun
    INTEGER,          INTENT(IN)                  :: NUNIT
    INTEGER,          INTENT(IN)                  :: NJXX
    CHARACTER(LEN=*), INTENT(IN)                  :: NAMFIL
    CHARACTER(LEN=6), INTENT(IN), DIMENSION(NJXX) :: TITLEJX
!
! !OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(OUT)                 :: RC
!
! !REMARKS:
!  Now flag special reactions that are to be adjusted for FlexChem later.
!
! !REVISION HISTORY:
!  28 Mar 2013 - S. D. Eastham - Copied from Fast-JX v7.0
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    LOGICAL            :: FileExists
    INTEGER            :: J, JJ, K
#if defined( MODEL_CESM )
    INTEGER            :: ierr
#endif
    REAL(fp)           :: F_FJX

    ! Strings
    CHARACTER(LEN=6  ) :: T_FJX
    CHARACTER(LEN=50 ) :: T_REACT
    CHARACTER(LEN=50 ) :: TEXT
    CHARACTER(LEN=120) :: CLINE
    CHARACTER(LEN=255) :: ErrMsg, ThisLoc, FileMsg

    ! String arrays
    CHARACTER(LEN=6)   :: JMAP(JVN_)
#if defined( MODEL_CESM )
    CHARACTER(LEN=*), PARAMETER :: subname = 'rd_js_jx'
#endif

    !========================================================================
    ! RD_JS_JX begins here!
    !========================================================================

    ! Initialize
    RC       = GC_SUCCESS
    ErrMsg   = ''
    ThisLoc  = ' -> at Rd_Js_Jx (in module GeosCore/fast_jx_mod.F90)'

    !========================================================================
    ! In dry-run mode, print file path to dryrun log and exit.
    ! Otherwise, print file path to stdout and continue.
    !========================================================================

    ! Test if the file exists
    INQUIRE( FILE=TRIM( NamFil ), EXIST=FileExists )

    ! Test if the file exists and define an output string
    IF ( FileExists ) THEN
       FileMsg = 'FAST-JX (RD_JS_JX): Opening'
    ELSE
       FileMsg = 'FAST-JX (RD_JS_JX): REQUIRED FILE NOT FOUND'
    ENDIF

    ! Write to stdout for both regular and dry-run simulations
    IF ( amIRoot ) THEN
       WRITE( 6, 300 ) TRIM( FileMsg ), TRIM( NamFil )
300    FORMAT( a, ' ', a )
    ENDIF

    ! For dry-run simulations, return to calling program.
    ! For regular simulations, throw an error if we can't find the file.
    IF ( dryrun ) THEN
       RETURN
    ELSE
       IF ( .not. FileExists ) THEN
          WRITE( ErrMsg, 300 ) TRIM( FileMsg ), TRIM( NamFil )
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
    ENDIF

    !========================================================================
    ! Read the FJX_j2j.dat file to map model specific J's onto fast-JX J's
    ! The chemistry code title describes fully the reaction (a50)
    ! Blank (unfilled) chemistry J's are unmapped
    ! The number NRATJ is the last JJ readin that is .le. JVN
    ! include fractional quantum yield for the fast-JX J's
    !========================================================================

    JLABEL(:) = '------'
    JMAP(:)   = '------'
    JFACTA(:) = 0.e+0_fp

#if defined( MODEL_CESM )
    ! Only read file on root thread if using CESM
    IF ( amIRoot ) THEN
#endif

    ! Open file
    open (NUNIT,file=NAMFIL,status='old',form='formatted')

    read (NUNIT,'(a)') CLINE
    IF ( amIRoot ) THEN
       write(6,'(a)') CLINE
    ENDIF
    do J = 1,JVN_
       read (NUNIT,'(i4,1x,a50,4x,f5.3,2x,a6)') JJ,T_REACT,F_FJX,T_FJX
       IF (JJ.gt.JVN_) THEN
          IF ( JJ .eq. 9999 ) THEN
             close(NUNIT)
             exit
          ELSE
             ErrMsg = 'Number of reactions in FJX_j2j.dat exceeds JVN_.' //&
                      'Adjust JVN_ in CMN_FJX_mod.F90 to get past error.'
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN
          ENDIF
       ENDIF

       JLABEL(JJ) = T_REACT
       JFACTA(JJ) = F_FJX
       JMAP(JJ) = T_FJX
       NRATJ = JJ
       ! SDE 03/31/13: Check number of branches
       ! Note that the order of the branches in
       ! globchem.dat must match the order in
       ! FJX_j2j.dat
       READ (T_REACT(1:10),"(a10)") RNAMES(JJ)
       RNAMES(JJ) = TRIM(RNAMES(JJ))
       BRANCH(JJ) = 1
       DO K=1,(JJ-1)
          IF (RNAMES(JJ) == RNAMES(K)) THEN
             BRANCH(JJ) = BRANCH(K) + 1
          ENDIF
       ENDDO
    enddo

20  close(NUNIT)

#if defined( MODEL_CESM )
    ENDIF

    CALL MPI_BCAST( JLABEL, JVN_*50, mpi_character, masterprocid, mpirun )
    IF ( ierr /= mpi_success ) CALL endrun(subname//': MPI_BCAST ERROR: JLABEL')
    CALL MPI_BCAST( JFACTA, JVN_,    mpi_real8,     masterprocid, mpirun )
    IF ( ierr /= mpi_success ) CALL endrun(subname//': MPI_BCAST ERROR: JFACTA')
    CALL MPI_BCAST( JMAP,   JVN_*6,  mpi_character, masterprocid, mpirun )
    IF ( ierr /= mpi_success ) CALL endrun(subname//': MPI_BCAST ERROR: JMAP')
    CALL MPI_BCAST( NRATJ,  1,       mpi_integer,   masterprocid, mpirun )
    IF ( ierr /= mpi_success ) CALL endrun(subname//': MPI_BCAST ERROR: NRATJ')
    CALL MPI_BCAST( RNAMES, JVN_*10, mpi_character, masterprocid, mpirun )
    IF ( ierr /= mpi_success ) CALL endrun(subname//': MPI_BCAST ERROR: RNAMES')
    CALL MPI_BCAST( BRANCH, JVN_,    mpi_integer,   masterprocid, mpirun )
    IF ( ierr /= mpi_success ) CALL endrun(subname//': MPI_BCAST ERROR: BRANCH')
#endif

    ! Zero / Set index arrays that map Jvalue(j) onto rates
    do K = 1,NRATJ
       JIND(K) = 0
       do J = 1,NJXX
          T_FJX = TITLEJX(J)
          if (JMAP(K) .eq. TITLEJX(J)) then
             JIND(K)=J
          endif
       enddo
    enddo

    IF ( amIRoot ) THEN
       write(6,'(a,i4,a)')'Photochemistry Scheme with',NRATJ,' J-values'
    ENDIF
    do K=1,NRATJ
       if (JMAP(K) .ne. '------' ) then
          J = JIND(K)
          IF ( amIRoot ) THEN
             if (J.eq.0) then
                write(6,'(i5,1x,a50,f6.3,a,1x,a6)') K,JLABEL(K),JFACTA(K), &
                     ' no mapping onto fast-JX',JMAP(K)
             else
                write(6,'(i5,1x,a50,f6.3,a,i4,1x,a6)') K,JLABEL(K),JFACTA(K), &
                     ' mapped to FJX:',J,TITLEJX(J)
             endif
          ENDIF
       endif
    enddo

  END SUBROUTINE RD_JS_JX
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: blkslv
!
! !DESCRIPTION: Subroutine BLKSLV solves the block tri-diagonal system
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE BLKSLV(FJ,POMEGA,FZ,ZTAU,ZFLUX,RFL,PM,PM0,FJTOP,FJBOT,ND)
!
! !USES:
!
  USE CMN_FJX_Mod, ONLY : EMU, M_, M2_, N_, W_, WT
!
! !INPUT PARAMETERS:
!
    INTEGER,  INTENT(IN)  :: ND
    REAL(fp), INTENT(IN)  :: POMEGA(M2_,N_,W_)
    REAL(fp), INTENT(IN)  :: FZ(N_,W_)
    REAL(fp), INTENT(IN)  :: ZTAU(N_,W_)
    REAL(fp), INTENT(IN)  :: PM(M_,M2_)
    REAL(fp), INTENT(IN)  :: PM0(M2_)
    REAL(fp), INTENT(IN)  :: RFL(W_)
    REAL(fp), INTENT(IN)  :: ZFLUX(W_)
!
! !OUTPUT PARAMETERS:
!
    REAL(fp), INTENT(OUT) :: FJ(N_,W_)
    REAL(fp), INTENT(OUT) :: FJTOP(W_)
    REAL(fp), INTENT(OUT) :: FJBOT(W_)
!
! !REMARKS:
! The block tri-diagonal system:
!       A(I)*X(I-1) + B(I)*X(I) + C(I)*X(I+1) = H(I)
!
! !REVISION HISTORY:
!  27 Mar 2013 - S. D. Eastham - Copied from GEOS-Chem v9-01-03
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    REAL(fp), DIMENSION(M_,N_,W_)    ::  A,C,H,   RR
    REAL(fp), DIMENSION(M_,M_,N_,W_) ::  B,AA,CC,  DD
    REAL(fp), DIMENSION(M_,M_)       ::  E
    REAL(fp)  SUMB,SUMBX,SUMT
    INTEGER I, J, K, L

    !=================================================================
    ! BLKSLV begins here!
    !=================================================================

    do K = 1,W_
       call GEN_ID (POMEGA(1,1,K),FZ(1,K),ZTAU(1,K),ZFLUX(K),RFL(K), &
                    PM,PM0, B(1,1,1,K),CC(1,1,1,K),AA(1,1,1,K),      &
                    A(1,1,K),H(1,1,K),C(1,1,K), ND)
    enddo

    do K = 1,W_
       ! UPPER BOUNDARY L=1
       L = 1
       do J = 1,M_
       do I = 1,M_
          E(I,J) = B(I,J,1,K)
       enddo
       enddo

       ! setup L & U matrices
       E(2,1) = E(2,1)/E(1,1)
       E(2,2) = E(2,2)-E(2,1)*E(1,2)
       E(2,3) = E(2,3)-E(2,1)*E(1,3)
       E(2,4) = E(2,4)-E(2,1)*E(1,4)
       E(3,1) = E(3,1)/E(1,1)
       E(3,2) = (E(3,2)-E(3,1)*E(1,2))/E(2,2)
       E(3,3) = E(3,3)-E(3,1)*E(1,3)-E(3,2)*E(2,3)
       E(3,4) = E(3,4)-E(3,1)*E(1,4)-E(3,2)*E(2,4)
       E(4,1) = E(4,1)/E(1,1)
       E(4,2) = (E(4,2)-E(4,1)*E(1,2))/E(2,2)
       E(4,3) = (E(4,3)-E(4,1)*E(1,3)-E(4,2)*E(2,3))/E(3,3)
       E(4,4) = E(4,4)-E(4,1)*E(1,4)-E(4,2)*E(2,4)-E(4,3)*E(3,4)
       ! invert L
       E(4,3) = -E(4,3)
       E(4,2) = -E(4,2)-E(4,3)*E(3,2)
       E(4,1) = -E(4,1)-E(4,2)*E(2,1)-E(4,3)*E(3,1)
       E(3,2) = -E(3,2)
       E(3,1) = -E(3,1)-E(3,2)*E(2,1)
       E(2,1) = -E(2,1)
       ! invert U
       E(4,4) = 1.e+0_fp/E(4,4)
       E(3,4) = -E(3,4)*E(4,4)/E(3,3)
       E(3,3) = 1.e+0_fp/E(3,3)
       E(2,4) = -(E(2,3)*E(3,4)+E(2,4)*E(4,4))/E(2,2)
       E(2,3) = -E(2,3)*E(3,3)/E(2,2)
       E(2,2) = 1.e+0_fp/E(2,2)
       E(1,4) = -(E(1,2)*E(2,4)+E(1,3)*E(3,4)+E(1,4)*E(4,4))/E(1,1)
       E(1,3) = -(E(1,2)*E(2,3)+E(1,3)*E(3,3))/E(1,1)
       E(1,2) = -E(1,2)*E(2,2)/E(1,1)
       E(1,1) = 1.e+0_fp/E(1,1)
       ! multiply U-invers * L-inverse
       E(1,1) = E(1,1)+E(1,2)*E(2,1)+E(1,3)*E(3,1)+E(1,4)*E(4,1)
       E(1,2) = E(1,2)+E(1,3)*E(3,2)+E(1,4)*E(4,2)
       E(1,3) = E(1,3)+E(1,4)*E(4,3)
       E(2,1) = E(2,2)*E(2,1)+E(2,3)*E(3,1)+E(2,4)*E(4,1)
       E(2,2) = E(2,2)+E(2,3)*E(3,2)+E(2,4)*E(4,2)
       E(2,3) = E(2,3)+E(2,4)*E(4,3)
       E(3,1) = E(3,3)*E(3,1)+E(3,4)*E(4,1)
       E(3,2) = E(3,3)*E(3,2)+E(3,4)*E(4,2)
       E(3,3) = E(3,3)+E(3,4)*E(4,3)
       E(4,1) = E(4,4)*E(4,1)
       E(4,2) = E(4,4)*E(4,2)
       E(4,3) = E(4,4)*E(4,3)

       do J = 1,M_
          do I = 1,M_
             DD(I,J,1,K) = -E(I,1)*CC(1,J,1,K)-E(I,2)*CC(2,J,1,K) &
                           -E(I,3)*CC(3,J,1,K)-E(I,4)*CC(4,J,1,K)
          enddo
          RR(J,1,K) = E(J,1)*H(1,1,K)+E(J,2)*H(2,1,K) &
                    + E(J,3)*H(3,1,K)+E(J,4)*H(4,1,K)
       enddo

       ! CONTINUE THROUGH ALL DEPTH POINTS ID=2 TO ID=ND-1
       do L = 2,ND-1

          do J = 1,M_
             do I = 1,M_
                B(I,J,L,K) = B(I,J,L,K) + A(I,L,K)*DD(I,J,L-1,K)
             enddo
             H(J,L,K) = H(J,L,K) - A(J,L,K)*RR(J,L-1,K)
          enddo

          do J = 1,M_
             do I = 1,M_
                E(I,J) = B(I,J,L,K)
             enddo
          enddo

          ! setup L & U matrices
          E(2,1) = E(2,1)/E(1,1)
          E(2,2) = E(2,2)-E(2,1)*E(1,2)
          E(2,3) = E(2,3)-E(2,1)*E(1,3)
          E(2,4) = E(2,4)-E(2,1)*E(1,4)
          E(3,1) = E(3,1)/E(1,1)
          E(3,2) = (E(3,2)-E(3,1)*E(1,2))/E(2,2)
          E(3,3) = E(3,3)-E(3,1)*E(1,3)-E(3,2)*E(2,3)
          E(3,4) = E(3,4)-E(3,1)*E(1,4)-E(3,2)*E(2,4)
          E(4,1) = E(4,1)/E(1,1)
          E(4,2) = (E(4,2)-E(4,1)*E(1,2))/E(2,2)
          E(4,3) = (E(4,3)-E(4,1)*E(1,3)-E(4,2)*E(2,3))/E(3,3)
          E(4,4) = E(4,4)-E(4,1)*E(1,4)-E(4,2)*E(2,4)-E(4,3)*E(3,4)
          ! invert L
          E(4,3) = -E(4,3)
          E(4,2) = -E(4,2)-E(4,3)*E(3,2)
          E(4,1) = -E(4,1)-E(4,2)*E(2,1)-E(4,3)*E(3,1)
          E(3,2) = -E(3,2)
          E(3,1) = -E(3,1)-E(3,2)*E(2,1)
          E(2,1) = -E(2,1)
          ! invert U
          E(4,4) = 1.e+0_fp/E(4,4)
          E(3,4) = -E(3,4)*E(4,4)/E(3,3)
          E(3,3) = 1.e+0_fp/E(3,3)
          E(2,4) = -(E(2,3)*E(3,4)+E(2,4)*E(4,4))/E(2,2)
          E(2,3) = -E(2,3)*E(3,3)/E(2,2)
          E(2,2) = 1.e+0_fp/E(2,2)
          E(1,4) = -(E(1,2)*E(2,4)+E(1,3)*E(3,4)+E(1,4)*E(4,4))/E(1,1)
          E(1,3) = -(E(1,2)*E(2,3)+E(1,3)*E(3,3))/E(1,1)
          E(1,2) = -E(1,2)*E(2,2)/E(1,1)
          E(1,1) = 1.e+0_fp/E(1,1)
          ! multiply U-invers * L-inverse
          E(1,1) = E(1,1)+E(1,2)*E(2,1)+E(1,3)*E(3,1)+E(1,4)*E(4,1)
          E(1,2) = E(1,2)+E(1,3)*E(3,2)+E(1,4)*E(4,2)
          E(1,3) = E(1,3)+E(1,4)*E(4,3)
          E(2,1) = E(2,2)*E(2,1)+E(2,3)*E(3,1)+E(2,4)*E(4,1)
          E(2,2) = E(2,2)+E(2,3)*E(3,2)+E(2,4)*E(4,2)
          E(2,3) = E(2,3)+E(2,4)*E(4,3)
          E(3,1) = E(3,3)*E(3,1)+E(3,4)*E(4,1)
          E(3,2) = E(3,3)*E(3,2)+E(3,4)*E(4,2)
          E(3,3) = E(3,3)+E(3,4)*E(4,3)
          E(4,1) = E(4,4)*E(4,1)
          E(4,2) = E(4,4)*E(4,2)
          E(4,3) = E(4,4)*E(4,3)

          do J = 1,M_
             do I = 1,M_
                DD(I,J,L,K) = - E(I,J)*C(J,L,K)
             enddo
             RR(J,L,K) = E(J,1)*H(1,L,K)+E(J,2)*H(2,L,K) &
                       + E(J,3)*H(3,L,K)+E(J,4)*H(4,L,K)
          enddo

       enddo

       ! FINAL DEPTH POINT: L=ND
       L = ND
       do J = 1,M_
          do I = 1,M_
             B(I,J,L,K) = B(I,J,L,K) &
                   + AA(I,1,L,K)*DD(1,J,L-1,K) + AA(I,2,L,K)*DD(2,J,L-1,K) &
                   + AA(I,3,L,K)*DD(3,J,L-1,K) + AA(I,4,L,K)*DD(4,J,L-1,K)
          enddo
          H(J,L,K) = H(J,L,K) &
                   - AA(J,1,L,K)*RR(1,L-1,K) - AA(J,2,L,K)*RR(2,L-1,K) &
                   - AA(J,3,L,K)*RR(3,L-1,K) - AA(J,4,L,K)*RR(4,L-1,K)
       enddo

       do J = 1,M_
          do I = 1,M_
             E(I,J) = B(I,J,L,K)
          enddo
       enddo

       ! setup L & U matrices
       E(2,1) = E(2,1)/E(1,1)
       E(2,2) = E(2,2)-E(2,1)*E(1,2)
       E(2,3) = E(2,3)-E(2,1)*E(1,3)
       E(2,4) = E(2,4)-E(2,1)*E(1,4)
       E(3,1) = E(3,1)/E(1,1)
       E(3,2) = (E(3,2)-E(3,1)*E(1,2))/E(2,2)
       E(3,3) = E(3,3)-E(3,1)*E(1,3)-E(3,2)*E(2,3)
       E(3,4) = E(3,4)-E(3,1)*E(1,4)-E(3,2)*E(2,4)
       E(4,1) = E(4,1)/E(1,1)
       E(4,2) = (E(4,2)-E(4,1)*E(1,2))/E(2,2)
       E(4,3) = (E(4,3)-E(4,1)*E(1,3)-E(4,2)*E(2,3))/E(3,3)
       E(4,4) = E(4,4)-E(4,1)*E(1,4)-E(4,2)*E(2,4)-E(4,3)*E(3,4)
       ! invert L
       E(4,3) = -E(4,3)
       E(4,2) = -E(4,2)-E(4,3)*E(3,2)
       E(4,1) = -E(4,1)-E(4,2)*E(2,1)-E(4,3)*E(3,1)
       E(3,2) = -E(3,2)
       E(3,1) = -E(3,1)-E(3,2)*E(2,1)
       E(2,1) = -E(2,1)
       ! invert U
       E(4,4) = 1.e+0_fp/E(4,4)
       E(3,4) = -E(3,4)*E(4,4)/E(3,3)
       E(3,3) = 1.e+0_fp/E(3,3)
       E(2,4) = -(E(2,3)*E(3,4)+E(2,4)*E(4,4))/E(2,2)
       E(2,3) = -E(2,3)*E(3,3)/E(2,2)
       E(2,2) = 1.e+0_fp/E(2,2)
       E(1,4) = -(E(1,2)*E(2,4)+E(1,3)*E(3,4)+E(1,4)*E(4,4))/E(1,1)
       E(1,3) = -(E(1,2)*E(2,3)+E(1,3)*E(3,3))/E(1,1)
       E(1,2) = -E(1,2)*E(2,2)/E(1,1)
       E(1,1) = 1.e+0_fp/E(1,1)
       ! multiply U-invers * L-inverse
       E(1,1) = E(1,1)+E(1,2)*E(2,1)+E(1,3)*E(3,1)+E(1,4)*E(4,1)
       E(1,2) = E(1,2)+E(1,3)*E(3,2)+E(1,4)*E(4,2)
       E(1,3) = E(1,3)+E(1,4)*E(4,3)
       E(2,1) = E(2,2)*E(2,1)+E(2,3)*E(3,1)+E(2,4)*E(4,1)
       E(2,2) = E(2,2)+E(2,3)*E(3,2)+E(2,4)*E(4,2)
       E(2,3) = E(2,3)+E(2,4)*E(4,3)
       E(3,1) = E(3,3)*E(3,1)+E(3,4)*E(4,1)
       E(3,2) = E(3,3)*E(3,2)+E(3,4)*E(4,2)
       E(3,3) = E(3,3)+E(3,4)*E(4,3)
       E(4,1) = E(4,4)*E(4,1)
       E(4,2) = E(4,4)*E(4,2)
       E(4,3) = E(4,4)*E(4,3)

       do J = 1,M_
          RR(J,L,K) = E(J,1)*H(1,L,K)+E(J,2)*H(2,L,K) &
                    + E(J,3)*H(3,L,K)+E(J,4)*H(4,L,K)
       enddo

       ! BACK SOLUTION
       do L = ND-1,1,-1
          do J = 1,M_
             RR(J,L,K) = RR(J,L,K) &
                       + DD(J,1,L,K)*RR(1,L+1,K) + DD(J,2,L,K)*RR(2,L+1,K) &
                       + DD(J,3,L,K)*RR(3,L+1,K) + DD(J,4,L,K)*RR(4,L+1,K)
          enddo
       enddo

       ! mean J & H
       do L = 1,ND,2
          FJ(L,K) = RR(1,L,K)*WT(1) + RR(2,L,K)*WT(2) &
                  + RR(3,L,K)*WT(3) + RR(4,L,K)*WT(4)
       enddo
       do L = 2,ND,2
          FJ(L,K) = RR(1,L,K)*WT(1)*EMU(1) + RR(2,L,K)*WT(2)*EMU(2) &
                  + RR(3,L,K)*WT(3)*EMU(3) + RR(4,L,K)*WT(4)*EMU(4)
       enddo

       ! FJTOP = scaled diffuse flux out top-of-atmosphere (limit = mu0)
       ! FJBOT = scaled diffuse flux onto surface:
       ! ZFLUX = reflect/(1 + reflect) * mu0 * Fsolar(lower boundary)
       ! SUMBX = flux from Lambert reflected I+
       SUMT = RR(1, 1,K)*WT(1)*EMU(1) + RR(2, 1,K)*WT(2)*EMU(2) &
            + RR(3, 1,K)*WT(3)*EMU(3) + RR(4, 1,K)*WT(4)*EMU(4)
       SUMB = RR(1,ND,K)*WT(1)*EMU(1) + RR(2,ND,K)*WT(2)*EMU(2) &
            + RR(3,ND,K)*WT(3)*EMU(3) + RR(4,ND,K)*WT(4)*EMU(4)
       SUMBX = 4.e+0_fp*SUMB*RFL(K)/(1.0e+0_fp + RFL(K)) + ZFLUX(K)

       FJTOP(K) = 4.e+0_fp*SUMT
       FJBOT(K) = 4.e+0_fp*SUMB - SUMBX

    enddo

  END SUBROUTINE BLKSLV
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: gen_id
!
! !DESCRIPTION: Subroutine GEN generates coefficient matrices for the block
!  tri-diagonal system described in BLKSLV.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE GEN_ID(POMEGA,FZ,ZTAU,ZFLUX,RFL,PM,PM0,B,CC,AA,A,H,C,ND)
!
! !USES:
! 
    USE CMN_FJX_Mod, ONLY : W_, WT, EMU, M_, M2_, N_
!
! !INPUT PARAMETERS:
!
    INTEGER,  INTENT(IN)  :: ND
    REAL(fp), INTENT(IN)  :: POMEGA(M2_,N_)
    REAL(fp), INTENT(IN)  :: PM(M_,M2_)
    REAL(fp), INTENT(IN)  :: PM0(M2_)
    REAL(fp), INTENT(IN)  :: ZFLUX,RFL
    REAL(fp), INTENT(IN), DIMENSION(N_) :: FZ,ZTAU
!
! !OUTPUT PARAMETERS:
!
    REAL(fp), INTENT(OUT),DIMENSION(M_,M_,N_) :: B,AA,CC
    REAL(fp), INTENT(OUT),DIMENSION(M_,N_)    :: A,C,H

!
! !REVISION HISTORY:
!  28 Mar 2013 - S. D. Eastham - Copied from Fast-JX v7.0
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER I, J, K, L1,L2,LL
    REAL(fp)  SUM0, SUM1, SUM2, SUM3
    REAL(fp)  DELTAU, D1, D2, SURFAC

    REAL(fp), DIMENSION(M_,M_) :: S,T,U,V,W

    !=================================================================
    ! GEN_ID begins here!
    !=================================================================

    ! upper boundary:  2nd-order terms
    L1 = 1
    L2 = 2
    do I = 1,M_
       SUM0 = POMEGA(1,L1)*PM(I,1)*PM0(1) + POMEGA(3,L1)*PM(I,3)*PM0(3) &
            + POMEGA(5,L1)*PM(I,5)*PM0(5) + POMEGA(7,L1)*PM(I,7)*PM0(7)
       SUM2 = POMEGA(1,L2)*PM(I,1)*PM0(1) + POMEGA(3,L2)*PM(I,3)*PM0(3) &
            + POMEGA(5,L2)*PM(I,5)*PM0(5) + POMEGA(7,L2)*PM(I,7)*PM0(7)
       SUM1 = POMEGA(2,L1)*PM(I,2)*PM0(2) + POMEGA(4,L1)*PM(I,4)*PM0(4) &
            + POMEGA(6,L1)*PM(I,6)*PM0(6) + POMEGA(8,L1)*PM(I,8)*PM0(8)
       SUM3 = POMEGA(2,L2)*PM(I,2)*PM0(2) + POMEGA(4,L2)*PM(I,4)*PM0(4) &
            + POMEGA(6,L2)*PM(I,6)*PM0(6) + POMEGA(8,L2)*PM(I,8)*PM0(8)
       H(I,L1) = 0.5e+0_fp*(SUM0*FZ(L1) + SUM2*FZ(L2))
       A(I,L1) = 0.5e+0_fp*(SUM1*FZ(L1) + SUM3*FZ(L2))
    enddo

    do I = 1,M_
    do J = 1,I
       SUM0 = POMEGA(1,L1)*PM(I,1)*PM(J,1) + POMEGA(3,L1)*PM(I,3)*PM(J,3) &
            + POMEGA(5,L1)*PM(I,5)*PM(J,5) + POMEGA(7,L1)*PM(I,7)*PM(J,7)
       SUM2 = POMEGA(1,L2)*PM(I,1)*PM(J,1) + POMEGA(3,L2)*PM(I,3)*PM(J,3) &
            + POMEGA(5,L2)*PM(I,5)*PM(J,5) + POMEGA(7,L2)*PM(I,7)*PM(J,7)
       SUM1 = POMEGA(2,L1)*PM(I,2)*PM(J,2) + POMEGA(4,L1)*PM(I,4)*PM(J,4) &
            + POMEGA(6,L1)*PM(I,6)*PM(J,6) + POMEGA(8,L1)*PM(I,8)*PM(J,8)
       SUM3 = POMEGA(2,L2)*PM(I,2)*PM(J,2) + POMEGA(4,L2)*PM(I,4)*PM(J,4) &
            + POMEGA(6,L2)*PM(I,6)*PM(J,6) + POMEGA(8,L2)*PM(I,8)*PM(J,8)
       S(I,J) = - SUM2*WT(J)
       S(J,I) = - SUM2*WT(I)
       T(I,J) = - SUM1*WT(J)
       T(J,I) = - SUM1*WT(I)
       V(I,J) = - SUM3*WT(J)
       V(J,I) = - SUM3*WT(I)
       B(I,J,L1) = - 0.5e+0_fp*(SUM0 + SUM2)*WT(J)
       B(J,I,L1) = - 0.5e+0_fp*(SUM0 + SUM2)*WT(I)
    enddo
    enddo

    do I = 1,M_
       S(I,I)   = S(I,I)    + 1.0e+0_fp
       T(I,I)   = T(I,I)    + 1.0e+0_fp
       V(I,I)   = V(I,I)    + 1.0e+0_fp
       B(I,I,L1)= B(I,I,L1) + 1.0e+0_fp

       C(I,L1)= S(I,1)*A(1,L1)/EMU(1) + S(I,2)*A(2,L1)/EMU(2) &
              + S(I,3)*A(3,L1)/EMU(3) + S(I,4)*A(4,L1)/EMU(4)
    enddo

    do I = 1,M_
    do J = 1,M_
       W(J,I) = S(J,1)*T(1,I)/EMU(1) + S(J,2)*T(2,I)/EMU(2) &
              + S(J,3)*T(3,I)/EMU(3) + S(J,4)*T(4,I)/EMU(4)
       U(J,I) = S(J,1)*V(1,I)/EMU(1) + S(J,2)*V(2,I)/EMU(2) &
              + S(J,3)*V(3,I)/EMU(3) + S(J,4)*V(4,I)/EMU(4)
    enddo
    enddo
    ! upper boundary, 2nd-order, C-matrix is full (CC)
    DELTAU = ZTAU(L2) - ZTAU(L1)
    D2 = 0.25e+0_fp*DELTAU
    do I = 1,M_
       do J = 1,M_
          B(I,J,L1) = B(I,J,L1) + D2*W(I,J)
          CC(I,J,L1) = D2*U(I,J)
       enddo
       H(I,L1) = H(I,L1) + 2.0e+0_fp*D2*C(I,L1)
       A(I,L1) = 0.0e+0_fp
    enddo
    do I = 1,M_
       D1 = EMU(I)/DELTAU
       B(I,I,L1)  = B(I,I,L1) + D1
       CC(I,I,L1) = CC(I,I,L1) - D1
    enddo

    ! intermediate points:  can be even or odd, A & C diagonal
    ! mid-layer h-points, Legendre terms 2,4,6,8
    do LL=2,ND-1,2
       DELTAU = ZTAU(LL+1) - ZTAU(LL-1)
       do I = 1,M_
          A(I,LL) = EMU(I)/DELTAU
          C(I,LL) = -A(I,LL)
          H(I,LL) = FZ(LL)*( &
               POMEGA(2,LL)*PM(I,2)*PM0(2) + POMEGA(4,LL)*PM(I,4)*PM0(4) &
             + POMEGA(6,LL)*PM(I,6)*PM0(6) + POMEGA(8,LL)*PM(I,8)*PM0(8))
       enddo
       do I = 1,M_
       do J=1,I
          SUM0 = POMEGA(2,LL)*PM(I,2)*PM(J,2) + POMEGA(4,LL)*PM(I,4)*PM(J,4) &
               + POMEGA(6,LL)*PM(I,6)*PM(J,6) + POMEGA(8,LL)*PM(I,8)*PM(J,8)
          B(I,J,LL) =  - SUM0*WT(J)
          B(J,I,LL) =  - SUM0*WT(I)
       enddo
       enddo
       do I = 1,M_
          B(I,I,LL) = B(I,I,LL) + 1.0e+0_fp
       enddo
    enddo

    ! odd-layer j-points, Legendre terms 1,3,5,7
    do LL=3,ND-2,2
       DELTAU = ZTAU(LL+1) - ZTAU(LL-1)
       do I = 1,M_
          A(I,LL) = EMU(I)/DELTAU
          C(I,LL) = -A(I,LL)
          H(I,LL) = FZ(LL)*( &
               POMEGA(1,LL)*PM(I,1)*PM0(1) + POMEGA(3,LL)*PM(I,3)*PM0(3) &
             + POMEGA(5,LL)*PM(I,5)*PM0(5) + POMEGA(7,LL)*PM(I,7)*PM0(7))
       enddo
       do I = 1,M_
       do J=1,I
          SUM0 = POMEGA(1,LL)*PM(I,1)*PM(J,1) + POMEGA(3,LL)*PM(I,3)*PM(J,3) &
               + POMEGA(5,LL)*PM(I,5)*PM(J,5) + POMEGA(7,LL)*PM(I,7)*PM(J,7)
          B(I,J,LL) =  - SUM0*WT(J)
          B(J,I,LL) =  - SUM0*WT(I)
       enddo
       enddo
       do I = 1,M_
          B(I,I,LL) = B(I,I,LL) + 1.0e+0_fp
       enddo
    enddo

    ! lower boundary:  2nd-order terms
    L1 = ND
    L2 = ND-1
    do I = 1,M_
       SUM0 = POMEGA(1,L1)*PM(I,1)*PM0(1) + POMEGA(3,L1)*PM(I,3)*PM0(3) &
            + POMEGA(5,L1)*PM(I,5)*PM0(5) + POMEGA(7,L1)*PM(I,7)*PM0(7)
       SUM2 = POMEGA(1,L2)*PM(I,1)*PM0(1) + POMEGA(3,L2)*PM(I,3)*PM0(3) &
            + POMEGA(5,L2)*PM(I,5)*PM0(5) + POMEGA(7,L2)*PM(I,7)*PM0(7)
       SUM1 = POMEGA(2,L1)*PM(I,2)*PM0(2) + POMEGA(4,L1)*PM(I,4)*PM0(4) &
            + POMEGA(6,L1)*PM(I,6)*PM0(6) + POMEGA(8,L1)*PM(I,8)*PM0(8)
       SUM3 = POMEGA(2,L2)*PM(I,2)*PM0(2) + POMEGA(4,L2)*PM(I,4)*PM0(4) &
            + POMEGA(6,L2)*PM(I,6)*PM0(6) + POMEGA(8,L2)*PM(I,8)*PM0(8)
       H(I,L1) = 0.5e+0_fp*(SUM0*FZ(L1) + SUM2*FZ(L2))
       A(I,L1) = 0.5e+0_fp*(SUM1*FZ(L1) + SUM3*FZ(L2))
    enddo

    do I = 1,M_
    do J = 1,I
       SUM0 = POMEGA(1,L1)*PM(I,1)*PM(J,1) + POMEGA(3,L1)*PM(I,3)*PM(J,3) &
            + POMEGA(5,L1)*PM(I,5)*PM(J,5) + POMEGA(7,L1)*PM(I,7)*PM(J,7)
       SUM2 = POMEGA(1,L2)*PM(I,1)*PM(J,1) + POMEGA(3,L2)*PM(I,3)*PM(J,3) &
            + POMEGA(5,L2)*PM(I,5)*PM(J,5) + POMEGA(7,L2)*PM(I,7)*PM(J,7)
       SUM1 = POMEGA(2,L1)*PM(I,2)*PM(J,2) + POMEGA(4,L1)*PM(I,4)*PM(J,4) &
            + POMEGA(6,L1)*PM(I,6)*PM(J,6) + POMEGA(8,L1)*PM(I,8)*PM(J,8)
       SUM3 = POMEGA(2,L2)*PM(I,2)*PM(J,2) + POMEGA(4,L2)*PM(I,4)*PM(J,4) &
            + POMEGA(6,L2)*PM(I,6)*PM(J,6) + POMEGA(8,L2)*PM(I,8)*PM(J,8)
       S(I,J) = - SUM2*WT(J)
       S(J,I) = - SUM2*WT(I)
       T(I,J) = - SUM1*WT(J)
       T(J,I) = - SUM1*WT(I)
       V(I,J) = - SUM3*WT(J)
       V(J,I) = - SUM3*WT(I)
       B(I,J,L1) = - 0.5e+0_fp*(SUM0 + SUM2)*WT(J)
       B(J,I,L1) = - 0.5e+0_fp*(SUM0 + SUM2)*WT(I)
    enddo
    enddo

    do I = 1,M_
       S(I,I)   = S(I,I)   + 1.0e+0_fp
       T(I,I)   = T(I,I)   + 1.0e+0_fp
       V(I,I)   = V(I,I)   + 1.0e+0_fp
       B(I,I,L1)= B(I,I,L1) + 1.0e+0_fp

       C(I,L1)= S(I,1)*A(1,L1)/EMU(1) + S(I,2)*A(2,L1)/EMU(2) &
              + S(I,3)*A(3,L1)/EMU(3) + S(I,4)*A(4,L1)/EMU(4)
    enddo

    do I = 1,M_
    do J = 1,M_
       W(J,I) = S(J,1)*T(1,I)/EMU(1) + S(J,2)*T(2,I)/EMU(2) &
              + S(J,3)*T(3,I)/EMU(3) + S(J,4)*T(4,I)/EMU(4)
       U(J,I) = S(J,1)*V(1,I)/EMU(1) + S(J,2)*V(2,I)/EMU(2) &
              + S(J,3)*V(3,I)/EMU(3) + S(J,4)*V(4,I)/EMU(4)
    enddo
    enddo

    ! lower boundary, 2nd-order, A-matrix is full (AA)
    DELTAU = ZTAU(L1) - ZTAU(L2)
    D2 = 0.25e+0_fp*DELTAU
    SURFAC = 4.0e+0_fp*RFL/(1.0e+0_fp + RFL)
    do I = 1,M_
       D1 = EMU(I)/DELTAU
       SUM0 = D1 + D2*(W(I,1)+W(I,2)+W(I,3)+W(I,4))
       SUM1 = SURFAC*SUM0
       do J = 1,M_
          AA(I,J,L1) = - D2*U(I,J)
          B(I,J,L1) = B(I,J,L1) + D2*W(I,J) - SUM1*EMU(J)*WT(J)
       enddo
       H(I,L1) = H(I,L1) - 2.0e+0_fp*D2*C(I,L1) + SUM0*ZFLUX
    enddo

    do I = 1,M_
       D1 = EMU(I)/DELTAU
       AA(I,I,L1) = AA(I,I,L1) + D1
       B(I,I,L1)  = B(I,I,L1) + D1
       C(I,L1) = 0.0e+0_fp
    enddo

  END SUBROUTINE GEN_ID
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: jratet
!
! !DESCRIPTION: Subroutine JRATET calculates temperature-dependent J-rates.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE JRATET(PPJ,TTJ,FFF,VALJL,LCTM,LCHEM,NJXU)
!
! !USES:
!
  USE CMN_FJX_Mod, ONLY : JXL1_, W_, NJX, X_, LQQ, QQQ, SQQ, TQQ
  USE CMN_FJX_Mod, ONLY : QO2, QO3, Q1D
!
! !INPUT PARAMETERS:
!
    integer,  intent(in)    :: LCTM,LCHEM,NJXU
    real(fp), intent(in)    ::  PPJ(JXL1_+1),TTJ(JXL1_+1)
!
! !INPUT/OUTPUT PARAMETERS:
!
    real(fp), intent(inout) ::  FFF(W_,LCTM)
!
! !OUTPUT VARIABLES:
!
    real(fp), intent(out), dimension(LCTM,NJXU) ::  VALJL
!
! !REMARKS:
!
! !REVISION HISTORY:
!  28 Mar 2013 - S. D. Eastham - Copied from Fast-JX v7.0
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    real(fp)  VALJ(X_)
    real(fp)  QO2TOT, QO3TOT, QO31DY, QO31D, QQQT, TFACT
    real(fp)  TT,PP,DD,TT200,TFACA,TFAC0,TFAC1,TFAC2
    real(fp)  QQQA,QQ2,QQ1A,QQ1B
    integer   J,K,L, IV

    !=================================================================
    ! JRATET begins here!
    !=================================================================

    if (NJXU .lt. NJX) then
       call EXITC(' JRATET:  CTM has not enough J-values dimensioned')
    endif
    do L = 1,LCTM
       ! need temperature, pressure, and density at mid-layer
       ! (for some quantum yields):
       TT   = TTJ(L)
       if (L .eq. 1) then
          PP = PPJ(1)
       else
          PP  = (PPJ(L)+PPJ(L+1))*0.5e+0_fp
       endif
       DD = 7.24e18*PP/TT

       ! if W_=18/12, must zero bin-11/5 below 100 hPa, since O2 e-fold is
       ! too weak and does not represent the decay of 215.5-221.5 nm sunlight.
       if (PP .gt. 100.e+0_fp) then
          if (W_ .eq. 18) then
             FFF(11,L) = 0.e+0_fp
          elseif (W_ .eq. 12) then
             FFF(5,L) = 0.e+0_fp
          endif
       endif

       do J = 1,NJXU
          VALJ(J) = 0.e+0_fp
       enddo

       do K = 1,W_
          call X_interp (TT,QO2TOT, TQQ(1,1),QO2(K,1), &
                         TQQ(2,1),QO2(K,2), TQQ(3,1),QO2(K,3), LQQ(1))
          call X_interp (TT,QO3TOT, TQQ(1,2),QO3(K,1), &
                         TQQ(2,2),QO3(K,2), TQQ(3,2),QO3(K,3), LQQ(2))
          call X_interp (TT,QO31DY, TQQ(1,3),Q1D(K,1), &
                         TQQ(2,3),Q1D(K,2), TQQ(3,3),Q1D(K,3), LQQ(3))
          QO31D  = QO31DY*QO3TOT
          VALJ(1) = VALJ(1) + QO2TOT*FFF(K,L)
          VALJ(2) = VALJ(2) + QO3TOT*FFF(K,L)
          VALJ(3) = VALJ(3) + QO31D*FFF(K,L)
       enddo

       do J = 4,NJXU
       do K = 1,W_
          ! also need to allow for Pressure interpolation if SQQ(J) = 'p'
          if (SQQ(J) .eq.'p') then
             call X_interp (PP,QQQT, TQQ(1,J),QQQ(K,1,J), &
                            TQQ(2,J),QQQ(K,2,J), TQQ(3,J),QQQ(K,3,J), LQQ(J))
          else
             call X_interp (TT,QQQT, TQQ(1,J),QQQ(K,1,J), &
                            TQQ(2,J),QQQ(K,2,J), TQQ(3,J),QQQ(K,3,J), LQQ(J))
          endif
          VALJ(J) = VALJ(J) + QQQT*FFF(K,L)
       enddo
       enddo

       do J=1,NJXU
          VALJL(L,J) = VALJ(J)
       enddo

    enddo

    ! Zero non-chemistry layers
    if (LCHEM.lt.LCTM) then
       do L=(LCTM+1),LCHEM
       do J=1,NJXU
          VALJL(L,J) = 0.e+0_fp
       enddo
       enddo
    endif

  END SUBROUTINE JRATET
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: x_interp
!
! !DESCRIPTION: Subroutine X\_INTERP is an up-to-three-point linear interp.
!  function for cross-sections.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE X_INTERP (TINT,XINT,T1,X1,T2,X2,T3,X3,L123)
!
! !USES:
!
!
! !INPUT PARAMETERS:
!
    REAL(fp), INTENT(IN)  ::  TINT,T1,T2,T3, X1,X2,X3
    INTEGER,  INTENT(IN)  ::  L123
!
! !OUTPUT VARIABLES:
!
    REAL(fp), INTENT(OUT) ::  XINT
!
! !REMARKS:
!
! !REVISION HISTORY:
!  28 Mar 2013 - S. D. Eastham - Copied from Fast-JX v7.0
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    REAL(fp)  TFACT

    !=================================================================
    ! X_INTERP begins here!
    !=================================================================

    if (L123 .le. 1) then
       XINT = X1
    elseif (L123 .eq. 2) then
       TFACT = max(0.e+0_fp,min(1.e+0_fp,(TINT-T1)/(T2-T1) ))
       XINT = X1 + TFACT*(X2 - X1)
    else
       if (TINT.le. T2) then
          TFACT = max(0.e+0_fp,min(1.e+0_fp,(TINT-T1)/(T2-T1) ))
          XINT = X1 + TFACT*(X2 - X1)
       else
          TFACT = max(0.e+0_fp,min(1.e+0_fp,(TINT-T2)/(T3-T2) ))
          XINT = X2 + TFACT*(X3 - X2)
       endif
    endif

  END SUBROUTINE X_INTERP
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: sphere2
!
! !DESCRIPTION: Subroutine SPHERE2 is an AMF2.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE SPHERE2 (U0,ZHL,AMF2,L1U,LJX1U)
!
! !USES:
!
    USE CMN_FJX_Mod, ONLY: RAD
!
! !INPUT PARAMETERS:
!
    INTEGER,  INTENT(IN)  ::   L1U, LJX1U
    REAL(fp), INTENT(IN)  ::   U0,ZHL(L1U+1)
!
! !OUTPUT VARIABLES:
!
    REAL(fp), INTENT(OUT) ::   AMF2(2*LJX1U+1,2*LJX1U+1)
!
! !REMARKS:
! Quoting from the original:
!  New v6.2: does AirMassFactors for mid-layer, needed for SZA ~ 90
!  This new AMF2 does each of the half-layers of the CTM separately,
!     whereas the original, based on the pratmo code did the whole layers
!     and thus calculated the ray-path to the CTM layre edges, NOT the middle.
!  Since fast-JX is meant to calculate the intensity at the mid-layer, the
!     solar beam at low sun (interpolated between layer edges) was incorrect.
!  This new model does make some approximations of the geometry of the layers:
!     the CTM layer is split evenly in mass (good) and in height (approx).
!                                                                             .
!  Calculation of spherical geometry; derive tangent heights, slant path
!  lengths and air mass factor for each layer. Not called when
!  SZA > 98 degrees.  Beyond 90 degrees, include treatment of emergent
!  beam (where tangent height is below altitude J-value desired at).
!                                                                             .
!  ---------------------------------------------------------------------
!  Inputs:
!     U0      cos(solar zenith angle)
!     RAD  radius of Earth mean sea level (cm)
!     ZHL(L)  height (cm) of the bottom edge of CTM level L
!     ZZHT    scale height (cm) used above top of CTM (ZHL(L_+1))
!     L1U     dimension of CTM = levels +1 (L+1 = above-CTM level)
!  Outputs:
!     AMF2(I,J) = air mass factor for CTM level I for sunlight reaching J
!         ( these are calculated for both layer middle and layer edge)
!  ---------------------------------------------------------------------
!
! !REVISION HISTORY:
!  28 Mar 2013 - S. D. Eastham - Copied from Fast-JX v7.0
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER, PARAMETER  ::  LSPH_ = 200

    ! RZ      Distance from centre of Earth to each point (cm)
    ! RQ      Square of radius ratios
    ! SHADHT  Shadow height for the current SZA
    ! XL      Slant path between points
    INTEGER  I, J, K, II, L2
    REAL(fp)   XMU1,XMU2,XL,DIFF,SHADHT,RZ(LSPH_+1)
    REAL(fp)   RZ2(2*LSPH_+1),RQ2(2*LSPH_+1)

    !=================================================================
    ! SPHERE2 begins here!
    !=================================================================

    ! must have top-of-atmos (NOT top-of-CTM) defined
    !      ZHL(L1U+1) = ZHL(L1U) + ZZHT

    if (L1U .gt. LSPH_) then
       call EXITC(' SPHERE2: temp arrays not large enough')
    endif

    RZ(1) = RAD + ZHL(1)
    do II = 2,L1U+1
       RZ(II)   = RAD + ZHL(II)
    enddo

    ! calculate heights for edges of split CTM-layers
    L2 = 2*L1U
    do II = 2,L2,2
       I = II/2
       RZ2(II-1) = RZ(I)
       RZ2(II) = 0.5e+0_fp*(RZ(I)+RZ(I+1))
    enddo
    RZ2(L2+1) = RZ(L1U+1)
    do II = 1,L2
       RQ2(II) = (RZ2(II)/RZ2(II+1))**2
    enddo

    ! shadow height for SZA > 90
    if (U0 .lt. 0.0e+0_fp)  then
       SHADHT = RZ2(1)/sqrt(1.0e+0_fp - U0**2)
    else
       SHADHT = 0.e+0_fp
    endif

    ! up from the surface calculating the slant paths between each level
    ! and the level above, and deriving the appropriate Air Mass Factor
    AMF2(:,:) = 0.e+0_fp

    do 16 J = 1,2*L1U+1

       ! Air Mass Factors all zero if below the tangent height
       if (RZ2(J) .lt. SHADHT) goto 16

       ! Ascend from layer J calculating AMF2s
       XMU1 = abs(U0)
       do I = J,2*L1U
          XMU2     = sqrt(1.0e+0_fp - RQ2(I)*(1.0e+0_fp-XMU1**2))
          XL       = RZ2(I+1)*XMU2 - RZ2(I)*XMU1
          AMF2(I,J) = XL / (RZ2(I+1)-RZ2(I))
          XMU1     = XMU2
       enddo

       ! fix above top-of-atmos (L=L1U+1), must set DTAU(L1U+1)=0
       AMF2(2*L1U+1,J) = 1.e+0_fp

       ! Twilight case - Emergent Beam, calc air mass factors below layer
       if (U0 .ge. 0.0e+0_fp) goto 16

       ! Descend from layer J
       XMU1       = abs(U0)
       do II = J-1,1,-1
          DIFF        = RZ2(II+1)*sqrt(1.0e+0_fp-XMU1**2)-RZ2(II)
          if (II.eq.1)  DIFF = max(DIFF,0.e+0_fp)   ! filter

          ! Tangent height below current level - beam passes through twice
          if (DIFF .lt. 0.0e+0_fp)  then
             XMU2      = sqrt(1.0e+0_fp - (1.0e+0_fp-XMU1**2)/RQ2(II))
             XL        = abs(RZ2(II+1)*XMU1-RZ2(II)*XMU2)
             AMF2(II,J) = 2.e+0_fp*XL/(RZ2(II+1)-RZ2(II))
             XMU1      = XMU2

          ! Lowest level intersected by emergent beam
          else
             XL        = RZ2(II+1)*XMU1*2.0e+0_fp
             AMF2(II,J) = XL/(RZ2(II+1)-RZ2(II))
             goto 16
          endif
       enddo

16  continue

  END SUBROUTINE SPHERE2
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: extral
!
! !DESCRIPTION: Subroutine EXTRAL adds sub-layers to thick cloud/aerosol layers
!  using log-spacing for sub-layers of increasing thickness ATAU.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE EXTRAL (Input_Opt,State_Diag,DTAUX,L1X,L2X,NX,JXTRA,ILON,ILAT)
!
! !USES:
!
    USE CMN_FJX_Mod, ONLY : ATAU, ATAU0, JTAUMX
    USE Input_Opt_Mod,  ONLY : OptInput
    USE State_Diag_Mod, ONLY : DgnState
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN) :: Input_Opt   ! Input options
    INTEGER,        INTENT(IN) :: L1X,L2X     !index of cloud/aerosol
    integer,        intent(in) :: NX          !Mie scattering array size
    real(fp),       intent(in) :: DTAUX(L1X)  !cloud+3aerosol OD in each layer
    integer,        intent(in) :: ILON, ILAT  !lon,lat index
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(DgnState), INTENT(INOUT) :: State_Diag  ! Diagnostics State object
!
! !OUTPUT VARIABLES:
!
    integer,        intent(out):: JXTRA(L2X+1)!number of sub-layers to be added
!
! !REMARKS:
!     DTAUX(L=1:L1X) = Optical Depth in layer L (generally 600 nm OD)
!        This can be just cloud or cloud+aerosol, it is used only to set
!        the number in levels to insert in each layer L
!        Set for log-spacing of tau levels, increasing top-down.
!                                                                             .
!     N.B. the TTAU, etc calculated here are NOT used elsewhere
!                                                                             .
!   The log-spacing parameters have been tested for convergence and chosen
!     to be within 0.5% for ranges OD=1-500, rflect=0-100%, mu0=0.1-1.0
!     use of ATAU = 1.18 and min = 0.01, gives at most +135 pts for OD=100
!     ATAU = 1.12 now recommended for more -accurate heating rates (not J's)
!
! !REVISION HISTORY:
!  28 Mar 2013 - S. D. Eastham - Copied from Fast-JX v7.0
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER JTOTL,I,L,L2
    REAL(fp)  TTAU(L2X+1),DTAUJ, ATAU1,ATAULN,ATAUM,ATAUN1

#ifdef MODEL_GEOS
    ! ckeller, 5/21/18
    LOGICAL  :: failed
    INTEGER  :: N, NMAX
    REAL(fp) :: ATAULOC
#endif

    !=================================================================
    ! EXTRAL begins here!
    !=================================================================

#ifdef MODEL_GEOS
    ! This routine now repeats the extra layer computation with an
    ! increased heating rate (ATAU) if there is an array overfloat.
    ! This is repeated maximum 5 times. The diagnostics arrays
    ! EXTRAL_NLEVS and EXTRAL_NITER archive the number of extra layers
    ! and the number of iterations needed to converge to that solution.
    ! Ideally, NITER is 1 and no adjustments to the heating rate are
    ! needed (ckeller, 5/22/18).
    NMAX = MAX(1,Input_Opt%FJX_EXTRAL_ITERMAX)
    DO N=1,NMAX
       ! local heating rate
       ATAULOC = ATAU + (0.06*(N-1))
#endif

       ! Reinitialize arrays
       TTAU(:)  = 0.e+0_fp
       JXTRA(:) = 0

       ! combine these edge- and mid-layer points into grid of size:
       !     L2X+1 = 2*L1X+1 = 2*L_+3
       ! calculate column optical depths above each level, TTAU(1:L2X+1)
       !     note that TTAU(L2X+1)=0 and TTAU(1)=total OD
       !
       ! Divide thick layers to achieve better accuracy in the scattering code
       ! In the original fast-J, equal sub-layers were chosen, this is wasteful
       ! and this new code (ver 5.3) uses log-scale:
       !     Each succesive layer (down) increase thickness by ATAU > 1
       !     e.g., if ATAU = 2, a layer with OD = 15 could be divided into
       !     4 sub-layers with ODs = 1 - 2 - 4 - 8
       ! The key parameters are:
       !     ATAU = factor increase from one layer to the next
       !     ATAUMN = the smallest OD layer desired
       !     JTAUMX = maximum number of divisions (i.e., may not get to ATAUMN)
       ! These are set in CMN_FJX_MOD, and have been tested/optimized

#if defined( MODEL_GEOS )
       ATAU1  = ATAULOC - 1.e+0_fp
       ATAULN = log(ATAULOC)
#else
       ATAU1  = ATAU - 1.e+0_fp
       ATAULN = log(ATAU)
#endif
       TTAU(L2X+1)  = 0.0e+0_fp

       do L2 = L2X,1,-1
          L         = (L2+1)/2
          DTAUJ     = 0.5e+0_fp * DTAUX(L)
          TTAU(L2)  = TTAU(L2+1) + DTAUJ
          ! Now compute the number of log-spaced sub-layers to be added in
          ! the interval TTAU(L2) > TTAU(L2+1)
          ! The objective is to have successive TAU-layers increasing by factor
          ! ATAU >1 the number of sub-layers + 1
          if (TTAU(L2) .lt. ATAU0) then
             JXTRA(L2) = 0
          else
             ATAUM    = max(ATAU0, TTAU(L2+1))
             ATAUN1 = log(TTAU(L2)/ATAUM) / ATAULN
             JXTRA(L2) = min(JTAUMX, max(0, int(ATAUN1 - 0.5e+0_fp)))
          endif
       enddo

       ! check on overflow of arrays, cut off JXTRA at lower L if too many
       ! levels
#ifdef MODEL_GEOS
       failed   = .FALSE.
       JTOTL    = L2X + 2
       do L2 = L2X,1,-1
          JTOTL  = JTOTL + JXTRA(L2)
          if (JTOTL .gt. NX/2)  then
             failed = .TRUE.
             exit
          endif
       enddo

       ! exit loop if not failed
       if ( .not. failed ) exit
    enddo

    ! print error and cut off JXTRAL at lower L if too many levels
    if ( failed ) then
       IF ( Input_Opt%FJX_EXTRAL_ERR ) THEN
          write(6,'(A,7I5)') 'N_/L2_/L2-cutoff JXTRA:',ILON,ILAT,NX,L2X,L2,JXTRA(L2),JTOTL
       ENDIF
       do L = L2,1,-1
          JXTRA(L) = 0
       enddo
       !go to 10
    endif
    !enddo
    !10 continue

    ! Fill diagnostics arrays
    IF ( State_Diag%Archive_EXTRALNLEVS ) THEN
       State_Diag%EXTRALNLEVS(ILON,ILAT) = SUM(JXTRA(:))
    ENDIF
    IF ( State_Diag%Archive_EXTRALNITER ) THEN
       State_Diag%EXTRALNITER(ILON,ILAT) = N
    ENDIF
#else
    JTOTL    = L2X + 2
    do L2 = L2X,1,-1
       JTOTL  = JTOTL + JXTRA(L2)
       if (JTOTL .gt. NX/2)  then
          write(6,'(A,7I5)') 'N_/L2_/L2-cutoff JXTRA:',ILON,ILAT,NX,L2X,L2,JXTRA(L2),JTOTL
          do L = L2,1,-1
             JXTRA(L) = 0
          enddo
          go to 10
       endif
    enddo
10  continue
#endif

  END SUBROUTINE EXTRAL
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: exitc
!
! !DESCRIPTION: Subroutine EXITC forces an error in GEOS-Chem and quits.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE EXITC (T_EXIT)
!
! !USES:
!
    USE ERROR_MOD, ONLY : ERROR_STOP
!
! !INPUT PARAMETERS:
!
    CHARACTER(LEN=*), INTENT(IN) ::  T_EXIT
!
! !REVISION HISTORY:
!  28 Mar 2013 - S. D. Eastham - Copied from Fast-JX v7.0
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    CALL ERROR_STOP( T_EXIT, 'fast_jx_mod.F90' )

  END SUBROUTINE EXITC
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: opmie
!
! !DESCRIPTION: Subroutine OPMIE is a core Fast-JX scattering subroutine,
!  specifically for Mie scattering.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE OPMIE(DTAUX,POMEGAX,U0,RFL,AMF2,JXTRA, &
                   FJACT,FJTOP,FJBOT,FSBOT,FJFLX,FLXD,FLXD0,LU)
!
! !USES:
!
    USE CMN_FJX_Mod, ONLY : ATAU0, JXL_, JXL1_, JXL2_, M2_, N_, W_
!
! !INPUT PARAMETERS:
!
    REAL(fp), INTENT(IN)  ::  DTAUX(JXL1_,W_),POMEGAX(8,JXL1_,W_)
    REAL(fp), INTENT(IN)  ::  AMF2(2*JXL1_+1,2*JXL1_+1)
    REAL(fp), INTENT(IN)  ::  U0,RFL(W_)
    INTEGER,  INTENT(IN)  ::  JXTRA(JXL2_+1), LU
!
! !OUTPUT VARIABLES:
!
    REAL(fp), INTENT(OUT) ::  FJACT(JXL_,W_),FJTOP(W_)
    REAL(fp), INTENT(OUT) ::  FJBOT(W_),FSBOT(W_)
    REAL(fp), INTENT(OUT) ::  FJFLX(JXL_,W_),FLXD(JXL1_,W_),FLXD0(W_)
!
! !REVISION HISTORY:
!  28 Mar 2013 - S. D. Eastham - Copied from Fast-JX v7.0
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER JNDLEV(JXL_),JNELEV(JXL1_)
    INTEGER JADDLV(JXL2_+1),JADDTO(JXL2_+1),L2LEV(JXL2_+1)
    INTEGER JTOTL,I,II,J,K,L,LL,IX,JK,   L2,L2L,L22,LZ,LZZ,ND
    INTEGER L1U,L2U,   LZ0,LZ1,LZMID
    REAL(fp)   SUMT,SUMJ

    REAL(fp)  DTAU(JXL1_+1,W_),POMEGAJ(M2_,JXL2_+1,W_)
    REAL(fp)  TTAU(JXL2_+1,W_)
    REAL(fp)  FTAU2(JXL2_+1,W_),POMEGAB(M2_,W_)
    REAL(fp)  ATAUA,ATAUZ,XLTAU,TAUDN,TAUUP,DTAUJ,FJFLX0
    REAL(fp), DIMENSION(W_) :: TAUBTM,TAUTOP,FBTM,FTOP,ZFLUX

    !  variables used in mie code-----------------------------------------
    REAL(fp), DIMENSION(W_)         :: FJT,FJB
    REAL(fp), DIMENSION(N_,W_)      :: FJ,FZ,ZTAU
    REAL(fp), DIMENSION(M2_,N_,W_)  :: POMEGA
    REAL(fp), DIMENSION(2*JXL1_,W_)  :: FLXD2

    !=================================================================
    ! OPMIE begins here!
    !=================================================================

    ! there is a parallel correspondence:
    !  dimension of JX arrays JXL_ .ge. dimension that CTM is using = L_
    !  but calculation is done for L_=LU, L1_=L1U, L2_=L2U lengths of CTM
    !
    !  fast-J Mie code for J_s, only uses 8-term expansion, 4-Gauss pts
    !
    ! in:
    !     DTAUX(1:L1_,1:W_) = optical depth of each layer
    !     POMEGAX(1:8,1:L1_,1:W_) = scattering phase fn (multiplied by s-s albedo)
    !     U0  = cos (SZA)
    !     RFL(1:W_) = Lambertian albedo of surface
    !     AMF2(1:2*L1_+1,1:2*L1_+1) = air mass factor (I,L)=wt of layer-I to layer-L
    !        AMF2 now does both edges and middle of CTM layers
    !     JXTRA(1:L1_) = number 0:J = no. of additional levels to be inserted
    ! out:
    !     FJACT(1:L_,1:W_) = mean actinic flux(diff+direct) at std CTM levels(mid-lyr)
    !  (new ver 5.7 diagnostics for fluxes, deposition)  fluxes 'down' are <0
    !     FJTOP(1:W_) = diffuse flux out top-of-atmosphere (TAU=0 above top model lyr)
    !     FJBOT(1:W_) = diffuse flux onto surface (<0 by definition)
    !     FSBOT(1:W_) = direct/solar flux onto surface  (<0 by definition)
    !     FJFLX(1:L_,1:W_) = diffuse flux across top of model layer L
    !        this connects with FJBOT = FJFLX(0) & FJTOP = FJFLX(L_+1) (not dim!!)
    !     FLXD(1:L_+1,1:W_) = solar flux deposited in layer L (includes lyr above CTM)
    !        this should take into account sphericity, and is not just = mu0
    !     FLXD0(1:W_) = sum of solar flux deposited in atmos
    !        does NOT include flux on lower surface, does NOT mean absorbed!
    ! --------------------------------------------------------------------
    !
    !     DTAU     Local optical depth of each CTM level
    !     TTAU     Optical depth of air vertically above each point (to top of atm)
    !     FTAU2     Attenuation of solar beam
    !     POMEGAJ  Scattering phase function
    !
    ! new ver 5.3 code adds sub-layers (# = JXTRA(L2)) using ATAU as the
    !   factor increase from sub-layer to sub-layer
    !
    ! ------------------SET UP FOR MIE CODE-------------------------------
    !
    ! --------------wavelength independent--------------------------------
    !
    !  Transpose the ascending TTAU grid to a descending ZTAU grid.
    !  Double the resolution - TTAU points become the odd points on the
    !  ZTAU grid, even points needed for asymm phase fn soln, contain 'h'.
    !  Odd point added at top of grid for unattenuated beam   (Z='inf')
    !
    !  The following mapping holds for JADDLV=0
    !        Surface:   TTAU(1)    ==> ZTAU(2*L2_+1)
    !        Top:       TTAU(L2_)  ==> ZTAU(3)
    !        Infinity:     0.0     ==> ZTAU(1)
    !        index: 2*(L2_+1-L2)+1 ==> LZ
    !
    !  Mie scattering code only used from surface to level L2_
    ! ---------------------------------------------------------------------
    !
    ! ---------------------------------------------------------------------
    !  Insert new levels, working downwards from the top of the atmosphere
    !  to the surface (down in 'LZ', up in 'L2'). This allows ztau and pomega
    !  to be incremented linearly, and the flux fz to be attenuated top-down
    !    (avoiding problems where lower level fluxes are zero).
    ! ---------------------------------------------------------------------
    !
    !  Ascend through atmosphere transposing grid and adding extra points
    !  remember L2=1 is surface of CTM, but last layer (LZ) in scattering code.
    !  there are twice the number of layers in the LZ arrays (2*L2_ + 2*JADDTO + 1)
    !    because we need to insert the intermediate layers (even LZ) for the
    !    asymmetric scattering code.
    !
    !  Transfer the L2=1:L2_+1 values (TTAU,FTAU2,POMEGAJ) onto the reverse
    !    order, expanded, doubled-level scatter grid.
    !    Note that we need to deal with the expansion by JADD levels (L2L).
    !      These JADDLV levels are skipped and need to be interpolated later.
    !    Note that only odd LZ levels are filled,
    !
    ! -------------------re-grid data---------------------------------------------
    !  Calculate cumulative total and define levels we want J-values at.
    !  Sum upwards for levels, and then downwards for Mie code readjustments.
    !
    !     JXTRA(L2)  Number of new levels to add between (L2) and (L2+1)
    !           ***JXTRA(1:L2_+1) is calculated based on the aerosol+cloud OD_s
    !     JADDLV(L2)  Number of new levels actually added at each wavelength
    !            where JADDLV = 0 when there is effectively no FTAU2
    !     JADDTO(L2)   Total number of new levels to add to and above level (L2)
    !     JNDLEV(L) = L2 index that maps on CTM mid-layer L
    !
    ! JADDLV(L2=1:L2_) = number of levels to add between TTAU2(L2) and TTAU(L2+1)
    !     JADDLV is taken from JXTRA, which is based on visible OD.
    !     JADDTO(L2=1:L2_+1) is the cumulative number of levels to be added
    ! these should be fixed for all wavelengths to lock-in the array sizes

    if (LU .gt. JXL_) then
       call EXITC (' OPMIE:  JXL_ .lt. L_')
    endif

    L1U = LU + 1
    L2U = 2*LU + 2

    do L2 = 1,L2U,1
       JADDLV(L2) = JXTRA(L2)
    enddo
    JADDTO(L2U+1) = 0
    do L2 = L2U,1,-1
       JADDTO(L2) = JADDTO(L2+1) + JADDLV(L2)
    enddo

    ! expanded grid now included CTM edge and mid layers plus expanded
    !     grid to allow for finer delta-tau at tops of clouds.
    !     DIM of new grid = L2U + JADDTO(1) + 1

    ! L2LEV(L2) = L2-index for old level L2 in expanded J-grid (w/JADDLV)
    !     in absence of JADDLV, L2LEV(L2) = L2
    L2LEV(1)  = 1
    do L2 = 2,L2U+1
       L2LEV(L2) = L2LEV(L2-1) + 1 + JADDLV(L2-1)
    enddo

    ! JNDLEV(L=1:L_) = L2-index in expanded grid for CTM mid-layer L
    ! JNELEV(L=1:L_) = L2-index for top of layer L
    do L = 1,LU
       JNDLEV(L) = L2LEV(2*L)
       JNELEV(L) = L2LEV(2*L+1)
    enddo
    JNELEV(LU+1) = 0  !need to set this to top-of-atmosphere

    ND = 2*L2U + 2*JADDTO(1) + 1

    if(ND .gt. N_) then
       call EXITC (' overflow of scatter arrays: ND > N_')
    endif

    ! -------------begin wavelength dependent set up---------------------------

    ! Reinitialize arrays
    ZTAU(:,:)     = 0.e+0_fp
    FZ(:,:)       = 0.e+0_fp
    POMEGA(:,:,:) = 0.e+0_fp

    do K=1,W_

       ! Set up optical depth DTAU(L)
       do L = 1,L1U
          DTAU(L,K) = DTAUX(L,K)
       enddo
       DTAU(L1U+1,K) = 0.e+0_fp

       ! Define the total scattering phase fn for each CTM layer L=1:L_+1
       !    from a DTAU-wt_d mix of aerosols, cloud & Rayleigh
       ! No. of quadrature pts fixed at 4(M_), expansion of phase fn @ 8
       do L = 1,L1U
          do I = 1,M2_
             POMEGAJ(I,L,K) = POMEGAX(I,L,K)
          enddo
       enddo

       ! Calculate attenuated incident beam exp(-TTAU/U0 = DTAU * AirMassFactor)
       !       at the middle & edges of the CTM layers L=1:2*L1_+1
       !   L1_ is top-edge of CTM (ie, L=38 = 2 hPa) which has TAU > 0
       !   note that DTAU(L1_) is optical depth in the FULL CTM layer just above
       FTAU2(:,:) = 0.e+0_fp
       FTAU2(L2U+1,:) = 1.0e+0_fp
       do LL = 1,2*L1U+1
          L = (LL+1)/2
          if (AMF2(LL,LL) .gt. 0.0e+0_fp) then
             XLTAU = 0.0e+0_fp
             do II = 1,2*L1U+1
                I = (II+1)/2
                XLTAU = XLTAU + 0.5e+0_fp*DTAU(I,K)*AMF2(II,LL)
             enddo
             if (XLTAU .lt. 76.e+0_fp) then   ! zero out flux at 1e-33
                FTAU2(LL,K) = exp(-XLTAU)
             endif
          endif
       enddo

       ! calculate direct solar flux deposited in each CTM half-layer: L=1:L2_
       !      use FSBOT for surface flux, cannot do layer above CTM (L_+1)
       FLXD2(:,:) = 0.e+0_fp
       do LL = 1,2*L1U
          if (AMF2(LL,LL) .gt. 0.e+0_fp) then
             FLXD2(LL,K) = (FTAU2(LL+1,K) - FTAU2(LL,K))/AMF2(LL,LL)
          endif
       enddo
       if (AMF2(1,1) .gt. 0.e+0_fp) then
          FSBOT(K) = FTAU2(1,K)/AMF2(1,1)
       else
          FSBOT(K) = 0.e+0_fp
       endif

       do LL = 2,2*L1U,2
          L=LL/2
          FLXD(L,K) = FLXD2(LL,K)+FLXD2(LL-1,K)
       enddo

       ! integrate solar flux depositied in CTM layers L=1:L_, cannot do top layer
       !   note FLXD0 .ne. (1.e+0_fp - FTAU(L_+1))/AMF(L_+1,L_+1) with spherical atmos
       FLXD0(K) = 0.e+0_fp
       if (AMF2(2*L1U,2*L1U) .gt. 0.e+0_fp) then
          do L=1,L1U
             FLXD0(K) = FLXD0(K) + FLXD(L,K)
          enddo
       endif

       ! ---------------------------------------------------------------------
       !  Take optical properties on CTM layers and convert to a photolysis
       !  level grid corresponding to layer centres and boundaries. This is
       !  required so that J-values can be calculated for the centre of CTM
       !  layers; the index of these layers is kept in the JNDLEV array.
       ! ---------------------------------------------------------------------
       ! Now combine the CTM layer edges (1:L_+2) with the CTM mid-layer
       !     points (1:L_) plus 1 for the mid point of added top layer.
       ! combine these edge- and mid-layer points into grid of size:
       !               L2_+1 = 2*L1_+1 = 2*L_+3
       ! calculate column optical depths above each level, TTAU(1:L2_+1)
       !       note that TTAU(L2_+1)=0 and TTAU(1)=total OD

       TTAU(L2U+1,K) = 0.0e+0_fp
       do L2 = L2U,1,-1
          L          = (L2+1)/2
          DTAUJ      = 0.5e+0_fp * DTAU(L,K)
          TTAU(L2,K)   = TTAU(L2+1,K) + DTAUJ
       enddo

       ! -solar flux incident on lower boundary & Lambertian reflect factor:
       if (FSBOT(K) .gt. 0.e+0_fp) then
          ZFLUX(K) = FSBOT(K)*RFL(K)/(1.e+0_fp+RFL(K))
       else
          ZFLUX(K) = 0.e+0_fp
       endif

       !  Calculate scattering properties, level centres then level boundaries
       !>>>>>be careful of order, we are overwriting/shifting the 'POMEGAJ'
       !     upward in index
       do L2 = L2U,2,-2
          L   = L2/2
          do I = 1,M2_
             POMEGAJ(I,L2,K) = POMEGAJ(I,L,K)
          enddo
       enddo
       ! lower boundary value is set (POMEGAJ(I,1)), but set upper:
       do I = 1,M2_
          POMEGAJ(I,L2U+1,K) = POMEGAJ(I,L2U,K)
       enddo
       ! now have POMEGAJ filled at even points from L2=3:L2_-1
       ! use inverse interpolation for correct tau-weighted values at edges
       do L2 = 3,L2U-1,2
          TAUDN = TTAU(L2-1,K)-TTAU(L2,K)
          TAUUP = TTAU(L2,K)-TTAU(L2+1,K)
          do I = 1,M2_
             POMEGAJ(I,L2,K) = (POMEGAJ(I,L2-1,K)*TAUDN + &
                                POMEGAJ(I,L2+1,K)*TAUUP) / (TAUDN+TAUUP)
          enddo
       enddo

       ! at this point FTAU2(1:L2_+1) and POMEAGJ(1:8, 1:L2_+1)
       !     where FTAU2(L2_+1) = 1.0 = top-of-atmos, FTAU2(1) = surface

       do L2 = 1,L2U+1          ! L2 = index of CTM edge- and mid-layers
          L2L = L2LEV(L2)        ! L2L = index for L2 in expanded scale(JADD)
          LZ  = ND + 2 - 2*L2L  ! LZ = index for L2 in scatt arrays
          ZTAU(LZ,K) = TTAU(L2,K)
          FZ(LZ,K)   = FTAU2(L2,K)
          do I=1,M2_
             POMEGA(I,LZ,K) = POMEGAJ(I,L2,K)
          enddo
       enddo

       ! Now go thru the pairs of L2 levels to see if we need JADD levels
       do L2 = 1,L2U             ! L2 = index of CTM edge- and mid-layers
          L2L = L2LEV(L2)        ! L2L = index for L2 in expanded scale(JADD)
          LZ  = ND + 2 - 2*L2L   ! LZ = index for L2 in scatt arrays
          L22 = L2LEV(L2+1) - L2LEV(L2) - 1   ! L22 = 0 if no added levels

          if (L22 .gt. 0) then
             TAUBTM(K) = TTAU(L2,K)
             TAUTOP(K) = TTAU(L2+1,K)
             FBTM(K)   = FTAU2(L2,K)
             FTOP(K)   = FTAU2(L2+1,K)
             do I = 1,M2_
                POMEGAB(I,K) = POMEGAJ(I,L2,K)
             enddo

             ! to fit L22 new layers between TAUBOT > TAUTOP, calculate new
             ! 1/ATAU factor
             ! such that TAU(just above TAU-btm) = ATUAZ * TAUBTM < TAUBTM
             ATAUZ = exp(-log(TAUBTM(K)/max(TAUTOP(K),ATAU0))/float(L22+1))
             do L = 1,L22      ! add odd levels between L2LEV(L2) & L2LEV(L2+1)
                LZZ = LZ - 2*L ! LZZ = index(odd) of added level in scatt arrays
                ZTAU(LZZ,K) = TAUBTM(K) * ATAUZ

                ! fraction from TAUBTM=>TAUTOP
                ATAUA=(TAUBTM(K)-ZTAU(LZZ,K))/(TAUBTM(K)-TAUTOP(K))
                ! solar flux at interp-levels: use exp(TAU/U0) if U0>0.02
                ! (89 deg), else scale by TAU
                if (U0 .gt. 0.02e+0_fp) then
                   FZ(LZZ,K) = FTOP(K) * exp((TAUTOP(K)-ZTAU(LZZ,K))/U0)
                else
                   if (FBTM(K) .lt. 1.d-32) then
                      FZ(LZZ,K) = 0.e+0_fp
                   else
                      FZ(LZZ,K) = FBTM(K) * (FTOP(K)/FBTM(K))**ATAUA
                   endif
                endif
                do I = 1,M2_
                   POMEGA(I,LZZ,K) = POMEGAB(I,K) + &
                        ATAUA*(POMEGAJ(I,L2+1,K)-POMEGAB(I,K))
                enddo
                TAUBTM(K)    = ZTAU(LZZ,K)
                FBTM(K)      = FZ(LZZ,K)
                do I = 1,M2_
                   POMEGAB(I,K) = POMEGA(I,LZZ,K)
                enddo
             enddo
          endif
       enddo

       ! Now fill in the even points with simple interpolation in scatter arrays
       do LZ = 2,ND-1,2
          ZTAU(LZ,K) = 0.5e+0_fp*(ZTAU(LZ-1,K)+ZTAU(LZ+1,K))
          FZ(LZ,K)   = sqrt(FZ(LZ-1,K)*FZ(LZ+1,K))
          do I=1,M2_
             POMEGA(I,LZ,K) = 0.5e+0_fp*(POMEGA(I,LZ-1,K)+POMEGA(I,LZ+1,K))
          enddo
       enddo

    enddo  ! wavelength loop!

    ! --------------------------------------------------------------------
    call MIESCT(FJ,FJT,FJB,POMEGA,FZ,ZTAU,ZFLUX,RFL,U0,ND)
    ! --------------------------------------------------------------------

    ! Move mean intensity from scatter array FJ(LZ=1:ND)
    !               to CTM mid-level array FJACT(L=1:L_)

    do K=1,W_

       ! mean intensity at mid-layer:  4*<I> + solar
       !do L = 1,LU
       ! L2L = JNDLEV(L)
       ! LZ  = ND+2 - 2*L2L
       ! FJACT(L,K) = 4.e+0_fp*FJ(LZ,K) + FZ(LZ,K)
       !enddo

       ! mean intensity averaged throughout layer:
       do L = 1,LU
          LZ0 = ND+2 - 2*JNELEV(L)
          if (L .gt. 1) then
             LZ1 = ND+2 - 2*JNELEV(L-1)
          else
             LZ1 = ND
          endif
          SUMJ = (4.e+0_fp*FJ(LZ0,K)+FZ(LZ0,K))*(ZTAU(LZ0+2,K)-ZTAU(LZ0,K)) &
               + (4.e+0_fp*FJ(LZ1,K)+FZ(LZ1,K))*(ZTAU(LZ1,K)-ZTAU(LZ1-2,K))
          SUMT = ZTAU(LZ0+2,K)-ZTAU(LZ0,K) + ZTAU(LZ1,K)-ZTAU(LZ1-2,K)

          do LZ = LZ0+2,LZ1-2,2
             SUMJ =SUMJ+(4.e+0_fp*FJ(LZ,K)+FZ(LZ,K))*(ZTAU(LZ+2,K)-ZTAU(LZ-2,K))
             SUMT =SUMT + ZTAU(LZ+2,K)-ZTAU(LZ-2,K)
          enddo
          FJACT(L,K) = SUMJ/SUMT

       enddo

       ! mean diffuse flux:  4<I*mu> (not solar) at top of layer L
       !       average (tau-wtd) the h's just above and below the L-edge
       do L = 1,LU
          L2L = JNELEV(L)
          LZ  = ND+2 - 2*L2L
          FJFLX0 = (ZTAU(LZ+1,K)-ZTAU(LZ,K))/(ZTAU(LZ+1,K)-ZTAU(LZ-1,K))
          FJFLX(L,K)=4.e+0_fp*(FJ(LZ-1,K)*FJFLX0 +FJ(LZ+1,K)*(1.e+0_fp-FJFLX0))
       enddo

       ! diffuse fluxes reflected at top, incident at bottom
       FJTOP(K) = FJT(K)
       FJBOT(K) = FJB(K)

    enddo  ! wavelength loop!

  END SUBROUTINE OPMIE
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: miesct
!
! !DESCRIPTION: Subroutine MIESCT is an adaptation of the Prather radiative
!  transfer code (mjp, 10/95).
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE MIESCT(FJ,FJT,FJB, POMEGA,FZ,ZTAU,ZFLUX,RFL,U0,ND)
!
! !USES:
!
    USE CMN_FJX_Mod, ONLY : EMU, M_, M2_, N_, W_
!
!
! !INPUT PARAMETERS:
!
    INTEGER,  INTENT(IN)  ::  ND
    REAL(fp), INTENT(IN)  ::  POMEGA(M2_,N_,W_),FZ(N_,W_),ZTAU(N_,W_), &
                              RFL(W_),U0,ZFLUX(W_)
!
! !OUTPUT VARIABLES:
!
    REAL(fp), INTENT(OUT) ::  FJ(N_,W_),FJT(W_),FJB(W_)
!
! !REMARKS:
!     Prather, 1974, Astrophys. J. 192, 787-792.
!         Solution of inhomogeneous Rayleigh scattering atmosphere.
!         (original Rayleigh w/ polarization)
!     Cochran and Trafton, 1978, Ap.J., 219, 756-762.
!         Raman scattering in the atmospheres of the major planets.
!         (first use of anisotropic code)
!     Jacob, Gottlieb and Prather, 1989, J.Geophys.Res., 94, 12975-13002.
!         Chemistry of a polluted cloudy boundary layer,
!         (documentation of extension to anisotropic scattering)
!                                                                             .
!    takes atmospheric structure and source terms from std J-code
!    ALSO limited to 4 Gauss points, only calculates mean field! (M=1)
!
! !REVISION HISTORY:
!  28 Mar 2013 - S. D. Eastham - Copied from Fast-JX v7.0
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    REAL(fp)  PM(M_,M2_),PM0(M2_)
    INTEGER I, IM  ,K

    !=================================================================
    ! MIESCT begins here!
    !=================================================================

    do I = 1,M_
       call LEGND0 (EMU(I),PM0,M2_)
       do IM = 1,M2_
          PM(I,IM) = PM0(IM)
       enddo
    enddo

    call LEGND0 (-U0,PM0,M2_)
    do IM=1,M2_
       PM0(IM) = 0.25e+0_fp*PM0(IM)
    enddo

    ! BLKSLV now called with all the wavelength arrays (K=1:W_)
    call BLKSLV(FJ,POMEGA,FZ,ZTAU,ZFLUX,RFL,PM,PM0,FJT,FJB, ND)

  END SUBROUTINE MIESCT
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: legnd0
!
! !DESCRIPTION: Subroutine LEGND0 calculates ordinary Legendre functions
!  of X (real) from $P[0] = PL(1) = 1, P[1] = X, \dots, P[N-1] = PL(N)$
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE LEGND0 (X,PL,N)
!
! !USES:
!
!
! !INPUT PARAMETERS:
!
    INTEGER,  INTENT(IN)  :: N
    REAL(fp), INTENT(IN)  :: X
!
! !OUTPUT VARIABLES:
!
    REAL(fp), INTENT(OUT) :: PL(N)
!
! !REMARKS:
!
! !REVISION HISTORY:
!  28 Mar 2013 - S. D. Eastham - Copied from Fast-JX v7.0
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER I
    REAL(fp)  DEN

    !=================================================================
    ! LEGND0 begins here!
    !=================================================================

    ! Always does PL(2) = P[1]
    PL(1) = 1.e+0_fp
    PL(2) = X
    do I = 3,N
       DEN = (I-1)
       PL(I) = PL(I-1)*X*(2.e+0_fp-1.0/DEN) - PL(I-2)*(1.e+0_fp-1.e+0_fp/DEN)
    enddo

  END SUBROUTINE LEGND0
!EOC

END MODULE FJX_MOD
