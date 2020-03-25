!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: hemcox_dustginoux_mod.F90
!
! !DESCRIPTION: Paul GINOUX dust source function.  This subroutine updates
!  the surface mixing ratio of dust aerosols for NDSTBIN size bins.  The
!  uplifting of dust depends in space on the source function, and in time
!  and space on the soil moisture and surface wind speed (10 meters).  Dust
!  is uplifted if the wind speed is greater than a threshold velocity which
!  is calculated with the formula of Marticorena et al.  (JGR, v.102,
!  pp 23277-23287, 1997).  To run this subroutine you need the source
!  function which can be obtained by contacting Paul Ginoux at
!  ginoux@rondo.gsfc.nasa.gov/  If you are not using GEOS DAS met fields,
!  you will most likely need to adapt the adjusting parameter.
!\\
!\\
! This is a HEMCO extension module that uses many of the HEMCO core
! utilities.
!\\
!\\
! References:
!
! \begin{enumerate}
! \item Ginoux, P., M. Chin, I. Tegen, J. Prospero, B. Hoben, O. Dubovik,
!        and S.-J. Lin, "Sources and distributions of dust aerosols simulated
!        with the GOCART model", J. Geophys. Res., 2001
! \item Chin, M., P. Ginoux, S. Kinne, B. Holben, B. Duncan, R. Martin,
!        J. Logan, A. Higurashi, and T. Nakajima, "Tropospheric aerosol
!        optical thickness from the GOCART model and comparisons with
!        satellite and sunphotometers measurements", J. Atmos Sci., 2001.
! \end{enumerate}
!
! !AUTHOR:
!  Paul Ginoux (ginoux@rondo.gsfc.nasa.gov)
!
! !INTERFACE:
!
MODULE HCOX_DustGinoux_Mod
!
! !USES:
!
  USE HCO_Error_Mod
  USE HCO_Diagn_Mod
  USE HCO_State_Mod,  ONLY : HCO_State
  USE HCOX_State_Mod, ONLY : Ext_State

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC :: HcoX_DustGinoux_Run
  PUBLIC :: HcoX_DustGinoux_Init
  PUBLIC :: HcoX_DustGinoux_Final
  PUBLIC :: HcoX_DustGinoux_GetChDust
!
! !REVISION HISTORY:
!  08 Apr 2004 - T. D. Fairlie - Initial version
!  (1 ) Added OpenMP parallelization (bmy, 4/8/04)
!  (2 ) Now references DATA_DIR from "directory_mod.f" (bmy, 7/20/04)
!  25 Aug 2010 - R. Yantosca - Added ProTeX headers
!  01 Mar 2012 - R. Yantosca - Now use GET_AREA_M2(I,J,L) from grid_mod.F90
!  01 Aug 2012 - R. Yantosca - Add reference to findFreeLUN from inqure_mod.F90
!  03 Aug 2012 - R. Yantosca - Move calls to findFreeLUN out of DEVEL block
!  09 Nov 2012 - M. Payer    - Replaced all met field arrays with State_Met
!                              derived type object
!  26 Feb 2013 - R. Yantosca - Now accept Input_Opt via the arg list
!  11 Dec 2013 - C. Keller   - Now a HEMCO extension.
!  29 Sep 2014 - R. Yantosca - Now make NBINS a variable and not a parameter
!  29 Sep 2014 - R. Yantosca - Now use F90 free-format indentation
!  08 Jul 2015 - M. Sulprizio- Now include dust alkalinity source (tdf 04/10/08)
!  14 Oct 2016 - C. Keller   - Now use HCO_EvalFld instead of HCO_GetPtr.
!  25 Jan 2019 - M. Sulprizio- Add instance wrapper
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !PRIVATE TYPES:
!
  TYPE :: MyInst

   ! Quantities related to dust bins
   INTEGER              :: Instance
   INTEGER              :: NBINS
   INTEGER              :: ExtNr    = -1     ! Extension number  for DustGinoux
   INTEGER              :: ExtNrAlk = -1     ! Extension number  for DustAlk
   INTEGER, ALLOCATABLE :: HcoIDs    (:)     ! HEMCO species IDs for DustGinoux
   INTEGER, ALLOCATABLE :: HcoIDsAlk (:)     ! HEMCO species IDs for DustAlk
   INTEGER,  POINTER    :: IPOINT    (:)     ! 1=sand, 2=silt, 3=clay
   REAL,     POINTER    :: FRAC_S    (:)     !
   REAL,     POINTER    :: DUSTDEN   (:)     ! dust density     [kg/m3]
   REAL,     POINTER    :: DUSTREFF  (:)     ! effective radius [um]
   REAL(hp), POINTER    :: FLUX(:,:,:)
   REAL(hp), POINTER    :: FLUX_ALK(:,:,:)

   ! Source functions (get from HEMCO core)
   REAL(hp), POINTER    :: SRCE_SAND(:,:) => NULL()
   REAL(hp), POINTER    :: SRCE_SILT(:,:) => NULL()
   REAL(hp), POINTER    :: SRCE_CLAY(:,:) => NULL()

   ! Transfer coefficient (grid-dependent)
   REAL(dp)             :: CH_DUST

   TYPE(MyInst), POINTER :: NextInst => NULL()
  END TYPE MyInst

  ! Pointer to instances
  TYPE(MyInst), POINTER  :: AllInst => NULL()

CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCOX_DustGinoux_Run
!
! !DESCRIPTION: Subroutine HcoX\_DustGinoux\_Run is the driver routine
! for the Paul Ginoux dust source function HEMCO extension.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HcoX_DustGinoux_Run( ExtState, HcoState, RC )
!
! !USES:
!
    USE HCO_Calc_Mod,     ONLY : HCO_EvalFld
    USE HCO_EmisList_Mod, ONLY : HCO_GetPtr
    USE HCO_FluxArr_Mod,  ONLY : HCO_EmisAdd
    USE HCO_Clock_Mod,    ONLY : HcoClock_First
!
! !INPUT PARAMETERS:
!
    TYPE(Ext_State), POINTER        :: ExtState    ! Options for this ext
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(HCO_State), POINTER        :: HcoState    ! HEMCO state object
    INTEGER,         INTENT(INOUT)  :: RC          ! Success or failure?
!
! !REMARKS:
!    SRCE_FUNK Source function                               (-)
!              for 1: Sand, 2: Silt, 3: Clay
!                                                                             .
!    DUSTDEN   Dust density                                  (kg/m3)
!    DUSTREFF  Effective radius                              (um)
!    AD        Air mass for each grid box                    (kg)
!    NTDT      Time step                                     (s)
!    W10m      Velocity at the anemometer level (10meters)   (m/s)
!    GWET      Surface wetness                               (-)
!                                                                             .
!  Dust properties used in GOCART
!                                                                             .
!  Size classes: 01-1, 1-1.8, 1.8-3, 3-6 (um)
!  Radius: 0.7, 1.5, 2.5, 4  (um)
!  Density: 2500, 2650, 2650, 2650 (kg/m3)
!
! !REVISION HISTORY:
!  08 Apr 2004 - T. D. Fairlie - Initial version
!  (1 ) Added OpenMP parallelization (bmy, 4/8/04)
!  (2 ) Now references DATA_DIR from "directory_mod.f" (bmy, 7/20/04)
!  25 Aug 2010 - R. Yantosca - Added ProTeX headers
!  01 Mar 2012 - R. Yantosca - Now use GET_AREA_M2(I,J,L) from grid_mod.F90
!  01 Aug 2012 - R. Yantosca - Add reference to findFreeLUN from inqure_mod.F90
!  03 Aug 2012 - R. Yantosca - Move calls to findFreeLUN out of DEVEL block
!  09 Nov 2012 - M. Payer    - Replaced all met field arrays with State_Met
!                              derived type object
!  26 Feb 2013 - R. Yantosca - Now accept Input_Opt via the arg list
!  11 Dec 2013 - C. Keller   - Now a HEMCO extension
!  29 Sep 2014 - R. Yantosca - Bug fix: SRCE_CLAY should have been picked when
!                              M=3 but was picked when M=2.  Now corrected.
!  26 Jun 2015 - E. Lundgren - Add L. Zhang new dust size distribution scheme
!  08 Jul 2015 - M. Sulprizio- Now include dust alkalinity source (tdf 04/10/08)
!  26 Oct 2016 - R. Yantosca - Don't nullify local ptrs in declaration stmts
!  07 Jul 2017 - R. Yantosca - Bug fix: Skip DustAlk IF block unless that
!                              extension has been turned on in the config file
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETER:
!
    REAL*8, PARAMETER :: RHOA     = 1.25d-3

!
! !LOCAL VARIABLES:
!
    ! SAVED scalars
!    LOGICAL, SAVE     :: FIRST = .TRUE.

    ! Scalars
    INTEGER           :: I, J, N, M, tmpID
    LOGICAL           :: ERR
    REAL*8            :: W10M,   DEN,    DIAM,   U_TS0, U_TS
    REAL*8            :: SRCE_P, REYNOL, ALPHA,  BETA
    REAL*8            :: GAMMA,  CW,     DTSRCE, A_M2,  G
    REAL              :: DSRC
    CHARACTER(LEN=63) :: MSG

    ! Arrays
    REAL*8            :: DUST_EMI_TOTAL(HcoState%NX, HcoState%NY)

    ! Pointers
    TYPE(MyInst), POINTER :: Inst
    REAL(hp),     POINTER :: Arr2D(:,:)

    !=======================================================================
    ! HCOX_DUSTGINOUX_RUN begins here!
    !=======================================================================

    ! Return if extension is disabled
    IF ( ExtState%DustGinoux <= 0 ) RETURN

    ! Enter
    CALL HCO_ENTER(HcoState%Config%Err,'HCOX_DustGinoux_Run (hcox_dustginoux_mod.F90)',RC)
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Get instance
    Inst   => NULL()
    CALL InstGet ( ExtState%DustGinoux, Inst, RC )
    IF ( RC /= HCO_SUCCESS ) THEN
       WRITE(MSG,*) 'Cannot find DustGinoux instance Nr. ', ExtState%DustGinoux
       CALL HCO_ERROR(HcoState%Config%Err,MSG,RC)
       RETURN
    ENDIF

    ! Set gravity at earth surface (cm/s^2)
    G       = HcoState%Phys%g0 * 1.0d2

    ! Emission timestep [s]
    DTSRCE  = HcoState%TS_EMIS

    ! Initialize total dust emissions array [kg/m2/s]
    DUST_EMI_TOTAL = 0.0d0

    ! Error check
    ERR     = .FALSE.

    ! Init
    Arr2D    => NULL()

    !=================================================================
    ! Point to DUST source functions
    !=================================================================
    !IF ( HcoClock_First(HcoState%Clock,.TRUE.) ) THEN

       ! Sand
       CALL HCO_EvalFld( HcoState, 'GINOUX_SAND', Inst%SRCE_SAND, RC )
       IF ( RC /= HCO_SUCCESS ) RETURN

       ! Silt
       CALL HCO_EvalFld( HcoState, 'GINOUX_SILT', Inst%SRCE_SILT, RC )
       IF ( RC /= HCO_SUCCESS ) RETURN

       ! Clay
       CALL HCO_EvalFld( HcoState, 'GINOUX_CLAY', Inst%SRCE_CLAY, RC )
       IF ( RC /= HCO_SUCCESS ) RETURN
    !ENDIF

    !=================================================================
    ! Compute dust emisisons
    !=================================================================
!$OMP PARALLEL DO                                             &
!$OMP DEFAULT( SHARED )                                       &
!$OMP PRIVATE( I,      J,     M,      N,      DEN,   DIAM   ) &
!$OMP PRIVATE( REYNOL, ALPHA, BETA,   GAMMA,  U_TS0, U_TS   ) &
!$OMP PRIVATE( CW,     W10M,  SRCE_P, RC                    ) &
!$OMP SCHEDULE( DYNAMIC )
    DO N = 1, Inst%NBINS

       !====================================================================
       ! Threshold velocity as a function of the dust density and the
       ! diameter from Bagnold (1941), valid for particles larger
       ! than 10 um.
       !
       ! u_ts0 = 6.5*sqrt(dustden(n)*g0*2.*dustreff(n))
       !
       ! Threshold velocity from Marticorena and Bergametti
       ! Convert units to fit dimensional parameters
       !====================================================================
       DEN    = Inst%DUSTDEN(N) * 1.d-3                   ! [g/cm3]
       DIAM   = 2d0 * Inst%DUSTREFF(N) * 1.d2             ! [cm in diameter]
       REYNOL = 1331.d0 * DIAM**(1.56d0) + 0.38d0         ! [Reynolds number]
       ALPHA  = DEN * G * DIAM / RHOA
       BETA   = 1d0 + ( 6.d-3 / ( DEN * G * DIAM**(2.5d0) ) )
       GAMMA  = ( 1.928d0 * REYNOL**(0.092d0) ) - 1.d0

       !====================================================================
       ! I think the 129.d-5 is to put U_TS in m/sec instead of cm/sec
       ! This is a threshold friction velocity!       from M&B
       ! i.e. Ginoux uses the Gillette and Passi formulation
       ! but has substituted Bagnold's Ut with M&B's U*t.
       ! This appears to be a problem.  (tdf, 4/2/04)
       !====================================================================

       ! [m/s]
       U_TS0  = 129.d-5 * SQRT( ALPHA ) * SQRT( BETA ) / SQRT( GAMMA )

       ! Index used to select the source function (1=sand, 2=silt, 3=clay)
       M = Inst%IPOINT(N)

       ! Loop over grid boxes
       DO J = 1, HcoState%NY
       DO I = 1, HcoState%NX

          ! Fraction of emerged surfaces
          ! (subtract lakes, coastal ocean,...)
          CW = 1.d0

          ! Case of surface dry enough to erode
          IF ( ExtState%GWETTOP%Arr%Val(I,J) < 0.2d0 ) THEN

             U_TS = U_TS0 *( 1.2d0 + 0.2d0 * &
                    LOG10( MAX(1.d-3,ExtState%GWETTOP%Arr%Val(I,J))))
             U_TS = MAX( 0.d0, U_TS )

          ELSE

             ! Case of wet surface, no erosion
             U_TS = 100.d0

          ENDIF

          ! 10m wind speed squared [m2/s2]
          W10M = ExtState%U10M%Arr%Val(I,J)**2 &
               + ExtState%V10M%Arr%Val(I,J)**2

          ! Get source function
          SELECT CASE( M )
             CASE( 1 )
                SRCE_P = Inst%SRCE_SAND(I,J)
             CASE( 2 )
                SRCE_P = Inst%SRCE_SILT(I,J)
             CASE( 3 )
                SRCE_P = Inst%SRCE_CLAY(I,J)
          END SELECT

          ! Units are m2
          SRCE_P = Inst%FRAC_S(N) * SRCE_P !* A_M2

          ! Dust source increment [kg/m2/s]
          Inst%FLUX(I,J,N) = CW           * Inst%CH_DUST * SRCE_P * W10M &
                           * ( SQRT(W10M) - U_TS )

          ! Not less than zero
          IF ( Inst%FLUX(I,J,N) < 0.d0 ) Inst%FLUX(I,J,N) = 0.d0

          ! Increment total dust emissions [kg/m2/s] (L. Zhang, 6/26/15)
          DUST_EMI_TOTAL(I,J) = DUST_EMI_TOTAL(I,J) + Inst%FLUX(I,J,N)

          ! Include DUST Alkalinity SOURCE, assuming an alkalinity
          ! of 4% by weight [kg].                  !tdf 05/10/08
          !tdf 3% Ca + equ 1% Mg = 4% alkalinity
          IF ( Inst%ExtNrAlk > 0 ) THEN
             Inst%FLUX_ALK(I,J,N) = 0.04 * Inst%FLUX(I,J,N)
          ENDIF

       ENDDO
       ENDDO
    ENDDO
!$OMP END PARALLEL DO

    ! Error check
    IF ( ERR ) THEN
       RC = HCO_FAIL
       RETURN
    ENDIF

    ! Redistribute dust emissions across bins (L. Zhang, 6/26/15)
!$OMP PARALLEL DO                                           &
!$OMP DEFAULT( SHARED )                                     &
!$OMP PRIVATE( I, J, N )                                    &
!$OMP SCHEDULE( DYNAMIC )
     DO N=1,Inst%NBINS
     DO J=1,HcoState%NY
     DO I=1,HcoState%NX
        SELECT CASE( N )
           CASE( 1 )
              Inst%FLUX(I,J,N) = DUST_EMI_TOTAL(I,J) * 0.0766d0
           CASE( 2 )
              Inst%FLUX(I,J,N) = DUST_EMI_TOTAL(I,J) * 0.1924d0
           CASE( 3 )
              Inst%FLUX(I,J,N) = DUST_EMI_TOTAL(I,J) * 0.3491d0
           CASE( 4 )
              Inst%FLUX(I,J,N) = DUST_EMI_TOTAL(I,J) * 0.3819d0
        END SELECT
     ENDDO
     ENDDO
     ENDDO
!$OMP END PARALLEL DO

    !=======================================================================
    ! PASS TO HEMCO STATE AND UPDATE DIAGNOSTICS
    !=======================================================================
    DO N = 1, Inst%NBINS
       IF ( Inst%HcoIDs(N) > 0 ) THEN

          ! Add flux to emission array
          CALL HCO_EmisAdd( HcoState,       Inst%FLUX(:,:,N), &
                            Inst%HcoIDs(N), RC,       ExtNr=Inst%ExtNr   )
          IF ( RC /= HCO_SUCCESS ) THEN
             WRITE(MSG,*) 'HCO_EmisAdd error: dust bin ', N
             CALL HCO_ERROR(HcoState%Config%Err,MSG, RC )
             RETURN
          ENDIF

       ENDIF

       ! This block is only relevant if the DustAlk extension
       ! has been turned on.  Skip othewrise. (bmy, 7/7/17)
       IF ( Inst%ExtNrAlk > 0 ) THEN
          IF ( Inst%HcoIDsAlk(N) > 0 ) THEN

             ! Add flux to emission array
             CALL HCO_EmisAdd( HcoState,          Inst%FLUX_Alk(:,:,N), &
                               Inst%HcoIDsAlk(N), RC, ExtNr=Inst%ExtNrAlk)
             IF ( RC /= HCO_SUCCESS ) THEN
                WRITE(MSG,*) 'HCO_EmisAdd error: dust alkalinity bin ', N
                CALL HCO_ERROR(HcoState%Config%Err,MSG, RC )
                RETURN
             ENDIF
          ENDIF
       ENDIF

    ENDDO

    !=======================================================================
    ! Cleanup & quit
    !=======================================================================

    ! Nullify pointers
    Inst    => NULL()

    ! Leave w/ success
    CALL HCO_LEAVE( HcoState%Config%Err,RC )

  END SUBROUTINE HcoX_DustGinoux_Run
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCOX_DustGinoux_Init
!
! !DESCRIPTION: Subroutine HcoX\_DustGinoux\_Init initializes the HEMCO
! DUSTGINOUX extension.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HcoX_DustGinoux_Init( HcoState, ExtName, ExtState, RC )
!
! !USES:
!
    USE HCO_ExtList_Mod, ONLY : GetExtNr, GetExtOpt
    USE HCO_State_Mod,   ONLY : HCO_GetExtHcoID
!
! !INPUT PARAMETERS:
!
    TYPE(HCO_State),  POINTER        :: HcoState   ! HEMCO State object
    CHARACTER(LEN=*), INTENT(IN   )  :: ExtName    ! Extension name
    TYPE(Ext_State),  POINTER        :: ExtState   ! Extension options
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(INOUT)  :: RC         ! Success or failure?
!
! !REVISION HISTORY:
!  11 Dec 2013 - C. Keller   - Now a HEMCO extension
!  26 Sep 2014 - R. Yantosca - Updated for TOMAS
!  29 Sep 2014 - R. Yantosca - Now initialize NBINS from HcoState%N_DUST_BINS
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER                        :: N, AS, nSpc, nSpcAlk, ExtNr
    CHARACTER(LEN=255)             :: MSG
    REAL(dp)                       :: Mp, Rp, TmpScal
    LOGICAL                        :: FOUND

    ! Arrays
    CHARACTER(LEN=31), ALLOCATABLE :: SpcNames(:)
    CHARACTER(LEN=31), ALLOCATABLE :: SpcNamesAlk(:)

    ! Pointers
    TYPE(MyInst), POINTER          :: Inst

    !=======================================================================
    ! HCOX_DUSTGINOUX_INIT begins here!
    !=======================================================================

    ! Extension Nr.
    ExtNr = GetExtNr( HcoState%Config%ExtList, TRIM(ExtName) )
    IF ( ExtNr <= 0 ) RETURN

    ! Create Instance
    Inst => NULL()
    CALL InstCreate ( ExtNr, ExtState%DustGinoux, Inst, RC )
    IF ( RC /= HCO_SUCCESS ) THEN
       CALL HCO_ERROR ( HcoState%Config%Err, 'Cannot create DustGinoux instance', RC )
       RETURN
    ENDIF
    ! Also fill Inst%ExtNr
    Inst%ExtNr = ExtNr

    ! Check for dust alkalinity option
    Inst%ExtNrAlk = GetExtNr( HcoState%Config%ExtList, 'DustAlk' )

    ! Enter
    CALL HCO_ENTER(HcoState%Config%Err,'HCOX_DustGinoux_Init (hcox_dustginoux_mod.F90)',RC)
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Get the expected number of dust species
    Inst%NBINS = HcoState%nDust

    ! Get the actual number of dust species defined for DustGinoux extension
    CALL HCO_GetExtHcoID( HcoState, Inst%ExtNr, Inst%HcoIDs, &
                          SpcNames, nSpc, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Get the dust alkalinity species defined for DustAlk option
    IF ( Inst%ExtNrAlk > 0 ) THEN
       CALL HCO_GetExtHcoID( HcoState,    Inst%ExtNrAlk, Inst%HcoIDsAlk, &
                             SpcNamesAlk, nSpcAlk,  RC)
       IF ( RC /= HCO_SUCCESS ) RETURN
    ENDIF

    ! Make sure the # of dust species is as expected
    IF ( nSpc /= Inst%NBINS ) THEN
       WRITE( MSG, 100 ) Inst%NBINS, nSpc
 100   FORMAT( 'Expected ', i3, ' DustGinoux species but only found ', i3, &
               ' in the HEMCO configuration file!  Exiting...' )
       CALL HCO_ERROR(HcoState%Config%Err,MSG, RC )
       RETURN
    ENDIF

    ! Set scale factor: first try to read from configuration file. If
    ! not specified, call wrapper function which sets teh scale factor
    ! based upon compiler switches.
    CALL GetExtOpt( HcoState%Config, Inst%ExtNr, 'Mass tuning factor', &
                     OptValDp=TmpScal, Found=FOUND, RC=RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Set parameter FLX_MSS_FDG_FCT to specified tuning factor. Get from
    ! wrapper routine if not defined in configuration file
    IF ( FOUND ) THEN
       Inst%CH_DUST = TmpScal
    ELSE
       ! Get global mass flux tuning factor
       Inst%CH_DUST = HcoX_DustGinoux_GetCHDust( Inst, HcoState )
       IF ( Inst%CH_DUST < 0.0_dp ) THEN
          RC = HCO_FAIL
          RETURN
       ENDIF
    ENDIF

    ! Verbose mode
    IF ( HcoState%amIRoot ) THEN
       MSG = 'Use Ginoux dust emissions (extension module)'
       CALL HCO_MSG(HcoState%Config%Err,MSG )

       IF ( Inst%ExtNrAlk > 0 ) THEN
          MSG = 'Use dust alkalinity option'
          CALL HCO_MSG(HcoState%Config%Err,MSG, SEP1='-' )
       ENDIF

       MSG = 'Use the following species (Name: HcoID):'
       CALL HCO_MSG(HcoState%Config%Err,MSG)
       DO N = 1, nSpc
          WRITE(MSG,*) TRIM(SpcNames(N)), ':', Inst%HcoIDs(N)
          CALL HCO_MSG(HcoState%Config%Err,MSG)
       ENDDO
       IF ( Inst%ExtNrAlk > 0 ) THEN
          DO N = 1, nSpcAlk
             WRITE(MSG,*) TRIM(SpcNamesAlk(N)), ':', Inst%HcoIDsAlk(N)
             CALL HCO_MSG(HcoState%Config%Err,MSG)
          ENDDO
       ENDIF

       WRITE(MSG,*) 'Global mass flux tuning factor: ', Inst%CH_DUST
       CALL HCO_MSG(HcoState%Config%Err,MSG,SEP2='-')
    ENDIF

    ! Allocate vectors holding bin-specific informations
    ALLOCATE ( Inst%IPOINT  (Inst%NBINS) )
    ALLOCATE ( Inst%FRAC_S  (Inst%NBINS) )
    ALLOCATE ( Inst%DUSTDEN (Inst%NBINS) )
    ALLOCATE ( Inst%DUSTREFF(Inst%NBINS) )
    ALLOCATE ( Inst%FLUX    (HcoState%NX,HcoState%NY,Inst%NBINS) )
    ALLOCATE ( Inst%FLUX_ALK(HcoState%NX,HcoState%NY,Inst%NBINS) )

    ! Allocate arrays
    ALLOCATE ( Inst%SRCE_SAND ( HcoState%NX, HcoState%NY ), &
               Inst%SRCE_SILT ( HcoState%NX, HcoState%NY ), &
               Inst%SRCE_CLAY ( HcoState%NX, HcoState%NY ), &
               STAT = AS )
    IF ( AS /= 0 ) THEN
       CALL HCO_ERROR(HcoState%Config%Err,'Allocation error', RC )
       RETURN
    ENDIF

    ! Init
    Inst%FLUX      = 0.0_hp
    Inst%FLUX_ALK  = 0.0_hp
    Inst%SRCE_SAND = 0.0_hp
    Inst%SRCE_SILT = 0.0_hp
    Inst%SRCE_CLAY = 0.0_hp


    !=======================================================================
    ! Setup for simulations that use 4 dust bins (w/ or w/o TOMAS)
    !=======================================================================

    ! Fill bin-specific information
    IF ( Inst%NBINS == 4 ) THEN

       Inst%IPOINT  (1:Inst%NBINS) = (/ 3,       2,       2,       2       /)
       Inst%FRAC_S  (1:Inst%NBINS) = (/ 0.095d0, 0.3d0,   0.3d0,   0.3d0   /)
       Inst%DUSTDEN (1:Inst%NBINS) = (/ 2500.d0, 2650.d0, 2650.d0, 2650.d0 /)
       Inst%DUSTREFF(1:Inst%NBINS) = (/ 0.73d-6, 1.4d-6,  2.4d-6,  4.5d-6  /)

    ELSE

#if !defined( TOMAS )
       MSG = 'Cannot have > 4 GINOUX dust bins unless you are using TOMAS!'
       CALL HCO_ERROR(HcoState%Config%Err,MSG, RC )
       RETURN
#endif

    ENDIF

#if defined( TOMAS )

    !=======================================================================
    ! Setup for TOMAS simulations using more than 4 dust bins
    !
    ! from Ginoux:
    ! The U.S. Department of Agriculture (USDA) defines particles
    ! with a radius between 1 um and 25 um as silt, and below 1 um
    ! as clay [Hillel, 1982]. Mineralogical silt particles are mainly
    ! composed of quartz, but they are often coated with strongly
    ! adherent clay such that their physicochemical properties are
    ! similar to clay [Hillel, 1982].
    !
    ! SRCE_FUNC Source function
    ! for 1: Sand, 2: Silt, 3: Clay
    !=======================================================================
    IF ( Inst%NBINS == HcoState%MicroPhys%nBins ) THEN

       !--------------------------------------------------------------------
       ! Define the IPOINT array based on particle size
       !--------------------------------------------------------------------

       ! Loop over # of TOMAS bins
       DO N = 1, HcoState%MicroPhys%nBins

          ! Compute particle mass and radius
          Mp = 1.4 * HcoState%MicroPhys%BinBound(N)
          Rp = ( ( Mp /2500. ) * (3./(4.*HcoState%Phys%PI)))**(0.333)

          ! Pick the source function based on particle size
          IF ( Rp < 1.d-6 ) THEN
             Inst%IPOINT(N) = 3
          ELSE
             Inst%IPOINT(N) = 2
          END IF
       END DO

       !--------------------------------------------------------------------
       ! Set up dust density (DUSTDEN) array
       !--------------------------------------------------------------------
       DO N = 1, HcoState%MicroPhys%nBins
          IF ( HcoState%MicroPhys%BinBound(N) < 4.0D-15 ) THEN
             Inst%DUSTDEN(N)  = 2500.d0
          ELSE
             Inst%DUSTDEN(N)  = 2650.d0
          ENDIF
       ENDDO

       !--------------------------------------------------------------------
       ! Set up dust density (DUSTDEN) array
       !--------------------------------------------------------------------
       DO N = 1, HcoState%MicroPhys%nBins
          Inst%DUSTREFF(N) = 0.5d0                                    &
                      * ( SQRT( HcoState%MicroPhys%BinBound(N) *      &
                                HcoState%MicroPhys%BinBound(N+1) )    &
                      /   Inst%DUSTDEN(N) * 6.d0/HcoState%Phys%PI )**( 0.333d0 )
       ENDDO

       !--------------------------------------------------------------------
       ! Set up the FRAC_S array
       !--------------------------------------------------------------------

       ! Initialize
       Inst%FRAC_S( 1:HcoState%MicroPhys%nBins )           = 0d0

# if  defined( TOMAS12 ) || defined( TOMAS15 )

       !---------------------------------------------------
       ! TOMAS simulations with 12 or 15 size bins
       !---------------------------------------------------
       Inst%FRAC_S( HcoState%MicroPhys%nActiveModeBins + 1  )  = 7.33E-10
       Inst%FRAC_S( HcoState%MicroPhys%nActiveModeBins + 2  )  = 2.032E-08
       Inst%FRAC_S( HcoState%MicroPhys%nActiveModeBins + 3  )  = 3.849E-07
       Inst%FRAC_S( HcoState%MicroPhys%nActiveModeBins + 4  )  = 5.01E-06
       Inst%FRAC_S( HcoState%MicroPhys%nActiveModeBins + 5  )  = 4.45E-05
       Inst%FRAC_S( HcoState%MicroPhys%nActiveModeBins + 6  )  = 2.714E-04
       Inst%FRAC_S( HcoState%MicroPhys%nActiveModeBins + 7  )  = 1.133E-03
       Inst%FRAC_S( HcoState%MicroPhys%nActiveModeBins + 8  )  = 3.27E-03
       Inst%FRAC_S( HcoState%MicroPhys%nActiveModeBins + 9  )  = 6.81E-03
       Inst%FRAC_S( HcoState%MicroPhys%nActiveModeBins + 10 )  = 1.276E-02
       Inst%FRAC_S( HcoState%MicroPhys%nActiveModeBins + 11 )  = 2.155E-01
       Inst%FRAC_S( HcoState%MicroPhys%nActiveModeBins + 12 )  = 6.085E-01

# else

       !---------------------------------------------------
       ! TOMAS simulations with 30 or 40 size bins
       !---------------------------------------------------
       Inst%FRAC_S( HcoState%MicroPhys%nActiveModeBins +  1 )  = 1.05d-10
       Inst%FRAC_S( HcoState%MicroPhys%nActiveModeBins +  2 )  = 6.28d-10
       Inst%FRAC_S( HcoState%MicroPhys%nActiveModeBins +  3 )  = 3.42d-09
       Inst%FRAC_S( HcoState%MicroPhys%nActiveModeBins +  4 )  = 1.69d-08
       Inst%FRAC_S( HcoState%MicroPhys%nActiveModeBins +  5 )  = 7.59d-08
       Inst%FRAC_S( HcoState%MicroPhys%nActiveModeBins +  6 )  = 3.09d-07
       Inst%FRAC_S( HcoState%MicroPhys%nActiveModeBins +  7 )  = 1.15d-06
       Inst%FRAC_S( HcoState%MicroPhys%nActiveModeBins +  8 )  = 3.86d-06
       Inst%FRAC_S( HcoState%MicroPhys%nActiveModeBins +  9 )  = 1.18d-05
       Inst%FRAC_S( HcoState%MicroPhys%nActiveModeBins + 10 )  = 3.27d-05
       Inst%FRAC_S( HcoState%MicroPhys%nActiveModeBins + 11 )  = 8.24d-05
       Inst%FRAC_S( HcoState%MicroPhys%nActiveModeBins + 12 )  = 1.89d-04
       Inst%FRAC_S( HcoState%MicroPhys%nActiveModeBins + 13 )  = 3.92d-04
       Inst%FRAC_S( HcoState%MicroPhys%nActiveModeBins + 14 )  = 7.41d-04
       Inst%FRAC_S( HcoState%MicroPhys%nActiveModeBins + 15 )  = 1.27d-03
       Inst%FRAC_S( HcoState%MicroPhys%nActiveModeBins + 16 )  = 2.00d-03
       Inst%FRAC_S( HcoState%MicroPhys%nActiveModeBins + 17 )  = 2.89d-03
       Inst%FRAC_S( HcoState%MicroPhys%nActiveModeBins + 18 )  = 3.92d-03
       Inst%FRAC_S( HcoState%MicroPhys%nActiveModeBins + 19 )  = 5.26d-03
       Inst%FRAC_S( HcoState%MicroPhys%nActiveModeBins + 20 )  = 7.50d-03
       Inst%FRAC_S( HcoState%MicroPhys%nActiveModeBins + 21 )  = 1.20d-02
       Inst%FRAC_S( HcoState%MicroPhys%nActiveModeBins + 22 )  = 2.08d-02
       Inst%FRAC_S( HcoState%MicroPhys%nActiveModeBins + 23 )  = 3.62d-02
       Inst%FRAC_S( HcoState%MicroPhys%nActiveModeBins + 24 )  = 5.91d-02
       Inst%FRAC_S( HcoState%MicroPhys%nActiveModeBins + 25 )  = 8.74d-02
       Inst%FRAC_S( HcoState%MicroPhys%nActiveModeBins + 26 )  = 1.15d-01
       Inst%FRAC_S( HcoState%MicroPhys%nActiveModeBins + 27 )  = 1.34d-01
       Inst%FRAC_S( HcoState%MicroPhys%nActiveModeBins + 28 )  = 1.37d-01
       Inst%FRAC_S( HcoState%MicroPhys%nActiveModeBins + 29 )  = 1.24d-01
       Inst%FRAC_S( HcoState%MicroPhys%nActiveModeBins + 30 )  = 9.85d-02

# endif

    ELSE

       ! Stop w/ error message
       CALL HCO_ERROR( HcoState%Config%Err, 'Wrong number of TOMAS dust bins!', RC )

    ENDIF

#endif

    !=====================================================================
    ! Activate fields in ExtState used by Ginoux dust
    !=====================================================================

    ! Activate met. fields required by this module
    ExtState%U10M%DoUse    = .TRUE.
    ExtState%V10M%DoUse    = .TRUE.
    ExtState%GWETTOP%DoUse = .TRUE.

    !=======================================================================
    ! Leave w/ success
    !=======================================================================
    IF ( ALLOCATED(SpcNames) ) DEALLOCATE(SpcNames)

    ! Nullify pointers
    Inst    => NULL()

    CALL HCO_LEAVE( HcoState%Config%Err,RC )

  END SUBROUTINE HcoX_DustGinoux_Init
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCOX_DustGinoux_Final
!
! !DESCRIPTION: Subroutine HcoX\_DustGinoux\_Final finalizes the HEMCO
! DUSTGINOUX extension.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HcoX_DustGinoux_Final( ExtState )
!
! !INPUT PARAMETERS:
!
    TYPE(Ext_State),  POINTER       :: ExtState   ! Module options
!
! !REVISION HISTORY:
!  11 Dec 2013 - C. Keller - Now a HEMCO extension
!EOP
!------------------------------------------------------------------------------
!BOC

    !=======================================================================
    ! HCOX_DUSTGINOUX_FINAL begins here!
    !=======================================================================

    CALL InstRemove ( ExtState%DustGinoux )



  END SUBROUTINE HcoX_DustGinoux_Final
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCOX_DustGinoux_GetChDust
!
! !DESCRIPTION: Function HCOX\_DustGinoux\_GetChDust returns the CH\_DUST
! parameter for the current simulation type.
!\\
!\\
! !INTERFACE:
!
  FUNCTION HCOX_DustGinoux_GetChDust( Inst, HcoState ) RESULT( CH_DUST )
!
! !INPUT PARAMETERS:
!
    TYPE(MyInst),    POINTER        :: Inst      ! Instance
    TYPE(HCO_State), POINTER        :: HcoState  ! Hemco state
!
! !RETURN VALUE:
!
    REAL*8 :: CH_DUST
!
! !REMARKS:
!  The logic in the #ifdefs may need to be cleaned up later on.  We have
!  just replicated the existing code in pre-HEMCO versions of dust_mod.F.
!
! !REVISION HISTORY:
!  11 Dec 2013 - C. Keller   - Initial version
!  25 Sep 2014 - R. Yantosca - Updated for TOMAS
!  24 Aug 2017 - M. Sulprizio- Remove support for GRID1x1
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Transfer coeff for type natural source  (kg*s2/m5)
    ! Emission reduction factor for China-nested grid domain (win, 4/27/08)

    IF ( TRIM(HcoState%Config%GridRes)  == '4.0x5.0'  ) THEN

       !-----------------------------------------------------------------------
       ! All 4x5 simulations (including TOMAS)
       !-----------------------------------------------------------------------
       Inst%CH_DUST  = 9.375d-10

    ELSE

       !-----------------------------------------------------------------------
       ! All other resolutions
       !-----------------------------------------------------------------------

       ! Start w/ same value as for 4x5
       Inst%CH_DUST  = 9.375d-10

#if defined( TOMAS )
       ! KLUDGE: For TOMAS simulations at grids higher than 4x5 (e.g. 2x25),
       ! then multiplyCH_DUST by 0.75.  (Sal Farina)
       Inst%CH_DUST  = Inst%CH_DUST * 0.75d0
#endif

    ENDIF

  END FUNCTION HCOX_DustGinoux_GetChDust
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: InstGet
!
! !DESCRIPTION: Subroutine InstGet returns a poiner to the desired instance.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE InstGet ( Instance, Inst, RC, PrevInst )
!
! !INPUT PARAMETERS:
!
    INTEGER                             :: Instance
    TYPE(MyInst),     POINTER           :: Inst
    INTEGER                             :: RC
    TYPE(MyInst),     POINTER, OPTIONAL :: PrevInst
!
! !REVISION HISTORY:
!  18 Feb 2016 - C. Keller   - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
    TYPE(MyInst),     POINTER    :: PrvInst

    !=================================================================
    ! InstGet begins here!
    !=================================================================

    ! Get instance. Also archive previous instance.
    PrvInst => NULL()
    Inst    => AllInst
    DO WHILE ( ASSOCIATED(Inst) )
       IF ( Inst%Instance == Instance ) EXIT
       PrvInst => Inst
       Inst    => Inst%NextInst
    END DO
    IF ( .NOT. ASSOCIATED( Inst ) ) THEN
       RC = HCO_FAIL
       RETURN
    ENDIF

    ! Pass output arguments
    IF ( PRESENT(PrevInst) ) PrevInst => PrvInst

    ! Cleanup & Return
    PrvInst => NULL()
    RC = HCO_SUCCESS

  END SUBROUTINE InstGet
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: InstCreate
!
! !DESCRIPTION: Subroutine InstCreate creates a new instance.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE InstCreate ( ExtNr, Instance, Inst, RC )
!
! !INPUT PARAMETERS:
!
    INTEGER,       INTENT(IN)       :: ExtNr
!
! !OUTPUT PARAMETERS:
!
    INTEGER,       INTENT(  OUT)    :: Instance
    TYPE(MyInst),  POINTER          :: Inst
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,       INTENT(INOUT)    :: RC
!
! !REVISION HISTORY:
!  18 Feb 2016 - C. Keller   - Initial version
!  26 Oct 2016 - R. Yantosca - Don't nullify local ptrs in declaration stmts
!EOP
!------------------------------------------------------------------------------
!BOC
    TYPE(MyInst), POINTER          :: TmpInst
    INTEGER                        :: nnInst

    !=================================================================
    ! InstCreate begins here!
    !=================================================================

    ! ----------------------------------------------------------------
    ! Generic instance initialization
    ! ----------------------------------------------------------------

    ! Initialize
    Inst => NULL()

    ! Get number of already existing instances
    TmpInst => AllInst
    nnInst = 0
    DO WHILE ( ASSOCIATED(TmpInst) )
       nnInst  =  nnInst + 1
       TmpInst => TmpInst%NextInst
    END DO

    ! Create new instance
    ALLOCATE(Inst)
    Inst%Instance = nnInst + 1
    Inst%ExtNr    = ExtNr

    ! Attach to instance list
    Inst%NextInst => AllInst
    AllInst       => Inst

    ! Update output instance
    Instance = Inst%Instance

    ! ----------------------------------------------------------------
    ! Type specific initialization statements follow below
    ! ----------------------------------------------------------------

    ! Return w/ success
    RC = HCO_SUCCESS

  END SUBROUTINE InstCreate
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!BOP
!
! !IROUTINE: InstRemove
!
! !DESCRIPTION: Subroutine InstRemove creates a new instance.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE InstRemove ( Instance )
!
! !INPUT PARAMETERS:
!
    INTEGER                         :: Instance
!
! !REVISION HISTORY:
!  18 Feb 2016 - C. Keller   - Initial version
!  26 Oct 2016 - R. Yantosca - Don't nullify local ptrs in declaration stmts
!EOP
!------------------------------------------------------------------------------
!BOC
    INTEGER                     :: RC
    TYPE(MyInst), POINTER       :: PrevInst
    TYPE(MyInst), POINTER       :: Inst

    !=================================================================
    ! InstRemove begins here!
    !=================================================================

    ! Init
    PrevInst => NULL()
    Inst     => NULL()

    ! Get instance. Also archive previous instance.
    CALL InstGet ( Instance, Inst, RC, PrevInst=PrevInst )

    ! Instance-specific deallocation
    IF ( ASSOCIATED(Inst) ) THEN

       ! Pop off instance from list
       IF ( ASSOCIATED(PrevInst) ) THEN

          ! Free pointer
          IF ( ASSOCIATED( Inst%SRCE_SAND ) ) DEALLOCATE( Inst%SRCE_SAND )
          IF ( ASSOCIATED( Inst%SRCE_SILT ) ) DEALLOCATE( Inst%SRCE_SILT )
          IF ( ASSOCIATED( Inst%SRCE_CLAY ) ) DEALLOCATE( Inst%SRCE_CLAY )

          ! Cleanup option object
          IF ( ASSOCIATED( Inst%IPOINT    ) ) DEALLOCATE( Inst%IPOINT    )
          IF ( ASSOCIATED( Inst%FRAC_S    ) ) DEALLOCATE( Inst%FRAC_S    )
          IF ( ASSOCIATED( Inst%DUSTDEN   ) ) DEALLOCATE( Inst%DUSTDEN   )
          IF ( ASSOCIATED( Inst%DUSTREFF  ) ) DEALLOCATE( Inst%DUSTREFF  )
          IF ( ASSOCIATED( Inst%FLUX      ) ) DEALLOCATE( Inst%FLUX      )
          IF ( ASSOCIATED( Inst%FLUX_ALK  ) ) DEALLOCATE( Inst%FLUX_ALK  )
          IF ( ALLOCATED ( Inst%HcoIDs    ) ) DEALLOCATE( Inst%HcoIDs    )
          IF ( ALLOCATED ( Inst%HcoIDsALK ) ) DEALLOCATE( Inst%HcoIDsALK )

          PrevInst%NextInst => Inst%NextInst
       ELSE
          AllInst => Inst%NextInst
       ENDIF
       DEALLOCATE(Inst)
       Inst => NULL()
    ENDIF

   END SUBROUTINE InstRemove
!EOC
END MODULE HCOX_DustGinoux_Mod
