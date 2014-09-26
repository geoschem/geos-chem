!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
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
!
! !INTERFACE:
!
      MODULE HCOX_DUSTGINOUX_MOD 
!
! !USES:
!
      USE HCO_ERROR_MOD
      USE HCO_DIAGN_MOD
      USE HCOX_State_MOD,    ONLY : Ext_State
      USE HCO_STATE_MOD,     ONLY : HCO_State 

      IMPLICIT NONE
      PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
      PUBLIC :: HcoX_DustGinoux_Run
      PUBLIC :: HcoX_DustGinoux_Init
      PUBLIC :: HcoX_DustGinoux_Final
      PUBLIC :: HcoX_DustGinoux_GETCHDUST
!
! !REVISION HISTORY:
!  (1 ) Added parallel DO loop in GET_ORO (bmy, 4/14/04)
!  (2 ) Now references "directory_mod.f" (bmy, 7/20/04)
!  (3 ) Fixed typo in ORO_IS_LND for PGI compiler (bmy, 3/1/05)
!  (4 ) Modified for GEOS-5 and GCAP met fields (swu, bmy, 8/16/05)
!  (5 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (6 ) Now uses GOCART source function (tdf, bmy, 1/25/07)
!  (7 ) Modifications for 0.5 x 0.667 grid (yxw, dan, bmy, 11/6/08)
!  (8 ) Updates for nested grids (amv, bmy, 12/18/09)
!  01 Mar 2012 - R. Yantosca - Now reference new grid_mod.F90
!  11 Dec 2013 - C. Keller - Now a HEMCO extension
!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
      INTEGER, PARAMETER   :: NBINS = 4         ! # of dust bins 
!
! !MODULE VARIABLES:
!
      ! Parameters related to dust bins
      INTEGER              :: ExtNr = -1
      INTEGER, ALLOCATABLE :: HcoIDs  (:)       ! HEMCO species IDs
      INTEGER, ALLOCATABLE :: IPOINT  (:)       ! 1=sand, 2=silt, 3=clay 
      REAL,    ALLOCATABLE :: FRAC_S  (:)       !  
      REAL,    ALLOCATABLE :: DUSTDEN (:)       ! dust density     [kg/m3] 
      REAL,    ALLOCATABLE :: DUSTREFF(:)       ! effective radius [um] 

      ! Source functions (get from HEMCO core) 
      REAL(hp), POINTER    :: SRCE_SAND(:,:) => NULL()
      REAL(hp), POINTER    :: SRCE_SILT(:,:) => NULL()
      REAL(hp), POINTER    :: SRCE_CLAY(:,:) => NULL()

      ! Transfer coeff
      REAL*8               :: CH_DUST 

      CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
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
      SUBROUTINE HcoX_DustGinoux_Run( am_I_Root, ExtState, HcoState, RC )
!
! !USES:
!
      USE HCO_EMISLIST_MOD,  ONLY : HCO_GetPtr
      USE HCO_FLUXARR_MOD,   ONLY : HCO_EmisAdd 
!
! !ARGUMENTS:
!
      LOGICAL,         INTENT(IN   )  :: am_I_Root
      TYPE(Ext_State), POINTER        :: ExtState    ! Module options
      TYPE(HCO_State), POINTER        :: HcoState   ! Hemco state 
      INTEGER,         INTENT(INOUT)  :: RC 
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
!  11 Dec 2013 - C. Keller - Now a HEMCO extension
!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER                :: I, J, N, M, tmpID
      LOGICAL                :: ERR
      LOGICAL, SAVE          :: FIRST = .TRUE. 
      REAL*8                 :: W10M,   DEN,    DIAM,   U_TS0, U_TS
      REAL*8                 :: SRCE_P, REYNOL, ALPHA,  BETA
      REAL*8                 :: GAMMA,  CW,     DTSRCE, A_M2,  G
      REAL                   :: DSRC

      REAL*8, PARAMETER      :: RHOA     = 1.25d-3

      ! Flux array
      REAL(hp), TARGET       :: FLUX(HcoState%NX,HcoState%NY,NBINS)

      ! For diagnostics
      REAL(hp), POINTER      :: Arr2D(:,:) => NULL()

#if defined( TOMAS )
      ! Quantities for the TOMAS microphysics simulation
      REAL*8                 :: Mp, Rp
#endif

      !=================================================================
      ! HCOX_DUSTGINOUX_RUN begins here!
      !=================================================================

      ! Return if extension disabled 
      IF ( .NOT. ExtState%DustGinoux ) RETURN

      ! Enter
      CALL HCO_ENTER('HCOX_DustGinoux_Run (hcox_dustginoux_mod.F90)',RC)
      IF ( RC /= HCO_SUCCESS ) RETURN

      ! Set gravity at earth surface (cm/s^2)
      G = HcoState%Phys%g0 * 1.0d2 

      ! Emission timestep [s]
      DTSRCE = GET_TS_EMIS() * 60d0

      !=================================================================
      ! Point to DUST source functions 
      !=================================================================
      IF ( FIRST ) THEN
         CALL HCO_GetPtr ( am_I_Root, 'GINOUX_SAND', SRCE_SAND, RC )
         IF ( RC /= HCO_SUCCESS ) RETURN
         CALL HCO_GetPtr ( am_I_Root, 'GINOUX_SILT', SRCE_SILT, RC )
         IF ( RC /= HCO_SUCCESS ) RETURN
         CALL HCO_GetPtr ( am_I_Root, 'GINOUX_CLAY', SRCE_CLAY, RC )
         IF ( RC /= HCO_SUCCESS ) RETURN
         FIRST = .FALSE.
      ENDIF

      ! Error check
      ERR = .FALSE.

      ! Init
      FLUX = 0.0_hp

!$OMP PARALLEL DO                                            &
!$OMP DEFAULT( SHARED )                                      &
!$OMP PRIVATE( I,      J,     M,    N,      DEN,   DIAM    ) &
!$OMP PRIVATE( REYNOL, ALPHA, BETA, GAMMA,  U_TS0, A_M2    ) &
!$OMP PRIVATE( CW,     U_TS,  W10M, SRCE_P, tmpID, RC      ) &
!$OMP SCHEDULE( DYNAMIC )
      ! Loop over size bins
      DO N = 1, NBINS

         !==============================================================
         ! Threshold velocity as a function of the dust density and the 
         ! diameter from Bagnold (1941), valid for particles larger 
         ! than 10 um.
         ! 
         ! u_ts0 = 6.5*sqrt(dustden(n)*g0*2.*dustreff(n))
         !
         ! Threshold velocity from Marticorena and Bergametti
         ! Convert units to fit dimensional parameters
         !==============================================================
         DEN    = DUSTDEN(N) * 1.d-3            ! [g/cm3]
         DIAM   = 2d0 * DUSTREFF(N) * 1.d2      ! [cm in diameter]
         REYNOL = 1331.d0 * DIAM**(1.56d0) + 0.38d0    ! [Reynolds number]
         ALPHA  = DEN * G * DIAM / RHOA
         BETA   = 1d0 + ( 6.d-3 / ( DEN * G * DIAM**(2.5d0) ) )
         GAMMA  = ( 1.928d0 * REYNOL**(0.092d0) ) - 1.d0

         !==============================================================
         ! I think the 129.d-5 is to put U_TS in m/sec instead of cm/sec
         ! This is a threshold friction velocity!       from M&B
         ! i.e. Ginoux uses the Gillette and Passi formulation
         ! but has substituted Bagnold's Ut with M&B's U*t.
         ! This appears to be a problem.  (tdf, 4/2/04)
         !==============================================================

         ! [m/s] 
         U_TS0  = 129.d-5 * SQRT( ALPHA ) * SQRT( BETA ) / SQRT( GAMMA )
         M = IPOINT(N)

         ! Loop over grid boxes 
         DO J = 1, HcoState%NY 
            DO I = 1, HcoState%NX

!               ! Get grid box surface area [m2]
!               A_M2 = HcoState%Grid%AREA_M2( I, J )

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

               ! 10m wind speed [m/s]
               W10M = ExtState%U10M%Arr%Val(I,J)**2 + &
                      ExtState%V10M%Arr%Val(I,J)**2

               ! Get source function
               IF ( M == 1 ) THEN
                  SRCE_P = SRCE_SAND(I,J)
               ELSEIF ( M == 2 ) THEN
                  SRCE_P = SRCE_SILT(I,J)
               ELSEIF ( M == 2 ) THEN
                  SRCE_P = SRCE_CLAY(I,J)
               ENDIF

               ! Units are m2
               SRCE_P = FRAC_S(N) * SRCE_P !* A_M2

               ! Dust source increment [kg/m2/s]
               FLUX(I,J,N) = CW           * CH_DUST * SRCE_P * W10M &
                           * ( SQRT(W10M) - U_TS )

               ! Not less than zero
               IF ( FLUX(I,J,N) < 0.d0 ) FLUX(I,J,N) = 0.d0

               !========================================================
               ! ND06 diagnostics: dust emissions [kg/timestep]
               !========================================================
!               IF ( ND06 > 0 ) THEN
!                  AD06(I,J,N) = AD06(I,J,N) + ( DSRC * A_M2 )
!               ENDIF
            ENDDO
         ENDDO
      ENDDO
!$OMP END PARALLEL DO

      ! Error check
      IF ( ERR ) THEN
         RC = HCO_FAIL
         RETURN 
      ENDIF

      !=================================================================
      ! PASS TO HEMCO STATE AND UPDATE DIAGNOSTICS 
      !=================================================================
      DO N = 1, NBINS
         IF ( HcoIDs(N) > 0 ) THEN

            ! Add flux to emission array
            CALL HCO_EmisAdd( HcoState, FLUX(:,:,N), HcoIDs(N), RC)
            IF ( RC /= HCO_SUCCESS ) RETURN 

            ! Eventually update diagnostics
            IF ( Diagn_AutoFillLevelDefined(2) ) THEN
               Arr2D => FLUX(:,:,N)
               CALL Diagn_Update( am_I_Root, HcoState, ExtNr=ExtNr, &
                                  Cat=-1, Hier=-1, HcoID=HcoIDs(N), &
                                  AutoFill=1, Array2D=Arr2D, RC=RC   )
               IF ( RC /= HCO_SUCCESS ) RETURN 
               Arr2D => NULL() 
            ENDIF
         ENDIF
      ENDDO !N

      ! Leave w/ success
      CALL HCO_LEAVE ( RC )

      END SUBROUTINE HcoX_DustGinoux_Run
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
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
      SUBROUTINE HcoX_DustGinoux_Init( am_I_Root, HcoState, ExtName, &
                                       ExtState,    RC                )
!
! !USES:
!
      USE HCO_ExtList_Mod,     ONLY : GetExtNr
      USE HCO_STATE_MOD,       ONLY : HCO_GetExtHcoID
!
! !ARGUMENTS:
!
      LOGICAL,          INTENT(IN   )  :: am_I_Root
      TYPE(HCO_State),  POINTER        :: HcoState   ! Hemco state 
      CHARACTER(LEN=*), INTENT(IN   )  :: ExtName    ! Extension name
      TYPE(Ext_State),  POINTER        :: ExtState     ! Module options
      INTEGER,          INTENT(INOUT)  :: RC 

! !REVISION HISTORY:
!  11 Dec 2013 - C. Keller - Now a HEMCO extension
!
! !NOTES: 
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      CHARACTER(LEN=255)             :: MSG
      INTEGER                        :: N, nSpc
      CHARACTER(LEN=31), ALLOCATABLE :: SpcNames(:)
      LOGICAL                        :: verb

      !=================================================================
      ! HCOX_DUSTGINOUX_INIT begins here!
      !=================================================================

      ! Extension Nr.
      ExtNr = GetExtNr( TRIM(ExtName) )
      IF ( ExtNr <= 0 ) RETURN

      ! Enter
      CALL HCO_ENTER('HCOX_DustGinoux_Init (hcox_dustginoux_mod.F90)',RC)
      IF ( RC /= HCO_SUCCESS ) RETURN

      CALL HCO_GetExtHcoID( HcoState, ExtNr, HcoIDs, SpcNames, nSpc, RC)
      IF ( RC /= HCO_SUCCESS ) RETURN

      ! Sanity check
      IF ( nSpc /= NBINS ) THEN
         MSG = 'Ginoux dust model does not have four species!'
         CALL HCO_ERROR ( MSG, RC )
         RETURN
      ENDIF

      ! Get global mass flux tuning factor
      CH_DUST = HcoX_DustGinoux_GetCHDust()

      ! Verbose mode
      IF ( verb ) THEN
         MSG = 'Use Ginoux dust emissions (extension module)'
         CALL HCO_MSG( MSG )

         MSG = 'Use the following species (Name: HcoID):'
         CALL HCO_MSG(MSG)
         DO N = 1, nSpc
            WRITE(MSG,*) TRIM(SpcNames(N)), ':', HcoIDs(N)
            CALL HCO_MSG(MSG)
         ENDDO

         WRITE(MSG,*) 'Global mass flux tuning factor: ', CH_DUST
            CALL HCO_MSG(MSG,SEP2='-')
      ENDIF

      ! Allocate vectors holding bin-specific informations 
      ALLOCATE ( IPOINT  (NBINS) ) 
      ALLOCATE ( FRAC_S  (NBINS) ) 
      ALLOCATE ( DUSTDEN (NBINS) ) 
      ALLOCATE ( DUSTREFF(NBINS) ) 

#if defined( TOMAS )
      !=================================================================
      ! Setup for TOMAS simulations 
      !
      !from ginoux:
      !The U.S. Department of Agriculture (USDA) defines 
      !particles with a radius between 1 um and 25 um as silt, 
      !and below 1 um as clay [Hillel, 1982]. Mineralogical silt 
      !particles are mainly composed of quartz, but they are 
      !often coated with strongly adherent clay such that their 
      !physicochemical properties are similar to clay [Hillel, 
      !1982].
      !
      !SRCE_FUNC Source function                              
      !for 1: Sand, 2: Silt, 3: Clay
      !=================================================================

      ! Give values for IPOINT array (win, 3/11/08)
      IF ( NBINS == ExtState%IBINS ) THEN
         DO N = 1, ExtState%IBINS
            Mp = 1.4*ExtState%Xk(N)
            Rp= ( ( Mp /2500. ) * (3./(4.*HcoState%HcoPhys%PI)))**(0.333)
            IF ( Rp < 1.d-6 ) THEN
               IPOINT(N) = 3
            ELSE
               IPOINT(N) = 2
            END IF
         END DO
      ELSE
         PRINT *,' SRC_DUST_GINOUX: ',NDSTBIN,' bins not supported'
         PRINT *,' SRC_DUST_GINOUX: in this executable.'
         CALL ERROR_STOP('Need to change hard-wired array IPOINT',
     &                   'SRC_DUST_GINOUX: dust_mod.f' )
      ENDIF

      !-------------------------------
      ! Set up dust density
      !-------------------------------
      IF ( NDSTBIN == 4 ) THEN
         Input_Opt%DUSTDEN(1:4) = (/ 2500.d0, 2650.d0, 
     &                               2650.d0, 2650.d0 /)

      ELSE IF ( NDSTBIN == IBINS ) THEN

         DO I = 1, IBINS
            IF ( ExtState%Xk(I) < 4.0D-15 ) THEN
               Input_Opt%DUSTDEN(I)  = 2500.d0
            ELSE
               Input_Opt%DUSTDEN(I)  = 2650.d0
            ENDIF
         ENDDO

      ENDIF

      !--------------------------------
      ! Set up dust effective radius
      !--------------------------------
      IF ( NDSTBIN == 4 ) THEN
         DUSTREFF(1) = 0.73d-6
         DUSTREFF(2) =  1.4d-6
         DUSTREFF(3) =  2.4d-6
         DUSTREFF(4) =  4.5d-6

      ELSE IF ( NDSTBIN == IBINS ) THEN

         ! TOMAS dust Reff (win, 7/17/09)
         DO I = 1, IBINS
            DUSTREFF(I) = 
     &           0.5d0 * ( SQRT(Xk(I) * Xk(I+1)) /
     &           DUSTDEN(I) * 6.d0/PI  )**( 0.333d0 )
         ENDDO 
        
      ENDIF

      !----------------------------------
      ! Set up FRAC_S (only for Ginoux)
      !----------------------------------

      !--------------------------------
      ! Set up dust effective radius
      !--------------------------------
      IF ( NDSTBIN == 4 ) THEN

         ! 4 dust bins
         FRAC_S(1) = 0.095d0
         FRAC_S(2) =   0.3d0
         FRAC_S(3) =   0.3d0
         FRAC_S(4) =   0.3d0

      ELSE IF ( NDSTBIN == ExtState%IBINS ) THEN

         DO I = 1, IBINS
            FRAC_S( I )  = 0.00d-00  !scf initialize
                                     !first few bins in TOMAS15/40 frac_s = 0.0d0
         ENDDO

# if  defined( TOMAS12 ) || defined( TOMAS15 )

         !---------------------------------------------------
         ! TOMAS simulations with 12 or 15 size bins
         !---------------------------------------------------
         FRAC_S( ExtState%ACTMODEBINS + 1  )  = 7.33E-10
         FRAC_S( ExtState%ACTMODEBINS + 2  )  = 2.032E-08
         FRAC_S( ExtState%ACTMODEBINS + 3  )  = 3.849E-07
         FRAC_S( ExtState%ACTMODEBINS + 4  )  = 5.01E-06
         FRAC_S( ExtState%ACTMODEBINS + 5  )  = 4.45E-05
         FRAC_S( ExtState%ACTMODEBINS + 6  )  = 2.714E-04
         FRAC_S( ExtState%ACTMODEBINS + 7  )  = 1.133E-03
         FRAC_S( ExtState%ACTMODEBINS + 8  )  = 3.27E-03
         FRAC_S( ExtState%ACTMODEBINS + 9  )  = 6.81E-03
         FRAC_S( ExtState%ACTMODEBINS + 10 )  = 1.276E-02
         FRAC_S( ExtState%ACTMODEBINS + 11 )  = 2.155E-01
         FRAC_S( ExtState%ACTMODEBINS + 12 )  = 6.085E-01

# else

         !---------------------------------------------------
         ! TOMAS simulations with 30 or 40 size bins
         !---------------------------------------------------        
         FRAC_S( ExtState%ACTMODEBINS +  1 )  = 1.05d-10
         FRAC_S( ExtState%ACTMODEBINS +  2 )  = 6.28d-10
         FRAC_S( ExtState%ACTMODEBINS +  3 )  = 3.42d-09
         FRAC_S( ExtState%ACTMODEBINS +  4 )  = 1.69d-08
         FRAC_S( ExtState%ACTMODEBINS +  5 )  = 7.59d-08
         FRAC_S( ExtState%ACTMODEBINS +  6 )  = 3.09d-07
         FRAC_S( ExtState%ACTMODEBINS +  7 )  = 1.15d-06
         FRAC_S( ExtState%ACTMODEBINS +  8 )  = 3.86d-06
         FRAC_S( ExtState%ACTMODEBINS +  9 )  = 1.18d-05
         FRAC_S( ExtState%ACTMODEBINS + 10 )  = 3.27d-05
         FRAC_S( ExtState%ACTMODEBINS + 11 )  = 8.24d-05
         FRAC_S( ExtState%ACTMODEBINS + 12 )  = 1.89d-04
         FRAC_S( ExtState%ACTMODEBINS + 13 )  = 3.92d-04
         FRAC_S( ExtState%ACTMODEBINS + 14 )  = 7.41d-04
         FRAC_S( ExtState%ACTMODEBINS + 15 )  = 1.27d-03
         FRAC_S( ExtState%ACTMODEBINS + 16 )  = 2.00d-03
         FRAC_S( ExtState%ACTMODEBINS + 17 )  = 2.89d-03
         FRAC_S( ExtState%ACTMODEBINS + 18 )  = 3.92d-03
         FRAC_S( ExtState%ACTMODEBINS + 19 )  = 5.26d-03
         FRAC_S( ExtState%ACTMODEBINS + 20 )  = 7.50d-03
         FRAC_S( ExtState%ACTMODEBINS + 21 )  = 1.20d-02
         FRAC_S( ExtState%ACTMODEBINS + 22 )  = 2.08d-02
         FRAC_S( ExtState%ACTMODEBINS + 23 )  = 3.62d-02
         FRAC_S( ExtState%ACTMODEBINS + 24 )  = 5.91d-02
         FRAC_S( ExtState%ACTMODEBINS + 25 )  = 8.74d-02
         FRAC_S( ExtState%ACTMODEBINS + 26 )  = 1.15d-01
         FRAC_S( ExtState%ACTMODEBINS + 27 )  = 1.34d-01
         FRAC_S( ExtState%ACTMODEBINS + 28 )  = 1.37d-01
         FRAC_S( ExtState%ACTMODEBINS + 29 )  = 1.24d-01
         FRAC_S( ExtState%ACTMODEBINS + 30 )  = 9.85d-02

# endif

      ENDIF

#else
      !=================================================================
      ! Setup for simulations without TOMAS
      !=================================================================

      ! Fill bin-specific information
      IPOINT(1:NBINS)   = (/ 3, 2, 2, 2 /)
      FRAC_S(1:NBINS)   = (/ 0.095d0, 0.3d0, 0.3d0,   0.3d0   /)
      DUSTDEN(1:NBINS)  = (/ 2500.d0, 2650.d0, 2650.d0, 2650.d0 /)
      DUSTREFF(1:NBINS) = (/ 0.73d-6, 1.4d-6, 2.4d-6,  4.5d-6  /)

      ! Activate met. fields required by this module
      ExtState%U10M%DoUse    = .TRUE.
      ExtState%V10M%DoUse    = .TRUE.
      ExtState%GWETTOP%DoUse = .TRUE.

      ! Activate this module
      ExtState%DustGinoux = .TRUE.

#if   defined( TOMAS )


      ! Get the # of activation mode bins
      ActModeBins = HCOX_DustGinoux_GetActModeBins()

#endif

      ! Leave w/ success
      IF ( ALLOCATED(SpcNames) ) DEALLOCATE(SpcNames)
      CALL HCO_LEAVE ( RC ) 

      END SUBROUTINE HcoX_DustGinoux_Init
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
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
      SUBROUTINE HcoX_DustGinoux_Final()
!
! !REVISION HISTORY:
!  11 Dec 2013 - C. Keller - Now a HEMCO extension
!
! !NOTES: 
!EOP
!------------------------------------------------------------------------------
!BOC

      !=================================================================
      ! HCOX_DUSTGINOUX_FINAL begins here!
      !=================================================================
 
      ! Free pointer
      SRCE_SAND => NULL()
      SRCE_SILT => NULL()
      SRCE_CLAY => NULL()

      ! Cleanup option object
      IF ( ALLOCATED(IPOINT  ) ) DEALLOCATE(IPOINT   )
      IF ( ALLOCATED(FRAC_S  ) ) DEALLOCATE(FRAC_S   )
      IF ( ALLOCATED(DUSTDEN ) ) DEALLOCATE(DUSTDEN  )
      IF ( ALLOCATED(DUSTREFF) ) DEALLOCATE(DUSTREFF )

      END SUBROUTINE HcoX_DustGinoux_Final
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
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
      FUNCTION HCOX_DustGinoux_GetChDust() RESULT( CH_DUST )
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
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      ! Transfer coeff for type natural source  (kg*s2/m5)
      ! Emission reduction factor for China-nested grid domain (win, 4/27/08)

#if defined( GRID4x5 )

      !---------------------------------------------------------------------
      ! All 4x5 simulations (including TOMAS)
      !---------------------------------------------------------------------
      CH_DUST  = 9.375d-10

#elif defined( GRID1x1 ) && defined( NESTED_CH )

      !---------------------------------------------------------------------
      ! Note: monthly emission over the China nested-grid domain is about
      !       2 times higher compared to the same domain in 4x5 resolution
      !       Thus applying 1/2  factor to correct the emission.
      !
      !%%% NOTE: MAY NEED TO UPDATE THIS STATEMENT FOR HIGHER RESOLUTION
      !%%% NESTED GRIDS.  THIS WAS ORIGINALLY DONE FOR THE GEOS-3 1x1
      !%%% NESTED GRID.  LOOK INTO THIS LATER.  (bmy, 9/25/14)
      !---------------------------------------------------------------------
      CH_DUST  = 9.375d-10 * 0.5d0

#else

      !---------------------------------------------------------------------
      ! All other resolutions
      !---------------------------------------------------------------------

      ! Start w/ same value as for 4x5
      CH_DUST  = 9.375d-10

#if defined( TOMAS )
      ! KLUDGE: For TOMAS simulations at grids higher than 4x5 (e.g. 2x25),
      ! then multiplyCH_DUST by 0.75.  (Sal Farina)
      CH_DUST  = CH_DUST * 0.75d0
#endif

#endif
      
      END FUNCTION HCOX_DustGinoux_GetChDust
!EOC
      END MODULE HCOX_DustGinoux_Mod
