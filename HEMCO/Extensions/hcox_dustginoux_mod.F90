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
! !MODULE VARIABLES:
!

      ! Parameters related to dust bins
      INTEGER                   :: ExtNr = -1
      INTEGER, PARAMETER        :: NBINS = 4     ! # of dust bins 
      INTEGER, ALLOCATABLE      :: HcoIDs  (:)   ! HEMCO species IDs
      INTEGER, ALLOCATABLE      :: IPOINT  (:)   ! 1=sand, 2=silt, 3=clay 
      REAL,    ALLOCATABLE      :: FRAC_S  (:)   !  
      REAL,    ALLOCATABLE      :: DUSTDEN (:)   ! dust density     [kg/m3] 
      REAL,    ALLOCATABLE      :: DUSTREFF(:)   ! effective radius [um] 

      ! Transfer coeff
      REAL*8                    :: CH_DUST 

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
      USE HCO_EMISLIST_MOD,  ONLY : EmisList_GetDataArr
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
      REAL*8                 :: W10M,   DEN,    DIAM,   U_TS0, U_TS
      REAL*8                 :: SRCE_P, REYNOL, ALPHA,  BETA
      REAL*8                 :: GAMMA,  CW,     DTSRCE, A_M2,  G
      REAL                   :: DSRC

      ! Source functions (get from HEMCO core) 
      REAL(hp), POINTER      :: SRCE_SAND(:,:) => NULL()
      REAL(hp), POINTER      :: SRCE_SILT(:,:) => NULL()
      REAL(hp), POINTER      :: SRCE_CLAY(:,:) => NULL()

      REAL*8, PARAMETER      :: RHOA     = 1.25d-3

      ! Flux array
      REAL(hp), TARGET       :: FLUX(HcoState%NX,HcoState%NY,NBINS)

      ! For diagnostics
      REAL(hp), POINTER      :: Arr2D(:,:) => NULL()

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

      !=================================================================
      ! Point to DUST source functions 
      !=================================================================

      CALL EmisList_GetDataArr ( am_I_Root, 'GINOUX_SAND', SRCE_SAND, RC )
      IF ( RC /= HCO_SUCCESS ) RETURN
      CALL EmisList_GetDataArr ( am_I_Root, 'GINOUX_SILT', SRCE_SILT, RC )
      IF ( RC /= HCO_SUCCESS ) RETURN
      CALL EmisList_GetDataArr ( am_I_Root, 'GINOUX_CLAY', SRCE_CLAY, RC )
      IF ( RC /= HCO_SUCCESS ) RETURN

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
         M      = IPOINT(N)

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

      ! Free pointer
      SRCE_SAND => NULL()
      SRCE_SILT => NULL()
      SRCE_CLAY => NULL()

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
      FUNCTION HCOX_DustGinoux_GetChDust RESULT ( CH_DUST )
!
! !ARGUMENTS:
!
      REAL*8 :: CH_DUST
!
! !REVISION HISTORY:
!  11 Dec 2013 - C. Keller - Initial version 
!
! !NOTES: 
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      ! Transfer coeff for type natural source  (kg*s2/m5)
!Prior 4/27/08
!      REAL*8, PARAMETER      :: CH_DUST  = 9.375d-10
!Emission reduction factor for China-nested grid domain (win, 4/27/08)
#if   defined( GRID4x5  )
      CH_DUST  = 9.375d-10

#elif defined( GRID1x1  ) && defined( NESTED_CH )
      CH_DUST  = 9.375d-10 * 0.5d0
      ! Note: monthly emission over the China nested-grid domain is about
      !       2 times higher compared to the same domain in 4x5 resolution
      !       Thus applying 1/2  factor to correct the emission.
#else

      CH_DUST  = 9.375d-10
#endif
      
      END FUNCTION HCOX_DustGinoux_GetChDust
!EOC
      END MODULE HCOX_DustGinoux_Mod
!EOM
