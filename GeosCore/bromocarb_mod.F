!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: bromocarb_mod
!
! !DESCRIPTION: Module BROMOCARB\_MOD contains variables and routines for the
!  GEOS-CHEM bromocarbon simulation.
!\\
!\\
! !INTERFACE:
!
      MODULE BROMOCARB_MOD
!
! !USES:
!
      IMPLICIT NONE
      PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
!      PUBLIC :: EMIS_CHBr3
!      PUBLIC :: EMIS_CH2Br2
      PUBLIC :: SET_CH3Br
      PUBLIC :: SET_BrO
!      PUBLIC :: SEA_SURFACE
      PUBLIC :: INIT_BROMOCARB
      PUBLIC :: CLEANUP_BROMOCARB
!
! !PUBLIC DATA MEMBERS:
!
      ! For scaling bromine emissions (mpayer, 5/15/12)
!      REAL*8, PUBLIC :: Br_SCALING
!
! !REVISION HISTORY:
!  13 Aug 2007 - J. Parrella - Initial version
!  22 May 2012 - M. Payer    - Added ProTeX headers 
!  27 Aug 2012 - M. Payer    - Now parallelize DO loops
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
      ! Molecular weight of CHBr3  [kg/mol]
      REAL*8,  PARAMETER   :: MWT_CHBr3  = 2.53d-1

      ! Molecular weight of CH2Br2 [kg/mol]
      REAL*8,  PARAMETER   :: MWT_CH2Br2 = 1.74d-1

      ! Number of elements for TOT_OC_AREA
      INTEGER, PARAMETER   :: NUM_REGIONS = 2
!
! PRIVATE TYPES:
!
      ! For TOT_OC_AREA:
      ! If Region = 1, band between 20S and 20N
      !    Region = 2, band between 20  and 50  degrees, NH + SH
!      REAL*8,  ALLOCATABLE :: TOT_OC_AREA(:)
!      REAL*8,  ALLOCATABLE :: FRAC_IN_ab(:,:,:)  ! (region, I, J)

      ! Surface area of the box; only depends on latitude
      REAL*8               :: A_M2

      ! Qing Liang's emissions variables
!      REAL*8,  ALLOCATABLE :: chbr3_emiss(:,:)
!      REAL*8,  ALLOCATABLE :: ch2br2_emiss(:,:)

      LOGICAL, SAVE        :: FIRST_INIT = .TRUE.

      CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: set_ch3br
!
! !DESCRIPTION: Subroutine SET\_CH3Br set CH3Br Concentrations in the planetary
!  boundary layer. Based on latitude bands (1) 90-55N, (2) 55N-0, (3) 0-55S,
!  (4) 55-90S. Values for setting pbl flux were determined by surface
!  measurements from NOAA 2006 data.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE SET_CH3Br( N_TRACERS, TCVV, STT, unit_flag, State_Met )
! 
! !USES:
!
      USE GIGC_State_Met_Mod, ONLY : MetState
      USE GRID_MOD,           ONLY : GET_YMID
      USE PBL_MIX_MOD,        ONLY : GET_FRAC_UNDER_PBLTOP
      USE TRACER_MOD,         ONLY : TRACER_NAME

      USE CMN_SIZE_MOD             ! Size parameters
!
! !INPUT PARAMETERS:
!
      LOGICAL,        INTENT(IN)    :: unit_flag
      INTEGER,        INTENT(IN)    :: N_TRACERS
      REAL*8,         INTENT(IN)    :: TCVV(N_TRACERS) ! [Air MW / Tracer MW]
      TYPE(MetState), INTENT(IN)    :: State_Met       ! Met State object
!
! !INPUT/OUTPUT PARAMETERS:
!
      REAL*8,         INTENT(INOUT) :: STT(IIPAR,JJPAR,LLPAR,N_TRACERS)
!
! !REMARKS: 
!   1) STT is converted back and forth between units of [kg] and
!      [v/v]. Placement of the call to SET_CH3Br in main.f (it's
!      with the emissions) means that it should be in [kg].
!
! !REVISION HISTORY: 
!  12 Feb 2008 - J. Parrella - Initial version
!  22 May 2012 - M. Payer    - Added ProTeX headers
!  27 Aug 2012 - M. Payer    - Added parallel DO loop
!  28 Aug 2012 - M. Payer    - Add error check for CH3Br to avoid OOB error
!  09 Nov 2012 - M. Payer    - Replaced all met field arrays with State_Met
!                              derived type object
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
      ! CH3Br values ( from pptv )
      REAL*8, PARAMETER :: gt55N     = 8.35d-12 
      REAL*8, PARAMETER :: gt0_lt55N = 8.27d-12
      REAL*8, PARAMETER :: lt0_gt55S = 6.94d-12
      REAL*8, PARAMETER :: lt55S     = 6.522d-12
!
! !LOCAL VARIABLES:
!
      INTEGER           :: I,  J,  L,  N
      INTEGER           :: CH3Br_sel
      REAL*8            :: CH3Br_conc
      REAL*8            :: LAT_MID

      ! for testing
      LOGICAL           :: FIRST_COUNT = .TRUE.

      !===================================
      ! Initializing Variables
      !===================================

      LAT_MID    = 0d0
      CH3Br_conc = 0d0

      ! get ID # for CH3Br
      DO N = 1, N_TRACERS
         IF(TRACER_NAME(N) == 'CH3Br') THEN
            CH3Br_sel = N
         ENDIF
      ENDDO

      ! Return if CH3Br is not found
      IF ( CH3Br_sel <= 0 ) THEN
         WRITE( 6, '(a)' ) REPEAT( '=', 79 )
         WRITE( 6, '(a)' ) 'SET_CH3Br: CH3Br not found, so do not'
         WRITE( 6, '(a)' ) 'set concentrations in STT'
         WRITE( 6, '(a)' ) REPEAT( '=', 79 )
         RETURN 
      ENDIF

      ! -----------------------------------------
      ! If we aren't using bromocarbons, then
      ! set the CH3Br equal to zero.
      ! -----------------------------------------
!      IF ( .not. LWARWICK_VSLS ) THEN
!         STT(:,:,:,CH3Br_sel) = 0.d0
!         RETURN
!      ENDIF

      !=================================================================
      ! SET_CH3Br begins here!
      !=================================================================

      ! Loop over grid boxes
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, LAT_MID, CH3Br_conc )
      DO L = 1, LLPAR
      DO J = 1, JJPAR
      DO I = 1, IIPAR
         IF ( GET_FRAC_UNDER_PBLTOP( I, J, L ) > 0d0 ) THEN

            ! base lat band selection on midpoint
            ! latitude of the box
            LAT_MID = GET_YMID( I, J, L )

            ! Selecting the latitude bands:
            IF ( LAT_MID > 55d0 ) THEN
               CH3Br_conc = gt55N
            ELSEIF ( (LAT_MID >= 0d0) .and. (LAT_MID <= 55d0) ) THEN
               CH3Br_conc = gt0_lt55N
            ELSEIF ( (LAT_MID < 0d0) .and. (LAT_MID >= -55d0) ) THEN
               CH3Br_conc = lt0_gt55S
            ELSEIF ( LAT_MID < -55d0 ) THEN
               CH3Br_conc = lt55S
            ENDIF

            ! Make sure we're using the correct units
            IF ( unit_flag ) THEN
               ! if the flag is true, then STT has been
               ! converted from kg/box to v/v mixing ratio.
               ! so we must supply v/v
               STT(I,J,L,CH3Br_sel) = CH3Br_conc
            ELSE
               ! Now convert the [v/v] units to [kg]
               ! as in convert_units subroutine in dao_mod.f
               STT(I,J,L,CH3Br_sel) = CH3Br_conc * State_Met%AD(I,J,L) /
     &                                TCVV(CH3Br_sel)
            ENDIF

!            IF ( FIRST_COUNT ) THEN
!               print *, 'New STT for CH3Br'
!               print *, STT(I,J,L,CH3Br_sel), ' [kg] '
!               print *, CH3Br_conc, ' [v/v] '
!               print *, 'LAT_MID = ', LAT_MID
!               FIRST_COUNT = .FALSE.
!            ENDIF

         ENDIF  ! end selection of PBL boxes

      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      END SUBROUTINE SET_CH3Br
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: set_bro
!
! !DESCRIPTION: Subroutine SET\_BRO sets BrO concentrations in the planetary
!  boundary layer. Based on latitude bands (1) 90-55N, (2) 55N-0, (3) 0-55S,
!  (4) 55-90S. Values for setting pbl flux were determined by surface
!  measurements from NOAA 2006 data.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE SET_BRO( am_I_Root, Input_Opt, State_Met, 
     &                    unit_flag, State_Chm, RC         )

! 
! !USES:
!
      USE CMN_SIZE_MOD
      USE DAO_MOD,            ONLY : IS_WATER
      USE GIGC_ErrCode_Mod
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
      USE GIGC_State_Met_Mod, ONLY : MetState
      USE LOGICAL_MOD,        ONLY : LFIX_PBL_BRO
      USE PBL_MIX_MOD,        ONLY : GET_FRAC_UNDER_PBLTOP
      USE TRACER_MOD,         ONLY : TRACER_NAME

!
! !INPUT PARAMETERS:
!
      LOGICAL,        INTENT(IN)    :: am_I_Root
      TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
      TYPE(MetState), INTENT(IN)    :: State_Met        ! Met State object
      LOGICAL,        INTENT(IN)    :: unit_flag

      REAL*8,         INTENT(IN)    :: TCVV(N_TRACERS)  ! [Air MW / Tracer MW]
!
! !INPUT/OUTPUT PARAMETERS:
!
      REAL*8,         INTENT(INOUT) :: STT(IIPAR,JJPAR,LLPAR,N_TRACERS)
!
! !OUTPUT PARAMETERS:
!
      INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
!
! !REMARKS:
!   1) STT is converted back and forth between units of [kg] and
!      [v/v]. Placement of the call to SET_Bro in main.f (it's
!      with the emissions) means that it should be in [kg].
!
! !REVISION HISTORY:
!  12 Feb 2008 - J. Parrella - Initial version
!  22 May 2012 - M. Payer    - Added ProTeX headers
!  27 Aug 2012 - M. Payer    - Added parallel DO loop
!  28 Aug 2012 - M. Payer    - Add error check for BrO to avoid OOB error
!  09 Nov 2012 - M. Payer    - Replaced all met field arrays with State_Met
!                              derived type object
!  27 Nov 2012 - R. Yantosca - Replace SUNCOS with State_Met%SUNCOS
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER :: I,  J,  L,  N, N_TRACERS
      INTEGER :: IJLOOP
      INTEGER :: BrO_sel
      REAL*8  :: BrO_conc


      REAL*8          :: TCVV(Input_Opt%N_TRACERS)
      REAL*8, POINTER :: STT(:,:,:,:)

      ! Assume success
      RC = GIGC_SUCCESS

      ! -----------------------------------------
      ! If we aren't using this 1pptv experiment
      ! then return without updating STT array.
      ! -----------------------------------------
      IF ( .not. Input_Opt%LFIX_PBL_BRO ) THEN
         RETURN
      ENDIF

      ! Reference fields in derived-type objects
      N_TRACERS = Input_Opt%N_TRACERS
      TCVV      => Input_Opt%TCVV(1:N_TRACERS)
      STT       => State_Chm%TRACERS

      ! ------------------------------------------------------
      ! Set the BrO concentration to 1 pptv inside the PBL.
      ! ------------------------------------------------------
      BrO_conc = 1.0d-12

      ! get ID # for Bro
      DO N = 1, N_TRACERS
         IF(TRACER_NAME(N) == 'BrO') THEN
            BrO_sel = N
         ENDIF
      ENDDO

      ! Return if BrO is not found
      IF ( BrO_sel <= 0 ) THEN
         WRITE( 6, '(a)' ) REPEAT( '=', 79 )
         WRITE( 6, '(a)' ) 'SET_BrO: BrO not found, so do not'
         WRITE( 6, '(a)' ) 'set concentrations in STT'
         WRITE( 6, '(a)' ) REPEAT( '=', 79 )
         RETURN 
      ENDIF

      !=================================================================
      ! SET_BrO begins here!
      !=================================================================

      ! Loop over grid boxes
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, IJLOOP, BrO_conc )
      DO J = 1, JJPAR
         lon_loop: DO I = 1, IIPAR

            ! -----------------------------------------
            ! 1. Determine if we're in the marine
            !    boundary layer. If so, procede,
            !    otherwise, skip. (note, we should NOT
            !    0 the concentration... it can be
            !    contributed from other timesteps and
            !    sources.
            ! -----------------------------------------
            IF ( .not. IS_WATER( I, J, State_Met) ) THEN
               CYCLE lon_loop
            ENDIF

            ! -----------------------------------------
            ! 2. Get the cosine of the SZA to determine
            !    if there's available sunlight for
            !    activation of bromine-chemistry.
            !    If so, set the PBL BrO to 1ppt.
            ! -----------------------------------------
            IJLOOP = ( (J-1) * IIPAR ) + I

            IF ( State_Met%SUNCOS(I,J) > 0.d0 ) THEN
               BrO_conc = 1.0d-12 ! 1pptv if daytime
            ELSE
               BrO_conc = 0.d0    ! 0 otherwise
            ENDIF

            DO L = 1, LLPAR
               IF ( GET_FRAC_UNDER_PBLTOP( I, J, L ) > 0d0 ) THEN

                  ! Make sure we're using the correct units
                  IF ( unit_flag ) THEN
                     ! If the flag is true, then STT has been
                     ! converted from kg/box to v/v mixing ratio.
                     ! so we must supply v/v
                     STT(I,J,L,BrO_sel) = BrO_conc
                  ELSE
                     ! Now convert the [v/v] units to [kg]
                     ! as in convert_units subroutine in dao_mod.f
                     STT(I,J,L,BrO_sel) = BrO_conc * 
     &                                    State_Met%AD(I,J,L) /
     &                                    TCVV(BrO_sel)
                  ENDIF

               ENDIF            ! end selection of PBL boxes

            END DO              !END l-loop over altitudes
         END DO lon_loop        !END i-loop over longitudes
      END DO                    !END j-loop over latitudes
!$OMP END PARALLEL DO

      ! Free pointer
      NULLIFY( STT )

      RETURN

      END SUBROUTINE SET_BRO
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_bromocarb
!
! !DESCRIPTION: Subroutine INIT\_BROMOCARB allocates and zeroes BROMOCARB 
!  module arrays.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE INIT_BROMOCARB( State_Met )
!
! !USES:
!
!      USE DIRECTORY_MOD,      ONLY : DATA_DIR, DATA_DIR_1x1
      USE ERROR_MOD,          ONLY : ALLOC_ERR, DEBUG_MSG
      USE GIGC_State_Met_Mod, ONLY : MetState
!      USE GRID_MOD,           ONLY : GET_AREA_M2
!      USE LOGICAL_MOD,        ONLY : LPRT,      LWARWICK_VSLS
      USE LOGICAL_MOD,        ONLY : LPRT
!      USE REGRID_A2A_MOD

      USE CMN_SIZE_MOD             ! Size parameters
      USE COMODE_LOOP_MOD          ! Avogadro's #, called 'AVG' = 6.02252d+23

!      USE m_netcdf_io_open
!      USE m_netcdf_io_read
!      USE m_netcdf_io_close
!
! !INPUT PARAMETERS: 
!
      TYPE(MetState), INTENT(IN) :: State_Met   ! Meteorology State object
!
! !REVISION HISTORY:
!  13 Aug 2007 - J. Parrella - Initial version
!  22 May 2012 - M. Payer    - Added ProTeX headers
!  24 Aug 2012 - R. Yantosca - DO_REGRID_A2A now reads netCDF input file
!  27 Aug 2012 - M. Payer    - Added parallel DO loop
!  03 Jan 2013 - M. Payer    - Renamed PERAREA to IS_MASS in DO_REGRID_A2A
!  05 Jun 2013 - R. Yantosca - Add st2d, ct2d arrays to avoid having the
!                              code create array temporaries in call to NcRd
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      ! Scalars
      INTEGER            :: I, J, fileID
      INTEGER            :: AS
!      CHARACTER(len=60)  :: vsl_file, fmt, vsldir
!      CHARACTER(len=10)  :: cnlon
!      CHARACTER(len=255) :: LLFILENAME, FILENAME

!      ! Arrays
!      INTEGER            :: st2d(2)
!      INTEGER            :: ct2d(2)
!      REAL*8             :: ARRAY2x25(144,91)

      !=================================================================
      ! INIT_BROMOCARB begins here!
      !=================================================================

!      ! Start and count arrays for netCDF reads
!      st2d = (/ 1,   1  /)
!!      ct2d = (/ 144, 91 /)
!
!      !--------------------------
!      ! Allocate and zero arrays
!      !--------------------------
!
!      ! Total regional areas
!      ALLOCATE( TOT_OC_AREA( NUM_REGIONS ), STAT=AS )
!      IF ( AS /= 0 ) CALL ALLOC_ERR( 'TOT_OC_AREA' )
!      TOT_OC_AREA = 0d0
!
!      ! Fraction of box in a or b regions
!      ALLOCATE( FRAC_IN_ab( NUM_REGIONS, IIPAR, JJPAR ), STAT=AS )
!      IF ( AS /= 0 ) CALL ALLOC_ERR( 'FRAC_IN_ab' )
!      FRAC_IN_ab = 0d0
!
!      ! CHBr3 emissions [kg(CHBr3)/m2/s]
!      ALLOCATE( chbr3_emiss(IIPAR, JJPAR), STAT=AS)
!      IF ( AS /= 0 ) CALL ALLOC_ERR( 'chbr3_emiss' )
!      chbr3_emiss = 0d0
!
!      ! CH2Br2 emissions [kg(CH2Br2)/m2/s]
!      ALLOCATE( ch2br2_emiss(IIPAR, JJPAR), STAT=AS)
!      IF ( AS /= 0 ) CALL ALLOC_ERR( 'ch2br2_emiss' )
!      ch2br2_emiss = 0d0
!
!      ! Do not read bromine files if emissions are off or if Br_SCALING = 0
!      IF ( .not. LWARWICK_VSLs               ) RETURN
!      IF ( .not. ( ABS( Br_SCALING ) > 0d0 ) ) RETURN
!
!      !------------------------------------------
!      ! Call Sea Surface to start filling arrays
!      !------------------------------------------
!
!      CALL SEA_SURFACE( State_Met )
!      IF ( LPRT ) CALL DEBUG_MSG( '### EMISSBROMOCARB: SEA_SURFACE' )
!
!      ! ----------------------------------------------------
!      ! Now read in Qing Liang's CHBr3 and CH2Br2 Emissions
!      ! **
!      !  These emissions are in 2 x 2.5 (lat x lon) and 
!      !  regridded to current model resolution
!      ! ----------------------------------------------------
!
!      ! File with lat/lon edges for regridding
!      LLFILENAME = TRIM( DATA_DIR_1x1) //
!     &     'MAP_A2A_Regrid_201203/MAP_A2A_latlon_geos2x25.nc'
!      
!      ! Input file name
!      FILENAME = TRIM( DATA_DIR_1x1 ) //
!     &     'bromine_201205/VSL_emissions/SLBromine_Emission.2x25.nc'
!      
!      ! Read in CHBr3 emission file ([kg(CHBr3)/m2/s]), regrid
!      call NcOp_Rd( fileID, TRIM( FILENAME ) )
!      call NcRd( ARRAY2x25, fileID, 'CHBr3_emission', st2d, ct2d )
!      CALL DO_REGRID_A2A( LLFILENAME, 144,         91, 
!     &                    ARRAY2x25,  CHBr3_EMISS, IS_MASS=0,
!     &                    netCDF=.TRUE.                       )
!
!
!      ! Read in CH2Br2 emission file ([kg(CH2Br2)/m2/s]), regrid
!      call NcOp_Rd( fileID, TRIM( FILENAME ) )
!      call NcRd( ARRAY2x25, fileID, 'CH2Br2_emission', st2d, ct2d )
!      CALL DO_REGRID_A2A( LLFILENAME, 144,          91,
!     &                    ARRAY2x25,  CH2Br2_EMISS, IS_MASS=0,
!     &                    netCDF=.TRUE.                        )
!
!      ! Close emission file
!      call NcCl( fileID )
!
!      ! Calculating the CHBr3 emission rate [molecules/box/s]
!      ! from Qing Liang's emissions distribution...
!!$OMP PARALLEL DO
!!$OMP+DEFAULT( SHARED )
!!$OMP+PRIVATE( I, J, A_M2 )
!      DO J = 1, JJPAR
!      DO I = 1, IIPAR
!
!         ! Grid box surface area [m2]-- it's only a fn' of latitude
!         A_M2 = GET_AREA_M2( I, J, 1 )
!
!         ! Conversions:
!         ! kg/m2/s ---> molecules/box/second
!         chbr3_emiss(I,J) = chbr3_emiss(I,J) / MWT_CHBr3
!     &                    * AVG * A_M2
!
!         ch2br2_emiss(I,J) = ch2br2_emiss(I,J) / MWT_CH2Br2
!     &                     * AVG * A_M2
!
!      ENDDO
!      ENDDO
!!$OMP END PARALLEL DO

      IF ( LPRT ) CALL DEBUG_MSG( ' ### Post-Reading VSL files' )

      END SUBROUTINE INIT_BROMOCARB
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: cleanup_bromocarb
!
! !DESCRIPTION: Subroutine CLEANUP\_BROMOCARB deallocates all BROMOCARB 
!  module arrays.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE CLEANUP_BROMOCARB
!
! !REVISION HISTORY:
!  13 Aug 2007 - J. Parrella - Initial version
!  22 May 2012 - M. Payer    - Added ProTeX headers 
!EOP
!------------------------------------------------------------------------------
!BOC
      !=================================================================
      ! CLEANUP_BROMOCARB begins here!
      !=================================================================
!      IF ( ALLOCATED( TOT_OC_AREA   ) ) DEALLOCATE( TOT_OC_AREA   )
!      IF ( ALLOCATED( FRAC_IN_ab    ) ) DEALLOCATE( FRAC_IN_ab    )
!      IF ( ALLOCATED( chbr3_emiss   ) ) DEALLOCATE( chbr3_emiss   )
!      IF ( ALLOCATED( ch2br2_emiss  ) ) DEALLOCATE( ch2br2_emiss  )

      END SUBROUTINE CLEANUP_BROMOCARB
!EOC
      END MODULE BROMOCARB_MOD
