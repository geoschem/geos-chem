!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: bromocarb_mod
!
! !DESCRIPTION: Module BROMOCARB\_MOD contains variables and routines for the
!  GEOS-Chem bromocarbon simulation.
!\\
!\\
! !INTERFACE:
!
      MODULE BROMOCARB_MOD
!
! !USES:
!
      USE PRECISION_MOD

      IMPLICIT NONE
      PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
      PUBLIC :: SET_CH3Br
      PUBLIC :: SET_BrO
!
! !REMARKS:
!  HEMCO has made most of this routine obsolete.  HEMCO reads in the
!  data files in netCDF format and tracks them in its data structure.
!
! !REVISION HISTORY:
!  13 Aug 2007 - J. Parrella - Initial version
!  22 May 2012 - M. Payer    - Added ProTeX headers 
!  27 Aug 2012 - M. Payer    - Now parallelize DO loops
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
!  24 Jun 2014 - R. Yantosca - Removed code made obsolete by HEMCO
!  24 Jun 2014 - R. Yantosca - Removed INIT_BROMOCARB, CLEANUP_BROMOCARB
!  06 Nov 2015 - M. Yannetti - Added PRECISION_MOD
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
      ! Molecular weight of CHBr3  [kg/mol]
      REAL(fp),  PARAMETER   :: MWT_CHBr3  = 2.53e-1_fp

      ! Molecular weight of CH2Br2 [kg/mol]
      REAL(fp),  PARAMETER   :: MWT_CH2Br2 = 1.74e-1_fp

      ! Number of elements for TOT_OC_AREA
      INTEGER, PARAMETER   :: NUM_REGIONS = 2
!
! PRIVATE TYPES:
!
      ! Surface area of the box; only depends on latitude
      REAL(fp)               :: A_M2

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
      SUBROUTINE SET_CH3Br( am_I_Root, Input_Opt, State_Met,
     &                      unit_flag, State_Chm, RC         )
! 
! !USES:
!
      USE CMN_SIZE_MOD
      USE GIGC_ErrCode_Mod
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
      USE GIGC_State_Chm_Mod, ONLY : ChmState
      USE GIGC_State_Met_Mod, ONLY : MetState
      USE GRID_MOD,           ONLY : GET_YMID
      USE PBL_MIX_MOD,        ONLY : GET_FRAC_UNDER_PBLTOP
!
! !INPUT PARAMETERS:
!
      LOGICAL,        INTENT(IN)    :: am_I_Root   ! Are we on the root CPU?
      TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
      TYPE(MetState), INTENT(IN)    :: State_Met   ! Meteorology State object
      LOGICAL,        INTENT(IN)    :: unit_flag   ! Convert kg -> v/v?
!
! !INPUT/OUTPUT PARAMETERS:
!
      TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
!
! !OUTPUT PARAMETERS:
!
      INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
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
!  24 Jun 2014 - R. Yantosca - Now accept am_I_Root, Input_Opt, State_Met,
!                              State_Chm.  Get fields from these objects.
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
      ! CH3Br values ( from pptv )
      REAL(fp), PARAMETER :: gt55N     = 8.35e-12_fp 
      REAL(fp), PARAMETER :: gt0_lt55N = 8.27e-12_fp
      REAL(fp), PARAMETER :: lt0_gt55S = 6.94e-12_fp
      REAL(fp), PARAMETER :: lt55S     = 6.522e-12_fp
!
! !LOCAL VARIABLES:
!
      INTEGER           :: I,  J,  L,  N
      INTEGER           :: CH3Br_sel
      REAL(fp)            :: CH3Br_conc
      REAL(fp)            :: LAT_MID

      ! For values from Input_Opt
      INTEGER           :: N_TRACERS
      REAL(fp)            :: TCVV(Input_Opt%N_TRACERS)
      REAL(fp),   POINTER :: STT(:,:,:,:)

      !===================================
      ! Initializing Variables
      !===================================

      ! Get fields from derived-type objects
      N_TRACERS  =  Input_Opt%N_TRACERS
      TCVV       =  Input_Opt%TCVV(1:N_TRACERS)
      STT        => State_Chm%TRACERS

      ! Zero
      LAT_MID    = 0e+0_fp
      CH3Br_conc = 0e+0_fp

      ! get ID # for CH3Br
      DO N = 1, N_TRACERS
         IF( Input_Opt%TRACER_NAME(N) == 'CH3Br') THEN
            CH3Br_sel = N
         ENDIF
      ENDDO

      ! Return if CH3Br is not found
      IF ( CH3Br_sel <= 0 .and. am_I_Root ) THEN
         WRITE( 6, '(a)' ) REPEAT( '=', 79 )
         WRITE( 6, '(a)' ) 'SET_CH3Br: CH3Br not found, so do not'
         WRITE( 6, '(a)' ) 'set concentrations in STT'
         WRITE( 6, '(a)' ) REPEAT( '=', 79 )
         RETURN 
      ENDIF

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
         IF ( GET_FRAC_UNDER_PBLTOP( I, J, L ) > 0e+0_fp ) THEN

            ! base lat band selection on midpoint
            ! latitude of the box
            LAT_MID = GET_YMID( I, J, L )

            ! Selecting the latitude bands:
            IF ( LAT_MID > 55e+0_fp ) THEN
               CH3Br_conc = gt55N
            ELSEIF ( (LAT_MID >= 0e+0_fp) .and. (LAT_MID <= 55e+0_fp) ) 
     &       THEN
               CH3Br_conc = gt0_lt55N
            ELSEIF ( (LAT_MID < 0e+0_fp) .and. (LAT_MID >= -55e+0_fp) )
     &       THEN
               CH3Br_conc = lt0_gt55S
            ELSEIF ( LAT_MID < -55e+0_fp ) THEN
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

         ENDIF  ! end selection of PBL boxes

      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      ! Free pointer
      NULLIFY( STT )

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
      USE GIGC_State_Chm_Mod, ONLY : ChmState
      USE GIGC_State_Met_Mod, ONLY : MetState
      USE PBL_MIX_MOD,        ONLY : GET_FRAC_UNDER_PBLTOP
!
! !INPUT PARAMETERS:
!
      LOGICAL,        INTENT(IN)    :: am_I_Root   ! Are we on the root CPU?
      TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
      TYPE(MetState), INTENT(IN)    :: State_Met   ! Met State object
      LOGICAL,        INTENT(IN)    :: unit_flag   ! Convert kg -> v/v ?
!
! !INPUT/OUTPUT PARAMETERS:
!
      TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
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
!  24 Jun 2014 - R. Yantosca - Remove IJLOOP, it's not used anymroe
!  06 Nov 2014 - M. Yannetti - Changed REAL*8 to REAL(fp)
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      ! Scalars
      INTEGER         :: I,  J,  L,  N, N_TRACERS, BrO_sel
      REAL(fp)          :: BrO_conc

      ! Pointers
      REAL(fp)          :: TCVV(Input_Opt%N_TRACERS)
      REAL(fp), POINTER :: STT(:,:,:,:)

      !====================================================================
      ! SET_BrO begins here!
      !====================================================================

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
      N_TRACERS =  Input_Opt%N_TRACERS
      TCVV      =  Input_Opt%TCVV(1:N_TRACERS)
      STT       => State_Chm%TRACERS

      ! ------------------------------------------------------
      ! Set the BrO concentration to 1 pptv inside the PBL.
      ! ------------------------------------------------------
      BrO_conc = 1.0d-12

      ! get ID # for Bro
      DO N = 1, N_TRACERS
         IF( Input_Opt%TRACER_NAME(N) == 'BrO') THEN
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
!$OMP+PRIVATE( I, J, L, BrO_conc )
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
      END MODULE BROMOCARB_MOD
