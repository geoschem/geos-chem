!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: sfcvmr_mod.F90
!
! !DESCRIPTION: Module sfcvmr\_mod.F90 is a simple module which forces
!  surface concentrations of relevant species to pre-determined values.
!\\
!\\
! !INTERFACE:
!
MODULE SFCVMR_MOD
!
! !USES:
!
    USE PhysConstants       ! Physical constants
    USE inquireMod,    ONLY : findFreeLUN
    USE PRECISION_MOD       ! For GEOS-Chem Precision (fp)

    IMPLICIT NONE
    PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
    PUBLIC :: fixSfcVMR
    PUBLIC :: CLEANUP_SfcVMR
    PUBLIC :: INIT_SfcVMR
!
! !REVISION HISTORY:
!  24 Dec 2016 - S. D. Eastham - Initial version.
!  13 Mar 2019 - T. Sherwen - Updates to use CMIP6 monthly Organo-halogens
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !PRIVATE TYPES:
!
    ! Pointers to fields in the HEMCO data structure.
    ! These need to be declared REAL(f4), aka REAL*4.
    REAL(f4), POINTER :: HCO_CH3BR     (:,:) => NULL()
    REAL(f4), POINTER :: HCO_CH3Cl     (:,:) => NULL()
    REAL(f4), POINTER :: HCO_CH2Cl2    (:,:) => NULL()
    REAL(f4), POINTER :: HCO_CHCl3     (:,:) => NULL()

    CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: FIXSFCVMR
!
! !DESCRIPTION: Subroutine FIXSFCVMR fixes the VMR of selected species
! throughout the PBL to observed values.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE fixSfcVMR( am_I_Root, Input_Opt, State_Met, State_Chm, RC )
!
! !USES:
!
    USE CMN_SIZE_MOD
    USE ErrCode_Mod
    USE ERROR_MOD,          ONLY : ERROR_STOP
    USE Input_Opt_Mod,      ONLY : OptInput
    USE State_Met_Mod,      ONLY : MetState
    USE State_Chm_Mod,      ONLY : ChmState
    USE State_Chm_Mod,      ONLY : Ind_
    USE GC_GRID_MOD,        ONLY : GET_YMID
    USE TIME_MOD,           ONLY : Get_Month
    USE HCO_EmisList_Mod,   ONLY : HCO_GetPtr
    USE HCO_INTERFACE_MOD,  ONLY : HcoState
    USE HCO_STATE_MOD,      ONLY : HCO_GetHcoID

    ! Needed for the new CHxCly boundary condition
    USE PBL_MIX_MOD,        ONLY : GET_FRAC_UNDER_PBLTOP
    Use PhysConstants,      Only : AirMW
!
! !INPUT PARAMETERS:
!
    LOGICAL,          INTENT(IN   )  :: am_I_Root  ! root CPU?
    TYPE(MetState),   INTENT(IN   )  :: State_Met  ! Met state
    TYPE(ChmState),   INTENT(INOUT)  :: State_Chm  ! Chemistry State object
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(OptInput),   INTENT(INOUT)  :: Input_Opt  ! Input opts
    INTEGER,          INTENT(INOUT)  :: RC         ! Failure or success
!
! !REVISION HISTORY:
!  27 Aug 2014 - C. Keller   - Initial version
!  16 Jun 2016 - J. Sheng    - Added tracer index retriever
!  20 Jun 2016 - R. Yantosca - Now define species IDs only in the INIT phase
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Loop indices
    Integer            :: I, J, L, MONTH
    ! Species index
    Integer            :: id_Spc
    ! Local variables for quantities from Input_Opt
    LOGICAL            :: LHCodedOrgHal
    LOGICAL            :: LCMIP6OrgHal
    ! Saved (main) scalars
    LOGICAL,          SAVE      :: FIRST         = .TRUE.
    INTEGER,          SAVE      :: id_CH3Br      = -1
    INTEGER,          SAVE      :: id_CH3Cl      = -1
    INTEGER,          SAVE      :: id_CH2Cl2     = -1
    INTEGER,          SAVE      :: id_CHCl3      = -1
    ! Saved (HEMCO) scalars
    INTEGER,          SAVE      :: IDCH3Br       = -1
    INTEGER,          SAVE      :: IDCH3Cl       = -1
    INTEGER,          SAVE      :: IDCH2Cl2      = -1
    INTEGER,          SAVE      :: IDCHCl3       = -1

    Real(fp)          :: A3090S (12), A0030S (12), A0030N(12)
    Real(fp)          :: A3060N (12), A6090N (12)
    Real(fp)          :: YLAT, PPT

    ! Strings
    CHARACTER(LEN=255)         :: LOC = 'fixSfcVMR (sfcvmr_mod.F90)'

    ! Pointer to the species array
    Real(fp),    Pointer       :: Spc(:,:,:,:)
    ! Other pointers
!     REAL(fp), PUBLIC, ALLOCATABLE   :: FIXED_CH3BR(   :,: )
!     REAL(fp), PUBLIC, ALLOCATABLE   :: FIXED_CH3Cl(   :,: )
!     REAL(fp), PUBLIC, ALLOCATABLE   :: FIXED_CH2Cl2(  :,: )
!     REAL(fp), PUBLIC, ALLOCATABLE   :: FIXED_CHCl3(   :,: )

    ! Check that species units are in kg/kg dry air (ewl, 9/10/15)
    IF ( TRIM( State_Chm%Spc_Units ) /= 'kg/kg dry' ) THEN
        CALL GC_Error( 'Incorrect species units: ' //         &
                      State_Chm%Spc_Units, RC,                &
                      'Routine SET_CH3Br in bromocarb_mod.F' )
        RETURN
    ENDIF



    !=================================================================
    ! FIXSFCVMR begins here!
    !=================================================================

    ! Assume success
    RC        = GC_SUCCESS

    ! Get a pointer to the species array
    Spc => State_Chm%Species

    ! Copy booleans from INPUT_OPT
    LHCodedOrgHal     = Input_Opt%LHCodedOrgHal
    LCMIP6OrgHal      = Input_Opt%LCMIP6OrgHal

    ! Import emissions from HEMCO (through HEMCO state)
    IF ( .NOT. ASSOCIATED(HcoState) ) THEN
        CALL ERROR_STOP ( 'HcoState not defined!', LOC )
    ENDIF

    WRITE(*,*) 'TMS CHECK - LCMIP6OrgHal', Input_Opt%LCMIP6OrgHal
    WRITE(*,*) 'TMS CHECK - LHCodedOrgHal', Input_Opt%LHCodedOrgHal

    !=================================================================
    ! Tracer ID setup and error checks (first-time only)
    !=================================================================
    IF ( FIRST ) THEN

        !--------------------------------------------------------------
        ! Query main cod for all of the tracer IDs
        !--------------------------------------------------------------
        id_CH3Br  = Ind_('CH3Br')
        id_CH3Cl  = Ind_('CH3Cl')
        id_CH2Cl2 = Ind_('CH2Cl2')
        id_CHCl3  = Ind_('CHCl3')
        ! Check these...
        WRITE(*,*) 'TMS CHECK - ID id_CH3Br:', id_CH3Br
        WRITE(*,*) 'TMS CHECK - ID id_CH3Cl:', id_CH3Cl
        WRITE(*,*) 'TMS CHECK - ID id_CH2Cl2:', id_CH2Cl2
        WRITE(*,*) 'TMS CHECK - ID id_CHCl3:', id_CHCl3

        ! Return if CH3Br is not found
        IF ( id_CH3Br <= 0 .and. am_I_Root ) THEN
            WRITE( 6, '(a)' ) REPEAT( '=', 79 )
            WRITE( 6, '(a)' ) 'SET_CH3Br: CH3Br not found, so do not'
            WRITE( 6, '(a)' ) 'set concentrations in Spc'
            WRITE( 6, '(a)' ) REPEAT( '=', 79 )
            RETURN
        ENDIF

        !--------------------------------------------------------------
        ! Query HEMCO for all of the tracer IDs
        !--------------------------------------------------------------
        IDCH3Br   = HCO_GetHcoID( 'CH3Br',   HcoState )
        IDCH3Cl   = HCO_GetHcoID( 'CH3Cl',   HcoState )
        IDCH2Cl2  = HCO_GetHcoID( 'CH2Cl2',  HcoState )
        IDCHCl3   = HCO_GetHcoID( 'CHCl3',   HcoState )
        ! Check these...
        WRITE(*,*) 'TMS CHECK - HcoID IDCH3Br:', IDCH3Br
        WRITE(*,*) 'TMS CHECK - HcoID IDCH3Cl:', IDCH3Cl
        WRITE(*,*) 'TMS CHECK - HcoID IDCH2Cl2:', IDCH2Cl2
        WRITE(*,*) 'TMS CHECK - HcoID IDCHCl3:', IDCHCl3

        !--------------------------------------------------------------
        ! Add error checks.  Check to see that the tracer index
        ! and the corresponding field in the HEMCO state are defined.
        ! Move this outside the parallel loop to avoid problems when
        ! exiting while w/in a parallel region. (bmy, 4/17/15)
        !--------------------------------------------------------------
        IF ( Input_Opt%LCMIP6OrgHal ) THEN
            IF ( IDCH3Br < 1 ) THEN
                CALL ERROR_STOP( 'IDCH3Br not defined', LOC )
            ENDIF
!             IF ( .not. ASSOCIATED(HcoState%Spc(IDCH3Br)%Emis%Val) ) THEN
!                 CALL ERROR_STOP('fixed CH3Br emissions not defined', LOC )
!             ENDIF
        ENDIF

        IF ( Input_Opt%LCMIP6OrgHal ) THEN
            IF ( IDCH3Cl < 1 ) THEN
                CALL ERROR_STOP( 'IDCH3Cl not defined', LOC )
            ENDIF
!             IF ( .not. ASSOCIATED(HcoState%Spc(IDCH3Cl)%Emis%Val) ) THEN
!                 CALL ERROR_STOP('fixed CH3Cl emissions not defined', LOC )
!             ENDIF
        ENDIF

        IF ( Input_Opt%LCMIP6OrgHal ) THEN
            IF ( IDCH2Cl2 < 1 ) THEN
                CALL ERROR_STOP( 'IDCH2Cl2 not defined', LOC )
            ENDIF
!             IF ( .not. ASSOCIATED(HcoState%Spc(IDCH2Cl2)%Emis%Val) ) THEN
!                 CALL ERROR_STOP('fixed CH2Cl2 emissions not defined', LOC )
!             ENDIF
        ENDIF

        IF ( Input_Opt%LCMIP6OrgHal ) THEN
            IF ( IDCHCl3 < 1 ) THEN
                CALL ERROR_STOP( 'IDCHCl3 not defined', LOC )
            ENDIF
!             IF ( .not. ASSOCIATED(HcoState%Spc(IDCHCl3)%Emis%Val) ) THEN
!                 CALL ERROR_STOP('fixed CHCl3 emissions not defined', LOC )
!             ENDIF
        ENDIF

        !--------------------------------------------------------------
        ! Also test if the biofuel & biomass tracers are on
        ! (i.e. if the pointers in the HEMCO state are associated).
        ! We can do this just once, outside of the parallel loop
        !--------------------------------------------------------------
        !     BIOFUEL_ON = ASSOCIATED( HcoState%Spc(IDbf)%Emis%Val )
        !     BIOMASS_ON = ASSOCIATED( HcoState%Spc(IDbb)%Emis%Val )

        !--------------------------------------------------------------
        ! Error check: For now, the emission grid must be
        ! on the simulation grid.
        !--------------------------------------------------------------
        IF ( HcoState%NX /= IIPAR .OR.   &
             Hcostate%NY /= JJPAR .OR.   &
             Hcostate%NZ /= LLPAR      ) THEN
            CALL ERROR_STOP( 'HEMCO grid not same as sim. grid!', LOC )
        ENDIF

        ! Set first-time flag to false
        FIRST = .FALSE.
    ENDIF


    ! Now set fixed "emissions" (mixing ratios) in arrays
    IF ( Input_Opt%LHCodedOrgHal ) THEN

        WRITE(*,*) 'TMS CHECK - Using hardcoded emiss. for fixed CHyXz'

        ! Current month
        MONTH = Get_Month()

        YLAT = 0.0
        PPT = 0.0_fp

        ! ---------------------------------------------------
        ! JAS, 9/17/15:
        ! Set mixing ratio of CH3Cl, CH2Cl2, and CHCl3 in PBL
        ! SDE 2016-12-14: This now replaces the UCX routine
        ! which previously gave only CH3Cl.
        ! ---------------------------------------------------
        ! Set CH3Cl mixing ratio in PBL
        id_Spc = Ind_('CH3Cl')

        A3090S = (/514,512,512,509,512,521,530,533,536,534,531,527/)
        A0030S = (/549,549,553,564,574,581,587,575,564,556,549,544/)
        A0030N = (/605,610,609,613,607,580,554,537,539,557,570,579/)
        A3060N = (/573,572,592,598,590,572,545,527,512,532,547,554/)
        A6090N = (/532,547,574,582,572,539,495,464,466,488,503,515/)

        IF ( id_Spc > 0 ) THEN
!$OMP PARALLEL DO                                                 &
!$OMP DEFAULT( SHARED )                                           &
!$OMP PRIVATE( I, J, L, YLAT, PPT)
           DO L = 1, LLPAR
           DO J = 1, JJPAR
           DO I = 1, IIPAR
              YLAT = GET_YMID( I, J, L )
              IF ( YLAT < -30.0_fp ) THEN
                 PPT = A3090S(MONTH)
              ELSE IF ( YLAT >= -30.0_fp .and. YLAT < 0.0_fp  ) THEN
                 PPT = A0030S(MONTH)
              ELSE IF ( YLAT >=   0.0_fp .and. YLAT < 30.0_fp ) THEN
                 PPT = A0030N(MONTH)
              ELSE IF ( YLAT >=  30.0_fp .and. YLAT < 60.0_fp ) THEN
                 PPT = A3060N(MONTH)
              ELSE
                 PPT = A6090N(MONTH)
              ENDIF
              IF (GET_FRAC_UNDER_PBLTOP(I,J,L)>0e+0_fp) THEN
                 Spc(I,J,L,id_Spc) = PPT*1e-12_fp / ( AIRMW / &
                    State_Chm%SpcData(id_Spc)%Info%emMW_g )
              ENDIF  ! end selection of PBL boxes
           ENDDO
           ENDDO
           ENDDO
!$OMP END PARALLEL DO
        ENDIF

        ! Set CH2Cl2 mixing ratio in PBL
        id_Spc = Ind_('CH2Cl2')

        A3090S = (/10,10,10,11,12,12,13,14,13,13,12,11/)
        A0030S = (/15,17,17,18,18,18,19,20,19,18,17,18/)
        A0030N = (/45,45,44,39,45,42,39,34,31,29,35,46/)
        A3060N = (/54,56,56,58,59,52,46,41,41,45,50,54/)
        A6090N = (/65,66,63,63,61,56,51,47,44,47,54,61/)

        IF ( id_Spc > 0 ) THEN
!$OMP PARALLEL DO                                                 &
!$OMP DEFAULT( SHARED )                                           &
!$OMP PRIVATE( I, J, L, YLAT, PPT)
           DO L = 1, LLPAR
           DO J = 1, JJPAR
           DO I = 1, IIPAR
              YLAT = GET_YMID( I, J, L )
              IF ( YLAT < -30.0_fp ) THEN
                 PPT = A3090S(MONTH)
              ELSE IF ( YLAT >= -30.0_fp .and. YLAT < 0.0_fp  ) THEN
                 PPT = A0030S(MONTH)
              ELSE IF ( YLAT >=   0.0_fp .and. YLAT < 30.0_fp ) THEN
                 PPT = A0030N(MONTH)
              ELSE IF ( YLAT >=  30.0_fp .and. YLAT < 60.0_fp ) THEN
                 PPT = A3060N(MONTH)
              ELSE
                 PPT = A6090N(MONTH)
              ENDIF
              IF (GET_FRAC_UNDER_PBLTOP(I,J,L)>0e+0_fp) THEN
                 Spc(I,J,L,id_Spc) = PPT*1e-12_fp / ( AIRMW / &
                    State_Chm%SpcData(id_Spc)%Info%emMW_g )
              ENDIF  ! end selection of PBL boxes
           ENDDO
           ENDDO
           ENDDO
!$OMP END PARALLEL DO
        ENDIF

        ! Set CHCl3 mixing ratio in PBL
        id_Spc = Ind_('CHCl3')

        A3090S = (/5, 5, 5, 5, 6, 6, 7, 7, 7, 6, 6, 6 /)
        A0030S = (/5, 5, 5, 6, 6, 6, 6, 6, 6, 5, 5, 5 /)
        A0030N = (/10,10,10,8, 9, 9, 9, 8, 8, 7, 8, 10/)
        A3060N = (/14,14,14,14,14,12,11,12,14,14,15,14/)
        A6090N = (/16,15,15,14,14,13,13,13,13,14,15,16/)

        IF ( id_Spc > 0 ) THEN
!$OMP PARALLEL DO                                                 &
!$OMP DEFAULT( SHARED )                                           &
!$OMP PRIVATE( I, J, L, YLAT )
           DO L = 1, LLPAR
           DO J = 1, JJPAR
           DO I = 1, IIPAR
              YLAT = GET_YMID( I, J, L )
              IF ( YLAT < -30.0_fp ) THEN
                 PPT = A3090S(MONTH)
              ELSE IF ( YLAT >= -30.0_fp .and. YLAT < 0.0_fp  ) THEN
                 PPT = A0030S(MONTH)
              ELSE IF ( YLAT >=   0.0_fp .and. YLAT < 30.0_fp ) THEN
                 PPT = A0030N(MONTH)
              ELSE IF ( YLAT >=  30.0_fp .and. YLAT < 60.0_fp ) THEN
                 PPT = A3060N(MONTH)
              ELSE
                 PPT = A6090N(MONTH)
              ENDIF
              IF (GET_FRAC_UNDER_PBLTOP(I,J,L)>0e+0_fp) THEN
                 Spc(I,J,L,id_Spc) = PPT*1e-12_fp / ( AIRMW / &
                    State_Chm%SpcData(id_Spc)%Info%emMW_g )
              ENDIF  ! end selection of PBL boxes
           ENDDO
           ENDDO
           ENDDO
!$OMP END PARALLEL DO
        ENDIF

    ! Use the
    ELSE IF ( Input_Opt%LCMIP6OrgHal ) THEN

        WRITE(*,*) 'TMS CHECK - Using CMIP6 emissions from NetCDF'

        !=================================================================
        ! Read in monthly fixed mixing ratio fields
        !=================================================================

        ! Fill an allocatable array with CMIP6 fixed CH3Br fixed field read in by HEMCO
        ! This listed in the NON-EMISSIONS DATA section of the HEMCO config. file.
!         CALL HCO_GetPtr( am_I_Root, HcoState, 'CMIP6_CH3Br_FIXED', &
!                         FIXED_CH3Br, RC )
!         IF ( RC /= GC_SUCCESS ) THEN
!             CALL ERROR_STOP ( 'CMIP6 CH3Br fixed conc not defined', LOC )
!         ENDIF
!
!         ! Get a pointer to the CMIP6 fixed CH3Cl fixed field read in by HEMCO
!         ! This listed in the NON-EMISSIONS DATA section of the HEMCO config. file.
!         CALL HCO_GetPtr( am_I_Root, HcoState, 'CMIP6_CH3Cl_FIXED', &
!                         FIXED_CH3Cl, RC )
!         IF ( RC /= GC_SUCCESS ) THEN
!             CALL ERROR_STOP ( 'CMIP6 CH3Cl fixed conc not defined', LOC )
!         ENDIF
!
!         ! Get a pointer to the CMIP6 fixed CH2Cl2 fixed field read in by HEMCO
!         ! This listed in the NON-EMISSIONS DATA section of the HEMCO config. file.
!         CALL HCO_GetPtr( am_I_Root, HcoState, 'CMIP6_CH2Cl2_FIXED', &
!                          FIXED_CH2Cl2, RC )
!         IF ( RC /= GC_SUCCESS ) THEN
!             CALL ERROR_STOP ( 'CMIP6 CH2Cl2 fixed conc not defined', LOC )
!         ENDIF
!
!         ! Get a pointer to the CMIP6 fixed CH2Cl2 fixed field read in by HEMCO
!         ! This listed in the NON-EMISSIONS DATA section of the HEMCO config. file.
!         CALL HCO_GetPtr( am_I_Root, HcoState, 'CMIP6_CH2Cl2_FIXED', &
!                         FIXED_CH2Cl2, RC )
!         IF ( RC /= GC_SUCCESS ) THEN
!             CALL ERROR_STOP ( 'CMIP6 CH2Cl2 fixed conc not defined', LOC )
!         ENDIF


        ! Get a pointer to the CMIP6 fixed CH3Br fixed field read in by HEMCO
        ! This listed in the NON-EMISSIONS DATA section of the HEMCO config. file.
        CALL HCO_GetPtr( am_I_Root, HcoState, 'CMIP6_CH3Br_FIXED', &
                        HCO_CH3BR, RC )
        IF ( RC /= GC_SUCCESS ) THEN
            CALL ERROR_STOP ( 'CMIP6 CH3Br fixed conc not defined', LOC )
        ENDIF

        ! Get a pointer to the CMIP6 fixed CH3Cl fixed field read in by HEMCO
        ! This listed in the NON-EMISSIONS DATA section of the HEMCO config. file.
        CALL HCO_GetPtr( am_I_Root, HcoState, 'CMIP6_CH3Cl_FIXED', &
                        HCO_CH3Cl, RC )
        IF ( RC /= GC_SUCCESS ) THEN
            CALL ERROR_STOP ( 'CMIP6 CH3Cl fixed conc not defined', LOC )
        ENDIF

        ! Get a pointer to the CMIP6 fixed CH2Cl2 fixed field read in by HEMCO
        ! This listed in the NON-EMISSIONS DATA section of the HEMCO config. file.
        CALL HCO_GetPtr( am_I_Root, HcoState, 'CMIP6_CH2Cl2_FIXED', &
                         HCO_CH2Cl2, RC )
        IF ( RC /= GC_SUCCESS ) THEN
            CALL ERROR_STOP ( 'CMIP6 CH2Cl2 fixed conc not defined', LOC )
        ENDIF

        ! Get a pointer to the CMIP6 fixed CH2Cl2 fixed field read in by HEMCO
        ! This listed in the NON-EMISSIONS DATA section of the HEMCO config. file.
        CALL HCO_GetPtr( am_I_Root, HcoState, 'CMIP6_CHCl3_FIXED', &
                        HCO_CHCl3, RC )
        IF ( RC /= GC_SUCCESS ) THEN
            CALL ERROR_STOP ( 'CMIP6 CHCl3 fixed conc not defined', LOC )
        ENDIF

      !=================================================================
      ! Update concentration in the PBL.
      !=================================================================
      ! Loop over grid boxes
      DO L = 1, LLPAR
      DO J = 1, JJPAR
      DO I = 1, IIPAR
         IF ( GET_FRAC_UNDER_PBLTOP( I, J, L ) > 0e+0_fp ) THEN

            ! --- Use pointers to get HEMCO values
            ! Convert the [v/v] CH3Br conc units (dry air) to [kg/kg]
            ! when setting species concentration. Note: also pptv=>v/v (/1E12).
            Spc(I,J,L,id_CH3Br) = HCO_CH3Br(I,J)  / 1e12_fp / ( AIRMW / &
                              State_Chm%SpcData(id_CH3Br)%Info%emMW_g )

            ! Convert the [v/v] CH3Cl conc units (dry air) to [kg/kg]
            ! when setting species concentration. Note: also pptv=>v/v (/1E12).
            Spc(I,J,L,id_CH3Cl) = HCO_CH3Cl(I,J) / 1e12_fp / ( AIRMW / &
                              State_Chm%SpcData(id_CH3Cl)%Info%emMW_g )

            ! Convert the [v/v] CH3Cl conc units (dry air) to [kg/kg]
            ! when setting species concentration. Note: also pptv=>v/v (/1E12).
            Spc(I,J,L,id_CH2Cl2) = HCO_CH2Cl2(I,J) / 1e12_fp / ( AIRMW / &
                              State_Chm%SpcData(id_CH2Cl2)%Info%emMW_g )

            ! Convert the [v/v] CH3Cl conc units (dry air) to [kg/kg]
            ! when setting species concentration. Note: also pptv=>v/v (/1E12).
            Spc(I,J,L,id_CHCl3) = HCO_CHCl3(I,J) / 1e12_fp / ( AIRMW / &
                              State_Chm%SpcData(id_CHCl3)%Info%emMW_g )
            ! --- Use allocatable arrays to hold HEMCO values

!             Convert the [v/v] CH3Br conc units (dry air) to [kg/kg]
!             when setting species concentration. Note: also pptv=>v/v (/1E12).
!             Spc(I,J,L,id_CH3Br) = FIXED_CH3Br(I,J)  / 1e12_fp / ( AIRMW / &
!                               State_Chm%SpcData(id_CH3Br)%Info%emMW_g )
!
!             Convert the [v/v] CH3Cl conc units (dry air) to [kg/kg]
!             when setting species concentration. Note: also pptv=>v/v (/1E12).
!             Spc(I,J,L,id_CH3Cl) = FIXED_CH3Cl(I,J) / 1e12_fp / ( AIRMW / &
!                               State_Chm%SpcData(id_CH3Cl)%Info%emMW_g )
!
!             Convert the [v/v] CH3Cl conc units (dry air) to [kg/kg]
!             when setting species concentration. Note: also pptv=>v/v (/1E12).
!             Spc(I,J,L,id_CH2Cl2) = FIXED_CH2Cl2(I,J) / 1e12_fp / ( AIRMW / &
!                               State_Chm%SpcData(id_CH2Cl2)%Info%emMW_g )
!
!             Convert the [v/v] CH3Cl conc units (dry air) to [kg/kg]
!             when setting species concentration. Note: also pptv=>v/v (/1E12).
!             Spc(I,J,L,id_CHCl3) = FIXED_CHCl3(I,J) / 1e12_fp / ( AIRMW / &
!                               State_Chm%SpcData(id_CHCl3)%Info%emMW_g )

         ENDIF  ! end selection of PBL boxes

      ENDDO
      ENDDO
      ENDDO

    ELSE
        WRITE(*,*) 'Please select either hardcored of CMIP monthly emissions!'
        WRITE(*,*) 'LHCodedOrgHal = ', LHCodedOrgHal, Input_Opt%LHCodedOrgHal
        WRITE(*,*) 'LCMIP6OrgHal = ', LCMIP6OrgHal, Input_Opt%LCMIP6OrgHal
    END IF

    ! Free pointer
    Spc => NULL()

    IF ( RC/=GC_SUCCESS ) RETURN

  END SUBROUTINE fixSfcVMR
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_co2
!
! !DESCRIPTION: Subroutine INIT\_SfcVMR allocates memory to module arrays and
!  reads in monthly emissions.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE INIT_SfcVMR( am_I_Root, Input_Opt, RC )
!
! !USES:
!
      ! References to F90 modules
      USE CMN_SIZE_MOD
!      USE GIGC_ErrCode_Mod
      USE ErrCode_Mod
      USE Input_Opt_Mod,      ONLY : OptInput
      USE ERROR_MOD,          ONLY : ALLOC_ERR
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
!  13 Mar 2019 - T. Sherwen - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL, SAVE      :: IS_INIT = .FALSE.
      INTEGER            :: AS

      ! For values from Input_Opt
      LOGICAL            :: LHCodedOrgHal
      LOGICAL            :: LCMIP6OrgHal

      !=================================================================
      ! INIT_CO2 begins here!
      !=================================================================

      ! Return success
      RC          = GC_SUCCESS

      ! Exit if we have already initialised
      IF ( IS_INIT ) RETURN

      ! Copy values from Input_Opt
      LHCodedOrgHal     = Input_Opt%LHCodedOrgHal
      LCMIP6OrgHal      = Input_Opt%LCMIP6OrgHal

      ! Check flags at INIT
      WRITE(*,*) 'TMS CHECK INIT_SfcVMR '
      WRITE(*,*) 'TMS CHECK INIT - LCMIP6OrgHal', Input_Opt%LCMIP6OrgHal
      WRITE(*,*) 'TMS CHECK INIT - LHCodedOrgHal', Input_Opt%LHCodedOrgHal
!       Array for Fossil Fuel regions
!       ALLOCATE( FIXED_CH3Br( IIPAR, JJPAR ), STAT=AS )
!       IF ( AS /= 0 ) CALL ALLOC_ERR( 'FIXED_CH3Br' )
!       FIXED_CH3Br = 0e+0_fp
!
!       Array for Biospheric regions
!       ALLOCATE( FIXED_CH3Cl( IIPAR, JJPAR ), STAT=AS )
!       IF ( AS /= 0 ) CALL ALLOC_ERR( 'FIXED_CH3Cl' )
!       FIXED_CH3Cl = 0e+0_fp
!
!       Array for Ocean Regions
!       ALLOCATE( FIXED_CH2Cl2( IIPAR, JJPAR ), STAT=AS )
!       IF ( AS /= 0 ) CALL ALLOC_ERR( 'FIXED_CH2Cl2' )
!       FIXED_CH2Cl2 = 0e+0_fp
!
!       Array for Biospheric regions
!       ALLOCATE( FIXED_CHCl3( IIPAR, JJPAR ), STAT=AS )
!       IF ( AS /= 0 ) CALL ALLOC_ERR( 'FIXED_CHCl3' )
!       FIXED_CHCl3 = 0e+0_fp

      ! Reset IS_INIT flag
      IS_INIT = .TRUE.

      END SUBROUTINE INIT_SfcVMR
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: cleanup_SfcVMR
!
! !DESCRIPTION: Subroutine CLEANUP\_SfcVMR deallocates all module arrays.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE CLEANUP_SfcVMR
!
! !REVISION HISTORY:
!  13 Mar 2019 - T. Sherwen - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
      !=================================================================
      ! CLEANUP_CO2 begins here!
      !=================================================================
!       IF ( ALLOCATED( FIXED_CH3Br  ) ) DEALLOCATE( FIXED_CH3Br )
!       IF ( ALLOCATED( FIXED_CH3Cl  ) ) DEALLOCATE( FIXED_CH3Cl )
!       IF ( ALLOCATED( FIXED_CH2Cl2 ) ) DEALLOCATE( FIXED_CH2Cl2 )
!       IF ( ALLOCATED( FIXED_CHCl3  ) ) DEALLOCATE( FIXED_CHCl3  )

      ! Free pointers
      HCO_CH3BR      => NULL()
      HCO_CH3Cl      => NULL()
      HCO_CH2Cl2     => NULL()
      HCO_CH2Cl2     => NULL()


      END SUBROUTINE CLEANUP_SfcVMR

!EOC
END MODULE SFCVMR_MOD
