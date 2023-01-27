#include "MAPL_Generic.h"
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: gchp_chunk_mod
!
! !DESCRIPTION: Module GC\_CHUNK\_MOD is the module that contains the init,
!  and run methods for the ESMF interface to GEOS-Chem.
!\\
!\\
! !INTERFACE:
!
MODULE GCHP_Chunk_Mod
!
! !USES:
!
  USE MAPL_MOD
  USE ESMF
  USE ErrCode_Mod
  USE Precision_Mod

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC :: GCHP_Chunk_Init
  PUBLIC :: GCHP_Chunk_Run
!
! !PRIVATE MEMBER FUNCTIONS:
!
  INTEGER  ::  MemDebugLevel
!
! !REVISION HISTORY:
!  22 Jun 2009 - R. Yantosca & P. Le Sager - Chunkized & cleaned up.
!  See https://github.com/geoschem/geos-chem for history
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
! !IROUTINE: gchp_chunk_init
!
! !DESCRIPTION: Subroutine GCHP\_CHUNK\_INIT is the ESMF init method for
!  GEOS-Chem.  This routine calls routines within core GEOS-Chem to allocate
!  arrays and read input files.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE GCHP_Chunk_Init( nymdB,         nhmsB,      nymdE,           &
                              nhmsE,         tsChem,     tsDyn,           &
                              tsRad,         lonCtr,     latCtr,          &
#if !defined( MODEL_GEOS )
                              GC,            EXPORT,                      &
#endif
                              Input_Opt,     State_Chm,  State_Diag,      &
                              State_Grid,    State_Met,  HcoConfig,       &
                              HistoryConfig, RC )
!
! !USES:
!
    USE Chemistry_Mod,           ONLY : Init_Chemistry
    USE Emissions_Mod,           ONLY : Emissions_Init
    USE GC_Environment_Mod
    USE GC_Grid_Mod,             ONLY : SetGridFromCtr
    USE GCHP_HistoryExports_Mod, ONLY : HistoryConfigObj
    USE HCO_Types_Mod,           ONLY : ConfigObj
    USE Input_Mod,               ONLY : Read_Input_File
    USE Input_Opt_Mod,           ONLY : OptInput, Set_Input_Opt
    USE Linear_Chem_Mod,         ONLY : Init_Linear_Chem
    USE Linoz_Mod,               ONLY : Linoz_Read
    USE PhysConstants,           ONLY : PI_180
    USE Pressure_Mod,            ONLY : Init_Pressure
    USE Roundoff_Mod,            ONLY : RoundOff
    USE State_Chm_Mod,           ONLY : ChmState, Ind_
    USE State_Diag_Mod,          ONLY : DgnState
    USE State_Grid_Mod,          ONLY : GrdState, Init_State_Grid
    USE State_Met_Mod,           ONLY : MetState
!#if defined( MODEL_GEOS )
!    USE Tendencies_Mod,          ONLY : TEND_INIT
!#endif
    USE Time_Mod,                ONLY : Set_Timesteps
    USE UCX_MOD,                 ONLY : INIT_UCX
    USE UnitConv_Mod,            ONLY : Convert_Spc_Units
#ifdef ADJOINT
    USE Charpak_Mod,             ONLY : To_UpperCase
#endif
#if defined( RRTMG )
    USE RRTMG_RAD_TRANSFER_MOD,  ONLY : Init_RRTMG_Rad_Transfer
    USE RRTMG_LW_Init,           ONLY : RRTMG_LW_Ini
    USE RRTMG_SW_Init,           ONLY : RRTMG_SW_Ini
#endif
!
! !INPUT PARAMETERS:
!
    INTEGER,            INTENT(IN)    :: nymdB       ! YYYYMMDD @ start of run
    INTEGER,            INTENT(IN)    :: nhmsB       ! hhmmss   @ start of run
    INTEGER,            INTENT(IN)    :: nymdE       ! YYYYMMDD @ end of run
    INTEGER,            INTENT(IN)    :: nhmsE       ! hhmmss   @ end of run
    REAL,               INTENT(IN)    :: tsChem      ! Chemistry timestep [s]
    REAL,               INTENT(IN)    :: tsDyn       ! Chemistry timestep [s]
    REAL,               INTENT(IN)    :: tsRad       ! Chemistry timestep [s]
    REAL(ESMF_KIND_R4), INTENT(IN)    :: lonCtr(:,:) ! Lon centers [radians]
    REAL(ESMF_KIND_R4), INTENT(IN)    :: latCtr(:,:) ! Lat centers [radians]
!
! !INPUT/OUTPUT PARAMETERS:
!
#if !defined( MODEL_GEOS )
    TYPE(ESMF_State),    INTENT(INOUT), TARGET :: EXPORT ! Export state object
    TYPE(ESMF_GridComp), INTENT(INOUT)         :: GC     ! Ref to this GridComp
#endif
    TYPE(OptInput),      INTENT(INOUT) :: Input_Opt      ! Input Options object
    TYPE(ChmState),      INTENT(INOUT) :: State_Chm      ! Chem State object
    TYPE(DgnState),      INTENT(INOUT) :: State_Diag     ! Diag State object
    TYPE(GrdState),      INTENT(INOUT) :: State_Grid     ! Grid State object
    TYPE(MetState),      INTENT(INOUT) :: State_Met      ! Met State object
    TYPE(ConfigObj),     POINTER       :: HcoConfig      ! HEMCO config obj
    TYPE(HistoryConfigObj), POINTER    :: HistoryConfig  ! History config obj
!
! !OUTPUT PARAMETERS:
!
    INTEGER,             INTENT(OUT)   :: RC             ! Success or failure?
!
! !REMARKS:
!  Need to add better error checking
!
! !REVISION HISTORY:
!  18 Jul 2011 - M. Long     - Initial Version
!  See https://github.com/geoschem/geos-chem for history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER                        :: I, J, L, STATUS
    CHARACTER(LEN=ESMF_MAXSTR)     :: Iam
    TYPE(ESMF_Config)              :: CF            ! Grid comp config object

#ifdef ADJOINT
    ! Adoint variables
    ! Local Finite Difference variables
    REAL(fp)                       :: FD_LAT, FD_LON
    INTEGER                        :: FD_STEP
    CHARACTER(LEN=ESMF_MAXSTR)     :: FD_SPEC
    REAL(fp)                       :: d, dmin
    INTEGER                        :: imin, jmin, NFD, LFD
    INTEGER                        :: IFD, JFD
    CHARACTER(LEN=ESMF_MAXSTR)     :: FD_TYPE

    ! At present, we are unable to load cube-sphere files through ExtData
    ! so we will define the cost function region thusly in GCHP.rc
    INTEGER                        :: CF_IMIN, CF_IMAX
    INTEGER                        :: CF_JMIN, CF_JMAX
    INTEGER                        :: CF_LMIN, CF_LMAX

    ! Need to get gloabl grid information for some FD spot tests
    TYPE(ESMF_Grid)                :: grid           ! ESMF Grid object
    INTEGER                        :: IL_PET, IU_PET ! Global lon bounds on this PET
    INTEGER                        :: JL_PET, JU_PET ! Global lat bounds on this PET

    ! Model phase: fwd, TLM, ADJOINT
    CHARACTER(LEN=ESMF_MAXSTR)     :: ModelPhase
#endif

    !=======================================================================
    ! GCHP_CHUNK_INIT begins here
    !=======================================================================

    ! Error trap
    Iam = 'GCHP_CHUNK_INIT (gchp_chunk_mod.F90)'

    ! Assume success
    RC = GC_SUCCESS

#if !defined( MODEL_GEOS )
    ! Get memory debug level
    call ESMF_GridCompGet ( GC, config=CF, RC=STATUS )
    _VERIFY(STATUS)
    call ESMF_ConfigGetAttribute(CF, MemDebugLevel, &
                                 Label="MEMORY_DEBUG_LEVEL:" , RC=STATUS)
    _VERIFY(STATUS)
#endif

    ! Update Input_Opt with timing fields
    ! We will skip defining these in READ_INPUT_FILE
    Input_Opt%NYMDb   = nymdB           ! YYYYMMDD @ start of simulation
    Input_Opt%NHMSb   = nhmsB           ! hhmmss   @ end   of simulation
    Input_Opt%NYMDe   = nymdE           ! YYYYMMDD @ start of simulation
    Input_Opt%NHMSe   = nhmsE           ! hhmmss   @ end   of simulation
    Input_Opt%TS_CHEM = INT( tsChem )   ! Chemistry timestep [sec]
    Input_Opt%TS_EMIS = INT( tsChem )   ! Chemistry timestep [sec]
    Input_Opt%TS_DYN  = INT( tsDyn  )   ! Dynamic   timestep [sec]
    Input_Opt%TS_CONV = INT( tsDyn  )   ! Dynamic   timestep [sec]
    Input_Opt%TS_RAD  = INT( tsRad  )   ! RRTMG     timestep [sec]

    ! Read geoschem_config.yml at very beginning of simulation on every CPU
    CALL Read_Input_File( Input_Opt, State_Grid, RC )
    _ASSERT(RC==GC_SUCCESS, 'Error calling Read_Input_File')

    ! Initialize GEOS-Chem horizontal grid structure
    CALL GC_Init_Grid( Input_Opt, State_Grid, RC )
    _ASSERT(RC==GC_SUCCESS, 'Error calling GC_Init_Grid')

    ! Set maximum number of levels in the chemistry grid
    State_Grid%MaxChemLev  = State_Grid%MaxStratLev

    ! In the ESMF/MPI environment, we can get the total overhead ozone
    ! either from the met fields (GCHPsa) or from the Import State (GEOS-5)
    Input_Opt%USE_O3_FROM_MET = .TRUE.

    ! Read LINOZ climatology
    IF ( Input_Opt%LLINOZ ) THEN
       CALL Linoz_Read( Input_Opt, RC )
       _ASSERT(RC==GC_SUCCESS, 'Error calling Linoz_Read')
    ENDIF

    ! Allocate all lat/lon arrays
    CALL GC_Allocate_All( Input_Opt, State_Grid, RC )
    _ASSERT(RC==GC_SUCCESS, 'Error calling GC_Allocate_All')

    ! Set grid based on passed mid-points
    CALL SetGridFromCtr( Input_Opt, State_Grid, lonCtr, latCtr, RC )
    _ASSERT(RC==GC_SUCCESS, 'Error calling "SetGridFromCtr"')

    ! Set GEOS-Chem timesteps on all CPUs
    CALL Set_Timesteps( Input_Opt  = Input_Opt,                              &
                        Chemistry  = Input_Opt%TS_CHEM,                      &
                        Convection = Input_Opt%TS_CONV,                      &
                        Dynamics   = Input_Opt%TS_DYN,                       &
                        Emission   = Input_Opt%TS_EMIS,                      &
                        Radiation  = Input_Opt%TS_RAD,                       &
                        Unit_Conv  = MAX( Input_Opt%TS_DYN,                  &
                                          Input_Opt%TS_CONV ),               &
                        Diagnos    = Input_Opt%TS_DIAG         )

    ! Initialize derived-type objects for met, chem, and diag
    CALL GC_Init_StateObj( HistoryConfig%DiagList,                           &
                           HistoryConfig%TaggedDiagList, Input_Opt,          &
                           State_Chm, State_Diag, State_Grid, State_Met, RC )
    _ASSERT(RC==GC_SUCCESS, 'Error calling GC_Init_StateObj')


#ifdef ADJOINT
    ! Are we running the adjoint?
    call ESMF_ConfigGetAttribute(CF, ModelPhase,            &
                                 Label="MODEL_PHASE:" ,         &
                                 Default="FORWARD",  RC=STATUS)
    _VERIFY(STATUS)
    call WRITE_PARALLEL('Checking if this is adjoint. Model phase = "' // trim(ModelPhase) // '"')
    input_opt%IS_ADJOINT = .FALSE.
    if (TRIM(ModelPhase) .eq. 'ADJOINT') THEN
       call WRITE_PARALLEL('Yes! Setting IS_ADJOINT to true.')
       input_opt%IS_ADJOINT = .TRUE.
    endif

    call ESMF_ConfigGetAttribute(CF, FD_TYPE, &
         Label="FD_TYPE:" , Default='NONE', RC=STATUS)
    _VERIFY(STATUS)

    Input_Opt%IS_FD_GLOBAL = TRIM(To_UpperCase(FD_TYPE(1:4))) == 'GLOB'
    Input_Opt%IS_FD_SPOT   = TRIM(To_UpperCase(FD_TYPE(1:4))) == 'SPOT'
    IF (MAPL_Am_I_Root()) THEN
       WRITE(*,1091) TRIM(FD_TYPE), Input_Opt%IS_FD_GLOBAL, Input_Opt%IS_FD_SPOT
    ENDIF
1091   FORMAT('FD_TYPE = ', a6, ', FD_GLOB = ', L1, ', FD_SPOT = ', L1)



    call ESMF_ConfigGetAttribute(CF, FD_STEP, &
         Label="FD_STEP:" , Default=-1, RC=STATUS)
    _VERIFY(STATUS)

    IF (Input_Opt%IS_FD_GLOBAL .or. Input_Opt%IS_FD_SPOT)  THEN
       _ASSERT(FD_STEP /= -1, 'FD_GLOB or FD_SPOT require FD_STEP')
    ENDIF


    if (.not. FD_STEP == -1 .or. input_opt%IS_ADJOINT) THEN
       Input_Opt%FD_STEP = FD_STEP

       call ESMF_ConfigGetAttribute(CF, FD_SPEC, &
            Label="FD_SPEC:", default="", RC=STATUS)
       _VERIFY(STATUS)
       IF (TRIM(FD_SPEC ) == "") THEN
          NFD = -1
       ELSE
          NFD = Ind_(FD_SPEC)
       ENDIF

       call ESMF_ConfigGetAttribute(CF, IFD, &
            Label="IFD:", default=-1, RC=STATUS)
       _VERIFY(STATUS)

       call ESMF_ConfigGetAttribute(CF, JFD, &
            Label="JFD:", default=-1, RC=STATUS)
       _VERIFY(STATUS)

       ! Get the ESMF grid attached to this gridded component
       CALL ESMF_GridCompGet( GC, grid=Grid, __RC__ )

       ! Get the upper and lower bounds of on each PET using MAPL
       CALL MAPL_GridGetInterior( Grid, IL_PET, IU_PET, JL_PET, JU_PET )

       ! See if we specified IFD and JFD in GCHP.rc
       IF ( IFD > 0 .and. JFD > 0 ) THEN

          if (IL_PET .le. IFD .and. IFD .le. IU_PET .and. &
               JL_PET .le. JFD .and. JFD .le. JU_PET) THEN
             Input_Opt%IS_FD_SPOT_THIS_PET = .true.
             Input_opt%IFD = IFD - IL_PET + 1
             Input_Opt%JFD = JFD - JL_PET + 1

             ! set these for debug printing
             DMIN = 0.0
             IMIN = Input_Opt%IFD
             JMIN = Input_Opt%JFD
          ENDIF

       ELSE

          call ESMF_ConfigGetAttribute(CF, FD_LAT, &
               Label="FD_LAT:", default=-999.0d0, RC=STATUS)
          _VERIFY(STATUS)

          call ESMF_ConfigGetAttribute(CF, FD_LON, &
               Label="FD_LON:", default=-999.0d0, RC=STATUS)
          _VERIFY(STATUS)

          _ASSERT( FD_LAT .ne. -999.0d0 .and. FD_LON .ne. -999.0d0, 'FD_SPOT requires either IFD and JFD or FD_LAT and FD_LON be set in GCHP.rc')


          dmin = 99999.9
          imin = -1
          jmin = -1
          ! try to find lat lon grid cell closest to 44.65, -63.58 (Halifax, NS)
          DO I = 1, state_grid%nx
             DO J = 1, state_grid%ny
                d = sqrt((state_grid%XMID(I,J) - FD_LON)**2 + &
                     (state_grid%YMID(I,J) - FD_LAT)**2)
                if (d < dmin) then
                   dmin = d
                   imin = i
                   jmin = j
                endif
             enddo
          enddo
          ! this is terrible. We need a better way to figure out if we're really in
          ! a grid cell, bbut I don't know how to do that. For now we're just hardcoding
          ! to the value for C24 and hoping for no points near cubed-sphere face
          ! boundaries
          if (dmin < 3.2) then
             ! getting the global grid offset is possible, see Chem_GridCompMod.F90:Extract_
             Input_Opt%IS_FD_SPOT_THIS_PET = .true.
             Input_Opt%IFD = IMIN
             Input_Opt%JFD = JMIN

          end if
       ENDIF

       Input_Opt%NFD = NFD

       call ESMF_ConfigGetAttribute(CF, LFD, &
            Label="LFD:", RC=STATUS)
       _VERIFY(STATUS)

       Input_Opt%LFD = LFD

       ! Read in cost function region

       call ESMF_ConfigGetAttribute(CF, CF_IMIN, &
            Label="CF_IMIN:", default=-1, RC=STATUS)
       _VERIFY(STATUS)

       CF_IMIN = CF_IMIN - IL_PET + 1

       call ESMF_ConfigGetAttribute(CF, CF_IMAX, &
            Label="CF_IMAX:", default=-1, RC=STATUS)
       _VERIFY(STATUS)

       CF_IMAX = CF_IMAX - IL_PET + 1

       call ESMF_ConfigGetAttribute(CF, CF_JMIN, &
            Label="CF_JMIN:", default=-1, RC=STATUS)
       _VERIFY(STATUS)

       CF_JMIN = CF_JMIN - JL_PET + 1

       call ESMF_ConfigGetAttribute(CF, CF_JMAX, &
            Label="CF_JMAX:", default=-1, RC=STATUS)
       _VERIFY(STATUS)

       CF_JMAX = CF_JMAX - JL_PET + 1

       call ESMF_ConfigGetAttribute(CF, CF_LMIN, &
            Label="CF_LMIN:", default=-1, RC=STATUS)
       _VERIFY(STATUS)

       call ESMF_ConfigGetAttribute(CF, CF_LMAX, &
            Label="CF_LMAX:", default=-1, RC=STATUS)
       _VERIFY(STATUS)

       IF (CF_IMIN < 1 .OR. CF_IMIN > State_Grid%NX .OR. &
            CF_IMAX < 1 .OR. CF_IMAX > State_Grid%NX .OR. &
            CF_JMIN < 1 .OR. CF_JMIN > State_Grid%NY .OR. &
            CF_JMAX < 1 .OR. CF_JMAX > State_Grid%NY) THEN
       WRITE(*,1028) Input_Opt%thisCPU,   &
            Input_Opt%CF_IMIN, Input_Opt%CF_IMAX, &
            Input_Opt%CF_JMIN, Input_Opt%CF_JMAX, &
            Input_Opt%CF_LMIN, Input_Opt%CF_LMAX
1028   FORMAT('Pre-CF on Pet ', i3, ' I = (', i3, ', ', i3, ') &
             J = ( ', i3, ', ', i3, ') &
             L = (', i3, ', ', i3, ')')

          CF_IMIN = -1
          CF_IMAX = -1
          CF_JMIN = -1
          CF_JMAX = -1
          CF_LMIN = -1
          CF_LMAX = -1
       ENDIF

       _ASSERT(CF_IMIN * CF_IMAX > 0, 'Please define both max and min for CF_I')
       _ASSERT(CF_JMIN * CF_JMAX > 0, 'Please define both max and min for CF_J')
       _ASSERT(CF_LMIN * CF_LMAX > 0, 'Please define both max and min for CF_L')

       _ASSERT(CF_LMIN * CF_IMIN > 0, 'If CF_I: is defined, please define CF_L')
       _ASSERT(CF_JMIN * CF_IMIN > 0, 'If CF_I: is defined, please define CF_J')
       
       ! At this point, they should all be set or all be negative (probably -1)
       IF (CF_IMIN > 0) THEN
          Input_Opt%CF_IMIN = CF_IMIN
          Input_Opt%CF_IMAX = CF_IMAX
          Input_Opt%CF_JMIN = CF_JMIN
          Input_Opt%CF_JMAX = CF_JMAX
          Input_Opt%CF_LMIN = CF_LMIN
          Input_Opt%CF_LMAX = CF_LMAX
       ELSEIF (Input_Opt%IS_FD_SPOT_THIS_PET) THEN
          Input_Opt%CF_IMIN = Input_Opt%IFD
          Input_Opt%CF_IMAX = Input_Opt%IFD
          Input_Opt%CF_JMIN = Input_Opt%JFD
          Input_Opt%CF_JMAX = Input_Opt%JFD
          Input_Opt%CF_LMIN = Input_Opt%LFD
          Input_Opt%CF_LMAX = Input_Opt%LFD
       WRITE(*,1027) Input_Opt%thisCPU,   &
            Input_Opt%CF_IMIN, Input_Opt%CF_IMAX, &
            Input_Opt%CF_JMIN, Input_Opt%CF_JMAX, &
            Input_Opt%CF_LMIN, Input_Opt%CF_LMAX
1027   FORMAT('CF on Pet ', i3, ' I = (', i3, ', ', i3, ') &
             J = ( ', i3, ', ', i3, ') &
             L = (', i3, ', ', i3, ')')

       ELSE
          Input_Opt%CF_IMIN = -1
          Input_Opt%CF_IMAX = -1
          Input_Opt%CF_JMIN = -1
          Input_Opt%CF_JMAX = -1
          Input_Opt%CF_LMIN = -1
          Input_Opt%CF_LMAX = -1
       ENDIF

       IF ( Input_Opt%IS_FD_SPOT_THIS_PET ) THEN
          write (*,1011) Input_Opt%thisCPU, dmin, imin, jmin, &
               state_grid%YMID(IMIN,JMIN), state_grid%XMID(IMIN,JMIN)

#ifdef DEBUG
          ! Get the ESMF grid attached to this gridded component
          CALL ESMF_GridCompGet( GC, grid=Grid, __RC__ )

          ! Get the upper and lower bounds of on each PET using MAPL
          CALL MAPL_GridGetInterior( Grid, IL_PET, IU_PET, JL_PET, JU_PET )
          WRITE(*,1013) IL_PET, IU_PET
          WRITE(*,1014) JL_PET, JU_PET
#endif
          
         ENDIF

1011   FORMAT('Found FD_SPOT on PET ', i5, ' ', f7.2, &
            ' degrees from cell ', i3, ', ', i3, ' (', f7.2, ', ', f7.2, ')')
1012   FORMAT('Did not find FD_SPOT on PET ', i5, ' ', f7.2,&
            ' degrees from cell ', i3, ', ', i3, ' (', f7.2, ', ', f7.2, ')')
1013   FORMAT('   XminOffset = ', i3, '     XmaxOffset = ', i3)
1014   FORMAT('   YminOffset = ', i3, '     YmaxOffset = ', i3)
1015   FORMAT('   GlobalXMid(', i3, ', ', i3, ') = (', f7.2, ', ' f7.2, ')')
1016   FORMAT('       SPC(', a10, ', FD_SPOT) = ', e22.10)
1019   FORMAT('   SPC_ADJ(', a10, ', FD_SPOT) = ', e22.10)
    ENDIF
#endif

    ! Initialize other GEOS-Chem modules
    CALL GC_Init_Extra( HistoryConfig%DiagList, Input_Opt,    &
                        State_Chm, State_Diag, State_Grid, RC )
    _ASSERT(RC==GC_SUCCESS, 'Error calling GC_Init_Extra')


    ! Set initial species units to internal state units, the same
    ! units as the restart file values. Note that species concentrations
    ! are all still zero at this point since internal state values are not
    ! copied to State_Chm%Species%Conc until Run (post-initialization).
# if defined( MODEL_GEOS )
    State_Chm%Spc_Units = 'kg/kg total'
#else
    State_Chm%Spc_Units = 'v/v dry'
#endif

    ! Initialize chemistry mechanism
    IF ( Input_Opt%ITS_A_FULLCHEM_SIM      .or. &
         Input_Opt%ITS_AN_AEROSOL_SIM      .or. &
         Input_Opt%ITS_A_MERCURY_SIM       .or. &
         Input_Opt%ITS_A_CARBON_SIM      ) THEN
       CALL INIT_CHEMISTRY ( Input_Opt,  State_Chm, State_Diag, &
                             State_Grid, RC )
       _ASSERT(RC==GC_SUCCESS, 'Error calling INIT_CHEMISTRY')
    ENDIF

#if defined( RRTMG )
       ! RRTMG initialization
    IF ( Input_Opt%LRAD ) THEN
       CALL Init_RRTMG_Rad_Transfer( Input_Opt, State_Diag, State_Grid, RC )
       _ASSERT(RC==GC_SUCCESS, 'Error calling "Init_RRTMG_Rad_Transfer"!')
       CALL Rrtmg_Lw_Ini()
       CALL Rrtmg_Sw_Ini()
       State_Chm%RRTMG_iCld  = 0
       State_Chm%RRTMG_iSeed = 10
    ENDIF
#endif

    ! Initialize HEMCO
    CALL EMISSIONS_INIT( Input_Opt, State_Chm, State_Grid, State_Met, RC, &
                         HcoConfig=HcoConfig )
    _ASSERT(RC==GC_SUCCESS, 'Error calling EMISSIONS_INIT')

    ! Initialize UCX routines
    CALL INIT_UCX( Input_Opt, State_Chm, State_Diag, State_Grid )

#if defined( MODEL_GEOS )
    ! Keep commented out line as a GEOS-5 option reminder
    !IF ( Input_Opt%LINEAR_CHEM .AND. Input_Opt%LLSTRAT < value_LM ) THEN
#endif
     IF ( Input_Opt%LINEAR_CHEM ) THEN
       CALL INIT_LINEAR_CHEM( Input_Opt, State_Chm, State_Met, State_Grid, RC )
       _ASSERT(RC==GC_SUCCESS, 'Error calling INIT_LINEAR_CHEM')
    ENDIF

    !-------------------------------------------------------------------------
    ! Diagnostics and tendencies
    !-------------------------------------------------------------------------

!#if defined( MODEL_GEOS )
!    ! The GEOS-Chem diagnostics list, stored in HistoryConfig, is initialized
!    ! during GCHP_INIT_SIMULATION, and corresponding arrays in State_Diag are
!    ! allocated accordingly when initializing State_Diag. Here, we thus
!    ! only need to initialize the tendencies, which have not been initialized
!    ! yet (ckeller, 11/29/17).
!    CALL Tend_Init ( Input_Opt, State_Chm, State_Grid, State_Met, RC )
!    _ASSERT(RC==GC_SUCCESS, 'Error calling Tend_Init')
!#endif

    ! Return success
    RC = GC_Success

  END SUBROUTINE GCHP_Chunk_Init
!EOC

!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: gchp_chunk_run
!
! !DESCRIPTION: Subroutine GCHP\_CHUNK\_RUN is the ESMF run method for
!  GEOS-Chem.
!
! !INTERFACE:
!
  SUBROUTINE GCHP_Chunk_Run( GC,                                             &
                             nymd,       nhms,       year,       month,      &
                             day,        dayOfYr,    hour,       minute,     &
                             second,     utc,        hElapsed,   Input_Opt,  &
                             State_Chm,  State_Diag, State_Grid, State_Met,  &
                             Phase,      IsChemTime, IsRadTime,              &
#if defined( MODEL_GEOS )
                             FrstRewind, &
#endif
#if defined( ADJOINT )
                             IsStarttime, &
#endif
                             RC )
!
! !USES:
!
    ! GEOS-Chem state objects
    USE Input_Opt_Mod,      ONLY : OptInput
    USE State_Chm_Mod,      ONLY : ChmState
    USE State_Diag_Mod
    USE State_Grid_Mod,     ONLY : GrdState
    USE State_Met_Mod,      ONLY : MetState

    ! GEOS-Chem components
    USE Chemistry_Mod,      ONLY : Do_Chemistry, Recompute_OD
    USE Convection_Mod,     ONLY : Do_Convection
    USE DryDep_Mod,         ONLY : Do_DryDep
    USE Emissions_Mod,      ONLY : Emissions_Run
    USE Mixing_Mod,         ONLY : Do_Tend, Do_Mixing
    USE WetScav_Mod,        ONLY : Setup_WetScav, Do_WetDep

    ! HEMCO components (eventually moved to a separate GridComp?)
    USE HCO_State_GC_Mod,   ONLY : HcoState, ExtState
    USE HCO_Interface_Common, ONLY : SetHcoTime
    USE HCO_Interface_GC_Mod, ONLY : Compute_Sflx_For_Vdiff

    ! Specialized subroutines
    USE Calc_Met_Mod,       ONLY : AirQnt
    USE Calc_Met_Mod,       ONLY : Set_Dry_Surface_Pressure
    USE Calc_Met_Mod,       ONLY : Set_Clock_Tracer
    USE Calc_Met_Mod,       ONLY : GCHP_Cap_Tropopause_Prs
    USE Set_Global_CH4_Mod, ONLY : Set_CH4
    USE MODIS_LAI_Mod,      ONLY : Compute_XLAI
    USE PBL_Mix_Mod,        ONLY : Compute_PBL_Height
    USE Pressure_Mod,       ONLY : Set_Floating_Pressures
    USE TOMS_Mod,           ONLY : Compute_Overhead_O3
    USE UCX_Mod,            ONLY : Set_H2O_Trac
    USE Vdiff_Mod,          ONLY : Max_PblHt_for_Vdiff

    ! Utilities
    USE ErrCode_Mod
    USE HCO_Error_Mod
    USE MAPL_MemUtilsMod
    USE Pressure_Mod,       ONLY : Accept_External_Pedge
    USE State_Chm_Mod,      ONLY : IND_
    USE Time_Mod,           ONLY : Accept_External_Date_Time
    USE UnitConv_Mod,       ONLY : Convert_Spc_Units, Print_Global_Species_Kg

    ! Diagnostics
    USE Diagnostics_Mod,    ONLY : Zero_Diagnostics_StartofTimestep
    USE Diagnostics_Mod,    ONLY : Set_Diagnostics_EndofTimestep
#ifdef ADJOINT
    USE PhysConstants,      ONLY : AIRMW
    USE Diagnostics_Mod,    ONLY :  Set_SpcAdj_Diagnostic
#endif
    USE Aerosol_Mod,        ONLY : Set_AerMass_Diagnostic

#if defined( RRTMG )
    USE RRTMG_RAD_TRANSFER_MOD,  ONLY : Do_RRTMG_Rad_Transfer
    USE RRTMG_RAD_TRANSFER_MOD,  ONLY : Set_SpecMask
#endif

    USE Calc_Met_Mod,           ONLY : GET_COSINE_SZA
#if defined( MODEL_GEOS )
    USE HCO_Interface_GC_Mod,   ONLY : HCOI_GC_WriteDiagn
#endif
    USE Species_Mod,   ONLY : Species

!
! !INPUT PARAMETERS:
!
    INTEGER,        INTENT(IN)    :: nymd        ! YYYY/MM/DD @ current time
    INTEGER,        INTENT(IN)    :: nhms        ! hh:mm:ss   @ current time
    INTEGER,        INTENT(IN)    :: year        ! UTC year
    INTEGER,        INTENT(IN)    :: month       ! UTC month
    INTEGER,        INTENT(IN)    :: day         ! UTC day
    INTEGER,        INTENT(IN)    :: dayOfYr     ! UTC day of year
    INTEGER,        INTENT(IN)    :: hour        ! UTC hour
    INTEGER,        INTENT(IN)    :: minute      ! UTC minute
    INTEGER,        INTENT(IN)    :: second      ! UTC second
    REAL*4,         INTENT(IN)    :: utc         ! UTC time [hrs]
    REAL*4,         INTENT(IN)    :: hElapsed    ! Elapsed hours
    INTEGER,        INTENT(IN)    :: Phase       ! Run phase (-1, 1 or 2)
    LOGICAL,        INTENT(IN)    :: IsChemTime  ! Time for chemistry? 
    LOGICAL,        INTENT(IN)    :: IsRadTime   ! Time for RRTMG? 
#if defined( MODEL_GEOS )
    LOGICAL,        INTENT(IN)    :: FrstRewind  ! Is it the first rewind?
#endif
#if defined ( ADJOINT )
    LOGICAL,        INTENT(IN)    :: IsStarttime ! Have we reached the start time
                                                 ! in an adjoint run
#endif 
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ESMF_GridComp), INTENT(INOUT) :: GC          ! Ref to this GridComp
    TYPE(OptInput),      INTENT(INOUT) :: Input_Opt   ! Input Options obj
    TYPE(ChmState),      INTENT(INOUT) :: State_Chm   ! Chemistry State obj
    TYPE(DgnState),      INTENT(INOUT) :: State_Diag  ! Diagnostics State obj
    TYPE(GrdState),      INTENT(INOUT) :: State_Grid  ! Grid State obj
    TYPE(MetState),      INTENT(INOUT) :: State_Met   ! Meteorology State obj
!
! !OUTPUT PARAMETERS:
!
    INTEGER,             INTENT(OUT)   :: RC          ! Return code
!
! !REMARKS:
!
! !REVISION HISTORY:
!  18 Jul 2011 - M. Long     - Initial Version
!  See https://github.com/geoschem/geos-chem for history
!EOP
!------------------------------------------------------------------------------
!BOC
    TYPE(ESMF_STATE)               :: INTSTATE
    TYPE(MAPL_MetaComp), POINTER   :: STATE
    TYPE(ESMF_VM)                  :: VM            ! ESMF VM object
    TYPE(ESMF_Field)               :: IntField
    REAL*8                         :: DT
    CHARACTER(LEN=ESMF_MAXSTR)     :: Iam, OrigUnit
    INTEGER                        :: STATUS, HCO_PHASE, RST
#if defined( MODEL_GEOS )
    INTEGER                        :: I, J, L
#endif

    ! Local logicals to turn on/off individual components
    ! The parts to be executed are based on the input options,
    ! the time step and the phase.
    LOGICAL                        :: DoConv
    LOGICAL                        :: DoDryDep
    LOGICAL                        :: DoEmis
    LOGICAL                        :: DoTend
    LOGICAL                        :: DoTurb
    LOGICAL                        :: DoChem
    LOGICAL                        :: DoWetDep
    LOGICAL                        :: DoRad

    ! First call?
    LOGICAL, SAVE                  :: FIRST    = .TRUE.
    LOGICAL, SAVE                  :: FIRST_RT = .TRUE. ! RRTMG

    ! # of times this routine has been called. Only temporary for printing
    ! processes on the first 10 calls.
    INTEGER, SAVE                  :: NCALLS = 0

    ! Strat. H2O settings
    LOGICAL                        :: SetStratH2O
#if defined( MODEL_GEOS )
    LOGICAL, SAVE                  :: LSETH2O_orig
#endif

    ! For RRTMG
    INTEGER                        :: N

    ! Whether to scale mixing ratio with meteorology update in AirQnt
    LOGICAL, SAVE                  :: scaleMR = .FALSE.

    ! Debug variables
    INTEGER, parameter             :: I_DBG = 6, J_DBG = 5, L_DBG=1
#ifdef ADJOINT
    ! Adjoint Finitie Difference Variables
    INTEGER                        :: IFD, JFD, LFD, NFD
    INTEGER                        :: I, J, L
    REAL*8                         :: CFN
    CHARACTER(len=ESMF_MAXSTR)     :: FD_SPEC, TRACNAME
    TYPE(Species),       POINTER   :: ThisSpc
#endif

    !=======================================================================
    ! GCHP_CHUNK_RUN begins here
    !=======================================================================

    ! Error trap
    Iam = 'GCHP_CHUNK_RUN (gchp_chunk_mod.F90)'

    ! Assume success
    RC = GC_SUCCESS

    ! Get state object (needed for timers)
    CALL MAPL_GetObjectFromGC(GC, STATE, __RC__)

    ! Get the VM for optional memory prints (level >= 2)
    !-----------------------------------
    if ( MemDebugLevel > 0 ) THEN
       call ESMF_VmGetCurrent(VM, RC=STATUS)
       _VERIFY(STATUS)
    endif

    !=======================================================================
    ! Define processes to be covered in this phase
    !
    ! In the standard GEOS-Chem, the following operator sequence is used:
    ! 1. DryDep (kg)
    ! 2. Emissions (kg)
    ! 3. Turbulence (v/v)
    ! 4. Convection (v/v)
    ! 5. Chemistry (kg)
    ! 6. Wetdep (kg)
    !
    ! The GEOS-5 operator sequence is:
    ! 1. Gravity wave drag
    ! 2. Moist (convection)
    ! 3. Chemistry 1 (drydep and emissions)
    ! 4. Surface 1
    ! 5. Turbulence 1
    ! 6. Surface 2
    ! 7. Turbulence 2
    ! 8. Chemistry 2 (chemistry and wet deposition)
    ! 9. Radiation
    !
    ! Here, we use the following operator sequence:
    !
    ! 1.  Convection (v/v) --> Phase 1
    ! 2.  DryDep (kg)      --> Phase 1
    ! 3.  Emissions (kg)   --> Phase 1
    ! 4a. Tendencies (v/v) --> Phase 1
    ! -------------------------------
    ! 4b. Turbulence (v/v) --> Phase 2
    ! 5.  Chemistry (kg)   --> Phase 2
    ! 6.  WetDep (kg)      --> Phase 2
    !
    ! Any of the listed processes is only executed if the corresponding switch
    ! in the geoschem_config.yml file is enabled. If the physics component
    ! already covers convection or turbulence, they should not be applied here!
    ! The tendencies are only applied if turbulence is not done within
    ! GEOS-Chem (ckeller, 10/14/14).
    !
    ! The standard number of phases in GCHP is 1, set in GCHP.rc, which
    ! results in Phase -1 in gchp_chunk_run. This results in executing
    ! all GEOS-Chem components in a single run rather than splitting up
    ! across two runs as is done in GEOS-5. (ewl, 10/26/18)
    !=======================================================================

    ! By default, do processes as defined in geoschem_config.yml. DoTend
    ! defined below.
    DoConv   = Input_Opt%LCONV                    ! dynamic time step
    DoDryDep = Input_Opt%LDRYD .AND. IsChemTime   ! chemistry time step
    DoEmis   = IsChemTime                         ! chemistry time step
#if defined( MODEL_GEOS )
    DoTurb   = Input_Opt%LTURB .AND. IsChemTime   ! dynamic time step
#else
    DoTurb   = Input_Opt%LTURB                    ! dynamic time step
#endif
    DoChem   = Input_Opt%LCHEM .AND. IsChemTime   ! chemistry time step
    DoWetDep = Input_Opt%LWETD                    ! dynamic time step 
    DoRad    = Input_Opt%LRAD  .AND. IsRadTime    ! radiation time step

    ! If Phase is not -1, only do selected processes for given phases:
    ! Phase 1: disable turbulence, chemistry and wet deposition.
    IF ( Phase == 1 ) THEN
       DoTurb   = .FALSE.
       DoChem   = .FALSE.
       DoWetDep = .FALSE.

    ! Phase 2: disable convection, drydep and emissions.
    ELSEIF ( Phase == 2 ) THEN
       DoConv   = .FALSE.
       DoDryDep = .FALSE.
       DoEmis   = .FALSE.
    ENDIF

    ! Check if tendencies need be applied. The drydep and emission calls
    ! only calculates the emission / drydep rates, but do not apply the
    ! tendencies to the tracer array yet. If turbulence is done as part of
    ! GEOS-5, we need to make sure that these tendencies are applied to the
    ! tracer array. If turbulence is explicitly covered by GEOS-Chem,
    ! however, the tendencies become automatically applied within the PBL
    ! mixing routines (DO_MIXING), so we should never apply the tendencies
    ! in this case.
    DoTend = ( DoEmis .OR. DoDryDep ) .AND. .NOT. Input_Opt%LTURB

    ! testing only
    IF ( Input_Opt%AmIRoot .and. NCALLS < 10 ) THEN
       write(*,*) 'GEOS-Chem phase ', Phase, ':'
       write(*,*) 'DoConv   : ', DoConv
       write(*,*) 'DoDryDep : ', DoDryDep
       write(*,*) 'DoEmis   : ', DoEmis
       write(*,*) 'DoTend   : ', DoTend
       write(*,*) 'DoTurb   : ', DoTurb
       write(*,*) 'DoChem   : ', DoChem
       write(*,*) 'DoWetDep : ', DoWetDep
       write(*,*) ' '
    ENDIF

    !-------------------------------------------------------------------------
    ! Pre-Run assignments
    !-------------------------------------------------------------------------

    ! Zero out certain State_Diag arrays. This should not be done in a phase 2
    ! call since this can erase diagnostics filled during phase 1 (e.g., drydep) 
    ! (ckeller, 1/21/2022).
    IF ( Phase /= 2 ) THEN
       CALL Zero_Diagnostics_StartOfTimestep( Input_Opt, State_Diag, RC )
    ENDIF

    ! Pass time values obtained from the ESMF environment to GEOS-Chem
    CALL Accept_External_Date_Time( value_NYMD     = nymd,       &
                                    value_NHMS     = nhms,       &
                                    value_YEAR     = year,       &
                                    value_MONTH    = month,      &
                                    value_DAY      = day,        &
                                    value_DAYOFYR  = dayOfYr,    &
                                    value_HOUR     = hour,       &
                                    value_MINUTE   = minute,     &
                                    value_HELAPSED = hElapsed,   &
                                    value_UTC      = utc,        &
                                    RC             = RC         )

    ! Pass time values obtained from the ESMF environment to HEMCO
#if !defined( MODEL_GEOS )
    CALL SetHcoTime ( HcoState,   ExtState,   year,    month,   day,   &
                      dayOfYr,    hour,       minute,  second,  DoEmis,  RC )
#endif

    ! Calculate MODIS leaf area indexes needed for dry deposition
    CALL Compute_XLAI( Input_Opt, State_Grid, State_Met, RC )

    ! Set the pressure at level edges [hPa] from the ESMF environment
    CALL Accept_External_Pedge( State_Met  = State_Met,   &
                                State_Grid = State_Grid,  &
                                RC         = RC          )

    ! Set dry surface pressure (PS1_DRY) from State_Met%PS1_WET
    CALL SET_DRY_SURFACE_PRESSURE( State_Grid, State_Met, 1 )

    ! Set dry surface pressure (PS2_DRY) from State_Met%PS2_WET
    CALL SET_DRY_SURFACE_PRESSURE( State_Grid, State_Met, 2 )

    ! Initialize surface pressures to match the post-advection pressures
    State_Met%PSC2_WET = State_Met%PS1_WET
    State_Met%PSC2_DRY = State_Met%PS1_DRY
    CALL SET_FLOATING_PRESSURES( State_Grid, State_Met, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    ! Define airmass and related quantities
#if defined( MODEL_GEOS )
    CALL AirQnt( Input_Opt, State_Chm, State_Grid, State_Met, RC, .FALSE. )
#else
    ! Scale mixing ratio with changing met only if FV advection is off.
    ! Only do this the first timestep if DELP_DRY found in restart.
    IF ( FIRST .and. .not. Input_Opt%LTRAN ) THEN
       CALL MAPL_Get ( STATE, INTERNAL_ESMF_STATE=INTSTATE, __RC__ )
       CALL ESMF_StateGet( INTSTATE, 'DELP_DRY', IntField, RC=STATUS )
       _VERIFY(STATUS)
       CALL ESMF_AttributeGet( IntField, NAME="RESTART", VALUE=RST, RC=STATUS )
       _VERIFY(STATUS)
       IF ( .not. ( RST == MAPL_RestartBootstrap .OR. &
                    RST == MAPL_RestartSkipInitial ) ) scaleMR = .TRUE.
       CALL AirQnt( Input_Opt, State_Chm, State_Grid, State_Met, RC, scaleMR )
       scaleMR = .TRUE.
    ELSE
       CALL AirQnt( Input_Opt, State_Chm, State_Grid, State_Met, RC, scaleMR )
    ENDIF
#endif

    ! Initialize/reset wetdep after air quantities computed
    IF ( DoConv .OR. DoChem .OR. DoWetDep ) THEN
       CALL SETUP_WETSCAV( Input_Opt, State_Chm, State_Grid, State_Met, RC )
    ENDIF

    ! Cap the polar tropopause pressures at 200 hPa, in order to avoid
    ! tropospheric chemistry from happening too high up (cf. J. Logan)
    CALL GCHP_Cap_Tropopause_Prs( Input_Opt      = Input_Opt,  &
                                  State_Grid     = State_Grid, &
                                  State_Met      = State_Met,  &
                                  RC             = RC         )

    ! Update clock tracer if relevant
    IF (  IND_('CLOCK','A') > 0 ) THEN
       CALL Set_Clock_Tracer( State_Chm, State_Grid )
    ENDIF

    ! Call PBL quantities. Those are always needed
    CALL Compute_Pbl_Height( Input_Opt, State_Grid, State_Met, RC )
    _ASSERT(RC==GC_SUCCESS, 'Error calling COMPUTE_PBL_HEIGHT')

    ! Convert to dry mixing ratio
    CALL Convert_Spc_Units ( Input_Opt, State_Chm, State_Grid, State_Met, &
                             'kg/kg dry', RC, OrigUnit=OrigUnit )
    _ASSERT(RC==GC_SUCCESS, 'Error calling CONVERT_SPC_UNITS')

    ! SDE 05/28/13: Set H2O to STT if relevant
    IF ( IND_('H2O','A') > 0 ) THEN
       SetStratH2O = .FALSE.
#if defined( MODEL_GEOS )
       !=======================================================================
       ! Tropospheric H2O is always prescribed (using GEOS Q). For strat H2O
       ! there are three options, controlled by toggles 'set initial global MR'
       ! in geoschem_config.yml and 'Prescribe_strat_H2O' in GEOSCHEMchem_GridComp.rc:
       ! (A) never prescribe strat H2O -> both toggles off
       ! (B) prescribe strat H2O on init time step -> toggle in input.goes on
       ! (C) always prescribe strat H2O -> toggle in GEOSCHEMchem_GridComp.rc on
       !=======================================================================
       IF ( FIRST ) THEN
          LSETH2O_orig = Input_Opt%LSETH2O
       ENDIF
       IF ( FIRST .OR. FrstRewind ) THEN
          Input_Opt%LSETH2O = LSETH2O_orig
       ENDIF
       IF ( Input_Opt%LSETH2O .OR. Input_Opt%AlwaysSetH2O ) THEN
#else
       IF ( Input_Opt%LSETH2O ) THEN
#endif
          SetStratH2O = .TRUE.
       ENDIF
       CALL SET_H2O_TRAC( SetStratH2O, Input_Opt, State_Chm, &
                          State_Grid,  State_Met, RC )
       _ASSERT(RC==GC_SUCCESS, 'Error calling SET_H2O_TRAC')

      ! Only force strat once
       IF ( Input_Opt%LSETH2O ) Input_Opt%LSETH2O = .FALSE.
    ENDIF

    ! Compute the cosine of the solar zenith angle array:
    !    State_Met%SUNCOS     => COS(SZA) at the current time
    !    State_Met%SUNCOSmid  => COS(SZA) at the midpt of the chem timestep
    !    COS(SZA) at the midpt of the chem timestep 5hrs ago is now
    !    calculated elsewhere, in the HEMCO PARANOx extension
    CALL GET_COSINE_SZA( Input_Opt, State_Grid, State_Met, RC )
    _ASSERT(RC==GC_SUCCESS, 'Error calling GET_COSINE_SZA')

#ifdef ADJOINT
    if (.not. first) &
         CALL Print_Global_Species_Kg( I_DBG, J_DBG, L_DBG,           &
                                       'CO2', Input_Opt, State_Chm,   &
                                       State_Grid, State_Met, trim(Iam) // &
                                       ' before first unit conversion', RC)
    CALL GCHP_PRINT_MET( I_DBG, J_DBG, L_DBG, Input_Opt,&
         State_Grid, State_Met, trim(Iam) // ' before first unit conversion.', RC)
    
    IF (first .and. Input_Opt%IS_FD_SPOT_THIS_PET .and.  Input_Opt%IS_FD_SPOT) THEN
       FD_SPEC = transfer(state_chm%SpcData(Input_Opt%NFD)%Info%Name, FD_SPEC)
       IFD = Input_Opt%IFD
       JFD = Input_Opt%JFD
       LFD = Input_Opt%LFD
       NFD = Input_Opt%NFD
       WRITE (*, 1017) TRIM(FD_SPEC), state_chm%Species(NFD)%Conc(IFD,JFD,LFD)
       IF (Input_Opt%IS_ADJOINT) THEN
          WRITE(*,*) ' Computing final cost function'
          CFN = 0d0
          state_chm%SpeciesAdj(:,:,:,NFD) = 0d0
          DO L = 1,State_Grid%NZ
          DO J = 1,State_Grid%NY
          DO I = 1,State_Grid%NX
             if (State_chm%CostFuncMask(I,J,L) > 0d0) THEN
                WRITE (*, 1047) I, J, L, State_Chm%Species(NFD)%Conc(I,J,L)
                state_chm%SpeciesAdj(I,J,L, NFD) = 1.0d0
                CFN = CFN + State_Chm%Species(NFD)%Conc(I,J,L)
             endif
          ENDDO
          ENDDO
          ENDDO
          WRITE(*,'(a7, e22.10)') ' CFN = ', CFN
1047      FORMAT('  SPC(', i2, ', ', i2, ', ', i2, ') = ', e22.10)
       ELSE
          IF (Input_Opt%FD_STEP .eq. 0) THEN
             WRITE(*, *) '    Not perturbing'
          ELSEIF (Input_Opt%FD_STEP .eq. 1) THEN
             WRITE(*, *) '    Perturbing +0.1'
             State_Chm%Species(NFD)%Conc(IFD,JFD,LFD) = State_Chm%Species(NFD)%Conc(IFD,JFD,LFD) * 1.1d0
          ELSEIF (Input_Opt%FD_STEP .eq. 2) THEN
             WRITE(*, *) '    Perturbing -0.1'
             State_Chm%Species(NFD)%Conc(IFD,JFD,LFD) = State_Chm%Species(NFD)%Conc(IFD,JFD,LFD) * 0.9d0
          ELSE
             WRITE(*, *) '    FD_STEP = ', Input_Opt%FD_STEP, ' NOT SUPPORTED!'
          ENDIF
          WRITE (*, 1017) TRIM(FD_SPEC), State_Chm%Species(NFD)%Conc(IFD,JFD,LFD)
       ENDIF
    ENDIF

    IF (first .and. Input_Opt%IS_FD_GLOBAL) THEN
       FD_SPEC = transfer(state_chm%SpcData(Input_Opt%NFD)%Info%Name, FD_SPEC)
       NFD = Input_Opt%NFD
       LFD = Input_Opt%LFD
       IF (Input_Opt%IS_FD_SPOT_THIS_PET) THEN
          IFD = Input_Opt%IFD
          JFD = Input_Opt%JFD
          WRITE (*, 1017) TRIM(FD_SPEC), State_Chm%Species(NFD)%Conc(IFD,JFD,LFD)
          IF (Input_Opt%Is_Adjoint) &
               WRITE (*, 1018) TRIM(FD_SPEC), state_chm%SpeciesAdj(IFD, JFD, LFD, NFD)
       ENDIF
       IF (.not. Input_Opt%IS_ADJOINT) THEN
          IF (Input_Opt%FD_STEP .eq. 0) THEN
             WRITE(*, *) '    Not perturbing'
          ELSEIF (Input_Opt%FD_STEP .eq. 1) THEN
             WRITE(*, *) '    Perturbing +0.1'
             State_Chm%Species(NFD)%Conc = State_Chm%Species(NFD)%Conc(:,:,:) * 1.1d0
          ELSEIF (Input_Opt%FD_STEP .eq. 2) THEN
             WRITE(*, *) '    Perturbing -0.1'
             State_Chm%Species(NFD)%Conc = State_Shm%Species(NFD)%Conc(:,:,:) * 0.9d0
          ELSE
             WRITE(*, *) '    FD_STEP = ', Input_Opt%FD_STEP, ' NOT SUPPORTED!'
          ENDIF
          IF (Input_Opt%IS_FD_SPOT_THIS_PET) &
               WRITE (*, 1017) TRIM(FD_SPEC), State_Chm%Species(NFD)%Conc(IFD,JFD,LFD)
       ELSE
          state_chm%SpeciesAdj(:,:,:,:) = 0d0
          IF (NFD > 0) THEN
             IF (LFD > 0) THEN
                IF (Input_opt%amIRoot) THEN
                   WRITE(*,*) ' Setting Level ', LFD, ' forcing to 1'
                ENDIF
                state_chm%SpeciesAdj(:,:,LFD,NFD) = 1d0
             ELSE
                IF (Input_opt%amIRoot) THEN
                   WRITE(*,*) ' Setting all forcing to 1'
                ENDIF
                state_chm%SpeciesAdj(:,:,:,NFD) = 1d0
             ENDIF
          ENDIF
       ENDIF
    ENDIF

1017 FORMAT('       SPC(', a10, ', FD_SPOT) = ', e22.10)
1018   FORMAT('   SPC_ADJ(', a10, ', FD_SPOT) = ', e22.10)
    IF (Input_Opt%IS_FD_SPOT_THIS_PET ) THEN
       FD_SPEC = transfer(state_chm%SpcData(Input_Opt%NFD)%Info%Name, FD_SPEC)
       NFD = Input_Opt%NFD
       IFD = Input_Opt%IFD
       JFD = Input_Opt%JFD
       LFD = Input_Opt%LFD
       WRITE(*,1017) TRIM(FD_SPEC), State_Chm%Species(NFD)%Conc(IFD,JFD,LFD)
       IF (Input_Opt%Is_Adjoint) &
            WRITE (*, 1018) TRIM(FD_SPEC), state_chm%SpeciesAdj(IFD, JFD, LFD, NFD)
    ENDIF
#endif

    !=======================================================================
    ! EMISSIONS. Pass HEMCO Phase 1 which only updates the HEMCO clock
    ! and the HEMCO data list. Should be called every time to make sure
    ! that the HEMCO clock and the HEMCO data list are up to date.
    !=======================================================================
    HCO_PHASE = 1
    CALL EMISSIONS_RUN( Input_Opt, State_Chm, State_Diag, &
                        State_Grid, State_Met, DoEmis, HCO_PHASE, RC  )
    _ASSERT(RC==GC_SUCCESS, 'Error calling EMISSIONS_RUN')

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!                                PHASE 1 or -1                           !!!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !=======================================================================
    ! 1. Convection
    !
    ! Call GEOS-Chem internal convection routines if convection is enabled
    ! in geoschem_config.yml. This should only be done if convection is not
    ! covered by another gridded component and/or the GC species are not made
    ! friendly to this component!!
    !=======================================================================
    IF ( DoConv ) THEN
       if(Input_Opt%AmIRoot.and.NCALLS<10) write(*,*) ' --- Do convection now'
       CALL MAPL_TimerOn( STATE, 'GC_CONV' )

       CALL DO_CONVECTION ( Input_Opt, State_Chm, State_Diag, &
                            State_Grid, State_Met, RC )
       _ASSERT(RC==GC_SUCCESS, 'Error calling DO_CONVECTION')

       CALL MAPL_TimerOff( STATE, 'GC_CONV' )
       if(Input_Opt%AmIRoot.and.NCALLS<10) write(*,*) ' --- Convection done!'
    ENDIF

    !=======================================================================
    ! 2. Dry deposition
    !
    ! Calculates the deposition rates in [s-1].
    !=======================================================================
    IF ( DoDryDep ) THEN
       if(Input_Opt%AmIRoot.and.NCALLS<10) THEN
          write(*,*) ' --- Do drydep now'
          write(*,*) '     Use FULL PBL: ', Input_Opt%PBL_DRYDEP
       endif
       CALL MAPL_TimerOn( STATE, 'GC_DRYDEP' )

       ! Do dry deposition
       CALL Do_DryDep ( Input_Opt, State_Chm, State_Diag, &
                        State_Grid, State_Met, RC )
       _ASSERT(RC==GC_SUCCESS, 'Error calling Do_DryDep')

       CALL MAPL_TimerOff( STATE, 'GC_DRYDEP' )
       if(Input_Opt%AmIRoot.and.NCALLS<10) write(*,*) ' --- Drydep done!'
    ENDIF

    !=======================================================================
    ! 3. Emissions (HEMCO)
    !
    ! HEMCO must be called on first time step to make sure that the HEMCO
    ! data lists are all properly set up.
    !=======================================================================
    IF ( DoEmis ) THEN
#if !defined( MODEL_GEOS )
       ! Optional memory prints (level >= 3)
       if ( MemDebugLevel > 0 ) THEN
          call ESMF_VMBarrier(VM, RC=STATUS)
          _VERIFY(STATUS)
          call MAPL_MemUtilsWrite(VM, &
                  'gchp_chunk_run, before Emissions_Run', RC=STATUS )
          _VERIFY(STATUS)
       endif
#endif

       if(Input_Opt%AmIRoot.and.NCALLS<10) write(*,*) ' --- Do emissions now'
       CALL MAPL_TimerOn( STATE, 'GC_EMIS' )

       ! Do emissions. Pass HEMCO Phase 2 which performs the emissions
       ! calculations.
       HCO_PHASE = 2
       CALL EMISSIONS_RUN( Input_Opt, State_Chm, State_Diag, &
                           State_Grid, State_Met, DoEmis, HCO_PHASE, RC )
       _ASSERT(RC==GC_SUCCESS, 'Error calling EMISSIONS_RUN')

       CALL MAPL_TimerOff( STATE, 'GC_EMIS' )
       if(Input_Opt%AmIRoot.and.NCALLS<10) write(*,*) ' --- Emissions done!'

       ! Optional memory prints (level >= 3)
       if ( MemDebugLevel > 0 ) THEN
          call ESMF_VMBarrier(VM, RC=STATUS)
          _VERIFY(STATUS)
          call MAPL_MemUtilsWrite(VM,&
                  'gchp_chunk_run, after  Emissions_Run', RC=STATUS )
          _VERIFY(STATUS)
       endif

    ENDIF

    !=======================================================================
    ! If physics covers turbulence, simply add the emission and dry
    ! deposition fluxes calculated above to the tracer array, without caring
    ! about the vertical distribution. The tracer tendencies are only added
    ! to the tracers array after emissions, drydep. So we need to use the
    ! emissions time step here.
    !=======================================================================
    IF ( DoTend ) THEN
       if(Input_Opt%AmIRoot.and.NCALLS<10) write(*,*)   &
                           ' --- Add emissions and drydep to tracers'
       CALL MAPL_TimerOn( STATE, 'GC_FLUXES' )

       ! Get emission time step [s].
       _ASSERT(ASSOCIATED(HcoState), 'Error: HcoState not associated')
       DT = HcoState%TS_EMIS

       ! Apply tendencies over entire PBL. Use emission time step.
       CALL DO_TEND( Input_Opt, State_Chm, State_Diag, &
                     State_Grid, State_Met, .FALSE., RC, DT=DT )
       _ASSERT(RC==GC_SUCCESS, 'Error calling DO_TEND')

       ! testing only
       if(Input_Opt%AmIRoot.and.NCALLS<10) write(*,*)   &
                                 '     Tendency time step [s]: ', DT

       CALL MAPL_TimerOff( STATE, 'GC_FLUXES' )
       if(Input_Opt%AmIRoot.and.NCALLS<10) write(*,*)   &
                                 ' --- Fluxes applied to tracers!'
    ENDIF ! Tendencies

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!                              PHASE 2 or -1                             !!!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !=======================================================================
    ! 4. Turbulence
    !
    ! Call GEOS-Chem internal turbulence routines if turbulence is enabled
    ! in geoschem_config.yml. This should only be done if turbulence is not
    ! covered by another gridded component and/or the GC species are not made
    ! friendly to this component!!
    !=======================================================================
    IF ( DoTurb ) THEN
       if(Input_Opt%AmIRoot.and.NCALLS<10) write(*,*) ' --- Do turbulence now'
       CALL MAPL_TimerOn( STATE, 'GC_TURB' )

       ! Only do the following for the non-local PBL mixing
       IF ( Input_Opt%LNLPBL ) THEN

          ! Once the initial met fields have been read in, we need to find
          ! the maximum PBL level for the non-local mixing algorithm.
          ! This only has to be done once. (bmy, 5/28/20)
          IF ( FIRST ) THEN
             CALL Max_PblHt_For_Vdiff( Input_Opt, State_Grid, State_Met, RC )
             _ASSERT(RC==GC_SUCCESS, 'Error calling MAX_PBLHT_FOR_VDIFF')
          ENDIF

          ! Compute the surface flux for the non-local mixing,
          ! (which means getting emissions & drydep from HEMCO)
          ! and store it in State_Chm%Surface_Flux
          CALL Compute_Sflx_For_Vdiff( Input_Opt,  State_Chm, State_Diag,    &
                                       State_Grid, State_Met, RC            )
          _ASSERT(RC==GC_SUCCESS, 'Error calling COMPUTE_SFLX_FOR_VDIFF')
       ENDIF

       ! Do mixing and apply tendencies. This will use the dynamic time step,
       ! which is fine since this call will be executed on every time step.
       CALL DO_MIXING ( Input_Opt, State_Chm, State_Diag,                    &
                        State_Grid, State_Met, RC                           )
       _ASSERT(RC==GC_SUCCESS, 'Error calling DO_MIXING')

       CALL MAPL_TimerOff( STATE, 'GC_TURB' )
       if(Input_Opt%AmIRoot.and.NCALLS<10) write(*,*) ' --- Turbulence done!'
    ENDIF

    ! Set tropospheric CH4 concentrations and fill species array with
    ! current values.
#if defined( MODEL_GEOS )
    IF ( DoTurb .OR. DoTend ) THEN
#else
    IF ( Phase /= 2 .AND. Input_Opt%ITS_A_FULLCHEM_SIM  &
         .AND. IND_('CH4','A') > 0 ) THEN
#endif
       CALL SET_CH4 ( Input_Opt, State_Chm, State_Diag, &
                      State_Grid, State_Met, RC )
       _ASSERT(RC==GC_SUCCESS, 'Error calling SET_CH4')
    ENDIF

    !=======================================================================
    ! 5. Chemistry
    !=======================================================================
    IF ( DoChem ) THEN
       if(Input_Opt%AmIRoot.and.NCALLS<10) write(*,*) ' --- Do chemistry now'
       CALL MAPL_TimerOn( STATE, 'GC_CHEM' )

       IF ( Input_Opt%ITS_A_FULLCHEM_SIM ) THEN
          ! Calculate TOMS O3 overhead. For now, always use it from the
          ! Met field. State_Met%TO3 is imported from PCHEM (ckeller, 10/21/2014).
          CALL COMPUTE_OVERHEAD_O3( Input_Opt, State_Grid, State_Chm, DAY, &
                                    .TRUE., State_Met%TO3, RC )
       ENDIF

#if !defined( MODEL_GEOS )
       ! Set H2O to species value if H2O is advected
       IF ( IND_('H2O','A') > 0 ) THEN
          CALL SET_H2O_TRAC( .FALSE., Input_Opt, &
                             State_Chm, State_Grid, State_Met, RC )
       ENDIF
#endif

       ! Optional memory prints (level >= 3)
       if ( MemDebugLevel > 0 ) THEN
          call ESMF_VMBarrier(VM, RC=STATUS)
          _VERIFY(STATUS)
          call MAPL_MemUtilsWrite(VM, &
                  'gchp_chunk_run:, before Do_Chemistry', RC=STATUS )
          _VERIFY(STATUS)
       endif

       ! Do chemistry
       CALL Do_Chemistry( Input_Opt, State_Chm, State_Diag, &
                          State_Grid, State_Met, RC )
       _ASSERT(RC==GC_SUCCESS, 'Error calling Do_Chemistr')

       CALL MAPL_TimerOff( STATE, 'GC_CHEM' )
       if(Input_Opt%AmIRoot.and.NCALLS<10) write(*,*) ' --- Chemistry done!'

       ! Optional memory prints (level >= 3)
       if ( MemDebugLevel > 0 ) THEN
          call ESMF_VMBarrier(VM, RC=STATUS)
          _VERIFY(STATUS)
          call MAPL_MemUtilsWrite(VM, &
                  'gchp_chunk_run, after  Do_Chemistry', RC=STATUS )
          _VERIFY(STATUS)
       endif

    ENDIF

    !=======================================================================
    ! 6. Wet deposition
    !=======================================================================
    IF ( DoWetDep ) THEN
       if(Input_Opt%AmIRoot.and.NCALLS<10) write(*,*) ' --- Do wetdep now'
       CALL MAPL_TimerOn( STATE, 'GC_WETDEP' )

       ! Do wet deposition
       CALL DO_WETDEP( Input_Opt, State_Chm, State_Diag, &
                       State_Grid, State_Met, RC )
       _ASSERT(RC==GC_SUCCESS, 'Error calling DO_WETDEP')

       CALL MAPL_TimerOff( STATE, 'GC_WETDEP' )
       if(Input_Opt%AmIRoot.and.NCALLS<10) write(*,*) ' --- Wetdep done!'
    ENDIF

    !=======================================================================
    ! Diagnostics
    !=======================================================================

    !==============================================================
    !      ***** U P D A T E  O P T I C A L  D E P T H *****
    !==============================================================
    ! Recalculate the optical depth at the wavelength(s) specified
    ! in the Radiation Menu. This must be done before the call to any
    ! diagnostic and only on a chemistry timestep.
    ! (skim, 02/05/11)
    IF ( DoChem ) THEN
       CALL RECOMPUTE_OD ( Input_Opt, State_Chm, State_Diag, &
                           State_Grid, State_Met, RC )
       _ASSERT(RC==GC_SUCCESS, 'Error calling RECOMPUTE_OD')
    ENDIF

#if defined( RRTMG )
    ! RRTMG diagnostics
    IF ( DoRad ) THEN

       CALL MAPL_TimerOn( STATE, 'GC_RAD' )

       IF ( Input_Opt%amIRoot .AND. FIRST_RT ) THEN
             WRITE( 6, '(a)' ) REPEAT( '#', 79 )
             WRITE( 6, 500 ) 'R R T M G : Radiative Transfer Model (by AER)'
500          FORMAT( '#####', 12x, a, 12x, '#####' )
             WRITE( 6, '(a)' ) REPEAT( '#', 79 )
       ENDIF

       State_Chm%RRTMG_iSeed = State_Chm%RRTMG_iSeed + 15

       !-----------------------------------------------------------
       ! Determine if we are doing clear-sky or all-sky.
       ! Clear-sky is output with all-sky, so we just need
       ! to run once regardless of whether both are required
       ! or just one.
       !-----------------------------------------------------------
       IF (Input_Opt%LSKYRAD(2) ) Then
          State_Chm%RRTMG_iCld = 1
       ELSE
          State_Chm%RRTMG_iCld = 0      !clouds are on
       ENDIF

       !-----------------------------------------------------------
       ! Calculation for each of the potential output types
       ! See: wiki.geos-chem.org/Coupling_GEOS-Chem_with_RRTMG
       !
       ! RRTMG outputs (scheduled in HISTORY.rc):
       !  0-BA  1=O3  2=ME  3=SU   4=NI  5=AM
       !  6=BC  7=OA  8=SS  9=DU  10=PM  11=ST
       !
       ! State_Diag%RadOutInd(1) will ALWAYS correspond to BASE due
       ! to how it is populated from HISTORY.rc diaglist_mod.F90.
       ! BASE is always calculated first since its flux is used to calculate
       ! other RRTMG flux diagnostics.
       !-----------------------------------------------------------

       ! Calculate BASE first
       N = 1

       ! Echo info
       IF ( Input_Opt%amIRoot ) THEN
          PRINT *, 'Calling RRTMG to compute fluxes and optics'
          IF ( FIRST_RT ) THEN
             WRITE( 6, 520 ) State_Diag%RadOutName(N), State_Diag%RadOutInd(N)
          ENDIF
       ENDIF

       ! Generate mask for species in RT
       CALL Set_SpecMask( State_Diag%RadOutInd(N) )

       ! Compute radiative fluxes for the given output
       CALL Do_RRTMG_Rad_Transfer( ThisDay    = Day,                     &
                                   ThisMonth  = Month,                   &
                                   iCld       = State_Chm%RRTMG_iCld,    &
                                   iSpecMenu  = State_Diag%RadOutInd(N), &
                                   iNcDiag    = N,                       &
                                   iSeed      = State_Chm%RRTMG_iSeed,   &
                                   Input_Opt  = Input_Opt,               &
                                   State_Chm  = State_Chm,               &
                                   State_Diag = State_Diag,              &
                                   State_Grid = State_Grid,              &
                                   State_Met  = State_Met,               &
                                   RC         = RC                     )

       ! Trap potential errors
       _ASSERT(RC==GC_SUCCESS, 'Error encounted in Do_RRTMG_Rad_Transfer' )

       ! Calculate for rest of outputs, if any
       DO N = 2, State_Diag%nRadOut
          IF ( Input_Opt%amIRoot .AND. FIRST_RT ) THEN
             WRITE( 6, 520 ) State_Diag%RadOutName(N), State_Diag%RadOutInd(N)
          ENDIF
          CALL Set_SpecMask( State_Diag%RadOutInd(N) )
          CALL Do_RRTMG_Rad_Transfer( ThisDay    = Day,                    &
                                      ThisMonth  = Month,                  &
                                      iCld       = State_Chm%RRTMG_iCld,   &
                                      iSpecMenu  = State_Diag%RadOutInd(N),&
                                      iNcDiag    = N,                      &
                                      iSeed      = State_Chm%RRTMG_iSeed,  &
                                      Input_Opt  = Input_Opt,              &
                                      State_Chm  = State_Chm,              &
                                      State_Diag = State_Diag,             &
                                      State_Grid = State_Grid,             &
                                      State_Met  = State_Met,              &
                                      RC         = RC          )
          _ASSERT(RC==GC_SUCCESS, 'Error encounted in Do_RRTMG_Rad_Transfer')
       ENDDO

520    FORMAT( 5x, '- ', &
                  a4, ' (Index = ', i2.2, ')' )

       IF ( FIRST_RT ) THEN
          FIRST_RT = .FALSE.
       ENDIF

       CALL MAPL_TimerOff( STATE, 'GC_RAD' )
    ENDIF
#endif

    if(Input_Opt%AmIRoot.and.NCALLS<10) write(*,*) ' --- Do diagnostics now'
    CALL MAPL_TimerOn( STATE, 'GC_DIAGN' )

    ! Set certain diagnostics dependent on state at end of step. This
    ! includes species concentration and dry deposition flux.
    ! For GEOS, this is now done in Chem_GridCompMod.F90. This makes sure
    ! that the diagnostics include any post-run updates (e.g., if assimilation
    ! increments are being applied (ckeller, 2/7/22).
#if !defined( MODEL_GEOS )
    CALL Set_Diagnostics_EndofTimestep( Input_Opt,  State_Chm, State_Diag, &
                                        State_Grid, State_Met, RC )
    _ASSERT(RC==GC_SUCCESS, 'Error calling Set_Diagnostics_EndofTimestep')
#endif

    ! Archive aerosol mass and PM2.5 diagnostics
    IF ( State_Diag%Archive_AerMass ) THEN
       CALL Set_AerMass_Diagnostic( Input_Opt,  State_Chm, State_Diag, &
                                    State_Grid, State_Met, RC )
       _ASSERT(RC==GC_SUCCESS, 'Error calling Set_AerMass_Diagnostic')
    ENDIF

#if defined( MODEL_GEOS )
    ! Save specific humidity and dry air mass for total mixing ratio
    ! adjustment in next timestep, if needed (ewl, 11/8/18)
    State_Met%SPHU_PREV = State_Met%SPHU
#endif

#ifdef ADJOINT
       if (Input_Opt%IS_FD_SPOT_THIS_PET .and. Input_opt%IFD > 0) THEN
       DO N = 1, State_Chm%nSpecies
          ThisSpc => State_Chm%SpcData(N)%Info
          write(*,*) 'SpcAdj(', TRIM(thisSpc%Name), ') = ',  &
               State_Chm%SpeciesAdj(Input_Opt%IFD,Input_Opt%JFD,Input_Opt%LFD,N)
       ENDDO
       ENDIF
    !=======================================================================
    ! If this is an adjoint run, we need to check for the final (first)
    ! timestep and multiply the scaling factor adjoint by the initial concs
    !=======================================================================
    IF (Input_Opt%IS_ADJOINT .and. IsStarttime) THEN
       if (Input_opt%amIRoot) WRITE(*,*) '   Adjoint multiplying SF_ADJ by ICS'
       DO N = 1, State_Chm%nSpecies
          ThisSpc => State_Chm%SpcData(N)%Info

          ! Find the non-adjoint variable or this
          TRACNAME = ThisSpc%Name

          State_Chm%SpeciesAdj(:,:,:,N) = State_Chm%SpeciesAdj(:,:,:,N) * State_Chm%Species(N)%Conc(:,:,:) * &
               ( AIRMW / State_Chm%SpcData(N)%Info%MW_g )

          if (Input_Opt%IS_FD_SPOT_THIS_PET .and. Input_Opt%IFD > 0) THEN
             write(*,*) 'After conversion ',  &
                  State_Chm%SpeciesAdj(Input_Opt%IFD,Input_Opt%JFD,Input_Opt%LFD,N)
          ENDIF
       ENDDO

       CALL Set_SpcAdj_Diagnostic( Input_Opt,  State_Chm, State_Diag,        &
                                   State_Grid, State_Met, RC                )
    ENDIF
#endif

    CALL MAPL_TimerOff( STATE, 'GC_DIAGN' )
    if(Input_Opt%AmIRoot.and.NCALLS<10) write(*,*) ' --- Diagnostics done!'

    !=======================================================================
    ! Convert State_Chm%Species units
    !=======================================================================
    CALL Convert_Spc_Units ( Input_Opt, State_Chm, State_Grid, State_Met, &
                             OrigUnit, RC )
    _ASSERT(RC==GC_SUCCESS, 'Error calling CONVERT_SPC_UNITS')

    !=======================================================================
    ! Clean up
    !=======================================================================

    ! testing only
    IF ( PHASE /= 1 .AND. NCALLS < 10 ) NCALLS = NCALLS + 1

    ! First call is done
    FIRST = .FALSE.

    ! Return success
    RC = GC_SUCCESS

  END SUBROUTINE GCHP_Chunk_Run
!EOC

!BOP
  SUBROUTINE GCHP_PRINT_MET(I, J, L,         &
       Input_Opt, State_Grid, State_Met, LOC, RC )

    !
    ! !USES:
    !
    USE State_Met_Mod,        ONLY : MetState
    USE Input_Opt_Mod,        ONLY : OptInput
    USE State_Grid_Mod,       ONLY : GrdState

    !
    ! !INPUT PARAMETERS: 
    !
    INTEGER,          INTENT(IN)    :: I         ! Grid cell lat index
    INTEGER,          INTENT(IN)    :: J         ! Grid cell lon index
    INTEGER,          INTENT(IN)    :: L         ! Grid cell lev index
    CHARACTER(LEN=*), INTENT(IN)    :: LOC       ! Call location string
    TYPE(OptInput),   INTENT(IN)    :: Input_Opt ! Input Options object
    TYPE(GrdState),   INTENT(IN)    :: State_Grid! Grid State object
    TYPE(MetState),   INTENT(IN)    :: State_Met ! Meteorology State object
    !
    ! !INPUT/OUTPUT PARAMETERS: 
    !

    !
    ! !OUTPUT PARAMETERS:
    !
    INTEGER,          INTENT(OUT)   :: RC        ! Success or failure?! 
    ! !REMARKS:
    !
    ! !REVISION HISTORY: 
    !EOP
    !------------------------------------------------------------------------------
    !BOC
    !
    ! !LOCAL VARIABLES:
    !     
    CHARACTER(LEN=255) :: ErrorMsg, ThisLoc


    !=========================================================================
    ! GCHP_PRINT_MET begins here!
    !=========================================================================

    ErrorMsg  = ''
    ThisLoc   = ' -> at GCHP_Print_Met (in module ' // &
         'Interfaces/GCHP/gchp_chunk_mod.F)'

    ! Assume success
    RC = GC_SUCCESS

    ! Echo info
    IF ( Input_Opt%amIRoot ) THEN
       WRITE( 6, 100 ) TRIM( LOC )
       WRITE( 6, 113 ) State_Grid%YMid(I,J), State_Grid%XMid(I,J)
    ENDIF
100 FORMAT( /, '%%%%% GCHP_PRINT_MET at ', a )
113 FORMAT( 'Lat: ', f5.1, '   Lon: ', f5.1 )

    ! Write formatted output
    IF ( Input_Opt%amIRoot ) THEN
       ! 2-D Fields
       WRITE( 6, 114 ) 'PBLH',     State_Met%PBLH(I,J),     I, J
       WRITE( 6, 114 ) 'PSC2_WET', State_Met%PSC2_WET(I,J), I, J
       WRITE( 6, 114 ) 'PSC2_DRY', State_Met%PSC2_DRY(I,J), I, J
       WRITE( 6, 114 ) 'PS1_WET',  State_Met%PS1_WET(I,J), I, J
       WRITE( 6, 114 ) 'PS1_DRY',  State_Met%PS1_DRY(I,J), I, J
       WRITE( 6, 114 ) 'PS2_WET',  State_Met%PS2_WET(I,J), I, J
       WRITE( 6, 114 ) 'PS2_DRY',  State_Met%PS2_DRY(I,J), I, J
       WRITE( 6, 114 ) 'TS',       State_Met%TS(I,J),       I, J
       WRITE( 6, 114 ) 'U10M',     State_Met%U10M(I,J),     I, J
       ! 3-D Fields
       WRITE( 6, 115 ) 'CLDF',     State_Met%CLDF(I,J,L),      I, J, L
       WRITE( 6, 115 ) 'OMEGA',    State_Met%OMEGA(I,J,L),     I, J, L
       WRITE( 6, 115 ) 'PEDGE',    State_Met%PEDGE(I,J,L),     I, J, L
       WRITE( 6, 115 ) 'T',        State_Met%T(I,J,L),         I, J, L
       WRITE( 6, 115 ) 'U',        State_Met%U(I,J,L),         I, J, L
       WRITE( 6, 115 ) 'V',        State_Met%V(I,J,L),         I, J, L
       WRITE( 6, 115 ) 'AD',       State_Met%AD(I,J,L),        I, J, L
       WRITE( 6, 115 ) 'PREVSPHU', State_Met%SPHU_PREV(I,J,L), I, J, L
       WRITE( 6, 115 ) 'SPHU',     State_Met%SPHU(I,J,L),      I, J, L
       ! terminator
       WRITE( 6, 120 )
    ENDIF
114 FORMAT( 'Grid cell  for ', a8, ' = ', es24.16, ', I,J  = ',2I4 )
115 FORMAT( 'Grid cell  for ', a8, ' = ', es24.16, ', I,J,L= ',3I4 )
120 FORMAT( / )


  END SUBROUTINE GCHP_PRINT_MET
!EOC
END MODULE GCHP_Chunk_Mod
