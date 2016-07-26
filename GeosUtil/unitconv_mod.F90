!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: unitconv_mod.F90
!
! !DESCRIPTION: Module UNITCONV\_MOD contains routines which are used to 
!  convert the units of tracer concentrations between mass mixing ratio 
!  [kg/kg air], mass per grid box per area [kg/m2], molar mixing ratio 
!  [vol/vol], and molecular number density [molecules/cm3]. There are 
!  different conversion routines for dry air and total (wet) air mixing 
!  ratios. Conversions involving column area will be phased out for
!  grid-independent GEOS-Chem. 
!\\  
!\\
! !INTERFACE: 
!
MODULE UnitConv_Mod
!
! !USES:
!
  ! GEOS-Chem Modules
  USE CMN_SIZE_MOD          ! Size parameters
  USE PRECISION_MOD         ! GEOS-Chem Flexible Precision (fp)
  USE GIGC_ErrCode_Mod
  USE ERROR_MOD
                    
  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! State_Chm%TRACERS: KG/KG <-> V/V
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  ! kg/kg dry air <-> v/v dry air
  ! Used in DO_TEND in mixing
  PUBLIC  :: Convert_KgKgDry_to_VVDry
  PUBLIC  :: Convert_VVDry_to_KgKgDry

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! State_Chm%TRACERS: MIXING RATIO <-> AREA DENSITY
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  ! kg/kg dry air <-> kg/m2
  ! Used for wet deposition, DO_TEND in mixing,
  ! and around AIRQNT and SET_H2O_TRAC in main
  PUBLIC  :: Convert_KgKgDry_to_Kgm2
  PUBLIC  :: Convert_kgm2_to_KgKgDry

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! AREA-DEPENDENT (temporary routines)
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  ! v/v dry air <-> kg/grid box (State_Chm%TRACERS)
  ! Temporarily replaces legacy CONVERT_UNITS
  ! Used in strat_chem_mod and sulfate_mod
  PUBLIC  :: Convert_VVDry_to_Kg
  PUBLIC  :: Convert_Kg_to_VVDry

  ! kg/kg dry air <-> kg/grid box (State_Chm%TRACERS)
  ! Used in aerosol_mod, tomas_mod, emissions_mod,
  ! strat_chem_mod, exchange_mod, rrtmg_rad_transfer_mod,
  ! chemistry_mod, sulfate_mod, and carbon_mod
  ! This is since RRTMG, TOMAS, exchange_mod, chemistry,
  ! and EMISSMERCURY are still in [kg]
  PUBLIC  :: Convert_KgKgDry_to_Kg
  PUBLIC  :: Convert_Kg_to_KgKgDry

  ! kg/kg dry air <-> kg/grid box (State_Chm%TRACERS, single box only)
  ! Used for TOMAS compatibility in WASHOUT
  PUBLIC  :: ConvertBox_KgKgDry_to_Kg
  PUBLIC  :: ConvertBox_Kg_to_KgKgDry

  ! kg <-> kg/m2 (State_Chm%TRACERS, single box only)
  ! Used for TOMAS compatibility in WASHOUT within wetscav_mod
  PUBLIC  :: ConvertBox_Kg_to_Kgm2
  PUBLIC  :: ConvertBox_Kgm2_to_Kg 
!
! !REMARKS:
!  The routines in this module are used to convert the units of tracer 
!  concentrations in various GEOS-Chem routines.
!
! !REVISION HISTORY:
!  23 Jun 2015 - E. Lundgren - Initial version
!  13 Aug 2015 - E. Lundgren - Add tracer unit error handling
!  29 Sep 2015 - E. Lundgren - Adjust some of the unit conversions to/from kg
!                              to be for a single grid box for TOMAS
!  26 Jul 2016 - E. Lundgren - Remove unused conversions and use "Box" in 
!                              TOMAS-specific unit conversions
!EOP
!------------------------------------------------------------------------------
!BOC
CONTAINS
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: convert_kgkgdry_to_vvdry
!
! !DESCRIPTION: Subroutine Convert\_KgKgDry\_to\_VVDry converts the units of 
!  tracer concentrations from mass mixing ratio (KGKG) [kg/kg] to 
!  volume ratio (VR) [vol/vol] (same as molar ratio [mol/mol]). 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Convert_KgKgDry_to_VVDry( am_I_Root, Input_Opt, &
                                       State_Chm, RC ) 
!
! USES: 
!
    USE GIGC_Input_Opt_Mod, ONLY : OptInput
    USE GIGC_State_Chm_Mod, ONLY : ChmState
!
! !INPUT PARAMETERS: 
!
    LOGICAL,        INTENT(IN)    :: am_I_Root   ! Are we on the root CPU?
    TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
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
!
! !REVISION HISTORY: 
!  08 Jan 2015 - E. Lundgren - Initial version
!  30 Sep 2015 - E. Lundgren - Remove N_TRACERS from arguments list
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER :: I, J, L, N
    CHARACTER(LEN=255) :: MSG, LOC
    
    !====================================================================
    ! Convert_KgKgDry_to_VVDry begins here!
    !====================================================================

    ! Assume success
    RC        =  GIGC_SUCCESS

    ! Verify correct initial units. If current units are unexpected,
    ! write error message and location to log, then pass failed RC
    ! to calling routine. 
    IF ( TRIM( State_Chm%Trac_Units ) /= 'kg/kg dry' ) THEN
       MSG = 'Incorrect initial units: ' // TRIM( State_Chm%Trac_Units )
       LOC = 'Routine Convert_KgKgDry_to_VVDry in unitconv_mod.F90'
       CALL GIGC_Error( MSG, RC, LOC )
       RETURN
    ENDIF

    !====================================================================
    !
    !  The conversion is as follows:
    !
    !   kg tracer(N)    g dry air      mol tracer(N)    
    !   -----------  * ----------  *  -------------  
    !     kg air       mol air         g tracer(N)          
    !
    !   = mass mixing ratio * ratio of air to tracer molecular weights  
    !   
    !   = molar ratio
    !
    ! Therefore, with:
    !
    !  TCVV(N) = dry air molecular wt / tracer molecular wt 
    !     
    ! the conversion is:
    ! 
    !  Tracers(I,J,L,N) [vol/vol]
    !
    !    = Tracers(I,J,L,N) [kg/kg] * TCVV(N)
    !                   
    !====================================================================
 
    !$OMP PARALLEL DO           &
    !$OMP DEFAULT( SHARED     ) &
    !$OMP PRIVATE( I, J, L, N ) 
    DO N = 1, Input_Opt%N_TRACERS
    DO L = 1, LLPAR
    DO J = 1, JJPAR
    DO I = 1, IIPAR
       State_Chm%Tracers(I,J,L,N) = State_Chm%Tracers(I,J,L,N)  &
                                  * Input_Opt%TCVV(N)
    ENDDO
    ENDDO
    ENDDO
    ENDDO
    !$OMP END PARALLEL DO

    ! Update tracer units
    State_Chm%Trac_Units = 'v/v dry'

  END SUBROUTINE Convert_KgKgDry_to_VVDry
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: convert_vvdry_to_kgkgdry
!
! !DESCRIPTION: Subroutine Convert\_VVDry\_to\_KgKgDry converts the units of 
!  tracer concentrations from volume ratio (VR) [vol/vol] (same 
!  as molar mixing ratio [mol/mol]) to mass mixing ratio [kg/kg]. 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Convert_VVDry_to_KgKgDry( am_I_Root, Input_Opt, &
                                       State_Chm, RC ) 
!
! USES: 
!
    USE GIGC_Input_Opt_Mod, ONLY : OptInput
    USE GIGC_State_Chm_Mod, ONLY : ChmState
!
! !INPUT PARAMETERS: 
!
    LOGICAL,        INTENT(IN)    :: am_I_Root   ! Are we on the root CPU?
    TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
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
!
! !REVISION HISTORY: 
!  08 Jan 2015 - E. Lundgren - Initial version
!  30 Sep 2015 - E. Lundgren - Remove N_TRACERS from arguments list
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER :: I, J, L, N
    CHARACTER(LEN=255) :: MSG, LOC

    !====================================================================
    ! Convert_VVDry_to_KgKgDry begins here!
    !=================================================================

      ! Assume success
      RC        =  GIGC_SUCCESS

    ! Verify correct initial units. If current units are unexpected,
    ! write error message and location to log, then pass failed RC
    ! to calling routine. 
    IF ( TRIM( State_Chm%Trac_Units ) /= 'v/v dry' ) THEN
       MSG = 'Incorrect initial units: ' // TRIM( State_Chm%Trac_Units )
       LOC = 'Routine Convert_VVDry_to_KgKgDry in unitconv_mod.F90'
       CALL GIGC_Error( MSG, RC, LOC )
       RETURN
    ENDIF

         !==============================================================
         !
         !  The conversion is as follows:
         !
         !   mol tracer(N)  mol dry air     g tracer(N)         
         !   -----------  * -----------  *  -------------  
         !     mol air       g dry air      mol tracer(N)           
         !
         !   = volume ratio / ratio of air to tracer molecular wts  
         !   
         !   = mass mixing ratio ([g/g] is equivalent to [kg/kg])
         !
         ! Therefore, with:
         !
         !  TCVV(N) = dry air molecular wt / tracer molecular wt 
         !     
         ! the conversion is:
         ! 
         !  Tracers(I,J,L,N) [vol/vol]
         !
         !    = Tracers(I,J,L,N) [kg/kg] / TCVV(N)
         !                   
         !==============================================================

      !$OMP PARALLEL DO           &
      !$OMP DEFAULT( SHARED     ) &
      !$OMP PRIVATE( I, J, L, N ) 
      DO N = 1, Input_Opt%N_TRACERS
      DO L = 1, LLPAR
      DO J = 1, JJPAR
      DO I = 1, IIPAR
        State_Chm%Tracers(I,J,L,N) = State_Chm%Tracers(I,J,L,N)   &
                                        / Input_Opt%TCVV(N)
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      !$OMP END PARALLEL DO

    ! Update tracer units
    State_Chm%Trac_Units = 'kg/kg dry'

    END SUBROUTINE Convert_VVDry_to_KgKgDry
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: convert_kgkgdry_to_kgm2
!
! !DESCRIPTION: Subroutine Convert\_kgkgdry\_to\_kgm2 converts the units of 
!  a 3D array from dry mass mixing ratio [kg/kg dry air] to area density 
!  [kg/m2].  
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Convert_KgKgDry_to_Kgm2( am_I_Root, Input_Opt, State_Met,  &
                                      State_Chm, RC             )
!
! !USES:
!
    USE GIGC_Input_Opt_Mod, ONLY : OptInput
    USE GIGC_State_Chm_Mod, ONLY : ChmState
    USE GIGC_State_Met_Mod, ONLY : MetState
    USE PHYSCONSTANTS

!
! !INPUT PARAMETERS: 
!
    LOGICAL,        INTENT(IN)    :: am_I_Root     ! Are we on the root CPU?
    TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
    TYPE(MetState), INTENT(IN)    :: State_Met     ! Meteorology state object

!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry state object 
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC            ! Success or failure?
!
! !REMARKS:
!
! !REVISION HISTORY: 
!  14 Aug 2015 - E. Lundgren - Initial version
!  30 Sep 2015 - E. Lundgren - Remove N_TRACERS from arguments list
!  06 Jul 2016 - E. Lundgren - Replace DELP and SPHU with DELP_DRY
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER            :: I, J, L, N
    CHARACTER(LEN=255) :: MSG, LOC
!    REAL(fp)           :: SPHU_kgkg
!    REAL(fp)           :: temp

    !====================================================================
    ! Convert_KgKgDry_to_Kgm2 begins here!
    !====================================================================

    ! Assume success
    RC        =  GIGC_SUCCESS

    ! Verify correct initial units. If current units are unexpected,
    ! write error message and location to log, then pass failed RC
    ! to calling routine. 
    IF ( TRIM( State_Chm%Trac_Units ) /= 'kg/kg dry' ) THEN
       MSG = 'Incorrect initial units: ' // TRIM( State_Chm%Trac_Units )
       LOC = 'Routine Convert_KgKgDry_to_Kgm2 in unitconv_mod.F90'
       CALL GIGC_Error( MSG, RC, LOC )
       RETURN
    ENDIF

    !====================================================================
    !
    !  The conversion is as follows:
    !
! MOISTURE FIX
!    !   kg tracer      Delta P [hPa]   100 [Pa]     kg dry air 
!    !   -----------  * ------------- * -------- * --------------     
!    !   kg dry air     g [m/s2]        [hPa]      kg total air  
    !   kg tracer      Delta dry P [hPa]   100 [Pa]
    !   -----------  * ----------------- * --------
    !   kg dry air     g [m/s2]            [hPa]   
    !
    !   = kg / m2
    !
    ! where:
    !
    !  Delta dry P = dry pressure difference across level as derived
    !                from the dry surface pressure with A and B params
    !  g = acceleration due to gravity
    !  kg dry air / kg total air  = 1 - specific humidity
    !     
    !====================================================================

    !$OMP PARALLEL DO        &
    !$OMP DEFAULT( SHARED  ) &
    !$OMP PRIVATE( I, J, L ) 
    DO N = 1, Input_Opt%N_TRACERS
    DO L = 1, LLPAR
    DO J = 1, JJPAR
    DO I = 1, IIPAR

! Prior to 7/6/2016:
!       ! Convert specific humidity from [g/kg] to [kg/kg]
!       SPHU_kgkg = State_Met%SPHU(I,J,L) * 1.0e-3_fp 
!
!       State_Chm%TRACERS(I,J,L,N) = State_Chm%TRACERS(I,J,L,N)          &
!                                    * ( ( 1.0e+0_fp - SPHU_kgkg )       &
!                                    * ( g0_100                          &
!                                    * State_Met%DELP(I,J,L) ) )           
! MOISTURE FIX:
       ! MOISTURE FIX - use new State_Met variable DELP_DRY from sfc P
       ! Area-independent conversion
       State_Chm%TRACERS(I,J,L,N) = State_Chm%TRACERS(I,J,L,N)          &
                                    * ( g0_100                          &
                                    * State_Met%DELP_DRY(I,J,L) )           

    ENDDO
    ENDDO
    ENDDO
    ENDDO
    !$OMP END PARALLEL DO

    ! Update tracer units
    State_Chm%Trac_Units = 'kg/m2'

  END SUBROUTINE Convert_KgKgDry_to_Kgm2
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: convert_kgm2_to_kgkgdry
!
! !DESCRIPTION: Subroutine Convert\_Kgm2\_to\_kgkgdry converts the units of 
!  tracer concentrations from area density [kg/m2] to dry mass mixing ratio 
!  [kg/kg dry air].  
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Convert_Kgm2_to_KgKgDry( am_I_Root, Input_Opt, State_Met, &
                                      State_Chm, RC          )
!
! !USES:
!
    USE GIGC_Input_Opt_Mod, ONLY : OptInput
    USE GIGC_State_Met_Mod, ONLY : MetState
    USE GIGC_State_Chm_Mod, ONLY : ChmState
    USE PHYSCONSTANTS
!
! !INPUT PARAMETERS: 
!
    LOGICAL,        INTENT(IN)    :: am_I_Root   ! Are we on the root CPU?
    TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
    TYPE(MetState), INTENT(IN)    :: State_Met   ! Meteorology state object
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry state object 
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
!
! !REMARKS:
!
! !REVISION HISTORY: 
!  14 Aug 2015 - E. Lundgren - Initial version
!  30 Sep 2015 - E. Lundgren - Replace N_TRACERS with Input_Opt in args list
!  06 Jul 2016 - E. Lundgren - Replace DELP and SPHU with DELP_DRY
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER            :: I, J, L, N
    CHARACTER(LEN=255) :: MSG, LOC
!    REAL(fp)           :: SPHU_kgkg
!    REAL(fp)           :: temp

    !====================================================================
    ! Convert_Kgm2_to_KgKgDry begins here!
    !====================================================================

    ! Assume success
    RC        =  GIGC_SUCCESS

    ! Verify correct initial units. If current units are unexpected,
    ! write error message and location to log, then pass failed RC
    ! to calling routine. 
    IF ( TRIM( State_Chm%Trac_Units ) /= 'kg/m2' ) THEN
       MSG = 'Incorrect initial units: ' // TRIM( State_Chm%Trac_Units )
       LOC = 'Routine Convert_Kgm2_to_KgKgDry in unitconv_mod.F90'
       CALL GIGC_Error( MSG, RC, LOC )
       RETURN
    ENDIF

    !====================================================================
    !
    !  The conversion is as follows:
    !
!    !   kg tracer(N)   g [m/s2]        [hPa]      kg total air   
!    !   -----------  * ------------- * -------- * --------------     
!    !        m2        Delta P [hPa]   100 [Pa]    kg dry air
    !   kg tracer(N)   g [m/s2]            [hPa]
    !   -----------  * ----------------- * --------
    !        m2        Delta dry P [hPa]   100 [Pa]
    !
    !   = kg tracer(N) / kg dry air
    !
    ! where:
    !
    !  Delta dry P = dry pressure difference across level as derived
    !                from the dry surface pressure with A and B params
    !  g = acceleration due to gravity
    !  kg dry air / kg total air  = 1 - specific humidity
    !     
    !====================================================================

    !$OMP PARALLEL DO           &
    !$OMP DEFAULT( SHARED     ) &
    !$OMP PRIVATE( I, J, L, N ) 
    DO N = 1, Input_Opt%N_TRACERS
    DO L = 1, LLPAR
    DO J = 1, JJPAR
    DO I = 1, IIPAR

! Prior to 7/6/2016:
!       ! Convert specific humidity from [g/kg] to [kg/kg]
!       SPHU_kgkg = State_Met%SPHU(I,J,L) * 1.0e-3_fp 
!
!       ! Area-independent conversion
!       State_Chm%TRACERS(I,J,L,N) = State_Chm%TRACERS(I,J,L,N)          &
!                                    * ( 1.0e+0_fp                       &      
!                                    / ( ( 1.0e+0_fp - SPHU_kgkg )       &
!                                    * ( g0_100                          &
!                                    * State_Met%DELP(I,J,L) ) ) )          
! MOISTURE FIX:
       ! Area-independent conversion
       State_Chm%TRACERS(I,J,L,N) = State_Chm%TRACERS(I,J,L,N)          &
                                    * ( 1.0e+0_fp                       &      
                                    / ( g0_100                          &
                                    * State_Met%DELP_DRY(I,J,L) ) )          

    ENDDO
    ENDDO
    ENDDO
    ENDDO
    !$OMP END PARALLEL DO

    ! Update tracer units
    State_Chm%Trac_Units = 'kg/kg dry'

  END SUBROUTINE Convert_Kgm2_to_KgKgDry
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: convert_vvdry_to_kg
!
! !DESCRIPTION: Subroutine Convert\_VVDry\_to\_Kg converts the units of 
!  tracer concentrations from dry volume mixing ratio 
!  [mol tracer/mol dry air] to tracer mass per grid box [kg]. 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Convert_VVDry_to_Kg( am_I_Root, Input_Opt,    &
                                    State_Met, State_Chm, RC   ) 
!
! !USES:
!
    USE GIGC_Input_Opt_Mod, ONLY : OptInput
    USE GIGC_State_Met_Mod, ONLY : MetState
    USE GIGC_State_Chm_Mod, ONLY : ChmState
!
! !INPUT PARAMETERS: 
!
    LOGICAL,        INTENT(IN)    :: am_I_Root   ! Are we on the root CPU?
    TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
    TYPE(MetState), INTENT(IN)    :: State_Met   ! Meteorology state object
!
! !INPUT/OUTPUT PARAMETERS:
!
    ! Object containing tracer concentration
    TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry state object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
!
! !REMARKS:
! This routine replaces legacy routine CONVERT_UNITS and will be removed
! once GEOS-Chem is entirely area independent
!
! !REVISION HISTORY: 
!  08 Jan 2015 - E. Lundgren - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER :: I, J, L, N
    CHARACTER(LEN=255) :: MSG, LOC

    !====================================================================
    ! Convert_KgKgDry_to_Kg begins here!
    !====================================================================

    ! Assume success
    RC        =  GIGC_SUCCESS

    ! Verify correct initial units. If current units are unexpected,
    ! write error message and location to log, then pass failed RC
    ! to calling routine. 
    IF ( TRIM( State_Chm%Trac_Units ) /= 'v/v dry' ) THEN
       MSG = 'Incorrect initial units: ' // TRIM( State_Chm%Trac_Units )
       LOC = 'Routine Convert_VVDry_to_Kg in unitconv_mod.F90'
       CALL GIGC_Error( MSG, RC, LOC )
       RETURN
    ENDIF

    !====================================================================
    !
    !  The conversion is as follows:
    !
    !   mol tracer(N)                g/mol tracer(N)
    !   -----------  *  kg dry air * --------------      
    !   mol dry air                  g/mol dry air
    !
    !   = volume mixing ratio * dry air mass * MW tracer(N) / MW dry air 
    !   
    !   = kg tracer(N)
    !
    ! Therefore, with:
    !
    !  AD(I,J,L)   = grid box dry air mass [kg]
    !  TCVV(N)     = ratio of dry air and tracer molecular weights
    !     
    ! the conversion is:
    ! 
    !  TRACERS(I,J,L,N) [kg]
    !
    !    = TRACERS(I,J,L,N) [v/v] * AD(I,J,L) / TCVV (N)
    !                   
    !====================================================================

    !$OMP PARALLEL DO           &
    !$OMP DEFAULT( SHARED     ) &
    !$OMP PRIVATE( I, J, L, N ) 
    DO N = 1, Input_Opt%N_TRACERS
    DO L = 1, LLPAR
    DO J = 1, JJPAR
    DO I = 1, IIPAR
       State_Chm%TRACERS(I,J,L,N) = State_Chm%TRACERS(I,J,L,N)  &
                                  * State_Met%AD(I,J,L)         &
                                  / Input_Opt%TCVV(N)
    ENDDO
    ENDDO
    ENDDO
    ENDDO
    !$OMP END PARALLEL DO

    ! Update tracer units
    State_Chm%Trac_Units = 'kg'

  END SUBROUTINE Convert_VVDry_to_Kg
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: convert_kg_to_vvdry
!
! !DESCRIPTION: Subroutine Convert\_Kg\_to\_VVDry converts the units of 
!  tracer concentrations from tracer mass per grid box [kg] to dry volume 
!  mixing ratio [mol tracer/mol dry air].  
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Convert_Kg_to_VVDry( am_I_Root, Input_Opt,   &
                                    State_Met, State_Chm, RC   ) 
!
! !USES:
!
    USE GIGC_Input_Opt_Mod, ONLY : OptInput
    USE GIGC_State_Met_Mod, ONLY : MetState
    USE GIGC_State_Chm_Mod, ONLY : ChmState
!
! !INPUT PARAMETERS: 
!
    LOGICAL,        INTENT(IN)    :: am_I_Root   ! Are we on the root CPU?
    TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
    TYPE(MetState), INTENT(IN)    :: State_Met   ! Meteorology state object
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry state object 
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
!
! !REMARKS:
! This routine replaces legacy routine CONVERT_UNITS and will be removed
! once GEOS-Chem is entirely area independent
!
! !REVISION HISTORY: 
!  08 Jan 2015 - E. Lundgren - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER :: I, J, L, N
    INTEGER :: N_TRACERS
    CHARACTER(LEN=255) :: MSG, LOC

    !====================================================================
    ! Convert_Kg_to_VVDry begins here!
    !====================================================================

    ! Assume success
    RC        =  GIGC_SUCCESS

    ! Verify correct initial units. If current units are unexpected,
    ! write error message and location to log, then pass failed RC
    ! to calling routine. 
    IF ( TRIM( State_Chm%Trac_Units ) /= 'kg' ) THEN
       MSG = 'Incorrect initial units: ' // TRIM( State_Chm%Trac_Units )
       LOC = 'Routine Convert_Kg_to_VVDry in unitconv_mod.F90'
       CALL GIGC_Error( MSG, RC, LOC )
       RETURN
    ENDIF

    !====================================================================
    !
    !  The conversion is as follows:
    !
    !                         1           g/mol dry air
    !   kg tracer(N)  * --------------  * --------------    
    !                    kg dry air       g/mol tracer(N)      
    !
    !   = kg tracer(N) / dry air mass * MW dry air / MW tracer(N)
    !   
    !   = volume mixing ratio
    !
    ! Therefore, with:
    !
    !  AD(I,J,L)    = grid box dry air mass [kg]
    !  TCVV(N)      = ratio of dry air and tracer molecular weights
    !     
    ! the conversion is:
    ! 
    !  TRACERS(I,J,L,N) [v/v]
    !
    !    = TRACERS(I,J,L,N) [kg] * TCVV(N) / AD(I,J,L) 
    !                   
    !====================================================================

    !$OMP PARALLEL DO           &
    !$OMP DEFAULT( SHARED     ) &
    !$OMP PRIVATE( I, J, L, N ) 
    DO N = 1, Input_Opt%N_TRACERS
    DO L = 1, LLPAR
    DO J = 1, JJPAR
    DO I = 1, IIPAR
       State_Chm%TRACERS(I,J,L,N) = State_Chm%TRACERS(I,J,L,N)  &
                                  * Input_Opt%TCVV(N)           &
                                  / State_Met%AD(I,J,L)
    ENDDO
    ENDDO
    ENDDO
    ENDDO
    !$OMP END PARALLEL DO

    ! Update tracer units
    State_Chm%Trac_Units = 'v/v dry'

  END SUBROUTINE Convert_Kg_to_VVDry 
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: convert_kgkgdry_to_kg
!
! !DESCRIPTION: Subroutine Convert\_KgKgDry\_to\_Kg converts the units of 
!  tracer concentrations from dry mass mixing ratio 
!  [kg tracer/kg dry air] to tracer mass per grid box [kg]. 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Convert_KgKgDry_to_Kg( am_I_Root, Input_Opt,    &
                                    State_Met, State_Chm, RC   ) 
!
! !USES:
!
    USE GIGC_Input_Opt_Mod, ONLY : OptInput
    USE GIGC_State_Met_Mod, ONLY : MetState
    USE GIGC_State_Chm_Mod, ONLY : ChmState
!
! !INPUT PARAMETERS: 
!
    LOGICAL,        INTENT(IN)    :: am_I_Root   ! Are we on the root CPU?
    TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
    TYPE(MetState), INTENT(IN)    :: State_Met   ! Meteorology state object
!
! !INPUT/OUTPUT PARAMETERS:
!
    ! Object containing tracer concentration
    TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry state object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
!
! !REMARKS:
!
! !REVISION HISTORY: 
!  08 Jan 2015 - E. Lundgren - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER :: I, J, L, N
    CHARACTER(LEN=255) :: MSG, LOC

    !====================================================================
    ! Convert_KgKgDry_to_Kg begins here!
    !====================================================================

    ! Assume success
    RC        =  GIGC_SUCCESS

    ! Verify correct initial units. If current units are unexpected,
    ! write error message and location to log, then pass failed RC
    ! to calling routine. 
    IF ( TRIM( State_Chm%Trac_Units ) /= 'kg/kg dry' ) THEN
       MSG = 'Incorrect initial units: ' // TRIM( State_Chm%Trac_Units )
       LOC = 'Routine Convert_KgKgDry_to_Kg in unitconv_mod.F90'
       CALL GIGC_Error( MSG, RC, LOC )
       RETURN
    ENDIF

    !====================================================================
    !
    !  The conversion is as follows:
    !
    !   kg tracer(N)            
    !   -----------  *  kg dry air       
    !   kg dry air                   
    !
    !   = mass mixing ratio * dry air mass  
    !   
    !   = kg tracer(N)
    !
    ! Therefore, with:
    !
    !  AD(I,J,L)   = grid box dry air mass [kg]
    !     
    ! the conversion is:
    ! 
    !  TRACERS(I,J,L,N) [kg]
    !
    !    = TRACERS(I,J,L,N) [kg/kg] * AD(I,J,L)
    !                   
    !====================================================================

    !$OMP PARALLEL DO           &
    !$OMP DEFAULT( SHARED     ) &
    !$OMP PRIVATE( I, J, L, N ) 
    DO N = 1, Input_Opt%N_TRACERS
    DO L = 1, LLPAR
    DO J = 1, JJPAR
    DO I = 1, IIPAR
       State_Chm%TRACERS(I,J,L,N) = State_Chm%TRACERS(I,J,L,N) &
                                  * State_Met%AD(I,J,L)
    ENDDO
    ENDDO
    ENDDO
    ENDDO
    !$OMP END PARALLEL DO

    ! Update tracer units
    State_Chm%Trac_Units = 'kg'

  END SUBROUTINE Convert_KgKgDry_to_Kg
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: convert_kg_to_kgkgdry
!
! !DESCRIPTION: Subroutine Convert\_Kg\_to\_KgKgDry converts the units of 
!  tracer concentrations from tracer mass per grid box [kg] to dry mass 
!  mixing ratio [kg tracer/kg dry air].  
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Convert_Kg_to_KgKgDry( am_I_Root, Input_Opt,    &
                                    State_Met, State_Chm, RC   ) 
!
! !USES:
!
    USE GIGC_Input_Opt_Mod, ONLY : OptInput
    USE GIGC_State_Met_Mod, ONLY : MetState
    USE GIGC_State_Chm_Mod, ONLY : ChmState
!
! !INPUT PARAMETERS: 
!
    LOGICAL,        INTENT(IN)    :: am_I_Root   ! Are we on the root CPU?
    TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
    TYPE(MetState), INTENT(IN)    :: State_Met   ! Meteorology state object
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry state object 
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
!
! !REMARKS:
!
! !REVISION HISTORY: 
!  08 Jan 2015 - E. Lundgren - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER :: I, J, L, N
    CHARACTER(LEN=255) :: MSG, LOC

    !====================================================================
    ! Convert_Kg_to_KgKgDry begins here!
    !====================================================================

    ! Assume success
    RC        =  GIGC_SUCCESS

    ! Verify correct initial units. If current units are unexpected,
    ! write error message and location to log, then pass failed RC
    ! to calling routine. 
    IF ( TRIM( State_Chm%Trac_Units ) /= 'kg' ) THEN
       MSG = 'Incorrect initial units: ' // TRIM( State_Chm%Trac_Units )
       LOC = 'Routine Convert_Kg_to_KgKgDry in unitconv_mod.F90'
       CALL GIGC_Error( MSG, RC, LOC )
       RETURN
    ENDIF

    !====================================================================
    !
    !  The conversion is as follows:
    !
    !                         1          
    !   kg tracer(N)  * --------------      
    !                    kg dry air              
    !
    !   = kg tracer(N) / dry air mass
    !   
    !   = mass mixing ratio
    !
    ! Therefore, with:
    !
    !  AD(I,J,L)    = grid box dry air mass [kg]
    !     
    ! the conversion is:
    ! 
    !  TRACERS(I,J,L,N) [kg/kg]
    !
    !    = TRACERS(I,J,L,N) [kg] / AD(I,J,L) 
    !                   
    !====================================================================

    !$OMP PARALLEL DO           &
    !$OMP DEFAULT( SHARED     ) &
    !$OMP PRIVATE( I, J, L, N ) 
    DO N = 1, Input_Opt%N_TRACERS
    DO L = 1, LLPAR
    DO J = 1, JJPAR
    DO I = 1, IIPAR
       State_Chm%TRACERS(I,J,L,N) = State_Chm%TRACERS(I,J,L,N) &
                                  / State_Met%AD(I,J,L)
    ENDDO
    ENDDO
    ENDDO
    ENDDO
    !$OMP END PARALLEL DO

    ! Update tracer units
    State_Chm%Trac_Units = 'kg/kg dry'

  END SUBROUTINE Convert_Kg_to_KgKgDry 
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: convertbox_kgkgdry_to_kg
!
! !DESCRIPTION: Subroutine ConvertBox\_KgKgDry\_to\_Kg converts the units of 
!  tracer concentrations from dry mass mixing ratio [kg tracer/kg dry air]
!  to tracer mass per grid box [kg] for a single grid box. This routine is 
!  temporary during the unit transition of TOMAS to area-independence.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE ConvertBox_KgKgDry_to_Kg( am_I_Root, I, J, L, Input_Opt,  &
                                       State_Met, State_Chm, RC   ) 
!
! !USES:
!
    USE GIGC_Input_Opt_Mod, ONLY : OptInput
    USE GIGC_State_Met_Mod, ONLY : MetState
    USE GIGC_State_Chm_Mod, ONLY : ChmState
!
! !INPUT PARAMETERS: 
!
    LOGICAL,        INTENT(IN)    :: am_I_Root   ! Are we on the root CPU?
    INTEGER,        INTENT(IN)    :: I, J, L     ! Grid box indexes
    TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
    TYPE(MetState), INTENT(IN)    :: State_Met   ! Meteorology state object
!
! !INPUT/OUTPUT PARAMETERS:
!
    ! Object containing tracer concentration
    TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry state object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
!
! !REMARKS:
!  This routine is temporary and is only used for local conversion of tracer 
!  concentrations for use in TOMAS within wetscav_mod routine WASHOUT. 
!  That routine is called within a parallel do loop and therefore units can
!  only be converted per grid box to avoid excessive computation time. Also,
!  State_Chm%Trac_Units cannot be changed within the parallel do loop without
!  causing problems. It is therefore left out of this routine.
!
! !REVISION HISTORY: 
!  26 Jul 2016 - E. Lundgren - Initial version, an adaptation of 
!                              convert_kgkgtotal_to_kg
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER            :: N
    CHARACTER(LEN=255) :: MSG, LOC

    !=================================================================
    ! ConvertBox_KgKgDry_to_Kg begins here!
    !=================================================================

    ! Assume success
    RC        =  GIGC_SUCCESS

    !==============================================================
    !
    !  The conversion is as follows:
    !
    !   kg tracer(N)            
    !   -----------  *  kg dry air       
    !   kg dry air                   
    !
    !   = mass mixing ratio * dry air mass  
    !   
    !   = kg tracer(N)
    !
    ! Therefore, with:
    !
    !  AD(I,J,L)   = grid box dry air mass [kg] (wet)
    !     
    ! the conversion is:
    ! 
    !  TRACERS(I,J,L,N) [kg]
    !
    !    = TRACERS(I,J,L,N) [kg/kg] * AD(I,J,L)
    !                   
    !==============================================================

    !$OMP PARALLEL DO           &
    !$OMP DEFAULT( SHARED     ) &
    !$OMP PRIVATE( N ) 
    DO N = 1, Input_Opt%N_TRACERS
       State_Chm%TRACERS(I,J,L,N) = State_Chm%TRACERS(I,J,L,N)  &
                                  * State_Met%AD(I,J,L)
    ENDDO
    !$OMP END PARALLEL DO

  END SUBROUTINE ConvertBox_KgKgDry_to_Kg
!EOC

!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: convertbox_kg_to_kgkgdry
!
! !DESCRIPTION: Subroutine ConvertBox\_Kg\_to\_KgKgDry converts the units
!  of tracer concentrations from tracer mass per grid box [kg] to mass 
!  mixing ratio [kg tracer/kg dry air] for a single grid box.
!  This routine is temporary during the unit transition of TOMAS to 
!  area-independence.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE ConvertBox_Kg_to_KgKgDry( am_I_Root, I, J, L, Input_Opt, &
                                       State_Met, State_Chm, RC   ) 
!
! !USES:
!
    USE GIGC_Input_Opt_Mod, ONLY : OptInput
    USE GIGC_State_Met_Mod, ONLY : MetState
    USE GIGC_State_Chm_Mod, ONLY : ChmState
!
! !INPUT PARAMETERS: 
!
    LOGICAL,        INTENT(IN)    :: am_I_Root   ! Are we on the root CPU?
    INTEGER,        INTENT(IN)    :: I, J, L     ! Grid box indexes
    TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
    TYPE(MetState), INTENT(IN)    :: State_Met   ! Meteorology state object
!
! !INPUT/OUTPUT PARAMETERS:
!
    ! Object containing tracer concentration
    TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry state object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
!
! !REMARKS:
!  This routine is temporary and is only used for local conversion of tracer 
!  concentrations for use in TOMAS within wetscav_mod routine WASHOUT. 
!  That routine is called within a parallel do loop and therefore units can
!  only be converted per grid box to avoid excessive computation time. Also,
!  State_Chm%Trac_Units cannot be changed within the parallel do loop without
!  causing problems. It is therefore left out of this routine.
!
! !REVISION HISTORY: 
!  26 Jul 2016 - E. Lundgren - Initial version, an adaptation of 
!                              convert_kg_to_kgkgtotal
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER :: N
    CHARACTER(LEN=255) :: MSG, LOC

    !====================================================================
    ! Convert_Kg_to_KgKgDry begins here!
    !====================================================================

    ! Assume success
    RC        =  GIGC_SUCCESS

    !====================================================================
    !
    !  The conversion is as follows:
    !
    !                         1          
    !   kg tracer(N)  * --------------      
    !                    kg dry air              
    !
    !   = kg tracer(N) / dry air mass
    !   
    !   = mass mixing ratio
    !
    ! Therefore, with:
    !
    !  AD(I,J,L)    = grid box dry air mass [kg]
    !     
    ! the conversion is:
    ! 
    !  TRACERS(I,J,L,N) [kg/kg]
    !
    !    = TRACERS(I,J,L,N) [kg] / AD(I,J,L) 
    !                   
    !====================================================================

    !$OMP PARALLEL DO           &
    !$OMP DEFAULT( SHARED     ) &
    !$OMP PRIVATE( N ) 
    DO N = 1, Input_Opt%N_TRACERS
       State_Chm%TRACERS(I,J,L,N) = State_Chm%TRACERS(I,J,L,N) &
                                  / State_Met%AD(I,J,L)
    ENDDO
    !$OMP END PARALLEL DO

  END SUBROUTINE ConvertBox_Kg_to_KgKgDry
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: convertbox_kg_to_kgm2
!
! !DESCRIPTION: Subroutine ConvertBox\_Kg\_to\_kgm2 converts the units of 
!  mass [kg] to area density [kg/m2] for a single grid box.  This routine is
!  temporary during the unit transition of TOMAS to area-independence. 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE ConvertBox_Kg_to_Kgm2( am_I_Root, I, J, L, Input_Opt,      &
                                    State_Met, State_Chm, RC        )
!
! !USES:
!
    USE GIGC_Input_Opt_Mod, ONLY : OptInput
    USE GIGC_State_Chm_Mod, ONLY : ChmState
    USE GIGC_State_Met_Mod, ONLY : MetState
    USE PHYSCONSTANTS
!
! !INPUT PARAMETERS: 
!
    LOGICAL,        INTENT(IN)    :: am_I_Root     ! Are we on the root CPU?
    INTEGER,        INTENT(IN)    :: I, J, L       ! Grid box indexes
    TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
    TYPE(MetState), INTENT(IN)    :: State_Met     ! Meteorology state object
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry state object 
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC            ! Success or failure?
!
! !REMARKS:
!  This routine is temporary and is only used for local conversion of tracer 
!  concentrations for use in TOMAS within wetscav_mod routine WASHOUT. 
!  That routine is called within a parallel do loop and therefore units can
!  only be converted per grid box to avoid excessive computation time. Also,
!  State_Chm%Trac_Units cannot be changed within the parallel do loop without
!  causing problems. It is therefore left out of this routine.
!
! !REVISION HISTORY: 
!  14 Aug 2015 - E. Lundgren - Initial version
!  29 Sep 2015 - E. Lundgren - Now convert for a single grid box only
!  26 Jul 2016 - E. Lundgren - Rename from Convert_Kg_to_Kgm2
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER :: N
    CHARACTER(LEN=255) :: MSG, LOC

    !====================================================================
    ! ConvertBox_Kg_to_Kgm2 begins here!
    !====================================================================

    ! Assume success
    RC        =  GIGC_SUCCESS

    !====================================================================
    !
    !  The conversion is as follows:
    !
    !====================================================================

    !$OMP PARALLEL DO        &
    !$OMP DEFAULT( SHARED  ) &
    !$OMP PRIVATE( N ) 
    DO N = 1, Input_Opt%N_TRACERS
       State_Chm%TRACERS(I,J,L,N) = State_Chm%TRACERS(I,J,L,N)    &
                                     / State_Met%AREA_M2(I,J,1)
    ENDDO
    !$OMP END PARALLEL DO

  END SUBROUTINE ConvertBox_Kg_to_Kgm2
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: convertbox_kgm2_to_kg
!
! !DESCRIPTION: Subroutine ConvertBox\_Kgm2\_to\_Kg converts the units of area 
!  density [kg/m2] to mass [kg] for a single grid box. This routine is
!  temporary during the unit transition of TOMAS to area-independence.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE ConvertBox_Kgm2_to_Kg( am_I_Root, I, J, L, Input_Opt,   &
                                 State_Met, State_Chm, RC        )
!
! !USES:
!
    USE GIGC_Input_Opt_Mod, ONLY : OptInput
    USE GIGC_State_Chm_Mod, ONLY : ChmState
    USE GIGC_State_Met_Mod, ONLY : MetState
    USE PHYSCONSTANTS
!
! !INPUT PARAMETERS: 
!
    LOGICAL,        INTENT(IN)    :: am_I_Root     ! Are we on the root CPU?
    INTEGER,        INTENT(IN)    :: I, J, L       ! Grid box indexes
    TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
    TYPE(MetState), INTENT(IN)    :: State_Met     ! Meteorology state object
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry state object 
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC            ! Success or failure?
!
! !REMARKS:
!  This routine is temporary and is only used for local conversion of tracer 
!  concentrations for use in TOMAS within wetscav_mod routine WASHOUT. 
!  That routine is called within a parallel do loop and therefore units can
!  only be converted per grid box to avoid excessive computation time. Also,
!  State_Chm%Trac_Units cannot be changed within the parallel do loop without
!  causing problems. It is therefore left out of this routine.
!
! !REVISION HISTORY: 
!  14 Aug 2015 - E. Lundgren - Initial version
!  29 Sep 2015 - E. Lundgren - Now convert for a single grid box only
!  26 Jul 2016 - E. Lundgren - Rename from Convert_Kgm2_to_Kg
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER            :: N
    CHARACTER(LEN=255) :: MSG, LOC

    !====================================================================
    ! ConvertBox_Kgm2_to_Kg begins here!
    !====================================================================

    ! Assume success
    RC        =  GIGC_SUCCESS

    !====================================================================
    !
    !  The conversion is as follows:
    !
    !====================================================================

    !$OMP PARALLEL DO        &
    !$OMP DEFAULT( SHARED  ) &
    !$OMP PRIVATE( N ) 
    DO N = 1, Input_Opt%N_TRACERS
       State_Chm%TRACERS(I,J,L,N) = State_Chm%TRACERS(I,J,L,N)    &
                                     * State_Met%AREA_M2(I,J,1)
    ENDDO
    !$OMP END PARALLEL DO

  END SUBROUTINE ConvertBox_Kgm2_to_Kg
!EOC
END MODULE UnitConv_Mod
