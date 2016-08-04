!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: unitconv_mod.F90
!
! !DESCRIPTION: Module UNITCONV\_MOD contains routines which are used to 
!  convert the units of tracer or species concentrations between mass 
!  mixing ratio [kg/kg air], mass per grid box per area [kg/m2], molar 
!  mixing ratio [vol/vol], and molecular number density [molecules/cm3]. 
!  There are different conversion routines for dry air and total (wet) 
!  air mixing ratios. Conversions involving column area will be phased
!  out for grid-independent GEOS-Chem. 
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
  USE PHYSCONSTANTS
  USE GIGC_ErrCode_Mod
  USE ERROR_MOD
                    
  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! KG/KG TOTAL <-> KG/KG DRY 
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  ! Tracers
  ! Used for transport and convection
  PUBLIC  :: Convert_KgKgDry_to_KgKgTotal
  PUBLIC  :: Convert_KgKgTotal_to_KgKgDry

  ! Species
  PUBLIC  :: ConvertSpc_KgKgDry_to_KgKgTotal
  PUBLIC  :: ConvertSpc_KgKgTotal_to_KgKgDry

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! KG/KG DRY <-> V/V DRY
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  ! Tracers
  ! Used in DO_TEND in mixing
  PUBLIC  :: Convert_KgKgDry_to_VVDry
  PUBLIC  :: Convert_VVDry_to_KgKgDry

  ! Species
  PUBLIC  :: ConvertSpc_KgKgDry_to_VVDry
  PUBLIC  :: ConvertSpc_VVDry_to_KgKgDry

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! KG/KG DRY <-> KG/M2
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  ! Tracers
  ! Used for wet deposition, DO_TEND in mixing,
  ! and around AIRQNT and SET_H2O_TRAC in main
  PUBLIC  :: Convert_KgKgDry_to_Kgm2
  PUBLIC  :: Convert_kgm2_to_KgKgDry

  ! Species
  PUBLIC  :: ConvertSpc_KgKgDry_to_Kgm2
  PUBLIC  :: ConvertSpc_kgm2_to_KgKgDry

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! KG/KG DRY <-> MOLEC/CM3
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  ! Species
  PUBLIC  :: ConvertSpc_KgKgDry_to_MND
  PUBLIC  :: ConvertSpc_MND_to_KgKgDry

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! AREA-DEPENDENT (TO AND FROM KG)
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  ! Tracers
  ! v/v dry air <-> kg/grid box
  ! Temporarily replaces legacy CONVERT_UNITS
  ! Used in strat_chem_mod and sulfate_mod
  PUBLIC  :: Convert_VVDry_to_Kg
  PUBLIC  :: Convert_Kg_to_VVDry

  ! Species
  ! v/v dry air <-> kg/grid box
  PUBLIC  :: ConvertSpc_VVDry_to_Kg
  PUBLIC  :: ConvertSpc_Kg_to_VVDry

  ! Tracers
  ! kg/kg dry air <-> kg/grid box
  ! Used in aerosol_mod, tomas_mod, emissions_mod,
  ! strat_chem_mod, exchange_mod, rrtmg_rad_transfer_mod,
  ! chemistry_mod, sulfate_mod, and carbon_mod
  ! This is since RRTMG, TOMAS, exchange_mod, chemistry,
  ! and EMISSMERCURY are still in [kg]
  PUBLIC  :: Convert_KgKgDry_to_Kg
  PUBLIC  :: Convert_Kg_to_KgKgDry

  ! Species
  ! kg/kg dry air <-> kg/grid box
  PUBLIC  :: ConvertSpc_KgKgDry_to_Kg
  PUBLIC  :: ConvertSpc_Kg_to_KgKgDry

  ! Tracers
  ! kg/kg total air <-> kg/grid box (single box only)
  ! Used for TOMAS compatibility in WASHOUT
  PUBLIC  :: Convert_KgKgTotal_to_Kg
  PUBLIC  :: Convert_Kg_to_KgKgTotal

  ! Species
  ! kg/kg total air <-> kg/grid box (single box only)
  PUBLIC  :: ConvertSpc_KgKgTotal_to_Kg
  PUBLIC  :: ConvertSpc_Kg_to_KgKgTotal

  !Tracers
  ! kg <-> kg/m2 (single box only)
  ! Used for TOMAS compatibility in WASHOUT within wetscav_mod
  PUBLIC  :: Convert_Kg_to_Kgm2
  PUBLIC  :: Convert_Kgm2_to_Kg 

  !Species
  ! kg <-> kg/m2 (single box only)
  PUBLIC  :: ConvertSpc_Kg_to_Kgm2
  PUBLIC  :: ConvertSpc_Kgm2_to_Kg 

  ! Species
  ! kg <-> molec/cm3 dry air
  PUBLIC  :: ConvertSpc_Kg_to_MND
  PUBLIC  :: ConvertSpc_MND_to_Kg

!
! !REMARKS:
!  The routines in this module are used to convert the units of tracer 
!  or species concentrations in various GEOS-Chem routines.
!
! !REVISION HISTORY:
!  23 Jun 2015 - E. Lundgren - Initial version
!  13 Aug 2015 - E. Lundgren - Add tracer unit error handling
!  29 Sep 2015 - E. Lundgren - Adjust some of the unit conversions to/from kg
!                              to be for a single grid box for TOMAS
!  21 Jul 2016 - E. Lundgren - Add species unit conversion routines
!EOP
!------------------------------------------------------------------------------
!BOC
CONTAINS
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: convert_kgkgdry_to_kgkgtotal
!
! !DESCRIPTION: Subroutine Convert\_KgKgDry\_to\_KgKgTotal converts the units 
!  of tracer mass mixing ratio from tracer mass per dry air mass [kg tracer/
!  kg dry air] to tracer mass per moist air mass [kg tracer/kg moist air]. 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Convert_KgKgDry_to_KgKgTotal( am_I_Root, Input_Opt,  &
                                           State_Met, State_Chm, RC ) 
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
!  16 Apr 2015 - E. Lundgren - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER :: I, J, L, N
    CHARACTER(LEN=255) :: MSG, LOC

    !====================================================================
    ! Convert_KgKgDry_to_KgKgTotal begins here!
    !====================================================================

    ! Assume success
    RC        =  GIGC_SUCCESS

    ! Verify correct initial units. If current units are unexpected,
    ! write error message and location to log, then pass failed RC
    ! to calling routine. 
    IF ( TRIM( State_Chm%Trac_Units ) /= 'kg/kg dry' ) THEN
       MSG = 'Incorrect initial units: ' // TRIM( State_Chm%Trac_Units )
       LOC = 'Routine Convert_KgKgDry_to_KgKgTotal in unitconv_mod.F90'
       CALL GIGC_Error( MSG, RC, LOC )
    ENDIF

    !====================================================================
    !
    !  The conversion is as follows:
    !
    !   kg tracer(N)     kg dry air     
    !   ------------  *  ------------    
    !   kg dry air       kg moist air            
    !
    !   = kg tracer(N) / kg moist air
    !
    ! Therefore, with kg dry air / kg moist air defined as the
    !  complement of specific humidity (kg water vapor / kg moist air),
    !  the conversion is:
    ! 
    ! TRACERS(I,J,L,N) [kg/kg moist]
    !
    !    = TRACERS(I,J,L,N) [kg/kg dry] * ( 1 - SPHU(I,J,L) )
    !
    ! Note that State_Met%SPHU is in units of [g/kg] and so must
    ! be converted to [kg/kg]. 
    !                   
    !====================================================================

    !$OMP PARALLEL DO           &
    !$OMP DEFAULT( SHARED     ) &
    !$OMP PRIVATE( I, J, L, N ) 
    DO N = 1, Input_Opt%N_TRACERS
    DO L = 1, LLPAR
    DO J = 1, JJPAR
    DO I = 1, IIPAR
       State_Chm%TRACERS(I,J,L,N) = State_Chm%TRACERS(I,J,L,N)        &
                                  * ( 1.0e+0_fp - State_Met%SPHU(I,J,L) &
                                  * 1.e-3_fp )
    ENDDO
    ENDDO
    ENDDO
    ENDDO
    !$OMP END PARALLEL DO

    ! Update tracer units
    State_Chm%Trac_Units = 'kg/kg total'

  END SUBROUTINE Convert_KgKgDry_to_KgKgTotal
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: convert_kgkgtotal_to_kgkgdry
!
! !DESCRIPTION: Subroutine Convert\_KgKgTotal\_to\_KgKgDry converts the units 
!  of tracer mass mixing ratio from tracer mass per moist air mass [kg tracer/
!  kg moist air] to tracer mass per dry air mass [kg tracer/kg dry air]. 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Convert_KgKgTotal_to_KgKgDry( am_I_Root, Input_Opt, &
                                           State_Met, State_Chm, RC ) 
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
!  16 Apr 2015 - E. Lundgren - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER :: I, J, L, N
    CHARACTER(LEN=255) :: MSG, LOC

    !====================================================================
    ! Convert_KgKgTotal_to_KgKgDry begins here!
    !====================================================================

    ! Assume success
    RC        =  GIGC_SUCCESS

    ! Verify correct initial units. If current units are unexpected,
    ! write error message and location to log, then pass failed RC
    ! to calling routine. 
    IF ( TRIM( State_Chm%Trac_Units ) /= 'kg/kg total' ) THEN
       MSG = 'Incorrect initial units: ' // TRIM( State_Chm%Trac_Units )
       LOC = 'Routine Convert_KgKgTotal_to_KgKgDry in unitconv_mod.F90'
       CALL GIGC_Error( MSG, RC, LOC )
       RETURN
    ENDIF

    !====================================================================
    !
    !  The conversion is as follows:
    !
    !   kg tracer(N)     kg moist air     
    !   ------------  *  ------------    
    !   kg moist air     kg dry air            
    !
    !   = kg tracer(N) / kg dry air
    !
    ! Therefore, with kg dry air / kg moist air defined as the
    !  complement of specific humidity (kg water vapor / kg moist air),
    !  the conversion is:
    ! 
    ! TRACERS(I,J,L,N) [kg/kg dry]
    !
    !    = TRACERS(I,J,L,N) [kg/kg moist] / ( 1 - SPHU(I,J,L) )
    !      
    ! Note that State_Met%SPHU is in units of [g/kg] and so must
    ! be converted to [kg/kg]. 
    !                                
    !====================================================================

    !$OMP PARALLEL DO           &
    !$OMP DEFAULT( SHARED     ) &
    !$OMP PRIVATE( I, J, L, N ) 
    DO N = 1, Input_Opt%N_TRACERS
    DO L = 1, LLPAR
    DO J = 1, JJPAR
    DO I = 1, IIPAR
       State_Chm%TRACERS(I,J,L,N) = State_Chm%TRACERS(I,J,L,N)          &
                                  / ( 1.0e+0_fp - State_Met%SPHU(I,J,L) &
                                  * 1.e-3_fp )
    ENDDO
    ENDDO
    ENDDO
    ENDDO
    !$OMP END PARALLEL DO

    ! Update tracer units
    State_Chm%Trac_Units = 'kg/kg dry'
  
  END SUBROUTINE Convert_KgKgTotal_to_KgKgDry
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: convertspc_kgkgdry_to_kgkgtotal
!
! !DESCRIPTION: Subroutine ConvertSpc\_KgKgDry\_to\_KgKgTotal converts the 
!  units of species mass mixing ratio from species mass per dry air mass 
!  [kg species/kg dry air] to species mass per moist air mass 
!  [kg species/kg moist air]. 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE ConvertSpc_KgKgDry_to_KgKgTotal( am_I_Root,  State_Met, &
                                              State_Chm, RC ) 
!
! !USES:
!
    USE GIGC_State_Met_Mod, ONLY : MetState
    USE GIGC_State_Chm_Mod, ONLY : ChmState
!
! !INPUT PARAMETERS: 
!
    LOGICAL,        INTENT(IN)    :: am_I_Root   ! Are we on the root CPU?
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
!  21 Jul 2016 - E. Lundgren - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER :: I, J, L, N
    CHARACTER(LEN=255) :: MSG, LOC

    !====================================================================
    ! ConvertSpc_KgKgDry_to_KgKgTotal begins here!
    !====================================================================

    ! Assume success
    RC        =  GIGC_SUCCESS

    ! Verify correct initial units. If current units are unexpected,
    ! write error message and location to log, then pass failed RC
    ! to calling routine. 
    IF ( TRIM( State_Chm%Spc_Units ) /= 'kg/kg dry' ) THEN
       MSG = 'Incorrect initial units: ' // TRIM( State_Chm%Spc_Units )
       LOC = 'Routine ConvertSpc_KgKgDry_to_KgKgTotal in unitconv_mod.F90'
       CALL GIGC_Error( MSG, RC, LOC )
    ENDIF

    !====================================================================
    !
    !  The conversion is as follows:
    !
    !   kg species(N)     kg dry air     
    !   ------------  *  ------------    
    !   kg dry air       kg moist air            
    !
    !   = kg species(N) / kg moist air
    !
    ! Therefore, with kg dry air / kg moist air defined as the
    !  complement of specific humidity (kg water vapor / kg moist air),
    !  the conversion is:
    ! 
    ! Species(I,J,L,N) [kg/kg moist]
    !
    !    = Species(I,J,L,N) [kg/kg dry] * ( 1 - SPHU(I,J,L) )
    !
    ! Note that State_Met%SPHU is in units of [g/kg] and so must
    ! be converted to [kg/kg]. 
    !                   
    !====================================================================

    !$OMP PARALLEL DO           &
    !$OMP DEFAULT( SHARED     ) &
    !$OMP PRIVATE( I, J, L, N ) 
    DO N = 1, State_Chm%nSpecies
    DO L = 1, LLPAR
    DO J = 1, JJPAR
    DO I = 1, IIPAR
       State_Chm%Species(I,J,L,N) = State_Chm%Species(I,J,L,N)        &
                                  * ( 1.0e+0_fp - State_Met%SPHU(I,J,L) &
                                  * 1.e-3_fp )
    ENDDO
    ENDDO
    ENDDO
    ENDDO
    !$OMP END PARALLEL DO

    ! Update species units
    State_Chm%Spc_Units = 'kg/kg total'

  END SUBROUTINE ConvertSpc_KgKgDry_to_KgKgTotal
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: convertspc_kgkgtotal_to_kgkgdry
!
! !DESCRIPTION: Subroutine ConvertSpc\_KgKgTotal\_to\_KgKgDry converts the 
!  units of species mass mixing ratio from species mass per moist air mass 
!  [kg species/ kg moist air] to species mass per dry air mass 
!  [kg species/kg dry air]. 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE ConvertSpc_KgKgTotal_to_KgKgDry( am_I_Root, State_Met, &
                                              State_Chm, RC ) 
!
! !USES:
!
    USE GIGC_State_Met_Mod, ONLY : MetState
    USE GIGC_State_Chm_Mod, ONLY : ChmState
!
! !INPUT PARAMETERS: 
!
    LOGICAL,        INTENT(IN)    :: am_I_Root   ! Are we on the root CPU?
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
!  21 Jul 2016 - E. Lundgren - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER :: I, J, L, N
    CHARACTER(LEN=255) :: MSG, LOC

    !====================================================================
    ! Convert_KgKgTotal_to_KgKgDry begins here!
    !====================================================================

    ! Assume success
    RC        =  GIGC_SUCCESS

    ! Verify correct initial units. If current units are unexpected,
    ! write error message and location to log, then pass failed RC
    ! to calling routine. 
    IF ( TRIM( State_Chm%Spc_Units ) /= 'kg/kg total' ) THEN
       MSG = 'Incorrect initial units: ' // TRIM( State_Chm%Spc_Units )
       LOC = 'Routine ConvertSpc_KgKgTotal_to_KgKgDry in unitconv_mod.F90'
       CALL GIGC_Error( MSG, RC, LOC )
       RETURN
    ENDIF

    !====================================================================
    !
    !  The conversion is as follows:
    !
    !   kg species(N)    kg moist air     
    !   ------------  *  ------------    
    !   kg moist air     kg dry air            
    !
    !   = kg species(N) / kg dry air
    !
    ! Therefore, with kg dry air / kg moist air defined as the
    !  complement of specific humidity (kg water vapor / kg moist air),
    !  the conversion is:
    ! 
    ! Species(I,J,L,N) [kg/kg dry]
    !
    !    = Species(I,J,L,N) [kg/kg moist] / ( 1 - SPHU(I,J,L) )
    !      
    ! Note that State_Met%SPHU is in units of [g/kg] and so must
    ! be converted to [kg/kg]. 
    !                                
    !====================================================================

    !$OMP PARALLEL DO           &
    !$OMP DEFAULT( SHARED     ) &
    !$OMP PRIVATE( I, J, L, N ) 
    DO N = 1, State_Chm%nSpecies
    DO L = 1, LLPAR
    DO J = 1, JJPAR
    DO I = 1, IIPAR
       State_Chm%Species(I,J,L,N) = State_Chm%Species(I,J,L,N)          &
                                  / ( 1.0e+0_fp - State_Met%SPHU(I,J,L) &
                                  * 1.e-3_fp )
    ENDDO
    ENDDO
    ENDDO
    ENDDO
    !$OMP END PARALLEL DO

    ! Update species units
    State_Chm%Spc_Units = 'kg/kg dry'
  
  END SUBROUTINE ConvertSpc_KgKgTotal_to_KgKgDry
!EOC
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
! !IROUTINE: convertspc_kgkgdry_to_vvdry
!
! !DESCRIPTION: Subroutine ConvertSpc\_KgKgDry\_to\_VVDry converts the 
!  units of species concentrations from mass mixing ratio (KGKG) [kg/kg] to 
!  volume ratio (VR) [vol/vol] (same as molar ratio [mol/mol]). 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE ConvertSpc_KgKgDry_to_VVDry( am_I_Root, State_Chm, RC ) 
!
! USES: 
!
    USE GIGC_State_Chm_Mod, ONLY : ChmState
    USE Species_Mod,        ONLY : Species
!
! !INPUT PARAMETERS: 
!
    LOGICAL,        INTENT(IN)    :: am_I_Root   ! Are we on the root CPU?
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
!  21 Jul 2016 - E. Lundgren - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER                :: I, J, L, N
    CHARACTER(LEN=255)     :: MSG, LOC
    REAL(fp)               :: MolecRatio, MW_G
    TYPE(Species), POINTER :: ThisSpc
    
    !====================================================================
    ! ConvertSpc_KgKgDry_to_VVDry begins here!
    !====================================================================

    ! Assume success
    RC = GIGC_SUCCESS

    ! Initialize pointer
    ThisSpc => NULL()

    ! Verify correct initial units. If current units are unexpected,
    ! write error message and location to log, then pass failed RC
    ! to calling routine. 
    IF ( TRIM( State_Chm%Spc_Units ) /= 'kg/kg dry' ) THEN
       MSG = 'Incorrect initial units: ' // TRIM( State_Chm%Spc_Units )
       LOC = 'Routine ConvertSpc_KgKgDry_to_VVDry in unitconv_mod.F90'
       CALL GIGC_Error( MSG, RC, LOC )
       RETURN
    ENDIF

    !====================================================================
    !
    !  The conversion is as follows:
    !
    !   kg species(N)   g dry air      mol species(N)    
    !   ------------- * ----------  *  -------------  
    !     kg air        mol air         g species(N)          
    !
    !   = mass mixing ratio * ratio of air to species molecular weights  
    !   
    !   = molar ratio
    !
    ! Therefore, with:
    !
    !  AIRMW   = dry air molecular wt [g/mol]
    !  MW_G(N) = species molecular wt [g/mol]
    !     
    ! the conversion is:
    ! 
    !  Species(I,J,L,N) [vol/vol]
    !
    !    = Species(I,J,L,N) [kg/kg] * ( AIRMW / MW_G(N) )
    !                   
    !====================================================================
 
    !$OMP PARALLEL DO           &
    !$OMP DEFAULT( SHARED     ) &
    !$OMP PRIVATE( I, J, L, N, MW_G ) 
    DO N = 1, State_Chm%nSpecies

       ! Get info about this species from the species database
       ThisSpc => State_Chm%SpcData(N)%Info
       MolecRatio = State_Chm%SpcData(N)%Info%MolecRatio ! mol C / mol spc?
       MW_G       = State_Chm%SpcData(N)%Info%emMW_g
    
       DO L = 1, LLPAR
       DO J = 1, JJPAR
       DO I = 1, IIPAR
          State_Chm%Species(I,J,L,N) = State_Chm%Species(I,J,L,N)  &
                                     * ( AIRMW / MW_G )
       ENDDO
       ENDDO
       ENDDO

       ! Free pointer
       ThisSpc => NULL()

    ENDDO
    !$OMP END PARALLEL DO

    ! Update species units
    State_Chm%Spc_Units = 'v/v dry'

  END SUBROUTINE ConvertSpc_KgKgDry_to_VVDry
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: convertspc_vvdry_to_kgkgdry
!
! !DESCRIPTION: Subroutine ConvertSpc\_VVDry\_to\_KgKgDry converts the 
!  units of species concentrations from volume ratio (VR) [vol/vol] (same 
!  as molar mixing ratio [mol/mol]) to mass mixing ratio [kg/kg]. 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE ConvertSpc_VVDry_to_KgKgDry( am_I_Root, State_Chm, RC ) 
!
! USES: 
!
    USE GIGC_State_Chm_Mod, ONLY : ChmState
    USE Species_Mod,        ONLY : Species
!
! !INPUT PARAMETERS: 
!
    LOGICAL,        INTENT(IN)    :: am_I_Root   ! Are we on the root CPU?
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
!  21 Jul 2016 - E. Lundgren - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER                :: I, J, L, N
    CHARACTER(LEN=255)     :: MSG, LOC
    REAL(fp)               :: MolecRatio, MW_G
    TYPE(Species), POINTER :: ThisSpc

    !====================================================================
    ! ConvertSpc_VVDry_to_KgKgDry begins here!
    !=================================================================

    ! Assume success
    RC = GIGC_SUCCESS

    ! Initialize pointer
    ThisSpc => NULL()

    ! Verify correct initial units. If current units are unexpected,
    ! write error message and location to log, then pass failed RC
    ! to calling routine. 
    IF ( TRIM( State_Chm%Spc_Units ) /= 'v/v dry' ) THEN
       MSG = 'Incorrect initial units: ' // TRIM( State_Chm%Spc_Units )
       LOC = 'Routine ConvertSpc_VVDry_to_KgKgDry in unitconv_mod.F90'
       CALL GIGC_Error( MSG, RC, LOC )
       RETURN
    ENDIF

    !==============================================================
    !
    !  The conversion is as follows:
    !
    !   mol species(N)  mol dry air     g species(N)         
    !   -----------  * -----------  *  -------------  
    !     mol air       g dry air      mol species(N)           
    !
    !   = volume ratio / ratio of air to species molecular wts  
    !   
    !   = mass mixing ratio ([g/g] is equivalent to [kg/kg])
    !
    ! Therefore, with:
    !
    !  AIRMW   = dry air molecular wt [g/mol]
    !  MW_G(N) = species molecular wt [g/mol]
    !     
    ! the conversion is:
    ! 
    !  Species(I,J,L,N) [vol/vol]
    !
    !    = Species(I,J,L,N) [kg/kg] / ( AIRMW / MW_G(N) )
    !                   
    !==============================================================

    !$OMP PARALLEL DO           &
    !$OMP DEFAULT( SHARED     ) &
    !$OMP PRIVATE( I, J, L, N, MW_G ) 
    DO N = 1, State_Chm%nSpecies

       ! Get info about this species from the species database
       ThisSpc => State_Chm%SpcData(N)%Info
       MolecRatio = State_Chm%SpcData(N)%Info%MolecRatio ! mol C / mol spc?
       MW_G       = State_Chm%SpcData(N)%Info%emMW_g

       DO L = 1, LLPAR
       DO J = 1, JJPAR
       DO I = 1, IIPAR
         State_Chm%Species(I,J,L,N) = State_Chm%Species(I,J,L,N)   &
                                      / ( AIRMW / MW_G )
       ENDDO
       ENDDO
       ENDDO

       ! Free pointer
       ThisSpc => NULL()

    ENDDO
    !$OMP END PARALLEL DO

    ! Update species units
    State_Chm%Spc_Units = 'kg/kg dry'

    END SUBROUTINE ConvertSpc_VVDry_to_KgKgDry
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
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER            :: I, J, L, N
    CHARACTER(LEN=255) :: MSG, LOC
    REAL(fp)           :: SPHU_kgkg

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
    !   kg tracer      Delta P [hPa]   100 [Pa]     kg dry air 
    !   -----------  * ------------- * -------- * --------------     
    !   kg dry air     g [m/s2]        [hPa]      kg total air  
    !
    !   = kg / m2
    !
    ! where:
    !
    !  Delta P = edge pressure difference across level
    !  g = acceleration due to gravity
    !  kg dry air / kg total air  = 1 - specific humidity
    !     
    !====================================================================

    !$OMP PARALLEL DO        &
    !$OMP DEFAULT( SHARED  ) &
    !$OMP PRIVATE( I, J, L, SPHU_kgkg ) 
    DO N = 1, Input_Opt%N_TRACERS
    DO L = 1, LLPAR
    DO J = 1, JJPAR
    DO I = 1, IIPAR

       ! Convert specific humidity from [g/kg] to [kg/kg]
       SPHU_kgkg = State_Met%SPHU(I,J,L) * 1.0e-3_fp 

       ! Area-independent conversion
       State_Chm%TRACERS(I,J,L,N) = State_Chm%TRACERS(I,J,L,N)          &
                                    * ( ( 1.0e+0_fp - SPHU_kgkg )       &
                                    * ( g0_100                          &
                                    * State_Met%DELP(I,J,L) ) )           

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
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER            :: I, J, L, N
    CHARACTER(LEN=255) :: MSG, LOC
    REAL(fp)           :: SPHU_kgkg

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
    !   kg tracer(N)   g [m/s2]        [hPa]      kg total air   
    !   -----------  * ------------- * -------- * --------------     
    !        m2        Delta P [hPa]   100 [Pa]    kg dry air
    !
    !   = kg tracer(N) / kg dry air
    !
    ! where:
    !
    !  Delta P = edge pressure difference across level
    !  g = acceleration due to gravity
    !  kg dry air / kg total air  = 1 - specific humidity
    !     
    !====================================================================

    !$OMP PARALLEL DO           &
    !$OMP DEFAULT( SHARED     ) &
    !$OMP PRIVATE( I, J, L, N, SPHU_kgkg ) 
    DO N = 1, Input_Opt%N_TRACERS
    DO L = 1, LLPAR
    DO J = 1, JJPAR
    DO I = 1, IIPAR

       ! Convert specific humidity from [g/kg] to [kg/kg]
       SPHU_kgkg = State_Met%SPHU(I,J,L) * 1.0e-3_fp 

       ! Area-independent conversion
       State_Chm%TRACERS(I,J,L,N) = State_Chm%TRACERS(I,J,L,N)          &
                                    * ( 1.0e+0_fp                       &      
                                    / ( ( 1.0e+0_fp - SPHU_kgkg )       &
                                    * ( g0_100                          &
                                    * State_Met%DELP(I,J,L) ) ) )          

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
! !IROUTINE: convertspc_kgkgdry_to_kgm2
!
! !DESCRIPTION: Subroutine ConvertSpc\_kgkgdry\_to\_kgm2 converts the units of 
!  a 3D array from dry mass mixing ratio [kg/kg dry air] to area density 
!  [kg/m2].  
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE ConvertSpc_KgKgDry_to_Kgm2( am_I_Root, State_Met,  &
                                         State_Chm, RC             )
!
! !USES:
!
    USE GIGC_State_Chm_Mod, ONLY : ChmState
    USE GIGC_State_Met_Mod, ONLY : MetState
!
! !INPUT PARAMETERS: 
!
    LOGICAL,        INTENT(IN)    :: am_I_Root     ! Are we on the root CPU?
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
!  21 Jul 2016 - E. Lundgren - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER            :: I, J, L, N
    CHARACTER(LEN=255) :: MSG, LOC
    REAL(fp)           :: SPHU_kgkg

    !====================================================================
    ! ConvertSpc_KgKgDry_to_Kgm2 begins here!
    !====================================================================

    ! Assume success
    RC        =  GIGC_SUCCESS

    ! Verify correct initial units. If current units are unexpected,
    ! write error message and location to log, then pass failed RC
    ! to calling routine. 
    IF ( TRIM( State_Chm%Spc_Units ) /= 'kg/kg dry' ) THEN
       MSG = 'Incorrect initial units: ' // TRIM( State_Chm%Spc_Units )
       LOC = 'Routine ConvertSpc_KgKgDry_to_Kgm2 in unitconv_mod.F90'
       CALL GIGC_Error( MSG, RC, LOC )
       RETURN
    ENDIF

    !====================================================================
    !
    !  The conversion is as follows:
    !
    !   kg species     Delta P [hPa]   100 [Pa]     kg dry air 
    !   -----------  * ------------- * -------- * --------------     
    !   kg dry air     g [m/s2]        [hPa]      kg total air  
    !
    !   = kg / m2
    !
    ! where:
    !
    !  Delta P = edge pressure difference across level
    !  g = acceleration due to gravity
    !  kg dry air / kg total air  = 1 - specific humidity
    !     
    !====================================================================

    !$OMP PARALLEL DO        &
    !$OMP DEFAULT( SHARED  ) &
    !$OMP PRIVATE( I, J, L, N, SPHU_kgkg ) 
    DO N = 1, State_Chm%nSpecies
    DO L = 1, LLPAR
    DO J = 1, JJPAR
    DO I = 1, IIPAR

       ! Convert specific humidity from [g/kg] to [kg/kg]
       SPHU_kgkg = State_Met%SPHU(I,J,L) * 1.0e-3_fp 

       ! Area-independent conversion
       State_Chm%Species(I,J,L,N) = State_Chm%Species(I,J,L,N)          &
                                    * ( ( 1.0e+0_fp - SPHU_kgkg )       &
                                    * ( g0_100                          &
                                    * State_Met%DELP(I,J,L) ) )           

    ENDDO
    ENDDO
    ENDDO
    ENDDO
    !$OMP END PARALLEL DO

    ! Update species units
    State_Chm%Spc_Units = 'kg/m2'

  END SUBROUTINE ConvertSpc_KgKgDry_to_Kgm2
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: convertspc_kgm2_to_kgkgdry
!
! !DESCRIPTION: Subroutine ConvertSpc\_Kgm2\_to\_kgkgdry converts the units of 
!  species concentrations from area density [kg/m2] to dry mass mixing ratio 
!  [kg/kg dry air].  
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE ConvertSpc_Kgm2_to_KgKgDry( am_I_Root, State_Met, &
                                         State_Chm, RC          )
!
! !USES:
!
    USE GIGC_State_Met_Mod, ONLY : MetState
    USE GIGC_State_Chm_Mod, ONLY : ChmState
!
! !INPUT PARAMETERS: 
!
    LOGICAL,        INTENT(IN)    :: am_I_Root   ! Are we on the root CPU?
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
!  21 Jul 2016 - E. Lundgren - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER            :: I, J, L, N
    CHARACTER(LEN=255) :: MSG, LOC
    REAL(fp)           :: SPHU_kgkg

    !====================================================================
    ! ConvertSpc_Kgm2_to_KgKgDry begins here!
    !====================================================================

    ! Assume success
    RC        =  GIGC_SUCCESS

    ! Verify correct initial units. If current units are unexpected,
    ! write error message and location to log, then pass failed RC
    ! to calling routine. 
    IF ( TRIM( State_Chm%Spc_Units ) /= 'kg/m2' ) THEN
       MSG = 'Incorrect initial units: ' // TRIM( State_Chm%Spc_Units )
       LOC = 'Routine ConvertSpc_Kgm2_to_KgKgDry in unitconv_mod.F90'
       CALL GIGC_Error( MSG, RC, LOC )
       RETURN
    ENDIF

    !====================================================================
    !
    !  The conversion is as follows:
    !
    !   kg species(N)   g [m/s2]        [hPa]      kg total air   
    !   -----------  * ------------- * -------- * --------------     
    !        m2        Delta P [hPa]   100 [Pa]    kg dry air
    !
    !   = kg species(N) / kg dry air
    !
    ! where:
    !
    !  Delta P = edge pressure difference across level
    !  g = acceleration due to gravity
    !  kg dry air / kg total air  = 1 - specific humidity
    !     
    !====================================================================

    !$OMP PARALLEL DO           &
    !$OMP DEFAULT( SHARED     ) &
    !$OMP PRIVATE( I, J, L, N, SPHU_kgkg ) 
    DO N = 1, State_Chm%nSpecies
    DO L = 1, LLPAR
    DO J = 1, JJPAR
    DO I = 1, IIPAR

       ! Convert specific humidity from [g/kg] to [kg/kg]
       SPHU_kgkg = State_Met%SPHU(I,J,L) * 1.0e-3_fp 

       ! Area-independent conversion
       State_Chm%Species(I,J,L,N) = State_Chm%Species(I,J,L,N)          &
                                    * ( 1.0e+0_fp                       &      
                                    / ( ( 1.0e+0_fp - SPHU_kgkg )       &
                                    * ( g0_100                          &
                                    * State_Met%DELP(I,J,L) ) ) )          

    ENDDO
    ENDDO
    ENDDO
    ENDDO
    !$OMP END PARALLEL DO

    ! Update species units
    State_Chm%Spc_Units = 'kg/kg dry'

  END SUBROUTINE ConvertSpc_Kgm2_to_KgKgDry
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
! !IROUTINE: convertspc_vvdry_to_kg
!
! !DESCRIPTION: Subroutine ConvertSpc\_VVDry\_to\_Kg converts the units of 
!  species concentrations from dry volume mixing ratio 
!  [mol species/mol dry air] to species mass per grid box [kg]. 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE ConvertSpc_VVDry_to_Kg( am_I_Root, State_Met,  &
                                     State_Chm, RC   ) 
!
! !USES:
!
    USE GIGC_State_Met_Mod, ONLY : MetState
    USE GIGC_State_Chm_Mod, ONLY : ChmState
    USE Species_Mod,        ONLY : Species
!
! !INPUT PARAMETERS: 
!
    LOGICAL,        INTENT(IN)    :: am_I_Root   ! Are we on the root CPU?
    TYPE(MetState), INTENT(IN)    :: State_Met   ! Meteorology state object
!
! !INPUT/OUTPUT PARAMETERS:
!
    ! Object containing species concentration
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
!  21 Jul 2016 - E. Lundgren - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER                :: I, J, L, N
    CHARACTER(LEN=255)     :: MSG, LOC
    REAL(fp)               :: MW_G
    TYPE(Species), POINTER :: ThisSpc


    !====================================================================
    ! ConvertSpc_KgKgDry_to_Kg begins here!
    !====================================================================

    ! Assume success
    RC = GIGC_SUCCESS

    ! Verify correct initial units. If current units are unexpected,
    ! write error message and location to log, then pass failed RC
    ! to calling routine. 
    IF ( TRIM( State_Chm%Spc_Units ) /= 'v/v dry' ) THEN
       MSG = 'Incorrect initial units: ' // TRIM( State_Chm%Spc_Units )
       LOC = 'Routine ConvertSpc_VVDry_to_Kg in unitconv_mod.F90'
       CALL GIGC_Error( MSG, RC, LOC )
       RETURN
    ENDIF

    !====================================================================
    !
    !  The conversion is as follows:
    !
    !   mol species(N)                g/mol species(N)
    !   -------------  * kg dry air * ----------------      
    !   mol dry air                   g/mol dry air
    !
    !   = volume mixing ratio * dry air mass * MW species(N) / MW dry air 
    !   
    !   = kg species(N)
    !
    ! Therefore, with:
    !
    !  AD(I,J,L) = grid box dry air mass [kg]
    !  AIRMW     = dry air molecular wt [g/mol]
    !  MW_G(N)   = species molecular wt [g/mol]
    !     
    ! the conversion is:
    ! 
    !  SPECIES(I,J,L,N) [kg]
    !
    !    = SPECIES(I,J,L,N) [v/v] * AD(I,J,L) /  ( AIRMW / MW_G(N) )
    !                   
    !====================================================================

    !$OMP PARALLEL DO           &
    !$OMP DEFAULT( SHARED     ) &
    !$OMP PRIVATE( I, J, L, N, MW_G ) 
    DO N = 1, State_Chm%nSpecies

       ! Get info about this species from the species database
       ThisSpc    => State_Chm%SpcData(N)%Info
       MW_G      = State_Chm%SpcData(N)%Info%emMW_g

       DO L = 1, LLPAR
       DO J = 1, JJPAR
       DO I = 1, IIPAR
          State_Chm%Species(I,J,L,N) = State_Chm%Species(I,J,L,N)  &
                                       * State_Met%AD(I,J,L)       &
                                       / ( AIRMW / MW_G )
       ENDDO
       ENDDO
       ENDDO

    ENDDO
    !$OMP END PARALLEL DO

    ! Update species units
    State_Chm%Spc_Units = 'kg'

  END SUBROUTINE ConvertSpc_VVDry_to_Kg
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: convertspc_kg_to_vvdry
!
! !DESCRIPTION: Subroutine ConvertSpc\_Kg\_to\_VVDry converts the units of 
!  species concentrations from species mass per grid box [kg] to dry volume 
!  mixing ratio [mol species/mol dry air].  
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE ConvertSpc_Kg_to_VVDry( am_I_Root, State_Met, &
                                     State_Chm, RC   ) 
!
! !USES:
!
    USE GIGC_State_Met_Mod, ONLY : MetState
    USE GIGC_State_Chm_Mod, ONLY : ChmState
    USE Species_Mod,        ONLY : Species
!
! !INPUT PARAMETERS: 
!
    LOGICAL,        INTENT(IN)    :: am_I_Root   ! Are we on the root CPU?
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
!  21 Jul 2016 - E. Lundgren - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER                :: I, J, L, N
    CHARACTER(LEN=255)     :: MSG, LOC
    REAL(fp)               :: MW_G
    TYPE(Species), POINTER :: ThisSpc

    !====================================================================
    ! ConvertSpc_Kg_to_VVDry begins here!
    !====================================================================

    ! Assume success
    RC = GIGC_SUCCESS

    ! Verify correct initial units. If current units are unexpected,
    ! write error message and location to log, then pass failed RC
    ! to calling routine. 
    IF ( TRIM( State_Chm%Spc_Units ) /= 'kg' ) THEN
       MSG = 'Incorrect initial units: ' // TRIM( State_Chm%Spc_Units )
       LOC = 'Routine ConvertSpc_Kg_to_VVDry in unitconv_mod.F90'
       CALL GIGC_Error( MSG, RC, LOC )
       RETURN
    ENDIF

    !====================================================================
    !
    !  The conversion is as follows:
    !
    !                         1           g/mol dry air
    !   kg species(N)  * -------------- * ----------------    
    !                    kg dry air       g/mol species(N)      
    !
    !   = kg species(N) / dry air mass * MW dry air / MW species(N)
    !   
    !   = volume mixing ratio
    !
    ! Therefore, with:
    !
    !  AD(I,J,L) = grid box dry air mass [kg]
    !  AIRMW     = dry air molecular wt [g/mol]
    !  MW_G(N)   = species molecular wt [g/mol]
    !     
    ! the conversion is:
    ! 
    !  Species(I,J,L,N) [v/v]
    !
    !    = Species(I,J,L,N) [kg] * ( AIRMW / MW_G(N) ) / AD(I,J,L) 
    !                   
    !====================================================================

    !$OMP PARALLEL DO           &
    !$OMP DEFAULT( SHARED     ) &
    !$OMP PRIVATE( I, J, L, N, MW_G ) 
    DO N = 1, State_Chm%nSpecies

       ! Get info about this species from the species database
       ThisSpc    => State_Chm%SpcData(N)%Info
       MW_G      = State_Chm%SpcData(N)%Info%emMW_g

       DO L = 1, LLPAR
       DO J = 1, JJPAR
       DO I = 1, IIPAR
          State_Chm%Species(I,J,L,N) = State_Chm%Species(I,J,L,N)  &
                                       *  ( AIRMW / MW_G )         &
                                       / State_Met%AD(I,J,L)
       ENDDO
       ENDDO
       ENDDO

    ENDDO
    !$OMP END PARALLEL DO

    ! Update species units
    State_Chm%Spc_Units = 'v/v dry'

  END SUBROUTINE ConvertSpc_Kg_to_VVDry 
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
! !IROUTINE: convertspc_kgkgdry_to_kg
!
! !DESCRIPTION: Subroutine ConvertSpc\_KgKgDry\_to\_Kg converts the units of 
!  species concentrations from dry mass mixing ratio 
!  [kg species/kg dry air] to species mass per grid box [kg]. 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE ConvertSpc_KgKgDry_to_Kg( am_I_Root, State_Met, &
                                       State_Chm, RC   ) 
!
! !USES:
!
    USE GIGC_State_Met_Mod, ONLY : MetState
    USE GIGC_State_Chm_Mod, ONLY : ChmState
!
! !INPUT PARAMETERS: 
!
    LOGICAL,        INTENT(IN)    :: am_I_Root   ! Are we on the root CPU?
    TYPE(MetState), INTENT(IN)    :: State_Met   ! Meteorology state object
!
! !INPUT/OUTPUT PARAMETERS:
!
    ! Object containing species concentration
    TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry state object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
!
! !REMARKS:
!
! !REVISION HISTORY: 
!  21 Jul 2016 - E. Lundgren - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER :: I, J, L, N
    CHARACTER(LEN=255) :: MSG, LOC

    !====================================================================
    ! ConvertSpc_KgKgDry_to_Kg begins here!
    !====================================================================

    ! Assume success
    RC = GIGC_SUCCESS

    ! Verify correct initial units. If current units are unexpected,
    ! write error message and location to log, then pass failed RC
    ! to calling routine. 
    IF ( TRIM( State_Chm%Spc_Units ) /= 'kg/kg dry' ) THEN
       MSG = 'Incorrect initial units: ' // TRIM( State_Chm%Spc_Units )
       LOC = 'Routine ConvertSpc_KgKgDry_to_Kg in unitconv_mod.F90'
       CALL GIGC_Error( MSG, RC, LOC )
       RETURN
    ENDIF

    !====================================================================
    !
    !  The conversion is as follows:
    !
    !   kg species(N)            
    !   -----------  *  kg dry air       
    !   kg dry air                   
    !
    !   = mass mixing ratio * dry air mass  
    !   
    !   = kg species(N)
    !
    ! Therefore, with:
    !
    !  AD(I,J,L)   = grid box dry air mass [kg]
    !     
    ! the conversion is:
    ! 
    !  Species(I,J,L,N) [kg]
    !
    !    = Species(I,J,L,N) [kg/kg] * AD(I,J,L)
    !                   
    !====================================================================

    !$OMP PARALLEL DO           &
    !$OMP DEFAULT( SHARED     ) &
    !$OMP PRIVATE( I, J, L, N ) 
    DO N = 1, State_Chm%nSpecies
    DO L = 1, LLPAR
    DO J = 1, JJPAR
    DO I = 1, IIPAR
       State_Chm%Species(I,J,L,N) = State_Chm%Species(I,J,L,N) &
                                  * State_Met%AD(I,J,L)
    ENDDO
    ENDDO
    ENDDO
    ENDDO
    !$OMP END PARALLEL DO

    ! Update species units
    State_Chm%Spc_Units = 'kg'

  END SUBROUTINE ConvertSpc_KgKgDry_to_Kg
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: convertspc_kg_to_kgkgdry
!
! !DESCRIPTION: Subroutine ConvertSpc\_Kg\_to\_KgKgDry converts the units of 
!  species concentrations from species mass per grid box [kg] to dry mass 
!  mixing ratio [kg species/kg dry air].  
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE ConvertSpc_Kg_to_KgKgDry( am_I_Root, State_Met, &
                                       State_Chm, RC   ) 
!
! !USES:
!
    USE GIGC_State_Met_Mod, ONLY : MetState
    USE GIGC_State_Chm_Mod, ONLY : ChmState
!
! !INPUT PARAMETERS: 
!
    LOGICAL,        INTENT(IN)    :: am_I_Root   ! Are we on the root CPU?
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
!  21 Jul 2016 - E. Lundgren - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER :: I, J, L, N
    CHARACTER(LEN=255) :: MSG, LOC

    !====================================================================
    ! ConvertSpc_Kg_to_KgKgDry begins here!
    !====================================================================

    ! Assume success
    RC = GIGC_SUCCESS

    ! Verify correct initial units. If current units are unexpected,
    ! write error message and location to log, then pass failed RC
    ! to calling routine. 
    IF ( TRIM( State_Chm%Spc_Units ) /= 'kg' ) THEN
       MSG = 'Incorrect initial units: ' // TRIM( State_Chm%Spc_Units )
       LOC = 'Routine ConvertSpc_Kg_to_KgKgDry in unitconv_mod.F90'
       CALL GIGC_Error( MSG, RC, LOC )
       RETURN
    ENDIF

    !====================================================================
    !
    !  The conversion is as follows:
    !
    !                         1          
    !   kg species(N)  * --------------      
    !                      kg dry air              
    !
    !   = kg species(N) / dry air mass
    !   
    !   = mass mixing ratio
    !
    ! Therefore, with:
    !
    !  AD(I,J,L)    = grid box dry air mass [kg]
    !     
    ! the conversion is:
    ! 
    !  Species(I,J,L,N) [kg/kg]
    !
    !    = Species(I,J,L,N) [kg] / AD(I,J,L) 
    !                   
    !====================================================================

    !$OMP PARALLEL DO           &
    !$OMP DEFAULT( SHARED     ) &
    !$OMP PRIVATE( I, J, L, N ) 
    DO N = 1, State_Chm%nSpecies
    DO L = 1, LLPAR
    DO J = 1, JJPAR
    DO I = 1, IIPAR
       State_Chm%Species(I,J,L,N) = State_Chm%Species(I,J,L,N) &
                                    / State_Met%AD(I,J,L)
    ENDDO
    ENDDO
    ENDDO
    ENDDO
    !$OMP END PARALLEL DO

    ! Update species units
    State_Chm%Spc_units = 'kg/kg dry'

  END SUBROUTINE ConvertSpc_Kg_to_KgKgDry 
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: convert_kgkgtotal_to_kg
!
! !DESCRIPTION: Subroutine Convert\_KgKgTotal\_to\_Kg converts the units of 
!  tracer concentrations from moist mass mixing ratio [kg tracer/kg total 
!  (wet) air] to tracer mass per grid box [kg] for a single grid box. 
!  This routine is temporary during the unit transition of TOMAS to 
!  area-independence.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Convert_KgKgTotal_to_Kg( am_I_Root, I, J, L, Input_Opt,    &
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
!  08 Jan 2015 - E. Lundgren - Initial version
!  29 Sep 2015 - E. Lundgren - Now convert for a single grid box 
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER            :: N
    CHARACTER(LEN=255) :: MSG, LOC

    !=================================================================
    ! Convert_KgKgTotal_to_Kg begins here!
    !=================================================================

    ! Assume success
    RC        =  GIGC_SUCCESS

    !==============================================================
    !
    !  The conversion is as follows:
    !
    !   kg tracer(N)            
    !   -----------  *  kg total air       
    !   kg total air                   
    !
    !   = mass mixing ratio * total air mass  
    !   
    !   = kg tracer(N)
    !
    ! Therefore, with:
    !
    !  ADMOIST(I,J,L)   = grid box total air mass [kg] (wet)
    !     
    ! the conversion is:
    ! 
    !  TRACERS(I,J,L,N) [kg]
    !
    !    = TRACERS(I,J,L,N) [kg/kg] * ADMOIST(I,J,L)
    !                   
    !==============================================================

    !$OMP PARALLEL DO           &
    !$OMP DEFAULT( SHARED     ) &
    !$OMP PRIVATE( N ) 
    DO N = 1, Input_Opt%N_TRACERS
       State_Chm%TRACERS(I,J,L,N) = State_Chm%TRACERS(I,J,L,N)  &
                                  * State_Met%ADMOIST(I,J,L)
    ENDDO
    !$OMP END PARALLEL DO

  END SUBROUTINE Convert_KgKgTotal_to_Kg
!EOC

!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: convert_kg_to_kgkgtotal
!
! !DESCRIPTION: Subroutine Convert\_Kg\_to\_KgKgTotal converts the units of 
!  tracer concentrations from tracer mass per grid box [kg] to mass 
!  mixing ratio [kg tracer/kg total (wet) air] for a single grid box.
!  This routine is temporary during the unit transition of TOMAS to 
!  area-independence.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Convert_Kg_to_KgKgTotal( am_I_Root, I, J, L, Input_Opt,    &
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
!  08 Jan 2015 - E. Lundgren - Initial version
!  29 Sep 2015 - E. Lundgren - Now convert for a single grid box
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER :: N
    CHARACTER(LEN=255) :: MSG, LOC

    !====================================================================
    ! Convert_Kg_to_KgKgTotal begins here!
    !====================================================================

    ! Assume success
    RC        =  GIGC_SUCCESS

    !====================================================================
    !
    !  The conversion is as follows:
    !
    !                         1          
    !   kg tracer(N)  * --------------      
    !                    kg total air              
    !
    !   = kg tracer(N) / total air mass
    !   
    !   = mass mixing ratio
    !
    ! Therefore, with:
    !
    !  ADMOIST(I,J,L)    = grid box total air mass [kg] (wet)
    !     
    ! the conversion is:
    ! 
    !  TRACERS(I,J,L,N) [kg/kg]
    !
    !    = TRACERS(I,J,L,N) [kg] / ADMOIST(I,J,L) 
    !                   
    !====================================================================

    !$OMP PARALLEL DO           &
    !$OMP DEFAULT( SHARED     ) &
    !$OMP PRIVATE( N ) 
    DO N = 1, Input_Opt%N_TRACERS
       State_Chm%TRACERS(I,J,L,N) = State_Chm%TRACERS(I,J,L,N) &
                                  / State_Met%ADMOIST(I,J,L)
    ENDDO
    !$OMP END PARALLEL DO

  END SUBROUTINE Convert_Kg_to_KgKgTotal
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: convertspc_kgkgtotal_to_kg
!
! !DESCRIPTION: Subroutine ConvertSpc\_KgKgTotal\_to\_Kg converts the units of 
!  species concentrations from moist mass mixing ratio [kg species/kg total 
!  (wet) air] to species mass per grid box [kg] for a single grid box. 
!  This routine is temporary during the unit transition of TOMAS to 
!  area-independence.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE ConvertSpc_KgKgTotal_to_Kg( am_I_Root, I, J, L,    &
                                         State_Met, State_Chm, RC   ) 
!
! !USES:
!
    USE GIGC_State_Met_Mod, ONLY : MetState
    USE GIGC_State_Chm_Mod, ONLY : ChmState
!
! !INPUT PARAMETERS: 
!
    LOGICAL,        INTENT(IN)    :: am_I_Root   ! Are we on the root CPU?
    INTEGER,        INTENT(IN)    :: I, J, L     ! Grid box indexes
    TYPE(MetState), INTENT(IN)    :: State_Met   ! Meteorology state object
!
! !INPUT/OUTPUT PARAMETERS:
!
    ! Object containing species concentration
    TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry state object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
!
! !REMARKS:
!  This routine is temporary and is only used for local conversion of species 
!  concentrations for use in TOMAS within wetscav_mod routine WASHOUT. 
!  That routine is called within a parallel do loop and therefore units can
!  only be converted per grid box to avoid excessive computation time. Also,
!  State_Chm%Spc_Units cannot be changed within the parallel do loop without
!  causing problems. It is therefore left out of this routine.
!
! !REVISION HISTORY: 
!  21 Jul 2015 - E. Lundgren - Initial version - convert for a single grid box 
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER            :: N
    CHARACTER(LEN=255) :: MSG, LOC

    !=================================================================
    ! ConvertSpc_KgKgTotal_to_Kg begins here!
    !=================================================================

    ! Assume success
    RC        =  GIGC_SUCCESS

    !==============================================================
    !
    !  The conversion is as follows:
    !
    !   kg species(N)            
    !   -----------  *  kg total air       
    !   kg total air                   
    !
    !   = mass mixing ratio * total air mass  
    !   
    !   = kg species(N)
    !
    ! Therefore, with:
    !
    !  ADMOIST(I,J,L)   = grid box total air mass [kg] (wet)
    !     
    ! the conversion is:
    ! 
    !  Species(I,J,L,N) [kg]
    !
    !    = Species(I,J,L,N) [kg/kg] * ADMOIST(I,J,L)
    !                   
    !==============================================================

    !$OMP PARALLEL DO           &
    !$OMP DEFAULT( SHARED     ) &
    !$OMP PRIVATE( N ) 
    DO N = 1, State_Chm%nSpecies
       State_Chm%Species(I,J,L,N) = State_Chm%Species(I,J,L,N)  &
                                  * State_Met%ADMOIST(I,J,L)
    ENDDO
    !$OMP END PARALLEL DO

  END SUBROUTINE ConvertSpc_KgKgTotal_to_Kg
!EOC

!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: convertspc_kg_to_kgkgtotal
!
! !DESCRIPTION: Subroutine ConvertSpc\_Kg\_to\_KgKgTotal converts the units of 
!  species concentrations from species mass per grid box [kg] to mass 
!  mixing ratio [kg species/kg total (wet) air] for a single grid box.
!  This routine is temporary during the unit transition of TOMAS to 
!  area-independence.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE ConvertSpc_Kg_to_KgKgTotal( am_I_Root, I, J, L,    &
                                         State_Met, State_Chm, RC   ) 
!
! !USES:
!
 
    USE GIGC_State_Met_Mod, ONLY : MetState
    USE GIGC_State_Chm_Mod, ONLY : ChmState
!
! !INPUT PARAMETERS: 
!
    LOGICAL,        INTENT(IN)    :: am_I_Root   ! Are we on the root CPU?
    INTEGER,        INTENT(IN)    :: I, J, L     ! Grid box indexes
    TYPE(MetState), INTENT(IN)    :: State_Met   ! Meteorology state object
!
! !INPUT/OUTPUT PARAMETERS:
!
    ! Object containing species concentration
    TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry state object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
!
! !REMARKS:
!  This routine is temporary and is only used for local conversion of species 
!  concentrations for use in TOMAS within wetscav_mod routine WASHOUT. 
!  That routine is called within a parallel do loop and therefore units can
!  only be converted per grid box to avoid excessive computation time. Also,
!  State_Chm%Spc_Units cannot be changed within the parallel do loop without
!  causing problems. It is therefore left out of this routine.
!
! !REVISION HISTORY: 
!  21 Jul 2016 - E. Lundgren - Initial version - convert for a single grid box
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER :: N
    CHARACTER(LEN=255) :: MSG, LOC

    !====================================================================
    ! ConvertSpc_Kg_to_KgKgTotal begins here!
    !====================================================================

    ! Assume success
    RC        =  GIGC_SUCCESS

    !====================================================================
    !
    !  The conversion is as follows:
    !
    !                         1          
    !   kg species(N)  * --------------      
    !                    kg total air              
    !
    !   = kg species(N) / total air mass
    !   
    !   = mass mixing ratio
    !
    ! Therefore, with:
    !
    !  ADMOIST(I,J,L)    = grid box total air mass [kg] (wet)
    !     
    ! the conversion is:
    ! 
    !  Species(I,J,L,N) [kg/kg]
    !
    !    = Species(I,J,L,N) [kg] / ADMOIST(I,J,L) 
    !                   
    !====================================================================

    !$OMP PARALLEL DO           &
    !$OMP DEFAULT( SHARED     ) &
    !$OMP PRIVATE( N ) 
    DO N = 1, State_Chm%nSpecies
       State_Chm%Species(I,J,L,N) = State_Chm%Species(I,J,L,N) &
                                  / State_Met%ADMOIST(I,J,L)
    ENDDO
    !$OMP END PARALLEL DO

  END SUBROUTINE ConvertSpc_Kg_to_KgKgTotal
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: convert_kg_to_kgm2
!
! !DESCRIPTION: Subroutine Convert\_Kg\_to\_kgm2 converts the units of 
! mass [kg] to area density [kg/m2] for a single grid box.  This routine is
!  temporary during the unit transition of TOMAS to area-independence. 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Convert_Kg_to_Kgm2( am_I_Root, I, J, L, Input_Opt,      &
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
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER :: N
    CHARACTER(LEN=255) :: MSG, LOC

    !====================================================================
    ! Convert_Kg_to_Kgm2 begins here!
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

  END SUBROUTINE Convert_Kg_to_Kgm2
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: convert_kgm2_to_kg
!
! !DESCRIPTION: Subroutine Convert\_Kgm2\_to\_Kg converts the units of area 
!  density [kg/m2] to mass [kg] for a single grid box. This routine is
!  temporary during the unit transition of TOMAS to area-independence.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Convert_Kgm2_to_Kg( am_I_Root, I, J, L, Input_Opt,   &
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
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER            :: N
    CHARACTER(LEN=255) :: MSG, LOC

    !====================================================================
    ! Convert_Kgm2_to_Kg begins here!
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

  END SUBROUTINE Convert_Kgm2_to_Kg
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: convertspc_kg_to_kgm2
!
! !DESCRIPTION: Subroutine ConvertSpc\_Kg\_to\_kgm2 converts the units of 
! mass [kg] to area density [kg/m2] for a single grid box.  This routine is
!  temporary during the unit transition of TOMAS to area-independence. 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE ConvertSpc_Kg_to_Kgm2( am_I_Root, I, J, L,      &
                                    State_Met, State_Chm, RC        )
!
! !USES:
!
    USE GIGC_State_Chm_Mod, ONLY : ChmState
    USE GIGC_State_Met_Mod, ONLY : MetState
    USE PHYSCONSTANTS
!
! !INPUT PARAMETERS: 
!
    LOGICAL,        INTENT(IN)    :: am_I_Root     ! Are we on the root CPU?
    INTEGER,        INTENT(IN)    :: I, J, L       ! Grid box indexes
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
!  This routine is temporary and is only used for local conversion of species 
!  concentrations for use in TOMAS within wetscav_mod routine WASHOUT. 
!  That routine is called within a parallel do loop and therefore units can
!  only be converted per grid box to avoid excessive computation time. Also,
!  State_Chm%Spc_Units cannot be changed within the parallel do loop without
!  causing problems. It is therefore left out of this routine.
!
! !REVISION HISTORY: 
!  21 Jul 2016 - E. Lundgren - Initial version - convert single grid box only
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER :: N
    CHARACTER(LEN=255) :: MSG, LOC

    !====================================================================
    ! ConvertSpc_Kg_to_Kgm2 begins here!
    !====================================================================

    ! Assume success
    RC        =  GIGC_SUCCESS

    !$OMP PARALLEL DO        &
    !$OMP DEFAULT( SHARED  ) &
    !$OMP PRIVATE( N ) 
    DO N = 1, State_Chm%nSpecies
       State_Chm%Species(I,J,L,N) = State_Chm%Species(I,J,L,N)    &
                                    / State_Met%AREA_M2(I,J,1)
    ENDDO
    !$OMP END PARALLEL DO

  END SUBROUTINE ConvertSpc_Kg_to_Kgm2
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: convertspc_kgm2_to_kg
!
! !DESCRIPTION: Subroutine ConvertSpc\_Kgm2\_to\_Kg converts the units of area 
!  density [kg/m2] to mass [kg] for a single grid box. This routine is
!  temporary during the unit transition of TOMAS to area-independence.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE ConvertSpc_Kgm2_to_Kg( am_I_Root, I, J, L,   &
                                    State_Met, State_Chm, RC        )
!
! !USES:
!
    USE GIGC_State_Chm_Mod, ONLY : ChmState
    USE GIGC_State_Met_Mod, ONLY : MetState
    USE PHYSCONSTANTS
!
! !INPUT PARAMETERS: 
!
    LOGICAL,        INTENT(IN)    :: am_I_Root     ! Are we on the root CPU?
    INTEGER,        INTENT(IN)    :: I, J, L       ! Grid box indexes
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
!  This routine is temporary and is only used for local conversion of species 
!  concentrations for use in TOMAS within wetscav_mod routine WASHOUT. 
!  That routine is called within a parallel do loop and therefore units can
!  only be converted per grid box to avoid excessive computation time. Also,
!  State_Chm%Spc_Units cannot be changed within the parallel do loop without
!  causing problems. It is therefore left out of this routine.
!
! !REVISION HISTORY: 
!  21 Jul 2016 - E. Lundgren - Initial version - convert single grid box only
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER            :: N
    CHARACTER(LEN=255) :: MSG, LOC

    !====================================================================
    ! ConvertSpc_Kgm2_to_Kg begins here!
    !====================================================================

    ! Assume success
    RC        =  GIGC_SUCCESS

    !$OMP PARALLEL DO        &
    !$OMP DEFAULT( SHARED  ) &
    !$OMP PRIVATE( N ) 
    DO N = 1, State_Chm%nSpecies
       State_Chm%Species(I,J,L,N) = State_Chm%Species(I,J,L,N)    &
                                     * State_Met%AREA_M2(I,J,1)
    ENDDO
    !$OMP END PARALLEL DO

  END SUBROUTINE ConvertSpc_Kgm2_to_Kg
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: convertspc_kgkgdry_to_mnd
!
! !DESCRIPTION: Subroutine ConvertSpc\_KgKgDry\_to\_MND converts the units of 
!  species concentrations from dry mass mixing ratio [kg/kg dry air] to 
!  molecular number density (MND) [molecules/cm3].  
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE ConvertSpc_KgKgDry_to_MND( am_I_Root, State_Met, &
                                        State_Chm, RC )
!
! !USES:
!
    USE GIGC_State_Met_Mod, ONLY : MetState
    USE GIGC_State_Chm_Mod, ONLY : ChmState    
    USE Species_Mod,        ONLY : Species
!
! !INPUT PARAMETERS: 
!
    LOGICAL,        INTENT(IN)    :: am_I_Root   ! Are we on the root CPU?
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
!  21 Jul 2016 - E. Lundgren - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER                :: I, J, L, N
    CHARACTER(LEN=255)     :: MSG, LOC
    REAL(fp)               :: MolecRatio, MW_KG
    TYPE(Species), POINTER :: ThisSpc

    !====================================================================
    ! ConvertSpc_KgKgDry_to_MND begins here!
    !====================================================================

    ! Assume success
    RC        =  GIGC_SUCCESS

    ! Initialize pointer
    ThisSpc => NULL()

    ! Verify correct initial units. If current units are unexpected,
    ! write error message and location to log, then pass failed RC
    ! to calling routine. 
    IF ( TRIM( State_Chm%Spc_Units ) /= 'kg/kg dry' ) THEN
       MSG = 'Incorrect initial units: ' // TRIM( State_Chm%Spc_Units )
       LOC = 'Routine ConvertSpc_KgKgDry_to_MND in unitconv_mod.F90'
       CALL GIGC_Error( MSG, RC, LOC )
       RETURN
    ENDIF

    !====================================================================
    !
    ! The conversion is as follows:
    !
    !   kg species(N)    kg air     molec   mol species(N)     m3    
    !   -----------   * --------  * ----- * -------------  * -------   
    !   kg dry air         m3        mol        kg           1E6 cm3     
    !
    !   = mixing ratio * air density * Avogadro's # / MW * conversion factors
    !   
    !   = molecules per cm3
    !
    ! Therefore, with:
    !
    !  AVO             = Avogadro's #
    !  AIRDEN(I,J,L)   = grid box dry air density [kg/m3]
    !  MW_KG           = molecules species / kg species
    !     
    ! the conversion is:
    ! 
    !  Spcies(I,J,L,N) [molecules/cm3]
    !
    !    = Species(I,J,L,N) [kg/kg] * AIRDEN(I,J,L) * AVO / MW_KG / 1e6
    !
    ! NOTES: 
    !   (1) Also divide by mol C / mol species (equal to one for most spc)
    !   (2) Use AD/AIRVOL instead of AIRDEN to preserve legacy method
    !====================================================================

    !$OMP PARALLEL DO           &
    !$OMP DEFAULT( SHARED     ) &
    !$OMP PRIVATE( I, J, L, N, MolecRatio, MW_KG )  
    DO N = 1, State_Chm%nSpecies 

       ! Get info about this species from the species database
       ThisSpc    => State_Chm%SpcData(N)%Info
       MolecRatio = State_Chm%SpcData(N)%Info%MolecRatio ! mol C / mol spc
       MW_KG      = State_Chm%SpcData(N)%Info%emMW_g * 1.e-3_fp

       ! Loop over grid boxes
       DO L = 1, LLPAR
       DO J = 1, JJPAR
       DO I = 1, IIPAR
          State_Chm%Species(I,J,L,N) = State_Chm%Species(I,J,L,N)           &
                                     * ( State_Met%AD(I,J,L) /              &
                                         State_Met%AIRVOL(I,J,L) )          &
                                     * ( AVO / MW_KG )                      &
                                     / ( 1e+6_fp * MolecRatio )  

       ENDDO
       ENDDO
       ENDDO

       ! Free pointer
       ThisSpc => NULL()

    ENDDO
    !$OMP END PARALLEL DO

    ! Update species units
    State_Chm%Spc_Units = 'molec/cm3'

  END SUBROUTINE ConvertSpc_KgKgDry_to_MND
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: convertspc_mnd_to_kgkgdry
!
! !DESCRIPTION: Subroutine ConvertSpc\_MND\_to\_KgKgDry converts the units of 
!  species concentrations from molecular number density (MND) [molecules/cm3]
!  to dry mass mixing ratio [kg/kg dry air].
!  .  
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE ConvertSpc_MND_to_KgKgDry( am_I_Root, State_Met, &
                                        State_Chm, RC )
!
! !USES:
!
    USE GIGC_State_Met_Mod, ONLY : MetState
    USE GIGC_State_Chm_Mod, ONLY : ChmState
    USE Species_Mod,        ONLY : Species
!
! !INPUT PARAMETERS: 
!
    LOGICAL,        INTENT(IN)    :: am_I_Root   ! Are we on the root CPU?
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
!  21 Jul 2016 - E. Lundgren - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER                :: I, J, L, N
    CHARACTER(LEN=255)     :: MSG, LOC
    REAL(fp)               :: MolecRatio, MW_KG
    TYPE(Species), POINTER :: ThisSpc

    !====================================================================
    ! ConvertSpc_MND_to_KgKgDry begins here!
    !====================================================================

    ! Assume success
    RC        =  GIGC_SUCCESS

    ! Initialize pointer
    ThisSpc => NULL()

    ! Verify correct initial units. If current units are unexpected,
    ! write error message and location to log, then pass failed RC
    ! to calling routine. 
    IF ( TRIM( State_Chm%Spc_Units ) /= 'molec/cm3' ) THEN
       MSG = 'Incorrect initial units: ' // TRIM( State_Chm%Spc_Units )
       LOC = 'Routine ConvertSpc_MND_to_KgKgDry in unitconv_mod.F90'
       CALL GIGC_Error( MSG, RC, LOC )
       RETURN
    ENDIF

    !====================================================================
    !
    ! The conversion is as follows:
    !
    !   molec species(N)   mol     kg species(N)      m3      1E6 cm3     
    !   ---------------- * ----- * -------------- * ------  * -------
    !       cm3            molec   mol species(N)   kg air      m3    

    !
    !   = # density / Avogadro's # * MW / air density * conversion factors
    !   
    !   = kg species / kg dry air
    !
    ! Therefore, with:
    !
    !  AVO             = Avogadro's #
    !  AIRDEN(I,J,L)   = grid box dry air density [kg/m3]
    !  MW_KG           = molecules species / kg species
    !     
    ! the conversion is:
    ! 
    !  Spcies(I,J,L,N) [kg/kg dry air]
    !
    !    = Species(I,J,L,N) [molecules/cm3] * AIRDEN(I,J,L) * AVO / MW_KG / 1e6
    !
    ! NOTES: 
    !  (1) Also multiply by mol C / mol species (equal to one for most spc)
    !  (2) Use exact reverse of the mixing ratio -> # density conversion to
    !      avoid numerical noise differences
    !  (3) Use AD/AIRVOL instead of AIRDEN to preserve legacy method
    !
    !====================================================================

    !$OMP PARALLEL DO           &
    !$OMP DEFAULT( SHARED     ) &
    !$OMP PRIVATE( I, J, L, N, MolecRatio, MW_KG )  
    DO N = 1, State_Chm%nSpecies 

       ! Get info about this species from the species database
       ThisSpc => State_Chm%SpcData(N)%Info
       MolecRatio = State_Chm%SpcData(N)%Info%MolecRatio ! mol C / mol spc?
       MW_KG      = State_Chm%SpcData(N)%Info%emMW_g * 1.e-3_fp

       ! Loop over grid boxes
       DO L = 1, LLPAR
       DO J = 1, JJPAR
       DO I = 1, IIPAR
          State_Chm%Species(I,J,L,N) = State_Chm%Species(I,J,L,N)           &
                                     * ( 1e+6_fp * MolecRatio )             &
                                     / ( AVO / MW_KG )                      &
                                     / ( State_Met%AD(I,J,L) /              &
                                         State_Met%AIRVOL(I,J,L) )              
       ENDDO
       ENDDO
       ENDDO

       ! Free pointer
       ThisSpc => NULL()

    ENDDO
    !$OMP END PARALLEL DO

    ! Update species units
    State_Chm%Spc_Units = 'kg/kg dry'

  END SUBROUTINE ConvertSpc_MND_to_KgKgDry
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: convertspc_kg_to_mnd
!
! !DESCRIPTION: Subroutine ConvertSpc\_Kg\_to\_MND converts the units of 
!  species concentrations from mass per grid box [kg] to molecular 
!  number density (MND) [molecules/cm3].  
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE ConvertSpc_Kg_to_MND( am_I_Root, State_Met, State_Chm, RC )
!
! !USES:
!
    USE GIGC_State_Met_Mod, ONLY : MetState
    USE GIGC_State_Chm_Mod, ONLY : ChmState
    USE Species_Mod,        ONLY : Species
!
! !INPUT PARAMETERS: 
!
    LOGICAL,        INTENT(IN)    :: am_I_Root   ! Are we on the root CPU?
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
!  21 Jul 2016 - E. Lundgren - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER                :: I, J, L, N
    CHARACTER(LEN=255)     :: MSG, LOC
    REAL(fp)               :: MolecRatio, MW_KG
    TYPE(Species), POINTER :: ThisSpc

    !====================================================================
    ! ConvertSpc_Kg_to_MND begins here!
    !====================================================================

    ! Assume success
    RC        =  GIGC_SUCCESS

    ! Initialize pointer
    ThisSpc => NULL()

    ! Verify correct initial units. If current units are unexpected,
    ! write error message and location to log, then pass failed RC
    ! to calling routine. 
    IF ( TRIM( State_Chm%Spc_Units ) /= 'kg' ) THEN
       MSG = 'Incorrect initial units: ' // TRIM( State_Chm%Spc_Units )
       LOC = 'Routine ConvertSpc_KgKgDry_to_MND in unitconv_mod.F90'
       CALL GIGC_Error( MSG, RC, LOC )
       RETURN
    ENDIF

    !====================================================================
    !
    ! The conversion is as follows:
    !
    !                   molec   mol species(N)      1       m3    
    !   kg species(N) * ----- * -------------- * ------ * -------   
    !                    mol    kg species(N)    box m3   1E6 cm3     
    !
    !   = species mass * Avogadro's # / MW / box volume * conversion factors
    !   
    !   = molecules per cm3
    !
    ! Therefore, with:
    !
    !  AVO             = Avogadro's #
    !  AIRVOL(I,J,L)   = grid box volume [m3]
    !  MW_KG           = molecules species / kg species
    !     
    ! the conversion is:
    ! 
    !  Spcies(I,J,L,N) [molecules/cm3]
    !
    !    = Species(I,J,L,N) [kg] / AIRVOL(I,J,L) * AVO / MW_KG / 1e6
    !
    ! NOTE: Also divide by mol C / mol species (equal to one for most spc)
    !                   
    !====================================================================

    !$OMP PARALLEL DO           &
    !$OMP DEFAULT( SHARED     ) &
    !$OMP PRIVATE( I, J, L, N, MolecRatio, MW_KG )  
    DO N = 1, State_Chm%nSpecies 

       ! Get info about this species from the species database
       ThisSpc    => State_Chm%SpcData(N)%Info
       MolecRatio = State_Chm%SpcData(N)%Info%MolecRatio ! mol C / mol spc
       MW_KG      = State_Chm%SpcData(N)%Info%emMW_g * 1.e-3_fp

       ! Loop over grid boxes
       DO L = 1, LLPAR
       DO J = 1, JJPAR
       DO I = 1, IIPAR
          State_Chm%Species(I,J,L,N) = State_Chm%Species(I,J,L,N)           &
                                     / State_Met%AIRVOL(I,J,L)              &
                                     * ( AVO / MW_KG )                      &
                                     / ( 1e+6_fp * MolecRatio )   
       ENDDO
       ENDDO
       ENDDO

       ! Free pointer
       ThisSpc => NULL()

    ENDDO
    !$OMP END PARALLEL DO

    ! Update species units
    State_Chm%Spc_Units = 'molec/cm3'

  END SUBROUTINE ConvertSpc_Kg_to_MND
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: convertspc_mnd_to_kg
!
! !DESCRIPTION: Subroutine ConvertSpc\_MND\_to\_Kg converts the units of 
!  species concentrations from molecular number density (MND) 
!  [molecules/cm3] to mass per grid box [kg].
!  .  
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE ConvertSpc_MND_to_Kg( am_I_Root, State_Met, State_Chm, RC )
!
! !USES:
!
    USE GIGC_State_Met_Mod, ONLY : MetState
    USE GIGC_State_Chm_Mod, ONLY : ChmState
    USE Species_Mod,        ONLY : Species
!
! !INPUT PARAMETERS: 
!
    LOGICAL,        INTENT(IN)    :: am_I_Root   ! Are we on the root CPU?
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
!  21 Jul 2016 - E. Lundgren - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER                :: I, J, L, N
    CHARACTER(LEN=255)     :: MSG, LOC
    REAL(fp)               :: MolecRatio, MW_KG
    TYPE(Species), POINTER :: ThisSpc

    !====================================================================
    ! ConvertSpc_MND_to_Kg begins here!
    !====================================================================

    ! Assume success
    RC        =  GIGC_SUCCESS

    ! Initialize pointer
    ThisSpc => NULL()

    ! Verify correct initial units. If current units are unexpected,
    ! write error message and location to log, then pass failed RC
    ! to calling routine. 
    IF ( TRIM( State_Chm%Spc_Units ) /= 'molec/cm3' ) THEN
       MSG = 'Incorrect initial units: ' // TRIM( State_Chm%Spc_Units )
       LOC = 'Routine ConvertSpc_MND_to_KgKgDry in unitconv_mod.F90'
       CALL GIGC_Error( MSG, RC, LOC )
       RETURN
    ENDIF

    !====================================================================
    !
    ! The conversion is as follows:
    !
    !   molec species(N)   mol     kg species(N)    box m3    1E6 cm3     
    !   ---------------- * ----- * -------------- * ------  * -------
    !       cm3            molec   mol species(N)      1        m3    

    !
    !   = # density / Avogadro's # * MW / box volume * conversion factors
    !   
    !   = kg species
    !
    ! Therefore, with:
    !
    !  AVO             = Avogadro's #
    !  AIRVOL(I,J,L)   = grid box volume [m3]
    !  MW_KG           = molecules species / kg species
    !     
    ! the conversion is:
    ! 
    !  Species(I,J,L,N) [kg]
    !
    !    = Species(I,J,L,N) [molecules/cm3] * AIRVOL(I,J,L) * AVO / MW_KG / 1e6
    !
    ! NOTES: 
    !  (1) Also multiply by mol C / mol species (equal to one for most spc)
    !  (2) Use exact reverse of the species mass -> # density conversion to
    !      avoid numerical noise differences
    !
    !====================================================================

    !$OMP PARALLEL DO           &
    !$OMP DEFAULT( SHARED     ) &
    !$OMP PRIVATE( I, J, L, N, MolecRatio, MW_KG )  
    DO N = 1, State_Chm%nSpecies 

       ! Get info about this species from the species database
       ThisSpc => State_Chm%SpcData(N)%Info
       MolecRatio = State_Chm%SpcData(N)%Info%MolecRatio ! mol C / mol spc?
       MW_KG      = State_Chm%SpcData(N)%Info%emMW_g * 1.e-3_fp

       ! Loop over grid boxes
       DO L = 1, LLPAR
       DO J = 1, JJPAR
       DO I = 1, IIPAR
          State_Chm%Species(I,J,L,N) = State_Chm%Species(I,J,L,N)           &
                                     * ( 1e+6_fp * MolecRatio )             &
                                     / ( AVO / MW_KG )                      &
                                     / State_Met%AIRVOL(I,J,L)              
       ENDDO
       ENDDO
       ENDDO

       ! Free pointer
       ThisSpc => NULL()

    ENDDO
    !$OMP END PARALLEL DO

    ! Update species units
    State_Chm%Spc_Units = 'kg'

  END SUBROUTINE ConvertSpc_MND_to_Kg
!EOC
END MODULE UnitConv_Mod
