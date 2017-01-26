!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: unitconv_mod.F90
!
! !DESCRIPTION: Module UNITCONV\_MOD contains routines which are used to 
!  convert the units of species concentrations between mass 
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
  USE ErrCode_Mod
  USE ERROR_MOD
  USE PHYSCONSTANTS
  USE PRECISION_MOD         ! GEOS-Chem Flexible Precision (fp)
                    
  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! Wrapper routine 
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  PUBLIC :: Convert_Units

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! KG/KG DRY <-> V/V DRY
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! kg/kg dry air <-> v/v dry air
  ! Used in DO_TEND in mixing
  PUBLIC  :: ConvertSpc_KgKgDry_to_VVDry
  PUBLIC  :: ConvertSpc_VVDry_to_KgKgDry

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! KG/KG DRY <-> KG/M2
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! kg/kg dry air <-> kg/m2
  ! Used for wet deposition, DO_TEND in mixing,
  ! and around AIRQNT and SET_H2O_TRAC in main
  PUBLIC  :: ConvertSpc_KgKgDry_to_Kgm2
  PUBLIC  :: ConvertSpc_kgm2_to_KgKgDry

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! KG/KG DRY <-> MOLEC/CM3
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  PUBLIC  :: ConvertSpc_KgKgDry_to_MND
  PUBLIC  :: ConvertSpc_MND_to_KgKgDry

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! AREA-DEPENDENT (temporary routines)
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  ! v/v dry air <-> kg/grid box
  ! Temporarily replaces legacy CONVERT_UNITS
  ! Used in strat_chem_mod and sulfate_mod
  PUBLIC  :: ConvertSpc_VVDry_to_Kg
  PUBLIC  :: ConvertSpc_Kg_to_VVDry

  ! kg/kg dry air <-> kg/grid box
  ! Used in aerosol_mod, tomas_mod, emissions_mod,
  ! strat_chem_mod, exchange_mod, rrtmg_rad_transfer_mod,
  ! chemistry_mod, sulfate_mod, and carbon_mod
  ! This is since RRTMG, TOMAS, exchange_mod, chemistry,
  ! and EMISSMERCURY are still in [kg]
  PUBLIC  :: ConvertSpc_KgKgDry_to_Kg
  PUBLIC  :: ConvertSpc_Kg_to_KgKgDry

  ! molec/cm3 dry air <-> kg/gridbox
  PUBLIC  :: ConvertSpc_MND_to_Kg
  PUBLIC  :: ConvertSpc_Kg_to_MND

  ! kg/kg dry air <-> kg/grid box (single box only)
  ! Used for TOMAS compatibility in WASHOUT
  PUBLIC  :: ConvertBox_KgKgDry_to_Kg
  PUBLIC  :: ConvertBox_Kg_to_KgKgDry

  ! kg <-> kg/m2 (single box only)
  ! Used for TOMAS compatibility in WASHOUT within wetscav_mod
  PUBLIC  :: ConvertBox_Kgm2_to_Kg
  PUBLIC  :: ConvertBox_Kg_to_Kgm2
!
! !REMARKS:
!  The routines in this module are used to convert the units of
!  species concentrations in various GEOS-Chem routines.
!
! !REVISION HISTORY:
!  23 Jun 2015 - E. Lundgren - Initial version
!  13 Aug 2015 - E. Lundgren - Add tracer unit error handling
!  29 Sep 2015 - E. Lundgren - Adjust some of the unit conversions to/from kg
!                              to be for a single grid box for TOMAS
!  21 Jul 2016 - E. Lundgren - Add species unit conversion routines
!  26 Jul 2016 - E. Lundgren - Remove unused conversions and use "Box" in 
!                              TOMAS-specific unit conversions
!  23 Aug 2016 - M. Sulprizio- Remove tracer unit conversion routines, only
!                              species unit conversion routines remain
!EOP
!------------------------------------------------------------------------------
!BOC
CONTAINS


!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Convert_Units
!
! !DESCRIPTION: Subroutine Convert\_Units is a wrapper function to convert
!  the tracer input array to a desired unit. This routine is currently only
!  used by tendences\_mod.F90 and incomplete. 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Convert_Units ( am_I_Root, Input_Opt, State_Met, &
                             State_Chm, OutUnit,   RC, InUnit ) 
!
! !USES:
!
    USE Input_Opt_Mod, ONLY : OptInput
    USE State_Met_Mod, ONLY : MetState
    USE State_Chm_Mod, ONLY : ChmState
!
! !INPUT PARAMETERS: 
!
    LOGICAL,          INTENT(IN)            :: am_I_Root   ! Are we on the root CPU?
    TYPE(OptInput),   INTENT(IN)            :: Input_Opt   ! Input Options object
    TYPE(MetState),   INTENT(IN)            :: State_Met   ! Meteorology state object
    CHARACTER(LEN=*), INTENT(IN)            :: OutUnit     ! Desired output unit
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ChmState),   INTENT(INOUT)         :: State_Chm   ! Chemistry state object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(OUT)           :: RC          ! Success or failure?
    CHARACTER(LEN=*), INTENT(OUT), OPTIONAL :: InUnit      ! Units of input data 
!
! !REMARKS:
!
! !REVISION HISTORY: 
!  14 Apr 2016 - C. Keller   - Initial version
!  10 Oct 2016 - C. Keller   - Update to v11-01h 
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    CHARACTER(LEN=255) :: MSG, LOC
    LOGICAL            :: DONE

    !====================================================================
    ! Convert_Units begins here!
    !====================================================================

    ! Assume success
    RC =  GC_SUCCESS

    ! Archive units of input data
    IF ( PRESENT(InUnit) ) InUnit = State_Chm%Spc_Units 

    ! Initialize flag to indicate if we are done or not
    DONE = .FALSE.

    ! Leave here if data is already on desired units
    IF ( TRIM(OutUnit) == TRIM(State_Chm%Spc_Units) ) DONE = .TRUE. 

    !====================================================================
    ! Convert to kg/kg dry
    !====================================================================
    IF ( .NOT. DONE .AND. TRIM(OutUnit) == 'kg/kg dry' ) THEN 

       IF ( TRIM(State_Chm%Spc_Units) == 'v/v dry' ) THEN
          CALL ConvertSpc_VVDry_to_KgKgDry( am_I_Root, State_Chm, RC ) 
          IF ( RC /= GC_SUCCESS ) RETURN
          DONE = .TRUE.

       ELSEIF ( TRIM(State_Chm%Spc_Units) == 'kg' ) THEN
          CALL ConvertSpc_Kg_to_KgKgDry( am_I_Root, State_Met, State_Chm, RC ) 
          IF ( RC /= GC_SUCCESS ) RETURN
          DONE = .TRUE.

       ELSEIF ( TRIM(State_Chm%Spc_Units) == 'kg/m2' ) THEN
          CALL ConvertSpc_Kgm2_to_KgKgDry( am_I_Root, State_Met, State_Chm, RC ) 
          IF ( RC /= GC_SUCCESS ) RETURN
          DONE = .TRUE.
       ENDIF
    ENDIF

    !====================================================================
    ! Convert to v/v dry
    !====================================================================
    IF ( .NOT. DONE .AND. TRIM(OutUnit) == 'v/v dry' ) THEN 

       IF ( TRIM(State_Chm%Spc_Units) == 'kg/kg dry' ) THEN
          CALL ConvertSpc_KgKgDry_to_VVDry( am_I_Root, State_Chm, RC ) 
          IF ( RC /= GC_SUCCESS ) RETURN
          DONE = .TRUE.
       ENDIF
    ENDIF

    !====================================================================
    ! Convert to kg
    !====================================================================
    IF ( .NOT. DONE .AND. TRIM(OutUnit) == 'kg' ) THEN 

       IF ( TRIM(State_Chm%Spc_Units) == 'kg/kg dry' ) THEN
          CALL ConvertSpc_KgKgDry_to_Kg( am_I_Root, State_Met, State_Chm, RC ) 
          IF ( RC /= GC_SUCCESS ) RETURN
          DONE = .TRUE.
       ENDIF
    ENDIF

    !====================================================================
    ! Convert to kg/m2
    !====================================================================
    IF ( .NOT. DONE .AND. TRIM(OutUnit) == 'kg/m2' ) THEN 

       IF ( TRIM(State_Chm%Spc_Units) == 'kg/kg dry' ) THEN
          CALL ConvertSpc_KgKgDry_to_Kgm2( am_I_Root, State_Met, State_Chm, RC ) 
          IF ( RC /= GC_SUCCESS ) RETURN
          DONE = .TRUE.
       ENDIF
    ENDIF

    !====================================================================
    ! Error check
    !====================================================================
    IF ( .NOT. DONE ) THEN 
       MSG = 'Cannot convert from '//TRIM(InUnit)//' to '//TRIM(OutUnit)
       LOC = 'Routine Convert_Units in unitconv_mod.F90'
       CALL GC_Error( MSG, RC, LOC )
       RETURN
    ENDIF

  END SUBROUTINE Convert_Units
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ConvertSpc_kgkgdry_to_vvdry
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
    USE State_Chm_Mod,      ONLY : ChmState
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
    REAL(fp)               :: MW_g
    CHARACTER(LEN=255)     :: MSG, LOC
    
    !====================================================================
    ! ConvertSpc_KgKgDry_to_VVDry begins here!
    !====================================================================

    ! Assume success
    RC =  GC_SUCCESS

    ! Verify correct initial units. If current units are unexpected,
    ! write error message and location to log, then pass failed RC
    ! to calling routine. 
    IF ( TRIM( State_Chm%Spc_Units ) /= 'kg/kg dry' ) THEN
       MSG = 'Incorrect initial units: ' // TRIM( State_Chm%Spc_Units )
       LOC = 'Routine ConvertSpc_KgKgDry_to_VVDry in unitconv_mod.F90'
       CALL GC_Error( MSG, RC, LOC )
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
    
    ! Loop over all species
    !$OMP PARALLEL DO                 &
    !$OMP DEFAULT( SHARED           ) &
    !$OMP PRIVATE( I, J, L, N, MW_g ) 
    DO N = 1, State_Chm%nSpecies

       ! (Emitted) molecular weight for the species [g]
       ! NOTE: Non-advected species will have a MW of -1, which will
       ! make the species concentration negative.  This can be used to
       ! flag that the species should not be used.  The inverse unit
       ! conversion will flip the sign back to positive (ewl, bmy, 8/4/16)
       MW_g = State_Chm%SpcData(N)%Info%emMW_g
    
       ! Loop over grid boxes and do unit conversion
       DO L = 1, LLPAR
       DO J = 1, JJPAR
       DO I = 1, IIPAR
          State_Chm%Species(I,J,L,N) = State_Chm%Species(I,J,L,N)  &
                                     * ( AIRMW / MW_G )
       ENDDO
       ENDDO
       ENDDO
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
! !IROUTINE: ConvertSpc_vvdry_to_kgkgdry
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
    USE State_Chm_Mod,      ONLY : ChmState
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
    REAL(fp)               :: MW_g
    CHARACTER(LEN=255)     :: MSG, LOC

    !====================================================================
    ! ConvertSpc_VVDry_to_KgKgDry begins here!
    !=================================================================

    ! Assume success
    RC = GC_SUCCESS

    ! Verify correct initial units. If current units are unexpected,
    ! write error message and location to log, then pass failed RC
    ! to calling routine. 
    IF ( TRIM( State_Chm%Spc_Units ) /= 'v/v dry' ) THEN
       MSG = 'Incorrect initial units: ' // TRIM( State_Chm%Spc_Units )
       LOC = 'Routine ConvertSpc_VVDry_to_KgKgDry in unitconv_mod.F90'
       CALL GC_Error( MSG, RC, LOC )
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

    ! Loop over all species
    !$OMP PARALLEL DO                 &
    !$OMP DEFAULT( SHARED           ) &
    !$OMP PRIVATE( I, J, L, N, MW_g ) 
    DO N = 1, State_Chm%nSpecies

       ! (Emitted) molecular weight for the species [g]
       ! NOTE: Non-advected species will have a MW of -1, which will
       ! make the species concentration negative.  This can be used to
       ! flag that the species should not be used.  The inverse unit
       ! conversion will flip the sign back to positive (ewl, bmy, 8/4/16)
       MW_g = State_Chm%SpcData(N)%Info%emMW_g

       ! Loop over grid boxes and do the unit conversion
       DO L = 1, LLPAR
       DO J = 1, JJPAR
       DO I = 1, IIPAR
         State_Chm%Species(I,J,L,N) = State_Chm%Species(I,J,L,N)   &
                                    / ( AIRMW / MW_G )
       ENDDO
       ENDDO
       ENDDO

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
! !IROUTINE: ConvertSpc_kgkgdry_to_kgm2
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
    USE State_Chm_Mod,      ONLY : ChmState
    USE State_Met_Mod,      ONLY : MetState
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
!  16 Sep 2016 - E. Lundgren - Replace DELP and SPHU with DELP_DRY
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER            :: I, J, L, N
    CHARACTER(LEN=255) :: MSG, LOC

    !====================================================================
    ! ConvertSpc_KgKgDry_to_Kgm2 begins here!
    !====================================================================

    ! Assume success
    RC = GC_SUCCESS

    ! Verify correct initial units. If current units are unexpected,
    ! write error message and location to log, then pass failed RC
    ! to calling routine. 
    IF ( TRIM( State_Chm%Spc_Units ) /= 'kg/kg dry' ) THEN
       MSG = 'Incorrect initial units: ' // TRIM( State_Chm%Spc_Units )
       LOC = 'Routine ConvertSpc_KgKgDry_to_Kgm2 in unitconv_mod.F90'
       CALL GC_Error( MSG, RC, LOC )
       RETURN
    ENDIF

    !====================================================================
    !
    !  The conversion is as follows:
    !
    !   kg species     Delta dry P [hPa]   100 [Pa]
    !   -----------  * ----------------- * --------
    !   kg dry air     g [m/s2]            [hPa]   
    !
    !   = kg species / m2
    !
    ! where:
    !
    !  Delta dry P = dry pressure difference across level as derived
    !                from the dry surface pressure with A and B params
    !  g = acceleration due to gravity
    !  kg dry air / kg total air  = 1 - specific humidity
    !     
    !====================================================================
    !$OMP PARALLEL DO            &
    !$OMP DEFAULT( SHARED      ) &
    !$OMP PRIVATE( I, J, L, N  ) 
    DO N = 1, State_Chm%nSpecies
    DO L = 1, LLPAR
    DO J = 1, JJPAR
    DO I = 1, IIPAR
       State_Chm%Species(I,J,L,N) = State_Chm%Species(I,J,L,N)          &
                                    * ( g0_100                          &
                                    * State_Met%DELP_DRY(I,J,L) ) 
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
! !IROUTINE: ConvertSpc_kgm2_to_kgkgdry
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
    USE State_Chm_Mod,      ONLY : ChmState
    USE State_Met_Mod,      ONLY : MetState
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
!  16 Sep 2016 - E. Lundgren - Replace DELP and SPHU with DELP_DRY
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER            :: I, J, L, N
    CHARACTER(LEN=255) :: MSG, LOC

    !====================================================================
    ! ConvertSpc_Kgm2_to_KgKgDry begins here!
    !====================================================================

    ! Assume success
    RC = GC_SUCCESS

    ! Verify correct initial units. If current units are unexpected,
    ! write error message and location to log, then pass failed RC
    ! to calling routine. 
    IF ( TRIM( State_Chm%Spc_Units ) /= 'kg/m2' ) THEN
       MSG = 'Incorrect initial units: ' // TRIM( State_Chm%Spc_Units )
       LOC = 'Routine ConvertSpc_Kgm2_to_KgKgDry in unitconv_mod.F90'
       CALL GC_Error( MSG, RC, LOC )
       RETURN
    ENDIF

    !====================================================================
    !
    !  The conversion is as follows:
    !
    !   kg species(N)   g [m/s2]            [hPa]
    !   -----------  * ----------------- * --------
    !        m2        Delta dry P [hPa]   100 [Pa]
    !
    !   = kg species(N) / kg dry air
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
    DO N = 1, State_Chm%nSpecies
    DO L = 1, LLPAR
    DO J = 1, JJPAR
    DO I = 1, IIPAR
       State_Chm%Species(I,J,L,N) = State_Chm%Species(I,J,L,N)          &
                                    * ( 1.0e+0_fp                       &      
                                    / ( g0_100                          &
                                    * State_Met%DELP_DRY(I,J,L) ) )
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
! !IROUTINE: ConvertSpc_kgkgdry_to_mnd
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
    USE State_Chm_Mod,      ONLY : ChmState    
    USE State_Met_Mod,      ONLY : MetState
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
    REAL(fp)           :: MolecRatio, MW_kg
    CHARACTER(LEN=255) :: MSG, LOC

    !====================================================================
    ! ConvertSpc_KgKgDry_to_MND begins here!
    !====================================================================

    ! Assume success
    RC = GC_SUCCESS

    ! Verify correct initial units. If current units are unexpected,
    ! write error message and location to log, then pass failed RC
    ! to calling routine. 
    IF ( TRIM( State_Chm%Spc_Units ) /= 'kg/kg dry' ) THEN
       MSG = 'Incorrect initial units: ' // TRIM( State_Chm%Spc_Units )
       LOC = 'Routine ConvertSpc_KgKgDry_to_MND in unitconv_mod.F90'
       CALL GC_Error( MSG, RC, LOC )
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

    ! Loop over all species
    !$OMP PARALLEL DO                              &
    !$OMP DEFAULT( SHARED                        ) &
    !$OMP PRIVATE( I, J, L, N, MolecRatio, MW_kg )  
    DO N = 1, State_Chm%nSpecies 

       ! Moles C / moles species
       MolecRatio = State_Chm%SpcData(N)%Info%MolecRatio

       ! (Emitted) molecular weight for the species [kg]
       ! NOTE: Non-advected species will have a MW of -1, which will
       ! make the species concentration negative.  This can be used to
       ! flag that the species should not be used.  The inverse unit
       ! conversion will flip the sign back to positive (ewl, bmy, 8/4/16)
       MW_kg      = State_Chm%SpcData(N)%Info%emMW_g * 1.e-3_fp

       ! Loop over grid boxes and do the unit conversion
       DO L = 1, LLPAR
       DO J = 1, JJPAR
       DO I = 1, IIPAR
          State_Chm%Species(I,J,L,N) = State_Chm%Species(I,J,L,N)           &
                                     * State_Met%AIRDEN(I,J,L)              &
                                     * ( AVO / MW_kg )                      &
                                     / ( 1e+6_fp * MolecRatio )  

       ENDDO
       ENDDO
       ENDDO

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
! !IROUTINE: ConvertSpc_mnd_to_kgkgdry
!
! !DESCRIPTION: Subroutine ConvertSpc\_MND\_to\_KgKgDry converts the units of 
!  species concentrations from molecular number density (MND) [molecules/cm3]
!  to dry mass mixing ratio [kg/kg dry air].
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE ConvertSpc_MND_to_KgKgDry( am_I_Root, State_Met, &
                                        State_Chm, RC )
!
! !USES:
!
    USE State_Chm_Mod,      ONLY : ChmState
    USE State_Met_Mod,      ONLY : MetState
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
    REAL(fp)           :: MolecRatio, MW_kg
    CHARACTER(LEN=255) :: MSG, LOC

    !====================================================================
    ! ConvertSpc_MND_to_KgKgDry begins here!
    !====================================================================

    ! Assume success
    RC = GC_SUCCESS

    ! Verify correct initial units. If current units are unexpected,
    ! write error message and location to log, then pass failed RC
    ! to calling routine. 
    IF ( TRIM( State_Chm%Spc_Units ) /= 'molec/cm3' ) THEN
       MSG = 'Incorrect initial units: ' // TRIM( State_Chm%Spc_Units )
       LOC = 'Routine ConvertSpc_MND_to_KgKgDry in unitconv_mod.F90'
       CALL GC_Error( MSG, RC, LOC )
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

    ! Loop over species
    !$OMP PARALLEL DO                              &
    !$OMP DEFAULT( SHARED                        ) &
    !$OMP PRIVATE( I, J, L, N, MolecRatio, MW_kg )  
    DO N = 1, State_Chm%nSpecies 

       ! Moles C / moles species
       MolecRatio = State_Chm%SpcData(N)%Info%MolecRatio

       ! (Emitted) molecular weight for the species [kg]
       ! NOTE: Non-advected species will have a MW of -1, which will
       ! make the species concentration negative.  This can be used to
       ! flag that the species should not be used.  The inverse unit
       ! conversion will flip the sign back to positive (ewl, bmy, 8/4/16)
       MW_kg      = State_Chm%SpcData(N)%Info%emMW_g * 1.e-3_fp

       ! Loop over grid boxes and do the unit conversion
       DO L = 1, LLPAR
       DO J = 1, JJPAR
       DO I = 1, IIPAR
          State_Chm%Species(I,J,L,N) = State_Chm%Species(I,J,L,N)           &
                                     * ( 1e+6_fp * MolecRatio )             &
                                     / ( AVO / MW_kg )                      &
                                     / State_Met%AIRDEN(I,J,L)
       ENDDO
       ENDDO
       ENDDO

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
! !IROUTINE: ConvertSpc_vvdry_to_kg
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
    USE State_Chm_Mod,      ONLY : ChmState
    USE State_Met_Mod,      ONLY : MetState
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
    INTEGER            :: I, J, L, N
    REAL(fp)           :: MW_g
    CHARACTER(LEN=255) :: MSG, LOC

    !====================================================================
    ! ConvertSpc_KgKgDry_to_Kg begins here!
    !====================================================================

    ! Assume success
    RC = GC_SUCCESS

    ! Verify correct initial units. If current units are unexpected,
    ! write error message and location to log, then pass failed RC
    ! to calling routine. 
    IF ( TRIM( State_Chm%Spc_Units ) /= 'v/v dry' ) THEN
       MSG = 'Incorrect initial units: ' // TRIM( State_Chm%Spc_Units )
       LOC = 'Routine ConvertSpc_VVDry_to_Kg in unitconv_mod.F90'
       CALL GC_Error( MSG, RC, LOC )
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

    ! Loop over all species
    !$OMP PARALLEL DO                 &
    !$OMP DEFAULT( SHARED           ) &
    !$OMP PRIVATE( I, J, L, N, MW_g ) 
    DO N = 1, State_Chm%nSpecies

       ! (Emitted) molecular weight for the species [g]
       ! NOTE: Non-advected species will have a MW of -1, which will
       ! make the species concentration negative.  This can be used to
       ! flag that the species should not be used.  The inverse unit
       ! conversion will flip the sign back to positive (ewl, bmy, 8/4/16)
       MW_g = State_Chm%SpcData(N)%Info%emMW_g

       ! Loop over grid boxes and do the unit conversion
       DO L = 1, LLPAR
       DO J = 1, JJPAR
       DO I = 1, IIPAR
          State_Chm%Species(I,J,L,N) = State_Chm%Species(I,J,L,N)  &
                                       * State_Met%AD(I,J,L)       &
                                       / ( AIRMW / MW_g )
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
! !IROUTINE: ConvertSpc_kg_to_vvdry
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
    USE State_Chm_Mod,      ONLY : ChmState
    USE State_Met_Mod,      ONLY : MetState
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
    INTEGER            :: I, J, L, N
    REAL(fp)           :: MW_g
    CHARACTER(LEN=255) :: MSG, LOC

    !====================================================================
    ! ConvertSpc_Kg_to_VVDry begins here!
    !====================================================================

    ! Assume success
    RC = GC_SUCCESS

    ! Verify correct initial units. If current units are unexpected,
    ! write error message and location to log, then pass failed RC
    ! to calling routine. 
    IF ( TRIM( State_Chm%Spc_Units ) /= 'kg' ) THEN
       MSG = 'Incorrect initial units: ' // TRIM( State_Chm%Spc_Units )
       LOC = 'Routine ConvertSpc_Kg_to_VVDry in unitconv_mod.F90'
       CALL GC_Error( MSG, RC, LOC )
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

    ! Loop over all species
    !$OMP PARALLEL DO                 &
    !$OMP DEFAULT( SHARED           ) &
    !$OMP PRIVATE( I, J, L, N, MW_g ) 
    DO N = 1, State_Chm%nSpecies

       ! (Emitted) molecular weight for the species [g]
       ! NOTE: Non-advected species will have a MW of -1, which will
       ! make the species concentration negative.  This can be used to
       ! flag that the species should not be used.  The inverse unit
       ! conversion will flip the sign back to positive (ewl, bmy, 8/4/16)
       MW_g = State_Chm%SpcData(N)%Info%emMW_g

       ! Loop over grid boxes and do the unit conversion
       DO L = 1, LLPAR
       DO J = 1, JJPAR
       DO I = 1, IIPAR
          State_Chm%Species(I,J,L,N) = State_Chm%Species(I,J,L,N)  &
                                       *  ( AIRMW / MW_g )         &
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
! !IROUTINE: ConvertSpc_kgkgdry_to_kg
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
    USE State_Chm_Mod,      ONLY : ChmState
    USE State_Met_Mod,      ONLY : MetState
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
    INTEGER            :: I, J, L, N
    CHARACTER(LEN=255) :: MSG, LOC

    !====================================================================
    ! ConvertSpc_KgKgDry_to_Kg begins here!
    !====================================================================

    ! Assume success
    RC = GC_SUCCESS

    ! Verify correct initial units. If current units are unexpected,
    ! write error message and location to log, then pass failed RC
    ! to calling routine. 
    IF ( TRIM( State_Chm%Spc_Units ) /= 'kg/kg dry' ) THEN
       MSG = 'Incorrect initial units: ' // TRIM( State_Chm%Spc_Units )
       LOC = 'Routine ConvertSpc_KgKgDry_to_Kg in unitconv_mod.F90'
       CALL GC_Error( MSG, RC, LOC )
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
! !IROUTINE: ConvertSpc_kg_to_kgkgdry
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
    USE State_Chm_Mod,      ONLY : ChmState
    USE State_Met_Mod,      ONLY : MetState
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

    !====================================================================
    ! ConvertSpc_Kg_to_KgKgDry begins here!
    !====================================================================

    ! Assume success
    RC = GC_SUCCESS

    ! Verify correct initial units. If current units are unexpected,
    ! write error message and location to log, then pass failed RC
    ! to calling routine. 
    IF ( TRIM( State_Chm%Spc_Units ) /= 'kg' ) THEN
       MSG = 'Incorrect initial units: ' // TRIM( State_Chm%Spc_Units )
       LOC = 'Routine ConvertSpc_Kg_to_KgKgDry in unitconv_mod.F90'
       CALL GC_Error( MSG, RC, LOC )
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
! !IROUTINE: ConvertSpc_mnd_to_kg
!
! !DESCRIPTION: Subroutine ConvertSpc\_MND\_to\_Kg converts the units of 
!  species concentrations from molecular number density (MND) 
!  [molecules/cm3] to mass per grid box [kg].
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE ConvertSpc_MND_to_Kg( am_I_Root, State_Met, State_Chm, RC )
!
! !USES:
!
    USE State_Chm_Mod,      ONLY : ChmState
    USE State_Met_Mod,      ONLY : MetState
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
    ! Scalars
    INTEGER                :: I, J, L, N
    REAL(fp)               :: MolecRatio, MW_kg

    ! Strings
    CHARACTER(LEN=255)     :: MSG, LOC

    !====================================================================
    ! ConvertSpc_MND_to_Kg begins here!
    !====================================================================

    ! Assume success
    RC = GC_SUCCESS

    ! Verify correct initial units. If current units are unexpected,
    ! write error message and location to log, then pass failed RC
    ! to calling routine. 
    IF ( TRIM( State_Chm%Spc_Units ) /= 'molec/cm3' ) THEN
       MSG = 'Incorrect initial units: ' // TRIM( State_Chm%Spc_Units )
       LOC = 'Routine ConvertSpc_MND_to_KgKgDry in unitconv_mod.F90'
       CALL GC_Error( MSG, RC, LOC )
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

    ! Loop over all species
    !$OMP PARALLEL DO                              &
    !$OMP DEFAULT( SHARED                        ) &
    !$OMP PRIVATE( I, J, L, N, MolecRatio, MW_kg )  
    DO N = 1, State_Chm%nSpecies 

       ! Moles C / moles species
       MolecRatio = State_Chm%SpcData(N)%Info%MolecRatio

       ! (Emitted) molecular weight for the species [g]
       ! NOTE: Non-advected species will have a MW of -1, which will
       ! make the species concentration negative.  This can be used to
       ! flag that the species should not be used.  The inverse unit
       ! conversion will flip the sign back to positive (ewl, bmy, 8/4/16)
       MW_kg      = State_Chm%SpcData(N)%Info%emMW_g * 1.e-3_fp

       ! Loop over grid boxes and do the unit conversion
       DO L = 1, LLPAR
       DO J = 1, JJPAR
       DO I = 1, IIPAR
          State_Chm%Species(I,J,L,N) = State_Chm%Species(I,J,L,N)           &
                                     / ( AVO / MW_kg )                      &
                                     * (  State_Met%AIRVOL(I,J,L)           &
                                          * 1e+6_fp * MolecRatio )          
       ENDDO
       ENDDO
       ENDDO

    ENDDO
    !$OMP END PARALLEL DO

    ! Update species units
    State_Chm%Spc_Units = 'kg'

  END SUBROUTINE ConvertSpc_MND_to_Kg
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ConvertSpc_kg_to_mnd
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
    USE State_Chm_Mod,      ONLY : ChmState
    USE State_Met_Mod,      ONLY : MetState
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
    ! Scalars
    INTEGER            :: I, J, L, N
    REAL(fp)           :: MolecRatio, MW_kg

    ! Strings
    CHARACTER(LEN=255) :: MSG, LOC

    !====================================================================
    ! ConvertSpc_Kg_to_MND begins here!
    !====================================================================

    ! Assume success
    RC = GC_SUCCESS

    ! Verify correct initial units. If current units are unexpected,
    ! write error message and location to log, then pass failed RC
    ! to calling routine. 
    IF ( TRIM( State_Chm%Spc_Units ) /= 'kg' ) THEN
       MSG = 'Incorrect initial units: ' // TRIM( State_Chm%Spc_Units )
       LOC = 'Routine ConvertSpc_KgKgDry_to_MND in unitconv_mod.F90'
       CALL GC_Error( MSG, RC, LOC )
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

    ! Loop over all species
    !$OMP PARALLEL DO                              &
    !$OMP DEFAULT( SHARED                        ) &
    !$OMP PRIVATE( I, J, L, N, MolecRatio, MW_kg )  
    DO N = 1, State_Chm%nSpecies 

       ! Moles C / moles species
       MolecRatio = State_Chm%SpcData(N)%Info%MolecRatio

       ! (Emitted) molecular weight for the species [kg]
       ! NOTE: Non-advected species will have a MW of -1, which will
       ! make the species concentration negative.  This can be used to
       ! flag that the species should not be used.  The inverse unit
       ! conversion will flip the sign back to positive (ewl, bmy, 8/4/16)
       MW_kg      = State_Chm%SpcData(N)%Info%emMW_g * 1.e-3_fp

       ! Loop over grid boxes and do the unit conversion
       DO L = 1, LLPAR
       DO J = 1, JJPAR
       DO I = 1, IIPAR
          State_Chm%Species(I,J,L,N) = State_Chm%Species(I,J,L,N)           &
                                     * ( AVO / MW_kg )                      &
                                     / ( State_Met%AIRVOL(I,J,L)            &
                                         * 1e+6_fp * MolecRatio )   
       ENDDO
       ENDDO
       ENDDO

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
! !IROUTINE: ConvertBox_kgkgdry_to_kg
!
! !DESCRIPTION: Subroutine ConvertBox\_KgKgDry\_to\_Kg converts the units of 
!  species concentrations from dry mass mixing ratio [kg tracer/kg dry air]
!  to tracer mass per grid box [kg] for a single grid box. This routine is 
!  temporary during the unit transition of TOMAS to area-independence.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE ConvertBox_KgKgDry_to_Kg( am_I_Root, I, J, L,       &
                                       State_Met, State_Chm, RC ) 
!
! !USES:
!
    USE State_Chm_Mod,      ONLY : ChmState
    USE State_Met_Mod,      ONLY : MetState
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
!  16 Sep 2016 - E. Lundgren - Initial version, an adaptation of 
!                              convertspc_kgkgdry_to_kg
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER            :: N
    CHARACTER(LEN=255) :: MSG, LOC

    !====================================================================
    ! ConvertBox_KgKgDry_to_Kg begins here!
    !====================================================================

    ! Assume success
    RC = GC_SUCCESS

    !$OMP PARALLEL DO           &
    !$OMP DEFAULT( SHARED     ) &
    !$OMP PRIVATE( N ) 
    DO N = 1, State_Chm%nSpecies
       State_Chm%Species(I,J,L,N) = State_Chm%Species(I,J,L,N) &
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
! !IROUTINE: ConvertBox_kg_to_kgkgdry
!
! !DESCRIPTION: Subroutine ConvertBox\_Kg\_to\_KgKgDry converts the units of 
!  species concentrations from species mass per grid box [kg] to mass 
!  mixing ratio [kg tracer/kg dry air] for a single grid box.
!  This routine is temporary during the unit transition of TOMAS to 
!  area-independence.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE ConvertBox_Kg_to_KgKgDry( am_I_Root, I, J, L,         &
                                       State_Met, State_Chm, RC   ) 
!
! !USES:
!
    USE State_Chm_Mod,      ONLY : ChmState
    USE State_Met_Mod,      ONLY : MetState
!
! !INPUT PARAMETERS: 
!
    LOGICAL,        INTENT(IN)    :: am_I_Root   ! Are we on the root CPU?
    INTEGER,        INTENT(IN)    :: I, J, L     ! Grid box indexes
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
!  This routine is temporary and is only used for local conversion of species
!  concentrations for use in TOMAS within wetscav_mod routine WASHOUT. 
!  That routine is called within a parallel do loop and therefore units can
!  only be converted per grid box to avoid excessive computation time. Also,
!  State_Chm%Spc_Units cannot be changed within the parallel do loop without
!  causing problems. It is therefore left out of this routine.
!
! !REVISION HISTORY: 
!  16 Sep 2016 - E. Lundgren - Initial version, an adaptation of 
!                              convertspc_kg_to_kgkgdry
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER :: N
    CHARACTER(LEN=255) :: MSG, LOC

    !====================================================================
    ! ConvertSpc_Kg_to_KgKgDry begins here!
    !====================================================================

    ! Assume success
    RC = GC_SUCCESS

    !$OMP PARALLEL DO           &
    !$OMP DEFAULT( SHARED     ) &
    !$OMP PRIVATE( N ) 
    DO N = 1, State_Chm%nSpecies
       State_Chm%Species(I,J,L,N) = State_Chm%Species(I,J,L,N) &
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
! !IROUTINE: ConvertBox_kgm2_to_kg
!
! !DESCRIPTION: Subroutine ConvertBox\_Kgm2\_to\_Kg converts the units of area 
!  density [kg/m2] to mass [kg] for a single grid box. This routine is
!  temporary during the unit transition of TOMAS to area-independence.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE ConvertBox_Kgm2_to_Kg( am_I_Root, I, J, L,       &
                                    State_Met, State_Chm, RC )
!
! !USES:
!
    USE State_Chm_Mod,      ONLY : ChmState
    USE State_Met_Mod,      ONLY : MetState
!
! !INPUT PARAMETERS: 
!
    LOGICAL,        INTENT(IN)    :: am_I_Root     ! Are we on the root CPU?
    INTEGER,        INTENT(IN)    :: I, J, L       ! Grid box indexes
    TYPE(MetState), INTENT(IN)    :: State_Met     ! Meteorology state object
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ChmState), INTENT(INOUT) :: State_Chm     ! Chemistry state object 
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
!  16 Sep 2016 - E. Lundgren - Rename from ConvertSpc_Kgm2_to_Kg
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
    RC = GC_SUCCESS

    !$OMP PARALLEL DO       &
    !$OMP DEFAULT( SHARED ) &
    !$OMP PRIVATE( N      ) 
    DO N = 1, State_Chm%nSpecies
       State_Chm%Species(I,J,L,N) = State_Chm%Species(I,J,L,N)    &
                                  * State_Met%AREA_M2(I,J,1)
    ENDDO
    !$OMP END PARALLEL DO

  END SUBROUTINE ConvertBox_Kgm2_to_Kg
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ConvertBox_kg_to_kgm2
!
! !DESCRIPTION: Subroutine ConvertBox\_Kg\_to\_kgm2 converts the units of 
! mass [kg] to area density [kg/m2] for a single grid box.  This routine is
!  temporary during the unit transition of TOMAS to area-independence. 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE ConvertBox_Kg_to_Kgm2( am_I_Root, I, J, L,       &
                                    State_Met, State_Chm, RC )
!
! !USES:
!
    USE State_Chm_Mod,      ONLY : ChmState
    USE State_Met_Mod,      ONLY : MetState
!
! !INPUT PARAMETERS: 
!
    LOGICAL,        INTENT(IN)    :: am_I_Root     ! Are we on the root CPU?
    INTEGER,        INTENT(IN)    :: I, J, L       ! Grid box indexes
    TYPE(MetState), INTENT(IN)    :: State_Met     ! Meteorology state object
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ChmState), INTENT(INOUT) :: State_Chm     ! Chemistry state object 
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
!  16 Sep 2016 - E. Lundgren - Rename from ConvertSpc_Kg_to_Kgm2
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER            :: N
    CHARACTER(LEN=255) :: MSG, LOC

    !====================================================================
    ! ConvertBox_Kg_to_Kgm2 begins here!
    !====================================================================

    ! Assume success
    RC = GC_SUCCESS

    !$OMP PARALLEL DO       &
    !$OMP DEFAULT( SHARED ) &
    !$OMP PRIVATE( N      ) 
    DO N = 1, State_Chm%nSpecies
       State_Chm%Species(I,J,L,N) = State_Chm%Species(I,J,L,N)    &
                                  / State_Met%AREA_M2(I,J,1)
    ENDDO
    !$OMP END PARALLEL DO

  END SUBROUTINE ConvertBox_Kg_to_Kgm2
!EOC
END MODULE UnitConv_Mod
