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
  USE ErrCode_Mod
  USE Error_Mod
  USE PhysConstants
  USE Precision_Mod
  USE Input_Opt_Mod,  ONLY : OptInput
  USE State_Chm_Mod,  ONLY : ChmState
  USE State_Chm_Mod,  ONLY : Ind_
  USE State_Grid_Mod, ONLY : GrdState
  USE State_Met_Mod,  ONLY : MetState

  IMPLICIT NONE
  PRIVATE
!
! !DEFINED PARAMETERS:
!
  ! Unit flags (use integers; faster comparisons than strings)
  ! NOTE: mol/mol are equivalent units to v/v
  INTEGER, PARAMETER, PUBLIC :: KG_SPECIES                        = 1
  INTEGER, PARAMETER, PUBLIC :: KG_SPECIES_PER_KG_DRY_AIR         = 2
  INTEGER, PARAMETER, PUBLIC :: KG_SPECIES_PER_KG_TOTAL_AIR       = 3
  INTEGER, PARAMETER, PUBLIC :: KG_SPECIES_PER_M2                 = 4
  INTEGER, PARAMETER, PUBLIC :: MOLECULES_SPECIES_PER_CM3         = 5
  INTEGER, PARAMETER, PUBLIC :: MOLES_SPECIES_PER_MOLES_DRY_AIR   = 6
  INTEGER, PARAMETER, PUBLIC :: MOLES_SPECIES_PER_MOLES_TOTAL_AIR = 7

  ! Labels corresponding to each integer unit flag defined above.
  ! This array is private to this module, whereas the flags are public.
  CHARACTER(LEN=13), PARAMETER, PUBLIC :: UNIT_STR(7) =                   (/ &
     'kg           ', 'kg/kg dry    ', 'kg/kg total  ', 'kg/m2        ',     &
     'molec/cm3    ', 'mol/mol dry  ', 'mol/mol total'                       &
                                                                         /)
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC :: Check_Units
  PUBLIC :: Convert_Spc_Units
  PUBLIC :: Print_Global_Species_Kg
  PUBLIC :: Print_Species_Units

  ! kg/kg dry air <-> kg/grid box (single box only)
  ! Used for TOMAS compatibility in WASHOUT
  PUBLIC :: ConvertBox_KgKgDry_to_Kg
  PUBLIC :: ConvertBox_Kg_to_KgKgDry

  ! kg <-> kg/m2 (single box only)
  ! Used for TOMAS compatibility in WASHOUT within wetscav_mod
  PUBLIC :: ConvertBox_Kgm2_to_Kg
  PUBLIC :: ConvertBox_Kg_to_Kgm2
!
! !PRIVATE MEMBER FUNCTIONS:
!
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! KG/KG DRY <-> V/V DRY
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! kg/kg dry air <-> v/v dry air
  ! Used in DO_TEND in mixing
  PRIVATE :: ConvertSpc_KgKgDry_to_VVDry
  PRIVATE :: ConvertSpc_VVDry_to_KgKgDry

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! KG/KG DRY <-> KG/KG TOTAL
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! kg/kg dry air <-> kg/kg total air
  ! Used in DO_TEND in mixing
  PRIVATE :: ConvertSpc_KgKgDry_to_KgKgTotal
  PRIVATE :: ConvertSpc_KgKgTotal_to_KgKgDry

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! KG/KG DRY <-> KG/M2
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! kg/kg dry air <-> kg/m2
  ! Used for wet deposition, DO_TEND in mixing,
  ! and around AIRQNT and SET_H2O_TRAC in main
  PRIVATE :: ConvertSpc_KgKgDry_to_Kgm2
  PRIVATE :: ConvertSpc_kgm2_to_KgKgDry

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! KG/KG DRY <-> MOLEC/CM3
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  PRIVATE :: ConvertSpc_KgKgDry_to_MND
  PRIVATE :: ConvertSpc_MND_to_KgKgDry

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! AREA-DEPENDENT (temporary routines)
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  ! v/v dry air <-> kg/grid box
  ! Temporarily replaces legacy CONVERT_UNITS
  ! Used in strat_chem_mod and sulfate_mod
  PRIVATE :: ConvertSpc_VVDry_to_Kg
  PRIVATE :: ConvertSpc_Kg_to_VVDry

  ! kg/kg dry air <-> kg/grid box
  ! Used in aerosol_mod, tomas_mod, emissions_mod,
  ! strat_chem_mod, exchange_mod, rrtmg_rad_transfer_mod,
  ! chemistry_mod, sulfate_mod, and carbon_mod
  ! This is since RRTMG, TOMAS, exchange_mod, chemistry,
  ! and EMISSMERCURY are still in [kg]
  PRIVATE :: ConvertSpc_KgKgDry_to_Kg
  PRIVATE :: ConvertSpc_Kg_to_KgKgDry

  ! molec/cm3 dry air <-> kg/gridbox
  PRIVATE :: ConvertSpc_MND_to_Kg
  PRIVATE :: ConvertSpc_Kg_to_MND
!
! !REMARKS:
!  The routines in this module are used to convert the units of
!  species concentrations in various GEOS-Chem routines.
!
! !REVISION HISTORY:
!  23 Jun 2015 - E. Lundgren - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
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
! !IROUTINE: Convert_Spc_Units
!
! !DESCRIPTION: Subroutine Convert\_Spc\_Units is a wrapper function to convert
!  the species input array to a desired unit.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Convert_Spc_Units( Input_Opt, State_Chm,      State_Grid,       &
                                State_Met, new_units,      RC,               &
                                mapping,   previous_units                   )
!
! !USES:
!
    USE TIMERS_MOD
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN)        :: Input_Opt      ! Input Options object
    TYPE(GrdState), INTENT(IN)        :: State_Grid     ! Grid state object
    TYPE(MetState), INTENT(IN)        :: State_Met      ! Met State object
    INTEGER,        INTENT(IN)        :: new_units      ! Units to convert to
    INTEGER,        OPTIONAL, POINTER :: mapping(:)     ! Spc ID -> modelId
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ChmState), INTENT(INOUT)     :: State_Chm      ! Chem State object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)       :: RC             ! Success or failure?
    INTEGER,        OPTIONAL          :: previous_units ! Previous units
!
! !CALLING SEQUENCE:
!
!    ! Convert units from mol/mol dry to kg
!    CALL Convert_Spc_Units(                    &
!         Input_Opt      = Input_Opt,           &
!         State_Chm      = State_Chm,           &
!         State_Grid     = State_Grid,          &
!         State_Met      = State_Met,           &
!         mapping        = State_Chm%Map_XXXXX, & ! Uses Map_All if omitted
!         new_units      = KG_SPECIES,          &
!         previous_units = previous_units,      &
!         RC             = RC                  )
!
!    ...computation...
!
!    ! Convert back to original units
!    CALL Convert_Spc_Units(                    &
!         Input_Opt  = Input_Opt,               &
!         State_Chm  = State_Chm,               &
!         State_Grid = State_Grid,              &
!         State_Met  = State_Met,               &
!         mapping    = State_Chm%Map_XXXXX,     & ! Uses Map_All if omitted
!         new_units  = previous_units,          &
!         RC         = RC                      )
!
! !REVISION HISTORY:
!  14 Apr 2016 - C. Keller    - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    LOGICAL            :: isAdjoint
    INTEGER            :: current_units

    ! Pointers
    INTEGER, POINTER   :: theMapping(:)

    ! Strings
    CHARACTER(LEN=255) :: errNoIn, errNoOut, errMsg, errUnits, thisLoc

    !====================================================================
    ! Convert_Spc_Units begins here!
    !====================================================================
    IF ( Input_Opt%useTimers ) THEN
       CALL Timer_Start( "Unit conversions", RC )
    ENDIF

    ! Convert units for all species unless mapping is passed
    ! NOTE: Avoid an ELSE statement here, which can be a bottleneck.
    theMapping => State_Chm%Map_All
    IF ( PRESENT( mapping ) ) theMapping => mapping

    ! Error check the mapping argument
    IF ( SIZE(theMapping) < 1 ) THEN
       errMsg = 'The "mapping" argument has zero elements!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    ! Initialize
    RC             = GC_SUCCESS
    current_units  = State_Chm%Species(theMapping(1))%Units
    isAdjoint      = .FALSE.
    thisLoc        = ' -> at Convert_Spc_Units (in GeosUtil/unitconv_mod.F90)'

    ! Error messages
    errNoOut  = 'Conversion to ' // TRIM( UNIT_STR( new_units ) )         // &
                ' is not defined!'
    errNoIn   = 'Conversion from ' // TRIM( UNIT_STR( current_units ) )   // &
                ' is not defined!'
    errMsg    = 'Error in conversion from ' // &
                 TRIM( UNIT_STR( current_units ) ) // ' to '              // &
                 TRIM( UNIT_STR( new_units     ) ) //  '!'
    errUnits  = ''

    ! Debugging print
    IF ( Input_Opt%Verbose ) THEN
       WRITE( 6, 100 ) TRIM( UNIT_STR( current_units ) ),                    &
                       TRIM( UNIT_STR( new_units     ) )
 100   FORMAT( '     ### Species Unit Conversion: ', a, ' -> ', a )
    ENDIF

    ! Exit if in and out units are the same
    IF ( new_units == current_units ) THEN
       IF ( PRESENT( previous_units ) ) previous_units = new_units
       IF ( Input_Opt%useTimers ) THEN
          CALL Timer_End( "Unit conversions", RC )
       ENDIF
       RETURN
    ENDIF

    ! Make sure all species have consistent starting units
    IF ( .not. Check_Units( State_Chm, current_units, theMapping ) ) THEN
       errMsg = 'All species do not have consistent starting units!'
       theMapping => NULL()
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

#ifdef ADJOINT
    ! Set a flag if we have compiled for the adjoint
    isAdjoint = Input_Opt%is_Adjoint
#endif

    ! Convert based on input and output units
    SELECT CASE ( current_units )

       !================================================================
       ! Convert from kg/kg dry
       !================================================================
       CASE ( KG_SPECIES_PER_KG_DRY_AIR )

          SELECT CASE ( new_units )

             CASE ( MOLES_SPECIES_PER_MOLES_DRY_AIR )
                CALL ConvertSpc_KgKgDry_to_VVDry(                            &
                     State_Chm,  State_Grid, theMapping,                     &
                     isAdjoint,  RC                                         )

             CASE ( KG_SPECIES_PER_KG_TOTAL_AIR )
                CALL ConvertSpc_KgKgDry_to_KgKgTotal(                        &
                     State_Chm,  State_Grid, State_Met,                      &
                     theMapping, isAdjoint,  RC                              )

             CASE ( KG_SPECIES )
                CALL ConvertSpc_KgKgDry_to_Kg(                               &
                     State_Chm,  State_Grid, State_Met,                      &
                     theMapping, isAdjoint,  RC                             )

             CASE ( KG_SPECIES_PER_M2 )
                CALL ConvertSpc_KgKgDry_to_Kgm2(                             &
                     State_Chm,  State_Grid, State_Met,                      &
                     theMapping, isAdjoint,  RC                             )

             CASE ( MOLECULES_SPECIES_PER_CM3 )
                CALL ConvertSpc_KgKgDry_to_MND(                              &
                     State_Chm,  State_Grid, State_Met,                      &
                     theMapping, isAdjoint,  RC                             )

             CASE DEFAULT
                CALL GC_Error( errNoOut, RC, thisLoc )

          END SELECT

       !================================================================
       ! Convert from kg/kg total
       !================================================================
       CASE ( KG_SPECIES_PER_KG_TOTAL_AIR )

          SELECT CASE ( new_units )

             CASE ( KG_SPECIES_PER_KG_DRY_AIR )
                CALL ConvertSpc_KgKgTotal_to_KgKgDry(                        &
                     State_Chm,  State_Grid, State_Met,                      &
                     theMapping, isAdjoint,  RC                             )

             CASE ( KG_SPECIES )
                CALL ConvertSpc_KgKgTotal_to_KgKgDry(                        &
                     State_Chm,  State_Grid, State_Met,                      &
                     theMapping, isAdjoint,  RC                             )
                CALL ConvertSpc_KgKgDry_to_Kg(                               &
                     State_Chm,  State_Grid, State_Met,                      &
                     theMapping, isAdjoint,  RC                             )

             CASE ( MOLECULES_SPECIES_PER_CM3 )
                CALL ConvertSpc_KgKgTotal_to_KgKgDry(                        &
                     State_Chm,  State_Grid, State_Met,                      &
                     theMapping, isAdjoint,  RC                             )
                CALL ConvertSpc_KgKgDry_to_MND(                              &
                     State_Chm,  State_Grid, State_Met,                      &
                     theMapping, isAdjoint,  RC                             )

             CASE DEFAULT
                CALL GC_Error( errNoOut, RC, thisLoc )
          END SELECT

       !====================================================================
       ! Convert from v/v dry
       !====================================================================
       CASE ( MOLES_SPECIES_PER_MOLES_DRY_AIR )

          SELECT CASE ( new_units )

             CASE ( KG_SPECIES_PER_KG_DRY_AIR )
                CALL ConvertSpc_VVDry_to_KgKgDry(                            &
                     State_Chm,  State_Grid,                                 &
                     theMapping, isAdjoint,  RC                             )

             CASE ( KG_SPECIES )
                CALL ConvertSpc_VVDry_to_Kg(                                 &
                     State_Chm,  State_Grid, State_Met,                      &
                     theMapping, isAdjoint,  RC                             )

             CASE ( KG_SPECIES_PER_M2 )
                CALL ConvertSpc_VVDry_to_KgKgDry(                            &
                     State_Chm,  State_Grid,                                 &
                     theMapping, isAdjoint,  RC                             )
                CALL ConvertSpc_KgKgDry_to_Kgm2(                             &
                     State_Chm,  State_Grid, State_Met,                      &
                     theMapping, isAdjoint,  RC                             )

             CASE DEFAULT
                CALL GC_Error( errNoOut, RC, thisLoc )

          END SELECT

       !====================================================================
       ! Convert from kg
       !====================================================================
       CASE ( KG_SPECIES )

          SELECT CASE ( new_units )

            CASE ( KG_SPECIES_PER_KG_DRY_AIR )
                CALL ConvertSpc_Kg_to_KgKgDry(                               &
                     State_Chm,  State_Grid, State_Met,                      &
                     theMapping, isAdjoint,  RC                             )

             CASE ( KG_SPECIES_PER_KG_TOTAL_AIR )
                CALL ConvertSpc_Kg_to_KgKgDry(                               &
                     State_Chm,  State_Grid, State_Met,                      &
                     theMapping, isAdjoint,  RC                             )
                CALL ConvertSpc_KgKgDry_to_KgKgTotal(                        &
                     State_Chm,  State_Grid, State_Met,                      &
                     theMapping, isAdjoint,  RC                             )

             CASE ( MOLES_SPECIES_PER_MOLES_DRY_AIR )
                CALL ConvertSpc_Kg_to_VVDry(                                 &
                     State_Chm,  State_Grid, State_Met,                      &
                     theMapping, isAdjoint,  RC                             )

             CASE ( MOLECULES_SPECIES_PER_CM3 )
                CALL ConvertSpc_Kg_to_MND(                                   &
                     State_Chm,  State_Grid, State_Met,                      &
                     theMapping, isAdjoint,  RC                             )

             CASE DEFAULT
                CALL GC_Error( errNoOut, RC, thisLoc )

          END SELECT

       !====================================================================
       ! Convert from kg/m2
       !====================================================================
       CASE ( KG_SPECIES_PER_M2 )

          SELECT CASE ( new_units )

             CASE( KG_SPECIES_PER_KG_DRY_AIR )
                CALL ConvertSpc_Kgm2_to_KgKgDry(                             &
                     State_Chm,  State_Grid, State_Met,                      &
                     theMapping, isAdjoint,  RC                             )

             CASE ( MOLES_SPECIES_PER_MOLES_DRY_AIR )
                CALL ConvertSpc_Kgm2_to_KgKgDry(                             &
                     State_Chm,  State_Grid, State_Met,                      &
                     theMapping, isAdjoint,  RC                             )
                CALL ConvertSpc_KgKgDry_to_VVDry(                            &
                     State_Chm,  State_Grid,                                 &
                     theMapping, isAdjoint,  RC                             )

             CASE DEFAULT
                CALL GC_Error( errNoOut, RC, thisLoc )

          END SELECT

       !====================================================================
       ! Convert from molecular number density (MND)
       !====================================================================
       CASE ( MOLECULES_SPECIES_PER_CM3 )

          SELECT CASE ( new_units )

             CASE ( KG_SPECIES )
                CALL ConvertSpc_MND_to_Kg(                                   &
                     State_Chm,  State_Grid, State_Met,                      &
                     theMapping, isAdjoint,  RC                             )

             CASE ( KG_SPECIES_PER_KG_DRY_AIR )
                CALL ConvertSpc_MND_to_KgKgDry(                              &
                     State_Chm,  State_Grid, State_Met,                      &
                     theMapping, isAdjoint,  RC                             )

             CASE ( KG_SPECIES_PER_KG_TOTAL_AIR )
                CALL ConvertSpc_MND_to_KgKgDry(                              &
                     State_Chm,  State_Grid, State_Met,                      &
                     theMapping, isAdjoint,  RC                             )
                CALL ConvertSpc_KgKgDry_to_KgKgTotal(                        &
                     State_Chm,  State_Grid, State_Met,                      &
                     theMapping, isAdjoint,  RC                             )

             CASE DEFAULT
                CALL GC_Error( errNoOut, RC, thisLoc )

          END SELECT

       ! Error if input units not found
       CASE DEFAULT
          CALL GC_Error( errNoIn, RC, thisLoc )

    END SELECT

    !========================================================================
    ! Cleanup and quit
    !========================================================================

    ! Return the previous units (if necessary)
    IF ( PRESENT( previous_units ) ) previous_units = current_units

    ! Free pointer
    theMapping => NULL()

    ! Error if problem within called conversion routine
    IF ( RC /= GC_SUCCESS ) THEN
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    IF ( Input_Opt%useTimers ) THEN
       CALL Timer_End( "Unit conversions", RC )
    ENDIF

  END SUBROUTINE Convert_Spc_Units
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Print_Global_Species_Kg
!
! !DESCRIPTION: Subroutine Print\_Global\_Species\_Kg prints the
!   global and grid box (I,J,L) mass for species N to log. Species
!   units can be any unit for which conversion to kg is defined in
!   the unit conversion module unitconv\_mod.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Print_Global_Species_Kg( I,          J,         L,              &
                                      Spc,        Input_Opt, State_Chm,      &
                                      State_Grid, State_Met, thisLoc,        &
                                      RC                                    )
!
! !INPUT PARAMETERS:
!
    INTEGER,          INTENT(IN)    :: I           ! Grid cell lat index
    INTEGER,          INTENT(IN)    :: J           ! Grid cell lon index
    INTEGER,          INTENT(IN)    :: L           ! Grid cell lev index
    CHARACTER(LEN=*), INTENT(IN)    :: Spc         ! Species abbrev string
    CHARACTER(LEN=*), INTENT(IN)    :: thisLoc     ! Call location string
    TYPE(OptInput),   INTENT(IN)    :: Input_Opt   ! Input Options object
    TYPE(GrdState),   INTENT(IN)    :: State_Grid  ! Grid State object
    TYPE(MetState),   INTENT(IN)    :: State_Met   ! Meteorology State object
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ChmState),   INTENT(INOUT) :: State_Chm   ! Chemistry State object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(OUT)   :: RC          ! Success or failure?
!
! !REMARKS:
!  This routine is for debugging purposes to helptrace where species
!  mass is not conserved.
!
! !REVISION HISTORY:
!  22 Jun 2016 - E. Lundgren - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER            :: N
    INTEGER            :: previous_units
    REAL(fp)           :: SpcTotal

    ! Arrays
    INTEGER, TARGET    :: mapping(1)

    ! Pointers
    INTEGER, POINTER   :: theMapping(:)

    ! Strings
    CHARACTER(LEN=12)  :: SpcName
    CHARACTER(LEN=255) :: errMsg, errLoc

    !========================================================================
    ! Print_Global_Species_Kg begins here!
    !========================================================================

    RC     = GC_SUCCESS
    errMsg = ''
    errLoc = &
     ' -> at Print_Global_Species_Kg (in GeosUtil/unitconv_mod.F90)'

    ! Get species index
    N = Ind_(spc)

    ! Create a pointer for the mapping indices
    mapping(1) =  N
    theMapping => mapping

    ! Convert species conc units to kg
    CALL Convert_Spc_Units(                                                  &
         Input_Opt      = Input_Opt,                                         &
         State_Chm      = State_Chm,                                         &
         State_Grid     = State_Grid,                                        &
         State_Met      = State_Met,                                         &
         mapping        = theMapping,                                        &
         new_units      = KG_SPECIES,                                        &
         previous_units = previous_units,                                    &
         RC             = RC                                                )

    ! Trap potential errors
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Unit conversion error!'
       CALL GC_Error( errMsg, RC, errLoc )
       RETURN
    ENDIF

    ! Echo info
    IF ( Input_Opt%amIRoot ) THEN
       WRITE( 6, 100 ) TRIM( thisLoc )
    ENDIF
100 FORMAT( /, '%%%%% PRINT_GLOBAL_SPECIES_KG at ', a )

    ! Compute global sum
    SpcTotal = SUM( State_Chm%Species(N)%Conc(:,:,:) )

    ! Get species name from the species database
    SpcName = TRIM( State_Chm%SpcData(N)%Info%Name )

    ! Write formatted output
    IF ( Input_Opt%amIRoot ) THEN
       WRITE( 6, 110 ) SpcName, SpcTotal
       WRITE( 6, 115 ) SpcName, State_Chm%Species(N)%Conc(I,J,L), I, J, L
       WRITE( 6, 115 ) 'AD', State_Met%AD(I,J,L), I, J, L
       WRITE( 6, 115 ) 'PREVSPHU', State_Met%SPHU_PREV(I,J,L), I, J, L
       WRITE( 6, 115 ) 'SPHU', State_Met%SPHU(I,J,L), I, J, L
       WRITE( 6, 120 )
    ENDIF
110 FORMAT( 'Global sum [kg] for ', a8, ' = ', es24.16 )
115 FORMAT( 'Grid cell  [kg] for ', a8, ' = ', es24.16, ', I,J,L= ',3I4 )
120 FORMAT( / )

    ! Convert back to original units
    CALL Convert_Spc_Units(                                                  &
         Input_Opt  = Input_Opt,                                             &
         State_Chm  = State_Chm,                                             &
         State_Grid = State_Grid,                                            &
         State_Met  = State_Met,                                             &
         mapping    = theMapping,                                            &
         new_units  = previous_units,                                        &
         RC         = RC                                                    )

    ! Free pointer
    theMapping => NULL()

    ! Trap potential errors
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Unit conversion error!'
       CALL GC_Error( errMsg, RC, errLoc )
       RETURN
    ENDIF

  END SUBROUTINE Print_Global_Species_Kg
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
  SUBROUTINE ConvertSpc_KgKgDry_to_VVDry( State_Chm, State_Grid,             &
                                          mapping,   isAdjoint,  RC         )
!
! !INPUT PARAMETERS:
!
    TYPE(GrdState), INTENT(IN)    :: State_Grid  ! Grid State object
    INTEGER,        INTENT(IN)    :: mapping(:)  ! Species map to modelId
    LOGICAL,        INTENT(IN)    :: isAdjoint   ! Is this reverse integration?
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
!
! !REVISION HISTORY:
!  21 Jul 2016 - E. Lundgren - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER                :: N,      S
    REAL(fp)               :: const

    ! Strings
    CHARACTER(LEN=255)     :: errMsg, thisLoc

    !========================================================================
    ! ConvertSpc_KgKgDry_to_VVDry begins here!
    !========================================================================

    ! Initialize
    RC      = GC_SUCCESS
    errMsg  = ''
    thisLoc = &
     ' -> at ConvertSpc_KgKgDry_to_VVDry (in GeosUtil/unitconv_mod.F90)'


    !========================================================================
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
    !========================================================================

    ! Loop over species
    !$OMP PARALLEL DO                                                        &
    !$OMP DEFAULT( SHARED                                                   )&
    !$OMP PRIVATE( S, N, const                                              )
    DO S = 1, SIZE( mapping )

       ! Get the modelId from the mapping array
       N = mapping(S)

       ! Compute this constant term only once
       const = AIRMW / State_Chm%SpcData(N)%Info%MW_g

       ! Convert species concentration units
       State_Chm%Species(N)%Conc =                                           &
       State_Chm%Species(N)%Conc * const

#ifdef ADJOINT
          IF ( isAdjoint ) THEN
             State_Chm%SpeciesAdj(:,:,:,N) =                                 &
             State_Chm%SpeciesAdj(:,:,:,N) * const
          ENDIF
#endif

       ! Update units metadata
       State_Chm%Species(N)%Units = MOLES_SPECIES_PER_MOLES_DRY_AIR

    ENDDO
    !$OMP END PARALLEL DO

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
  SUBROUTINE ConvertSpc_VVDry_to_KgKgDry( State_Chm, State_Grid,             &
                                          mapping,   isAdjoint,  RC         )
!
! !INPUT PARAMETERS:
!
    TYPE(GrdState), INTENT(IN)    :: State_Grid  ! Grid State object
    INTEGER,        INTENT(IN)    :: mapping(:)  ! Species map to modelId
    LOGICAL,        INTENT(IN)    :: isAdjoint   ! Is this reverse integration?
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
!
! !REVISION HISTORY:
!  21 Jul 2016 - E. Lundgren - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER                :: N,      S
    REAL(fp)               :: const

    ! Strings
    CHARACTER(LEN=255)     :: errMsg, thisLoc

    !========================================================================
    ! ConvertSpc_VVDry_to_KgKgDry begins here!
    !========================================================================

    ! Initialize
    RC      = GC_SUCCESS
    errMsg  = ''
    thisLoc = &
     ' -> at ConvertSpc_VVDry_to_KgKgDry (in GeosUtil/unitconv_mod.F90)'

    !========================================================================
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
    !========================================================================

    ! Loop over all species
    !$OMP PARALLEL DO                                                        &
    !$OMP DEFAULT( SHARED                                                   )&
    !$OMP PRIVATE( S, N, const                                              )
    DO S = 1, SIZE( mapping )

       ! Get the modelId from the mapping array
       N = mapping(S)

       ! Compute this constant only once
       const = AIRMW / State_Chm%SpcData(N)%Info%MW_g

       ! Convert species concentration units
       State_Chm%Species(N)%Conc =                                           &
       State_Chm%Species(N)%Conc / const

#ifdef ADJOINT
       IF ( isAdjoint ) THEN
          State_Chm%SpeciesAdj(:,:,:,N) =                                    &
          State_Chm%SpeciesAdj(:,:,:,N) / const
       ENDIF
#endif

       ! Update units metadata
       State_Chm%Species(N)%Units = KG_SPECIES_PER_KG_DRY_AIR

    ENDDO
    !$OMP END PARALLEL DO

    END SUBROUTINE ConvertSpc_VVDry_to_KgKgDry
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ConvertSpc_kgkgdry_to_kgkgtotal
!
! !DESCRIPTION: Subroutine ConvertSpc\_KgKgDry\_to\_KgKgTotal converts the
!  units of species concentrations from dry mass mixing ratio (KGKG) [kg/kg] to
!  total mass mixing ratio [kg/kg] including moisture.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE ConvertSpc_KgKgDry_to_KgKgTotal( State_Chm, State_Grid,         &
                                              State_Met, mapping,            &
                                              isAdjoint, RC                 )
!
! !INPUT PARAMETERS:
!
    TYPE(GrdState), INTENT(IN)    :: State_Grid  ! Grid State object
    TYPE(MetState), INTENT(IN)    :: State_Met   ! Meteorology State object
    INTEGER,        INTENT(IN)    :: mapping(:)  ! Species map to modelId
    LOGICAL,        INTENT(IN)    :: isAdjoint   ! Is this reverse integration?
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
!
! !REVISION HISTORY:
!  07 Nov 2018 - E. Lundgren - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER                :: N,      S

    ! Strings
    CHARACTER(LEN=255)     :: errMsg, thisLoc

    !========================================================================
    ! ConvertSpc_KgKgDry_to_KgKgTotal begins here!
    !========================================================================

    ! Initialize
    RC      = GC_SUCCESS
    errMsg  = ''
    thisLoc = &
     ' -> at ConvertSpc_KgKgDry_to_KgKgTotal (in GeosUtil/unitconv_mod.F90)'

    !========================================================================
    !
    !  The conversion is as follows:
    !
    !========================================================================

    ! Loop over all species
    !$OMP PARALLEL DO                                                        &
    !$OMP DEFAULT( SHARED                                                   )&
    !$OMP PRIVATE( S, N                                                     )
    DO S = 1, SIZE( mapping )

       ! Get the modelId from the mapping array
       N = mapping(S)

       ! Convert species concentration units
       State_Chm%Species(N)%Conc =                                           &
       State_Chm%Species(N)%Conc *                                           &
          ( 1.0_fp - ( State_Met%SPHU * 1e-3_fp ) )

#ifdef ADJOINT
       IF ( isAdjoint ) THEN
          State_Chm%SpeciesAdj(:,:,:,N) =                                    &
          State_Chm%SpeciesAdj(:,:,:,N) *                                    &
             ( 1.0_fp - ( State_Met%SPHU * 1e-3_fp ) )
       ENDIF
#endif

       ! Update units metadata
       State_Chm%Species(N)%Units = KG_SPECIES_PER_KG_TOTAL_AIR
    ENDDO
    !$OMP END PARALLEL DO

  END SUBROUTINE ConvertSpc_KgKgDry_to_KgKgTotal
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ConvertSpc_kgkgtotal_to_kgkgdry
!
! !DESCRIPTION: Subroutine ConvertSpc\_KgKgTotal\_to\_KgKgDry converts the
!  units of species concentrations from total mass mixing ratio [kg/kg]
!  (includes moisture) to dry mass mixing ratio [kg/kg].
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE ConvertSpc_KgKgTotal_to_KgKgDry( State_Chm, State_Grid,         &
                                              State_Met, mapping,            &
                                              isAdjoint, RC                 )
!
! !INPUT PARAMETERS:
!
    TYPE(GrdState), INTENT(IN)    :: State_Grid  ! Grid State object
    TYPE(MetState), INTENT(IN)    :: State_Met   ! Meteorology State object
    INTEGER,        INTENT(IN)    :: mapping(:)  ! Species map to modelId
    LOGICAL,        INTENT(IN)    :: isAdjoint   ! Is this reverse integration?
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
!
! !REVISION HISTORY:
!  07 Nov 2018 - E. Lundgren - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER                :: N,      S

    ! Strings
    CHARACTER(LEN=255)     :: errMsg, thisLoc

    !========================================================================
    ! ConvertSpc_KgKgTotal_to_KgKgDry begins here!
    !========================================================================

    ! Initialize
    RC      = GC_SUCCESS
    errMsg  = ''
    thisLoc = &
     ' -> at ConvertSpc_KgKgTotal_to_KgKgDry (in GeosUtil/unitconv_mod.F90)'

    !========================================================================
    !
    !  The conversion is as follows:
    !
    !========================================================================

    ! Loop over all species
    !$OMP PARALLEL DO                                                        &
    !$OMP DEFAULT( SHARED                                                   )&
    !$OMP PRIVATE( S, N                                                     )
    DO S = 1, SIZE( mapping )

       ! Get the modelId from the mapping array
       N = mapping(S)

       ! Convert species concentration units
       State_Chm%Species(N)%Conc =                                           &
       State_Chm%Species(N)%Conc /                                           &
          ( 1e0_fp - ( State_Met%SPHU * 1e-3_fp ) )

#ifdef ADJOINT
       IF ( isAdjoint ) THEN
          State_Chm%SpeciesAdj(:,:,:,N) =                                    &
          State_Chm%SpeciesAdj(:,:,:,N) /                                    &
             ( 1e0_fp - ( State_Met%SPHU * 1e-3_fp ) )
       ENDIF
#endif

       ! Update units metadata
       State_Chm%Species(N)%Units = KG_SPECIES_PER_KG_DRY_AIR

    ENDDO
    !$OMP END PARALLEL DO

    END SUBROUTINE ConvertSpc_KgKgTotal_to_KgKgDry
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
  SUBROUTINE ConvertSpc_KgKgDry_to_Kgm2( State_Chm, State_Grid, State_Met,   &
                                         mapping,   isAdjoint,  RC          )
!
! !INPUT PARAMETERS:
!
    TYPE(GrdState), INTENT(IN)    :: State_Grid  ! Grid State object
    TYPE(MetState), INTENT(IN)    :: State_Met   ! Meteorology state object
    INTEGER,        INTENT(IN)    :: mapping(:)  ! Species map to modelId
    LOGICAL,        INTENT(IN)    :: isAdjoint   ! Is this reverse integration?
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry state object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
!
! !REVISION HISTORY:
!  21 Jul 2016 - E. Lundgren - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER            :: N,      S

    ! Strings
    CHARACTER(LEN=255) :: errMsg, thisLoc

    !========================================================================
    ! ConvertSpc_KgKgDry_to_Kgm2 begins here!
    !========================================================================

    ! Initialize
    RC      = GC_SUCCESS
    errMsg  = ''
    thisLoc = &
     ' -> at ConvertSpc_KgKgDry_to_Kgm2 (in GeosUtil/unitconv_mod.F90)'

    !========================================================================
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
    !========================================================================
    !$OMP PARALLEL DO                                                        &
    !$OMP DEFAULT( SHARED                                                   )&
    !$OMP PRIVATE( S, N                                                     )
    DO S = 1, SIZE( mapping )

       ! Get the modelId from the mapping array
       N = mapping(S)

       ! Convert species concentration units
       State_Chm%Species(N)%Conc =                                           &
       State_Chm%Species(N)%Conc * ( g0_100 * State_Met%DELP_DRY )

#ifdef ADJOINT
       IF ( isAdjoint ) THEN
          State_Chm%SpeciesAdj(:,:,:,N) =                                    &
          State_Chm%SpeciesAdj(:,:,:,N) * ( g0_100 * State_Met%DELP_DRY )
       ENDIF
#endif

       ! Update units metadata
       State_Chm%Species(N)%Units = KG_SPECIES_PER_M2

    ENDDO
    !$OMP END PARALLEL DO

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
  SUBROUTINE ConvertSpc_Kgm2_to_KgKgDry( State_Chm, State_Grid, State_Met,   &
                                         mapping,   isAdjoint,  RC          )
!
! !INPUT PARAMETERS:
!
    TYPE(GrdState), INTENT(IN)    :: State_Grid  ! Grid State object
    TYPE(MetState), INTENT(IN)    :: State_Met   ! Meteorology state object
    INTEGER,        INTENT(IN)    :: mapping(:)  ! Species map to modelId
    LOGICAL,        INTENT(IN)    :: isAdjoint   ! Is this reverse integration?
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry state object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
!
! !REVISION HISTORY:
!  21 Jul 2016 - E. Lundgren - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER            :: N,      S

    ! Strings
    CHARACTER(LEN=255) :: errMsg, thisLoc

    !========================================================================
    ! ConvertSpc_Kgm2_to_KgKgDry begins here!
    !========================================================================

    ! Initialize
    RC      = GC_SUCCESS
    errMsg  = ''
    thisLoc = &
     ' -> at ConvertSpc_Kgm2_to_KgKgDry (in GeosUtil/unitconv_mod.F90)'

    !========================================================================
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
    !========================================================================

    !$OMP PARALLEL DO                                                        &
    !$OMP DEFAULT( SHARED                                                   )&
    !$OMP PRIVATE( S, N                                                     )
    DO S = 1, SIZE( mapping )

       ! Get the modelId from the mapping array
       N = mapping(S)

       ! Convert species concentration units
       State_Chm%Species(N)%Conc =                                           &
       State_Chm%Species(N)%Conc / ( g0_100  * State_Met%DELP_DRY )

#ifdef ADJOINT
       IF ( isAdjoint ) THEN
          State_Chm%SpeciesAdj(:,:,:,N) =                                    &
          State_Chm%SpeciesAdj(:,:,:,N) / ( g0_100  * State_Met%DELP_DRY )
       ENDIF
#endif

       ! Update units metadata
       State_Chm%Species(N)%Units = KG_SPECIES_PER_KG_DRY_AIR

    ENDDO
    !$OMP END PARALLEL DO

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
  SUBROUTINE ConvertSpc_KgKgDry_to_MND( State_Chm, State_Grid, State_Met,    &
                                        mapping,   isAdjoint,  RC           )
!
! !INPUT PARAMETERS:
!
    TYPE(GrdState), INTENT(IN)    :: State_Grid  ! Grid State object
    TYPE(MetState), INTENT(IN)    :: State_Met   ! Meteorology state object
    INTEGER,        INTENT(IN)    :: mapping(:)  ! Species map to modelId
    LOGICAL,        INTENT(IN)    :: isAdjoint   ! Is this reverse integration?
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry state object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
!
! !REVISION HISTORY:
!  21 Jul 2016 - E. Lundgren - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER            :: N,       S
    REAL(fp)           :: const,   MW_kg

    ! Strings
    CHARACTER(LEN=255) :: errMsg,  thisLoc

    !========================================================================
    ! ConvertSpc_KgKgDry_to_MND begins here!
    !========================================================================

    ! Initialize
    RC      = GC_SUCCESS
    errMsg  = ''
    thisLoc = &
     ' -> at ConvertSpc_KgKgDry_to_MND (in GeosUtil/unitconv_mod.F90)'

    !========================================================================
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
    !   (1) Use AD/AIRVOL instead of AIRDEN to preserve legacy method
    !========================================================================

    ! Loop over all species
    !$OMP PARALLEL DO                                                        &
    !$OMP DEFAULT( SHARED                                                   )&
    !$OMP PRIVATE( S, N, MW_kg, const                                       )
    DO S = 1, SIZE( mapping )

       ! Get the modelId from the mapping array
       N = mapping(S)

       ! Molecular weight for the species [kg]
       MW_kg = State_Chm%SpcData(N)%Info%MW_g * 1.e-3_fp

       ! Compute this constant term once
       const = ( AVO / MW_kg ) / 1.0e+6_fp

       ! Convert species concentration units
       State_Chm%Species(N)%Conc =                                           &
       State_Chm%Species(N)%Conc * State_Met%AIRDEN * const

#ifdef ADJOINT
       IF ( isAdjoint ) THEN
          State_Chm%SpeciesAdj(:,:,:,N) =                                    &
          State_Chm%SpeciesAdj(:,:,:,N) * State_Met%AIRDEN * const
       ENDIF
#endif

       ! Update units metadata
       State_Chm%Species(N)%Units = MOLECULES_SPECIES_PER_CM3

    ENDDO
    !$OMP END PARALLEL DO

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
  SUBROUTINE ConvertSpc_MND_to_KgKgDry( State_Chm, State_Grid, State_Met,    &
                                        mapping,   isAdjoint, RC            )
!
! !INPUT PARAMETERS:
!
    TYPE(GrdState), INTENT(IN)    :: State_Grid  ! Grid State object
    TYPE(MetState), INTENT(IN)    :: State_Met   ! Meteorology state object
    INTEGER,        INTENT(IN)    :: mapping(:)  ! Species map to modelId
    LOGICAL,        INTENT(IN)    :: isAdjoint   ! Is this reverse integration?
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry state object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
!
! !REVISION HISTORY:
!  21 Jul 2016 - E. Lundgren - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER            :: N,      S
    REAL(fp)           :: const,  MW_kg

    ! Strings
    CHARACTER(LEN=255) :: errMsg, thisLoc

    !========================================================================
    ! ConvertSpc_MND_to_KgKgDry begins here!
    !========================================================================

    ! Initialize
    RC      = GC_SUCCESS
    errMsg  = ''
    thisLoc = &
     ' -> at ConvertSpc_MND_to_KgKgDry (in GeosUtil/unitconv_mod.F90)'

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
    !  (1) Use exact reverse of the mixing ratio -> # density conversion to
    !      avoid numerical noise differences
    !  (2) Use AD/AIRVOL instead of AIRDEN to preserve legacy method
    !
    !========================================================================

    ! Loop over species
    !$OMP PARALLEL DO                                                        &
    !$OMP DEFAULT( SHARED                                                   )&
    !$OMP PRIVATE( S, N, MW_kg, const                                       )
    DO S = 1, SIZE( mapping )

       ! Get the modelId from the mapping array
       N = mapping(S)

       ! Molecular weight for the species [kg]
       MW_kg = State_Chm%SpcData(N)%Info%MW_g * 1.e-3_fp

       ! Compute this constant term once
       const = 1.0e+6_fp / ( AVO / MW_kg )

       ! Convert species concentration units
       State_Chm%Species(N)%Conc =                                           &
       State_Chm%Species(N)%Conc * const / State_Met%AIRDEN

#ifdef ADJOINT
       IF ( isAdjoint ) THEN
          State_Chm%SpeciesAdj(:,:,:,N) =                                    &
          State_Chm%SpeciesAdj(:,:,:,N) * const / State_Met%AIRDEN
       ENDIF
#endif

       ! Update units metadata
       State_Chm%Species(N)%Units = KG_SPECIES_PER_KG_DRY_AIR

    ENDDO
    !$OMP END PARALLEL DO

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
  SUBROUTINE ConvertSpc_VVDry_to_Kg( State_Chm, State_Grid, State_Met,       &
                                     mapping,   isAdjoint,  RC              )
!
! !INPUT PARAMETERS:
!
    TYPE(GrdState), INTENT(IN)    :: State_Grid  ! Grid State object
    TYPE(MetState), INTENT(IN)    :: State_Met   ! Meteorology state object
    INTEGER,        INTENT(IN)    :: mapping(:)  ! Species map to modelId
    LOGICAL,        INTENT(IN)    :: isAdjoint   ! Is this reverse integration?
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
!  This routine replaces legacy routine CONVERT_UNITS and will be removed
!  once GEOS-Chem is entirely area independent
!
! !REVISION HISTORY:
!  21 Jul 2016 - E. Lundgren - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER            :: N,      S
    REAL(fp)           :: const

    ! Strings
    CHARACTER(LEN=255) :: errMsg, thisLoc

    !========================================================================
    ! ConvertSpc_VVDry_to_Kg begins here!
    !========================================================================

    ! Initialize
    RC      = GC_SUCCESS
    errMsg  = ''
    thisLoc = &
     ' -> at ConvertSpc_VVDry_to_Kg (in GeosUtil/unitconv_mod.F90)'

    !========================================================================
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
    !========================================================================

    ! Loop over all species
    !$OMP PARALLEL DO                                                        &
    !$OMP DEFAULT( SHARED                                                   )&
    !$OMP PRIVATE( S, N, const                                              )
    DO S = 1, SIZE( mapping )

       ! Get the modelId from the mapping array
       N = mapping(S)

       ! Compute this constant term only once
       const = AIRMW / State_Chm%SpcData(N)%Info%MW_g

       ! Convert species concentration units
       State_Chm%Species(N)%Conc =                                           &
       State_Chm%Species(N)%Conc * State_Met%AD / const

#ifdef ADJOINT
       IF ( isAdjoint ) THEN
          State_Chm%SpeciesAdj(:,:,:,N) =                                    &
          State_Chm%SpeciesAdj(:,:,:,N) * State_Met%AD / const
       ENDIF
#endif

       ! Update units metadata
       State_Chm%Species(N)%Units = KG_SPECIES

    ENDDO
    !$OMP END PARALLEL DO

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
  SUBROUTINE ConvertSpc_Kg_to_VVDry( State_Chm, State_Grid, State_Met,       &
                                     mapping,   isAdjoint,  RC              )
!
! !INPUT PARAMETERS:
!
    TYPE(GrdState), INTENT(IN)    :: State_Grid  ! Grid State object
    TYPE(MetState), INTENT(IN)    :: State_Met   ! Meteorology state object
    INTEGER,        INTENT(IN)    :: mapping(:)  ! Species map to modelId
    LOGICAL,        INTENT(IN)    :: isAdjoint   ! Is this reverse integration?
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
!  This routine replaces legacy routine CONVERT_UNITS and will be removed
!  once GEOS-Chem is entirely area independent
!
! !REVISION HISTORY:
!  21 Jul 2016 - E. Lundgren - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER            :: N,      S
    REAL(fp)           :: const

    ! Strings
    CHARACTER(LEN=255) :: errMsg, thisLoc

    !========================================================================
    ! ConvertSpc_Kg_to_VVDry begins here!
    !========================================================================

    ! Initialize
    RC      = GC_SUCCESS
    errMsg  = ''
    thisLoc = &
     ' -> at ConvertSpc_Kg_to_VVDry (in GeosUtil/unitconv_mod.F90)'

    !========================================================================
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
    !========================================================================

    ! Loop over all species
    !$OMP PARALLEL DO                                                        &
    !$OMP DEFAULT( SHARED                                                   )&
    !$OMP PRIVATE( S, N, const                                              )
    DO S = 1, SIZE( mapping )

       ! Get the modelId from the mapping array
       N = mapping(S)

       ! Compute this constant only once
       const = AIRMW / State_Chm%SpcData(N)%Info%MW_g

       ! Convert species concentration units
       State_Chm%Species(N)%Conc =                                           &
       State_Chm%Species(N)%Conc * const / State_Met%AD

#ifdef ADJOINT
       IF ( isAdjoint ) THEN
          State_Chm%SpeciesAdj(:,:,:,N) =                                    &
          State_Chm%SpeciesAdj(:,:,:,N) * const / State_Met%AD
       ENDIF
#endif

       ! Update units metadata
       State_Chm%Species(N)%Units = MOLES_SPECIES_PER_MOLES_DRY_AIR

    ENDDO
    !$OMP END PARALLEL DO

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
  SUBROUTINE ConvertSpc_KgKgDry_to_Kg( State_Chm, State_Grid, State_Met,     &
                                       mapping,   isAdjoint, RC             )
!
! !INPUT PARAMETERS:
!
    TYPE(GrdState), INTENT(IN)    :: State_Grid  ! Grid State object
    TYPE(MetState), INTENT(IN)    :: State_Met   ! Meteorology state object
    INTEGER,        INTENT(IN)    :: mapping(:)  ! Species map to modelId
    LOGICAL,        INTENT(IN)    :: isAdjoint   ! Is this reverse integration?
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry state object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
!
! !REVISION HISTORY:
!  21 Jul 2016 - E. Lundgren - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER            :: N,      S

    ! Strings
    CHARACTER(LEN=255) :: errMsg, thisLoc

    !========================================================================
    ! ConvertSpc_KgKgDry_to_Kg begins here!
    !========================================================================

    ! Initialize
    RC      = GC_SUCCESS
    errMsg  = ''
    thisLoc = &
     ' -> at ConvertSpc_KgKgDry_to_Kg (in GeosUtil/unitconv_mod.F90)'

    !========================================================================
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
    !========================================================================

    !$OMP PARALLEL DO                                                        &
    !$OMP DEFAULT( SHARED                                                   )&
    !$OMP PRIVATE( S, N                                                     )
    DO S = 1, SIZE( mapping )

       ! Get the modelId from the mapping array
       N = mapping(S)

       ! Convert species concentration units
       State_Chm%Species(N)%Conc =                                           &
       State_Chm%Species(N)%Conc * State_Met%AD

#ifdef ADJOINT
       IF ( isAdjoint ) THEN
          State_Chm%SpeciesAdj(:,:,:,N) =                                 &
          State_Chm%SpeciesAdj(:,:,:,N) * State_Met%AD
       ENDIF
#endif

       ! Update units metadata
       State_Chm%Species(N)%Units = KG_SPECIES

    ENDDO
    !$OMP END PARALLEL DO

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
  SUBROUTINE ConvertSpc_Kg_to_KgKgDry( State_Chm, State_Grid, State_Met,     &
                                       mapping,   isAdjoint,  RC            )
!
! !INPUT PARAMETERS:
!
    TYPE(GrdState), INTENT(IN)    :: State_Grid  ! Grid State object
    TYPE(MetState), INTENT(IN)    :: State_Met   ! Meteorology state object
    INTEGER,        INTENT(IN)    :: mapping(:)  ! Species map to modelId
    LOGICAL,        INTENT(IN)    :: isAdjoint   ! Is this reverse integration?
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry state object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
!
! !REVISION HISTORY:
!  21 Jul 2016 - E. Lundgren - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER            :: N,      S

    ! Strings
    CHARACTER(LEN=255) :: errMsg, thisLoc

    !========================================================================
    ! ConvertSpc_Kg_to_KgKgDry begins here!
    !========================================================================

    ! Initialize
    RC      = GC_SUCCESS
    errMsg  = ''
    thisLoc = &
     ' -> at ConvertSpc_Kg_to_KgKgDry (in GeosUtil/unitconv_mod.F90)'

    !========================================================================
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
    !========================================================================

    !$OMP PARALLEL DO                                                        &
    !$OMP DEFAULT( SHARED                                                   )&
    !$OMP PRIVATE( S, N                                                     )
    DO S = 1, SIZE( mapping )

       ! Get the modelId from the mapping array
       N = mapping(S)

       ! Convert species concentration units
       State_Chm%Species(N)%Conc =                                           &
       State_Chm%Species(N)%Conc / State_Met%AD

#ifdef ADJOINT
       IF ( isAdjoint ) THEN
          State_Chm%SpeciesAdj(:,:,:,N) =                                  &
          State_Chm%SpeciesAdj(:,:,:,N) / State_Met%AD
       ENDIF
#endif

       ! Update units metadata
       State_Chm%Species(N)%Units = KG_SPECIES_PER_KG_DRY_AIR

    ENDDO
    !$OMP END PARALLEL DO

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
  SUBROUTINE ConvertSpc_MND_to_Kg( State_Chm, State_Grid, State_Met,         &
                                   mapping,   isAdjoint,  RC                )
!
! !INPUT PARAMETERS:
!
    TYPE(GrdState), INTENT(IN)    :: State_Grid  ! Grid State object
    TYPE(MetState), INTENT(IN)    :: State_Met   ! Meteorology state object
    INTEGER,        INTENT(IN)    :: mapping(:)  ! Species map to modelId
    LOGICAL,        INTENT(IN)    :: isAdjoint   ! Is this reverse integration?
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry state object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
!
! !REVISION HISTORY:
!  21 Jul 2016 - E. Lundgren - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER                :: N,      S
    REAL(fp)               :: const,  MW_kg

    ! Strings
    CHARACTER(LEN=255)     :: errMsg, thisLoc

    !========================================================================
    ! ConvertSpc_MND_to_Kg begins here!
    !========================================================================

    ! Initialize
    RC      = GC_SUCCESS
    errMsg  = ''
    thisLoc = &
     ' -> at ConvertSpc_MND_to_Kg (in GeosUtil/unitconv_mod.F90)'

    !========================================================================
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
    !    = Species(I,J,L,N) [molec/cm3] * AIRVOL(I,J,L) / AVO * MW_KG * 1e6
    !    = Species(I,J,L,N) [molec/cm3] * [1e6 / (AVO / MW_KG)] * AIRVOL(I,J,L)
    !
    ! NOTES:
    !  (1) Use exact reverse of the species mass -> # density conversion to
    !      avoid numerical noise differences
    !
    !========================================================================

    ! Loop over all species
    !$OMP PARALLEL DO                                                        &
    !$OMP DEFAULT( SHARED                                                   )&
    !$OMP PRIVATE( S, N, MW_kg, const                                       )
    DO S = 1, SIZE( mapping )

       ! Get the modelId from the mapping array
       N = mapping(S)

       ! Molecular weight for the species [g]
       MW_kg = State_Chm%SpcData(N)%Info%MW_g * 1.e-3_fp

       ! Define this constant term only once
       const = 1.0e6_fp / ( AVO / MW_kg )

       ! Convert species concentration units
       State_Chm%Species(N)%Conc =                                           &
       State_Chm%Species(N)%Conc * const * State_Met%AIRVOL

#ifdef ADJOINT
       IF ( isAdjoint ) THEN
          State_Chm%SpeciesAdj(:,:,:,N) =                                    &
          State_Chm%SpeciesAdj(:,:,:,N) * const * State_Met%AIRVOL
       ENDIF
#endif

       ! Update units metadata
       State_Chm%Species(N)%Units = KG_SPECIES

    ENDDO
    !$OMP END PARALLEL DO

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
  SUBROUTINE ConvertSpc_Kg_to_MND( State_Chm, State_Grid, State_Met,         &
                                   mapping,   isAdjoint,  RC                )
!
! !INPUT PARAMETERS:
!
    TYPE(GrdState), INTENT(IN)    :: State_Grid  ! Grid State object
    TYPE(MetState), INTENT(IN)    :: State_Met   ! Meteorology state object
    INTEGER,        INTENT(IN)    :: mapping(:)  ! Species map to modelId
    LOGICAL,        INTENT(IN)    :: isAdjoint   ! Is this reverse integration?
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry state object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
!
! !REVISION HISTORY:
!  21 Jul 2016 - E. Lundgren - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER            :: N,      S
    REAL(fp)           :: const,  MW_kg

    ! Strings
    CHARACTER(LEN=255) :: errMsg, thisLoc

    !========================================================================
    ! ConvertSpc_Kg_to_MND begins here!
    !========================================================================

    ! Initialize
    RC      = GC_SUCCESS
    errMsg  = ''
    thisLoc = &
     ' -> at ConvertSpc_Kg_to_MND (in GeosUtil/unitconv_mod.F90)'

    !========================================================================
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
    !    = Species(I,J,L,N) [kg] * [ AVO / MW_KG / 1e6 ] / AIRVOL(I,J,L)
    !
    !========================================================================

    ! Loop over all species
    !$OMP PARALLEL DO                                                        &
    !$OMP DEFAULT( SHARED                                                   )&
    !$OMP PRIVATE( S, N, MW_kg, const                                       )
    DO N = 1, State_Chm%nSpecies

       ! Molecular weight for the species [kg]
       MW_kg = State_Chm%SpcData(N)%Info%MW_g * 1.e-3_fp

       ! Compute this constant term only once
       const = ( AVO / MW_kg ) / 1.0e6_fp

       ! Convert species concentration units
       State_Chm%Species(N)%Conc =                                           &
       State_Chm%Species(N)%Conc * const / State_Met%AIRVOL

#ifdef ADJOINT
       IF ( isAdjoint ) THEN
          State_Chm%SpeciesAdj(:,:,:,N) =                                    &
          State_Chm%SpeciesAdj(:,:,:,N) * const / State_Met%AIRVOL
       ENDIF
#endif

       ! Update units metadata
       State_Chm%Species(N)%Units = MOLECULES_SPECIES_PER_CM3

    ENDDO
    !$OMP END PARALLEL DO

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
  SUBROUTINE ConvertBox_KgKgDry_to_Kg( I,         J,         L,              &
                                       State_Met, State_Chm, isAdjoint      )
!
! !INPUT PARAMETERS:
!
    INTEGER,           INTENT(IN)    :: I, J, L     ! Grid box indices
    TYPE(MetState),    INTENT(IN)    :: State_Met   ! Meteorology state object
    LOGICAL,           INTENT(IN)    :: isAdjoint   ! Reverse integration?
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ChmState),    INTENT(INOUT) :: State_Chm   ! Chemistry state object
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
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! LOCAL VARIABLES:
!
    INTEGER :: N

    !========================================================================
    ! ConvertBox_KgKgDry_to_Kg begins here!
    !========================================================================
    DO N = 1, State_Chm%nSpecies

       ! Convert species concentration units
       State_Chm%Species(N)%Conc(I,J,L) =                                    &
       State_Chm%Species(N)%Conc(I,J,L) * State_Met%AD(I,J,L)

#ifdef ADJOINT
       IF ( isAdjoint ) THEN
          State_Chm%SpeciesAdj(I,J,L,N) =                                    &
          State_Chm%SpeciesAdj(I,J,L,N) * State_Met%AD(I,J,L)
       ENDIF
#endif

       ! Update units metadata
       State_Chm%Species(N)%Units = KG_SPECIES
    ENDDO

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
  SUBROUTINE ConvertBox_Kg_to_KgKgDry( I,         J,         L,              &
                                       State_Met, State_Chm, isAdjoint      )
!
! !INPUT PARAMETERS:
!
    INTEGER,           INTENT(IN)    :: I, J, L     ! Grid box indices
    TYPE(MetState),    INTENT(IN)    :: State_Met   ! Meteorology state object
    LOGICAL,           INTENT(IN)    :: isAdjoint   ! Reverse integration?
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ChmState),    INTENT(INOUT) :: State_Chm   ! Chemistry state object
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
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! LOCAL VARIABLES:
!
    INTEGER :: N

    !========================================================================
    ! ConvertBox_Kg_to_KgKgDry begins here!
    !========================================================================
    DO N = 1, State_Chm%nSpecies

       ! Convert species concentration units
       State_Chm%Species(N)%Conc(I,J,L) =                                    &
       State_Chm%Species(N)%Conc(I,J,L) / State_Met%AD(I,J,L)

#ifdef ADJOINT
       IF ( isAdjoint ) THEN
          State_Chm%SpeciesAdj(I,J,L,N) =                                    &
          State_Chm%SpeciesAdj(I,J,L,N) / State_Met%AD(I,J,L)
       ENDIF
#endif

       ! Update units metadata
       State_Chm%Species(N)%Units = KG_SPECIES_PER_KG_DRY_AIR
    ENDDO

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
  SUBROUTINE ConvertBox_Kgm2_to_Kg( I,         J,          L,                &
                                    State_Chm, State_Grid, isAdjoint        )
!
! !INPUT PARAMETERS:
!
    INTEGER,           INTENT(IN)    :: I, J, L     ! Grid box indices
    TYPE(GrdState),    INTENT(IN)    :: State_Grid  ! Grid State object
    LOGICAL,           INTENT(IN)    :: isAdjoint   ! Reverse integration?
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ChmState),    INTENT(INOUT) :: State_Chm   ! Chemistry state object
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
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER :: N

    !========================================================================
    ! ConvertBox_Kgm2_to_Kg begins here!
    !========================================================================
    DO N = 1, State_Chm%nSpecies

       ! Convert species concentration units
       State_Chm%Species(N)%Conc(I,J,L) =                                    &
       State_Chm%Species(N)%Conc(I,J,L) * State_Grid%Area_M2(I,J)

#ifdef ADJOINT
       IF ( isAdjoint ) THEN
          State_Chm%SpeciesAdj(I,J,L,N) =                                    &
          State_Chm%SpeciesAdj(I,J,L,N) * State_Grid%Area_M2(I,J)
       ENDIF
#endif

       ! Update units metadata
       State_Chm%Species(N)%Units = KG_SPECIES
    ENDDO

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
  SUBROUTINE ConvertBox_Kg_to_Kgm2( I,         J,          L,                &
                                    State_Chm, State_Grid, isAdjoint        )
!
! !INPUT PARAMETERS:
!
    INTEGER,           INTENT(IN)    :: I, J, L     ! Grid box indexes
    TYPE(GrdState),    INTENT(IN)    :: State_Grid  ! Grid State object
    LOGICAL,           INTENT(IN)    :: isAdjoint   ! Reverse integration?
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ChmState),    INTENT(INOUT) :: State_Chm   ! Chemistry state object
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
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER :: N

    !========================================================================
    ! ConvertBox_Kg_to_Kgm2 begins here!
    !========================================================================
    DO N = 1, State_Chm%nSpecies

       ! Convert species concentration units
       State_Chm%Species(N)%Conc(I,J,L) =                                    &
       State_Chm%Species(N)%Conc(I,J,L) / State_Grid%Area_M2(I,J)

#ifdef ADJOINT
       IF ( isAdjoint ) THEN
          State_Chm%SpeciesAdj(I,J,L,N) =                                    &
          State_Chm%SpeciesAdj(I,J,L,N) / State_Grid%Area_M2(I,J)
       ENDIF
#endif

       ! Update units metadata
       State_Chm%Species(N)%Units = KG_SPECIES_PER_M2
    ENDDO

  END SUBROUTINE ConvertBox_Kg_to_Kgm2
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Check_Units
!
! !DESCRIPTION: Returns .TRUE. if all species have the same units, or
!  .FALSE. if not.
!\\
! !INTERFACE:
!
  FUNCTION Check_Units( State_Chm, units, mapping ) RESULT( same )
!
! !INPUT PARAMETERS:
!
    INTEGER,           INTENT(IN) :: units      ! Input units flag
    INTEGER, OPTIONAL, POINTER    :: mapping(:) ! Species ID -> modelId

!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ChmState),    INTENT(IN) :: State_Chm  ! Chemistry State object
!
! !RETURN VALUE:
!
    LOGICAL                       :: same       ! All species in same units?
!
! !REVISION HISTORY:
!  30 Nov 2023 - R. Yantosca - Initial version
!  See the subsequent Git history with the gitk browser!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER, POINTER :: theMapping(:)

    !========================================================================
    ! Check_Units begins here!
    !========================================================================

    ! Point to the mapping array (or use all species if not passed)
    theMapping => State_Chm%Map_All
    IF ( PRESENT( mapping ) ) theMapping => mapping

    ! Are all species in the same units?
    same = ALL( State_Chm%Species(theMapping)%Units == units )

    ! Free pointer
    theMapping => NULL()

  END FUNCTION Check_Units
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Print_Species_Units
!
! !DESCRIPTION: Prints each species name and its units.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Print_Species_Units( State_Chm, mapping )
!
! !USES:
!
    USE State_Chm_Mod, ONLY : ChmState
!
! !INPUT PARAMETERS:
!
    TYPE(ChmState), INTENT(IN) :: State_Chm    ! Chemistry state object
    INTEGER, OPTIONAL, POINTER :: mapping(:)   ! Species Id -> modelId
!
! !REVISION HISTORY:
!  23 Feb 2024 - R. Yantosca - Initial version
!  See the subsequent Git history with the gitk browser!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER :: N, S

    ! Pointers
    INTEGER, POINTER :: theMapping(:)

    !========================================================================
    ! Print_Species_Units begins here!
    !========================================================================

    ! Assume all species will be printed if mapping is not passed
    theMapping => State_Chm%Map_All
    IF ( PRESENT( mapping ) ) theMapping => mapping

    ! Loop over species
    DO S = 1, SIZE( theMapping )

       ! Get the modelId for each species
       N = theMapping(S)

       WRITE( 6, 100 ) N, ADJUSTL( State_Chm%SpcData(N)%Info%Name ),         &
                          UNIT_STR( State_Chm%Species(N)%Units )
 100   FORMAT( i5, 1x, a20, 1x, a )
    ENDDO

    ! Free pointer
    theMapping => NULL()

  END SUBROUTINE Print_Species_Units
!EOC
END MODULE UnitConv_Mod
