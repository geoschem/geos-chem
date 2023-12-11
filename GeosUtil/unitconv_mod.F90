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
  PUBLIC :: Convert_Spc_Units
  PUBLIC :: Print_Global_Species_Kg

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
  PRIVATE :: Check_Units

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
  SUBROUTINE Convert_Spc_Units( Input_Opt, State_Chm,     State_Grid,        &
                                State_Met, mapping,       new_units,         &
                                RC,        previous_units                   )
!
! !USES:
!
    USE TIMERS_MOD
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput),   INTENT(IN)    :: Input_Opt      ! Input Options object
    TYPE(GrdState),   INTENT(IN)    :: State_Grid     ! Grid state object
    TYPE(MetState),   INTENT(IN)    :: State_Met      ! Meteorology state object
    INTEGER,          INTENT(IN)    :: mapping(:)     ! Species map to modelId
    INTEGER,          INTENT(IN)    :: new_units      ! Units to convert to
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ChmState),   INTENT(INOUT) :: State_Chm      ! Chemistry state object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(OUT)   :: RC             ! Success or failure?
    INTEGER,          OPTIONAL      :: previous_units ! Previous units
!
! !REMARKS:
!  The purpose of optional output argument origUnit is to enable conversion
!  back to the original units in a second call to Convert_Spc_Units.
!  For example:
!
!      CALL Convert_Spc_Units( Input_Opt, State_Chm, State_Grid,       &
!                              State_Met, unit,      mapping,     RC  )
!      ...computation...
!      CALL Convert_Spc_Units( Input_Opt, State_Chm, State_Grid,       &
!                              State_Met, mapping,   RC               )
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
    INTEGER            :: in_units

    ! Strings
    CHARACTER(LEN=255) :: errNoIn, errNoOut, errMsg, errUnits, thisLoc

    !====================================================================
    ! Convert_Spc_Units begins here!
    !====================================================================
    IF ( Input_Opt%useTimers ) THEN
       CALL Timer_Start( "Unit conversions", RC )
    ENDIF

    ! Initialize
    RC        = GC_SUCCESS
    in_units  = State_Chm%Species(mapping(1))%Units
    isAdjoint = .FALSE.
    thisLoc   = ' -> at Convert_Spc_Units (in GeosUtil/unitconv_mod.F90)'
    errNoOut  = 'Conversion to '            // TRIM( UNIT_STR(newUnits) )  // &
                ' not defined!'
    errNoIn   = 'Conversion from '          // TRIM( UNIT_STR(inUnit ) )  // &
                ' not defined!'
    errMsg    = 'Error in conversion from ' // TRIM( UNIT_STR(inUnit ) )  // &
                ' to '                      // TRIM( UNIT_STR(newUnits) )  // &
                '!'
    errUnits  = ''

    ! TODO: Re-enable debug print
    ! Debugging print
    IF ( Input_Opt%Verbose ) THEN
       WRITE(6,'(a)') '     ### Species Unit Conversion: ' //                &
                      TRIM( UNIT_STR(inUnit ) )            // ' -> ' //      &
                      TRIM( UNIT_STR(newUnits) )            // ' ###'
    ENDIF

    ! Exit if in and out units are the same
    IF ( newUnits == inUnit ) THEN
       IF ( Input_Opt%useTimers ) THEN
          CALL Timer_End( "Unit conversions", RC )
       ENDIF
       RETURN
    ENDIF

    ! Make sure all species have consistent starting units
    CALL Check_Units( State_Chm, mapping, in_units, RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error encountered in routine "Check_Units"!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

#ifdef ADJOINT
    ! Set a flag if we have compiled for the adjoint
    isAdjoint = Input_Opt%is_Adjoint
#endif

    ! Convert based on input and output units
    SELECT CASE ( in_units )

       !================================================================
       ! Convert from kg/kg dry
       !================================================================
       CASE ( KG_SPECIES_PER_KG_DRY_AIR )

          SELECT CASE ( new_units )

             CASE ( MOLES_SPECIES_PER_MOLES_DRY_AIR )
                CALL ConvertSpc_KgKgDry_to_VVDry(                            &
                     State_Chm,  State_Grid,  mapping,                       &
                     isAdjoint,  RC                                         )

             CASE ( KG_SPECIES_PER_KG_TOTAL_AIR )
                CALL ConvertSpc_KgKgDry_to_KgKgTotal(                        &
                     State_Chm, State_Grid, State_Met,                       &
                     mapping,   isAdjoint,  RC                              )

             CASE ( KG_SPECIES )
                CALL ConvertSpc_KgKgDry_to_Kg(                               &
                     State_Chm, State_Grid, State_Met,                       &
                     mapping,   isAdjoint,  RC                              )

             CASE ( KG_SPECIES_PER_M2 )
                CALL ConvertSpc_KgKgDry_to_Kgm2(                             &
                     State_Chm, State_Grid, State_Met,                       &
                     mapping,   isAdjoint,  RC                              )

             CASE ( MOLECULES_SPECIES_PER_CM3 )
                CALL ConvertSpc_KgKgDry_to_MND(                              &
                     State_Chm, State_Grid, State_Met,                       &
                     mapping,   isAdjoint,  RC                              )

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
                     State_Chm, State_Grid, State_Met,                       &
                     mapping,   isAdjoint,  RC                              )

             CASE ( KG_SPECIES )
                CALL ConvertSpc_KgKgTotal_to_KgKgDry(                        &
                     State_Chm, State_Grid, State_Met,                       &
                     mapping,   isAdjoint,  RC                              )
                CALL ConvertSpc_KgKgDry_to_Kg(                               &
                     State_Chm, State_Grid, State_Met,                       &
                     mapping,   isAdjoint,  RC                              )

             CASE ( MOLECULES_SPECIES_PER_CM3 )
                CALL ConvertSpc_KgKgTotal_to_KgKgDry(                        &
                     State_Chm, State_Grid, State_Met,                       &
                     mapping,   isAdjoint,  RC                              )
                CALL ConvertSpc_KgKgDry_to_MND(                              &
                     State_Chm, State_Grid, State_Met,                       &
                     mapping,   isAdjoint,  RC                              )

             CASE DEFAULT
                CALL GC_Error( errNoOut, RC, thisLoc )
          END SELECT

       !====================================================================
       ! Convert from v/v dry
       !====================================================================
       CASE ( MOLES_SPECIES_PER_MOLES_DRY_AIR )

          SELECT CASE ( new_units )

             CASE ( KG_SPECIES_PER_KG_DRY_AIR )
                CALL ConvertSpc_VVDry_to_KgKgDry(                           &
                     State_Chm, State_Grid,                                 &
                     mapping,   isAdjoint,  RC                             )

             CASE ( KG_SPECIES )
                CALL ConvertSpc_VVDry_to_Kg(                                &
                     State_Chm, State_Grid, State_Met,                      &
                     mapping,   isAdjoint,  RC                             )

             CASE ( KG_SPECIES_PER_M2 )
                CALL ConvertSpc_VVDry_to_KgKgDry(                           &
                     State_Chm, State_Grid,                                 &
                     mapping,   isAdjoint,  RC                             )
                CALL ConvertSpc_KgKgDry_to_Kgm2(                            &
                     State_Chm, State_Grid, State_Met,                      &
                     mapping,   isAdjoint,  RC                             )

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
                     State_Chm, State_Grid, State_Met,                       &
                     mapping,   isAdjoint,  RC                              )

             CASE ( KG_SPECIES_PER_KG_TOTAL_AIR )
                CALL ConvertSpc_Kg_to_KgKgDry(                               &
                     State_Chm, State_Grid, State_Met,                       &
                     mapping,   isAdjoint,  RC                              )
                CALL ConvertSpc_KgKgDry_to_KgKgTotal(                        &
                     State_Chm, State_Grid, State_Met,                       &
                     mapping,   isAdjoint,  RC                              )

             CASE ( MOLES_SPECIES_PER_MOLES_DRY_AIR )
                CALL ConvertSpc_Kg_to_VVDry(                                 &
                     State_Chm, State_Grid, State_Met,                       &
                     mapping,   isAdjoint,  RC                              )

             CASE ( MOLECULES_SPECIES_PER_CM3 )
                CALL ConvertSpc_Kg_to_MND(                                   &
                     State_Chm, State_Grid, State_Met,                       &
                     mapping,   isAdjoint,  RC                              )

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
                     State_Chm, State_Grid, State_Met,                       &
                     mapping,   isAdjoint,  RC                              )

             CASE ( MOLES_SPECIES_PER_MOLES_DRY_AIR )
                CALL ConvertSpc_Kgm2_to_KgKgDry(                             &
                     State_Chm, State_Grid, State_Met,                       &
                     mapping,   isAdjoint,  RC                              )
                CALL ConvertSpc_KgKgDry_to_VVDry(                            &
                     State_Chm, State_Grid,                                  &
                     mapping,   isAdjoint,  RC )

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
                     State_Chm, State_Grid, State_Met,                       &
                     mapping,   isAdjoint,  RC                              )

             CASE ( KG_SPECIES_PER_KG_DRY_AIR )
                CALL ConvertSpc_MND_to_KgKgDry(                              &
                     State_Chm, State_Grid, State_Met,                       &
                     mapping,   isAdjoint,  RC                              )

             CASE ( KG_SPECIES_PER_KG_TOTAL_AIR )
                CALL ConvertSpc_MND_to_KgKgDry(                              &
                     State_Chm, State_Grid, State_Met,                       &
                     mapping,   isAdjoint,  RC                              )
                CALL ConvertSpc_KgKgDry_to_KgKgTotal(                        &
                     State_Chm, State_Grid, State_Met,                       &
                     mapping,   isAdjoint,  RC                              )

             CASE DEFAULT
                CALL GC_Error( errNoOut, RC, thisLoc )

          END SELECT

       ! Error if input units not found
       CASE DEFAULT
          CALL GC_Error( errNoIn, RC, thisLoc )

    END SELECT

    !========================================================================
    ! Additional error checks
    !========================================================================

    ! Make sure that all species have consistent "previous_units" values
    IF ( PRESENT( previous_units ) ) THEN
       CALL Check_Previous_Units( State_Chm, mapping, in_units, RC )
       IF ( RC /= GC_SUCCESS ) THEN
          errUnits = 'Error encountered in "Check_Previous_Units!"'
          CALL GC_Error( errUnits, RC, thisLoc )
          RETURN
       ENDIF
       previous_units = in_units
    ENDIF

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
    INTEGER            :: N
    INTEGER            :: previous_units
    REAL(fp)           :: SpcTotal
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

    ! Convert species conc units to kg
    CALL Convert_Spc_Units(                                                  &
         Input_Opt      = Input_Opt,                                         &
         State_Chm      = State_Chm,                                         &
         State_Grid     = State_Grid,                                        &
         State_Met      = State_Met,                                         &
         mapping        = (/ N /),                                           &
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
         mapping    = (/ N /),                                               &
         new_units  = previous_units,                                        &
         RC         = RC                                                    )

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
    INTEGER                :: I,      J
    INTEGER                :: L,      N,        S
    REAL(fp)               :: MW_g,   MW_ratio

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
    !$OMP PRIVATE( I, J, L, N, S, MW_g, MW_ratio                            )
    DO S = 1, SIZE( mapping )

       ! Get the modelId from the mapping array
       N = mapping(S)

       ! Molecular weight for the species [g]
       MW_g = State_Chm%SpcData(N)%Info%MW_g

       ! Compute the ratio (MW air / MW species) outside of the IJL loop
       MW_ratio = ( AIRMW / MW_g )

       ! Loop over grid boxes and do unit conversion
       DO L = 1, State_Grid%NZ
       DO J = 1, State_Grid%NY
       DO I = 1, State_Grid%NX

          State_Chm%Species(N)%Conc(I,J,L) =                                 &
          State_Chm%Species(N)%Conc(I,J,L) * MW_ratio

#ifdef ADJOINT
          IF ( isAdjoint ) THEN
             State_Chm%SpeciesAdj(I,J,L,N) =                                 &
             State_Chm%SpeciesAdj(I,J,L,N) * MW_ratio
          ENDIF
#endif
       ENDDO
       ENDDO
       ENDDO

       ! Update units metadata
       State_Chm%Species(N)%Previous_Units = State_Chm%Species(N)%Units
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
    INTEGER                :: I, J, L, N, S
    REAL(fp)               :: MW_g,   MW_ratio

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
    !$OMP PRIVATE( I, J, L, N, S, MW_g, MW_ratio                            )
    DO S = 1, SIZE( mapping )

       ! Get the modelId from the mapping array
       N = mapping(S)

       ! Molecular weight for the species [g]
       MW_g = State_Chm%SpcData(N)%Info%MW_g

       ! Compute the ratio (MW air / MW species) outside of the IJL loop
       MW_ratio = ( AIRMW / MW_g )

       ! Loop over grid boxes and do the unit conversion
       DO L = 1, State_Grid%NZ
       DO J = 1, State_Grid%NY
       DO I = 1, State_Grid%NX

          State_Chm%Species(N)%Conc(I,J,L) =                                 &
          State_Chm%Species(N)%Conc(I,J,L) / MW_ratio

#ifdef ADJOINT
          IF ( isAdjoint ) THEN
             State_Chm%SpeciesAdj(I,J,L,N) =                                 &
             State_Chm%SpeciesAdj(I,J,L,N) / MW_ratio
          ENDIF
#endif
       ENDDO
       ENDDO
       ENDDO

       ! Update units metadata
       State_Chm%Species(N)%Previous_Units = State_Chm%Species(N)%Units
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
    INTEGER                :: I, J, L, N, S
    REAL(fp)               :: convFac

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
    !$OMP PRIVATE( I, J, L, N, S, convFac                                   )
    DO S = 1, SIZE( mapping )

       ! Get the modelId from the mapping array
       N = mapping(S)

       ! Loop over grid boxes and do unit conversion
       DO L = 1, State_Grid%NZ
       DO J = 1, State_Grid%NY
       DO I = 1, State_Grid%NX

          convFac = ( 1.0_fp - ( State_Met%SPHU(I,J,L) * 1.0e-3_fp ) )

          State_Chm%Species(N)%Conc(I,J,L) =                                 &
          State_Chm%Species(N)%Conc(I,J,L) * convFac

#ifdef ADJOINT
          IF ( isAdjoint ) THEN
             State_Chm%SpeciesAdj(I,J,L,N) =                                 &
             State_Chm%SpeciesAdj(I,J,L,N) * convFac
          ENDIF
#endif

       ENDDO
       ENDDO
       ENDDO

       ! Update units metadata
       State_Chm%Species(N)%Previous_Units = State_Chm%Species(N)%Units
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
    INTEGER                :: I, J, L, N, S
    REAL(fp)               :: convFac

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
    !$OMP PRIVATE( I, J, L, N, S, convFac                                   )
    DO S = 1, SIZE( mapping )

       ! Get the modelId from the mapping array
       N = mapping(S)

       ! Loop over grid boxes and convert units
       DO L = 1, State_Grid%NZ
       DO J = 1, State_Grid%NY
       DO I = 1, State_Grid%NX

          convFac = ( 1.0_fp - ( State_Met%SPHU(I,J,L) * 1.0e-3_fp ) )

          State_Chm%Species(N)%Conc(I,J,L) =                                 &
          State_Chm%Species(N)%Conc(I,J,L) / convFac

#ifdef ADJOINT
          IF ( isAdjoint ) THEN
             State_Chm%SpeciesAdj(I,J,L,N) =                                 &
             State_Chm%SpeciesAdj(I,J,L,N) / convFac
          ENDIF
#endif
       ENDDO
       ENDDO
       ENDDO

       ! Update units metadata
       State_Chm%Species(N)%Previous_Units = State_Chm%Species(N)%Units
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
                                         mapping    isAdjoint,  RC          )
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
    INTEGER            :: I, J, L, N, S
    REAL(fp)           :: convFac

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
    !$OMP PRIVATE( I, J, L, N, S, convFac                                   )
    DO S = 1, SIZE( mapping )

       ! Get the modelId from the mapping array
       N = mapping(S)

       ! Loop over grid boxes and convert units
       DO L = 1, State_Grid%NZ
       DO J = 1, State_Grid%NY
       DO I = 1, State_Grid%NX

          convFac = ( g0_100 * State_Met%DELP_DRY(I,J,L) )

          State_Chm%Species(N)%Conc(I,J,L) =                                 &
          State_Chm%Species(N)%Conc(I,J,L) * convFac

#ifdef ADJOINT
          IF ( isAdjoint ) THEN
             State_Chm%SpeciesAdj(I,J,L,N) =                                 &
             State_Chm%SpeciesAdj(I,J,L,N) * convFac
          ENDIF
#endif
       ENDDO
       ENDDO
       ENDDO

       ! Update units metadata
       State_Chm%Species(N)%Previous_Units = State_Chm%Species(N)%Units
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
    INTEGER            :: I, J, L, N, S
    REAL(fp)           :: convFac

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
    !$OMP PRIVATE( I, J, L, N, S, convFac                                   )
    DO S = 1, SIZE( mapping )

       ! Get the modelId from the mapping array
       N = mapping(S)

       DO L = 1, State_Grid%NZ
       DO J = 1, State_Grid%NY
       DO I = 1, State_Grid%NX

          convFac = ( 1.0_fp / ( g0_100 * State_Met%DELP_DRY(I,J,L) )

          State_Chm%Species(N)%Conc(I,J,L) =                                 &
          State_Chm%Species(N)%Conc(I,J,L) * convFac

#ifdef ADJOINT
          IF ( isAdjoint ) THEN
             State_Chm%SpeciesAdj(I,J,L,N) =                                 &
             State_Chm%SpeciesAdj(I,J,L,N) * convFac
          ENDIF
#endif

       ENDDO
       ENDDO
       ENDDO

       ! Update units metadata
       State_Chm%Species(N)%Previous_Units = State_Chm%Species(N)%Units
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
    INTEGER            :: I, J, L, N, S
    REAL(fp)           :: avoTerm, convFac, MW_kg

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
    !$OMP PRIVATE( I, J, L, N, S, MW_kg, convFac                            )
    DO S = 1, SIZE( mapping )

       ! Get the modelId from the mapping array
       N = mapping(S)

       ! Molecular weight for the species [kg]
       MW_kg = State_Chm%SpcData(N)%Info%MW_g * 1.e-3_fp

       ! Compute AVO/ MW_kg / 1e6 term outside the IJL loop
       avoTerm = ( AVO / MW_kg ) / 1.0e+6_fp

       ! Loop over grid boxes and do the unit conversion
       DO L = 1, State_Grid%NZ
       DO J = 1, State_Grid%NY
       DO I = 1, State_Grid%NX

          convFac = ( State_Met%AIRDEN(I,J,L) * avoTerm )

          State_Chm%Species(N)%Conc(I,J,L) =                                 &
          State_Chm%Species(N)%Conc(I,J,L) * convFac

#ifdef ADJOINT
          IF ( isAdjoint ) THEN
             State_Chm%SpeciesAdj(I,J,L,N) =                                 &
             State_Chm%SpeciesAdj(I,J,L,N) * convFac
          ENDIF
#endif
       ENDDO
       ENDDO
       ENDDO

       ! Update units metadata
       State_Chm%Species(N)%Previous_Units = State_Chm%Species(N)%Units
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
                                        mappingm,  isAdjoint, RC            )
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
    INTEGER            :: I, J, L, N, S
    REAL(fp)           :: avoTerm, convFac, MW_kg
    CHARACTER(LEN=255) :: errMsg,  thisLoc

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
    !$OMP PRIVATE( I, J, L, N, S, MW_kg, avoTerm, convFac                   )
    DO S = 1, SIZE( mapping )

       ! Get the modelId from the mapping array
       N = mapping(S)

       ! Molecular weight for the species [kg]
       MW_kg = State_Chm%SpcData(N)%Info%MW_g * 1.e-3_fp

       ! Compute this term outside of the IJL loop
       avoTerm = ( 1.0e+6_fp / ( AVO / MW_kg ) )

       ! Loop over grid boxes and do the unit conversion
       DO L = 1, State_Grid%NZ
       DO J = 1, State_Grid%NY
       DO I = 1, State_Grid%NX

          convFac = ( avoTerm / State_Met%AIRDEN(I,J,L) )

          State_Chm%Species(N)%Conc(I,J,L) =                                 &
          State_Chm%Species(N)%Conc(I,J,L) * convFac
                                    
#ifdef ADJOINT
          IF ( isAdjoint ) THEN
             State_Chm%SpeciesAdj(I,J,L,N) =                                 &
             State_Chm%SpeciesAdj(I,J,L,N) * convFac
          ENDIF
#endif
       ENDDO
       ENDDO
       ENDDO

       ! Update units metadata
       State_Chm%Species(N)%Previous_Units = State_Chm%Species(N)%Units
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
    INTEGER            :: I, J, L, N, S
    REAL(fp)           :: convFac, MW_g,   MW_ratio
    CHARACTER(LEN=255) :: errMsg,  thisLoc

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
    !$OMP PRIVATE( I, J, L, N, S, MW_g, MW_ratio, convFac                   )
    DO S = 1, SIZE( mapping )

       ! Get the modelId from the mapping array
       N = mapping(S)

       ! Molecular weight for the species [g]
       MW_g = State_Chm%SpcData(N)%Info%MW_g

       ! Compute this term outside of the IJL loop
       MW_ratio = ( AIRMW / MW_g )

       ! Loop over grid boxes convert units
       DO L = 1, State_Grid%NZ
       DO J = 1, State_Grid%NY
       DO I = 1, State_Grid%NX

          convFac = ( State_Met%AD(I,J,L) / MW_ratio )

          State_Chm%Species(N)%Conc(I,J,L) =                                 &
          State_Chm%Species(N)%Conc(I,J,L) * convFac
                                      
#ifdef ADJOINT
          IF ( isAdjoint ) THEN
             State_Chm%SpeciesAdj(I,J,L,N) =                                 &
             State_Chm%SpeciesAdj(I,J,L,N) * convFac
          ENDIF
#endif
       ENDDO
       ENDDO
       ENDDO

       ! Update units metadata
       State_Chm%Species(N)%Previous_Units = State_Chm%Species(N)%Units
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
    INTEGER            :: I, J, L, N, S
    REAL(fp)           :: convFac, MW_g, MW_ratio

    ! Strings
    CHARACTER(LEN=255) :: errMsg,  thisLoc

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
    !$OMP PRIVATE( I, J, L, N, S, MW_g, MW_ratio, convFac                   )
    DO S = 1, SIZE( mapping )

       ! Get the modelId from the mapping array
       N = mapping(S)

       ! Molecular weight for the species [g]
       MW_g = State_Chm%SpcData(N)%Info%MW_g

       ! Compute this term outside the IJL loop
       MW_ratio = ( AIRMW / MW_g )

       ! Loop over grid boxes and do the unit conversion
       DO L = 1, State_Grid%NZ
       DO J = 1, State_Grid%NY
       DO I = 1, State_Grid%NX

          convFac = ( MW_ratio / State_Met%AD(I,J,L) )

          State_Chm%Species(N)%Conc(I,J,L) =                                 &
          State_Chm%Species(N)%Conc(I,J,L) * convFac
                                       
#ifdef ADJOINT
          IF ( isAdjoint ) THEN
             State_Chm%SpeciesAdj(I,J,L,N) =                                 &
             State_Chm%SpeciesAdj(I,J,L,N) * convFac
#endif
       ENDDO
       ENDDO
       ENDDO

       ! Update units metadata
       State_Chm%Species(N)%Previous_Units = State_Chm%Species(N)%Units
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
    INTEGER            :: I, J, L, N, S
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
    !$OMP PRIVATE( I, J, L, N, S                                            )
    DO S = 1, SIZE( mapping )

       ! Get the modelId from the mapping array
       N = mapping(S)

       ! Loop over grid boxes and convert units
       DO L = 1, State_Grid%NZ
       DO J = 1, State_Grid%NY
       DO I = 1, State_Grid%NX

          State_Chm%Species(N)%Conc(I,J,L) =                                 &
          State_Chm%Species(N)%Conc(I,J,L) * State_Met%AD(I,J,L)

#ifdef ADJOINT
          IF ( isAdjoint ) THEN
             State_Chm%SpeciesAdj(I,J,L,N) =                                 &
             State_Chm%SpeciesAdj(I,J,L,N) * State_Met%AD(I,J,L)
          ENDIF
#endif
       ENDDO
       ENDDO
       ENDDO

       ! Update units metadata
       State_Chm%Species(N)%Previous_Units = State_Chm%Species(N)%Units
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
    INTEGER            :: I, J, L, N, S

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
    !$OMP PRIVATE( I, J, L, N, S                                            )
    DO S = 1, SIZE( mapping )

       ! Get the modelId from the mapping array
       N = mapping(S)

       DO L = 1, State_Grid%NZ
       DO J = 1, State_Grid%NY
       DO I = 1, State_Grid%NX

          State_Chm%Species(N)%Conc(I,J,L) =                                 &
          State_Chm%Species(N)%Conc(I,J,L) / State_Met%AD(I,J,L)

#ifdef ADJOINT
          IF ( isAdjoint ) THEN
            State_Chm%SpeciesAdj(I,J,L,N) =                                  &
            State_Chm%SpeciesAdj(I,J,L,N) / State_Met%AD(I,J,L)
         ENDIF
#endif
       ENDDO
       ENDDO
       ENDDO

       ! Update units metadata
       State_Chm%Species(N)%Previous_Units = State_Chm%Species(N)%Units
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
    INTEGER                :: I, J, L, N, S
    REAL(fp)               :: avoTerm, MW_kg

    ! Strings
    CHARACTER(LEN=255)     :: errMsg,  thisLoc

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
    !    = Species(I,J,L,N) [molecules/cm3] * AIRVOL(I,J,L) / AVO * MW_KG * 1e6
    !
    ! NOTES:
    !  (1) Use exact reverse of the species mass -> # density conversion to
    !      avoid numerical noise differences
    !
    !========================================================================

    ! Loop over all species
    !$OMP PARALLEL DO                                                        &
    !$OMP DEFAULT( SHARED                                                   )&
    !$OMP PRIVATE( I, J, L, N, S, MW_kg, avoTerm                            )
    DO S = 1, SIZE( mapping )

       ! Get the modelId from the mapping array
       N = mapping(S)

       ! Molecular weight for the species [g]
       MW_kg = State_Chm%SpcData(N)%Info%MW_g * 1.e-3_fp

       ! Compute this term ouside the IJL loop
       avoTerm = ( AVO / MW_kg )

       ! Loop over grid boxes and do the unit conversion
       DO L = 1, State_Grid%NZ
       DO J = 1, State_Grid%NY
       DO I = 1, State_Grid%NX

          convFac = ( State_Met%AIRVOL(I,J,L) * 1.0e+6_fp ) / avoTerm

          State_Chm%Species(N)%Conc(I,J,L) =                                 &
          State_Chm%Species(N)%Conc(I,J,L) / convFac

#ifdef ADJOINT
          IF ( isAdjoint ) THEN
               State_Chm%SpeciesAdj(I,J,L,N) =                               &
               State_Chm%SpeciesAdj(I,J,L,N) * convFac
          ENDIF
#endif
       ENDDO
       ENDDO
       ENDDO

       ! Update units metadata
       State_Chm%Species(N)%Previous_Units = State_Chm%Species(N)%Units
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
    INTEGER            :: I, J, L, N, S
    REAL(fp)           :: avoTerm, MW_kg

    ! Strings
    CHARACTER(LEN=255) :: errMsg,  thisLoc

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
    !
    !========================================================================

    ! Loop over all species
    !$OMP PARALLEL DO                                                        &
    !$OMP DEFAULT( SHARED                                                   )&
    !$OMP PRIVATE( I, J, L, N, MW_kg, avoTerm                               )
    DO N = 1, State_Chm%nSpecies

       ! Molecular weight for the species [kg]
       MW_kg = State_Chm%SpcData(N)%Info%MW_g * 1.e-3_fp

       ! Compute this term outside the IJL loop
       avoTerm = ( AVO / MW_kg )

       ! Loop over grid boxes and do the unit conversion
       DO L = 1, State_Grid%NZ
       DO J = 1, State_Grid%NY
       DO I = 1, State_Grid%NX

          convFac = avoTerm / ( State_Met%AIRVOL(I,J,L) * 1.0e+6_fp )

          State_Chm%Species(N)%Conc(I,J,L) =                                 &
          State_Chm%Species(N)%Conc(I,J,L) * convFac

#ifdef ADJOINT
          IF ( isAdjoint ) THEN
             State_Chm%SpeciesAdj(I,J,L,N) =                                 &
             State_Chm%SpeciesAdj(I,J,L,N) * convFac
         ENDIF
#endif
       ENDDO
       ENDDO
       ENDDO

       ! Update units metadata
       State_Chm%Species(N)%Previous_Units = State_Chm%Species(N)%Units
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
  SUBROUTINE ConvertBox_KgKgDry_to_Kg( I,         J,         L,             &
                                       State_Met, State_Chm, mapping,       &
                                       isAdjoint, RC                       )
!
! !INPUT PARAMETERS:
!
    INTEGER,        INTENT(IN)    :: I, J, L     ! Grid box indexes
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
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER            :: N, S

    ! Strings
    CHARACTER(LEN=255) :: errMsg, thisLoc

    !========================================================================
    ! ConvertBox_KgKgDry_to_Kg begins here!
    !========================================================================

    ! Initialize
    RC      = GC_SUCCESS
    errMsg  = ''
    thisLoc = &
     ' -> at ConvertBox_KgKgDry_to_Kg (in GeosUtil/unitconv_mod.F90)'

    ! Loop over species
    !$OMP PARALLEL DO                                                        &
    !$OMP DEFAULT( SHARED                                                   )&
    !$OMP PRIVATE( S, N                                                     )
    DO S = 1, SIZE( mapping )

       ! Get the modelId from the mapping array
       N = mapping(S)

       State_Chm%Species(N)%Conc(I,J,L) =                                    &
       State_Chm%Species(N)%Conc(I,J,L) * State_Met%AD(I,J,L)

#ifdef ADJOINT
       IF ( isAdjoint ) THEN
          State_Chm%SpeciesAdj(I,J,L,N) =                                    &
          State_Chm%SpeciesAdj(I,J,L,N) * State_Met%AD(I,J,L)
       ENDIF
#endif

       ! Update units metadata
       State_Chm%Species(N)%Previous_Units = State_Chm%Species(N)%Units
       State_Chm%Species(N)%Units = KG_SPECIES

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
  SUBROUTINE ConvertBox_Kg_to_KgKgDry( I,        J,         L,               &
                                      State_Met, State_Chm, mapping,         &
                                      isAdjoint, RC                         )
!
! !INPUT PARAMETERS:
!
    INTEGER,        INTENT(IN)    :: I, J, L     ! Grid box indexes
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
! !LOCAL VARIABLES:
!
    INTEGER :: N
    CHARACTER(LEN=255) :: errMsg, thisLoc

    !========================================================================
    ! ConvertSpc_Kg_to_KgKgDry begins here!
    !========================================================================

    ! Initialize
    RC      = GC_SUCCESS
    errMsg  = ''
    thisLoc = &
     ' -> at ConvertBox_Kg_to_KgKgDry (in GeosUtil/unitconv_mod.F90)'

    !$OMP PARALLEL DO                                                        &
    !$OMP DEFAULT( SHARED                                                   )&
    !$OMP PRIVATE( S, N                                                     )
    DO S = 1, SIZE( mapping )

       ! Get the modelId from the mapping array
       N = mapping(S)

       State_Chm%Species(N)%Conc(I,J,L) =                                    &
       State_Chm%Species(N)%Conc(I,J,L) / State_Met%AD(I,J,L)

#ifdef ADJOINT
       IF ( isAdjoint ) THEN
          State_Chm%SpeciesAdj(I,J,L,N) =                                    &
          State_Chm%SpeciesAdj(I,J,L,N) / State_Met%AD(I,J,L)
       ENDIF
#endif

       ! Update units metadata
       State_Chm%Species(N)%Previous_Units = State_Chm%Species(N)%Units
       State_Chm%Species(N)%Units = KG_SPECIES_PER_KG_DRY_AIR

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
  SUBROUTINE ConvertBox_Kgm2_to_Kg( I,         J,          L,                &
                                    State_Chm, State_Grid, mapping,          &
                                    isAdjoint, RC                           )
!
! !INPUT PARAMETERS:
!
    INTEGER,        INTENT(IN)    :: I, J, L     ! Grid box indexes
    TYPE(GrdState), INTENT(IN)    :: State_Grid  ! Grid State object
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
    INTEGER            :: N,      S
    CHARACTER(LEN=255) :: errMsg, thisLoc

    !========================================================================
    ! ConvertBox_Kgm2_to_Kg begins here!
    !========================================================================

    ! Initialize
    RC      = GC_SUCCESS
    errMsg  = ''
    thisLoc = &
     ' -> at ConvertBox_Kgm2_to_Kg (in GeosUtil/unitconv_mod.F90)'

    ! Loop over species
    !$OMP PARALLEL DO                                                        &
    !$OMP DEFAULT( SHARED                                                   )&
    !$OMP PRIVATE( N, S                                                     )
    DO S = 1, SIZE( mapping )

       ! Get the modelId from the mapping array
       N = mapping(S)

       State_Chm%Species(N)%Conc(I,J,L) =                                    &
       State_Chm%Species(N)%Conc(I,J,L) * State_Grid%Area_M2(I,J)

#ifdef ADJOINT
       IF ( isAdjoint ) THEN
          State_Chm%SpeciesAdj(I,J,L,N) =                                    &
          State_Chm%SpeciesAdj(I,J,L,N) * State_Grid%Area_M2(I,J)
       ENDIF
#endif

       ! Update units metadata
       State_Chm%Species(N)%Previous_Units = State_Chm%Species(N)%Units
       State_Chm%Species(N)%Units = KG_SPECIES

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
  SUBROUTINE ConvertBox_Kg_to_Kgm2( I,         J,          L,                &
                                    State_Chm, State_Grid, mapping,          &
                                    isAdjoint, RC                           )
!
! !INPUT PARAMETERS:
!
    INTEGER,        INTENT(IN)    :: I, J, L     ! Grid box indexes
    TYPE(GrdState), INTENT(IN)    :: State_Grid  ! Grid State object
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
    INTEGER            :: N,      S
    CHARACTER(LEN=255) :: errMsg, thisLoc

    !========================================================================
    ! ConvertBox_Kg_to_Kgm2 begins here!
    !========================================================================

    ! Initialize
    RC      = GC_SUCCESS
    errMsg  = ''
    thisLoc = &
     ' -> at ConvertBox_Kg_to_Kgm2 (in GeosUtil/unitconv_mod.F90)'

    ! Loop over species
    !$OMP PARALLEL DO                                                        &
    !$OMP DEFAULT( SHARED                                                   )&
    !$OMP PRIVATE( N, S                                                     )
    DO S = 1, SIZE( mapping )

       ! Get the modelId from the mapping array
       N = mapping(S)

       State_Chm%Species(N)%Conc(I,J,L) =                                    &
       State_Chm%Species(N)%Conc(I,J,L) / State_Grid%Area_M2(I,J)

#ifdef ADJOINT
       IF ( isAdjoint ) THEN
          State_Chm%SpeciesAdj(I,J,L,N) =                                    &
          State_Chm%SpeciesAdj(I,J,L,N) / State_Grid%Area_M2(I,J)
       ENDIF
#endif

       ! Update units metadata
       State_Chm%Species(N)%Previous_Units = State_Chm%Species(N)%Units
       State_Chm%Species(N)%Units = KG_SPECIES_PER_M2

    ENDDO
    !$OMP END PARALLEL DO

  END SUBROUTINE ConvertBox_Kg_to_Kgm2
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Check_Units
!
! !DESCRIPTION: Ensures that all species have the same Units value before
!  unit conversion is done.
!\\
! !INTERFACE:
!
  SUBROUTINE Check_Units( State_Chm, mapping, unit, RC )
!
! !INPUT PARAMETERS:
!
    INTEGER,        INTENT(IN)    :: mapping(:)  ! Mapping to species Id
    INTEGER,        INTENT(IN)    :: unit        ! Input unit flag
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure
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
    ! Scalars
    INTEGER            :: N,      S

    ! Strings
    CHARACTER(LEN=255) :: errMsg, thisLoc

    !========================================================================
    ! Check_Units begins here!
    !========================================================================

    ! Initialize
    RC      = GC_SUCCESS
    errMsg  = ''
    thisLoc = ' -> at Check_Units (in module GeosUtil/unitconv_mod.F90)'

    ! Make sure all species start with the proper unit, or throw an error
    IF ( .not. ALL( State_Chm%Species(mapping)%Units == unit ) ) THEN
       errMsg = 'All species do not have Units = '                        // &
                 TRIM( UNIT_STR( unit ) )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

  END SUBROUTINE Check_Units
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Check_Previous_Units
!
! !DESCRIPTION:  Ensures that all species have the same Previous_Units value
!  before unit conversion is done.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Check_Previous_Units( State_Chm, mapping, unit, RC )
!
! !INPUT PARAMETERS:
!
    INTEGER,        INTENT(IN)    :: mapping(:)  ! Mapping to species Id
    INTEGER,        INTENT(IN)    :: unit        ! Input unit flag
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure
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
    ! Scalars
    INTEGER            :: N,      S

    ! Strings
    CHARACTER(LEN=255) :: errMsg, thisLoc

    !========================================================================
    ! Check_Previous_Units begins here!
    !========================================================================

    ! Initialize
    RC      = GC_SUCCESS
    errMsg  = ''
    thisLoc = &
     ' -> at Check_PreviousUnits (in module GeosUtil/unitconv_mod.F90)'

    ! Make sure all species start with the proper unit, or throw an error
    IF ( .not. ALL( State_Chm%Species(mapping)%Previous_Units == unit ) ) THEN
       errMsg = 'All species do not have Previous_Units = '               // &
                TRIM( UNIT_STR( unit ) )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

  END SUBROUTINE Check_Previous_Units
!EOC
END MODULE UnitConv_Mod
