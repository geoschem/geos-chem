!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: print_mod.F90
!
! !DESCRIPTION: Module PRINT\_MOD contains routines which are used as a
!  general utility to print various quantities to log for debugging
!  or informational purposes.
!\\
!\\
! !INTERFACE:
!
MODULE Print_Mod
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
  USE UnitConv_Mod

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC :: Print_Species_Min_Max_Sum
  PUBLIC :: Print_Species_Global_Mass
  PUBLIC :: Print_Species_Global_Mass_From_VVDry
!
! !REMARKS:
!
!
! !REVISION HISTORY:
!  15 Oct 20324 - E. Lundgren - Initial version
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
! !IROUTINE: Print_Species_Min_Max_Sum
!
! !DESCRIPTION: Subroutine Print\_Species\_Min\_Max\_Sum prints the
!   minimum, maximum, and sum of species concentrations on the root thread
!   to log. The default is to write all species. Arguments can be passed to
!   to specify start index and stop index of State_Chm%Species array to
!   limit species to one species or a consecutive sequence.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Print_Species_Min_Max_Sum( msg, Input_Opt, State_Chm,  &
                                        RC,  nStart,    nStop        )
!
! !INPUT PARAMETERS:
!

    CHARACTER(LEN=*), INTENT(IN)           :: msg       ! Message to print
    TYPE(OptInput),   INTENT(IN)           :: Input_Opt ! Input Options object
    TYPE(ChmState),   INTENT(IN)           :: State_Chm ! Chemistry State object
    INTEGER,          INTENT(IN), OPTIONAL :: nStart    ! Index of 1st species to print
    INTEGER,          INTENT(IN), OPTIONAL :: nStop     ! Index of last species to print
!
! !OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(OUT)   :: RC     ! Success or failure?
!
! !REMARKS:
!
! !REVISION HISTORY:
!  07 Oct 2024 - E. Lundgren - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER            :: N, N_Start, N_Stop

    ! Strings
    CHARACTER(LEN=255) :: errMsg, errLoc, units

    !========================================================================
    ! Print_Species_Min_Max_Sum begins here!
    !========================================================================

    RC     = GC_SUCCESS
    errMsg = ''
    errLoc = ' -> at Print_Species_Min_Max_Sum (in GeosUtil/print_mod.F90)'

    ! Set defaults
    N_START = 1
    N_STOP = State_Chm%nSpecies

    ! Override with optional args
    IF ( PRESENT(nStart) ) N_START = nStart
    IF ( PRESENT(nStop ) ) N_STOP  = nStop

    ! Write to log
    IF ( Input_Opt%amIRoot ) THEN
       WRITE(6,*) TRIM(msg) // ' (' // TRIM(UNIT_STR(State_Chm%Species(1)%Units)) // ')'
       DO N = N_START, N_STOP
          WRITE( 6, 120 ) N, TRIM( State_Chm%SpcData(N)%Info%Name ), &
               MINVAL( State_Chm%Species(N)%Conc(:,:,:) ), &
               MAXVAL( State_Chm%Species(N)%Conc(:,:,:) ), &
               SUM ( State_Chm%Species(N)%Conc(:,:,:) )
       ENDDO
    ENDIF
120 FORMAT( '   Species ', i3, ', ', a8, ': Min = ', es15.9, &
         '  Max = ',es15.9, '  Sum = ',es15.9)

  END SUBROUTINE Print_Species_Min_Max_Sum
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Print_Species_Global_Mass
!
! !DESCRIPTION: Subroutine Print\_Species\_Global\_Mass prints the
!   global sum of species mass in kg on the root thread to log.
!   The default is to write all species. Arguments can be passed to
!   to specify start index and stop index of State_Chm%Species array to
!   limit species to one species or a consecutive sequence.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Print_Species_Global_Mass( msg,        Input_Opt,   State_Chm,  &
                                        State_Met,  State_Grid,  RC,         &
                                        nStart,     nStop                   )
!
! !INPUT PARAMETERS:
!

    CHARACTER(LEN=*), INTENT(IN)           :: msg        ! Message to print
    TYPE(OptInput),   INTENT(IN)           :: Input_Opt  ! Input Options
    TYPE(MetState),   INTENT(IN)           :: State_Met  ! Meteorology State
    TYPE(GrdState),   INTENT(IN)           :: State_Grid ! Grid State object
    INTEGER,          INTENT(IN), OPTIONAL :: nStart     ! Index of 1st spc
    INTEGER,          INTENT(IN), OPTIONAL :: nStop      ! Index of last spc
!
! INPUT/OUTPUT PARAMETERS:
!

    TYPE(ChmState),   INTENT(INOUT)        :: State_Chm  ! Chemistry State
!
! !OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(OUT)          :: RC         ! Success or failure?
!
! !REMARKS:
!  You may see very small numerical differences when comparing results
!  from a simulation where Print_Species_Global_Mass has been called to a
!  simulation where it has not been called.  This is because there is an
!  additional unit conversion in Print_Species_Global_Mass, which may
!  cause numerical noise in the species concentration array.
!
! !REVISION HISTORY:
!  07 Oct 2024 - E. Lundgren - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER            :: N, N_Start, N_Stop, previous_units

    ! Strings
    CHARACTER(LEN=255) :: errMsg, errLoc, units

    !========================================================================
    ! Print_Species_Global_Mass begins here!
    !========================================================================

    RC     = GC_SUCCESS
    errMsg = ''
    errLoc = ' -> at Print_Species_Global_Mass (in GeosUtil/print_mod.F90)'

    ! Set defaults
    N_START = 1
    N_STOP = State_Chm%nSpecies

    ! Override with optional args
    IF ( PRESENT( nStart) ) N_START = nStart
    IF ( PRESENT( nStop ) ) N_STOP  = nStop

    ! Convert species to kg if needed
    CALL Convert_Spc_Units(                                                  &
         Input_Opt      = Input_Opt,                                         &
         State_Chm      = State_Chm,                                         &
         State_Grid     = State_Grid,                                        &
         State_Met      = State_Met,                                         &
         new_units      = KG_SPECIES,                                        &
         previous_units = previous_units,                                    &
         RC             = RC                                                )

    ! Write to log
    IF ( Input_Opt%amIRoot ) THEN
       WRITE( 6, '(a)' ) 'Global mass sum of each species [kg]'
       IF ( LEN_TRIM ( msg ) > 0 ) WRITE( 6, '(a)' ) TRIM( msg )
       DO N = N_START, N_STOP
          WRITE( 6, 130 ) N, TRIM( State_Chm%SpcData(N)%Info%Name ),         &
                             SUM( State_Chm%Species(N)%Conc )
       ENDDO
    ENDIF

    ! Convert species to original units
    CALL Convert_Spc_Units(                                                  &
         Input_Opt      = Input_Opt,                                         &
         State_Chm      = State_Chm,                                         &
         State_Grid     = State_Grid,                                        &
         State_Met      = State_Met,                                         &
         new_units      = previous_units,                                    &
         RC             = RC                                                )

130 FORMAT( '   Species ', i3, ', ', a9, ': Global mass = ', es15.9 )

  END SUBROUTINE Print_Species_Global_Mass
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Print_Species_Global_Mass
!
! !DESCRIPTION: Subroutine Print\_Species\_Global\_Mass prints the
!   global sum of species mass in kg on the root thread to log.
!   The default is to write all species. Arguments can be passed to
!   to specify start index and stop index of State_Chm%Species array to
!   limit species to one species or a consecutive sequence.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Print_Species_Global_Mass_From_VVDry( msg,        Input_Opt,    &
                                                   State_Chm,  State_Met,    &
                                                   State_Grid, RC,           &
                                                   nStart,     nStop        )
!
! !INPUT PARAMETERS:
!

    CHARACTER(LEN=*), INTENT(IN)    :: msg        ! Message to print
    TYPE(OptInput),   INTENT(IN)    :: Input_Opt  ! Input Options object
    TYPE(MetState),   INTENT(IN)    :: State_Met  ! Meteorology State object
    TYPE(GrdState),   INTENT(IN)    :: State_Grid ! Grid State object
    INTEGER,          OPTIONAL      :: nStart     ! Index of 1st species
    INTEGER,          OPTIONAL      :: nStop      ! Index of last species
!
! INPUT/OUTPUT PARAMETERS:
!
    TYPE(ChmState),   INTENT(INOUT) :: State_Chm  ! Chemistry State object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(OUT)   :: RC         ! Success or failure?
!
! !REMARKS:
!  This routine prints species global masses without modifying the
!  State_Chm%Species array.  This should prevent very small numerical
!  differences caused by roundoff.
!
! !REVISION HISTORY:
!  05 Mar 2026 - R. Yantosca - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER            :: N,       N_Start, N_Stop

    ! Arrays
    REAL(fp)           :: conv(State_Grid%NX,State_Grid%NY,State_Grid%NZ)

    !========================================================================
    ! Print_Species_Global_Mass_From_VVDry begins here!
    !========================================================================

    ! Set defaults
    N_START = 1
    N_STOP = State_Chm%nSpecies

    ! Override with optional args
    IF ( PRESENT( nStart ) ) N_START = nStart
    IF ( PRESENT( nStop  ) ) N_STOP  = nStop

    ! Write to log
    IF ( Input_Opt%amIRoot ) THEN
       WRITE( 6, '(a)' ) 'Global mass sum of each species [kg]'
       IF ( LEN_TRIM ( msg ) > 0 ) WRITE( 6, '(a)' ) TRIM( msg )

       ! Loop over species
       DO N = N_START, N_STOP

          ! Conversion from [mol/mol dry] -> [kg]
          conv = State_Met%AD * ( State_Chm%SpcData(N)%Info%MW_g / AIRMW )

          ! Print species mass in kg
          WRITE( 6, 10 ) N, TRIM( State_Chm%SpcData(N)%Info%Name  ),         &
                            SUM( State_Chm%Species(N)%Conc * conv )
       ENDDO
    ENDIF

 10 FORMAT( '   Species ', i3, ', ', a9, ': Global mass = ', es15.9 )

  END SUBROUTINE Print_Species_Global_Mass_From_VVDry
!EOC
END MODULE Print_Mod
