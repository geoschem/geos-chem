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
END MODULE Print_Mod
