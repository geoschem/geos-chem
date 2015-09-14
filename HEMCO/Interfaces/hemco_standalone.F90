!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: hemco_standalone.F90 
!
! !DESCRIPTION: Program HEMCO\_StandAlone is the driver for HEMCO in
!  in stand-alone mode. It receives the configuration file name as 
!  input argument and then calls the routine HCOI\_StandAlone\_Run 
!  (from module file hcoi\_standalone\_mod.F90) to repeatedly call HEMCO 
!  for one or more timesteps on the predefined grid.
!\\
!\\
! !INTERFACE:
!
PROGRAM HEMCO_StandAlone
!
! !USES:
!
  USE HCOI_STANDALONE_MOD, ONLY : HCOI_Standalone_Run

  IMPLICIT NONE
  INTRINSIC :: TRIM
!
! !REVISION HISTORY:
!  16 Jul 2014 - R. Yantosca - Initial version
!  21 Jul 2013 - C. Keller   - Now pass configuration file name as argument
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
  CHARACTER(LEN=63)  :: ProgramName
  CHARACTER(LEN=255) :: ConfigFile
  INTEGER            :: NARG

  !------------------------------------------------
  ! HEMCO_StandAlone begins here!
  !------------------------------------------------

  ! Get configuration file (passed as argument).
  NARG = COMMAND_ARGUMENT_COUNT()
  CALL GET_COMMAND_ARGUMENT(0,ProgramName)

  IF ( NARG > 1) THEN
     WRITE(*,*) 'HEMCO_StandAlone takes only 1 argument!'
     STOP
  ENDIF

  IF ( NARG == 0) THEN
     WRITE(*,*) 'Please provide HEMCO configuration file!'
     STOP
  ENDIF

  CALL GET_COMMAND_ARGUMENT(1,ConfigFile)

  ! Run HEMCO in standalone configuration for the given parameters
  ! as specified in the various user-defined input files
  CALL HCOI_StandAlone_Run( TRIM(ConfigFile) )

  ! End of simulation
  WRITE(*,*) 'HEMCO_STANDALONE FINISHED!'

!EOC
END PROGRAM HEMCO_StandAlone
