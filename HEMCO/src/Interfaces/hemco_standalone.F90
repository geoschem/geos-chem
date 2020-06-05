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
  USE HCO_Error_Mod
  USE HCOI_StandAlone_Mod, ONLY : HCOI_StandAlone_Run

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
  ! Scalars
  LOGICAL            :: IsDryRun
  LOGICAL            :: Exists
  INTEGER            :: ArgLen
  INTEGER            :: nArg
  INTEGER            :: RC

  ! Strings
  CHARACTER(LEN=63)  :: ProgramName
  CHARACTER(LEN=255) :: ArgVal
  CHARACTER(LEN=255) :: ConfigFile

  !=========================================================================
  ! Initialize
  !=========================================================================
  RC         = HCO_SUCCESS
  nArg       = 0
  ArgLen     = 0
  ArgVal     = ''
  ConfigFile = ''
  IsDryRun   = .FALSE.

  !=========================================================================
  ! The first argument is always the name of the program, so skip ahead
  !=========================================================================
  nArg = nArg + 1
  CALL Get_Command_Argument( nArg, ArgVal, ArgLen )

  !=========================================================================
  ! Parse remaining arguments to determine if this is a dry-run
  !=========================================================================
  DO

     ! Initialize for next argument
     ArgLen  = 0
     ArgVal  = ''

     ! Get the next argument
     CALL Get_Command_Argument( nArg, ArgVal, ArgLen )
     IF ( ArgLen == 0 ) EXIT

     ! Parse the arguments
     SELECT CASE( TRIM( ArgVal ) )

        ! Test for the configuration file
        CASE( '--config-file', '--config', '-c' )

           ! Look for the config file following the the delimiter.
           ! Otherwise use the default dry run logfile name.
           ! Error check for bad input.
           nArg = nArg + 1
           CALL Get_Command_Argument( nArg, ArgVal, ArgLen )
           SELECT CASE( TRIM( ArgVal ) )
              CASE( '--dryrun' )
                 WRITE( 6, '(a)' ) 'Please provide a HEMCO configuration ' // &
                       'file or remove config file argument to use the '   // &
                       'default HEMCO_sa_Config.rc.'
                 STOP
              CASE DEFAULT
                 IF ( ArgLen > 0 ) THEN
                    ConfigFile = TRIM( ArgVal )
                 ENDIF
           END SELECT

        ! Test for the dry-run switch
        CASE( '--dryrun' )
           IsDryRun  = .TRUE.

        CASE DEFAULT
           ! Pass

     END SELECT

     ! Increment the argument counter
     nArg = nArg + 1
  ENDDO

  ! If no configuration file is provided, use HEMCO_sa_Config.rc
  IF ( LEN_TRIM( ConfigFile ) == 0 ) THEN
     ConfigFile='HEMCO_sa_Config.rc'
  ENDIF

  ! Test if the config file exists
  INQUIRE( FILE=TRIM( ConfigFile ), EXIST=Exists )

  ! Test if the file exists and define an output string
  IF ( .not. Exists ) THEN
     ! Write message to stdout
     WRITE( 6, '(a)' ) TRIM( ConfigFile ) // &
                       TRIM( ' NOT FOUND! Check the path of your HEMCO ' // &
                             'configuration file. (hemco_standalone.F90)' )
     STOP
  ENDIF

  !=========================================================================
  ! Run HEMCO in standalone configuration for the given parameters
  ! as specified in the various user-defined input files
  !=========================================================================
  CALL HCOI_StandAlone_Run( ConfigFile = TRIM( ConfigFile ),                 &
                            IsDryRun   = IsDryRun,                           &
                            RC         = RC                                 )

  ! Trap potential errors
  IF ( RC /= HCO_SUCCESS ) THEN
     WRITE( 6, '(a)' ) 'HEMCO_STANDALONE EXITED WITH ERROR!'
     STOP
  ENDIF

  ! End of simulation
  WRITE( 6, '(a)' ) 'HEMCO_STANDALONE FINISHED!'
!EOC
END PROGRAM HEMCO_StandAlone
