!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: kpp_standalone_interface.F90
!
! !DESCRIPTION: Contains routines to print the full chemical state in fullchem, which can be used as input to the KPP Standalone.
!\\
!\\
! !INTERFACE:
!
MODULE KPP_Standalone_Interface
!
! !USES:
!
  USE PRECISION_MOD             ! For GEOS-Chem Precision (fp)
  USE HCO_ERROR_MOD             ! For real precisions (hp)
  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBERS:
!
  PUBLIC  :: Check_Domain
  PUBLIC  :: Check_ActiveCell
  PUBLIC  :: Config_KPP_Standalone
  PUBLIC  :: Write_Samples
  PUBLIC  :: Cleanup_KPP_Standalone
!
! !DERIVED TYPES:
!
  TYPE, PRIVATE :: KPP_Standalone_Interface_Type
     ! Scalars
     INTEGER :: NLOC
     LOGICAL :: Active_Cell
     LOGICAL :: SkipIt
     
     ! Strings
     CHARACTER(LEN=255) :: Active_Cell_Name
     CHARACTER(LEN=255) :: Output_Directory

     ! Allocatable arrays
     CHARACTER(LEN=255), DIMENSION(:), ALLOCATABLE :: LocationName
     REAL(hp),           DIMENSION(:), ALLOCATABLE :: LocationLons
     REAL(hp),           DIMENSION(:), ALLOCATABLE :: LocationLats
     INTEGER,            DIMENSION(:), ALLOCATABLE :: IDX
     INTEGER,            DIMENSION(:), ALLOCATABLE :: JDX
     INTEGER,            DIMENSION(:), ALLOCATABLE :: Levels
  END TYPE KPP_Standalone_Interface_Type

  TYPE, PRIVATE :: KPP_Standalone_ActiveCell_Type
     ! Scalars
     LOGICAL            :: Active_Cell
     CHARACTER(LEN=255) :: Active_Cell_Name
  END TYPE KPP_Standalone_ActiveCell_Type
!
! !PRIVATE DATA MEMBERS:
!
  TYPE(KPP_Standalone_Interface_Type),  PRIVATE :: KPP_Standalone_YAML
  TYPE(KPP_Standalone_ActiveCell_Type), PRIVATE :: KPP_Standalone_ActiveCell
  !$OMP THREADPRIVATE( KPP_Standalone_ActiveCell )

! !REVISION HISTORY:
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
! !IROUTINE: check_domain
!
! !DESCRIPTION: Subroutine Check_Domain is used to identify if a
! specified latitude and longitude falls within a grid cell on the
! current CPU. Multiple lat/lon pairs can be checked simultaneously.
! Obin Sturm (psturm@usc.edu) 2023/12/29
!\\
!\\
! !INTERFACE:
!
   SUBROUTINE Check_Domain( RC )

! !USES:
     USE HCO_GeoTools_Mod,         ONLY:  HCO_GetHorzIJIndex
     USE HCO_State_GC_Mod,         ONLY : HcoState
     USE HCO_ERROR_MOD             ! For real precisions (hp)
! !OUTPUT PARAMETERS
    integer, intent(out)       :: RC
   
    
    ! Early exit if no locations
    IF ( KPP_Standalone_YAML%SkipIt ) THEN
       RETURN
    END IF

    CALL HCO_GetHorzIJIndex( HcoState,                      &
                             KPP_Standalone_YAML%NLOC,         &
                             KPP_Standalone_YAML%LocationLons, &
                             KPP_Standalone_YAML%LocationLats, &
                             KPP_Standalone_YAML%IDX,          &
                             KPP_Standalone_YAML%JDX,          &
                             RC)
   END SUBROUTINE Check_Domain
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: check_activecell
!
! !DESCRIPTION: Subroutine Check_ActiveCell is used to identify if a grid cell
! is within a specified latitude and longitude to print the full chemical state
! (all concentrations, reaction rates, rate constants, and meteo metadata).
! Obin Sturm (psturm@usc.edu) 2024/03/11
!\\
!\\
! !INTERFACE:
!
   SUBROUTINE Check_ActiveCell( I, J, L, State_Grid )

! !USES:
     USE State_Grid_Mod,   ONLY : GrdState
! !INPUT PARAMETERS:
     INTEGER, INTENT(IN)   :: I,J,L        ! Grid Indices
     TYPE(GrdState), INTENT(IN)     :: State_Grid ! Grid State object
! !LOCAL VARIABLES
    INTEGER                :: K

    ! Early exit if there was no YAML file or no active cells
    IF ( KPP_Standalone_YAML%SkipIt ) RETURN

    KPP_Standalone_ActiveCell%Active_Cell = .FALSE.
    KPP_Standalone_ActiveCell%Active_Cell_Name = ''
    
    IF ( ANY(L == KPP_Standalone_YAML%Levels) ) THEN
       DO K = 1,KPP_Standalone_YAML%NLOC
          IF ( KPP_Standalone_YAML%IDX(K) == I .AND. &
               KPP_Standalone_YAML%JDX(K) == J ) THEN
             KPP_Standalone_ActiveCell%Active_Cell = .TRUE.
             KPP_Standalone_ActiveCell%Active_Cell_Name =  &
                  KPP_Standalone_YAML%LocationName(K)
             !write(*,*) trim(KPP_Standalone_YAML%Active_Cell_Name), " LatLon: " , State_Grid%YMid(I,J), State_Grid%XMid(I,J)
          ENDIF
       ENDDO
    ENDIF
  END SUBROUTINE Check_ActiveCell
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Config_KPP_Standalone
!
! !DESCRIPTION: Subroutine Config_KPP_Standalone reads a set of gridcells to be sampled
! and the full chemical state printed.
! Obin Sturm (psturm@usc.edu) 2024/03/11
!\\
!\\
! !INTERFACE:
!
   SUBROUTINE Config_KPP_Standalone( Input_Opt, RC ) 
      USE QfYaml_Mod
      USE ErrCode_Mod
      USE Input_Opt_Mod, ONLY : OptInput
      USE RoundOff_Mod,  ONLY : Cast_and_RoundOff
      USE inquireMod,    ONLY : findFreeLUN
! !INPUT PARAMETERS:
!
      TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
!
! !OUTPUT PARAMETERS:
!
      INTEGER,        INTENT(OUT)   :: RC          ! Success or failure 
!
! !LOCAL VARIABLES:
!
      ! Scalars
      INTEGER                      :: I, N
      INTEGER                      :: IU_FILE     ! Available unit for writing
      INTEGER                      :: path_exists
      LOGICAL                      :: file_exists
      LOGICAL                      :: v_bool
 
      ! Strings
      CHARACTER(LEN=255)           :: thisLoc
      CHARACTER(LEN=512)           :: errMsg
      CHARACTER(LEN=QFYAML_NamLen) :: key
      CHARACTER(LEN=QFYAML_StrLen) :: v_str

      ! Objects
      TYPE(QFYAML_t)               :: Config, ConfigAnchored

      ! Arrays
      INTEGER                      :: a_int(QFYAML_MaxArr)

      ! String arrays
      CHARACTER(LEN=QFYAML_NamLen) :: a_str(QFYAML_MaxArr)
  
      ! YAML configuration file name to be read
      CHARACTER(LEN=30), PARAMETER :: configFile = &
           './kpp_standalone_interface.yml'
  
      ! Inquire if YAML interface exists -- if not, skip initializing 
      KPP_Standalone_YAML%SkipIt = .FALSE.
      INQUIRE( FILE=configFile, EXIST=file_exists ) 
      IF ( .NOT. file_exists ) THEN
         KPP_Standalone_YAML%SkipIt = .TRUE.
         IF ( Input_Opt%amIRoot ) &
            write(*,*) "Config file ", configFile, " not found, skipping KPP Standalone interface"
         RETURN
      END IF
      
      ! Assume success
      RC      = GC_SUCCESS
      errMsg  = ''
      thisLoc = ' -> at Config_KPP_Standalone (in module GeosCore/kpp_standalone_interface.F90)'     
 
      !========================================================================
      ! Read the YAML file into the Config object
      !========================================================================
      CALL QFYAML_Init( configFile, Config, ConfigAnchored, RC )
      IF ( RC /= GC_SUCCESS ) THEN
         errMsg = 'Error reading configuration file: ' // TRIM( configFile )
         CALL GC_Error( errMsg, RC, thisLoc )
         RETURN
      ENDIF

      !========================================================================
      ! Read the main on/off switch; Exit if the switch is turned off
      !========================================================================
      key = "settings%activate"
      v_bool = MISSING_BOOL
      CALL QFYAML_Add_Get( Config, TRIM( key ), v_bool, "", RC )
      IF ( RC /= GC_SUCCESS ) THEN
         errMsg = 'Error parsing ' // TRIM( key ) // '!'
         CALL GC_Error( errMsg, RC, thisLoc )
         RETURN
      ENDIF
      KPP_Standalone_YAML%SkipIt = ( .not. v_bool )

      !========================================================================
      ! Read the list of active cells
      !========================================================================
      key = "active_cells"
      a_str = MISSING_STR
      CALL QFYAML_Add_Get( Config, key, a_str, "", RC, dynamic_size=.TRUE. )
      IF ( RC /= GC_SUCCESS ) THEN
         errMsg = 'Error parsing ' // TRIM( key ) // '!'
         CALL GC_Error( errMsg, RC, thisLoc )
         RETURN
      ENDIF

      !========================================================================
      ! Get the number of active cells (if 0, return) and the list of names
      !========================================================================
      KPP_Standalone_YAML%NLOC = Find_Number_of_Locations( a_str )
      IF ( KPP_Standalone_YAML%NLOC .eq. 0 ) THEN
         ! Set SkipIt flag to short circuit other subroutines
         KPP_Standalone_YAML%SkipIt = .TRUE.
         IF ( Input_Opt%amIRoot ) &
            write(*,*) "No active cells for box modeling in kpp_standalone_interface.yml"
         RETURN
      END IF
      ALLOCATE( KPP_Standalone_YAML%LocationName( KPP_Standalone_YAML%NLOC ), STAT=RC )
      CALL GC_CheckVar( 'KPP_Standalone_YAML%LocationName', 0, RC )
      IF ( RC /= GC_SUCCESS ) RETURN
      DO I = 1,KPP_Standalone_YAML%NLOC
         KPP_Standalone_YAML%LocationName(I) = TRIM( a_str(I) )
         print*, trim(KPP_Standalone_YAML%LocationName(I))
      END DO

      !========================================================================
      ! Read latitude and longitude of active cells
      !========================================================================

      ! Allocate number of locations for lats and lons
      ALLOCATE( KPP_Standalone_YAML%LocationLons( KPP_Standalone_YAML%NLOC ), STAT=RC )
      CALL GC_CheckVar( 'KPP_Standalone_YAML%LocationLons', 0, RC )
      IF ( RC /= GC_SUCCESS ) RETURN
      
      ALLOCATE( KPP_Standalone_YAML%LocationLats( KPP_Standalone_YAML%NLOC ), STAT=RC )
      CALL GC_CheckVar( 'KPP_Standalone_YAML%LocationLats', 0, RC )
      IF ( RC /= GC_SUCCESS ) RETURN
     
      ! Read coordinates
      DO I = 1,KPP_Standalone_YAML%NLOC
      ! Read longitudes
         key = "locations%"//TRIM( KPP_Standalone_YAML%LocationName(I) )//"%longitude"
         v_str = MISSING_STR
         CALL QFYAML_Add_Get( Config, TRIM( key ), v_str, "", RC )
         IF ( RC /= GC_SUCCESS ) THEN
            errMsg = 'Error parsing ' // TRIM( key ) // '!'
            CALL GC_Error( errMsg, RC, thisLoc )
            RETURN
         ENDIF
         KPP_Standalone_YAML%LocationLons( I ) = Cast_and_RoundOff( TRIM( v_str ), places=-1 )
      ! Read latitudes
         key = "locations%"//TRIM( KPP_Standalone_YAML%LocationName(I) )//"%latitude"
         v_str = MISSING_STR
         CALL QFYAML_Add_Get( Config, TRIM( key ), v_str, "", RC )
         IF ( RC /= GC_SUCCESS ) THEN
            errMsg = 'Error parsing ' // TRIM( key ) // '!'
            CALL GC_Error( errMsg, RC, thisLoc )
            RETURN
         ENDIF
         KPP_Standalone_YAML%LocationLats( I ) = Cast_and_RoundOff( TRIM( v_str ), places=-1 ) 
      END DO
      
      ! Allocate IDX and JDX (masks for whether a location is on the CPU)
      ALLOCATE( KPP_Standalone_YAML%IDX( KPP_Standalone_YAML%NLOC ), STAT=RC )
      CALL GC_CheckVar( 'KPP_Standalone_YAML%IDX', 0, RC )
      IF ( RC /= GC_SUCCESS ) RETURN

      ALLOCATE( KPP_Standalone_YAML%JDX( KPP_Standalone_YAML%NLOC ), STAT=RC )
      CALL GC_CheckVar( 'KPP_Standalone_YAML%JDX', 0, RC )
      IF ( RC /= GC_SUCCESS ) RETURN

      KPP_Standalone_YAML%IDX(:) = -1
      KPP_Standalone_YAML%JDX(:) = -1

      !========================================================================
      ! Get the list of levels and number of levels
      !========================================================================
      ! Note: could add capability for location specific levels
      key = "settings%levels"
      a_int = MISSING_INT
      CALL QFYAML_Add_Get( Config, key, a_int, "", RC, dynamic_size=.TRUE. )
      IF ( RC /= GC_SUCCESS ) THEN
         errMsg = 'Error parsing ' // TRIM( key ) // '!'
         CALL GC_Error( errMsg, RC, thisLoc )
         RETURN
      ENDIF
      N = Find_Number_of_Levels( a_int )
      ! if no specified levels, print the surface
      IF ( N .eq. 0 ) THEN
         N = 1
         a_int(1) = 1
      END IF
      ALLOCATE( KPP_Standalone_YAML%Levels( N ), STAT=RC )
      CALL GC_CheckVar( 'KPP_Standalone_YAML%Levels', 0, RC )
      IF ( RC /= GC_SUCCESS ) RETURN
      DO I = 1,N
         KPP_Standalone_YAML%Levels(I) = a_int(I)
      END DO

      !========================================================================
      ! Set the output directory
      !========================================================================
      ! Get that value
      key = "settings%output_directory"
      v_str = MISSING_STR
      CALL QFYAML_Add_Get( Config, TRIM( key ), v_str, "", RC )
      IF ( RC /= GC_SUCCESS ) THEN
         errMsg = 'Error parsing ' // TRIM( key ) // '!'
         CALL GC_Error( errMsg, RC, thisLoc )
         RETURN
      ENDIF
      ! Check to see if the directory exists
      ! Do this in a portable way that works across compilers
      ! The directory specifier in inquire might be specific to ifort
      ! So instead try to open a test file within the output directory
      ! Check ./OutputDir (which exists for GEOS-Chem and GCHP) as backup
      IU_FILE = findFreeLUN()
      open(IU_FILE,FILE=trim(v_str)//'/.test_directory_existence', &
           action = "WRITE",iostat=path_exists,access='SEQUENTIAL')

      ! If the specified folder doesn't exist, try OutputDir
      IF ( path_exists /= 0 ) THEN 
         open(IU_FILE,FILE='./OutputDir/.test_directory_existence', &
              action = "WRITE",iostat=path_exists,access='SEQUENTIAL')
         KPP_Standalone_YAML%Output_Directory = "./OutputDir"
         IF ( Input_Opt%amIRoot )                                                        &
             write(*,*) "KPP Standalone Interface warning: Specified output directory ", &
                         trim(v_str), " was not found, trying default output path './OutputDir' "
         ! If OutputDir doesn't exist, write to the current directory
         IF ( (path_exists /= 0) ) THEN
            IF ( Input_Opt%amIRoot )                                                       &
               write(*,*) "KPP Standalone Interface warning: Specified output directory ", &
                           trim(v_str), " and default output directory './OutputDir' " //  &
                           "were not found, writing output to the current directory './'" 
            KPP_Standalone_YAML%Output_Directory = "./"
         ENDIF
      ELSE 
         KPP_Standalone_YAML%Output_Directory = trim(v_str)
         close(IU_FILE)
      END IF
       
   END SUBROUTINE Config_KPP_Standalone
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Write_Samples
!
! !DESCRIPTION: Subroutine Write_Samples writes the full chemical state
! (concentrations, reaction rates and rate constants, meteorological conditions).
! Obin Sturm (psturm@usc.edu) 2024/03/11
!\\
!\\
! !INTERFACE:
!
   SUBROUTINE Write_Samples( I, J, L, initC, localRCONST, initHvalue, exitHvalue,  &
                             State_Grid,   State_Chm, State_Met, Input_Opt,        &
                             KPP_TotSteps, RC, FORCE_WRITE, CELL_NAME )
      USE ErrCode_Mod
      USE State_Grid_Mod,           ONLY : GrdState
      USE State_Chm_Mod,            ONLY : ChmState
      USE State_Met_Mod,            ONLY : MetState
      USE Input_Opt_Mod,            ONLY : OptInput
      USE GcKpp_Function
      USE GcKpp_Parameters,         ONLY : NSPEC, NREACT, NVAR 
      USE TIME_MOD,                 ONLY : GET_TS_CHEM
      USE TIME_MOD,                 ONLY : TIMESTAMP_STRING
      USE TIME_MOD,                 ONLY : Get_Minute
      USE TIME_MOD,                 ONLY : Get_Hour
      USE TIME_MOD,                 ONLY : Get_Day
      USE TIME_MOD,                 ONLY : Get_Month
      USE TIME_MOD,                 ONLY : Get_Year
      USE Pressure_Mod,             ONLY : Get_Pcenter
      USE inquireMod,               ONLY : findFreeLUN
! !INPUT PARAMETERS:
!
      INTEGER,        INTENT(IN)    :: I            ! Longitude index
      INTEGER,        INTENT(IN)    :: J            ! Latitude index
      INTEGER,        INTENT(IN)    :: L            ! GEOS-Chem vertical level
      INTEGER,        INTENT(IN)    :: KPP_TotSteps ! Total KPP integrator steps

      TYPE(GrdState), INTENT(IN)    :: State_Grid                 ! Grid State object
      TYPE(ChmState), INTENT(IN)    :: State_Chm                  ! Chemistry State object
      TYPE(MetState), INTENT(IN)    :: State_Met                  ! Meteorology State object
      TYPE(OptInput), INTENT(IN)    :: Input_Opt                  ! Input Options object
      REAL(dp), INTENT(IN)          :: initC(NSPEC)               ! Initial concentrations
      REAL(dp), INTENT(IN)          :: localRCONST(NREACT)        ! Rate constants
      REAL(dp)                      :: initHvalue                 ! Initial timestep
      REAL(dp)                      :: exitHvalue                 ! Final timestep, RSTATE(Nhexit)

! !OPTIONAL INPUT PARAMETER
      LOGICAL, OPTIONAL             :: FORCE_WRITE  ! Write even if not in an active cell
      CHARACTER(LEN=*), OPTIONAL    :: CELL_NAME    ! Customize the name of this file 
!
! !AUXILLIARY LOCAL PARAMETERS (pass the aux bc Fortran doesn't have defaults for kwargs)
      LOGICAL            :: FORCE_WRITE_AUX   ! Write even if not in an active cell
      CHARACTER(LEN=255) :: CELL_NAME_AUX     ! Customize the name of this file 
!
! !OUTPUT PARAMETERS:
!
      INTEGER,        INTENT(OUT)   :: RC          ! Success or failure
!
! !LOCAL VARIABLES:
      ! Integers
      INTEGER                :: N          ! Loop index
      INTEGER                :: IU_FILE    ! Available unit for writing
      INTEGER                :: SpcID      ! Mapping from State_Chm and KPP
      REAL(fp)               :: DT         ! Chemistry operator timestep
      ! Strings
      CHARACTER(LEN=255)     :: YYYYMMDD_hhmmz
      CHARACTER(LEN=255)     :: level_string
      CHARACTER(LEN=512)     :: errMsg, filename
      
      ! Arrays
      REAL(dp)               :: Vloc(NVAR), Aout(NREACT)  ! For KPP reaction rate diagnostics


      ! Did a user want to write the chemical state even if 
      ! not in an active cell?
      FORCE_WRITE_AUX = .FALSE.
      IF ( PRESENT( FORCE_WRITE ) ) FORCE_WRITE_AUX = FORCE_WRITE

      ! Quit early if there's no writing to be done
      IF ( .not. KPP_Standalone_ActiveCell%Active_Cell  .AND.               &
           .not. FORCE_WRITE_AUX                       ) THEN
         RETURN
      END IF

      ! Did the call include an optional cell name?
      CELL_NAME_AUX = ''
      IF ( PRESENT( CELL_NAME ) ) CELL_NAME_AUX = CELL_NAME

      ! Get KPP state
      CALL Fun( V       = initC(1:NVAR),                                     &
                F       = initC(NVAR+1:NSPEC),                               &
                RCT     = localRCONST,                                       &
                Vdot    = Vloc,                                              &
                Aout    = Aout                                          )

      DT = GET_TS_CHEM()

      !========================================================================
      ! Write the file
      !========================================================================

      ! Find a free file LUN
      IU_FILE = findFreeLUN()
      write(level_string,'(I0)') L
      write(YYYYMMDD_hhmmz,'(I0.4,I0.2,I0.2,a,I0.2,I0.2)' ) &
            Get_Year(), Get_Month(), Get_Day(),'_', Get_Hour(), Get_Minute()

      ! Filename for output
      filename = TRIM( KPP_Standalone_YAML%Output_Directory             ) // &
                 '/'                                                      // &
                 TRIM( Cell_Name_Aux                                    ) // &
                 TRIM( KPP_Standalone_ActiveCell%Active_Cell_Name       ) // &
                 '_L'                                                     // &
                 trim( level_string                                     ) // &
                 '_'                                                      // &
                 TRIM( YYYYMMDD_hhmmz                                   ) // &
                 '.txt'

      ! Open the file
      open( IU_FILE,  FILE=TRIM(filename), ACTION="WRITE",                   &
                      IOSTAT=RC,           ACCESS='SEQUENTIAL')
      IF ( RC /= GC_SUCCESS ) THEN
         errMsg = 'Error writing chemical state to KPP Standalone file' 
         CALL GC_Error( errMsg, RC, '' )
         RETURN
      ENDIF

      ! Write header to file
      write(IU_FILE, '(a)') '48                                                                         '
      write(IU_FILE, '(a)') '==========================================================================='
      write(IU_FILE, '(a)') '                                                                           '
      write(IU_FILE, '(a)') '                  KPP Standalone Atmospheric Chemical State                '
      write(IU_FILE, '(a)') 'File Description:                                                          '
      write(IU_FILE, '(a)') 'This file contains model output of the atmospheric chemical state          '
      write(IU_FILE, '(a)') 'as simulated by the GEOS-Chem chemistry module in a 3D setting.            '
      write(IU_FILE, '(a)') 'Each grid cell represents the chemical state of an individual location,    '
      write(IU_FILE, '(a)') 'suitable for input into a separate KPP Standalone program which will       '
      write(IU_FILE, '(a)') 'replicate the chemical evolution of that grid cell for mechanism analysis. '
      write(IU_FILE, '(a)') 'Note that the KPP Standalone will only use concentrations, rate constants, '
      write(IU_FILE, '(a)') 'and KPP-specific fields. All other fields are for reference. The first line'
      write(IU_FILE, '(a)') 'contains the number of lines in this header. If wanting to use this output '
      write(IU_FILE, '(a)') 'for other analysis, a Python class to read these fields is available by    '
      write(IU_FILE, '(a)') 'request, contact Obin Sturm (psturm@usc.edu).                              '
      write(IU_FILE, '(a)') '                                                                           '
      write(IU_FILE, '(a)') 'Generated by the GEOS-Chem Model                                           '
      write(IU_FILE, '(a)') '       (https://geos-chem.org/)                                            '
      write(IU_FILE, '(a)') 'Using the KPP Standalone Interface                                         '
      write(IU_FILE, '(a)') 'github.com/GEOS-ESM/geos-chem/tree/feature/psturm/kpp_standalone_interface '
      write(IU_FILE, '(a)') '     With contributions from:                                              '
      write(IU_FILE, '(a)') '        Obin Sturm (psturm@usc.edu)                                        '
      write(IU_FILE, '(a)') '        Christoph Keller                                                   '
      write(IU_FILE, '(a)') '        Michael Long                                                       '
      write(IU_FILE, '(a)') '        Sam Silva                                                          '
      write(IU_FILE, '(a)') '                                                                           '
      ! Write the grid cell metadata as part of the header
      write(IU_FILE,'(a)'      ) 'Meteorological and general grid cell metadata    '
      write(IU_FILE,'(a,a)'    ) 'Location:                                        ', trim(CELL_NAME_AUX)//trim(KPP_Standalone_YAML%ACTIVE_CELL_NAME)
      write(IU_FILE,'(a,a)'    ) 'Timestamp:                                       ', TIMESTAMP_STRING()
      write(IU_FILE,'(a,F11.4)') 'Longitude (degrees):                             ', State_Grid%XMid(I,J)
      write(IU_FILE,'(a,F11.4)') 'Latitude (degrees):                              ', State_Grid%YMid(I,J)
      write(IU_FILE,'(a,i6)'   ) 'GEOS-Chem Vertical Level:                        ', L
      write(IU_FILE,'(a,F11.4)') 'Pressure (hPa):                                  ', Get_Pcenter(I,J,L)
      write(IU_FILE,'(a,F11.2)') 'Temperature (K):                                 ', State_Met%T(I,J,L)
      write(IU_FILE,'(a,e11.4)') 'Dry air density (molec/cm3):                     ', State_Met%AIRNUMDEN(I,J,L)
      write(IU_FILE,'(a,e11.4)') 'Water vapor mixing ratio (vol H2O/vol dry air):  ', State_Met%AVGW(I,J,L)
      write(IU_FILE,'(a,e11.4)') 'Cloud fraction:                                  ', State_Met%CLDF(I,J,L)
      write(IU_FILE,'(a,e11.4)') 'Cosine of solar zenith angle:                    ', State_Met%SUNCOSmid(I,J)
      write(IU_FILE,'(a)'      ) 'KPP Integrator-specific parameters               '
      write(IU_FILE,'(a,F11.4)') 'Init KPP Timestep (seconds):                     ', initHvalue
      write(IU_FILE,'(a,F11.4)') 'Exit KPP Timestep (seconds):                     ', exitHvalue
      write(IU_FILE,'(a,F11.4)') 'Chemistry operator timestep (seconds):           ', DT
      write(IU_FILE,'(a,i6)'   ) 'Number of internal timesteps:                    ', KPP_TotSteps
      write(IU_File,'(a)'      ) 'CSV data of full chemical state, including species concentrations,    '
      write(IU_File,'(a)'      ) 'rate constants (R) and instantaneous reaction rates (A).              '
      write(IU_File,'(a)'      ) 'All concentration units are in molecules/cc and rates in molec/cc/s.  '
      write(IU_FILE, '(a)') '                                                                           '
      write(IU_FILE, '(a)') '==========================================================================='
      write(IU_FILE, '(a)') 'Name,   Value                                                              '
      DO N=1,NSPEC
         SpcID = State_Chm%Map_KppSpc(N)
         IF ( SpcID <= 0 ) THEN
            write(IU_FILE,'(A,I0,A,E25.16E3)') "C",N,",",initC(N)
            CYCLE
         ENDIF
         write(IU_FILE,'(A,A,E25.16E3)') trim(State_Chm%SpcData(SpcID)%Info%Name),',',initC(N)
      ENDDO
      DO N=1,NREACT
         write(IU_FILE,'(A,I0,A,E25.16E3)') 'R',N,',', localRCONST(N)
      ENDDO
      DO N=1,NREACT
         write(IU_FILE,'(A,I0,A,E25.16E3)') 'A',N,',', Aout(N)
      ENDDO
      close(IU_FILE)

   END SUBROUTINE Write_Samples
!EOC
! !INPUT PARAMETERS:
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: cleanup_kpp_standalone
!
! !DESCRIPTION: Deallocates module variables that may have been allocated at run time
!               and unnecessary files required during the process
!\\
!\\
! !INTERFACE:
   SUBROUTINE Cleanup_KPP_Standalone( RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE inquireMod,      ONLY : findFreeLUN
!
! !OUTPUT PARAMETERS:
!
    INTEGER, INTENT(OUT) :: RC          ! Success or failure?
!
! !REVISION HISTORY:
!  11 Mar 2024 - Obin Sturm - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
    ! Assume success
    RC = GC_SUCCESS

    IF ( ALLOCATED( KPP_Standalone_YAML%LocationName ) ) THEN
       DEALLOCATE( KPP_Standalone_YAML%LocationName, STAT=RC  )
       CALL GC_CheckVar( 'kpp_standalone_interface.F90:KPP_Standalone_YAML%LocationName', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ALLOCATED( KPP_Standalone_YAML%LocationLons ) ) THEN
       DEALLOCATE( KPP_Standalone_YAML%LocationLons, STAT=RC  )
       CALL GC_CheckVar( 'kpp_standalone_interface.F90:KPP_Standalone_YAML%LocationLons', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ALLOCATED( KPP_Standalone_YAML%LocationLats ) ) THEN
       DEALLOCATE( KPP_Standalone_YAML%LocationLats, STAT=RC  )
       CALL GC_CheckVar( 'kpp_standalone_interface.F90:KPP_Standalone_YAML%LocationLats', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ALLOCATED( KPP_Standalone_YAML%IDX ) ) THEN
       DEALLOCATE( KPP_Standalone_YAML%IDX, STAT=RC  )
       CALL GC_CheckVar( 'kpp_standalone_interface.F90:KPP_Standalone_YAML%IDX', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ALLOCATED( KPP_Standalone_YAML%JDX ) ) THEN
       DEALLOCATE( KPP_Standalone_YAML%JDX, STAT=RC  )
       CALL GC_CheckVar( 'kpp_standalone_interface.F90:KPP_Standalone_YAML%JDX', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ALLOCATED( KPP_Standalone_YAML%Levels ) ) THEN
       DEALLOCATE( KPP_Standalone_YAML%Levels, STAT=RC  )
       CALL GC_CheckVar( 'kpp_standalone_interface.F90:KPP_Standalone_YAML%Levels', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF
    
   END SUBROUTINE Cleanup_KPP_Standalone
!EOC
! !INPUT PARAMETERS:
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Find_Number_of_Locations
!
! !DESCRIPTION: Searches a string array containing location names and returns
!  the number of valid locations (i.e. char that do not match MISSING_STR).
!  Assumes all the valid locations will be listed contiguously at the front
!  of the array. Taken from Get_Number_of_Species from input_mod.F90
!\\
!\\
! !INTERFACE:
   FUNCTION Find_Number_of_Locations( a_str ) RESULT( n_valid )
!
! !INPUT PARAMETERS:
!
    CHARACTER(LEN=*), INTENT(IN) :: a_str(:)
!
! !RETURN VALUE:
!
    INTEGER                      :: n_valid
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER :: N

    ! Return the number of valid locations
    n_valid = 0
    DO N = 1, SIZE( a_str )
       IF ( TRIM( a_str(N) ) == MISSING_STR ) EXIT
       n_valid = n_valid + 1
    ENDDO

  END FUNCTION Find_Number_of_Locations
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Find_Number_of_Levels
!
! !DESCRIPTION: Searches an integer array containing location names and returns
!  the number of valid levels (i.e. int that do not match MISSING_INT).
!  Assumes all the valid levels will be listed contiguously at the front
!  of the array. Taken from Get_Number_of_Species from input_mod.F90
!\\
!\\
! !INTERFACE:
   FUNCTION Find_Number_of_Levels( a_int ) RESULT( n_valid )
!
! !INPUT PARAMETERS:
!
    INTEGER, INTENT(IN) :: a_int(:)
!
! !RETURN VALUE:
!
    INTEGER             :: n_valid
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER :: N

    ! Return the number of valid locations
    n_valid = 0
    DO N = 1, SIZE( a_int )
       IF ( a_int(N) == MISSING_INT ) EXIT
       n_valid = n_valid + 1
    ENDDO

  END FUNCTION Find_Number_of_Levels
!EOC
END MODULE KPP_Standalone_Interface
