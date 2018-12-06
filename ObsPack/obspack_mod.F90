!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: OBSPACK_MOD
!
! !DESCRIPTION: Module OBSPACK\_MOD contains variables and routines
!  which to sample a GEOS-Chem model simulation for in situ observations
!  contained in an ObsPack file.
!\\
!\\
! !INTERFACE:
!
MODULE ObsPack_Mod
!
! !USES:
!
  USE Precision_Mod    ! For GEOS-Chem Precision (fp)

  IMPLICIT NONE
  PRIVATE

#include "netcdf.inc"
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC  :: ObsPack_Cleanup
  PUBLIC  :: ObsPack_Init
  PUBLIC  :: ObsPack_Sample
  PUBLIC  :: ObsPack_Write_Output
!
! !PRIVATE MEMBER FUNCTIONS:
!
  PRIVATE :: ObsPack_Read_Input
  PRIVATE :: ObsPack_Get_Indices
  PRIVATE :: Seconds_Since_1970
!
! !AUTHOR:
!  Andrew Jacobson (NOAA); andy.jacobson@noaa.gov
!  Andrew Schuh (Colorado State University), aschuh@colostate.edu
!
! !REMARKS:
!  NOTE: All ObsPack variables are saved as fields of State_Diag,
!  which will also facilitate using ObsPack in other model contexts.
!  Right now this has only been validated with GEOS-Chem "Classic", though.
!
!  TODO list
!  Read variables/tracers to output from input file?  (See
!  planeflight_mod.F)
!  
!  Write variable/tracer names to output file?  Compare with
!  user_output_flask.F90
!
!  RC error code passing
!  
!  Output surface height?
!
! !REVISION HISTORY:
!  04 Jun 2015 - A. Jacobson - Adapted from v10.1 planeflight_mod.f, following
!                              similar work done in v9.2 by Andrew Schuh.
!  05 Dec 2018 - R. Yantosca - Implemented into the standard GEOS-Chem code
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
! !IROUTINE: obspack_init
!
! !DESCRIPTION: Subroutine OBSPACK\_INIT reads information from
!  the input ObsPack file in order to initialize the obspack
!  diagnostic.  It allocates memory in the obs array.  If this
!  structure is already allocated, the routine presumes that the
!  current samples need to be written out before the new ObsPack
!  information is read.  In that case, OBSPACK\_WRITE\_OUTPUT is
!  called.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE ObsPack_Init( am_I_Root,  yyyymmdd,   hhmmss,                   &
                           Input_Opt,  State_Diag, RC                       )
!
! !USES:
!
    USE ErrCode_Mod
    USE File_Mod,       ONLY : File_Exists
    USE Input_Opt_Mod,  ONLY : OptInput 
    USE State_Diag_Mod, ONLY : DgnState
    USE TIME_MOD,       ONLY : Expand_Date
    USE TIME_MOD,       ONLY : Get_Time_Ahead
!
! !INPUT PARAMETERS:
!
    LOGICAL,        INTENT(IN)    :: am_I_Root   ! Is this the root core?
    INTEGER,        INTENT(IN)    :: yyyymmdd    ! Current date
    INTEGER,        INTENT(IN)    :: hhmmss      ! Current time
    TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(DgnState), INTENT(INOUT) :: State_Diag   ! Chemistry State object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
!
! !REVISION HISTORY
!  05 Dec 2018 - R. Yantosca - Implemented into the standard GEOS-Chem code
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    LOGICAL            :: prtDebug

    ! Strings
    CHARACTER(LEN=255) :: ErrMsg, ThisLoc

    !=======================================================================
    ! OBSPACK_INIT begins here
    !=======================================================================

    ! Initialize
    RC       =  GC_SUCCESS
    prtDebug = ( am_I_Root .and. Input_Opt%LPRT )
    ErrMsg   = ''
    ThisLoc  = ' -> at ObsPack_Init (in module ObsPack/obspack_mod.F90' 

    ! Assume that there are ObsPack data for today
    State_Diag%Do_ObsPack = .TRUE.

    WRITE (6,'("obspack_init")') 
    FLUSH(6)
    IF ( prtDebug ) THEN
       CALL DEBUG_MSG( '### OBSPACK_INIT: starting' )
    ENDIF

    !-----------------------------------------------------------------------
    ! If obs array is already allocated, assume that we need to
    ! write out existing results before reading new input.  This
    ! could happen in a multi-day run with daily input files.
    !-----------------------------------------------------------------------
    IF ( ASSOCIATED( State_Diag%ObsPack_Id ) ) THEN

       ! write output
       CALL ObsPack_Write_Output( am_I_Root, Input_Opt, State_Diag, RC )
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Error encountered in "ObsPack_Write_Output"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

       ! Deallocate arrays
       CALL ObsPack_Cleanup( am_I_Root, State_Diag, RC )
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Error encountered in "ObsPack_Cleanup"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
    ENDIF

    ! Reset the number of observations
    State_Diag%ObsPack_nObs = 0

    !-----------------------------------------------------------------------
    ! Form input and output file names with daily timestamps.
    ! Example: flask_input.2011063000_2011070100.nc and
    ! flask_output.2011063000_2011070100.nc. The current
    ! implementation uses today's time for the first timestamp, and
    ! today's time plus 24 hours for the second.  We could also use
    ! the IVALb and IVALe times, if the input files have been
    ! processed accordingly
    !-----------------------------------------------------------------------
    State_Diag%ObsPack_InFile  = TRIM( Input_Opt%ObsPackInputFile  )
    State_Diag%ObsPack_OutFile = TRIM( Input_Opt%ObsPackOutputFile )

    ! Replace YYYYMMDD with date and time
    CALL EXPAND_DATE( State_Diag%ObsPack_InFile,  yyyymmdd, hhmmss )
    CALL EXPAND_DATE( State_Diag%ObsPack_OutFile, yyyymmdd, hhmmss )

    WRITE(6,*) 'OBSPACK INPUT FILE  : ', TRIM( State_Diag%ObsPack_InFile  )
    WRITE(6,*) 'OBSPACK OUTPUT FILE : ', TRIM( State_Diag%ObsPack_OutFile )

    ! If we can't find a ObsPack file for today's date, return
    IF ( .NOT. FILE_EXISTS( TRIM( State_Diag%ObsPack_InFile ) ) ) THEN 
       State_Diag%Do_ObsPack = .FALSE.
       WRITE(6,*) 'INPUT FILE DOES NOT EXIST!'
       RETURN
    ENDIF

    ! Get number of tracers to sample
    State_Diag%ObsPack_nTracers = Input_Opt%N_ADVECT

    !-----------------------------------------------------------------------
    ! Get the list of lon/lat/alt at which to save out GEOS-Chem data
    !-----------------------------------------------------------------------
    CALL ObsPack_Read_Input( am_I_Root, Input_Opt, State_Diag, RC )
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Error encountered in "ObsPack_Write_Output"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    WRITE(6,*) 'OBSPACK NTRACERS: ', State_Diag%ObsPack_nTracers
    WRITE(6,*) 'OBSPACK NOBS: ', State_Diag%ObsPack_nObs

  END SUBROUTINE OBSPACK_INIT
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: obspack_read_input
!
! !DESCRIPTION: Subroutine OBSPACK\_READ\_INPUT allocates space
!  for variables in the input file, and reads those data into the
!  module arrays.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE ObsPack_Read_Input( am_I_root, Input_Opt, State_Diag, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE File_Mod,               ONLY : File_Exists
    USE Input_Opt_Mod,          ONLY : OptInput 
    USE State_Diag_Mod,         ONLY : DgnState
    USE m_netcdf_io_open,       ONLY : Ncop_Rd
    USE m_netcdf_io_get_dimlen, ONLY : Ncget_Dimlen
    USE m_netcdf_io_read
    USE m_netcdf_io_close,      ONLY : Nccl
!
! !INPUT PARAMETERS:
!
    LOGICAL,        INTENT(IN)    :: am_I_Root   ! Is this the root core
    TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
!
! !OUTPUT PARAMETERS:
!
    TYPE(DgnState), INTENT(INOUT) :: State_Diag  ! Diagnostics State object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
! 
! !REMARKS:
!  We read in all data available in the input file,
!  but it is possible that there are observations that
!  fall outside the IVALb-IVALe time period
!
! !REVISION HISTORY: 
!  05 Jun 2015 - A. Jacobson - first version
!  06 Dec 2018 - R. Yantosca - Implemented into the standard GEOS_Chem code
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER              :: fId, N, nObs, nTracers

    ! Arrays
    INTEGER, ALLOCATABLE :: central_time(:,:)
    INTEGER              :: st1d(1), ct1d(1)
    INTEGER              :: st2d(2), ct2d(2)

    ! Strings
    CHARACTER(LEN=30 )   :: varName
    CHARACTER(LEN=512)   :: ErrMsg
    CHARACTER(LEN=255)   :: ThisLoc

    !=======================================================================
    ! OBSPACK_READ_INPUT begins here
    !=======================================================================

    ! Initialize
    RC      = GC_SUCCESS
    ErrMsg  = ''
    ThisLoc = ' -> at ObsPack_Read_Input (in module ObsPack/obspack_mod.F90)'

    !=======================================================================
    ! If obs array is already allocated, assume that we need to
    ! write out existing results before reading new input.  This
    ! could happen in a multi-day run with daily input files.
    !=======================================================================
    IF ( ASSOCIATED( State_Diag%ObsPack_Id ) ) THEN

       ! write output
       CALL ObsPack_Write_Output( am_I_Root, Input_Opt, State_Diag, RC )
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Error encountered in "ObsPack_Write_Output"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

       ! Deallocate arrays
       CALL ObsPack_Cleanup( am_I_Root, State_Diag, RC )
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Error encountered in "ObsPack_Cleanup"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
    ENDIF

    !=======================================================================
    ! Exit if we are not on the root core
    !
    ! Or for MPI (e.g. WRF-GC), gather inforamation here from all cores
    !=======================================================================
    IF ( .not. am_I_Root ) RETURN

    !=======================================================================
    ! Get the number of observations in the input netCDF file
    ! (or exit if the input netCDF file does not exist)
    !=======================================================================

    ! First test if the input file exists
    IF ( .NOT. File_Exists( State_Diag%ObsPack_InFile ) ) THEN
       State_Diag%Do_ObsPack   = .FALSE.
       State_Diag%ObsPack_nObs = 0
       WRITE(6,*) 'INPUT FILE DOES NOT EXIST', State_Diag%ObsPack_InFile
       RETURN
    ENDIF

    ! If the file exists, open it.  Save the file ID to State_Diag.
    CALL Ncop_Rd( fId, State_Diag%ObsPack_InFile )
    State_Diag%ObsPack_fId = fId

    ! Get the # of observations in the file (and save it to State_Diag).
    CALL Ncget_Dimlen( fId, 'obs', nObs )
    State_Diag%ObsPack_nObs = nObs
    WRITE (6,'("[obspack_read_input] ",i10," obs in input file.")') nObs

    !=======================================================================
    ! Allocate the relevant fields of State_Diag
    !=======================================================================
    ALLOCATE( State_Diag%ObsPack_ID( nObs ), STAT=RC ) 
    CALL GC_CheckVar( 'State_Diag%ObsPack_Id', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
   
    ALLOCATE( State_Diag%ObsPack_nSamples( nObs ), STAT=RC )
    CALL GC_CheckVar( 'State_Diag%ObsPack_nSamples', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Diag%ObsPack_nSamples = 0

    ALLOCATE( State_Diag%ObsPack_Strategy( nObs ), STAT=RC )
    CALL GC_CheckVar( 'State_Diag%ObsPack_Strategy', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    
    ALLOCATE( State_Diag%ObsPack_Latitude( nObs ), STAT=RC ) 
    CALL GC_CheckVar( 'State_Diag%ObsPack_Latitude', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Diag%ObsPack_Latitude = 0.0_f4

    ALLOCATE( State_Diag%ObsPack_Longitude( nObs ), STAT=RC )
    CALL GC_CheckVar( 'State_Diag%ObsPack_Longitude', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Diag%ObsPack_Longitude = 0.0_f4

    ALLOCATE( State_Diag%ObsPack_Altitude( nObs ), STAT=RC )
    CALL GC_CheckVar( 'State_Diag%ObsPack_Altitude', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Diag%ObsPack_Altitude = 0.0_f4

    ALLOCATE( State_Diag%ObsPack_Ival_Start( nObs ), STAT=RC )
    CALL GC_CheckVar( 'State_Diag%ObsPack_Ival_Start', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Diag%ObsPack_Ival_Start = 0.0_f8

    ALLOCATE( State_Diag%ObsPack_Ival_Center( nObs ), STAT=RC )
    CALL GC_CheckVar( 'State_Diag%ObsPack_Ival_Center', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Diag%ObsPack_Ival_Center = 0.0_f8

    ALLOCATE( State_Diag%ObsPack_Ival_End( nObs ), STAT=RC )
    CALL GC_CheckVar( 'State_Diag%ObsPack_Ival_End', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Diag%ObsPack_Ival_End = 0.0_f8

    ALLOCATE( State_Diag%ObsPack_U( nObs ), STAT=RC )
    CALL GC_CheckVar( 'State_Diag%ObsPack_U', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Diag%ObsPack_U = 0.0_f4

    ALLOCATE( State_Diag%ObsPack_V( nObs ), STAT=RC )
    CALL GC_CheckVar( 'State_Diag%ObsPack_V', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Diag%ObsPack_V = 0.0_f4

    ALLOCATE( State_Diag%ObsPack_BLH( nObs ), STAT=RC )
    CALL GC_CheckVar( 'State_Diag%ObsPack_BLH', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Diag%ObsPack_BLH = 0.0_f4

    ALLOCATE( State_Diag%ObsPack_Q( nObs ), STAT=RC )
    CALL GC_CheckVar( 'State_Diag%ObsPack_Q', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Diag%ObsPack_Q = 0.0_f4

    ALLOCATE( State_Diag%ObsPack_Pressure( nObs ), STAT=RC )
    CALL GC_CheckVar( 'State_Diag%ObsPack_Pressure', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Diag%ObsPack_Pressure = 0.0_f4

    ALLOCATE( State_Diag%ObsPack_Temperature( nObs ), STAT=RC )
    CALL GC_CheckVar( 'State_Diag%ObsPack_Temperature', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Diag%ObsPack_temperature = 0.0_f4

    nTracers = State_Diag%ObsPack_nTracers
    ALLOCATE( State_Diag%ObsPack_Tracers( nObs, nTracers ), STAT=RC )
    CALL GC_CheckVar( 'State_Diag%ObsPack_nTracers', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Diag%ObsPack_tracers = 0.0_f4

    ! central_time is a local work array
    ALLOCATE( central_time( 6, nObs ), STAT=RC )
    CALL GC_CheckVar( 'State_Diag%ObsPack_Id', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    !=======================================================================
    ! Read input information from the netCDF file:
    !=======================================================================

    !----------------------------
    ! Read coordinate arrays
    !----------------------------
    st1d = (/ 1    /)
    ct1d = (/ nobs /)

    varName = 'latitude'
    CALL NcRd( State_Diag%ObsPack_Latitude,  fId, TRIM(varName), st1d, ct1d )
    
    varName = 'longitude'
    CALL NcRd( State_Diag%ObsPack_Longitude, fId, TRIM(varName), st1d, ct1d )

    varName = 'altitude'
    CALL NcRd( State_Diag%ObsPack_Altitude,  fId, TRIM(varName), st1d, ct1d )

    varName = 'sampling_strategy'
    CALL NcRd( State_Diag%ObsPack_Strategy,  fId, TRIM(varName), st1d, ct1d )  

    !----------------------------
    ! Read ID string
    !----------------------------
    st2d = (/ 1, 1      /)
    ct2d = (/ 100, nObs /)

    varName = 'ObsPack_id'
    CALL NcRd( State_Diag%ObsPack_Id, fId, TRIM(varName), st2d, ct2d )

    !----------------------------
    ! Read time components
    ! (Y, M, D, hr, min, sec)
    !----------------------------
    st2d = (/ 1, 1 /)
    ct2d = (/ 6, nobs /)

    varName = 'time_components'
    CALL NcRd( Central_Time, fId, TRIM(varName), st2d, ct2d )

    ! Close input netCDF file
    CALL NcCl( fId )
    State_Diag%ObsPack_fId = fId

    !=======================================================================
    ! Fill sampling window start and end times
    !=======================================================================
    DO N = 1, nObs

       !-----------------------
       ! Test for valid data
       !-----------------------
       IF ( State_Diag%ObsPack_Longitude(N) <   - 180.0_f4 .or.              &
            State_Diag%ObsPack_Longitude(N) >     180.0_f4 .or.              &
            State_Diag%ObsPack_Latitude(N)  <     -90.0_f4 .or.              &
            State_Diag%ObsPack_Latitude(N)  >      90.0_f4 .or.              &
            State_Diag%ObsPack_Altitude(N)  <   -1000.0_f4 .or.              &
            State_Diag%ObsPack_Altitude(N)  >  200000.0_f4 ) THEN

          ! Write an error message if the data is bad
          ErrMsg = '"WARNING!  Bad coordinates on obs with ObsPack_id: '
          WRITE( 6, '(2a)' ) TRIM( ErrMsg                      ),            &
                             TRIM( State_Diag%ObsPack_id(N) )

          WRITE (6,'("  -  lon ",f12.4,", lat ",f12.4,", alt ",f12.4,".")')              &
             State_Diag%ObsPack_Longitude(N),                                &
             State_Diag%ObsPack_Latitude(N),                                 &
             State_Diag%ObsPack_Altitude(N)

          ! Mark that the data will be skipped
          State_Diag%ObsPack_strategy(N) = 0
       ENDIF
       
       ! Compute the center time of the observation as seconds since 1970
       State_Diag%ObsPack_Ival_Center(N) =                                   &
           Seconds_Since_1970( central_time(1,N),                            &
                               central_time(2,N),                            &
                               central_time(3,N),                            &
                               central_time(4,N),                            &
                               central_time(5,N),                            & 
                               central_time(6,N)  )

       ! Pick the start and end time of the averaging interval
       ! depending on the averaging strategy listed in the file
       SELECT CASE ( State_Diag%ObsPack_Strategy(N) )

          !------------------
          ! DO-NOT-SAMPLE
          !------------------
          CASE( 0 ) 
             State_Diag%ObsPack_Ival_Start(N) =                              &
                State_Diag%ObsPack_Ival_Center(N) + 1

             State_Diag%ObsPack_Ival_End(N) =                                &
                State_Diag%ObsPack_Ival_Center(N) - 1

          !------------------
          ! 4-hour window
          !------------------
          CASE( 1 )
             State_Diag%ObsPack_Ival_Start(N) =                              &
                State_Diag%ObsPack_Ival_Center(N) - 7200.0_f8

             State_Diag%ObsPack_Ival_End(N) =                                &
                State_Diag%ObsPack_Ival_center(N) + 7200.0_f8

          !------------------
          ! 1-hour window
          !------------------
          CASE( 2 ) 
             State_Diag%ObsPack_Ival_Start(N) =                              &
                State_Diag%ObsPack_Ival_Center(N) - 1800.0_f8

             State_Diag%ObsPack_Ival_end(N) =                                &
                State_Diag%ObsPack_Ival_center(N) + 1800.0_f8

          !------------------
          ! 1-hour window
          !------------------
          CASE( 3 )
             State_Diag%ObsPack_Ival_Start(N) =                              &
                State_Diag%ObsPack_Ival_Center(N) - 2700.0_f8

             State_Diag%ObsPack_Ival_End(N) =                                &
                State_Diag%ObsPack_Ival_Center(N) + 2700.0_f8

          !------------------
          ! Exit w/ error
          !------------------
          CASE DEFAULT
             ErrMsg = 'Observation with ObsPack_id: '                     // &
                      TRIM( ADJUSTL( State_Diag%ObsPack_id(N) ) )         // &
                      ' has an unknown or invalid sampling strategy!'
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN

       END SELECT

    ENDDO

    WRITE(6,*) 'ObsPack No. OBS',nobs

    ! Deallocate the central time array
    IF ( ALLOCATED( central_time ) ) THEN
       DEALLOCATE( central_time, STAT=RC )
       CALL GC_CheckVar( 'Central_Time', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

  END SUBROUTINE ObsPack_READ_INPUT
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ObsPack_cleanup
!
! !DESCRIPTION: Subroutine ObsPack\_CLEANUP deallocates all ObsPack
!  fields of State_Diag.
!\\
!\\
! !INTERFACE:
!
 SUBROUTINE ObsPack_Cleanup( am_I_Root, State_Diag, RC )
!
! !USES:
!     
   USE ErrCode_Mod
   USE State_Diag_Mod, ONLY : DgnState
!
! !INPUT PARAMETERS:
!
    LOGICAL,        INTENT(IN)    :: am_I_Root   ! Are we on the root core?
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(DgnState), INTENT(INOUT) :: State_Diag  ! Diagnostics State object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,         INTENT(OUT)  :: RC          ! Success or failure?
! 
! !REVISION HISTORY: 
!  05 Jun 2015 - A. Jacobson - first version
!  05 Dec 2018 - R. Yantosca - Implemented into GEOS-Chem standard code
!EOP
!------------------------------------------------------------------------------
!BOC

    ! Assume success
    RC = GC_SUCCESS

    !=======================================================================
    ! Deallocate ObsPack variables
    !=======================================================================
    IF ( ASSOCIATED( State_Diag%ObsPack_Id ) ) THEN
       DEALLOCATE( State_Diag%ObsPack_Id, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%ObsPack_Id', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%ObsPack_Id => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%ObsPack_nSamples ) ) THEN
       DEALLOCATE( State_Diag%ObsPack_nSamples, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%ObsPack_nSamples', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%ObsPack_nSamples => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%ObsPack_Strategy ) ) THEN
       DEALLOCATE( State_Diag%ObsPack_Strategy, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%ObsPack_Strategy', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%ObsPack_Strategy => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%ObsPack_Latitude ) ) THEN
       DEALLOCATE( State_Diag%ObsPack_Latitude, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%ObsPack_Latitude', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%ObsPack_Latitude => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%ObsPack_Longitude ) ) THEN
       DEALLOCATE( State_Diag%ObsPack_Longitude, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%ObsPack_Longitude', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%ObsPack_Longitude => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%ObsPack_Altitude ) ) THEN
       DEALLOCATE( State_Diag%ObsPack_Altitude, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%ObsPack_Altitude', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%ObsPack_Altitude => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%ObsPack_Ival_Start ) ) THEN
       DEALLOCATE( State_Diag%ObsPack_Ival_Start, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%ObsPack_Ival_Start', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%ObsPack_Ival_Start => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%ObsPack_Ival_Center ) ) THEN
       DEALLOCATE( State_Diag%ObsPack_Ival_Center, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%ObsPack_Ival_Center', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%ObsPack_Ival_Center => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%ObsPack_Ival_End ) ) THEN
       DEALLOCATE( State_Diag%ObsPack_Ival_End, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%ObsPack_Ival_End', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%ObsPack_Ival_End => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%ObsPack_Accum_Weight ) ) THEN
       DEALLOCATE( State_Diag%ObsPack_Accum_Weight, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%ObsPack_Accum_Weight', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%ObsPack_Accum_Weight => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%ObsPack_U ) ) THEN
       DEALLOCATE( State_Diag%ObsPack_U, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%ObsPack_U', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%ObsPack_U => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%ObsPack_V ) ) THEN
       DEALLOCATE( State_Diag%ObsPack_V, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%ObsPack_V', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%ObsPack_V => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%ObsPack_BLH ) ) THEN
       DEALLOCATE( State_Diag%ObsPack_BLH, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%ObsPack_BLH', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%ObsPack_BLH => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%ObsPack_Q ) ) THEN
       DEALLOCATE( State_Diag%ObsPack_Q, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%ObsPack_Q', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%ObsPack_Q => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%ObsPack_Pressure ) ) THEN
       DEALLOCATE( State_Diag%ObsPack_Pressure, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%ObsPack_Pressure', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%ObsPack_Pressure => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%ObsPack_Temperature ) ) THEN
       DEALLOCATE( State_Diag%ObsPack_Temperature, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%ObsPack_Temperature', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%ObsPack_Temperature => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%ObsPack_Tracers ) ) THEN
       DEALLOCATE( State_Diag%ObsPack_Tracers, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%ObsPack_Tracers', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%ObsPack_Tracers => NULL()
    ENDIF

  END SUBROUTINE ObsPack_Cleanup
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: obspack_write_output
!
! !DESCRIPTION: Subroutine ObsPack\_WRITE\_OUTPUT computes window averages
!  and writes data to output file
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE ObsPack_Write_Output( am_I_Root, Input_Opt, State_Diag, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE Input_Opt_Mod,  ONLY : OptInput
    USE State_Diag_Mod, ONLY : DgnState
    USE Time_Mod,       ONLY : Get_NHMSb
    USE Time_Mod,       ONLY : Get_NHMSe
    USE Time_Mod,       ONLY : Get_NYMDb
    USE Time_Mod,       ONLY : Get_NYMDe
    USE Time_Mod,       ONLY : System_Timestamp
    USE Time_Mod,       ONLY : Ymd_Extract
    
    ! NetCDF modules
    USE m_netcdf_io_define
    USE m_netcdf_io_create
    USE m_netcdf_io_write
    USE m_netcdf_io_close     
!
! !INPUT PARAMETERS:
!
    LOGICAL,        INTENT(IN)    :: am_I_Root   ! Is this the root CPU?
    TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(DgnState), INTENT(INOUT) :: State_Diag  ! Diagnostic State object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
! 
! !REVISION HISTORY: 
!  05 Jun 2015 - A. Jacobson - First version
!  06 Dec 2018 - R. Yantosca - Implemented into the standard GEOS-Chem code
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES
!
    ! Scalars
    INTEGER               :: fId,          nObs,          nTracers
    INTEGER               :: ymd,          hms,           omode
    INTEGER               :: yr,           mo,            da
    INTEGER               :: hr,           mn,            sc
    INTEGER               :: dimid_obs,    dimid_tracer
    INTEGER               :: dimid_tnmlen, dimid_char100  
    INTEGER               :: vid,          N

    ! Arrays
    INTEGER               :: st1d(1),      ct1d(1),       dims_1d(1)
    INTEGER               :: st2d(2),      ct2d(2),       dims_2d(2)
    REAL(f8), ALLOCATABLE :: avetime(:)

    ! Strings
    CHARACTER(LEN=16)     :: stamp
    CHARACTER(LEN=30)     :: varName
    CHARACTER(LEN=255)    :: attstring
    CHARACTER(LEN=255)    :: ThisLoc
    CHARACTER(LEN=512)    :: ErrMsg

    !=======================================================================
    ! ObsPack_Write_Output begins here!
    !=======================================================================

    ! Initialize
    RC       =  GC_SUCCESS
    ErrMsg   = ''
    ThisLoc  = ' -> at ObsPack_Write_Output (in module Obspack/obspack_mod.F90)'
    nObs     = State_Diag%ObsPack_nObs
!    nTracers = State_Diag%ObsPack_nTracers

    ! Exit if there are no observations
    IF ( nObs == 0 ) RETURN
 
    !=======================================================================
    ! Compute averages
    !=======================================================================
    DO N = 1, nObs

       ! Only compute averages if there are samples
       IF ( State_Diag%ObsPack_nSamples(N) > 0 ) THEN
          State_Diag%ObsPack_U(N)           =                               &
          State_Diag%ObsPack_U(N)           / State_Diag%ObsPack_nSamples(N)

          State_Diag%ObsPack_V(N)           =                               &
          State_Diag%ObsPack_V(N)           / State_Diag%ObsPack_nSamples(N)

          State_Diag%ObsPack_BLH(N)         =                               &
          State_Diag%ObsPack_BLH(N)         / State_Diag%ObsPack_nSamples(N)

          State_Diag%ObsPack_Q(N)           =                               &
          State_Diag%ObsPack_Q(N)           / State_Diag%ObsPack_nSamples(N)

          State_Diag%ObsPack_Temperature(N) =                               &
          State_Diag%ObsPack_temperature(N) / State_Diag%ObsPack_nSamples(N)

          State_Diag%ObsPack_Pressure(N)    =                               &
          State_Diag%ObsPack_Pressure(N)    / State_Diag%ObsPack_nSamples(N)

          State_Diag%ObsPack_Tracers(N,:)   =                               &
          State_Diag%ObsPack_tracers(N,:)   / State_Diag%ObsPack_nSamples(N)
       ENDIF
    ENDDO

    !=======================================================================
    ! Open file for output
    !=======================================================================
    WRITE (6,'("[obspack_write_output] Creating ",a,"...")') State_Diag%ObsPack_outfile
   
    ! Create netCDF file and save the file ID in State_Diag
    CALL NcCr_Wr( fId, State_Diag%ObsPack_outfile )
    State_Diag%ObsPack_fId = fId

    ! Turn filling off
    CALL NcSetFill( fId, NF_NOFILL, omode )

    ! Define dimensions
    varName = 'obs'
    CALL NcDef_Dimension( fId, TRIM(varName), NF_UNLIMITED,   dimid_obs     )

    varName = 'tracer'
    CALL NcDef_Dimension( fId, TRIM(varName), nTracers,       dimid_tracer  )

    varName = 'tracer_name_len'
! What is this???
!    CALL NcDef_Dimension( fId, TRIM(varName), tracer_namelen, dimid_tnmlen  )
    CALL NcDef_Dimension( fId, TRIM(varName), 31,             dimid_tnmlen  )

    varName = 'char100'
    CALL NcDef_Dimension( fId, TRIM(varName), 100,            dimid_char100 )

    !=======================================================================
    ! Set global attributes
    !=======================================================================
    stamp = SYSTEM_TIMESTAMP()
    attstring = 'GEOS-Chem simulation at ' // stamp
    CALL NcDef_Glob_Attributes( fId, 'History',           TRIM(attstring) )
    CALL NcDef_Glob_Attributes( fId, 'Conventions',       'CF-1.4'  )

    ymd = GET_NYMDb()
    hms = GET_NHMSb()
    CALL Ymd_Extract( ymd, yr, mo, da)
    CALL Ymd_Extract( hms, hr, mn, sc)
    WRITE( attstring, 100 ) yr, mo, da, hr, mn, sc
100 FORMAT( i4.4,"/",i2.2,"/",i2.2," ",i2.2,":",i2.2,":",i2.2, " UTC" )
    CALL NcDef_Glob_Attributes( fId, 'model_start_date',  TRIM(attstring) )

    ymd = GET_NYMDe()
    hms = GET_NHMSe()
    CALL YMD_EXTRACT( ymd, yr, mo, da )
    CALL YMD_EXTRACT( hms, hr, mn, sc )
    WRITE( attstring, 100 ) yr,mo,da,hr,mn,sc
    CALL NcDef_Glob_Attributes( fId, 'model_end_date',    TRIM(attstring)   )

    !=======================================================================
    ! Define variables and attributes
    !=======================================================================
    dims_2d = (/ dimid_char100, dimid_obs /)
    vid     = 0
    CALL NcDef_Variable      ( fId, 'obspack_id', NF_CHAR,  2, dims_2d, vid  )
    CALL NcDef_Var_Attributes( fId, vid, 'long_name',      'obspack_id'      )
    CALL NcDef_Var_Attributes( fId, vid, 'units',          'unitless'                ) 

    vid  = vid + 1
    dims_2d = (/ dimid_tracer, dimid_obs /)
    CALL NcDef_Variable( fId, 'flask', NF_DOUBLE, 2, dims_2d, vid )
    CALL NcDef_Var_Attributes( fId, vid, 'long_name', 'mole_fraction_of_trace_gas_in_air' )
    CALL NcDef_Var_Attributes( fId, vid, 'units',     'mol mol-1'      )
    CALL NcDef_Var_Attributes( fId, vid, '_FillValue', -1e34          )

    vid  = vid + 1
    dims_1d = (/ dimid_obs /)
    CALL NcDef_Variable( fId, 'nsamples', NF_INT, 1, dims_1d, vid )
    CALL NcDef_Var_Attributes( fId, vid, 'long_name', 'no. of model samples')
    CALL NcDef_Var_Attributes( fId, vid, 'units',     'unitless' )
    CALL NcDef_Var_Attributes( fId, vid, 'comment',     'Number of discrete model samples in average.' )

    vid  = vid + 1
    dims_1d = (/ dimid_obs /)
    CALL NcDef_Variable( fId, 'averaging_time', NF_INT, 1, dims_1d, vid )
    CALL NcDef_Var_Attributes( fId, vid, 'long_name', 'averaging time')
    CALL NcDef_Var_Attributes( fId, vid, 'units',     'seconds' )
    CALL NcDef_Var_Attributes( fId, vid, 'comment',     'Amount of model time over which this sample is averaged.' )

    vid  = vid + 1
    dims_1d = (/ dimid_obs /)
    CALL NcDef_Variable( fId, 'u', NF_DOUBLE, 1, dims_1d, vid )
    CALL NcDef_Var_Attributes( fId, vid, 'long_name', 'u-wind')
    CALL NcDef_Var_Attributes( fId, vid, 'units',     'm s^-1' )

    vid  = vid + 1
    dims_1d = (/ dimid_obs /)
    CALL NcDef_Variable( fId, 'v', NF_DOUBLE, 1, dims_1d, vid )
    CALL NcDef_Var_Attributes( fId, vid, 'long_name', 'v-wind')
    CALL NcDef_Var_Attributes( fId, vid, 'units',     'm s^-1' )

    vid  = vid + 1
    dims_1d = (/ dimid_obs /)
    CALL NcDef_Variable( fId, 'blh', NF_DOUBLE, 1, dims_1d, vid )
    CALL NcDef_Var_Attributes( fId, vid, 'long_name', 'v-wind')
    CALL NcDef_Var_Attributes( fId, vid, 'units',     'm s^-1' )

    vid  = vid + 1
    dims_1d = (/ dimid_obs /)
    CALL NcDef_Variable( fId, 'q', NF_DOUBLE, 1, dims_1d, vid )
    CALL NcDef_Var_Attributes( fId, vid, 'long_name', 'mass_fraction_of_water_inair')
    CALL NcDef_Var_Attributes( fId, vid, 'units',     'kg water (kg air)^-1' )

    vid  = vid + 1
    dims_1d = (/ dimid_obs /)
    CALL NcDef_Variable( fId, 'pressure', NF_DOUBLE, 1, dims_1d, vid )
    CALL NcDef_Var_Attributes( fId, vid, 'long_name', 'pressure')
    CALL NcDef_Var_Attributes( fId, vid, 'units',     'Pa' )

    vid  = vid + 1
    dims_1d = (/ dimid_obs /)
    CALL NcDef_Variable( fId, 'temperature', NF_DOUBLE, 1, dims_1d, vid )
    CALL NcDef_Var_Attributes( fId, vid, 'long_name', 'temperature')
    CALL NcDef_Var_Attributes( fId, vid, 'units',     'K' )

    ! End the definition section
    CALL NcEnd_def( fId )
 
    !=======================================================================
    ! Write variables to disk
    !=======================================================================
    varName = 'obspack_id'
    st2d    = (/ 1,   1    /)
    ct2d    = (/ 100, nObs /)
    CALL NcWr( State_Diag%ObsPack_Id, fId, TRIM(varName), st2d, ct2d )

    varName = 'flask'
    st2d    = (/ 1,        1    /)
    ct2d    = (/ ntracers, nObs /)
    CALL NcWr( State_Diag%ObsPack_Tracers, fId, TRIM(varName), st2d, ct2d )

    varName = 'nsamples'
    st1d    = (/ 1    /)
    ct1d    = (/ nObs /)
    CALL NcWr( State_Diag%ObsPack_nSamples, fId, TRIM(varName), st1d, ct1d )

    varname = 'averaging_time'
    st1d    = (/ 1    /)
    ct1d    = (/ nObs /)
    ALLOCATE( avetime( nObs ))
    aveTime = State_Diag%ObsPack_Ival_End - State_Diag%ObsPack_Ival_Start
    CALL NcWr( aveTime, fId, TRIM(varName),  st1d, ct1d )
    DEALLOCATE( aveTime )

    varName = 'u'
    st1d    = (/ 1    /)
    ct1d    = (/ nObs /)
    CALL NcWr( State_Diag%ObsPack_V, fId, TRIM(varName), st1d, ct1d )

    varName = 'v'
    st1d    = (/ 1    /)
    ct1d    = (/ nObs /)
    CALL NcWr( State_Diag%ObsPack_V, fId, TRIM(varName), st1d, ct1d )
    
    varName = 'blh'
    st1d    = (/ 1    /)
    ct1d    = (/ nObs /)
    CALL NcWr( State_Diag%ObsPack_BLH, fId, TRIM(varName), st1d, ct1d )

    varName = 'q'
    st1d    = (/ 1    /)
    ct1d    = (/ nObs /)
    CALL NcWr( State_Diag%ObsPack_Q, fId, TRIM(varName), st1d, ct1d )

    varName = 'pressure'
    st1d    = (/ 1    /)
    ct1d    = (/ nObs /)
    CALL NcWr( State_Diag%ObsPack_Pressure, fId, TRIM(varName), st1d, ct1d )

    varName = 'temperature'
    st1d = (/ 1      /)
    ct1d = (/ nObs   /)
    CALL NcWr( State_Diag%ObsPack_Temperature, fId, TRIM(varName), st1d, ct1d )

    ! Close the netCDF file
    CALL NcCl( fId )

    ! Cleanup arrays
    CALL ObsPack_Cleanup( am_I_Root, State_Diag, RC )
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Error encountered in "ObsPack_Cleanup!"'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

  END SUBROUTINE ObsPack_Write_Output
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ObsPack\_SAMPLE
!
! !DESCRIPTION: Subroutine ObsPack\_SAMPLE performs the model sampling
!  and saves concentrations to locations corresponding to a flight track.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE ObsPack_Sample( am_I_Root, yyyymmdd,  hhmmss,     Input_Opt,    &
                             State_Met, State_Chm, State_Diag, RC           )
!
! !USES:
!
    USE ErrCode_Mod,
    USE Input_Opt_Mod,  ONLY : OptInput
    USE State_Chm_Mod,  ONLY : ChmState
    USE State_Met_Mod,  ONLY : MetState
    USE State_Diag_Mod, ONLY : DgnState
    USE Time_Mod,       ONLY : Ymd_Extract
    USE UnitConv_Mod,   ONLY : Convert_Spc_Units
!
! !INPUT PARAMETERS: 
!
    LOGICAL,        INTENT(IN)    :: am_I_Root   ! Are we on the root CPU?
    INTEGER,        INTENT(IN)    :: yyyymmdd    ! Current date
    INTEGER,        INTENT(IN)    :: hhmmss      ! Current time
    TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
    TYPE(MetState), INTENT(IN)    :: State_Met   ! Meteorology State object
!
! !INPUT/OUTPUT PARAMETERS: 
!
    TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
    TYPE(DgnState), INTENT(INOUT) :: State_Diag  ! Chemistry State object

!
! !OUTPUT PARAMETERS: 
!
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
!
! !REMARKS:
!
! !REVISION HISTORY: 
!  08 Jun 2015 - A. Jacobson, A. Schuh - imported from Andrew Schuh's
!                                        ct_mod.F, itself modified from
!                                        planeflight_mod.F
!  03 Mar 2017 - A. Jacobson - Update to v11 (get species in "v/v dry")
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    LOGICAL             :: prtDebug
    INTEGER             :: I,  J,  L,  T,  N
    INTEGER             :: Yr, Mo, Da, Hr, Mn, Sc
    REAL(f8)            :: Elapsed, LastFileWrite, this_taue, this_taus

    ! Strings
    CHARACTER(LEN=255)  :: PriorUnit, ErrMsg, ThisLoc

    !=================================================================
    ! ObsPack_Sample begins here
    !=================================================================

    ! Initialize
    RC       =  GC_SUCCESS
    ErrMsg   = ''
    ThisLoc  = ' -> at ObsPack_Sample (in module ObsPack/obspack_mod.F90)'
    prtDebug = ( am_I_Root .and. Input_Opt%LPRT)

    ! Return if ObsPack sampling is turned off (perhaps
    ! because there are no data at this time).
    IF ( .not. State_Diag%DO_ObsPack ) RETURN

    IF ( prtDebug ) THEN
       CALL DEBUG_MSG( '### ObsPack_SAMPLE: starting' )
    ENDIF

    ! Ensure that units of species are "v/v dry", which is dry
    ! air mole fraction.  Capture the InUnit value, this is
    ! what the units are prior to this call.  After we sample
    ! the species, we'll call this again requesting that the
    ! species are converted back to the InUnit values.
    CALL Convert_Spc_Units( am_I_root, Input_Opt, State_Met,                 &
                            State_Chm, "v/v dry", RC,       PriorUnit       )

    ! Trap potential errors
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Error encountered in "Convert_Units"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    !=======================================================================
    !=======================================================================

    ! Extract date and time into components
    CALL Ymd_Extract( yyyymmdd, Yr, Mo, Da )
    CALL Ymd_Extract( hhmmss,   Hr, Mn, Sc )
 
    ! Compute elapsed seconds since 1970
    Elapsed = Seconds_Since_1970( Yr, Mo, Da, Hr, Mn, Sc )
    This_taue = elapsed

    ! Don't use TAU values
    THIS_TAUS = THIS_TAUE !!( GET_TS_DIAG() / 60d0 )
  
    !=======================================================================

    DO N = 1, State_Diag%ObsPack_nObs

       IF ( State_Diag%ObsPack_Strategy(N) == 0 ) CYCLE

       IF ( State_Diag%ObsPack_Ival_Start(N) <= THIS_TAUS .and.          &
            State_Diag%ObsPack_Ival_End(N)   >= THIS_TAUE ) THEN

           WRITE (6,'("   ")')
           FLUSH(6)
           WRITE (6,'("   ")')
           FLUSH(6)
           WRITE (6,'("sampling obs ",i10,", obspack_id: ",a)') &
                N, trim(State_Diag%ObsPack_id(N))
           FLUSH(6)

          ! Return grid box indices for the chemistry region
          CALL ObsPack_Get_Indices( N, State_Met, State_Diag, I, J, L, RC )

          ! Trap potential errors
          IF ( RC /= GC_SUCCESS ) THEN
             ErrMsg = 'Error encountered in "ObsPack_Get_Indices"!'
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN
          ENDIF

          ! NOTE: NEED A MAPPING TO ONLY SELECT A SUBSET OF TRACERS!
          DO T = 1, State_Diag%ObsPack_nTracers
             State_Diag%ObsPack_Tracers(N,T) =                               &
             State_Diag%ObsPack_Tracers(N,T) +  State_Chm%SPECIES(I,J,L,T)
          ENDDO

          State_Diag%ObsPack_U(N)           =                                &
          State_Diag%ObsPack_U(N)           + State_Met%U(I,J,L) 

          State_Diag%ObsPack_V(N)           =                                &
          State_Diag%ObsPack_V(N)           + State_Met%V(I,J,L) 

          State_Diag%ObsPack_BLH(N)         =                                &
          State_Diag%ObsPack_BLH(N)         + State_Met%PBLH(I,J) 

          State_Diag%ObsPack_Q(N)           =                                &
          State_Diag%ObsPack_Q(N)           + State_Met%SPHU(I,J,L) 

          State_Diag%ObsPack_pressure(N)    =                                &
          State_Diag%ObsPack_pressure(N)    + State_Met%PMID(I,J,L) 

          State_Diag%ObsPack_temperature(N) =                                &
          State_Diag%ObsPack_temperature(N) + State_Met%T(I,J,L) 

          State_Diag%ObsPack_nSamples(N)    =                                &
          State_Diag%ObsPack_nSamples(N)    + 1

       ENDIF
    ENDDO

    ! Return State_Chm%SPECIES to whatever units they had
    ! coming into this routine
    call Convert_Spc_Units( am_I_root, Input_Opt, State_Met,                 &
                            State_Chm, PriorUnit, RC                        )

    ! Trap potential errors
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Error encountered in "Convert_Units"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

  END SUBROUTINE ObsPack_Sample
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: 
!
! !DESCRIPTION: Subroutine ObsPack\_GET\_GRID\_INDICES returns the
!  grid box indices (I, J, L) corresponding to the input point
!  defined by longitude, latitude, altitude.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE ObsPack_Get_Indices( iObs, State_Met, State_Diag, I, J, L, RC )
!
! !USES:
!
    USE GC_GRID_MOD,    ONLY : GET_XOFFSET
    USE GC_GRID_MOD,    ONLY : GET_YOFFSET
    USE CMN_SIZE_MOD,   ONLY : LLPAR, DISIZE, DJSIZE, IIPAR
    USE ErrCode_Mod
    USE State_Met_Mod,  ONLY : MetState
    USE State_Diag_Mod, ONLY : DgnState 
!
! !INPUT PARAMETERS: 
!
    INTEGER,        INTENT(IN)  :: iObs         ! ObsPack Observation number
    TYPE(MetState), INTENT(IN)  :: State_Met    ! Meteorology State object
    TYPE(DgnState), INTENT(IN)  :: State_Diag   ! Diagnostics State object
!
! !OUTPUT PARAMETERS: 
!
    INTEGER,        INTENT(OUT) :: I, J, L      ! Lon, lat, level indices
    INTEGER,        INTENT(OUT) :: RC           ! Success or failure?
!
! !REMARKS:
!
! !REVISION HISTORY:
!  5 Jun 2015 - A. Jacobson - First version
!  3 Mar 2017 - A. Jacobson - Update to v11 (use State_Met%BXHEIGHT instead of my own hypsometry)
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER             :: idx, I0, J0
    REAL(f8)            :: Z

    ! 
    CHARACTER(LEN=255)  :: ErrMsg, ThisLoc

    !========================================================================
    ! ObsPack_Get_Indices begins here!
    !========================================================================

    ! Initialize
    RC      = GC_SUCCESS
    ErrMsg  = ''
    ThisLoc = ' -> at ObsPack_Get_Indices (in module ObsPack/obspack_mod.F90'

    !-----------------------------------------------------------------------
    ! Exit with error if the observations have invalid lon/lat/alt values
    !-----------------------------------------------------------------------
    IF ( State_Diag%ObsPack_Longitude(iObs) <    -180.0_f8   .or.            & 
         State_Diag%ObsPack_Longitude(iObs) >     180.0_f8   .or.            &
         State_Diag%ObsPack_Latitude(iObs)  <     -90.0_f8   .or.            &
         State_Diag%ObsPack_Latitude(iObs)  >      90.0_f8   .or.            &
         State_Diag%ObsPack_Altitude(iObs)  <   -1000.0_f8   .or.            &
         State_Diag%ObsPack_Altitude(iObs)  >  200000.0_f8 ) THEN

       ! Print error message and exit
       I = -1
       J = -1
       L = -1
       CALL GC_Error( ErrMsg, Rc, ThisLoc )
       RETURN
    ENDIF

    ! Added correct definitions for I and J based on nested regions 
    ! (lds, 8/25/11)
    I0 = GET_XOFFSET( GLOBAL=.TRUE. )
    J0 = GET_YOFFSET( GLOBAL=.TRUE. )

    !-----------------------------------------------------------------------
    ! Find I corresponding to the ObsPack longitude value
    !-----------------------------------------------------------------------
    I = INT( ( State_Diag%ObsPack_Longitude(iObs) + 180.0_f8                 &
                                                  - ( I0 * DISIZE ) )        &
                                                  / DISIZE + 1.5d0    )

    ! Handle date line correctly (bmy, 4/23/04)
    IF ( I > IIPAR ) I = I - IIPAR

    !-----------------------------------------------------------------------
    ! Find J corresponding to the ObsPack latitude value
    !-----------------------------------------------------------------------
    J = INT( ( State_Diag%ObsPack_Latitude(iObs)  +  90.0_f8                 &
                                                  - ( J0 * DJSIZE ) )        &
                                                  / DJSIZE + 1.5d0    )

    !-----------------------------------------------------------------------
    ! Find L corresponding to the ObsPack altitude value
    !-----------------------------------------------------------------------
    Z = 0.0_f8
    DO L = 1, LLPAR
       Z = Z + State_Met%BXHEIGHT(I,J,L)
       IF ( Z >= State_Diag%ObsPack_Altitude(iObs) ) RETURN 
    ENDDO

    !========================================================================
    ! Issue error message if we get here.
    !========================================================================
    !WRITE (6,*) 'At I,J =',i,j
    !FLUSH(6)

    DO L = 1,LLPAR
       Z=Z+State_Met%BXHEIGHT(I,J,L)
       WRITE (6,*) L,' bxheight and z:', State_Met%BXHEIGHT(I,J,L),Z
    ENDDO
    WRITE( ErrMsg, 100 ) State_Diag%ObsPack_Altitude(iObs), Z,               &
                         State_Diag%ObsPack_Latitude(iObs),                  &
                         State_Diag%ObsPack_Longitude
100 FORMAT( "Altitude ",f11.4, "m exceeds TOA of ",f11.4,"m at ",f11.4,"N,",f11.4,"E" )

    CALL GC_Error( ErrMsg, RC, ThisLoc )
    
  END SUBROUTINE ObsPack_Get_Indices
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Seconds_Since_1970
!
! !DESCRIPTION: Computes the seconds that have elapsed since 1970/01/01
!  at 00:00 UTC.
!\\
!\\
! !INTERFACE:
!
  FUNCTION Seconds_Since_1970( Year, Month,  Day,                            &
                               Hour, Minute, Second ) RESULT( Elapsed_Sec )
!
! !USES:
!
    USE Julday_Mod, ONLY : Julday
!
! !INPUT PARAMETERS: 
!
    INTEGER, INTENT(IN) :: Year, Month, Day       ! Current date
    INTEGER, INTENT(IN) :: Hour, Minute, Second   ! Current time
!
! !RETURN VALUE: 
!
    REAL(f8)            :: Elapsed_Sec            ! Elapsed seconds since
                                                  ! 1970/01/01 00:00 UTC
!
! !REMARKS:
!
! !REVISION HISTORY:
!  06 Dec 2018 - R. Yantosca - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    REAL(f8) :: FracDay, Jsec

    !=======================================================================
    ! Seconds_Since_1970 begins here!
    !=======================================================================

    ! Compute the fractional day
    FracDay = DBLE( Day ) + ( DBLE( Hour   ) /    24.0_f8 )  +               & 
                            ( DBLE( Minute ) /  3600.0_f8 )  +               &
                            ( DBLE( Second ) / 86400.0_f8 ) 

    ! Compute the Astronomical Julian Date (in decimal days)
    ! and then convert the result to seconds.
    JSec = JulDay( Year, Month, FracDay ) * 86400.0_f8

    ! Compute the Elapsed seconds since 1970/01/01 by subtracting the
    ! Astronomical Julian Date at 1970/01/01 00:00 UTC converted to seconds.
    ! Just keep the integer part, since we are dealing in integral seconds.
    ! NINT ensures that we round up in case there is underflow.
    Elapsed_Sec = NINT( JSec - 210866760000.0_f8 )

  END FUNCTION Seconds_Since_1970
!EOC
END MODULE ObsPack_Mod


