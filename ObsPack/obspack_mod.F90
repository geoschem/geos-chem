!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: obspack_mod.F90
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
  PUBLIC  :: ObsPack_Init
  PUBLIC  :: ObsPack_Cleanup
  PUBLIC  :: ObsPack_Sample
  PUBLIC  :: ObsPack_Write_Output
  PUBLIC  :: ObsPack_SpeciesMap_Init
  PUBLIC  :: ObsPack_SpeciesMap_Cleanup
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
! !REVISION HISTORY:
!  04 Jun 2015 - A. Jacobson - Adapted from v10.1 planeflight_mod.f, following
!                              similar work done in v9.2 by Andrew Schuh.
!  05 Dec 2018 - R. Yantosca - Implemented into the standard GEOS-Chem code
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
  INTEGER, PARAMETER :: CHAR_LEN_OBS  = 200  ! Length of ObsPack ID strings
  INTEGER, PARAMETER :: CHAR_LEN_SPEC = 31   ! Length of species names

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
  SUBROUTINE ObsPack_Init( yyyymmdd, hhmmss, Input_Opt, State_Diag, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE Error_Mod,      ONLY : Debug_Msg
    USE File_Mod,       ONLY : File_Exists
    USE Input_Opt_Mod,  ONLY : OptInput 
    USE State_Diag_Mod, ONLY : DgnState
    USE Time_Mod,       ONLY : Expand_Date
!
! !INPUT PARAMETERS:
!
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
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    LOGICAL            :: prtDebug
    INTEGER            :: N

    ! Strings
    CHARACTER(LEN=255) :: ErrMsg, ThisLoc

    !=======================================================================
    ! OBSPACK_INIT begins here
    !=======================================================================

    ! Initialize
    RC       =  GC_SUCCESS
    prtDebug = ( Input_Opt%LPRT .and. Input_Opt%amIRoot )
    ErrMsg   = ''
    ThisLoc  = ' -> at ObsPack_Init (in module ObsPack/obspack_mod.F90' 

    ! Assume that there are ObsPack data for today
    State_Diag%Do_ObsPack = .TRUE.

    IF ( prtDebug ) THEN
       CALL DEBUG_MSG( '### OBSPACK_INIT: starting' )
    ENDIF

    !-----------------------------------------------------------------------
    ! If obs array is already allocated, assume that we need to
    ! write out existing results before reading new input.  This
    ! could happen in a multi-day run with daily input files.
    !-----------------------------------------------------------------------
    IF ( ASSOCIATED( State_Diag%ObsPack_Id ) ) THEN

       ! Write any remaining ObsPack data to disk, and immediately
       ! thereafter free the ObsPack pointer fields of State_Diag
       CALL ObsPack_Write_Output( Input_Opt, State_Diag, RC )
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Error encountered in "ObsPack_Write_Output"!'
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
    State_Diag%ObsPack_InFile  = TRIM( Input_Opt%ObsPack_InputFile  )
    State_Diag%ObsPack_OutFile = TRIM( Input_Opt%ObsPack_OutputFile )

    ! Replace YYYYMMDD with date and time
    CALL Expand_Date( State_Diag%ObsPack_InFile,  yyyymmdd, hhmmss )
    CALL Expand_Date( State_Diag%ObsPack_OutFile, yyyymmdd, hhmmss )

    ! If we can't find a ObsPack file for today's date, return
    IF ( .NOT. FILE_EXISTS( TRIM( State_Diag%ObsPack_InFile ) ) ) THEN 
       State_Diag%Do_ObsPack = .FALSE.
       State_Diag%ObsPack_nObs = 0
       RETURN
    ENDIF

    !-----------------------------------------------------------------------
    ! Get the list of lon/lat/alt at which to save out GEOS-Chem data
    !-----------------------------------------------------------------------
    CALL ObsPack_Read_Input( Input_Opt, State_Diag, RC )
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Error encountered in "ObsPack_Write_Output"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    !-----------------------------------------------------------------------
    ! Print info about diagnostics that will be saved out
    !-----------------------------------------------------------------------
    IF ( Input_Opt%amIRoot ) THEN

       ! Print info
       WRITE( 6, '(/,a)' ) REPEAT( '=', 79 )
       WRITE( 6, 100     ) yyyymmdd
       WRITE( 6, 110     ) TRIM( State_Diag%ObsPack_InFile  )
       WRITE( 6, 120     ) TRIM( State_Diag%ObsPack_OutFile )
       WRITE( 6, 130     ) State_Diag%ObsPack_nObs
       WRITE( 6, 140     ) State_Diag%ObsPack_nSpecies
       WRITE( 6, '(a,/)' ) REPEAT( '=', 79 )

       ! FORMAT statements
 100   FORMAT( 'OBSPACK for date ',     i8  )
 110   FORMAT( '-> Input file     : ', a   )
 120   FORMAT( '-> Output file    : ', a   )
 130   FORMAT( '-> # observations : ', i10 )
 140   FORMAT( '-> # species/obs  : ', i10 )

    ENDIF

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
  SUBROUTINE ObsPack_Read_Input( Input_Opt, State_Diag, RC )
!
! !USES:
!
    USE CharPak_Mod,            ONLY : CStrip
    USE ErrCode_Mod
    USE File_Mod,               ONLY : File_Exists
    USE Input_Opt_Mod,          ONLY : OptInput 
    USE State_Diag_Mod,         ONLY : DgnState
    USE m_netcdf_io_open,       ONLY : Ncop_Rd
    USE m_netcdf_io_get_dimlen, ONLY : Ncget_Dimlen
    USE m_netcdf_io_read
    USE m_netcdf_io_close,      ONLY : Nccl
    USE m_netcdf_io_checks,     ONLY : NcDoes_Var_Exist
!
! !INPUT PARAMETERS:
!
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
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    LOGICAL              :: It_Exists
    INTEGER              :: fId, N, nObs, nSpecies

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
       CALL ObsPack_Write_Output( Input_Opt, State_Diag, RC )
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Error encountered in "ObsPack_Write_Output"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
    ENDIF

    !=======================================================================
    ! Exit if we are not on the root core
    !
    ! Or for MPI (e.g. WRF-GC), gather inforamation here from all cores
    !=======================================================================
    IF ( .not. Input_Opt%amIRoot ) RETURN

    !=======================================================================
    ! Get the number of observations in the input netCDF file
    !=======================================================================

    ! If the file exists, open it.  Save the file ID to State_Diag.
    CALL Ncop_Rd( fId, State_Diag%ObsPack_InFile )
    State_Diag%ObsPack_fId = fId

    ! Get the # of observations in the file (and save it to State_Diag).
    CALL Ncget_Dimlen( fId, 'obs', nObs )
    State_Diag%ObsPack_nObs = nObs

    !=======================================================================
    ! Allocate the relevant fields of State_Diag
    !=======================================================================
    ALLOCATE( State_Diag%ObsPack_ID( nObs ), STAT=RC ) 
    CALL GC_CheckVar( 'State_Diag%ObsPack_Id', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Diag%ObsPack_Id = ''
   
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

    ALLOCATE( State_Diag%ObsPack_P( nObs ), STAT=RC )
    CALL GC_CheckVar( 'State_Diag%ObsPack_P', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Diag%ObsPack_P = 0.0_f4

    ALLOCATE( State_Diag%ObsPack_T( nObs ), STAT=RC )
    CALL GC_CheckVar( 'State_Diag%ObsPack_T', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Diag%ObsPack_T = 0.0_f4

    nSpecies = State_Diag%ObsPack_nSpecies
    ALLOCATE( State_Diag%ObsPack_Species( nObs, nSpecies ), STAT=RC )
    CALL GC_CheckVar( 'State_Diag%ObsPack_nSpecies', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Diag%ObsPack_Species = 0.0_f4

    ! central_time is a local work array
    ALLOCATE( central_time( 6, nObs ), STAT=RC )
    CALL GC_CheckVar( 'obspack_mod.F:central_time', 0, RC )
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

    ! First check if the "CT_sampling_strategy" variable exists.
    ! If it does not, assume hourly sampling (strategy value = 2).
    varName = 'CT_sampling_strategy'
    It_Exists = NcDoes_Var_Exist( fId, varName )
    IF ( It_Exists ) THEN
       CALL NcRd( State_Diag%ObsPack_Strategy,  fId, TRIM(varName), st1d, ct1d)
    ELSE
       ErrMsg = 'Could not find "CT_sampling_strategy" in file: '         // &
                TRIM( State_Diag%ObsPack_InFile  )                        // &
                '.  Will use hourly sampling by default.'
       CALL GC_Warning( ErrMsg, RC, ThisLoc )
       State_Diag%ObsPack_Strategy = 2
    ENDIF

    !----------------------------
    ! Read ID string
    !----------------------------
    st2d = (/ 1,            1    /)
    ct2d = (/ CHAR_LEN_OBS, nObs /)
    varName = 'obspack_id'
    CALL NcRd( State_Diag%ObsPack_Id, fId, TRIM(varName), st2d, ct2d )

    ! Strip white space (e.g. TABS) from the ObsPack ID strings
    DO N = 1, nObs
       CALL CStrip( State_Diag%ObsPack_Id(N), KeepSpaces=.TRUE. )
    ENDDO

    !----------------------------
    ! Read time components
    ! (Y, M, D, hr, min, sec)
    !----------------------------
    st2d = (/ 1, 1    /)
    ct2d = (/ 6, nObs /)

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
          ! 90-minute window
          !------------------
          CASE( 3 )
             State_Diag%ObsPack_Ival_Start(N) =                              &
                State_Diag%ObsPack_Ival_Center(N) - 2700.0_f8

             State_Diag%ObsPack_Ival_End(N) =                                &
                State_Diag%ObsPack_Ival_Center(N) + 2700.0_f8

          !---------------------
          ! Instaneous Sampling
          !---------------------
          CASE( 4 )
             State_Diag%ObsPack_Ival_Start(N) =                              &
                State_Diag%ObsPack_Ival_Center(N)
                                                              
             State_Diag%ObsPack_Ival_End(N) =                                &
                State_Diag%ObsPack_Ival_Center(N) 


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

    ! Deallocate the central time array
    IF ( ALLOCATED( central_time ) ) THEN
       DEALLOCATE( central_time, STAT=RC )
       CALL GC_CheckVar( 'Central_Time', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

  END SUBROUTINE ObsPack_Read_Input
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ObsPack_Cleanup
!
! !DESCRIPTION: Subroutine ObsPack\_CLEANUP deallocates all ObsPack
!  fields of State_Diag.
!\\
!\\
! !INTERFACE:
!
 SUBROUTINE ObsPack_Cleanup( Input_Opt, State_Diag, RC )
!
! !USES:
!     
   USE ErrCode_Mod
   USE Input_Opt_Mod,  ONLY : OptInput
   USE State_Diag_Mod, ONLY : DgnState
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
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
!  See https://github.com/geoschem/geos-chem for complete history
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

    IF ( ASSOCIATED( State_Diag%ObsPack_P ) ) THEN
       DEALLOCATE( State_Diag%ObsPack_P, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%ObsPack_P', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%ObsPack_P => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%ObsPack_T ) ) THEN
       DEALLOCATE( State_Diag%ObsPack_T, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%ObsPack_T', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%ObsPack_T => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%ObsPack_Species ) ) THEN
       DEALLOCATE( State_Diag%ObsPack_Species, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%ObsPack_Species', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%ObsPack_Species => NULL()
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
!  and writes data to an output file.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE ObsPack_Write_Output( Input_Opt, State_Diag, RC )
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
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES
!
    ! Scalars
    INTEGER               :: fId,        nObs,       nSpecies
    INTEGER               :: ymd,        hms,        omode
    INTEGER               :: yr,         mo,         da
    INTEGER               :: hr,         mn,         sc
    INTEGER               :: vid,        nSamples
    INTEGER               :: N,          S
    INTEGER               :: dId_obs,    dId_spec
    INTEGER               :: dId_obslen, dId_speclen  

    ! Arrays
    INTEGER               :: st1d(1),    ct1d(1),    dims_1d(1)
    INTEGER               :: st2d(2),    ct2d(2),    dims_2d(2)
    INTEGER,  ALLOCATABLE :: aveStart(:)
    INTEGER,  ALLOCATABLE :: aveEnd(:)
    REAL(f8), ALLOCATABLE :: aveTime(:)

    ! Pointers
    REAL(f4), POINTER     :: ptr1d(:)

    ! Strings
    CHARACTER(LEN=16)     :: stamp
    CHARACTER(LEN=31)     :: varName
    CHARACTER(LEN=255)    :: attVal
    CHARACTER(LEN=255)    :: ThisLoc
    CHARACTER(LEN=512)    :: ErrMsg

    !=======================================================================
    ! ObsPack_Write_Output begins here!
    !=======================================================================

    ! Initialize
    RC       =  GC_SUCCESS
    ErrMsg   = ''
    ThisLoc  = ' -> at ObsPack_Write_Output (in module Obspack/obspack_mod.F90)'
    fId      = 0
    vId      = 0
    ptr1d    => NULL()

    ! Exit if no ObsPack observations are found
    nObs = State_Diag%ObsPack_nObs
    IF ( nObs == 0 ) RETURN

    ! Number of species to save out per observation
    nSpecies = State_Diag%ObsPack_nSpecies

    !=======================================================================
    ! Allocate temporary arrays
    !=======================================================================

    ! Cast averaging interval start to integer
    ALLOCATE( aveStart( nObs ), STAT=RC )
    CALL GC_CheckVar( 'obspack_mod.F90:aveStart', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    aveStart = State_Diag%ObsPack_Ival_Start

    ! Cast averaging interval end to integer
    ALLOCATE( aveEnd( nObs ), STAT=RC )
    CALL GC_CheckVar( 'obspack_mod.F90:aveEnd', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    aveEnd = State_Diag%ObsPack_Ival_Start

    ! Compute averaging interval
    ALLOCATE( aveTime( nObs ), STAT=RC )
    CALL GC_CheckVar( 'obspack_mod.F90:aveTime', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    aveTime = State_Diag%ObsPack_Ival_End - State_Diag%ObsPack_Ival_Start
 
    !=======================================================================
    ! Compute averages
    !=======================================================================
    !$OMP PARALLEL DO               &
    !$OMP DEFAULT( SHARED         ) &
    !$OMP PRIVATE( N, S, nSamples )
    DO N = 1, nObs

       ! Number of GEOS-Chem samples for each observation
       nSamples = State_Diag%ObsPack_nSamples(N)

       ! Compute averages of all GEOS-Chem samples for each observation
       ! Only compute averages if there are samples
       IF ( nSamples > 0 ) THEN

          ! Met field quantities
          State_Diag%ObsPack_U(N)   = State_Diag%ObsPack_U(N)   / nSamples
          State_Diag%ObsPack_V(N)   = State_Diag%ObsPack_V(N)   / nSamples
          State_Diag%ObsPack_BLH(N) = State_Diag%ObsPack_BLH(N) / nSamples
          State_Diag%ObsPack_Q(N)   = State_Diag%ObsPack_Q(N)   / nSamples
          State_Diag%ObsPack_T(N)   = State_Diag%ObsPack_T(N)   / nSamples
          State_Diag%ObsPack_P(N)   = State_Diag%ObsPack_P(N)   / nSamples
          
          ! Species concentrations
          DO S = 1, State_Diag%ObsPack_nSpecies
             State_Diag%ObsPack_Species(N,S) =                               &
             State_Diag%ObsPack_Species(N,S) / nSamples
          ENDDO
       ENDIF
    ENDDO
    !$OMP END PARALLEL DO

    !=======================================================================
    ! Open netCDF file for output
    !=======================================================================

    ! Print info
    IF ( Input_Opt%amIRoot ) THEN
       WRITE( 6, 100 ) TRIM( State_Diag%ObsPack_OutFile )
100    FORMAT( '     - OBSPACK: Writing file ', a ) 
    ENDIF

    ! Create netCDF file and save the file ID in State_Diag
    CALL NcCr_Wr( fId, State_Diag%ObsPack_OutFile, WRITE_NC4=.TRUE. )
    State_Diag%ObsPack_fId = fId

    ! Trap potential errors
    IF ( State_Diag%ObsPack_fId < 0 ) THEN
       ErrMsg = 'Invalid netCDF file Id!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    ! Turn filling off
    CALL NcSetFill( fId, NF_NOFILL, omode )

    !=======================================================================
    ! Define dimensions of netCDF variables
    !=======================================================================

    ! Number of observations
    varName = 'obs'
    CALL NcDef_Dimension( fId, TRIM(varName), NF_UNLIMITED,  dId_obs     )

    ! Number of species
    varName = 'species'
    CALL NcDef_Dimension( fId, TRIM(varName), nSpecies,      dId_spec    )

    ! Character length of ObsPack Id strings
    varName = 'char_len_obs'
    CALL NcDef_Dimension( fId, TRIM(varName), CHAR_LEN_OBS,  dId_obslen  )

    !=======================================================================
    ! Set global attributes
    !=======================================================================

    ! History
    stamp   = SYSTEM_TIMESTAMP()
    attVal = 'GEOS-Chem simulation at ' // stamp
    CALL NcDef_Glob_Attributes( fId, 'history', TRIM(attVal) )

    ! Conventions
    CALL NcDef_Glob_Attributes( fId, 'conventions', 'CF-1.4'  )

    ! Reference
    attVal= 'www.geos-chem.org; wiki.geos-chem.org'
    CALL NcDef_Glob_Attributes( fId, 'references', TRIM(attVal)  )

    ! Model start date
    ymd = GET_NYMDb()
    hms = GET_NHMSb()
    CALL Ymd_Extract( ymd, yr, mo, da )
    CALL Ymd_Extract( hms, hr, mn, sc )
    WRITE( attVal, 150 ) yr, mo, da, hr, mn, sc
150 FORMAT( i4.4,"/",i2.2,"/",i2.2," ",i2.2,":",i2.2,":",i2.2, " UTC" )
    CALL NcDef_Glob_Attributes( fId, 'model_start_date',  TRIM(attVal) )

    ! Model end date
    ymd = GET_NYMDe()
    hms = GET_NHMSe()
    CALL Ymd_Extract( ymd, yr, mo, da )
    CALL Ymd_Extract( hms, hr, mn, sc )
    WRITE( attVal, 150 ) yr,mo,da,hr,mn,sc
    CALL NcDef_Glob_Attributes( fId, 'model_end_date', TRIM(attVal) )

    !=======================================================================
    ! Define variables and attributes
    ! NOTE: Dimension order is row-major (i.e. the reverse of Fortran)
    !=======================================================================

    ! Dimension arrays
    dims_1d = (/             dId_obs /)
    dims_2d = (/ dId_obslen, dId_obs /)

    ! ID
    varName = 'obspack_id'
    CALL NcDef_Variable      ( fId, TRIM(varName), NF_CHAR, 2, dims_2d, vId )
    CALL NcDef_Var_Attributes( fId, vId, 'long_name', 'obspack_id'          )
    CALL NcDef_Var_Attributes( fId, vId, 'units',     '1'                   )

    ! # of samples
    varName = 'nsamples'
    CALL NcDef_Variable( fId, TRIM(varName), NF_INT, 1, dims_1d, vId        )
    attVal = 'no. of model samples'
    CALL NcDef_Var_Attributes( fId, vId, 'long_name', TRIM(attVal)          )
    CALL NcDef_Var_Attributes( fId, vId, 'units',     '1'                   )
    attVal = 'Number of discrete model samples in average'
    CALL NcDef_Var_Attributes( fId, vId, 'comment',   TRIM(attVal)          )

    ! Averaging interval
    varName =  'averaging_interval'
    CALL NcDef_Variable( fId, TRIM(varName), NF_INT, 1, dims_1d, vId        )
    attVal = 'Amount of model time over which this observation is averaged'
    CALL NcDef_Var_Attributes( fId, vId, 'long_name', TRIM(attVal)          )
    CALL NcDef_Var_Attributes( fId, vId, 'units',     'seconds'             )

    ! Averaging interval start time
    varName = 'averaging_interval_start'
    CALL NcDef_Variable( fId, TRIM(varName), NF_INT, 1, dims_1d, vId        )
    attVal = 'Start of averaging interval'
    CALL NcDef_Var_Attributes( fId, vId, 'long_name', TRIM(attVal)          )
    attVal = 'seconds since 1970-01-01 00:00:00 UTC'
    CALL NcDef_Var_Attributes( fId, vId, 'units',     TRIM(attVal)          )
    CALL NcDef_Var_Attributes( fId, vId, 'calendar',  'standard'            )

    ! Averaging interval end time
    varName = 'averaging_interval_end'
    CALL NcDef_Variable( fId, TRIM(varName), NF_INT, 1, dims_1d, vId        )
    attVal = 'End of averaging interval'
    CALL NcDef_Var_Attributes( fId, vId, 'long_name', TRIM(attVal)          )
    attVal = 'seconds since 1970-01-01 00:00:00 UTC'
    CALL NcDef_Var_Attributes( fId, vId, 'units',     TRIM(attVal)          )
    CALL NcDef_Var_Attributes( fId, vId, 'calendar',  'standard'            )

    ! Longitude
    CALL NcDef_Variable( fId, 'lon', NF_FLOAT, 1, dims_1d, vId              )
    CALL NcDef_Var_Attributes( fId, vId, 'long_name', 'longitude'           )
    CALL NcDef_Var_Attributes( fId, vId, 'units',     'degrees_east'        )

    ! Latitude
    CALL NcDef_Variable( fId, 'lat', NF_FLOAT, 1, dims_1d, vId              )
    CALL NcDef_Var_Attributes( fId, vId, 'long_name', 'latitude'            )
    CALL NcDef_Var_Attributes( fId, vId, 'units',     'degrees_north'       )

!    ! Altitude
!    CALL NcDef_Variable( fId, 'height', NF_FLOAT, 1, dims_1d, vId           )
!    CALL NcDef_Var_Attributes( fId, vId, 'long_name', 'longitude'           )
!    CALL NcDef_Var_Attributes( fId, vId, 'units',     'degrees_east'        )

    ! U-wind
    CALL NcDef_Variable( fId, 'u', NF_FLOAT, 1, dims_1d, vId                )
    attVal = 'Zonal component of wind'
    CALL NcDef_Var_Attributes( fId, vId, 'long_name', TRIM(attVal)          )
    CALL NcDef_Var_Attributes( fId, vId, 'units',     'm s^-1'              )

    ! V-wind
    CALL NcDef_Variable( fId, 'v', NF_FLOAT, 1, dims_1d, vId                )
    attVal = 'Meridional component of wind'
    CALL NcDef_Var_Attributes( fId, vId, 'long_name', TRIM(attVal)          )
    CALL NcDef_Var_Attributes( fId, vId, 'units',     'm s^-1'              )

    ! Boundary layer height
    CALL NcDef_Variable( fId, 'blh', NF_FLOAT, 1, dims_1d, vId              )
    attVal = 'Boundary layer height'
    CALL NcDef_Var_Attributes( fId, vId, 'long_name', TRIM(attVal)          )
    CALL NcDef_Var_Attributes( fId, vId, 'units',     'm'                   )

    ! Specific humidity
    CALL NcDef_Variable( fId, 'q', NF_FLOAT, 1, dims_1d, vId                )
    attVal = 'mass_fraction_of_water_inair'
    CALL NcDef_Var_Attributes( fId, vId, 'long_name', TRIM(attVal)          )
    attVal = 'kg water (kg air)^-1'
    CALL NcDef_Var_Attributes( fId, vId, 'units',     TRIM(attVal)          )

    ! Pressure
    CALL NcDef_Variable( fId, 'pressure', NF_FLOAT, 1, dims_1d, vId         )
    CALL NcDef_Var_Attributes( fId, vId, 'long_name', 'pressure'            )
    CALL NcDef_Var_Attributes( fId, vId, 'units',     'Pa'                  )

    ! Temperature
    CALL NcDef_Variable( fId, 'temperature', NF_FLOAT, 1, dims_1d, vId      )
    CALL NcDef_Var_Attributes( fId, vId, 'long_name', 'temperature'         )
    CALL NcDef_Var_Attributes( fId, vId, 'units',     'K'                   )

    ! Species concentration
    DO S = 1, State_Diag%ObsPack_nSpecies
       varName  = State_Diag%ObsPack_Species_Name(S)
       attVal = State_Diag%ObsPack_Species_LName(S)
       CALL NcDef_Variable( fId, TRIM(varName), NF_FLOAT, 1, dims_1d, vId   )
       CALL NcDef_Var_Attributes( fId, vId, 'long_name', TRIM(attVal)       )
       CALL NcDef_Var_Attributes( fId, vId, 'units',     'mol mol-1'        )
       CALL NcDef_Var_Attributes( fId, vId, '_FillValue', -1e34             )
       
    ENDDO

    ! End the definition section
    CALL NcEnd_def( fId )
 
    !=======================================================================
    ! Write variables to disk
    ! NOTE: Dimension order is row-major (i.e. the reverse of Fortran)
    !=======================================================================

    !---------------------------------------
    ! Write 2-D variables 
    !---------------------------------------
    st2d    = (/ 1,            1    /)
    ct2d    = (/ CHAR_LEN_OBS, nObs /)

    varName = 'obspack_id'
    CALL NcWr( State_Diag%ObsPack_Id,       fId, TRIM(varName), st2d, ct2d  )

    !---------------------------------------
    ! Write 1-D variables 
    !---------------------------------------
    st1d = (/ 1    /)
    ct1d = (/ nObs /)

    varName = 'nsamples'
    CALL NcWr( State_Diag%ObsPack_nSamples, fId, TRIM(varName), st1d, ct1d  )

    varname = 'averaging_interval'
    CALL NcWr( aveTime,                     fId, TRIM(varName), st1d, ct1d  )

    varname = 'averaging_interval_start'
    CALL NcWr( aveStart,                    fId, TRIM(varName), st1d, ct1d  )

    varname = 'averaging_interval_end'
    CALL NcWr( aveEnd,                      fId, TRIM(varName), st1d, ct1d  )

    varName = 'lon'
    CALL NcWr( State_Diag%ObsPack_Longitude,fId, TRIM(varName), st1d, ct1d  )

    varName = 'lat'
    CALL NcWr( State_Diag%ObsPack_Latitude, fId, TRIM(varName), st1d, ct1d  )

    varName = 'u'
    CALL NcWr( State_Diag%ObsPack_U,        fId, TRIM(varName), st1d, ct1d  )

    varName = 'v'
    CALL NcWr( State_Diag%ObsPack_V,        fId, TRIM(varName), st1d, ct1d  )
    
    varName = 'blh'
    CALL NcWr( State_Diag%ObsPack_BLH,      fId, TRIM(varName), st1d, ct1d  )

    varName = 'q'
    CALL NcWr( State_Diag%ObsPack_Q,        fId, TRIM(varName), st1d, ct1d  )

    varName = 'pressure'
    CALL NcWr( State_Diag%ObsPack_P,        fId, TRIM(varName), st1d, ct1d  )

    varName = 'temperature'
    CALL NcWr( State_Diag%ObsPack_T,        fId, TRIM(varName), st1d, ct1d  )

    !---------------------------------------
    ! Write species concentrations
    !---------------------------------------
    DO S = 1, State_Diag%ObsPack_nSpecies
       varName =  State_Diag%ObsPack_Species_Name(S)
       ptr1d   => State_Diag%ObsPack_Species(:,S)
       CALL NcWr( ptr1d, fId, TRIM(varName), st1d, ct1d  )
       ptr1d   => NULL()
    ENDDO

    ! Close the netCDF file
    CALL NcCl( fId )

    !=======================================================================
    ! Cleanup and quit
    !=======================================================================
    IF ( ALLOCATED( aveStart ) ) THEN
       DEALLOCATE( aveStart, STAT=RC )
       CALL GC_CheckVar( 'obspack_mod.F90:aveStart', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ALLOCATED( aveEnd ) ) THEN
       DEALLOCATE( aveEnd, STAT=RC )
       CALL GC_CheckVar( 'obspack_mod.F90:aveEnd', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ALLOCATED( aveTime ) ) THEN
       DEALLOCATE( aveTime, STAT=RC )
       CALL GC_CheckVar( 'obspack_mod.F90:aveTime', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    CALL ObsPack_Cleanup( Input_Opt, State_Diag, RC )
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
! !IROUTINE: ObsPack_Sample
!
! !DESCRIPTION: Subroutine ObsPack\_Sample performs the model sampling
!  and saves concentrations to locations corresponding to a flight track.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE ObsPack_Sample( yyyymmdd,   hhmmss,     Input_Opt, State_Chm, &
                             State_Diag, State_Grid, State_Met, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE Error_Mod,      ONLY : Debug_Msg
    USE Input_Opt_Mod,  ONLY : OptInput
    USE State_Chm_Mod,  ONLY : ChmState
    USE State_Diag_Mod, ONLY : DgnState
    USE State_Grid_Mod, ONLY : GrdState
    USE State_Met_Mod,  ONLY : MetState
    USE Time_Mod,       ONLY : Ymd_Extract
    USE UnitConv_Mod,   ONLY : Convert_Spc_Units
!
! !INPUT PARAMETERS: 
!
    INTEGER,        INTENT(IN)    :: yyyymmdd    ! Current date
    INTEGER,        INTENT(IN)    :: hhmmss      ! Current time
    TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
    TYPE(GrdState), INTENT(IN)    :: State_Grid  ! Grid State object
    TYPE(MetState), INTENT(IN)    :: State_Met   ! Meteorology State object
!
! !INPUT/OUTPUT PARAMETERS: 
!
    TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
    TYPE(DgnState), INTENT(INOUT) :: State_Diag  ! Diagnostics State object

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
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    LOGICAL             :: prtLog,  prtDebug, doSample
    INTEGER             :: I,       J,        L,  N,  R,  S
    INTEGER             :: Yr,      Mo,       Da, Hr, Mn, Sc
    REAL(f8)            :: TsStart, TsEnd

    ! Strings
    CHARACTER(LEN=255)  :: PriorUnit, ErrMsg, ThisLoc

    !=================================================================
    ! ObsPack_Sample begins here
    !=================================================================

    ! Initialize
    RC       =  GC_SUCCESS
    ErrMsg   = ''
    ThisLoc  = ' -> at ObsPack_Sample (in module ObsPack/obspack_mod.F90)'
    prtLog   = (Input_Opt%amIRoot .and. ( .not. Input_Opt%ObsPack_Quiet ) )
    prtDebug = (Input_Opt%amIRoot .and. Input_Opt%LPRT                    )

    ! Return if ObsPack sampling is turned off (perhaps
    ! because there are no data at this time).
    IF ( .not. State_Diag%Do_ObsPack ) RETURN

    ! Ensure that units of species are "v/v dry", which is dry
    ! air mole fraction.  Capture the InUnit value, this is
    ! what the units are prior to this call.  After we sample
    ! the species, we'll call this again requesting that the
    ! species are converted back to the InUnit values.
    CALL Convert_Spc_Units( Input_Opt, State_Chm, State_Grid,  State_Met, &
                            'v/v dry', RC, PriorUnit )

    ! Trap potential errors
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Error encountered in "Convert_Units"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    !=======================================================================
    ! Sample the GEOS-CHem data corresponding to this location & time
    !=======================================================================

    ! Extract date and time into components
    CALL Ymd_Extract( yyyymmdd, Yr, Mo, Da )
    CALL Ymd_Extract( hhmmss,   Hr, Mn, Sc )
 
    ! Compute elapsed seconds since 1970
    TsEnd   = Seconds_Since_1970( Yr, Mo, Da, Hr, Mn, Sc )
       
    TsStart = TsEnd - Input_Opt%TS_DYN

    ! Logfile header
    IF ( prtLog ) THEN
       WRITE( 6, '(/,a)' ) REPEAT( '=', 79 )
       WRITE( 6, 100     ) Yr, Mo, Da, Hr, Mn, Sc
 100   FORMAT( 'OBSPACK SAMPLING at ', i4,   '/', i2.2, '/', i2.2,           &
               ' ',                    i2.2, ':', i2.2, ':', i2.2           )
       WRITE( 6, '(/,a)' ) ' Obs # ID string'
       WRITE( 6, '(a)'   ) '------ -----------------'

    ENDIF

    ! Loop over observations
    DO N = 1, State_Diag%ObsPack_nObs
       
       !initializing flag for whether sampling should occur at this timestep
       doSample = .false.

       SELECT CASE ( State_Diag%ObsPack_Strategy(N) )
          CASE ( 0 )
             ! Skip observation if the sampling strategy says to do so
             CYCLE
          CASE ( 1:3 )
             ! If the sample covers the entire dynamic timestep, then...
             IF ( State_Diag%ObsPack_Ival_Start(N) <= TsStart .and.                &
                  State_Diag%ObsPack_Ival_End(N)   >= TsEnd ) doSample = .true.
          CASE ( 4 )
             ! If Instantaneous sampling choose the closest timestep
             IF ( (TsEnd - State_Diag%ObsPack_Ival_Center(N)) <= (Input_Opt%TS_DYN/2.) .and.    &
                  (State_Diag%ObsPack_Ival_Center(N) - TsEnd) < (Input_Opt%TS_DYN/2.) ) doSample =.true.
          CASE DEFAULT
             ErrMsg = "Sample Strategy not implemented in ObsPack_Sample Subroutine"
             CALL GC_Error( ErrMsg, RC, ThisLoc )     
             RETURN
       END SELECT
                
                

       ! If sampling strategy time-step conditions are met, sample at these times
       IF ( doSample ) THEN 
          
          ! Print the observations that are sampled here
          IF ( prtLog ) THEN
             WRITE( 6, '(i6,1x,a)' ) N, TRIM( State_Diag%ObsPack_Id(N) )
          ENDIF

          ! Return grid box indices for the chemistry region
          CALL ObsPack_Get_Indices( N, State_Diag, State_Grid, State_Met,    &
                                    I, J, L, RC )

          ! Trap potential errors
          IF ( RC /= GC_SUCCESS ) THEN
             ErrMsg = 'Error encountered in "ObsPack_Get_Indices"!'
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN
          ENDIF

          !------------------------
          ! Update sample counter
          !------------------------
          State_Diag%ObsPack_nSamples(N) = State_Diag%ObsPack_nSamples(N) + 1

          !-----------------------
          ! Archive species
          !-----------------------

          ! Loop over the # of requested species
          DO R = 1, State_Diag%ObsPack_nSpecies

             ! Get the proper index for State_Chm%Species
             S = State_Diag%ObsPack_Species_Ind(R)
             
             ! Add the species concentration
             !%%% TODO: Maybe do an in-line unit conversion to avoid
             !%%% converting units for all species if we are only saving
             !%%% out a few.
             State_Diag%ObsPack_Species(N,R) =                               &
             State_Diag%ObsPack_Species(N,R) + State_Chm%SPECIES(I,J,L,S)
          ENDDO

          !-----------------------
          ! Archive met fields
          !-----------------------
          State_Diag%ObsPack_U(N)   = State_Diag%ObsPack_U(N)                &
                                    + State_Met%U(I,J,L) 

          State_Diag%ObsPack_V(N)   = State_Diag%ObsPack_V(N)                &
                                    + State_Met%V(I,J,L) 

          State_Diag%ObsPack_BLH(N) = State_Diag%ObsPack_BLH(N)              &
                                    + State_Met%PBLH(I,J) 

          State_Diag%ObsPack_Q(N)   = State_Diag%ObsPack_Q(N)                &
                                    + State_Met%SPHU(I,J,L) 

          State_Diag%ObsPack_P(N)   = State_Diag%ObsPack_P(N)                &
                                    + State_Met%PMID(I,J,L) 

          State_Diag%ObsPack_T(N)   = State_Diag%ObsPack_T(N)                &
                                    + State_Met%T(I,J,L) 

       ENDIF
    ENDDO

    ! Logfile footer
    IF ( prtLog ) THEN
       WRITE( 6, '(a,/)' ) REPEAT( '=', 79 )
    ENDIF

    !=======================================================================
    ! Cleanup and quit
    !=======================================================================

    ! Return State_Chm%SPECIES to whatever units they had
    ! coming into this routine
    call Convert_Spc_Units( Input_Opt, State_Chm, State_Grid, State_Met, &
                            PriorUnit, RC )

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
! !IROUTINE: ObsPack_Get_Grid_Indices
!
! !DESCRIPTION: Subroutine ObsPack\_GET\_GRID\_INDICES returns the
!  grid box indices (I, J, L) corresponding to the input point
!  defined by longitude, latitude, altitude.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE ObsPack_Get_Indices( iObs, State_Diag, State_Grid, State_Met, &
                                  I, J, L, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE State_Diag_Mod, ONLY : DgnState 
    USE State_Grid_Mod, ONLY : GrdState 
    USE State_Met_Mod,  ONLY : MetState
!
! !INPUT PARAMETERS: 
!
    INTEGER,        INTENT(IN)  :: iObs         ! ObsPack Observation number
    TYPE(DgnState), INTENT(IN)  :: State_Diag   ! Diagnostics State object
    TYPE(GrdState), INTENT(IN)  :: State_Grid   ! Grid State object
    TYPE(MetState), INTENT(IN)  :: State_Met    ! Meteorology State object
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
!  See https://github.com/geoschem/geos-chem for complete history
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
    I0 = State_Grid%XMinOffset
    J0 = State_Grid%YMinOffset

    !-----------------------------------------------------------------------
    ! Find I corresponding to the ObsPack longitude value
    !-----------------------------------------------------------------------
    I = INT( ( State_Diag%ObsPack_Longitude(iObs) + 180.0_f8                 &
                                                  - ( I0 * State_Grid%DX ) ) &
                                                  / State_Grid%DX + 1.5d0  )

    ! Handle date line correctly (bmy, 4/23/04)
    IF ( I > State_Grid%NX ) I = I - State_Grid%NX

    !-----------------------------------------------------------------------
    ! Find J corresponding to the ObsPack latitude value
    !-----------------------------------------------------------------------
    J = INT( ( State_Diag%ObsPack_Latitude(iObs)  +  90.0_f8                 &
                                                  - ( J0 * State_Grid%DY ) ) &
                                                  / State_Grid%DY + 1.5d0  )

    !-----------------------------------------------------------------------
    ! Find L corresponding to the ObsPack altitude value
    !-----------------------------------------------------------------------
    Z = 0.0_f8
    DO L = 1, State_Grid%NZ
       Z = Z + State_Met%BXHEIGHT(I,J,L)
       IF ( Z >= State_Diag%ObsPack_Altitude(iObs) ) RETURN 
    ENDDO

    !========================================================================
    ! Issue error message if we get here.
    !========================================================================
    !WRITE (6,*) 'At I,J =',i,j
    !FLUSH(6)

    DO L = 1,State_Grid%NZ
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
!  See https://github.com/geoschem/geos-chem for complete history
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
                            ( DBLE( Minute ) /  1440.0_f8 )  +               &
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
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ObsPack_SpeciesMap_Init
!
! !DESCRIPTION: Gets the index, name, and long name for each requested
!  ObsPack species from the GEOS-Chem Species Database.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE ObsPack_SpeciesMap_Init( Input_Opt, State_Chm, State_Diag, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE Input_Opt_Mod,  ONLY : OptInput
    USE Species_Mod,    ONLY : Species
    USE State_Chm_Mod,  ONLY : ChmState, Ind_
    USE State_Diag_Mod, ONLY : DgnState
!
! !INPUT PARAMETERS: 
!
    TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
!
! !INPUT/OUTPUT PARAMETERS: 
!
    TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
    TYPE(DgnState), INTENT(INOUT) :: State_Diag  ! Diagnostics State object
!
! !OUTPUT PARAMETERS: 
!
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
!
! !REMARKS:
!  This routine can be called at the model initialization stage.
!  (unlike ObsPack_Init, which is called every time a new input file is read).
!
! !REVISION HISTORY:
!  04 Jan 2019 - R. Yantosca - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    LOGICAL                :: WildCard 
    INTEGER                :: N,       nSpc,     R

    ! Strings
    CHARACTER(LEN=255)     :: ErrMsg,  ThisLoc,  SpcName

    ! Objects
    TYPE(Species), POINTER :: ThisSpc

    !=======================================================================
    ! ObsPack_SpeciesMap_Innit begins here!
    !=======================================================================
 
    ! Initialize
    RC      = GC_SUCCESS
    ErrMsg  =''
    ThisLoc =                                                                &
     ' -> at ObsPack_SpeciesMap_Init (in module Obspack/obspack_mod.F90)'
 
    !=======================================================================
    ! Allocate appropriate fields of State_Diag 
    !=======================================================================

    ! Number of species that will be saved to ObsPack output
    ! (accounting for any wildcards)
    IF ( TRIM( Input_Opt%ObsPack_SpcName(1) ) == '?ADV?' ) THEN
       State_Diag%ObsPack_nSpecies = State_Chm%nAdvect
       WildCard                    = .TRUE.
    ELSE
       State_Diag%ObsPack_nSpecies = Input_Opt%ObsPack_nSpc
       WildCard                    = .FALSE. 
    ENDIF

    ! Save in a shadow variable
    nSpc = State_Diag%ObsPack_nSpecies

    ! Trap potential errors
    IF ( nSpc < 1 ) THEN
       ErrMsg = 'OBSPACK: No valid species found!  Exiting...'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    ! ObsPack species index
    ALLOCATE( State_Diag%ObsPack_Species_Ind( nSpc ), STAT=RC ) 
    CALL GC_CheckVar( 'State_Diag%ObsPack_Species_Ind', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Diag%ObsPack_Species_Ind = 0

    ! ObsPack species name
    ALLOCATE( State_Diag%ObsPack_Species_Name( nSpc ), STAT=RC ) 
    CALL GC_CheckVar( 'State_Diag%ObsPack_Species_Name', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Diag%ObsPack_Species_Name = ''

    ! ObsPack species long-name
    ALLOCATE( State_Diag%ObsPack_Species_LName( nSpc ), STAT=RC ) 
    CALL GC_CheckVar( 'State_Diag%ObsPack_Species_LName', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Diag%ObsPack_Species_LName = ''

    !=======================================================================
    ! Map species: for each ObsPack requested species, find the 
    ! corresponding entry in the GEOS-Chem Species Database.
    !=======================================================================
    IF ( WildCard ) THEN

       !--------------------------------------------------------------
       ! Take all advected species if the ?ADV? wildcard is found
       !--------------------------------------------------------------

       ! Loop over all advected species
       DO R = 1, State_Chm%nAdvect

          ! Get the overall species index
          N = State_Chm%Map_Advect(R)

          ! Point to the relevant entry in the species database
          ThisSpc => State_Chm%SpcData(N)%Info

          ! Save metadata from species database to ObsPack arrays
          State_Diag%ObsPack_Species_Ind(R)   = N
          State_Diag%ObsPack_Species_Name(R)  = TRIM( ThisSpc%Name     )
          State_Diag%ObsPack_Species_LName(R) = TRIM( ThisSpc%FullName )

          ! Free pointer
          ThisSpc => NULL()
       ENDDO
       
    ELSE

       !--------------------------------------------------------------
       ! Otherwise look up individual species names
       !--------------------------------------------------------------
       DO R = 1, State_Diag%ObsPack_nSpecies

          ! Get the species index 
          N = Ind_( Input_Opt%ObsPack_SpcName(R) )

          ! If the species isn't found, exit with error
          IF ( N < 0 ) THEN 
             ErrMsg = 'Could not find species database info for the '     // &
                      'ObsPack requested species: '                       // &
                      TRIM( Input_Opt%ObsPack_SpcName(R) )
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN
          ENDIF

          ! Point to the relevant entry in the species database
          ThisSpc => State_Chm%SpcData(N)%Info

          ! Save metadata from species database to ObsPack arrays
          State_Diag%ObsPack_Species_Ind(R)   = N
          State_Diag%ObsPack_Species_Name(R)  = TRIM( ThisSpc%Name     )
          State_Diag%ObsPack_Species_LName(R) = TRIM( ThisSpc%FullName )

          ! Free pointer
          ThisSpc => NULL()

       ENDDO
    ENDIF

    !=======================================================================
    ! Print output of Obspack requested species names
    !=======================================================================
    IF ( Input_Opt%amIRoot .and. ( .not. Input_Opt%ObsPack_Quiet ) ) THEN

       ! Header
       WRITE( 6, '(/,a)' ) REPEAT( '=', 79 )
       WRITE( 6, '(a,/)' ) 'OBSPACK: Species that will be saved as output'
       WRITE( 6, 100     ) 
       WRITE( 6, 110     )

       ! Species info
       DO R = 1, State_Diag%ObsPack_nSpecies
          WRITE( 6, 120 ) State_Diag%ObsPack_Species_Ind(R),                 &
                          State_Diag%ObsPack_Species_Name(R)(1:31),          &
                          State_Diag%ObsPack_Species_Lname(R)(1:45)
       ENDDO

       ! FORMAT statements
 100   FORMAT( ' Index  Name', 29x, 'Long name' )
 110   FORMAT( '------  ----', 29x, '---------' )   
 120   FORMAT( i6, 2x, a31, 2x, a40 )

       ! Footer
       WRITE( 6, '(a,/)'   ) REPEAT( '=', 79 )
    ENDIF

  END SUBROUTINE ObsPack_SpeciesMap_Init
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ObsPack_SpeciesMap_Cleanup
!
! !DESCRIPTION: Deallocates and frees fields of State\_Diag that are
!  relevant to the ObsPack species mapping.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE ObsPack_SpeciesMap_Cleanup( Input_Opt, State_Diag, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE Input_Opt_Mod,  ONLY : OptInput
    USE State_Chm_Mod,  ONLY : ChmState
    USE State_Diag_Mod, ONLY : DgnState
!
! !INPUT PARAMETERS: 
!
    TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
!
! !INPUT/OUTPUT PARAMETERS: 
!
    TYPE(DgnState), INTENT(INOUT) :: State_Diag  ! Diagnostics State object
!
! !OUTPUT PARAMETERS: 
!
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
!
! !REMARKS:
!  This routine can be called at the overall model finalization stage
!  (unlike ObsPack_Cleanup, which is called every time a file is written).
!
! !REVISION HISTORY:
!  04 Jan 2019 - R. Yantosca - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
    
    !=======================================================================
    ! ObsPack_SpeciesMap_Cleanup begins here!
    !=======================================================================
    
    ! Assume success
    RC = GC_SUCCESS
    
    ! Deallocate arrays 
    IF ( ASSOCIATED( State_Diag%ObsPack_Species_Ind ) ) THEN
       DEALLOCATE( State_Diag%ObsPack_Species_Ind, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%ObsPack_Species_Ind', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%ObsPack_Species_Ind => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%ObsPack_Species_Name ) ) THEN
       DEALLOCATE( State_Diag%ObsPack_Species_Name, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%ObsPack_Species_Name', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%ObsPack_Species_Name => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Diag%ObsPack_Species_Lname ) ) THEN
       DEALLOCATE( State_Diag%ObsPack_Species_LName, STAT=RC )
       CALL GC_CheckVar( 'State_Diag%ObsPack_Species_LName', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%ObsPack_Species_Name => NULL()
    ENDIF
   
  END SUBROUTINE ObsPack_SpeciesMap_Cleanup
!EOC
END MODULE ObsPack_Mod


