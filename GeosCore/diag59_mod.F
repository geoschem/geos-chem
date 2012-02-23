!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: diag59_mod
!
! !DESCRIPTION: Module DIAG59\_MOD contains variables and routines to save 
!  out the fraction of NOx remaining and integrated OPE to disk
!  (gvinken, 25/02/11) 
!\\
!\\
! !INTERFACE: 
!
      MODULE DIAG59_MOD
!
! !USES:
!
      IMPLICIT NONE
      PRIVATE
!
! !PUBLIC DATA MEMBERS:
!
      LOGICAL, PUBLIC :: DO_SAVE_DIAG59 
!
! !PUBLIC MEMBER FUNCTIONS:
!
      PUBLIC :: DIAG59
      PUBLIC :: ITS_TIME_FOR_DIAG59
      PUBLIC :: INIT_DIAG59
!
! !REMARKS:
!  ND59 tracer numbers:
!  ============================================================================
!  1  : Fraction of NOx remaining  [unitless]
!  2  : Integrated OPE             [molec O3 produced / molec NOx lost]
!
! !REVISION HISTORY:
!  25 Feb 2011 - G. Vinken   - Initial version based on the orig. diag49_mod.f
!  07 Feb 2012 - M. Payer    - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !PRIVATE TYPES:
!
      !=================================================================
      ! MODULE VARIABLES
      !
      ! DO_SAVE_DIAG59   : Switch to turn ND59 timeseries on/off 
      ! I0               : Lon offset between global & nested grid
      ! J0               : Lat offset between global & nested grid
      ! IOFF             : Offset between relative & absolute lon
      ! JOFF             : Offset between relative & absolute lat
      ! ND59_IMIN        : Minimum longitude index
      ! ND59_IMAX        : Maximum latitude  index
      ! ND59_JMIN        : Minimum longitude index
      ! ND59_JMAX        : Maximum longitude index
      ! ND59_FREQ        : Frequency which to save to disk [min]
      ! ND59_N_TRACERS   : Number of tracers for ND59 timeseries
      ! ND59_OUTPUT_FILE : Name of timeseries output file
      ! ND59_TRACERS     : Array w/ tracer #'s to save to disk
      ! HALFPOLAR        : Used for binary punch file write
      ! CENTER180        : Used for binary punch file write
      ! LONRES           : Used for binary punch file write
      ! LATRES           : Used for binary punch file write
      ! RESERVED         : Used for binary punch file write
      ! MODELNAME        : Used for binary punch file write
      ! TITLE            : Used for binary punch file write 
      !=================================================================

      INTEGER            :: IOFF,           JOFF
      INTEGER            :: I0,             J0
      INTEGER            :: ND59_N_TRACERS, ND59_TRACERS(120)
      INTEGER            :: ND59_IMIN,      ND59_IMAX
      INTEGER            :: ND59_JMIN,      ND59_JMAX
      INTEGER            :: ND59_FREQ,      ND59_NI
      INTEGER            :: ND59_NJ,        ND59_NL
      INTEGER            :: HALFPOLAR
      INTEGER, PARAMETER :: CENTER180=1 
      REAL*4             :: LONRES,         LATRES
      CHARACTER(LEN=20)  :: MODELNAME
      CHARACTER(LEN=40)  :: RESERVED = ''
      CHARACTER(LEN=80)  :: TITLE
      CHARACTER(LEN=255) :: ND59_OUTPUT_FILE

      CONTAINS
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: diag59 
!
! !DESCRIPTION: Subroutine DIAG59 produces time series (instantaneous fields) 
!  for a geographical domain from the information read in timeseries.dat.  
!  Output will be in binary punch (BPCH) format.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE DIAG59
!
! !USES:
!
      USE BPCH2_MOD,    ONLY : BPCH2,   OPEN_BPCH2_FOR_WRITE
      USE DAO_MOD,      ONLY : AD,      AIRDEN, BXHEIGHT, CLDF 
      USE DAO_MOD,      ONLY : CLDTOPS, OPTD,   RH,       SLP     
      USE DAO_MOD,      ONLY : T,       UWND,   VWND
      USE DAO_MOD,      ONLY : TS
      USE DAO_MOD,      ONLY : PARDF, PARDR
      USE LAI_MOD,      ONLY : ISOLAI
      USE FILE_MOD,     ONLY : IU_ND59
      USE GRID_MOD,     ONLY : GET_XOFFSET,        GET_YOFFSET
      USE TIME_MOD,     ONLY : EXPAND_DATE
      USE TIME_MOD,     ONLY : GET_NYMD,           GET_NHMS
      USE TIME_MOD,     ONLY : GET_NYMD_DIAG,      GET_TS_DIAG
      USE TIME_MOD,     ONLY : GET_TAU,            GET_HOUR
      USE TIME_MOD,     ONLY : ITS_A_NEW_DAY,      TIMESTAMP_STRING
      USE PBL_MIX_MOD,  ONLY : GET_PBL_TOP_L,      GET_PBL_TOP_m
      USE TRACER_MOD,   ONLY : ITS_A_FULLCHEM_SIM, N_TRACERS
      USE TRACER_MOD,   ONLY : STT,                TCVV
      USE TRACER_MOD,   ONLY : XNUMOLAIR
      USE PRESSURE_MOD, ONLY : GET_PEDGE
      USE TRACERID_MOD, ONLY : IDTHNO3, IDTHNO4, IDTN2O5, IDTNOX  
      USE TRACERID_MOD, ONLY : IDTPAN,  IDTPMN,  IDTPPN,  IDTOX   
      USE TRACERID_MOD, ONLY : IDTR4N2, IDTSALA, IDTSALC 
      USE DIAG_MOD,     ONLY : AD59     

      USE CMN_FJ_MOD         ! FAST-J stuff, includes CMN_SIZE
      USE JV_CMN_MOD         ! ODAER, QAA, QAA_AOD (clh)
      USE CMN_O3_MOD	     ! Pure O3, SAVENO2
      USE CMN_GCTM_MOD       ! XTRA2
!
! !REVISION HISTORY:
!  25 Feb 2011 - G. Vinken   - Initial version based on DIAG49 
!  07 Feb 2012 - M. Payer    - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL, SAVE            :: FIRST  = .TRUE.
      LOGICAL, SAVE            :: IS_FULLCHEM, IS_NOx,     IS_Ox 
      LOGICAL, SAVE            :: IS_NOy,      IS_CLDTOPS, IS_OPTD
      LOGICAL, SAVE            :: IS_SEASALT,  IS_SLP
      INTEGER                  :: IOS,  GMTRC, GMNL, I, J, K, L 
      INTEGER                  :: N,    R,     H,    W, X, Y
      INTEGER                  :: NHMS, TS_DIAG
      REAL*8                   :: TAU, TMP,   SCALEAODnm
      REAL*8                   :: Q( ND59_NI, ND59_NJ, 1)
      CHARACTER(LEN=16)        :: STAMP
      CHARACTER(LEN=40)        :: CATEGORY
      CHARACTER(LEN=40)        :: UNIT
      CHARACTER(LEN=255)       :: FILENAME

      !=================================================================
      ! DIAG59 begins here!
      !=================================================================

      ! Set logical flags on first timestep
      IF ( FIRST ) THEN
         IS_CLDTOPS  = ALLOCATED( CLDTOPS )
         IS_OPTD     = ALLOCATED( OPTD    )
         IS_SLP      = ALLOCATED( SLP     )
         IS_FULLCHEM = ITS_A_FULLCHEM_SIM()
         IS_SEASALT  = ( IDTSALA > 0 .and. IDTSALC > 0 )
         IS_Ox       = ( IS_FULLCHEM .and. IDTOX   > 0 )
         IS_NOx      = ( IS_FULLCHEM .and. IDTNOX  > 0 )
         IS_NOy      = ( IS_FULLCHEM .and. 
     &                   IDTNOX  > 0 .and. IDTPAN  > 0 .and.
     &                   IDTHNO3 > 0 .and. IDTPMN  > 0 .and.
     &                   IDTPPN  > 0 .and. IDTR4N2 > 0 .and.
     &                   IDTN2O5 > 0 .and. IDTHNO4 > 0 ) 
         FIRST       = .FALSE.
      ENDIF

      !=================================================================
      ! If it's a new day, open a new BPCH file and write file header
      ! We need to check if it's a new day + 1 ND59 time step (ccc, 8/12/09)
      !=================================================================
!--- Previous to (ccc, 8/12/09)
!      IF ( ITS_A_NEW_DAY() ) THEN
      NHMS    = GET_NHMS()
      TS_DIAG = ND59_FREQ

      ! To change TS_DIAG to NHMS format
      TS_DIAG = TS_DIAG/60 * 10000 + (TS_DIAG - (TS_DIAG/60)*60) * 100  

      IF ( NHMS == TS_DIAG ) THEN     ! It's a new day for diagnostics.

         ! Expand date tokens in the file name
         FILENAME = TRIM( ND59_OUTPUT_FILE )
         CALL EXPAND_DATE( FILENAME, GET_NYMD(), GET_NHMS() )

         ! Echo info
         WRITE( 6, 100 ) TRIM( FILENAME )
 100     FORMAT( '     - DIAG59: Opening file ', a )
        
         ! Open bpch file and write top-of-file header
         CALL OPEN_BPCH2_FOR_WRITE( IU_ND59, FILENAME, TITLE )
      ENDIF

      !=================================================================
      ! Save tracers to timeseries file
      !=================================================================

      ! Echo info
      STAMP = TIMESTAMP_STRING()
      WRITE( 6, 110 ) STAMP
 110  FORMAT( '     - DIAG59: Saving timeseries at ', a )

      ! Time for BPCH file
      TAU  = GET_TAU()

      ! Loop over tracers
      DO W = 1, ND59_N_TRACERS

         ! Zero summing array
         Q = 0d0

         !-------------------------------------
         ! SHIP-diagnostic tracers [unitless]
         !-------------------------------------
         CATEGORY = 'SHIP-$$$'
         UNIT     = ''           ! Let GAMAP pick the unit
         GMNL     = ND59_NL
         GMTRC    = W
         
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, X, Y)
            DO Y = 1, ND59_NJ
               J = JOFF + Y
            DO X = 1, ND59_NI
               I = GET_I( X )
               Q(X,Y,1) = AD59(I,J,W)
            ENDDO
            ENDDO
!$OMP END PARALLEL DO

         !==============================================================
         ! Save this data block to the ND59 timeseries file
         !==============================================================
         CALL BPCH2( IU_ND59,      MODELNAME,    LONRES,   
     &               LATRES,       HALFPOLAR,    CENTER180, 
     &               CATEGORY,     GMTRC,        UNIT,      
     &               TAU,          TAU,          RESERVED,  
     &               ND59_NI,      ND59_NJ,      1,  
     &               ND59_IMIN+I0, ND59_JMIN+J0, 1, 
     &               REAL( Q(1:ND59_NI, 1:ND59_NJ, 1) ) )
      ENDDO
            
      !=================================================================
      ! Close the file at the proper time
      !=================================================================
      IF ( ITS_TIME_TO_CLOSE_FILE() ) THEN

         ! Expand date tokens in the file name
         FILENAME = TRIM( ND59_OUTPUT_FILE )
!--- Previous to (ccc, 8/12/09)
!         CALL EXPAND_DATE( FILENAME, GET_NYMD(), GET_NHMS() )
         CALL EXPAND_DATE( FILENAME, GET_NYMD_DIAG(), GET_NHMS() )

         ! Echo info
         WRITE( 6, 120 ) TRIM( FILENAME )
 120     FORMAT( '     - DIAG59: Closing file : ', a )

         ! Close file
         CLOSE( IU_ND59 ) 
      ENDIF

      END SUBROUTINE DIAG59
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: its_time_to_close_file
!
! !DESCRIPTION: Function ITS\_TIME\_TO\_CLOSE\_FILE returns TRUE if it's 
!  time to close the ND59 bpch file before the end of the day.
!\\
!\\
! !INTERFACE:
!
      FUNCTION ITS_TIME_TO_CLOSE_FILE() RESULT( ITS_TIME )
!
! !USES:
!
      USE TIME_MOD, ONLY : GET_HOUR
      USE TIME_MOD, ONLY : GET_MINUTE
! 
! !REVISION HISTORY: 
!  20 Jul 2004 - R. Yantosca - Initial version
!  (1 ) The time is already updated to the next time step (ccc, 8/12/09)
!  02 Dec 2010 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL :: ITS_TIME
      REAL*8  :: HR1

      !=================================================================
      ! ITS_TIME_TO_CLOSE_FILE begins here!
      !=================================================================

      ! Current hour
      HR1      = GET_HOUR() + ( GET_MINUTE() / 60d0 )

!--- Previous to (ccc, 8/12/09)
!      ! Hour at the next dynamic timestep
!      HR2      = HR1        + ( ND59_FREQ / 60d0 )

       ! If the next dyn step is the start of a new day, return TRUE
!--- Previous to (ccc, 11/11/10)
!       HR1 varies between 00 and 23:59. So compares to 00 not 24 anymore.
!      ITS_TIME = ( INT( HR1 ) == 24 )
      ITS_TIME = ( INT( HR1 ) == 00 )
 
       END FUNCTION ITS_TIME_TO_CLOSE_FILE
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: its_time_for_diag59
!
! !DESCRIPTION: Function ITS\_TIME\_FOR\_DIAG59 returns TRUE if ND59 is 
!  turned on and it is time to call DIAG59 -- or FALSE otherwise.
!\\
!\\
! !INTERFACE:
!
      FUNCTION ITS_TIME_FOR_DIAG59() RESULT( ITS_TIME )
!
! !USES:
!
      USE TIME_MOD,  ONLY : GET_ELAPSED_MIN
      USE TIME_MOD,  ONLY : GET_TS_DIAG
      USE ERROR_MOD, ONLY : GEOS_CHEM_STOP
!
! !RETURN VALUE:
!
      LOGICAL :: ITS_TIME
! 
! !REVISION HISTORY: 
!  25 Feb 2011 - G. Vinken   - Initial version based on ITS_TIME_FOR_DIAG49
!  07 Feb 2012 - M. Payer    - Added ProTeX headers 
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER       :: XMIN, TS_DIAG
      LOGICAL, SAVE :: FIRST = .TRUE.

      !=================================================================
      ! ITS_TIME_FOR_DIAG59 begins here!
      !=================================================================

      IF ( DO_SAVE_DIAG59 ) THEN
         IF ( FIRST ) THEN
            TS_DIAG = GET_TS_DIAG()
            
            ! Check if ND59_FREQ is a multiple of TS_DIAG
            IF ( MOD( ND59_FREQ, TS_DIAG ) /= 0 ) THEN
               WRITE( 6, 100 ) 'ND59', ND59_FREQ, TS_DIAG
 100           FORMAT( 'The ',a,' output frequency must be a multiple '
     &              'of the largest time step:', i5, i5 )
               CALL GEOS_CHEM_STOP
            ENDIF
            FIRST = .FALSE.
         ENDIF
         
         ! Time already elapsed in this run
         XMIN     = GET_ELAPSED_MIN()
         
         ! Is the elapsed time a multiple of ND59_FREQ?
         ITS_TIME = ( DO_SAVE_DIAG59 .and. MOD( XMIN, ND59_FREQ ) == 0 )
      ELSE
         ITS_TIME = DO_SAVE_DIAG59
      ENDIF
            
      END FUNCTION ITS_TIME_FOR_DIAG59
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_i
!
! !DESCRIPTION: Function GET\_I returns the absolute longitude index (I), 
!  given the relative longitude index (X).
!\\
!\\
! !INTERFACE:
!

      FUNCTION GET_I( X ) RESULT( I )
!
! !USES:
!
      USE CMN_SIZE_MOD         ! Size parameters
!
! !INPUT PARAMETERS: 
!
      INTEGER, INTENT(IN) :: X   ! Relative longitude index (used by Q array)
!
! !RETURN VALUE:
!
      INTEGER             :: I   ! Absolute longitude index
! 
! !REVISION HISTORY: 
!  20 Jul 2004 - R. Yantosca - Initial version
!  02 Dec 2010 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC      
      !=================================================================
      ! GET_I begins here!
      !=================================================================

      ! Add the offset to X to get I  
      I = IOFF + X

      ! Handle wrapping around the date line, if necessary
      IF ( I > IIPAR ) I = I - IIPAR

      END FUNCTION GET_I
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_diag59
!
! !DESCRIPTION: Subroutine INIT\_DIAG59 allocates and zeroes all module 
!  arrays.  It also gets values for module variables from "input\_mod.f". 
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE INIT_DIAG59( DO_ND59, N_ND59, TRACERS, IMIN,    
     &                        IMAX,    JMIN,   JMAX,    FREQ,   
     &                        FILE )
!
! !USES:
!
      USE BPCH2_MOD, ONLY : GET_MODELNAME
      USE BPCH2_MOD, ONLY : GET_HALFPOLAR
      USE GRID_MOD,  ONLY : GET_XOFFSET
      USE GRID_MOD,  ONLY : GET_YOFFSET
      USE GRID_MOD,  ONLY : ITS_A_NESTED_GRID
      USE ERROR_MOD, ONLY : ERROR_STOP

      USE CMN_SIZE_MOD    ! Size parameters
!
! !INPUT PARAMETERS: 
!
      ! DO_ND59 : Switch to turn on ND59 timeseries diagnostic
      ! N_ND59  : Number of ND59 read by "input_mod.f"
      ! TRACERS : Array w/ ND59 tracer #'s read by "input_mod.f"
      ! IMIN    : Min longitude index read by "input_mod.f"
      ! IMAX    : Max longitude index read by "input_mod.f" 
      ! JMIN    : Min latitude index read by "input_mod.f" 
      ! JMAX    : Min latitude index read by "input_mod.f" 
      ! FREQ    : Frequency for saving to disk [min]
      ! FILE    : ND59 output file name read by "input_mod.f"
      LOGICAL,            INTENT(IN) :: DO_ND59
      INTEGER,            INTENT(IN) :: N_ND59, TRACERS(100)
      INTEGER,            INTENT(IN) :: IMIN,   IMAX 
      INTEGER,            INTENT(IN) :: JMIN,   JMAX      
      INTEGER,            INTENT(IN) :: FREQ
      CHARACTER(LEN=255), INTENT(IN) :: FILE
! 
! !REVISION HISTORY:
!  25 Feb 2011 - G. Vinken   - Initial version based on INIT_DIAG49
!  07 Feb 2012 - M. Payer    - Added ProTeX headers 
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      CHARACTER(LEN=255)             :: LOCATION
      
      !=================================================================
      ! INIT_DIAG59 begins here!
      !=================================================================

      ! Initialize
      LOCATION               = 'INIT_DIAG59 ("diag59_mod.f")'
      ND59_TRACERS(:)        = 0

      ! Get values from "input_mod.f"
      DO_SAVE_DIAG59         = DO_ND59 
      ND59_N_TRACERS         = N_ND59
      ND59_TRACERS(1:N_ND59) = TRACERS(1:N_ND59)
      ND59_IMIN              = IMIN
      ND59_IMAX              = IMAX
      ND59_JMIN              = JMIN
      ND59_JMAX              = JMAX
      ND59_FREQ              = FREQ
      ND59_OUTPUT_FILE       = FILE
     
      ! Return if we are not saving ND59 diagnostics
      IF ( .not. DO_SAVE_DIAG59 ) RETURN

      !=================================================================
      ! Compute lon, lat, alt extents and check for errors
      !=================================================================

      ! Get grid offsets for error checking
      IF ( ITS_A_NESTED_GRID() ) THEN
         I0 = GET_XOFFSET()
         J0 = GET_YOFFSET()
      ELSE
         I0 = GET_XOFFSET( GLOBAL=.TRUE. )
         J0 = GET_YOFFSET( GLOBAL=.TRUE. )
      ENDIF

      !-----------
      ! Longitude
      !-----------

      ! Error check ND59_IMIN
      IF ( ND59_IMIN+I0 < 1 .or. ND59_IMIN+I0 > IGLOB ) THEN
         CALL ERROR_STOP( 'Bad ND59_IMIN value!', LOCATION )
      ENDIF

      ! Error check ND59_IMAX
      IF ( ND59_IMAX+I0 < 1 .or. ND59_IMAX+I0 > IGLOB ) THEN
         CALL ERROR_STOP( 'Bad ND59_IMAX value!', LOCATION )
      ENDIF

      ! Compute longitude limits to write to disk 
      ! Also handle wrapping around the date line
      IF ( ND59_IMAX >= ND59_IMIN ) THEN
         ND59_NI = ( ND59_IMAX - ND59_IMIN ) + 1
      ELSE 
         ND59_NI = ( IIPAR - ND59_IMIN ) + 1 + ND59_IMAX
         WRITE( 6, '(a)' ) 'We are wrapping over the date line!'
      ENDIF

      ! Make sure that ND59_NI <= IIPAR
      IF ( ND59_NI > IIPAR ) THEN
         CALL ERROR_STOP( 'Too many longitudes!', LOCATION )
      ENDIF

      !-----------
      ! Latitude
      !-----------
      
      ! Error check JMIN_AREA
      IF ( ND59_JMIN+J0 < 1 .or. ND59_JMIN+J0 > JGLOB ) THEN
         CALL ERROR_STOP( 'Bad ND59_JMIN value!', LOCATION)
      ENDIF
     
      ! Error check JMAX_AREA
      IF ( ND59_JMAX+J0 < 1 .or.ND59_JMAX+J0 > JGLOB ) THEN
         CALL ERROR_STOP( 'Bad ND59_JMAX value!', LOCATION)
      ENDIF

      ! Compute latitude limits to write to disk (bey, bmy, 3/16/99)
      IF ( ND59_JMAX >= ND59_JMIN ) THEN      
         ND59_NJ = ( ND59_JMAX - ND59_JMIN ) + 1
      ELSE
         CALL ERROR_STOP( 'ND59_JMAX < ND59_JMIN!', LOCATION )
      ENDIF     
  
      !-----------
      ! Altitude
      !-----------

      ! # of levels to save in ND53 timeseries
      ND59_NL = 1


      !-----------
      ! Offsets
      !-----------
      IOFF      = ND59_IMIN - 1
      JOFF      = ND59_JMIN - 1

      !-----------
      ! For bpch
      !-----------
      TITLE     = 'GEOS-CHEM DIAG59 instantaneous timeseries'
      LONRES    = DISIZE
      LATRES    = DJSIZE
      MODELNAME = GET_MODELNAME()
      HALFPOLAR = GET_HALFPOLAR()
      
      ! Reset grid offsets to global values for bpch write
      I0        = GET_XOFFSET( GLOBAL=.TRUE. )
      J0        = GET_YOFFSET( GLOBAL=.TRUE. )      

      END SUBROUTINE INIT_DIAG59
!EOC
      END MODULE DIAG59_MOD
