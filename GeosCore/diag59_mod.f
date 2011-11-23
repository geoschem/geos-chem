! $Id: diag59_mod.f,v 1.4 2011/02/25 gvinken Exp $
      MODULE DIAG59_MOD
!
!******************************************************************************
!  Module DIAG59_MOD contains variables and routines to save out the fraction 
!  of NOx remaining and integrated OPE to disk (gvinken, 25/02/11)
!      Based on the original diag49_mod.f (bmy, 7/20/04, 10/13/09)
!
!  Module Variables:
!  ============================================================================
!  (1 ) DO_SAVE_DIAG59   (LOGICAL ) : Switch to turn ND59 timeseries on/off 
!  (2 ) I0               (INTEGER ) : Lon offset between global & nested grid
!  (3 ) J0               (INTEGER ) : Lat offset between global & nested grid
!  (4 ) IOFF             (INTEGER ) : Offset between relative & absolute lon
!  (5 ) JOFF             (INTEGER ) : Offset between relative & absolute lat
!  (6 ) ND59_IMIN        (INTEGER ) : Minimum longitude index
!  (7 ) ND59_IMAX        (INTEGER ) : Maximum latitude  index
!  (8 ) ND59_JMIN        (INTEGER ) : Minimum longitude index
!  (9 ) ND59_JMAX        (INTEGER ) : Maximum longitude index
!  (10) ND59_FREQ        (INTEGER ) : Frequency which to save to disk [min]
!  (11) ND59_N_TRACERS   (INTEGER ) : Number of tracers for ND49 timeseries
!  (12) ND59_OUTPUT_FILE (CHAR*255) : Name of timeseries output file
!  (13) ND59_TRACERS     (INTEGER ) : Array w/ tracer #'s to save to disk
!  (14) HALFPOLAR        (INTEGER ) : Used for binary punch file write
!  (15) CENTER180        (INTEGER ) : Used for binary punch file write
!  (16) LONRES           (REAL*4  ) : Used for binary punch file write
!  (17) LATRES           (REAL*4  ) : Used for binary punch file write
!  (18) RESERVED         (CHAR*40 ) : Used for binary punch file write
!  (19) MODELNAME        (CHAR*20 ) : Used for binary punch file write
!  (20) TITLE            (CHAR*80 ) : Used for binary punch file write 
!
!  Module Routines:
!  ============================================================================
!  (1 ) DIAG59                 : Main driver routine
!  (2 ) ITS_TIME_TO_CLOSE_FILE : Returns TRUE if it's time to close ND49 file
!  (3 ) ITS_TIME_FOR_DIAG59    : Returns TRUE if it's time to save to disk
!  (4 ) GET_I                  : Converts relative longitude index to absolute
!  (5 ) INIT_DIAG59            : Gets variable values from "input_mod.f"
!
!  GEOS-CHEM modules referenced by "diag59_mod.f" 
!  ============================================================================
!  (1 ) bpch2_mod.f            : Module w/ routines for binary punch file I/O
!  (2 ) dao_mod.f              : Module w/ arrays for DAO met fields
!  (3 ) file_mod.f             : Module w/ file unit numbers & error checks
!  (4 ) grid_mod.f             : Module w/ horizontal grid information   
!  (5 ) pbl_mix_mod.f          : Module w/ routines for PBL height & mixing
!  (6 ) pressure_mod.f         : Module w/ routines to compute P(I,J,L)
!  (7 ) time_mod.f             : Module w/ routines for computing time & date 
!  (8 ) tracer_mod.f           : Module w/ GEOS-CHEM tracer array STT etc.  
!  (9 ) tracerid_mod.f         : Module w/ pointers to tracers & emissions
!
!  ND59 tracer numbers:
!  ============================================================================
!  1             : Fraction of NOx remaining [unitless]
!  2             : Integrated OPE (molecules of O3 produced per NOx molecule lost)
! 
!  NOTES:
!  (1 ) 
!******************************************************************************
!
      IMPLICIT NONE

      !=================================================================
      ! MODULE PRIVATE DECLARATIONS -- keep certain internal variables 
      ! and routines from being seen outside "diag59_mod.f"
      !=================================================================

      ! Make everything PRIVATE ...
      PRIVATE 

      ! ... except these variables ...
      PUBLIC :: DO_SAVE_DIAG59 

      ! ... except these routines 
      PUBLIC :: DIAG59
      PUBLIC :: ITS_TIME_FOR_DIAG59
      PUBLIC :: INIT_DIAG59

      !=================================================================
      ! MODULE VARIABLES
      !=================================================================
      LOGICAL            :: DO_SAVE_DIAG59
      INTEGER            :: IOFF,           JOFF
      INTEGER            :: I0,             J0
      ! Increased to 120 from 100 (mpb,2009)
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

      !=================================================================
      ! MODULE ROUTINES -- follow below the "CONTAINS" statement
      !=================================================================
      CONTAINS

!------------------------------------------------------------------------------

      SUBROUTINE DIAG59
! 
!******************************************************************************
!  Subroutine DIAG59 (based on DIAG49 from diag49_mod.f) produces time series
!  (instantaneous fields) for a geographical domain from the information read 
!  in timeseries.dat.  Output will be in binary punch (BPCH) format. 
!  (bey, bmy, rvm, 4/9/99, 10/7/08)
!
!  NOTES:
!  (1 ) Now bundled into "diag49_mod.f".  Now reference STT from 
!        "tracer_mod.f".  Now scale aerosol & dust OD's to 400 nm.  
!        (bmy, rvm, aad, 7/9/04)
!  (2 ) Updated tracer # for NO2 (bmy, 10/25/04)
!  (3 ) Remove reference to "CMN".  Also now get PBL heights in meters and 
!        model layers from GET_PBL_TOP_m and GET_PBL_TOP_L of "pbl_mix_mod.f".
!        (bmy, 2/16/05)
!  (4 ) Now reference CLDF and BXHEIGHT from "dao_mod.f".  Now save 3-D cloud 
!        fraction as tracer #79 and box height as tracer #93.  Now remove 
!        reference to PBL from "dao_mod.f"(bmy, 4/20/05)
!  (5 ) Remove references to TRCOFFSET because it is always zero (bmy, 6/24/05)
!  (6 ) Now do not save SLP data if it is not allocated (bmy, 8/2/05)
!  (7 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (8 ) Now references XNUMOLAIR from "tracer_mod.f".  Bug fix: now must sum
!        aerosol OD's over all RH bins.  Also zero Q array. (bmy, 11/1/05)
!  (9 ) Bug fix: accumulate into Q(X,Y,K) for dust OD (qli, bmy, 4/30/07)
!  (10) Bug fix: UNIT should be "levels" for tracer 77.  Also RH should be
!        tracer #17 under "TIME-SER" category. (cdh, bmy, 2/11/08)
!  (11) Bug fix: replace "PS-PTOP" with "PEDGE-$" (bmy, phs, 10/7/08)
!  (12) Change the new day condition to open a new file. (ccc, 8/12/09)
!  (13) Change the timestamp for the filename when closing (ccc, 8/12/09)
!  (14) Add outputs for EMISS_BVOC (10 tracers), TS, PARDR, PARDF and ISOLAI
!        (mpb, 11/19/09)
!******************************************************************************
!
      ! References to F90 modules
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

#     include "cmn_fj.h"        ! FAST-J stuff, includes CMN_SIZE
#     include "jv_cmn.h"        ! ODAER, QAA, QAA_AOD (clh)
#     include "CMN_O3"		! Pure O3, SAVENO2
#     include "CMN_GCTM"        ! XTRA2

      ! Local variables
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
      ! We need to check if it's a new day + 1 ND49 time step (ccc, 8/12/09)
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

      ! Return to calling program
      END SUBROUTINE DIAG59

!------------------------------------------------------------------------------

      FUNCTION ITS_TIME_TO_CLOSE_FILE() RESULT( ITS_TIME )
!
!******************************************************************************
!  Function ITS_TIME_TO_CLOSE_FILE returns TRUE if it's time to close the
!  ND59 bpch file before the end of the day. (bmy, 7/20/04)
!
!  NOTES:
!  (1 ) The time is already updated to the next time step (ccc, 8/12/09)
!  (2 ) The time varies now between 00 and 23:59, so we need to compare HR1 to
!       00 instead of 24. (ccc, 11/11/10)
!******************************************************************************
!
      ! References to F90 modules
      USE TIME_MOD, ONLY : GET_HOUR, GET_MINUTE

      ! Local variables
      LOGICAL :: ITS_TIME
      REAL*8  :: HR1

      !=================================================================
      ! ITS_TIME_TO_CLOSE_FILE begins here!
      !=================================================================

      ! Current hour
      HR1      = GET_HOUR() + ( GET_MINUTE() / 60d0 )

!--- Previous to (ccc, 8/12/09)
!      ! Hour at the next dynamic timestep
!      HR2      = HR1        + ( ND49_FREQ / 60d0 )

       ! If the next dyn step is the start of a new day, return TRUE
!--- Previous to (ccc, 11/11/10)
!       HR1 varies between 00 and 23:59. So compares to 00 not 24 anymore.
!      ITS_TIME = ( INT( HR1 ) == 24 )
      ITS_TIME = ( INT( HR1 ) == 00 )
 
       ! Return to calling program
       END FUNCTION ITS_TIME_TO_CLOSE_FILE

!------------------------------------------------------------------------------

      FUNCTION ITS_TIME_FOR_DIAG59() RESULT( ITS_TIME )
!
!******************************************************************************
!  Function ITS_TIME_FOR_DIAG59 returns TRUE if ND59 is turned on and it is 
!  time to call DIAG59 -- or FALSE otherwise. (bmy, 7/20/04)
!
!  NOTES:
!  (1 ) Add a check on the output frequency for validity compared to time 
!        steps used. (ccc, 5/21/09)
!******************************************************************************
!
      ! References to F90 modules
      USE TIME_MOD,  ONLY : GET_ELAPSED_MIN, GET_TS_DIAG
      USE ERROR_MOD, ONLY : GEOS_CHEM_STOP

      ! Local variables
      INTEGER :: XMIN, TS_DIAG
      LOGICAL :: ITS_TIME
      LOGICAL, SAVE :: FIRST = .TRUE.

      !=================================================================
      ! ITS_TIME_FOR_DIAG59 begins here!
      !=================================================================

      IF ( DO_SAVE_DIAG59 ) THEN
         IF ( FIRST ) THEN
            TS_DIAG = GET_TS_DIAG()
            
            ! Check if ND59_FREQ is a multiple of TS_DIAG
            IF ( MOD( ND59_FREQ, TS_DIAG ) /= 0 ) THEN
               WRITE( 6, 100 ) 'ND49', ND59_FREQ, TS_DIAG
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
            
      ! Return to calling program
      END FUNCTION ITS_TIME_FOR_DIAG59

!------------------------------------------------------------------------------

      FUNCTION GET_I( X ) RESULT( I )
!
!******************************************************************************
!  Function GET_I returns the absolute longitude index (I), given the 
!  relative longitude index (X).  (bmy, 7/20/04)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) X (INTEGER) : Relative longitude index (used by Q)
!
!  NOTES:
!******************************************************************************
!      
#     include "CMN_SIZE"   ! Size parameters

      ! Arguments
      INTEGER, INTENT(IN) :: X

      ! Local variables
      INTEGER             :: I
      
      !=================================================================
      ! GET_I begins here!
      !=================================================================

      ! Add the offset to X to get I  
      I = IOFF + X

      ! Handle wrapping around the date line, if necessary
      IF ( I > IIPAR ) I = I - IIPAR

      ! Return to calling program
      END FUNCTION GET_I

!------------------------------------------------------------------------------

      SUBROUTINE INIT_DIAG59( DO_ND59, N_ND59, TRACERS, IMIN,    
     &                        IMAX,    JMIN,   JMAX, FREQ,   FILE )
!
!******************************************************************************
!  Subroutine INIT_DIAG59 allocates and zeroes all module arrays.  
!  It also gets values for module variables from "input_mod.f". Based on 
!  INIT_DIAG59 (bmy, 7/20/04, 11/30/06)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) DO_ND59 (LOGICAL ) : Switch to turn on ND59 timeseries diagnostic
!  (2 ) N_ND59  (INTEGER ) : Number of ND59 read by "input_mod.f"
!  (3 ) TRACERS (INTEGER ) : Array w/ ND59 tracer #'s read by "input_mod.f"
!  (4 ) IMIN    (INTEGER ) : Min longitude index read by "input_mod.f"
!  (5 ) IMAX    (INTEGER ) : Max longitude index read by "input_mod.f" 
!  (6 ) JMIN    (INTEGER ) : Min latitude index read by "input_mod.f" 
!  (7 ) JMAX    (INTEGER ) : Min latitude index read by "input_mod.f" 
!  (8 ) LMIN    (INTEGER ) : Min level index read by "input_mod.f" 
!  (9 ) LMAX    (INTEGER ) : Min level index read by "input_mod.f" 
!  (10) FREQ    (INTEGER ) : Frequency for saving to disk [min]
!  (11) FILE    (CHAR*255) : ND59 output file name read by "input_mod.f"
! 
!  NOTES:
!  (1 ) Now get I0 and J0 correctly for nested grid simulations (bmy, 11/9/04)
!  (2 ) Now call GET_HALFPOLAR from "bpch2_mod.f" to get the HALFPOLAR flag 
!        value for GEOS or GCAP grids. (bmy, 6/28/05)
!  (3 ) Now allow ND49_IMIN to be equal to ND49_IMAX and ND49_JMIN to be
!        equal to ND49_JMAX.  This will allow us to save out longitude
!        or latitude transects.  (cdh, bmy, 11/30/06)
!******************************************************************************
!      
      ! References to F90 modules
      USE BPCH2_MOD, ONLY : GET_MODELNAME, GET_HALFPOLAR
      USE GRID_MOD,  ONLY : GET_XOFFSET, GET_YOFFSET, ITS_A_NESTED_GRID
      USE ERROR_MOD, ONLY : ERROR_STOP

#     include "CMN_SIZE" ! Size parameters

      ! Arguments
      LOGICAL,            INTENT(IN) :: DO_ND59
      INTEGER,            INTENT(IN) :: N_ND59, TRACERS(100)
      INTEGER,            INTENT(IN) :: IMIN,   IMAX 
      INTEGER,            INTENT(IN) :: JMIN,   JMAX      
      INTEGER,            INTENT(IN) :: FREQ
      CHARACTER(LEN=255), INTENT(IN) :: FILE

      ! Local variables
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

      ! Error check ND49_IMAX
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

      ! Return to calling program
      END SUBROUTINE INIT_DIAG59

!------------------------------------------------------------------------------

      ! End of module
      END MODULE DIAG59_MOD
