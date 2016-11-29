!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: diag63_mod.F
!
! !DESCRIPTION: Module DIAG63\_MOD contains variables and routines to save 
!  out the fraction of NOx remaining and integrated OPE to disk
!  (gvinken, 25/02/11) 
!\\
!\\
! !INTERFACE: 
!
      MODULE DIAG63_MOD
!
! !USES:
!
      USE PRECISION_MOD    ! For GEOS-Chem Precision (fp)

      IMPLICIT NONE
      PRIVATE
!
! !PUBLIC DATA MEMBERS:
!
      LOGICAL, PUBLIC :: DO_SAVE_DIAG63 
!
! !PUBLIC MEMBER FUNCTIONS:
!
      PUBLIC :: DIAG63
      PUBLIC :: ITS_TIME_FOR_DIAG63
      PUBLIC :: INIT_DIAG63
!
! !REMARKS:
!  ND63 tracer numbers:
!  ============================================================================
!  1  : Fraction of NOx remaining  [unitless]
!  2  : Integrated OPE             [molec O3 produced / molec NOx lost]
!
! !REVISION HISTORY:
!  25 Feb 2011 - G. Vinken   - Initial version based on the orig. diag49_mod.f
!  07 Feb 2012 - M. Payer    - Added ProTeX headers
!  24 Feb 2012 - M. Payer    - Rename module from diag59_mod to diag63_mod.
!                              Diag59 is used by TOMAS. Fix this throughout.
!  06 Aug 2012 - R. Yantosca - Now make IU_ND63 a local module variable
!  10 Nov 2014 - M. Yannetti - Added PRECISION_MOD
!  29 Nov 2016 - R. Yantosca - grid_mod.F90 is now gc_grid_mod.F90
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !PRIVATE TYPES:
!
      !=================================================================
      ! MODULE VARIABLES
      !
      ! DO_SAVE_DIAG63   : Switch to turn ND63 timeseries on/off 
      ! I0               : Lon offset between global & nested grid
      ! J0               : Lat offset between global & nested grid
      ! IOFF             : Offset between relative & absolute lon
      ! JOFF             : Offset between relative & absolute lat
      ! ND63_IMIN        : Minimum longitude index
      ! ND63_IMAX        : Maximum latitude  index
      ! ND63_JMIN        : Minimum longitude index
      ! ND63_JMAX        : Maximum longitude index
      ! ND63_FREQ        : Frequency which to save to disk [min]
      ! ND63_N_TRACERS   : Number of tracers for ND63 timeseries
      ! ND63_OUTPUT_FILE : Name of timeseries output file
      ! ND63_TRACERS     : Array w/ tracer #'s to save to disk
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
      INTEGER            :: ND63_N_TRACERS, ND63_TRACERS(120)
      INTEGER            :: ND63_IMIN,      ND63_IMAX
      INTEGER            :: ND63_JMIN,      ND63_JMAX
      INTEGER            :: ND63_FREQ,      ND63_NI
      INTEGER            :: ND63_NJ,        ND63_NL
      INTEGER            :: HALFPOLAR
      INTEGER, PARAMETER :: CENTER180=1 
      REAL*4             :: LONRES,         LATRES
      CHARACTER(LEN=20)  :: MODELNAME
      CHARACTER(LEN=40)  :: RESERVED = ''
      CHARACTER(LEN=80)  :: TITLE
      CHARACTER(LEN=255) :: ND63_OUTPUT_FILE

      ! LUN for ND63 diagnostic file
      INTEGER            :: IU_ND63

      CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: diag63 
!
! !DESCRIPTION: Subroutine DIAG63 produces time series (instantaneous fields) 
!  for a geographical domain from the information read in timeseries.dat.  
!  Output will be in binary punch (BPCH) format.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE DIAG63 ( am_I_Root, Input_Opt, State_Chm, RC )
!
! !USES:
!
      USE BPCH2_MOD,          ONLY : BPCH2,   OPEN_BPCH2_FOR_WRITE
      USE CMN_O3_MOD	      ! Pure O3, SAVENO2
      USE ErrCode_Mod
      USE ERROR_MOD,          ONLY : SAFE_DIV, ERROR_STOP
      USE GC_GRID_MOD,        ONLY : GET_XOFFSET,     GET_YOFFSET
      USE HCO_DIAGN_MOD
      USE HCO_ERROR_MOD
      USE HCO_TYPES_MOD,      ONLY : DiagnCont
      USE HCO_INTERFACE_MOD,  ONLY : HcoState
      USE Input_Opt_Mod,      ONLY : OptInput
      USE inquireMod,         ONLY : findFreeLUN
      USE PhysConstants
      USE PBL_MIX_MOD,        ONLY : GET_PBL_TOP_L,   GET_PBL_TOP_m
      USE State_Chm_Mod,      ONLY : ChmState
      USE State_Chm_Mod,      ONLY : Ind_
      USE TIME_MOD,           ONLY : EXPAND_DATE
      USE TIME_MOD,           ONLY : GET_NYMD,        GET_NHMS
      USE TIME_MOD,           ONLY : GET_NYMD_DIAG,   GET_TS_DIAG
      USE TIME_MOD,           ONLY : GET_TAU,         GET_HOUR
      USE TIME_MOD,           ONLY : ITS_A_NEW_DAY,   TIMESTAMP_STRING
      USE TIME_MOD,           ONLY : GET_TAUb
!
! !INPUT ARGUMENTS
!
      LOGICAL,        INTENT(IN   )       :: am_I_Root
      TYPE(OptInput), INTENT(IN   )       :: Input_Opt  ! Input opts 
!
! !INPUT/OUTPUT ARGUMENTS
!
      TYPE(ChmState), INTENT(INOUT)       :: State_Chm  ! Chemistry State object

      INTEGER,        INTENT(INOUT)       :: RC
!
! !REVISION HISTORY:
!  25 Feb 2011 - G. Vinken   - Initial version based on DIAG49 
!  07 Feb 2012 - M. Payer    - Added ProTeX headers
!  11 Apr 2012 - R. Yantosca - Remove reference to lai_mod.F, it's not needed
!  03 Aug 2012 - R. Yantosca - Move calls to findFreeLUN out of DEVEL block
!  17 Dec 2014 - R. Yantosca - Leave time/date variables as 8-byte
!  26 Feb 2015 - E. Lundgren - Remove pressure_mod from use list (not used)
!  29 Apr 2016 - R. Yantosca - Don't initialize pointers in declaration stmts
!  31 May 2016 - E. Lundgren - Use TCVV instead of XNUMOL for molecular wts
!  22 Jun 2016 - M. Yannetti - Replaced references to TCVV, now takes in
!                              State_Chm for spec db access
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER                       :: IOS,  GMTRC, GMNL, I, J, K, L 
      INTEGER                       :: N,    R,     H,    W, X, Y
      INTEGER                       :: NHMS, TS_DIAG
      REAL(f8)                      :: TAU1, TAU0
      REAL(fp)                      :: Q( ND63_NI, ND63_NJ, 1)
      REAL(fp)                      :: MW_NO
      CHARACTER(LEN=16)             :: STAMP
      CHARACTER(LEN=40)             :: CATEGORY
      CHARACTER(LEN=40)             :: UNIT
      CHARACTER(LEN=255)            :: FILENAME

      ! Connect to HEMCO diagnostics
      INTEGER                       :: FLAG
      INTEGER                       :: AF
      CHARACTER(LEN= 63)            :: DiagnName
      CHARACTER(LEN=255)            :: MSG
      TYPE(DiagnCont), POINTER      :: DiagnCnt
      REAL(f8)                      :: FACTOR

      ! Species ID flag
      INTEGER                       :: id_NO
!
! !DEFINED PARAMETERS:
!
      REAL(f8), PARAMETER           :: CM2PERM2 = 10000e+0_f8
      CHARACTER(LEN=255), PARAMETER :: LOC = 'DIAG63 (DIAG63_MOD.F)' 

      !=================================================================
      ! DIAG63 begins here!
      !=================================================================

      ! Initialize
      DiagnCnt => NULL()

      !=================================================================
      ! If it's a new day, open a new BPCH file and write file header
      ! We need to check if it's a new day + 1 ND63 time step (ccc, 8/12/09)
      !=================================================================
!--- Previous to (ccc, 8/12/09)
!      IF ( ITS_A_NEW_DAY() ) THEN
      NHMS    = GET_NHMS()
      TS_DIAG = ND63_FREQ

      ! To change TS_DIAG to NHMS format
      TS_DIAG = TS_DIAG/60 * 10000 + (TS_DIAG - (TS_DIAG/60)*60) * 100  

      ! Check if it's a new day for diagnostics, 
      ! or if the model run just started
      IF ( (NHMS == TS_DIAG) .OR. 
     &     (GET_TAU() == ( GET_TAUb() + ND63_FREQ/60e+0_fp )) ) THEN

         ! Find a free file LUN
         IU_ND63  = findFreeLUN()

         ! Expand date tokens in the file name
         FILENAME = TRIM( ND63_OUTPUT_FILE )
         CALL EXPAND_DATE( FILENAME, GET_NYMD(), GET_NHMS() )

         ! Echo info
         WRITE( 6, 100 ) TRIM( FILENAME )
 100     FORMAT( '     - DIAG63: Opening file ', a, ' on unit ', i6 )
        
         ! Open bpch file and write top-of-file header
         CALL OPEN_BPCH2_FOR_WRITE( IU_ND63, FILENAME, TITLE )
      ENDIF

      !=================================================================
      ! Save tracers to timeseries file
      !=================================================================

      ! Get the species ID of NO
      id_NO = Ind_('NO')
      IF ( id_NO < 0 ) THEN
         MSG = 'PARANOX is active but NO is not defined!'
         CALL ERROR_STOP ( MSG, LOC )
      ENDIF

      ! Get MW from species database
      MW_NO = State_Chm%SpcData(id_NO)%Info%emMW_g

      ! Echo info
      STAMP = TIMESTAMP_STRING()
      WRITE( 6, 110 ) STAMP
 110  FORMAT( '     - DIAG63: Saving timeseries at ', a )

      ! Time for BPCH file
      TAU1  = GET_TAU()
      
      ! Time for start of averaging period
      TAU0  = TAU1 - (ND63_FREQ / 60e+0_fp) 

      ! Loop over tracers
      DO W = 1, ND63_N_TRACERS

         ! Zero summing array
         Q = 0e+0_fp

         !-------------------------------------
         ! SHIP-diagnostic tracers [unitless]
         !-------------------------------------
         CATEGORY = 'SHIP-$$$'
         GMNL     = ND63_NL
         GMTRC    = W

         ! Pick variables to get diagnostics from HEMCO
         SELECT CASE ( W )
            CASE ( 1 ) 
               DiagnName = 'PARANOX_NOXFRAC_REMAINING'
               UNIT      = 'UNITLESS'
               AF        = 0
               FACTOR    = 1.0d0
            CASE ( 2 ) 
               DiagnName = 'PARANOX_OPE'
               UNIT      = 'UNITLESS'
               AF        = 0
               FACTOR    = 1.0d0
            CASE ( 3 ) 
               DiagnName = 'PARANOX_OX_PRODUCTION'
               UNIT      = 'molec/cm2/s'
               AF        = 0
               FACTOR    = AVO * ( AIRMW / MW_NO ) / 
     &                     ( CM2PERM2 * AIRMW * 1.e-3_fp ) 
            CASE ( 4 ) 
               DiagnName = 'PARANOX_NOX_REMAINING'
               UNIT      = 'molec/cm2/s'
               AF        = 0
               FACTOR    = AVO * ( AIRMW / MW_NO ) / 
     &                     ( CM2PERM2 * AIRMW * 1.e-3_fp ) 
            CASE ( 5 ) 
               DiagnName = 'SHIP_NO'
               UNIT      = 'molec/cm2/s'
               AF        = 1
               FACTOR    = AVO * ( AIRMW / MW_NO ) /
     &                     ( CM2PERM2 * AIRMW * 1.e-3_fp )  
            CASE DEFAULT
               MSG = 'Only 5 diagnostics expected!'
               CALL ERROR_STOP ( MSG, LOC )
         END SELECT

         ! Get diagnostics from HEMCO
         DiagnCnt => NULL()
         CALL Diagn_Get( am_I_Root, HcoState, .FALSE., DiagnCnt,
     &                   FLAG,  RC, cName=TRIM(DiagnName),
     &                   AutoFill=AF, 
     &                   COL=HcoState%Diagn%HcoDiagnIDManual )
         IF ( RC /= HCO_SUCCESS ) THEN
            MSG = 'Cannot find diagnostics ' // TRIM(DiagnName)
            CALL ERROR_STOP ( MSG, LOC )
         ENDIF
        
         ! Save into ARRAY. Convert to desired output units
         IF ( FLAG == HCO_SUCCESS ) THEN
            Q(:,:,1) = DiagnCnt%Arr2D%Val(:,:) * FACTOR
         ELSE
            MSG = 'No diagnostics returned: ' // TRIM(DiagnName)
            MSG = TRIM(MSG) // ' - will write zeros!'
            CALL HCO_WARNING ( MSG, RC, THISLOC=LOC )
         ENDIF
 
         !==============================================================
         ! Save this data block to the ND63 timeseries file
         !==============================================================
         CALL BPCH2( IU_ND63,      MODELNAME,    LONRES,   
     &               LATRES,       HALFPOLAR,    CENTER180, 
     &               CATEGORY,     GMTRC,        UNIT,      
     &               TAU0,         TAU1,         RESERVED,  
     &               ND63_NI,      ND63_NJ,      1,  
     &               ND63_IMIN+I0, ND63_JMIN+J0, 1, 
     &               REAL( Q(1:ND63_NI, 1:ND63_NJ, 1) ) )
      ENDDO
 
      !=================================================================
      ! Close the file at the proper time
      !=================================================================
      IF ( ITS_TIME_TO_CLOSE_FILE() ) THEN

         ! Expand date tokens in the file name
         FILENAME = TRIM( ND63_OUTPUT_FILE )
!--- Previous to (ccc, 8/12/09)
!         CALL EXPAND_DATE( FILENAME, GET_NYMD(), GET_NHMS() )
         CALL EXPAND_DATE( FILENAME, GET_NYMD_DIAG(), GET_NHMS() )

         ! Echo info
         WRITE( 6, 120 ) TRIM( FILENAME )
 120     FORMAT( '     - DIAG63: Closing file : ', a )

         ! Close file
         CLOSE( IU_ND63 ) 
      ENDIF

      ! Return w/ success
      DiagnCnt => NULL()
      RC = GC_SUCCESS

      END SUBROUTINE DIAG63
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: its_time_to_close_file
!
! !DESCRIPTION: Function ITS\_TIME\_TO\_CLOSE\_FILE returns TRUE if it's 
!  time to close the ND63 bpch file before the end of the day.
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
      REAL(fp)  :: HR1

      !=================================================================
      ! ITS_TIME_TO_CLOSE_FILE begins here!
      !=================================================================

      ! Current hour
      HR1      = GET_HOUR() + ( GET_MINUTE() / 60e+0_fp )

!--- Previous to (ccc, 8/12/09)
!      ! Hour at the next dynamic timestep
!      HR2      = HR1        + ( ND63_FREQ / 60d0 )

       ! If the next dyn step is the start of a new day, return TRUE
!--- Previous to (ccc, 11/11/10)
!       HR1 varies between 00 and 23:59. So compares to 00 not 24 anymore.
!      ITS_TIME = ( INT( HR1 ) == 24 )
      ITS_TIME = ( INT( HR1 ) == 00 )
 
       END FUNCTION ITS_TIME_TO_CLOSE_FILE
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: its_time_for_diag63
!
! !DESCRIPTION: Function ITS\_TIME\_FOR\_DIAG63 returns TRUE if ND63 is 
!  turned on and it is time to call DIAG63 -- or FALSE otherwise.
!\\
!\\
! !INTERFACE:
!
      FUNCTION ITS_TIME_FOR_DIAG63() RESULT( ITS_TIME )
!
! !USES:
!
      USE ERROR_MOD, ONLY : GEOS_CHEM_STOP
      USE TIME_MOD,  ONLY : GET_ELAPSED_MIN
      USE TIME_MOD,  ONLY : GET_TS_DIAG
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
      ! ITS_TIME_FOR_DIAG63 begins here!
      !=================================================================

      IF ( DO_SAVE_DIAG63 ) THEN
         IF ( FIRST ) THEN
            TS_DIAG = GET_TS_DIAG()
            
            ! Check if ND63_FREQ is a multiple of TS_DIAG
            IF ( MOD( ND63_FREQ, TS_DIAG ) /= 0 ) THEN
               WRITE( 6, 100 ) 'ND63', ND63_FREQ, TS_DIAG
 100           FORMAT( 'The ',a,' output frequency must be a multiple '
     &              'of the largest time step:', i5, i5 )
               CALL GEOS_CHEM_STOP
            ENDIF
            FIRST = .FALSE.
         ENDIF
         
         ! Time already elapsed in this run
         XMIN     = GET_ELAPSED_MIN()
         
         ! Is the elapsed time a multiple of ND63_FREQ?
         ITS_TIME = ( DO_SAVE_DIAG63 .and. MOD( XMIN, ND63_FREQ ) == 0 )
      ELSE
         ITS_TIME = DO_SAVE_DIAG63
      ENDIF
            
      END FUNCTION ITS_TIME_FOR_DIAG63
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
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
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_diag63
!
! !DESCRIPTION: Subroutine INIT\_DIAG63 allocates and zeroes all module 
!  arrays.  It also gets values for module variables from "input\_mod.f". 
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE INIT_DIAG63( DO_ND63, N_ND63, TRACERS, IMIN,    
     &                        IMAX,    JMIN,   JMAX,    FREQ,   
     &                        FILE )
!
! !USES:
!
      USE BPCH2_MOD,   ONLY : GET_MODELNAME
      USE BPCH2_MOD,   ONLY : GET_HALFPOLAR
      USE CMN_SIZE_MOD      ! Size parameters
      USE GC_GRID_MOD, ONLY : GET_XOFFSET
      USE GC_GRID_MOD, ONLY : GET_YOFFSET
      USE GC_GRID_MOD, ONLY : ITS_A_NESTED_GRID
      USE ERROR_MOD,   ONLY : ERROR_STOP
!
! !INPUT PARAMETERS: 
!
      ! DO_ND63 : Switch to turn on ND63 timeseries diagnostic
      ! N_ND63  : Number of ND63 read by "input_mod.f"
      ! TRACERS : Array w/ ND63 tracer #'s read by "input_mod.f"
      ! IMIN    : Min longitude index read by "input_mod.f"
      ! IMAX    : Max longitude index read by "input_mod.f" 
      ! JMIN    : Min latitude index read by "input_mod.f" 
      ! JMAX    : Min latitude index read by "input_mod.f" 
      ! FREQ    : Frequency for saving to disk [min]
      ! FILE    : ND63 output file name read by "input_mod.f"
      LOGICAL,            INTENT(IN) :: DO_ND63
      INTEGER,            INTENT(IN) :: N_ND63, TRACERS(100)
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
      ! INIT_DIAG63 begins here!
      !=================================================================

      ! Initialize
      LOCATION               = 'INIT_DIAG63 ("diag63_mod.f")'
      ND63_TRACERS(:)        = 0

      ! Get values from "input_mod.f"
      DO_SAVE_DIAG63         = DO_ND63 
      ND63_N_TRACERS         = N_ND63
      ND63_TRACERS(1:N_ND63) = TRACERS(1:N_ND63)
      ND63_IMIN              = IMIN
      ND63_IMAX              = IMAX
      ND63_JMIN              = JMIN
      ND63_JMAX              = JMAX
      ND63_FREQ              = FREQ
      ND63_OUTPUT_FILE       = FILE
     
      ! Return if we are not saving ND63 diagnostics
      IF ( .not. DO_SAVE_DIAG63 ) RETURN

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

      ! Error check ND63_IMIN
      IF ( ND63_IMIN+I0 < 1 .or. ND63_IMIN+I0 > IGLOB ) THEN
         CALL ERROR_STOP( 'Bad ND63_IMIN value!', LOCATION )
      ENDIF

      ! Error check ND63_IMAX
      IF ( ND63_IMAX+I0 < 1 .or. ND63_IMAX+I0 > IGLOB ) THEN
         CALL ERROR_STOP( 'Bad ND63_IMAX value!', LOCATION )
      ENDIF

      ! Compute longitude limits to write to disk 
      ! Also handle wrapping around the date line
      IF ( ND63_IMAX >= ND63_IMIN ) THEN
         ND63_NI = ( ND63_IMAX - ND63_IMIN ) + 1
      ELSE 
         ND63_NI = ( IIPAR - ND63_IMIN ) + 1 + ND63_IMAX
         WRITE( 6, '(a)' ) 'We are wrapping over the date line!'
      ENDIF

      ! Make sure that ND63_NI <= IIPAR
      IF ( ND63_NI > IIPAR ) THEN
         CALL ERROR_STOP( 'Too many longitudes!', LOCATION )
      ENDIF

      !-----------
      ! Latitude
      !-----------
      
      ! Error check JMIN_AREA
      IF ( ND63_JMIN+J0 < 1 .or. ND63_JMIN+J0 > JGLOB ) THEN
         CALL ERROR_STOP( 'Bad ND63_JMIN value!', LOCATION)
      ENDIF
     
      ! Error check JMAX_AREA
      IF ( ND63_JMAX+J0 < 1 .or.ND63_JMAX+J0 > JGLOB ) THEN
         CALL ERROR_STOP( 'Bad ND63_JMAX value!', LOCATION)
      ENDIF

      ! Compute latitude limits to write to disk (bey, bmy, 3/16/99)
      IF ( ND63_JMAX >= ND63_JMIN ) THEN      
         ND63_NJ = ( ND63_JMAX - ND63_JMIN ) + 1
      ELSE
         CALL ERROR_STOP( 'ND63_JMAX < ND63_JMIN!', LOCATION )
      ENDIF     
  
      !-----------
      ! Altitude
      !-----------

      ! # of levels to save in ND53 timeseries
      ND63_NL = 1


      !-----------
      ! Offsets
      !-----------
      IOFF      = ND63_IMIN - 1
      JOFF      = ND63_JMIN - 1

      !-----------
      ! For bpch
      !-----------
      TITLE     = 'GEOS-CHEM DIAG63 instantaneous timeseries'
      LONRES    = DISIZE
      LATRES    = DJSIZE
      MODELNAME = GET_MODELNAME()
      HALFPOLAR = GET_HALFPOLAR()
      
      ! Reset grid offsets to global values for bpch write
      I0        = GET_XOFFSET( GLOBAL=.TRUE. )
      J0        = GET_YOFFSET( GLOBAL=.TRUE. )      

      END SUBROUTINE INIT_DIAG63
!EOC
      END MODULE DIAG63_MOD
