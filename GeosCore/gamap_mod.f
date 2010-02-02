! $Id: gamap_mod.f,v 1.7 2010/02/02 16:57:53 bmy Exp $
      MODULE GAMAP_MOD
!
!******************************************************************************
!  Module GAMAP_MOD contains routines to create GAMAP "tracerinfo.dat" and
!  "diaginfo.dat" files which are customized to each particular GEOS-Chem
!  simulation. (bmy, 5/3/05, 12/18/09)
! 
!  Module Variables:
!  ============================================================================
!  (1 ) MAXCAT    (INTEGER ) : Maximum # of GAMAP diagnostic categories
!  (2 ) SPACING   (INTEGER ) : Spacing between GAMAP diagnostic categories
!  (3 ) NCATS     (INTEGER ) : Counter for number of defined GAMAP categories
!  (4 ) OFFSET    (INTEGER ) : GAMAP offset for each diagnostic category
!  (5 ) CATEGORY  (CHAR*40 ) : GAMAP category name
!  (6 ) DESCRIPT  (CHAR*40 ) : GAMAP category description
!  (7 ) DFILE     (CHAR*255) : Path name of the GAMAP "diaginfo.dat" file
!  (8 ) MAXDIAG   (INTEGER ) : Maximum # of GEOS-CHEM diagnostics
!  (9 ) MAXTRACER (INTEGER ) : Maximum # of tracers per GEOS-CHEM diagnostic
!  (10) NTRAC     (INTEGER ) : Number of tracers for each GEOS-CHEM diagnostic
!  (11) INDEX     (INTEGER ) : Diagnostic tracer numbers
!  (12) MOLC      (INTEGER ) : Ratio of (moles C) / (moles tracer)
!  (13) MWT       (REAL*4  ) : Tracer molecular weights
!  (14) SCALE     (REAL*4  ) : GAMAP scale factors (e.g. 1e9 for v/v -> ppbv)
!  (15) NAME      (CHAR*40 ) : Tracer names
!  (16) FNAME     (CHAR*40 ) : Full tracer names
!  (17) UNIT      (CHAR*40 ) : Unit string for each tracer
!  (18) TFILE     (CHAR*255) : Path name of the GAMAP "tracerinfo.dat" file
!  (19) STAMP     (CHAR*16 ) : Timestamp w/ system date & time
!  (20) SIM_NAME  (CHAR*40 ) : Name of the GEOS-CHEM simulation 
!
!  Module Routines:
!  ============================================================================
!  (1 ) DO_GAMAP             : Driver routine
!  (2 ) CREATE_DINFO         : Writes customized "diaginfo.dat" file
!  (3 ) CREATE_TINFO         : Writes customized "tracerinfo.dat" file
!  (4 ) WRITE_TINFO          : Writes one line to disk for "tracerinfo.dat" 
!  (5 ) WRITE_SEPARATOR      : Writes separator blocks to "tracerinfo.dat"
!  (6 ) INIT_DIAGINFO        : Initializes arrays for diaginfo.dat file
!  (7 ) INIT_TRACERINFO      : Initializes arrays for tracerinfo.dat file
!  (8 ) INIT_GAMAP           : Allocates and initializes all module arrays
!  (9 ) CLEANUP_GAMAP        : Deallocates all module arrays
!
!  GEOS-CHEM modules referenced by "gamap_mod.f"
!  ============================================================================
!  (1 ) diag03_mod.f         : Module w/ routines for mercury diagnostics
!  (2 ) diag41_mod.f         : Module w/ routines for afternoon PBL diag's
!  (3 ) diag48_mod.f         : Module w/ routines for station timeseries
!  (4 ) diag49_mod.f         : Module w/ routines for inst timeseries
!  (5 ) diag50_mod.f         : Module w/ routines for 24hr avg timeseries
!  (6 ) diag51_mod.f         : Module w/ routines for morning/aft timeseries
!  (7 ) diag_pl_mod.f        : Module w/ routines for saving family P & L
!  (8 ) drydep_mod.f         : Module w/ GEOS-CHEM drydep routines
!  (9 ) error_mod.f          : Module w/ error and NaN check routines
!  (10) file_mod.f           : Module w/ file unit numbers & I/O error checks
!  (11) time_mod.f           : Module w/ routines for computing time & date
!  (12) tracer_mod.f         : Module w/ GEOS-CHEM tracer array STT etc.
!  (13) tracerid_mod.f       : Module w/ GEOS-CHEM tracer ID flags
!  (14) wetscav_mod.f        : Module w/ routines for wetdep/scavenging
!
!  References:
!  ============================================================================
!  (1 ) For more information, please see the GAMAP Online Users' Manual:
!        http://www-as.harvard.edu/chemistry/trop/gamap/documentation/
!
!  NOTES:
!  (1 ) Minor bug fix for Rn/Pb/Be simulations (bmy, 5/11/05)
!  (2 ) Added ND09 diagnostic for HCN/CH3CN simulation. (bmy, 6/30/05)
!  (3 ) Added ND04 diagnostic for CO2 simulation (bmy, 7/25/05)
!  (4 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (5 ) Add MBO to ND46 diagnostic (tmf, bmy, 10/20/05)
!  (6 ) Updated for tagged Hg simulation (cdh, bmy, 4/6/06)
!  (7 ) Updated for ND56 lightning flash diagnostics (ltm, bmy, 5/5/06)
!  (8 ) Updated for ND42 SOA concentration diagnostics (dkh, bmy, 5/22/06)
!  (9 ) Updated for ND36 CH3I simulation diagnostics (bmy, 7/25/06)
!  (10) Remove support for GEOS-1 and GEOS-STRAT met fields (bmy, 8/4/06)
!  (11) Add routines INIT_DIAGINFO, INIT_TRACERINFO for clarity.  Added new
!        entries for biomass burning (ND28) and time in tropopshere (ND54)
!        in INIT_DIAGINFO and INIT_TRACERINFO. (phs, bmy, 10/17/06)
!  (12) Now write GPROD & APROD info to diaginfo.dat, tracerinfo.dat files,
!        for the SOA restart files (tmf, havala, bmy, 2/6/07)
!  (13) Added ND10 diagnostic for H2/HD simulation. (phs, 9/18/07)
!  (14) Change category name for ND31 diagnostic (bmy, 11/16/07)
!  (15) Add to tracerinfo.dat file for timeseries and Rn-Pb-Be (bmy, 2/22/08)
!  (16) Added ND52 diagnostic for gamma HO2 (jaegle 02/26/09)
!  (17) Add gamap info for dicarbonyl simulation (tmf, 3/10/09)
!  (18) Add C2H4 in ND46 (ccc, 3/10/09)
!  (19) Add EFLUX to ND67 (lin, ccc, 5/29/09)
!  (20) Minor bug fixes (bmy, phs, 10/9/09)
!  (20) Minor bug fixes (dkh, bmy, 11/19/09)
!  (21) Include second satellite overpass diagnostic.  Adjust AOD name to 550 
!        nm from 400 nm.  Add additional dust AOD bins.  Output values to 
!        hdf_mod. (amv, bmy, 12/1
!******************************************************************************
!
      IMPLICIT NONE

      !=================================================================
      ! MODULE PRIVATE DECLARATIONS -- keep certain internal variables 
      ! and routines from being seen outside "gamap_mod.f"
      !=================================================================

      ! Make everything PRIVATE ...
      PRIVATE

      ! ... except these routines
      PUBLIC :: DO_GAMAP

      !=================================================================
      ! MODULE VARIABLES
      !=================================================================

      ! For "diaginfo.dat"
      INTEGER,           PARAMETER   :: MAXCAT    = 150
      INTEGER,           PARAMETER   :: SPACING   = 1000
      INTEGER                        :: NCATS
      INTEGER,           ALLOCATABLE :: OFFSET(:)
      CHARACTER(LEN=40), ALLOCATABLE :: CATEGORY(:)
      CHARACTER(LEN=40), ALLOCATABLE :: DESCRIPT(:)
      CHARACTER(LEN=255)             :: DFILE

      ! For "tracerinfo.dat"
      INTEGER,           PARAMETER   :: MAXDIAG   = 70   
      INTEGER,           PARAMETER   :: MAXTRACER = 120  
      INTEGER,           ALLOCATABLE :: NTRAC(:)
      INTEGER,           ALLOCATABLE :: INDEX(:,:)
      INTEGER,           ALLOCATABLE :: MOLC(:,:)
      REAL*4,            ALLOCATABLE :: MWT(:,:)
      REAL*4,            ALLOCATABLE :: SCALE(:,:)
      CHARACTER(LEN=40), ALLOCATABLE :: NAME(:,:)
      CHARACTER(LEN=40), ALLOCATABLE :: FNAME(:,:)
      CHARACTER(LEN=40), ALLOCATABLE :: UNIT(:,:)
      CHARACTER(LEN=255)             :: TFILE
      
      ! Other variables
      CHARACTER(LEN=16)              :: STAMP
      CHARACTER(LEN=40)              :: SIM_NAME

      !=================================================================
      ! MODULE ROUTINES -- follow below the "CONTAINS" statement
      !=================================================================
      CONTAINS

!------------------------------------------------------------------------------

      SUBROUTINE DO_GAMAP( DIAGINFO, TRACERINFO )  
!
!******************************************************************************
!  Subroutine DO_GAMAP is the driver program for creating the customized GAMAP
!  files "diaginfo.dat" and "tracerinfo.dat". (bmy, 5/3/05)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) DIAGINFO   (CHARACTER) : Path name of the GAMAP "diaginfo.dat"   file
!  (2 ) TRACERINFO (CHARACTER) : Path name of the GAMAP "tracerinfo.dat" file
!
!  NOTES: 
!******************************************************************************
!
      ! References to F90 modules
      USE TIME_MOD,   ONLY : SYSTEM_TIMESTAMP
      USE TRACER_MOD, ONLY : GET_SIM_NAME 

      ! Arguments
      CHARACTER(LEN=255), INTENT(IN) :: DIAGINFO  
      CHARACTER(LEN=255), INTENT(IN) :: TRACERINFO

      !=================================================================
      ! DO_GAMAP begins here!
      !=================================================================

      ! Allocate and initialize variables
      CALL INIT_GAMAP( DIAGINFO, TRACERINFO )

      ! Create a timestamp with the system date & time
      STAMP    = SYSTEM_TIMESTAMP()

      ! Get simulation name
      SIM_NAME = GET_SIM_NAME()

      ! Write "diaginfo.dat" file
      CALL CREATE_DINFO

      ! Write "tracerinfo.dat" file
      CALL CREATE_TINFO

      ! Deallocate variables
      CALL CLEANUP_GAMAP

      ! Return to calling program
      END SUBROUTINE DO_GAMAP

!------------------------------------------------------------------------------

      SUBROUTINE CREATE_DINFO
!
!******************************************************************************
!  Subroutine CREATE_DINFO  writes information about diagnostic categories
!  to a customized "diaginfo.dat" file. (bmy, 5/3/05)
!
!  NOTES:
!******************************************************************************
!                                
      ! References to F90 modules
      USE FILE_MOD, ONLY : IOERROR, IU_FILE

      ! Local variables
      INTEGER           :: IOS, N

      !=================================================================
      ! CREATE_DINFO begins here!
      !=================================================================

      ! Open "diaginfo.dat" file for output
      OPEN( IU_FILE, FILE=TRIM( DFILE ), STATUS='UNKNOWN', IOSTAT=IOS )
      IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_FILE, 'create_dinfo:1' )

      ! Write file header
      WRITE( IU_FILE, '(a)' ) '#' // REPEAT( '=', 78 )
      WRITE( IU_FILE,  100  ) STAMP
      WRITE( IU_FILE,  105  ) TRIM( SIM_NAME )
      WRITE( IU_FILE,  110  )
      WRITE( IU_FILE,  115  )
      WRITE( IU_FILE, '(a)' ) '# File Format:'
      WRITE( IU_FILE, '(a)' ) '# ' // REPEAT( '-', 77 )
      WRITE( IU_FILE,  120  )
      WRITE( IU_FILE,  125  )
      WRITE( IU_FILE,  130  )
      WRITE( IU_FILE,  135  )
      WRITE( IU_FILE,  125  )
      WRITE( IU_FILE, '(a)' ) '#'
      WRITE( IU_FILE,  140  ) SPACING
      WRITE( IU_FILE, '(a)' ) '#' // REPEAT( '=', 78 )

      ! Loop over categories
      DO N = 1, NCATS

         ! Write one line to "diaginfo.dat" file
         WRITE( IU_FILE, 145, IOSTAT=IOS ) 
     &        OFFSET(N), ADJUSTL( CATEGORY(N) ), ADJUSTL( DESCRIPT(N) )

         ! Error check
         IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_FILE, 'create_dinfo:1' )
      ENDDO

      ! FORMAT strings
 100  FORMAT( '# diaginfo.dat: Created by GEOS-CHEM at ', a,     /,'#' )
 105  FORMAT( '# ****** CUSTOMIZED FOR ', a, ' SIMULATION *****',/,'#' )
 110  FORMAT( '# This file contains category names and the offsets',
     &        ' which they are stored'                                 ) 
 115  FORMAT( '# in file "tracerinfo.dat".  This file is read into',
     &        ' GAMAP by routine', /, '# "ctm_diaginfo.pro".',   /,'#' )
 120  FORMAT( '# OFFSET    (I8 )  Constant to add to tracer numbers', 
     &        ' in order to distinguish', /, '#', 18x, 'for the given',
     &        ' diagnostic category, as stored in file', /, '#', 18x, 
     &        '"tracerinfo.dat".  OFFSET may be up to 8 digits long.'  )
 125  FORMAT( '#  --       (1X )  1-character spacer'                  )    
 130  FORMAT( '# CATEGORY  (A40)  Category name for CTM diagnostics.', 
     &        ' NOTE: The category name', /, '#', 18x, 
     &        'can be up to 40 chars long, but historically the', 
     &        ' GEOS-CHEM', /,'#', 18x, 'and GISS models have used an', 
     &        ' 8-character category name.'                            )
 135  FORMAT( '# COMMENT   (A  )  Descriptive comment string',   /,'#' )
 140  FORMAT('##### SPACING BETWEEN DIAGNOSTIC CATEGORY OFFSETS = ',i8 )
 145  FORMAT( i8, 1x, a40, 1x, a )

      ! Close file
      CLOSE( IU_FILE )

      ! Return to calling program
      END SUBROUTINE CREATE_DINFO

!------------------------------------------------------------------------------

      SUBROUTINE CREATE_TINFO
!
!******************************************************************************
!  Subroutine CREATE_TINFO writes information about tracers to a customized
!  "tracerinfo.dat" file. (bmy, 4/21/05, 2/6/07)
!
!  NOTES:
!  (1 ) Now write out tracers in ug/m3 (dkh, bmy, 5/22/06)
!  (2 ) Now write out GPROD & APROD info (tmf, havala, bmy, 2/6/07)
!******************************************************************************
!      
      ! References to F90 modules
      USE FILE_MOD,    ONLY : IOERROR, IU_FILE
      USE LOGICAL_MOD, ONLY : LSOA

#     include "CMN_SIZE"    ! Size parameters
#     include "CMN_DIAG"    ! NDxx flags

      ! Local variables
      INTEGER              :: D, IOS, N, T
      REAL*4               :: SCALE_NEW
      CHARACTER(LEN=2)     :: C
      CHARACTER(LEN=40)    :: UNIT_NEW, NAME_NEW

      !=================================================================
      ! CREATE_TINFO begins here!
      !=================================================================

      ! Open "tracerinfo.dat" file for output
      OPEN( IU_FILE, FILE=TRIM( TFILE ), STATUS='UNKNOWN', IOSTAT=IOS )
      IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_FILE, 'create_tinfo:1' )

      ! Write file header
      WRITE( IU_FILE, '(a)' ) '#' // REPEAT( '=', 78 )
      WRITE( IU_FILE,  100  ) STAMP
      WRITE( IU_FILE,  105  ) TRIM( SIM_NAME )
      WRITE( IU_FILE,  110  ) 
      WRITE( IU_FILE,  115  )
      WRITE( IU_FILE, '(a)' ) '# File Format:' 
      WRITE( IU_FILE, '(a)' ) '# ' // REPEAT( '-', 77 )    
      WRITE( IU_FILE,  120  )
      WRITE( IU_FILE,  125  )
      WRITE( IU_FILE,  130  )
      WRITE( IU_FILE,  135  ) 
      WRITE( IU_FILE,  140  )
      WRITE( IU_FILE,  145  )
      WRITE( IU_FILE,  150  )
      WRITE( IU_FILE,  125  )
      WRITE( IU_FILE,  155  )

      ! FORMAT strings
 100  FORMAT( '# tracerinfo.dat: Created by GEOS-CHEM at ', a,   /,'#' )
 105  FORMAT( '# ****** CUSTOMIZED FOR ', a, ' SIMULATION *****',/,'#' )
 110  FORMAT( '# This file contains name weight and index', 
     &        ' information about GEOS-CHEM'                           )
 115  FORMAT( '# tracers.  It is read by routine ',
     &        '"ctm_tracerinfo.pro" of the GAMAP package.',      /,'#' )
 120  FORMAT( '# NAME     (A8   )  Tracer name (up to 8 chars)'        )
 125  FORMAT( '#  --      (1X   )  1-character spacer'                 )
 130  FORMAT( '# FULLNAME (A30  )  Full tracer name (up to 30 chars)'  )
 135  FORMAT( '# MOLWT    (E10.0)  Molecular weight (kg/mole)'         )
 140  FORMAT( '# C        (I3   )  For HC''s: # moles C/moles tracer;',
     &        ' otherwise set=1'                                       )
 145  FORMAT( '# TRACER   (I9   )  Tracer number (up to 9 digits)'     )
 150  FORMAT( '# SCALE    (E10.3)  Standard scale factor to convert', 
     &        ' to unit given below'                                   )
 155  FORMAT( '# UNIT     (A40  )  Unit string',                 /,'#' )

      !-----------------------------------
      ! 0: Tracers [ppbv]
      !-----------------------------------

      ! Write separator line
      CALL WRITE_SEPARATOR( 0 )

      ! Loop over tracers
      DO T = 1, NTRAC(45)

         ! GAMAP tracer number
         N = ( SPACING * 0 ) + T

         ! Write tracers [ppbv] to "tracerinfo.dat" file
         CALL WRITE_TINFO( NAME(T,45), FNAME(T,45), MWT(T,45), 
     &                     MOLC(T,45), SCALE(T,45), UNIT(T,45), N )
      ENDDO

      !-----------------------------------
      ! SPACING*1: Tracers [molec/cm2/s]
      !-----------------------------------

      ! Write separator line
      CALL WRITE_SEPARATOR( 100 )

      ! Loop over tracers
      DO T = 1, NTRAC(45)
         
         ! GAMAP tracer number
         N = ( SPACING * 1 ) + T

         ! New scale
         SCALE_NEW = 1.0e0
         
         ! New unit
         IF ( TRIM( UNIT(T,45) ) == 'ppbC' ) THEN
            UNIT_NEW = 'atoms C/cm2/s'
         ELSE
            UNIT_NEW = 'molec/cm2/s'
         ENDIF

         ! Write tracers [molec/cm2/s] to "tracerinfo.dat"
         CALL WRITE_TINFO( NAME(T,45), FNAME(T,45), MWT(T,45), 
     &                     MOLC(T,45), SCALE_NEW,   UNIT_NEW, N )
      ENDDO

      !-----------------------------------
      ! SPACING*2: Tracers [molec/cm2]
      !-----------------------------------

      ! Write separator line
      CALL WRITE_SEPARATOR( 200 )

      ! Loop over tracers
      DO T = 1, NTRAC(45)
         
         ! GAMAP tracer number
         N = ( SPACING * 2 ) + T

         ! New scale
         SCALE_NEW = 1.0e0
         
         ! New unit
         IF ( TRIM( UNIT(T,45) ) == 'ppbC' ) THEN
            UNIT_NEW = 'atoms C/cm2'
         ELSE
            UNIT_NEW = 'molec/cm2'
         ENDIF

         ! Write tracers [molec/cm2] to "tracerinfo.dat"
         CALL WRITE_TINFO( NAME(T,45), FNAME(T,45), MWT(T,45), 
     &                     MOLC(T,45), SCALE_NEW,   UNIT_NEW, N )
      ENDDO

      !-----------------------------------
      ! SPACING*3: Tracers [kg/s]
      !-----------------------------------

      ! Write separator line
      CALL WRITE_SEPARATOR( 300 )

      ! Loop over tracers
      DO T = 1, NTRAC(45)

         ! GAMAP tracer number
         N = ( SPACING * 3 ) + T

         ! New scale
         SCALE_NEW = 1.0e0

         ! New unit
         IF ( TRIM( UNIT(T,45) ) == 'ppbC' ) THEN
            UNIT_NEW = 'kg C/s'
         ELSE
            UNIT_NEW = 'kg/s'
         ENDIF

         ! Write tracers [kg/s] to "tracerinfo.dat"
         CALL WRITE_TINFO( NAME(T,45), FNAME(T,45), MWT(T,45), 
     &                     MOLC(T,45), SCALE_NEW,   UNIT_NEW, N )
      ENDDO

      !-----------------------------------
      ! SPACING*4: Tracers [kg]
      !-----------------------------------

      ! Write separator line
      CALL WRITE_SEPARATOR( 400 )

      ! Loop over tracers
      DO T = 1, NTRAC(45)

         ! GAMAP tracer number
         N = ( SPACING * 4 ) + T

         ! New scale
         SCALE_NEW = 1.0e0

         ! New unit
         IF ( TRIM( UNIT(T,45) ) == 'ppbC' ) THEN
            UNIT_NEW = 'kg C'
         ELSE
            UNIT_NEW = 'kg'
         ENDIF

         ! Write tracers [kg] to "tracerinfo.dat"
         CALL WRITE_TINFO( NAME(T,45), FNAME(T,45), MWT(T,45), 
     &                     MOLC(T,45), SCALE_NEW,   UNIT_NEW, N )
      ENDDO

! Add this later on (bmy, 5/22/06)
!      !-----------------------------------
!      ! SPACING*5: Tracers [ug/m3]
!      !-----------------------------------
!
!      ! Write separator line
!      CALL WRITE_SEPARATOR( 500 )
!
!      ! Loop over tracers
!      DO T = 1, NTRAC(45)
!
!         ! GAMAP tracer number
!         N = ( SPACING * 5 ) + T
!
!         ! New scale
!         SCALE_NEW = 1.0e0
!
!         ! New unit
!         IF ( TRIM( UNIT(T,45) ) == 'ppbC' ) THEN
!            UNIT_NEW = 'ug C/m3'
!         ELSE
!            UNIT_NEW = 'ug/m3'
!         ENDIF
!
!         ! Write tracers [kg] to "tracerinfo.dat"
!         CALL WRITE_TINFO( NAME(T,45), FNAME(T,45), MWT(T,45), 
!     &                     MOLC(T,45), SCALE_NEW,   UNIT_NEW, N )
!      ENDDO

      !-----------------------------------
      ! SPACING*6: GPROD & APROD [kg/kg]
      !-----------------------------------
      IF ( LSOA ) THEN

         ! Write separator line
         CALL WRITE_SEPARATOR( 600 )

         ! Loop over tracers
         DO T = 1, 18
            
            ! GAMAP tracer number
            N = ( SPACING * 6 ) + T
            
            ! Make a character 
            WRITE( C, '(i2.2)' ) T

            ! Tracer name
            NAME_NEW = 'PROD' // C
            
            ! Write tracers [kg] to "tracerinfo.dat"
            CALL WRITE_TINFO( NAME_NEW, NAME_NEW,  1e0, 
     &                        1,        1e0,      'kg/kg', N )
         ENDDO
      ENDIF

      !------------------------------
      ! All other diagnostics 
      !------------------------------
      DO D = 1, MAXDIAG

         ! If tracers are defined then...
         IF ( NTRAC(D) > 0 ) THEN

            ! Skip ND45, we already wrote tracers above
            IF ( D == 45 ) CYCLE

            ! Write separator
            CALL WRITE_SEPARATOR( D )

            ! Write tracers to file
            DO T = 1, NTRAC(D)
               CALL WRITE_TINFO( NAME(T,D), FNAME(T,D), MWT(T,D), 
     &                           MOLC(T,D), SCALE(T,D), UNIT(T,D),
     &                           INDEX(T,D) )
            ENDDO
         ENDIF
      ENDDO

      ! Return to calling program
      END SUBROUTINE CREATE_TINFO

!------------------------------------------------------------------------------

      SUBROUTINE WRITE_TINFO( NAME, FNAME, MWT, MOLC, SCALE, UNIT, N )
!
!******************************************************************************
!  Subroutine WRITE_TINFO writes one line to the customized "tracerinfo.dat" 
!  file. (bmy, 5/3/05)
! 
!  Arguments as Input:
!  ============================================================================
!  (1 ) NAME  (CHARACTER) : GAMAP short tracer name
!  (2 ) FNAME (CHARACTER) : GAMAP long tracer name
!  (3 ) MWT   (REAL*8   ) : Molecular weight [kg/mole]
!  (4 ) MOLC  (REAL*8   ) : Moles C per mole tracer (for hydrocarbons)
!  (5 ) SCALE (REAL*8   ) : GAMAP scale factor (e.g. 1e+9 scales v/v --> ppbv)
!  (6 ) UNIT  (CHARACTER) : Unit string
!  (7 ) N     (INTEGER  ) : Tracer number
!
!  NOTES:
!******************************************************************************
!
      ! References to F90 modules
      USE FILE_MOD, ONLY : IU_FILE, IOERROR

      ! Argumetns
      INTEGER,          INTENT(IN) :: MOLC, N
      REAL*4,           INTENT(IN) :: MWT,  SCALE
      CHARACTER(LEN=*), INTENT(IN) :: NAME, FNAME, UNIT

      ! Local variables
      INTEGER                      :: IOS

      !=================================================================
      ! WRITE_TINFO begins here!
      !=================================================================

      ! Write one line to "tracerinfo.dat" file
      WRITE( IU_FILE, 100, IOSTAT=IOS ) 
     &     ADJUSTL( NAME ), ADJUSTL( FNAME ), MWT, 
     &     MOLC,            N,                SCALE, TRIM( UNIT )

      ! Error check
      IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_FILE, 'write_tinfo:1' )

      ! FORMAT string
 100  FORMAT( a8, 1x, a30, es10.3, i3, i9, es10.3, 1x, a )

      ! Return to calling program
      END SUBROUTINE WRITE_TINFO

!------------------------------------------------------------------------------

      SUBROUTINE WRITE_SEPARATOR( DIAG )
!
!******************************************************************************
!  Subroutine WRITE_SEPARATOR writes a separator block to the customized
!  "tracerinfo.dat" file. (bmy, 5/3/05, 2/6/07)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) DIAG (INTEGER) : GEOS-CHEM diagnostic number 
!
!  NOTES:
!  (1 ) Added new header for GPROD & APROD info (bmy, 2/6/07)
!******************************************************************************
!
      ! References to F90 modules
      USE FILE_MOD, ONLY : IU_FILE, IOERROR

      ! Arguments
      INTEGER, INTENT(IN) :: DIAG

      ! Local variables
      INTEGER             :: IOS
      CHARACTER(LEN=79)   :: SEPARATOR

      !=================================================================
      ! WRITE_SEPARATOR begins here! 
      !=================================================================

      ! Create separator string
      SEPARATOR = '#' // REPEAT( '=', 78 )

      ! Write separator string
      WRITE( IU_FILE, '(a)', IOSTAT=IOS ) SEPARATOR
      IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_FILE, 'write_separator:1' )

      ! Write the appropriate message
      SELECT CASE( DIAG )
         CASE( 0 ) 
            WRITE( IU_FILE, 100, IOSTAT=IOS )
         CASE( 100 ) 
            WRITE( IU_FILE, 110, IOSTAT=IOS )
         CASE( 200 ) 
            WRITE( IU_FILE, 120, IOSTAT=IOS )
         CASE( 300 ) 
            WRITE( IU_FILE, 130, IOSTAT=IOS )
         CASE( 400 ) 
            WRITE( IU_FILE, 140, IOSTAT=IOS )
         CASE( 500 ) 
            WRITE( IU_FILE, 150, IOSTAT=IOS )
         CASE( 600 ) 
            WRITE( IU_FILE, 160, IOSTAT=IOS )
         CASE DEFAULT
            WRITE( IU_FILE, 170, IOSTAT=IOS ) DIAG
      END SELECT

      ! Error check
      IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_FILE, 'write_separator:2' )

      ! Write separator string
      WRITE( IU_FILE, '(a)', IOSTAT=IOS ) SEPARATOR
      IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_FILE, 'write_separator:3' )

      ! FORMAT strings
 100  FORMAT( '# GEOS-CHEM tracers [ppbv]'           )
 110  FORMAT( '# GEOS-CHEM tracers [molec/cm2/s]'    )
 120  FORMAT( '# GEOS-CHEM tracers [molec/cm2]'      )
 130  FORMAT( '# GEOS-CHEM tracers [kg/s]'           )
 140  FORMAT( '# GEOS-CHEM tracers [kg]'             )
 150  FORMAT( '# GEOS-CHEM tracers [ug/m3]'          )
 160  FORMAT( '# SOA GPROD & APROD [kg/kg]'          )
 170  FORMAT( '# ND', i2.2, ' diagnostic quantities' )
      
      ! Return to calling program
      END SUBROUTINE WRITE_SEPARATOR

!------------------------------------------------------------------------------

      SUBROUTINE INIT_DIAGINFO
!
!******************************************************************************
!  Subroutine INIT_DIAGINFO initializes the CATEGORY, DESCRIPT, and OFFSET
!  variables, which are used to define the "diaginfo.dat" file for GAMAP.
!  (bmy, 10/17/06, 11/16/07)
!
!  NOTES:
!  (1 ) Split this code off from INIT_GAMAP, for clarity.  Now declare biomass
!        burning emissions w/ offset of 45000.  Now declare time in the 
!        troposphere diagnostic with offset of 46000. (phs, bmy, 10/17/06)
!  (2 ) Now add IJ-GPROD & IJ-APROD w/ offset of SPACING*6, for the SOA
!        GPROD & APROD restart file. (tmf, havala, bmy, 2/6/07)
!  (3 ) Now declare H2-HD sources w/ offset of 48000. Now declare H2-HD
!        production/loss w/ offset of 47000. (phs, 9/18/07)
!  (4 ) Change diagnostic category for ND31 diagnostic from "PS-PTOP" 
!        to "PEDGE-$" (bmy, 11/16/07)
!  (5 ) Add categories CH4-LOSS, CH4-EMISS and WET-FRAC (kjw, 8/18/09)
!******************************************************************************
!
      ! Local variables 
      INTEGER :: N

      !=================================================================
      ! INIT_DIAGINFO begins here!
      !=================================================================

      N           = 1
      CATEGORY(N) = 'IJ-AVG-$'
      DESCRIPT(N) = 'Tracer concentration'
      OFFSET(N)   = SPACING * 0

      N           = N + 1
      CATEGORY(N) = 'IJ-24H-$'
      DESCRIPT(N) = '24-hr avg tracer conc.'
      OFFSET(N)   = SPACING * 0
         
      N           = N + 1
      CATEGORY(N) = 'INST-MAP'
      DESCRIPT(N) = 'Instantaneous tracer'
      OFFSET(N)   = SPACING * 0

      N           = N + 1
      CATEGORY(N) = 'ANTHSRCE'
      DESCRIPT(N) = 'Anthropogenic emissions'
      OFFSET(N)   = SPACING * 1

      N           = N + 1
      CATEGORY(N) = 'BIOFSRCE'
      DESCRIPT(N) = 'Biofuel emissions'
      OFFSET(N)   = SPACING * 1

      N           = N + 1
      CATEGORY(N) = 'NOX-AC-$'
      DESCRIPT(N) = 'Aircraft NOx'
      OFFSET(N)   = SPACING * 1

      N           = N + 1
      CATEGORY(N) = 'NOX-AN-$'
      DESCRIPT(N) = 'Anthropogenic NOx'
      OFFSET(N)   = SPACING * 1

      N           = N + 1
      CATEGORY(N) = 'NOX-BIOB'
      DESCRIPT(N) = 'Biomass NOx'
      OFFSET(N)   = SPACING * 1

      N           = N + 1
      CATEGORY(N) = 'NOX-BIOF'
      DESCRIPT(N) = 'Biofuel NOx'
      OFFSET(N)   = SPACING * 1

      N           = N + 1
      CATEGORY(N) = 'NOX-LI-$'
      DESCRIPT(N) = 'Lightning NOx'
      OFFSET(N)   = SPACING * 1

      N           = N + 1
      CATEGORY(N) = 'NOX-SOIL'
      DESCRIPT(N) = 'Soil NOx'
      OFFSET(N)   = SPACING * 1

      N           = N + 1
      CATEGORY(N) = 'NOX-FERT'
      DESCRIPT(N) = 'Fertilizer NOx'
      OFFSET(N)   = SPACING * 1

      N           = N + 1
      CATEGORY(N) = 'NOX-STRT'
      DESCRIPT(N) = 'Stratopsheric NOx'
      OFFSET(N)   = SPACING * 1
 
      N           = N + 1
      CATEGORY(N) = 'INST_COL'
      DESCRIPT(N) = 'Instantaneous columns'
      OFFSET(N)   = SPACING * 2

      N           = N + 1
      CATEGORY(N) = 'CV-FLX-$'
      DESCRIPT(N) = 'Convective mass flux'
      OFFSET(N)   = SPACING * 3

      N           = N + 1
      CATEGORY(N) = 'TURBMC-$'
      DESCRIPT(N) = 'PBL mixing mass flux'
      OFFSET(N)   = SPACING * 3

      N           = N + 1
      CATEGORY(N) = 'EW-FLX-$'
      DESCRIPT(N) = 'E/W transport flux'
      OFFSET(N)   = SPACING * 3

      N           = N + 1
      CATEGORY(N) = 'NS-FLX-$'
      DESCRIPT(N) = 'N/S transport flux'
      OFFSET(N)   = SPACING * 3

      N           = N + 1
      CATEGORY(N) = 'UP-FLX-$'
      DESCRIPT(N) = 'Up/down transport flux'
      OFFSET(N)   = SPACING * 3

      N           = N + 1
      CATEGORY(N) = 'STRT-FLX'
      DESCRIPT(N) = 'Flux from stratosphere'
      OFFSET(N)   = SPACING * 3

      N           = N + 1
      CATEGORY(N) = 'RN--SRCE'
      DESCRIPT(N) = 'Rn-Pb-Be source'
      OFFSET(N)   = SPACING * 3

      N           = N + 1
      CATEGORY(N) = 'RN-DECAY'
      DESCRIPT(N) = 'Rn-Pb-Be loss'
      OFFSET(N)   = SPACING * 3

      N           = N + 1
      CATEGORY(N) = 'WETDCV-$'
      DESCRIPT(N) = 'Conv wet scavenging'
      OFFSET(N)   = SPACING * 3

      N           = N + 1
      CATEGORY(N) = 'WETDLS-$'
      DESCRIPT(N) = 'Wet deposition'
      OFFSET(N)   = SPACING * 3

      N           = N + 1
      CATEGORY(N) = 'DMS-BIOG'
      DESCRIPT(N) = 'Biogenic DMS'
      OFFSET(N)   = SPACING * 4

      N           = N + 1
      CATEGORY(N) = 'DUSTSRCE'
      DESCRIPT(N) = 'Dust emission'
      OFFSET(N)   = SPACING * 4

      N           = N + 1
      CATEGORY(N) = 'NVOCSRCE'
      DESCRIPT(N) = 'NVOC emissions'
      OFFSET(N)   = SPACING * 4

      N           = N + 1
      CATEGORY(N) = 'SALTSRCE'
      DESCRIPT(N) = 'Seasalt emission'
      OFFSET(N)   = SPACING * 4

      N           = N + 1
      CATEGORY(N) = 'SO2-AC-$'
      DESCRIPT(N) = 'Aircraft SO2 emissions'
      OFFSET(N)   = SPACING * 4

      N           = N + 1
      CATEGORY(N) = 'SO2-AN-$'
      DESCRIPT(N) = 'Anthro SO2 emissions'
      OFFSET(N)   = SPACING * 4

      N           = N + 1
      CATEGORY(N) = 'SO2-BIOB'
      DESCRIPT(N) = 'Biomass SO2 emissions'
      OFFSET(N)   = SPACING * 4

      N           = N + 1
      CATEGORY(N) = 'SO2-BIOF'
      DESCRIPT(N) = 'Biofuel SO2 emissions'
      OFFSET(N)   = SPACING * 4 

      N           = N + 1
      CATEGORY(N) = 'SO2-EV-$'
      DESCRIPT(N) = 'Erup. Volcano SO2'
      OFFSET(N)   = SPACING * 4

      N           = N + 1
      CATEGORY(N) = 'SO2-NV-$'
      DESCRIPT(N) = 'Non-Erup. Volcano SO2'
      OFFSET(N)   = SPACING * 4

      N           = N + 1
      CATEGORY(N) = 'SO2-SHIP'
      DESCRIPT(N) = 'SO2 from ship exhaust'
      OFFSET(N)   = SPACING * 4

      N           = N + 1
      CATEGORY(N) = 'SO4-AN-$'
      DESCRIPT(N) = 'Anthro SO4 emissions'
      OFFSET(N)   = SPACING * 4

      N           = N + 1
      CATEGORY(N) = 'SO4-BIOF'
      DESCRIPT(N) = 'Biofuel SO4 emissions'
      OFFSET(N)   = SPACING * 4

      N           = N + 1
      CATEGORY(N) = 'NH3-ANTH'
      DESCRIPT(N) = 'Anthro NH3 emissions'
      OFFSET(N)   = SPACING * 4

      N           = N + 1
      CATEGORY(N) = 'NH3-NATU'
      DESCRIPT(N) = 'Natural NH3 emissions'
      OFFSET(N)   = SPACING * 4

      N           = N + 1
      CATEGORY(N) = 'NH3-BIOB'
      DESCRIPT(N) = 'Biomass NH3 emissions'
      OFFSET(N)   = SPACING * 4

      N           = N + 1
      CATEGORY(N) = 'NH3-BIOF'
      DESCRIPT(N) = 'Biofuel NH3 emissions'
      OFFSET(N)   = SPACING * 4

      N           = N + 1
      CATEGORY(N) = 'TROPO-AV'
      DESCRIPT(N) = 'Trop avg''d tracer'
      OFFSET(N)   = SPACING * 4

      N           = N + 1
      CATEGORY(N) = 'TCMASS-$'
      DESCRIPT(N) = 'Tracer mass (kg)'
      OFFSET(N)   = SPACING * 4

      N           = N + 1
      CATEGORY(N) = 'CV-FLX-$'
      DESCRIPT(N) = 'Upward flux from wet conv'
      OFFSET(N)   = SPACING * 4

      N           = N + 1
      CATEGORY(N) = 'TURBMC-$'
      DESCRIPT(N) = 'Upward flux from PBL mixing'
      OFFSET(N)   = SPACING * 4

      N           = N + 1
      CATEGORY(N) = 'EW-FLX-$'
      DESCRIPT(N) = 'E/W transport flux'
      OFFSET(N)   = SPACING * 4

      N           = N + 1
      CATEGORY(N) = 'NS-FLX-$'
      DESCRIPT(N) = 'N/S transport flux'
      OFFSET(N)   = SPACING * 4

      N           = N + 1
      CATEGORY(N) = 'UP-FLX-$'
      DESCRIPT(N) = 'Up/down transport flux'
      OFFSET(N)   = SPACING * 4

      N           = N + 1
      CATEGORY(N) = 'IJ-GPROD'
      DESCRIPT(N) = 'SOA GPROD restart'
      OFFSET(N)   = SPACING * 6

      N           = N + 1
      CATEGORY(N) = 'IJ-APROD'
      DESCRIPT(N) = 'SOA APROD restart'
      OFFSET(N)   = SPACING * 6

      N           = N + 1
      CATEGORY(N) = 'PEDGE-$'
      DESCRIPT(N) = 'Pressure at level edges'
      OFFSET(N)   = SPACING * 10 

      N           = N + 1
      CATEGORY(N) = 'DAO-FLDS'
      DESCRIPT(N) = 'GMAO 2-D met fields'
      OFFSET(N)   = SPACING * 11

      N           = N + 1
      CATEGORY(N) = 'DAO-3D-$'
      DESCRIPT(N) = 'GMAO 3-D met fields'
      OFFSET(N)   = SPACING * 12

      N           = N + 1
      CATEGORY(N) = 'JV-MAP-$'
      DESCRIPT(N) = 'Photolysis rates'
      OFFSET(N)   = SPACING * 13

      N           = N + 1
      CATEGORY(N) = 'OD-MAP-$'
      DESCRIPT(N) = 'Optical Depths'
      OFFSET(N)   = SPACING * 14

      N           = N + 1
      CATEGORY(N) = 'LANDMAP'
      DESCRIPT(N) = 'Land type map'
      OFFSET(N)   = SPACING * 15

      N           = N + 1
      CATEGORY(N) = 'CHEM-L=$'
      DESCRIPT(N) = 'Chemical Prod/Loss'
      OFFSET(N)   = SPACING * 16

      N           = N + 1
      CATEGORY(N) = 'PORL-L=$'
      DESCRIPT(N) = 'ND65 P/L family diagnostics'
      OFFSET(N)   = SPACING * 17

      N           = N + 1
      CATEGORY(N) = 'PL-SUL=$'
      DESCRIPT(N) = 'P/L of sulfur species'
      OFFSET(N)   = SPACING * 18

      N           = N + 1
      CATEGORY(N) = 'TIME-SER'
      DESCRIPT(N) = 'Timeseries quantities'
      OFFSET(N)   = SPACING * 19

      N           = N + 1
      CATEGORY(N) = 'CO--SRCE'
      DESCRIPT(N) = 'CO Source diagnostic'
      OFFSET(N)   = SPACING * 20

      N           = N + 1
      CATEGORY(N) = 'BIOGSRCE'
      DESCRIPT(N) = 'Biogenic emissions'
      OFFSET(N)   = SPACING * 21

      N           = N + 1
      CATEGORY(N) = 'ACETSRCE'
      DESCRIPT(N) = 'Acetone emissions'
      OFFSET(N)   = SPACING * 22

      N           = N + 1
      CATEGORY(N) = 'EMDIS-BL'
      DESCRIPT(N) = 'Emissions in PBL'
      OFFSET(N)   = SPACING * 23

      N           = N + 1
      CATEGORY(N) = 'BXHGHT-$'
      DESCRIPT(N) = 'Boxheight, airmass, etc'
      OFFSET(N)   = SPACING * 24

      N           = N + 1
      CATEGORY(N) = 'DXYP'
      DESCRIPT(N) = 'Surface area'
      OFFSET(N)   = SPACING * 25

      N           = N + 1
      CATEGORY(N) = 'TR-PAUSE'
      DESCRIPT(N) = 'Annual mean tropopause'
      OFFSET(N)   = SPACING * 26

      N           = N + 1
      CATEGORY(N) = 'PBLDEPTH'
      DESCRIPT(N) = 'Afternoon PBL height'
      OFFSET(N)   = SPACING * 27

      N           = N + 1
      CATEGORY(N) = 'WD-FRC-$'
      DESCRIPT(N) = 'Wet dep fraction'
      OFFSET(N)   = SPACING * 28

      N           = N + 1
      CATEGORY(N) = 'WD-LSR-$'
      DESCRIPT(N) = 'Large-scale rainout'
      OFFSET(N)   = SPACING * 29

      N           = N + 1
      CATEGORY(N) = 'WD-CVR-$'
      DESCRIPT(N) = 'Convective rainout'
      OFFSET(N)   = SPACING * 29

      N           = N + 1
      CATEGORY(N) = 'WD-LSW-$'
      DESCRIPT(N) = 'Large-scale washout'
      OFFSET(N)   = SPACING * 29

      N           = N + 1
      CATEGORY(N) = 'WD-CVW-$'
      DESCRIPT(N) = 'Convective washout'
      OFFSET(N)   = SPACING * 29

      N           = N + 1
      CATEGORY(N) = 'MC-FRC-$'
      DESCRIPT(N) = 'Moist conv fraction'
      OFFSET(N)   = SPACING * 30

      N           = N + 1
      CATEGORY(N) = 'COBUDGET'
      DESCRIPT(N) = 'bnd CO-OH budget'
      OFFSET(N)   = SPACING * 31

      N           = N + 1
      CATEGORY(N) = 'CH4-LOSS'
      DESCRIPT(N) = 'CH4 Loss by OH'
      OFFSET(N)   = SPACING * 32

      N           = N + 1
      CATEGORY(N) = 'BC-ANTH'
      DESCRIPT(N) = 'Anthro BC emission'
      OFFSET(N)   = SPACING * 33

      N           = N + 1
      CATEGORY(N) = 'BC-BIOB'
      DESCRIPT(N) = 'Biomass BC emission'
      OFFSET(N)   = SPACING * 33

      N           = N + 1
      CATEGORY(N) = 'BC-BIOF'
      DESCRIPT(N) = 'Biofuel BC emission'
      OFFSET(N)   = SPACING * 33

      N           = N + 1
      CATEGORY(N) = 'OC-ANTH'
      DESCRIPT(N) = 'Anthro OC emission'
      OFFSET(N)   = SPACING * 33

      N           = N + 1
      CATEGORY(N) = 'OC-BIOB'
      DESCRIPT(N) = 'Biomass OC emission'
      OFFSET(N)   = SPACING * 33

      N           = N + 1
      CATEGORY(N) = 'OC-BIOF'
      DESCRIPT(N) = 'Biofuel OC emission'
      OFFSET(N)   = SPACING * 33

      N           = N + 1
      CATEGORY(N) = 'OC-BIOG'
      DESCRIPT(N) = 'Biogenic OC emission'
      OFFSET(N)   = SPACING * 33

      N           = N + 1
      CATEGORY(N) = 'OC-ALPH'
      DESCRIPT(N) = 'Biogenic ALPH emission'
      OFFSET(N)   = SPACING * 33

      N           = N + 1
      CATEGORY(N) = 'OC-LIMO'
      DESCRIPT(N) = 'Biogenic ALPH emission'
      OFFSET(N)   = SPACING * 33

      N           = N + 1
      CATEGORY(N) = 'OC-TERP'
      DESCRIPT(N) = 'Biogenic TERP emission'
      OFFSET(N)   = SPACING * 33

      N           = N + 1
      CATEGORY(N) = 'OC-ALCO'
      DESCRIPT(N) = 'Biogenic ALCO emission'
      OFFSET(N)   = SPACING * 33

      N           = N + 1
      CATEGORY(N) = 'OC-SESQ'
      DESCRIPT(N) = 'Biogenic SESQ emission'
      OFFSET(N)   = SPACING * 33

      N           = N + 1
      CATEGORY(N) = 'PL-BC=$'
      DESCRIPT(N) = 'H-philic from H-phobic BC'
      OFFSET(N)   = SPACING * 33

      N           = N + 1
      CATEGORY(N) = 'PL-OC=$'
      DESCRIPT(N) = 'H-philic from H-phobic OC'
      OFFSET(N)   = SPACING * 33

      N           = N + 1
      CATEGORY(N) = 'SALT-SR$'
      DESCRIPT(N) = 'Sea salt emission'
      OFFSET(N)   = SPACING * 33

      N           = N + 1
      CATEGORY(N) = 'BLKC-SR$'
      DESCRIPT(N) = 'Black carbon emission'
      OFFSET(N)   = SPACING * 33

      N           = N + 1
      CATEGORY(N) = 'ORGC-SR$'
      DESCRIPT(N) = 'Organic Carbon emission'
      OFFSET(N)   = SPACING * 33

      N           = N + 1
      CATEGORY(N) = 'SOAGM=$'
      DESCRIPT(N) = 'Dicarbonyl SOA production'
      OFFSET(N)   = SPACING * 33

      N           = N + 1
      CATEGORY(N) = 'HG-SRCE'
      DESCRIPT(N) = 'Hg emissions'
      OFFSET(N)   = SPACING * 34

      N           = N + 1
      CATEGORY(N) = 'PL-HG2-$'
      DESCRIPT(N) = 'Prod / loss of Hg2'
      OFFSET(N)   = SPACING * 35

      N           = N + 1
      CATEGORY(N) = 'DRYD-FLX'
      DESCRIPT(N) = 'Drydep fluxes'
      OFFSET(N)   = SPACING * 36

      N           = N + 1
      CATEGORY(N) = 'DRYD-VEL'
      DESCRIPT(N) = 'Drydep velocities'
      OFFSET(N)   = SPACING * 37

      N           = N + 1
      CATEGORY(N) = 'HCN-PL-$'
      DESCRIPT(N) = 'HCN & CH3CN sinks'
      OFFSET(N)   = SPACING * 38

      N           = N + 1
      CATEGORY(N) = 'HCN-SRCE'
      DESCRIPT(N) = 'HCN & CH3CN sources'
      OFFSET(N)   = SPACING * 39

      N           = N + 1
      CATEGORY(N) = 'CO2-SRCE'
      DESCRIPT(N) = 'CO2 fluxes'
      OFFSET(N)   = SPACING * 40

      N           = N + 1
      CATEGORY(N) = 'OCEAN-HG'
      DESCRIPT(N) = 'Oceanic Hg emissions'
      OFFSET(N)   = SPACING * 41

      N           = N + 1
      CATEGORY(N) = 'LFLASH-$'
      DESCRIPT(N) = 'Lightning flash rates'
      OFFSET(N)   = SPACING * 42

      N           = N + 1
      CATEGORY(N) = 'IJ-SOA-$'
      DESCRIPT(N) = 'SOA concentrations'
      OFFSET(N)   = SPACING * 43

      N           = N + 1
      CATEGORY(N) = 'CH3ISRCE'
      DESCRIPT(N) = 'CH3I emissions'
      OFFSET(N)   = SPACING * 44

      N           = N + 1
      CATEGORY(N) = 'BIOBSRCE'
      DESCRIPT(N) = 'Biomass emissions'
      OFFSET(N)   = SPACING * 45

      N           = N + 1
      CATEGORY(N) = 'TIME-TPS'
      DESCRIPT(N) = 'Fraction of time in troposphere'
      OFFSET(N)   = SPACING * 46

      N           = N + 1
      CATEGORY(N) = 'PL-H2HD-'
      DESCRIPT(N) = 'Prod / loss of H2-HD'
      OFFSET(N)   = SPACING * 47

      N           = N + 1
      CATEGORY(N) = 'H2HD-SRC'
      DESCRIPT(N) = 'H2 HD emissions'
      OFFSET(N)   = SPACING * 48

      ! New ND52 diagnostic for gamma(HO2). (jaegle, 2/29/09)
      N           = N + 1
      CATEGORY(N) = 'GAMMA'
      DESCRIPT(N) = 'gamma HO2'
      OFFSET(N)   = SPACING * 49

      N           = N + 1
      CATEGORY(N) = 'CH4-EMIS'
      DESCRIPT(N) = 'CH4 Emissions'
      OFFSET(N)   = SPACING * 50

      N           = N + 1
      CATEGORY(N) = 'WET-FRAC'
      DESCRIPT(N) = 'Wetland Fraction'
      OFFSET(N)   = SPACING * 51

! NEED TO DO THAT AT ONE POINT      
!      ! New for CSPEC
!      N           = N + 1
!      CATEGORY(N) = 'IJ-CHK-$'
!      DESCRIPT(N) = 'species concentrations after chemistry'
!      OFFSET(N)   = SPACING * 50
!      
      ! Number of categories
      NCATS = N
      
      ! Return to calling program
      END SUBROUTINE INIT_DIAGINFO

!------------------------------------------------------------------------------

      SUBROUTINE INIT_TRACERINFO
!
!******************************************************************************
!  Subroutine INIT_TRACERINFO initializes the NAME, FNAME, MWT, MOLC, INDEX,
!  MOLC, UNIT arrays which are used to define the "tracerinfo.dat" file.
!  (bmy, phs, 10/17/06, 11/19/09)
!
!  NOTES:
!  (1 ) Split this code off from INIT_GAMAP, for clarity.  Also now declare
!        biomass burning emissions w/ offset of 45000.  Bug fix: write out
!        26 tracers for ND48, ND49, ND50, ND51 timeseries.  Also define
!        ND54 diagnostic with offset of 46000. (bmy, 10/17/06)
!  (2 ) Modifications for H2/HD in ND10, ND44 diagnostics (phs, 9/18/07)
!  (3 ) Now write out PBLDEPTH diagnostic information to "tracerinfo.dat" if 
!        any of ND41, ND48, ND49, ND50, ND51 are turned on.  Also set the
!        unit to "kg/s" for the Rn-Pb-Be ND44 drydep diag. (cdh, bmy, 2/22/08)
!  (4 ) Added C2H4 in ND46 (ccc, 2/2/09)
!  (5 ) Add EFLUX to ND67 (lin, ccc, 5/29/08)
!  (6 ) Bug fix in ND28: ALD2 should have 2 carbons, not 3.  Also bug fix
!        in ND66 to print out the name of ZMMU correctly. (dbm, bmy, 10/9/09)
!  (7 ) Previous bug fix was erroneous; now corrected (dkh, bmy, 11/19/09)
!  (8 ) Include second satellite overpass diagnostic.  Adjust AOD name to 550 
!        nm from 400 nm.  Add additional dust AOD bins (amv, bmy, 12/18/09)
!******************************************************************************
!
      ! References to F90 modules
      USE DIAG03_MOD,   ONLY : ND03
      USE DIAG04_MOD,   ONLY : ND04
      USE DIAG41_MOD,   ONLY : ND41
      USE DIAG42_MOD,   ONLY : ND42
      USE DIAG48_MOD,   ONLY : DO_SAVE_DIAG48
      USE DIAG49_MOD,   ONLY : DO_SAVE_DIAG49
      USE DIAG50_MOD,   ONLY : DO_SAVE_DIAG50
      USE DIAG51_MOD,   ONLY : DO_SAVE_DIAG51
      USE DIAG51b_MOD,  ONLY : DO_SAVE_DIAG51b
      USE DIAG56_MOD,   ONLY : ND56
      USE DIAG_PL_MOD,  ONLY : DO_SAVE_PL,  GET_NFAM
      USE DIAG_PL_MOD,  ONLY : GET_FAM_MWT, GET_FAM_NAME
      USE DRYDEP_MOD,   ONLY : DEPNAME,     NUMDEP,    NTRAIND
      USE LOGICAL_MOD,  ONLY : LSOA
      USE TRACER_MOD,   ONLY : ITS_A_CO2_SIM,    ITS_A_H2HD_SIM
      USE TRACER_MOD,   ONLY : ITS_A_CH3I_SIM,   ITS_A_FULLCHEM_SIM
      USE TRACER_MOD,   ONLY : ITS_A_HCN_SIM,    ITS_A_MERCURY_SIM
      USE TRACER_MOD,   ONLY : ITS_A_RnPbBe_SIM, ITS_A_TAGOX_SIM
      USE TRACER_MOD,   ONLY : N_TRACERS,        TRACER_COEFF
      USE TRACER_MOD,   ONLY : TRACER_MW_KG,     TRACER_NAME
      USE TRACERID_MOD, ONLY : IDTBCPI, IDTOCPI, IDTALPH, IDTLIMO
      USE TRACERID_MOD, ONLY : IDTSOA1, IDTSOA2, IDTSOA3, NEMANTHRO
      USE TRACERID_MOD, ONLY : IDTSOA4, IDTSOAM, IDTSOAG
      USE WETSCAV_MOD,  ONLY : GET_WETDEP_IDWETD, GET_WETDEP_NSOL

#     include "CMN_SIZE"     ! Size parameters
#     include "CMN_DIAG"     ! NDxx flags

      ! Local variables
      INTEGER               :: N, NN, NYMDb, NHMSb, T
      LOGICAL               :: DO_TIMESERIES

      !=================================================================
      ! INIT_TRACERINFO begins here!
      !=================================================================

      ! Set a flag if any timeseries diagnostics are turned on
      DO_TIMESERIES = ( DO_SAVE_DIAG48 .or. DO_SAVE_DIAG49 .or.
     &                  DO_SAVE_DIAG50 .or. DO_SAVE_DIAG51 .or.
     &                  DO_SAVE_DIAG51b )

      !----------------------------------
      ! General tracer information
      !----------------------------------
      DO T = 1, N_TRACERS

         ! Store quantities for each tracer
         NAME (T,45) = TRACER_NAME(T)
         FNAME(T,45) = TRIM( NAME(T,45) ) // ' tracer'
         MOLC (T,45) = INT( TRACER_COEFF(T,1) )
         SCALE(T,45) = 1.0e+9
         INDEX(T,45) = T   !changed from N (phs, 3/19/03)

         ! We usually report NOx as Tg N
         IF ( TRIM( NAME(T,45) ) == 'NOx'  .or. 
     &        TRIM( NAME(T,45) ) == 'NOX' ) THEN
            MWT(T,45) = 14e-3
         ELSE
            MWT(T,45) = TRACER_MW_KG(T) 
         ENDIF

         ! Special handling for hydrocarbons
         IF ( MOLC(T,45) > 1 ) THEN
            UNIT(T,45) = 'ppbC'
         ELSE
            UNIT(T,45) = 'ppbv'
         ENDIF

         ! Special handling for Rn-Pb-Be simulation (bmy, 5/11/05)
         IF ( ITS_A_RnPbBe_SIM() ) THEN
            SELECT CASE( T ) 
               CASE( 1 )
                  UNIT (T,45) = 'pCi/SCM'
                  SCALE(T,45) = 1.5243e21
               CASE( 2 ) 
                  UNIT (T,45) = 'fCi/SCM'
                  SCALE(T,45) = 7.0651e20
               CASE( 3 )
                  UNIT (T,45) = 'fCi/SCM'
                  SCALE(T,45) = 1.0938e23
               CASE DEFAULT
                  ! Nothing
            END SELECT
         ENDIF

         ! For mercury, print as pptv (bmy, 1/24/06)
         IF ( ITS_A_MERCURY_SIM() ) THEN
            UNIT(T,45)  = 'pptv'
            SCALE(T,45) = 1.0e+12
         ENDIF
      ENDDO

      ! Number of ND45 tracers 
      NTRAC(45) = N_TRACERS

      !-------------------------------------
      ! FULL-CHEMISTRY RUNS ONLY: 
      ! Pure O3 is stored as N_TRACERS+1
      !-------------------------------------
      IF ( ITS_A_FULLCHEM_SIM() ) THEN
         NTRAC(45)   = N_TRACERS + 1
         T           = N_TRACERS + 1
         NAME (T,45) = 'O3'
         FNAME(T,45) = 'Pure O3 tracer'
         UNIT (T,45) = 'ppbv'
         MWT  (T,45) = TRACER_MW_KG(2)           ! Ox is tracer 2
         MOLC (T,45) = INT( TRACER_COEFF(2,1) )  ! Ox is tracer 2
         SCALE(T,45) = 1.0e+9
         INDEX(T,45) = T
      ENDIF

      !-------------------------------------
      ! Hg source, production & loss (ND03)
      !
      ! Updated for tagged Hg simulation
      ! (cdh, bmy, 1/9/06)
      !-------------------------------------
      IF ( ND03 > 0 ) THEN
         
         ! Number of tracers
         NTRAC(03) = 16 + N_TRACERS

         ! Loop over tracers for HG-SRCE, PL-HG2-$, OCEAN-HG
         DO T = 1, NTRAC(03)

            ! Define quantities
            UNIT (T,03) = 'kg'
            MOLC (T,03) = 1
            MWT  (T,03) = 201e-3
            SCALE(T,03) = 1e0

            ! Get name, long-name, index, and new units
            SELECT CASE( T )
               CASE( 1  )
                  NAME (T,03) = 'Hg0_an'
                  FNAME(T,03) = 'Anthro elemental Hg'
                  INDEX(T,03) = T + ( SPACING * 34 )
               CASE( 2  )
                  NAME (T,03) = 'Hg0_aq'
                  FNAME(T,03) = 'Ocean mass of elemental Hg'
                  INDEX(T,03) = T + ( SPACING * 34 )
               CASE( 3  )
                  NAME (T,03) = 'Hg0_oc'
                  FNAME(T,03) = 'Ocean-emitted elemental Hg'
                  INDEX(T,03) = T + ( SPACING * 34 )
               CASE( 4  )
                  NAME (T,03) = 'Hg0_ln'
                  FNAME(T,03) = 'Land re-emitted elemental Hg'
                  INDEX(T,03) = T + ( SPACING * 34 )
               CASE( 5  )
                  NAME (T,03) = 'Hg0_na'
                  FNAME(T,03) = 'Natural land source'
                  INDEX(T,03) = T + ( SPACING * 34 )
               CASE( 6  )
                  NAME (T,03) = 'Hg2_an'
                  FNAME(T,03) = 'Anthro divalent Hg'
                  INDEX(T,03) = T + ( SPACING * 34 )
               CASE( 7  )
                  NAME (T,03) = 'Hg2_aq'
                  FNAME(T,03) = 'Ocean mass of divalent Hg'
                  INDEX(T,03) = T + ( SPACING * 34 )
               CASE( 8  )
                  NAME (T,03) = 'Hg2_sk'
                  FNAME(T,03) = 'Mass of Hg2 sunk in ocean'
                  INDEX(T,03) = T + ( SPACING * 34 )
               CASE( 9  )
                  NAME (T,03) = 'HgP_an'
                  FNAME(T,03) = 'Anthro particulate Hg'
                  INDEX(T,03) = T + ( SPACING * 34 )
               CASE( 10 )
                  NAME (T,03) = 'KwHg'
                  FNAME(T,03) = 'Henry''s Law exchange constant'
                  UNIT (T,03) = 'cm/h'
                  INDEX(T,03) = T + ( SPACING * 34 )
                  MWT  (T,03) = 0e0
               CASE( 11 )
                  NAME (T,03) = 'HgC'
                  FNAME(T,03) = 'Hg in Colloidal phase'
                  INDEX(T,03) = T + ( SPACING * 34 )
               CASE( 12 )
                  NAME (T,03) = 'Hg_to_C'
                  FNAME(T,03) = 'Hg converted to colloidal'
                  UNIT (T,03) = 'kg/m2/s'
                  INDEX(T,03) = T + ( SPACING * 34 )
               CASE( 13 )
                  NAME (T,03) = 'Hg2_Hg0'
                  FNAME(T,03) = 'Prod of Hg2 from Hg0'
                  INDEX(T,03) = ( T - 12 ) + ( SPACING * 35 )
               CASE( 14 )
                  NAME (T,03) = 'Hg2_OH'
                  FNAME(T,03) = 'Prod of Hg2 from OH'
                  INDEX(T,03) = ( T - 12 ) + ( SPACING * 35 )
               CASE( 15 )
                  NAME (T,03) = 'Hg2_O3'
                  FNAME(T,03) = 'Prod of Hg2 from O3'
                  INDEX(T,03) = ( T - 12 ) + ( SPACING * 35 )
               CASE( 16 )
                  NAME (T,03) = 'Hg2_SS'
                  FNAME(T,03) = 'Loss of Hg2 from sea salt'
                  INDEX(T,03) = ( T - 12 ) + ( SPACING * 35 )
               CASE ( 17: )
                  NAME (T,03) = TRACER_NAME(T-16)

                  ! Tracer 3 should be "HgC" instead of "HgP"
                  IF ( TRIM( NAME(T,03) ) == 'HgP' ) THEN
                     NAME(T,03) = 'HgC'
                  ENDIF

                  FNAME(T,03) = 'Oceanic ' // TRIM( NAME(T,03) )
                  INDEX(T,03) = ( T - 16 ) + ( SPACING * 41 )
               CASE DEFAULT
                  ! Nothing
            END SELECT
         ENDDO
      ENDIF

      !-------------------------------------
      ! CO2 sources & fluxes
      !-------------------------------------
      IF ( ND04 > 0 ) THEN
         
         ! Number of tracers
         NTRAC(04) = 6

         ! Loop over tracers: HG-SRCE
         DO T = 1, NTRAC(04)

            ! Define quantities
            UNIT (T,04) = 'atoms C/cm2/s'
            INDEX(T,04) = T + ( SPACING * 40 )
            MOLC (T,04) = 1
            MWT  (T,04) = 12e-3
            SCALE(T,04) = 1e0

            ! Get name, long-name
            SELECT CASE( T )
               CASE( 1  )
                  NAME (T,04) = 'CO2ff'
                  FNAME(T,04) = 'CO2 fossil fuel emiss'
               CASE( 2  )
                  NAME (T,04) = 'CO2oc'
                  FNAME(T,04) = 'CO2 ocean emissions'
               CASE( 3  )
                  NAME (T,04) = 'CO2bal'
                  FNAME(T,04) = 'CO2 balanced biosphere'
               CASE( 4  )
                  NAME (T,04) = 'CO2bb'
                  FNAME(T,04) = 'CO2 biomass burning emiss'
               CASE( 5  )
                  NAME (T,04) = 'CO2bf'
                  FNAME(T,04) = 'CO2 biofuel emission'
               CASE( 6  )
                  NAME (T,04) = 'CO2net'
                  FNAME(T,04) = 'CO2 net terr exchange'
               CASE DEFAULT
                  ! Nothing
            END SELECT
         ENDDO
      ENDIF

      !-------------------------------------
      ! Sulfur production & loss (ND05)
      !-------------------------------------
      IF ( ND05 > 0 ) THEN

         ! Number of tracers
         NTRAC(05) = 10

         ! Loop over tracers
         DO T = 1, NTRAC(05)

            ! Get name, long-name, unit, molecular weight
            SELECT CASE( T )
               CASE( 1  )
                  NAME (T,05) = 'SO2dms'
                  FNAME(T,05) = 'P(SO2) from DMS+OH'
                  UNIT (T,05) = 'kg S'
                  MWT  (T,05) = 32e-3
               CASE( 2  )
                  NAME (T,05) = 'SO2no3'
                  FNAME(T,05) = 'P(SO2) from DMS+NO3'
                  UNIT (T,05) = 'kg S'
                  MWT  (T,05) = 32e-3
               CASE( 3  )
                  NAME (T,05) = 'SO2tot'
                  FNAME(T,05) = 'Total P(SO2)'
                  UNIT (T,05) = 'kg S'
                  MWT  (T,05) = 32e-3
               CASE( 4  )
                  NAME (T,05) = 'MSAdms'
                  FNAME(T,05) = 'P(MSA) from DMS'
                  UNIT (T,05) = 'kg S'
                  MWT  (T,05) = 32e-3
               CASE( 5  )
                  NAME (T,05) = 'SO4gas'
                  FNAME(T,05) = 'P(SO4) by gas phase'
                  UNIT (T,05) = 'kg S'
                  MWT  (T,05) = 32e-3
               CASE( 6  )
                  NAME (T,05) = 'SO4h2o2'
                  FNAME(T,05) = 'P(SO4) by cloud H2O2'
                  UNIT (T,05) = 'kg S'
                  MWT  (T,05) = 32e-3
               CASE( 7  )
                  NAME (T,05) = 'SO4o3'
                  FNAME(T,05) = 'P(SO4) by cloud O3'
                  UNIT (T,05) = 'kg S'
                  MWT  (T,05) = 32e-3
               CASE( 8  )
                  NAME (T,05) = 'SO4ss'
                  FNAME(T,05) = 'P(SO4) by seasalt O3'
                  UNIT (T,05) = 'kg S'
                  MWT  (T,05) = 32e-3
               CASE( 9  )
                  NAME (T,05) = 'LOH'
                  FNAME(T,05) = 'L(OH) by DMS'
                  UNIT (T,05) = 'kg OH'
                  MWT  (T,05) = 17e-3
               CASE( 10 )
                  NAME (T,05) = 'LNO3'
                  FNAME(T,05) = 'L(NO3) by DMS'
                  UNIT (T,05) = 'kg NO3'
                  MWT  (T,05) = 62e-3
               CASE DEFAULT
                  ! Nothing
            END SELECT

            ! Define the rest of the quantities
            INDEX(T,05) = T + ( SPACING * 18 )
            MOLC (T,05) = 1
            SCALE(T,05) = 1e0
         ENDDO
      ENDIF

      !-------------------------------------
      ! Carbon aerosol emissions (ND07)
      !-------------------------------------
      IF ( ND07 > 0 ) THEN

         ! Number of tracers
         IF ( LSOA ) THEN
            NTRAC(07) = 15
         ELSE
            NTRAC(07) = 2
         ENDIF

         ! Loop over tracers
         DO T = 1, NTRAC(07)

            ! Define the default quantities
            UNIT (T,07) = 'kg'
            MOLC (T,07) = 1
            MWT  (T,07) = 12e-3
            SCALE(T,07) = 1e0
            ! Get name, long-name, tracern umber
            SELECT CASE( T )
               CASE( 1 )
                  NAME (T,07) = 'BLKC'
                  FNAME(T,07) = 'Black (Elemental) Carbon'
                  INDEX(T,07) = IDTBCPI + ( SPACING * 33 )
               CASE( 2 )
                  NAME (T,07) = 'ORGC'
                  FNAME(T,07) = 'Organic Carbon'
                  INDEX(T,07) = IDTOCPI + ( SPACING * 33 )
               CASE( 3 )
                  NAME (T,07) = 'ALPH'
                  FNAME(T,07) = 'Alpha-Pinene'
                  MWT  (T,07) = 136e-3
                  INDEX(T,07) = IDTALPH + ( SPACING * 33 )
               CASE( 4 )
                  NAME (T,07) = 'LIMO'
                  FNAME(T,07) = 'Limonene'
                  MWT  (T,07) = 136e-3
                  INDEX(T,07) = IDTLIMO + ( SPACING * 33 )
               CASE( 5 )
                  NAME (T,07) = 'TERP'
                  FNAME(T,07) = 'Terpenes'
                  MWT  (T,07) = 136e-3
                  INDEX(T,07) = IDTLIMO + 1 + ( SPACING * 33 )
               CASE( 6 )
                  NAME (T,07) = 'ALCO'
                  FNAME(T,07) = 'Alcohols'
                  MWT  (T,07) = 142e-3     ! Actually carbon_mod.f uses 154.25
                  INDEX(T,07) = IDTLIMO + 2 + ( SPACING * 33 )
               CASE( 7 )
                  NAME (T,07) = 'SESQ'
                  FNAME(T,07) = 'Sesqterpene'
                  MWT  (T,07) = 204e-3
                  INDEX(T,07) = IDTLIMO + 3 + ( SPACING * 33 )
               CASE( 8 )
                  NAME (T,07) = 'SOA1'
                  FNAME(T,07) = 'Aer prods of ALPH+LIMO+TERP ox'
                  MWT  (T,07) = 150e-3
                  INDEX(T,07) = IDTSOA1 + ( SPACING * 33 )
               CASE( 9 )
                  NAME (T,07) = 'SOA2'
                  FNAME(T,07) = 'Aerosol prod of ALCO ox'
                  MWT  (T,07) = 160e-3
                  INDEX(T,07) = IDTSOA2 + ( SPACING * 33 )
               CASE( 10 )
                  NAME (T,07) = 'SOA3'
                  FNAME(T,07) = 'Aerosol prod of SESQ ox'
                  MWT  (T,07) = 220e-3
                  INDEX(T,07) = IDTSOA3 + ( SPACING * 33 )

               CASE( 11 )
                  NAME (T,07) = 'SOA4'
                  FNAME(T,07) = 'Aerosol prod of ISOP ox'
                  MWT  (T,07) = 130e-3
                  INDEX(T,07) = IDTSOA4 + ( SPACING * 33 )
               CASE( 12 )
                  NAME (T,07) = 'SOAG'
                  FNAME(T,07) = 'SOAG production in aq. aerosol'
                  MWT  (T,07) = 58e-3
                  INDEX(T,07) = IDTSOAG + ( SPACING * 33 )

               CASE( 13 )
                  NAME (T,07) = 'SOAM'
                  FNAME(T,07) = 'SOAM production in aq. aerosol'
                  MWT  (T,07) = 72e-3
                  INDEX(T,07) = IDTSOAM + ( SPACING * 33 )

               CASE( 14 )
                  NAME (T,07) = 'SOAG'
                  FNAME(T,07) = 'SOAG production in clouds'
                  MWT  (T,07) = 58e-3
                  INDEX(T,07) = 91 + ( SPACING * 33 )

               CASE( 15 )
                  NAME (T,07) = 'SOAM'
                  FNAME(T,07) = 'SOAM production in clouds'
                  MWT  (T,07) = 72e-3
                  INDEX(T,07) = 92 + ( SPACING * 33 )

               CASE DEFAULT
                  ! Nothing
            END SELECT

         ENDDO
      ENDIF

      !-------------------------------------
      ! HCN & CH3CN sources & sinks (ND09)
      !-------------------------------------
      IF ( ND09 > 0 ) THEN
         
         ! Number of tracers
         NTRAC(09) = N_TRACERS + 6

         ! Loop over tracers: HCN-PL-$
         DO T = 1, NTRAC(09)-6
            NAME (T,09) = TRACER_NAME(T)
            FNAME(T,09) = TRIM( TRACER_NAME(T) ) // ' lost via OH rxn'
            UNIT (T,09) = 'kg'
            INDEX(T,09) = T + ( SPACING * 38 )
            MOLC (T,09) = 1
            MWT  (T,09) = TRACER_MW_KG(T)
            SCALE(T,09) = 1e0
         ENDDO

         ! Loop over tracers: HCN-SRCE
         DO T = NTRAC(09)-6+1, NTRAC(09)

            ! Get tracer number, C number, scale factor, unit
            INDEX(T,09) = ( T - N_TRACERS ) + ( SPACING * 39 )
            MOLC (T,09) = 1
            SCALE(T,09) = 1e0
            UNIT (T,09) = 'molec/cm2/s'
              
            ! Get name, longname, unit, mwt
            IF ( T == N_TRACERS+1 ) THEN
               NAME (T,09) = 'HCNbb'
               FNAME(T,09) = 'HCN biomass emissions'
               MWT  (T,09) = 27d-3
            ELSE IF ( T == N_TRACERS+2 ) THEN
               NAME (T,09) = 'CH3CNbb'
               FNAME(T,09) = 'CH3CN biomass emissions'
               MWT  (T,09) = 41d-3
            ELSE IF ( T == N_TRACERS+3 ) THEN
               NAME (T,09) = 'HCNdf'
               FNAME(T,09) = 'HCN domestic FF emiss'
               MWT  (T,09) = 27d-3
            ELSE IF ( T == N_TRACERS+4 ) THEN
               NAME (T,09) = 'CH3CNdf'
               FNAME(T,09) = 'CH3CN domestic FF emiss'
               MWT  (T,09) = 41d-3
            ELSE IF ( T == N_TRACERS+5 ) THEN
               NAME (T,09) = 'HCNoc'
               FNAME(T,09) = 'HCN ocean uptake'
               MWT  (T,09) = 27d-3
            ELSE IF ( T == N_TRACERS+6 ) THEN
               NAME (T,09) = 'CH3CNoc'
               FNAME(T,09) = 'CH3CN ocean uptake'
               MWT  (T,09) = 41d-3
            ENDIF
         ENDDO
      ENDIF
    
      !-------------------------------------      
      ! H2-HD source, production & loss (ND10)
      !-------------------------------------      
      IF ( ND10 > 0 ) THEN
         
         ! Number of tracers
         NTRAC(10) = 20

         ! Loop over tracers 
         DO T = 1, NTRAC(10)

            ! Define quantities
            UNIT (T,10) = 'molec/cm3/s' ! overwrite below when needed
            MOLC (T,10) = 1   
            SCALE(T,10) = 1e0

            ! Get name, long-name, index, and new units
            SELECT CASE( T )
               CASE( 1  )
                  NAME (T,10) = 'H2oh'
                  FNAME(T,10) = 'Loss of H2 by OH'
                  MWT  (T,10) = 2e-3
                  INDEX(T,10) = T + ( SPACING * 47 )
               CASE( 2  )
                  NAME (T,10) = 'H2iso'
                  FNAME(T,10) = 'Prod of H2 from Isoprene'
                  MWT  (T,10) = 2e-3
                  INDEX(T,10) = T + ( SPACING * 47 )
               CASE( 3  )
                  NAME (T,10) = 'H2ch4'
                  FNAME(T,10) = 'Prod of H2 from CH4'
                  MWT  (T,10) = 2e-3
                  INDEX(T,10) = T + ( SPACING * 47 )
               CASE( 4  )
                  NAME (T,10) = 'H2ch3oh'
                  FNAME(T,10) = 'Prod of H2 from CH3OH'
                  MWT  (T,10) = 2e-3
                  INDEX(T,10) = T + ( SPACING * 47 )
               CASE( 5  )
                  NAME (T,10) = 'H2mono'
                  FNAME(T,10) = 'Prod of H2 from MONO'
                  MWT  (T,10) = 2e-3
                  INDEX(T,10) = T + ( SPACING * 47 )
               CASE( 6  )
                  NAME (T,10) = 'H2acet'
                  FNAME(T,10) = 'Prod of H2 from ACET'
                  MWT  (T,10) = 2e-3
                  INDEX(T,10) = T + ( SPACING * 47 )
               CASE( 7  )
                  NAME (T,10) = 'H2o1d'
                  FNAME(T,10) = 'Loss of H2 by strat O1D'
                  MWT  (T,10) = 2e-3
                  INDEX(T,10) = T + ( SPACING * 47 )
               CASE( 8  )
                  NAME (T,10) = 'HDoh'
                  FNAME(T,10) = 'Loss of HD by OH'
                  UNIT (T,10) = 's-1'
                  MWT  (T,10) = 3e-3
                  INDEX(T,10) = T + ( SPACING * 47 )
               CASE( 9  )
                  NAME (T,10) = 'HDiso'
                  FNAME(T,10) = 'Prod of HD from Isoprene'
                  MWT  (T,10) = 3e-3
                  INDEX(T,10) = T + ( SPACING * 47 )
               CASE( 10 )
                  NAME (T,10) = 'HDch4'
                  FNAME(T,10) = 'Prod of HD from CH4'
                  MWT  (T,10) = 3e-3
                  INDEX(T,10) = T + ( SPACING * 47 )
               CASE( 11 )
                  NAME (T,10) = 'HDch3oh'
                  FNAME(T,10) = 'Prod of HD from CH3OH'
                  MWT  (T,10) = 3e-3
                  INDEX(T,10) = T + ( SPACING * 47 )
               CASE( 12 )
                  NAME (T,10) = 'HDmono'
                  FNAME(T,10) = 'Prod of HD from MONO'
                  MWT  (T,10) = 3e-3
                  INDEX(T,10) = T + ( SPACING * 47 )
               CASE( 13 )
                  NAME (T,10) = 'HDacet'
                  FNAME(T,10) = 'Prod of HD from ACET'
                  MWT  (T,10) = 3e-3
                  INDEX(T,10) = T + ( SPACING * 47 )
               CASE( 14 )
                  NAME (T,10) = 'HDo1d'
                  FNAME(T,10) = 'Loss of HD by strat O1D'
                  MWT  (T,10) = 3e-3
                  INDEX(T,10) = T + ( SPACING * 47 )
               CASE( 15 )
                  NAME (T,10) = 'Alpha'
                  FNAME(T,10) = 'Ratio of OH k rates kHD/kH2'
                  UNIT (T,10) = 'unitless'
                  INDEX(T,10) = T + ( SPACING * 47 )
               CASE( 16 )
                  NAME (T,10) = 'H2anth'
                  FNAME(T,10) = 'Fossil Fuel H2 (anthropogenic)'
                  UNIT (T,10) = 'molec/cm2/s'
                  MWT  (T,10) = 2e-3
                  INDEX(T,10) = ( T - 15 ) + ( SPACING * 48 )
               CASE( 17 )
                  NAME (T,10) = 'H2bb'
                  FNAME(T,10) = 'Biomass Burning of H2'
                  UNIT (T,10) = 'molec/cm2/s'
                  MWT  (T,10) = 2e-3
                  INDEX(T,10) = ( T - 15 ) + ( SPACING * 48 )

               CASE( 18 )
                  NAME (T,10) = 'H2bf'
                  FNAME(T,10) = 'Biofuel Burning of H2'
                  UNIT (T,10) = 'molec/cm2/s'
                  MWT  (T,10) = 2e-3
                  INDEX(T,10) = ( T - 15 ) + ( SPACING * 48 )

               CASE( 19 )
                  NAME (T,10) = 'H2ocean'
                  FNAME(T,10) = 'Ocean emissions of H2'
                  UNIT (T,10) = 'molec/cm2/s'
                  MWT  (T,10) = 2e-3
                  INDEX(T,10) = ( T - 15 ) + ( SPACING * 48 )
               CASE( 20 )
                  NAME (T,10) = 'HDocean'
                  FNAME(T,10) = 'Ocean emissions of HD'
                  UNIT (T,10) = 'molec/cm2/s'
                  MWT  (T,10) = 3e-3
                  INDEX(T,10) = ( T - 15 ) + ( SPACING * 48 )

               CASE DEFAULT
                  ! Nothing
            END SELECT
         ENDDO
      ENDIF

      !-------------------------------------
      ! Acetone sources & sinks (ND11)
      !-------------------------------------
      IF ( ND11 > 0 ) THEN

         ! Number of tracers
         NTRAC(11) = 7

         ! Loop over tracers
         DO T = 1, NTRAC(11)

            ! Get name and long-name
            SELECT CASE( T )
               CASE( 1 )
                  NAME (T,11) = 'ACETmo'
                  FNAME(T,11) = 'ACET from monoterpenes'
               CASE( 2 )
                  NAME (T,11) = 'ACETmb'
                  FNAME(T,11) = 'ACET from methyl butenol'
               CASE( 3 )
                  NAME (T,11) = 'ACETbg'
                  FNAME(T,11) = 'ACET from direct emissions'
               CASE( 4 )
                  NAME (T,11) = 'ACETdl'
                  FNAME(T,11) = 'ACET from dry leaf matter'
               CASE( 5 )
                  NAME (T,11) = 'ACETgr'
                  FNAME(T,11) = 'ACET from grasslands'
               CASE( 6 )
                  NAME (T,11) = 'ACETop'
                  FNAME(T,11) = 'ACET from ocean source'
               CASE( 7 )
                  NAME (T,11) = 'ACETol'
                  FNAME(T,11) = 'ACET from ocean sink'
               CASE DEFAULT
                  ! Nothing
            END SELECT

            ! Define the rest of the quantities
            UNIT (T,11) = 'atoms C/cm2/s'
            INDEX(T,11) = T + ( SPACING * 22 )
            MOLC (T,11) = 3
            MWT  (T,11) = 12e-3
            SCALE(T,11) = 1e0
         ENDDO
      ENDIF

      !-------------------------------------      
      ! Distrib of emissions in PBL (ND12)
      !-------------------------------------      
      IF ( ND12 > 0 ) THEN
         T           = 1
         NTRAC(12)   = T
         NAME (T,12) = 'BLfrac'
         FNAME(T,12) = 'Boundary layer fraction'
         UNIT (T,12) = 'unitless'
         INDEX(T,12) = T + ( SPACING * 23 )
         MOLC (T,12) = 1
         MWT  (T,12) = 0e0
         SCALE(T,12) = 1e0
      ENDIF

      !-------------------------------------      
      ! Fraction of grid box experiencing
      ! LS or conv precipitation (ND16)
      !-------------------------------------   
      IF ( ND16 > 0 ) THEN

         ! Number of tracers
         NTRAC(16) = 2

         ! Loop over tracers
         DO T = 1, NTRAC(16)

            ! Get name, longname
            SELECT CASE( T )
               CASE( 1 )
                  NAME (T,16) = 'LS-f'
                  FNAME(T,16) = 'Large scale fraction'
               CASE( 2 )
                  NAME (T,16) = 'CV-f'
                  FNAME(T,16) = 'Convective fraction' 
               CASE DEFAULT
                  ! Nothing
            END SELECT

            ! Define the rest of the quantities
            UNIT (T,16) = 'unitless'
            INDEX(T,16) = T + ( SPACING * 28 )
            MOLC (T,16) = 1
            MWT  (T,16) = 0e0
            SCALE(T,16) = 1e0
         ENDDO
      ENDIF

      !-------------------------------------      
      ! F/f ratios for LS & convective
      ! rainout (ND17) & washout (ND18)
      !-------------------------------------
      IF ( ND17 + ND18 > 0 ) THEN
         
         ! Number of tracers
         NTRAC(17) = GET_WETDEP_NSOL()

         ! Loop over tracers
         DO N = 1, NTRAC(17)

            ! Tracer number for Nth wetdep species
            T = GET_WETDEP_IDWETD(N)
          
            ! Define quantities
            NAME (N,17) = 'Ff_' // TRIM( TRACER_NAME(T) )
            FNAME(N,17) = 'F/f ratio for ' // TRIM( TRACER_NAME(T) )
            UNIT (N,17) = 'unitless'
            INDEX(N,17) = T + ( SPACING * 29 )
            MWT  (N,17) = TRACER_MW_KG(T)
            MOLC (N,17) = INT( TRACER_COEFF(T,1) )
            SCALE(N,17) = 1e0
         ENDDO
      ENDIF

      !-------------------------------------      
      ! CH4 Loss due to reaction w/ OH
      !-------------------------------------
      IF ( ND19 > 0 ) THEN
         
         ! Number of tracers
         NTRAC(19) = 1
         T = 1
          
         ! Define quantities
         NAME (T,19) = 'CH4Loss'
         FNAME(T,19) = 'CH4 Loss from Reaction with OH'
         UNIT (T,19) = 'kg'
         INDEX(T,19) = T + ( SPACING * 32 )
         MWT  (T,19) = 1.6e-3
         MOLC (T,19) = 1
         SCALE(T,19) = 1e0

      ENDIF

      !-------------------------------------      
      ! Optical depths (ND21)
      !-------------------------------------
      IF ( ND21 > 0 .or. DO_TIMESERIES ) THEN

         ! Number of tracers
         NTRAC(21) = 27

         ! Loop over tracers
         DO T = 1, NTRAC(21)

            ! Define quantities
            UNIT (T,21) = 'unitless'
            INDEX(T,21) = T + ( SPACING * 14 )
            MOLC (T,21) = 1
            MWT  (T,21) = 0e0
            SCALE(T,21) = 1e0

            ! Get name long-name (and sometimes, unit)
            SELECT CASE( T )
               CASE( 1  )
                  NAME (T,21) = 'OPTD'
                  FNAME(T,21) = 'Cloud optical depth'
               CASE( 2  )
                  ! GEOS-3, GEOS-4: CLDTOT
                  NAME (T,21) = 'CLDTOT'
                  FNAME(T,21) = '3-D cloud frc'
               CASE( 4  )
                  NAME (T,21) = 'OPD'
                  FNAME(T,21) = 'Mineral dust opt depth'
               CASE( 5  )
                  NAME (T,21) = 'SD'
                  FNAME(T,21) = 'Mineral dust surface area'
                  UNIT (T,21) = 'cm2/cm3'
               CASE( 6  )
                  NAME (T,21) = 'OPSO4'
                  FNAME(T,21) = 'Sulfate opt depth (550 nm)'
               CASE( 7  )
                  NAME (T,21) = 'HGSO4'
                  FNAME(T,21) = 'Hygr growth of SO4'
               CASE( 8  )
                  NAME (T,21) = 'SSO4'
                  FNAME(T,21) = 'Sulfate surface area'
                  UNIT (T,21) = 'cm2/cm3'
               CASE( 9  )
                  NAME (T,21) = 'OPBC'
                  FNAME(T,21) = 'Black carbon opt depth (550 nm)'
               CASE( 10 )
                  NAME (T,21) = 'HGBC'
                  FNAME(T,21) = 'Hygr growth of BC'
               CASE( 11 )
                  NAME (T,21) = 'SBC'
                  FNAME(T,21) = 'BC surface area'
                  UNIT (T,21) = 'cm2/cm3'
               CASE( 12 )
                  NAME (T,21) = 'OPOC'
                  FNAME(T,21) = 'Org carbon opt depth (550 nm)'
               CASE( 13 )
                  NAME (T,21) = 'HGOC'
                  FNAME(T,21) = 'Hygr growth of OC'
               CASE( 14 )
                  NAME (T,21) = 'SOC'
                  FNAME(T,21) = 'OC surface area'
                  UNIT (T,21) = 'cm2/cm3'
               CASE( 15 )
                  NAME (T,21) = 'OPSSa'
                  FNAME(T,21) = 'Seasalt (accum) opt depth (550 nm)'
               CASE( 16 ) 
                  NAME (T,21) = 'HGSSa'
                  FNAME(T,21) = 'Hygr growth of seasalt (accum)'
               CASE( 17 )
                  NAME (T,21) = 'SSSa'
                  FNAME(T,21) = 'Seasalt (accum) surface area'
                  UNIT (T,21) = 'cm2/cm3'
               CASE( 18 )
                  NAME (T,21) = 'OPSSc'
                  FNAME(T,21) = 'Seasalt (coarse) opt depth (550 nm)'
               CASE( 19 ) 
                  NAME (T,21) = 'HGSSc'
                  FNAME(T,21) = 'Hygr growth of seasalt (coarse)'
               CASE( 20 )
                  NAME (T,21) = 'SSSc'
                  FNAME(T,21) = 'Seasalt (coarse) surface area'
                  UNIT (T,21) = 'cm2/cm3'
               CASE( 21 )
                  NAME (T,21) = 'OPD1'
                  FNAME(T,21) = 'Dust bin 1 AOD (550 nm)'
               CASE( 22 )
                  NAME (T,21) = 'OPD2'
                  FNAME(T,21) = 'Dust bin 2 AOD (550 nm)'
               CASE( 23 )
                  NAME (T,21) = 'OPD3'
                  FNAME(T,21) = 'Dust bin 3 AOD (550 nm)'
               CASE( 24 )
                  NAME (T,21) = 'OPD4'
                  FNAME(T,21) = 'Dust bin 4 AOD (550 nm)'
               CASE( 25 )
                  NAME (T,21) = 'OPD5'
                  FNAME(T,21) = 'Dust bin 5 AOD (550 nm)'
               CASE( 26 )
                  NAME (T,21) = 'OPD6'
                  FNAME(T,21) = 'Dust bin 6 AOD (550 nm)'
               CASE( 27 )
                  NAME (T,21) = 'OPD7'
                  FNAME(T,21) = 'Dust bin 7 AOD (550 nm)'
               CASE DEFAULT
                  ! Nothing
            END SELECT
         ENDDO
      ENDIF

      !-------------------------------------      
      ! Optical depths (ND22)
      !-------------------------------------
      IF ( ND22 > 0 ) THEN

         IF ( ITS_A_FULLCHEM_SIM() ) THEN 

            !------------------
            ! Full-chem run
            !------------------

            ! Number of tracers
            NTRAC(22) = 8

            ! Loop over tracers
            DO T = 1, NTRAC(22)

               ! Get name
               SELECT CASE( T )
                  CASE( 1  )
                     NAME(T,22) = 'JNO2'
                  CASE( 2  )
                     NAME(T,22) = 'JHNO3'
                  CASE( 3  )
                     NAME(T,22) = 'JH2O2'
                  CASE( 4  )
                     NAME(T,22) = 'JCH2O'
                  CASE( 5  )
                     NAME(T,22) = 'JO3'
                  CASE( 6  )
                     NAME(T,22) = 'POH'

                  CASE( 7  )
                     NAME(T,22) = 'JGLYX'
                  CASE( 8  )
                     NAME(T,22) = 'JMGLY'
                 CASE DEFAULT
                  ! Nothing
               END SELECT

               ! Define quantities
               FNAME(T,22) = TRIM( NAME(T,22) ) // ' photolysis rate'
               UNIT (T,22) = 's-1'
               INDEX(T,22) = T + ( SPACING * 13 )
               MOLC (T,22) = 1
               MWT  (T,22) = 0e0
               SCALE(T,22) = 1e0
            ENDDO
         
         ELSE IF ( ITS_A_CH3I_SIM() ) THEN

            !------------------
            ! CH3I run
            !------------------
            NTRAC(22)   = 1 
            NAME (T,22) = 'JCH3I'
            FNAME(T,22) = TRIM( NAME(T,22) ) // ' photolysis rate'
            UNIT (T,22) = 's-1'
            INDEX(T,22) = T + ( SPACING * 13 )
            MOLC (T,22) = 1
            MWT  (T,22) = 0e0
            SCALE(T,22) = 1e0            

         ELSE IF ( ITS_A_HCN_SIM() ) THEN

            !------------------
            ! HCN run
            !------------------
            NTRAC(22)   = 1 
            NAME (T,22) = 'JHCN'
            FNAME(T,22) = TRIM( NAME(T,22) ) // ' photolysis rate'
            UNIT (T,22) = 's-1'
            INDEX(T,22) = T + ( SPACING * 13 )
            MOLC (T,22) = 1
            MWT  (T,22) = 0e0
            SCALE(T,22) = 1e0    

         ENDIF
      ENDIF

      !-------------------------------------      
      ! Biomass emissions (ND28)
      !-------------------------------------
      IF ( ND28 > 0 ) THEN

         IF ( ITS_A_CO2_SIM() ) THEN

            !---------------------------
            ! CO2 simulation only
            !---------------------------

            ! Number of tracers
            NTRAC(28) = 1
            
            ! Define quantities
            NAME (T,28) = 'CO2'
            FNAME(T,28) = TRIM( NAME(T,28) ) // ' biomass'
            INDEX(T,28) = 1 + ( SPACING * 45 )
            MWT  (T,28) = 44e-3
            MOLC (T,28) = 1
            UNIT (T,28) = 'molec/cm2/s'
            SCALE(T,28) = 1e0

         ELSE

            !---------------------------
            ! Full-chemistry simulation
            !---------------------------

            ! Number of tracers
            ! Add GLYX, MGLY, GLYC, HAC, C2H2,BENZ, TOLU, XYLE, C2H4
            NTRAC(28) = 23

            ! Loop over tracers
            DO T = 1, NTRAC(28)

               ! Case statement
               SELECT CASE( T )
                  CASE( 1  )
                     NAME (T,28) = 'NOx'
                     INDEX(T,28) = 1 + ( SPACING * 45 )
                     MWT  (T,28) = 14e-3
                     MOLC (T,28) = 1
                     UNIT (T,28) = 'molec/cm2/s'
                  CASE( 2  )
                     NAME (T,28) = 'CO'
                     INDEX(T,28) = 4 + ( SPACING * 45 )
                     MWT  (T,28) = 28e-3
                     MOLC (T,28) = 1
                     UNIT (T,28) = 'molec/cm2/s'
                  CASE( 3  )
                     NAME (T,28) = 'ALK4'
                     INDEX(T,28) = 5 + ( SPACING * 45 )
                     MWT  (T,28) = 12e-3
                     MOLC (T,28) = 4
                     UNIT (T,28) = 'atoms C/cm2/s'
                  CASE( 4  )
                     NAME (T,28) = 'ACET'
                     INDEX(T,28) = 9 + ( SPACING * 45 )
                     MWT  (T,28) = 12e-3
                     MOLC (T,28) = 3
                     UNIT (T,28) = 'atoms C/cm2/s'
                  CASE( 5  )
                     NAME (T,28) = 'MEK'
                     INDEX(T,28) = 10 + ( SPACING * 45 )
                     MWT  (T,28) = 12e-3
                     MOLC (T,28) = 4
                     UNIT (T,28) = 'atoms C/cm2/s'
                  CASE( 6  )
                     NAME (T,28) = 'ALD2'
                     INDEX(T,28) = 11 + ( SPACING * 45 )
                     MWT  (T,28) = 12e-3
                     MOLC (T,28) = 3
                     UNIT (T,28) = 'atoms C/cm2/s'
                  CASE( 7  )
                     NAME (T,28) = 'PRPE'
                     INDEX(T,28) = 18 + ( SPACING * 45 )
                     MWT  (T,28) = 12e-3
                     MOLC (T,28) = 3
                     UNIT (T,28) = 'atoms C/cm2/s'
                  CASE( 8  )
                     NAME (T,28) = 'C3H8'
                     INDEX(T,28) = 19 + ( SPACING * 45 )
                     MWT  (T,28) = 12e-3
                     MOLC (T,28) = 3
                     UNIT (T,28) = 'atoms C/cm2/s'
                  CASE( 9  )
                     NAME (T,28) = 'CH2O'
                     INDEX(T,28) = 20 + ( SPACING * 45 )
                     MWT  (T,28) = 30e-3
                     MOLC (T,28) = 1
                     UNIT (T,28) = 'molec/cm2/s'
                  CASE( 10 )
                     NAME (T,28) = 'C2H6'
                     INDEX(T,28) = 21 + ( SPACING * 45 )
                     MWT  (T,28) = 12e-3
                     MOLC (T,28) = 2
                     UNIT (T,28) = 'atoms C/cm2/s'
                  CASE( 11 )
                     NAME (T,28) = 'SO2'
                     INDEX(T,28) = 26 + ( SPACING * 45 )
                     MWT  (T,28) = 32e-3
                     MOLC (T,28) = 1
                     UNIT (T,28) = 'atoms S/cm2/s'
                  CASE( 12 )
                     NAME (T,28) = 'NH3'
                     INDEX(T,28) = 30 + ( SPACING * 45 )
                     MWT  (T,28) = 17e-3
                     MOLC (T,28) = 1
                     UNIT (T,28) = 'molec/cm2/s'
                  CASE( 13 )
                     NAME (T,28) = 'BC'
                     INDEX(T,28) = 34 + ( SPACING * 45 )
                     MWT  (T,28) = 12e-3
                     MOLC (T,28) = 1
                     UNIT (T,28) = 'atoms C/cm2/s'
                  CASE( 14 )
                     NAME (T,28) = 'OC'
                     INDEX(T,28) = 35 + ( SPACING * 45 )
                     MWT  (T,28) = 12e-3
                     MOLC (T,28) = 1
                     UNIT (T,28) = 'atoms C/cm2/s'
                  CASE( 15 )
                     NAME (T,28) = 'GLYX'
                     INDEX(T,28) = 55 + ( SPACING * 45 )
                     MWT  (T,28) = 58e-3
                     MOLC (T,28) = 1
                     UNIT (T,28) = 'molec/cm2/s'
                  CASE( 16 )
                     NAME (T,28) = 'MGLY'
                     INDEX(T,28) = 56 + ( SPACING * 45 )
                     MWT  (T,28) = 72e-3
                     MOLC (T,28) = 1
                     UNIT (T,28) = 'molec/cm2/s'
                  CASE( 17 )
                     NAME (T,28) = 'BENZ'
                     INDEX(T,28) = 57 + ( SPACING * 45 )
                     MWT  (T,28) = 12e-3
                     MOLC (T,28) = 6
                     UNIT (T,28) = 'atoms C/cm2/s'
                  CASE( 18 )
                     NAME (T,28) = 'TOLU'
                     INDEX(T,28) = 58 + ( SPACING * 45 )
                     MWT  (T,28) = 12e-3
                     MOLC (T,28) = 7
                     UNIT (T,28) = 'atoms C/cm2/s'
                  CASE( 19 )
                     NAME (T,28) = 'XYLE'
                     INDEX(T,28) = 59 + ( SPACING * 45 )
                     MWT  (T,28) = 12e-3
                     MOLC (T,28) = 8
                     UNIT (T,28) = 'atoms C/cm2/s'
                  CASE( 20 )
                     NAME (T,28) = 'C2H4'
                     INDEX(T,28) = 63 + ( SPACING * 45 )
                     MWT  (T,28) = 12e-3
                     MOLC (T,28) = 2
                     UNIT (T,28) = 'atoms C/cm2/s'
                  CASE( 21 )
                     NAME (T,28) = 'C2H2'
                     INDEX(T,28) = 64 + ( SPACING * 45 )
                     MWT  (T,28) = 12e-3
                     MOLC (T,28) = 2
                     UNIT (T,28) = 'atoms C/cm2/s'
                  CASE( 22 )
                     NAME (T,28) = 'GLYC'
                     INDEX(T,28) = 66 + ( SPACING * 45 )
                     MWT  (T,28) = 60e-3
                     MOLC (T,28) = 1
                     UNIT (T,28) = 'molec/cm2/s'
                  CASE( 23 )
                     NAME (T,28) = 'HAC'
                     INDEX(T,28) = 67 + ( SPACING * 45 )
                     MWT  (T,28) = 74e-3
                     MOLC (T,28) = 1
                     UNIT (T,28) = 'molec/cm2/s'
                  CASE DEFAULT
                     ! Nothing
               END SELECT

               ! Define other quantities
               FNAME(T,28) = TRIM( NAME(T,28) ) // ' biomass'
               SCALE(T,28) = 1e0
            ENDDO
         ENDIF
      ENDIF

      !-------------------------------------      
      ! CO emissions (ND29)
      !-------------------------------------
      IF ( ND29 > 0 ) THEN

         ! Number of tracers
         NTRAC(29) = 5

         ! Loop over tracers
         DO T = 1, NTRAC(29)

            ! Get name, longname
            SELECT CASE( T )
               CASE( 1 )
                  NAME (T,29) = 'COanth'
                  FNAME(T,29) = 'Anthropogenic CO'
               CASE( 2 )
                  NAME (T,29) = 'CObb'
                  FNAME(T,29) = 'Biomass CO'
               CASE( 3 )
                  NAME (T,29) = 'CObf'
                  FNAME(T,29) = 'Biofuel CO'
               CASE( 4 )
                  NAME (T,29) = 'COmeth'
                  FNAME(T,29) = 'CO from methanol'
               CASE( 5 )
                  NAME (T,29) = 'COmono'
                  FNAME(T,29) = 'CO from monoterpenes'
               CASE DEFAULT
                  ! Nothing
            END SELECT

            ! Define the rest of the quantities
            UNIT (T,29) = 'molec/cm2/s'
            INDEX(T,29) = T + ( SPACING * 20 )
            MOLC (T,29) = 1
            MWT  (T,29) = 28e-3
            SCALE(T,29) = 1e0
         ENDDO
      ENDIF

      !-------------------------------------      
      ! Land map (ND30)
      !-------------------------------------      
      IF ( ND30 > 0 ) THEN
         T           = 1
         NTRAC(30)   = T
         NAME (T,30) = 'LWI'
         FNAME(T,30) = 'GMAO Land-water indices'
         UNIT (T,30) = 'unitless'
         INDEX(T,30) = T + ( SPACING * 15 )
         MOLC (T,30) = 1
         MWT  (T,30) = 0e0
         SCALE(T,30) = 1e0
      ENDIF

      !-------------------------------------      
      ! Surface pressure (ND31)
      !-------------------------------------      
      IF ( ND31 > 0 ) THEN
         T           = 1
         NTRAC(31)   = T
         NAME (T,31) = 'PSURF'
         FNAME(T,31) = 'Surface pressure'
         UNIT (T,31) = 'hPa'
         INDEX(T,31) = T + ( SPACING * 10 )
         MOLC (T,31) = 1
         MWT  (T,31) = 0e0
         SCALE(T,31) = 1e0
      ENDIF

      !-------------------------------------      
      ! CH3I emissions sources (ND36)
      !-------------------------------------      
      IF ( ND36 > 0 .and. ITS_A_CH3I_SIM() ) THEN

         ! Number of tracers
         NTRAC(36) = NEMANTHRO

         ! The first few are tracers
         DO T = 1, NTRAC(36)
            IF ( T <= N_TRACERS ) THEN
               NAME (T,36) = TRACER_NAME(T)
               FNAME(T,36) = TRIM( NAME(T,36) ) // ' tracer'
            ELSE IF ( T == 6 ) THEN
               NAME (T,36) = 'CH3Iof'
               FNAME(T,36) = 'CH3I ocean flux??'
            ELSE IF ( T == 7 ) THEN
               NAME (T,36) = 'CH3Igf'
               FNAME(T,36) = 'CH3I gas flux??'
            ELSE IF ( T == 8 ) THEN
               NAME (T,36) = 'CH3Itf'
               FNAME(T,36) = 'CH3I total flux??'
            ENDIF

            ! Define other quantities
            MOLC (T,36) = 1
            SCALE(T,36) = 1e0
            INDEX(T,36) = T + ( SPACING * 44 )
            MWT  (T,36) = 141.9e-3
            UNIT (T,36) = 'ng/m2/s'
         ENDDO
      ENDIF

      !-------------------------------------      
      ! Afternoon-average boundary 
      ! layer heights (ND41) + timeseries
      !-------------------------------------      
      IF ( ND41 > 0 .or. DO_TIMESERIES ) THEN

         ! Number of tracers
         NTRAC(41) = 2

         ! Loop over tracers
         DO T = 1, NTRAC(41)

            ! Get name and unit for each met field
            SELECT CASE( T )
               CASE( 1 )
                  NAME(T,41) = 'PBL-M'
                  UNIT(T,41) = 'm'
               CASE( 2 )
                  NAME(T,41) = 'PBL-L'
                  UNIT(T,41) = 'levels'
               CASE DEFAULT
                  ! Nothing
            END SELECT

            ! Define the rest of the quantities
            FNAME(T,41) = 'PBL depth'
            INDEX(T,41) = T + ( SPACING * 27 )
            MOLC (T,41) = 1
            MWT  (T,41) = 0e0
            SCALE(T,41) = 1e0
         ENDDO
      ENDIF

      !-------------------------------------      
      ! SOA concentrations (ND42)
      !-------------------------------------    
      IF ( ND42 > 0 ) THEN

         ! Number of tracers
         NTRAC(42) = 8

         ! Loop over tracers
         DO T = 1, NTRAC(42)

            ! Get name, long name, unit, and mol wt for each field
            SELECT CASE( T )
               CASE( 1 )
                  NAME (T,42) = 'SOA1'
                  FNAME(T,42) = 'Aer prods of ALPH+LIMO+TERP ox'
                  UNIT (T,42) = 'ug/m3'
                  MWT  (T,42) = 150e-3
               CASE( 2 )
                  NAME (T,42) = 'SOA2'
                  FNAME(T,42) = 'Aerosol prod of ALCO ox'
                  UNIT (T,42) = 'ug/m3'
                  MWT  (T,42) = 160e-3
               CASE( 3 )
                  NAME (T,42) = 'SOA3'
                  FNAME(T,42) = 'Aerosol prod of SESQ ox'
                  UNIT (T,42) = 'ug/m3'
                  MWT  (T,42) = 220e-3
               CASE( 4 )
                  NAME (T,42) = 'SOA4'
                  FNAME(T,42) = 'Aerosol prod of ISOP ox'
                  UNIT (T,42) = 'ug/m3'
                  MWT  (T,42) = 130e-3
               CASE( 5 )
                  NAME (T,42) = 'SOA1-3'
                  FNAME(T,42) = 'SOA1 + SOA2 + SOA3'
                  UNIT (T,42) = 'ug/m3'
                  MWT  (T,42) = 177d-3
               CASE( 6 )
                  NAME (T,42) = 'SOA1-4'
                  FNAME(T,42) = 'SOA1 + SOA2 + SOA3 + SOA4'
                  UNIT (T,42) = 'ug/m3'
                  MWT  (T,42) = 165e-3
               CASE( 7 )
                  NAME (T,42) = 'sumOC'
                  FNAME(T,42) = 'Sum of organic carbon'
                  UNIT (T,42) = 'ug C/m3'
                  MWT  (T,42) = 12e-3
               CASE( 8 )
                  NAME (T,42) = 'sumOCstp'
                  FNAME(T,42) = 'Sum of organic carbon @ STP'
                  UNIT (T,42) = 'ug C/sm3'
                  MWT  (T,42) = 12e-3
               CASE DEFAULT
                  ! Nothing
            END SELECT

            ! Define the rest of the quantities
            INDEX(T,42) = T + ( SPACING * 43 )
            MOLC (T,42) = 1
            SCALE(T,42) = 1e0
         ENDDO
      ENDIF
      
      !-------------------------------------      
      ! Chemically-produced OH, etc (ND43)
      !-------------------------------------      
      IF ( ND43 > 0 ) THEN 

         ! Number of tracers
         NTRAC(43) = 5

         ! Loop over tracers
         DO T = 1, NTRAC(43)

            ! Get name and unit for each met field
            SELECT CASE( T )
               CASE( 1 )
                  NAME(T,43) = 'OH'
                  UNIT(T,43) = 'molec/cm3'
                  MWT (T,43) = 17e-3
               CASE( 2 )
                  NAME(T,43) = 'NO'
                  UNIT(T,43) = 'v/v'
                  MWT (T,43) = 30e-3
               CASE( 3 )
                  NAME(T,43) = 'HO2'
                  UNIT(T,43) = 'v/v'
                  MWT (T,43) = 33e-3
               CASE( 4 )
                  NAME(T,43) = 'NO2'
                  UNIT(T,43) = 'v/v'
                  MWT (T,43) = 46e-3
               CASE( 5 )
                  NAME(T,43) = 'NO3'
                  UNIT(T,43) = 'v/v'
                  MWT (T,43) = 62e-3
               CASE DEFAULT
                  ! Nothing
            END SELECT

            ! Define the rest of the quantities
            FNAME(T,43) = 'Chemically produced ' // TRIM( NAME(T,43) )
            INDEX(T,43) = T + ( SPACING * 16 )
            MOLC (T,43) = 1
            SCALE(T,43) = 1e0
         ENDDO
      ENDIF

      !-------------------------------------      
      ! Dry deposition fluxes     (ND44)
      ! Dry deposition velocities (ND44)
      !-------------------------------------      
      IF ( ND44 > 0 ) THEN

         ! Special handling for tagged simulations
         IF ( ITS_A_TAGOX_SIM()   .or. 
     &        ITS_A_MERCURY_SIM() .or. ITS_A_H2HD_SIM()  ) THEN

            !----------------------------------------------------
            ! Tagged runs: Save drydep flux for all tracers
            ! but drydep velocity for only NUMDEP species
            !----------------------------------------------------

            ! Drydep flux (all tracers)
            DO T = 1, N_TRACERS
               NAME (T,44)  = TRIM( TRACER_NAME(T) ) // 'df'        
               FNAME(T,44)  = TRIM( TRACER_NAME(T) ) // ' drydep flux'
               UNIT (T,44)  = 'molec/cm2/s'
               MWT  (T,44)  = TRACER_MW_KG(T) 
               MOLC (T,44)  = INT( TRACER_COEFF(T,1) )
               SCALE(T,44)  = 1.0e0
               INDEX(T,44)  = T + ( SPACING * 36 )
               NTRAC(44)    = NTRAC(44) + 1
            ENDDO

            ! Drydep velocity (only deposited species)
            DO N = 1, NUMDEP
               T            = NTRAIND(N)
               NN           = N + N_TRACERS
               NAME (NN,44) = TRIM( DEPNAME(N) ) // 'dv'        
               FNAME(NN,44) = TRIM( DEPNAME(N) ) // ' drydep velocity'
               UNIT (NN,44) = 'cm/s'
               MWT  (NN,44) = TRACER_MW_KG(T) 
               MOLC (NN,44) = INT( TRACER_COEFF(T,1) )
               SCALE(NN,44) = 1.0e0
               INDEX(NN,44) = T + ( SPACING * 37 )
               NTRAC(44)    = NTRAC(44) + 1
            ENDDO

         ELSE

            !----------------------------------------------------
            ! Otherwise save drydep flux for NUMDEP species only
            !----------------------------------------------------

            ! Loop over drydep species
            DO N = 1, NUMDEP
               
               ! GEOS-CHEM tracer #
               T            = NTRAIND(N)         

               ! Drydep flux (deposited species only)
               NAME (N,44)  = TRIM( DEPNAME(N) ) // 'df'        
               FNAME(N,44)  = TRIM( DEPNAME(N) ) // ' drydep flux'
               UNIT (N,44)  = 'molec/cm2/s'
               MWT  (N,44)  = TRACER_MW_KG(T) 
               MOLC (N,44)  = INT( TRACER_COEFF(T,1) )
               SCALE(N,44)  = 1.0e0
               INDEX(N,44)  = T + ( SPACING * 36 )
               NTRAC(44)    = NTRAC(44) + 1

               ! For the Rn simulation, unit is kg/s (bmy, 2/22/08)
               IF ( ITS_A_RnPbBe_SIM() ) UNIT(N,44) = 'kg/s'

               ! Drydep velocity (deposited species only)
               NN           = N + NUMDEP
               NAME (NN,44) = TRIM( DEPNAME(N) ) // 'dv'        
               FNAME(NN,44) = TRIM( DEPNAME(N) ) // ' drydep velocity'
               UNIT (NN,44) = 'cm/s'
               MWT  (NN,44) = TRACER_MW_KG(T) 
               MOLC (NN,44) = INT( TRACER_COEFF(T,1) )
               SCALE(NN,44) = 1.0e0
               INDEX(NN,44) = T + ( SPACING * 37 )
               NTRAC(44)    = NTRAC(44) + 1
            ENDDO
         ENDIF
      ENDIF

      !-------------------------------------      
      ! Biogenic emissions (ND46)
      !-------------------------------------      
      IF ( ND46 > 0 ) THEN 

         ! Number of tracers
         NTRAC(46) = 13 ! was 6 (mpb,2009)

         ! Loop over tracers
         DO T = 1, NTRAC(46)

            ! Get name and unit for each met field
            SELECT CASE( T )
               CASE( 1 )
                  NAME(T,46) = 'ISOP'
                  MOLC(T,46) = 5
               CASE( 2 )
                  NAME(T,46) = 'ACET'
                  MOLC(T,46) = 3
               CASE( 3 )
                  NAME(T,46) = 'PRPE'
                  MOLC(T,46) = 3
               CASE( 4 )
                  NAME(T,46) = 'MONOT'
                  MOLC(T,46) = 10
               CASE( 5 )
                  NAME(T,46) = 'MBO'
                  MOLC(T,46) = 5
               CASE( 6 )
                  NAME(T,46) = 'C2H4'
                  MOLC(T,46) = 2
               CASE( 7 )
                  NAME(T,46) = 'APINE'
                  MOLC(T,46) = 10
               CASE( 8 )
                  NAME(T,46) = 'BPINE'
                  MOLC(T,46) = 10
               CASE( 9 )
                  NAME(T,46) = 'LIMON'
                  MOLC(T,46) = 10
               CASE( 10)
                  NAME(T,46) = 'SABIN'
                  MOLC(T,46) = 10
               CASE( 11)
                  NAME(T,46) = 'MYRCN'
                  MOLC(T,46) = 10
               CASE( 12)
                  NAME(T,46) = 'CAREN'
                  MOLC(T,46) = 10
               CASE( 13)
                  NAME(T,46) = 'OCIMN'
                  MOLC(T,46) = 10
               CASE DEFAULT
                  ! Nothing
            END SELECT

            ! Define the rest of the quantities
            INDEX(T,46) = T + ( SPACING * 21 )
            FNAME(T,46) = TRIM( NAME(T,46) ) // ' emissions'
            UNIT (T,46) = 'atoms C/cm2/s'
            MWT  (T,46) = 12e-3
            SCALE(T,46) = 1e0
         ENDDO
      ENDIF

      !-------------------------------------      
      ! Time series diags (ND{48,49,50,51})
      !-------------------------------------      
      IF ( DO_TIMESERIES ) THEN 

         ! Number of tracers
         ! Increased from 26 to 32 (mpb,2009)
         NTRAC(48) = 32

         ! Loop over tracers
         DO T = 1, NTRAC(48)

            ! Define quantities
            INDEX(T,48) = T + ( SPACING * 19 )
            MOLC (T,48) = 1
            MWT  (T,48) = 0e0
            SCALE(T,48) = 1e0

            ! Get name, long-name (and others where necessary)
            SELECT CASE( T )
               CASE( 1  )
                  NAME (T,48) = 'O3'
                  FNAME(T,48) = 'O3 time series'
                  UNIT (T,48) = 'ppbv'
                  MWT  (T,48) = 48e-3
                  SCALE(T,48) = 1.0e+9
               CASE( 2  )
                  NAME (T,48) = 'OH'
                  FNAME(T,48) = 'OH time series'
                  UNIT (T,48) = 'molec/cm3'
               CASE( 3  )
                  NAME (T,48) = 'NOy'
                  FNAME(T,48) = 'NOy time series'
                  UNIT (T,48) = 'ppbv'
                  SCALE(T,48) = 1.0e+9
               CASE( 4  )
                  NAME (T,48) = 'NO2dv'
                  FNAME(T,48) = 'NO2 drydep vel'
                  UNIT (T,48) = 'cm/s'
               CASE( 5  )
                  NAME (T,48) = 'O3dv'
                  FNAME(T,48) = 'O3 drydep vel'
                  UNIT (T,48) = 'cm/s'
               CASE( 6  )
                  NAME (T,48) = 'PANdv'
                  FNAME(T,48) = 'PAN drydep vel'
                  UNIT (T,48) = 'cm/s'
               CASE( 7  )
                  NAME (T,48) = 'HNO3dv'
                  FNAME(T,48) = 'HNO3 drydep vel'
                  UNIT (T,48) = 'cm/s'
               CASE( 8  )
                  NAME (T,48) = 'H2O2dv'
                  FNAME(T,48) = 'H2O2 drydep vel'
                  UNIT (T,48) = 'cm/s'
               CASE( 9  )
                  NAME (T,48) = 'NO'
                  FNAME(T,48) = 'NO time series'
                  UNIT (T,48) = 'ppbv'
                  MWT  (T,48) = 30e-3
                  SCALE(T,48) = 1.0e+9
               CASE( 10  )
                  NAME (T,48) = 'PBLTOP'
                  FNAME(T,48) = 'PBL top pressure'
                  UNIT (T,48) = 'hPa'
               CASE( 11  )
                  NAME (T,48) = 'LT'
                  FNAME(T,48) = 'Station local time'
                  UNIT (T,48) = 'h'
               CASE( 12 )
                  NAME (T,48) = 'ISOPf'
                  FNAME(T,48) = 'ISOP emission flux'
                  UNIT (T,48) = 'molec/cm2/s'
                  MWT  (T,48) = 12e-3
               CASE( 13 )
                  NAME (T,48) = 'MOLEN'
                  FNAME(T,48) = 'Monin-Obhukov length'
                  UNIT (T,48) = 'm'
               CASE( 14 )
                  NAME (T,48) = 'NO3'
                  FNAME(T,48) = 'NO3 time series'
                  UNIT (T,48) = 'ppbv'
                  MWT  (T,48) = 62e-3
                  SCALE(T,48) = 1.0e+9
               CASE( 15 )
                  NAME (T,48) = 'PS'
                  FNAME(T,48) = 'Surface pressure'
                  UNIT (T,48) = 'hPa'
               CASE( 16 )
                  NAME (T,48) = 'AVGW'
                  FNAME(T,48) = 'Water vapor mixing ratio'
                  UNIT (T,48) = 'ppbv'
                  MWT  (T,48) = 33e-3
                  SCALE(T,48) = 1.0e+9
               CASE( 17 )
                  NAME (T,48) = 'RH'
                  FNAME(T,48) = 'Relative humidity'
                  UNIT (T,48) = '%'
               CASE( 18 )
                  NAME (T,48) = 'HCNcol'
                  FNAME(T,48) = 'HCN column density'
                  UNIT (T,48) = 'molec/cm2'
               CASE( 19 )
                  NAME (T,48) = 'CF'
                  FNAME(T,48) = 'Cloud fraction'
                  UNIT (T,48) = 'unitless'
               CASE( 20 )
                  NAME (T,48) = 'ODCOL'
                  FNAME(T,48) = 'Cloud optical depth'
                  UNIT (T,48) = 'unitless'
               CASE( 21 )
                  NAME (T,48) = 'CThgt'
                  FNAME(T,48) = 'Cloud top height'
                  UNIT (T,48) = 'hPa'
               CASE( 22 )
                  NAME (T,48) = 'AIRDEN'
                  FNAME(T,48) = 'Air density'
                  UNIT (T,48) = 'molec/cm3'
               CASE( 23 )
                  NAME (T,48) = 'TOTDUST'
                  FNAME(T,48) = 'Total dust opt depth'
                  UNIT (T,48) = 'unitless'
               CASE( 24 )
                  NAME (T,48) = 'TOTSALT'
                  FNAME(T,48) = 'Total seasalt tracer'
                  UNIT (T,48) = 'ppbv'
                  MWT  (T,48) = 36e-3
                  SCALE(T,48) = 1.0e+9
               CASE( 25 )
                  NAME (T,48) = 'NO2'
                  FNAME(T,48) = 'NO2 concentration'
                  UNIT (T,48) = 'ppbv'
                  MWT  (T,48) = 46e-3
                  SCALE(T,48) = 1.0e+9                  
               CASE( 26 )
                  NAME (T,48) = 'TOTAOD'
                  FNAME(T,48) = 'Total aerosol opt depth'
                  UNIT (T,48) = 'unitless'
                  ! +++++++++++++++++++++++++++++++++++++
                  ! Add in gamma factors (mpb,2008)
                  ! +++++++++++++++++++++++++++++++++++++
               CASE( 27 )
                  NAME (T,48) = 'G_T'
                  FNAME(T,48) = 'Temperature activity factor'
                  UNIT (T,48) = 'unitless'
               CASE( 28 )
                  NAME (T,48) = 'G_PAR'
                  FNAME(T,48) = 'Light activity factor'
                  UNIT (T,48) = 'unitless'
               CASE( 29 )              
                  NAME (T,48) = 'G_AGE'
                  FNAME(T,48) = 'Leaf-age activity factor'
                  UNIT (T,48) = 'unitless'
               CASE( 30 )
                  NAME (T,48) = 'G_LAI'
                  FNAME(T,48) = 'LAI activity factor'
                  UNIT (T,48) = 'unitless'
               CASE( 31 )
                  NAME (T,48) = 'G_SM'
                  FNAME(T,48) = 'Soil-moisture activity factor'
                  UNIT (T,48) = 'unitless'
               CASE( 32 )
                  NAME (T,48) = 'D_LAI'
                  FNAME(T,48) = 'Daily LAI'
                  UNIT (T,48) = 'm2/m2'
               CASE DEFAULT
                  ! Nothing
            END SELECT

         ENDDO
      ENDIF


      !-----------------------------------
      ! gamma HO2 (ND52) jaegle (02/26/09)
      !-----------------------------------
      IF ( ND52 > 0 ) THEN
         T           = 1
         NTRAC(52)   = T
         NAME (T,52) = 'GAMMAHO2'
         FNAME(T,52) = 'Gamma HO2'
         UNIT (T,52) = 'unitless'
         INDEX(T,52) = T + ( SPACING * 49 )
         MOLC (T,52) = 1
         MWT  (T,52) = 0e0
         SCALE(T,52) = 1e0
      ENDIF

      !-------------------------------------      
      ! Time in the troposphere (ND54)
      !-------------------------------------      
      IF ( ND54 > 0 ) THEN
         T           = 1
         NTRAC(54)   = T
         NAME (T,54) = 'TIMETROP'
         FNAME(T,54) = 'Time in the troposphere'
         UNIT (T,54) = 'unitless'
         INDEX(T,54) = T + ( SPACING * 46 )
         MOLC (T,54) = 1
         MWT  (T,54) = 0e0
         SCALE(T,54) = 1e0
      ENDIF

      !-------------------------------------      
      ! Tropopause quantities (ND55)
      !-------------------------------------      
      IF ( ND55 > 0 ) THEN 

         ! Number of tracers
         NTRAC(55) = 3

         ! Loop over tracers
         DO T = 1, NTRAC(55)

            ! Get name and unit for each met field
            SELECT CASE( T )
               CASE( 1 )
                  NAME (T,55) = 'TP-LEVEL'
                  FNAME(T,55) = 'Tropopause level'
                  UNIT (T,55) = 'level'
               CASE( 2 )
                  NAME (T,55) = 'TP-HGHT'
                  FNAME(T,55) = 'Tropopause height'
                  UNIT (T,55) = 'km'
               CASE( 3 )
                  NAME (T,55) = 'TP-PRESS'
                  FNAME(T,55) = 'Tropopause pressure'
                  UNIT (T,55) = 'hPa'
               CASE DEFAULT
                  ! Nothing
            END SELECT

            ! Define the rest of the quantities
            INDEX(T,55) = T + ( SPACING * 26 )
            MOLC (T,55) = 1
            MWT  (T,55) = 0e0
            SCALE(T,55) = 1e0
         ENDDO
      ENDIF

      !-------------------------------------      
      ! Lightning flash rates (ND56)
      !-------------------------------------      
      IF ( ND56 > 0 ) THEN 

         ! Number of tracers
         NTRAC(56) = 3

         ! Loop over tracers
         DO T = 1, NTRAC(56)

            ! Get name and long name for each met field
            SELECT CASE( T )
               CASE( 1 )
                  NAME (T,56) = 'L-RATE'
                  FNAME(T,56) = 'Total lightning flash rate'
               CASE( 2 )
                  NAME (T,56) = 'IC-RATE'
                  FNAME(T,56) = 'Intra-cloud flash rate'
               CASE( 3 )
                  NAME (T,56) = 'CG-RATE'
                  FNAME(T,56) = 'Cloud-ground flash rate'
               CASE DEFAULT
                  ! Nothing
            END SELECT

            ! Define the rest of the quantities
            UNIT (T,56) = 'flashes/min/km2'
            INDEX(T,56) = T + ( SPACING * 42 )
            MOLC (T,56) = 1
            MWT  (T,56) = 0e0
            SCALE(T,56) = 1e0
         ENDDO
      ENDIF

      !-------------------------------------
      ! CH4 Emissions
      !-------------------------------------
      IF ( ND58 > 0 ) THEN

         ! Number of tracers
         NTRAC(58) = 12

         ! Loop over tracers
         DO T = 1, NTRAC(58)

            ! Get name and long name for each tracer
            SELECT CASE( T )
               CASE( 1 )
                  NAME (T,58) = 'CH4-TOT'
                  FNAME(T,58) = 'CH4 Emissions total (w/o sab)'
               CASE( 2 )
                  NAME (T,58) = 'CH4-GAO'
                  FNAME(T,58) = 'CH4 Emissions gas & oil'
               CASE( 3 )
                  NAME (T,58) = 'CH4-COL'
                  FNAME(T,58) = 'CH4 Emissions coal'
               CASE( 4 )
                  NAME (T,58) = 'CH4-LIV'
                  FNAME(T,58) = 'CH4 Emissions livestock'
               CASE( 5 )
                  NAME (T,58) = 'CH4-WST'
                  FNAME(T,58) = 'CH4 Emissions waste'
               CASE( 6 )
                  NAME (T,58) = 'CH4-BFL'
                  FNAME(T,58) = 'CH4 Emissions biofuel'
               CASE( 7 )
                  NAME (T,58) = 'CH4-RIC'
                  FNAME(T,58) = 'CH4 Emissions rice'
               CASE( 8 )
                  NAME (T,58) = 'CH4-OTA'
                  FNAME(T,58) = 'CH4 Emissions other anthro'
               CASE( 9 )
                  NAME (T,58) = 'CH4-BBN'
                  FNAME(T,58) = 'CH4 Emissions bioburn'
               CASE( 10 )
                  NAME (T,58) = 'CH4-WTL'
                  FNAME(T,58) = 'CH4 Emissions wetlands'
               CASE( 11 )
                  NAME (T,58) = 'CH4-SAB'
                  FNAME(T,58) = 'CH4 Emissions soil abs'
               CASE( 12 )
                  NAME (T,58) = 'CH4-OTN'
                  FNAME(T,58) = 'CH4 Emissions other natural'
               CASE( 13 )
                  NAME (T,58) = 'CH4-Loss'
                  FNAME(T,58) = 'CH4 Loss by Reaction w/ OH'
               CASE DEFAULT
                  ! Nothing
            END SELECT

            ! Define the rest of the quantities
            UNIT (T,58) = 'kg'
            INDEX(T,58) = T + ( SPACING * 50 )
            MOLC (T,58) = 1
            MWT  (T,58) = 1.6e-3
            SCALE(T,58) = 1e0
         ENDDO
      ENDIF

      !-------------------------------------      
      ! Wetland Fraction (ND60)
      !-------------------------------------      
      IF ( ND60 > 0 ) THEN

         ! Number of tracers
         NTRAC(60) = 1
         T         = 1

         NAME (T,60) = 'WETFRAC'
         FNAME(T,60) = 'Wetland Fraction'
         INDEX(T,60) = T + ( SPACING * 51 )
         MOLC (T,60) = 0 
         MWT  (T,60) = 0e0
         SCALE(T,60) = 1e0
         UNIT(T,60)  = 'unitless'
      ENDIF

      !-------------------------------------      
      ! Family prod & loss (ND65)
      !-------------------------------------      
      IF ( DO_SAVE_PL ) THEN

         ! Number of P/L families
         NTRAC(65) = GET_NFAM()

         ! Loop over each P/L family
         DO T = 1, NTRAC(65)
            NAME (T,65) = GET_FAM_NAME( T ) 
            FNAME(T,65) = TRIM( NAME(T,65) ) // ' P/L family'
            INDEX(T,65) = T + ( SPACING * 17 )
            MOLC (T,65) = 1 
            MWT  (T,65) = GET_FAM_MWT( T )
            SCALE(T,65) = 1e0
            
            ! Unit for Tag Ox is kg/s, otherwise molec/cm3/s 
            IF ( ITS_A_TAGOX_SIM() ) THEN
               UNIT(T,65) = 'kg/s'
            ELSE
               UNIT(T,65) = 'molec/cm3/s'
            ENDIF
         ENDDO
      ENDIF

      !-------------------------------------      
      ! 3-D GMAO met fields (ND66)
      ! also for timeseries diagnostics
      !-------------------------------------      
      IF ( ND66 > 0 .or. DO_TIMESERIES ) THEN

         ! Number of tracers
         NTRAC(66) = 6

         ! Loop over tracers
         DO T = 1, NTRAC(66)

            ! Get name and unit for each met field
            SELECT CASE( T )
               CASE( 1  )
                  NAME(T,66) = 'UWND'
                  UNIT(T,66) = 'm/s'
               CASE( 2  )
                  NAME(T,66) = 'VWND'
                  UNIT(T,66) = 'm/s'
               CASE( 3  )
                  NAME(T,66) = 'TMPU'
                  UNIT(T,66) = 'K'
               CASE( 4 )
                  NAME(T,66) = 'SPHU'
                  UNIT(T,66) = 'g/kg'
               CASE( 5 )
#if   defined( GEOS_4 )
                  NAME(T,66) = 'ZMMU'
                  UNIT(T,66) = 'Pa/s'
#elif defined( GEOS_5 )
                  NAME(T,66) = 'CMFMC'
                  UNIT(T,66) = 'kg/m2/s'
#else
                  NAME(T,66) = 'CLDMAS'
                  UNIT(T,66) = 'kg/m2/s'
#endif
               CASE( 6 )
                  NAME(T,66) = 'DTRAIN'
                  UNIT(T,66) = 'kg/m2/s'
               CASE DEFAULT
                  ! Nothing
            END SELECT

            ! Define the rest of the quantities
            FNAME(T,66) = 'GMAO ' // TRIM( NAME(T,66) ) // ' field'
            INDEX(T,66) = T + ( SPACING * 12 )
            MOLC (T,66) = 1
            MWT  (T,66) = 0e0
            SCALE(T,66) = 1e0
         ENDDO
      ENDIF

      !-------------------------------------      
      ! 2-D GMAO met fields (ND67)
      !-------------------------------------      
      IF ( ND67 > 0 .or. DO_TIMESERIES ) THEN 

         ! Number of tracers
         NTRAC(67) = 23 ! (Lin, 03/31/09)

         ! Loop over tracers
         DO T = 1, NTRAC(67)

            ! Get name and unit for each met field
            SELECT CASE( T )
               CASE( 1  )
                  NAME(T,67) = 'HFLUX'
                  UNIT(T,67) = 'W/m2'
               CASE( 2  )
                  NAME(T,67) = 'RADSWG'
                  UNIT(T,67) = 'W/m2'
               CASE( 3  )
                  NAME(T,67) = 'PREACC'
                  UNIT(T,67) = 'W/m2'
               CASE( 4  )
                  NAME(T,67) = 'PRECON'
                  UNIT(T,67) = 'W/m2'
               CASE( 5  )
                  NAME(T,67) = 'TS'
                  UNIT(T,67) = 'K'
               CASE( 6  )
                  NAME(T,67) = 'RADSWT'
                  UNIT(T,67) = 'W/m2'
               CASE( 7  )
                  NAME(T,67) = 'USTAR'
                  UNIT(T,67) = 'm/s'
               CASE( 8  )
                  NAME(T,67) = 'Z0'
                  UNIT(T,67) = 'm'
               CASE( 9  )
                  NAME(T,67) = 'PBL'
#if   defined( GEOS_4 )
                  UNIT(T,67) = 'm'
#else
                  UNIT(T,67) = 'hPa'
#endif
               CASE( 10 )
                  NAME(T,67) = 'CLDFRC'
                  UNIT(T,67) = 'unitless'
               CASE( 11 )
                  NAME(T,67) = 'U10M'
                  UNIT(T,67) = 'm/s'
               CASE( 12 )
                  NAME(T,67) = 'V10M'
                  UNIT(T,67) = 'm/s'
               CASE( 13 )
                  NAME(T,67) = 'PS-PBL'
                  UNIT(T,67) = 'hPa'
               CASE( 14 )
                  NAME(T,67) = 'ALBD'
                  UNIT(T,67) = 'unitless'
               CASE( 15 )
                  NAME(T,67) = 'PHIS'
                  UNIT(T,67) = 'm'
               CASE( 16 )
                  NAME(T,67) = 'CLDTOP'
                  UNIT(T,67) = 'level'
               CASE( 17 )
                  NAME(T,67) = 'TROPP'
                  UNIT(T,67) = 'hPa'
               CASE( 18 )
                  NAME(T,67) = 'SLP'
                  UNIT(T,67) = 'hPa'
               CASE( 19 )
                  NAME(T,67) = 'TSKIN'
                  UNIT(T,67) = 'K'
               CASE( 20 )
                  NAME(T,67) = 'PARDF'
                  UNIT(T,67) = 'W/m2'
               CASE( 21 )
                  NAME(T,67) = 'PARDR'
                  UNIT(T,67) = 'W/m2'
               CASE( 22 )
                  NAME(T,67) = 'GWET'
                  UNIT(T,67) = 'unitless'
               CASE( 23 )
                  ! Add EFLUX (Lin, 05/16/08)
                  NAME(T,67) = 'EFLUX'
                  UNIT(T,67) = 'W/m2'
               CASE DEFAULT
                  ! Nothing
            END SELECT

            ! Define the rest of the quantities
            FNAME(T,67) = 'GMAO ' // TRIM( NAME(T,67) ) // ' field'
            INDEX(T,67) = T + ( SPACING * 11 )
            MOLC (T,67) = 1
            MWT  (T,67) = 0e0
            SCALE(T,67) = 1e0
         ENDDO
      ENDIF

      !-------------------------------------      
      ! Grid box heights and related 
      ! quantities (ND68) + timeseries
      !-------------------------------------      
      IF ( ND68 > 0 .or. DO_TIMESERIES ) THEN

         ! Number of tracers
         NTRAC(68) = 4

         ! Loop over tracers
         DO T = 1, NTRAC(68)

            ! Get name and unit for each met field
            SELECT CASE( T )
               CASE( 1 )
                  NAME (T,68) = 'BXHEIGHT'
                  FNAME(T,68) = 'Grid box height'
                  UNIT (T,68) = 'm'
               CASE( 2 )
                  NAME (T,68) = 'AD'
                  FNAME(T,68) = 'Air mass in grid box'
                  UNIT (T,68) = 'kg'
               CASE( 3 )
                  NAME (T,68) = 'AVGW'
                  FNAME(T,68) = 'Mixing ratio of H2O vapor'
                  UNIT (T,68) = 'v/v'
               CASE( 4 )
                  NAME (T,68) = 'N(AIR)'
                  FNAME(T,68) = 'Number density of air'
                  UNIT (T,68) = 'molec/m3'
               CASE DEFAULT
                  ! Nothing
            END SELECT

            ! Define the rest of the quantities
            INDEX(T,68) = T + ( SPACING * 24 )
            MOLC (T,68) = 1
            MWT  (T,68) = 0e0
            SCALE(T,68) = 1e0
         ENDDO
      ENDIF

      !-------------------------------------      
      ! Grid box surface area (ND69)
      !-------------------------------------      
      IF ( ND69 > 0 ) THEN 
         T           = 1
         NTRAC(69)   = T
         NAME (T,69) = 'DXYP'
         FNAME(T,69) = 'Grid box surface area'
         UNIT (T,69) = 'm2'
         INDEX(T,69) = T + ( SPACING * 25 )
         MOLC (T,69) = 1
         MWT  (T,69) = 0e0
         SCALE(T,69) = 1e0
      ENDIF

      ! Return to calling program
      END SUBROUTINE INIT_TRACERINFO

!------------------------------------------------------------------------------

      SUBROUTINE INIT_GAMAP( DIAGINFO, TRACERINFO )
!
!******************************************************************************
!  Subroutine INIT_GAMAP allocates and initializes all module variables.  
!  (bmy, 4/22/05, 8/4/06)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) DIAGINFO   (CHARACTER) : Path name of the GAMAP "diaginfo.dat"   file
!  (2 ) TRACERINFO (CHARACTER) : Path name of the GAMAP "tracerinfo.dat" file
!
!  NOTES:
!  (1 ) Now add proper UNIT & SCALE for Rn/Pb/Be simulations (bmy, 5/11/05)
!  (2 ) Added HCN & CH3CN source & sink info for ND09 (bmy, 6/27/05)
!  (3 ) Bug fix: removed duplicate category names.  Updated for CO2-SRCE 
!        diagnostic.  Now references ND04 from "diag04_mod.f. 
!        (pns, bmy, 7/25/05)
!  (4 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (5 ) Now save MBO as tracer #5 for ND46 (tmf, bmy, 10/20/05)
!  (6 ) Now add categories CV-FLX-$, TURBMC-$, EW-FLX-$, NS-FLX-$, UP-FLX-$
!        which had been inadvertently omitted.  Also add OCEAN-HG category.
!        Rewrote do loop and case statement to add new diagnostics to ND03. 
!        Now make units of Hg tracers "pptv", not "ppbv".  Now remove 
!        restriction on printing out cloud mass flux in GEOS-4 for the ND66 
!        diagnostic.  Added new sea salt category. (cdh, eck, bmy, 4/6/06)
!  (7 ) Now references ND56 from "diag56_mod.f" (ltm, bmy, 5/5/06) 
!  (8 ) Now references ND42 from "diag42_mod.f".  Also updated for extra SOA
!        tracers in ND07 diagnostic. (dkh, bmy, 5/22/06)
!  (9 ) Updated ND36 for CH3I simulation (bmy, 7/25/06)
!  (10) Remove support for GEOS-1 and GEOS-STRAT met fields (bmy, 8/4/06)
!  (11) Split into INIT_DIAGINFO, INIT_TRACERINFO for clarity (bmy, 9/28/06)
!  (12) Save output to HDF_MOD (amv, bmy, 12/18/09)
!******************************************************************************
!
      ! References to F90 modules
      USE ERROR_MOD,   ONLY : ALLOC_ERR
      USE TIME_MOD,    ONLY : EXPAND_DATE, GET_NHMSb, GET_NYMDb
      USE LOGICAL_MOD, ONLY : LND50_HDF, LND51_HDF, LND51b_HDF

      USE HDF_MOD, ONLY : INIT_HDF
      USE HDF_MOD, ONLY : HDFCATEGORY
      USE HDF_MOD, ONLY : HDFDESCRIPT
      USE HDF_MOD, ONLY : HDFNAME
      USE HDF_MOD, ONLY : HDFFNAME
      USE HDF_MOD, ONLY : HDFUNIT
      USE HDF_MOD, ONLY : HDFMOLC
      USE HDF_MOD, ONLY : HDFMWT
      USE HDF_MOD, ONLY : HDFSCALE

#     include "CMN_SIZE"    ! Size parameters
#     include "CMN_DIAG"    ! NDxx flags

      ! Arguments
      CHARACTER(LEN=255), INTENT(IN) :: DIAGINFO  
      CHARACTER(LEN=255), INTENT(IN) :: TRACERINFO

      ! Local variables
      INTEGER                        :: AS, NYMDb, NHMSb

      !=================================================================
      ! INIT_GAMAP begins here!
      !=================================================================

      ! Save from arguments to module variables
      DFILE = DIAGINFO
      TFILE = TRACERINFO

      ! Get starting date & time
      NYMDb = GET_NYMDb()
      NHMSb = GET_NHMSb()
 
      ! Replace any date/time tokens in the file names
      CALL EXPAND_DATE( DFILE, NYMDb, NHMSb )
      CALL EXPAND_DATE( TFILE, NYMDb, NHMSb )

      !=================================================================
      ! Allocate module arrays
      !=================================================================

      ALLOCATE( OFFSET( MAXCAT ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'OFFSET' )
      OFFSET = 0

      ALLOCATE( CATEGORY( MAXCAT ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'CATEGORY' )
      CATEGORY = ''

      ALLOCATE( DESCRIPT( MAXCAT ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'DESCRIPT' )
      DESCRIPT = ''

      ALLOCATE( NTRAC( MAXDIAG ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'NTRAC' )
      NTRAC = 0

      ALLOCATE( INDEX( MAXTRACER, MAXDIAG ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'INDEX' )      
      INDEX = 0

      ALLOCATE( MOLC( MAXTRACER, MAXDIAG ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'MOLC' )      
      MOLC = 0

      ALLOCATE( MWT( MAXTRACER, MAXDIAG ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'MWT' )      
      MWT = 0.0

      ALLOCATE( SCALE( MAXTRACER, MAXDIAG ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'SCALE' )      
      SCALE = 0.0

      ALLOCATE( NAME( MAXTRACER, MAXDIAG ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'NAME' )      
      NAME = ''

      ALLOCATE( FNAME( MAXTRACER, MAXDIAG ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'FNAME' )      
      FNAME = ''

      ALLOCATE( UNIT( MAXTRACER, MAXDIAG ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'UNIT' )      
      UNIT = ''

      !=================================================================
      ! Initialize arrays for "diaginfo.dat" & "tracerinfo.dat" files
      !=================================================================

      ! Initialize arrays for "diaginfo.dat"
      CALL INIT_DIAGINFO

      ! Initialize arrays for "tracerinfo.dat"
      CALL INIT_TRACERINFO

      ! Store values for hdf output
      IF ( LND50_HDF .or. LND51_HDF .or. LND51b_HDF ) THEN

         CALL INIT_HDF(MAXCAT, MAXTRACER, MAXDIAG)

         HDFCATEGORY(:) = CATEGORY(:)
         HDFDESCRIPT(:) = DESCRIPT(:)
         HDFNAME(:,:) = NAME(:,:)
         HDFFNAME(:,:) = FNAME(:,:)
         HDFUNIT(:,:) = UNIT(:,:)
         HDFMOLC(:,:) = MOLC(:,:)
         HDFMWT(:,:) = MWT(:,:)
         HDFSCALE(:,:) = SCALE(:,:)
         
      ENDIF

      ! Return to calling program
      END SUBROUTINE INIT_GAMAP

!------------------------------------------------------------------------------

      SUBROUTINE CLEANUP_GAMAP
!
!******************************************************************************
!  Subroutine CLEANUP_GAMAP deallocates all module arrays (bmy, 4/25/05)
! 
!  NOTES:
!******************************************************************************
!
      !===================================================================
      ! CLEANUP_GAMAP begins here!
      !===================================================================
      IF ( ALLOCATED( CATEGORY ) ) DEALLOCATE( CATEGORY )
      IF ( ALLOCATED( DESCRIPT ) ) DEALLOCATE( DESCRIPT )
      IF ( ALLOCATED( FNAME    ) ) DEALLOCATE( FNAME    )
      IF ( ALLOCATED( INDEX    ) ) DEALLOCATE( INDEX    )
      IF ( ALLOCATED( MOLC     ) ) DEALLOCATE( MOLC     )
      IF ( ALLOCATED( MWT      ) ) DEALLOCATE( MWT      )
      IF ( ALLOCATED( NAME     ) ) DEALLOCATE( NAME     )
      IF ( ALLOCATED( NTRAC    ) ) DEALLOCATE( NTRAC    )
      IF ( ALLOCATED( OFFSET   ) ) DEALLOCATE( OFFSET   )
      IF ( ALLOCATED( SCALE    ) ) DEALLOCATE( SCALE    )
      IF ( ALLOCATED( UNIT     ) ) DEALLOCATE( UNIT     )

      ! Return to calling program
      END SUBROUTINE CLEANUP_GAMAP

!------------------------------------------------------------------------------

      ! End of module
      END MODULE GAMAP_MOD
