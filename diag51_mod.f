! $Id: diag51_mod.f,v 1.7 2004/04/19 15:09:52 bmy Exp $
      MODULE DIAG51_MOD
!
!******************************************************************************
!  Module DIAG51_MOD contains variables and routines to generate save 
!  timeseries data over the United States where the local time is between 
!  two user-defined limits. (amf, bey, bdf, pip, bmy, 11/30/00, 4/16/04)
!
!  Module Variables:
!  ============================================================================
!  (1 ) BLEVCT     (INTEGER ) : Counter to accompany BLLEVS array
!  (2 ) GOOD       (INTEGER ) : Index for grid boxes w/in a local time range
!  (3 ) PBLTOT     (INTEGER ) : Accumulator for PBL [m]
!  (4 ) WTOTH      (INTEGER ) : Indicates whe
!  (5 ) ISET       (INTEGER ) : Longitude offset
!  (6 ) JSET       (INTEGER ) : Latitude offset
!  (7 ) LSET       (INTEGER ) : Altitude offset
!  (8 ) BLARR      (REAL*8  ) : Afternoon points for PBL [m] 
!  (9 ) BLLEVS     (REAL*8  ) : Accumulator for PBL [model layers]
!  (10) STT_TEMP51 (REAL*8  ) : Accumulator for various quantities
!  (11) XLOCTM     (REAL*8  ) : Array of local times for each longitude
!  (12) ZTAU1      (REAL*8  ) : Starting TAU used to index the bpch file
!  (13) ZTAU2      (REAL*8  ) : Ending TAU used to index the bpch file
!  (14) NHMS_WRITE (INTEGER ) : Time to write average data to the bpch file
!  (15) HR1        (REAL*8  ) : Starting hour of user-defined LT interval
!  (16) HR2        (REAL*8  ) : Ending hour of user-defined LT interval
!  (17) FILENAME   (CHAR*255) : Name of bpch file containing timeseries data
!
!  Module Procedures:
!  ============================================================================
!  (1 ) DIAG51            : Driver subroutine for US grid timeseries 
!  (2 ) GET_LOCAL_TIME    : Computes the local times at each grid box
!  (3 ) WRITE_DIAG51      : Writes timeseries data to a binary punch file
!  (4 ) ACCUMULATE_DIAG51 : Accumulates data over the US for later averaging
!  (5 ) INIT_DIAG51       : Allocates and zeroes all module arrays 
!  (6 ) CLEANUP_DIAG51    : Deallocates all module arrays
!
!  GEOS-CHEM modules referenced by diag51_mod.f
!  ============================================================================
!  (1 ) bpch2_mod.f    : Module containing routines for binary punch file I/O
!  (2 ) dao_mod.f      : Module containing arrays for DAO met fields
!  (3 ) error_mod.f    : Module containing NaN and other error check routines
!  (4 ) file_mod.f     : Module containing file unit numbers and error checks
!  (5 ) grid_mod.f     : Module containing horizontal grid information
!  (6 ) pressure_mod.f : Module containing routines to compute P(I,J,L) 
!  (7 ) time_mod.f     : Module containing routines to compute date & time
!  (8 ) tracerid_mod.f : Module containing pointers to tracers & emissions
!
!  DIAG51 Tracers (specified in "timeseries.dat")
!  ============================================================================
!   1-31 : CTM Tracers               [v/v]
!   32   : Pure O3 (full chem only)  [v/v]
!   33   : Pure NO (full chem only)  [v/v]
!   34   : NOy (full chem only)      [v/v]
!   35   : OH  (full chem only)      [molec/cm3/s]
!   36   : Pure NO2 (full chem only) [v/v]
!   37   : PBL depth                 [m]
!   38   : PBL depth                 [# of model layers]
!   39   : Air density               [molec/cm3/s]
!   42   : Cloud fraction            [unitless]
!   43   : Column Optical Depths     [unitless]
!   44   : Cloud Top Heights         [mb]
!   98   : PS - Ptop                 [mb]
!   99   : Temperature               [K]
!
!  NOTES:
!  (1 ) "PRIVATE" means that the module variables will not be "seen" anywhere
!         outside of this module.  This prevents overwriting. (bmy, 12/1/00)
!  (2 ) Removed obsolete code from WRITE_DIAG51 (bmy, 12/18/00)
!  (3 ) Updated comments (bmy, 9/4/01)
!  (4 ) XTRA2(IREF,JREF,5) is now XTRA2(I,J) (bmy, 9/25/01)
!  (5 ) Replace PW with P (bmy, 10/3/01)
!  (6 ) Removed obsolete code from 9/01 (bmy, 10/23/01)
!  (7 ) Now specify the first box correctly for window regions smaller
!        than global size (bmy, 12/18/01)
!  (8 ) Eliminate obsolete code from 12/01 (bmy, 2/27/02)
!  (9 ) Now divide module header into MODULE PRIVATE, MODULE VARIABLES, and
!        MODULE ROUTINES sections.  Updated comments (bmy, 5/28/02)
!  (10) Now reference "file_mod.f" (bmy, 6/26/02)
!  (11) Replaced references to P(I,J) with call to GET_PEDGE(I,J,1) from 
!        "pressure_mod.f".  Also eliminated obsolete, commented-out code
!         from 6/02. (dsa, bdf, bmy, 8/20/02)
!  (12) Now reference AD, BXHEIGHT, and T from "dao_mod.f".  Also removed
!        obsolete code from various routines.  Now references "error_mod.f".
!        Now references "tracerid_mod.f". (bmy, 10/15/02)
!  (13) Change tracer #'s for sulfate tracers beyond 24.  Also updated
!        comments (rjp, bmy, 3/23/03)
!  (14) Now references "time_mod.f" and "grid_mod.f" (bmy, 3/27/03)
!  (15) Bug fix for LINUX in calls to TIMESTAMP_STRING (bmy, 9/29/03)
!  (16) Bug fix in WRITE_DIAG51.  Also multiply AIRDEN by 1d-6 to convert
!        from m3 to cm3. (bmy, 3/24/04)
!  (17) Update tracer #'s for 41 regular CTM tracers (bmy, 4/16/04)
!******************************************************************************
!
      IMPLICIT NONE

      !=================================================================
      ! MODULE PRIVATE DECLARATIONS -- keep certain internal variables 
      ! and routines from being seen outside "diag51_mod.f"
      !=================================================================

      ! PRIVATE module variables
      PRIVATE :: BLEVCT,     GOOD,       PBLTOT
      PRIVATE :: WTOTH,      ISET,       JSET
      PRIVATE :: LSET,       BLARR,      BLLEVS
      PRIVATE :: STT_TEMP51, XLOCTM,     ZTAU1
      PRIVATE :: ZTAU2,      NHMS_WRITE, HR1  
      PRIVATE :: HR2,        FILENAME 

      ! PRIVATE module routines
      PRIVATE :: GET_LOCAL_TIME    
      PRIVATE :: WRITE_DIAG51      
      PRIVATE :: ACCUMULATE_DIAG51 
      PRIVATE :: INIT_DIAG51       

      !=================================================================
      ! MODULE VARIABLES
      !=================================================================
      INTEGER, ALLOCATABLE :: BLEVCT(:,:)
      INTEGER, ALLOCATABLE :: GOOD(:)
      INTEGER, ALLOCATABLE :: PBLTOT(:,:)
      INTEGER, ALLOCATABLE :: WTOTH(:,:)
      INTEGER              :: ISET, JSET, LSET
                           
      REAL*8,  ALLOCATABLE :: BLARR(:,:)
      REAL*8,  ALLOCATABLE :: BLLEVS(:,:)
      REAL*8,  ALLOCATABLE :: STT_TEMP51(:,:,:,:)
      REAL*8,  ALLOCATABLE :: XLOCTM(:) 
      REAL*8               :: ZTAU1, ZTAU2

      !=================================================================
      ! For timeseries between 1300 and 1700 LT, uncomment this code:
      !
      ! Need to write to the bpch file at 12 GMT, since this covers
      ! an entire day over the US grid (amf, bmy, 12/1/00)
      !
      INTEGER, PARAMETER   :: NHMS_WRITE = 120000
      REAL*8,  PARAMETER   :: HR1        = 13d0
      REAL*8,  PARAMETER   :: HR2        = 17d0
      CHARACTER(LEN=255)   :: FILENAME   = 'ts1_4pm.bpch'
      !=================================================================
      ! For timeseries between 1000 and 1200 LT, uncomment this code:
      !
      ! Between 10 and 12 has been chosen because the off-polar orbit 
      ! of GOME traverses (westward) through local times between 12 
      ! and 10 over North America, finally crossing the equator at 
      ! 10.30 (local time).
      !
      ! Need to write to the bpch file at 00 GMT, since we will be 
      ! interested in the whole northern hemisphere (pip, 12/1/00)
      !
      !INTEGER, PARAMETER   :: NHMS_WRITE = 000000
      !REAL*8,  PARAMETER   :: HR1        = 10d0
      !REAL*8,  PARAMETER   :: HR2        = 12d0
      !CHARACTER(LEN=255)   :: FILENAME   ='ts10_12pm.bpch'
      !=================================================================

      !=================================================================
      ! MODULE ROUTINES -- follow below the "CONTAINS" statement 
      !=================================================================
      CONTAINS

!------------------------------------------------------------------------------

      SUBROUTINE DIAG51
!
!******************************************************************************
!  Subroutine DIAG51 generates time series (averages from 10am - 12pm LT 
!  or 1pm - 4pm LT) for the US grid area.  Output is to binary punch files.
!  (amf, bey, bdf, pip, bmy, 11/15/99, 2/27/02)
!
!  *** NEW *** DIAG51 was not archiving for the proper day for the west coast 
!  because the end-of-day test was based on GMT midnight, which was still only
!  3 or 4 pm in california.  Provide user the option of using GMT noon (for
!  US grid) or GMT midnight (elsewhere) as the end-of-day test criteria.
!  (amf, bey, bdf, 11/30/00)
!
!  NOTES:
!  (1 ) Added to "diag51_mod.f" (bmy, 11/30/00)
!  (2 ) Split off sections into separate module subroutines (bmy, 11/30/00)
!  (3 ) TAU is referenced to GMT.  For US grid box, now shift the punch
!        file write times to 12 GMT, to cover entire day.  However, still
!        keep TAU in the bpch file as GMT, for consistency with GAMAP.
!        (amf, bmy, 12/1/00)
!  (4 ) Make sure NOx is a defined tracer before archiving NO2 concentrations 
!        to the timeseries file.  Now save out pure NO2 mixing ratio, by using
!        the FRACNO2 array -- the fraction of NOx that is NO2.  Updated
!        comments (rvm, bmy, 2/27/02)
!  (5 ) Now use functions GET_NHMS, GET_TAUb, GET_TAUe, GET_TAU from
!        "time_mod.f".  Removed NHMS from arg list. (bmy, 3/27/03)
!******************************************************************************
!
      ! References to F90 modules
      USE TIME_MOD, ONLY : GET_NHMS, GET_TAUb, GET_TAUe, GET_TAU

#     include "CMN_SIZE"     ! Size parameters
#     include "CMN"          ! TOFDAY, TAU, TAUI, TAUE
#     include "CMN_TIMES"    ! IMIN_AREA, JMIN_AREA, LMIN_AREA

      ! Local variables
      LOGICAL, SAVE       :: FIRSTDIAG51 = .TRUE.

      !=================================================================
      ! DIAG51 begins here!
      !
      ! Initialization -- first time only
      !=================================================================
      IF ( FIRSTDIAG51 ) THEN
         
         ! Allocate and zero accumulating arrays
         CALL INIT_DIAG51

         ! STARTING TAU for the first punch file write time
         ZTAU1 = GET_TAUb()

         ! Calculate offset to correctly access where 
         ! STT_TEMP51 is in the global STT array:
         ISET = IMIN_AREA - 1
         JSET = JMIN_AREA - 1
         LSET = LMIN_AREA - 1
      ENDIF  

      !=================================================================
      ! Compute local times
      !=================================================================
      CALL GET_LOCAL_TIME

      !=================================================================
      ! If it is a new day, call WRITE_DIAG51 to:
      ! * compute the average values of the selected quantities 
      ! * write the average values in the timeseries file
      ! * reset variables to zero
      !
      ! NOTE: Do not write to disk if this is the first timestep,
      !       since nothing will have been archived yet. 
      !=================================================================
      IF ( ( GET_NHMS() == NHMS_WRITE   .or. 
     &       GET_TAU()  == GET_TAUe() ) .and. 
     &     ( .not. FIRSTDIAG51 ) ) THEN
    
         ! Set ENDING TAU for this bpch file write time
         ZTAU2 = GET_TAU()

         ! Compute averages and write timeseries data to bpch file
         CALL WRITE_DIAG51

         ! Set STARTING TAU for the next bpch file write time
         ZTAU1 = GET_TAU()
      ENDIF

      !=================================================================
      ! Accumulate quantities where local time is between HR1 and HR2
      !=================================================================
      CALL ACCUMULATE_DIAG51

      !=================================================================
      ! We are no longer in the first ND51 step -- reset FIRSTDIAG51
      !=================================================================
      FIRSTDIAG51 = .FALSE.

      ! Return to calling program
      END SUBROUTINE DIAG51

!------------------------------------------------------------------------------

      SUBROUTINE GET_LOCAL_TIME
!
!******************************************************************************
!  Subroutine GET_LOCAL_TIME computes the local time and returns an array 
!  of points where the local time is between two user-defined limits. 
!  (bmy, 11/29/00, 3/27/03)
!
!  NOTES:
!  (1 ) The 1d-3 in the computation of XLOCTM is to remove roundoff ambiguity 
!        if a the local time should fall exactly on an hour boundary.
!        (bmy, 11/29/00)
!  (2 ) Bug fix: XMID(I) should be XMID(II).  Also updated comments.
!        (bmy, 7/6/01)
!  (3 ) Updated comments (rvm, bmy, 2/27/02)
!  (4 ) Now uses function GET_LOCALTIME of "time_mod.f" (bmy, 3/27/03) 
!******************************************************************************
!
      ! References to F90 modules
      USE TIME_MOD, ONLY : GET_LOCALTIME

#     include "CMN_SIZE"   ! Size parameters
#     include "CMN"        ! TOFDAY, XMID

      ! Local variables 
      INTEGER              :: II
      INTEGER, PARAMETER   :: ITEST = 2

      !=================================================================
      ! GET_LOCAL_TIME begins here!
      !
      ! Compute local times for all locations
      !=================================================================
      DO II = 1, IIPAR

         ! Local time = GMT + ( longitude in degrees / 15 degrees per hour )
         XLOCTM(II) = GET_LOCALTIME(II)

         ! GOOD indicates which boxes have local times between HR1 and HR2
         IF ( XLOCTM(II) >= HR1 .and. XLOCTM(II) <= HR2 ) THEN
            GOOD(II) = 1
         ELSE
            GOOD(II) = 0
         ENDIF
      ENDDO

      ! Return to calling program
      END SUBROUTINE GET_LOCAL_TIME

!------------------------------------------------------------------------------

      SUBROUTINE WRITE_DIAG51
!
!******************************************************************************
!  Subroutine WRITE_DIAG51 (bmy, 12/1/00, 4/16/04) does the following:
!  
!  (1) Divides quantities by the number of times it was between
!       HR1 and HR2 local time in each grid box  
!  (2) Writes average data to a binary punch file
!  (3) Zeroes accumulating and counter arrays 
!
!  NOTES:
!  (1 ) Split off from the old "diag51.f" and added to F90 module
!       "diag51_mod.f" (bmy, 12/1/00)
!  (2 ) Now allocate and deallocate the ZZ array within WRITE_DIAG51.
!  (3 ) Readjust GAMAP tracer number for NO2 (bmy, 12/7/00)
!  (4 ) Removed obsolete code from 12/7/00 (bmy, 12/18/00)
!  (5 ) Now pass IMIN_AREA+I0 as IFIRST and JMIN_AREA+J0 as JFIRST in call 
!        to BPCH2.  This will take care of window regions smaller than the 
!        globe. (yxw, bmy, 12/18/01)
!  (6 ) Eliminate obsolete code from 12/01 (bmy, 2/27/02)
!  (7 ) Now reference IU_ND51 from "file_mod.f".  Now use IU_ND51 instead of 
!        IUT as the file unit #. Also call routine OPEN_BPCH2_FOR_WRITE to 
!        start writing to the output file. (bmy, 6/26/02)
!  (8 ) Deleted obsolete, commented-out code (bmy, 8/20/02)
!  (9 ) Now reference ALLOC_ERR from "error_mod.f" (bmy, 10/15/02)
!  (10) Change tracer #'s for sulfate tracers beyond 24.  Also updated
!        comments (rjp, bmy, 3/23/03)
!  (11) Now uses functions GET_XOFFSET and GET_YOFFSET from "grid_mod.f". 
!        I0 and J0 are now local variables.  Now uses TIMESTAMP_STRING and
!        ITS_TIME_FOR_EXIT from "time_mod.f" (bmy, 3/27/03)     
!  (12) LINUX has a problem putting a function call w/in a WRITE statement.  
!        Now save output from TIMESTAMP_STRING to STAMP and print that.
!        (bmy, 9/29/03)
!  (13) Bug fix: GMNL should be NL for tracer #39, not 1.  Also multiply 
!        AIRDEN by 1d-6 to convert from m3 to cm3. (bmy, 3/23/04)
!  (14) Readjust tracer numbers for other fields to accomodate 41 regular CTM 
!        tracers -- QUICK FIX (bmy, 4/16/04)
!******************************************************************************
!
      ! Reference to F90 modules
      USE BPCH2_MOD
      USE ERROR_MOD, ONLY : ALLOC_ERR
      USE FILE_MOD,  ONLY : IU_ND51
      USE GRID_MOD,  ONLY : GET_XOFFSET, GET_YOFFSET
      USE TIME_MOD,  ONLY : ITS_TIME_FOR_EXIT, TIMESTAMP_STRING

#     include "CMN_SIZE"  ! Size Parameters
#     include "CMN"       ! TAU, TOFDAY, STT, etc.
#     include "CMN_DIAG"  ! TRCOFFSET
#     include "CMN_TIMES" ! NI, NJ, NL, etc.

      ! Local variables
      LOGICAL             :: FIRST = .TRUE.

      INTEGER             :: I, J, L, N, GMNL, GMTRC, IOS, XTRAC, I0, J0
      INTEGER, PARAMETER  :: HALFPOLAR=1, CENTER180=1
      
      REAL*4              :: LONRES, LATRES
      REAL*4, ALLOCATABLE :: ZZ(:,:,:)

      CHARACTER(LEN=20)   :: MODELNAME
      CHARACTER(LEN=40)   :: CATEGORY
      CHARACTER(LEN=40)   :: UNIT     = '' 
      CHARACTER(LEN=40)   :: RESERVED = ''
      CHARACTER(LEN=80)   :: TITLE

      ! For LINUX fix (bmy, 9/29/03)
      CHARACTER(LEN=16)   :: STAMP

      !=================================================================
      ! WRITE_DIAG51 begins here!
      !=================================================================

      ! Get nested-grid offsets
      I0 = GET_XOFFSET( GLOBAL=.TRUE. )
      J0 = GET_YOFFSET( GLOBAL=.TRUE. )

      ! Allocate ZZ as a local array to this subroutine
      ALLOCATE( ZZ( NI, NJ, NL ), STAT=IOS )
      IF ( IOS /= 0 ) CALL ALLOC_ERR( 'ZZ' )
      ZZ(:,:,:) = 0e0

      !=================================================================
      ! Information for the header of the timeseries file
      !=================================================================
      LONRES    = DISIZE
      LATRES    = DJSIZE
      MODELNAME = GET_MODELNAME()

      !=================================================================
      ! If this is the first call to WRITE_DIAG51, then open 
      ! the binary timeseries file and write the file header.
      !=================================================================
      IF ( FIRST ) THEN

         ! Write header 
         TITLE = 'GEOS-CTM time series for geographical domain'

         ! Echo info
         WRITE( 6, 100 ) TRIM( FILENAME )
 100     FORMAT( '     - DIAG51: Opening file ', a ) 

         ! Open output file
         CALL OPEN_BPCH2_FOR_WRITE( IU_ND51, FILENAME, TITLE )

         ! Reset FIRST flag
         FIRST = .FALSE.
      ENDIF

      ! Echo info
      STAMP = TIMESTAMP_STRING()
      WRITE( 6, 110 ) STAMP
 110  FORMAT( '     - DIAG51: Saving to disk at ', a ) 

      !=================================================================
      ! Compute averages of quantities where it is between HR1 and
      ! HR2 local time -- Use a parallel DO-loop for efficiency
      ! Also avoid divide-by-zero errors!
      !=================================================================
!$OMP PARALLEL DO 
!$OMP+DEFAULT( SHARED ) 
!$OMP+PRIVATE( I, J, L, N )
      DO N = 1, NTRAC_AREA
      DO L = 1, NL
      DO J = 1, NJ
      DO I = 1, NI

         ! Averages for 3-D quantities
         IF ( WTOTH(I,J) > 0 ) THEN
            STT_TEMP51(I,J,L,N) = STT_TEMP51(I,J,L,N) / WTOTH(I,J) 
         ELSE
            STT_TEMP51(I,J,L,N) = 0d0
         ENDIF

         ! Averages for 2-D boundary layer arrays
         IF ( N == 1 .and. L == 1 ) THEN 

            ! Boundary layer [m]
            IF ( PBLTOT(I,J) > 0 ) THEN 
               BLARR(I,J) = BLARR(I,J) / PBLTOT(I,J)               
            ELSE     
               BLARR(I,J) = 0d0
            ENDIF

            ! Boundary layer [model layers]
            IF ( BLEVCT(I,J) > 0 ) THEN
               BLLEVS(I,J) = BLLEVS(I,J) / BLEVCT(I,J)               
            ELSE
               BLLEVS(I,J) = 0d0
            ENDIF
         ENDIF
      ENDDO
      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      !=================================================================
      ! Write each tracer from "timeseries.dat" to the timeseries file
      !=================================================================
      DO N = 1, NTRAC_AREA
         XTRAC = TRAC_AREA(N)

         ! Add correct tracer number (GMTRC) for GAMAP
         ! Save the proper number of levels (GMNL) to the punch file
         SELECT CASE ( XTRAC ) 

            ! CTM Tracers [v/v] (3-D field)
            !-------------------
            ! Prior to 4/16/04:
            !CASE ( 1:31 )
            !-------------------
            CASE ( 1:41 )
               ZZ(1:NI, 1:NJ, 1:NL) = STT_TEMP51(1:NI, 1:NJ, 1:NL, N)
               CATEGORY             = 'IJ-AVG-$'
               GMNL                 = NL
               GMTRC                = XTRAC + TRCOFFSET

            ! Pure O3 (3-D field)
            !-------------------
            ! Prior to 4/16/04:
            !CASE ( 32 )
            !-------------------
            CASE ( 71 )
               ZZ(1:NI, 1:NJ, 1:NL) = STT_TEMP51(1:NI, 1:NJ, 1:NL, N)
               CATEGORY             = 'IJ-AVG-$'
               GMNL                 = NL
               GMTRC                = XTRAC + TRCOFFSET

            ! Pure NO (3-D field)
            !-------------------
            ! Prior to 4/16/04:
            !CASE ( 33 )
            !-------------------
            CASE ( 72 )
               ZZ(1:NI, 1:NJ, 1:NL) = STT_TEMP51(1:NI, 1:NJ, 1:NL, N)
               CATEGORY             = 'TIME-SER'
               GMNL                 = NL
               GMTRC                = 9

            ! NOy (3-D field)
            !-------------------
            ! Prior to 4/16/04:
            !CASE ( 34 )
            !-------------------
            CASE ( 73 )
               ZZ(1:NI, 1:NJ, 1:NL) = STT_TEMP51(1:NI, 1:NJ, 1:NL, N)
               CATEGORY             = 'TIME-SER'
               GMNL                 = NL
               GMTRC                = 3

            ! OH (3-D field)
            !-------------------
            ! Prior to 4/16/04:
            !CASE ( 35 )
            !-------------------
            CASE ( 74 )
               ZZ(1:NI, 1:NJ, 1:NL) = STT_TEMP51(1:NI, 1:NJ, 1:NL, N)
               CATEGORY             = 'TIME-SER'
               GMNL                 = NL
               GMTRC                = 2

            ! NO2 (3-D field)
            !-------------------
            ! Prior to 4/16/04:
            !CASE ( 36 )
            !-------------------
            CASE ( 75 )
               ZZ(1:NI, 1:NJ, 1:NL) = STT_TEMP51(1:NI, 1:NJ, 1:NL, N)
               CATEGORY             = 'TIME-SER'
               GMNL                 = NL
               GMTRC                = 23

            ! PBL Height -- meters (1-D field)
            !-------------------
            ! Prior to 4/16/04:
            !CASE ( 37 )
            !-------------------
            CASE ( 76 )
               ZZ(1:NI, 1:NJ, 1)    = BLARR(1:NI, 1:NJ)
               CATEGORY             = 'PBLDEPTH'
               GMNL                 = 1
               GMTRC                = 1

            ! PBL Height -- model layers (1-D field)
            !-------------------
            ! Prior to 4/16/04:
            !CASE ( 38 )
            !-------------------
            CASE ( 77 )
               ZZ(1:NI, 1:NJ, 1)    = BLLEVS(1:NI, 1:NJ)
               CATEGORY             = 'PBLDEPTH'
               GMNL                 = 1
               GMTRC                = 2

            ! Air Density (3-D field)
            !-------------------
            ! Prior to 4/16/04:
            !CASE ( 39 )
            !-------------------
            CASE ( 78 )
               ZZ(1:NI, 1:NJ, 1:NL) = STT_TEMP51(1:NI, 1:NJ, 1:NL, N)
               CATEGORY             = 'TIME-SER'
               GMNL                 = NL
               GMTRC                = 22

            ! Cloud fractions (1-D field)
            !-------------------
            ! Prior to 4/16/04:
            !CASE ( 42 )
            !-------------------
            CASE ( 79 )
               ZZ(1:NI, 1:NJ, 1)    = STT_TEMP51(1:NI, 1:NJ, 1, N)
               CATEGORY             = 'TIME-SER'
               GMNL                 = 1
               GMTRC                = 19

            ! Optical depths (1-D field)
            !-------------------
            ! Prior to 4/16/04:
            !CASE ( 43 ) 
            !-------------------
            CASE ( 80 ) 
               ZZ(1:NI, 1:NJ, 1)    = STT_TEMP51(1:NI, 1:NJ, 1, N)
               CATEGORY             = 'TIME-SER'
               GMNL                 = 1
               GMTRC                = 20
                    
            ! Cloud top heights (1-D field)
            !-------------------
            ! Prior to 4/16/04:
            !CASE ( 44 ) 
            !-------------------
            CASE ( 81 ) 
               ZZ(1:NI, 1:NJ, 1)    = STT_TEMP51(1:NI, 1:NJ, 1, N)
               CATEGORY             = 'TIME-SER'
               GMNL                 = 1
               GMTRC                = 21

            ! Psurface - PTOP (1-D field)
            CASE ( 98 )
               ZZ(1:NI, 1:NJ, 1)    = STT_TEMP51(1:NI, 1:NJ, 1, N)
               CATEGORY             = 'PS-PTOP'
               GMNL                 = 1
               GMTRC                = 1

            ! Temperature (3-D field)
            CASE ( 99 )
               ZZ(1:NI, 1:NJ, 1:NL) = STT_TEMP51(1:NI, 1:NJ, 1:NL, N)
               CATEGORY             = 'DAO-3D-$'
               GMNL                 = NL
               GMTRC                = 3
            
            ! Everything else is an invalid tracer -- skip!
            CASE DEFAULT
               CYCLE
            
         END SELECT

         ! Save each field to the binary punch file
         CALL BPCH2( IU_ND51,      MODELNAME,    LONRES,   
     &               LATRES,       HALFPOLAR,    CENTER180, 
     &               CATEGORY,     GMTRC,        UNIT,      
     &               ZTAU1,        ZTAU2,        RESERVED,  
     &               NI,           NJ,           GMNL,     
     &               IMIN_AREA+I0, JMIN_AREA+J0, LMIN_AREA, 
     &               ZZ(1:NI, 1:NJ, 1:GMNL) )
      ENDDO

      !=================================================================
      ! Flush the buffer to disk after each write
      ! If this is the ending time, close the file
      !=================================================================
      CALL FLUSH( IU_ND51 )

      IF ( ITS_TIME_FOR_EXIT() ) THEN
         CLOSE( IU_ND51 )
         WRITE( 6, 120 ) TRIM( FILENAME )
 120     FORMAT( '     - DIAG51: Closing file ', a )
      ENDIF

      !=================================================================
      ! Zero diagnostic and counter arrays after writing to disk
      ! Use a parallel DO-loop for efficiency
      !=================================================================
      STAMP = TIMESTAMP_STRING()
      WRITE( 6, 130 ) STAMP
 130  FORMAT( '     - DIAG51: Zeroing arrays at ', a )

!$OMP PARALLEL DO 
!$OMP+DEFAULT( SHARED ) 
!$OMP+PRIVATE( I, J, L, N )
      DO N = 1, NTRAC_AREA
      DO L = 1, NL
      DO J = 1, NJ
      DO I = 1, NI

         ! Zero accumulating array for tracer
         STT_TEMP51(I,J,L,N) = 0d0

         ! Zero 2-D accumulating and counter arrays
         IF ( N == 1 ) THEN
            IF ( L == 1 ) THEN 
               BLARR(I,J)  = 0d0
               BLLEVS(I,J) = 0d0
               BLEVCT(I,J) = 0
               PBLTOT(I,J) = 0
               WTOTH(I,J)  = 0
            ENDIF
         ENDIF
      ENDDO
      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO
      
      !=================================================================
      ! Deallocate the local ZZ array 
      !=================================================================
      DEALLOCATE( ZZ )
      
      ! Return to calling program
      END SUBROUTINE WRITE_DIAG51
      
!------------------------------------------------------------------------------

      SUBROUTINE ACCUMULATE_DIAG51
!
!******************************************************************************
!  Subroutine ACCUMULATE_DIAG51 accumulates tracers into the STT_TEMP51
!  array. (bmy, 8/20/02, 9/29/03)
!
!  NOTES:
!  (1 ) Now reference routine GET_PEDGE from "pressure_mod.f", which
!        returns the correct "floating" pressure. (dsa, bdf, bmy, 8/20/02)
!  (2 ) Now reference AD, BXHEIGHT, T from "dao_mod.f".  Now reference
!        F90 module "tracerid_mod.f" (bmy, 11/6/02)
!  (3 ) Change tracer #'s for sulfate tracers beyond 24.  Also updated
!        comments (rjp, bmy, 3/23/03)
!  (4 ) Now uses routine TIMESTAMP_STRING from "time_mod.f" (bmy, 3/14/03)
!  (5 ) LINUX has a problem putting a function call w/in a WRITE statement.  
!        Now save output from TIMESTAMP_STRING to STAMP and print that.
!        (bmy, 9/29/03)
!******************************************************************************
!
      ! Reference to F90 modules
      USE DAO_MOD,      ONLY : AD,     AIRDEN, BXHEIGHT, CLDTOPS, 
     &                         CLMOSW, CLROSW, OPTD,     T
      USE PRESSURE_MOD, ONLY : GET_PEDGE
      USE TIME_MOD,     ONLY : TIMESTAMP_STRING
      USE TRACERID_MOD

#     include "CMN_SIZE"  ! Size parameters 
#     include "CMN"       ! NSRCX, STT, TCVV, etc.
#     include "CMN_TIMES" ! NI, NJ, NL
#     include "CMN_O3"    ! FRACO3, FRACNO, SAVEO3, SAVENO2, SAVEHO2, FRACNO2
#     include "CMN_SETUP" ! HR1_NO, HR2_NO
#     include "CMN_DIAG"  ! TRCOFFSET

      ! Local variables
      INTEGER             :: I, J, L, M, N, II, JJ, LL, PBLINT, XTRAC
      REAL*8              :: C1, C2, PBLDEC, TEMPBL, XX

      ! For LINUX fix
      CHARACTER(LEN=16)   :: STAMP

      !=================================================================
      ! ACCUMULATE_DIAG51 begins here!
      !
      ! Echo time information to the screen
      !=================================================================
      STAMP = TIMESTAMP_STRING()
      WRITE( 6, 100 ) STAMP
 100  FORMAT( '     - DIAG51: Accumulation at ', a )

      !=================================================================
      ! Archive WTOTH outside the TRACER LOOP 
      ! WTOTH is the cumulative counter of afternoon points
      !=================================================================
      DO J = 1, NJ
         JJ = JSET + J
         DO I = 1, NI
            II = ISET + I                 
            WTOTH(I,J) = WTOTH(I,J) + GOOD(II)
         ENDDO
      ENDDO 
      
      !=================================================================
      ! Archive tracers from "timeseries.dat" into accumulating array
      !
      ! Use parallel DO-loops to archive tracers, for efficiency
      !=================================================================
!$OMP PARALLEL DO 
!$OMP+DEFAULT( SHARED ) 
!$OMP+PRIVATE( I, J, L, M, N, II, JJ, LL, XX, XTRAC )
!$OMP+PRIVATE( PBLDEC, PBLINT, TEMPBL, C1, C2  )
!$OMP+SCHEDULE( DYNAMIC )
      DO N = 1, NTRAC_AREA
         XTRAC = TRAC_AREA(N)

         ! Case statement for tracers
         SELECT CASE ( XTRAC )

            !===========================================================
            ! Archive CHEMICAL TRACERS [v/v]
            !===========================================================
            CASE ( 1:31 )

               DO L = 1, NL
                  LL = LSET + L
               DO J = 1, NJ
                  JJ = JSET + J
               DO I = 1, NI
                  II = ISET + I

                  ! Archive afternoon points in ND51 
                  STT_TEMP51(I,J,L,N) = STT_TEMP51(I,J,L,N) + 
     &                 ( STT(II,JJ,LL,XTRAC) * TCVV(XTRAC) / 
     &                   AD(II,JJ,LL)        * GOOD(II) )
               ENDDO
               ENDDO
               ENDDO

            !===========================================================
            ! Archive OZONE [v/v] -- only for full chemistry
            !==========================================================
            CASE ( 32 )
               
               ! Make sure Ox is a tracer
               IF ( NSRCX == 3 .and. IDTOX > 0 ) THEN
                  DO L = 1, NL
                     LL = LSET + L
                  DO J = 1, NJ
                     JJ = JSET + J
                  DO I = 1, NI
                     II = ISET + I       

                     ! Pure O3 mixing ratio
                     XX = STT(II,JJ,LL,IDTOX) * FRACO3(II,JJ,LL) *
     &                    TCVV(IDTOX)         / AD(II,JJ,LL)     * 
     &                    GOOD(II)

                     ! Archive the afternoon points in STT_TEMP51
                     STT_TEMP51(I,J,L,N) = STT_TEMP51(I,J,L,N) + XX
                  ENDDO
                  ENDDO
                  ENDDO
               ENDIF


            !===========================================================
            ! Archive NO [v/v] -- only for full chemistry
            !===========================================================
            CASE ( 33 )
                  
               ! Make sure NOx is a tracer
               IF ( NSRCX == 3 .and. IDTNOX > 0 ) THEN
                  DO L = 1, NL
                     LL = LSET + L
                  DO J = 1, NJ
                     JJ = JSET + J
                  DO I = 1, NI
                     II = ISET + I     
  
                     ! Pure NO mixing ratio
                     XX = STT(II,JJ,LL,IDTNOX) * FRACNO(II,JJ,LL) *
     &                    TCVV(IDTNOX)         / AD(II,JJ,LL)     * 
     &                    GOOD(II)
                                             
                     ! Archive the afternoon points in STT_TEMP51
                     STT_TEMP51(I,J,L,N) = STT_TEMP51(I,J,L,N) + XX
                  ENDDO
                  ENDDO
                  ENDDO
               ENDIF


            !===========================================================
            ! Archive NOy [v/v] -- only for full chemistry
            !===========================================================
            CASE ( 34 )

               ! Full chemistry only...
               IF ( NSRCX == 3 ) THEN 
                  DO L = 1, NL
                     LL = LSET + L
                  DO J = 1, NJ
                     JJ = JSET + J
                  DO I = 1, NI
                     II = ISET + I   
  
                     ! Temp variable for accumulation
                     XX = 0d0
            
                     ! NOx
                     IF ( IDTNOX > 0 ) THEN
                        XX = XX + ( TCVV(IDTNOX)         * 
     &                              STT(II,JJ,LL,IDTNOX) / 
     &                              AD(II,JJ,LL) )
                     ENDIF 

                     ! PAN
                     IF ( IDTPAN /= 0 ) THEN
                        XX = XX + ( TCVV(IDTPAN)         * 
     &                              STT(II,JJ,LL,IDTPAN) / 
     &                              AD(II,JJ,LL) )
                     ENDIF

                     ! HNO3
                     IF ( IDTHNO3 /= 0 ) THEN
                        XX = XX + ( TCVV(IDTHNO3)         *
     &                              STT(II,JJ,LL,IDTHNO3) / 
     &                              AD(II,JJ,LL) )
                     ENDIF
            
                     ! PMN
                     IF ( IDTPMN /= 0 ) THEN
                        XX = XX + ( TCVV(IDTPMN)         *
     &                              STT(II,JJ,LL,IDTPMN) / 
     &                              AD(II,JJ,LL) )
                     ENDIF

                     ! PPN
                     IF ( IDTPPN /= 0 ) THEN
                        XX = XX + ( TCVV(IDTPPN)         * 
     &                              STT(II,JJ,LL,IDTPPN) / 
     &                              AD(II,JJ,LL) )
                     ENDIF 

                     ! ISN2
                     IF ( IDTISN2 /= 0 ) THEN 
                        XX = XX + ( TCVV(IDTISN2)         *
     &                              STT(II,JJ,LL,IDTISN2) / 
     &                              AD(II,JJ,LL) )
                     ENDIF
            
                     ! R4N2
                     IF ( IDTR4N2 /= 0 ) THEN
                        XX = XX + ( TCVV(IDTR4N2)         *
     &                              STT(II,JJ,LL,IDTR4N2) / 
     &                              AD(II,JJ,LL) )
                     ENDIF
            
                     ! N2O5
                     IF ( IDTN2O5 /= 0 ) THEN
                        XX = XX + ( 2 * TCVV(IDTN2O5)     * 
     &                              STT(II,JJ,LL,IDTN2O5) / 
     &                              AD(II,JJ,LL) )
                     ENDIF
                        
                     ! HNO4
                     IF ( IDTHNO4 /= 0 ) THEN
                        XX = XX + ( TCVV(IDTHNO4)         * 
     &                              STT(II,JJ,LL,IDTHNO4) / 
     &                              AD(II,JJ,LL) )
                     ENDIF

                     STT_TEMP51(I,J,L,N) = STT_TEMP51(I,J,L,N) + XX 
                  ENDDO
                  ENDDO
                  ENDDO
               ENDIF
                        
            !===========================================================
            ! Archive OH [molec/cm3/s] -- only for full chemistry
            ! OH is stored in the "SAVEOH" array from "CMN_O3"
            !===========================================================
            CASE ( 35 )

               ! Full chemistry only...
               IF ( NSRCX == 3 ) THEN
                  DO L = 1, NL
                     LL = LSET + L
                  DO J = 1, NJ
                     JJ = JSET + J
                  DO I = 1, NI
                     II = ISET + I  

                     ! Archive afternoon points in STT_TEMP51
                     STT_TEMP51(I,J,L,N) = STT_TEMP51(I,J,L,N) + 
     &                    ( SAVEOH(II,JJ,LL) * GOOD(II) )
                  ENDDO
                  ENDDO
                  ENDDO
               ENDIF
              
            !===========================================================
            ! Archive NO2 [molec/cm3/s] -- only for full chemistry
            ! OH is stored in the "SAVENO2" array from "CMN_O3"
            !===========================================================
            CASE ( 36 )

               ! Full chemistry only...
               IF ( NSRCX == 3 .and. IDTNOX > 0 ) THEN
                  DO L = 1, NL
                     LL = LSET + L
                  DO J = 1, NJ
                     JJ = JSET + J
                  DO I = 1, NI
                     II = ISET + I  

                     ! Pure NO2 mixing ratio (rvm, bmy, 2/27/02)
                     XX = STT(II,JJ,LL,IDTNOX) * FRACNO2(II,JJ,LL) *
     &                    TCVV(IDTNOX)         / AD(II,JJ,LL)      *
     &                    GOOD(II)
                     
                     ! Archive points in STT_TEMP51 (rvm, bmy, 2/27/02)
                     STT_TEMP51(I,J,L,N) = STT_TEMP51(I,J,L,N) + XX


                  ENDDO
                  ENDDO
                  ENDDO
               ENDIF
 
            !========================================================
            ! Archive PBL HEIGHTS [m]
            !========================================================
            CASE ( 37 )

               DO J = 1, NJ
                  JJ = JSET + J
               DO I = 1, NI
                  II = ISET + I  

                  ! Temp variable
                  TEMPBL = 0.d0

                  ! Get full boxes & add on box height.  For those, find 
                  ! fraction of highest box & add on corresponding box 
                  ! height for that.  Store total & number of entries to 
                  ! calculate average in diag3.f
                  PBLINT = FLOOR( XTRA2(II,JJ) )
                  PBLDEC = MOD( XTRA2(II,JJ),1d0 )

                  IF ( PBLINT > 0 ) THEN
                     DO M = 1, PBLINT
                        TEMPBL = TEMPBL + BXHEIGHT(II,JJ,M)
                     ENDDO
                  ENDIF

                  TEMPBL = TEMPBL + 
     &                     ( PBLDEC * BXHEIGHT(II,JJ,PBLINT+1) )

                  ! Save afternoon points in BLARR
                  ! Save count of afternoon points in PBLTOT
                  BLARR(I,J)  = BLARR(I,J)  + ( TEMPBL * GOOD(II) )
                  PBLTOT(I,J) = PBLTOT(I,J) + ( GOOD(II) )
               ENDDO
               ENDDO

            !===========================================================
            ! Archive PBL HEIGHTS [model levels]
            ! This is stored in XTRA2(i,j,5) (amf, 8/18/00)
            !===========================================================
            CASE ( 38 )

               DO J = 1, NJ
                  JJ = JSET + J
               DO I = 1, NI
                  II = ISET + I  

                  ! Archive into accumulating and counter arrays
                  BLLEVS(I,J) = BLLEVS(I,J) + XTRA2(II,JJ)
                  BLEVCT(I,J) = BLEVCT(I,J) + 1
               ENDDO
               ENDDO

            !===========================================================
            ! Archive AIR DENSITY [molec/cm3]
            !===========================================================
            CASE ( 39 )

               IF ( ALLOCATED( AIRDEN ) ) THEN
                  DO L = 1, NL
                     LL = LSET + L
                  DO J = 1, NJ
                     JJ = JSET + J
                  DO I = 1, NI
                     II = ISET + I  

                     ! Archive the afternoon points into STT_TEMP51
                     STT_TEMP51(I,J,L,N) = STT_TEMP51(I,J,L,N) +
     &                    ( AIRDEN(LL,II,JJ) *XNUMOLAIR*1d-6*GOOD(II) )
                  ENDDO
                  ENDDO
                  ENDDO
               ENDIF

            !===========================================================
            ! Archive CLOUD FRACTION [unitless] 
            ! Code taken from make_cldfrc.f
            ! NOTE: This has to be done outside the I-J-L loop above
            !===========================================================
            CASE ( 42 )

               ! Make sure cloud fraction arrays are allocated
               IF ( ALLOCATED( CLMOSW )  .and. 
     &              ALLOCATED( CLROSW ) ) THEN
                  DO J = 1, NJ
                     JJ = JSET + J
                  DO I = 1, NI
                     II = ISET + I

                     C1 = CLMOSW(1,II,JJ)       
                     C2 = 1.0d0 - CLROSW(1,II,JJ) 

                     DO L = 2, NL
                        LL = LSET + L
                        
                        IF ( CLMOSW(LL,II,JJ) > 
     &                       CLMOSW(LL-1,II,JJ) ) THEN
                           C1 = CLMOSW(LL,II,JJ)
                        ENDIF

                        C2 = C2 * ( 1.0d0 - CLROSW(LL,II,JJ) )
                     ENDDO

                     C1 = 1.0d0 - C1
                     
                     STT_TEMP51(I,J,1,N) = STT_TEMP51(I,J,1,N) + 
     &                    ( 1.0d0 - (C1 * C2) ) * GOOD(II)
                  ENDDO
                  ENDDO
               ENDIF

            !===========================================================
            ! Archive COLUMN OPTICAL DEPTHS [unitless]
            !===========================================================
            CASE ( 43 )

               ! Make sure OPTD array has been allocated
               IF ( ALLOCATED( OPTD ) ) THEN
                  DO L = 1, NL
                     LL = LSET + L
                  DO J = 1, NJ
                     JJ = JSET + J
                  DO I = 1, NI
                     II = ISET + I

                     ! Store afternoon points in STT_TEMP51
                     STT_TEMP51(I,J,1,N) = STT_TEMP51(I,J,1,N) + 
     &                                     ( OPTD(LL,II,JJ) * GOOD(II) )

                  ENDDO
                  ENDDO
                  ENDDO
               ENDIF

            !===========================================================
            ! Archive CLOUD TOP HEIGHTS [mb]
            !
            ! Assimilated meteorology is in the form of a model 
            ! level number. Need to convert this into a pressure.
            !===========================================================
            CASE ( 44 )
               
               ! Make sure CLDTOPS array has been allocated
               IF ( ALLOCATED( CLDTOPS ) ) THEN 
                  DO J = 1, NJ
                     JJ = JSET + J
                  DO I = 1, NI
                     II = ISET + I

                     ! Pressure at the cloud top
                     XX = GET_PEDGE( II, JJ, CLDTOPS(II,JJ) )

                     ! Store afternoon points in STT_TEMP51
                     STT_TEMP51(I,J,1,N) = STT_TEMP51(I,J,1,N) + 
     &                                     ( XX * GOOD(II) )
                  ENDDO
                  ENDDO
               ENDIF

            !===========================================================
            ! Archive PS - PTOP [mb]
            !===========================================================
            CASE ( 98 )

               DO J = 1, NJ
                  JJ = JSET + J
               DO I = 1, NI
                  II = ISET + I       

                  ! Archive the afternoon points in STT_TEMP51
                  STT_TEMP51(I,J,1,N) = STT_TEMP51(I,J,1,N) +
     &                 ( ( GET_PEDGE(II,JJ,1) - PTOP ) * GOOD(II) )

               ENDDO
               ENDDO

            !===========================================================
            ! Archive TEMPERATURE [K]
            !===========================================================
            CASE( 99 )

               DO L = 1, NL
                  LL = LSET + L
               DO J = 1, NJ
                  JJ = JSET + J
               DO I = 1, NI
                  II = ISET + I       

                  ! Archive the afternoon points in STT_TEMP51
                  STT_TEMP51(I,J,L,N) = STT_TEMP51(I,J,L,N) +
     &                                  ( T(II,JJ,LL) * GOOD(II) )
               ENDDO
               ENDDO
               ENDDO

            ! Everything else is an invalid tracer -- skip!
            CASE DEFAULT  
               CYCLE

         END SELECT

      ENDDO 
!$OMP END PARALLEL DO

      ! Return to calling program
      END SUBROUTINE ACCUMULATE_DIAG51

!------------------------------------------------------------------------------

      SUBROUTINE INIT_DIAG51
!
!******************************************************************************
!  Subroutine INIT_DIAG51 allocates and zeroes all module arrays.
!  (bmy, 11/29/00, 10/15/02)
! 
!  NOTES:
!  (1 ) Now references ALLOC_ERR from "error_mod.f" (bmy, 10/15/02)
!******************************************************************************
!    
      ! References to F90 modules
      USE ERROR_MOD, ONLY : ALLOC_ERR
  
#     include "CMN_SIZE"
#     include "CMN_TIMES"

      ! Local variables
      INTEGER :: AS

      !=================================================================
      ! INIT_DIAG51 begins here
      !=================================================================
      ALLOCATE( BLARR( NI, NJ ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'BLARR' )       
      BLARR(:,:) = 0d0

      ALLOCATE( BLLEVS( NI, NJ ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'BLLEVS' )
      BLLEVS(:,:) = 0

      ALLOCATE( BLEVCT( NI, NJ ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'BLEVCT' )
      BLEVCT(:,:) = 0

      ALLOCATE( GOOD( IIPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'GOOD' )
      GOOD(:) = 0

      ALLOCATE( PBLTOT( NI ,NJ ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'PBLTOT' )
      PBLTOT(:,:) = 0

      ALLOCATE( STT_TEMP51( NI, NJ, NL, NTRAC_AREA ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'STT_TEMP51' )
      STT_TEMP51(:,:,:,:) = 0d0

      ALLOCATE( WTOTH( NI, NJ ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'WTOTH' )
      WTOTH(:,:) = 0

      ALLOCATE( XLOCTM( IIPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'XLOCTM' )
      XLOCTM(:) = 0
         
      ! Return to calling program
      END SUBROUTINE INIT_DIAG51

!------------------------------------------------------------------------------

      SUBROUTINE CLEANUP_DIAG51
!
!******************************************************************************
!  Subroutine INIT_DIAG51 deallocates all module arrays. (bmy, 11/29/00)
!******************************************************************************
! 
      IF ( ALLOCATED( BLARR      ) ) DEALLOCATE( BLARR      )       
      IF ( ALLOCATED( BLLEVS     ) ) DEALLOCATE( BLLEVS     )
      IF ( ALLOCATED( BLEVCT     ) ) DEALLOCATE( BLEVCT     )
      IF ( ALLOCATED( GOOD       ) ) DEALLOCATE( GOOD       )
      IF ( ALLOCATED( PBLTOT     ) ) DEALLOCATE( PBLTOT     )
      IF ( ALLOCATED( STT_TEMP51 ) ) DEALLOCATE( STT_TEMP51 )
      IF ( ALLOCATED( WTOTH      ) ) DEALLOCATE( WTOTH      )
      IF ( ALLOCATED( XLOCTM     ) ) DEALLOCATE( XLOCTM     )

      ! Return to calling program
      END SUBROUTINE CLEANUP_DIAG51

!------------------------------------------------------------------------------

      END MODULE DIAG51_MOD







