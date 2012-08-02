! $Id: benchmark_mod.f,v 1.1 2009/09/16 14:06:41 bmy Exp $
      MODULE BENCHMARK_MOD
!
!******************************************************************************
!  Module BENCHMARK_MOD contains routines to save out initial and final
!  tracer masses which are needed for GEOS-CHEM benchmark diagnostics.
!  (bmy, 7/20/04, 10/3/05)
!
!  Module Variables:
!  ============================================================================
!  (1 ) INITIAL_FILE (CHAR*255) : Name of file w/ initial tracer mass
!  (2 ) FINAL_FILE   (CHAR*255) : Name of file w/ final tracer mass
!
!  Module Routines:
!  ============================================================================
!  (1 ) STDRUN : Saves initial or final tracer mass to bpch file format
!
!  GEOS-CHEM modules referenced by biomass_mod.f
!  ============================================================================
!  (1 ) bpch2_mod.f     : Module containing routines for binary punch file I/O
!  (2 ) file_mod.f          : Contains file unit numbers and error checks
!  (3 ) logical_mod.f   : Module containing GEOS-CHEM logical switches
!  (4 ) time_mod.f      : Module containing routines for computing time & date
!  (5 ) tracer_mod.f    : Module containing GEOS-CHEM tracer array STT etc.
!  (6 ) tracerid_mod.f  : Module containing pointers to tracers & emissions
!
!  NOTES:
!  (1 ) Now expand date & time tokens in filenames (bmy, 1/31/05)
!  (2 ) Now modified for GCAP grid (swu, bmy, 6/28/05)
!  (3 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  01 Aug 2012 - R. Yantosca - Add reference to findFreeLUN from inqure_mod.F90
!******************************************************************************
!
      IMPLICIT NONE
#     include "define.h"

      !=================================================================
      ! MODULE VARIABLES
      !=================================================================
      CHARACTER(LEN=255) :: INITIAL_FILE
      CHARACTER(LEN=255) :: FINAL_FILE

      !=================================================================
      ! MODULE ROUTINES -- Follow below the "CONTAINS" statement
      !=================================================================
      CONTAINS

!------------------------------------------------------------------------------

      SUBROUTINE STDRUN( LBEGIN )
!
!******************************************************************************
!  Subroutine STDRUN dumps the mass of either Ox [kg] or 222Rn, 210Pb, and 7Be
!  [kg] at the start & end of each run.  This is necessary for GEOS-CHEM
!  benchmarking.  (bmy, 8/12/02, 10/3/05)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) LBEGIN (LOGICAL) : TRUE  denotes beginning of the run;
!                          FALSE denotes the end of the run
!
!  NOTES:
!  (1 ) Changed name from STDRUN_Ox to STDRUN, since we now can also save out 
!        Rn/Pb/Be for NSRCX==1.  Also deleted obsolete code from 6/02.  Added 
!        LBEGIN as an argument to determine if this is the start or end of the 
!        run.  (bmy, 8/12/02)
!  (2 ) Bundled into "benchmark_mod.f" (bmy, 7/20/04)
!  (3 ) Now expand date tokens in the filename (bmy, 1/31/05)
!  (4 ) Now call GET_HALFPOLAR from "bpch2_mod.f" to get the HALFPOLAR flag
!        value for GEOS or GCAP grids .  Also removed references to CMN_DIAG
!        and TRCOFFSET. (bmy, 6/28/05)
!  (5 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  01 Aug 2012 - R. Yantosca - Add reference to findFreeLUN from inqure_mod.F90
!******************************************************************************
!
      ! References to F90 modules
      USE BPCH2_MOD,    ONLY : OPEN_BPCH2_FOR_WRITE, BPCH2             
      USE BPCH2_MOD,    ONLY : GET_HALFPOLAR,        GET_MODELNAME
      USE FILE_MOD,     ONLY : IOERROR
      USE TIME_MOD,     ONLY : EXPAND_DATE,          GET_NYMD
      USE TIME_MOD,     ONLY : GET_NHMS,             GET_TAU
      USE TRACER_MOD,   ONLY : ITS_A_FULLCHEM_SIM,   ITS_A_RnPbBe_SIM
      USE TRACER_MOD,   ONLY : STT,                  N_TRACERS
      USE TRACERID_MOD, ONLY : IDTOX
      USE CMN_SIZE_MOD       ! Size parameters
#if   defined( DEVEL )
      USE inquireMod,   ONLY : findFreeLUN
#else
      USE FILE_MOD,     ONLY : IU_FILE
#endif

      ! Arguments
      LOGICAL, INTENT(IN) :: LBEGIN 

      ! Local variables
#if   defined( DEVEL )
      INTEGER             :: IU_FILE
#endif
      INTEGER             :: N,        NYMD,     NHMS
      INTEGER, PARAMETER  :: IFIRST=1, JFIRST=1, LFIRST=1
      INTEGER, PARAMETER  :: CENTER180=1
      INTEGER             :: HALFPOLAR
      REAL*4              :: ARRAY(IIPAR,JJPAR,LLPAR)
      REAL*4              :: LONRES, LATRES
      REAL*8              :: TAU
      CHARACTER(LEN=20)   :: MODELNAME 
      CHARACTER(LEN=40)   :: CATEGORY, RESERVED, UNIT
      CHARACTER(LEN=80)   :: TITLE
      CHARACTER(LEN=255)  :: FILENAME

      !=================================================================
      ! STDRUN begins here!
      !=================================================================

      ! Return if we are not doing either a radon or fullchem stdrun
      IF ( ( .not. ITS_A_FULLCHEM_SIM() ) .and. 
     &     ( .not. ITS_A_RnPbBe_SIM() ) ) RETURN

#if   defined( DEVEL )
      ! Get the next free file unit
      IU_FILE   = findFreeLUN()
#endif

      ! Define variables for binary punch file
      MODELNAME = GET_MODELNAME()
      HALFPOLAR = GET_HALFPOLAR()
      CATEGORY  = 'TCMASS-$'
      UNIT      = 'kg'
      RESERVED  = ''      
      LONRES    = DISIZE
      LATRES    = DJSIZE
      NYMD      = GET_NYMD()
      NHMS      = GET_NHMS()
      TAU       = GET_TAU()

      ! Define filename for beginning or end of benchmark run
      IF ( LBEGIN ) THEN
         TITLE    = 'GEOS-CHEM Benchmark: Initial Tracer Mass'
         FILENAME = INITIAL_FILE
      ELSE
         TITLE    = 'GEOS-CHEM Benchmark: Final Tracer Mass'
         FILENAME = FINAL_FILE
      ENDIF
           
      ! Expand any date tokens in the filename
      CALL EXPAND_DATE( FILENAME, NYMD, NHMS )

      !=================================================================
      ! Save the mass of 222Rn, 210Pb, 7Be to a file
      !=================================================================
      IF ( ITS_A_RnPbBE_SIM() ) THEN

         ! Open binary punch file for writing
         CALL OPEN_BPCH2_FOR_WRITE( IU_FILE, FILENAME, TITLE )

         ! Loop over tracers
         DO N = 1, N_TRACERS

            ! Save Rn, Pb, Be as REAL*4
            ARRAY(:,:,:) = STT(:,:,:,N)

            ! Write Rn, Pb, Be to binary punch file
            CALL BPCH2( IU_FILE,   MODELNAME, LONRES,    LATRES,
     &                  HALFPOLAR, CENTER180, CATEGORY,  N,    
     &                  UNIT,      TAU,       TAU,       RESERVED,   
     &                  IIPAR,     JJPAR,     LLPAR,     IFIRST,     
     &                  JFIRST,    LFIRST,    ARRAY(:,:,:) )

         ENDDO

      !=================================================================
      ! Save the mass of Ox to a file
      !=================================================================
      ELSE IF ( ITS_A_FULLCHEM_SIM() .and. IDTOX > 0 ) THEN

         ! Open binary punch file for writing
         CALL OPEN_BPCH2_FOR_WRITE( IU_FILE, FILENAME, TITLE )
        
         ! Save Ox as REAL*4
         ARRAY(:,:,:) = STT(:,:,:,IDTOX)

         ! Write Ox to binary punch file
         CALL BPCH2( IU_FILE,   MODELNAME, LONRES,    LATRES,
     &               HALFPOLAR, CENTER180, CATEGORY,  IDTOX,    
     &               UNIT,      GET_TAU(), GET_TAU(), RESERVED,   
     &               IIPAR,     JJPAR,     LLPAR,     IFIRST,     
     &               JFIRST,    LFIRST,    ARRAY(:,:,:) )
               
      ENDIF

      ! Close file
      CLOSE( IU_FILE )

      ! Return to MAIN program
      END SUBROUTINE STDRUN

!------------------------------------------------------------------------------

      ! End of module
      END MODULE BENCHMARK_MOD
