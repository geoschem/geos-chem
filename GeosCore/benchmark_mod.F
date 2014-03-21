!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: benchmark_mod
!
! !DESCRIPTION: Module BENCHMARK\_MOD contains routines to save out initial 
!  and final tracer masses which are needed for GEOS-Chem 1-month benchmark 
!  simulations
!\\
!\\
! !INTERFACE: 
!
      MODULE BENCHMARK_MOD
!
! !USES:
!
      IMPLICIT NONE
      PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
      PUBLIC :: STDRUN
!
! !PUBLIC DATA MEMBERS:
!
      CHARACTER(LEN=255), PUBLIC :: INITIAL_FILE ! File w/ initial tracer mass
      CHARACTER(LEN=255), PUBLIC :: FINAL_FILE   ! File w/ final   tracer mass
!
! !REVISION HISTORY:
!  (1 ) Now expand date & time tokens in filenames (bmy, 1/31/05)
!  (2 ) Now modified for GCAP grid (swu, bmy, 6/28/05)
!  (3 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  01 Aug 2012 - R. Yantosca - Add reference to findFreeLUN from inqure_mod.F90
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
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
! !ROUTINE: stdrun
!
! !DESCRIPTION: Subroutine STDRUN dumps the mass of either O3 [kg] or 222Rn, 
!  210Pb, and 7Be [kg] at the start \& end of each run.  This is necessary
!  for GEOS-CHEM benchmarking.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE STDRUN( Input_Opt, State_Chm, LBEGIN )
!
! !USES:
!
      USE BPCH2_MOD,          ONLY : BPCH2             
      USE BPCH2_MOD,          ONLY : OPEN_BPCH2_FOR_WRITE
      USE BPCH2_MOD,          ONLY : GET_HALFPOLAR
      USE BPCH2_MOD,          ONLY : GET_MODELNAME
      USE CMN_SIZE_MOD
      USE FILE_MOD,           ONLY : IOERROR
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
      USE GIGC_State_Chm_Mod, ONLY : ChmState
      USE inquireMod,         ONLY : findFreeLUN
      USE TIME_MOD,           ONLY : EXPAND_DATE
      USE TIME_MOD,           ONLY : GET_NYMD
      USE TIME_MOD,           ONLY : GET_NHMS
      USE TIME_MOD,           ONLY : GET_TAU
      USE TRACERID_MOD,       ONLY : IDTO3
!
! !INPUT PARAMETERS: 
!
      TYPE(OptInput), INTENT(IN) :: Input_Opt   ! Input Options object
      TYPE(ChmState), INTENT(IN) :: State_Chm   ! Chemistry State object
      LOGICAL,        INTENT(IN) :: LBEGIN      ! =T denotes start of run
                                                ! =F denotes end of run
! 
! !REVISION HISTORY: 
!  12 Aug 2002 - R. Yantosca - Initial version
!  03 Aug 2012 - R. Yantosca - Added ProTeX headers
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
!  03 Aug 2012 - R. Yantosca - Now use findFreeLUN to obtain file unit #'s
!  14 Mar 2013 - M. Payer    - Replace Ox with O3 as part of removal of 
!                              NOx-Ox partitioning
!  25 Mar 2013 - M. Payer    - Now pass State_Chm object via the arg list
!  30 May 2013 - R. Yantosca - Now pass Input_Opt object via the arg list
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER             :: IU_FILE
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

      ! Pointers
      ! We need to define local arrays to hold corresponding values 
      ! from the Chemistry State (State_Chm) object. (mpayer, 12/6/12)
      REAL*8, POINTER :: STT(:,:,:,:)

      !=================================================================
      ! STDRUN begins here!
      !=================================================================

      ! Return if we are not doing either a radon or fullchem stdrun
      IF ( ( .not. Input_Opt%ITS_A_FULLCHEM_SIM )  .and. 
     &     ( .not. Input_Opt%ITS_A_RnPbBe_SIM   ) ) RETURN

      ! Find a free file LUN
      IU_FILE   = findFreeLUN()

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

      ! Initialize GEOS-Chem tracer array [kg] from Chemistry State object
      ! (mpayer, 12/6/12)
      STT       => State_Chm%Tracers

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
      IF ( Input_Opt%ITS_A_RnPbBE_SIM ) THEN

         ! Open binary punch file for writing
         CALL OPEN_BPCH2_FOR_WRITE( IU_FILE, FILENAME, TITLE )

         ! Loop over tracers
         DO N = 1, Input_Opt%N_TRACERS

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
      ! Save the mass of O3 to a file
      !=================================================================
      ELSE IF ( Input_Opt%ITS_A_FULLCHEM_SIM .and. IDTO3 > 0 ) THEN

         ! Open binary punch file for writing
         CALL OPEN_BPCH2_FOR_WRITE( IU_FILE, FILENAME, TITLE )
        
         ! Save O3 as REAL*4
         ARRAY(:,:,:) = STT(:,:,:,IDTO3)

         ! Write O3 to binary punch file
         CALL BPCH2( IU_FILE,   MODELNAME, LONRES,    LATRES,
     &               HALFPOLAR, CENTER180, CATEGORY,  IDTO3,
     &               UNIT,      GET_TAU(), GET_TAU(), RESERVED,
     &               IIPAR,     JJPAR,     LLPAR,     IFIRST,
     &               JFIRST,    LFIRST,    ARRAY(:,:,:) )
               
      ENDIF

      ! Close file
      CLOSE( IU_FILE )

      ! Free pointer
      NULLIFY( STT )

      END SUBROUTINE STDRUN
!EOC
      END MODULE BENCHMARK_MOD
