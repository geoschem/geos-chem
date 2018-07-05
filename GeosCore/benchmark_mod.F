!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: benchmark_mod.F
!
! !DESCRIPTION: Module BENCHMARK\_MOD contains routines to save out initial 
!  and final species masses which are needed for GEOS-Chem 1-month benchmark
!  simulations
!\\
!\\
! !INTERFACE: 
!
      MODULE BENCHMARK_MOD
!
! !USES:
!
      USE PRECISION_MOD

      IMPLICIT NONE
      PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
      PUBLIC :: STDRUN
!
! !PUBLIC DATA MEMBERS:
!
!
! !REVISION HISTORY:
!  (1 ) Now expand date & time tokens in filenames (bmy, 1/31/05)
!  (2 ) Now modified for GCAP grid (swu, bmy, 6/28/05)
!  (3 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  01 Aug 2012 - R. Yantosca - Add reference to findFreeLUN from inqure_mod.F90
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
!  05 Nov 2014 - M. Yannetti - PRECISION_MOD Changed REAL*8 to REAL(fp)
!  19 Mar 2015 - E. Lundgren - Change tracer units from kg to kg/kg
!  29 Nov 2016 - M. Sulprizio- Convert species locally from kg/kg to kg
!  30 Jan 2018 - M. Sulprizio- Remove INITIAL_FILE and FINAL_FILE and get from
!                              Input_Opt instead
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
      SUBROUTINE STDRUN( am_I_Root, Input_Opt,  State_Met, 
     &                   State_Chm, State_Diag, RC,        LBEGIN )
!
! !USES:
!
#if defined( BPCH_DIAG )
      USE BPCH2_MOD,      ONLY : BPCH2
      USE BPCH2_MOD,      ONLY : OPEN_BPCH2_FOR_WRITE
      USE BPCH2_MOD,      ONLY : GET_HALFPOLAR
      USE BPCH2_MOD,      ONLY : GET_MODELNAME
#endif
      USE CMN_SIZE_MOD
      USE ErrCode_Mod
      USE FILE_MOD,       ONLY : IOERROR
      USE Input_Opt_Mod,  ONLY : OptInput
      USE inquireMod,     ONLY : findFreeLUN
      USE State_Chm_Mod,  ONLY : ChmState
      USE State_Chm_Mod,  ONLY : Ind_
      USE State_Diag_Mod, ONLY : DgnState
      USE State_Met_Mod,  ONLY : MetState
      USE TIME_MOD,       ONLY : EXPAND_DATE
      USE TIME_MOD,       ONLY : GET_NYMD
      USE TIME_MOD,       ONLY : GET_NHMS
      USE TIME_MOD,       ONLY : GET_TAU
      USE UnitConv_Mod,   ONLY : Convert_Spc_Units
!
! !INPUT PARAMETERS: 
!
      LOGICAL,        INTENT(IN)    :: am_I_Root   ! Is this the root CPU?
      TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
      TYPE(MetState), INTENT(IN)    :: State_Met   ! Meteorology State object
      LOGICAL,        INTENT(IN)    :: LBEGIN      ! =T denotes start of run
                                                   ! =F denotes end of run
!
! !INPUT/OUTPUT PARAMETERS:
!
      TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
      TYPE(DgnState), INTENT(INOUT) :: State_Diag  ! Diagnostics State object
!
! !OUTPUT PARAMETERS:
!
      INTEGER,        INTENT(OUT)   :: RC          ! Success or failure

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
!  17 Dec 2014 - R. Yantosca - Leave time/date variables as 8-byte
!  17 Dec 2014 - R. Yantosca - Now use State_Chm%TRACERS instead of STT
!  19 Mar 2015 - E. Lundgren - Change tracer units from kg to kg/kg
!  16 Jun 2016 - R. Yantosca - Now use Ind_ to get species ID for O3
!  30 Jun 2016 - R. Yantosca - Remove instances of STT.  Now get the advected
!                              species ID from State_Chm%Map_Advect.
!  11 Aug 2016 - R. Yantosca - Remove temporary tracer-removal code
!  29 Nov 2016 - M. Sulprizio- Convert species locally from kg/kg to kg
!  26 Jun 2017 - R. Yantosca - GC_ERROR is now contained in errcode_mod.F90
!  28 Sep 2017 - E. Lundgren - Simplify unit conversions using wrapper routine
!  03 Nov 2017 - R. Yantosca - Now accept State_Diag as an argument
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER             :: IU_FILE
      INTEGER             :: id_O3,    NA,       nAdvect
      INTEGER             :: N,        NYMD,     NHMS
      INTEGER, PARAMETER  :: IFIRST=1, JFIRST=1, LFIRST=1
      INTEGER, PARAMETER  :: CENTER180=1
      INTEGER             :: HALFPOLAR
      REAL(f4)            :: ARRAY(IIPAR,JJPAR,LLPAR)
      REAL(f4)            :: LONRES, LATRES
      REAL(f8)            :: TAU
      CHARACTER(LEN=20)   :: MODELNAME 
      CHARACTER(LEN=40)   :: CATEGORY, RESERVED, UNIT
      CHARACTER(LEN=80)   :: TITLE
      CHARACTER(LEN=255)  :: FILENAME
      CHARACTER(LEN=63)   :: OrigUnit

      !=================================================================
      ! STDRUN begins here!
      !=================================================================

      ! Assume success
      RC = GC_SUCCESS

#if defined( BPCH_DIAG )

      ! Number of advected species
      nAdvect = State_Chm%nAdvect

      ! Return if we are not doing either a radon, tagO3, or fullchem stdrun
      IF ( ( .not. Input_Opt%ITS_A_FULLCHEM_SIM )  .and. 
     &     ( .not. Input_Opt%ITS_A_TAGO3_SIM    )  .and.
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

      ! Get species ID for O3
      id_O3     = Ind_('O3')

      ! Define filename for beginning or end of benchmark run
      IF ( LBEGIN ) THEN
         TITLE    = 'GEOS-CHEM Benchmark: Initial Species Mass'
         FILENAME = Input_Opt%STDRUN_INIT_FILE
      ELSE
         TITLE    = 'GEOS-CHEM Benchmark: Final Species Mass'
         FILENAME = Input_Opt%STDRUN_FINAL_FILE
      ENDIF
           
      ! Expand any date tokens in the filename
      CALL EXPAND_DATE( FILENAME, NYMD, NHMS )

      ! Convert species from to [kg]
      CALL Convert_Spc_Units( am_I_Root, Input_Opt, State_Met, 
     &                        State_Chm, 'kg', RC, OrigUnit=OrigUnit )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error('Unit conversion error', RC,
     &                 'Start of STDRUN in benchmark_mod.F')
         RETURN
      ENDIF

      !=================================================================
      ! Save the mass of 222Rn, 210Pb, 7Be to a file
      !=================================================================
      IF ( Input_Opt%ITS_A_RnPbBE_SIM ) THEN

         ! Open binary punch file for writing
         CALL OPEN_BPCH2_FOR_WRITE( IU_FILE, FILENAME, TITLE )

         ! Loop over only the advected species
         DO NA = 1, nAdvect

            ! Get the species ID from the advected species ID
            N = State_Chm%Map_Advect(NA)

            ! Save Rn, Pb, Be as REAL*4
            ARRAY(:,:,:) = State_Chm%Species(:,:,:,N)

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
      ELSE IF ( ( Input_Opt%ITS_A_FULLCHEM_SIM   .or. 
     &            Input_OPt%ITS_A_TAGO3_SIM    ) .and. id_O3 > 0 ) THEN

         ! Open binary punch file for writing
         CALL OPEN_BPCH2_FOR_WRITE( IU_FILE, FILENAME, TITLE )
        
         ! Save O3 as REAL*4
         ARRAY(:,:,:) = State_Chm%Species(:,:,:,id_O3)

         ! Write O3 to binary punch file
         CALL BPCH2( IU_FILE,   MODELNAME, LONRES,    LATRES,
     &               HALFPOLAR, CENTER180, CATEGORY,  id_O3,
     &               UNIT,      TAU,       TAU,       RESERVED,
     &               IIPAR,     JJPAR,     LLPAR,     IFIRST,
     &               JFIRST,    LFIRST,    ARRAY(:,:,:) )
               
      ENDIF

      ! Close file
      CLOSE( IU_FILE )

      ! Convert species back to original units
      CALL Convert_Spc_Units( am_I_Root, Input_Opt, State_Met, 
     &                        State_Chm, OrigUnit, RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error('Unit conversion error', RC,
     &                 'End of STDRUN in benchmark_mod.F')
         RETURN
      ENDIF
#endif

      END SUBROUTINE STDRUN
!EOC
      END MODULE BENCHMARK_MOD
