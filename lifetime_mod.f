      MODULE LIFETIME_MOD

!
!******************************************************************************
!  Module LIFETIME_MOD is written to caclulate the lifetimes of desired
!  species you select. The routine for retrieving the rates of the
!  reactions follows the machinary in place inside PLANEFLIGHT_MOD (the
!  subroutine ARCHIVE_RXNS_FOR_PF( JO1D, N2O5 ) which is called inside of
!  calcrate.f
!  (jpp, 8/20/08)
!
!  Module Variables:
!  ============================================================================
!
!******************************************************************************
!


      IMPLICIT NONE

      !=================================================================
      ! MODULE PRIVATE DECLARATIONS -- keep certain internal variables 
      ! and routines from being seen outside "planeflight_mod.f"
      !=================================================================

      ! Make everything PRIVATE ...
      PRIVATE

      ! ... except these routines
      PUBLIC :: INIT_LIFETIME
      PUBLIC :: SET_VARS_4LT
      PUBLIC :: ARCHIVE_RXNS_NEW
!      PUBLIC :: GET_CONC_4LT
      PUBLIC :: CLEANUP_LIFETIME
      PUBLIC :: save_rxn_rates
      ! jpp, debugging
      PUBLIC :: test_lt

      ! make these variables PUBLIC
      ! for access from diag3.f
!      PUBLIC :: TAU_LT
!      PUBLIC :: NTRA_4LT
!      PUBLIC :: NRXN_4LT
!      PUBLIC :: NLOSS
      PUBLIC :: LOSS_NUM
      PUBLIC :: rate_const_count
      !=================================================================
      ! MODULE VARIABLES 
      !=================================================================

      ! Parameters
      INTEGER,           PARAMETER   :: MAXVARS   = 95
      INTEGER,           PARAMETER   :: MAXREAC   = 50
      INTEGER,           PARAMETER   :: MAXRO2    = 20

      ! For specifying date/time
!     INTEGER,           ALLOCATABLE :: PDATE(:)
!     INTEGER,           ALLOCATABLE :: PTIME(:)
!     REAL*4,            ALLOCATABLE :: PTAU(:)

      ! jpp: hardwiring for now... read in later
      ! For specifying how many reaction rates and how
      ! many tracers to keep track of in calculating
      ! lifetimes:
!      INTEGER, PARAMETER             :: NRXN_4LT = 24
!      INTEGER                        :: NPREAC = 24
      INTEGER, PARAMETER             :: NRXN_4LT = 32 ! jpp replaced 26 ! jpp, replaced 25, 3/16/2010
      INTEGER                        :: NPREAC = 32 ! jpp replaced 26 (3/8/2011)
      INTEGER, PARAMETER             :: NTRA_4LT = 10

      ! using RXNVAR to store the indices for reaction rates
      INTEGER,           ALLOCATABLE :: RXNVAR(:)
      ! using TRAVAR to store the indices for tracer #s
      INTEGER,           ALLOCATABLE :: TRAVAR(:)
!      CHARACTER(LEN=10), ALLOCATABLE :: PNAME(:)

      ! For specifying SMVGEAR rxns to save at each flight point

      INTEGER,     ALLOCATABLE :: PREAC(:) 
      REAL*8,      ALLOCATABLE :: PRRATE(:,:,:,:) ! rxn constant
      REAL*8,      ALLOCATABLE :: RXN_RATE(:,:)
      ! jpp: count the number of times you store a new
      !      reaction rate, so that you can average
      !      properly
      INTEGER,  allocatable :: rate_const_count(:,:,:,:)
      INTEGER,           ALLOCATABLE :: RXN_COUNT(:,:,:)
      INTEGER,           ALLOCATABLE :: SPC_COUNT(:,:,:)
      INTEGER,    ALLOCATABLE :: SP2_COUNT(:,:)
      ! jpp: saving the mass:
      REAL*8,            ALLOCATABLE :: MASS(:,:)
      ! jpp: saving the concentrations,
      !      averaged over time:
      REAL*8,            ALLOCATABLE :: CONCS(:,:)
      REAL*8,            ALLOCATABLE :: TAU_LT(:,:,:,:,:)

      ! the molecular weights:
      REAL*8  :: MWEIGHT(NTRA_4LT)

      ! loop integer
      INTEGER  :: KKLOOP

      ! logical mask for collapsing TAU_LT(I,J,L,NRXN,NTRA)
      ! into the AD(I,J,L,SUM(NLOSS))
      LOGICAL, SAVE  :: D52MASK(NRXN_4LT, NTRA_4LT) = .FALSE.

      ! Logical for first call to the program
      LOGICAL, SAVE  :: IF_FIRST = .TRUE.

      ! logical hardwired to TRUE, allowing the program
      ! to run and store reactions. Eventually make it
      ! a common variable that can be selected within
      ! the input.geos file.
      LOGICAL, SAVE  :: DO_LT = .TRUE.

      ! tracers used here for reaction rates
      REAL*8, ALLOCATABLE  :: O3_lt(:), HO2_lt(:), NO3_lt(:)
      REAL*8, ALLOCATABLE  :: CH2O_lt(:), NO2_lt(:), NO_lt(:)
      REAL*8, ALLOCATABLE  :: OH_lt(:)


      ! now load the molecular weights (ordered as in input.geos)
      data MWEIGHT/.160d0,80.d-3,96.d-3,97.d-3,8.d-3,
     &     126.d-3,142.d-3,253.d-3,174.d-3,95.d-3/

      ! Integer for storing the number of loss rates for
      ! each chemical species
      INTEGER :: NLOSS(NTRA_4LT)
      data NLOSS/2,6,7,1,1,1,2,1,1,2/
      ! jpp, hardwired this... didn't work inside allocation
!      INTEGER :: LOSS_NUM(26) ! jpp, 3/16/2010: added space for BrNO3 hydrolysis
      INTEGER :: LOSS_NUM(32) ! jpp, added 4 new spots for the
                              !      Br + VOCs or hydrocarbons
                              ! and + 2 for HOBr + HBr + aer. psuedo rxns

      ! jpp, 3/8/2011: add space for 4 additional reactions
      !                some extra Br + hydrocarbons and VOCs
      data LOSS_NUM/301,374,    ! Br2 loss rxns
     &     293,295,300,303,304,305, 306, 307, 308, 309,  ! Br
     &     294,297,298,299,302,311,375, ! BrO
     &     317,376, !HOBr !includes the HOBr + HBr pseudo reaction
     &     296,318, !HBr  !includes the HOBr + HBr pseudo reaction
     &     379, !BrNO2
     &     300, 316, 377, 378, !BrNO3
     &     313, 380, !CHBr3
     &     314, !CH2Br2
     &     315/ !CH3Br

!jpt      data LOSS_NUM/304,375,            ! Br2 loss rxns
!jpt     &     296,298,303,306,307,308, 309, 310, 311, 312,  ! Br
!jpt     &     297,300,301,302,305,314,376, ! BrO
!jpt     &     377, !HOBr !includes the HOBr + HBr pseudo reaction
!jpt     &     299, !HBr  !includes the HOBr + HBr pseudo reaction
!jpt     &     380, !BrNO2
!jpt     &     303, 319, 378, 379, !BrNO3
!jpt     &     316, 381, !CHBr3
!jpt     &     317, !CH2Br2
!jpt     &     318/ !CH3Br


!      data LOSS_NUM/307,372,299,301,306,309,310,311,
!     &     300,303,304,305,308,313,373,374,302,377,306
!     &     375,376,315,378,316,317/

      !=================================================================
      ! MODULE ROUTINES -- follow below the "CONTAINS" statement 
      !=================================================================
      CONTAINS

!------------------------------------------------------------------------------

      SUBROUTINE GET_CONC_4LT( SUNCOS )
!
!******************************************************************************
! This subroutine returns the mass grid corresponding to the
! specified compounds. Right now it's just hard-wired
! for each run to generate the results. Eventually we might
! want to streamline this so that everyone can use it.
! (jpp, 5/21/08)
!
! ---------------------------------------------------------------
! NOTES:
!   (1 ) change so that conc's are in mass: use VOLUME(JLOOP) to
!        get the volume of the box. And then convert with Av and
!        the appropriate molar masses. this will weight lifetime
!        to represent the lower troposphere... molecules/cm3 will
!        weight to the upper troposphere. jpp, 7/30/08
!
!
!******************************************************************************
!

      ! References to F90 modules
      USE COMODE_MOD,     ONLY : JLOP, IXSAVE
      USE COMODE_MOD,     ONLY : IYSAVE, IZSAVE
      USE COMODE_MOD,     ONLY : CSPEC, VOLUME
      USE TROPOPAUSE_MOD, ONLY : ITS_IN_THE_TROP
      USE ERROR_MOD,      ONLY : GEOS_CHEM_STOP
      USE TRACER_MOD,     ONLY : N_TRACERS
      USE TRACERID_MOD,   ONLY : IDO3, IDHO2, IDOH
      USE TRACERID_MOD,   ONLY : IDNO, IDNO2, IDNO3
      USE TRACERID_MOD,   ONLY : IDCH2O

#     include "CMN_SIZE"  ! size parameters
#     include "comode.h"  ! for using NTLOOP, and avagadro's # (Avg)

      ! Arguments 
      REAL*8, INTENT(IN) :: SUNCOS( MAXIJ )

      ! Local Variables
      REAL*8          :: tnum ! temporary molecules/cm3
                              ! taken from CSPEC array
                              ! directly after gasconc.f
                              ! is called in chemdr.f
      ! (a) loop control
      INTEGER         :: I, J, L, JLOOP, N
      INTEGER         :: IJWINDOW
      ! (b) trick for previous count
      REAL*8     :: prev_ct

      ! species loop integer
      INTEGER         :: V
      ! jpp, debugging
      INTEGER :: COUNT

      ! reactant names
      character(len=10) :: rctnames(7), br_select(10)
      data rctnames/'O3','HO2','NO3','CH2O','NO2',
     &     'NO','OH'/
      data br_select/'Br2','Br','BrO','HOBr','HBr',
     &     'BrNO2','BrNO3','CHBr3','CH2Br2','CH3Br'/

      !------------------
      ! begin here
      !------------------
      print *, 'NTSPECGAS=', NTSPECGAS
      print *, IDO3
      print *, IDOH
      ! retrieve the necessary concentrations of 
      ! co-reactant tracers at the proper jloop
      ! locations. These are in units of molecules/cm3
      DO N = 1, NTSPECGAS  ! loop through the number of active gases


         DO JLOOP = 1, NTLOOP
            ! get box info
            I = IXSAVE( JLOOP )
            J = IYSAVE( JLOOP )
            L = IZSAVE( JLOOP )

            IJWINDOW  = (J-1)*IIPAR + I

            ! make sure it's in the trop and that
            ! it's in sunlight
            IF ( (SUNCOS(IJWINDOW) > 0) .AND. 
     &           (ITS_IN_THE_TROP(I,J,L)) ) THEN

               IF ( N == IDO3 ) THEN

                  ! start counting for averaging
                  SP2_COUNT(JLOOP, 1) = SP2_COUNT(JLOOP, 1) + 1

                  !previous count
                  prev_ct = DBLE( MAX0( 
     &                 SP2_COUNT(JLOOP, 1)-1, 1 ) )
                  ! get the concs
                  IF (SP2_COUNT(JLOOP,1) == 1) THEN
                     O3_lt(JLOOP)  = CSPEC(JLOOP,N)
                  ELSE
                     ! get concentrations in each box
                     ! and average them along the way
                     O3_lt(JLOOP) = O3_lt(JLOOP) *
     &                    prev_ct/ DBLE(SP2_COUNT(JLOOP,1)) +
     &                    CSPEC(JLOOP, N) /
     &                    DBLE(SP2_COUNT(JLOOP,1))
                  ENDIF

               ELSE IF ( N == IDHO2 ) THEN
!               CASE( 'HO2' )
                  ! start counting for averaging
                  SP2_COUNT(JLOOP, 2) = SP2_COUNT(JLOOP, 2) + 1

                  !previous count
                  prev_ct = DBLE( MAX0( 
     &                 SP2_COUNT(JLOOP, 2)-1, 1 ) )

                  ! get the concentrations
                  IF ( SP2_COUNT(JLOOP, 2) == 1 ) THEN
                     HO2_lt(JLOOP)  = CSPEC(JLOOP,N)
                  ELSE
                     ! get concentrations in each box
                     ! and average them along the way
                     HO2_lt(JLOOP) = HO2_lt(JLOOP) *
     &                    prev_ct/ DBLE(SP2_COUNT(JLOOP,2)) +
     &                    CSPEC(JLOOP, N) /
     &                    DBLE(SP2_COUNT(JLOOP,2))
                  ENDIF

               ELSE IF ( N == IDNO3 ) THEN
!               CASE( 'NO3' )
                  ! start counting for averaging
                  SP2_COUNT(JLOOP, 3) = SP2_COUNT(JLOOP, 3) + 1

                  !previous count
                  prev_ct = DBLE( MAX0( 
     &                 SP2_COUNT(JLOOP, 3)-1, 1 ) )

                  ! get the concentrations
                  IF ( SP2_COUNT(JLOOP, 3) == 1 ) THEN
                     NO3_lt(JLOOP)  = CSPEC(JLOOP,N)
                  ELSE
                     ! get concentrations in each box
                     ! and average them along the way
                     NO3_lt(JLOOP) = NO3_lt(JLOOP) *
     &                    prev_ct/ DBLE(SP2_COUNT(JLOOP,3)) +
     &                    CSPEC(JLOOP, N) /
     &                    DBLE(SP2_COUNT(JLOOP,3))
                  ENDIF

               ELSE IF ( N == IDCH2O ) THEN
!               CASE( 'CH2O' )
                  ! start counting for averaging
                  SP2_COUNT(JLOOP, 4) = SP2_COUNT(JLOOP, 4) + 1

                  !previous count
                  prev_ct = DBLE( MAX0( 
     &                 SP2_COUNT(JLOOP, 4)-1, 1 ) )

                  ! get the concentrations
                  IF ( SP2_COUNT(JLOOP, 4) == 1 ) THEN
                     CH2O_lt(JLOOP)  = CSPEC(JLOOP,N)
                  ELSE
                     ! get concentrations in each box
                     ! and average them along the way
                     CH2O_lt(JLOOP) = CH2O_lt(JLOOP) *
     &                    prev_ct/ DBLE(SP2_COUNT(JLOOP,4)) +
     &                    CSPEC(JLOOP, N) /
     &                    DBLE(SP2_COUNT(JLOOP,4))
                  ENDIF

               ELSE IF ( N == IDNO2 ) THEN
!               CASE( 'NO2' )
                  ! start counting for averaging
                  SP2_COUNT(JLOOP, 5) = SP2_COUNT(JLOOP, 5) + 1

                  !previous count
                  prev_ct = DBLE( MAX0( 
     &                 SP2_COUNT(JLOOP, 5)-1, 1 ) )

                  ! get the concentrations
                  IF ( SP2_COUNT(JLOOP, 5) == 1 ) THEN
                     NO2_lt(JLOOP)  = CSPEC(JLOOP,N)
                  ELSE
                     ! get concentrations in each box
                     ! and average them along the way
                     NO2_lt(JLOOP) = NO2_lt(JLOOP) *
     &                    prev_ct/ DBLE(SP2_COUNT(JLOOP,5)) +
     &                    CSPEC(JLOOP, N) /
     &                    DBLE(SP2_COUNT(JLOOP,5))
                  ENDIF

               ELSE IF ( N == IDNO ) THEN
!               CASE( 'NO' )
                  ! start counting for averaging
                  SP2_COUNT(JLOOP, 6) = SP2_COUNT(JLOOP, 6) + 1

                  !previous count
                  prev_ct = DBLE( MAX0( 
     &                 SP2_COUNT(JLOOP, 6)-1, 1 ) )

                  ! get the concentrations
                  IF ( SP2_COUNT(JLOOP, 6) == 1 ) THEN
                     NO_lt(JLOOP)  = CSPEC(JLOOP,N)
                  ELSE
                     ! get concentrations in each box
                     ! and average them along the way
                     NO_lt(JLOOP) = NO_lt(JLOOP) *
     &                    prev_ct/ DBLE(SP2_COUNT(JLOOP,6)) +
     &                    CSPEC(JLOOP, N) /
     &                    DBLE(SP2_COUNT(JLOOP,6))
                  ENDIF

               ELSE IF ( N == IDOH ) THEN
!               CASE( 'OH' )
                  ! start counting for averaging
                  SP2_COUNT(JLOOP, 7) = SP2_COUNT(JLOOP, 7) + 1

                  !previous count
                  prev_ct = DBLE( MAX0( 
     &                 SP2_COUNT(JLOOP, 7)-1, 1 ) )

                  ! get the concentrations
                  IF ( SP2_COUNT(JLOOP, 7) == 1 ) THEN
                     OH_lt(JLOOP)  = CSPEC(JLOOP,N)
                  ELSE
                     ! get concentrations in each box
                     ! and average them along the way
                     OH_lt(JLOOP) = OH_lt(JLOOP) *
     &                    prev_ct/ DBLE(SP2_COUNT(JLOOP,7)) +
     &                    CSPEC(JLOOP, N) /
     &                    DBLE(SP2_COUNT(JLOOP,7))
                  ENDIF

               ENDIF
!               END SELECT

            ENDIF
         ENDDO ! end loop over troposphere boxes
      ENDDO    ! end loop over species names





      ! find and store the masses for
      ! Bromine compounds
      DO N = 1, NTRA_4LT
      !jpp, debugging
      count=0

      DO JLOOP = 1, NTLOOP
         I = IXSAVE(JLOOP)
         J = IYSAVE(JLOOP)
         L = IZSAVE(JLOOP)
         IJWINDOW  = (J-1)*IIPAR + I

         ! make sure it's in the trop and that
         ! it's in sunlight
         IF ( (SUNCOS(IJWINDOW) > 0) .AND. 
     &        (ITS_IN_THE_TROP(I,J,L)) ) THEN

            ! begin the counting for subsequent
            ! averaging
            SPC_COUNT(I, J, L) =
     &           SPC_COUNT(I, J, L) + 1

            ! now store the previous count if
            ! we're beyond count = 1
            ! store it as double precision for
            ! calculating the average
            prev_ct = DBLE( MAX0( 
     &           SPC_COUNT(I, J, L)-1, 1 ) )

            ! avoid factoring 0 concentrations
            ! from the initialization of the CONCS
            ! array into the averaging scheme
            IF ( SPC_COUNT(I,J,L) == 1 ) THEN
               CONCS(JLOOP,N) = CSPEC(JLOOP, TRAVAR(N))
            ELSE
               ! get concentrations in each box
               ! and average them along the way
               CONCS(JLOOP,N) = CONCS(JLOOP, N) *
     &              prev_ct/ DBLE(SPC_COUNT(I,J,L)) +
     &              CSPEC(JLOOP, TRAVAR(N)) /
     &              DBLE(SPC_COUNT(I,J,L))
            ENDIF

            ! now calculate the mass: kg
            ! VOLUME() is in units cm3... CONCS is
            ! in units of molecules/cm3... so no conversion
            ! necessary
            MASS(JLOOP, N) = CONCS(JLOOP, N) * 
     &           VOLUME(JLOOP) / AVG * MWEIGHT(N)

!                  DO V = 1, NTRA_4LT
!                     ! adding up the masses
!                     ! from the different species
!                     ! of interest...e.g. Br + BrO = Brx 
!                     tnum = tnum + CSPEC(JLOOP, TRAVAR(V))
!                  ENDDO

         ENDIF

!      !jpp, debugging:
!      print *, 'N =', N
!      print *, 'count=',count
      ENDDO
      ENDDO

      ! jpp, debugging:
      print *, 'Max(BrO)=', MAXVAL(CONCS(:,3))
      print *, 'Min(BrO)=', MINVAL(CONCS(:,3))
      print *, 'llpar=', LLPAR
      print *, 'stopping run'
!      call flush(6)
!      call geos_chem_stop


      ! jpp, debugging
      ! still might be some problems... several
      ! species have ridiculously small maximum
      ! values for CSPEC, jpp 7/5/08
!      do N = 1, NTRA_4LT
!         print *, 'max(conc) = ', MAXVAL(CONCS(:,N))
!      enddo
!      print *,'jpp: worked... not an index problem'
!      call flush(6)
!      CALL GEOS_CHEM_STOP


      ! Return to calling program
      END SUBROUTINE GET_CONC_4LT

!------------------------------------------------------------------------------



      SUBROUTINE SET_VARS_4LT
!
!******************************************************************************
!  Subroutine READ_VARIABLES reads the list of variables (SMVGEAR species,
!  SMVGEAR rxn rates, GMAO met fields, or GEOS-Chem tracers) to be printed
!  out and sorts the information into the appropriate module variables.
!  (mje, bmy, 7/30/02, 10/16/06)
!
!  NOTES:
!  (1 ) Now references GEOS_CHEM_STOP from "error_mod.f", which frees all
!        allocated memory before stopping the run. (bmy, 10/15/02)
!  (2 ) Bug fix: replace missing commas in FORMAT statement (bmy, 3/23/03)
!  (3 ) Bug fix: replace NAMESPEC w/ NAMEGAS for SMVGEAR II (lyj, bmy, 7/9/09)
!  (4 ) Now locate reordered rxn numbers for SMVGEAR II. (mje, bmy, 8/1/03)
!  (5 ) Now flag N2O5 hydrolysis rxn as a special case (bmy, 8/8/03)
!  (6 ) Changed variable name prefix "DAO" to "GMAO".  Also added aerosol
!        optical depths w/ tracer offset 2000. (bmy, 4/23/04)
!  (7 ) Now references N_TRACERS & ITS_A_FULLCHEM_SIM from "tracer_mod.f"
!        (bmy, 7/20/04)
!  (8 ) Bug fix: extract tracer # when reading rxn rates (bmy, 1/7/05)
!  (9 ) Now computes column AOD's and AOD's below plane (bmy, 10/24/05)
!  (10) We need to trim NAMEGAS before comparing to LINE so that comparisons 
!        for species like "O3" will work.  Also set NCS=NCSURBAN at the top
!        of the subroutine, to avoid out of bounds error. (dbm, bmy, 10/16/06)
!******************************************************************************
!
      ! References to F90 modules
      USE ERROR_MOD,  ONLY : GEOS_CHEM_STOP
      USE FILE_MOD,   ONLY : IU_FILE,   IOERROR
      USE TRACER_MOD, ONLY : N_TRACERS, ITS_A_FULLCHEM_SIM

#     include "CMN_SIZE"  ! Size parameters
#     include "comode.h"  ! NAMEGAS, NSPEC, MAPPL
      ! jpp: MAPPL maps the old species number to
      !      the new chemistry solver one...

      ! Local variables
      LOGICAL             :: IS_FULLCHEM
      INTEGER             :: M, N, NUM, R, IK, IOS
      CHARACTER(LEN=255)  :: LINE
      !----------------------------------
      ! The following variables are
      ! simply hardwired for now:
      ! eventually have them read in from
      ! some other file (jpp, 5/21/08)
      !----------------------------------
      ! number of reactions to save... NRXN_4LT and NTRA_4LT defined
      ! above... so all local subroutines can use
!      INTEGER, PARAMETER  :: NRXN_4LT = 7 ! number of rxn rates
      CHARACTER(LEN=255)  :: RXNNAME(NRXN_4LT) ! the rxn #'s in smv2.log
      ! number of species concentrations desired
!      INTEGER, PARAMETER  :: NTRA_4LT = 2 ! Br + BrO = Brx
      CHARACTER(LEN=255)  :: TRANAME(NTRA_4LT) ! the species number in input.geos
      CHARACTER(LEN=255)  :: rntemp! temporary rxnname entry
      CHARACTER(LEN=255)  :: trtemp! temporary traname entry

      !=================================================================
      ! READ_VARIABLES begins here!
      !=================================================================

      ! Reset NCS to NCSURBAN for safety's sake (dbm, bmy, 10/16/06)
      NCS = NCSURBAN

      ! Test if this is a fullchem run
      IS_FULLCHEM = ITS_A_FULLCHEM_SIM()


      
      ! Zero reaction counter
      R = 0
      ! Loop over all the desired reactions
      DO N = 1, NRXN_4LT

            !===========================================================
            ! SMVGEAR rxn rate: listed as "REA_001", etc.
            ! PVAR offset: 10000
            !---------------------------------------------
            ! Get Reaction Rate #'s
            !===========================================================
            ! jpp: here is the reaction rate selection case

         ! Skip if not SMVGEAR!
         IF ( IS_FULLCHEM ) THEN 
               
            ! Increment rxn counter
            R = R + 1

            !==================================================
            ! NOTE: the reaction numbers listed in smv2.log 
            ! aren't really used to index SMVGEAR II rxns.  The 
            ! rxns get reordered.  Find the right rxn number, 
            ! which is stored in NOLDFNEW.  We assume only one 
            ! chemistry scheme. (mje, bmy, 8/1/03)
            !==================================================

            ! Extract tracer # from the string
            NUM = LOSS_NUM(N) ! take the proper rxn number

            ! Initialize
            RXNVAR(N)  = -999
            PREAC(R)   = -999

            ! Search for proper rxn number
            DO IK = 1, NMTRATE 

               ! Offset other reaction rates by 10000
               IF ( NOLDFNEW(IK,1) == NUM ) THEN 
                  RXNVAR(N)  = IK
                  PREAC(R)   = IK
                  EXIT
               ENDIF

            ENDDO

            ! Stop w/ error 
            IF ( RXNVAR(N) == -999 ) THEN 
               WRITE (6,*) 'Cant match up reaction number'
               WRITE (6,*) NUM
               WRITE (6,*) 'Is it the second line of the'
               WRITE (6,*) 'Three body reaction'
               WRITE (6,*) 'Stopping'
               CALL GEOS_CHEM_STOP
            ENDIF

         ENDIF


         ! Echo species names/numbers to screen
        WRITE( 6, *) ' inside READ_VARIABLES '
!     WRITE( 6, 120 ) N, TRIM( RXNNAME(N) ), RXNVAR(N)
        WRITE( 6, 120 ) N, '; preac = ', PREAC(R)
120     FORMAT( i4, 1x, a10, 1x, i10 )

      ENDDO

      !jpp debuggin:
      ! did it read all of the species
!      CALL FLUSH ( 6 )
!      CALL GEOS_CHEM_STOP



      RETURN
      ! Return to calling program
      END SUBROUTINE SET_VARS_4LT

!------------------------------------------------------------------------------

      SUBROUTINE ARCHIVE_RXNS_NEW( SUNCOS )

      ! References to F90 modules
      USE COMODE_MOD, ONLY     : IXSAVE, IYSAVE, IZSAVE
      USE TROPOPAUSE_MOD, ONLY : ITS_IN_THE_TROP
      USE TIME_MOD,       only : GET_LOCALTIME
      USE ERROR_MOD,  ONLY     : GEOS_CHEM_STOP
      USE DIAG_MOD,   ONLY     : LT_Br
      USE DIAG_MOD,       ONLY : AD53

#     include "CMN_SIZE"  ! Size parameters
C!!#     include "CMN_DIAG"  ! ND40 switch
#     include "comode.h"  ! TRATE,RRATE, JLOOPLO, KBLOOP
 
      ! Arguments 
      REAL*8, INTENT(IN) :: SUNCOS( MAXIJ )

      ! internal variables
      INTEGER :: KLOOP, JLOOP, I, J, L
      INTEGER :: IJWINDOW, R
      REAL*8  :: RATE
      real*8  :: loctime

      ! Smallest, largest REAL*4 #'s representable on this machine
      REAL*4, PARAMETER  :: SMALLEST=TINY(1e0), LARGEST=HUGE(1e0)



      ! loop over smgear reactions
      DO R = 1, NPREAC

         ! Store rate in PRRATE
         DO KLOOP = 1, KTLOOP
            ! LREORDER() array is now necessary for
            ! using the 1D array indexing in v.8-02-02. I
            ! am updating this for lifetime_mod. (jpp, 7/10/09)
            JLOOP = LREORDER(JLOOPLO + KLOOP)
!jpt            JLOOP = JLOOPLO + KLOOP
            ! get box info
            I = IXSAVE( JLOOP )
            J = IYSAVE( JLOOP )
            L = IZSAVE( JLOOP )

            ! store the actual rate constant from
            ! calcrate.f using RRATE common variable
            RATE  = RRATE(KLOOP,PREAC(R))

            ! get the local time:
            loctime = get_localtime(J)

            ! now only select reaction rates between
            ! 09:00 and 17:00 (btwn 8am and 5pm local time)
!            IF ( LT_Br(I,J) > 0 ) then
            if ( (loctime >= 9.0) .and. 
     &           (loctime <= 17.0) ) then
!
!               if ( its_in_the_trop(i,j,l) ) then

                  ! Avoid overflow/underflow
                  IF ( RATE < SMALLEST ) RATE = 0e0
                  IF ( RATE > LARGEST  ) RATE = LARGEST

                  ! Double check on LREORDER... what exactly
                  ! does it do. Should I use it here or not.
                  PRRATE(I,J,L,R) = RATE

                  AD53(I,J,L,R) = AD53(I,J,L,R) 
     &                 + PRRATE(I,J,L,R)

                  ! count how many times you save a conc
                  rate_const_count(I,J,L,R) =  1 + 
     &                 rate_const_count(I,J,L,R)

!               endif

!               ENDIF            ! troposphere selection
!
            endif               ! Time selection

         ENDDO
      ENDDO



      ! return to calling program
      END SUBROUTINE ARCHIVE_RXNS_NEW

!------------------------------------------------------------------------------

!----------------------------------------------------------------
      subroutine save_rxn_rates

! this is going to take the place of calculating lifetimes here
! will just save out rate constants
      USE COMODE_MOD,     ONLY : IXSAVE, IYSAVE, IZSAVE
      USE COMODE_MOD,     ONLY : CSPEC
      use error_mod,      only : geos_chem_stop
      USE DIAG_MOD,       ONLY : AD53
      USE TROPOPAUSE_MOD, ONLY : ITS_IN_THE_TROP
      USE DAO_MOD,        ONLY : SUNCOS
      USE TIME_MOD,       only : GET_LOCALTIME

#     include "CMN_SIZE"
#     include "comode.h"

      ! internal variables
      INTEGER :: KLOOP, JLOOP, ix, iy, iz !, I, J, L
      INTEGER :: IJWINDOW, R
      real*8  :: loctime

      ! begin here
      DO R = 1, NPREAC
         ! -----------------------------------
         ! changed this loop to ensure proper
         ! indexing for v.8-02-02, using
         ! LREORDER. (jpp, 7/10/09)
         ! -----------------------------------
         do ix = 1, IIPAR
            do iy = 1, JJPAR

               ! get the local time:
!               loctime = get_localtime(iy)

               do iz = 1, LLPAR

                  ! index for suncosine
!                  IJWINDOW  = (iy-1)*IIPAR + ix

                  ! now only select reaction rates between
                  ! 09:00 and 17:00 (btwn 8am and 5pm local time)
!                  if ( (loctime >= 9.0) .and. 
!     &                 (loctime <= 17.0) ) then

!                     if (ITS_IN_THE_TROP(ix,iy,iz)) THEN
                        ! now save the rate constants
                  AD53(ix,iy,iz,R) = AD53(ix,iy,iz,R) 
     &                 + PRRATE(ix,iy,iz,R)

                  ! count how many times you save a conc
                  rate_const_count(ix,iy,iz,R) =  1 + 
     &                 rate_const_count(ix,iy,iz,R)

                        ! jpp, debugging
!jpt                        write(6,*) 'jpp inside archive_rxns_new'
!jpt                        if ( (iy == 20 ) .and. (ix <= 55) .and.
!jpt     &                       (ix > 45) .and. (iz == 10) ) then
!jpt                           write(6,*) ' PRRATE = ', PRRATE(ix,iy,iz,R)
!jpt                        endif
!jpt                        print*, 'R =', R
!jpt                        call flush(6)
!jpt                        call geos_chem_stop

!                     endif
!                  endif


!jpt                  ! Treat daylight-only species
!jpt                  ! as all but the VSLs and CH3Br loss rates.
!jpt                  select case (R)
!jpt                  CASE(1,2,  ! Br2
!jpt     &                    3,4,5,6,7,8, ! Br
!jpt     &                    9,10,11,12,13,14,15, ! BrO
!jpt     &                    16,   ! HOBr
!jpt     &                    17,   ! HBr
!jpt     &                    18,   ! BrNO2
!jpt     &                    19, 20, 21, 22) ! BrNO3
!jpt                     ! make sure it's in the trop and that
!jpt                     ! it's in sunlight
!jpt                     IF ( (SUNCOS(IJWINDOW) > 0) .AND. 
!jpt     &                 (ITS_IN_THE_TROP(ix,iy,iz)) ) THEN
!jpt
!jpt                        ! now save the rate constants
!jpt                        AD53(ix,iy,iz,R) = AD53(ix,iy,iz,R) 
!jpt     &                       + PRRATE(ix,iy,iz,R)
!jpt
!jpt                        ! count how many times you save a conc
!jpt                        rate_const_count(ix,iy,iz,R) =  1 + 
!jpt     &                       rate_const_count(ix,iy,iz,R)
!jpt                     endif
!jpt                  CASE(23:26) ! VSLs and CHBr3
!jpt                     IF (ITS_IN_THE_TROP(ix,iy,iz)) THEN
!jpt                        ! now save the rate constants
!jpt                        AD53(ix,iy,iz,R) = AD53(ix,iy,iz,R) 
!jpt     &                       + PRRATE(ix,iy,iz,R)
!jpt
!jpt                        ! count how many times you save a conc
!jpt                        rate_const_count(ix,iy,iz,R) =  1 + 
!jpt     &                       rate_const_count(ix,iy,iz,R)
!jpt                     endif
!jpt                  end select

               enddo
            enddo
         enddo
      enddo                     ! loop over all reactions

      ! return to calling program
      return

      END SUBROUTINE save_rxn_rates
!----------------------------------------------------------------



!------------------------------------------------------------------------------

      !-------------------------------
      ! jpp, debugging
      SUBROUTINE TEST_LT(flag)

      !-------------------------------
      USE COMODE_MOD,     ONLY : IXSAVE, IYSAVE, IZSAVE
      USE COMODE_MOD,     ONLY : CSPEC
      use error_mod,      only : geos_chem_stop
      USE DIAG_MOD,       ONLY : AD53
!      USE DIRECTORY_MOD,  ONLY : RUN_DIR
      USE GRID_MOD,       ONLY : get_xmid, get_ymid
      USE DAO_MOD,        ONLY : suncos

      ! include statements
#     include "CMN_SIZE"
#     include "comode.h"

      ! arguments
      logical, intent(IN), optional :: flag
      character(len=23) :: log_dir = '/home/jpp/testrun/logs/'
      character(len=20) :: names(ntra_4lt)

      ! internal variables
      INTEGER :: I, J, K, N, JLOOP, IJWINDOW
      real*8  :: lat, lon
      ! control integer for selecting longitude
      integer,parameter :: lonsel = 10

      ! set the name data for header
      data names/'Br2','Br','BrO','HOBr','HBr','BrNO2',
     &     'BrNO3', 'CHBr3', 'CH2Br2', 'CH3Br'/

      ! begin here

      open(unit=123,file=log_dir//'test_ltmod.dat', 
     &     status='unknown', action='write')

      ! now calculations recording...
      lon = get_xmid( lonsel )
      write(123,*) 'longitude =', lon
      write(123,113) 'lat ', 'suncos',(trim(adjustl(names(N))), 
     &     N = 1, NTRA_4LT)

      ! begin writing here
      DO JLOOP = 1, NTLOOP
         I = IXSAVE(JLOOP)
         J = IYSAVE(JLOOP)
         K = IZSAVE(JLOOP)
         ! record the locations in lat x lon
         lat = GET_YMID( J )
         ! are we in the sun?
         IJWINDOW  = (J-1)*IIPAR + I
         IF ( I == lonsel .and. K == 1) then
            write(123,114) lat, suncos(ijwindow),
     &           (AD53(I, J, K, N), N = 1, 32)
!            write(123,114) lat, suncos(ijwindow),
!     &           (RXN_RATE(JLOOP, N), N = 1, 25)
!            write(123,114) lat, suncos(ijwindow),
!     &           (PRRATE(JLOOP, N), N = 1, NRXN_4LT)
!            write(123, 112) lat, suncos(ijwindow),
!     &           (MASS(JLOOP, N), N = 1,NTRA_4LT)
!     &           (CSPEC(JLOOP, TRAVAR(N)), N = 1,NTRA_4LT)
!            write(123, 115) lat, suncos(ijwindow),
!     &           O3_lt(JLOOP), OH(JLOOP), HO2(JLOOP), 
!     &           NO(JLOOP), NO2_lt(JLOOP),
!     &           NO3_lt(JLOOP), CH2O_lt(JLOOP)

         endif
      enddo
 112  format(f6.2, 1x, f5.3, 10(1x, es14.6))
 113  format(a3,11(6x,a6))
! 114  format(f6.2, 1x, f5.3, 24(1x, es14.6))
 114  format(f6.2, 1x, f5.3, 25(1x, es11.3))
 115  format(f6.2, 1x, f5.3, 7(1x, es14.6))

      call flush(123)
      close( 123 )

      ! default is that the flag is false and we don't
      ! quit the run
      if ( present(flag) ) then
!         IF ( (MAXVAL(PRRATE(:,2)) > 0.d0) .AND. flag ) THEN
!            print*, 'stoppping in TEST_LT'
         call flush(6)
         call geos_chem_stop
      !   endif
      else
         return
      endif


      RETURN
      END SUBROUTINE TEST_LT
      !-------------------------------
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------

      SUBROUTINE INIT_LIFETIME
!
!******************************************************************************
!  Subroutine INIT_LIFETIM... initialize lifetime arrays
!  (jpp, 7/9/08)
!
!  NOTES:
!  
!******************************************************************************
!
      ! References to F90 modules
      USE ERROR_MOD, ONLY : ALLOC_ERR, GEOS_CHEM_STOP
      USE FILE_MOD,  ONLY : IU_FILE,   IOERROR
      ! initializing AD53 here, since it's hardwired
      USE DIAG_MOD,  ONLY : AD53
      

#     include "CMN_SIZE" ! Size Parameters
#     include "comode.h" ! NTLOOP, KTLOOP, IIPAR, etc.

      ! Local variables 
      LOGICAL            :: IS_INIT = .FALSE.
      INTEGER            :: N, AS, IOS
      CHARACTER(LEN=20)  :: LINE

      !=================================================================
      ! INIT_LIFETIME begins here!
      !=================================================================

      !-------------------------
      ! Arrays of size NPREAC
      !-------------------------      
      ALLOCATE( PREAC( NRXN_4LT ), STAT=AS )
      IF ( AS /= 0 ) THEN
         print *, 'AS = ', AS
         CALL FLUSH(6)
         CALL ALLOC_ERR( 'PREAC' )
      ENDIF


      ALLOCATE( PRRATE( IIPAR, JJPAR, LLPAR,
     &     NRXN_4LT ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'PRRATE' )

      !--------------------------
      ! Arrays of size MAXVARS
      ! jpp: eventually, taylor this to
      !      NRXN_4LT and NTRA_4LT to get dimensions
      !--------------------------
      ALLOCATE( RXNVAR( NRXN_4LT ), STAT=AS )
      IF ( AS /= 0 ) THEN
         print *, 'AS = ', AS
         CALL FLUSH(6)
         CALL ALLOC_ERR( 'RXNVAR' )
      ENDIF

      ALLOCATE( RXN_RATE( NTLOOP, NRXN_4LT ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'RXN_RATE' )

      ALLOCATE( TRAVAR( NTRA_4LT ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'TRAVAR' )

      ! jpp: storing the mass in each grid
      ALLOCATE( MASS( NTLOOP, NTRA_4LT ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'MASS' )

      ! concentration array of relevent species
      ALLOCATE( CONCS(NTLOOP, NTRA_4LT), STAT=AS)
      IF (AS /= 0) CALL ALLOC_ERR('CONCS')

      ! final array containing spatially resolved
      ! liftimes for each reaction and each bromine
      ! species... has 5 dimensions as a result.
      ALLOCATE( TAU_LT( IIPAR, JJPAR, LLPAR, NRXN_4LT,
     &     NTRA_4LT), STAT=AS )
      IF (AS /= 0 ) CALL ALLOC_ERR( 'TAU_LT' )

      !------------------------------
      ! jpp: counting array used in
      !      taking the average over
      !      all of the lifetimes
      !------------------------------
      ALLOCATE( RXN_COUNT(IIPAR, JJPAR, LLPAR), STAT=AS)
      IF (AS /= 0) CALL ALLOC_ERR('RXN_COUNT')

      ALLOCATE( SPC_COUNT(IIPAR, JJPAR, LLPAR), STAT=AS)
      IF (AS /= 0) CALL ALLOC_ERR('SPC_COUNT')

      ALLOCATE( rate_const_count(IIPAR, JJPAR, LLPAR,
     &     NRXN_4LT), STAT=AS)
      IF (AS /= 0) CALL ALLOC_ERR('rate_const_count')

      !----------------------------------------------
      ! allocating the AD53 diagnostic array.
      ! eventually move to ndxx_setup.f using
      ! a new lifetime diagnostic menu in input.geos
      !----------------------------------------------
      ALLOCATE( AD53( IIPAR, JJPAR, LLPAR, 
     &     nrxn_4lt ), STAT=AS )
      IF (AS /= 0) CALL ALLOC_ERR('lifetime_mod: AD53')

      ! ------------------------------------
      ! new arrays for reactant concs
      ! ------------------------------------
      ! tracers used here for reaction rates
      ALLOCATE( O3_lt(NTLOOP), STAT=AS )
      IF (AS /= 0) CALL ALLOC_ERR('O3_lt')

      ALLOCATE( HO2_lt(NTLOOP), STAT=AS )
      IF (AS /= 0) CALL ALLOC_ERR('HO2_lt')

      ALLOCATE( NO3_lt(NTLOOP), STAT=AS )
      IF (AS /= 0) CALL ALLOC_ERR('NO3_lt')

      ALLOCATE( CH2O_lt(NTLOOP), STAT=AS )
      IF (AS /= 0) CALL ALLOC_ERR('CH2O_lt')

      ALLOCATE( NO2_lt(NTLOOP), STAT=AS )
      IF (AS /= 0) CALL ALLOC_ERR('NO2_lt')

      ALLOCATE( NO_lt(NTLOOP), STAT=AS )
      IF (AS /= 0) CALL ALLOC_ERR('NO_lt')

      ALLOCATE( OH_lt(NTLOOP), STAT=AS )
      IF (AS /= 0) CALL ALLOC_ERR('OH_lt')

      ! note, this is a hack... 7 is for the number
      ! of co-reactant species that are contributing
      ! to the loss reactions of interest.
      ALLOCATE( SP2_COUNT( NTLOOP, 7 ), STAT=AS )
      IF (AS /= 0) CALL ALLOC_ERR('SP2_COUNT')

      !=================================================================
      ! Initialize arrays 
      !=================================================================
      ! initializing tau_lt to -999 to avoid mis-matching it
      ! with the legitimate 0's that will come about from
      ! zero concentrations...
      TAU_LT(:,:,:,:,:) = -999.d0
      ! initializing AD53 here for now, but once
      ! its allocation is moved to ndxx_setup.f, shift
      ! the initialization to initialize.f. jpp 7/09/08
      AD53(:,:,:,:) = 0.d0
      RXN_COUNT(:,:,:) = 0
      rate_const_count(:,:,:,:) = 0 
      RXN_RATE(:,:) = 0.d0
      SPC_COUNT(:,:,:)=0
      SP2_COUNT(:, :) =0
      CONCS(:,:) = 0.d0
      MASS(:,:) = 0.d0
      PREAC  = 0
      PRRATE(:,:,:,:) = 0d0
      RXNVAR  = 0
      TRAVAR  = 0
      O3_lt = 0d0
      HO2_lt = 0d0
      NO3_lt = 0d0
      CH2O_lt = 0d0
      NO2_lt = 0d0
      NO_lt  = 0d0
      OH_lt = 0d0

      ! Return to calling program
      END SUBROUTINE INIT_LIFETIME

!------------------------------------------------------------------------------

      SUBROUTINE CLEANUP_LIFETIME
!
!******************************************************************************
!  Subroutine CLEANUP_PLANEFLIGHT deallocates all allocatable module arrays.
!  (mje, bmy, 7/1/02, 4/1/03)
!
!  NOTES:
!  (1 ) Renamed PRATE to PRRATE to avoid conflict w/ SMVGEAR II (bmy, 4/1/03)
!******************************************************************************
!

      ! use the following modules

      print*, 'lt 1'
      IF ( ALLOCATED( RXNVAR    ) ) DEALLOCATE( RXNVAR    )
      IF ( ALLOCATED( RXN_RATE  ) ) DEALLOCATE( RXN_RATE  )
      IF ( ALLOCATED( TRAVAR    ) ) DEALLOCATE( TRAVAR    )
      IF ( ALLOCATED( PREAC     ) ) DEALLOCATE( PREAC     )
      print*, 'lt 2'
      IF ( ALLOCATED( PRRATE    ) ) DEALLOCATE( PRRATE    )
      IF ( ALLOCATED( CONCS     ) ) DEALLOCATE( CONCS     )
      IF ( ALLOCATED( MASS      ) ) DEALLOCATE( MASS      )
      IF ( ALLOCATED( RXN_COUNT ) ) DEALLOCATE( RXN_COUNT )
      IF ( ALLOCATED( TAU_LT    ) ) DEALLOCATE( TAU_LT    )
      print*, 'lt 3'
      IF ( ALLOCATED( O3_lt     ) ) DEALLOCATE( O3_lt     )
      IF ( ALLOCATED( HO2_lt    ) ) DEALLOCATE( HO2_lt    )
      IF ( ALLOCATED( NO3_lt    ) ) DEALLOCATE( NO3_lt    )
      IF ( ALLOCATED( CH2O_lt   ) ) DEALLOCATE( CH2O_lt   )
      IF ( ALLOCATED( NO2_lt    ) ) DEALLOCATE( NO2_lt    )
      print*, 'lt 4'
      IF ( ALLOCATED( NO_lt     ) ) DEALLOCATE( NO_lt     )
      print*, 'lt5'
      IF ( ALLOCATED( OH_lt     ) ) DEALLOCATE( OH_lt     )
      print*, 'lt6'
      call flush(6)
      IF ( ALLOCATED( SP2_COUNT ) ) then
         print *, 'sp2_count is allocated'
         call flush(6)
!         DEALLOCATE( SP2_COUNT )
      endif
      print*, 'lt7'
      call flush(6)
      ! Return to calling program
      END SUBROUTINE CLEANUP_LIFETIME

!------------------------------------------------------------------------------

!------------------------------------------------------------------------------

!      SUBROUTINE WRITE_4LT( IND, VARI )
!
!******************************************************************************
!  Subroutine WRITE_VARS_TO_FILE writes the values of all the variables for
!  a given flight track point to the output file. (mje, bmy, 7/8/02. 3/25/05)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) IND  (INTEGER) : Number of the flight track point to print to file 
!  (2 ) VARI (REAL*4 ) : Array holding variable values to print to file
!
!  NOTES:
!  (1 ) The max line length for output seems to be 1024 characters.  Adjust
!        MAXVARS accordingly so that we don't exceed this. (bmy, 7/8/02)
!  (2 ) Now do not write file header -- this is now done in subroutine
!        SETUP_PLANEFLIGHT at the start of each day (bmy, 3/25/05)
!******************************************************************************
!
      ! References to F90 modules
!      USE FILE_MOD, ONLY : IOERROR
!      USE COMODE_MOD, ONLY : IXSAVE, IYSAVE, IZSAVE
!      USE BPCH2_MOD ! use all the bpch output tools
!
!#     include comode.h
!
!      ! Arguments
!      INTEGER, INTENT(IN) :: IND, tmpnl
!      REAL*8,  INTENT(IN) :: VARI(NPVAR)
!     
!      ! Local variables
!      LOGICAL, SAVE       :: FIRST = .TRUE.
!      INTEGER             :: I, IOS
!
!      CHARACTER(LEN=100)  :: title
!      ! this I/O unit is currently hardwired. eventually
!      ! build it into the I/O unit module.
!      INTEGER             :: IU_LIFETIME=132
!
!      ! parameters for bpch writing
!      INTEGER, PARAMETER :: CENTER180 = 1
!
!      ! loop control
!      INTEGER :: NN, BOXLOOP, jrxn, I, J, L
!
!      !=================================================================
!      ! WRITE_VARS_TO_FILE begins here!
!      !=================================================================
!
!      ! open bpch file for writing
!      title = 'GEOS-Chem binary punch file v 2.0: specific'//
!     &     ' for the lifetimes module output'
!      title = TRIM( ADJUSTL( title ) )
!      OPEN_BPCH2_FOR_WRITE( IU_LIFETIME, 'gc.lifetimes.bpch',title)
!
!      ! get model and output info from BPCH2_MOD to aid in
!      ! printing to the binary punch file
!      MODELNAME = GET_MODELNAME()
!      HALFPOLAR = GET_HALFPOLAR()
!      LONRES    = DISIZE
!      LATRES    = DJSIZE
!
!      ! write to bpch file
!
!      
!!      CALL BPCH2( IU_LIFETIME, MODELNAME, LONRES, LATRES,
!!     &     HALFPOLAR, CENTER180, IIPAR, JJPAR, NRXN_4LT, NTRA_4LT,
!!     &     
!
!      ! Write data to file
!      DO NN = 1, NTRA_4LT ! loop over the number of tracers
!      DO BOXLOOP = 1, NTLOOP ! loop over all the boxes
!         ! store the number of loss reactions
!         ! to keep track of
!         tmpnl = NLOSS(NN)
!
!         ! retrieve 3D counting variables
!         ! for storing the lifetimes
!         I = IXSAVE(BOXLOOP)
!         J = IYSAVE(BOXLOOP)
!         L = IZSAVE(BOXLOOP)
!
!         WRITE( IU_LIFETIME, 110, IOSTAT=IOS )
!         &     TAU_LT(I,J,L,jrxn,NN), jrxn=1,tmpnl )
!
!      enddo
!      enddo
!
!      ! Format string
! 110  FORMAT(I5,X,A5,X,I8.8,X,I4.4,X,F7.2,X,F7.2,X,F7.2,X,95(es10.3,x))
!
!      ! Error check
!      IF ( IOS /= 0 ) CALL IOERROR( IOS,IU_PLANE,'write_vars_to_file:1')
!
!      ! Flush the file to disk
!      CALL FLUSH( IU_LIFETIME )
!
!      ! Close Device
!      CLOSE(IU_LIFETIME)

      ! Return to calling program
!      END SUBROUTINE WRITE_4LT

!------------------------------------------------------------------------------

      END MODULE LIFETIME_MOD
