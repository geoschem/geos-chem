!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: passive_tracer_mod.F90
!
! !DESCRIPTION: Module passive\_tracer\_mod.F90 contains variables and routines
!  for using passive tracers in GEOS-Chem. Passive tracers are tracers that are 
!  passively transported by GEOS-Chem, with a simple first order loss rate 
!  being applied to each tracer. The number of passive tracers, corresponding 
!  tracer properties as well as loss rates and default initial concentrations 
!  can be specified by the user via the GEOS-Chem input file (input.geos).
! 
!  The passive tracer module is designed to work in combination with any
!  existing GEOS-Chem simulation type, even though it has only been tested
!  with the Radon simulation at this point.
!
! !REMARKS:  
!  Passive tracers are defined in input.geos in the PASSIVE TRACERS menu. For
!  instance, to use the Radon simulation with two passive tracers ('Rn_ps' and 
!  'Dummy') with atmospheric lifetimes of 3.8 days and 1 hour, respectively, 
!  add the following entries to input.geos:
!
!  %%% PASSIVE TRACERS %%% :
!  Number of pass. tracers : 2
!  Passive tracer #1       : Rn_ps 328320.0 1.0e-20
!  Passive tracer #2       : Dummy 3600.0   1.0e-20
!
!  The 3rd column of the tracer definition denotes the default initial
!  concentration (in v/v) of the species of interest (1.0e-20 v/v in this
!  case). This value will be used if the GEOS-Chem tracer restart file has
!  no concentration field for the given tracer. 
!
!  There must be a matching entry in the tracers menu for every passive tracer
!  defined in the passive tracers menu:
!
!  %%% TRACER MENU %%%     :
!  Type of simulation      : 1
!  Number of Tracers       : 4
!  Tracer Entries -------> : TR#   Name  g/mole   Tracer Members; () = emitted
!  Tracer #1               :   1   Rn     222.0
!  Tracer #2               :   2   Pb     210.0
!  Tracer #3               :   3   Be7      7.0
!  Tracer #4               :   4   Rn_ps  222.0
!  Tracer #5               :   5   Dummy  100.0
! 
!  In this example, tracers 1-3 are the default tracers for this simulation type
!  while tracers 4-5 are the user-specific passive tracers.
!
!  As for regular GEOS-Chem tracers, emissions can be assigned to passive 
!  tracers via the HEMCO configuration file. For example, to assign a uniform 
!  flux of 0.1 kg/m2/s to passive tracer 'Dummy', add the following line to 
!  section base emissions of your HEMCO_Config.rc:
!
!  0 DUMMY_EMIS 0.1 - - - xy kg/m2/s Dummy - 1 1
!\\
!\\
! !INTERFACE:
!
MODULE PASSIVE_TRACER_MOD 
!
! !USES:
!
  USE PRECISION_MOD

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC :: INIT_PASSIVE_TRACER
  PUBLIC :: ADD_PASSIVE_TRACER
  PUBLIC :: PASSIVE_TRACER_GETRATE
  PUBLIC :: PASSIVE_TRACER_INQUIRE
  PUBLIC :: CLEANUP_PASSIVE_TRACER
!
! !PRIVATE MEMBER FUNCTIONS:
!
! !REVISION HISTORY:
!  04 Sep 2015 - C. Keller   - Initial version. 
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !MODULE VARIABLES:
!
  INTEGER,                        PUBLIC :: NPASSIVE = 0
  CHARACTER(LEN=63), ALLOCATABLE, PUBLIC :: PASSIVE_NAME    (:)
  INTEGER,  ALLOCATABLE                  :: PASSIVE_ID      (:)
  REAL(fp), ALLOCATABLE                  :: PASSIVE_MW      (:)
  REAL(fp), ALLOCATABLE                  :: PASSIVE_TAU     (:)
  REAL(fp), ALLOCATABLE                  :: PASSIVE_INITCONC(:)

  !============================================================================
  ! MODULE ROUTINES -- follow below the "CONTAINS" statement
  !============================================================================
CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Init_Passive_Tracer
!
! !DESCRIPTION: Subroutine INIT\_PASSIVE\_TRACER initializes the passive
! tracers arrays.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE INIT_PASSIVE_TRACER ( am_I_Root, NPT, RC ) 
!
! !USES:
!
    USE ErrCode_Mod
!
! !INPUT PARAMETERS:
!
    LOGICAL,          INTENT(IN   )  :: am_I_Root  ! root CPU?
    INTEGER,          INTENT(IN   )  :: NPT        ! # of passive tracers 
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(INOUT)  :: RC         ! Failure or success
!
! !REMARKS:
!
! !REVISION HISTORY: 
!  04 Sep 2015 - C. Keller    - Initial version 
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER         :: AS

    !=================================================================
    ! INIT_PASSIVE_TRACER begins here!
    !=================================================================

    ! Assume success
    RC = GC_SUCCESS
  
    ! Set module variable 
    NPASSIVE = NPT
    
    ! Allocate arrays
    IF ( NPASSIVE > 0 ) THEN
       ALLOCATE( PASSIVE_ID(NPASSIVE),       PASSIVE_TAU(NPASSIVE), &
                 PASSIVE_INITCONC(NPASSIVE), PASSIVE_MW(NPASSIVE),  &
                 PASSIVE_NAME(NPASSIVE),     STAT=AS )
       IF ( AS /= 0 ) THEN
          WRITE(*,*) 'Cannot allocate passive tracers arrays'
          RC = GC_FAILURE
          RETURN
       ENDIF

       PASSIVE_ID(:)       = -999
       PASSIVE_NAME(:)     = ''
       PASSIVE_MW(:)       = 0.0_fp
       PASSIVE_TAU(:)      = 0.0_fp
       PASSIVE_INITCONC(:) = 0.0_fp
    ENDIF

  END SUBROUTINE INIT_PASSIVE_TRACER
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Add_Passive_Tracer
!
! !DESCRIPTION: Subroutine ADD\_PASSIVE\_TRACER registers a passive tracer 
!  based on the passed input arguments.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE ADD_PASSIVE_TRACER ( am_I_Root,   TrcName, TrcTau, &
                                  TrcInitConc, TrcMW,   RC ) 
!
! !USES:
!
    USE ErrCode_Mod
!
! !INPUT PARAMETERS:
!
    LOGICAL,          INTENT(IN   )  :: am_I_Root   ! root CPU?
    CHARACTER(LEN=*), INTENT(IN   )  :: TrcName     ! Tracer name
    REAL(fp),         INTENT(IN   )  :: TrcTau      ! Tracer lifetime (s) 
    REAL(fp),         INTENT(IN   )  :: TrcInitConc ! Tracer default init conc (v/v) 
    REAL(fp),         INTENT(IN   )  :: TrcMW       ! Tracer molec. weight (g/mol) 
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(INOUT)  :: RC          ! Failure or success
!
! !REMARKS:
!
! !REVISION HISTORY: 
!  04 Sep 2015 - C. Keller    - Initial version 
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER             :: I, IDX, TRCID
    CHARACTER(LEN=255)  :: MSG
    CHARACTER(LEN=255)  :: LOC = 'ADD_PASSIVE_TRACER (passive_tracer_mod.F90)'

    !=================================================================
    ! ADD_PASSIVE_TRACER begins here!
    !=================================================================

    ! Assume success
    RC = GC_SUCCESS

    ! Error check: cannot define passive tracer if # of passive tracers is 0.
    IF ( NPASSIVE <= 0 ) THEN
       WRITE(*,*) 'Cannot add passive tracer ', TRIM(TrcName), &
             ': # of passive tracers is smaller than 1!'
       RC = GC_FAILURE
       RETURN 
    ENDIF 

    ! Find empty slot
    IDX = -1
    DO I = 1, NPASSIVE
       IF ( PASSIVE_ID(I) < 0 ) THEN
          IDX = I
          EXIT
       ENDIF
    ENDDO 
    
    IF ( IDX <= 0 ) THEN
       WRITE(*,*) 'Cannot add passive tracer ', TRIM(TrcName), &
                    ': all ', NPASSIVE, ' tracer slots are already being used.'
       RC = GC_FAILURE
       RETURN 
    ENDIF 

!    ! Find GEOS-Chem tracer ID for this tracer (by name)
!    TRCID = ind_( TRIM(TrcName) )

!    ! Return w/ error if this tracer is not defined as GEOS-Chem tracer
!    IF ( TRCID <= 0 ) THEN
!       WRITE(MSG,*) 'Cannot add passive tracer ', TRIM(TrcName), &
!                    ': this is not a GEOS-Chem tracer.'
!       CALL ERROR_STOP ( TRIM(MSG), TRIM(LOC) )
!       RC = GC_FAILURE
!       RETURN 
!    ENDIF 

    ! Register tracer
    PASSIVE_ID(IDX)       = IDX
    PASSIVE_NAME(IDX)     = TRIM(TrcName)
    PASSIVE_TAU(IDX)      = TrcTau 
    PASSIVE_INITCONC(IDX) = TrcInitConc 
    PASSIVE_MW(IDX)       = TrcMW

    ! Verbose
    IF ( am_I_Root ) THEN
       WRITE( 6, '(a)' ) 'Added passive tracer: '
       WRITE( 6, 110   ) ' - Tracer name                 : ', PASSIVE_NAME(IDX) 
       WRITE( 6, 120   ) ' - Molec. weight [g/mol]       : ', PASSIVE_MW(IDX)
       WRITE( 6, 120   ) ' - Lifetime [s]                : ', PASSIVE_TAU(IDX)
       WRITE( 6, 130   ) ' - Default concentration [v/v] : ', PASSIVE_INITCONC(IDX)
    ENDIF

110 FORMAT( A, A )
120 FORMAT( A, F10.2  )
130 FORMAT( A, ES10.2 )
 
    ! Return w/ success
    RC = GC_SUCCESS

  END SUBROUTINE ADD_PASSIVE_TRACER
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Passive_Tracer_Getrate
!
! !DESCRIPTION: Subroutine PASSIVE\_TRACER\_GETRATE returns the unitless decay
!  rate for the given tracer and chemistry time step. 
!  on all passive tracers.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE PASSIVE_TRACER_GETRATE ( am_I_Root, TrcName, DT, Rate, RC ) 
!
! !USES:
!
    USE ErrCode_Mod
!
! !INPUT PARAMETERS:
!
    LOGICAL,          INTENT(IN   )  :: am_I_Root  ! root CPU?
    CHARACTER(LEN=*), INTENT(IN   )  :: TrcName    ! Passive tracer name 
    REAL(fp),         INTENT(IN   )  :: DT         ! Time step in s
!
! !OUTPUT PARAMETERS:
!
    REAL(fp),         INTENT(  OUT)  :: Rate       ! Decay rate (unitless 
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(INOUT)  :: RC         ! Failure or success
!
! !REMARKS:
!
! !REVISION HISTORY: 
!  04 Sep 2015 - C. Keller    - Initial version 
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER             :: N, ID
    REAL(fp)            :: Decay
    CHARACTER(LEN=255)  :: MSG
    CHARACTER(LEN=255)  :: LOC = "PASSIVE_TRACER_GETRATE (PASSIVE_TRACER_MOD.F90)"

    REAL(fp), PARAMETER :: ln2 = 0.693147181E+00_fp

    !=================================================================
    ! PASSIVE_TRACER_GETRATE begins here!
    !=================================================================

    ! Assume success
    RC = GC_SUCCESS

    ! Find the index
    ID = -1
    DO N = 1, NPASSIVE

       IF ( TRIM(TrcName) == TRIM(PASSIVE_NAME(N)) ) THEN
          ID = N
          EXIT
       ENDIF
    ENDDO

    ! Error check
    IF ( ID <= 0 ) THEN
       WRITE(*,*) 'This is not a passive tracer ', TRIM(TrcName)
       RC = GC_FAILURE
       RETURN 
    ENDIF 

    ! No loss needed if tau is zero or negative
    IF ( PASSIVE_TAU(N) <= 0.0_fp ) THEN
       Rate = 1.0

    ! Calculate decay rate (unitless)
    ELSE 
       Decay    = ln2 / PASSIVE_TAU(N)
       Rate     = EXP( - DT * Decay )

    ENDIF

    ! Return w/ success 
    RC = GC_SUCCESS
 
  END SUBROUTINE PASSIVE_TRACER_GETRATE
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Passive_Tracer_Inquire
!
! !DESCRIPTION: Function PASSIVE\_TRACER\_INQUIRE is a wrapper routine to 
!  inquire information about a passive tracer. 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE PASSIVE_TRACER_INQUIRE ( TrcName, IsPassive, MW, InitConc ) 
!
! !USES:
!
    USE ErrCode_Mod
!
! !INPUT PARAMETERS:
!
    CHARACTER(LEN=*),   INTENT(IN   )  :: TrcName    ! GC tracer name 
!
! !OUTPUT PARAMETERS:
!
    LOGICAL,  OPTIONAL, INTENT(  OUT)  :: IsPassive  ! Is TrcID a passive tracer?
    REAL(fp), OPTIONAL, INTENT(  OUT)  :: MW         ! Molecular weight (g/mol) 
    REAL(fp), OPTIONAL, INTENT(  OUT)  :: InitConc   ! Initial concentration (v/v) 
!
! !REMARKS:
!
! !REVISION HISTORY: 
!  04 Sep 2015 - C. Keller    - Initial version 
!  09 Mar 2017 - C. Keller    - Bug fix: initialize molw to default value of 0.0
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER    :: N, ID
    LOGICAL    :: IsPass
    REAL(fp)   :: InConc, molw

    !=================================================================
    ! PASSIVE_TRACER_INQUIRE begins here!
    !=================================================================

    ! Init
    IsPass = .FALSE.
    InConc = 0.0_fp 
    molw   = 0.0_fp 

    ! Nothing to do if no passive tracers defined
    IF ( NPASSIVE > 0 ) THEN 

       ! Loop over all passive tracers
       DO N = 1, NPASSIVE
          IF ( TRIM(TrcName) == TRIM(PASSIVE_NAME(N)) ) THEN
             IsPass = .TRUE.
             InConc = PASSIVE_INITCONC(N)
             molw   = PASSIVE_MW(N)
             EXIT
          ENDIF
       ENDDO
    ENDIF

    ! Pass to output
    IF ( PRESENT(IsPassive) ) IsPassive = IsPass
    IF ( PRESENT(InitConc ) ) InitConc  = InConc
    IF ( PRESENT(MW       ) ) MW        = molw    

  END SUBROUTINE PASSIVE_TRACER_INQUIRE
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Cleanup_Passive_Tracer
!
! !DESCRIPTION: Subroutine CLEANUP\_PASSIVE\_TRACER finalizes the passive
! tracers arrays.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE CLEANUP_PASSIVE_TRACER ( am_I_Root, RC ) 
!
! !USES:
!
    USE ErrCode_Mod
!
! !INPUT PARAMETERS:
!
    LOGICAL,          INTENT(IN   )  :: am_I_Root  ! root CPU?
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(INOUT)  :: RC         ! Failure or success
!
! !REMARKS:
!
! !REVISION HISTORY: 
!  04 Sep 2015 - C. Keller    - Initial version 
!EOP
!------------------------------------------------------------------------------
!BOC

    !=================================================================
    ! CLEANUP_PASSIVE_TRACER begins here!
    !=================================================================

    ! Deallocate arrays
    IF ( ALLOCATED( PASSIVE_ID       ) ) DEALLOCATE( PASSIVE_ID       )
    IF ( ALLOCATED( PASSIVE_TAU      ) ) DEALLOCATE( PASSIVE_TAU      )
    IF ( ALLOCATED( PASSIVE_INITCONC ) ) DEALLOCATE( PASSIVE_INITCONC )
    IF ( ALLOCATED( PASSIVE_MW       ) ) DEALLOCATE( PASSIVE_MW       )
    IF ( ALLOCATED( PASSIVE_NAME     ) ) DEALLOCATE( PASSIVE_NAME     )

    ! Return w/ success
    RC = GC_SUCCESS

  END SUBROUTINE CLEANUP_PASSIVE_TRACER
!EOC
END MODULE PASSIVE_TRACER_MOD 
