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
! !EXAMPLE:  
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
  PUBLIC :: CHEM_PASSIVE_TRACER
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
! !DEFINED PARAMETERS:
!
!
! !MODULE VARIABLES:
!
  INTEGER                :: NPASSIVE = 0
  INTEGER,  ALLOCATABLE  :: PASSIVE_ID      (:)
  REAL(fp), ALLOCATABLE  :: PASSIVE_TAU     (:)
  REAL(fp), ALLOCATABLE  :: PASSIVE_INITCONC(:)

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
! !IROUTINE: INIT_PASSIVE_TRACER
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
    USE GIGC_ErrCode_Mod
    USE ERROR_MOD,          ONLY : ALLOC_ERR 
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
! !REMARKS
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
    RC = GIGC_SUCCESS
  
    ! Set module variable 
    NPASSIVE = NPT
    
    ! Allocate arrays
    IF ( NPASSIVE > 0 ) THEN
       ALLOCATE( PASSIVE_ID(NPASSIVE),       PASSIVE_TAU(NPASSIVE), &
                 PASSIVE_INITCONC(NPASSIVE), STAT=AS )
       IF ( AS /= 0 ) THEN
          CALL ALLOC_ERR( 'passive tracer arrays' )
          RC = GIGC_FAILURE
          RETURN
       ENDIF

       PASSIVE_ID(:)       = -999
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
! !IROUTINE: ADD_PASSIVE_TRACER
!
! !DESCRIPTION: Subroutine ADD\_PASSIVE\_TRACER registers a passive tracer 
!  based on the passed input arguments.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE ADD_PASSIVE_TRACER ( am_I_Root, Input_Opt,   TrcName, &
                                  TrcTau,    TrcInitConc, RC        ) 
!
! !USES:
!
    USE GIGC_ErrCode_Mod
    USE GIGC_Input_Opt_Mod, ONLY : OptInput
    USE ERROR_MOD,          ONLY : ERROR_STOP
    USE GIGC_State_Chm_Mod, ONLY : Get_Indx
!
! !INPUT PARAMETERS:
!
    LOGICAL,          INTENT(IN   )  :: am_I_Root   ! root CPU?
    TYPE(OptInput),   INTENT(IN   )  :: Input_Opt   ! Input opts
    CHARACTER(LEN=*), INTENT(IN   )  :: TrcName     ! Tracer name
    REAL(fp),         INTENT(IN   )  :: TrcTau      ! Tracer lifetime (s) 
    REAL(fp),         INTENT(IN   )  :: TrcInitConc ! Tracer default init conc (v/v) 
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(INOUT)  :: RC          ! Failure or success
!
! !REMARKS
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
    RC = GIGC_SUCCESS

    ! Error check: cannot define passive tracer if # of passive tracers is 0.
    IF ( NPASSIVE <= 0 ) THEN
       MSG = 'Cannot add passive tracer ' // TRIM(TrcName) // &
             ': # of passive tracers is smaller than 1!'
       CALL ERROR_STOP ( TRIM(MSG), TRIM(LOC) )
       RC = GIGC_FAILURE
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
       WRITE(MSG,*) 'Cannot add passive tracer ', TRIM(TrcName), &
                    ': all ', NPASSIVE, ' tracer slots are already being used.'
       CALL ERROR_STOP ( TRIM(MSG), TRIM(LOC) )
       RC = GIGC_FAILURE
       RETURN 
    ENDIF 

    ! Find GEOS-Chem tracer ID for this tracer (by name)
    TRCID = Get_Indx( TRIM(TrcName), Input_Opt%ID_TRACER, Input_Opt%TRACER_NAME ) 

    ! Return w/ error if this tracer is not defined as GEOS-Chem tracer
    IF ( TRCID <= 0 ) THEN
       WRITE(MSG,*) 'Cannot add passive tracer ', TRIM(TrcName), &
                    ': this is not a GEOS-Chem tracer.'
       CALL ERROR_STOP ( TRIM(MSG), TRIM(LOC) )
       RC = GIGC_FAILURE
       RETURN 
    ENDIF 

    ! Register tracer
    PASSIVE_ID(IDX)       = TRCID
    PASSIVE_TAU(IDX)      = TrcTau 
    PASSIVE_INITCONC(IDX) = TrcInitConc 

    ! Verbose
    IF ( am_I_Root ) THEN
       WRITE( 6, '(a)' ) 'Added passive tracer: '
       WRITE( 6, 110   ) ' - Tracer name                 : ', TRIM(TrcName)
       WRITE( 6, 120   ) ' - Lifetime [s]                : ', PASSIVE_TAU(IDX)
       WRITE( 6, 130   ) ' - Default concentration [v/v] : ', PASSIVE_INITCONC(IDX)
    ENDIF

110 FORMAT( A, A )
120 FORMAT( A, F10.2  )
130 FORMAT( A, ES10.2 )
 
    ! Return w/ success
    RC = GIGC_SUCCESS

  END SUBROUTINE ADD_PASSIVE_TRACER
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: CHEM_PASSIVE_TRACER
!
! !DESCRIPTION: Subroutine RUN\_PASSIVE\_TRACER performs loss chemistry 
!  on all passive tracers.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE CHEM_PASSIVE_TRACER ( am_I_Root, Input_Opt, &
                                   State_Met, State_Chm, RC ) 
!
! !USES:
!
    USE GIGC_ErrCode_Mod
    USE GIGC_Input_Opt_Mod,     ONLY : OptInput
    USE GIGC_State_Chm_Mod,     ONLY : ChmState
    USE GIGC_State_Met_Mod,     ONLY : MetState
    USE TIME_MOD,               ONLY : GET_TS_CHEM
!
! !INPUT PARAMETERS:
!
    LOGICAL,          INTENT(IN   )  :: am_I_Root  ! root CPU?
    TYPE(OptInput),   INTENT(IN   )  :: Input_Opt  ! Input options object
    TYPE(MetState),   INTENT(IN   )  :: State_Met  ! Meteorology state object
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ChmState),   INTENT(IN   )  :: State_Chm  ! Chemistry state object 
    INTEGER,          INTENT(INOUT)  :: RC         ! Failure or success
!
! !REMARKS
!
! !REVISION HISTORY: 
!  04 Sep 2015 - C. Keller    - Initial version 
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER           :: I, J, L, N, GCID
    REAL(fp)          :: DTCHEM, Decay, Rate

    REAL(fp), PARAMETER :: ln2 = 0.693147181E+00_fp

    !=================================================================
    ! CHEM_PASSIVE_TRACER begins here!
    !=================================================================

    ! Assume success
    RC = GIGC_SUCCESS

    ! Return if there are no passive tracers
    IF ( NPASSIVE <= 0 ) RETURN

    ! Get chemistry timestep [s]
    DTCHEM = GET_TS_CHEM() * 60.0_fp

    ! Do for every passive tracer
    DO N = 1, NPASSIVE

       ! Sanity check (GCID should never be negative) 
       GCID = PASSIVE_ID(N)
       IF ( GCID <= 0 ) CYCLE

       ! No loss needed if tau is zero or negative
       IF ( PASSIVE_TAU(N) <= 0.0_fp ) CYCLE 

       ! Calculate decay rate (unitless)
       Decay    = ln2 / PASSIVE_TAU(N)
       Rate     = EXP( - DTCHEM * Decay )

       ! Apply loss
       State_Chm%Tracers(:,:,:,GCID) = State_Chm%Tracers(:,:,:,GCID) * Rate

    ENDDO

    ! Return w/ success 
    RC = GIGC_SUCCESS
 
  END SUBROUTINE CHEM_PASSIVE_TRACER
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: PASSIVE_TRACER_INQUIRE
!
! !DESCRIPTION: Function PASSIVE\_TRACER\_INQUIRE is a wrapper routine to 
!  inquire information about a passive tracer. 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE PASSIVE_TRACER_INQUIRE ( am_I_Root, TrcID, RC, IsPassive, InitConc ) 
!
! !USES:
!
    USE GIGC_ErrCode_Mod
!
! !INPUT PARAMETERS:
!
    LOGICAL,            INTENT(IN   )  :: am_I_Root  ! Root CPU? 
    INTEGER,            INTENT(IN   )  :: TrcID      ! GC tracer ID 
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,            INTENT(INOUT)  :: RC         ! Return code 
!
! !OUTPUT PARAMETERS:
!
    LOGICAL,  OPTIONAL, INTENT(  OUT)  :: IsPassive  ! Is TrcID a passive tracer?
    REAL(fp), OPTIONAL, INTENT(  OUT)  :: InitConc   ! Initial concentration (v/v) 
!
! !REMARKS
!
! !REVISION HISTORY: 
!  04 Sep 2015 - C. Keller    - Initial version 
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER    :: N
    LOGICAL    :: IsPass
    REAL(fp)   :: InConc

    !=================================================================
    ! PASSIVE_TRACER_INQUIRE begins here!
    !=================================================================

    ! Init
    IsPass = .FALSE.
    InConc = 0.0_fp 

    ! Nothing to do if no passive tracers defined
    IF ( NPASSIVE > 0 ) THEN 

       ! Loop over all passive tracers
       DO N = 1, NPASSIVE
          IF ( PASSIVE_ID(N) == TrcID ) THEN
             IsPass = .TRUE.
             InConc = PASSIVE_INITCONC(N)
             EXIT
          ENDIF
       ENDDO
    ENDIF

    ! Pass to output
    IF ( PRESENT(IsPassive) ) IsPassive = IsPass
    IF ( PRESENT(InitConc ) ) InitConc  = InConc

    ! Return w/ success
    RC = GIGC_SUCCESS

  END SUBROUTINE PASSIVE_TRACER_INQUIRE
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: CLEANUP_PASSIVE_TRACER
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
    USE GIGC_ErrCode_Mod
!
! !INPUT PARAMETERS:
!
    LOGICAL,          INTENT(IN   )  :: am_I_Root  ! root CPU?
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(INOUT)  :: RC         ! Failure or success
!
! !REMARKS
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

    ! Return w/ success
    RC = GIGC_SUCCESS

  END SUBROUTINE CLEANUP_PASSIVE_TRACER
!EOC
END MODULE PASSIVE_TRACER_MOD 
