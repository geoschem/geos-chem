#if defined (ESMF_)
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: gc_initialization_mod
!
! !DESCRIPTION: Module GC\_INITIALIZATION\_MOD is the module that contains 
!  the GEOS-Chem initialize methods for the ESMF framework.
!\\
!\\
! !INTERFACE: 
!      
MODULE GC_INITIALIZATION_MOD
!
! !USES:
!      
  USE GC_TYPE_MOD,        ONLY : GC_MET_LOCAL
  USE GC_TYPE2_MOD,       ONLY : CHEMSTATE

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!      

  PUBLIC :: GC_CHEMINIT
  PUBLIC :: GC_SETENV
  PUBLIC :: GC_INITRUN
  PUBLIC :: GC_GETOPTS
  PUBLIC :: GC_INIT_DIMENSIONS

  ! Not sure if we need the SAVE statement here (bmy, 10/15/12)
!        SAVE 
!        TYPE(GC_MET_LOCAL) GC_MET
!        TYPE(CHEMSTATE)    GC_STATE

!BOC
CONTAINS
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: gc_cheminit
!
! !DESCRIPTION: Routine GC\_CHEMINIT initializes the GEOS-Chem chemistry
!  module so that it can connect to the ESMF interface.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE GC_CHEMINIT( PLON, PLAT )
!
! !USES:
!
    USE LOGICAL_MOD,   ONLY : DO_DIAG_WRITE
    USE CMN_SIZE_MOD,  ONLY : JJPAR, IIPAR
!
! !INPUT PARAMETERS: 
!
    INTEGER, INTENT(IN) :: PLON   ! Number of
    INTEGER, INTENT(IN) :: PLAT
! 
! !REVISION HISTORY: 
!  15 Oct 2012 - M. Long     - Initial version
!  15 Oct 2012 - R. Yantosca - Added ProTeX Headers, use F90 format/indents
!EOP
!------------------------------------------------------------------------------
!BOC

    ! Don't write diagnostic files
    DO_DIAG_WRITE = .FALSE. 

    !### Debug
    ! SET MODEL DIMENSIONS
    !IIPAR = 1
    !JJPAR = PCOLS

    !=======================================================================
    ! Read GEOS-Chem emissions data ets
    !=======================================================================
   
    ! For now comment this out
    !CALL GC_EMIS_INTI( PLONL, PLATL, 1 )

  END SUBROUTINE GC_CHEMINIT
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: gc_emis_inti
!
! !DESCRIPTION: Routine GC\_EMIS\_INTI establishes the emissions fields for
!   each set of geos-chem emissions used. 
!  component.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE GC_EMIS_INTI( PLONL, PLATL, PPLON )
!
! !USES:
!
!    USE GC_EMIS_EDGAR, ONLY: GC_EDGAR_EMIS_INTI
!
! !INPUT PARAMETERS: 
!
    INTEGER, INTENT(IN) :: PLONL
    INTEGER, INTENT(IN) :: PLATL
    INTEGER, INTENT(IN) :: PPLON
! 
! !REMARKS:
!  This is a stub routine for now.  Will eventually connect with the Emissions
!  Component as being developed by Christoph Keller.
!
! !REVISION HISTORY: 
!  15 Oct 2012 - M. Long     - Initial version
!  15 Oct 2012 - R. Yantosca - Added ProTeX Headers, use F90 format/indents
!EOP
!------------------------------------------------------------------------------
!BOC

    ! Comment out for now
    !CALL GC_EDGAR_EMIS_INTI( PLONL, PLATL, PPLON )

  END SUBROUTINE GC_EMIS_INTI
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: gc_initrun
!
! !DESCRIPTION: Routine GC\_INITRUN is the Initialize method for the ESMF
!  interface that connects GEOS-Chem to the GEOS-5 GCM.  Calls to the various
!  GEOS-Chem init routines (which allocate arrays, etc.) are made from here.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE GC_INITRUN( LOCAL_MET, CHEM_STATE, tsChem,    &
                         nymd,      nhms,       am_I_Root )
      
!
! !USES:
!
    USE GC_TYPE2_MOD,       ONLY : INIT_CHEMSTATE      
    USE CMN_SIZE_MOD,       ONLY : IIPAR
    USE CMN_SIZE_MOD,       ONLY : JJPAR
    USE CMN_SIZE_MOD,       ONLY : LLTROP
    USE COMODE_MOD,         ONLY : INIT_COMODE
    USE COMODE_LOOP_MOD
    USE GC_ENVIRONMENT_MOD, ONLY : ALLOCATE_ALL
    USE GC_ENVIRONMENT_MOD, ONLY : INIT_ALL
    !USE GC_ENVIRONMENT_MOD, ONLY : TRACER_INDEX
    !USE GC_ENVIRONMENT_MOD, ONLY : TRACER_NAMES
    USE GCKPP_COMODE_MOD,   ONLY : INIT_GCKPP_COMODE
    USE GRID_MOD,           ONLY : INIT_GRID
    USE DAO_MOD,            ONLY : INIT_DAO
    USE LOGICAL_MOD,        ONLY : LKPP
    USE LOGICAL_MOD,        ONLY : LPRT
    USE LOGICAL_MOD,        ONLY : LEMIS
    USE LOGICAL_MOD,        ONLY : LCHEM
    USE PBL_MIX_MOD,        ONLY : INIT_PBL_MIX
    USE PRESSURE_MOD,       ONLY : INIT_PRESSURE
    USE TRACER_MOD,         ONLY : ITS_A_FULLCHEM_SIM
    USE TRACER_MOD,         ONLY : ITS_AN_AEROSOL_SIM
    USE TRACER_MOD,         ONLY : ID_TRACER
    USE TRACER_MOD,         ONLY : TRACER_NAME
    USE TRACER_MOD,         ONLY : N_TRACERS
    USE TOMS_MOD,           ONLY : TO3_DAILY
    USE WETSCAV_MOD,        ONLY : INIT_WETSCAV

    ! Comment these out for now (bmy, 10/15/12)
    !USE TIME_MANAGER,       ONLY : GET_STEP_SIZE
    !USE GC_GRIDAREA
!
! !INPUT PARAMETERS: 
!
    REAL,               INTENT(IN)  :: tsChem       ! Chemistry timestep [s]
    INTEGER,            INTENT(IN)  :: nymd         ! GMT date (YYYY/MM/DD)
    INTEGER,            INTENT(IN)  :: nhms         ! GMT time (hh:mm:ss)
    LOGICAL,            INTENT(IN)  :: am_I_Root    ! Is this the root CPU?
!
! !OUTPUT PARAMETERS:
!
    TYPE(CHEMSTATE),    INTENT(OUT) :: CHEM_STATE   ! Obj for chemistry
    TYPE(GC_MET_LOCAL), INTENT(OUT) :: LOCAL_MET    ! Obj for meteorology
!
! !REMARKS
!  Add other calls to GEOS-Chem init routines as necessary
!
! !REVISION HISTORY: 
!  15 Oct 2012 - M. Long     - Initial version
!  15 Oct 2012 - R. Yantosca - Added ProTeX Headers, use F90 format/indents
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
! 
    INTEGER :: DTIME          ! CHEMISTY TIME STEP
    INTEGER :: K, RC, AS
    INTEGER :: PLON, PLAT

    !=======================================================================
    ! GC_INITRUN begins here!
    !=======================================================================

    ! Initialize
    RC     = 0
    DTIME = tsChem
       
    print*, '### am_I_Root: ', am_I_root

    ! Initialize the timing routines
    CALL GC_INIT_TIMEINTERFACE( DTIME, nymd, nhms, am_I_Root )
    IF ( LPRT .and. am_I_Root ) THEN
       print*, '### after gc_init_timeinterface'
    endif

    ! Allocate alll 
    CALL ALLOCATE_ALL

    ! Allocate
    CALL ALLOCATE_INTERFACE
              
    ! Read options from the GEOS-Chem input file "input.geos"
    CALL GC_GETOPTS( am_I_Root )

    ! Initialize derived-type objects for meteorology & chemistry states
    CALL INIT_ALL( LOCAL_MET, CHEM_STATE )

    ! Save tracer names and ID's into CHEM_STATE
    CHEM_STATE%TRAC_NAME(1:N_TRACERS) = TRACER_NAME
    CHEM_STATE%TRAC_ID(1:N_TRACERS)   = ID_TRACER

    ! Allocate and zero GEOS-Chem diagnostic arrays
    CALL NDXX_SETUP
       
    ! Initialize 
    CALL GC_CHEMINIT( PLON, PLAT )

    ! Allocate and initialize met field arrays
    CALL INIT_DAO
    
    ! Initialize the GEOS-Chem pressure module (set Ap & Bp)
    CALL INIT_PRESSURE( am_I_Root )
    ! IF ( LPRT ) CALL ENDRUN( '### MAIN: A INIT_PRESSURE' )

! We should block these off with IF statements
!    ! Initialize the PBL mixing module
!    CALL INIT_PBL_MIX
!       
!    ! Initialize the  
!    CALL INIT_WETSCAV

      ! INITIALIZE ALLOCATABLE SMVGEAR/KPP ARRAYS
    IF ( LEMIS .OR. LCHEM ) THEN
       IF ( ITS_A_FULLCHEM_SIM() ) CALL INIT_COMODE( am_I_Root )
       IF ( ITS_AN_AEROSOL_SIM() ) CALL INIT_COMODE( am_I_Root )
       IF ( LKPP ) THEN
          CALL INIT_GCKPP_COMODE( am_I_Root, IIPAR,   JJPAR, LLTROP,  &
                                  ITLOOP,    NMTRATE, IGAS,  RC      )
       ENDIF
    ENDIF
    
    ! Allocate array of overhead O3 columns for TOMS
    ALLOCATE( TO3_DAILY( IIPAR, JJPAR ), STAT=AS )
    TO3_DAILY = 0d0

    ! Initialize dry deposition (work in progress), add here

  END SUBROUTINE GC_INITRUN
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: gc_getopts
!
! !DESCRIPTION: Routine GC\_GETOPTS reads options for a GEOS-Chem simulation
!  from the input.geos.rc input file.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE GC_GETOPTS( am_I_Root )
!
! !USES:
!
    USE INPUT_MOD, ONLY : READ_INPUT_FILE
!
! !INPUT PARAMETERS: 
!
    LOGICAL, INTENT(IN) :: am_I_Root   ! Are we on the root CPU?
! 
! !REVISION HISTORY: 
!  15 Oct 2012 - M. Long     - Initial version
!  15 Oct 2012 - R. Yantosca - Added ProTeX Headers, use F90 format/indents
!EOP
!------------------------------------------------------------------------------
!BOC
    
    ! Read the GEOS-Chem input file here
    CALL READ_INPUT_FILE( am_I_Root )

  END SUBROUTINE GC_GETOPTS
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: gc_init_timeinterface
!
! !DESCRIPTION: Routine GC\_INIT\_TIMEINTERFACE intializes the starting date
!  of the run as well as the dynamic timestep
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE GC_INIT_TIMEINTERFACE( DT, nymd, nhms, am_I_Root )
!
! !USES:
!
    USE JULDAY_MOD,   ONLY : JULDAY, CALDATE
    USE TIME_MOD
!
! !INPUT PARAMETERS: 
!
    INTEGER, INTENT(IN) :: DT          ! Timestep [seconds]
    INTEGER, INTENT(IN) :: nymd        ! GMT date (YYYY/MM/DD)
    INTEGER, INTENT(IN) :: nhms        ! GMT time (hh:mm:ss)
    LOGICAL, INTENT(IN) :: am_I_Root   ! Are we on the root CPU
! 
! !REVISION HISTORY: 
!  15 Oct 2012 - M. Long     - Initial version
!  15 Oct 2012 - R. Yantosca - Added ProTeX Headers, use F90 format/indents
!  15 Oct 2012 - R. Yantosca - Need to pass am_I_Root to INITIALIZE from here
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER :: NYMDB, NHMSB, Y, M, D, TOD, H, MN, S
    INTEGER :: DTM
    REAL*8  :: FRAC_DAY
       
    ! Convert timestep from seconds to minutes
    DTM =  DT / 60

    ! Set GEOS-Chem timesteps (for now, set all to the same value)
    CALL SET_TIMESTEPS( am_I_Root, DTM, DTM, DTM, DTM, DTM, DTM )
       
    ! SET NYMDB & NHMSB
    Y    = nymd / 10000                    ! Year
    M    = ( nymd - Y*10000 ) /100         ! Month
    D    = nymd - Y*10000 - M*100          ! Day
                                           
    H    = nhms/10000                      ! Hour
    MN   = ( nhms - H*10000 ) /100         ! Minute
    S    = nhms - H*10000 - MN*100         ! Seconds
                                           
    TOD = H*3600 + MN*60 + S               ! Seconds since start of day


    ! Fraction of the day that has elapsed
    IF ( TOD .GT. 0 ) THEN
       FRAC_DAY = DBLE( D ) + DBLE( TOD ) / DBLE( 86400 )   
    ELSE
       FRAC_DAY = DBLE( D )
    ENDIF
    
    ! Find the starting date
    CALL CALDATE( JULDAY( Y, M, FRAC_DAY ), NYMDB, NHMSB )

    ! Set the starting date & time of the GEOS-Chem simulation
    CALL SET_BEGIN_TIME( NYMDB, NHMSB )

    ! Set the current time to the starting time
    CALL SET_CURRENT_TIME()

    ! INITIALIZE AND SET TIMESTEP COUNTS USING
    ! GEOS-CHEM'S INITIALIZE (SEE GEOS-CHEM'S MAIN.F)
    CALL INITIALIZE( 2, am_I_Root )
    CALL INITIALIZE( 3, am_I_Root )

  END SUBROUTINE GC_INIT_TIMEINTERFACE
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: gc_init_dimensions
!
! !DESCRIPTION: Routine GC\_INIT\_DIMENSIONS
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE GC_INIT_DIMENSIONS( NI, NJ, NL )
!
! !USES:
!
    USE CMN_SIZE_MOD, ONLY: IIPAR, JJPAR, LLPAR, IGLOB, JGLOB, LGLOB
!
! !INPUT PARAMETERS: 
!
    INTEGER, INTENT(IN) :: NI    ! Size of local I dimension
    INTEGER, INTENT(IN) :: NJ    ! Size of local J dimension
    INTEGER, INTENT(IN) :: NL    ! Number of levels
! 
! !REVISION HISTORY: 
!  15 Oct 2012 - M. Long     - Initial version
!  15 Oct 2012 - R. Yantosca - Added ProTeX Headers, use F90 format/indents
!EOP
!------------------------------------------------------------------------------
!BOC
    
    ! Set GEOS-Chem size variables from the locally-defined 
    ! dimensions returned by the ESMF framework
    IGLOB = NI
    JGLOB = NJ
    LGLOB = NL
    IIPAR = NI
    JJPAR = NJ
    LLPAR = NL

  END SUBROUTINE GC_INIT_DIMENSIONS
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ALLOCATE_INTERFACE
!
! !DESCRIPTION: Routine ALLOCATE\_INTERFACE allocates GEOS-Chem arrays that
!  are not allocated elsewhere.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE ALLOCATE_INTERFACE
!
! !USES:
!
    USE CMN_SIZE_MOD, ONLY : IIPAR, JJPAR, LLPAR
    USE ERROR_MOD,    ONLY : ALLOC_ERR
    USE GRID_MOD,     ONLY : INIT_GRID
    USE UVALBEDO_MOD, ONLY : UVALBEDO
    USE DAO_MOD,      ONLY : TO3
! 
! !REVISION HISTORY: 
!  15 Oct 2012 - M. Long     - Initial version
!  15 Oct 2012 - R. Yantosca - Added ProTeX Headers, use F90 format/indents
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER :: AS

    ! UV albedo for photolysis
    ALLOCATE( UVALBEDO( IIPAR, JJPAR), STAT=AS  )

  END SUBROUTINE ALLOCATE_INTERFACE
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: gc_setenv
!
! !DESCRIPTION: GC\_SETENV
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE GC_SETENV
!
! !REMARKS:
!  This is a stub routine for now

! !REVISION HISTORY: 
!  15 Oct 2012 - M. Long     - Initial version
!  15 Oct 2012 - R. Yantosca - Added ProTeX Headers, use F90 format/indents
!EOP
!------------------------------------------------------------------------------
!BOC
  END SUBROUTINE GC_SETENV
!EOC      
END MODULE GC_INITIALIZATION_MOD
#endif
