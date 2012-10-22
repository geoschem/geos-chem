#if defined (ESMF_)
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: gc_chunk_mod
!
! !DESCRIPTION: Module GC\_CHUNK\_MOD is the module that contains 
!  the GEOS-Chem chunk code init, run and finalize methods.
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

        PUBLIC :: GC_CHEMINIT, GC_SETENV
        PUBLIC :: GC_INITRUN, GC_GETOPTS, GC_INIT_DIMENSIONS
!        PUBLIC :: GC_MET, GC_STATE

        SAVE 

!        TYPE(GC_MET_LOCAL) GC_MET
!        TYPE(CHEMSTATE)    GC_STATE

      CONTAINS

        SUBROUTINE GC_CHEMINIT( PLON, PLAT )
!-----------------------------------------------------------------------
! 	... CHEMISTRY MODULE INTIALIZATION
!-----------------------------------------------------------------------

        USE LOGICAL_MOD,   ONLY : DO_DIAG_WRITE
        USE CMN_SIZE_MOD,  ONLY : JJPAR, IIPAR

        IMPLICIT NONE

!-----------------------------------------------------------------------
! 	... DUMMY ARGUMENTS
!-----------------------------------------------------------------------
        INTEGER, INTENT(IN) :: PLON, PLAT

!-----------------------------------------------------------------------
!       INITIALIZE GEOS-CHEM RUNTIME ENVIRONMENT
!-----------------------------------------------------------------------

        DO_DIAG_WRITE = .FALSE. ! DON'T WANT TO WRITE DIAGNOSTICS FILES!

        ! SET MODEL DIMENSIONS
!        IIPAR = 1
!        JJPAR = PCOLS

!-----------------------------------------------------------------------
! 	... READ GEOS-CHEM EMISSIONS DATASETS 
!-----------------------------------------------------------------------
!        CALL GC_EMIS_INTI( PLONL, PLATL, 1 )
!-----------------------------------------------------------------------
      END SUBROUTINE GC_CHEMINIT
      
      SUBROUTINE GC_EMIS_INTI( PLONL, PLATL, PPLON )
        
!      	USE GC_EMIS_EDGAR, ONLY: GC_EDGAR_EMIS_INTI
        
      	IMPLICIT NONE

      	INTEGER, INTENT(IN) :: PLONL, PLATL, PPLON
!----------------------------------------------------------------------
! THIS SECTION WILL ESTABLISH THE EMISSIONS FIELDS FOR EACH SET
! OF GEOS-CHEM EMISSIONS USED. <<USE GC NAMELIST OR MODIFY BCC'S?>>
!----------------------------------------------------------------------
!      	CALL GC_EDGAR_EMIS_INTI( PLONL, PLATL, PPLON )

     END SUBROUTINE GC_EMIS_INTI
     
     SUBROUTINE GC_INITRUN(LOCAL_MET, CHEM_STATE, tsChem, nymd, nhms, am_I_Root)
       
       USE CMN_SIZE_MOD,       ONLY : IIPAR, JJPAR, LLTROP
       USE GC_ENVIRONMENT_MOD, ONLY : ALLOCATE_ALL, INIT_ALL!, TRACER_INDEX, TRACER_NAMES
       USE GRID_MOD,           ONLY : INIT_GRID !AREA_M2, AREA_CM2
!       USE TIME_MANAGER,       ONLY : GET_STEP_SIZE
!       USE GC_GRIDAREA
       USE DAO_MOD,            ONLY : INIT_DAO
       USE COMODE_MOD,         ONLY : INIT_COMODE
       USE PRESSURE_MOD,       ONLY : INIT_PRESSURE
       USE GCKPP_COMODE_MOD,   ONLY : INIT_GCKPP_COMODE
       USE LOGICAL_MOD,        ONLY : LKPP, LPRT, LEMIS, LCHEM
       USE TRACER_MOD,         ONLY : ITS_A_FULLCHEM_SIM, ITS_AN_AEROSOL_SIM, &
                                      ID_TRACER, TRACER_NAME, N_TRACERS
       USE PBL_MIX_MOD,        ONLY : INIT_PBL_MIX

       USE GC_TYPE2_MOD,       ONLY : INIT_CHEMSTATE
       
       USE COMODE_LOOP_MOD

       USE TOMS_MOD,           ONLY : TO3_DAILY
       
       USE WETSCAV_MOD,        ONLY : INIT_WETSCAV

       IMPLICIT NONE
       

       INTEGER             :: DTIME ! CHEMISTY TIME STEP
       INTEGER             :: K, RC, AS
       INTEGER             :: PLON, PLAT

       REAL, INTENT(IN)    :: tsChem
       INTEGER, INTENT(IN) :: nymd, nhms
       LOGICAL, INTENT(IN) :: am_I_Root
       TYPE(CHEMSTATE), INTENT(OUT)    :: CHEM_STATE
       TYPE(GC_MET_LOCAL), INTENT(OUT) :: LOCAL_MET
       
!---------------------------------------------------------------------

       RC = 0
       DTIME = tsChem
       
       CALL GC_INIT_TIMEINTERFACE(DTIME, nymd, nhms, am_I_Root)
       
       CALL ALLOCATE_ALL

       CALL ALLOCATE_INTERFACE
              
       CALL GC_GETOPTS(am_I_Root)

       CALL INIT_ALL(LOCAL_MET, CHEM_STATE)

       CHEM_STATE%TRAC_NAME(1:N_TRACERS) = TRACER_NAME
       CHEM_STATE%TRAC_ID(1:N_TRACERS)   = ID_TRACER

       CALL NDXX_SETUP

       CALL GC_CHEMINIT( PLON, PLAT )

      ! GEOS-CHEM INITIALIZATION ROUTINES
       CALL INIT_DAO

      ! INITIALIZE THE NEW HYBRID PRESSURE MODULE.  DEFINE AP AND BP.
       CALL INIT_PRESSURE(am_I_Root)
!       IF ( LPRT ) CALL ENDRUN( '### MAIN: A INIT_PRESSURE' )

       CALL INIT_PBL_MIX

       CALL INIT_WETSCAV

      ! INITIALIZE ALLOCATABLE SMVGEAR/KPP ARRAYS
       IF ( LEMIS .OR. LCHEM ) THEN
          IF ( ITS_A_FULLCHEM_SIM() ) CALL INIT_COMODE(am_I_Root)
          IF ( ITS_AN_AEROSOL_SIM() ) CALL INIT_COMODE(am_I_Root)
          IF ( LKPP ) THEN
             CALL INIT_GCKPP_COMODE( am_I_Root, IIPAR,   JJPAR, LLTROP,  &
                                     ITLOOP,    NMTRATE, IGAS,  RC      )
          ENDIF
       ENDIF

       ALLOCATE( TO3_DAILY( IIPAR, JJPAR ), STAT=AS )
       TO3_DAILY = 0d0

!---------------------------------------------------------------------
!     	... INITIALIZE DEPOSITION VELOCITIES
!--------------------------------------------------------------------- 

     END SUBROUTINE GC_INITRUN

     SUBROUTINE GC_GETOPTS(am_I_Root)
       
       USE INPUT_MOD,   ONLY : READ_INPUT_FILE
       LOGICAL, INTENT(IN) :: am_I_Root

       CALL READ_INPUT_FILE(am_I_Root)

     END SUBROUTINE GC_GETOPTS

     SUBROUTINE GC_INIT_TIMEINTERFACE(DT, nymd, nhms, am_I_Root)

       USE JULDAY_MOD,   ONLY : JULDAY, CALDATE
       USE TIME_MOD
       
       IMPLICIT NONE
       
       INTEGER             :: NYMDB, NHMSB, Y, M, D, TOD, H, MN, S
       INTEGER, INTENT(IN) :: DT, nymd, nhms
       INTEGER             :: DTM
       REAL*8              :: FRAC_DAY
       LOGICAL, INTENT(IN) :: am_I_Root
       
        ! SET GEOS-CHEM TIMESTEPS
       DTM = DT/60 ! TIMESTEP FROM SEC -> MIN
       CALL SET_TIMESTEPS(am_I_Root, DTM, DTM, DTM, DTM, DTM, DTM)
       
        ! SET NYMDB & NHMSB

       Y = nymd/10000
       M = (nymd-Y*10000)/100
       D = nymd-Y*10000-M*100

       H = nhms/10000
       MN = (nhms-H*10000)/100
       S = nhms-H*10000-MN*100

       TOD = H*3600+MN*60+S

       IF ( TOD .GT. 0 ) THEN
          FRAC_DAY = DBLE(D) + DBLE(TOD)/DBLE(86400)
       ELSE
          FRAC_DAY = DBLE(D)
       ENDIF

       CALL CALDATE( JULDAY( Y, M, FRAC_DAY ), NYMDB, NHMSB )

        ! SET BEGIN & CURRENT TIME VALUES FROM CAM VALUES
       
       CALL SET_BEGIN_TIME( NYMDB, NHMSB )

       CALL SET_CURRENT_TIME()

        ! INITIALIZE AND SET TIMESTEP COUNTS USING
        ! GEOS-CHEM'S INITIALIZE (SEE GEOS-CHEM'S MAIN.F)
       CALL INITIALIZE( 2 )
       CALL INITIALIZE( 3 )

     END SUBROUTINE GC_INIT_TIMEINTERFACE

     SUBROUTINE GC_INIT_DIMENSIONS(NI,NJ,NL)
       
       USE CMN_SIZE_MOD, ONLY: IIPAR,JJPAR,IGLOB,JGLOB, &
            LGLOB, LLPAR

       INTEGER, INTENT(IN) :: NI ! Size of local I dimension
       INTEGER, INTENT(IN) :: NJ ! Size of local J dimension
       INTEGER, INTENT(IN) :: NL ! Number of levels
       ! SET MODEL DIMENSIONS

       IGLOB = NI
       JGLOB = NJ
       LGLOB = NL
       
       IIPAR = NI
       JJPAR = NJ
       LLPAR = NL

     END SUBROUTINE GC_INIT_DIMENSIONS

     SUBROUTINE ALLOCATE_INTERFACE

       USE CMN_SIZE_MOD, ONLY : IIPAR, JJPAR, LLPAR
       USE ERROR_MOD,    ONLY : ALLOC_ERR
       USE GRID_MOD,     ONLY : INIT_GRID!AREA_M2, AREA_CM2, INIT_GRID
       USE UVALBEDO_MOD, ONLY : UVALBEDO
       USE DAO_MOD,      ONLY : TO3
       
       IMPLICIT NONE

       INTEGER AS

       ALLOCATE(UVALBEDO(IIPAR,JJPAR))
!       ALLOCATE(TO3(IIPAR,JJPAR))

     END  SUBROUTINE ALLOCATE_INTERFACE

     SUBROUTINE GC_SETENV
     END SUBROUTINE GC_SETENV
      
   END MODULE GC_INITIALIZATION_MOD
#endif
