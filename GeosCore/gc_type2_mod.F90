#if defined( DEVEL ) && ! defined ( ESMF_ )
! $Id: gc_type2_mod.F
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: gc_type2_mod
!
! !DESCRIPTION: Module GC\_TYPE2\_MOD defines the derived types that will represent
!               the "socket" with external modules. These types will be implemented
!               thru GEOS-Chem in a top-down format to tightly control memory
!               allocation and access both within and outside of GEOS-Chem.
!
!\\
!\\
!  NOTE: This is mostly for testing the grid-independent code in the current 
!  GEOS-Chem.  Many of these inputs will come from the GEOS-5 interface. It will
!  remain in DEVEL state for some time.
!\\
!\\
! !INTERFACE: 
!
      MODULE GC_TYPE2_MOD
!
! USES:
!
      IMPLICIT NONE
#     include "define.h"
      PUBLIC
      
      TYPE :: CHEMSTATE

        INTEGER, POINTER, DIMENSION(:)           :: &
             TRAC_ID, &  ! Tracer ID's set in TRACER_MOD
             SMVG_ID     ! Smvgear ID's set in TRACER_MOD
        CHARACTER(LEN=14), POINTER, DIMENSION(:) :: &
             TRAC_NAME ! Tracer names set in TRACER_MOD

        REAL*8, POINTER, DIMENSION(:,:,:,:)      :: &
             TRAC_TEND,  & ! Tracer  Tendency (<<units>>)
             TRAC_BTEND, & ! Biomass Tendency (<<units>>)
             TRACERS       ! Tracer concentration (Kg (per grid))

      END TYPE CHEMSTATE

 !-------------------------------------------------------------------------------

      TYPE PHYSTATE

         INTEGER                                   :: &
              BEGIN_I, &
              END_I,   &
              BEGIN_J, &
              END_J,   &
              NINDX,   &
              NJNDX
         REAL*8, ALLOCATABLE, DIMENSION(:)         :: &
              LAT,     &  ! LATITUDE (DEG)
              LON,     &  ! LONGITUDE (DEG)
              PS,      &  ! SURFACE PRESSURE
              PSDRY,   &  ! DRY SURFACE PRESSURE
              PHIS,    &  ! SURFACE GEOPOTENTIAL
              ULAT,    &  ! UNIQUE LATITUDES  (DEG)
              ULON        ! UNIQUE LONGITUDES (DEG)
         REAL*8, ALLOCATABLE, DIMENSION(:,:)       :: &
              T,       &  ! TEMPERATURE (K)
              U,       &  ! ZONAL WIND (M/S)
              V,       &  ! MERIDIONAL WIND (M/S)
              OMEGA,   &  ! VERTICAL PRESSURE VELOCITY (PA/S) 
              PMID,    &  ! MIDPOINT PRESSURE (PA) 
              PMIDDRY, &  ! MIDPOINT PRESSURE DRY (PA) 
              PDEL,    &  ! LAYER THICKNESS (PA)
              PDELDRY, &  ! LAYER THICKNESS DRY (PA)
              RPDEL,   &  ! RECIPROCAL OF LAYER THICKNESS (PA)
              RPDELDRY,&  ! RECIPRICOL LAYER THICKNESS DRY (PA)
              UZM,     &  ! ZONAL WIND FOR QBO (M/S)
              ZM          ! GEOPOTENTIAL HEIGHT ABOVE SURFACE AT MIDPOINTS (M)
         
!         REAL*8, DIMENSION(PCOLS,PVER+1)           :: &
!              PINT,    &  ! INTERFACE PRESSURE (PA)
!              PINTDRY, &  ! INTERFACE PRESSURE DRY (PA) 
!              LNPINT,  &  ! LN(PINT)
!              LNPINTDRY,& ! LOG INTERFACE PRESSURE DRY (PA) 
!              ZI          ! GEOPOTENTIAL HEIGHT ABOVE SURFACE AT INTERFACES (M)
         
!         INTEGER, DIMENSION(PCOLS) :: &
!              LATMAPBACK, &! MAP FROM COLUMN TO UNIQUE LAT FOR THAT COLUMN
!              LONMAPBACK, &! MAP FROM COLUMN TO UNIQUE LON FOR THAT COLUMN
!              CID         ! UNIQUE COLUMN ID
!         INTEGER :: ULATCNT, & ! NUMBER OF UNIQUE LATS IN CHUNK
!              ULONCNT     ! NUMBER OF UNIQUE LONS IN CHUNK

      END TYPE PHYSTATE
 !-----------------------------------------------------------------------

      TYPE(CHEMSTATE) :: CHEM_STATE
!      INTEGER         :: NULL

      REAL*8, PUBLIC, ALLOCATABLE :: EXT_STRATOH(:,:)
      REAL*8, PUBLIC, ALLOCATABLE :: EXT_SJVALUE(:,:,:)
      REAL*8, PUBLIC, ALLOCATABLE :: EXT_COPROD(:,:)
      REAL*8, PUBLIC, ALLOCATABLE :: EXT_COLOSS(:,:)

      CONTAINS

      INTEGER FUNCTION GET_TRACER_ID ( SPC_NAME )

 !-----------------------------------------------------------------------
 !     ... return overall species index associated with spc_name                                                                                                                
 !-----------------------------------------------------------------------   

      USE TRACER_MOD, ONLY : N_TRACERS 

      IMPLICIT NONE

    !-----------------------------------------------------------------------
    !     ... dummy arguments
    !-----------------------------------------------------------------------                                                                                                      

      CHARACTER(LEN=*), INTENT(IN) :: SPC_NAME

    !-----------------------------------------------------------------------
    !     ... local variables
    !-----------------------------------------------------------------------                                                                                                      

      INTEGER :: M

      GET_TRACER_ID = -1
      DO M = 1,N_TRACERS
!         IF( TRIM( SPC_NAME ) == TRIM( CHEM_STATE%TRAC_NAME(M) ) ) THEN
!            GET_TRACER_ID = M
!            EXIT
!         END IF
      END DO

      END FUNCTION GET_TRACER_ID

      SUBROUTINE INIT_CHEMSTATE(CHEM_STATE)
      
        USE CMN_SIZE_MOD,    ONLY : IIPAR, JJPAR, LLPAR, NBIOMAX, NNPAR
        USE COMODE_LOOP_MOD, ONLY : IGAS
        USE TRACER_MOD,      ONLY : N_TRACERS

        IMPLICIT NONE

        TYPE(CHEMSTATE), INTENT(out) :: CHEM_STATE

        ! 1-D ALLOCATIONS
        ALLOCATE( CHEM_STATE%TRAC_ID(N_TRACERS+1) )
        ALLOCATE( CHEM_STATE%SMVG_ID(IGAS) )

        ! 3-D ALLOCATIONS
        ALLOCATE( CHEM_STATE%TRAC_TEND(IIPAR,JJPAR,LLPAR,N_TRACERS+1) )
        ALLOCATE( CHEM_STATE%TRAC_BTEND(IIPAR,JJPAR,LLPAR,NBIOMAX) )
        ALLOCATE( CHEM_STATE%TRACERS(IIPAR,JJPAR,LLPAR,N_TRACERS+1) )

!        NULL = N_TRACERS+1 ! Dummy index just-in-case (e.g. CO2)

        CHEM_STATE%TRAC_ID = 0
        CHEM_STATE%SMVG_ID = 0

        CHEM_STATE%TRAC_TEND  = 0.
        CHEM_STATE%TRAC_BTEND = 0.
        CHEM_STATE%TRACERS    = 0.

        ALLOCATE(EXT_STRATOH(JJPAR, LLPAR))
        ALLOCATE(EXT_SJVALUE(JJPAR, LLPAR, 13)) ! NSPHOTO = 13
        ALLOCATE(EXT_COPROD( JJPAR, LLPAR))
        ALLOCATE(EXT_COLOSS( JJPAR, LLPAR))
        
      END SUBROUTINE INIT_CHEMSTATE

      END MODULE GC_TYPE2_MOD
#endif
!
!-----------------------------------------------------------------------
! FULLCHEM Simulation Emissions
! Done | Units? | Routine
! Y/P  --> Yes or Partially(fix needed)
!         good  --> Verifed that they are in Kg/s
!-----------------------------------------------------------------------
! Y    |        | CALL COMPUTE_BIOMASS_EMISSIONS
! Y    |        | CALL EMISS_STREETS_ANTHRO_05x0666
! Y    |        | CALL EMISS_STREETS_ANTHRO
! Y    |        | CALL EMISS_EDGAR( YEAR, MONTH )
! Y    |  good  | CALL EMISS_RETRO
! Y    |  good* | CALL EMISS_EPA_NEI
! Y    |        | CALL EMISS_VISTAS_ANTHRO
! Y    |        | CALL EMISS_BRAVO
! Y    |        | CALL EMISS_EMEP_05x0666
! Y    |        | CALL EMISS_EMEP
! Y    |        | CALL EMISS_CAC_ANTHRO_05x0666
! Y    |        | CALL EMISS_CAC_ANTHRO
! Y    |        | CALL EMISS_EPA_NEI
! Y    |        | CALL EMISS_NEI2005_ANTHRO_05x0666
! Y    |        | CALL EMISS_NEI2005_ANTHRO
! Y    |        | CALL EMISS_ARCTAS_SHIP( YEAR )
! Y    |        | CALL EMISS_ICOADS_SHIP
!      |        | CALL EMISSDR
! Y    |        | CALL EMISSSEASALT
! Y    |        | CALL EMISSSULFATE --> Be sure there's no PBL mixing
! Y    |        | CALL EMISSCARBON --> Be sure there's no PBL mixing
! Y    |        | CALL EMISSDUST --> Be sure there's no PBL mixing
! Y    |        | AIRCRAFT_NOX
! Y    |        | LIGHTNING_NOX
! Y    |        | SOIL_NOX
!      |        | BIOFUEL_BURN (NOx and CO)
!-----------------------------------------------------------------------
! Notes:
! LNLPBL Switch --> NEEDS TO BE ON (>=1)
!               --> But, does VDIFF need to be turned off? 
! NOT ALL EMISSIONS ARE JUST AT SFC, e.g. SO2
! LFUTURE --> How to deal with this? e.g. EDGAR emissions
! REGIONAL EMISSIONS OVERWRITE GLOBAL!!!!! DEAL WITH THIS!
!
! STT<-->CSPEC mapped in PARTITION
!
! KEEP EMISSIONS FROM UPDATING STT DIRECTLY
!
! NEI EMISSIONS: BIOFUEL EMISSIONS ARE NOT 'REALLY' BIOFUEL.
!                AS IN THERE'S NO IDBF'SPEC' INDEX.
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
! FULLCHEM Simulation Chemistry Routines
! Done | Units? | Routine
! Y/P  --> Yes or Partially(fix needed)
!         good  --> Verifed that they are in Kg
!-----------------------------------------------------------------------
!      |        | CALL CHEMDR
!      |        | CALL CHEMSEASALT
!      |        | CALL CHEMSULFATE
!      |        | CALL DO_ISOROPIAII
!      |        | CALL DO_RPMARES
!      |        | CALL CHEMCARBON
!      |        | CALL CHEMDUST
!      |        | CALL DRYFLX     
!      |        | CALL DIAGOH
!      |        | CALL OCEAN_SINK_ACET( STT(:,:,1,IDTACET) ) 
!-----------------------------------------------------------------------
! Notes:
! (1) If STT IS TIGHTLY LINKED TO CHEM_STATE, THEN THE ONLY CHANGE
!     NEEDED IN "DO_CHEMISTRY" IS TO PASS "CHEM_STATE" IN AND OUT
! (2) 1st PASS, WORKS ONLY WITH FULLCHEM, NOT APM OR ANY ADD-ON SIM
!     OPTIONS.
!-----------------------------------------------------------------------
