! $Id: biomass_mod.f,v 1.13 2006/10/17 17:51:07 bmy Exp $
      MODULE BIOMASS_MOD
!
!******************************************************************************
!  Module BIOMASS_MOD is a "wrapper" module, which allows us to select either
!  GFED2 biomass burning emissions, or the default GEOS-Chem biomass burning
!  emissions (based on Bryan Duncan et al).  (psk, bmy, 4/5/06, 9/28/06)
!
!  GEOS-Chem has the following biomass burning gas-phase species:
!
!  Species   Index   G-C Tracer #          Units
!  ----------------------------------------------------------------------------
!  GAS PHASE SPECIES (contained in both GFED2 & Duncan et al 2001)
!
!   NOx        1          1          [molec NOx /cm2/s]
!   CO         2          4          [molec CO  /cm2/s]
!   ALK4       3          5          [atoms C   /cm2/s]
!   ACET       4          9          [atoms C   /cm2/s]
!   MEK        5          10         [atoms C   /cm2/s]
!   ALD2       6          11         [atoms C   /cm2/s]
!   PRPE       7          18         [atoms C   /cm2/s]
!   C3H8       8          19         [atoms C   /cm2/s]
!   CH2O       9          20         [molec CH2O/cm2/s]
!   C2H6       10         21         [atoms C   /cm2/s]
!
!  ----------------------------------------------------------------------------
!  AEROSOL SPECIES (contained in GFED2; read separately in Duncan et al 2001)  
!
!   SO2        11         26         [molec SO2 /cm2/s]
!   NH3        12         32         [molec NH3 /cm2/s]
!   BC         13         34         [atoms C   /cm2/s]
!   OC         14         35         [atoms C   /cm2/s]
!
!  ----------------------------------------------------------------------------
!  FOR CO2 SIMULATION ONLY
!
!   CO2        15         1          [molec CO2 /cm2/s]
!
!
!  Module Variables:
!  ============================================================================
!  (1 ) BIOMASS      (REAL*8 )    : Biomass emissions [molec/cm3/s]
!  (2 ) BIOMASS_SAVE (REAL*8 )    : Internal array for biomass emissions 
!  (3 ) BIOTRCE      (INTEGER)    : Index array tracer #'s for biomass species
!  (4 ) IDBNOX       (INTEGER)    : Index for NOx  in BIOMASS, BIOMASS_SAVE
!  (5 ) IDBCO        (INTEGER)    : Index for CO   in BIOMASS,c BIOMASS_SAVE
!  (6 ) IDBC2H6      (INTEGER)    : Index for C2H6 in BIOMASS, BIOMASS_SAVE
!  (7 ) NBIOMAX      (INTEGER)    : Number of biomass burning species
!  (8 ) NBIOMAX_GAS  (INTEGER)    : Number of gas-phase biomass burning species
!
!  Module Routines:
!  ============================================================================
!  (1 ) COMPUTE_BIOMASS_EMISSIONS : Gets biomass emissions; updates diagnostics
!  (2 ) INIT_BIOMASS              : Allocates & zeroes module arrays
!  (3 ) CLEANUP_BIOMASS           : Deallocates module arrays
! 
!  GEOS-Chem modules referenced by "biomass_mod.f"
!  ============================================================================
!  (1 ) bpch2_mod.f               : Module w/ routines for bpch file I/O
!  (2 ) dao_mod.f                 : Module w/ arrays for DAO met fields
!  (3 ) diag_mod.f                : Module w/ GEOS-CHEM diagnostic arrays
!  (4 ) directory_mod.f           : Module w/ GEOS-CHEM data & met field dirs
!  (5 ) error_mod.f               : Module w/ I/O error and NaN check routines
!  (6 ) gc_biomass_mod.f          : Module w/ routines for default G-C biomass
!  (7 ) gfed2_biomass_mod.f       : Module w/ routines for GFED2 biomass 
!  (8 ) grid_mod.f                : Module w/ horizontal grid information
!  (9 ) logical_mod.f             : Module w/ GEOS-CHEM logical switches
!  (10) time_mod.f                : Module w/ routines for computing time/ date
!
!  NOTES:  
!  (1 ) Rewrote so that all 15 biomass species (from either GFED2 or Duncan
!        et al 2001) are contained in the BIOMASS array.  Also removed the
!        BIOMASS_SAVE array because we no longer need to convert the data
!        to [molec/cm3/s] on each timestep (bmy, 9/28/06)
!******************************************************************************
!
      IMPLICIT NONE

      !=================================================================
      ! MODULE PRIVATE DECLARATIONS -- keep certain internal variables 
      ! and routines from being seen outside "biomass_mod.f"
      !=================================================================

      ! Make everything PRIVATE ...
      PRIVATE

      ! ... except these variables
      PUBLIC :: NBIOMAX
      PUBLIC :: NBIOMAX_GAS
      PUBLIC :: BIOMASS
      PUBLIC :: BIOTRCE
      PUBLIC :: IDBBC
      PUBLIC :: IDBCO
      PUBLIC :: IDBCO2
      PUBLIC :: IDBC2H6
      PUBLIC :: IDBNH3
      PUBLIC :: IDBNOX
      PUBLIC :: IDBOC
      PUBLIC :: IDBSO2

      ! ... and these routines
      PUBLIC :: CLEANUP_BIOMASS
      PUBLIC :: COMPUTE_BIOMASS_EMISSIONS

      !=================================================================
      ! MODULE VARIABLES
      !=================================================================

      ! Parameters
      INTEGER, PARAMETER   :: NBIOMAX     = 15
      INTEGER, PARAMETER   :: NBIOMAX_GAS = 10
      INTEGER, PARAMETER   :: IDBNOX      = 1
      INTEGER, PARAMETER   :: IDBCO       = 2
      INTEGER, PARAMETER   :: IDBC2H6     = 10
      INTEGER, PARAMETER   :: IDBSO2      = 11
      INTEGER, PARAMETER   :: IDBNH3      = 12
      INTEGER, PARAMETER   :: IDBBC       = 13
      INTEGER, PARAMETER   :: IDBOC       = 14
      INTEGER, PARAMETER   :: IDBCO2      = 15
      
      ! Arrays
      INTEGER              :: BIOTRCE(NBIOMAX)
      REAL*8,  ALLOCATABLE :: BIOMASS(:,:,:)
      !-------------------------------------------------------------------
      ! Prior to 9/28/06:
      ! BIOMASS_SAVE is redundant now.  We only used to keep this here
      ! because we had to preserve biomass in [molec/cm2/s] locally, so
      ! that we could compute [molec/cm3/s] on each timestep.  But now
      ! we no longer do that unit conversion, so get rid of this array.
      ! (bmy, 9/28/06)
      !REAL*8,  ALLOCATABLE :: BIOMASS_SAVE(:,:,:)
      !-------------------------------------------------------------------

      !=================================================================
      ! MODULE ROUTINES -- follow below the "CONTAINS" statement 
      !=================================================================
      CONTAINS

!------------------------------------------------------------------------------

      SUBROUTINE COMPUTE_BIOMASS_EMISSIONS( YEAR, MONTH )
!
!******************************************************************************
!  Subroutine COMPUTE_BIOMASS_EMISSIONS is a wrapper which allows us to select
!  either the GFED2 biomass burning emissions, or the regular GEOS-Chem
!  biomass burning emissions (Duncan et al 2001). (psk, bmy, 4/5/06, 9/28/06)
!
!  This routine is called on each timestep.  At the start of a new month,
!  new biomass burning emissions are read from disk.  The ND28, ND29, ND32
!  diagnostics are updated on each timestep.
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) YEAR  (INTEGER) : Current year  
!  (2 ) MONTH (INTEGER) : Current month (1-12)
!
!  NOTES:
!  (1 ) Now store all biomass species in BIOMASS, from GFED2 or Duncan et al 
!        2001.  Also remove obsolete BIOMASS_SAVE array. (bmy, 9/28/06)
!******************************************************************************
!
      ! References to F90 modules
      USE DIAG_MOD,          ONLY : AD28, AD29, AD32_bb
      USE GC_BIOMASS_MOD,    ONLY : GC_COMPUTE_BIOMASS
      USE GC_BIOMASS_MOD,    ONLY : GC_READ_BIOMASS_BCOC
      USE GC_BIOMASS_MOD,    ONLY : GC_READ_BIOMASS_CO2
      USE GC_BIOMASS_MOD,    ONLY : GC_READ_BIOMASS_NH3
      USE GC_BIOMASS_MOD,    ONLY : GC_READ_BIOMASS_SO2
      USE GFED2_BIOMASS_MOD, ONLY : GFED2_COMPUTE_BIOMASS
      USE LOGICAL_MOD,       ONLY : LBIOMASS, LGFED2BB
      USE TIME_MOD,          ONLY : ITS_A_NEW_MONTH
      USE TRACER_MOD,        ONLY : ITS_A_CO2_SIM
      USE TRACERID_MOD,      ONLY : IDTBCPO, IDTNH3, IDTOCPO, IDTSO2

#     include "CMN_SIZE"          ! Size parameters
#     include "CMN_DIAG"          ! Diagnostic flags

      ! Arguments
      INTEGER, INTENT(IN)        :: YEAR, MONTH

      ! Local variables
      LOGICAL, SAVE              :: FIRST = .TRUE.
      LOGICAL                    :: DO_ND28, DO_ND29, DO_ND32
      INTEGER                    :: I,       J,       N,      N_BIOB
      REAL*8                     :: BXHT_CM, DTSRCE
      
      !=================================================================
      ! COMPUTE_BIOMASS_EMISSIONS begins here!
      !=================================================================

      ! If there are biomass emissions ...
      IF ( LBIOMASS ) THEN

         ! First-time initialization
         IF ( FIRST ) THEN
            CALL INIT_BIOMASS
            FIRST = .FALSE.
         ENDIF
         
         ! Define diagnostic flags
         DO_ND28 = ( ND28 > 0 )
         DO_ND29 = ( ND29 > 0 )
         DO_ND32 = ( ND32 > 0 )

         !==============================================================
         ! Read biomass emissions at the start of a new month
         !==============================================================
         IF ( ITS_A_NEW_MONTH() ) THEN

            ! Zero the array for biomass burning
            BIOMASS = 0d0

            ! Test for type of biomass emissions
            IF ( LGFED2BB ) THEN
            
               !---------------------------------
               ! GFED2 biomass inventory for
               ! gas-phase, aerosols, and CO2
               !---------------------------------

               ! Get emissions [molec/cm2/s] or [atoms C/cm2/s]
               CALL GFED2_COMPUTE_BIOMASS( YEAR, MONTH, BIOMASS )

            ELSE

               ! Test if it's a CO2 simulation
               IF ( ITS_A_CO2_SIM() ) THEN

                  !------------------------------
                  ! CO2 emissions (based on 
                  ! Duncan et al 2001 CO)
                  !------------------------------
                  
                  ! Get CO2 emissions [molec/cm2/s]
                  CALL GC_READ_BIOMASS_CO2( YEAR, MONTH,
     &                                      BIOMASS(:,:,IDBCO2) )
               ELSE

                  !------------------------------
                  ! Default GEOS-Chem inventory
                  ! (Bryan Duncan et al 2001)
                  !------------------------------

                  ! Get emissions of gas-phase species
                  ! in [molec/cm2/s] or [atoms C/cm2/s]
                  CALL GC_COMPUTE_BIOMASS( YEAR, MONTH, 
     &                                     BIOMASS(:,:,1:NBIOMAX_GAS) )

                  ! Get biomass SO2 [molec/cm2/s]
                  IF ( IDTSO2 > 0 ) THEN
                     CALL GC_READ_BIOMASS_SO2( YEAR, MONTH, 
     &                                         BIOMASS(:,:,IDBSO2) )
                  ENDIF

                  ! Get biomass NH3 [molec/cm2/s]
                  IF ( IDTNH3 > 0 ) THEN
                     CALL GC_READ_BIOMASS_NH3( YEAR, MONTH,
     &                                         BIOMASS(:,:,IDBNH3) )
                  ENDIF

                  ! Get biomass BC & OC [molec/cm2/s]
                  IF ( IDTBCPO > 0 .and. IDTOCPO > 0 ) THEN
                     CALL GC_READ_BIOMASS_BCOC( YEAR, MONTH,
     &                                          BIOMASS(:,:,IDBBC), 
     &                                          BIOMASS(:,:,IDBOC) ) 
                  ENDIF
               ENDIF
            ENDIF
         ENDIF

         !==============================================================
         ! Do the following on every timestep:
         !
         ! ND28, ND29, ND32 diags [molec/cm2/s] or [atoms C/cm2/s] 
         !==============================================================
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, N )
         DO N = 1, NBIOMAX
         DO J = 1, JJPAR
         DO I = 1, IIPAR

            ! ND28: biomass emissions in [molec/cm2/s]
            IF ( DO_ND28 ) THEN 
               AD28(I,J,N)  = AD28(I,J,N)  + BIOMASS(I,J,N)
            ENDIF
            
            ! ND29: CO biomass emissions [molec/cm2/s]
            IF ( DO_ND29 .and. N == IDBCO ) THEN 
               AD29(I,J,2)  = AD29(I,J,2)  + BIOMASS(I,J,IDBCO)
            ENDIF
            
            ! ND32: NOx biomass emissions in [molec/cm2/s]
            IF ( DO_ND32 .and. N == IDBNOx ) THEN
               AD32_bb(I,J) = AD32_bb(I,J) + BIOMASS(I,J,IDBNOx)
            ENDIF

         ENDDO
         ENDDO
         ENDDO
!$OMP END PARALLEL DO
      ENDIF

      ! Return to calling program
      END SUBROUTINE COMPUTE_BIOMASS_EMISSIONS

!------------------------------------------------------------------------------

      SUBROUTINE INIT_BIOMASS
!
!******************************************************************************
!  Subroutine INIT_BIOMASS allocates and zeroes the module arrays.
!  (bmy, 4/5/06, 9/28/06)
!
!  NOTES:
!  (1 ) Now set BIOTRCE for 15 biomass species (bmy, 9/28/06)
!  (2 ) Now remove BIOMASS_SAVE array, it's redundant (bmy, 9/28/06)
!******************************************************************************
!
      ! References to F90 modules
      USE ERROR_MOD,   ONLY : ALLOC_ERR
      USE LOGICAL_MOD, ONLY : LBIOMASS

#     include "CMN_SIZE"    ! Size parameters

      INTEGER              :: AS

      !=================================================================
      ! INIT_BIOMASS begins here!
      !=================================================================

      ! If there are biomass emissions ...
      IF ( LBIOMASS ) THEN

         ! Tracer numbers for each biomass species (CO2 is last)
         BIOTRCE(:) = (/ 1,  4,  5,  9,  10, 11, 18, 
     &                   19, 20, 21, 26, 30, 34, 35, 1/)
         
         ! Allocate array to hold monthly biomass emissions
         ALLOCATE( BIOMASS( IIPAR, JJPAR, NBIOMAX ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'BIOMASS' )
         BIOMASS = 0d0
         
         !-------------------------------------------------------------------
         ! Prior to 9/28/06:
         !! Allocate array to hold monthly biomass emissions
         !ALLOCATE( BIOMASS_SAVE( IIPAR, JJPAR, NBIOMAX ), STAT=AS )
         !IF ( AS /= 0 ) CALL ALLOC_ERR( 'BIOMASS_SAVE' )
         !BIOMASS_SAVE = 0d0
         !-------------------------------------------------------------------

      ENDIF
      
      ! Return to calling program
      END SUBROUTINE INIT_BIOMASS

!------------------------------------------------------------------------------

      SUBROUTINE CLEANUP_BIOMASS
!
!******************************************************************************
!  Subroutine CLEANUP_BIOMASS deallocates all module arrays.
!  (psk, bmy, 4/5/06, 9/28/06)
!
!  NOTES:
!  (1 ) Now remove BIOMASS_SAVE array, it's redundant (bmy, 9/28/06)
!******************************************************************************
!     
      !=================================================================
      ! CLEANUP_BIOMASS begins here!
      !=================================================================
      IF ( ALLOCATED( BIOMASS ) ) DEALLOCATE( BIOMASS )
      !----------------------------------------------------------------------
      ! Prior to 9/28/06:
      !IF ( ALLOCATED( BIOMASS_SAVE ) ) DEALLOCATE( BIOMASS_SAVE )      
      !----------------------------------------------------------------------

      ! Return to calling program
      END SUBROUTINE CLEANUP_BIOMASS

!------------------------------------------------------------------------------
      
      ! End of module
      END MODULE BIOMASS_MOD
