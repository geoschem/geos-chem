! $Id: biomass_mod.f,v 1.12 2006/04/21 15:39:51 bmy Exp $
      MODULE BIOMASS_MOD
!
!******************************************************************************
!  Module BIOMASS_MOD is a "wrapper" module, which allows us to select either
!  GFED2 biomass burning emissions, or the default GEOS-Chem biomass burning
!  emissions (based on Bryan Duncan et al).  (psk, bmy, 4/5/06)
!
!  GEOS-CHEM has the following biomass burning gas-phase species:
!
!     Species   Index   G-C Tracer #   
!     -------------------------------
!      NOx        1          1          
!      CO         2          4         
!      ALK4       3          5         
!      ACET       4          9         
!      MEK        5          10        
!      ALD2       6          11        
!      PRPE       7          18        
!      C3H8       8          19        
!      CH2O       9          20        
!      C2H6       10         21    
!
!  NOTE: Aerosol biomass species are read in separately in "carbon_mod.f"
!        and "sulfate_mod.f"
!
!  Module Variables:
!  ============================================================================
!  (1 ) BIOMASS      (REAL*8 )    : Biomass emissions [molec/cm3/s]
!  (2 ) BIOMASS_SAVE (REAL*8 )    : Internal array for biomass emissions 
!  (3 ) BIOTRCE      (INTEGER)    : Index array tracer #'s for biomass species
!  (4 ) IDBNOX       (INTEGER)    : Index for NOx  in BIOMASS, BIOMASS_SAVE
!  (5 ) IDBCO        (INTEGER)    : Index for CO   in BIOMASS, BIOMASS_SAVE
!  (6 ) IDBC2H6      (INTEGER)    : Index for C2H6 in BIOMASS, BIOMASS_SAVE
!  (7 ) NBIOMAX      (INTEGER)    : Number of biomass burning species
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
      PUBLIC :: BIOMASS
      PUBLIC :: BIOTRCE
      PUBLIC :: IDBNOX
      PUBLIC :: IDBCO
      PUBLIC :: IDBC2H6

      ! ... and these routines
      PUBLIC :: CLEANUP_BIOMASS
      PUBLIC :: COMPUTE_BIOMASS_EMISSIONS

      !=================================================================
      ! MODULE VARIABLES
      !=================================================================

      ! Parameters
      INTEGER, PARAMETER   :: NBIOMAX = 10
      INTEGER, PARAMETER   :: IDBNOX  = 1
      INTEGER, PARAMETER   :: IDBCO   = 2
      INTEGER, PARAMETER   :: IDBC2H6 = 10
      
      ! Arrays
      INTEGER              :: BIOTRCE(NBIOMAX)
      REAL*8,  ALLOCATABLE :: BIOMASS(:,:,:)
      REAL*8,  ALLOCATABLE :: BIOMASS_SAVE(:,:,:)

      !=================================================================
      ! MODULE ROUTINES -- follow below the "CONTAINS" statement 
      !=================================================================
      CONTAINS

!------------------------------------------------------------------------------

      SUBROUTINE COMPUTE_BIOMASS_EMISSIONS( YEAR, MONTH )
!
!******************************************************************************
!  Subroutine COMPUTE_BIOMASS_EMISSIONS is a wrapper which allows us to select
!  either the GFED2 biomass burning emissions, or the regular GEOS-CHEM
!  biomass burning emissions (from Bryan Duncan).  (psk, bmy, 4/5/06)
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
!******************************************************************************
!
      ! References to F90 modules
      USE DAO_MOD,           ONLY : BXHEIGHT
      USE DIAG_MOD,          ONLY : AD28, AD29, AD32_bb
      USE GC_BIOMASS_MOD,    ONLY : GC_COMPUTE_BIOMASS
      USE GFED2_BIOMASS_MOD, ONLY : GFED2_COMPUTE_BIOMASS
      USE LOGICAL_MOD,       ONLY : LBIOMASS, LGFED2BB
      USE TIME_MOD,          ONLY : ITS_A_NEW_MONTH, GET_TS_EMIS

#     include "CMN_SIZE"          ! Size parameters
#     include "CMN_DIAG"          ! Diagnostic flags

      ! Arguments
      INTEGER, INTENT(IN)        :: YEAR, MONTH

      ! Local variables
      LOGICAL, SAVE              :: FIRST = .TRUE.
      LOGICAL                    :: DO_ND28, DO_ND29, DO_ND32
      INTEGER                    :: I,       J,       N
      REAL*8                     :: BXHT_CM, DTSRCE
      
      !=================================================================
      ! COMPUTE_BIOMASS_EMISSIONS begins here!
      !=================================================================

      ! Emission timestep [s]
      DTSRCE = GET_TS_EMIS() * 60d0

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

         ! Read biomass emissions once per month
         IF ( ITS_A_NEW_MONTH() ) THEN

            ! Test for type of biomass emissions
            IF ( LGFED2BB ) THEN
            
               !------------------------------
               ! GFED2 biomass inventory
               !------------------------------

               ! Get emissions [molec/cm2/s] or [atoms C/cm2/s]
               CALL GFED2_COMPUTE_BIOMASS( YEAR, MONTH, BIOMASS_SAVE )

            ELSE

               !------------------------------
               ! Default GEOS-Chem inventory
               ! (based on Bryan Duncan)
               !------------------------------

               ! Get emissions [molec/cm2/s] or [atoms C/cm2/s]               
               CALL GC_COMPUTE_BIOMASS( YEAR, MONTH, BIOMASS_SAVE )

            ENDIF

         ENDIF

         !==============================================================
         ! Do the following on every timestep:
         !
         ! (1) Diagnostics     : ND28, ND29, ND32
         ! (2) Unit conversion : [  molec/cm2/s] --> [  molec/cm3/s] or
         !                       [atoms C/cm2/s] --> [atoms C/cm3/s]
         !==============================================================
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, N, BXHT_CM )
!$OMP+SCHEDULE( DYNAMIC )
         DO N = 1, NBIOMAX
         DO J = 1, JJPAR
         DO I = 1, IIPAR

            !--------------------------------------------------
            ! Archive Diagnostics: ND28, ND29, ND32
            ! NOTE: Store in molec/cm2, convert to molec/cm2/s
            !--------------------------------------------------

            ! ND28: biomass emissions in [molec/cm2/s]
            IF ( DO_ND28 ) THEN 
               AD28(I,J,N)  = AD28(I,J,N)  + BIOMASS_SAVE(I,J,N)
            ENDIF
            
            ! ND29: CO biomass emissions [molec/cm2/s]
            IF ( DO_ND29 .and. N == IDBCO ) THEN 
               AD29(I,J,2)  = AD29(I,J,2)  + BIOMASS_SAVE(I,J,IDBCO)
            ENDIF
            
            ! ND32: NOx biomass emissions in [molec/cm2/s]
            IF ( DO_ND32 .and. N == IDBNOx ) THEN
               AD32_bb(I,J) = AD32_bb(I,J) + BIOMASS_SAVE(I,J,IDBNOx)
            ENDIF

            !---------------------------------------------------
            ! Convert units to [molec/cm3/s] or [atoms C/cm3/s]
            !---------------------------------------------------

            ! Grid box height [cm]
            BXHT_CM         = BXHEIGHT(I,J,1) * 100d0

            ! Save [molec/cm3/s] or [atoms C/cm3/s] in BIOMASS
            ! for use in other routines
            BIOMASS(I,J,N)  = BIOMASS_SAVE(I,J,N) / BXHT_CM

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
!  (bmy, 4/5/06)
!
!  NOTES:
!******************************************************************************
!
      ! References to F90 modules
      USE ERROR_MOD,   ONLY : ALLOC_ERR
      USE LOGICAL_MOD, ONLY : LBIOMASS, LGFED2BB

#     include "CMN_SIZE"    ! Size parameters

      INTEGER              :: AS

      !=================================================================
      ! INIT_BIOMASS begins here!
      !=================================================================

      ! If there are biomass emissions ...
      IF ( LBIOMASS ) THEN

         ! Define BIOTRCE for backwards compatibility
         BIOTRCE(:) = (/ 1, 4, 5, 9, 10, 11, 18, 19, 20, 21/)

         ! Allocate array to hold monthly biomass emissions
         ALLOCATE( BIOMASS( IIPAR, JJPAR, NBIOMAX ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'BIOMASS' )
         BIOMASS = 0d0

         ! Allocate array to hold monthly biomass emissions
         ALLOCATE( BIOMASS_SAVE( IIPAR, JJPAR, NBIOMAX ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'BIOMASS_SAVE' )
         BIOMASS_SAVE = 0d0

      ENDIF
      
      ! Return to calling program
      END SUBROUTINE INIT_BIOMASS

!------------------------------------------------------------------------------

      SUBROUTINE CLEANUP_BIOMASS
!
!******************************************************************************
!  Subroutine CLEANUP_BIOMASS deallocates all module arrays (psk, bmy, 4/5/06)
!
!  NOTES:
!******************************************************************************
!     
      !=================================================================
      ! CLEANUP_BIOMASS begins here!
      !=================================================================
      IF ( ALLOCATED( BIOMASS      ) ) DEALLOCATE( BIOMASS      )
      IF ( ALLOCATED( BIOMASS_SAVE ) ) DEALLOCATE( BIOMASS_SAVE )      

      ! Return to calling program
      END SUBROUTINE CLEANUP_BIOMASS

!------------------------------------------------------------------------------
      
      ! End of module
      END MODULE BIOMASS_MOD
