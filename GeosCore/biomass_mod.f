! $Id: biomass_mod.f,v 1.1 2009/09/16 14:06:39 bmy Exp $
      MODULE BIOMASS_MOD
!
!******************************************************************************
!  Module BIOMASS_MOD is a "wrapper" module, which allows us to select either
!  GFED2 biomass burning emissions, or the default GEOS-Chem biomass burning
!  emissions (based on Bryan Duncan et al).  (psk, bmy, 4/5/06, 9/18/07)
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
!   CO2        24         1          [molec CO2 /cm2/s]
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
!  (2 ) SCALE_BIOMASS_CO          : applies scale factors to CO for VOC 
!                                   oxidation
!  (3 ) INIT_BIOMASS              : Allocates & zeroes module arrays
!  (4 ) CLEANUP_BIOMASS           : Deallocates module arrays
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
!  (2 ) Modification for H2/HD simulation (phs, 9/18/07)
!  (3 ) Added 9 gaseous emissions from biomass burning: BENZ, TOLU, XYLE
!        C2H2, C2H4, GLYX, MGLY, GLYC, HAC  (tmf, 1/8/08)
!  (4 ) Hard-wired IDBCO2 and BIOTRCE (tmf, 7/30/08)
!  (5 ) Add CO scaling for VOC production. Routine SCALE_BIOMASS_CO 
!        transfered from gc_biomass_mod.f (jaf, mak, 2/6/09)
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
      INTEGER, PARAMETER   :: NBIOMAX     = 24
      INTEGER, PARAMETER   :: NBIOMAX_GAS = 19
      INTEGER, PARAMETER   :: IDBNOX      = 1
      INTEGER, PARAMETER   :: IDBCO       = 2
      INTEGER, PARAMETER   :: IDBC2H6     = 10
      INTEGER, PARAMETER   :: IDBSO2      = 11
      INTEGER, PARAMETER   :: IDBNH3      = 12
      INTEGER, PARAMETER   :: IDBBC       = 13
      INTEGER, PARAMETER   :: IDBOC       = 14
      INTEGER, PARAMETER   :: IDBCO2      = 24
      
      ! Arrays
      INTEGER              :: BIOTRCE(NBIOMAX)
      REAL*8,  ALLOCATABLE :: BIOMASS(:,:,:)

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
!  biomass burning emissions (Duncan et al 2001). (psk, bmy, 4/5/06, 9/18/07)
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
!  (2 ) Reference ITS_A_H2HD_SIM from "tracer_mod.f" to deal with ND29
!        (phs, 9/18/07)
!  (3 ) Now make a more general call to GFED2 reader to account for all
!        four options (phs, 17/12/08)
!  (4 ) Add CO scaling for VOC production (jaf, mak, 2/6/09)
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
      USE LOGICAL_MOD,       ONLY : L8DAYBB,  LSYNOPBB, L3HRBB
      USE TIME_MOD,          ONLY : ITS_A_NEW_MONTH
      USE TRACER_MOD,        ONLY : ITS_A_CO2_SIM
      USE TRACER_MOD,        ONLY : ITS_A_H2HD_SIM
      USE TRACER_MOD,        ONLY : ITS_A_FULLCHEM_SIM
      USE TRACER_MOD,        ONLY : ITS_A_TAGCO_SIM
      USE TRACERID_MOD,      ONLY : IDTBCPO, IDTNH3, IDTOCPO, IDTSO2
      USE TRACERID_MOD,      ONLY : IDTCO

#     include "CMN_SIZE"          ! Size parameters
#     include "CMN_DIAG"          ! Diagnostic flags

      ! Arguments
      INTEGER, INTENT(IN)        :: YEAR, MONTH

      ! Local variables
      LOGICAL, SAVE              :: FIRST = .TRUE.
      LOGICAL                    :: DO_ND28, DO_ND29, DO_ND32
      LOGICAL, SAVE              :: USE_GFED
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
            FIRST    = .FALSE.
            USE_GFED = LGFED2BB .or. L8DAYBB .or. LSYNOPBB .or. L3HRBB
         ENDIF
         
         ! Define diagnostic flags
         DO_ND28 = ( ND28 > 0 )
         DO_ND29 = ( ND29 > 0 )
         DO_ND32 = ( ND32 > 0 )

         !==============================================================
         ! GFED2 updates BIOMASS if needed (phs, 12/17/08)
         !==============================================================
         IF ( USE_GFED ) THEN

            ! Get emissions [molec/cm2/s] or [atoms C/cm2/s]
            CALL GFED2_COMPUTE_BIOMASS( YEAR, MONTH, BIOMASS )
               
            
         !==============================================================
         ! Read GC biomass emissions at the start of a new month
         !==============================================================
         ELSE IF ( ITS_A_NEW_MONTH() ) THEN

            ! Zero the array for biomass burning
            BIOMASS = 0d0

            ! Test if it's a CO2 simulation
            IF ( ITS_A_CO2_SIM() ) THEN

               !------------------------------
               ! CO2 emissions (based on 
               ! Duncan et al 2001 CO)
               !------------------------------
               
               ! Get CO2 emissions [molec/cm2/s]
               CALL GC_READ_BIOMASS_CO2( YEAR, MONTH,
     &                                   BIOMASS(:,:,IDBCO2) )
            ELSE

               !------------------------------
               ! Default GEOS-Chem inventory
               ! (Bryan Duncan et al 2001)
               !------------------------------

               ! Get emissions of gas-phase species
               ! in [molec/cm2/s] or [atoms C/cm2/s]
               CALL GC_COMPUTE_BIOMASS( YEAR, MONTH, 
     &                                  BIOMASS(:,:,1:NBIOMAX_GAS) )

               ! Get biomass SO2 [molec/cm2/s]
               IF ( IDTSO2 > 0 ) THEN
                  CALL GC_READ_BIOMASS_SO2( YEAR, MONTH, 
     &                                      BIOMASS(:,:,IDBSO2) )
               ENDIF

               ! Get biomass NH3 [molec/cm2/s]
               IF ( IDTNH3 > 0 ) THEN
                  CALL GC_READ_BIOMASS_NH3( YEAR, MONTH,
     &                                      BIOMASS(:,:,IDBNH3) )
               ENDIF

               ! Get biomass BC & OC [molec/cm2/s]
               IF ( IDTBCPO > 0 .and. IDTOCPO > 0 ) THEN
                  CALL GC_READ_BIOMASS_BCOC( YEAR, MONTH,
     &                                       BIOMASS(:,:,IDBBC), 
     &                                       BIOMASS(:,:,IDBOC) ) 
               ENDIF
            ENDIF
!            ENDIF
         ENDIF

         ! Irrespective of inventory type, we need to scale biomass
         ! CO to account for CO production from VOC's that are not
         ! explicitly carried in the chemistry mechanisms. This used
         ! to be done in gc_biomass_mod.f but then is not used for 
         ! GFED2, FLAMBE, etc. (jaf, mak, 2/6/09)
         IF ( ITS_A_FULLCHEM_SIM() ) THEN
            BIOMASS(:,:,IDBCO) = BIOMASS(:,:,IDBCO)*1.05d0
         ELSE IF ( ITS_A_TAGCO_SIM() ) THEN
            BIOMASS(:,:,IDBCO) = BIOMASS(:,:,IDBCO)*1.11d0
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

               IF ( ITS_A_H2HD_SIM() .and. (.not. USE_GFED ) ) THEN
                  AD29(I,J,2) = AD29(I,J,2) + 
     &                          BIOMASS(I,J,IDBCO) * 1.11d0
               ELSE
                  AD29(I,J,2) = AD29(I,J,2) + BIOMASS(I,J,IDBCO)
               ENDIF
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

      SUBROUTINE SCALE_BIOMASS_CO( BBARRAY )
!
!******************************************************************************
!  Subroutine SCALE_BIOMASS_CO multiplies the CO biomass emissions by scale 
!  factors to account for CO production from VOC's that are not explicitly 
!  carried in the chemistry mechanisms. (bnd, bmy, 8/21/01, 7/20/04)
!  
!  Arguments as Input:
!  ============================================================================
!  (1 ) BBARRAY (REAL*8) : Array containing biomass burning CO emissions
!
!  NOTES:
!  (1 ) Scale factors were determined by Jennifer Logan (jal@io.harvard.edu),
!       Bryan Duncan (bnd@io.harvard.edu) and Daniel Jacob (djj@io.harvard.edu)
!  (2 ) Scale factors have been corrected to 5% and 11% (bnd, bmy, 8/21/01)
!  (3 ) BBARRAY is now dimensioned (IIPAR,JJPAR) (bmy, 9/28/01)
!  (4 ) Removed obsolete code from 9/01 (bmy, 10/23/01)
!  (5 ) Now references ITS_A_FULLCHEM_SIM, ITS_A_TAGCO_SIM from "tracer_mod.f"
!        (bmy, 7/20/04)
!******************************************************************************
!
      ! References to F90 modules
      USE TRACER_MOD, ONLY : ITS_A_FULLCHEM_SIM, ITS_A_TAGCO_SIM

#     include "CMN_SIZE"    ! Size parameters

      ! Arguments
      REAL*8, INTENT(INOUT) :: BBARRAY(IIPAR,JJPAR) 

      !=================================================================
      ! SCALE_BIOMASS_CO begins here!
      !=================================================================
      IF ( ITS_A_FULLCHEM_SIM() ) THEN

         ! Full chemistry w/ SMVGEAR  -- enhance by 5%
         BBARRAY = BBARRAY * 1.05d0
         
      ELSE IF ( ITS_A_TAGCO_SIM() ) THEN

         ! Tagged CO -- enhance by 11%
         BBARRAY = BBARRAY * 1.11d0

      ENDIF

      ! Return to calling program  
      END SUBROUTINE SCALE_BIOMASS_CO

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
!  (3 ) Now set BIOTRCE for 24 biomass species (tmf, 7/30/08)
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
     &                   19, 20, 21, 26, 30, 34, 35, 
     &                   55, 56, 57, 58, 59, 63, 64, 
     &                   66, 67, 1/)
         ! Allocate array to hold monthly biomass emissions
         ALLOCATE( BIOMASS( IIPAR, JJPAR, NBIOMAX ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'BIOMASS' )
         BIOMASS = 0d0

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

      ! Return to calling program
      END SUBROUTINE CLEANUP_BIOMASS

!------------------------------------------------------------------------------
      
      ! End of module
      END MODULE BIOMASS_MOD
