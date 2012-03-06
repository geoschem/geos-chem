! $Id: biomass_mod.f,v 1.3 2010/03/15 19:33:25 ccarouge Exp $
      MODULE BIOMASS_MOD
!
!******************************************************************************
!  Module BIOMASS_MOD is a "wrapper" module, which allows us to select either
!  GFED2 biomass burning emissions, or the default GEOS-Chem biomass burning
!  emissions (based on Bryan Duncan et al).  (psk, bmy, 4/5/06, 11/6/09)
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
!  ----------------------------------------------------------------------------
!  FOR CH4 SIMULATION ONLY (kjw)
!
!   CH4        25         1          [molec CH4 /cm2/s]!
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
!  (8 ) gfed3_biomass_mod.f       : Module w/ routines for GFED3 biomass 
!  (9 ) grid_mod.f                : Module w/ horizontal grid information
!  (10 ) logical_mod.f             : Module w/ GEOS-CHEM logical switches
!  (11) time_mod.f                : Module w/ routines for computing time/ date
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
!  (6 ) Now always scale biomass CO regardless of inventory (jaf, mak, 11/6/09)
!  (7 ) Updates to remove all hard-wired order. (fp, 2/2/10)
!  (8)  Updates to include GFED3 (psk, 1/5/11)
!  14 Feb 2012 - M. Payer    - Add modifications for CH4 simulation (K. Wecht)
!******************************************************************************
!
      USE CMN_SIZE_MOD               ! Size parameters

      IMPLICIT NONE
#     include "define.h"

      !=================================================================
      ! MODULE PRIVATE DECLARATIONS -- keep certain internal variables 
      ! and routines from being seen outside "biomass_mod.f"
      !=================================================================

      ! Make everything PRIVATE ...
      PRIVATE

      ! ... except these variables
      ! FP: NBIOMAX now specified in CMN_SIZE_mod and is the max
      ! number of BB species allowed, used for array allocation (hotp 7/30/09)
      ! if additional BB species are added, increase NBIOMAX in CMN_SIZE_mod
      !PUBLIC :: NBIOMAX
      !PUBLIC :: NBIOMAX_GAS ! no longer needed (hotp 8/3/09)
      !FP_ISOP
      ! NBIOTRCE is the number of BB species in a simulation
      ! and is determined online FP (hotp 7/30/09)
      PUBLIC :: NBIOTRCE
      PUBLIC :: BIOMASS
      PUBLIC :: BIOTRCE
      ! IDBxxx are now determined online and stored in
      ! BIOTRCE array by FP (hotp 7/30/09)
      ! IDBxxx are part of tracerid_mod by FP (hotp 7/31/09)
      !PUBLIC :: IDBBC
      !PUBLIC :: IDBCO
      !PUBLIC :: IDBCO2
      !PUBLIC :: IDBC2H6
      !PUBLIC :: IDBNH3
      !PUBLIC :: IDBNOX
      !PUBLIC :: IDBOC
      !PUBLIC :: IDBSO2

      ! Logical array to determine if species is treated in setemis 
      ! (hotp 8/3/09)
      PUBLIC :: BIOBGAS

      ! ... and these routines
      PUBLIC :: CLEANUP_BIOMASS
      PUBLIC :: COMPUTE_BIOMASS_EMISSIONS
      !FP_ISOP
      ! FP SET_BIOTRCE is a new routine for setting the IDBxxx in BIOTRCE (hotp 7/30/09)
      PUBLIC :: SET_BIOTRCE

      !=================================================================
      ! MODULE VARIABLES
      !=================================================================

      ! FP Removed hard-wired IDBxxx and move NBIOMAX to CMN_SIZE_mod (hotp 7/30/09)      ! Parameters
!      INTEGER, PARAMETER   :: NBIOMAX     = 24
!      INTEGER, PARAMETER   :: NBIOMAX_GAS = 19
!      INTEGER, PARAMETER   :: IDBNOX      = 1
!      INTEGER, PARAMETER   :: IDBCO       = 2
!      INTEGER, PARAMETER   :: IDBC2H6     = 10
!      INTEGER, PARAMETER   :: IDBSO2      = 11
!      INTEGER, PARAMETER   :: IDBNH3      = 12
!      INTEGER, PARAMETER   :: IDBBC       = 13
!      INTEGER, PARAMETER   :: IDBOC       = 14
!      INTEGER, PARAMETER   :: IDBCO2      = 24
      
      !FP_ISOP
      ! FP: NBIOTRCE is number of biomass burning tracers (hotp 7/30/09)
      INTEGER              :: NBIOTRCE

      ! Arrays
      INTEGER              :: BIOTRCE(NBIOMAX)
      REAL*8,  ALLOCATABLE :: BIOMASS(:,:,:)

      ! Logical array to determine if an emission should be treated
      ! in setemis or elsewhere (hotp 8/3/09)
      ! make sure it saves it's value
      LOGICAL, SAVE        :: BIOBGAS(NBIOMAX)

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
!  biomass burning emissions (Duncan et al 2001). (psk, bmy, 4/5/06, 11/6/09)
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
!  (5 ) Irrespective of inventory type, we need to scale biomass CO to account 
!        for CO production from VOC's that are not explicitly carried in the 
!        chemistry mechanisms. This used to be done in gc_biomass_mod.f but 
!        then is not used for GFED2, FLAMBE, etc. (jaf, mak, 11/6/09)
!  27 Jul 2011 - R. Yantosca - Make sure that IDBCO > 0 before using it
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
      USE GFED3_BIOMASS_MOD, ONLY : GFED3_COMPUTE_BIOMASS
      USE LOGICAL_MOD,       ONLY : LBIOMASS, LGFED2BB, LGFED3BB
      USE LOGICAL_MOD,       ONLY : L8DAYBB,  LSYNOPBB, L3HRBB
      USE LOGICAL_MOD,       ONLY : L8DAYBB3,  LSYNOPBB3, L3HRBB3
      USE TIME_MOD,          ONLY : ITS_A_NEW_MONTH
      USE TRACER_MOD,        ONLY : ITS_A_CO2_SIM
      USE TRACER_MOD,        ONLY : ITS_A_H2HD_SIM
      USE TRACER_MOD,        ONLY : ITS_A_FULLCHEM_SIM
      USE TRACER_MOD,        ONLY : ITS_A_TAGCO_SIM
      USE TRACERID_MOD,      ONLY : IDTBCPO, IDTNH3, IDTOCPO, IDTSO2
      !FP_ISOP
      ! FP: IDBxxxs will need to used here (hotp 7/30/09)
      USE TRACERID_MOD,      ONLY : IDBBC, IDBNH3, IDBOC, IDBSO2
      USE TRACERID_MOD,      ONLY : IDBCO, IDBNOx
      USE TRACERID_MOD,      ONLY : IDBCO2


      USE CMN_SIZE_MOD          ! Size parameters
      USE CMN_DIAG_MOD          ! Diagnostic flags

      ! Arguments
      INTEGER, INTENT(IN)        :: YEAR, MONTH

      ! Local variables
      LOGICAL, SAVE              :: FIRST = .TRUE.
      LOGICAL                    :: DO_ND28, DO_ND29, DO_ND32
      LOGICAL, SAVE              :: USE_GFED
      LOGICAL, SAVE              :: USE_GFED2
      LOGICAL, SAVE              :: USE_GFED3
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
            USE_GFED2 = LGFED2BB .or. L8DAYBB .or. LSYNOPBB .or. L3HRBB
            USE_GFED3 = LGFED3BB .or.
     &                  L8DAYBB3 .or. LSYNOPBB3 .or. L3HRBB3
            USE_GFED = USE_GFED2. or. USE_GFED3
         ENDIF
         
         ! Define diagnostic flags
         DO_ND28 = ( ND28 > 0 )
         DO_ND29 = ( ND29 > 0 )
         DO_ND32 = ( ND32 > 0 )

         !==============================================================
         ! GFED2 updates BIOMASS if needed (phs, 12/17/08)
         !==============================================================
         IF ( USE_GFED2 ) THEN

            ! Get emissions [molec/cm2/s] or [atoms C/cm2/s]
            CALL GFED2_COMPUTE_BIOMASS( YEAR, MONTH, BIOMASS )
               
            ! Irrespective of inventory type, we need to scale biomass
            ! CO to account for CO production from VOC's that are not
            ! explicitly carried in the chemistry mechanisms. This used
            ! to be done in gc_biomass_mod.f but then is not used for 
            ! GFED2, FLAMBE, etc. (jaf, mak, 2/6/09)
            !
            ! Now make sure IDBCO is defined before using it.  This will
            ! avoid an out-of-bounds error. (bmy, 7/29/11)
            IF ( IDBCO > 0 ) THEN
               IF ( ITS_A_FULLCHEM_SIM() ) THEN
                  BIOMASS(:,:,IDBCO) = BIOMASS(:,:,IDBCO) * 1.05d0
               ELSE IF ( ITS_A_TAGCO_SIM() ) THEN
                  BIOMASS(:,:,IDBCO) = BIOMASS(:,:,IDBCO) * 1.11d0
               ENDIF
            ENDIF

         !==============================================================
         ! GFED3 updates BIOMASS if needed (psk, 1/5/11)
         !==============================================================
         ELSE IF ( USE_GFED3 ) THEN

            ! Get emissions [molec/cm2/s] or [atoms C/cm2/s]
            CALL GFED3_COMPUTE_BIOMASS( YEAR, MONTH, BIOMASS )

            ! Irrespective of inventory type, we need to scale biomass
            ! CO to account for CO production from VOC's that are not
            ! explicitly carried in the chemistry mechanisms. This used
            ! to be done in gc_biomass_mod.f but then is not used for
            ! GFED3, FLAMBE, etc. (jaf, mak, 2/6/09)
            IF ( ITS_A_FULLCHEM_SIM() ) THEN
               BIOMASS(:,:,IDBCO) = BIOMASS(:,:,IDBCO) * 1.05d0
            ELSE IF ( ITS_A_TAGCO_SIM() ) THEN
               BIOMASS(:,:,IDBCO) = BIOMASS(:,:,IDBCO) * 1.11d0
            ENDIF

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
!               CALL GC_COMPUTE_BIOMASS( YEAR, MONTH, 
!     &                                  BIOMASS(:,:,1:NBIOMAX_GAS) )
               ! Now send entire BIOMASS array since gas and aerosol
               ! species might be intermixed (hotp 8/3/09)
               ! Only the standard 10 gas phase species will be updated
               CALL GC_COMPUTE_BIOMASS( YEAR, MONTH, BIOMASS )

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

            ! Irrespective of inventory type, we need to scale biomass
            ! CO to account for CO production from VOC's that are not
            ! explicitly carried in the chemistry mechanisms. This used
            ! to be done in gc_biomass_mod.f but then is not used for 
            ! GFED2, FLAMBE, etc. (jaf, mak, 2/6/09)
            !
            ! Now make sure IDBCO is defined before using it.  This will
            ! avoid an out-of-bounds error. (bmy, 7/29/11)
            IF ( IDBCO > 0 ) THEN
               IF ( ITS_A_FULLCHEM_SIM() ) THEN
                  BIOMASS(:,:,IDBCO) = BIOMASS(:,:,IDBCO) * 1.05d0
               ELSE IF ( ITS_A_TAGCO_SIM() ) THEN
                  BIOMASS(:,:,IDBCO) = BIOMASS(:,:,IDBCO) * 1.11d0
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
      SUBROUTINE SET_BIOTRCE

!******************************************************************************
!  Routine created by FP to set IDBxxxs according to the tracers defined
!  for the run. The routine provides greater flexibility when ordering
!  tracers in input.geos (FP 6/2009, hotp 7/30/09)
!
!  NOTES:
!******************************************************************************

!FP_ISOP
!AS BF_TRACE

      ! References to F90 modules
      ! ERROR_STOP to prevent array alloc issues (hotp 7/31/09)
      USE ERROR_MOD, ONLY    : ERROR_STOP

      ! need both IDBxxss and IDTxxxs
      USE TRACERID_MOD, ONLY : IDBNOX,  IDBCO,   IDBALK4, IDBACET 
      USE TRACERID_MOD, ONLY : IDBMEK,  IDBALD2, IDBPRPE, IDBC3H8
      USE TRACERID_MOD, ONLY : IDBCH2O, IDBC2H6
      USE TRACERID_MOD, ONLY : IDBSO2,  IDBNH3
      USE TRACERID_MOD, ONLY : IDBBC,   IDBOC
      USE TRACERID_MOD, ONLY : IDBCO2
      USE TRACERID_MOD, ONLY : IDBCH4   !kjw
      USE TRACERID_MOD, ONLY : IDBBENZ, IDBTOLU, IDBXYLE
      ! Add dicarbonyls
      USE TRACERID_MOD, ONLY : IDBGLYX, IDBMGLY, IDBC2H4
      USE TRACERID_MOD, ONLY : IDBC2H2, IDBGLYC, IDBHAC

      USE TRACERID_MOD, ONLY : IDTNOX,  IDTCO,   IDTALK4, IDTACET 
      USE TRACERID_MOD, ONLY : IDTMEK,  IDTALD2, IDTPRPE, IDTC3H8
      USE TRACERID_MOD, ONLY : IDTCH2O, IDTC2H6
      USE TRACERID_MOD, ONLY : IDTSO2,  IDTNH3
      USE TRACERID_MOD, ONLY : IDTBCPI, IDTOCPI
      USE TRACERID_MOD, ONLY : IDTBENZ, IDTTOLU, IDTXYLE
      USE TRACERID_MOD, ONLY : IDTHAC, IDTGLYC, IDTMGLY, IDTGLYX
      USE TRACERID_MOD, ONLY : IDTC2H2, IDTC2H4

      USE TRACER_MOD,   ONLY : ITS_A_CH4_SIM ! (mpayer, 2/14/12)

      !=================================================================
      ! SET_BIOTRCE begins here!
      !=================================================================

      ! Initialize
      NBIOTRCE = 0
      
      ! Increment NBIOTRCE for each turned on biomass burning tracer
      IF ( IDBNOX   /= 0 ) NBIOTRCE = NBIOTRCE + 1
      IF ( IDBCO    /= 0 ) NBIOTRCE = NBIOTRCE + 1 
      IF ( IDBALK4  /= 0 ) NBIOTRCE = NBIOTRCE + 1 
      IF ( IDBACET  /= 0 ) NBIOTRCE = NBIOTRCE + 1 
      IF ( IDBMEK   /= 0 ) NBIOTRCE = NBIOTRCE + 1 
      IF ( IDBALD2  /= 0 ) NBIOTRCE = NBIOTRCE + 1 
      IF ( IDBPRPE  /= 0 ) NBIOTRCE = NBIOTRCE + 1 
      IF ( IDBC3H8  /= 0 ) NBIOTRCE = NBIOTRCE + 1 
      IF ( IDBCH2O  /= 0 ) NBIOTRCE = NBIOTRCE + 1 
      IF ( IDBC2H6  /= 0 ) NBIOTRCE = NBIOTRCE + 1 
      !
      IF ( IDBSO2   /= 0 ) NBIOTRCE = NBIOTRCE + 1 
      IF ( IDBNH3   /= 0 ) NBIOTRCE = NBIOTRCE + 1 
      IF ( IDBBC    /= 0 ) NBIOTRCE = NBIOTRCE + 1 
      IF ( IDBOC    /= 0 ) NBIOTRCE = NBIOTRCE + 1 
      !

      IF ( IDBXYLE  /= 0 ) NBIOTRCE = NBIOTRCE + 1 
      IF ( IDBBENZ  /= 0 ) NBIOTRCE = NBIOTRCE + 1 
      IF ( IDBTOLU  /= 0 ) NBIOTRCE = NBIOTRCE + 1

      IF ( IDBGLYX  /= 0 ) NBIOTRCE = NBIOTRCE + 1 
      IF ( IDBMGLY  /= 0 ) NBIOTRCE = NBIOTRCE + 1 
      IF ( IDBC2H4  /= 0 ) NBIOTRCE = NBIOTRCE + 1 
      IF ( IDBC2H2  /= 0 ) NBIOTRCE = NBIOTRCE + 1 
      IF ( IDBGLYC  /= 0 ) NBIOTRCE = NBIOTRCE + 1 
      IF ( IDBHAC   /= 0 ) NBIOTRCE = NBIOTRCE + 1 

      ! Fill BIOTRCE w/ appropriate TRACER ID #'s
      IF ( IDBNOX   /= 0 ) BIOTRCE(IDBNOX ) = IDTNOX
      IF ( IDBCO    /= 0 ) BIOTRCE(IDBCO  ) = IDTCO
      IF ( IDBALK4  /= 0 ) BIOTRCE(IDBALK4) = IDTALK4
      IF ( IDBACET  /= 0 ) BIOTRCE(IDBACET) = IDTACET
      IF ( IDBMEK   /= 0 ) BIOTRCE(IDBMEK ) = IDTMEK
      IF ( IDBALD2  /= 0 ) BIOTRCE(IDBALD2) = IDTALD2
      IF ( IDBPRPE  /= 0 ) BIOTRCE(IDBPRPE) = IDTPRPE
      IF ( IDBC3H8  /= 0 ) BIOTRCE(IDBC3H8) = IDTC3H8
      IF ( IDBCH2O  /= 0 ) BIOTRCE(IDBCH2O) = IDTCH2O  
      IF ( IDBC2H6  /= 0 ) BIOTRCE(IDBC2H6) = IDTC2H6 
      !
      IF ( IDBSO2   /= 0 ) BIOTRCE(IDBSO2)  = IDTSO2 
      IF ( IDBNH3   /= 0 ) BIOTRCE(IDBNH3)  = IDTNH3 
      IF ( IDBBC    /= 0 ) BIOTRCE(IDBBC)   = IDTBCPI

      ! IDBOC may be associated with either OCPI or POA 
      ! depending on if POA is semivolatile (hotp 8/23/09)
      IF ( IDBOC    /= 0 ) THEN 
          IF ( IDTOCPI  /= 0 ) BIOTRCE(IDBOC) = IDTOCPI
      ENDIF

      IF ( IDBXYLE  /= 0 ) BIOTRCE(IDBXYLE) = IDTXYLE 
      IF ( IDBTOLU  /= 0 ) BIOTRCE(IDBTOLU) = IDTTOLU
      IF ( IDBBENZ  /= 0 ) BIOTRCE(IDBBENZ) = IDTBENZ

      ! Dicarbonyls 
      IF ( IDBGLYX  /= 0 ) BIOTRCE(IDBGLYX) = IDTGLYX
      IF ( IDBMGLY  /= 0 ) BIOTRCE(IDBMGLY) = IDTMGLY
      IF ( IDBC2H4  /= 0 ) BIOTRCE(IDBC2H4) = IDTC2H4
      IF ( IDBC2H2  /= 0 ) BIOTRCE(IDBC2H2) = IDTC2H2
      IF ( IDBGLYC  /= 0 ) BIOTRCE(IDBGLYC) = IDTGLYC
      IF ( IDBHAC   /= 0 ) BIOTRCE(IDBHAC)  = IDTHAC


      ! Set T/F of BIOBGAS() (hotp 8/3/09)
      ! BIOBGAS(IDBxx) is true if IDBxx should be handled by
      ! SMVGEAR/setemis
      ! Initialize to false
      BIOBGAS = .FALSE.
      ! Original gas phase emissions handled in setemis
      IF ( IDBNOX   /= 0 ) BIOBGAS(IDBNOX ) = .TRUE.
      IF ( IDBCO    /= 0 ) BIOBGAS(IDBCO  ) = .TRUE.
      IF ( IDBALK4  /= 0 ) BIOBGAS(IDBALK4) = .TRUE.
      IF ( IDBACET  /= 0 ) BIOBGAS(IDBACET) = .TRUE.
      IF ( IDBMEK   /= 0 ) BIOBGAS(IDBMEK ) = .TRUE.
      IF ( IDBALD2  /= 0 ) BIOBGAS(IDBALD2) = .TRUE.
      IF ( IDBPRPE  /= 0 ) BIOBGAS(IDBPRPE) = .TRUE.
      IF ( IDBC3H8  /= 0 ) BIOBGAS(IDBC3H8) = .TRUE.
      IF ( IDBCH2O  /= 0 ) BIOBGAS(IDBCH2O) = .TRUE.  
      IF ( IDBC2H6  /= 0 ) BIOBGAS(IDBC2H6) = .TRUE. 
      ! SO2 and NH3 handled by sulfate_mod.f
      IF ( IDBSO2   /= 0 ) BIOBGAS(IDBSO2)  = .FALSE. 
      IF ( IDBNH3   /= 0 ) BIOBGAS(IDBNH3)  = .FALSE. 
      ! BC and OC handled by carbon_mod.f
      IF ( IDBBC    /= 0 ) BIOBGAS(IDBBC)   = .FALSE.
      IF ( IDBOC    /= 0 ) BIOBGAS(IDBOC)   = .FALSE.
      ! Aromatics handled by setemis
      IF ( IDBXYLE  /= 0 ) BIOBGAS(IDBXYLE) = .TRUE. 
      IF ( IDBTOLU  /= 0 ) BIOBGAS(IDBTOLU) = .TRUE.
      IF ( IDBBENZ  /= 0 ) BIOBGAS(IDBBENZ) = .TRUE.

      ! Dicarbonyls handled by setemis
      IF ( IDBGLYX  /= 0 ) BIOBGAS(IDBGLYX) = .TRUE.
      IF ( IDBMGLY  /= 0 ) BIOBGAS(IDBMGLY) = .TRUE.
      IF ( IDBC2H4  /= 0 ) BIOBGAS(IDBC2H4) = .TRUE.
      IF ( IDBC2H2  /= 0 ) BIOBGAS(IDBC2H2) = .TRUE.
      IF ( IDBGLYC  /= 0 ) BIOBGAS(IDBGLYC) = .TRUE.
      IF ( IDBHAC   /= 0 ) BIOBGAS(IDBHAC ) = .TRUE.



      !FP_ISOP
      !FOR CO2
      ! IDBCO2 is placed at the end after all full chem tracers
      IDBCO2            = NBIOTRCE + 1
      BIOTRCE( IDBCO2 ) = 1
      BIOBGAS( IDBCO2 ) = .FALSE.

      ! For CH4 (kjw)
      !  IDBCH4 is placed at end after all full chem tracers and CO2
      IF ( ITS_A_CH4_SIM() ) THEN
         IDBCH4            = NBIOTRCE + 1
         BIOTRCE( IDBCH4 ) = 1
         BIOBGAS( IDBCH4 ) = .FALSE.  ! not treated in setemis
      ENDIF

      ! Error check 
      ! NBIOTRCE should be less than or equal to NBIOMAX (hotp 7/31/09)
      !IF ( NBIOTRCE > NBIOMAX ) THEN
      ! Need to account for the fact that CO2 is placed last in full
      ! chem sim (hotp 8/3/09)
      ! if NBIOTRCE and NBIOMAX are equal, BIOTRCE( IDBCO2 ) will have
      ! indexed out of array dimensions
      IF ( NBIOTRCE .GE. NBIOMAX ) THEN
          CALL ERROR_STOP ( 'NBIOTRCE too large!', '(biomass_mod)' )
          ! If you get this error, an appropriate fix is to increase
          ! NBIOMAX in CMN_SIZE_mod (hotp 7/31/09)
      ENDIF

      ! Echo biomass burning tracer information
      WRITE( 6, 100 ) BIOTRCE( 1:NBIOTRCE )
 100  FORMAT( 'TRACERID: Biomass burning tracers        :', 20i3 )
      
      ! Return to calling program
      END SUBROUTINE SET_BIOTRCE

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

      USE CMN_SIZE_MOD    ! Size parameters

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

      USE CMN_SIZE_MOD    ! Size parameters

      INTEGER              :: AS

      !=================================================================
      ! INIT_BIOMASS begins here!
      !=================================================================

      ! If there are biomass emissions ...
      IF ( LBIOMASS ) THEN

         ! Tracer numbers for each biomass species (CO2 is last)
         ! BIOTRCE is now computed online in SET_BIOTRCE. (fp, 2/2/10)
!         BIOTRCE(:) = (/ 1,  4,  5,  9,  10, 11, 18, 
!     &                   19, 20, 21, 26, 30, 34, 35, 
!     &                   55, 56, 57, 58, 59, 63, 64, 
!     &                   66, 67, 1/)
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
