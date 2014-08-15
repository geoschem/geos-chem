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
!  References:
!  ============================================================================
!
!  (1 ) Andreae, M.O., and P. Merlet, "Emissions of trace gases and aerosols
!        from biomass burning", Global Biogeochemical Cycles, Vol 15, pp
!        955-966, 2001.
!  (2 ) Hays, M.D., C.D. Geron, K.J. Linna, N.D. Smith, and J.J. Schauer, 
!        "Speciation of gas-phase and fine particle emissions from burning of
!        foliar fuels", Environ. Sci. Technol., Vol 36, pp 2281-2295, 2002.
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
!  14 Mar 2013 - M. Payer    - Replace NOx emissions with NO emissions as part
!                              of removal of NOx-Ox partitioning
!  13 Aug 2013 - M. Sulprizio- Add modifications for SOA + semivol POA (H. Pye)
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
!  20 Jun 2014 - R. Yantosca - Remove obsolete code now replaced by HEMCO
!******************************************************************************
!
      USE CMN_SIZE_MOD               ! Size parameters

      IMPLICIT NONE

      !=================================================================
      ! MODULE PRIVATE DECLARATIONS -- keep certain internal variables 
      ! and routines from being seen outside "biomass_mod.f"
      !=================================================================

      ! Make everything PRIVATE ...
      PRIVATE

      ! NBIOTRCE is the number of BB species in a simulation
      ! and is determined online FP (hotp 7/30/09)
      PUBLIC :: NBIOTRCE
      PUBLIC :: BIOMASS
      PUBLIC :: BIOTRCE

      ! Logical array to determine if species is treated in setemis 
      ! (hotp 8/3/09)
      PUBLIC :: BIOBGAS

      ! ... and these routines
      PUBLIC :: CLEANUP_BIOMASS
!      PUBLIC :: COMPUTE_BIOMASS_EMISSIONS
      !FP_ISOP
      ! FP SET_BIOTRCE is a new routine for setting the IDBxxx in BIOTRCE (hotp 7/30/09)
      PUBLIC :: SET_BIOTRCE

      !=================================================================
      ! MODULE VARIABLES
      !=================================================================
      
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

      SUBROUTINE SET_BIOTRCE( am_I_Root, Input_Opt )

!******************************************************************************
!  Routine created by FP to set IDBxxxs according to the tracers defined
!  for the run. The routine provides greater flexibility when ordering
!  tracers in input.geos (FP 6/2009, hotp 7/30/09)
!
!  NOTES:
!  30 Jul 2012 - R. Yantosca - Now accept am_I_Root as an argument when
!                              running with the traditional driver main.F
!  14 Mar 2013 - M. Payer    - Replace NOx emissions with NO emissions as part
!                              of removal of NOx-Ox partitioning
!  13 Aug 2013 - M. Sulprizio- Add NAP and POA biomass burning emissions for 
!                              SOA + semivolatile POA (H. Pye)
!  20 Feb 2014 - M. Sulprizio- Now pass Input_Opt argument
!  25 Jun 2014 - R. Yantosca - Removed references to tracer_mod.F
!******************************************************************************

!
! !USES:
!
      USE ERROR_MOD,          ONLY : ERROR_STOP
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
      USE TRACERID_MOD
!
! !INPUT PARAMETERS:
!
      LOGICAL,        INTENT(IN) :: am_I_Root   ! Is this the root CPU?
      TYPE(OptInput), INTENT(IN) :: Input_Opt   ! Input options
!
! !LOCAL VARIABLES:
!
      ! Local variables for quantities from Input_Opt
      LOGICAL :: LUCX

      !=================================================================
      ! SET_BIOTRCE begins here!
      !=================================================================

      ! Copy fields from INPUT_OPT
      LUCX = Input_Opt%LUCX

      ! Initialize
      NBIOTRCE = 0
      
      ! Increment NBIOTRCE for each turned on biomass burning tracer
      ! Replaced NOx with NO (mpayer, 3/14/13)
      IF ( IDBNO    /= 0 ) NBIOTRCE = NBIOTRCE + 1
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
      ! SOAupdate: Add NAP (hotp 8/23/09)
      IF ( IDBNAP   /= 0 ) NBIOTRCE = NBIOTRCE + 1

      ! Fill BIOTRCE w/ appropriate TRACER ID #'s
      ! Replaced NOx with NO (mpayer, 3/14/13)
      IF ( IDBNO    /= 0 ) BIOTRCE(IDBNO  ) = IDTNO
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
         IF ( IDTPOA1  /= 0 ) BIOTRCE(IDBOC) = IDTPOA1
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
      ! SOAupdate: Add NAP (hotp 8/23/09)
      IF ( IDBNAP   /= 0 ) BIOTRCE(IDBNAP)  = IDTNAP

      ! Set T/F of BIOBGAS() (hotp 8/3/09)
      ! BIOBGAS(IDBxx) is true if IDBxx should be handled by
      ! SMVGEAR/setemis
      ! Initialize to false
      BIOBGAS = .FALSE.
      ! Original gas phase emissions handled in setemis
      ! Replaced NOx with NO (mpayer, 3/14/13)
      IF ( IDBNO    /= 0 ) BIOBGAS(IDBNO  ) = .TRUE.
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
      ! SOAupdate: NAP handled by setemis (hotp 8/23/09)
      IF ( IDBNAP   /= 0 ) BIOBGAS(IDBNAP)  = .TRUE.

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
#if defined( UCX )
      IF (LUCX) THEN
         IF ( IDBCH4 /= 0 ) THEN
            NBIOTRCE          = NBIOTRCE + 1
            BIOTRCE( IDBCH4 ) = IDTCH4
            BIOBGAS( IDBCH4 ) = .TRUE. ! Handled in setemis
         ENDIF
      ELSEIF ( Input_Opt%ITS_A_CH4_SIM ) THEN
         IDBCH4            = NBIOTRCE + 1
         BIOTRCE( IDBCH4 ) = 1
         BIOBGAS( IDBCH4 ) = .FALSE.  ! not treated in setemis
      ENDIF
#else
      IF ( Input_Opt%ITS_A_CH4_SIM ) THEN
         IDBCH4            = NBIOTRCE + 1
         BIOTRCE( IDBCH4 ) = 1
         BIOBGAS( IDBCH4 ) = .FALSE.  ! not treated in setemis
      ENDIF
#endif

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
      IF ( am_I_Root ) THEN
         WRITE( 6, 100 ) BIOTRCE( 1:NBIOTRCE )
      ENDIF
 100  FORMAT( 'TRACERID: Biomass burning tracers        :', 20i3 )
      
      ! Return to calling program
      END SUBROUTINE SET_BIOTRCE

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
!      USE LOGICAL_MOD, ONLY : LBIOMASS

      USE CMN_SIZE_MOD    ! Size parameters

      INTEGER              :: AS

      !=================================================================
      ! INIT_BIOMASS begins here!
      !=================================================================

      ! If there are biomass emissions ...
!      IF ( LBIOMASS ) THEN

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

!      ENDIF
      
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
