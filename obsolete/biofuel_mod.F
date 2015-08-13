! $Id: biofuel_mod.f,v 1.1 2009/09/16 14:06:40 bmy Exp $
      MODULE BIOFUEL_MOD
!
!******************************************************************************
!  Module BIOFUEL_MOD contains arrays and routines to compute yearly
!  biofuel emissions for NOx, CO, ALK4, ACET, MEK, ALD2, PRPE, C3H8, 
!  CH2O, and C2H6. (bmy, 9/12/00, 9/18/07)
!
!  Module Variables:
!  ============================================================================
!  (1 ) NBFMAX             : Maximum # of biofuel burning species
!  (2 ) NBFTRACE           : # of emitted biofuel burning species (<= NBFMAX)
!  (3 ) BFTRACE            : Array of tracer #'s for emitted biofuel species
!  (4 ) BIOFUEL            : array containing biofuel emissions
!
!  Module Routines:
!  ============================================================================
!  (1 ) BIOFUEL_BURN       : Reads data from disk & computes biofuel emissions
!  (2 ) SCALE_BIOFUEL_CO   : Scales biofuel CO to account for VOC oxidation
!  (3 ) SCALE_BIOFUEL_ACET : Scales biofuel ACET to match a posteriori source
!  (4 ) SCALE_FUTURE       : Applies future scale factors to biofuel emissions
!  (5 ) SET_BFTRACE        : Initializes NBFTRACE counter and BFTRACE array
!  (6 ) INIT_BIOFUEL       : Initializes the BIOFUEL array
!  (7 ) CLEANUP_BIOFUEL    : Deallocates the BIOFUEL array
!
!  GEOS-CHEM modules referenced by biofuel_mod.f
!  ============================================================================
!  (1 ) bpch2_mod.f        : Module w/ routines for binary punch file I/O
!  (2 ) dao_mod.f          : Module w/ DAO met field arrays
!  (3 ) diag_mod.f         : Module w/ GEOS-CHEM diagnostic arrays
!  (4 ) directory_mod.f    : Module w/ GEOS-CHEM data & met field dirs
!  (5 ) epa_nei_mod.f      : Module w/ routines to read EPA/NEI99 data
!  (6 ) error_mod.f        : Module w/ NaN and other error check routines
!  (7 ) logical_mod.f      : Module w/ GEOS-CHEM logical switches
!  (8 ) tracer_mod.f       : Module w/ GEOS-CHEM tracer array etc.
!  (9 ) tracerid_mod.f     : Module w/ pointers to tracers & emissions 
!  (10) transfer_mod.f     : Module w/ routines to cast & resize arrays
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
!  (1 ) Now account for extra production of CO from VOC's for Tagged CO
!        and CO-OH simulations (bmy, 1/3/01)
!  (2 ) Now read NBIOFUEL=10 biofuel species.  Also archive biofuel emissions
!        in the ND34 diagnostic. (bmy, 4/17/01)
!  (3 ) Now dimension BFTRACE arrays to be of size NBFMAX instead of having 
!        them be made allocatable.  Also updated comments. (bmy, 4/17/01)
!  (4 ) Bug fix: now make sure to index biofuel tracers w/ the correct tracer
!        number, even when there are less than the maximum species being
!        requested (bmy, 8/24/01)
!  (5 ) Bug fix: now index biofuel CH2O correctly (bmy, 8/28/01)
!  (6 ) Now scale biofuel ACET by 0.82, in order to match the a posteriori
!        acetone source from Jacob et al 2001.  Also updated comments.
!        (bdf, bmy, 9/10/01)
!  (7 ) BIOFUEL is now declared (NBFTRACE,IIPAR,JJPAR).  Now use TRANSFER_2D
!        from "transfer_mod.f" to copy data into BIOFUEL. (bmy, 9/28/01)
!  (8 ) Deleted obsolete code from 9/01 (bmy, 11/15/01)
!  (9 ) Now do unit conversion every time step.  Also added private
!        array BIOFUEL_KG to hold emissions in kg over the entire
!        month. (bmy, 5/9/02)
!  (10) Now divide module header into MODULE PRIVATE, MODULE VARIABLES, and
!        MODULE ROUTINES sections.  Updated comments.  BIOMASS_KG is now
!        an allocatable module array instead of a local array in routine
!        "biofuel_burn.f". (bmy, 5/28/02)  
!  (11) Now reference BXHEIGHT from "dao_mod.f". Now references "error_mod.f".
!        Also deleted obsolete code from various routines.  Also references
!        "tracerid_mod.f"  Added routine SET_NBFTRACE. (bmy, 11/6/02)
!  (12) Now call READ_BPCH2 with QUIET=.TRUE. to suppress output (bmy, 3/14/03)
!  (13) Now references "directory_mod.f" (bmy, 7/19/04)
!  (14) Now references "time_mod.f" and "epa_nei_mod.f" (bmy, 11/5/04)
!  (15) Now can read data for both GEOS and GCAP grids (bmy, 8/16/05)
!  (16) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (17) Rewrite IF statements to avoid seg fault errors when LNEI99 is turned 
!        off (bmy, 2/1/06)
!  (18) Modified for IPCC future emissions scale factors.  Added private
!        routine SCALE_FUTURE. (swu, bmy, 5/30/06)
!  (19) Modified for VOC-scaling of CO emissions for H2/HD sim (phs, 5/16/07)
!  (20) Added 9 gaseous biofuel emissions: GLYX, MGLY, BENZ, 
!        TOLU, XYLE, C2H4, C2H2, GLYC, HAC. (tmf, 1/7/09)
!  (21) Emissions for these 9 tracers are scaled from CO emissions. (tmf, 1/7/09)
!  14 Mar 2013 - M. Payer    - Replace NOx emissions with NO emissions as part
!                              of removal of NOx-Ox partitioning
!  12 Jun 2013 - R. Yantosca - Now declare BIOFUEL_KG w/ the TARGET attribute
!  14 Jun 2013 - R. Yantosca - Now determine weekday/weekend with respect to
!                              the local time at each grid box.  (Formerly,
!                              this had been done w/r/t the GMT time.
!  13 Aug 2013 - M. Sulprizio- Add NAP biofuel emissions for SOA + semivolatile
!                              POA simulation (H. Pye)
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
!  20 Jun 2014 - R. Yantosca - Removed BIOFUEL_BURN etc; HEMCO replaces this
!******************************************************************************
!
      IMPLICIT NONE

      !=================================================================
      ! MODULE PRIVATE DECLARATIONS -- keep certain internal variables 
      ! and routines from being seen outside "biofuel_mod.f"
      !=================================================================

      ! Make everything PRIVATE ...
      PRIVATE

      ! ... except these variables ...
      PUBLIC :: NBFMAX
      PUBLIC :: NBFTRACE
      PUBLIC :: BFTRACE
      PUBLIC :: BIOFUEL

      ! ... and these routines
!      PUBLIC :: BIOFUEL_BURN
      PUBLIC :: CLEANUP_BIOFUEL
      PUBLIC :: INIT_BIOFUEL
      PUBLIC :: SET_BFTRACE

      !=================================================================     
      ! MODULE VARIABLES
      !=================================================================

      ! SOAupdate: increase NBFMAX for NAP (hotp 6/25/09)
      INTEGER, PARAMETER          :: NBFMAX = 20

      INTEGER                     :: NBFTRACE
      INTEGER                     :: BFTRACE(NBFMAX) 

      REAL*8, ALLOCATABLE         :: BIOFUEL(:,:,:)
      REAL*8, ALLOCATABLE, TARGET :: BIOFUEL_KG(:,:,:)

      !=================================================================
      ! MODULE ROUTINES -- follow below the "CONTAINS" statement 
      !=================================================================
      CONTAINS
!EO
      SUBROUTINE SET_BFTRACE( am_I_Root )
!
!******************************************************************************
!  Subroutine SET_NBFTRACE sets the NBFTRACE variable with the number of
!  biofuel tracers that are turned on.  This was split off from "tracerid.f"
!  in order to prevent circular module references. (bmy, 11/6/02, 10/3/05)
!
!  NOTES:
!  (1 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05
!  30 Jul 2012 - R. Yantosca - Now accept am_I_Root as an argument when
!                              running with the traditional driver main.F
!  14 Mar 2013 - M. Payer    - Replace NOx emissions with NO emissions as part
!                              of removal of NOx-Ox partitioning
!  13 Aug 2013 - M. Sulprizio- Add NAP for SOA + semivolatile POA (H. Pye)
!******************************************************************************
!
      ! References to F90 modules
      USE TRACERID_MOD, ONLY : IDBFACET, IDBFALD2, IDBFALK4, IDBFC2H6 
      USE TRACERID_MOD, ONLY : IDBFC3H8, IDBFCH2O, IDBFCO,   IDBFMEK  
      USE TRACERID_MOD, ONLY : IDBFNO,   IDBFPRPE, IDTACET,  IDTALD2 
      USE TRACERID_MOD, ONLY : IDTALK4,  IDTC2H6,  IDTC3H8,  IDTCH2O
      USE TRACERID_MOD, ONLY : IDTCO,    IDTMEK,   IDTNO,    IDTPRPE 

      USE TRACERID_MOD, ONLY : IDBFGLYX, IDBFMGLY, IDBFBENZ, IDBFTOLU
      USE TRACERID_MOD, ONLY : IDBFXYLE, IDBFC2H4, IDBFC2H2, IDBFGLYC
      USE TRACERID_MOD, ONLY : IDBFHAC
      USE TRACERID_MOD, ONLY : IDTGLYX,  IDTMGLY,  IDTBENZ,  IDTTOLU
      USE TRACERID_MOD, ONLY : IDTXYLE,  IDTC2H4,  IDTC2H2,  IDTGLYC
      USE TRACERID_MOD, ONLY : IDTHAC
      ! SOAupdate: for gas phase NAP chemistry, NAP biofuel emis (hotp 6/24/09)
      USE TRACERID_MOD, ONLY : IDTNAP
      USE TRACERID_MOD, ONLY : IDBFNAP

      LOGICAL, INTENT(IN) :: am_I_Root

      !=================================================================
      ! SET_BFTRACE begins here!
      !=================================================================

      ! Initialize
      NBFTRACE = 0
      
      ! Increment NBFTRACE for each turned on biofuel tracer
      ! Replaced NOx with NO (mpayer, 3/14/13)
      IF ( IDBFNO   /= 0 ) NBFTRACE = NBFTRACE + 1
      IF ( IDBFCO   /= 0 ) NBFTRACE = NBFTRACE + 1 
      IF ( IDBFALK4 /= 0 ) NBFTRACE = NBFTRACE + 1 
      IF ( IDBFACET /= 0 ) NBFTRACE = NBFTRACE + 1 
      IF ( IDBFMEK  /= 0 ) NBFTRACE = NBFTRACE + 1 
      IF ( IDBFALD2 /= 0 ) NBFTRACE = NBFTRACE + 1 
      IF ( IDBFPRPE /= 0 ) NBFTRACE = NBFTRACE + 1 
      IF ( IDBFC3H8 /= 0 ) NBFTRACE = NBFTRACE + 1 
      IF ( IDBFCH2O /= 0 ) NBFTRACE = NBFTRACE + 1 
      IF ( IDBFC2H6 /= 0 ) NBFTRACE = NBFTRACE + 1 
      IF ( IDBFGLYX /= 0 ) NBFTRACE = NBFTRACE + 1 
      IF ( IDBFMGLY /= 0 ) NBFTRACE = NBFTRACE + 1 
      IF ( IDBFBENZ /= 0 ) NBFTRACE = NBFTRACE + 1 
      IF ( IDBFTOLU /= 0 ) NBFTRACE = NBFTRACE + 1 
      IF ( IDBFXYLE /= 0 ) NBFTRACE = NBFTRACE + 1 
      IF ( IDBFC2H4 /= 0 ) NBFTRACE = NBFTRACE + 1 
      IF ( IDBFC2H2 /= 0 ) NBFTRACE = NBFTRACE + 1 
      IF ( IDBFGLYC /= 0 ) NBFTRACE = NBFTRACE + 1 
      IF ( IDBFHAC  /= 0 ) NBFTRACE = NBFTRACE + 1 
      ! SOAupdate: for gas phase NAP chemistry, NAP bf emiss (hotp 6/24/09)
      IF ( IDBFNAP  /= 0 ) NBFTRACE = NBFTRACE + 1

      ! Fill BFTRACE w/ appropriate TRACER ID #'s
      ! Replaced NOx with NO (mpayer, 3/14/13)
      IF ( IDBFNO   /= 0 ) BFTRACE(IDBFNO  ) = IDTNO
      IF ( IDBFCO   /= 0 ) BFTRACE(IDBFCO  ) = IDTCO
      IF ( IDBFALK4 /= 0 ) BFTRACE(IDBFALK4) = IDTALK4
      IF ( IDBFACET /= 0 ) BFTRACE(IDBFACET) = IDTACET
      IF ( IDBFMEK  /= 0 ) BFTRACE(IDBFMEK ) = IDTMEK
      IF ( IDBFALD2 /= 0 ) BFTRACE(IDBFALD2) = IDTALD2
      IF ( IDBFPRPE /= 0 ) BFTRACE(IDBFPRPE) = IDTPRPE
      IF ( IDBFC3H8 /= 0 ) BFTRACE(IDBFC3H8) = IDTC3H8
      IF ( IDBFCH2O /= 0 ) BFTRACE(IDBFCH2O) = IDTCH2O  
      IF ( IDBFC2H6 /= 0 ) BFTRACE(IDBFC2H6) = IDTC2H6  
      IF ( IDBFGLYX /= 0 ) BFTRACE(IDBFGLYX) = IDTGLYX
      IF ( IDBFMGLY /= 0 ) BFTRACE(IDBFMGLY) = IDTMGLY
      IF ( IDBFBENZ /= 0 ) BFTRACE(IDBFBENZ) = IDTBENZ
      IF ( IDBFTOLU /= 0 ) BFTRACE(IDBFTOLU) = IDTTOLU
      IF ( IDBFXYLE /= 0 ) BFTRACE(IDBFXYLE) = IDTXYLE
      IF ( IDBFC2H4 /= 0 ) BFTRACE(IDBFC2H4) = IDTC2H4
      IF ( IDBFC2H2 /= 0 ) BFTRACE(IDBFC2H2) = IDTC2H2
      IF ( IDBFGLYC /= 0 ) BFTRACE(IDBFGLYC) = IDTGLYC
      IF ( IDBFHAC  /= 0 ) BFTRACE(IDBFHAC ) = IDTHAC
      ! SOAupdate: for gas phase NAP chemistry, NAP bf emiss (hotp 6/24/09)
      IF ( IDBFNAP  /= 0 ) BFTRACE(IDBFNAP ) = IDTNAP

      ! Echo biofuel tracer information
      IF ( am_I_Root ) THEN
         WRITE( 6, 100 ) BFTRACE( 1:NBFTRACE )
      ENDIF
 100  FORMAT( 'TRACERID: Biofuel burning tracers        :', 20i3 )
      
      END SUBROUTINE SET_BFTRACE

!------------------------------------------------------------------------------

      SUBROUTINE INIT_BIOFUEL
!
!******************************************************************************
!  Subroutine INIT_BIOFUEL allocates and zeroes the BIOFUEL array. 
!  (bmy, 9/12/00, 10/15/02)
!
!  NOTES:
!  (1 ) Increase BIOFUEL array from 2 to NBIOFUEL=10 elements (bmy, 3/15/01)
!  (2 ) Make sure NBFTRACE > 0 before allocating BIOFUEL (bmy, 4/17/01)
!  (3 ) BIOFUEL is now declared (NBFTRACE,IIPAR,JJPAR) (bmy, 9/28/01)
!  (4 ) Deleted obsolete code from 9/01 (bmy, 11/15/01)
!  (5 ) Now references ALLOC_ERR from "error_mod.f".  Also references IDBFNOX,
!        IDBFCO, etc from "tracerid_mod.f" (bmy, 11/6/02)
!  (6 ) Replace LWOODCO w/ LBIOFUEL from "logical_mod.f" (bmy, 7/19/04)
!  (7 ) Remove reference to TRACERID_MOD, it's obsolete (bmy, 10/3/05)
!******************************************************************************
!
      ! References to F90 modules
      USE ERROR_MOD,    ONLY : ALLOC_ERR
!      USE LOGICAL_MOD,  ONLY : LBIOFUEL

      USE CMN_SIZE_MOD  ! Size parameters, etc

      ! Local variables
      INTEGER :: AS

      !=================================================================
      ! INIT_BIOFUEL begins here!
      !=================================================================
!      IF ( LBIOFUEL .and. NBFTRACE > 0 ) THEN
      IF ( NBFTRACE > 0 ) THEN
         ALLOCATE( BIOFUEL( NBFTRACE, IIPAR, JJPAR ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'BIOFUEL' )
         BIOFUEL = 0d0

         ! This is a local array to hold biofuel in kg
         ALLOCATE( BIOFUEL_KG( NBFTRACE, IIPAR, JJPAR ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'BIOFUEL_KG' )
         BIOFUEL_KG = 0d0
      ENDIF

      ! Return to calling program
      END SUBROUTINE INIT_BIOFUEL

!------------------------------------------------------------------------------

      SUBROUTINE CLEANUP_BIOFUEL
!
!******************************************************************************
!  Subroutine CLEANUP_BIOFUEL deallocates the BIOFUEL array (bmy, 9/11/00)
!******************************************************************************
!      
      ! CLEANUP_BIOFUEL begins here!
      IF ( ALLOCATED( BIOFUEL    ) ) DEALLOCATE( BIOFUEL    )
      IF ( ALLOCATED( BIOFUEL_KG ) ) DEALLOCATE( BIOFUEL_KG )

      ! Return to calling program
      END SUBROUTINE CLEANUP_BIOFUEL

!------------------------------------------------------------------------------
      
      END MODULE BIOFUEL_MOD
