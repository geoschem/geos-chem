!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: ndxx_setup.F
!
! !DESCRIPTION: Subroutine NDXX\_SETUP dynamically allocates memory for 
!  certain diagnostic arrays that  are declared allocatable in "diag\_mod.f". 
!\\
!\\
!  This allows us to reduce the amount of memory that needs to be declared 
!  globally.  We only allocate memory for arrays if the corresponding 
!  diagnostic is turned on.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE NDXX_SETUP( am_I_Root, Input_Opt, State_Chm, RC )
!
! !USES:
!
#if defined( BPCH_DIAG ) || defined( MODEL_GEOS )
      USE CMN_DIAG_MOD
      USE CMN_SIZE_MOD
      USE DIAG_MOD
#endif
      USE DIAG_OH_MOD,        ONLY : INIT_DIAG_OH
      USE ErrCode_Mod
      USE ERROR_MOD,          ONLY : ALLOC_ERR,   ERROR_STOP
      USE Input_Opt_Mod,      ONLY : OptInput
      USE PLANEFLIGHT_MOD,    ONLY : SETUP_PLANEFLIGHT
      USE State_Chm_Mod,      ONLY : ChmState
      USE State_Chm_Mod,      ONLY : Ind_
      USE State_Met_Mod,      ONLY : MetState
#if   defined( TOMAS )
      USE TOMAS_MOD,          ONLY : IBINS, ICOMP, IDIAG   !(win, 7/9/09)
#endif

      IMPLICIT NONE
!
! !INPUT PARAMETERS:
!
      LOGICAL,        INTENT(IN)     :: am_I_Root   ! Are we on the root CPU?
!
! !INPUT/OUTPUT PARAMETERS:
!
      TYPE(OptInput), INTENT(INOUT)  :: Input_Opt   ! Input Options object
      TYPE(ChmState), INTENT(INOUT)  :: State_Chm   ! Chemistry state object
!
! !OUTPUT PARAMETERS:
!
      INTEGER,        INTENT(OUT)    :: RC          ! Success or failure?
! 
! !REVISION HISTORY:
!  16 Jun 1998 - I. Bey, R. Yantosca - Initial version
!  (1 ) This subroutine was split off from subroutine INPUT, for clarity
!  (2 ) Added call to READ49 (bey, 2/99)
!  (3 ) Eliminate GISS-Specific code, and AIJ, AIL diagnostics (bmy, 3/15/99)
!  (4 ) Define tracer offset TRCOFFSET for "alternate chemistry" runs.
!  (5 ) Multi-level diagnostics ND21, ND22, ND43, ND45, ND66, and ND68 have 
!        now been split off from the AIJ arrays (bmy, 3/29/99)
!  (6 ) Added code for ND14 and ND15.  Also eliminated obsolete code
!        and updated comments (bmy, 11/10/99)
!  (7 ) Added new ND41 and ND51 diagnostics (from amf).  Freed up obsolete
!        diagnostics ND34. ND37, and ND42 and updated comments. (bmy, 11/15/99)
!        Also note: ND41 uses allocatable array AD41. (bmy, 12/6/99)
!  (8 ) The following diagnostic arrays are now declared allocatable
!        in "diag_mod.f": AD21, AD22, AD38, AD39, AD43, AD45, AD47, 
!        AD66, AD68, CONVFLUP, TURBFLUP, MASSFLEW, MASSFLNS, MASSFLUP, TCOBOX
!        Allocate memory for these arrays only if their respective
!        diagnostic is turned on.  This will save memory. (bmy, 11/29/99)
!  (9 ) Added ND55 diagnostic for tropopause heights (hyl, bmy, 12/1/99)
!  (10) ND50 and ND20 now have dynamically allocatable arrays. (bmy, 1/5/00)  
!  (11) ND27 diagnostic now also turns on ND24, ND25, ND26 (bmy, 1/7/00)
!  (12) ND31, ND33, ND35, ND37, ND67, and ND69 now use dynamically 
!        allocatable arrays declared in "diag_mod.f". (bmy, 2/17/00)
!  (13) ND16, ND17, ND18 now use allocatable arrays.  Also now use internal
!        subroutine "alloc_err" to print error messages. (bmy, 3/14/00)
!  (14) AIJ is now obsolete.  All diagnostic variables now use allocatable
!        arrays (cf. "diag_mod.f").  This is necessary in order to keep the
!        size of the 2 x 2.5 executable within machine limits. (bmy, 3/28/00)
!  (15) Removed obsolete code.  Added TRCOFFSET of 3 for CO run
!        with parameterized OH.  Removed reference to KAIJPAR. (bmy, 4/19/00)
!  (16) Add TRCOFFSET of 50 for DMS/SO2/SO4/MSA.  Also added arrays for
!        ND13 diagnostic for sulfur emissions (bmy, 6/6/00)
!  (17) Add reference to F90 module "biomass_mod.f".  Also added array
!        AD32_bf for biofuel NOx. (bmy, 9/11/00)
!  (18) Use NTRACE + 2 prodloss families for Tagged CO for the
!        ND65 diagnostic (bmy, 10/6/00)
!  (19) Adjust TRCOFFSET for 10-tracer Tagged CO run.  Redimensioned 
!        AD45 and AD47 to save memory.  Renamed STATUS to AS. (bmy, 10/18/00)
!  (20) Removed obsolete code from 10/00.  Save out ND65 only to LLTROP 
!        levels for full chemistry.  Save out ND43 only to LLTROP levels 
!        for full chemistry.  Dimension DIAGCHLORO up to LLTROP for 
!        full chemistry (or LLPAR for CO/OH chemistry).  ND24, ND25, ND26 
!        can now save out less than LLPAR levels.  Eliminate dependence 
!        on PD35, PD37, PD39 parameters (bmy, 12/5/00)
!  (21) Only save out a maximum of LCONVM layers for ND14 (bmy, 12/7/00)
!  (22) Removed obsolete code from 7/00, 9/00, and 12/00 (bmy, 12/21/00)
!  (23) Increase to NTRACE + 4 prodloss families for Tagged CO (bmy, 1/2/01)
!  (24) Add TRCOFFSET of 54 for CH4 chemistry (NSRCX == 9) (bmy, 1/16/01)
!  (25) Now allocate DIAGCHLORO (ND23 diagnostic) for CH4 runs (bmy, 1/18/01)
!  (26) For ND43, save up to LLTROP for full chemistry, but save up to
!        LLPAR for Tagged CO or CO-OH chemistry (bmy, 2/12/01)
!  (27) Now allocate AD34 for biofuel burning emissions (bmy, 3/15/01)
!  (28) Add L(CH3I) to ND65 diagnostic (nad, bmy, 3/20/01)
!  (29) For full chemistry, we only need to save up to LLTROP levels
!        for the ND22 J-value diagnostic (bmy, 4/2/01)
!  (30) Remove reference to NBIOMAX from "biomass_mod.f" (bmy, 4/17/01)
!  (31) Eliminate obsolete commented-out code (bmy, 4/20/01)
!  (32) Now also allocate the AD12 diagnostic array (bdf, bmy, 6/15/01)
!  (33) Now assign TRCOFFSET = 40 for multi-tracer Ox run (when NSRCX = 6 
!        and LSPLIT = T).  Reference CMN_SETUP for LSPLIT.  Allocate AD44 
!        with NTRACE instead of NUMDEP for single or multi-tracer Ox runs 
!        (NSRCX = 6).  Now define NFAM as NTRACE*2 for single or multi-tracer 
!        Ox runs.  Updated comments & made cosmetic changes. (bmy, 7/3/01)
!  (34) Added AD11 diagnostic for acetone source.  Also removed obsolete
!        code from 7/01. (bmy, 9/4/01)
!  (35) Turn off ND23 unless NSRCX = 3, 5, or 9.  This prevents us from
!        referencing an unallocated DIAGCHLORO array.  Add error check for
!        ND65, make sure that NFAM > 0.  Also clean up the code that 
!        allocates AD65 and FAMPL arrays. (bmy, 1/14/02)
!  (36) Now set TRCOFFSET = 64 for tagged C2H6 chemistry (bmy, 1/25/02)
!  (37) Eliminate obsolete code from 1/02 and 2/02.  Also allocate LTNO2,
!        CTNO2, LTHO2, CTHO2 for the ND43 diagnostic. (bmy, 2/27/02)
!  (38) Call SETUP_PLANEFLIGHT to initialize the ND40 plane flight diagnostic
!        for non-SMVGEAR chemistry runs. (mje, bmy, 7/2/02)
!  (39) Now set up variables & arrays for ND01 and ND02 diagnostics (i.e.
!        Rn-Pb-Be emissions and decay).  (bmy, 9/20/02)
!  (40) Now allocate AD05 array.   Now allocate routines ALLOC_ERR and 
!        ERROR_STOP from "error_mod.f".  Now reference NEMANTHRO from F90
!        module "tracerid_mod.f" instead of "comtrid.h".  Also added array
!        AD13_SO2_bf for biofuel SO2. (bmy, 1/16/03)
!  (41) Now also allocate AD13_NH3_na array for ND13 (rjp, bmy, 3/23/03)
!  (42) Added ND03 diagnostic for Kr85 prod/loss.  Also removed special case
!        TRCOFFSET for single-tracer Ox. (jsw, bmy, 8/20/03)
!  (43) Now use GET_WETDEP_NMAX to get max # of soluble tracers for ND37,
!        ND18, and ND19.  Also set NFAM=NTRACE+5 for Tagged CO simulation. 
!        (3/18/04)
!  (44) Now initialize AD06 and AD07* arrays (rjp, tdf, bmy, 4/5/04)
!  (45) Now initialize AD08 array.  Reset TRCOFFSET for tagged CO from
!        84 to 80.  Also activate ND52 diagnostic for ICARTT.
!        (rjp, bec, stu, cas, bmy, 4/20/04)
!  (46) Now allocate AD13_SO2_sh array for ND13 (bec, bmy, 5/20/04)
!  (47) Now allocate AD07_HC array for ND07 (rjp, bmy, 7/13/04)
!  (48) Now references "tracer_mod.f" and "logical_mod.f" instead of "CMN"
!        and "CMN_SETUP".  Now references INIT_DIAG_OH from "diag_oh_mod.f"
!        Adjust TRCOFFSET for various aerosol simulations. (bmy, 7/20/04)
!  (49) Make sure ND21 only goes from 1-LLTROP (bmy, 9/28/04)
!  (50) Now allocate AD13_SO4_bf array (bmy, 11/17/04)
!  (51) Now allocate extra arrays for ND03 mercury diag.  Also set up for
!        mercury tracers in ND44 diagnostic. (bmy, 12/14/04)
!  (52) Added separate ND21 array for cryst sulfur tracers.  Now reinstated
!        AD03 array for mercury simulation.  Now move ND03 diagnostics into
!        a separate module.  Remove TCOBOX reference, it's obsolete.
!        (cas, sas, bmy, 1/21/05)
!  (53) Now remove references to AD41 & AFTTOT.  Now call SETUP_PLANEFLIGHT 
!        for non-full-chemistry runs in main.f -- this will allow it to look 
!        for flight files for each day (bmy, 3/24/05)
!  (54) Now use PD05=10 to dimension AD05 array (bmy, 4/13/05)
!  (55) Now also allocates AD09 and AD09_em (bmy, 6/27/05)
!  (56) Now allocates AD30 (bmy, 8/18/05)
!  (57) Removed duplicate variable declarations (bmy, 2/6/06)
!  (58) Now remove NBIOTRCE; it's obsolete.  Replace w/ NBIOMAX (bmy, 4/5/06)
!  (59) Now remove TRCOFFSET; it's obsolete (bmy, 5/16/06)
!  (60) Added the ND54 for time spend in the troposphere (phs, 10/17/06)
!  (61) Now allocate ND43 and ND45 counter arrays as 3-D (phs, 1/19/07)
!  (62) For ND20 diagnostic, reset ND65 diagnostic with LLTROP_FIX instead of 
!        LLTROP.  Added ND10 diagnostic setup.  Added modifications for H2-HD 
!        simulation. (phs, bmy, 9/18/07)
!  (63) Now save true pressure edges for ND31 diagnostic (bmy, 11/16/07)
!  (64) Now stop the run if ND20 is defined but ND65 isn't (bmy, 12/4/07)
!  (65) Allocate CTO3_24h (phs, 11/18/08)      
!  (66) We don't need to set LD65=1 here anymore, we now call NDXX_SETUP!
!        after DIAG_PL_MOD. (phs, bmy, 12/18/08)
!  (67) Added ND52 for GAMMA HO2 diagnostic. (ccc, jaegle, 2/26/09)
!  (68) Add AD07_SOAGM (tmf, 1/7/09) 
!  (67) Added ND52 for GAMMA HO2 diagnostic. (ccc, jaegle, 2/26/09)
!  (68) Add AD07_SOAGM (tmf, 1/7/09) 
!  (69) Now always allocate Mass Flux arrays (phs, 4/15/09)
!  (70) Allocate LTO3. (ccc, 7/20/09)
!  (71) Add AD19, AD58, AD60 (kjw, 8/18/09)
!  (72) Now AD13_SO2_an and AD13_SO4_an have NOXLEVELS levels to accomodate
!     NEI 2005 (amv, 10/9/09)
!  (73) AD13_NH3_an is 3D now (phs, 10/22/09)
!  (74) Add new diagnostic ND59, ND60, ND61 (win, 7/9/09) 
!  (75) Increase size for AD44 for TOMAS aerosol mass (win, 7/14/09)
!  (76) Initialize values for LD59, LD60, and LD61 (win, 8/10/09)
!  (77) NBIOMAX is now in CMN_SIZE. (fp, 2/26/10)
!  26 Aug 2010 - R. Yantosca - Added ProTeX headers
!  16 Feb 2011 - R. Yantosca - Add modifications for APM from G. Luo
!  09 Nov 2012 - R. Yantosca - Added GIGC-specific modifications
!  29 Mar 2013 - R. Yantosca - Pass objects to GET_WETDEP_NMAX
!  02 Apr 2013 - M. Payer    - Remove allocation of *NO, *NO2, and *NO3 arrays
!                              for ND43. These are no longer needed because NO,
!                              NO2, and NO3 are now tracers.
!  13 Aug 2013 - M. Sulprizio- Modify AD07_HC for updated SOA (H. Pye)
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
!  12 Sep 2013 - M. Sulprizio- Modify AD06 to include dust alkalinity tracers
!                              (T.D. Fairlie)
!  08 Nov 2013 - M. Sulprizio- Removed CTO3 and LTO3. They are no longer used
!                              because O3 is now a tracer.
!  15 Aug 2014 - R. Yantosca - Removed reference to biofuel_mod.F
!  15 Aug 2014 - R. Yantosca - Now reference all arrays in diag_mod.F
!  22 May 2015 - R. Yantosca - Removed variables made obsolete by HEMCO
!  03 Sep 2015 - R. Yantosca - Now pass State_Chm to this routine so that
!                              we can take advantage of the species database
!  20 Jan 2016 - M. Sulprizio- Remove code to allocate the AD06, AD07, and AD08
!                              arrays. These diagnostics are tracked by HEMCO.
!  23 Jun 2016 - R. Yantosca - Remove references to APM code; it is no longer
!                              compatible with the FlexChem implementation
!  20 Jul 2016 - R. Yantosca - Remove references to NNPAR (now use nAdvect)
!  20 Sep 2016 - R. Yantosca - Bug fix for Gfortran: Make sure the IT_IS_*
!                              variables are declared as LOGICAL, not INTEGER
!  03 Oct 2016 - R. Yantosca - Bug fix: don't allocate ND71 scalar variables
!  09 Mar 2018 - R. Yantosca - Now use LLPAR to size diagnostic arrays, 
!                              which facilitates bpch vs.netCDF comparisons
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER :: NMAX, AS, NEMISS, LMAX, nAdvect, NFAM
      LOGICAL :: IT_IS_A_CH3I_SIM
      LOGICAL :: IT_IS_A_FULLCHEM_SIM
      LOGICAL :: IT_IS_A_MERCURY_SIM
      LOGICAL :: IT_IS_A_TAGO3_SIM
      LOGICAL :: IT_IS_A_H2HD_SIM
      LOGICAL :: LDUST
      LOGICAL :: LCARB
      LOGICAL :: LSSALT
      LOGICAL :: LDRYD
      LOGICAL :: LGTMM

      !=================================================================
      ! NDXX_SETUP begins here! 
      !=================================================================
#if defined( BPCH_DIAG )

      nAdvect              = State_Chm%nAdvect
      LDUST                = Input_Opt%LDUST
      LCARB                = Input_Opt%LCARB
      LSSALT               = Input_Opt%LSSALT
      LDRYD                = Input_Opt%LDRYD
      LGTMM                = Input_Opt%LGTMM
      IT_IS_A_CH3I_SIM     = Input_Opt%ITS_A_CH3I_SIM
      IT_IS_A_FULLCHEM_SIM = Input_Opt%ITS_A_FULLCHEM_SIM
      IT_IS_A_MERCURY_SIM  = Input_Opt%ITS_A_MERCURY_SIM  
      IT_IS_A_TAGO3_SIM    = Input_Opt%ITS_A_TAGO3_SIM 
      IT_IS_A_H2HD_SIM     = Input_Opt%ITS_A_H2HD_SIM

      ! Assume success
      RC                   = GC_SUCCESS

      !=================================================================
      ! Initialize some multi-level variables       
      !=================================================================
      LD01 = 1
      LD02 = 1
      LD05 = 1
      LD09 = 1
      LD10 = 1
      LD07 = 1
      LD12 = 1
      LD13 = 1
      LD14 = 1
      LD15 = 1
      LD16 = 1
      LD17 = 1
      LD18 = 1
      LD19 = 1
      LD21 = 1
      LD22 = 1
      LD24 = 1
      LD25 = 1
      LD26 = 1
      LD31 = 1
      LD37 = 1
      LD38 = 1
      LD39 = 1
      LD43 = 1
      LD45 = 1
      LD47 = 1
      LD52 = 1
      LD54 = 1
      LD66 = 1
      LD68 = 1

      !=================================================================
      ! ND01: Rn, Pb, Be emissions
      !=================================================================
      IF ( ND01 > 0 ) THEN 
         LD01 = MIN( ND01, LLPAR )
         ALLOCATE( AD01( IIPAR, JJPAR, LD01 ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD01' )
      ENDIF

      !=================================================================
      ! ND02: Rn, Pb, Be decay
      !=================================================================
      IF ( ND02 > 0 ) THEN 
         LD02 = MIN( ND02, LLPAR )
         ALLOCATE( AD02( IIPAR, JJPAR, LD02, nAdvect ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD02' )
      ENDIF

      !=================================================================
      ! ND04: CO2 source - see ??
      ! 
      ! ND05: Sulfate Prod/loss
      !=================================================================
      IF ( ND05 > 0 ) THEN

         ! Now size ND05 for all vertical levels (bmy, 3/9/18)
         LD05 = MIN( ND05, LLPAR )

         ALLOCATE( AD05( IIPAR, JJPAR, LD05, PD05 ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD05' )
      ENDIF

      !=================================================================
      ! ND07: Carbonaceous aerosols emissions and chemical conversion
      !=================================================================
      IF ( ND07 > 0 .and. LCARB ) THEN 
         LD07 = MIN( ND07, LLPAR )

         ALLOCATE( AD07_BC( IIPAR, JJPAR, LD07 ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD07_BC' )

         ALLOCATE( AD07_OC( IIPAR, JJPAR, LD07 ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD07_OC' )

         ! SOAupdate (hotp 6/15/09)
         ! semivolpoa2: increase dimension for POA (hotp 3/3/09)
         ! semivolpoa4: increase dimension for OPOA(hotp 3/27/09)
         ! semivolpoa4: add dimension for POG ox to OPOG (hotp 3/28/09)
         ! NAPSOA: increase AD07_HC dimension by one to 9 for production
         ! of SOA from IVOCs (hotp 7/22/09)
         ALLOCATE( AD07_HC( IIPAR, JJPAR, LD07, 9 ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD07_HC' )
      ENDIF  

      !=================================================================
      ! ND09: HCN / CH3CN source & sink
      !=================================================================
      IF ( ND09 > 0 ) THEN
         LD09 = MIN( ND09, LLPAR )

!         ALLOCATE( AD09( IIPAR, JJPAR, LD09, nAdvect ), STAT=AS )
!         IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD09' )
!
!         ALLOCATE( AD09_em( IIPAR, JJPAR, PD09 ), STAT=AS )
!         IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD09_em' )
      ENDIF

      !=================================================================
      ! ND12: Distribution of emissions in boundary layer [fraction]
      !       --> uses AD12 array (allocatable)
      !=================================================================
      LD12 = MIN( ND12, LLPAR )

      IF ( ND12 > 0 ) THEN 
         ALLOCATE( AD12( IIPAR, JJPAR, LD12 ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'PL24H' )
      ENDIF

      !=================================================================
      ! ND14: Upward flux of from wet conv [kg/s] 
      !       --> uses CONVFLUP array (allocatable)
      !=================================================================
      IF ( ND14 > 0 ) THEN
!------------------------------------------------
! Prior to 3/19/18:
! Use LLPAR instead of LLCONVM (bmy, 3/19/18)
!         LD14 = MIN( ND14, LLCONVM )
!------------------------------------------------
         LD14 = MIN( ND14, LLPAR )
         NMAX = nAdvect

         ALLOCATE( CONVFLUP( IIPAR, JJPAR, LLPAR, NMAX ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'CONVFLUP' )
      ENDIF

      !=================================================================
      ! ND15: Mass change from BL-mixing [kg/s] 
      !       --> uses TURBFLUP array (allocatable)
      !=================================================================
      IF ( ND15 > 0 ) THEN
         LD15 = MIN( ND15, LLPAR )
         NMAX = nAdvect

         ALLOCATE( TURBFLUP( IIPAR, JJPAR, LLPAR, NMAX ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'TURBFLUP' )
      ENDIF

      !=================================================================
      ! ND16: Fraction of grid box experiencing large-scale and
      !       convective precipitation --> uses AD16 array (allocatable)
      !=================================================================
      IF ( ND16 > 0 ) THEN
         LD16 = MIN( ND16, LLPAR )

         ! Store LS Precip fraction only
         ALLOCATE( AD16( IIPAR, JJPAR, LD16, 1 ), STAT=AS ) 
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD16' )

         ! Counter array for AD16
         ALLOCATE( CT16( IIPAR, JJPAR, LD16, 1 ), STAT=AS ) 
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'CT16' )
      ENDIF

      !=================================================================
      ! ND17: Fraction of tracer lost to rainout (in both large-scale
      !       and conv precipitation) --> uses AD17 array (allocatable)
      !=================================================================
      IF ( ND17 > 0 ) THEN

         ! Get # of soluble tracers for this simulation
         NMAX = State_Chm%nWetDep

         ! Turn off ND17 if there are no wetdep species
         IF ( NMAX == 0 ) THEN 
            ND17 = 0
            IF ( am_I_Root ) THEN 
               print*, '### No wetdep species, so turn off ND17'
            ENDIF
         ENDIF

         ! Number of levels
         LD17 = MIN( ND17, LLPAR )

         ! Store both LS precip only
         ALLOCATE( AD17( IIPAR, JJPAR, LD17, NMAX, 1 ), STAT=AS ) 
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD17' )

         ! Counter array for AD17
         ALLOCATE( CT17( IIPAR, JJPAR, LD17, 1 ), STAT=AS ) 
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'CT17' )
      ENDIF

      !=================================================================
      ! ND18: Fraction of tracer lost to washout (in both large-scale
      !       and convective precipitation) --> uses AD18 array (alloc.)
      !=================================================================
      IF ( ND18 > 0 ) THEN

         ! Get # of soluble tracers for this simulation
         NMAX = State_Chm%nWetDep

         ! Turn off ND17 if there are no wetdep species
         IF ( NMAX == 0 ) THEN 
            ND18 = 0
            IF ( am_I_Root ) THEN 
               print*, '### No wetdep species, so turn off ND18'
            ENDIF
         ENDIF

         ! Number of levels
         LD18 = MIN( ND18, LLPAR )

         ! Store both LS precip fraction only
         ALLOCATE( AD18( IIPAR, JJPAR, LD18, NMAX, 1 ), STAT=AS) 
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD18' )

         ! Counter array for AD17
         ALLOCATE( CT18( IIPAR, JJPAR, LD18, 1 ), STAT=AS ) 
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'CT18' )
      ENDIF

      !=================================================================
      ! ND19: CH4 loss by OH 
      !=================================================================
      IF ( ND19 > 0 ) THEN 
         LD19 = MIN( ND19, LLPAR )
	 
         ALLOCATE( AD19( IIPAR, JJPAR, LD19, PD19 ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD19' )
      ENDIF 

      !=================================================================
      ! ND20: Save O3 P-L losses to disk for single-tracer O3 run
      !       in the PL24H array.  Also turn on ND65, since the P-L 
      !       rates are computed by ND65.  
      !=================================================================
      IF ( ND20 > 0 ) THEN 
         IF ( ND65 == 0 ) THEN
            CALL ERROR_STOP( 'ND65 must be turned on for ND20 output!',
     &                       'ndxx_setup.f'  )
         ENDIF
      ENDIF
#endif
#if defined( BPCH_DIAG ) || defined( MODEL_GEOS )
      !=================================================================
      ! ND21: Optical depths and cloud fractions [unitless]
      !       --> uses AD21 array (allocatable) 
      !=================================================================
      IF ( Input_Opt%ND21 > 0 ) THEN
         ! Now size diagnostic arrays w/ all vertical levels (bmy, 3/9/18)
         Input_Opt%LD21 = MIN( Input_Opt%ND21, LLCHEM )

         ! For regular 
         ALLOCATE( AD21( IIPAR, JJPAR, Input_Opt%LD21, PD21 ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD21' )
      ENDIF
#endif
#if defined( BPCH_DIAG )

      !=================================================================
      ! ND22: J-value diagnostics [s^-1] 
      !       --> uses AD22 array (allocatable)
      !=================================================================
      IF ( ND22 > 0 ) THEN

         ! Now size J-value diagnostic for all vertical levels (bmy, 3/9/18)
         LD22 = MIN( ND22, LLPAR  )

         ! Accumulating diagnostic array
         ALLOCATE( AD22( IIPAR, JJPAR, LD22, PD22 ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD22' )

         ! Number of times where LT is between HR1_JV and HR2_JV
         ALLOCATE( CTJV( IIPAR, JJPAR ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'CTJV' )         
      ENDIF

      ! Locations where LT is between HR1_JV and HR2_JV
      ! NOTE: LTJV is currently used for netCDF diagnostics in flexchem_mod.F90
      ! regardless of the ND22 setting, so we always need to allocate it here
      ! (mps, 1/5/18)
      ALLOCATE( LTJV( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'LTJV' )         

      !=================================================================
      ! ND27: Flux of Ox across the annual mean tropopause [kg/s]
      !       ND27 will also turn on ND24, ND25, ND26 diagnostics
      !=================================================================
      IF ( ND27 > 0 ) THEN
         ND24 = LLPAR
         ND25 = LLPAR
         ND26 = LLPAR
      ENDIF

      !=================================================================
      ! ND24: Eastward mass flux from transport [kg/s] 
      !       --> uses MASSFLEW array (allocatable)
      !=================================================================
      IF ( ND24 > 0 ) THEN
         LD24 = MIN( ND24, LLPAR )
         NMAX = nAdvect
         ALLOCATE( MASSFLEW( IIPAR, JJPAR, LLPAR, NMAX ), STAT=AS) 
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'MASSFLEW' )
      ELSE
         ALLOCATE( MASSFLEW( 1, 1, 1, 1 ), STAT=AS) 
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'MASSFLEW' )
      ENDIF

      !=================================================================
      ! ND25: Northward mass flux from transport [kg/s] 
      !       --> uses MASSFLNS array (allocatable)
      !=================================================================
      IF ( ND25 > 0 ) THEN
         LD25 = MIN( ND25, LLPAR )
         NMAX = nAdvect
         ALLOCATE( MASSFLNS( IIPAR, JJPAR, LLPAR, NMAX ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'MASSFLNS' ) 
      ELSE
         ALLOCATE( MASSFLNS( 1, 1, 1, 1 ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'MASSFLNS' ) 
      ENDIF

      !=================================================================
      ! ND26: Vertical mass flux from transport [kg/s] 
      !       --> uses MASSFLUP array (allocatable)
      !=================================================================
      IF ( ND26 > 0 ) THEN
         LD26 = MIN( ND26, LLPAR )
         NMAX = nAdvect
         ALLOCATE( MASSFLUP( IIPAR, JJPAR, LLPAR, NMAX ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'MASSFLUP' )
      ELSE
         ALLOCATE( MASSFLUP( 1, 1, 1, 1 ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'MASSFLUP' )
      ENDIF

      !=================================================================
      ! ND30: Land/water/ice flags
      !=================================================================      
      IF ( ND30 > 0 ) THEN
         ALLOCATE( AD30( IIPAR, JJPAR ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD30' )
      ENDIF

      !=================================================================
      ! ND31: 3-D pressure edges [hPa] --> Uses AD31 array (allocatable)
      !=================================================================
      IF ( ND31 > 0 ) THEN
         LD31 = MIN( ND31, LLPAR+1 )
         ALLOCATE( AD31( IIPAR, JJPAR, LD31 ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD31' )
      ENDIF

      !=================================================================
      ! ND33: Column sum of tracer [kg]
      !       --> uses AD33 array (allocatable)!
      !=================================================================
      IF ( ND33 > 0 ) THEN
         ALLOCATE( AD33( IIPAR, JJPAR, nAdvect ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD33' )
      ENDIF

      !=================================================================
      ! ND35: Tracer at 500 mb [v/v] (this is ~ level 9 for GEOS-CHEM)
      !       --> uses AD35 array (allocatable)
      !=================================================================
      IF ( ND35 > 0 ) THEN
         ALLOCATE( AD35( IIPAR, JJPAR, nAdvect ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD35' ) 
      ENDIF

      !=================================================================
      ! ND37: Fraction of tracer scavenged in cloud updrafts 
      !       --> Uses AD37 array (allocatable)
      !=================================================================
      IF ( ND37 > 0 ) THEN
         
         ! Get # of soluble tracers for this simulation
         NMAX = State_Chm%nWetDep

         ! Turn off ND37 if there are no wetdep species
         IF ( NMAX == 0 ) THEN 
            ND37 = 0
            IF ( am_I_Root ) THEN 
               print*, '### No wetdep species, so turn off ND37'
            ENDIF
         ENDIF

         ! Number of levels
         LD37 = MIN( ND37, LLPAR )

         ! Allocate array accordingly
         IF ( NMAX > 0 ) THEN
            ALLOCATE( AD37( IIPAR, JJPAR, LD37, NMAX ), STAT=AS )
            IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD37' )
         ENDIF
      ENDIF

      !=================================================================
      ! ND38: Rainout of tracer (nfcldmx.f) 
      !       --> uses AD38 array (allocatable)
      !=================================================================
      ! Allocate AD38 for GTMM. (ccc, 10/29/09)
      IF ( ND38 > 0 .OR. LGTMM ) THEN

         ! Get # of soluble tracers for this simulation
         NMAX = State_Chm%nWetDep
         LD38 = MIN( ND38, LLPAR )

         ! If GTMM on and LD38 off, ND38 should be LLPAR, (ccc, 10/29/09)
         IF ( ND38 == 0 ) LD38 = LLPAR
         
         ! Allocate AD38 array accordingly
         IF ( NMAX > 0 ) THEN
            ALLOCATE( AD38( IIPAR, JJPAR, LD38, NMAX ), STAT=AS )
            IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD38' )
         ENDIF
      ENDIF

      !=================================================================
      ! ND39: Rainout of tracer (wetdep.f) 
      !       --> uses AD39 array (allocatable)
      !=================================================================
      ! Allocate AD39 for LGTMM. (ccc, 10/29/09)
      IF ( ND39 > 0 .OR. LGTMM ) THEN

         LD39 = MIN( ND39, LLPAR )
         ! If GTMM on and ND38 off, ND38 should be LLPAR, (ccc, 10/29/09)
         IF ( ND39 == 0 ) LD39 = LLPAR

         ! Get # of soluble tracers for this simulation
         NMAX = State_Chm%nWetDep

         ! Allocate AD39 array accordingly
         IF ( NMAX > 0 ) THEN
            ALLOCATE( AD39( IIPAR, JJPAR, LD39, NMAX ), STAT=AS )
            IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD39' )
         ENDIF
      ENDIF

      !=================================================================
      ! ND42: SOA concentration
      !
      ! ND43: Chemical diagnostics: OH [molec/cm3/s] and HO2 [v/v]
      !       --> uses AD43 array (allocatable)
      !=================================================================
      IF ( ND43 > 0 ) THEN
         
         ! Now size diagnostic arrays with all vertical levels (bmy, 3/9/18)
         LD43 = MIN( ND43, LLPAR  )

         ! Accumulating diagnostic array
         ALLOCATE( AD43( IIPAR, JJPAR, LD43, PD43 ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD43' )

         ! Locations where LT is between HR1_OH and HR2_OH
         ALLOCATE( CTOH( IIPAR, JJPAR, LD43 ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'CTOH' )

         ! Locations where LT is between HR1_OH and HR2_OH
         ALLOCATE( CTHO2( IIPAR, JJPAR, LD43 ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'CTHO2' )

         ! Locations where LT is between HR1_OH and HR2_OH
         ALLOCATE( CTO1D( IIPAR, JJPAR, LD43 ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'CTO1D' )

         ! Locations where LT is between HR1_OH and HR2_OH
         ALLOCATE( CTO3P( IIPAR, JJPAR, LD43 ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'CTO3P' )
      ENDIF

      ! The following are also used for netcdf diagnostics and therefore
      ! must come out of the BPCH_DIAG ifdef and not be predicated on
      ! an NDXX flag. Additional work should be done to only allocated
      ! these arrays if either the bpch or netcdf diagnostic is on 
      ! (ewl, 1/5/2018)

      ! Locations where LT is between HR1_OH and HR2_OH
      ALLOCATE( LTOH( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'LTOH' )

      ! Locations where LT is between HR1_OH and HR2_OH
      ALLOCATE( LTHO2( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'LTHO2' )

      ! Locations where LT is between HR1_OH and HR2_OH
      ALLOCATE( LTO1D( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'LTO1D' )

      ! Locations where LT is between HR1_OH and HR2_OH
      ALLOCATE( LTO3P( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'LTO3P' )

      !=================================================================
      ! ND44: Drydep fluxes [s-1] and drydep velocities [cm/s]
      !       --> uses AD44 arrays (allocatable)
      !=================================================================

      ! Turn off ND44 if drydep is turned off
      IF ( .not. LDRYD ) ND44 = 0

      IF ( ND44 > 0 ) THEN
#if   defined( TOMAS )
         ! add space in diag array for TOMAS aerosol mass (win, 7/14/09)
         ! Now use State_Chm%nDryDep for # dry depositing species
         ! (ewl, 10/14/15)
         IF ( Ind_('NK1') > 1 ) THEN
            NMAX = State_Chm%nDryDep + ( ICOMP - IDIAG )* IBINS
         ELSE
            NMAX = State_Chm%nDryDep
         ENDIF
#else
         ! Number of dry depositing species
         NMAX = State_Chm%nDryDep
#endif
         ! Allocate AD44 array
         ALLOCATE( AD44( IIPAR, JJPAR, NMAX, 2 ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD44' )

!------------------------------------------------------------------------------
! This diagnostic requires some fixes. Comment out for now (mps,3/2/18)
!         ALLOCATE( AD44b( IIPAR, JJPAR, NMAX), STAT=AS )
!         IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD44b' )
!
!	 ALLOCATE( AD_RA(IIPAR, JJPAR), STAT=AS )
!         IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD_RA')
!
!	 ALLOCATE( AD_ZL(IIPAR, JJPAR), STAT=AS )
!         IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD_ZL')
!------------------------------------------------------------------------------

      ENDIF
      
      !=================================================================
      ! ND45: Tracer concentrations [v/v] between HR1_OTH and HR2_OTH
      !       --> uses AD45 array (allocatable)
      !=================================================================
      IF ( ND45 > 0 ) THEN
         LD45 = MIN( ND45, LLPAR )
         NMAX = nAdvect

         ! Accumulating diagnostic array
         ! Resize to NMAX to save memory (bmy, 10/18/00)
         ALLOCATE( AD45( IIPAR, JJPAR, LD45, NMAX ), STAT=AS )
         IF ( AS > 0 ) CALL ALLOC_ERR( 'AD45' )

         ! Locations where LT is between HR1_OTH and HR2_OTH
         ALLOCATE( LTOTH( IIPAR, JJPAR ), STAT=AS )
         IF ( AS > 0 ) CALL ALLOC_ERR( 'LTOTH' )
    
         ! Number of times LT is between HR1_OTH and HR2_OTH
         ALLOCATE( CTOTH( IIPAR, JJPAR ), STAT=AS )
         IF ( AS > 0 ) CALL ALLOC_ERR( 'CTOTH' )
      ENDIF

      !=================================================================
      ! ND47: 24-h averaged tracer concentration [v/v]
      !       --> uses AD47 array (allocatable)
      !
      ! NOTE: ND47 is always a 24-h average field, while ND45 
      !       can be averaged over any arbitrary time period.
      !=================================================================
      IF ( ND47 > 0 ) THEN
         LD47 = MIN( ND47, LLPAR   )
         NMAX = nAdvect

         ! Resize to NMAX to save memory (bmy, 10/18/00)
         ALLOCATE( AD47( IIPAR, JJPAR, LD47, NMAX ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD47' )
      ENDIF

      !=================================================================
      ! ND52: gamma values
      !=================================================================
      IF ( ND52 > 0 ) THEN
         LD52 = MIN( ND52, LLPAR)

         ALLOCATE( AD52( IIPAR, JJPAR, LD52, PD52 ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD52' )
      ENDIF

      !=================================================================
      ! ND53: POPs emissions
      !
      ! ND54 - Time spend in the troposphere
      !        --> uses AD54 array (allocatable)
      !=================================================================
      IF ( ND54 > 0 ) THEN

         ! Now size diagnostic arrays w/ all vertical levels (bmy, 3/9/18)
!----------------------------------------------
! Prior to 3/19/18:
! Now set max value as LLPAR (bmy, 3/19/18)
!         LD54 = MIN( ND54, LLCHEM)
!----------------------------------------------
         LD54 = MIN( ND54, LLPAR )

         ALLOCATE( AD54( IIPAR, JJPAR, LD54 ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD54' )
      ENDIF

      !=================================================================
      ! ND55: Tropopause diagnostics [level, height, and pressure]
      !       --> uses AD55 array (allocatable)
      !=================================================================
      IF ( ND55 > 0 ) THEN
         ALLOCATE( AD55( IIPAR, JJPAR, PD55 ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD55' )
      ENDIF

      !=================================================================
      ! ND56: Lightning flash rate diagnostics
      !=================================================================

      !FP_ISOP (6/2009)
      !=================================================================
      ! ND57: theta
      !=================================================================
      IF ( ND57 > 0 ) THEN
         LD57 = MIN( ND57, LLPAR)
         ALLOCATE( AD57( IIPAR, JJPAR, LD57 ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD57' )
      ENDIF

#if   defined ( TOMAS )
      !=================================================================
      ! ND59: Size-resolved primary aerosol emissions      !(win, 7/9/09)
      !         Emissions to number, sulfate, sea-salt, carb, dust
      !      ----> save 3-D (I,J,1) or up to (I,J,2) 
      !=================================================================
      IF ( ND59 > 0 ) THEN
         LD59 = MIN( ND59, LLPAR )

         ! Number emission
         ALLOCATE( AD59_NUMB( IIPAR, JJPAR, 2, IBINS ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD59_NUMB' )

         ALLOCATE( AD59_SULF( IIPAR, JJPAR, 2, IBINS ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD59_SULF' )

         ALLOCATE( AD59_SALT( IIPAR, JJPAR, 2, IBINS ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD59_SALT' )

         ALLOCATE( AD59_ECIL( IIPAR, JJPAR, 2, IBINS ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD59_ECIL' )

         ALLOCATE( AD59_ECOB( IIPAR, JJPAR, 2, IBINS ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD59_ECOB' )

         ALLOCATE( AD59_OCIL( IIPAR, JJPAR, 2, IBINS ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD59_OCIL' )

         ALLOCATE( AD59_OCOB( IIPAR, JJPAR, 2, IBINS ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD59_OCOB' )

         ALLOCATE( AD59_DUST( IIPAR, JJPAR, 2, IBINS ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD59_DUST' )
      ENDIF

      !=================================================================
      ! ND60: TOMAS microphysical process rates (condensation, 
      !           coagulation, nucleation, aqueous oxidation, 
      !           error-fudging)
      !       ---> save 2-D (J,L) for 30-bin of each aerosol species
      !=================================================================
      IF ( ND60 > 0 ) THEN
         LD60 = MIN( ND60, LLPAR )

         ! Now the array dimension is IBINS*(ICOMP-IDIAG+1) because
         ! we need it for all prognostic mass species + 1 for number
         ! IDIAG = # of diagnostic species.  (win, 9/27/08)
         !Condensation rate
         ALLOCATE( AD60_COND( 1, JJPAR, LD60, IBINS*(ICOMP-IDIAG+1)), 
     &             STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD60_COND' )

         !Coagulation rate
         ALLOCATE( AD60_COAG( 1, JJPAR, LD60, IBINS*(ICOMP-IDIAG+1)), 
     &             STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD60_COAG' )

         !Nucleation rate
         ALLOCATE( AD60_NUCL( 1, JJPAR, LD60, IBINS*(ICOMP-IDIAG+1)), 
     &             STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD60_NUCL' )

         !Aqueous oxidation rate
         ALLOCATE( AD60_AQOX( 1, JJPAR, LD60, IBINS*(ICOMP-IDIAG+1)), 
     &             STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD60_AQOX' )

         !Accumulated error-fudging 
         ALLOCATE( AD60_ERROR( 1, JJPAR, LD60, IBINS*(ICOMP-IDIAG+1)), 
     &             STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD60_ERROR' )

         !SOA Condensation rate 
         ALLOCATE( AD60_SOA( 1, JJPAR, LD60, IBINS*(ICOMP-IDIAG+1)), 
     &             STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD60_SOA' )
      ENDIF

#else

      !=================================================================
      ! ND60: WETLAND FRACTION 
      !=================================================================
      IF ( ND60 > 0 ) THEN 
         ALLOCATE( AD60( IIPAR, JJPAR ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD60' )
      ENDIF
#endif


#if   defined( TOMAS )
      !=================================================================
      ! ND61: 3-D TOMAS process rate diagnostic
      !       --> Uses AD61 array (allocatable)
      !  NOTE: ND61 is used for 10-nm particle formation 
      !        rate and cluster-size nucleation rate.  So the array
      !        is declared for (IIPAR,JJPAR,LD61,2) (win, 10/6/08)
      !=================================================================
      IF ( ND61 > 0 ) THEN
         LD61 = MIN( ND61, LLPAR )

         ALLOCATE( AD61( IIPAR, JJPAR, LD61, PD61 ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD61' )
         ALLOCATE( AD61_INST( IIPAR, JJPAR, LD61, PD61 ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD61_INST' )
      ENDIF
#endif

      !=================================================================
      ! ND64: Radiative flux from FAST-JX
      !       --> Uses AD64 array (allocatable)
      !=================================================================
      IF ( ND64 > 0 ) THEN
         LD64 = MIN( ND64, LLPAR+1 )

         ALLOCATE( AD64( IIPAR, JJPAR, LD64, PD64, 3 ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD64' )
      ENDIF

      !=================================================================
      ! ND65: Production & loss rates from KPP
      !       --> Uses AD65 array (allocatable)
      !=================================================================
      IF ( ND65 > 0 ) THEN
         
         ! Now size diagnostic arrays with all vertical levels (bmy, 3/9/18)
         LD65 = MIN( ND65, LLPAR  )

         IF ( Input_Opt%DO_SAVE_PCO ) THEN
            NFAM = Input_Opt%NFAM + 2
         ELSE
            NFAM = Input_Opt%NFAM
         ENDIF

         ALLOCATE( AD65( IIPAR, JJPAR, LD65, NFAM ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD64' )
      ENDIF

      !=================================================================
      ! ND65 or ND47(O3): 24-h averaged tropospheric diagnostic 
      !=================================================================
      IF ( ND65 > 0 .OR. ND47 > 0 ) THEN

         ! Need both LD47 and LD65 set above
         LMAX = MAX( LD47, LD65 )

         ! Number of times in the troposphere 
         ALLOCATE( CTO3_24h( IIPAR, JJPAR, LMAX ), STAT=AS )
         IF ( AS > 0 ) CALL ALLOC_ERR( 'CTO3_24h' )
      ENDIF

      !=================================================================
      ! ND66: DAO 3-D fields (UWND, VWND, SPHU, TMPU, RH) 
      !       --> uses AD66 array (allocatable)
      !=================================================================
      IF ( ND66 > 0 ) THEN
         LD66 = MIN( ND66, LLPAR )
         ALLOCATE( AD66( IIPAR, JJPAR, LD66, PD66 ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD66' )
      ENDIF

      !=================================================================
      ! ND67: DAO A-3 and surface fields 
      !       --> Uses AD67 array (allocatable)
      !=================================================================
      IF ( ND67 > 0 ) THEN
         ALLOCATE( AD67( IIPAR, JJPAR, PD67 ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD67' )
      ENDIF

      !=================================================================
      ! ND68: Air mass diagnostics (BXHEIGHT, AD, AVGW, N_AIR)
      !       --> uses AD68 array (allocatable)
      !=================================================================
      IF ( ND68 > 0 ) THEN
         LD68 = MIN( ND68, LLPAR )
         ALLOCATE( AD68( IIPAR, JJPAR, LD68, PD68 ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD68' )
      ENDIF

      !=================================================================
      ! ND69: DXYP -- grid box surface areas [m^2] 
      !       --> uses AD69 array (allocatable)
      !=================================================================
      IF ( ND69 > 0 ) THEN
         ALLOCATE( AD69( IIPAR, JJPAR, PD69 ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD69' )
      ENDIF

      !=================================================================
      ! ND71: Maximum tracer concentration [v/v]
      !================================================================
      IF ( ND71 > 0 ) THEN
         NMAX = nAdvect+1

         ! Resize to NMAX to save memory (bmy, 10/18/00)
         ALLOCATE( AD71( IIPAR, JJPAR, NMAX ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD71' )

         ALLOCATE( AD71_HR( IIPAR, JJPAR, NMAX ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD71_HR' )

         ALLOCATE( AD71_DAY( IIPAR, JJPAR, NMAX ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD71_DAY' )

      ENDIF

      !=================================================================
      ! ND72: Radiative output diagnostics (TOASW, SRFSW, TOALW, SRFLW,
      !       GTOASW, GSRFSW, GTOALW, GSRFLW, ATOASW, ASRFSW, ATOALW,
      !       ASRFLW) for both clear and all sky (prefix CLR and ALL)
      !       --> uses AD72 array (allocatable)
      !=================================================================
      IF ( ND72 > 0 ) THEN
         ALLOCATE( AD72( IIPAR, JJPAR, PD72 ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD72' )
      ENDIF

      !=================================================================
      ! ND73: ISORROPIA pH and H+
      !=================================================================
      IF ( ND73 > 0 ) THEN 
         LD73 = MIN( ND73, LLPAR )

         ALLOCATE( AD73( IIPAR, JJPAR, LD73, PD73 ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD73' )
      ENDIF

      ! Save LD* fields back into the Input_Opt object (bmy, 11/9/12)
      Input_Opt%LD01  = LD01
      Input_Opt%LD02  = LD02
      Input_Opt%LD05  = LD05
      Input_Opt%LD07  = LD07
      Input_Opt%LD09  = LD09
      Input_Opt%LD10  = LD10
      Input_Opt%LD12  = LD12
      Input_Opt%LD13  = LD13
      Input_Opt%LD14  = LD14
      Input_Opt%LD15  = LD15
      Input_Opt%LD16  = LD16
      Input_Opt%LD17  = LD17
      Input_Opt%LD18  = LD18
      Input_Opt%LD19  = LD19
      Input_Opt%LD21  = LD21
      Input_Opt%LD22  = LD22
      Input_Opt%LD24  = LD24
      Input_Opt%LD25  = LD25
      Input_Opt%LD26  = LD26
      Input_Opt%LD31  = LD31
      Input_Opt%LD37  = LD37
      Input_Opt%LD38  = LD38
      Input_Opt%LD39  = LD39
      Input_Opt%LD43  = LD43
      Input_Opt%LD45  = LD45
      Input_Opt%LD47  = LD47
      Input_Opt%LD52  = LD52
      Input_Opt%LD54  = LD54
      Input_Opt%LD57  = LD57
      Input_Opt%LD59  = LD59
      Input_Opt%LD60  = LD60
      Input_Opt%LD61  = LD61
      Input_Opt%LD62  = LD62
      Input_Opt%LD64  = LD64
      Input_Opt%LD65  = LD65
      Input_Opt%LD66  = LD66
      Input_Opt%LD68  = LD68
      Input_Opt%LD73  = LD73
#endif

      ! Return w/ success
      RC = GC_SUCCESS

      END SUBROUTINE NDXX_SETUP
!EOC
