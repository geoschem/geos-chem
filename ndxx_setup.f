! $Id: ndxx_setup.f,v 1.13 2004/10/15 20:16:42 bmy Exp $
      SUBROUTINE NDXX_SETUP
!
!******************************************************************************
!  NDXX_SETUP dynamically allocates memory for certain diagnostic arrays that 
!  are declared allocatable in "diag_mod.f". (bmy, bey, 6/16/98, 9/28/04)
!
!  This allows us to reduce the amount of memory that needs to be declared 
!  globally.  We only allocate memory for arrays if the corresponding 
!  diagnostic is turned on.
!
!  NOTES:
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
!******************************************************************************
!
      ! References to F90 modules
      USE BIOMASS_MOD,     ONLY : NBIOTRCE
      USE BIOFUEL_MOD,     ONLY : NBFTRACE
      USE DIAG_MOD
      USE DIAG_OH_MOD,     ONLY : INIT_DIAG_OH
      USE DRYDEP_MOD,      ONLY : NUMDEP
      USE ERROR_MOD,       ONLY : ALLOC_ERR, ERROR_STOP
      USE LOGICAL_MOD      
      USE PLANEFLIGHT_MOD, ONLY : SETUP_PLANEFLIGHT
      USE TRACER_MOD 
      USE TRACERID_MOD,    ONLY : NEMANTHRO
      USE WETSCAV_MOD,     ONLY : GET_WETDEP_NMAX

      IMPLICIT NONE

#     include "CMN_SIZE"   ! Size parameters
#     include "CMN_DIAG"   ! Diagnostic switches & arrays

      ! Local variables
      INTEGER :: NMAX, AS

      !=================================================================
      ! NDXX_SETUP begins here! 
      !
      ! Initialize some multi-level variables       
      !=================================================================
      LD01 = 1
      LD02 = 1
      LD03 = 1
      LD05 = 1
      LD07 = 1
      LD12 = 1
      LD13 = 1
      LD14 = 1
      LD15 = 1
      LD16 = 1
      LD17 = 1
      LD18 = 1
      LD21 = 1
      LD22 = 1
      LD24 = 1
      LD25 = 1
      LD26 = 1
      LD37 = 1
      LD38 = 1
      LD39 = 1
      LD43 = 1
      LD45 = 1
      LD47 = 1
      LD64 = 1
      LD65 = 1
      LD66 = 1
      LD68 = 1

      !=================================================================
      ! Define TRCOFFSET for "alternate chemistry" runs
      ! TRCOFFSET will be used to offset punch file tracer indices
      !=================================================================

      IF ( ITS_A_RnPbBe_SIM() ) THEN
         TRCOFFSET = 60         
            
      ELSE IF ( ITS_A_CH3I_SIM() ) THEN
         TRCOFFSET = 70         

      ELSE IF ( ITS_A_HCN_SIM() ) THEN
         TRCOFFSET = 80         

      ELSE IF ( ITS_A_COPARAM_SIM() ) THEN
         TRCOFFSET = 3          

      ELSE IF ( ITS_A_TAGOX_SIM() ) THEN
         TRCOFFSET = 40         

      ELSE IF ( ITS_A_TAGCO_SIM() ) THEN
         TRCOFFSET = 80         

      ELSE IF ( ITS_A_C2H6_SIM() ) THEN
         TRCOFFSET = 64         

      ELSE IF ( ITS_A_CH4_SIM() ) THEN
         TRCOFFSET = 54         

      ELSE IF ( ITS_AN_AEROSOL_SIM() ) THEN
         TRCOFFSET = 50         

         IF ( LCARB ) THEN

            ! Tracer offset for carbon-aerosol-only simulation
            IF ( .not. ( LSULF .and. LDUST .and. LSSALT ) ) THEN
               TRCOFFSET = 58
            ENDIF
          
         ELSE IF ( LDUST ) THEN

            ! Tracer offset dust-aerosol-only simulation
            IF ( .not. ( LSULF .and. LCARB .and. LSSALT ) ) THEN
               TRCOFFSET = 62
            ENDIF            

         ELSE IF ( LSSALT ) THEN
         
            ! Tracer offset seasalt-aerosol-only simulation
            IF ( .not. ( LSULF .and. LCARB .and. LSSALT ) ) THEN
               TRCOFFSET = 66
            ENDIF    

         ENDIF

      ! Implement this more fully later...
      !ELSE IF ( ITS_A_Kr85_SIM() ) THEN
      !   TRCOFFSET = 75        
             
      ELSE
         TRCOFFSET = 0

      ENDIF

      !=================================================================
      ! ND01: Rn, Pb, Be emissions
      !=================================================================
      IF ( ND01 > 0 ) THEN 
         LD01 = MIN( ND01, LLPAR )

         ALLOCATE( AD01( IIPAR, JJPAR, LD01, N_TRACERS ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD01' )
      ENDIF

      !=================================================================
      ! ND02: Rn, Pb, Be decay
      !=================================================================
      IF ( ND02 > 0 ) THEN 
         LD02 = MIN( ND02, LLPAR )

         ALLOCATE( AD02( IIPAR, JJPAR, LD02, N_TRACERS ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD02' )
      ENDIF

      !=================================================================
      ! ND03: Kr85 prod/loss
      !=================================================================
      IF ( ND03 > 0 ) THEN 
         LD03 = MIN( ND03, LLPAR )

         ALLOCATE( AD03( IIPAR, JJPAR, LD03, PD03 ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD03' )

         ! Also need to turn on ND25 N/S flux diagnostic
         ND25 = LLPAR
      ENDIF

      !=================================================================
      ! ND04: Free diagnostic
      ! 
      ! ND05: Sulfate Prod/loss
      !=================================================================
      IF ( ND05 > 0 ) THEN
         LD05 = MIN( ND05, LLTROP )

         ALLOCATE( AD05( IIPAR, JJPAR, LD05, 10 ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD05' )
      ENDIF

      !================================================================
      ! ND06: Dust emissions
      !================================================================
      IF ( ND06 > 0 .and. LDUST ) THEN 
         ALLOCATE( AD06( IIPAR, JJPAR, NDSTBIN ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD06' )
      ENDIF

      !=================================================================
      ! ND07: Carbonaceous aerosols emissions and chemical conversion
      !=================================================================
      IF ( ND07 > 0 .and. LCARB ) THEN 
         LD07 = MIN( ND07, LLPAR )

         ALLOCATE( AD07( IIPAR, JJPAR, PD07 ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD07' )

         ALLOCATE( AD07_BC( IIPAR, JJPAR, LD07 ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD07_BC' )

         ALLOCATE( AD07_OC( IIPAR, JJPAR, LD07 ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD07_OC' )

         ALLOCATE( AD07_HC( IIPAR, JJPAR, LD07, 5 ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD07_HC' )
      ENDIF  

      !================================================================
      ! ND08: Dust emissions
      !================================================================
      IF ( ND08 > 0 .and. LSSALT ) THEN 
         ALLOCATE( AD08( IIPAR, JJPAR, PD08 ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD08' )
      ENDIF

      !=================================================================
      ! ND09 - ND10: Free diagnostics
      !
      ! ND11: Acetone source diagnostics [atoms C/cm2/s]
      !       --> uses AD11 array (allocatable)
      !=================================================================
      IF ( ND11 > 0 ) THEN 
         ALLOCATE( AD11( IIPAR, JJPAR, PD11 ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD11' )
      ENDIF
      
      !=================================================================
      ! ND12: Distribution of emissions in boundary layer [fraction]
      !       --> uses AD12 array (allocatable)
      !=================================================================
      LD12 = MIN( ND12, LLTROP )

      IF ( ND12 > 0 ) THEN 
         ALLOCATE( AD12( IIPAR, JJPAR, LD12 ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'PL24H' )
      ENDIF

      !=================================================================
      ! ND13: Sulfur emissions from DMS, SO2, and SO4 
      !=================================================================
      IF ( ND13 > 0 ) THEN 
         LD13 = MIN( ND13, LLPAR )

         ALLOCATE( AD13_DMS( IIPAR, JJPAR ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD13_DMS' )

         ALLOCATE( AD13_SO2_ac( IIPAR, JJPAR, LD13 ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD13_SO2_ac' )

         ALLOCATE( AD13_SO2_an( IIPAR, JJPAR, 2 ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD13_SO2_an' )

         ALLOCATE( AD13_SO2_bb( IIPAR, JJPAR ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD13_SO2_bb' )

         ALLOCATE( AD13_SO2_bf( IIPAR, JJPAR ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD13_SO2_bf' )

         ALLOCATE( AD13_SO2_ev( IIPAR, JJPAR, LD13 ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD13_SO2_ev' )

         ALLOCATE( AD13_SO2_nv( IIPAR, JJPAR, LD13 ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD13_SO2_nv' )

         ALLOCATE( AD13_SO4_an( IIPAR, JJPAR, 2 ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD13_SO4_an' )

         ALLOCATE( AD13_SO2_sh( IIPAR, JJPAR ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD13_SO4_sh' )

         ALLOCATE( AD13_NH3_an( IIPAR, JJPAR ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD13_NH3_an' )

         ALLOCATE( AD13_NH3_na( IIPAR, JJPAR ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD13_NH3_na' )

         ALLOCATE( AD13_NH3_bb( IIPAR, JJPAR ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD13_NH3_bb' )

         ALLOCATE( AD13_NH3_bf( IIPAR, JJPAR ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD13_NH3_bf' )
      ENDIF         

      !=================================================================
      ! ND14: Upward flux of from wet conv [kg/s] 
      !       --> uses CONVFLUP array (allocatable)
      !=================================================================
      IF ( ND14 > 0 ) THEN
         LD14 = MIN( ND14,      LLCONVM )
         NMAX = MIN( N_TRACERS, NNPAR   )

         ALLOCATE( CONVFLUP( IIPAR, JJPAR, LLCONVM, NMAX ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'CONVFLUP' )
      ENDIF

      !=================================================================
      ! ND15: Mass change from BL-mixing [kg/s] 
      !       --> uses TURBFLUP array (allocatable)
      !=================================================================
      IF ( ND15 > 0 ) THEN
         LD15 = MIN( ND15,      LLPAR )
         NMAX = MIN( N_TRACERS, NNPAR )

         ALLOCATE( TURBFLUP( IIPAR, JJPAR, LLPAR, NMAX ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'TURBFLUP' )
      ENDIF

      !=================================================================
      ! ND16: Fraction of grid box experiencing large-scale and
      !       convective precipitation --> uses AD16 array (allocatable)
      !=================================================================
      IF ( ND16 > 0 ) THEN
         LD16 = MIN( ND16, LLPAR )

         ! Store both LS and convective fractions
         ALLOCATE( AD16( IIPAR, JJPAR, LD16, 2 ), STAT=AS ) 
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD16' )

         ! Counter array for AD16
         ALLOCATE( CT16( IIPAR, JJPAR, LD16, 2 ), STAT=AS ) 
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'CT16' )
      ENDIF

      !=================================================================
      ! ND17: Fraction of tracer lost to rainout (in both large-scale
      !       and conv precipitation) --> uses AD17 array (allocatable)
      !=================================================================
      IF ( ND17 > 0 ) THEN
         LD17 = MIN( ND17, LLPAR )

         ! Get max # of soluble tracers for this simulation
         NMAX = GET_WETDEP_NMAX()

         ! Store both LS and convective rainout fractions
         ALLOCATE( AD17( IIPAR, JJPAR, LD17, NMAX, 2 ), STAT=AS ) 
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD17' )

         ! Counter array for AD17
         ALLOCATE( CT17( IIPAR, JJPAR, LD17, 2 ), STAT=AS ) 
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'CT17' )
      ENDIF

      !=================================================================
      ! ND18: Fraction of tracer lost to washout (in both large-scale
      !       and convective precipitation) --> uses AD18 array (alloc.)
      !=================================================================
      IF ( ND18 > 0 ) THEN
         LD18 = MIN( ND18, LLPAR )

         ! Get max # of soluble tracers for this simulation
         NMAX = GET_WETDEP_NMAX()

         ! Store both LS and convective rainout fractions         
         ALLOCATE( AD18( IIPAR, JJPAR, LD18, NMAX, 2 ), STAT=AS) 
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD18' )

         ! Counter array for AD17
         ALLOCATE( CT18( IIPAR, JJPAR, LD18, 2 ), STAT=AS ) 
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'CT18' )
      ENDIF

      !=================================================================
      ! ND19: Free diagnostic
      !
      ! ND20: Save O3 P-L losses to disk for single-tracer O3 run
      !       in the PL24H array.  Also turn on ND65, since the P-L 
      !       rates are computed by ND65.  
      !=================================================================
      IF ( ND20 > 0 ) THEN 
         IF ( ND65 == 0 ) ND65 = LLTROP
      ENDIF

      !=================================================================
      ! ND21: Optical depths and cloud fractions [unitless]
      !       --> uses AD21 array (allocatable) 
      !=================================================================
      IF ( ND21 > 0 ) THEN
         !-----------------------------
         ! Prior to 9/28/04:
         !LD21 = MIN( ND21, LLPAR )
         !-----------------------------
         LD21 = MIN( ND21, LLTROP )

         ALLOCATE( AD21( IIPAR, JJPAR, LD21, PD21 ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD21' )
      ENDIF

      !=================================================================
      ! ND22: J-value diagnostics [s^-1] 
      !       --> uses AD22 array (allocatable)
      !=================================================================
      IF ( ND22 > 0 ) THEN

         ! For full chemistry, we only consider boxes below the
         ! tropopause.  Cap LD22 at LLTROP (bmy, 4/2/01)
         IF ( ITS_A_FULLCHEM_SIM() ) THEN
            LD22 = MIN( ND22, LLTROP )
         ELSE
            LD22 = MIN( ND22, LLPAR  )
         ENDIF

         ! Accumulating diagnostic array
         ALLOCATE( AD22( IIPAR, JJPAR, LD22, PD22 ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD22' )

         ! Locations where LT is between HR1_JV and HR2_JV
         ALLOCATE( LTJV( IIPAR, JJPAR ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'LTJV' )         

         ! Number of times where LT is between HR1_JV and HR2_JV
         ALLOCATE( CTJV( IIPAR, JJPAR ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'CTJV' )         
      ENDIF

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
         LD24 = MIN( ND24,      LLPAR )
         NMAX = MIN( N_TRACERS, NNPAR )

         ALLOCATE( MASSFLEW( IIPAR, JJPAR, LLPAR, NMAX ), STAT=AS) 
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'MASSFLEW' )
      ENDIF

      !=================================================================
      ! ND25: Northward mass flux from transport [kg/s] 
      !       --> uses MASSFLNS array (allocatable)
      !=================================================================
      IF ( ND25 > 0 ) THEN
         LD25 = MIN( ND25,      LLPAR )
         NMAX = MIN( N_TRACERS, NNPAR )

         ALLOCATE( MASSFLNS( IIPAR, JJPAR, LLPAR, NMAX ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'MASSFLNS' ) 
      ENDIF

      !=================================================================
      ! ND26: Vertical mass flux from transport [kg/s] 
      !       --> uses MASSFLUP array (allocatable)
      !=================================================================
      IF ( ND26 > 0 ) THEN
         LD26 = MIN( ND26,      LLPAR )
         NMAX = MIN( N_TRACERS, NNPAR )

         ALLOCATE( MASSFLUP( IIPAR, JJPAR, LLPAR, NMAX ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'MASSFLUP' )
      ENDIF

      !=================================================================
      ! ND28: Biomass burning diagnostic [molec/cm2/s]
      !       (NO, CO, ALK4, ACET, MEK, ALD2, PRPE, C3H8, CH2O, C2H6) 
      !       --> uses AD28 array (allocatable)
      !=================================================================
      IF ( ND28 > 0 ) THEN
         ALLOCATE( AD28( IIPAR, JJPAR, NBIOTRCE ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD28' )
      ENDIF

      !=================================================================
      ! ND29: CO-SRCE diagnostic [molec/cm2/s]
      !       (anthro, biomass, biofuel, from monot., from methanol) 
      !        --> uses AD29 array (allocatable)
      !=================================================================
      IF ( ND29 > 0 ) THEN
         ALLOCATE( AD29( IIPAR, JJPAR, PD29 ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD29' )
      ENDIF

      !=================================================================
      ! ND30: Plot land map...does not use any allocatable arrays
      !
      ! ND31: PS - PTOP [mb] --> Uses AD31 array (allocatable)
      !=================================================================
      IF ( ND31 > 0 ) THEN
         ALLOCATE( AD31( IIPAR, JJPAR, PD31 ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD31' )
      ENDIF

      !=================================================================
      ! ND32: Sources of NOx [molec/cm2/s]
      !       (aircraft, biomass, biofuel, lightning, 
      !        stratosphere, soils, fertilizer, anthropogenic) 
      !       --> Uses AD32_xx arrays (allocatable)
      !=================================================================
      IF ( ND32 > 0 ) THEN

         ! For aircraft NOx
         ALLOCATE( AD32_ac( IIPAR, JJPAR, LLTROP ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD32_ac' )
         
         ! For anthropogenic NOx
         ALLOCATE( AD32_an( IIPAR, JJPAR, NOXEXTENT ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD32_an' ) 

         ! For biomass burning NOx
         ALLOCATE( AD32_bb( IIPAR, JJPAR ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD32_bb' )

         ! For biofuel NOx
         ALLOCATE( AD32_bf( IIPAR, JJPAR ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD32_bf' )

         ! For fertilizer NOx
         ALLOCATE( AD32_fe( IIPAR, JJPAR ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD32_fe' ) 

         ! For Lightning NOx
         !------------------------------------------------------
         ! Prior to 7/20/04:
         !ALLOCATE( AD32_li( IIPAR, JJPAR, LCONVM ), STAT=AS )
         !------------------------------------------------------
         ALLOCATE( AD32_li( IIPAR, JJPAR, LLCONVM ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD32_li' )

         ! For soil NOx
         ALLOCATE( AD32_so( IIPAR, JJPAR ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD32_so' ) 

         ! For stratospheric NOx
         ALLOCATE( AD32_ub( IIPAR, JJPAR ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD32_st' ) 
      ENDIF

      !=================================================================
      ! ND33: Column sum of tracer [kg]
      !       --> uses AD33 array (allocatable)!
      !=================================================================
      IF ( ND33 > 0 ) THEN
         ALLOCATE( AD33( IIPAR, JJPAR, PD33 ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD33' )
      ENDIF

      !=================================================================
      ! ND34: Biofuel burning emissions [molec/cm2/s]
      !       (NO, CO, ALK4, ACET, MEK, ALD2, PRPE, C3H8, CH2O, C2H6)  
      !       --> uses AD34 array (allocatable)
      !=================================================================
      IF ( ND34 > 0 ) THEN
         ALLOCATE( AD34( IIPAR, JJPAR, NBFTRACE ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD34' )
      ENDIF         

      !=================================================================
      ! ND35: Tracer at 500 mb [v/v] (this is ~ level 9 for GEOS-CHEM)
      !       --> uses AD35 array (allocatable)
      !=================================================================
      IF ( ND35 > 0 ) THEN
         ALLOCATE( AD35( IIPAR, JJPAR, N_TRACERS ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD35' ) 
      ENDIF

      !=================================================================
      ! ND36: Anthropogenic Emisisons [molec/cm2/s]
      !       (NOx, CO, ALK4, ACET, MEK, ALD2, PRPE, C3H8, C2H6)
      !       --> uses AD36 array (allocatable)
      !  
      ! NOTE: For a CH3I run, use ND36 for CH3I emission diagnostics....
      !=================================================================
      IF ( ITS_A_CH3I_SIM() ) NEMANTHRO = 8   ! for CH3I

      IF ( ND36 > 0 ) THEN
         ALLOCATE( AD36( IIPAR, JJPAR, NEMANTHRO ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD36' )
      ENDIF

      !=================================================================
      ! ND37: Fraction of tracer scavenged in cloud updrafts 
      !       --> Uses AD37 array (allocatable)
      !=================================================================
      IF ( ND37 > 0 ) THEN
         LD37 = MIN( ND37, LLPAR )
         
         ! Get max # of soluble tracers for this simulation
         NMAX = GET_WETDEP_NMAX()

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
      IF ( ND38 > 0 ) THEN
         LD38 = MIN( ND38, LLPAR )
         
         ! Get max # of soluble tracers for this simulation
         NMAX = GET_WETDEP_NMAX()

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
      IF ( ND39 > 0 ) THEN
         LD39 = MIN( ND39, LLPAR )

         ! Get max # of soluble tracers for this simulation
         NMAX = GET_WETDEP_NMAX()
         
         ! Allocate AD39 array accordingly
         IF ( NMAX > 0 ) THEN
            ALLOCATE( AD39( IIPAR, JJPAR, LD39, NMAX ), STAT=AS )
            IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD39' )
         ENDIF
      ENDIF

      !=================================================================
      ! ND40: Plane flight track diagnostic
      !
      ! NOTE: For SMVGEAR chemistry, we first have to read "chem.dat"
      !       before initializing the plane flight diagnostic.  Call
      !       routine SETUP_PLANEFLIGHT from "chemdr.f" in that case.
      !=================================================================
      IF ( ND40 > 0 .and. ( .not. ITS_A_FULLCHEM_SIM() ) ) THEN
         CALL SETUP_PLANEFLIGHT
      ENDIF

      !=================================================================
      ! ND41: PM (1200-1600 LT) boundary layer height [layers]
      !       --> uses AD41 array (allocatable)
      !=================================================================
      IF ( ND41 > 0 ) THEN

         ! Accumulating diagnostic array
         ALLOCATE( AD41( IIPAR, JJPAR, PD41 ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD41' )

         ! Locations where local time is between 1 and 4 PM 
         ALLOCATE( AFTTOT( IIPAR, JJPAR ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'AFTTOT' )

         ! SAVENO2 is not used for now...
         !ALLOCATE( SAVENO2( IIPAR, JJPAR, PD41 ), STAT=AS )
         !IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD41' )
      ENDIF

      !=================================================================
      ! ND42: Free diagnostic
      !
      ! ND43: Chemical diagnostics: OH [molec/cm3/s] and NO [v/v]
      !       --> uses AD43 array (allocatable)
      !=================================================================
      IF ( ND43 > 0 ) THEN
         
         ! For full chemistry, only save OH up to the tropopause.  
         ! For tagged CO or CO-OH, save OH everywhere (bmy, 2/12/01)
         IF ( ITS_A_FULLCHEM_SIM() ) THEN
            LD43 = MIN( ND43, LLTROP )
         ELSE
            LD43 = MIN( ND43, LLPAR  )
         ENDIF

         ! Accumulating diagnostic array
         ALLOCATE( AD43( IIPAR, JJPAR, LD43, PD43 ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD43' )

         ! Locations where LT is between HR1_NO and HR2_NO
         ALLOCATE( LTNO( IIPAR, JJPAR ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'LTNO' )
         
         ! Number of times LT was between HR1_NO and HR2_NO
         ALLOCATE( CTNO( IIPAR, JJPAR ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'CTNO' )

         ! Locations where LT is between HR1_OH and HR2_OH
         ALLOCATE( LTOH( IIPAR, JJPAR ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'LTOH' )

         ! Locations where LT is between HR1_OH and HR2_OH
         ALLOCATE( CTOH( IIPAR, JJPAR ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'CTOH' )

         ! Locations where LT is between HR1_OH and HR2_OH
         ALLOCATE( LTHO2( IIPAR, JJPAR ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'LTHO2' )

         ! Locations where LT is between HR1_OH and HR2_OH
         ALLOCATE( CTHO2( IIPAR, JJPAR ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'CTHO2' )

         ! Locations where LT is between HR1_OH and HR2_OH
         ALLOCATE( LTNO2( IIPAR, JJPAR ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'LTHO2' )

         ! Locations where LT is between HR1_OH and HR2_OH
         ALLOCATE( CTNO2( IIPAR, JJPAR ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'CTHO2' )

         ! Locations where LT is between HR1_OH and HR2_OH
         ALLOCATE( LTNO3( IIPAR, JJPAR ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'LTHO2' )

         ! Locations where LT is between HR1_OH and HR2_OH
         ALLOCATE( CTNO3( IIPAR, JJPAR ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'CTHO2' )
      ENDIF

      !=================================================================
      ! ND44: Drydep fluxes [s-1] and drydep velocities [cm/s]
      !       --> uses AD44 arrays (allocatable)
      !=================================================================
      IF ( ND44 > 0 ) THEN

         ! Get number of tracers for ND44
         IF ( ITS_A_TAGOX_SIM() ) THEN
            NMAX = N_TRACERS
         ELSE
            NMAX = NUMDEP
         ENDIF

         ! Allocate AD44 array
         ALLOCATE( AD44( IIPAR, JJPAR, NMAX, 2 ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD44' )
      ENDIF

      !=================================================================
      ! ND45: Tracer concentrations [v/v] between HR1_OTH and HR2_OTH
      !       --> uses AD45 array (allocatable)
      !=================================================================
      IF ( ND45 > 0 ) THEN
         LD45 = MIN( ND45,        LLPAR   )
         NMAX = MIN( N_TRACERS+1, NNPAR+1 )

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
      ! ND46: Biogenic emissions [molec/cm2/s]
      !       (ISOP, PRPE, ACET, and MONOTERPENES)
      !       --> uses AD46 array (allocatable)
      !=================================================================
      IF ( ND46 > 0 ) THEN
         ALLOCATE( AD46( IIPAR, JJPAR, PD46 ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD46' ) 
      ENDIF

      !=================================================================
      ! ND47: 24-h averaged tracer concentration [v/v]
      !       --> uses AD47 array (allocatable)
      !
      ! NOTE: ND47 is always a 24-h average field, while ND45 
      !       can be averaged over any arbitrary time period.
      !=================================================================
      IF ( ND47 > 0 ) THEN
         LD47 = MIN( ND47,        LLPAR   )
         NMAX = MIN( N_TRACERS+1, NNPAR+1 )

         ! Resize to NMAX to save memory (bmy, 10/18/00)
         ALLOCATE( AD47( IIPAR, JJPAR, LD47, NMAX ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD47' )
      ENDIF

      !=================================================================
      ! ND48: Station timeseries diagnostics (also see diag48.f)
      !       --> uses TCOBOX array (allocatable)
      !=================================================================
      IF ( ND48 > 0 ) THEN
         ALLOCATE( TCOBOX( MAXACC, NNSTA ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'TCOBOX' )
      ENDIF        

!-----------------------------------------------------------------------------
! Prior to 7/20/04:
!      !=================================================================
!      ! ND49: Geographic domain -- inst. time series of tracers
!      ! ND50: Geographic domain -- 24h average time series of tracers
!      ! ND51: Geographic domain -- 1-4pm time series for tracers
!      !
!      ! Call READ49 to read "timeseries.dat" for ND49, ND50, ND51
!      !=================================================================
!      !--------------------------------------------------------
!      ! Prior to 4/22/04: 
!      ! Now activate for ND52 diagnostic (bmy, 4/22/04)
!      !IF ( ND49 > 0 .or. ND50 .gt. 0 .or. ND51 > 0 ) THEN
!      !--------------------------------------------------------
!      IF ( ND49 > 0 .or. ND50 > 0 .or. ND51 > 0 .or. ND52 > 0 ) THEN
!         CALL READ49
!
!         ! For ND50 -- allocate space for the STT_TEMPO2 array. 
!         ! Use NTRACE tracers, and pure O3 and NO if necessary (bmy, 10/13/00)
!         IF ( ND50 > 0 ) THEN
!            ALLOCATE( STT_TEMPO2( IIPAR, JJPAR, LLPAR, NTRACE+2 ), 
!     &                            STAT=AS )
!
!            IF ( AS /= 0 ) CALL ALLOC_ERR( 'STT_TEMPO2' )
!         ENDIF
!      ENDIF
!-----------------------------------------------------------------------------

      !=================================================================
      ! ND52 - ND54: Free diagnostics
      ! 
      ! ND55: Tropopause diagnostics [level, height, and pressure]
      !       --> uses AD55 array (allocatable)
      !=================================================================
      IF ( ND55 > 0 ) THEN
         ALLOCATE( AD55( IIPAR, JJPAR, PD55 ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD55' )
      ENDIF

      !=================================================================
      ! ND56 - ND64: Free diagnostics
      !
!-----------------------------------------------------------------------------
! Prior to 7/20/04:
! We now do this in "diag_pl_mod.f"
!      ! ND65: Chemical family production and loss [molec/cm3/s]
!      !       Call setnfam to get # of families to keep track of (nfam)
!      !=================================================================
!      NFAM = 0
!
!
!      IF ( ND65 > 0 ) THEN
!
!         ! Select NFAM by simulation
!         SELECT CASE ( NSRCX )
!
!            ! CH3I: NTRACE families 
!            CASE ( 2 )
!               NFAM = NTRACE
!
!            ! Tagged Ox: NTRACE*2 families
!            CASE ( 6 )
!               NFAM = NTRACE * 2
!
!            ! Tagged CO: NTRACE+5 families
!            CASE ( 7 )
!               NFAM = NTRACE + 5
!
!            ! Other simulations: call SETNFAM
!            CASE DEFAULT
!               CALL SETNFAM
!
!         END SELECT
!
!         ! Error check: NFAM must be > 0! (bmy, 1/15/02)
!         IF ( NFAM == 0 ) THEN 
!            CALL ERROR_STOP( 'NFAM must be >0 for ND65', 'ndxx_setup.f')
!         ENDIF
!
!         ! Select by simulation
!         SELECT CASE ( NSRCX )
!
!            !----------------------------------------------
!            ! For full chem w/ SMVGEAR: save up to LLTROP
!            !----------------------------------------------
!            CASE ( 3 )
!
!               ! Number of vertical levels to print out
!               LD65 = MIN( ND65, LLTROP )
!
!               ! Allocate AD65
!               ALLOCATE( AD65( IIPAR, JJPAR, LD65, NFAM ), STAT=AS )
!               IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD65' )
!
!               ! Allocate FAMPL for "plsave.f", "diag20.f", etc.
!               ALLOCATE( FAMPL( IIPAR, JJPAR, LLTROP, NFAM ), STAT=AS )
!               IF ( AS /= 0 ) CALL ALLOC_ERR( 'FAMPL' )
!
!            !----------------------------------------------
!            ! All other simulations: save up to LLPAR
!            !----------------------------------------------
!            CASE DEFAULT
!
!               ! Number of vertical levels to print out
!               LD65 = MIN( ND65, LLPAR )
!
!               ! Allocate AD65
!               ALLOCATE( AD65( IIPAR, JJPAR, LD65, NFAM ), STAT=AS )
!               IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD65' )
!
!          END SELECT
!
!      ENDIF
!-----------------------------------------------------------------------------

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

!------------------------------------------------------------------------------
! Prior to 7/20/04:
!      !=================================================================
!      ! Write NDxx's to standard output
!      !=================================================================
! 975  FORMAT(' ND', I1, '-', 10I4)
!      WRITE(6,975) 0,ND01,ND02,ND03,ND04,ND05,ND06,ND07,ND08,ND09,ND10
!      WRITE(6,975) 1,ND11,ND12,ND13,ND14,ND15,ND16,ND17,ND18,ND19,ND20
!      WRITE(6,975) 2,ND21,ND22,ND23,ND24,ND25,ND26,ND27,ND28,ND29,ND30
!      WRITE(6,975) 3,ND31,ND32,ND33,ND34,ND35,ND36,ND37,ND38,ND39,ND40
!      WRITE(6,975) 4,ND41,ND42,ND43,ND44,ND45,ND46,ND47,ND48,ND49,ND50
!      WRITE(6,975) 5,ND51,ND52,ND53,ND54,ND55,ND56,ND57,ND58,ND59,ND60
!      WRITE(6,975) 6,ND61,ND62,ND63,ND64,ND65,ND66,ND67,ND68,ND69,ND70
!------------------------------------------------------------------------------

      ! Return to calling program
      END SUBROUTINE NDXX_SETUP
