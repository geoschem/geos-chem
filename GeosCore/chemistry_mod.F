!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: chemistry_mod
!
! !DESCRIPTION: Module CHEMISTRY\_MOD is used to call the proper chemistry 
!  subroutine for the various GEOS-Chem simulations. 
!\\
!\\
! !INTERFACE:
!
      MODULE CHEMISTRY_MOD
!
! !USES:
!
      IMPLICIT NONE
#     include "define.h"
      PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
      PUBLIC :: DO_CHEMISTRY
      PUBLIC :: GCKPP_DRIVER
      PUBLIC :: RECOMPUTE_OD
!
! !REVISION HISTORY: 
!  (1 ) Bug fix in DO_CHEMISTRY (bnd, bmy, 4/14/03)
!  (2 ) Now references DEBUG_MSG from "error_mod.f" (bmy, 8/7/03)
!  (3 ) Now references "tagged_ox_mod.f"(bmy, 8/18/03)
!  (4 ) Now references "Kr85_mod.f" (jsw, bmy, 8/20/03)
!  (5 ) Bug fix: Now also call OPTDEPTH for GEOS-4 (bmy, 1/27/04)
!  (6 ) Now references "carbon_mod.f" and "dust_mod.f" (rjp, tdf, bmy, 4/5/04)
!  (7 ) Now references "seasalt_mod.f" (rjp, bec, bmy, 4/20/04)
!  (8 ) Now references "logical_mod.f", "tracer_mod.f", "diag20_mod.f", and
!        "diag65_mod.f", and "aerosol_mod." (bmy, 7/20/04)
!  (9 ) Now references "mercury_mod.f" (bmy, 12/7/04)
!  (10) Updated for SO4s, NITs chemistry (bec, bmy, 4/13/05)
!  (11) Now call CHEM_HCN_CH3CN from "hcn_ch3cn_mod.f".  Also remove all
!        references to the obsolete CO-OH param simulation. (xyp, bmy, 6/24/05)
!  (12) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (13) Now call MAKE_RH from "main.f" (bmy, 3/16/06)
!  (14) Updated for SOA from isoprene (dkh, bmy, 6/1/06)
!  (15) Remove support for GEOS-1 and GEOS-STRAT met fields (bmy, 8/4/06)
!  (16) For now, replace use RPMARES instead of ISORROPIA. (bmy, 4/2/08)
!  (17) Added KPP chemistry driver subroutine (phs,ks,dhk, 09/15/09)
!  (18) Added public member function recompute_OD (skim, 02/03/11)
!  17 Dec 2009 - R. Yantosca - Added ProTeX headers
!  28 Jan 2010 - C. Carouge, R. Yantosca - Modified for ISORROPIA II
!  08 Aug 2012 - R. Yantosca - Now align IF statements better
!  10 Aug 2012 - R. Yantosca - Cosmetic changes
!EOP
!------------------------------------------------------------------------------
!BOC
      CONTAINS
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: do_chemistry
!
! !DESCRIPTION: Subroutine DO\_CHEMISTRY is the driver routine which calls 
!  the appropriate chemistry subroutine for the various GEOS-Chem simulations.
!\\
!\\
! !INTERFACE:
!
#if defined( DEVEL ) || defined( EXTERNAL_GRID ) || defined( EXTERNAL_FORCING )
      !-----------------------------------------------------------------
      !     %%%%% CONNECTING TO GEOS-5 GCM via ESMF INTERFACE %%%%%
      !
      ! Call CHEMDR with dimension values as well as the Input Options,
      ! Meteorology State,  and Chemistry State objects. (bmy, 11/8/12)
      !-----------------------------------------------------------------
      SUBROUTINE DO_CHEMISTRY( am_I_Root, NI,        NJ,        NL, 
     &                         Input_Opt, State_Chm, State_Met, RC )
#else
      !-----------------------------------------------------------------
      !              %%%%% TRADITIONAL GEOS-Chem %%%%%
      !
      ! Current practice is to call DO_CHEMISTRY with am_I_Root as the
      ! only argument. (bmy, 10/25/12)
      !-----------------------------------------------------------------
      SUBROUTINE DO_CHEMISTRY( am_I_Root )
#endif
!
! !USES:
!
      USE ACETONE_MOD,     ONLY : OCEAN_SINK_ACET
      USE AEROSOL_MOD,     ONLY : AEROSOL_CONC, AEROSOL_RURALBOX
      USE AEROSOL_MOD,     ONLY : RDAER,        SOILDUST
      USE C2H6_MOD,        ONLY : CHEMC2H6
      USE CARBON_MOD,      ONLY : CHEMCARBON
      USE CH3I_MOD,        ONLY : CHEMCH3I
      USE DAO_MOD,         ONLY : CLDF,    DELP
      USE DAO_MOD,         ONLY : OPTDEP,  OPTD,   T
      USE DRYDEP_MOD,      ONLY : DRYFLX, DRYFLXRnPbBe, DRYFLXH2HD
      USE DUST_MOD,        ONLY : CHEMDUST, RDUST_ONLINE
      USE ERROR_MOD,       ONLY : DEBUG_MSG
      USE GLOBAL_CH4_MOD,  ONLY : CHEMCH4
      USE H2_HD_MOD,       ONLY : CHEM_H2_HD
      USE HCN_CH3CN_MOD,   ONLY : CHEM_HCN_CH3CN
      USE ISOROPIAII_MOD,  ONLY : DO_ISOROPIAII
      USE MERCURY_MOD,     ONLY : CHEMMERCURY
      USE OPTDEPTH_MOD,    ONLY : OPTDEPTH
      USE RnPbBe_MOD,      ONLY : CHEMRnPbBe
      USE RPMARES_MOD,     ONLY : DO_RPMARES
      USE SEASALT_MOD,     ONLY : CHEMSEASALT
      USE SULFATE_MOD,     ONLY : CHEMSULFATE
      USE STRAT_CHEM_MOD,  ONLY : DO_STRAT_CHEM
      USE TAGGED_CO_MOD,   ONLY : CHEM_TAGGED_CO
      USE TAGGED_OX_MOD,   ONLY : CHEM_TAGGED_OX
      USE TIME_MOD,        ONLY : GET_ELAPSED_MIN, GET_TS_CHEM
      USE TRACER_MOD,      ONLY : STT  
      USE TRACERID_MOD,    ONLY : IDTACET, IDTISOP

      USE CMN_SIZE_MOD          ! Size parameters
      USE COMODE_LOOP_MOD       ! NPHOT
      USE CMN_DIAG_MOD          ! NDxx flags

#if defined( DEVEL ) || defined( EXTERNAL_GRID ) || defined( EXTERNAL_FORCING )
      !-----------------------------------------------------------------
      !     %%%%% CONNECTING TO GEOS-5 GCM via ESMF INTERFACE %%%%%
      !
      ! We need to reference some modules particular to the Grid-
      ! Independent GEOS-Chem implementation here.   In particular,
      ! the settings from the Input Options object will replace
      ! logical flags from logical_mod.F. (bmy, 11/8/12)
      !-----------------------------------------------------------------
      USE GIGC_ErrCode_Mod
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
      USE GIGC_State_Chm_Mod, ONLY : ChmState
      USE GIGC_State_Met_Mod, ONLY : MetState
#else
      !-----------------------------------------------------------------
      !              %%%%% TRADITIONAL GEOS-Chem %%%%%
      !
      ! Current practice in the standard GEOS-Chem is to reference
      ! logical flags directly from module GeosCore/logical_mod.F.
      ! We also reference the ITS_A_*_SIM functions directly from
      ! module GeosCore/tracer_mod.F.   (bmy, 11/8/12)
      !-----------------------------------------------------------------
      USE LOGICAL_MOD,     ONLY : LCARB, LCHEM,  LCRYST, LDUST, LSCHEM
      USE LOGICAL_MOD,     ONLY : LPRT,  LSSALT, LSULF,  LSOA,  LNLPBL 
      USE TRACER_MOD,      ONLY : N_TRACERS
      USE TRACER_MOD,      ONLY : ITS_A_C2H6_SIM
      USE TRACER_MOD,      ONLY : ITS_A_CH3I_SIM
      USE TRACER_MOD,      ONLY : ITS_A_CH4_SIM
      USE TRACER_MOD,      ONLY : ITS_A_FULLCHEM_SIM
      USE TRACER_MOD,      ONLY : ITS_A_H2HD_SIM
      USE TRACER_MOD,      ONLY : ITS_A_HCN_SIM
      USE TRACER_MOD,      ONLY : ITS_A_MERCURY_SIM
      USE TRACER_MOD,      ONLY : ITS_A_RnPbBe_SIM
      USE TRACER_MOD,      ONLY : ITS_A_TAGCO_SIM
      USE TRACER_MOD,      ONLY : ITS_A_TAGOX_SIM
      USE TRACER_MOD,      ONLY : ITS_AN_AEROSOL_SIM
      USE TRACER_MOD,      ONLY : ITS_NOT_COPARAM_OR_CH4
#endif
!
! !INPUT PARAMETERS:
!
      LOGICAL,        INTENT(IN)    :: am_I_Root   ! Is this the root CPU?
#if defined( DEVEL ) || defined( EXTERNAL_GRID ) || defined( EXTERNAL_FORCING )
      INTEGER,        INTENT(IN)    :: NI          ! # of lons on this PET
      INTEGER,        INTENT(IN)    :: NJ          ! # of lats on this PET 
      INTEGER,        INTENT(IN)    :: NL          ! # of levels on this PET
      TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
      TYPE(MetState), INTENT(IN)    :: State_Met   ! Meteorology State object
!
! !INPUT/OUTPUT PARAMETERS:
!
      TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
!
! !OUTPUT PARAMETERS:
!
      INTEGER,        INTENT(OUT)   :: RC          ! Success or failure
#endif
!
! !REMARKS:
!  (1) State_Chm%TRACERS needs to enter this routine w/ units of [kg].
!  (2) State_Chm%SPECIES needs to enter this routine w/ units of [molec/cm3].
!  (3) As of 25 Oct 2012, we have not replaced the met-field arrays from
!       dao_mod.F with the fields in State_Met.  This will be done in a
!       piecemeal fashion so as not to disrupt the existing G-C workflow.
!
! !REVISION HISTORY: 
!  (1 ) Now reference DELP, T from "dao_mod.f" since we need to pass this
!        to OPTDEPTH for GEOS-1 or GEOS-STRAT met fields (bnd, bmy, 4/14/03)
!  (2 ) Now references DEBUG_MSG from "error_mod.f" (bmy, 8/7/03)
!  (3 ) Removed call to CHEMO3, it's obsolete.  Now calls CHEM_TAGGED_OX !
!        from "tagged_ox_mod.f" when NSRCX==6.  Now calls Kr85 chemistry if 
!        NSRCX == 12 (jsw, bmy, 8/20/03)
!  (4 ) Bug fix: added GEOS-4 to the #if block in the call to OPTDEPTH.
!        (bmy, 1/27/04)
!  (5 ) Now calls CHEMCARBON and CHEMDUST to do carbon aerosol & dust 
!        aerosol chemistry (rjp, tdf, bmy, 4/2/04)
!  (6 ) Now calls CHEMSEASALT to do seasalt aerosol chemistry 
!        (rjp, bec, bmy, 4/20/04)
!  (7 ) Now references "logical_mod.f" & "tracer_mod.f".  Now references
!        AEROSOL_CONC, AEROSOL_RURALBOX, and RDAER from "aerosol_mod.f".  
!        Now includes "CMN_DIAG" and "comode.h".  Also call READER, READCHEM, 
!        and INPHOT to initialize the FAST-J arrays so that we can save out !
!        AOD's to the ND21 diagnostic for offline runs. (bmy, 7/20/04)
!  (8 ) Now call routine CHEMMERCURY from "mercury_mod.f" for an offline
!        Hg0/Hg2/HgP simulation. (eck, bmy, 12/7/04)
!  (9 ) Now do not call DO_RPMARES if we are doing an offline aerosol run
!        with crystalline sulfur & aqueous tracers (cas, bmy, 1/7/05)
!  (10) Now use ISOROPIA for aer thermodyn equilibrium if we have seasalt 
!        tracers defined, or RPMARES if not.  Now call CHEMSEASALT before
!        CHEMSULFATE.  Now do aerosol thermodynamic equilibrium before
!        aerosol chemistry for offline aerosol runs.  Now also reference 
!        CLDF from "dao_mod.f" (bec, bmy, 4/20/05)
!  (11) Now modified for GCAP met fields.  Now call CHEM_HCN_CH3CN from 
!        "hcn_ch3cn_mod.f".  Also remove allreferences to the obsolete 
!         CO-OH param simulation. (xyp, bmy, 6/23/05)
!  (12) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (13) Now call MAKE_RH from "main.f" (bmy, 3/16/06)
!  (14) Removed ISOP_PRIOR as a local variable (dkh, bmy, 6/1/06)
!  (15) Remove support for GEOS-1 and GEOS-STRAT met fields (bmy, 8/4/06)
!  (16) Now use DRYFLXH2HD and CHEM_H2_HD for H2/HD sim (lyj, phs, 9/18/07)
!  (17) Bug fix: now hardwired to use RPMARES since ISORROPIA can return very
!        unphysical values at low RH.  Wait for ISORROPIA II. (bmy, 4/2/08)
!  (18) The dry deposition diagnostic (ND44) is done in vdiff_mod if using non-
!        local PBL (lin, ccc, 5/29/09)
!  17 Dec 2009 - R. Yantosca - Added ProTeX headers
!  28 Jan 2010 - C. Carouge, R. Yantosca - Modified for ISORROPIA II
!  19 Mar 2012 - R. Yantosca - Add C-preprocessor switch to shut off 
!                              ISORROPIA to facilitate debugging
!  30 Jul 2012 - R. Yantosca - Now accept am_I_Root as an argument, and pass
!                              this down to lower-level chem routines for GIGC
!  08 Aug 2012 - R. Yantosca - Now align IF statements better
!  10 Aug 2012 - R. Yantosca - Cosmetic changes
!  18 Oct 2012 - R. Yantosca - Rename GC_MET argument to State_Met
!  18 Oct 2012 - R. Yantosca - Rename CHEM_STATE argument to State_Chem
!  19 Oct 2012 - R. Yantosca - Now reference gigc_state_chm_mod.F90
!  19 Oct 2012 - R. Yantosca - Now reference gigc_state_met_mod.F90
!  25 Oct 2012 - R. Yantosca - Add comments for GIGC #ifdefs
!  25 Oct 2012 - R. Yantosca - Add the RC output argument for the GIGC
!  08 Nov 2012 - R. Yantosca - Now pass Input_Opt argument for the GIGC and
!                              use fields of Input_Opt to replace logicals
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
! 
      INTEGER       :: N_TROP
      INTEGER       :: MONTH
      INTEGER       :: YEAR
      INTEGER       :: WAVELENGTH
      LOGICAL       :: IT_IS_A_C2H6_SIM
      LOGICAL       :: IT_IS_A_CH3I_SIM
      LOGICAL       :: IT_IS_A_CH4_SIM
      LOGICAL       :: IT_IS_A_FULLCHEM_SIM
      LOGICAL       :: IT_IS_A_H2HD_SIM
      LOGICAL       :: IT_IS_A_HCN_SIM
      LOGICAL       :: IT_IS_A_MERCURY_SIM
      LOGICAL       :: IT_IS_A_RnPbBe_SIM
      LOGICAL       :: IT_IS_A_TAGCO_SIM
      LOGICAL       :: IT_IS_A_TAGOX_SIM
      LOGICAL       :: IT_IS_AN_AEROSOL_SIM
      LOGICAL       :: IT_IS_NOT_COPARAM_OR_CH4

#if defined( DEVEL ) || defined( EXTERNAL_GRID ) || defined( EXTERNAL_FORCING )
      !-----------------------------------------------------------------
      !     %%%%% CONNECTING TO GEOS-5 GCM via ESMF INTERFACE %%%%%
      !
      ! We need to define local logical flags to hold corresponding
      ! values from the Input Options object.  We do not need to define
      ! these for the traditional GEOS-Chem since it will refer to 
      ! variables from the same name from GeosCore/logical_mod.F.
      ! Also need to shadow the N_TRACERS variable from module
      ! GeosCore/tracer_mod.F. (bmy, 11/9/12)
      !-----------------------------------------------------------------
      INTEGER       :: N_TRACERS
      LOGICAL       :: LCARB
      LOGICAL       :: LCHEM
      LOGICAL       :: LCRYST
      LOGICAL       :: LDUST
      LOGICAL       :: LSCHEM
      LOGICAL       :: LPRT
      LOGICAL       :: LSSALT
      LOGICAL       :: LSULF
      LOGICAL       :: LSOA
      LOGICAL       :: LNLPBL
#endif

      ! SAVEd variables
      LOGICAL, SAVE :: FIRST = .TRUE.

      !=================================================================
      ! DO_CHEMISTRY begins here!
      !=================================================================

#if defined( DEVEL ) || defined( EXTERNAL_GRID ) || defined( EXTERNAL_FORCING )
      !-----------------------------------------------------------------
      !     %%%%% CONNECTING TO GEOS-5 GCM via ESMF INTERFACE %%%%%
      !
      ! At present (25 Oct 2012), we have not yet completely replaced 
      ! the GEOS-Chem arrays with fields from the Chemistry State
      ! (i.e. State_Chm) and Meteorology State (i.e. State_Met) 
      ! objects.  In order to continue GIGC development without 
      ! disrupting the existing GEOS-Chem workflow, we shall copy the 
      ! tracer array State_Chm%TRACERS to STT before chemistry. 
      !
      ! Also note: State_Chm%TRACERS needs to enter this routine with
      ! units of [kg].  If we are calling this routine from the 
      ! traditional "main.F" driver, then the tracers are already in
      ! [kg] and no further unit conversion is necessary.  If on the
      ! other hand we are calling this routine from the ESMF interface,
      ! then we need to convert State_Chm%TRACERS from [v/v] to [kg]
      ! in the calling routine GIGC_Do_Chem (in ESMF/gigc_chemdr.F90)
      !
      ! To simplify the coding for both the ESMF interface and the
      ! traditional GEOS-Chem, define local logical variables for
      ! various flags and simulation types. (bmy, 11/8/12)
      !-----------------------------------------------------------------

      ! Assume success
      RC                       = GIGC_SUCCESS

      ! Copy logical fields from INPUT_OPT to local variables for use below
      N_TRACERS                = Input_Opt%N_TRACERS
      LCARB                    = Input_Opt%LCARB                        
      LCHEM                    = Input_Opt%LCHEM
      LCRYST                   = Input_Opt%LCRYST
      LDUST                    = Input_Opt%LDUST
      LSCHEM                   = Input_Opt%LSCHEM
      LPRT                     = Input_Opt%LPRT
      LSSALT                   = Input_Opt%LSSALT
      LSULF                    = Input_Opt%LSULF
      LSOA                     = Input_Opt%LSOA
      LNLPBL                   = Input_Opt%LNLPBL
      IT_IS_A_C2H6_SIM         = Input_Opt%ITS_A_C2H6_SIM
      IT_IS_A_CH3I_SIM         = Input_Opt%ITS_A_CH3I_SIM
      IT_IS_A_CH4_SIM          = Input_Opt%ITS_A_CH4_SIM 
      IT_IS_A_FULLCHEM_SIM     = Input_Opt%ITS_A_FULLCHEM_SIM
      IT_IS_A_H2HD_SIM         = Input_Opt%ITS_A_H2HD_SIM
      IT_IS_A_HCN_SIM          = Input_Opt%ITS_A_HCN_SIM
      IT_IS_A_MERCURY_SIM      = Input_Opt%ITS_A_MERCURY_SIM
      IT_IS_A_RnPbBe_SIM       = Input_Opt%ITS_A_TAGCO_SIM
      IT_IS_A_TAGCO_SIM        = Input_Opt%ITS_A_TAGCO_SIM
      IT_IS_A_TAGOX_SIM        = Input_Opt%ITS_A_TAGOX_SIM
      IT_IS_AN_AEROSOL_SIM     = Input_Opt%ITS_AN_AEROSOL_SIM
      IT_IS_NOT_COPARAM_OR_CH4 = Input_Opt%ITS_NOT_COPARAM_OR_CH4

      ! Initialize GEOS-Chem tracer array [kg] from Chemistry State object
      STT(:,:,:,1:N_TRACERS)   = State_Chm%TRACERS(:,:,:,1:N_TRACERS)

#else
      !--------------------------------------------------------------
      !               %%%%% TRADITIONAL GEOS-Chem %%%%%
      !
      ! Define local logical flags to store values from the various
      ! ITS_A_*_SIM routines from GeosCore/tracer_mod.F.  We do this
      ! in order to clarify the coding for both the ESMF interface
      ! as well as for the traditional GEOS-Chem.  (bmy, 11/8/12)
      !--------------------------------------------------------------
      IT_IS_A_C2H6_SIM         = ITS_A_C2H6_SIM()
      IT_IS_A_CH3I_SIM         = ITS_A_CH3I_SIM()
      IT_IS_A_CH4_SIM          = ITS_A_CH4_SIM()
      IT_IS_A_FULLCHEM_SIM     = ITS_A_FULLCHEM_SIM()
      IT_IS_A_H2HD_SIM         = ITS_A_H2HD_SIM()
      IT_IS_A_HCN_SIM          = ITS_A_HCN_SIM()
      IT_IS_A_MERCURY_SIM      = ITS_A_MERCURY_SIM()
      IT_IS_A_RnPbBe_SIM       = ITS_A_RnPbBe_SIM()
      IT_IS_A_TAGCO_SIM        = ITS_A_TAGCO_SIM()
      IT_IS_A_TAGOX_SIM        = ITS_A_TAGOX_SIM()
      IT_IS_AN_AEROSOL_SIM     = ITS_AN_AEROSOL_SIM()
      IT_IS_NOT_COPARAM_OR_CH4 = ITS_NOT_COPARAM_OR_CH4()
#endif

      IF ( .not. IT_IS_A_CH4_SIM ) THEN
!#if defined( DEVEL ) || defined( EXTERNAL_GRID ) || defined( EXTERNAL_FORCING )
!         !--------------------------------------------------------------
!         !   %%%%% CONNECTING TO GEOS-5 GCM via ESMF INTERFACE %%%%%
!         !
!         ! Now pass size dimensions, the Input Options object, and the
!         ! Meteorology State object to OPTDEPTH.  NOTE: The OPTD
!         ! and CLDF fields of the State_Met object have already
!         ! been populated in main.F. (bmy, 11/9/12)
!         !--------------------------------------------------------------
!
!         ! Compute optical depth diagnostics
!         CALL OPTDEPTH( am_I_Root, NI,        NJ, NL,
!     &                  Input_Opt, State_Met, RC      )
!#else
         !--------------------------------------------------------------
         !              %%%%% TRADITIONAL GEOS-Chem %%%%%
         !
         ! Current practice in the standard GEOS-Chem is to pass
         ! met field arrays as arguments to OPTDEPTH. (bmy, 11/9/12)
         !--------------------------------------------------------------

         ! Compute visible optical depths
         CALL OPTDEPTH( LLPAR, CLDF, OPTDEP, OPTD )
!#endif
      ENDIF

      !=================================================================
      ! If LCHEM=T then call the chemistry subroutines
      !=================================================================
      IF ( LCHEM ) THEN

         !---------------------------------
         ! NOx-Ox-HC (w/ or w/o aerosols) 
         !---------------------------------
         IF ( IT_IS_A_FULLCHEM_SIM ) THEN 

#if defined( DEVEL ) || defined( EXTERNAL_GRID ) || defined( EXTERNAL_FORCING )
            !-----------------------------------------------------------
            !  %%%%% CONNECTING TO GEOS-5 GCM via ESMF INTERFACE %%%%%
            !
            ! Call CHEMDR with dimension values as well as the Input 
            ! Options, Meteorology State, and Chemistry State objects. 
            ! (bmy, 11/8/12)
            !-----------------------------------------------------------

            ! Call gas-phase chemistry
            CALL CHEMDR( am_I_Root, NI,        NJ,        NL, 
     &                   Input_Opt, State_Met, State_Chm, RC  )

            ! Call linearized stratospheric scheme
            IF ( LSCHEM ) THEN 
               CALL DO_STRAT_CHEM( am_I_Root, Input_Opt, State_Chm, RC )
            ENDIF
#else
            !-----------------------------------------------------------
            !            %%%%% TRADITIONAL GEOS-Chem %%%%%
            !
            ! Current practice in the std GEOS-Chem is to call both
            ! CHEMDR and DO_STRAT_CHEM with just the am_I_Root
            ! argument. (bmy, 11/8/12)
            !-----------------------------------------------------------

            ! Call gas-phase chemistry
            CALL CHEMDR( am_I_Root )

            ! Call linearized stratospheric scheme
            IF ( LSCHEM ) CALL DO_STRAT_CHEM( am_I_Root )
#endif

            ! Do seasalt aerosol chemistry
            IF ( LSSALT ) CALL CHEMSEASALT( am_I_Root )

            ! Also do sulfate chemistry
            IF ( LSULF ) THEN

               ! Do sulfate chemistry
               CALL CHEMSULFATE( am_I_Root )
 
               ! Do aerosol thermodynamic equilibrium
               IF ( LSSALT ) THEN
                
#if   !defined( NO_ISORROPIA )
                  ! ISOROPIA takes Na+, Cl- into account
                  CALL DO_ISOROPIAII( am_I_Root )
#endif
                
               ELSE
 
                  ! RPMARES does not take Na+, Cl- into account
                  CALL DO_RPMARES( am_I_Root )
 
               ENDIF
                
            ENDIF
 
            ! Do carbonaceous aerosol chemistry
            IF ( LCARB ) CALL CHEMCARBON( am_I_Root )
 
            ! Do dust aerosol chemistry
            IF ( LDUST ) CALL CHEMDUST( am_I_Root )
 
            ! ND44 drydep fluxes
            ! The drydep fluxes diag. is done in vdiff_mod.f when non-local 
            ! PBL is used ( Lin, 03/31/09) 
            IF ( .NOT. LNLPBL ) CALL DRYFLX 

            ! ND43 chemical production
            CALL DIAGOH
 
            ! Remove acetone ocean sink
            IF ( IDTACET /= 0 ) THEN
               CALL OCEAN_SINK_ACET( STT(:,:,1,IDTACET) ) 
            ENDIF
             
         !---------------------------------
         ! Offline aerosol simulation
         !---------------------------------
         ELSE IF ( IT_IS_AN_AEROSOL_SIM ) THEN

            ! Define loop index and other SMVGEAR arrays
            ! N_TROP, the # of trop boxes, is returned
            CALL AEROSOL_RURALBOX( N_TROP, am_I_Root )
 
            ! Initialize FAST-J quantities for computing AOD's
            IF ( FIRST ) THEN

               ! Init chem mechanism and photolysis
               CALL READER  ( FIRST,         am_I_Root )
               CALL READCHEM(                am_I_Root )
               CALL INPHOT  ( LLTROP, NPHOT, am_I_Root )
 
               ! Reset NCS with NCSURBAN
               NCS     = NCSURBAN
 
               ! Reset NTLOOP and NTTLOOP after call to READER
               ! with the actual # of boxes w/in the ann mean trop
               NTLOOP  = N_TROP
               NTTLOOP = N_TROP
 
               ! Reset first-time flag
               FIRST = .FALSE.
            ENDIF
 
            ! Compute aerosol & dust concentrations [kg/m3]
            ! (NOTE: SOILDUST in "aerosol_mod.f" is computed here)
            CALL AEROSOL_CONC
 
            ! Compute AOD's and surface areas at 999 nm
            MONTH      = 0
            YEAR       = 0
            WAVELENGTH = 0
            CALL RDAER( am_I_Root, MONTH, YEAR, WAVELENGTH )
 
            !*** AEROSOL THERMODYNAMIC EQUILIBRIUM ***
            IF ( LSSALT ) THEN
             
#if   !defined( NO_ISORROPIA )
               ! ISOROPIA takes Na+, Cl- into account
               CALL DO_ISOROPIAII( am_I_Root )
#endif
             
            ELSE
 
               ! RPMARES does not take Na+, Cl- into account
               ! (skip for crystalline & aqueous offline run)
               IF ( .not. LCRYST ) CALL DO_RPMARES( am_I_Root )
 
            ENDIF
 
            !*** SEASALT AEROSOLS ***
            IF ( LSSALT ) CALL CHEMSEASALT( am_I_Root )
 
            !*** SULFATE AEROSOLS ***
            IF ( LSULF .or. LCRYST ) THEN
 
               ! Do sulfate chemistry
               CALL CHEMSULFATE( am_I_Root )
 
            ENDIF
                
            !*** CARBON AND 2NDARY ORGANIC AEROSOLS ***
            IF ( LCARB ) CALL CHEMCARBON( am_I_Root )
 
            !*** MINERAL DUST AEROSOLS ***
            IF ( LDUST ) THEN 
 
               ! Do dust aerosol chemsitry
               CALL CHEMDUST( am_I_Root )
 
               WAVELENGTH = 0

               ! Compute dust OD's & surface areas
               CALL RDUST_ONLINE( SOILDUST, WAVELENGTH )
            ENDIF
 
         !---------------------------------
         ! Rn-Pb-Be
         !---------------------------------                 
         ELSE IF ( IT_IS_A_RnPbBe_SIM ) THEN
            CALL CHEMRnPbBe 
            CALL DRYFLXRnPbBe
                   
         !---------------------------------
         ! CH3I
         !---------------------------------
         ELSE IF ( IT_IS_A_CH3I_SIM ) THEN
            CALL CHEMCH3I
 
         !---------------------------------            
         ! HCN
         !---------------------------------
         ELSE IF ( IT_IS_A_HCN_SIM ) THEN
            CALL CHEM_HCN_CH3CN( N_TRACERS, STT )
 
         !---------------------------------
         ! Tagged O3
         !---------------------------------
         ELSE IF ( IT_IS_A_TAGOX_SIM ) THEN 
            CALL CHEM_TAGGED_OX
 
#if defined( DEVEL )
            !-----------------------------------------------------------
            !  %%%%% CONNECTING TO GEOS-5 GCM via ESMF INTERFACE %%%%%
            ! 
            ! Call DO_STRAT_CHEM (the G-C chemistry driver routine) 
            ! with the Input Options and Chemistry State objects. 
            ! (bmy, 11/9/12)
            !-----------------------------------------------------------

            ! Call linearized stratospheric scheme
            IF ( LSCHEM ) THEN 
               CALL DO_STRAT_CHEM( am_I_Root, Input_Opt, State_Chm, RC )
            ENDIF
#else
            !-----------------------------------------------------------
            !             %%%%% TRADITIONAL GEOS-Chem %%%%%
            !
            ! Current practice is to call DO_STRAT_CHEM with am_I_Root 
            ! as the only argument. (bmy, 11/9/12)
            !-----------------------------------------------------------

            ! Call linearized stratospheric scheme
            IF ( LSCHEM ) CALL DO_STRAT_CHEM( am_I_Root )
#endif

         !---------------------------------
         ! Tagged CO
         !---------------------------------
         ELSE IF ( IT_IS_A_TAGCO_SIM ) THEN
            CALL CHEM_TAGGED_CO
 
         !---------------------------------
         ! C2H6
         !---------------------------------
         ELSE IF ( IT_IS_A_C2H6_SIM ) THEN
            CALL CHEMC2H6
 
         !---------------------------------
         ! CH4
         !---------------------------------
         ELSE IF ( IT_IS_A_CH4_SIM ) THEN
 
            ! Only call after the first 24 hours
            IF ( GET_ELAPSED_MIN() >= GET_TS_CHEM() ) THEN
               CALL CHEMCH4
            ENDIF
 
         !---------------------------------
         ! Mercury
         !---------------------------------
         ELSE IF ( IT_IS_A_MERCURY_SIM ) THEN
 
            ! Do Hg chemistry
            CALL CHEMMERCURY
               
         !---------------------------------
         ! Offline H2/HD
         !---------------------------------
         ELSE IF ( IT_IS_A_H2HD_SIM ) THEN
            CALL CHEM_H2_HD
            CALL DRYFLXH2HD

#if defined( DEVEL )
            !-----------------------------------------------------------
            !  %%%%% CONNECTING TO GEOS-5 GCM via ESMF INTERFACE %%%%%
            ! 
            ! Call CHEMDR (the G-C chemistry driver routine) with
            ! local dimension values and the the Meteorology and 
            ! Chemistry State objects.  Have CHEMDR return success or
            ! failure via the RC argument. (bmy, 10/25/12)
            !-----------------------------------------------------------

            ! Call linearized stratospheric scheme
            IF ( LSCHEM ) THEN 
               CALL DO_STRAT_CHEM( am_I_Root, Input_Opt, State_Chm, RC )
            ENDIF
#else
            !-----------------------------------------------------------
            !             %%%%% TRADITIONAL GEOS-Chem %%%%%
            !
            ! Current practice is to call DO_STRAT_CHEM with am_I_Root 
            ! as the only argument. (bmy, 10/25/12)
            !-----------------------------------------------------------

            ! Call linearized stratospheric scheme
            CALL DO_STRAT_CHEM( am_I_Root )
#endif  

         ENDIF

         !### Debug
         IF ( LPRT .and. am_I_Root ) THEN
            CALL DEBUG_MSG( '### MAIN: a CHEMISTRY' )
         ENDIF
      ENDIF
      
#if defined( DEVEL ) || defined( EXTERNAL_GRID ) || defined( EXTERNAL_FORCING )
      !-----------------------------------------------------------------
      !     %%%%% CONNECTING TO GEOS-5 GCM via ESMF INTERFACE %%%%%
      !
      ! At present (25 Oct 2012), we have not yet completely replaced
      ! the GEOS-Chem arrays with fields from the Chemistry State 
      ! (i.e. State_Chm) and Meteorology State (i.e. State_Met)
      ! objects.  In order to continue GIGC development without 
      ! disrupting the existing GEOS-Chem workflow, we shall copy the 
      ! tracer array STT back to State_Chm%TRACERS after chemistry. 
      !
      ! Also note: State_Chm%TRACERS needs to leave this routine with
      ! units of [kg].  If we are calling this routine from the 
      ! traditional "main.F" driver, then the tracers are already in
      ! [kg] and no further unit conversion is necessary.  If on the
      ! other hand we are calling this routine from the ESMF interface,
      ! then we need to convert State_Chm%TRACERS from [kg] to [v/v]
      ! in the calling routine GIGC_Do_Chem (in ESMF/gigc_chemdr.F90).
      ! (bmy, 10/24/12) 
      !-----------------------------------------------------------------

      ! Reset Chemistry State tracers [kg] from GEOS-Chem tracer array
      State_Chm%TRACERS(:,:,:,1:N_TRACERS) = STT(:,:,:,1:N_TRACERS)
#endif

      END SUBROUTINE DO_CHEMISTRY
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: gckpp_driver
!
! !DESCRIPTION: Subroutine GCKPP\_DRIVER is the driver routine to perform 
!  integration with the full KPP chemistry mechanism.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE GCKPP_DRIVER( KTLOOP, JLOOPLO, R_KPP, NSPEC_GC ) 
!
! !USES:
!
      USE COMODE_MOD,           ONLY : JLOP,   CSPEC
      USE COMODE_MOD,           ONLY : IXSAVE, IYSAVE,    IZSAVE
      USE GCKPP_COMODE_MOD,     ONLY : HSAVE_KPP 
      USE TIME_MOD,             ONLY : GET_TS_CHEM
      USE GCKPP_UTIL,           ONLY : SHUFFLE_KPP2USER
      USE GCKPP_INITIALIZE,     ONLY : INITIALIZE
      USE GCKPP_MODEL
      USE GCKPP_GLOBAL    
      USE GCKPP_RATES,          ONLY : UPDATE_RCONST
      USE GCKPP_MONITOR,        ONLY : SPC_NAMES
      USE GCKPP_FUNCTION
      USE ERROR_MOD,            ONLY : ERROR_STOP
      USE GCKPP_INTEGRATOR,     ONLY : NHNEW, NHEXIT, INTEGRATE
!
! !INPUT PARAMETERS: 
!
      INTEGER, INTENT(IN) :: KTLOOP       ! Local loop index
      INTEGER, INTENT(IN) :: JLOOPLO      ! JLOOPLO + KLOOP = JLOOP
      REAL*8,  INTENT(IN) :: R_KPP(:,:)   ! Array of reaction rates
      INTEGER, INTENT(IN) :: NSPEC_GC     ! # of active chemical species
!
! !REMARKS:
!  Variables used to pass the last/first step size b/w call 
!                                                                             .
!  For Rosenbrock:
!  ----------------
!      Nhexit=2, Nhnew = 3
!  OUT    
!      RSTATUS(2)  -> Hexit, last accepted step before exit
!      RSTATUS(3)  -> Hnew, last predicted step (not yet taken)
!      For multiple restarts, use Hnew as Hstart in the subsequent run
!  IN
!      RCNTRL(3)   -> Hstart, starting value for the integration step size
!                                                                             .
!                                                                             .
!  For LSODE:
!  ------------
!  OUT
!      RSTATUS(1)  -> Texit, the time corresponding to the
!                     computed Y upon return
!      RSTATUS(2)  -> Hexit, last predicted step before exit
!      For multiple restarts, use Hexit as Hstart in the following run
!
!  IN      
!      RCNTRL(3)  -> Hstart, starting value for the integration step size
!                                                                             .
!                                                                             .
!  For RADAU5:
!  ------------
!  OUT
!      RSTATUS(1)  -> final time     
!  IN      
!      RCNTRL(3)   -> not used
!                                                                             .
!                                                                             .
! For RUNGE-KUTTA
! ---------------
! OUT
!     same as Rosenbrock
! 
! !REVISION HISTORY: 
!  24 Jan 2008 - Kumaresh    - Based on Daven Henze's GCKPP_DRIVER.
!  16 Sep 2009 - R. Yantosca - Commented, and updated to call various 
!  03 Dec 2009 - C. Carouge  - Use CSPEC instead of CSPEC_FOR_KPP 
!                              to save memory space
!  17 Dec 2009 - R. Yantosca - Added ProTeX headers
!  20 Jan 2010 - C. Carouge  - Now call GCKPP_DRIVER from physproc.f to save 
!                              memory. 
!  20 Jan 2010 - C. Carouge  - Now use the # of active species from GC to 
!                              update CSPEC and not the of variable species 
!                              from KPP.
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      REAL*8        :: T, TIN, TOUT
      INTEGER       :: ICNTRL(20)
      REAL(kind=dp) :: RCNTRL(20)
      INTEGER       :: ISTATUS(20)
      INTEGER       :: I, J, L, N, KLOOP, IERR
      REAL(kind=dp) :: RSTATE(20)
      LOGICAL, SAVE :: FIRST = .TRUE.
      CHARACTER(LEN=255) :: ERR_MSG
      
      !=================================================================

      STEPMIN = 0.0d0
      STEPMAX = 0.0d0

      ! Suggested tolerance for Rosenbrock solvers
      DO i=1,NVAR
!         RTOL(i) = 2.d-1
         RTOL(i) = 1.d-1
         ATOL(i) = 1.d-2
      END DO

      ! for LSODE
!      DO i=1,NVAR
!         RTOL(i) = 1.0d-3
!         ATOL(i) = 1.0d6
!      END DO

      
      ! Set parameters to default. See comments in Integrator module for
      ! a list of the defaults.
      ICNTRL(:) = 0
      RCNTRL(:) = 0.d0

      ! Change some parameters from the default to new values
      ICNTRL(1) = 1             ! Autonomous
      ICNTRL(2) = 0             ! Nonautonomous

      ! Select Integrator
      !    ICNTRL(3)  -> selection of a particular method.
      !                  For Rosenbrock, options are:
      !        = 0 :  default method is Rodas3
      !        = 1 :  method is  Ros2
      !        = 2 :  method is  Ros3 
      !        = 3 :  method is  Ros4 
      !        = 4 :  method is  Rodas3
      !        = 5:   method is  Rodas4
      ICNTRL(3) = 4    
      ICNTRL(7) = 1 ! No adjoint

      IF(FIRST)THEN

         RSTATE(Nhexit) = 0d0

         FIRST = .FALSE. 

      ENDIF

      ! GET TS_CHEM and convert it to seconds. 
      DT = GET_TS_CHEM() * 60d0

      ! Set time parameters. 
      T = 0d0
      TIN = T
      TOUT = T + DT

      !=================================================================
      ! Solve Chemistry
      !=================================================================
 100  format('No. of function calls:', i6, /,
     &       'No. of jacobian calls:', i6, /,
     &       'No. of steps:         ', i6, /,
     &       'No. of accepted steps:', i6, /,
     &       'No. of rejected steps ', i6, /,
     $       '       (except at very beginning)',          /,
     &       'No. of LU decompositions:             ', i6, /,
     &       'No. of forward/backward substitutions:', i6, /,
     &       'No. of singular matrix decompositions:', i6, /,
     &       /,
     &       'Texit, the time corresponding to the      ',        /,
     &       '       computed Y upon return:            ', f11.4, /,
     &       'Hexit, last accepted step before exit:    ', f11.4, /,
     &       'Hnew, last predicted step (not yet taken):', f11.4 )


      DO KLOOP = 1, KTLOOP

         ! 1D grid box index on the global grid
         JLOOP         = JLOOPLO+KLOOP
           
         ! Get 3D coords from SMVGEAR's 1D coords
         I = IXSAVE(JLOOP)
         J = IYSAVE(JLOOP)
         L = IZSAVE(JLOOP)

         ! Pass tracer concentrations from CSPEC to KPP working vector V_CSPEC
         DO N =1, NVAR
            V_CSPEC(N) = CSPEC(JLOOP,N)
         END DO

         ! Pass tracer concentrations from V_CSPEC to KPP working vector VAR.
         ! This also initializes the constant rate constants FIX.
         CALL Initialize()


         ! starting value for integration time step
         RCNTRL(3) = HSAVE_KPP(I,J,L)

         ! 1D grid box index in the KTLOOP subset.
         ! R_KPP is only defined on KTLOOP boxes (ccc, 12/3/09)
         JLOOP = KLOOP
         CALL Update_RCONST(R_KPP)
         
         ! 1D grid box index on the global grid
         JLOOP = JLOOPLO+KLOOP


         ! Integrate FWD -- phs --
         call integrate(TIN, TOUT, ICNTRL, RCNTRL, ISTATUS,
     $        RSTATE, IERR)
         
         ! Try another time if it failed
         IF ( IERR < 0 ) THEN
            write(6,*) ''
            write(6,100) istatus(1:8),rstate(1:3)
            write(6,*) ''
            write(ERR_MSG,'(a, i3)') 'Integrator error code :',IERR
            write(6,*)  'JLOOP, I, J, L ', JLOOP, I, J, L

            ! Reset first time step and start concentrations
            RCNTRL(3)  = 0d0
            CALL Initialize( )  ! v2.1 

            ! 1D grid box index in the KTLOOP subset.
            ! R_KPP is only defined on KTLOOP boxes (ccc, 12/3/09)
            JLOOP = KLOOP
            CALL Update_RCONST(R_KPP)

            ! 1D grid box index on the global grid
            JLOOP = JLOOPLO+KLOOP

            call integrate(TIN, TOUT, ICNTRL, RCNTRL, ISTATUS,
     $           RSTATE, IERR)
            IF ( IERR < 0 ) THEN              
               print*, 'failed twice !!! '
               CALL ERROR_STOP(ERR_MSG, 'INTEGRATE_KPP')
            ENDIF
                            
         ENDIF

!--------- using the adjoint integrator (limited to Rosenbrock)           
!adj!           ! using adjoint integrator           
!adj!           CALL INTEGRATE_ADJ( 1, VAR, VAR_ADJ, TIN, TOUT,ATOL, 
!adj!     &        RTOL, ICNTRL, RCNTRL, ISTATUS, RSTATE)
!adj!
!adj!
!adj!           ! Update the printout and functionality (dkh, 07/28/09)
!adj!           IF ( ISTATUS(NIERR) < 0 ) THEN 
!adj!              print*, 'KPP Integrator failed.  Trying again' 
!adj!              print*, 'IERR = ', RSTATE(20)
!adj!              print*, 'RSTAT = ', RSTATE
!adj!              print*, 'ISTAT = ', ISTATUS
!adj!              print*, 'RCNTRL = ', RCNTRL
!adj!              print*, 'ICNTRL = ', ICNTRL
!adj!              print*, 'JLOOP, I, J, L ', JLOOP, I, J, L
!adj!              rcntrl(3)  = 0d0
!adj!              CALL Initialize( ) ! v2.1 
!adj!              CALL Update_RCONST()
!adj!              CALL INTEGRATE_ADJ( 1, VAR, VAR_ADJ, TIN, TOUT,
!adj!     &             ATOL, RTOL, ICNTRL, RCNTRL, ISTATUS, RSTATE)
!adj!              IF ( ISTATUS(NIERR) < 0 ) THEN
!adj!                 print*, 'failed twice !!! '
!adj!                 CALL ERROR_STOP('IERR < 0 ', 'INTEGRATE_ADJ')
!adj!              ENDIF
!adj!           ENDIF
!--------- using the adjoint integrator (limited to Rosenbrock)           


         ! Set negative values to SMAL2
         DO N = 1, NVAR
            VAR(N) = MAX(VAR(N),SMAL2)
         ENDDO

         ! save next integration time step
         HSAVE_KPP(I,J,L) = RSTATE(Nhnew)

         ! Pass KPP concentrations from VAR to V_CSPEC
         CALL Shuffle_kpp2user(VAR,V_CSPEC)  

         ! Pass KPP concentrations V_CSPEC to geos-chem CSPEC
         DO N =1, NSPEC_GC
            CSPEC(JLOOP,N) = V_CSPEC(N)
         ENDDO

      ENDDO

      ! Return to calling program
      END SUBROUTINE GCKPP_DRIVER
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: recompute_od
!
! !DESCRIPTION: Subroutine RECOMPUTE\_OD will update the optical depth values 
!  before accumulating or writing the diagnostics.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE RECOMPUTE_OD( am_I_Root )
!
! !USES:
!
      ! References to F90 modules
      USE AEROSOL_MOD,     ONLY : AEROSOL_CONC
      USE AEROSOL_MOD,     ONLY : RDAER, SOILDUST
      USE DUST_MOD,        ONLY : RDUST_ONLINE, RDUST_OFFLINE
      USE ERROR_MOD,       ONLY : DEBUG_MSG
      USE LOGICAL_MOD,     ONLY : LCARB, LCHEM, LCRYST, LDUST
      USE LOGICAL_MOD,     ONLY : LPRT, LSSALT, LSULF,  LSOA
      USE TIME_MOD,        ONLY : GET_MONTH, GET_YEAR
      USE TRACER_MOD,      ONLY : ITS_A_FULLCHEM_SIM
      USE TRACER_MOD,      ONLY : ITS_AN_AEROSOL_SIM
!
! !INPUT PARAMETERS:
!
      LOGICAL, INTENT(IN) :: am_I_Root   ! Is this the root CPU?
!
! !REVISION HISTORY: 
!  03 Fev 2011 - Adapted from chemdr.f by skim
!  30 Jul 2012 - R. Yantosca - Now accept am_I_Root as an argument when
!                              running with the traditional driver main.F
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!

      INTEGER :: MONTH, YEAR, WAVELENGTH
 
      !=================================================================
      ! RECOMPUTE_OD begins here!
      !=================================================================

      ! Get month and year
      MONTH = GET_MONTH()
      YEAR  = GET_YEAR()

      IF ( LCHEM ) THEN

         IF ( ITS_A_FULLCHEM_SIM() .or. ITS_AN_AEROSOL_SIM() ) THEN
            IF ( LSULF .or. LCARB .or. LDUST .or. LSSALT ) THEN

               ! Skip this section if all of these are turned off
               CALL AEROSOL_CONC

               !===============================================================
               ! Call RDAER -- computes aerosol optical depths
               !===============================================================

               ! Calculate the AOD at the wavelength specified in jv_spec_aod
               WAVELENGTH = 1
               CALL RDAER( am_I_Root, MONTH, YEAR, WAVELENGTH )

               !### Debug
               IF ( LPRT .and. am_I_Root ) THEN 
                  CALL DEBUG_MSG( '### RECOMPUTE_OD: after RDAER' )
               ENDIF

               !===============================================================
               ! If LDUST is turned on, then we have online dust aerosol in
               ! GEOS-CHEM...so just pass SOILDUST to RDUST_ONLINE in order to
               ! compute aerosol optical depth for FAST-J, etc.
               !
               ! If LDUST is turned off, then we don't have online dust aerosol
               ! in GEOS-CHEM...so read monthly-mean dust files from disk.
               ! (rjp, tdf, bmy, 4/1/04)
               !===============================================================
               IF ( LDUST ) THEN
                  CALL RDUST_ONLINE( SOILDUST, WAVELENGTH )
               ELSE
                  CALL RDUST_OFFLINE( MONTH,      YEAR, 
     &                                WAVELENGTH, am_I_Root )
               ENDIF

               !### Debug
               IF ( LPRT .and. am_I_Root ) THEN
                  CALL DEBUG_MSG( '### RECOMPUTE_OD: after RDUST' )
               ENDIF

            ENDIF
         ENDIF

      ENDIF

      END SUBROUTINE RECOMPUTE_OD
!EOC
      END MODULE CHEMISTRY_MOD
