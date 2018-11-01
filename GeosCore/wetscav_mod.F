!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!     
! !MODULE: wetscav_mod.F
!     
! !DESCRIPTION: Module WETSCAV\_MOD contains routines and variables used in 
!  the wet scavenging of species in cloud updrafts, rainout, and washout. 
!\\   
!\\   
! !INTERFACE: 
!
      MODULE WETSCAV_MOD
!
! !USES:
!
      USE PRECISION_MOD    ! For GEOS-Chem Precision (fp)

      IMPLICIT NONE
      PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
      PUBLIC  :: CLEANUP_WETSCAV
      PUBLIC  :: COMPUTE_F
      PUBLIC  :: DO_WETDEP
      PUBLIC  :: INIT_WETSCAV
      PUBLIC  :: SETUP_WETSCAV
      PUBLIC  :: WASHOUT
      PUBLIC  :: LS_K_RAIN
      PUBLIC  :: LS_F_PRIME
!
! !PRIVATE MEMBER FUNCTIONS:
!
      PRIVATE :: COMPUTE_L2G
      PRIVATE :: CONV_F_PRIME
      PRIVATE :: E_ICE
      PRIVATE :: RAINOUT
      PRIVATE :: GET_RAINFRAC
      PRIVATE :: SAFETY
      PRIVATE :: WASHFRAC_FINE_AEROSOL
      PRIVATE :: WASHFRAC_COARSE_AEROSOL
      PRIVATE :: WASHFRAC_LIQ_GAS
      PRIVATE :: WASHFRAC_HNO3
      PRIVATE :: GET_VUD
!
! !REMARKS:
!  References:
!  ============================================================================
!  (1 ) Liu,H., D.J. Jacob, I. Bey and R.M. Yantosca, "Constraints from 210Pb 
!        and 7Be on wet deposition and transport in a global three-dimensional
!        chemical tracer model driven by assimilated meteorological fields", 
!        JGR, Vol 106, pp 12109-12128, 2001.
!  (2 ) D.J. Jacob, H. Liu, C. Mari, and R. M. Yantosca, "Harvard wet 
!        deposition scheme for GMI", Harvard Atmospheric Chemistry Modeling 
!        Group, March 2000.
!  (3 ) Chin, M., D.J. Jacob, G.M. Gardner, M.S. Foreman-Fowler, and P.A. 
!        Spiro, "A global three-dimensional model of tropospheric sulfate", 
!        J. Geophys. Res., 101, 18667-18690, 1996.
!  (4 ) Balkanski, Y  D.J. Jacob, G.M. Gardner, W.C. Graustein, and K.K.
!        Turekian, "Transport and Residence Times of Tropospheric Aerosols
!        from a Global Three-Dimensional Simulation of 210Pb", JGR, Vol 98, 
!        (D11) pp 20573-20586, 1993.  
!  (5 ) Giorgi, F, & W.L. Chaimedes, "Rainout Lifetimes of Highly Soluble
!        Aerosols and Gases as Inferred from Simulations With a General
!        Circulation Model", JGR, Vol 86 (D13) pp 14367-14376, 1986.  
!
! !REVISION HISTORY:
!  (1 ) Now trap allocation errors with routine ALLOC_ERR. (bmy, 7/11/00)
!  (2 ) Moved routine MAKE_QQ here from "dao_mod.f" (bmy, 10/12/00)
!  (3 ) Reordered arguments in INIT_PRECIP (bmy, 10/12/00)
!  (4 ) Updated comments (bmy, 9/4/01)
!  (5 ) Bug fix in MAKE_QQ: BXHEIGHT is sized IIPAR,JJPAR,LLPAR (bmy, 10/4/01)
!  (6 ) Removed obsolete, commented-out code from 10/01 (bmy, 11/26/01)
!  (7 ) Now divide module header into MODULE PRIVATE, MODULE VARIABLES, and
!        MODULE ROUTINES sections.  Updated comments (bmy, 5/28/02)
!  (8 ) Now zero allocatable arrays (bmy, 8/5/02)
!  (9 ) Bug fix: ND39 diagnostic now closes the budget.  Also bundled several
!        standalone routines into this module.  Now references F90 module
!        "tracerid_mod.f".  Also set NSOLMAX=10 since we now have sulfate
!        tracers for wetdep.   Now prevent out-of-bounds errors in routine
!        WETDEP.  Added GET_WETDEP_NMAX function to return max # of soluble
!        tracers for allocating diagnostic arrays.  Added functions 
!        GET_WETDEP_NSOL and GET_WETDEP_IDWETD.  Now init H2O2s and SO2s
!        to the initial H2O2 and SO2 from STT.  Updated comments. 
!        (qli, bmy, 1/14/03)
!  (10) Improvements for SO2/SO4 scavenging (rjp, bmy, 3/23/03)
!  (11) Now references "time_mod.f".  Added driver routine DO_WETDEP to
!        remove cumbersome calling sequence from MAIN program.  Also declared
!        WETDEP and MAKE_QQ PRIVATE to this module. (bmy, 3/27/03)
!  (11) Add parallelization to routine WETDEP (bmy, 3/17/04)
!  (12) Added carbon and dust aerosol tracers (rjp, tdf, bmy, 4/5/04)
!  (13) Added seasalt aerosol tracers (rjp, bec, bmy, 4/20/04)
!  (14) Added secondary organic aerosol tracers (rjp, bmy, 7/13/04)
!  (15) Now references "logical_mod.f" and "tracer_mod.f".  Now move all 
!        internal routines to the module and pass arguments explicitly in
!        order to facilitate parallelization on the Altix. (bmy, 7/20/04)
!  (16) Updated for mercury aerosol tracers (eck, bmy, 12/9/04)
!  (17) Updated for AS, AHS, LET, NH4aq, SO4aq.  Also now pass Hg2 wetdep loss
!        to "ocean_mercury_mod.f". (cas, sas, bmy, 1/20/05)
!  (18) Bug fix to avoid numerical blowup in WETDEP.  Now use analytical
!        function for E_ICE(T). (bmy, 3/7/05)
!  (19) Added SO4s, NITs.  Increased NSOLMAX to 31.  Also block out 
!        parallel loop in WETDEP for SGI MIPS compiler. (bec, bmy, 5/5/05)
!  (20) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (21) Bug fixes: do not over-deplete H2O2s.  Also include updates for
!        tagged Hg simulation. (dkh, rjp, eck, cdh, bmy, 1/6/06)
!  (22) Now wet deposit SOG4, SOA4. Remove unnecessary variables in WETDEP.
!        (dkh, bmy, 5/18/06)
!  (23) Bug fixes in COMPUTE_F (bmy, 7/26/06)
!  (24) Resize DSTT array in WETDEP to save memory.  Added fixes for GEOS-5
!        wet deposition per Hongyu Liu's suggestions. (bmy, 3/5/08)
!  (25) Add wet scavenging of GLYX, MGLY, GLYC, SOAG, SOAM (tmf, 1/7/09)
!  (26) Effective Henry's law constant and coefficient from 
!       Sander, R, 1999, Compilation of Henry's Law Constants for 
!          Inorganic and Organic Species of Potential Importance in 
!          Environmental Chemistry.
!          http://www.mpch-mainz.mpg.de/~sander/res/henry.html
!       (tmf, 1/7/09)
!  (27) Remove support for SGI compiler.  Bug fix in RAINOUT. (bmy, 7/20/09)
!  (28) Update mercury simulation. (ccc, 5/17/10)
!  (29) Add LGTMM as condition to output AD39. (ccc, 11/18/09)
!  (30) Add snow scavenging, different washout/rainout ratio 
!       (wqq, ccc, 7/13/10)
!  13 Aug 2010 - R. Yantosca - Add modifications for MERRA (treat like GEOS-5)
!  16 Sep 2010 - R. Yantosca - Added ProteX headers
!  20 Sep 2010 - H. Amos, R. Yantosca - Implement new algorithms for MERRA
!  08 Oct 2010 - H. Amos     - WASHFRAC_LIQ_GAS is now a subroutine
!  08 Oct 2010 - H. Amos     - Various other modifications in WETDEP_MERRA
!  01 Aug 2011 - H. Amos     - Bug fix for function WASHFRAC_LIQ_GAS
!  01 Aug 2011 - H. Amos     - Updated comments
!  09 Feb 2012 - R. Yantosca - Add modifications for GEOS-5.7.x
!  21 Jun 2012 - R. Yantosca - Declare H2O2s, SO2s as TARGETs for pointers
!  23 Apr 2013 - R. Yantosca - Bug fix, eliminate white space from #if block
!  23 Apr 2013 - R. Yantosca - Remove LTOMAS logical, since we now invoke TOMAS
!                              with either TOMAS=yes or TOMAS40=yes
!  13 Aug 2013 - M. Sulprizio- Add modifications for updated SOA and SOA + 
!                              semivolatile POA simulations (H. Pye)
!  12 Sep 2013 - M. Sulprizio- Add modifications for acid uptake on dust
!                              aerosols (T.D. Fairlie)
!  24 Nov 2014 - M. Yannetti - Added PRECISION_MOD
!  12 Feb 2015 - C. Keller   - Added GET_F and GET_VUD. VUD is now stored in
!                              State_Met. In an ESMF environment, it's values
!                              are directly taken from the GEOS-5 model. In 
!                              the traditional GEOS-Chem, the updraft velocity
!                              is still calculated in the old manner (5 m/s 
!                              over water, 10 m/s over land/ice).
!  23 Jun 2015 - M. Sulprizio- Add impaction scavenging for hydrophobic BC and
!                              homogeneous IN removal from Qiaoqiao Wang
!  08 Jul 2015 - E. Lundgren - Add marine organic aerosols (B.Gantt, M.Johnson)
!  03 Sep 2015 - R. Yantosca - Remove NSOLMAX, we now get the # of wetdep
!                              species from State_Chm%nWetDep
!  03 Sep 2015 - R. Yantosca - INIT_WETSCAV now just allocates arrays;
!                              SETUP_WETSCAV now just sets up quantities
!  03 Sep 2015 - R. Yantosca - Remove WETDEPID routine, it's not needed
!  22 Sep 2015 - R. Yantosca - Re-introduce GET_WETDEP_NSOL routine
!  25 Sep 2015 - R. Yantosca - Use species database to collapse down routines
!                              COMPUTE_F, RAINOUT, WASHOUT
!  16 Jun 2016 - L. Hu       - Replaced IDTXXX with Ind_('XXX')
!  20 Jun 2016 - R. Yantosca - Now define species ID's only in the INIT phase
!  20 Jun 2016 - R. Yantosca - Renamed species ID's IDTxxxx to id_XXXX
!  06 Jul 2016 - R. Yantosca - Removed GET_WETDEP_NSOL and GET_WETDEP_IDWETD,
!                              these are now replaced by State_Chm fields
!  17 Mar 2017 - R. Yantosca - Finally remove NSOL and IDWETD; they're unused
!  24 Aug 2017 - M. Sulprizio- Remove support for GCAP, GEOS-4, GEOS-5 and MERRA
!                              and rename routine WETDEP_MERRA to WETDEP (the
!                              old WETDEP routine for GEOS-5 met was removed)
!  13 Nov 2017 - R. Yantosca - Added arrays for State_Diag diagnostics
!  11 Oct 2018 - M. Sulprizio- Move H2O2s and SO2s to State_Chm and rename as
!                              H2O2AfterChem and SO2AfterChem
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
      ! TINY number
      REAL(fp), PARAMETER           :: TINY_FP = TINY(1.0_fp)

      ! Use condensed water content (L+W) =1.0e-6 cm3/cm3 (= 1.0e-3 kg/m3).
      ! which was recommended by Qiaoqiao Wang et al (2014).
      REAL(fp), PARAMETER           :: COND_WATER_CONTENT = 1.0e-6_fp
!
! !LOCAL VARIABLES:
!
      ! Arrays
      REAL(fp), ALLOCATABLE, TARGET :: C_H2O (:,:,:) ! Mix ratio of H2O [v/v]
      REAL(fp), ALLOCATABLE, TARGET :: CLDICE(:,:,:) ! Cloud ice mixing ratio
                                                     !  [cm3 ice/cm3 air]
      REAL(fp), ALLOCATABLE, TARGET :: CLDLIQ(:,:,:) ! Cloud liquid water
                                                     !  mixing ratio
                                                     !  [cm3 H2O/cm3 air]
      REAL(fp), ALLOCATABLE         :: PDOWN (:,:,:) ! Precipitation thru the
                                                     !  bottom of the grid box
                                                     !  [cm3 H2O/cm2 area/s]
      REAL(fp), ALLOCATABLE         :: QQ    (:,:,:) ! Rate of new precip
                                                     !  formation
                                                     !  [cm3 H2O/cm3 air/s]
      REAL*8,  ALLOCATABLE          :: REEVAP(:,:,:) ! Rate of precip 
                                                     !  reevaporation
                                                     !  [cm3 H2O/cm3 air/s]

      ! Define local shadow variables for values in Input_Opt
      LOGICAL                       :: LPRT
      LOGICAL                       :: LGTMM
      LOGICAL                       :: LSOILNOX
      LOGICAL                       :: LDYNOCEAN
      LOGICAL                       :: ITS_A_MERCURY_SIM
      LOGICAL                       :: ITS_A_POPS_SIM

      ! Species ID flags
      INTEGER                       :: id_DUST1
      INTEGER                       :: id_H2O2
      INTEGER                       :: id_NK1
      INTEGER                       :: id_SF1
      INTEGER                       :: id_SO2
      INTEGER                       :: id_SO4

      CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: do_wetdep
!
! !DESCRIPTION: Subroutine DO\_WETDEP is a driver for the wet deposition code, 
!  called from the MAIN program.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE DO_WETDEP( am_I_Root, Input_Opt,  State_Met,
     &                      State_Chm, State_Diag, RC         )
!
! !USES:
!
      USE CMN_SIZE_MOD
#if defined( NC_DIAG )
      USE Diagnostics_Mod, ONLY : Compute_Budget_Diagnostics
      USE Diagnostics_Mod, ONLY : Compute_Column_Mass
#endif
      USE ErrCode_Mod
      USE Error_Mod,       ONLY : Debug_MSG
      USE Input_Opt_Mod,   ONLY : OptInput
      USE PhysConstants
      USE State_Chm_Mod,   ONLY : ChmState
      USE State_Diag_Mod,  ONLY : DgnState
      USE State_Met_Mod,   ONLY : MetState
#if defined( USE_TEND )
      USE TENDENCIES_MOD
#endif
      USE TIME_MOD,       ONLY : GET_TS_DYN
      USE UnitConv_Mod
!
! !INPUT PARAMETERS:
!
      LOGICAL,        INTENT(IN)    :: am_I_Root   ! Are we on the root CPU?
      TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
      TYPE(MetState), INTENT(IN)    :: State_Met   ! Meteorology State object
!
! !INPUT/OUTPUT PARAMETERS:
!
      TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
      TYPE(DgnState), INTENT(INOUT) :: State_Diag  ! Diagnostics State object
!
! !OUTPUT PARAMETERS:
!
      INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
!
! !REVISION HISTORY: 
!  27 Mar 2003 - R. Yantosca - Initial version
!  (1 ) Now references LPRT from "logical_mod.f" (bmy, 7/20/04)
!  (2 ) Don't do rainout/washout for conv precip for GEOS-5 (hyl, bmy, 3/5/08)
!  13 Aug 2010 - R. Yantosca - Treat GEOS-5 like MERRA
!  16 Sep 2010 - R. Yantosca - Added ProTeX headers
!  20 Sep 2010 - R. Yantosca - Rewrote #if block structure for clarity
!  09 Feb 2012 - R. Yantosca - Treat GEOS-5.7.x in the same way as MERRA
!  25 Mar 2013 - R. Yantosca - Now accept am_I_Root, Input_Opt, State_Chm, RC
!  02 Jul 2013 - R. Yantosca - Bug fix: add State_Chm to WETDEP calls for
!                              met other than MERRA, GEOS-5.7 or GEOS-5
!  26 Sep 2013 - R. Yantosca - Renamed GEOS_57 Cpp switch to GEOS_FP
!  08 Aug 2015 - E. Lundgren - Tracer units are now input as [kg/kg]
!  12 Aug 2015 - R. Yantosca - Add support for MERRA2 meteorology
!  29 Apr 2016 - R. Yantosca - Don't initialize pointers in declaration stmts
!  19 Jul 2016 - R. Yantosca - Now bracket tendency calls with #ifdef USE_TEND
!  09 Nov 2017 - R. Yantosca - Now accept State_Diag as an argument
!  09 Nov 2017 - R. Yantosca - Return error condition to calling program
!  27 Aug 2018 - E. Lundgren - Implement budget diagnostics
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      ! Scalars
      LOGICAL                 :: doPrt
      INTEGER                 :: I, J, L
      REAL(fp)                :: DT_TEND

      ! Strings
      CHARACTER(LEN=255)      :: ErrMsg, ThisLoc

#if defined( NC_DIAG )
      REAL(fp)                :: DT_Dyn
#endif

      !=================================================================
      ! Initialize
      !=================================================================

      ! Initialize
      RC      = GC_SUCCESS
      ErrMsg  = ''
      ThisLoc = ' -> at Do_WetDep (in module GeosCore/wetscav_mod.F)'

      ! Copy values from Input_Opt to module shadow variables
      LPRT              = Input_Opt%LPRT
      LGTMM             = Input_Opt%LGTMM
      LSOILNOX          = Input_Opt%LSOILNOX
      LDYNOCEAN         = Input_Opt%LDYNOCEAN
      ITS_A_MERCURY_SIM = Input_Opt%ITS_A_MERCURY_SIM
      ITS_A_POPS_SIM    = Input_Opt%ITS_A_POPS_SIM

      ! Only print output on the root CPU
      doPrt = ( LPRT .and. am_I_Root )

#if defined( NC_DIAG )
      !----------------------------------------------------------
      ! Wet deposition budget diagnostics - Part 1 of 2
      !----------------------------------------------------------
      IF ( State_Diag%Archive_BudgetWetDep ) THEN
         ! Get initial column masses
         CALL Compute_Column_Mass( am_I_Root, 
     &                             Input_Opt, State_Met, State_Chm, 
     &                             State_Chm%Map_WetDep,        
     &                             State_Diag%Archive_BudgetWetDepFull, 
     &                             State_Diag%Archive_BudgetWetDepTrop, 
     &                             State_Diag%Archive_BudgetWetDepPBL,  
     &                             State_Diag%BudgetMass1, 
     &                             RC )
      ENDIF
#endif

#if defined( USE_TEND )
      !=================================================================
      ! Archive species concentrations for tendencies (ckeller,7/15/2015)
      !=================================================================
      CALL TEND_STAGE1( am_I_Root, Input_Opt, State_Met, 
     &                  State_Chm, 'WETD',    RC         )

      ! Trap potential errors
      IF ( RC /= GC_SUCCESS ) THEN
         ErrMsg = 'Error encountered in "TEND_STAGE1"!'
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
#endif

      !=================================================================
      ! Only do wet deposition for large-scale + anvil precip
      !=================================================================

#if defined( NC_DIAG )
      !------------------------------------------
      ! Zero diagnostic arrays in State_Diag
      ! at the start of a new wetdep cycle
      !
      ! *FracLs diagnostics might need work...
      !------------------------------------------
      IF ( State_Diag%Archive_PrecipFracLS ) 
     &                              State_Diag%PrecipFracLS = 0.0_f4
      IF ( State_Diag%Archive_RainFracLS   ) 
     &                              State_Diag%RainFracLS   = 0.0_f4
      IF ( State_Diag%Archive_WashFracLS   ) 
     &                              State_Diag%WashFracLS   = 0.0_f4
      IF ( State_Diag%Archive_WetLossLS    ) 
     &                              State_Diag%WetLossLS    = 0.0_f4
#endif


      !------------------------------------------
      ! Create precip fields
      ! (Always assume large-scale precip)
      !------------------------------------------
      CALL MAKE_QQ( am_I_Root = am_I_Root,
     &              State_Met = State_Met,
     &              LS        =.TRUE.,
     &              RC        = RC         )

       ! Trap potential errors
      IF ( RC /= GC_SUCCESS ) THEN
         ErrMsg = 'Error encountered in "Make_QQ"!'
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF

      ! Debug print
      IF ( doPrt ) CALL DEBUG_MSG( '### DO_WETDEP: before LS wetdep' )

      !-----------------------------------------
      ! Do wet deposition
      ! (Always assume large-scale precip)
      !-----------------------------------------
      CALL WETDEP( am_I_Root  = am_I_Root,
     &             Input_Opt  = Input_Opt, 
     &             State_Met  = State_Met,
     &             State_Chm  = State_Chm,
     &             State_Diag = State_Diag,
     &             LS         = .TRUE.,
     &             RC         = RC          )

      ! Trap potential errors
      IF ( RC /= GC_SUCCESS ) THEN
         ErrMsg = 'Error encountered in "Wetdep"!'
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF

      ! Debug print
      IF ( doPrt ) CALL DEBUG_MSG( '### DO_WETDEP: after LS wetdep' )

#if defined( USE_TEND )
      !=================================================================
      ! Calculate tendencies and write to diagnostics 
      ! (ckeller,7/15/2015)
      !=================================================================
      DT_TEND = GET_TS_DYN()
      CALL TEND_STAGE2( am_I_Root, Input_Opt, State_Met, 
     &                  State_Chm, 'WETD',    DT_TEND,   RC )

      ! Trap potential errors
      IF ( RC /= GC_SUCCESS ) THEN
         ErrMsg = 'Error encountered in "TEND_STAGE1"!'
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
#endif

#if defined( NC_DIAG )
      !----------------------------------------------------------
      ! Wet deposition budget diagnostics - Part 2 of 2
      !----------------------------------------------------------
      IF ( State_Diag%Archive_BudgetWetDep ) THEN
         ! Get final masses and compute diagnostics
         CALL Compute_Column_Mass( am_I_Root, 
     &                          Input_Opt, State_Met, State_Chm, 
     &                          State_Chm%Map_WetDep,        
     &                          State_Diag%Archive_BudgetWetDepFull, 
     &                          State_Diag%Archive_BudgetWetDepTrop, 
     &                          State_Diag%Archive_BudgetWetDepPBL,  
     &                          State_Diag%BudgetMass2, 
     &                          RC ) 
         DT_Dyn = Get_Ts_Dyn()
         CALL Compute_Budget_Diagnostics( am_I_Root, 
     &                          State_Chm%Map_WetDep,
     &                          DT_Dyn,                                
     &                          State_Diag%Archive_BudgetWetDepFull, 
     &                          State_Diag%Archive_BudgetWetDepTrop, 
     &                          State_Diag%Archive_BudgetWetDepPBL,  
     &                          State_Diag%BudgetWetDepFull,         
     &                          State_Diag%BudgetWetDepTrop,         
     &                          State_Diag%BudgetWetDepPBL,          
     &                          State_Diag%BudgetMass1,          
     &                          State_Diag%BudgetMass2,          
     &                          RC )
         IF ( RC /= GC_SUCCESS ) THEN
            ErrMsg = 'Wetdep budget diagnostics error 2'
            CALL GC_Error( ErrMsg, RC, ThisLoc )
            RETURN
         ENDIF
      ENDIF
#endif

      END SUBROUTINE DO_WETDEP
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: make_qq
!
! !DESCRIPTION: Subroutine MAKE\_QQ computes the large-scale or convective 
!  precipitation fields for use with WETDEP
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE MAKE_QQ( am_I_Root, State_Met, LS, RC )
!
! !USES:
!
      USE CMN_SIZE_mod
      USE ErrCode_Mod
      USE State_Met_Mod, ONLY : MetState
!
! !INPUT PARAMETERS: 
!
      LOGICAL,        INTENT(IN)  :: am_I_Root ! Are we on the root CPU
      TYPE(MetState), INTENT(IN)  :: State_Met ! Meteorology State object
      LOGICAL,        INTENT(IN)  :: LS        ! =T, denotes large scale precip
                                               ! =F, denotes convective precip
!
! !OUTPUT PARAMETERS:
!
      INTEGER,        INTENT(OUT) :: RC        ! Success or failure
! 
! !REMARKS:
!  Construct QQ and PDOWN directly from met fields.
!                                                                             .
!  This only applies to large-scale precip, as the #if defined
!  block in routine DO_WETDEP prevents the wet deposition
!  routines from being called if it is convective precip.
!                                                                             .
!  Met fields:
!  =================
!  DQRLSAN   = 3-D precip production rate  (LS+anvil) [kg/kg/s]
!  PFILSAN   = Dwnwd flux of ice precip    (LS+anvil) [kg/m2/s]
!  PFLLSAN   = Dwnwd flux of liquid precip (LS+anvil) [kg/m2/s]
!  REEVAPLS  = Evap of precip'ing LS+anvil condensate [kg/kg/s]
!                                                                             .
!  Unit conversion for QQ:
!  =======================
! 
!      kg H2O   |   m^3 H2O   | AIRDEN kg air       m^3 H2O
!   ------------+-------------+--------------- = -------------   
!    kg air * s | 1000 kg H2O |    m^3 air        m^3 air * s
! 
!  and [m^3 H2O/m3 air] = [cm^3 H2O/cm3 air] because the same conversion 
!  factor from m^3 -> cm^3 is in both the numerator and the denominator.
!                                                                             .
!  Unit conversion for PDOWN:
!  ==========================
!                                                                             .
!      kg H2O |   m^3 H2O   | 1e6 cm^3 |  m^2       
!   ----------+-------------+----------+--------- +
!     m^2 * s | 1000 kg H2O |   m^3    | 1e4 cm2 
!                                                                             .
!      kg ice |   m^3 ice   | 1e6 cm^3 |  m^2
!   ----------+-------------+----------+---------
!     m^2 * s |  917 kg ice |   m^3    | 1e4 cm2 
!                                                                             .
!  = [ (PFILSAN/1000) * 100 ] + [ (PFILSAN/1000) * 100]
!
!
! !REMARKS:
!  The PFILSAN and PFLLSAN fields are defined on level edges.
!  Therefore, we must use L+1 to index them.
!
! !REVISION HISTORY: 
!  29 Feb 2000 - H. Liu, R. Yantosca - Initial version
!  (1 ) Now we partition MOISTQ into large-scale and convective parts, using
!        total precipitation PREACC and convective precipitation PRECON (both
!        are vertical integral amounts). The precipitation field at altitudes
!        (PDOWN) is also made (hyl, djj, 10/17/98).
!  (2 ) MAKE_QQ is written in Fixed-Form Fortran 90. (bmy, 4/2/99)!
!  (3 ) AIRDEN, MOISTQ, QQ, and PDOWN are dimensioned (LLPAR,IIPAR,JJPAR) 
!       in order to maximize loop efficiency when processing an (I,J) 
!       column layer by layer. (bmy, 3/14/00)
!  (4 ) MOISTQ is originally [g H2O/kg air/day], and is converted in
!        READ_A6 to [kg H2O/kg air/s]. (bmy, 3/14/00)
!  (5 ) Now reference PREACC, PRECON from "dao_mod.f" instead of from
!        common block header file "CMN_PRECIP" (bmy, 6/26/00)
!  (6 ) Now pass BXHEIGHT as an argument.  Also added to "dao_mod.f". 
!        (bmy, 6/26/00)
!  (7 ) Moved from "dao_mod.f" to "wetscav_mod.f".  Also made PREACC
!        and PRECON into arguments. (bmy, 10/12/00)
!  (8 ) Updated comments (bmy, 9/4/01)
!  (9 ) BXHEIGHT is now sized (IIPAR,JJPAR,LLPAR) (bmy, 10/4/01)
!  (10) Removed obsolete, commented-out code from 10/01 (bmy, 11/26/01)
!  (11) Now reference met field arrays directly from "dao_mod.f" (bmy, 11/8/02)
!  16 Sep 2010 - R. Yantosca - Added ProTeX headers
!  16 Sep 2010 - R. Yantosca - Compute QQ and PDOWN from MERRA met fields
!  09 Feb 2012 - R. Yantosca - Added modifications for GEOS-5.7.x met fields
!  09 Nov 2012 - M. Payer    - Replaced all met field arrays with State_Met
!                              derived type object
!  26 Sep 2013 - R. Yantosca - Renamed GEOS_57 Cpp switch to GEOS_FP
!  06 Nov 2014 - R. Yantosca - Now use State_Met%AIRDEN(I,J,L)
!  06 Nov 2014 - R. Yantosca - Now use State_Met%MOISTQ(I,J,L)
!  29 Apr 2015 - E. Lundgren - Now use State_Met%MAIRDEN instead of AIRDEN
!                              since AIRDEN is now dry air density and MAIRDEN
!                              is moist air density, needed for use with MOISTQ
!  12 Aug 2015 - R. Yantosca - Add support for MERRA2 meteorology
!  09 Nov 2017 - R. Yantosca - Move allocations to INIT_WETSCAV
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER :: I, J, L

      !=================================================================
      ! Loop over surface grid boxes
      !=================================================================

      ! Initialize
      RC = GC_SUCCESS

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L )
      DO J = 1, JJPAR
      DO I = 1, IIPAR
      DO L = 1, LLPAR

         ! Rate of new precipitation formation in grid box (I,J,L)
         ! [cm3 H2O/cm3 air/s]
         QQ(L,I,J)     = ( State_Met%DQRLSAN(I,J,L)                 )
     &                 * ( State_Met%MAIRDEN(I,J,L)     / 1000.0_fp )

         ! Rate of re-evaporation in grid box (I,J,L)
         ! [cm3 H2O/cm3 air/s]
         REEVAP(L,I,J) = ( State_Met%REEVAPLS(I,J,L)                )
     &                 * ( State_Met%AIRDEN(I,J,L)      / 1000.0_fp )


            ! Column precipitation [cm3 H2O/cm2 air/s]
         PDOWN(L,I,J)  = ( ( State_Met%PFLLSAN(I,J,L+1) / 1000.0_fp )
     &                 +   ( State_Met%PFILSAN(I,J,L+1) /  917.0_fp ) )
     &                 * 100.0_fp

      ENDDO
      ENDDO  
      ENDDO
!$OMP END PARALLEL DO

      END SUBROUTINE MAKE_QQ
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: e_ice
!
! !DESCRIPTION: Subroutine E\_ICE computes Eice(T), the saturation vapor 
!  pressure of ice at a given Celsius temperature. 
!\\
!\\
! !INTERFACE:
!
      FUNCTION E_ICE( TK ) RESULT( VALUE )
!
! !INPUT PARAMETERS: 
!
      REAL(fp), INTENT(IN) :: TK      ! Temperature [K] 
!
! !RETURN VALUE:
!
      REAL(fp)             :: VALUE   ! Saturation vapor pressure [hPa]
!
! !REMARKS:
!  Marti & Mauersberber (GRL '93) formulation of saturation 
!  vapor pressure of ice [Pa] is: log P = A/TK + B
! 
! !REVISION HISTORY: 
!  08 Feb 2005 - R. Yantosca - Initial version
!  (1 ) Now use the same analytic function as the Goddard CTM (bmy, 2/8/05)
!  16 Sep 2010 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
      REAL(fp), PARAMETER  :: A = -2663.5e+0_fp
      REAL(fp), PARAMETER  :: B =  12.537e+0_fp

      !=================================================================
      ! E_ICE begins here!
      !=================================================================
      
      ! Saturation vap press of Ice [Pa] -- divide by 100 for [hPa]
      IF ( TK <= TINY_FP ) THEN
         VALUE = 0.0_fp
      ELSE
         VALUE = ( 10e+0_fp**( A/TK + B ) ) / 100e+0_fp 
      ENDIF

      END FUNCTION E_ICE
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: compute_l2g
!
! !DESCRIPTION: Subroutine COMPUTE\_L2G computes the ratio L2G = Cliq / Cgas, 
!  which is the mixing ratio of species in the liquid phase, divided by the 
!  mixing ratio of species in the gas phase.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE COMPUTE_L2G( K0, CR, pKa, TK, H2OLIQ, L2G )
!
! !USES:
!
      USE Henry_Mod, ONLY : Calc_KH
      USE Henry_Mod, ONLY : Calc_Heff
!     
! !INPUT PARAMETERS: 
!
      REAL(f8), INTENT(IN)  :: K0     ! Henry's solubility constant [M/atm]
      REAL(f8), INTENT(IN)  :: CR     ! Henry's volatility constant [K] 
      REAL(f8), INTENT(IN)  :: pKa    ! Henry's pH correction factor [1]
      REAL(fp), INTENT(IN)  :: TK     ! Temperature [K]
      REAL(fp), INTENT(IN)  :: H2OLIQ ! Liquid water content [cm3 H2O/cm3 air]
!                                     
! !OUTPUT PARAMETERS:                 
!                                     
      REAL(fp), INTENT(OUT) :: L2G    ! Cliq/Cgas ratio [1]
!
! !REMARKS:
!  The ratio Cliq / Cgas is obtained via Henry's law.  The appropriate 
!  values of Kstar298 and H298_R must be supplied for each species.  
!  (cf Jacob et al 2000, p. 3)
! 
! !REVISION HISTORY: 
!  23 Feb 2000 - R. Yantosca - Initial version
!  (1 ) Bundled into "wetscav_mod.f" (bmy, 11/8/02)
!  16 Sep 2010 - R. Yantosca - Added ProTeX headers
!  10-Jan-2011 - H. Amos - Corrected the units on KStar298 from moles/atm
!                          to M/atm
!  15-May-2013 - F. Paulot - Fix R constant
!  08 Dec 2015 - R. Yantosca - Now use functions from henry_mod.F
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER  :: RC
      REAL(f8) :: HEFF, KH, pH, TK_8

      !=================================================================
      ! COMPUTE_L2G begins here!
      !=================================================================

      ! Cast temperature to REAL*8
      TK_8 = TK

      ! For wetdep, we assume a pH of 4.5 for rainwater
      pH = 4.5_f8

      ! Calculate the Henry's law constant
      CALL CALC_KH( K0, CR, TK_8, KH, RC )

      ! Calculate effective Henry's law constant, corrected for pH
      ! (for those species that have a defined pKa value)
      CALL CALC_HEFF( pKa, pH, KH, HEFF, RC )

      ! Use Henry's Law to get the ratio:
      ! [ mixing ratio in liquid phase / mixing ratio in gas phase ]
      L2G   = HEFF * H2OLIQ

      END SUBROUTINE COMPUTE_L2G
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: compute_f 
!
! !DESCRIPTION: Subroutine COMPUTE\_F computes F, the fraction of soluble 
!  species lost by scavenging in convective cloud updrafts.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE COMPUTE_F( am_I_Root, N,         F,         ISOL, 
     &                      Input_Opt, State_Met, State_Chm, RC    )
!
! !USES:
!
      USE CMN_Size_Mod
      USE ErrCode_Mod
      USE Error_Mod
      USE Input_Opt_Mod,      ONLY : OptInput
      USE Species_Mod,        ONLY : Species
      USE State_Chm_Mod,      ONLY : ChmState
      USE State_Met_Mod,      ONLY : MetState
#if defined( TOMAS )
      USE Tomas_Mod,          ONLY : GetFraction
#endif
!
! !INPUT PARAMETERS: 
!
      LOGICAL,        INTENT(IN)    :: am_I_Root  ! Are we on the root CPU?
      INTEGER,        INTENT(IN)    :: N          ! Species ID
      TYPE(OptInput), INTENT(IN)    :: Input_Opt  ! Input Options object
      TYPE(MetState), INTENT(IN)    :: State_Met  ! Met State object
!
! !INPUT/OUTPUT PARAMETERS:
!
      TYPE(ChmState), INTENT(INOUT) :: State_Chm  ! Chemistry State
!
! !OUTPUT PARAMETERS:
!
      INTEGER,        INTENT(OUT)   :: ISOL       ! Index for ND38 diag
      REAL(fp),       INTENT(OUT)   :: F(:,:,:)   ! Soluble fraction of species
      INTEGER,        INTENT(OUT)   :: RC         ! Success or failure?
!
! !REVISION HISTORY: 
!  23 Feb 2000 - H. Liu, R. Yantosca - Initial version
!  (1 ) Currently works computes scavenging fractions for either full
!        chemistry simulation (NSRCX == 3) or Rn-Pb-Be chemistry simulation
!        (NSRCX == 1).  Set the scavenging fraction to zero for other
!        simulations which do not carry soluble tracers. (bmy, 3/2/00)
!  (2 ) Need to call INIT_SCAV to initialize the Vud, C_H2O, CLDLIQ, 
!        and CLDICE fields once per timestep. (bmy, 2/23/00)
!  (3 ) For aerosols only: now apply Eq. 2 for all temperatures.  Also
!        use the distance between the grid box centers in Eq. 2.  Updated
!        comments and made some cosmetic changes (hyl, bmy, 6/18/01)
!  (4 ) Remove IREF, JREF -- these are obsolete.  T is now dimensioned
!        (IIPAR,JJPAR,LLPAR).  T(IREF,JREF,L) is now T(I,J,L). (bmy, 9/27/01)
!  (5 ) Removed obsolete code from 9/01 (bmy, 10/23/01)
!  (6 ) Fix 2 bugs for aerosol scavenging in Rn-Pb-Be simulation: 
!        (a) set F(:,:,1) = 0 since we don't do any scavenging there.  
!        (b) DO L = 2, LLPAR to avoid any subscript range out of bounds 
!        errors (rjp, hyl, bmy, 1/10/02)
!  (7 ) Now set F=0 in the first level for all tracers.  Also now
!        compute the distance between grid box centers and use that in 
!        in Eq. 10 from Jacob et al, 2000 to compute F. (hyl, bmy, 1/24/02)
!  (8 ) Eliminated obsolete code from 1/02 (bmy, 2/27/02)
!  (9 ) Now reference T from "dao_mod.f" instead of from "CMN".  Also reference
!        BXHEIGHT from "dao_mod.f" instead of from "CMN_NOX".  Now bundled
!        into "wetscav_mod.f".  Now references IDTHNO3, IDTH2O2, etc, from
!        F90 module "tracerid_mod.f".  Added internal routines F_AEROSOL
!        and GET_ISOL.  Rewritten so that we don't duplicate code for 
!        different chemistry simulations. (bmy, 1/17/03)
!  (10) Now compute F for SO2 in the same way for both fullchem and offline 
!        simulations (rjp, bmy, 3/23/03)
!  (11) Added slots for carbon aerosol & dust tracers.  Now modified internal
!        routine GET_ISOL so it's not hardwired anymore. (rjp, bmy, 4/5/04)
!  (12) Added slots for sea salt aerosol tracers (rjp, bec, bmy, 4/20/04)
!  (13) Added slots for secondary organic aerosol tracers (rjp, bmy, 7/13/04)
!  (14) Remove reference to CMN, it's not needed.  Made internal routine
!        F_AEROSOL a module procedure rather than an internal routine to
!        COMPUTE_F in order to facilitate parallelization on the Altix.  Also
!        now pass all arguments explicitly to F_AEROSOL. (bmy, 7/20/04)
!  (15) Now wet scavenge mercury aerosol tracers (eck, bmy, 12/9/04)
!  (16) Updated for AS, AHS, LET, NH4aq, SO4aq.  Also condensed the IF
!        statement by combining branches for aerosols. (cas, bmy, 12/20/04)
!  (17) Updated for SO4s, NITs (bec, bmy, 4/25/05)
!  (18) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (19) Bug fix: Now do not over-deplete H2O2s.  Also change Henry's law
!        constant for Hg2 to 1.0d+14. Now use functions IS_Hg2 and IS_HgP to 
!        determine if a tracer is an Hg2 or HgP tagged tracer. 
!        (dkh, rjp, eck, cdh, bmy, 1/6/06)
!  (20) Updated for SOG4 and SOA4 (dkh, bmy, 5/18/06)
!  (21) Bug fix: now use separate conversion factors for H2O2 and NH3.
!        (havala, bmy, 7/26/06)
!  16 Sep 2010 - R. Yantosca - Added ProTeX headers
!  10-Jan-2011 - H.Amos      - Changed Hg2 Henry's law constant from 1.0d14 (no 
!                              citation) to 1.4d6 M/atm (HgCl2, Lindqvist &
!                              Rhode, 1985). Henry's law constant in
!                              wetscav_mod.f is now consistent with what's used
!                              in mercury_mod.f
!  27 Sep 2011 - H. Amos     - remove LHg_WETDasHNO3 logical, it's obsolete
!  09 Nov 2012 - M. Payer    - Replaced all met field arrays with State_Met
!                              derived type object
!  06 Mar 2013 - H. Amos     - merge C. Friedman's POP code
!  31 May 2013 - R. Yantosca - Now accept State_Chm, and pass it to TOMAS code
!  13 Aug 2013 - M. Sulprizio- Add modifications for updated SOA and SOA + 
!                              semivolatile POA simulations (H. Pye)
!  12 Sep 2013 - M. Sulprizio- Add modifications for acid uptake on dust
!                              aerosols (T.D. Fairlie)
!  18 Apr 2014 - R. Yantosca - Now make the F argument an assumed-shape array
!                              so that we can pass a pointer array slice to it
!  25 Jun 2014 - R. Yantosca - Now accept Input_Opt via the arg list
!  19 Dec 2014 - M. Sulprizio- Add bug fixes for CONV_H2O2 and CONV_NH3 from
!                              Duncan Fairlie. The square root term should be
!                              inverted to be consistent with Mari et al (2000).
!  04 Jun 2015 - E. Lundgren - Now accept am_I_Root and RC as arguments
!  23 Jun 2015 - M. Sulprizio- Add impaction scavenging for hydrophobic BC and
!                              homogeneous IN removal from Qiaoqiao Wang
!  04 Sep 2015 - R. Yantosca - Use species database to radically simplify
!                              the structure of this routine.
!  24 Sep 2015 - E. Lundgren - Convert kg/kg total air <-> kg for TOMAS
!  25 Sep 2015 - R. Yantosca - Move the computation of the Ki rate (Eq. 1 from
!                              Jacob et al 2000) into routine COMPUTE_Ki
!  25 Sep 2015 - R. Yantosca - Remove T pointer, it's not needed here
!  30 Sep 2015 - R. Yantosca - Now use ThisSpc%WD_Is_HNO3 and ThisSpc%WD_Is_SO2
!                              to flag the special cases of HNO3 & SO2 wetdep
!  05 Oct 2015 - R. Yantosca - Need to make ThisSpc !$OMP THREADPRIVATE
!  06 Oct 2015 - R. Yantosca - Add missing variables for TOMAS
!  16 Oct 2015 - E. Lundgren - Consolidate remaining embedded ifelse blocks
!  01 Jul 2016 - R. Yantosca - Now rename species DB object ThisSpc to SpcInfo
!  07 Jul 2016 - R. Yantosca - Now set ISOL = -1 if not a wetdep species
!  07 Jul 2016 - R. Yantosca - Remove pointer nullification from declarations
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      ! Scalars
      INTEGER                  :: I, J, L
      REAL(fp)                 :: SO2LOSS,   Ki
#if defined( TOMAS )
      REAL(fp)                 :: SOLFRAC,   XFRAC
#endif

      ! Arrays
      REAL(fp)                 :: KcScale(3)

      ! Pointers
      REAL(fp),      POINTER   :: p_C_H2O
      REAL(fp),      POINTER   :: p_CLDICE
      REAL(fp),      POINTER   :: p_CLDLIQ
      REAL(fp),      POINTER   :: p_T
      REAL(fp),      POINTER   :: H2O2s(:,:,:)
      REAL(fp),      POINTER   :: SO2s(:,:,:)

      ! Objects
      TYPE(Species), POINTER   :: SpcInfo
!
! !DEFINED PARAMETERS:
!
      ! Kc is the conversion rate from cloud condensate to precip [s^-1]
      REAL(fp),      PARAMETER :: Kc = 5e-3_fp

      !=================================================================
      ! COMPUTE_F begins here!
      !=================================================================

      ! Assume success
      RC         =  GC_SUCCESS

      ! Initialize
      F          =  0.0_fp
      p_C_H2O    => NULL()
      p_CLDICE   => NULL()
      p_CLDLIQ   => NULL()
      p_T        => NULL()
      H2O2s      => State_Chm%H2O2AfterChem
      SO2s       => State_Chm%SO2AfterChem
      SpcInfo    => State_Chm%SpcData(N)%Info

      ! ISOL is the wetdep ID (will be -999 if not a wetdep species)
      ISOL       =  SpcInfo%WetDepId

      ! Exit with F=0, ISOL=0 if this is not a wetdep species
      IF ( .not. SpcInfo%Is_WetDep ) THEN
         SpcInfo => NULL()
         RETURN
      ENDIF

      ! Temperature-dependent scale factors for the KC rate
      ! (conversion of cloud condensate -> precipitation)
      KcScale    = SpcInfo%WD_KcScaleFac

      !=================================================================
      ! %%% SPECIAL CASE %%%
      ! SO2 scavenges like an aerosol although it is considered
      ! to be a gas-phase species elsewhere (e.g. dry deposition)
      !=================================================================
      IF ( SpcInfo%WD_Is_SO2 ) THEN

         ! Compute fraction of SO2 scavenged
         CALL F_AEROSOL( KC, KcScale, State_Met, Input_Opt, F )

         !--------------------------------------------------------------
         ! Coupled full chemistry/aerosol simulation:
         ! Use the wet scavenging formula of Chin et al [1996], 
         ! such that a soluble fraction of SO2 is limited by the
         ! availability of H2O2 in the precipitating grid box. 
         ! Scavenge the soluble SO2 at the same rate as the sulfate.
         ! Update H2O2_sav and SO2_sav for use in RAINOUT, WASHOUT
         !--------------------------------------------------------------
         DO L = 2, LLPAR
         DO J = 1, JJPAR
         DO I = 1, IIPAR

            ! Make sure to deplete H2O2s the same as SO2s. 
            ! (dkh, rjp, bmy, 11/17/05)
            IF ( SO2s(I,J,L) > TINY_FP ) THEN
       
               ! Limit F
               SO2LOSS      = MIN( H2O2s(I,J,L), SO2s(I,J,L) )
               F(I,J,L)     = F(I,J,L) * SO2LOSS / SO2s(I,J,L)
               F(I,J,L)     = MAX(F(I,J,L), 0e+0_fp)
        
               ! Update saved H2O2 concentration
               H2O2s(I,J,L) = H2O2s(I,J,L) - ( SO2s(I,J,L) * F(I,J,L) )
               H2O2s(I,J,L) = MAX( H2O2s(I,J,L), TINY_FP )
        
            ELSE

               ! Set F = 0 if SO2s < EPSILON (dkh, rjp, bmy, 11/17/05) 
               F(I,J,L)     = 0e+0_fp
        
            ENDIF

            ! Update SO2
            SO2s(I,J,L)     = SO2s(I,J,L) * ( 1e+0_fp - F(I,J,L) )
            SO2s(I,J,L)     = MAX( SO2s(I,J,L), TINY_FP )
              
         ENDDO
         ENDDO
         ENDDO

      !=================================================================
      ! For all other species, compute the fraction of species 
      ! scavenged in updrafts.
      !=================================================================

      !--------------------------------------------------------------
      ! Soluble gas-phase species
      !
      ! NOTE: HNO3 scavenges like an aerosol, although it is
      ! considered a gas-phase species elsewhere.  Compute the
      ! fraction of HNO3 scavenged out of the column further down
      ! in the last ELSE block.
      !--------------------------------------------------------------
      ELSE IF ( SpcInfo%Is_Gas .and. ( .not. SpcInfo%WD_Is_HNO3 ) ) THEN

         ! No scavenging at the surface
         F(:,:,1) = 0.0_fp

         ! Start scavenging at level 2
         DO L = 2, LLPAR
         DO J = 1, JJPAR
         DO I = 1, IIPAR

            ! Set pointers
            p_C_H2O  => C_H2O (I,J,L)
            p_CLDICE => CLDICE(I,J,L)
            p_CLDLIQ => CLDLIQ(I,J,L)
            p_T      => State_Met%T(I,J,L)

            ! Compute Ki, the loss rate of a gas-phase species from
            ! the convective updraft (Eq. 1, Jacob et al, 2000)
            CALL COMPUTE_Ki( SpcInfo,  p_C_H2O, p_CLDICE,
     &                       p_CLDLIQ, Kc,      p_T,      Ki )

            ! Free pointers
            p_C_H2O  => NULL()
            p_CLDICE => NULL()
            p_CLDLIQ => NULL()
            p_T      => NULL()

            ! Compute F, the fraction of scavenged H2O2.
            ! (Eq. 2, Jacob et al, 2000)
            F(I,J,L) = GET_F( State_Met, Input_Opt, I, J, L, Ki )

         ENDDO
         ENDDO
         ENDDO
         
      !-----------------------------------------------------------
      ! Size-resolved soluble aerosol species
      ! (Microphysics simulations only)
      !-----------------------------------------------------------
      ELSE IF ( SpcInfo%MP_SizeResAer ) THEN

#if defined( TOMAS )
         ! Get the fraction of species scavenged in updrafts
         ! NOTE: The surface layer F(:,:,1) will be returned
         ! as zero, to shut off scavenging at the surface
         CALL F_AEROSOL( KC, KcScale, State_Met, Input_Opt, F )

         ! Adjust F for size-resolved aerosol (multiply by XFRAC)
         DO L = 2, LLPAR
         DO J = 1, JJPAR
         DO I = 1, IIPAR
            CALL GETFRACTION( I, J, L, N, .FALSE., State_Met,
     &                        XFRAC,      SOLFRAC, State_Chm )
            F(I,J,L) = XFRAC * F(I,J,L)
         ENDDO
         ENDDO
         ENDDO

#endif

      !-----------------------------------------------------------
      ! Size-resolved aerosol number
      ! (Microphysics simulations only)
      !-----------------------------------------------------------
      ELSE IF ( SpcInfo%MP_SizeResNum ) THEN

#if defined( TOMAS )
         ! Get the fraction of species scavenged in updrafts
         ! NOTE: The surface layer F(:,:,1) will be returned
         ! as zero, to shut off scavenging at the surface
         CALL F_AEROSOL( KC, KcScale, State_Met, Input_Opt, F )

         ! Adjust F for size-resolved aerosol number
         ! (multiply by XFRAC * SOLFRAC)
         DO L = 2, LLPAR
         DO J = 1, JJPAR
         DO I = 1, IIPAR
            CALL GETFRACTION( I, J, L, N, .FALSE., State_Met,
     &                              XFRAC,      SOLFRAC, State_Chm )
            F(I,J,L) = XFRAC * SOLFRAC * F(I,J,L)
         ENDDO
         ENDDO
         ENDDO

#endif

      !-----------------------------------------------------------
      ! Soluble aerosol species (non-size-resolved)
      ! including the special case of HNO3 as well as H2SO4
      ! if using TOMAS microphysics
      !-----------------------------------------------------------
      ELSE

         ! Get the fraction of species scavenged in updrafts
         ! NOTE: The surface layer F(:,:,1) will be returned
         ! as zero, to shut off scavenging at the surface
         CALL F_AEROSOL( KC, KcScale, State_Met, Input_Opt, F )

         ! Multiply by the aerosol scavenging efficiency
         ! For most species this is 1.0
         ! For SOA species this is usually 0.8
         IF ( SpcInfo%WD_AerScavEff > 0.0_fp ) THEN
            F = F * SpcInfo%WD_AerScavEff
         ENDIF

      ENDIF

      ! Nullify pointers
      p_C_H2O   => NULL()
      p_CLDICE  => NULL()
      p_CLDLIQ  => NULL()
      p_T       => NULL()
      H2O2s     => NULL()
      SO2s      => NULL()
      SpcInfo   => NULL()

      END SUBROUTINE COMPUTE_F
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: compute_ki
!
! !DESCRIPTION: Subroutine COMPUTE\_Ki computes the loss of species
!  by scavenging according to Jacob et al 2000, eq. 1.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE COMPUTE_Ki( SpcInfo, C_H2O, CLDICE, CLDLIQ, Kc, T, Ki )
!
! !USES:
!
      USE Species_Mod, ONLY : Species
!
! !INPUT PARAMETERS: 
!
      TYPE(Species), INTENT(IN)  :: SpcInfo  ! Species database object
      REAL(fp),      INTENT(IN)  :: C_H2O    ! Mixing ratio of H2O [v/v]
      REAL(fp),      INTENT(IN)  :: CLDICE   ! Cloud ice mixing ratio
                                             !  [cm3 ice/cm3 air]
      REAL(fp),      INTENT(IN)  :: CLDLIQ   ! Cloud liquid water mix ratio
                                             !  [cm3 H2O/cm3 air]
      REAL(fp),      INTENT(IN)  :: Kc       ! Rate for conversion of cloud
                                             !  condensate -> precip [1/s]
      REAL(fp),      INTENT(IN)  :: T        ! Temperature [K]
!
! !OUTPUT PARAMETERS:
!
      REAL(fp),      INTENT(OUT) :: Ki       ! Loss of species from updraft
                                             !  (cf Eq. 1, Jacob et al, 2000)
!
! !REMARKS:
!  This routine centralizes computations that are used in routines
!  COMPUTE_F and RAINOUT.
!
! !REVISION HISTORY:
!  25 Sep 2015 - R. Yantosca - Initial version
!  08 Dec 2015 - R. Yantosca - Make K0, CR, pKa 8-byte variables for
!                              compatibility w/ CALC_KH, CALC_HEFF routines
!  01 Jul 2016 - R. Yantosca - Now rename species DB object ThisSpc to SpcInfo
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      ! Scalars
      REAL(fp) :: L2G, I2G, C_TOT, F_L, F_I
      REAL(f8) :: K0, CR, pKa
      
      !=================================================================
      ! COMPUTE_Ki begins here!
      !=================================================================

      ! Get Henry's law parameters
      K0   = SpcInfo%Henry_K0
      CR   = SpcInfo%Henry_CR
      pKa  = SpcInfo%Henry_pKa

      IF ( SpcInfo%WD_LiqAndGas ) THEN

         ! For species that consider ice and liquid phases
         ! in wet deposition, compute ice to gas ratio for by
         ! co-condensation (Eq. 9, Jacob et al, 2000)
         IF ( C_H2O > 0.0_fp ) THEN
            I2G = ( CLDICE / C_H2O ) * SpcInfo%WD_ConvFacI2G
         ELSE
            I2G = 0.0_fp
         ENDIF

      ELSE

         ! For all other species, set the ice/gas ratio to zero
         I2G = 0.0_fp

      ENDIF

      ! Compute liquid to gas ratio for using 
      ! the appropriate parameters for Henry's law
      ! (Eqs. 7, 8, and Table 1, Jacob et al, 2000)
      CALL COMPUTE_L2G( K0, CR, pKa, T, CLDLIQ, L2G )

      ! Fraction of species in liquid & ice phases
      ! (Eqs. 4, 5, 6, Jacob et al, 2000)
      C_TOT = 1.0_fp + L2G + I2G
      F_L   = L2G / C_TOT
      F_I   = I2G / C_TOT

      ! Compute the rate constant Ki for loss of species from
      ! convective updraft scavenging (Eq. 1, Jacob et al, 2000)
      IF ( T >= 268.0_fp ) THEN
         Ki = KC * ( F_L + F_I )

      ELSE IF ( T > 248.0_fp  .and. T < 268.0_fp ) THEN
         Ki = KC * ( ( SpcInfo%WD_RetFactor * F_L ) + F_I )

      ELSE
         Ki = KC * F_I

      ENDIF

      END SUBROUTINE COMPUTE_Ki
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: f_aerosol
!
! !DESCRIPTION: Subroutine F\_AEROSOL returns the fraction of aerosol 
!  scavenged in updrafts
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE F_AEROSOL( KC, KcScale, State_Met, Input_Opt, F )
!
! !USES:
!
      USE CMN_SIZE_MOD
      USE State_Met_Mod,      ONLY: MetState
      USE Input_Opt_Mod,      ONLY: OptInput
!
!
! !INPUT PARAMETERS: 
!
      REAL(fp),       INTENT(IN)  :: KC                   ! Cloud condensate to
                                                          !  precipitation rate 
                                                          !  [1/s]
      REAL(fp),       INTENT(IN)  :: KcScale(3)           ! Scale factors for Kc
                                                          !  for 3 temperature
                                                          !  regimes
      TYPE(MetState), INTENT(IN)  :: State_Met            ! Meteorology State
      TYPE(OptInput), INTENT(IN)  :: Input_Opt            ! Input Options object
!
! !OUTPUT PARAMETERS:
!
      REAL(fp),       INTENT(OUT) :: F(IIPAR,JJPAR,LLPAR) ! Fraction of aerosol
                                                          !  scavenged in 
                                                          !  convective updrafts
! 
! !REVISION HISTORY: 
!  07 Nov 2002 - R. Yantosca - Initial version
!  16 Sep 2010 - R. Yantosca - Added ProTeX headers
!  09 Nov 2012 - M. Payer    - Replaced all met field arrays with State_Met
!                              derived type object
!  29 May 2013 - R. Yantosca - Segregate TOMAS-specific code with #ifdefs
!  10 Nov 2014 - C. Keller   - Now also apply TINY check in ESMF environment.
!  23 Jun 2015 - M. Sulprizio- Add impaction scavenging for hydrophobic BC and
!                              homogeneous IN removal from Qiaoqiao Wang
!  25 Sep 2015 - R. Yantosca - Rewrite this routine to avoid testing on the
!                              tracer number.  Remove the N argument.
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      ! Scalars
      INTEGER  :: I, J, L
      REAL(fp) :: TMP, FF, Scaled_KC

      !=================================================================
      ! F_AEROSOL begins here!
      !
      ! Aerosol species are 100% in the cloud condensate phase, so 
      ! we set K = Kc, and compute F accordingly (cf Jacob et al 2000 )    
      !=================================================================

      ! Turn off scavenging in the first level by setting F = 0
      F(:,:,1) = 0.0_fp

      ! Apply scavenging in levels 2 and higher
      DO L = 2, LLPAR
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! Apply temperature-dependent scale factors to the KC rate for ..
         IF ( State_Met%T(I,J,L) < 237.0_fp ) THEN

            ! Ice: T < 237 K:
            Scaled_KC = KC * KcScale(1)

         ELSE IF ( State_Met%T(I,J,L) >= 237.0_fp  .and. 
     &             State_Met%T(I,J,L) <  258.0_fp ) THEN
            
            ! Snow: 237 K <= T < 258 K
            Scaled_KC = KC * KcScale(2)

         ELSE

            ! Rain: T > 258 K
            Scaled_KC = KC * KcScale(3) 

         ENDIF

         ! (Eq. 2, Jacob et al, 2000, with K = Kc)
         ! Kc now has been scaled for impaction scavenging (bmy, 9/24/15)
         F(I,J,L) = GET_F( State_Met, Input_Opt, I, J, L, Scaled_KC )
            
      ENDDO
      ENDDO
      ENDDO
      
      END SUBROUTINE F_AEROSOL 
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: rainout
!
! !DESCRIPTION: Subroutine RAINOUT computes RAINFRAC, the fraction of soluble 
!  species lost to rainout events in precipitation.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE RAINOUT( I,         J,         L,         N,
     &                    K_RAIN,    DT,        F,         RAINFRAC,
     &                    Input_Opt, State_Met, State_Chm, RC        )
!
! !USES:
!
      USE CMN_Size_Mod
      USE ErrCode_Mod
      USE Input_Opt_Mod, ONLY : OptInput
      USE Species_Mod,   ONLY : Species
      USE State_Chm_Mod, ONLY : ChmState
      USE State_Met_Mod, ONLY : MetState
!
! !INPUT PARAMETERS: 
!
      INTEGER,        INTENT(IN)  :: I          ! Longitude index
      INTEGER,        INTENT(IN)  :: J          ! Latitude index
      INTEGER,        INTENT(IN)  :: L          ! Level index
      INTEGER,        INTENT(IN)  :: N          ! Species number
      REAL(fp),       INTENT(IN)  :: K_RAIN     ! Rainout rate constant [1/s]
      REAL(fp),       INTENT(IN)  :: DT         ! Timestep for rainout event [s]
      REAL(fp),       INTENT(IN)  :: F          ! Fraction of grid box that is
                                                !  precipitating [unitless]
      TYPE(OptInput), INTENT(IN)  :: Input_Opt  ! Input options
      TYPE(MetState), INTENT(IN)  :: State_Met  ! Meteorology State object
!
! !INPUT/OUTPUT PARAMETERS:
!
      TYPE(ChmState), INTENT(IN)  :: State_Chm  ! Chemistry State object

!
! !OUTPUT PARAMETERS:
!
      REAL(fp),       INTENT(OUT) :: RAINFRAC   ! Fraction of species lost 
                                                !  to rainout [unitless]
      INTEGER,        INTENT(OUT) :: RC         ! Success or failure?
! 
! !REVISION HISTORY: 
!  28 Feb 2000 - R. Yantosca - Initial version
!  (1 ) Currently works for either full chemistry simulation (NSRCX == 3) 
!        or Rn-Pb-Be chemistry simulation (NSRCX == 1).  Other simulations
!        do not carry soluble tracer, so set RAINFRAC = 0. (bmy, 2/28/00)
!  (2 ) Need to call INIT_SCAV to initialize the Vud, C_H2O, CLDLIQ, 
!        and CLDICE fields once per dynamic timestep. (bmy, 2/28/00)
!  (3 ) K_RAIN, the rainout rate constant, and F, the areal fraction of the 
!        grid box undergoing precipitiation, are computed according to 
!        Giorgi & Chaimedes, as described in Jacob et al, 2000.
!  (4 ) Now no longer suppress scavenging of HNO3 and aerosol below 258K.
!        Updated comments, cosmetic changes.  Now set TK = T(I,J,L) since
!        T is now sized (IIPAR,JJPAR,LLPAR) in "CMN". (djj, hyl, bmy, 1/24/02)
!  (5 ) Eliminated obsolete code (bmy, 2/27/02)
!  (6 ) Now reference T from "dao_mod.f".  Updated comments.  Now bundled 
!        into "wetscav_mod.f". Now refererences "tracerid_mod.f".  Also 
!        removed reference to CMN since we don't need NSRCX. (bmy, 11/8/02)
!  (7 ) Now updated for carbon & dust aerosol tracers (rjp, bmy, 4/5/04)
!  (8 ) Now updated for seasalt aerosol tracers (rjp, bec, bmy, 4/20/04)
!  (9 ) Now updated for secondary aerosol tracers (rjp, bmy, 7/13/04)
!  (10) Now treat rainout of mercury aerosol tracers (eck, bmy, 12/9/04)
!  (11) Updated for AS, AHS, LET, NH4aq, SO4aq.  Also condensed the IF
!        statement by grouping blocks together. (cas, bmy, 12/20/04)
!  (12) Updated for SO4s, NITs (bec, bmy, 4/25/05)
!  (13) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (14) Change Henry's law constant for Hg2 to 1.0d+14.  Now use functions
!        IS_Hg2 and IS_HgP to determine if the tracer is a tagged Hg0 or
!        HgP tracer. (eck, cdh, bmy, 1/6/06)
!  (15) Updated for SOG4 and SOA4 (dkh, bmy, 5/18/06)
!  (16) For GEOS-5, suppress rainout when T < 258K (hyl, bmy, 3/5/08)
!  (17) Bug fix: need to use separate conversion parameters for H2O2 and
!        NH3.  This was the same fix as in COMPUTE_F but until now we had
!        overlooked this. (havala, bmy, 7/20/09)
!  25 Aug 2010 - R. Yantosca - Treat MERRA in the same way as GEOS-5
!  16 Sep 2010 - R. Yantosca - Added ProTeX headers
!  27 Sep 2011 - H. Amos     - remove LHg_WETDasHNO3 logical, it's obsolete
!  09 Feb 2012 - R. Yantosca - Treat GEOS-5.7.x in the same way as MERRA
!  09 Nov 2012 - M. Payer    - Replaced all met field arrays with State_Met
!                              derived type object
!  13 Aug 2013 - M. Sulprizio- Add modifications for updated SOA and SOA + 
!                              semivolatile POA simulations (H. Pye)
!  12 Sep 2013 - M. Sulprizio- Add modifications for acid uptake on dust
!                              aerosols (T.D. Fairlie)
!  26 Sep 2013 - R. Yantosca - Renamed GEOS_57 Cpp switch to GEOS_FP
!  25 Aug 2014 - M. Sulprizio- Now accept Input_Opt as an argument
!  19 Dec 2014 - M. Sulprizio- Add bug fixes for CONV_H2O2 and CONV_NH3 from
!                              Duncan Fairlie. The square root term should be
!                              inverted to be consistent with Mari et al (2000).
!  23 Jun 2015 - M. Sulprizio- Add impaction scavenging for hydrophobic BC and
!                              homogeneous IN removal from Qiaoqiao Wang
!  08 Jul 2015 - E. Lundgren - Add marine organic aerosols (B.Gantt, M.Johnson)
!  23 Sep 2015 - R. Yantosca - Use the species database to eliminate many
!                              IF blocks to simplify the code logic
!  25 Sep 2015 - R. Yantosca - Now call routine COMPUTE_Ki to compute the
!                              rate of scavenging (Eq. 1, Jacob et al 2000)
!  25 Sep 2015 - R. Yantosca - Now define THREADPRIVATE pointers
!  28 Sep 2015 - R. Yantosca - Needed to add special case for HNO3
!  30 Sep 2015 - R. Yantosca - Now use ThisSpc%WD_Is_HNO3 and ThisSpc%WD_Is_SO2
!                              to flag the special cases of HNO3 & SO2 wetdep
!  16 Oct 2015 - E. Lundgren - Use ThisSpc%WD_Is_H2SO4 to treat H2SO4 as
!                              aerosol for wetdep (for microphysics)
!  01 Jul 2016 - R. Yantosca - Now rename species DB object ThisSpc to SpcInfo
!  07 Jul 2016 - R. Yantosca - Now remove pointer threadprivate declarations
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      ! Scalars
      REAL(fp)               :: Ki, SO2LOSS

      ! Pointers
      REAL(fp),      POINTER :: p_C_H2O
      REAL(fp),      POINTER :: p_CLDICE
      REAL(fp),      POINTER :: p_CLDLIQ
      REAL(fp),      POINTER :: p_T
      REAL(fp),      POINTER :: H2O2s(:,:,:)
      REAL(fp),      POINTER :: SO2s(:,:,:)

      ! Objects
      TYPE(Species), POINTER :: SpcInfo

      !==================================================================
      ! RAINOUT begins here!
      !
      ! For aerosols, set Kc = K_RAIN and compute RAINFRAC according
      ! to Eq. 10 of Jacob et al 2000.  Call function GET_RAINFRAC.
      !==================================================================

      ! Initialize
      RC       =  GC_SUCCESS

      ! Set pointers
      p_C_H2O  => C_H2O(I,J,L)
      p_CLDICE => CLDICE(I,J,L)
      p_CLDLIQ => CLDLIQ(I,J,L)
      p_T      => State_Met%T(I,J,L)
      H2O2s    => State_Chm%H2O2AfterChem
      SO2s     => State_Chm%SO2AfterChem
      SpcInfo  => State_Chm%SpcData(N)%Info

      !=================================================================
      ! %%% SPECIAL CASE %%%
      ! SO2 scavenges like an aerosol although it is considered
      ! to be a gas-phase species elsewhere (e.g. dry deposition)
      !=================================================================
      IF ( SpcInfo%WD_Is_SO2 ) THEN
  
         ! Update SO2 and H2O2
         IF ( SO2s(I,J,L) > TINY_FP ) THEN 

            ! Treat SO2 as an aerosol
            RAINFRAC = GET_RAINFRAC( K_RAIN, F, DT )

            ! Apply temperature-dependent rainout efficiencies
            ! This accounts for impaction scavenging of certain aerosols
            CALL APPLY_RAINOUT_EFF( p_T, SpcInfo, RAINFRAC )

            ! Limit RAINFRAC 
            SO2LOSS      = MIN( SO2s(I,J,L), H2O2s(I,J,L) )
            RAINFRAC     = SO2LOSS * RAINFRAC / SO2s(I,J,L)
            RAINFRAC     = MAX( RAINFRAC, 0e+0_fp )
         
            ! Update saved H2O2 concentration
            H2O2s(I,J,L) = H2O2s(I,J,L) - ( SO2s(I,J,L) * RAINFRAC )
            H2O2s(I,J,L) = MAX( H2O2s(I,J,L), TINY_FP )
           
         ELSE

            ! If SO2s is not defined (i.e. if wetdep and convection
            ! are turned off), then set 
            RAINFRAC     = 0.0_fp
           
         ENDIF

         ! Update saved SO2 concentration
         SO2s(I,J,L)     = SO2s(I,J,L) * ( 1.0_fp - RAINFRAC )
         SO2s(I,J,L)     = MAX( SO2s(I,J,L), TINY_FP )

      !=================================================================
      ! Compute rainout fraction for soluble gas-phase species
      ! (except for HNO3 and H2SO4 which scavenge like aerosols)
      !=================================================================
      ELSE IF ( SpcInfo%Is_Gas                 .and. 
     &          ( .not. SpcInfo%WD_Is_HNO3 )   .and. 
     &          ( .not. SpcInfo%WD_Is_H2SO4 ) ) THEN  

         ! Compute Ki, the loss rate of a gas-phase species from
         ! the convective updraft (Eq. 1, Jacob et al, 2000)
         CALL COMPUTE_Ki( SpcInfo,  p_C_H2O, p_CLDICE,
     &                    p_CLDLIQ, K_RAIN,  p_T,      Ki )


         ! Compute RAINFRAC, the fraction of rained-out H2O2
         ! (Eq. 10, Jacob et al, 2000)
         RAINFRAC = GET_RAINFRAC( Ki, F, DT )

      !=================================================================
      ! Compute rainout fraction for aerosol species
      ! (including HNO3 and H2SO4 which scavenge like aerosols)
      !=================================================================
      ELSE

         ! Compute rainout fraction for aerosol tracres
         RAINFRAC = GET_RAINFRAC( K_RAIN, F, DT )

         ! Apply temperature-dependent rainout efficiencies
         ! This accounts for impaction scavenging of certain aerosols
         CALL APPLY_RAINOUT_EFF( p_T, SpcInfo, RAINFRAC )

      ENDIF

      ! Free pointers
      p_C_H2O  => NULL()
      p_CLDICE => NULL()
      p_CLDLIQ => NULL()
      p_T      => NULL()
      H2O2s   => NULL()
      SO2s    => NULL()
      SpcInfo  => NULL()

      END SUBROUTINE RAINOUT
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: apply_rainout_eff
!
! !DESCRIPTION: Subroutine APPLY\_RAINOUT\_EFF multiplies the rainout fraction
!  computed by RAINOUT with the rainout efficiency for one of 3 temperature
!  ranges: (1) T < 237 K; (2) 237 K <= T < 258 K; (3) T > 258 K. The rainout
!  efficiencies for each aerosol species are defined in the species database
!  object (i.e. State\_Chm%SpcData(:)%Info).
!\\
!\\
!  This allows us to apply the impaction scavenging of certain aerosol species
!  (BC, dust, HNO3) as implemented by Qiaoqiao Wang, while also suppressing 
!  rainout for other aerosol species.  The prior code achieved this by using
!  a large and confusing IF statement, whose logic was hard to understand.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE APPLY_RAINOUT_EFF( TK, SpcInfo, RainFrac )
!
! !USES:
!
      USE Species_Mod, ONLY : Species
!
! !INPUT PARAMETERS: 
!
      REAL(fp),      INTENT(IN)    :: TK         ! Temperature [K]
      TYPE(Species), INTENT(IN)    :: SpcInfo    ! Species Database object
!
! !INPUT/OUTPUT PARAMETERS: 
!
      REAL(fp),      INTENT(INOUT) :: RainFrac   ! Rainout fraction
!
! !REVISION HISTORY:
!  06 Jan 2015 - R. Yantosca - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:

      ! Apply temperature-dependent rainout efficiencies
      ! This accounts for impaction scavenging of certain aerosols
      IF ( TK < 237.0_fp ) THEN 

         ! Ice: T < 237 K
         RainFrac = RainFrac * SpcInfo%WD_RainoutEff(1)

      ELSE IF ( TK >= 237.0_fp .and. TK < 258.0_fp ) THEN

         ! Snow: 237 K <= T < 258 K
         RainFrac = RainFrac * SpcInfo%WD_RainoutEff(2)

      ELSE
            
         ! Liquid rain: T > 258 K
         RainFrac = RainFrac * SpcInfo%WD_RainoutEff(3)

      ENDIF

      END SUBROUTINE APPLY_RAINOUT_EFF
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_rainfrac
!
! !DESCRIPTION: Function GET\_RAINFRAC computes the fraction of species 
!  lost to rainout according to Jacob et al 2000.
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_RAINFRAC( K, F, DT ) RESULT( RAINFRAC )
!
! !INPUT PARAMETERS: 
!
      REAL(fp), INTENT(IN) :: K          ! Rainout rate constant [1/s]
      REAL(fp), INTENT(IN) :: F          ! Timestep for rainout event [s]
      REAL(fp), INTENT(IN) :: DT         ! Fraction of grid box that is
                                         !  undergoing precipitation [unitless]
!
! !RETURN VALUE:
!
      REAL(fp)             :: RAINFRAC   ! Fraction of species lost to rainout
!
! !REVISION HISTORY: 
!  08 Nov 2002 - R. Yantosca - Initial version
!  (1 ) Now move internal routines GET_RAINFRAC to the module and pass all 
!        arguments explicitly.  This facilitates parallelization on the 
!        Altix platform (bmy, 7/20/04) 
!  16 Sep 2010 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
      !=================================================================
      ! GET_RAINFRAC begins here!
      !=================================================================

      ! (Eq. 10, Jacob et al, 2000 ) 
      RAINFRAC = F * ( 1 - EXP( -K * DT ) )

      END FUNCTION GET_RAINFRAC
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: washout
!
! !DESCRIPTION: Subroutine WASHOUT computes WASHFRAC, the fraction of 
!  soluble species lost to washout events in precipitation.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE WASHOUT( am_I_Root, I,        J,         L,  
     &                    N,         BXHEIGHT, TK,        PP, 
     &                    DT,        F,        H2O2s,     SO2s, 
     &                    WASHFRAC,  KIN,      Input_Opt, State_Met, 
     &                    State_Chm, RC                              )
!
! !USES:
!
      USE CMN_SIZE_MOD
      USE ErrCode_Mod
      USE Input_Opt_Mod, ONLY : OptInput
      USE Species_Mod,   ONLY : Species
      USE State_Chm_Mod, ONLY : ChmState
      USE State_Met_Mod, ONLY : MetState
#if defined( TOMAS )
      USE ERROR_MOD
      USE TOMAS_MOD,     ONLY : IBINS,    ICOMP
      USE UnitConv_Mod
#endif
!
! !INPUT PARAMETERS: 
!
      LOGICAL,        INTENT(IN)    :: am_I_Root  ! Are we on the root CPU?
      INTEGER,        INTENT(IN)    :: I          ! Longitude index
      INTEGER,        INTENT(IN)    :: J          ! Latitude index
      INTEGER,        INTENT(IN)    :: L          ! Level index
      INTEGER,        INTENT(IN)    :: N          ! Species number
      REAL(fp),       INTENT(IN)    :: BXHEIGHT   ! Grid box height [m]
      REAL(fp),       INTENT(IN)    :: TK         ! Temperature [K]
      REAL(fp),       INTENT(IN)    :: PP         ! Precip rate thru bottom 
                                                  !  of grid (I,J,L)  
                                                  !  [cm3 H2O/cm2 air/s]
      REAL(fp),       INTENT(IN)    :: DT         ! Timestep [s]
      REAL(fp),       INTENT(IN)    :: F          ! Fraction of grid box that 
                                                  !  is precipitating [1]

      TYPE(OptInput), INTENT(IN)    :: Input_Opt  ! Input Options
      TYPE(MetState), INTENT(IN)    :: State_Met  ! Meteorology State

! !INPUT/OUTPUT PARAMETERS:
!
      TYPE(ChmState), INTENT(INOUT) :: State_Chm  ! Chemistry State object
      REAL(fp),       INTENT(INOUT) :: H2O2s      ! H2O2 [v/v] and SO2 [v/v] 
      REAL(fp),       INTENT(INOUT) :: SO2s       ! conc's after aqueous rxns 
                                                  ! are applied.  These are 
                                                  ! computed in the sulfate 
                                                  ! chemistry module and
                                                  ! passed here as arguments. 
!
! !OUTPUT PARAMETERS:
!
      REAL(fp),       INTENT(OUT)   :: WASHFRAC   ! Fraction of species lost
                                                  !  to washout [1]
      LOGICAL,        INTENT(OUT)   :: KIN        ! =T washout is a 
                                                  !    kinetic process
                                                  ! =F washout is an 
                                                  !    equilibrium process
      INTEGER,        INTENT(OUT)   :: RC         ! Success or failure?       
!
! !REVISION HISTORY: 
!  28 Feb 2000 - R. Yantosca - Initial version
!  (1 ) Currently works for either full chemistry simulation (NSRCX == 3) 
!        or Rn-Pb-Be chemistry simulation (NSRCX == 1).  Other simulations
!        do not carry soluble tracers, so set WASHFRAC = 0. 
!  (2 ) K_WASH, the rainout rate constant, and F, the areal fraction of the 
!        grid box undergoing precipitiation, are computed according to 
!        Giorgi & Chaimedes, as described in Jacob et al, 2000.
!  (3 ) Washout is only done for T >= 268 K, when the cloud condensate is
!        in the liquid phase. 
!  (4 ) T(I+I0,J+J0,L) is now T(I,J,L).  Removed IREF, JREF -- these are 
!        obsolete.  Updated comments. (bmy, 9/27/01)
!  (5 ) Removed obsolete commented out code from 9/01 (bmy, 10/24/01)
!  (6 ) Now reference BXHEIGHT, T from "dao_mod.f".  Also remove reference
!        to "CMN_NOX".  Updated comments.  Now bundled into "wetscav_mod.f".
!        Now also references "tracerid_mod.f".  Added internal routines
!        WASHFRAC_AEROSOL and WASHFRAC_LIQ_GAS.  Also removed reference to
!        CMN since we don't need to use NSRCX here. (bmy, 11/6/02)
!  (7 ) Updated for carbon aerosol and dust tracers (rjp, bmy, 4/5/04)
!  (8 ) Updated for seasalt aerosol tracers (rjp, bec, bmy, 4/20/04)
!  (9 ) Updated for secondary organic aerosol tracers (rjp, bmy, 7/13/04)
!  (10) Now move internal routines WASHFRAC_AEROSOL and WASHFRAC_LIQ_GAS
!        to the module and pass all arguments explicitly.  This facilitates
!        parallelization on the Altix platform (bmy, 7/20/04)
!  (11) Now handle washout of mercury aerosol tracers (eck, bmy, 12/9/04)
!  (13) Updated for AS, AHS, LET, NH4aq, SO4aq.  Also condensed the IF
!        statement by grouping blocks together (cas, bmy, 12/20/04)
!  (14) Updated for SO4s, NITs (bec, bmy, 4/25/05)
!  (15) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (16) Bug fix: Deplete H2O2s the same as SO2s.  Also change Henry's law
!        constant for Hg2 to 1.0d+14. Now use functions IS_Hg2 and IS_HgP to 
!        determine if a tracer is a tagged Hg0 or HgP tracer.
!        (dkh, rjp, eck, cdh, bmy, 1/6/06)
!  (17) Updated for SOG4 and SOA4 (bmy, 5/18/06)
!  16 Sep 2010 - R. Yantosca - Added ProTeX headers
!  30 Sep 2010 - H. Amos     - WASHFRAC_LIQ_GAS now a subroutine (was an 
!                              external function)
!  14 Oct 2010 - H. Amos     - Remove dependence on I, J. That means removing
!                              I, J as input arguments and adding T, BXHEIGHT,
!                              H2O2s, and SO4s and input arguments.
!  16 Aug 2011 - H. Amos     - move K_WASH to WASHFRAC_AEROSOL, WASHFRAC_HNO3,
!                              and WASHFRAC_LIQ_GAS 
!  16 Aug 2011 - H. Amos     - Replace logical AER with KIN. Serves the same 
!                              purpose in the code, but emphasizes to the user
!                              that the difference in washout isn't whether or
!                              not the tracer is an aerosol, it's whether or not
!                              washout is modeled as a kinetic vs equilibrium 
!                              process.
!  27 Sep 2011 - H. Amos     - remove LHg2_WETDasHNO3 logical, it's obsolete
!  20 Jan 2012 - H. Amos     - WASHFRAC for aerosol is now either computed by
!                              WASHFRAC_FINE_AEROSOL or WASHFRAC_COARSE_AEROSOL
!  31 May 2013 - R. Yantosca - Now accept State_Chm, and pass it to TOMAS code
!  13 Aug 2013 - M. Sulprizio- Add modifications for updated SOA and SOA + 
!                              semivolatile POA simulations (H. Pye)
!  12 Sep 2013 - M. Sulprizio- Add modifications for acid uptake on dust
!                              aerosols (T.D. Fairlie)
!  27 Sep 2013 - M. Sulprizio- DST2-DST4 are now considered coarse mode aerosol
!                              following the recommendation of T.D. Fairlie.
!                              NITd, SO4d, and DALK are similarly updated.
!  25 Aug 2014 - M. Sulprizio- Now accept Input_Opt as an argument
!  12 May 2015 - E. Lundgren - Change tracer units from v/v -> kg for TOMAS if
!                              not already in kg (ie. called from convection)
!  04 Jun 2015 - E. Lundgren - Now accept am_I_Root and RC as arguments
!  23 Jun 2015 - E. Lundgren - Adjust 5/12/15 TOMAS bug fix to convert
!                              kg/kg total air <-> kg for new convection units
!  08 Jul 2015 - E. Lundgren - Add marine organic aerosols (B.Gantt, M.Johnson)
!  22 Sep 2015 - E. Lundgren - Add kg/m2 <-> kg conversion for TOMAS
!  30 Sep 2015 - R. Yantosca - Now use ThisSpc%WD_Is_HNO3 and ThisSpc%WD_Is_SO2
!                              to flag the special cases of HNO3 & SO2 wetdep
!  16 Oct 2015 - E. Lundgren - Use ThisSpc%WD_Is_H2SO4 to treat H2SO4 as aerosol
!                              for wetdep, and consolidate ifelse blocks
!  27 Jul 2016 - E. Lundgren - Now expect dry mixing ratio not total
!  22 Dec 2016 - R. Yantosca - Bug fix: K0, CR, pKa have to be REAL(f8)
!  09 Nov 2017 - R. Yantosca - Return error condition to calling program
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      ! Scalars
#if defined( TOMAS )
      LOGICAL                :: UNITCHANGE_KGKG ! flag for STT units kg/kg
      LOGICAL                :: UNITCHANGE_KGM2 ! flag for STT units kg/m2 
#endif
      REAL(fp)               :: L2G, DZ, SO2LOSS
      REAL(f8)               :: K0,  CR, pKa

      ! Strings
      CHARACTER(LEN=255)     :: ErrMsg, ThisLoc

      ! Pointers
      TYPE(Species), POINTER :: SpcInfo

      !=================================================================
      ! WASHOUT begins here!
      !
      ! Call either WASHFRAC_FINE_AEROSOL, WASHFRAC_COARSE_AEROSOL,
      ! or WASHFRAC_LIQ_GAS to compute the fraction of species lost to 
      ! washout according to Jacob et al 2000
      !=================================================================

      ! Initialize
      RC      = GC_SUCCESS
      ErrMsg  = ''
      ThisLoc = ' -> at Washout (in module GeosCore/wetscav_mod.F)'

#if defined( TOMAS )
      !-----------------------------------------------------------------
      ! TOMAS MICROPHYSICS ONLY 
      !
      ! Convert species concentration units to [kg] if not already
      ! since TOMAS functions and routines expect [kg]. Units are
      ! kg/kg total air if WASHOUT is called from convection and are
      ! kg/m2 is called from DO_WASHOUT_ONLY. Since WASHOUT is called
      ! within an (I,J,L) loop, only convert units for a single grid
      ! box. Otherwise, run will take too long (ewl, 9/22/15)
      !-----------------------------------------------------------------
      UNITCHANGE_KGKG = .FALSE.
      UNITCHANGE_KGM2 = .FALSE.

      IF ( TRIM( State_Chm%Spc_Units ) .eq. 'kg/kg dry' ) THEN
         UNITCHANGE_KGKG = .TRUE.
         CALL ConvertBox_KgKgDry_to_Kg( am_I_Root, I, J, L,  
     &                                  State_Met, State_Chm, RC )

         ! Trap potential errors
         IF ( RC /= GC_SUCCESS ) THEN
            ErrMsg = 'Error encountered in "ConvertBox_KgKgDry_to_Kg"!'
            CALL GC_Error( ErrMsg, RC, ThisLoc )
            RETURN
         ENDIF

      ELSE IF ( TRIM( State_Chm%Spc_Units ) .eq. 'kg/m2' ) THEN
         UNITCHANGE_KGM2 = .TRUE.
         CALL ConvertBox_Kgm2_to_Kg( am_I_Root, I, J, L, 
     &                               State_Met, State_Chm, RC )

         ! Trap potential errors
         IF ( RC /= GC_SUCCESS ) THEN
            ErrMsg = 'Error encountered in "ConvertBox_KgM2_to_Kg"!'
            CALL GC_Error( ErrMsg, RC, ThisLoc )
            RETURN
         ENDIF

      ELSE

         ! Exit if units are not as expected
         ErrMsg = 'Incorrect initial species units:' // 
     &             TRIM( State_Chm%Spc_Units )
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN

      ENDIF
#endif

      ! DZ is the height of the grid box in cm
      DZ      =  BXHEIGHT * 1e+2_fp

      ! Get info about Nth species from the species database
      SpcInfo => State_Chm%SpcData(N)%Info
      
      !=================================================================
      ! %%% SPECIAL CASE %%%
      ! HNO3 scavenges like an aerosol although it is considered
      ! to be a gas-phase species elsewhere (e.g. dry deposition)
      !=================================================================
      IF ( SpcInfo%WD_Is_HNO3 ) THEN

         ! Washout is a kinetic process
         KIN      = .TRUE.

         ! Get washout fraction
         WASHFRAC = WASHFRAC_HNO3( DT, F, PP, TK )

      !=================================================================
      ! %%% SPECIAL CASE %%%
      ! SO2 scavenges like an aerosol although it is considered
      ! to be a gas-phase species elsewhere (e.g. dry deposition)
      !=================================================================
      ELSE IF ( SpcInfo%WD_Is_SO2 ) THEN

         ! NOTE: Even though SO2 is not an aerosol we treat it as SO4 in
         ! wet scavenging.  When evaporation occurs, it returns to SO4.
         KIN      = .TRUE.
         WASHFRAC = WASHFRAC_FINE_AEROSOL( DT, F, PP, TK )

         ! Use the wet-scavenging following [Chin et al, 1996] such 
         ! that a soluble fraction of SO2 is limited by the availability 
         ! of H2O2 in the precipitating grid box.  Then scavenge the 
         ! soluble SO2 at the same rate as sulfate.
         IF ( TK >= 268e+0_fp .AND. SO2s > TINY_FP ) THEN
         
            ! Adjust WASHFRAC
            SO2LOSS  = MIN( SO2s, H2O2s )
            WASHFRAC = SO2LOSS * WASHFRAC / SO2s
            WASHFRAC = MAX( WASHFRAC, 0e+0_fp )
                  
            ! Deplete H2O2s the same as SO2s (dkh, rjp, bmy, 11/17/05)
            H2O2s = H2O2s - ( SO2s * WASHFRAC )
            H2O2s = MAX( H2O2s, TINY_FP )

         ELSE
            WASHFRAC = 0e+0_fp
         
         ENDIF
         
         ! Update saved SO2 concentration 
         SO2s = SO2s * ( 1e+0_fp - WASHFRAC )
         SO2s = MAX( SO2s, TINY_FP ) 

      !=================================================================
      ! All other species
      !=================================================================

      !-----------------------------------------------------------------
      ! Washout for gas-phase species 
      ! (except H2SO4, NO3, and SO2, which scavenge like aerosols; 
      ! NO3 and SO2 are handled above and H2SO4 is handled in TOMAS 
      ! block further below)
      !-----------------------------------------------------------------
      ELSE IF ( SpcInfo%Is_Gas .and.
     &        ( .not. SpcInfo%WD_Is_H2SO4 ) ) THEN

         ! Get Henry's law parameters
         K0  = SpcInfo%Henry_K0
         CR  = SpcInfo%Henry_CR
         pKa = SpcInfo%Henry_pKa

         ! Washout is an equilibrium process
         KIN = .FALSE.
         
         ! Get the washout fraction for this species
         CALL WASHFRAC_LIQ_GAS( K0, CR, pKa, PP,       DT,
     &                          F,  DZ, TK,  WASHFRAC, KIN )

      !-----------------------------------------------------------------
      ! Washout for size-resolved aerosol species or
      ! size-resolved aerosol number (e.g. from TOMAS)
      ! NOTE: treat H2SO4 as an aerosol for wetdep in TOMAS
      !-----------------------------------------------------------------
      ELSE IF ( SpcInfo%MP_SizeResAer .or. SpcInfo%MP_SizeResNum ) THEN

#if defined( TOMAS )
         ! Washout is a kinetic process
         KIN      = .TRUE.

         ! Compute washout fraction
         CALL WASHFRAC_SIZE_AEROSOL( DT, F, PP, TK, N, I, J, L,
     &                             State_Met, State_Chm, WASHFRAC, RC )
#endif

      !-----------------------------------------------------------------
      ! Washout for coarse aerosol species (Reff >= 1um)
      !-----------------------------------------------------------------
      ELSE IF ( SpcInfo%WD_CoarseAer ) THEN
         
         ! Washout is a kinetic process
         KIN      = .TRUE.

         ! Compute washout fraction
         WASHFRAC = WASHFRAC_COARSE_AEROSOL( DT, F, PP, TK)

      !-----------------------------------------------------------------
      ! Washout for fine-aerosol species (Reff < 1um)
      !-----------------------------------------------------------------
      ELSE

         ! Washout is a kinetic process
         KIN      = .TRUE.

         ! Compute washout fraction
         WASHFRAC = WASHFRAC_FINE_AEROSOL( DT, F, PP, TK )

      ENDIF

#if defined( TOMAS )
      !-----------------------------------------------------------------
      ! TOMAS MICROPHYSICS ONLY 
      ! Convert State_Chm%Species units back to original units 
      ! if conversion occurred at start of WASHOUT (ewl, 5/12/15)
      !-----------------------------------------------------------------
      IF ( UNITCHANGE_KGKG ) THEN
         CALL ConvertBox_Kg_to_KgKgDry( am_I_Root, I, J, L, 
     &                                  State_Met, State_Chm, RC )

         ! Trap potential errors
         IF ( RC /= GC_SUCCESS ) THEN
            ErrMsg = 'Error encountered in "ConvertBox_Kg_to_KgKgDry"!'
            CALL GC_Error( ErrMsg, RC, ThisLoc )
            RETURN
         ENDIF

      ELSE IF ( UNITCHANGE_KGM2 ) THEN
         CALL ConvertBox_Kg_to_Kgm2( am_I_Root, I, J, L,
     &                               State_Met, State_Chm, RC )

         ! Trap potential errors
         IF ( RC /= GC_SUCCESS ) THEN
            ErrMsg = 'Error encountered in "ConvertBox_Kg_to_KgM2"!'
            CALL GC_Error( ErrMsg, RC, ThisLoc )
            RETURN
         ENDIF

      ENDIF

      ! Check that species units are as expected (ewl, 9/29/15)
      IF ( TRIM( State_Chm%Spc_Units ) /= 'kg/kg dry' .AND.
     &     TRIM( State_Chm%Spc_Units ) /= 'kg/m2' ) THEN
         ErrMsg = 'Incorrect final species units:' // 
     &             TRIM( State_Chm%Spc_Units )
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
#endif

      ! Free pointer
      SpcInfo => NULL()

      END SUBROUTINE WASHOUT
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: washfrac_fine_aerosol
!
! !DESCRIPTION: Function WASHFRAC\_FINE\_AEROSOL returns the fraction of
!  soluble aerosol species lost to washout.
!\\
!\\
! !INTERFACE:
!
      FUNCTION WASHFRAC_FINE_AEROSOL( DT, F, PP, TK ) 
     &         RESULT( WASHFRAC )
!
! !USES:
!
!
! !INPUT PARAMETERS: 
!
      REAL(fp), INTENT(IN) :: DT        ! Timestep of washout event [s]
      REAL(fp), INTENT(IN) :: F         ! Fraction of grid box that is
                                        !  precipitating [unitless]
      REAL(fp), INTENT(IN) :: PP        ! Precip rate thru bottom of grid
                                        !  box (I,J,L)  [cm3 H2O/cm2 air/s]
      REAL(fp), INTENT(IN) :: TK        ! Temperature in grid box [K]
!
! !RETURN VALUE:
!
      REAL(fp)             :: WASHFRAC  ! Fraction of soluble species
                                        !  lost to washout [1]
! 
! !REVISION HISTORY: 
!  08 Nov 2002 - R. Yantosca - Initial version
!  (1 ) WASHFRAC_AEROSOL used to be an internal function to subroutine WASHOUT.
!        This caused NaN's in the parallel loop on Altix, so we moved it to
!        the module and now pass Iall arguments explicitly (bmy, 7/20/04)
!  16 Sep 2010 - R. Yantosca - Added ProTeX headers
!  21 Jan 2011 - J. Fisher & Q. Wang - Update to account for time-dependent
!        shift in aerosol size distribution that slows washout as a rain
!        event proceeds (see e.g. Feng et al., 2007, 2009). 
!  16 Aug 2011 - H Amos      - Remove K_WASH from input list, make a defined 
!                              parameter.
!  20 Jan 2012 - H Amos      - rename WASHFRAC_FINE_AEROSOL to distinguish 
!                              this function from WASHFRAC_COARSE_AEROSOL
!  04 Sep 2013 - R. Yantosca - Bug fix: Prevent div-by-zero if F=0.  Because F 
!                              multiplies the whole expression for WASHFRAC,
!                              WASHFRAC=0 whenever F=0 anyway.
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETER:
!
      ! Washout rate constant for aerosols: aP^b (p: mm h^-1)
      ! K_WASH for aerosols in accumulation mode (qq,10/11/2011)
      REAL(fp), PARAMETER :: K_WASH = 1.06e-3_fp

      !=================================================================
      ! WASHFRAC_FINE_AEROSOL begins here!
      !=================================================================
      IF ( ( TK >= 268e+0_fp ) .OR. ITS_A_POPS_SIM ) THEN

         !---------------------------------
         ! T >= 268K (or POPS simulation)  
         !---------------------------------
         IF ( F > 0e+0_fp ) THEN
            WASHFRAC = F *(1e+0_fp - EXP( -K_WASH * 
     &              (PP / F*3.6e+4_fp )**0.61e+0_fp * DT / 3.6e+3_fp ))
         ELSE
            WASHFRAC = 0e+0_fp
         ENDIF

      ELSE

         !---------------------------------
         ! T < 268K
         !---------------------------------
         IF ( F > 0e+0_fp ) THEN 
            WASHFRAC = F *(1e+0_fp - EXP( -2.6e+1_fp*K_WASH *
     &              (PP / F*3.6e+4_fp )**0.96e+0_fp * DT / 3.6e+3_fp ))
         ELSE
            WASHFRAC = 0e+0_fp
         ENDIF

      ENDIF

      END FUNCTION WASHFRAC_FINE_AEROSOL
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: washfrac_coarse_aerosol
!
! !DESCRIPTION: Function WASHFRAC\_COARSE\_AEROSOL returns the fraction of 
!  soluble aerosol species lost to washout.
!\\
!\\
! !INTERFACE:
!
      FUNCTION WASHFRAC_COARSE_AEROSOL( DT, F, PP, TK ) 
     &         RESULT( WASHFRAC )
!
! !USES:
!
!
! !INPUT PARAMETERS: 
!
      REAL(fp), INTENT(IN) :: DT         ! Timestep of washout event [s]
      REAL(fp), INTENT(IN) :: F          ! Fraction of grid box that is
                                         !  precipitating [unitless]
      REAL(fp), INTENT(IN) :: PP         ! Precip rate thru bottom of grid 
                                         !  box (I,J,L)  [cm3 H2O/cm2 air/s]
      REAL(fp), INTENT(IN) :: TK         ! Temperature in grid box [K]
!
! !RETURN VALUE:
!
      REAL(fp)             :: WASHFRAC   ! Fraction of soluble species 
                                         !  lost to washout
! 
! !REVISION HISTORY: 
!  08 Nov 2002 - R. Yantosca - Initial version
!  (1 ) WASHFRAC_AEROSOL used to be an internal function to subroutine WASHOUT.
!        This caused NaN's in the parallel loop on Altix, so we moved it to
!        the module and now pass Iall arguments explicitly (bmy, 7/20/04)
!  16 Sep 2010 - R. Yantosca - Added ProTeX headers
!  16 Aug 2011 - H Amos      - Remove K_WASH from input list, make a defined 
!                              parameter.
!  20 Jan 2012 - H Amos      - WASHFRAC_COARSE_AEROSOL created to handle
!                              SALC and DST4
!  04 Sep 2013 - R. Yantosca - Bug fix: Prevent div-by-zero if F=0.  Because F
!                              multiplies the whole expression for WASHFRAC,
!                              WASHFRAC=0 whenever F=0 anyway.
!EOP
!------------------------------------------------------------------------------
!BOC
!
      !=================================================================
      ! WASHFRAC_COARSE_AEROSOL begins here!
      !=================================================================

      IF ( TK >= 268e+0_fp .OR. ITS_A_POPS_SIM ) THEN

         !---------------------------------
         ! T >= 268K (or POPS simulation)  
         !---------------------------------
         IF ( F > 0e+0_fp ) THEN
            WASHFRAC = F*(1.0_fp - EXP(-0.92_fp * (PP / F*3.6e+4_fp)
     &                ** 0.79_fp * DT / 3.6e+3_fp ))
         ELSE
            WASHFRAC = 0.0_fp
         ENDIF

      ELSE

         !---------------------------------
         ! T < 268K
         !---------------------------------
         IF ( F > 0e+0_fp ) THEN
            WASHFRAC = F *(1.0_fp - EXP( -1.57_fp *
     &                (PP / F*3.6e+4_fp)**0.96_fp * DT / 3.6e+3_fp ))
         ELSE
            WASHFRAC = 0.0_fp
         ENDIF

      ENDIF   

      END FUNCTION WASHFRAC_COARSE_AEROSOL
!EOC
#if defined( TOMAS )
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: washfrac_size_aerosol
!
! !DESCRIPTION: Subroutine WASHFRAC\_SIZE\_AEROSOL retrieves fraction of 
!  soluble aerosol species lost to washout. Size resolved version for TOMAS.
!\\
!\\
! !INTERFACE:
!

      SUBROUTINE WASHFRAC_SIZE_AEROSOL( DT, F, PP, TK, N, I, J, L,
     &                                  State_Met, State_Chm, WASHFRAC,
     &                                  RC ) 
!
! !USES:
!
      USE ErrCode_Mod
      USE ERROR_MOD
      USE State_Chm_Mod,      ONLY : ChmState
      USE State_Met_Mod,      ONLY : MetState
      USE TOMAS_MOD,          ONLY : IBINS, GETDP, STRATSCAV
!
! !INPUT PARAMETERS:
! 
      REAL(fp),       INTENT(IN)    :: DT         ! Dynamic timestep [s]
      REAL(fp),       INTENT(IN)    :: F          ! Fraction of grid box
                                                  !  that is precipitating
      REAL(fp),       INTENT(IN)    :: PP         ! Precip rate thru bottom
                                                  !  of grid box (I,J,L)
                                                  !  [cm3 H2O/cm2 air/s]
      REAL(fp),       INTENT(IN)    :: TK         ! Temperature [K]
      INTEGER,        INTENT(IN)    :: I          ! Longitude index
      INTEGER,        INTENT(IN)    :: J          ! Latitude index
      INTEGER,        INTENT(IN)    :: L          ! Level index
      INTEGER,        INTENT(IN)    :: N          ! Species index

      TYPE(MetState), INTENT(IN)    :: State_Met  ! Meteorology State object
!
! !INPUT/OUTPUT PARAMETERS:
!
      TYPE(ChmState), INTENT(INOUT) :: State_Chm  ! Chemistry State object
!
! !OUTPUT PARAMETERS:
!
      REAL(fp),       INTENT(OUT)   :: WASHFRAC   ! Fraction of species 
                                                  !  lost to washout [1]
      INTEGER,        INTENT(OUT)   :: RC         ! Success or failure?
! 
! !REVISION HISTORY: 
!  31 May 2013 - R. Yantosca - Now accept State_Met, State_Chm as arguments
!  04 Sep 2013 - R. Yantosca - Bug fix: Prevent div-by-zero if F=0.  Because F 
!                              multiplies the whole expression for WASHFRAC,
!                              WASHFRAC=0 whenever F=0 anyway.
!  24 Sep 2015 - E. Lundgren - Function is now a subroutine and outputs RC      
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      REAL(fp)       :: DPAERO          ! Average diameter of particle 
      REAL(fp)       :: SCAVR !Below-cloud scavenging coefficient (per cm rain)
      REAL(fp), SAVE :: SCAVRSAVE(IBINS)
      INTEGER        :: BIN 
      
      !=================================================================
      ! WASHFRAC_SIZE_AEROSOL begins here!
      !=================================================================

      IF ( TK >= 268e+0_fp ) THEN

         !-------------
         ! T >= 268K 
         !-------------

!--------------------------------------------------------------
!!sfarina - This contruct assumes species are dealt with sequentially,
!!          but wetdep parallelizes over species
!!          It could be possible to calculated the lookup table and save
!!          in an I,J,L,BIN array but for now we will calculate redundantly.
!         ! For aerosol number, get Dp and calculate scavr 
!         IF ( N < id_NK1 + IBINS ) THEN
!            DPAERO = GETDP( I, J, L, N, State_Met, State_Chm )
!            ! External function stratscav returns the scavenging rate (mm^-1)
!            ! Let scavr has a unit of cm^-1
!            SCAVR = 10.e+0_fp* STRATSCAV( DPAERO )
!            SCAVRSAVE(N-id_NK1+1) = scavr
!         ELSE
!            BIN = MOD( N - id_NK1 + 1, IBINS )
!            IF( BIN == 0 ) BIN = IBINS
!            SCAVR = SCAVRSAVE(BIN)
!         ENDIF
!---------------------------------------------------------------

            DPAERO = GETDP( I, J, L, N, State_Met, State_Chm, RC )
            ! External function stratscav returns the scavenging rate (mm^-1)
            ! Let scavr has a unit of cm^-1
            SCAVR = 10.e+0_fp* STRATSCAV( DPAERO )

         ! Prevent div by zero (bmy, 9/4/13)
         IF ( F > 0e+0_fp ) THEN
            WASHFRAC = F * ( 1e+0_fp - EXP( -SCAVR * ( PP / F ) * DT ) )
         ELSE
            WASHFRAC = 0e+0_fp
         ENDIF

      ELSE

         !-------------
         ! T < 268K 
         !-------------
         WASHFRAC = 0e+0_fp

      ENDIF

      ! Return to calling program
      END SUBROUTINE WASHFRAC_SIZE_AEROSOL
!EOC
#endif
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: washfrac_hno3
!
! !DESCRIPTION: Function WASHFRAC\_HNO3 returns the fraction of HNO3 
!               species lost to washout.
!\\
!\\
! !INTERFACE:
!
      FUNCTION WASHFRAC_HNO3( DT, F, PP, TK ) RESULT( WASHFRAC )
!
! !USES:
!
!
! !INPUT PARAMETERS: 
!
      REAL(fp), INTENT(IN) :: DT        ! Timestep of washout event [s]
      REAL(fp), INTENT(IN) :: F         ! Fraction of grid box that is
                                        !  precipitating [unitless]
      REAL(fp), INTENT(IN) :: PP        ! Precip rate thru bottom of grid 
                                        !  box (I,J,L)  [cm3 H2O/cm2 air/s]
      REAL(fp), INTENT(IN) :: TK        ! Temperature in grid box [K]
!
! !RETURN VALUE:
!
      REAL(fp)             :: WASHFRAC  ! Fraction of soluble species 

! 
! !REVISION HISTORY: 
!  13 Aug 2011, H Amos: Initial version, modeled after WASHFRAC_AEROSOL.
!                       Seperate function created to emphasize that the new,
!                       updated washout coefficients from Feng et al (2007;
!                       2009) should only be applied to aerosol species. It
!                       was a coincidence before that the original washout
!                       coefficients for aerosols and HNO3 were the same.
!  16 Aug 2011, H Amos: Remove K_WASH from input list, now a defined parameter
!  04 Sep 2013 - R. Yantosca - Bug fix: Prevent div-by-zero if F=0.  Because F 
!                              multiplies the whole expression for WASHFRAC,
!                              WASHFRAC=0 whenever F=0 anyway.
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETER:
!
      REAL(fp), PARAMETER :: K_WASH = 1.0_fp  ! First order washout rate
                                              ! constant [cm^-1].

      !=================================================================
      ! WASHFRAC_HNO3 begins here!
      !=================================================================
      IF ( TK >= 268e+0_fp  ) THEN

         !------------------------
         ! T >= 268K: Do washout
         !------------------------
         IF ( F > 0e+0_fp ) THEN
            WASHFRAC = F * (1e+0_fp - EXP(-K_WASH * (PP / F) * DT))
         ELSE
            WASHFRAC = 0e+0_fp
         ENDIF

      ELSE

         !------------------------
         ! T < 268K: No washout
         !------------------------
         WASHFRAC = 0e+0_fp

      ENDIF

      END FUNCTION WASHFRAC_HNO3
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: washfrac_liq_gas
!!
! !DESCRIPTION: Subroutine WASHFRAC\_LIQ\_GAS returns the fraction of soluble 
!  liquid/gas phase species lost to washout.
!\\
!\\
! !INTERFACE:
! 
      SUBROUTINE WASHFRAC_LIQ_GAS( K0, CR, pKa, PP,       DT,
     &                             F,  DZ, TK,  WASHFRAC, KIN ) 
!     
! !INPUT PARAMETERS: 
!
      REAL(f8), INTENT(IN)  :: K0        ! Henry's solubility constant [M/atm]
      REAL(f8), INTENT(IN)  :: CR        ! Henry's volatility constant [K]
      REAL(f8), INTENT(IN)  :: pKa       ! Henry's pH correction [1]
      REAL(fp), INTENT(IN)  :: PP        ! Precip rate thru bottom of the
                                         !  grid box [cm3 H2O/cm2 air/s]
      REAL(fp), INTENT(IN)  :: DT        ! Timestep for washout event [s]
      REAL(fp), INTENT(IN)  :: F         ! Fraction of grid box that is
                                         !  precipitating [unitless]
      REAL(fp), INTENT(IN)  :: DZ        ! Height of grid box [cm]
      REAL(fp), INTENT(IN)  :: TK        ! Temperature in grid box [K]
!
! !OUTPUT PARAMETERS:
!
      REAL(fp), INTENT(OUT) :: WASHFRAC  ! Fraction of species lost to washout
      LOGICAL,  INTENT(OUT) :: KIN       ! T = washout is a kinetic process
                                         ! F = washout is an equilibrium process
!
! !REVISION HISTORY: 
!  20 Jul 2004 - R. Yantosca - Initial version
!  (1 ) WASHFRAC_LIQ_GAS used to be an internal function to subroutine WASHOUT.
!        This caused NaN's in the parallel loop on Altix, so we moved it to
!        the module and now pass all arguments explicitly (bmy, 7/20/04)
!  16 Sep 2010 - R. Yantosca - Added ProTeX headers 
!  10 Jan 2011 - H. Amos     - Remove AER from the argument list
!  03 Jun 2011 - H. Amos     - convert from a function to a subroutine and
!                              add AER to the argument list
!  16 Aug 2011 - H. Amos     - remove K_WASH from input list, now a defined 
!                              parameter
!  16 Aug 2911 - H. Amos     - rename AER logical KIN to emphasize that washout
!                              is either a kinetic or equilibrium process
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      REAL(fp)            :: L2G, LP, WASHFRAC_F_14
!
! !DEFINED PARAMETERS
!
      REAL(fp), PARAMETER :: K_WASH = 1e+0_fp  ! First order washout rate 
                                               ! constant [cm^-1].

      !=================================================================
      ! WASHFRAC_LIQ_GAS begins here!
      !=================================================================

      ! Start with the assumption that washout will be an 
      ! equilibrium process (H Amos, 03 Jun 2011)
      KIN = .FALSE.

      ! Suppress washout below 268 K
      IF ( TK >= 268e+0_fp ) THEN

         !------------------------
         ! T >= 268K: Do washout
         !------------------------

         ! Rainwater content in the grid box (Eq. 17, Jacob et al, 2000)
         LP = ( PP * DT ) / ( F * DZ ) 

         ! Compute liquid to gas ratio for H2O2, using the appropriate 
         ! parameters for Henry's law -- also use rainwater content Lp
         ! (Eqs. 7, 8, and Table 1, Jacob et al, 2000)
         CALL COMPUTE_L2G( K0, CR, pKa, TK, LP, L2G )

         ! Washout fraction from Henry's law (Eq. 16, Jacob et al, 2000)
         WASHFRAC = L2G / ( 1e+0_fp + L2G )

         ! Washout fraction / F from Eq. 14, Jacob et al, 2000
         ! Note: WASHFRAC_F_14 should match what's used for HNO3 (hma, 13aug2011)
         WASHFRAC_F_14 = 1e+0_fp - EXP( -K_WASH * ( PP / F ) * DT )

         ! Do not let the Henry's law washout fraction exceed
         ! that of HNO3 -- this is a cap
         IF ( WASHFRAC > WASHFRAC_F_14 ) THEN
            WASHFRAC = F * WASHFRAC_F_14
            KIN = .TRUE. ! washout is a kinetic process
         ENDIF
            
      ELSE

         !------------------------
         ! T < 268K: No washout
         !------------------------
         WASHFRAC = 0e+0_fp
            
      ENDIF

      END SUBROUTINE WASHFRAC_LIQ_GAS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: wetdep
!
! !DESCRIPTION: Subroutine WETDEP computes the downward mass flux of 
!  species due to washout and rainout of aerosols and soluble species in a 
!  column.  This subroutine implements an algorithm in which the 
!  precipitation fields come directly from the met archive.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE WETDEP( am_I_Root, Input_Opt,  State_Met,
     &                   State_Chm, State_Diag, LS,        RC )
!
! !USES:
!
      USE CMN_SIZE_MOD
#if defined( BPCH_DIAG )
      USE CMN_DIAG_MOD
#endif
      USE DEPO_MERCURY_MOD, ONLY : ADD_Hg2_WD
      USE DEPO_MERCURY_MOD, ONLY : ADD_HgP_WD
      USE DEPO_MERCURY_MOD, ONLY : ADD_Hg2_SNOWPACK
      USE ErrCode_Mod
      USE ERROR_MOD
      USE Input_Opt_Mod,    ONLY : OptInput
      USE State_Chm_Mod,    ONLY : ChmState
      USE State_Diag_Mod,   ONLY : DgnState
      USE State_Met_Mod,    ONLY : MetState
      USE TIME_MOD,         ONLY : GET_TS_DYN
      USE Species_Mod,      ONLY : Species
      USE UnitConv_Mod,     ONLY : Convert_Spc_Units
!
! !INPUT PARAMETERS: 
!
      LOGICAL,        INTENT(IN)    :: am_I_Root   ! Are we on the root CPU?
      LOGICAL,        INTENT(IN)    :: LS          ! =T for large-scale precip
                                                   ! =F for convective precip
      TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input options
      TYPE(MetState), INTENT(IN)    :: State_Met   ! Meteorology State object
!
! !INPUT/OUTPUT PARAMETERS: 
!
      TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
      TYPE(DgnState), INTENT(INOUT) :: State_Diag  ! Diagnostics State object
!
! !OUTPUT PARAMETERS:
!
      INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
!
! !REMARKS:
!  Precipitation fields:
!  =====================
!                                                                             .
!       Layer        Formation of       Precipitation
!                     New Precip        falling down
!      ==================================================== Top of Atm.
!        LM           QQ(L,I,J)         PDOWN(LM,I,J)
!                         |                   |
!      ----------------------------------------------------
!        LM-1         QQ(L,I,J)         PDOWN(LM-1,I,J)    
!                         |                   |
!      -------------------V-------------------V------------
!                        ...                 ...     
!                         |                   |
!      -------------------V-------------------V------------
!        3            QQ(L,I,J)         PDOWN(3,I,J)
!                         |                   |
!      -------------------V--------------------------------
!        2            QQ(L,I,J)         PDOWN(2,I,J) 
!                         |                   |
!      ----------------------------------------------------
!        1            QQ(L,I,J)         PDOWN(1,I,J) 
!                         |                   |
!      ===================V===================V============ Ground
!                                                                             .
!  Where:
!    (a) New formation forming in grid box (I,J,L) = QQ(L,I,J)
!    (b) Precip coming in  thru top    of layer L  = PDOWN(L+1,I,J)
!    (c) Precip going  out thru bottom of layer L  = PDOWN(L,  I,J) 
!                                                                             .
!  Rainout:
!  ========   
!  Rainout occurs when there is more precipitation in grid box (I,J,L) than
!  in grid box (I,J,L+1).  In other words, rainout occurs when the amount of 
!  rain falling through the bottom of grid box (I,J,L) is more than the amount
!  of rain coming in through the top of grid box (I,J,L). 
!                                                                             .
!  Soluble gases/aerosols are incorporated into the raindrops and are 
!  completely removed from grid box (I,J,LLPAR).  There is no evaporation 
!  and "resuspension" of aerosols during a rainout event.
!                                                                             .
!  For large-scale (a.k.a. stratiform) precipitation, the first order rate 
!  constant for rainout in the grid box (I,J,L=LLPAR) (cf. Eq. 12, Jacob 
!  et al, 2000) is given by:
!                                                                             .
!                           Q        
!       K_RAIN = K_MIN + -------    [units: s^-1]
!                         L + W    
!                                                                             .
!  and the areal fraction of grid box (I,J,L=LLPAR) that is actually 
!  experiencing large-scale precipitation (cf. Eq. 11, Jacob et al, 2000) 
!  is given by: 
!                                                                             .
!                         Q               
!       F'     =  -------------------   [unitless]
!                  K_RAIN * ( L + W )    
!                                                                             .
!  Where:
!                                                                             .
!       K_MIN  = minimum value for K_RAIN         
!              = 1.0e-4 [s^-1]
!                                                                             .
!       L + W  = condensed water content in cloud 
!              = 1.5e-6 [cm3 H2O/cm3 air]
!                                                                             .
!       Q = QQ = rate of precipitation formation 
!                [ cm3 H2O / cm3 air / s ]
!                                                                             .
!  For convective precipitation, K_RAIN = 5.0e-3 [s^-1], and the expression 
!  for F' (cf. Eq. 13, Jacob et al, 2000) becomes:
!                                                                             .
!                                       { DT        }
!                         FMAX * Q * MIN{ --- , 1.0 }
!                                       { TAU       }
!       F' = ------------------------------------------------------
!                    { DT        }
!             Q * MIN{ --- , 1.0 }  +  FMAX * K_RAIN * ( L + W )
!                    { TAU       } 
!                                                                             .
!  Where:
!                                                                             .
!       Q = QQ = rate of precipitation formation 
!              [cm3 H2O/cm3 air/s]
!                                                                             .
!       FMAX   = maximum value for F' 
!              = 0.3
!                                                                             .
!       DT     = dynamic time step from the CTM [s]
!                                                                             .
!       TAU    = duration of rainout event 
!              = 1800 s (30 min)
!                                                                             .
!       L + W  = condensed water content in cloud 
!              = 2.0e-6 [cm3 H2O/cm3 air]
!                                                                             .
!  K_RAIN and F' are needed to compute the fraction of species in grid box 
!  (I,J,L=LLPAR) lost to rainout.  This is done in module routine RAINOUT.
!                                                                             .
!  Washout:
!  ========   
!  Washout occurs when we have evaporation (or no precipitation at all) at 
!  grid box (I,J,L), but have rain coming down from grid box (I,J,L+1).
!                                                                             .
! !REVISION HISTORY: 
!  20 Sep 2010 - R. Yantosca - Initial version, based on WETDEP
!  28 Sep 2010 - H. Amos     - Now define Q, QDOWN directly from MERRA met
!  08 Oct 2010 - R. Yantosca - Adjusted OpenMP do loop
!  09-Dec-2010 - H. Amos     - Added PDOWN(L+1) > 0 to criterion for IS_WASHOUT
!  09-Dec-2010 - H. Amos     - SAFETY now prints PDOWN(L+1) instead of PDOWN(L)
!  31-Dec-2010 - H. Amos     - Clean up code, remove obsolete code
!  31-Dec-2010 - H. Amos     - Added comments
!  27 May 2011 - R. Yantosca - Now pass F_RAINOUT to DO_WASHOUT_ONLY
!  25 Aug 2014 - M. Sulprizio- Now accept Input_Opt as an argument
!  01 Apr 2015 - L. Zhang    - Don't do wetdep in nested-grid buffer zone
!  02 Apr 2015 - E. Lundgren - Move tracer unit conversion from kg/kg to kg
!                              to within this routine
!  04 Jun 2015 - E. Lundgren - Now accept am_I_Root and RC as arguments
!  09 Jun 2015 - R. Yantosca - Now deposit Hg2, HgP to snowpack regardless of
!                              whether the dynamic ocean is used
!  04 Mar 2016 - C. Keller   - Added option to scale convective fraction of
!                              large-scale precip.
!  22 Apr 2016 - R. Yantosca - Now get Is_Hg2, Is_HgP from species database
!  25 Apr 2016 - R. Yantosca - Now get the Hg category # from species database
!  29 Apr 2016 - R. Yantosca - Don't initialize pointers in declaration stmts
!  30 Jun 2016 - R. Yantosca - Remove instances of STT.  Now get the advected
!                              species ID from State_Chm%Map_Advect.
!  30 Jun 2016 - R. Yantosca - Replace STT with Spc and DSTT with DSpc
!  01 Jul 2016 - R. Yantosca - Now rename species DB object ThisSpc to SpcInfo
!  05 Jul 2016 - R. Yantosca - Now replace IDWETD with State_Chm%Map_Wetdep
!  05 Aug 2016 - R. Yantosca - Remove temporary tracer-removal code
!  05 Aug 2016 - R. Yantosca - Remove reference to Spc, it's not needed
!  17 Mar 2017 - R. Yantosca - Replace NSOL with State_Chm%nWetDep
!  24 Aug 2017 - M. Sulprizio- Rename routine from WETDEP_MERRA to WETDEP
!  28 Sep 2017 - E. Lundgren - Simplify unit conversions using wrapper routine
!  09 Nov 2017 - R. Yantosca - Now accept State_Diag as an argument
!  09 Nov 2017 - R. Yantosca - Return error condition to calling program
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      ! SAVEd Scalars
      LOGICAL, SAVE          :: FIRST = .TRUE.

      ! Scalars
      LOGICAL                :: IS_Hg
      LOGICAL                :: KIN
      LOGICAL                :: IS_RAINOUT, IS_WASHOUT, IS_BOTH
      INTEGER                :: I,     IDX,    J,         L
      INTEGER                :: N,     NW,     Hg_Cat,    EC
      REAL(fp)               :: Q,     QDOWN,  DT,        DT_OVER_TAU
      REAL(fp)               :: K,     K_MIN,  K_RAIN,    RAINFRAC
      REAL(fp)               :: F,     FTOP,   F_PRIME,   WASHFRAC
      REAL(fp)               :: LOST,  GAINED, MASS_WASH, MASS_NOWASH
      REAL(fp)               :: ALPHA, ALPHA2, WETLOSS,   TMP
      REAL(fp)               :: F_RAINOUT,     F_WASHOUT 
      REAL(fp)               :: DEP_HG
      REAL(fp)               :: CNVSCL
      CHARACTER(LEN=63)      :: OrigUnit

      ! Arrays
      ! DSpc is the accumulator array of rained-out 
      ! soluble species for a given (I,J) column
      REAL(fp)               :: DSpc(State_Chm%nWetDep,
     &                               LLPAR,IIPAR,JJPAR)
 
      ! Strings
      CHARACTER(LEN=255)     :: ErrMsg, ErrorMsg, ThisLoc

      ! Pointers
      TYPE(Species), POINTER :: SpcInfo

      !=================================================================
      ! (1)  I n i t i a l i z e   V a r i a b l e s
      !=================================================================

      ! Initialize
      RC        = GC_SUCCESS
      ErrorMsg  = ''
      ThisLoc   = ' -> at WetDep (in module GeosCore/wetscav_mod.F)'

      ! Is this a mercury simulation?
      IS_Hg = ITS_A_MERCURY_SIM

      ! Initialize pointers
      SpcInfo => NULL()

      ! Convert species concentration to mass per unit area (kg/m2) for 
      ! wet deposition since computation is done per column (ewl, 9/8/15)
      CALL Convert_Spc_Units( am_I_Root, Input_Opt, State_Met, 
     &                        State_Chm, 'kg/m2', RC, 
     &                        OrigUnit=OrigUnit )

      ! Trap potential errors
      IF ( RC /= GC_SUCCESS ) THEN
         ErrorMsg = 'Unit conversion error at start of WETDEP!'
         CALL GC_Error( ErrorMsg, RC, ThisLoc )
         RETURN
      ENDIF  

      ! Dynamic timestep [s]
      DT  = GET_TS_DYN()
      
      ! Select index for diagnostic arrays -- will archive either
      ! large-scale or convective rainout/washout fractions
      !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      !%%%% NOTE: FOR GEOS-FP AND MERRA-2, WE ALWAYS ASSUME THAT   %%%%%
      !%%%% WETDEP OPERATES ONLY ON LARGE-SCALE PRECIP.  THUS, LS  %%%%%
      !%%%% IS ALWAYS SET TO "TRUE" IN THE CALLING ROUTINE, AND    %%%%%
      !%%%% IDX WILL ALWAYS BE 1 HERE AND IN SUBSEQUENT ROUTINES.  %%%%%  
      !%%%% (bmy, 11/13/17)                                        %%%%%
      !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      IF ( LS ) THEN
         IDX = 1
      ELSE
         IDX = 2
      ENDIF

      ! CNVSCL is the scale factor to be applied to convective large
      ! scale precipitation: a fraction of model-resolved vertical 
      ! transport is of convective nature (especially at high model
      ! resolution). Treating all of it as large-scale precipiation 
      ! may overestimate total washout. CNVSCL is a scale factor to
      ! adjust the convective fraction of the precipitation. If set
      ! to 1.0, it will yield the same result as without correction.
      ! If set to 0.0, all convective large-scale precipitation is 
      ! suppressed (ckeller, 3/4/16).
      CNVSCL = -1.0 
      IF ( ASSOCIATED(State_Met%CNV_FRC) ) THEN
         CNVSCL = Input_Opt%WETD_CONV_SCAL
         CNVSCL = MIN(MAX(CNVSCL,0.0_fp),1.0_fp)
      ENDIF

      !=================================================================
      ! (2)  L o o p   O v e r   (I, J)   S u r f a c e   B o x e s
      !=================================================================
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I,           J,          FTOP,        L          )    
!$OMP+PRIVATE( NW,          ERRMSG,     F,           F_PRIME    )
!$OMP+PRIVATE( F_RAINOUT,   F_WASHOUT,  K_RAIN,      Q          )       
!$OMP+PRIVATE( QDOWN,       IS_RAINOUT, IS_WASHOUT,  N          ) 
!$OMP+PRIVATE( DEP_HG,      SpcInfo,    Hg_Cat,      EC         )
!$OMP+SCHEDULE( DYNAMIC )   
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! Initialize PRIVATE error-handling variables
         EC        = GC_SUCCESS
         ErrorMsg  = ''

         ! Don't do wetdep in nested-grid buffer zone (lzh, 4/1/15)
         IF ( Input_Opt%ITS_A_NESTED_GRID ) THEN
           IF ( J <=          Input_Opt%NESTED_J0W ) CYCLE
           IF ( J >   JJPAR - Input_Opt%NESTED_J0E ) CYCLE
           IF ( I <=          Input_Opt%NESTED_I0W ) CYCLE
           IF ( I >   IIPAR - Input_Opt%NESTED_I0E ) CYCLE
        ENDIF

         ! Zero FTOP
         FTOP = 0e+0_fp

         ! Zero accumulator array
         DO L  = 1, LLPAR
         DO NW = 1, State_Chm%nWetDep
            DSpc(NW,L,I,J) = 0e+0_fp
         ENDDO
         ENDDO

         !==============================================================
         ! (3)  R a i n o u t   F r o m   T o p   L a y e r  (L = LLPAR) 
         !==============================================================

         ! Zero variables for this level
         ERRMSG      = 'RAINOUT: Top of atm'
         F           = 0e+0_fp
         F_PRIME     = 0e+0_fp
         F_RAINOUT   = 0e+0_fp      
         F_WASHOUT   = 0e+0_fp      
         K_RAIN      = 0e+0_fp
         Q           = 0e+0_fp     

         ! Start at the top of the atmosphere
         L = LLPAR

         ! If precip forms at (I,J,L), assume it all rains out
         IF ( QQ(L,I,J) > 0e+0_fp ) THEN

            ! Q is the new precip that is forming within grid box (I,J,L)
            Q = QQ(L,I,J)

            ! Compute K_RAIN and F_RAINOUT for large-scale 
            ! precipitation (cf. Eqs. 11-13, Jacob et al, 2000) 
            K_RAIN  = LS_K_RAIN( Q )
            F_RAINOUT = LS_F_PRIME( Q, K_RAIN )
            
            ! Set F = F_RAINOUT, since there is no FTOP at L = LLPAR
            F = F_RAINOUT

            ! Adjust convective large-scale precip (ckeller, 3/4/16)
            IF ( CNVSCL >= 0.0_fp ) THEN
               F = ( ( 1.0_fp - State_Met%CNV_FRC(I,J) ) * F )
     &           + ( ( CNVSCL * State_Met%CNV_FRC(I,J) ) * F )
            ENDIF       

            ! Only compute rainout if F > 0. 
            ! This helps to eliminate unnecessary CPU cycles.
            IF ( F > 0e+0_fp ) THEN 
               CALL DO_RAINOUT_ONLY( am_I_Root  = am_I_Root,
     &                               LS         = LS,
     &                               I          = I, 
     &                               J          = J,
     &                               L          = L,
     &                               IDX        = IDX,        
     &                               ERRMSG     = ERRMSG,
     &                               F_RAINOUT  = F,
     &                               K_RAIN     = K_RAIN,     
     &                               DT         = DT,
     &                               DSpc       = DSpc,
     &                               Input_Opt  = Input_Opt,
     &                               State_Met  = State_Met,
     &                               State_Chm  = State_Chm,
     &                               State_Diag = State_Diag, 
     &                               RC         = EC          )

               ! Trap potential errors
               IF ( EC /= GC_SUCCESS ) THEN
                  RC       = EC
                  ErrorMsg = 
     &               'Error encountered in "Do_Rainout_Only (3)"!'
               ENDIF
            ENDIF

            ! Save FTOP for the next lower level 
            FTOP = F

         ENDIF         

         !==============================================================
         ! (4)  R a i n o u t   a n d   W a s h o u t 
         !      i n   t h e   M i d d l e   L e v e l s
         !==============================================================
         DO L = LLPAR-1, 2, -1

            ! Zero variables for each level
            F           = 0e+0_fp
            F_PRIME     = 0e+0_fp
            F_RAINOUT   = 0e+0_fp 
            F_WASHOUT   = 0e+0_fp 
            K_RAIN      = 0e+0_fp
            Q           = 0e+0_fp
            QDOWN       = 0e+0_fp

            ! If there is new precip forming w/in the grid box ...
            IF ( QQ(L,I,J) > 0e+0_fp ) THEN

               ! Compute K_RAIN and F_RAINOUT for large-scale 
               ! precipitation (cf. Eqs. 11-13, Jacob et al, 2000) 
               ! F_RAINOUT is the fraction of grid box (I,J,L) 
               ! experiencing rainout.
               K_RAIN  = LS_K_RAIN( QQ(L,I,J) )
               F_PRIME = LS_F_PRIME( QQ(L,I,J), K_RAIN )
               
            ELSE
               
               F_PRIME = 0e+0_fp
 
            ENDIF

            ! The following block implements Qiaoqiao's changes
            ! Calculate the fractional areas subjected to rainout and
            ! washout. If PDOWN = 0, then all dissolved species returns
            ! to the atmosphere. (cdh, 7/13/10)
            IF ( PDOWN(L,I,J) > 0e+0_fp ) THEN
               F_RAINOUT = F_PRIME
               ! Washout occurs where there is no rainout
               F_WASHOUT = MAX( FTOP - F_RAINOUT, 0e+0_fp )
            ELSE
               F_RAINOUT = 0e+0_fp
               F_WASHOUT = 0e+0_fp
            ENDIF

            ! Adjust convective large-scale precip (ckeller, 3/4/16)
            IF ( CNVSCL >= 0.0_fp ) THEN
               F_RAINOUT = 
     &             ( ( 1.0_fp - State_Met%CNV_FRC(I,J) ) * F_RAINOUT )
     &           + ( ( CNVSCL * State_Met%CNV_FRC(I,J) ) * F_RAINOUT )
               F_WASHOUT = 
     &             ( ( 1.0_fp - State_Met%CNV_FRC(I,J) ) * F_WASHOUT )
     &           + ( ( CNVSCL * State_Met%CNV_FRC(I,J) ) * F_WASHOUT )
            ENDIF       


            IF ( F_WASHOUT > 0e+0_fp ) THEN    

               ! QDOWN is the precip leaving thru the bottom of box (I,J,L)
               ! Q     is the precip that is evaporating within box (I,J,L) 
               QDOWN = PDOWN(L,I,J)               
               Q     = REEVAP(L,I,J)                         
	       Q     = MAX( Q, 0e+0_fp ) ! Negative values are unphysical
               
               !  Define PDOWN and p
               IF ( F_RAINOUT > 0e+0_fp ) THEN

                  ! The precipitation causing washout 
                  ! is the precip entering thru the top
                  QDOWN = PDOWN(L+1,I,J)            
                   
                  !** GEOS-FP and MERRA-2 distinguish between rates of
                  ! new precipitation (field DQRLSAN) and evaporation of
                  ! precipitation (field REEVAPLS).
                  ! So the assumption below is not required. ** VPS (6/21/16)
                  
                  ! The amount of precipitating water entering from above 
                  ! which evaporates. If there is rainout (new precip
                  ! forming) then we have no way to estimate this, so assume
                  ! zero for now. Consequently there will be no resuspended
                  ! aerosol.
                  !Q = 0e+0_fp
               ENDIF
            ENDIF 

            !-----------------------------------------------------------
            ! Determine if we have the following conditions:
            !
            ! (a) Rainout 
            ! (b) Washout 
            !
            ! Note that rainout and washout can happen in the same 
            ! grid box.
            !-----------------------------------------------------------

            ! If a non-zero fraction of the grid box is 
            ! experiencing rainout...
            IS_RAINOUT = ( F_RAINOUT > 0e+0_fp )

            ! If a non-zero fraction of the grid box is 
            ! experiencing washout...
            IS_WASHOUT = ( F_WASHOUT > 0e+0_fp )

            IF ( IS_RAINOUT ) THEN

               !--------------------------------------------------------
               ! RAINOUT 
               !--------------------------------------------------------

               ! Error msg for stdout
               ERRMSG = 'RAINOUT'

               ! Do rainout if we meet the above criteria
               CALL DO_RAINOUT_ONLY( am_I_Root  = am_I_Root,
     &                               LS         = LS,
     &                               I          = I, 
     &                               J          = J,
     &                               L          = L,
     &                               IDX        = IDX,        
     &                               ERRMSG     = ERRMSG,
     &                               F_RAINOUT  = F_RAINOUT,
     &                               K_RAIN     = K_RAIN,     
     &                               DT         = DT,
     &                               DSpc       = DSpc,
     &                               Input_Opt  = Input_Opt,
     &                               State_Met  = State_Met,
     &                               State_Chm  = State_Chm,
     &                               State_Diag = State_Diag, 
     &                               RC         = EC          )

               ! Trap potential errors
               IF ( EC /= GC_SUCCESS ) THEN
                  RC     = EC
                  ErrorMsg = 
     &               'Error encountered in "Do_Rainout_Only (4)!'
               ENDIF

            ENDIF  

            IF ( IS_WASHOUT ) THEN 

               !--------------------------------------------------------
               ! WASHOUT ONLY
               !--------------------------------------------------------

               ! Error msg for stdout
               ERRMSG = 'WASHOUT'

               ! Do the washout
               CALL DO_WASHOUT_ONLY( am_I_Root  = am_I_Root, 
     &                               LS         = LS,        
     &                               I          = I,
     &                               J          = J,         
     &                               L          = L,         
     &                               IDX        = IDX,
     &                               ERRMSG     = ERRMSG,    
     &                               QDOWN      = QDOWN,     
     &                               Q          = Q,
     &                               F_WASHOUT  = F_WASHOUT, 
     &                               F_RAINOUT  = F_RAINOUT,
     &                               DT         = DT,
     &                               PDOWN      = PDOWN,     
     &                               DSpc       = DSpc,      
     &                               Input_Opt  = Input_Opt,  
     &                               State_Met  = State_Met,
     &                               State_Chm  = State_Chm,
     &                               State_Diag = State_Diag,
     &                               RC         = EC          )

               ! Trap potential errors
               IF ( EC /= GC_SUCCESS ) THEN
                  RC       = EC
                  ErrorMsg = 
     &               'Error encountered in "Do_Washout_Only (4)!'
               ENDIF
            ENDIF

            !===========================================================
            ! (6)  N o   D o w n w a r d   P r e c i p i t a t i o n 
            !
            ! If there is no precipitation leaving grid box (I,J,L), 
            ! then  set F, the effective area of precipitation in grid 
            ! box (I,J,L), to zero.
            !
            ! Also, all of the previously rained-out species that is now 
            ! coming down from grid box (I,J,L+1) will evaporate and 
            ! re-enter the atmosphere in the gas phase in grid box 
            ! (I,J,L).  This is called "resuspension".
            !===========================================================

            ! Check if there is precip entering grid box, but not
            ! leaving grid box
            IF ( F_WASHOUT == 0e+0_fp .and. F_RAINOUT == 0e+0_fp ) THEN

               ! No precipitation at grid box (I,J,L), thus F = 0
               F = 0e+0_fp
              
               ! Error message
               ERRMSG = 'RESUSPENSION in middle levels'

               ! Re-evaporate all of the rain
               CALL DO_COMPLETE_REEVAP( am_I_Root  = am_I_Root,
     &                                  LS         = LS,       
     &                                  I          = I, 
     &                                  J          = J,
     &                                  L          = L,
     &                                  IDX        = IDX,       
     &                                  ERRMSG     = ERRMSG,
     &                                  DT         = DT,
     &                                  DSpc       = DSpc,      
     &                                  Input_Opt  = Input_Opt,
     &                                  State_Met  = State_Met,
     &                                  State_Chm  = State_Chm, 
     &                                  State_Diag = State_Diag,
     &                                  RC         = EC         )

               ! Trap potential errors
               IF ( EC /= GC_SUCCESS ) THEN
                  RC       = EC
                  ErrorMsg = 
     &               'Error encountered in "Do_Complete_Reevap" (6)!'
               ENDIF
            ENDIF 

            ! Save FTOP for next level
            FTOP = F_RAINOUT + F_WASHOUT

         ENDDO               

         !==============================================================
         ! (7)  W a s h o u t   i n   L e v e l   1
         !==============================================================

         ! Zero variables for this level
         ERRMSG  = 'WASHOUT: at surface'
         F       = 0e+0_fp
         F_PRIME = 0e+0_fp
         K_RAIN  = 0e+0_fp
         Q       = 0e+0_fp
         QDOWN   = 0e+0_fp
         
         ! We are at the surface, set L = 1
         L = 1

         ! Washout at level 1 criteria
         IF ( PDOWN(L+1,I,J) > 0e+0_fp ) THEN

            ! QDOWN is the precip leaving thru the bottom of box (I,J,L+1)
            QDOWN = PDOWN(L+1,I,J)

            ! Since no precipitation is forming within grid box (I,J,L),
            ! F' = 0, and F = MAX( F', FTOP ) reduces to F = FTOP.
            F = FTOP

            ! Only compute washout if F > 0.
            IF ( F > 0e+0_fp ) THEN
               CALL DO_WASHOUT_AT_SFC( am_I_Root  = am_I_Root,
     &                                 LS         = LS,
     &                                 I          = I,
     &                                 J          = J,
     &                                 L          = L,
     &                                 IDX        = IDX,
     &                                 ERRMSG     = ERRMSG,
     &                                 QDOWN      = QDOWN,
     &                                 F          = F,
     &                                 DT         = DT,
     &                                 DSpc       = DSpc,
     &                                 Input_Opt  = Input_Opt,
     &                                 State_Met  = State_Met,
     &                                 State_Chm  = State_Chm, 
     &                                 State_Diag = State_Diag, 
     &                                 RC         = EC          )

               ! Trap potential errors
               IF ( EC /= GC_SUCCESS ) THEN
                  RC       = EC
                  ErrorMsg = 
     &               'Error encountered in "Do_Washout_at_Sfc (7)!'
               ENDIF
            ENDIF    
         ENDIF

         !==============================================================
         ! (8)  M e r c u r y   S i m u l a t i o n   O n l y 
         !
         ! For the mercury simulation, we need to archive the amt of 
         ! Hg2 [kg] that is scavenged out of the column.  Also applies
         ! to the tagged Hg simulation.
         !
         ! NOTES:
         ! (a) Now moved outside the loop above for clarity and to 
         !      fix a bug where HgP scavenging was not recorded. 
         ! (b) The values of DSpc in the first layer accumulates all 
         !      scavenging and washout in the column
         ! (c) Updates from cdh. (ccc, 5/17/10)
         !==============================================================
         IF ( IS_Hg ) THEN

            ! Loop over soluble species and/or aerosol species
            DO NW = 1, State_Chm%nWetDep

               ! Get the species index from the wetdep index
               N       = State_Chm%Map_WetDep(NW)

               ! Amount of Hg wet-deposited out of the column
               DEP_HG  = DSpc(NW,1,I,J)

               ! Point to the Species Database entry for species N
               SpcInfo => State_Chm%SpcData(N)%Info

               ! Check if it is a gaseous Hg2 tag
               IF ( SpcInfo%Is_Hg2 ) THEN
                  
                  ! Get the category # for this Hg2 species
                  Hg_Cat = SpcInfo%Hg_Cat

                  ! Archive wet-deposited Hg2
                  CALL ADD_Hg2_WD      ( I, J, Hg_Cat,    DEP_HG     )
                  CALL ADD_Hg2_SNOWPACK( I, J, Hg_Cat,    DEP_HG, 
     &                              State_Met, State_Chm, State_Diag )

               ! Check if it is a HgP tag
               ELSE IF ( SpcInfo%Is_HgP ) THEN
                  
                  ! Get the category # for this HgP species
                  Hg_Cat = SpcInfo%Hg_Cat

                  ! Archive wet-deposited HgP
                  CALL ADD_HgP_WD      ( I, J, Hg_Cat,    DEP_HG     )
                  CALL ADD_Hg2_SNOWPACK( I, J, Hg_Cat,    DEP_HG, 
     &                              State_Met, State_Chm, State_Diag )
                  
               ENDIF

               ! Free pointer
               SpcInfo => NULL()
            ENDDO
            
         ENDIF
         
      ENDDO	
      ENDDO	
!$OMP END PARALLEL DO

      ! Exit with error condition 
      IF ( RC /= GC_SUCCESS ) THEN
         CALL GC_Error( ErrorMsg, RC, ThisLoc )
         RETURN
      ENDIF

      ! Convert species concentration back to original unit (ewl, 9/8/15)
      CALL Convert_Spc_Units( am_I_Root, Input_Opt, State_Met, 
     &                        State_Chm, OrigUnit,  RC )

      ! Trap potential errors
      IF ( RC /= GC_SUCCESS ) THEN
         ErrMsg = 'Unit conversion error at end of WETDEP!'
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF  

      END SUBROUTINE WETDEP
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ls_k_rain
!
! !DESCRIPTION: Function LS\_K\_RAIN computes K\_RAIN, the first order 
!  rainout rate constant for large-scale (a.k.a. stratiform) precipitation.
!\\
!\\
! !INTERFACE:
!
      FUNCTION LS_K_RAIN( Q ) RESULT( K_RAIN )
!
! !INPUT PARAMETERS: 
!
      REAL(fp), INTENT(IN) :: Q        ! Rate of precipitation formation 
                                       !  [cm3 H2O/cm3 air/s]
!
! !RETURN VALUE:
!
      REAL(fp)             :: K_RAIN   ! 1st order rainout rate constant [1/s]
!
! !REVISION HISTORY: 
!  18 Mar 2004 - R. Yantosca - Initial version
!  (1 ) Now made into a MODULE routine since we cannot call internal routines
!        from w/in a parallel loop.  Updated comments. (bmy, 3/18/04)
!  16 Sep 2010 - R. Yantosca - Added ProTeX headers
!  23 Oct 2015 - R. Yantosca - Now use parameter COND_WATER_CONTENT
!EOP
!------------------------------------------------------------------------------
!BOC
      !==================================================================
      ! LS_K_RAIN begins here!
      !==================================================================

      ! Compute rainout rate constant K in s^-1 (Eq. 12, Jacob et al, 2000).
      ! (1) 1.0d-4             = K_MIN, a minimum value for K_RAIN
      ! (2) COND_WATER_CONTENT = L + W [cm3/cm3], the condensed water
      !                          content (liquid + ice)  in the cloud, 
      !                          (cf. Jacob et al 2000, Eq. 12)
      K_RAIN = 1.0e-4_fp + ( Q / COND_WATER_CONTENT )
      
      END FUNCTION LS_K_RAIN
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ls_f_prime
!
! !DESCRIPTION: Function LS\_F\_PRIME computes F', the fraction of the 
!  grid box that is precipitating during large scale (a.k.a. stratiform) 
!  precipitation.
!\\
!\\
! !INTERFACE:
!
      FUNCTION LS_F_PRIME( Q, K_RAIN ) RESULT( F_PRIME )
!
! !INPUT PARAMETERS: 
!
      REAL(fp), INTENT(IN) :: Q        ! Rate of precipitation formation 
                                       !  [cm3 H2O/cm3 air/s]
      REAL(fp), INTENT(IN) :: K_RAIN   ! 1st order rainout rate constant [1/s]
!
! !REMARKS:
! 
      REAL(fp)             :: F_PRIME  ! Fraction of grid box undergoing 
                                       !  large-scale precipitation [unitless]
! 
! !REVISION HISTORY: 
!  18 Mar 2004 - R. Yantosca - Initial version
!  (1 ) Now made into a MODULE routine since we cannot call internal routines
!        from w/in a parallel loop.  Updated comments. (bmy, 3/18/04)
!  16 Sep 2010 - R. Yantosca - Added ProTeX headers
!  23 Oct 2015 - R. Yantosca - Now use COND_WATER_CONTENT parameter
!EOP
!------------------------------------------------------------------------------
!BOC
      !=================================================================
      ! LS_F_PRIME begins here!
      !=================================================================

      ! Compute F', the area of the grid box undergoing precipitation
      ! COND_WATER_CONTENT = L + W [cm3/cm3], the condensed water
      !                      content (liquid + ice)  in the cloud 
      !                      (cf. Jacob et al 2000, Eq. 12)
      F_PRIME = Q / ( K_RAIN * COND_WATER_CONTENT )

      END FUNCTION LS_F_PRIME
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: conv_f_prime
!
! !DESCRIPTION: Function CONV\_F\_PRIME computes F', the fraction of the 
!  grid box that is precipitating during convective precipitation.
!\\
!\\
! !INTERFACE:
!
      FUNCTION CONV_F_PRIME( Q, K_RAIN, DT ) RESULT( F_PRIME )
!
! !INPUT PARAMETERS: 
!
      REAL(fp), INTENT(IN) :: Q         ! Rate of precipitation formation 
                                      !  [cm3 H2O/cm3 air/s]
      REAL(fp), INTENT(IN) :: K_RAIN    ! 1st order rainout rate constant [1/s]
      REAL(fp), INTENT(IN) :: DT        ! Wet deposition timestep [s]
!
! !RETURN VALUE:
!
      REAL(fp)             :: F_PRIME   ! Frac. of grid box undergoing
                                      !  convective precipitation [unitless]
! 
! !REVISION HISTORY: 
!  18 Mar 2004 - R. Yantosca - Initial version
!  (1 ) Now made into a MODULE routine since we cannot call internal routines
!        from w/in a parallel loop.  Updated comments. (bmy, 3/18/04)
!  16 Sep 2010 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      REAL(fp) :: TIME

      !=================================================================
      ! CONV_F_PRIME begins here!
      !=================================================================
      
      ! Assume the rainout event happens in 30 minutes (1800 s)
      ! Compute the minimum of DT / 1800s and 1.0
      TIME = MIN( DT / 1800e+0_fp, 1e+0_fp )

      ! Compute F' for convective precipitation (Eq. 13, Jacob et al, 2000)
      ! 0.3  = FMAX, the maximum value of F' for convective precip
      ! 2d-6 = L + W, the condensed water content [cm3 H2O/cm3 air]
      F_PRIME = ( 0.3e+0_fp * Q * TIME ) / 
     &          ( ( Q * TIME ) + ( 0.3e+0_fp * K_RAIN * 2e-6_fp ) )

      END FUNCTION CONV_F_PRIME
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: do_rainout_only
!
! !DESCRIPTION: Subroutine DO\_RAINOUT\_ONLY removes species by rainout.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE DO_RAINOUT_ONLY( am_I_Root,  LS,        I,
     &                            J,          L,         IDX,
     &                            ERRMSG,     F_RAINOUT, K_RAIN,
     &                            DT,         DSpc,      Input_Opt,  
     &                            State_Met,  State_Chm, State_Diag,
     &                            RC                                 )
!
! !USES:
!
      USE CMN_SIZE_MOD                                 ! Size parameters
#if defined( BPCH_DIAG )
      USE CMN_DIAG_MOD                                 ! Diagnostic flags
      USE DIAG_MOD,       ONLY : AD16                  ! ND16 diag array
      USE DIAG_MOD,       ONLY : AD17                  ! ND17 diag array
      USE DIAG_MOD,       ONLY : AD39                  ! ND39 diag array
      USE DIAG_MOD,       ONLY : CT16                  ! ND16 diag counter
      USE DIAG_MOD,       ONLY : CT17                  ! ND17 diag counter
#endif
      USE ErrCode_Mod
      USE ERROR_MOD,      ONLY : IT_IS_NAN             ! Test for NaN
      USE GET_NDEP_MOD,   ONLY : SOIL_WETDEP           ! Wet deposited species
      USE Input_Opt_Mod,  ONLY : OptInput              ! Input options type
      USE State_Chm_Mod,  ONLY : ChmState              ! Chemistry State type
      USE State_Met_Mod,  ONLY : MetState              ! Met State type
      USE State_Diag_Mod, ONLY : DgnState
#if defined( TOMAS )
      USE TOMAS_MOD,      ONLY : IBINS, ICOMP, AQOXID
      USE TOMAS_MOD,      ONLY : GETFRACTION
#if defined( BPCH_DIAG )
      USE DIAG_MOD,       ONLY : AD05
#endif
#endif
!
! !INPUT PARAMETERS: 
!     
      LOGICAL,          INTENT(IN)    :: am_I_Root     ! Are we on root CPU?    
      LOGICAL,          INTENT(IN)    :: LS            ! =T denotes LS precip
      INTEGER,          INTENT(IN)    :: I             ! Longitude index
      INTEGER,          INTENT(IN)    :: J             ! Latitude index
      INTEGER,          INTENT(IN)    :: L             ! Level index
      INTEGER,          INTENT(IN)    :: IDX           ! ND38 index
      REAL(fp),         INTENT(IN)    :: F_RAINOUT     ! Fraction of grid box 
                                                       !  undergoing rainout
      REAL(fp),         INTENT(IN)    :: K_RAIN        ! Rainout constant
      REAL(fp),         INTENT(IN)    :: DT            ! Rainout timestep [s]
      CHARACTER(LEN=*), INTENT(IN)    :: ERRMSG        ! Error message
      TYPE(OptInput),   INTENT(IN)    :: Input_Opt     ! Input options
      TYPE(MetState),   INTENT(IN)    :: State_Met     ! Met State object
      
!
! !INPUT/OUTPUT PARAMETERS: 
!
      REAL(fp),         INTENT(INOUT) :: DSpc(:,:,:,:) ! Accumulator array [kg]
      TYPE(ChmState),   INTENT(INOUT) :: State_Chm     ! Chemistry State object
      TYPE(DgnState),   INTENT(INOUT) :: State_Diag    ! Diags State object
!
! !OUTPUT PARAMETERS:
!
      INTEGER,          INTENT(OUT)   :: RC            ! Success or failure?    
! 
! !REMARKS:
!  An IF statement in WETDEP decides if this rainout is to be done (and thus
!  if this routine will be called.  The criteria for rainout is:
!                                                                             .
!     There is rainout if there is new precip formation in the grid box 
!     (i.e. DQRLSAN(I,J,L) > 0) and the fraction of the grid box experiencing 
!     rainout (i.e. F_RAINOUT) is greater than or equal to the fraction of 
!     the grid box directly overhead experiencing precip (i.e. FTOP).
!        -- Helen Amos (9/10/10)
!
! !REVISION HISTORY: 
!  16 Sep 2010 - R. Yantosca - Initial version
!  25 Aug 2014 - M. Sulprizio- Now accept Input_Opt as an argument
!  22 Sep 2015 - E. Lundgren - Input STT is now in kg/m2 (prev kg)
!  24 Sep 2015 - E. Lundgren - Now pass am_I_Root and RC as arguments 
!  29 Sep 2015 - E. Lundgren - Now pass State_Chm as argument and use 
!                              local pointer for STT array [kg/m2] 
!  30 Jun 2016 - R. Yantosca - Replace STT with Spc and DSTT with DSpc
!  05 Jul 2016 - R. Yantosca - Now replace IDWETD with State_Chm%Map_WetDep
!  09 Nov 2017 - R. Yantosca - Now accept State_Diag as an argument
!  09 Nov 2017 - R. Yantosca - Return error condition to calling program
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      ! Scalars
      INTEGER            :: N,        NW
      REAL(fp)           :: RAINFRAC, WETLOSS

      ! Strings
      CHARACTER(LEN=255) :: ErrorMsg, ThisLoc

      ! Pointers
      REAL(fp), POINTER  :: Spc(:,:,:,:)

#if defined( TOMAS )
      ! Scavenging fraction of 30-bin aerosols (win, 7/16/09)
      REAL(fp)           :: TOM_SC_FRACTION(IBINS)
      REAL(fp)           :: SOLFRAC, XFRAC
#endif

      !=================================================================
      ! DO_RAINOUT_ONLY begins here!
      !=================================================================

      ! Initialize
      RC        = GC_SUCCESS
      ErrorMsg  = ''
      ThisLoc   = 
     &   ' -> at DO_RAINOUT_ONLY (in module GeosCore/wetscav_mod.F)'

      ! Point to the chemical species array [kg/m2]
      Spc => State_Chm%Species

#if defined( BPCH_DIAG )
      !-----------------------------------------------------------------
      ! ND16 (bpch) diagnostic
      !-----------------------------------------------------------------

      ! ND16 diagnostic...save F 
      IF ( ND16 > 0 .and. L <= LD16 ) THEN
         AD16(I,J,L,IDX) = AD16(I,J,L,IDX) + F_RAINOUT
         CT16(I,J,L,IDX) = CT16(I,J,L,IDX) + 1 
      ENDIF

      ! ND17 diagnostic...increment counter
      IF ( ND17 > 0 .and. L <= LD17 ) THEN
         CT17(I,J,L,IDX) = CT17(I,J,L,IDX) + 1
      ENDIF
#endif

#if defined( NC_DIAG )
      !-----------------------------------------------------------------
      ! HISTORY (aka netCDF diagnostics)
      !
      ! Archive the fraction of the box that is undergoing large-scale 
      ! precipitation (PrecipFracLs).  This includes contributions
      ! from both rainout and washout.  Here we add the contribution
      ! from rainout.
      !
      ! NOTE: We always assume large-scale precipitation, because
      ! the LS flag is always set to TRUE in the calling routine 
      ! for both GEOS-FP and MERRA-2 meteorology.
      !-----------------------------------------------------------------

      ! NOTE: This diagnostic may need some work
      ! Units: [1]
      IF ( State_Diag%Archive_PrecipFracLS ) THEN
         State_Diag%PrecipFracLS(I,J,L) = 
     &   State_Diag%PrecipFracLS(I,J,L) + F_Rainout
      ENDIF
#endif

      !-----------------------------------------------------------------
      ! Loop over all wet deposition species
      !-----------------------------------------------------------------
      DO NW = 1, State_Chm%nWetDep

         ! Get the wetdep ID from the species ID
         N = State_Chm%Map_WetDep(NW)
                  
         ! Call subroutine RAINOUT to comptue the fraction
         ! of species lost to rainout in grid box (I,J,L) 
         CALL RAINOUT( I,         J,         L,         N, 
     &                 K_RAIN,    DT,        F_RAINOUT, RAINFRAC,
     &                 Input_Opt, State_Met, State_Chm, RC        )

         ! Trap potential errors
         IF ( RC /= GC_SUCCESS ) THEN
            Spc    => NULL()
            ErrorMsg = 'Error encountered in "Rainout"!'
            CALL GC_Error( ErrorMsg, RC, ThisLoc )
            RETURN
         ENDIF

#if defined( TOMAS )
         IF ( id_NK1 > 0 ) THEN
            IF ( N >= id_NK1 .and. N < id_NK1 + IBINS ) THEN  

               CALL GETFRACTION( I, J, L, N, LS, State_Met, XFRAC, 
     &                           SOLFRAC, State_Chm )
 
               RAINFRAC = RAINFRAC * XFRAC * SOLFRAC
             ELSE IF ( N >= id_SF1             .and. 
     &                 N <  id_DUST1 + IBINS ) THEN

               CALL GETFRACTION( I, J, L, N, LS, State_Met, XFRAC,  
     &                           SOLFRAC, State_Chm ) 

               RAINFRAC = RAINFRAC * XFRAC
             ENDIF
          ENDIF
#endif
         ! WETLOSS is the amount of species in grid box per unit area
         ! (I,J,L) that is lost to rainout [kg/m2]
         WETLOSS = Spc(I,J,L,N) * RAINFRAC

         ! Subtract the rainout loss in grid box (I,J,L) from Spc [kg/m2]
         Spc(I,J,L,N) = Spc(I,J,L,N) - WETLOSS

         IF ( L == LLPAR ) THEN
           
            ! DSpc is an accumulator array for rained-out species.  
            ! The species in DSpc are in the liquid phase and will 
            ! precipitate to the levels below until a washout occurs.
            ! Initialize DSpc at (I,J,L=LLPAR) with WETLOSS.
            DSpc(NW,L,I,J) = WETLOSS

         ELSE
                                 
            ! Add to DSpc the species lost to rainout in grid box        
            ! (I,J,L) plus the species lost to rainout from grid box 
            ! (I,J,L+1), which has by now precipitated down into 
            ! grid box (I,J,L).  DSpc will continue to accumulate 
            ! rained out species in this manner until a washout 
            ! event occurs.
            DSpc(NW,L,I,J) = DSpc(NW,L+1,I,J) + WETLOSS

         ENDIF

#if defined( BPCH_DIAG )
         !--------------------------------------------------------------
         ! ND17 and ND39 (bpch) diagnostics
         !--------------------------------------------------------------

         ! ND17 diagnostic...rainout fractions [unitless]
         IF ( ND17 > 0 .and. L <= LD17 ) THEN
            AD17(I,J,L,NW,IDX) = 
     &           AD17(I,J,L,NW,IDX) + RAINFRAC / F_RAINOUT
         ENDIF

         ! ND39 diag -- save rainout losses in [kg/s] 
         ! NOTE: eventually change this to kg/m2/s (ewl, 11/16/15)
         ! Add LGTMM in condition for AD39 (ccc, 11/18/09)
         IF ( ( ND39 > 0 .or. LGTMM ) .and. L <= LD39 ) THEN
            AD39(I,J,L,NW) = AD39(I,J,L,NW) + WETLOSS / DT
     &                       * State_Met%AREA_M2(I,J,1)
         ENDIF
#endif

#if defined( NC_DIAG )
         !--------------------------------------------------------------
         ! HISTORY (aka netCDF diagnostics)
         !
         ! (1) Archive the fraction of soluble species lost to rainout
         !     in large-scale precipitation (RainFracLS)
         !     
         ! (2) Archive the amount of soluble species lost to large-
         !     scale precipitation (WetLossLS).  This includes 
         !     contributions from rainout, washout, and reevaporation.
         !     Here we add the component from rainout.
         !--------------------------------------------------------------

         ! Units: [1]
         IF ( State_Diag%Archive_RainFracLs .and. 
     &        F_Rainout > 0.0_fp ) THEN 
            State_Diag%RainFracLs(I,J,L,NW) = RainFrac / F_Rainout
         ENDIF

         ! Units: [kg/s], but eventually consider changing to [kg/m2/s]
         IF ( State_Diag%Archive_WetLossLs ) THEN
            State_Diag%WetLossLs(I,J,L,NW) =
     &      State_Diag%WetLossLs(I,J,L,NW) + ( WetLoss / DT )
     &                                     *  State_Met%Area_M2(I,J,1)
         ENDIF
#endif

         ! Archive wet loss in kg/m2/s
         IF ( LSOILNOX ) THEN
            CALL SOIL_WETDEP ( I, J, L, N, WETLOSS / DT, State_Chm )
         ENDIF

         !--------------------------------------------------------------
         ! Error checks
         !--------------------------------------------------------------
         IF ( IT_IS_NAN( Spc(I,J,L,N) )  .or.
     &        Spc(I,J,L,N)   < 0e+0_fp   .or.
     &        DSpc(NW,L,I,J) < 0e+0_fp        ) THEN

            ! Print error message
            CALL SAFETY( I, J, L, N, ERRMSG,
     &                   LS          = LS, 
     &                   PDOWN       = PDOWN(L,I,J),
     &                   QQ          = QQ(L,I,J),
     &                   ALPHA       = 0e+0_fp,
     &                   ALPHA2      = 0e+0_fp,
     &                   RAINFRAC    = RAINFRAC,
     &                   WASHFRAC    = 0e+0_fp,
     &                   MASS_WASH   = 0e+0_fp,
     &                   MASS_NOWASH = 0e+0_fp,
     &                   WETLOSS     = WETLOSS,
     &                   GAINED      = 0e+0_fp,
     &                   LOST        = 0e+0_fp,
     &                   DSpc        = DSpc(NW,:,I,J),
     &                   Spc         = Spc(I,J,:,N),
     &                   RC          = RC              )
            
            ! Trap potential errors
            IF ( RC /= GC_SUCCESS ) THEN
               Spc      => NULL()
               ErrorMsg = 'Error encountered in "Safety"!'
               CALL GC_Error( ErrorMsg, RC, ThisLoc )
               RETURN
            ENDIF

         ENDIF
      ENDDO

      ! Free pointer
      Spc => NULL()

      END SUBROUTINE DO_RAINOUT_ONLY
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: do_washout_only
!
! !DESCRIPTION: Subroutine DO\_WASHOUT\_ONLY removes species by washout.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE DO_WASHOUT_ONLY( am_I_Root,  LS,        I, J, L,
     &                            IDX,        ERRMSG,    QDOWN, 
     &                            Q,          F_WASHOUT, F_RAINOUT,
     &                            DT,         PDOWN,     DSpc, 
     &                            Input_Opt,  State_Met, State_Chm,
     &                            State_Diag, RC,        REEVAP     )
!
! !USES:
!
      USE CMN_SIZE_MOD                                 ! Size parameters
#if defined( BPCH_DIAG )
      USE CMN_DIAG_MOD                                 ! Diagnostic flags
      USE DIAG_MOD,       ONLY : AD05
      USE DIAG_MOD,       ONLY : AD16                  ! ND16 diag array
      USE DIAG_MOD,       ONLY : AD17                  ! ND17 diag array
      USE DIAG_MOD,       ONLY : AD18                  ! ND18 diag array
      USE DIAG_MOD,       ONLY : AD39                  ! ND39 diag array
      USE DIAG_MOD,       ONLY : CT16                  ! ND16 diag counter
      USE DIAG_MOD,       ONLY : CT17                  ! ND17 diag counter
      USE DIAG_MOD,       ONLY : CT18                  ! ND18 diag counter
#endif
      USE ErrCode_Mod                                  
      USE ERROR_MOD
      USE GET_NDEP_MOD,   ONLY : SOIL_WETDEP           ! Wet deposited species
      USE Input_Opt_Mod,  ONLY : OptInput              ! Input options
      USE State_Chm_Mod,  ONLY : ChmState              ! Chemistry State object
      USE State_Diag_Mod, ONLY : DgnState
      USE State_Met_Mod,  ONLY : MetState              ! Met State object
#if defined( TOMAS )      
      USE TOMAS_MOD,      ONLY : IBINS, ICOMP, AQOXID
#endif
!
! !INPUT PARAMETERS: 
!
      LOGICAL,          INTENT(IN)    :: am_I_Root     ! Are we on root CPU?  
      LOGICAL,OPTIONAL, INTENT(IN)    :: REEVAP        ! Do re-evaporation?
      LOGICAL,          INTENT(IN)    :: LS            ! =T denotes LS precip
      INTEGER,          INTENT(IN)    :: I             ! Longitude index
      INTEGER,          INTENT(IN)    :: J             ! Latitude index
      INTEGER,          INTENT(IN)    :: L             ! Level index
      INTEGER,          INTENT(IN)    :: IDX           ! ND38 index
      CHARACTER(LEN=*), INTENT(IN)    :: ERRMSG        ! Error message
      REAL(fp),         INTENT(IN)    :: QDOWN         ! Precip leaving thru
                                                       !  bottom of box (I,J,L)
      REAL(fp),         INTENT(IN)    :: Q             ! Precip forming or
                                                       ! evaporating 
                                                       ! in box (I,J,L)
      REAL(fp),         INTENT(IN)    :: F_WASHOUT     ! Fraction of grid box 
                                                       !  undergoing washout
      REAL(fp),         INTENT(IN)    :: F_RAINOUT     ! Fraction of grid box 
                                                       !  undergoing rainout
      REAL(fp),         INTENT(IN)    :: DT            ! Rainout timestep [s]
      REAL(fp),         INTENT(IN)    :: PDOWN(:,:,:)  ! Precip 
      TYPE(OptInput),   INTENT(IN)    :: Input_Opt     ! Input options
      TYPE(MetState),   INTENT(IN)    :: State_Met     ! Met State object
!
! !INPUT/OUTPUT PARAMETERS: 
!
      REAL(fp),         INTENT(INOUT) :: DSpc(:,:,:,:) ! Accumulator array 
                                                       ! [kg/m2]
      TYPE(ChmState),   INTENT(INOUT) :: State_Chm     ! Chemistry State object
      TYPE(DgnState),   INTENT(INOUT) :: State_Diag    ! Diags State object

!
! !OUTPUT PARAMETERS:
!
      INTEGER,          INTENT(OUT)   :: RC            ! Success or failure?
! 
! !REMARKS:
!  A fraction ALPHA of the raindrops falling down from grid 
!  box (I,J,L+1) to grid box (I,J,L) will evaporate along the 
!  way.  ALPHA is given by:
!   
!             precip leaving (I,J,L+1) - precip leaving (I,J,L)
!   ALPHA = ---------------------------------------------------
!                      precip leaving (I,J,L+1)
! 
! 
!                     -QQ(L,I,J) * DZ(I,J,L)
!         =         --------------------------
!                         PDOWN(L+1,I,J)
! 
!  We assume that a fraction ALPHA2 = 0.5 * ALPHA of the 
!  previously rained-out aerosols and HNO3 coming down from 
!  level (I,J,L+1) will evaporate and re-enter the atmosphere 
!  in the gas phase in grid box (I,J,L).  This process is 
!  called "resuspension".  
! 
!  For non-aerosol species, the amount of previously rained 
!  out mass coming down from grid box (I,J,L+1) to grid box 
!  (I,J,L) is figured into the total mass available for 
!  washout in grid box (I,J,L).  We therefore do not have to
!  use the fraction ALPHA2 to compute the resuspension.
! 
!  NOTE from Hongyu Liu about ALPHA (hyl, 2/29/00)
!  =============================================================
!  If our QQ field was perfect, the evaporated amount in grid 
!  box (I,J,L) would be at most the total rain amount coming 
!  from above (i.e. PDOWN(I,J,L+1) ). But this is not true for 
!  the MOISTQ field we are using.  Sometimes the evaporation in 
!  grid box (I,J,L) can be more than the rain amount from above.  
!  The reason is our "evaporation" also includes the effect of 
!  cloud detrainment.  For now we cannot find a way to 
!  distinguish betweeen the two. We then decided to release 
!  aerosols in both the detrained air and the evaporated air. 
! 
!  Therefore, we should use this term in the numerator:
!  
!                 -QQ(I,J,L) * BXHEIGHT(I,J,L) 
! 
!  instead of the term:
!  
!                 PDOWN(L+1)-PDOWN(L)
! 
!  Recall that in make_qq.f we have restricted PDOWN to 
!  positive values, otherwise, QQ would be equal to 
!  PDOWN(L+1)-PDOWN(L).           
!  =============================================================  
!  Update (V. Shah 6/29/16)
!  For GEOS-FP and MERRA2 met fields we use the 
!  following term in the numerator instead:
!		REEVAP(L,I,J) * BXHEIGHT(I,J,L) 
!  =============================================================
!
! !REVISION HISTORY: 
!  16 Sep 2010 - R. Yantosca - Initial version
!  20 Sep 2010 - R. Yantosca - Update definition of ALPHA if we are doing
!                              partial re-evaporation.
!  28 Sep 2010 - H. Amos     - Now check for NaN's with function IT_IS_NAN
!  31 Dec 2010 - H. Amos     - new variable, TK, for temperature
!  26 May 2011 - R. Yantosca - Bug fix: Only apply the error trap for the
!                              condition WASHFRAC < 1d-3 for MERRA met
!  25 May 2011 - Q. Wang     - new variable, TF, for total precip fraction
!  25 May 2011 - Q. Wang     - Also pass F_RAINOUT via the arg list
!  25 May 2011 - Q. Wang     - Correct the washfrac to make sure that washout 
!                              is applied the area of F_washout instead of 
!                              total area of (F_washout+F_rainout)
!  27 May 2011 - R. Yantosca - Added comments, readjusted IF statements to 
!                              avoid floating-point problems
!  16 Aug 2011 - H. Amos     - Replace logical AER with KIN to emphasize that
!                              washout is either modeled as a kinetic process
!                              or an equilibrium process
!  09 Nov 2012 - M. Payer    - Replaced all met field arrays with State_Met
!                              derived type object
!  26 Sep 2013 - R. Yantosca - Renamed GEOS_57 Cpp switch to GEOS_FP
!  25 Aug 2014 - M. Sulprizio- Now accept Input_Opt as an argument
!  04 Feb 2015 - M. Sulprizio- Fix calculation of WETLOSS for non-aerosol 
!                              tracers (C. Friedman)
!  20 May 2015 - M. Sulprizio- Remove WASHFRAC < 1D-3 error trap for MERRA
!                              following recommendation by Viral Shah
!  04 Jun 2015 - E. Lundgren - Now accept am_I_Root and RC as arguments
!  22 Sep 2015 - E. Lundgren - Incoming tracer units are now [kg/m2]
!  29 Sep 2015 - E. Lundgren - Now use local pointer for STT array [kg/m2]
!  30 Jun 2016 - R. Yantosca - Replace STT with Spc and DSTT with DSpc
!  05 Jul 2016 - R. Yantosca - Replace IDWETD with State_Chm%Map_WetDep
!  09 Nov 2017 - R. Yantosca - Now accept State_Diag as an argument
!  09 Nov 2017 - R. Yantosca - Return error condition to calling program
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      ! Scalars
      LOGICAL            :: KIN,       DO_REEVAP
      INTEGER            :: N,         NW
      REAL(fp)           :: ALPHA,     ALPHA2,      GAINED,   LOST
      REAL(fp)           :: MASS_WASH, MASS_NOWASH, WASHFRAC, WETLOSS
      REAL(fp)           :: TK,        TF  

      ! Strings
      CHARACTER(LEN=255) :: ErrorMsg,  ThisLoc

      ! Pointers
      ! We need to define local array to hold species concentraiton values 
      ! from the Chemistry State (State_Chm) object [kg/m2] (ewl, 9/29/15)
      REAL(fp), POINTER  :: Spc(:,:,:,:)

#if defined( TOMAS )
      REAL(fp)           :: REEVAPSO2  !(win, 7/16/09)
      INTEGER            :: KMIN       !(win, 7/16/09)
#endif

      !=================================================================
      ! DO_WASHOUT_ONLY begins here!
      !=================================================================

      ! Initialize
      RC        = GC_SUCCESS
      ErrorMsg  = ''
      ThisLoc   = 
     &   ' -> at Do_Washout_Only (in module GeosCore/wetscav_mod.F)'

      ! Point to the chemical species array [kg/m2]
      Spc => State_Chm%Species

#if defined( BPCH_DIAG )
      !-----------------------------------------------------------------
      ! ND16 (bpch) diagnostic
      !-----------------------------------------------------------------

      ! ND16 diagnostic...save F (fraction of grid box raining)
      IF ( ND16 > 0e+0_fp .and. L <= LD16 ) THEN
         AD16(I,J,L,IDX) = AD16(I,J,L,IDX) + F_WASHOUT
         CT16(I,J,L,IDX) = CT16(I,J,L,IDX) + 1 
      ENDIF

      ! ND18 diagnostic...increment counter
      IF ( ND18 > 0 .and. L <= LD18 ) THEN
         CT18(I,J,L,IDX) = CT18(I,J,L,IDX) + 1
      ENDIF
#endif

#if defined( NC_DIAG )
      !-----------------------------------------------------------------
      ! HISTORY (aka netCDF diagnostics)
      !
      ! Archive the fraction of the box that is undergoing large-scale 
      ! precipitation (PrecipFracLs).  This includes contributions
      ! from both rainout and washout.  Here we add the contribution
      ! from washout.
      !
      ! NOTE: We always assume large-scale precipitation, because
      ! the LS flag is always set to TRUE in the calling routine 
      ! for both GEOS-FP and MERRA-2 meteorology.
      !-----------------------------------------------------------------

      ! NOTE: This diagnostic may need some work
      ! Units: [1]
      IF ( State_Diag%Archive_PrecipFracLS ) THEN
         State_Diag%PrecipFracLS(I,J,L) = 
     &   State_Diag%PrecipFracLS(I,J,L) + F_Washout
      ENDIF
#endif

      ! air temperature [K]
      TK  = State_Met%T(I,J,L)

      ! TOTAL precipitation fraction
      TF  = F_WASHOUT + F_RAINOUT

      !-----------------------------------------------------------------
      ! Loop over all wet deposition species
      !-----------------------------------------------------------------
      DO NW = 1, State_Chm%nWetDep

         ! Get the species ID from the wetdep ID
         N           = State_Chm%Map_WetDep(NW)

         ! zero local variables
         ALPHA       = 0e+0_fp
         ALPHA2      = 0e+0_fp
         WASHFRAC    = 0e+0_fp
         MASS_WASH   = 0e+0_fp
         MASS_NOWASH = 0e+0_fp
         WETLOSS     = 0e+0_fp
         GAINED      = 0e+0_fp
         LOST        = 0e+0_fp

         ! Call WASHOUT to compute the fraction of 
         ! species lost to washout in grid box (I,J,L)
         CALL WASHOUT( am_I_Root, I, J, L, N,
     &                 State_Met%BXHEIGHT(I,J,L),
     &                 TK,        QDOWN,        DT,
     &                 TF,        
     &                 State_Chm%H2O2AfterChem(I,J,L),
     &                 State_Chm%SO2AfterChem(I,J,L),
     &                 WASHFRAC,  KIN,          Input_Opt,
     &                 State_Met, State_Chm,    RC           )

         ! Trap potential errors
         IF ( RC /= GC_SUCCESS ) THEN
            ErrorMsg = 'Error encountered in "Washout"!'
            CALL GC_Error( ErrorMsg, RC, ThisLoc )
            RETURN
         ENDIF

         !%%% BUG FIX (bmy, hamos, 5/26/11, 8/12/15)
         !
         ! Check if WASHFRAC = NaN or WASHFRAC < 0.1 %
         ! 
         ! If WASHFRAC = NaN, then DSpc = NaN and SAFETY trips because 
         ! species concentrations must be finite.WASHFRAC = NaN when F = 0. 
         ! When less than 0.1% of a soluble species is available for washout
         ! DSpc < 0 and SAFETY trips.  (Helen Amos, 20100928)
         !
         ! Viral Shah wrote:
         !  The condition requiring that the WASHFRAC>1D-3 or any number 
         !  greater than 0 for washout to occur, prevents the partial 
         !  resuspension of the dissolved mass falling from above, and leads
         !  to an overestimate in the wet deposited mass... I recommend that 
         !  the condition on WASHFRAC should be removed.  If  this leads to 
         !  negative values of DSpc, we should add a few lines to restrict
         !  DSpc to a minimum of 0. (mps, 5/20/15)
         !
         IF ( IT_IS_NAN( WASHFRAC ) ) THEN
            CYCLE
         ENDIF

         ! Adjust WASHFRAC accordingly for aerosols.  NOTE: TF is always 
         ! > 0 since DO_WASHOUT_ONLY is only called if F_WASHOUT > 0.  
         ! We will never get a div-by-zero error here. (bmy, 5/27/11)
         IF ( KIN ) THEN
           WASHFRAC = WASHFRAC / TF * F_WASHOUT
         ENDIF
            
         !==============================================================      
         ! Washout of aerosol species -- 
         ! this is modeled as a kinetic process
         !==============================================================      
         IF ( KIN ) THEN

            ! Define ALPHA, the fraction of the raindrops that 
            ! re-evaporate when falling from (I,J,L+1) to (I,J,L)
            ALPHA = ( ABS( Q ) * State_Met%BXHEIGHT(I,J,L) * 100e+0_fp )
     &            / ( PDOWN(L+1,I,J)                               )

            ! Restrict ALPHA to be less than 1 (>1 is unphysical)
            ! (hma, 24-Dec-2010)
            IF ( ALPHA > 1e+0_fp ) THEN 
               ALPHA = 1e+0_fp
            ENDIF

            ! ALPHA2 is the fraction of the rained-out aerosols
            ! that gets resuspended in grid box (I,J,L)
            ALPHA2  = 0.5e+0_fp * ALPHA

            ! GAINED is the rained out aerosol coming down from 
            ! grid box (I,J,L+1) that will evaporate and re-enter 
            ! the atmosphere in the gas phase in grid box (I,J,L).
            GAINED  = DSpc(NW,L+1,I,J) * ALPHA2

            ! Amount of aerosol lost to washout in grid box
            ! (qli, bmy, 10/29/02)
            WETLOSS = Spc(I,J,L,N) * WASHFRAC - GAINED

            ! Remove washout losses in grid box (I,J,L) from Spc.
            ! Add the aerosol that was reevaporated in (I,J,L).
            ! SO2 in sulfate chemistry is wet-scavenged on the
            ! raindrop and converted to SO4 by aqeuous chem.
            ! If evaporation occurs then SO2 comes back as SO4
            ! (rjp, bmy, 3/23/03)
            IF ( N == id_SO2 ) THEN
               Spc(I,J,L,id_SO4) = Spc(I,J,L,id_SO4) 
     &                           + GAINED * 96e+0_fp / 64e+0_fp

               Spc(I,J,L,N)      = Spc(I,J,L,N) *
     &                                          ( 1e+0_fp - WASHFRAC )


#if defined( TOMAS )
!added for TOMAS (win, 7/16/09)
#if defined( BPCH_DIAG )
            ! Save the amout of SO4 [kg S/m2] added via aqueous 
            ! chem to ND05(6) diagnostic assuming it's all 
            ! by reacting with H2O2 (win, 7/16/09)
            IF ( ND05 > 0 .and. L <= LD05 ) 
     &         AD05(I,J,L,6) = AD05(I,J,L,6) +
     &         ( GAINED * 32e+0_fp / 64e+0_fp )
#endif

            ! Re-evaporated portion get distributed onto
            ! size-resolved sulfate by AQOXID (win, 7/16/09)
            IF ( GAINED > 0e+0_fp ) THEN 
               IF ( LS ) THEN
! JKodros (6/2/15 - allow for different TOMAS bin lengths)
# if  defined( TOMAS12)
                  KMIN = 5
# elif  defined( TOMAS15)
                  KMIN = 8
# elif  defined( TOMAS30)
                  KMIN = 10
# else
                  KMIN = 20
# endif

               ELSE
# if  defined( TOMAS12)
                  KMIN = 3
# elif  defined( TOMAS15)
                  KMIN = 6
# elif  defined( TOMAS30)
                  KMIN = 6
# else
                  KMIN = 16
# endif
               ENDIF

               ! ***NOTE*** Species concentration units are currently in
               ! [kg/m2] which is incompatible with TOMAS. Units are
               ! therefore converted to [kg] locally within AQOXID.
               ! GAINED is now [kg/m2] ans so is multiplied
               ! by area prior to passing REEVAPSO2 to AQOXID (ewl, 9/30/15)
               IF ( TRIM( State_Chm%Spc_Units ) .eq. 'kg/m2' ) THEN
                  REEVAPSO2 = GAINED * 96e+0_fp / 64e+0_fp
     &                        * State_Met%AREA_M2(I,J,1)
               ELSE
                  ErrorMsg= 'Unexpected species units: ' // 
     &                       TRIM( State_Chm%Spc_Units )
                  CALL GC_Error( ErrorMsg, RC, ThisLoc )
                  RETURN
               ENDIF
               CALL AQOXID( am_I_Root, REEVAPSO2, KMIN, I, J, L, 
     &                      Input_Opt, State_Met, State_Chm, RC  )
            ENDIF
!end -added for TOMAS  (win, 7/16/09)
#endif

            ELSE
               Spc(I,J,L,N)      = Spc(I,J,L,N) - WETLOSS
            ENDIF

            ! LOST is the rained out aerosol coming down from
            ! grid box (I,J,L+1) that will remain in the liquid
            ! phase in grid box (I,J,L) and will NOT re-evaporate.
            LOST = DSpc(NW,L+1,I,J) - GAINED

            ! Add the washed out species from grid box (I,J,L) to 
            ! DSpc.  Also add the amount of species coming down
            ! from grid box (I,J,L+1) that does NOT re-evaporate.
            IF ( F_RAINOUT > 0e+0_fp ) THEN 
               DSpc(NW,L,I,J) = DSpc(NW,L,  I,J) + WETLOSS
            ELSE
               DSpc(NW,L,I,J) = DSpc(NW,L+1,I,J) + WETLOSS
            ENDIF

#if defined( BPCH_DIAG )
            !-----------------------------------------------------------
            ! ND18 (bpch) diagnostic
            !-----------------------------------------------------------

            ! ND18 diagnostic...divide washout fraction by F
            IF ( ND18 > 0 .and. L <= LD18 ) THEN
               AD18(I,J,L,NW,IDX) = 
     &         AD18(I,J,L,NW,IDX) + ( WASHFRAC / F_WASHOUT )
            ENDIF
#endif

#if defined( NC_DIAG )
            !-----------------------------------------------------------
            ! HISTORY (aka netCDF diagnostics)
            !
            ! Archive the fraction of soluble species lost to washout
            ! in large-scale precipitation (WashFracLS)
            !
            ! Here we only handle the soluble aerosol species
            !-----------------------------------------------------------
            
            ! Units: [1]
            IF ( State_Diag%Archive_WashFracLS .and. 
     &           F_Washout > 0.0_fp ) THEN
               State_Diag%WashFracLS(I,J,L,NW) = WashFrac / F_Washout
            ENDIF
#endif

         !==============================================================      
         ! Washout of non-aerosol species
         ! This is modeled as an equilibrium process
         !==============================================================      
         ELSE
                  
            ! MASS_NOWASH is the amount of non-aerosol species in 
            ! grid box (I,J,L) that is NOT available for washout.
            MASS_NOWASH = ( 1e+0_fp - F_WASHOUT ) * Spc(I,J,L,N)
               
            ! MASS_WASH is the total amount of non-aerosol species
            ! that is available for washout in grid box (I,J,L).
            ! It consists of the mass in the precipitating
            ! part of box (I,J,L), plus the previously rained-out
            ! species coming down from grid box (I,J,L+1).
            ! (Eq. 15, Jacob et al, 2000).
            MASS_WASH = ( F_WASHOUT*Spc(I,J,L,N) ) + DSpc(NW,L+1,I,J)

            ! WETLOSS is the amount of species mass in 
            ! grid box (I,J,L) that is lost to washout.
            ! (Eq. 16, Jacob et al, 2000)
            WETLOSS = ( MASS_WASH - DSpc(NW,L+1,I,J) ) * WASHFRAC

            ! The species left in grid box (I,J,L) is what was
            ! in originally in the non-precipitating fraction 
            ! of the box, plus MASS_WASH, less WETLOSS. 
            Spc(I,J,L,N) = Spc(I,J,L,N) - WETLOSS  
            
            ! Add washout losses in grid box (I,J,L) to DSpc 
            IF ( F_RAINOUT > 0e+0_fp ) THEN
               DSpc(NW,L,I,J) = DSpc(NW,L,  I,J) + WETLOSS
            ELSE
               DSpc(NW,L,I,J) = DSpc(NW,L+1,I,J) + WETLOSS
            ENDIF

#if defined( BPCH_DIAG )
            !-----------------------------------------------------------
            ! ND18 (bpch) diagnostic
            !-----------------------------------------------------------

            ! ND18 diagnostic...we don't have to divide the
            ! washout fraction by F since this is accounted for.
            IF ( ND18 > 0 .and. L <= LD18 ) THEN
               AD18(I,J,L,NW,IDX) = 
     &         AD18(I,J,L,NW,IDX) + WASHFRAC
            ENDIF
#endif

#if defined( NC_DIAG )
            !-----------------------------------------------------------
            ! HISTORY (aka netCDF diagnostics)
            !
            ! Archive the fraction of soluble species lost to washout
            ! in large-scale precipitation (WashFracLS)
            !
            ! Here we handle the non-aerosol soluble species.  We don't 
            ! have to divide by F_Washout, since this has already been 
            ! accounted for in the equations above.
            !-----------------------------------------------------------
            
            ! Units: [1]
            IF ( State_Diag%Archive_WashFracLS ) THEN
               State_Diag%WashFracLS(I,J,L,NW) = WashFrac
            ENDIF
#endif

         ENDIF

#if defined( BPCH_DIAG )
         !--------------------------------------------------------------
         ! ND39 (bpch) diagnostic
         !--------------------------------------------------------------

         ! ND39 diag -- save washout losses in [kg/s]
         ! NOTE: Eventually this will be [kg/m2/s] (ewl, 11/16/15)
         ! Add LGTMM in condition for AD39 (ccc, 11/18/09)
         IF ( ( ND39 > 0 .or. LGTMM ) .and. L <= LD39 ) THEN
            AD39(I,J,L,NW) = AD39(I,J,L,NW) + WETLOSS / DT
     &                       * State_Met%AREA_M2(I,J,1)

         ENDIF
#endif

#if defined( NC_DIAG )
         !--------------------------------------------------------------
         ! HISTORY (aka netCDF diagnostics)
         !
         ! Archive the amount of soluble species lost to large-scale 
         ! precipitation (WetLossLS).  This includes contributions 
         ! from rainout, washout, and reevaporation.  Here we add the 
         ! component from washout.
         !--------------------------------------------------------------

         ! Units: [kg/s], but eventually consider changing to [kg/m2/s]
         IF ( State_Diag%Archive_WetLossLs ) THEN
            State_Diag%WetLossLS(I,J,L,NW) = 
     &      State_Diag%WetLossLS(I,J,L,NW) + ( WetLoss / DT )
     &                                     *  State_Met%Area_M2(I,J,1)
         ENDIF
#endif

         ! Archive wet loss in kg/m2/s
         IF ( LSOILNOX ) THEN
            CALL SOIL_WETDEP ( I, J, L, N, WETLOSS / DT, State_Chm )
         ENDIF

         !--------------------------------------------------------------
         ! Error checks
         !--------------------------------------------------------------
         IF ( IT_IS_NAN( Spc(I,J,L,N) )       .or.
     &        Spc(I,J,L,N)   < 0e+0_fp        .or. 
     &        DSpc(NW,L,I,J) < 0e+0_fp      ) THEN

            ! Print error message and stop simulaton
            CALL SAFETY( I, J, L, N, ERRMSG,
     &                   LS          = LS,
     &                   PDOWN       = PDOWN(L+1,I,J),
     &                   QQ          = QQ(L,I,J),
     &                   ALPHA       = ALPHA,
     &                   ALPHA2      = ALPHA2,
     &                   RAINFRAC    = 0e+0_fp,
     &                   WASHFRAC    = WASHFRAC,
     &                   MASS_WASH   = MASS_WASH,
     &                   MASS_NOWASH = MASS_NOWASH,
     &                   WETLOSS     = WETLOSS,
     &                   GAINED      = GAINED,
     &                   LOST        = LOST,
     &                   DSpc        = DSpc(NW,:,I,J),
     &                   Spc         = Spc(I,J,:,N),
     &                   RC          = RC              )

            ! Trap potential errors
            IF ( RC /= GC_SUCCESS ) THEN
               Spc    => NULL()
               ErrorMsg = 'Error encountered in "Safety"!'
               CALL GC_Error( ErrorMsg, RC, ThisLoc )
               RETURN
            ENDIF
         ENDIF
      ENDDO  

      ! Free pointer
      Spc => NULL()
 
      END SUBROUTINE DO_WASHOUT_ONLY
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: do_complete_reevap
!
! !DESCRIPTION: Subroutine DO\_COMPLETE\_REEVAP re-evaporates all of the
!  soluble species back into the atmosphere.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE DO_COMPLETE_REEVAP( am_I_Root, LS,         I, J, L,
     &                               IDX,       ERRMSG,     DT,     
     &                               DSpc,      Input_Opt,  State_Met,
     &                               State_Chm, State_Diag, RC         ) 
!     
! !USES:
!
      USE CMN_SIZE_MOD                                 ! Size parameters
#if defined( BPCH_DIAG )                              
      USE CMN_DIAG_MOD                                 ! Diagnostic flags
      USE DIAG_MOD,       ONLY : AD05
      USE DIAG_MOD,       ONLY : AD16                  ! ND16 diag array
      USE DIAG_MOD,       ONLY : AD17                  ! ND17 diag array
      USE DIAG_MOD,       ONLY : AD18                  ! ND18 diag array
      USE DIAG_MOD,       ONLY : AD39                  ! ND39 diag array
      USE DIAG_MOD,       ONLY : CT16                  ! ND16 diag counter
      USE DIAG_MOD,       ONLY : CT17                  ! ND17 diag counter
      USE DIAG_MOD,       ONLY : CT18                  ! ND18 diag counter
#endif                                                
      USE ErrCode_Mod
      USE Error_Mod,      ONLY : It_Is_Nan
      USE GET_NDEP_MOD,   ONLY : SOIL_WETDEP           ! Wet deposited species
      USE Input_Opt_Mod,  ONLY : OptInput
      USE State_Chm_Mod,  ONLY : ChmState
      USE State_Diag_Mod, ONLY : DgnState
      USE State_Met_Mod,  ONLY : MetState
#if defined( TOMAS )
      USE TOMAS_MOD,      ONLY : IBINS, ICOMP, AQOXID
#endif
!
! !INPUT PARAMETERS: 
!
      LOGICAL,          INTENT(IN)    :: am_I_Root     ! Are we on root CPU?
      LOGICAL,          INTENT(IN)    :: LS            ! =T denotes LS precip
      INTEGER,          INTENT(IN)    :: I             ! Longitude index
      INTEGER,          INTENT(IN)    :: J             ! Latitude index
      INTEGER,          INTENT(IN)    :: L             ! Level index
      INTEGER,          INTENT(IN)    :: IDX           ! ND38 index
      CHARACTER(LEN=*), INTENT(IN)    :: ERRMSG        ! Error message
      REAL(fp),         INTENT(IN)    :: DT            ! Rainout timestep [s]
      TYPE(OptInput),   INTENT(IN)    :: Input_Opt     ! Input options
      TYPE(MetState),   INTENT(IN)    :: State_Met     ! Meteorology State obj
!
! !INPUT/OUTPUT PARAMETERS: 
!
      TYPE(ChmState),   INTENT(INOUT) :: State_Chm     ! Chemistry State object
      TYPE(DgnState),   INTENT(INOUT) :: State_Diag    ! Diagnoss State object

      REAL(fp),         INTENT(INOUT) :: DSpc(:,:,:,:) ! Accumulator array 
                                                       ! [kg/m2] 
!
! !OUTPUT PARAMETERS:
!
      INTEGER,          INTENT(OUT)   :: RC            ! Success or failure?
!
! !REVISION HISTORY: 
!  16 Sep 2010 - R. Yantosca - Initial version
!  22 Sep 2015 - E. Lundgren - Incoming tracer units are now kg/m2 (prev kg)
!  24 Sep 2015 - E. Lundgren - Now pass am_I_Root and RC as arguments and  
!                              local pointer for STT array [kg/m2] 
!  30 Sep 2015 - E. Lundgren - Now pass Input_Opt for use in TOMAS AQOXID
!  30 Jun 2016 - R. Yantosca - Replace STT with Spc and DSTT with DSpc
!  05 Jul 2016 - R. Yantosca - Replace IDWETD with State_Chm%Map_WetDep
!  09 Nov 2017 - R. Yantosca - Now accept State_Diag as an argument
!  09 Nov 2017 - R. Yantosca - Return error condition to calling program
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      ! Scalars
      INTEGER            :: N, NW
      REAL(fp)           :: WETLOSS

      ! Strings
      CHARACTER(LEN=255) :: ErrorMsg, ThisLoc

      ! Pointers
      REAL(fp), POINTER  :: Spc(:,:,:,:)

#if defined( TOMAS )
      REAL(fp)           :: REEVAPSO2  !(win, 7/16/09)
      INTEGER            :: KMIN       !(win, 7/16/09)
#endif

      !=================================================================
      ! DO_COMPLETE_REEVAP begins here!
      !=================================================================

      ! Initialize
      RC        = GC_SUCCESS
      ErrorMsg  = ''
      ThisLoc   = 
     &   ' -> at Do_Complete_Reevap (in module GeosCore/wetscav_mod.F)'

      ! Point to chemical species array [kg/m2]
      Spc => State_Chm%Species

      ! Loop over wetdep species
      DO NW = 1, State_Chm%nWetDep
 
         ! Get the species ID from the wetdep ID
         N = State_Chm%Map_WetDep(NW)

         ! WETLOSS is the amount of species in grid box (I,J,L) per area
         ! that is lost to rainout. (qli, bmy, 10/29/02)
         WETLOSS = -DSpc(NW,L+1,I,J)

         ! All of the rained-out species coming from grid box
         ! (I,J,L+1) goes back into the gas phase at (I,J,L)
         ! In evap, SO2 comes back as SO4 (rjp, bmy, 3/23/03)
         IF ( N == id_SO2 ) THEN
            Spc(I,J,L,id_SO4) = Spc(I,J,L,id_SO4) 
     &                        - ( WETLOSS * 96e+0_fp / 64e+0_fp )


#if defined( TOMAS )
!added for TOMAS (win, 7/16/09)
                     ! Save the amout of SO4 [kg S] added via aqueous 
                     ! chem to ND05(6) diagnostic assuming it's all 
                     ! by reacting with H2O2 (win, 7/16/09)

      !=================================================================
      ! sfarina - commenting out this DIAG for now... unclear if this 
      !           was even correct in the older verions because of an
      !           error in the calculation of GAINED
      !=================================================================
!                     IF ( ND05 > 0 .and. L <= LD05 ) 
!     &                    AD05(I,J,L,6) = AD05(I,J,L,6) +
!     &                                    ( GAINED * 32D0 / 64D0 )
                            
                       ! Re-evaporated portion get distributed onto
                     ! size-resolved sulfate by AQOXID (win, 7/16/09)
                     IF ( ABS(WETLOSS) > 0e+0_fp ) THEN 
                        IF ( LS ) THEN
! JKodros (6/2/15) - Allow for different TOMAS bin lengths
# if  defined( TOMAS12)
                           KMIN = 5
# elif  defined( TOMAS15)
                           KMIN = 8
# elif  defined( TOMAS30)
                           KIMIN = 10
# else
                           KMIN = 20
# endif

                        ELSE
# if  defined( TOMAS12)
                           KMIN = 3
# elif  defined( TOMAS15)
                           KMIN = 6
# elif  defined( TOMAS30)
                           KIMIN = 6
# else
                           KMIN = 16
# endif
                        ENDIF

               ! ***NOTE*** Species concentration units are currently in
               ! [kg/m2] which is incompatible with TOMAS. Units are
               ! therefore converted to [kg] locally within AQOXID.
               ! WETLOSS is now [kg/m2] and so is multiplied
               ! by area prior to passing REEVAPSO2 to AQOXID (ewl, 9/30/15)
                        IF ( TRIM( State_Chm%Spc_Units ) .eq.
     &                                          'kg/m2' ) THEN
                           REEVAPSO2 =  - ( WETLOSS * 96e+0_fp 
     &                                  / 64e+0_fp )
     &                                  * State_Met%AREA_M2(I,J,1)
                        ELSE
                           Spc    => NULL()
                           ErrorMsg = 'Unexpected species units: ' 
     &                              // TRIM(State_Chm%Spc_Units)
                           CALL GC_Error( ErrorMsg, RC, ThisLoc )
                           RETURN
                        ENDIF
                        CALL AQOXID( am_I_Root, REEVAPSO2, KMIN, I, J, 
     &                               L, Input_Opt, State_Met, 
     &                               State_Chm, RC )
                     ENDIF
!end- added for TOMAS (win, 7/16/09)
#endif


         ELSE
            Spc(I,J,L,N)      = Spc(I,J,L,N) - WETLOSS
         ENDIF

         ! There is nothing rained out/washed out in grid box
         ! (I,J,L), so set DSpc at grid box (I,J,L) to zero.
         DSpc(NW,L,I,J) = 0e+0_fp
         
#if defined( BPCH_DIAG )
         !--------------------------------------------------------------
         ! ND39 (bpch) diagnostic
         !--------------------------------------------------------------

         ! ND39 diag -- save rainout losses in [kg/s]
         ! NOTE: Eventually this will be [kg/m2/s] (ewl, 11/16/15)
         ! Add LGTMM in condition for AD39 (ccc, 11/18/09)
         IF ( ( ND39 > 0 .or. LGTMM ) .and. L <= LD39 ) THEN
            AD39(I,J,L,NW) = AD39(I,J,L,NW) + ( WETLOSS / DT )
     &                       * State_Met%AREA_M2(I,J,1)

         ENDIF
#endif

#if defined( NC_DIAG )
         !--------------------------------------------------------------
         ! HISTORY (aka netCDF diagnostics)
         !
         ! Archive the amount of soluble species lost to large-scale 
         ! precipitation (WetLossLs).  This includes contributions 
         ! from rainout, washout, and reevaporation.  Here we add the 
         ! component from reevaporation (which is negative).
         !--------------------------------------------------------------

         ! Units: [kg/s], but eventually consider changing to [kg/m2/s]
         IF ( State_Diag%Archive_WetLossLs ) THEN
            State_Diag%WetLossLs(I,J,L,NW) = 
     &      State_Diag%WetLossLs(I,J,L,NW) + ( WetLoss / DT )
     &                                     *  State_Met%Area_M2(I,J,1)
         ENDIF
#endif

         ! Archive wet loss in kg/m2/s
         IF ( LSOILNOX ) THEN
            CALL SOIL_WETDEP ( I, J, L, N, WETLOSS / DT, State_Chm )
         ENDIF

         !--------------------------------------------------------------
         ! Error checks
         !--------------------------------------------------------------
         IF ( IT_IS_NAN( Spc(I,J,L,N) )   .or.
     &        Spc(I,J,L,N)   < 0e+0_fp        .or. 
     &        DSpc(NW,L,I,J) < 0e+0_fp      ) THEN

            ! Print error message and stop simulaton
            CALL SAFETY( I, J, L, N, ERRMSG,
     &                   LS          = LS,
     &                   PDOWN       = 0e+0_fp,
     &                   QQ          = 0e+0_fp,
     &                   ALPHA       = 0e+0_fp,
     &                   ALPHA2      = 0e+0_fp,
     &                   RAINFRAC    = 0e+0_fp,     
     &                   WASHFRAC    = 0e+0_fp,
     &                   MASS_WASH   = 0e+0_fp,
     &                   MASS_NOWASH = 0e+0_fp,
     &                   WETLOSS     = WETLOSS,
     &                   GAINED      = 0e+0_fp,
     &                   LOST        = 0e+0_fp,
     &                   DSpc        = DSpc(NW,:,I,J),
     &                   Spc         = Spc(I,J,:,N),
     &                   RC          = RC              )

            ! Trap potential errors
            IF ( RC /= GC_SUCCESS ) THEN
               Spc      => NULL()
               ErrorMsg = 'Error encountered in "Safety"!'
               CALL GC_Error( ErrorMsg, RC, ThisLoc )
               RETURN
            ENDIF
         ENDIF
      ENDDO

      ! Free pointer
      Spc => NULL()

      END SUBROUTINE DO_COMPLETE_REEVAP
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: do_washout_at_sfc
!
! !DESCRIPTION: Subroutine DO\_WASHOUT\_AT\_SFC washes out the species 
!  at the surface.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE DO_WASHOUT_AT_SFC( am_I_Root,  LS,        I, J, L,
     &                              IDX,        ERRMSG,    QDOWN,
     &                              F,          DT,        DSpc,
     &                              Input_Opt,  State_Met, State_Chm,
     &                              State_Diag, RC                    )
!
! !USES:
!
      USE CMN_SIZE_MOD                                 ! Size parameters
#if defined( BPCH_DIAG )                               
      USE CMN_DIAG_MOD                                 ! Diagnostic flags
      USE DIAG_MOD,       ONLY : AD16                  ! ND16 diag array
      USE DIAG_MOD,       ONLY : AD17                  ! ND17 diag array
      USE DIAG_MOD,       ONLY : AD18                  ! ND18 diag array
      USE DIAG_MOD,       ONLY : AD39                  ! ND39 diag array
      USE DIAG_MOD,       ONLY : CT16                  ! ND16 diag counter
      USE DIAG_MOD,       ONLY : CT17                  ! ND17 diag counter
      USE DIAG_MOD,       ONLY : CT18                  ! ND18 diag counter
#endif                                                 
      USE ErrCode_Mod                                  
      USE ERROR_MOD,      ONLY : IT_IS_NAN             ! Test for NaN
      USE GET_NDEP_MOD,   ONLY : SOIL_WETDEP           ! Wet deposited species
      USE Input_Opt_Mod,  ONLY : OptInput              ! Input options
      USE State_Chm_Mod,  ONLY : ChmState              ! Chm State object
      USE State_Diag_Mod, ONLY : DgnState
      USE State_Met_Mod,  ONLY : MetState              ! Met State object
!
! !INPUT PARAMETERS: 
!     
      LOGICAL,          INTENT(IN)    :: am_I_Root     ! Are we on the root CPU?
      LOGICAL,          INTENT(IN)    :: LS            ! =T denotes LS precip
      INTEGER,          INTENT(IN)    :: I             ! Longitude index
      INTEGER,          INTENT(IN)    :: J             ! Latitude index
      INTEGER,          INTENT(IN)    :: L             ! Level index
      INTEGER,          INTENT(IN)    :: IDX           ! ND38 index
      CHARACTER(LEN=*), INTENT(IN)    :: ERRMSG        ! Error message
      REAL(fp),         INTENT(IN)    :: QDOWN         ! Precip leaving thru
                                                       !  bottom of box (I,J,L)
      REAL(fp),         INTENT(IN)    :: F             ! Fraction of grid box 
                                                       !  undergoing precip
      REAL(fp),         INTENT(IN)    :: DT            ! Rainout timestep [s]
      TYPE(OptInput),   INTENT(IN)    :: Input_Opt     ! Input options
      TYPE(MetState),   INTENT(IN)    :: State_Met     ! Met State object
!
! !INPUT/OUTPUT PARAMETERS: 
!
      REAL(fp),         INTENT(INOUT) :: DSpc(:,:,:,:) ! Accumulator array 
                                                       ! [kg/m2]
      TYPE(ChmState),   INTENT(INOUT) :: State_Chm     ! Chemistry State object
      TYPE(DgnState),   INTENT(INOUT) :: State_Diag    ! Diags State object
! 
! !OUTPUT PARAMETERS:
!
      INTEGER,          INTENT(OUT)   :: RC            ! Success or failure?
!
! !REMARKS:
!  Assume all of the species precipitating down from grid box (I,J,L=2) to 
!  grid box (I,J,L=1) gets washed out in grid box (I,J,L=1).
!
! !REVISION HISTORY: 
!  16 Sep 2010 - R. Yantosca - Initial version
!  09 Nov 2012 - M. Payer    - Replaced all met field arrays with State_Met
!                              derived type object
!  25 Aug 2014 - M. Sulprizio- Now accept Input_Opt as an argument
!  04 Jun 2015 - E. Lundgren - Now accept am_I_Root and RC as arguments
!  22 Sep 2015 - E. Lundgren - Tracer units are now kg/m2 (prev kg)
!  29 Sep 2015 - E. Lundgren - Now use local pointer for STT array [kg/m2]
!  30 Jun 2016 - R. Yantosca - Replace STT with Spc and DSTT with DSpc
!  05 Jul 2016 - R. Yantosca - Now replace IDWETD with State_Chm%Map_Wetdep
!  09 Nov 2017 - R. Yantosca - Now accept State_Diag as an argument
!  09 Nov 2017 - R. Yantosca - Return error condition to calling program
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!     
      ! Scalars
      LOGICAL            :: KIN
      INTEGER            :: N,        NW
      REAL(fp)           :: WASHFRAC, WETLOSS, TMP, TK

      ! Strings
      CHARACTER(LEN=255) :: ErrorMsg, ThisLoc

      ! Pointers
      REAL(fp), POINTER  :: Spc(:,:,:,:)

      !=================================================================
      ! DO_WASHOUT_AT_SFC begins here!
      !=================================================================

      ! Initialize
      RC        = GC_SUCCESS
      ErrorMsg  = ''
      ThisLoc   = 
     &   ' -> at Do_Washout_at_Sfc (in module GeosCore/wetscav_mod.F)'

      ! Point to the chemical species array [kg/m2]
      Spc => State_Chm%Species

#if defined( BPCH_DIAG )
      !-----------------------------------------------------------------
      ! ND16 (bpch) diagnostic
      !-----------------------------------------------------------------

      ! ND16 diagnostic...save F 
      IF ( ND16 > 0 .and. L <= LD16 ) THEN
         AD16(I,J,L,IDX) = AD16(I,J,L,IDX) + F
         CT16(I,J,L,IDX) = CT16(I,J,L,IDX) + 1
      ENDIF

      ! ND18 diagnostic...increment counter
      IF ( ND18 > 0 .and. L <= LD18 ) THEN
         CT18(I,J,L,IDX) = CT18(I,J,L,IDX) + 1
      ENDIF
#endif

#if defined( NC_DIAG )
      !-----------------------------------------------------------------
      ! HISTORY (aka netCDF diagnostics)
      !
      ! Archive the fraction of the box that is undergoing large-scale 
      ! precipitation (PrecipFracLs).  This includes contributions
      ! from both rainout and washout.  Here we add the contribution
      ! from washout.
      !
      ! NOTE: We always assume large-scale precipitation, because
      ! the LS flag is always set to TRUE in the calling routine 
      ! for both GEOS-FP and MERRA-2 meteorology.
      !-----------------------------------------------------------------

      ! NOTE: This diagnostic may need some work
      ! Units: [1]
      IF ( State_Diag%Archive_PrecipFracLS ) THEN
         State_Diag%PrecipFracLS(I,J,L) =
     &   State_Diag%PrecipFracLS(I,J,L) + F

      ENDIF
#endif

      ! air temperature [K]
      TK = State_Met%T(I,J,L)

      !-----------------------------------------------------------------
      ! Loop over all wet deposition species
      !-----------------------------------------------------------------
      DO NW = 1, State_Chm%nWetDep

         ! Get species ID from wetdep ID
         N = State_Chm%Map_WetDep(NW)

         ! Call WASHOUT to compute the fraction of species 
         ! in grid box (I,J,L) that is lost to washout.  
         CALL WASHOUT( am_I_Root, I, J, L, N,
     &                 State_Met%BXHEIGHT(I,J,L), 
     &                 TK ,       QDOWN,        DT,
     &                 F,
     &                 State_Chm%H2O2AfterChem(I,J,L),
     &                 State_Chm%SO2AfterChem(I,J,L),
     &                 WASHFRAC,  KIN,          Input_Opt,
     &                 State_Met, State_Chm,    RC           )

         ! Trap potential errors
         IF ( RC /= GC_SUCCESS ) THEN
            Spc    => NULL()
            ErrorMsg = 'Error encountered in "Washout"!'
            CALL GC_Error( ErrMsg, RC, ThisLoc )
            RETURN
         ENDIF

         ! NOTE: for HNO3 and aerosols, there is an F factor
         ! already present in WASHFRAC.  For other soluble
         ! gases, we need to multiply by the F (hyl, bmy, 10/27/00)
         IF ( KIN ) THEN
            WETLOSS = Spc(I,J,L,N) * WASHFRAC
         ELSE
            WETLOSS = Spc(I,J,L,N) * WASHFRAC * F
         ENDIF

         ! Subtract WETLOSS from Spc [kg/m2]
         Spc(I,J,L,N) = Spc(I,J,L,N) - WETLOSS     
      
         ! Add washout losses in grid box (I,J,L=1) to DSpc [kg/m2] 
         ! (added cdh, 4/14/2009)
         DSpc(NW,L,I,J) = DSpc(NW,L+1,I,J) + WETLOSS

#if defined( BPCH_DIAG )
         !--------------------------------------------------------------
         ! ND18 and ND39 (bpch) diagnostics
         !--------------------------------------------------------------

         ! ND18 diagnostic...LS and conv washout fractions [unitless]
         IF ( ND18 > 0 .and. L <= LD18 ) THEN

            ! Only divide WASHFRAC by F for aerosols, since
            ! for non-aerosols this is already accounted for
            IF ( KIN ) THEN
               TMP = WASHFRAC / F
            ELSE
               TMP = WASHFRAC
            ENDIF
            
            AD18(I,J,L,NW,IDX) = AD18(I,J,L,NW,IDX) + TMP
         ENDIF

         ! ND39 diag -- save washout loss in [kg/s]
         ! NOTE: Eventually this will be [kg/m2/s] (ewl, 11/16/15)
         ! Add LGTMM in condition for AD39 (ccc, 11/18/09)
         IF ( ( ND39 > 0 .or. LGTMM ) .and. L <= LD39 ) THEN
            AD39(I,J,L,NW) = AD39(I,J,L,NW) + WETLOSS / DT
     &                       * State_Met%AREA_M2(I,J,1)
         ENDIF
#endif

#if defined( NC_DIAG )
         !--------------------------------------------------------------
         ! HISTORY (aka netCDF diagnostics)
         !
         ! Archive the fraction of soluble species lost to washout
         ! in large-scale precipitation (WashFracLS)
         !--------------------------------------------------------------
         IF ( State_Diag%Archive_WashFracLS ) THEN

            ! Only divide WASHFRAC by F for aerosols, since
            ! for non-aerosols this is already accounted for
            IF ( KIN ) THEN
               Tmp = WashFrac / F
            ELSE
               TMP = WashFrac
            ENDIF

            ! Units: [1]
            State_Diag%WashFracLS(I,J,L,NW) = Tmp
         ENDIF

         !--------------------------------------------------------------
         ! HISTORY (aka netCDF diagnostics)
         !
         ! Archive the amount of soluble species lost to large-scale 
         ! precipitation (WetLossLS).  This includes contributions 
         ! from rainout, washout, and reevaporation.  Here we add the 
         ! component from washout.
         !--------------------------------------------------------------

         ! Units: [kg/s], but eventually consider changing to [kg/m2/s]
         IF ( State_Diag%Archive_WetLossLS ) THEN
            State_Diag%WetLossLS(I,J,L,NW) = 
     &      State_Diag%WetLossLS(I,J,L,NW) + ( WetLoss / DT )
     &                                     *  State_Met%Area_M2(I,J,1)
         ENDIF

#endif

         ! Archive wet loss in kg/m2/s (check source code for this routine
         ! - ewl )
!         IF ( LSOILNOX ) THEN
            CALL SOIL_WETDEP ( I, J, L, N, WETLOSS / DT, State_Chm )
!         ENDIF

         !-----------------------------------------------------
         ! Dirty kludge to prevent wet deposition from removing 
         ! stuff from stratospheric boxes -- this can cause 
         ! negative species (rvm, bmy, 6/21/00)
         !
         IF ( Spc(I,J,L,N) < 0e+0_fp .and. L > 23 ) THEN
             WRITE ( 6, 101 ) I, J, L, N, 7
 101         FORMAT( 'WETDEP - Spc < 0 at ', 3i4, 
     &               ' for species ', i4, 'in area ', i4 )
             PRINT*, 'Spc:', Spc(I,J,:,N)
             Spc(I,J,L,N) = 0e+0_fp
         ENDIF
         !-----------------------------------------------------

         !--------------------------------------------------------------
         ! Error checks
         !--------------------------------------------------------------
         IF ( IT_IS_NAN( Spc(I,J,L,N) )   .or.
     &        Spc(I,J,L,N)   < 0e+0_fp        .or. 
     &        DSpc(NW,L,I,J) < 0e+0_fp      ) THEN

            !PRINT*, 'WASHFRAC = ', WASHFRAC
            !PRINT*, 'F        = ', F

            ! Print error message and stop simulaton
            CALL SAFETY( I, J, L, N, ERRMSG,
     &                   LS          = LS,
     &                   PDOWN       = 0e+0_fp,
     &                   QQ          = 0e+0_fp,
     &                   ALPHA       = 0e+0_fp,
     &                   ALPHA2      = 0e+0_fp,
     &                   RAINFRAC    = 0e+0_fp,
     &                   WASHFRAC    = 0e+0_fp, 
     &                   MASS_WASH   = 0e+0_fp,
     &                   MASS_NOWASH = 0e+0_fp,
     &                   WETLOSS     = WETLOSS,
     &                   GAINED      = 0e+0_fp,
     &                   LOST        = 0e+0_fp,
     &                   DSpc        = DSpc(NW,:,I,J),
     &                   Spc         = Spc(I,J,:,N),
     &                   RC          = RC              )

            ! Trap potential errors
            IF ( RC /= GC_SUCCESS ) THEN
               Spc      => NULL()
               ErrorMsg = 'Error encountered in "Safety"!'
               CALL GC_Error( ErrorMsg, RC, ThisLoc )
               RETURN
            ENDIF
         ENDIF
      ENDDO

      ! Free pointer
      Spc => NULL()

      END SUBROUTINE DO_WASHOUT_AT_SFC
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: safety
!
! !DESCRIPTION: Subroutine SAFETY stops the run with debug output and an 
!  error message if negative species are found. 
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE SAFETY( I,         J,           L,        N,
     &                   A,         LS,          PDOWN,    QQ,
     &                   ALPHA,     ALPHA2,      RAINFRAC, WASHFRAC,
     &                   MASS_WASH, MASS_NOWASH, WETLOSS,  GAINED, 
     &                   LOST,      DSpc,        Spc,      RC        )

!
! !USES:
!
      USE CMN_SIZE_MOD
      USE ErrCode_Mod
!
! !INPUT PARAMETERS: 
!
      ! Arguments
      LOGICAL,          INTENT(IN)  :: LS            !
      INTEGER,          INTENT(IN)  :: I             !
      INTEGER,          INTENT(IN)  :: J             !
      INTEGER,          INTENT(IN)  :: L             !
      INTEGER,          INTENT(IN)  :: N             !
      CHARACTER(LEN=*), INTENT(IN)  :: A             !
      REAL(fp),         INTENT(IN)  :: PDOWN         !
      REAL(fp),         INTENT(IN)  :: QQ            !
      REAL(fp),         INTENT(IN)  :: ALPHA         !
      REAL(fp),         INTENT(IN)  :: ALPHA2        !
      REAL(fp),         INTENT(IN)  :: RAINFRAC      !
      REAL(fp),         INTENT(IN)  :: WASHFRAC      ! 
      REAL(fp),         INTENT(IN)  :: MASS_WASH     !
      REAL(fp),         INTENT(IN)  :: MASS_NOWASH   !
      REAL(fp),         INTENT(IN)  :: WETLOSS       !
      REAL(fp),         INTENT(IN)  :: GAINED        !
      REAL(fp),         INTENT(IN)  :: LOST          !
      REAL(fp),         INTENT(IN)  :: DSpc(LLPAR)   !
      REAL(fp),         INTENT(IN)  :: Spc(LLPAR)    !
!
! ! OUTPUT PARAMETERS:
!
      INTEGER,          INTENT(OUT) :: RC            ! Success or failure?
! 
! !REVISION HISTORY: 
!  18 Mar 2004 - R. Yantosca - Initial version
!  (1 ) Now made into a MODULE routine since we cannot call internal routines
!        from w/in a parallel loop.  Updated comments. (bmy, 3/18/04)
!  16 Sep 2010 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      ! Strings
      CHARACTER(LEN=255) :: ErrMsg, ThisLoc

      !=================================================================
      ! SAFETY begins here!
      !=================================================================
      
      ! Initialize
      RC      = GC_SUCCESS
      ErrMsg  = 'Error encountered in wet deposition!'
      ThisLoc = ' -> at SAFETY (in module GeosCore/wetscav_mod.F)'

      ! Print line
      WRITE( 6, '(a)' ) REPEAT( '=', 79 )

      ! Write error message and stop the run
      WRITE ( 6, 100 ) I, J, L, N, TRIM( A )
 100  FORMAT( 'WETDEP: ERROR at ', 3i4, ' for species ', i4, 
     &        ' in area ', a )

      PRINT*, 'LS          : ', LS
      PRINT*, 'PDOWN       : ', PDOWN
      PRINT*, 'QQ          : ', QQ
      PRINT*, 'ALPHA       : ', ALPHA
      PRINT*, 'ALPHA2      : ', ALPHA2
      PRINT*, 'RAINFRAC    : ', RAINFRAC
      PRINT*, 'WASHFRAC    : ', WASHFRAC
      PRINT*, 'MASS_WASH   : ', MASS_WASH
      PRINT*, 'MASS_NOWASH : ', MASS_NOWASH
      PRINT*, 'WETLOSS     : ', WETLOSS
      PRINT*, 'GAINED      : ', GAINED
      PRINT*, 'LOST        : ', LOST
      PRINT*, 'DSpc(NW,:)  : ', DSpc(:)
      PRINT*, 'Spc(I,J,:N) : ', Spc(:)

      ! Print line
      WRITE( 6, '(a)' ) REPEAT( '=', 79 )

      ! Return with error
      CALL GC_Error( ErrMsg, RC, ThisLoc )

      END SUBROUTINE SAFETY
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_vud
!
! !DESCRIPTION: Function GET\_VUD returns the vertical updraft velocity in
!  m/s at location I, J, L. 
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_VUD( State_Met, Input_Opt, I, J, L ) RESULT( VUD )
!
! !USES:
!
      USE DAO_MOD,            ONLY : IS_WATER
      USE State_Met_Mod,      ONLY : MetState
      USE Input_Opt_Mod,      ONLY : OptInput 
!
! !INPUT PARAMETERS: 
!
      TYPE(MetState), INTENT(IN) :: State_Met   ! Meteorology State object
      TYPE(OptInput), INTENT(IN) :: Input_Opt   ! Input Options object
      INTEGER,        INTENT(IN) :: I, J, L     ! Location 
!
! !RETURN VALUE:
!
      REAL(fp)                   :: VUD    ! Vertical updraft velocity in m/s. 
! 
! !REVISION HISTORY:
!  12 Feb 2015 - C. Keller   - Initial version
!  26 Feb 2015 - E. Lundgren - Remove dependency on pressure_mod since not used. 
!  27 Jun 2016 - C. Keller   - Use VUD_DEFINED instead of ESMF_ compiler switch.
!  13 Mar 2017 - C. Keller   - Bug fix: make sure VUD_DEFINED is updated after
!                              defining VUD.
!  01 Jun 2017 - C. Keller   - Set VUD_DEFINED only if UPDVVEL has valid values. 
!  05 Jun 2018 - C. Keller   - Add flag for online VUD calculation
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL :: VUD_DEFINED

      !=================================================================
      ! GET_VUD begins here!
      !=================================================================

      ! Init
      VUD_DEFINED = .FALSE.

      !=================================================================
      ! Use vertical updraft velocity from State_Met
      ! Convert hPa/s to m/s here
      !=================================================================
#if defined( MODEL_GEOS )
      IF ( Input_Opt%UseOnlineVUD ) THEN
         IF ( ASSOCIATED( State_Met%UPDVVEL ) ) THEN
            IF ( State_Met%DELP(I,J,L)    > TINY_FP .AND. 
     &           State_Met%UPDVVEL(I,J,L) > 0.0_fp       ) THEN
    
               ! Compute VUD 
               VUD = State_Met%UPDVVEL (I,J,L) 
     &             * State_Met%BXHEIGHT(I,J,L) 
     &             / State_Met%DELP    (I,J,L) 
   
               ! VUD is now defined
               VUD_DEFINED = .TRUE.
            ENDIF
         ENDIF
      ENDIF
#else
      IF ( ASSOCIATED( State_Met%UPDVVEL ) ) THEN
         IF ( State_Met%DELP(I,J,L)    > TINY_FP .AND. 
     &        State_Met%UPDVVEL(I,J,L) > 0.0_fp       ) THEN
 
            ! Compute VUD 
            VUD = State_Met%UPDVVEL (I,J,L) 
     &          * State_Met%BXHEIGHT(I,J,L) 
     &          / State_Met%DELP    (I,J,L) 

            ! VUD is now defined
            VUD_DEFINED = .TRUE.
          ENDIF
       ENDIF
#endif

      !=================================================================
      ! Traditional GEOS-Chem:
      ! Compute Vud -- 5 m/s over oceans, 10 m/s over land (or ice?)
      ! Assume Vud is the same at all altitudes; the array can be 2-D
      !=================================================================
      IF ( .NOT. VUD_DEFINED ) THEN
         IF ( IS_WATER( I, J, State_Met ) ) THEN
            VUD = 5e+0_fp
         ELSE
            VUD = 10e+0_fp
         ENDIF
      ENDIF

      END FUNCTION GET_VUD
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_f
!
! !DESCRIPTION: Function GET\_F returns the scavenged fraction at location
! I, J, L and for the given rate constant K.
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_F( State_Met, Input_Opt, I, J, L, K ) RESULT( F )
!
! !USES:
!
      USE DAO_MOD,            ONLY : IS_WATER
      USE State_Met_Mod,      ONLY : MetState
      USE Input_Opt_Mod,      ONLY : OptInput 
!
! !INPUT PARAMETERS: 
!
      TYPE(MetState), INTENT(IN) :: State_Met  ! Meteorology State object
      TYPE(OptInput), INTENT(IN) :: Input_Opt   ! Input Options object
      INTEGER,        INTENT(IN) :: I, J, L    ! Lon, lat, level indices
      REAL(fp),       INTENT(IN) :: K          ! Rate constant
!
! !RETURN VALUE:
!
      REAL(fp)                   :: F          ! Fraction of species scavenged
                                               !  out of the updraft [1]
! 
! !REVISION HISTORY:
!  12 Feb 2015 - C. Keller   - Initial version 
!  26 Feb 2015 - E. Lundgren - Remove dependency on pressure_mod since not used.
!  25 Jun 2015 - M. Sulprizio- Now pass FF as an optional argument for aerosols
!  25 Sep 2015 - R. Yantosca - Remove FF; we scale Kc in F_AEROSOL now.
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      REAL(fp)    :: TMP, VUD

      !=================================================================
      ! GET_F begins here!
      !=================================================================

      ! Distance between grid box centers [m]
      TMP = 0.5_fp * ( State_Met%BXHEIGHT(I,J,L-1) + 
     &                 State_Met%BXHEIGHT(I,J,L) )
               
      ! Vertical updraft velocity [m/s]
      Vud = GET_VUD( State_Met, Input_Opt, I, J, L )

      ! Compute F (avoid div-by-zero errors)
      IF ( Vud > TINY_FP ) THEN
         F = 1.0_fp - EXP( -K * TMP / Vud )
      ELSE
         F = 0.0_fp
      ENDIF

      END FUNCTION GET_F
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_wetscav
!
! !DESCRIPTION: Subroutine INIT\_WETSCAV initializes updraft velocity, cloud 
!  liquid water content, cloud ice content, and mixing ratio of water fields, 
!  which are used in the wet scavenging routines.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE INIT_WETSCAV( am_I_Root, Input_Opt, 
     &                         State_Chm, State_Diag, RC )
!
! !USES:
!
      USE CMN_Size_Mod
      USE ErrCode_Mod
      USE Input_Opt_Mod,  ONLY : OptInput
      USE Species_Mod,    ONLY : Species
      USE State_Chm_Mod,  ONLY : ChmState
      USE State_Chm_Mod,  ONLY : Ind_
      USE State_Diag_Mod, ONLY : DgnState
!
! !INPUT PARAMETERS:
!
      LOGICAL,        INTENT(IN)  :: am_I_Root   ! Are we on the root CPU?
      TYPE(OptInput), INTENT(IN)  :: Input_Opt   ! Input Options object
      TYPE(ChmState), INTENT(IN)  :: State_Chm   ! Chemistry State object
      TYPE(DgnState), INTENT(IN)  :: State_Diag  ! Diagnostics State object
!
! !OUTPUT PARAMETERS:
!
      INTEGER,        INTENT(OUT) :: RC          ! Success or failure?
! 
! !REVISION HISTORY: 
!  03 Sep 2015 - R. Yantosca - Split off from the old INIT_WETSCAV.  We now
!                              just allocate arrays here.  Arrays are assigned
!                              values from State_Met in SETUP_WETSCAV.
!  03 Sep 2015 - R. Yantosca - Also initialize IDWETD here; this makes routine
!                              WETDEPID obsolete, so we can remove it.
!  22 Sep 2015 - R. Yantosca - Bug fix: only define the IDWETD array if
!                              wetdep or drydep is turned on.  Also only
!                              print info if wetdep is turned on.
!  29 Apr 2016 - R. Yantosca - Don't initialize pointers in declaration stmts
!  20 Jun 2016 - R. Yantosca - Now define species IDs in the init phase
!  09 Nov 2017 - R. Yantosca - Now accept State_Diag as an argument
!  09 Nov 2017 - R. Yantosca - Return error condition to calling program
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      ! Scalars
      INTEGER                :: N, N_WD
      LOGICAL, SAVE          :: FIRST = .TRUE.      
            
      ! Strings
      CHARACTER(LEN=9)       :: K0,  CR,  pKA
      CHARACTER(LEN=80)      :: LINE

      ! Pointers
      TYPE(Species), POINTER :: SpcInfo

      !=================================================================
      ! INIT_WETSCAV begins here!
      !=================================================================

      ! Assume success
      RC       = GC_SUCCESS

      ! Exit if we have already executed this routine
      IF ( .not. FIRST ) RETURN

      ! Initialize
      SpcInfo  => NULL()
      
      ! Define species ID flags
      id_NK1   = Ind_('NK1'  )
      id_SF1   = Ind_('SF1'  )
      id_DUST1 = Ind_('DUST1')
      id_SO2   = Ind_('SO2'  )
      id_SO4   = Ind_('SO4'  )
      id_H2O2  = Ind_('H2O2' )

      !-----------------------------------------------------------------
      ! Only allocate arrays if wetdep or convection is turned on
      !-----------------------------------------------------------------
      IF ( Input_Opt%LWETD .or. Input_Opt%LCONV ) THEN

         ! Allocate PDOWN on first call
         ALLOCATE( PDOWN( LLPAR, IIPAR, JJPAR ), STAT=RC )
         CALL GC_CheckVar( 'wetscav_mod.F:PDOWN', 0, RC )
         IF ( RC /= GC_SUCCESS ) RETURN
         PDOWN = 0e+0_fp
      
         ! Allocate QQ on first call
         ALLOCATE( QQ( LLPAR, IIPAR, JJPAR ), STAT=RC )
         CALL GC_CheckVar( 'wetscav_mod.F:QQ', 0, RC )
         IF ( RC /= GC_SUCCESS ) RETURN
         QQ = 0e+0_fp

         ! Allocate REEVAP on first call
         ALLOCATE( REEVAP( LLPAR, IIPAR, JJPAR ), STAT=RC )
         CALL GC_CheckVar( 'wetscav_mod.F:REEVAP', 0, RC )
         IF ( RC /= GC_SUCCESS ) RETURN
         REEVAP = 0d0

         ! Allocate C_H2O on first call
         ALLOCATE( C_H2O( IIPAR, JJPAR, LLPAR ), STAT=RC )
         CALL GC_CheckVar( 'wetscav_mod.F:C_H2O', 0, RC )
         IF ( RC /= GC_SUCCESS ) RETURN
         C_H2O = 0.0_fp
            
         ! Allocate CLDLIQ on first call
         ALLOCATE( CLDLIQ( IIPAR, JJPAR, LLPAR ), STAT=RC )
         CALL GC_CheckVar( 'wetscav_mod.F:CLDLIQ', 0, RC )
         IF ( RC /= GC_SUCCESS ) RETURN
         CLDLIQ = 0.0_fp
            
         ! Allocate CLDICE on first call
         ALLOCATE( CLDICE( IIPAR, JJPAR, LLPAR ), STAT=RC )
         CALL GC_CheckVar( 'wetscav_mod.F:CLDICE', 0, RC )
         IF ( RC /= GC_SUCCESS ) RETURN
         CLDICE = 0.0_fp

      ENDIF

      !=================================================================
      ! Print information about wet-depositing species
      !=================================================================
      IF ( am_I_Root .and. State_Chm%nWetDep > 0 ) THEN

         ! Title
         LINE = 'INIT_WETSCAV: List of soluble species: '
         WRITE( 6, '(/,a,/)' ) TRIM( LINE )

         ! 1st line
         LINE = '  #             Name  Species Mol Wt ' // 
     &          ' Henry K0  Henry CR  Henry pKa'
         WRITE( 6, '(a)'     ) TRIM( LINE )

         ! 2nd line
         LINE = 'ID    g/mol    M/atm      K         1'
         WRITE( 6, '(24x,a)' ) TRIM( LINE )

         ! Separator line
         WRITE( 6, '(a)'     ) REPEAT( '-', 70 )

         ! Loop over all wet-depositing species
         DO N_WD = 1, State_Chm%nWetDep

            ! Get the corresponding species ID
            N = State_Chm%Map_WetDep(N_WD)

            ! Get physical parameters from the species database object
            SpcInfo => State_Chm%SpcData(N)%Info 

            ! Convert Henry's law K0 parameter to string
            IF ( SpcInfo%Henry_K0 > 0.0_f8 ) THEN
               WRITE( K0, '(es9.2)' ) SpcInfo%Henry_K0
            ELSE
               K0 = '    -    '
            ENDIF

            ! Convert Henry's law CR parameter to string
            IF ( SpcInfo%Henry_K0 > 0.0_f8 ) THEN
               WRITE( CR, '(es9.2)' ) SpcInfo%Henry_CR
            ELSE
               CR = '    -    '
            ENDIF

            ! Convert Henry's law pKa parameter to string
            IF ( SpcInfo%Henry_pKa > 0.0_f8 ) THEN
               WRITE( pKa, '(es9.2)' ) SpcInfo%Henry_pKa
            ELSE
               pKa = '    -    '
            ENDIF

            ! Write info to stdout
            WRITE( 6, 100 ) N_WD,           
     &                      TRIM( SpcInfo%Name ),  
     &                      N,    
     &                      SpcInfo%EmMW_g, 
     &                      K0, CR, pKA

 100        FORMAT( i3,3x,a14,3x,i3,3x,f6.1,3(1x,a9) )

            ! Free pointer
            SpcInfo => NULL()

         ENDDO
      ENDIF


      ! Reset flag
      FIRST = .FALSE. 

      END SUBROUTINE INIT_WETSCAV
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: setup_wetscav
!
! !DESCRIPTION: Subroutine SETUP\_WETSCAV initializes updraft velocity, cloud 
!  liquid water content, cloud ice content, and mixing ratio of water fields, 
!  which are used in the wet scavenging routines.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE SETUP_WETSCAV( am_I_Root, Input_Opt, 
     &                          State_Met, State_Chm, RC )
!
! !USES:
!
      USE CMN_Size_Mod
      USE ErrCode_Mod
      USE Error_Mod,          ONLY : Alloc_Err
      USE Input_Opt_Mod,      ONLY : OptInput
      USE PhysConstants,      ONLY : AIRMW
      USE State_Chm_Mod,      ONLY : ChmState
      USE State_Met_Mod,      ONLY : MetState
      USE UnitConv_Mod
!
! !INPUT PARAMETERS:
!
      LOGICAL,        INTENT(IN)    :: am_I_Root   ! Are we on the root CPU?
      TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
      TYPE(MetState), INTENT(IN)    :: State_Met   ! Meteorology State object
!
! !INPUT/OUTPUT PARAMETERS:
!
      TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
!
! !OUTPUT PARAMETERS:
!
      INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
! 
! !REVISION HISTORY: 
!  23 Feb 2000 - R. Yantosca - Initial version
!  (1 ) References "e_ice.f" -- routine to compute Eice(T).
!  (2 ) Vud, CLDLIQ, CLDICE, C_H2O are all independent of tracer, so we
!        can compute them once per timestep, before calling the cloud 
!        convection and wet deposition routines.
!  (3 ) Set C_H2O = 0 below -120 Celsius.  E_ICE(T) has a lower limit of
!        -120 Celsius, so temperatures lower than this will cause a stop
!        with an error message. (bmy, 6/15/00)
!  (4 ) Replace {IJL}GLOB with IIPAR,JJPAR,LLPAR.  Also rename PW to P.
!        Remove IREF, JREF, these are obsolete.  Now reference IS_WATER
!        from "dao_mod.f" to determine water boxes. 
!  (5 ) Removed obsolete code from 9/01.  Updated comments and made
!        cosmetic changes. (bmy, 10/24/01)
!  (6 ) Now use routine GET_PCENTER from "pressure_mod.f" to compute the
!        pressure at the midpoint of grid box (I,J,L).  Also removed P and
!        SIG from the argument list (dsa, bdf, bmy, 8/20/02)
!  (7 ) Now reference T from "dao_mod.f".  Updated comments.  Now allocate
!        Vud, C_H2O, CLDLIQ and CLDICE here on the first call.  Now references
!        ALLOC_ERR from "error_mod.f".  Now set H2O2s and SO2s to the initial
!        values from for the first call to COMPUTE_F .  Now call WETDEPID
!        on the first call to initialize the wetdep index array. (bmy, 1/27/03)
!  (8 ) Now references STT from "tracer_mod.f".  Also now we call WETDEPID
!        from "input_mod.f" (bmy, 7/20/04)
!  (9 ) Now references new function E_ICE, which is an analytic function of 
!        Kelvin temperature instead of Celsius. (bmy, 3/7/05)
!  16 Sep 2010 - R. Yantosca - Added ProTeX headers
!  18 Oct 2012 - R. Yantosca - Removed DEVEL from #ifdef statement, now use
!                              EXTERNAL_GRID || EXTERNAL_FORCING
!  09 Nov 2012 - M. Payer    - Replaced all met field arrays with State_Met
!                              derived type object
!  26 Feb 2015 - E. Lundgren - Replace GET_PCENTER with State_Met%PMID_DRY.
!                              Remove dependency on pressure_mod.
!  28 Apr 2015 - E. Lundgren - Change PMID_DRY to PMID for Dalton's Law.
!  13 Aug 2015 - E. Lundgren - Incoming tracer units are now kg/kg dry air
!  03 Sep 2015 - R. Yantosca - Now allocate IDWETD
!  03 Sep 2015 - R. Yantosca - Now move array allocation into INIT_WETSCAV
!  22 Sep 2015 - R. Yantosca - Bug fix: now initialize H2O2s and SO2s on
!                              the first call.  This was done in INIT_WETSCAV,
!                              but that is now called before the restart file
!                              is read from disk.
!  23 Sep 2015 - R. Yantosca - Now accept State_Chm as an argument
!  29 Apr 2016 - R. Yantosca - Don't initialize pointers in declaration stmts
!  23 Jun 2016 - R. Yantosca - Remove reference to tracerid_mod.F
!  30 Jun 2016 - R. Yantosca - Remove instances of STT and State_Chm%TRACERS
!  04 Jul 2016 - M. Yannetti - Replace TCVV with species db MW and phys constant
!  05 Aug 2016 - R. Yantosca - Remove temporary tracer-removal code
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      ! Scalars
      LOGICAL, SAVE     :: FIRST = .TRUE.
      INTEGER           :: I, J, L

      ! Pointers
      REAL(fp), POINTER :: PL
      REAL(fp), POINTER :: TK

      !=====================================================================
      ! Compute Vud, CLDLIQ, CLDICE, C_H2O, following Jacob et al, 2000.
      !=====================================================================

      ! Initialize
      PL => NULL()
      TK => NULL()
      RC = GC_SUCCESS

      ! Only do computation if wetdep or convection is turned on
      IF ( Input_Opt%LWETD .or. Input_Opt%LCONV ) THEN

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, TK, PL )
!$OMP+SCHEDULE( DYNAMIC )
         DO L = 1, LLPAR
         DO J = 1, JJPAR
         DO I = 1, IIPAR
         
            ! Point to Temp [K] and Pressure [hPa]
            TK => State_Met%T(I,J,L)
            PL => State_Met%PMID(I,J,L)
         
            !--------------------------------------------------------------
            ! CLDLIQ, the cloud liquid water content [cm3 H2O/cm3 air], 
            ! is a function of the local Kelvin temperature:
            ! Tunable parameter and use 1e-6 here (qq,10/14/2011) 
            !    CLDLIQ = 1e-6                    [     T >= 268 K    ]
            !    CLDLIQ = 1e-6 * ((T - 248) / 20) [ 248 K < T < 268 K ]
            !    CLDLIQ = 0                       [     T <= 248 K    ]
            !
            !--------------------------------------------------------------
            IF ( TK >= 268.0_fp ) THEN
               CLDLIQ(I,J,L) = 1e-6_fp
         
            ELSE IF ( TK > 248.0_fp .and. TK < 268.0_fp ) THEN
               CLDLIQ(I,J,L) = 1e-6_fp * ((TK - 248.0_fp) / 20.0_fp )
         
            ELSE
               CLDLIQ(I,J,L) = 0.0_fp
               
            ENDIF
              
            !--------------------------------------------------------------
            ! CLDICE, the cloud ice content [cm3 ice/cm3 air] is given by:
            !
            !    CLDICE = 1.0e-6 - CLDLIQ
            !--------------------------------------------------------------
            CLDICE(I,J,L) = 1e-6_fp - CLDLIQ(I,J,L)
         
            !--------------------------------------------------------------
            ! C_H2O is given by Dalton's Law as:
            !
            !       C_H2O = Eice( Tk(I,J,L) ) / P(I,J,L)
            !
            ! where P(L) = pressure in grid box (I,J,L)
            !
            ! and   Tk(I,J,L) is the Kelvin temp. of grid box (I,J,L).
            !
            ! and   Eice( Tk(I,J,L) ) is the saturation vapor pressure 
            !       of ice [hPa] at temperature Tk(I,J,L) -- computed in 
            !       routine E_ICE above.
            !--------------------------------------------------------------
            IF ( PL <= TINY_FP ) THEN
               C_H2O(I,J,L) = 0.0_fp 
            ELSE
               C_H2O(I,J,L) = E_ICE( TK ) / PL
            ENDIF         

            ! Free pointers
            PL => NULL()
            TK => NULL()

         ENDDO
         ENDDO
         ENDDO
!$OMP END PARALLEL DO
      ENDIF

      !====================================================================
      ! We need to initialize the H2O2s and SO2s arrays to the values
      ! of State_Chm%Species for H2O2 and SO2.  This only needs to be
      ! done the first time this routine is called (which happens after
      ! the restart file is read from disk).
      !
      ! Prior to GEOS-Chem v11-01, this was done in routine INIT_WETSCAV, 
      ! which also allocated the arrays.  In GEOS-Chem v11-01 and higher
      ! versions, we do array allocations only in INIT_WETSCAV, which is
      ! called in the init phase (before restart file I/O).  We now
      ! call SETUP_WETSCAV to compute the CLDLIQ, CLDICE, and CH2O 
      ! arrays on each dynamic timestep.  (bmy, 9/22/15)
      !====================================================================
      IF ( FIRST ) THEN 

         ! Set H2O2s to the initial H2O2 from the species array, so that we will
         ! have nonzero values for the first call to COMPUTE_F (bmy, 1/14/03)
         ! While State_Chm%Species are now in units of [kg/kg dry air] for
         ! call to SETUP_WETSCAV, store H2O2s in legacy units [v/v dry air] 
         ! for now for use in sulfate_mod and WASHOUT (ewl, 10/15/15)
         IF ( id_H2O2 > 0 ) THEN
            State_Chm%H2O2AfterChem = State_Chm%Species(:,:,:,id_H2O2)
     &            * ( AIRMW / State_Chm%SpcData(id_H2O2)%Info%emMW_g )
         ENDIF

         ! Set SO2s to the initial SO2 from the species array, so that we will
         ! have nonzero values for the first call to COMPUTE_F (bmy, 1/14/03)
         ! While State_Chm%Species are now in units of [kg/kg dry air] for 
         ! call to SETUP_WETSCAV, store SO2s in units [v/v dry air] for now 
         ! for use in sulfate_mod and WASHOUT (ewl, 10/15/15)
         IF ( id_SO2 > 0 ) THEN
            State_Chm%SO2AfterChem = State_Chm%Species(:,:,:,id_SO2)
     &           * ( AIRMW / State_Chm%SpcData(id_SO2)%Info%emMW_g )
         ENDIF

         ! Reset first-time flag
         FIRST = .FALSE.
      ENDIF

      END SUBROUTINE SETUP_WETSCAV
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: cleanup_wetscav
!
! !DESCRIPTION: Subroutine CLEANUP\_WETSCAV deallocates all module arrays.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE CLEANUP_WETSCAV( am_I_Root, RC )
!
! !USES:
!
      USE ErrCode_Mod
!
! !INPUT_PARAMETERS:
! 
      LOGICAL, INTENT(IN)  :: am_I_Root   ! Are we on the root CPU?
!
! !OUTPUT PARAMETERS:
!
      INTEGER, INTENT(OUT) :: RC          ! Success or failure?
!
! !REVISION HISTORY: 
!  23 Feb 2000 - R. Yantosca - Initial version
!  16 Sep 2010 - R. Yantosca - Added ProTeX headers
!  03 Sep 2015 - R. Yantosca - Now deallocate IDWETD
!  17 Mar 2017 - R. Yantosca - Now remove IDWETD, it's obsolete
!  09 Nov 2017 - R. Yantosca - Return error condition to calling program
!EOP
!------------------------------------------------------------------------------
!BOC
!
      !=================================================================
      ! CLEANUP_WETSCAV begins here!
      !=================================================================
      IF ( ALLOCATED( C_H2O  ) ) DEALLOCATE( C_H2O,  STAT=RC )
      CALL GC_CheckVar( 'wetscav_mod.F:C_H2O', 2, RC )
      IF ( RC /= GC_SUCCESS ) RETURN

      IF ( ALLOCATED( CLDLIQ ) ) DEALLOCATE( CLDLIQ, STAT=RC )
      CALL GC_CheckVar( 'wetscav_mod.F:CLDLIQ', 2, RC )
      IF ( RC /= GC_SUCCESS ) RETURN

      IF ( ALLOCATED( CLDICE ) ) DEALLOCATE( CLDICE, STAT=RC )
      CALL GC_CheckVar( 'wetscav_mod.F:CLDICE', 2, RC )
      IF ( RC /= GC_SUCCESS ) RETURN

      IF ( ALLOCATED( PDOWN  ) ) DEALLOCATE( PDOWN,  STAT=RC )
      CALL GC_CheckVar( 'wetscav_mod.F:PDOWN', 2, RC )
      IF ( RC /= GC_SUCCESS ) RETURN

      IF ( ALLOCATED( QQ     ) ) DEALLOCATE( QQ,     STAT=RC )
      CALL GC_CheckVar( 'wetscav_mod.F:QQ', 2, RC )
      IF ( RC /= GC_SUCCESS ) RETURN

      IF ( ALLOCATED( REEVAP ) ) DEALLOCATE( REEVAP, STAT=RC )
      CALL GC_CheckVar( 'wetscav_mod.F:REEVAP', 2, RC )
      IF ( RC /= GC_SUCCESS ) RETURN

      END SUBROUTINE CLEANUP_WETSCAV
!EOC
      END MODULE WETSCAV_MOD
