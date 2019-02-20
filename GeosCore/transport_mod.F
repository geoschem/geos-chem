!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: transport_mod.F
!
! !DESCRIPTION: Module TRANSPORT\_MOD is used to call the proper version of 
!  the TPCORE advection scheme for different met field datasets and their
!  nested or global grids.
!\\
!\\
! !INTERFACE: 
!
      MODULE TRANSPORT_MOD
! 
! !USES:
!
      USE PRECISION_MOD      ! For GEOS-Chem Precision (fp)
      USE PRESSURE_MOD

      IMPLICIT NONE
      PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
      PUBLIC  :: CLEANUP_TRANSPORT
      PUBLIC  :: DO_TRANSPORT
      PUBLIC  :: INIT_TRANSPORT
      PUBLIC  :: INIT_WINDOW_TRANSPORT
!
! !PRIVATE MEMBER FUNCTIONS:
!
      PRIVATE :: DO_GLOBAL_ADV      
      PRIVATE :: DO_WINDOW_TRANSPORT
!
! !REVISION HISTORY:
!  10 Mar 2003 - Y. Wang, R. Yantosca - Initial version
!  (1 ) Now can select transport scheme for GEOS-3 winds.  Added code for PJC 
!        pressure fixer. (bdf, bmy, 5/8/03)
!  (2 ) Now delete DSIG array, it's obsolete.  Also added new PRIVATE function 
!        GET_AIR_MASS to compute air masses from the input/output pressures
!        from the new GEOS-4/fvDAS TPCORE. (bmy, 6/24/03)
!  (3 ) Now references DEBUG_MSG from "error_mod.f". (bmy, 8/7/03)
!  (4 ) Bug fix in DO_GLOBAL_TRANSPORT (bmy, 10/21/03)
!  (5 ) IORD, JORD, KORD are now module variables.  Now references 
!        "logical_mod.f" and "tracer_mod.f" (bmy, 7/20/04)
!  (6 ) Add mass-flux diagnostics to TPCORE_FVDAS (bdf, bmy, 9/28/04)
!  (7 ) Now references "diag_mod.f" (bmy, 9/28/04)
!  (8 ) Now modified for GEOS-5 and GCAP met fields (swu, bmy, 5/25/05)
!  (9 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (10) Now flip arrays in call to TPCORE_FVDAS (bmy, 6/16/06)
!  (11) Added modifications for SUN compiler (bmy, 7/12/06)
!  (12) Bug fixes in DO_GLOBAL_TRANSPORT (bmy, 11/29/06)
!  (13) Split off GCAP, GEOS-3, GEOS-4/GEOS-5 specific calling sequences
!        into separate subroutines.  Also removed some obsolete module
!        variables. (bmy, 10/30/07)
!  (14) Modifications for GEOS-5 nested grid (yxw, dan, bmy, 11/6/08)
!  (15) Bug fix in mass balance in GCAP_GLOBAL_ADV and GEOS4_GEOS5_GLOBAL_ADV.
!        (ccc, 2/17/09)
!  26 Feb 2010 - R. Yantosca - Removed references to obsolete LEMBED switch
!  26 Feb 2010 - R. Yantosca - Added ProTex Headers
!  08 Mar 2010 - C. Carouge  - Modify call to tpcore_fvdas. We do not re-order 
!                              mass fluxes diagnostics anymore.
!  28 Feb 2012 - R. Yantosca - Removed support for GEOS-3
!  01 Mar 2012 - R. Yantosca - Now reference new grid_mod.F90
!  21 Jun 2012 - R. Yantosca - Now use pointers to flip arrays in vertical
!  21 Jun 2012 - R. Yantosca - Comment out GEOS-3 window subroutine
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
!  26 Sep 2013 - R. Yantosca - Renamed "GEOS57" to "GEOSFP" in routine names
!  24 Nov 2014 - M. Yannetti - Added PRECISION_MOD
!  12 Feb 2015 - E. Lundgren - Added functionality for writing diags to netcdf
!  15 Apr 2015 - R. Yantosca - Add TARGET spec to A_M2 array
!  24 Jun 2015 - E. Lundgren - Now use kg/kg total air as tracer units
!                              (previously v/v) for advection
!  24 Jun 2015 - E. Lundgren - Remove post-advection residual mass correction
!                              and GET_AIR_MASS function since no longer needed
!                              when using kg/kg total air in advection
!  13 Aug 2015 - E. Lundgren - Tracer units are now input as kg/kg dry air 
!                              (previously v/v)
!  19 Jan 2016 - E. Lundgren - Move DiagnUpdate_Transport to diagnostics_mod
!                              and consolidate netcdf and bpch diagnostic code
!  23 Jun 2016 - R. Yantosca - Remove references to APM code; it is no longer
!                              compatible with the FlexChem implementation
!  29 Nov 2016 - R. Yantosca - grid_mod.F90 is now gc_grid_mod.F90
!  24 Aug 2017 - M. Sulprizio- Remove support for GCAP, GEOS-4, GEOS-5 and MERRA
!                              and rename routine GEOS4_GEOS5_GLOBAL_ADV to
!                              DO_GLOBAL_ADV
!EOP
!------------------------------------------------------------------------------
!BOC
      !=================================================================
      ! MODULE VARIABLES:
      !
      ! (1 ) Ap     (REAL(fp) ) : Vertical coordinate array for TPCORE
      ! (2 ) A_M2   (REAL(fp) ) : Grid box surface areas [m2]
      ! (3 ) Bp     (REAL(fp) ) : Vertical coordinate array for TPCORE
      ! (7 ) JLAST  (INTEGER)   : For fvDAS TPCORE
      ! (8 ) MG     (INTEGER)   : For fvDAS TPCORE
      ! (9 ) NG     (INTEGER)   : For fvDAS TPCORE
      ! (10) N_ADJ  (INTEGER)   : For fvDAS TPCORE
      !=================================================================
      INTEGER                       :: JFIRST 
      INTEGER                       :: JLAST, NG,   MG,   N_ADJ
      REAL(fp), ALLOCATABLE         :: Ap(:)
      REAL(fp), ALLOCATABLE         :: Bp(:)
      REAL(fp), ALLOCATABLE, TARGET :: A_M2(:)

      CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: do_transport
!
! !DESCRIPTION: Subroutine DO\_TRANSPORT is the driver routine for the proper 
!  TPCORE program for GEOS-3, GEOS-4/GEOS-5, or window simulations.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE DO_TRANSPORT( am_I_Root, Input_Opt,  State_Met,
     &                         State_Chm, State_Diag, RC         )
!
! !USES:
!
      USE CMN_SIZE_MOD
      USE Diagnostics_Mod, ONLY : Compute_Column_Mass
      USE Diagnostics_Mod, ONLY : Compute_Budget_Diagnostics
      USE ErrCode_Mod
      USE GC_GRID_MOD,    ONLY : ITS_A_NESTED_GRID
      USE Input_Opt_Mod,  ONLY : OptInput
      USE State_Chm_Mod,  ONLY : ChmState
      USE State_Diag_Mod, ONLY : DgnState
      USE State_Met_Mod,  ONLY : MetState
      USE TPCORE_BC_MOD,  ONLY : INIT_TPCORE_BC
      USE TIME_MOD,       ONLY : GET_TS_DYN
#if defined( USE_TEND )
      USE TENDENCIES_MOD
#endif
!
! !INPUT PARAMETERS:
!
      LOGICAL,        INTENT(IN)    :: am_I_Root   ! Are we on the root CPU?
      TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
!
! !INPUT/OUTPUT PARAMETERS:
!
      TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
      TYPE(DgnState), INTENT(INOUT) :: State_Diag  ! Diagnostics State object
      TYPE(MetState), INTENT(INOUT) :: State_Met   ! Meteorology State object
!
! !OUTPUT PARAMETERS:
!
      INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
! 
! !REVISION HISTORY: 
!  10 Mar 2003 - R. Yantosca - Initial version
!  (1 ) Removed IORD, JORD, KORD from the arg list.  Also now removed
!        reference to CMN, it's not needed. (bmy, 7/20/04)
!  (2 ) Now call separate routines for different met fields. (bmy, 10/30/07)
!  (3 ) Now references subroutine INIT_TPCORE_BC from tpcore_bc_mod.f and
!        DO_GEOS5_FVDAS_WINDOW_TRANSPORT from 
!        "tpcore_geos5_fvdas_window_mod.f90". (yxw, dan, bmy, 11/6/08)
!  26 Feb 2010 - R. Yantosca - Removed references to obsolete LEMBED switch
!  26 Feb 2010 - R. Yantosca - Added ProTeX headers
!  06 Oct 2010 - R. Yantosca - Treat MERRA in the same way as GEOS-5.
!  28 Feb 2012 - R. Yantosca - Treat GEOS-5.7 in the same way as MERRA
!  28 Feb 2012 - R. Yantosca - Removed support for GEOS-3
!  25 Mar 2013 - R. Yantosca - Now accept am_I_Root, Input_Opt, State_Chm, RC
!  31 Mar 2015 - E. Lundgren - Make State_Met INTENT(INOUT) for AIRQNT calls
!  19 Jul 2016 - R. Yantosca - Now bracket DO_TEND calls with #ifdef USE_TEND
!  03 Nov 2017 - R. Yantosca - Now accept State_Diag as an argument
!  09 Nov 2017 - R. Yantosca - Return error condition to calling program
!  27 Aug 2018 - E. Lundgren - Implement budget diagnostics
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      ! Scalars
      LOGICAL, SAVE      :: FIRST = .TRUE.
      CHARACTER(LEN=255) :: ErrMsg, ThisLoc
      REAL(fp)           :: DT_Dyn

      !=================================================================
      ! DO_TRANSPORT begins here!
      !=================================================================

      ! Initialize
      RC      = GC_SUCCESS
      ErrMsg  = ''
      ThisLoc = 
     &   ' -> at Do_Transport  (in module GeosCore/transport_mod.F)'

      !----------------------------------------------------------
      ! Transport (advection) budget diagnostics - Part 1 of 2
      !----------------------------------------------------------
      IF ( State_Diag%Archive_BudgetTransport ) THEN
         ! Get initial column masses
         CALL Compute_Column_Mass( am_I_Root, 
     &                          Input_Opt, State_Met, State_Chm, 
     &                          State_Chm%Map_Advect,        
     &                          State_Diag%Archive_BudgetTransportFull, 
     &                          State_Diag%Archive_BudgetTransportTrop, 
     &                          State_Diag%Archive_BudgetTransportPBL,  
     &                          State_Diag%BudgetMass1, 
     &                          RC )
         IF ( RC /= GC_SUCCESS ) THEN
            ErrMsg = 'Transport budget diagnostics error 1'
            CALL GC_Error( ErrMsg, RC, ThisLoc )
            RETURN
         ENDIF
      ENDIF

      ! First-time initialization
      IF ( FIRST ) THEN 

#if defined( GRID05x0625 ) || defined( GRID025x03125 )

         ! All nested grid simulations
         CALL INIT_WINDOW_TRANSPORT( am_I_Root,  Input_Opt, 
     &                               State_Diag, RC         )

         ! Trap potential errors
         IF ( RC /= GC_SUCCESS ) THEN
            ErrMsg = 'Error encountered in "Init_Window_Transport"!'
            CALL GC_Error( ErrMsg, RC, ThisLoc )
            RETURN
         ENDIF

#else

         ! All global simulations
         CALL INIT_TRANSPORT( am_I_Root, Input_Opt, State_Diag, RC )

         ! Trap potential errors
         IF ( RC /= GC_SUCCESS ) THEN
            ErrMsg = 'Error encountered in "Init_Transport"!'
            CALL GC_Error( ErrMsg, RC, ThisLoc )
            RETURN
         ENDIF
#endif

         FIRST = .FALSE.

      ENDIF

#if defined( USE_TEND )
      !=================================================================
      ! Archive tracer concentrations for tendencies (ckeller,7/15/2015)
      !=================================================================
      CALL TEND_STAGE1( am_I_Root, Input_Opt, State_Met, 
     &                  State_Chm, 'ADV',     RC         )

      ! Trap potential errors
      IF ( RC /= GC_SUCCESS ) THEN
         ErrMsg = 'Error encountered in "Tend_Stage1"!'
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
#endif

      !=================================================================
      ! Choose the proper version of TPCORE for the nested-grid window 
      ! region (usually 1x1 grids) or for the entire globe
      !=================================================================
      IF ( ITS_A_NESTED_GRID() ) THEN

#if defined( GRID05x0625 ) || defined( GRID025x03125 )

         ! Nested-grid simulation with GEOS-FP/MERRA2 met
         CALL DO_WINDOW_TRANSPORT( am_I_Root, Input_Opt,  State_Met,
     &                             State_Chm, State_Diag, RC         )

         ! Trap potential errors
         IF ( RC /= GC_SUCCESS ) THEN
            ErrMsg = 'Error encountered in "Do_Window_Transport"!'
            CALL GC_Error( ErrMsg, RC, ThisLoc )
            RETURN
         ENDIF

#endif

      !=================================================================
      ! Choose the proper version of TPCORE for global simulations
      !=================================================================
      ELSE

         ! Call TPCORE w/ proper settings for the GEOS-FP/MERRA2 met fields
         CALL DO_GLOBAL_ADV( am_I_Root, Input_Opt,  State_Met,
     &                       State_Chm, State_Diag, RC         )

         ! Trap potential errors
         IF ( RC /= GC_SUCCESS ) THEN
            ErrMsg = 'Error encountered in "Do_Global_Adv"!'
            CALL GC_Error( ErrMsg, RC, ThisLoc )
            RETURN
         ENDIF

      ENDIF

      ! Dynamic timestep [s]
      DT_Dyn = GET_TS_DYN()

#if defined( USE_TEND )
      !=================================================================
      ! Calculate tendencies and write to diagnostics 
      ! (ckeller,7/15/2015)
      !=================================================================
      CALL TEND_STAGE2( am_I_Root, Input_Opt, State_Met, 
     &                  State_Chm, 'ADV', DT_Dyn, RC )

      ! Trap potential errors
      IF ( RC /= GC_SUCCESS ) THEN
         ErrMsg = 'Error encountered in "Tend_Stage2"!'
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
#endif

      !----------------------------------------------------------
      ! Transport (advection) budget diagnostics - Part 2 of 2
      !----------------------------------------------------------
      IF ( State_Diag%Archive_BudgetTransport ) THEN
         ! Get final column masses and compute diagnostics
         CALL Compute_Column_Mass( am_I_Root, 
     &                          Input_Opt, State_Met, State_Chm, 
     &                          State_Chm%Map_Advect,        
     &                          State_Diag%Archive_BudgetTransportFull, 
     &                          State_Diag%Archive_BudgetTransportTrop, 
     &                          State_Diag%Archive_BudgetTransportPBL,  
     &                          State_Diag%BudgetMass2, 
     &                          RC ) 
         CALL Compute_Budget_Diagnostics( am_I_Root, 
     &                          State_Chm%Map_Advect,
     &                          DT_Dyn,                                
     &                          State_Diag%Archive_BudgetTransportFull, 
     &                          State_Diag%Archive_BudgetTransportTrop, 
     &                          State_Diag%Archive_BudgetTransportPBL,  
     &                          State_Diag%BudgetTransportFull,         
     &                          State_Diag%BudgetTransportTrop,         
     &                          State_Diag%BudgetTransportPBL,          
     &                          State_Diag%BudgetMass1,          
     &                          State_Diag%BudgetMass2,          
     &                          RC )
         IF ( RC /= GC_SUCCESS ) THEN
            ErrMsg = 'Transport budget diagnostics error 2'
            CALL GC_Error( ErrMsg, RC, ThisLoc )
            RETURN
         ENDIF
      ENDIF

      ! Return to calling program
      END SUBROUTINE DO_TRANSPORT
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: do_global_adv
!
! !DESCRIPTION: Subroutine DO\_GLOBAL\_ADV is the driver routine 
!  for TPCORE with the GMAO GEOS-FP or MERRA-2 met fields.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE DO_GLOBAL_ADV( am_I_Root, Input_Opt,  State_Met,
     &                          State_Chm, State_Diag, RC )
!
! !USES:
!
      USE CMN_SIZE_MOD             ! Size parameters        
      USE DAO_MOD,            ONLY : AIRQNT
#if defined( BPCH_DIAG )
      USE CMN_DIAG_MOD             ! NDxx flags
      USE DIAG_MOD,           ONLY : MASSFLEW,  MASSFLNS,  MASSFLUP
#endif
      USE ErrCode_Mod
      USE ERROR_MOD
      USE Input_Opt_Mod,      ONLY : OptInput
      USE State_Chm_Mod,      ONLY : ChmState
      USE State_Diag_Mod,     ONLY : DgnState

      USE State_Met_Mod,      ONLY : MetState
      USE State_Met_Mod,      ONLY : MetState
      USE PhysConstants            ! Physical constants
      USE PJC_PFIX_MOD,       ONLY : DO_PJC_PFIX
      USE TIME_MOD,           ONLY : GET_TS_DYN
      USE TPCORE_FVDAS_MOD,   ONLY : TPCORE_FVDAS
      USE UnitConv_Mod
!
! !INPUT PARAMETERS:
!
      LOGICAL,        INTENT(IN)    :: am_I_Root   ! Are we on the root CPU?
      TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
!
! !INPUT/OUTPUT PARAMETERS:
!
      TYPE(MetState), INTENT(INOUT) :: State_Met   ! Meteorology State object
      TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
      TYPE(DgnState), INTENT(INOUT) :: State_Diag  ! Diagnostics State object
!
! !OUTPUT PARAMETERS:
!
      INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
! 
! !REMARKS:
!  As of July 2016, we assume that all of the advected species are listed
!  first in the species database.  This is the easiest way to pass a slab
!  to the TPCORE routine.  This may change in the future. (bmy, 7/13/16)
!                                                                             .
!  Note: the mass flux diagnostic arrays (MASSFLEW, MASSFLNS and MASSFLUP)
!  are incremented upside-down (level 1 = top of the atmosphere).
!  The levels order is reversed only when written out to diagnostic output.
!
! !REVISION HISTORY: 
!  30 Oct 2007 - R. Yantosca - Initial version
!  (1 ) Split off the GEOS-4 & GEOS-5 relevant parts from the previous 
!        routine DO_GLOBAL_TRANSPORT (bmy, 10/30/07)
!  (2 ) Activate the call to SAVE_GLOBAL_TPCORE_BC (yxw, dan, bmy, 11/6/08)
!  (3 ) Bug fix in mass balance: only account for cells of STT with non-zero
!        concentrations when doing the computation (ccc, bmy, 2/17/09)
!  26 Feb 2010 - R. Yantosca - Removed references to obsolete LEMBED switch
!  26 Feb 2010 - R. Yantosca - Added ProTeX headers
!  16 Feb 2011 - R. Yantosca - Add modifications for APM microphysics (G. Luo)
!  21 Jun 2012 - R. Yantosca - Now use pointers to flip indices in vertical
!  09 Nov 2012 - M. Payer    - Replaced all met field arrays with State_Met
!                              derived type object
!  25 Mar 2013 - R. Yantosca - Now accept am_I_Root, Input_Opt, State_Chm, RC
!  04 Jun 2013 - R. Yantosca - Replace pointer STT with State_Chm%Tracers
!  26 Feb 2015 - E. Lundgren - Replace GET_PEDGE with State_Met%PEDGE.
!  31 Mar 2015 - E. Lundgren - Move AIRQNT call from main to within transport
!  09 Jun 2015 - E. Lundgren - Remove dependency on tpcore_bc_mod since not used
!  15 Jun 2015 - E. Lundgren - Tracer units in advection are now kg/kg total air
!                              (previously v/v)
!  15 Jun 2015 - E. Lundgren - Remove mass residual adjustment to tracer
!                              concentration following advection
!  13 Aug 2015 - E. lundgren - Remove v/v <-> kg/kg conversion
!  29 Apr 2016 - R. Yantosca - Don't initialize pointers in declaration stmts
!  30 Jun 2016 - R. Yantosca - Remove instances of STT.  Now get the advected
!                              species ID from State_Chm%Map_Advect.
!  20 Jul 2016 - R. Yantosca - Add minor bug fixes for netCDF diagnostics
!  20 Jul 2016 - R. Yantosca - Now replace NNPAR with State_Chm%nAdvect
!  20 Jul 2016 - R. Yantosca - Now pass State_Chm to DiagnUpdate_Transport_Flux
!  03 Aug 2016 - R. Yantosca - Removed temporary tracer-removal code
!  03 Nov 2017 - R. Yantosca - Now accept State_Diag as an argument
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      ! Scalars
      LOGICAL            :: LFILL
      LOGICAL            :: LPRT
      LOGICAL            :: LWINDO
      INTEGER            :: IORD, JORD, KORD
      INTEGER            :: I, J, L, L2, N, N_DYN, NA, nAdvect
      INTEGER            :: ND24x, ND25x, ND26x
      REAL(fp)           :: D_DYN

      ! Arrays
      REAL(fp)           :: P_TP1      (IIPAR,JJPAR      )
      REAL(fp)           :: P_TP2      (IIPAR,JJPAR      )
      REAL(fp)           :: P_TEMP     (IIPAR,JJPAR      )
      REAL(fp),  TARGET  :: XMASS      (IIPAR,JJPAR,LLPAR) 
      REAL(fp),  TARGET  :: YMASS      (IIPAR,JJPAR,LLPAR) 

      ! Pointers
      REAL(fp),  POINTER :: p_UWND (:,:,:  )
      REAL(fp),  POINTER :: p_VWND (:,:,:  )
      REAL(fp),  POINTER :: p_XMASS(:,:,:  )
      REAL(fp),  POINTER :: p_YMASS(:,:,:  )
      REAL(fp),  POINTER :: p_MFLEW(:,:,:,:)
      REAL(fp),  POINTER :: p_MFLNS(:,:,:,:)
      REAL(fp),  POINTER :: p_MFLUP(:,:,:,:)
      REAL(fp),  POINTER :: p_Spc  (:,:,:,:)

      !=================================================================
      ! DO_GLOBAL_ADV begins here!
      !=================================================================

      ! Assume success
      RC          =  GC_SUCCESS

      ! Initialize
      LFILL       =  Input_Opt%LFILL
      LPRT        =  Input_Opt%LPRT
      LWINDO      =  Input_Opt%LWINDO
      IORD        =  Input_Opt%TPCORE_IORD
      JORD        =  Input_Opt%TPCORE_JORD
      KORD        =  Input_Opt%TPCORE_KORD
      nAdvect     =  State_Chm%nAdvect
      p_MFLEW     => NULL()
      p_MFLNS     => NULL()
      p_MFLUP     => NULL()

      ! Dynamic timestep [s]
      N_DYN       =  GET_TS_DYN()
      D_DYN       =  DBLE( N_DYN )

      ! Define shadow variables for ND24, ND25, ND26
#if defined( BPCH_DIAG )
      ND24x       = ND24
      ND25x       = ND25
      ND26x       = ND26
#else
      ND24x       = 0
      ND25x       = 0
      ND26x       = 0
#endif

      !=================================================================
      ! Prepare variables for calls to PJC pressure-fixer and TPCORE
      !
      ! For hybrid grids, the pressure at the 
      ! bottom edge of grid box (I,J,L) is given by:
      !
      !    P(I,J,L) = Ap(L) + [ Bp(L) * Psurface(I,J) ]
      !
      ! where Psurface is the true surface pressure (i.e. not PS-PTOP).
      ! and Ap(L), Bp(L) define the vertical grid (see pressure_mod.f)
      !=================================================================

      !### DEBUG: Print a few global species sums
      IF ( LPRT ) THEN
         CALL Print_Global_Species_Kg( State_Chm, State_Met, 
     &                                "pre-advection")
      ENDIF
    
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J )
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! Set true dry sfc pressure at midpoint of dynamic timestep [hPa]
         P_TP1(I,J) = GET_PEDGE_DRY(I,J,1)

         ! Set true dry sfc pressure at end of dynamic timestep [hPa]
         P_TP2(I,J) = State_Met%PSC2_DRY(I,J)    
   
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      !=================================================================
      ! Call the PJC/LLNL pressure fixer to get the adjusted air 
      ! masses XMASS and YMASS.  XMASS and YMASS need to be passed to 
      ! TPCORE_FVDAS in order to ensure mass conservation.
      !=================================================================

      ! NOTE: P_TP1 and P_TP2 are the true surface pressures!
      CALL DO_PJC_PFIX( D_DYN,        
     &                  P_TP1,        P_TP2,
     &                  State_Met%U,  State_Met%V,
     &                  XMASS,        YMASS )

      !=================================================================
      ! Call TPCORE_FVDAS to perform the advection
      !=================================================================

      ! Flip array indices in the vertical using pointer storage
      p_UWND    => State_Met%U      ( :, :, LLPAR:1:-1            )
      p_VWND    => State_Met%V      ( :, :, LLPAR:1:-1            )
      p_XMASS   => XMASS            ( :, :, LLPAR:1:-1            )
      p_YMASS   => YMASS            ( :, :, LLPAR:1:-1            )

      ! NOTE: For now, so as to avoid having to rewrite the internals
      ! of the TPCORE routines, just point to 1:nAdvect entries of
      ! State_Chm%Species.  This is OK for now, as of July 2016, all of
      ! the advected species are listed first.  This may change in the
      ! future, but we'll worry about that later.  The units of p_SPC
      ! will be converted to [kg/kg moist air] below. (bmy, 7/13/16)
      p_Spc     => State_Chm%Species( :, :, LLPAR:1:-1, 1:nAdvect )

      ! Diagnostic arrays are vertically flipped at writeout
#if defined( BPCH_DIAG )
      ! Point to global diag_mod arrays for bpch
      p_MFLEW   => MASSFLEW
      p_MFLNS   => MASSFLNS
      p_MFLUP   => MASSFLUP
#endif

      ! Do the advection
      CALL TPCORE_FVDAS( D_DYN,     Re,       IIPAR,    JJPAR,
     &                   LLPAR,     JFIRST,   JLAST,    NG,
     &                   MG,        nAdvect,  Ap,       Bp,
     &                   p_UWND,    p_VWND,   P_TP1,    P_TP2,
     &                   P_TEMP,    p_Spc,    IORD,     JORD,      
     &                   KORD,      N_ADJ,    p_XMASS,  p_YMASS,
     &                   LFILL,
#if defined( BPCH_DIAG )
     &                   p_MFLEW,   p_MFLNS,  p_MFLUP,
#endif
     &                   A_M2,      ND24x,    ND25x,    ND26x,
     &                   State_Diag                            )

      ! Free pointer memory
      p_UWND  => NULL()
      p_VWND  => NULL()
      p_XMASS => NULL()
      p_YMASS => NULL()
      p_Spc   => NULL()
      p_MFLEW => NULL()
      p_MFLNS => NULL()  
      p_MFLUP => NULL()

      !=================================================================
      ! Reset surface pressure and ensure mass conservation
      !=================================================================

      ! Update dry and wet floating pressures to the most recently 
      ! interpolated values (State_Met%PSC2_DRY and State_Met%PSC2)
      ! (ewl, 7/6/16)
      CALL SET_FLOATING_PRESSURES(  am_I_Root, State_Met, RC )

      ! Update State_Met air quantities with new pressures.
      ! Do not update tracer mixing ratio because after advection
      ! the mixing ratio values reflect the new air pressure (ewl, 3/31/15)
      CALL AIRQNT( am_I_Root, Input_opt, State_Met, State_Chm, RC, 
     &             update_mixing_ratio=.FALSE. )

      !### DEBUG: Print a few global species sums
      IF ( LPRT ) THEN
         CALL Print_Global_Species_Kg( State_Chm, State_Met,
     &                                "post-airqnt")
         CALL DEBUG_MSG( '### G4_G5_GLOB_ADV: a TPCORE' ) 
      ENDIF

      END SUBROUTINE DO_GLOBAL_ADV
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: do_window_transport
!
! !DESCRIPTION: Subroutine DO\_WINDOW\_TRANSPORT is the driver program 
!  for the proper TPCORE program for the GEOS-FP/MERRA2 nested-grid
!  simulations. 
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE DO_WINDOW_TRANSPORT( am_I_Root, Input_Opt,  State_Met,
     &                                State_Chm, State_Diag, RC        )
!
! !USES:
!
      USE CMN_SIZE_MOD               ! Size parameters
      USE DAO_MOD,              ONLY : AIRQNT
#if defined( BPCH_DIAG )
      USE CMN_DIAG_MOD               ! NDxx flags
      USE DIAG_MOD,             ONLY : MASSFLEW, MASSFLNS, MASSFLUP
#endif
      USE ErrCode_Mod
      USE ERROR_MOD
      USE GC_GRID_MOD,          ONLY : GET_XOFFSET,   GET_YOFFSET
      USE Input_Opt_Mod,        ONLY : OptInput
      USE State_Chm_Mod,        ONLY : ChmState
      USE State_Met_Mod,        ONLY : MetState
      USE State_Diag_Mod,       ONLY : DgnState
      USE PhysConstants              ! Physical constants
      USE PJC_PFIX_WINDOW_MOD,  ONLY : DO_PJC_PFIX_WINDOW
      USE TIME_MOD,             ONLY : GET_TS_DYN
      USE TPCORE_BC_MOD,        ONLY : I0_W, J0_W, I1_W, J1_W
      USE TPCORE_BC_MOD,        ONLY : I2_W, J2_W, IM_W, JM_W, IGZD
      USE TPCORE_BC_MOD,        ONLY : DO_WINDOW_TPCORE_BC 
      USE TPCORE_WINDOW_MOD,    ONLY : TPCORE_WINDOW
      USE UnitConv_Mod
!
! !INPUT PARAMETERS:
!
      LOGICAL,        INTENT(IN)    :: am_I_Root   ! Are we on the root CPU?
      TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
!
! !INPUT/OUTPUT PARAMETERS:
!
      TYPE(MetState), INTENT(INOUT) :: State_Met   ! Meteorology State object
      TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
      TYPE(DgnState), INTENT(INOUT) :: State_Diag  ! Diagnostics State object
!
! !OUTPUT PARAMETERS:
!
      INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
!
! !REMARKS:
!  As of July 2016, we assume that all of the advected species are listed
!  first in the species database.  This is the easiest way to pass a slab
!  to the TPCORE routine.  This may change in the future. (bmy, 7/13/16)
!                                                                             .
!  Note: the mass flux diagnostic arrays (MASSFLEW, MASSFLNS and MASSFLUP)
!  are incremented upside-down (level 1 = top of the atmosphere).
!  The levels order is reversed only when written out to diagnostic output.
!
! !REVISION HISTORY: 
!  10 Mar 2003 - R. Yantosca - Initial version
!  26 Feb 2010 - R. Yantosca - Removed references to obsolete LEMBED switch
!  26 Feb 2010 - R. Yantosca - Added ProTeX headers
!  16 Feb 2011 - R. Yantosca - Add modifications for APM from G. Luo
!  21 Jun 2012 - R. Yantosca - Now use pointers to flip arrays in the vertical
!  04 Jun 2013 - R. Yantosca - Replace pointer STT with State_Chm%Tracers
!  26 Feb 2015 - E. Lundgren - Replace GET_PEDGE with State_Met%PEDGE.
!  31 Mar 2015 - E. Lundgren - Move AIRQNT call from main to within transport
!  01 Apr 2015 - L. Zhang    - Add updates to skip the window region
!  15 Apr 2015 - R. Yantosca - Add more pointers to avoid array temporaries
!  15 Jun 2015 - E. Lundgren - Tracer units in advection are now kg/kg total air
!                              (previously v/v)
!  15 Jun 2015 - E. Lundgren - Remove mass residual adjustment to tracer
!                              concentration prior to advection
!  13 Aug 2015 - E. lundgren - Remove v/v <-> kg/kg conversion
!  04 Nov 2015 - M. Sulprizio- Rename from DO_GEOSFP_WINDOW_TRANSPORT to
!                              DO_WINDOW_TRANSPORT for use with all nested grids
!  29 Apr 2016 - R. Yantosca - Don't initialize pointers in declaration stmts
!  30 Jun 2016 - R. Yantosca - Remove instances of STT.  Now get the advected
!                              species ID from State_Chm%Map_Advect.
!  20 Jul 2016 - R. Yantosca - Add minor bug fixes for netCDF diagnostics
!  20 Jul 2016 - R. Yantosca - Now replace NNPAR with State_Chm%nAdvect
!  20 Jul 2016 - R. Yantosca - Now pass State_Chm to DiagnUpdate_Transport_Flux
!  03 Aug 2016 - R. Yantosca - Removed temporary tracer-removal code
!  11 Jul 2017 - R. Yantosca - Prevent segfault if ND24, ND25, or ND26 are off
!  03 Nov 2017 - R. Yantosca - Now accept State_Diag as an argument
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      ! Scalars
      LOGICAL            :: LFILL
      LOGICAL            :: LPRT
      LOGICAL            :: LWINDO
      INTEGER            :: IORD,  JORD,  KORD
      INTEGER            :: I0,    J0,    NA,    nAdvect, N_DYN 
      INTEGER            :: IM_W1, JM_W1, I0_W1, J0_W1,   BUFF_SIZE
      INTEGER            :: I,     J,     L,     L2,      N
      INTEGER            :: ND24x, ND25x, ND26x
      REAL(fp)           :: D_DYN

      ! Arrays
      REAL(fp),  TARGET  :: P_TP1      (IIPAR,JJPAR      )
      REAL(fp),  TARGET  :: P_TP2      (IIPAR,JJPAR      )
      REAL(fp),  TARGET  :: P_TEMP     (IIPAR,JJPAR      )
      REAL(fp),  TARGET  :: XMASS      (IIPAR,JJPAR,LLPAR)
      REAL(fp),  TARGET  :: YMASS      (IIPAR,JJPAR,LLPAR)

      ! Pointers
      REAL(fp),  POINTER :: p_A_M2  (  :    )
      REAL(fp),  POINTER :: p_P_TP1 (:,:    )
      REAL(fp),  POINTER :: p_P_TP2 (:,:    )
      REAL(fp),  POINTER :: p_P_TEMP(:,:    )
      REAL(fp),  POINTER :: p_UWND  (:,:,:  )
      REAL(fp),  POINTER :: p_VWND  (:,:,:  )
      REAL(fp),  POINTER :: p_XMASS (:,:,:  )
      REAL(fp),  POINTER :: p_YMASS (:,:,:  )
      REAL(fp),  POINTER :: p_MFLEW (:,:,:,:)
      REAL(fp),  POINTER :: p_MFLNS (:,:,:,:)
      REAL(fp),  POINTER :: p_MFLUP (:,:,:,:)
      REAL(fp),  POINTER :: p_Spc   (:,:,:,:)

      !=================================================================
      ! DO_FVDAS_WINDOW_TRANSPORT begins here!
      !=================================================================

      ! Assume success
      RC          =  GC_SUCCESS

      ! Copy values from Input_Opt
      LFILL       =  Input_Opt%LFILL
      LPRT        =  Input_Opt%LPRT
      LWINDO      =  Input_Opt%LWINDO
      IORD        =  Input_Opt%TPCORE_IORD
      JORD        =  Input_Opt%TPCORE_JORD
      KORD        =  Input_Opt%TPCORE_KORD
      nAdvect     =  State_Chm%nAdvect

      ! Initialize pointers
      p_A_M2      => NULL()
      p_P_TP1     => NULL()
      p_P_TP2     => NULL()
      p_P_TEMP    => NULL()
      p_UWND      => NULL()
      p_VWND      => NULL()
      p_XMASS     => NULL()
      p_YMASS     => NULL()
      p_MFLEW     => NULL()
      p_MFLNS     => NULL()
      p_MFLUP     => NULL()
      p_Spc       => NULL()

      ! Get nested-grid lon/lat offsets [# boxes]
      I0          =  GET_XOFFSET( GLOBAL=.TRUE. )
      J0          =  GET_YOFFSET( GLOBAL=.TRUE. )

      ! Dynamic timestep [s]
      N_DYN       =  GET_TS_DYN()
      D_DYN       =  DBLE( N_DYN )

      ! (lzh, 09/01/2014)
      BUFF_SIZE   =  2
      IM_W1       =  IM_W + 2 * BUFF_SIZE
      JM_W1       =  JM_W + 2 * BUFF_SIZE
      I0_W1       =  I0_W     - BUFF_SIZE
      J0_W1       =  J0_W     - BUFF_SIZE

      ! Define shadow variables for ND24, ND25, ND26
#if defined( BPCH_DIAG )
      ND24x       = ND24
      ND25x       = ND25
      ND26x       = ND26
#else
      ND24x       = 0
      ND25x       = 0
      ND26x       = 0
#endif

      !=================================================================
      ! Prepare variables for calls to PJC pressure-fixer and TPCORE
      !
      ! For hybrid grids, the pressure at the
      ! bottom edge of grid box (I,J,L) is given by:
      !
      !    P(I,J,L) = Ap(L) + [ Bp(L) * Psurface(I,J) ]
      !
      ! where Psurface is the true surface pressure (i.e. not PS-PTOP).
      ! and Ap(L), Bp(L) define the vertical grid (see pressure_mod.f)
      !=================================================================

      IF ( LPRT ) THEN
         CALL Print_Global_Species_Kg( State_Chm, State_Met, 
     &                                "pre-advection")
      ENDIF

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J )
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! Set true dry sfc pressure at midpoint of dynamic timestep [hPa]
         P_TP1(I,J) = GET_PEDGE_DRY(I,J,1)

         ! Set true dry sfc pressure at end of dynamic timestep [hPa]
         P_TP2(I,J) = State_Met%PSC2_DRY(I,J)    

      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      !=================================================================
      ! Call the PJC/LLNL pressure fixer to get the adjusted air
      ! masses XMASS and YMASS.  XMASS and YMASS need to be passed to
      ! TPCORE_FVDAS in order to ensure mass conservation.
      !=================================================================
      XMASS = 0e+0_fp !(dan)
      YMASS = 0e+0_fp
      ! NOTE: P_TP1 and P_TP2 are the true surface pressures!  
      CALL DO_PJC_PFIX_WINDOW( D_DYN,
     &                         P_TP1,        P_TP2,
     &                         State_Met%U,  State_Met%V,
     &                         XMASS,        YMASS )

      IF ( LPRT ) CALL DEBUG_MSG( '### FVDAS_WINDOW: a PJC_PFIX_WINDOW')

      ! Impose TPCORE boundary conditions @ edges of nested grid
      CALL DO_WINDOW_TPCORE_BC( am_I_Root, Input_Opt, State_Chm, RC )

      ! Flip array indices in the vertical using pointer storage

      ! Exclude the buffer zone (lzh, 4/1/2015)
      p_UWND  => State_Met%U(       I0_W1+1 : I0_W1+IM_W1, 
     &                              J0_W1+1 : J0_W1+JM_W1, 
     &                              LLPAR:1:-1                   )
                                    
      p_VWND  => State_Met%V(       I0_W1+1 : I0_W1+IM_W1,
     &                              J0_W1+1 : J0_W1+JM_W1, 
     &                              LLPAR:1:-1                   )

      ! NOTE: For now, so as to avoid having to rewrite the internals
      ! of the TPCORE routines, just point to 1:nAdvect entries of
      ! State_Chm%Species.  This is OK for now, as of July 2016, all of
      ! the advected species are listed first.  This may change in the
      ! future, but we'll worry about that later. (bmy, 7/13/16)
      p_Spc   => State_Chm%Species( I0_W1+1 : I0_W1+IM_W1,
     &                              J0_W1+1 : J0_W1+JM_W1, 
     &                              LLPAR:1:-1,                
     &                              1:nAdvect                    )             

      p_XMASS => XMASS(             I0_W1+1 : I0_W1+IM_W1,
     &                              J0_W1+1 : J0_W1+JM_W1,
     &                              LLPAR:1:-1                   )

      p_YMASS => YMASS(             I0_W1+1 : I0_W1+IM_W1, 
     &                              J0_W1+1 : J0_W1+JM_W1,
     &                              LLPAR:1:-1                   )
                                    
      p_P_TP1 => P_TP1(             I0_W1+1 : I0_W1+IM_W1, 
     &                              J0_W1+1 : J0_W1+JM_W1        )
                                    
      p_P_TP2 => P_TP2(             I0_W1+1 : I0_W1+IM_W1, 
     &                              J0_W1+1 : J0_W1+JM_W1        )

      p_P_TEMP => P_TEMP(           I0_W1+1 : I0_W1+IM_W1,
     &                              J0_W1+1 : J0_W1+JM_W1        )
                                   
#if defined( BPCH_DIAG )
      !---------------------------------------------------------------------
      ! NOTE: If ND24 is off, then MASSFLEW is sized (1,1,1,1), and
      !       If ND25 is off, then MASSFLNS is sized (1,1,1,1), and
      !       If ND26 is off, then MASSFLUP is sized (1,1,1,1).
      !
      ! Using the window offsets (which are greater than 1) will cause
      ! an array-out-of-bounds error.  For now if any of ND24, ND25,
      ! or ND26 is off, just point to the whole array and don't try
      ! to subset it
      !
      ! NOTE: BPCH diagnostics are going away really soon!!! 
      ! (bmy, 7/11/17)
      !---------------------------------------------------------------------
      IF ( Input_Opt%ND24 > 0 ) THEN
         p_MFLEW => MASSFLEW(       I0_W1+1 : I0_W1+IM_W1,
     &                              J0_W1+1 : J0_W1+JM_W1,  
     &                              :,                         : )
      ELSE
         p_MFLEW => MASSFLEW
      ENDIF

      IF ( Input_Opt%ND25 > 0 ) THEN
         p_MFLNS => MASSFLNS(       I0_W1+1 : I0_W1+IM_W1,
     &                              J0_W1+1 : J0_W1+JM_W1,  
     &                              :,                         : )
      ELSE
         p_MFLNS => MASSFLNS
      ENDIF

      IF ( Input_Opt%ND26 > 0 ) THEN
         p_MFLUP => MASSFLUP(       I0_W1+1 : I0_W1+IM_W1,
     &                              J0_W1+1 : J0_W1+JM_W1,  
     &                              :,                         : )
      ELSE
         p_MFLUP => MASSFLUP
      ENDIF
#endif


      p_A_M2  => A_M2( J0_W1+1 : J0_W1+JM_W1                     )

      ! Do the advection
      CALL TPCORE_WINDOW( 
     &      D_DYN,   Re,       IM_W1,   JM_W1,   LLPAR, 
     &      JFIRST,  JLAST,    NG,      MG,      nAdvect, 
     &      Ap,      Bp,       p_UWND,  p_VWND,  p_P_TP1, 
     &      p_P_TP2, p_P_TEMP, p_Spc,   IORD,    JORD, 
     &      KORD,    N_ADJ,    p_XMASS, p_YMASS,
#if defined( BPCH_DIAG )
     &      p_MFLEW, p_MFLNS,  p_MFLUP,
#endif
     &      p_A_M2,  ND24x,    ND25x,   ND26x                    ) 

      ! Free pointer memory
      p_UWND   => NULL()  
      p_VWND   => NULL()
      p_Spc    => NULL()
      p_XMASS  => NULL()
      p_YMASS  => NULL()
      p_P_TP1  => NULL()
      p_P_TP2  => NULL()
      p_P_TEMP => NULL()
      p_A_M2   => NULL()
      p_MFLEW  => NULL()
      p_MFLNS  => NULL() 
      p_MFLUP  => NULL() 

      !=================================================================
      ! Reset surface pressure and ensure mass conservation
      !=================================================================

      ! Update dry and wet floating pressures to the most recently 
      ! interpolated values (State_Met%PSC2_DRY and State_Met%PSC2)
      ! (ewl, 7/6/16)
      CALL SET_FLOATING_PRESSURES(  am_I_Root, State_Met, RC )

      ! Update State_Met air quantities with new pressures.
      ! Do not update tracer mixing ratio because after advection
      ! the mixing ratio values reflect the new air pressure (ewl, 3/31/15)
      CALL AIRQNT( am_I_Root, Input_opt, State_Met, State_Chm, RC, 
     &             update_mixing_ratio=.FALSE. )

      !### Debug
      IF ( LPRT ) THEN
         CALL Print_Global_Species_Kg( State_Chm, State_Met,
     &                                "post-airqnt")
         CALL DEBUG_MSG( '### NESTED_ADV: a TPCORE' ) 
      ENDIF

      END SUBROUTINE DO_WINDOW_TRANSPORT
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_transport
!
! !DESCRIPTION: Subroutine INIT\_TRANSPORT initializes all module variables 
!  and arrays.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE INIT_TRANSPORT( am_I_Root, Input_Opt, State_Diag, RC )
!
! !USES:
!
      USE CMN_SIZE_MOD           ! Size parameters
      USE ErrCode_Mod
      USE ERROR_MOD,        ONLY : ALLOC_ERR
      USE GC_GRID_MOD,      ONLY : GET_AREA_M2, GET_YMID_R
      USE Input_Opt_Mod,    ONLY : OptInput
      USE State_Diag_Mod,   ONLY : DgnState
      USE PhysConstants          ! Re
      USE TIME_MOD,         ONLY : GET_TS_DYN
      USE TPCORE_FVDAS_MOD, ONLY : INIT_TPCORE
!
! !INPUT PARAMETERS:
!
      LOGICAL,        INTENT(IN)  :: am_I_Root   ! Are we on the root CPU?
      TYPE(OptInput), INTENT(IN)  :: Input_Opt   ! Input Options object
      TYPE(DgnState), INTENT(IN)  :: State_Diag  ! Diagnostics State object
!
! !OUTPUT PARAMETERS:
!
      INTEGER,        INTENT(OUT) :: RC          ! Success or failure?
! 
! !REVISION HISTORY: 
!  10 Mar 2003 - R. Yantosca - Initial version
!  (1 ) Now references GET_TS_DYN from "time_mod.f", INIT_TPCORE_FVDAS from
!        "tpcore_fvdas_mod.f90", and GET_YMID_R from "grid_mod.f".  Now also
!        include "CMN_SETUP".  (bdf, bmy, 4/28/03)
!  (2 ) Remove reference to DSIG, it's obsolete. (bmy, 6/24/03)
!  (3 ) Now references LEMBED & LTPFV from "logical_mod.f".  Now references
!        N_TRACERS from "tracer_mod.f". (bmy, 7/20/04)
!  (4 ) Now modified for GEOS-5 and GCAP met fields (swu, bmy, 5/25/05)
!  (5 ) Removed reference to USE_GEOS_4_TRANSPORT, STT_I1, STT_I2, STT_J1,
!        STT_J2, variables (bmy, 10/30/07)
!  (6 ) Deleted reference to CMN, it's not needed anymore (bmy, 11/6/08)
!  26 Feb 2010 - R. Yantosca - Removed references to obsolete LEMBED switch
!  26 Feb 2010 - R. Yantosca - Added ProTeX headers
!  28 Feb 2012 - R. Yantosca - Removed support for GEOS-3
!  01 Mar 2012 - R. Yantosca - Now use GET_AREA_M2(I,J,L) from grid_mod.F90
!  01 Mar 2012 - R. Yantosca - Now use GET_YMID_R_W(I,J,L) from grid_mod.F90
!  25 Mar 2013 - R. Yantosca - Now accept am_I_Root, Input_Opt, RC
!  05 Dec 2014 - M. Yannetti - Added REAL_N_DYN
!  03 Nov 2017 - R. Yantosca - Now accept State_Diag as an argument
!  09 Nov 2017 - R. Yantosca - Return error condition to calling program
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      ! Scalars
      LOGICAL            :: LTRAN
      INTEGER            :: J, K, L, N_DYN
      REAL(fp)           :: YMID_R(JJPAR)
      REAL(fp)           :: REAL_N_DYN

      ! Strings
      CHARACTER(LEN=255) :: ErrMsg, ThisLoc

      !=================================================================
      ! Initialize
      !=================================================================
      RC      = GC_SUCCESS
      ErrMsg  = ''
      ThisLoc = 
     &   ' -> at Init_Transport  (in module GeosCore/transport_mod.F)'

      !=================================================================
      ! Allocate arrays for TPCORE vertical coordinates 
      !
      ! For fvDAS TPCORE with for GEOS-FP or MERRA-2 met fields:
      !
      !    P(I,J,L) = Ap(L) + ( Bp(L) * Psurf(I,J) )
      !
      ! Also here Ap, Bp will be flipped since both TPCORE versions
      ! index levels from the atm. top downwards (bdf, bmy, 10/30/07)
      !=================================================================
      ALLOCATE( Ap( LLPAR+1 ), STAT=RC )
      CALL GC_CheckVar( 'transport_mod.F:Ap', 0, RC )
      IF ( RC /= GC_SUCCESS ) RETURN
 
      ALLOCATE( Bp( LLPAR+1 ), STAT=RC )
      CALL GC_CheckVar( 'transport_mod.F:Bp', 0, RC )
      IF ( RC /= GC_SUCCESS ) RETURN

      ! Flip Ap and Bp for TPCORE
      DO L = 1, LLPAR+1 

         ! As L runs from the surface up, 
         ! K runs from the top down
         K = ( LLPAR + 1 ) - L + 1

         Ap(L) = GET_AP(K)          ! Ap(L) is in [hPa] 
         Bp(L) = GET_BP(K)
      ENDDO

      !=================================================================
      ! Allocate arrays for surface area and layer thickness
      !=================================================================
      ALLOCATE( A_M2( JJPAR ), STAT=RC )
      CALL GC_CheckVar( 'transport_mod.F:A_m2', 0, RC )
      IF ( RC /= GC_SUCCESS ) RETURN

      ! Surface area [m2]
      DO J = 1, JJPAR
         A_M2(J) = GET_AREA_M2( 1, J, 1 )
      ENDDO

      !=================================================================
      ! Additional setup for the fvDAS version of TPCORE
      !=================================================================

      ! Initialize
      N_DYN = GET_TS_DYN()
      N_ADJ = 0
      NG    = 0
      MG    = 0

      ! YMID_R is latitude of grid box center [radian]
      DO J = 1,JJPAR
         YMID_R(J) = GET_YMID_R( 1, J, 1 )
      ENDDO

      REAL_N_DYN = N_DYN

      ! Call INIT routine from "tpcore_fvdas_mod.f" 
      CALL INIT_TPCORE( IIPAR, JJPAR,  LLPAR,      JFIRST, 
     &                  JLAST, NG,     MG,         REAL_N_DYN, 
     &                  Re,    YMID_R, State_Diag, RC          )

      ! Trap potential errors
      IF ( RC /= GC_SUCCESS ) THEN
         ErrMsg = 'Error encountered in "Init_Tpcore"!'
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF

      END SUBROUTINE INIT_TRANSPORT
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_window_transport
!
! !DESCRIPTION: Subroutine INIT\_WINDOW\_TRANSPORT initializes all 
!  module variables and arrays for the GEOS-FP/MERRA2 nested grid
!  simulation.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE INIT_WINDOW_TRANSPORT( am_I_Root,  Input_Opt,
     &                                  State_Diag, RC         )
!
! !USES:    
!
      USE CMN_SIZE_MOD             ! Size parameters
      USE ErrCode_Mod
      USE ERROR_MOD,          ONLY : ALLOC_ERR
      USE GC_GRID_MOD,        ONLY : GET_AREA_M2
      USE GC_GRID_MOD,        ONLY : GET_YMID_R_W
      USE Input_Opt_Mod,      ONLY : OptInput
      USE PhysConstants            ! Re
      USE State_Diag_Mod,     ONLY : DgnState
      USE TIME_MOD,           ONLY : GET_TS_DYN
      USE TPCORE_FVDAS_MOD,   ONLY : INIT_TPCORE
      USE TPCORE_BC_MOD,      ONLY : I0_W, J0_W, I1_W, J1_W
      USE TPCORE_BC_MOD,      ONLY : I2_W, J2_W, IM_W, JM_W
      USE TPCORE_BC_MOD,      ONLY : IGZD, INIT_TPCORE_BC
      USE TPCORE_WINDOW_MOD,  ONLY : INIT_WINDOW
!
! !INPUT PARAMETERS:
!
      LOGICAL,        INTENT(IN)  :: am_I_Root   ! Are we on the root CPU?
      TYPE(OptInput), INTENT(IN)  :: Input_Opt   ! Input Options object
      TYPE(DgnState), INTENT(IN)  :: State_Diag  ! Diagnostics State object
!
! !OUTPUT PARAMETERS:
!
      INTEGER,        INTENT(OUT) :: RC          ! Success or failure?
!
! !REVISION HISTORY: 
!  06 Jun 2008 - D. Chen & R. Yantosca - Initial version
!  26 Feb 2010 - R. Yantosca - Removed references to obsolete LEMBED switch
!  26 Feb 2010 - R. Yantosca - Added ProTeX headers
!  01 Mar 2012 - R. Yantosca - Now use GET_AREA_CM2(I,J,L) from grid_mod.F90
!  01 Mar 2012 - R. Yantosca - Now use GET_YMID_R_W(I,J,L) from grid_mod.F90
!  25 Mar 2013 - R. Yantosca - Now accept am_I_Root, Input_Opt, State_Chm, RC
!  26 Sep 2013 - R. Yantosca - Renamed to INIT_GEOSFP_WINDOW_TRANSPORT
!  05 Dec 2014 - M. Yannetti - Added REAL_N_DYN
!  04 Nov 2015 - M. Sulprizio- Renamed from INIT_GEOSFP_WINDOW_TRANSPORT to
!                              INIT_WINDOW_TRANSPORT for use with all nested
!                              grids
!  26 Jun 2017 - R. Yantosca - GC_ERROR is now contained in errcode_mod.F90
!  03 Nov 2017 - R. Yantosca - Now accept State_Diag as an argument
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      ! Scalars
      LOGICAL  :: LTRAN
      INTEGER  :: BUFF_SIZE
      INTEGER  :: J,     K,     L,     N_DYN
      INTEGER  :: IM_W1, JM_W1, I0_W1, J0_W1
      REAL(fp) :: REAL_N_DYN

      ! Arrays
      REAL(fp) :: YMID_R_W(0:JJPAR+1)

      !=================================================================
      ! Initialize
      !=================================================================

      ! Assume success
      RC        =  GC_SUCCESS

      ! Copy values from Input_Opt
      LTRAN     = Input_Opt%LTRAN

      ! Cast N_DYN to flexible precision
#if defined( USE_REAL8 )
      REAL_N_DYN = DBLE( N_DYN )
#else
      REAL_N_DYN = FLOAT( N_DYN )
#endif

      !=================================================================
      ! Allocate arrays for TPCORE vertical coordinates
      ! GEOS-FP/MERRA2 nested grid simulation only!!!
      !
      ! For fvDAS TPCORE with for GEOS-FP/MERRA2 met fields:
      !
      !    P(I,J,L) = Ap(L) + ( Bp(L) * Psurf(I,J) )
      !
      ! Also here Ap, Bp will be flipped since both TPCORE versions
      ! index levels from the atm. top downwards (bdf, bmy, 10/30/07)
      !=================================================================
      ALLOCATE( Ap( LLPAR+1 ), STAT=RC )
      IF ( RC /= 0 ) CALL ALLOC_ERR( 'Ap' )

      ALLOCATE( Bp( LLPAR+1 ), STAT=RC )
      IF ( RC /= 0 ) CALL ALLOC_ERR( 'Bp' )

      ! Flip Ap and Bp for TPCORE
      DO L = 1, LLPAR+1

         ! As L runs from the surface up,
         ! K runs from the top down
         K = ( LLPAR + 1 ) - L + 1

         Ap(L) = GET_AP(K)
         Bp(L) = GET_BP(K)
      ENDDO

      !=================================================================
      ! Allocate arrays for surface area and layer thickness
      !=================================================================
      ALLOCATE( A_M2( JJPAR ), STAT=RC )
      IF ( RC /= 0 ) CALL ALLOC_ERR( 'A_M2' )

      ! Surface area [m2]
      DO J = 1, JJPAR
         A_M2(J) = GET_AREA_M2( 1, J, 1 )
      ENDDO

      !=================================================================
      ! Additional setup for the fvDAS version of TPCORE
      !=================================================================

      ! Initialize
      N_DYN = GET_TS_DYN()
      N_ADJ = 0
      NG    = 0
      MG    = 0

      ! (lzh, 4/1/2015)
      BUFF_SIZE = 2
      IM_W1     = IM_W + 2 * BUFF_SIZE
      JM_W1     = JM_W + 2 * BUFF_SIZE
      I0_W1     = I0_W - BUFF_SIZE
      J0_W1     = J0_W - BUFF_SIZE

      ! YMID_R is latitude of grid box center [radians]
      DO J =0, JJPAR+1
         YMID_R_W(J) = GET_YMID_R_W( 1, J, 1 )
      ENDDO

      REAL_N_DYN = N_DYN

      ! Call INIT routine from "tpcore_window_mod.F90"
      CALL INIT_WINDOW( IM_W1, 
     &                  JM_W1,    
     &                  LLPAR, 
     &                  JFIRST, 
     &                  JLAST, 
     &                  NG,       
     &                  MG,    
     &                  REAL_N_DYN,
     &                  Re,  
     &                  YMID_R_W( J0_W1:(J0_W1+JM_W1+1) ) )

      END SUBROUTINE INIT_WINDOW_TRANSPORT
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: cleanup_transport
!
! !DESCRIPTION: Subroutine CLEANUP\_TRANSPORT deallocates all module arrays. 
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE CLEANUP_TRANSPORT
! 
! !REVISION HISTORY: 
!  10 Mar 2003 - R. Yantosca - Initial version
!  (1 ) Remove reference to DSIG, it's obsolete. (bmy, 6/24/03)
!  (2 ) Remove obsolete embedded chemistry arrays (bmy, 10/30/07)
!  26 Feb 2010 - R. Yantosca - Added ProTeX headers
!  12 Feb 2015 - E. Lundgren - Added new diagnostics arrays for netcdf output
!EOP
!------------------------------------------------------------------------------
!BOC
!
      IF ( ALLOCATED( Ap     ) ) DEALLOCATE( Ap     )
      IF ( ALLOCATED( A_M2   ) ) DEALLOCATE( A_M2   )
      IF ( ALLOCATED( Bp     ) ) DEALLOCATE( Bp     )

      END SUBROUTINE CLEANUP_TRANSPORT
!EOC

      END MODULE TRANSPORT_MOD
