! $Id: chemistry_mod.f,v 1.32 2009/09/15 15:51:48 phs Exp $
      MODULE CHEMISTRY_MOD
!
!******************************************************************************
!  Module CHEMISTRY_MOD is used to call the proper chemistry subroutine
!  for the various GEOS-CHEM simulations. (bmy, 4/14/03, 4/2/08)
! 
!  Module Routines:
!  ============================================================================
!  (1 ) DO_CHEMISTRY       : Driver which calls various chemistry routines
!  (2 ) GCKPP_DRIVER       : KPP chemistry driver subroutine
!
!  GEOS-CHEM modules referenced by chemistry_mod.f
!  ============================================================================
!  (1 ) acetone_mod.f      : Module w/ routines for ACET chemistry
!  (2 ) c2h6_mod.f         : Module w/ routines for C2H6 chemistry
!  (3 ) carbon_mod.f       : Module w/ routines for carbon arsl chem.
!  (4 ) ch3i_mod.f         : Module w/ routines for CH3I chemistry
!  (5 ) dao_mod.f          : Module w/ arrays for DAO met fields
!  (6 ) diag_pl_mod.f      : Module w/ routines to save P(Ox), L(Ox)
!  (7 ) drydep_mod.f       : Module w/ GEOS-CHEM drydep routines
!  (8 ) dust_mod.f         : Module w/ routines for dust arsl chem.
!  (9 ) error_mod.f        : Module w/ NaN and error checks
!  (10) global_ch4_mod.f   : Module w/ routines for CH4 chemistry
!  (11) hcn_ch3cn_mod.f    : Module w/ routines for HCN and CH3CN chemistry
!  (12) Kr85_mod.f         : Module w/ routines for Kr85 chemistry
!  (13) logical_mod.f      : Module w/ GEOS-CHEM logical switches
!  (14) RnPbBe_mod.f       : Module w/ routines for Rn-Pb-Be chemistry
!  (15) rpmares_mod.f      : Module w/ routines for arsl phase equilib.
!  (16) seasalt_mod.f      : Module w/ routines for seasalt chemistry
!  (17) sulfate_mod.f      : Module w/ routines for sulfate chemistry
!  (18) tagged_co_mod.f    : Module w/ routines for Tagged CO chemistry
!  (19) tagged_ox_mod.f    : Module w/ routines for Tagged Ox chemistry
!  (20) time_mod.f         : Module w/ routines to compute time & date
!  (21) tracer_mod.f       : Module w/ GEOS-CHEM tracer array STT etc. 
!  (22) tracerid_mod.f     : Module w/ pointers to tracers & emissions
!
!  NOTES:
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
!******************************************************************************
!
      IMPLICIT NONE

      !=================================================================
      ! MODULE ROUTINES -- follow below the "CONTAINS" statement
      !=================================================================
      CONTAINS

!------------------------------------------------------------------------------
      
      SUBROUTINE DO_CHEMISTRY
!
!******************************************************************************
!  Subroutine DO_CHEMISTRY is the driver routine which calls the appropriate
!  chemistry subroutine for the various GEOS-CHEM simulations. 
!  (bmy, 2/11/03, 9/18/07)
!
!  NOTES:
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
!******************************************************************************
!
      ! References to F90 modules
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
      USE ISOROPIA_MOD,    ONLY : DO_ISOROPIA
      USE Kr85_MOD,        ONLY : CHEMKr85
      USE LOGICAL_MOD,     ONLY : LCARB, LCHEM,  LCRYST, LDUST
      USE LOGICAL_MOD,     ONLY : LPRT,  LSSALT, LSULF,  LSOA
      USE MERCURY_MOD,     ONLY : CHEMMERCURY
      USE OPTDEPTH_MOD,    ONLY : OPTDEPTH
      USE RnPbBe_MOD,      ONLY : CHEMRnPbBe
      USE RPMARES_MOD,     ONLY : DO_RPMARES
      USE SEASALT_MOD,     ONLY : CHEMSEASALT
      USE SULFATE_MOD,     ONLY : CHEMSULFATE
      USE TAGGED_CO_MOD,   ONLY : CHEM_TAGGED_CO
      USE TAGGED_OX_MOD,   ONLY : CHEM_TAGGED_OX
      USE TIME_MOD,        ONLY : GET_ELAPSED_MIN, GET_TS_CHEM
      USE TRACER_MOD,      ONLY : N_TRACERS,       STT  
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
      USE TRACERID_MOD,    ONLY : IDTACET, IDTISOP
      USE LOGICAL_MOD,     ONLY : LNLPBL ! (Lin, 03/31/09)

#     include "CMN_SIZE"        ! Size parameters
#     include "CMN_DIAG"        ! NDxx flags
#     include "comode.h"        ! NPHOT

      ! Local variables
      LOGICAL, SAVE            :: FIRST = .TRUE.
      INTEGER                  :: N_TROP

      !=================================================================
      ! DO_CHEMISTRY begins here!
      !=================================================================

      ! Compute optical depths (except for CH4 simulation)
      IF ( .not. ITS_A_CH4_SIM() ) THEN
         CALL OPTDEPTH( LLPAR, CLDF, OPTDEP, OPTD )
      ENDIF

      !=================================================================
      ! If LCHEM=T then call the chemistry subroutines
      !=================================================================
      IF ( LCHEM ) THEN 

         !---------------------------------
         ! NOx-Ox-HC (w/ or w/o aerosols) 
         !---------------------------------
         IF ( ITS_A_FULLCHEM_SIM() ) THEN 

            ! Call SMVGEAR routines
            CALL CHEMDR

             ! Do seasalt aerosol chemistry
             IF ( LSSALT ) CALL CHEMSEASALT
 
             ! Also do sulfate chemistry
             IF ( LSULF ) THEN
 
                ! Do sulfate chemistry
                CALL CHEMSULFATE
 
                ! Do aerosol thermodynamic equilibrium
                !------------------------------------------------------------
                ! Prior to 4/2/08:
                ! Bug fix: ISORROPIA can return very unphysical values when
                ! RH is very low.  We will replace the current version of
                ! ISORROPIA with ISORROPIA II.  In the meantime, we shall
                ! use RPMARES to do the ATE computations. (bmy, 4/2/08)
                !IF ( LSSALT ) THEN
                !
                !   ! ISOROPIA takes Na+, Cl- into account
                !   CALL DO_ISOROPIA
                !
                !ELSE
 
                   ! RPMARES does not take Na+, Cl- into account
                   CALL DO_RPMARES
 
                !ENDIF
                !------------------------------------------------------------
                
             ENDIF
 
             ! Do carbonaceous aerosol chemistry
             IF ( LCARB ) CALL CHEMCARBON
 
             ! Do dust aerosol chemistry
             IF ( LDUST ) CALL CHEMDUST
 
             ! ND44 drydep fluxes
             ! The drydep fluxes diag. is done in vdiff_mod.f when non-local 
             ! PBL is used ( Lin, 03/31/09) 
             IF (.NOT. LNLPBL) CALL DRYFLX     
 
             ! ND43 chemical production
             CALL DIAGOH
 
             ! Remove acetone ocean sink
             IF ( IDTACET /= 0 ) THEN
                CALL OCEAN_SINK_ACET( STT(:,:,1,IDTACET) ) 
             ENDIF
 
          !---------------------------------
          ! Offline aerosol simulation
          !---------------------------------
          ELSE IF ( ITS_AN_AEROSOL_SIM() ) THEN
 
             ! Define loop index and other SMVGEAR arrays
             ! N_TROP, the # of trop boxes, is returned
             CALL AEROSOL_RURALBOX( N_TROP )
 
             ! Initialize FAST-J quantities for computing AOD's
             IF ( FIRST ) THEN
                CALL READER( FIRST )
                CALL READCHEM
                CALL INPHOT( LLTROP, NPHOT )
 
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
 
             ! Compute AOD's and surface areas
             CALL RDAER
 
             !*** AEROSOL THERMODYNAMIC EQUILIBRIUM ***
             !-------------------------------------------------------------
             ! Prior to 4/2/08:
             ! Bug fix: ISORROPIA can return very unphysical values when
             ! RH is very low.  We will replace the current version of
             ! ISORROPIA with ISORROPIA II.  In the meantime, we shall
             ! use RPMARES to do the ATE computations. (bmy, 4/2/08)
             !IF ( LSSALT ) THEN
             !
             !   ! ISOROPIA takes Na+, Cl- into account
             !   CALL DO_ISOROPIA
             !
             !ELSE
 
                ! RPMARES does not take Na+, Cl- into account
                ! (skip for crystalline & aqueous offline run)
                IF ( .not. LCRYST ) CALL DO_RPMARES
 
             !ENDIF
             !-------------------------------------------------------------
 
             !*** SEASALT AEROSOLS ***
             IF ( LSSALT ) CALL CHEMSEASALT
 
             !*** SULFATE AEROSOLS ***
             IF ( LSULF .or. LCRYST ) THEN
 
                ! Do sulfate chemistry
                CALL CHEMSULFATE
 
             ENDIF
                
             !*** CARBON AND 2NDARY ORGANIC AEROSOLS ***
             IF ( LCARB ) CALL CHEMCARBON
 
             !*** MINERAL DUST AEROSOLS ***
             IF ( LDUST ) THEN 
 
                ! Do dust aerosol chemsitry
                CALL CHEMDUST
 
                ! Compute dust OD's & surface areas
                CALL RDUST_ONLINE( SOILDUST )
             ENDIF
 
          !---------------------------------
          ! Rn-Pb-Be
          !---------------------------------                 
          ELSE IF ( ITS_A_RnPbBe_SIM() ) THEN
             CALL CHEMRnPbBe 
             CALL DRYFLXRnPbBe
                   
          !---------------------------------
          ! CH3I
          !---------------------------------
          ELSE IF ( ITS_A_CH3I_SIM() ) THEN
             CALL CHEMCH3I
 
          !---------------------------------            
          ! HCN
          !---------------------------------
          ELSE IF ( ITS_A_HCN_SIM() ) THEN
             CALL CHEM_HCN_CH3CN( N_TRACERS, STT )
 
          !---------------------------------
          ! Tagged O3
          !---------------------------------
          ELSE IF ( ITS_A_TAGOX_SIM() ) THEN 
             CALL CHEM_TAGGED_OX
 
          !---------------------------------
          ! Tagged CO
          !---------------------------------
          ELSE IF ( ITS_A_TAGCO_SIM() ) THEN
             CALL CHEM_TAGGED_CO
 
          !---------------------------------
          ! C2H6
          !---------------------------------
          ELSE IF ( ITS_A_C2H6_SIM() ) THEN
             CALL CHEMC2H6
 
          !---------------------------------
          ! CH4
          !---------------------------------
          ELSE IF ( ITS_A_CH4_SIM() ) THEN
 
             ! Only call after the first 24 hours
             IF ( GET_ELAPSED_MIN() >= GET_TS_CHEM() ) THEN
                CALL CHEMCH4
             ENDIF
 
          !---------------------------------
          ! Mercury
          !---------------------------------
          ELSE IF ( ITS_A_MERCURY_SIM() ) THEN
 
             ! Do Hg chemistry
             CALL CHEMMERCURY
               
          !---------------------------------
          ! Offline H2/HD
          !---------------------------------
          ELSE IF ( ITS_A_H2HD_SIM() ) THEN
             CALL CHEM_H2_HD
             CALL DRYFLXH2HD
  
 !-----------------------------------------------------------------------------
 ! Prior to 7/19/04:
 ! Fully install Kr85 run later (bmy, 7/19/04)
 !            !---------------------------------
 !            ! Kr85   
 !            !---------------------------------
 !            CASE ( 12 )
 !               CALL CHEMKr85
 !-----------------------------------------------------------------------------
          ENDIF

         !### Debug
         IF ( LPRT ) CALL DEBUG_MSG( '### MAIN: a CHEMISTRY' )
      ENDIF
         
      ! Return to calling program
      END SUBROUTINE DO_CHEMISTRY

!------------------------------------------------------------------------------
      
      SUBROUTINE GCKPP_DRIVER( ) 
!
!******************************************************************************
!     Driver routine to perform integration of the full KPP chemistry mechanism.
!     Based on Daven Henze's GCKPP_DRIVER.           (Kumaresh, 01/24/2008)
!     Updated for       
!******************************************************************************
!
      ! Reference to f90 modules
      USE COMODE_MOD,           ONLY : JLOP,   CSPEC,    CSPEC_FOR_KPP
      USE COMODE_MOD,           ONLY : IXSAVE, IYSAVE,   IZSAVE
      USE COMODE_MOD,           ONLY : R_KPP,  HSAVE_KPP
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
      
!======================================================================
! Variables used to pass the last/first step size b/w call 
!      
! For Rosenbrock:
!----------------
!      Nhexit=2, Nhnew = 3
! OUT    
!     RSTATUS(2)  -> Hexit, last accepted step before exit
!     RSTATUS(3)  -> Hnew, last predicted step (not yet taken)
!     For multiple restarts, use Hnew as Hstart 
!     in the subsequent run
! IN
!    RCNTRL(3)  -> Hstart, starting value for the integration step size
!
!      
! For LSODE:
!------------
! OUT
!     RSTATUS(1)  -> Texit, the time corresponding to the
!     computed Y upon return
!     RSTATUS(2)  -> Hexit, last predicted step before exit
!     For multiple restarts, use Hexit as Hstart in the following run
!
! IN      
!     RCNTRL(3)  -> Hstart, starting value for the integration step size
!
!
! For RADAU5:
!------------
! OUT
!     RSTATUS(1)  -> final time     
! IN      
!     RCNTRL(3)  -> not used
!
!
! For RUNGE-KUTTA
! ---------------
! OUT
!     same as Rosenbrock
!
!======================================================================
      
      
      ! Local variables
      REAL*8        :: T, TIN, TOUT
      INTEGER       :: ICNTRL(20)
      REAL(kind=dp) :: RCNTRL(20)
      INTEGER       :: ISTATUS(20)
      INTEGER       :: I, J, L, N, JJLOOP, IERR
      REAL(kind=dp) :: RSTATE(20)
      LOGICAL, SAVE :: FIRST = .TRUE.
      CHARACTER(LEN=255) :: ERR_MSG
      
      !=================================================================

      STEPMIN = 0.0d0
      STEPMAX = 0.0d0

      ! Suggested tolerance for Rosenbrock solvers
      DO i=1,NVAR
         RTOL(i) = 2.d-1
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

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( JJLOOP, I, J, L, N, RSTATE, ISTATUS )
!$OMP+FIRSTPRIVATE( RCNTRL, ICNTRL )
!$OMP+COPYIN( TIME )
!$OMP+SCHEDULE( DYNAMIC )
        DO JJLOOP = 1,NTT

           JLOOP = JJLOOP
           
           ! Get 3D coords from SMVGEAR's 1D coords
           I = IXSAVE(JJLOOP)
           J = IYSAVE(JJLOOP)
           L = IZSAVE(JJLOOP)

           ! Pass tracer concentrations from CSPEC_FOR_KPP to KPP working vector V_CSPEC
           DO N =1, NVAR
              V_CSPEC(N) = CSPEC_FOR_KPP(JLOOP,N)
           END DO

           ! Pass tracer concentrations from V_CSPEC to KPP working vector VAR.
           ! This also initializes the constant rate constants FIX.
           CALL Initialize()


           ! starting value for integration time step
           RCNTRL(3) = HSAVE_KPP(I,J,L)

           ! Recalculate rate constants
           CALL Update_RCONST()


           ! Integrate FWD -- phs --
           call integrate(TIN, TOUT, ICNTRL, RCNTRL, ISTATUS,
     $          RSTATE, IERR)
         
           ! Try another time if it failed
           IF ( IERR < 0 ) THEN
              write(6,*), ''
              write(6,100) istatus(1:8),rstate(1:3)
              write(6,*), ''
              write(ERR_MSG,'(a, i3)') 'Integrator error code :',IERR
              write(6,*),   'JLOOP, I, J, L ', JLOOP, I, J, L

              ! Reset first time step and start concentrations
              RCNTRL(3)  = 0d0
              CALL Initialize( ) ! v2.1 
              CALL Update_RCONST()
              call integrate(TIN, TOUT, ICNTRL, RCNTRL, ISTATUS,
     $                       RSTATE, IERR)
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
           DO N =1, NVAR
              CSPEC(JLOOP,N) = V_CSPEC(N)
           END DO

        ENDDO
!$OMP END PARALLEL DO

        ! Return to calling program
      END SUBROUTINE GCKPP_DRIVER

!------------------------------------------------------------------------------

      ! End of module
      END MODULE CHEMISTRY_MOD
