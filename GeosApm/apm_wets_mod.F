! $Id: wetscav_mod.f,v 1.30 2008/08/08 17:20:37 bmy Exp $
      MODULE APM_WETS_MOD
!
! Modified from GEOS-Chem wetscav_mod.f to consider wetscav for APM tracers 
! (GanLuo, 2010)
!
!******************************************************************************
!  Module APM_WETS_MOD contains arrays for used in the wet scavenging of
!  tracer in cloud updrafts, rainout, and washout. (bmy, 2/28/00, 3/5/08)
!
!  Module Variables:
!  ============================================================================
!  (1 ) NSOLMAX (INTEGER) : Max # of soluble tracers       [unitless]
!  (2 ) NSOL    (INTEGER) : Actual # of soluble tracers    [unitless]
!  (3 ) IDWETD  (INTEGER) : Index array for WETDEPBIN routine [unitless]
!  (4 ) Vud     (REAL*8 ) : Array for updraft velocity     [m/s]
!  (5 ) CLDLIQ  (REAL*8 ) : Array for cloud liquid water   [cm3 H2O/cm3 air]
!  (6 ) CLDICE  (REAL*8 ) : Array for cloud ice content    [cm3 ice/cm3 air]
!  (7 ) C_H2O   (REAL*8 ) : Array for Mixing ratio of ,      
!                            water, computed from Eice(T)  [v/v]
!  (8 ) PDOWN   (REAL*8 ) : Precip thru bottom of grid box [cm3 H2O/cm2 area/s]
!  (9 ) QQ      (REAL*8 ) : Rate of new precip formation   [cm3 H2O/cm3 air/s]
!  (10) EPSILON (REAL*8 ) : A very small positive number   [unitless]
!  (11) H2O2s   (REAL*8 ) : Array to save H2O2 for wetdep  [v/v]
!  (12) SO2s    (REAL*8 ) : Array to save SO2 for wetdep   [v/v]
!
!  Module Routines:
!  ============================================================================
!  (1 ) MAKE_QQ           : Constructs the QQ field (precipitable water)
!  (2 ) E_ICE             : Computes saturation vapor pressure for ice 
!  (3 ) COMPUTE_L2G       : Computes the ratio [v/v liquid] / [v/v gas] 
!  (4 ) COMPUTE_F         : Computes fraction of tracer lost in cloud updrafts
!  (5 ) F_AEROSOL         : Computes fraction of tracer scavenged in updrafts
!  (6 ) GET_ISOL          : Returns correct index for ND38 diagnostic
!  (7 ) RAINOUT           : Computes fraction of soluble tracer lost to rainout
!  (8 ) GET_RAINFRAC      : Computes rainout fraction -- called by RAINOUT
!  (9 ) WASHOUT           : Computes fraction of soluble tracer lost to washout
!  (10) WASHFRAC_AEROSOL  : Computes fraction of aerosol lost to washout
!  (11) WASHFRAC_LIQ_GAS  : Computes fraction of soluble gases lost to washout
!  (12) WETDEPBIN            : Driver routine for computing wet deposition losses
!  (13) LS_K_RAIN         : Computes K_RAIN (for LS precipitation)
!  (14) LS_F_PRIME        : Computes F_PRIME (for LS precipitation)
!  (15) CONV_F_PRIME      : Computes F_PRIME (for convective precipitation)
!  (16) SAFETY            : Stops WETDEPBIN w/ error msg if negative tracer found
!  (17) WETDEPBINID          : Initalizes the IDWETD array for routine WETDEPBIN
!  (18) GET_WETDEPBIN_NMAX   : Returns max # of soluble tracers per simulation
!  (19) GET_WETDEPBIN_NSOL   : Returns actual # of soluble tracers per simulation
!  (20) GET_WETDEPBIN_IDWETD : Returns CTM tracer # of for a given wetdep species 
!  (21) INIT_WETSCAVBIN      : Initializes fields used for computing wetdep losses
!  (22) CLEANUP_WETSCAVBIN   : Deallocates all allocatable module arrays
!
!  GEOS-CHEM modules referenced by wetscav_mod.f
!  ============================================================================
!  (1 ) dao_mod.f      : Module containing arrays for DAO met fields
!  (2 ) diag_mod.f     : Module containing GEOS-CHEM diagnostic arrays
!  (3 ) error_mod.f    : Module containing NaN and other error check routines
!  (4 ) logical_mod.f  : Module containing GEOS-CHEM logical switches
!  (5 ) pressure_mod.f : Module containing routines to compute P(I,J,L)
!  (6 ) tracer_mod.f   : Module containing GEOS-CHEM tracer array STT etc.
!  (7 ) tracerid_mod.f : Module containing pointers to tracers and emissions
!
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
!  NOTES:
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
!        WETDEPBIN.  Added GET_WETDEPBIN_NMAX function to return max # of soluble
!        tracers for allocating diagnostic arrays.  Added functions 
!        GET_WETDEPBIN_NSOL and GET_WETDEPBIN_IDWETD.  Now init H2O2s and SO2s
!        to the initial H2O2 and SO2 from STT.  Updated comments. 
!        (qli, bmy, 1/14/03)
!  (10) Improvements for SO2/SO4 scavenging (rjp, bmy, 3/23/03)
!  (11) Now references "time_mod.f".  Added driver routine DO_WETDEPBIN to
!        remove cumbersome calling sequence from MAIN program.  Also declared
!        WETDEPBIN and MAKE_QQ PRIVATE to this module. (bmy, 3/27/03)
!  (11) Add parallelization to routine WETDEPBIN (bmy, 3/17/04)
!  (12) Added carbon and dust aerosol tracers (rjp, tdf, bmy, 4/5/04)
!  (13) Added seasalt aerosol tracers (rjp, bec, bmy, 4/20/04)
!  (14) Added secondary organic aerosol tracers (rjp, bmy, 7/13/04)
!  (15) Now references "logical_mod.f" and "tracer_mod.f".  Now move all 
!        internal routines to the module and pass arguments explicitly in
!        order to facilitate parallelization on the Altix. (bmy, 7/20/04)
!  (16) Updated for mercury aerosol tracers (eck, bmy, 12/9/04)
!  (17) Updated for AS, AHS, LET, NH4aq, SO4aq.  Also now pass Hg2 wetdep loss
!        to "ocean_mercury_mod.f". (cas, sas, bmy, 1/20/05)
!  (18) Bug fix to avoid numerical blowup in WETDEPBIN.  Now use analytical
!        function for E_ICE(T). (bmy, 3/7/05)
!  (19) Added SO4s, NITs.  Increased NSOLMAX to 31.  Also block out 
!        parallel loop in WETDEPBIN for SGI MIPS compiler. (bec, bmy, 5/5/05)
!  (20) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (21) Bug fixes: do not over-deplete H2O2s.  Also include updates for
!        tagged Hg simulation. (dkh, rjp, eck, cdh, bmy, 1/6/06)
!  (22) Now wet deposit SOG4, SOA4. Remove unnecessary variables in WETDEPBIN.
!        (dkh, bmy, 5/18/06)
!  (23) Bug fixes in COMPUTE_F (bmy, 7/26/06)
!  (24) Resize DSTT array in WETDEPBIN to save memory.  Added fixes for GEOS-5
!        wet deposition per Hongyu Liu's suggestions. (bmy, 3/5/08)
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
!******************************************************************************
!
      IMPLICIT NONE

      !=================================================================
      ! MODULE PRIVATE DECLARATIONS
      !=================================================================

      ! Make everything PRIVATE ...
      PRIVATE

      ! ... and these routines 
      PUBLIC :: CLEANUP_WETSCAVBIN
      PUBLIC :: DO_WETDEPBIN
      PUBLIC :: GET_WETDEPBIN_IDWETD
      PUBLIC :: INIT_WETSCAVBIN
      PUBLIC :: WETDEPBINID

      !=================================================================
      ! MODULE VARIABLES
      !=================================================================

      ! Parameters
      INTEGER, PARAMETER   :: NSOLMAX = 94
      REAL*8,  PARAMETER   :: EPSILON = 1d-32

      ! Scalars
      INTEGER              :: NSOL 

      ! Arrays
      INTEGER              :: IDWETD(NSOLMAX)
      REAL*8,  ALLOCATABLE :: Vud(:,:)
      REAL*8,  ALLOCATABLE :: C_H2O(:,:,:)
      REAL*8,  ALLOCATABLE :: CLDLIQ(:,:,:)
      REAL*8,  ALLOCATABLE :: CLDICE(:,:,:)
      REAL*8,  ALLOCATABLE :: PDOWN(:,:,:)
      REAL*8,  ALLOCATABLE :: QQ(:,:,:)

      !=================================================================
      ! MODULE ROUTINES -- follow below the "CONTAINS" statement 
      !=================================================================
      CONTAINS

!------------------------------------------------------------------------------

      SUBROUTINE DO_WETDEPBIN( State_Met, State_Chm )
!
!******************************************************************************
!  Subroutine DO_WETDEPBIN is a driver for the wet deposition code, called
!  from the MAIN program. (bmy, 3/27/03, 3/5/08)
!
!  NOTES:
!  (1 ) Now references LPRT from "logical_mod.f" (bmy, 7/20/04)
!  (2 ) Don't do rainout/washout for conv precip for GEOS-5 (hyl, bmy, 3/5/08)
!  17 Oct 2013 - M. Sulprizio- Now pass State_Chm as argument
!******************************************************************************
!
      ! References to F90 modules
      USE ERROR_MOD,          ONLY : DEBUG_MSG
      USE GIGC_State_Met_Mod, ONLY : MetState
      USE GIGC_State_Chm_Mod, ONLY : ChmState
      USE LOGICAL_MOD,        ONLY : LPRT

      USE CMN_SIZE_MOD             ! Size parameters
!
! !INPUT PARAMETERS: 
!
      TYPE(MetState), INTENT(IN)    :: State_Met   ! Meteorology State object
!
! !INPUT/OUTPUT PARAMETERS:
!
      TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object

      !==================================================================
      ! DO_WETDEPBIN begins here!
      !==================================================================

      ! Wetdep by large-scale (stratiform) precip
      CALL MAKE_QQ( State_Met, .TRUE. )
      IF ( LPRT ) CALL DEBUG_MSG( '### DO_WETDEPBIN: before LS wetdep' )
      CALL WETDEPBIN( State_Met, State_Chm, .TRUE. )
      IF ( LPRT ) CALL DEBUG_MSG( '### DO_WETDEPBIN: after LS wetdep' )

#if   !defined( GEOS_5 )

      !------------------------------------------------------------------
      ! NOTE FROM HONGYU LIU (hyl@nianet.org) -- 3/5/08
      !
      ! Rainout and washout from convective precipitation for previous
      ! GEOS archives were intended to represent precipitation from 
      ! cloud anvils [Liu et al., 2001]. For GEOS-5 (as archived at 
      ! Harvard), the cloud anvil precipitation was already included 
      ! in the large-scale precipitation. 
      !
      ! Therefore, we insert a #if block to ensure that call MAKE_QQ
      ! and WETDEPBIN are not called for convective precip in GEOS-5.
      ! (hyl, bmy, 3/5/08)
      !------------------------------------------------------------------

      ! Wetdep by convective precip
      CALL MAKE_QQ( State_Met, .FALSE. )
      IF ( LPRT ) CALL
     &DEBUG_MSG( '### DO_WETDEPBIN: before conv wetdep' )
      CALL WETDEPBIN( State_Met, State_Chm, .FALSE. )
      IF ( LPRT ) CALL
     &DEBUG_MSG( '### DO_WETDEPBIN: after conv wetdep' )

#endif

      ! Return to calling program
      END SUBROUTINE DO_WETDEPBIN

!------------------------------------------------------------------------------

      SUBROUTINE MAKE_QQ( State_Met, LS )
!
!*****************************************************************************
!  Subroutine MAKE_QQ computes the large-scale or convective precipitation
!  fields for use with wetdep.f. (hyl, bmy, 2/29/00, 11/8/02)
!
!  References (see above for full citations):
!  ===========================================================================
!  (1 ) Liu et al, 2000
!  (2 ) Jacob et al, 2000
!
!  NOTES:
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
!  09 Nov 2012 - M. Payer    - Replaced all met field arrays with State_Met
!                              derived type object
!******************************************************************************
! 
      ! References to F90 modules
      USE ERROR_MOD,          ONLY : ALLOC_ERR
      USE GIGC_State_Met_Mod, ONLY : MetState

      USE CMN_SIZE_MOD             ! Size parameters

      ! Arguments
      LOGICAL,        INTENT(IN)  :: LS          ! = T for Large-scale precip
      TYPE(MetState), INTENT(IN)  :: State_Met   ! Meteorology State object

      ! Local variables
      INTEGER              :: I, J, L, AS
      REAL*8               :: PTEMP, FRAC
      LOGICAL              :: FIRST = .TRUE.

      !=================================================================
      ! MAKE_QQ begins here!
      !=================================================================
      IF ( FIRST ) THEN

         ! Allocate PDOWN on first call
         ALLOCATE( PDOWN( LLPAR, IIPAR, JJPAR ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'PDOWN' )
         PDOWN = 0d0
      
         ! Allocate QQ on first call
         ALLOCATE( QQ( LLPAR, IIPAR, JJPAR ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'QQ' )
         QQ = 0d0
         
         ! Reset flag
         FIRST = .FALSE.
      ENDIF

      !=================================================================
      ! Loop over surface grid boxes
      !=================================================================
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, FRAC, L, PTEMP )
!$OMP+SCHEDULE( DYNAMIC )
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         !==============================================================
         ! If there is total precipitation in the (I,J) column, then:
         ! 
         ! (1) Compute FRAC, the large scale fraction (if LS = .TRUE.) 
         !     or convective fraction (if LS = .FALSE.) total 
         !     precipitation.  FRAC is computed from PRECTOT and PRECCON.
         !
         ! (2) Compute QQ, the rate of formation of precipitation 
         !     [cm3 H2O/cm3 air/s].  From MOISTQ [kg H2O/kg air/s], 
         !     the unit conversion is: 
         !
         !     kg H2O   |   m^3 H2O   | AIRDEN kg air         m^3 H2O
         !  ------------+-------------+--------------- ==> -------------   
         !   kg air * s | 1000 kg H2O |    m^3 air          m^3 air * s
         !
         ! and
         !
         !         m^3 H2O                         cm^3 H2O
         !      -------------  is equivalent to  -------------- 
         !       m^3 air * s                      cm^3 air * s!
         !       
         ! since the same conversion factor (10^6 cm^3/m^3) is in both
         ! the numerator and the denominator.
         !
         ! Therefore, the equation for QQ is:
         !
         !   QQ(L,I,J) = FRAC * MOISTQ(L,I,J) * AIRDEN(L,I,J) / 1000.0
         !     
         ! (3) Compute PDOWN, the column precipitation 
         !     [cm3 H2O/cm2 air/s], by multiplying QQ(L,I,J) by 
         !     BXHEIGHT(I,J,L) * 100 cm.  
         !
         ! (4) The reason why we do not force PTEMP to be positive is 
         !     that PRECTOT is the integral of the MOISTQ field.  MOISTQ 
         !     contains both negative (evap) and positive (precip) 
         !     values.  If we forced PTEMP to be positive, then we would
         !     be adding extra precipitation to PDOWN (hyl, bmy, 3/6/99).
         !==============================================================
         IF ( State_Met%PRECTOT(I,J) > 0d0 ) THEN

            ! Large scale or convective fraction of precipitation
            IF ( LS ) THEN
               FRAC = ( State_Met%PRECTOT(I,J) - State_Met%PRECCON(I,J))
     &                / State_Met%PRECTOT(I,J) 
            ELSE
               FRAC =   State_Met%PRECCON(I,J) / State_Met%PRECTOT(I,J)
            ENDIF

            ! Start at the top of the atmosphere
            L = LLPAR

            ! Compute QQ and PDOWN.  Keep PTEMP for the next level
            QQ(L,I,J)    = FRAC * State_Met%MOISTQ(L,I,J) *
     &                     State_Met%AIRDEN(L,I,J) / 1d3
            PTEMP        = QQ(L,I,J) * State_Met%BXHEIGHT(I,J,L) * 1d2
            PDOWN(L,I,J) = PTEMP

            ! PDOWN cannot be negative
            IF ( PDOWN(L,I,J) < 0d0 ) PDOWN(L,I,J) = 0.d0

            ! Loop down from LLPAR to the surface
            DO L = LLPAR-1, 1, -1
               
               ! Compute QQ and PDOWN.  Keep PTEMP for the next level.
               QQ(L,I,J)    = FRAC  * State_Met%MOISTQ(L,I,J) *
     &                        State_Met%AIRDEN(L,I,J) / 1d3
               PDOWN(L,I,J) = PTEMP + QQ(L,I,J) *
     &                        State_Met%BXHEIGHT(I,J,L) * 1d2  
               PTEMP        = PDOWN(L,I,J)

               ! PDOWN cannot be negative
               IF ( PDOWN(L,I,J) < 0.0d0 ) PDOWN(L,I,J) = 0.d0
            ENDDO
  
         !==============================================================
         ! If there is no precipitation reaching the surface in the 
         ! (I,J) column, then assume any precipitation at altitude to 
         ! be large-scale.
         ! 
         ! (1) Assume the large scale fraction = 1d0, 
         !                convective fraction  = 0d0
         ! (2) Compute QQ as described above
         ! (3) Compute PDOWN as described above
         !==============================================================
         ELSE

            ! Assume large-scale precipitation!
            IF ( LS ) THEN
               FRAC = 1d0
            ELSE         
               FRAC = 0d0
            ENDIF

            ! Start at the top of the atmosphere
            L = LLPAR

            ! Compute QQ and PDOWN.  Keep PTEMP for the next level
            QQ(L,I,J)    = FRAC * State_Met%MOISTQ(L,I,J) *
     &                     State_Met%AIRDEN(L,I,J) / 1d3
            PTEMP        = QQ(L,I,J) * State_Met%BXHEIGHT(I,J,L) * 1d2
            PDOWN(L,I,J) = PTEMP
           
            ! PDOWN cannot be negative
            IF( PDOWN(L,I,J) < 0d0 ) PDOWN(L,I,J) = 0.d0

            ! Loop down from LLPAR to the surface
            DO L = LLPAR-1, 1, -1
              
               ! Compute QQ and PDOWN.  Keep PTEMP for the next level
               QQ(L,I,J)    = FRAC  * State_Met%MOISTQ(L,I,J) *
     &                        State_Met%AIRDEN(L,I,J) / 1d3
               PDOWN(L,I,J) = PTEMP + QQ(L,I,J) *
     &                        State_Met%BXHEIGHT(I,J,L) * 1d2 
               PTEMP        = PDOWN(L,I,J)

               ! PDOWN cannot be negative
               IF ( PDOWN(L,I,J) < 0.0d0 ) PDOWN(L,I,J) = 0.d0
            ENDDO
         ENDIF
      ENDDO  ! J
      ENDDO  ! I
!$OMP END PARALLEL DO

      ! Return to calling program
      END SUBROUTINE MAKE_QQ

!------------------------------------------------------------------------------

      FUNCTION E_ICE( TK ) RESULT( VALUE )
! 
!******************************************************************************
!  Subroutine E_ICE computes Eice(T), the saturation vapor pressure of ice
!  at a given Celsius temperature. (bmy, 2/8/05)
!  
!  Arguments as Input:
!  ============================================================================
!  (1 ) TK (REAL*8) : Ambient temperature [K] 
! 
!  References:
!  ============================================================================
!  (1 ) Marti & Mauersberber (GRL '93) formulation of saturation 
!        vapor pressure of ice [Pa] is: log P = A/TK + B
!
!  NOTES:
!  (1 ) Now use the same analytic function as the Goddard CTM (bmy, 2/8/05)
!******************************************************************************
!
      ! Arguments as Input
      REAL*8, INTENT(IN) :: TK

      ! Return value
      REAL*8             :: VALUE

      ! Parameters
      REAL*8, PARAMETER  :: A = -2663.5d0
      REAL*8, PARAMETER  :: B =  12.537d0

      !=================================================================
      ! E_ICE begins here!
      !=================================================================
      
      ! Saturation vap press of Ice [Pa] -- divide by 100 for [hPa]
      VALUE = ( 10d0**( A/TK + B ) ) / 100d0 

      ! Return to calling program
      END FUNCTION E_ICE

!------------------------------------------------------------------------------

      SUBROUTINE COMPUTE_L2G( Kstar298, H298_R, TK, H2OLIQ, L2G )
!
!******************************************************************************
!  Subroutine COMPUTE_L2G computes the ratio L2G = Cliq / Cgas, which is 
!  the mixing ratio of tracer in the liquid phase, divided by the mixing 
!  ratio of tracer in the gas phase.  (bmy, 2/23/00, 11/8/02)
!
!  The ratio Cliq / Cgas is obtained via Henry's law.  The appropriate 
!  values of Kstar298 and H298_R must be supplied for each tracer.  
!  (cf Jacob et al 2000, p. 3)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) Kstar298 (REAL*8) : Eff. Henry's law constant @ 298 K   [moles/atm]
!  (2 ) H298_R   (REAL*8) : Molar heat of formation @ 298 K / R [K]
!  (3 ) TK       (REAL*8) : Temperature at grid box (I,J,L)     [K]
!  (4 ) H2OLIQ   (REAL*8) : Liquid water content at (I,J,L)     [cm3 H2O/cm3 air]
!
!  Arguments as Output:
!  ============================================================================
!  (5 ) L2G      (REAL*8) : Cliq/Cgas ratio for given tracer  [unitless]
!
!  References (see above for full citations):
!  ===========================================================================
!  (1 ) Jacob et al, 2000
!
!  NOTES:
!  (1 ) Bundled into "wetscav_mod.f" (bmy, 11/8/02)
!******************************************************************************
!
      ! Arguments
      REAL*8, INTENT(IN)  :: KStar298, H298_R, TK, H2OLIQ
      REAL*8, INTENT(OUT) :: L2G
      
      ! Local variables
      REAL*8              :: Kstar

      ! R = universal gas constant [atm/moles/K]
      REAL*8, PARAMETER   :: R = 8.32d-2

      ! INV_T0 = 1/298 K
      REAL*8, PARAMETER   :: INV_T0 = 1d0 / 298d0

      !=================================================================
      ! COMPUTE_L2G begins here!
      !=================================================================

      ! Get Kstar, the effective Henry's law constant for temperature TK
      Kstar = Kstar298 * EXP( -H298_R * ( ( 1d0 / TK ) - INV_T0 ) )

      ! Use Henry's Law to get the ratio:
      ! [ mixing ratio in liquid phase / mixing ratio in gas phase ]
      L2G   = Kstar * H2OLIQ * R * TK

      ! Return to calling program
      END SUBROUTINE COMPUTE_L2G

!------------------------------------------------------------------------------

      SUBROUTINE RAINOUT( I, J, L, N, K_RAIN, DT, F, RAINFRAC,
     &                    State_Met )
!
!******************************************************************************
!  Subroutine RAINOUT computes RAINFRAC, the fraction of soluble tracer
!  lost to rainout events in precipitation. (djj, bmy, 2/28/00, 3/5/08)
!
!  Arguments as Input:
!  ============================================================================
!  (1-3) I, J, L  (INTEGER) : Grid box lon-lat-alt indices
!  (4  ) N        (INTEGER) : Tracer number
!  (5  ) K_RAIN   (REAL*8 ) : Rainout rate constant for tracer N [s^-1]
!  (6  ) DT       (REAL*8 ) : Timestep for rainout event         [s]
!  (7  ) F        (REAL*8 ) : Fraction of grid box precipitating [unitless]
!
!  Arguments as Output:
!  ============================================================================
!  (8  ) RAINFRAC (REAL*8)  : Fraction of tracer lost to rainout [unitless]
!
!  References (see above for full citations):
!  ============================================================================
!  (1 ) Jacob et al, 2000
!  (2 ) Chin et al, 1996
!
!  NOTES:
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
!  09 Nov 2012 - M. Payer    - Replaced all met field arrays with State_Met
!                              derived type object
!******************************************************************************
!
      ! References to F90 modules
      USE APM_INIT_MOD,       ONLY : NGCOND,    NSO4
      USE APM_INIT_MOD,       ONLY : NSEA,      NDSTB
      USE APM_INIT_MOD,       ONLY : NCTSO4,    NCTBCOC
      USE APM_INIT_MOD,       ONLY : NCTDST,    NCTSEA
      USE APM_INIT_MOD,       ONLY : NBCPIF,    NBCPIO
      USE APM_INIT_MOD,       ONLY : NOCPIF,    NOCPIO
      USE APM_INIT_MOD,       ONLY : NBCPOF,    NBCPOO
      USE APM_INIT_MOD,       ONLY : NOCPOF,    NOCPOO
      USE APM_INIT_MOD,       ONLY : IDTSO4G,   IDTSO4BIN1, IDTDSTBIN1
      USE APM_INIT_MOD,       ONLY : IDTCTSO4,  IDTCTBCOC,  IDTCTDST
      USE APM_INIT_MOD,       ONLY : IDTCTSEA,  IDTSEABIN1
      USE APM_INIT_MOD,       ONLY : IDTBCPIFF, IDTBCPIBB
      USE APM_INIT_MOD,       ONLY : IDTOCPIFF, IDTOCPIBB
      USE APM_INIT_MOD,       ONLY : IDTBCPOFF, IDTBCPOBB
      USE APM_INIT_MOD,       ONLY : IDTOCPOFF, IDTOCPOBB
      USE GIGC_State_Met_Mod, ONLY : MetState
      USE ERROR_MOD,          ONLY : ERROR_STOP

      USE CMN_SIZE_MOD             ! Size parameters

      IMPLICIT NONE

      ! Arguments
      INTEGER,        INTENT(IN)  :: I, J, L, N
      REAL*8,         INTENT(IN)  :: K_RAIN, DT, F
      REAL*8,         INTENT(OUT) :: RAINFRAC
      TYPE(MetState), INTENT(IN)  :: State_Met   ! Meteorology State object

      ! Local variables 
      REAL*8               :: L2G, I2G, C_TOT, F_L, F_I, K, TK, SO2LOSS

      ! CONV = 0.6 * SQRT( 1.9 ), used for the ice to gas ratio for H2O2
      REAL*8, PARAMETER    :: CONV = 8.27042925126d-1

      !==================================================================
      ! RAINOUT begins here!
      !
      ! For aerosols, set K = K_RAIN and compute RAINFRAC according
      ! to Eq. 10 of Jacob et al 2000.  Call function GET_RAINFRAC.
      !==================================================================

      ! Save the local temperature in TK for convenience
      TK = State_Met%T(I,J,L)

#if   defined( GEOS_5 )
      !------------------------------------------------------------------
      ! NOTE FROM HONGYU LIU (hyl@nianet.org) -- 3/5/08
      !
      ! Lead-210 (210Pb) and Beryllium-7 (7Be) simulations indicate 
      ! that we can improve the GEOS-5 simulation by (1) turning off
      ! rainout/washout for convective precip (see DO_WETDEPBIN) 
      ! and (2) suppressing rainout for large-scale precip at  
      ! temperatures below 258K.
      !
      ! Place an #if block here to set RAINFRAC=0 when T < 258K for 
      ! GEOS-5 met.  This will suppress rainout. (hyl, bmy, 3/5/08)
      !-------------------------------------------------------------------   
      IF ( TK < 258d0 ) THEN
         RAINFRAC = 0d0
         RETURN
      ENDIF
#endif

      !------------------------------
      ! SO4G (aerosol)
      !------------------------------
      IF ( N >=IDTSO4G .and. (N<(IDTSO4G+NGCOND)) ) THEN
         RAINFRAC = GET_RAINFRAC( K_RAIN, F, DT )

      !------------------------------
      ! SO2
      !------------------------------
      ! SO4 and SO4aq (aerosol)
      !----------------------------
      ELSE IF ( N >=IDTSO4BIN1 .and. (N<(IDTSO4BIN1+NSO4)) ) THEN
         RAINFRAC = GET_RAINFRAC( K_RAIN, F, DT )

      !-------------------------------
      ! Sulfate and SOM (aerosol)
      !-------------------------------
      ELSE IF ( N >= IDTCTSO4 .and. N<(IDTCTSO4+NCTSO4) ) THEN
         RAINFRAC = GET_RAINFRAC( K_RAIN, F, DT )

      !------------------------------
      ! BC HYDROPHILIC (aerosol) or
      ! OC HYDROPHILIC (aerosol)
      !------------------------------
      ELSE IF ( (N >= IDTCTBCOC .and. N<(IDTCTBCOC+NCTBCOC))
     &          .or. N == IDTBCPIFF .or. N == IDTBCPIBB
     &          .or. N == IDTOCPIFF .or. N == IDTOCPIBB ) THEN
         RAINFRAC = GET_RAINFRAC( K_RAIN, F, DT )

      !-------------------------------
      ! BC HYDROPHOBIC (aerosol) or
      ! OC HYDROPHOBIC (aerosol)
      !-------------------------------
      ELSE IF ( N == IDTBCPOFF .or. N == IDTBCPOBB
     &          .or. N == IDTOCPOFF .or. N == IDTOCPOBB ) THEN

         ! No rainout 
         RAINFRAC = 0.0D0                  

      !-------------------------------
      ! DUST all size bins (aerosol)
      !-------------------------------
      ELSE IF ( (N >= IDTCTDST .and. N<(IDTCTDST+NCTDST)) .or.
     &          (N >= IDTDSTBIN1 .and. N<(IDTDSTBIN1+NDSTB))  ) THEN
         RAINFRAC = GET_RAINFRAC( K_RAIN, F, DT )

      !------------------------------
      ! Accum  seasalt (aerosol) or
      ! Coarse seasalt (aerosol)
      !------------------------------
      ELSE IF ( (N >= IDTCTSEA .and. N<(IDTCTSEA+NCTSEA)) .or.
     &          (N >= IDTSEABIN1 .and. N<(IDTSEABIN1+NSEA)) ) THEN
         RAINFRAC = GET_RAINFRAC( K_RAIN, F, DT )      

      !------------------------------
      ! ERROR: insoluble tracer!
      !------------------------------
      ELSE
         CALL ERROR_STOP( 'Invalid tracer!', 'RAINOUT (wetscav_mod.f)' )

      ENDIF
      
      ! Return to calling program
      END SUBROUTINE RAINOUT

!------------------------------------------------------------------------------

      FUNCTION GET_RAINFRAC( K, F, DT ) RESULT( RAINFRAC )
!
!******************************************************************************
!  Function GET_RAINFRAC computes the fraction of tracer lost to rainout 
!  according to Jacob et al 2000. (bmy, 11/8/02, 7/20/04)
!
!  Arguments as Input:
!  =========================================================================== 
!  (1 ) K  (REAL*8) : Rainout rate constant              [1/s]
!  (2 ) DT (REAL*8) : Timestep for rainout event         [s]
!  (3 ) F  (REAL*8) : Fraction of grid box precipitating [unitless]
!
!  NOTES:
!  (1 ) Now move internal routines GET_RAINFRAC to the module and pass all 
!        arguments explicitly.  This facilitates parallelization on the 
!        Altix platform (bmy, 7/20/04) 
!******************************************************************************
!
      ! Arguments
      REAL*8, INTENT(IN) :: K, F, DT
      
      ! Local variables
      REAL*8             :: RAINFRAC

      !=================================================================
      ! GET_RAINFRAC begins here!
      !=================================================================

      ! (Eq. 10, Jacob et al, 2000 ) 
      RAINFRAC = F * ( 1 - EXP( -K * DT ) )

      ! Return to RAINOUT
      END FUNCTION GET_RAINFRAC

!------------------------------------------------------------------------------

      SUBROUTINE WASHOUT( I, J, L, N, PP, DT, F, WASHFRAC, AER,
     &                    State_Met )
!
!******************************************************************************
!  Subroutine WASHOUT computes WASHFRAC, the fraction of soluble tracer
!  lost to washout events in precipitation. (djj, bmy, 2/28/00, 5/18/06)
!
!  Arguments as Input:
!  ============================================================================
!  (1-3) I, J, L  (INTEGER) : Grid box lon-lat-alt indices [unitless]
!  (4  ) N        (INTEGER) : Tracer number                [unitless]
!  (5  ) PP       (REAL*8 ) : Precip rate thru at bottom    
!                             of grid box (I,J,L)          [cm3 H2O/cm2 air/s]
!  (6  ) DT       (REAL*8 ) : Timestep for rainout event   [s]
!  (7  ) F        (REAL*8 ) : Fraction of grid box 
!                             that is precipitating        [unitless]
!
!  Arguments as Output:
!  ============================================================================
!  (8  ) WASHFRAC (REAL*8)  : Fraction of tracer lost to rainout [unitless]
!  (9  ) AER      (LOGICAL) : = T if the tracer is an aerosol, =F otherwise
!
!  Reference (see above for full citations):
!  ============================================================================
!  (1  ) Jacob et al, 2000
!
!  NOTES:
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
!******************************************************************************
!
      ! References to F90 modules
      USE APM_INIT_MOD,       ONLY : NGCOND,     NSO4
      USE APM_INIT_MOD,       ONLY : NSEA,       NDSTB
      USE APM_INIT_MOD,       ONLY : NCTSO4,     NCTBCOC
      USE APM_INIT_MOD,       ONLY : NCTDST,     NCTSEA
      USE APM_INIT_MOD,       ONLY : NBCPIF,     NBCPIO
      USE APM_INIT_MOD,       ONLY : NOCPIF,     NOCPIO
      USE APM_INIT_MOD,       ONLY : NBCPOF,     NBCPOO
      USE APM_INIT_MOD,       ONLY : NOCPOF,     NOCPOO
      USE APM_INIT_MOD,       ONLY : IDTSO4G,    IDTSO4BIN1
      USE APM_INIT_MOD,       ONLY : IDTCTBCOC,  IDTCTDST
      USE APM_INIT_MOD,       ONLY : IDTCTSO4,   IDTCTSEA
      USE APM_INIT_MOD,       ONLY : IDTSEABIN1, IDTDSTBIN1
      USE APM_INIT_MOD,       ONLY : IDTBCPIFF,  IDTBCPIBB
      USE APM_INIT_MOD,       ONLY : IDTOCPIFF,  IDTOCPIBB
      USE APM_INIT_MOD,       ONLY : IDTBCPOFF,  IDTBCPOBB
      USE APM_INIT_MOD,       ONLY : IDTOCPOFF,  IDTOCPOBB
      USE ERROR_MOD,          ONLY : ERROR_STOP
      USE GIGC_State_Met_Mod, ONLY : MetState

      USE CMN_SIZE_MOD             ! Size parameters

      ! Arguments
      INTEGER,        INTENT(IN)  :: I, J, L, N
      REAL*8,         INTENT(IN)  :: PP, DT, F
      REAL*8,         INTENT(OUT) :: WASHFRAC
      LOGICAL,        INTENT(OUT) :: AER
      TYPE(MetState), INTENT(IN)  :: State_Met  ! Meteorology State object

      ! Local variables 
      REAL*8               :: L2G, DZ, TK, SO2LOSS

      ! First order washout rate constant for HNO3, aerosols = 1 cm^-1
      REAL*8, PARAMETER    :: K_WASH = 1d0

      !=================================================================
      ! WASHOUT begins here!
      !
      ! Call either WASHFRAC_AEROSOL or WASHFRAC_LIQ_GAS to compute the
      ! fraction of tracer lost to washout according to Jacob et al 2000
      !=================================================================

      ! TK is Kelvin temperature 
      TK = State_Met%T(I,J,L)

      ! DZ is the height of the grid box in cm
      DZ = State_Met%BXHEIGHT(I,J,L) * 1d2

      !------------------------------
      ! HNO3 (aerosol)
      !------------------------------
      IF ( N >= IDTSO4G .and. N<(IDTSO4G+NGCOND) ) THEN
         AER      = .TRUE.
         WASHFRAC = WASHFRAC_AEROSOL( DT, F, K_WASH, PP, TK )

      !------------------------------
      ! SO2 (aerosol treatment)
      !------------------------------
      ! SO4 and SO4aq (aerosol)
      !------------------------------
      ELSE IF ( N >= IDTSO4BIN1 .and. N < (IDTSO4BIN1+NSO4) ) THEN
         AER      = .TRUE.
         WASHFRAC = WASHFRAC_AEROSOL( DT, F, K_WASH, PP, TK )

      !------------------------------
      ! Sulfate and SOM (aerosol)
      !------------------------------
      ELSE IF ( N >= IDTCTSO4 .and. N<(IDTCTSO4+NCTSO4) ) THEN
         AER      = .TRUE.
         WASHFRAC = WASHFRAC_AEROSOL( DT, F, K_WASH, PP, TK )

      !------------------------------
      ! BC HYDROPHILIC (aerosol) or
      ! OC HYDROPHILIC (aerosol) or
      ! BC HYDROPHOBIC (aerosol) or
      ! OC HYDROPHOBIC (aerosol) 
      !------------------------------
      ELSE IF ( (N >= IDTCTBCOC .and. N<(IDTCTBCOC+NCTBCOC))
     &          .or. N == IDTBCPIFF .or. N == IDTBCPIBB
     &          .or. N == IDTOCPIFF .or. N == IDTOCPIBB
     &          .or. N == IDTBCPOFF .or. N == IDTBCPOBB
     &          .or. N == IDTOCPOFF .or. N == IDTOCPOBB ) THEN
         AER      = .TRUE.
         WASHFRAC = WASHFRAC_AEROSOL( DT, F, K_WASH, PP, TK )

      !------------------------------
      ! DUST all size bins (aerosol)
      !------------------------------
      ELSE IF ( (N >= IDTCTDST .and. N<(IDTCTDST+NCTDST)) .or.
     &          (N >= IDTDSTBIN1 .and. N<(IDTDSTBIN1+NDSTB))  ) THEN
         AER      = .TRUE.
         WASHFRAC = WASHFRAC_AEROSOL( DT, F, K_WASH, PP, TK )

      !------------------------------
      ! Accum  seasalt (aerosol) or
      ! Coarse seasalt (aerosol)
      !------------------------------
      ELSE IF ( (N >= IDTCTSEA .and. N<(IDTCTSEA+NCTSEA)) .or.
     &          (N >= IDTSEABIN1 .and. N<(IDTSEABIN1+NSEA)) ) THEN
         AER      = .TRUE.
         WASHFRAC = WASHFRAC_AEROSOL( DT, F, K_WASH, PP, TK )

      !------------------------------
      ! ERROR: Insoluble tracer
      !------------------------------
      ELSE 
         CALL ERROR_STOP('Invalid tracer!','WASHOUT(wetscavbin_mod.f)')

      ENDIF

      ! Return to calling program
      END SUBROUTINE WASHOUT

!------------------------------------------------------------------------------

      FUNCTION WASHFRAC_AEROSOL( DT, F, K_WASH, PP, TK ) 
     &         RESULT( WASHFRAC )
!
!******************************************************************************
!  Function WASHFRAC_AEROSOL returns the fraction of soluble aerosol tracer 
!  lost to washout.  (bmy, 11/8/02, 7/20/04)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) TK       (REAL*8 ) : Temperature in grid box            [K]
!  (2 ) F        (REAL*8 ) : Fraction of grid box 
!                             that is precipitating             [unitless]
!  (3 ) K_WASH   (REAL*8 ) : 1st order washout rate constant    [1/cm]
!  (3 ) PP       (REAL*8 ) : Precip rate thru at bottom    
!                             of grid box (I,J,L)           [cm3 H2O/cm2 air/s]
!
!  NOTES:
!  (1 ) WASHFRAC_AEROSOL used to be an internal function to subroutine WASHOUT.
!        This caused NaN's in the parallel loop on Altix, so we moved it to
!        the module and now pass Iall arguments explicitly (bmy, 7/20/04)
!******************************************************************************
!   
      ! Arguments
      REAL*8, INTENT(IN) :: DT, F, K_WASH, PP, TK

      ! Function value
      REAL*8             :: WASHFRAC

      !=================================================================
      ! WASHFRAC_AEROSOL begins here!
      !=================================================================

      ! Washout only happens at or above 268 K
      IF ( TK >= 268d0 ) THEN
         WASHFRAC = F * ( 1d0 - EXP( -K_WASH * ( PP / F ) * DT ) )
      ELSE
         WASHFRAC = 0d0
      ENDIF

      ! Return to calling program
      END FUNCTION WASHFRAC_AEROSOL

!------------------------------------------------------------------------------

      FUNCTION WASHFRAC_LIQ_GAS( Kstar298, H298_R, PP, DT, 
     &                           F,        DZ,     TK, K_WASH ) 
     &         RESULT( WASHFRAC )
!
!******************************************************************************
!  Function WASHFRAC_LIQ_GAS returns the fraction of soluble liquid/gas phase 
!  tracer lost to washout. (bmy, 11/8/02, 7/20/04)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) Kstar298 (REAL*8 ) : Eff. Henry's law constant @ 298 K  [moles/atm]
!  (2 ) H298_R   (REAL*8 ) : Henry's law coefficient            [K]
!  (3 ) PP       (REAL*8 ) : Precip rate thru at bottom    
!                             of grid box (I,J,L)           [cm3 H2O/cm2 air/s]
!  (4 ) DT       (REAL*8 ) : Dynamic timestep                   [s]
!  (5 ) F        (REAL*8 ) : Fraction of grid box 
!                             that is precipitating             [unitless]
!  (6 ) DZ       (REAL*8 ) : Height of grid box                 [cm]
!  (7 ) TK       (REAL*8 ) : Temperature in grid box            [K]
!  (8 ) K_WASH   (REAL*8 ) : 1st order washout rate constant    [1/cm]
!
!  NOTES:
!  (1 ) WASHFRAC_LIQ_GAS used to be an internal function to subroutine WASHOUT.
!        This caused NaN's in the parallel loop on Altix, so we moved it to
!        the module and now pass all arguments explicitly (bmy, 7/20/04)
!******************************************************************************
!
      ! Arguments
      REAL*8, INTENT(IN) :: Kstar298, H298_R, PP, DT, F, DZ, TK, K_WASH

      ! Local variables
      REAL*8             :: L2G, LP, WASHFRAC, WASHFRAC_F_14

      !=================================================================
      ! WASHFRAC_LIQ_GAS begins here!
      !=================================================================

      ! Suppress washout below 268 K
      IF ( TK >= 268d0 ) THEN

         ! Rainwater content in the grid box (Eq. 17, Jacob et al, 2000)
         LP = ( PP * DT ) / ( F * DZ ) 

         ! Compute liquid to gas ratio for H2O2, using the appropriate 
         ! parameters for Henry's law -- also use rainwater content Lp
         ! (Eqs. 7, 8, and Table 1, Jacob et al, 2000)
         CALL COMPUTE_L2G( Kstar298, H298_R, TK, LP, L2G )

         ! Washout fraction from Henry's law (Eq. 16, Jacob et al, 2000)
         WASHFRAC = L2G / ( 1d0 + L2G )

         ! Washout fraction / F from Eq. 14, Jacob et al, 2000
         WASHFRAC_F_14 = 1d0 - EXP( -K_WASH * ( PP / F ) * DT )

         ! Do not let the Henry's law washout fraction exceed
         ! ( washout fraction / F ) from Eq. 14 -- this is a cap
         IF ( WASHFRAC > WASHFRAC_F_14 ) WASHFRAC = WASHFRAC_F_14
            
      ELSE
         WASHFRAC = 0d0
            
      ENDIF

      ! Return to calling program
      END FUNCTION WASHFRAC_LIQ_GAS

!------------------------------------------------------------------------------

      SUBROUTINE WETDEPBIN( State_Met, State_Chm, LS )
!
!******************************************************************************
!  Subroutine WETDEPBIN computes the downward mass flux of tracer due to washout 
!  and rainout of aerosols and soluble tracers in a column.  The timestep is 
!  the dynamic timestep. (hyl, bey, bmy, djj, 4/2/99, 5/24/06)
!
!  The precip fields through the bottom of each level are indexed as follows:
!
!       Layer          GISS-CTM II         GEOS-CTM
!
!      ------------------------------------------------- Top of Atm.
!        LM            PSSW4(I,J,LM-1)   PDOWN(LM,I,J)
!                          |                  |
!      ====================V==================V========= Max Extent 
!        LM-1          PSSW4(I,J,LM)     PDOWN(LM-1,I,J)  of Clouds
!                          |                  |
!      --------------------V------------------V---------
!                         ...                ...             
!
!      -------------------------------------------------
!        4             PSSW4(I,J,3)      PDOWN(4,I,J)
!                          |                  |
!      --------------------V------------------V----------
!        3             PSSW4(I,J,2)      PDOWN(3,I,J)
!                          |                  |
!      --------------------V------------------V--------- Cloud base
!        2             PSSW4(I,J,1)      PDOWN(2,I,J) 
!                          |                  |
!      -  -  -  -  -  -  - V -  -   -   -   - V -  -  - 
!        1                               PDOWN(1,I,J) 
!                                             |
!      =======================================V========= Ground
!
!  From the diagram, we have the following for layer L:
!     
!  GISS-CTM:
!  (a) Precip coming in  thru top    of layer L = PSSW4(I,J,L  )
!  (b) Precip going  out thru bottom of layer L = PSSW4(I,J,L-1)
!
!  GEOS-CHEM
!  (a) Precip coming in  thru top    of layer L = PDOWN(L+1,I,J)
!  (b) Precip going  out thru bottom of layer L = PDOWN(L,  I,J) 
!
!  Thus: Precip coming in:  PSSW4(I,J,L  ) is analogous to PDOWN(L+1,I,J).
!        Precip going  out: PSSW4(I,J,L-1) is analogous to PDOWN(L,I,J  ).
!
!  References (see above for full citations):
!  ============================================================================
!  (1 ) Jacob et al, 2000
!  (2 ) Balkanski et al, 1993
!  (3 ) Giorgi & Chaimedes, 1986
!
!  NOTES: 
!  (1 ) WETDEPBIN should be called twice, once with LS = .TRUE. and once
!        with LS = .FALSE.  This will handle both large-scale and
!        convective precipitation. (bmy, 2/28/00)
!  (2 ) Call subroutine MAKE_QQ to construct the QQ and PDOWN precipitation
!        fields before calling WETDEPBIN. (bmy, 2/28/00)
!  (3 ) Since we are working with an (I,J) column, the ordering of the
!        loops goes J - I - L - N.  Dimension arrays DSTT, PDOWN, QQ
!        to take advantage of this optimal configuration (bmy, 2/28/00)
!  (4 ) Use double-precision exponents to force REAL*8 accuracy
!        (e.g. 1d0, bmy, 2/28/00)
!  (5 ) Diagnostics ND16, ND17, ND18, and ND39 use allocatable arrays 
!        from "diag_mod.f"  (bmy, bey, 3/14/00)
!  (6 ) WETDEPBIN only processes soluble tracers and/or aerosols, as are
!        defined in the NSOL and IDWETD arrays (bmy, 3/14/00)
!  (7 ) Add kludge to prevent wet deposition in the stratosphere (bmy, 6/21/00)
!  (8 ) Removed obsolete code from 10/27/00 (bmy, 12/21/00)
!  (9 ) Remove IREF, JREF -- they are obsolete (bmy, 9/27/01)
!  (10) Removed obsolete commented out code from 9/01 (bmy, 10/24/01)
!  (11) Replaced all instances of IM with IIPAR and JM with JJPAR, in order
!        to prevent namespace confusion for the new TPCORE (bmy, 6/25/02)
!  (12) Now reference BXHEIGHT from "dao_mod.f".  Also references routine
!        GEOS_CHEM_STOP from "error_mod.f".  Also fix ND39 diagnostic so that
!        the budget of tracer lost to wetdep is closed.  Now bundled into
!        "wetscav_mod.f".  Now only save to AD16, AD17, AD18, AD39 if L<=LD16,
!        L<=LD17, L<=LD18, and L<=LD39 respectively; this avoids out-of-bounds
!        array errors. Updated comments, cosmetic changes. (qli, bmy, 11/26/02)
!  (13) References IDTSO2, IDTSO4 from "tracerid_mod.f". SO2 in sulfate 
!        chemistry is wet-scavenged on the raindrop and converted to SO4 by 
!        aqueous chem. If evaporation occurs then SO2 comes back as SO4.
!        (rjp, bmy, 3/23/03)  
!  (14) Now use function GET_TS_DYN() from "time_mod.f" (bmy, 3/27/03)
!  (15) Now parallelize over outermost J-loop.  Also move internal routines
!        LS_K_RAIN, LS_F_PRIME, CONV_F_PRIME, and SAFETY to the module, since
!        we cannot call internal routines from w/in a parallel loop. 
!        (bmy, 3/18/04)
!  (16) Now references STT & N_TRACERS from "tracer_mod.f".  Also now make
!        DSTT a 4-d internal array so as to facilitate -C checking on the
!        SGI platform. (bmy, 7/20/04)
!  (17) Now references IDTHg2 from "tracerid_mod.f".  Now pass the amt of
!        Hg2 wet scavenged out of the column to "ocean_mercury_mod.f" via
!        routine ADD_Hg2_WD. (sas, bmy, 1/19/05)
!  (18) Bug fix: replace line that can cause numerical blowup with a safer
!        analytical expression. (bmy, 2/23/05)
!  (19) Block out parallel loop with #ifdef statements for SGI_MIPS compiler.
!        For some reason this causes an error. (bmy, 5/5/05)
!  (20) Now use function IS_Hg2 to determine if a tracer is a tagged Hg2 
!        tracer.  Now also pass N to ADD_Hg2_WD.  Now references LDYNOCEAN
!        from "logical_mod.f".  Now do not call ADD_Hg2_WD if we are not
!        using the dynamic ocean model. (eck, sas, cdh, bmy, 2/27/06)
!  (21) Eliminate unnecessary variables XDSTT, L_PLUS_W.  Also zero all 
!        unused variables for each grid box. (bmy, 5/24/06)
!  (22) Redimension DSTT with NSOL instead of NSOLMAX. In many cases, NSOL is
!        less than NSOLMAX and this will help to save memory especially when
!        running at 2x25 or greater resolution. (bmy, 1/31/08)
!  09 Nov 2012 - M. Payer    - Replaced all met field arrays with State_Met
!                              derived type object
!  17 Oct 2013 - M. Sulprizio- Now pass STT via the State_Chm object
!******************************************************************************
!
      ! References to F90 modules
      USE APM_INIT_MOD,       ONLY : NGCOND,    NSO4
      USE APM_INIT_MOD,       ONLY : NSEA,      NDSTB
      USE APM_INIT_MOD,       ONLY : NCTSO4,    NCTBCOC
      USE APM_INIT_MOD,       ONLY : NCTDST,    NCTSEA
      USE APM_INIT_MOD,       ONLY : NBCPIF,    NBCPIO
      USE APM_INIT_MOD,       ONLY : NOCPIF,    NOCPIO
      USE APM_INIT_MOD,       ONLY : NBCPOF,    NBCPOO
      USE APM_INIT_MOD,       ONLY : NOCPOF,    NOCPOO
      USE APM_INIT_MOD,       ONLY : IDTSO4G,   IDTSO4BIN1, IDTDSTBIN1
      USE APM_INIT_MOD,       ONLY : IDTCTSO4,  IDTCTBCOC,  IDTCTDST
      USE APM_INIT_MOD,       ONLY : IDTCTSEA,  IDTSEABIN1
      USE APM_INIT_MOD,       ONLY : IDTBCPIFF, IDTBCPIBB
      USE APM_INIT_MOD,       ONLY : IDTOCPIFF, IDTOCPIBB
      USE APM_INIT_MOD,       ONLY : IDTBCPOFF, IDTBCPOBB
      USE APM_INIT_MOD,       ONLY : IDTOCPOFF, IDTOCPOBB
      USE APM_DRIV_MOD,       ONLY : IACT1, IACT3, FCLOUD,  GFTOT3D
      USE APM_INIT_MOD,       ONLY : RDRY,  RSALT, IACTSS1, IACTSS3
      USE DIAG_MOD,           ONLY : AD16,  AD17,  AD18
      USE DIAG_MOD,           ONLY : CT16,  CT17,  CT18,    AD39 
      USE ERROR_MOD,          ONLY : GEOS_CHEM_STOP, IT_IS_NAN
      USE GIGC_State_Met_Mod, ONLY : MetState
      USE GIGC_State_Chm_Mod, ONLY : ChmState
      USE LOGICAL_MOD,        ONLY : LDYNOCEAN
      USE TIME_MOD,           ONLY : GET_TS_DYN
      USE WETSCAV_MOD,        ONLY : SO2GAINED, SO2WETLOSS

      USE CMN_SIZE_MOD             ! Size parameters
      USE CMN_DIAG_MOD             ! Diagnostic arrays and switches 

      IMPLICIT NONE
!
! !INPUT PARAMETERS: 
!
      LOGICAL,        INTENT(IN)    :: LS          ! =T for Large-Scale precip
      TYPE(MetState), INTENT(IN)    :: State_Met   ! Meteorology State object
!
! !INPUT/OUTPUT PARAMETERS:
!
      TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
!
! !LOCAL VARIABLES:
!
      LOGICAL, SAVE       :: FIRST = .TRUE.
      LOGICAL             :: AER

      INTEGER             :: I, IDX, J, L, N, NN
      
      REAL*8              :: Q,     QDOWN,  DT,        DT_OVER_TAU
      REAL*8              :: K,     K_MIN,  K_RAIN,    RAINFRAC
      REAL*8              :: F,     FTOP,   F_PRIME,   WASHFRAC
      REAL*8              :: LOST,  GAINED, MASS_WASH, MASS_NOWASH
      REAL*8              :: ALPHA, ALPHA2, WETLOSS,   TMP
      REAL*8              :: RSIZEIN, IACT, IACTSS

      ! DSTT is the accumulator array of rained-out 
      ! soluble tracer for a given (I,J) column
      REAL*8              :: DSTT(NSOL,LLPAR,IIPAR,JJPAR)
      REAL*8              :: MASS1, MASS2, MASS3, MASS4
 
      ! Make a pointer to the tracer array
      REAL*8, POINTER :: STT(:,:,:,:)

      !=================================================================
      ! WETDEPBIN begins here!
      !
      ! (1)  I n i t i a l i z e   V a r i a b l e s
      !=================================================================

      ! Point to STT
      STT => State_Chm%TRACERS

      ! Dynamic timestep [s]
      DT    = GET_TS_DYN() * 60d0
      
      ! Select index for diagnostic arrays -- will archive either
      ! large-scale or convective rainout/washout fractions
      IF ( LS ) THEN
         IDX = 1
      ELSE
         IDX = 2
      ENDIF

      !=================================================================
      ! (2)  L o o p   O v e r   (I, J)   S u r f a c e   B o x e s
      !
      ! Process rainout / washout by columns.
      !=================================================================

#if   !defined( SGI_MIPS )
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I,       J,           FTOP,      ALPHA              )
!$OMP+PRIVATE( ALPHA2,  F,           F_PRIME,   GAINED,   K_RAIN   )
!$OMP+PRIVATE( LOST,    MASS_NOWASH, MASS_WASH, RAINFRAC, WASHFRAC )
!$OMP+PRIVATE( WETLOSS, L,           Q,         NN,       N        )
!$OMP+PRIVATE( QDOWN,   AER,         TMP                           )
!$OMP+PRIVATE( IACT,    IACTSS, RSIZEIN                            )
!$OMP+PRIVATE( MASS1, MASS2, MASS3, MASS4                          )
!$OMP+SCHEDULE( DYNAMIC )
#endif
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! Zero FTOP
         FTOP = 0d0

         ! Zero accumulator array
         DO L  = 1, LLPAR
         DO NN = 1, NSOL
            DSTT(NN,L,I,J) = 0d0
         ENDDO
         ENDDO

         !==============================================================
         ! (3)  R a i n o u t   F r o m   T o p   L a y e r  (L = LLPAR) 
         !
         ! Assume that rainout is happening in the top layer if 
         ! QQ(LLPAR,I,J) > 0.  In other words, if any precipitation 
         ! forms in grid box (I,J,LLPAR), assume that all of it falls 
         ! down to lower levels.
         !
         ! Soluble gases/aerosols are incorporated into the raindrops 
         ! and are completely removed from grid box (I,J,LLPAR).  There 
         ! is no evaporation and "resuspension" of aerosols during a 
         ! rainout event.
         !
         ! For large-scale (a.k.a. stratiform) precipitation, the first 
         ! order rate constant for rainout in the grid box (I,J,L=LLPAR) 
         ! (cf. Eq. 12, Jacob et al, 2000) is given by:
         !
         !                        Q        
         !    K_RAIN = K_MIN + -------    [units: s^-1]
         !                      L + W    
         !          
         ! and the areal fraction of grid box (I,J,L=LLPAR) that 
         ! is actually experiencing large-scale precipitation 
         ! (cf. Eq. 11, Jacob et al, 2000) is given by: 
         ! 
         !                  Q               
         !    F' =  -------------------   [unitless]
         !           K_RAIN * ( L + W )    
         !
         ! Where:
         !
         !    K_MIN  = minimum value for K_RAIN         
         !           = 1.0e-4 [s^-1]
         !
         !    L + W  = condensed water content in cloud 
         !           = 1.5e-6 [cm3 H2O/cm3 air]
         !
         !    Q = QQ = rate of precipitation formation 
         !             [ cm3 H2O / cm3 air / s ]
         !
         ! For convective precipitation, K_RAIN = 5.0e-3 [s^-1], and the
         ! expression for F' (cf. Eq. 13, Jacob et al, 2000) becomes:
         !
         !                                  { DT        }
         !                    FMAX * Q * MIN{ --- , 1.0 }
         !                                  { TAU       }
         !  F' = ------------------------------------------------------
         !               { DT        }
         !        Q * MIN{ --- , 1.0 }  +  FMAX * K_RAIN * ( L + W )
         !               { TAU       } 
         !
         ! Where:
         !
         !    Q = QQ = rate of precipitation formation 
         !             [cm3 H2O/cm3 air/s]
         !
         !    FMAX   = maximum value for F' 
         !           = 0.3
         !
         !    DT     = dynamic time step from the CTM [s]
         !
         !    TAU    = duration of rainout event 
         !           = 1800 s (30 min)
         !
         !    L + W  = condensed water content in cloud 
         !           = 2.0e-6 [cm3 H2O/cm3 air]
         !
         ! K_RAIN and F' are needed to compute the fraction of tracer
         ! in grid box (I,J,L=LLPAR) lost to rainout.  This is done in 
         ! module routine RAINOUT.
         !==============================================================

         ! Zero variables for this level
         ALPHA       = 0d0
         ALPHA2      = 0d0
         F           = 0d0
         F_PRIME     = 0d0
         GAINED      = 0d0
         K_RAIN      = 0d0
         LOST        = 0d0
         Q           = 0d0
         QDOWN       = 0d0
         MASS_NOWASH = 0d0
         MASS_WASH   = 0d0
         RAINFRAC    = 0d0
         WASHFRAC    = 0d0
         WETLOSS     = 0d0

         ! Start at the top of the atmosphere
         L = LLPAR

         ! If precip forms at (I,J,L), assume it all rains out
         IF ( QQ(L,I,J) > 0d0 ) THEN

            ! Q is the new precip that is forming within grid box (I,J,L)
            Q = QQ(L,I,J)

            ! Compute K_RAIN and F' for either large-scale or convective
            ! precipitation (cf. Eqs. 11-13, Jacob et al, 2000) 
            IF ( LS ) THEN
               K_RAIN  = LS_K_RAIN( Q )
               F_PRIME = LS_F_PRIME( Q, K_RAIN )
               IACT = IACT3(I,J,L)
               IACTSS = IACTSS3
            ELSE
               K_RAIN  = 1.5d-3
               F_PRIME = CONV_F_PRIME( Q, K_RAIN, DT )
               IACT = IACT1(I,J,L)
               IACTSS = IACTSS1
            ENDIF
            
            ! Set F = F', since there is no FTOP at L = LLPAR
            F = F_PRIME

            ! Only compute rainout if F > 0. 
            ! This helps to eliminate unnecessary CPU cycles.
            IF ( F > 0d0 ) THEN 

               MASS1=SUM(STT(I,J,L,IDTSO4BIN1:(IDTSO4BIN1+NSO4-1)))
               MASS3=SUM(STT(I,J,L,IDTSEABIN1:(IDTSEABIN1+NSEA-1)))

               ! Loop over soluble tracers and/or aerosol tracers    
               DO NN = 1, NSOL
                  N = IDWETD(NN)
                  if(n<idtdstbin1.or.n>=(idtdstbin1+ndstb))then

                  ! Call subroutine RAINOUT to compute the fraction
                  ! of tracer lost to rainout in grid box (I,J,L=LLPAR)
                  CALL RAINOUT( I, J, L, N, K_RAIN, DT, F, RAINFRAC,
     &                          State_Met )

                  IF(.NOT.((N>=IDTSO4BIN1.and.N<(IDTSO4BIN1+IACT)).or.
     &               (N>=IDTSEABIN1.and.N<(IDTSEABIN1+IACTSS)).or.
     &               (N>=IDTCTSO4.and.N<(IDTCTSO4+NCTSO4)).or.
     &               (N>=IDTCTSEA.and.N<(IDTCTSEA+NCTSEA))))THEN
                  ! WETLOSS is the amount of soluble tracer 
                  ! lost to rainout in grid box (I,J,L=LLPAR)
                  WETLOSS = STT(I,J,L,N) * RAINFRAC

                  ! Remove rainout losses in grid box (I,J,L=LLPAR) from STT
                  STT(I,J,L,N) = STT(I,J,L,N) - WETLOSS

                  ! DSTT is an accumulator array for rained-out tracers.  
                  ! The tracers in DSTT are in the liquid phase and will 
                  ! precipitate to the levels below until a washout occurs.
                  ! Initialize DSTT at (I,J,L=LLPAR) with WETLOSS.
                  DSTT(NN,L,I,J) = WETLOSS

                  ! Negative tracer...call subroutine SAFETY
                  IF ( STT(I,J,L,N) < 0d0 ) THEN
                     CALL SAFETY( I, J, L, N, 3, 
     &                            LS,             PDOWN(L,I,J), 
     &                            QQ(L,I,J),      ALPHA,      
     &                            ALPHA2,         RAINFRAC,    
     &                            WASHFRAC,       MASS_WASH,    
     &                            MASS_NOWASH,    WETLOSS,    
     &                            GAINED,         LOST,        
     &                            DSTT(NN,:,I,J), STT(I,J,:,N) )
                  ENDIF

                  ENDIF
                  endif
               ENDDO

               MASS2=SUM(STT(I,J,L,IDTSO4BIN1:(IDTSO4BIN1+NSO4-1)))
               MASS4=SUM(STT(I,J,L,IDTSEABIN1:(IDTSEABIN1+NSEA-1)))

               ! Loop over soluble tracers and/or aerosol tracers
               DO NN = 1, NCTSO4
                  N = IDTCTSO4+NN-1
                  IF(MASS1>1.D-30)THEN
                    DSTT((N-IDTSO4G+1),L,I,J)=
     &              STT(I,J,L,N)*(1.D0-(MASS2/MASS1))
                    STT(I,J,L,N)=STT(I,J,L,N)*(MASS2/MASS1)
                  ENDIF
                  ! Negative tracer...call subroutine SAFETY
                  IF ( STT(I,J,L,N) < 0d0 ) THEN
                     CALL SAFETY( I, J, L, N, 3,
     &                            LS,             PDOWN(L,I,J),
     &                            QQ(L,I,J),      ALPHA,
     &                            ALPHA2,         RAINFRAC,
     &                            WASHFRAC,       MASS_WASH,
     &                            MASS_NOWASH,    WETLOSS,
     &                            GAINED,         LOST,
     &                            DSTT(NN,:,I,J), STT(I,J,:,N) )
                  ENDIF
               ENDDO
               DO NN = 1, NCTSEA
                  N = IDTCTSEA+NN-1
                  IF(MASS3>1.D-30)THEN
                    DSTT((N-IDTSO4G+1),L,I,J)=
     &              STT(I,J,L,N)*(1.D0-(MASS4/MASS3))
                    STT(I,J,L,N)=STT(I,J,L,N)*(MASS4/MASS3)
                  ENDIF
                  ! Negative tracer...call subroutine SAFETY
                  IF ( STT(I,J,L,N) < 0d0 ) THEN
                     CALL SAFETY( I, J, L, N, 3,
     &                            LS,             PDOWN(L,I,J),
     &                            QQ(L,I,J),      ALPHA,
     &                            ALPHA2,         RAINFRAC,
     &                            WASHFRAC,       MASS_WASH,
     &                            MASS_NOWASH,    WETLOSS,
     &                            GAINED,         LOST,
     &                            DSTT(NN,:,I,J), STT(I,J,:,N) )
                  ENDIF
               ENDDO

            ENDIF

            ! Save FTOP for the next lower level 
            FTOP = F
         ENDIF

         !==============================================================
         ! (4)  R a i n o u t   i n   t h e   M i d d l e   L e v e l s
         ! 
         ! Rainout occurs when there is more precipitation in grid box 
         ! (I,J,L) than in grid box (I,J,L+1).  In other words, rainout 
         ! occurs when the amount of rain falling through the bottom of 
         ! grid box (I,J,L) is more than the amount of rain coming in 
         ! through the top of grid box (I,J,L). 
         !
         ! Thus ( PDOWN(L,I,J) > 0 and QQ(L,I,J) > 0 ) is the 
         ! criterion for Rainout.
         !
         ! Soluble gases/aerosols are incorporated into the raindrops 
         ! and are completely removed from grid box (I,J,L).  There is 
         ! no evaporation and "resuspension" of aerosols during a 
         ! rainout event.
         !
         ! Compute K_RAIN and F' for grid box (I,J,L) exactly as 
         ! described above in Section (4).  K_RAIN and F' depend on 
         ! whether we have large-scale or convective precipitation.
         !
         ! F' is the areal fraction of grid box (I,J,L) that is 
         ! precipitating.  However, the effective area of precipitation
         ! that layer L sees (cf. Eqs. 11-13, Jacob et al, 2000) is 
         ! given by:
         !
         !                   F = MAX( F', FTOP )
         !
         ! where FTOP = F' at grid box (I,J,L+1), that is, for the grid
         ! box immediately above the current grid box.  
         !
         ! Therefore, the effective area of precipitation in grid box
         ! (I,J,L) depends on the area of precipitation in the grid 
         ! boxes above it.
         !
         ! Having computed K_RAIN and F for grid box (I,J,L), call 
         ! routine RAINOUT to compute the fraction of tracer lost to 
         ! rainout conditions.
         !==============================================================
         DO L = LLPAR-1, 2, -1

            ! Zero variables for each level
            ALPHA       = 0d0
            ALPHA2      = 0d0
            F           = 0d0
            F_PRIME     = 0d0
            GAINED      = 0d0
            K_RAIN      = 0d0
            LOST        = 0d0
            MASS_NOWASH = 0d0
            MASS_WASH   = 0d0
            Q           = 0d0
            QDOWN       = 0d0
            RAINFRAC    = 0d0
            WASHFRAC    = 0d0
            WETLOSS     = 0d0

            ! Rainout criteria
            IF ( PDOWN(L,I,J) > 0d0 .and. QQ(L,I,J) > 0d0 ) THEN

               ! Q is the new precip that is forming within grid box (I,J,L)
               Q = QQ(L,I,J)

               ! Compute K_RAIN and F' for either large-scale or convective
               ! precipitation (cf. Eqs. 11-13, Jacob et al, 2000) 
               IF ( LS ) THEN
                  K_RAIN  = LS_K_RAIN( Q )
                  F_PRIME = LS_F_PRIME( Q, K_RAIN )
                  IACT = IACT3(I,J,L)
                  IACTSS = IACTSS3
               ELSE
                  K_RAIN  = 1.5d-3
                  F_PRIME = CONV_F_PRIME( Q, K_RAIN, DT )
                  IACT = IACT1(I,J,L)
                  IACTSS = IACTSS1
               ENDIF

               ! F is the effective area of precip seen by grid box (I,J,L) 
               F = MAX( F_PRIME, FTOP )

               ! Only compute rainout if F > 0. 
               ! This helps to eliminate unnecessary CPU cycles. 
               IF ( F > 0d0 ) THEN

               MASS1=SUM(STT(I,J,L,IDTSO4BIN1:(IDTSO4BIN1+NSO4-1)))
               MASS3=SUM(STT(I,J,L,IDTSEABIN1:(IDTSEABIN1+NSEA-1)))

                  ! Loop over soluble tracers and/or aerosol tracers    
                  DO NN = 1, NSOL
                     N = IDWETD(NN)
                     if(n<idtdstbin1.or.n>=(idtdstbin1+ndstb))then

                     ! Call subroutine RAINOUT to comptue the fraction
                     ! of tracer lost to rainout in grid box (I,J,L) 
                     CALL RAINOUT( I, J, L, N, K_RAIN, DT, F, RAINFRAC,
     &                             State_Met )

                  IF(.NOT.((N>=IDTSO4BIN1.and.N<(IDTSO4BIN1+IACT)).or.
     &               (N>=IDTSEABIN1.and.N<(IDTSEABIN1+IACTSS)).or.
     &               (N>=IDTCTSO4.and.N<(IDTCTSO4+NCTSO4)).or.
     &               (N>=IDTCTSEA.and.N<(IDTCTSEA+NCTSEA))))THEN
                     ! WETLOSS is the amount of tracer in grid box 
                     ! (I,J,L) that is lost to rainout.
                     WETLOSS = STT(I,J,L,N) * RAINFRAC

                     ! Subtract the rainout loss in grid box (I,J,L) from STT
                     STT(I,J,L,N) = STT(I,J,L,N) - WETLOSS

                     ! Add to DSTT the tracer lost to rainout in grid box 
                     ! (I,J,L) plus the tracer lost to rainout from grid box 
                     ! (I,J,L+1), which has by now precipitated down into 
                     ! grid box (I,J,L).  DSTT will continue to accumulate 
                     ! rained out tracer in this manner until a washout 
                     ! event occurs.
                     DSTT(NN,L,I,J) = DSTT(NN,L+1,I,J) + WETLOSS

                     ! Negative tracer...call subroutine SAFETY
                     IF ( STT(I,J,L,N) < 0d0 .or.
     &                    IT_IS_NAN( STT(I,J,L,N) ) ) THEN
                        CALL SAFETY( I, J, L, N, 4, 
     &                               LS,             PDOWN(L,I,J),  
     &                               QQ(L,I,J),      ALPHA,        
     &                               ALPHA2,         RAINFRAC,     
     &                               WASHFRAC,       MASS_WASH,    
     &                               MASS_NOWASH,    WETLOSS,      
     &                               GAINED,         LOST,         
     &                               DSTT(NN,:,I,J), STT(I,J,:,N) )
                     ENDIF

                     ENDIF
                     endif
                  ENDDO

               MASS2=SUM(STT(I,J,L,IDTSO4BIN1:(IDTSO4BIN1+NSO4-1)))
               MASS4=SUM(STT(I,J,L,IDTSEABIN1:(IDTSEABIN1+NSEA-1)))

               ! Loop over soluble tracers and/or aerosol tracers
               DO NN = 1, NCTSO4
                  N = IDTCTSO4+NN-1
                  IF(MASS1>1.D-30)THEN
                    DSTT((N-IDTSO4G+1),L,I,J)=
     &              DSTT((N-IDTSO4G+1),L+1,I,J)+
     &              STT(I,J,L,N)*(1.D0-(MASS2/MASS1))
                    STT(I,J,L,N)=STT(I,J,L,N)*(MASS2/MASS1)
                  ENDIF
                  ! Negative tracer...call subroutine SAFETY
                  IF ( STT(I,J,L,N) < 0d0 ) THEN
                     CALL SAFETY( I, J, L, N, 4,
     &                            LS,             PDOWN(L,I,J),
     &                            QQ(L,I,J),      ALPHA,
     &                            ALPHA2,         RAINFRAC,
     &                            WASHFRAC,       MASS_WASH,
     &                            MASS_NOWASH,    WETLOSS,
     &                            GAINED,         LOST,
     &                            DSTT(NN,:,I,J), STT(I,J,:,N) )
                  ENDIF
               ENDDO
               DO NN = 1, NCTSEA
                  N = IDTCTSEA+NN-1
                  IF(MASS3>1.D-30)THEN
                    DSTT((N-IDTSO4G+1),L,I,J)=
     &              DSTT((N-IDTSO4G+1),L+1,I,J)+
     &              STT(I,J,L,N)*(1.D0-(MASS4/MASS3))
                    STT(I,J,L,N)=STT(I,J,L,N)*(MASS4/MASS3)
                  ENDIF
                  ! Negative tracer...call subroutine SAFETY
                  IF ( STT(I,J,L,N) < 0d0 ) THEN
                     CALL SAFETY( I, J, L, N, 4,
     &                            LS,             PDOWN(L,I,J),
     &                            QQ(L,I,J),      ALPHA,
     &                            ALPHA2,         RAINFRAC,
     &                            WASHFRAC,       MASS_WASH,
     &                            MASS_NOWASH,    WETLOSS,
     &                            GAINED,         LOST,
     &                            DSTT(NN,:,I,J), STT(I,J,:,N) )
                  ENDIF
               ENDDO

               ENDIF

               ! Save FTOP for next level
               FTOP = F 

            !==============================================================
            ! (5)  W a s h o u t   i n   t h e   m i d d l e   l e v e l s
            ! 
            ! Washout occurs when we have evaporation (or no precipitation 
            ! at all) at grid box (I,J,L), but have rain coming down from 
            ! grid box (I,J,L+1).
            !
            ! Thus PDOWN(L,I,J) > 0 and QQ(L,I,J) <= 0 is the criterion 
            ! for Washout.  Also recall that QQ(L,I,J) < 0 denotes 
            ! evaporation and not precipitation.
            !
            ! A fraction ALPHA of the raindrops falling down from grid 
            ! box (I,J,L+1) to grid box (I,J,L) will evaporate along the 
            ! way.  ALPHA is given by:
            !  
            !            precip leaving (I,J,L+1) - precip leaving (I,J,L)
            !  ALPHA = ---------------------------------------------------
            !                     precip leaving (I,J,L+1)
            !
            !
            !                    -QQ(L,I,J) * DZ(I,J,L)
            !        =         --------------------------
            !                        PDOWN(L+1,I,J)
            !
            ! We assume that a fraction ALPHA2 = 0.5 * ALPHA of the 
            ! previously rained-out aerosols and HNO3 coming down from 
            ! level (I,J,L+1) will evaporate and re-enter the atmosphere 
            ! in the gas phase in grid box (I,J,L).  This process is 
            ! called "resuspension".  
            !
            ! For non-aerosol species, the amount of previously rained 
            ! out mass coming down from grid box (I,J,L+1) to grid box 
            ! (I,J,L) is figured into the total mass available for 
            ! washout in grid box (I,J,L).  We therefore do not have to
            ! use the fraction ALPHA2 to compute the resuspension.
            !
            ! NOTE from Hongyu Liu about ALPHA (hyl, 2/29/00)
            ! =============================================================
            ! If our QQ field was perfect, the evaporated amount in grid 
            ! box (I,J,L) would be at most the total rain amount coming 
            ! from above (i.e. PDOWN(I,J,L+1) ). But this is not true for 
            ! the MOISTQ field we are using.  Sometimes the evaporation in 
            ! grid box (I,J,L) can be more than the rain amount from above.  
            ! The reason is our "evaporation" also includes the effect of 
            ! cloud detrainment.  For now we cannot find a way to 
            ! distinguish betweeen the two. We then decided to release 
            ! aerosols in both the detrained air and the evaporated air. 
            !
            ! Therefore, we should use this term in the numerator:
            ! 
            !                -QQ(I,J,L) * BXHEIGHT(I,J,L) 
            !
            ! instead of the term:
            ! 
            !                PDOWN(L+1)-PDOWN(L)
            !
            ! Recall that in make_qq.f we have restricted PDOWN to 
            ! positive values, otherwise, QQ would be equal to 
            ! PDOWN(L+1)-PDOWN(L).           
            !==============================================================
            ELSE IF ( PDOWN(L,I,J) > 0d0 .and. QQ(L,I,J) <= 0d0 ) THEN

               ! QDOWN is the precip leaving thru the bottom of box (I,J,L)
               ! Q     is the new precip that is forming within box (I,J,L)
               QDOWN = PDOWN(L,I,J)
               Q     = QQ(L,I,J)

               ! Since no precipitation is forming within grid box (I,J,L),
               ! F' = 0, and F = MAX( F', FTOP ) reduces to F = FTOP.
               F = FTOP
 
               ! Only compute washout if F > 0.
               ! This helps to eliminate needless CPU cycles.
               IF ( F > 0d0 ) THEN

               MASS1=SUM(STT(I,J,L,IDTSO4BIN1:(IDTSO4BIN1+NSO4-1)))
               MASS3=SUM(STT(I,J,L,IDTSEABIN1:(IDTSEABIN1+NSEA-1)))

                  ! Loop over soluble tracers and/or aerosol tracers    
                  DO NN = 1, NSOL
                     N  = IDWETD(NN)
                     if(n<idtdstbin1.or.n>=(idtdstbin1+ndstb))then

                  IF(.NOT.((N>=IDTCTSO4.and.N<(IDTCTSO4+NCTSO4)).or.
     &               (N>=IDTCTSEA.and.N<(IDTCTSEA+NCTSEA))))THEN
                     ! Call WASHOUT to compute the fraction of 
                     ! tracer lost to washout in grid box (I,J,L)
                     CALL WASHOUT( I,     J,  L, N, 
     &                             QDOWN, DT, F, WASHFRAC, AER,
     &                             State_Met )
                  
                     IF( N >= IDTSO4BIN1 .and.
     &                   N < (IDTSO4BIN1+NSO4) )THEN
                       IF ( State_Met%T(I,J,L) > 268.D0 ) THEN
                       RSIZEIN=GFTOT3D(I,J,L,1)*RDRY(N-IDTSO4BIN1+1)
                       CALL SIZEWASHFRAC( RSIZEIN, QDOWN, WASHFRAC )
                       WASHFRAC=(1.D0-EXP(-WASHFRAC*DT))
                       ELSE
                       WASHFRAC=0.D0
                       ENDIF
                     ENDIF

                     IF( N >= IDTSEABIN1 .and.
     &                   N < (IDTSEABIN1+NSEA) )THEN
                       IF ( State_Met%T(I,J,L) > 268.D0 ) THEN
                       RSIZEIN=GFTOT3D(I,J,L,2)*RSALT(N-IDTSEABIN1+1)
                       CALL SIZEWASHFRAC( RSIZEIN, QDOWN, WASHFRAC )
                       WASHFRAC=(1.D0-EXP(-WASHFRAC*DT))
                       ELSE
                       WASHFRAC=0.D0
                       ENDIF
                     ENDIF

                     !=====================================================
                     ! Washout of aerosol tracers -- 
                     ! this is modeled as a kinetic process
                     !=====================================================
                     IF ( AER ) THEN

                        ! Amount of aerosol lost to washout in grid box
                        ! (qli, bmy, 10/29/02)
                        WETLOSS = STT(I,J,L,N) * WASHFRAC

                        ! Remove washout losses in grid box (I,J,L) from STT.
                        ! Add the aerosol that was reevaporated in (I,J,L).
                        ! SO2 in sulfate chemistry is wet-scavenged on the
                        ! raindrop and converted to SO4 by aqeuous chem.
                        ! If evaporation occurs then SO2 comes back as SO4
                        ! (rjp, bmy, 3/23/03)
                        STT(I,J,L,N) = STT(I,J,L,N) - WETLOSS

                        ! Add the washed out tracer from grid box (I,J,L) to 
                        ! DSTT.  Also add the amount of tracer coming down
                        ! from grid box (I,J,L+1) that does NOT re-evaporate.
                        DSTT(NN,L,I,J) = DSTT(NN,L+1,I,J) + WETLOSS

                     !=====================================================
                     ! Washout of non-aerosol tracers
                     ! This is modeled as an equilibrium process
                     !=====================================================
                     ELSE
                  
                        ! MASS_NOWASH is the amount of non-aerosol tracer in 
                        ! grid box (I,J,L) that is NOT available for washout.
                        MASS_NOWASH = ( 1d0 - F ) * STT(I,J,L,N)
                     
                        ! MASS_WASH is the total amount of non-aerosol tracer
                        ! that is available for washout in grid box (I,J,L).
                        ! It consists of the mass in the precipitating
                        ! part of box (I,J,L), plus the previously rained-out
                        ! tracer coming down from grid box (I,J,L+1).
                        ! (Eq. 15, Jacob et al, 2000).
                        MASS_WASH = ( F*STT(I,J,L,N) ) +DSTT(NN,L+1,I,J)

                        ! WETLOSS is the amount of tracer mass in 
                        ! grid box (I,J,L) that is lost to washout.
                        ! (Eq. 16, Jacob et al, 2000)
                        WETLOSS = MASS_WASH * WASHFRAC -DSTT(NN,L+1,I,J)

                        ! The tracer left in grid box (I,J,L) is what was
                        ! in originally in the non-precipitating fraction 
                        ! of the box, plus MASS_WASH, less WETLOSS. 
                        STT(I,J,L,N) = STT(I,J,L,N) - WETLOSS  
                  
                        ! Add washout losses in grid box (I,J,L) to DSTT 
                        DSTT(NN,L,I,J) = DSTT(NN,L+1,I,J) + WETLOSS

                     ENDIF

                     ! Negative tracer...call subroutine SAFETY
                     IF ( STT(I,J,L,N) < 0d0 .or. 
     &                    IT_IS_NAN( STT(I,J,L,N) ) ) THEN
                        CALL SAFETY( I, J, L, N, 5, 
     &                               LS,             PDOWN(L,I,J), 
     &                               QQ(L,I,J),      ALPHA,        
     &                               ALPHA2,         RAINFRAC,     
     &                               WASHFRAC,       MASS_WASH,    
     &                               MASS_NOWASH,    WETLOSS,      
     &                               GAINED,         LOST,         
     &                               DSTT(NN,:,I,J), STT(I,J,:,N) )
                     ENDIF

                  ENDIF

                  endif
                  ENDDO

               MASS2=SUM(STT(I,J,L,IDTSO4BIN1:(IDTSO4BIN1+NSO4-1)))
               MASS4=SUM(STT(I,J,L,IDTSEABIN1:(IDTSEABIN1+NSEA-1)))

               ! Loop over soluble tracers and/or aerosol tracers
               DO NN = 1, NCTSO4
                  N = IDTCTSO4+NN-1
                  IF(MASS1>1.D-30)THEN
                    DSTT((N-IDTSO4G+1),L,I,J)=
     &              DSTT((N-IDTSO4G+1),L+1,I,J)+
     &              STT(I,J,L,N)*(1.D0-(MASS2/MASS1))
                    STT(I,J,L,N)=STT(I,J,L,N)*(MASS2/MASS1)
                  ENDIF
                     ! Negative tracer...call subroutine SAFETY
                     IF ( STT(I,J,L,N) < 0d0 .or.
     &                    IT_IS_NAN( STT(I,J,L,N) ) ) THEN
                        CALL SAFETY( I, J, L, N, 5,
     &                               LS,             PDOWN(L,I,J),
     &                               QQ(L,I,J),      ALPHA,
     &                               ALPHA2,         RAINFRAC,
     &                               WASHFRAC,       MASS_WASH,
     &                               MASS_NOWASH,    WETLOSS,
     &                               GAINED,         LOST,
     &                               DSTT(NN,:,I,J), STT(I,J,:,N) )
                     ENDIF
               ENDDO
               DO NN = 1, NCTSEA
                  N = IDTCTSEA+NN-1
                  IF(MASS3>1.D-30)THEN
                    DSTT((N-IDTSO4G+1),L,I,J)=
     &              DSTT((N-IDTSO4G+1),L+1,I,J)+
     &              STT(I,J,L,N)*(1.D0-(MASS4/MASS3))
                    STT(I,J,L,N)=STT(I,J,L,N)*(MASS4/MASS3)
                  ENDIF
                     ! Negative tracer...call subroutine SAFETY
                     IF ( STT(I,J,L,N) < 0d0 .or.
     &                    IT_IS_NAN( STT(I,J,L,N) ) ) THEN
                        CALL SAFETY( I, J, L, N, 5,
     &                               LS,             PDOWN(L,I,J),
     &                               QQ(L,I,J),      ALPHA,
     &                               ALPHA2,         RAINFRAC,
     &                               WASHFRAC,       MASS_WASH,
     &                               MASS_NOWASH,    WETLOSS,
     &                               GAINED,         LOST,
     &                               DSTT(NN,:,I,J), STT(I,J,:,N) )
                     ENDIF
               ENDDO

               ! Loop over soluble tracers and/or aerosol tracers
               DO NN = 1, NSOL
                  N  = IDWETD(NN)
                  if(n<idtdstbin1.or.n>=(idtdstbin1+ndstb))then

                  IF ( AER ) THEN

                     ! ALPHA is the fraction of the raindrops that
                     ! re-evaporate when falling from (I,J,L+1) to (I,J,L)
                     ALPHA = ( ABS( Q ) * State_Met%BXHEIGHT(I,J,L) *
     &                         100d0 )  / PDOWN(L+1,I,J)

                     ! ALPHA2 is the fraction of the rained-out aerosols
                     ! that gets resuspended in grid box (I,J,L)
                     ALPHA2 = 0.5d0 * ALPHA

                     ! GAINED is the rained out aerosol coming down from
                     ! grid box (I,J,L+1) that will evaporate and re-enter
                     ! the atmosphere in the gas phase in grid box (I,J,L).
                     GAINED = DSTT(NN,L+1,I,J) * ALPHA2

                     ! Amount of aerosol lost to washout in grid box
                     ! (qli, bmy, 10/29/02)
                     WETLOSS = - GAINED

                     ! Remove washout losses in grid box (I,J,L) from STT.
                     ! Add the aerosol that was reevaporated in (I,J,L).
                     ! SO2 in sulfate chemistry is wet-scavenged on the
                     ! raindrop and converted to SO4 by aqeuous chem.
                     ! If evaporation occurs then SO2 comes back as SO4
                     ! (rjp, bmy, 3/23/03)
                     STT(I,J,L,N)      = STT(I,J,L,N) - WETLOSS

                     ! Add the washed out tracer from grid box (I,J,L) to
                     ! DSTT.  Also add the amount of tracer coming down
                     ! from grid box (I,J,L+1) that does NOT re-evaporate.
                     DSTT(NN,L,I,J) = DSTT(NN,L+1,I,J) + WETLOSS

                     ! Negative tracer...call subroutine SAFETY
                     IF ( STT(I,J,L,N) < 0d0 .or.
     &                    IT_IS_NAN( STT(I,J,L,N) ) ) THEN
                        CALL SAFETY( I, J, L, N, 5,
     &                               LS,             PDOWN(L,I,J),
     &                               QQ(L,I,J),      ALPHA,
     &                               ALPHA2,         RAINFRAC,
     &                               WASHFRAC,       MASS_WASH,
     &                               MASS_NOWASH,    WETLOSS,
     &                               GAINED,         LOST,
     &                               DSTT(NN,:,I,J), STT(I,J,:,N) )
                     ENDIF

                  ENDIF

                  endif
               ENDDO

               !********************************************************
               ! deal with SO2GAINED
               DO N = IDTSO4BIN1,(IDTSO4BIN1+NSO4-1)
                  STT(I,J,L,N) = STT(I,J,L,N)
     &                         + SO2GAINED(I,J,L)
     &                         * FCLOUD(I,J,L,(N-IDTSO4BIN1+1))
                  ! Negative tracer...call subroutine SAFETY
                  IF ( STT(I,J,L,N) < 0d0 ) THEN
                     CALL SAFETY( I, J, L, N, 51,
     &                            LS,             PDOWN(L,I,J),
     &                            QQ(L,I,J),      ALPHA,
     &                            ALPHA2,         RAINFRAC,
     &                            WASHFRAC,       MASS_WASH,
     &                            MASS_NOWASH,    WETLOSS,
     &                            GAINED,         LOST,
     &                            DSTT(NN,:,I,J), STT(I,J,:,N) )
                  ENDIF
               ENDDO
               IF(NCTSO4>1)THEN
               DO N=IDTCTSO4,(IDTCTSO4+NCTSO4-1)
                  STT(I,J,L,N) = STT(I,J,L,N)
     &                         + SO2GAINED(I,J,L)
     &                         * SUM(FCLOUD(I,J,L,1:NSO4))
                  ! Negative tracer...call subroutine SAFETY
                  IF ( STT(I,J,L,N) < 0d0 ) THEN
                     CALL SAFETY( I, J, L, N, 51,
     &                            LS,             PDOWN(L,I,J),
     &                            QQ(L,I,J),      ALPHA,
     &                            ALPHA2,         RAINFRAC,
     &                            WASHFRAC,       MASS_WASH,
     &                            MASS_NOWASH,    WETLOSS,
     &                            GAINED,         LOST,
     &                            DSTT(NN,:,I,J), STT(I,J,:,N) )
                  ENDIF
               ENDDO
               ENDIF
               DO N = IDTCTBCOC,IDTCTBCOC
                  STT(I,J,L,N) = STT(I,J,L,N)
     &                         + SO2GAINED(I,J,L)
     &                         * FCLOUD(I,J,L,(NSO4+1))
                  ! Negative tracer...call subroutine SAFETY
                  IF ( STT(I,J,L,N) < 0d0 ) THEN
                     CALL SAFETY( I, J, L, N, 51,
     &                            LS,             PDOWN(L,I,J),
     &                            QQ(L,I,J),      ALPHA,
     &                            ALPHA2,         RAINFRAC,
     &                            WASHFRAC,       MASS_WASH,
     &                            MASS_NOWASH,    WETLOSS,
     &                            GAINED,         LOST,
     &                            DSTT(NN,:,I,J), STT(I,J,:,N) )
                  ENDIF
               ENDDO
               DO N = (IDTCTBCOC+1),(IDTCTBCOC+1)
                  STT(I,J,L,N) = STT(I,J,L,N)
     &                         + SO2GAINED(I,J,L)
     &                         * FCLOUD(I,J,L,(NSO4+2))
                  ! Negative tracer...call subroutine SAFETY
                  IF ( STT(I,J,L,N) < 0d0 ) THEN
                     CALL SAFETY( I, J, L, N, 51,
     &                            LS,             PDOWN(L,I,J),
     &                            QQ(L,I,J),      ALPHA,
     &                            ALPHA2,         RAINFRAC,
     &                            WASHFRAC,       MASS_WASH,
     &                            MASS_NOWASH,    WETLOSS,
     &                            GAINED,         LOST,
     &                            DSTT(NN,:,I,J), STT(I,J,:,N) )
                  ENDIF
               ENDDO
               DO N=IDTCTDST,IDTCTDST
                  STT(I,J,L,N) = STT(I,J,L,N)
     &                         + SO2GAINED(I,J,L)
     &                         * FCLOUD(I,J,L,(NSO4+3))
                  ! Negative tracer...call subroutine SAFETY
                  IF ( STT(I,J,L,N) < 0d0 ) THEN
                     CALL SAFETY( I, J, L, N, 51,
     &                            LS,             PDOWN(L,I,J),
     &                            QQ(L,I,J),      ALPHA,
     &                            ALPHA2,         RAINFRAC,
     &                            WASHFRAC,       MASS_WASH,
     &                            MASS_NOWASH,    WETLOSS,
     &                            GAINED,         LOST,
     &                            DSTT(NN,:,I,J), STT(I,J,:,N) )
                  ENDIF
               ENDDO
               DO N=IDTCTSEA,IDTCTSEA
                  STT(I,J,L,N) = STT(I,J,L,N)
     &                         + SO2GAINED(I,J,L)
     &                         * FCLOUD(I,J,L,(NSO4+4))
                  ! Negative tracer...call subroutine SAFETY
                  IF ( STT(I,J,L,N) < 0d0 ) THEN
                     CALL SAFETY( I, J, L, N, 51,
     &                            LS,             PDOWN(L,I,J),
     &                            QQ(L,I,J),      ALPHA,
     &                            ALPHA2,         RAINFRAC,
     &                            WASHFRAC,       MASS_WASH,
     &                            MASS_NOWASH,    WETLOSS,
     &                            GAINED,         LOST,
     &                            DSTT(NN,:,I,J), STT(I,J,:,N) )
                  ENDIF
               ENDDO

               ENDIF

               ! Save FTOP for next level
               FTOP = F   

            !===========================================================
            ! (6)  N o   D o w n w a r d   P r e c i p i t a t i o n 
            !
            ! If there is no precipitation leaving grid box (I,J,L), 
            ! then  set F, the effective area of precipitation in grid 
            ! box (I,J,L), to zero.
            !
            ! Also, all of the previously rained-out tracer that is now 
            ! coming down from grid box (I,J,L+1) will evaporate and 
            ! re-enter the atmosphere in the gas phase in grid box 
            ! (I,J,L).  This is called "resuspension".
            !===========================================================
            ELSE IF ( ABS( PDOWN(L,I,J) ) < 1d-30 ) THEN

               ! No precipitation at grid box (I,J,L), thus F = 0
               F = 0d0

               DO N = IDTSO4BIN1,(IDTSO4BIN1+NSO4-1)
                  STT(I,J,L,N) = STT(I,J,L,N)
     &                         - SO2WETLOSS(I,J,L)
     &                         * FCLOUD(I,J,L,(N-IDTSO4BIN1+1))
                  ! Negative tracer...call subroutine SAFETY
                  IF ( STT(I,J,L,N) < 0d0 ) THEN
                     CALL SAFETY( I, J, L, N, 6,
     &                            LS,             PDOWN(L,I,J),
     &                            QQ(L,I,J),      ALPHA,
     &                            ALPHA2,         RAINFRAC,
     &                            WASHFRAC,       MASS_WASH,
     &                            MASS_NOWASH,    WETLOSS,
     &                            GAINED,         LOST,
     &                            DSTT(NN,:,I,J), STT(I,J,:,N) )
                  ENDIF
               ENDDO
               IF(NCTSO4>1)THEN
               DO N=IDTCTSO4,(IDTCTSO4+NCTSO4-1)
                  STT(I,J,L,N) = STT(I,J,L,N)
     &                         - SO2WETLOSS(I,J,L)
     &                         * SUM(FCLOUD(I,J,L,1:NSO4))
                  ! Negative tracer...call subroutine SAFETY
                  IF ( STT(I,J,L,N) < 0d0 ) THEN
                     CALL SAFETY( I, J, L, N, 6,
     &                            LS,             PDOWN(L,I,J),
     &                            QQ(L,I,J),      ALPHA,
     &                            ALPHA2,         RAINFRAC,
     &                            WASHFRAC,       MASS_WASH,
     &                            MASS_NOWASH,    WETLOSS,
     &                            GAINED,         LOST,
     &                            DSTT(NN,:,I,J), STT(I,J,:,N) )
                  ENDIF
               ENDDO
               ENDIF
               DO N = IDTCTBCOC,IDTCTBCOC
                  STT(I,J,L,N) = STT(I,J,L,N)
     &                         - SO2WETLOSS(I,J,L)
     &                         * FCLOUD(I,J,L,(NSO4+1))
                  ! Negative tracer...call subroutine SAFETY
                  IF ( STT(I,J,L,N) < 0d0 ) THEN
                     CALL SAFETY( I, J, L, N, 6,
     &                            LS,             PDOWN(L,I,J),
     &                            QQ(L,I,J),      ALPHA,
     &                            ALPHA2,         RAINFRAC,
     &                            WASHFRAC,       MASS_WASH,
     &                            MASS_NOWASH,    WETLOSS,
     &                            GAINED,         LOST,
     &                            DSTT(NN,:,I,J), STT(I,J,:,N) )
                  ENDIF
               ENDDO
               DO N = (IDTCTBCOC+1),(IDTCTBCOC+1)
                  STT(I,J,L,N) = STT(I,J,L,N)
     &                         - SO2WETLOSS(I,J,L)
     &                         * FCLOUD(I,J,L,(NSO4+2))
                  ! Negative tracer...call subroutine SAFETY
                  IF ( STT(I,J,L,N) < 0d0 ) THEN
                     CALL SAFETY( I, J, L, N, 6,
     &                            LS,             PDOWN(L,I,J),
     &                            QQ(L,I,J),      ALPHA,
     &                            ALPHA2,         RAINFRAC,
     &                            WASHFRAC,       MASS_WASH,
     &                            MASS_NOWASH,    WETLOSS,
     &                            GAINED,         LOST,
     &                            DSTT(NN,:,I,J), STT(I,J,:,N) )
                  ENDIF
               ENDDO
               DO N=IDTCTDST,IDTCTDST
                  STT(I,J,L,N) = STT(I,J,L,N)
     &                         - SO2WETLOSS(I,J,L)
     &                         * FCLOUD(I,J,L,(NSO4+3))
                  ! Negative tracer...call subroutine SAFETY
                  IF ( STT(I,J,L,N) < 0d0 ) THEN
                     CALL SAFETY( I, J, L, N, 6,
     &                            LS,             PDOWN(L,I,J),
     &                            QQ(L,I,J),      ALPHA,
     &                            ALPHA2,         RAINFRAC,
     &                            WASHFRAC,       MASS_WASH,
     &                            MASS_NOWASH,    WETLOSS,
     &                            GAINED,         LOST,
     &                            DSTT(NN,:,I,J), STT(I,J,:,N) )
                  ENDIF
               ENDDO
               DO N=IDTCTSEA,IDTCTSEA
                  STT(I,J,L,N) = STT(I,J,L,N)
     &                         - SO2WETLOSS(I,J,L)
     &                         * FCLOUD(I,J,L,(NSO4+4))
                  ! Negative tracer...call subroutine SAFETY
                  IF ( STT(I,J,L,N) < 0d0 ) THEN
                     CALL SAFETY( I, J, L, N, 6,
     &                            LS,             PDOWN(L,I,J),
     &                            QQ(L,I,J),      ALPHA,
     &                            ALPHA2,         RAINFRAC,
     &                            WASHFRAC,       MASS_WASH,
     &                            MASS_NOWASH,    WETLOSS,
     &                            GAINED,         LOST,
     &                            DSTT(NN,:,I,J), STT(I,J,:,N) )
                  ENDIF
               ENDDO

               ! Loop over soluble tracers and/or aerosol tracers
               DO NN = 1, NSOL
                  N = IDWETD(NN)
                  if(n<idtdstbin1.or.n>=(idtdstbin1+ndstb))then

                  IF(.NOT.((N>=IDTCTSO4.and.N<(IDTCTSO4+NCTSO4)).or.
     &               (N>=IDTCTSEA.and.N<(IDTCTSEA+NCTSEA))))THEN
                  ! WETLOSS is the amount of tracer in grid box (I,J,L) 
                  ! that is lost to rainout. (qli, bmy, 10/29/02)
                  WETLOSS = -DSTT(NN,L+1,I,J)

                  STT(I,J,L,N) = STT(I,J,L,N) - WETLOSS

                  ! There is nothing rained out/washed out in grid box
                  ! (I,J,L), so set DSTT at grid box (I,J,L) to zero.
                  DSTT(NN,L,I,J) = 0d0
                  
                  ! Negative tracer...call subroutine SAFETY
                  IF ( STT(I,J,L,N) < 0d0 ) THEN
                     CALL SAFETY( I, J, L, N, 61, 
     &                            LS,             PDOWN(L,I,J), 
     &                            QQ(L,I,J),      ALPHA,      
     &                            ALPHA2,         RAINFRAC,    
     &                            WASHFRAC,       MASS_WASH,    
     &                            MASS_NOWASH,    WETLOSS,        
     &                            GAINED,         LOST,        
     &                            DSTT(NN,:,I,J), STT(I,J,:,N) )
                  ENDIF

                  ENDIF

                  endif
               ENDDO

               ! Loop over soluble tracers and/or aerosol tracers
               DO NN = 1, NCTSO4
                  N = IDTCTSO4+NN-1
                    STT(I,J,L,N) = STT(I,J,L,N) + 
     &              DSTT((N-IDTSO4G+1),L+1,I,J)
                    DSTT((N-IDTSO4G+1),L,I,J) = 0d0
                  ! Negative tracer...call subroutine SAFETY
                  IF ( STT(I,J,L,N) < 0d0 ) THEN
                     CALL SAFETY( I, J, L, N, 61,
     &                            LS,             PDOWN(L,I,J),
     &                            QQ(L,I,J),      ALPHA,
     &                            ALPHA2,         RAINFRAC,
     &                            WASHFRAC,       MASS_WASH,
     &                            MASS_NOWASH,    WETLOSS,
     &                            GAINED,         LOST,
     &                            DSTT(NN,:,I,J), STT(I,J,:,N) )
                  ENDIF
               ENDDO
               DO NN = 1, NCTSEA
                  N = IDTCTSEA+NN-1
                    STT(I,J,L,N) = STT(I,J,L,N) +
     &              DSTT((N-IDTSO4G+1),L+1,I,J)
                    DSTT((N-IDTSO4G+1),L,I,J) = 0d0
                  ! Negative tracer...call subroutine SAFETY
                  IF ( STT(I,J,L,N) < 0d0 ) THEN
                     CALL SAFETY( I, J, L, N, 61,
     &                            LS,             PDOWN(L,I,J),
     &                            QQ(L,I,J),      ALPHA,
     &                            ALPHA2,         RAINFRAC,
     &                            WASHFRAC,       MASS_WASH,
     &                            MASS_NOWASH,    WETLOSS,
     &                            GAINED,         LOST,
     &                            DSTT(NN,:,I,J), STT(I,J,:,N) )
                  ENDIF
               ENDDO

               ! Save FTOP for next level
               FTOP = F
            ENDIF
         ENDDO

         !==============================================================
         ! (7)  W a s h o u t   i n   L e v e l   1
         !
         ! Assume all of the tracer precipitating down from grid box 
         ! (I,J,L=2) to grid box (I,J,L=1) gets washed out in grid box 
         ! (I,J,L=1).
         !==============================================================

         ! Zero variables for this level
         ALPHA       = 0d0
         ALPHA2      = 0d0
         F           = 0d0
         F_PRIME     = 0d0
         GAINED      = 0d0
         K_RAIN      = 0d0
         LOST        = 0d0
         MASS_NOWASH = 0d0
         MASS_WASH   = 0d0
         Q           = 0d0
         QDOWN       = 0d0
         RAINFRAC    = 0d0
         WASHFRAC    = 0d0
         WETLOSS     = 0d0
         
         ! We are at the surface, set L = 1
         L = 1

         IF ( PDOWN(L+1,I,J) > 0d0 ) THEN

            ! QDOWN is the precip leaving thru the bottom of box (I,J,L+1)
            QDOWN = PDOWN(L+1,I,J)

            ! Since no precipitation is forming within grid box (I,J,L),
            ! F' = 0, and F = MAX( F', FTOP ) reduces to F = FTOP.
            F = FTOP

            ! Only compute washout if F > 0.
            ! This helps to eliminate unnecessary CPU cycles.
            IF ( F > 0d0 ) THEN

               MASS1=SUM(STT(I,J,L,IDTSO4BIN1:(IDTSO4BIN1+NSO4-1)))
               MASS3=SUM(STT(I,J,L,IDTSEABIN1:(IDTSEABIN1+NSEA-1)))

               ! Loop over soluble tracers and/or aerosol tracers
               DO NN = 1, NSOL
                  N = IDWETD(NN)
                  if(n<idtdstbin1.or.n>=(idtdstbin1+ndstb))then

                  IF(.NOT.((N>=IDTCTSO4.and.N<(IDTCTSO4+NCTSO4)).or.
     &               (N>=IDTCTSEA.and.N<(IDTCTSEA+NCTSEA))))THEN
                  ! Call WASHOUT to compute the fraction of tracer 
                  ! in grid box (I,J,L) that is lost to washout.  
                  CALL WASHOUT( I,     J,  L, N, 
     &                          QDOWN, DT, F, WASHFRAC, AER,
     &                          State_Met )

                  IF( N >= IDTSO4BIN1 .and.
     &                N < (IDTSO4BIN1+NSO4) )THEN
                    IF ( State_Met%T(I,J,L) > 268.D0 ) THEN
                    RSIZEIN=GFTOT3D(I,J,L,1)*RDRY(N-IDTSO4BIN1+1)
                    CALL SIZEWASHFRAC( RSIZEIN, QDOWN, WASHFRAC )
                    WASHFRAC=(1.D0-EXP(-WASHFRAC*DT))
                    ELSE
                    WASHFRAC=0.D0
                    ENDIF
                  ENDIF

                  IF( N >= IDTSEABIN1 .and.
     &                N < (IDTSEABIN1+NSEA) )THEN
                    IF ( State_Met%T(I,J,L) > 268.D0 ) THEN
                    RSIZEIN=GFTOT3D(I,J,L,2)*RSALT(N-IDTSEABIN1+1)
                    CALL SIZEWASHFRAC( RSIZEIN, QDOWN, WASHFRAC )
                    WASHFRAC=(1.D0-EXP(-WASHFRAC*DT))
                    ELSE
                    WASHFRAC=0.D0
                    ENDIF
                  ENDIF

                  ! NOTE: for HNO3 and aerosols, there is an F factor
                  ! already present in WASHFRAC.  For other soluble
                  ! gases, we need to multiply by the F (hyl, bmy, 10/27/00)
                  IF ( AER ) THEN
                     WETLOSS = STT(I,J,L,N) * WASHFRAC
                  ELSE
                     WETLOSS = STT(I,J,L,N) * WASHFRAC * F
                  ENDIF

                  ! Subtract WETLOSS from STT
                  STT(I,J,L,N) = STT(I,J,L,N) - WETLOSS     
              
                  !-----------------------------------------------------
                  ! Dirty kludge to prevent wet deposition from removing 
                  ! stuff from stratospheric boxes -- this can cause 
                  ! negative tracer (rvm, bmy, 6/21/00)
                  !
                  IF ( STT(I,J,L,N) < 0d0 .and. L > 23 ) THEN
                      WRITE ( 6, 101 ) I, J, L, N, 7
 101                  FORMAT( 'WETDEPBIN - STT < 0 at ', 3i4, 
     &                        ' for tracer ', i4, 'in area ', i4 )
                      PRINT*, 'STT:', STT(I,J,:,N)
                      STT(I,J,L,N) = 0d0
                  ENDIF
                  !-----------------------------------------------------

                  ! Negative tracer...call subroutine SAFETY
                  IF ( STT(I,J,L,N) < 0d0 ) THEN
                     CALL SAFETY( I, J, L, N, 7, 
     &                            LS,             PDOWN(L,I,J), 
     &                            QQ(L,I,J),      ALPHA,           
     &                            ALPHA2,         RAINFRAC,    
     &                            WASHFRAC,       MASS_WASH,    
     &                            MASS_NOWASH,    WETLOSS,    
     &                            GAINED,         LOST,        
     &                            DSTT(NN,:,I,J), STT(I,J,:,N) )
                  ENDIF

               ENDIF

               endif
               ENDDO

               MASS2=SUM(STT(I,J,L,IDTSO4BIN1:(IDTSO4BIN1+NSO4-1)))
               MASS4=SUM(STT(I,J,L,IDTSEABIN1:(IDTSEABIN1+NSEA-1)))

               ! Loop over soluble tracers and/or aerosol tracers
               DO NN = 1, NCTSO4
                  N = IDTCTSO4+NN-1
                  IF(MASS1>1.D-30)THEN
                    STT(I,J,L,N)=STT(I,J,L,N)*(MASS2/MASS1)
                  ENDIF
                  ! Negative tracer...call subroutine SAFETY
                  IF ( STT(I,J,L,N) < 0d0 ) THEN
                     CALL SAFETY( I, J, L, N, 7,
     &                            LS,             PDOWN(L,I,J),
     &                            QQ(L,I,J),      ALPHA,
     &                            ALPHA2,         RAINFRAC,
     &                            WASHFRAC,       MASS_WASH,
     &                            MASS_NOWASH,    WETLOSS,
     &                            GAINED,         LOST,
     &                            DSTT(NN,:,I,J), STT(I,J,:,N) )
                  ENDIF
               ENDDO
               DO NN = 1, NCTSEA
                  N = IDTCTSEA+NN-1
                  IF(MASS3>1.D-30)THEN
                    STT(I,J,L,N)=STT(I,J,L,N)*(MASS4/MASS3)
                  ENDIF
                  ! Negative tracer...call subroutine SAFETY
                  IF ( STT(I,J,L,N) < 0d0 ) THEN
                     CALL SAFETY( I, J, L, N, 7,
     &                            LS,             PDOWN(L,I,J),
     &                            QQ(L,I,J),      ALPHA,
     &                            ALPHA2,         RAINFRAC,
     &                            WASHFRAC,       MASS_WASH,
     &                            MASS_NOWASH,    WETLOSS,
     &                            GAINED,         LOST,
     &                            DSTT(NN,:,I,J), STT(I,J,:,N) )
                  ENDIF
               ENDDO

            ENDIF    
         ENDIF

      ENDDO
      ENDDO
#if   !defined( SGI_MIPS )
!$OMP END PARALLEL DO
#endif

      ! Clear the pointer
      NULLIFY( STT )

      ! Return to calling program
      END SUBROUTINE WETDEPBIN

!------------------------------------------------------------------------------

      FUNCTION LS_K_RAIN( Q ) RESULT( K_RAIN )
!
!******************************************************************************
!  Function LS_K_RAIN computes K_RAIN, the first order rainout rate constant 
!  for large-scale (a.k.a. stratiform) precipitation (bmy, 3/18/04)
! 
!  Arguments as Input:
!  ============================================================================
!  (1 ) Q      (REAL*8) : Rate of precip formation [cm3 H2O/cm3 air/s]
!
!  Function value:
!  ============================================================================
!  (2 ) K_RAIN (REAL*8) : 1st order rainout rate constant [s-1]
!
!  NOTES:
!  (1 ) Now made into a MODULE routine since we cannot call internal routines
!        from w/in a parallel loop.  Updated comments. (bmy, 3/18/04)
!******************************************************************************
!
      ! Arguments
      REAL*8, INTENT(IN) :: Q
      
      ! Function value
      REAL*8             :: K_RAIN

      !==================================================================
      ! LS_K_RAIN begins here!
      !==================================================================

      ! Compute rainout rate constant K in s^-1 (Eq. 12, Jacob et al, 2000).
      ! 1.0d-4 = K_MIN, a minimum value for K_RAIN 
      ! 1.5d-6 = L + W, the condensed water content (liq + ice) in the cloud
      K_RAIN = 1.0d-4 + ( Q / 1.5d-6 ) 
      
      ! Return to WETDEPBIN
      END FUNCTION LS_K_RAIN

!------------------------------------------------------------------------------

      FUNCTION LS_F_PRIME( Q, K_RAIN ) RESULT( F_PRIME )
!
!******************************************************************************
!  Function LS_F_PRIME computes F', the fraction of the grid box that is 
!  precipitating during large scale (a.k.a. stratiform) precipitation.
!  (bmy, 3/18/04)
! 
!  Arguments as Input:
!  ============================================================================
!  (1 ) Q       (REAL*8) : Rate of precip formation [cm3 H2O/cm3 air/s]
!  (2 ) K_RAIN  (REAL*8) : 1st order rainout rate constant [s-1]
!
!  Function value:
!  ============================================================================
!  (3 ) F_PRIME (REAL*8) : Fraction of grid box undergoing LS precip [unitless]
!
!  NOTES:
!  (1 ) Now made into a MODULE routine since we cannot call internal routines
!        from w/in a parallel loop.  Updated comments. (bmy, 3/18/04)
!******************************************************************************
!
      ! Arguments
      REAL*8, INTENT(IN) :: Q, K_RAIN

      ! Function value
      REAL*8             :: F_PRIME

      !=================================================================
      ! LS_F_PRIME begins here!
      !=================================================================

      ! Compute F', the area of the grid box undergoing precipitation
      ! 1.5d-6 = L + W, the condensed water content [cm3 H2O/cm3 air]
      F_PRIME = Q / ( K_RAIN * 1.5d-6 )

      ! Return to WETDEPBIN
      END FUNCTION LS_F_PRIME

!------------------------------------------------------------------------------

      FUNCTION CONV_F_PRIME( Q, K_RAIN, DT ) RESULT( F_PRIME )
!
!******************************************************************************
!  Function CONV_F_PRIME computes F', the fraction of the grid box that is 
!  precipitating during convective precipitation. (bmy, 3/18/04)
! 
!  Arguments as Input:
!  ============================================================================
!  (1 ) Q       (REAL*8) : Rate of precip formation        [cm3 H2O/cm3 air/s]
!  (2 ) K_RAIN  (REAL*8) : 1st order rainout rate constant [s-1]
!  (3 ) DT      (REAL*8) : Wet deposition timestep         [s]
!
!  Function value:
!  ============================================================================
!  (4 ) F_PRIME (REAL*8) : Frac. of grid box undergoing CONV precip [unitless]
!
!  NOTES:
!  (1 ) Now made into a MODULE routine since we cannot call internal routines
!        from w/in a parallel loop.  Updated comments. (bmy, 3/18/04)
!******************************************************************************
!
      ! Arguments
      REAL*8, INTENT(IN) :: Q, K_RAIN, DT
      
      ! Local variables
      REAL*8             :: TIME

      ! Function value
      REAL*8             :: F_PRIME

      !=================================================================
      ! CONV_F_PRIME begins here!
      !=================================================================
      
      ! Assume the rainout event happens in 30 minutes (1800 s)
      ! Compute the minimum of DT / 1800s and 1.0
      TIME = MIN( DT / 1800d0, 1d0 )

      ! Compute F' for convective precipitation (Eq. 13, Jacob et al, 2000)
      ! 0.3  = FMAX, the maximum value of F' for convective precip
      ! 2d-6 = L + W, the condensed water content [cm3 H2O/cm3 air]
      F_PRIME = ( 0.3d0 * Q * TIME ) / 
     &          ( ( Q * TIME ) + ( 0.3d0 * K_RAIN * 2d-6 ) )

      ! Return to WETDEPBIN
      END FUNCTION CONV_F_PRIME

!------------------------------------------------------------------------------

      SUBROUTINE SAFETY( I,         J,           L,        N,     
     &                   A,         LS,          PDOWN,    QQ,       
     &                   ALPHA,     ALPHA2,      RAINFRAC, WASHFRAC, 
     &                   MASS_WASH, MASS_NOWASH, WETLOSS,  GAINED,   
     &                   LOST,      DSTT,        STT )
!
!******************************************************************************
!  Subroutine SAFETY stops the run with debug output and an error message 
!  if negative tracers are found. (bmy, 3/18/04)
! 
!  Arguments as Input:
!  ============================================================================
!  (1 ) Q       (REAL*8) : Rate of precip formation        [cm3 H2O/cm3 air/s]
!  (2 ) K_RAIN  (REAL*8) : 1st order rainout rate constant [s-1]
!  (3 ) DT      (REAL*8) : Wet deposition timestep         [s]
!
!  Function value:
!  ============================================================================
!  (4 ) F_PRIME (REAL*8) : Frac. of grid box undergoing CONV precip [unitless]
!
!  NOTES:
!  (1 ) Now made into a MODULE routine since we cannot call internal routines
!        from w/in a parallel loop.  Updated comments. (bmy, 3/18/04)
!******************************************************************************
!
      ! References to F90 modules
      USE ERROR_MOD, ONLY : GEOS_CHEM_STOP

      USE CMN_SIZE_MOD    ! Size parameters

      ! Arguments
      LOGICAL, INTENT(IN) :: LS
      INTEGER, INTENT(IN) :: I, J, L, N, A
      REAL*8,  INTENT(IN) :: PDOWN,    QQ,       ALPHA,     ALPHA2
      REAL*8,  INTENT(IN) :: RAINFRAC, WASHFRAC, MASS_WASH, MASS_NOWASH
      REAL*8,  INTENT(IN) :: WETLOSS,  GAINED,   LOST,      DSTT(LLPAR)
      REAL*8,  INTENT(IN) :: STT(LLPAR)

      !=================================================================
      ! SAFETY begins here!
      !=================================================================
      
      ! Print line
      WRITE( 6, '(a)' ) REPEAT( '=', 79 )

      ! Write error message and stop the run
      WRITE ( 6, 100 ) I, J, L, N, A
 100  FORMAT( 'WETDEPBIN - STT < 0 at ', 3i4, ' for tracer ', i4, 
     &        ' in area ', i4 )

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
      PRINT*, 'DSTT(NN,:)  : ', DSTT(:)
      PRINT*, 'STT(I,J,:N) : ', STT(:)

      ! Print line
      WRITE( 6, '(a)' ) REPEAT( '=', 79 )

      ! Deallocate memory and stop
      CALL GEOS_CHEM_STOP

      ! Return to WETDEPBIN
      END SUBROUTINE SAFETY

!------------------------------------------------------------------------------

      SUBROUTINE SIZEWASHFRAC( RIN, PP, WASHRATE )
      REAL*8, INTENT(IN) :: RIN, PP
      REAL*8, INTENT(INOUT) :: WASHRATE

      INTEGER            :: N
      REAL*8             :: RINUM,PAR1,PAR2,PAR3,RATE
      REAL*8             :: PHOUR

       integer,parameter :: resol=100 ! number of differential scav
coefs used
       real*8,parameter,dimension(resol) :: 
     & A0   = (/
     & 1.2421339D-03,  1.1542652D-03,  1.0697728D-03,  9.8944578D-04,
     & 9.1371929D-04,  8.4278590D-04,  7.7667052D-04,  7.1528548D-04,
     & 6.5846901D-04,  6.0600966D-04,  5.5766889D-04,  5.1319298D-04,
     & 4.7232405D-04,  4.3480799D-04,  4.0039657D-04,  3.6885168D-04,
     & 3.3994846D-04,  3.1347378D-04,  2.8923032D-04,  2.6703325D-04,
     & 2.4671193D-04,  2.2810873D-04,  2.1107841D-04,  1.9548827D-04,
     & 1.8121649D-04,  1.6815161D-04,  1.5619207D-04,  1.4524521D-04,
     & 1.3522667D-04,  1.2605986D-04,  1.1767508D-04,  1.1000906D-04,
     & 1.0300445D-04,  9.6609200D-05,  9.0776336D-05,  8.6302263D-05,
     & 8.1774526D-05,  7.7650854D-05,  7.3907784D-05,  7.0525143D-05,
     & 6.7486259D-05,  6.4778164D-05,  6.2440762D-05,  6.0625356D-05,
     & 5.8598430D-05,  5.7528951D-05,  5.6584556D-05,  5.5783209D-05,
     & 5.5543436D-05,  5.5543436D-05,  5.5543436D-05,  5.7475162D-05,
     & 5.8920602D-05,  6.1073444D-05,  6.4497544D-05,  6.8131062D-05,
     & 7.2558585D-05,  7.7855845D-05,  8.4144327D-05,  9.1530451D-05,
     & 1.0014767D-04,  1.1015226D-04,  1.2173151D-04,  1.3511098D-04,
     & 1.3144316D+00,  2.1080957D+01,  6.4889491D+01,  1.9500575D+02,
     & 7.2548248D+02,  9.6941400D+02,  1.1577594D+03,  1.3120109D+03,
     & 1.4584424D+03,  1.5450455D+03,  1.6161885D+03,  1.6828528D+03,
     & 1.7346022D+03,  2.3845894D+03,  2.6876754D+03,  2.7387250D+03,
     & 2.7742642D+03,  2.8069166D+03,  2.8392922D+03,  2.7104748D+03,
     & 2.0122139D+03,  1.4870979D+03,  1.1876294D+03,  8.0640117D+02,
     & 3.2910188D+02,  7.1894308D+01,  2.1724393D+00,  1.3292166D-01,
     & 1.2247745D-01,  1.1591321D-01,  1.1317314D-01,  1.1448342D-01,
     & 1.2051768D-01,  1.3249278D-01,  1.5257908D-01,  1.8452159D-01 /)

       real*8,parameter,dimension(resol) ::
     & A1   = (/
     & 5.0580353D-03,  4.7109454D-03,  4.4035710D-03,  4.1284255D-03,
     & 3.8799211D-03,  3.6537793D-03,  3.4466705D-03,  3.2559540D-03,
     & 3.0795037D-03,  2.9155923D-03,  2.7627973D-03,  2.6199358D-03,
     & 2.4860212D-03,  2.3602081D-03,  2.2417832D-03,  2.1301353D-03,
     & 2.0247372D-03,  1.9251300D-03,  1.8309154D-03,  1.7417413D-03,
     & 1.6573005D-03,  1.5773210D-03,  1.5015609D-03,  1.4298067D-03,
     & 1.3618641D-03,  1.2975605D-03,  1.2367403D-03,  1.1792643D-03,
     & 1.1250047D-03,  1.0738467D-03,  1.0256848D-03,  9.8042322D-04,
     & 9.3797448D-04,  8.9825859D-04,  8.6120289D-04,  8.2707145D-04,
     & 7.9140109D-04,  7.5859975D-04,  7.2858162D-04,  7.0127807D-04,
     & 6.7664220D-04,  6.5465171D-04,  6.3491029D-04,  6.1917999D-04,
     & 6.0180213D-04,  5.9237683D-04,  5.8410135D-04,  5.7713998D-04,
     & 5.7521646D-04,  5.7521646D-04,  5.7521646D-04,  5.9573825D-04,
     & 6.1128457D-04,  6.3433325D-04,  6.7077811D-04,  7.0902428D-04,
     & 7.5496154D-04,  8.0921016D-04,  8.7210055D-04,  9.4417636D-04,
     & 1.0260809D-03,  1.1184538D-03,  1.2219490D-03,  1.3372226D-03,
     & 3.3669307D-06,  6.8382524D-07,  4.0853339D-07,  2.0483149D-07,
     & 7.4624009D-08,  7.0603692D-08,  7.1071154D-08,  7.2576727D-08,
     & 7.3354584D-08,  7.6011076D-08,  7.8321576D-08,  7.9914898D-08,
     & 8.1448083D-08,  6.1700605D-08,  5.6634320D-08,  5.7221333D-08,
     & 5.7962467D-08,  5.8669217D-08,  5.9349052D-08,  6.3637204D-08,
     & 8.7872021D-08,  1.2220386D-07,  1.5786685D-07,  2.4108266D-07,
     & 6.1647029D-07,  2.9680970D-06,  1.0439819D-04,  1.8745902D-03,
     & 2.2149189D-03,  2.5806136D-03,  2.9535205D-03,  3.3070670D-03,
     & 3.6059673D-03,  3.8132840D-03,  3.8954592D-03,  3.8302728D-03 /)

       real*8,parameter,dimension(resol) ::
     & A2   = (/
     & 6.6424127D-01,  6.6576376D-01,  6.6713995D-01,  6.6839236D-01,
     & 6.6953812D-01,  6.7059130D-01,  6.7156306D-01,  6.7246240D-01,
     & 6.7329683D-01,  6.7407304D-01,  6.7479637D-01,  6.7547196D-01,
     & 6.7610345D-01,  6.7669458D-01,  6.7724858D-01,  6.7776851D-01,
     & 6.7825605D-01,  6.7871442D-01,  6.7914446D-01,  6.7954857D-01,
     & 6.7992783D-01,  6.8028357D-01,  6.8061731D-01,  6.8092889D-01,
     & 6.8121941D-01,  6.8148927D-01,  6.8173855D-01,  6.8196686D-01,
     & 6.8217380D-01,  6.8235808D-01,  6.8251825D-01,  6.8265207D-01,
     & 6.8275644D-01,  6.8282760D-01,  6.8286013D-01,  6.8282405D-01,
     & 6.8277795D-01,  6.8266725D-01,  6.8247922D-01,  6.8219803D-01,
     & 6.8180393D-01,  6.8127402D-01,  6.8056216D-01,  6.7977710D-01,
     & 6.7843974D-01,  6.7732850D-01,  6.7585803D-01,  6.7349287D-01,
     & 6.7115454D-01,  6.7115454D-01,  6.7115454D-01,  6.6188326D-01,
     & 6.5857402D-01,  6.5475758D-01,  6.5009930D-01,  6.4629666D-01,
     & 6.4264392D-01,  6.3921057D-01,  6.3599722D-01,  6.3311369D-01,
     & 6.3051789D-01,  6.2819509D-01,  6.2612153D-01,  6.2426768D-01,
     & 4.9996833D-01,  6.0034315D-01,  6.6171289D-01,  6.9929362D-01,
     & 7.2455931D-01,  7.4315159D-01,  7.5769135D-01,  7.6946300D-01,
     & 7.7915407D-01,  7.8717177D-01,  7.9378675D-01,  7.9919713D-01,
     & 8.0356028D-01,  8.0700554D-01,  8.0963902D-01,  8.1154737D-01,
     & 8.1280391D-01,  8.1342701D-01,  8.1346568D-01,  8.1291308D-01,
     & 8.1175022D-01,  8.0993808D-01,  8.0742008D-01,  8.0412217D-01,
     & 7.9996035D-01,  7.9483587D-01,  7.8800543D-01,  7.6975396D-01,
     & 7.6003299D-01,  7.4935955D-01,  7.3794038D-01,  7.2605591D-01,
     & 7.1404544D-01,  7.0226337D-01,  6.9104364D-01,  6.8066128D-01 /)

       real*8,parameter,dimension(resol) ::
     & radresol = (/
     & 1.1220191D-03,  1.2589269D-03,  1.4125401D-03,  1.5848971D-03,
     & 1.7782848D-03,  1.9952694D-03,  2.2387304D-03,  2.5118983D-03,
     & 2.8183979D-03,  3.1622965D-03,  3.5481569D-03,  3.9810999D-03,
     & 4.4668699D-03,  5.0119134D-03,  5.6234631D-03,  6.3096331D-03,
     & 7.0795286D-03,  7.9433667D-03,  8.9126099D-03,  1.0000118D-02,
     & 1.1220323D-02,  1.2589417D-02,  1.4125566D-02,  1.5849154D-02,
     & 1.7783055D-02,  1.9952927D-02,  2.2387566D-02,  2.5119277D-02,
     & 2.8184310D-02,  3.1623334D-02,  3.5481986D-02,  3.9811466D-02,
     & 4.4669226D-02,  5.0119724D-02,  5.6235291D-02,  6.3097067D-02,
     & 7.0796117D-02,  7.9434596D-02,  8.9127131D-02,  1.0000235D-01,
     & 1.1220455D-01,  1.2589565D-01,  1.4125732D-01,  1.5849341D-01,
     & 1.7783263D-01,  1.9953163D-01,  2.2387829D-01,  2.5119573D-01,
     & 2.8184637D-01,  3.1623703D-01,  3.5482404D-01,  3.9811930D-01,
     & 4.4669750D-01,  5.0120312D-01,  5.6235945D-01,  6.3097805D-01,
     & 7.0796943D-01,  7.9435527D-01,  8.9128178D-01,  1.0000352D+00,
     & 1.1220586D+00,  1.2589712D+00,  1.4125898D+00,  1.5849527D+00,
     & 1.7783473D+00,  1.9953395D+00,  2.2388091D+00,  2.5119867D+00,
     & 2.8184969D+00,  3.1624076D+00,  3.5482817D+00,  3.9812400D+00,
     & 4.4670272D+00,  5.0120902D+00,  5.6236610D+00,  6.3098550D+00,
     & 7.0797777D+00,  7.9436460D+00,  8.9129219D+00,  1.0000469D+01,
     & 1.1220717D+01,  1.2589860D+01,  1.4126063D+01,  1.5849712D+01,
     & 1.7783680D+01,  1.9953630D+01,  2.2388355D+01,  2.5120161D+01,
     & 2.8185301D+01,  3.1624447D+01,  3.5483231D+01,  3.9812866D+01,
     & 4.4670795D+01,  5.0121487D+01,  5.6237267D+01,  6.3099289D+01,
     & 7.0798607D+01,  7.9437386D+01,  8.9130272D+01,  1.0000587D+02 /)


      !=================================================================
      ! WASHFRAC begins here!
      !=================================================================
      !PP cm/s
      PHOUR=PP*3.6d6 !mm/hour

      RINUM=RIN*1.D6 !m to um

      DO N=2,resol
      IF(RINUM<=radresol(1))THEN
        PAR1=A0(1)
        PAR2=A1(1)
        PAR3=A2(1)
      ELSE IF(RINUM>radresol(resol))THEN
        PAR1=A0(resol)
        PAR2=A1(resol)
        PAR3=A2(resol)
      ELSE IF((RINUM>radresol(N-1)).AND.(RINUM<=radresol(N)))THEN
        RATE=(RINUM-radresol(N-1))/(radresol(N)-radresol(N-1))
        PAR1=A0(N-1)+RATE*(A0(N)-A0(N-1))
        PAR2=A1(N-1)+RATE*(A1(N)-A1(N-1))
        PAR3=A2(N-1)+RATE*(A2(N)-A2(N-1))
      ENDIF
      ENDDO

      WASHRATE=PAR1*(EXP(PAR2*(PHOUR**PAR3))-1.D0) !s-1

      END SUBROUTINE SIZEWASHFRAC
!------------------------------------------------------------------------------

      SUBROUTINE WETDEPBINID
!
!******************************************************************************
!  Subroutine WETDEPBINID sets up the index array of soluble tracers used in
!  the WETDEPBIN routine above (bmy, 11/8/02, 5/18/06)
! 
!  NOTES:
!  (1 ) Now references "tracerid_mod.f".  Also references "CMN" in order to
!        pass variables NSRCX and NTRACE. (bmy, 11/8/02)
!  (2 ) Updated for carbon aerosol & dust tracers (rjp, bmy, 4/5/04)
!  (3 ) Updated for seasalt aerosol tracers.  Also added fancy output.
!        (rjp, bec, bmy, 4/20/04)
!  (4 ) Updated for secondary organic aerosol tracers (bmy, 7/13/04)
!  (5 ) Now references N_TRACERS, TRACER_NAME, TRACER_MW_KG from
!        "tracer_mod.f".  Removed reference to NSRCX.  (bmy, 7/20/04)
!  (6 ) Updated for mercury aerosol tracers (eck, bmy, 12/9/04)
!  (7 ) Updated for AS, AHS, LET, NH4aq, SO4aq (cas, bmy, 12/20/04)
!  (8 ) Updated for SO4s, NITs (bec, bmy, 4/25/05)
!  (9 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (10) Now use IS_Hg2 and IS_HgP to determine if a tracer is a tagged Hg2
!        or HgP tracer (bmy, 1/6/06)
!  (11) Now added SOG4 and SOA4 (dkh, bmy, 5/18/06)
!******************************************************************************
!
      ! References To F90 modules
      USE ERROR_MOD,    ONLY : ERROR_STOP
      USE TRACER_MOD,   ONLY : N_TRACERS, N_APMTRA

      USE CMN_SIZE_MOD        ! Size parameters

      ! Local variables
      INTEGER :: N, NN

      !=================================================================
      ! WETDEPBINID begins here!
      !=================================================================

      ! Zero NSOL
      NSOL = 0

      ! Sort soluble tracers into IDWETD
      DO N = N_TRACERS+1, N_TRACERS+N_APMTRA

         !-----------------------------
         ! aerosize tracers
         !-----------------------------
         NSOL         = NSOL + 1
         IDWETD(NSOL) = N

      ENDDO

      ! Error check: Make sure that NSOL is less than NSOLMAX
      IF ( NSOL > NSOLMAX ) THEN
         WRITE(*,*)NSOL,'>',NSOLMAX
         CALL ERROR_STOP( 'NSOL > NSOLMAX!',
     &        'WETDEPBINID (wetscavbin_mod.f)')
      ENDIF
      
      ! Return to calling program
      END SUBROUTINE WETDEPBINID

!------------------------------------------------------------------------------

      FUNCTION GET_WETDEPBIN_IDWETD( NWET ) RESULT( N )
!
!******************************************************************************
!  Function GET_WETDEPBIN_IDWETD returns the tracer number of wet deposition 
!  species  NWET.  This is meant to be called outside of WETSCAVBIN_MOD so that 
!  IDWETD can be kept as a PRIVATE variable. (bmy, 1/10/03)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) NWET (INTEGER) : Wet deposition species N
!
!  NOTES:
!******************************************************************************
!
      ! References to F90 modules
      USE ERROR_MOD, ONLY : ERROR_STOP

      ! Arguments
      INTEGER, INTENT(IN) :: NWET
      
      ! Function value
      INTEGER             :: N

      !=================================================================
      ! GET_WETDEPBIN_IDWETD begins here!
      !=================================================================

      ! Make sure NWET is valid
      IF ( NWET < 1 .or. NWET > NSOLMAX ) THEN
         CALL ERROR_STOP( 'Invalid value of NWET!', 
     &                    'GET_N_WETDEPBIN (wetscav_mod.f)' )
      ENDIF

      ! Get the tracer # for wet deposition species N
      N = IDWETD(NWET)
     
      ! Return to calling program
      END FUNCTION GET_WETDEPBIN_IDWETD

!------------------------------------------------------------------------------

      SUBROUTINE INIT_WETSCAVBIN( State_Met )
!
!******************************************************************************
!  Subroutine INIT_WETSCAVBIN initializes updraft velocity, cloud liquid water
!  content, cloud ice content, and mixing ratio of water fields, which
!  are used in the wet scavenging routines. (bmy, 2/23/00, 3/7/05)
!
!  NOTES:
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
!        values from for the first call to COMPUTE_F .  Now call WETDEPBINID
!        on the first call to initialize the wetdep index array. (bmy, 1/27/03)
!  (8 ) Now references STT from "tracer_mod.f".  Also now we call WETDEPBINID
!        from "input_mod.f" (bmy, 7/20/04)
!  (9 ) Now references new function E_ICE, which is an analytic function of 
!        Kelvin temperature instead of Celsius. (bmy, 3/7/05)
!  09 Nov 2012 - M. Payer    - Replaced all met field arrays with State_Met
!                              derived type object
!******************************************************************************
!
      ! References to F90 modules
      USE DAO_MOD,            ONLY : IS_WATER
      USE ERROR_MOD,          ONLY : ALLOC_ERR
      USE GIGC_State_Met_Mod, ONLY : MetState
      USE PRESSURE_MOD,       ONLY : GET_PCENTER

      USE CMN_SIZE_MOD             ! Size parameters

      ! Arguments
      TYPE(MetState), INTENT(IN)  :: State_Met   ! Meteorology State object

      ! Local variables
      INTEGER             :: I, J, L, AS
      REAL*8              :: PL, TK
      LOGICAL, SAVE       :: FIRST = .TRUE.      
            
      !=================================================================
      ! INIT_WETSCAVBIN begins here!
      !=================================================================
      IF ( FIRST ) THEN

         ! Allocate Vud on first call
         ALLOCATE( Vud( IIPAR, JJPAR ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'Vud' )
         Vud = 0d0

         ! Allocate C_H2O on first call
         ALLOCATE( C_H2O( IIPAR, JJPAR, LLPAR ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'C_H2O' )
         C_H2O = 0d0

         ! Allocate CLDLIQ on first call
         ALLOCATE( CLDLIQ( IIPAR, JJPAR, LLPAR ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'CLDLIQ' )
         CLDLIQ = 0d0

         ! Allocate CLDICE on first call
         ALLOCATE( CLDICE( IIPAR, JJPAR, LLPAR ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'CLDICE' )
         CLDICE = 0d0

         ! Reset flag
         FIRST = .FALSE. 
      ENDIF

      !=================================================================
      ! Compute Vud, CLDLIQ, CLDICE, C_H2O, following Jacob et al, 2000.
      !=================================================================
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, TK, PL )
!$OMP+SCHEDULE( DYNAMIC )
      DO L = 1, LLPAR
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! Compute Temp [K] and Pressure [hPa]
         TK = State_Met%T(I,J,L)
         PL = GET_PCENTER(I,J,L)

         !==============================================================
         ! Compute Vud -- 5 m/s over oceans, 10 m/s over land (or ice?)
         ! Assume Vud is the same at all altitudes; the array can be 2-D
         !==============================================================
         IF ( L == 1 ) THEN
            IF ( IS_WATER( I, J, State_Met ) ) THEN
               Vud(I,J) = 5d0
            ELSE
               Vud(I,J) = 10d0
            ENDIF
         ENDIF

         !==============================================================
         ! CLDLIQ, the cloud liquid water content [cm3 H2O/cm3 air], 
         ! is a function of the local Kelvin temperature:
         !  
         !    CLDLIQ = 2e-6                    [     T >= 268 K    ]
         !    CLDLIQ = 2e-6 * ((T - 248) / 20) [ 248 K < T < 268 K ]
         !    CLDLIQ = 0                       [     T <= 248 K    ]
         !==============================================================
         IF ( TK >= 268d0 ) THEN
            CLDLIQ(I,J,L) = 2d-6

         ELSE IF ( TK > 248d0 .and. TK < 268d0 ) THEN
            CLDLIQ(I,J,L) = 2d-6 * ( ( TK - 248d0 ) / 20d0 )

         ELSE
            CLDLIQ(I,J,L) = 0d0
            
         ENDIF
           
         !=============================================================
         ! CLDICE, the cloud ice content [cm3 ice/cm3 air] is given by:
         !
         !    CLDICE = 2e-6 - CLDLIQ
         !=============================================================
         CLDICE(I,J,L) = 2d-6 - CLDLIQ(I,J,L)

         !=============================================================
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
         !==============================================================
         C_H2O(I,J,L) = E_ICE( TK ) / PL

      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      ! Return to calling program
      END SUBROUTINE INIT_WETSCAVBIN

!------------------------------------------------------------------------------

      SUBROUTINE CLEANUP_WETSCAVBIN

      !=================================================================
      ! Subroutine CLEANUP_WETSCAVBIN deallocates arrays for 
      ! wet scavenging / wet deposition
      !=================================================================
      IF ( ALLOCATED( Vud    ) ) DEALLOCATE( Vud    )
      IF ( ALLOCATED( C_H2O  ) ) DEALLOCATE( C_H2O  )
      IF ( ALLOCATED( CLDLIQ ) ) DEALLOCATE( CLDLIQ )
      IF ( ALLOCATED( CLDICE ) ) DEALLOCATE( CLDICE )
      IF ( ALLOCATED( PDOWN  ) ) DEALLOCATE( PDOWN  )
      IF ( ALLOCATED( QQ     ) ) DEALLOCATE( QQ     )

      ! Return to calling program
      END SUBROUTINE CLEANUP_WETSCAVBIN

!-----------------------------------------------------------------------------

      ! End of module
      END MODULE APM_WETS_MOD
