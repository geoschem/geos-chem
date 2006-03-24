! $Id: mercury_mod.f,v 1.9 2006/03/24 20:22:54 bmy Exp $
      MODULE MERCURY_MOD
!
!******************************************************************************
!  Module MERCURY_MOD contains variables and routines for the GEOS-CHEM 
!  mercury simulation. (eck, bmy, 12/14/04, 3/16/06)
!
!  Module Variables:
!  ============================================================================
!  (1 ) AN_Hg0   (INTEGER) : Tracer index array for tagged anth Hg(0)  regions
!  (2 ) AN_Hg2   (INTEGER) : Tracer index array for tagged anth Hg(II) regions
!  (3 ) AN_HgP   (INTEGER) : Tracer index array for tagged anth HgP    regions
!  (4 ) COSZM    (REAL*8 ) : Max daily value of COS( S. Z. Angle ) [unitless]
!  (5 ) DRYHg2   (INTEGER) : Index for Hg2 in DEPSAV array (drydep freq)
!  (6 ) DRYHgP   (INTEGER) : Index for HgP in DEPSAV array (drydep freq)
!  (7 ) EHg0_an  (REAL*8 ) : Anthropogenic Hg0 emissions [kg/s]
!  (8 ) EHg2_an  (REAL*8 ) : Anthropogenic Hg2 emissions [kg/s]
!  (9 ) EHgP_an  (REAL*8 ) : Anthropogenic HgP emissions [kg/s]
!  (10) EHg0_oc  (REAL*8 ) : Hg0 emissions from oceans [kg/s]
!  (11) EHg0_ln  (REAL*8 ) : Re-emissions of Hg0 from land [kg/s]
!  (12) EHg0_nt  (REAL*8 ) : Hg0 emissions from natural land sources [kg/s] 
!  (13) TCOSZ    (REAL*8 ) : Sum of COS( Solar Zenith Angle ) [unitless]
!  (14) TTDAY    (REAL*8 ) : Total daylight time at location (I,J) [minutes]
!  (15) T44      (REAL*4 ) : Local array for drydep diagnostic
!  (15) N_HgTAGS (INTEGER) : Number of tagged sources (1 or 8)
!  (16) ZERO_DVEL(REAL*8 ) : Array with zero dry deposition velocity
!  (17) ANTHRO_Hg_YEAR(INT): Anthropogenic Hg emissions year (1995 or 2000)
!
!  Module Routines:
!  ============================================================================
!  (1 ) CHEMMERCURY        : Chemistry routine for Hg
!  (2 ) CHEM_Hg0_Hg2       : Chemistry for Hg0, Hg2 and drydep of Hg2
!  (3 ) RXN_Hg0_Hg2        : Conversion of Hg(0) --> Hg(II) via reduction rxn
!  (4 ) RXN_Hg0            : Conversion of Hg(0) --> Hg(II) via oxidation rxns
!  (5 ) RXN_Hg2_DRYD       : Prod of Hg(II) from Hg(0); also Hg(II) drydep loss
!  (6 ) RXN_Hg2            : Prod of Hg(II) from Hg(0) 
!  (7 ) CHEM_HGP           : Chemistry (via drydep loss) for HgP
!  (8 ) RXN_HgP_DRYD       : Loss of HgP via drydep
!  (9 ) EMISSMERCURY       : Emissions for Hg
!  (10) SRCHG0             : Applies emissions of Hg0
!  (11) SRCHG2             : Applies emissions of Hg2
!  (12) SRCHGP             : Applies emissions of HgP
!  (13) MERCURY_READYR     : Reads mercury emissions and converts to [kg/s]
!  (14) GET_LWC            : Computes liquid water content as a function of T
!  (15) GET_VCLDF          : Computes volume cloud fraction 
!  (16) GET_O3             : Returns monthly mean O3 field
!  (17) GET_OH             : Returns monthly mean OH field (w/ diurnal scaling)
!  (18) OHNO3TIME          : Computes diurnal scaling for monthly mean OH
!  (19) DEFINE_TAGGED_Hg   : Defines tracer number for tagged Hg tracers 
!  (20) INIT_MERCURY       : Allocates and zeroes module arrays
!  (21) CLEANUP_MERCURY    : Deallocates module arrays
!
!  GEOS-CHEM modules referenced by mercury_mod.f
!  ============================================================================
!  (1 ) bpch2_mod.f        : Module w/ routines for binary pch file I/O
!  (2 ) comode_mod.f       : Module w/ SMVGEAR allocatable arrays
!  (3 ) dao_mod.f          : Module w/ arrays for DAO met fields
!  (4 ) diag_mod.f         : Module w/ GEOS-CHEM diagnostic arrays
!  (5 ) drydep_mod.f       : Module w/ GEOS-CHEM dry deposition routines
!  (6 ) error_mod.f        : Module w/ NaN, other error check routines
!  (7 ) global_o3_mod.f    : Module w/ routines to read 3-D O3 field
!  (8 ) global_oh_mod.f    : Module w/ routines to read 3-D OH field
!  (9 ) grid_mod.f         : Module w/ horizontal grid information
!  (10) logical_mod.f      : Module w/ GEOS-CHEM logical switches
!  (11) pbl_mix_mod.f      : Module w/ routines for PBL height & mixing
!  (12) pressure_mod.f     : Module w/ routines to compute P(I,J,L)
!  (13) regrid_1x1_mod.f   : Module w/ routines to regrid 1x1 data
!  (13) time_mod.f         : Module w/ routines to compute date & time
!  (14) tracer_mod.f       : Module w/ GEOS-CHEM tracer array STT etc.
!  (15) tracerid_mod.f     : Module w/ pointers to tracers & emissions
!  (16) transfer_mod.f     : Module w/ routines to cast & resize arrays
!
!  Nomenclature: 
!  ============================================================================
!  (1 ) Hg(0)  a.k.a. Hg0  : Elemental   mercury
!  (2 ) Hg(II) a.k.a. Hg2  : Divalent    mercury
!  (3 ) HgP                : Particulate mercury
!
!  Mercury Tracers (1-3 are always defined; 4-21 are defined for tagged runs)
!  ============================================================================
!  (1 ) Hg(0)              : Hg(0)  - total tracer
!  (2 ) Hg(II)             : Hg(II) - total tracer 
!  (3 ) HgP                : HgP    - total tracer
!  ------------------------+---------------------------------------------------
!  (4 ) Hg(0)_an           : Hg(0)  - North American anthropogenic
!  (5 ) Hg(0)_ae           : Hg(0)  - European Anthropogenic
!  (6 ) Hg(0)_aa           : Hg(0)  - Asian anthropogenic
!  (7 ) Hg(0)_ar           : Hg(0)  - Rest of World Anthropogenic
!  (8 ) Hg(0)_oc           : Hg(0)  - Ocean emission
!  (9 ) Hg(0)_ln           : Hg(0)  - Land reemission
!  (10) Hg(0)_nt           : Hg(0)  - Land natural emission
!  ------------------------+---------------------------------------------------
!  (11) Hg(II)_an          : Hg(II) - North American anthropogenic
!  (12) Hg(II)_ae          : Hg(II) - European Anthropogenic
!  (13) Hg(II)_aa          : Hg(II) - Asian anthropogenic
!  (14) Hg(II)_ar          : Hg(II) - Rest of World Anthropogenic
!  (15) Hg(II)_oc          : Hg(II) - Ocean emission
!  (16) Hg(II)_ln          : Hg(II) - Land reemission
!  (17) Hg(II)_nt          : Hg(II) - Land natural emission
!  ------------------------+---------------------------------------------------
!  (18) HgP_an             : HgP    - North American anthropogenic
!  (19) HgP_ae             : HgP    - European anthropogenic
!  (20) HgP_aa             : HgP    - Asian anthropogenic
!  (21) HgP_ar             : HgP    - Rest of world anthropogenic
!  ------------------------+---------------------------------------------------
!  (22) HgP_oc             : HgP    - Ocean emission        (FOR FUTURE USE)
!  (23) HgP_ln             : HgP    - Land reemission       (FOR FUTURE USE)
!  (24) HgP_nt             : HgP    - Land natural emission (FOR FUTURE USE)
!
!  References:
!  ============================================================================
!  (1 ) Hall, B. (1995). "The gas phase oxidation of elemental mercury by 
!        ozone.", Water, Air, and Soil Pollution 80: 301-315.
!  (2 ) Sommar, J., et al. (2001). "A kinetic study of the gas-phase 
!        reaction between the hydroxyl radical and atomic mercury." 
!        Atmospheric Environment 35: 3049-3054.
!
!  NOTES:
!  (1 ) Updated for reduction rxn and online Hg0 ocean flux.  Now use 
!        diagnostic arrays from "diag03_mod.f".  (eck, sas, bmy, 1/21/05)
!  (2 ) Now references "pbl_mix_mod.f".  Remove FEMIS array and routine
!        COMPUTE_FEMIS. (bmy, 2/15/05)
!  (3 ) Now can read data for both GEOS and GCAP grids (bmy, 8/16/05)
!  (4 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (5 ) Now references XNUMOL from "tracer_mod.f" (bmy, 10/25/05)
!  (6 ) Various updates added for tagged Hg sim. (eck, sas, cdh, bmy, 3/16/06)
!******************************************************************************
!
      IMPLICIT NONE

      !=================================================================
      ! MODULE PRIVATE DECLARATIONS -- keep certain internal variables 
      ! and routines from being seen outside "mercury_mod.f"
      !=================================================================

      ! Make everything PRIVATE ...
      PRIVATE

      ! ... except these routines
      PUBLIC :: CHEMMERCURY
      PUBLIC :: CLEANUP_MERCURY
      PUBLIC :: INIT_MERCURY
      PUBLIC :: EMISSMERCURY

      !=================================================================
      ! MODULE VARIABLES
      !=================================================================

      ! Scalars 
      !-----------------------------------------------
      ! Prior to 2/24/06:
      ! This is now in "logical_mod.f" (bmy, 2/24/06)
      !LOGICAL              :: LDYNOCEAN
      !-----------------------------------------------
      INTEGER              :: ANTHRO_Hg_YEAR
      INTEGER              :: DRYHg2
      INTEGER              :: DRYHgP

      ! Parameters
      REAL*8,  PARAMETER   :: SMALLNUM = 1d-32

      ! Arrays
      INTEGER, ALLOCATABLE :: AN_Hg0(:,:)
      INTEGER, ALLOCATABLE :: AN_Hg2(:,:)
      INTEGER, ALLOCATABLE :: AN_HgP(:,:)
      REAL*8,  ALLOCATABLE :: COSZM(:,:)
      REAL*8,  ALLOCATABLE :: EHg0_an(:,:)
      REAL*8,  ALLOCATABLE :: EHg2_an(:,:)
      REAL*8,  ALLOCATABLE :: EHgP_an(:,:)
      REAL*8,  ALLOCATABLE :: EHg0_oc(:,:,:)
      REAL*8,  ALLOCATABLE :: EHg0_ln(:,:)
      REAL*8,  ALLOCATABLE :: EHg0_nt(:,:)
      REAL*4,  ALLOCATABLE :: T44(:,:,:,:)
      REAL*8,  ALLOCATABLE :: TCOSZ(:,:)
      REAL*8,  ALLOCATABLE :: TTDAY(:,:)
      REAL*8,  ALLOCATABLE :: ZERO_DVEL(:,:)

      !=================================================================
      ! MODULE ROUTINES -- follow below the "CONTAINS" statement
      !=================================================================
      CONTAINS

!------------------------------------------------------------------------------

      SUBROUTINE CHEMMERCURY
!
!******************************************************************************
!  Subroutine CHEMMERCURY is the driver routine for mercury chemistry
!  in the GEOS-CHEM module. (eck, bmy, 12/6/04, 1/9/06)
!
!  NOTES:
!  (1 ) Now references routine GET_PBL_MAX_L from "pbl_mix_mod.f".  Now
!        references AD44 from "diag_mod.f".  Now sum the levels from T44 into 
!        the AD44 array.  Now references N_TRACERS from "tracer_mod.f".
!        (bmy, 2/24/05)
!  (2 ) Bug fix: Set T44 to 0e0 for single precision.  Now allow for zero
!        dry deposition velocity.  Now call INIT_MERCURY from "input_mod.f"
!        (bmy, 2/24/06)
!******************************************************************************
!
      ! References to F90 modules
      USE DIAG_MOD,      ONLY : AD44
      USE DRYDEP_MOD,    ONLY : DEPSAV
      USE ERROR_MOD,     ONLY : DEBUG_MSG
      USE GLOBAL_O3_MOD, ONLY : GET_GLOBAL_O3
      USE GLOBAL_OH_MOD, ONLY : GET_GLOBAL_OH
      USE PBL_MIX_MOD,   ONLY : GET_PBL_MAX_L
      USE LOGICAL_MOD,   ONLY : LPRT
      USE TIME_MOD,      ONLY : GET_MONTH, ITS_A_NEW_MONTH
      USE TRACER_MOD,    ONLY : N_TRACERS

#     include "CMN_SIZE"      ! Size parameters
#     include "CMN_DIAG"      ! ND44

      ! Local variables
      LOGICAL, SAVE          :: FIRST = .TRUE.
      INTEGER                :: I, J, L, MONTH, N, PBL_MAX

      !=================================================================
      ! CHEMMERCURY begins here!
      !=================================================================

      !-----------------------------------------------------------------------
      ! Prior to 2/24/06:
      ! Now call INIT_MERCURY from "input_mod.f" (bmy, 2/24/06)
      !! First-time initialization
      !IF ( FIRST ) THEN
      !
      !   ! Allocate arrays & locate drydep indices (if not already done)
      !   !CALL INIT_MERCURY
      !
      !   ! Reset first-time flag
      !   FIRST = .FALSE.
      !ENDIF
      !-----------------------------------------------------------------------

      !=================================================================
      ! Read monthly mean OH and O3 fields
      !=================================================================
      IF ( ITS_A_NEW_MONTH() ) THEN 

         ! Get the current month
         MONTH = GET_MONTH()

         ! Read monthly mean OH and O3 from disk
         CALL GET_GLOBAL_OH( MONTH )
         IF ( LPRT ) CALL DEBUG_MSG( '### CHEMMERC: a GET_GLOBAL_OH' )

         CALL GET_GLOBAL_O3( MONTH )
         IF ( LPRT ) CALL DEBUG_MSG( '### CHEMMERC: a GET_GLOBAL_O3' )

      ENDIF
      
      !=================================================================
      ! Perform chemistry on Hg tracers 
      !=================================================================
      
      ! Compute diurnal scaling for OH
      CALL OHNO3TIME
      IF ( LPRT ) CALL DEBUG_MSG( 'CHEMMERCURY: a OHNO3TIME' )

      ! Zero dry deposition tmp array
      IF ( ND44 > 0 ) THEN
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, N )
         DO N = 1, N_TRACERS
         DO L = 1, LLTROP
         DO J = 1, JJPAR
         DO I = 1, IIPAR
            T44(I,J,L,N) = 0e0
         ENDDO
         ENDDO
         ENDDO
         ENDDO
!$OMP END PARALLEL DO
      ENDIF

      !-------------------------
      ! Hg0 and Hg2 chemistry
      !-------------------------
      IF ( DRYHg2 > 0 ) THEN
         
         ! If DRYHG2 > 0 then drydep is active; 
         ! pass drydep frequency to CHEM_Hg0_Hg2
         CALL CHEM_Hg0_Hg2( DEPSAV(:,:,DRYHg2) )
         IF ( LPRT ) CALL DEBUG_MSG( 'CHEMMERCURY: a CHEM_Hg0_Hg2' )
      
      ELSE

         ! Otherwise pass zero drydep frequency
         CALL CHEM_Hg0_Hg2( ZERO_DVEL )
         IF ( LPRT ) CALL DEBUG_MSG( 'CHEMMERCURY: a CHEM_Hg0_Hg2' )
         
      ENDIF
         
      !--------------------------
      ! HgP chemistry
      !--------------------------
      IF ( DRYHgP > 0 ) THEN

         ! If DRYHgP > 0, then drydep is active;
         ! Pass drydep frequency to CHEM_HgP
         CALL CHEM_HgP( DEPSAV(:,:,DRYHgP) )
         IF ( LPRT ) CALL DEBUG_MSG( 'CHEMMERCURY: a CHEM_HgP' )

      ELSE
         
         ! Otherwise pass zero drydep frequency
         CALL CHEM_HgP( ZERO_DVEL )
         IF ( LPRT ) CALL DEBUG_MSG( 'CHEMMERCURY: a CHEM_HgP' )

      ENDIF

      ! Archive drydep fluxes into ND44
      IF ( ND44 > 0 ) THEN
         
         ! Model layers where the PBL top occurs
         PBL_MAX = GET_PBL_MAX_L()

         ! Sum levels into AD44 
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, N )
         DO N = 1, N_TRACERS
         DO J = 1, JJPAR
         DO I = 1, IIPAR
         DO L = 1, PBL_MAX
            AD44(I,J,N,1) = AD44(I,J,N,1) + T44(I,J,L,N)
         ENDDO
         ENDDO
         ENDDO
         ENDDO
!$OMP END PARALLEL DO
      ENDIF

      ! Return to calling program
      END SUBROUTINE CHEMMERCURY

!------------------------------------------------------------------------------

      SUBROUTINE CHEM_Hg0_Hg2( V_DEP_Hg2 )
!
!******************************************************************************
!  Subroutine CHEM_Hg0_Hg2 is the chemistry subroutine for the oxidation/
!  reduction of Hg0/Hg(II).  For tracers with dry deposition, the loss rate 
!  of dry dep is combined in the chemistry loss term. 
!  (eck, bmy, 12/6/04, 1/9/06)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) V_DEP_Hg2 (REAL*8) : Dry deposition velocity for Hg(II) [cm/s]
!
!  Description of the chemistry mechanism:
!  ============================================================================
!  (1 ) Conversion from Hg(0) to Hg(II): Oxidation by O3 and OH
!
!         Hg(0)(g)+ O3(g) --> Hg(II) ,  k  = 3.0e-20 cm3 molec-1 s-1
!                                       Source: Hall, 1995
!       
!         Hg(0)(g)+ OH(g) --> Hg(II) ,  k  = 8.7e-14 cm3 s-1 
!                                       Source: Sommar et al. 2001
!           
!  (2 ) Aqueous-phase photochemical reduction of Hg(II) is included based
!        on estimate of rate constant and scaled to OH concentrations.
! 
!  (3 ) Hg(II) is dry-deposited,        kd = Dvel/DELZ (sec-1)        
!     
!  NOTES:
!  (1 ) Updated for reduction reaction.  Now use diagnostic arrays from
!        "diag03_mod.f" (eck, bmy, 1/21/05)
!  (2 ) Now references GET_FRAC_UNDER_PBLTOP from "pbl_mix_mod.f".  Now
!        performs drydep for all levels in the PBL.  Changed Kred to 2.1e-10
!        on advice of Noelle Eckley Selin. (bmy, 2/24/05)
!  (3 ) Now references XNUMOL from "tracer_mod.f" (bmy, 10/25/05)
!  (4 ) Now prevent divide-by-zero error.  Now use ID_Hg0 and ID_Hg2 index
!        arrays from "tracerid_mod.f".  Also modified for updated ocean
!        mercury module.  Updated some constants. (eck, cdh, sas, bmy, 1/9/06) 
!******************************************************************************
!
      ! References to F90 modules
      USE DAO_MOD,      ONLY : T, AD
      USE DIAG03_MOD,   ONLY : AD03_Hg2_Hg0, AD03_Hg2_O3
      USE DIAG03_MOD,   ONLY : AD03_Hg2_OH,  LD03,  ND03
      USE LOGICAL_MOD,  ONLY : LSPLIT
      USE PBL_MIX_MOD,  ONLY : GET_FRAC_UNDER_PBLTOP
      USE TIME_MOD,     ONLY : GET_TS_CHEM
      USE TRACER_MOD,   ONLY : STT,    XNUMOL
      USE TRACERID_MOD, ONLY : ID_Hg0, ID_Hg2, ID_Hg_tot, N_Hg_CATS

#     include "CMN_SIZE"     ! Size parameters

      ! Arguments
      REAL*8, INTENT(IN)    :: V_DEP_Hg2(IIPAR,JJPAR)

      ! Local variables
      INTEGER               :: I, J, L, N, NN
      REAL*8                :: AREA_CM2,    C_O3,        C_OH
      REAL*8                :: DRYDEP,      DTCHEM,      E_RKT
      REAL*8                :: E_K1,        E_K2,        FA
      REAL*8                :: FC,          K1,          K2          
      REAL*8                :: Ko,          Kr,          Kt
      REAL*8                :: Kr_Kt,       NK1,         NK2
      REAL*8                :: LOST_Hg0_O3, LOST_Hg0_OH, LWC
      REAL*8                :: F_UNDER_TOP, RKT,         E_SALT
      REAL*8                :: LOST_Hg0(N_Hg_CATS)
      REAL*8                :: Ox_Hg0(N_Hg_CATS)

      ! K for reaction Hg0 + O3 in cm3 molec-1 s-1 (Source: Hall '95)
      REAL*8, PARAMETER     :: K       = 3.0d-20

      ! K for reaction Hg2 + OH in cm3 s-1 
      !-------------------------------------------------
      ! Prior to 1/9/06:
      ! Update value (eck, bmy, 3/16/06)
      !REAL*8, PARAMETER     :: K_HG_OH = 8.7d-14
      !-------------------------------------------------
      REAL*8, PARAMETER     :: K_HG_OH = 7.8d-14

      ! Henry's Law constant for Hg2
      REAL*8, PARAMETER     :: HL      = 1.4d6

      ! Gas constant??
      REAL*8, PARAMETER     :: R       = 8.2d-2

      ! K for reduction (scaled to budget and OH conc)
      ! eck may change this later
      !-------------------------------------------------
      ! Prior to 1/9/06:
      ! Update value (eck, bmy, 3/16/06)
      !REAL*8, PARAMETER     :: Kred    = 2.1d-10
      !-------------------------------------------------
      REAL*8, PARAMETER     :: Kred    = 8.4d-10
     
      ! K for sea salt (eck, bmy, 1/9/06)
      REAL*8, PARAMETER     :: KSalt   = 3.8d-5

      ! External functions
      REAL*8,  EXTERNAL     :: BOXVL

      !=================================================================
      ! CHEM_Hg0_Hg2 begins here!
      !=================================================================

      ! Chemistry timestep [s]
      DTCHEM = GET_TS_CHEM() * 60d0

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I,           J,      L,           NN,          LOST_Hg0 )
!$OMP+PRIVATE( F_UNDER_TOP, C_O3,   C_OH,        FC,          LWC      )
!$OMP+PRIVATE( FA,          K1,     K2,          Ko,          Kr       )
!$OMP+PRIVATE( Kt,          Kr_Kt,  RKT,         E_RKT,       E_K1     )
!$OMP+PRIVATE( E_K2,        E_SALT, LOST_Hg0_O3, LOST_Hg0_OH, Ox_Hg0   )
      DO L = 1, LLPAR
      DO J = 1, JJPAR
      DO I = 1, IIPAR
        
         ! Initialize LOST_Hg0, Ox_Hg0 arrays
         DO NN = 1, N_Hg_CATS
            LOST_Hg0(NN) = 0d0
            Ox_Hg0(NN)   = 0d0
         ENDDO
         
         ! Also initialize for ND03 diag
         LOST_Hg0_O3 = 0d0
         LOST_Hg0_OH = 0d0

         ! Fraction of box (I,J,L) underneath the PBL top [unitless]
         F_UNDER_TOP = GET_FRAC_UNDER_PBLTOP( I, J, L )

         ! Monthly mean O3 and OH concentrations [molec/cm3]
         C_O3        = GET_O3( I, J, L )
         C_OH        = GET_OH( I, J, L )
                
         ! Get fraction of cloud in the gridbox [unitless]
         FC          = GET_VCLDF( I, J, L )

         ! Get liquid water content in cloudy area of gridbox (m3/m3)
         LWC         = GET_LWC( T(I,J,L) ) * FC 

         ! Define fraction of Hg(II) which is in air
         FA          = ( HL * R * T(I,J,L) * LWC )
         FA          = FA / ( 1d0 + FA )

         ! Define K's for the oxidation reactions
         K1          = K       * C_O3
         K2          = K_HG_OH * C_OH
         Ko          = K1      + K2
         
         ! Define K for the reduction reaction.  Include the fraction of 
         ! Hg(II) in air within the Kreduction & scale to OH concentration
         Kr          = Kred * FA * C_OH 
         
         ! Total rate constant: Oxidation + Reduction
         Kt          = Ko + Kr

         ! Ratio of Kr / Kt 
         IF ( Kt > SMALLNUM ) THEN
            Kr_Kt    = Kr / Kt
         ELSE
            Kr_Kt    = SMALLNUM
         ENDIF
         
         ! Kt * timestep (i.e. the argument for EXP) [unitless]
         RKT         = Kt * DTCHEM
        
         ! Compute the exponential terms for use below
         E_RKT       = EXP( -RKT            )
         E_K1        = EXP( -K1    * DTCHEM )
         E_K2        = EXP( -K2    * DTCHEM )
         E_SALT      = EXP( -Ksalt * DTCHEM )

         !==============================================================
         ! Hg(0) Chemistry: Conversion from Hg(0) to Hg(II)
         !
         ! CASE 1: REDUCTION REACTION
         ! -------------------------------------------------------------
         ! Aqueous chem occurs if Kr>0, it's cloudy, and T > -15 C.
         ! In this case we have the reaction:
         !
         !    New Hg(0) = {  ( Old Hg(0) * EXP( -Kt * DT ) )          } +
         !                {  ( Old Hg(0) + Old Hg(II)      ) * Kr/Kt * 
         !                   (         1 - EXP( -Kt * DT ) )          }    
         !
         ! where: Kr         = K of the reduction rxn
         !        Kt         = K of the total rxn (oxidation + reduction)
         !        DT         = Chemistry timestep [s]
         !        Old Hg(0)  = Hg(0)  at start of this timestep
         !        Old Hg(II) = Hg(II) as start of this timestep
         !
         ! The amount of Hg(0) lost in this rxn becomes Hg(II):
         !
         !    Hg(0) converted to Hg(II) = Old Hg(0) - New Hg(0)
         !
         ! CASE 2: OXIDATION REACTIONS
         ! -------------------------------------------------------------
         ! If aqueous chemistry does not happen, then Hg(0) is converted 
         ! to Hg(II) via the oxidation rxns with rate constants K1, K2:
         !
         !    New Hg(0) = Old_Hg(0) * EXP( -K1 * DT ) * EXP( -K2 * DT )
         !
         ! The amt of Hg(0) lost in these rxns rxn becomes Hg(II):
         !
         !    Hg(0) converted to Hg(II) = Old Hg(0) - New Hg(0)
         !==============================================================
         IF ( Kr > 0d0 .and. FC > 0d0 .and. T(I,J,L) > 258d0 ) THEN
                     
            !--------------------------------------------
            ! CASE 1: Total Hg(II) tracer, aq chem
            !--------------------------------------------
            
            ! Loop over all Hg0 tracers
            DO NN = 1, N_Hg_CATS

               ! Conversion of total Hg(0) --> Hg(II) 
!--------------------------------------------------------------------------
! Prior to 2/17/06:  
! cdh modification 
!               CALL RXN_Hg0_Hg2( I,          J,          L,     
!     &                           ID_Hg0(NN), ID_Hg2(NN), RKT, 
!     &                           E_RKT,      Kr_Kt,      LOST_Hg0(NN) )
!--------------------------------------------------------------------------
               CALL RXN_Hg0_Hg2( I,            J,          L,     
     &                           ID_Hg0(NN),   ID_Hg2(NN), RKT, 
     &                           E_RKT,        Kr_Kt,      Kt,
     &                           Ko,           Kr,         DTCHEM,
     &                           LOST_Hg0(NN), Ox_Hg0(NN) )

            ENDDO

         ELSE

            !--------------------------------------------
            ! CASE 2: Total Hg(II) tracer, non-aq chem
            !--------------------------------------------
            
            ! Loop over all Hg0 tracers
            DO NN = 1, N_Hg_CATS

               ! Total Hg(0) loss by oxidation rxn
               CALL RXN_Hg0( I,          J,    L, 
     &                       ID_Hg0(NN), E_K1, E_K2, LOST_Hg0(NN) )

               ! Gross oxidation flux Hg(0) -> Hg(II)
               Ox_Hg0(NN) = LOST_Hg0(NN)

            ENDDO

         ENDIF
            
         !==============================================================
         ! Compute Hg(II) production from OH and O3 rxns for diagnostic
         !==============================================================
         !------------------------------------------
         ! Prior to 2/16/06: 
         ! cdh modification
         !IF ( LOST_Hg0(ID_Hg_tot) > 0d0 ) THEN
         !------------------------------------------
         IF ( Ox_Hg0(ID_Hg_tot) > 0d0 ) THEN
         
            ! This is for the diagnostics OH and O3 prod of Hg(II).
            ! They are messed up a little since adding the reduction 
            ! reaction haven't fixed them yet. (eck, 12/7/04)
            !----------------------------------------------------------------
            ! Prior to 1/9/06:
            ! Rewrite code below to prevent underflow in ND03 diagnostic
            ! Also note LOST_Hg0_O3 and LOST_Hg0_OH are set to zero at the
            ! top of the loop and will remain so, unless reset below.
            ! (bmy, 1/9/06)
            !! Define new k's here to avoid NaN error in denominator.
            !! (eck, 12/7/04)
            !NK1         = MAX( K1, SMALLNUM )
            !NK2         = MAX( K2, SMALLNUM )
            !
            ! Production of Hg(II) from O3 rxn [kg]
            !LOST_Hg0_O3 = ( NK1 / ( NK1 + NK2 ) ) * LOST_Hg0(ID_Hg_tot)
            !
            ! Production of Hg(II) from OH rxn [kg]
            !LOST_Hg0_OH = ( NK2 / ( NK1 + NK2 ) ) * LOST_Hg0(ID_Hg_tot)
            !----------------------------------------------------------------

            ! Avoid division by zero
            IF ( ( K1 + K2 ) > 0d0 ) THEN

               ! Production of Hg(II) from O3 rxn [kg]
               !-------------------------------------------------------------
               ! Prior to 2/16/06:
               ! cdh modification
               !LOST_Hg0_O3 = ( K1 / ( K1 + K2 ) ) * LOST_Hg0(ID_Hg_tot)
               !-------------------------------------------------------------
               LOST_Hg0_O3 = ( K1 / ( K1 + K2 ) ) * Ox_Hg0(ID_Hg_tot)

               ! Production of Hg(II) from OH rxn [kg]
               !-------------------------------------------------------------
               ! Prior to 2/16/06:
               ! cdh modification
               !LOST_Hg0_OH = ( K2 / ( K1 + K2 ) ) * LOST_Hg0(ID_Hg_tot)
               !-------------------------------------------------------------
               LOST_Hg0_OH = ( K2 / ( K1 + K2 ) ) * Ox_Hg0(ID_Hg_tot)

            ENDIF

         ENDIF

         !==============================================================
         ! Hg(II) chemistry: Conversion from Hg(0) and drydep loss
         !
         ! CASE 1: WITHIN THE PLANETARY BOUNDARY LAYER (PBL)
         ! -------------------------------------------------------------
         ! At the surface we have both dry deposition of Hg(II) plus
         ! conversion of Hg(0) into Hg(II).  In this case we use the
         ! following rxns:
         !
         !    CASE 1a: If Conv Hg(0) > 0:
         !    ---------------------------
         !    New Hg(II) = ( Old Hg(II)     *       EXP( -RKT ) ) +
         !                 ( Conv Hg(0)/RKT * ( 1 - EXP( -RKT ) )
         !               
         !
         !    CASE 1b: If Conv Hg(0) <= 0:
         !    ----------------------------
         !    New Hg(II) = ( Old Hg(II) + Conv Hg(0) ) * EXP( -RKT )
         !
         ! Where: RKT         = DTCHEM * Drydep Vel of Hg(II)
         !       "Conv Hg(0)" = Amt of Hg(0) that became Hg(II), which
         !                      is archived from the Hg(0) rxns above
         !
         ! CASE 2: OUTSIDE THE PLANETARY BOUNDARY LAYER (PBL)
         ! -------------------------------------------------------------
         ! At levels higher than the surface, we do not have drydep.
         ! Therefore we only have conversion of Hg(0) into Hg(II), and
         ! we use this rxn:
         !
         !     New Hg(II) = Old Hg(II) + Conv Hg(0)
         !
         !
         ! NOTE: ND44 diagnostics are archived in RXN_Hg2_DRYD.
         !==============================================================

         ! If we are in the PBL and there is nonzero drydep velocity ...
         IF ( F_UNDER_TOP > 0d0 .and. V_DEP_Hg2(I,J) > 0d0 ) THEN
         
            ! Hg2 drydep frequency [1/s] -- F_UNDER_TOP accounts for the 
            ! fraction of box (I,J,L) that is located beneath the PBL top
            RKT   = V_DEP_Hg2(I,J) * DTCHEM * F_UNDER_TOP

            ! Pre-compute exponential term 
            E_RKT = EXP( -RKT )
               
            !-----------------------------------------------
            ! CASE 1: Total Hg(II) tracer, in PBL
            !-----------------------------------------------
            
            ! Loop over all Hg2 tracers
            DO NN = 1, N_Hg_CATS

               ! Compute new Hg(II) concentration in PBL
               CALL RXN_Hg2_DRYD( I,          J,            L,        
     &                            ID_Hg2(NN), RKT,          E_RKT, 
     &                            E_SALT,     LOST_Hg0(NN), DTCHEM )
          
            ENDDO

         ELSE

            !--------------------------------------------
            ! CASE 2: Total Hg(II) tracer, outside PBL
            !--------------------------------------------

            ! Loop over all Hg2 tracers
            DO NN = 1, N_Hg_CATS

               ! Compute new concentration of total Hg(II) outside PBL [kg]
               CALL RXN_Hg2( I, J, L, ID_Hg2(NN), LOST_Hg0(NN) )
            
            ENDDO

         ENDIF
     
         !==============================================================
         ! ND03 diagnostic: Hg(II) production [kg]
         !==============================================================
         IF ( ND03 > 0 .and. L <= LD03 ) THEN
            NN = ID_Hg_tot
            AD03_Hg2_Hg0(I,J,L) = AD03_Hg2_Hg0(I,J,L) + LOST_Hg0(NN)
            AD03_Hg2_OH(I,J,L)  = AD03_Hg2_OH(I,J,L)  + LOST_Hg0_OH
            AD03_Hg2_O3(I,J,L)  = AD03_Hg2_O3(I,J,L)  + LOST_Hg0_O3
         ENDIF

      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO
      
      ! Return to calling program
      END SUBROUTINE CHEM_Hg0_Hg2

!------------------------------------------------------------------------------

      SUBROUTINE RXN_Hg0_Hg2( I,   J,      L,        N_Hg0,  N_Hg2, 
     &                        RKT, E_RKT,  Kr_Kt,    Kt,     Ko, 
     &                        Kr,  DTCHEM, LOST_Hg0, Ox_Hg0 )
!
!******************************************************************************
!  Subroutine RXN_Hg0_Hg2 computes the conversion of Hg(0) to Hg(II) via
!  an aqueous chemistry reduction reaction.  The formula used below is a s
!  solution of the 2-box model equation. (eck, bmy, 12/14/04, 2/17/06)
! 
!  Arguments as Input:
!  ============================================================================
!  (1-3) I, J, L  (INTEGER) : GEOS-CHEM lon, lat, alt grid box indices
!  (4  ) N_Hg0    (INTEGER) : Index for Hg(0)  total or tagged tracers
!  (5  ) N_Hg2    (INTEGER) : Index for Hg(II) total or tagged tracers
!  (6  ) RKT      (REAL*8 ) : Value of R * k * Temp for this rxn [unitless]
!  (7  ) E_RKT    (REAL*8 ) : Value of EXP( - RKT ) [unitless]
!  (8  ) Kr_Kt    (REAL*8 ) : Ratio of Kr / Kt (redux K / total K) [unitless]
!
!  Arguments as Output:
!  ============================================================================
!  (9  ) LOST_Hg0 (REAL*8 ) : Loss term: Hg(0) before - Hg(0) after [kg]
!
!  NOTES:
!  (1  ) Changed equation to reflect reduction rxn.  Also modified to output
!         the gross oxidation flux. (eck, cdh, bmy, 2/17/06)
!******************************************************************************
!
      ! References to F90 modules
      USE TRACER_MOD, ONLY  : STT

      ! Arguments
      INTEGER, INTENT(IN)  :: I,        J,     L,     N_Hg0, N_Hg2
      REAL*8,  INTENT(IN)  :: RKT,      E_RKT, Kr_Kt, Kt,    Kr  
      REAL*8,  INTENT(IN)  :: Ko,       DTCHEM
      REAL*8,  INTENT(OUT) :: LOST_Hg0, Ox_Hg0

      ! Local variables
      REAL*8               :: OLD_Hg0, OLD_Hg2, NEW_Hg0

      !=================================================================
      ! RXN_Hg0_Hg2 begins here!
      !=================================================================

      ! Error check tracer number
      IF ( N_Hg0 < 1 .or. N_Hg2 < 1 ) RETURN

      ! Initial concentrations of Hg(0) and Hg(II) [kg]
      OLD_Hg0          = MAX( STT(I,J,L,N_Hg0), SMALLNUM )
      OLD_Hg2          = MAX( STT(I,J,L,N_Hg2), SMALLNUM )

      ! New concentration of Hg(0) [kg]
      NEW_Hg0 = (             OLD_Hg0   *                   E_RKT   ) + 
     &          ((( OLD_Hg2 + OLD_Hg0 ) * Kr_Kt ) * ( 1d0 - E_RKT ) )

      ! Gross oxidation flux Hg(0) -> Hg(II) 
      Ox_Hg0  = ( Kt - Kr ) / Kt**2 *
     &        ( ( OLD_Hg0 + OLD_Hg2 ) * Kr *
     &          ( E_RKT - 1d0 + Kt * DTCHEM )
     &          + OLD_Hg0 * Kt * ( 1d0 -E_RKT ) )
               
      ! Set a floor of zero
      OX_Hg0  = MAX( OX_Hg0, 0d0 )

      ! Save back into STT array [kg]
      NEW_Hg0          = MAX( NEW_Hg0, SMALLNUM )
      STT(I,J,L,N_Hg0) = NEW_Hg0 
      
      ! Compute amount of Hg(0) which has now become Hg(II) [kg]
      LOST_Hg0         = OLD_Hg0 - NEW_Hg0

      ! Return to calling program
      END SUBROUTINE RXN_Hg0_Hg2
      
!------------------------------------------------------------------------------

      SUBROUTINE RXN_Hg0( I, J, L, N, E_K1, E_K2, LOST_Hg0 )
!
!******************************************************************************
!  Subroutine RXN_Hg0 computes the loss of Hg(0) by oxidation reactions with
!  rate constants K1 and K2. (eck, bmy, 12/14/04)
! 
!  Arguments as Input:
!  ============================================================================
!  (1-3) I, J, L  (INTEGER) : GEOS-CHEM lon, lat, alt grid box indices
!  (4  ) N        (INTEGER) : Index for Hg(0) total or tagged tracers
!  (5  ) E_K1     (REAL*8 ) : Value of EXP( -K1 ) [unitless]
!  (6  ) E_K2     (REAL*8 ) : Value of EXP( -K2 ) [unitless]
!
!  Arguments as Output:
!  ============================================================================
!  (7  ) LOST_Hg0 (REAL*8 ) : Loss term: Hg(0) before - Hg(0) after [kg]
!
!  NOTES:
!******************************************************************************
!
      ! References to F90 modules
      USE TRACER_MOD, ONLY : STT

      ! Arguments
      INTEGER, INTENT(IN)  :: I, J, L, N
      REAL*8,  INTENT(IN)  :: E_K1, E_K2
      REAL*8,  INTENT(OUT) :: LOST_Hg0

      ! Local variables
      REAL*8               :: OLD_Hg0, NEW_Hg0

      !=================================================================
      ! RXN_Hg0 begins here!
      !=================================================================

      ! Error check tracer number
      IF ( N < 1 ) RETURN

      ! Initial concentration of Hg(0) [kg]
      OLD_Hg0      = MAX( STT(I,J,L,N), SMALLNUM )

      ! New concentration of Hg0 after oxidation rxns [kg]
      NEW_Hg0      = MAX( ( OLD_Hg0 * E_K1 * E_K2 ), SMALLNUM )

      ! Save back into STT array [kg]
      STT(I,J,L,N) = NEW_Hg0

      ! Compute amount of Hg(0) which has been oxidized into Hg(II) [kg]
      LOST_Hg0     = OLD_Hg0 - NEW_Hg0

      ! Return to calling program
      END SUBROUTINE RXN_Hg0

!------------------------------------------------------------------------------

      SUBROUTINE RXN_Hg2_DRYD( I,      J,        L,      
     &                         N,      RKT,      E_RKT, 
     &                         E_SALT, LOST_Hg0, DTCHEM )
!
!******************************************************************************
!  Subroutine RXN_Hg2_DRYD computes the new concentration of Hg(II) from the
!  converted Hg(0) plus the drydep of Hg(II). (eck, bmy, 12/14/04, 2/27/06)
! 
!  Arguments as Input:
!  ============================================================================
!  (1-3) I, J, L  (INTEGER) : GEOS-CHEM lon, lat, alt grid box indices
!  (4  ) N        (INTEGER) : Index for Hg(II) total or tagged tracers
!  (5  ) RKT      (REAL*8 ) : Value of R * k * Temp for drydep [unitless]
!  (6  ) E_RKT    (REAL*8 ) : Value of EXP( - RKT ) [unitless]
!  (7  ) E_SALT   (REAL*8 ) : Value of EXP( -KSALT * DTCHEM ) 
!  (8  ) LOST_Hg0 (REAL*8 ) : Amount of Hg(0) that became Hg(II) [kg]
!  (9  ) DTCHEM   (INTEGER) : Chemistry timestep [s]
!
!  NOTES:
!  (1  ) Now use 2 different solutions, depending on whether or not LOST_Hg0 
!         (amt of converted Hg0 -> Hg2) is positive.  Also now archive the 
!         amount of total Hg(II) tracer lost to drydep for the ocean flux 
!         routines in "ocean_mercury_mod.f". (sas, bmy, 1/19/05)
!  (2  ) Remove references to "diag_mod.f" and "CMN_DIAG".  Now save drydep
!         fluxes into T44 array. (bmy, 2/24/05)
!  (3  ) Now references XNUMOL from "tracer_mod.f" (bmy, 10/25/05)
!  (4  ) Now add E_SALT to account for uptake by sea salt.  Now references
!         IS_WATER from "dao_mod.f".  Now uses ID_Hg2 and ID_Hg_tot from
!         "tracerid_mod.f".  Now references LDYNOCEAN from "logical_mod.f".
!         Now do not call ADD_Hg2_DD if we are not using the dynamic ocean
!         module. (eck, cdh, bmy, 2/27/06)
!******************************************************************************
!
      ! References to F90 modules
      USE DAO_MOD,           ONLY : IS_WATER
      USE GRID_MOD,          ONLY : GET_AREA_CM2
      USE LOGICAL_MOD,       ONLY : LDYNOCEAN
      USE OCEAN_MERCURY_MOD, ONLY : ADD_Hg2_DD
      USE TRACER_MOD,        ONLY : STT, XNUMOL
      USE TRACERID_MOD,      ONLY : IS_Hg2

#     include "CMN_SIZE"          ! Size parameters
#     include "CMN_DIAG"          ! ND44

      ! Arguments
      INTEGER, INTENT(IN)        :: I,   J,     L,      N
      REAL*8,  INTENT(IN)        :: RKT, E_RKT, E_SALT, LOST_Hg0, DTCHEM

      ! Local variables
      INTEGER                    :: NN
      REAL*8                     :: AREA_CM2, DRYDEP, OLD_Hg2, NEW_Hg2

      !=================================================================
      ! RXN_Hg2_DRYD begins here!
      !=================================================================

      ! Error check tracer number
      IF ( N < 1 ) RETURN

      ! Initial concentration of Hg(II) [kg]
      OLD_Hg2 = MAX( STT(I,J,L,N), SMALLNUM )

      IF ( LOST_Hg0 > 0d0 ) THEN

         !--------------------------------------------------------------
         ! CASE 1a: Amt of Hg(0) converted to Hg(II) is positive
         !          Use Solution #1 of the differential equation
         !--------------------------------------------------------------

         ! Concentration of Hg(II) after drydep [kg]
         ! Factor in drydep at this time
         NEW_Hg2 = ( OLD_Hg2      *         E_RKT   )  + 
     &             ( LOST_Hg0/RKT * ( 1d0 - E_RKT ) )
 
      ELSE

         !--------------------------------------------------------------
         ! CASE 1b: Amt of Hg(0) converted to Hg(II) is negative or zero
         !          Use Solution #2 of the differential equation
         !--------------------------------------------------------------

         ! New concentration of Hg(II) [kg]
         ! Assume drydep takes place after chemistry 
         NEW_Hg2 = ( OLD_Hg2 + LOST_Hg0 ) * E_RKT 

      ENDIF

      ! Also account for uptake of Hg(II) by sea salt aerosol
      IF ( IS_WATER( I, J ) ) NEW_Hg2 = NEW_Hg2 * E_SALT

      ! Save back into STT array [kg]
      NEW_Hg2      = MAX( NEW_Hg2, SMALLNUM )
      STT(I,J,L,N) = NEW_Hg2

      !=================================================================
      ! Compute amount of Hg(II) lost to drydep
      !=================================================================

      ! Amount of Hg(II) lost to dry deposition [kg]
      DRYDEP   = ( OLD_Hg2 - NEW_Hg2 ) + LOST_Hg0

      ! Archive Hg(II) lost to drydep [kg] for the ocean mercury flux routines
      ! in "ocean_mercury_mod.f" if necessary.  Do not call ADD_Hg2_DD if the 
      ! dynamic ocean model is turned off. (sas, bmy, 2/27/06)
      !---------------------------------------------------------------------
      ! Prior to 2/27/06:
      !IF ( IS_Hg2(N) ) CALL ADD_Hg2_DD( I, J, N, DRYDEP )
      !---------------------------------------------------------------------
      IF ( IS_Hg2(N) .and. LDYNOCEAN ) THEN
         CALL ADD_Hg2_DD( I, J, N, DRYDEP )
      ENDIF

      !=================================================================
      ! ND44 diagnostic: drydep flux of Hg(II) [molec/cm2/s]
      !=================================================================
      IF ( ND44 > 0 ) THEN
      
         ! Grid box surface area [cm2]
         AREA_CM2     = GET_AREA_CM2( J )
      
         ! Amt of Hg(II) lost to drydep [molec/cm2/s]
         DRYDEP       = DRYDEP * XNUMOL(N) / ( AREA_CM2 * DTCHEM )
      
         ! Archive Hg(II) drydep flux in T44 array [molec/cm2/s]
         T44(I,J,L,N) = T44(I,J,L,N) + DRYDEP

      ENDIF
          
      ! Return to calling program
      END SUBROUTINE RXN_Hg2_DRYD

!------------------------------------------------------------------------------

      SUBROUTINE RXN_Hg2( I, J, L, N, LOST_Hg0 )
!
!******************************************************************************
!  Subroutine RXN_Hg2 computes the new concentration of Hg(II) which is the 
!  old concentration plus the amount of converted Hg(0) (eck, bmy, 12/14/04)
! 
!  Arguments as Input:
!  ============================================================================
!  (1-3) I, J, L  (INTEGER) : GEOS-CHEM lon, lat, alt grid box indices
!  (4  ) N        (INTEGER) : Index for Hg(II) total or tagged tracers
!  (5  ) LOST_Hg0 (REAL*8 ) : Amount of Hg(0) that has become Hg(II) [kg]
!
!  NOTES:
!******************************************************************************
!
      ! References to F90 modules
      USE TRACER_MOD, ONLY : STT

      ! Arguments
      INTEGER, INTENT(IN)  :: I, J, L, N
      REAL*8,  INTENT(IN)  :: LOST_Hg0

      ! Local variables
      REAL*8               :: OLD_Hg2, NEW_Hg2
      
      !=================================================================
      ! RXN_Hg2 begins here!
      !=================================================================

      ! Error check tracer number
      IF ( N < 1 ) RETURN

      ! Initial concentration of Hg(II) [kg]
      OLD_Hg2      = MAX( STT(I,J,L,N), SMALLNUM )

      ! New concentration of Hg(II) is the old concentration
      ! plus the amount of Hg(0) which was converted to Hg(II) [kg]
      NEW_Hg2      = MAX( ( OLD_Hg2 + LOST_Hg0 ), SMALLNUM )
      
      ! Save new concentration of Hg(II) to STT array [kg]
      STT(I,J,L,N) = NEW_Hg2

      ! Return to calling program
      END SUBROUTINE RXN_Hg2

!------------------------------------------------------------------------------

      SUBROUTINE CHEM_HgP( V_DEP_HgP )
!
!******************************************************************************
!  Subroutine CHEM_HgP is the chemistry subroutine for HgP (particulate
!  mercury.  HgP is lost via dry deposition. (eck, bmy, 12/7/04, 1/9/06)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) V_DEP_HgP (REAL*8) : Dry deposition velocity for Hg(II) [cm/s]
!     
!  NOTES:
!  (1 ) Removed references to AD44 and "CMN_DIAG".  Now compute drydep for all
!        levels within the PBL.  Now references ND03, LD03 from "diag03_mod.f".
!        (bmy, 2/24/05)
!  (2 ) Now references XNUMOL & XNUMOLAIR from "tracer_mod.f" (bmy, 10/25/05)
!  (3 ) Now use ID_HgP index array from "tracerid_mod.f". (cdh, bmy, 1/9/06)
!******************************************************************************
!
      ! Refernces to F90 modules
      USE DIAG03_MOD,   ONLY : ND03, LD03
      USE GRID_MOD,     ONLY : GET_AREA_CM2
      USE LOGICAL_MOD,  ONLY : LSPLIT
      USE PBL_MIX_MOD,  ONLY : GET_FRAC_UNDER_PBLTOP, GET_PBL_MAX_L
      USE TIME_MOD,     ONLY : GET_TS_CHEM
      USE TRACER_MOD,   ONLY : STT,       XNUMOL
      USE TRACERID_MOD, ONLY : ID_HgP,    N_Hg_CATS

#     include "CMN_SIZE"  ! Size parameters

      ! Arguments
      REAL*8, INTENT(IN) :: V_DEP_HgP(IIPAR,JJPAR)

      ! Local variables
      INTEGER            :: I,      J,    L,          NN,  PBL_MAX
      REAL*8             :: DTCHEM, E_KT, F_UNDER_TOP

      !=================================================================
      ! CHEM_HgP begins here!
      !=================================================================

      ! Chemistry timestep [s]
      DTCHEM  = GET_TS_CHEM() * 60d0

      ! Maximum extent of the PBL
      PBL_MAX = GET_PBL_MAX_L()

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, F_UNDER_TOP, E_KT, NN )
!$OMP+SCHEDULE( DYNAMIC )
      DO L = 1, PBL_MAX
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! Fraction of box (I,J,L) underneath the PBL top [unitless]
         F_UNDER_TOP = GET_FRAC_UNDER_PBLTOP( I, J, L )

         ! If we are in the PBL and there is nonzero drydep vel ...
         IF ( F_UNDER_TOP > 0d0 .and. V_DEP_HgP(I,J) > 0d0 ) THEN

            ! Pre-compute the exponential term for use below [unitless]
            E_KT = EXP( -V_DEP_HgP(I,J) * DTCHEM * F_UNDER_TOP  )

            ! Compute new conc of HgP tracers after drydep [kg]
            DO NN = 1, N_Hg_CATS
               
               ! Do dry deposition (skip undefined HgP tagged tracers)
               IF ( ID_HgP(NN) > 0 ) THEN
                  CALL RXN_HgP_DRYD( I, J, L, ID_HgP(NN), E_KT, DTCHEM )
               ENDIF

            ENDDO
         ENDIF
      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      ! Return to callwing program
      END SUBROUTINE CHEM_HgP

!------------------------------------------------------------------------------

      SUBROUTINE RXN_HgP_DRYD( I, J, L, N, E_KT, DTCHEM )
!
!******************************************************************************
!  Subroutine RXN_Hg0_DRYD computes the new concentration of HgP 
!  after dry deposition.  ND44 diagnostics are also archived. 
!  (eck, bmy, 12/14/04, 1/9/06)
! 
!  Arguments as Input:
!  ============================================================================
!  (1-3) I, J, L  (INTEGER) : GEOS-CHEM lon, lat, alt grid box indices
!  (4  ) N        (INTEGER) : Index for HgP total or tagged tracers
!  (5  ) E_KT     (REAL*8 ) : Value of EXP( -KT ) for drydep rxn [unitless]
!  (6  ) DTCHEM   (INTEGER) : Chemistry timestep [s]
!
!  NOTES:
!  (1  ) Remove references to "diag_mod.f" and "CMN_DIAG".  Now save drydep
!         fluxes into T44 array. (bmy, 2/24/05)
!  (2  ) Now references XNUMOL from "tracer_mod.f" (bmy, 10/25/05)
!  (3  ) Now uses ID_HgP and ID_Hg_tot from "tracerid_mod.f" (cdh, bmy, 1/9/06)
!******************************************************************************
!
      ! References to F90 modules
      USE GRID_MOD,   ONLY : GET_AREA_CM2
      USE TRACER_MOD, ONLY : STT, XNUMOL

#     include "CMN_SIZE"   ! Size parameters
#     include "CMN_DIAG"   ! ND44 

      ! Arguments
      INTEGER, INTENT(IN) :: I, J, L, N
      REAL*8,  INTENT(IN) :: DTCHEM, E_KT

      ! Local variables
      INTEGER             :: NN
      REAL*8              :: AREA_CM2, DRYDEP, OLD_HgP, NEW_HgP

      !=================================================================
      ! RXN_HgP_DRYD begins here!
      !=================================================================

      ! Error check tracer number
      IF ( N < 1 ) RETURN

      ! Initial concentration of HgP [kg]
      OLD_HgP      = MAX( STT(I,J,L,N), SMALLNUM )

      ! New concentration of HgP after drydep [kg]
      NEW_HgP      = MAX( ( OLD_HgP * E_KT ), SMALLNUM )

      ! Save new concentration of HgP in STT [kg]
      STT(I,J,L,N) = NEW_HgP

      !=================================================================
      ! ND44 diagnostic: drydep flux of Hg(II) [molec/cm2/s]
      !=================================================================
      IF ( ND44 > 0 ) THEN
         
         ! Grid box surface area [cm2]
         AREA_CM2     = GET_AREA_CM2( J )

         ! Amt of Hg(II) lost to drydep [molec/cm2/s]
         DRYDEP       = OLD_HgP - NEW_HgP  
         DRYDEP       = ( DRYDEP * XNUMOL(N) ) / ( AREA_CM2 * DTCHEM )

         ! Archive Hg(II) drydep flux in T44 array [molec/cm2/s]
         T44(I,J,L,N) = T44(I,J,L,N) + DRYDEP

      ENDIF

      ! Return to calling program
      END SUBROUTINE RXN_HgP_DRYD

!------------------------------------------------------------------------------

      SUBROUTINE EMISSMERCURY
!
!******************************************************************************
!  Subroutine EMISSMERCURY is the driver routine for mercury emissions.
!  (eck, bmy, 12/7/04, 2/24/06)
! 
!  NOTES:
!  (1 ) Now call OCEAN_MERCURY_FLUX from "ocean_mercury_mod.f" to compute 
!        the emissions of Hg0 from the ocean instead of reading it from disk.
!        (sas, bmy, 1/20/05)
!  (2 ) Now no longer call COMPUTE_FEMIS, since we can get the same information
!        from routine GET_FRAC_OF_PBL in "pbl_mix_mod.f" (bmy, 2/22/05)
!  (3 ) Now modified for new ocean mercury module. (cdh, sas, bmy, 2/24/06)
!******************************************************************************
!
      ! References to F90 modules
      USE ERROR_MOD,         ONLY : DEBUG_MSG
      USE LOGICAL_MOD,       ONLY : LPRT, LDYNOCEAN
      USE OCEAN_MERCURY_MOD, ONLY : OCEAN_MERCURY_FLUX
      USE TIME_MOD,          ONLY : GET_MONTH
      USE TRACER_MOD,        ONLY : STT
      
#     include "CMN_SIZE"          ! Size parameters

      ! Local variables
      LOGICAL, SAVE :: FIRST = .TRUE. 
      INTEGER       :: THISMONTH
     
      !=================================================================
      ! EMISSMERCURY begins here!
      !=================================================================

      ! First-time initialization
      IF ( FIRST ) THEN

         !-------------------------------------------------
         ! Prior to 2/24/06:
         ! Now call this from "input_mod.f" (bmy, 2/24/06)
         !! Allocate arrays (if not done before)
         !CALL INIT_MERCURY
         !-------------------------------------------------

         ! Read anthro, ocean, land emissions of Hg from disk
         CALL MERCURY_READYR

         ! Reset first-time flag
         FIRST = .FALSE.
      ENDIF

      !=================================================================
      ! Call emission routines for Hg(0), Hg(II), and Hg(P)
      !=================================================================

      ! Ocean flux of Hg(0)
      IF ( LDYNOCEAN ) THEN
         CALL OCEAN_MERCURY_FLUX( EHg0_oc )
         IF ( LPRT ) CALL DEBUG_MSG( '### EMISSMERCURY: a OCEAN_FLUX' )
      ENDIF

      ! Add Hg(0) source into STT [kg]
      CALL SRCHg0
      IF ( LPRT ) CALL DEBUG_MSG( '### EMISSMERCURY: a SRCHg0' )

      ! Add Hg(II) source into STT [kg]
      CALL SRCHg2
      IF ( LPRT ) CALL DEBUG_MSG( '### EMISSMERCURY: a SRCHg2' )

      ! Add HgP source into STT [kg]
      CALL SRCHgP
      IF ( LPRT ) CALL DEBUG_MSG( '### EMISSMERCURY: a SRCHgP' )

      ! Return to calling program
      END SUBROUTINE EMISSMERCURY

!------------------------------------------------------------------------------

      SUBROUTINE SRCHg0
!
!******************************************************************************
!  Subroutine SRCHg0 is the subroutine for Hg(0) emissions.  
!  Emissions of Hg(0) will be distributed throughout the boundary layer. 
!  (eck, cdh, bmy, 1/21/05, 1/9/06)
! 
!  Arguments as Input/Output:
!  ============================================================================
!  (1 ) TC (REAL*8) : Tracer concentration of Hg(0) [kg]
!
!  NOTES:
!  (1 ) Now use diagnostic arrays from "diag03_mod.f" (bmy, 1/21/05)
!  (2 ) Now references GET_FRAC_OF_PBL and GET_PBL_MAX_L from "pbl_mix_mod.f".
!        Remove reference to FEMIS. (bmy, 2/22/05)
!  (3 ) EHg0_an is now a 2-D array.  Modified for new ocean mercury module.
!        Now use ID_Hg0 index array from "tracerid_mod.f" (cdh, bmy, 1/9/06).
!******************************************************************************
!
      ! Reference to diagnostic arrays
      USE DIAG03_MOD,   ONLY : AD03, ND03
      USE ERROR_MOD,    ONLY : ERROR_STOP
      USE LOGICAL_MOD,  ONLY : LSPLIT
      USE PBL_MIX_MOD,  ONLY : GET_FRAC_OF_PBL, GET_PBL_MAX_L
      USE TIME_MOD,     ONLY : GET_TS_EMIS
      USE TRACER_MOD,   ONLY : STT
      USE TRACERID_MOD, ONLY : ID_Hg_tot, ID_Hg_na, ID_Hg_eu, ID_Hg_as
      USE TRACERID_MOD, ONLY : ID_Hg_rw,  ID_Hg_oc, ID_Hg_ln, ID_Hg_nt
      USE TRACERID_MOD, ONLY : ID_Hg0
      
#     include "CMN_SIZE"     ! Size parameters

      ! Local variables
      INTEGER               :: I,      J,    L,        N,    PBL_MAX
      REAL*8                :: DTSRCE, E_Hg, F_OF_PBL, T_Hg, T_Hg_An

      !=================================================================
      ! SRCHg0 begins here!
      !=================================================================

      ! Emission timestep [s]
      DTSRCE  = GET_TS_EMIS() * 60d0

      ! Maximum extent of the PBL [model levels]
      PBL_MAX = GET_PBL_MAX_L() 

      ! Loop over grid boxes
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, N, T_Hg_An, T_Hg, F_OF_PBL, E_Hg )
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! Compute total anthropogenic Hg(0) emissions
         T_Hg_An = EHg0_an(I,J)

         ! Compute total Hg(0) emissions (anthro+oceans+land+natural)
         T_Hg    = T_Hg_An        + EHg0_oc(I,J,ID_Hg_tot)   + 
     &             EHg0_ln(I,J)   + EHg0_nt(I,J)

         !==============================================================
         ! Partition Hg0 throughout PBL; store into STT [kg]
         !==============================================================

         ! Loop up to max PBL level
         DO L = 1, PBL_MAX

            ! Fraction of box (I,J,L) w/in the PBL [unitless]
            F_OF_PBL        = GET_FRAC_OF_PBL( I, J, L )

            !-----------------
            ! Total Hg tracer
            !-----------------
            N               = ID_Hg0(ID_Hg_tot)
            E_Hg            = F_OF_PBL     * T_Hg
            STT(I,J,L,N)    = STT(I,J,L,N) + ( E_Hg * DTSRCE )

            !-----------------
            ! Tagged tracers
            !-----------------
            IF ( LSPLIT ) THEN 

               !--------------------
               ! Primary emissions
               !--------------------

               ! Anthro Hg0 by region
               N            = AN_Hg0(I,J)
               E_Hg         = F_OF_PBL     * T_Hg_an
               STT(I,J,L,N) = STT(I,J,L,N) + ( E_Hg * DTSRCE )

               ! Land re-emissions of Hg0
               N            = ID_Hg0(ID_Hg_ln)
               E_Hg         = F_OF_PBL     * EHg0_ln(I,J)
               STT(I,J,L,N) = STT(I,J,L,N) + ( E_Hg * DTSRCE ) 

               ! Natural land sources of Hg0
               N            = ID_Hg0(ID_Hg_nt)
               E_Hg         = F_OF_PBL     * EHg0_nt(I,J)
               STT(I,J,L,N) = STT(I,J,L,N) + ( E_Hg * DTSRCE )

               !--------------------
               ! Ocean re-emissions
               !--------------------

               ! Anthro re-emission from ocean in N. AMERICA
               N            = ID_Hg0(ID_Hg_na)
               E_Hg         = F_OF_PBL     * EHg0_oc(I,J,ID_Hg_na)
               STT(I,J,L,N) = STT(I,J,L,N) + ( E_Hg * DTSRCE )

               ! Anthro re-emission from ocean in EUROPE
               N            = ID_Hg0(ID_Hg_eu)
               E_Hg         = F_OF_PBL     * EHg0_oc(I,J,ID_Hg_eu)
               STT(I,J,L,N) = STT(I,J,L,N) + ( E_Hg * DTSRCE )

               ! Anthro re-emission from ocean in ASIA
               N            = ID_Hg0(ID_Hg_as)
               E_Hg         = F_OF_PBL     * EHg0_oc(I,J,ID_Hg_as)
               STT(I,J,L,N) = STT(I,J,L,N) + ( E_Hg * DTSRCE )

               ! Anthro re-emission from ocean in REST OF WORLD
               N            = ID_Hg0(ID_Hg_rw)
               E_Hg         = F_OF_PBL     * EHg0_oc(I,J,ID_Hg_rw)
               STT(I,J,L,N) = STT(I,J,L,N) + ( E_Hg * DTSRCE )

               ! Re-emission from ocean in OCEAN 
               N            = ID_Hg0(ID_Hg_oc)
               E_Hg         = F_OF_PBL     * EHg0_oc(I,J,ID_Hg_oc)
               STT(I,J,L,N) = STT(I,J,L,N) + ( E_Hg * DTSRCE )

               ! Re-emission from ocean in LAND REEMISSION
               N            = ID_Hg0(ID_Hg_ln)
               E_Hg         = F_OF_PBL     * EHg0_oc(I,J,ID_Hg_ln)
               STT(I,J,L,N) = STT(I,J,L,N) + ( E_Hg * DTSRCE )

               ! Re-emission from ocean in NATURAL
               N            = ID_Hg0(ID_Hg_nt)
               E_Hg         = F_OF_PBL     * EHg0_oc(I,J,ID_Hg_nt)
               STT(I,J,L,N) = STT(I,J,L,N) + ( E_Hg * DTSRCE )
               
            ENDIF
         ENDDO
        
         !==============================================================
         ! ND03 diagnostic: Total Hg(0) emissions [kg]
         ! 1=anthro; 3=from ocean; 4=land re-emission; 5=natural src
         !==============================================================
         IF ( ND03 > 0 ) THEN
            N = ID_Hg_tot
            AD03(I,J,1) = AD03(I,J,1) + ( T_Hg_An        * DTSRCE )
            AD03(I,J,3) = AD03(I,J,3) + ( EHg0_oc(I,J,N) * DTSRCE )
            AD03(I,J,4) = AD03(I,J,4) + ( EHg0_ln(I,J)   * DTSRCE )
            AD03(I,J,5) = AD03(I,J,5) + ( EHg0_nt(I,J)   * DTSRCE )
         ENDIF
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      ! Return to calling program
      END SUBROUTINE SRCHg0

!-----------------------------------------------------------------------------

      SUBROUTINE SRCHg2
!
!******************************************************************************
!  Subroutine SRCHg2 is the subroutine for Hg(II) emissions.  
!  Emissions of Hg(II) will be distributed throughout the boundary layer. 
!  (eck, bmy, 12/7/04, 1/9/06)
! 
!  Arguments as Input/Output:
!  ============================================================================
!  (1 ) TC (REAL*8) : Tracer concentration of Hg(II) [kg]
!
!  NOTES:
!  (1 ) Now use diagnostic arrays from "diag03_mod.f" (bmy, 1/21/05)
!  (2 ) Now references GET_FRAC_OF_PBL and GET_PBL_MAX_L from "pbl_mix_mod.f".
!        Remove reference to FEMIS. (bmy, 2/22/05)
!  (3 ) EHg2_an is now a 2-D array.  Now use ID_Hg2 index array from 
!        "tracerid_mod.f". (eck, cdh, bmy, 1/9/06)
!******************************************************************************
!
      ! Reference to F90 modules
      USE DIAG03_MOD,   ONLY : AD03, ND03
      USE ERROR_MOD,    ONLY : ERROR_STOP
      USE LOGICAL_MOD,  ONLY : LSPLIT
      USE PBL_MIX_MOD,  ONLY : GET_FRAC_OF_PBL, GET_PBL_MAX_L
      USE TIME_MOD,     ONLY : GET_TS_EMIS
      USE TRACER_MOD,   ONLY : STT
      USE TRACERID_MOD, ONLY : ID_Hg_tot, ID_Hg2

#     include "CMN_SIZE"     ! Size parameters

      ! Local variables
      INTEGER               :: I,      J,        L,    N,   PBL_MAX
      REAL*8                :: DTSRCE, F_OF_PBL, E_Hg 

      !=================================================================
      ! SRCHg2 begins here!
      !=================================================================
      
      ! Emission timestep [s]
      DTSRCE  = GET_TS_EMIS() * 60d0

      ! Maximum extent of the PBL [model levels]
      PBL_MAX = GET_PBL_MAX_L()

      ! Loop over grid boxes
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, F_OF_PBL, E_Hg, N )
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! Loop up to the max PBL layer
         DO L = 1, PBL_MAX

            ! Fraction of box (I,J,L) w/in the PBL [unitless]
            F_OF_PBL           = GET_FRAC_OF_PBL( I, J, L )

            ! Partition total Hg2 into box (I,J,L) [kg]
            E_Hg               = F_OF_PBL * EHg2_an(I,J) * DTSRCE

            !---------------------------
            ! Total anthro Hg(II) [kg]
            !---------------------------
            N               = ID_Hg2(ID_Hg_tot)
            STT(I,J,L,N)    = STT(I,J,L,N) + E_Hg

            !---------------------------
            ! Tagged anthro Hg(II) [kg]
            !---------------------------
            IF ( LSPLIT ) THEN 
               N            = AN_Hg2(I,J)
               STT(I,J,L,N) = STT(I,J,L,N) + E_Hg
            ENDIF
         ENDDO
            
         !-------------------------------
         ! ND03 diag: Anthro Hg(II) [kg]
         !-------------------------------
         IF ( ND03 > 0 ) THEN
            AD03(I,J,6) = AD03(I,J,6) + ( EHg2_an(I,J) * DTSRCE )
         ENDIF
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      ! Return to calling program
      END SUBROUTINE SRCHg2

!-----------------------------------------------------------------------------

      SUBROUTINE SRCHgP
!
!******************************************************************************
!  Subroutine SRCHgP is the subroutine for HgP emissions.  
!  Emissions of HgP will be distributed throughout the boundary layer. 
!  (eck, bmy, 12/7/04, 1/9/06)
! 
!  Arguments as Input/Output:
!  ============================================================================
!  (1 ) TC (REAL*8) : Tracer concentration of Hg(II) [kg]
!
!  NOTES:
!  (1 ) Now use diagnostic arrays from "diag03_mod.f" (bmy, 1/21/05)
!  (2 ) Now references GET_FRAC_OF_PBL and GET_PBL_MAX_L from "pbl_mix_mod.f".
!        Remove reference to FEMIS. (bmy, 2/22/05)
!  (3 ) EHgP_an is now a 2-D array.  Now use ID_HgP index array from 
!        "tracerid_mod.f". (eck, cdh, bmy, 1/9/06) 
!******************************************************************************
!
      ! Reference to diagnostic arrays
      USE DIAG03_MOD,   ONLY : AD03, ND03
      USE ERROR_MOD,    ONLY : ERROR_STOP
      USE LOGICAL_MOD,  ONLY : LSPLIT
      USE PBL_MIX_MOD,  ONLY : GET_FRAC_OF_PBL, GET_PBL_MAX_L
      USE TIME_MOD,     ONLY : GET_TS_EMIS 
      USE TRACER_MOD,   ONLY : STT
      USE TRACERID_MOD, ONLY : ID_Hg_tot, ID_HgP

#     include "CMN_SIZE"     ! Size paramters

      ! Local variables
      INTEGER               :: I,      J,        L,    N,   PBL_MAX
      REAL*8                :: DTSRCE, F_OF_PBL, E_Hg 

      !=================================================================
      ! SRCHgP begins here!
      !=================================================================

      ! Chemistry timestep [s]
      DTSRCE  = GET_TS_EMIS() * 60d0
      
      ! Maximum extent of the PBL [model levels]
      PBL_MAX = GET_PBL_MAX_L()

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, F_OF_PBL, E_Hg, N )
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! Loop up to PBL top layer
         DO L = 1, PBL_MAX

            ! Fraction of box (I,J,L) w/in the PBL [unitless]
            F_OF_PBL           = GET_FRAC_OF_PBL( I, J, L )

            ! Partition HgP into box (I,J,L) [kg]
            E_Hg               = F_OF_PBL * EHgP_an(I,J) * DTSRCE

            !------------------------
            ! Total anthro HgP [kg]
            !------------------------
            N               = ID_HgP(ID_Hg_tot)
            STT(I,J,L,N)    = STT(I,J,L,N) + E_Hg

            !------------------------
            ! Tagged anthro HgP [kg]
            !------------------------
            IF ( LSPLIT ) THEN
               N            = AN_HgP(I,J)
               STT(I,J,L,N) = STT(I,J,L,N) + E_Hg
            ENDIF
         ENDDO

         !----------------------------
         ! ND03 diag: Anthro HgP [kg]
         !----------------------------
         IF ( ND03 > 0 ) THEN
            AD03(I,J,9)  = AD03(I,J,9) + ( EHgP_an(I,J) * DTSRCE )
         ENDIF
         
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      ! Return to calling program
      END SUBROUTINE SRCHgP

!------------------------------------------------------------------------------

      SUBROUTINE MERCURY_READYR
!
!******************************************************************************
!  Subroutine MERCURY_READYR reads the year-invariant emissions for Mercury
!  from anthropogenic, ocean, and land sources. (eck, bmy, 12/6/04, 3/16/06)
!  
!  NOTES:
!  (1 ) Now read data from mercury_200501 subdirectory.  Now compute oceanic 
!        Hg(0) emissions w/ ocean flux module instead of reading them from 
!        disk.  Now use 1985 TAU values. (sas, bmy, 1/20/05) 
!  (2 ) Now can read data for both GEOS and GCAP grids (bmy, 8/16/05)
!  (3 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (4 ) Now read anthro emissions on GEOS 1x1 grid in DATA_1x1_DIR.  Also
!        keep 2x25 and 4x5 files together in DATA_1x1_DIR.  Also now use new
!        land re-emissions files from Noelle Selin. (eck, bmy, 3/16/06)
!******************************************************************************
!
      ! References to F90 modules
      USE BPCH2_MOD,      ONLY : GET_NAME_EXT_2D, GET_RES_EXT
      USE BPCH2_MOD,      ONLY : GET_TAU0,        READ_BPCH2
      USE DIRECTORY_MOD,  ONLY : DATA_DIR_1x1
      USE LOGICAL_MOD,    ONLY : LDYNOCEAN
      USE REGRID_1x1_MOD, ONLY : DO_REGRID_1x1
      USE TIME_MOD,       ONLY : EXPAND_DATE
      USE TRANSFER_MOD,   ONLY : TRANSFER_2D
      
#     include "CMN_SIZE"       ! Size parameters

      ! Local variables
      INTEGER                 :: NYMD
      REAL*4                  :: ARRAY(IGLOB,JGLOB,1)
      REAL*4                  :: ARRAY1(I1x1,J1x1,1)
      REAL*8                  :: XTAU
      REAL*8, PARAMETER       :: SEC_PER_YR = 365.25d0 * 86400d0
      CHARACTER(LEN=255)      :: FILENAME

      !=================================================================
      ! MERCURY_READYR begins here!
      !
      ! Read annual anthropogenic mercury emissions [kg/s]
      !=================================================================

      ! Mercury data is either for 1995 or 2000
      NYMD     = ( ANTHRO_Hg_YEAR * 10000 ) + 0101
      XTAU     = GET_TAU0( 1, 1, ANTHRO_Hg_YEAR )

      !---------------------------
      ! Hg(0) emissions [kg/s]
      !---------------------------

      ! Filename for anthropogenic mercury source
      FILENAME = TRIM( DATA_DIR_1x1 ) // 
     &           'mercury_200511/GEIA_Hg0.geos.1x1.YYYY'

      ! Add year to the filename
      CALL EXPAND_DATE( FILENAME, NYMD, 000000 )

      ! Echo info
      WRITE( 6, 100 ) TRIM( FILENAME )
 100  FORMAT( '     - MERCURY_READYR: Reading ', a )

      ! Read data in [kg/yr]
      CALL READ_BPCH2( FILENAME, 'HG-SRCE', 1, 
     &                 XTAU,      I1x1,     J1x1,    
     &                 1,         ARRAY1,   QUIET=.TRUE. )

      ! Regrid from 1x1 to the current grid
      CALL DO_REGRID_1x1( 'kg', ARRAY1, EHg0_an )

      ! Convert from [kg/yr] to [kg/s]
      EHg0_an = EHg0_an / SEC_PER_YR

      !---------------------------
      ! Hg(II) emissions [kg/s]
      !---------------------------

      ! Filename for anthropogenic mercury source
      FILENAME = TRIM( DATA_DIR_1x1 ) // 
     &           'mercury_200511/GEIA_Hg2.geos.1x1.YYYY'

      ! Add year to the filename
      CALL EXPAND_DATE( FILENAME, NYMD, 000000 )

      ! Echo info
      WRITE( 6, 100 ) TRIM( FILENAME )

      ! Read data in [kg/yr]
      CALL READ_BPCH2( FILENAME, 'HG-SRCE', 6, 
     &                 XTAU,      I1x1,     J1x1,    
     &                 1,         ARRAY1,   QUIET=.TRUE. )

      ! Regrid from 1x1 to the current grid
      CALL DO_REGRID_1x1( 'kg', ARRAY1, EHg2_an )

      ! Convert from [kg/yr] to [kg/s]
      EHg2_an = EHg2_an / SEC_PER_YR

      !---------------------------
      ! HgP emissions [kg/s]
      !---------------------------

      ! Filename for anthropogenic mercury source
      FILENAME = TRIM( DATA_DIR_1x1 ) // 
     &           'mercury_200511/GEIA_HgP.geos.1x1.YYYY'

      ! Add year to the filename
      CALL EXPAND_DATE( FILENAME, NYMD, 000000 )

      ! Echo info
      WRITE( 6, 100 ) TRIM( FILENAME )

      ! Read data in [kg/yr]
      CALL READ_BPCH2( FILENAME, 'HG-SRCE', 9, 
     &                 XTAU,      I1x1,     J1x1,    
     &                 1,         ARRAY1,   QUIET=.TRUE. )

      ! Regrid from 1x1 to the current grid
      CALL DO_REGRID_1x1( 'kg', ARRAY1, EHgP_an )

      ! Convert from [kg/yr] to [kg/s]
      EHgP_an = EHgP_an / SEC_PER_YR

      !=================================================================
      ! Read annual emissions of anthropogenic Hg(0) which is
      ! re-emitted from the land [kg/s]
      !=================================================================

      ! Use "generic" year 1985 for TAU values
      XTAU     = GET_TAU0( 1, 1, 1985 )

      ! Filename for re-emitted anthropogenic mercury
      FILENAME = TRIM( DATA_DIR_1x1 ) // 
!-----------------------------------------------------------------------------
! Prior to 3/14/06:
! Now use better land-re-emissions files from Noelle (bmy, 3/16/06)
!     &           'mercury_200511/Hg_land.' // GET_NAME_EXT_2D() //
!     &           '.'                       // GET_RES_EXT()   
!-----------------------------------------------------------------------------
     &           'mercury_200511/Hg_land_reemission.' // 
     &            GET_NAME_EXT_2D()   // '.' // GET_RES_EXT()   

      ! Echo info
      WRITE( 6, 100 ) TRIM( FILENAME )

      ! Read data in [kg/yr]
      CALL READ_BPCH2( FILENAME, 'HG-SRCE',     4,  
     &                 XTAU,      IGLOB,        JGLOB,    
     &                 1,         ARRAY(:,:,1), QUIET=.TRUE. )

      ! Cast to REAL*8 and resize
      CALL TRANSFER_2D( ARRAY(:,:,1), EHg0_ln )

      ! Quick fix: change total amt from 2000 to 1500
      EHg0_ln = EHg0_ln * 0.75d0
      
      ! Convert from [kg/yr] to [kg/s]
      EHg0_ln = EHg0_ln / SEC_PER_YR

      !=================================================================
      ! Read annual emissions of Hg(0) from natural land sources [kg/s]
      !=================================================================

      ! Use "generic" year 1985 for TAU values
      XTAU     = GET_TAU0( 1, 1, 1985 )

      ! Filename for natural land-source mercury
      FILENAME = TRIM( DATA_DIR_1x1 )         // 
     &           'mercury_200511/Hg_natural.' // GET_NAME_EXT_2D() //
     &           '.'                          // GET_RES_EXT() 

      ! Echo info
      WRITE( 6, 100 ) TRIM( FILENAME )

      ! TAU value corresponding to the data 
      XTAU   = GET_TAU0( 1, 1, 1995 )

      ! Read data in [kg/yr]
      CALL READ_BPCH2( FILENAME, 'HG-SRCE',     5, 
     &                 XTAU,      IGLOB,        JGLOB,    
     &                 1,         ARRAY(:,:,1), QUIET=.TRUE. )

      ! Cast to REAL*8 and resize
      CALL TRANSFER_2D( ARRAY(:,:,1), EHg0_nt )

      ! Convert from [kg/yr] to [kg/s]
      EHg0_nt = EHg0_nt / SEC_PER_YR

      !=================================================================
      ! Read offline ocean Hg0 emissions (if LDYNOCEAN = .FALSE.)
      !=================================================================
      IF ( .not. LDYNOCEAN ) THEN

         ! File name
         FILENAME = TRIM( DATA_DIR_1x1 )            // 
     &              'mercury_200511/Hg_ocean.geos.' // GET_RES_EXT()    

         ! Echo info
         WRITE( 6, 100 ) TRIM( FILENAME )

         ! TAU value corresponding to the data 
         XTAU  = GET_TAU0( 1, 1, 1985 )

         ! Read data in [kg/yr]
         CALL READ_BPCH2( FILENAME, 'HG-SRCE',     3,  
     &                    XTAU,      IGLOB,        JGLOB,    
     &                    1,         ARRAY(:,:,1), QUIET=.TRUE. )

         ! Cast to REAL*8 and resize
         CALL TRANSFER_2D( ARRAY(:,:,1), EHg0_oc(:,:,1) )

         ! Convert from [kg/yr] to [kg/s]
         EHg0_oc = EHg0_oc / SEC_PER_YR

      ENDIF

      !=================================================================
      ! Print totals to the screen in [Gg/yr]
      !=================================================================
      WRITE( 6, '(a)' ) REPEAT( '=', 79 )
      WRITE( 6, 110   )
      WRITE( 6, '(a)' )
      WRITE( 6, 111   ) SUM( EHg0_an ) * SEC_PER_YR * 1d-6
      WRITE( 6, 113   ) SUM( EHg0_ln ) * SEC_PER_YR * 1d-6
      WRITE( 6, 114   ) SUM( EHg0_nt ) * SEC_PER_YR * 1d-6

      ! Only write ocean total if we are doing offline ocean
      IF ( .not. LDYNOCEAN ) THEN
         WRITE( 6, 117   ) SUM( EHg0_oc ) * SEC_PER_YR * 1d-6 
      ENDIF

      WRITE( 6, 115   ) SUM( EHg2_an ) * SEC_PER_YR * 1d-6
      WRITE( 6, 116   ) SUM( EHgP_an ) * SEC_PER_YR * 1d-6
      WRITE( 6, '(a)' ) REPEAT( '=', 79 )

      ! FORMAT strings
 110  FORMAT( 'M E R C U R Y   E M I S S I O N S' )
 111  FORMAT( 'Total Anthro     Hg(0)  : ', f7.3, ' [Gg/yr]' )
 113  FORMAT( 'Total Re-Emitted Hg(0)  : ', f7.3, ' [Gg/yr]' )
 114  FORMAT( 'Total Natural    Hg(0)  : ', f7.3, ' [Gg/yr]' )
 115  FORMAT( 'Total Anthro     Hg(II) : ', f7.3, ' [Gg/yr]' )
 116  FORMAT( 'Total Anthro     HgP    : ', f7.3, ' [Gg/yr]' )
 117  FORMAT( 'Total Ocean      Hg(0)  : ', f7.3, ' [Gg/yr]' )

      ! Return to calling program
      END SUBROUTINE MERCURY_READYR

!------------------------------------------------------------------------------

      FUNCTION GET_LWC( T ) RESULT( LWC )
!
!******************************************************************************
!  Function GET_LWC returns the cloud liquid water content at a GEOS-CHEM
!  grid box as a function of temperature. (rjp, bmy, 10/31/02, 12/7/04)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) T (REAL*8) : Temperature value at a GEOS-CHEM grid box [K]
!
!  NOTES:
!******************************************************************************
!
      ! Arguments
      REAL*8, INTENT(IN) :: T

      ! Function value
      REAL*8             :: LWC

      !=================================================================
      ! GET_LWC begins here!
      !=================================================================

      ! Compute Liquid water content in [g/m3]
      IF ( T >= 280.d0 .AND. T <= 293.d0 ) THEN
         LWC = 0.32d0 - 0.0060d0 * ( T - 273.D0 ) 
 
      ELSE IF ( T >= 248.d0 .AND. T < 280.d0 ) THEN
         LWC = 0.23d0 + 0.0065d0 * ( T - 273.D0 )

      ELSE IF ( T < 248.d0 ) THEN
         LWC = 0.07d0

      ENDIF

      ! Convert from [g/m3] to [m3/m3]
      LWC = LWC * 1.D-6         

      ! Return to calling program
      END FUNCTION GET_LWC

!------------------------------------------------------------------------------

      FUNCTION GET_VCLDF( I, J, L ) RESULT( VCLDF )
! 
!******************************************************************************
!  Subroutine GET_VCLDF computes the volume cloud fraction for Hg0 and Hg2
!  chemistry. (eck, bmy, 12/6/04)
!
!  Arguments as Input:
!  ============================================================================
!  (1-3) I, J, L (INTEGER) : GEOS-CHEM lon, lat, alt indices 
!
!  References:
!  ============================================================================
!  (1  ) Sundqvist et al. [1989]
!
!  NOTES:
!  (1  ) Copied from "sulfate_mod.f" but was made into a function since we are
!         already looping over (I,J,L) in CHEM_Hg0_Hg2 (eck, bmy, 12/6/04)
!******************************************************************************
!
      ! References to F90 modules 
      USE DAO_MOD,      ONLY : RH
      USE PRESSURE_MOD, ONLY : GET_PCENTER, GET_PEDGE
      
      ! Arguments
      INTEGER, INTENT(IN)  :: I, J, L

      ! Local variables
      REAL*8               :: PRES, PSFC, RH2, R0, B0
      REAL*8,  PARAMETER   :: ZRT = 0.60d0, ZRS = 0.99d0
	
      ! Function value
      REAL*8               :: VCLDF
	
      !=================================================================
      ! GET_VCLDF begins here!
      !=================================================================
     
      ! Surface pressure
      PSFC = GET_PEDGE(I,J,1)         

      ! Pressure at center of grid box L
      PRES = GET_PCENTER(I,J,L)
        
      ! RH (from "dao_mod.f") is relative humidity [%]
      ! Convert to fraction and store in RH2
      RH2  = RH(I,J,L) * 1.0d-2
        
      ! Terms from Sundqvist ???
      R0   = ZRT + ( ZRS - ZRT ) * EXP( 1d0 - ( PSFC / PRES )**2.5 )
      B0   = ( RH2 - R0 ) / ( 1d0 - R0 )
	
      ! Force B0 into the range 0-1
      IF ( RH2 < R0  ) B0 = 0d0
      IF ( B0  > 1d0 ) B0 = 1d0

      ! Volume cloud fraction
      VCLDF = 1d0 - SQRT( 1d0 - B0 )

      ! Return to calling program
      END FUNCTION GET_VCLDF

!------------------------------------------------------------------------------

      FUNCTION GET_O3( I, J, L ) RESULT( O3_MOLEC_CM3 )
!
!******************************************************************************
!  Function GET_O3 returns monthly mean O3 for offline sulfate aerosol
!  simulations. (bmy, 12/16/02)
!
!  Arguments as Input:
!  ============================================================================
!  (1-3) I, J, L   (INTEGER) : Grid box indices for lon, lat, vertical level
!
!  NOTES:
!  (1 ) We assume SETTRACE has been called to define IDO3. (bmy, 12/16/02)
!  (2 ) Now reference inquiry functions from "tracer_mod.f" (bmy, 7/20/04)
!******************************************************************************
!
      ! References to F90 modules
      USE DAO_MOD,       ONLY : AD
      USE GLOBAL_O3_MOD, ONLY : O3

#     include "CMN_SIZE"  ! Size parameters

      ! Arguments
      INTEGER, INTENT(IN) :: I, J, L

      ! Local variables
      REAL*8              :: O3_MOLEC_CM3

      ! External functions
      REAL*8, EXTERNAL    :: BOXVL
      
      !=================================================================
      ! GET_O3 begins here!
      !=================================================================

      ! Get tropospheric O3 [v/v] for this gridbox & month
      ! and convert to [molec/cm3] (eck, 12/2/04)
      IF ( L <= LLTROP ) THEN
         O3_MOLEC_CM3 = O3(I,J,L) * ( 6.022d23 / 28.97d-3 ) * 
     &                  AD(I,J,L)  /  BOXVL(I,J,L)
      ELSE
         O3_MOLEC_CM3 = 0d0
      ENDIF

      ! Return to calling program
      END FUNCTION GET_O3

!------------------------------------------------------------------------------

      FUNCTION GET_OH( I, J, L ) RESULT( OH_MOLEC_CM3 )
!
!******************************************************************************
!  Function GET_OH returns monthly mean OH and imposes a diurnal variation.
!  (eck, bmy, 12/7/04)
!
!  Arguments as Input:
!  ============================================================================
!  (1-3) I, J, L (INTEGER) : Grid box indices for lon, lat, vertical level
!
!  NOTES:
!******************************************************************************
!
      ! References to F90 modules
      USE DAO_MOD,       ONLY : SUNCOS
      USE GLOBAL_OH_MOD, ONLY : OH
      USE TIME_MOD,      ONLY : GET_TS_CHEM

#     include "CMN_SIZE"  ! Size parameters

      ! Arguments
      INTEGER, INTENT(IN) :: I, J, L

      ! Local variables
      INTEGER             :: JLOOP
      REAL*8              :: OH_MOLEC_CM3
 
      !=================================================================
      ! GET_OH begins here!
      !=================================================================

      ! 1-D grid box index for SUNCOS
      JLOOP = ( (J-1) * IIPAR ) + I

      ! Test for sunlight...
      IF ( SUNCOS(JLOOP) > 0d0 .and. TCOSZ(I,J) > 0d0 ) THEN

         ! Impose a diurnal variation on OH during the day
         OH_MOLEC_CM3 = OH(I,J,L)                      *           
     &                  ( SUNCOS(JLOOP) / TCOSZ(I,J) ) *
     &                  ( 1440d0        / GET_TS_CHEM() )

         ! Make sure OH is not negative
         OH_MOLEC_CM3 = MAX( OH_MOLEC_CM3, 0d0 )
               
      ELSE

         ! At night, OH goes to zero
         OH_MOLEC_CM3 = 0d0

      ENDIF

      ! Return to calling program
      END FUNCTION GET_OH

!------------------------------------------------------------------------------

      SUBROUTINE OHNO3TIME
!
!******************************************************************************
!  Subroutine OHNO3TIME computes the sum of cosine of the solar zenith
!  angle over a 24 hour day, as well as the total length of daylight. 
!  This is needed to scale the offline OH and NO3 concentrations.
!  (rjp, bmy, 12/16/02, 12/8/04)
!
!  NOTES:
!  (1 ) Copy code from COSSZA directly for now, so that we don't get NaN
!        values.  Figure this out later (rjp, bmy, 1/10/03)
!  (2 ) Now replace XMID(I) with routine GET_XMID from "grid_mod.f".  
!        Now replace RLAT(J) with routine GET_YMID_R from "grid_mod.f". 
!        Removed NTIME, NHMSb from the arg list.  Now use GET_NHMSb,
!        GET_ELAPSED_SEC, GET_TS_CHEM, GET_DAY_OF_YEAR, GET_GMT from 
!        "time_mod.f". (bmy, 3/27/03)
!  (3 ) Now store the peak SUNCOS value for each surface grid box (I,J) in 
!        the COSZM array. (rjp, bmy, 3/30/04)
!  (4 ) Also added parallel loop over grid boxes (eck, bmy, 12/8/04)
!******************************************************************************
!
      ! References to F90 modules
      USE GRID_MOD, ONLY : GET_XMID,    GET_YMID_R
      USE TIME_MOD, ONLY : GET_NHMSb,   GET_ELAPSED_SEC
      USE TIME_MOD, ONLY : GET_TS_CHEM, GET_DAY_OF_YEAR, GET_GMT

#     include "CMN_SIZE"  ! Size parameters
#     include "CMN_GCTM"  ! Physical constants

      ! Local variables
      LOGICAL, SAVE       :: FIRST = .TRUE.
      INTEGER             :: I, IJLOOP, J, L, N, NT, NDYSTEP
      REAL*8              :: A0, A1, A2, A3, B1, B2, B3
      REAL*8              :: LHR0, R, AHR, DEC, TIMLOC, YMID_R
      REAL*8              :: SUNTMP(MAXIJ)
      
      !=================================================================
      ! OHNO3TIME begins here!
      !=================================================================

      !  Solar declination angle (low precision formula, good enough for us):
      A0 = 0.006918
      A1 = 0.399912
      A2 = 0.006758
      A3 = 0.002697
      B1 = 0.070257
      B2 = 0.000907
      B3 = 0.000148
      R  = 2.* PI * float( GET_DAY_OF_YEAR() - 1 ) / 365.

      DEC = A0 - A1*cos(  R) + B1*sin(  R)
     &         - A2*cos(2*R) + B2*sin(2*R)
     &         - A3*cos(3*R) + B3*sin(3*R)

      LHR0 = int(float( GET_NHMSb() )/10000.)

      ! Only do the following at the start of a new day
      IF ( FIRST .or. GET_GMT() < 1e-5 ) THEN 
      
         ! Zero arrays
         TTDAY(:,:) = 0d0
         TCOSZ(:,:) = 0d0
         COSZM(:,:) = 0d0

         ! NDYSTEP is # of chemistry time steps in this day
         NDYSTEP = ( 24 - INT( GET_GMT() ) ) * 60 / GET_TS_CHEM()         

         ! NT is the elapsed time [s] since the beginning of the run
         NT = GET_ELAPSED_SEC()

         ! Loop forward through NDYSTEP "fake" timesteps for this day 
         DO N = 1, NDYSTEP
            
            ! Zero SUNTMP array
            SUNTMP(:) = 0d0

            ! Loop over surface grid boxes
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, YMID_R, IJLOOP, TIMLOC, AHR )
            DO J = 1, JJPAR

               ! Grid box latitude center [radians]
               YMID_R = GET_YMID_R( J )

            DO I = 1, IIPAR

               ! Increment IJLOOP
               IJLOOP = ( (J-1) * IIPAR ) + I
               TIMLOC = real(LHR0) + real(NT)/3600.0 + GET_XMID(I)/15.0
         
               DO WHILE (TIMLOC .lt. 0)
                  TIMLOC = TIMLOC + 24.0
               ENDDO

               DO WHILE (TIMLOC .gt. 24.0)
                  TIMLOC = TIMLOC - 24.0
               ENDDO

               AHR = abs(TIMLOC - 12.) * 15.0 * PI_180

            !===========================================================
            ! The cosine of the solar zenith angle (SZA) is given by:
            !     
            !  cos(SZA) = sin(LAT)*sin(DEC) + cos(LAT)*cos(DEC)*cos(AHR) 
            !                   
            ! where LAT = the latitude angle, 
            !       DEC = the solar declination angle,  
            !       AHR = the hour angle, all in radians. 
            !
            ! If SUNCOS < 0, then the sun is below the horizon, and 
            ! therefore does not contribute to any solar heating.  
            !===========================================================

               ! Compute Cos(SZA)
               SUNTMP(IJLOOP) = sin(YMID_R) * sin(DEC) +
     &                          cos(YMID_R) * cos(DEC) * cos(AHR)

               ! TCOSZ is the sum of SUNTMP at location (I,J)
               ! Do not include negative values of SUNTMP
               TCOSZ(I,J) = TCOSZ(I,J) + MAX( SUNTMP(IJLOOP), 0d0 )

               ! COSZM is the peak value of SUMTMP during a day at (I,J)
               ! (rjp, bmy, 3/30/04)
               COSZM(I,J) = MAX( COSZM(I,J), SUNTMP(IJLOOP) )

               ! TTDAY is the total daylight time at location (I,J)
               IF ( SUNTMP(IJLOOP) > 0d0 ) THEN
                  TTDAY(I,J) = TTDAY(I,J) + DBLE( GET_TS_CHEM() )
               ENDIF
            ENDDO
            ENDDO
!$OMP END PARALLEL DO

            ! Increment elapsed time [sec]
            NT = NT + ( GET_TS_CHEM() * 60 )             
         ENDDO

         ! Reset first-time flag
         FIRST = .FALSE.
      ENDIF

      ! Return to calling program
      END SUBROUTINE OHNO3TIME

!------------------------------------------------------------------------------

      SUBROUTINE DEFINE_TAGGED_Hg
!
!******************************************************************************
!  Subroutine DEFINE_TAGGED_Hg defines the tagged tracer numbers for
!  anthropogenic (by geographic region) Hg0, Hg2, and HgP.  The position of 
!  Hg2 and HgP in the DEPSAV array is also computed for future use.  This 
!  routine only has to be called once at the start of the simulation. 
!  (eck, cdh, bmy, 12/15/04, 1/9/06)
!
!  NOTES:
!  (1 ) Now only define AN_Hg0, AN_Hg2, AN_HgP.  Now use ID_Hg0, ID_Hg2, and
!        ID_HgP index arrays from "tracerid_mod.f". (eck, cdh, bmy, 1/9/06)
!******************************************************************************
!
      ! References to F90 modules
      USE GRID_MOD,     ONLY : GET_XMID, GET_YMID
      USE ERROR_MOD,    ONLY : ERROR_STOP
      USE TRACER_MOD,   ONLY : N_TRACERS
      USE TRACERID_MOD, ONLY : ID_Hg0,    ID_Hg2,   ID_HgP
      USE TRACERID_MOD, ONLY : ID_Hg_tot, ID_Hg_na, ID_Hg_eu, ID_Hg_as
      USE TRACERID_MOD, ONLY : ID_Hg_rw,  ID_Hg_oc, ID_Hg_ln, ID_Hg_nt

#     include "CMN_SIZE"     ! Size parameters

      ! Local variables
      INTEGER               :: I, J
      REAL*8                :: X, Y
      CHARACTER(LEN=255)    :: LOCATION

      !=================================================================
      ! DEFINE_TAGGED_Hg begins here!
      !=================================================================

      ! Location string for ERROR_STOP
      LOCATION = 'DEFINE_TAGGED_Hg ("mercury_mod.f")'

      !---------------------------------
      ! Anthropogenic tracer indices
      !---------------------------------
!$OMP PARALLEL DO 
!$OMP+DEFAULT( SHARED ) 
!$OMP+PRIVATE( I, J, X, Y )

      ! Loop over latitudes
      DO J = 1, JJPAR

         ! Grid-box latitude [degrees]
         Y = GET_YMID( J )
         
         ! Loop over longitudes
         DO I = 1, IIPAR

            ! Grid box longitude [degrees]
            X = GET_XMID( I )

            ! North American Anthro Hg 
            IF ( ( X >= -172.5 .and. X < -17.5 )   .and. 
     &           ( Y >=   24.0 .and. Y <  88.0 ) ) THEN
               AN_Hg0(I,J) = ID_Hg0(ID_Hg_na)
               AN_Hg2(I,J) = ID_Hg2(ID_Hg_na)
               AN_HgP(I,J) = ID_HgP(ID_Hg_na)
            
            ! European Anthro Hg (1st sub-box)
            ELSE IF ( ( X >= -17.5 .and. X < 72.5 )  .and. 
     &                ( Y >=  36.0 .and. Y < 45.0 ) ) THEN
               AN_Hg0(I,J) = ID_Hg0(ID_Hg_eu)
               AN_Hg2(I,J) = ID_Hg2(ID_Hg_eu)
               AN_HgP(I,J) = ID_HgP(ID_Hg_eu)

            ! European Anthro Hg (2nd sub-box)
            ELSE IF ( ( X >= -17.5 .and. X < 172.5 )  .and. 
     &                ( Y >=  45.0 .and. Y <  88.0 ) ) THEN
               AN_Hg0(I,J) = ID_Hg0(ID_Hg_eu)
               AN_Hg2(I,J) = ID_Hg2(ID_Hg_eu)
               AN_HgP(I,J) = ID_HgP(ID_Hg_eu)

            ! Asian Anthro Hg 
            ELSE IF ( ( X >= 70.0 .and. X < 152.5 )  .and.
     &                ( Y >=  8.0 .and. Y <  45.0 ) ) THEN
               AN_Hg0(I,J) = ID_Hg0(ID_Hg_as)
               AN_Hg2(I,J) = ID_Hg2(ID_Hg_as)
               AN_HgP(I,J) = ID_HgP(ID_Hg_as)

            ! Rest-of-world Anthro Hg 
            ELSE
               AN_Hg0(I,J) = ID_Hg0(ID_Hg_rw)
               AN_Hg2(I,J) = ID_Hg2(ID_Hg_rw)
               AN_HgP(I,J) = ID_HgP(ID_Hg_rw)

            ENDIF
         ENDDO
      ENDDO
!$OMP END PARALLEL DO

      !=================================================================
      ! Error check tracers: make sure they are not zero
      ! since this can cause array-out-of-bounds errors
      !=================================================================

      ! Tagged Hg0
      IF ( ANY( AN_Hg0 == 0 ) ) THEN
         CALL ERROR_STOP( 'AN_Hg0 tracers are undefined!', LOCATION )   
      ENDIF
      
      ! Tagged Hg2
      IF ( ANY( AN_Hg2 == 0 ) ) THEN
         CALL ERROR_STOP( 'AN_Hg2 tracers are undefined!', LOCATION )   
      ENDIF

      ! Tagged Hg2
      IF ( ANY( AN_HgP == 0 ) ) THEN
         CALL ERROR_STOP( 'AN_HgP tracer are undefined!', LOCATION )   
      ENDIF

      !---------------------------------
      ! Error check # of tracers
      !---------------------------------
      IF ( N_TRACERS < 21 ) THEN 
         CALL ERROR_STOP( 'Too few Hg tagged tracers!',
     &                    'DEFINE_TAGGED_Hg ("mercury_mod.f")' )
      ENDIF

      ! Return to calling program
      END SUBROUTINE DEFINE_TAGGED_Hg

!------------------------------------------------------------------------------

      SUBROUTINE INIT_MERCURY( THIS_ANTHRO_Hg_YEAR )
!
!******************************************************************************
!  Subroutine INIT_MERCURY allocates and zeroes all module arrays.
!  (eck, cdh, sas, bmy, 12/2/04, 1/9/06)
!  
!  NOTES:
!  (1 ) Removed reference to FEMIS array.  Now also allocates and zeroes
!        the T44 array.  Added reference to CMN_DIAG.  Now references 
!        N_TRACERS from "tracer_mod.f". (bmy, 2/24/05)
!  (2 ) EHg0_an, EHg2_an, EHgP_an are now 2-D arrays.  Now modified for 
!        updated ocean mercury module. (eck, cdh, sas, bmy, 1/9/06)
!******************************************************************************
!
      ! References to F90 modules
      USE DRYDEP_MOD,   ONLY : DEPNAME,   NUMDEP
      USE ERROR_MOD,    ONLY : ALLOC_ERR, ERROR_STOP
      USE LOGICAL_MOD,  ONLY : LSPLIT,    LDRYD
      USE TRACER_MOD,   ONLY : N_TRACERS
      USE TRACERID_MOD, ONLY : N_Hg_CATS

#     include "CMN_SIZE"     ! Size parameters
#     include "CMN_DIAG"     ! ND44

      ! Arguments
      INTEGER, INTENT(IN)   :: THIS_ANTHRO_Hg_YEAR

      ! Local variables
      LOGICAL, SAVE         :: IS_INIT = .FALSE. 
      INTEGER               :: AS, N
      CHARACTER(LEN=255)    :: LOCATION

      !=================================================================
      ! INIT_MERCURY begins here!
      !=================================================================

      ! Return if we have already allocated arrays
      IF ( IS_INIT ) RETURN

      !---------------------------------------------
      ! Prior to 2/24/06
      ! This is now in "logical_mod.f"
      !! For now, define dynamic ocean flag here
      !LDYNOCEAN = .FALSE.
      !---------------------------------------------

      ! Anthropogenic Hg emissions year
      ANTHRO_Hg_YEAR = THIS_ANTHRO_Hg_YEAR

      ! Location string for ERROR_STOP
      LOCATION       = 'DEFINE_TAGGED_Hg ("mercury_mod.f")'

      !=================================================================
      ! Allocate arrays
      !=================================================================
      ALLOCATE( COSZM( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'COSZM' )
      COSZM = 0d0

      ALLOCATE( EHg0_an( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'EHg0_an' )
      EHg0_an = 0d0

      ALLOCATE( EHg2_an( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'EHg2_an' )
      EHg2_an = 0d0

      ALLOCATE( EHgP_an( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'EHgP_an' )
      EHgP_an = 0d0

      ALLOCATE( EHg0_oc( IIPAR, JJPAR, N_Hg_CATS ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'EHg0_oc' )
      EHg0_oc = 0d0

      ALLOCATE( EHg0_ln( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'EHg0_ln' )
      EHg0_ln = 0d0

      ALLOCATE( EHg0_nt( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'EHg0_nt' )
      EHg0_nt = 0d0

      ALLOCATE( TCOSZ( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'TCOSZ' )
      TCOSZ = 0d0

      ALLOCATE( TTDAY( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'TTDAY' )
      TTDAY = 0d0

      ! Allocate ZERO_DVEL if drydep is turned off 
      IF ( .not. LDRYD ) THEN
         ALLOCATE( ZERO_DVEL( IIPAR, JJPAR ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'ZERO_DVEL' )
         ZERO_DVEL = 0d0
      ENDIF

      !=================================================================
      ! Allocate & initialize arrays for tagged tracer indices
      !=================================================================
      IF ( LSPLIT ) THEN

         ! Tracer indices for tagged anthro regions
         ALLOCATE( AN_Hg0( IIPAR, JJPAR ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'AN_Hg0' )
         AN_Hg0 = 0d0

         ALLOCATE( AN_Hg2( IIPAR, JJPAR ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'AN_Hg2' )
         AN_Hg2 = 0d0

         ALLOCATE( AN_HgP( IIPAR, JJPAR ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'AN_HgP' )
         AN_HgP = 0d0

         ! Define the tagged tracer indices 
         CALL DEFINE_TAGGED_Hg

      ENDIF

      !=================================================================
      ! Locate the drydep species w/in the DEPSAV array (for ND44)
      !=================================================================

      ! Initialize flags
      DRYHg2 = 0d0 
      DRYHgP = 0d0 

      ! If drydep is turned on ...
      IF ( LDRYD ) THEN
         
         ! Loop over drydep species
         DO N = 1, NUMDEP

            ! Locate by DEPNAME
            SELECT CASE ( TRIM( DEPNAME(N) ) )
               CASE( 'Hg2' )
                  DRYHg2 = N
               CASE( 'HgP' )
                  DRYHgP = N
               CASE DEFAULT
                  ! nothing
            END SELECT
         ENDDO

         ! Temporary array for ND44 diagnosti
         !----------------------------------------------------------------
         ! Prior to 2/27/06:
         ! For now, comment this out so that MERCURY MENU can be placed
         ! above the DIAGNOSTIC MENU in "input.geos" (bmy, 2/27/06)
         !IF ( ND44 > 0 ) THEN
         !----------------------------------------------------------------
            ALLOCATE( T44( IIPAR, JJPAR, LLTROP, N_TRACERS ), STAT=AS )
            IF ( AS /= 0 ) CALL ALLOC_ERR( 'T44' )
            T44 = 0e0
         !----------------------------------------------------------------
         ! Prior to 2/27/06:
         ! For now, comment this out so that MERCURY MENU can be placed
         ! above the DIAGNOSTIC MENU in "input.geos" (bmy, 2/27/06)
         !ENDIF
         !----------------------------------------------------------------
      ENDIF

      ! Reset IS_INIT, since we have already allocated arrays
      IS_INIT = .TRUE.

      ! Return to calling program
      END SUBROUTINE INIT_MERCURY

!------------------------------------------------------------------------------

      SUBROUTINE CLEANUP_MERCURY
!
!******************************************************************************
!  Subroutine CLEANUP_MERCURY deallocates all module arrays.
!  (eck, bmy, 12/6/04, 2/24/05)
!
!  NOTES:
!  (1 ) Now deallocate MLD, NPP, RAD (sas, bmy, 1/18/05)
!  (2 ) Now deallocate T44 (bmy, 2/24/05)
!******************************************************************************
!
      !=================================================================
      ! CLEANUP_MERCURY begins here!
      !=================================================================
      IF ( ALLOCATED( AN_Hg0  ) ) DEALLOCATE( AN_Hg0  )  
      IF ( ALLOCATED( AN_Hg2  ) ) DEALLOCATE( AN_Hg2  )      
      IF ( ALLOCATED( AN_HgP  ) ) DEALLOCATE( AN_HgP  )      
      IF ( ALLOCATED( COSZM   ) ) DEALLOCATE( COSZM   )      
      IF ( ALLOCATED( EHg0_an ) ) DEALLOCATE( EHg0_an )
      IF ( ALLOCATED( EHg2_an ) ) DEALLOCATE( EHg2_an )
      IF ( ALLOCATED( EHgP_an ) ) DEALLOCATE( EHgP_an )
      IF ( ALLOCATED( EHg0_oc ) ) DEALLOCATE( EHg0_oc )
      IF ( ALLOCATED( EHg0_ln ) ) DEALLOCATE( EHg0_ln )
      IF ( ALLOCATED( EHg0_nt ) ) DEALLOCATE( EHg0_nt )
      IF ( ALLOCATED( TCOSZ   ) ) DEALLOCATE( TCOSZ   )
      IF ( ALLOCATED( T44     ) ) DEALLOCATE( T44     )
      IF ( ALLOCATED( TTDAY   ) ) DEALLOCATE( TTDAY   )

      ! Return to calling program
      END SUBROUTINE CLEANUP_MERCURY

!------------------------------------------------------------------------------

      ! End of module
      END MODULE MERCURY_MOD
