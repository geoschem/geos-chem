! $Id: mercury_mod.f,v 1.1 2004/12/16 16:52:45 bmy Exp $
      MODULE MERCURY_MOD
!
!******************************************************************************
!  Module MERCURY_MOD contains variables and routines for the GEOS-CHEM 
!  mercury simulation. (eck, bmy, 12/14/04)
!
!  Module Variables:
!  ============================================================================
!  (1 ) AN_Hg0  (INTEGER) : Tracer index array for tagged anthro Hg(0)  regions
!  (2 ) AN_Hg2  (INTEGER) : Tracer index array for tagged anthro Hg(II) regions
!  (3 ) AN_HgP  (INTEGER) : Tracer index array for tagged anthro HgP    regions
!  (4 ) COSZM   (REAL*8 ) : Max daily value of COS( S. Z. Angle ) [unitless]
!  (5 ) DRYHg2  (INTEGER) : Index for Hg2 in DEPSAV array (drydep freq)
!  (6 ) DRYHgP  (INTEGER) : Index for HgP in DEPSAV array (drydep freq)
!  (7 ) EHg0_an (REAL*8 ) : Anthropogenic Hg0 emissions [kg/s]
!  (8 ) EHg2_an (REAL*8 ) : Anthropogenic Hg2 emissions [kg/s]
!  (9 ) EHgP_an (REAL*8 ) : Anthropogenic HgP emissions [kg/s]
!  (10) EHg0_oc (REAL*8 ) : Hg0 emissions from oceans [kg/s]
!  (11) EHg0_ln (REAL*8 ) : Re-emissions of Hg0 from land [kg/s]
!  (12) EHg0_nt (REAL*8 ) : Hg0 emissions from natural land sources [kg/s] 
!  (13) FEMIS   (REAL*8 ) : Fraction of emissions in layer [unitless]
!  (14) LN_Hg0  (INTEGER) : Tracer index for tagged land re-emission of Hg(0)
!  (15) LN_Hg2  (INTEGER) : Tracer index for tagged land re-emission of Hg(II)
!  (16) LN_HgP  (INTEGER) : Tracer index for tagged land re-emission of HgP
!  (17) NT_Hg0  (INTEGER) : Tracer index for natural land emissions of Hg(0)
!  (18) NT_Hg2  (INTEGER) : Tracer index for natural land emissions of Hg(II)
!  (19) NT_HgP  (INTEGER) : Tracer index for natural land emissions of HgP
!  (20) OC_Hg0  (INTEGER) : Tracer index for oceanic emissions of Hg(0)
!  (21) OC_Hg2  (INTEGER) : Tracer index for oceanic emissions of Hg(II)
!  (22) OC_HgP  (INTEGER) : Tracer index for oceanic emissions of HgP 
!  (23) TCOSZ   (REAL*8 ) : Sum of COS( Solar Zenith Angle ) [unitless]
!  (24) TTDAY   (REAL*8 ) : Total daylight time at location (I,J) [minutes]
!
!  Module Routines:
!  ============================================================================
!  (1 ) CHEMMERCURY       : Chemistry routine for Hg
!  (2 ) CHEM_Hg0_Hg2      : Chemistry for Hg0, Hg2 and drydep of Hg2
!  (3 ) RXN_Hg0_Hg2       : Conversion of Hg(0) --> Hg(II) via reduction rxn
!  (4 ) RXN_Hg0           : Conversion of Hg(0) --> Hg(II) via oxidation rxns
!  (5 ) RXN_Hg2_DRYD      : Prod of Hg(II) from Hg(0); also Hg(II) drydep loss
!  (6 ) RXN_Hg2           : Prod of Hg(II) from Hg(0) 
!  (7 ) CHEM_HGP          : Chemistry (via drydep loss) for HgP
!  (8 ) RXN_HgP_DRYD      : Loss of HgP via drydep
!  (9 ) EMISSMERCURY      : Emissions for Hg
!  (10) SRCHG0            : Applies emissions of Hg0
!  (11) SRCHG2            : Applies emissions of Hg2
!  (12) SRCHGP            : Applies emissions of HgP
!  (13) COMPUTE_FEMIS     : Compute fraction of emissions in layer
!  (14) MERCURY_READYR    : Reads mercury emissions and converts to [kg/s]
!  (15) GET_LWC           : Computes liquid water content as a function of T
!  (16) GET_VCLDF         : Computes volume cloud fraction for SO2 chemistry 
!  (17) GET_O3            : Returns monthly mean O3 values
!  (18) GET_OH            : Returns monthly mean OH values (w/ diurnal scaling)
!  (19) OHNO3TIME         : Computes diurnal scaling for monthly mean OH
!  (20) DEFINE_TAGGED_Hg  : Defines tracer number for tagged Hg tracers 
!  (20) INIT_MERCURY      : Allocates and zeroes module arrays
!  (21) CLEANUP_MERCURY   : Deallocates module arrays
!
!  GEOS-CHEM modules referenced by mercury_mod.f
!  ============================================================================
!  (1 ) bpch2_mod.f       : Module containing routines for binary pch file I/O
!  (2 ) comode_mod.f      : Module containing SMVGEAR allocatable arrays
!  (3 ) dao_mod.f         : Module containing arrays for DAO met fields
!  (4 ) diag_mod.f        : Module containing GEOS-CHEM diagnostic arrays
!  (5 ) drydep_mod.f      : Module containing GEOS-CHEM dry deposition routines
!  (6 ) error_mod.f       : Module containing NaN, other error check routines
!  (7 ) global_o3_mod.f   : Module containing routines to read 3-D O3 field
!  (8 ) global_oh_mod.f   : Module containing routines to read 3-D OH field
!  (9 ) grid_mod.f        : Module containing horizontal grid information
!  (10) logical_mod.f     : Module containing GEOS-CHEM logical switches
!  (11) pressure_mod.f    : Module containing routines to compute P(I,J,L)
!  (12) time_mod.f        : Module containing routines to compute date & time
!  (13) tracer_mod.f      : Module containing GEOS-CHEM tracer array STT etc.
!  (14) tracerid_mod.f    : Module containing pointers to tracers & emissions
!  (15) transfer_mod.f    : Module containing routines to cast & resize arrays
!
!  Nomenclature: 
!  ============================================================================
!  (1 ) Hg(0)  a.k.a. Hg0 : Elemental   mercury
!  (2 ) Hg(II) a.k.a. Hg2 : Divalent    mercury
!  (3 ) HgP               : Particulate mercury
!
!  Mercury Tracers (1-3 are always defined; 4-21 are defined for tagged runs)
!  ============================================================================
!  (1 ) Hg(0)             : Hg(0)  - total tracer
!  (2 ) Hg(II)            : Hg(II) - total tracer 
!  (3 ) HgP               : HgP    - total tracer
!  -----------------------+----------------------------------------------------
!  (4 ) Hg(0)_an_na       : Hg(0)  - North American anthropogenic
!  (5 ) Hg(0)_an_eu       : Hg(0)  - European Anthropogenic
!  (6 ) Hg(0)_an_as       : Hg(0)  - Asian anthropogenic
!  (7 ) Hg(0)_an_rw       : Hg(0)  - Rest of World Anthropogenic
!  (8 ) Hg(0)_oc          : Hg(0)  - Ocean emission
!  (9 ) Hg(0)_ln          : Hg(0)  - Land reemission
!  (10) Hg(0)_nt          : Hg(0)  - Land natural emission
!  (11) Hg(II)_an_na      : Hg(II) - North American anthropogenic
!  (12) Hg(II)_an_eu      : Hg(II) - European Anthropogenic
!  (13) Hg(II)_an_as      : Hg(II) - Asian anthropogenic
!  (14) Hg(II)_an_rw      : Hg(II) - Rest of World Anthropogenic
!  (15) Hg(II)_oc         : Hg(II) - Ocean emission
!  (16) Hg(II)_ln         : Hg(II) - Land reemission
!  (17) Hg(II)_nt         : Hg(II) - Land natural emission
!  (18) HgP_an_na         : HgP    - North American anthropogenic
!  (19) HgP_an_eur        : HgP    - European anthropogenic
!  (20) HgP_an_as         : HgP    - Asian anthropogenic
!  (21) HgP_an_rw         : HgP    - Rest of world anthropogenic
!  (22) HgP_oc            : HgP    - Ocean emission        (FOR FUTURE USE)
!  (23) HgP_ln            : HgP    - Land reemission       (FOR FUTURE USE)
!  (24) HgP_nt            : HgP    - Land natural emission (FOR FUTURE USE)
!
!  References
!  ============================================================================
!  (1 ) Hall, 1995
!  (2 ) Sommar et al. 2001
!
!  NOTES:
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
      PUBLIC :: EMISSMERCURY

      !=================================================================
      ! MODULE VARIABLES
      !=================================================================

      ! Scalars
      INTEGER              :: DRYHg2, DRYHgP
      INTEGER, PARAMETER   :: Hg_LEVELS = 3
      INTEGER              :: OC_Hg0, LN_Hg0, NT_Hg0
      INTEGER              :: OC_Hg2, LN_Hg2, NT_Hg2
      INTEGER              :: OC_HgP, LN_HgP, NT_HgP
      REAL*8,  PARAMETER   :: SMALLNUM  = 1d-32

      ! Arrays
      INTEGER, ALLOCATABLE :: AN_Hg0(:,:)
      INTEGER, ALLOCATABLE :: AN_Hg2(:,:)
      INTEGER, ALLOCATABLE :: AN_HgP(:,:)
      REAL*8,  ALLOCATABLE :: COSZM(:,:)
      REAL*8,  ALLOCATABLE :: EHg0_an(:,:,:)
      REAL*8,  ALLOCATABLE :: EHg2_an(:,:,:)
      REAL*8,  ALLOCATABLE :: EHgP_an(:,:,:)
      REAL*8,  ALLOCATABLE :: EHg0_oc(:,:)
      REAL*8,  ALLOCATABLE :: EHg0_ln(:,:)
      REAL*8,  ALLOCATABLE :: EHg0_nt(:,:)
      REAL*8,  ALLOCATABLE :: FEMIS(:,:,:)
      REAL*8,  ALLOCATABLE :: TCOSZ(:,:)
      REAL*8,  ALLOCATABLE :: TTDAY(:,:)

      !=================================================================
      ! MODULE ROUTINES -- follow below the "CONTAINS" statement
      !=================================================================
      CONTAINS

!------------------------------------------------------------------------------

      SUBROUTINE CHEMMERCURY
!
!******************************************************************************
!  Subroutine CHEMMERCURY is the driver routine for mercury chemistry
!  in the GEOS-CHEM module. (eck, bmy, 12/6/04)
!
!  NOTES:
!******************************************************************************
!
      ! References to F90 modulees
      USE DRYDEP_MOD,    ONLY : DEPSAV
      USE ERROR_MOD,     ONLY : DEBUG_MSG
      USE GLOBAL_O3_MOD, ONLY : GET_GLOBAL_O3
      USE GLOBAL_OH_MOD, ONLY : GET_GLOBAL_OH
      USE LOGICAL_MOD,   ONLY : LPRT
      USE TIME_MOD,      ONLY : GET_MONTH, ITS_A_NEW_MONTH

#     include "CMN_SIZE"      ! Size parameters

      ! Local variables
      LOGICAL, SAVE          :: FIRST = .TRUE.
      INTEGER                :: N, MONTH

      !=================================================================
      ! CHEMMERCURY begins here!
      !=================================================================

      ! First-time initialization
      IF ( FIRST ) THEN

         ! Allocate arrays & locate drydep indices (if not already done)
         CALL INIT_MERCURY

         ! Reset first-time flag
         FIRST = .FALSE.
      ENDIF

      !=================================================================
      ! Read monthly mean OH and O3 fields
      !=================================================================
      IF ( ITS_A_NEW_MONTH() ) THEN 

         ! Get the current month
         MONTH = GET_MONTH()

         ! Read monthly mean OH and O3 from disk
         CALL GET_GLOBAL_OH( MONTH )
         CALL GET_GLOBAL_O3( MONTH )

      ENDIF
      
      !=================================================================
      ! Perform chemistry on Hg tracers 
      !=================================================================
      
      ! Compute diurnal scaling for OH
      CALL OHNO3TIME

      ! Hg0 and Hg2 
      IF ( DRYHg2 > 0 ) THEN
         CALL CHEM_Hg0_Hg2( DEPSAV(:,:,DRYHg2) )
         IF ( LPRT ) CALL DEBUG_MSG( 'CHEMMERCURY: a CHEM_Hg0_Hg2' )
      ENDIF

      ! HgP
      IF ( DRYHgP > 0 ) THEN
         CALL CHEM_HgP( DEPSAV(:,:,DRYHgP) )
         IF ( LPRT ) CALL DEBUG_MSG( 'CHEMMERCURY: a CHEM_HgP' )
      ENDIF

      ! Return to calling program
      END SUBROUTINE CHEMMERCURY

!------------------------------------------------------------------------------

      SUBROUTINE CHEM_Hg0_Hg2( V_DEP_Hg2 )
!
!******************************************************************************
!  Subroutine CHEM_Hg0_Hg2 is the chemistry subroutine for the oxidation/
!  reduction of Hg0/Hg(II).  For tracers with dry deposition, the loss rate 
!  of dry dep is combined in the chemistry loss term. (eck, bmy, 12/6/04)
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
!******************************************************************************
!
      ! References to F90 modules
      USE DIAG_MOD,     ONLY : AD03_Hg2_Hg0, AD03_Hg2_O3, AD03_Hg2_OH
      USE DAO_MOD,      ONLY : T, AD
      USE LOGICAL_MOD,  ONLY : LSPLIT
      USE TIME_MOD,     ONLY : GET_TS_CHEM
      USE TRACER_MOD,   ONLY : STT
      USE TRACERID_MOD, ONLY : IDTHg0, IDTHg2

#     include "CMN_SIZE"     ! Size parameters
#     include "CMN_O3"       ! XNUMOL
#     include "CMN_DIAG"     ! ND03, LD03

      ! Arguments
      REAL*8, INTENT(IN)     :: V_DEP_Hg2(IIPAR,JJPAR)

      ! Local variables
      INTEGER                :: I, J, L, N
      REAL*8                 :: AREA_CM2,    C_O3,        C_OH
      REAL*8                 :: DRYDEP,      DTCHEM,      E_RKT
      REAL*8                 :: E_K1,        E_K2,        FA
      REAL*8                 :: FC,          K1,          K2          
      REAL*8                 :: KD,          KO,          KR
      REAL*8                 :: KT,          LOST_Hg0,    LOST_Hg0_an
      REAL*8                 :: LOST_Hg0_oc, LOST_Hg0_ln, LOST_Hg0_nt
      REAL*8                 :: LOST_Hg0_O3, LOST_Hg0_OH, LWC
      REAL*8                 :: NK1,         NK2,         RKT

      ! K for reaction Hg0 + O3 in cm3 molec-1 s-1 (Source: Hall '95)
      REAL*8, PARAMETER      :: K       = 3.0d-20

      ! K for reaction Hg2 + OH in cm3 s-1 
      REAL*8, PARAMETER      :: K_HG_OH = 8.7d-14

      ! Henry's Law constant for Hg2
      REAL*8, PARAMETER      :: HL      = 1.4D6

      ! R
      REAL*8, PARAMETER      :: R       = 8.2D-2

      ! K for reduction (scaled to budget and OH conc)
      ! testing in progress
      !-----------------------------------------------------------------
      ! Prior to 12/10/04:
      ! Set KRED to zero for testing.  We still need to pin down the
      ! value of the reduction rate constant (eck, bmy, 12/10/04)
      !REAL*8, PARAMETER      :: KRED    = 1.5D-17
      !-----------------------------------------------------------------
      REAL*8, PARAMETER      :: KRED    = 0d0
     
      ! External functions
      REAL*8,  EXTERNAL      :: BOXVL

      !=================================================================
      ! CHEM_Hg0_Hg2 begins here!
      !=================================================================

      ! Chemistry timestep [s]
      DTCHEM = GET_TS_CHEM() * 60d0

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I,           J,           L,           LOST_Hg0    )
!$OMP+PRIVATE( LOST_Hg0_an, LOST_Hg0_oc, LOST_Hg0_ln, LOST_Hg0_nt )
!$OMP+PRIVATE( LOST_Hg0_O3, LOST_Hg0_OH, C_O3,        C_OH        )
!$OMP+PRIVATE( FC,          LWC,         FA,          K1          )
!$OMP+PRIVATE( K2,          KO,          KR,          KT          )
!$OMP+PRIVATE( RKT,         E_RKT,       E_K1,        E_K2        )
!$OMP+PRIVATE( NK1,         NK2                                   )
      DO L = 1, LLPAR
      DO J = 1, JJPAR
      DO I = 1, IIPAR
        
         ! Initialize 
         LOST_Hg0    = 0d0
         LOST_Hg0_an = 0d0
         LOST_Hg0_oc = 0d0
         LOST_Hg0_ln = 0d0
         LOST_Hg0_nt = 0d0
         LOST_Hg0_O3 = 0d0
         LOST_Hg0_OH = 0d0
   
         ! Monthly mean O3 and OH concentrations [molec/cm3]
         C_O3        = GET_O3(I,J,L)
         C_OH        = GET_OH(I,J,L)
                
         ! Get fraction of cloud in the gridbox [unitless]
         FC          = GET_VCLDF(I,J,L)

         ! Get liquid water content in cloudy area of gridbox (m3/m3)
         LWC         = GET_LWC( T(I,J,L) ) * FC 

         ! Define fraction of Hg(II) which is in air
         FA          = ( HL * R * T(I,J,L) * LWC )
         FA          = FA / ( 1d0 + FA )

         ! Define K's for the oxidation reactions
         K1          = K       * C_O3
         K2          = K_HG_OH * C_OH
         KO          = K1      + K2
         
         ! Define K for the reduction reaction.  Include the fraction of 
         ! Hg(II) in air within the Kreduction & scale to OH concentration
         KR          = KRED * FA * C_OH 
         
         ! Total rate constant: Oxidation + Reduction
         KT          = KO + KR
         RKT         = KT * DTCHEM
            
         ! Compute the exponential terms for use below
         E_RKT       = EXP( -RKT         )
         E_K1        = EXP( -K1 * DTCHEM )
         E_K2        = EXP( -K2 * DTCHEM )

         !==============================================================
         ! Hg(0) Chemistry: Conversion from Hg(0) to Hg(II)
         !
         ! CASE 1: REDUCTION REACTION
         ! -------------------------------------------------------------
         ! Aqueous chem occurs if Kr>0, it's cloudy, and T > -15 C.
         ! In this case we have the reaction:
         !
         !    New Hg(0) =  ( Old Hg(0)     *   EXP( -KT )         ) + 
         !                 ( Old Hg(II)/KT * ( 1d0 - EXP( -KT ) ) )
         !
         ! where: KT         = K of the reduction rxn * Temperature
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
         !    New Hg(0) = Old_Hg(0) * EXP( -K1 * T ) * EXP( -K2 * T )
         !
         ! The amt of Hg(0) lost in these rxns rxn becomes Hg(II):
         !
         !    Hg(0) converted to Hg(II) = Old Hg(0) - New Hg(0)
         !==============================================================
         IF ( KR > 0d0 .and. FC > 0d0 .and. T(I,J,L) > 258d0 ) THEN
                     
            !--------------------------------------------
            ! CASE 1: Total Hg(II) tracer, aq chem
            !--------------------------------------------

            ! Conversion of total Hg(0) --> Hg(II) 
            CALL RXN_Hg0_Hg2( I,      J,   L,     IDTHG0, 
     &                        IDTHG2, RKT, E_RKT, LOST_Hg0 )

            !--------------------------------------------
            ! CASE 1: Tagged Hg(II) tracers, aq chem
            !--------------------------------------------
            IF ( LSPLIT ) THEN
               
               ! Conversion of tagged anthro Hg(0) --> Hg(II) 
               CALL RXN_Hg0_Hg2( I,           J,   L,     AN_Hg0(I,J),
     &                           AN_Hg2(I,J), RKT, E_RKT, LOST_Hg0_an )

               ! Conversion of tagged oceanic Hg(0) --> Hg(II) 
               CALL RXN_Hg0_Hg2( I,           J,   L,     OC_Hg0, 
     &                           OC_Hg2,      RKT, E_RKT, LOST_Hg0_oc )

               ! Conversion of tagged land re-emission Hg(0) --> Hg(II) 
               CALL RXN_Hg0_Hg2( I,           J,   L,     LN_Hg0, 
     &                           LN_Hg2,      RKT, E_RKT, LOST_Hg0_ln )

               ! Conversion of tagged natural source Hg(0) --> Hg(II)
               CALL RXN_Hg0_Hg2( I,           J,   L,     NT_Hg0, 
     &                           NT_Hg2,      RKT, E_RKT, LOST_Hg0_nt )
            ENDIF

         ELSE

            !--------------------------------------------
            ! CASE 2: Total Hg(II) tracer, non-aq chem
            !--------------------------------------------

            ! Total Hg(0) loss by oxidation rxn
            CALL RXN_Hg0( I, J, L, IDTHg0, E_K1, E_K2, LOST_Hg0 )

            !--------------------------------------------
            ! CASE 2: Tagged Hg(II) tracers, non-aq chem
            !--------------------------------------------
            IF ( LSPLIT ) THEN

               ! Tagged anthro Hg(0) loss by oxidation rxn
               CALL RXN_Hg0( I,           J,    L, 
     &                       AN_Hg0(I,J), E_K1, E_K2, LOST_Hg0_an )

               ! Tagged oceanic Hg(0) loss by oxidation rxn
               CALL RXN_Hg0( I,           J,    L,
     &                       OC_Hg0,      E_K1, E_K2, LOST_Hg0_oc )

               ! Tagged land re-emission Hg(0) loss by oxidation rxn
               CALL RXN_Hg0( I,           J,    L,
     &                       LN_Hg0,      E_K1, E_K2, LOST_Hg0_ln )

               ! Tagged natural source Hg(0) loss by oxidation rxn
               CALL RXN_Hg0( I,           J,    L,
     &                       NT_Hg0,      E_K1 ,E_K2, LOST_Hg0_nt )
            ENDIF

         ENDIF 

         !==============================================================
         ! Compute Hg(II) production from OH and O3 rxns for diagnostic
         !==============================================================
         IF ( LOST_Hg0 > 0d0 ) THEN

            ! This is for the diagnostics OH and O3 prod of Hg(II).
            ! They are messed up a little since adding the reduction 
            ! reaction haven't fixed them yet.  
            ! Define new k's here to avoid NaN error in denominator.
            ! (eck, 12/7/04)
            NK1      = MAX( K1, SMALLNUM )
            NK2      = MAX( K2, SMALLNUM )

            ! Production of Hg(II) from O3 rxn [kg]
            LOST_Hg0_O3 = ( NK1 / ( NK1 + NK2 ) ) * LOST_Hg0

            ! Production of Hg(II) from OH rxn [kg]
            LOST_Hg0_OH = ( NK2 / ( NK1 + NK2 ) ) * LOST_Hg0

         ENDIF

         !==============================================================
         ! Hg(II) chemistry: Conversion from Hg(0) and drydep loss
         !
         ! CASE 1: SURFACE LEVEL
         ! -------------------------------------------------------------
         ! At the surface we have both dry deposition of Hg(II) plus
         ! conversion of Hg(0) into Hg(II).  In this case we use the
         ! following rxn:
         !
         !     New Hg(II) = ( Old Hg(II)         *       EXP( -KT ) ) +
         !                  ( Converted Hg(0)/KT * ( 1 - EXP( -KT ) )
         !               
         ! Where: KT               = DTCHEM * Drydep Vel of Hg(II)
         !       "Converted Hg(0)" = Amt of Hg(0) that became Hg(II)
         !                           archived from the Hg(0) rxns above
         !
         ! CASE 2: ALL OTHER LEVELS
         ! --------------------------------------------------------------
         ! At levels higher than the surface, we do not have drydep.
         ! Therefore we only have conversion of Hg(0) into Hg(II), and
         ! we use this rxn:
         !
         !     New Hg(II) = Old Hg(II) + Converted Hg(0)
         !
         !
         ! NOTE: ND44 diagnostics are archived in RXN_Hg2_DRYD.
         !==============================================================
         IF ( L == 1 ) THEN

            ! Compute exponential terms from drydep rxns for use below
            RKT      = V_DEP_Hg2(I,J) * DTCHEM
            E_RKT    = EXP( -RKT )
            
            !--------------------------------------------
            ! CASE 1: Total Hg(II) tracer, surface
            !--------------------------------------------

            ! Compute new Hg(II) concentration at surface
            CALL RXN_Hg2_DRYD( I,   J,     1,        IDTHG2,   
     &                         RKT, E_RKT, LOST_Hg0, DTCHEM )
          
            !--------------------------------------------
            ! CASE 2: Tagged Hg(II ) tracers: surface
            !--------------------------------------------
            IF ( LSPLIT ) THEN

               ! Compute new conc of tagged anthro Hg(II) at surface
               CALL RXN_Hg2_DRYD( I,   J,     1,           AN_Hg2(I,J), 
     &                            RKT, E_RKT, LOST_Hg0_an, DTCHEM )

               ! Compute new conc of tagged oceanic Hg(II) at surface
               CALL RXN_Hg2_DRYD( I,   J,     1,           OC_Hg2,      
     &                            RKT, E_RKT, LOST_Hg0_oc, DTCHEM )

               ! Compute new conc of land re-emission Hg(II) at surface
               CALL RXN_Hg2_DRYD( I,   J,     1,           LN_Hg2,      
     &                            RKT, E_RKT, LOST_Hg0_ln, DTCHEM )

               ! Compute new conc of tagged natural source Hg(II) at surface
               CALL RXN_Hg2_DRYD( I,   J,     1,           NT_Hg2,      
     &                            RKT, E_RKT, LOST_Hg0_nt, DTCHEM )

            ENDIF
                   
         ELSE

            !--------------------------------------------
            ! CASE 2: Total Hg(II) tracer, other levels
            !--------------------------------------------

            ! Compute new concentration of total Hg(II) for other levels [kg]
            CALL RXN_Hg2( I, J, L, IDTHg2, LOST_Hg0 )
            
            !--------------------------------------------
            ! CASE 2: Tagged Hg(II) tracers, other levels
            !--------------------------------------------
            IF ( LSPLIT ) THEN

               ! Compute new concentration of tagged Hg(II) tracers
               ! at levels higher than the surface [kg]
               CALL RXN_Hg2( I, J, L, AN_Hg2(I,J), LOST_Hg0_an )
               CALL RXN_Hg2( I, J, L, OC_Hg2,      LOST_Hg0_oc )
               CALL RXN_Hg2( I, J, L, LN_Hg2,      LOST_Hg0_ln )
               CALL RXN_Hg2( I, J, L, NT_Hg2,      LOST_Hg0_nt )

            ENDIF

         ENDIF
     
         !==============================================================
         ! ND03 diagnostic: Hg(II) production [kg]
         !==============================================================
         IF ( ND03 > 0 .and. L <= LD03 ) THEN
            AD03_Hg2_Hg0(I,J,L) = AD03_Hg2_Hg0(I,J,L) + LOST_Hg0
            AD03_Hg2_OH(I,J,L)  = AD03_Hg2_OH(I,J,L)  + LOST_Hg0_OH
            AD03_Hg2_O3(I,J,L)  = AD03_Hg2_O3(I,J,L)  + LOST_Hg0_O3
         ENDIF

      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO
      
      ! Return to calling program
      END SUBROUTINE CHEM_HG0_HG2

!------------------------------------------------------------------------------

      SUBROUTINE RXN_Hg0_Hg2( I,     J,   L,     N_Hg0, 
     &                        N_Hg2, RKT, E_RKT, LOST_Hg0 )
!
!******************************************************************************
!  Subroutine RXN_Hg0_Hg2 computes the conversion of Hg(0) to Hg(II) via
!  an aqueous chemistry reduction reaction.  The formula used below is a s
!  solution of the 2-box model equation. (eck, bmy, 12/14/04)
! 
!  Arguments as Input:
!  ============================================================================
!  (1-3) I, J, L  (INTEGER) : GEOS-CHEM lon, lat, alt grid box indices
!  (4  ) N_Hg0    (INTEGER) : Index for Hg(0)  total or tagged tracers
!  (5  ) N_Hg2    (INTEGER) : Index for Hg(II) total or tagged tracers
!  (6  ) RKT      (REAL*8 ) : Value of R * k * Temp for this rxn [unitless]
!  (7  ) E_RKT    (REAL*8 ) : Value of EXP( - RKT ) [unitless]
!
!  Arguments as Output:
!  ============================================================================
!  (8  ) LOST_Hg0 (REAL*8 ) : Loss term: Hg(0) before - Hg(0) after [kg]
!
!  NOTES:
!******************************************************************************
!
      ! References to F90 modules
      USE TRACER_MOD, ONLY : STT

      ! Arguments
      INTEGER, INTENT(IN)  :: I, J, L, N_Hg0, N_Hg2
      REAL*8,  INTENT(IN)  :: RKT, E_RKT
      REAL*8,  INTENT(OUT) :: LOST_Hg0

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
      NEW_Hg0          = ( OLD_Hg0      *   E_RKT         ) + 
     &                   ( OLD_Hg2/RKT  * ( 1d0 - E_RKT ) )

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

      SUBROUTINE RXN_Hg2_DRYD( I,   J,     L,        N,    
     &                         RKT, E_RKT, LOST_Hg0, DTCHEM )
!
!******************************************************************************
!  Subroutine RXN_Hg0_DRYD computes the new concentration of Hg(II) from the
!  converted Hg(0) plus the drydep of Hg(II). (eck, bmy, 12/14/04)
! 
!  Arguments as Input:
!  ============================================================================
!  (1-3) I, J, L  (INTEGER) : GEOS-CHEM lon, lat, alt grid box indices
!  (4  ) N        (INTEGER) : Index for Hg(II) total or tagged tracers
!  (5  ) RKT      (REAL*8 ) : Value of R * k * Temp for drydep [unitless]
!  (6  ) E_RKT    (REAL*8 ) : Value of EXP( - RKT ) [unitless]
!  (7  ) LOST_Hg0 (REAL*8 ) : Amount of Hg(0) that became Hg(II) [kg]
!  (8  ) DTCHEM   (INTEGER) : Chemistry timestep [s]
!
!  NOTES:
!******************************************************************************
!
      ! References to F90 modules
      USE DIAG_MOD,     ONLY : AD44
      USE GRID_MOD,     ONLY : GET_AREA_CM2
      USE TRACER_MOD,   ONLY : STT
      USE TRACERID_MOD, ONLY : IDTHG2

#     include "CMN_SIZE"    ! Size parameters
#     include "CMN_DIAG"    ! ND44
#     include "CMN_O3"      ! XNUMOL

      ! Arguments
      INTEGER, INTENT(IN)  :: I,   J,     L,        N
      REAL*8,  INTENT(IN)  :: RKT, E_RKT, LOST_Hg0, DTCHEM

      ! Local variables
      REAL*8               :: AREA_CM2, DRYDEP, OLD_Hg2, NEW_Hg2

      !=================================================================
      ! RXN_Hg2_DRYD begins here!
      !=================================================================

      ! Error check tracer number
      IF ( N < 1 ) RETURN

      ! Initial concentration of Hg(II) [kg]
      OLD_Hg2      = MAX( STT(I,J,L,N), SMALLNUM )

      ! Concentration of Hg(II) after drydep: also accounts for
      ! the amount of Hg(0) which was converted to Hg(II) [kg]
      NEW_Hg2      = ( OLD_Hg2  *       E_RKT         )  + 
     &               ( LOST_Hg0/RKT * ( 1d0 - E_RKT ) )
 
      ! Save back into STT array [kg]
      NEW_Hg2      = MAX( NEW_Hg2, SMALLNUM )
      STT(I,J,L,N) = NEW_Hg2

      !=================================================================
      ! ND44 diagnostic: drydep flux of Hg(II) [molec/cm2/s]
      !=================================================================
      IF ( ND44 > 0 ) THEN
      
         ! Grid box surface area [cm2]
         AREA_CM2 = GET_AREA_CM2( J )
      
         ! Amt of Hg(II) lost to drydep [molec/cm2/s]
         DRYDEP   = ( OLD_Hg2 - NEW_Hg2 ) + LOST_Hg0
         DRYDEP   = DRYDEP * XNUMOL(IDTHG2) / ( AREA_CM2 * DTCHEM )
      
         ! Archive Hg(II) drydep flux in AD44 array [molec/cm2/s]
         AD44(I,J,N,1) = AD44(I,J,N,1) + DRYDEP
      ENDIF
          
      ! Return to calling program
      END SUBROUTINE RXN_HG2_DRYD

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
!  mercury.  HgP is lost via dry deposition. (eck, bmy, 12/7/04)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) V_DEP_HgP (REAL*8) : Dry deposition velocity for Hg(II) [cm/s]
!     
!  NOTES:
!******************************************************************************
!
      ! Refernces to F90 modules
      USE DIAG_MOD,     ONLY : AD44 
      USE GRID_MOD,     ONLY : GET_AREA_CM2
      USE TIME_MOD,     ONLY : GET_TS_CHEM
      USE LOGICAL_MOD,  ONLY : LSPLIT
      USE TRACER_MOD,   ONLY : STT
      USE TRACERID_MOD, ONLY : IDTHGP,    IDTHGP_NA,
     &                         IDTHGP_EU, IDTHGP_AS, IDTHGP_RW

#     include "CMN_SIZE"  ! Size parameters
#     include "CMN_O3"    ! XNUMOL
#     include "CMN_DIAG"  ! ND03, LD03

      ! Arguments
      REAL*8, INTENT(IN) :: V_DEP_HgP(IIPAR,JJPAR)

      ! Local variables
      INTEGER            :: I,       J
      REAL*8             :: DTCHEM, E_KT

      !=================================================================
      ! CHEM_HgP begins here!
      !=================================================================

      ! Chemistry timestep [s]
      DTCHEM = GET_TS_CHEM() * 60d0

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, E_KT )
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! Pre-compute the exponential term for use below [unitless]
         E_KT = EXP( -V_DEP_HgP(I,J) * DTCHEM )

         !--------------------
         ! Total HgP
         !--------------------

         ! Compute new concentration of HgP after drydep [kg]
         CALL RXN_HgP_DRYD( I, J, 1, IDTHgP, E_KT, DTCHEM )
         
         !--------------------
         ! Tagged HgP tracers
         !--------------------
         IF ( LSPLIT ) THEN

            ! Compute new conc of tagged HgP tracers after drydep [kg]
            CALL RXN_HgP_DRYD( I, J, 1, AN_HgP(I,J), E_KT, DTCHEM )

!------------------------------------------------------------------------------
! NOTE: Leave for future expansion for other HgP tagged tracers
! (eck, bmy, 12/14/04)
!            CALL RXN_HgP_DRYD( I, J, 1, OC_HgP, E_KT, DTCHEM )
!            CALL RXN_HgP_DRYD( I, J, 1, LN_HgP, E_KT, DTCHEM )
!            CALL RXN_HgP_DRYD( I, J, 1, NT_HgP, E_KT, DTCHEM )
!------------------------------------------------------------------------------
            
         ENDIF

      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      ! Return to calling program
      END SUBROUTINE CHEM_HGP

!------------------------------------------------------------------------------

      SUBROUTINE RXN_HgP_DRYD( I, J, L, N, E_KT, DTCHEM )
!
!******************************************************************************
!  Subroutine RXN_Hg0_DRYD computes the new concentration of HgP after
!  dry deposition.  ND44 diagnostics are also archived. (eck, bmy, 12/14/04)
! 
!  Arguments as Input:
!  ============================================================================
!  (1-3) I, J, L  (INTEGER) : GEOS-CHEM lon, lat, alt grid box indices
!  (4  ) N        (INTEGER) : Index for HgP total or tagged tracers
!  (5  ) E_KT     (REAL*8 ) : Value of EXP( -KT ) for drydep rxn [unitless]
!  (6  ) DTCHEM   (INTEGER) : Chemistry timestep [s]
!
!  NOTES:
!******************************************************************************
!
      ! References to F90 modules
      USE DIAG_MOD,     ONLY : AD44
      USE GRID_MOD,     ONLY : GET_AREA_CM2
      USE TRACER_MOD,   ONLY : STT
      USE TRACERID_MOD, ONLY : IDTHGP
    
#     include "CMN_SIZE"   ! Size parameters
#     include "CMN_DIAG"   ! ND44
#     include "CMN_O3"     ! XNUMOL

      ! Arguments
      INTEGER, INTENT(IN)  :: I, J, L, N
      REAL*8,  INTENT(IN)  :: DTCHEM, E_KT

      ! Local variables
      REAL*8               :: AREA_CM2, DRYDEP, OLD_HgP, NEW_HgP

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
         AREA_CM2 = GET_AREA_CM2( J )

         ! Amt of Hg(II) lost to drydep [molec/cm2/s]
         DRYDEP   = OLD_HgP - NEW_HgP  
         DRYDEP   = DRYDEP * XNUMOL(IDTHGP) / ( AREA_CM2 * DTCHEM )

         ! Archive Hg(II) drydep flux in AD44 array [molec/cm2/s]
         AD44(I,J,N,1) = AD44(I,J,N,1) + DRYDEP

      ENDIF

      ! Return to calling program
      END SUBROUTINE RXN_HgP_DRYD

!------------------------------------------------------------------------------

      SUBROUTINE EMISSMERCURY
!
!******************************************************************************
!  Subroutine EMISSMERCURY is the driver routine for mercury emissions.
!  (eck, bmy, 12/7/04)
! 
!  NOTES:
!******************************************************************************
!
      ! References to F90 modules
      USE ERROR_MOD,    ONLY : DEBUG_MSG
      USE LOGICAL_MOD,  ONLY : LPRT
      USE TRACER_MOD,   ONLY : STT
      USE TRACERID_MOD, ONLY : IDTHg0, IDTHg2, IDTHgP
      
#     include "CMN_SIZE"

      ! Local variables
      LOGICAL, SAVE      :: FIRST = .TRUE. 
     
      !=================================================================
      ! EMISSMERCURY begins here!
      !=================================================================

      ! First-time initialization
      IF ( FIRST ) THEN

         ! Allocate arrays (if not done before)
         CALL INIT_MERCURY

         ! Read anthro, ocean, land emissions of Hg from disk
         CALL MERCURY_READYR

         ! Reset first-time flag
         FIRST = .FALSE.
      ENDIF

      !=================================================================
      ! Call emission routines for Hg(0), Hg(II), and Hg(P)
      !=================================================================

      ! Compute fraction of each layer in the PBL
      CALL COMPUTE_FEMIS

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
!  Subroutine SRCHg0 is the subroutine for Hg(0) emissions.  Emissions of
!  Hg(0) will be distributed throughout the boundary layer. (eck, bmy, 12/7/04)
! 
!  Arguments as Input/Output:
!  ============================================================================
!  (1 ) TC (REAL*8) : Tracer concentration of Hg(0) [kg]
!
!  NOTES:
!******************************************************************************
!
      ! Reference to diagnostic arrays
      USE DIAG_MOD,     ONLY : AD03_Hg0_an, AD03_Hg0_oc, 
     &                         AD03_Hg0_ln, AD03_Hg0_nt
      USE ERROR_MOD,    ONLY : ERROR_STOP
      USE LOGICAL_MOD,  ONLY : LSPLIT
      USE TIME_MOD,     ONLY : GET_TS_EMIS
      USE TRACER_MOD,   ONLY : STT
      USE TRACERID_MOD, ONLY : IDTHG0

#     include "CMN_SIZE"     ! Size parameters
#     include "CMN_DIAG"     ! ND03

      ! Local variables
      INTEGER               :: I,      J,    L,    N
      REAL*8                :: DTSRCE, E_Hg, T_Hg, T_Hg_An

      !=================================================================
      ! SRCHg0 begins here!
      !=================================================================

      ! Emission timestep [s]
      DTSRCE = GET_TS_EMIS() * 60d0

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, N, T_Hg_An, T_Hg, E_Hg )
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         !==============================================================
         ! Handle total Hg(0) emissions
         !==============================================================

         ! Compute total anthropogenic Hg(0) emissions
         T_Hg_An = EHg0_an(I,J,1) + EHg0_an(I,J,2) + EHg0_an(I,J,3)

         ! Compute total Hg(0) emissions (anthro+oceans+land+natural)
         T_Hg    = T_Hg_An        + EHg0_oc(I,J)   + 
     &             EHg0_ln(I,J)   + EHg0_nt(I,J)

         ! Partition total Hg(0) throughout PBL; store into STT [kg]
         DO L = 1, LLTROP
            E_Hg              = FEMIS(L,I,J)      * T_Hg
            STT(I,J,L,IDTHG0) = STT(I,J,L,IDTHG0) + ( E_Hg * DTSRCE )
         ENDDO

         !==============================================================
         ! Handle tagged Hg(0) emissions
         !==============================================================
         IF ( LSPLIT ) THEN

            ! Partition tagged Hg(0) throughout PBL; store into STT [kg]
            DO L = 1, LLTROP

               ! Anthro Hg0 by region
               N            = AN_Hg0(I,J)
               E_Hg         = FEMIS(L,I,J) * T_Hg_an
               STT(I,J,L,N) = STT(I,J,L,N) + ( E_Hg * DTSRCE )

               ! Oceanic Hg0 
               N            = OC_Hg0
               E_Hg         = FEMIS(L,I,J) * EHg0_oc(I,J)
               STT(I,J,L,N) = STT(I,J,L,N) + ( E_Hg * DTSRCE )

               ! Land re-emissions of Hg0
               N            = LN_Hg0
               E_Hg         = FEMIS(L,I,J) * EHg0_ln(I,J)
               STT(I,J,L,N) = STT(I,J,L,N) + ( E_Hg * DTSRCE ) 

               ! Natural land sources
               N            = NT_Hg0
               E_Hg         = FEMIS(L,I,J) * EHg0_nt(I,J)
               STT(I,J,L,N) = STT(I,J,L,N) + ( E_Hg * DTSRCE )

            ENDDO
         ENDIF
        
         !==============================================================
         ! ND03 diagnostic: Total Hg(0) emissions [kg]
         !==============================================================
         IF ( ND03 > 0 ) THEN

            ! Anthro emissions
            AD03_Hg0_an(I,J) = AD03_Hg0_an(I,J) 
     &                       + ( T_Hg_An      * DTSRCE )

            ! Ocean emissions
            AD03_Hg0_oc(I,J) = AD03_Hg0_oc(I,J) 
     &                       + ( EHg0_oc(I,J) * DTSRCE )

            ! Re-emission from land
            AD03_Hg0_ln(I,J) = AD03_Hg0_ln(I,J) 
     &                       + ( EHg0_ln(I,J) * DTSRCE )

            ! Natural land source emissions
            AD03_Hg0_nt(I,J) = AD03_Hg0_nt(I,J) 
     &                       + ( EHg0_nt(I,J) * DTSRCE )
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
!  (eck, bmy, 12/7/04)
! 
!  Arguments as Input/Output:
!  ============================================================================
!  (1 ) TC (REAL*8) : Tracer concentration of Hg(II) [kg]
!
!  NOTES:
!******************************************************************************
!
      ! Reference to F90 modules
      USE DIAG_MOD,     ONLY : AD03_HG2_AN
      USE ERROR_MOD,    ONLY : ERROR_STOP
      USE LOGICAL_MOD,  ONLY : LSPLIT
      USE TIME_MOD,     ONLY : GET_TS_EMIS
      USE TRACER_MOD,   ONLY : STT
      USE TRACERID_MOD, ONLY : IDTHg2

#     include "CMN_SIZE"      ! Size parameters
#     include "CMN_DIAG"      ! ND03 (for now)

      ! Local variables
      INTEGER                :: I,      J,    L,   N
      REAL*8                 :: DTSRCE, T_Hg, E_Hg

      !=================================================================
      ! SRCHg2 begins here!
      !=================================================================
      
      ! Emission timestep [s]
      DTSRCE = GET_TS_EMIS() * 60d0

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, N, T_Hg, E_Hg )
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         !==============================================================
         ! Handle total Hg(II) emissions
         !==============================================================

         ! Compute total Hg(II) emission over all layers
         T_Hg  = EHg2_an(I,J,1) + EHg2_an(I,J,2) + EHg2_an(I,J,3)

         ! Partition total Hg(II) throughout PBL; store in STT [kg]
         DO L  = 1, LLTROP
            E_Hg              = FEMIS(L,I,J)      * T_Hg
            STT(I,J,L,IDTHg2) = STT(I,J,L,IDTHg2) + ( E_Hg * DTSRCE )
         ENDDO

         !==============================================================
         ! Handle tagged Hg(II) emissions
         !==============================================================
         IF ( LSPLIT ) THEN 

            ! Partition tagged Hg(II) throughout PBL; store in STT [kg]
            DO L = 1, LLTROP
               N            = AN_Hg2(I,J)
               E_Hg         = FEMIS(L,I,J) * T_Hg
               STT(I,J,L,N) = STT(I,J,L,N) + ( E_Hg * DTSRCE ) 
            ENDDO

         ENDIF

         !==============================================================
         ! ND03 diagnostic: Total anthro Hg(II) emissions [kg]
         !==============================================================
         IF ( ND03 > 0 ) THEN
            AD03_Hg2_an(I,J) = AD03_Hg2_an(I,J) + ( T_Hg * DTSRCE )
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
!  (eck, bmy, 12/7/04)
! 
!  Arguments as Input/Output:
!  ============================================================================
!  (1 ) TC (REAL*8) : Tracer concentration of Hg(II) [kg]
!
!  NOTES:
!******************************************************************************
!
      ! Reference to diagnostic arrays
      USE DIAG_MOD,     ONLY : AD03_HGP_AN
      USE ERROR_MOD,    ONLY : ERROR_STOP
      USE LOGICAL_MOD,  ONLY : LSPLIT
      USE TIME_MOD,     ONLY : GET_TS_EMIS 
      USE TRACER_MOD,   ONLY : STT
      USE TRACERID_MOD, ONLY : IDTHGP

#     include "CMN_SIZE"      ! Size paramters
#     include "CMN_DIAG"      ! ND03 (for now)

      ! Local variables
      INTEGER                :: I,      J,    L,   N
      REAL*8                 :: DTSRCE, T_Hg, E_Hg 

      !=================================================================
      ! SRCHgP begins here!
      !=================================================================

      ! Chemistry timestep [s]
      DTSRCE = GET_TS_EMIS() * 60d0

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, N, T_Hg, E_Hg )
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         !==============================================================
         ! Handle total HgP emissions
         !==============================================================

         ! Compute total HgP emissions over all levels
         T_Hg  = EHgP_an(I,J,1) + EHgP_an(I,J,2) + EHgP_an(I,J,3) 

         ! Partition total HgP throughout PBL; store into STT [kg]
         DO L  = 1, LLTROP
            E_Hg              = FEMIS(L,I,J)      * T_Hg
            STT(I,J,L,IDTHGP) = STT(I,J,L,IDTHGP) + ( E_Hg * DTSRCE )
         ENDDO

         !==============================================================
         ! Handle total HgP emissions
         !==============================================================
         IF ( LSPLIT ) THEN

            ! Partition total HgP throughout PBL; store into STT [kg]
            DO L  = 1, LLTROP
               N            = AN_HgP(I,J)
               E_Hg         = FEMIS(L,I,J) * T_Hg
               STT(I,J,L,N) = STT(I,J,L,N) + ( E_Hg * DTSRCE )
            ENDDO
            
         ENDIF

         !==============================================================
         ! ND03 diagnostic: Total anthro HgP emissions [kg]
         !==============================================================
         IF ( ND03 > 0 ) THEN
            AD03_HgP_an(I,J) = AD03_HgP_an(I,J) + ( T_Hg * DTSRCE )
         ENDIF

      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      ! Return to calling program
      END SUBROUTINE SRCHgP

!------------------------------------------------------------------------------

      SUBROUTINE COMPUTE_FEMIS
!
!******************************************************************************
!  Subroutine COMPUTE_FEMIS computes the fraction of the entire boundary
!  layer which is occupied by each vertical level. (eck, bmy, 12/7/04)
!
!  NOTES:
!******************************************************************************
!
      ! Reference to diagnostic arrays
      USE DAO_MOD,      ONLY : PBL
      USE ERROR_MOD,    ONLY : ERROR_STOP
      USE PRESSURE_MOD, ONLY : GET_PEDGE

#     include "CMN_SIZE"     ! Size parameters
#     include "CMN_GCTM"     ! Physical constants

      ! Local variables 
      INTEGER               :: I,     J,      L
      REAL*8                :: BLTOP, BLTHIK, DELP, FTOT, PB, PT

      !=================================================================
      ! COMPUTE_FEMIS begins here
      !=================================================================
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, FTOT, BLTOP, BLTHIK, PB, PT, DELP )
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! Zero FTOT summing variable
         FTOT = 0d0

         ! Zero FEMIS array
         DO L = 1, LLTROP
            FEMIS(L,I,J) = 0d0
         ENDDO

#if   defined( GEOS_4 )

         ! BLTOP = pressure at PBL top [hPa]
         ! Use barometric law since PBL is in [m]
         BLTOP  = GET_PEDGE(I,J,1) * EXP( -PBL(I,J) / SCALE_HEIGHT )

         ! BLTHIK is PBL thickness [hPa]
         BLTHIK = GET_PEDGE(I,J,1) - BLTOP

#else

         ! BLTOP = pressure of PBL top [hPa]
         BLTOP  = GET_PEDGE(I,J,1) - MAX( PBL(I,J), 1.D0 )

         ! BLTHIK is PBL thickness [hPa]
         BLTHIK = MAX( PBL(I,J), 1.D0 )

#endif

         !==============================================================
         ! Loop thru tropospheric levels and compute the fraction
         ! of the PBL occupied by each level.  Store in FEMIS array.
         !==============================================================
         DO L  = 1, LLTROP

            ! Pressure at bottom & top edges of grid box (I,J,L) [hPa]
            PB   = GET_PEDGE(I,J,L)
            PT   = GET_PEDGE(I,J,L+1)

            ! Thickness of grid box (I,J,L) [hPa]
            DELP = PB - PT

            ! FEMIS is the fraction of the PBL 
            ! which is occupied by this level L
            IF ( BLTOP <= PT )  THEN
               FEMIS(L,I,J) = DELP / BLTHIK

            ELSEIF ( BLTOP > PT .AND. BLTOP <= PB ) THEN
               FEMIS(L,I,J) = ( PB - BLTOP ) / BLTHIK

            ELSEIF ( BLTOP > PB ) THEN
               CYCLE      
 
            ENDIF
            
            ! Fraction of data partitioned into 
            ! each level should sum to 1.0
            FTOT = FTOT + FEMIS(L,I,J)

         ENDDO

         ! Error check
         IF ( ABS( FTOT - 1.d0 ) > 1.d-3 ) THEN
            CALL ERROR_STOP( 'Check vertical. distribution!',
     &                       'COMPUTE_FEMIS ("mercury_mod.f")' )
         ENDIF

      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      ! Return to calling program
      END SUBROUTINE COMPUTE_FEMIS

!------------------------------------------------------------------------------

      SUBROUTINE MERCURY_READYR
!
!******************************************************************************
!  Subroutine MERCURY_READYR reads the year-invariant emissions for Mercury
!  from anthropogenic, ocean, and land sources. (eck, bmy, 12/6/04)
!  
!  NOTES:
!******************************************************************************
!
      ! References to F90 modules
      USE BPCH2_MOD
      USE DIRECTORY_MOD, ONLY : DATA_DIR
      USE TRANSFER_MOD,  ONLY : TRANSFER_2D

#     include "CMN_SIZE"   ! Size parameters

      ! Local variables
      INTEGER             :: L
      REAL*4              :: ARRAY(IGLOB,JGLOB,Hg_LEVELS)
      REAL*8              :: XTAU
      REAL*8, PARAMETER   :: SEC_PER_YR = 365.25d0 * 86400d0
      CHARACTER(LEN=255)  :: FILENAME

      !=================================================================
      ! MERCURY_READYR begins here!
      !
      ! Read annual anthropogenic mercury emissions [kg/s]
      !=================================================================

      ! Filename for anthropogenic mercury source
      FILENAME = TRIM( DATA_DIR )                  // 
     &           'mercury_200412/Hg_anthro.geos.'  // GET_RES_EXT()

      ! Echo info
      WRITE( 6, 100 ) TRIM( FILENAME )
 100  FORMAT( 'MERCURY_READYR: Reading ', a )

      ! TAU value corresponding to the data
      XTAU = GET_TAU0( 1, 1, 1995 )

      !---------------------------
      ! Hg(0) emissions [kg/s]
      !---------------------------

      ! Read data in [kg/yr]
      CALL READ_BPCH2( FILENAME, 'HG-SRCE', 1, 
     &                 XTAU,      IGLOB,    JGLOB,    
     &                 3,         ARRAY,    QUIET=.TRUE. )

      ! Cast to REAL*8 and convert from [kg/yr] to [kg/s]
      DO L = 1, Hg_LEVELS
         CALL TRANSFER_2D( ARRAY(:,:,L), EHg0_an(:,:,L) )
         EHg0_an(:,:,L) = EHg0_an(:,:,L) / SEC_PER_YR
      ENDDO

      !---------------------------
      ! Hg(II) emissions [kg/s]
      !---------------------------

      ! Read data in [kg/yr]
      CALL READ_BPCH2( FILENAME, 'HG-SRCE', 5, 
     &                 XTAU,      IGLOB,    JGLOB,    
     &                 3,         ARRAY,    QUIET=.TRUE. )

      ! Cast to REAL*8 and convert from [kg/yr] to [kg/s]
      DO L = 1, Hg_LEVELS
         CALL TRANSFER_2D( ARRAY(:,:,L), EHg2_an(:,:,L) )
         EHg2_an(:,:,L) = EHg2_an(:,:,L) / SEC_PER_YR
      ENDDO

      !---------------------------
      ! HgP emissions [kg/s]
      !---------------------------

      ! Read data in [kg/yr]
      CALL READ_BPCH2( FILENAME, 'HG-SRCE', 6, 
     &                 XTAU,      IGLOB,    JGLOB,    
     &                 3,         ARRAY,    QUIET=.TRUE. )

      ! Cast to REAL*8 and convert from [kg/yr] to [kg/s]
      DO L = 1, Hg_LEVELS
         CALL TRANSFER_2D( ARRAY(:,:,L), EHgP_an(:,:,L) )
         EHgP_an(:,:,L) = EHgP_an(:,:,L) / SEC_PER_YR
      ENDDO
  
      !=================================================================
      ! Read annual emissions of oceanic Hg(0) [kg/s]
      !=================================================================

      ! Filename for ocean source
      FILENAME = TRIM( DATA_DIR )                // 
     &           'mercury_200412/Hg_ocean.geos.' // GET_RES_EXT()    

      ! Echo info
      WRITE( 6, 100 ) TRIM( FILENAME )

      ! TAU value corresponding to the data 
      XTAU  = GET_TAU0( 1, 1, 1990 )

      ! Read data in [kg/yr]
      CALL READ_BPCH2( FILENAME, 'HG-SRCE',     2,  
     &                 XTAU,      IGLOB,        JGLOB,    
     &                 1,         ARRAY(:,:,1), QUIET=.TRUE. )

      ! Cast to REAL*8 and resize
      CALL TRANSFER_2D( ARRAY(:,:,1), EHg0_oc )

      ! Convert from [kg/yr] to [kg/s]
      EHg0_oc = EHg0_oc / SEC_PER_YR
 
      !=================================================================
      ! Read annual emissions of anthropogenic Hg(0) which is
      ! re-emitted from the land [kg/s]
      !=================================================================

      ! Filename for re-emitted anthropogenic mercury
      FILENAME = TRIM( DATA_DIR )               // 
     &           'mercury_200412/Hg_land.geos.' // GET_RES_EXT()   

      ! Echo info
      WRITE( 6, 100 ) TRIM( FILENAME )

      ! TAU value corresponding to the data 
      XTAU  = GET_TAU0( 1, 1, 1990 )

      ! Read data in [kg/yr]
      CALL READ_BPCH2( FILENAME, 'HG-SRCE',     3,  
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

      ! Filename for natural land-source mercury
      FILENAME = TRIM( DATA_DIR )                  // 
     &           'mercury_200412/Hg_natural.geos.' // GET_RES_EXT() 
      ! Echo info
      WRITE( 6, 100 ) TRIM( FILENAME )

      ! TAU value corresponding to the data 
      XTAU   = GET_TAU0( 1, 1, 1995 )

      ! Read data in [kg/yr]
      CALL READ_BPCH2( FILENAME, 'HG-SRCE',     4, 
     &                 XTAU,      IGLOB,        JGLOB,    
     &                 1,         ARRAY(:,:,1), QUIET=.TRUE. )

      ! Cast to REAL*8 and resize
      CALL TRANSFER_2D( ARRAY(:,:,1), EHg0_nt )

      ! Convert from [kg/yr] to [kg/s]
      EHg0_nt = EHg0_nt / SEC_PER_YR

      !=================================================================
      ! Print totals to the screen in [Gg/yr]
      !=================================================================
      WRITE( 6, '(a)' ) REPEAT( '=', 79 )
      WRITE( 6, 110   )
      WRITE( 6, '(a)' )
      WRITE( 6, 111   ) SUM( EHg0_an ) * SEC_PER_YR * 1d-6
      WRITE( 6, 112   ) SUM( EHg0_oc ) * SEC_PER_YR * 1d-6
      WRITE( 6, 113   ) SUM( EHg0_ln ) * SEC_PER_YR * 1d-6
      WRITE( 6, 114   ) SUM( EHg0_nt ) * SEC_PER_YR * 1d-6
      WRITE( 6, 115   ) SUM( EHg2_an ) * SEC_PER_YR * 1d-6
      WRITE( 6, 116   ) SUM( EHgP_an ) * SEC_PER_YR * 1d-6
      WRITE( 6, '(a)' ) REPEAT( '=', 79 )

      ! FORMAT strings
 110  FORMAT( 'M E R C U R Y   E M I S S I O N S' )
 111  FORMAT( 'Total Anthro     Hg(0)  : ', f7.3, ' [Gg/yr]' )
 112  FORMAT( 'Total Oceanic    Hg(0)  : ', f7.3, ' [Gg/yr]' )
 113  FORMAT( 'Total Re-Emitted Hg(0)  : ', f7.3, ' [Gg/yr]' )
 114  FORMAT( 'Total Natural    Hg(0)  : ', f7.3, ' [Gg/yr]' )
 115  FORMAT( 'Total Anthro     Hg(II) : ', f7.3, ' [Gg/yr]' )
 116  FORMAT( 'Total Anthro     HgP    : ', f7.3, ' [Gg/yr]' )

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
      USE TIME_MOD, ONLY : GET_NHMSb,   GET_ELAPSED_SEC, 
     &                     GET_TS_CHEM, GET_DAY_OF_YEAR, GET_GMT

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
!  anthropogenic (by geographic region), oceanic, land re-emission, and natural
!  Hg0, Hg2, and HgP.  The position of Hg2 and HgP in the DEPSAV array is
!  also computed for future use.  This routine only has to be called once
!  at the start of the simulation. (eck, bmy, 12/15/04)
!
!  NOTES:
!******************************************************************************
!
      ! References to F90 modules
      USE GRID_MOD,   ONLY : GET_XMID, GET_YMID
      USE ERROR_MOD,  ONLY : ERROR_STOP
      USE TRACER_MOD, ONLY : N_TRACERS
      USE TRACERID_MOD

#     include "CMN_SIZE"   ! Size parameters

      ! Local variables
      INTEGER             :: I, J
      REAL*8              :: X, Y
      CHARACTER(LEN=255)  :: LOCATION

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
               AN_Hg0(I,J) = IDTHg0_NA
               AN_Hg2(I,J) = IDTHg2_NA
               AN_HgP(I,J) = IDTHgP_NA
            
            ! European Anthro Hg (1st sub-box)
            ELSE IF ( ( X >= -17.5 .and. X < 72.5 )  .and. 
     &                ( Y >=  36.0 .and. Y < 45.0 ) ) THEN
               AN_Hg0(I,J) = IDTHg0_EU
               AN_Hg2(I,J) = IDTHg2_EU
               AN_HgP(I,J) = IDTHgP_EU

            ! European Anthro Hg (2nd sub-box)
            ELSE IF ( ( X >= -17.5 .and. X < 172.5 )  .and. 
     &                ( Y >=  45.0 .and. Y <  88.0 ) ) THEN
               AN_Hg0(I,J) = IDTHg0_EU 
               AN_Hg2(I,J) = IDTHg2_EU
               AN_HgP(I,J) = IDTHgP_EU

            ! Asian Anthro Hg 
            ELSE IF ( ( X >= 70.0 .and. X < 152.5 )  .and.
     &                ( Y >=  8.0 .and. Y <  45.0 ) ) THEN
               AN_Hg0(I,J) = IDTHg0_AS
               AN_Hg2(I,J) = IDTHg2_AS
               AN_HgP(I,J) = IDTHgP_AS

            ! Rest-of-world Anthro Hg 
            ELSE
               AN_Hg0(I,J) = IDTHg0_RW
               AN_Hg2(I,J) = IDTHg2_RW
               AN_HgP(I,J) = IDTHgP_RW

            ENDIF
         ENDDO
      ENDDO
!$OMP END PARALLEL DO

      !---------------------------------
      ! Oceanic tracer indices
      !---------------------------------
      OC_Hg0 = IDTHg0_OC
      OC_Hg2 = IDTHg2_OC
      OC_HgP = IDTHgP_OC

      !---------------------------------
      ! Land re-emission tracer indices
      !---------------------------------
      LN_Hg0 = IDTHg0_LN
      LN_Hg2 = IDTHg2_LN
      LN_HgP = IDTHgP_LN

      !---------------------------------
      ! Natural emission tracer indices
      !---------------------------------
      NT_Hg0 = IDTHg0_NT
      NT_Hg2 = IDTHg2_NT
      NT_HgP = IDTHgP_NT
    
      !=================================================================
      ! Error check tracers: make sure they are not zero
      ! since this can cause array-out-of-bounds errors
      !=================================================================

      !-------------
      ! Tagged Hg0
      !-------------
      IF ( ANY( AN_Hg0 == 0 ) ) THEN
         CALL ERROR_STOP( 'Hg0_an_xx tracers are undefined!', LOCATION )   
      ENDIF
      
      IF ( OC_Hg0 == 0 ) THEN
         CALL ERROR_STOP( 'Hg0_oc tracer is undefined!', LOCATION )
      ENDIF

      IF ( LN_Hg0 == 0 ) THEN
         CALL ERROR_STOP( 'Hg0_ln tracer is undefined!', LOCATION )
      ENDIF

      IF ( NT_Hg0 == 0 ) THEN
         CALL ERROR_STOP( 'Hg0_nt tracer is undefined!', LOCATION )
      ENDIF

      !-------------
      ! Tagged Hg2
      !-------------
      IF ( ANY( AN_Hg2 == 0 ) ) THEN
         CALL ERROR_STOP( 'Hg2_an_xx tracers are undefined!', LOCATION )   
      ENDIF

      IF ( OC_Hg2 == 0 ) THEN
         CALL ERROR_STOP( 'Hg2_oc tracer is undefined!', LOCATION )
      ENDIF

      IF ( LN_Hg2 == 0 ) THEN
         CALL ERROR_STOP( 'Hg2_ln tracer is undefined!', LOCATION )
      ENDIF

      IF ( NT_Hg2 == 0 ) THEN
         CALL ERROR_STOP( 'Hg2_nt tracer is undefined!', LOCATION )
      ENDIF

      !-------------
      ! Tagged Hg2
      !-------------
      IF ( ANY( AN_HgP == 0 ) ) THEN
         CALL ERROR_STOP( 'HgP_an_xx tracer are undefined!', LOCATION )   
      ENDIF

      !------------------------------------------------------------------
      ! NOTE: These tracers have not yet been implemented but we
      ! will leave this here for future use
      !IF ( OC_HgP == 0 ) THEN
      !   CALL ERROR_STOP( 'HgP tracer is undefined!', LOCATION )
      !ENDIF
      !
      !IF ( LN_HgP == 0 ) THEN
      !   CALL ERROR_STOP( 'HgP_ln tracer is undefined!', LOCATION )
      !ENDIF
      !
      !IF ( NT_HgP == 0 ) THEN
      !   CALL ERROR_STOP( 'HgP_nt tracer is undefined!', LOCATION )
      !ENDIF
      !------------------------------------------------------------------

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

      SUBROUTINE INIT_MERCURY
!
!******************************************************************************
!  Subroutine INIT_MERCURY allocates all module arrays (eck, bmy, 12/2/04)
!  
!  NOTES:
!******************************************************************************
!
      ! References to F90 modules
      USE DRYDEP_MOD,   ONLY : DEPNAME,   NUMDEP
      USE ERROR_MOD,    ONLY : ALLOC_ERR, ERROR_STOP
      USE LOGICAL_MOD,  ONLY : LSPLIT,    LDRYD
      USE TRACERID_MOD, ONLY : IDTHg0,    IDTHg2,   IDTHgP

#     include "CMN_SIZE"     ! Size parameters

      ! Local variables
      LOGICAL, SAVE         :: IS_INIT = .FALSE. 
      INTEGER               :: AS, N
      CHARACTER(LEN=255)    :: LOCATION

      !=================================================================
      ! INIT_MERCURY begins here!
      !=================================================================

      ! Return if we have already allocated arrays
      IF ( IS_INIT ) RETURN

      ! Location string for ERROR_STOP
      LOCATION = 'DEFINE_TAGGED_Hg ("mercury_mod.f")'

      ! Error check Hg0 tracer
      IF ( IDTHg0 == 0 ) THEN
         CALL ERROR_STOP( 'Total Hg0 tracer is undefined!', LOCATION )
      ENDIF

      ! Error check Hg2 tracer
      IF ( IDTHg2 == 0 ) THEN
         CALL ERROR_STOP( 'Total Hg2 tracer is undefined!', LOCATION )
      ENDIF

      ! Error check HgP tracer
      IF ( IDTHgP == 0 ) THEN
         CALL ERROR_STOP( 'Total HgP tracer is undefined!', LOCATION )
      ENDIF

      !=================================================================
      ! Allocate arrays
      !=================================================================
      ALLOCATE( COSZM( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'COSZM' )
      COSZM = 0d0

      ALLOCATE( EHg0_an( IIPAR, JJPAR, 3 ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'EHg0_an' )
      EHg0_an = 0d0

      ALLOCATE( EHg2_an( IIPAR, JJPAR, 3 ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'EHg2_an' )
      EHg2_an = 0d0

      ALLOCATE( EHgP_an( IIPAR, JJPAR, 3 ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'EHgP_an' )
      EHgP_an = 0d0

      ALLOCATE( EHg0_oc( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'EHg0_oc' )
      EHg0_oc = 0d0

      ALLOCATE( EHg0_ln( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'EHg0_ln' )
      EHg0_ln = 0d0

      ALLOCATE( EHg0_nt( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'EHg0_nt' )
      EHg0_nt = 0d0

      ALLOCATE( FEMIS( LLTROP, IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'FEMIS' )
      FEMIS = 0d0

      ALLOCATE( TCOSZ( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'TCOSZ' )
      TCOSZ = 0d0

      ALLOCATE( TTDAY( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'TTDAY' )
      TTDAY = 0d0

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

         ! Tracer indices for land re-emission
         LN_Hg0 = 0d0
         LN_Hg2 = 0d0
         LN_HgP = 0d0

         ! Tracer indices for natural land emissions
         NT_Hg0 = 0d0
         NT_Hg2 = 0d0
         NT_HgP = 0d0

         ! Tracer indices for oceanic emissions
         OC_Hg0 = 0d0
         OC_Hg2 = 0d0
         OC_HgP = 0d0

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
      ENDIF

      ! Reset IS_INIT, since we have already allocated arrays
      IS_INIT = .TRUE.

      ! Return to calling program
      END SUBROUTINE INIT_MERCURY

!------------------------------------------------------------------------------

      SUBROUTINE CLEANUP_MERCURY
!
!******************************************************************************
!  Subroutine CLEANUP_MERCURY deallocates all module arrays (eck, bmy, 12/6/04)
!
!  NOTES:
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
      IF ( ALLOCATED( FEMIS   ) ) DEALLOCATE( FEMIS   )
      IF ( ALLOCATED( TCOSZ   ) ) DEALLOCATE( TCOSZ   )
      IF ( ALLOCATED( TTDAY   ) ) DEALLOCATE( TTDAY   )

      ! Return to calling program
      END SUBROUTINE CLEANUP_MERCURY

!------------------------------------------------------------------------------

      END MODULE MERCURY_MOD
