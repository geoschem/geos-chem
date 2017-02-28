! $Id: acetone_mod.f,v 1.8 2006/09/08 19:20:50 bmy Exp $
      MODULE ALD2_MOD
!
!******************************************************************************
!  F90 module ALD2_MOD contains some of the acetaldehyde subroutines.
!  Note that some sources are handled directly in the appropriate locations
!  e.g., biogenic emissions in megan_mod, biomass burning emissions, etc.
!  dbm, 08/11/2008

      !=================================================================
      ! MODULE PRIVATE DECLARATIONS -- keep certain internal variables 
      ! and routines from being seen outside "acetone_mod.f"
      !=================================================================

      ! PRIVATE module variables
  
      !=================================================================
      ! MODULE VARIABLES
      !=================================================================

      ! Molecules C / kg C 
      REAL*8,  PARAMETER   :: XNUMOL_C = 6.022d+23 / 12d-3 

      !=================================================================
      ! MODULE ROUTINES -- follow below the "CONTAINS" statement 
      !=================================================================
      CONTAINS
!------------------------------------------------------------------------------

      SUBROUTINE EMISS_ALD2_SEN( I,    J,    ACETALDEHYDE )
!
!******************************************************************************
!  Subroutine EMISS_ALD2_SEN computes the biogenic emissions of ACETALDEHYDE
!  from decaying plants (dbm, 08/11/2008)
!
!  Arguments as Input:
!  ============================================================================
!  (1-2) I, J    (INTEGER) : Lon/lat grid box indices
!
!  Arguments as Input/Output:
!  ============================================================================
!  (2  ) ACETALDEHYDE (REAL*8 ) : Total biogenic acetaldehyde emissions from
!        decaying plants [atoms C/box/s]
!
!******************************************************************************
!
      ! References to F90 modules
      USE GRID_MOD, ONLY : GET_AREA_M2
      USE ACETONE_MOD, ONLY: XRESP!, XNUMOL_C

      USE CMN_SIZE_MOD

      ! Arguments
      INTEGER, INTENT(IN)    :: I, J 
      REAL*8,  INTENT(OUT)   :: ACETALDEHYDE

      ! Local variables
      REAL*8                 :: YIELD_RESP,  ALD2_RESP,    AREA_M2

      !=================================================================
      ! EMISS_ALD2_SEN begins here!
      !=================================================================

      ! Grid box areas in [m2] and [cm2]
      AREA_M2  = GET_AREA_M2(I, J,1 )

      !=================================================================
      ! (1) BIOGENIC ACETALDEHYDE FROM DRY LEAF MATTER / DEAD PLANTS
      !
      ! Literature survey shows acetaldehyde yields 3-80 ug/gDW.
      ! Here we use 40ug/gDW.
      ! Assume organic matter has molecular formula CH2O (MW 30.03), which
      ! means 0.4 gC / g plant
      ! Acetaldehyde is 44.05 g/mole, which means 0.55 gC/g acetaldehyde
      ! Therefore 40ug (acetaldehyde) / g (plant)
      ! ==> 55 ug C in acetaldehyde / gC plant
      ! The monthly values of dry leaf matter 
      ! comes from estimates of resp, heterotrophic respiration, from 
      ! Parvada Suntharalingham.  XRESP is in units of g C/m2/s.
      !=================================================================

      ! Convert from [g C in dry leaves/m2/s] to [g C in ACETALDEHYDE/m2/s]
      YIELD_RESP = 55d-6
      ALD2_RESP  = XRESP(I,J) * YIELD_RESP

      ! Convert [g C in acetaldehyde/m2/s] to [atoms C/box/s] 
      ! The 1000 is for changing [g C] to [kg C]
      ACETALDEHYDE =  ALD2_RESP * AREA_M2 * XNUMOL_C / 1000d0

      ! Return to calling program
      END SUBROUTINE EMISS_ALD2_SEN

!------------------------------------------------------------------------------

      SUBROUTINE EMISS_EOH_SEN( I,    J,    ETHANOL )
!
!******************************************************************************
!  Subroutine EMISS_EOH_SEN computes the biogenic emissions of ETHANOL
!  from decaying plants (dbm, 09/29/2008)
!
!  Arguments as Input:
!  ============================================================================
!  (1-2) I, J    (INTEGER) : Lon/lat grid box indices
!
!  Arguments as Input/Output:
!  ============================================================================
!  (2  ) ETHANOL (REAL*8 ) : Total biogenic ethanol emissions from
!        decaying plants [atoms C/box/s]
!
!******************************************************************************
!
      ! References to F90 modules
      USE GRID_MOD, ONLY : GET_AREA_M2
      USE ACETONE_MOD        ! XRESP, XNUMOL_C

      USE CMN_SIZE_MOD

      ! Arguments
      INTEGER, INTENT(IN)    :: I, J 
      REAL*8,  INTENT(OUT)   :: ETHANOL

      ! Local variables
      REAL*8                 :: YIELD_RESP,  EOH_RESP,    AREA_M2

      !=================================================================
      ! EMISS_EOH_SEN begins here!
      !=================================================================

      ! Grid box areas in [m2] and [cm2]
      AREA_M2  = GET_AREA_M2(I, J, 1)

      !=================================================================
      ! (1) BIOGENIC ETHANOL FROM DRY LEAF MATTER / DEAD PLANTS
      !
      ! Use 40ug/gDW.
      ! Assume organic matter has molecular formula CH2O (MW 30.03), which
      ! means 0.4 gC / g plant
      ! Ethanol is 46 g/mole ==> 0.52 gC/g eoh
      ! Therefore 40ug (eoh) / g (plant)
      ! ==> 52 ug C in ethanol / gC plant
      ! The monthly values of dry leaf matter 
      ! comes from estimates of resp, heterotrophic respiration, from 
      ! Parvada Suntharalingham.  XRESP is in units of g C/m2/s.
      !=================================================================

      ! Convert from [g C in dry leaves/m2/s] to [g C in ETHANOL/m2/s]
      YIELD_RESP = 52d-6
      EOH_RESP  = XRESP(I,J) * YIELD_RESP

      ! Convert [g C in ethanol/m2/s] to [atoms C/box/s] 
      ! The 1000 is for changing [g C] to [kg C]
      ETHANOL =  EOH_RESP * AREA_M2 * XNUMOL_C / 1000d0

      ! Return to calling program
      END SUBROUTINE EMISS_EOH_SEN

!------------------------------------------------------------------------------

!------------------------------------------------------------------------------

      SUBROUTINE EMISS_ALD2_OCN( I,    J,    ACETALDEHYDE, State_Met)
!
!******************************************************************************
!  Subroutine EMISS_ALD2_OCN computes the ocean source ACETALDEHYDE
!  from photodegradation of DOM (dbm, 12/08/2008)
!
!  Arguments as Input:
!  ============================================================================
!  (1-2) I, J    (INTEGER) : Lon/lat grid box indices
!
!  Arguments as Input/Output:
!  ============================================================================
!  (2  ) ACETALDEHYDE (REAL*8 ) : Total oceanic acetaldehyde emissions
!        [atoms C/box/s]
!
!***************************************************************************
!  ALD2 air-sea exchange
!  -------------------------------
!  The flux of gas (ALD2) through the air-sea interface can be expressed as:
!
!  Liss & Slater, Nature 247, 181-184, 1974:
!  F = KKl * ( Cg / H - Cl)
!
!  where
!  (a)KKl    :  total transfer velocity expressed on a liquid phase (seawater) 
!               basis, 
!
!  (b) Cg, Cl :  gas (ALD2) concentration in gas and liquid phase, respectively 
!  (c) H      :  dimensionless Henry's law constant 
!                 [mass MOH/volume air]/[mass MOH/volume H2O]
!
!  At equilibrium, Cl = Cg / H
!
!  Gross ocean source: F = KKl * Cl
!  Gross ocean sink  : F = KKl * Cg/H
!
!  Parameters H, KKL are computed in subroutine GET_ALD2_OCN_PARAMS.
!
!  (d)Reference:
!
!             Liss, P.S., and P.G. Slater, Flux of gases across the air-
!                sea interface, Nature, 247, 181-184, 1974
!                Over the Ocean, J. Geophys. Res., 97, 7373-7382, 1992
!
!  -- dbm 12/08/2008
!*****************************************************************************
!
      ! References to F90 modules
      USE GRID_MOD,      ONLY : GET_AREA_CM2
      !evf, replace DAO_MOD with MetState object, 5/21/13
      USE GIGC_State_Met_Mod, ONLY : MetState
      !USE DAO_MOD,       ONLY : ALBD, evf

      USE TIME_MOD,      ONLY : GET_MONTH, ITS_A_NEW_MONTH
      USE BPCH2_MOD,     ONLY : GET_TAU0,        READ_BPCH2, GET_RES_EXT
      USE TRANSFER_MOD, ONLY : TRANSFER_2D

      USE CMN_SIZE_MOD
      !USE CMN_DEP_MOD, evf obsolete

      ! Arguments
      INTEGER, INTENT(IN)    :: I, J 
      REAL*8,  INTENT(OUT)   :: ACETALDEHYDE

            !evf, replace DAO_MOD with MetState object, 5/21/13
      TYPE(MetState), INTENT(IN)  :: State_Met   ! Meteorology State object

      ! Local variables
      REAL*8                 :: AREA_CM2, FOCEAN
      REAL*8                 :: KKL, HENCONST, XTAU
      REAL*8                 :: THIS_ALD2_OML_CONC
      REAL*4                 :: ARRAY(IIPAR,JJPAR,LLTROP)
      INTEGER                :: MONTH
      CHARACTER(LEN=255)     :: FILENAME

      ! Saved variables
      !REAL*8, SAVE           :: ALD2_OML_PROD(IIPAR,JJPAR)
      REAL*8                  :: ALD2_OML_PROD(IIPAR,JJPAR)

      ! molecular weights (kg/mole)
      REAL*8,  PARAMETER     :: MW_ALD2  = 44.0526d-3

      ! Avogaddro's number
      REAL*8,  PARAMETER     :: AVO      = 6.0223D23

      !=================================================================
      ! EMISS_ALD2_OCN begins here!
      !=================================================================
      
      ! Read in ocean mixed layer ALD2 production fields from file
      IF ( ITS_A_NEW_MONTH() .AND. I == 1 .AND. J == 1 ) THEN

         ! Read in the ALD2 values (nM/h)
         filename='/home/efischer' //
     &        '/Acetaldehyde/input_files_for_emily/' //
     &        'ALD2_OML_PROD_2006_nM_PER_HOUR.' //
     &        GET_RES_EXT() // '.bpch'

         MONTH   = GET_MONTH()
         ! Get TAU0 value to index the punch file -- 
         ! use generic year 2006
         XTAU = GET_TAU0( MONTH, 1, 2006 )

        ! read ALD2 production in units of nM/h
         write(6,*) '% Reading ALD2 OML production ',filename
         call READ_BPCH2(filename,'BIOGSRCE', 99, XTAU, 
     &        IIPAR,JJPAR,1,ARRAY(:,:,1))
         call TRANSFER_2D(ARRAY(:,:,1),ALD2_OML_PROD)

      ENDIF

      ! Initialize
      ACETALDEHYDE = 0d0

      ! Grid cell area
      AREA_CM2 = GET_AREA_CM2(I, J,1 )

      ! Fraction of grid box that is ocean
      FOCEAN = 1D0 - State_Met%FRCLND(I,J)

      !===========================================================
      ! Only compute ocean source if there is more than 50% ocean
      ! in the grid box, and if it is not ice (albedo > 0.4)
      ! (mje, rvm, bmy, 11/26/01)
      !===========================================================
      IF ( FOCEAN > 0.5d0 .AND.State_Met%ALBD(I,J) <= 0.4d0 ) THEN

         ! Calculate Henry constant and transfer velocity
         CALL GET_ALD2_OCN_PARAMS( I, J, HENCONST, KKL, State_Met)

         !-------------------------------------------------------
         ! Compute gross ocean source
         ! sea to air FLUX(kg/cm2/s) = KKl(cm/s) * Cl(kg/cm3)
         !-------------------------------------------------------

         ! Assume an OML lifetime for ALD2.  Then the steady state
         ! concentration is given by (production nM/h) * (lifetime h) = nM
         !Write(6,*) '% Assuming ALD2 OML lifetime = 1.25 h'
         ! sensitivity run to examine higher ocean source
         THIS_ALD2_OML_CONC = 0.5d0 * ALD2_OML_PROD( I , J )
         !THIS_ALD2_OML_CONC = 1.05d0 * ALD2_OML_PROD( I , J )
      
         ! Convert to kg/cm3
         Cl = THIS_ALD2_OML_CONC *1d-9 * MW_ALD2 / 1000.0d0

         ! Correct for the fraction of the grid cell that is ocean
         ! and compute flux (kg/cm3 * cm/s = kg/cm2/s)
         ACETALDEHYDE = KKL * FOCEAN * Cl
            
         ! Convert to atomsC / box / second
         ! kg/cm2/s * cm2/box * mole/kg * atomsC/molecule * molec/mole = atomsC/box/s
         ACETALDEHYDE = ACETALDEHYDE * AREA_CM2 * AVO / MW_ALD2 * 2d0

      ENDIF                     ! IF FOCEAN

      ! Return to calling program
      END SUBROUTINE EMISS_ALD2_OCN

!------------------------------------------------------------------------------
      SUBROUTINE ALD2_OCN_SINK(State_Met)
!
!******************************************************************************
!  Subroutine ALD2_OCN_SINK computes the ocean sink of ACETALDEHYDE
!  It is called from chemistry_mod.f (dbm, 12/08/
!
!***************************************************************************
!  ALD2 air-sea exchange
!  -------------------------------
!  The flux of gas (ALD2) through the air-sea interface can be expressed as:
!
!  Liss & Slater, Nature 247, 181-184, 1974:
!  F = KKl * ( Cg / H - Cl)
!
!  where
!  (a)KKl    :  total transfer velocity expressed on a liquid phase (seawater) 
!               basis, 
!
!  (b) Cg, Cl :  gas (ALD2) concentration in gas and liquid phase, respectively 
!  (c) H      :  dimensionless Henry's law constant 
!                 [mass MOH/volume air]/[mass MOH/volume H2O]
!
!  At equilibrium, Cl = Cg / H
!
!  Gross ocean source: F = KKl * Cl
!  Gross ocean sink  : F = KKl * Cg/H
!
!  Parameters H, KKL are computed in subroutine GET_ALD2_OCN_PARAMS.
!
!  (d)Reference:
!
!             Liss, P.S., and P.G. Slater, Flux of gases across the air-
!                sea interface, Nature, 247, 181-184, 1974
!                Over the Ocean, J. Geophys. Res., 97, 7373-7382, 1992
!
!  -- dbm 12/08/2008
!*****************************************************************************
!
      ! References to F90 modules
      USE GRID_MOD,     ONLY : GET_AREA_CM2

      !USE DAO_MOD,      ONLY : ALBD, BXHEIGHT
      !evf, replace DAO_MOD with MetState object, 5/21/13
      USE GIGC_State_Met_Mod, ONLY : MetState
      USE PBL_MIX_MOD,  ONLY : GET_FRAC_UNDER_PBLTOP, GET_PBL_TOP_m
      USE PBL_MIX_MOD,  ONLY : GET_PBL_MAX_L
      USE TRACER_MOD,   ONLY : STT
      USE TRACERID_MOD, ONLY : IDTALD2
      USE TIME_MOD,     ONLY : GET_TS_CHEM
      USE DIAG_MOD,     ONLY : AD46

      USE CMN_SIZE_MOD
      !USE CMN_DEP_MOD, evf
      USE CMN_DIAG_MOD

            !evf, replace DAO_MOD with MetState object, 5/21/13
      TYPE(MetState), INTENT(IN)  :: State_Met   ! Meteorology State object

      ! Local variables
      INTEGER                :: I, J, L, PBL_MAX
      REAL*8                 :: AREA_CM2, FOCEAN
      REAL*8                 :: KKL, HENCONST, THIK, DV, LF
      REAL*8                 :: STT0, F_UNDER_TOP, FREQRUR
      REAL*8                 :: E_RKT, DTCHEM

      ! molecular weights (kg/mole)
      REAL*8,  PARAMETER     :: MW_ALD2  = 44.0526d-3
      ! Avogaddro's number
      REAL*8,  PARAMETER     :: AVO      = 6.0223D23

      !=================================================================
      ! ALD2_OCN_SINK begins here!
      !=================================================================

      write(6,*) '% ALD2_MOD: Computing ocean sink'

      ! Max PBL height for globe
      PBL_MAX = GET_PBL_MAX_L()
      
      ! DTCHEM is the chemistry timestep in seconds
      DTCHEM = GET_TS_CHEM() * 60d0 

      DO J = 1, JJPAR

      ! Grid cell area

      DO I = 1, IIPAR
      AREA_CM2 = GET_AREA_CM2(I, J,1)
         ! Fraction of grid box that is ocean
         FOCEAN = 1D0 - State_Met%FRCLND(I,J)

         !===========================================================
         ! Only compute ocean sink if there is more than 50% ocean
         ! in the grid box, and if it is not ice (albedo > 0.4)
         ! (mje, rvm, bmy, 11/26/01)
         !===========================================================
         IF ( FOCEAN > 0.5d0 .AND. State_Met%ALBD(I,J) <= 0.4d0 ) THEN

            ! Calculate Henry constant and transfer velocity
            CALL GET_ALD2_OCN_PARAMS( I, J, HENCONST, KKL, State_Met )
            !-------------------------------------------------------
            ! Compute gross ocean sink
            ! air to sea FLUX = KKl * Cg/H
            !-------------------------------------------------------

            ! Net air-to-sea deposition velocity
            ! cm/s * unitless/ unitless = cm/s
            DV = KKL * FOCEAN / HENCONST

            ! Mixed layer depth (m).  
            ! Make sure that it is greater than or equal 
            ! to the thickness of the 1st layer
            THIK    = MAX(GET_PBL_TOP_m(I, J), 
     &                    State_Met%BXHEIGHT(I, J, 1))

            ! Loss frequency averaged over MBL (s-1)
            LF = ( DV / 100d0 ) / THIK 
                   
            ! Distribute loss through MBL
            DO L = 1, PBL_MAX           

               ! Initial concentration of ald2
               STT0 = STT(I,J,L,IDTALD2)

               ! Fraction under the PBL top
               F_UNDER_TOP = GET_FRAC_UNDER_PBLTOP( I, J, L )              
                  
               ! Account for frac w/ in the PBL 
               FREQRUR = LF * F_UNDER_TOP 

               ! Skip if no loss occurs
               IF ( FREQRUR > 0d0 ) THEN 

                  ! Pre-compute the exponential
                  E_RKT = EXP( -FREQRUR * DTCHEM )

                  ! M = Mo * e^(-kt)
                  STT(I,J,L,IDTALD2) = STT(I,J,L,IDTALD2) * E_RKT

                  !==================================================
                  ! Use ND46 for this
                  ! Ocean uptake is process #19
                  ! kg * (molec/mole) * (atomsC/molec) 
                  !    / (cm2 * s * kg / mole) = atomsC/cm2/s
                  !
                  ! bug here - should multiply by (1/12) instead of 
                  ! (2/MW_ALD2) since STT is kgC.
                  !
                  ! Now do it correctly:
                  ! kgC * (1 moleC/12E-3 kgC) * (atomsC/moleC) / (cm2 * s)
                  ! = atomsC/cm2/s
                  !==================================================

                  IF ( ND46 > 0 ) THEN
                     AD46(I,J,19) = AD46(I,J,19) + 
     &                    (STT0 - STT (I,J,L,IDTALD2) ) * 
     &                    AVO * 1d0/(AREA_CM2 * DTCHEM * 12D-3)

                  ENDIF
                           
               ENDIF            ! IF FREQRUR > 0
            ENDDO               ! L
         ENDIF                  ! IF FOCEAN and ALBD           
      ENDDO                     ! I
      ENDDO                     ! J
      END SUBROUTINE ALD2_OCN_SINK

!------------------------------------------------------------------------------

      SUBROUTINE GET_ALD2_OCN_PARAMS( I, J, HENCONST, KKL, State_Met )

!******************************************************************************
!  Subroutine GET_ALD2_OCN_PARAMS calculates the Henry constant and air-to-sea
!  transfer velocity, H and KKL. (dbm, 12/08/2008)
!
!  NOTES:
!
!  Liss & Slater, Nature 247, 181-184, 1974:
!  F = KKl * ( Cg / H - Cl)
!
!  where
!  (a)KKl    :  total transfer velocity expressed on a liquid phase (seawater) 
!               basis, 
!
!                 1      1        1
!               ----- = ---  + --------
!                KKl     kl     H * kg
!
!               where 
!                (a.1) kl : transfer velocity in liquid phase (seawater) (cm/hr)
!                (a.2) kg : transfer velocity in gas phase
!
!                      KL, KG are computed according to the recommendations of Nightingale and Asher.
!                      See below for details and references.
!
!  Henry constant computed here is dimensionless [conc. in air] / [conc. in water].
!                        1
!               H = ------------- 
!                    KH * Rg * T
!
!               where
!                (c.1) KH  :  Henry's law constant 
!                              [X]aq / P(X) (moles / litr / atm) 
!                (c.2) Rg  :  Universal gas constant
!                (c.3) T   :  Temperature
!
!  dbm, 01/15/2009: note that in the aqueous phase acetaldehyde is partially present
!        in the hydrated form as the gem-diol CH3CH(OH)2.  Therefore the appropriate
!        KH is the apparent Henry's law constant, where
!        [X]aq is the total aldehyde in the aqueous phase, hydrated + unhydrated.
!
!  dbm, 01/15/2009: also note that the hydration effect increases the rate of aqueous
!        diffusion.  This means that:
!
!                 1          1            1
!               ----- =  ---------  +   ------
!                KKl    (alpha)(kl)     H * kg
!
!       Where alpha is the enhancement factor.  Use alpha = 2.0 following
!       Zhou and Mopper [1997].
!******************************************************************************
!
      ! References to F90 modules
      !USE DAO_MOD,      ONLY : TS ! evf, replace dao_mod, 5/21/13
      USE GIGC_State_Met_Mod, ONLY : MetState 

      IMPLICIT NONE

      ! Arguments
      INTEGER, INTENT(IN)  :: I, J
      REAL*8,  INTENT(OUT) :: HENCONST
      REAL*8,  INTENT(OUT) :: KKL

            !evf, replace DAO_MOD with MetState object, 5/21/13
      TYPE(MetState), INTENT(IN)  :: State_Met   ! Meteorology State object

      ! Local variables

      ! Coefficients for fitting the Schmidt number for acetaldehyde
      ! Derived using polynomial fit (code provided by qli, same as used
      ! for acetone, methanol)
      ! and partial molal volume of acetaldehyde at its normal boiling
      ! temperature (51.8 cm3/g/mole) calculated using Le Bas method
      ! see "The Properties of Gases and Liquids", Reid, Prausnitz, Sherwood.
 
      REAL*8,  PARAMETER     :: A0 = 2581.709d0
      REAL*8,  PARAMETER     :: A1 = -106.9671d0
      REAL*8,  PARAMETER     :: A2 = 1.73263d0
      REAL*8,  PARAMETER     :: A3 = -0.0110773d0

      REAL*8                 :: KL, KL600, KG, KH
      REAL*8                 :: SC, U, TC

      ! External functions
      REAL*8, EXTERNAL       :: SFCWINDSQR



      !=================================================================
      !  GET_ALD2_OCN_PARAMS begins here!
      !=================================================================

      !-------------------------------------------------------
      ! Compute H
      !-------------------------------------------------------

      ! From Sander compilation:
      ! kH (M/atm) = 13 * exp ( 5750 * ( 1/T - 1/298. ) )
      ! Multiply by R * T to get dimensionless Henry's constant
      !     HENCONST in (g / cm3 water) / (g / cm3 air)
      !KH       = 1.3d1 * EXP( 5750d0 * ( 1d0/TS(I,J) - 1d0/298d0 ) )

      ! Use value from Betterton and Hoffmann - apparent Henry's law
      ! coefficient
      KH       = 1.1d1 * EXP( 6300d0 *  
     &                    ( 1d0/State_Met%TS(I,J) - 1d0/298d0 ) )
      HENCONST = 0.08206d0 * 298d0 * KH

      ! For Liss formulation, want Henry constant exressed as
      ! concentration in air/concentration in water.
      ! Take reciprocal.  
      HENCONST = 1d0/HENCONST
                  
      !-------------------------------------------------------
      ! Compute schmidt number
      !-------------------------------------------------------
      ! Magnitude of surface wind [m/s]
      ! SFCWINDSQR(I,J) is needed since this will compute the 
      ! square of the surface wind correctly for all GEOS models 
      U     = SQRT( SFCWINDSQR( I, J, State_Met%U10M, State_Met%V10M ) )

      ! TC is temperature in Celsius
      TC    = State_Met%TS(I,J) - 273.15d0
      ! Bug fix: make sure -5 <= TC <= 30 (mje, rvm, bmy, 11/26/01) 
      TC    = MIN( MAX( TC, -5d0 ), 30d0 )

      ! SC is Schmidt # for ALD2 [unitless]
      SC    = A0 + A1*TC + A2*TC**2 + A3*TC**3

      !-------------------------------------------------------
      ! Liquid-phase transfer velocity
      ! Nightingale et al., Global Biogeochem. Cyc., 14, 373-387 (2000)
      ! 
      ! KL = KL600 * ( Sc / 600 ) ^-0.5
      !
      ! Parameterized relationship for available observations:
      !        KL600 = 0.222*U^2 + 0.333*U 
      !                U = wind speed @ 10m (m/s)
      !
      ! Units of KL and KL600 are cm/h
      !-------------------------------------------------------

      KL600 = ( 0.222d0 * U * U ) + ( 0.333D0 * U )
      KL    = KL600 * ( SC / 600d0 )**( -.5d0 )
            
      !-------------------------------------------------------
      ! Compute kg (conductance for mass transfer in gas phase)
      !
      ! Asher: The sea-surface microlayer and its effect on global air-sea
      !        gas transfer (1997)
      ! In   : The sea surface and global change, P.S. Liss & R. A. Duce, Eds., 1997.
      !
      ! Kg = (5.2E-5 + 3.2E-3 * u) * SQRT(WTH2O / WTALD2)
      !                  
      !                  where WTH2O =18.0 and WTMALD2 =44.1 are the molecular weight 
      !                  of H2O and ALD2 in (g/mole), respectively.
      ! Units of Kg are m/s.
      !-------------------------------------------------------
      KG    = SQRT( 18.015d0 / 44.053d0 ) * ( 5.2d-5 + 3.2d-3 * U )
      ! Convert from [m/s] to [cm/hr]
      KG    = KG * 3600d0 * 100d0

      !-------------------------------------------------------
      ! Compute KKl (total air to sea transfer velocity expressed on a 
      ! liquid phase basis)
      !                 1          1            1
      !               ----- =  ---------  +  -------
      !                KKl    (alpha)(kl)     H * kg
      !-------------------------------------------------------

      ! KKL is the air-to-sea transfer velocity (Liss and Slater 1974)
      ! Divide KKL by 3600 to convert from [cm/hr] to [cm/s]
!      KKL   = 1d0 / ( 1d0/KL + 1d0/( HENCONST * KG ) )
      KKL   = 1d0 / ( 1d0/(2d0*KL) + 1d0/( HENCONST * KG ) )
      KKL   = KKL / 3600d0


      ! Return to calling program
      END SUBROUTINE GET_ALD2_OCN_PARAMS
!------------------------------------------------------------------------------
      END MODULE ALD2_MOD
