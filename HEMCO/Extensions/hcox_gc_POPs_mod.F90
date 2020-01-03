!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: hcox_gc_POPs_mod.F90
!
! !DESCRIPTION: Defines the HEMCO extension for the GEOS-Chem persistent
!   organic pollutants (POPs) specialty simulation.
!\\
!\\
! !INTERFACE:
!
MODULE HCOX_GC_POPs_Mod
!
! !USES:
!
  USE HCO_Error_Mod
  USE HCO_Diagn_Mod
  USE HCO_State_Mod,  ONLY : HCO_State   ! Derived type for HEMCO state
  USE HCOX_State_Mod, ONLY : Ext_State   ! Derived type for External state

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC  :: HcoX_GC_POPs_Run
  PUBLIC  :: HcoX_GC_POPs_Init
  PUBLIC  :: HcoX_Gc_POPs_Final
!
! !PRIVATE MEMBER FUNCTIONS:
!
  PRIVATE :: VEGEMISPOP
  PRIVATE :: LAKEEMISPOP
  PRIVATE :: SOILEMISPOP
  PRIVATE :: IS_LAND
  PRIVATE :: IS_ICE
!
! !REMARKS:
!
!  References:
!  ============================================================================
!  (1 ) Zhang, Y., and S. Tao, Global atmospheric emission inventory of
!        polycyclic aromatic hydrocarbons (PAHs) for 2004. Atm Env, 43, 812-819,
!        2009.
!  (2 ) Friedman, C.L, and N.E. Selin, Long-Range Atmospheric Transport of
!        Polycyclic Aromatic Hydrocarbons: A Global 3-D Model Analysis
!        Including Evaluation of Arctic Sources, Environ. Sci. Technol., 46(17),
!        9501-9510, 2012.
!  (3 ) Friedman, C.L., Y. Zhang, and N.E. Selin, Climate change and
!        emissions impacts on atmospheric PAH transport to the Arctic, Environ.
!        Sci. Technol., 48, 429-437, 2014.
!  (4 ) Friedman, C.L., J.R. Pierce, and N.E. Selin, Assessing the influence of
!        secondary organic versus primary carbonaceous aerosols on long-range
!        atmospheric polycyclic aromatic hydrocarbon transport, Environ. Sci.
!        Technol., 48(6), 3293-3302, 2014.
!
! !REVISION HISTORY:
!  20 Sep 2010 - N.E. Selin    - Initial Version
!  04 Jan 2011 - C.L. Friedman - Expansion on initial version
!  19 Aug 2014 - M. Sulprizio  - Now a HEMCO extension
!  18 Aug 2015 - M. Sulprizio  - Add VEGEMISPOP, LAKEEMISPOP, and SOILEMISPOP
!                                routines from new land_pops_mod.F written by
!                                C.L. Friedman.
!  24 Aug 2017 - M. Sulprizio  - Remove support for GCAP
!  25 Jan 2019 - M. Sulprizio  - Add instance wrapper
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
  REAL(hp), PARAMETER           :: SMALLNUM    = 1e-20_hp
!
! !PRIVATE TYPES:
!
  TYPE :: MyInst
   ! Fields required by module
   INTEGER               :: Instance
   INTEGER               :: ExtNr            ! HEMCO Extension number
   INTEGER               :: IDTPOPPOCPO      ! Index # for POPPOC tracer
   INTEGER               :: IDTPOPPBCPO      ! Index # for POPPBC tracer
   INTEGER               :: IDTPOPG          ! Index # for POPG   tracer

   ! Pointers to emission arrays read from disk
   REAL(sp), POINTER     :: POP_TOT_EM(:,:) => NULL() ! [kg/m2/s]
   REAL(sp), POINTER     :: POP_SURF(:,:)   => NULL() ! [kg]
   REAL(sp), POINTER     :: C_OC(:,:,:)     => NULL() ! [kg]
   REAL(sp), POINTER     :: C_BC(:,:,:)     => NULL() ! [kg]
   REAL(sp), POINTER     :: F_OC_SOIL(:,:)  => NULL() ! [kg/m2]

   ! Calculated emissions of OC-phase, BC-phase, and gas-phase POPs [kg/m2/s]
   REAL(hp), POINTER     :: EPOP_G    (:,:,:)
   REAL(hp), POINTER     :: EPOP_OC   (:,:,:)
   REAL(hp), POINTER     :: EPOP_BC   (:,:,:)
   REAL(hp), POINTER     :: EPOP_VEG  (:,:)
   REAL(hp), POINTER     :: EPOP_LAKE (:,:)
   REAL(hp), POINTER     :: EPOP_SOIL (:,:)
   REAL(hp), POINTER     :: EPOP_OCEAN(:,:)
   REAL(hp), POINTER     :: EPOP_SNOW (:,:)

   ! For diagnostics
   CHARACTER(LEN=63)     :: DiagnName
   REAL(hp), POINTER     :: SUM_OC_EM    (:,:)
   REAL(hp), POINTER     :: SUM_BC_EM    (:,:)
   REAL(hp), POINTER     :: SUM_G_EM     (:,:)
   REAL(hp), POINTER     :: SUM_OF_ALL   (:,:)
   REAL(hp), POINTER     :: EMIS_SOIL    (:,:)
   REAL(hp), POINTER     :: FLUX_SOIL2AIR(:,:)
   REAL(hp), POINTER     :: FLUX_AIR2SOIL(:,:)
   REAL(hp), POINTER     :: FUG_SOILAIR  (:,:)
   REAL(hp), POINTER     :: EMIS_LAKE    (:,:)
   REAL(hp), POINTER     :: FLUX_LAKE2AIR(:,:)
   REAL(hp), POINTER     :: FLUX_AIR2LAKE(:,:)
   REAL(hp), POINTER     :: FUG_LAKEAIR  (:,:)
   REAL(hp), POINTER     :: EMIS_LEAF    (:,:)
   REAL(hp), POINTER     :: FLUX_LEAF2AIR(:,:)
   REAL(hp), POINTER     :: FLUX_AIR2LEAF(:,:)
   REAL(hp), POINTER     :: FUG_LEAFAIR  (:,:)

   TYPE(MyInst), POINTER :: NextInst => NULL()
  END TYPE MyInst

  ! Pointer to instances
  TYPE(MyInst), POINTER  :: AllInst => NULL()

CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCOX_GC_POPs_run
!
! !DESCRIPTION: Subroutine HcoX\_Gc\_POPs\_Run computes emissions of OC-phase,
!  BC-phase, and gas-phase POPs for the GEOS-Chem POPs specialty simulation.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCOX_GC_POPs_Run( ExtState, HcoState, RC )
!
! !USES:
!
    ! HEMCO modules
    USE HCO_EmisList_Mod,  ONLY : HCO_GetPtr
    USE HCO_FluxArr_Mod,   ONLY : HCO_EmisAdd
    USE HCO_Clock_Mod,     ONLY : HcoClock_First
!
! !INPUT PARAMETERS:
!
    TYPE(Ext_State),  POINTER       :: ExtState    ! Options for POPs sim
    TYPE(HCO_State),  POINTER       :: HcoState    ! HEMCO state
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(INOUT) :: RC          ! Success or failure?
!
! !REMARKS:
!  This code is based on routine EMISSPOPS in prior versions of GEOS-Chem.
!
! !REVISION HISTORY:
!  20 Sep 2010 - N.E. Selin  - Initial Version based on EMISSMERCURY
!  29 Nov 2012 - M. Payer    - Replaced all met field arrays with State_Met
!                              derived type object
!  13 Dec 2012 - R. Yantosca - Remove reference to obsolete CMN_DEP_mod.F
!  25 Mar 2013 - R. Yantosca - Now accept am_I_Root, Input_Opt, State_Chm, RC
!  14 Apr 2014 - R. Yantosca - Prevent div-by-zero error w/ SUM_OF_ALL
!  19 Aug 2014 - M. Sulprizio- Now a HEMCO extension
!  07 Jan 2016 - E. Lundgren - Update molar gas constant to NIST 2014
!  26 Oct 2016 - R. Yantosca - Don't nullify local ptrs in declaration stmts
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
    ! Universal gas constant for adjusting KOA for temp: 8.3144598 [J/mol/K]
    REAL(hp), PARAMETER :: R          = 8.3144598e+0_hp

    ! Density of octanol, needed for partitioning into OC: 820 [kg/m^3]
    REAL(hp), PARAMETER :: DENS_OCT   = 82e+1_hp

    ! Density of BC, needed for partitioning onto BC: 1 [kg/L] or 1000 [kg/m^3]
    ! From Lohmann and Lammel, Environ. Sci. Technol., 2004, 38:3793-3803.
    REAL(hp), PARAMETER :: DENS_BC    = 1e+3_hp

!
! !LOCAL VARIABLES:
!
    INTEGER             :: I, J, L
    INTEGER             :: PBL_MAX
    INTEGER             :: MONTH,            YEAR
    REAL(hp)            :: F_OF_PBL,         TK
    REAL(hp)            :: T_POP
    REAL(hp)            :: C_OC1,            C_BC1
    REAL(hp)            :: C_OC2,            C_BC2
    REAL(hp)            :: F_POP_OC,         F_POP_BC
    REAL(hp)            :: F_POP_G,          AIR_VOL
    REAL(hp)            :: KOA_T,            KBC_T
    REAL(hp)            :: KOC_BC_T,         KBC_OC_T
    REAL(hp)            :: VR_OC_AIR,        VR_OC_BC
    REAL(hp)            :: VR_BC_AIR,        VR_BC_OC
    REAL(hp)            :: SUM_F
    REAL(hp)            :: OC_AIR_RATIO,     OC_BC_RATIO
    REAL(hp)            :: BC_AIR_RATIO,     BC_OC_RATIO
    REAL(hp)            :: FRAC_SNOW_OR_ICE, FRAC_SNOWFREE_LAND
    REAL(hp)            :: FRAC_LEAF, FRAC_LAKE, FRAC_SOIL
!    LOGICAL, SAVE       :: FIRST = .TRUE.
    LOGICAL             :: aIR
    LOGICAL             :: IS_SNOW_OR_ICE,   IS_LAND_OR_ICE
    CHARACTER(LEN=255)  :: MSG

    ! Delta H for POP [kJ/mol]. Delta H is enthalpy of phase transfer
    ! from gas phase to OC. For now we use Delta H for phase transfer
    ! from the gas phase to the pure liquid state.
    ! For PHENANTHRENE:
    ! this is taken as the negative of the Delta H for phase transfer
    ! from the pure liquid state to the gas phase (Schwarzenbach,
    ! Gschwend, Imboden, 2003, pg 200, Table 6.3), or -74000 [J/mol].
    ! For PYRENE:
    ! this is taken as the negative of the Delta H for phase transfer
    ! from the pure liquid state to the gas phase (Schwarzenbach,
    ! Gschwend, Imboden, 2003, pg 200, Table 6.3), or -87000 [J/mol].
    ! For BENZO[a]PYRENE:
    ! this is also taken as the negative of the Delta H for phase transfer
    ! from the pure liquid state to the gas phase (Schwarzenbach,
    ! Gschwend, Imboden, 2003, pg 452, Prob 11.1), or -110,000 [J/mol]
    REAL(hp)            :: DEL_H

    ! KOA_298 for partitioning of gas phase POP to atmospheric OC
    ! KOA_298 = Cpop in octanol/Cpop in atmosphere at 298 K
    ! For PHENANTHRENE:
    ! log KOA_298 = 7.64, or 4.37*10^7 [unitless]
    ! For PYRENE:
    ! log KOA_298 = 8.86, or 7.24*10^8 [unitless]
    ! For BENZO[a]PYRENE:
    ! log KOA_298 = 11.48, or 3.02*10^11 [unitless]
    ! (Ma et al., J. Chem. Eng. Data, 2010, 55:819-825).
    REAL(hp)            :: KOA_298

    ! KBC_298 for partitioning of gas phase POP to atmospheric BC
    ! KBC_298 = Cpop in black carbon/Cpop in atmosphere at 298 K
    ! For PHENANTHRENE:
    ! log KBC_298 = 10.0, or 1.0*10^10 [unitless]
    ! For PYRENE:
    ! log KBC_298 = 11.0, or 1.0*10^11 [unitless]
    ! For BENZO[a]PYRENE:
    ! log KBC_298 = 13.9, or 7.94*10^13 [unitless]
    ! (Lohmann and Lammel, EST, 2004, 38:3793-3802)
    REAL(hp)            :: KBC_298

    ! Pointers
    REAL(hp),     POINTER :: Arr3D(:,:,:)
    TYPE(MyInst), POINTER :: Inst

    !=======================================================================
    ! HCOX_GC_POPs_RUN begins here!
    !=======================================================================

    ! Return if extension not turned on
    IF ( ExtState%GC_POPs <= 0 ) RETURN

    ! Enter
    CALL HCO_ENTER( HcoState%Config%Err, 'HCOX_GC_POPs_Run (hcox_gc_POPs_mod.F90)', RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Get instance
    Inst => NULL()
    CALL InstGet ( ExtState%GC_POPs, Inst, RC )
    IF ( RC /= HCO_SUCCESS ) THEN
       WRITE(MSG,*) 'Cannot find GC_POPs instance Nr. ', ExtState%GC_POPs
       CALL HCO_ERROR(HcoState%Config%Err,MSG,RC)
       RETURN
    ENDIF

    DEL_H   = ExtState%POP_DEL_H
    KOA_298 = ExtState%POP_KOA
    KBC_298 = ExtState%POP_KBC
    Arr3D   => NULL()

    !=======================================================================
    ! Get pointers to gridded data imported through config. file
    !=======================================================================
    IF ( HcoClock_First(HcoState%Clock,.TRUE.) ) THEN

       CALL HCO_GetPtr( HcoState, 'TOT_POP',     Inst%POP_TOT_EM, RC )
       IF ( RC /= HCO_SUCCESS ) RETURN

       CALL HCO_GetPtr( HcoState, 'GLOBAL_OC',   Inst%C_OC,       RC )
       IF ( RC /= HCO_SUCCESS ) RETURN

       CALL HCO_GetPtr( HcoState, 'GLOBAL_BC',   Inst%C_BC,       RC )
       IF ( RC /= HCO_SUCCESS ) RETURN

       CALL HCO_GetPtr( HcoState, 'SURF_POP',    Inst%POP_SURF,   RC )
       IF ( RC /= HCO_SUCCESS ) RETURN

       CALL HCO_GetPtr( HcoState, 'SOIL_CARBON', Inst%F_OC_SOIL,  RC )
       IF ( RC /= HCO_SUCCESS ) RETURN

       ! Convert F_OC_SOIL from kg/m2 to fraction
       DO J=1, HcoState%NY
       DO I=1, HcoState%NX

          ! Assume most of carbon mass extends to 5 cm and calculate
          ! concentration in kg/kg
          ! For now, assume a mean soil bulk density of 1300 kg/m3 similar to
          ! McLachlan 2002 to calculate a dry weight fraction
          Inst%F_OC_SOIL(I,J) = Inst%F_OC_SOIL(I,J) / 30e-2_hp / 13e+2_hp

       ENDDO
       ENDDO

!       FIRST = .FALSE.

    ENDIF

    ! Maximum extent of the PBL [model level]
    IF ( .NOT. ASSOCIATED(ExtState%PBL_MAX) ) THEN
       CALL HCO_ERROR ( HcoState%Config%Err, 'PBL_MAX not defined in ExtState!', RC )
       RETURN
    ELSE
       PBL_MAX = DBLE( ExtState%PBL_MAX )
    ENDIF

    !=================================================================
    ! Call emissions routines for revolatilization fluxes from surfaces
    ! Assume all re-emited POPs are in the gas phase until partitioning
    ! to ambient OC and BC in the boundary layer.
    ! Re-emission flux/mass depends on type of POP
    ! First draft, CLF, 28 Aug 2012
    !=================================================================
    CALL SOILEMISPOP( Inst%POP_SURF, Inst%F_OC_SOIL, Inst%EPOP_SOIL, &
                      HcoState, ExtState, Inst )
    CALL LAKEEMISPOP( Inst%POP_SURF, Inst%EPOP_LAKE,                 &
                      HcoState, ExtState, Inst )
    CALL VEGEMISPOP ( Inst%POP_SURF, Inst%EPOP_VEG,                  &
                      HcoState, ExtState, Inst )

    Inst%EPOP_SNOW  = 0e+0_hp
    Inst%EPOP_OCEAN = 0e+0_hp

    ! Loop over grid boxes
    DO J = 1, HcoState%Ny
    DO I = 1, HcoState%Nx

       F_OF_PBL = 0e+0_hp
       T_POP    = 0e+0_hp

       ! Here, save the total from the emissions array
       ! into the T_POP variable [kg/m2/s]
       T_POP = Inst%POP_TOT_EM(I,J)

       ! Now add revolatilization (secondary) emissions to primary [kg/m2/s]
       T_POP = T_POP + Inst%EPOP_VEG(I,J) + Inst%EPOP_LAKE(I,J) + &
                       Inst%EPOP_SOIL(I,J) !+ &
!                      Inst%EPOP_SNOW(I,J) + Inst%EPOP_OCEAN(I,J)

       !====================================================================
       ! Apportion total POPs emitted to gas phase, OC-bound, and BC-bound
       ! emissions (clf, 2/1/2011)
       ! Then partition POP throughout PBL; store into STT [kg]
       ! Now make sure STT does not underflow (cdh, bmy, 4/6/06; eck 9/20/10)
       !====================================================================

       ! Loop up to max PBL level
       DO L = 1, PBL_MAX

          ! Get temp [K]
          TK = ExtState%TK%Arr%Val(I,J,L)

          ! Define temperature-dependent partition coefficients:
          ! KOA_T, the octanol-air partition coeff at temp T [unitless]
          KOA_T = KOA_298 * EXP((-DEL_H/R) * ((1e+0_hp/TK) - (1e+0_hp/298e+0_hp)))

          ! Define KBC_T, the BC-air partition coeff at temp T [unitless]
          ! TURN OFF TEMPERATURE DEPENDENCY FOR SENSITIVITY ANALYSIS
          KBC_T = KBC_298 * EXP((-DEL_H/R) * ((1e+0_hp/TK) - (1e+0_hp/298e+0_hp)))

          ! Define KOC_BC_T, theoretical OC-BC part coeff at temp T [unitless]
          KOC_BC_T = KOA_T / KBC_T

          ! Define KBC_OC_T, theoretical BC_OC part coeff at temp T [unitless]
          KBC_OC_T = 1d0 / KOC_BC_T

          ! Get monthly mean OC and BC concentrations [kg/box]
          C_OC1    = Inst%C_OC(I,J,L)
          C_BC1    = Inst%C_BC(I,J,L)

          ! Make sure OC is not negative
          C_OC1    = MAX( C_OC1, 0e+0_hp )

          ! Convert C_OC and C_BC units to volume per box
          ! [m^3 OC or BC/box]
          C_OC2    = C_OC1 / DENS_OCT
          C_BC2    = C_BC1 / DENS_BC

          ! Get air volume (m^3)
          AIR_VOL  = ExtState%AIRVOL%Arr%Val(I,J,L)

          ! Define volume ratios:
          ! VR_OC_AIR = volume ratio of OC to air [unitless]
          VR_OC_AIR   = C_OC2 / AIR_VOL

          ! VR_OC_BC  = volume ratio of OC to BC [unitless]
          VR_OC_BC    = C_OC2 / C_BC2

          ! VR_BC_AIR = volume ratio of BC to air [unitless]
          VR_BC_AIR   = VR_OC_AIR / VR_OC_BC

          ! VR_BC_OC  = volume ratio of BC to OC [unitless]
          !VR_BC_OC(I,J,L)    = 1d0 / VR_OC_BC(I,J,L)
          VR_BC_OC    = 1d0 / VR_OC_BC

          ! Redefine fractions of total POPs in box (I,J,L) that are OC-phase,
          ! BC-phase, and gas phase with new time step (should only change if
          ! temp changes or OC/BC concentrations change)
          OC_AIR_RATIO = 1e+0_hp / (KOA_T    * VR_OC_AIR)
          OC_BC_RATIO  = 1e+0_hp / (KOC_BC_T * VR_OC_BC)

          BC_AIR_RATIO = 1e+0_hp / (KBC_T    * VR_BC_AIR)
          BC_OC_RATIO  = 1e+0_hp / (KBC_OC_T * VR_BC_OC)

          ! If there are zeros in OC or BC concentrations, make sure they
          ! don't cause problems with phase fractions
          IF ( C_OC1 > SMALLNUM .and. C_BC1 > SMALLNUM ) THEN
             F_POP_OC  = 1e+0_hp / (1e+0_hp + OC_AIR_RATIO + OC_BC_RATIO)
             F_POP_BC  = 1e+0_hp / (1e+0_hp + BC_AIR_RATIO + BC_OC_RATIO)

          ELSE IF ( C_OC1 > SMALLNUM .and. C_BC1 .le. SMALLNUM ) THEN
             F_POP_OC  = 1e+0_hp / (1e+0_hp + OC_AIR_RATIO)
             F_POP_BC  = SMALLNUM

          ELSE IF ( C_OC1 .le. SMALLNUM .and. C_BC1 > SMALLNUM ) THEN
             F_POP_OC  = SMALLNUM
             F_POP_BC  = 1e+0_hp / (1e+0_hp + BC_AIR_RATIO)

          ELSE IF ( C_OC1 .le. SMALLNUM .and. C_BC1 .le. SMALLNUM) THEN
             F_POP_OC = SMALLNUM
             F_POP_BC = SMALLNUM
          ENDIF

          ! Gas-phase:
          F_POP_G   = 1e+0_hp - F_POP_OC - F_POP_BC

          ! Check that sum of fractions equals 1
          SUM_F = F_POP_OC + F_POP_BC + F_POP_G

          ! Fraction of PBL that box (I,J,L) makes up [unitless]
          F_OF_PBL = ExtState%FRAC_OF_PBL%Arr%Val(I,J,L)

          ! Calculate rates of POP emissions in each phase [kg/m2/s]
          ! OC-phase:
          Inst%EPOP_OC(I,J,L) = F_POP_OC * F_OF_PBL * T_POP

          ! BC-phase
          Inst%EPOP_BC(I,J,L) = F_POP_BC * F_OF_PBL * T_POP

          ! Gas-phase
          Inst%EPOP_G(I,J,L)  = F_POP_G  * F_OF_PBL * T_POP

       ENDDO

!-----------------------------------------------------------------------------
! This code is not actually used (mps, 8/24/15)
!       !==================================================================
!       ! Sum different POPs emissions phases (OC, BC, and gas phase)
!       ! through bottom layer to top of PBL for storage in ND53 diagnostic
!       !==================================================================
!
!       Inst%SUM_OC_EM(I,J) =  SUM(Inst%EPOP_OC(I,J,1:PBL_MAX))
!       Inst%SUM_BC_EM(I,J) =  SUM(Inst%EPOP_BC(I,J,1:PBL_MAX))
!       Inst%SUM_G_EM(I,J)  =  SUM(Inst%EPOP_G(I,J,1:PBL_MAX))
!
!       Inst%SUM_OF_ALL(I,J) = Inst%SUM_OC_EM(I,J) + Inst%SUM_BC_EM(I,J) + &
!                              Inst%SUM_G_EM(I,J)
!
!       ! Check that sum thru PBL is equal to original emissions array
!       ! NOTE: Prevent div-by-zero floating point error (bmy, 4/14/14)
!       IF ( Inst%SUM_OF_ALL(I,J) > 0e+0_hp ) THEN
!          Inst%SUM_OF_ALL(I,J) = Inst%POP_TOT_EM(I,J) / Inst%SUM_OF_ALL(I,J)
!       ENDIF
!-----------------------------------------------------------------------------

    ENDDO
    ENDDO

    !=======================================================================
    ! Add POPs emissions to HEMCO data structure & diagnostics
    !=======================================================================

    !----------------------
    ! OC-PHASE EMISSIONS
    !----------------------
    IF ( Inst%IDTPOPPOCPO > 0 ) THEN

       ! Add flux to emissions array
       Arr3D => Inst%EPOP_OC(:,:,:)
       CALL HCO_EmisAdd( HcoState, Arr3D, Inst%IDTPOPPOCPO, &
                         RC, ExtNr=Inst%ExtNr )
       Arr3D => NULL()
       IF ( RC /= HCO_SUCCESS ) THEN
          CALL HCO_ERROR( HcoState%Config%Err, 'HCO_EmisAdd error: EPOP_OC', RC )
          RETURN
       ENDIF
    ENDIF

    !----------------------
    ! BC-PHASE EMISSIONS
    !----------------------
    IF ( Inst%IDTPOPPBCPO > 0 ) THEN

       ! Add flux to emissions array
       Arr3D => Inst%EPOP_BC(:,:,:)
       CALL HCO_EmisAdd( HcoState, Arr3D, Inst%IDTPOPPBCPO, &
                         RC, ExtNr=Inst%ExtNr )
       Arr3D => NULL()
       IF ( RC /= HCO_SUCCESS ) THEN
          CALL HCO_ERROR( HcoState%Config%Err, 'HCO_EmisAdd error: EPOP_BC', RC )
          RETURN
       ENDIF
    ENDIF

    !----------------------
    ! GASEOUS EMISSIONS
    !----------------------
    IF ( Inst%IDTPOPG > 0 ) THEN

       ! Add flux to emissions array
       Arr3D => Inst%EPOP_G(:,:,:)
       CALL HCO_EmisAdd( HcoState, Arr3D, Inst%IDTPOPG, &
                         RC, ExtNr=Inst%ExtNr )
       Arr3D => NULL()
       IF ( RC /= HCO_SUCCESS ) THEN
          CALL HCO_ERROR( HcoState%Config%Err, 'HCO_EmisAdd error: EPOP_G', RC )
          RETURN
       ENDIF

    ENDIF

    !----------------------
    ! Manual diagnostics
    !----------------------

    Inst%DiagnName = 'AD53_POPG_SOIL'
    CALL Diagn_Update( HcoState, ExtNr=Inst%ExtNr, &
                       cName=TRIM(Inst%DiagnName), &
                       Array2D=Inst%EMIS_SOIL, RC=RC)
    IF ( RC /= HCO_SUCCESS ) RETURN

    Inst%DiagnName = 'AD53_POPG_LAKE'
    CALL Diagn_Update( HcoState, ExtNr=Inst%ExtNr, &
                       cName=TRIM(Inst%DiagnName), &
                       Array2D=Inst%EMIS_LAKE, RC=RC)
    IF ( RC /= HCO_SUCCESS ) RETURN

    Inst%DiagnName = 'AD53_POPG_LEAF'
    CALL Diagn_Update( HcoState, ExtNr=Inst%ExtNr, &
                       cName=TRIM(Inst%DiagnName), &
                       Array2D=Inst%EMIS_LEAF, RC=RC)
    IF ( RC /= HCO_SUCCESS ) RETURN

    Inst%DiagnName = 'AD53_SOIL2AIR'
    CALL Diagn_Update( HcoState, ExtNr=Inst%ExtNr, &
                       cName=TRIM(Inst%DiagnName), &
                       Array2D=Inst%FLUX_SOIL2AIR, RC=RC)
    IF ( RC /= HCO_SUCCESS ) RETURN

    Inst%DiagnName = 'AD53_AIR2SOIL'
    CALL Diagn_Update( HcoState, ExtNr=Inst%ExtNr, &
                       cName=TRIM(Inst%DiagnName), &
                       Array2D=Inst%FLUX_AIR2SOIL, RC=RC)
    IF ( RC /= HCO_SUCCESS ) RETURN

    Inst%DiagnName = 'AD53_LAKE2AIR'
    CALL Diagn_Update( HcoState, ExtNr=Inst%ExtNr, &
                       cName=TRIM(Inst%DiagnName), &
                       Array2D=Inst%FLUX_LAKE2AIR, RC=RC)
    IF ( RC /= HCO_SUCCESS ) RETURN

    Inst%DiagnName = 'AD53_AIR2LAKE'
    CALL Diagn_Update( HcoState, ExtNr=Inst%ExtNr, &
                       cName=TRIM(Inst%DiagnName), &
                       Array2D=Inst%FLUX_AIR2LAKE, RC=RC)
    IF ( RC /= HCO_SUCCESS ) RETURN

    Inst%DiagnName = 'AD53_LEAF2AIR'
    CALL Diagn_Update( HcoState, ExtNr=Inst%ExtNr, &
                       cName=TRIM(Inst%DiagnName), &
                       Array2D=Inst%FLUX_LEAF2AIR, RC=RC)
    IF ( RC /= HCO_SUCCESS ) RETURN

    Inst%DiagnName = 'AD53_AIR2LEAF'
    CALL Diagn_Update( HcoState, ExtNr=Inst%ExtNr, &
                       cName=TRIM(Inst%DiagnName), &
                       Array2D=Inst%FLUX_AIR2LEAF, RC=RC)
    IF ( RC /= HCO_SUCCESS ) RETURN

    Inst%DiagnName = 'AD53_SOILAIR_FUG'
    CALL Diagn_Update( HcoState, ExtNr=Inst%ExtNr, &
                       cName=TRIM(Inst%DiagnName), &
                       Array2D=Inst%FUG_SOILAIR, RC=RC)
    IF ( RC /= HCO_SUCCESS ) RETURN

    Inst%DiagnName = 'AD53_LAKEAIR_FUG'
    CALL Diagn_Update( HcoState, ExtNr=Inst%ExtNr, &
                       cName=TRIM(Inst%DiagnName), &
                       Array2D=Inst%FUG_LAKEAIR, RC=RC)
    IF ( RC /= HCO_SUCCESS ) RETURN

    Inst%DiagnName = 'AD53_LEAFAIR_FUG'
    CALL Diagn_Update( HcoState, ExtNr=Inst%ExtNr, &
                       cName=TRIM(Inst%DiagnName), &
                       Array2D=Inst%FUG_LEAFAIR, RC=RC)
    IF ( RC /= HCO_SUCCESS ) RETURN

    !=======================================================================
    ! Cleanup & quit
    !=======================================================================

    ! Nullify pointers
    Inst    => NULL()

    ! Return w/ success
    CALL HCO_LEAVE( HcoState%Config%Err,RC )

  END SUBROUTINE HCOX_GC_POPs_Run
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: soilemispop
!
! !DESCRIPTION: Subroutine SOILEMISPOP is the subroutine for secondary
!  POP emissions from soils.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE SOILEMISPOP( POP_SURF, F_OC_SOIL, EPOP_SOIL, &
                              HcoState, ExtState,  Inst )
!
! !INPUT PARAMETERS:
!
      REAL(sp), DIMENSION(:,:), INTENT(IN)  :: POP_SURF   ! POP sfc conc [kg]
      REAL(sp), DIMENSION(:,:), INTENT(IN)  :: F_OC_SOIL  ! Frac C in soil [g/g]
      TYPE(HCO_STATE),          POINTER     :: HcoState   ! Hemco state
      TYPE(Ext_State),          POINTER     :: ExtState   ! Module options
      TYPE(MyInst),             POINTER     :: Inst       ! Instance
!
! !OUTPUT PARAMETERS:
!
      REAL(hp), DIMENSION(:,:), INTENT(OUT) :: EPOP_SOIL  ! POP emissions from
                                                          ! soil [kg/m2/s]
!
! !REMARKS:
!
!
!
! !REVISION HISTORY:
!  21 Aug 2012 - C.L. Friedman - Initial version based on LAND_MERCURY_MOD
!  25 Aug 2015 - M. Sulprizio  - Moved to hcox_gc_POPs_mod.F90
!  02 Oct 2015 - E. Lundgren   - ExtState%POPG is now kg/kg dry air (prev kg)
!  07 Jan 2016 - E. Lundgren   - Update molar gas constant to NIST 2014
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL  :: IS_SNOW_OR_ICE
      LOGICAL  :: IS_LAND_OR_ICE
      INTEGER  :: I, J, L
      REAL(hp) :: POPG
      REAL(hp) :: TK_SURF
      REAL(hp) :: KSA_T, FLUX, F_OC
      REAL(hp) :: SOIL_CONC, DTSRCE, KOA_T
      REAL(hp) :: DIFF, FSOIL, FAIR, DS
      REAL(hp) :: DSA, DAD, DWD, PL
      REAL(hp) :: ZSOIL, ZAIR, TK
      REAL(hp) :: DTCHEM, NEWSOIL, E_KDEG
      REAL(hp) :: FRAC_LAKE, FRAC_SOIL
      REAL(hp) :: FUG_R
      REAL(hp) :: AREA_M2
!
! !DEFINED PARAMETERS:
!
      ! Delta H for POP [kJ/mol]. Delta H is enthalpy of phase transfer
      ! from gas phase to OC. For now we use Delta H for phase transfer
      ! from the gas phase to the pure liquid state.
      ! For PHENANTHRENE:
      ! this is taken as the negative of the Delta H for phase transfer
      ! from the pure liquid state to the gas phase (Schwarzenbach,
      !  Gschwend, Imboden, 2003, pg 200, Table 6.3), or -74000 [J/mol].
      ! For PYRENE:
      ! this is taken as the negative of the Delta H for phase transfer
      ! from the pure liquid state to the gas phase (Schwarzenbach,
      ! Gschwend, Imboden, 2003, pg 200, Table 6.3), or -87000 [J/mol].
      ! For BENZO[a]PYRENE:
      ! this is also taken as the negative of the Delta H for phase transfer
      ! from the pure liquid state to the gas phase (Schwarzenbach,
      ! Gschwend, Imboden, 2003, pg 452, Prob 11.1), or -110,000 [J/mol]
      REAL(hp)            :: DEL_H

      ! R = universal gas constant for adjusting KOA for temp:
      ! 8.3144598 [J/mol/K OR m3*Pa/K/mol]
      REAL(hp), PARAMETER :: R = 8.3144598d0

      ! Molecular weight
      ! For phe, 0.17823 kg/mol
      ! For pyr, 0.20225 kg/mol
      ! For BaP, 0,25231 kg/mol
      REAL(hp)            :: MW

      ! Molecular weight of air
      REAL(hp), PARAMETER :: MWAIR  = 28.97d0 ! g/mol

      ! For PHENANTHRENE:
      ! log KOA_298 = 7.64, or 4.37*10^7 [unitless]
      ! For PYRENE:
      ! log KOA_298 = 8.86, or 7.24*10^8 [unitless]
      ! For BENZO[a]PYRENE:
      ! log KOA_298 = 11.48, or 3.02*10^11 [unitless]
      ! (Ma et al., J. Chem. Eng. Data, 2010, 55:819-825).
      REAL(hp)            :: KOA_298

      ! Set transfer velocity and diffusion coefficient values
      REAL(hp), PARAMETER :: KSA  = 1d0     ![m/h]
      REAL(hp), PARAMETER :: BA   = 0.04d0  ![m2/h]
      REAL(hp), PARAMETER :: BW   = 4d-6    ![m2/h]

      ! Set soil degradation rate
      REAL(hp), PARAMETER :: DEGR = 3.5d-5  ![/h]

      !=================================================================
      ! SOILEMISPOP begins here!
      !=================================================================

!      IF (.NOT. LSECEMISPOP) THEN
!         EPOP_SOIL = 0e+0_hp
!      ELSE

      ! Copy values from ExtState
      DEL_H   = ExtState%POP_DEL_H
      KOA_298 = ExtState%POP_KOA
      MW      = ExtState%POP_XMW

      ! Emission timestep [s]
      DTSRCE  = HcoState%TS_EMIS

      ! Chemistry timestep [h]
      DTCHEM = HcoState%TS_CHEM / 60e+0_hp

      DO J=1, HcoState%NY
      DO I=1, HcoState%NX

         ! Set logicals
         ! Is grid box covered by land/ice or by water? (IS_LAND_OR_ICE)
         ! IS_LAND will return non-ocean boxes but may still contain lakes
         ! If land, is it covered by snow/ice? (IS_SNOW_OR_ICE)
         IS_LAND_OR_ICE = ( ( IS_LAND(I,J,ExtState) ) .OR.  &
                            ( IS_ICE (I,J,ExtState) ) )
         IS_SNOW_OR_ICE = ( ( IS_ICE (I,J,ExtState) ) .OR.  &
                            ( IS_LAND(I,J,ExtState)   .AND. &
                              ExtState%SNOWHGT%Arr%Val(I,J) > 10e+0_hp ) )

         ! Do soils routine only if we are on land that is not covered with
         ! snow or ice
         IF ((IS_LAND_OR_ICE) .AND. .NOT. ( IS_SNOW_OR_ICE )) THEN

            ! Get fraction of grid box covered by lake surface area
            FRAC_LAKE = ExtState%FRLAKE%Arr%Val(I,J)

            ! Get fraction of land remaining
            ! Assume the remaining land is soil and get OC content.
            ! If remaining land is not soil (e.g., desert), there
            ! should be a characteristically low OC content
            ! that will have little capacity to store POPs
            ! ONLY SUBTRACT FRAC LAKE NOW
            FRAC_SOIL = MAX(1e+0_hp - FRAC_LAKE, 0e+0_hp)

            ! Get surface skin temp [K]
            TK_SURF = ExtState%TSKIN%Arr%Val(I,J)

            ! Get air temp [K]
            TK = ExtState%TK%Arr%Val(I,J,1)

            ! Get gas phase air POP concentration at surface in mol/m3
            ! ExtState%POPG is now in units of kg/kg dry air (ewl, 10/2/15)
            POPG = MAX( ExtState%POPG%Arr%Val(I,J,1), SMALLNUM )

! old
!            ! kg / (0.178 kg/mol) /m3 in gridbox
!            POPG = POPG / MW / ExtState%AIRVOL%Arr%Val(I,J,1) ! mol/m3
!            !WRITE(6,*) 'POPG (mol/m3) =', POPG
! new
!            ! (kg trc/kg dry air) / (0.178 kg trc/mol) * (kg dry air/m3)
            POPG = POPG / MW * ExtState%AIRDEN%Arr%Val(I,J,1) ! mol/m3

            ! Grid box surface area [m2]
            AREA_M2   = HcoState%Grid%AREA_M2%Val(I,J)

            ! Get soil concentration in top 5 cm of soil
            ! (following Howsam et al 2000)
            ! From Howsam et al, soil burdens are equal to
            ! 2.6 years deposition for PHE,
            ! 10 years for PYR, and 9.4 years for BaP
            ! Convert to mol/m3
            ! 2.6 yrs * kg deposited to soil in 1 yr * / 0.178 kg/mol
            ! / area grid box (m2) / 0.05 m
            SOIL_CONC = 10e+0_hp * POP_SURF(I,J) / MW / &
                        AREA_M2 / 5e-2_hp ! mol/m3

            ! Get rid of mass due to degradation
            ! Use rate constant for BaP from Mackay and Paterson 1991:
            ! 3.5*10^-5 /h
            ! Calculate exponential factor
            E_KDEG = EXP (-DEGR * DTCHEM)

            ! Adjust conc
            NEWSOIL = SOIL_CONC * E_KDEG

            ! Get foc from GTMM saved files
            F_OC = F_OC_SOIL(I,J)

            ! Define temperature-dependent KOA:
            KOA_T = KOA_298 * EXP((-DEL_H/R) * ( ( 1e+0_hp / TK_SURF ) - &
                    ( 1e+0_hp / 298e+0_hp ) ))

            ! Dimensionless coefficient (mol/m3 soil / mol/m3 air)
            ! KSA = 1.5 (fTOC)*Koa
            KSA_T = 1.5 * F_OC * KOA_T
            KSA_T = MAX( KSA_T, SMALLNUM )

            ! Calculate fugacities from concentrations by dividing by "Z"
            ! values, or the fugacity capacity in mol/m3*Pa following
            ! Mackay and Paterson, 1991

            ! fsoil = Csoil * R * T / KSA [Pa]
            ! where Csoil is in mol/m3, R is in Pa * m3 / mol K
            ! T is in K and KSA is dimensionless
            FSOIL = NEWSOIL * R * TK_SURF / KSA_T

            ! fair = Cair * R * T [Pa]
            ! where Cair is in mol/m3, R is in Pa * m3 / mol K and T is in K
            FAIR = POPG * R * TK

            ! Calculate the fugacity gradient [Pa]
            ! If the gradient is negative, fair is larger and the POP will
            ! diffuse from air to soil
            ! If the gradient is positive, fsoil is larger and POP will
            ! diffuse from soil to air
            DIFF  = FSOIL - FAIR
            FUG_R = FSOIL/FAIR

            ! Calculate "Z" values from fugacities.
            ! Z is the fugacity capacity in mol/m3*Pa. C = Z*f, so Z = C/f
            ZAIR  = POPG / FAIR ! (mol/m3) / (Pa)
            ZSOIL = NEWSOIL / FSOIL ! (mol/m3) / (Pa)

            ! Calculate the "D" value, or the transfer coefficient that
            ! describes the movement of POP between phases (Mackay and
            ! Paterson, 1991). [mol/h*Pa]
            ! The D value for soil-air diffusion is given by
            ! Ds = 1 / (1/Dsa + 1/(Dad + Dwd))
            ! Dsa is the air-side boundary layer diffusion parameter [mol/h*Pa]
            ! Dad is the diffusion parameter between soil particles and
            !  "soil air" [mol/h*Pa]
            ! Dwd is the diffusion parameter between soil particles and
            !  porewater [mol/h*Pa]
            ! Dsa is in series with soil-air and soil-water diffusion,
            !  which are in parallel

            ! Need to define each D value
            ! DSA = kSA * Zair
            ! where kSA is a mass transfer coefficient [m/h],
            ! Zair is the air fugacity capacity [mol/m3*Pa]
            ! ***********
            ! DAD = BA * Zair
            ! where BA is the molecular diffusivity in air [m2/h]
            ! ***********
            ! DWD = BW * Zsoil
            ! where BW is the molecular diffusivity in water [m2/h]
            ! **** PL = the soil diffusion pathlength, set to half the
            ! soil depth (0.025 m)
            DSA = KSA * ZAIR  ! (m/h)  * (mol/m3*Pa) = mol/m2*h*Pa
            DAD = BA* ZAIR    ! (m2/h) * (mol/m3*Pa) = mol/m*h*Pa
            DWD = BW * ZSOIL  ! (m2/h) * (mol/m3*Pa) = mol/m*h*Pa
            PL  = 0.025e+0_hp
            DS  = 1e+0_hp / ( 1e+0_hp/DSA + PL/(DAD+DWD) ) ! mol/(m2*h*Pa) [* m3*Pa/K/mol = m/h/K

            ! Calculate Flux in mol/m2/h
            FLUX = DS * DIFF

            ! Change to units of ng/m2/d for storage
            FLUX = FLUX * 24e+0_hp * MW * 1e+12_hp

            ! Kludge soil emissions from poles for now
            ! Bug somewhere that allows GCAP versions to think some high polar
            ! boxes during some months are land rather than ice - results in
            ! extremely high fluxes
            IF ( HcoState%Grid%YMID%Val(I,J) > 60    .OR. &
                 HcoState%Grid%YMID%Val(I,J) < -60 ) THEN
               FLUX = 0e+0_hp
               DIFF = 0e+0_hp
            ENDIF

            ! Convert to an emission rate in kg/m2/s for returning to
            ! HcoX_GC_POPs_Run
            Inst%EPOP_SOIL(I,J) = MAX(FLUX / 24e+0_hp / 3600e+0_hp / 1e+12_hp, &
                                      0e+0_hp )

            ! Multiply the mass emitted by the fraction of land that is soil
            Inst%EPOP_SOIL(I,J) = FRAC_SOIL * Inst%EPOP_SOIL(I,J)

            ! If the flux is positive, then the direction will be from the
            ! soil to the air.
            ! Store this in a diagnostic array.
            ! If the flux is zero or negative, store it in a separate array.
            IF ( FLUX > 0e+0_hp ) THEN

               ! Store total mass emitted from soil [kg] for ND53 diagnostic.
               ! (We don't care about the mass in the other direction right
               ! now)
               Inst%EMIS_SOIL(I,J)     = Inst%EPOP_SOIL(I,J) * AREA_M2 * DTSRCE

               ! Store positive flux
               Inst%FLUX_SOIL2AIR(I,J) = FLUX

               ! Make sure negative flux diagnostic has nothing added to it
               Inst%FLUX_AIR2SOIL(I,J) = 0e+0_hp

               ! Store the soil/air fugacity ratio
               Inst%FUG_SOILAIR(I,J)   = FUG_R

            ELSE IF ( FLUX <= 0e+0_hp ) THEN

               ! Store the negative flux
               Inst%FLUX_AIR2SOIL(I,J) = FLUX

               ! Add nothing to positive flux or mass diagnostics
               Inst%EMIS_SOIL(I,J)     = 0e+0_hp
               Inst%FLUX_SOIL2AIR(I,J) = 0e+0_hp

               ! Continue to store the fugacity ratio
               Inst%FUG_SOILAIR(I,J)   = FUG_R

            ENDIF

         ELSE

            ! We are not on land or the land is covered with ice or snow
            FLUX                    = 0e+0_hp
            Inst%EPOP_SOIL(I,J)     = 0e+0_hp
            Inst%EMIS_SOIL(I,J)     = 0e+0_hp
            Inst%FLUX_SOIL2AIR(I,J) = 0e+0_hp
            Inst%FLUX_AIR2SOIL(I,J) = 0e+0_hp
            Inst%FUG_SOILAIR(I,J)   = 0e+0_hp

         ENDIF

      ENDDO
      ENDDO

      END SUBROUTINE SOILEMISPOP
!EOC

!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: lakeemispop
!
! !DESCRIPTION: Subroutine LAKEEMISPOP is the subroutine for secondary
!  POP emissions from lakes.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE LAKEEMISPOP( POP_SURF, EPOP_LAKE, HcoState, ExtState, Inst )
!
! !USES:
!
!
! !INPUT PARAMETERS:
!
      REAL(sp), DIMENSION(:,:), INTENT(IN)  :: POP_SURF   ! POP sfc conc [kg]
      TYPE(HCO_STATE),          POINTER     :: HcoState   ! Hemco state
      TYPE(Ext_State),          POINTER     :: ExtState   ! Module options
      TYPE(MyInst),             POINTER     :: Inst       ! Instance
!
! !OUTPUT PARAMETERS:
!
      REAL(hp), DIMENSION(:,:), INTENT(OUT) :: EPOP_LAKE  ! POP emissions from
                                                          ! lakes [kg/m2/s]
!
! !REMARKS:
!
! !REVISION HISTORY:
!  21 Aug 2012 - C.L. Friedman - Initial version based on LAND_MERCURY_MOD
!  25 Aug 2015 - M. Sulprizio  - Moved to hcox_gc_POPs_mod.F90
!  02 Oct 2015 - E. Lundgren   - ExtState%POPG is now kg/kg dry air (prev kg)
!  07 Jan 2016 - E. Lundgren   - Update molar gas constant to NIST 2014

!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!

      INTEGER  :: I, J, L
      REAL(hp) :: TK_SURF
      REAL(hp) :: KAW_T, FLUX, KOL_T
      REAL(hp) :: DTSRCE
      REAL(hp) :: KA_H2O, KA_POP, KW_CO2, KW_POP
      REAL(hp) :: DA_H2O, DA_POP, DW_CO2, DW_POP
      REAL(hp) :: SCH_CO2, SCH_POP
      REAL(hp) :: TK, PRESS
      REAL(hp) :: C_DISS
      REAL(hp) :: POPG, U10M, ALPHA
      REAL(hp) :: FRAC_LAKE, VISC_H2O
      REAL(hp) :: SFCWINDSQR
      REAL(hp) :: AREA_M2
      LOGICAL  :: IS_SNOW_OR_ICE, IS_LAND_OR_ICE
!
! !DEFINED PARAMETERS:
!
      ! Delta H for POP:
      ! For PHENANTHRENE:
      ! this is the Delta H for phase transfer
      ! from air to water (Schwarzenbach,
      !  Gschwend, Imboden, 2003, pg 200, Table 6.3), or 47 [kJ/mol].
      ! For PYRENE: 43000 [kJ/mol]
      REAL(hp)            :: DEL_HW

      ! R = universal gas constant for adjusting KOA for temp:
      ! 8.3144598d-3 [kJ/mol/K]
      REAL(hp), PARAMETER :: R = 8.3144598d-3

      ! Molecular weight
      ! For phe, 0.17823 kg/mol
      REAL(hp)            :: MWPOP

      ! Molecular weight of air
      REAL(hp), PARAMETER :: MWAIR  = 28.97d0 ! g/mol

      ! Molecular weight of water
      REAL(hp), PARAMETER :: MWH2O  = 18.1d0 ! g/mol

      ! Molar volumes calculated following Abraham and McGowan 1987 as
      ! summarized by Schwarzenbach et al. 2003.
      ! Each element is assigned a characteristic atomic volume, and an atomic
      ! volume of 6.56 cm3/mol is subtracted for each bond, no matter whether
      ! single, double, or triple
      ! C = 16.35,  H = 8.71,  O = 12.43, N = 14.39, P = 24.87, F = 10.48
      ! Br = 26.21, I = 34.53, S = 22.91, Si = 26.83

      ! Molar volume of water
      ! 2*(8.71) + 12.43 - 2*(6.56) = 16.73
      REAL(hp), PARAMETER :: V_H2O  = 16.73d0 ! cm3/mol

      ! Molar volume of CO2
      ! 16.35 + 2*(12.43) - 2*(6.56) = 28.1
      REAL(hp), PARAMETER :: V_CO2   = 28.1d0  ! cm3/mol

      ! Molar volume of POP
      ! For PHE (C16H10):
      ! 16*(16.35) + 10*(8.71) - 29*(6.56) =
      REAL(hp), PARAMETER :: V_POP  = 538.94d0 ! cm3/mol

      ! Molar volume of air - average of gases in air
      REAL(hp), PARAMETER :: V_AIR  = 20.1d0     ! cm3/mol

      ! For PHENANTHRENE:
      ! log KAW_298 = -2.76, or 1.74*10-3 [unitless]
      ! For PYRENE:
      ! log KAW_298 = -3.27, or 5.37*10-4 [unitless]
      REAL(hp)            :: KAW_298

      ! Set the kinematic viscosity of freshwater at 20C
!      REAL(hp), PARAMETER :: VISC_H2O  ! = 1d0    ![cm2/s]

      ! Set aqueous degradation rate
!      REAL(hp), PARAMETER :: DEGR      != 3.5d-5  ![/h]

      !=================================================================
      ! LAKEEMISPOP begins here!
      !=================================================================

!      IF (.NOT. LSECEMISPOP) THEN
!         Inst%EPOP_LAKE = 0e+0_hp
!      ELSE

      ! Do lake emissions routine:

      ! Copy values from ExtState
      DEL_HW  = ExtState%POP_DEL_HW
      MWPOP   = ExtState%POP_XMW
      KAW_298 = ExtState%POP_HSTAR

      ! Emission timestep [s]
      DTSRCE  = HcoState%TS_EMIS

      DO J=1, HcoState%NY
      DO I=1, HcoState%NX

         ! Set logicals
         ! Is grid box covered by land/ice or by water? (IS_LAND_OR_ICE)
         ! IS_LAND will return non-ocean boxes but may still contain lakes
         ! If land, is it covered by snow/ice? (IS_SNOW_OR_ICE)
         IS_LAND_OR_ICE = ( ( IS_LAND(I,J,ExtState) ) .OR.  &
                            ( IS_ICE (I,J,ExtState) ) )
         IS_SNOW_OR_ICE = ( ( IS_ICE (I,J,ExtState) ) .OR.  &
                            ( IS_LAND(I,J,ExtState)   .AND. &
                              ExtState%SNOWHGT%Arr%Val(I,J) > 10e+0_hp ) )

         ! Do soils routine only if we are on land that is not covered with
         ! snow or ice
         IF ((IS_LAND_OR_ICE) .AND. .NOT. ( IS_SNOW_OR_ICE )) THEN

            ! Get fraction of grid box covered by lake surface area
            FRAC_LAKE = ExtState%FRLAKE%Arr%Val(I,J)

            IF ( FRAC_LAKE > 0e+0_hp ) THEN

               ! Get surface skin temp [K]
               TK_SURF = ExtState%TSKIN%Arr%Val(I,J)

               ! Get air temp [K]
               TK = ExtState%TK%Arr%Val(I,J,1)

               ! Get surface pressure at end of dynamic time step [hPa]
               PRESS = ExtState%PSC2_WET%Arr%Val(I,J)

               ! Convert to units of atm
               PRESS = PRESS / 1013.25e+0_hp

               ! Get gas phase air POP concentration at surface in mol/m3
               ! ExtState%POPG is now in units of kg/kg dry air (ewl, 10/2/15)
               POPG = MAX( ExtState%POPG%Arr%Val(I,J,1), SMALLNUM )

! old
!               ! kg / (kg/mol) /m3 in gridbox
!               POPG = POPG / MWPOP / ExtState%AIRVOL%Arr%Val(I,J,1) ! mol/m3
! new
               ! (kg trc/kg dry air) / (kg trc/mol) * (kg dry air/m3)
               POPG = POPG / MWPOP * ExtState%AIRDEN%Arr%Val(I,J,1) ! mol/m3


               ! Grid box surface area [m2]
               AREA_M2   = HcoState%Grid%AREA_M2%Val(I,J)

               ! Get the dissolved POP concentration at lake surface in mol/m3
               ! Distribute the total deposited mass to a volume that best
               ! matches observed dissolved concentrations
               ! Future versions should consider aqueous particle concentrations
               ! and sinking rates and photolytic/microbial degradation

               ! Start with 1 m - scale by 100
               C_DISS = POP_SURF(I,J) / MWPOP / AREA_M2 / 100e+0_hp

               ! Wind speed at 10m altitude [m/s]
               SFCWINDSQR = ExtState%U10M%Arr%Val(I,J)**2 &
                          + ExtState%V10M%Arr%Val(I,J)**2
               U10M       = SQRT( SFCWINDSQR )

               ! Need to calculate water-side and air-side mass transfer
               ! coefficients
               ! Start with air-side
               ! Relate air-side MTC of POP to that of H2O
               ! First, calculate air-side MTC for water following
               ! Schwarzenbach, Gschwend, Imboden 2003
               KA_H2O = 0.2e+0_hp * U10M + 0.3e+0_hp   ! cm/s

               ! Relate air-side MTC for water to that of POP via diffusivities
               ! following Schwarzenbach et al 2003
               ! Calculate temperature-dependent diffusivities in air first
               ! folling Fuller et al. 1966 (summarized by Schwarzenbach
               ! et al. 2003)
               DA_POP = 1e-3_hp * TK**1.75e+0_hp * ( (1e+0_hp/MWAIR) + &
                        (1e+0_hp / (MWPOP*1e+3_hp) ) )**0.5e+0_hp /    &
                        ( PRESS * ( V_AIR**(1e+0_hp/3e+0_hp) +         &
                        V_POP**(1e+0_hp/3e+0_hp))**2e+0_hp )  ! cm2/s

               DA_H2O = 1e-3_hp * TK**1.75e+0_hp * ( (1e+0_hp/MWAIR) + &
                        (1e+0_hp / MWH2O ) )**0.5e+0_hp /              &
                        ( PRESS * ( V_AIR**(1e+0_hp/3e+0_hp) +         &
                        V_H2O**(1e+0_hp/3e+0_hp))**2e+0_hp )  ! cm2/s

               ! Relate POP and H2O air-side MTCs
               KA_POP = KA_H2O * ( DA_POP / DA_H2O )**( 0.67e+0_hp ) ! cm/s

               ! Now calculate water-side MTCs
               ! Start with calculating the water side MTC of CO2
               ! This depends on wind speed (Schwarzenbach et al 2003)
               ! Three different scenarios for the water surface under
               ! different wind speeds are considered:
               ! Smooth Surface Regime (SSR): u10 <= 4.2 m/s
               ! Rough Surface Regime (RSR):  4.2 m/s < u10 <= 13 m/s
               ! Breaking Wave Regime (BWR):  u10 > 13 m/s
               ! Alpha, the exponent in the relationship for the water side MTC,
               ! also depends on the wind speed. Set this as well
               IF ( U10M <= 4.2e+0_hp ) THEN
                  KW_CO2 = 0.65e-3_hp       ! cm/s
                  ALPHA  = 0.67e+0_hp
               ELSE IF ( U10M > 4.2e+0_hp .AND. U10M <= 13e+0_hp ) THEN
                  KW_CO2 = ( 0.79e+0_hp * U10M - 2.68e+0_hp ) * 1e-3_hp
                  ALPHA  = 0.5e+0_hp
               ELSE IF (U10M > 13e+0_hp) THEN
                  KW_CO2 = ( 1.64e+0_hp * U10M - 13.69e+0_hp ) * 1e-3_hp
                  ALPHA  = 0.50e+0_hp
               ENDIF

               ! Get the temperature-dependent kinematic viscosity of water
               IF (TK_SURF <= 273.15 ) THEN
                  VISC_H2O = 1.787e-2_hp      ! [cm2/s]
               ELSE IF (TK_SURF > 273.15 .AND. TK_SURF <= 278.15 ) THEN
                  VISC_H2O = 1.518e-2_hp
               ELSE IF (TK_SURF > 258.15 .AND. TK_SURF <= 283.15 ) THEN
                  VISC_H2O = 1.307e-2_hp
               ELSE IF (TK_SURF > 283.15 .AND. TK_SURF <= 287.15 ) THEN
                  VISC_H2O = 1.139e-2_hp
               ELSE IF (TK_SURF > 287.15 .AND. TK_SURF <= 293.15 ) THEN
                  VISC_H2O = 1.002e-2_hp
               ELSE IF (TK_SURF > 293.15 .AND. TK_SURF <= 298.15 ) THEN
                  VISC_H2O = 0.89e-2_hp
               ELSE IF (TK_SURF > 298.15 ) THEN
                  VISC_H2O = 0.797e-2_hp
               ENDIF

               ! Calculate the diffusivites of CO2 and POP in water
               DW_CO2 = ( 13.26 * 1e-5_hp ) /               &
                        (( VISC_H2O*1e+2_hp )**1.14e+0_hp * &
                        (V_CO2)**0.589e+0_hp )   ! [cm2/s]

               DW_POP = ( 13.26 * 1e-5_hp ) /               &
                        (( VISC_H2O*1e+2_hp )**1.14e+0_hp * &
                        (V_POP)**0.589e+0_hp )   ! [cm2/s]

               ! Calculate the Schmidt numbers for CO2 and POP
               SCH_CO2 = VISC_H2O / DW_CO2    ! [unitless]
               SCH_POP = VISC_H2O / DW_POP    ! [unitless]

               ! Calculate the water-side MTC for POP
               KW_POP = KW_CO2 * ( SCH_POP / SCH_CO2 ) ** (-ALPHA) ! [cm/s]

               ! Calculate the temperature-dependent dimensionless Henry's Law
               ! constant
               KAW_T = KAW_298 * EXP((-DEL_HW/R) * ((1e+0_hp/TK_SURF) - &
                        (1e+0_hp/298e+0_hp)))       ! [unitless]

               ! Now calculate the overall air-water MTC
               KOL_T = 1e+0_hp / ( 1e+0_hp/KW_POP + &
                       1e+0_hp / (KA_POP*KAW_T) ) ! [cm/s]

               ! Calculate Flux in ng/m2/day !
               FLUX = KOL_T * 3600e+0_hp * 24e+0_hp * ( C_DISS - &
                      POPG/KAW_T ) * MWPOP * 1e+12_hp / 100e+0_hp

               ! Convert to an emission rate in kg/m2/s for returning to
               ! HcoX_GC_POPs_Run
               ! Only return it if it's positive
               EPOP_LAKE(I,J) = MAX(FLUX / 24e+0_hp / 3600e+0_hp / 1e+12_hp, &
                                    0e+0_hp )

               ! Multiply the mass emitted by the fraction of land that is water
               EPOP_LAKE(I,J) = FRAC_LAKE * EPOP_LAKE(I,J)

               ! If the flux is positive, then the direction will be from the
               ! soil to the air.
               ! Store this in a diagnostic array.
               ! If the flux is zero or negative, store it in a separate array.
               IF ( FLUX > 0e+0_hp ) THEN

                  ! Store total mass emitted from soil [kg] in ND53 diagnostic.
                  ! (We don't care about the mass in the other direction right
                  ! now)
                  Inst%EMIS_LAKE(I,J)     = EPOP_LAKE(I,J) * AREA_M2 * DTSRCE

                  ! Store positive flux
                  Inst%FLUX_LAKE2AIR(I,J) = + FLUX

                  ! Make sure negative flux diagnostic has nothing added to it
                  Inst%FLUX_AIR2LAKE(I,J) = 0e+0_hp

                     ! Store the soil/air fugacity ratio
                  Inst%FUG_LAKEAIR(I,J)   = C_DISS / (POPG/KAW_T)

               ELSE IF ( FLUX <= 0e+0_hp ) THEN

                  ! Store the negative flux
                  Inst%FLUX_AIR2LAKE(I,J) = FLUX

                  ! Add nothing to positive flux or mass diagnostics
                  Inst%EMIS_LAKE(I,J)     = 0e+0_hp
                  Inst%FLUX_LAKE2AIR(I,J) = 0e+0_hp

                  ! Continue to store the fugacity ratio
                  Inst%FUG_LAKEAIR(I,J)   = C_DISS / (POPG/KAW_T)

               ENDIF

            ENDIF

         ELSE

            ! We are not on land or the land is covered with ice or snow
            ! or we are land but there is no water
            FLUX                    = 0e+0_hp
            Inst%EPOP_LAKE(I,J)     = 0e+0_hp
            Inst%EMIS_LAKE(I,J)     = 0e+0_hp
            Inst%FLUX_LAKE2AIR(I,J) = 0e+0_hp
            Inst%FLUX_AIR2LAKE(I,J) = 0e+0_hp
            Inst%FUG_LAKEAIR(I,J)   = 0e+0_hp

         ENDIF

      ENDDO
      ENDDO

      END SUBROUTINE LAKEEMISPOP
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: vegemispop
!
! !DESCRIPTION: Subroutine VEGEMISPOP is the subroutine for secondary
!  POP emissions from soils.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE VEGEMISPOP( POP_SURF,  EPOP_VEG, HcoState, ExtState, Inst )
!
! !USES:
!
!
! !INPUT PARAMETERS:
!
      REAL(sp), DIMENSION(:,:), INTENT(IN)  :: POP_SURF   ! POP sfc conc [kg]
      TYPE(HCO_State),          POINTER     :: HcoState   ! Hemco state
      TYPE(Ext_State),          POINTER     :: ExtState   ! Module options
      TYPE(MyInst),             POINTER     :: Inst       ! Instance
!
! !OUTPUT PARAMETERS:
!
      REAL(hp), DIMENSION(:,:), INTENT(OUT) :: EPOP_VEG   ! POP emissions from
                                                          ! leaves [kg/m2/s]
!
! !REMARKS:
!
!
!
! !REVISION HISTORY:
!  21 Aug 2012 - C.L. Friedman - Initial version based on LAND_MERCURY_MOD
!  25 Aug 2015 - M. Sulprizio  - Moved to hcox_gc_POPs_mod.F90
!  02 Oct 2015 - E. Lundgren   - ExtState%POPG is now kg/kg dry air (prev kg)
!  07 Jan 2016 - E. Lundgren   - Update molar gas constant to NIST 2014
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!

      INTEGER  :: I, J, L
      REAL(hp) :: POPG, POPG_GL
      REAL(hp) :: FRAC_SNOWFREE_LAND, TK_SURF
      REAL(hp) :: KLA_T, FLUX, K_MT
      REAL(hp) :: LEAF_CONC, DTSRCE, KOA_T
      REAL(hp) :: DIFF, FLEAF, FAIR, DS
      REAL(hp) :: DAB_F, DC, DAL, PC, UC
      REAL(hp) :: ZLEAF, ZAIR, TK
      REAL(hp) :: NEWLEAF
      REAL(hp) :: LAI, KAW_T, KOW_T
      REAL(hp) :: FUG_R, NEW_LEAF, DLA
      REAL(hp) :: AREA_M2
      LOGICAL  :: IS_SNOW_OR_ICE, IS_LAND_OR_ICE
      LOGICAL  :: IS_SNOWFREE_LAND

!
! !DEFINED PARAMETERS:
!
      ! Delta H for POP [kJ/mol]. Delta H is enthalpy of phase transfer
      ! from gas phase to OC. For now we use Delta H for phase transfer
      ! from the gas phase to the pure liquid state.
      ! For PHENANTHRENE:
      ! this is taken as the negative of the Delta H for phase transfer
      ! from the pure liquid state to the gas phase (Schwarzenbach,
      !  Gschwend, Imboden, 2003, pg 200, Table 6.3), or -74000 [J/mol].
      ! For PYRENE:
      ! this is taken as the negative of the Delta H for phase transfer
      ! from the pure liquid state to the gas phase (Schwarzenbach,
      ! Gschwend, Imboden, 2003, pg 200, Table 6.3), or -87000 [J/mol].
      ! For BENZO[a]PYRENE:
      ! this is also taken as the negative of the Delta H for phase transfer
      ! from the pure liquid state to the gas phase (Schwarzenbach,
      ! Gschwend, Imboden, 2003, pg 452, Prob 11.1), or -110,000 [J/mol]
      REAL(hp)             :: DEL_H

      ! Delta H for POP [kJ/mol].
      ! For PHENANTHRENE:
      ! this is the Delta H for phase transfer
      ! from air to water (Schwarzenbach,
      !  Gschwend, Imboden, 2003, pg 200, Table 6.3), or 47000 [J/mol].
      ! For PYRENE: 43000 [J/mol]
      REAL(hp)            :: DEL_HW

      ! R = universal gas constant for adjusting KOA for temp:
      ! 8.314 [J/mol/K OR m3*Pa/K/mol]
      REAL(hp), PARAMETER :: R = 8.3144598e+0_hp

      ! Molecular weight
      ! For phe, 0.17823 kg/mol
      REAL(hp)            :: MW

      ! Molecular weight of air
      REAL(hp), PARAMETER :: MWAIR  = 28.97d0 ! g/mol

      ! For PHENANTHRENE:
      ! log KOA_298 = 7.64, or 4.37*10^7 [unitless]
      ! For PYRENE:
      ! log KOA_298 = 8.86, or 7.24*10^8 [unitless]
      ! For BENZO[a]PYRENE:
      ! log KOA_298 = 11.48, or 3.02*10^11 [unitless]
      ! (Ma et al., J. Chem. Eng. Data, 2010, 55:819-825).
      REAL(hp)            :: KOA_298

      ! For PHENANTHRENE:
      ! log KAW_298 = -2.76, or 1.74*10-3 [unitless]
      ! For PYRENE:
      ! log KAW_298 = -3.27, or 5.37*10-4 [unitless]
      REAL(hp)            :: KAW_298

      ! Set volume fractions of octanol and water in surface and reservoir
      ! leaf compartments [unitless]
      REAL(hp), PARAMETER :: OCT_SURF = 0.8e+0_hp
      REAL(hp), PARAMETER :: OCT_RES  = 0.02e+0_hp
      REAL(hp), PARAMETER :: H2O_RES  = 0.7e+0_hp

      ! Set thickness of different leaf compartments. Volumes calculated by
      ! multiplying thicknesses by leaf area index
      REAL(hp), PARAMETER :: SURF_THICK = 2e-6_hp   ! m
      REAL(hp), PARAMETER :: RES_THICK  = 250e-6_hp ! m

      ! Set transfer velocity and diffusion coefficient values
      REAL(hp), PARAMETER :: UAB_F = 9e+0_hp   ![m/h]

      ! Set soil degradation rate
      REAL(hp), PARAMETER :: DEGR  = 3.5e-5_hp ![/h]

      !=================================================================
      ! VEGEMISPOP begins here!
      !=================================================================

      ! Copy values from ExtState
      DEL_H   = ExtState%POP_DEL_H
      KOA_298 = ExtState%POP_KOA
      DEL_HW  = ExtState%POP_DEL_Hw
      KAW_298 = ExtState%POP_HSTAR
      MW      = ExtState%POP_XMW

      ! Emission timestep [s]
      DTSRCE  = HcoState%TS_EMIS

      DO J=1, HcoState%NY
      DO I=1, HcoState%NX

         ! Set logicals
         ! Is grid box covered by land/ice or by water? (IS_LAND_OR_ICE)
         ! IS_LAND will return non-ocean boxes but may still contain lakes
         ! If land, is it covered by snow/ice? (IS_SNOW_OR_ICE)
         IS_LAND_OR_ICE = ( (IS_LAND(I,J,ExtState)) .OR.  &
                            (IS_ICE (I,J,ExtState)) )
         IS_SNOW_OR_ICE = ( (IS_ICE (I,J,ExtState)) .OR.  &
                            (IS_LAND(I,J,ExtState)  .AND. &
                             ExtState%SNOWHGT%Arr%Val(I,J) > 10e+0_hp ) )

         ! Do soils routine only if we are on land that is not covered with
         ! snow or ice
         IF ((IS_LAND_OR_ICE) .AND. .NOT. ( IS_SNOW_OR_ICE )) THEN

            ! Get fraction of grid box covered by leaf surface area
            ! Do not consider different vegetation types for now
            LAI = ExtState%LAI%Arr%Val(I,J)

            IF ( LAI > 0e+0_hp ) THEN

               ! Get surface skin temp [K]
               TK_SURF = ExtState%TSKIN%Arr%Val(I,J)

               ! Get air temp [K]
               TK = ExtState%TK%Arr%Val(I,J,1)

               ! Get gas phase air POP concentration at surface in mol/m3
               ! ExtState%POPG is now in units of kg/kg dry air (ewl, 10/2/15)
               POPG = MAX( ExtState%POPG%Arr%Val(I,J,1), SMALLNUM )

! old
!               ! kg / (0.178 kg/mol) /m3 in gridbox
!               POPG = POPG / MW / ExtState%AIRVOL%Arr%Val(I,J,1) ! mol/m3
! new
               ! (kg trc/kg dry air) / (kg trc/mol) * (kg dry air/m3)
               POPG = POPG / MW * ExtState%AIRDEN%Arr%Val(I,J,1) ! mol/m3


               ! Grid box surface area [m2]
               AREA_M2   = HcoState%Grid%AREA_M2%Val(I,J)

               ! Only consider partitioning into leaf surface (not reservoir)
               ! for now following Mackay et al 2006 Environ Sci & Pollut Res
               ! Include reservoir when land-atm models become dynamic

               ! Assume that all leaf surfaces contain an average lipid content
               ! of 80% (Mackay et al 2006)

               ! Get leaf concentration
               ! Convert to mol/m3
               ! kg deposited to leaf in 1 yr * / 0.178 kg/mol
               ! / area grid box (m2) / surface thickness m
               LEAF_CONC = POP_SURF(I,J) / MW / AREA_M2 / SURF_THICK ! mol/m3

               ! Check concentration in leaves by assuming a density similar
               ! to water (1 g/cm3)

               ! No degradation/metabolism for now. Just scale leaf
               ! concentrations to match flux observations
               NEWLEAF = LEAF_CONC/1e+4_hp  !SCALING FACTOR

               ! Define temperature-dependent KOA:
               KOA_T = KOA_298 * EXP((-DEL_H/R) * ((1e+0_hp/TK_SURF) - &
                       (1e+0_hp/298e+0_hp)))

               ! Calculate the temperature-dependent dimensionless Henry's Law
               ! constant
               KAW_T = KAW_298 * EXP((-DEL_HW/R) * ((1e+0_hp/TK_SURF) - &
                       (1e+0_hp/298e+0_hp)))       ! [unitless]

               ! Estimate the temperature-dependent dimensionless octanol-water
               ! constant
               KOW_T = KOA_T * KAW_T

               ! Define dimensionless leaf surface-air partition coefficient
               ! (mol/m3 leaf / mol/m3 air)
               ! KLA = foct_surf * Koa
               KLA_T = OCT_SURF * KOA_T
!               KLA_T = MAX( KLA_T, SMALLNUM )

               ! Calculate fugacities from concentrations by dividing by "Z"
               ! values, or the fugacity capacity in mol/m3*Pa following Mackay
               ! and Paterson, 1991

               ! fleaf = Cleaf * R * T / KSA [Pa]
               ! where Cleaf is in mol/m3, R is in Pa * m3 / mol K
               ! T is in K and KLA is dimensionless
               FLEAF = NEWLEAF * R * TK_SURF / KLA_T

               ! fair = Cair * R * T [Pa]
               ! where Cair is in mol/m3, R is in Pa * m3 / mol K and T is in K
               FAIR = POPG * R * TK

               ! Calculate the fugacity gradient [Pa]
               ! If the gradient is negative, fair is larger and the POP will
               ! diffuse from air to soil
               ! If the gradient is positive, fsoil is larger and POP will
               ! diffuse from soil to air
               DIFF = FLEAF - FAIR
               FUG_R = FLEAF/FAIR

               ! Calculate "Z" values from fugacities.
               ! Z is the fugacity capacity in mol/m3*Pa. C = Z*f, so Z = C/f
               ZAIR  = POPG / FAIR ! (mol/m3) / (Pa)
               ZLEAF = NEWLEAF / FLEAF ! (mol/m3) / (Pa)

               ! Calculate the "D" value, or the transfer coefficient that
               ! describes the movement of POP between phases (Mackay and
               ! Paterson, 1991, Cousins and Mackay 2000, internal report).
               ! [mol/h*Pa]
               ! The D value for leaf surface-air gas diffusion is given by
               ! Dla = 1 / (1/Dc + 1/(Dab-f)) [mol/(Pa*h)]
               ! where Dab-f is the boundary layer diffusion  [mol/h*Pa]
               ! given by Dab-f = As * L * Uab-f * Za
               ! where As is the area of the land surface [m2], L is the leaf
               ! area index [m2/m2],
               ! Uab-f is a mass transfer coefficient for surface-air boundary
               ! layer diffusion [m/h],
               ! and Za is the fugacity capacity of the air [mol/(m3*Pa)]
               ! Dc is the cuticle diffusion, given by
               ! Dc = As * L * Uc * Zf
               ! where As and L are as above, Uc is the cuticle mass transfer
               ! coefficient [m/h],
               ! and Zf is the fugacity capacity of the leaf surface
               ! (mol/(m3*Pa))

               ! Uc is given by
               ! Uc = 3600 * Pc * 1/Kaw
               ! where Pc is the cuticle permeance (m/s) and Kaw is the
               ! dimensionless air-water partition coefficient.
               ! Pc is given by
               ! Log Pc = ((0.704 * log Kow - 11.2) +
               !          (-3.47 - 2.79 * logMW + 0.970 log Kow)) / 2
               ! (an average of two equations)

               ! Need to define each D value
               ! DAB_F:
               !  m/h * mol/(m3*Pa)  =  (mol/h*Pa*m2)
               DAB_F = UAB_F * ZAIR  !  mol/(h*Pa*m2)

               ! Calculate PC and then Uc in order to calculate Dc
               ! PC, UC are calculated according to Cousins and Mackay,
               ! Chemosphere, 2001, Table 2
               PC = 10** (( 0.704e+0_hp * LOG (KOW_T) - 11.2e+0_hp ) +    &
                    ( -3.47e+0_hp -2.79e+0_hp* LOG(MW*1000d0) + 0.97e+0_hp &
                    * LOG(KOW_T))/2e+0_hp) ![m/s]

               UC = 3600e+0_hp * PC * 1e+0_hp/KAW_T ! [m/h]

               DC = UC * ZLEAF    ! mol/(h*Pa*m2)

               ! Now calculate overall transfer  ! mol/(h*Pa*m2)
               DLA = 1e+0_hp / (1e+0_hp/DC + 1e+0_hp/DAB_F)

               ! Calculate Flux in mol/h/m2
               FLUX = DLA * DIFF

               ! Change to units of ng/m2/d for storage
               FLUX = FLUX * 24e+0_hp * MW * 1e+12_hp

               ! Convert to an emission rate in kg/m2/s for returning to
               ! HcoX_GC_POPs_Run
               ! Only want to add rates that are positive
               Inst%EPOP_VEG(I,J) = MAX(FLUX * LAI / 24e+0_hp / 3600e+0_hp / &
                                   1e+12_hp, 0e+0_hp)

               ! If the flux is positive, then the direction will be from the
               ! soil to the air.
               ! Store this in a diagnostic array.
               ! If the flux is zero or negative, store it in a separate array.
               IF ( FLUX > 0e+0_hp ) THEN

                  ! Store total mass emitted from soil [kg] in ND53 diagnostic.
                  ! (We don't care about the mass in the other direction right
                  ! now)
                  Inst%EMIS_LEAF(I,J)     = EPOP_VEG(I,J) * AREA_M2 * DTSRCE

                  ! Store positive flux
                  Inst%FLUX_LEAF2AIR(I,J) = FLUX

                  ! Make sure negative flux diagnostic has nothing added to it
                  Inst%FLUX_AIR2LEAF(I,J) = 0e+0_hp

                  ! Store the soil/air fugacity ratio
                  Inst%FUG_LEAFAIR(I,J) = FUG_R

               ELSE IF ( FLUX <= 0e+0_hp ) THEN

                  ! Store the negative flux
                  Inst%FLUX_AIR2LEAF(I,J) = FLUX

                  ! Add nothing to positive flux or mass diagnostics
                  Inst%EMIS_LEAF(I,J)     = 0e+0_hp
                  Inst%FLUX_LEAF2AIR(I,J) = 0e+0_hp

                  ! Continue to store the fugacity ratio
                  Inst%FUG_LEAFAIR(I,J)   = FUG_R

               ENDIF

            ELSE

               ! We are not on land or the land is covered with ice or snow
               FLUX                    = 0e+0_hp
               Inst%EPOP_VEG(I,J)      = 0e+0_hp
               Inst%EMIS_LEAF(I,J)     = 0e+0_hp
               Inst%FLUX_LEAF2AIR(I,J) = 0e+0_hp
               Inst%FLUX_AIR2LEAF(I,J) = 0e+0_hp
               Inst%FUG_LEAFAIR(I,J)   = 0e+0_hp

            ENDIF

         ENDIF

      ENDDO
      ENDDO

      END SUBROUTINE VEGEMISPOP
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: is_land
!
! !DESCRIPTION: Function IS\_LAND returns TRUE if surface grid box (I,J) is
!  a land box.
!\\
!\\
! !INTERFACE:
!
      FUNCTION IS_LAND( I, J, ExtState ) RESULT ( LAND )
!
! !USES:
!
!
! !INPUT PARAMETERS:
!
      INTEGER,         INTENT(IN) :: I           ! Longitude index of grid box
      INTEGER,         INTENT(IN) :: J           ! Latitude  index of grid box
      TYPE(Ext_State), POINTER    :: ExtState    ! Module options
!
! !RETURN VALUE:
!
      LOGICAL                     :: LAND        ! =T if it is a land box
!
! !REVISION HISTORY:
!  26 Jun 2000 - R. Yantosca - Initial version
!  (1 ) Now use ALBEDO field to determine land or land ice boxes for GEOS-3.
!        (bmy, 4/4/01)
!  (2 ) For 4x5 data, regridded albedo field can cause small inaccuracies
!        near the poles (bmy, 4/4/01)
!  (3 ) Add references to CMN_SIZE and CMN, so that we can use the JYEAR
!        variable to get the current year.  Also, for 1998, we need to compute
!        if is a land box or not from the surface albedo, since for this
!        year the LWI/SURFTYPE field is not given.  For other years than 1998,
!        we use LWI(I,J) < 50 as our land box criterion.  Deleted obsolete
!        code and updated comments.(mje, bmy, 1/9/02)
!  (4 ) Deleted GEOS-2 #ifdef statement.  GEOS-2 met fields never really
!        materialized, we use GEOS-3 instead. (bmy, 9/18/02)
!  (5 ) Now uses function GET_YEAR from "time_mod.f".  Removed reference
!        to CMN header file. (bmy, 3/11/03)
!  (6 ) Added code to determine land boxes for GEOS-4 (bmy, 6/18/03)
!  (7 ) Now modified for GEOS-5 and GCAP met fields (swu, bmy, 5/25/05)
!  (8 ) Now return TRUE only for land boxes (w/ no ice) (bmy, 8/10/05)
!  (9 ) Now use NINT to round LWI for GEOS-4/GEOS-5 (ltm, bmy, 5/9/06)
!  (10) Remove support for GEOS-1 and GEOS-STRAT met fields (bmy, 8/4/06)
!  16 Aug 2010 - R. Yantosca - Added ProTeX headers
!  25 Aug 2010 - R. Yantosca - Treat MERRA in the same way as GEOS-5
!  06 Feb 2012 - R. Yantosca - Treat GEOS-5.7.x in the same way as MERRA/GEOS-5
!  28 Feb 2012 - R. Yantosca - Removed support for GEOS-3
!  09 Nov 2012 - M. Payer    - Replaced all met field arrays with State_Met
!                              derived type object
!EOP
!------------------------------------------------------------------------------
!BOC

      ! LWI=1 and ALBEDO less than 69.5% is a LAND box
      LAND = ( NINT( ExtState%WLI%Arr%Val(I,J) ) == 1   .and. &
                     ExtState%ALBD%Arr%Val(I,J)  <  0.695e+0_hp )

      END FUNCTION IS_LAND
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: is_ice
!
! !DESCRIPTION: Function IS\_ICE returns TRUE if surface grid box (I,J)
!  contains either land-ice or sea-ice.
!\\
!\\
! !INTERFACE:
!
      FUNCTION IS_ICE( I, J, ExtState ) RESULT ( ICE )
!
! !USES:
!
!
! !INPUT PARAMETERS:
!
      INTEGER,         INTENT(IN) :: I           ! Longitude index of grid box
      INTEGER,         INTENT(IN) :: J           ! Latitude  index of grid box
      TYPE(Ext_State), POINTER    :: ExtState    ! Module options
!
! !RETURN VALUE:
!
      LOGICAL                     :: ICE         ! =T if this is an ice box
!
!
! !REVISION HISTORY:
!  09 Aug 2005 - R. Yantosca - Initial version
!  (1 ) Remove support for GEOS-1 and GEOS-STRAT met fields (bmy, 8/4/06)
!  16 Aug 2010 - R. Yantosca - Added ProTeX headers
!  25 Aug 2010 - R. Yantosca - Treat MERRA in the same way as GEOS-5
!  06 Feb 2012 - R. Yantosca - Treat GEOS-5.7.x in the same way as MERRA/GEOS-5
!  28 Feb 2012 - R. Yantosca - Removed support for GEOS-3
!EOP
!------------------------------------------------------------------------------
!BOC

      ! LWI=2 or ALBEDO > 69.5% is ice
      ICE = ( NINT( ExtState%WLI%Arr%Val(I,J) ) == 2       .or. &
                    ExtState%ALBD%Arr%Val(I,J)  >= 0.695e+0_hp )

      END FUNCTION IS_ICE
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCOX_GC_POPs_Init
!
! !DESCRIPTION: Subroutine HcoX\_GC\_POPs\_Init initializes the HEMCO
! GC\_POPs extension.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCOX_GC_POPs_Init( HcoState, ExtName, ExtState, RC )
!
! !USES:
!
    USE HCO_ExtList_Mod,   ONLY : GetExtNr
    USE HCO_STATE_MOD,     ONLY : HCO_GetExtHcoID
!
! !INPUT PARAMETERS:
!
    CHARACTER(LEN=*), INTENT(IN   )  :: ExtName     ! Extension name
    TYPE(Ext_State),  POINTER        :: ExtState    ! Module options
    TYPE(HCO_State),  POINTER        :: HcoState    ! Hemco state
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(INOUT)  :: RC

! !REVISION HISTORY:
!  19 Aug 2014 - M. Sulprizio- Initial version
!  01 May 2015 - R. Yantosca - Bug fix: need to zero arrays after allocating
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER                        :: N, nSpc, ExtNr
    CHARACTER(LEN=255)             :: MSG

    ! Arrays
    INTEGER,           ALLOCATABLE :: HcoIDs(:)
    CHARACTER(LEN=31), ALLOCATABLE :: SpcNames(:)

    ! Pointers
    TYPE(MyInst), POINTER          :: Inst

    !=======================================================================
    ! HCOX_GC_POPs_INIT begins here!
    !=======================================================================

    ! Get the extension number
    ExtNr = GetExtNr( HcoState%Config%ExtList, TRIM(ExtName) )
    IF ( ExtNr <= 0 ) RETURN

    ! Enter HEMCO
    CALL HCO_ENTER( HcoState%Config%Err, 'HcoX_GC_POPs_Init (hcox_gc_POPs_mod.F90)', RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Create Instance
    Inst => NULL()
    CALL InstCreate ( ExtNr, ExtState%GC_POPs, Inst, RC )
    IF ( RC /= HCO_SUCCESS ) THEN
       CALL HCO_ERROR ( HcoState%Config%Err, 'Cannot create GC_POPs instance', RC )
       RETURN
    ENDIF
    ! Also fill ExtNrSS - this is the same as the parent ExtNr
    Inst%ExtNr = ExtNr

    ! Set species IDs
    CALL HCO_GetExtHcoID( HcoState, Inst%ExtNr, HcoIDs, SpcNames, nSpc, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Verbose mode
    IF ( HcoState%amIRoot ) THEN
       MSG = 'Use GC_POPs emissions module (extension module)'
       CALL HCO_MSG(HcoState%Config%Err,MSG )

       MSG = 'Use the following species (Name: HcoID):'
       CALL HCO_MSG(HcoState%Config%Err,MSG)
       DO N = 1, nSpc
          WRITE(MSG,*) TRIM(SpcNames(N)), ':', HcoIDs(N)
          CALL HCO_MSG(HcoState%Config%Err,MSG)
       ENDDO
    ENDIF

    ! Set up tracer indices
    DO N = 1, nSpc
       SELECT CASE( TRIM( SpcNames(N) ) )
          CASE( 'POPG' )
             Inst%IDTPOPG     = HcoIDs(N)
          CASE( 'POPPOCPO' )
             Inst%IDTPOPPOCPO = HcoIDs(N)
          CASE( 'POPPBCPO' )
             Inst%IDTPOPPBCPO = HcoIDs(N)
          CASE DEFAULT
             ! Do nothing
       END SELECT
    ENDDO

    ! ERROR: POPG tracer is not found!
    IF ( Inst%IDTPOPG <= 0 ) THEN
       RC = HCO_FAIL
       CALL HCO_ERROR( HcoState%Config%Err, 'Cannot find POPG tracer in list of species!', RC )
       RETURN
    ENDIF

    ! ERROR! POPPOCPO tracer is not found
    IF ( Inst%IDTPOPPOCPO <= 0 ) THEN
       RC = HCO_FAIL
       CALL HCO_ERROR( HcoState%Config%Err, 'Cannot find POPPOCPO tracer in list of species!', RC )
       RETURN
    ENDIF

    ! ERROR! POPPBCPO tracer is not found
    IF ( Inst%IDTPOPPBCPO <= 0 ) THEN
       RC = HCO_FAIL
       CALL HCO_ERROR( HcoState%Config%Err, 'Cannot find POPPBCPO tracer in list of species!', RC )
       RETURN
    ENDIF

    !=======================================================================
    ! Activate this module and the fields of ExtState that it uses
    !=======================================================================

    ! Activate met fields required by this extension
    ExtState%POPG%DoUse        = .TRUE.
    ExtState%ALBD%DoUse        = .TRUE.
    ExtState%AIRVOL%DoUse      = .TRUE.
    ExtState%AIRDEN%DoUse      = .TRUE.
    ExtState%FRAC_OF_PBL%DoUse = .TRUE.
    ExtState%FRLAKE%DoUse      = .TRUE.
    ExtState%LAI%DoUse         = .TRUE.
    ExtState%PSC2_WET%DoUse    = .TRUE.
    ExtState%SNOWHGT%DoUse     = .TRUE.
    ExtState%TK%DoUse          = .TRUE.
    ExtState%TSKIN%DoUse       = .TRUE.
    ExtState%U10M%DoUse        = .TRUE.
    ExtState%V10M%DoUse        = .TRUE.
    ExtState%WLI%DoUse         = .TRUE.

    !=======================================================================
    ! Initialize data arrays
    !=======================================================================

    ALLOCATE( Inst%EPOP_G ( HcoState%NX, HcoState%NY, HcoState%NZ ), STAT=RC )
    IF ( RC /= 0 ) THEN
       CALL HCO_ERROR ( HcoState%Config%Err, 'Cannot allocate EPOP_G', RC )
       RETURN
    ENDIF
    Inst%EPOP_G = 0.0e0_hp

    ALLOCATE( Inst%EPOP_OC( HcoState%NX, HcoState%NY, HcoState%NZ ), STAT=RC )
    IF ( RC /= 0 ) THEN
       CALL HCO_ERROR ( HcoState%Config%Err, 'Cannot allocate EPOP_OC', RC )
       RETURN
    ENDIF
    Inst%EPOP_OC = 0.0e0_hp

    ALLOCATE( Inst%EPOP_BC( HcoState%NX, HcoState%NY, HcoState%NZ ), STAT=RC )
    IF ( RC /= 0 ) THEN
       CALL HCO_ERROR ( HcoState%Config%Err, 'Cannot allocate EPOP_BC', RC )
       RETURN
    ENDIF
    Inst%EPOP_BC = 0.0e0_hp

    ALLOCATE( Inst%EPOP_VEG( HcoState%NX, HcoState%NY ), STAT=RC )
    IF ( RC /= 0 ) THEN
       CALL HCO_ERROR ( HcoState%Config%Err, 'Cannot allocate EPOP_VEG', RC )
       RETURN
    ENDIF
    Inst%EPOP_VEG = 0.0e0_hp

    ALLOCATE( Inst%EPOP_LAKE( HcoState%NX, HcoState%NY ), STAT=RC )
    IF ( RC /= 0 ) THEN
       CALL HCO_ERROR ( HcoState%Config%Err, 'Cannot allocate EPOP_LAKE', RC )
       RETURN
    ENDIF
    Inst%EPOP_LAKE = 0.0e0_hp

    ALLOCATE( Inst%EPOP_SOIL( HcoState%NX, HcoState%NY ), STAT=RC )
    IF ( RC /= 0 ) THEN
       CALL HCO_ERROR ( HcoState%Config%Err, 'Cannot allocate EPOP_SOIL', RC )
       RETURN
    ENDIF
    Inst%EPOP_SOIL = 0.0e0_hp

    ALLOCATE( Inst%EPOP_OCEAN( HcoState%NX, HcoState%NY ), STAT=RC )
    IF ( RC /= 0 ) THEN
       CALL HCO_ERROR ( HcoState%Config%Err, 'Cannot allocate EPOP_OCEAN', RC )
       RETURN
    ENDIF
    Inst%EPOP_OCEAN = 0.0e0_hp

    ALLOCATE( Inst%EPOP_SNOW( HcoState%NX, HcoState%NY ), STAT=RC )
    IF ( RC /= 0 ) THEN
       CALL HCO_ERROR ( HcoState%Config%Err, 'Cannot allocate EPOP_SNOW', RC )
       RETURN
    ENDIF
    Inst%EPOP_SNOW = 0.0e0_hp

    ALLOCATE( Inst%SUM_OC_EM( HcoState%NX, HcoState%NY ), STAT=RC )
    IF ( RC /= 0 ) THEN
       CALL HCO_ERROR ( HcoState%Config%Err, 'Cannot allocate SUM_OC_EM', RC )
       RETURN
    ENDIF
    Inst%SUM_OC_EM = 0.0e0_hp

    ALLOCATE( Inst%SUM_BC_EM( HcoState%NX, HcoState%NY ), STAT=RC )
    IF ( RC /= 0 ) THEN
       CALL HCO_ERROR ( HcoState%Config%Err, 'Cannot allocate SUM_BC_EM', RC )
       RETURN
    ENDIF
    Inst%SUM_BC_EM = 0.0e0_hp

    ALLOCATE( Inst%SUM_G_EM( HcoState%NX, HcoState%NY ), STAT=RC )
    IF ( RC /= 0 ) THEN
       CALL HCO_ERROR ( HcoState%Config%Err, 'Cannot allocate SUM_G_EM', RC )
       RETURN
    ENDIF
    Inst%SUM_G_EM = 0.0e0_hp

    ALLOCATE( Inst%SUM_OF_ALL( HcoState%NX, HcoState%NY ), STAT=RC )
    IF ( RC /= 0 ) THEN
       CALL HCO_ERROR ( HcoState%Config%Err, 'Cannot allocate SUM_OF_ALL', RC )
       RETURN
    ENDIF
    Inst%SUM_OF_ALL = 0.0e0_hp

    ALLOCATE( Inst%EMIS_SOIL( HcoState%NX, HcoState%NY ), STAT=RC )
    IF ( RC /= 0 ) THEN
       CALL HCO_ERROR ( HcoState%Config%Err, 'Cannot allocate EMIS_SOIL', RC )
       RETURN
    ENDIF
    Inst%EMIS_SOIL = 0.0e0_hp

    ALLOCATE( Inst%FLUX_SOIL2AIR( HcoState%NX, HcoState%NY ), STAT=RC )
    IF ( RC /= 0 ) THEN
       CALL HCO_ERROR ( HcoState%Config%Err, 'Cannot allocate FLUX_SOIL2AIR', RC )
       RETURN
    ENDIF
    Inst%FLUX_SOIL2AIR = 0.0e0_hp

    ALLOCATE( Inst%FLUX_AIR2SOIL( HcoState%NX, HcoState%NY ), STAT=RC )
    IF ( RC /= 0 ) THEN
       CALL HCO_ERROR ( HcoState%Config%Err, 'Cannot allocate FLUX_AIR2SOIL', RC )
       RETURN
    ENDIF
    Inst%FLUX_AIR2SOIL = 0.0e0_hp

    ALLOCATE( Inst%FUG_SOILAIR( HcoState%NX, HcoState%NY ), STAT=RC )
    IF ( RC /= 0 ) THEN
       CALL HCO_ERROR ( HcoState%Config%Err, 'Cannot allocate FUG_SOILAIR', RC )
       RETURN
    ENDIF
    Inst%FUG_SOILAIR = 0.0e0_hp

    ALLOCATE( Inst%EMIS_LAKE( HcoState%NX, HcoState%NY ), STAT=RC )
    IF ( RC /= 0 ) THEN
       CALL HCO_ERROR ( HcoState%Config%Err, 'Cannot allocate EMIS_LAKE', RC )
       RETURN
    ENDIF
    Inst%EMIS_LAKE = 0.0e0_hp

    ALLOCATE( Inst%FLUX_LAKE2AIR( HcoState%NX, HcoState%NY ), STAT=RC )
    IF ( RC /= 0 ) THEN
       CALL HCO_ERROR ( HcoState%Config%Err, 'Cannot allocate FLUX_LAKE2AIR', RC )
       RETURN
    ENDIF
    Inst%FLUX_LAKE2AIR = 0.0e0_hp

    ALLOCATE( Inst%FLUX_AIR2LAKE( HcoState%NX, HcoState%NY ), STAT=RC )
    IF ( RC /= 0 ) THEN
       CALL HCO_ERROR ( HcoState%Config%Err, 'Cannot allocate FLUX_AIR2LAKE', RC )
       RETURN
    ENDIF
    Inst%FLUX_AIR2LAKE = 0.0e0_hp

    ALLOCATE( Inst%FUG_LAKEAIR( HcoState%NX, HcoState%NY ), STAT=RC )
    IF ( RC /= 0 ) THEN
       CALL HCO_ERROR ( HcoState%Config%Err, 'Cannot allocate FUG_LAKEAIR', RC )
       RETURN
    ENDIF
    Inst%FUG_LAKEAIR = 0.0e0_hp

    ALLOCATE( Inst%EMIS_LEAF( HcoState%NX, HcoState%NY ), STAT=RC )
    IF ( RC /= 0 ) THEN
       CALL HCO_ERROR ( HcoState%Config%Err, 'Cannot allocate EMIS_LEAF', RC )
       RETURN
    ENDIF
    Inst%EMIS_LEAF = 0.0e0_hp

    ALLOCATE( Inst%FLUX_LEAF2AIR( HcoState%NX, HcoState%NY ), STAT=RC )
    IF ( RC /= 0 ) THEN
       CALL HCO_ERROR ( HcoState%Config%Err, 'Cannot allocate FLUX_LEAF2AIR', RC )
       RETURN
    ENDIF
    Inst%FLUX_LEAF2AIR = 0.0e0_hp

    ALLOCATE( Inst%FLUX_AIR2LEAF( HcoState%NX, HcoState%NY ), STAT=RC )
    IF ( RC /= 0 ) THEN
       CALL HCO_ERROR ( HcoState%Config%Err, 'Cannot allocate FLUX_AIR2LEAF', RC )
       RETURN
    ENDIF
    Inst%FLUX_AIR2LEAF = 0.0e0_hp

    ALLOCATE( Inst%FUG_LEAFAIR  ( HcoState%NX, HcoState%NY ), STAT=RC )
    IF ( RC /= 0 ) THEN
       CALL HCO_ERROR ( HcoState%Config%Err, 'Cannot allocate FUG_LEAFAIR', RC )
       RETURN
    ENDIF
    Inst%FUG_LEAFAIR = 0.0e0_hp

    !=======================================================================
    ! Leave w/ success
    !=======================================================================
    IF ( ALLOCATED( HcoIDs   ) ) DEALLOCATE( HcoIDs   )
    IF ( ALLOCATED( SpcNames ) ) DEALLOCATE( SpcNames )

    ! Nullify pointers
    Inst           => NULL()

    CALL HCO_LEAVE( HcoState%Config%Err,RC )

  END SUBROUTINE HCOX_GC_POPs_Init
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCOX_GC_POPs_Final
!
! !DESCRIPTION: Subroutine HcoX\_GC\_POPs\_Final finalizes the HEMCO
!  extension for the GEOS-Chem POPs specialty simulation.  All module
!  arrays will be deallocated.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCOX_GC_POPs_Final( ExtState )
!
! !INPUT PARAMETERS:
!
    TYPE(Ext_State),  POINTER       :: ExtState   ! Module options
!
! !REVISION HISTORY:
!  19 Aug 2014 - M. Sulprizio- Initial version
!EOP
!------------------------------------------------------------------------------
!BOC

    !=======================================================================
    ! HCOX_GC_POPs_FINAL begins here!
    !=======================================================================

    CALL InstRemove( ExtState%GC_POPs )

  END SUBROUTINE HCOX_GC_POPs_Final
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: InstGet
!
! !DESCRIPTION: Subroutine InstGet returns a poiner to the desired instance.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE InstGet ( Instance, Inst, RC, PrevInst )
!
! !INPUT PARAMETERS:
!
    INTEGER                             :: Instance
    TYPE(MyInst),     POINTER           :: Inst
    INTEGER                             :: RC
    TYPE(MyInst),     POINTER, OPTIONAL :: PrevInst
!
! !REVISION HISTORY:
!  18 Feb 2016 - C. Keller   - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
    TYPE(MyInst),     POINTER    :: PrvInst

    !=================================================================
    ! InstGet begins here!
    !=================================================================

    ! Get instance. Also archive previous instance.
    PrvInst => NULL()
    Inst    => AllInst
    DO WHILE ( ASSOCIATED(Inst) )
       IF ( Inst%Instance == Instance ) EXIT
       PrvInst => Inst
       Inst    => Inst%NextInst
    END DO
    IF ( .NOT. ASSOCIATED( Inst ) ) THEN
       RC = HCO_FAIL
       RETURN
    ENDIF

    ! Pass output arguments
    IF ( PRESENT(PrevInst) ) PrevInst => PrvInst

    ! Cleanup & Return
    PrvInst => NULL()
    RC = HCO_SUCCESS

  END SUBROUTINE InstGet
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: InstCreate
!
! !DESCRIPTION: Subroutine InstCreate creates a new instance.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE InstCreate ( ExtNr, Instance, Inst, RC )
!
! !INPUT PARAMETERS:
!
    INTEGER,       INTENT(IN)       :: ExtNr
!
! !OUTPUT PARAMETERS:
!
    INTEGER,       INTENT(  OUT)    :: Instance
    TYPE(MyInst),  POINTER          :: Inst
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,       INTENT(INOUT)    :: RC
!
! !REVISION HISTORY:
!  18 Feb 2016 - C. Keller   - Initial version
!  26 Oct 2016 - R. Yantosca - Don't nullify local ptrs in declaration stmts
!EOP
!------------------------------------------------------------------------------
!BOC
    TYPE(MyInst), POINTER          :: TmpInst
    INTEGER                        :: nnInst

    !=================================================================
    ! InstCreate begins here!
    !=================================================================

    ! ----------------------------------------------------------------
    ! Generic instance initialization
    ! ----------------------------------------------------------------

    ! Initialize
    Inst => NULL()

    ! Get number of already existing instances
    TmpInst => AllInst
    nnInst = 0
    DO WHILE ( ASSOCIATED(TmpInst) )
       nnInst  =  nnInst + 1
       TmpInst => TmpInst%NextInst
    END DO

    ! Create new instance
    ALLOCATE(Inst)
    Inst%Instance = nnInst + 1
    Inst%ExtNr    = ExtNr

    ! Attach to instance list
    Inst%NextInst => AllInst
    AllInst       => Inst

    ! Update output instance
    Instance = Inst%Instance

    ! ----------------------------------------------------------------
    ! Type specific initialization statements follow below
    ! ----------------------------------------------------------------

    ! Return w/ success
    RC = HCO_SUCCESS

  END SUBROUTINE InstCreate
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!BOP
!
! !IROUTINE: InstRemove
!
! !DESCRIPTION: Subroutine InstRemove creates a new instance.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE InstRemove ( Instance )
!
! !INPUT PARAMETERS:
!
    INTEGER                         :: Instance
!
! !REVISION HISTORY:
!  18 Feb 2016 - C. Keller   - Initial version
!  26 Oct 2016 - R. Yantosca - Don't nullify local ptrs in declaration stmts
!EOP
!------------------------------------------------------------------------------
!BOC
    INTEGER                     :: RC
    TYPE(MyInst), POINTER       :: PrevInst
    TYPE(MyInst), POINTER       :: Inst

    !=================================================================
    ! InstRemove begins here!
    !=================================================================

    ! Init
    PrevInst => NULL()
    Inst     => NULL()

    ! Get instance. Also archive previous instance.
    CALL InstGet ( Instance, Inst, RC, PrevInst=PrevInst )

    ! Instance-specific deallocation
    IF ( ASSOCIATED(Inst) ) THEN

       ! Pop off instance from list
       IF ( ASSOCIATED(PrevInst) ) THEN

          IF ( ASSOCIATED(Inst%EPOP_G       ) ) DEALLOCATE(Inst%EPOP_G       )
          IF ( ASSOCIATED(Inst%EPOP_OC      ) ) DEALLOCATE(Inst%EPOP_OC      )
          IF ( ASSOCIATED(Inst%EPOP_BC      ) ) DEALLOCATE(Inst%EPOP_BC      )
          IF ( ASSOCIATED(Inst%SUM_OC_EM    ) ) DEALLOCATE(Inst%SUM_OC_EM    )
          IF ( ASSOCIATED(Inst%SUM_BC_EM    ) ) DEALLOCATE(Inst%SUM_BC_EM    )
          IF ( ASSOCIATED(Inst%SUM_G_EM     ) ) DEALLOCATE(Inst%SUM_G_EM     )
          IF ( ASSOCIATED(Inst%SUM_OF_ALL   ) ) DEALLOCATE(Inst%SUM_OF_ALL   )
          IF ( ASSOCIATED(Inst%EMIS_SOIL    ) ) DEALLOCATE(Inst%EMIS_SOIL    )
          IF ( ASSOCIATED(Inst%FLUX_SOIL2AIR) ) DEALLOCATE(Inst%FLUX_SOIL2AIR)
          IF ( ASSOCIATED(Inst%FLUX_AIR2SOIL) ) DEALLOCATE(Inst%FLUX_AIR2SOIL)
          IF ( ASSOCIATED(Inst%FUG_SOILAIR  ) ) DEALLOCATE(Inst%FUG_SOILAIR  )
          IF ( ASSOCIATED(Inst%EMIS_LAKE    ) ) DEALLOCATE(Inst%EMIS_LAKE    )
          IF ( ASSOCIATED(Inst%FLUX_LAKE2AIR) ) DEALLOCATE(Inst%FLUX_LAKE2AIR)
          IF ( ASSOCIATED(Inst%FLUX_AIR2LAKE) ) DEALLOCATE(Inst%FLUX_AIR2LAKE)
          IF ( ASSOCIATED(Inst%FUG_LAKEAIR  ) ) DEALLOCATE(Inst%FUG_LAKEAIR  )
          IF ( ASSOCIATED(Inst%EMIS_LEAF    ) ) DEALLOCATE(Inst%EMIS_LEAF    )
          IF ( ASSOCIATED(Inst%FLUX_LEAF2AIR) ) DEALLOCATE(Inst%FLUX_LEAF2AIR)
          IF ( ASSOCIATED(Inst%FLUX_AIR2LEAF) ) DEALLOCATE(Inst%FLUX_AIR2LEAF)
          IF ( ASSOCIATED(Inst%FUG_LEAFAIR  ) ) DEALLOCATE(Inst%FUG_LEAFAIR  )

          PrevInst%NextInst => Inst%NextInst
       ELSE
          AllInst => Inst%NextInst
       ENDIF
       DEALLOCATE(Inst)
       Inst => NULL()
    ENDIF

   END SUBROUTINE InstRemove
!EOC
END MODULE HCOX_GC_POPs_Mod
