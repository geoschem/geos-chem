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
! !REMARKS:
!  POPs Tracers
!  ============================================================================
!  (1 ) POPG   : Gaseous POP    - total tracer  
!  (2 ) POPPOC : OC-sorbed POP  - total tracer
!  (3 ) POPPBC : BC-sorbed POP  - total tracer
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
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !MODULE VARIABLES:
!
  ! Fields required by module
  INTEGER                       :: ExtNr     = -1   ! HEMCO Extension number
  INTEGER                       :: IDTPOPPOC = -1   ! Index # for POPPOC tracer
  INTEGER                       :: IDTPOPPBC = -1   ! Index # for POPPBC tracer
  INTEGER                       :: IDTPOPG   = -1   ! Index # for POPG   tracer

  REAL*8,  PARAMETER            :: SMALLNUM = 1D-20

  ! Pointers to emission arrays read from disk
  REAL(sp), POINTER             :: POP_TOT_EM(:,:) => NULL()
  REAL(sp), POINTER             :: C_OC(:,:,:)     => NULL()
  REAL(sp), POINTER             :: C_BC(:,:,:)     => NULL()

  ! Calculated emissions of OC-phase, BC-phase, and gas-phase POPs
  REAL(hp), ALLOCATABLE, TARGET :: EPOP_OC(:,:,:)
  REAL(hp), ALLOCATABLE, TARGET :: EPOP_BC(:,:,:)
  REAL(hp), ALLOCATABLE, TARGET :: EPOP_G (:,:,:)

  ! For diagnostics
  REAL(hp), ALLOCATABLE, TARGET :: SUM_OC_EM (:,:)
  REAL(hp), ALLOCATABLE, TARGET :: SUM_BC_EM (:,:)
  REAL(hp), ALLOCATABLE, TARGET :: SUM_G_EM  (:,:)
  REAL(hp), ALLOCATABLE, TARGET :: SUM_OF_ALL(:,:)

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
  SUBROUTINE HCOX_GC_POPs_Run( am_I_Root, ExtState, HcoState, RC )
!
! !USES:
!
    ! HEMCO modules
    USE HCO_EmisList_Mod,  ONLY : HCO_GetPtr
    USE HCO_FluxArr_Mod,   ONLY : HCO_EmisAdd
!
! !INPUT PARAMETERS:
!
    LOGICAL,          INTENT(IN   ) :: am_I_Root   ! Are we on the root CPU?
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
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
    ! Universal gas constant for adjusting KOA for temp: 8.3145 [J/mol/K]
    REAL*8, PARAMETER :: R          = 8.31d0  

    ! Density of octanol, needed for partitioning into OC: 820 [kg/m^3]
    REAL*8, PARAMETER :: DENS_OCT   = 82d1

    ! Density of BC, needed for partitioning onto BC: 1 [kg/L] or 1000 [kg/m^3]
    ! From Lohmann and Lammel, Environ. Sci. Technol., 2004, 38:3793-3803.
    REAL*8, PARAMETER :: DENS_BC    = 1d3

!
! !LOCAL VARIABLES:
!
    INTEGER           :: I, J, L
    INTEGER           :: PBL_MAX
    INTEGER           :: MONTH,            YEAR
    REAL*8            :: F_OF_PBL,         TK
    REAL*8            :: T_POP
    REAL*8            :: C_OC1,            C_BC1
    REAL*8            :: C_OC2,            C_BC2
    REAL*8            :: F_POP_OC,         F_POP_BC 
    REAL*8            :: F_POP_G,          AIR_VOL
    REAL*8            :: KOA_T,            KBC_T
    REAL*8            :: KOC_BC_T,         KBC_OC_T
    REAL*8            :: VR_OC_AIR,        VR_OC_BC
    REAL*8            :: VR_BC_AIR,        VR_BC_OC
    REAL*8            :: SUM_F
    REAL*8            :: OC_AIR_RATIO,     OC_BC_RATIO
    REAL*8            :: BC_AIR_RATIO,     BC_OC_RATIO
    LOGICAL, SAVE     :: FIRST = .TRUE.
    LOGICAL           :: aIR

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
    REAL*8            :: DEL_H

    ! KOA_298 for partitioning of gas phase POP to atmospheric OC
    ! KOA_298 = Cpop in octanol/Cpop in atmosphere at 298 K 
    ! For PHENANTHRENE:
    ! log KOA_298 = 7.64, or 4.37*10^7 [unitless]
    ! For PYRENE:
    ! log KOA_298 = 8.86, or 7.24*10^8 [unitless]
    ! For BENZO[a]PYRENE:
    ! log KOA_298 = 11.48, or 3.02*10^11 [unitless]
    ! (Ma et al., J. Chem. Eng. Data, 2010, 55:819-825).
    REAL*8            :: KOA_298

    ! KBC_298 for partitioning of gas phase POP to atmospheric BC
    ! KBC_298 = Cpop in black carbon/Cpop in atmosphere at 298 K
    ! For PHENANTHRENE:
    ! log KBC_298 = 10.0, or 1.0*10^10 [unitless]
    ! For PYRENE:
    ! log KBC_298 = 11.0, or 1.0*10^11 [unitless]
    ! For BENZO[a]PYRENE:
    ! log KBC_298 = 13.9, or 7.94*10^13 [unitless]
    ! (Lohmann and Lammel, EST, 2004, 38:3793-3802)
    REAL*8            :: KBC_298

    ! Pointers for diagnostics
    REAL(hp), POINTER :: Arr3D(:,:,:) => NULL()

    !=======================================================================
    ! HCOX_GC_POPs_RUN begins here!
    !=======================================================================

    ! Return if extension not turned on
    IF ( .NOT. ExtState%GC_POPs ) RETURN

    ! Enter
    CALL HCO_ENTER( 'HCOX_GC_POPs_Run (hcox_gc_POPs_mod.F90)', RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! am I root? 
    aIR = am_I_Root

    DEL_H   = ExtState%POP_DEL_H
    KOA_298 = ExtState%POP_KOA
    KBC_298 = ExtState%POP_KBC

    !=======================================================================
    ! Get pointers to gridded data imported through config. file
    !=======================================================================
    IF ( FIRST ) THEN

       CALL HCO_GetPtr( aIR, 'TOT_POP', POP_TOT_EM, RC )
       IF ( RC /= HCO_SUCCESS ) RETURN

       CALL HCO_GetPtr( aIR, 'GLOBAL_OC', C_OC, RC )
       IF ( RC /= HCO_SUCCESS ) RETURN

       CALL HCO_GetPtr( aIR, 'GLOBAL_BC', C_BC, RC )
       IF ( RC /= HCO_SUCCESS ) RETURN

       FIRST = .FALSE.

    ENDIF

    ! Maximum extent of the PBL [model level]
    IF ( .NOT. ASSOCIATED(ExtState%PBL_MAX) ) THEN
       CALL HCO_ERROR ( 'PBL_MAX not defined in ExtState!', RC )
       RETURN
    ELSE
       PBL_MAX = DBLE( ExtState%PBL_MAX )
    ENDIF

    ! Loop over grid boxes
    DO J = 1, HcoState%Ny
    DO I = 1, HcoState%Nx

       F_OF_PBL = 0d0 
       T_POP = 0d0       

       ! Here, save the total from the emissions array
       ! into the T_POP variable [kg/m2/s]
       T_POP = POP_TOT_EM(I,J)

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
          KOA_T = KOA_298 * EXP((-DEL_H/R) * ((1d0/TK) - (1d0/298d0)))

          ! Define KBC_T, the BC-air partition coeff at temp T [unitless]
          ! TURN OFF TEMPERATURE DEPENDENCY FOR SENSITIVITY ANALYSIS
          KBC_T = KBC_298 * EXP((-DEL_H/R) * ((1d0/TK) - (1d0/298d0)))

          ! Define KOC_BC_T, theoretical OC-BC part coeff at temp T [unitless]
          KOC_BC_T = KOA_T / KBC_T

          ! Define KBC_OC_T, theoretical BC_OC part coeff at temp T [unitless]
          KBC_OC_T = 1d0 / KOC_BC_T

          ! Get monthly mean OC and BC concentrations [kg/box]
          C_OC1    = C_OC(I,J,L)
          C_BC1    = C_BC(I,J,L)
           
          ! Make sure OC is not negative
          C_OC1 = MAX( C_OC1, 0d0 )

          ! Convert C_OC and C_BC units to volume per box 
          ! [m^3 OC or BC/box]
          !C_OC(I,J,L)        = GET_OC(I,J,L) / DENS_OCT
          !C_BC(I,J,L)        = GET_BC(I,J,L) / DENS_BC
          C_OC2    = C_OC1 / DENS_OCT
          C_BC2    = C_BC1 / DENS_BC

          ! Get air volume (m^3)
          AIR_VOL   = ExtState%AIRVOL%Arr%Val(I,J,L) 

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
          OC_AIR_RATIO = 1d0 / (KOA_T    * VR_OC_AIR) 
          OC_BC_RATIO  = 1d0 / (KOC_BC_T * VR_OC_BC) 
  
          BC_AIR_RATIO = 1d0 / (KBC_T * VR_BC_AIR) 
          BC_OC_RATIO  = 1d0 / (KBC_OC_T * VR_BC_OC)

          ! If there are zeros in OC or BC concentrations, make sure they
          ! don't cause problems with phase fractions
          IF ( C_OC1 > SMALLNUM .and. C_BC1 > SMALLNUM ) THEN
             F_POP_OC  = 1d0 / (1d0 + OC_AIR_RATIO + OC_BC_RATIO) 
             F_POP_BC  = 1d0 / (1d0 + BC_AIR_RATIO + BC_OC_RATIO)
         
          ELSE IF ( C_OC1 > SMALLNUM .and. C_BC1 .le. SMALLNUM ) THEN
             F_POP_OC  = 1d0 / (1d0 + OC_AIR_RATIO)
             F_POP_BC  = SMALLNUM           

          ELSE IF ( C_OC1 .le. SMALLNUM .and. C_BC1 > SMALLNUM ) THEN
             F_POP_OC  = SMALLNUM
             F_POP_BC  = 1d0 / (1d0 + BC_AIR_RATIO)

          ELSE IF ( C_OC1 .le. SMALLNUM .and. C_BC1 .le. SMALLNUM) THEN
             F_POP_OC = SMALLNUM
             F_POP_BC = SMALLNUM
          ENDIF

          ! Gas-phase:
          F_POP_G   = 1d0 - F_POP_OC - F_POP_BC

          ! Check that sum of fractions equals 1
          SUM_F = F_POP_OC + F_POP_BC + F_POP_G                
            
          ! Fraction of PBL that box (I,J,L) makes up [unitless]
          F_OF_PBL = ExtState%FRAC_OF_PBL%Arr%Val(I,J,L) 

          ! Calculate rates of POP emissions in each phase [kg/s]
          ! OC-phase:
          EPOP_OC(I,J,L) = F_POP_OC * F_OF_PBL * T_POP 

          ! BC-phase
          EPOP_BC(I,J,L) = F_POP_BC * F_OF_PBL * T_POP

          ! Gas-phase
          EPOP_G(I,J,L)  = F_POP_G  * F_OF_PBL * T_POP

       ENDDO

       !==================================================================
       ! Sum different POPs emissions phases (OC, BC, and gas phase)
       ! through bottom layer to top of PBL for storage in ND53 diagnostic
       !==================================================================

       SUM_OC_EM(I,J) =  SUM(EPOP_OC(I,J,1:PBL_MAX))  
       SUM_BC_EM(I,J) =  SUM(EPOP_BC(I,J,1:PBL_MAX))
       SUM_G_EM(I,J)  =  SUM(EPOP_G(I,J,1:PBL_MAX))           
       
       SUM_OF_ALL(I,J) = SUM_OC_EM(I,J) + SUM_BC_EM(I,J) + &
                         SUM_G_EM(I,J)

       ! Check that sum thru PBL is equal to original emissions array
       ! NOTE: Prevent div-by-zero floating point error (bmy, 4/14/14)
       IF ( SUM_OF_ALL(I,J) > 0d0 ) THEN 
          SUM_OF_ALL(I,J) = POP_TOT_EM(I,J) / SUM_OF_ALL(I,J)
       ENDIF

    ENDDO
    ENDDO

    !=======================================================================
    ! Add POPs emissions to HEMCO data structure & diagnostics
    !=======================================================================

    !----------------------
    ! OC-PHASE EMISSIONS
    !----------------------
    IF ( IDTPOPPOC > 0 ) THEN

       ! Add flux to emissions array
       Arr3D => EPOP_OC(:,:,:)
       CALL HCO_EmisAdd( am_I_Root, HcoState, Arr3D, IDTPOPPOC, RC, ExtNr=ExtNr )
       Arr3D => NULL()
       IF ( RC /= HCO_SUCCESS ) THEN
          CALL HCO_ERROR( 'HCO_EmisAdd error: EPOP_OC', RC )
          RETURN 
       ENDIF
    ENDIF

    !----------------------
    ! BC-PHASE EMISSIONS
    !----------------------
    IF ( IDTPOPPBC > 0 ) THEN

       ! Add flux to emissions array
       Arr3D => EPOP_BC(:,:,:)
       CALL HCO_EmisAdd( am_I_Root, HcoState, Arr3D, IDTPOPPBC, RC, ExtNr=ExtNr )
       Arr3D => NULL()
       IF ( RC /= HCO_SUCCESS ) THEN
          CALL HCO_ERROR( 'HCO_EmisAdd error: EPOP_BC', RC )
          RETURN 
       ENDIF
    ENDIF

    !----------------------
    ! GASEOUS EMISSIONS
    !----------------------
    IF ( IDTPOPG > 0 ) THEN

       ! Add flux to emissions array
       Arr3D => EPOP_G(:,:,:)
       CALL HCO_EmisAdd( am_I_Root, HcoState, Arr3D, IDTPOPG, RC, ExtNr=ExtNr )
       Arr3D => NULL()
       IF ( RC /= HCO_SUCCESS ) THEN
          CALL HCO_ERROR( 'HCO_EmisAdd error: EPOP_G', RC )
          RETURN 
       ENDIF

    ENDIF

    !=======================================================================
    ! Cleanup & quit
    !=======================================================================

    ! Return w/ success
    CALL HCO_LEAVE ( RC )

  END SUBROUTINE HCOX_GC_POPs_Run
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
  SUBROUTINE HCOX_GC_POPs_Init( am_I_Root, HcoState, ExtName, ExtState, RC )
!
! !USES:
!
    USE HCO_ExtList_Mod,   ONLY : GetExtNr
    USE HCO_STATE_MOD,     ONLY : HCO_GetExtHcoID
!
! !INPUT PARAMETERS:
!
    LOGICAL,          INTENT(IN   )  :: am_I_Root
    CHARACTER(LEN=*), INTENT(IN   )  :: ExtName     ! Extension name
    TYPE(Ext_State),  POINTER        :: ExtState    ! Module options      
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(HCO_State),  POINTER        :: HcoState    ! Hemco state 
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
    INTEGER                        :: N, nSpc
    CHARACTER(LEN=255)             :: MSG 

    ! Arrays
    INTEGER,           ALLOCATABLE :: HcoIDs(:)
    CHARACTER(LEN=31), ALLOCATABLE :: SpcNames(:)

    !=======================================================================
    ! HCOX_GC_POPs_INIT begins here!
    !=======================================================================

    ! Get the extension number
    ExtNr = GetExtNr( TRIM( ExtName ) )
    IF ( ExtNr <= 0 ) RETURN

    ! Enter HEMCO
    CALL HCO_ENTER( 'HcoX_GC_POPs_Init (hcox_gc_POPs_mod.F90)', RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Set species IDs      
    CALL HCO_GetExtHcoID( HcoState, ExtNr, HcoIDs, SpcNames, nSpc, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Verbose mode
    IF ( am_I_Root ) THEN
       MSG = 'Use gc_POPs emissions module (extension module)'
       CALL HCO_MSG( MSG )

       MSG = 'Use the following species (Name: HcoID):'
       CALL HCO_MSG(MSG)
       DO N = 1, nSpc
          WRITE(MSG,*) TRIM(SpcNames(N)), ':', HcoIDs(N)
          CALL HCO_MSG(MSG)
       ENDDO
    ENDIF

    ! Set up tracer indices
    DO N = 1, nSpc
       SELECT CASE( TRIM( SpcNames(N) ) )
          CASE( 'POPG' )
             IDTPOPG   = HcoIDs(N)
          CASE( 'POPPOC' )
             IDTPOPPOC = HcoIDs(N)
          CASE( 'POPPBC' )
             IDTPOPPBC = HcoIDs(N)
          CASE DEFAULT
             ! Do nothing
       END SELECT
    ENDDO

    ! ERROR: POPG tracer is not found!
    IF ( IDTPOPG <= 0 ) THEN
       RC = HCO_FAIL
       CALL HCO_ERROR( 'Cannot find POPG tracer in list of species!', RC )
       RETURN
    ENDIF
    
    ! ERROR! POPPOC tracer is not found
    IF ( IDTPOPPOC <= 0 ) THEN
       RC = HCO_FAIL
       CALL HCO_ERROR( 'Cannot find POPPOC tracer in list of species!', RC )
       RETURN
    ENDIF

    ! ERROR! POPPBC tracer is not found
    IF ( IDTPOPPBC <= 0 ) THEN
       RC = HCO_FAIL
       CALL HCO_ERROR( 'Cannot find POPPBC tracer in list of species!', RC )
       RETURN
    ENDIF

    ! Activate met fields required by this extension
    ExtState%AIRVOL%DoUse      = .TRUE. 
    ExtState%FRAC_OF_PBL%DoUse = .TRUE. 
    ExtState%TK%DoUse          = .TRUE. 

    ! Activate this extension
    ExtState%GC_POPs           = .TRUE.

    !=======================================================================
    ! Initialize data arrays
    !=======================================================================

    ALLOCATE( EPOP_G ( HcoState%NX, HcoState%NY, HcoState%NZ ), STAT=RC )
    IF ( RC /= 0 ) THEN
       CALL HCO_ERROR ( 'Cannot allocate EPOP_G', RC )
       RETURN
    ENDIF 
    EPOP_G = 0.0e0_hp

    ALLOCATE( EPOP_OC( HcoState%NX, HcoState%NY, HcoState%NZ ), STAT=RC )
    IF ( RC /= 0 ) THEN
       CALL HCO_ERROR ( 'Cannot allocate EPOP_OC', RC )
       RETURN
    ENDIF 
    EPOP_OC = 0.0e0_hp

    ALLOCATE( EPOP_BC( HcoState%NX, HcoState%NY, HcoState%NZ ), STAT=RC )
    IF ( RC /= 0 ) THEN
       CALL HCO_ERROR ( 'Cannot allocate EPOP_BC', RC )
       RETURN
    ENDIF 
    EPOP_BC = 0.0e0_hp

    ALLOCATE( SUM_OC_EM( HcoState%NX, HcoState%NY ), STAT=RC )
    IF ( RC /= 0 ) THEN
       CALL HCO_ERROR ( 'Cannot allocate SUM_OC_EM', RC )
       RETURN
    ENDIF 
    SUM_OC_EM = 0.0e0_hp

    ALLOCATE( SUM_BC_EM( HcoState%NX, HcoState%NY ), STAT=RC )
    IF ( RC /= 0 ) THEN
       CALL HCO_ERROR ( 'Cannot allocate SUM_BC_EM', RC )
       RETURN
    ENDIF 
    SUM_BC_EM = 0.0e0_hp

    ALLOCATE( SUM_G_EM( HcoState%NX, HcoState%NY ), STAT=RC )
    IF ( RC /= 0 ) THEN
       CALL HCO_ERROR ( 'Cannot allocate SUM_G_EM', RC )
       RETURN
    ENDIF     
    SUM_G_EM = 0.0e0_hp

    ALLOCATE( SUM_OF_ALL( HcoState%NX, HcoState%NY ), STAT=RC )
    IF ( RC /= 0 ) THEN
       CALL HCO_ERROR ( 'Cannot allocate SUM_OF_ALL', RC )
       RETURN
    ENDIF 
    SUM_OF_ALL = 0.0e0_hp

    !=======================================================================
    ! Leave w/ success
    !=======================================================================
    IF ( ALLOCATED( HcoIDs   ) ) DEALLOCATE( HcoIDs   )
    IF ( ALLOCATED( SpcNames ) ) DEALLOCATE( SpcNames )

    CALL HCO_LEAVE ( RC ) 

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
  SUBROUTINE HCOX_GC_POPs_Final()
!
! !REVISION HISTORY:
!  19 Aug 2014 - M. Sulprizio- Initial version
!EOP
!------------------------------------------------------------------------------
!BOC

    !=======================================================================
    ! HCOX_GC_POPs_FINAL begins here!
    !=======================================================================
    IF ( ALLOCATED( EPOP_G     ) ) DEALLOCATE( EPOP_G     )
    IF ( ALLOCATED( EPOP_OC    ) ) DEALLOCATE( EPOP_OC    )
    IF ( ALLOCATED( EPOP_BC    ) ) DEALLOCATE( EPOP_BC    )
    IF ( ALLOCATED( SUM_OC_EM  ) ) DEALLOCATE( SUM_OC_EM  )
    IF ( ALLOCATED( SUM_BC_EM  ) ) DEALLOCATE( SUM_BC_EM  )
    IF ( ALLOCATED( SUM_G_EM   ) ) DEALLOCATE( SUM_G_EM   )
    IF ( ALLOCATED( SUM_OF_ALL ) ) DEALLOCATE( SUM_OF_ALL )

  END SUBROUTINE HCOX_GC_POPs_Final
!EOC
END MODULE HCOX_GC_POPs_Mod
