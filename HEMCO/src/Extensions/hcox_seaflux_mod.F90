!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: hcox_seaflux_mod.F90
!
! !DESCRIPTION: Module HCOX\_SeaFlux\_Mod contains routines to calculate
! the oceanic emissions of a number of defined species.
! The oceanic flux is parameterized according to Liss and Slater, 1974:
! F = Kg * ( Cair - H Cwater )
! where F is the net flux, Kg is the exchange velocity, Cair and Cwater
! are the air and aqueous concentrations, respectively, and H is the
! dimensionless air over water Henry constant.
!\\
!\\
! This module calculates the source and sink terms separately. The source
! is given as flux, the sink as deposition rate:
! source = Kg * H * Cwater     [kg m-2 s-1]
! sink   = Kg / DEPHEIGHT      [s-1]
!
! The deposition rate is obtained by dividing the exchange velocity Kg
! by the deposition height DEPHEIGHT, e.g. the height over which
! deposition occurs. This can be either the first grid box only, or the
! entire planetary boundary layer. The HEMCO option 'PBL\_DRYDEP' determines
! which option is being used.
!\\
!\\
! Kg is calculated following Johnson, 2010, which is largely based on
! the work of Nightingale et al., 2000a/b.
! The salinity and seawater pH are currently set to constant global values
! of 35 ppt and 8.0, respectively.
! Since Kg is only little sensitive to these variables, this should not
! introduce a notable error.
!\\
!\\
! This is a HEMCO extension module that uses many of the HEMCO core
! utilities.
!\\
!\\
! Air-sea exchange is calculated for all species defined during
! extension initialization. For each species, the following parameter
! must be specified: species name, model species ID (i.e. ID of this
! species in the external model), parameterization type of Schmidt
! number in water, liquid molar volume of species, and the name of the
! field containing species sea-water concentrations. See initialization
! routine for more details.
! To add new species to this module, the abovementioned arrays have to
! be extended accordingly.
!\\
!\\
! References:
! \begin{itemize}
! \item Johnson, M.: A numerical scheme to calculate temperature and salinity
!    dependent air-water transfer velocities for any gas, Ocean Science, 6,
!    2010.
! \item Liss and Slater: Flux of gases across the air-sea interface, Nature,
!    247, 1974.
! \item Nightingale et al.: In situ evaluation of air-sea gas exchange
!    parameterizations using novel conservative and volatile tracers,
!    Global Biogeochemical Cycles, 14, 2000a.
! \item Nightingale et al.: Measurements of air-sea gas transfer during an
!    open ocean algal bloom, Geophys. Res. Lett., 27, 2000b.
! \item Saltzman et al.: Experimental determination of the diffusion
!    coefficient of dimethylsulfide in water, J. Geophys. Res., 98, 1993.
! \end{itemize}
!
! !INTERFACE:
!
MODULE HCOX_SeaFlux_Mod
!
! !USES:
!
  USE HCO_Error_MOD
  USE HCO_Diagn_MOD
  USE HCO_State_MOD,  ONLY : HCO_State
  USE HCOX_State_MOD, ONLY : Ext_State

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC  :: HCOX_SeaFlux_Init
  PUBLIC  :: HCOX_SeaFlux_Run
  PUBLIC  :: HCOX_SeaFlux_Final
!
! !PRIVATE MEMBER FUNCTIONS:
!
  PRIVATE :: Calc_SeaFlux
!
! !REVISION HISTORY:
!  16 Apr 2013 - C. Keller   - Initial version
!  01 Oct 2013 - C. Keller   - Now a HEMCO extension module
!  11 Dec 2013 - C. Keller   - Now define container name during initialization
!  01 Jul 2014 - R. Yantosca - Now use F90 free-format indentation
!  01 Jul 2014 - R. Yantosca - Cosmetic changes in ProTeX headers
!  06 Nov 2015 - C. Keller   - Now use land type definitions instead of FRCLND
!  14 Oct 2016 - C. Keller   - Now use HCO_EvalFld instead of HCO_GetPtr.
!  10 Mar 2017 - M. Sulprizio- Add fix for acetone parameterization of Schmidt
!                              number - use SCWPAR = 3 instead of 1
!  11 Sep 2018 - C. Keller   - Added instances wrapper
!  08 May 2019 - J. Fisher   - Add C1-C2 alkyl nitrates (MENO3, ETNO3)
!EOP
!------------------------------------------------------------------------------
!
! !PRIVATE TYPES:
!
  ! Ocean species
  TYPE :: OcSpec
     INTEGER            :: HcoID        ! HEMCO species ID
     CHARACTER(LEN=31)  :: OcSpcName    ! oc. species name
     CHARACTER(LEN=31)  :: OcDataName   ! seawater conc. field name
     REAL*8             :: LiqVol       ! liq. molecular volume
     INTEGER            :: SCWPAR       ! Schmidt # parameterization type
  END TYPE OcSpec

  TYPE :: MyInst
   ! Tracer IDs
   INTEGER                :: Instance
   ! Variables carrying information about ocean species
   INTEGER                :: ExtNr
   INTEGER                :: nOcSpc            ! # of ocean species
   TYPE(OcSpec), POINTER  :: OcSpecs(:)
   TYPE(MyInst), POINTER  :: NextInst => NULL()
  END TYPE MyInst

  ! Pointer to instances
  TYPE(MyInst), POINTER   :: AllInst => NULL()

CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCOX_SeaFlux_Run
!
! !DESCRIPTION: Subroutine HcoX\_SeaFlux\_Run is the run routine to
! calculate oceanic emissions for the current time step.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCOX_SeaFlux_Run( ExtState, HcoState, RC )
!
! !USES:
!
    USE HCO_FLUXARR_MOD,  ONLY : HCO_EmisAdd
    USE HCO_FLUXARR_MOD,  ONLY : HCO_DepvAdd
    USE HCO_CALC_MOD,     ONLY : HCO_EvalFld
!    USE HCO_EMISLIST_MOD, ONLY : HCO_GetPtr
!
! !INPUT PARAMETERS:
!
    TYPE(HCO_State), POINTER       :: HcoState   ! Output obj
    TYPE(Ext_State), POINTER       :: ExtState  ! Module options
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,         INTENT(INOUT) :: RC         ! Success or failure?
!
! !REVISION HISTORY:
!  16 Apr 2013 - C. Keller - Initial version
!  26 Oct 2016 - R. Yantosca - Don't nullify local ptrs in declaration stmts
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    TYPE(MyInst), POINTER :: Inst
    INTEGER               :: OcID, HcoID
    REAL(hp), TARGET      :: SOURCE(HcoState%NX,HcoState%NY)
    REAL(hp), TARGET      :: SINK  (HcoState%NX,HcoState%NY)
    REAL(hp), TARGET      :: SeaConc(HcoState%NX,HcoState%NY)
    CHARACTER(LEN=255)    :: ContName
    CHARACTER(LEN=255)    :: MSG
    LOGICAL               :: VERBOSE

    ! Pointers
    REAL(hp), POINTER     :: Arr2D(:,:)

    !=================================================================
    ! HCOX_SeaFlux_Run begins here!
    !=================================================================

    ! Enter
    CALL HCO_ENTER( HcoState%Config%Err, 'HCOX_SeaFlux_Run (hcox_seaflux_mod.F90)', RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Return if extension disabled
    IF ( ExtState%SeaFlux <= 0 ) RETURN

    ! Verbose?
    verbose = HCO_IsVerb(HcoState%Config%Err,1)

    ! Nullify
    Arr2D => NULL()

    ! Get instance
    Inst => NULL()
    CALL InstGet ( ExtState%SeaFlux, Inst, RC )
    IF ( RC /= HCO_SUCCESS ) THEN
       WRITE(MSG,*) 'Cannot find SeaFlux instance Nr. ', ExtState%SeaFlux
       CALL HCO_ERROR(HcoState%Config%Err,MSG,RC)
       RETURN
    ENDIF

    ! ---------------------------------------------------------------
    ! Calculate emissions
    ! ---------------------------------------------------------------

    ! Loop over all model species
    DO OcID = 1, Inst%nOcSpc

       ! Get HEMCO species ID
       HcoID = Inst%OcSpecs(OcID)%HcoID

       ! Skip this species if it has no corresponding HEMCO and/or
       ! model species
       IF ( HcoID                     < 0 ) CYCLE
       IF ( HcoState%Spc(HcoID)%ModID < 0 ) CYCLE

       IF ( verbose ) THEN
          WRITE(MSG,'(A40,I5)') &
               'Calculate air-sea flux for HEMCO species', HcoID
          CALL HCO_MSG(HcoState%Config%Err,MSG)
          WRITE(MSG,*) 'Module species name: ', &
                        TRIM(Inst%OcSpecs(OcID)%OcSpcName)
          CALL HCO_MSG(HcoState%Config%Err,MSG)
       ENDIF

       ! Get seawater concentration of given compound (from HEMCO core).
       ContName = TRIM(Inst%OcSpecs(OcID)%OcDataName)
       CALL HCO_EvalFld ( HcoState, ContName, SeaConc, RC )
       IF ( RC /= HCO_SUCCESS ) RETURN

       ! Calculate oceanic source (kg/m2/s) as well as the deposition
       ! velocity (1/s).
       CALL Calc_SeaFlux ( HcoState, ExtState, Inst,    &
                           SOURCE,   SINK,     SeaConc, &
                           OcID,     HcoID,    RC       )
       IF ( RC /= HCO_SUCCESS ) RETURN

       ! Set flux in HEMCO object [kg/m2/s]
       CALL HCO_EmisAdd ( HcoState, SOURCE, HcoID, RC, ExtNr=Inst%ExtNr )
       IF ( RC /= HCO_SUCCESS ) THEN
          MSG = 'HCO_EmisAdd error: ' // TRIM(Inst%OcSpecs(OcID)%OcSpcName)
          CALL HCO_ERROR(HcoState%Config%Err,MSG, RC )
          RETURN
       ENDIF
       IF ( RC /= HCO_SUCCESS ) RETURN

       ! Set deposition velocity in HEMCO object [1/s]
       CALL HCO_DepvAdd ( HcoState, SINK, HcoID, RC )
       IF ( RC /= HCO_SUCCESS ) RETURN

       ! Free pointers
       !SeaConc => NULL()

       ! Eventually add to dry deposition diagnostics
       ContName = 'DRYDEP_VEL_' // TRIM(HcoState%Spc(HcoID)%SpcName)
       Arr2D    => SINK
       CALL Diagn_Update( HcoState,                 &
                          cName   = TRIM(ContName), &
                          Array2D = Arr2D,          &
                          COL     = -1,             &
                          RC      = RC              )
       Arr2D => NULL()
    ENDDO !SpcID

    ! Cleanup
    Inst => NULL()

    ! Leave w/ success
    CALL HCO_LEAVE( HcoState%Config%Err,RC )

  END SUBROUTINE HCOX_SeaFlux_Run
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Calc_SeaFlux
!
! !DESCRIPTION: Subroutine CALC\_SEAFLUX calculates oceanic emissions
! of the specified tracer using the parameterization described in
! Johnson, 2010.
!\\
!\\
! The net emission flux is given by F = - Kg ( Cg - Caq*H ). Here, we
! calculate the source term ( Kg * H * Caq ) in units of kg/m2/s as
! well as the deposition velocity Kg in m/s.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Calc_SeaFlux( HcoState, ExtState,           &
                           Inst,     SOURCE,   SINK,     &
                           SeaConc,  OcID,     HcoID, RC )
!
! !USES:
!
    USE Ocean_ToolBox_Mod,  ONLY : CALC_KG
    USE Henry_Mod,          ONLY : CALC_KH, CALC_HEFF
    USE HCO_CALC_MOD,       ONLY : HCO_CheckDepv
    USE HCO_GeoTools_Mod,   ONLY : HCO_LANDTYPE
!
! !INPUT PARAMETERS:
!
    INTEGER,         INTENT(IN   ) :: OcID                ! ocean species ID
    INTEGER,         INTENT(IN   ) :: HcoID               ! HEMCO species ID
    TYPE(HCO_State), POINTER       :: HcoState            ! Output obj
    TYPE(Ext_State), POINTER       :: ExtState
    TYPE(MyInst),    POINTER       :: Inst
!
! !OUTPUT PARAMETERS:
!
    REAL(hp),        INTENT(  OUT) :: SOURCE(HcoState%NX,HcoState%NY )
    REAL(hp),        INTENT(  OUT) :: SINK  (HcoState%NX,HcoState%NY )
!
! !INPUT/OUTPUT PARAMETERS:
!
    REAL(hp),        INTENT(INOUT) :: SeaConc(HcoState%NX,HcoState%NY )
    INTEGER,         INTENT(INOUT) :: RC                 ! Error stat

!
! !REMARKS:
!  For now, the salinity and pH of seawater are prescribed to 35ppt and 8.0,
!  respectively.  The oceanic flux is not expected to be sensitive to these
!  parameters (which have only little variations anyway), but we may use
!  climatologies for these parameter at some point nevertheless!
!
! !REVISION HISTORY:
!  16 Apr 2013 - C. Keller   - Initial version
!  15 Aug 2014 - C. Keller   - Now restrict calculations to temperatures above
!                              10 deg C.
!  03 Oct 2014 - C. Keller   - Added surface temperature limit of 45 degrees C
!                              to avoid negative Schmidt numbers.
!  07 Oct 2014 - C. Keller   - Now use skin temperature instead of air temperature
!  06 Mar 2015 - C. Keller   - Now calculate deposition rate over entire PBL.
!  14 Oct 2015 - R. Yantosca - Pulled variables MW, VB, SCW out of the parallel
!                              loop.
!  06 Nov 2015 - C. Keller   - Now use HCO_LANDTYPE instead of FRCLND
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER             :: I, J, L, N
    REAL*8              :: IJSRC
    INTEGER             :: SCW
    REAL*8              :: P, V, VB, MW, KG
    REAL*8              :: K0, CR, PKA
    REAL*8              :: KH, HEFF
    REAL*8              :: TK, TC
    REAL(hp)            :: DEP_HEIGHT
    INTEGER             :: OLDWARN
    INTEGER             :: PBL_MAX
    INTEGER, SAVE       :: WARN = 0

    ! For now, hardcode salinity
    REAL(dp), PARAMETER :: S = 35.0_dp

    ! Set seawater PH to constant value of 8
    REAL(dp), PARAMETER :: PH = 8.0_dp

    ! Maximum allowed temperature (to avoid neg. Schmidt number)
    ! Set to 45 C (= 318.15 K)
    REAL(dp), PARAMETER :: TMAX = 318.15_dp

    ! Error handling
    CHARACTER(LEN=255)  :: MSG

    !=================================================================
    ! CALC_SEAFLUX begins here!
    !=================================================================

    ! Enter
    CALL HCO_ENTER( HcoState%Config%Err, 'Calc_SeaFlux (hcox_seaflux_mod.F90)', RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Init
    SOURCE  = 0.0_hp
    SINK    = 0.0_hp

    ! Extract Henry coefficients
    K0      = HcoState%Spc(HcoID)%HenryK0
    CR      = HcoState%Spc(HcoID)%HenryCR
    PKA     = HcoState%Spc(HcoID)%HenryPKA

    ! molecular weight [g/mol]
    ! Use real species molecular weight and not the emitted
    ! molecular weight. The molecular weight is only needed to
    ! calculate the air-side Schmidt number, which should be
    ! using the actual species MW.
    MW      = HcoState%Spc(HcoID)%MW_g

    ! Liquid molar volume at boiling point [cm3/mol]
    VB      = Inst%OcSpecs(OcID)%LiqVol

    ! Get parameterization type for Schmidt number in water
    SCW     = Inst%OcSpecs(OcID)%SCWPAR

    ! Model surface layer
    L       = 1

    ! Write out original warning status
    OLDWARN = WARN

    ! Loop over all grid boxes. Only emit into lowest layer

!$OMP PARALLEL DO                                     &
!$OMP DEFAULT( SHARED )                               &
!$OMP PRIVATE( I,  J,     N,        TK,        TC   ) &
!$OMP PRIVATE( P,  V,     KH,       RC,        HEFF ) &
!$OMP PRIVATE( KG, IJSRC, PBL_MAX,  DEP_HEIGHT      ) &
!$OMP SCHEDULE( DYNAMIC )

    DO J = 1, HcoState%NY
    DO I = 1, HcoState%NX

       ! Make sure we have no negative seawater concentrations
       IF ( SeaConc(I,J) < 0.0_sp ) SeaConc(I,J) = 0.0_sp

       ! Assume no air-sea exchange over snow/ice (ALBEDO > 0.4)
       IF ( ExtState%ALBD%Arr%Val(I,J) > 0.4_hp ) CYCLE

       ! Do only over the ocean:
       IF ( HCO_LANDTYPE( ExtState%WLI%Arr%Val(I,J), &
                          ExtState%ALBD%Arr%Val(I,J) ) == 0 ) THEN

          !-----------------------------------------------------------
          ! Get grid box and species specific quantities
          !-----------------------------------------------------------

          ! skin surface temp in K
          TK = ExtState%TSKIN%Arr%Val(I,J)

          ! Error check: the Schmidt number may become negative for
          ! very high temperatures - hence cap temperature at specified
          ! limit
          IF ( TK > TMAX ) THEN
             WARN = 1
             TK   = TMAX
          ENDIF

          ! Temperature in C
          TC = TK - 273.15d0

          ! Assume no air-sea exchange for temperatures below -10 deg C.
          ! This is rather arbitrary, but seawater should be frozen at
          ! that temperature anyways. Also, this ensures that the cal-
          ! culation of KG doesn't produce an overflow error, which occurs
          ! at temperatures of -10.7 to -10.9 deg C.
          IF ( TC < -10.0d0 ) CYCLE

          ! surface pressure [Pa]
          P = HcoState%Grid%PEDGE%Val(I,J,L)

          ! 10-m wind speed [m/s]
          V = ExtState%U10M%Arr%Val(I,J)**2 + &
              ExtState%V10M%Arr%Val(I,J)**2
          V = SQRT(V)

          ! Henry gas over liquid dimensionless constant and
          ! effective Henry constant [both unitless].
          CALL CALC_KH ( K0, CR, TK, KH, RC )  ! liquid over gas
          ! Exit here if error. Use error flags from henry_mod.F!
          IF ( RC /= 0 ) THEN
             RC  = HCO_FAIL
             WRITE(MSG,*) 'Cannot calculate KH: ', K0, CR, TK
             EXIT
          ENDIF
          CALL CALC_HEFF ( PKA, PH, KH, HEFF, RC )  ! liquid over gas
          ! Exit here if error. Use error flags from henry_mod.F!
          IF ( RC /= 0 ) THEN
             RC  = HCO_FAIL
             WRITE(MSG,*) 'Cannot calculate HEFF: ', PKA, PH, KH
             EXIT
          ENDIF

          ! Gas over liquid
          KH   = 1d0 / KH
          HEFF = 1d0 / HEFF

          !-----------------------------------------------------------
          ! Calculate exchange velocity KG in [m s-1]
          !-----------------------------------------------------------

          ! Get exchange velocity KG (m/s) following Johnson, 2010.
          ! Kg is defined as 1 / (1/k_air + H/k_water). Note that Kg
          ! is denoted Ka in Johnson, 2010!
          ! Use effective Henry constant here to account for
          ! hydrolysis!
          CALL CALC_KG( TC, P, V, S, HEFF, VB, MW, SCW, KG, RC )
          IF ( RC /= 0 ) THEN
             RC = HCO_FAIL
             WRITE(MSG,*) 'Cannot calculate KG: ', TC, P, V, S, HEFF
             EXIT
          ENDIF

          !-----------------------------------------------------------
          ! Calculate flux from the ocean (kg m-2 s-1):
          !-----------------------------------------------------------

          ! Fwa = KG * Cwater * H (Liss and Slater, 1974)
          ! Oceanic concentration is im [kg m-3], H is
          ! dimensionless, and KG is [m s-1], so IJSRC is
          ! [kg m-2 s-1].
          ! OcArr already accounts for pH effects, so apply the
          ! 'regular' Henry constant H here.
          IJSRC = KG * KH * SeaConc(I,J)

          ! Pass to flux array
          SOURCE(I,J) = IJSRC

          !-----------------------------------------------------------
          ! Calculate deposition rate to the ocean (s-1):
          !-----------------------------------------------------------

          ! Determine deposition height based on HEMCO option regarding
          ! the deposition length scale.
          IF ( HcoState%Options%PBL_DRYDEP ) THEN
             DO N = HcoState%NZ, 1, -1
                IF ( ExtState%FRAC_OF_PBL%Arr%Val(I,J,N) > 0.0_hp ) THEN
                   PBL_MAX = N
                   EXIT
                ENDIF
             ENDDO
          ELSE
             PBL_MAX = 1
          ENDIF
          DEP_HEIGHT = SUM(HcoState%Grid%BXHEIGHT_M%Val(I,J,1:PBL_MAX))

          ! Now calculate deposition rate from velocity and deposition
          ! height: [s-1] = [m s-1] / [m].
          SINK(I,J) = KG / DEP_HEIGHT

          ! Check validity of value
          CALL HCO_CheckDepv( HcoState, SINK(I,J), RC )

       ENDIF !Over ocean
    ENDDO !I
    ENDDO !J
!$OMP END PARALLEL DO


    ! Check exit status
    IF ( RC /= HCO_SUCCESS ) THEN
       CALL HCO_ERROR(HcoState%Config%Err,MSG, RC )
       RETURN
    ENDIF

    ! Warning?
    IF ( WARN /= OLDWARN ) THEN
       WRITE(MSG,*) 'Temperature limited to ', TMAX, 'K'
       CALL HCO_WARNING(HcoState%Config%Err, MSG, RC )
    ENDIF

    ! Leave w/ success
    CALL HCO_LEAVE( HcoState%Config%Err,RC )

  END SUBROUTINE Calc_SeaFlux
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCOX_SeaFlux_Init
!
! !DESCRIPTION: Subroutine HCOX\_SeaFlux\_Init initializes all module
! variables, including all species - specific parameter such as
! the liquid molar volume (Vb), the parameterization type for the
! Schmidt number in water (SCWPAR) and the name of the field containing
! oceanic concentrations.
!\\
!\\
! LiqVol is the liquid molar volume [cm3/mol]. If not stated otherwise,
! it is calculated using the Schroeder additive method as described in
! Johnson, 2010. Note that experimental values for LiqVol should be used
! if available!
!\\
!\\
! Table 3 of Johnson, 2010: Schroeder additive method for calculating
! Vb. For all atoms/structural items a molecule contains, the sum of the
! incre- ments will give the molar volume. e.g. CH2=CH2 contains 2 car-
! bon atoms, 4 hydrogen atoms and 1 double bond so the Schroeder
! Vb is 2x7 + 4x7 + 7 = 49cm3mol-1. * applies to all kinds of cyclic features
! and is applied only once to ring-containing compounds irrespective
! of the number of rings present.
!
! \begin{itemize}
! \item Atom/feature Increment/cm3mole-1
! \item Carbon       7.0
! \item Hydrogen     7.0
! \item Oxygen       7.0
! \item Nitrogen     7.0
! \item Bromine     31.5
! \item Chlorine    24.5
! \item Fluorine    10.5
! \item Iodine      38.5
! \item Sulfur      21.0
! \item Ring*       -7.0
! \item Double bond  7.0
! \item Triple bond 14.0
! \end{itemize}
!
! SCWPAR denotes which parameterization will be used to calculate the
! Schmidt number in water (in ocean\_toolbox\_mod). The following
! parameterizations are currently supported:
!
! \begin{enumerate}
! \item Parameterization as in Johnson, 2010 (default).
! \item Parameterization for DMS according to Saltzman et al., 1993.
! \item Parameterization for Acetone as in former acetone\_mod.F in GC.
! \item Parameterization for Acetaldehyde as in ald2\_mod.F from D. Millet
! \item Parameterization for MENO3, ETNO3 as in Fisher et al., 2018
! \end{enumerate}

! The oceanic surface concentrations of all species are obtained from
! external fields. These field names are specified in array OcDataName.
! For now, we obtain these concentrations from netCDF-files through the
! HEMCO core module, i.e. for each species there need to be a
! corresponding seawater concentration data file specified in the HEMCO
! configuration file. Once we use a coupled (ESMF) system, these names
! may be used to refer to the names of the concentration fields imported
! from the ocean model component.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCOX_SeaFlux_Init( HcoState, ExtName, ExtState, RC )
!
! !USES:
!
    USE HCO_ExtList_Mod,        ONLY : GetExtNr
    USE HCO_STATE_MOD,          ONLY : HCO_GetExtHcoID
!
! !INPUT PARAMETERS:
!
    TYPE(HCO_State),  POINTER        :: HcoState     ! Hemco State obj.
    CHARACTER(LEN=*), INTENT(IN   )  :: ExtName      ! Extension name
    TYPE(Ext_State),  POINTER        :: ExtState       ! Ext. obj.
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(INOUT)  :: RC           ! Return status
!
! !REVISION HISTORY:
!  16 Apr 2013 - C. Keller - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES
!
    ! Scalars
    TYPE(MyInst), POINTER          :: Inst
    INTEGER                        :: ExtNr, I, J, nSpc
    CHARACTER(LEN=255)             :: NAME_OC, MSG, ERR

    ! Arrays
    INTEGER,           ALLOCATABLE :: HcoIDs(:)
    CHARACTER(LEN=31), ALLOCATABLE :: SpcNames(:)

    !=================================================================
    ! HCOX_SeaFlux_Init begins here!
    !=================================================================

    ! Extension Nr.
    ExtNr = GetExtNr( HcoState%Config%ExtList, TRIM(ExtName) )
    IF ( ExtNr <= 0 ) RETURN

    ! Enter
    CALL HCO_ENTER( HcoState%Config%Err, 'HCOX_SeaFlux_Init (hcox_seaflux_mod.F90)', RC )
    IF ( RC /= HCO_SUCCESS ) RETURN
    ERR = 'nOcSpc too low!'

    ! Create instance for this simulation
    CALL InstCreate ( ExtNr, ExtState%SeaFlux, Inst, RC )
    IF ( RC /= HCO_SUCCESS ) THEN
       CALL HCO_ERROR ( HcoState%Config%Err, 'Cannot create SeaFlux instance', RC )
       RETURN
    ENDIF

    ! Verbose mode
    IF ( HcoState%amIRoot ) THEN
       MSG = 'Use air-sea flux emissions (extension module)'
       CALL HCO_MSG(HcoState%Config%Err,MSG,SEP1='-' )
       MSG = '   - Use species:'
       CALL HCO_MSG(HcoState%Config%Err,MSG )
    ENDIF

    ! ----------------------------------------------------------------------
    ! Get species IDs and settings
    ! ----------------------------------------------------------------------

    ! # of species for which air-sea exchange will be calculated
    Inst%nOcSpc = 7  ! updated to include MENO3, ETNO3, MOH

    ! Initialize vector w/ species information
    ALLOCATE ( Inst%OcSpecs(Inst%nOcSpc) )
    DO I = 1, Inst%nOcSpc
       Inst%OcSpecs(I)%HcoID      = -1
       Inst%OcSpecs(I)%OcSpcName  = ''
       Inst%OcSpecs(I)%OcDataName = ''
       Inst%OcSpecs(I)%LiqVol     = 0d0
       Inst%OcSpecs(I)%SCWPAR     = 1
    ENDDO

    ! Counter
    I = 0

    ! ----------------------------------------------------------------------
    ! CH3I:
    ! ----------------------------------------------------------------------

    I = I + 1
    IF ( I > Inst%nOcSpc ) THEN
       CALL HCO_ERROR ( HcoState%Config%Err, ERR, RC )
       RETURN
    ENDIF

    Inst%OcSpecs(I)%OcSpcName  = 'CH3I'
    Inst%OcSpecs(I)%OcDataName = 'CH3I_SEAWATER'
    Inst%OcSpecs(I)%LiqVol     = 1d0*7d0 + 3d0*7d0 + 1d0*38.5d0 ! Johnson, 2010
    Inst%OcSpecs(I)%SCWPAR     = 1 ! Schmidt number following Johnson, 2010

    ! ----------------------------------------------------------------------
    ! DMS:
    ! ----------------------------------------------------------------------

    I = I + 1
    IF ( I > Inst%nOcSpc ) THEN
       CALL HCO_ERROR ( HcoState%Config%Err, ERR, RC )
       RETURN
    ENDIF

    Inst%OcSpecs(I)%OcSpcName  = 'DMS'
    Inst%OcSpecs(I)%OcDataName = 'DMS_SEAWATER'
    Inst%OcSpecs(I)%LiqVol     = 2d0*7d0 + 6d0*7d0 + 1d0*21.0d0 ! Johnson, 2010
    Inst%OcSpecs(I)%SCWPAR     = 2 ! Schmidt number following Saltzman et al., 1993

    ! ----------------------------------------------------------------------
    ! Acetone:
    ! ----------------------------------------------------------------------

    I = I + 1
    IF ( I > Inst%nOcSpc ) THEN
       CALL HCO_ERROR ( HcoState%Config%Err, ERR, RC )
       RETURN
    ENDIF

    Inst%OcSpecs(I)%OcSpcName  = 'ACET'
    Inst%OcSpecs(I)%OcDataName = 'ACET_SEAWATER'
    Inst%OcSpecs(I)%LiqVol     = 3d0*7d0 + 6d0*7d0 + 1d0*7d0 + 1d0*7d0 ! Johnson, 2010
    Inst%OcSpecs(I)%SCWPAR     = 3 ! Schmidt number of acetone

    ! ----------------------------------------------------------------------
    ! Methanol:
    ! ----------------------------------------------------------------------

    I = I + 1
    IF ( I > Inst%nOcSpc ) THEN
       CALL HCO_ERROR ( HcoState%Config%Err, ERR, RC )
       RETURN
    ENDIF

    Inst%OcSpecs(I)%OcSpcName  = 'MOH'
    Inst%OcSpecs(I)%OcDataName = 'MOH_SEAWATER'
    Inst%OcSpecs(I)%LiqVol     = 1d0*7d0 + 4d0*7d0 + 1d0*7d0 ! Johnson, 2010
    Inst%OcSpecs(I)%SCWPAR     = 1 ! Schmidt number of methanol

    ! ----------------------------------------------------------------------
    ! Acetaldehyde:
    ! ----------------------------------------------------------------------

    I = I + 1
    IF ( I > Inst%nOcSpc ) THEN
       CALL HCO_ERROR ( ERR, RC )
       RETURN
    ENDIF

    Inst%OcSpecs(I)%OcSpcName  = 'ALD2'
    Inst%OcSpecs(I)%OcDataName = 'ALD2_SEAWATER'
    Inst%OcSpecs(I)%LiqVol     = 2d0*7d0 + 4d0*7d0 + 1d0*7d0 + 1d0*7d0 ! Johnson, 2010
    Inst%OcSpecs(I)%SCWPAR     = 4 ! Schmidt number of acetaldehyde

    ! ----------------------------------------------------------------------
    ! Methyl nitrate:
    ! ----------------------------------------------------------------------

    I = I + 1
    IF ( I > Inst%nOcSpc ) THEN
       CALL HCO_ERROR ( ERR, RC )
       RETURN
    ENDIF

    Inst%OcSpecs(I)%OcSpcName  = 'MENO3'
    Inst%OcSpecs(I)%OcDataName = 'MENO3_SEAWATER'
    Inst%OcSpecs(I)%LiqVol     = 64d0 ! Kornilov & Klselev 2015
    Inst%OcSpecs(I)%SCWPAR     = 1

    ! ----------------------------------------------------------------------
    ! Ethyl nitrate:
    ! ----------------------------------------------------------------------

    I = I + 1
    IF ( I > Inst%nOcSpc ) THEN
       CALL HCO_ERROR ( ERR, RC )
       RETURN
    ENDIF

    Inst%OcSpecs(I)%OcSpcName  = 'ETNO3'
    Inst%OcSpecs(I)%OcDataName = 'ETNO3_SEAWATER'
    Inst%OcSpecs(I)%LiqVol     = 82.2d0 ! Kornilov & Klselev 2015
    Inst%OcSpecs(I)%SCWPAR     = 1

    ! ----------------------------------------------------------------------
    ! Match module species with species assigned to this module in config.
    ! file
    ! ----------------------------------------------------------------------

    ! HEMCO species IDs of species names defined in config. file
    CALL HCO_GetExtHcoID( HcoState, ExtNr, HcoIDs, SpcNames, nSpc, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Set information in module variables
    DO I = 1, Inst%nOcSpc

       ! Append ocean tag '__OC' to this species name to make sure
       ! that we will also register non-tagged species.
       NAME_OC = TRIM(Inst%OcSpecs(I)%OcSpcName) // '__OC'

       DO J = 1, nSpc

          ! Compare model species names against defined module species.
          ! Also accept species names without the tag __OC, e.g.
          ! 'ACET' only instead of 'ACET__OC'.
          IF ( TRIM(SpcNames(J)) == TRIM(Inst%OcSpecs(I)%OcSpcName) .OR. &
               TRIM(SpcNames(J)) == TRIM(NAME_OC)              ) THEN
             Inst%OcSpecs(I)%HcoID = HcoIDs(J)
             EXIT
          ENDIF
       ENDDO !J

       ! verbose
       IF ( Inst%OcSpecs(I)%HcoID > 0 .AND. HcoState%amIRoot ) THEN
          WRITE(MSG,*) '   - ', &
               TRIM(Inst%OcSpecs(I)%OcSpcName), Inst%OcSpecs(I)%HcoID
          CALL HCO_MSG(HcoState%Config%Err,MSG)
       ENDIF
    ENDDO !I

    ! Set met fields
    ExtState%U10M%DoUse        = .TRUE.
    ExtState%V10M%DoUse        = .TRUE.
    ExtState%TSKIN%DoUse       = .TRUE.
    ExtState%ALBD%DoUse        = .TRUE.
    ExtState%WLI%DoUse         = .TRUE.
    IF ( HcoState%Options%PBL_DRYDEP ) THEN
       ExtState%FRAC_OF_PBL%DoUse = .TRUE.
    ENDIF
!    ExtState%FRCLND%DoUse      = .TRUE.

    ! Enable extensions
    !ExtState%SeaFlux = .TRUE.

    ! Return w/ success
    IF ( ALLOCATED(HcoIDs  ) ) DEALLOCATE(HcoIDs  )
    IF ( ALLOCATED(SpcNames) ) DEALLOCATE(SpcNames)
    CALL HCO_LEAVE( HcoState%Config%Err,RC )

  END SUBROUTINE HCOX_SeaFlux_Init
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCOX_SeaFlux_Final
!
! !DESCRIPTION: Subroutine HCOX\_SeaFlux\_Final deallocates
!  all module arrays.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCOX_SeaFlux_Final( ExtState )
!
! !INPUT PARAMETERS:
!
    TYPE(Ext_State),  POINTER       :: ExtState   ! Module options
!
! !REVISION HISTORY:
!  16 Apr 2013 - C. Keller - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
    !=================================================================
    ! HCOX_SeaFlux_Final begins here!
    !=================================================================
    CALL InstRemove( ExtState%SeaFlux )

    !IF ( ASSOCIATED( OcSpecs )) DEALLOCATE( OcSpecs )

  END SUBROUTINE HCOX_SeaFlux_Final
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
          IF ( ASSOCIATED( Inst%OcSpecs )) DEALLOCATE( Inst%OcSpecs )
          PrevInst%NextInst => Inst%NextInst
       ELSE
          AllInst => Inst%NextInst
       ENDIF
       DEALLOCATE(Inst)
       Inst => NULL()
    ENDIF

   END SUBROUTINE InstRemove
!EOC
END MODULE HCOX_SeaFlux_Mod
!EOM
