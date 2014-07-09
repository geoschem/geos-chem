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

  ! Variables carrying information about ocean species 
  INTEGER                     :: ExtNr
  INTEGER                     :: nOcSpc            ! # of ocean species
  TYPE(OcSpec), POINTER       :: OcSpecs(:)

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
  SUBROUTINE HCOX_SeaFlux_Run( am_I_Root, ExtState, HcoState, RC )
!
! !USES:
!
    USE HCO_FLUXARR_MOD,  ONLY : HCO_EmisAdd
    USE HCO_FLUXARR_MOD,  ONLY : HCO_DepvAdd
    USE HCO_EMISLIST_MOD, ONLY : EmisList_GetDataArr 
!
! !INPUT PARAMETERS:
!
    LOGICAL,         INTENT(IN   ) :: am_I_Root  ! root CPU?
    TYPE(HCO_State), POINTER       :: HcoState   ! Output obj
    TYPE(Ext_State), POINTER       :: ExtState  ! Module options  
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,         INTENT(INOUT) :: RC         ! Success or failure?
!
! !REVISION HISTORY: 
!  16 Apr 2013 - C. Keller - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER            :: OcID, HcoID
    REAL(hp), TARGET   :: SOURCE(HcoState%NX,HcoState%NY) 
    REAL(hp), TARGET   :: SINK  (HcoState%NX,HcoState%NY)
    CHARACTER(LEN=255) :: ContName
    CHARACTER(LEN=255) :: MSG
    LOGICAL            :: VERBOSE

    ! Pointers
    REAL(hp), POINTER  :: Arr2D(:,:)   => NULL() 
    REAL(hp), POINTER  :: SeaConc(:,:) => NULL()

    !=================================================================
    ! HCOX_SeaFlux_Run begins here!
    !=================================================================

    ! Enter
    CALL HCO_ENTER( 'HCOX_SeaFlux_Run (hcox_seaflux_mod.F90)', RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Return if extension disabled 
    IF ( .NOT. ExtState%SeaFlux ) RETURN

    ! Verbose?
    VERBOSE = am_I_Root .AND. HCO_VERBOSE_CHECK() 

    ! ---------------------------------------------------------------
    ! Calculate emissions
    ! ---------------------------------------------------------------

    ! Loop over all model species 
    DO OcID = 1, nOcSpc

       ! Get HEMCO species ID 
       HcoID = OcSpecs(OcID)%HcoID

       ! Skip this species if it has no corresponding HEMCO and/or
       ! model species 
       IF ( HcoID                         < 0 ) CYCLE
       IF ( HcoState%Spc(HcoID)%ModID < 0 ) CYCLE

       IF ( verbose ) THEN
          WRITE(MSG,'(A40,I5)') & 
               'Calculate air-sea flux for HEMCO species', HcoID
          CALL HCO_MSG(MSG)
          WRITE(MSG,*) 'Module species name: ', &
                        TRIM(OcSpecs(OcID)%OcSpcName)
          CALL HCO_MSG(MSG)
       ENDIF

       ! Get seawater concentration of given compound (from HEMCO core).
       ContName = TRIM(OcSpecs(OcID)%OcDataName)
       CALL EmisList_GetDataArr ( am_I_Root, ContName, SeaConc, RC )
       IF ( RC /= HCO_SUCCESS ) RETURN

       ! Calculate oceanic source (kg/m2/s) as well as the deposition 
       ! velocity (m/s).
       CALL Calc_SeaFlux ( am_I_Root, HcoState, ExtState, &
                           SOURCE,    SINK,     SeaConc,   &
                           OcID,      HcoID,    RC          )
       IF ( RC /= HCO_SUCCESS ) RETURN

       ! Set flux in HEMCO object [kg/m2/s]
       CALL HCO_EmisAdd ( HcoState, SOURCE, HcoID, RC )
       IF ( RC /= HCO_SUCCESS ) RETURN
      
       ! Set deposition velocity in HEMCO object [m/s]
       CALL HCO_DepvAdd ( HcoState, SINK, HcoID, RC )  
       IF ( RC /= HCO_SUCCESS ) RETURN
      
       ! Free pointers
       SeaConc => NULL()

       ! Eventually update diagnostics
       IF ( Diagn_AutoFillLevelDefined(2) ) THEN
          Arr2D => SOURCE 
          CALL Diagn_Update( am_I_Root, HcoState, ExtNr=ExtNr, &
                             Cat=-1, Hier=-1, HcoID=HcoID,     &
                             AutoFill=1, Array2D=Arr2D, RC=RC   )
          IF ( RC /= HCO_SUCCESS ) RETURN 
          Arr2D => NULL() 
       ENDIF

       ! TODO: update depostion diagnostics
 
    ENDDO !SpcID

    ! Leave w/ success
    CALL HCO_LEAVE ( RC ) 

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
  SUBROUTINE Calc_SeaFlux( am_I_Root, HcoState, ExtState, & 
                           SOURCE,    SINK,     SeaConc,   &
                           OcID,      HcoID,    RC          )
!
! !USES:
! 
    USE Ocean_ToolBox_Mod,  ONLY : CALC_KG
    USE Henry_Mod,          ONLY : CALC_KH, CALC_HEFF
!
! !INPUT PARAMETERS:
!
    LOGICAL,         INTENT(IN   ) :: am_I_Root           ! root CPU?
    INTEGER,         INTENT(IN   ) :: OcID                ! ocean species ID 
    INTEGER,         INTENT(IN   ) :: HcoID               ! HEMCO species ID
    TYPE(HCO_State), POINTER       :: HcoState            ! Output obj
    TYPE(Ext_State), POINTER       :: ExtState     
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
!  16 Apr 2013 - C. Keller - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER            :: I, J, L
    REAL*8             :: IJSRC, IJSINK
    INTEGER            :: SCW
    REAL*8             :: P, V, S, VB, MW, KG
    REAL*8             :: K0, CR, PKA
    REAL*8             :: KH, HEFF, PH
    REAL*8             :: TK, TC
                       
    ! Error handling   
    CHARACTER(LEN=255) :: MSG

    !=================================================================
    ! CALC_SEAFLUX begins here!
    !=================================================================

    ! Enter
    CALL HCO_ENTER ( 'Calc_SeaFlux (hcox_seaflux_mod.F90)', RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Init
    SOURCE(:,:) = 0d0
    SINK  (:,:) = 0d0

    ! Extract Henry coefficients
    K0  = HcoState%Spc(HcoID)%HenryK0
    CR  = HcoState%Spc(HcoID)%HenryCR
    PKA = HcoState%Spc(HcoID)%HenryPKA

    ! Model surface layer
    L = 1

    ! Loop over all grid boxes. Only emit into lowest layer

!$OMP PARALLEL DO                                                   &
!$OMP DEFAULT( SHARED )                                             &
!$OMP PRIVATE( I,           J,        PH,       TK                ) &
!$OMP PRIVATE( TC,          P,        MW,       VB,     S         ) &
!$OMP PRIVATE( V,           KH,       RC,       HEFF,   SCW       ) &
!$OMP PRIVATE( KG,          IJSRC                                 ) &
!$OMP SCHEDULE( DYNAMIC )

    DO J = 1, HcoState%NY
    DO I = 1, HcoState%NX
      
       ! Make sure we have no negative seawater concentrations 
       IF ( SeaConc(I,J) < 0d0 ) SeaConc(I,J) = 0d0

       ! Do only over the ocean, i.e. if land fraction is less
       ! than 0.8
       IF ( ExtState%FRCLND%Arr%Val(I,J) < 0.8 ) THEN

          !-----------------------------------------------------------
          ! Get grid box and species specific quantities
          !-----------------------------------------------------------

          ! pH of sea water. For the moment, set to 8
          PH = 8d0
 
          ! surface air temp in K and C
          TK = ExtState%TSURFK%Arr%Val(I,J) 
          TC = TK - 273.15d0
 
          ! surface pressure [Pa]
          P = ExtState%PSURF%Arr%Val(I,J) * 100d0

          ! molecular weight [g/mol]
          MW = HcoState%Spc(HcoID)%MW_g

          ! Liquid molar volume at boiling point [cm3/mol]
          VB = OcSpecs(OcID)%LiqVol

          ! Salinity [ppt]
          ! Set to constant value for now!
          S = 35d0 

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
             MSG = 'Cannot calculate KH'
             EXIT
          ENDIF
          CALL CALC_HEFF ( PKA, PH, KH, HEFF, RC )  ! liquid over gas 
          ! Exit here if error. Use error flags from henry_mod.F!
          IF ( RC /= 0 ) THEN
             RC  = HCO_FAIL
             MSG = 'Cannot calculate HEFF'
             EXIT
          ENDIF

          ! Gas over liquid
          KH   = 1d0 / KH
          HEFF = 1d0 / HEFF

          ! Get parameterization type for Schmidt number in water 
          SCW = OcSpecs(OcID)%SCWPAR

          !-----------------------------------------------------------
          ! Calculate exchange velocity KG in [m s-1]
          !-----------------------------------------------------------

          ! Assume no air-sea exchange over snow/ice (ALBEDO > 0.4)
          IF ( ExtState%ALBD%Arr%Val(I,J) > 0.4d0 ) THEN
             KG = 0d0

          ! Get exchange velocity KG (m/s) following Johnson, 2010.
          ! Kg is defined as 1 / (1/k_air + H/k_water). Note that Kg 
          ! is denoted Ka in Johnson, 2010!
          ! Use effective Henry constant here to account for
          ! hydrolysis!
          ELSE
             CALL CALC_KG( TC, P, V, S, HEFF, VB, MW, SCW, KG )
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
          ! Calculate deposition velocity to the ocean (m s-1):
          !-----------------------------------------------------------

          ! Pass to deposition array
          SINK(I,J) = KG 

       ENDIF !Over ocean
    ENDDO !I
    ENDDO !J
!$OMP END PARALLEL DO

    ! Check exit status
    IF ( RC /= HCO_SUCCESS ) THEN
       CALL HCO_ERROR ( MSG, RC )
       RETURN
    ENDIF

    ! Leave w/ success
    CALL HCO_LEAVE ( RC ) 

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
! \item Parameteriaztion for DMS according to Saltzman et al., 1993.
! \item Parameteriaztion for Acetone as in former acetone\_mod.F in GC. 
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
  SUBROUTINE HCOX_SeaFlux_Init( am_I_Root, HcoState, ExtName, ExtState, RC )
!
! !USES:
!
    USE HCOX_ExtList_Mod,       ONLY : GetExtNr, GetExtHcoID
!
! !INPUT PARAMETERS:
!
    LOGICAL,          INTENT(IN   )  :: am_I_Root    ! root CPU?
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
    INTEGER                        :: I, J, nSpc
    CHARACTER(LEN=255)             :: NAME_OC, MSG

    ! Arrays
    INTEGER,           ALLOCATABLE :: HcoIDs(:)
    CHARACTER(LEN=31), ALLOCATABLE :: SpcNames(:)

    !=================================================================
    ! HCOX_SeaFlux_Init begins here!
    !=================================================================

    ! Extension Nr.
    ExtNr = GetExtNr( TRIM(ExtName) )
    IF ( ExtNr <= 0 ) RETURN
 
    ! Enter 
    CALL HCO_ENTER (  'HCOX_SeaFlux_Init (hcox_seaflux_mod.F90)', RC )
    IF ( RC /= HCO_SUCCESS ) RETURN
    MSG = 'nOcSpc too low!'

    ! Verbose mode
    MSG = 'Use air-sea flux emissions (extension module)'
    CALL HCO_MSG( MSG,SEP1='-' )
    MSG = '   - Use species:'
    CALL HCO_MSG( MSG )

    ! ---------------------------------------------------------------------- 
    ! Get species IDs and settings 
    ! ---------------------------------------------------------------------- 
    
    ! # of species for which air-sea exchange will be calculated
    nOcSpc = 3

    ! Initialize vector w/ species information
    ALLOCATE ( OcSpecs(nOcSpc) ) 
    DO I = 1, nOcSpc
       OcSpecs(I)%HcoID      = -1
       OcSpecs(I)%OcSpcName  = ''
       OcSpecs(I)%OcDataName = ''
       OcSpecs(I)%LiqVol     = 0d0
       OcSpecs(I)%SCWPAR     = 1
    ENDDO

    ! Counter
    I = 0

    ! ----------------------------------------------------------------------
    ! CH3I:
    ! ----------------------------------------------------------------------

    I = I + 1
    IF ( I > nOcSpc ) THEN
       CALL HCO_ERROR ( MSG, RC )
       RETURN
    ENDIF

    OcSpecs(I)%OcSpcName  = 'CH3I'
    OcSpecs(I)%OcDataName = 'CH3I_SEAWATER'
    OcSpecs(I)%LiqVol     = 1d0*7d0 + 3d0*7d0 + 1d0*38.5d0 ! Johnson, 2010
    OcSpecs(I)%SCWPAR     = 1

    ! ----------------------------------------------------------------------
    ! DMS:
    ! ----------------------------------------------------------------------

    I = I + 1
    IF ( I > nOcSpc ) THEN
       CALL HCO_ERROR ( MSG, RC )
       RETURN
    ENDIF

    OcSpecs(I)%OcSpcName  = 'DMS'
    OcSpecs(I)%OcDataName = 'DMS_SEAWATER'
    OcSpecs(I)%LiqVol     = 2d0*7d0 + 6d0*7d0 + 1d0*21.0d0 ! Johnson, 2010
    OcSpecs(I)%SCWPAR     = 2

    ! ----------------------------------------------------------------------
    ! Acetone:
    ! ----------------------------------------------------------------------

    I = I + 1
    IF ( I > nOcSpc ) THEN
       CALL HCO_ERROR ( MSG, RC )
       RETURN
    ENDIF

    OcSpecs(I)%OcSpcName  = 'ACET'
    OcSpecs(I)%OcDataName = 'ACET_SEAWATER'
    OcSpecs(I)%LiqVol     = 3d0*7d0 + 6d0*7d0 + 1d0*7d0 + 1d0*7d0 ! Johnson, 2010
    OcSpecs(I)%SCWPAR     = 1

    ! ----------------------------------------------------------------------
    ! Match module species with species assigned to this module in config.
    ! file 
    ! ----------------------------------------------------------------------

    ! HEMCO species IDs of species names defined in config. file 
    CALL GetExtHcoID( HcoState, ExtNr, HcoIDs, SpcNames, nSpc, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Set information in module variables
    DO I = 1, nOcSpc 

       ! Append ocean tag '__OC' to this species name to make sure
       ! that we will also register non-tagged species.
       NAME_OC = TRIM(OcSpecs(I)%OcSpcName) // '__OC'

       DO J = 1, nSpc

          ! Compare model species names against defined module species. 
          ! Also accept species names without the tag __OC, e.g.
          ! 'ACET' only instead of 'ACET__OC'. 
          IF ( TRIM(SpcNames(J)) == TRIM(OcSpecs(I)%OcSpcName) .OR. &
               TRIM(SpcNames(J)) == TRIM(NAME_OC)              ) THEN
             OcSpecs(I)%HcoID = HcoIDs(J)
             EXIT
          ENDIF
       ENDDO !J

       ! verbose
       IF ( OcSpecs(I)%HcoID > 0 ) THEN
          WRITE(MSG,*) '   - ', &
               TRIM(OcSpecs(I)%OcSpcName), OcSpecs(I)%HcoID
          CALL HCO_MSG(MSG)
       ENDIF
    ENDDO !I

    ! Set met fields
    ExtState%U10M%DoUse   = .TRUE.
    ExtState%V10M%DoUse   = .TRUE.
    ExtState%PSURF%DoUse  = .TRUE.
    ExtState%TSURFK%DoUse = .TRUE.
    ExtState%ALBD%DoUse   = .TRUE.
    ExtState%FRCLND%DoUse = .TRUE.
    
    ! Enable extensions
    ExtState%SeaFlux = .TRUE.

    ! Return w/ success
    IF ( ALLOCATED(HcoIDs  ) ) DEALLOCATE(HcoIDs  )
    IF ( ALLOCATED(SpcNames) ) DEALLOCATE(SpcNames)
    CALL HCO_LEAVE ( RC )

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
  SUBROUTINE HCOX_SeaFlux_Final()
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

    IF ( ASSOCIATED( OcSpecs )) DEALLOCATE( OcSpecs ) 

  END SUBROUTINE HCOX_SeaFlux_Final
!EOC
END MODULE HCOX_SeaFlux_Mod
!EOM
