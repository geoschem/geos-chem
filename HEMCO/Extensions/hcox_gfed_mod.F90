!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: hcox_gfed_mod.F90
!
! !DESCRIPTION: Module HCOX\_GFED\_MOD contains routines to calculate
! GFED4 biomass burning emissions in HEMCO.
!
! !INTERFACE:
!
MODULE HCOX_GFED_MOD
!
! !USES:
!
  USE HCO_ERROR_MOD
  USE HCO_DIAGN_MOD
  USE HCOX_TOOLS_MOD
  USE HCO_STATE_MOD,  ONLY : HCO_State
  USE HCOX_State_MOD, ONLY : Ext_State

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC :: HCOX_GFED_Init
  PUBLIC :: HCOX_GFED_Run
  PUBLIC :: HCOX_GFED_Final
!
! !REMARKS:
!  Monthly emissions of DM are read from disk,
!  multiplied by daily and 3hourly fractions (if necessary), and then
!  multiplied by the appropriate emission factors to produce biomass
!  burning emissions.
!
! All species to be used must be listed in the settings section of the HEMCO
! configuration file. For every listed species, individual scale factors as
! well as masks can be defined. For example, to scale FINN CO emissions by a
! factor of 1.05 and restrict them to North America, as well as to scale NO
! emissions by a factor of 1.5:
!
!111     GFED              : on  NO/CO/ALK4/ACET/MEK/ALD2/PRPE/C3H8/CH2O/C2H6/SO2/NH3/BC/OC/GLYC/MGLY/BENZ/TOLU/XYLE/C2H4/C2H2/GLYC/CO2/CH4/HCOOH/DMS/ISOP/LIMO/MOH/EOH/ACTA/GLYX/HAC
!    --> GFED4             :       true
!    --> GFED_daily        :       false
!    --> GFED_3hourly      :       false
!    --> hydrophilic BC    :       0.2
!    --> hydrophilic OC    :       0.5
!    --> Mask_CO           :       NAMASK
!    --> Scaling_CO        :       1.05
!    --> Scaling_NO        :       1.5
!
! Field NAMASK must be defined in section mask of the HEMCO configuration file.
!
! For SOA_SVPOA mechanism:
! * If tracers POG1 and POG2 are specified, emissions are calculated from OC,
!   multiplied by a POG scale factor (Scaling_POG1, Scaling_POG2) that must be
!   specified in the HEMCO configuration file.
! * If tracer NAP is specified, emissions are calculated from CO, multiplied
!   by a NAP scale factor (Scaling_NAP) that must be specified in the HEMCO
!   configuration file.
!
!  References:
!  ============================================================================
!  (1 ) Original GFED3 database from Guido van der Werf
!        http://www.falw.vu/~gwerf/GFED/GFED3/emissions/
!  (2 ) Giglio, L., Randerson, J. T., van der Werf, G. R., Kasibhatla, P. S.,
!       Collatz, G. J., Morton, D. C., and DeFries, R. S.: Assessing
!       variability and long-term trends in burned area by merging multiple
!       satellite fire products, Biogeosciences, 7, 1171-1186,
!       doi:10.5194/bg-7-1171-2010, 2010.
!  (3 ) van der Werf, G. R., Randerson, J. T., Giglio, L., Collatz, G. J.,
!       Mu, M., Kasibhatla, P. S., Morton, D. C., DeFries, R. S., Jin, Y.,
!       and van Leeuwen, T. T.: Global fire emissions and the contribution of
!       deforestation, savanna, forest, agricultural, and peat fires
!       (1997â~@~S2009), Atmos. Chem. Phys., 10, 11707-11735,
!       doi:10.5194/acp-10-11707-2010, 2010.
!
! !REVISION HISTORY:
!  07 Sep 2011 - P. Kasibhatla - Initial version, based on GFED2
!  07 Sep 2011 - R. Yantosca   - Added ProTeX headers
!  14 Feb 2012 - M. Payer      - Add modifications for CH4 (K. Wecht)
!  01 Mar 2012 - R. Yantosca   - Now reference new grid_mod.F90
!  06 Mar 2012 - P. Kasibhatla - Final version
!  01 Aug 2012 - R. Yantosca - Add reference to findFreeLUN from inqure_mod.F90
!  03 Aug 2012 - R. Yantosca - Move calls to findFreeLUN out of DEVEL block
!  14 Mar 2013 - M. Payer    - Replace NOx emissions with NO emissions as part
!                              of removal of NOx-Ox partitioning
!  15 Dec 2013 - C. Keller   - Now a HEMCO extension. Emissions in kg/m2/s,
!                              emission factors in kg/kgDM.
!  01 Jul 2014 - R. Yantosca - Now use F90 free-format indentation
!  01 Jul 2014 - R. Yantosca - Cosmetic changes in ProTeX headers
!  08 Aug 2014 - R. Yantosca - Now avoid ASCII file reads for ESMF
!  23 Sep 2014 - C. Keller   - Increase N_SPEC to 26 (+Hg0)
!  12 Mar 2015 - C. Keller / P. Kasibhatla - Added GFED-4.
!  03 Jun 2015 - C. Keller / P. Kasibhatla - GFED-4 update: now use GFED-4
!                                            specific emission factors and DM data.
!  14 Oct 2016 - C. Keller    - Now use HCO_EvalFld instead of HCO_GetPtr.
!  11 Feb 2017 - S. Farina    - Increase N_SPEC to 27 (+SOAP)
!  23 Mar 2017 - M. Sulprizio - Increase N_SPEC to 29 (+EOH+MTPA)
!  29 Mar 2018 - K. Travis    - Update GFED4 emission factors, increase to 34 species
!  29 Mar 2018 - K. Travis    - Remove GFED3
!  12 Sep 2018 - C. Keller    - Added instance wrapper
!EOP
!------------------------------------------------------------------------------
!
! !DEFINED PARAMETERS:
!
  !=================================================================
  ! MODULE PARAMETERS
  !
  ! N_EMFAC : Number of emission factors per species
  ! N_SPEC  : Max. number of species
  !=================================================================
  INTEGER,           PARAMETER :: N_EMFAC = 6
  INTEGER,           PARAMETER :: N_SPEC  = 35 ! increase from 34 (v12.5.0 default)
                                               ! to 35 for MOH
!
! !PRIVATE TYPES:
!
  TYPE :: MyInst
   !=================================================================
   ! HEMCO VARIABLES
   !
   ! ExtNr   : Extension number
   ! DoDay   : TRUE if dialy scale factors are used
   ! Do3Hr   : TRUE if 3-hourly scale factors are used
   !=================================================================
   INTEGER                       :: Instance
   INTEGER                       :: ExtNr
   LOGICAL                       :: DoDay
   LOGICAL                       :: Do3Hr
   LOGICAL                       :: IsGFED4

   !=================================================================
   ! SPECIES VARIABLES
   !
   ! nSpc     : Number of GFED species (specified in config. file)
   ! SpcNames : Names of all used GFED species
   ! HcoIDs   : HEMCO species IDs of all used GFED species
   ! gfedIDs  : Index of used GFED species in scale factor table
   ! SpcScal  : Additional scaling factors assigned to species through
   !            the HEMCO configuration file (e.g. Scaling_CO).
   !=================================================================
   INTEGER                    :: nSpc
   CHARACTER(LEN=31), POINTER :: SpcNames(:) => NULL()
   CHARACTER(LEN=61), POINTER :: SpcScalFldNme(:) => NULL()
   INTEGER,           POINTER :: HcoIDs(:) => NULL()
   INTEGER,           POINTER :: GfedIDs(:) => NULL()
   REAL(sp),          POINTER :: SpcScal(:) => NULL()

   !=================================================================
   ! SCALE FACTORS
   !
   ! GFED_EMFAC:  emission scale factors for each species and
   !              emission factor type. The filename of the emissions
   !              emissions factor table is specified in the HEMCO
   !              configuration file. All scale factors in kg/kgDM.
   ! OCPIfrac   : Fraction of OC that converts into hydrophilic OC.
   !              Can be set in HEMCO configuration file (default=0.5)
   ! BCPIfrac   : Fraction of BC that converts into hydrophilic BC.
   !              Can be set in HEMCO configuration file (default=0.2)
   ! POG1frac   : Fraction of SVOC that is assigned to POG1.
   !              Can be set in HEMCO configuration file (default=0.49)
   !=================================================================
   REAL(hp), POINTER              :: GFED4_EMFAC(:,:) => NULL()
   REAL(hp), POINTER              :: GFED_EMFAC (:,:) => NULL()
   REAL(sp)                       :: OCPIfrac
   REAL(sp)                       :: BCPIfrac
   REAL(sp)                       :: POG1frac
   REAL(sp)                       :: SOAPfrac

   !=================================================================
   ! DATA ARRAY POINTERS
   !
   ! These are the pointers to the 6 input data specified in the
   ! the configuration file
   !=================================================================
   REAL(hp), POINTER           :: GFED_SAVA(:,:) => NULL()
   REAL(hp), POINTER           :: GFED_BORF(:,:) => NULL()
   REAL(hp), POINTER           :: GFED_TEMP(:,:) => NULL()
   REAL(hp), POINTER           :: GFED_DEFO(:,:) => NULL()
   REAL(hp), POINTER           :: GFED_PEAT(:,:) => NULL()
   REAL(hp), POINTER           :: GFED_AGRI(:,:) => NULL()
   REAL(hp), POINTER           :: DAYSCAL (:,:) => NULL()
   REAL(hp), POINTER           :: HRSCAL  (:,:) => NULL()

   TYPE(MyInst), POINTER       :: NextInst => NULL()
  END TYPE MyInst

  ! Pointer to instances
  TYPE(MyInst), POINTER        :: AllInst => NULL()

CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCOX_GFED_Run
!
! !DESCRIPTION: Subroutine HcoX\_GFED\_Run is the driver run routine to
! calculate seasalt emissions in HEMCO.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCOX_GFED_Run( ExtState, HcoState, RC )
!
! !USES:
!
    USE HCO_Calc_Mod,     ONLY : HCO_EvalFld
    USE HCO_EmisList_Mod, ONLY : HCO_GetPtr
    USE HCO_FluxArr_MOD,  ONLY : HCO_EmisAdd
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(HCO_State), POINTER        :: HcoState   ! Output obj
    TYPE(Ext_State), POINTER        :: ExtState  ! Module options
    INTEGER,         INTENT(INOUT)  :: RC         ! Success or failure?
!
! !REVISION HISTORY:
!  07 Sep 2011 - P. Kasibhatla - Initial version, based on GFED2
!  15 Dec 2013 - C. Keller     - Now a HEMCO extension
!  03 Apr 2015 - C. Keller     - Humid tropical forest mask is not binary
!                                any more but fraction (0.0 - 1.0).
!  21 Sep 2016 - R. Yantosca   - Bug fix: move WHERE statement for HUMTROP
!                                into the GFED3 block to avoid segfault
!  10 Mar 2017 - M. Sulprizio  - Add SpcArr3D for emitting 65% of biomass
!                                burning emissions into the PBL and 35% into the
!                                free troposphere, following code from E.Fischer
!  24 Apr 2017 - M. Sulprizio  - Comment out vertical distribution of biomass
!                                burning emissions for now.
!  12 May 2017 - M. Sulprizio  - Comment out partitioning of NO directly to PAN
!                                and HNO3 for now.
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    LOGICAL, SAVE       :: FIRST = .TRUE.
    INTEGER             :: N, M
    REAL(hp), POINTER   :: TmpPtr(:,:)
    CHARACTER(LEN=63)   :: MSG

    REAL(hp), TARGET    :: SpcArr(HcoState%NX,HcoState%NY)
    REAL(hp), TARGET    :: TypArr(HcoState%NX,HcoState%NY)

    TYPE(MyInst), POINTER :: Inst

!==============================================================================
! This code is required for the vertical distribution of biomass burning emiss.
! We will keep it here for a future implementation. (mps, 4/24/17)
!    INTEGER             :: I, J, L, N, M
!    INTEGER             :: PBL_MAX
!    REAL(hp)            :: PBL_FRAC, F_OF_PBL, F_OF_FT
!    REAL(hp)            :: DELTPRES, TOTPRESFT
!    REAL(hp), TARGET    :: SpcArr3D(HcoState%NX,HcoState%NY,HcoState%NZ)
!==============================================================================

    !=================================================================
    ! HCOX_GFED_Run begins here!
    !=================================================================

    ! Return if extension disabled
    IF ( ExtState%GFED <= 0 ) RETURN

    ! Enter
    CALL HCO_ENTER( HcoState%Config%Err, 'HCOX_GFED_Run (hcox_gfed_mod.F90)', RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Get instance
    Inst => NULL()
    CALL InstGet ( ExtState%GFED, Inst, RC )
    IF ( RC /= HCO_SUCCESS ) THEN
       WRITE(MSG,*) 'Cannot find GFED instance Nr. ', ExtState%GFED
       CALL HCO_ERROR(HcoState%Config%Err,MSG,RC)
       RETURN
    ENDIF

!==============================================================================
! This code is required for the vertical distribution of biomass burning emiss.
! We will keep it here for a future implementation. (mps, 4/24/17)
!    ! Add only 65% biomass burning source to boundary layer, the
!    ! rest is emitted into the free troposphere (mps from evf+tjb, 3/10/17)
!    PBL_FRAC = 0.65_hp
!==============================================================================

    !-----------------------------------------------------------------
    ! Get pointers to data arrays
    !-----------------------------------------------------------------
    !IF ( FIRST ) THEN

    IF ( Inst%IsGFED4 ) THEN
          CALL HCO_EvalFld ( HcoState, 'GFED_SAVA', Inst%GFED_SAVA, RC )
          IF ( RC /= HCO_SUCCESS ) RETURN
          CALL HCO_EvalFld ( HcoState, 'GFED_BORF', Inst%GFED_BORF, RC )
          IF ( RC /= HCO_SUCCESS ) RETURN
          CALL HCO_EvalFld ( HcoState, 'GFED_TEMP', Inst%GFED_TEMP, RC )
          IF ( RC /= HCO_SUCCESS ) RETURN
          CALL HCO_EvalFld ( HcoState, 'GFED_DEFO', Inst%GFED_DEFO, RC )
          IF ( RC /= HCO_SUCCESS ) RETURN
          CALL HCO_EvalFld ( HcoState, 'GFED_PEAT', Inst%GFED_PEAT, RC )
          IF ( RC /= HCO_SUCCESS ) RETURN
          CALL HCO_EvalFld ( HcoState, 'GFED_AGRI', Inst%GFED_AGRI, RC )
          IF ( RC /= HCO_SUCCESS ) RETURN
       ENDIF

       ! Also point to scale factors if needed
       IF ( Inst%DoDay ) THEN
          CALL HCO_EvalFld ( HcoState, 'GFED_FRAC_DAY', &
                             Inst%DAYSCAL,   RC )
          IF ( RC /= HCO_SUCCESS ) RETURN
       ENDIF
       IF ( Inst%Do3Hr ) THEN
          CALL HCO_EvalFld ( HcoState, 'GFED_FRAC_3HOUR', &
                             Inst%HRSCAL,    RC )
          IF ( RC /= HCO_SUCCESS ) RETURN
       ENDIF

       FIRST = .FALSE.
    !ENDIF

    !-----------------------------------------------------------------
    ! Calculate emissions for defined species
    !-----------------------------------------------------------------
    DO N = 1, Inst%nSpc

       ! Continue if species not defined
       IF ( Inst%HcoIDs(N)  < 0 ) CYCLE
       IF ( Inst%GfedIDs(N) < 0 ) CYCLE

       ! SpcArr are the total biomass burning emissions for this
       ! species. TypArr are the emissions from a given source type.
       SpcArr   = 0.0_hp
!==============================================================================
! This code is required for the vertical distribution of biomass burning emiss.
! We will keep it here for a future implementation. (mps, 4/24/17)
!       SpcArr3D = 0.0_hp
!==============================================================================

       ! Calculate emissions for all source types
       DO M = 1, N_EMFAC

          ! Point to the emission factor array for each source type
          SELECT CASE ( M )
             CASE( 1 )
                TMPPTR => Inst%GFED_SAVA
             CASE( 2 )
                TMPPTR => Inst%GFED_BORF
             CASE( 3 )
                TMPPTR => Inst%GFED_TEMP
             CASE( 4 )
                TMPPTR => Inst%GFED_DEFO
             CASE( 5 )
                TMPPTR => Inst%GFED_PEAT
             CASE( 6 )
                TMPPTR => Inst%GFED_AGRI
             CASE DEFAULT
                CALL HCO_ERROR ( HcoState%Config%Err, 'Undefined emission factor', RC )
                RETURN
          END SELECT

          ! Calculate emissions for this type. The emission factors
          ! per type are in kgDM/m2/s, and the GFED_EMFAC scale factors
          ! are in kg/kgDM (or kgC/kgDM for VOCs). This gives us TypArr
          ! in kg/m2/s.
          ! Use woodland emission factors for 'deforestation' outside
          ! humid tropical forest.
          ! Deforestation emissions now use the weighted sum of
          ! deforestation and woodland scale factors, based on the value
          ! of the humid tropical forest mask. This makes the calculation
          ! less dependent on model resolution. (ckeller, 4/3/15)
          TypArr = TmpPtr * Inst%GFED_EMFAC(Inst%GfedIDs(N),M)

          ! Eventually add daily / 3-hourly scale factors. These scale
          ! factors are unitless.
          IF ( Inst%DoDay ) THEN
             !IF ( ASSOCIATED(DAYSCAL) ) THEN
                TypArr = TypArr * Inst%DAYSCAL
             !ENDIF
          ENDIF
          IF ( Inst%Do3Hr ) THEN
             !IF ( ASSOCIATED(HRSCAL) ) THEN
                TypArr = TypArr * Inst%HRSCAL
             !ENDIF
          ENDIF

          ! Add to output array
          SpcArr = SpcArr + TypArr

          ! Nullify pointer
          TmpPtr  => NULL()

       ENDDO !M

       ! Apply species specific scale factors
       SpcArr = SpcArr * Inst%SpcScal(N)

       SELECT CASE ( Inst%SpcNames(N) )
          CASE ( 'OCPI' )
             SpcArr = SpcArr * Inst%OCPIfrac
          CASE ( 'OCPO' )
             SpcArr = SpcArr * (1.0_sp - Inst%OCPIfrac)
          CASE ( 'BCPI' )
             SpcArr = SpcArr * Inst%BCPIfrac
          CASE ( 'BCPO' )
             SpcArr = SpcArr * (1.0_sp - Inst%BCPIfrac)
          CASE ( 'POG1' )
             SpcArr = SpcArr * Inst%POG1frac
          CASE ( 'POG2' )
             SpcArr = SpcArr * (1.0_sp - Inst%POG1frac)
          CASE ( 'SOAP' )
             SpcArr = SpcArr * Inst%SOAPfrac
!==============================================================================
! This code is required for partitioning NOx emissions directly to PAN and HNO3.
! We will keep it here as an option for users focusing on North American fires.
! (mps, 5/12/17)
!          ! Put 40% of NOx Biomass emissions into PAN
!          ! and 20% into HNO3 (evf, 9/9/11, 9/15/11)
!          ! Sensitivity study with Hudman 2007 recommendation
!          ! of 80% of NOX as PAN. (evf, 4/25/12)
!          CASE ( 'NO' )
!             SpcArr = SpcArr * 0.40_sp
!          CASE ( 'PAN' )
!             SpcArr = SpcArr * 0.40_sp
!          CASE ( 'HNO3' )
!             SpcArr = SpcArr * 0.20_sp
!==============================================================================
       END SELECT

       ! Check for masking
       CALL HCOX_SCALE( HcoState, SpcArr, TRIM(Inst%SpcScalFldNme(N)), RC )
       IF ( RC /= HCO_SUCCESS ) RETURN

!==============================================================================
! This code is required for the vertical distribution of biomass burning emiss.
! We will keep it here for a future implementation. (mps, 4/24/17)
!
!       !--------------------------------------------------------------------
!       ! For grid boxes with emissions, distribute 65% to PBL and 35% to FT
!       !--------------------------------------------------------------------
!       DO J = 1, HcoState%Ny
!       DO I = 1, HcoState%Nx
!
!          IF ( SpcArr(I,J) > 0e+0_hp ) THEN
!
!             ! Initialize
!             PBL_MAX  = 1
!             F_OF_PBL = 0e+0_hp
!             F_OF_FT  = 0e+0_hp
!             DELTPRES = 0e+0_hp
!
!             ! Determine PBL height
!             DO L = HcoState%NZ, 1, -1
!                IF ( ExtState%FRAC_OF_PBL%Arr%Val(I,J,L) > 0.0_hp ) THEN
!                   PBL_MAX = L
!                   EXIT
!                ENDIF
!             ENDDO
!
!             ! Loop over the boundary layer
!             DO L = 1, PBL_MAX
!
!                ! Fraction of PBL that box (I,J,L) makes up [unitless]
!                F_OF_PBL = ExtState%FRAC_OF_PBL%Arr%Val(I,J,L)
!
!                ! Add only 65% biomass burning source to PBL
!                ! Distribute emissions thru the entire boundary layer
!                ! (mps from evf+tjb, 3/10/17)
!                SpcArr3D(I,J,L) = SpcArr(I,J) * PBL_FRAC * F_OF_PBL
!
!             ENDDO
!
!
!             ! Total thickness of the free troposphere [hPa]
!             ! (considered here to be 10 levels above the PBL)
!             TOTPRESFT = HcoState%Grid%PEDGE%Val(I,J,PBL_MAX+1) - &
!                         HcoState%Grid%PEDGE%Val(I,J,PBL_MAX+11)
!
!
!             ! Loop over the free troposphere
!             DO L = PBL_MAX+1, PBL_MAX+10
!
!                ! Thickness of level L [hPa]
!                DELTPRES = HcoState%Grid%PEDGE%Val(I,J,L) - &
!                           HcoState%Grid%PEDGE%Val(I,J,L+1)
!
!                ! Fraction of FT that box (I,J,L) makes up [unitless]
!                F_OF_FT = DELTPRES / TOTPRESFT
!
!                ! Add 35% of biomass burning source to free troposphere
!                ! Distribute emissions thru 10 model levels above the BL
!                ! (mps from evf+tjb, 3/10/17)
!                SpcArr3D(I,J,L) = SpcArr(I,J) * (1.0-PBL_FRAC) * F_OF_FT
!
!             ENDDO
!
!          ENDIF
!
!       ENDDO
!       ENDDO
!
!       ! Add flux to HEMCO emission array
!       ! Now 3D flux (mps, 3/10/17)
!       CALL HCO_EmisAdd( HcoState,   SpcArr3D, HcoIDs(N), &
!                         RC,        ExtNr=ExtNr )
!==============================================================================

       ! Add flux to HEMCO emission array
       CALL HCO_EmisAdd( HcoState, SpcArr, Inst%HcoIDs(N), RC, ExtNr=Inst%ExtNr )
       IF ( RC /= HCO_SUCCESS ) THEN
          MSG = 'HCO_EmisAdd error: ' // TRIM(HcoState%Spc(Inst%HcoIDs(N))%SpcName)
          CALL HCO_ERROR(HcoState%Config%Err,MSG, RC )
          RETURN
       ENDIF

    ENDDO !N

    ! Nullify pointers for safety's sake
    TmpPtr  => NULL()
    Inst    => NULL()

    ! Leave w/ success
    CALL HCO_LEAVE( HcoState%Config%Err,RC )

  END SUBROUTINE HCOX_GFED_Run
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCOX_GFED_Init
!
! !DESCRIPTION: Subroutine HcoX\_GFED\_Init initializes all
!  extension variables.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCOX_GFED_Init ( HcoState, ExtName, ExtState, RC )
!
! !USES:
!
    USE HCO_STATE_MOD,          ONLY : HCO_GetHcoID
    USE HCO_STATE_MOD,          ONLY : HCO_GetExtHcoID
    USE HCO_ExtList_Mod,        ONLY : GetExtNr, GetExtOpt
    USE HCO_ExtList_Mod,        ONLY : GetExtSpcVal
!
! !INPUT/OUTPUT PARAMETERS:
!
    CHARACTER(LEN=*), INTENT(IN   )  :: ExtName     ! Extension name
    TYPE(Ext_State),  POINTER        :: ExtState    ! Options object
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(HCO_State),  POINTER        :: HcoState    ! HEMCO state object
    INTEGER,          INTENT(INOUT)  :: RC          ! Return status
!
! !REVISION HISTORY:
!  07 Sep 2011 - P. Kasibhatla - Initial version, based on GFED2
!  15 Dec 2013 - C. Keller     - Now a HEMCO extension
!  08 Aug 2014 - R. Yantosca   - Now include hcox_gfed_include.H, which defines
!                                GFED_SPEC_NAME and GFED_EMFAC arrays
!  11 Nov 2014 - C. Keller     - Now get hydrophilic fractions via config file
!  22 Apr 2015 - R. Yantosca   - Now explicitly test for "POA scale factor"
!                                and "NAP scale factor" to avoid search errors
!  07 Jan 2016 - M. Sulprizio  - Change 'POA1' to 'POG1' to better reflect that
!                                SVOC emissions are added to the gas-phase
!                                species in carbon_mod.F
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    CHARACTER(LEN=255) :: MSG, ScalFile
    INTEGER            :: ExtNr, tmpNr, AS, IU_FILE, IOS
    INTEGER            :: nSpc, N, M, NDUM, NCHAR
    CHARACTER(LEN=31)  :: tmpName
    CHARACTER(LEN=31)  :: SpcName
    LOGICAL            :: FOUND, Matched
    REAL(sp)           :: ValSp
    TYPE(MyInst), POINTER :: Inst

    CHARACTER(LEN=255), POINTER :: GFED_SPEC_NAME (:) => NULL()
    CHARACTER(LEN=255), TARGET  :: GFED4_SPEC_NAME(N_SPEC)

    CHARACTER(LEN=31), ALLOCATABLE :: SpcNames(:)
    CHARACTER(LEN=61), ALLOCATABLE :: SpcScalFldNme(:)
    INTEGER,           ALLOCATABLE :: HcoIDs(:)
    REAL(sp),          ALLOCATABLE :: SpcScal(:)

   !=================================================================

    !=================================================================
    ! HCOX_GFED_Init begins here!
    !=================================================================

    ! Extension Nr.
    ExtNr = GetExtNr( HcoState%Config%ExtList, TRIM(ExtName) )
    IF ( ExtNr <= 0 ) RETURN

    ! Enter
    CALL HCO_ENTER( HcoState%Config%Err, 'HCOX_GFED_Init (hcox_gfed_mod.F90)', RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Create local instance for this simulation
    Inst => NULL()
    CALL InstCreate ( ExtNr, ExtState%GFED, Inst, RC )
    IF ( RC /= HCO_SUCCESS ) THEN
       CALL HCO_ERROR ( HcoState%Config%Err, 'Cannot create GFED instance', RC )
       RETURN
    ENDIF

    ! Check if this is GFED4
    CALL GetExtOpt( HcoState%Config, Inst%ExtNr, 'GFED4', &
                     OptValBool=Inst%IsGFED4, FOUND=FOUND, RC=RC )
    IF ( .NOT. FOUND ) THEN
       Inst%IsGFED4 = .FALSE.
    ENDIF

    ! Error checks
    IF ( .NOT. Inst%IsGFED4  ) THEN
       MSG = 'GFED is enabled but no GFED version is selected. ' // &
             'Please set GFED4 in HEMCO configuration file.'
       CALL HCO_ERROR(HcoState%Config%Err,MSG, RC )
       RETURN
    ENDIF

    ! ----------------------------------------------------------------------
    ! Get settings
    ! The speciation of carbon aerosols into hydrophilic and hydrophobic
    ! fractions can be specified in the configuration file, e.g.:
    ! 100     GFED           : on    NO/CO/OCPI/OCPO/BCPI/BCPO
    !    --> hydrophilic BC  :       0.2
    !    --> hydrophilic OC  :       0.5
    !
    ! Setting these values is optional and default values are applied if
    ! they are not specified. The values only take effect if the
    ! corresponding species (CO, BCPI/BCPO, OCPI/OCPO) are listed as species
    ! to be used.
    ! ----------------------------------------------------------------------

    ! Try to read hydrophilic fractions of BC. Defaults to 0.2.
    CALL GetExtOpt( HcoState%Config, Inst%ExtNr, 'hydrophilic BC', &
                     OptValSp=ValSp, FOUND=FOUND, RC=RC )
    IF ( RC /= HCO_SUCCESS ) RETURN
    IF ( .NOT. FOUND ) THEN
       Inst%BCPIfrac = 0.2
    ELSE
       Inst%BCPIfrac = ValSp
    ENDIF

    ! Try to read hydrophilic fractions of OC. Defaults to 0.5.
    CALL GetExtOpt( HcoState%Config, Inst%ExtNr, 'hydrophilic OC', &
                     OptValSp=ValSp, FOUND=FOUND, RC=RC )
    IF ( RC /= HCO_SUCCESS ) RETURN
    IF ( .NOT. FOUND ) THEN
       Inst%OCPIfrac = 0.5
    ELSE
       Inst%OCPIfrac = ValSp
    ENDIF

    ! Try to read POG1 fraction of SVOC. Defaults to 0.49.
    CALL GetExtOpt ( HcoState%Config, ExtNr, 'fraction POG1', &
                     OptValSp=ValSp, FOUND=FOUND, RC=RC )
    IF ( RC /= HCO_SUCCESS ) RETURN
    IF ( .NOT. FOUND ) THEN
       Inst%POG1frac = 0.49
    ELSE
       Inst%POG1frac = ValSp
    ENDIF

    CALL GetExtOpt( HcoState%Config, ExtNr, 'CO to SOAP', &
                     OptValSp=ValSp, FOUND=FOUND, RC=RC )
    IF ( RC /= HCO_SUCCESS ) RETURN
    IF ( .NOT. FOUND ) THEN
       Inst%SOAPfrac = 0.0
    ELSE
       Inst%SOAPfrac = ValSp
    ENDIF

    ! Error check: OCPIfrac, BCPIfrac, and POG1frac must be between 0 and 1
    IF ( Inst%OCPIfrac < 0.0_sp .OR. Inst%OCPIfrac > 1.0_sp .OR. &
         Inst%BCPIfrac < 0.0_sp .OR. Inst%BCPIfrac > 1.0_sp .OR. &
         Inst%SOAPfrac < 0.0_sp .OR. Inst%SOAPfrac > 1.0_sp .OR. &
         Inst%POG1frac < 0.0_sp .OR. Inst%POG1frac > 1.0_sp     ) THEN
       WRITE(MSG,*) 'fractions must be between 0-1: ', &
          Inst%OCPIfrac, Inst%BCPIfrac, Inst%POG1frac, Inst%SOAPfrac
       CALL HCO_ERROR(HcoState%Config%Err,MSG, RC )
       RETURN
    ENDIF

    ! Use daily scale factors?
    CALL GetExtOpt( HcoState%Config, Inst%ExtNr, 'GFED_daily', &
                     OptValBool=Inst%DoDay, FOUND=FOUND, RC=RC )
    IF ( RC /= HCO_SUCCESS ) RETURN
    IF ( .NOT. FOUND ) THEN
       Inst%DoDay = .FALSE.
    ENDIF

    ! Use 3-hourly scale factors?
    CALL GetExtOpt( HcoState%Config, ExtNr, 'GFED_3hourly', &
                     OptValBool=Inst%Do3Hr, FOUND=FOUND, RC=RC )
    IF ( RC /= HCO_SUCCESS ) RETURN
    IF ( .NOT. FOUND ) THEN
       Inst%Do3Hr = .FALSE.
    ENDIF

    !-----------------------------------------------------------------------
    ! Initialize GFED scale factors
    !-----------------------------------------------------------------------

    ! Allocate scale factors table
    ALLOCATE ( Inst%GFED4_EMFAC ( N_SPEC, N_EMFAC ), STAT=AS )
    IF ( AS/=0 ) THEN
       CALL HCO_ERROR( HcoState%Config%Err, 'Cannot allocate GFED_EMFAC', RC )
       RETURN
    ENDIF
    Inst%GFED4_EMFAC = 0.0_hp

    ALLOCATE( Inst%GFED_SAVA(HcoState%NX,HcoState%NY) )
    ALLOCATE( Inst%GFED_BORF(HcoState%NX,HcoState%NY) )
    ALLOCATE( Inst%GFED_TEMP(HcoState%NX,HcoState%NY) )
    ALLOCATE( Inst%GFED_DEFO(HcoState%NX,HcoState%NY) )
    ALLOCATE( Inst%GFED_PEAT(HcoState%NX,HcoState%NY) )
    ALLOCATE( Inst%GFED_AGRI(HcoState%NX,HcoState%NY) )
    ALLOCATE( Inst%DAYSCAL (HcoState%NX,HcoState%NY) )
    ALLOCATE( Inst%HRSCAL  (HcoState%NX,HcoState%NY) )

    ! Now get definitions for GFED_EMFAC and GFED_SPEC_NAME from an include
    ! file.  This avoids ASCII file reads in the ESMF environment.  To update
    ! the emission factors, one just needs to modify the include file.
    ! This can be done with the script HEMCO/Extensions/Preprocess/gfed.pl,
    ! (bmy, 8/14/14)
#include "hcox_gfed_include_gfed4.H"

    ! Set working pointers
    IF ( Inst%IsGFED4 ) THEN
       Inst%GFED_EMFAC => Inst%GFED4_EMFAC
       GFED_SPEC_NAME  => GFED4_SPEC_NAME
    ENDIF

    !-----------------------------------------------------------------------
    ! Match specified species with GFED species
    ! The species to be used are specified in the HEMCO configuration file.
    ! Match these species with the ones found in the scale factors table.
    !-----------------------------------------------------------------------

    ! Prompt to log file
    IF ( HcoState%amIRoot ) THEN
       MSG = 'Use GFED extension'
       CALL HCO_MSG(HcoState%Config%Err,MSG, SEP1='-' )
       WRITE(MSG,*) '   - Use GFED-4              : ', Inst%IsGFED4
       CALL HCO_MSG(HcoState%Config%Err,MSG )
       WRITE(MSG,*) '   - Use daily scale factors : ', Inst%DoDay
       CALL HCO_MSG(HcoState%Config%Err,MSG )
       WRITE(MSG,*) '   - Use hourly scale factors: ', Inst%Do3Hr
       CALL HCO_MSG(HcoState%Config%Err,MSG )
       WRITE(MSG,*) '   - Hydrophilic OC fraction : ', Inst%OCPIfrac
       CALL HCO_MSG(HcoState%Config%Err,MSG )
       WRITE(MSG,*) '   - Hydrophilic BC fraction : ', Inst%BCPIfrac
       CALL HCO_MSG(HcoState%Config%Err,MSG )
       WRITE(MSG,*) '   - POG1 fraction           : ', Inst%POG1frac
       CALL HCO_MSG(HcoState%Config%Err,MSG )
       WRITE(MSG,*) '   - SOAP fraction           : ', Inst%SOAPfrac
       CALL HCO_MSG(HcoState%Config%Err,MSG )
    ENDIF

    ! Get HEMCO species IDs of all species specified in configuration file
    CALL HCO_GetExtHcoID( HcoState, Inst%ExtNr, HcoIDs, SpcNames, Inst%nSpc, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN
    IF ( Inst%nSpc == 0 ) THEN
       MSG = 'No GFED species specified'
       CALL HCO_ERROR(HcoState%Config%Err,MSG, RC )
       RETURN
    ENDIF
    ALLOCATE(Inst%HcoIDs(Inst%nSpc),Inst%SpcNames(Inst%nSpc))
    Inst%HcoIDs   = HcoIDs
    Inst%SpcNames = SpcNames
    DEALLOCATE(HcoIDs,SpcNames)

    ! Get species scale factors
    CALL GetExtSpcVal( HcoState%Config, Inst%ExtNr, Inst%nSpc, &
                       Inst%SpcNames, 'Scaling', 1.0_sp, SpcScal, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Get species mask fields
    CALL GetExtSpcVal( HcoState%Config, Inst%ExtNr, Inst%nSpc, &
                       Inst%SpcNames, 'ScaleField', HCOX_NOSCALE, SpcScalFldNme, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Pass to instance
    nSpc = Inst%nSpc
    ALLOCATE(Inst%SpcScal(nSpc),Inst%SpcScalFldNme(nSpc))
    Inst%SpcScal       = SpcScal
    Inst%SpcScalFldNme = SpcScalFldNme
    DEALLOCATE(SpcScal,SpcScalFldNme)

    ! Error trap: in previous versions, CO, POA and NAP scale factor were given as
    ! 'CO scale factor', etc. Make sure those attributes do not exist any more!
    CALL GetExtOpt( HcoState%Config, Inst%ExtNr, 'CO scale factor', &
                     OptValSp=ValSp, FOUND=FOUND, RC=RC )
    IF ( .NOT. FOUND ) THEN
       CALL GetExtOpt( HcoState%Config, Inst%ExtNr, 'POA scale factor', &
                        OptValSp=ValSp, FOUND=FOUND, RC=RC )
    ENDIF
    IF ( .NOT. FOUND ) THEN
       CALL GetExtOpt( HcoState%Config, Inst%ExtNr, 'NAP scale factor', &
                        OptValSp=ValSp, FOUND=FOUND, RC=RC )
    ENDIF
    IF ( FOUND ) THEN
       MSG = 'Found old definition of CO, POA and/or NAP scale factor! '  // &
             'This version of HEMCO expects species scale factors to be ' // &
             'set as `Scaling_XX` instead of `XX scale factor`. '         // &
             'Please update the GFED settings section accordingly.'
       CALL HCO_ERROR(HcoState%Config%Err,MSG, RC )
       RETURN
    ENDIF

    ! GFEDIDS are the matching indeces of the HEMCO species in GFED_EMFAC.
    ALLOCATE ( Inst%GfedIDs(Inst%nSpc), STAT=AS )
    IF ( AS/=0 ) THEN
       CALL HCO_ERROR( HcoState%Config%Err, 'Cannot allocate GfedIDs', RC )
       RETURN
    ENDIF
    Inst%GfedIDs = -1

    ! Find matching GFED index for each specified species
    DO N = 1, Inst%nSpc
       IF ( Inst%HcoIDs(N) < 0 ) CYCLE

       ! SpcName is the GFED species name to be searched. Adjust
       ! if necessary.
       SpcName = Inst%SpcNames(N)
       NCHAR   = LEN(SpcName)
       IF ( NCHAR > 3 ) THEN
          IF ( SpcName(1:3) == 'CO2' ) THEN
             SpcName = 'CO2'
          ELSEIF ( SpcName(1:3) == 'CH4' ) THEN
             SpcName = 'CH4'
          ELSEIF ( SpcName(1:3) == 'CO_' ) THEN
             SpcName = 'CO'
          ELSEIF ( SpcName(1:2) == 'BC' ) THEN
             SpcName = 'BC'
          ELSEIF ( SpcName(1:2) == 'OC' ) THEN
             SpcName = 'OC'
          ENDIF
       ENDIF
       IF ( TRIM(SpcName) == 'POG1' ) SpcName = 'OC'
       IF ( TRIM(SpcName) == 'POG2' ) SpcName = 'OC'
       IF ( TRIM(SpcName) == 'NAP'  ) SpcName = 'CO'
!==============================================================================
! This code is required for partitioning NOx emissions directly to PAN and HNO3.
! We will keep it here as an option for users focusing on North American fires.
! (mps, 5/12/17)
!       IF ( TRIM(SpcName) == 'PAN'  ) SpcName = 'NO'
!       IF ( TRIM(SpcName) == 'HNO3' ) SpcName = 'NO'
!==============================================================================

       ! adjust SOAP scale factor by CO scale factor (SOAP co-emitted with CO)
       IF ( TRIM(SpcName) == 'CO' ) THEN
         Inst%SOAPfrac = Inst%SOAPfrac * Inst%SpcScal(N)
       END IF

       ! Search for matching GFED species by name
       Matched = .FALSE.
       DO M = 1, N_SPEC

          IF ( TRIM(SpcName) == TRIM(GFED_SPEC_NAME(M)) ) THEN
             Inst%GfedIDs(N) = M
             Matched    = .TRUE.

             ! Verbose
             IF ( HcoState%amIRoot ) THEN
                MSG = '   - Emit GFED species ' // TRIM(GFED_SPEC_NAME(M)) // &
                      '     as model species ' // TRIM(Inst%SpcNames(N))
                CALL HCO_MSG(HcoState%Config%Err,MSG )
                WRITE(MSG,*) '     --> Will use scale factor: ', Inst%SpcScal(N)
                CALL HCO_MSG(HcoState%Config%Err,MSG )
                WRITE(MSG,*) '     --> Will use scale field : ', TRIM(Inst%SpcScalFldNme(N))
                CALL HCO_MSG(HcoState%Config%Err,MSG )
             ENDIF
             EXIT ! go to next species
          ENDIF
       ENDDO
       IF ( .NOT. Matched ) THEN
          MSG = 'Species '// TRIM(SpcName) //' not found in GFED'
          CALL HCO_ERROR(HcoState%Config%Err,MSG, RC )
          RETURN
       ENDIF
    ENDDO !N

    !=======================================================================
    ! Activate this module and the fields of ExtState that it uses
    !=======================================================================

!==============================================================================
! This code is required for the vertical distribution of biomass burning emiss.
! We will keep it here for a future implementation. (mps, 4/24/17)
!    ! Activate met fields required by this extension
!    ExtState%FRAC_OF_PBL%DoUse = .TRUE.
!==============================================================================

    ! Enable module
    !ExtState%GFED = .TRUE.

    ! Cleanup
    GFED_SPEC_NAME => NULL()
    Inst           => NULL()

    ! Return w/ success
    CALL HCO_LEAVE( HcoState%Config%Err,RC )

  END SUBROUTINE HCOX_GFED_Init
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCOX_GFED_Final
!
! !DESCRIPTION: Subroutine HcoX\_GFED\_Final deallocates
!  all module arrays.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCOX_GFED_Final ( ExtState )
!
! !INPUT PARAMETERS:
!
    TYPE(Ext_State),  POINTER       :: ExtState   ! Module options
!
! !REVISION HISTORY:
!  07 Sep 2011 - P. Kasibhatla - Initial version, based on GFED2
!  15 Dec 2013 - C. Keller     - Now a HEMCO extension
!EOP
!------------------------------------------------------------------------------
!BOC
!
    !=================================================================
    ! HCOX_GFED_Final begins here!
    !=================================================================

    CALL InstRemove ( ExtState%GFED )

  END SUBROUTINE HCOX_GFED_Final
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
          ! Free pointers
          Inst%GFED_EMFAC => NULL()

          DEALLOCATE( Inst%GFED_SAVA)
          DEALLOCATE( Inst%GFED_BORF)
          DEALLOCATE( Inst%GFED_TEMP)
          DEALLOCATE( Inst%GFED_DEFO)
          DEALLOCATE( Inst%GFED_PEAT)
          DEALLOCATE( Inst%GFED_AGRI)
          DEALLOCATE( Inst%DAYSCAL )
          DEALLOCATE( Inst%HRSCAL  )

          ! Cleanup module arrays
          IF ( ASSOCIATED( Inst%GFED4_EMFAC  ) ) DEALLOCATE( Inst%GFED4_EMFAC  )
          IF ( ASSOCIATED( Inst%GfedIDs      ) ) DEALLOCATE( Inst%GfedIds      )
          IF ( ASSOCIATED( Inst%HcoIDs       ) ) DEALLOCATE( Inst%HcoIDs       )
          IF ( ASSOCIATED( Inst%SpcNames     ) ) DEALLOCATE( Inst%SpcNames     )
          IF ( ASSOCIATED( Inst%SpcScal      ) ) DEALLOCATE( Inst%SpcScal      )
          IF ( ASSOCIATED( Inst%SpcScalFldNme) ) DEALLOCATE( Inst%SpcScalFldNme)

          PrevInst%NextInst => Inst%NextInst
       ELSE
          AllInst => Inst%NextInst
       ENDIF
       DEALLOCATE(Inst)
       Inst => NULL()
    ENDIF

   END SUBROUTINE InstRemove
!EOC
END MODULE HCOX_GFED_MOD
