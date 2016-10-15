!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: hcox_gfed_mod.F90
!
! !DESCRIPTION: Module HCOX\_GFED\_MOD contains routines to calculate
! GFED biomass burning emissions in HEMCO. This can be GFED-3 or GFED-4,
! depending on the input data selected in the HEMCO configuration file. 
!
! !NOTES:
!
! !REFERENCES:
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
!111     GFED              : on    NO/CO/ALK4/ACET/MEK/ALD2/PRPE/C3H8/CH2O/C2H6/SO2/NH3/BCPI/BCPO/OCPI/OCPO
!    --> GFED3             :       false
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
!  14 Oct 2016 - C. Keller   - Now use HCO_EvalFld instead of HCO_GetPtr.
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
  INTEGER,           PARAMETER :: N_SPEC  = 26
!
! !PRIVATE TYPES:
!
  !=================================================================
  ! HEMCO VARIABLES 
  !
  ! ExtNr   : Extension number 
  ! DoDay   : TRUE if dialy scale factors are used 
  ! Do3Hr   : TRUE if 3-hourly scale factors are used 
  !=================================================================
  INTEGER                       :: ExtNr
  LOGICAL                       :: DoDay
  LOGICAL                       :: Do3Hr
  LOGICAL                       :: IsGFED3
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
  INTEGER                        :: nSpc
  CHARACTER(LEN=31), ALLOCATABLE :: SpcNames(:)
  CHARACTER(LEN=61), ALLOCATABLE :: SpcScalFldNme(:)
  INTEGER,           ALLOCATABLE :: HcoIDs(:)
  INTEGER,           ALLOCATABLE :: GfedIDs(:)
  REAL(sp),          ALLOCATABLE :: SpcScal(:)

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
  REAL(hp), ALLOCATABLE, TARGET  :: GFED3_EMFAC(:,:)
  REAL(hp), ALLOCATABLE, TARGET  :: GFED4_EMFAC(:,:)
  REAL(hp),              POINTER :: GFED_EMFAC (:,:) => NULL()
  REAL(sp)                       :: OCPIfrac 
  REAL(sp)                       :: BCPIfrac
  REAL(sp)                       :: POG1frac

  !=================================================================
  ! DATA ARRAY POINTERS 
  !
  ! These are the pointers to the 6 input data specified in the 
  ! the configuration file
  !=================================================================
!  REAL(sp), POINTER   :: GFED_WDL(:,:) => NULL()
!  REAL(sp), POINTER   :: GFED_SAV(:,:) => NULL()
!  REAL(sp), POINTER   :: GFED_PET(:,:) => NULL()
!  REAL(sp), POINTER   :: GFED_FOR(:,:) => NULL()
!  REAL(sp), POINTER   :: GFED_AGW(:,:) => NULL()
!  REAL(sp), POINTER   :: GFED_DEF(:,:) => NULL()
!  REAL(sp), POINTER   :: HUMTROP (:,:) => NULL()
!  REAL(sp), POINTER   :: DAYSCAL (:,:) => NULL()
!  REAL(sp), POINTER   :: HRSCAL  (:,:) => NULL()
 
  REAL(hp), ALLOCATABLE, TARGET   :: GFED_WDL(:,:)
  REAL(hp), ALLOCATABLE, TARGET   :: GFED_SAV(:,:)
  REAL(hp), ALLOCATABLE, TARGET   :: GFED_PET(:,:)
  REAL(hp), ALLOCATABLE, TARGET   :: GFED_FOR(:,:)
  REAL(hp), ALLOCATABLE, TARGET   :: GFED_AGW(:,:)
  REAL(hp), ALLOCATABLE, TARGET   :: GFED_DEF(:,:)
  REAL(hp), ALLOCATABLE, TARGET   :: HUMTROP (:,:)
  REAL(hp), ALLOCATABLE, TARGET   :: DAYSCAL (:,:)
  REAL(hp), ALLOCATABLE, TARGET   :: HRSCAL  (:,:)




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
  SUBROUTINE HCOX_GFED_Run( am_I_Root, ExtState, HcoState, RC )
!
! !USES:
!
    USE HCO_Calc_Mod,     ONLY : HCO_EvalFld
    USE HCO_EmisList_Mod, ONLY : HCO_GetPtr
    USE HCO_FluxArr_MOD,  ONLY : HCO_EmisAdd
!
! !INPUT PARAMETERS:
!
    LOGICAL,         INTENT(IN   )  :: am_I_Root  ! root CPU?
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
   
    !=================================================================
    ! HCOX_GFED_Run begins here!
    !=================================================================

    ! Return if extension disabled 
    IF ( .NOT. ExtState%GFED ) RETURN

    ! Enter 
    CALL HCO_ENTER( HcoState%Config%Err, 'HCOX_GFED_Run (hcox_gfed_mod.F90)', RC ) 
    IF ( RC /= HCO_SUCCESS ) RETURN

    !-----------------------------------------------------------------
    ! Get pointers to data arrays 
    !-----------------------------------------------------------------
    !IF ( FIRST ) THEN

       ! Get pointers to GFED3 data 
       IF ( IsGFED3 ) THEN
          CALL HCO_EvalFld ( am_I_Root, HcoState, 'GFED_WDL', GFED_WDL, RC )
          IF ( RC /= HCO_SUCCESS ) RETURN
          CALL HCO_EvalFld ( am_I_Root, HcoState, 'GFED_SAV', GFED_SAV, RC )
          IF ( RC /= HCO_SUCCESS ) RETURN
          CALL HCO_EvalFld ( am_I_Root, HcoState, 'GFED_PET', GFED_PET, RC )
          IF ( RC /= HCO_SUCCESS ) RETURN
          CALL HCO_EvalFld ( am_I_Root, HcoState, 'GFED_FOR', GFED_FOR, RC )
          IF ( RC /= HCO_SUCCESS ) RETURN
          CALL HCO_EvalFld ( am_I_Root, HcoState, 'GFED_AGW', GFED_AGW, RC )
          IF ( RC /= HCO_SUCCESS ) RETURN
          CALL HCO_EvalFld ( am_I_Root, HcoState, 'GFED_DEF', GFED_DEF, RC )
          IF ( RC /= HCO_SUCCESS ) RETURN
          CALL HCO_EvalFld ( am_I_Root, HcoState, 'GFED_HUMTROP', HUMTROP, RC )
          IF ( RC /= HCO_SUCCESS ) RETURN

          ! Make sure HUMTROP does not exceed one
          WHERE ( HUMTROP > 1.0_sp ) 
             HUMTROP = 1.0_sp
          END WHERE

       ! Get pointers to GFED4 data
       ELSEIF ( IsGFED4 ) THEN
          CALL HCO_EvalFld ( am_I_Root, HcoState, 'GFED_TEMP', GFED_WDL, RC )
          IF ( RC /= HCO_SUCCESS ) RETURN
          CALL HCO_EvalFld ( am_I_Root, HcoState, 'GFED_SAVA', GFED_SAV, RC )
          IF ( RC /= HCO_SUCCESS ) RETURN
          CALL HCO_EvalFld ( am_I_Root, HcoState, 'GFED_PEAT', GFED_PET, RC )
          IF ( RC /= HCO_SUCCESS ) RETURN
          CALL HCO_EvalFld ( am_I_Root, HcoState, 'GFED_BORF', GFED_FOR, RC )
          IF ( RC /= HCO_SUCCESS ) RETURN
          CALL HCO_EvalFld ( am_I_Root, HcoState, 'GFED_AGRI', GFED_AGW, RC )
          IF ( RC /= HCO_SUCCESS ) RETURN
          CALL HCO_EvalFld ( am_I_Root, HcoState, 'GFED_DEFO', GFED_DEF, RC )
          IF ( RC /= HCO_SUCCESS ) RETURN
       ENDIF

       ! Also point to scale factors if needed
       IF ( DoDay ) THEN
          CALL HCO_EvalFld ( am_I_Root, HcoState, 'GFED_FRAC_DAY', &
                             DAYSCAL,   RC )
          IF ( RC /= HCO_SUCCESS ) RETURN
       ENDIF
       IF ( Do3Hr ) THEN
          CALL HCO_EvalFld ( am_I_Root, HcoState, 'GFED_FRAC_3HOUR', &
                             HRSCAL,    RC )
          IF ( RC /= HCO_SUCCESS ) RETURN
       ENDIF

       FIRST = .FALSE.
    !ENDIF

    !-----------------------------------------------------------------
    ! Calculate emissions for defined species
    !-----------------------------------------------------------------
    DO N = 1, nSpc

       ! Continue if species not defined
       IF ( HcoIDs(N)  < 0 ) CYCLE
       IF ( GfedIDs(N) < 0 ) CYCLE

       ! SpcArr are the total biomass burning emissions for this
       ! species. TypArr are the emissions from a given source type. 
       SpcArr = 0.0_hp

       ! Calculate emissions for all source types
       DO M = 1, N_EMFAC

          ! Point to the emission factor array for each source type
          SELECT CASE ( M ) 
             CASE( 1 )
                TMPPTR => GFED_AGW
             CASE( 2 )
                TMPPTR => GFED_DEF
             CASE( 3 )
                TMPPTR => GFED_FOR
             CASE( 4 )
                TMPPTR => GFED_PET
             CASE( 5 )
                TMPPTR => GFED_SAV
             CASE( 6 )
                TMPPTR => GFED_WDL
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
          IF ( M == 2 .AND. IsGFED3 ) THEN
             TypArr = TmpPtr *         HUMTROP  * GFED_EMFAC(GfedIDs(N),M) &
                    + TmpPtr * (1.0_sp-HUMTROP) * GFED_EMFAC(GfedIDs(N),6)
          ELSE
             TypArr = TmpPtr * GFED_EMFAC(GfedIDs(N),M)
          ENDIF

          ! Eventually add daily / 3-hourly scale factors. These scale
          ! factors are unitless.
          IF ( DoDay ) THEN
             !IF ( ASSOCIATED(DAYSCAL) ) THEN
                TypArr = TypArr * DAYSCAL
             !ENDIF
          ENDIF
          IF ( Do3Hr ) THEN
             !IF ( ASSOCIATED(HRSCAL) ) THEN
                TypArr = TypArr * HRSCAL
             !ENDIF
          ENDIF

          ! Add to output array
          SpcArr = SpcArr + TypArr

          ! Nullify pointer
          TmpPtr  => NULL()

       ENDDO !M

       ! Apply species specific scale factors
       SpcArr = SpcArr * SpcScal(N)

       SELECT CASE ( SpcNames(N) ) 
          CASE ( 'OCPI' )
             SpcArr = SpcArr * OCPIfrac
          CASE ( 'OCPO' )
             SpcArr = SpcArr * (1.0_sp - OCPIfrac)
          CASE ( 'BCPI' )
             SpcArr = SpcArr * BCPIfrac
          CASE ( 'BCPO' )
             SpcArr = SpcArr * (1.0_sp - BCPIfrac)
          CASE ( 'POG1' )
             SpcArr = SpcArr * POG1frac
          CASE ( 'POG2' )
             SpcArr = SpcArr * (1.0_sp - POG1frac)
       END SELECT

       ! Check for masking
       CALL HCOX_SCALE( am_I_Root, HcoState, SpcArr, TRIM(SpcScalFldNme(N)), RC ) 
       IF ( RC /= HCO_SUCCESS ) RETURN

       ! Add flux to HEMCO emission array
       CALL HCO_EmisAdd( am_I_Root, HcoState, SpcArr, HcoIDs(N), RC, ExtNr=ExtNr ) 
       IF ( RC /= HCO_SUCCESS ) THEN
          MSG = 'HCO_EmisAdd error: ' // TRIM(HcoState%Spc(HcoIDs(N))%SpcName)
          CALL HCO_ERROR(HcoState%Config%Err,MSG, RC )
          RETURN 
       ENDIF

    ENDDO !N

    ! Nullify pointers for safety's sake
    TmpPtr  => NULL()

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
  SUBROUTINE HCOX_GFED_Init ( am_I_Root, HcoState, ExtName, &
                              ExtState,  RC                  ) 
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
    LOGICAL,          INTENT(IN   )  :: am_I_Root   ! root CPU?
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
    INTEGER            :: tmpNr, AS, IU_FILE, IOS
    INTEGER            :: N, M, NDUM, NCHAR
    CHARACTER(LEN=31)  :: tmpName
    CHARACTER(LEN=31)  :: SpcName
    LOGICAL            :: FOUND, Matched
    REAL(sp)           :: ValSp

    CHARACTER(LEN=255), POINTER :: GFED_SPEC_NAME (:) => NULL()
    CHARACTER(LEN=255), TARGET  :: GFED3_SPEC_NAME(N_SPEC)
    CHARACTER(LEN=255), TARGET  :: GFED4_SPEC_NAME(N_SPEC)

    !=================================================================
    ! HCOX_GFED_Init begins here!
    !=================================================================

    ! Extension Nr.
    ExtNr = GetExtNr( HcoState%Config%ExtList, TRIM(ExtName) )
    IF ( ExtNr <= 0 ) RETURN

    ! Enter 
    CALL HCO_ENTER( HcoState%Config%Err, 'HCOX_GFED_Init (hcox_gfed_mod.F90)', RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Check if this is GFED3 or GFED4
    CALL GetExtOpt( HcoState%Config, ExtNr, 'GFED3', &
                     OptValBool=IsGFED3, FOUND=FOUND, RC=RC )
    IF ( .NOT. FOUND ) THEN
       IsGFED3 = .FALSE.
    ENDIF
    CALL GetExtOpt( HcoState%Config, ExtNr, 'GFED4', &
                     OptValBool=IsGFED4, FOUND=FOUND, RC=RC )
    IF ( .NOT. FOUND ) THEN
       IsGFED4 = .FALSE.
    ENDIF

    ! Error checks
    IF ( .NOT. IsGFED4 .AND. .NOT. IsGFED3 ) THEN
       MSG = 'GFED is enabled but no GFED version is selected. ' // &
             'Please set GFED3 or GFED4 in HEMCO configuration file.'
       CALL HCO_ERROR(HcoState%Config%Err,MSG, RC )
       RETURN
    ENDIF
    IF ( IsGFED4 .AND. IsGFED3 ) THEN
       MSG = 'Cannot use GFED3 and GFED4 together! Please select ' // &
             'only one model version in the HEMCO configuration file.'
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
    CALL GetExtOpt( HcoState%Config, ExtNr, 'hydrophilic BC', &
                     OptValSp=ValSp, FOUND=FOUND, RC=RC )
    IF ( RC /= HCO_SUCCESS ) RETURN
    IF ( .NOT. FOUND ) THEN
       BCPIfrac = 0.2
    ELSE
       BCPIfrac = ValSp
    ENDIF

    ! Try to read hydrophilic fractions of OC. Defaults to 0.5.
    CALL GetExtOpt( HcoState%Config, ExtNr, 'hydrophilic OC', &
                     OptValSp=ValSp, FOUND=FOUND, RC=RC )
    IF ( RC /= HCO_SUCCESS ) RETURN
    IF ( .NOT. FOUND ) THEN
       OCPIfrac = 0.5
    ELSE
       OCPIfrac = ValSp
    ENDIF

    ! Try to read POG1 fraction of SVOC. Defaults to 0.49.
    CALL GetExtOpt ( HcoState%Config, ExtNr, 'fraction POG1', &
                     OptValSp=ValSp, FOUND=FOUND, RC=RC )
    IF ( RC /= HCO_SUCCESS ) RETURN
    IF ( .NOT. FOUND ) THEN
       POG1frac = 0.49
    ELSE
       POG1frac = ValSp
    ENDIF

    ! Error check: OCPIfrac, BCPIfrac, and POG1frac must be between 0 and 1
    IF ( OCPIfrac < 0.0_sp .OR. OCPIfrac > 1.0_sp .OR. &
         BCPIfrac < 0.0_sp .OR. BCPIfrac > 1.0_sp .OR. &
         POG1frac < 0.0_sp .OR. POG1frac > 1.0_sp     ) THEN
       WRITE(MSG,*) 'fractions must be between 0-1: ', &
          OCPIfrac, BCPIfrac, POG1frac
       CALL HCO_ERROR(HcoState%Config%Err,MSG, RC )
       RETURN
    ENDIF

    ! Use daily scale factors?
    CALL GetExtOpt( HcoState%Config, ExtNr, 'GFED_daily', &
                     OptValBool=DoDay, FOUND=FOUND, RC=RC )
    IF ( RC /= HCO_SUCCESS ) RETURN
    IF ( .NOT. FOUND ) THEN
       DoDay = .FALSE.
    ENDIF 

    ! Use 3-hourly scale factors?
    CALL GetExtOpt( HcoState%Config, ExtNr, 'GFED_3hourly', &
                     OptValBool=Do3Hr, FOUND=FOUND, RC=RC )
    IF ( RC /= HCO_SUCCESS ) RETURN
    IF ( .NOT. FOUND ) THEN
       Do3Hr = .FALSE.
    ENDIF 

    !----------------------------------------------------------------------- 
    ! Initialize GFED scale factors
    !----------------------------------------------------------------------- 

    ! Allocate scale factors table
    ALLOCATE ( GFED3_EMFAC ( N_SPEC, N_EMFAC ),        &
               GFED4_EMFAC ( N_SPEC, N_EMFAC ), STAT=AS )
    IF ( AS/=0 ) THEN
       CALL HCO_ERROR( HcoState%Config%Err, 'Cannot allocate GFED_EMFAC', RC )
       RETURN
    ENDIF
    GFED4_EMFAC = 0.0_hp
    GFED3_EMFAC = 0.0_hp

    ALLOCATE( GFED_WDL(HcoState%NX,HcoState%NY) )
    ALLOCATE( GFED_SAV(HcoState%NX,HcoState%NY) )
    ALLOCATE( GFED_PET(HcoState%NX,HcoState%NY) )
    ALLOCATE( GFED_FOR(HcoState%NX,HcoState%NY) )
    ALLOCATE( GFED_AGW(HcoState%NX,HcoState%NY) )
    ALLOCATE( GFED_DEF(HcoState%NX,HcoState%NY) )
    ALLOCATE( HUMTROP (HcoState%NX,HcoState%NY) )
    ALLOCATE( DAYSCAL (HcoState%NX,HcoState%NY) )
    ALLOCATE( HRSCAL  (HcoState%NX,HcoState%NY) )

    ! Now get definitions for GFED_EMFAC and GFED_SPEC_NAME from an include 
    ! file.  This avoids ASCII file reads in the ESMF environment.  To update
    ! the emission factors, one just needs to modify the include file.
    ! This can be done with the script HEMCO/Extensions/Preprocess/gfed.pl,
    ! (bmy, 8/14/14)
#include "hcox_gfed_include_gfed3.H"
#include "hcox_gfed_include_gfed4.H"

    ! Set working pointers
    IF ( IsGFED3 ) THEN
       GFED_EMFAC     => GFED3_EMFAC
       GFED_SPEC_NAME => GFED3_SPEC_NAME
    ELSEIF ( IsGFED4 ) THEN
       GFED_EMFAC     => GFED4_EMFAC
       GFED_SPEC_NAME => GFED4_SPEC_NAME
    ENDIF

    !----------------------------------------------------------------------- 
    ! Match specified species with GFED species
    ! The species to be used are specified in the HEMCO configuration file.
    ! Match these species with the ones found in the scale factors table.
    !----------------------------------------------------------------------- 

    ! Prompt to log file
    IF ( am_I_Root ) THEN
       MSG = 'Use GFED extension'
       CALL HCO_MSG(HcoState%Config%Err,MSG, SEP1='-' )
       WRITE(MSG,*) '   - Use GFED-3              : ', IsGFED3 
       CALL HCO_MSG(HcoState%Config%Err,MSG )
       WRITE(MSG,*) '   - Use GFED-4              : ', IsGFED4 
       CALL HCO_MSG(HcoState%Config%Err,MSG )
       WRITE(MSG,*) '   - Use daily scale factors : ', DoDay 
       CALL HCO_MSG(HcoState%Config%Err,MSG )
       WRITE(MSG,*) '   - Use hourly scale factors: ', Do3Hr
       CALL HCO_MSG(HcoState%Config%Err,MSG )
       WRITE(MSG,*) '   - Hydrophilic OC fraction : ', OCPIfrac
       CALL HCO_MSG(HcoState%Config%Err,MSG )
       WRITE(MSG,*) '   - Hydrophilic BC fraction : ', BCPIfrac
       CALL HCO_MSG(HcoState%Config%Err,MSG )
       WRITE(MSG,*) '   - POG1 fraction           : ', POG1frac
       CALL HCO_MSG(HcoState%Config%Err,MSG )
    ENDIF

    ! Get HEMCO species IDs of all species specified in configuration file
    CALL HCO_GetExtHcoID( HcoState, ExtNr, HcoIDs, SpcNames, nSpc, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN
    IF ( nSpc == 0 ) THEN
       MSG = 'No GFED species specified'
       CALL HCO_ERROR(HcoState%Config%Err,MSG, RC ) 
       RETURN
    ENDIF

    ! Get species scale factors
    CALL GetExtSpcVal( HcoState%Config, ExtNr, nSpc, &
                       SpcNames, 'Scaling', 1.0_sp, SpcScal, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Get species mask fields
    CALL GetExtSpcVal( HcoState%Config, ExtNr, nSpc, &
                       SpcNames, 'ScaleField', HCOX_NOSCALE, SpcScalFldNme, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Error trap: in previous versions, CO, POA and NAP scale factor were given as
    ! 'CO scale factor', etc. Make sure those attributes do not exist any more!
    CALL GetExtOpt( HcoState%Config, ExtNr, 'CO scale factor', &
                     OptValSp=ValSp, FOUND=FOUND, RC=RC )
    IF ( .NOT. FOUND ) THEN
       CALL GetExtOpt( HcoState%Config, ExtNr, 'POA scale factor', &
                        OptValSp=ValSp, FOUND=FOUND, RC=RC )
    ENDIF
    IF ( .NOT. FOUND ) THEN
       CALL GetExtOpt( HcoState%Config, ExtNr, 'NAP scale factor', &
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
    ALLOCATE ( GfedIDs(nSpc), STAT=AS )
    IF ( AS/=0 ) THEN
       CALL HCO_ERROR( HcoState%Config%Err, 'Cannot allocate GfedIDs', RC )
       RETURN
    ENDIF
    GfedIDs = -1

    ! Find matching GFED index for each specified species
    DO N = 1, nSpc
       IF ( HcoIDs(N) < 0 ) CYCLE

       ! SpcName is the GFED species name to be searched. Adjust
       ! if necessary.
       SpcName = SpcNames(N)
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

       ! Search for matching GFED species by name
       Matched = .FALSE.
       DO M = 1, N_SPEC

          IF ( TRIM(SpcName) == TRIM(GFED_SPEC_NAME(M)) ) THEN
             GfedIDs(N) = M
             Matched    = .TRUE.

             ! Verbose
             IF ( am_I_Root ) THEN
                MSG = '   - Emit GFED species ' // TRIM(GFED_SPEC_NAME(M)) // &
                      '     as model species ' // TRIM(SpcNames(N))
                CALL HCO_MSG(HcoState%Config%Err,MSG )
                WRITE(MSG,*) '     --> Will use scale factor: ', SpcScal(N)
                CALL HCO_MSG(HcoState%Config%Err,MSG )
                WRITE(MSG,*) '     --> Will use scale field : ', TRIM(SpcScalFldNme(N))
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

    ! Enable module
    ExtState%GFED = .TRUE.

    ! Cleanup
    GFED_SPEC_NAME => NULL()

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
  SUBROUTINE HCOX_GFED_Final
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

    ! Free pointers
!    GFED_WDL   => NULL()
!    GFED_SAV   => NULL()
!    GFED_PET   => NULL()
!    GFED_FOR   => NULL()
!    GFED_AGW   => NULL()
!    GFED_DEF   => NULL()
!    HUMTROP    => NULL()
!    DAYSCAL    => NULL()
!    HRSCAL     => NULL()
    GFED_EMFAC => NULL()
  
    DEALLOCATE( GFED_WDL)
    DEALLOCATE( GFED_SAV)
    DEALLOCATE( GFED_PET)
    DEALLOCATE( GFED_FOR)
    DEALLOCATE( GFED_AGW)
    DEALLOCATE( GFED_DEF)
    DEALLOCATE( HUMTROP )
    DEALLOCATE( DAYSCAL )
    DEALLOCATE( HRSCAL  )

    ! Cleanup module arrays
    IF ( ALLOCATED( GFED3_EMFAC  ) ) DEALLOCATE( GFED3_EMFAC  )
    IF ( ALLOCATED( GFED4_EMFAC  ) ) DEALLOCATE( GFED4_EMFAC  )
    IF ( ALLOCATED( GfedIDs      ) ) DEALLOCATE( GfedIds      )
    IF ( ALLOCATED( HcoIDs       ) ) DEALLOCATE( HcoIDs       )
    IF ( ALLOCATED( SpcNames     ) ) DEALLOCATE( SpcNames     )
    IF ( ALLOCATED( SpcScal      ) ) DEALLOCATE( SpcScal      )
    IF ( ALLOCATED( SpcScalFldNme) ) DEALLOCATE( SpcScalFldNme)

  END SUBROUTINE HCOX_GFED_Final
!EOC
END MODULE HCOX_GFED_MOD
