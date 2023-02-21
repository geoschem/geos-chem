!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: species_database_mod.F90
!
! !DESCRIPTION: Module SPECIES\_DATABASE\_MOD contains routines to set up
!  a database object containing physical properties for each GEOS-Chem
!  species.  This allows us to consolidate all species properties into a
!  single data structure, for convenience.
!\\
!\\
! !INTERFACE:
!
MODULE Species_Database_Mod
!
! !USES:
!
  USE Precision_Mod

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC  :: Init_Species_Database
  PUBLIC  :: Cleanup_Species_Database

#if defined ( EXTERNAL_GRID ) || defined( EXTERNAL_FORCING )
  !-----------------------------------------------------------------
  !         %%%%%%% GEOS-Chem HP (with ESMF & MPI) %%%%%%%
  !
  ! Cleanup routines for restoring the internal state of this
  ! module are exposed, so the DB can be reset from an external
  ! interface to perform multiple initializations of
  ! chemistry states. (hplin, 6/4/18)
  !-----------------------------------------------------------------
  PUBLIC  :: Cleanup_Work_Arrays
#endif
!
! !PRIVATE MEMBER FUNCTIONS:
!
  PRIVATE :: TranUc
!
! !REVISION HISTORY:
!  28 Aug 2015 - R. Yantosca - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC

  ! Work array to hold the list of species names, which combines the advected
  ! species from geoschem_config.yml with the KPP species names (and removes
  ! duplicates)
  CHARACTER(LEN=31), ALLOCATABLE :: Species_Names(:)

  ! Work array to hold the list of all KPP species indices
  ! (Non-KPP species are given missing values)
  INTEGER,           ALLOCATABLE :: KppSpcId(:)

  ! Work array to hold the list of KPP fixed species indices
  ! (Non-KPP species are given missing values)
  INTEGER,           ALLOCATABLE :: KppFixId(:)

  ! Work array to hold the unique list of KPP variable species indices
  ! (Non-KPP species are given missing values)
  INTEGER,           ALLOCATABLE :: KppVarId(:)

CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Init_Species_Database
!
! !DESCRIPTION: Initializes the GEOS-Chem Species database from
!  YAML file format input.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Init_Species_Database( Input_Opt, SpcData, SpcCount, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE Input_Opt_Mod, ONLY : OptInput
    USE QFYAML_Mod
    USE RoundOff_Mod,  ONLY : Cast_and_Roundoff
    USE Species_Mod
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN)    :: Input_Opt    ! Input Options object
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(SpcPtr),   POINTER       :: SpcData(:)   ! Species database object
    TYPE(SpcIndCt), INTENT(INOUT) :: SpcCount     ! Species index counters
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC           ! Success/failure
!
! !REMARKS:
!  Uses the QFYAML parser, see: https://github.com/yantosca/qfyaml
!
! !REVISION HISTORY:
!  23 Apr 2020 - R. Yantosca - Initial version
!  See the subsequent Git history with the gitk browser!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    LOGICAL                     :: found_dd_dvzaersnow_luo
    LOGICAL                     :: found_dd_dvzminval_luo
    LOGICAL                     :: found_henry_cr_luo
    LOGICAL                     :: found_henry_k0_luo
    LOGICAL                     :: found_wd_convfaci2g_luo
    LOGICAL                     :: found_wd_kcscalefac_luo
    LOGICAL                     :: found_wd_liqandgas_luo
    LOGICAL                     :: found_wd_rainouteff_luo
    LOGICAL                     :: found_wd_retfactor_luo
    LOGICAL                     :: no_luo
    LOGICAL                     :: prtDebug
    LOGICAL                     :: v_bool
    LOGICAL                     :: wd_liqandgas_luo
    INTEGER                     :: v_int
    INTEGER                     :: nSpecies
    INTEGER                     :: N
    INTEGER                     :: S
    REAL(f4)                    :: v_real
    REAL(f4)                    :: dd_dvzaersnow_luo
    REAL(f4)                    :: henry_cr_luo
    REAL(f4)                    :: henry_k0_luo
    REAL(f4)                    :: wd_convfaci2g_luo
    REAL(f4)                    :: wd_retfactor_luo

    ! Strings
    CHARACTER(LEN=17)           :: tag
    CHARACTER(LEN=31)           :: spc
    CHARACTER(LEN=255)          :: v_str
    CHARACTER(LEN=255)          :: key
    CHARACTER(LEN=255)          :: thisLoc
    CHARACTER(LEN=512)          :: errMsg

    ! Arrays
    REAL(f4)                    :: a_real_2(2)
    REAL(f4)                    :: a_real_3(3)
    REAL(f4)                    :: dd_dvzminval_luo(2)
    REAL(f4)                    :: wd_kcscalefac_luo(3)
    REAL(f4)                    :: wd_rainouteff_luo(3)

    ! String arrays
    CHARACTER(LEN=17)           :: tags(48)

    ! Objects
    TYPE(QFYAML_t)              :: yml
    TYPE(Species),    POINTER   :: ThisSpc

    !=======================================================================
    ! Init_Species_Database begins here!
    !=======================================================================

    ! Initialize
    RC         = GC_SUCCESS
    prtDebug   = ( Input_Opt%LPRT .and. Input_Opt%amIRoot )
    errMsg     = ""
    thisLoc    = &
    " -> at Init_Species_Database (in module Headers/species_database_mod.F90"

    ! Zero counters
    SpcCount%nAdvect  = 0
    SpcCount%nAeroSpc = 0
    SpcCount%nDryAlt  = 0
    SpcCount%nDryDep  = 0
    SpcCount%nGasSpc  = 0
    SpcCount%nHygGrth = 0
    SpcCount%nKppVar  = 0
    SpcCount%nKppFix  = 0
    SpcCount%nKppSpc  = 0
    SpcCount%nOmitted = 0
    SpcCount%nPhotol  = 0
    SpcCount%nRadNucl = 0
    SpcCount%nRealSpc = 0
    SpcCount%nWetDep  = 0
    SpcCount%nHg0     = 0
    SpcCount%nHg2     = 0
    SpcCount%nHgP     = 0

    ! Species database tags to match
    tags = (/"Background_VV    ",  &
             "DD_AeroDryDep    ",  &
             "DD_DustDryDep    ",  &
             "DD_DvzAerSnow    ",  &
             "DD_DvzAerSnow_Luo",  &
             "DD_DvzMinVal     ",  &
             "DD_DvzMinVal_Luo ",  &
             "DD_F0            ",  &
             "DD_Hstar         ",  &
             "DD_KOA           ",  &
             "Density          ",  &
             "Formula          ",  &
             "FullName         ",  &
             "Is_Aerosol       ",  &
             "Is_DryAlt        ",  &
             "Is_DryDep        ",  &
             "Is_HygroGrowth   ",  &
             "Is_Gas           ",  &
             "Is_Hg0           ",  &
             "Is_Hg2           ",  &
             "Is_HgP           ",  &
             "Is_Photolysis    ",  &
             "Is_RadioNuclide  ",  &
             "Is_WetDep        ",  &
             "Henry_CR         ",  &
             "Henry_CR_Luo     ",  &
             "Henry_K0         ",  &
             "Henry_K0_Luo     ",  &
             "Henry_pKa        ",  &
             "MP_SizeResAer    ",  &
             "MP_SizeResNum    ",  &
             "MW_g             ",  &
             "Radius           ",  &
             "WD_AerScavEff    ",  &
             "WD_CoarseAer     ",  &
             "WD_ConvFacI2G    ",  &
             "WD_ConvFacI2G_Luo",  &
             "WD_KcScaleFac    ",  &
             "WD_KcScaleFac_Luo",  &
             "WD_Is_H2SO4      ",  &
             "WD_Is_HNO3       ",  &
             "WD_Is_SO2        ",  &
             "WD_LiqAndGas     ",  &
             "WD_LiqAndGas_Luo ",  &
             "WD_RainoutEff    ",  &
             "WD_RainoutEff_Luo",  &
             "WD_RetFactor     ",  &
             "WD_RetFactor_Luo "   /)

    !=======================================================================
    ! Store the list unique GEOS-Chem species names in work arrays for use
    ! below. This is the combined list of advected species (from
    ! geoschem_config.yml) plus KPP species (from SPC_NAMES in
    ! gckpp_Monitor.F90), with all duplicates removed. Also stores the
    ! corresponding indices in the KPP VAR and FIX arrays.  For simulations
    ! that do not use KPP, the unique species list is the list of advected
    ! species from geoschem_config.yml.
    !=======================================================================
    CALL Unique_Species_Names( Input_Opt, nSpecies, RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = "Could not determine species names!"
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    ! Initialize the species database vector and
    ! set all tags for each species to missing values
    CALL SpcData_Init( Input_Opt, nSpecies, SpcData, RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = "Could not initialize the species database object!"
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !=======================================================================
    ! Read the species metadata from YAML files into a QFYAML object
    !=======================================================================
    CALL Read_Species_Database( Input_Opt, yml, RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error encountered in routine "Read_Species_Database"!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !=======================================================================
    ! Extract species metadata and store in the Species Database object
    !=======================================================================

    ! Loop over the number of species
    DO S = 1, nSpecies

       ! Species name
       spc = species_names(S)

       !--------------------------------------------------------------------
       ! If the species is a "dummy" species (i.e. used for bookkeeping in
       ! KPP rxns), then flag it as such and skip to the next species
       !--------------------------------------------------------------------
       v_bool = MISSING_BOOL
       key    =  TRIM( spc ) // "%Is_Omitted"
       CALL QFYAML_Add_Get( yml, key, v_bool, "", RC )
       IF ( v_bool ) THEN
          SpcCount%nOmitted = SpcCount%nOmitted + 1
          CYCLE
       ENDIF

       !--------------------------------------------------------------------
       ! Look up this species in the database, and assign its name to the
       ! modelId field.  Subtract # of omitted species from the modelId.
       !--------------------------------------------------------------------
       N                 =  S - SpcCount%nOmitted
       ThisSpc           => SpcData(N)%Info
       ThisSpc%ModelId   =  N
       ThisSpc%Name      =  TRIM( spc )
       SpcCount%nRealSpc =  SpcCount%nRealSpc + 1

       !--------------------------------------------------------------------
       ! Set the Is_Advected tag (check against Input_Opt%AdvecSpc list)
       !-------------------------------------------------------------------
       v_bool = ANY( Input_Opt%AdvectSpc_Name == spc )
       IF ( v_bool ) THEN
          SpcCount%nAdvect    = SpcCount%nAdvect + 1
          ThisSpc%AdvectId    = SpcCount%nAdvect
          ThisSpc%Is_Advected = v_bool
       ENDIF

       !--------------------------------------------------------------------
       ! Set tags for species in the KPP mechanism
       !-------------------------------------------------------------------

       ! Is this species in the KPP mechanism?
       IF ( KppSpcId(S) > 0 ) THEN
          SpcCount%nKppSpc = SpcCount%nKppSpc + 1
          ThisSpc%KppSpcId = KppSpcId(S)
       ENDIF

       ! Is this species an active KPP species?
       IF ( KppVarId(S) > 0 ) THEN
          SpcCount%nKppVar = SpcCount%nKppVar + 1
          ThisSpc%KppVarId = KppVarId(S)
       ENDIF

       ! Is this species a fixed KPP species?
       IF ( KppFixId(S) > 0 ) THEN
          SpcCount%nKppFix = SpcCount%nKppFix + 1
          ThisSpc%KppFixId = KppFixId(S)
       ENDIF

       ! Is the species part of the KPP chemical mechanism?
       ThisSpc%Is_Kpp = ( ThisSpc%KppVarId > 0  .or. ThisSpc%KppFixId > 0 )

       ! Is the species an active or fixed species in the chemical mechanism?
       ThisSpc%Is_ActiveChem = ( ThisSpc%KppVarId >  0 .and.                 &
                                 ThisSpc%KppFixId <= 0                      )
       ThisSpc%Is_FixedChem  = ( ThisSpc%KppFixId >  0                      )

       !--------------------------------------------------------------------
       ! Initialize found flags
       !-------------------------------------------------------------------
       found_dd_dvzaersnow_luo = .FALSE.
       found_dd_dvzminval_luo  = .FALSE.
       found_henry_cr_luo      = .FALSE.
       found_henry_k0_luo      = .FALSE.
       found_wd_convfaci2g_luo = .FALSE.
       found_wd_kcscalefac_luo = .FALSE.
       found_wd_liqandgas_luo  = .FALSE.
       found_wd_rainouteff_luo = .FALSE.
       found_wd_retfactor_luo  = .FALSE.

       !--------------------------------------------------------------------
       ! Loop over the remaining tags in the species database and
       ! copy values from the QFYAML object to the SpcData object
       !--------------------------------------------------------------------
       DO N = 1, SIZE( tags )

          ! Set intial values to default "missing" values
          ! If the tag isn't found for a given species, then
          ! it will be given the appropriate missing value.
          a_real_2 = MISSING_REAL
          a_real_3 = MISSING_REAL
          v_bool   = MISSING_BOOL
          v_int    = MISSING_INT
          v_real   = MISSING_REAL
          v_str    = MISSING_STR

          ! Create search key for each variable
          key = TRIM( spc ) // '%' // TRIM( tags(N) )

          ! Set a flag if "Luo" is not found in the key
          no_luo = ( INDEX( key, "Luo" ) <= 0 )

          ! Save into the proper field of the species database
          ! NOTE: Attempt to round off values to 2 decimal places,
          ! unless the values can be either too large or too small
          ! for the roundoff algorithm.
          IF ( INDEX( key, "%Background_VV" ) > 0 ) THEN
             v_real = MISSING_VV
             CALL QFYAML_Add_Get( yml, key, v_real, "", RC )
             IF ( RC /= GC_SUCCESS ) GOTO 999
             ThisSpc%BackgroundVV = DBLE( v_real )   ! Don't round off

          ELSE IF ( INDEX( key, "%DD_AeroDryDep" ) > 0 ) THEN
             CALL QFYAML_Add_Get( yml, key, v_bool, "", RC )
             IF ( RC /= GC_SUCCESS ) GOTO 999
             ThisSpc%DD_AeroDryDep = v_bool

          ELSE IF ( INDEX( key, "%DD_DustDryDep" ) > 0 ) THEN
             CALL QFYAML_Add_Get( yml, key, v_bool, "", RC )
             IF ( RC /= GC_SUCCESS ) GOTO 999
             ThisSpc%DD_DustDryDep = v_bool

          ELSE IF ( INDEX( key, "%DD_DvzAerSnow" ) >  0  .and. no_luo ) THEN
             CALL QFYAML_Add_Get( yml, key, v_real, "", RC )
             IF ( RC /= GC_SUCCESS ) GOTO 999
             ThisSpc%DD_DvzAerSnow = Cast_and_RoundOff( v_real, 2 )

          ELSE IF ( INDEX( key, "%DD_DvzAerSnow_Luo" ) > 0 ) THEN
             CALL QFYAML_Add_Get( yml, key, v_real, "", RC )
             IF ( RC /= GC_SUCCESS ) GOTO 999
             dd_dvzaersnow_luo = Cast_and_RoundOff( v_real, 2 )
             IF ( dd_dvzaersnow_luo /= MISSING_REAL ) THEN
                found_dd_dvzaersnow_luo = .TRUE.
             ENDIF

          ELSE IF ( INDEX( key, "%DD_DvzMinVal" ) > 0 .and. no_luo ) THEN
             CALL QFYAML_Add_Get( yml, key, a_real_2, "", RC )
             IF ( RC /= GC_SUCCESS ) GOTO 999
             ThisSpc%DD_DvzMinVal(1) = Cast_and_RoundOff( a_real_2(1), 2 )
             ThisSpc%DD_DvzMinVal(2) = Cast_and_RoundOff( a_real_2(2), 2 )

          ELSE IF ( INDEX( key, "%DD_DvzMinVal_Luo" ) > 0 ) THEN
             CALL QFYAML_Add_Get( yml, key, a_real_2, "", RC )
             IF ( RC /= GC_SUCCESS ) GOTO 999
             dd_dvzminval_luo(1) = Cast_and_RoundOff( a_real_2(1), 2 )
             dd_dvzminval_luo(2) = Cast_and_RoundOff( a_real_2(2), 2 )
             IF ( dd_dvzminval_luo(1) /= MISSING_REAL ) THEN
                found_dd_dvzminval_luo = .TRUE.
             ENDIF

          ELSE IF ( INDEX( key, "%DD_F0" ) > 0 ) THEN
             CALL QFYAML_Add_Get( yml, key, v_real, "", RC )
             IF ( RC /= GC_SUCCESS ) GOTO 999
             ThisSpc%DD_F0 = DBLE( v_real )          ! Don't round off

          ELSE IF ( INDEX( key, "%DD_Hstar" ) > 0 ) THEN
             CALL QFYAML_Add_Get( yml, key, v_real, "", RC )
             IF ( RC /= GC_SUCCESS ) GOTO 999
             ThisSpc%DD_Hstar = DBLE( v_real )       ! Don't round off

          ELSE IF ( INDEX( key, "%DD_KOA" ) > 0 ) THEN
             CALL QFYAML_Add_Get( yml, key, v_real, "", RC )
             IF ( RC /= GC_SUCCESS ) GOTO 999
             ThisSpc%DD_KOA = DBLE( v_real )       ! Don't round off

          ELSE IF ( INDEX( key, "%Density" ) > 0 ) THEN
             CALL QFYAML_Add_Get( yml, key, v_real, "", RC )
             IF ( RC /= GC_SUCCESS ) GOTO 999
             ThisSpc%Density = Cast_and_RoundOff( v_real, 2 )

          ELSE IF ( INDEX( key, "%Formula" ) > 0 ) THEN
             CALL QFYAML_Add_Get( yml, key, v_str, "", RC )
             IF ( RC /= GC_SUCCESS ) GOTO 999
             ThisSpc%Formula = TRIM( v_str )

          ELSE IF ( INDEX( key, "%FullName" ) > 0 ) THEN
             CALL QFYAML_Add_Get( yml, key, v_str, "", RC )
             IF ( RC /= GC_SUCCESS ) GOTO 999
             ThisSpc%FullName = TRIM( v_str )

          ELSE IF ( INDEX( key, "%Henry_CR" ) > 0 .and. no_luo ) THEN
             CALL QFYAML_Add_Get( yml, key, v_real, "", RC )
             IF ( RC /= GC_SUCCESS ) GOTO 999
             ThisSpc%Henry_CR = DBLE( v_real )       ! Don't round off

          ELSE IF ( INDEX( key, "%Henry_CR_Luo" ) > 0 ) THEN
             CALL QFYAML_Add_Get( yml, key, v_real, "", RC )
             IF ( RC /= GC_SUCCESS ) GOTO 999
             henry_cr_luo = DBLE( v_real )           ! Don't round off
             IF ( henry_cr_luo /= MISSING_REAL ) THEN
                found_henry_cr_luo = .TRUE.
             ENDIF

          ELSE IF ( INDEX( key, "%Henry_K0" ) > 0 .and. no_luo ) THEN
             CALL QFYAML_Add_Get( yml, key, v_real, "", RC )
             IF ( RC /= GC_SUCCESS ) GOTO 999
             ThisSpc%Henry_K0 = DBLE( v_real )       ! Don't round off

          ELSE IF ( INDEX( key, "%Henry_K0_Luo" ) > 0 ) THEN
             CALL QFYAML_Add_Get( yml, key, v_real, "", RC )
             IF ( RC /= GC_SUCCESS ) GOTO 999
             henry_k0_luo = DBLE( v_real )           ! Don't round off
             IF ( henry_k0_luo /= MISSING_REAL ) THEN
                found_henry_k0_luo = .TRUE.
             ENDIF

          ELSE IF ( INDEX( key, "%Is_Aerosol" ) > 0  ) THEN
             CALL QFYAML_Add_Get( yml, key, v_bool, "", RC )
             IF ( RC /= GC_SUCCESS ) GOTO 999
             IF ( v_bool ) THEN
                SpcCount%nAeroSpc  = SpcCount%nAeroSpc + 1
                ThisSpc%AerosolId  = SpcCount%nAeroSpc
                ThisSpc%Is_Aerosol = v_bool
             ENDIF

          ELSE IF ( INDEX( key, "%Is_DryAlt" ) > 0 ) THEN
             CALL QFYAML_Add_Get( yml, key, v_bool, "", RC )
             IF ( RC /= GC_SUCCESS ) GOTO 999
             IF ( v_bool ) THEN
                SpcCount%nDryAlt  = SpcCount%nDryAlt + 1
                ThisSpc%DryAltId  = SpcCount%nDryAlt
                ThisSpc%Is_DryAlt = v_bool
             ENDIF

          ELSE IF ( INDEX( key, "%Is_DryDep" ) > 0 ) THEN
             CALL QFYAML_Add_Get( yml, key, v_bool, "", RC )
             IF ( RC /= GC_SUCCESS ) GOTO 999
             IF ( v_bool .AND. ThisSpc%Is_Advected ) THEN
                SpcCount%nDryDep  = SpcCount%nDryDep + 1
                ThisSpc%DryDepId  = SpcCount%nDryDep
                ThisSpc%Is_DryDep = v_bool
             ENDIF

          ELSE IF ( INDEX( key, "%Is_HygroGrowth" ) > 0 ) THEN
             CALL QFYAML_Add_Get( yml, key, v_bool, "", RC )
             IF ( RC /= GC_SUCCESS ) GOTO 999
             IF ( v_bool ) THEN
                SpcCount%nHygGrth      = SpcCount%nHygGrth + 1
                ThisSpc%HygGrthId      = SpcCount%nHygGrth
                ThisSpc%Is_HygroGrowth = v_bool
             ENDIF

          ELSE IF ( INDEX( key, "%Is_Gas" ) > 0 ) THEN
             CALL QFYAML_Add_Get( yml, key, v_bool, "", RC )
             IF ( RC /= GC_SUCCESS ) GOTO 999
             IF ( v_bool ) THEN
                SpcCount%nGasSpc = SpcCount%nGasSpc + 1
                ThisSpc%GasSpcId = SpcCount%nGasSpc
                ThisSpc%Is_Gas   = v_bool
             ENDIF

          ELSE IF ( INDEX( key, "%Is_Hg0" ) > 0 ) THEN
             CALL QFYAML_Add_Get( yml, key, v_bool, "", RC )
             IF ( RC /= GC_SUCCESS ) GOTO 999
             IF ( v_bool ) THEN
                SpcCount%nHg0  = SpcCount%nHg0 + 1
                ThisSpc%Is_Hg0 = v_bool
             ENDIF

          ELSE IF ( INDEX( key, "%Is_Hg2" ) > 0 ) THEN
             CALL QFYAML_Add_Get( yml, key, v_bool, "", RC )
             IF ( RC /= GC_SUCCESS ) GOTO 999
             IF ( v_bool ) THEN
                SpcCount%nHg2  = SpcCount%nHg2 + 1
                ThisSpc%Is_Hg2 = v_bool
             ENDIF

          ELSE IF ( INDEX( key, "%Is_HgP" ) > 0 ) THEN
             CALL QFYAML_Add_Get( yml, key, v_bool, "", RC )
             IF ( RC /= GC_SUCCESS ) GOTO 999
             IF ( v_bool ) THEN
                SpcCount%nHgP  = SpcCount%nHgP + 1
                ThisSpc%Is_HgP = v_bool
             ENDIF

          ELSE IF ( INDEX( key, "%Is_Photolysis" ) > 0 ) THEN
             CALL QFYAML_Add_Get( yml, key, v_bool, "", RC )
             IF ( RC /= GC_SUCCESS ) GOTO 999
             IF ( v_bool ) THEN
                SpcCount%nPhotol      = SpcCount%nPhotol + 1
                ThisSpc%PhotolId      = SpcCount%nPhotol
                ThisSpc%Is_Photolysis = v_bool
             ENDIF

          ELSE IF ( INDEX( key, "%Is_RadioNuclide" ) > 0 ) THEN
             CALL QFYAML_Add_Get( yml, key, v_bool, "", RC )
             IF ( RC /= GC_SUCCESS ) GOTO 999
             IF ( v_bool ) THEN
                SpcCount%nRadNucl       = SpcCount%nRadNucl + 1
                ThisSpc%RadNuclId       = SpcCount%nRadNucl
                ThisSpc%Is_RadioNuclide = v_bool
             ENDIF

          ELSE IF ( INDEX( key, "%Is_WetDep" ) > 0 ) THEN
             CALL QFYAML_Add_Get( yml, key, v_bool, "", RC )
             IF ( RC /= GC_SUCCESS ) GOTO 999
             IF ( v_bool ) THEN
                SpcCount%nWetDep  = SpcCount%nWetDep + 1
                ThisSpc%WetDepID  = SpcCount%nWetDep
                ThisSpc%Is_WetDep = v_bool
             ENDIF

          ELSE IF ( INDEX( key, "%MP_SizeResAer" ) > 0 ) THEN
             CALL QFYAML_Add_Get( yml, key, v_bool, "", RC )
             IF ( RC /= GC_SUCCESS ) GOTO 999
             ThisSpc%MP_SizeResAer = v_bool

          ELSE IF ( INDEX( key, "%MP_SizeResNum" ) > 0 ) THEN
             CALL QFYAML_Add_Get( yml, key, v_bool, "", RC )
             IF ( RC /= GC_SUCCESS ) GOTO 999
             ThisSpc%MP_SizeResNum = v_bool

          ELSE IF ( INDEX( key, "%MW_g" ) > 0 ) THEN
             CALL QFYAML_Add_Get( yml, key, v_real, "", RC )
             IF ( RC /= GC_SUCCESS ) GOTO 999
             ThisSpc%MW_g = Cast_and_RoundOff( v_real, 2 )

          ELSE IF ( INDEX( key, "%Radius" ) > 0 ) THEN
             CALL QFYAML_Add_Get( yml, key, v_real, "", RC )
             IF ( RC /= GC_SUCCESS ) GOTO 999
             ThisSpc%Radius = DBLE( v_real )         ! Don't round off

          ELSE IF ( INDEX( key, "%WD_AerScavEff" ) > 0 ) THEN
             CALL QFYAML_Add_Get( yml, key, v_real, "", RC )
             IF ( RC /= GC_SUCCESS ) GOTO 999
             ThisSpc%WD_AerScavEff = Cast_and_RoundOff( v_real, 2 )

          ELSE IF ( INDEX( key, "%WD_CoarseAer" ) > 0 ) THEN
             CALL QFYAML_Add_Get( yml, key, v_bool, "", RC )
             IF ( RC /= GC_SUCCESS ) GOTO 999
             ThisSpc%WD_CoarseAer = v_bool

          ELSE IF ( INDEX( key, "%WD_ConvFacI2G" ) > 0 .and. no_luo ) THEN
             CALL QFYAML_Add_Get( yml, key, v_real, "", RC )
             IF ( RC /= GC_SUCCESS ) GOTO 999
             ThisSpc%WD_ConvFacI2G = DBLE( v_real )  ! Don't round off

          ELSE IF ( INDEX( key, "%WD_ConvFacI2G_Luo" ) > 0 ) THEN
             CALL QFYAML_Add_Get( yml, key, v_real, "", RC )
             IF ( RC /= GC_SUCCESS ) GOTO 999
             wd_convfaci2g_luo = DBLE( v_real )      ! Don't round off
             IF ( wd_convfaci2g_luo /= MISSING_REAL ) THEN
                found_wd_convfaci2g_luo = .TRUE.
             ENDIF

          ELSE IF ( INDEX( key, "%WD_KcScaleFac" ) > 0  .and. no_luo ) THEN
             CALL QFYAML_Add_Get( yml, key, a_real_3, "", RC )
             IF ( RC /= GC_SUCCESS ) GOTO 999
             ThisSpc%WD_KcScaleFac(1) = Cast_and_RoundOff( a_real_3(1), 2 )
             ThisSpc%WD_KcScaleFac(2) = Cast_and_RoundOff( a_real_3(2), 2 )
             ThisSpc%WD_KcScaleFac(3) = Cast_and_RoundOff( a_real_3(3), 2 )

          ELSE IF ( INDEX( key, "%WD_KcScaleFac_Luo" ) > 0 ) THEN
             CALL QFYAML_Add_Get( yml, key, a_real_3, "", RC )
             IF ( RC /= GC_SUCCESS ) GOTO 999
             wd_kcscalefac_luo(1) = Cast_and_RoundOff( a_real_3(1), 2 )
             wd_kcscalefac_luo(2) = Cast_and_RoundOff( a_real_3(2), 2 )
             wd_kcscalefac_luo(3) = Cast_and_RoundOff( a_real_3(3), 2 )
             IF ( wd_kcscalefac_luo(1) /= MISSING_REAL ) THEN
                found_wd_kcscalefac_luo = .TRUE.
             ENDIF

          ELSE IF ( INDEX( key, "%WD_Is_H2SO4" ) > 0 ) THEN
             CALL QFYAML_Add_Get( yml, key, v_bool, "", RC )
             IF ( RC /= GC_SUCCESS ) GOTO 999
             ThisSpc%WD_Is_H2SO4 = v_bool

          ELSE IF ( INDEX( key, "%WD_Is_HNO3" ) > 0 ) THEN
             CALL QFYAML_Add_Get( yml, key, v_bool, "", RC )
             IF ( RC /= GC_SUCCESS ) GOTO 999
             ThisSpc%WD_Is_HNO3 = v_bool

          ELSE IF ( INDEX( key, "%WD_Is_SO2" ) > 0 ) THEN
             CALL QFYAML_Add_Get( yml, key, v_bool, "", RC )
             IF ( RC /= GC_SUCCESS ) GOTO 999
             ThisSpc%WD_Is_SO2 = v_bool

          ELSE IF ( INDEX( key, "%WD_LiqAndGas" ) > 0 .and. no_luo ) THEN
             CALL QFYAML_Add_Get( yml, key, v_bool, "", RC )
             IF ( RC /= GC_SUCCESS ) GOTO 999
             ThisSpc%WD_LiqAndGas = v_bool

          ELSE IF ( INDEX( key, "%WD_LiqAndGas_Luo" ) > 0 ) THEN
             CALL QFYAML_Add_Get( yml, key, v_bool, "", RC )
             IF ( RC /= GC_SUCCESS ) GOTO 999
             wd_liqandgas_luo = v_bool
             found_wd_liqandgas_luo = wd_liqandgas_luo

          ELSE IF ( INDEX( key, "%WD_RainoutEff" ) > 0 .and. no_luo ) THEN
             CALL QFYAML_Add_Get( yml, key, a_real_3, "", RC )
             IF ( RC /= GC_SUCCESS ) GOTO 999
             ThisSpc%WD_RainoutEff(1) = Cast_and_RoundOff( a_real_3(1), 2 )
             ThisSpc%WD_RainoutEff(2) = Cast_and_RoundOff( a_real_3(2), 2 )
             ThisSpc%WD_RainoutEff(3) = Cast_and_RoundOff( a_real_3(3), 2 )

          ELSE IF ( INDEX( key, "%WD_RainoutEff_Luo" ) > 0 ) THEN
             CALL QFYAML_Add_Get( yml, key, a_real_3, "", RC )
             IF ( RC /= GC_SUCCESS ) GOTO 999
             wd_rainouteff_luo(1) = Cast_and_RoundOff( a_real_3(1), 2 )
             wd_rainouteff_luo(2) = Cast_and_RoundOff( a_real_3(2), 2 )
             wd_rainouteff_luo(3) = Cast_and_RoundOff( a_real_3(3), 2 )
             IF ( wd_rainouteff_luo(1) /= MISSING_REAL ) THEN
                found_wd_rainouteff_luo = .TRUE.
             ENDIF

          ELSE IF ( INDEX( key, "%WD_RetFactor" ) > 0 .and. no_luo ) THEN
             CALL QFYAML_Add_Get( yml, key, v_real, "", RC )
             IF ( RC /= GC_SUCCESS ) GOTO 999
             ThisSpc%WD_RetFactor = Cast_and_RoundOff( v_real, 2 )

          ELSE IF ( INDEX( key, "%WD_RetFactor_Luo" ) > 0 ) THEN
             CALL QFYAML_Add_Get( yml, key, v_real, "", RC )
             IF ( RC /= GC_SUCCESS ) GOTO 999
             wd_retfactor_luo = Cast_and_RoundOff( v_real, 2 )
             IF ( wd_retfactor_luo /= MISSING_REAL ) THEN
                found_wd_retfactor_luo = .TRUE.
             ENDIF

          ELSE
             ! Pass

          ENDIF

       ENDDO

       !--------------------------------------------------------------------
       ! SANITY CHECKS
       !--------------------------------------------------------------------

       ! Is_Gas and Is_Aero tags cannot both be TRUE at the same time
       IF ( ThisSpc%Is_Gas .and. ThisSpc%Is_Aerosol ) THEN
          errMsg = "Is_Gas and Is_Aerosol are both TRUE for species "     // &
                   TRIM( spc ) // "!"
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       ! Is_Gas and Is_Aero tags cannot both be FALSE at the same time
       IF ( ( .not. ThisSpc%Is_Gas     )   .and.                             &
            ( .not. ThisSpc%Is_Aerosol )   .and.                             &
            ( .not. ThisSpc%Is_Omitted ) ) THEN

          ! Check if this is a KPP species, is so set Is_Gas to TRUE, otherwise
          ! return with an error. This will account for P/L families not
          ! defined in the species database.
          IF ( ThisSpc%Is_Kpp ) THEN
             ThisSpc%Is_Gas = .TRUE.
          ELSE
             errMsg = "Is_Gas and Is_Aerosol are both FALSE for species " // &
                      TRIM( spc ) // "!"
             CALL GC_Error( errMsg, RC, thisLoc )
             RETURN
          ENDIF
       ENDIF

       ! If the species is a gas, set all aerosol fields to missing values
       IF ( ThisSpc%Is_Gas ) THEN

          SELECT CASE( TRIM( spc ) )
             CASE( 'H2SO4' )
                ! H2SO4 are gases that wetdep like aerosols,
                ! so keep both all gas and aerosol properties.
             CASE( 'HNO3', 'SO2' )
                ! HNO3 and SO2 drydep like gases but wetdep like fine
                ! aerosols, so set certain fields to missing values.
                ThisSpc%DD_DvzAerSnow = MISSING
                ThisSpc%MP_SizeResAer = MISSING_BOOL
                ThisSpc%MP_SizeResNum = MISSING_BOOL
                ThisSpc%WD_CoarseAer  = MISSING_BOOL
             CASE DEFAULT
                ! For all other gas-phase species, set all
                ! aerosol fields to missing values
                ThisSpc%DD_DvzAerSnow = MISSING
                ThisSpc%MP_SizeResAer = MISSING_BOOL
                ThisSpc%MP_SizeResNum = MISSING_BOOL
                ThisSpc%WD_CoarseAer  = MISSING_BOOL
                ThisSpc%WD_AerScavEff = MISSING
                ThisSpc%WD_KcScaleFac = MISSING
                ThisSpc%WD_RainoutEff = MISSING
          END SELECT
       ENDIF

       ! If the species is an aerosol, set all gas fields to missing values
       IF ( ThisSpc%Is_Aerosol ) THEN
          ThisSpc%WD_ConvFacI2G = MISSING
          ThisSpc%WD_RetFactor  = MISSING
          ThisSpc%WD_LiqAndGas  = MISSING_BOOL
       ENDIF

#ifdef LUO_WETDEP
       !--------------------------------------------------------------------
       ! For Luo et al 2020 wetdep
       ! Overwrite with special values if present in file
       !--------------------------------------------------------------------
       IF ( found_dd_dvzaersnow_luo ) THEN
          ThisSpc%DD_DvzAerSnow = dd_dvzaersnow_luo
       ENDIF

       IF ( found_dd_dvzminval_luo ) THEN
          ThisSpc%DD_DvzMinVal(1) = dd_dvzminval_luo(1)
          ThisSpc%DD_DvzMinVal(2) = dd_dvzminval_luo(2)
       ENDIF

       IF ( found_henry_cr_luo ) THEN
          ThisSpc%Henry_CR = henry_cr_luo
       ENDIF

       IF ( found_henry_k0_luo ) THEN
          ThisSpc%Henry_K0 = henry_k0_luo
       ENDIF

       IF ( found_wd_convfaci2g_luo ) THEN
          ThisSpc%WD_ConvFacI2G = wd_convfaci2g_luo
       ENDIF

       IF ( found_wd_liqandgas_luo ) THEN
          ThisSpc%WD_LiqAndGas = wd_liqandgas_luo
       ENDIF

       IF ( found_wd_kcscalefac_luo ) THEN
          ThisSpc%WD_KcScaleFac(1) = wd_kcscalefac_luo(1)
          ThisSpc%WD_KcScaleFac(2) = wd_kcscalefac_luo(2)
          ThisSpc%WD_KcScaleFac(3) = wd_kcscalefac_luo(3)
       ENDIF

       IF ( found_wd_rainouteff_luo ) THEN
          ThisSpc%WD_RainoutEff(1) = wd_rainouteff_luo(1)
          ThisSpc%WD_RainoutEff(2) = wd_rainouteff_luo(2)
          ThisSpc%WD_RainoutEff(3) = wd_rainouteff_luo(3)
       ENDIF

       IF ( found_wd_retfactor_luo ) THEN
          ThisSpc%WD_RetFactor = wd_retfactor_luo
       ENDIF
#endif

       ! Free pointer
       ThisSpc => NULL()
    ENDDO

    ! FORMAT statements
10  FORMAT( a30, " | ", a      )
20  FORMAT( a30, " | ", L10    )
30  FORMAT( a30, " | ", f10.2  )
31  FORMAT( a30, " | ", 2f10.2 )
32  FORMAT( a30, " | ", 3f10.2 )
40  FORMAT( a30, " | ", i10    )

    !=======================================================================
    ! Print metadata for only the species that are defined in this
    ! simulation (but not the entire species database) to a YAML file.
    ! This file may be used for pre-processing files in other models
    ! when updating GEOS-Chem versions, such as in WRF and CESM. It
    ! should not be generated when running those models. Output file is
    ! set in simulation%species_metadata_output_file in geoschem_config.yml.
    !=======================================================================
    IF ( LEN(TRIM( Input_Opt%SpcMetaDataOutFile )) > 0 ) THEN
       IF ( Input_Opt%amIRoot ) THEN
          CALL QFYAML_Print( yml        = yml,                               &
                             fileName   = Input_Opt%SpcMetaDataOutFile,      &
                             searchKeys = species_names,                     &
                             RC         = RC                                )

       ENDIF
    ENDIF

    !=======================================================================
    ! Normal exit
    !=======================================================================

    ! Free objects and arrays, then return
    ThisSpc => NULL()
    CALL QFYAML_CleanUp( yml )
    CALL Cleanup_Work_Arrays()

    !### Uncomment this to stop here when debugging species info
    !STOP

    RETURN

    !=======================================================================
    ! Abnormal exit
    !=======================================================================
999 CONTINUE

    ! Free objects and arrays
    ThisSpc => NULL()
    CALL QFYAML_CleanUp( yml )
    CALL Cleanup_Work_Arrays()

    ! Exit with error
    errMsg = 'Could not read species database variable: ' // TRIM( key )
    CALL GC_Error( errMsg, RC, thisLoc )

  END SUBROUTINE Init_Species_Database
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Read_Species_Database
!
! !DESCRIPTION: Reads the metadata for each species into a QFYAML object.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Read_Species_Database( Input_Opt, yml, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE Input_Opt_Mod, ONLY : OptInput
    USE QFYAML_Mod
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN) :: Input_Opt   ! Input options
!
! !OUTPUT PARAMETERS:
!
    TYPE(QFYAML_t), INTENT(OUT) :: yml
    INTEGER,        INTENT(OUT) :: RC
!
! !REVISION HISTORY:
!  28 Apr 2020 - R. Yantosca - Initial version
!  See the subsequent Git history with the gitk browser!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Strings
    CHARACTER(LEN=255) :: fileName
    CHARACTER(LEN=255) :: thisLoc
    CHARACTER(LEN=512) :: errMsg

    ! Objects
    TYPE(QFYAML_t)     :: yml_1
    TYPE(QFYAML_t)     :: yml_2
    TYPE(QFYAML_t)     :: yml_anchored

    !=========================================================================
    ! Read_Species_Database begins here!
    !=========================================================================

    ! Initialize
    RC      = GC_SUCCESS
    errMsg  = ""
    thisLoc = &
    " -> at Read Species Database (in module Headers/species_database_mod.F90)"

    !=======================================================================
    ! Read metadata for GEOS-Chem species
    !========================================================================
    fileName = TRIM( Input_Opt%SpcDatabaseFile )
    CALL QFYAML_Init( fileName, yml, yml_anchored, RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = "Error reading " // TRIM( fileName )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    CALL QFYAML_CleanUp( yml_anchored )

  END SUBROUTINE Read_Species_Database
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Cleanup_Species_Database
!
! !DESCRIPTION: Finalizes the vector with species information.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Cleanup_Species_Database( SpcData, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE Species_Mod
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(SpcPtr),   POINTER     :: SpcData(:)   ! Species database object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT) :: RC           ! Success or failure?
!
! !REMARKS:
!
! !REVISION HISTORY:
!  22 Jul 2015 - R. Yantosca - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Assume success
    RC = GC_SUCCESS

    ! Deallocate the species database object
    CALL SpcData_Cleanup( SpcData )

  END SUBROUTINE Cleanup_Species_Database
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: TranUc
!
! !DESCRIPTION: Tranlate a character variable to all upper case letters.
!  Non-alphabetic characters are not affected.  The original "text" is
!  destroyed.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE TranUc( text )
!
! !INPUT/OUTPUT PARAMETERS:
!
    CHARACTER(LEN=*), INTENT(INOUT) :: text
!
! !AUTHOR:
!  Robert D. Stewart, May 19, 1992 (part of CHARPAK)
!
! !REMARKS:
!  Keep a private shadow copy of this routine here so as not to
!  incur a dependency with GeosUtil/charpak_mod.F90.  This lets us
!  keep species_datbase_mod.F90 in the Headers/ folder together
!  with state_chm_mod.F90 and species_mod.F90.
!
! !REVISION HISTORY:
!  06 Jan 2015 - R. Yantosca - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER :: iasc, i, ilen

    ilen = LEN(text)
    DO i=1,ilen
       iasc = ICHAR(text(i:i))
       IF ((iasc.GT.96).AND.(iasc.LT.123)) THEN
          text(i:i) = CHAR(iasc-32)
       ENDIF
    ENDDO

  END SUBROUTINE TranUc
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Unique_Species_Names
!
! !DESCRIPTION: Stores the list of unique species names (i.e. removing
!  duplicates from the list of advected species and the the list of KPP
!  species) for later use.  Also computes the corresponding indices for
!  the KPP variable and fixed species arrays (VAR and FIX, respectively).
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Unique_Species_Names( Input_Opt, nSpecies, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE Input_Opt_Mod,    ONLY : OptInput
    USE GcKpp_Monitor,    ONLY : Spc_Names
    USE GcKpp_Parameters, ONLY : NFIX, NSPEC, NVAR
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN)  :: Input_Opt   ! Input Options object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT) :: nSpecies    ! Number of unique species
    INTEGER,        INTENT(OUT) :: RC          ! Success or failure
!
! !REMARKS:
!  This may not be the fastest search algorithm (because it relies on string
!  comparisons).  But it is only executed at startup so we can live with it.
!  We could make it faster by hashing but that seems like overkill.
!
! !REVISION HISTORY:
!  09 May 2016 - R. Yantosca - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER                        :: nAdvect, K, S

    ! Strings
    CHARACTER(LEN=255)             :: errMsg
    CHARACTER(LEN=255)             :: thisLoc

    ! Arrays
    CHARACTER(LEN=31), ALLOCATABLE :: Tmp(:)
    CHARACTER(LEN=31)              :: SpcName

    !=======================================================================
    ! UNIQUE_SPECIES_NAMES begins here!
    !=======================================================================

    ! Assume success
    RC       = GC_SUCCESS
    errMsg   = ''
    thisLoc  = &
    ' -> at Unique_Species_Names (in module Headers/species_database_mod.F90)'

    ! Number of advected species listed in geoschem_config.yml
    nAdvect  = Input_Opt%N_Advect

    ! First set the # of species to the # of advected species
    nSpecies = nAdvect

    !=======================================================================
    ! For KPP-based simulations, get the list of all of
    ! species names in the KPP mechanism, and their indices
    !=======================================================================
    IF ( Input_Opt%ITS_A_FULLCHEM_SIM      .or.                             &
         Input_Opt%ITS_A_MERCURY_SIM       .or.                             &
         Input_Opt%ITS_A_CARBON_SIM      ) THEN

       ! Allocate a temporary array large enough to hold all of the
       ! advected species listed in geoschem_config.yml as well as all of the
       ! KPP species names (listed in SPC_NAMES of gckpp_Monitor.F90)
       ALLOCATE( Tmp( nAdvect + NSPEC ), STAT=RC )
       CALL GC_CheckVar( 'species_database_mod.F90:Tmp', 0 , RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       Tmp = ''

       !--------------------------------------------------------------------
       ! First determine the unique list of species in the KPP mechanism
       ! (so that we don't duplicate storage for advected & chemical species)
       !--------------------------------------------------------------------

       ! First, store advected species (from geoschem_config.yml) in the
       ! TMP array
       DO S = 1, nSpecies
          Tmp(S) = Input_Opt%AdvectSpc_Name(S)
       ENDDO

       ! Loop over KPP species
       DO K = 1, NSPEC

          ! Skip dummy RR species for prod/loss diagnostic (mps, 8/23/16)
          SpcName = ADJUSTL( Spc_Names(K) )
          IF ( SpcName(1:2) == 'RR' ) CYCLE

          ! Next, add to the TMP array those KPP species that aren't already
          ! listed as advected species.  nSpecies is the # of unique species.
          IF ( .not. ANY( Input_Opt%AdvectSpc_Name == Spc_Names(K) ) ) THEN
             nSpecies      = nSpecies + 1
             Tmp(nSpecies) = Spc_Names(K)
          ENDIF

       ENDDO

       ! Allocate the species names array precisely of length nSpecies
       ALLOCATE( Species_Names( nSpecies ), STAT=RC )
       CALL GC_CheckVar( 'species_database_mod.F90:Species_Names', 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       Species_Names = Tmp(1:nSpecies )

       ! Free temporary array
       IF ( ALLOCATED( Tmp ) ) DEALLOCATE( Tmp )

       !--------------------------------------------------------------------
       ! Now determine the KPP indices for each unique species name
       !--------------------------------------------------------------------

       ! Work array to hold the list of all KPP species indices
       ALLOCATE( KppSpcId( nSpecies ), STAT=RC )
       CALL GC_CheckVar( 'species_database_mod.F90:KppSpcId', 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       KppSpcId = MISSING_INT

       ! Work array to hold the list of KPP fixed species indices
       ALLOCATE( KppFixId( nSpecies ), STAT=RC )
       CALL GC_CheckVar( 'species_database_mod.F90:KppFixId', 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       KppFixId = MISSING_INT

       ! Work array to hold the list of KPP variable species indices
       ALLOCATE( KppVarId( nSpecies ), STAT=RC )
       CALL GC_CheckVar( 'species_database_mod.F90:KppVarId', 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       KppVarId = MISSING_INT

       ! Loop through the list of unique species names
       DO S = 1, nSpecies

          ! Loop through the list of KPP species (stored in SPC_NAMES)
          DO K = 1, NSPEC

             ! Skip dummy RR species for prod/loss diagnostic (mps, 8/23/16)
             SpcName = ADJUSTL( Spc_Names(K) )
             IF ( SpcName(1:2) == 'RR' ) CYCLE

             ! Test the unique species names (stored in SPECIES_NAMES)
             ! against the list of KPP species (in SPC_NAMES).  The K
             ! index corresponds to the location of the species in the
             ! KPP chemical mechanism:  1..NSPEC = [ 1..NVAR, 1..NFIX].
             IF ( Species_Names(S) == Spc_Names(K) ) THEN

                ! KPP species index (1..NSPEC).  These
                ! are used to index species in the KPP "C" array.
                ! These include both variable and fixed species.
                KppSpcId(S) = K

                IF ( K <= NVAR ) THEN

                   ! KPP variable species index (1..NVAR).  These
                   ! are used to index species in the KPP "C" array
                   ! (as well as the "VAR" array).
                   KppVarId(S) = K

                ELSE

                   ! KPP fixed species also have entries (1..NFIX).  These
                   ! are used to index species in the KPP "FIX" array.
                   KppFixId(S) = K - NVAR

                ENDIF

                ! Skip to next species
                EXIT
             ENDIF
          ENDDO
       ENDDO

    !=======================================================================
    ! For specialty simulations, we do not have KPP species.  Thus, the
    ! of species is just the list of advected species from geoschem_config.yml
    !=======================================================================
    ELSE

       ! Initialize the species names array from Input_Opt
       ALLOCATE( Species_Names( nSpecies ), STAT=RC )
       CALL GC_CheckVar( 'species_database_mod.F90:Species_Names', 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       Species_Names = Input_Opt%AdvectSpc_Name(1:nSpecies)

       ! Set KppSpcId to missing value
       ALLOCATE( KppSpcId( nSpecies ), STAT=RC )
       CALL GC_CheckVar( 'species_database_mod.F90:KppSpcId', 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       KppSpcId = MISSING_INT

       ! Set KppFixId to missing value
       ALLOCATE( KppFixId( nSpecies ), STAT=RC )
       CALL GC_CheckVar( 'species_database_mod.F90:KppFixId', 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       KppFixId = MISSING_INT

       ! Set KppVarId to missing value
       ALLOCATE( KppVarId( nSpecies ), STAT=RC )
       CALL GC_CheckVar( 'species_database_mod.F90:KppVarId', 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       KppVarId = MISSING_INT

    ENDIF

  END SUBROUTINE Unique_Species_Names
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Cleanup_Work_Arrays
!
! !DESCRIPTION: Cleans working (temporary) arrays used by this module,
!  restoring them to an unused state. It is called at the end of
!  Init\_Species\_Database or by an external module when needed to
!  reinitialize the species DB.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Cleanup_Work_Arrays()
!
! !REMARKS:
!  This routine allows Species_Database_Mod to be initialized more than once
!  in the same CPU, if called externally before re-initializing a State_Chm
!  derived type object.
!
! !REVISION HISTORY:
!  06 May 2016 - R. Yantosca - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Deallocate arrays
    IF ( ALLOCATED( Species_Names ) ) DEALLOCATE( Species_Names )
    IF ( ALLOCATED( KppFixId      ) ) DEALLOCATE( KppFixId      )
    IF ( ALLOCATED( KppVarId      ) ) DEALLOCATE( KppVarId      )
    IF ( ALLOCATED( KppSpcId      ) ) DEALLOCATE( KppSpcId      )

  END SUBROUTINE Cleanup_Work_Arrays
!EOC
END MODULE Species_Database_Mod
