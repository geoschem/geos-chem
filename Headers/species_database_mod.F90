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
  ! species from input.geos with the KPP species names (and removes duplicates)
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
  LOGICAL                     :: prtDebug
  LOGICAL                     :: v_bool
  INTEGER                     :: v_int
  INTEGER                     :: nSpecies
  INTEGER                     :: N
  INTEGER                     :: S
  REAL(f4)                    :: v_real
  REAL(f4)                    :: mw_g

  ! Strings
  CHARACTER(LEN=14)           :: tag
  CHARACTER(LEN=14)           :: spc
  CHARACTER(LEN=255)          :: v_str
  CHARACTER(LEN=255)          :: key
  CHARACTER(LEN=255)          :: fileDir
  CHARACTER(LEN=255)          :: fileName
  CHARACTER(LEN=255)          :: thisLoc
  CHARACTER(LEN=512)          :: errMsg

  ! Arrays
  REAL(f4)                    :: a_real_2(2)
  REAL(f4)                    :: a_real_3(3)

  ! String arrays
  CHARACTER(LEN=17)           :: tags(45)

  ! Objects
  TYPE(QFYAML_t)              :: yml
  TYPE(QFYAML_t)              :: yml_anchored
  TYPE(Species),    POINTER   :: ThisSpc

  !=========================================================================
  ! Init_Species_Database begins here!
  !=========================================================================

  ! Initialize
  RC         = GC_SUCCESS
  mw_g       = MISSING_R4
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
  SpcCount%nPhotol  = 0
  SpcCount%nWetDep  = 0
  SpcCount%nHg0     = 0
  SpcCount%nHg2     = 0
  SpcCount%nHgP     = 0  

  ! Species database tags to match
  tags = (/ "BackgroundVV     ", "DD_AeroDryDep    ", "DD_DustDryDep    ",   &
            "DD_DvzAerSnow    ", "DD_DvzMinVal     ", "DD_F0            ",   &
            "DD_Hstar         ", "DD_KOA           ", "Density          ",   &
            "Formula          ", "Fullname         ", "Is_ActiveChem    ",   &
            "Is_Advected      ", "Is_Aero          ", "Is_DryAlt        ",   &
            "Is_DryDep        ", "Is_FixedChem     ", "Is_HygroGrowth   ",   &
            "Is_Kpp           ", "Is_Gas           ", "Is_Hg0           ",   &
            "Is_Hg2           ", "Is_HgP           ", "Is_Photolysis    ",   &
            "Is_WetDep        ", "Henry_CR         ", "Henry_K0         ",   &
            "Henry_pKa        ", "MP_SizeResAer    ", "MP_SizeResNum    ",   &
            "MolecRatio       ", "MW_g             ", "EmMw_g           ",   &
            "Radius           ", "WD_AerScavEff    ", "WD_CoarseAer     ",   &
            "WD_ConvFacI2G    ", "WD_KcScaleFac    ", "WD_KcScaleFac_Luo",   &
            "WD_Is_H2SO4      ", "WD_Is_HNO3       ", "WD_Is_SO2        ",   &
            "WD_LiqAndGas     ", "WD_RainoutEff    ", "WD_RainoutEff_Luo"  /)

  !=========================================================================
  ! Store the list unique GEOS-Chem species names in work arrays for use
  ! below.  This is the combined list of advected species (from input.geos)
  ! plus KPP species (from SPC_NAMES in gckpp_Monitor.F90), with all
  ! duplicates removed.  Also stores the corresponding indices in the
  ! KPP VAR and FIX arrays.  For simulations that do not use KPP, the
  ! unique species list is the list of advected species from input.geos.
  !=========================================================================
  CALL Unique_Species_Names( Input_Opt, nSpecies, RC )

  ! Initialize the species vector
  CALL SpcData_Init( Input_Opt, nSpecies, SpcData, RC )
  IF ( RC /= GC_SUCCESS ) THEN
     errMsg = "Could not initialize the species database object!"
     CALL GC_Error( errMsg, RC, thisLoc )
     RETURN
  ENDIF

  print*, '### size spcdata:', nSpecies

  !=========================================================================
  ! Read the species database from a YAML file
  ! Store variable information into the yml* configuration objects
  !=========================================================================

  !-------------------------------------------------------------------------
  ! Set the directory for the species database file
  ! USE THIS PATH FOR GEOS-Chem 12.9.0
  fileDir = "./CodeDir/Headers/"
  !-------------------------------------------------------------------------
  ! Set the directory for the species database file
  ! USE THIS PATH FOR GEOS-Chem 13.0.0 AND LATER
  !fileDir = "./CodeDir/src/GEOS-Chem/Headers/"
  !-------------------------------------------------------------------------

  ! Set the path to the species database file
#if defined(TOMAS)
  fileName = TRIM( fileDir ) // "species_database_tomas.yml"
#elif defined(APM)
  fileName = TRIM( fileDir ) // "species_database_apm.yml"
#else
  fileName = TRIM( fileDir ) // "species_database.yml"
#endif

  ! Read species metadata into configuration objects
  CALL QFYAML_Init( fileName, yml, yml_anchored, RC )
  IF ( RC /= GC_SUCCESS ) THEN
     errMsg = "Error reading " // TRIM( fileName )
     CALL GC_Error( errMsg, RC, thisLoc )
     RETURN
  ENDIF
  
  !=========================================================================
  ! Extract species metadata into the species database object
  !=========================================================================

  ! Loop over the number of species
  DO S = 1, nSpecies

     ! Species name
     spc = species_names(S)

     ! Point to the corresponding species database entry
     ThisSpc => SpcData(S)%Info

     ! Model Id = overall species index
     ThisSpc%ModelId = S

     ! Short species name
     ThisSpc%Name = TRIM( spc )

     ! Loop over the number of tags in the species database
     DO N = 1, SIZE( tags )

        ! Set intial values to default "missing" values
        ! If the tag isn't found for a given species, then
        ! it will be given the appropriate missing value.
        a_real_2 = MISSING_R4
        a_real_3 = MISSING_R4
        v_bool   = MISSING_BOOL
        v_int    = MISSING_INT
        v_real   = MISSING_R4
        v_str    = MISSING_STR

        ! Create search key for each variable
        key = TRIM( spc ) // '%' // TRIM( tags(N) )

        ! Save into the proper field of the species database
        IF ( INDEX( key, "%BackgroundVV" ) > 0 ) THEN
           CALL QFYAML_Add_Get( yml, key, v_real, "", RC )
           IF ( RC /= GC_SUCCESS ) RETURN
           IF ( prtDebug ) WRITE( 6, 20 ) TRIM( key ), v_real
           ThisSpc%BackgroundVV = v_real

        ELSE IF ( INDEX( key, "%DD_AeroDryDep" ) > 0 ) THEN
           CALL QFYAML_Add_Get( yml, key, v_bool, "", RC )
           IF ( RC /= GC_SUCCESS ) RETURN
           IF ( prtDebug ) WRITE( 6, 20 ) TRIM( key ), v_bool
           ThisSpc%DD_AeroDryDep = v_bool

        ELSE IF ( INDEX( key, "%DD_DustDryDep" ) > 0 ) THEN
           CALL QFYAML_Add_Get( yml, key, v_bool, "", RC )
           IF ( RC /= GC_SUCCESS ) RETURN
           IF ( prtDebug ) WRITE( 6, 20 ) TRIM( key ), v_bool
           ThisSpc%DD_DustDryDep = v_bool

        ELSE IF ( INDEX( key, "%DD_DvzAerSnow" ) > 0 ) THEN
           CALL QFYAML_Add_Get( yml, key, v_real, "", RC )
           IF ( RC /= GC_SUCCESS ) RETURN
           IF ( prtDebug ) WRITE( 6, 30 ) TRIM( key ), v_real
           ThisSpc%DD_DvzAerSnow = v_real

        ELSE IF ( INDEX( key, "%DD_DvzMinVal" ) > 0 ) THEN
           CALL QFYAML_Add_Get( yml, key, a_real_2, "", RC )
           IF ( RC /= GC_SUCCESS ) RETURN
           IF ( prtDebug ) WRITE( 6, 31 ) TRIM( key ), a_real_2
           ThisSpc%DD_DvzMinVal = a_real_2

        ELSE IF ( INDEX( key, "%DD_F0" ) > 0 ) THEN
           CALL QFYAML_Add_Get( yml, key, v_real, "", RC )
           IF ( RC /= GC_SUCCESS ) RETURN
           IF ( prtDebug ) WRITE( 6, 30 ) TRIM( key ), v_real
           ThisSpc%DD_F0 = v_real

        ELSE IF ( INDEX( key, "%DD_Hstar" ) > 0 ) THEN
           CALL QFYAML_Add_Get( yml, key, v_real, "", RC )
           IF ( RC /= GC_SUCCESS ) RETURN
           IF ( prtDebug ) WRITE( 6, 30 ) TRIM( key ), v_real
           ThisSpc%DD_Hstar = v_real

        ELSE IF ( INDEX( key, "%Density" ) > 0 ) THEN
           CALL QFYAML_Add_Get( yml, key, v_real, "", RC )
           IF ( RC /= GC_SUCCESS ) RETURN
           IF ( prtDebug ) WRITE( 6, 30 ) TRIM( key ), v_real
           ThisSpc%Density = v_real

        ELSE IF ( INDEX( key, "%Formula" ) > 0 ) THEN
           CALL QFYAML_Add_Get( yml, key, v_str, "", RC )
           IF ( RC /= GC_SUCCESS ) RETURN
           IF ( prtDebug ) WRITE( 6, 10 ) TRIM( key ), TRIM( v_str )
           ThisSpc%Formula = TRIM( v_str )

        ELSE IF ( INDEX( key, "%Fullname" ) > 0 ) THEN
           CALL QFYAML_Add_Get( yml, key, v_str, "", RC )
           IF ( RC /= GC_SUCCESS ) RETURN
           IF ( prtDebug ) WRITE( 6, 10 ) TRIM( key ), TRIM( v_str )
           ThisSpc%Fullname = TRIM( v_str )

        ELSE IF ( INDEX( key, "%Is_Advected" ) > 0 ) THEN
           CALL QFYAML_Add_Get( yml, key, v_bool, "", RC )
           IF ( RC /= GC_SUCCESS ) RETURN
           IF ( prtDebug ) WRITE( 6, 20 ) TRIM( key ), v_bool
           ThisSpc%Is_Advected = v_bool
           SpcCount%nAdvect    = SpcCount%nAdvect + 1

        ELSE IF ( INDEX( key, "%Is_Aero" ) > 0  ) THEN
           CALL QFYAML_Add_Get( yml, key, v_bool, "", RC )
           IF ( RC /= GC_SUCCESS ) RETURN
           IF ( prtDebug ) WRITE( 6, 20 ) TRIM( key ), v_bool
           SpcCount%nAeroSpc = SpcCount%nAeroSpc + 1
           ThisSpc%AeroId    = SpcCount%nAeroSpc
           ThisSpc%Is_Aero   = v_bool

        ELSE IF ( INDEX( key, "%Is_DryAlt" ) > 0 ) THEN
           CALL QFYAML_Add_Get( yml, key, v_bool, "", RC )
           IF ( RC /= GC_SUCCESS ) RETURN
           IF ( prtDebug ) WRITE( 6, 20 ) TRIM( key ), v_bool
           SpcCount%nDryAlt  = SpcCount%nDryAlt + 1
           ThisSpc%DryAltId  = SpcCount%nDryAlt
           ThisSpc%Is_DryAlt = v_bool

        ELSE IF ( INDEX( key, "%Is_DryDep" ) > 0 ) THEN
           CALL QFYAML_Add_Get( yml, key, v_bool, "", RC )
           IF ( RC /= GC_SUCCESS ) RETURN
           IF ( prtDebug ) WRITE( 6, 20 ) TRIM( key ), v_bool
           SpcCount%nDryDep  = SpcCount%nDryDep + 1
           ThisSpc%DryDepId  = SpcCount%nDryDep
           ThisSpc%Is_DryDep = v_bool

        ELSE IF ( INDEX( key, "%Is_HygroGrowth" ) > 0 ) THEN
           CALL QFYAML_Add_Get( yml, key, v_bool, "", RC )
           IF ( RC /= GC_SUCCESS ) RETURN
           IF ( prtDebug ) WRITE( 6, 20 ) TRIM( key ), v_bool
           SpcCount%nHygGrth      = SpcCount%nHygGrth + 1
           ThisSpc%HygGrthId      = SpcCount%nHygGrth
           ThisSpc%Is_HygroGrowth = v_bool

        ELSE IF ( INDEX( key, "%Is_Gas" ) > 0 ) THEN
           CALL QFYAML_Add_Get( yml, key, v_bool, "", RC )
           IF ( RC /= GC_SUCCESS ) RETURN
           IF ( prtDebug ) WRITE( 6, 20 ) TRIM( key ), v_bool
           SpcCount%nGasSpc = SpcCount%nGasSpc + 1
           ThisSpc%GasSpcId = SpcCount%nGasSpc
           ThisSpc%Is_Gas   = v_bool

        ELSE IF ( INDEX( key, "%Is_Hg0" ) > 0 ) THEN
           CALL QFYAML_Add_Get( yml, key, v_bool, "", RC )
           IF ( RC /= GC_SUCCESS ) RETURN
           IF ( prtDebug ) WRITE( 6, 20 ) TRIM( key ), v_bool
           ThisSpc%Is_Hg0 = v_bool
           SpcCount%nHg0  = SpcCount%nHg0 + 1

        ELSE IF ( INDEX( key, "%Is_Hg2" ) > 0 ) THEN
           CALL QFYAML_Add_Get( yml, key, v_bool, "", RC )
           IF ( RC /= GC_SUCCESS ) RETURN
           IF ( prtDebug ) WRITE( 6, 20 ) TRIM( key ), v_bool
           ThisSpc%Is_Hg2 = v_bool
           SpcCount%nHg2  = SpcCount%nHg2 + 1

        ELSE IF ( INDEX( key, "%Is_HgP" ) > 0 ) THEN
           CALL QFYAML_Add_Get( yml, key, v_bool, "", RC )
           IF ( RC /= GC_SUCCESS ) RETURN
           IF ( prtDebug ) WRITE( 6, 20 ) TRIM( key ), v_bool
           ThisSpc%Is_HgP = v_bool
           SpcCount%nHgP  = SpcCount%nHgP + 1

        ELSE IF ( INDEX( key, "%Is_Photolysis" ) > 0 ) THEN
           CALL QFYAML_Add_Get( yml, key, v_bool, "", RC )
           IF ( RC /= GC_SUCCESS ) RETURN
           IF ( prtDebug ) WRITE( 6, 20 ) TRIM( key ), v_bool
           SpcCount%nPhotol      = SpcCount%nPhotol + 1
           ThisSpc%PhotolId      = SpcCount%nPhotol
           ThisSpc%Is_Photolysis = v_bool

        ELSE IF ( INDEX( key, "%Is_WetDep" ) > 0 ) THEN
           CALL QFYAML_Add_Get( yml, key, v_bool, "", RC )
           IF ( RC /= GC_SUCCESS ) RETURN
           IF ( prtDebug ) WRITE( 6, 20 ) TRIM( key ), v_bool
           SpcCount%nWetDep  = SpcCount%nWetDep + 1
           ThisSpc%WetDepID  = SpcCount%nWetDep
           ThisSpc%Is_WetDep = v_bool

        ELSE IF ( INDEX( key, "%Henry_CR" ) > 0 ) THEN
           CALL QFYAML_Add_Get( yml, key, v_real, "", RC )
           IF ( RC /= GC_SUCCESS ) RETURN
           IF ( prtDebug ) WRITE( 6, 30 ) TRIM( key ), v_real
           ThisSpc%Henry_CR = v_real

        ELSE IF ( INDEX( key, "%Henry_K0" ) > 0 ) THEN
           CALL QFYAML_Add_Get( yml, key, v_real, "", RC )
           IF ( RC /= GC_SUCCESS ) RETURN
           IF ( prtDebug ) WRITE( 6, 30 ) TRIM( key ), v_real
           ThisSpc%Henry_K0 = v_real

        ELSE IF ( INDEX( key, "%Henry_pKa" ) > 0 ) THEN
           CALL QFYAML_Add_Get( yml, key, v_real, "", RC )
           IF ( RC /= GC_SUCCESS ) RETURN
           IF ( prtDebug ) WRITE( 6, 30 ) TRIM( key ), v_real
           ThisSpc%Henry_pKa = v_real

        ELSE IF ( INDEX( key, "%MP_SizeResAer" ) > 0 ) THEN
           CALL QFYAML_Add_Get( yml, key, v_bool, "", RC )
           IF ( RC /= GC_SUCCESS ) RETURN
           IF ( prtDebug ) WRITE( 6, 20 ) TRIM( key ), v_bool
           ThisSpc%MP_SizeResAer = v_bool

        ELSE IF ( INDEX( key, "%MP_SizeResNum" ) > 0 ) THEN
           CALL QFYAML_Add_Get( yml, key, v_bool, "", RC )
           IF ( RC /= GC_SUCCESS ) RETURN
           IF ( prtDebug ) WRITE( 6, 20 ) TRIM( key ), v_bool
           ThisSpc%MP_SizeResNum = v_bool

        ELSE IF ( INDEX( key, "%MolecRatio" ) > 0 ) THEN
           v_real = ONE_R4                               ! Set default to 1
           CALL QFYAML_Add_Get( yml, key, v_real, "", RC )
           IF ( RC /= GC_SUCCESS ) RETURN
           IF ( prtDebug ) WRITE( 6, 30 ) TRIM( key ), v_real
           ThisSpc%MolecRatio = v_real

        ELSE IF ( INDEX( key, "%MW_g" ) > 0 ) THEN
           CALL QFYAML_Add_Get( yml, key, v_real, "", RC )
           IF ( RC /= GC_SUCCESS ) RETURN
           IF ( prtDebug ) WRITE( 6, 30 ) TRIM( key ), v_real
           ThisSpc%MW_g = v_real
           mw_g         = v_real                         ! For missing EmMw_g

        ELSE IF ( INDEX( key, "%EmMW_g" ) > 0 ) THEN
           v_real = mw_g                            
           CALL QFYAML_Add_Get( yml, key, v_real, "", RC )
           IF ( RC /= GC_SUCCESS ) RETURN
           IF ( prtDebug ) WRITE( 6, 30 ) TRIM( key ), v_real
           ThisSpc%EmMw_G = v_real
           mw_g           = MISSING_R4                   ! Reset for next spc

        ELSE IF ( INDEX( key, "%Radius" ) > 0 ) THEN
           CALL QFYAML_Add_Get( yml, key, v_real, "", RC )
           IF ( RC /= GC_SUCCESS ) RETURN
           IF ( prtDebug ) WRITE( 6, 30 ) TRIM( key ), v_real
           ThisSpc%Radius = v_real

        ELSE IF ( INDEX( key, "%WD_AerScavEff" ) > 0 ) THEN
           CALL QFYAML_Add_Get( yml, key, v_real, "", RC )
           IF ( RC /= GC_SUCCESS ) RETURN
           IF ( prtDebug ) WRITE( 6, 30 ) TRIM( key ), v_real
           ThisSpc%WD_AerScavEff = v_real

        ELSE IF ( INDEX( key, "%WD_CoarseAer" ) > 0 ) THEN
           CALL QFYAML_Add_Get( yml, key, v_bool, "", RC )
           IF ( RC /= GC_SUCCESS ) RETURN
           IF ( prtDebug ) WRITE( 6, 20 ) TRIM( key ), v_bool
           ThisSpc%WD_CoarseAer = v_bool

        ELSE IF ( INDEX( key, "%WD_ConvFacI2G" ) > 0 ) THEN
           CALL QFYAML_Add_Get( yml, key, v_real, "", RC )
           IF ( RC /= GC_SUCCESS ) RETURN
           IF ( prtDebug ) WRITE( 6, 30 ) TRIM( key ), v_real
           ThisSpc%WD_ConvFacI2G = v_real

#ifdef LUO_WETDEP
        ELSE IF ( INDEX( key, "%WD_KcScaleFac_Luo" ) > 0 ) THEN
           CALL QFYAML_Add_Get( yml, key, a_real_3, "", RC )
           IF ( RC /= GC_SUCCESS ) RETURN
           IF ( prtDebug ) WRITE( 6, 32 ) TRIM( key ), a_real_3
           ThisSpc%WD_KcScaleFac = a_real_3
#else
        ELSE IF ( INDEX( key, "%WD_KcScaleFac" ) > 0 ) THEN
           CALL QFYAML_Add_Get( yml, key, a_real_3, "", RC )
           IF ( RC /= GC_SUCCESS ) RETURN
           IF ( prtDebug ) WRITE( 6, 32 ) TRIM( key ), a_real_3
           ThisSpc%WD_KcScaleFac = a_real_3
#endif

        ELSE IF ( INDEX( key, "%WD_Is_H2SO4" ) > 0 ) THEN
           CALL QFYAML_Add_Get( yml, key, v_bool, "", RC )
           IF ( RC /= GC_SUCCESS ) RETURN
           IF ( prtDebug ) WRITE( 6, 20 ) TRIM( key ), v_bool
           ThisSpc%WD_Is_H2SO4 = v_bool

        ELSE IF ( INDEX( key, "%WD_Is_HNO3" ) > 0 ) THEN
           CALL QFYAML_Add_Get( yml, key, v_bool, "", RC )
           IF ( RC /= GC_SUCCESS ) RETURN
           IF ( prtDebug ) WRITE( 6, 20 ) TRIM( key ), v_bool
           ThisSpc%WD_Is_HNO3 = v_bool

        ELSE IF ( INDEX( key, "%WD_Is_SO2" ) > 0 ) THEN
           CALL QFYAML_Add_Get( yml, key, v_bool, "", RC )
           IF ( RC /= GC_SUCCESS ) RETURN
           IF ( prtDebug ) WRITE( 6, 20 ) TRIM( key ), v_bool
           ThisSpc%WD_Is_SO2 = v_bool

        ELSE IF ( INDEX( key, "%WD_LiqAndGas" ) > 0 ) THEN
           CALL QFYAML_Add_Get( yml, key, v_bool, "", RC )
           IF ( RC /= GC_SUCCESS ) RETURN
           IF ( prtDebug ) WRITE( 6, 20 ) TRIM( key ), v_bool
           ThisSpc%WD_LiqAndGas = v_bool

#ifdef LUO WETDEP
        ELSE IF ( INDEX( key, "%WD_RainoutEff_Luo" ) > 0 ) THEN
           CALL QFYAML_Add_Get( yml, key, a_real_3, "", RC )
           IF ( RC /= GC_SUCCESS ) RETURN
           IF ( prtDebug ) WRITE( 6, 32 ) TRIM( key ), a_real_3
           ThisSpc%WD_RainoutEff = a_real_3
#else
        ELSE IF ( INDEX( key, "%WD_RainoutEff" ) > 0 ) THEN
           CALL QFYAML_Add_Get( yml, key, a_real_3, "", RC )
           IF ( RC /= GC_SUCCESS ) RETURN
           IF ( prtDebug ) WRITE( 6, 32 ) TRIM( key ), a_real_3
           ThisSpc%WD_RainoutEff = a_real_3
#endif

        ELSE
           ! Pass

        ENDIF

     ENDDO

     ! Free pointer
     ThisSpc => NULL()

     ! Debug output
     IF ( prtDebug ) WRITE( 6, '(a)' )

     !---------------------------------------------------------------------
     ! Manually update some fields in anchored variables
     !---------------------------------------------------------------------
     SELECT CASE( TRIM( spc ) ) 
        CASE( 'Be7' ) 
           key   = "Be7%Fullname"
           v_str = "Beryllium-7 isotope"
           CALL QFYAML_Update( yml, key, v_str )
           IF ( prtDebug ) WRITE( 6, 10 ) TRIM( key ), TRIM( v_str )
           ThisSpc%Fullname = TRIM( v_str )

           key   = "Be7%Formula"
           v_str = "Be"
           CALL QFYAML_Update( yml, key, v_str )
           IF ( prtDebug ) WRITE( 6, 10 ) TRIM( key ), TRIM( v_str )
           ThisSpc%Fullname = TRIM( v_str )

        CASE( 'Be7Strat' ) 
           key   = "Be7%Fullname"
           v_str = "Beryllium-7 isotope in stratosphere"
           CALL QFYAML_Update( yml, key, v_str )
           IF ( prtDebug ) WRITE( 6, 10 ) TRIM( key ), TRIM( v_str )
           ThisSpc%Fullname = TRIM( v_str )

           key   = "Be7%Formula"
           v_str = "Be"
           CALL QFYAML_Update( yml, key, v_str )
           IF ( prtDebug ) WRITE( 6, 10 ) TRIM( key ), TRIM( v_str )
           ThisSpc%Fullname = TRIM( v_str )

        CASE( 'Be10' ) 
           key   = "Be7%Fullname"
           v_str = "Beryllium-10 isotope"
           CALL QFYAML_Update( yml, key, v_str )
           IF ( prtDebug ) WRITE( 6, 10 ) TRIM( key ), TRIM( v_str )
           ThisSpc%Fullname = TRIM( v_str )

           key   = "Be10%Formula"
           v_str = "Be10"
           CALL QFYAML_Update( yml, key, v_str )
           IF ( prtDebug ) WRITE( 6, 10 ) TRIM( key ), TRIM( v_str )
           ThisSpc%Fullname = TRIM( v_str )

        CASE( 'Be10Strat' ) 
           key   = "Be7%Fullname"
           v_str = "Beryllium-7 isotope in stratosphere"
           CALL QFYAML_Update( yml, key, v_str )
           IF ( prtDebug ) WRITE( 6, 10 ) TRIM( key ), TRIM( v_str )
           ThisSpc%Fullname = TRIM( v_str )

           key   = "Be10%Formula"
           v_str = "Be10"
           CALL QFYAML_Update( yml, key, v_str )
           IF ( prtDebug ) WRITE( 6, 10 ) TRIM( key ), TRIM( v_str )
           ThisSpc%Fullname = TRIM( v_str )

        CASE DEFAULT
           ! Pass
     END SELECT

  ENDDO

  ! FORMAT statements
10 FORMAT( a30, " | ", a      )
20 FORMAT( a30, " | ", L10    )
30 FORMAT( a30, " | ", f10.2  )
31 FORMAT( a30, " | ", 2f10.2 )
32 FORMAT( a30, " | ", 3f10.2 )
40 FORMAT( a30, " | ", i10    )

  !=========================================================================
  ! Finalize the config objects
  !=========================================================================
  CALL QFYAML_CleanUp( yml          )
  CALL QFYAML_CleanUp( yml_anchored )
  stop

  END SUBROUTINE Init_Species_Database
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
    USE Species_Mod,      ONLY : MISSING_INT
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

    ! Arrays
    CHARACTER(LEN=15), ALLOCATABLE :: Tmp(:)
    CHARACTER(LEN=15)              :: SpcName

    !=======================================================================
    ! UNIQUE_SPECIES_NAMES begins here!
    !=======================================================================

    ! Assume success
    RC       = GC_SUCCESS

    ! Number of advected species listed in input.geos
    nAdvect  = Input_Opt%N_Advect

    ! First set the # of species to the # of advected species
    nSpecies = nAdvect

    !=======================================================================
    ! For full-chemistry simulations with KPP, get the list of all of
    ! species names in the KPP mechanism, and their indices
    !=======================================================================
    IF ( Input_Opt%ITS_A_FULLCHEM_SIM ) THEN

       ! Allocate a temporary array large enough to hold all of the
       ! advected species listed in input.geos as well as all of the
       ! KPP species names (listed in SPC_NAMES of gckpp_Monitor.F90)
       ALLOCATE( Tmp( nAdvect + NSPEC ), STAT=RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       Tmp = ''

       !--------------------------------------------------------------------
       ! First determine the unique list of species in the KPP mechanism
       ! (so that we don't duplicate storage for advected & chemical species)
       !--------------------------------------------------------------------

       ! First, store advected species (from input.geos) in the TMP array
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
       ALLOCATE( Species_Names( nSpecies ) )
       Species_Names = Tmp(1:nSpecies )

       ! Free temporary array
       IF ( ALLOCATED( Tmp ) ) DEALLOCATE( Tmp )

       !--------------------------------------------------------------------
       ! Now determine the KPP indices for each unique species name
       !--------------------------------------------------------------------

       ! Work array to hold the list of all KPP species indices
       ALLOCATE( KppSpcId( nSpecies ), STAT=RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       KppSpcId = MISSING_INT

       ! Work array to hold the list of KPP fixed species indices
       ALLOCATE( KppFixId( nSpecies ), STAT=RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       KppFixId = MISSING_INT

       ! Work array to hold the list of KPP variable species indices
       ALLOCATE( KppVarId( nSpecies ), STAT=RC )
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
    ! of species is just the list of advected species from input.geos
    !=======================================================================
    ELSE

       ! Initialize the species names array from Input_Opt
       ALLOCATE( Species_Names( nSpecies ), STAT=RC )
       Species_Names = Input_Opt%AdvectSpc_Name(1:nSpecies)

       ! Set KppSpcId to missing value
       ALLOCATE( KppSpcId( nSpecies ), STAT=RC )
       KppSpcId = MISSING_INT

       ! Set KppFixId to missing value
       ALLOCATE( KppFixId( nSpecies ), STAT=RC )
       KppFixId = MISSING_INT

       ! Set KppVarId to missing value
       ALLOCATE( KppVarId( nSpecies ), STAT=RC )
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
