!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: species_mod.F90
!
! !DESCRIPTION: Module SPECIES\_MOD contains types and routines to define
!  the GEOS-Chem species object.
!\\
!\\
! !INTERFACE:
!
MODULE Species_Mod
!
! USES:
!
  USE Precision_Mod

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC :: SpcData_Init
  PUBLIC :: SpcData_Cleanup
  PUBLIC :: Spc_Create
  PUBLIC :: Spc_GetNumSpecies
  PUBLIC :: Spc_Print
!
! !PUBLIC TYPES:
!
  !=========================================================================
  ! Counters for the species indices
  !=========================================================================
  INTEGER, PRIVATE :: AdvectCount  = 0  ! Counter of advected species
  INTEGER, PRIVATE :: AeroCount    = 0  ! Counter of aerosol species
  INTEGER, PRIVATE :: DryAltCount  = 0  ! Counter of dry-dep species to save
                                        !  at a user-defined altitude
  INTEGER, PRIVATE :: DryDepCount  = 0  ! Counter of dry-deposited species
  INTEGER, PRIVATE :: GasSpcCount  = 0  ! Counter of gas-phase species
  INTEGER, PRIVATE :: HygGrthCount = 0  ! Counter of hygroscopic growth spc
  INTEGER, PRIVATE :: KppVarCount  = 0  ! Counter of variable KPP species
  INTEGER, PRIVATE :: KppFixCount  = 0  ! Counter of fixed KPP species
  INTEGER, PRIVATE :: KppSpcCount  = 0  ! Counter of species in KPP matrix
  INTEGER, PRIVATE :: PhotolCount  = 0  ! Counter of photolysis species
  INTEGER, PRIVATE :: WetDepCount  = 0  ! Counter of wet-deposited species
  INTEGER, PRIVATE :: Hg0Count     = 0  ! Number  of Hg0 tracers
  INTEGER, PRIVATE :: Hg2Count     = 0  ! Number  of Hg2 tracers
  INTEGER, PRIVATE :: HgPCount     = 0  ! Number  of HgP tracers

  !=========================================================================
  ! Type for the Species Database object (vector of type Species)
  !=========================================================================
  TYPE, PUBLIC :: SpcPtr
     TYPE(Species), POINTER :: Info         ! Single entry of Species Database
  END TYPE SpcPtr

  !=========================================================================
  ! Type for individual species information
  ! (i.e. this is a single entry in the Species Database)
  !=========================================================================
  TYPE, PUBLIC :: Species

     ! Indices
     INTEGER            :: ModelID          ! Model species ID
     INTEGER            :: AdvectID         ! Advection index
     INTEGER            :: AeroID           ! Aerosol species index
     INTEGER            :: DryAltId         ! Dry dep species at altitude ID
     INTEGER            :: DryDepID         ! Dry deposition index
     INTEGER            :: GasSpcID         ! Gas-phase species index
     INTEGER            :: HygGrthID        ! Hygroscopic growth species index
     INTEGER            :: KppVarId         ! KPP variable species index
     INTEGER            :: KppFixId         ! KPP fixed spcecies index
     INTEGER            :: KppSpcID         ! KPP species index
     INTEGER            :: PhotolID         ! Photolysis index
     INTEGER            :: WetDepID         ! Wet deposition index

     ! Names
     CHARACTER(LEN=31)  :: Name             ! Short name
     CHARACTER(LEN=80)  :: FullName         ! Long name
     CHARACTER(LEN=80)  :: Formula          ! Chemical formula

     ! Logical switches
     LOGICAL            :: Is_Advected      ! Is it advected?
     LOGICAL            :: Is_Aero          ! Is it an aerosol species?
     LOGICAL            :: Is_DryAlt        ! Is it a dry-dep species that we
                                            !  want to save at a given altitude?
     LOGICAL            :: Is_DryDep        ! Is it dry-deposited?
     LOGICAL            :: Is_Gas           ! Is it a gas?  If not, aerosol.
     LOGICAL            :: Is_HygroGrowth   ! Does it have hygroscropic growth?
     LOGICAL            :: Is_ActiveChem    ! Is it an active chemical species?
     LOGICAL            :: Is_FixedChem     ! Is it a fixed chemical species?
     LOGICAL            :: Is_Kpp           ! Is it in the KPP mechanism?
     LOGICAL            :: Is_Photolysis    ! Is it an photolysis species?
     LOGICAL            :: Is_WetDep        ! Is it wet-deposited?
     LOGICAL            :: Is_InRestart     ! Is it in the restart file?

     ! Molecular weights
     REAL(fp)           :: MW_g             ! Species molecular weight [g/mol]
     REAL(fp)           :: EmMW_g           ! Emitted molecular weight [g/mol]
     REAL(fp)           :: MolecRatio       ! Mol carbon / mol species [1    ]

     ! Default background concentration
     REAL(fp)           :: BackgroundVV     ! Background conc [v/v]

     ! Density and radius
     REAL(fp)           :: Density          ! Density [kg/m3]
     REAL(fp)           :: Radius           ! Radius  [m]

     ! Henry's law parameters
     REAL(f8)           :: Henry_K0         ! Liq./gas Henry const [M/atm ]
     REAL(f8)           :: Henry_CR         ! d(ln K0) / d(1/T)    [K     ]
     REAL(f8)           :: Henry_PKA        ! pKa for Henry const. correction

     ! Drydep parameters
     LOGICAL            :: DD_AeroDryDep    ! Use AERO_SFCRSII for drydep?
     LOGICAL            :: DD_DustDryDep    ! Use DUST_SFCRSII for drydep?
     REAL(fp)           :: DD_DvzAerSnow    ! Vd for aerosols on snow [cm/s]
     REAL(fp)           :: DD_DvzMinVal(2)  ! Min Vd for aerosols [cm/s]
     REAL(fp)           :: DD_F0            ! F0 (reactivity) factor [1]
     REAL(fp)           :: DD_KOA           ! KOA factor for POPG
     !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     !%%% NOTE: We will eventually replace this with the common Henry's law
     !%%% parameters.  But in order to replicate the prior behavior,
     !%%% we will need to supply the dry deposition code with the same
     !%%% HSTAR values that are currently set in INIT_DRYDEP.  Therefore,
     !%%% add this field as a temporary placeholder for the Hstar quantity
     !%%% from drydep_mod.F90.  We will remove this later on. (bmy, 8/24/15)
     !%%%
     REAL(fp)           :: DD_Hstar_Old     ! HSTAR value in drydep_mod [M/atm]
     !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

     ! Wetdep parameters, gas-phase species
     LOGICAL            :: WD_LiqAndGas     ! Consider liquid and gas phases?
     REAL(fp)           :: WD_ConvFacI2G    ! Conv. factor for ice/gas ratio
     REAL(fp)           :: WD_RetFactor     ! Retention factor [1]

     ! Wetdep parameters, aerosol-phase species
     LOGICAL            :: WD_Is_H2SO4      ! Flag to denote H2SO4 wetdep
     LOGICAL            :: WD_Is_HNO3       ! Flag to denote HNO3 wetdep
     LOGICAL            :: WD_Is_SO2        ! Flag to denote SO2 wetdep
     LOGICAL            :: WD_CoarseAer     ! T=coarse aerosol; F=fine aerosol
     REAL(fp)           :: WD_AerScavEff    ! Aerosol scavenging efficiency
     REAL(fp)           :: WD_KcScaleFac(3) ! Temperature-dependent scale
                                            !  factors to multiply Kc rate
                                            !  (conv of condensate -> precip)
                                            !  in F_AEROSOL (wetscav_mod.F90)
     REAL(fp)           :: WD_RainoutEff(3) ! Temperature-dependent scale
                                            !  factors for rainout efficiency

     ! Microphysics parameters
     LOGICAL            :: MP_SizeResAer    ! T=size-resolved aerosol (TOMAS)
     LOGICAL            :: MP_SizeResNum    ! T=size-resolved aerosol number

     ! Tagged mercury parameters
     LOGICAL            :: Is_Hg0           ! T=total or tagged Hg0 species
     LOGICAL            :: Is_Hg2           ! T=total or tagged Hg2 species
     LOGICAL            :: Is_HgP           ! T=total or tagged HgP species
     INTEGER            :: Hg_Cat           ! Tagged Hg category number

  END TYPE Species
!
! !DEFINED PARAMETERS
!
  !=========================================================================
  ! Missing value parameters
  !=========================================================================
  INTEGER,  PARAMETER, PUBLIC :: MISSING_INT = -999         ! Integer
  REAL(fp), PARAMETER, PUBLIC :: MISSING     = -999e+0_fp   ! Flexible precision
  REAL(f8), PARAMETER, PUBLIC :: MISSING_R8  = -999e+0_f8   ! 8-byte precision
  REAL(fp), PARAMETER, PUBLIC :: ZERO        =  0.0e+0_fp   ! Flexible precision
  REAL(f8), PARAMETER, PUBLIC :: ZERO_R8     =  0.0e+0_f8   ! 8-byte precision
  REAL(fp), PARAMETER, PUBLIC :: MISSING_MW  = 1.0_fp       ! Missing MW values

  !=========================================================================
  ! Missing species concentration value if not in restart file and special
  ! background value not defined
  !=========================================================================
  REAL(fp), PARAMETER, PUBLIC :: MISSING_VV  = 1.0e-20_fp ! Missing spc conc
!
! !REMARKS:
! (1) The emission molecular weight is the molecular weight of the emitted
!      compound. This value is only different to MW_g if the emitted compound
!      does not correspond to the transported species, e.g. if emissions are
!      in kg C4H10 but the corresponding species is transported as mass Carbon.
! (2) MolecRatio is the ratio between # of species molecules per emitted
!      molecule, e.g. 4 if emissions are kg C4H10 but model species are kg C.
!
! !REVISION HISTORY:
!  28 Feb 2014 - C. Keller   - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: SpcData_Init
!
! !DESCRIPTION: Routine SpcData\_Init initializes species database object.
!  This is an array where each element is of type Species.  This object holds
!  the metadata for each species (name, molecular weight, Henry's law
!  constants, drydep info, wetdep info, etc.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE SpcData_Init( Input_Opt, nSpecies, SpecDB, RC )
!
! !USES:
!
    USE Input_Opt_Mod, ONLY : OptInput
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput),       INTENT(IN)    :: Input_Opt    ! Input Options object
    INTEGER,              INTENT(IN)    :: nSpecies     ! # of species
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(SpcPtr),         POINTER       :: SpecDB(:)    ! Species database
    INTEGER,              INTENT(INOUT) :: RC           ! Return code
!
! !REVISION HISTORY:
!  20 Aug 2013 - C. Keller   - Adapted from gigc_state_chm_mod.F90
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER :: N, AS

    !=====================================================================
    ! SpcData_Init begins here!
    !=====================================================================

    ! Check if already allocated
    IF ( ASSOCIATED( SpecDB ) ) THEN
       CALL SpcData_Cleanup( SpecDB )
    ENDIF

    ! Allocate collection type
    ALLOCATE ( SpecDB( nSpecies ), STAT=AS )
    IF ( AS/=0 ) THEN
       RC = -999
       RETURN
    ENDIF

    ! Nullify each species pointer
    DO N = 1, nSpecies
       NULLIFY( SpecDB(N)%Info )
    ENDDO

  END SUBROUTINE SpcData_Init
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: SpcData_Cleanup
!
! !DESCRIPTION: Routine SpcData\_Cleanup cleans up the passed species
! collection object
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE SpcData_Cleanup( SpecDB )
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(SpcPtr), POINTER :: SpecDB(:)  ! Species database object
!
! !REVISION HISTORY:
!  20 Aug 2013 - C. Keller   - Adapted from gigc_state_chm_mod.F90
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER :: N, nSpecies

    !=====================================================================
    ! SpcData_Cleanup begins here!
    !=====================================================================

    ! Check if already allocated
    IF ( ASSOCIATED( SpecDB ) ) THEN

       ! First get the size of the SpecDb object
       nSpecies = SIZE( SpecDb )

       ! If there are more than 0 elements ...
       IF ( nSpecies > 0 ) THEN

          ! Nullify each entry in the species database
          DO N = 1, nSpecies
             IF( ASSOCIATED( SpecDB(N)%Info ) ) THEN
                DEALLOCATE( SpecDB(N)%Info )
             ENDIF
          ENDDO

          ! And free the object's memory
          DEALLOCATE( SpecDB )
       ENDIF
    ENDIF

  END SUBROUTINE SpcData_Cleanup
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Spc_Create
!
! !DESCRIPTION: Routine Spc\_Create creates a new object that holds
!  information about a given species, and assigns values to it.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Spc_Create( Input_Opt,      ThisSpc,       ModelID,             &
                         DryDepID,       Name,          FullName,            &
                         Formula,                                            &
                         MW_g,           EmMW_g,        MolecRatio,          &
                         BackgroundVV,   Henry_K0,      Henry_CR,            &
                         Henry_PKA,      Density,       Radius,              &
                         DD_AeroDryDep,  DD_DustDryDep, DD_DvzAerSnow,       &
                         DD_DvzMinVal,   DD_F0,         DD_KOA,              &
                         DD_HStar_Old,   MP_SizeResAer, MP_SizeResNum,       &
                         WD_RetFactor,   WD_LiqAndGas,  WD_ConvFacI2G,       &
                         WD_AerScavEff,  WD_KcScaleFac, WD_RainoutEff,       &
                         WD_CoarseAer,   Is_Advected,   Is_DryAlt,           &
                         Is_Drydep,      Is_Gas,        Is_HygroGrowth,      &
                         Is_Photolysis,  Is_Wetdep,     Is_InRestart,        &
                         Is_Hg0,         Is_Hg2,        Is_HgP,              &
                         KppSpcId,       KppVarId,      KppFixId,            &
                         RC                                                 )
!
! !USES:
!
    USE CharPak_Mod,        ONLY : To_UpperCase
    USE Input_Opt_Mod,      ONLY : OptInput
    USE PhysConstants,      ONLY : AIRMW,      AVO
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput),   INTENT(IN)  :: Input_Opt        ! Input Options object
    INTEGER,          OPTIONAL    :: ModelID          ! Model ID number
    INTEGER,          OPTIONAL    :: DryDepID         ! Drydep ID number
    CHARACTER(LEN=*), OPTIONAL    :: Name             ! Short name of species
    CHARACTER(LEN=*), OPTIONAL    :: FullName         ! Long name of species
    CHARACTER(LEN=*), OPTIONAL    :: Formula          ! Chemical formula
    REAL(fp),         OPTIONAL    :: MW_g             ! Molecular weight [g]
    REAL(fp),         OPTIONAL    :: EmMW_g           ! Emissions mol. wt [g]
    REAL(fp),         OPTIONAL    :: MolecRatio       ! Molec ratio
    REAL(fp),         OPTIONAL    :: BackgroundVV     ! Background conc [v/v]
    REAL(f8),         OPTIONAL    :: Henry_K0         ! Henry's law K0 [M/atm]
    REAL(f8),         OPTIONAL    :: Henry_CR         ! Henry's law CR [K]
    REAL(f8),         OPTIONAL    :: Henry_PKA        ! Henry's law pKa [1]
    REAL(fp),         OPTIONAL    :: Density          ! Density [kg/m3]
    REAL(fp),         OPTIONAL    :: Radius           ! Radius [m]
    LOGICAL,          OPTIONAL    :: DD_AeroDryDep    ! Use AERO_SFCRSII?
    LOGICAL,          OPTIONAL    :: DD_DustDryDep    ! Use DUST_SFCRSII?
    REAL(fp),         OPTIONAL    :: DD_DvzAerSnow    ! Vd for aerosols
                                                      !  on snow/ice [cm/s]
    REAL(fp),         OPTIONAL    :: DD_DvzMinVal(2)  ! Min Vd for aerosols
                                                      !  (cf GOCART) [cm/s]
    REAL(fp),         OPTIONAL    :: DD_F0            ! Drydep reactivity [1]
    REAL(fp),         OPTIONAL    :: DD_KOA           ! Drydep KOA parameter
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%% NOTE: We will eventually replace this with the common Henry's law
    !%%% parameters.  But in order to replicate the prior behavior,
    !%%% we will need to supply the dry deposition code with the same
    !%%% HSTAR values that are currently set in INIT_DRYDEP.  Therefore,
    !%%% add this field as a temporary placeholder for the Hstar quantity
    !%%% from drydep_mod.F90.  We will remove this later on. (bmy, 8/24/15)
    !%%%
    REAL(fp),         OPTIONAL    :: DD_Hstar_Old     ! HSTAR, drydep [M/atm]
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    LOGICAL,          OPTIONAL    :: MP_SizeResAer    ! Size resolved aerosol?
    LOGICAL,          OPTIONAL    :: MP_SizeResNum    ! Size resolved aer #?
    REAL(fp),         OPTIONAL    :: WD_RetFactor     ! Wetdep retention factor
    LOGICAL,          OPTIONAL    :: WD_LiqAndGas     ! Liquid and gas phases?
    REAL(fp),         OPTIONAL    :: WD_ConvFacI2G    ! Factor for ice/gas ratio
    REAL(fp),         OPTIONAL    :: WD_AerScavEff    ! Aerosol scavenging eff.
    REAL(fp),         OPTIONAL    :: WD_KcScaleFac(3) ! Factor to multiply Kc
                                                      !  rate in F_AEROSOL
    REAL(fp),         OPTIONAL    :: WD_RainoutEff(3) ! Rainout efficiency
    LOGICAL,          OPTIONAL    :: WD_CoarseAer     ! Coarse aerosol?
    LOGICAL,          OPTIONAL    :: Is_Advected      ! Is it advected?
    LOGICAL,          OPTIONAL    :: Is_DryAlt        ! Is it a drydep species
                                                      !  to save at a given alt?
    LOGICAL,          OPTIONAL    :: Is_Drydep        ! Is it dry deposited?
    LOGICAL,          OPTIONAL    :: Is_Gas           ! Gas (T) or aerosol (F)?
    LOGICAL,          OPTIONAL    :: Is_HygroGrowth   ! Is hygroscopic growth?
    LOGICAL,          OPTIONAL    :: Is_Photolysis    ! Is it photolysis spc?
    LOGICAL,          OPTIONAL    :: Is_Wetdep        ! Is it wet deposited?
    LOGICAL,          OPTIONAL    :: Is_InRestart     ! Is it in restart file?
    LOGICAL,          OPTIONAL    :: Is_Hg0           ! Denotes Hg0 species
    LOGICAL,          OPTIONAL    :: Is_Hg2           ! Denotes Hg2 species
    LOGICAL,          OPTIONAL    :: Is_HgP           ! Denotes HgP species
    INTEGER,          OPTIONAL    :: KppSpcId         ! KPP species ID
    INTEGER,          OPTIONAL    :: KppVarId         ! KPP variable species ID
    INTEGER,          OPTIONAL    :: KppFixId         ! KPP fixed species ID
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(Species),    POINTER     :: ThisSpc       ! Object w/ species info
!
! !OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(OUT) :: RC            ! Return code
!
! !REMARKS:
!  (1 ) If Fullname is not specified, it will use the value assigned to Name.
!  (2 ) If EmMw_g is not specified, it will use the value assigned to MW_g.
!  (3 ) If MolecRatio is not specified, it will be set to 1.
!  (4 ) WD_Is_HNO3 is automatically set according to Name.
!  (5 ) WD_Is_SO2 is automatically set according to Name.
!  (6 ) All other fields, if not specified, will be set to MISSING
!  (7 ) If Is_Gas = T, aerosol-specific fields will be set to MISSING.
!        (except for HNO3 and SO2, which wet scavenge as aerosols).
!  (8 ) If Is_Gas = F, gas-phase specific-fields will be set to MISSING.
!  (9 ) If Is_Advected = T, this will automatically update AdvectId.
!  (10) If Is_Drydep = T, this will automatically update DryDepId.
!  (11) If Is_Wetdep = T, this will automatically update WetDepId.
!
! !REVISION HISTORY:
!  20 Aug 2013 - C. Keller - Adapted from gigc_state_chm_mod.F90
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    CHARACTER(LEN=14) :: Name14

    !=====================================================================
    ! Spc_Create begins here!
    !=====================================================================

    ! Check if already allocated
    IF ( ASSOCIATED( ThisSpc ) ) DEALLOCATE( ThisSpc )

    ! Allocate pointer
    ALLOCATE( ThisSpc )

    !---------------------------------------------------------------------
    ! Model ID #
    !---------------------------------------------------------------------
    IF ( PRESENT( ModelID ) ) THEN
       ThisSpc%ModelID = ModelID
    ELSE
       ThisSpc%ModelID = MISSING_INT
    ENDIF

    !---------------------------------------------------------------------
    ! Short name.  Also construct a hash value from the short name.
    !---------------------------------------------------------------------
    IF ( PRESENT( Name ) ) THEN
       ThisSpc%Name     = Name
    ELSE
       ThisSpc%Name     = ''
    ENDIF

    !---------------------------------------------------------------------
    ! Long name (default to short name if not passed)
    !---------------------------------------------------------------------
    IF ( PRESENT( FullName ) ) THEN
       ThisSpc%FullName = FullName
    ELSE
       IF ( PRESENT( Name ) ) THEN
          ThisSpc%FullName = Name
       ELSE
          ThisSpc%FullName = ''
       ENDIF
    ENDIF

    !---------------------------------------------------------------------
    ! Chemical formula
    !---------------------------------------------------------------------
    IF ( PRESENT( Formula ) ) THEN
       ThisSpc%Formula = Formula
    ELSE
       ThisSpc%Formula = ''
    ENDIF

    !---------------------------------------------------------------------
    ! Molecular weight [g]
    !---------------------------------------------------------------------
    IF ( PRESENT( MW_g ) ) THEN
       ThisSpc%MW_g = MW_g
    ELSE
       ThisSpc%MW_g = MISSING_MW
    ENDIF

    !---------------------------------------------------------------------
    ! Emission molecular weight [g]
    ! (Defaults to molecular weight MW_g if not specified)
    !---------------------------------------------------------------------
    IF ( PRESENT( EmMW_g ) ) THEN
       ThisSpc%EmMW_g = EmMW_g
    ELSE
       IF ( PRESENT( MW_g ) ) THEN
          ThisSpc%EmMW_g = MW_g
       ELSE
          ThisSpc%EmMW_g = MISSING_MW
       ENDIF
    ENDIF

    !---------------------------------------------------------------------
    ! Molecule ratio (i.e. moles carbon per moles species)
    ! (Defaults to 1.0 if not specified)
    !---------------------------------------------------------------------
    IF ( PRESENT( MolecRatio ) ) THEN
       ThisSpc%MolecRatio = MolecRatio
    ELSE
       ThisSpc%MolecRatio = 1e+0_fp
    ENDIF

    !---------------------------------------------------------------------
    ! Radius [m]
    !---------------------------------------------------------------------
    IF ( PRESENT( Radius ) ) THEN
       ThisSpc%Radius = Radius
    ELSE
       ThisSpc%Radius = MISSING
    ENDIF

    !---------------------------------------------------------------------
    ! Density [kg/m3]
    !---------------------------------------------------------------------
    IF ( PRESENT( Density ) ) THEN
       ThisSpc%Density = Density
    ELSE
       ThisSpc%Density = MISSING
    ENDIF

    !---------------------------------------------------------------------
    ! DD_AeroDryDep: call routine AERO_SFCRSII to do drydep
    ! (i.e. special drydep handling for sea salt species)
    !---------------------------------------------------------------------
    IF ( PRESENT( DD_AeroDryDep ) ) THEN
       ThisSpc%DD_AeroDryDep = DD_AeroDryDep
    ELSE
       ThisSpc%DD_AeroDryDep = .FALSE.
    ENDIF

    !---------------------------------------------------------------------
    ! DD_DustDryDep: call routine DUST_SFCRSII to do drydep
    ! (i.e. special drydep handling for dust species)
    !---------------------------------------------------------------------
    IF ( PRESENT( DD_DustDryDep ) ) THEN
       ThisSpc%DD_DustDryDep = DD_DustDryDep
    ELSE
       ThisSpc%DD_DustDryDep = .FALSE.
    ENDIF

    !---------------------------------------------------------------------
    ! F0 (stickiness) parameter for drydep
    !---------------------------------------------------------------------
    IF ( PRESENT( DD_F0 ) ) THEN
       ThisSpc%DD_F0 = DD_F0
    ELSE
       ThisSpc%DD_F0 = MISSING
    ENDIF

    !---------------------------------------------------------------------
    ! Dry deposition velocity for aerosols on snow [cm/s]
    !---------------------------------------------------------------------
    IF ( PRESENT( DD_DvzAerSnow ) ) THEN
       ThisSpc%DD_DvzAerSnow = DD_DvzAerSnow
    ELSE
       ThisSpc%DD_DvzAerSnow = MISSING
    ENDIF

    !---------------------------------------------------------------------
    ! Mimimum value for drydep velocity (cf GOCART model) [cm/s]
    !---------------------------------------------------------------------
    IF ( PRESENT( DD_DvzMinVal ) ) THEN
       ThisSpc%DD_DvzMinVal(:) = DD_DvzMinVal(:)
    ELSE
       ThisSpc%DD_DvzMinVal(:) = MISSING
    ENDIF

    !---------------------------------------------------------------------
    ! KOA parameter for drydep (POPs species only)
    !---------------------------------------------------------------------
    IF ( PRESENT( DD_KOA ) ) THEN
       ThisSpc%DD_KOA = DD_KOA
    ELSE
       ThisSpc%DD_KOA = MISSING
    ENDIF

    !---------------------------------------------------------------------
    ! Old HSTAR parameter from drydep_mod
    !---------------------------------------------------------------------
    IF ( PRESENT( DD_Hstar_Old ) ) THEN
       ThisSpc%DD_Hstar_Old = DD_Hstar_old
    ELSE
       ThisSpc%DD_Hstar_Old = MISSING
    ENDIF

    !---------------------------------------------------------------------
    ! Henry's law K0 parameter (aka Hcp)
    !---------------------------------------------------------------------
    IF ( PRESENT( Henry_K0 ) ) THEN
       ThisSpc%Henry_K0 = Henry_K0
    ELSE
       ThisSpc%Henry_K0 = MISSING_R8
    ENDIF

    !---------------------------------------------------------------------
    ! Henry's law CR parameter
    !---------------------------------------------------------------------
    IF ( PRESENT( Henry_CR ) ) THEN
       ThisSpc%Henry_CR = Henry_CR
    ELSE
       ThisSpc%Henry_CR = MISSING_R8
    ENDIF

    !---------------------------------------------------------------------
    ! Henry's law pKA parameter [1]
    !---------------------------------------------------------------------
    IF ( PRESENT( Henry_PKA ) ) THEN
       ThisSpc%Henry_PKA = Henry_PKA
    ELSE
       ThisSpc%Henry_PKA = MISSING_R8
    ENDIF

    !---------------------------------------------------------------------
    ! Retention factor for wetdep (gas-phase species only)
    !---------------------------------------------------------------------
    IF ( PRESENT( WD_RetFactor ) ) THEN
       ThisSpc%WD_RetFactor = WD_RetFactor
    ELSE
       ThisSpc%WD_RetFactor = MISSING
    ENDIF

    !---------------------------------------------------------------------
    ! Use liquid and gas phases for gas-phase species wetdep?
    !---------------------------------------------------------------------
    IF ( PRESENT( WD_LiqAndGas ) ) THEN
       ThisSpc%WD_LiqAndGas = WD_LiqAndGas
    ELSE
       ThisSpc%WD_LiqAndGas = .FALSE.
    ENDIF

    !---------------------------------------------------------------------
    ! Conversion factor for computing the ice/gas ratio for wetdep
    ! (gas-phase species only)
    !---------------------------------------------------------------------
    IF ( PRESENT( WD_ConvFacI2G ) ) THEN
       ThisSpc%WD_ConvFacI2G = WD_ConvFacI2G
    ELSE
       ThisSpc%WD_ConvFacI2G = MISSING
    ENDIF

    !---------------------------------------------------------------------
    ! Scavenging efficiency for wetdep (aerosol species only)
    !---------------------------------------------------------------------
    IF ( PRESENT( WD_AerScavEff ) ) THEN
       ThisSpc%WD_AerScavEff = WD_AerScavEff
    ELSE
       ThisSpc%WD_AerScavEff = MISSING
    ENDIF

    !---------------------------------------------------------------------
    ! Scale factor used to multiply the Kc rate (condensate -> precip)
    ! in routine F_AEROSOL in wetscav_mod.F90.  This implments the
    ! impaction scavenging for aerosol species.
    !---------------------------------------------------------------------
    IF ( PRESENT( WD_KcScaleFac ) ) THEN
       ThisSpc%WD_KcScaleFac(:) = WD_KcScaleFac(:)
    ELSE
       ThisSpc%WD_KcScaleFac(:) = MISSING
    ENDIF

    !---------------------------------------------------------------------
    ! Rainout efficiency for wetdep (aerosol species only)
    !---------------------------------------------------------------------
    IF ( PRESENT( WD_RainOutEff ) ) THEN
       ThisSpc%WD_RainOutEff(:) = WD_RainOutEff(:)
    ELSE
       ThisSpc%WD_RainOutEff(:) = MISSING
    ENDIF

    !---------------------------------------------------------------------
    ! Is it advected?
    !---------------------------------------------------------------------
    IF ( PRESENT( Is_Advected ) ) THEN

       ! Increment the count of advected species
       ThisSpc%Is_Advected = Is_Advected

       ! Update count & index of advected species
       IF ( Is_Advected ) THEN
          AdvectCount         = AdvectCount + 1
          ThisSpc%AdvectID    = AdvectCount
       ELSE
          ThisSpc%Is_Advected = .FALSE.
          ThisSpc%AdvectID    = MISSING_INT
       ENDIF

    ELSE
       ThisSpc%Is_Advected    = .FALSE.
       ThisSpc%AdvectID       = MISSING_INT
    ENDIF

    !---------------------------------------------------------------------
    ! Is it dry-deposited?  If TRUE, then update drydep species index.
    !---------------------------------------------------------------------
    IF ( PRESENT( Is_Drydep ) ) THEN
       ThisSpc%Is_Drydep       = Is_Drydep

       IF ( Is_Drydep ) THEN

          ! Increment count of drydep'd species
          DryDepCount          = DryDepCount + 1

          ! If the dry deposition ID # is passed, then use it;
          ! Otherwise increment the index of drydep'd species
          IF ( PRESENT( DryDepID ) ) THEN
             ThisSpc%DryDepID  = DryDepID
          ELSE
             ThisSpc%DryDepID  = DryDepCount
          ENDIF

       ELSE
          ThisSpc%Is_DryDep    = .FALSE.
          ThisSpc%DryDepID     = MISSING_INT
       ENDIF

    ELSE
       ThisSpc%Is_Drydep       = .FALSE.
       ThisSpc%DryDepID        = MISSING_INT
    ENDIF


    !---------------------------------------------------------------------
    ! Is it a drydep species that we want to save at a given altitude
    ! above the surface?
    !---------------------------------------------------------------------
    IF ( PRESENT( Is_DryAlt ) ) THEN
       ThisSpc%Is_DryAlt = Is_DryAlt

       ! Update count & index
       IF ( Is_DryAlt ) THEN
          DryAltCount      = DryAltCount + 1
          ThisSpc%DryAltID = DryAltCount
       ELSE
          ThisSpc%Is_DryAlt = .FALSE.
          ThisSpc%DryAltID  = MISSING_INT
       ENDIF
    ELSE
       ThisSpc%Is_DryAlt = .FALSE.
       ThisSpc%DryAltID  = MISSING_INT
    ENDIF

    !---------------------------------------------------------------------
    ! Is it a gas?  (If FALSE, then it's an aerosol)
    !---------------------------------------------------------------------
    IF ( PRESENT( Is_Gas ) ) THEN
       ThisSpc%Is_Gas     = Is_Gas

       ! Update count & index of gas species
       IF ( Is_Gas ) THEN
          GasSpcCount      = GasSpcCount + 1
          ThisSpc%GasSpcID = GasSpcCount
       ELSE
          ThisSpc%Is_Gas   = .FALSE.
          ThisSpc%GasSpcID = MISSING_INT
       ENDIF
    ELSE
       ThisSpc%Is_Gas = .FALSE.
       ThisSpc%GasSpcID = MISSING_INT
    ENDIF

    !---------------------------------------------------------------------
    ! Is it an aerosol? (TRUE if not a gas)
    !---------------------------------------------------------------------
    ThisSpc%Is_Aero = .NOT. ThisSpc%Is_Gas
    IF ( ThisSpc%Is_Aero ) THEN
       ! Update count & index of aerosol species
       AeroCount      = AeroCount + 1
       ThisSpc%AeroID = AeroCount
    ELSE
       ThisSpc%AeroID = MISSING_INT
    ENDIF

    !---------------------------------------------------------------------
    ! Is it an aerosol with hygroscopic growth?
    !---------------------------------------------------------------------
    IF ( PRESENT( Is_HygroGrowth ) ) THEN
       ThisSpc%Is_HygroGrowth = Is_HygroGrowth

       ! Increment count & index of hygroscopic growth species
       IF ( Is_HygroGrowth ) THEN
          HygGrthCount = HygGrthCount + 1
          ThisSpc%HygGrthID   = HygGrthCount
       ELSE
          ThisSpc%Is_HygroGrowth = .FALSE.
          ThisSpc%HygGrthID   = MISSING_INT
       ENDIF
    ELSE
       ThisSpc%Is_HygroGrowth = .FALSE.
       ThisSpc%HygGrthID   = MISSING_INT
    ENDIF

    !---------------------------------------------------------------------
    ! Is it a photolysis species in the KPP chemical mechanism?
    !---------------------------------------------------------------------
    IF ( PRESENT( Is_Photolysis ) ) THEN
       ThisSpc%Is_Photolysis = Is_Photolysis

       ! Increment count & index of wet deposited species
       IF ( Is_Photolysis ) THEN
          PhotolCount        = PhotolCount + 1
          ThisSpc%PhotolID   = PhotolCount
       ELSE
          ThisSpc%Is_Photolysis = .FALSE.
          ThisSpc%PhotolID      = MISSING_INT
       ENDIF

    ELSE
       ThisSpc%Is_Photolysis = .FALSE.
       ThisSpc%PhotolID      = MISSING_INT
    ENDIF

    !---------------------------------------------------------------------
    ! Is it wet deposited?  If TRUE, then update wetdep species index.
    !---------------------------------------------------------------------
    IF ( PRESENT( Is_Wetdep ) ) THEN
       ThisSpc%Is_Wetdep    = Is_Wetdep

       ! Increment count & index of wet deposited species
       IF ( Is_WetDep ) THEN
          WetDepCount       = WetDepCount + 1
          ThisSpc%WetDepID  = WetDepCount
       ELSE
          ThisSpc%Is_Wetdep = .FALSE.
          ThisSpc%WetDepID  = MISSING_INT
       ENDIF

    ELSE
       ThisSpc%Is_Wetdep    = .FALSE.
       ThisSpc%WetDepID     = MISSING_INT
    ENDIF

    !---------------------------------------------------------------------
    ! Is there a default background concentration for this species
    ! [mol spc/mol dry air]? If not, use a default value.
    !---------------------------------------------------------------------
    IF ( PRESENT( BackgroundVV ) ) THEN
       ThisSpc%BackgroundVV = BackgroundVV
    ELSE
       ThisSpc%BackgroundVV = MISSING_VV
    ENDIF

    !---------------------------------------------------------------------
    ! Is it a species in the KPP chemical mechanism?
    !---------------------------------------------------------------------
    IF ( PRESENT( KppSpcId ) ) THEN
       KppSpcCount      = KppSpcCount + 1
       ThisSpc%KppSpcId = KppSpcId
    ELSE
       ThisSpc%KppSpcId = MISSING_INT
    ENDIF

    !---------------------------------------------------------------------
    ! Is it a variable species in the KPP chemical mechanism?
    !---------------------------------------------------------------------
    IF ( PRESENT( KppVarId ) ) THEN
       KppVarCount      = KppVarCount + 1
       ThisSpc%KppVarId = KppVarId
    ELSE
       ThisSpc%KppVarId = MISSING_INT
    ENDIF

    !---------------------------------------------------------------------
    ! Is it a fixed species in the KPP chemical mechanism?
    !---------------------------------------------------------------------
    IF ( PRESENT( KppFixId ) ) THEN
       KppFixCount      = KppFixCount + 1
       ThisSpc%KppFixId = KppFixId
    ELSE
       ThisSpc%KppFixId = MISSING_INT
    ENDIF

    !---------------------------------------------------------------------
    ! Is it stored in the restart file?
    !---------------------------------------------------------------------
    IF ( PRESENT( Is_InRestart ) ) THEN

       ! Assume presence in the restart file until proven otherwise
       ! within READ_GC_RESTART
       ! NOTE: Eventually this logical flag can be removed since all
       ! species will be stored in the restart file (7/25/16)
       ThisSpc%Is_InRestart = .TRUE.
    ELSE
       ThisSpc%Is_InRestart = .FALSE.
    ENDIF

    !---------------------------------------------------------------------
    ! Is it a coarse aerosol?
    !---------------------------------------------------------------------
    IF ( PRESENT( WD_CoarseAer ) ) THEN
       ThisSpc%WD_CoarseAer = WD_CoarseAer
    ELSE
       ThisSpc%WD_CoarseAer = .FALSE.
    ENDIF

    !---------------------------------------------------------------------
    ! Is it a size-resolved aerosol species for microphysics?
    !---------------------------------------------------------------------
    IF ( PRESENT( MP_SizeResAer ) ) THEN
       ThisSpc%MP_SizeResAer = MP_SizeResAer
    ELSE
       ThisSpc%MP_SizeResAer = .FALSE.
    ENDIF

    !---------------------------------------------------------------------
    ! Is it a size-resolved bin number?
    !---------------------------------------------------------------------
    IF ( PRESENT( MP_SizeResNum ) ) THEN
       ThisSpc%MP_SizeResNum = MP_SizeResNum
    ELSE
       ThisSpc%MP_SizeResNum = .FALSE.
    ENDIF

    !---------------------------------------------------------------------
    ! Is it a Hg0 species (total or tagged)?
    !---------------------------------------------------------------------
    IF ( PRESENT( Is_Hg0 ) ) THEN
       ThisSpc%Is_Hg0 = Is_Hg0

       ! Increment count and index of Hg0 categories
       IF ( Is_Hg0 ) THEN
          Hg0Count       = Hg0Count + 1
          ThisSpc%Hg_Cat = Hg0Count
       ELSE
          ThisSpc%Is_Hg0 = .FALSE.
       ENDIF

    ELSE
       ThisSpc%Is_Hg0    = .FALSE.
    ENDIF

    !---------------------------------------------------------------------
    ! Is it a Hg2 species (total or tagged)?
    !---------------------------------------------------------------------
    IF ( PRESENT( Is_Hg2 ) ) THEN
       ThisSpc%Is_Hg2 = Is_Hg2

       ! Increment count of Hg2 species
       IF ( Is_Hg2 ) THEN
          Hg2Count       = Hg2Count + 1
          ThisSpc%Hg_Cat = Hg2Count
       ELSE
          ThisSpc%Is_Hg2 = .FALSE.
       ENDIF

    ELSE
       ThisSpc%Is_Hg2    = .FALSE.
    ENDIF

    !---------------------------------------------------------------------
    ! Is it a HgP species (total or tagged)?
    !---------------------------------------------------------------------
    IF ( PRESENT( Is_HgP ) ) THEN
       ThisSpc%Is_HgP = Is_HgP

       ! Increment count of HgP species
       IF ( Is_HgP ) THEN
          HgPCount       = HgPCount + 1
          ThisSpc%Hg_Cat = HgPCount
       ELSE
          ThisSpc%Is_HgP = .FALSE.
       ENDIF

    ELSE
       ThisSpc%Is_HgP    = .FALSE.
    ENDIF

    !---------------------------------------------------------------------
    ! Sanity checks
    !---------------------------------------------------------------------

    ! Is the species part of the KPP chemical mechanism?
    ThisSpc%Is_Kpp = ( ThisSpc%KppVarId > 0 .or. ThisSpc%KppFixId > 0 )

    ! Is the species an active or fixed species in the chemical mechanism?
    ThisSpc%Is_ActiveChem = ( ThisSpc%KppVarId >  0 .and. &
                              ThisSpc%KppFixId <= 0 )
    ThisSpc%Is_FixedChem  = ( ThisSpc%KppFixId >  0 )

    ! Assume the species is not H2SO4, HNO3, or SO2 at first
    ThisSpc%WD_Is_H2SO4 = .FALSE.
    ThisSpc%WD_Is_HNO3  = .FALSE.
    ThisSpc%WD_Is_SO2   = .FALSE.

    IF ( ThisSpc%Is_Gas ) THEN

       ! If this is a gas, then zero out all aerosol fields ...
       ! ... except for those species that wetdep like aerosols
       SELECT CASE( TRIM( ThisSpc%Name ) )
          CASE( 'H2SO4' )
             ThisSpc%WD_Is_H2SO4      = .TRUE.   ! Set flag for H2SO4
          CASE( 'HNO3' )
             ThisSpc%DD_DvzAerSnow    = MISSING
             ThisSpc%MP_SizeResAer    = .FALSE.
             ThisSpc%MP_SizeResNum    = .FALSE.
             ThisSpc%WD_CoarseAer     = .FALSE.
             ThisSpc%WD_Is_HNO3       = .TRUE.   ! Set flag for HNO3
          CASE( 'SO2' )
             ThisSpc%DD_DvzAerSnow    = MISSING
             ThisSpc%MP_SizeResAer    = .FALSE.
             ThisSpc%MP_SizeResNum    = .FALSE.
             ThisSpc%WD_CoarseAer     = .FALSE.
             ThisSpc%WD_Is_SO2        = .TRUE.   ! Set flag for SO2
          CASE DEFAULT
             ThisSpc%DD_DvzAerSnow    = MISSING
             ThisSpc%MP_SizeResAer    = .FALSE.
             ThisSpc%MP_SizeResNum    = .FALSE.
             ThisSpc%WD_CoarseAer     = .FALSE.
             ThisSpc%WD_AerScavEff    = MISSING
             ThisSpc%WD_KcScaleFac(:) = MISSING
             ThisSpc%WD_RainoutEff(:) = MISSING
       END SELECT

    ELSE

       ! If this species is an aerosol, zero out gas-phase fields
       ThisSpc%WD_RetFactor  = MISSING
       ThisSpc%WD_LiqAndGas  = .FALSE.
       ThisSpc%WD_ConvFacI2G = MISSING

    ENDIF

  END SUBROUTINE Spc_Create
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Spc_Print
!
! !DESCRIPTION: Routine Spc\_Create prints the fields of the species object.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Spc_Print( Input_Opt, ThisSpc, RC )
!
! !USES:
!
    USE Input_Opt_Mod, ONLY : OptInput
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput),   INTENT(IN)    :: Input_Opt    ! Input Options object
    TYPE(Species),    POINTER       :: ThisSpc      ! Object w/ species info
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(INOUT) :: RC           ! Return code
!
! !REMARKS:
!  Optional fields are not printed out if they are not defined (i.e. if they
!  have a "missing data value" of -999).
!
! !REVISION HISTORY:
!  27 Jul 2015 - R. Yantosca - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC

    !=====================================================================
    ! Spc_Create begins here!
    !=====================================================================
    IF ( Input_Opt%amIRoot ) THEN

       !-------------------------
       ! Print general info
       !-------------------------
       WRITE( 6, '(a)' ) REPEAT( '=', 79 )
       WRITE( 6, 100 ) 'Species ID            ',  ThisSpc%ModelID
       WRITE( 6, 110 ) 'Name                  ',  TRIM( ThisSpc%Name     )
       WRITE( 6, 110 ) 'FullName              ',  TRIM( ThisSpc%FullName )
       WRITE( 6, 110 ) 'Formula               ',  TRIM( ThisSpc%Formula  )
       WRITE( 6, 120 ) 'Molecular weight [g]  ',  ThisSpc%MW_g
       WRITE( 6, 120 ) 'Emitted mol. wt [g]   ',  ThisSpc%EmMW_g
       WRITE( 6, 120 ) 'Molecular ratio       ',  ThisSpc%MolecRatio
       WRITE( 6, 130 ) 'Is it a gas?          ',  ThisSpc%Is_Gas
       WRITE( 6, 120 ) 'Default bckgrnd [v/v] ',  ThisSpc%BackgroundVV

       !-------------------------
       ! Print density & radius
       !-------------------------
       IF ( ThisSpc%Density > ZERO ) THEN
          WRITE( 6, 120 ) 'Density [kg/m3]      ',       ThisSpc%Density
       ENDIF

       IF ( ThisSpc%Radius > ZERO ) THEN
          WRITE( 6, 120 ) 'Radius [m]           ',       ThisSpc%Radius
       ENDIF

       !-------------------------
       ! Print Henry's Law info
       !-------------------------
       IF ( ThisSpc%Henry_K0 > ZERO_R8 ) THEN
          WRITE( 6, 120 ) 'Henry''s law K0       ',      ThisSpc%Henry_K0
       ENDIF

       IF ( ThisSpc%Henry_CR > ZERO_R8 ) THEN
          WRITE( 6, 120 ) 'Henry''s law CR       ',      ThisSpc%Henry_CR
       ENDIF

       IF ( ThisSpc%Henry_pKa > ZERO_R8 ) THEN
          WRITE( 6, 120 ) 'Henry''s law pKa      ',      ThisSpc%Henry_pKa
       ENDIF

       !-------------------------
       ! Print advected ID
       !-------------------------
       WRITE( 6, 130 ) 'Is it advected?      ',          ThisSpc%Is_Advected

       IF ( ThisSpc%Is_Advected ) THEN
          WRITE( 6, 100 )    ' -> Advected index   ',    ThisSpc%AdvectId
       ENDIF

       !-------------------------
       ! Print KPP Id's
       !-------------------------
       WRITE( 6, 130 ) 'Is it a KPP species? ',          ThisSpc%Is_Kpp

       IF ( ThisSpc%Is_Kpp ) THEN
          IF ( ThisSpc%KppSpcId > 0 ) THEN
             WRITE( 6, 100 )    ' -> ID in C   array  ', ThisSpc%KppSpcId
          ENDIF

          WRITE( 6, 130 ) 'Is it an active spc? ',       ThisSpc%Is_ActiveChem
          IF ( ThisSpc%KppVarId > 0 ) THEN
             WRITE( 6, 100 )    ' -> ID in VAR array  ', ThisSpc%KppVarId
          ENDIF

          WRITE( 6, 130 ) 'Is it a fixed spc?   ',       ThisSpc%Is_FixedChem
          IF ( ThisSpc%KppFixId > 0 ) THEN
             WRITE( 6, 100 )    ' -> ID in FIX array  ', ThisSpc%KppFixId
          ENDIF
       ENDIF

       !-------------------------
       ! Print photolysis info
       !-------------------------
       WRITE( 6, 130 ) 'Is it a photol. spc? ',          ThisSpc%Is_Photolysis

       !------------------------------
       ! Print hygroscopic growth info
       !------------------------------
       WRITE( 6, 130 ) 'Is it a hygro. spc? ',          ThisSpc%Is_HygroGrowth

       !-------------------------
       ! Print drydep info
       !-------------------------
       WRITE( 6, 130 ) 'Is it dry deposited? ',          ThisSpc%Is_DryDep

       IF ( ThisSpc%Is_DryDep ) THEN
          WRITE( 6, 100 ) ' -> Drydep index     ',       ThisSpc%DryDepID

          IF ( ThisSpc%DD_F0 > ZERO ) THEN
             WRITE( 6, 120 ) ' -> F0 parameter     ',    ThisSpc%DD_F0
          ENDIF

          IF ( ThisSpc%DD_KOA > ZERO ) THEN
             WRITE( 6, 120 ) ' -> KOA parameter    ',    ThisSpc%DD_KOA
          ENDIF

          IF ( ThisSpc%DD_Hstar_old > ZERO ) THEN
             WRITE( 6, 120 ) ' -> Old HSTAR value  ',    ThisSpc%DD_Hstar_Old
          ENDIF
       ENDIF

       !-------------------------
       ! Print wetdep info
       !-------------------------
       WRITE( 6, 130 ) 'Is it wet deposited? ',          ThisSpc%Is_WetDep

       IF ( ThisSpc%Is_WetDep ) THEN
          WRITE( 6, 100 ) ' -> Wetdep index:    ',       ThisSpc%WetDepID

          IF ( ThisSpc%WD_LiqAndGas ) THEN
             WRITE( 6, 130 ) ' -> Liq & gas phases ',    ThisSpc%WD_LiqAndGas
             WRITE( 6, 120 ) ' -> Conv factor I2G  ',    ThisSpc%WD_ConvFacI2G
          ENDIF

          IF ( ThisSpc%WD_RetFactor > ZERO ) THEN
             WRITE( 6, 120 ) ' -> Ret. Factor      ',    ThisSpc%WD_RetFactor
          ENDIF

          IF ( ThisSpc%WD_CoarseAer ) THEN
             WRITE( 6, 130 ) ' -> Coarse aerosol?',      ThisSpc%WD_CoarseAer
          ENDIF

          IF ( ThisSpc%WD_AerScavEff > ZERO ) THEN
             WRITE( 6, 120 ) ' -> Scav. Effeciency ',    ThisSpc%WD_AerScavEff
          ENDIF

          IF ( SUM( ThisSpc%WD_KcScaleFac ) > ZERO ) THEN
             WRITE( 6, 140 ) ' -> KcScale factor   ',    ThisSpc%WD_KcScaleFac
          ENDIF

          IF ( SUM( ThisSpc%WD_RainoutEff ) > ZERO ) THEN
             WRITE( 6, 140 ) ' -> Rainout effic''ncy',   ThisSpc%WD_RainoutEff
          ENDIF

       ENDIF

       !-------------------------
       ! Print microphys info
       !-------------------------

       ! Format statements
 100   FORMAT( a, ' : ', i5          )
 110   FORMAT( a, ' : ', a           )
 120   FORMAT( a, ' : ', es13.6      )
 130   FORMAT( a, ' : ', L1          )
 140   FORMAT( a, ' : ', 3(f5.1, 1x) )
    ENDIF

  END SUBROUTINE Spc_Print
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Spc_GetNumSpecies
!
! !DESCRIPTION: Routine Spc\_GetNumSpecies returns the number of advected,
!  dry-deposited, and wet-deposited species to an external routine.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Spc_GetNumSpecies( nAdvect,  nAeroSpc, nDryAlt, nDryDep,        &
                                nGasSpc,  nHygGrth, nKppVar, nKppFix,        &
                                nKppSpc,  nPhotol,  nWetDep, nHg0Cats,       &
                                nHg2Cats, nHgPCats                          )
!
! !OUTPUT PARAMETERS:
!
    INTEGER, INTENT(OUT) :: nAdvect     ! # of advected species
    INTEGER, INTENT(OUT) :: nAeroSpc    ! # of aerosol species
    INTEGER, INTENT(OUT) :: nDryAlt     ! # of dry-dep species to save at a
                                        !  user-defined altitude above sfc.
    INTEGER, INTENT(OUT) :: nDryDep     ! # of dry-deposited species
    INTEGER, INTENT(OUT) :: nGasSpc     ! # of gas-phase species
    INTEGER, INTENT(OUT) :: nHygGrth    ! # of species with hygroscopic growth
    INTEGER, INTENT(OUT) :: nKppVar     ! # of variable KPP species
    INTEGER, INTENT(OUT) :: nKppFix     ! # of fixed KPP species
    INTEGER, INTENT(OUT) :: nKppSpc     ! # of species in KPP mechanism
    INTEGER, INTENT(OUT) :: nPhotol     ! # of photolysis species
    INTEGER, INTENT(OUT) :: nWetDep     ! # of wet-deposited species
    INTEGER, INTENT(OUT) :: nHg0Cats    ! # of Hg0 categories
    INTEGER, INTENT(OUT) :: nHg2Cats    ! # of Hg0 categories
    INTEGER, INTENT(OUT) :: nHgPCats    ! # of Hg0 categories
!
! !REVISION HISTORY:
!  02 Sep 2015 - R. Yantosca - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC

    ! Return module variables
    nAdvect  = AdvectCount
    nAeroSpc = AeroCount
    nDryAlt  = DryAltCount
    nDryDep  = DryDepCount
    nGasSpc  = GasSpcCount
    nHygGrth = HygGrthCount
    nKppVar  = KppVarCount
    nKppFix  = KppFixCount
    nKppSpc  = KppSpcCount
    nPhotol  = PhotolCount
    nWetDep  = WetDepCount
    nHg0Cats = Hg0Count
    nHg2Cats = Hg2Count
    nHgPCats = HgPCount

  END SUBROUTINE Spc_GetNumSpecies
!EOC
END MODULE Species_Mod

