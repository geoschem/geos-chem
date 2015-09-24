!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: species_mod 
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
  PUBLIC :: SpcData_GetIndx
  PUBLIC :: SpcData_Cleanup
  PUBLIC :: Spc_Create
  PUBLIC :: Spc_Print
  PUBLIC :: Spc_GetNumSpecies
!
! !PUBLIC TYPES: 
!
  !=========================================================================
  ! Counters for the species indices
  !=========================================================================
  INTEGER, PRIVATE :: AdvectCount = 0    ! Counter of advected species
  INTEGER, PRIVATE :: DryDepCount = 0    ! Counter of dry-deposited species
  INTEGER, PRIVATE :: WetDepCount = 0    ! Counter of wet-deposited species

  !=========================================================================
  ! Type for ASCII sums (fast-species lookup algorithm)
  !=========================================================================
  TYPE, PUBLIC :: AsciiSum
     INTEGER,           POINTER :: I(:)  !
     CHARACTER(len=14), POINTER :: N(:)  !
  END Type AsciiSum

  !=========================================================================
  ! Type for the species database object (vector of type Species)
  !=========================================================================
  TYPE, PUBLIC :: SpcPtr
     TYPE(Species),     POINTER :: Info  ! Species type
  END TYPE SpcPtr

  !=========================================================================
  ! Type for individual species information
  !=========================================================================
  TYPE, PUBLIC :: Species

     ! Indices
     INTEGER            :: ModelID       ! Model species ID
     INTEGER            :: AdvectID      ! Advection index
     INTEGER            :: DryDepID      ! Dry deposition index
     INTEGER            :: WetDepID      ! Wet deposition index

     ! Names
     CHARACTER(LEN=31)  :: Name          ! Short name
     CHARACTER(LEN=80)  :: FullName      ! Long name

     ! Logical switches
     LOGICAL            :: Is_Gas        ! Is it a gas?  If not, aerosol.
     LOGICAL            :: Is_Advected   ! Is it advected?
     LOGICAL            :: Is_DryDep     ! Is it dry-deposited?
     LOGICAL            :: Is_WetDep     ! Is it wet-deposited?

     ! Molecular weights and conversion factors
     REAL(fp)           :: MW_g          ! Species molecular weight  [g/mole]
     REAL(fp)           :: EmMW_g        ! Emitted molecular weight  [g/mole]
     REAL(fp)           :: MolecRatio    ! Molecule emission ratio   [1     ] 
     !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     !%%% NOTE: These parameters are kept for backwards compatibility 
     !%%% but will be removed when the new unit conversion modifications
     !%%% are introduced into GEOS-Chem. (bmy, 9/2/15)
     !%%%
     REAL(fp)           :: TCVV          ! Mol. Wt. dry air / Mol. wt. species
     REAL(fp)           :: XNUMOL        ! Molecules species / kg species
     !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

     ! Henry's law parameters
     REAL(f8)           :: Henry_K0      ! Liq. over gas Henry const [M/atm ]
     REAL(f8)           :: Henry_CR      ! d(ln K0) / d(1/T)         [K     ] 
     REAL(f8)           :: Henry_PKA     ! pKa for Henry const. correction

     ! Drydep parameters
     REAL(fp)           :: DD_A_Density  ! Aerosol density [kg/m3]
     REAL(fp)           :: DD_A_Radius   ! Aerosol radius  [um]
     REAL(fp)           :: DD_F0         ! F0 (reactivity) factor [1]
     REAL(fp)           :: DD_KOA        ! KOA factor
     !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     !%%% NOTE: We will eventually replace this with the common Henry's law
     !%%% parameters.  But in order to replicate the prior behavior,
     !%%% we will need to supply the dry deposition code with the same
     !%%% HSTAR values that are currently set in INIT_DRYDEP.  Therefore,
     !%%% add this field as a temporary placeholder for the Hstar quantity
     !%%% from drydep_mod.F.  We will remove this later on. (bmy, 8/24/15)
     !%%%
     REAL(fp)           :: DD_Hstar_Old     ! HSTAR value in drydep_mod [M/atm]
     !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     
     ! Wetdep parameters, gas-phase species
     REAL(fp)           :: WD_RetFactor     ! Retention factor [1]
     LOGICAL            :: WD_LiqAndGas     ! Consider liquid and gas phases?
     REAL(fp)           :: WD_ConvFacI2G    ! Conv. factor for ice/gas ratio 

     ! Wetdep parameters, aerosol-phase species
     REAL(fp)           :: WD_AerScavEff    ! Aerosol scavenging efficiency
     LOGICAL            :: WD_CoarseAer     ! T=coarse aerosol; F=fine aerosol
     REAL(fp)           :: WD_KcScaleFac(3) ! Factor used to multiply Kc rate
                                            !  (conv of condensate -> precip)
                                            !  in F_AEROSOL (wetscav_mod.F)
     REAL(fp)           :: WD_RainoutEff(3) ! Updraft scavenging efficiency
     LOGICAL            :: WD_SizeResAer    ! T=size-resolved aerosol (TOMAS)

  END TYPE Species
!
! !DEFINED PARAMETERS
!
  !=========================================================================
  ! Missing value parameters
  !=========================================================================
  INTEGER,  PARAMETER  :: MISSING_INT = -999         ! Integer 
  REAL(fp), PARAMETER  :: MISSING     = -999e+0_fp   ! Flexible precision
  REAL(f8), PARAMETER  :: MISSING_R8  = -999e+0_f8   ! 8-byte precision
  REAL(fp), PARAMETER  :: ZERO        =  0.0e+0_fp   ! Flexible precision
  REAL(f8), PARAMETER  :: ZERO_R8     =  0.0e+0_f8   ! 8-byte precision
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
!  22 Jul 2015 - R. Yantosca - Updated and cleaned up a bit 
!  18 Aug 2015 - R. Yantosca - Added indices for drydep, wetdep, transport
!  18 Aug 2015 - R. Yantosca - Added missing value parameters
!  31 Aug 2015 - R. Yantosca - Add AdvectId
!  24 Sep 2015 - R. Yantosca - Make WD_RainoutEff a 3-element vector:
!                              (1) T < 237K; (2) 237K < T < 258 K (3) T > 258K
!  24 Sep 2015 - R. Yantosca - Rename WD_ConvFactor to WD_ConvFacI2G

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
  SUBROUTINE SpcData_Init( am_I_Root, nSpecies, SpecDB, RC )
!
! !INPUT PARAMETERS:
! 
    LOGICAL,              INTENT(IN)    :: am_I_Root    ! root CPU?
    INTEGER,              INTENT(IN)    :: nSpecies     ! # of species 
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(SpcPtr),         POINTER       :: SpecDB(:)    ! Species database
    INTEGER,              INTENT(INOUT) :: RC           ! Return code
! 
! !REVISION HISTORY: 
!  20 Aug 2013 - C. Keller   - Adapted from gigc_state_chm_mod.F90
!  22 Jul 2015 - R. Yantosca - Cosmetic changes
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER :: N, AS

    !=====================================================================
    ! SpcPtr_Init begins here!
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
! !IROUTINE: SpcData_GetIndx 
!
! !DESCRIPTION: Function SpcData\_GetIndx returns the index of a given 
!  species in the species data base object.  You can search by the short
!  name or the full name of the species.
!\\
!\\
! !INTERFACE:
!
  FUNCTION SpcData_GetIndx( name, SpecDB, fullname ) RESULT( Indx )
!
! !INPUT PARAMETERS:
!
    CHARACTER(LEN=*),     INTENT(IN) :: name       ! Species or tracer name
    TYPE(SpcPtr),         POINTER    :: SpecDB(:)  ! Species database
    LOGICAL, INTENT(IN),  OPTIONAL   :: fullname   ! Search in fullname?
!
! !RETURN VALUE:
!
    INTEGER                          :: Indx       ! Index of this species 
!
! !REVISION HISTORY: 
!  09 Oct 2012 - M. Long     - Initial version, based on gc_esmf_utils_mod.F90
!  22 Jul 2015 - R. Yantosca - Cosmetic changes
!EOP
!------------------------------------------------------------------------------
!BOC
    INTEGER :: M
    LOGICAL :: do_fullname

    ! Initialize
    Indx= -1
    IF ( PRESENT(fullname) ) THEN
       do_fullname = fullname
    ELSE
       do_fullname = .FALSE.
    ENDIF

    ! Loop over all species names
    DO M = 1, SIZE(SpecDB)

       ! Return the index of the sought-for species
       IF ( do_fullname ) THEN
          IF( TRIM( name ) == TRIM( SpecDB(M)%Info%FullName ) ) THEN
             Indx = SpecDB(M)%Info%ModelID
             EXIT
          ENDIF

       ELSE
          IF( TRIM( name ) == TRIM( SpecDB(M)%Info%Name ) ) THEN
             Indx = SpecDB(M)%Info%ModelID
             EXIT
          ENDIF
       ENDIF
    ENDDO

  END FUNCTION SpcData_GetIndx 
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
    TYPE(SpcPtr), POINTER :: SpecDB(:)  ! Species database
! 
! !REVISION HISTORY: 
!  20 Aug 2013 - C. Keller   - Adapted from gigc_state_chm_mod.F90
!  22 Jul 2015 - R. Yantosca - Cosmetic changes
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER :: N

    !=====================================================================
    ! SpcData_Cleanup begins here!
    !=====================================================================

    ! Check if already allocated
    IF ( ASSOCIATED( SpecDB ) ) THEN

       ! Nullify each species pointer
       DO N = 1, SIZE( SpecDB )
          IF( ASSOCIATED( SpecDB(N)%Info ) ) THEN
             DEALLOCATE( SpecDB(N)%Info )
          ENDIF
       ENDDO
       DEALLOCATE( SpecDB )
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
  SUBROUTINE Spc_Create( am_I_Root,     ThisSpc,       ModelID,        &
                         DryDepID,      Name,          FullName,       &
                         MW_g,          EmMW_g,        MolecRatio,     &
                         Henry_K0,      Henry_CR,      Henry_PKA,      &
                         DD_A_Density,  DD_A_Radius,   DD_F0,          &
                         DD_KOA,        DD_HStar_Old,  WD_RetFactor,   &
                         WD_LiqAndGas,  WD_ConvFacI2G, WD_AerScavEff,  &
                         WD_KcScaleFac, WD_RainoutEff, WD_CoarseAer,   &
                         WD_SizeResAer, Is_Advected,   Is_Gas,         &
                         Is_Drydep,     Is_Wetdep,     RC             )
!
! !USES:
!
    USE CMN_GCTM_Mod, ONLY : AIRMW, AVO              ! Physical constants
!
! !INPUT PARAMETERS:
! 
    LOGICAL,          INTENT(IN)  :: am_I_Root        ! Are we on the root CPU?
    INTEGER,          OPTIONAL    :: ModelID          ! Model ID number
    INTEGER,          OPTIONAL    :: DryDepID         ! Drydep ID number
    CHARACTER(LEN=*), OPTIONAL    :: Name             ! Short name of species
    CHARACTER(LEN=*), OPTIONAL    :: FullName         ! Long name of species
    REAL(fp),         OPTIONAL    :: MW_g             ! Molecular weight [g]
    REAL(fp),         OPTIONAL    :: EmMW_g           ! Emissions mol. wt [g]
    REAL(fp),         OPTIONAL    :: MolecRatio       ! Molec ratio
    REAL(f8),         OPTIONAL    :: Henry_K0         ! Henry's law K0 [M/atm]
    REAL(f8),         OPTIONAL    :: Henry_CR         ! Henry's law CR [K]
    REAL(f8),         OPTIONAL    :: Henry_PKA        ! Henry's law pKa
    REAL(fp),         OPTIONAL    :: DD_A_Density     ! Aerosol density [kg/m3]
    REAL(fp),         OPTIONAL    :: DD_A_Radius      ! Aerosol radius [m]
    REAL(fp),         OPTIONAL    :: DD_F0            ! Drydep reactivity
    REAL(fp),         OPTIONAL    :: DD_KOA           ! Drydep KOA parameter
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%% NOTE: We will eventually replace this with the common Henry's law
    !%%% parameters.  But in order to replicate the prior behavior,
    !%%% we will need to supply the dry deposition code with the same
    !%%% HSTAR values that are currently set in INIT_DRYDEP.  Therefore,
    !%%% add this field as a temporary placeholder for the Hstar quantity
    !%%% from drydep_mod.F.  We will remove this later on. (bmy, 8/24/15)
    !%%%
    REAL(fp),         OPTIONAL    :: DD_Hstar_Old     ! HSTAR, drydep [M/atm]
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    REAL(fp),         OPTIONAL    :: WD_RetFactor     ! Wetdep retention factor
    LOGICAL,          OPTIONAL    :: WD_LiqAndGas     ! Liquid and gas phases?
    REAL(fp),         OPTIONAL    :: WD_ConvFacI2G    ! Factor for ice/gas ratio
    REAL(fp),         OPTIONAL    :: WD_AerScavEff    ! Aerosol scavenging eff.
    REAL(fp),         OPTIONAL    :: WD_KcScaleFac(3) ! Factor to multiply Kc
                                                      !  rate in F_AEROSOL
    REAL(fp),         OPTIONAL    :: WD_RainoutEff(3) ! Rainout efficiency
    LOGICAL,          OPTIONAL    :: WD_CoarseAer     ! Coarse aerosol?
    LOGICAL,          OPTIONAL    :: WD_SizeResAer    ! Size resolved aerosol?
    LOGICAL,          OPTIONAL    :: Is_Advected      ! Is it advected?
    LOGICAL,          OPTIONAL    :: Is_Gas           ! Gas (T) or aerosol (F)?
    LOGICAL,          OPTIONAL    :: Is_Drydep        ! Is it dry deposited?
    LOGICAL,          OPTIONAL    :: Is_Wetdep        ! Is it wet deposited?
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
!  (1) If Fullname   is not specified, it will use the value assigned to Name.
!  (2) If EmMw_g     is not specified, it will use the value assigned to MW_g.
!  (3) If MolecRatio is not specified, it will be set to 1.
!  (4) If WD_RainoutEff is not specified, it will default to WD_AerScavEff.
!  (4) All other fields, if not specified, will be set to -999 (missing value).
!
! !REVISION HISTORY: 
!  20 Aug 2013 - C. Keller   - Adapted from gigc_state_chm_mod.F90
!  22 Jul 2015 - R. Yantosca - Added RetFactor and drydep parameters
!  31 Aug 2015 - R. Yantosca - Now also compute AdvectId
!  04 Sep 2015 - R. Yantosca - Add arguments WD_RainoutEff, WD_CoarseAer,
!                              and WD_SizeResAer
!  24 Sep 2015 - R. Yantosca - Added WD_KcScaleFac argument
!EOP
!------------------------------------------------------------------------------
!BOC

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
       ThisSpc%ModelID = -1
    ENDIF

    !---------------------------------------------------------------------
    ! Short name
    !---------------------------------------------------------------------
    IF ( PRESENT( Name ) ) THEN
       ThisSpc%Name = Name
    ELSE
       ThisSpc%Name = ''
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
    ! Molecular weight [g]
    !---------------------------------------------------------------------
    IF ( PRESENT( MW_g ) ) THEN
       ThisSpc%MW_g = MW_g 
    ELSE
       ThisSpc%MW_g = MISSING
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
          ThisSpc%EmMW_g = MISSING
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
    ! Aerosol radius for drydep
    !---------------------------------------------------------------------
    IF ( PRESENT( DD_A_Radius ) ) THEN
       ThisSpc%DD_A_Radius = DD_A_Radius
    ELSE
       ThisSpc%DD_A_Radius = MISSING
    ENDIF

    !---------------------------------------------------------------------
    ! Aerosol density for drydep
    !---------------------------------------------------------------------
    IF ( PRESENT( DD_A_Density ) ) THEN
       ThisSpc%DD_A_Density = DD_A_Density
    ELSE
       ThisSpc%DD_A_Density = MISSING
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
    ! in routine F_AEROSOL in wetscav_mod.F.  This implments the 
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
    ! Is it a gas?  (If FALSE, then it's an aerosol)
    !---------------------------------------------------------------------
    IF ( PRESENT( Is_Gas ) ) THEN
       ThisSpc%Is_Gas = Is_Gas
    ELSE
       ThisSpc%Is_Gas = .FALSE.
    ENDIF

    !---------------------------------------------------------------------
    ! Is it advected?
    !---------------------------------------------------------------------
    IF ( PRESENT( Is_Advected ) ) THEN

       ! Increment the count of advected species
       ThisSpc%Is_Advected = Is_Advected
       
       ! Update count & index of advected species
       IF ( Is_Advected ) THEN
          AdvectCount      = AdvectCount + 1
          ThisSpc%AdvectID = AdvectCount
       ENDIF

    ELSE
       ThisSpc%Is_Advected = .FALSE.
    ENDIF

    !---------------------------------------------------------------------
    ! Is it dry-deposited?  If TRUE, then update drydep species index.
    !---------------------------------------------------------------------
    IF ( PRESENT( Is_Drydep ) ) THEN
       ThisSpc%Is_Drydep       = Is_Drydep

       IF( Is_Drydep ) THEN

          ! Increment count of drydep'd species
          DryDepCount          = DryDepCount + 1

          ! If the dry deposition ID # is passed, then use it;
          ! Otherwise increment the index of drydep'd species
          IF ( PRESENT( DryDepID ) ) THEN 
             ThisSpc%DryDepID  = DryDepID
          ELSE
             ThisSpc%DryDepID  = DryDepCount
          ENDIF

          ! ISOPN dry deposits as ISOPND + ISOPNB
          ! MMN   dry deposits as MACRN  + MVKN
          ! So we need to increment the drydep counter to
          ! leave space for the next species
          SELECT CASE( TRIM( ThisSpc%Name ) ) 
             CASE( 'ISOPN', 'MMN' )
                DryDepCount    = DryDepCount + 1
             CASE DEFAULT
                ! Nothing
          END SELECT
       ENDIF
    ELSE
       ThisSpc%Is_Drydep       = .FALSE.
       ThisSpc%DryDepID        = MISSING_INT
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
       ENDIF

    ELSE
       ThisSpc%Is_Wetdep    = .FALSE.
       ThisSpc%WetDepID     = MISSING_INT
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
    ! Is it a size-resolved aerosol?
    !---------------------------------------------------------------------
    IF ( PRESENT( WD_SizeResAer ) ) THEN
       ThisSpc%WD_SizeResAer = WD_SizeResAer
    ELSE
       ThisSpc%WD_SizeResAer = .FALSE.
    ENDIF

    !---------------------------------------------------------------------
    ! Conversion factors based on emitted molecular weight  
    ! These are mostly for backwards compatibility w/ existing code
    !
    ! NOTE: We have to use the emitted molecular weight EmMw_g instead
    ! of the actual molecular weight.  This will provide the proper
    ! unit conversion for species like ISOP which are emitted and
    ! transported as equivalent carbons.
    !
    ! These will be removed from GEOS-Chem once the new unit conversion
    ! routines are brought into the the code. (bmy, 9/2/15)
    !---------------------------------------------------------------------
    
    ! Ratio of MW of dry air per emitted MW Of species
    ThisSpc%TCVV   = AIRMW / ThisSpc%EmMw_g

    ! Molecules species per kg of species
    ThisSpc%XNUMOL = AVO   / ( ThisSpc%EmMw_g * 1e-3_fp ) 

    !---------------------------------------------------------------------
    ! Sanity checks
    !---------------------------------------------------------------------
    IF ( ThisSpc%Is_Gas ) THEN

       ! If this is a gas, then zero out all aerosol fields ...
       ThisSpc%WD_CoarseAer  = .FALSE.
       ThisSpc%WD_SizeResAer = .FALSE.
       
       ! ... except for those species that wetdep like aerosols
       SELECT CASE( TRIM( ThisSpc%Name ) )
          CASE( 'SO2', 'HNO3' )
             ! Do nothing
          CASE DEFAULT 
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
  SUBROUTINE Spc_Print( am_I_Root, ThisSpc, RC )
!
! !INPUT PARAMETERS:
! 
    LOGICAL,          INTENT(IN)    :: am_I_Root    ! Are we on the root CPU?
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
!EOP
!------------------------------------------------------------------------------
!BOC
 
    !=====================================================================
    ! Spc_Create begins here!
    !=====================================================================
    IF ( am_I_Root ) THEN

       !-----------------------
       ! Print general info
       !-----------------------
       WRITE( 6, '(a)' ) REPEAT( '=', 79 )  
       WRITE( 6, 100 ) 'Species ID           ',  ThisSpc%ModelID
       WRITE( 6, 110 ) 'Name                 ',  TRIM( ThisSpc%Name     )
       WRITE( 6, 110 ) 'FullName             ',  TRIM( ThisSpc%FullName )
       WRITE( 6, 120 ) 'Molecular weight [g] ',  ThisSpc%MW_g
       WRITE( 6, 120 ) 'Emitted mol. wt [g]  ',  ThisSpc%EmMW_g
       WRITE( 6, 120 ) 'Molecular ratio      ',  ThisSpc%MolecRatio
       WRITE( 6, 130 ) 'Is it a gas?         ',  ThisSpc%Is_Gas
       WRITE( 6, 130 ) 'Is it advected?      ',  ThisSpc%Is_Advected

       !-----------------------
       ! Print Henry's info
       !-----------------------
       IF ( ThisSpc%Henry_K0 > ZERO_R8 ) THEN
          WRITE( 6, 120 ) 'Henry''s law K0       ', ThisSpc%Henry_K0
       ENDIF

       IF ( ThisSpc%Henry_CR > ZERO_R8 ) THEN
          WRITE( 6, 120 ) 'Henry''s law CR       ', ThisSpc%Henry_CR
       ENDIF

       IF ( ThisSpc%Henry_pKa > ZERO_R8 ) THEN
          WRITE( 6, 120 ) 'Henry''s law pKa      ', ThisSpc%Henry_pKa
       ENDIF

       !-----------------------
       ! Print drydep info
       !-----------------------
       WRITE( 6, 130 ) 'Is it dry deposited? ',  ThisSpc%Is_DryDep

       IF ( ThisSpc%Is_DryDep ) THEN
          WRITE( 6, 100 ) ' -> Drydep index:    ',  ThisSpc%DryDepID

          IF ( ThisSpc%DD_A_Density > ZERO ) THEN
             WRITE( 6, 120 ) ' -> A_DENSITY        ',  ThisSpc%DD_A_Density
          ENDIF

          IF ( ThisSpc%DD_A_Radius > ZERO ) THEN
             WRITE( 6, 120 ) ' -> A_RADIUS         ',  ThisSpc%DD_A_Radius
          ENDIF

          IF ( ThisSpc%DD_F0 > ZERO ) THEN
             WRITE( 6, 120 ) ' -> F0 parameter     ',  ThisSpc%DD_F0
          ENDIF

          IF ( ThisSpc%DD_KOA > ZERO ) THEN
             WRITE( 6, 120 ) ' -> KOA parameter    ',  ThisSpc%DD_KOA
          ENDIF

          IF ( ThisSpc%DD_Hstar_old > ZERO ) THEN
             WRITE( 6, 120 ) ' -> Old HSTAR value  ',  ThisSpc%DD_Hstar_Old
          ENDIF
       ENDIF

       !-----------------------
       ! Print wetdep info
       !-----------------------
       WRITE( 6, 130 ) 'Is it wet deposited? ',  ThisSpc%Is_WetDep

       IF ( ThisSpc%Is_WetDep ) THEN
          WRITE( 6, 100 ) ' -> Wetdep index:    ',  ThisSpc%WetDepID

          IF ( ThisSpc%WD_RetFactor > ZERO ) THEN
             WRITE( 6, 120 ) ' -> Ret. Factor      ', ThisSpc%WD_RetFactor
          ENDIF

          IF ( ThisSpc%WD_AerScavEff > ZERO ) THEN
             WRITE( 6, 120 ) ' -> Scav. Effeciency ', ThisSpc%WD_AerScavEff
          ENDIF
       ENDIF

       ! Format statements
 100   FORMAT( a, ' : ', i5     )
 110   FORMAT( a, ' : ', a      )
 120   FORMAT( a, ' : ', es13.6 )
 130   FORMAT( a, ' : ', L1     )
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
  SUBROUTINE Spc_GetNumSpecies( nAdvect, nDryDep, nWetDep )
!
! !OUTPUT PARAMETERS:
!
    INTEGER, INTENT(OUT) :: nAdvect   ! # of advected species
    INTEGER, INTENT(OUT) :: nDryDep   ! # of dry-deposited species
    INTEGER, INTENT(OUT) :: nWetDep   ! # of wet-deposited species
! 
! !REVISION HISTORY: 
!   2 Sep 2015 - R. Yantosca - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC

    ! Return module variables
    nAdvect = AdvectCount
    nDryDep = DryDepCount
    nWetDep = WetDepCount
    
  END SUBROUTINE Spc_GetNumSpecies
!EOC
END MODULE Species_Mod

