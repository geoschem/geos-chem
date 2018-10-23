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
  PUBLIC  :: Spc_Info

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

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%% Uncomment this if you want to use the new henry's law constants
!%%% compiled by Katie Travis and posted on this wiki page:
!%%% http://wiki.geos-chem.org/Physical_properties_of_GEOS-Chem_species
!#define NEW_HENRY_CONSTANTS 1
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! !REMARKS:
!
! !REVISION HISTORY:
!  28 Aug 2015 - R. Yantosca - Initial version
!  02 Aug 2016 - M. Sulprizio- Add KppSpcId to store all KPP species incices.
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
! !DESCRIPTION: Initializes the GEOS-Chem species database object.  You can
!  add information about new species to this routine.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Init_Species_Database( am_I_Root, Input_Opt, SpcData, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE Input_Opt_Mod,       ONLY : OptInput
    USE Species_Mod
!
! !INPUT PARAMETERS:
!
    LOGICAL,        INTENT(IN)  :: am_I_Root    ! Are we on the root CPU?
    TYPE(OptInput), INTENT(IN)  :: Input_Opt    ! Input Options object
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(SpcPtr),   POINTER     :: SpcData(:)   ! Vector with species info
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT) :: RC           ! Success or failure?
!
! !REMARKS:
!  For detailed information about the species database (and the physical
!  properties that are specified there), please see the following GEOS-Chem
!  wiki pages:
!
!    (1) wiki.geos-chem.org/GEOS_Chem_species_database
!    (2) wiki.geos-chem.org/Physical_properties_of_GEOS-Chem_species
!
!  References for the new Henry's law constants:
!    (1) Sander et al [2015]: http://henrys-law.org
!
!  The Hcp (aka K0) parameter listed on the wiki page:
!
!     http://wiki.geos-chem.org/Physical_properties_of_GEOS-Chem_species
!
!  have units of [mol m-3 Pa-1].  But the GEOS-Chem species object expects
!  this parameter to have units of [M atm-1].  Therefore, when we pass the
!  Hcp paramter to routine Spc_Create via the HENRY_K0 argument, we will
!  multiply by the proper conversion factor (9.86923e-3) to convert
!  [mol m-3 Pa-1] to [M atm-1].
!
! !REVISION HISTORY:
!  22 Jul 2015 - R. Yantosca - Initial version
!  01 Sep 2015 - R. Yantosca - Add Henry K0, CR constants for DMS, ACET
!  02 Sep 2015 - R. Yantosca - Corrected typos for some SOA species
!  24 Sep 2015 - R. Yantosca - WD_RainoutEff is now a 3-element vector
!  24 Sep 2015 - R. Yantosca - Add WD_KcScaleFAc, a 3-element vector
!  30 Sep 2015 - R. Yantosca - DD_A_Density is renamed to Density
!  30 Sep 2015 - R. Yantosca - DD_A_Radius is renamed to Radius
!  01 Oct 2015 - R. Yantosca - Added DD_DvzMinVal field to put a minimum
!                              deposition velocity for sulfate aerosols
!  14 Oct 2015 - E. Lundgren - Treat H2SO4 as an aerosol for TOMAS
!  18 Nov 2015 - M. Sulprizio- Add passive tracers PASV to RnPbBe simulation
!  16 Dec 2015 - R. Yantosca - Use MW_g = 31.4 g/mol for SO4s and NITs
!  15 Mar 2016 - R. Yantosca - Added tagged CO tracer names
!  22 Apr 2016 - R. Yantosca - Now define Is_Hg0, Is_Hg2, Is_HgP fields
!  19 May 2016 - R. Yantosca - Remove DryDepId_PAN and DryDepId_HNO3; we shall
!                              now explicitly compute a drydep velocity for
!                              all GEOS-Chem species.
!  21 Jun 2016 - M. Sulprizio- Set Is_Photolysis to T for all species included
!                              in FAST-JX photolysis. Also added new species
!                              that are in FAST-JX photolysis but not already
!                              defined here.
!  18 Jul 2016 - M. Sulprizio- Remove family tracers ISOPN, MMN, CFCX, HCFCX
!                              and replace with their constituents.
!  02 Aug 2016 - M. Sulprizio- Add KppSpcId as argument passed to Spc_Create
!  11 Aug 2016 - E. Lundgren - Define special background conc for some species
!  22 Nov 2016 - M. Sulprizio- Move aerosol densities for BC, OC, and SO4 here
!                              from aerosol_mod.F
!  23 Feb 2017 - M. Sulprizio- Change MolecRatio for ALK4 from 4 to 4.3
!                              (B. Henderson)
!  13 Jun 2017 - M. Sulprizio- Add species for mechanistic isoprene SOA
!                              (E. Marais)
!  27 Nov 2017 - E. Lundgren - Add SALA, SALC, OCPO/OCPI, BCPO/BCPI, and SO4
!                              as hygroscopic growth species for cloud diags
!  14 Sep 2018 - C. Keller   - Now get species info from Spc_Info.
!  23 Oct 2018 - R. Yantosca - Cosmetic changes
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER             :: C,      N,   P,       nSpecies
    REAL(fp)            :: Radius, KOA, MW_g,    BackgroundVV, HStar
    REAL(f8)            :: K0,     CR

    ! Strings
    CHARACTER(LEN=31)   :: NameAllCaps
    CHARACTER(LEN=31)   :: Name
    CHARACTER(LEN=80)   :: FullName
    CHARACTER(LEN=80)   :: Formula

    ! For Tagged Hg species
    INTEGER             :: Hg0_CAT
    INTEGER             :: Hg2_CAT
    INTEGER             :: HgP_CAT

    ! For values from Input_Opt
    LOGICAL             :: Is_Advected
    LOGICAL             :: prtDebug

    ! For passive species
!    LOGICAL             :: IsPassive

    ! Species information
    REAL(fp)               :: EmMW_g           ! Emissions mol. wt [g]
    REAL(fp)               :: MolecRatio       ! Molec ratio
    REAL(f8)               :: Henry_K0         ! Henry's law K0 [M/atm]
    REAL(f8)               :: Henry_CR         ! Henry's law CR [K]
    REAL(f8)               :: Henry_PKA        ! Henry's law pKa [1]
    REAL(fp)               :: Density          ! Density [kg/m3]
    LOGICAL                :: DD_AeroDryDep    ! Use AERO_SFCRSII?
    LOGICAL                :: DD_DustDryDep    ! Use DUST_SFCRSII?
    REAL(fp)               :: DD_DvzAerSnow    ! Vd for aerosols
    REAL(fp)               :: DD_DvzMinVal(2)  ! Min Vd for aerosols
    REAL(fp)               :: DD_F0            ! Drydep reactivity [1]
    REAL(fp)               :: DD_KOA           ! Drydep KOA parameter
    REAL(fp)               :: DD_Hstar_Old     ! HSTAR, drydep [M/atm]
    LOGICAL                :: MP_SizeResAer    ! Size resolved aerosol?
    LOGICAL                :: MP_SizeResNum    ! Size resolved aer #?
    REAL(fp)               :: WD_RetFactor     ! Wetdep retention factor
    LOGICAL                :: WD_LiqAndGas     ! Liquid and gas phases?
    REAL(fp)               :: WD_ConvFacI2G    ! Factor for ice/gas ratio
    REAL(fp)               :: WD_AerScavEff    ! Aerosol scavenging eff.
    REAL(fp)               :: WD_KcScaleFac(3) ! Factor to multiply Kc
    REAL(fp)               :: WD_RainoutEff(3) ! Rainout efficiency
    LOGICAL                :: WD_CoarseAer     ! Coarse aerosol?
    LOGICAL                :: Is_Drydep        ! Is it dry deposited?
    LOGICAL                :: Is_Gas           ! Gas (T) or aerosol (F)?
    LOGICAL                :: Is_HygroGrowth   ! Is hygroscopic growth?
    LOGICAL                :: Is_Photolysis    ! Is it photolysis spc?
    LOGICAL                :: Is_Wetdep        ! Is it wet deposited?
    LOGICAL                :: Is_InRestart     ! Is it in restart file?
    LOGICAL                :: Is_Hg0           ! Denotes Hg0 species
    LOGICAL                :: Is_Hg2           ! Denotes Hg2 species
    LOGICAL                :: Is_HgP           ! Denotes HgP species
!
! !DEFINED PARAMETERS
!

    !=======================================================================
    ! Init_Species_Database begins here!
    !=======================================================================

    ! Initialize
    Hg0_CAT       = 0
    Hg2_CAT       = 0
    HgP_CAT       = 0

    ! Copy values from Input_Opt
    prtDebug      = ( Input_Opt%LPRT .and. am_I_Root )

    ! Store the list unique GEOS-Chem species names in work arrays for use
    ! below.  This is the combined list of advected species (from input.geos)
    ! plus KPP species (from SPC_NAMES in gckpp_Monitor.F90), with all
    ! duplicates removed.  Also stores the corresponding indices in the
    ! KPP VAR and FIX arrays.  For simulations that do not use KPP, the
    ! unique species list is the list of advected species from input.geos.
    CALL Unique_Species_Names( am_I_Root, Input_Opt, nSpecies, RC )

    ! Initialize the species vector
    CALL SpcData_Init( am_I_Root, nSpecies, SpcData, RC )
    IF ( RC /= GC_SUCCESS ) THEN
       PRINT*, '### Could not initialize species vector!'
       CALL EXIT( -999 )
    ENDIF

    ! Loop over all species
    DO N = 1, nSpecies

       ! It's advected if it's in the advected species list in input.geos
       Is_Advected = ANY( Input_Opt%AdvectSpc_Name == Species_Names(N) )

       ! Translate species name to uppercase
       NameAllCaps = TRIM( Species_Names(N) )
       CALL TranUc( NameAllCaps )

       ! Get species info, now write into local variable
       CALL Spc_Info( am_I_Root       = am_I_Root,                           &
                      iName           = TRIM(NameAllCaps),                   &
                      Input_Opt       = Input_Opt,                           &
                      KppSpcId        = KppSpcId(N),                         &
                      oFullName       = FullName,                            &
                      oFormula        = Formula,                             &
                      oMW_g           = MW_g,                                & 
                      oEmMW_g         = EmMW_g,                              &
                      oMolecRatio     = MolecRatio,                          &
                      oBackgroundVV   = BackgroundVV,                        &
                      oHenry_K0       = Henry_K0,                            &
                      oHenry_CR       = Henry_CR,                            &
                      oHenry_PKA      = Henry_PKA,                           &
                      oDensity        = Density,                             &
                      oRadius         = Radius,                              &
                      oDD_AeroDryDep  = DD_AeroDryDep,                       &
                      oDD_DustDryDep  = DD_DustDryDep,                       &
                      oDD_DvzAerSnow  = DD_DvzAerSnow,                       &
                      oDD_DvzMinVal   = DD_DvzMinVal,                        &
                      oDD_F0          = DD_F0,                               &
                      oDD_KOA         = DD_KOA,                              &
                      oDD_Hstar_Old   = DD_Hstar_Old,                        &
                      oMP_SizeResAer  = MP_SizeResAer,                       &
                      oMP_SizeResNum  = MP_SizeResNum,                       &
                      oWD_RetFactor   = WD_RetFactor,                        &
                      oWD_LiqAndGas   = WD_LiqAndGas,                        &
                      oWD_ConvFacI2G  = WD_ConvFacI2G,                       &
                      oWD_AerScavEff  = WD_AerScavEff,                       &
                      oWD_KcScaleFac  = WD_KcScaleFac,                       &
                      oWD_RainoutEff  = WD_RainoutEff,                       &
                      oWD_CoarseAer   = WD_CoarseAer,                        &
                      oIs_Drydep      = Is_Drydep,                           & 
                      oIs_Gas         = Is_Gas,                              &
                      oIs_HygroGrowth = Is_HygroGrowth,                      &
                      oIs_Photolysis  = Is_Photolysis,                       &
                      oIs_Wetdep      = Is_Wetdep,                           &
                      oIs_InRestart   = Is_InRestart,                        &
                      oIs_Hg0         = Is_Hg0,                              &
                      oIs_Hg2         = Is_Hg2,                              &  
                      oIs_HgP         = Is_HgP,                              &
                      RC              = RC                                  )  
       ! Error
       IF ( RC /= GC_SUCCESS ) THEN
          PRINT*, '### Could not get species info!',TRIM(NameAllCaps)
          CALL EXIT( -999 )
       ENDIF

       ! Create species
       IF ( EmMW_g == MISSING_MW ) THEN
          EmMW_g = MW_g
       ENDIF
       Name = TRIM(Species_Names(N))
       ! Handle exceptions
       IF ( TRIM(Name) == 'POx' ) Name = 'POX'
       IF ( TRIM(Name) == 'LOx' ) Name = 'LOX'
       IF ( TRIM(FullName) == '' ) FullName = TRIM(Name) 
       CALL Spc_Create( am_I_Root      = am_I_Root,                          &
                        ThisSpc        = SpcData(N)%Info,                    &
                        ModelID        = N,                                  &
                        KppSpcId       = KppSpcId(N),                        & 
                        KppVarId       = KppVarId(N),                        & 
                        KppFixId       = KppFixId(N),                        &
                        Name           = TRIM(Name),                         & 
                        FullName       = FullName,                           &
                        Formula        = Formula,                            &
                        MW_g           = MW_g,                               & 
                        EmMW_g         = EmMW_g,                             &
                        MolecRatio     = MolecRatio,                         &
                        BackgroundVV   = BackgroundVV,                       &
                        Henry_K0       = Henry_K0,                           &
                        Henry_CR       = Henry_CR,                           &
                        Henry_PKA      = Henry_PKA,                          &
                        Density        = Density,                            &
                        Radius         = Radius,                             &
                        DD_AeroDryDep  = DD_AeroDryDep,                      & 
                        DD_DustDryDep  = DD_DustDryDep,                      & 
                        DD_DvzAerSnow  = DD_DvzAerSnow,                      &
                        DD_DvzMinVal   = DD_DvzMinVal,                       &
                        DD_F0          = DD_F0,                              &
                        DD_KOA         = DD_KOA,                             &
                        DD_Hstar_Old   = DD_Hstar_Old,                       &
                        MP_SizeResAer  = MP_SizeResAer,                      &
                        MP_SizeResNum  = MP_SizeResNum,                      &
                        WD_RetFactor   = WD_RetFactor,                       &
                        WD_LiqAndGas   = WD_LiqAndGas,                       &
                        WD_ConvFacI2G  = WD_ConvFacI2G,                      &
                        WD_AerScavEff  = WD_AerScavEff,                      &
                        WD_KcScaleFac  = WD_KcScaleFac,                      &
                        WD_RainoutEff  = WD_RainoutEff,                      &
                        WD_CoarseAer   = WD_CoarseAer,                       &
                        Is_Advected    = Is_Advected,                        & 
                        Is_Drydep      = Is_Drydep,                          & 
                        Is_Gas         = Is_Gas,                             &
                        Is_HygroGrowth = Is_HygroGrowth,                     &
                        Is_Photolysis  = Is_Photolysis,                      &
                        Is_Wetdep      = Is_Wetdep,                          &
                        Is_InRestart   = Is_InRestart,                       &
                        Is_Hg0         = Is_Hg0,                             &
                        Is_Hg2         = Is_Hg2,                             &  
                        Is_HgP         = Is_HgP,                             &
                        RC             = RC                                 )  
       ! Error
       IF ( RC /= GC_SUCCESS ) THEN
          PRINT*, '### Could not initialize species vector!',TRIM(NameAllCaps)
          CALL EXIT( -999 )
       ENDIF

       ! Print info about each species
       ! testing only
       CALL Spc_Print( am_I_Root, SpcData(N)%Info, RC )
       IF ( prtDebug ) CALL Spc_Print( am_I_Root, SpcData(N)%Info, RC )

    ENDDO

    ! Deallocate temporary work arrays
    CALL Cleanup_Work_Arrays()

  END SUBROUTINE Init_Species_Database
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Spc_Info
!
! !DESCRIPTION: Routine Spc\_Info is a helper function that returns the
!  properties of a species. 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Spc_Info( am_I_Root,       iName,          Input_Opt,           &
                       KppSpcId,        oFullName,      oFormula,            &
                       oMW_g,           oEmMW_g,        oMolecRatio,         &
                       oBackgroundVV,   oHenry_K0,      oHenry_CR,           &
                       oHenry_PKA,      oDensity,       oRadius,             &
                       oDD_AeroDryDep,  oDD_DustDryDep, oDD_DvzAerSnow,      &
                       oDD_DvzMinVal,   oDD_F0,         oDD_KOA,             &
                       oDD_HStar_Old,   oMP_SizeResAer, oMP_SizeResNum,      &
                       oWD_RetFactor,   oWD_LiqAndGas,  oWD_ConvFacI2G,      &
                       oWD_AerScavEff,  oWD_KcScaleFac, oWD_RainoutEff,      &
                       oWD_CoarseAer,                                        &
                       oIs_Drydep,      oIs_Gas,        oIs_HygroGrowth,     &
                       oIs_Photolysis,  oIs_Wetdep,     oIs_InRestart,       &
                       oIs_Hg0,         oIs_Hg2,        oIs_HgP,             &
                       Found,           Underscores,    RC                  )
!
! !USES:
!
    USE Input_Opt_Mod,           ONLY : OptInput
    USE SPECIES_MOD,             ONLY : MISSING_INT,   MISSING,    MISSING_R8
    USE SPECIES_MOD,             ONLY : ZERO, ZERO_R8, MISSING_MW, MISSING_VV
!
! !INPUT PARAMETERS:
! 
    LOGICAL,           INTENT(IN)  :: am_I_Root        ! Are we on the root CPU?
    CHARACTER(LEN=*),  INTENT(IN)  :: iName            ! Short name of species
    TYPE(OptInput),    OPTIONAL    :: Input_Opt        ! Input Options object
    INTEGER,           INTENT(IN)  :: KppSpcId         ! KPP ID 
!
! !OUTPUT PARAMETERS:
! 
    CHARACTER(LEN=*), INTENT(OUT), OPTIONAL   :: oFullName         ! Long name of species
    CHARACTER(LEN=*), INTENT(OUT), OPTIONAL   :: oFormula          ! Chemical formula
    REAL(fp),INTENT(OUT), OPTIONAL    :: oMW_g             ! Molecular weight [g]
    REAL(fp),INTENT(OUT), OPTIONAL    :: oEmMW_g           ! Emissions mol. wt [g]
    REAL(fp),INTENT(OUT), OPTIONAL    :: oMolecRatio       ! Molec ratio
    REAL(fp),INTENT(OUT), OPTIONAL    :: oBackgroundVV     ! Background conc [v/v]
    REAL(f8),INTENT(OUT), OPTIONAL    :: oHenry_K0         ! Henry's law K0 [M/atm]
    REAL(f8),INTENT(OUT), OPTIONAL    :: oHenry_CR         ! Henry's law CR [K]
    REAL(f8),INTENT(OUT), OPTIONAL    :: oHenry_PKA        ! Henry's law pKa [1]
    REAL(fp),INTENT(OUT), OPTIONAL    :: oDensity          ! Density [kg/m3]
    REAL(fp),INTENT(OUT), OPTIONAL    :: oRadius           ! Radius [m]
    LOGICAL, INTENT(OUT), OPTIONAL    :: oDD_AeroDryDep    ! Use AERO_SFCRSII?
    LOGICAL, INTENT(OUT), OPTIONAL    :: oDD_DustDryDep    ! Use DUST_SFCRSII?
    REAL(fp),INTENT(OUT), OPTIONAL    :: oDD_DvzAerSnow    ! Vd for aerosols
                                                      !  on snow/ice [cm/s]
    REAL(fp),INTENT(OUT), OPTIONAL    :: oDD_DvzMinVal(2)  ! Min Vd for aerosols
                                                      !  (cf GOCART) [cm/s]
    REAL(fp),INTENT(OUT), OPTIONAL    :: oDD_F0            ! Drydep reactivity [1]
    REAL(fp),INTENT(OUT), OPTIONAL    :: oDD_KOA           ! Drydep KOA parameter
    REAL(fp),INTENT(OUT), OPTIONAL    :: oDD_Hstar_Old     ! HSTAR, drydep [M/atm]
    LOGICAL, INTENT(OUT), OPTIONAL    :: oMP_SizeResAer    ! Size resolved aerosol?
    LOGICAL, INTENT(OUT), OPTIONAL    :: oMP_SizeResNum    ! Size resolved aer #?
    REAL(fp),INTENT(OUT), OPTIONAL    :: oWD_RetFactor     ! Wetdep retention factor
    LOGICAL, INTENT(OUT), OPTIONAL    :: oWD_LiqAndGas     ! Liquid and gas phases?
    REAL(fp),INTENT(OUT), OPTIONAL    :: oWD_ConvFacI2G    ! Factor for ice/gas ratio
    REAL(fp),INTENT(OUT), OPTIONAL    :: oWD_AerScavEff    ! Aerosol scavenging eff.
    REAL(fp),INTENT(OUT), OPTIONAL    :: oWD_KcScaleFac(3) ! Factor to multiply Kc
                                                           !  rate in F_AEROSOL
    REAL(fp),INTENT(OUT), OPTIONAL    :: oWD_RainoutEff(3) ! Rainout efficiency
    LOGICAL, INTENT(OUT), OPTIONAL    :: oWD_CoarseAer     ! Coarse aerosol?
    LOGICAL, INTENT(OUT), OPTIONAL    :: oIs_Drydep        ! Is it dry deposited?
    LOGICAL, INTENT(OUT), OPTIONAL    :: oIs_Gas           ! Gas (T) or aerosol (F)?
    LOGICAL, INTENT(OUT), OPTIONAL    :: oIs_HygroGrowth   ! Is hygroscopic growth?
    LOGICAL, INTENT(OUT), OPTIONAL    :: oIs_Photolysis    ! Is it photolysis spc?
    LOGICAL, INTENT(OUT), OPTIONAL    :: oIs_Wetdep        ! Is it wet deposited?
    LOGICAL, INTENT(OUT), OPTIONAL    :: oIs_InRestart     ! Is it in restart file?
    LOGICAL, INTENT(OUT), OPTIONAL    :: oIs_Hg0           ! Denotes Hg0 species
    LOGICAL, INTENT(OUT), OPTIONAL    :: oIs_Hg2           ! Denotes Hg2 species
    LOGICAL, INTENT(OUT), OPTIONAL    :: oIs_HgP           ! Denotes HgP species
    INTEGER, INTENT(OUT)              :: RC                ! Return code
    LOGICAL, INTENT(OUT), OPTIONAL    :: Found             ! Species found? If arg present, 
                                                           ! no error if not found 
    LOGICAL, INTENT(IN),  OPTIONAL    :: Underscores       ! Replace blanks with underscores 
! 
! !REMARKS:
!
! !REVISION HISTORY: 
!  14 Sep 2018 - C. Keller   - Created standalone subroutine so that species
!                              info can be queried independently.
!  23 Oct 2018 - R. Yantosca - Cosmetic changes (consistent indentation)
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! local species information
    CHARACTER(LEN=31)   :: Name
    CHARACTER(LEN=80)   :: FullName
    CHARACTER(LEN=80)   :: Formula
    REAL(fp)            :: MW_g             ! Molecular weight [g]
    REAL(fp)            :: EmMW_g           ! Emissions mol. wt [g]
    REAL(fp)            :: MolecRatio       ! Molec ratio
    REAL(fp)            :: BackgroundVV     ! Background conc [v/v]
    REAL(f8)            :: Henry_K0         ! Henry's law K0 [M/atm]
    REAL(f8)            :: Henry_CR         ! Henry's law CR [K]
    REAL(f8)            :: Henry_PKA        ! Henry's law pKa [1]
    REAL(fp)            :: Density          ! Density [kg/m3]
    REAL(fp)            :: Radius           ! Radius [m]
    LOGICAL             :: DD_AeroDryDep    ! Use AERO_SFCRSII?
    LOGICAL             :: DD_DustDryDep    ! Use DUST_SFCRSII?
    REAL(fp)            :: DD_DvzAerSnow    ! Vd for aerosols
    REAL(fp)            :: DD_DvzMinVal(2)  ! Min Vd for aerosols
    REAL(fp)            :: DD_F0            ! Drydep reactivity [1]
    REAL(fp)            :: DD_KOA           ! Drydep KOA parameter
    REAL(fp)            :: DD_Hstar_Old     ! HSTAR, drydep [M/atm]
    LOGICAL             :: MP_SizeResAer    ! Size resolved aerosol?
    LOGICAL             :: MP_SizeResNum    ! Size resolved aer #?
    REAL(fp)            :: WD_RetFactor     ! Wetdep retention factor
    LOGICAL             :: WD_LiqAndGas     ! Liquid and gas phases?
    REAL(fp)            :: WD_ConvFacI2G    ! Factor for ice/gas ratio
    REAL(fp)            :: WD_AerScavEff    ! Aerosol scavenging eff.
    REAL(fp)            :: WD_KcScaleFac(3) ! Factor to multiply Kc
    REAL(fp)            :: WD_RainoutEff(3) ! Rainout efficiency
    LOGICAL             :: WD_CoarseAer     ! Coarse aerosol?
    LOGICAL             :: Is_Advected      ! Is it advected?
    LOGICAL             :: Is_Drydep        ! Is it dry deposited?
    LOGICAL             :: Is_Gas           ! Gas (T) or aerosol (F)?
    LOGICAL             :: Is_HygroGrowth   ! Is hygroscopic growth?
    LOGICAL             :: Is_Photolysis    ! Is it photolysis spc?
    LOGICAL             :: Is_Wetdep        ! Is it wet deposited?
    LOGICAL             :: Is_InRestart     ! Is it in restart file?
    LOGICAL             :: Is_Hg0           ! Denotes Hg0 species
    LOGICAL             :: Is_Hg2           ! Denotes Hg2 species
    LOGICAL             :: Is_HgP           ! Denotes HgP species

    ! Scalars
    REAL(f8)            :: K0,  CR
    REAL(fp)            :: KOA, HStar
    INTEGER             :: P, IDX, CNT
    LOGICAL             :: IsPassive
    LOGICAL             :: Uscore    

    ! Arrays
    REAL(fp)            :: DvzMinVal(2)
    REAL(fp)            :: KcScale(3)
    REAL(fp)            :: RainEff(3)
!
! !DEFINED PARAMETERS:
!
    ! Local parameter
    LOGICAL,  PARAMETER :: T        = .TRUE.         ! Yes
    LOGICAL,  PARAMETER :: F        = .FALSE.        ! No
    REAL(f8), PARAMETER :: To_M_atm = 9.86923e-3_f8  ! mol/m3/Pa -> M/atm

    !=====================================================================
    ! Spc_Info begins here!
    !=====================================================================

    ! Default values
    FullName         = ''
    Formula          = ''
    MW_g             = MISSING_MW
    EmMW_g           = MISSING_MW
    MolecRatio       = 1.0e+0_fp
    Radius           = MISSING
    Density          = MISSING
    DD_AeroDryDep    = .FALSE.
    DD_DustDryDep    = .FALSE.
    DD_F0            = MISSING
    DD_DvzAerSnow    = MISSING
    DD_DvzMinVal(:)  = MISSING
    DD_KOA           = MISSING
    DD_Hstar_Old     = MISSING
    Henry_K0         = MISSING_R8 
    Henry_CR         = MISSING_R8
    Henry_PKA        = MISSING_R8
    WD_RetFactor     = MISSING
    WD_LiqAndGas     = .FALSE.
    WD_ConvFacI2G    = MISSING
    WD_AerScavEff    = MISSING
    WD_KcScaleFac(:) = MISSING
    WD_RainOutEff(:) = MISSING
    Is_Drydep        = .FALSE.
    Is_Gas           = .FALSE.
    Is_HygroGrowth   = .FALSE.
    Is_Photolysis    = .FALSE.
    Is_Wetdep        = .FALSE.
    BackgroundVV     = MISSING_VV
    Is_InRestart     = .FALSE.
    WD_CoarseAer     = .FALSE.
    MP_SizeResAer    = .FALSE.
    MP_SizeResNum    = .FALSE.
    Is_Hg0           = .FALSE.
    Is_Hg2           = .FALSE.
    Is_HgP           = .FALSE.

    Uscore = .FALSE.
    IF ( PRESENT( Underscores ) ) Uscore = Underscores
    IF ( PRESENT( Found       ) ) Found  = .TRUE.

    ! Make sure its all caps
    Name = TRIM( iName )
    CALL TranUc( Name )

    ! Test for species name
    SELECT CASE( TRIM ( Name ) )

          !==================================================================
          ! Species for the various "full-chemistry" mechanisms
          !
          ! (1) benchmark  (2) SOA         (3) UCX
          ! (4) tropchem   (5) aciduptake  (6) marine OC
          !
          ! CH4 is also contained here, as it is part of the benchmark
          ! and UCX mechanisms.
          !==================================================================

          CASE( 'ACET' )
             FullName      = 'Acetone'
             Formula       = 'CH3C(O)CH3'
             MW_g          = 58.08_fp
             EmMW_g        = 12.00_fp
             MolecRatio    = 3.0_fp
             Is_Advected   = Is_Advected
             Is_Gas        = T
             Is_Drydep     = T
             Is_Wetdep     = F
             Is_Photolysis = T
             DD_F0         = 1.0_fp
#if defined( NEW_HENRY_CONSTANTS )
             Henry_K0      = 2.7e-1_f8 * To_M_atm
             Henry_CR      = 5500.0_f8
#else
             DD_Hstar_Old  = 1e5_fp
             Henry_K0      = 2.7e+1_f8
             Henry_CR      = 5300.0_f8
#endif

          CASE( 'ACTA' )
             FullName      = 'Acetic acid'
             Formula       = 'CH3C(O)OH'
             MW_g          = 60.00_fp
             Is_Gas        = T
             Is_Drydep     = T
             Is_Wetdep     = T
             DD_F0         = 1.0_fp
             DD_Hstar_Old  = 4.1e+3_fp
             Henry_K0      = 4.1e+3_f8
             Henry_CR      = 6300.0_f8
             WD_RetFactor  = 2.0e-2_fp


          CASE( 'ALD2' )
             FullName      = 'Acetaldehyde'
             Formula       = 'CH3CHO'
             MW_g          = 44.05_fp
             EmMW_g        = 12.0_fp
             MolecRatio    = 2.0_fp
             Is_Gas        = T
             Is_Drydep     = T
             Is_Wetdep     = T
             Is_Photolysis = T
             DD_F0         = 1.0_fp
#if defined( NEW_HENRY_CONSTANTS )
             Henry_K0      = 1.30e-01_f8 * To_M_atm
             Henry_CR      = 5900.0_f8
#else
             DD_Hstar_Old  = 1.1e+1_fp
             Henry_K0      = 1.1e+1_f8
             Henry_CR      = 6300.0_f8
#endif
             WD_RetFactor  = 2.0e-2_fp

          CASE( 'ALK4' )
             FullName      = 'Lumped >= C4 Alkanes'
             Formula       = ''
             MW_g          = 58.12_fp
             EmMW_g        = 12.0_fp
             MolecRatio    = 4.3_fp
             Is_Gas        = T
             Is_Drydep     = F
             Is_Wetdep     = F
#if defined( NEW_HENRY_CONSTANTS )
             Henry_K0      = 1.20e-5_f8 * To_M_atm
             Henry_CR      = 3100.0_f8
#endif

          CASE( 'ASOA1', 'ASOA2', 'ASOA3', 'ASOAN' )
             FullName = 'Lumped non-volatile aerosol products of light aromatics + IVOCs'

             ! Halve the Kc (cloud condensate -> precip) rate
             ! for the temperature range 237 K <= T < 258 K.
             KcScale = (/ 1.0_fp, 0.5_fp, 1.0_fp /)

             ! Turn off rainout only when 237 K <= T < 258K.
             ! NOTE: Rainout efficiency is 0.8 because these are SOA species.
             RainEff = (/ 0.8_fp, 0.0_fp, 0.8_fp /)

             Formula       = ''
             MW_g          = 150.0_fp
             Is_Gas        = F
             Is_Drydep     = T
             Is_Wetdep     = T
             DD_DvzAerSnow = 0.03_fp
             DD_F0         = 0.0_fp
             DD_Hstar_Old  = 0.0_fp
             WD_AerScavEff = 0.8_fp
             WD_KcScaleFac = KcScale
             WD_RainoutEff = RainEff

          CASE( 'ASOG1', 'ASOG2', 'ASOG3' )
             FullName = 'Lumped non-volatile gas products of light aromatics + IVOCs'
             Formula       = ''
             MW_g          = 150.0_fp
             Is_Gas        = T
             Is_Drydep     = T
             Is_Wetdep     = T
             DD_F0         = 0.0_fp
#if defined( NEW_HENRY_CONSTANTS )
             Henry_K0      = 1.00e+5_f8
             Henry_CR      = 6039.0_f8
#else
                              DD_Hstar_Old  = 1.00e+5_fp
                              Henry_K0      = 1.00e+5_f8
                              Henry_CR      = 6039.0_f8
#endif
                              WD_RetFactor  = 2.0e-2_fp

          CASE( 'BCPI' )
             FullName = 'Hydrophilic black carbon aerosol'

             ! Halve the Kc (cloud condensate -> precip) rate
             ! for the temperature range 237 K <= T < 258 K.
             KcScale  = (/ 1.0_fp, 0.5_fp, 1.0_fp /)

             ! Turn off rainout only when 237 K <= T < 258K.
             RainEff  = (/ 1.0_fp, 0.0_fp, 1.0_fp /)

             Formula       = ''
             MW_g          = 12.01_fp
             EmMW_g        = 12.0_fp
             Is_Gas        = F
             Is_Drydep     = T
             Is_Wetdep     = T
             Is_HygroGrowth= T
             Density       = 1800.0_fp
             DD_DvzAerSnow = 0.03_fp
             DD_F0         = 0.0_fp
             DD_Hstar_Old  = 0.0_fp
             WD_AerScavEff = 1.0_fp
             WD_KcScaleFac = KcScale
             WD_RainoutEff = RainEff

          CASE( 'BCPO' )
             Fullname = 'Hydrophobic black carbon aerosol'

             ! Halve the Kc (cloud condensate -> precip) rate
             ! for the temperature range T > 258 K
             KcScale  = (/ 1.0_fp, 1.0_fp, 0.5_fp /)

             ! Allow rainout of BCPO when T < 258 K, because
             ! BCPO is considered to be IN.
             RainEff  = (/ 1.0_fp, 1.0_fp, 0.0_fp /)

             Formula       = ''
             MW_g          = 12.01_fp
             EmMW_g        = 12.0_fp
             Is_Gas        = F
             Is_Drydep     = T
             Is_Wetdep     = T
             Is_HygroGrowth= F
             Density       = 1800.0_fp
             DD_DvzAerSnow = 0.03_fp
             DD_F0         = 0.0_fp
             DD_Hstar_Old  = 0.0_fp
             WD_AerScavEff = 1.0_fp
             WD_KcScaleFac = KcScale
             WD_RainoutEff = RainEff

          CASE( 'BENZ' )
             FullName      = 'Benzene'
             Formula       = 'C6H6'
             MW_g          = 78.11_fp
             EmMW_g        = 12.0_fp
             MolecRatio    = 6.0_fp
             Is_Gas        = T
             Is_Drydep     = F
             Is_Wetdep     = F
             
          CASE( 'BR' )
             FullName      = 'Atomic bromine'
             Formula       = 'Br'
             MW_g          = 80.0_fp
             Is_Gas        = T
             Is_Drydep     = F
             Is_Wetdep     = F
#if defined( NEW_HENRY_CONSTANTS )
             Henry_K0      = 3.40e-4_f8 * To_M_atm
             Henry_CR      = 1800.0_f8
#endif

          CASE( 'BR2' )
             FullName      = 'Molecular Bromine'
             Formula       = 'Br2'
             MW_g          = 160.0_fp
             Is_Gas        = T
             Is_Drydep     = T
             Is_Wetdep     = T
             Is_Photolysis = T
             DD_F0         = 0.0_fp
#if defined( NEW_HENRY_CONSTANTS )
             Henry_K0      = 7.20e-3_f8 * To_M_atm
             Henry_CR      = 4400.0_f8
#else
             DD_Hstar_Old  = 7.60e-1_fp
             Henry_K0      = 7.60e-1_f8
             Henry_CR      = 3720.0_f8
#endif
             WD_RetFactor  = 0.0_fp

          CASE( 'BRCL' )
             FullName      = 'Bromine chloride'
             Formula       = 'BrCl'
             MW_g          = 115.0_fp
             Is_Gas        = T
             Is_Drydep     = F
             Is_Wetdep     = F
             Is_Photolysis = T

          CASE( 'BRNO2' )
             FullName      = 'Nitryl bromide'
             Formula       = 'BrNO2'
             MW_g          = 126.0_fp
             Is_Gas        = T
             Is_Drydep     = F
             Is_Wetdep     = F
             Is_Photolysis = T

          CASE( 'BRNO3' )
             FullName      = 'Bromine nitrate'
             Formula       = 'BrNO3'
             MW_g          = 142.0_fp
             Is_Gas        = T
             Is_Drydep     = T
             Is_Wetdep     = F
             Is_Photolysis = T
             DD_F0         = 0.0_fp
             DD_Hstar_Old  = 1.00e+20_fp

          CASE( 'BRO' )
             FullName      = 'Bromine monoxide'
             Formula       = 'BrO'
             MW_g          = 96.0_fp
             Is_Gas        = T
             Is_Drydep     = F
             Is_Wetdep     = F
             Is_Photolysis = T

          CASE( 'C2H6' )
             FullName      = 'Ethane'
             Formula       = 'C2H6'
             MW_g          = 30.07_fp
             EmMW_g        = 12.0_fp
             MolecRatio    = 2.0_fp
             Is_Gas        = T
             Is_Drydep     = F
             Is_Wetdep     = F
             DD_F0         = 1.0_fp
#if defined( NEW_HENRY_CONSTANTS )
             Henry_K0      = 1.90e+5_f8 * To_M_atm
             Henry_CR      = 2400.0_f8
#endif

          CASE( 'C3H8' )
             FullName      = 'Propane'
             Formula       = 'C3H8'
             MW_g          = 44.1_fp
             EmMW_g        = 12.0_fp
             MolecRatio    = 3.0_fp
             Is_Gas        = T
             Is_Drydep     = F
             Is_Wetdep     = F
#if defined( NEW_HENRY_CONSTANTS )
             Henry_K0      = 1.50e-5_f8 * To_M_atm
             Henry_CR      = 2700.0_f8
#endif

          CASE( 'CCL4' )
             FullName      = 'Carbon tetrachloride'
             Formula       = 'CCl4'
             MW_g          = 152.0_fp
             Is_Gas        = T
             Is_Drydep     = F
             Is_Wetdep     = F
             Is_Photolysis = T

          CASE( 'CFC11' )
             FullName      = 'CFC-11'
             Formula       = 'CCl3F'
             MW_g          = 137.0_fp
             Is_Gas        = T
             Is_Drydep     = F
             Is_Wetdep     = F
             Is_Photolysis = T
             
          CASE( 'CFC12' )
             FullName      = 'CFC-12'
             Formula       = 'CCl2F2'
             MW_g          = 121.0_fp
             Is_Gas        = T
             Is_Drydep     = F
             Is_Wetdep     = F
             Is_Photolysis = T

          CASE( 'CFC113' )
             FullName      = 'CFC-113'
             Formula       = 'C2Cl3F3'
             MW_g          = 187.0_fp
             Is_Gas        = T
             Is_Drydep     = F
             Is_Wetdep     = F
             Is_Photolysis = T

          CASE( 'CFC114' )
             FullName      = 'CFC-114'
             Formula       = 'C2Cl2F4'
             MW_g          = 187.0_fp
             Is_Gas        = T
             Is_Drydep     = F
             Is_Wetdep     = F
             Is_Photolysis = T

          CASE( 'CFC115' )
             FullName      = 'CFC-115'
             Formula       = 'C2ClF5'
             MW_g          = 187.0_fp
             Is_Gas        = T
             Is_Drydep     = F
             Is_Wetdep     = F
             Is_Photolysis = T

          CASE( 'CH2BR2' )
             FullName      = 'Dibromomethane'
             Formula       = 'CH2Br2'
             MW_g          = 174.0_fp
             Is_Gas        = T
             Is_Drydep     = F
             Is_Wetdep     = F
             Is_Photolysis = T

          CASE( 'CH2O', 'HCHO' )
             FullName      = 'Formaldehyde'
             Formula       = 'CH2O'
             MW_g          = 30.0_fp
             BackgroundVV  = 4.0e-15_fp
             Is_Gas        = T
             Is_Drydep     = T
             Is_Wetdep     = T
             Is_Photolysis = T
             DD_F0         = 1.0_fp
#if defined( NEW_HENRY_CONSTANTS )
             Henry_K0      = 3.2e+1_f8 * To_M_atm
             Henry_CR      = 6800.0_f8
#else
             DD_Hstar_Old  = 3.0e+3_fp
             Henry_K0      = 3.0e+3_f8
             Henry_CR      = 7200.0_f8
#endif
             WD_RetFactor  = 2.0e-2_fp

          CASE( 'CH3BR' )
             FullName      = 'Methyl bromide'
             Formula       = 'CH3Br'
             MW_g          = 95.0_fp
             Is_Gas        = T
             Is_Drydep     = F
             Is_Wetdep     = F
             Is_Photolysis = T

          CASE( 'CHCL3' )
             FullName      = 'Chloroform'
             Formula       = 'CHCl3'
             MW_g          = 119.0_fp
             Is_Gas        = T
             Is_Drydep     = F
             Is_Wetdep     = F
             Is_Photolysis = T

          CASE( 'CH2CL2' )
             FullName      = 'Dichloromethane'
             Formula       = 'CH2Cl2'
             MW_g          = 85.0_fp
             Is_Gas        = T
             Is_Drydep     = F
             Is_Wetdep     = F
             Is_Photolysis = T

          CASE( 'CH3CL' )
             FullName      = 'Chloromethane'
             Formula       = 'CH3Cl'
             MW_g          = 50.0_fp
             Is_Gas        = T
             Is_Drydep     = F
             Is_Wetdep     = F
             Is_Photolysis = T

          CASE( 'CH3CCL3' )
             FullName      = 'Methyl chloroform'
             Formula       = 'CH3CCl3'
             MW_g          = 133.0_fp
             Is_Gas        = T
             Is_Drydep     = F
             Is_Wetdep     = F
             Is_Photolysis = T

          ! Now include both total and tagged CH4 species (mps, 6/21/17)
          ! All of these have identical properties except for the names
          CASE( 'CH4',     'CH4_OIL', 'CH4_GAS', 'CH4_COL', 'CH4_LIV', &
                'CH4_LDF', 'CH4_WST', 'CH4_RIC', 'CH4_OTA', 'CH4_BBN', &
                'CH4_WTL', 'CH4_SEE', 'CH4_LAK', 'CH4_TER', 'CH4_SAB' )

             SELECT CASE( TRIM( Name ) )
                CASE( 'CH4_OIL' )
                   FullName = 'CH4 from oil emissions'
                CASE( 'CH4_GAS' )
                   FullName = 'CH4 from gas emissions'
                CASE( 'CH4_COL' )
                   FullName = 'CH4 from coal mining emissions'
                CASE( 'CH4_LIV' )
                   FullName = 'CH4 from livestock emissions'
                CASE( 'CH4_LDF' )
                   FullName = 'CH4 from landfill emissions'
                CASE( 'CH4_WST' )
                   FullName = 'CH4 from waste emissions'
                CASE( 'CH4_RIC' )
                   FullName = 'CH4 from rice emissions'
                CASE( 'CH4_OTA' )
                   FullName = 'CH4 from other anthropogenic emissions'
                CASE( 'CH4_BBN' )
                   FullName = 'CH4 from biomass burning emissions'
                CASE( 'CH4_WTL' )
                   FullName = 'CH4 from wetland emissions'
                CASE( 'CH4_SEE' )
                   FullName = 'CH4 from geological seep emissions'
                CASE( 'CH4_LAK' )
                   FullName = 'CH4 from lake emissions'
                CASE( 'CH4_TER' )
                   FullName = 'CH4 from termite emissions'
                CASE( 'CH4_SAB' )
                   FullName = 'CH4 from soil absorption emissions'
                CASE DEFAULT
                   FullName = 'Methane'
             END SELECT

             SELECT CASE( TRIM( Name ) )
                CASE( 'CH4' )
                   BackgroundVV = 1.7e-06_fp
                CASE DEFAULT
                   BackgroundVV = MISSING_VV
             END SELECT

             Formula       = 'CH4'
             MW_g          = 16.0_fp
             Is_Gas        = T
             Is_Drydep     = F
             Is_Wetdep     = F

          CASE( 'CHBR3' )
             FullName      = 'Bromoform'
             Formula       = 'CHBr3'
             MW_g          = 253.0_fp
             Is_Gas        = T
             Is_Drydep     = F
             Is_Wetdep     = F
             Is_Photolysis = T

          CASE( 'CL' )
             FullName      = 'Atomic chlorine'
             Formula       = 'Cl'
             MW_g          = 35.0_fp
             Is_Gas        = T
             Is_Drydep     = F
             Is_Wetdep     = F

          CASE( 'CL2' )
             FullName      = 'Molecular chlorine'
             Formula       = 'Cl2'
             MW_g          = 71.0_fp
             Is_Gas        = T
             Is_Drydep     = F
             Is_Wetdep     = F
             Is_Photolysis = T

          CASE( 'CL2O2' )
             FullName      = 'Dichlorine dioxide'
             Formula       = 'Cl2O2'
             MW_g          = 103.0_fp
             Is_Gas        = T
             Is_Drydep     = F
             Is_Wetdep     = F
             Is_Photolysis = T

          CASE( 'CLNO2' )
             FullName      = 'Nitryl chloride'
             Formula       = 'ClNO2'
             MW_g          = 81.0_fp
             Is_Gas        = T
             Is_Drydep     = F
             Is_Wetdep     = F
             Is_Photolysis = T

          CASE( 'CLNO3' )
             FullName      = 'Chlorine nitrate'
             Formula       = 'ClNO3'
             MW_g          = 97.0_fp
             Is_Gas        = T
             Is_Drydep     = T
             Is_Wetdep     = F
             Is_Photolysis = T
             DD_F0         = 0.0_fp
             DD_Hstar_old  = 1.00e+20_fp

          CASE( 'CLO' )
             FullName      = 'Chlorine monoxide'
             Formula       = 'ClO'
             MW_g          = 51.0_fp
             Is_Gas        = T
             Is_Drydep     = F
             Is_Wetdep     = F
             Is_Photolysis = T

          CASE( 'CLOO' )
             FullName      = 'Chlorine dioxide'
             Formula       = 'ClOO'
             MW_g          = 67.0_fp
             Is_Gas        = T
             Is_Drydep     = F
             Is_Wetdep     = F
             Is_Photolysis = T

          CASE( 'CO',     'COUS',    'COEUR',  'COASIA', 'COOTH',  &
                'COBBAM', 'COBBAF',  'COBBAS', 'COBBOC', 'COBBEU', &
                'COBBNA', 'COBBOTH', 'COCH4',  'COBIOF', 'COISOP', &
                'COMONO', 'COMEOH',  'COACET', 'CONMVOC'            )

             SELECT CASE( TRIM( Name ) )
                CASE( 'CO'     )
                   FullName = 'Carbon monoxide'
                CASE( 'COUS'   )
                   FullName = 'Anthropogenic + biofuel CO emitted over the USA'
                CASE( 'COEUR'  )
                   FullName = 'Anthropogenic + biofuel CO emitted over Europe'
                CASE( 'COASIA' )
                   FullName = 'Anthropogenic + biofuel CO emitted over Asia'
                CASE( 'COOTH'  )
                   FullName = 'Anthropogenic + biofuel CO emitted everywhere else'
                CASE( 'COBBAM' )
                   FullName = 'Biomass burning CO emitted over South America'
                CASE( 'COBBAF' )
                   FullName = 'Biomass burning CO emitted over Africa'
                CASE( 'COBBAS' )
                   FullName = 'Biomass burning CO emitted over Asia'
                CASE( 'COBBOC' )
                   FullName = 'Biomass burning CO emitted over Oceania'
                CASE( 'COBBEU' )
                   FullName = 'Biomass burning CO emitted over Europe'
                CASE( 'COBBOTH' )
                   FullName = 'Biomass burning CO emitted everywhere else'
                CASE( 'COCH4'  )
                   FullName = 'CO produced from methane oxidation'
                CASE( 'COBIOF' )
                   FullName = 'CO produced from biofuels (whole world)'
                CASE( 'CONMVOC' )
                   FullName = 'CO produced from NMVOC oxidation'
                CASE( 'COISOP' )
                   FullName = 'CO produced from isoprene oxidation'
                CASE( 'COMONO' )
                   FullName = 'CO produced from monterpenes oxidation'
                CASE( 'COMEOH' )
                   FullName = 'CO produced from methanol oxidation'
                CASE( 'COACET' )
                   FullName = 'CO produced from acetone oxidation'
             END SELECT

             ! Set special default background for CO
             SELECT CASE( TRIM( Name ) )
                CASE( 'CO'     )
                   BackgroundVV = 1.0e-07_fp
                CASE DEFAULT
                   BackgroundVV = MISSING_VV
             END SELECT

             Formula       = 'CO'
             MW_g          = 28.0_fp
             Is_Gas        = T
             Is_Drydep     = F
             Is_Wetdep     = F
#if defined( NEW_HENRY_CONSTANTS )
             Henry_K0      = 9.70e-6_f8 * To_M_atm
             Henry_CR      = 1300.0_f8
#endif

          CASE( 'DHDC' )
             FullName      = 'Dihydroxyperoxide dicarbonyl'
             Formula       = 'C5H8O6'
             MW_g          = 164.0_fp
             Is_Gas        = T
             Is_Drydep     = F
             Is_Wetdep     = F

          CASE( 'DHDN' )
             ! DHDN uses the same DD_F0 and DD_Hstar_old values as ISOPN
             ! so that we can compute its drydep velocity explicitly.
             FullName      = 'C5 dihydroxydinitrate'
             Formula       = 'C5H10O8N2'
             MW_g          = 226.0_fp
             Is_Gas        = T
             Is_Drydep     = T
             Is_Wetdep     = T
             Is_Photolysis = T
             DD_F0         = 1.0_fp
             DD_Hstar_Old  = 2.00e+6_fp
             Henry_K0      = 2.00e+6_f8
             Henry_CR      = 9200.0_f8
             WD_RetFactor  = 2.0e-2_fp

          CASE( 'DHPCARP' )
             FullName = 'Dihydroxy a-formyl peroxy radical'
             Formula       = 'C5H9O7'
             MW_g          = 181.0_fp
             Is_Gas        = T
             Is_Drydep     = F
             Is_Wetdep     = F

          CASE( 'DMS' )
             FullName      = 'Dimethyl sulfide'
             Formula       = '(CH3)2S'
             MW_g          = 62.0_fp
             Is_Gas        = F
             Is_Drydep     = F
             Is_Wetdep     = F
             Henry_K0      = 4.80e-1_f8
             Henry_CR      = 3100.0_f8

          CASE( 'DST1', 'DSTAL1', 'NITD1', 'SO4D1' )

             ! These have identical properties except for the names
             SELECT CASE( Name )
                CASE( 'DST1' )
                   FullName = 'Dust aerosol, Reff = 0.7 microns'
                CASE( 'DSTAL1' )
                   FullName = 'Dust alkalinity, Reff = 0.7 microns'
                CASE( 'NITD1' )
                   FullName = 'Nitrate on dust, Reff = 0.7 microns'
                CASE( 'SO4D1' )
                   FullName = 'Sulfate on dust, Reff = 0.7 microns'
             END SELECT

             ! Do not reduce the Kc (cloud condensate -> precip) rate
             KcScale = (/ 1.0_fp, 1.0_fp, 1.0_fp /)

             ! Allow rainout of dust when T < 258K, because dust
             ! is considered to be IN.
             RainEff = (/ 1.0_fp, 1.0_fp, 0.0_fp /)

             Formula       = ''
             MW_g          = 29.0_fp
             Is_Gas        = F
             Is_Drydep     = T
             Is_Wetdep     = T
             Density       = 2500.0_fp
             Radius        = 7.3e-7_fp
             DD_DustDryDep = T
             DD_F0         = 0.0_fp
             DD_Hstar_Old  = 0.0_fp
             WD_AerScavEff = 1.0_fp
             WD_KcScaleFac = KcScale
             WD_RainoutEff = RainEff

          CASE( 'DST2', 'DSTAL2', 'NITD2', 'SO4D2' )

             ! These have identical properties except for the names
             SELECT CASE( Name )
                CASE( 'DST2' )
                   FullName = 'Dust aerosol, Reff = 1.4 microns'
                CASE( 'DSTAL2' )
                   FullName = 'Dust alkalinity, Reff = 1.4 microns'
                CASE( 'NITD2' )
                   FullName = 'Nitrate on dust, Reff = 1.4 microns'
                CASE( 'SO4D2' )
                   FullName = 'Sulfate on dust, Reff = 1.4 microns'
             END SELECT

             ! Do not reduce the Kc (cloud condensate -> precip) rate
             KcScale       = (/ 1.0_fp, 1.0_fp, 1.0_fp /)

             ! Allow rainout of dust when T < 258K, because dust
             ! is considered to be IN.
             RainEff       = (/ 1.0_fp, 1.0_fp, 0.0_fp /)

             Formula       = ''
             MW_g          = 29.0_fp
             Is_Gas        = F
             Is_Drydep     = T
             Is_Wetdep     = T
             Density       = 2650.0_fp
             Radius        = 1.4e-6_fp
             DD_DustDryDep = T
             DD_F0         = 0.0_fp
             DD_Hstar_Old  = 0.0_fp
             WD_AerScavEff = 1.0_fp
             WD_CoarseAer  = T
             WD_KcScaleFac = KcScale
             WD_RainoutEff = RainEff

          CASE( 'DST3', 'DSTAL3', 'NITD3', 'SO4D3' )

             ! These have identical properties except for the names
             SELECT CASE( Name )
                CASE( 'DST3' )
                   FullName = 'Dust aerosol, Reff = 2.4 microns'
                CASE( 'DSTAL3' )
                   FullName = 'Dust alkalinity, Reff = 2.4 microns'
                CASE( 'NITD3' )
                   FullName = 'Nitrate on dust, Reff = 2.4 microns'
                CASE( 'SO4D3' )
                   FullName = 'Sulfate on dust, Reff = 2.4 microns'
             END SELECT

             ! Do not reduce the Kc (cloud condensate -> precip) rate
             KcScale       = (/ 1.0_fp, 1.0_fp, 1.0_fp /)

             ! Allow rainout of dust when T < 258K, because dust
             ! is considered to be IN.
             RainEff       = (/ 1.0_fp, 1.0_fp, 0.0_fp /)

             Formula       = ''
             MW_g          = 29.0_fp
             Is_Gas        = F
             Is_Drydep     = T
             Is_Wetdep     = T
             Density       = 2650.0_fp
             Radius        = 2.4e-6_fp
             DD_DustDryDep = T
             DD_F0         = 0.0_fp
             DD_Hstar_Old  = 0.0_fp
             WD_AerScavEff = 1.0_fp
             WD_CoarseAer  = T
             WD_KcScaleFac = KcScale
             WD_RainoutEff = RainEff


          CASE( 'DST4', 'DSTAL4', 'NITD4', 'SO4D4' )

             ! These have identical properties except for the names
             SELECT CASE( Name )
                CASE( 'DST4' )
                   FullName = 'Dust aerosol, Reff = 4.5 microns'
                CASE( 'DSTAL4' )
                   FullName = 'Dust alkalinity, Reff = 4.5 microns'
                CASE( 'NITD4' )
                   FullName = 'Nitrate on dust, Reff = 4.5 microns'
                CASE( 'SO4D4' )
             END SELECT

             ! Do not reduce the Kc (cloud condensate -> precip) rate
             KcScale       = (/ 1.0_fp, 1.0_fp, 1.0_fp /)

             ! Allow rainout of dust when T < 258K, because dust
             ! is considered to be IN.
             RainEff       = (/ 1.0_fp, 1.0_fp, 0.0_fp /)

             Formula       = ''
             MW_g          = 29.0_fp
             Is_Gas        = F
             Is_Drydep     = T
             Is_Wetdep     = T
             Density       = 2650.0_fp
             Radius        = 4.50e-6_fp
             DD_DustDryDep = T
             DD_F0         = 0.0_fp
             DD_Hstar_Old  = 0.0_fp
             WD_AerScavEff = 1.0_fp
             WD_CoarseAer  = T
             WD_KcScaleFac = KcScale
             WD_RainoutEff = RainEff

          CASE( 'EOH' )
             FullName      = 'Ethanol'
             Formula       = 'C2H5OH'
             MW_g          = 46.07_fp
             EmMW_g        = 12.0_fp
             MolecRatio    = 2.0_fp
             Is_Gas        = T
             Is_Drydep     = T
             Is_Wetdep     = T
             DD_F0         = 0.0_fp
             DD_Hstar_old  = 1.9e+2_fp
             Henry_K0      = 1.9e+2_f8
             Henry_CR      = 6600.0_f8
             WD_RetFactor  = 2.0e-2_fp

          CASE( 'ETHLN' )
             ! ETHLN uses the same DD_F0 and DD_Hstar_old values as ISOPN
             ! so that we can compute its drydep velocity explicitly.
             FullName      = 'Ethanol nitrate'
             Formula       = 'CHOCH2ONO2'
             MW_g          = 105.0_fp
             Is_Gas        = T
             Is_Drydep     = T
             Is_Wetdep     = T
             Is_Photolysis = T
             DD_F0         = 1.0_fp
             DD_Hstar_old  = 2.00e+6_fp
             Henry_K0      = 2.00e+6_f8
             Henry_CR      = 9200.0_f8
             WD_RetFactor  = 2.0e-2_fp


          CASE( 'GLYC' )
             FullName      = 'Glycoaldehyde'
             Formula       = 'HOCH2CHO'
             MW_g          = 60.0_fp
             Is_Gas        = T
             Is_Drydep     = T
             Is_Wetdep     = T
             Is_Photolysis = T
             DD_F0         = 1.0_fp
#if defined( NEW_HENRY_CONSTANTS )
             Henry_K0      = 4.10e+2_f8 * To_M_atm
             Henry_CR      = 4600.0_f8
#else
             DD_Hstar_old  = 4.10e+4_fp
             Henry_K0      = 4.10e+4_f8
             Henry_CR      = 4600.0_f8
#endif
             WD_RetFactor  = 2.0e-2_fp

          CASE( 'GLYX' )
             FullName      = 'Glyoxal'
             Formula       = 'CHOCHO'
             MW_g          = 58.0_fp
             Is_Gas        = T
             Is_Drydep     = T
             Is_Wetdep     = T
             Is_Photolysis = T
             DD_F0         = 1.0_fp
             DD_Hstar_old  = 3.6e+5_fp
             Henry_K0      = 3.6e+5_f8
             Henry_CR      = 7200.0_f8
             WD_RetFactor  = 2.0e-2_fp

          CASE( 'H1211' )
             FullName      = 'Halon 1211, Freon 12B1'
             Formula       = 'CBrClF2'
             MW_g          = 165.0_fp
             Is_Gas        = T
             Is_Drydep     = F
             Is_Wetdep     = F
             Is_Photolysis = T

          CASE( 'H1301' )
             FullName      = 'Halon 1301, Freon 13B1'
             Formula       = 'CBrF3'
             MW_g          = 149.0_fp
             Is_Gas        = T
             Is_Drydep     = F
             Is_Wetdep     = F
             Is_Photolysis = T

          CASE( 'H2402' )
             FullName      = 'Halon 2402'
             Formula       = 'C2Br2F4'
             MW_g          = 260.0_fp
             Is_Gas        = T
             Is_Drydep     = F
             Is_Wetdep     = F
             Is_Photolysis = T

          CASE( 'H2O' )
             FullName      = 'Water vapor'
             Formula       = 'H2O'
             MW_g          = 18.0_fp
             BackgroundVV  = 1.839e-02_fp
             Is_Gas        = T
             Is_Drydep     = F
             Is_Wetdep     = F

          CASE( 'H2O2' )
             FullName      = 'Hydrogen peroxide'
             Formula       = 'H2O2'
             MW_g          = 34.0_fp
             BackgroundVV  = 4.0e-15_fp
             Is_Gas        = T
             Is_Drydep     = T
             Is_Wetdep     = T
             Is_Photolysis = T
             DD_F0         = 1.0_fp
#if defined( NEW_HENRY_CONSTANTS )
             Henry_K0      = 4.93e+5_f8 * To_M_atm
             Henry_CR      = 6600.0_f8
             Henry_pKa     = 11.6_f8
#else
             DD_Hstar_old  = 5.00e+7_fp
             Henry_K0      = 8.30e+4_f8
             Henry_CR      = 7400.0_f8
#endif
             WD_RetFactor  = 5e-2_fp
                              WD_LiqAndGas  = T
                              WD_ConvFacI2G = 4.36564e-1_fp

          CASE( 'HAC' )
             FullName      = 'Hydroxyacetone'
             Formula       = 'HOCH2C(O)CH3'
             MW_g          = 74.0_fp
             Is_Gas        = T
             Is_Drydep     = T
             Is_Wetdep     = T
             Is_Photolysis = T
             DD_F0         = 1.0_fp
#if defined( NEW_HENRY_CONSTANTS )
             Henry_K0      = 7.70e+1_f8 * To_M_atm
             Henry_CR      = 0.0_f8
#else
             DD_Hstar_old  = 1.40e+6_fp
             Henry_K0      = 1.40e+6_f8
             Henry_CR      = 7200.0_f8
#endif
             WD_RetFactor  = 2.0e-2_fp

          CASE( 'HBR' )
             FullName      = 'Hypobromic acid'
             Formula       = 'HBr'
             MW_g          = 81.0_fp
             Is_Gas        = T
             Is_Drydep     = T
             Is_Wetdep     = T
             DD_F0         = 0.0_fp
#if defined( NEW_HENRY_CONSTANTS )
             Henry_K0      = 2.40e-1_f8 * To_M_atm
             Henry_CR      = 370.0_f8
#else
             DD_Hstar_old  = 7.10e+15_fp
             Henry_K0      = 7.10e+13_f8
             Henry_CR      = 10200.0_f8
#endif
             WD_RetFactor  = 1.0_fp

          CASE( 'HC187' )
             FullName = 'Epoxide oxidation product m/z 187-189'

             ! HC187 uses the same DD_F0 and DD_Hstar_old values as HNO3,
             ! so that we can compute its drydep velocity explicitly.
             Formula       = ''
             MW_g          = 187.0_fp
             Is_Gas        = T
             Is_Drydep     = T
             Is_Wetdep     = F
             DD_F0         = 0.0_fp
#if defined( NEW_HENRY_CONSTANTS )
             Henry_K0      = 2.10e-2_f8 * To_M_atm
             Henry_CR      = 3400.0_f8
#else
             DD_Hstar_old  = 1.0e+14_fp
#endif

          CASE( 'HCFC123' )
             FullName      = 'HCFC-123, Freon 123'
             Formula       = 'C2HCl2F3'
             MW_g          = 117.0_fp
             Is_Gas        = T
             Is_Drydep     = F
             Is_Wetdep     = F
             Is_Photolysis = T

          CASE( 'HCFC141B' )
             FullName      = 'HCFC-141b, Freon 141b'
             Formula       = 'C(CH3)Cl2F'
             MW_g          = 117.0_fp
             Is_Gas        = T
             Is_Drydep     = F
             Is_Wetdep     = F
             Is_Photolysis = T

          CASE( 'HCFC142B' )
             FullName      = 'HCFC-142b, Freon 142b'
             Formula       = 'C(CH3)ClF2'
             MW_g          = 117.0_fp
             Is_Gas        = T
             Is_Drydep     = F
             Is_Wetdep     = F
             Is_Photolysis = T

          CASE( 'HCFC22' )
             FullName      = 'HCFC-22, Freon 22'
             Formula       = 'CHClF2'
             MW_g          = 86.0_fp
             Is_Gas        = T
             Is_Drydep     = F
             Is_Wetdep     = F
             Is_Photolysis = T

          CASE( 'HCL' )
             FullName      = 'Hydrochloric acid'
             Formula       = 'HCl'
             MW_g          = 36.0_fp
             Is_Gas        = T
             Is_Drydep     = T
             Is_Wetdep     = T
             DD_F0         = 0.0_fp
#if defined( MODEL_GEOS )
             DD_Hstar_old  = 2.05e+6_fp
             Henry_K0      = 7.10e+15_f8
#else
             DD_Hstar_old  = 2.05e+13_fp
             Henry_K0      = 7.00e+10_f8
#endif
             Henry_CR      = 11000.0_f8
             WD_RetFactor  = 1.0_fp

          CASE( 'HCOOH' )
             FullName      = 'Formic acid'
             Formula       = 'HCOOH'
             MW_g          = 46.0_fp
             Is_Gas        = T
             Is_Drydep     = T
             Is_Wetdep     = T
             DD_F0         = 1.0_fp
             DD_Hstar_old  = 8.90e+3_fp
             Henry_K0      = 8.90e+3_f8
             Henry_CR      = 6100.0_f8
                              WD_RetFactor  = 2.0e-2_fp

          CASE( 'HNO2' )
             FullName      = 'Nitrous acid'
             Formula       = 'HNO2'
             MW_g          = 47.0_fp
             BackgroundVV  = 4.0e-15_fp
             Is_Gas        = T
             Is_Drydep     = F
             Is_Wetdep     = F
             Is_Photolysis = T

          CASE( 'HNO3' )

             !%%% NOTE: HNO3 dry-deposits like a gas, but wet-deposits
             !%%% like an aerosol.  Therefore we need to define HNO3
             !%%% with both gas-phase and aerosol parameters. (bmy, 9/28/15)

             ! Do not reduce the Kc (cloud condensate -> precip) rate
             KcScale       = (/ 1.0_fp, 1.0_fp, 1.0_fp /)

             ! Allow rainout of HNO3 when T < 258K, becasue HNO3
             ! is considered to be IN.
             RainEff       = (/ 1.0_fp, 1.0_fp, 1.0_fp /)

             FullName      = 'Nitric acid'
             Formula       = 'HNO3'
             MW_g          = 63.0_fp
             BackgroundVV  = 4.0e-15_fp
             Is_Gas        = T
             Is_Drydep     = T
             Is_Wetdep     = T
             Is_Photolysis = T
             DD_F0         = 0.0_fp
#if defined( NEW_HENRY_CONSTANTS )
             Henry_K0      = 2.10e+3_f8 * To_M_atm
             Henry_CR      = 8700.0_f8
#else
             DD_Hstar_old  = 1.0e+14_fp
             Henry_K0      = 8.3e+4_f8
             Henry_CR      = 7400.0_f8
#endif
                              WD_AerScavEff = 1.0_fp
                              WD_KcScaleFac = KcScale
                              WD_RainoutEff = RainEff

          CASE( 'HNO4' )
             FullName      = 'Peroxynitric acid'
             Formula       = 'HNO4'
             MW_g          = 79.0_fp
             BackgroundVV  = 4.0e-15_fp
             Is_Gas        = T
             Is_Drydep     = F
             Is_Wetdep     = F
             Is_Photolysis = T
#if defined( NEW_HENRY_CONSTANTS )
             !!!Henry_K0      = 3.90e1X_f8 * To_M_atm
             Henry_K0      = 3.90e-1_f8 * To_M_atm
             Henry_CR      = 8400.0_f8
#endif

          CASE( 'HOBR' )
             FullName      = 'Hypobromous acid'
             Formula       = 'HOBr'
             MW_g          = 97.0_fp
             Is_Gas        = T
             Is_Drydep     = T
             Is_Wetdep     = T
             Is_Photolysis = T
             DD_F0         = 0.0_fp
#if defined( NEW_HENRY_CONSTANTS )
             Henry_K0      = 1.30e+0_f8 * To_M_atm
             Henry_CR      = 4000.0_f8
#else
             DD_Hstar_old  = 6.10e+3_fp
             Henry_K0      = 6.10e+3_f8
             Henry_CR      = 6014.0_f8
#endif
             WD_RetFactor  = 0.0_fp

          CASE( 'HOCL' )
             FullName      = 'Hypochlorous acid'
             Formula       = 'HOCl'
             MW_g          = 52.0_fp
             Is_Gas        = T
             Is_Drydep     = T
             Is_Wetdep     = T
             Is_Photolysis = T
             DD_F0         = 0.0_fp
             DD_Hstar_old  = 6.50e+2_fp
             Henry_K0      = 6.50e+2_f8
             Henry_CR      = 5900.0_f8
             WD_RetFactor  = 0.0_fp

          CASE( 'HONIT' )
             FullName = '2nd gen monoterpene organic nitrate'

             ! HONIT uses the same DD_F0 and DD_Hstar_old values as ISOPN
             ! so that we can compute its drydep velocity explicitly.
             Formula       = ''
             MW_g          = 215.0_fp
             Is_Gas        = T
             Is_Drydep     = T
             Is_Wetdep     = T
             Is_Photolysis = T
             DD_F0         = 1.0_fp
             DD_Hstar_old  = 2.00e+6_fp
             Henry_K0      = 2.69e+13_f8
             Henry_CR      = 5487.0_f8
             WD_RetFactor  = 2.0e-2_fp

          CASE( 'HPALD' )
             FullName      = 'Hydroperoxyaldehydes'
             Formula       = 'C5H8O3'
             MW_g          = 116.0_fp
             Is_Gas        = T
             Is_Drydep     = T
             Is_Wetdep     = F
             Is_Photolysis = T
             DD_F0         = 0.0_fp
             DD_Hstar_old  = 4.0e+4_fp

          CASE( 'HPC52O2' )
             FullName      = 'HOC52O2'
             Formula       = 'OOCC(C(C=O)O)(O[O])C'
             MW_g          = 165.0_fp
             Is_Gas        = T
             Is_Drydep     = F
             Is_Wetdep     = F

          CASE( 'IEPOXA' )
             FullName      = 'trans-Beta isoprene epoxydiol'
             Formula       = 'C5H10O3'
             MW_g          = 118.0_fp
             Is_Gas        = T
             Is_Drydep     = T
             Is_Wetdep     = T
             DD_F0         = 1.0_fp
             DD_Hstar_old  = 8.00e+7_fp
             Henry_K0      = 8.00e+7_f8
             Henry_CR      = 0.0_f8
             WD_RetFactor  = 2.0e-2_fp

          CASE( 'IEPOXB' )
             FullName      = 'cis-Beta isoprene epoxydiol'
             Formula       = 'C5H10O3'
             MW_g          = 118.0_fp
             Is_Gas        = T
             Is_Drydep     = T
             Is_Wetdep     = T
             DD_F0         = 1.0_fp
             DD_Hstar_old  = 8.00e+7_fp
             Henry_K0      = 8.00e+7_f8
             Henry_CR      = 0.0_f8
             WD_RetFactor  = 2.0e-2_fp

          CASE( 'IEPOXD' )
             FullName      = 'Delta isoprene epoxydiol'
             Formula       = 'C5H10O3'
             MW_g          = 118.0_fp
             Is_Gas        = T
             Is_Drydep     = T
             Is_Wetdep     = T
             DD_F0         = 1.0_fp
             DD_Hstar_old  = 8.00e+7_fp
             Henry_K0      = 8.00e+7_f8
             Henry_CR      = 0.0_f8
             WD_RetFactor  = 2.0e-2_fp

          CASE( 'IMAE' )
             FullName = 'C4 epoxide from oxidation of MPAN (PMN)'
             Formula       = 'C4H6O3'
             MW_g          = 102.0_fp
             Is_Gas        = T
             Is_Drydep     = T
             Is_Wetdep     = T
             DD_F0         = 1.0_fp
             ! DD_Hstar based on Pye et al. (2013)
             DD_Hstar_old  = 1.20e+5_fp
             Henry_K0      = 1.20e+5_f8
             Henry_CR      = 7200.0_f8
             WD_RetFactor  = 2.0e-2_fp

          CASE( 'INDIOL' )
             FullName      = 'Generic aerosol-phase organonitrate hydrolysis product'

             ! Halve the Kc (cloud condensate -> precip) rate
             ! for the temperature range 237 K <= T < 258 K.
             KcScale       = (/ 1.0_fp, 0.5_fp, 1.0_fp /)

             ! Turn off rainout only when 237 K <= T < 258K.
             ! NOTE: Rainout efficiency is 0.8 because these are SOA species.
             RainEff       = (/ 0.8_fp, 0.0_fp, 0.8_fp /)

             Formula       = ''
             MW_g          = 102.0_fp
             Is_Gas        = F
             Is_Drydep     = T
             Is_Wetdep     = T
             DD_DvzAerSnow = 0.03_fp
             DD_F0         = 0.0_fp
             DD_HStar_old  = 0.0_fp
             WD_AerScavEff = 0.8_fp
             WD_KcScaleFac = KcScale
             WD_RainoutEff = RainEff

          CASE( 'IONITA' )
             FullName      = 'Aer-phase organic nitrate from isoprene precursors'

             ! Halve the Kc (cloud condensate -> precip) rate
             ! for the temperature range 237 K <= T < 258 K.
             KcScale       = (/ 1.0_fp, 0.5_fp, 1.0_fp /)

             ! Turn off rainout only when 237 K <= T < 258K.
             ! NOTE: Rainout efficiency is 0.8 because these are SOA species.
             RainEff       = (/ 0.8_fp, 0.0_fp, 0.8_fp /)

             Formula       = ''
             MW_g          = 14.0_fp
             Is_Gas        = F
             Is_Drydep     = T
             Is_Wetdep     = T
             DD_DvzAerSnow = 0.03_fp
             DD_F0         = 0.0_fp
             DD_Hstar_old  = 0.0_fp
             WD_AerScavEff = 0.8_fp
             WD_KcScaleFac = KcScale
             WD_RainoutEff = RainEff

          CASE( 'IPMN' )
             FullName = 'Peroxymethacroyl nitrate (PMN) from isoprene oxidation'

             ! PMN uses the same DD_F0 and DD_Hstar_old values as PAN
             ! so that we can compute its drydep velocity explicitly.
             ! (bmy, 5/19/16)
             Formula       = 'CH2=C(CH3)C(O)OONO2'
             MW_g          = 147.0_fp
             Is_Gas        = T
             Is_Drydep     = T
             Is_Wetdep     = F
             Is_Photolysis = T
             DD_F0         = 1.0_fp
#if defined( NEW_HENRY_CONSTANTS )
             Henry_K0      = 1.70e-2_f8 * To_M_atm
#else
             DD_Hstar_old  = 3.60_fp
#endif

          CASE( 'ISN1' )
             ! ISN1 uses the same DD_F0 and DD_Hstar_old values as ISOPN
             ! so that we can compute its drydep velocity explicitly.
             FullName      = 'Nighttime isoprene nitrate'
             Formula       = 'C5H8NO4'
             MW_g          = 147.0_fp
             Is_Gas        = T
             Is_Drydep     = T
             Is_Wetdep     = T
             Is_Photolysis = T
             DD_F0         = 1.0_fp
             DD_Hstar_old  = 2.00e+6_fp
             Henry_K0      = 2.00e+6_f8
             Henry_CR      = 9200.0_f8
             WD_RetFactor  = 2.0e-2_fp

          CASE( 'ISN1OA' )
             FullName      = 'Aer-phase 2nd-gen hydroxynitrates from ISOP+NO3'

             ! Halve the Kc (cloud condensate -> precip) rate
             ! for the temperature range 237 K <= T < 258 K.
             KcScale       = (/ 1.0_fp, 0.5_fp, 1.0_fp /)

             ! Turn off rainout only when 237 K <= T < 258K.
             ! NOTE: Rainout efficiency is 0.8 because these are SOA species.
             RainEff       = (/ 0.8_fp, 0.0_fp, 0.8_fp /)

             Formula       = ''
             MW_g          = 226.0_fp
             Is_Gas        = F
             Is_Drydep     = T
             Is_Wetdep     = T
             DD_DvzAerSnow = 0.03_fp
             DD_F0         = 0.0_fp
             DD_HStar_old  = 0.0_fp
             WD_AerScavEff = 0.8_fp
             WD_KcScaleFac = KcScale
             WD_RainoutEff = RainEff

          CASE( 'ISN1OG' )
             FullName      = 'Gas-phase 2nd-gen hydroxynitrates from ISOP+NO3'
             Formula       = ''
             MW_g          = 226.0_fp
             Is_Gas        = T
             Is_Drydep     = T
             Is_Wetdep     = T
             DD_F0         = 1.0_fp
             DD_Hstar_old  = 2.30e+4_fp
             Henry_K0      = 2.30e+4_f8
             Henry_CR      = 9200.0_f8
             WD_RetFactor  = 2.0e-2_fp

          CASE( 'ISOA1', 'ISOA2', 'ISOA3' )
             FullName      = 'Lumped semivolatile aer products of isoprene oxidation'

             ! Halve the Kc (cloud condensate -> precip) rate
             ! for the temperature range 237 K <= T < 258 K.
             KcScale       = (/ 1.0_fp, 0.5_fp, 1.0_fp /)

             ! Turn off rainout only when 237 K <= T < 258K.
             ! NOTE: Rainout efficiency is 0.8 because these are SOA species.
             RainEff       = (/ 0.8_fp, 0.0_fp, 0.8_fp /)

             Formula       = ''
             MW_g          = 150.0_fp
             Is_Gas        = F
             Is_Drydep     = T
             Is_Wetdep     = T
             DD_DvzAerSnow = 0.03_fp
             DD_F0         = 0.0_fp
             DD_HStar_old  = 0.0_fp
             WD_AerScavEff = 0.8_fp
             WD_KcScaleFac = KcScale
             WD_RainoutEff = RainEff

          CASE( 'ISOG1', 'ISOG2', 'ISOG3' )
             FullName      = 'Lumped semivolatile gas products of isoprene oxidation'
             Formula       = ''
             MW_g          = 150.0_fp
             Is_Gas        = T
             Is_Drydep     = T
             Is_Wetdep     = T
             DD_F0         = 0.0_fp
             DD_Hstar_old  = 1.00e+5_fp
             Henry_K0      = 1.00e+5_f8
             Henry_CR      = 6039.0_f8
             WD_RetFactor  = 2.0e-2_fp

          CASE( 'ISOP' )
             FullName      = 'Isoprene'
             Formula       = 'CH2=C(CH3)CH=CH2'
             MW_g          = 68.12_fp
             EmMW_g        = 12.0_fp
             MolecRatio    = 5e+0_fp
             Is_Gas        = T
             Is_Drydep     = F
             Is_Wetdep     = F
#if defined( NEW_HENRY_CONSTANTS )
             Henry_K0      = 3.40e-4_f8 * To_M_atm
             Henry_CR      = 4400.0_f8
#endif

          CASE( 'ISOPNB' )
             FullName      = 'Isoprene nitrate Beta'
             Formula       = 'C5H9NO4'
             MW_g          = 147.0_fp
             Is_Gas        = T
             Is_Drydep     = T
             Is_Wetdep     = T
             Is_Photolysis = T
             DD_F0         = 1.0_fp
#if defined( NEW_HENRY_CONSTANTS )
             Henry_K0      = 1.97e+4_f8 * To_M_atm
#else
             DD_Hstar_old  = 2.00e+6_fp
             Henry_K0      = 2.00e+6_f8
             Henry_CR      = 9200.0_f8
#endif
             WD_RetFactor  = 2.0e-2_fp

          CASE( 'ISOPND'  )
             FullName      = 'Isoprene nitrate Delta'
             Formula       = 'C5H9NO4'
             MW_g          = 147.0_fp
             Is_Gas        = T
             Is_Drydep     = T
             Is_Wetdep     = T
             Is_Photolysis = T
             DD_F0         = 1.0_fp
#if defined( NEW_HENRY_CONSTANTS )
             Henry_K0      = 1.97e+4_f8 * To_M_atm
#else
             DD_Hstar_old  = 2.00e+6_fp
             Henry_K0      = 2.00e+6_f8
             Henry_CR      = 9200.0_f8
#endif
             WD_RetFactor  = 2.0e-2_fp

          CASE( 'LIMO' )
             FullName      = 'Limonene'
             Formula       = 'C10H16'
             MW_g          = 136.23_fp
             Is_Gas        = T
             Is_Drydep     = T
             Is_Wetdep     = T
             DD_F0         = 0.0_fp
             DD_Hstar_old  = 7.00e-2_fp
             Henry_K0      = 7.00e-2_f8
             Henry_CR      = 0.0_f8
             WD_RetFactor  = 2.0e-2_fp

          CASE( 'LVOC' )
             FullName = 'Gas-phase low-volatility non-IEPOX product of RIP ox'
             Formula       = 'C5H14O5'
             MW_g          = 154.0_fp
             Is_Gas        = T
             Is_Drydep     = T
             Is_Wetdep     = T
             DD_F0         = 1.0_fp
             DD_Hstar_old  = 1.00e+8_fp
             Henry_K0      = 1.00e+8_f8
             Henry_CR      = 7200.0_f8
             WD_RetFactor  = 2.0e-2_fp

          CASE( 'LVOCOA' )
             FullName      = 'Aer-phase low-volatility non-IEPOX product of RIP ox'

             ! Halve the Kc (cloud condensate -> precip) rate
             ! for the temperature range 237 K <= T < 258 K.
             KcScale       = (/ 1.0_fp, 0.5_fp, 1.0_fp /)

             ! Turn off rainout only when 237 K <= T < 258K.
             ! NOTE: Rainout efficiency is 0.8 because these are SOA species.
             RainEff       = (/ 0.8_fp, 0.0_fp, 0.8_fp /)

             Formula       = 'C5H14O5'
             MW_g          = 154.0_fp
             Is_Gas        = F
             Is_Drydep     = T
             Is_Wetdep     = T
             DD_DvzAerSnow = 0.03_fp
             DD_F0         = 0.0_fp
             DD_HStar_old  = 0.0_fp
             WD_AerScavEff = 0.8_fp
             WD_KcScaleFac = KcScale
             WD_RainoutEff = RainEff

          CASE( 'MACR' )
             FullName      = 'Methacrolein'
             Formula       = 'CH2=C(CH3)CHO'
             MW_g          = 70.0_fp
             Is_Gas        = T
             Is_Drydep     = T
             Is_Wetdep     = F
             Is_Photolysis = T
             DD_F0         = 1.0_fp
#if defined( NEW_HENRY_CONSTANTS )
             Henry_K0      = 4.8e-2_f8 * To_M_atm
             Henry_CR      = 4300.0_f8
#else
             DD_Hstar_old  = 6.5e+0_fp
#endif

          CASE( 'MACRN' )
             FullName      = 'Nitrate from MACR'
             Formula       = 'HOCH2C(ONO2)(CH3)CHO'
             MW_g          = 149.0_fp
             Is_Gas        = T
             Is_Drydep     = T
             Is_Wetdep     = T
             Is_Photolysis = T
             DD_F0         = 1.0_fp
#if defined( NEW_HENRY_CONSTANTS )
             Henry_K0      = 1.97e+4_f8 * To_M_atm
#else
             DD_Hstar_old  = 2.00e+6_fp
             Henry_K0      = 2.00e+6_f8
             Henry_CR      = 9200.0_f8
#endif
             WD_RetFactor  = 2.0e-2_fp

          CASE( 'MAP' )
              FullName      = 'Peroxyacetic acid'
              Formula       = 'CH3C(O)OOH'
              MW_g          = 76.0_fp
              Is_Gas        = T
              Is_Drydep     = T
              Is_Wetdep     = T
              Is_Photolysis = T
              DD_F0         = 1.0_fp
#if defined( NEW_HENRY_CONSTANTS )
              Henry_K0      = 8.30e+0_f8 * To_M_atm
              Henry_CR      = 5300.0_f8
#else
              DD_Hstar_old  = 8.40e+2_fp
              Henry_K0      = 8.40e+2_f8
              Henry_CR      = 5300.0_f8
#endif
              WD_RetFactor  = 2.0e-2_fp

          CASE( 'MEK' )
             FullName      = 'Methyl Ethyl Ketone'
             Formula       = 'RC(O)R'
             MW_g          = 72.11_fp
             EmMW_g        = 12.0_fp
             MolecRatio    = 4.0_fp
             Is_Gas        = T
             Is_Drydep     = F
             Is_Wetdep     = F
             Is_Photolysis = T
#if defined( NEW_HENRY_CONSTANTS )
             Henry_K0      = 2.90e+02_f8 * To_M_atm
             Henry_CR      = 5700.0_f8
#endif

          CASE( 'MGLY' )
             FullName      = 'Methylglyoxal'
             Formula       = 'CH3COCHO'
             MW_g          = 72.0_fp
             Is_Gas        = T
             Is_Drydep     = T
             Is_Wetdep     = T
             Is_Photolysis = T
             DD_F0         = 1.0_fp
             DD_Hstar_old  = 3.7e+3_fp
             Henry_K0      = 3.7e+3_f8
             Henry_CR      = 7500.0_f8
             WD_RetFactor  = 2.0e-2_fp

          CASE( 'MOBA' )
             FullName      = '5C acid from isoprene'
             Formula       = 'HOC(=O)C(CH3)=CHCHO'
             MW_g          = 114.0_fp
             Is_Gas        = T
             Is_Drydep     = F
             Is_Wetdep     = T
#if defined( NEW_HENRY_CONSTANTS )
             Henry_K0      = 2.27e+2_f8 * To_M_atm
             Henry_CR      = 6300.0_f8
#else
             Henry_K0      = 2.30e+4_f8
             Henry_CR      = 6300.0_f8
#endif
             WD_RetFactor  = 2.0e-2_fp

          CASE( 'MONITA' )
             FullName      = 'Aer-phase organic nitrate from monoterpene precursors'

             ! Halve the Kc (cloud condensate -> precip) rate
             ! for the temperature range 237 K <= T < 258 K.
             KcScale       = (/ 1.0_fp, 0.5_fp, 1.0_fp /)

             ! Turn off rainout only when 237 K <= T < 258K.
             ! NOTE: Rainout efficiency is 0.8 because these are SOA species.
             RainEff       = (/ 0.8_fp, 0.0_fp, 0.8_fp /)

             Formula       = ''
             MW_g          = 14.0_fp
             Is_Gas        = F
             Is_Drydep     = T
             Is_Wetdep     = T
             DD_DvzAerSnow = 0.03_fp
             DD_F0         = 0.0_fp
             DD_Hstar_old  = 0.0_fp
             WD_AerScavEff = 0.8_fp
             WD_KcScaleFac = KcScale
             WD_RainoutEff = RainEff

          CASE( 'MONITS' )
             FullName = 'Saturated 1st gen monoterpene organic nitrate'

             ! MONITS uses the same DD_F0 and DD_Hstar_old values as ISOPN
             ! so that we can compute its drydep velocity explicitly.
             Formula       = ''
             MW_g          = 215.0_fp
             Is_Gas        = T
             Is_Drydep     = T
             Is_Wetdep     = T
             Is_Photolysis = T
             DD_F0         = 1.0_fp
             DD_Hstar_old  = 2.00e+6_fp
             Henry_K0      = 1.70e+4_f8
             Henry_CR      = 9200.0_f8
             WD_RetFactor  = 2.0e-2_fp

          CASE( 'MONITU' )
             FullName = 'Unsaturated 1st gen monoterpene organic nitrate'

             ! MONITU uses the same DD_F0 and DD_Hstar_old values as ISOPN
             ! so that we can compute its drydep velocity explicitly.
             FullName      = FullName
             Formula       = ''
             MW_g          = 215.0_fp
             Is_Gas        = T
             Is_Drydep     = T
             Is_Wetdep     = T
             Is_Photolysis = T
             DD_F0         = 1.0_fp
             DD_Hstar_old  = 2.00e+6_fp
             Henry_K0      = 1.70e+4_f8
             Henry_CR      = 9200.0_f8
             WD_RetFactor  = 2.0e-2_fp

          CASE( 'MOPI', 'MOPO' )

             ! MOPO is treated the same as OCPO and
             ! MOPI is treated the same as OCPI (msj, krt, 8/23/17).

             ! These have mostly identical properties
             ! Turn off rainout for hydrophobic OC, for all temperatures.
             SELECT CASE( Name )

                CASE( 'MOPI' )
                   FullName = 'Hydrophilic marine organic carbon aerosol'

                   ! Halve the Kc (cloud condensate -> precip) rate
                   ! for the temperature range 237 K <= T < 258 K.
                   KcScale  = (/ 1.0_fp, 0.5_fp, 1.0_fp /)

                   ! Turn off rainout only when 237 K <= T < 258K.
                   RainEff  = (/ 1.0_fp, 0.0_fp, 1.0_fp /)

                CASE( 'MOPO' )
                   Fullname = 'Hydrophobic marine organic carbon aerosol'

                   ! For all temperatures:
                   ! (1) Halve the Kc (cloud condensate -> precip) rate
                   ! (2) Turn off rainout (OCPO is hydrophobic)
                   KcScale  = (/ 0.5_fp, 0.5_fp, 0.5_fp /)
                   RainEff  = (/ 0.0_fp, 0.0_fp, 0.0_fp /)

             END SELECT

             Formula       = ''
             MW_g          = 12.01_fp
             EmMW_g        = 12.0_fp
             Is_Gas        = F
             Is_Drydep     = T
             Is_Wetdep     = T
             Density       = 1300.0_fp
             DD_DvzAerSnow = 0.03_fp
             DD_F0         = 0.0_fp
             DD_Hstar_Old  = 0.0_fp
             WD_AerScavEff = 1.0_fp
             WD_KcScaleFac = KcScale
             WD_RainoutEff = RainEff

          CASE( 'MP', 'CH3OOH' )
             FullName      = 'Methyl hydro peroxide'
             Formula       = 'CH3OOH'
             MW_g          = 48.0_fp
             BackgroundVV  = 4.0e-15_fp
             Is_Gas        = T
             Is_Drydep     = F
             Is_Wetdep     = T
             Is_Photolysis = T
#if defined( NEW_HENRY_CONSTANTS )
             Henry_K0      = 2.90e+0_f8 * To_M_atm
             Henry_CR      = 5200.0_f8
#else
             Henry_K0      = 3.10e+2_f8
             Henry_CR      = 5200.0_f8
#endif
             WD_RetFactor  = 2.0e-2_fp

          CASE( 'MPN' )
             FullName      = 'Methyl peroxy nitrate'
             Formula       = 'CH3O2NO2'
             MW_g          = 93.0_fp
             Is_Gas        = T
             Is_Drydep     = F
             Is_Wetdep     = F
             Is_Photolysis = T

          CASE( 'MSA' )

             ! Halve the Kc (cloud condensate -> precip) rate
             ! for the temperature range 237 K <= T < 258 K.
             KcScale       = (/ 1.0_fp, 0.5_fp, 1.0_fp /)

             ! Turn off rainout only when 237 K <= T < 258K.
             RainEff       = (/ 1.0_fp, 0.0_fp, 1.0_fp /)

             ! Enforce minimum dry deposition velocity (Vd) for MSA
             ! (cf. Mian Chin's GOCART model)
             ! Minimum Vd over snow/ice : 0.01 cm/s
             ! Minimum Vd over land     : 0.01 cm/s
             DvzMinVal     = (/ 0.01_fp, 0.01_fp /)

             FullName      = 'Methyl sulfonic acid'
             Formula       = 'CH4SO3'
             MW_g          = 96.0_fp
             Is_Gas        = F
             Is_Drydep     = T
             Is_Wetdep     = T
             DD_DvzAerSnow = 0.03_fp
             DD_DvzMinVal  = DvzMinVal
             DD_F0         = 0.0_fp
             DD_Hstar_old  = 0.0_fp
             WD_AerScavEff = 1.0_fp
             WD_KcScaleFac = KcScale
             WD_RainoutEff = RainEff

          CASE( 'MTPA' )
             FullName      = 'a-pinene, b-pinene, sabinene, carene'
             Formula       = ''
             MW_g          = 136.23_fp
             Is_Gas        = T
             Is_Drydep     = T
             Is_Wetdep     = T
             DD_F0         = 0.0_fp
             DD_Hstar_old  = 4.90e-2_fp
             Henry_K0      = 4.90e-2_f8
             Henry_CR      = 0.0_f8
             WD_RetFactor  = 2.0e-2_fp

          CASE( 'MTPO' )
             FullName = 'Terpinene, terpinolene, myrcene, ocimene, other monoterpenes'
             Formula       = ''
             MW_g          = 136.23_fp
             MolecRatio    = 1.0_fp
             Is_Gas        = T
             Is_Drydep     = T
             Is_Wetdep     = T
             DD_F0         = 0.0_fp
             DD_Hstar_old  = 4.90e-2_fp
             Henry_K0      = 4.90e-2_f8
             Henry_CR      = 0.0_f8
             WD_RetFactor  = 2.0e-2_fp

          CASE( 'MVK' )
             FullName      = 'Methyl vinyl ketone'
             Formula       = 'CH2=CHC(=O)CH3'
             MW_g          = 70.0_fp
             Is_Gas        = T
             Is_Drydep     = T
             Is_Wetdep     = F
             Is_Photolysis = T
             DD_F0         = 1.0_fp
#if defined( NEW_HENRY_CONSTANTS )
             Henry_K0      = 2.6e-1_f8 * To_M_atm
             Henry_CR      = 4800.0_f8
#else
             DD_Hstar_old  = 4.4e+1_fp
#endif

          CASE( 'MVKN' )
             FullName      = 'Nitrate from MVK'
             Formula       = 'HOCH2CH(ONO2)C(=O)CH3'
             MW_g          = 149.0_fp
             Is_Gas        = T
             Is_Drydep     = T
             Is_Wetdep     = T
             Is_Photolysis = T
             DD_F0         = 1.0_fp
#if defined( NEW_HENRY_CONSTANTS )
             Henry_K0      = 1.97e+4_f8 * To_M_atm
#else
             DD_Hstar_old  = 2.00e+6_fp
             Henry_K0      = 2.00e+6_f8
             Henry_CR      = 9200.0_f8
#endif
             WD_RetFactor  = 2.0e-2_fp
                              
          CASE( 'N2O' )
             FullName      = 'Nitrous oxide'
             Formula       = 'N2O'
             MW_g          = 44.0_fp
             BackgroundVV  = 3.0e-07_fp
             Is_Gas        = T
             Is_Drydep     = F
             Is_Wetdep     = F
             Is_Photolysis = T

          CASE( 'N2O5' )

             ! N2O5 uses the same DD_F0 and DD_Hstar_old values as HNO3,
             ! so that we can compute its drydep velocity explicitly.
             ! (bmy, 5/19/16)
             FullName      = 'Dinitrogen pentoxide'
             Formula       = 'N2O5'
             MW_g          = 108.0_fp
             BackgroundVV  = 4.0e-15_fp
             Is_Gas        = T
             Is_Drydep     = T
             Is_Wetdep     = F
             Is_Photolysis = T
             DD_F0         = 0.0_fp
#if defined( NEW_HENRY_CONSTANTS )
             Henry_K0      = 2.10e-2_f8 * To_M_atm
             Henry_CR      = 3400.0_f8
#else
             DD_Hstar_old  = 1.0e+14_fp
#endif


          CASE( 'NAP' )
              FullName      = 'Naphtalene/IVOC surrogate'
              Formula       = 'C10H8'
              MW_g          = 128.17_fp
              EmMw_g        = 12.0_fp
              MolecRatio    = 10.0_fp
              Is_Gas        = T
              Is_Drydep     = F
              Is_Wetdep     = F

          CASE( 'NH3' )

             ! Enforce minimum dry deposition velocity (Vd) for NH3
             ! (cf. Mian Chin's GOCART model)
             ! Minimum Vd over snow/ice : 0.2 cm/s
             ! Minimum Vd over land     : 0.3 cm/s
             DvzMinVal = (/ 0.2_fp, 0.3_fp /)

             FullName      = 'Ammonia'
             Formula       = 'NH3'
             MW_g          = 17.0_fp
             Is_Gas        = T
             Is_Drydep     = T
             Is_Wetdep     = T
             DD_DvzAerSnow = 0.03_fp
             DD_DvzMinVal  = DvzMinVal
             DD_F0         = 0.0_fp
#if defined( NEW_HENRY_CONSTANTS )
             Henry_K0      = 5.90e-1_f8 * To_M_atm
             Henry_CR      = 4200.0_f8
#else
             DD_Hstar_old  = 2.0e+4_fp
             Henry_K0      = 3.30e+6_f8
             Henry_CR      = 4100.0_f8
#endif
             WD_RetFactor  = 5.0e-2_fp
             WD_LiqAndGas  = T
             WD_ConvFacI2G = 6.17395e-1_fp

          CASE( 'NH4' )

             ! Halve the Kc (cloud condensate -> precip) rate
             ! for the temperature range 237 K <= T < 258 K.
             KcScale       = (/ 1.0_fp, 0.5_fp, 1.0_fp /)

             ! Turn off rainout only when 237 K <= T < 258K.
             RainEff       = (/ 1.0_fp, 0.0_fp, 1.0_fp /)

             ! Enforce minimum dry deposition velocity (Vd) for NH4
             ! (cf. Mian Chin's GOCART model)
             ! Minimum Vd over snow/ice : 0.01 cm/s
             ! Minimum Vd over land     : 0.01 cm/s
             DvzMinVal     = (/ 0.01_fp, 0.01_fp /)

             FullName      = 'Ammonium'
             Formula       = 'NH4'
             MW_g          = 18.0_fp
             Is_Gas        = F
             Is_Drydep     = T
             Is_Wetdep     = T
             DD_DvzAerSnow = 0.03_fp
             DD_DvzMinVal  = DvzMinVal
             DD_F0         = 0.0_fp
             DD_Hstar_old  = 0.0_fp
             WD_AerScavEff = 1.0_fp
             WD_KcScaleFac = KcScale
             WD_RainoutEff = RainEff

          CASE( 'NIT' )

             ! Halve the Kc (cloud condensate -> precip) rate
             ! for the temperature range 237 K <= T < 258 K.
             KcScale   = (/ 1.0_fp, 0.5_fp, 1.0_fp /)

             ! Turn off rainout only when 237 K <= T < 258K.
             RainEff   = (/ 1.0_fp, 0.0_fp, 1.0_fp /)

             ! Enforce minimum dry deposition velocity (Vd) for NIT
             ! (cf. Mian Chin's GOCART model)
             ! Minimum Vd over snow/ice : 0.01 cm/s
             ! Minimum Vd over land     : 0.01 cm/s
             DvzMinVal = (/ 0.01_fp, 0.01_fp /)

             FullName      = 'Inorganic nitrates'
             Formula       = ''
             MW_g          = 62.0_fp
             Is_Gas        = F
             Is_Drydep     = T
             Is_Wetdep     = T
             DD_DvzAerSnow = 0.03_fp
             DD_DvzMinVal  = DvzMinVal
             DD_F0         = 0.0_fp
             DD_Hstar_Old  = 0.0_fp
             WD_AerScavEff = 1.0_fp
             WD_KcScaleFac = KcScale
             WD_RainoutEff = RainEff

          CASE( 'NITS' )

             ! NOTE: NITs (NIT on coarse sea salt aerosol) needs to use
             ! the same molecular weight as coarse sea salt (= 31.4 g/mol)
             ! instead of NIT (= 62 g/mol).
             !
             ! Becky Alexander (beckya@atmos.washington.edu) wrote:
             !
             !   "The reason for using sea salt's MW for NITss is that ...
             !    [it is] ... essentially internally mixed with coarse sea
             !    salt aerosol (SALC).  As coarse sea salt aerosol likely
             !    dominates the mass of ... [NITs] ...  it is appropriate
             !    to use sea salt's MW.  Another explanation is that since
             !    NITs ... [is ] internally mixed with sea salt, ... [it]
             !    should be treated identically to SALC in the code for
             !    all processes.  (15 Dec 2015)

             Fullname = 'Inorganic nitrates on surface of seasalt aerosol'
             IF ( Present(oRadius) ) THEN
                IF ( .NOT. Present(Input_Opt) ) THEN
                   WRITE( 6, '(a)' ) REPEAT( '=', 79 )
                   WRITE( 6, * ) 'Error getting radius for species ',TRIM(Name)
                   WRITE( 6, * ) 'Input_Opt is missing!'
                   WRITE( 6, * ) 'In module Headers/species_database_mod.F90!'
                   RC = -1
                   RETURN
                ENDIF
                Radius     = ( Input_Opt%SALC_REDGE_um(1) +                  &
                               Input_Opt%SALC_REDGE_um(2)  ) * 0.5e-6_fp
             ENDIF

             ! Halve the Kc (cloud condensate -> precip) rate
             ! for the temperature range 237 K <= T < 258 K.
             KcScale       = (/ 1.0_fp, 0.5_fp, 1.0_fp /)

             ! Turn off rainout only when 237 K <= T < 258K.
             RainEff       = (/ 1.0_fp, 0.0_fp, 1.0_fp /)

             Formula       = ''
             MW_g          = 31.4_fp
             Is_Gas        = F
             Is_Drydep     = T
             Is_Wetdep     = T
             Density       = 2200.0_fp
             Radius        = Radius
             DD_AeroDryDep = T
             DD_F0         = 0.0_fp
             DD_Hstar_Old  = 0.0_fp
             WD_AerScavEff = 1.0_fp
             WD_KcScaleFac = KcScale
             WD_RainoutEff = RainEff

          CASE( 'NO' )
             FullName      = 'Nitrogen oxide'
             Formula       = 'NO'
             MW_g          = 30.0_fp
             BackgroundVV  = 4.0e-13_fp
             Is_Gas        = T
             Is_Drydep     = F
             Is_Wetdep     = F
             Is_Photolysis = T
#if defined( NEW_HENRY_CONSTANTS )                                          
             Henry_K0      = 1.90e-5_f8 * To_M_atm
             Henry_CR      = 1600.0_f8
#endif

          CASE( 'NO2' )
             FullName      = 'Nitrogen dioxide'
             Formula       = 'NO2'
             MW_g          = 46.0_fp
             BackgroundVV  = 4.0e-13_fp
             Is_Gas        = T
             Is_Drydep     = T
             Is_Wetdep     = F
             Is_Photolysis = T
             DD_F0         = 0.1_fp
#if defined( NEW_HENRY_CONSTANTS )
             Henry_K0      = 1.20e-4_f8 * To_M_atm
             Henry_CR      = 2400.0_f8
#else
             DD_Hstar_old  = 1.00e-2_fp
#endif

          CASE( 'NO3' )
             FullName      = 'Nitrate radical'
             Formula       = 'NO3'
             MW_g          = 62.0_fp
             BackgroundVV  = 4.0e-15_fp
             MolecRatio    = 1.0_fp
             Is_Gas        = T
             Is_Drydep     = F
             Is_Wetdep     = F
             Is_Photolysis = T
#if defined( NEW_HENRY_CONSTANTS )
             Henry_K0      = 3.80e-4_f8 * To_M_atm
             Henry_CR      = 1900.0_f8
#endif

          CASE( 'NPMN' )
             FullName      = 'Non-isoprene peroxymethacroyl nitrate (PMN)'

             ! PMN uses the same DD_F0 and DD_Hstar_old values as PAN
             ! so that we can compute its drydep velocity explicitly.
             ! (bmy, 5/19/16)
             Formula       = 'CH2=C(CH3)C(O)OONO2 '
             MW_g          = 147.0_fp
             Is_Gas        = T
             Is_Drydep     = T
             Is_Wetdep     = F
             DD_F0         = 1.0_fp
#if defined( NEW_HENRY_CONSTANTS )
             Henry_K0      = 1.70e-2_f8 * To_M_atm
#else
             DD_Hstar_old  = 3.60_fp
#endif

          CASE( 'O3',     'O3STRAT', 'O3UT',   'O3MT',   'O3ROW',           &
                'O3PCBL', 'O3NABL',  'O3ATBL', 'O3EUBL', 'O3AFBL',          &
                'O3ASBL', 'O3INIT',  'O3USA',  'O3STRT'            )

             ! Now include both total and tagged ozone species (bmy, 10/5/15)
             ! All of these have identical properties except for the names
             SELECT CASE( TRIM( Name ) )
                CASE( 'O3' )
                   FullName ='Ozone'
                CASE ( 'O3STRAT', 'O3STRT' )
                   FullName = 'Ozone produced in the stratosphere'
                CASE( 'O3UT' )
                   FullName = 'Ozone produced in the upper tropopshere'
                CASE( 'O3MT' )
                   FullName = 'Ozone produced in the middle troposphere'
                CASE( 'O3ROW' )
                   FullName = 'Ozone produced in the rest of the world'
                CASE( 'O3PCBL' )
                   FullName = 'Ozone produced in the Pacific Ocean boundary layer'
                CASE( 'O3NABL' )
                   FullName = 'Ozone produced in the North American boundary layer'
                CASE( 'O3ATBL' )
                   FullName = 'Ozone produced in the Atlantic Ocean boundary layer'
                CASE( 'O3EUBL' )
                   FullName = 'Ozone produced in the European boundary layer'
                CASE( 'O3AFBL' )
                   FullName = 'Ozone produced in the African boundary layer'
                CASE( 'O3ASBL' )
                   FullName = 'Ozone produced in the Asian boundary layer'
                CASE( 'O3INIT' )
                   FullName = 'Ozone from the initial condition'
                CASE( 'O3USA' )
                   FullName = 'Ozone produced over the United States in PBL'
             END SELECT

             ! Set special default background for O3
             SELECT CASE( TRIM( Name ) )
                CASE( 'O3' )
                   BackgroundVV = 2.0e-08_fp
                CASE DEFAULT
                   BackgroundVV = MISSING_VV
             END SELECT

             Formula       = 'O3'
             MW_g          = 48.0_fp
             Is_Gas        = T
             Is_Drydep     = T
             Is_Wetdep     = F
             Is_Photolysis = T
             DD_F0         = 1.0_fp
#if defined( NEW_HENRY_CONSTANTS )
             Henry_K0      = 1.90e-05_f8 * To_M_atm
             Henry_CR      = 1600.0_f8
#else
             DD_Hstar_old  = 1.0e-2_fp
#endif

          CASE( 'OCLO' )
             FullName      = 'Chlorine dioxide'
             Formula       = 'OClO'
             MW_g          = 67.0_fp
             Is_Gas        = T
             Is_Drydep     = F
             Is_Wetdep     = F
             Is_Photolysis = T

          CASE( 'OCPI' )
             FullName      = 'Hydrophilic organic carbon aerosol'

             ! Halve the Kc (cloud condensate -> precip) rate
             ! for the temperature range 237 K <= T < 258 K.
             KcScale       = (/ 1.0_fp, 0.5_fp, 1.0_fp /)

             ! Turn off rainout only when 237 K <= T < 258K.
             RainEff       = (/ 1.0_fp, 0.0_fp, 1.0_fp /)

             Formula       = ''
             MW_g          = 12.01_fp
             EmMW_g        = 12.0_fp
             Is_Gas        = F
             Is_Drydep     = T
             Is_Wetdep     = T
             Is_HygroGrowth= T
             Density       = 1300.0_fp
             DD_DvzAerSnow = 0.03_fp
             DD_F0         = 0.0_fp
             DD_Hstar_Old  = 0.0_fp
             WD_AerScavEff = 1.0_fp
             WD_KcScaleFac = KcScale
             WD_RainoutEff = RainEff

          CASE( 'OCPO' )
             Fullname      = 'Hydrophobic organic carbon aerosol'

             ! For all temperatures:
             ! (1) Halve the Kc (cloud condensate -> precip) rate
             ! (2) Turn off rainout (OCPO is hydrophobic)
             KcScale       = (/ 0.5_fp, 0.5_fp, 0.5_fp /)
             RainEff       = (/ 0.0_fp, 0.0_fp, 0.0_fp /)

             Formula       = ''
             MW_g          = 12.01_fp
             EmMW_g        = 12.0_fp
             Is_Gas        = F
             Is_Drydep     = T
             Is_Wetdep     = T
             Is_HygroGrowth= F
             Density       = 1300.0_fp
             DD_DvzAerSnow = 0.03_fp
             DD_F0         = 0.0_fp
             DD_Hstar_Old  = 0.0_fp
             WD_AerScavEff = 1.0_fp
             WD_KcScaleFac = KcScale
             WD_RainoutEff = RainEff

          CASE( 'OCS' )
             FullName      = 'Carbonyl sulfide'
             Formula       = 'COS'
             MW_g          = 60.0_fp
             BackgroundVV  = 9.0e-15_fp
             Is_Gas        = T
             Is_Drydep     = F
             Is_Wetdep     = F
             Is_Photolysis = T
             
          CASE( 'OPOA1', 'OPOA2' )
             FullName      = 'Lumped aerosol product of SVOC oxidation'

             ! Halve the Kc (cloud condensate -> precip) rate
             ! for the temperature range 237 K <= T < 258 K.
             KcScale       = (/ 1.0_fp, 0.5_fp, 1.0_fp /)

             ! Turn off rainout only when 237 K <= T < 258K.
             ! NOTE: Rainout efficiency is 0.8 because these are SOA species.
             RainEff       = (/ 0.8_fp, 0.0_fp, 0.8_fp /)

             Formula       = ''
             MW_g          = 12.01_fp
             EmMW_g        = 12.0_fp
             MolecRatio    = 1.0_fp
             Is_Gas        = F
             Is_Drydep     = T
             Is_Wetdep     = T
             DD_DvzAerSnow = 0.03_fp
             DD_F0         = 0.0_fp
             DD_Hstar_old  = 0.0_fp
             WD_AerScavEff = 0.8_fp
             WD_KcScaleFac = KcScale
             WD_RainoutEff = RainEff

          CASE( 'OPOG1', 'OPOG2' )
             FullName      = 'Lumped gas product of SVOC oxidation'
             Formula       = ''
             MW_g          = 12.01_fp
             EmMW_g        = 12.0_fp
             MolecRatio    = 1.0_fp
             Is_Gas        = T
             Is_Drydep     = T
             Is_Wetdep     = T
             DD_F0         = 0.0_fp
             DD_Hstar_old  = 1.00e+5_fp
             Henry_K0      = 1.00e+5_f8
             Henry_CR      = 6039.0_f8
             WD_RetFactor  = 2.0e-2_fp

          CASE( 'PAN' )
             FullName      = 'Peroxyacetyl nitrate'
             Formula       = 'CH3C(O)OONO2'
             MW_g          = 121.0_fp
             Is_Gas        = T
             Is_Drydep     = T
             Is_Wetdep     = F
             Is_Photolysis = T
             DD_F0         = 1.0e+0_fp
#if defined( NEW_HENRY_CONSTANTS )
             Henry_K0      = 2.90e+02_f8 * To_M_atm
             Henry_CR      = 5700.0_f8
#else
             DD_Hstar_old  = 3.60e+0_fp
#endif

          CASE( 'PFE' )

             ! Halve the Kc (cloud condensate -> precip) rate
             ! for the temperature range 237 K <= T < 258 K.
             KcScale       = (/ 1.0_fp, 0.5_fp, 1.0_fp /)

             ! Turn off rainout only when 237 K <= T < 258K.
             RainEff       = (/ 1.0_fp, 0.0_fp, 1.0_fp /)   

             FullName      = 'Anthropogenic iron'
             Formula       = 'Fe'
             MW_g          = 55.85_fp
             Is_Gas        = F
             Is_Drydep     = T
             Is_Wetdep     = T
             DD_F0         = 0.0e+0_fp
             DD_Hstar_old  = 0.0e+0_fp
             WD_AerScavEff = 1.0_fp
             WD_KcScaleFac = KcScale
             WD_RainoutEff = RainEff

          CASE( 'POA1' )

             ! For all temperatures:
             ! (1) Halve the Kc (cloud condensate -> precip) rate
             ! (2) Turn off rainout (these are hydrophobic species)
             KcScale       = (/ 0.5_fp, 0.5_fp, 0.5_fp /)
             RainEff       = (/ 0.0_fp, 0.0_fp, 0.0_fp /)

             FullName      ='Lumped aerosol primary SVOCs'
             Formula       = ''
             MW_g          = 12.01_fp
             EmMW_g        = 12.0_fp
             MolecRatio    = 1.0_fp
             Is_Gas        = F
             Is_Drydep     = T
             Is_Wetdep     = T
             Is_HygroGrowth= T
             Density       = 1300.0_fp
             DD_DvzAerSnow = 0.03_fp
             DD_F0         = 0.0_fp
             DD_Hstar_old  = 0.0_fp
             WD_AerScavEff = 1.0_fp
             WD_KcScaleFac = KcScale
             WD_RainoutEff = RainEff

          CASE( 'POA2' )

             ! For all temperatures:
             ! (1) Halve the Kc (cloud condensate -> precip) rate
             ! (2) Turn off rainout (these are hydrophobic species)
             KcScale       = (/ 0.5_fp, 0.5_fp, 0.5_fp /)
             RainEff       = (/ 0.0_fp, 0.0_fp, 0.0_fp /)

             FullName      ='Lumped aerosol primary SVOCs'
             Formula       = ''
             MW_g          = 12.01_fp
             EmMW_g        = 12.0_fp
             MolecRatio    = 1.0_fp
             Is_Gas        = F
             Is_Drydep     = T
             Is_Wetdep     = T
             Density       = 1300.0_fp
             DD_DvzAerSnow = 0.03_fp
             DD_F0         = 0.0_fp
             DD_Hstar_old  = 0.0_fp
             WD_AerScavEff = 1.0_fp
             WD_KcScaleFac = KcScale
             WD_RainoutEff = RainEff

          CASE( 'POG1', 'POG2' )
             FullName      = 'Lumped gas primary SVOCs'
             Formula       = ''
             MW_g          = 12.01_fp
             EmMW_g        = 12.0_fp
             MolecRatio    = 1.0_fp
             Is_Gas        = T
             Is_Drydep     = T
             Is_Wetdep     = T
             DD_F0         = 0.0_fp
             DD_Hstar_old  = 9.50e+0_fp
             Henry_K0      = 9.50e+0_f8
             Henry_CR      = 4700.0_f8
             WD_RetFactor  = 2.0e-2_fp

          CASE( 'PPN' )
             ! PPN uses the same DD_F0 and DD_Hstar_old values as PAN
             ! so that we can compute its drydep velocity explicitly.
             ! (bmy, 5/19/16)
             FullName = 'Lumped peroxypropionyl nitrate'
             Formula       = 'CH3CH2C(O)OONO2'
             MW_g          = 135.0_fp
             Is_Gas        = T
             Is_Drydep     = T
             Is_Wetdep     = F
             DD_F0         = 1.0_fp
#if defined( NEW_HENRY_CONSTANTS )
             Henry_K0      = 2.9e-2_fp * To_M_atm
             Henry_CR      = 0.0_f8
#else
             DD_Hstar_old  = 3.60_fp
#endif

          CASE( 'PRPE' )
             FullName      = 'Lumped >= C3 alkenes'
             Formula       = 'C3H6'
             MW_g          = 42.08_fp
             EmMW_g        = 12.0_fp
             MolecRatio    = 3.0_fp
             Is_Gas        = T
             Is_Drydep     = F
             Is_Wetdep     = F
#if defined( NEW_HENRY_CONSTANTS )
             Henry_K0      = 7.3e-5_f8 * To_M_atm
             Henry_CR      = 3400.0_f8
#endif

          CASE( 'PROPNN' )
             FullName      = 'Propanone nitrate'
             Formula       = 'CH3C(=O)CH2ONO2'
             MW_g          = 119.0_fp
             Is_Gas        = T
             Is_Drydep     = T
             Is_Wetdep     = T
             Is_Photolysis = T
             DD_F0         = 1.0_fp
#if defined( NEW_HENRY_CONSTANTS )
             Henry_K0      = 4.93e+3_f8 * To_M_atm
             Henry_CR      = 0.0_f8
#else
             DD_Hstar_old  = 5.00e+5_fp
             Henry_K0      = 5.00e+5_f8
             Henry_CR      = 0.0_f8
#endif
             WD_RetFactor  = 2.0e-2_fp

          CASE( 'R4N2' )
             ! Now drydep is like other alkyl nitrates (Ito et al., 2007)
             ! (eam, 2014)
             FullName      = 'Lumped alkyl nitrate'
             Formula       = 'RO2NO'
             MW_g          = 119.0_fp
             Is_Gas        = T
             Is_Drydep     = T
             Is_Wetdep     = F
             Is_Photolysis = T
             DD_F0         = 1.0_fp
#if defined( NEW_HENRY_CONSTANTS )
             Henry_K0      = 1.0e-2_f8 * To_M_atm
             Henry_CR      = 5800.0_f8
#else
             DD_Hstar_old  = 1.70e+4_fp
#endif

          CASE( 'RCHO' )
             FullName      = 'Lumped aldehyde >= C3'
             Formula       = 'CH3CH2CHO'
             MW_g          = 58.0_fp
             Is_Gas        = T
             Is_Drydep     = F
             Is_Wetdep     = F
             Is_Photolysis = T
#if defined( NEW_HENRY_CONSTANTS )
             Henry_K0      = 9.5e-2_f8 * To_M_atm
             Henry_CR      = 6200.0_f8
#else
             Henry_K0      = 4.20e+3_f8
             Henry_CR      = 0.0_f8
#endif
             WD_RetFactor  = 2.0e-2_fp

          CASE( 'RIPA' )
             FullName      = '1,2-ISOPOOH'
             Formula       = 'C5H10O3'
             MW_g          = 118.0_fp
             Is_Gas        = T
             Is_Drydep     = T
             Is_Wetdep     = T
             Is_Photolysis = T
             DD_F0         = 1.0_fp
             DD_Hstar_old  = 1.70e+6_fp
             Henry_K0      = 1.70e+6_f8
             Henry_CR      = 0.0_f8
             WD_RetFactor  = 2.0e-2_fp

          CASE( 'RIPB' )
             FullName      = '4,3-ISOPOOH'
             Formula       = 'C5H10O3'
             MW_g          = 118.0_fp
             Is_Gas        = T
             Is_Drydep     = T
             Is_Wetdep     = T
             Is_Photolysis = T
             DD_F0         = 1.0_fp
             DD_Hstar_old  = 1.70e+6_fp
             Henry_K0      = 1.70e+6_f8
             Henry_CR      = 0.0_f8
             WD_RetFactor  = 2.0e-2_fp

          CASE( 'RIPD' )
             FullName      = 'd(1,4 and 4,1)-ISOPOOH'
             Formula       = 'C5H10O3'
             MW_g          = 118.0_fp
             Is_Gas        = T
             Is_Drydep     = T
             Is_Wetdep     = T
             Is_Photolysis = T
             DD_F0         = 1.0_fp
             DD_Hstar_old  = 1.70e+6_fp
             Henry_K0      = 1.70e+6_f8
             Henry_CR      = 0.0_f8
             WD_RetFactor  = 2.0e-2_fp

          CASE( 'SALA' )
             FullName = 'Accumulation mode sea salt aerosol'
             IF ( Present(oRadius) ) THEN
                IF ( .NOT. Present(Input_Opt) ) THEN
                   WRITE( 6, '(a)' ) REPEAT( '=', 79 )
                   WRITE( 6, * ) 'Error getting radius for species ',TRIM(Name)
                   WRITE( 6, * ) 'Input_Opt is missing!'
                   WRITE( 6, * ) 'In module Headers/species_database_mod.F90!'
                   RC = -1
                   RETURN
                ENDIF
                Radius     = ( Input_Opt%SALA_REDGE_um(1) +                  & 
                               Input_Opt%SALA_REDGE_um(2)  ) * 0.5e-6_fp
             ENDIF

             ! Halve the Kc (cloud condensate -> precip) rate
             ! for the temperature range 237 K <= T < 258 K.
             KcScale       = (/ 1.0_fp, 0.5_fp, 1.0_fp /)

             ! Turn off rainout only when 237 K <= T < 258K.
             RainEff       = (/ 1.0_fp, 0.0_fp, 1.0_fp /)

             Formula       = ''
             MW_g          = 31.4_fp
             Is_Gas        = F
             Is_Drydep     = T
             Is_Wetdep     = T
             Is_HygroGrowth= T
             Density       = 2200.0_fp
             Radius        = Radius
             DD_AeroDryDep = T
             DD_F0         = 0.0_fp
             DD_Hstar_old  = 0.0_fp
             WD_AerScavEff = 1.0_fp
             WD_KcScaleFac = KcScale
             WD_RainoutEff = RainEff

          CASE( 'SALC' )
             FullName = 'Coarse mode sea salt aerosol'
             IF ( Present(oRadius) ) THEN
                IF ( .NOT. Present(Input_Opt) ) THEN
                   WRITE( 6, '(a)' ) REPEAT( '=', 79 )
                   WRITE( 6, * ) 'Error getting radius for species ',TRIM(Name)
                   WRITE( 6, * ) 'Input_Opt is missing!'
                   WRITE( 6, * ) 'In module Headers/species_database_mod.F90!'
                   RC = -1
                   RETURN
                ENDIF
                Radius     = ( Input_Opt%SALC_REDGE_um(1) +                  &
                               Input_Opt%SALC_REDGE_um(2)  ) * 0.5e-6_fp
             ENDIF

             ! Halve the Kc (cloud condensate -> precip) rate
             ! for the temperature range 237 K <= T < 258 K.
             KcScale       = (/ 1.0_fp, 0.5_fp, 1.0_fp /)

             ! Turn off rainout only when 237 K <= T < 258K.
             RainEff       = (/ 1.0_fp, 0.0_fp, 1.0_fp /)

             FullName      = Fullname
             Formula       = ''
             MW_g          = 31.4_fp
             Is_Gas        = F
             Is_Drydep     = T
             Is_Wetdep     = T
             Is_HygroGrowth= T
             Density       = 2200.0_fp
             Radius        = Radius
             DD_AeroDryDep = T
             DD_F0         = 0.0_fp
             DD_Hstar_old  = 0.0_fp
             WD_AerScavEff = 1.0_fp
             WD_CoarseAer  = T
             WD_KcScaleFac = KcScale
             WD_RainoutEff = RainEff

          CASE( 'SO2' )

             !%%% NOTE: SO2 dry-deposits like a gas but wet-deposits
             !%%% like an aerosol.  Therefore, we need to define SO2 with
             !%%% both gas-phase and aerosol parameters. (bmy, 9/28/15)

             ! Halve the Kc (cloud condensate -> precip) rate
             ! for the temperature range 237 K <= T < 258 K.
             KcScale       = (/ 1.0_fp, 0.5_fp, 1.0_fp /)

             ! Turn off rainout only when 237 K <= T < 258K.
             RainEff       = (/ 1.0_fp, 0.0_fp, 1.0_fp /)

             ! Enforce minimum dry deposition velocity (Vd) for SO2
             ! (cf. Mian Chin's GOCART model)
             ! Minimum Vd over snow/ice : 0.2 cm/s
             ! Minimum Vd over land     : 0.3 cm/s
             DvzMinVal     = (/ 0.2_fp, 0.3_fp /)

             FullName      = 'Sulfur dioxide'
             Formula       = 'SO2'
             MW_g          = 64.0_fp
             Is_Gas        = T
             Is_Drydep     = T
             Is_Wetdep     = T
             DD_DvzAerSnow = 0.03_fp
             DD_DvzMinVal  = DvzMinVal
             DD_F0         = 0.0_fp
#if defined( NEW_HENRY_CONSTANTS )
             Henry_K0      = 1.30e-2_f8 * To_M_atm
             Henry_CR      = 2900.0_f8
#else
             DD_Hstar_old  = 1.00e+5_fp
#endif
             WD_AerScavEff = 1.0_fp
             WD_KcScaleFac = KcScale
             WD_RainoutEff = RainEff

          CASE( 'SO4' )

             ! Halve the Kc (cloud condensate -> precip) rate
             ! for the temperature range 237 K <= T < 258 K.
             KcScale       = (/ 1.0_fp, 0.5_fp, 1.0_fp /)

             ! Turn off rainout only when 237 K <= T < 258K.
             RainEff       = (/ 1.0_fp, 0.0_fp, 1.0_fp /)

             ! Enforce minimum dry deposition velocity (Vd) for SO4
             ! (cf. Mian Chin's GOCART model)
             ! Minimum Vd over snow/ice : 0.01 cm/s
             ! Minimum Vd over land     : 0.01 cm/s
             DvzMinVal     = (/ 0.01_fp, 0.01_fp /)

             FullName      = 'Sulfate'
             Formula       = 'SO4'
             MW_g          = 96.0_fp
             Is_Gas        = F
             Is_Drydep     = T
             Is_Wetdep     = T
             Is_Photolysis = T
             Is_HygroGrowth= T
             Density       = 1700.0_fp
             DD_DvzAerSnow = 0.03_fp
             DD_DvzMinVal  = DvzMinVal
             DD_F0         = 0.0_fp
             DD_Hstar_Old  = 0.0_fp
             WD_AerScavEff = 1.0_fp
             WD_KcScaleFac = KcScale
             WD_RainoutEff = RainEff

          CASE( 'SO4S' )

             ! NOTE: SO4s (SO4 on coarse sea salt aerosol) needs to use
             ! the same molecular weight as coarse sea salt (= 31.4 g/mol)
             ! instead of SO4 (= 96 g/mol).
             !
             ! Becky Alexander (beckya@atmos.washington.edu) wrote:
             !
             !   "The reason for using sea salt's MW for SO4s is that ...
             !    [it is] ... essentially internally mixed with coarse sea
             !    salt aerosol (SALC).  As coarse sea salt aerosol likely
             !    dominates the mass of ... [SO4s] ...  it is appropriate
             !    to use sea salt's MW.  Another explanation is that since
             !    SO4s ... [is ] internally mixed with sea salt, ... [it]
             !    should be treated identically to SALC in the code for
             !    all processes.  (15 Dec 2015)
             !
             Fullname = 'Sulfate on surface of seasalt aerosol'
             IF ( Present(oRadius) ) THEN
                IF ( .NOT. Present(Input_Opt) ) THEN
                   WRITE( 6, '(a)' ) REPEAT( '=', 79 )
                   WRITE( 6, * ) 'Error getting radius for species ',TRIM(Name)
                   WRITE( 6, * ) 'Input_Opt is missing!'
                   WRITE( 6, * ) 'In module Headers/species_database_mod.F90!'
                   RC = -1
                   RETURN
                ENDIF
                Radius     = ( Input_Opt%SALC_REDGE_um(1) +                  &
                               Input_Opt%SALC_REDGE_um(2)  ) * 0.5e-6_fp
             ENDIF

             ! Halve the Kc (cloud condensate -> precip) rate
             ! for the temperature range 237 K <= T < 258 K.
             KcScale       = (/ 1.0_fp, 0.5_fp, 1.0_fp /)

             ! Turn off rainout only when 237 K <= T < 258K.
             RainEff       = (/ 1.0_fp, 0.0_fp, 1.0_fp /)

             Formula       = ''
             MW_g          = 31.4_fp
             Is_Gas        = F
             Is_Drydep     = T
             Is_Wetdep     = T
             Density       = 2200.0_fp
             Radius        = Radius
             DD_AeroDryDep = T
             DD_F0         = 0.0_fp
             DD_Hstar_Old  = 0.0_fp
             WD_AerScavEff = 1.0_fp
             WD_KcScaleFac = KcScale
             WD_RainoutEff = RainEff

          CASE( 'SOAIE' )

             ! Halve the Kc (cloud condensate -> precip) rate
             ! for the temperature range 237 K <= T < 258 K.
             KcScale       = (/ 1.0_fp, 0.5_fp, 1.0_fp /)

             ! Turn off rainout only when 237 K <= T < 258K.
             ! NOTE: Rainout efficiency is 0.8 because these are SOA species.
             RainEff       = (/ 0.8_fp, 0.0_fp, 0.8_fp /)

             FullName      = 'Aerosol-phase IEPOX'
             Formula       = 'C5H10O3'
             MW_g          = 118.0_fp
             Is_Gas        = F
             Is_Drydep     = T
             Is_Wetdep     = T
             DD_DvzAerSnow = 0.03_fp
             DD_F0         = 0.0_fp
             DD_HStar_old  = 0.0_fp
             WD_AerScavEff = 0.8_fp
             WD_KcScaleFac = KcScale
             WD_RainoutEff = RainEff

          CASE( 'SOAGX' )

             ! Halve the Kc (cloud condensate -> precip) rate
             ! for the temperature range 237 K <= T < 258 K.
             KcScale       = (/ 1.0_fp, 0.5_fp, 1.0_fp /)

             ! Turn off rainout only when 237 K <= T < 258K.
             ! NOTE: Rainout efficiency is 0.8 because these are SOA species.
             RainEff       = (/ 0.8_fp, 0.0_fp, 0.8_fp /)

             FullName      = 'Aerosol-phase glyoxal'
             Formula       = 'C2H2O2'
             MW_g          = 58.0_fp
             Is_Gas        = F
             Is_Drydep     = T
             Is_Wetdep     = T
             DD_DvzAerSnow = 0.03_fp
             DD_F0         = 0.0_fp
             DD_HStar_old  = 0.0_fp
             WD_AerScavEff = 0.8_fp
             WD_KcScaleFac = KcScale
             WD_RainoutEff = RainEff

          CASE( 'SOAME' )

             ! Halve the Kc (cloud condensate -> precip) rate
             ! for the temperature range 237 K <= T < 258 K.
             KcScale       = (/ 1.0_fp, 0.5_fp, 1.0_fp /)

             ! Turn off rainout only when 237 K <= T < 258K.
             ! NOTE: Rainout efficiency is 0.8 because these are SOA species.
             RainEff       = (/ 0.8_fp, 0.0_fp, 0.8_fp /)

             FullName      = 'Aerosol-phase IMAE'
             Formula       = 'C4H6O3'
             MW_g          = 102.0_fp
             Is_Gas        = F
             Is_Drydep     = T
             Is_Wetdep     = T
             DD_DvzAerSnow = 0.03_fp
             DD_F0         = 0.0_fp
             DD_HStar_old  = 0.0_fp
             WD_AerScavEff = 0.8_fp
             WD_KcScaleFac = KcScale
             WD_RainoutEff = RainEff

          CASE( 'SOAMG' )

             ! Halve the Kc (cloud condensate -> precip) rate
             ! for the temperature range 237 K <= T < 258 K.
             KcScale       = (/ 1.0_fp, 0.5_fp, 1.0_fp /)

             ! Turn off rainout only when 237 K <= T < 258K.
             ! NOTE: Rainout efficiency is 0.8 because these are SOA species.
             RainEff       = (/ 0.8_fp, 0.0_fp, 0.8_fp /)

             FullName      = 'Aerosol-phase methylglyoxal'
             Formula       = 'C3H4O2'
             MW_g          = 72.0_fp
             Is_Gas        = F
             Is_Drydep     = T
             Is_Wetdep     = T
             DD_DvzAerSnow = 0.03_fp
             DD_F0         = 0.0_fp
             DD_HStar_old  = 0.0_fp
             WD_AerScavEff = 0.8_fp
             WD_KcScaleFac = KcScale
             WD_RainoutEff = RainEff

          CASE( 'SOAP' )
             FullName      = 'SOA Precursor - lumped species for simplified SOA paramterization'

             !SOAPis not removed because it is a simple parameterization,
             !not a physical model
             FullName      = FullName
             Formula       = ''
             MW_g          = 150.0_fp
             Is_Gas        = T
             Is_Drydep     = F
             Is_Wetdep     = F

          CASE( 'SOAS' )
             FullName = 'SOA Simple - simplified non-volatile SOA parameterization'
             !Copy data from ISOA

             ! Halve the Kc (cloud condensate -> precip) rate
             ! for the temperature range 237 K <= T < 258 K.
             KcScale       = (/ 1.0_fp, 0.5_fp, 1.0_fp /)

             ! Turn off rainout only when 237 K <= T < 258K.
             ! NOTE: Rainout efficiency is 0.8 because these are SOA species.
             RainEff       = (/ 0.8_fp, 0.0_fp, 0.8_fp /)

             FullName      = FullName
             Formula       = ''
             MW_g          = 150.0_fp
             Is_Gas        = F
             Is_Drydep     = T
             Is_Wetdep     = T
             DD_DvzAerSnow = 0.03_fp
             DD_F0         = 0.0_fp
             DD_HStar_old  = 0.0_fp
             WD_AerScavEff = 0.8_fp
             WD_KcScaleFac = KcScale
             WD_RainoutEff = RainEff

          CASE( 'TOLU' )
             FullName      = 'Toluene'
             Formula       = 'C7H8'
             MW_g          = 92.14_fp
             EmMW_g        = 12.0_fp
             MolecRatio    = 7.0_fp
             Is_Gas        = T
             Is_Drydep     = F
             Is_Wetdep     = F

          CASE( 'TSOA0', 'TSOA1', 'TSOA2', 'TSOA3' )
             FullName = 'Lumped semivolatile aerosol products of monoterpene + sesquiterpene oxidation'

             ! Halve the Kc (cloud condensate -> precip) rate
             ! for the temperature range 237 K <= T < 258 K.
             KcScale       = (/ 1.0_fp, 0.5_fp, 1.0_fp /)

             ! Turn off rainout only when 237 K <= T < 258K.
             ! NOTE: Rainout efficiency is 0.8 because these are SOA species.
             RainEff       = (/ 0.8_fp, 0.0_fp, 0.8_fp /)

             FullName      = FullName
             Formula       = ''
             MW_g          = 150.0_fp
             Is_Gas        = F
             Is_Drydep     = T
             Is_Wetdep     = T
             DD_DvzAerSnow = 0.03_fp
             DD_F0         = 0.0_fp
             DD_Hstar_old  = 0.0_fp
             WD_AerScavEff = 0.8_fp
             WD_KcScaleFac = KcScale
             WD_RainoutEff = RainEff

          CASE( 'TSOG0', 'TSOG1', 'TSOG2', 'TSOG3' )
             FullName      = 'Lumped semivolatile gas products of monoterpene + sesquiterpene oxidation'
             Formula       = ''
             MW_g          = 150.0_fp
             Is_Gas        = T
             Is_Drydep     = T
             Is_Wetdep     = T
             DD_F0         = 0.0_fp
             DD_Hstar_old  = 1.00e+5_fp
             Henry_K0      = 1.00e+5_f8
             Henry_CR      = 6039.0_f8
             WD_RetFactor  = 2.0e-2_fp

          CASE( 'XYLE' )
             FullName      = 'Xylene'
             Formula       = 'C8H10'
             MW_g          = 106.16_fp
             EmMW_g        = 12.0_fp
             MolecRatio    = 8.0_fp
             Is_Gas        = T
             Is_Drydep     = F
             Is_Wetdep     = F

          !==================================================================
          ! Species for simulations including iodine
          !==================================================================

          CASE( 'CH3I' )
             FullName      = 'Methyl iodide'
             Formula       = 'CH3I'
             MW_g          = 142.0_fp
             Is_Gas        = T
             Is_Drydep     = F
             Is_Wetdep     = F
             Is_Photolysis = T
             DD_F0         = 0.0_fp
             DD_Hstar_old  = 0.0_fp
             Henry_K0      = 2.0e-3_f8 * To_M_atm
             Henry_CR      = 3.6e+3_f8

          CASE( 'CH2I2' )
             FullName      = 'Diiodomethane'
             Formula       = 'CH2I2'
             MW_g          = 268.0_fp
             Is_Gas        = T
             Is_Drydep     = F
             Is_Wetdep     = F
             Is_Photolysis = T

          CASE( 'CH2ICL' )
             FullName      = 'Chloroiodomethane'
             Formula       = 'CH2ICl'
             MW_g          = 167.0_fp
             Is_Gas        = T
             Is_Drydep     = F
             Is_Wetdep     = F
             Is_Photolysis = T

          CASE( 'CH2IBR' )
             FullName      = 'Bromoiodomethane'
             Formula       = 'CH2IBr'
             MW_g          = 221.0_fp
             Is_Gas        = T
             Is_Drydep     = F
             Is_Wetdep     = F
             Is_Photolysis = T

          CASE( 'HOI' )
             FullName      = 'Hypoiodous acid'
             Formula       = 'HOI'
             MW_g          = 144.0_fp
             Is_Gas        = T
             Is_Drydep     = T
             Is_Wetdep     = T
             Is_Photolysis = T
             DD_F0         = 0.0_fp
#if defined( NEW_HENRY_CONSTANTS )
             Henry_K0      = -1.0e+0_f8 * To_M_atm
             Henry_CR      = -1.0e+0_f8
#else
             DD_Hstar_old  = 1.54e+4_f8
             Henry_K0      = 1.54e+4_f8
             Henry_CR      = 8.371e+3_f8
#endif
             WD_RetFactor  = 2.0e-2_fp

          CASE( 'I2' )
             FullName      = 'Molecular iodine'
             Formula       = 'I2'
             MW_g          = 254.0_fp
             Is_Gas        = T
             Is_Drydep     = T
             Is_Wetdep     = T
             Is_Photolysis = T
             DD_F0         = 0.0_fp
#if defined( NEW_HENRY_CONSTANTS )
             Henry_K0      = -1.0e+0_f8 * To_M_atm
             Henry_CR      = -1.0e+0_f8
#else
             DD_Hstar_old  = 2.70e+0_fp
             Henry_K0      = 2.7e+0_f8
             Henry_CR      = 7.5074e+3_f8
#endif
             WD_RetFactor  = 0.0_fp

          CASE( 'IBR' )
             FullName      = 'Iodine monobromide'
             Formula       = 'IBr'
             MW_g          = 207.0_fp
             Is_Gas        = T
             Is_Drydep     = T
             Is_Wetdep     = T
             Is_Photolysis = T
             DD_F0         = 0.0_fp
#if defined( NEW_HENRY_CONSTANTS )
             Henry_K0      = -1.0e+0_f8 * To_M_atm
             Henry_CR      = -1.0e+0_f8
#else
             DD_Hstar_old  = 2.43e+1_fp
             Henry_K0      = 2.4e+1_f8
             Henry_CR      = 4.9167e+3_f8
#endif
             WD_RetFactor  = 0.0_fp

          CASE( 'ICL' )
             FullName      = 'Iodine monochloride'
             Formula       = 'ICl'
             MW_g          = 162.0_fp
             Is_Gas        = T
             Is_Drydep     = T
             Is_Wetdep     = T
             Is_Photolysis = T
             DD_F0         = 0.0_fp
#if defined( NEW_HENRY_CONSTANTS )
             Henry_K0      = -1.0e+0_f8 * To_M_atm
             Henry_CR      = -1.0e+0_f8
#else
             DD_Hstar_old  = 1.11e+2_f8
             Henry_K0      = 1.11e+2_f8
             Henry_CR      = 2.1055e+3_f8
#endif
             WD_RetFactor  = 0.0_fp

          CASE( 'I' )
             FullName      = 'Atomic iodine'
             Formula       = 'I'
             MW_g          = 127.0_fp
             Is_Gas        = T
             Is_Drydep     = F
             Is_Wetdep     = F

          CASE( 'IO' )
             FullName      = 'Iodine monoxide'
             Formula       = 'IO'
             MW_g          = 143.0_fp
             Is_Gas        = T
             Is_Drydep     = F
             Is_Wetdep     = F
             Is_Photolysis = T

          CASE( 'HI' )
             FullName      = 'Hydrogen iodide'
             MW_g          = 128.0_fp
             Formula       = 'HI'
             Is_Gas        = T
             Is_Drydep     = T
             Is_Wetdep     = T
             DD_F0         = 0.0_fp
#if defined( NEW_HENRY_CONSTANTS )
             Henry_K0      = -1.0e+0_f8 * To_M_atm
             Henry_CR      = -1.0e+0_f8
#else
             DD_Hstar_old  = 2.35e+16_fp
             Henry_K0      = 7.43e+13_f8
             Henry_CR      = 3.1872e+3_f8
#endif
                              WD_RetFactor  = 1.0_fp

          CASE( 'OIO' )
             FullName      = 'Iodine dioxide'
             Formula       = 'OIO'
             MW_g          = 159.0_fp
             Is_Gas        = T
             Is_Drydep     = F
             Is_Wetdep     = F
             Is_Photolysis = T

          CASE( 'INO' )
             FullName      = 'Nitrosyl iodide'
             Formula       = 'INO'
             MW_g          = 157.0_fp
             Is_Gas        = T
             Is_Drydep     = F
             Is_Wetdep     = F
             Is_Photolysis = T

          CASE( 'IONO' )
             FullName      = 'Nitryl iodide'
             Formula       = 'IONO'
             MW_g          = 173.0_fp
             Is_Gas        = T
             Is_Drydep     = T
             Is_Wetdep     = T
             Is_Photolysis = T
             DD_F0         = 0.0_fp
#if defined( NEW_HENRY_CONSTANTS )
             Henry_K0      = -1.0e+0_f8 * To_M_atm
             Henry_CR      = -1.0e+0_f8
#else
             DD_Hstar_old  = 3.0e-1_f8
             Henry_K0      = 3.0e-1_f8
             Henry_CR      = 7.2404e+3_f8
#endif
             WD_RetFactor  = 2.0e-2_fp

          CASE( 'IONO2' )
             FullName      = 'Iodine nitrate'
             Formula       = 'IONO2'
             MW_g          = 189.0_fp
             Is_Gas        = T
             Is_Drydep     = T
             Is_Wetdep     = T
             Is_Photolysis = T
             DD_F0         = 0.0_fp
#if defined( NEW_HENRY_CONSTANTS )
             Henry_K0      = -1.0e+0_f8 * To_M_atm
             Henry_CR      = -1.0e+0_f8
#else
             DD_Hstar_old  = 1.0e+20_f8
             Henry_K0      = 1.0e+20_f8
             Henry_CR      = 3.98e+3_f8
#endif
             WD_RetFactor  = 1.0_fp

          CASE( 'I2O2' )
             FullName      = 'Diiodine dioxide'
             Formula       = 'I2O2'
             MW_g          = 286.0_fp
             Is_Gas        = T
             Is_Drydep     = T
             Is_Wetdep     = T
             Is_Photolysis = T
             DD_F0         = 0.0_fp
#if defined( NEW_HENRY_CONSTANTS )
             Henry_K0      = -1.0e+0_f8 * To_M_atm
             Henry_CR      = -1.0e+0_f8
#else
             DD_Hstar_old  = 1.0e+20_f8
             Henry_K0      = 1.0e+20_f8
             Henry_CR      = 1.89e+4_f8
#endif
             WD_RetFactor  = 1.0_fp

          CASE( 'I2O3' )
             FullName      = 'Diiodine trioxide'
             Formula       = 'I2O3'
             MW_g          = 302.0_fp
             Is_Gas        = T
             Is_Drydep     = T
             Is_Wetdep     = T
             Is_Photolysis = T
             DD_F0         = 0.0_fp
#if defined( NEW_HENRY_CONSTANTS )
             Henry_K0      = -1.0e+0_f8 * To_M_atm
             Henry_CR      = -1.0e+0_f8
#else
             DD_Hstar_old  = 1.0e+20_f8
             Henry_K0      = 1.0e+20_f8
             Henry_CR      = 1.34e+4_f8
#endif
             WD_RetFactor  = 1.0_fp

          CASE( 'I2O4' )
             FullName      = 'Diiodine tetraoxide'
             Formula       = 'I2O4'
             MW_g          = 318.0_fp
             Is_Gas        = T
             Is_Drydep     = T
             Is_Wetdep     = T
             Is_Photolysis = T
             DD_F0         = 0.0_fp
#if defined( NEW_HENRY_CONSTANTS )
             Henry_K0      = -1.0e+0_f8 * To_M_atm
             Henry_CR      = -1.0e+0_f8
#else
             DD_Hstar_old  = 1.0e+20_f8
             Henry_K0      = 1.0e+20_f8
             Henry_CR      = 1.34e+4_f8
#endif
             WD_RetFactor  = 1.0_fp

          CASE( 'BRSALA' )

             IF ( Present(oRadius) ) THEN
                IF ( .NOT. Present(Input_Opt) ) THEN
                   WRITE( 6, '(a)' ) REPEAT( '=', 79 )
                   WRITE( 6, * ) 'Error getting radius for species ',TRIM(Name)
                   WRITE( 6, * ) 'Input_Opt is missing!'
                   WRITE( 6, * ) 'In module Headers/species_database_mod.F90!'
                   RC = -1
                   RETURN
                ENDIF
                ! Mimic SALA
                Radius     = ( Input_Opt%SALA_REDGE_um(1) +                  &
                               Input_Opt%SALA_REDGE_um(2)  ) * 0.5e-6_fp
             ENDIF

             ! Halve the Kc (cloud condensate -> precip) rate
             ! for the temperature range 237 K <= T < 258 K.
             KcScale       = (/ 1.0_fp, 0.5_fp, 1.0_fp /)

             ! Turn off rainout only when 237 K <= T < 258K.
             RainEff       = (/ 1.0_fp, 0.0_fp, 1.0_fp /)

             FullName      = 'Fine sea salt bromine'
             Formula       = 'Br'
             MW_g          = 80.0_fp
             Is_Gas        = F
             Is_Drydep     = T
             Is_Wetdep     = T
             Density       = 2200.0_fp
             Radius        = Radius
             DD_AeroDryDep = T
             DD_F0         = 0.0_fp
             DD_Hstar_Old  = 0.0_fp
             WD_AerScavEff = 1.0_fp
             WD_KcScaleFac = KcScale
             WD_RainoutEff = RainEff

          CASE( 'BRSALC' )
             IF ( Present(oRadius) ) THEN
                IF ( .NOT. Present(Input_Opt) ) THEN
                   WRITE( 6, '(a)' ) REPEAT( '=', 79 )
                   WRITE( 6, * ) 'Error getting radius for species ',TRIM(Name)
                   WRITE( 6, * ) 'Input_Opt is missing!'
                   WRITE( 6, * ) 'In module Headers/species_database_mod.F90!'
                   RC = -1
                   RETURN
                ENDIF
                ! Mimic SALC
                Radius     = ( Input_Opt%SALC_REDGE_um(1) +                  &
                               Input_Opt%SALC_REDGE_um(2)  ) * 0.5e-6_fp
             ENDIF

             ! Halve the Kc (cloud condensate -> precip) rate
             ! for the temperature range 237 K <= T < 258 K.
             KcScale       = (/ 1.0_fp, 0.5_fp, 1.0_fp /)

             ! Turn off rainout only when 237 K <= T < 258K.
             RainEff       = (/ 1.0_fp, 0.0_fp, 1.0_fp /)

             FullName      = 'Coarse sea salt bromine'
             Formula       = 'Br'
             MW_g          = 80.0_fp
             Is_Gas        = F
             Is_Drydep     = T
             Is_Wetdep     = T
             Density       = 2200.0_fp
             Radius        = Radius
             DD_AeroDryDep = T
             DD_F0         = 0.0_fp
             DD_Hstar_Old  = 0.0_fp
             WD_AerScavEff = 1.0_fp
             WD_CoarseAer  = T
             WD_KcScaleFac = KcScale
             WD_RainoutEff = RainEff

          CASE( 'ISALA' )
             IF ( Present(oRadius) ) THEN
                IF ( .NOT. Present(Input_Opt) ) THEN
                   WRITE( 6, '(a)' ) REPEAT( '=', 79 )
                   WRITE( 6, * ) 'Error getting radius for species ',TRIM(Name)
                   WRITE( 6, * ) 'Input_Opt is missing!'
                   WRITE( 6, * ) 'In module Headers/species_database_mod.F90!'
                   RC = -1
                   RETURN
                ENDIF
                ! Mimic SALA
                Radius     = ( Input_Opt%SALA_REDGE_um(1) +                  &
                               Input_Opt%SALA_REDGE_um(2)  ) * 0.5e-6_fp
             ENDIF

             ! Halve the Kc (cloud condensate -> precip) rate
             ! for the temperature range 237 K <= T < 258 K.
             KcScale       = (/ 1.0_fp, 0.5_fp, 1.0_fp /)

             ! Turn off rainout only when 237 K <= T < 258K.
             RainEff       = (/ 1.0_fp, 0.0_fp, 1.0_fp /)

             FullName      = 'Fine sea salt iodine'
             Formula       = 'I'
             MW_g          = 127.0_fp
             Is_Gas        = F
             Is_Drydep     = T
             Is_Wetdep     = T
             Density       = 2200.0_fp
             Radius        = Radius
             DD_AeroDryDep = T
             DD_F0         = 0.0_fp
             DD_Hstar_Old  = 0.0_fp
             WD_AerScavEff = 1.0_fp
             WD_KcScaleFac = KcScale
             WD_RainoutEff = RainEff

          CASE( 'ISALC' )
             IF ( Present(oRadius) ) THEN
                IF ( .NOT. Present(Input_Opt) ) THEN
                   WRITE( 6, '(a)' ) REPEAT( '=', 79 )
                   WRITE( 6, * ) 'Error getting radius for species ',TRIM(Name)
                   WRITE( 6, * ) 'Input_Opt is missing!'
                   WRITE( 6, * ) 'In module Headers/species_database_mod.F90!'
                   RC = -1
                   RETURN
                ENDIF
                ! Mimic SALC
                Radius     = ( Input_Opt%SALC_REDGE_um(1) +                  &
                               Input_Opt%SALC_REDGE_um(2)  ) * 0.5e-6_fp
             ENDIF

             ! Halve the Kc (cloud condensate -> precip) rate
             ! for the temperature range 237 K <= T < 258 K.
             KcScale       = (/ 1.0_fp, 0.5_fp, 1.0_fp /)

             ! Turn off rainout only when 237 K <= T < 258K.
             RainEff       = (/ 1.0_fp, 0.0_fp, 1.0_fp /)

             FullName      = 'Coarse sea salt iodine'
             Formula       = 'I'
             MW_g          = 127.0_fp
             Is_Gas        = F
             Is_Drydep     = T
             Is_Wetdep     = T
             Density       = 2200.0_fp
             Radius        = Radius
             DD_AeroDryDep = T
             DD_F0         = 0.0_fp
             DD_Hstar_Old  = 0.0_fp
             WD_AerScavEff = 1.0_fp
             WD_CoarseAer  = T
             WD_KcScaleFac = KcScale
             WD_RainoutEff = RainEff

          CASE( 'AERI' )
             ! Mimic SO4 (AERI is essentiall iodine dissolved in aerosol)

             ! Halve the Kc (cloud condensate -> precip) rate
             ! for the temperature range 237 K <= T < 258 K.
             KcScale   = (/ 1.0_fp, 0.5_fp, 1.0_fp /)

             ! Turn off rainout only when 237 K <= T < 258K.
             RainEff   = (/ 1.0_fp, 0.0_fp, 1.0_fp /)

             ! (cf. Mian Chin's GOCART model)
             ! Minimum Vd over snow/ice : 0.01 cm/s
             ! Minimum Vd over land     : 0.01 cm/s
             DvzMinVal = (/ 0.01_fp, 0.01_fp /)

             FullName      = 'iodine on aerosol'
             Formula       = 'I'
             MW_g          = 127.0_fp
             Is_Gas        = F
             Is_Drydep     = T
             Is_Wetdep     = T
             DD_DvzAerSnow = 0.03_fp
             DD_DvzMinVal  = DvzMinVal
             DD_F0         = 0.0_fp
             DD_Hstar_Old  = 0.0_fp
             WD_AerScavEff = 1.0_fp
             WD_KcScaleFac = KcScale
             WD_RainoutEff = RainEff

          !==================================================================
          ! Species for the Rn-Pb-Be specialty simulation
          !==================================================================

          CASE( 'RN', '222RN', 'RN222' )
             FullName      = 'Radon-222 isotope'
             Formula       = 'Rn'
             MW_g          = 222.0_fp
             Is_Gas        = F
             Is_Drydep     = F
             Is_Wetdep     = F

          CASE( 'PB', '210PB', 'PB210' )

             ! Halve the Kc (cloud condensate -> precip) rate
             ! for the temperature range 237 K <= T < 258 K.
             KcScale       = (/ 1.0_fp, 0.5_fp, 1.0_fp /)

             ! Turn off rainout only when 237 K <= T < 258K.
             RainEff       = (/ 1.0_fp, 0.0_fp, 1.0_fp /)

             FullName      = 'Lead-210 isotope'
             Formula       = 'Pb'
             MW_g          = 210.0_fp
             Is_Gas        = F
             Is_Drydep     = T
             Is_Wetdep     = T
             DD_DvzAerSnow = 0.03_fp
             DD_F0         = 0.0_fp
             DD_HStar_old  = 0.0_fp
             WD_AerScavEff = 1.0_fp
             WD_KcScaleFac = KcScale
             WD_RainoutEff = RainEff

          CASE( 'BE', '7BE', 'BE7' )

             ! Halve the Kc (cloud condensate -> precip) rate
             ! for the temperature range 237 K <= T < 258 K.
             KcScale       = (/ 1.0_fp, 0.5_fp, 1.0_fp /)

             ! Turn off rainout only when 237 K <= T < 258K.
             RainEff       = (/ 1.0_fp, 0.0_fp, 1.0_fp /)

             FullName      = 'Beryllium-7 isotope'
             Formula       = 'Be'
             MW_g          = 7.0_fp
             Is_Gas        = F
             Is_Drydep     = T
             Is_Wetdep     = T
             DD_DvzAerSnow = 0.03_fp
             DD_F0         = 0.0_fp
             DD_HStar_old  = 0.0_fp
             WD_AerScavEff = 1.0_fp
             WD_KcScaleFac = KcScale
             WD_RainoutEff = RainEff

          !==================================================================
          ! Species for the Hg specialty simulation
          !==================================================================

          CASE( 'HG0',     'HG0_CAN', 'HG0_USA', 'HG0_CAM', 'HG0_SAM',      &
                'HG0_WAF', 'HG0_EAF', 'HG0_SAF', 'HG0_NAF', 'HG0_EUR',      &
                'HG0_EEU', 'HG0_MDE', 'HG0_SOV', 'HG0_SAS', 'HG0_EAS',      &
                'HG0_SEA', 'HG0_JPN', 'HG0_OCE', 'HG0_SO',  'HG0_BB',       &
                'HG0_GEO', 'HG0_ATL', 'HG0_NAT', 'HG0_SAT', 'HG0_NPA',      &
                'HG0_ARC', 'HG0_ANT', 'HG0_OCN', 'HG0_STR'   )

             ! Standardize tagged Hg0 species names
             SELECT CASE( TRIM( Name ) )
                CASE( 'HG0'     )
                   FullName = 'Elemental mercury'
                CASE( 'HG0_CAN' )
                   FullName = 'Elemental mercury from Canada'
                CASE( 'HG0_USA' )
                   FullName = 'Elemental mercury from USA'
                CASE( 'HG0_CAM' )
                   FullName = 'Elemental mercury from Central America'
                CASE( 'HG0_SAM' )
                   FullName = 'Elemental mercury from South America'
                CASE( 'HG0_WAF' )
                   FullName = 'Elemental mercury from West Africa'
                CASE( 'HG0_EAF' )
                   FullName = 'Elemental mercury from East Africa'
                CASE( 'HG0_SAF' )
                   FullName = 'Elemental mercury from South Africa'
                CASE( 'HG0_NAF' )
                   FullName = 'Elemental mercury from North Africa'
                CASE( 'HG0_EUR' )
                   FullName = 'Elemental mercury from OECD Europe'
                CASE( 'HG0_EEU' )
                   FullName = 'Elemental mercury from Eastern Europe'
                CASE( 'HG0_MDE' )
                   FullName = 'Elemental mercury from Middle East'
                CASE( 'HG0_SOV' )
                   FullName = 'Elemental mercury from former USSR'
                CASE( 'HG0_SAS' )
                   FullName = 'Elemental mercury from South Asia'
                CASE( 'HG0_EAS' )
                   FullName = 'Elemental mercury from East Asia'
                CASE( 'HG0_SEA' )
                   FullName = 'Elemental mercury from Southeast Asia'
                CASE( 'HG0_JPN' )
                   FullName = 'Elemental mercury from Japan'
                CASE( 'HG0_OCE' )
                   FullName = 'Elemental mercury from Oceania'
                CASE( 'HG0_SO'  )
                   FullName = 'Elemental mercury from Organic Soil'
                CASE( 'HG0_BB'  )
                   FullName = 'Elemental mercury from Biomass Burning'
                CASE( 'HG0_GEO' )
                   FullName = 'Elemental mercury from Geogenic Sources'
                CASE( 'HG0_ATL' )
                   FullName = 'Elemental mercury from Midatlantic Subsurface Water'
                CASE( 'HG0_NAT' )
                   FullName = 'Elemental mercury from N. Atlantic Subsurface Water'
                CASE( 'HG0_SAT' )
                   FullName = 'Elemental mercury from S. Atlantic Subsurface Water'
                CASE( 'HG0_NPA' )
                   FullName = 'Elemental mercury from N. Pacific Subsurface Water'
                CASE( 'HG0_ARC' )
                   FullName = 'Elemental mercury from Arctic Subsurface Water'
                CASE( 'HG0_ANT' )
                   FullName = 'Elemental mercury from Antarctic Subsurface Water'
                CASE( 'HG0_OCN' )
                   FullName = 'Elemental mercury from Indo-Pacific Subsurface Water'
                CASE( 'HG0_STR' )
                   FullName = 'Elemental mercury from Stratosphere'
             END SELECT

             FullName      = FullName
             Formula       = 'Hg'
             MW_g          = 201.0_fp
             Is_Gas        = T
             Is_Drydep     = T
             Is_Wetdep     = F
             Is_Hg0        = T
             DD_F0         = 1.0e-5_fp
             DD_Hstar_old  = 0.11_fp

          CASE( 'HG2',     'HG2_CAN', 'HG2_USA', 'HG2_CAM', 'HG2_SAM',      &
                'HG2_WAF', 'HG2_EAF', 'HG2_SAF', 'HG2_NAF', 'HG2_EUR',      &
                'HG2_EEU', 'HG2_MDE', 'HG2_SOV', 'HG2_SAS', 'HG2_EAS',      &
                'HG2_SEA', 'HG2_JPN', 'HG2_OCE', 'HG2_SO',  'HG2_BB',       &
                'HG2_GEO', 'HG2_ATL', 'HG2_NAT', 'HG2_SAT', 'HG2_NPA',      &
                'HG2_ARC', 'HG2_ANT', 'HG2_OCN', 'HG2_STR'   )

             ! Standardize tagged Hg0 species names
             SELECT CASE( TRIM( Name ) )
                CASE( 'HG2'     )
                   FullName = 'Divalent mercury'
                CASE( 'HG2_CAN' )
                   FullName = 'Divalent mercury from Canada'
                CASE( 'HG2_USA' )
                   FullName = 'Divalent mercury from USA'
                CASE( 'HG2_CAM' )
                   FullName = 'Divalent mercury from Central America'
                CASE( 'HG2_SAM' )
                   FullName = 'Divalent mercury from South America'
                CASE( 'HG2_WAF' )
                   FullName = 'Divalent mercury from West Africa'
                CASE( 'HG2_EAF' )
                   FullName = 'Divalent mercury from East Africa'
                CASE( 'HG2_SAF' )
                   FullName = 'Divalent mercury from South Africa'
                CASE( 'HG2_NAF' )
                   FullName = 'Divalent mercury from North Africa'
                CASE( 'HG2_EUR' )
                   FullName = 'Divalent mercury from OECD Europe'
                CASE( 'HG2_EEU' )
                   FullName = 'Divalent mercury from Eastern Europe'
                CASE( 'HG2_MDE' )
                   FullName = 'Divalent mercury from Middle East'
                CASE( 'HG2_SOV' )
                   FullName = 'Divalent mercury from former USSR'
                CASE( 'HG2_SAS' )
                   FullName = 'Divalent mercury from South Asia'
                CASE( 'HG2_EAS' )
                   FullName = 'Divalent mercury from East Asia'
                CASE( 'HG2_SEA' )
                   FullName = 'Divalent mercury from Southeast Asia'
                CASE( 'HG2_JPN' )
                   FullName = 'Divalent mercury from Japan'
                CASE( 'HG2_OCE' )
                   FullName = 'Divalent mercury from Oceania'
                CASE( 'HG2_SO'  )
                   FullName = 'Divalent mercury from Organic Soil'
                CASE( 'HG2_BB'  )
                   FullName = 'Divalent mercury from Biomass Burning'
                CASE( 'HG2_GEO' )
                   FullName = 'Divalent mercury from Geogenic Sources'
                CASE( 'HG2_ATL' )
                   FullName = 'Divalent mercury from Midatlantic Subsurface Water'
                CASE( 'HG2_NAT' )
                   FullName = 'Divalent mercury from N. Atlantic Subsurface Water'
                CASE( 'HG2_SAT' )
                   FullName = 'Divalent mercury from S. Atlantic Subsurface Water'
                CASE( 'HG2_NPA' )
                   FullName = 'Divalent mercury from N. Pacific Subsurface Water'
                CASE( 'HG2_ARC' )
                   FullName = 'Divalent mercury from Arctic Subsurface Water'
                CASE( 'HG2_ANT' )
                   FullName = 'Divalent mercury from Antarctic Subsurface Water'
                CASE( 'HG2_OCN' )
                   FullName = 'Divalent mercury from Indo-Pacific Subsurface Water'
                CASE( 'HG2_STR' )
                   FullName = 'Divalent mercury from Stratosphere'
             END SELECT

             FullName      = FullName
             Formula       = 'Hg'
             MW_g          = 201.0_fp
             Is_Gas        = T
             Is_Drydep     = T
             Is_Wetdep     = T
             Is_Hg2        = T
             DD_F0         = 0.0_fp
#if defined( NEW_HENRY_CONSTANTS )
             Henry_K0      = 1.40e+4_f8 * To_M_atm
             Henry_CR      = 5300.0_f8
#else
             DD_Hstar_old  = 1.00e+14_fp
             Henry_K0      = 1.40e+6_f8
             Henry_CR      = 8400.0_f8
#endif
             WD_RetFactor  = 1.0_fp

          CASE( 'HGP',     'HGP_CAN', 'HGP_USA', 'HGP_CAM', 'HGP_SAM',      &
                'HGP_WAF', 'HGP_EAF', 'HGP_SAF', 'HGP_NAF', 'HGP_EUR',      &
                'HGP_EEU', 'HGP_MDE', 'HGP_SOV', 'HGP_SAS', 'HGP_EAS',      &
                'HGP_SEA', 'HGP_JPN', 'HGP_OCE', 'HGP_SO',  'HGP_BB',       &
                'HGP_GEO', 'HGP_ATL', 'HGP_NAT', 'HGP_SAT', 'HGP_NPA',      &
                'HGP_ARC', 'HGP_ANT', 'HGP_OCN', 'HGP_STR' )

             ! Standardize tagged HgP species names
             SELECT CASE( TRIM( Name ) )
                 CASE( 'HGP'     )
                   FullName = 'Particulate mercury'
                CASE( 'HGP_CAN' )
                   FullName = 'Particulate mercury from Canada'
                CASE( 'HGP_USA' )
                   FullName = 'Particulate mercury from USA'
                CASE( 'HGP_CAM' )
                   FullName = 'Particulate mercury from Central America'
                CASE( 'HGP_SAM' )
                   FullName = 'Particulate mercury from South America'
                CASE( 'HGP_WAF' )
                   FullName = 'Particulate mercury from West Africa'
                CASE( 'HGP_EAF' )
                   FullName = 'Particulate mercury from East Africa'
                CASE( 'HGP_SAF' )
                   FullName = 'Particulate mercury from South Africa'
                CASE( 'HGP_NAF' )
                   FullName = 'Particulate mercury from North Africa'
                CASE( 'HGP_EUR' )
                   FullName = 'Particulate mercury from OECD Europe'
                CASE( 'HGP_EEU' )
                   FullName = 'Particulate mercury from Eastern Europe'
                CASE( 'HGP_MDE' )
                   FullName = 'Particulate mercury from Middle East'
                CASE( 'HGP_SOV' )
                   FullName = 'Particulate mercury from former USSR'
                CASE( 'HGP_SAS' )
                   FullName = 'Particulate mercury from South Asia'
                CASE( 'HGP_EAS' )
                   FullName = 'Particulate mercury from East Asia'
                CASE( 'HGP_SEA' )
                   FullName = 'Particulate mercury from Southeast Asia'
                CASE( 'HGP_JPN' )
                   FullName = 'Particulate mercury from Japan'
                CASE( 'HGP_OCE' )
                   FullName = 'Particulate mercury from Oceania'
                CASE( 'HGP_SO'  )
                   FullName = 'Particulate mercury from Organic Soil'
                CASE( 'HGP_BB'  )
                   FullName = 'Particulate mercury from Biomass Burning'
                CASE( 'HGP_GEO' )
                   FullName = 'Particulate mercury from Geogenic Sources'
                CASE( 'HGP_ATL' )
                   FullName = 'Particulate mercury from Midatlantic Subsurface Water'
                CASE( 'HGP_NAT' )
                   FullName = 'Particulate mercury from N. Atlantic Subsurface Water'
                CASE( 'HGP_SAT' )
                   FullName = 'Particulate mercury from S. Atlantic Subsurface Water'
                CASE( 'HGP_NPA' )
                   FullName = 'Particulate mercury from N. Pacific Subsurface Water'
                CASE( 'HGP_ARC' )
                   FullName = 'Particulate mercury from Arctic Subsurface Water'
                CASE( 'HGP_ANT' )
                   FullName = 'Particulate mercury from Antarctic Subsurface Water'
                CASE( 'HGP_OCN' )
                   FullName = 'Particulate mercury from Indo-Pacific Subsurface Water'
                CASE( 'HGP_STR' )
                   FullName = 'Particulate mercury from Stratosphere'
             END SELECT

             !%%% NOTE: In the prior code, the rainout fraction for HgP
             !%%% was computed before the shunt that turned off rainout
             !%%% when 237 K <= T < 258 K.  Therefore, in order to
             !%%% replicate the behavior of the prior code, we need to
             !%%% set the rainout efficiency (field WD_RainoutEff) to
             !%%% 1.0 for all temperature regimes.
             !%%%
             !%%% But in the prior code, the Kc rate (rate of conversion
             !%%% of cloud condensate to precipitation) for HgP was
             !%%% multiplied by 0.5 (as is done for most other aerosols)
             !%%% in routine F_AEROSOL.  This is part of the update to
             !%%% allow scavenging by snow that was implemented by Qiaoqiao
             !%%% Wang.
             !%%%
             !%%% Therefore, we have to ask Team Hg if we should allow
             !%%% the rainout for Hg to be turned off AND the Kc rate
             !%%% to be multiplied by 0.5.  They may have intended to
             !%%% not turn off rainout for HgP, but may also have been
             !%%% unaware of the scaling of the Kc rate by 0.5 in the
             !%%% F_AEROSOL routine.
             !%%%
             !%%% For the time being, we shall replicate the behavior of
             !%%% the prior code.  Therefore, we shall allow rainout of
             !%%% HgP to occur for 237 <= T < 258 K AND ALSO multiply
             !%%% Kc by 0.5.
             !%%%
             !%%% (bmy, 9/28/15)

             ! When 237 K <= T < 258 K:
             ! (1) Halve the Kc (cloud condensate -> precip) rate
             ! (2) DON'T TURN OFF RAINOUT! (at least until we talk to Team Hg)
             KcScale       = (/ 1.0_fp, 0.5_fp, 1.0_fp /)
             RainEff       = (/ 1.0_fp, 1.0_fp, 1.0_fp /)

             FullName      = FullName
             Formula       = 'Hg'
             MW_g          = 201.0_fp
             Is_Gas        = F
             Is_Drydep     = T
             Is_Wetdep     = T
             Is_HgP        = T
             DD_DvzAerSnow = 0.03_fp
             DD_F0         = 0.0_fp
             DD_Hstar_old  = 0.0_fp
             WD_AerScavEff = 1.0_fp
             WD_KcScaleFac = KcScale
             WD_RainoutEff = RainEff

          !==================================================================
          ! Species for the POPs specialty simulation
          !==================================================================

          CASE( 'POPG' )

             !----------------------------------------------------------------
             ! Notes for DD_Hstar_old from the v11-01c drydep_mod.F
             ! (cf. Carey Friedman and Helen Amos
             !----------------------------------------------------------------
             ! HSTAR is Henry's Law in mol/L/atm.
             ! For PHENANTHRENE, log Kaw = -2.76
             !  so unitless Kaw = 1.73*10^-3 and Kwa = 1/Kaw
             !  Divide by R (0.0821 atm/M/K) and T (298 K) and get
             !  HSTAR = 23.5 M/atm
             ! For PYRENE, log Kaw = -3.27
             !  Using the same conversion, HSTAR = 76.1 M/atm
             ! For BENZO[a]PYRENE, log Kaw = -4.51
             !  Using the same conversion, HSTAR = 1.32d3 M/atm
             !  All log Kaws from Ma et al., J Chem Eng Data 2010, 55:819
             !
             !----------------------------------------------------------------
             ! Notes for DD_KOA from the v11-01c drydep_mod.F
             ! (cf. Carey Friedman and Helen Amos)
             !----------------------------------------------------------------
             ! Adding Koa (octanol-ar partition coefficient) for POPs to
             !  account for accumulation in leaf cuticles
             !  Needs to be in units of mol/liter/atm as with HSTAR
             !  Divide unitless Koa at 298 K by product of R (0.0821 atm/M/K)
             !  and T (298 K)
             ! For PHENANTHRENE, log Koa = 7.64
             !  use same conversion as for HSTAR to get 1.78d6 M/atm
             ! For PYRENE, log Koa = 8.86
             !  use same conversion to get 2.96d7 M/atm
             ! For BENZO[a]PYRENE, log Koa = 11.48
             !  use same conversion to get 1.23d10 M/atm
             ! All log Koas from Ma et al., J Chem Eng Data 2010, 55:819
             ! Now add factor of 0.8 to account for 80% vol content of octanol
             !
             !----------------------------------------------------------------
             ! Notes for Henry_K0, Henry_CR from the v11-01c wetscav_mod.F
             ! (cf. Carey Friedman and Helen Amos
             !----------------------------------------------------------------
             ! Cocmpute liquid to gas ratio for POPs using
             !  the appropriate parameters for Henry's Law (M/atm, unitless
             !  Kaw divided by R (in atm/M/K, or 8.21d-2) and T (T = 298 K))
             !  as first argument and negative enthalpy of water-air exchange
             !  (kJ/mol) divided by R (in kJ/mol/K, or 8.32d-3) as second
             !  argument.
             ! For PHENANTHRENE, HSTAR = 2.35d1 and del_H = -5.65d3 (HSTAR from
             !  Ma et al, 2010 J. Chem. Eng. Data, and del_H from Scharzenbach
             !  2003, p200)
             ! For PYRENE, HSTAR = 7.61d1 and del_H = -5.17d3 (HSTAR from Ma
             !  et al and del_H from Scharzenbach 2003, p200)
             ! For BENZO[a]PYRENE, HSTAR = 1.32d3 and del_H = -5.65d3
             !  (HSTAR from Ma et al and Del_H the same as pyrene for now)
             !----------------------------------------------------------------
             IF ( .NOT. Present(Input_Opt) ) THEN
                WRITE( 6, '(a)' ) REPEAT( '=', 79 )
                WRITE( 6, * ) 'Error getting info for species ',TRIM(Name)
                WRITE( 6, * ) 'Input_Opt is missing!'
                WRITE( 6, * ) 'In module Headers/species_database_mod.F90!'
                RC = -1
                RETURN
             ENDIF

             SELECT CASE( TRIM( Input_Opt%POP_TYPE ) )
                CASE( 'PHE' )
                   FullName = 'Phenanthrene (gas phase)'
                   Formula  = 'C14H10'
                   MW_g     = 178.23_fp
                   KOA      = 4.37e+7_fp  * 0.0409_fp  * 0.8_fp
                   Hstar    = 1.0_fp      / 1.74e-3_fp * 0.0409_fp
                   K0       = 1.0_f8      / 1.74e-3_f8 / 8.21e-2_f8 / 298.0_f8
                   CR       = 47.0_f8     / 8.32e-3_f8
                CASE( 'PYR' )
                   FullName = 'Pyrene (gas phase)'
                   Formula  = 'C16H10'
                   MW_g     = 202.25_fp
                   KOA      = 7.24e+8_fp  * 0.0409_fp  * 0.8_fp
                   Hstar    = 1.0_fp      / 5.37e-4_fp * 0.0409_fp
                   K0       = 1.0_f8      / 5.37e-4_f8 / 8.21e-2_f8 / 298.0_f8
                   CR       = 43.0_f8     / 8.32e-3_f8
                CASE( 'BaP' )
                   FullName = 'Benzo(a)pyrene (gas phase)'
                   Formula  = 'C20H12'
                   MW_g     = 252.31_fp
                   KOA      = 3.02e+11_fp * 0.0409_fp  * 0.8_fp
                   Hstar    = 1.0_fp      / 3.10e-5_fp * 0.0409_fp
                   K0       = 1.0_f8      / 3.10e-5_f8 / 8.21e-2_f8 / 298.0_f8
                   CR       = 43.0_f8     / 8.32e-3_f8
             END SELECT

             FullName      = FullName
             Formula       = Formula
             MW_g          = MW_g
             Is_Gas        = T
             Is_Drydep     = T
             Is_Wetdep     = T
             DD_F0         = 0.0_fp
             DD_Hstar_old  = Hstar
             DD_KOA        = KOA
             Henry_K0      = K0
             Henry_CR      = CR
             WD_RetFactor  = 0.0_fp

          CASE( 'POPPBCPO', 'POPPOCPO' )
             IF ( .NOT. Present(Input_Opt) ) THEN
                WRITE( 6, '(a)' ) REPEAT( '=', 79 )
                WRITE( 6, * ) 'Error getting info for species ',TRIM(Name)
                WRITE( 6, * ) 'Input_Opt is missing!'
                WRITE( 6, * ) 'In module Headers/species_database_mod.F90!'
                RC = -1
                RETURN
             ENDIF


             ! These have identical properties except for the mol weights ...
             SELECT CASE( TRIM( Input_Opt%POP_TYPE ) )
                CASE( 'PHE' )
                   MW_g     = 178.23_fp
                   FullName = 'Phenanthrene particles on'
                   Formula  = 'C14H10'
                CASE( 'PYR' )
                   MW_g     = 202.25_fp
                   FullName = 'Pyrene particles on'
                   Formula  = 'C16H10'
                CASE( 'BaP' )
                   MW_g     = 252.31_fp
                   FullName = 'Benzo(a)pyrene particles on'
                   Formula  = 'C20H12'
             END SELECT

             ! ... and the names and rainout efficiencies.
             SELECT CASE( TRIM( Name ) )

                CASE( 'POPPBCPO' )
                   FullName = TRIM( FullName ) // ' hydrophobic black carbon'
                   Formula  = ''

                   ! Halve the Kc (cloud condensate -> precip) rate
                   ! for the temperature range 237 K <= T < 258 K.
                   KcScale  = (/ 1.0_fp, 1.0_fp, 0.5_fp /)

                   ! Allow rainout of POPPBCPO when T < 258 K, because
                   ! POPPBCPO is considered to be IN.
                   RainEff  = (/ 1.0_fp, 1.0_fp, 0.0_fp /)

                CASE( 'POPPOCPO' )
                   FullName = TRIM( FullName ) // ' hydrophobic organic carbon'
                   Formula  = ''

                   ! For all temperatures:
                   ! (1) Halve the Kc (cloud condensate -> precip) rate
                   ! (2) Turn off rainout (it's hydrophobic)
                   KcScale  = (/ 0.5_fp, 0.5_fp, 0.5_fp /)
                   RainEff  = (/ 0.0_fp, 0.0_fp, 0.0_fp /)

             END SELECT

             FullName      = FullName
             Formula       = Formula
             MW_g          = MW_g
             Is_Gas        = F
             Is_Drydep     = T
             Is_Wetdep     = T
             DD_DvzAerSnow = 0.03_fp
             DD_F0         = 0.0_fp
             DD_Hstar_old  = 0.0_fp
             WD_AerScavEff = 1.0_fp
             WD_KcScaleFac = KcScale
             WD_RainoutEff = RainEff

          CASE( 'POPPBCPI', 'POPPOCPI' )

             IF ( .NOT. Present(Input_Opt) ) THEN
                WRITE( 6, '(a)' ) REPEAT( '=', 79 )
                WRITE( 6, * ) 'Error getting info for species ',TRIM(Name)
                WRITE( 6, * ) 'Input_Opt is missing!'
                WRITE( 6, * ) 'In module Headers/species_database_mod.F90!'
                RC = -1
                RETURN
             ENDIF

             ! These have identical properties except for the mol weights ...
             SELECT CASE( TRIM( Input_Opt%POP_TYPE ) )
                CASE( 'PHE' )
                   MW_g     = 178.23_fp
                   FullName = 'Phenanthrene particles on'
                   Formula  = 'C14H10'
                CASE( 'PYR' )
                   MW_g     = 202.25_fp
                   FullName = 'Pyrene particles on'
                   Formula  = 'C16H10'
                CASE( 'BaP' )
                   MW_g     = 252.31_fp
                   FullName = 'Benzo(a)pyrene particles on'
                   Formula  = 'C20H12'
             END SELECT

             ! ... and the names
             SELECT CASE( TRIM( Name ) )
                CASE( 'POPPBCPI' )
                   FullName = TRIM( FullName ) // ' hydrophilic black carbon'
                   Formula  = ''
                CASE( 'POPPOCPI' )
                   FullName = TRIM( FullName ) // ' hydrophilic organic carbon'
                   Formula  = ''
             END SELECT

             ! Halve the Kc (cloud condensate -> precip) rate
             ! for the temperature range 237 K <= T < 258 K.
             KcScale = (/ 1.0_fp, 0.5_fp, 1.0_fp /)

             ! Turn off rainout only when 237 K <= T < 258K.
             RainEff = (/ 1.0_fp, 0.0_fp, 1.0_fp /)

             FullName      = FullName
             Formula       = Formula
             MW_g          = MW_g
             Is_Gas        = F
             Is_Drydep     = T
             Is_Wetdep     = T
             DD_DvzAerSnow = 0.03_fp
             DD_F0         = 0.0_fp
             DD_Hstar_old  = 0.0_fp
             WD_AerScavEff = 1.0_fp
             WD_KcScaleFac = KcScale
             WD_RainoutEff = RainEff

          !==================================================================
          ! Species for the CO2 specialty simulation
          !==================================================================

          CASE( 'CO2',    'CO2FF', 'CO2OC', 'CO2BAL', 'CO2BB', 'CO2BF',     &
                'CO2NTE', 'CO2SE', 'CO2AV', 'CO2CH',  'CO2CORR'         )

             ! These all have identical properties except for the names
             ! Add TOMAS bin number to full name

             SELECT CASE( TRIM( Name ) )
                CASE( 'CO2FF' )
                   FullName = 'Carbon dioxide from fossil fuel emissions'
                CASE( 'CO2OC' )
                   FullName = 'Carbon dioxide from ocean emissions'
                CASE( 'CO2BAL' )
                   FullName = 'Carbon dioxide from balanced biosphere'
                CASE( 'CO2BB' )
                   FullName = 'Carbon dioxide  from biomass burning emissions'
                CASE( 'CO2BF' )
                   FullName = 'Carbon dioxide  from biofuel emissions'
                CASE( 'CO2NTE' )
                   FullName = 'Carbon dioxide from net terrestrial exchange'
                CASE( 'CO2SE' )
                   FullName = 'Carbon dioxide from ship emissions'
                CASE( 'CO2AV' )
                   FullName = 'Carbon dioxide from aviation emissions'
                CASE( 'CO2CH' )
                   FullName = 'Carbon dioxide from chemical sources'
                CASE( 'CO2CORR' )
                   FullName = 'Carbon dioxide chemical source surface correction'
                CASE DEFAULT
                   FullName = 'Carbon dioxide'
             END SELECT

             ! Set special default background for CO2
             SELECT CASE( TRIM( Name ) )
                CASE( 'CO2' )
                   BackgroundVV = 3.55e-04_fp
                CASE DEFAULT
                   BackgroundVV = MISSING_VV
             END SELECT

             FullName      = FullName
             Formula       = 'CO2'
             MW_g          = 44.0_fp
             BackgroundVV  = BackgroundVV
             Is_Gas        = T
             Is_Drydep     = F
             Is_Wetdep     = F

#if defined( TOMAS )
          !==================================================================
          ! Species for the TOMAS microphysics simulations
          !==================================================================

          CASE( 'AW1',  'AW2',  'AW3',  'AW4',  'AW5',  'AW6',  'AW7',      &
                'AW8',  'AW9',  'AW10', 'AW11', 'AW12', 'AW13', 'AW14',     &
                'AW15', 'AW16', 'AW17', 'AW18', 'AW19', 'AW20', 'AW21',     &
                'AW22', 'AW23', 'AW24', 'AW25', 'AW26', 'AW27', 'AW28',     &
                'AW29', 'AW30', 'AW31', 'AW32', 'AW33', 'AW34', 'AW35',     &
                'AW36', 'AW37', 'AW38', 'AW39', 'AW40'  )

             ! Add TOMAS bin number to full name
             FullName      = 'Aerosol water, size bin ='
             C             = LEN_TRIM( Name )
             IF ( C == 3 ) THEN
                FullName   = TRIM( FullName ) // ' ' // Name(C:C)
             ELSE
                FullName   = TRIM( FullName ) // ' ' // Name(C-1:C)
             ENDIF

             ! Halve the Kc (cloud condensate -> precip) rate
             ! for the temperature range 237 K <= T < 258 K.
             KcScale       = (/ 1.0_fp, 0.5_fp, 1.0_fp /)

             ! Turn off rainout only when 237 K <= T < 258K.
             RainEff       = (/ 1.0_fp, 0.0_fp, 1.0_fp /)

             Formula       = ''
             MW_g          = 18.0_fp
             Is_Gas        = F
             Is_Drydep     = F
             Is_Wetdep     = F
             DD_DvzAerSnow = 0.03_fp
             DD_F0         = 0.0_fp
             DD_Hstar_old  = 0.0_fp
             MP_SizeResAer = T
             WD_AerScavEff = 1.0_fp
             WD_KcScaleFac = KcScale
             WD_RainoutEff = RainEff

          CASE( 'DUST1',  'DUST2',  'DUST3',  'DUST4',  'DUST5', &
                'DUST6',  'DUST7',  'DUST8',  'DUST9',  'DUST10',&
                'DUST11', 'DUST12', 'DUST13', 'DUST14', 'DUST15',&
                'DUST16', 'DUST17', 'DUST18', 'DUST19', 'DUST20',&
                'DUST21', 'DUST22', 'DUST23', 'DUST24', 'DUST25',&
                'DUST26', 'DUST27', 'DUST28', 'DUST29', 'DUST30',&
                'DUST31', 'DUST32', 'DUST33', 'DUST34', 'DUST35',&
                'DUST36', 'DUST37', 'DUST38', 'DUST39', 'DUST40'  )

             ! Add TOMAS bin number to full name
             FullName      = 'Mineral dust, size bin ='
             C             = LEN_TRIM( Name )
             IF ( C == 5 ) THEN
                FullName   = TRIM( FullName ) // ' ' // Name(C:C)
             ELSE
                FullName   = TRIM( FullName ) // ' ' // Name(C-1:C)
             ENDIF

             ! Halve the Kc (cloud condensate -> precip) rate
             ! for the temperature range 237 K <= T < 258 K.
             KcScale       = (/ 1.0_fp, 0.5_fp, 1.0_fp /)

             ! Turn off rainout only when 237 K <= T < 258K.
             RainEff       = (/ 1.0_fp, 0.0_fp, 1.0_fp /)

             FullName      = FullName
             Formula       = ''
             MW_g          = 100.0_fp
             Is_Gas        = F
             Is_Drydep     = F
             Is_Wetdep     = T
             DD_DvzAerSnow = 0.03_fp
             DD_F0         = 0.0_fp
             DD_Hstar_old  = 0.0_fp
             MP_SizeResAer = T
             WD_AerScavEff = 1.0_fp
             WD_KcScaleFac = KcScale
             WD_RainoutEff = RainEff

          CASE( 'ECIL1',  'ECIL2',  'ECIL3',  'ECIL4',  'ECIL5'
                'ECIL6',  'ECIL7',  'ECIL8',  'ECIL9',  'ECIL10'
                'ECIL11', 'ECIL12', 'ECIL13', 'ECIL14', 'ECIL15'
                'ECIL16', 'ECIL17', 'ECIL18', 'ECIL19', 'ECIL20'
                'ECIL21', 'ECIL22', 'ECIL23', 'ECIL24', 'ECIL25'
                'ECIL26', 'ECIL27', 'ECIL28', 'ECIL29', 'ECIL30'
                'ECIL31', 'ECIL32', 'ECIL33', 'ECIL34', 'ECIL35'
                'ECIL36', 'ECIL37', 'ECIL38', 'ECIL39', 'ECIL40'  )

             ! Add TOMAS bin number to full name
             FullName      = 'Hydrophilic elemental carbon, size bin ='
             C             = LEN_TRIM( Name )
             IF ( C == 5 ) THEN
                FullName   = TRIM( FullName ) // ' ' // Name(C:C)
             ELSE
                FullName   = TRIM( FullName ) // ' ' // Name(C-1:C)
             ENDIF

             ! Halve the Kc (cloud condensate -> precip) rate
             ! for the temperature range 237 K <= T < 258 K.
             KcScale       = (/ 1.0_fp, 0.5_fp, 1.0_fp /)

             ! Turn off rainout only when 237 K <= T < 258K.
             RainEff       = (/ 1.0_fp, 0.0_fp, 1.0_fp /)

             Formula       = ''
             MW_g          = 12.0_fp
             Is_Gas        = F
             Is_Drydep     = F
             Is_Wetdep     = T
             DD_DvzAerSnow = 0.03_fp
             DD_F0         = 0.0_fp
             DD_Hstar_old  = 0.0_fp
             MP_SizeResAer = T
             WD_AerScavEff = 1.0_fp
             WD_KcScaleFac = KcScale
             WD_RainoutEff = RainEff

          CASE( 'ECOB1',  'ECOB2',  'ECOB3',  'ECOB4',  'ECOB5'
                'ECOB6',  'ECOB7',  'ECOB8',  'ECOB9',  'ECOB10'
                'ECOB11', 'ECOB12', 'ECOB13', 'ECOB14', 'ECOB15'
                'ECOB16', 'ECOB17', 'ECOB18', 'ECOB19', 'ECOB20'
                'ECOB21', 'ECOB22', 'ECOB23', 'ECOB24', 'ECOB25'
                'ECOB26', 'ECOB27', 'ECOB28', 'ECOB29', 'ECOB30'
                'ECOB31', 'ECOB32', 'ECOB33', 'ECOB34', 'ECOB35'
                'ECOB36', 'ECOB37', 'ECOB38', 'ECOB39', 'ECOB40'  )

             ! Add TOMAS bin number to full name
             FullName      = 'Hydrophobic elemental carbon, size bin ='
             C             = LEN_TRIM( Name )
             IF ( C == 5 ) THEN
                FullName   = TRIM( FullName ) // ' ' // Name(C:C)
             ELSE
                FullName   = TRIM( FullName ) // ' ' // Name(C-1:C)
             ENDIF

             ! Halve the Kc (cloud condensate -> precip) rate
             ! for the temperature range 237 K <= T < 258 K.
             KcScale       = (/ 1.0_fp, 0.5_fp, 1.0_fp /)

             ! Turn off rainout only when 237 K <= T < 258K.
             RainEff       = (/ 1.0_fp, 0.0_fp, 1.0_fp /)

             Formula       = ''
             MW_g          = 12.0_fp
             Is_Gas        = F
             Is_Drydep     = F
             Is_Wetdep     = F
             DD_DvzAerSnow = 0.03_fp
             DD_F0         = 0.0_fp
             DD_Hstar_old  = 0.0_fp
             MP_SizeResAer = T
             WD_AerScavEff = 1.0_fp
             WD_KcScaleFac = KcScale
             WD_RainoutEff = RainEff

          CASE( 'H2SO4' )

             !%%% NOTE: The TOMAS H2SO4 species dry-deposits like a gas,
             !%%% wet-deposits as an aerosol.  So we need to give this
             !%%% both gas and aerosol properties (ewl, bmy, 10/13/15)

             ! Halve the Kc (cloud condensate -> precip) rate
             ! for the temperature range 237 K <= T < 258 K.
             KcScale       = (/ 1.0_fp, 0.5_fp, 1.0_fp /)

             ! Turn off rainout only when 237 K <= T < 258K.
             RainEff       = (/ 1.0_fp, 0.0_fp, 1.0_fp /)

             FullName      = 'Sulfuric acid'
             Formula       = 'H2SO4'
             MW_g          = 98.0_fp
             MolecRatio    = 1.0_fp
             Is_Gas        = T
             Is_Drydep     = T
             Is_Wetdep     = T
             DD_F0         = 0.0_fp
             DD_Hstar_old  = 1.00e+5_fp
             MP_SizeResAer = T
             WD_AerScavEff = 1.0_fp
             WD_KcScaleFac = KcScale
             WD_RainoutEff = RainEff

          CASE( 'NK1',  'NK2',  'NK3',  'NK4',  'NK5',  'NK6',  'NK7',      &
                'NK8',  'NK9',  'NK10', 'NK11', 'NK12', 'NK13', 'NK14',     &
                'NK15', 'NK16', 'NK17', 'NK18', 'NK19', 'NK20', 'NK21',     &
                'NK22', 'NK23', 'NK24', 'NK25', 'NK26', 'NK27', 'NK28',     &
                'NK29', 'NK30', 'NK31', 'NK32', 'NK33', 'NK34', 'NK35',     &
                'NK36', 'NK37', 'NK38', 'NK39', 'NK40'                  )

             ! Add TOMAS bin number to full name
             FullName      = 'Aerosol number, size bin ='
             C             = LEN_TRIM( Name )
             IF ( C == 3 ) THEN
                FullName   = TRIM( FullName ) // ' ' // Name(C:C)
             ELSE
                FullName   = TRIM( FullName ) // ' ' // Name(C-1:C)
             ENDIF

             ! Halve the Kc (cloud condensate -> precip) rate
             ! for the temperature range 237 K <= T < 258 K.
             KcScale       = (/ 1.0_fp, 0.5_fp, 1.0_fp /)

             ! Turn off rainout only when 237 K <= T < 258K.
             RainEff       = (/ 1.0_fp, 0.0_fp, 1.0_fp /)

             Formula       = ''
             MW_g          = 1.0_fp
             Is_Gas        = F
             Is_Drydep     = T
             Is_Wetdep     = T
             DD_DvzAerSnow = 0.03_fp
             DD_DustDryDep = T
             DD_F0         = 0.0_fp
             DD_Hstar_old  = 0.0_fp
             MP_SizeResNum = T
             WD_AerScavEff = 1.0_fp
             WD_KcScaleFac = KcScale
             WD_RainoutEff = RainEff

          CASE( 'OCIL1',  'OCIL2',  'OCIL3',  'OCIL4',  'OCIL5',&
                'OCIL6',  'OCIL7',  'OCIL8',  'OCIL9',  'OCIL10',&
                'OCIL11', 'OCIL12', 'OCIL13', 'OCIL14', 'OCIL15',&
                'OCIL16', 'OCIL17', 'OCIL18', 'OCIL19', 'OCIL20',&
                'OCIL21', 'OCIL22', 'OCIL23', 'OCIL24', 'OCIL25',&
                'OCIL26', 'OCIL27', 'OCIL28', 'OCIL29', 'OCIL30',&
                'OCIL31', 'OCIL32', 'OCIL33', 'OCIL34', 'OCIL35',&
                'OCIL36', 'OCIL37', 'OCIL38', 'OCIL39', 'OCIL40'  )

             ! Add TOMAS bin number to full name
             FullName = 'Hydrophilic organic carbon, size bin ='
             C             = LEN_TRIM( Name )
             IF ( C == 5 ) THEN
                FullName = TRIM( FullName ) // ' ' // Name(C:C)
             ELSE
                FullName = TRIM( FullName ) // ' ' // Name(C-1:C)
             ENDIF

             ! Halve the Kc (cloud condensate -> precip) rate
             ! for the temperature range 237 K <= T < 258 K.
             KcScale       = (/ 1.0_fp, 0.5_fp, 1.0_fp /)

             ! Turn off rainout only when 237 K <= T < 258K.
             RainEff       = (/ 1.0_fp, 0.0_fp, 1.0_fp /)

             Formula       = ''
             MW_g          = 12.0_fp
             Is_Gas        = F
             Is_Drydep     = F
             Is_Wetdep     = T
             DD_DvzAerSnow = 0.03_fp
             DD_F0         = 0.0_fp
             DD_Hstar_old  = 0.0_fp
             MP_SizeResAer = T
             WD_AerScavEff = 1.0_fp
             WD_KcScaleFac = KcScale
             WD_RainoutEff = RainEff

          CASE( 'OCOB1',  'OCOB2',  'OCOB3',  'OCOB4',  'OCOB5',&
                'OCOB6',  'OCOB7',  'OCOB8',  'OCOB9',  'OCOB10',&
                'OCOB11', 'OCOB12', 'OCOB13', 'OCOB14', 'OCOB15',&
                'OCOB16', 'OCOB17', 'OCOB18', 'OCOB19', 'OCOB20',&
                'OCOB21', 'OCOB22', 'OCOB23', 'OCOB24', 'OCOB25',&
                'OCOB26', 'OCOB27', 'OCOB28', 'OCOB29', 'OCOB30',&
                'OCOB31', 'OCOB32', 'OCOB33', 'OCOB34', 'OCOB35',&
                'OCOB36', 'OCOB37', 'OCOB38', 'OCOB39', 'OCOB40'  )

             ! Add TOMAS bin number to full name
             FullName      = 'Hydrophobic organic carbon, size bin ='
             C             = LEN_TRIM( Name )
             IF ( C == 5 ) THEN
                FullName = TRIM( FullName ) // ' ' // Name(C:C)
             ELSE
                FullName = TRIM( FullName ) // ' ' // Name(C-1:C)
             ENDIF

             ! Halve the Kc (cloud condensate -> precip) rate
             ! for the temperature range 237 K <= T < 258 K.
             KcScale       = (/ 1.0_fp, 0.5_fp, 1.0_fp /)

             ! Turn off rainout only when 237 K <= T < 258K.
             RainEff       = (/ 1.0_fp, 0.0_fp, 1.0_fp /)

             FullName      = FullName
             Formula       = ''
             MW_g          = 12.0_fp
             Is_Gas        = F
             Is_Drydep     = F
             Is_Wetdep     = T
             DD_DvzAerSnow = 0.03_fp
             DD_F0         = 0.0_fp
             DD_Hstar_old  = 0.0_fp
             MP_SizeResAer = T
             WD_AerScavEff = 1.0_fp
             WD_KcScaleFac = KcScale
             WD_RainoutEff = RainEff

          CASE( 'SF1',  'SF2',  'SF3',  'SF4',  'SF5',  'SF6',  'SF7',      &
                'SF8',  'SF9',  'SF10', 'SF11', 'SF12', 'SF13', 'SF14',     &
                'SF15', 'SF16', 'SF17', 'SF18', 'SF19', 'SF20', 'SF21',     &
                'SF22', 'SF23', 'SF24', 'SF25', 'SF26', 'SF27', 'SF28',     &
                'SF29', 'SF30', 'SF31', 'SF32', 'SF33', 'SF34', 'SF35',     &
                'SF36', 'SF37', 'SF38', 'SF39', 'SF40'                  )

             ! Add TOMAS bin number to full name
             FullName      = 'Sulfate aerosol, size bin ='
             C             = LEN_TRIM( Name )
             IF ( C == 3 ) THEN
                FullName   = TRIM( FullName ) // ' ' // Name(C:C)
             ELSE
                FullName   = TRIM( FullName ) // ' ' // Name(C-1:C)
             ENDIF

             ! Halve the Kc (cloud condensate -> precip) rate
             ! for the temperature range 237 K <= T < 258 K.
             KcScale       = (/ 1.0_fp, 0.5_fp, 1.0_fp /)

             ! Turn off rainout only when 237 K <= T < 258K.
             RainEff       = (/ 1.0_fp, 0.0_fp, 1.0_fp /)

             Formula       = ''
             MW_g          = 96.0_fp
             Is_Gas        = F
             Is_Drydep     = F
             Is_Wetdep     = T
             DD_DvzAerSnow = 0.03_fp
             DD_F0         = 0.0_fp
             DD_Hstar_old  = 0.0_fp
             MP_SizeResAer = T
             WD_AerScavEff = 1.0_fp
             WD_KcScaleFac = KcScale
             WD_RainoutEff = RainEff

          CASE( 'SS1',  'SS2',  'SS3',  'SS4',  'SS5',  'SS6',  'SS7',      &
                'SS8',  'SS9',  'SS10', 'SS11', 'SS12', 'SS13', 'SS14',     &
                'SS15', 'SS16', 'SS17', 'SS18', 'SS19', 'SS20', 'SS21',     &
                'SS22', 'SS23', 'SS24', 'SS25', 'SS26', 'SS27', 'SS28',     &
                'SS29', 'SS30', 'SS31', 'SS32', 'SS33', 'SS34', 'SS35',     &
                'SS36', 'SS37', 'SS38', 'SS39', 'SS40'                  )

             ! Add TOMAS bin number to full name
             FullName      = 'Sea salt aerosol, size bin = '
             C             = LEN_TRIM( Name )
             IF ( C == 3 ) THEN
                FullName = TRIM( FullName ) // ' ' // Name(C:C)
             ELSE
                FullName = TRIM( FullName ) // ' ' // Name(C-1:C)
             ENDIF

             ! Halve the Kc (cloud condensate -> precip) rate
             ! for the temperature range 237 K <= T < 258 K.
             KcScale       = (/ 1.0_fp, 0.5_fp, 1.0_fp /)

             ! Turn off rainout only when 237 K <= T < 258K.
             RainEff       = (/ 1.0_fp, 0.0_fp, 1.0_fp /)

             FullName      = FullName
             Formula       = ''
             MW_g          = 58.5_fp
             Is_Gas        = F
             Is_Drydep     = F
             Is_Wetdep     = T
             DD_DvzAerSnow = 0.03_fp
             DD_F0         = 0.0_fp
             DD_Hstar_old  = 0.0_fp
             MP_SizeResAer = T
             WD_AerScavEff = 1.0_fp
             WD_KcScaleFac = KcScale
             WD_RainoutEff = RainEff

#endif

          !==================================================================
          ! Additional species for FAST-JX photolysis
          !==================================================================

          CASE( 'O2' )
             BackgroundVV  = 2.095e-01_fp
             FullName      = 'Molecular oxygen'
             Formula       = 'O2'
             Is_Gas        = T
             Is_Drydep     = F
             Is_Wetdep     = F
             Is_Photolysis = T

          CASE( 'INPN' )
             FullName      = 'Peroxide from INO2'
             Formula       = 'O2NOCH2C(OOH)(CH3)CH=CH2'
             Is_Gas        = T
             Is_Drydep     = F
             Is_Wetdep     = F
             Is_Photolysis = T

          CASE( 'PRPN' )
             FullName      = 'Peroxide from PRN1'
             Formula       = 'O2NOCH2CH(OOH)CH3'
             Is_Gas        = T
             Is_Drydep     = F
             Is_Wetdep     = F
             Is_Photolysis = T

          CASE( 'ETP' )
             FullName      = 'Ethylhydroperoxide'
             Formula       = 'CH3CH2OOH '
             Is_Gas        = T
             Is_Drydep     = F
             Is_Wetdep     = F
             Is_Photolysis = T

          CASE( 'RA3P' )
             FullName      = 'Peroxide from A3O2'
             Formula       = 'CH3CH2CH2OOH'
             Is_Gas        = T
             Is_Drydep     = F
             Is_Wetdep     = F
             Is_Photolysis = T

          CASE( 'RB3P' )
             FullName      = 'Peroxide from B3O2'
             Formula       = 'CH3CH(OOH)CH3'
             Is_Gas        = T
             Is_Drydep     = F
             Is_Wetdep     = F
             Is_Photolysis = T

          CASE( 'R4P' )
             FullName      = 'Peroxide from R4O2'
             Formula       = 'CH3CH2CH2CH2OOH'
             Is_Gas        = T
             Is_Drydep     = F
             Is_Wetdep     = F
             Is_Photolysis = T

          CASE( 'PP' )
             FullName      = 'Peroxide from PO2'
             Formula       = 'HOCH2CH(OOH)CH3'
             Is_Gas        = T
             Is_Drydep     = F
             Is_Wetdep     = F
             Is_Photolysis = T

          CASE( 'RP' )
             FullName      = 'Peroxide from RCO3'
             Formula       = 'CH3CH2C(O)OOH'
             Is_Gas        = T
             Is_Drydep     = F
             Is_Wetdep     = F
             Is_Photolysis = T

          CASE( 'IAP' )
             FullName      = 'Peroxide from IAO2'
             Formula       = 'HOCH2C(CH3)(OOH)CH(OH)CHO'
             Is_Gas        = T
             Is_Drydep     = F
             Is_Wetdep     = F
             Is_Photolysis = T

          CASE( 'ISNP' )
             FullName      = 'Isoprene nitrate'
             Formula       = 'HOCH2C(OOH)(CH3)CH(ONO2)CH2OH'
             Is_Gas        = T
             Is_Drydep     = F
             Is_Wetdep     = F
             Is_Photolysis = T

          CASE( 'VRP' )
             FullName      = 'Peroxide from VRO2'
             Formula       = 'HOCH2CH(OOH)C(O)CH3'
             Is_Gas        = T
             Is_Drydep     = F
             Is_Wetdep     = F
             Is_Photolysis = T

          CASE( 'MRP' )
             FullName      = 'Peroxide from MRO2'
             Formula       = 'HOCH2C(OOH)(CH3)CHO '
             Is_Gas        = T
             Is_Drydep     = F
             Is_Wetdep     = F
             Is_Photolysis = T
             
          CASE( 'MAOP' )
             FullName      = 'Peroxide from MAO3'
             Formula       = 'CH2=C(CH3)C(O)OOH'
             Is_Gas        = T
             Is_Drydep     = F
             Is_Wetdep     = F
             Is_Photolysis = T

          CASE( 'ATOOH' )
             FullName      = 'ATO2 peroxide'
             Formula       = 'CH3C(O)CH2OOH '
             Is_Gas        = T
             Is_Drydep     = F
             Is_Wetdep     = F
             Is_Photolysis = T

          CASE( 'PIP' )
             FullName      = 'Peroxide from MTPA'
             Formula       = ''
             Is_Gas        = T
             Is_Drydep     = F
             Is_Wetdep     = F
             Is_Photolysis = T

          CASE( 'HO2' )
             FullName      = 'Hydroperoxyl radical'
             Formula       = 'HO2'
             BackgroundVV  = 4.0e-15_fp
             Is_Gas        = T
             Is_Drydep     = F
             Is_Wetdep     = F

          CASE( 'MO2' )
             FullName      = 'Methylperoxy radical'
             Formula       = 'CH3O2'
             BackgroundVV  = 4.0e-15_fp
             Is_Gas        = T
             Is_Drydep     = F
             Is_Wetdep     = F

          CASE( 'OH' )
             FullName      = 'Hydroxyl radical'
             Formula       = 'OH'
             BackgroundVV  = 4.0e-15_fp
             Is_Gas        = T
             Is_Drydep     = F
             Is_Wetdep     = F

          CASE( 'H2' )
             FullName      = 'Molecular hydrogen'
             Formula       = 'H2'
             BackgroundVV  = 5.0e-07_fp
             Is_Gas        = T
             Is_Drydep     = F
             Is_Wetdep     = F

          CASE( 'N' )
             FullName      = 'Atomic nitrogen'
             Formula       = 'N'
             BackgroundVV  = 4.0e-20_fp
             Is_Gas        = T
             Is_Drydep     = F
             Is_Wetdep     = F

          CASE( 'N2' )
             FullName      = 'Molecular nitrogen'
             Formula       = 'N2'
             BackgroundVV  = 7.808e-1_fp
             Is_Gas        = T
             Is_Drydep     = F
             Is_Wetdep     = F

          CASE( 'O1D' )
             FullName      = 'Excited atomic oxygen (1D)'
             Formula       = 'O(1D)'
             BackgroundVV  = 1.0e-15_fp
             Is_Gas        = T
             Is_Drydep     = F
             Is_Wetdep     = F

          !==================================================================
          ! Special handling for species not found in the list above
          !==================================================================
          CASE DEFAULT

             ! Check if passive species
             IsPassive = .FALSE.
             !MW_g = 0.0_fp
             !BackgroundVV = 0.0_fp

             IF ( Present(Input_Opt) ) THEN
                IF ( Input_Opt%NPASSIVE > 0 ) THEN
   
                   ! Loop over all passive species
                   DO P = 1, Input_Opt%NPASSIVE
                      IF ( TRIM(Name) ==    &
                           TRIM(Input_Opt%PASSIVE_NAME(P)) ) THEN
                         IsPassive = .TRUE.
                         BackgroundVV = Input_Opt%PASSIVE_INITCONC(P)
                         MW_g   = Input_Opt%PASSIVE_MW(P)
                         EXIT
                      ENDIF
                   ENDDO
                ENDIF
             ENDIF

             !---------------------------------------------------------------
             ! Add passive species if it is listed in the optional passive
             ! species menu in input.geos. When reading the input file, all
             ! listed passive species quantities are written to local
             ! variables within passive_species_mod.F90. Pass these values
             ! here to the species database (ckeller, 11/3/16).
             !
             ! NOTE: EmMw_g will be set to MW_g by default and MolecRatio
             ! will be set to 1 by default, so we can omit setting these
             ! explicitly.  Also the passive species should probably be a gas
             ! instead of an aerosol (i.e., set Is_Gas = T). (bmy, 3/29/17)
             !---------------------------------------------------------------
             IF ( IsPassive ) THEN

                ! Define passive species
                MW_g          = MW_g
                BackgroundVV  = BackgroundVV
                Is_Gas        = F
                Is_Drydep     = F
                Is_Wetdep     = F
                Is_Photolysis = F

             ! Test if this is a non-advected chemical species
             ELSEIF ( KppSpcId > 0 ) THEN

                !------------------------------------------------------------
                ! If this is a non-advected KPP chemical species, then just
                ! create a basic default entry in the species database
                !------------------------------------------------------------
                Is_Gas        = T
                Is_Drydep     = F
                Is_Wetdep     = F
                
             ELSE

                IF ( Present(Found) ) THEN
                   Found = .FALSE.
                ELSE
                   !------------------------------------------------------------
                   ! If this species i not found, the exit with error!
                   ! create a default entry in the species database
                   !------------------------------------------------------------
                   WRITE( 6, '(a)' ) REPEAT( '=', 79 )
                   WRITE( 6, 100 ) TRIM( Name )
                   WRITE( 6, 110 )
                   WRITE( 6, 120 )
100                FORMAT( 'Species ', a, ' not found in the database!'       )
110                FORMAT( 'Please add the information to the CASE statement' )
120                FORMAT( 'in module Headers/species_database_mod.F90!'      )
                   RC = -1
                   RETURN
                ENDIF
             ENDIF

    END SELECT

    ! Pass args
    IF ( PRESENT(oFullName) ) THEN
       oFullName = FullName
       IF ( Uscore ) THEN
          CNT = 0
          IDX = INDEX(TRIM(oFullName),' ')
          DO WHILE ( IDX > 0 ) 
             CNT = CNT + 1
             IF ( CNT > 100 ) EXIT
             oFullName(IDX:IDX) = '_'
             IDX = INDEX(TRIM(oFullName),' ')
          END DO
       ENDIF
    ENDIF

    !=======================================================================
    ! Assign species database values to the optional output arguments
    !=======================================================================
    IF ( PRESENT( oFormula        ) ) oFormula          = Formula
    IF ( PRESENT( oMW_g           ) ) oMW_g             = MW_g
    IF ( PRESENT( oEmMW_g         ) ) oEmMW_g           = EmMW_g
    IF ( PRESENT( oMolecRatio     ) ) oMolecRatio       = MolecRatio
    IF ( PRESENT( oRadius         ) ) oRadius           = Radius
    IF ( PRESENT( oDensity        ) ) oDensity          = Density
    IF ( PRESENT( oDD_AeroDryDep  ) ) oDD_AeroDryDep    = DD_AeroDryDep
    IF ( PRESENT( oDD_DustDryDep  ) ) oDD_DustDryDep    = DD_DustDryDep
    IF ( PRESENT( oDD_F0          ) ) oDD_F0            = DD_F0
    IF ( PRESENT( oDD_DvzAerSnow  ) ) oDD_DvzAerSnow    = DD_DvzAerSnow
    IF ( PRESENT( oDD_DvzMinVal   ) ) oDD_DvzMinVal(:)  = DD_DvzMinVal(:)
    IF ( PRESENT( oDD_KOA         ) ) oDD_KOA           = DD_KOA
    IF ( PRESENT( oDD_Hstar_Old   ) ) oDD_Hstar_Old     = DD_Hstar_Old
    IF ( PRESENT( oHenry_K0       ) ) oHenry_K0         = Henry_K0
    IF ( PRESENT( oHenry_CR       ) ) oHenry_CR         = Henry_CR
    IF ( PRESENT( oHenry_PKa      ) ) oHenry_PKa        = Henry_PKa
    IF ( PRESENT( oWD_RetFactor   ) ) oWD_RetFactor     = WD_RetFactor
    IF ( PRESENT( oWD_LiqAndGas   ) ) oWD_LiqAndGas     = WD_LiqAndGas
    IF ( PRESENT( oWD_ConvFacI2G  ) ) oWD_ConvFacI2G    = WD_ConvFacI2G
    IF ( PRESENT( oWD_AerScavEff  ) ) oWD_AerScavEff    = WD_AerScavEff
    IF ( PRESENT( oWD_KcScaleFac  ) ) oWD_KcScaleFac(:) = WD_KcScaleFac(:)
    IF ( PRESENT( oWD_RainoutEff  ) ) oWD_RainoutEff(:) = WD_RainoutEff(:)
    IF ( PRESENT( oIs_DryDep      ) ) oIs_DryDep        = Is_DryDep
    IF ( PRESENT( oIs_Gas         ) ) oIs_Gas           = Is_Gas
    IF ( PRESENT( oIs_HygroGrowth ) ) oIs_HygroGrowth   = Is_HygroGrowth
    IF ( PRESENT( oIs_Photolysis  ) ) oIs_Photolysis    = Is_Photolysis
    IF ( PRESENT( oIs_WetDep      ) ) oIs_WetDep        = Is_WetDep
    IF ( PRESENT( oBackgroundVV   ) ) oBackgroundVV     = BackgroundVV
    IF ( PRESENT( oIs_InRestart   ) ) oIs_InRestart     = Is_InRestart
    IF ( PRESENT( oWD_CoarseAer   ) ) oWD_CoarseAer     = WD_CoarseAer
    IF ( PRESENT( oMP_SizeResAer  ) ) oMP_SizeResAer    = MP_SizeResAer
    IF ( PRESENT( oMP_SizeResNum  ) ) oMP_SizeResNum    = MP_SizeResNum
    IF ( PRESENT( oIs_Hg0         ) ) oIs_Hg0           = Is_Hg0
    IF ( PRESENT( oIs_Hg2         ) ) oIs_Hg2           = Is_Hg2
    IF ( PRESENT( oIs_HgP         ) ) oIs_HgP           = Is_HgP

  END SUBROUTINE Spc_Info 
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
  SUBROUTINE Cleanup_Species_Database( am_I_Root, SpcData, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE Species_Mod
!
! !INPUT PARAMETERS:
!
    LOGICAL,        INTENT(IN)  :: am_I_Root    ! Are we on the root CPU?
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
!  incur a dependency with GeosUtil/charpak_mod.F.  This lets us
!  keep species_datbase_mod.F90 in the Headers/ folder together
!  with state_chm_mod.F90 and species_mod.F90.
!
! !REVISION HISTORY:
!  06 Jan 2015 - R. Yantosca - Initial version
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
  SUBROUTINE Unique_Species_Names( am_I_Root, Input_Opt, nSpecies, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE Input_Opt_Mod, ONLY : OptInput
    USE GcKpp_Monitor,      ONLY : Spc_Names
    USE GcKpp_Parameters,   ONLY : NFIX, NSPEC, NVAR
!
! !INPUT PARAMETERS:
!
    LOGICAL,        INTENT(IN)  :: am_I_Root   ! Are we on the root CPU?
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
!  02 Aug 2016 - M. Sulprizio- Add KppSpcId; Also only set KppVarId if loop
!                              indexis <= NVAR.
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
!
! !DEFINED PARAMETERS:
!
    ! Missing value
    INTEGER,           PARAMETER   :: MISSING_INT = -999

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
!  05 Jul 2018 - H.P. Lin    - Add missing KppSpcId
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
