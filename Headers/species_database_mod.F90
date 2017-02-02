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
    USE Input_Opt_Mod, ONLY : OptInput
    USE Passive_Tracer_Mod, ONLY : PASSIVE_TRACER_INQUIRE
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
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER             :: C,      N,   nSpecies
    REAL(fp)            :: Radius, KOA, MW_g,    BackgroundVV, HStar
    REAL(f8)            :: K0,     CR

    ! Arrays
    REAL(fp)            :: DvzMinVal(2)
    REAL(fp)            :: KcScale(3)
    REAL(fp)            :: RainEff(3)

    ! Strings
    CHARACTER(LEN=31)   :: NameAllCaps
    CHARACTER(LEN=31)   :: Name
    CHARACTER(LEN=80)   :: FullName

    ! For Tagged Hg species
    INTEGER             :: Hg0_CAT
    INTEGER             :: Hg2_CAT
    INTEGER             :: HgP_CAT

    ! For values from Input_Opt
    LOGICAL             :: Is_Advected
    LOGICAL             :: prtDebug

    ! For passive tracers
    LOGICAL             :: IsPassive
!
! !DEFINED PARAMETERS
!
    LOGICAL,  PARAMETER :: T        = .TRUE.         ! Yes
    LOGICAL,  PARAMETER :: F        = .FALSE.        ! No
    REAL(f8), PARAMETER :: To_M_atm = 9.86923e-3_f8  ! mol/m3/Pa -> M/atm

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

       ! Test for species name
       SELECT CASE( TRIM ( NameAllCaps ) )

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
             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              KppSpcId      = KppSpcId(N),                  &
                              KppVarId      = KppVarId(N),                  &
                              KppFixId      = KppFixId(N),                  &
                              Name          = NameAllCaps,                  &
                              FullName      = 'Acetone',                    &
                              MW_g          = 58.08_fp,                     &
                              EmMW_g        = 12.00_fp,                     &
                              MolecRatio    = 3.0_fp,                       &
                              Is_Advected   = Is_Advected,                  &
                              Is_Gas        = T,                            &
                              Is_Drydep     = T,                            &
                              Is_Wetdep     = F,                            &
                              Is_Photolysis = T,                            &
                              DD_F0         = 1.0_fp,                       &
#if defined( NEW_HENRY_CONSTANTS )                                          
                              Henry_K0      = 2.7e-1_f8 * To_M_atm,         &
                              Henry_CR      = 5500.0_f8,                    &
#else                                                                       
                              DD_Hstar_Old  = 1e5_fp,                       &
                              Henry_K0      = 2.7e+1_f8,                    &
                              Henry_CR      = 5300.0_f8,                    &
#endif                                      
                              RC            = RC )

          CASE( 'ALD2' )
             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              KppSpcId      = KppSpcId(N),                  &
                              KppVarId      = KppVarId(N),                  &
                              KppFixId      = KppFixId(N),                  &
                              Name          = NameAllCaps,                  &
                              FullName      = 'Acetaldehyde',               &
                              MW_g          = 44.05_fp,                     &
                              EmMW_g        = 12.0_fp,                      &
                              MolecRatio    = 2.0_fp,                       &
                              Is_Advected   = Is_Advected,                  &
                              Is_Gas        = T,                            &
                              Is_Drydep     = T,                            &
                              Is_Wetdep     = F,                            &
                              Is_Photolysis = T,                            &
                              DD_F0         = 1.0_fp,                       &
#if defined( NEW_HENRY_CONSTANTS )                                          
                              Henry_K0      = 1.30e-01_f8 * To_M_atm,       &
                              Henry_CR      = 5900.0_f8,                    &
#else                                                                       
                              DD_Hstar_old  = 1.5e+1_fp,                    &
#endif                                      
                              RC            = RC )

          CASE( 'ALK4' )
             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              KppSpcId      = KppSpcId(N),                  &
                              KppVarId      = KppVarId(N),                  &
                              KppFixId      = KppFixId(N),                  &
                              Name          = NameAllCaps,                  &
                              FullName      = 'Lumped >= C4 Alkanes',       &
                              MW_g          = 58.12_fp,                     &
                              EmMW_g        = 12.0_fp,                      &
                              MolecRatio    = 4.0_fp,                       &
                              Is_Advected   = Is_Advected,                  &
                              Is_Gas        = T,                            &
                              Is_Drydep     = F,                            &
                              Is_Wetdep     = F,                            &
#if defined( NEW_HENRY_CONSTANTS )                                          
                              Henry_K0      = 1.20e-5_f8 * To_M_atm,        &
                              Henry_CR      = 3100.0_f8,                    &
#endif                                      
                              RC            = RC )


          CASE( 'ASOA1', 'ASOA2', 'ASOA3', 'ASOAN' )
             FullName = 'Lumped non-volatile aerosol products of light aromatics + IVOCs'
             
             ! Halve the Kc (cloud condensate -> precip) rate
             ! for the temperature range 237 K <= T < 258 K.
             KcScale = (/ 1.0_fp, 0.5_fp, 1.0_fp /)

             ! Turn off rainout only when 237 K <= T < 258K.
             ! NOTE: Rainout efficiency is 0.8 because these are SOA species.
             RainEff = (/ 0.8_fp, 0.0_fp, 0.8_fp /)

             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              KppSpcId      = KppSpcId(N),                  &
                              KppVarId      = KppVarId(N),                  &
                              KppFixId      = KppFixId(N),                  &
                              Name          = NameAllCaps,                  &
                              FullName      = FullName,                     &
                              MW_g          = 150.0_fp,                     &
                              Is_Advected   = Is_Advected,                  &
                              Is_Gas        = F,                            &
                              Is_Drydep     = T,                            &
                              Is_Wetdep     = T,                            &
                              DD_DvzAerSnow = 0.03_fp,                      &
                              DD_F0         = 0.0_fp,                       &
                              DD_Hstar_old  = 0.0_fp,                       &
                              WD_AerScavEff = 0.8_fp,                       &
                              WD_KcScaleFac = KcScale,                      &
                              WD_RainoutEff = RainEff,                      &
                              RC            = RC )

          CASE( 'ASOG1', 'ASOG2', 'ASOG3' )
             FullName = 'Lumped non-volatile gas products of light aromatics + IVOCs'
             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              KppSpcId      = KppSpcId(N),                  &
                              KppVarId      = KppVarId(N),                  &
                              KppFixId      = KppFixId(N),                  &
                              Name          = NameAllCaps,                  &
                              FullName      = FullName,                     &
                              MW_g          = 150.0_fp,                     &
                              Is_Advected   = Is_Advected,                  &
                              Is_Gas        = T,                            &
                              Is_Drydep     = T,                            &
                              Is_Wetdep     = T,                            &
                              DD_F0         = 0.0_fp,                       &
#if defined( NEW_HENRY_CONSTANTS )					    
                              Henry_K0      = 1.00e+5_f8,                   &
                              Henry_CR      = 6039.0_f8,                    &
#else									    
                              DD_Hstar_old  = 1.00e+5_fp,                   &
                              Henry_K0      = 1.00e+5_f8,                   &
                              Henry_CR      = 6039.0_f8,                    &
#endif									    
                              WD_RetFactor  = 2.0e-2_fp,                    &
                              RC            = RC )

          CASE( 'BCPI', 'BCPO' )
             
             ! These have mostly identical properties
             SELECT CASE( NameAllCaps )

                CASE( 'BCPI' )
                   FullName = 'Hydrophilic black carbon aerosol'

                   ! Halve the Kc (cloud condensate -> precip) rate
                   ! for the temperature range 237 K <= T < 258 K.
                   KcScale  = (/ 1.0_fp, 0.5_fp, 1.0_fp /)

                   ! Turn off rainout only when 237 K <= T < 258K.
                   RainEff  = (/ 1.0_fp, 0.0_fp, 1.0_fp /)

                CASE( 'BCPO' )
                   Fullname = 'Hydrophobic black carbon aerosol'

                   ! Halve the Kc (cloud condensate -> precip) rate
                   ! for the temperature range T > 258 K
                   KcScale  = (/ 1.0_fp, 1.0_fp, 0.5_fp /)

                   ! Allow rainout of BCPO when T < 258 K, because
                   ! BCPO is considered to be IN.
                   RainEff  = (/ 1.0_fp, 1.0_fp, 0.0_fp /)

             END SELECT

             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              KppSpcId      = KppSpcId(N),                  &
                              KppVarId      = KppVarId(N),                  &
                              KppFixId      = KppFixId(N),                  &
                              Name          = NameAllCaps,                  &
                              FullName      = FullName,                     &
                              MW_g          = 12.01_fp,                     &
                              EmMW_g        = 12.0_fp,                      &
                              Is_Advected   = Is_Advected,                  &
                              Is_Gas        = F,                            &
                              Is_Drydep     = T,                            &
                              Is_Wetdep     = T,                            &
                              Density       = 1800.0_fp,                    &
                              DD_DvzAerSnow = 0.03_fp,                      &
                              DD_F0         = 0.0_fp,                       &
                              DD_Hstar_Old  = 0.0_fp,                       &
                              WD_AerScavEff = 1.0_fp,                       &
                              WD_KcScaleFac = KcScale,                      &
                              WD_RainoutEff = RainEff,                      &
                              RC            = RC )

          CASE( 'BENZ' )
             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              KppSpcId      = KppSpcId(N),                  &
                              KppVarId      = KppVarId(N),                  &
                              KppFixId      = KppFixId(N),                  &
                              Name          = NameAllCaps,                  &
                              FullName      = 'Benzene',                    &
                              MW_g          = 78.11_fp,                     &
                              EmMW_g        = 12.0_fp,                      &
                              MolecRatio    = 6.0_fp,                       &
                              Is_Advected   = Is_Advected,                  &
                              Is_Gas        = T,                            &
                              Is_Drydep     = F,                            &
                              Is_Wetdep     = F,                            &
                              RC            = RC )

          CASE( 'BR' )
             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              KppSpcId      = KppSpcId(N),                  &
                              KppVarId      = KppVarId(N),                  &
                              KppFixId      = KppFixId(N),                  &
                              Name          = 'Br',                         &
                              FullName      = 'Atomic bromine',             &
                              MW_g          = 80.0_fp,                      &
                              Is_Advected   = Is_Advected,                  &
                              Is_Gas        = T,                            &
                              Is_Drydep     = F,                            &
                              Is_Wetdep     = F,                            &
#if defined( NEW_HENRY_CONSTANTS )					    
                              Henry_K0      = 3.40e-4_f8 * To_M_atm,        &
                              Henry_CR      = 1800.0_f8,                    &
#endif									    
                              RC            = RC )

          CASE( 'BR2' )
             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              KppSpcId      = KppSpcId(N),                  &
                              KppVarId      = KppVarId(N),                  &
                              KppFixId      = KppFixId(N),                  &
                              Name          = 'Br2',                        &
                              FullName      = 'Molecular Bromine',          &
                              MW_g          = 160.0_fp,                     &
                              Is_Advected   = Is_Advected,                  &
                              Is_Gas        = T,                            &
                              Is_Drydep     = T,                            &
                              Is_Wetdep     = T,                            &
                              Is_Photolysis = T,                            &
                              DD_F0         = 0.0_fp,                       &
#if defined( NEW_HENRY_CONSTANTS )					    
                              Henry_K0      = 7.20e-3_f8 * To_M_atm,        &
                              Henry_CR      = 4400.0_f8,                    &
#else									    
                              DD_Hstar_old  = 7.60e-1_fp,                   &
                              Henry_K0      = 7.60e-1_f8,                   &
                              Henry_CR      = 3720.0_f8,                    &
#endif									    
                              WD_RetFactor  = 0.0_fp,                       &
                              RC            = RC )

          CASE( 'BRCL' )
             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              KppSpcId      = KppSpcId(N),                  &
                              KppVarId      = KppVarId(N),                  &
                              KppFixId      = KppFixId(N),                  &
                              Name          = 'BrCl',                       &
                              FullName      = 'Bromine chloride',           &
                              MW_g          = 115.0_fp,                     &
                              Is_Advected   = Is_Advected,                  &
                              Is_Gas        = T,                            &
                              Is_Drydep     = F,                            &
                              Is_Wetdep     = F,                            &
                              Is_Photolysis = T,                            &
                              RC            = RC )

          CASE( 'BRNO2' )
             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              KppSpcId      = KppSpcId(N),                  &
                              KppVarId      = KppVarId(N),                  &
                              KppFixId      = KppFixId(N),                  &
                              Name          = 'BrNO2',                      &
                              FullName      = 'Nitryl bromide',             &
                              MW_g          = 126.0_fp,                     &
                              Is_Advected   = Is_Advected,                  &
                              Is_Gas        = T,                            &
                              Is_Drydep     = F,                            &
                              Is_Wetdep     = F,                            &
                              Is_Photolysis = T,                            &
                              RC            = RC )

          CASE( 'BRNO3' )
             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              KppSpcId      = KppSpcId(N),                  &
                              KppVarId      = KppVarId(N),                  &
                              KppFixId      = KppFixId(N),                  &
                              Name          = 'BrNO3',                      &
                              FullName      = 'Bromine nitrate',            &
                              MW_g          = 142.0_fp,                     &
                              Is_Advected   = Is_Advected,                  &
                              Is_Gas        = T,                            &
                              Is_Drydep     = T,                            &
                              Is_Wetdep     = F,                            &
                              Is_Photolysis = T,                            &
                              DD_F0         = 0.0_fp,                       &
                              DD_Hstar_old  = 1.00e+20_fp,                  &
                              RC            = RC )

          CASE( 'BRO' )
             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              KppSpcId      = KppSpcId(N),                  &
                              KppVarId      = KppVarId(N),                  &
                              KppFixId      = KppFixId(N),                  &
                              Name          = 'BrO',                        &
                              FullName      = 'Bromine monoxide',           &
                              MW_g          = 96.0_fp,                      &
                              Is_Advected   = Is_Advected,                  &
                              Is_Gas        = T,                            &
                              Is_Drydep     = F,                            &
                              Is_Wetdep     = F,                            &
                              Is_Photolysis = T,                            &
                              RC            = RC )

          CASE( 'C2H6' )
             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              KppSpcId      = KppSpcId(N),                  &
                              KppVarId      = KppVarId(N),                  &
                              KppFixId      = KppFixId(N),                  &
                              Name          = NameAllCaps,                  &
                              FullName      = 'Ethane',                     &
                              MW_g          = 30.07_fp,                     &
                              EmMW_g        = 12.0_fp,                      &
                              MolecRatio    = 2.0_fp,                       &
                              Is_Advected   = Is_Advected,                  &
                              Is_Gas        = T,                            &
                              Is_Drydep     = F,                            &
                              Is_Wetdep     = F,                            &
                              DD_F0         = 1.0_fp,                       &
#if defined( NEW_HENRY_CONSTANTS )					    
                              Henry_K0      = 1.90e+5_f8 * To_M_atm,        &
                              Henry_CR      = 2400.0_f8,                    &
#endif									    
                              RC            = RC )

          CASE( 'C3H8' )
             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              KppSpcId      = KppSpcId(N),                  &
                              KppVarId      = KppVarId(N),                  &
                              KppFixId      = KppFixId(N),                  &
                              Name          = NameAllCaps,                  &
                              FullName      = 'Propane',                    &
                              MW_g          = 44.1_fp,                      &
                              EmMW_g        = 12.0_fp,                      &
                              MolecRatio    = 3.0_fp,                       &
                              Is_Advected   = Is_Advected,                  &
                              Is_Gas        = T,                            &
                              Is_Drydep     = F,                            &
                              Is_Wetdep     = F,                            &
#if defined( NEW_HENRY_CONSTANTS )					    
                              Henry_K0      = 1.50e-5_f8 * To_M_atm,        &
                              Henry_CR      = 2700.0_f8,                    &
#endif									    
                              RC            = RC )

          CASE( 'CCL4' )
             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              KppSpcId      = KppSpcId(N),                  &
                              KppVarId      = KppVarId(N),                  &
                              KppFixId      = KppFixId(N),                  &
                              Name          = 'CCl4',                       &
                              FullName      = 'Carbon tetrachloride',       &
                              MW_g          = 152.0_fp,                     &
                              Is_Advected   = Is_Advected,                  &
                              Is_Gas        = T,                            &
                              Is_Drydep     = F,                            &
                              Is_Wetdep     = F,                            &
                              Is_Photolysis = T,                            &
                              RC            = RC )
             
          CASE( 'CFC11' )
             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              KppSpcId      = KppSpcId(N),                  &
                              KppVarId      = KppVarId(N),                  &
                              KppFixId      = KppFixId(N),                  &
                              Name          = NameAllCaps,                  &
                              FullName      = 'CFC-11',                     &
                              MW_g          = 137.0_fp,                     &
                              Is_Advected   = Is_Advected,                  &
                              Is_Gas        = T,                            &
                              Is_Drydep     = F,                            &
                              Is_Wetdep     = F,                            &
                              Is_Photolysis = T,                            &
                              RC            = RC )

          CASE( 'CFC12' )
             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              KppSpcId      = KppSpcId(N),                  &
                              KppVarId      = KppVarId(N),                  &
                              KppFixId      = KppFixId(N),                  &
                              Name          = NameAllCaps,                  &
                              FullName      = 'CFC-12',                     &
                              MW_g          = 121.0_fp,                     &
                              Is_Advected   = Is_Advected,                  &
                              Is_Gas        = T,                            &
                              Is_Drydep     = F,                            &
                              Is_Wetdep     = F,                            &
                              Is_Photolysis = T,                            &
                              RC            = RC )

          CASE( 'CFC113' )
             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              KppSpcId      = KppSpcId(N),                  &
                              KppVarId      = KppVarId(N),                  &
                              KppFixId      = KppFixId(N),                  &
                              Name          = NameAllCaps,                  &
                              FullName      = 'CFC-113',                    &
                              MW_g          = 187.0_fp,                     &
                              Is_Advected   = Is_Advected,                  &
                              Is_Gas        = T,                            &
                              Is_Drydep     = F,                            &
                              Is_Wetdep     = F,                            &
                              Is_Photolysis = T,                            &
                              RC            = RC )

          CASE( 'CFC114' )
             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              KppSpcId      = KppSpcId(N),                  &
                              KppVarId      = KppVarId(N),                  &
                              KppFixId      = KppFixId(N),                  &
                              Name          = NameAllCaps,                  &
                              FullName      = 'CFC-114',                    &
                              MW_g          = 187.0_fp,                     &
                              Is_Advected   = Is_Advected,                  &
                              Is_Gas        = T,                            &
                              Is_Drydep     = F,                            &
                              Is_Wetdep     = F,                            &
                              Is_Photolysis = T,                            &
                              RC            = RC )

          CASE( 'CFC115' )
             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              KppSpcId      = KppSpcId(N),                  &
                              KppVarId      = KppVarId(N),                  &
                              KppFixId      = KppFixId(N),                  &
                              Name          = NameAllCaps,                  &
                              FullName      = 'CFC-115',                    &
                              MW_g          = 187.0_fp,                     &
                              Is_Advected   = Is_Advected,                  &
                              Is_Gas        = T,                            &
                              Is_Drydep     = F,                            &
                              Is_Wetdep     = F,                            &
                              Is_Photolysis = T,                            &
                              RC            = RC )

          CASE( 'CH2BR2' )
             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              KppSpcId      = KppSpcId(N),                  &
                              KppVarId      = KppVarId(N),                  &
                              KppFixId      = KppFixId(N),                  &
                              Name          = 'CH2Br2',                     &
                              FullName      = 'Dibromomethane',             &
                              MW_g          = 174.0_fp,                     &
                              Is_Advected   = Is_Advected,                  &
                              Is_Gas        = T,                            &
                              Is_Drydep     = F,                            &
                              Is_Wetdep     = F,                            &
#if defined( UCX )
                              Is_Photolysis = T,                            &
#endif
                              RC            = RC )

          CASE( 'CH2O', 'HCHO' )
             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              KppSpcId      = KppSpcId(N),                  &
                              KppVarId      = KppVarId(N),                  &
                              KppFixId      = KppFixId(N),                  &
                              Name          = 'CH2O',                       &
                              FullName      = 'Formaldehyde',               &
                              MW_g          = 30.0_fp,                      &
                              BackgroundVV  = 4.0e-15_fp,                   &
                              Is_Advected   = Is_Advected,                  &
                              Is_Gas        = T,                            &
                              Is_Drydep     = T,                            &
                              Is_Wetdep     = T,                            &
                              Is_Photolysis = T,                            &
                              DD_F0         = 1.0_fp,                       &
#if defined( NEW_HENRY_CONSTANTS )					    
                              Henry_K0      = 3.2e+1_f8 * To_M_atm,         &
                              Henry_CR      = 6800.0_f8,                    &
#else									    
                              DD_Hstar_old  = 6.0e+3_fp,                    &
                              Henry_K0      = 3.0e+3_f8,                    &
                              Henry_CR      = 7200.0_f8,                    &
#endif									    
                              WD_RetFactor  = 2.0e-2_fp,                    &
                              RC            = RC )

          CASE( 'CH3BR' )
             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              KppSpcId      = KppSpcId(N),                  &
                              KppVarId      = KppVarId(N),                  &
                              KppFixId      = KppFixId(N),                  &
                              Name          = 'CH3Br',                      &
                              FullName      = 'Methyl bromide',             &
                              MW_g          = 95.0_fp,                      &
                              Is_Advected   = Is_Advected,                  &
                              Is_Gas        = T,                            &
                              Is_Drydep     = F,                            &
                              Is_Wetdep     = F,                            &
#if defined( UCX )
                              Is_Photolysis = T,                            &
#endif
                              RC            = RC )

          CASE( 'CH3CL' )
             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              KppSpcId      = KppSpcId(N),                  &
                              KppVarId      = KppVarId(N),                  &
                              KppFixId      = KppFixId(N),                  &
                              Name          = 'CH3Cl',                      &
                              FullName      = 'Chloromethane',              &
                              MW_g          = 50.0_fp,                      &
                              Is_Advected   = Is_Advected,                  &
                              Is_Gas        = T,                            &
                              Is_Drydep     = F,                            &
                              Is_Wetdep     = F,                            &
                              Is_Photolysis = T,                            &
                              RC            = RC )

          CASE( 'CH3CCL3' )
             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              KppSpcId      = KppSpcId(N),                  &
                              KppVarId      = KppVarId(N),                  &
                              KppFixId      = KppFixId(N),                  &
                              Name          = 'CH3CCl3',                    &
                              FullName      = 'Methyl chloroform',          &
                              MW_g          = 133.0_fp,                     &
                              Is_Advected   = Is_Advected,                  &
                              Is_Gas        = T,                            &
                              Is_Drydep     = F,                            &
                              Is_Wetdep     = F,                            &
                              Is_Photolysis = T,                            &
                              RC            = RC )

          CASE( 'CH4' )
             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              KppSpcId      = KppSpcId(N),                  &
                              KppVarId      = KppVarId(N),                  &
                              KppFixId      = KppFixId(N),                  &
                              Name          = NameAllCaps,                  &
                              FullName      = 'Methane',                    &
                              MW_g          = 16.0_fp,                      &
                              BackgroundVV  = 1.7e-06_fp,                   &
                              Is_Advected   = Is_Advected,                  &
                              Is_Gas        = T,                            &
                              Is_Drydep     = F,                            &
                              Is_Wetdep     = F,                            &
                              RC            = RC )

          CASE( 'CHBR3' )
             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              KppSpcId      = KppSpcId(N),                  &
                              KppVarId      = KppVarId(N),                  &
                              KppFixId      = KppFixId(N),                  &
                              Name          = 'CHBr3',                      &
                              FullName      = 'Bromoform',                  &
                              MW_g          = 253.0_fp,                     &
                              Is_Advected   = Is_Advected,                  &
                              Is_Gas        = T,                            &
                              Is_Drydep     = F,                            &
                              Is_Wetdep     = F,                            &
                              Is_Photolysis = T,                            &
                              RC            = RC )

          CASE( 'CL' )
             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              KppSpcId      = KppSpcId(N),                  &
                              KppVarId      = KppVarId(N),                  &
                              KppFixId      = KppFixId(N),                  &
                              Name          = 'Cl',                         &
                              FullName      = 'Atomic chlorine',            &
                              MW_g          = 35.0_fp,                      &
                              Is_Advected   = Is_Advected,                  &
                              Is_Gas        = T,                            &
                              Is_Drydep     = F,                            &
                              Is_Wetdep     = F,                            &
                              RC            = RC )

          CASE( 'CL2' )
             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              KppSpcId      = KppSpcId(N),                  &
                              KppVarId      = KppVarId(N),                  &
                              KppFixId      = KppFixId(N),                  &
                              Name          = 'Cl2',                        &
                              FullName      = 'Molecular chlorine',         &
                              MW_g          = 71.0_fp,                      &
                              Is_Advected   = Is_Advected,                  &
                              Is_Gas        = T,                            &
                              Is_Drydep     = F,                            &
                              Is_Wetdep     = F,                            &
                              Is_Photolysis = T,                            &
                              RC            = RC )

          CASE( 'CL2O2' )
             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              KppSpcId      = KppSpcId(N),                  &
                              KppVarId      = KppVarId(N),                  &
                              KppFixId      = KppFixId(N),                  &
                              Name          = 'Cl2O2',                      &
                              FullName      = 'Dichlorine dioxide',         &
                              MW_g          = 103.0_fp,                     &
                              Is_Advected   = Is_Advected,                  &
                              Is_Gas        = T,                            &
                              Is_Drydep     = F,                            &
                              Is_Wetdep     = F,                            &
                              Is_Photolysis = T,                            &
                              RC            = RC )

          CASE( 'CLNO2' )
             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              KppSpcId      = KppSpcId(N),                  &
                              KppVarId      = KppVarId(N),                  &
                              KppFixId      = KppFixId(N),                  &
                              Name          = 'ClNO2',                      &
                              FullName      = 'Nitryl chloride',            &
                              MW_g          = 81.0_fp,                      &
                              Is_Advected   = Is_Advected,                  &
                              Is_Gas        = T,                            &
                              Is_Drydep     = F,                            &
                              Is_Wetdep     = F,                            &
                              Is_Photolysis = T,                            &
                              RC            = RC )
          CASE( 'CLNO3' )
             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              KppSpcId      = KppSpcId(N),                  &
                              KppVarId      = KppVarId(N),                  &
                              KppFixId      = KppFixId(N),                  &
                              Name          = 'ClNO3',                      &
                              FullName      = 'Chlorine nitrate',           &
                              MW_g          = 97.0_fp,                      &
                              Is_Advected   = Is_Advected,                  &
                              Is_Gas        = T,                            &
                              Is_Drydep     = F,                            &
                              Is_Wetdep     = F,                            &
                              Is_Photolysis = T,                            &
                              RC            = RC )

          CASE( 'CLO' )
             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              KppSpcId      = KppSpcId(N),                  &
                              KppVarId      = KppVarId(N),                  &
                              KppFixId      = KppFixId(N),                  &
                              Name          = 'ClO',                        &
                              FullName      = 'Chlorine monoxide',          &
                              MW_g          = 51.0_fp,                      &
                              Is_Advected   = Is_Advected,                  &
                              Is_Gas        = T,                            &
                              Is_Drydep     = F,                            &
                              Is_Wetdep     = F,                            &
                              Is_Photolysis = T,                            &
                              RC            = RC )

          CASE( 'CLOO', 'OCLO' )

             ! These have identical properties except for the names
             SELECT CASE( NameAllCaps )
                CASE( 'CLOO' )
                   Name     = 'ClOO'
                CASE( 'OCLO' )
                   Name     = 'OClO'
             END SELECT

             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              KppSpcId      = KppSpcId(N),                  &
                              KppVarId      = KppVarId(N),                  &
                              KppFixId      = KppFixId(N),                  &
                              Name          = Name,                         &
                              FullName      = 'Chlorine dioxide',           &
                              MW_g          = 67.0_fp,                      &
                              Is_Advected   = Is_Advected,                  &
                              Is_Gas        = T,                            &
                              Is_Drydep     = F,                            &
                              Is_Wetdep     = F,                            &
                              Is_Photolysis = T,                            &
                              RC            = RC )

          CASE( 'CO',     'COUS',    'COEUR',  'COASIA', 'COOTH',           &
                'COBBAM', 'COBBAF',  'COBBAS', 'COBBOC', 'COBBEU',          &
                'COBBNA', 'COBBOTH', 'COCH4',  'COBIOF', 'COISOP',          &
                'COMONO', 'COMEOH',  'COACET'                       )

             ! Set Name and LongName for the various CO species
             SELECT CASE( TRIM( NameAllCaps ) )
                CASE( 'CO'     ) 
                   Name     = 'CO'
                   FullName = 'Carbon monoxide'
                CASE( 'COUS'   ) 
                   Name     = 'COus'
                   FullName = 'Anthropogenic + biofuel CO emitted over the USA'
                CASE( 'COEUR'  )
                   Name     = 'COeur'
                   FullName = 'Anthropogenic + biofuel CO emitted over Europe'
                CASE( 'COASIA' ) 
                   Name     = 'COasia'
                   FullName = 'Anthropogenic + biofuel CO emitted over Asia'
                CASE( 'COOTH'  )
                   Name     = 'COoth'
                   FullName = 'Anthropogenic + biofuel CO emitted everywhere else'
                CASE( 'COBBAM' )
                   Name     = 'CObbam'
                   FullName = 'Biomass burning CO emitted over South America'
                CASE( 'COBBAF' )
                   Name     = 'CObbaf'
                   FullName = 'Biomass burning CO emitted over Africa'
                CASE( 'COBBAS' )
                   Name     = 'CObbas'
                   FullName = 'Biomass burning CO emitted over Asia'
                CASE( 'COBBOC' )
                   Name     = 'CObboc'
                   FullName = 'Biomass burning CO emitted over Oceania'
                CASE( 'COBBEU' )
                   Name     = 'CObbeu'
                   FullName = 'Biomass burning CO emitted over Europe'
                CASE( 'COBBOTH' )
                   Name     = 'CObboth'
                   FullName = 'Biomass burning CO emitted everywhere else'
                CASE( 'COCH4'  )
                   Name     = 'COch4'
                   FullName = 'CO produced from methane oxidation'
                CASE( 'COBIOF' )
                   Name     = 'CObiof'
                   FullName = 'CO produced from biofuels (whole world)'
                CASE( 'COISOP' )
                   Name     = 'COisop'
                   FullName = 'CO produced from isoprene oxidation'
                CASE( 'COMONO' )
                   Name     = 'COmono'
                   FullName = 'CO produced from monterpenes oxidation'
                CASE( 'COMEOH' ) 
                   Name     = 'COmeoh'
                   FullName = 'CO produced from methanol oxidation'
                CASE( 'COACET' )
                   Name     = 'COacet'
                   FullName = 'CO produced from acetone oxidation'
             END SELECT

             ! Set special default background for CO
             SELECT CASE( TRIM( NameAllCaps ) )
                CASE( 'CO'     ) 
                   BackgroundVV = 1.0e-07_fp
                CASE DEFAULT
                   BackgroundVV = MISSING_VV
             END SELECT
                                
             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              KppSpcId      = KppSpcId(N),                  &
                              KppVarId      = KppVarId(N),                  &
                              KppFixId      = KppFixId(N),                  &
                              Name          = Name,                         &
                              FullName      = FullName,                     &
                              MW_g          = 28.0_fp,                      &
                              BackgroundVV  = BackgroundVV,                 &
                              Is_Advected   = Is_Advected,                  &
                              Is_Gas        = T,                            &
                              Is_Drydep     = F,                            &
                              Is_Wetdep     = F,                            &
#if defined( NEW_HENRY_CONSTANTS )                                          
                              Henry_K0      = 9.70e-6_f8 * To_M_atm,        &
                              Henry_CR      = 1300.0_f8,                    &
#endif                                      
                              RC            = RC )

          CASE( 'DMS' )
             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              KppSpcId      = KppSpcId(N),                  &
                              KppVarId      = KppVarId(N),                  &
                              KppFixId      = KppFixId(N),                  &
                              Name          = NameAllCaps,                  &
                              FullName      = 'Dimethyl sulfide',           &
                              MW_g          = 62.0_fp,                      &
                              Is_Advected   = Is_Advected,                  &
                              Is_Gas        = F,                            &
                              Is_Drydep     = F,                            &
                              Is_Wetdep     = F,                            &
                              Henry_K0      = 4.80e-1_f8,                   &
                              Henry_CR      = 3100.0_f8,                    &
                              RC            = RC )

          CASE( 'DST1', 'DSTAL1', 'NITD1', 'SO4D1' )
             
             ! These have identical properties except for the names
             SELECT CASE( NameAllCaps )
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

             ! Allow rainout of dust when T < 258K, becasue dust
             ! is considered to be IN.
             RainEff = (/ 1.0_fp, 1.0_fp, 1.0_fp /)

             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              KppSpcId      = KppSpcId(N),                  &
                              KppVarId      = KppVarId(N),                  &
                              KppFixId      = KppFixId(N),                  &
                              Name          = NameAllCaps,                  &
                              FullName      = FullName,                     &
                              MW_g          = 29.0_fp,                      &
                              Is_Advected   = Is_Advected,                  &
                              Is_Gas        = F,                            &
                              Is_Drydep     = T,                            &
                              Is_Wetdep     = T,                            &
                              Density       = 2500.0_fp,                    &
                              Radius        = 7.3e-7_fp,                    &
                              DD_DustDryDep = T,                            &
                              DD_F0         = 0.0_fp,                       &
                              DD_Hstar_Old  = 0.0_fp,                       &
                              WD_AerScavEff = 1.0_fp,                       &
                              WD_KcScaleFac = KcScale,                      &
                              WD_RainoutEff = RainEff,                      &
                              RC            = RC )

          CASE( 'DST2', 'DSTAL2', 'NITD2', 'SO4D2' )
             
             ! These have identical properties except for the names
             SELECT CASE( NameAllCaps )
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
             KcScale = (/ 1.0_fp, 1.0_fp, 1.0_fp /)

             ! Allow rainout of dust when T < 258K, becasue dust
             ! is considered to be IN.
             RainEff = (/ 1.0_fp, 1.0_fp, 1.0_fp /)

             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              KppSpcId      = KppSpcId(N),                  &
                              KppVarId      = KppVarId(N),                  &
                              KppFixId      = KppFixId(N),                  &
                              Name          = NameAllCaps,                  &
                              FullName      = FullName,                     &
                              MW_g          = 29.0_fp,                      &
                              Is_Advected   = Is_Advected,                  &
                              Is_Gas        = F,                            &
                              Is_Drydep     = T,                            &
                              Is_Wetdep     = T,                            &
                              Density       = 2650.0_fp,                    &
                              Radius        = 1.4e-6_fp,                    &
                              DD_DustDryDep = T,                            &
                              DD_F0         = 0.0_fp,                       &
                              DD_Hstar_Old  = 0.0_fp,                       &
                              WD_AerScavEff = 1.0_fp,                       &
                              WD_CoarseAer  = T,                            &
                              WD_KcScaleFac = KcScale,                      &
                              WD_RainoutEff = RainEff,                      &
                              RC            = RC )

          CASE( 'DST3', 'DSTAL3', 'NITD3', 'SO4D3' )
             
             ! These have identical properties except for the names
             SELECT CASE( NameAllCaps )
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
             KcScale = (/ 1.0_fp, 1.0_fp, 1.0_fp /)

             ! Allow rainout of dust when T < 258K, becasue dust
             ! is considered to be IN.
             RainEff = (/ 1.0_fp, 1.0_fp, 1.0_fp /)

             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              KppSpcId      = KppSpcId(N),                  &
                              KppVarId      = KppVarId(N),                  &
                              KppFixId      = KppFixId(N),                  &
                              Name          = NameAllCaps,                  &
                              FullName      = FullName,                     &
                              MW_g          = 29.0_fp,                      &
                              Is_Advected   = Is_Advected,                  &
                              Is_Gas        = F,                            &
                              Is_Drydep     = T,                            &
                              Is_Wetdep     = T,                            &
                              Density       = 2650.0_fp,                    &
                              Radius        = 2.4e-6_fp,                    &
                              DD_DustDryDep = T,                            &
                              DD_F0         = 0.0_fp,                       &
                              DD_Hstar_Old  = 0.0_fp,                       &
                              WD_AerScavEff = 1.0_fp,                       &
                              WD_CoarseAer  = T,                            &
                              WD_KcScaleFac = KcScale,                      &
                              WD_RainoutEff = RainEff,                      &
                              RC            = RC )


          CASE( 'DST4', 'DSTAL4', 'NITD4', 'SO4D4' )
             
             ! These have identical properties except for the names
             SELECT CASE( NameAllCaps )
                CASE( 'DST4' )
                   FullName = 'Dust aerosol, Reff = 4.5 microns'
                CASE( 'DSTAL4' )
                   FullName = 'Dust alkalinity, Reff = 4.5 microns'
                CASE( 'NITD4' )
                   FullName = 'Nitrate on dust, Reff = 4.5 microns'
                CASE( 'SO4D4' )
                   Name     = 'Sulfate on dust, Reff = 4.5 microns'
             END SELECT

             ! Do not reduce the Kc (cloud condensate -> precip) rate
             KcScale = (/ 1.0_fp, 1.0_fp, 1.0_fp /)

             ! Allow rainout of dust when T < 258K, becasue dust
             ! is considered to be IN.
             RainEff = (/ 1.0_fp, 1.0_fp, 1.0_fp /)

             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              KppSpcId      = KppSpcId(N),                  &
                              KppVarId      = KppVarId(N),                  &
                              KppFixId      = KppFixId(N),                  &
                              Name          = NameAllCaps,                  &
                              FullName      = FullName,                     &
                              MW_g          = 29.0_fp,                      &
                              Is_Advected   = Is_Advected,                  &
                              Is_Gas        = F,                            &
                              Is_Drydep     = T,                            &
                              Is_Wetdep     = T,                            &
                              Density       = 2650.0_fp,                    &
                              Radius        = 4.50e-6_fp,                   &
                              DD_DustDryDep = T,                            &
                              DD_F0         = 0.0_fp,                       &
                              DD_Hstar_Old  = 0.0_fp,                       &
                              WD_AerScavEff = 1.0_fp,                       &
                              WD_CoarseAer  = T,                            &
                              WD_KcScaleFac = KcScale,                      &
                              WD_RainoutEff = RainEff,                      &
                              RC            = RC )

          CASE( 'GLYC' )
             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              KppSpcId      = KppSpcId(N),                  &
                              KppVarId      = KppVarId(N),                  &
                              KppFixId      = KppFixId(N),                  &
                              Name          = NameAllCaps,                  &
                              FullName      = 'Glycoaldehyde',              &
                              MW_g          = 60.0_fp,                      &
                              Is_Advected   = Is_Advected,                  &
                              Is_Gas        = T,                            &
                              Is_Drydep     = T,                            &
                              Is_Wetdep     = T,                            &
                              Is_Photolysis = T,                            &
                              DD_F0         = 1.0_fp,                       &
#if defined( NEW_HENRY_CONSTANTS )					    
                              Henry_K0      = 4.10e+2_f8 * To_M_atm,        &
                              Henry_CR      = 4600.0_f8,                    &
#else									    
                              DD_Hstar_old  = 4.10e+4_fp,                   &
                              Henry_K0      = 4.10e+4_f8,                   &
                              Henry_CR      = 4600.0_f8,                    &
#endif									    
                              WD_RetFactor  = 2.0e-2_fp,                    &
                              RC            = RC )

          CASE( 'H2O' )
             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              KppSpcId      = KppSpcId(N),                  &
                              KppVarId      = KppVarId(N),                  &
                              KppFixId      = KppFixId(N),                  &
                              Name          = NameAllCaps,                  &
                              FullName      = 'Water vapor',                &
                              MW_g          = 18.0_fp,                      &
                              BackgroundVV  = 1.839e-02_fp,                 &
                              Is_Advected   = Is_Advected,                  &
                              Is_Gas        = T,                            &
                              Is_Drydep     = F,                            &
                              Is_Wetdep     = F,                            &
                              RC            = RC )

          CASE( 'H2O2' )
             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              KppSpcId      = KppSpcId(N),                  &
                              KppVarId      = KppVarId(N),                  &
                              KppFixId      = KppFixId(N),                  &
                              Name          = NameAllCaps,                  &
                              FullName      = 'Hydrogen peroxide',          &
                              MW_g          = 34.0_fp,                      &
                              BackgroundVV  = 4.0e-15_fp,                   &
                              Is_Advected   = Is_Advected,                  &
                              Is_Gas        = T,                            &
                              Is_Drydep     = T,                            &
                              Is_Wetdep     = T,                            &
                              Is_Photolysis = T,                            &
                              DD_F0         = 1.0_fp,                       &
#if defined( NEW_HENRY_CONSTANTS )                                          
                              Henry_K0      = 4.93e+5_f8 * To_M_atm,        &
                              Henry_CR      = 6600.0_f8,                    &
                              Henry_pKa     = 11.6_f8,                      &
#else                                                                       
                              DD_Hstar_old  = 1e+5_fp,                      &
                              Henry_K0      = 8.30e+4_f8,                   &
                              Henry_CR      = 7400.0_f8,                    &
#endif                                                                      
                              WD_RetFactor  = 5e-2_fp,                      &
                              WD_LiqAndGas  = T,                            &
                              WD_ConvFacI2G = 4.36564e-1_fp,                &
                              RC            = RC )

          CASE( 'HAC' )
             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              KppSpcId      = KppSpcId(N),                  &
                              KppVarId      = KppVarId(N),                  &
                              KppFixId      = KppFixId(N),                  &
                              Name          = NameAllCaps,                  &
                              FullName      = 'Hydroxyacetone',             &
                              MW_g          = 74.0_fp,                      &
                              Is_Advected   = Is_Advected,                  &
                              Is_Gas        = T,                            &
                              Is_Drydep     = T,                            &
                              Is_Wetdep     = F,                            &
                              Is_Photolysis = T,                            &
                              DD_F0         = 1.0_fp,                       &
#if defined( NEW_HENRY_CONSTANTS )					    
                              Henry_K0      = 7.70e+1_f8 * To_M_atm,        &
                              Henry_CR      = 0.0_f8,                       &
#else									    
                              DD_Hstar_old  = 2.90e+3_fp,                   &
#endif									    
                              RC            = RC )

          CASE( 'H1211' )
             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              KppSpcId      = KppSpcId(N),                  &
                              KppVarId      = KppVarId(N),                  &
                              KppFixId      = KppFixId(N),                  &
                              Name          = NameAllCaps,                  &
                              FullName      = 'H-1211',                     &
                              MW_g          = 165.0_fp,                     &
                              Is_Advected   = Is_Advected,                  &
                              Is_Gas        = T,                            &
                              Is_Drydep     = F,                            &
                              Is_Wetdep     = F,                            &
                              Is_Photolysis = T,                            &
                              RC            = RC )

          CASE( 'H1301' )
             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              KppSpcId      = KppSpcId(N),                  &
                              KppVarId      = KppVarId(N),                  &
                              KppFixId      = KppFixId(N),                  &
                              Name          = NameAllCaps,                  &
                              FullName      = 'H-1301',                     &
                              MW_g          = 149.0_fp,                     &
                              Is_Advected   = Is_Advected,                  &
                              Is_Gas        = T,                            &
                              Is_Drydep     = F,                            &
                              Is_Wetdep     = F,                            &
                              Is_Photolysis = T,                            &
                              RC            = RC )

          CASE( 'H2402' )
             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              Name          = NameAllCaps,                  &
                              KppSpcId      = KppSpcId(N),                  &
                              KppVarId      = KppVarId(N),                  &
                              KppFixId      = KppFixId(N),                  &
                              FullName      = 'H-2402',                     &
                              MW_g          = 260.0_fp,                     &
                              Is_Advected   = Is_Advected,                  &
                              Is_Gas        = T,                            &
                              Is_Drydep     = F,                            &
                              Is_Wetdep     = F,                            &
                              Is_Photolysis = T,                            &
                              RC            = RC )

          CASE( 'HCFC123' )

             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              KppSpcId      = KppSpcId(N),                  &
                              KppVarId      = KppVarId(N),                  &
                              KppFixId      = KppFixId(N),                  &
                              Name          = NameAllCaps,                  &
                              FullName      = 'HCFC-123',                   &
                              MW_g          = 117.0_fp,                     &
                              Is_Advected   = Is_Advected,                  &
                              Is_Gas        = T,                            &
                              Is_Drydep     = F,                            &
                              Is_Wetdep     = F,                            &
                              Is_Photolysis = T,                            &
                              RC            = RC )

          CASE( 'HCFC141B' )

             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              KppSpcId      = KppSpcId(N),                  &
                              KppVarId      = KppVarId(N),                  &
                              KppFixId      = KppFixId(N),                  &
                              Name          = 'HCFC141b',                   &
                              FullName      = 'HCFC-141b',                  &
                              MW_g          = 117.0_fp,                     &
                              Is_Advected   = Is_Advected,                  &
                              Is_Gas        = T,                            &
                              Is_Drydep     = F,                            &
                              Is_Wetdep     = F,                            &
                              Is_Photolysis = T,                            &
                              RC            = RC )

          CASE( 'HCFC142B' )

             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              KppSpcId      = KppSpcId(N),                  &
                              KppVarId      = KppVarId(N),                  &
                              KppFixId      = KppFixId(N),                  &
                              Name          = 'HCFC142b',                   &
                              FullName      = 'HCFC-142b',                  &
                              MW_g          = 117.0_fp,                     &
                              Is_Advected   = Is_Advected,                  &
                              Is_Gas        = T,                            &
                              Is_Drydep     = F,                            &
                              Is_Wetdep     = F,                            &
                              Is_Photolysis = T,                            &
                              RC            = RC )

          CASE( 'HCFC22' )
             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              KppSpcId      = KppSpcId(N),                  &
                              KppVarId      = KppVarId(N),                  &
                              KppFixId      = KppFixId(N),                  &
                              Name          = NameAllCaps,                  &
                              FullName      = 'HCFC-22',                    &
                              MW_g          = 86.0_fp,                      &
                              Is_Advected   = Is_Advected,                  &
                              Is_Gas        = T,                            &
                              Is_Drydep     = F,                            &
                              Is_Wetdep     = F,                            &
                              Is_Photolysis = T,                            &
                              RC            = RC )

          CASE( 'HCL' )
             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              KppSpcId      = KppSpcId(N),                  &
                              KppVarId      = KppVarId(N),                  &
                              KppFixId      = KppFixId(N),                  &
                              Name          = 'HCl',                        &
                              FullName      = 'Hydrochloric acid',          &
                              MW_g          = 36.0_fp,                      &
                              Is_Advected   = Is_Advected,                  &
                              Is_Gas        = T,                            &
                              Is_Drydep     = T,                            &
                              Is_Wetdep     = T,                            &
                              DD_F0         = 0.0_fp,                       &
#if defined( NEW_HENRY_CONSTANTS )					    
#else									    
                              DD_Hstar_old  = 2.05e+6_fp,                   &
                              Henry_K0      = 7.10e+15_f8,                  &
                              Henry_CR      = 11000.0_f8,                   &
#endif									    
                              WD_RetFactor  = 1.0_fp,                       &
                              RC            = RC )

          CASE( 'HNO2' )
             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              KppSpcId      = KppSpcId(N),                  &
                              KppVarId      = KppVarId(N),                  &
                              KppFixId      = KppFixId(N),                  &
                              Name          = NameAllCaps,                  &
                              FullName      = 'Nitrous acid',               &
                              MW_g          = 47.0_fp,                      &
                              BackgroundVV  = 4.0e-15_fp,                   &
                              Is_Advected   = Is_Advected,                  &
                              Is_Gas        = T,                            &
                              Is_Drydep     = F,                            &
                              Is_Wetdep     = F,                            &
                              Is_Photolysis = T,                            &
                              RC            = RC )

          CASE( 'HNO3' )

             !%%% NOTE: HNO3 dry-deposits like a gas, but wet-deposits
             !%%% like an aerosol.  Therefore we need to define HNO3 
             !%%% with both gas-phase and aerosol parameters. (bmy, 9/28/15)

             ! Do not reduce the Kc (cloud condensate -> precip) rate
             KcScale = (/ 1.0_fp, 1.0_fp, 1.0_fp /)

             ! Allow rainout of HNO3 when T < 258K, becasue HNO3
             ! is considered to be IN.
             RainEff = (/ 1.0_fp, 1.0_fp, 1.0_fp /)

             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              KppSpcId      = KppSpcId(N),                  &
                              KppVarId      = KppVarId(N),                  &
                              KppFixId      = KppFixId(N),                  &
                              Name          = NameAllCaps,                  &
                              FullName      = 'Nitric acid',                &
                              MW_g          = 63.0_fp,                      &
                              BackgroundVV  = 4.0e-15_fp,                   &
                              Is_Advected   = Is_Advected,                  &
                              Is_Gas        = T,                            &
                              Is_Drydep     = T,                            &
                              Is_Wetdep     = T,                            &
                              Is_Photolysis = T,                            &
                              DD_F0         = 0.0_fp,                       &
#if defined( NEW_HENRY_CONSTANTS )                                          
                              Henry_K0      = 2.10e+3_f8 * To_M_atm,        &
                              Henry_CR      = 8700.0_f8,                    &
#else                                                                       
                              DD_Hstar_old  = 1.0e+14_fp,                   &
                              Henry_K0      = 8.3e+4_f8,                    &
                              Henry_CR      = 7400.0_f8,                    &
#endif                                      
                              WD_AerScavEff = 1.0_fp,                       &
                              WD_KcScaleFac = KcScale,                      &
                              WD_RainoutEff = RainEff,                      &
                              RC            = RC )

          CASE( 'HNO4' )
             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              KppSpcId      = KppSpcId(N),                  &
                              KppVarId      = KppVarId(N),                  &
                              KppFixId      = KppFixId(N),                  &
                              Name          = NameAllCaps,                  &
                              FullName      = 'Pernitric acid',             &
                              MW_g          = 79.0_fp,                      &
                              BackgroundVV  = 4.0e-15_fp,                   &
                              Is_Advected   = Is_Advected,                  &
                              Is_Gas        = T,                            &
                              Is_Drydep     = F,                            &
                              Is_Wetdep     = F,                            &
                              Is_Photolysis = T,                            &
#if defined( NEW_HENRY_CONSTANTS )					    
                              Henry_K0      = 3.90e1X_f8 * To_M_atm,        &
                              Henry_CR      = 8400.0_f8,                    &
#endif									    
                              RC            = RC )

          CASE( 'HBR' )
             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              KppSpcId      = KppSpcId(N),                  &
                              KppVarId      = KppVarId(N),                  &
                              KppFixId      = KppFixId(N),                  &
                              Name          = 'HBr',                        &
                              FullName      = 'Hypobromic acid',            &
                              MW_g          = 81.0_fp,                      &
                              Is_Advected   = Is_Advected,                  &
                              Is_Gas        = T,                            &
                              Is_Drydep     = T,                            &
                              Is_Wetdep     = T,                            &
                              DD_F0         = 0.0_fp,                       &
#if defined( NEW_HENRY_CONSTANTS )					    
                              Henry_K0      = 2.40e-1_f8 * To_M_atm,        &
                              Henry_CR      = 370.0_f8,                     &
#else									    
                              DD_Hstar_old  = 7.10e+15_fp,                  &
                              Henry_K0      = 7.10e+13_f8,                  &
                              Henry_CR      = 10200.0_f8,                   &
#endif									    
                              WD_RetFactor  = 1.0_fp,                       &
                              RC            = RC )

          CASE( 'HOBR' )
             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              KppSpcId      = KppSpcId(N),                  &
                              KppVarId      = KppVarId(N),                  &
                              KppFixId      = KppFixId(N),                  &
                              Name          = 'HOBr',                       &
                              FullName      = 'Hypobromous acid',           &
                              MW_g          = 97.0_fp,                      &
                              Is_Advected   = Is_Advected,                  &
                              Is_Gas        = T,                            &
                              Is_Drydep     = T,                            &
                              Is_Wetdep     = T,                            &
                              Is_Photolysis = T,                            &
                              DD_F0         = 0.0_fp,                       &
#if defined( NEW_HENRY_CONSTANTS )					    
                              Henry_K0      = 1.30e+0_f8 * To_M_atm,        &
                              Henry_CR      = 4000.0_f8,                    &
#else									    
                              DD_Hstar_old  = 6.10e+3_fp,                   &
                              Henry_K0      = 6.10e+3_f8,                   &
                              Henry_CR      = 6014.0_f8,                    &
#endif									    
                              WD_RetFactor  = 0.0_fp,                       &
                              RC            = RC )

          CASE( 'HOCL' )
             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              KppSpcId      = KppSpcId(N),                  &
                              KppVarId      = KppVarId(N),                  &
                              KppFixId      = KppFixId(N),                  &
                              Name          = 'HOCl',                       &
                              FullName      = 'Hypochlorous acid',          &
                              MW_g          = 52.0_fp,                      &
                              Is_Advected   = Is_Advected,                  &
                              Is_Gas        = T,                            &
                              Is_Drydep     = F,                            &
                              Is_Wetdep     = F,                            &
                              Is_Photolysis = T,                            &
                              RC            = RC )

          CASE( 'IEPOX' )
             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              KppSpcId      = KppSpcId(N),                  &
                              KppVarId      = KppVarId(N),                  &
                              KppFixId      = KppFixId(N),                  &
                              Name          = NameAllCaps,                  &
                              FullName      = 'Isoprene epoxide',           &
                              MW_g          = 118.0_fp,                     &
                              Is_Advected   = Is_Advected,                  &
                              Is_Gas        = T,                            &
                              Is_Drydep     = T,                            &
                              Is_Wetdep     = T,                            &
                              DD_F0         = 1.0_fp,                       &
#if defined( NEW_HENRY_CONSTANTS )					    
                              Henry_K0      = 7.60e+5_f8 * To_M_atm,        &
                              Henry_CR      = 0.0_f8,                       &
#else									    
                              DD_Hstar_old  = 1.30e+8_fp,                   &
                              Henry_K0      = 1.30e+8_f8,                   &
                              Henry_CR      = 0.0_f8,                       &
#endif									    
                              WD_RetFactor  = 2.0e-2_fp,                    &
                              RC            = RC )

! Leave for future expansion (bmy, 5/19/16)
!          CASE( 'ISN1' )
!             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
!                              ThisSpc       = SpcData(N)%Info,              &
!                              ModelID       = N,                            &
!                              Name          = NameAllCaps,                  &
!                              FullName      = '',                           &
!                              MW_g          = __.0_fp,                      &
!                              MolecRatio    = 1.0_fp,                       &
!                              Is_Advected   = Is_Advected,                  &
!                              Is_Gas        = T,                            &
!                              Is_Drydep     = T,                            &
!                              Is_Wetdep     = T,                            &
!                              DD_F0         = 0.0_fp,                       &
!                              DD_Hstar_old  = 1.0e+14_fp,                  &
!                              RC            = RC )

          CASE( 'ISOA1', 'ISOA2', 'ISOA3' )
             FullName = 'Lumped semivolatile gas products of isoprene oxidation'

             ! Halve the Kc (cloud condensate -> precip) rate
             ! for the temperature range 237 K <= T < 258 K.
             KcScale = (/ 1.0_fp, 0.5_fp, 1.0_fp /)

             ! Turn off rainout only when 237 K <= T < 258K.
             ! NOTE: Rainout efficiency is 0.8 because these are SOA species.
             RainEff = (/ 0.8_fp, 0.0_fp, 0.8_fp /)

             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              KppSpcId      = KppSpcId(N),                  &
                              KppVarId      = KppVarId(N),                  &
                              KppFixId      = KppFixId(N),                  &
                              Name          = NameAllCaps,                  &
                              FullName      = FullName,                     &
                              MW_g          = 150.0_fp,                     &
                              Is_Advected   = Is_Advected,                  &
                              Is_Gas        = F,                            &
                              Is_Drydep     = T,                            &
                              Is_Wetdep     = T,                            &
                              DD_DvzAerSnow = 0.03_fp,                      &
                              DD_F0         = 0.0_fp,                       &
                              DD_HStar_old  = 0.0_fp,                       &
                              WD_AerScavEff = 0.8_fp,                       &
                              WD_KcScaleFac = KcScale,                      &
                              WD_RainoutEff = RainEff,                      &
                              RC            = RC )

          CASE( 'ISOG1', 'ISOG2', 'ISOG3' )
             FullName = 'Lumped semivolatile gas products of isoprene oxidation'

             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              KppSpcId      = KppSpcId(N),                  &
                              KppVarId      = KppVarId(N),                  &
                              KppFixId      = KppFixId(N),                  &
                              Name          = NameAllCaps,                  &
                              FullName      = FullName,                     &
                              MW_g          = 150.0_fp,                     &
                              Is_Advected   = Is_Advected,                  &
                              Is_Gas        = T,                            &
                              Is_Drydep     = T,                            &
                              Is_Wetdep     = T,                            &
                              DD_F0         = 0.0_fp,                       &
                              DD_Hstar_old  = 1.00e+5_fp,                   &
                              Henry_K0      = 1.00e+5_f8,                   &
                              Henry_CR      = 6039.0_f8,                    &
                              WD_RetFactor  = 2.0e-2_fp,                    &
                              RC            = RC )

          CASE( 'ISOP' )
             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              KppSpcId      = KppSpcId(N),                  &
                              KppVarId      = KppVarId(N),                  &
                              KppFixId      = KppFixId(N),                  &
                              Name          = NameAllCaps,                  &
                              FullName      = 'Isoprene',                   &
                              MW_g          = 68.12_fp,                     &
                              EmMW_g        = 12.0_fp,                      &
                              MolecRatio    = 5e+0_fp,                      &
                              Is_Advected   = Is_Advected,                  &
                              Is_Gas        = T,                            &
                              Is_Drydep     = F,                            &
                              Is_Wetdep     = F,                            &
#if defined( NEW_HENRY_CONSTANTS )                                          
                              Henry_K0      = 3.40e-4_f8 * To_M_atm,        &
                              Henry_CR      = 4400.0_f8,                    &
#endif                                      
                              RC            = RC )

          CASE( 'ISOPNB' )
             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              KppSpcId      = KppSpcId(N),                  &
                              KppVarId      = KppVarId(N),                  &
                              KppFixId      = KppFixId(N),                  &
                              Name          = NameAllCaps,                  &
                              FullName      = 'Isoprene nitrate Beta',      &
                              MW_g          = 147.0_fp,                     &
                              Is_Advected   = Is_Advected,                  &
                              Is_Gas        = T,                            &
                              Is_Drydep     = T,                            &
                              Is_Wetdep     = T,                            &
                              Is_Photolysis = T,                            &
                              DD_F0         = 1.0_fp,                       &
#if defined( NEW_HENRY_CONSTANTS )					    
                              Henry_K0      = 1.97e+4_f8 * To_M_atm,        &
#else									    
                              DD_Hstar_old  = 1.70e+4_fp,                   &
                              Henry_K0      = 1.70e+4_f8,                   &
                              Henry_CR      = 9200.0_f8,                    &
#endif									    
                              WD_RetFactor  = 2.0e-2_fp,                    &
                              RC            = RC )
        
          CASE( 'ISOPND'  )
             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              KppSpcId      = KppSpcId(N),                  &
                              KppVarId      = KppVarId(N),                  &
                              KppFixId      = KppFixId(N),                  &
                              Name          = NameAllCaps,                  &
                              FullName      = 'Isoprene nitrate Delta',     &
                              MW_g          = 147.0_fp,                     &
                              Is_Advected   = Is_Advected,                  &
                              Is_Gas        = T,                            &
                              Is_Drydep     = T,                            &
                              Is_Wetdep     = T,                            &
                              Is_Photolysis = T,                            &
                              DD_F0         = 1.0_fp,                       &
#if defined( NEW_HENRY_CONSTANTS )					    
                              Henry_K0      = 1.97e+4_f8 * To_M_atm,        &
#else									    
                              DD_Hstar_old  = 1.70e+4_fp,                   &
                              Henry_K0      = 1.70e+4_f8,                   &
                              Henry_CR      = 9200.0_f8,                    &
#endif									    
                              WD_RetFactor  = 2.0e-2_fp,                    &
                              RC            = RC )

          CASE( 'LIMO' )
             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              KppSpcId      = KppSpcId(N),                  &
                              KppVarId      = KppVarId(N),                  &
                              KppFixId      = KppFixId(N),                  &
                              Name          = NameAllCaps,                  &
                              FullName      = 'Limonene',                   &
                              MW_g          = 136.23_fp,                    &
                              MolecRatio    = 1.0_fp,                       &
                              Is_Advected   = Is_Advected,                  &
                              Is_Gas        = T,                            &
                              Is_Drydep     = T,                            &
                              Is_Wetdep     = T,                            &
                              DD_F0         = 0.0_fp,                       &
#if defined( NEW_HENRY_CONSTANTS )					    
#else									    
                              DD_Hstar_old  = 7.00e-2_fp,                   &
                              Henry_K0      = 7.00e-2_f8,                   &
                              Henry_CR      = 0.0_f8,                       &
#endif									    
                              WD_RetFactor  = 2.0e-2_fp,                    &
                              RC            = RC )
     
          CASE( 'MACR' )
             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              KppSpcId      = KppSpcId(N),                  &
                              KppVarId      = KppVarId(N),                  &
                              KppFixId      = KppFixId(N),                  &
                              Name          = NameAllCaps,                  &
                              FullName      = 'Methacrolein',               &
                              MW_g          = 70.0_fp,                      &
                              Is_Advected   = Is_Advected,                  &
                              Is_Gas        = T,                            &
                              Is_Drydep     = T,                            &
                              Is_Wetdep     = F,                            &
                              Is_Photolysis = T,                            &
                              DD_F0         = 1.0_fp,                       &
#if defined( NEW_HENRY_CONSTANTS )                                          
                              Henry_K0      = 4.8e-2_f8 * To_M_atm,         &
                              Henry_CR      = 4300.0_f8,                    &
#else                                                                       
                              DD_Hstar_old  = 6.5e+0_fp,                    &
#endif
                              RC            = RC )

          CASE( 'MACRN' )
             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              KppSpcId      = KppSpcId(N),                  &
                              KppVarId      = KppVarId(N),                  &
                              KppFixId      = KppFixId(N),                  &
                              Name          = NameAllCaps,                  &
                              FullName      = 'Nitrate from MACR',          &
                              MW_g          = 149.0_fp,                     &
                              Is_Advected   = Is_Advected,                  &
                              Is_Gas        = T,                            &
                              Is_Drydep     = T,                            &
                              Is_Wetdep     = T,                            &
                              Is_Photolysis = T,                            &
                              DD_F0         = 1.0_fp,                       &
#if defined( NEW_HENRY_CONSTANTS )					    
                              Henry_K0      = 1.97e+4_f8 * To_M_atm,        &
#else									    
                              DD_Hstar_old  = 1.70e+4_fp,                   &
                              Henry_K0      = 1.70e+4_f8,                   &
                              Henry_CR      = 9200.0_f8,                    &
#endif									    
                              WD_RetFactor  = 2.0e-2_fp,                    &
                              RC            = RC )

          CASE( 'MAP' )
             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              KppSpcId      = KppSpcId(N),                  &
                              KppVarId      = KppVarId(N),                  &
                              KppFixId      = KppFixId(N),                  &
                              Name          = NameAllCaps,                  &
                              FullName      = 'Peroxyacetic acid',          &
                              MW_g          = 76.0_fp,                      &
                              Is_Advected   = Is_Advected,                  &
                              Is_Gas        = T,                            &
                              Is_Drydep     = T,                            &
                              Is_Wetdep     = T,                            &
                              Is_Photolysis = T,                            &
                              DD_F0         = 1.0_fp,                       &
#if defined( NEW_HENRY_CONSTANTS )					    
                              Henry_K0      = 8.30e+0_f8 * To_M_atm,        &
                              Henry_CR      = 5300.0_f8,                    &
#else									    
                              DD_Hstar_old  = 8.40e+2_fp,                   &
                              Henry_K0      = 8.40e+2_f8,                   &
                              Henry_CR      = 5300.0_f8,                    &
#endif									    
                              WD_RetFactor  = 2.0e-2_fp,                    &
                              RC            = RC )

          CASE( 'MEK' )
             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              KppSpcId      = KppSpcId(N),                  &
                              KppVarId      = KppVarId(N),                  &
                              KppFixId      = KppFixId(N),                  &
                              Name          = NameAllCaps,                  &
                              FullName      = 'Methyl Ethyl Ketone',        &
                              MW_g          = 72.11_fp,                     &
                              EmMW_g        = 12.0_fp,                      &
                              MolecRatio    = 4.0_fp,                       &
                              Is_Advected   = Is_Advected,                  &
                              Is_Gas        = T,                            &
                              Is_Drydep     = F,                            &
                              Is_Wetdep     = F,                            &
                              Is_Photolysis = T,                            &
#if defined( NEW_HENRY_CONSTANTS )                                          
                              Henry_K0      = 2.90e+02_f8 * To_M_atm,       &
                              Henry_CR      = 5700.0_f8,                    &
#endif                                      
                              RC            = RC )

          CASE( 'MOBA' )
             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              KppSpcId      = KppSpcId(N),                  &
                              KppVarId      = KppVarId(N),                  &
                              KppFixId      = KppFixId(N),                  &
                              Name          = NameAllCaps,                  &
                              FullName      = '5C acid from isoprene',      &
                              MW_g          = 114.0_fp,                     &
                              Is_Advected   = Is_Advected,                  &
                              Is_Gas        = T,                            &
                              Is_Drydep     = F,                            &
                              Is_Wetdep     = T,                            &
#if defined( NEW_HENRY_CONSTANTS )					    
                              Henry_K0      = 2.27e+2_f8 * To_M_atm,        &
                              Henry_CR      = 6300.0_f8,                    &
#else									    
                              Henry_K0      = 2.30e+4_f8,                   &
                              Henry_CR      = 6300.0_f8,                    &
#endif									    
                              WD_RetFactor  = 2.0e-2_fp,                    &
                              RC            = RC )

          CASE( 'MOPI' )
             
             ! Halve the Kc (cloud condensate -> precip) rate
             ! for the temperature range 237 K <= T < 258 K.
             KcScale = (/ 1.0_fp, 0.5_fp, 1.0_fp /)

             ! Turn off rainout only when 237 K <= T < 258K.
             RainEff = (/ 1.0_fp, 0.0_fp, 1.0_fp /)

             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              KppSpcId      = KppSpcId(N),                  &
                              KppVarId      = KppVarId(N),                  &
                              KppFixId      = KppFixId(N),                  &
                              Name          = NameAllCaps,                  &
                              FullName      = 'Hydrophilic marine OC',      &
                              MW_g          = 12.01_fp,                     &
                              EmMW_g        = 12.0_fp,                      &
                              MolecRatio    = 1.0_fp,                       &
                              Is_Advected   = Is_Advected,                  &
                              Is_Gas        = T,                            &
                              Is_Drydep     = F,                            &
                              Is_Wetdep     = F,                            &
                              DD_DvzAerSnow = 0.03_fp,                      &
                              DD_Hstar_old  = 0.0_fp,                       &
                              WD_AerScavEff = 1.0_fp,                       &
                              WD_KcScaleFac = KcScale,                      &
                              WD_RainoutEff = RainEff,                      &
                              RC            = RC )

          CASE( 'MOPO' )

             ! Turn off rainout because MOPO is hydrophobic
             KcScale = (/ 1.0_fp, 1.0_fp, 1.0_fp /)
             RainEff = (/ 0.0_fp, 0.0_fp, 0.0_fp /)

             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              KppSpcId      = KppSpcId(N),                  &
                              KppVarId      = KppVarId(N),                  &
                              KppFixId      = KppFixId(N),                  &
                              Name          = NameAllCaps,                  &
                              FullName      = 'Hydrophobic marine OC',      &
                              MW_g          = 12.01_fp,                     &
                              EmMW_g        = 12.0_fp,                      &
                              MolecRatio    = 1.0_fp,                       &
                              Is_Advected   = Is_Advected,                  &
                              Is_Gas        = T,                            &
                              Is_Drydep     = F,                            &
                              Is_Wetdep     = F,                            &
                              DD_DvzAerSnow = 0.03_fp,                      &
                              DD_Hstar_old  = 0.0_fp,                       &
                              WD_AerScavEff = 0.0_fp,                       &
                              WD_RainOutEff = RainEff,                      &
                              RC            = RC )

          CASE( 'MP', 'CH3OOH' )
             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              KppSpcId      = KppSpcId(N),                  &
                              KppVarId      = KppVarId(N),                  &
                              KppFixId      = KppFixId(N),                  &
                              Name          = 'MP',                         &
                              FullName      = 'Methyl hydro peroxide',      &
                              MW_g          = 48.0_fp,                      &
                              BackgroundVV  = 4.0e-15_fp,                   &
                              Is_Advected   = Is_Advected,                  &
                              Is_Gas        = T,                            &
                              Is_Drydep     = F,                            &
                              Is_Wetdep     = T,                            &
                              Is_Photolysis = T,                            &
#if defined( NEW_HENRY_CONSTANTS )					    
                              Henry_K0      = 2,90e+0_f8 * To_M_atm,        &
                              Henry_CR      = 5200.0_f8,                    &
#else									    
                              Henry_K0      = 3.10e+2_f8,                   &
                              Henry_CR      = 5200.0_f8,                    &
#endif									    
                              WD_RetFactor  = 2.0e-2_fp,                    &
                              RC            = RC )

          CASE( 'MPN' )
             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              KppSpcId      = KppSpcId(N),                  &
                              KppVarId      = KppVarId(N),                  &
                              KppFixId      = KppFixId(N),                  &
                              Name          = NameAllCaps,                  &
                              FullName      = 'Methyl peroxy nitrate',      &
                              MW_g          = 93.0_fp,                      &
                              Is_Advected   = Is_Advected,                  &
                              Is_Gas        = T,                            &
                              Is_Drydep     = F,                            &
                              Is_Wetdep     = F,                            &
                              Is_Photolysis = T,                            &
                              RC            = RC )

          CASE( 'MSA' )

             ! Halve the Kc (cloud condensate -> precip) rate
             ! for the temperature range 237 K <= T < 258 K.
             KcScale   = (/ 1.0_fp, 0.5_fp, 1.0_fp /)

             ! Turn off rainout only when 237 K <= T < 258K.
             RainEff   = (/ 1.0_fp, 0.0_fp, 1.0_fp /)

             ! Enforce minimum dry deposition velocity (Vd) for MSA
             ! (cf. Mian Chin's GOCART model)
             ! Minimum Vd over snow/ice : 0.01 cm/s
             ! Minimum Vd over land     : 0.01 cm/s
             DvzMinVal = (/ 0.01_fp, 0.01_fp /) 

             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              KppSpcId      = KppSpcId(N),                  &
                              KppVarId      = KppVarId(N),                  &
                              KppFixId      = KppFixId(N),                  &
                              Name          = NameAllCaps,                  &
                              FullName      = 'Methyl sulfonic acid',       &
                              MW_g          = 96.0_fp,                      &
                              Is_Advected   = Is_Advected,                  &
                              Is_Gas        = F,                            &
                              Is_Drydep     = T,                            &
                              Is_Wetdep     = T,                            &
                              DD_DvzAerSnow = 0.03_fp,                      &
                              DD_DvzMinVal  = DvzMinVal,                    &
                              DD_F0         = 0.0_fp,                       &
                              DD_Hstar_old  = 0.0_fp,                       &
                              WD_AerScavEff = 1.0_fp,                       &
                              WD_KcScaleFac = KcScale,                      &
                              WD_RainoutEff = RainEff,                      &
                              RC            = RC )

          CASE( 'MTPA' )
             FullName = 'a-pinene, b-pinene, sabinene, carene'
             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              KppSpcId      = KppSpcId(N),                  &
                              KppVarId      = KppVarId(N),                  &
                              KppFixId      = KppFixId(N),                  &
                              Name          = NameAllCaps,                  &
                              FullName      = FullName,                     &
                              MW_g          = 136.23_fp,                    &
                              Is_Advected   = Is_Advected,                  &
                              Is_Gas        = T,                            &
                              Is_Drydep     = T,                            &
                              Is_Wetdep     = T,                            &
                              DD_F0         = 0.0_fp,                       &
#if defined( NEW_HENRY_CONSTANTS )					    
#else									    
                              DD_Hstar_old  = 4.90e-2_fp,                   &
                              Henry_K0      = 4.90e-2_f8,                   &
                              Henry_CR      = 0.0_f8,                       &
#endif									    
                              WD_RetFactor  = 2.0e-2_fp,                    &
                              RC            = RC )

          CASE( 'MTPO' )
             FullName = 'Terpinene, terpinolene, myrcene, ocimene, other monoterpenes'
             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              KppSpcId      = KppSpcId(N),                  &
                              KppVarId      = KppVarId(N),                  &
                              KppFixId      = KppFixId(N),                  &
                              Name          = NameAllCaps,                  &
                              FullName      = FullName,                     &
                              MW_g          = 136.23_fp,                    &
                              MolecRatio    = 1.0_fp,                       &
                              Is_Advected   = Is_Advected,                  &
                              Is_Gas        = T,                            &
                              Is_Drydep     = T,                            &
                              Is_Wetdep     = T,                            &
                              DD_F0         = 0.0_fp,                       &
                              DD_Hstar_old  = 4.90e-2_fp,                   &
                              Henry_K0      = 4.90e-2_f8,                   &
                              Henry_CR      = 0.0_f8,                       &
                              WD_RetFactor  = 2.0e-2_fp,                    &
                              RC            = RC )

          CASE( 'MVK' )
             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              KppSpcId      = KppSpcId(N),                  &
                              KppVarId      = KppVarId(N),                  &
                              KppFixId      = KppFixId(N),                  &
                              Name          = NameAllCaps,                  &
                              FullName      = 'Methyl vinyl ketone',        &
                              MW_g          = 70.0_fp,                      &
                              Is_Advected   = Is_Advected,                  &
                              Is_Gas        = T,                            &
                              Is_Drydep     = T,                            &
                              Is_Wetdep     = F,                            &
                              Is_Photolysis = T,                            &
                              DD_F0         = 1.0_fp,                       &
#if defined( NEW_HENRY_CONSTANTS )                                          
                              Henry_K0      = 2.6e-1_f8 * To_M_atm,         &
                              Henry_CR      = 4800.0_f8,                    &
#else                                                                       
                              DD_Hstar_old  = 4.4e+1_fp,                    &
#endif                                      
                              RC            = RC )

          CASE( 'MVKN' )
             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              KppSpcId      = KppSpcId(N),                  &
                              KppVarId      = KppVarId(N),                  &
                              KppFixId      = KppFixId(N),                  &
                              Name          = NameAllCaps,                  &
                              FullName      = 'Nitrate from MVK',           &
                              MW_g          = 149.0_fp,                     &
                              Is_Advected   = Is_Advected,                  &
                              Is_Gas        = T,                            &
                              Is_Drydep     = T,                            &
                              Is_Wetdep     = T,                            &
                              Is_Photolysis = T,                            &
                              DD_F0         = 1.0_fp,                       &
#if defined( NEW_HENRY_CONSTANTS )					    
                              Henry_K0      = 1.97e+4_f8 * To_M_atm,        &
#else									    
                              DD_Hstar_old  = 1.70e+4_fp,                   &
                              Henry_K0      = 1.70e+4_f8,                   &
                              Henry_CR      = 9200.0_f8,                    &
#endif									    
                              WD_RetFactor  = 2.0e-2_fp,                    &
                              RC            = RC )

          CASE( 'NAP' )
             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              KppSpcId      = KppSpcId(N),                  &
                              KppVarId      = KppVarId(N),                  &
                              KppFixId      = KppFixId(N),                  &
                              Name          = NameAllCaps,                  &
                              FullName      = 'Naphtalene/IVOC surrogate',  &
                              MW_g          = 128.27_fp,                    &
                              EmMw_g        = 12.0_fp,                      &
                              MolecRatio    = 10.0_fp,                      &
                              Is_Advected   = Is_Advected,                  &
                              Is_Gas        = T,                            &
                              Is_Drydep     = F,                            &
                              Is_Wetdep     = F,                            &
                              RC            = RC )

          CASE( 'N2O' )
             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              KppSpcId      = KppSpcId(N),                  &
                              KppVarId      = KppVarId(N),                  &
                              KppFixId      = KppFixId(N),                  &
                              Name          = NameAllCaps,                  &
                              FullName      = 'Nitrous oxide',              &
                              MW_g          = 44.0_fp,                      &
                              BackgroundVV  = 3.0e-07_fp,                   &
                              Is_Advected   = Is_Advected,                  &
                              Is_Gas        = T,                            &
                              Is_Drydep     = F,                            &
                              Is_Wetdep     = F,                            &
                              Is_Photolysis = T,                            &
                              RC            = RC )

          CASE( 'N2O5' )

             ! N2O5 uses the same DD_F0 and DD_Hstar_old values as HNO3,
             ! so that we can compute its drydep velocity explicitly.
             ! (bmy, 5/19/16)
             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              KppSpcId      = KppSpcId(N),                  &
                              KppVarId      = KppVarId(N),                  &
                              KppFixId      = KppFixId(N),                  &
                              Name          = NameAllCaps,                  &
                              FullName      = 'Dinitrogen pentoxide',       &
                              MW_g          = 108.0_fp,                     &
                              BackgroundVV  = 4.0e-15_fp,                   &
                              Is_Advected   = Is_Advected,                  &
                              Is_Gas        = T,                            &
                              Is_Drydep     = T,                            &
                              Is_Wetdep     = F,                            &
                              Is_Photolysis = T,                            &
                              DD_F0         = 0.0_fp,                       &
#if defined( NEW_HENRY_CONSTANTS )					    
                              Henry_K0      = 2.10e-2_f8 * To_M_atm,        &
                              Henry_CR      = 3400.0_f8,                    &
#else									    
                              DD_Hstar_old  = 1.0e+14_fp,                   &
#endif									    
                              RC            = RC )

          CASE( 'NH3' )

             ! Enforce minimum dry deposition velocity (Vd) for NH3
             ! (cf. Mian Chin's GOCART model)
             ! Minimum Vd over snow/ice : 0.2 cm/s
             ! Minimum Vd over land     : 0.3 cm/s
             DvzMinVal = (/ 0.2_fp, 0.3_fp /) 

             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              KppSpcId      = KppSpcId(N),                  &
                              KppVarId      = KppVarId(N),                  &
                              KppFixId      = KppFixId(N),                  &
                              Name          = NameAllCaps,                  &
                              FullName      = 'Ammonia',                    &
                              MW_g          = 17.0_fp,                      &
                              Is_Advected   = Is_Advected,                  &
                              Is_Gas        = T,                            &
                              Is_Drydep     = T,                            &
                              Is_Wetdep     = T,                            &
                              DD_DvzAerSnow = 0.03_fp,                      &
                              DD_DvzMinVal  = DvzMinVal,                    &
                              DD_F0         = 0.0_fp,                       &
#if defined( NEW_HENRY_CONSTANTS )					    
                              Henry_K0      = 5.90e-1_f8 * To_M_atm,        &
                              Henry_CR      = 4200.0_f8,                    &
#else									    
                              DD_Hstar_old  = 2.0e+4_fp,                    &
                              Henry_K0      = 3.30e+6_f8,                   &
                              Henry_CR      = 4100.0_f8,                    &
#endif									    
                              WD_RetFactor  = 5.0e-2_fp,                    &
                              WD_LiqAndGas  = T,                            &
                              WD_ConvFacI2G = 6.17395e-1_fp,                &
                              RC            = RC )

          CASE( 'NH4' )

             ! Halve the Kc (cloud condensate -> precip) rate
             ! for the temperature range 237 K <= T < 258 K.
             KcScale   = (/ 1.0_fp, 0.5_fp, 1.0_fp /)

             ! Turn off rainout only when 237 K <= T < 258K.
             RainEff   = (/ 1.0_fp, 0.0_fp, 1.0_fp /)

             ! Enforce minimum dry deposition velocity (Vd) for NH4
             ! (cf. Mian Chin's GOCART model)
             ! Minimum Vd over snow/ice : 0.01 cm/s
             ! Minimum Vd over land     : 0.01 cm/s
             DvzMinVal = (/ 0.01_fp, 0.01_fp /) 

             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              KppSpcId      = KppSpcId(N),                  &
                              KppVarId      = KppVarId(N),                  &
                              KppFixId      = KppFixId(N),                  &
                              Name          = NameAllCaps,                  &
                              FullName      = 'Ammonium',                   &
                              MW_g          = 18.0_fp,                      &
                              Is_Advected   = Is_Advected,                  &
                              Is_Gas        = F,                            &
                              Is_Drydep     = T,                            &
                              Is_Wetdep     = T,                            &
                              DD_DvzAerSnow = 0.03_fp,                      &
                              DD_DvzMinVal  = DvzMinVal,                    &
                              DD_F0         = 0.0_fp,                       &
                              DD_Hstar_old  = 0.0_fp,                       &
                              WD_AerScavEff = 1.0_fp,                       &
                              WD_KcScaleFac = KcScale,                      &
                              WD_RainoutEff = RainEff,                      &
                              RC            = RC )

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

             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              KppSpcId      = KppSpcId(N),                  &
                              KppVarId      = KppVarId(N),                  &
                              KppFixId      = KppFixId(N),                  &
                              Name          = NameAllCaps,                  &
                              FullName      = 'Inorganic nitrates',         &
                              MW_g          = 62.0_fp,                      &
                              Is_Advected   = Is_Advected,                  &
                              Is_Gas        = F,                            &
                              Is_Drydep     = T,                            &
                              Is_Wetdep     = T,                            &
                              DD_DvzAerSnow = 0.03_fp,                      &
                              DD_DvzMinVal  = DvzMinVal,                    &
                              DD_F0         = 0.0_fp,                       &
                              DD_Hstar_Old  = 0.0_fp,                       &
                              WD_AerScavEff = 1.0_fp,                       &
                              WD_KcScaleFac = KcScale,                      &
                              WD_RainoutEff = RainEff,                      &
                              RC            = RC )

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
             Radius   = ( Input_Opt%SALC_REDGE_um(1) +                      &
                          Input_Opt%SALC_REDGE_um(2)  ) * 0.5e-6_fp

             ! Halve the Kc (cloud condensate -> precip) rate
             ! for the temperature range 237 K <= T < 258 K.
             KcScale  = (/ 1.0_fp, 0.5_fp, 1.0_fp /)

             ! Turn off rainout only when 237 K <= T < 258K.
             RainEff  = (/ 1.0_fp, 0.0_fp, 1.0_fp /)

             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              KppSpcId      = KppSpcId(N),                  &
                              KppVarId      = KppVarId(N),                  &
                              KppFixId      = KppFixId(N),                  &
                              Name          = 'NITs',                       &
                              FullName      = FullName,                     &
                              MW_g          = 31.4_fp,                      &
                              Is_Advected   = Is_Advected,                  &
                              Is_Gas        = F,                            &
                              Is_Drydep     = T,                            &
                              Is_Wetdep     = T,                            &
                              Density       = 2200.0_fp,                    &
                              Radius        = Radius,                       &
                              DD_AeroDryDep = T,                            &
                              DD_F0         = 0.0_fp,                       &
                              DD_Hstar_Old  = 0.0_fp,                       &
                              WD_AerScavEff = 1.0_fp,                       &
                              WD_KcScaleFac = KcScale,                      &
                              WD_RainoutEff = RainEff,                      &
                              RC            = RC )

          CASE( 'NO' )
             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              KppSpcId      = KppSpcId(N),                  &
                              KppVarId      = KppVarId(N),                  &
                              KppFixId      = KppFixId(N),                  &
                              Name          = NameAllCaps,                  &
                              FullName      = 'Nitrogen oxide',             &
                              MW_g          = 30.0_fp,                      &
                              BackgroundVV  = 4.0e-13_fp,                   &
                              Is_Advected   = Is_Advected,                  &
                              Is_Gas        = T,                            &
                              Is_Drydep     = F,                            &
                              Is_Wetdep     = F,                            &
#if defined( UCX )
                              Is_Photolysis = T,                            &
#endif
#if defined( NEW_HENRY_CONSTANTS )                                          
                              Henry_K0      = 1.90e-5_f8 * To_M_atm,        &
                              Henry_CR      = 1600.0_f8,                    &
#endif                                      
                              RC            = RC )

          CASE( 'NO2' )
             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              KppSpcId      = KppSpcId(N),                  &
                              KppVarId      = KppVarId(N),                  &
                              KppFixId      = KppFixId(N),                  &
                              Name          = NameAllCaps,                  &
                              FullName      = 'Nitrogen dioxide',           &
                              MW_g          = 46.0_fp,                      &
                              BackgroundVV  = 4.0e-13_fp,                   &
                              Is_Advected   = Is_Advected,                  &
                              Is_Gas        = T,                            &
                              Is_Drydep     = T,                            &
                              Is_Wetdep     = F,                            &
                              Is_Photolysis = T,                            &
                              DD_F0         = 0.1_fp,                       &
#if defined( NEW_HENRY_CONSTANTS )					    
                              Henry_K0      = 1.20e-4_f8 * To_M_atm,        &
                              Henry_CR      = 2400.0_f8,                    &
#else									    
                              DD_Hstar_old  = 1.00e-2_fp,                   &
#endif									    
                              RC            = RC )

          CASE( 'NO3' )
             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              KppSpcId      = KppSpcId(N),                  &
                              KppVarId      = KppVarId(N),                  &
                              KppFixId      = KppFixId(N),                  &
                              Name          = NameAllCaps,                  &
                              FullName      = 'Nitrate radical',            &
                              MW_g          = 62.0_fp,                      &
                              BackgroundVV  = 4.0e-15_fp,                   &
                              MolecRatio    = 1.0_fp,                       &
                              Is_Advected   = Is_Advected,                  &
                              Is_Gas        = T,                            &
                              Is_Drydep     = F,                            &
                              Is_Wetdep     = F,                            &
                              Is_Photolysis = T,                            &
#if defined( NEW_HENRY_CONSTANTS )					    
                              Henry_K0      = 3.80e-4_f8 * To_M_atm,        &
                              Henry_CR      = 1900.0_f8,                    &
#endif									    
                              RC            = RC )

          CASE( 'O3',     'O3STRAT', 'O3UT',   'O3MT',   'O3ROW',           &
                'O3PCBL', 'O3NABL',  'O3ATBL', 'O3EUBL', 'O3AFBL',          &
                'O3ASBL', 'O3INIT',  'O3USA',  'O3STRT'            )

             ! Now include both total and tagged ozone species (bmy, 10/5/15)
             ! All of these have identical properties except for the names
             SELECT CASE( TRIM( NameAllCaps ) )
                CASE( 'O3' )
                   FullName ='Ozone'
                   Name     = 'O3'
                CASE ( 'O3STRAT', 'O3STRT' )
                   FullName = 'Ozone produced in the stratosphere'
                   Name     = 'O3Strat'
                CASE( 'O3UT' )
                   FullName = 'Ozone produced in the upper tropopshere'
                   Name     = 'O3UT'
                CASE( 'O3MT' )
                   FullName = 'Ozone produced in the middle troposphere'
                   Name     = 'O3MT'
                CASE( 'O3ROW' )
                   FullName = 'Ozone produced in the rest of the world'
                   Name     = 'O3ROW'
                CASE( 'O3PCBL' )
                   FullName = 'Ozone produced in the Pacific Ocean boundary layer'
                   Name     = 'O3PcBL'
                CASE( 'O3NABL' )
                   FullName = 'Ozone produced in the North American boundary layer'
                   Name     = 'O3NaBL'
                CASE( 'O3ATBL' )
                   FullName = 'Ozone produced in the Atlantic Ocean boundary layer'
                   Name     = 'O3AtBL'
                CASE( 'O3EUBL' )
                   FullName = 'Ozone produced in the European boundary layer'
                   Name     = 'O3EuBL'
                CASE( 'O3AFBL' )
                   FullName = 'Ozone produced in the African boundary layer'
                   Name     = 'O3AfBL'
                CASE( 'O3ASBL' )
                   FullName = 'Ozone produced in the Asian boundary layer'
                   Name     = ''
                CASE( 'O3INIT' )
                   FullName = 'Ozone from the initial condition'
                   Name     = 'O3Init'
                CASE( 'O3USA' )
                   FullName = 'Ozone produced over the United States'
                   Name     = 'O3USA'
             END SELECT

             ! Set special default background for O3
             SELECT CASE( TRIM( NameAllCaps ) )
                CASE( 'O3' ) 
                   BackgroundVV = 2.0e-08_fp
                CASE DEFAULT
                   BackgroundVV = MISSING_VV
             END SELECT

             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              KppSpcId      = KppSpcId(N),                  &
                              KppVarId      = KppVarId(N),                  &
                              KppFixId      = KppFixId(N),                  &
                              Name          = Name,                         &
                              FullName      = FullName,                     &
                              MW_g          = 48.0_fp,                      &
                              BackgroundVV  = BackgroundVV,                 &
                              Is_Advected   = Is_Advected,                  &
                              Is_Gas        = T,                            &
                              Is_Drydep     = T,                            &
                              Is_Wetdep     = F,                            &
                              Is_Photolysis = T,                            &
                              DD_F0         = 1.0_fp,                       &
#if defined( NEW_HENRY_CONSTANTS )                                          
                              Henry_K0      = 1.90e-05_f8 * To_M_atm,       &
                              Henry_CR      = 1600.0_f8,                    &
#else                                                                       
                              DD_Hstar_old  = 1.0e-2_fp,                    &
#endif                                      
                              RC            = RC )

          CASE( 'OCPI', 'OCPO' )
             
             ! These have mostly identical properties
             ! Turn off rainout for hydrophobic OC, for all temperatures.
             SELECT CASE( NameAllCaps )

                CASE( 'OCPI' )
                   FullName = 'Hydrophilic organic carbon aerosol'

                   ! Halve the Kc (cloud condensate -> precip) rate
                   ! for the temperature range 237 K <= T < 258 K.
                   KcScale  = (/ 1.0_fp, 0.5_fp, 1.0_fp /)

                   ! Turn off rainout only when 237 K <= T < 258K.
                   RainEff  = (/ 1.0_fp, 0.0_fp, 1.0_fp /)

                CASE( 'OCPO' )
                   Fullname = 'Hydrophobic organic carbon aerosol'

                   ! For all temperatures:
                   ! (1) Halve the Kc (cloud condensate -> precip) rate
                   ! (2) Turn off rainout (OCPO is hydrophobic)
                   KcScale  = (/ 0.5_fp, 0.5_fp, 0.5_fp /)
                   RainEff  = (/ 0.0_fp, 0.0_fp, 0.0_fp /)

             END SELECT

             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              KppSpcId      = KppSpcId(N),                  &
                              KppVarId      = KppVarId(N),                  &
                              KppFixId      = KppFixId(N),                  &
                              Name          = NameAllCaps,                  &
                              FullName      = FullName,                     &
                              MW_g          = 12.01_fp,                     &
                              EmMW_g        = 12.0_fp,                      &
                              Is_Advected   = Is_Advected,                  &
                              Is_Gas        = F,                            &
                              Is_Drydep     = T,                            &
                              Is_Wetdep     = T,                            &
                              Density       = 1300.0_fp,                    &
                              DD_DvzAerSnow = 0.03_fp,                      &
                              DD_F0         = 0.0_fp,                       &
                              DD_Hstar_Old  = 0.0_fp,                       &
                              WD_AerScavEff = 1.0_fp,                       &
                              WD_KcScaleFac = KcScale,                      &
                              WD_RainoutEff = RainEff,                      &
                              RC            = RC )

          CASE( 'OPOA1', 'OPOA2' )
             FullName = 'Lumped aerosol product of SVOC oxidation'

             ! Halve the Kc (cloud condensate -> precip) rate
             ! for the temperature range 237 K <= T < 258 K.
             KcScale  = (/ 1.0_fp, 0.5_fp, 1.0_fp /)

             ! Turn off rainout only when 237 K <= T < 258K.
             ! NOTE: Rainout efficiency is 0.8 because these are SOA species.
             RainEff  = (/ 0.8_fp, 0.0_fp, 0.8_fp /) 

             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              KppSpcId      = KppSpcId(N),                  &
                              KppVarId      = KppVarId(N),                  &
                              KppFixId      = KppFixId(N),                  &
                              Name          = NameAllCaps,                  &
                              FullName      = FullName,                     &
                              MW_g          = 12.01_fp,                     &
                              EmMW_g        = 12.0_fp,                      &
                              MolecRatio    = 1.0_fp,                       &
                              Is_Advected   = Is_Advected,                  &
                              Is_Gas        = F,                            &
                              Is_Drydep     = T,                            &
                              Is_Wetdep     = T,                            &
                              DD_DvzAerSnow = 0.03_fp,                      &
                              DD_F0         = 0.0_fp,                       &
                              DD_Hstar_old  = 0.0_fp,                       &
                              WD_AerScavEff = 0.8_fp,                       &
                              WD_KcScaleFac = KcScale,                      &
                              WD_RainoutEff = RainEff,                      &
                              RC            = RC )

          CASE( 'OPOG1', 'OPOG2' )
             FullName = 'Lumped gas product of SVOC oxidation'
             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              KppSpcId      = KppSpcId(N),                  &
                              KppVarId      = KppVarId(N),                  &
                              KppFixId      = KppFixId(N),                  &
                              Name          = NameAllCaps,                  &
                              FullName      = FullName,                     &
                              MW_g          = 12.01_fp,                     &
                              EmMW_g        = 12.0_fp,                      &
                              MolecRatio    = 1.0_fp,                       &
                              Is_Advected   = Is_Advected,                  &
                              Is_Gas        = T,                            &
                              Is_Drydep     = T,                            &
                              Is_Wetdep     = T,                            &
                              DD_F0         = 0.0_fp,                       &
                              DD_Hstar_old  = 1.00e+5_fp,                   &
                              Henry_K0      = 1.00e+5_f8,                   &
                              Henry_CR      = 6039.0_f8,                    &
                              WD_RetFactor  = 2.0e-2_fp,                    &
                              RC            = RC )

          CASE( 'OCS' )
             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              KppSpcId      = KppSpcId(N),                  &
                              KppVarId      = KppVarId(N),                  &
                              KppFixId      = KppFixId(N),                  &
                              Name          = NameAllCaps,                  &
                              FullName      = 'Carbonyl sulfide',           &
                              MW_g          = 60.0_fp,                      &
                              BackgroundVV  = 9.0e-15_fp,                   &
                              Is_Advected   = Is_Advected,                  &
                              Is_Gas        = T,                            &
                              Is_Drydep     = F,                            &
                              Is_Wetdep     = F,                            &
                              Is_Photolysis = T,                            &
                              RC            = RC )

          CASE( 'PAN' )
             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              KppSpcId      = KppSpcId(N),                  &
                              KppVarId      = KppVarId(N),                  &
                              KppFixId      = KppFixId(N),                  &
                              Name          = NameAllCaps,                  &
                              FullName      = 'Peroxyacetyl nitrate',       &
                              MW_g          = 121.0_fp,                     &
                              Is_Advected   = Is_Advected,                  &
                              Is_Gas        = T,                            &
                              Is_Drydep     = T,                            &
                              Is_Wetdep     = F,                            &
                              Is_Photolysis = T,                            &
                              DD_F0         = 1.0e+0_fp,                    &
#if defined( NEW_HENRY_CONSTANTS )                                          
                              Henry_K0      = 2.90e+02_f8 * To_M_atm,       &
                              Henry_CR      = 5700.0_f8,                    &
#else                                                                       
                              DD_Hstar_old  = 3.60e+0_fp,                   &
#endif                                      
                              RC            = RC )

          CASE( 'PMN' )
             ! PMN uses the same DD_F0 and DD_Hstar_old values as PAN
             ! so that we can compute its drydep velocity explicitly.
             ! (bmy, 5/19/16)
             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              KppSpcId      = KppSpcId(N),                  &
                              KppVarId      = KppVarId(N),                  &
                              KppFixId      = KppFixId(N),                  &
                              Name          = NameAllCaps,                  &
                              FullName      = 'Peroxymethacroyl nitrate',   &
                              MW_g          = 147.0_fp,                     &
                              Is_Advected   = Is_Advected,                  &
                              Is_Gas        = T,                            &
                              Is_Drydep     = T,                            &
                              Is_Wetdep     = F,                            &
                              DD_F0         = 1.0_fp,                       &
#if defined( NEW_HENRY_CONSTANTS )                                          
                              Henry_K0      = 1.70e-2_f8 * To_M_atm,        &
#else                                                                       
                              DD_Hstar_old  = 3.60_fp,                      &
#endif
                              RC            = RC )

          CASE( 'PPN' )
             ! PPN uses the same DD_F0 and DD_Hstar_old values as PAN
             ! so that we can compute its drydep velocity explicitly.
             ! (bmy, 5/19/16)
             FullName = 'Lumped peroxypropionyl nitrate'
             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              KppSpcId      = KppSpcId(N),                  &
                              KppVarId      = KppVarId(N),                  &
                              KppFixId      = KppFixId(N),                  &
                              Name          = NameAllCaps,                  &
                              FullName      = FULLNAME,                     &
                              MW_g          = 135.0_fp,                     &
                              Is_Advected   = Is_Advected,                  &
                              Is_Gas        = T,                            &
                              Is_Drydep     = T,                            &
                              Is_Wetdep     = F,                            &
                              DD_F0         = 1.0_fp,                       &
#if defined( NEW_HENRY_CONSTANTS )                                          
                              Henry_K0      = 2.9e-2_f8 * To_M_atm,         &
                              Henry_CR      = _f8,                          &
#else                                                                       
                              DD_Hstar_old  = 3.60_fp,                      &
#endif                                                                      
                              RC            = RC )

          CASE( 'POA1', 'POA2' )

             ! For all temperatures:
             ! (1) Halve the Kc (cloud condensate -> precip) rate
             ! (2) Turn off rainout (these are hydrophobic species)
             KcScale = (/ 0.5_fp, 0.5_fp, 0.5_fp /)
             RainEff = (/ 0.0_fp, 0.0_fp, 0.0_fp /)

             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              KppSpcId      = KppSpcId(N),                  &
                              KppVarId      = KppVarId(N),                  &
                              KppFixId      = KppFixId(N),                  &
                              Name          = NameAllCaps,                  &
                              FullName      ='Lumped aerosol primary SVOCs',&
                              MW_g          = 12.01_fp,                     &
                              EmMW_g        = 12.0_fp,                      &
                              MolecRatio    = 1.0_fp,                       &
                              Is_Advected   = Is_Advected,                  &
                              Is_Gas        = F,                            &
                              Is_Drydep     = T,                            &
                              Is_Wetdep     = T,                            &
                              Density       = 1300.0_fp,                    &
                              DD_DvzAerSnow = 0.03_fp,                      &
                              DD_F0         = 0.0_fp,                       &
                              DD_Hstar_old  = 0.0_fp,                       &
                              WD_AerScavEff = 1.0_fp,                       &
                              WD_KcScaleFac = KcScale,                      &
                              WD_RainoutEff = RainEff,                      &
                              RC            = RC )

          CASE( 'POG1', 'POG2' )
             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              KppSpcId      = KppSpcId(N),                  &
                              KppVarId      = KppVarId(N),                  &
                              KppFixId      = KppFixId(N),                  &
                              Name          = NameAllCaps,                  &
                              FullName      = 'Lumped gas primary SVOCs',   &
                              MW_g          = 12.01_fp,                     &
                              EmMW_g        = 12.0_fp,                      &
                              MolecRatio    = 1.0_fp,                       &
                              Is_Advected   = Is_Advected,                  &
                              Is_Gas        = T,                            &
                              Is_Drydep     = T,                            &
                              Is_Wetdep     = T,                            &
                              DD_F0         = 0.0_fp,                       &
                              DD_Hstar_old  = 9.50e+0_fp,                   &
                              Henry_K0      = 9.50e+0_f8,                   &
                              Henry_CR      = 4700.0_f8,                    &
                              WD_RetFactor  = 2.0e-2_fp,                    &
                              RC            = RC )

          CASE( 'PRPE' )
             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              KppSpcId      = KppSpcId(N),                  &
                              KppVarId      = KppVarId(N),                  &
                              KppFixId      = KppFixId(N),                  &
                              Name          = NameAllCaps,                  &
                              FullName      = 'Lumped >= C3 alkenes',       &
                              MW_g          = 42.08_fp,                     &
                              EmMW_g        = 12.0_fp,                      &
                              MolecRatio    = 3.0_fp,                       &
                              Is_Advected   = Is_Advected,                  &
                              Is_Gas        = T,                            &
                              Is_Drydep     = F,                            &
                              Is_Wetdep     = F,                            &
#if defined( NEW_HENRY_CONSTANTS )					    
                              Henry_K0      = 7.3e-5_f8 * To_M_atm,         &
                              Henry_CR      = 3400.0_f8,                    &
#endif									    
                              RC            = RC )

          CASE( 'PROPNN' )
             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              KppSpcId      = KppSpcId(N),                  &
                              KppVarId      = KppVarId(N),                  &
                              KppFixId      = KppFixId(N),                  &
                              Name          = NameAllCaps,                  &
                              FullName      = 'Propanone nitrate',          &
                              MW_g          = 119.0_fp,                     &
                              Is_Advected   = Is_Advected,                  &
                              Is_Gas        = T,                            &
                              Is_Drydep     = T,                            &
                              Is_Wetdep     = T,                            &
                              Is_Photolysis = T,                            &
                              DD_F0         = 1.0_fp,                       &
#if defined( NEW_HENRY_CONSTANTS )					    
                              Henry_K0      = 4.93e+3_f8 * To_M_atm,        &
                              Henry_CR      = 0.0_f8,                       &
#else									    
                              DD_Hstar_old  = 1.00e+3_fp,                   &
                              Henry_K0      = 1.00e+3_f8,                   &
                              Henry_CR      = 0.0_f8,                       &
#endif									    
                              WD_RetFactor  = 2.0e-2_fp,                    &
                              RC            = RC )

          CASE( 'R4N2' )
             ! R4N2 uses the same DD_F0 and DD_Hstar_old values as PAN
             ! so that we can compute its drydep velocity explicitly.
             ! (bmy, 5/19/16)
             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              KppSpcId      = KppSpcId(N),                  &
                              KppVarId      = KppVarId(N),                  &
                              KppFixId      = KppFixId(N),                  &
                              Name          = NameAllCaps,                  &
                              FullName      = 'Lumped alkyl nitrate',       &
                              MW_g          = 119.0_fp,                     &
                              Is_Advected   = Is_Advected,                  &
                              Is_Gas        = T,                            &
                              Is_Drydep     = T,                            &
                              Is_Wetdep     = F,                            &
                              Is_Photolysis = T,                            &
                              DD_F0         = 1.0_fp,                       &   
#if defined( NEW_HENRY_CONSTANTS )                                          
                              Henry_K0      = 1.0e-2_f8 * To_M_atm,         &
                              Henry_CR      = 5800.0_f8,                    &
#else                                                                       
                              DD_Hstar_old  = 3.60_fp,                      &
#endif
                              RC            = RC )


          CASE( 'RCHO' )
             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              KppSpcId      = KppSpcId(N),                  &
                              KppVarId      = KppVarId(N),                  &
                              KppFixId      = KppFixId(N),                  &
                              Name          = NameAllCaps,                  &
                              FullName      = 'Lumped aldehyde >= C3',      &
                              MW_g          = 58.0_fp,                      &
                              Is_Advected   = Is_Advected,                  &
                              Is_Gas        = T,                            &
                              Is_Drydep     = F,                            &
                              Is_Wetdep     = F,                            &
                              Is_Photolysis = T,                            &
#if defined( NEW_HENRY_CONSTANTS )                                          
                              Henry_K0      = 9.5e-2_f8 * To_M_atm,         &
                              Henry_CR      = 6200.0_f8,                    &
#endif
                              RC            = RC )

          CASE( 'RIP' )
             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              KppSpcId      = KppSpcId(N),                  &
                              KppVarId      = KppVarId(N),                  &
                              KppFixId      = KppFixId(N),                  &
                              Name          = NameAllCaps,                  &
                              FullName      = 'Peroxide from RIO2',         &
                              MW_g          = 118.0_fp,                     &
                              Is_Advected   = Is_Advected,                  &
                              Is_Gas        = T,                            &
                              Is_Drydep     = T,                            &
                              Is_Wetdep     = T,                            &
                              Is_Photolysis = T,                            &
                              DD_F0         = 1.0_fp,                       &
#if defined( NEW_HENRY_CONSTANTS )					    
                              Henry_K0      = 7.90e+5_f8 * To_M_atm,        &
                              Henry_CR      = 0.0_f8,                       &
#else									    
                              DD_Hstar_old  = 1.70e+6_fp,                   &
                              Henry_K0      = 1.70e+6_f8,                   &
                              Henry_CR      = 0.0_f8,                       &
#endif									    
                              WD_RetFactor  = 2.0e-2_fp,                    &
                              RC            = RC )

          CASE( 'SALA' )
             FullName = 'Accumulation mode sea salt aerosol'
             Radius   = ( Input_Opt%SALA_REDGE_um(1) +                      &
                          Input_Opt%SALA_REDGE_um(2)  ) * 0.5e-6_fp

             ! Halve the Kc (cloud condensate -> precip) rate
             ! for the temperature range 237 K <= T < 258 K.
             KcScale  = (/ 1.0_fp, 0.5_fp, 1.0_fp /)

             ! Turn off rainout only when 237 K <= T < 258K.
             RainEff  = (/ 1.0_fp, 0.0_fp, 1.0_fp /)   

             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              KppSpcId      = KppSpcId(N),                  &
                              KppVarId      = KppVarId(N),                  &
                              KppFixId      = KppFixId(N),                  &
                              Name          = NameAllCaps,                  &
                              FullName      = Fullname,                     &
                              MW_g          = 31.4_fp,                      &
                              Is_Advected   = Is_Advected,                  &
                              Is_Gas        = F,                            &
                              Is_Drydep     = T,                            &
                              Is_Wetdep     = T,                            &
                              Density       = 2200.0_fp,                    &
                              Radius        = Radius,                       &
                              DD_AeroDryDep = T,                            &
                              DD_F0         = 0.0_fp,                       &
                              DD_Hstar_old  = 0.0_fp,                       &
                              WD_AerScavEff = 1.0_fp,                       &
                              WD_KcScaleFac = KcScale,                      &
                              WD_RainoutEff = RainEff,                      &
                              RC            = RC )

          CASE( 'SALC' )
             FullName = 'Coarse mode sea salt aerosol'
             Radius   = ( Input_Opt%SALC_REDGE_um(1) +                      &
                          Input_Opt%SALC_REDGE_um(2)  ) * 0.5e-6_fp
 
             ! Halve the Kc (cloud condensate -> precip) rate
             ! for the temperature range 237 K <= T < 258 K.
             KcScale  = (/ 1.0_fp, 0.5_fp, 1.0_fp /)

             ! Turn off rainout only when 237 K <= T < 258K.
             RainEff  = (/ 1.0_fp, 0.0_fp, 1.0_fp /) 

             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              KppSpcId      = KppSpcId(N),                  &
                              KppVarId      = KppVarId(N),                  &
                              KppFixId      = KppFixId(N),                  &
                              Name          = NameAllCaps,                  &
                              FullName      = Fullname,                     &
                              MW_g          = 31.4_fp,                      &
                              Is_Advected   = Is_Advected,                  &
                              Is_Gas        = F,                            &
                              Is_Drydep     = T,                            &
                              Is_Wetdep     = T,                            &
                              Density       = 2200.0_fp,                    &
                              Radius        = Radius,                       &
                              DD_AeroDryDep = T,                            &
                              DD_F0         = 0.0_fp,                       &
                              DD_Hstar_old  = 0.0_fp,                       &
                              WD_AerScavEff = 1.0_fp,                       &
                              WD_CoarseAer  = T,                            &
                              WD_KcScaleFac = KcScale,                      &
                              WD_RainoutEff = RainEff,                      &
                              RC            = RC )

          CASE( 'SO2' )

             !%%% NOTE: SO2 dry-deposits like a gas but wet-deposits
             !%%% like an aerosol.  Therefore, we need to define SO2 with
             !%%% both gas-phase and aerosol parameters. (bmy, 9/28/15)

             ! Halve the Kc (cloud condensate -> precip) rate
             ! for the temperature range 237 K <= T < 258 K.
             KcScale   = (/ 1.0_fp, 0.5_fp, 1.0_fp /)
 
             ! Turn off rainout only when 237 K <= T < 258K.
             RainEff   = (/ 1.0_fp, 0.0_fp, 1.0_fp /)

             ! Enforce minimum dry deposition velocity (Vd) for SO2
             ! (cf. Mian Chin's GOCART model)
             ! Minimum Vd over snow/ice : 0.2 cm/s
             ! Minimum Vd over land     : 0.3 cm/s
             DvzMinVal = (/ 0.2_fp, 0.3_fp /) 

             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              KppSpcId      = KppSpcId(N),                  &
                              KppVarId      = KppVarId(N),                  &
                              KppFixId      = KppFixId(N),                  &
                              Name          = NameAllCaps,                  &
                              FullName      = 'Sulfur dioxide',             &
                              MW_g          = 64.0_fp,                      &
                              Is_Advected   = Is_Advected,                  &
                              Is_Gas        = T,                            &
                              Is_Drydep     = T,                            &
                              Is_Wetdep     = T,                            &
                              DD_DvzAerSnow = 0.03_fp,                      &
                              DD_DvzMinVal  = DvzMinVal,                    &
                              DD_F0         = 0.0_fp,                       &
#if defined( NEW_HENRY_CONSTANTS )					    
                              Henry_K0      = 1.30e-2_f8 * To_M_atm,        &
                              Henry_CR      = 2900.0_f8,                    &
#else									    
                              DD_Hstar_old  = 1.00e+5_fp,                   &
#endif									    
                              WD_AerScavEff = 1.0_fp,                       &
                              WD_KcScaleFac = KcScale,                      &
                              WD_RainoutEff = RainEff,                      &
                              RC            = RC )

          CASE( 'SO4' )

             ! Halve the Kc (cloud condensate -> precip) rate
             ! for the temperature range 237 K <= T < 258 K.
             KcScale   = (/ 1.0_fp, 0.5_fp, 1.0_fp /)

             ! Turn off rainout only when 237 K <= T < 258K.
             RainEff   = (/ 1.0_fp, 0.0_fp, 1.0_fp /)   

             ! Enforce minimum dry deposition velocity (Vd) for SO4
             ! (cf. Mian Chin's GOCART model)
             ! Minimum Vd over snow/ice : 0.01 cm/s
             ! Minimum Vd over land     : 0.01 cm/s
             DvzMinVal = (/ 0.01_fp, 0.01_fp /) 

             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              KppSpcId      = KppSpcId(N),                  &
                              KppVarId      = KppVarId(N),                  &
                              KppFixId      = KppFixId(N),                  &
                              Name          = 'SO4',                        &
                              FullName      = 'Sulfate',                    &
                              MW_g          = 96.0_fp,                      &
                              Is_Advected   = Is_Advected,                  &
                              Is_Gas        = F,                            &
                              Is_Drydep     = T,                            &
                              Is_Wetdep     = T,                            &
#if defined( UCX )
                              Is_Photolysis = T,                            &
#endif
                              Density       = 1700.0_fp,                    &
                              DD_DvzAerSnow = 0.03_fp,                      &
                              DD_DvzMinVal  = DvzMinVal,                    &
                              DD_F0         = 0.0_fp,                       &
                              DD_Hstar_Old  = 0.0_fp,                       &
                              WD_AerScavEff = 1.0_fp,                       &
                              WD_KcScaleFac = KcScale,                      &
                              WD_RainoutEff = RainEff,                      &
                              RC            = RC )

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
             Radius   = ( Input_Opt%SALC_REDGE_um(1) +                      &
                          Input_Opt%SALC_REDGE_um(2)  ) * 0.5e-6_fp

             ! Halve the Kc (cloud condensate -> precip) rate
             ! for the temperature range 237 K <= T < 258 K.
             KcScale  = (/ 1.0_fp, 0.5_fp, 1.0_fp /)

             ! Turn off rainout only when 237 K <= T < 258K.
             RainEff  = (/ 1.0_fp, 0.0_fp, 1.0_fp /)   

             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              KppSpcId      = KppSpcId(N),                  &
                              KppVarId      = KppVarId(N),                  &
                              KppFixId      = KppFixId(N),                  &
                              Name          = 'SO4s',                       &
                              FullName      = FullName,                     &
                              MW_g          = 31.4_fp,                      &
                              Is_Advected   = Is_Advected,                  &
                              Is_Gas        = F,                            &
                              Is_Drydep     = T,                            &
                              Is_Wetdep     = T,                            &
                              Density       = 2200.0_fp,                    &
                              Radius        = Radius,                       &
                              DD_AeroDryDep = T,                            &
                              DD_F0         = 0.0_fp,                       &
                              DD_Hstar_Old  = 0.0_fp,                       &
                              WD_AerScavEff = 1.0_fp,                       &
                              WD_KcScaleFac = KcScale,                      &
                              WD_RainoutEff = RainEff,                      &
                              RC            = RC )

          CASE( 'TOLU' )
             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              KppSpcId      = KppSpcId(N),                  &
                              KppVarId      = KppVarId(N),                  &
                              KppFixId      = KppFixId(N),                  &
                              Name          = NameAllCaps,                  &
                              FullName      = 'Toluene',                    &
                              MW_g          = 92.14_fp,                     &
                              EmMW_g        = 12.0_fp,                      &
                              MolecRatio    = 7.0_fp,                       &
                              Is_Advected   = Is_Advected,                  &
                              Is_Gas        = T,                            &
                              Is_Drydep     = F,                            &
                              Is_Wetdep     = F,                            &
                              RC            = RC )

          CASE( 'TSOA0', 'TSOA1', 'TSOA2', 'TSOA3' )
             FullName = 'Lumped semivolatile aerosol products of monoterpene + sesquiterpene oxidation'

             ! Halve the Kc (cloud condensate -> precip) rate
             ! for the temperature range 237 K <= T < 258 K.
             KcScale  = (/ 1.0_fp, 0.5_fp, 1.0_fp /)

             ! Turn off rainout only when 237 K <= T < 258K.
             ! NOTE: Rainout efficiency is 0.8 because these are SOA species.
             RainEff  = (/ 0.8_fp, 0.0_fp, 0.8_fp /)   

             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              KppSpcId      = KppSpcId(N),                  &
                              KppVarId      = KppVarId(N),                  &
                              KppFixId      = KppFixId(N),                  &
                              Name          = NameAllCaps,                  &
                              FullName      = FullName,                     &
                              MW_g          = 150.0_fp,                     &
                              Is_Advected   = Is_Advected,                  &
                              Is_Gas        = F,                            &
                              Is_Drydep     = T,                            &
                              Is_Wetdep     = T,                            &
                              DD_DvzAerSnow = 0.03_fp,                      &
                              DD_F0         = 0.0_fp,                       &
                              DD_Hstar_old  = 0.0_fp,                       &
                              WD_AerScavEff = 0.8_fp,                       &
                              WD_KcScaleFac = KcScale,                      &
                              WD_RainoutEff = RainEff,                      &
                              RC            = RC )

          CASE( 'TSOG0', 'TSOG1', 'TSOG2', 'TSOG3' )
             FullName = 'Lumped semivolatile gas products of monoterpene + sesquiterpene oxidation'
             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              KppSpcId      = KppSpcId(N),                  &
                              KppVarId      = KppVarId(N),                  &
                              KppFixId      = KppFixId(N),                  &
                              Name          = NameAllCaps,                  &
                              FullName      = FullName,                     &
                              MW_g          = 150.0_fp,                     &
                              Is_Advected   = Is_Advected,                  &
                              Is_Gas        = T,                            &
                              Is_Drydep     = T,                            &
                              Is_Wetdep     = T,                            &
                              DD_F0         = 0.0_fp,                       &
                              DD_Hstar_old  = 1.00e+5_fp,                   &
                              Henry_K0      = 1.00e+5_f8,                   &
                              Henry_CR      = 6039.0_f8,                    &
                              WD_RetFactor  = 2.0e-2_fp,                    &
                              RC            = RC )

          CASE( 'XYLE' )
             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              KppSpcId      = KppSpcId(N),                  &
                              KppVarId      = KppVarId(N),                  &
                              KppFixId      = KppFixId(N),                  &
                              Name          = NameAllCaps,                  &
                              FullName      = 'Xylene',                     &
                              MW_g          = 106.16_fp,                    &
                              EmMW_g        = 12.0_fp,                      &
                              MolecRatio    = 8.0_fp,                       &
                              Is_Advected   = Is_Advected,                  &
                              Is_Gas        = T,                            &
                              Is_Drydep     = F,                            &
                              Is_Wetdep     = F,                            &
                              RC            = RC )

          !==================================================================
          ! Species for the Rn-Pb-Be specialty simulation
          !==================================================================

          CASE( 'RN', '222RN', 'RN222' )
             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              Name          = 'Rn',                         &
                              FullName      = 'Radon-222 isotope',          &
                              MW_g          = 222.0_fp,                     &
                              Is_Advected   = Is_Advected,                  &
                              Is_Gas        = F,                            &
                              Is_Drydep     = F,                            &
                              Is_Wetdep     = F,                            &
                              RC            = RC )

          CASE( 'PB', '210PB', 'PB210' )

             ! Halve the Kc (cloud condensate -> precip) rate
             ! for the temperature range 237 K <= T < 258 K.
             KcScale = (/ 1.0_fp, 0.5_fp, 1.0_fp /)

             ! Turn off rainout only when 237 K <= T < 258K.
             RainEff = (/ 1.0_fp, 0.0_fp, 1.0_fp /)

             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              Name          = 'Pb',                         &
                              FullName      = 'Lead-210 isotope',           &
                              MW_g          = 210.0_fp,                     &
                              Is_Advected   = Is_Advected,                  &
                              Is_Gas        = F,                            &
                              Is_Drydep     = T,                            &
                              Is_Wetdep     = T,                            &
                              DD_DvzAerSnow = 0.03_fp,                      &
                              DD_F0         = 0.0_fp,                       &
                              DD_HStar_old  = 0.0_fp,                       &
                              WD_AerScavEff = 1.0_fp,                       &
                              WD_KcScaleFac = KcScale,                      &
                              WD_RainoutEff = RainEff,                      &
                              RC            = RC )

          CASE( 'BE', '7BE', 'BE7' )

             ! Halve the Kc (cloud condensate -> precip) rate
             ! for the temperature range 237 K <= T < 258 K.
             KcScale = (/ 1.0_fp, 0.5_fp, 1.0_fp /)

             ! Turn off rainout only when 237 K <= T < 258K.
             RainEff = (/ 1.0_fp, 0.0_fp, 1.0_fp /)

             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              Name          = 'Be7',                        &
                              FullName      = 'Beryllium-7 isotope',        &
                              MW_g          = 7.0_fp,                       &
                              Is_Advected   = Is_Advected,                  &
                              Is_Gas        = F,                            &
                              Is_Drydep     = T,                            &
                              Is_Wetdep     = T,                            &
                              DD_DvzAerSnow = 0.03_fp,                      &
                              DD_F0         = 0.0_fp,                       &
                              DD_HStar_old  = 0.0_fp,                       &
                              WD_AerScavEff = 1.0_fp,                       &
                              WD_KcScaleFac = KcScale,                      &
                              WD_RainoutEff = RainEff,                      &
                              RC            = RC )

          CASE( 'PASV' )
             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              Name          = 'PASV',                       &
                              FullName      = 'Passive species',            &
                              MW_g          = 1.0_fp,                       &
                              Is_Advected   = Is_Advected,                  &
                              Is_Gas        = T,                            &
                              Is_Drydep     = F,                            &
                              Is_Wetdep     = F,                            &
                              RC            = RC )

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
             SELECT CASE( TRIM( NameAllCaps ) ) 
                CASE( 'HG0'     )
                   Name     = 'Hg0'
                   FullName = 'Elemental mercury'
                CASE( 'HG0_CAN' )
                   Name     = 'Hg0_can'
                   FullName = 'Elemental mercury from Canada'
                CASE( 'HG0_USA' )
                   Name     = 'Hg0_usa'
                   FullName = 'Elemental mercury from USA'
                CASE( 'HG0_CAM' )
                   Name     = 'Hg0_cam'
                   FullName = 'Elemental mercury from Central America'
                CASE( 'HG0_SAM' )
                   Name     = 'Hg0_sam'
                   FullName = 'Elemental mercury from South America'
                CASE( 'HG0_WAF' )
                   Name     = 'Hg0_waf'
                   FullName = 'Elemental mercury from West Africa'
                CASE( 'HG0_EAF' ) 
                   Name     = 'Hg0_eaf'
                   FullName = 'Elemental mercury from East Africa'
                CASE( 'HG0_SAF' )
                   Name     = 'Hg0_saf'
                   FullName = 'Elemental mercury from South Africa'
                CASE( 'HG0_NAF' )
                   Name     = 'Hg0_naf'
                   FullName = 'Elemental mercury from North Africa'
                CASE( 'HG0_EUR' )
                   Name     = 'Hg0_eur'
                   FullName = 'Elemental mercury from OECD Europe'
                CASE( 'HG0_EEU' )
                   Name     = 'Hg0_eeu'
                   FullName = 'Elemental mercury from Eastern Europe'
                CASE( 'HG0_MDE' )
                   Name     = 'Hg0_mde'
                   FullName = 'Elemental mercury from Middle East'
                CASE( 'HG0_SOV' )
                   Name     = 'Hg0_sov'
                   FullName = 'Elemental mercury from former USSR'
                CASE( 'HG0_SAS' )
                   Name     = 'Hg0_sas'
                   FullName = 'Elemental mercury from South Asia'
                CASE( 'HG0_EAS' )
                   Name     = 'Hg0_eas'
                   FullName = 'Elemental mercury from East Asia'
                CASE( 'HG0_SEA' )
                   Name     = 'Hg0_sea'
                   FullName = 'Elemental mercury from Southeast Asia'
                CASE( 'HG0_JPN' )
                   Name     = 'Hg0_jpn'
                   FullName = 'Elemental mercury from Japan'
                CASE( 'HG0_OCE' )
                   Name     = 'Hg0_oce'
                   FullName = 'Elemental mercury from Oceania'
                CASE( 'HG0_SO'  )  
                   Name     = 'Hg0_so'
                   FullName = 'Elemental mercury from Organic Soil'
                CASE( 'HG0_BB'  )
                   Name     = 'Hg0_bb'
                   FullName = 'Elemental mercury from Biomass Burning'
                CASE( 'HG0_GEO' )
                   Name     = 'Hg0_geo'
                   FullName = 'Elemental mercury from Geogenic Sources'
                CASE( 'HG0_ATL' )
                   Name     = 'Hg0_atl'
                   FullName = 'Elemental mercury from Midatlantic Subsurface Water'
                CASE( 'HG0_NAT' )
                   Name     = 'Hg0_nat'
                   FullName = 'Elemental mercury from N. Atlantic Subsurface Water'
                CASE( 'HG0_SAT' )
                   Name     = 'Hg0_sat'
                   FullName = 'Elemental mercury from S. Atlantic Subsurface Water'
                CASE( 'HG0_NPA' )
                   Name     = 'Hg0_npa'
                   FullName = 'Elemental mercury from N. Pacific Subsurface Water'
                CASE( 'HG0_ARC' )
                   Name     = 'Hg0_arc'
                   FullName = 'Elemental mercury from Arctic Subsurface Water'
                CASE( 'HG0_ANT' ) 
                   Name     = 'Hg0_ant'
                   FullName = 'Elemental mercury from Antarctic Subsurface Water'
                CASE( 'HG0_OCN' )
                   Name     = 'Hg0_ocn'
                   FullName = 'Elemental mercury from Indo-Pacific Subsurface Water'
                CASE( 'HG0_STR' )
                   Name     = 'Hg0_str'
                   FullName = 'Elemental mercury from Stratosphere'
             END SELECT

             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              Name          = Name,                         &
                              FullName      = FullName,                     &
                              MW_g          = 201.0_fp,                     &
                              Is_Advected   = Is_Advected,                  &
                              Is_Gas        = T,                            &
                              Is_Drydep     = T,                            &
                              Is_Wetdep     = F,                            &
                              Is_Hg0        = T,                            &
                              DD_F0         = 1.0e-5_fp,                    &
                              DD_Hstar_old  = 0.11_fp,                      &
                              RC            = RC )

          CASE( 'HG2',     'HG2_CAN', 'HG2_USA', 'HG2_CAM', 'HG2_SAM',      &
                'HG2_WAF', 'HG2_EAF', 'HG2_SAF', 'HG2_NAF', 'HG2_EUR',      &
                'HG2_EEU', 'HG2_MDE', 'HG2_SOV', 'HG2_SAS', 'HG2_EAS',      &
                'HG2_SEA', 'HG2_JPN', 'HG2_OCE', 'HG2_SO',  'HG2_BB',       &
                'HG2_GEO', 'HG2_ATL', 'HG2_NAT', 'HG2_SAT', 'HG2_NPA',      &
                'HG2_ARC', 'HG2_ANT', 'HG2_OCN', 'HG2_STR'   )

             ! Standardize tagged Hg0 species names 
             SELECT CASE( TRIM( NameAllCaps ) ) 
                CASE( 'HG2'     )
                   Name     = 'Hg2'
                   FullName = 'Divalent mercury'
                CASE( 'HG2_CAN' )
                   Name     = 'Hg2_can'
                   FullName = 'Divalent mercury from Canada'
                CASE( 'HG2_USA' )
                   Name     = 'Hg2_usa'
                   FullName = 'Divalent mercury from USA'
                CASE( 'HG2_CAM' )
                   Name     = 'Hg2_cam'
                   FullName = 'Divalent mercury from Central America'
                CASE( 'HG2_SAM' )
                   Name     = 'Hg2_sam'
                   FullName = 'Divalent mercury from South America'
                CASE( 'HG2_WAF' )
                   Name     = 'Hg2_waf'
                   FullName = 'Divalent mercury from West Africa'
                CASE( 'HG2_EAF' ) 
                   Name     = 'Hg2_eaf'
                   FullName = 'Divalent mercury from East Africa'
                CASE( 'HG2_SAF' )
                   Name     = 'Hg2_saf'
                   FullName = 'Divalent mercury from South Africa'
                CASE( 'HG2_NAF' )
                   Name     = 'Hg2_naf'
                   FullName = 'Divalent mercury from North Africa'
                CASE( 'HG2_EUR' )
                   Name     = 'Hg2_eur'
                   FullName = 'Divalent mercury from OECD Europe'
                CASE( 'HG2_EEU' )
                   Name     = 'Hg2_eeu'
                   FullName = 'Divalent mercury from Eastern Europe'
                CASE( 'HG2_MDE' )
                   Name     = 'Hg2_mde'
                   FullName = 'Divalent mercury from Middle East'
                CASE( 'HG2_SOV' )
                   Name     = 'Hg2_sov'
                   FullName = 'Divalent mercury from former USSR'
                CASE( 'HG2_SAS' )
                   Name     = 'Hg2_sas'
                   FullName = 'Divalent mercury from South Asia'
                CASE( 'HG2_EAS' )
                   Name     = 'Hg2_eas'
                   FullName = 'Divalent mercury from East Asia'
                CASE( 'HG2_SEA' )
                   Name     = 'Hg2_sea'
                   FullName = 'Divalent mercury from Southeast Asia'
                CASE( 'HG2_JPN' )
                   Name     = 'Hg2_jpn'
                   FullName = 'Divalent mercury from Japan'
                CASE( 'HG2_OCE' )
                   Name     = 'Hg2_oce'
                   FullName = 'Divalent mercury from Oceania'
                CASE( 'HG2_SO'  )  
                   Name     = 'Hg2_so'
                   FullName = 'Divalent mercury from Organic Soil'
                CASE( 'HG2_BB'  )
                   Name     = 'Hg2_bb'
                   FullName = 'Divalent mercury from Biomass Burning'
                CASE( 'HG2_GEO' )
                   Name     = 'Hg2_geo'
                   FullName = 'Divalent mercury from Geogenic Sources'
                CASE( 'HG2_ATL' )
                   Name     = 'Hg2_atl'
                   FullName = 'Divalent mercury from Midatlantic Subsurface Water'
                CASE( 'HG2_NAT' )
                   Name     = 'Hg2_nat'
                   FullName = 'Divalent mercury from N. Atlantic Subsurface Water'
                CASE( 'HG2_SAT' )
                   Name     = 'Hg2_sat'
                   FullName = 'Divalent mercury from S. Atlantic Subsurface Water'
                CASE( 'HG2_NPA' )
                   Name     = 'Hg2_npa'
                   FullName = 'Divalent mercury from N. Pacific Subsurface Water'
                CASE( 'HG2_ARC' )
                   Name     = 'Hg2_arc'
                   FullName = 'Divalent mercury from Arctic Subsurface Water'
                CASE( 'HG2_ANT' ) 
                   Name     = 'Hg2_ant'
                   FullName = 'Divalent mercury from Antarctic Subsurface Water'
                CASE( 'HG2_OCN' )
                   Name     = 'Hg2_ocn'
                   FullName = 'Divalent mercury from Indo-Pacific Subsurface Water'
                CASE( 'HG2_STR' )
                   Name     = 'Hg2_str'
                   FullName = 'Divalent mercury from Stratosphere'
             END SELECT

             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              Name          = Name,                         &
                              FullName      = FullName,                     &
                              MW_g          = 201.0_fp,                     &
                              Is_Advected   = Is_Advected,                  &
                              Is_Gas        = T,                            &
                              Is_Drydep     = T,                            &
                              Is_Wetdep     = T,                            &
                              Is_Hg2        = T,                            &
                              DD_F0         = 0.0_fp,                       &
#if defined( NEW_HENRY_CONSTANTS )					    
                              Henry_K0      = 1.40e+4_f8 * To_M_atm,        &
                              Henry_CR      = 5300.0_f8,                    &
#else									    
                              DD_Hstar_old  = 1.00e+14_fp,                  &
                              Henry_K0      = 1.40e+6_f8,                   &
                              Henry_CR      = 8400.0_f8,                    &
#endif									    
                              WD_RetFactor  = 1.0_fp,                       &
                              RC            = RC )

          CASE( 'HGP',     'HGP_CAN', 'HGP_USA', 'HGP_CAM', 'HGP_SAM',      &
                'HGP_WAF', 'HGP_EAF', 'HGP_SAF', 'HGP_NAF', 'HGP_EUR',      &
                'HGP_EEU', 'HGP_MDE', 'HGP_SOV', 'HGP_SAS', 'HGP_EAS',      &
                'HGP_SEA', 'HGP_JPN', 'HGP_OCE', 'HGP_SO',  'HGP_BB',       &
                'HGP_GEO', 'HGP_ATL', 'HGP_NAT', 'HGP_SAT', 'HGP_NPA',      &
                'HGP_ARC', 'HGP_ANT', 'HGP_OCN', 'HGP_STR' )

             ! Standardize tagged HgP species names 
             SELECT CASE( TRIM( NameAllCaps ) )
                 CASE( 'HGP'     )
                   Name     = 'HgP'
                   FullName = 'Particulate mercury'
                CASE( 'HGP_CAN' )
                   Name     = 'HgP_can'
                   FullName = 'Particulate mercury from Canada'
                CASE( 'HGP_USA' )
                   Name     = 'HgP_usa'
                   FullName = 'Particulate mercury from USA'
                CASE( 'HGP_CAM' )
                   Name     = 'HgP_cam'
                   FullName = 'Particulate mercury from Central America'
                CASE( 'HGP_SAM' )
                   Name     = 'HgP_sam'
                   FullName = 'Particulate mercury from South America'
                CASE( 'HGP_WAF' )
                   Name     = 'HgP_waf'
                   FullName = 'Particulate mercury from West Africa'
                CASE( 'HGP_EAF' ) 
                   Name     = 'HgP_eaf'
                   FullName = 'Particulate mercury from East Africa'
                CASE( 'HGP_SAF' )
                   Name     = 'HgP_saf'
                   FullName = 'Particulate mercury from South Africa'
                CASE( 'HGP_NAF' )
                   Name     = 'HgP_naf'
                   FullName = 'Particulate mercury from North Africa'
                CASE( 'HGP_EUR' )
                   Name     = 'HgP_eur'
                   FullName = 'Particulate mercury from OECD Europe'
                CASE( 'HGP_EEU' )
                   Name     = 'HgP_eeu'
                   FullName = 'Particulate mercury from Eastern Europe'
                CASE( 'HGP_MDE' )
                   Name     = 'HgP_mde'
                   FullName = 'Particulate mercury from Middle East'
                CASE( 'HGP_SOV' )
                   Name     = 'HgP_sov'
                   FullName = 'Particulate mercury from former USSR'
                CASE( 'HGP_SAS' )
                   Name     = 'HgP_sas'
                   FullName = 'Particulate mercury from South Asia'
                CASE( 'HGP_EAS' )
                   Name     = 'HgP_eas'
                   FullName = 'Particulate mercury from East Asia'
                CASE( 'HGP_SEA' )
                   Name     = 'HgP_sea'
                   FullName = 'Particulate mercury from Southeast Asia'
                CASE( 'HGP_JPN' )
                   Name     = 'HgP_jpn'
                   FullName = 'Particulate mercury from Japan'
                CASE( 'HGP_OCE' )
                   Name     = 'HgP_oce'
                   FullName = 'Particulate mercury from Oceania'
                CASE( 'HGP_SO'  )  
                   Name     = 'HgP_so'
                   FullName = 'Particulate mercury from Organic Soil'
                CASE( 'HGP_BB'  )
                   Name     = 'HgP_bb'
                   FullName = 'Particulate mercury from Biomass Burning'
                CASE( 'HGP_GEO' )
                   Name     = 'HgP_geo'
                   FullName = 'Particulate mercury from Geogenic Sources'
                CASE( 'HGP_ATL' )
                   Name     = 'HgP_atl'
                   FullName = 'Particulate mercury from Midatlantic Subsurface Water'
                CASE( 'HGP_NAT' )
                   Name     = 'HgP_nat'
                   FullName = 'Particulate mercury from N. Atlantic Subsurface Water'
                CASE( 'HGP_SAT' )
                   Name     = 'HgP_sat'
                   FullName = 'Particulate mercury from S. Atlantic Subsurface Water'
                CASE( 'HGP_NPA' )
                   Name     = 'HgP_npa'
                   FullName = 'Particulate mercury from N. Pacific Subsurface Water'
                CASE( 'HGP_ARC' )
                   Name     = 'HgP_arc'
                   FullName = 'Particulate mercury from Arctic Subsurface Water'
                CASE( 'HGP_ANT' ) 
                   Name     = 'HgP_ant'
                   FullName = 'Particulate mercury from Antarctic Subsurface Water'
                CASE( 'HGP_OCN' )
                   Name     = 'HgP_ocn'
                   FullName = 'Particulate mercury from Indo-Pacific Subsurface Water'
                CASE( 'HGP_STR' )
                   Name     = 'HgP_str'
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
             KcScale = (/ 1.0_fp, 0.5_fp, 1.0_fp /)
             RainEff = (/ 1.0_fp, 1.0_fp, 1.0_fp /)

             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              Name          = Name,                         &
                              FullName      = FullName,                     &
                              MW_g          = 201.0_fp,                     &
                              Is_Advected   = Is_Advected,                  &
                              Is_Gas        = F,                            &
                              Is_Drydep     = T,                            &
                              Is_Wetdep     = T,                            &
                              Is_HgP        = T,                            &
                              DD_DvzAerSnow = 0.03_fp,                      &
                              DD_F0         = 0.0_fp,                       &
                              DD_Hstar_old  = 0.0_fp,                       &
                              WD_AerScavEff = 1.0_fp,                       &
                              WD_KcScaleFac = KcScale,                      &
                              WD_RainoutEff = RainEff,                      &
                              RC            = RC )

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

             SELECT CASE( TRIM( Input_Opt%POP_TYPE ) )
                CASE( 'PHE' )
                   FullName = 'Phenanthrene (gas phase)'
                   MW_g     = 178.23_fp
                   KOA      = 4.37e+7_fp  * 0.0409_fp  * 0.8_fp
                   Hstar    = 1.0_fp      / 1.74e-3_fp * 0.0409_fp
                   K0       = 1.0_f8      / 1.74e-3_f8 / 8.21e-2_f8 / 298.0_f8
                   CR       = 47.0_f8     / 8.32e-3_f8
                CASE( 'PYR' )
                   FullName = 'Pyrene (gas phase)'
                   MW_g     = 202.25_fp
                   KOA      = 7.24e+8_fp  * 0.0409_fp  * 0.8_fp
                   Hstar    = 1.0_fp      / 5.37e-4_fp * 0.0409_fp
                   K0       = 1.0_f8      / 5.37e-4_f8 / 8.21e-2_f8 / 298.0_f8
                   CR       = 43.0_f8     / 8.32e-3_f8
                CASE( 'BaP' )
                   FullName = 'Benzo(a)pyrene (gas phase)' 
                   MW_g     = 252.31_fp
                   KOA      = 3.02e+11_fp * 0.0409_fp  * 0.8_fp
                   Hstar    = 1.0_fp      / 3.10e-5_fp * 0.0409_fp
                   K0       = 1.0_f8      / 3.10e-5_f8 / 8.21e-2_f8 / 298.0_f8
                   CR       = 43.0_f8     / 8.32e-3_f8
             END SELECT

             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              Name          = NameAllCaps,                  &
                              FullName      = FullName,                     &
                              MW_g          = MW_g,                         &
                              Is_Advected   = Is_Advected,                  &
                              Is_Gas        = T,                            &
                              Is_Drydep     = T,                            &
                              Is_Wetdep     = T,                            &
                              DD_F0         = 0.0_fp,                       &
                              DD_Hstar_old  = Hstar,                        &
                              DD_KOA        = KOA,                          &
                              Henry_K0      = K0,                           &
                              Henry_CR      = CR,                           &
                              WD_RetFactor  = 0.0_fp,                       &
                              RC            = RC )

          CASE( 'POPPBCPO', 'POPPOCPO' )

             ! These have identical properties except for the mol weights ...
             SELECT CASE( TRIM( Input_Opt%POP_TYPE ) )
                CASE( 'PHE' )
                   MW_g     = 178.23_fp
                   FullName = 'Phenanthrene particles on'
                CASE( 'PYR' )
                   MW_g     = 202.25_fp
                   FullName = 'Pyrene particles on'
                CASE( 'BaP' )
                   MW_g     = 252.31_fp
                   FullName = 'Benzo(a)pyrene particles on'
             END SELECT

             ! ... and the names and rainout efficiencies.
             SELECT CASE( TRIM( NameAllCaps ) ) 

                CASE( 'POPPBCPO' ) 
                   FullName = TRIM( FullName ) // ' hydrophobic black carbon'

                   ! Halve the Kc (cloud condensate -> precip) rate
                   ! for the temperature range 237 K <= T < 258 K.
                   KcScale  = (/ 1.0_fp, 1.0_fp, 0.5_fp /)

                   ! Allow rainout of POPPBCPO when T < 258 K, because
                   ! POPPBCPO is considered to be IN.
                   RainEff  = (/ 1.0_fp, 1.0_fp, 0.0_fp /)

                CASE( 'POPPOCPO' )
                   FullName = TRIM( FullName ) // ' hydrophobic organic carbon'

                   ! For all temperatures:
                   ! (1) Halve the Kc (cloud condensate -> precip) rate
                   ! (2) Turn off rainout (it's hydrophobic)
                   KcScale  = (/ 0.5_fp, 0.5_fp, 0.5_fp /)
                   RainEff  = (/ 0.0_fp, 0.0_fp, 0.0_fp /)

             END SELECT

             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              Name          = NameAllCaps,                  &
                              FullName      = FullName,                     &
                              MW_g          = MW_g,                         &
                              Is_Advected   = Is_Advected,                  &
                              Is_Gas        = F,                            &
                              Is_Drydep     = T,                            &
                              Is_Wetdep     = T,                            &
                              DD_DvzAerSnow = 0.03_fp,                      &
                              DD_F0         = 0.0_fp,                       &
                              DD_Hstar_old  = 0.0_fp,                       &
                              WD_AerScavEff = 1.0_fp,                       &
                              WD_KcScaleFac = KcScale,                      &
                              WD_RainoutEff = RainEff,                      &
                              RC            = RC )

          CASE( 'POPPBCPI', 'POPPOCPI' )

             ! These have identical properties except for the mol weights ...
             SELECT CASE( TRIM( Input_Opt%POP_TYPE ) )
                CASE( 'PHE' )
                   MW_g     = 178.23_fp
                   FullName = 'Phenanthrene particles on'
                CASE( 'PYR' )
                   MW_g     = 202.25_fp
                   FullName = 'Pyrene particles on'
                CASE( 'BaP' )
                   MW_g     = 252.31_fp
                   FullName = 'Benzo(a)pyrene particles on'
             END SELECT

             ! ... and the names
             SELECT CASE( TRIM( NameAllCaps ) ) 
                CASE( 'POPPBCPI' ) 
                   FullName = TRIM( FullName ) // ' hydrophilic black carbon'
                CASE( 'POPPOCPI' )
                   FullName = TRIM( FullName ) // ' hydrophilic organic carbon'
             END SELECT
             
             ! Halve the Kc (cloud condensate -> precip) rate
             ! for the temperature range 237 K <= T < 258 K.
             KcScale = (/ 1.0_fp, 0.5_fp, 1.0_fp /)

             ! Turn off rainout only when 237 K <= T < 258K.
             RainEff = (/ 1.0_fp, 0.0_fp, 1.0_fp /)

             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              Name          = NameAllCaps,                  &
                              FullName      = FullName,                     &
                              MW_g          = MW_g,                         &
                              Is_Advected   = Is_Advected,                  &
                              Is_Gas        = F,                            &
                              Is_Drydep     = T,                            &
                              Is_Wetdep     = T,                            &
                              DD_DvzAerSnow = 0.03_fp,                      &
                              DD_F0         = 0.0_fp,                       &
                              DD_Hstar_old  = 0.0_fp,                       &
                              WD_AerScavEff = 1.0_fp,                       &
                              WD_KcScaleFac = KcScale,                      &
                              WD_RainoutEff = RainEff,                      &
                              RC            = RC )

          !==================================================================
          ! Species for the CO2 specialty simulation
          !==================================================================

          CASE( 'CO2',    'CO2FF', 'CO2OC', 'CO2BAL', 'CO2BB', 'CO2BF',     &
                'CO2NTE', 'CO2SE', 'CO2AV', 'CO2CH',  'CO2CORR'         )

             ! These all have identical properties except for the names 
             ! Add TOMAS bin number to full name
             
             SELECT CASE( TRIM( NameAllCaps ) ) 
                CASE( 'CO2FF' )
                   Name     = 'CO2ff'
                   FullName = 'Carbon dioxide from fossil fuel emissions'
                CASE( 'CO2OC' )
                   Name     = 'CO2oc'
                   FullName = 'Carbon dioxide from ocean emissions'
                CASE( 'CO2BAL' )
                   Name     = 'CO2bal'
                   FullName = 'Carbon dioxide from balanced biosphere'
                CASE( 'CO2BB' )
                   Name     = 'CO2bb'
                   FullName = 'Carbon dioxide  from biomass burning emissions'
                CASE( 'CO2BF' )
                   Name     = 'CO2bf'
                   FullName = 'Carbon dioxide  from biofuel emissions'
                CASE( 'CO2NTE' )
                   Name     = 'CO2nte'
                   FullName = 'Carbon dioxide from net terrestrial exchange'
                CASE( 'CO2SE' )
                   Name     = 'CO2se'
                   FullName = 'Carbon dioxide from ship emissions'
                CASE( 'CO2AV' )
                   Name     = 'CO2av'
                   FullName = 'Carbon dioxide from aviation emissions'
                CASE( 'CO2CH' )
                   Name     = 'CO2ch'
                   FullName = 'Carbon dioxide from chemical sources'
                CASE( 'CO2CORR' )
                   Name     = 'CO2corr'
                   FullName = 'Carbon dioxide chemical source surface correction'
                CASE DEFAULT
                   Name     = 'CO2'
                   FullName = 'Carbon dioxide'
             END SELECT

             ! Set special default background for CO2
             SELECT CASE( TRIM( NameAllCaps ) )
                CASE( 'CO2' ) 
                   BackgroundVV = 3.55e-04_fp
                CASE DEFAULT
                   BackgroundVV = MISSING_VV
             END SELECT

             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              KppSpcId      = KppSpcId(N),                  &
                              KppVarId      = KppVarId(N),                  &
                              KppFixId      = KppFixId(N),                  &
                              Name          = Name,                         &
                              FullName      = FullName,                     &
                              MW_g          = 44.0_fp,                      &
                              BackgroundVV  = BackgroundVV,                 &
                              Is_Advected   = Is_Advected,                  &
                              Is_Gas        = T,                            &
                              Is_Drydep     = F,                            &
                              Is_Wetdep     = F,                            &
                              RC            = RC )

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
             FullName = 'Aerosol water, size bin ='
             C        = LEN_TRIM( NameAllCaps )
             IF ( C == 3 ) THEN
                FullName = TRIM( FullName ) // ' ' // NameAllCaps(C:C)
             ELSE
                FullName = TRIM( FullName ) // ' ' // NameAllCaps(C-1:C)
             ENDIF

             ! Halve the Kc (cloud condensate -> precip) rate
             ! for the temperature range 237 K <= T < 258 K.
             KcScale = (/ 1.0_fp, 0.5_fp, 1.0_fp /)

             ! Turn off rainout only when 237 K <= T < 258K.
             RainEff = (/ 1.0_fp, 0.0_fp, 1.0_fp /)

             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              Name          = NameAllCaps,                  &
                              FullName      = FullName,                     &
                              MW_g          = 18.0_fp,                      &
                              Is_Advected   = Is_Advected,                  &
                              Is_Gas        = F,                            &
                              Is_Drydep     = F,                            &
                              Is_Wetdep     = F,                            &
                              DD_DvzAerSnow = 0.03_fp,                      &
                              DD_F0         = 0.0_fp,                       &
                              DD_Hstar_old  = 0.0_fp,                       &
                              MP_SizeResAer = T,                            &
                              WD_AerScavEff = 1.0_fp,                       &
                              WD_KcScaleFac = KcScale,                      &
                              WD_RainoutEff = RainEff,                      &
                              RC            = RC )

          CASE( 'DUST1',  'DUST2',  'DUST3',  'DUST4',  'DUST5',            &
                'DUST6',  'DUST7',  'DUST8',  'DUST9',  'DUST10',           &
                'DUST11', 'DUST12', 'DUST13', 'DUST14', 'DUST15',           &
                'DUST16', 'DUST17', 'DUST18', 'DUST19', 'DUST20',           &
                'DUST21', 'DUST22', 'DUST23', 'DUST24', 'DUST25',           &
                'DUST26', 'DUST27', 'DUST28', 'DUST29', 'DUST30',           &
                'DUST31', 'DUST32', 'DUST33', 'DUST34', 'DUST35',           &
                'DUST36', 'DUST37', 'DUST38', 'DUST39', 'DUST40'  )
 
             ! Add TOMAS bin number to full name
             FullName = 'Mineral dust, size bin ='
             C        = LEN_TRIM( NameAllCaps )
             IF ( C == 5 ) THEN
                FullName = TRIM( FullName ) // ' ' // NameAllCaps(C:C)
             ELSE
                FullName = TRIM( FullName ) // ' ' // NameAllCaps(C-1:C)
             ENDIF

             ! Halve the Kc (cloud condensate -> precip) rate
             ! for the temperature range 237 K <= T < 258 K.
             KcScale = (/ 1.0_fp, 0.5_fp, 1.0_fp /)

             ! Turn off rainout only when 237 K <= T < 258K.
             RainEff = (/ 1.0_fp, 0.0_fp, 1.0_fp /)   

             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              Name          = NameAllCaps,                  &
                              FullName      = FullName,                     &
                              MW_g          = 100.0_fp,                     &
                              Is_Advected   = Is_Advected,                  &
                              Is_Gas        = F,                            &
                              Is_Drydep     = F,                            &
                              Is_Wetdep     = T,                            &
                              DD_DvzAerSnow = 0.03_fp,                      &
                              DD_F0         = 0.0_fp,                       &
                              DD_Hstar_old  = 0.0_fp,                       &
                              MP_SizeResAer = T,                            &
                              WD_AerScavEff = 1.0_fp,                       &
                              WD_KcScaleFac = KcScale,                      &
                              WD_RainoutEff = RainEff,                      &
                              RC            = RC )

          CASE( 'ECIL1',  'ECIL2',  'ECIL3',  'ECIL4',  'ECIL5',            &
                'ECIL6',  'ECIL7',  'ECIL8',  'ECIL9',  'ECIL10',           &
                'ECIL11', 'ECIL12', 'ECIL13', 'ECIL14', 'ECIL15',           &
                'ECIL16', 'ECIL17', 'ECIL18', 'ECIL19', 'ECIL20',           &
                'ECIL21', 'ECIL22', 'ECIL23', 'ECIL24', 'ECIL25',           &
                'ECIL26', 'ECIL27', 'ECIL28', 'ECIL29', 'ECIL30',           &
                'ECIL31', 'ECIL32', 'ECIL33', 'ECIL34', 'ECIL35',           &
                'ECIL36', 'ECIL37', 'ECIL38', 'ECIL39', 'ECIL40'  )
 
             ! Add TOMAS bin number to full name
             FullName = 'Hydrophilic elemental carbon, size bin ='
             C        = LEN_TRIM( NameAllCaps )
             IF ( C == 5 ) THEN 
                FullName = TRIM( FullName ) // ' ' // NameAllCaps(C:C)
             ELSE
                FullName = TRIM( FullName ) // ' ' // NameAllCaps(C-1:C)
             ENDIF

             ! Halve the Kc (cloud condensate -> precip) rate
             ! for the temperature range 237 K <= T < 258 K.
             KcScale = (/ 1.0_fp, 0.5_fp, 1.0_fp /)

             ! Turn off rainout only when 237 K <= T < 258K.
             RainEff = (/ 1.0_fp, 0.0_fp, 1.0_fp /)   

             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              Name          = NameAllCaps,                  &
                              FullName      = FullName,                     &
                              MW_g          = 12.0_fp,                      &
                              Is_Advected   = Is_Advected,                  &
                              Is_Gas        = F,                            &
                              Is_Drydep     = F,                            &
                              Is_Wetdep     = T,                            &
                              DD_DvzAerSnow = 0.03_fp,                      &
                              DD_F0         = 0.0_fp,                       &
                              DD_Hstar_old  = 0.0_fp,                       &
                              MP_SizeResAer = T,                            &
                              WD_AerScavEff = 1.0_fp,                       &
                              WD_KcScaleFac = KcScale,                      &
                              WD_RainoutEff = RainEff,                      &
                              RC            = RC )

          CASE( 'ECOB1',  'ECOB2',  'ECOB3',  'ECOB4',  'ECOB5',            &
                'ECOB6',  'ECOB7',  'ECOB8',  'ECOB9',  'ECOB10',           &
                'ECOB11', 'ECOB12', 'ECOB13', 'ECOB14', 'ECOB15',           &
                'ECOB16', 'ECOB17', 'ECOB18', 'ECOB19', 'ECOB20',           &
                'ECOB21', 'ECOB22', 'ECOB23', 'ECOB24', 'ECOB25',           &
                'ECOB26', 'ECOB27', 'ECOB28', 'ECOB29', 'ECOB30',           &
                'ECOB31', 'ECOB32', 'ECOB33', 'ECOB34', 'ECOB35',           &
                'ECOB36', 'ECOB37', 'ECOB38', 'ECOB39', 'ECOB40'  )
 
             ! Add TOMAS bin number to full name
             FullName = 'Hydrophobic elemental carbon, size bin ='
             C        = LEN_TRIM( NameAllCaps )
             IF ( C == 5 ) THEN 
                FullName = TRIM( FullName ) // ' ' // NameAllCaps(C:C)
             ELSE
                FullName = TRIM( FullName ) // ' ' // NameAllCaps(C-1:C)
             ENDIF

             ! Halve the Kc (cloud condensate -> precip) rate
             ! for the temperature range 237 K <= T < 258 K.
             KcScale = (/ 1.0_fp, 0.5_fp, 1.0_fp /)

             ! Turn off rainout only when 237 K <= T < 258K.
             RainEff = (/ 1.0_fp, 0.0_fp, 1.0_fp /)   

             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              Name          = NameAllCaps,                  &
                              FullName      = FullName,                     &
                              MW_g          = 12.0_fp,                      &
                              Is_Advected   = Is_Advected,                  &
                              Is_Gas        = F,                            &
                              Is_Drydep     = F,                            &
                              Is_Wetdep     = F,                            &
                              DD_DvzAerSnow = 0.03_fp,                      &
                              DD_F0         = 0.0_fp,                       &
                              DD_Hstar_old  = 0.0_fp,                       &
                              MP_SizeResAer = T,                            &
                              WD_AerScavEff = 1.0_fp,                       &
                              WD_KcScaleFac = KcScale,                      &
                              WD_RainoutEff = RainEff,                      &
                              RC            = RC )

          CASE( 'H2SO4' )

             !%%% NOTE: The TOMAS H2SO4 species dry-deposits like a gas, 
             !%%% wet-deposits as an aerosol.  So we need to give this
             !%%% both gas and aerosol properties (ewl, bmy, 10/13/15)

             ! Halve the Kc (cloud condensate -> precip) rate
             ! for the temperature range 237 K <= T < 258 K.
             KcScale = (/ 1.0_fp, 0.5_fp, 1.0_fp /)

             ! Turn off rainout only when 237 K <= T < 258K.
             RainEff = (/ 1.0_fp, 0.0_fp, 1.0_fp /)

             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              Name          = NameAllCaps,                  &
                              FullName      = 'Sulfuric acid',              &
                              MW_g          = 98.0_fp,                      &
                              MolecRatio    = 1.0_fp,                       &
                              Is_Advected   = Is_Advected,                  &
                              Is_Gas        = T,                            &
                              Is_Drydep     = T,                            &
                              Is_Wetdep     = T,                            &
                              DD_F0         = 0.0_fp,                       &
                              DD_Hstar_old  = 1.00e+5_fp,                   &
                              MP_SizeResAer = T,                            &
                              WD_AerScavEff = 1.0_fp,                       &
                              WD_KcScaleFac = KcScale,                      &
                              WD_RainoutEff = RainEff,                      &
                              RC            = RC )

          CASE( 'NK1',  'NK2',  'NK3',  'NK4',  'NK5',  'NK6',  'NK7',      &
                'NK8',  'NK9',  'NK10', 'NK11', 'NK12', 'NK13', 'NK14',     &
                'NK15', 'NK16', 'NK17', 'NK18', 'NK19', 'NK20', 'NK21',     &
                'NK22', 'NK23', 'NK24', 'NK25', 'NK26', 'NK27', 'NK28',     &
                'NK29', 'NK30', 'NK31', 'NK32', 'NK33', 'NK34', 'NK35',     &
                'NK36', 'NK37', 'NK38', 'NK39', 'NK40'                  )

             ! Add TOMAS bin number to full name
             FullName = 'Aerosol number, size bin ='
             C        = LEN_TRIM( NameAllCaps )
             IF ( C == 3 ) THEN 
                FullName = TRIM( FullName ) // ' ' // NameAllCaps(C:C)
             ELSE
                FullName = TRIM( FullName ) // ' ' // NameAllCaps(C-1:C)
             ENDIF

             ! Halve the Kc (cloud condensate -> precip) rate
             ! for the temperature range 237 K <= T < 258 K.
             KcScale = (/ 1.0_fp, 0.5_fp, 1.0_fp /)

             ! Turn off rainout only when 237 K <= T < 258K.
             RainEff = (/ 1.0_fp, 0.0_fp, 1.0_fp /)   

             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              KppSpcId      = KppSpcId(N),                  &
                              KppVarId      = KppVarId(N),                  &
                              KppFixId      = KppFixId(N),                  &
                              Name          = NameAllCaps,                  &
                              FullName      = FullName,                     &
                              MW_g          = 1.0_fp,                       &
                              Is_Advected   = Is_Advected,                  &
                              Is_Gas        = F,                            &
                              Is_Drydep     = T,                            &
                              Is_Wetdep     = T,                            &
                              DD_DvzAerSnow = 0.03_fp,                      &
                              DD_DustDryDep = T,                            &
                              DD_F0         = 0.0_fp,                       &
                              DD_Hstar_old  = 0.0_fp,                       &
                              MP_SizeResNum = T,                            &
                              WD_AerScavEff = 1.0_fp,                       &
                              WD_KcScaleFac = KcScale,                      &
                              WD_RainoutEff = RainEff,                      &
                              RC            = RC )

          CASE( 'OCIL1',  'OCIL2',  'OCIL3',  'OCIL4',  'OCIL5',            &
                'OCIL6',  'OCIL7',  'OCIL8',  'OCIL9',  'OCIL10',           &
                'OCIL11', 'OCIL12', 'OCIL13', 'OCIL14', 'OCIL15',           &
                'OCIL16', 'OCIL17', 'OCIL18', 'OCIL19', 'OCIL20',           &
                'OCIL21', 'OCIL22', 'OCIL23', 'OCIL24', 'OCIL25',           &
                'OCIL26', 'OCIL27', 'OCIL28', 'OCIL29', 'OCIL30',           &
                'OCIL31', 'OCIL32', 'OCIL33', 'OCIL34', 'OCIL35',           &
                'OCIL36', 'OCIL37', 'OCIL38', 'OCIL39', 'OCIL40'  )
 
             ! Add TOMAS bin number to full name
             FullName = 'Hydrophilic organic carbon, size bin ='
             C        = LEN_TRIM( NameAllCaps )
             IF ( C == 5 ) THEN 
                FullName = TRIM( FullName ) // ' ' // NameAllCaps(C:C)
             ELSE
                FullName = TRIM( FullName ) // ' ' // NameAllCaps(C-1:C)
             ENDIF

             ! Halve the Kc (cloud condensate -> precip) rate
             ! for the temperature range 237 K <= T < 258 K.
             KcScale = (/ 1.0_fp, 0.5_fp, 1.0_fp /)

             ! Turn off rainout only when 237 K <= T < 258K.
             RainEff = (/ 1.0_fp, 0.0_fp, 1.0_fp /)   

             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              Name          = NameAllCaps,                  &
                              FullName      = FullName,                     &
                              MW_g          = 12.0_fp,                      &
                              Is_Advected   = Is_Advected,                  &
                              Is_Gas        = F,                            &
                              Is_Drydep     = F,                            &
                              Is_Wetdep     = T,                            &
                              DD_DvzAerSnow = 0.03_fp,                      &
                              DD_F0         = 0.0_fp,                       &
                              DD_Hstar_old  = 0.0_fp,                       &
                              MP_SizeResAer = T,                            &
                              WD_AerScavEff = 1.0_fp,                       &
                              WD_KcScaleFac = KcScale,                      &
                              WD_RainoutEff = RainEff,                      &
                              RC            = RC )

          CASE( 'OCOB1',  'OCOB2',  'OCOB3',  'OCOB4',  'OCOB5',            &
                'OCOB6',  'OCOB7',  'OCOB8',  'OCOB9',  'OCOB10',           &
                'OCOB11', 'OCOB12', 'OCOB13', 'OCOB14', 'OCOB15',           &
                'OCOB16', 'OCOB17', 'OCOB18', 'OCOB19', 'OCOB20',           &
                'OCOB21', 'OCOB22', 'OCOB23', 'OCOB24', 'OCOB25',           &
                'OCOB26', 'OCOB27', 'OCOB28', 'OCOB29', 'OCOB30',           &
                'OCOB31', 'OCOB32', 'OCOB33', 'OCOB34', 'OCOB35',           &
                'OCOB36', 'OCOB37', 'OCOB38', 'OCOB39', 'OCOB40'  )
 
             ! Add TOMAS bin number to full name
             FullName = 'Hydrophobic organic carbon, size bin ='
             C        = LEN_TRIM( NameAllCaps )
             IF ( C == 5 ) THEN 
                FullName = TRIM( FullName ) // ' ' // NameAllCaps(C:C)
             ELSE
                FullName = TRIM( FullName ) // ' ' // NameAllCaps(C-1:C)
             ENDIF

             ! Halve the Kc (cloud condensate -> precip) rate
             ! for the temperature range 237 K <= T < 258 K.
             KcScale = (/ 1.0_fp, 0.5_fp, 1.0_fp /)

             ! Turn off rainout only when 237 K <= T < 258K.
             RainEff = (/ 1.0_fp, 0.0_fp, 1.0_fp /)   

             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              Name          = NameAllCaps,                  &
                              FullName      = FullName,                     &
                              MW_g          = 12.0_fp,                      &
                              Is_Advected   = Is_Advected,                  &
                              Is_Gas        = F,                            &
                              Is_Drydep     = F,                            &
                              Is_Wetdep     = T,                            &
                              DD_DvzAerSnow = 0.03_fp,                      &
                              DD_F0         = 0.0_fp,                       &
                              DD_Hstar_old  = 0.0_fp,                       &
                              MP_SizeResAer = T,                            &
                              WD_AerScavEff = 1.0_fp,                       &
                              WD_KcScaleFac = KcScale,                      &
                              WD_RainoutEff = RainEff,                      &
                              RC            = RC )

          CASE( 'SF1',  'SF2',  'SF3',  'SF4',  'SF5',  'SF6',  'SF7',      &
                'SF8',  'SF9',  'SF10', 'SF11', 'SF12', 'SF13', 'SF14',     &
                'SF15', 'SF16', 'SF17', 'SF18', 'SF19', 'SF20', 'SF21',     &
                'SF22', 'SF23', 'SF24', 'SF25', 'SF26', 'SF27', 'SF28',     &
                'SF29', 'SF30', 'SF31', 'SF32', 'SF33', 'SF34', 'SF35',     &
                'SF36', 'SF37', 'SF38', 'SF39', 'SF40'                  )
 
             ! Add TOMAS bin number to full name
             FullName = 'Sulfate aerosol, size bin ='
             C        = LEN_TRIM( NameAllCaps )
             IF ( C == 3 ) THEN 
                FullName = TRIM( FullName ) // ' ' // NameAllCaps(C:C)
             ELSE
                FullName = TRIM( FullName ) // ' ' // NameAllCaps(C-1:C)
             ENDIF

             ! Halve the Kc (cloud condensate -> precip) rate
             ! for the temperature range 237 K <= T < 258 K.
             KcScale = (/ 1.0_fp, 0.5_fp, 1.0_fp /)

             ! Turn off rainout only when 237 K <= T < 258K.
             RainEff = (/ 1.0_fp, 0.0_fp, 1.0_fp /)   

             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              Name          = NameAllCaps,                  &
                              FullName      = FullName,                     &
                              MW_g          = 96.0_fp,                      &
                              Is_Advected   = Is_Advected,                  &
                              Is_Gas        = F,                            &
                              Is_Drydep     = F,                            &
                              Is_Wetdep     = T,                            &
                              DD_DvzAerSnow = 0.03_fp,                      &
                              DD_F0         = 0.0_fp,                       &
                              DD_Hstar_old  = 0.0_fp,                       &
                              MP_SizeResAer = T,                            &
                              WD_AerScavEff = 1.0_fp,                       &
                              WD_KcScaleFac = KcScale,                      &
                              WD_RainoutEff = RainEff,                      &
                              RC            = RC )

          CASE( 'SS1',  'SS2',  'SS3',  'SS4',  'SS5',  'SS6',  'SS7',      &
                'SS8',  'SS9',  'SS10', 'SS11', 'SS12', 'SS13', 'SS14',     &
                'SS15', 'SS16', 'SS17', 'SS18', 'SS19', 'SS20', 'SS21',     &
                'SS22', 'SS23', 'SS24', 'SS25', 'SS26', 'SS27', 'SS28',     &
                'SS29', 'SS30', 'SS31', 'SS32', 'SS33', 'SS34', 'SS35',     &
                'SS36', 'SS37', 'SS38', 'SS39', 'SS40'                  )
 
             ! Add TOMAS bin number to full name
             FullName = 'Sea salt aerosol, size bin = '
             C        = LEN_TRIM( NameAllCaps )
             IF ( C == 3 ) THEN 
                FullName = TRIM( FullName ) // ' ' // NameAllCaps(C:C)
             ELSE
                FullName = TRIM( FullName ) // ' ' // NameAllCaps(C-1:C)
             ENDIF

             ! Halve the Kc (cloud condensate -> precip) rate
             ! for the temperature range 237 K <= T < 258 K.
             KcScale = (/ 1.0_fp, 0.5_fp, 1.0_fp /)

             ! Turn off rainout only when 237 K <= T < 258K.
             RainEff = (/ 1.0_fp, 0.0_fp, 1.0_fp /)   

             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              Name          = NameAllCaps,                  &
                              FullName      = FullName,                     &
                              MW_g          = 58.5_fp,                      &
                              Is_Advected   = Is_Advected,                  &
                              Is_Gas        = F,                            &
                              Is_Drydep     = F,                            &
                              Is_Wetdep     = T,                            &
                              DD_DvzAerSnow = 0.03_fp,                      &
                              DD_F0         = 0.0_fp,                       &
                              DD_Hstar_old  = 0.0_fp,                       &
                              MP_SizeResAer = T,                            &
                              WD_AerScavEff = 1.0_fp,                       &
                              WD_KcScaleFac = KcScale,                      &
                              WD_RainoutEff = RainEff,                      &
                              RC            = RC )

#endif

          !==================================================================
          ! Additional species for FAST-JX photolysis
          !==================================================================

          CASE( 'O2' )
             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              BackgroundVV  = 2.095e-01_fp,                 &
                              KppSpcId      = KppSpcId(N),                  &
                              KppVarId      = KppVarId(N),                  &
                              KppFixId      = KppFixId(N),                  &
                              Name          = NameAllCaps,                  &
                              Is_Advected   = F,                            &
                              Is_Gas        = T,                            &
                              Is_Drydep     = F,                            &
                              Is_Wetdep     = F,                            &
                              Is_Photolysis = T,                            &
                              RC            = RC )

          CASE( 'GLYX' )
             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              KppSpcId      = KppSpcId(N),                  &
                              KppVarId      = KppVarId(N),                  &
                              KppFixId      = KppFixId(N),                  &
                              Name          = NameAllCaps,                  &
                              Is_Advected   = F,                            &
                              Is_Gas        = T,                            &
                              Is_Drydep     = F,                            &
                              Is_Wetdep     = F,                            &
                              Is_Photolysis = T,                            &
                              RC            = RC )

          CASE( 'MGLY' )
             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              KppSpcId      = KppSpcId(N),                  &
                              KppVarId      = KppVarId(N),                  &
                              KppFixId      = KppFixId(N),                  &
                              Name          = NameAllCaps,                  &
                              Is_Advected   = F,                            &
                              Is_Gas        = T,                            &
                              Is_Drydep     = F,                            &
                              Is_Wetdep     = F,                            &
                              Is_Photolysis = T,                            &
                              RC            = RC )

          CASE( 'INPN' )
             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              KppSpcId      = KppSpcId(N),                  &
                              KppVarId      = KppVarId(N),                  &
                              KppFixId      = KppFixId(N),                  &
                              Name          = NameAllCaps,                  &
                              Is_Advected   = F,                            &
                              Is_Gas        = T,                            &
                              Is_Drydep     = F,                            &
                              Is_Wetdep     = F,                            &
                              Is_Photolysis = T,                            &
                              RC            = RC )

          CASE( 'PRPN' )
             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              KppSpcId      = KppSpcId(N),                  &
                              KppVarId      = KppVarId(N),                  &
                              KppFixId      = KppFixId(N),                  &
                              Name          = NameAllCaps,                  &
                              Is_Advected   = F,                            &
                              Is_Gas        = T,                            &
                              Is_Drydep     = F,                            &
                              Is_Wetdep     = F,                            &
                              Is_Photolysis = T,                            &
                              RC            = RC )

          CASE( 'ETP' )
             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              KppSpcId      = KppSpcId(N),                  &
                              KppVarId      = KppVarId(N),                  &
                              KppFixId      = KppFixId(N),                  &
                              Name          = NameAllCaps,                  &
                              Is_Advected   = F,                            &
                              Is_Gas        = T,                            &
                              Is_Drydep     = F,                            &
                              Is_Wetdep     = F,                            &
                              Is_Photolysis = T,                            &
                              RC            = RC )

          CASE( 'RA3P' )
             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              KppSpcId      = KppSpcId(N),                  &
                              KppVarId      = KppVarId(N),                  &
                              KppFixId      = KppFixId(N),                  &
                              Name          = NameAllCaps,                  &
                              Is_Advected   = F,                            &
                              Is_Gas        = T,                            &
                              Is_Drydep     = F,                            &
                              Is_Wetdep     = F,                            &
                              Is_Photolysis = T,                            &
                              RC            = RC )

          CASE( 'RB3P' )
             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              KppSpcId      = KppSpcId(N),                  &
                              KppVarId      = KppVarId(N),                  &
                              KppFixId      = KppFixId(N),                  &
                              Name          = NameAllCaps,                  &
                              Is_Advected   = F,                            &
                              Is_Gas        = T,                            &
                              Is_Drydep     = F,                            &
                              Is_Wetdep     = F,                            &
                              Is_Photolysis = T,                            &
                              RC            = RC )

          CASE( 'R4P' )
             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              KppSpcId      = KppSpcId(N),                  &
                              KppVarId      = KppVarId(N),                  &
                              KppFixId      = KppFixId(N),                  &
                              Name          = NameAllCaps,                  &
                              Is_Advected   = F,                            &
                              Is_Gas        = T,                            &
                              Is_Drydep     = F,                            &
                              Is_Wetdep     = F,                            &
                              Is_Photolysis = T,                            &
                              RC            = RC )

          CASE( 'PP' )
             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              KppSpcId      = KppSpcId(N),                  &
                              KppVarId      = KppVarId(N),                  &
                              KppFixId      = KppFixId(N),                  &
                              Name          = NameAllCaps,                  &
                              Is_Advected   = F,                            &
                              Is_Gas        = T,                            &
                              Is_Drydep     = F,                            &
                              Is_Wetdep     = F,                            &
                              Is_Photolysis = T,                            &
                              RC            = RC )

          CASE( 'RP' )
             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              KppSpcId      = KppSpcId(N),                  &
                              KppVarId      = KppVarId(N),                  &
                              KppFixId      = KppFixId(N),                  &
                              Name          = NameAllCaps,                  &
                              Is_Advected   = F,                            &
                              Is_Gas        = T,                            &
                              Is_Drydep     = F,                            &
                              Is_Wetdep     = F,                            &
                              Is_Photolysis = T,                            &
                              RC            = RC )

          CASE( 'IAP' )
             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              KppSpcId      = KppSpcId(N),                  &
                              KppVarId      = KppVarId(N),                  &
                              KppFixId      = KppFixId(N),                  &
                              Name          = NameAllCaps,                  &
                              Is_Advected   = F,                            &
                              Is_Gas        = T,                            &
                              Is_Drydep     = F,                            &
                              Is_Wetdep     = F,                            &
                              Is_Photolysis = T,                            &
                              RC            = RC )

          CASE( 'ISNP' )
             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              KppSpcId      = KppSpcId(N),                  &
                              KppVarId      = KppVarId(N),                  &
                              KppFixId      = KppFixId(N),                  &
                              Name          = NameAllCaps,                  &
                              Is_Advected   = F,                            &
                              Is_Gas        = T,                            &
                              Is_Drydep     = F,                            &
                              Is_Wetdep     = F,                            &
                              Is_Photolysis = T,                            &
                              RC            = RC )

          CASE( 'VRP' )
             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              KppSpcId      = KppSpcId(N),                  &
                              KppVarId      = KppVarId(N),                  &
                              KppFixId      = KppFixId(N),                  &
                              Name          = NameAllCaps,                  &
                              Is_Advected   = F,                            &
                              Is_Gas        = T,                            &
                              Is_Drydep     = F,                            &
                              Is_Wetdep     = F,                            &
                              Is_Photolysis = T,                            &
                              RC            = RC )

          CASE( 'MRP' )
             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              KppSpcId      = KppSpcId(N),                  &
                              KppVarId      = KppVarId(N),                  &
                              KppFixId      = KppFixId(N),                  &
                              Name          = NameAllCaps,                  &
                              Is_Advected   = F,                            &
                              Is_Gas        = T,                            &
                              Is_Drydep     = F,                            &
                              Is_Wetdep     = F,                            &
                              Is_Photolysis = T,                            &
                              RC            = RC )

          CASE( 'MAOP' )
             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              KppSpcId      = KppSpcId(N),                  &
                              KppVarId      = KppVarId(N),                  &
                              KppFixId      = KppFixId(N),                  &
                              Name          = NameAllCaps,                  &
                              Is_Advected   = F,                            &
                              Is_Gas        = T,                            &
                              Is_Drydep     = F,                            &
                              Is_Wetdep     = F,                            &
                              Is_Photolysis = T,                            &
                              RC            = RC )

          CASE( 'ATOOH' )
             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              KppSpcId      = KppSpcId(N),                  &
                              KppVarId      = KppVarId(N),                  &
                              KppFixId      = KppFixId(N),                  &
                              Name          = NameAllCaps,                  &
                              Is_Advected   = F,                            &
                              Is_Gas        = T,                            &
                              Is_Drydep     = F,                            &
                              Is_Wetdep     = F,                            &
                              Is_Photolysis = T,                            &
                              RC            = RC )

          CASE( 'HO2' )
                CALL Spc_Create( am_I_Root     = am_I_Root,                 &
                                 ThisSpc       = SpcData(N)%Info,           &
                                 ModelID       = N,                         &
                                 KppSpcId      = KppSpcId(N),               &
                                 KppVarId      = KppVarId(N),               &
                                 KppFixId      = KppFixId(N),               &
                                 Name          = NameAllCaps,               &
                                 BackgroundVV  = 4.0e-15_fp,                &
                                 Is_Advected   = F,                         &
                                 Is_Gas        = T,                         &
                                 Is_Drydep     = F,                         &
                                 Is_Wetdep     = F,                         &
                                 RC            = RC )

          CASE( 'MO2' )
                CALL Spc_Create( am_I_Root     = am_I_Root,                 &
                                 ThisSpc       = SpcData(N)%Info,           &
                                 ModelID       = N,                         &
                                 KppSpcId      = KppSpcId(N),               &
                                 KppVarId      = KppVarId(N),               &
                                 KppFixId      = KppFixId(N),               &
                                 Name          = NameAllCaps,               &
                                 BackgroundVV  = 4.0e-15_fp,                &
                                 Is_Advected   = F,                         &
                                 Is_Gas        = T,                         &
                                 Is_Drydep     = F,                         &
                                 Is_Wetdep     = F,                         &
                                 RC            = RC )

          CASE( 'OH' )
                CALL Spc_Create( am_I_Root     = am_I_Root,                 &
                                 ThisSpc       = SpcData(N)%Info,           &
                                 ModelID       = N,                         &
                                 KppSpcId      = KppSpcId(N),               &
                                 KppVarId      = KppVarId(N),               &
                                 KppFixId      = KppFixId(N),               &
                                 Name          = NameAllCaps,               &
                                 BackgroundVV  = 4.0e-15_fp,                &
                                 Is_Advected   = F,                         &
                                 Is_Gas        = T,                         &
                                 Is_Drydep     = F,                         &
                                 Is_Wetdep     = F,                         &
                                 RC            = RC )

          CASE( 'H2' )
                CALL Spc_Create( am_I_Root     = am_I_Root,                 &
                                 ThisSpc       = SpcData(N)%Info,           &
                                 ModelID       = N,                         &
                                 KppSpcId      = KppSpcId(N),               &
                                 KppVarId      = KppVarId(N),               &
                                 KppFixId      = KppFixId(N),               &
                                 Name          = NameAllCaps,               &
                                 BackgroundVV  = 5.0e-07_fp,                &
                                 Is_Advected   = F,                         &
                                 Is_Gas        = T,                         &
                                 Is_Drydep     = F,                         &
                                 Is_Wetdep     = F,                         &
                                 RC            = RC )

          CASE( 'N' )
                CALL Spc_Create( am_I_Root     = am_I_Root,                 &
                                 ThisSpc       = SpcData(N)%Info,           &
                                 ModelID       = N,                         &
                                 KppSpcId      = KppSpcId(N),               &
                                 KppVarId      = KppVarId(N),               &
                                 KppFixId      = KppFixId(N),               &
                                 Name          = NameAllCaps,               &
                                 BackgroundVV  = 4.0e-20_fp,                &
                                 Is_Advected   = F,                         &
                                 Is_Gas        = T,                         &
                                 Is_Drydep     = F,                         &
                                 Is_Wetdep     = F,                         &
                                 RC            = RC )

          CASE( 'O1D' )
                CALL Spc_Create( am_I_Root     = am_I_Root,                 &
                                 ThisSpc       = SpcData(N)%Info,           &
                                 ModelID       = N,                         &
                                 KppSpcId      = KppSpcId(N),               &
                                 KppVarId      = KppVarId(N),               &
                                 KppFixId      = KppFixId(N),               &
                                 Name          = NameAllCaps,               &
#if defined( UCX )
                                 BackgroundVV  = 1.0e-15_fp,                &
#else
                                 BackgroundVV  = 4.0e-22_fp,                &
#endif
                                 Is_Advected   = F,                         &
                                 Is_Gas        = T,                         &
                                 Is_Drydep     = F,                         &
                                 Is_Wetdep     = F,                         &
                                 RC            = RC )

          !==================================================================
          ! Special handling for species not found in the list above
          !==================================================================
          CASE DEFAULT
  
             ! Test if this is a passive tracer
             CALL PASSIVE_TRACER_INQUIRE( NameAllCaps,          &
                                          IsPassive=IsPassive,  &
                                          MW=MW_g,              &
                                          InitConc=BackgroundVV  )

             ! Add as passive tracer if it is listed in passive_tracer_mod.F90.
             ! Passive tracers can be listed in the optional passive tracer menu
             ! in input_geos. When reading the input file, all listed passive 
             ! tracer quantities are written to local variables within module
             ! passive_tracer_mod.F90. Pass these values here to the species
             ! database (ckeller, 11/3/16).
             IF ( IsPassive ) THEN

                CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                                 ThisSpc       = SpcData(N)%Info,              &
                                 ModelID       = N,                            &
                                 Name          = NameAllCaps,                  &
                                 MW_g          = MW_g,                         &
                                 EmMW_g        = MW_g,                         &
                                 BackgroundVV  = BackgroundVV,                 &
                                 MolecRatio    = 1.0_fp,                       &
                                 Is_Advected   = T,                            &
                                 Is_Gas        = F,                            &
                                 Is_Drydep     = F,                            &
                                 Is_Wetdep     = F,                            &
                                 Is_Photolysis = F,                            &
                                 RC            = RC )
   
             ! Test if this is a non-advected chemical species
             ELSEIF ( KppSpcId(N) > 0 ) THEN 
                
                !------------------------------------------------------------
                ! If this is a non-advected KPP chemical species, then just
                ! create a basic default entry in the species database
                !------------------------------------------------------------
                CALL Spc_Create( am_I_Root     = am_I_Root,                 &
                                 ThisSpc       = SpcData(N)%Info,           &
                                 ModelID       = N,                         &
                                 KppSpcId      = KppSpcId(N),               &
                                 KppVarId      = KppVarId(N),               &
                                 KppFixId      = KppFixId(N),               &
                                 Name          = NameAllCaps,               &
                                 Is_Advected   = F,                         &
                                 Is_Gas        = T,                         &
                                 Is_Drydep     = F,                         &
                                 Is_Wetdep     = F,                         &
                                 RC            = RC )

             ELSE
                   
                !------------------------------------------------------------
                ! If this species i not found, the exit with error!
                ! create a default entry in the species database
                !------------------------------------------------------------
                WRITE( 6, '(a)' ) REPEAT( '=', 79 )
                WRITE( 6, 100 ) TRIM( NameAllCaps )
                WRITE( 6, 110 )
                WRITE( 6, 120 )
100             FORMAT( 'Species ', a, ' not found in the database!'       )
110             FORMAT( 'Please add the information to the CASE statement' )
120             FORMAT( 'in module Headers/species_database_mod.F90!'      )
                RC = -1
                RETURN
             ENDIF

       END SELECT

       ! Error
       IF ( RC /= GC_SUCCESS ) THEN
          PRINT*, '### Could not initialize species vector!'
          CALL EXIT( -999 ) 
       ENDIF

       ! Print info about each species
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
! !DESCRIPTION: Stores the list of unique species names (i.e. removing
!  duplicates from the list of advected species and the the list of KPP
!  species) for later use.  Also computes the indices for KPP variable
!  and fixed indices.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Cleanup_Work_Arrays()
!
! !REMARKS:
!  This may not be the fastest search algorithm, but it is only executed
!  once, at startup.
!
! !REVISION HISTORY:
!  06 May 2016 - R. Yantosca - Initial version
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

  END SUBROUTINE Cleanup_Work_ArrayS
!EOC
END MODULE Species_Database_Mod
