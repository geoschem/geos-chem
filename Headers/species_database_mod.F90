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
! !IROUTINE: init_species_database
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
    USE GIGC_ErrCode_Mod
    USE GIGC_Input_Opt_Mod, ONLY : OptInput
    USE Species_Mod
!
! !INPUT PARAMETERS: 
!
    LOGICAL,        INTENT(IN)  :: am_I_Root    ! Are we on the root CPU?
    TYPE(OptInput), INTENT(IN)  :: Input_Opt    ! Input Options object
!
! !INPUT/OUTPUT PARAMETERS: 
!
    TYPE(SpcPtr),   POINTER     :: SpcData(:)    ! Vector with species info
!
! !OUTPUT PARAMETERS: 
!
    INTEGER,        INTENT(OUT) :: RC           ! Success or failure?
!
! !REMARKS:
!  References for the new Henry's law constants:
!    (1) Sander et al [2015]: http://henrys-law.org
!    (2) http://wiki.geos-chem.org/Physical_properties_of_GEOS-Chem_species
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
!  ALSO NOTE: For now we are using this to define the tracers but we will
!  eventually extend this to the number of species. (bmy, 7/22/15)
!
! !REVISION HISTORY:
!  22 Jul 2015 - R. Yantosca - Initial version
!  01 Sep 2015 - R. Yantosca - Add Henry K0, CR constants for DMS, ACET
!  02 Sep 2015 - R. Yantosca - Corrected typos for some SOA species
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER             :: C,        N,   nSpecies
    REAL(fp)            :: A_Radius, KOA, MW_g   
    REAL(fp)            :: K0,       CR,  HStar

    ! Arrays
    REAL(fp)            :: RainEff(3)

    ! Strings
    CHARACTER(LEN=31)   :: NameAllCaps
    CHARACTER(LEN=31)   :: Name
    CHARACTER(LEN=80)   :: FullName

    ! Some species use the same drydep velocities as others, etc.
    INTEGER             :: DryDepID_PAN
    INTEGER             :: DryDepID_HNO3
!
! !DEFINED PARAMETERS
!
    LOGICAL,  PARAMETER :: T        = .TRUE.         ! Yes
    LOGICAL,  PARAMETER :: F        = .FALSE.        ! No
    REAL(f8), PARAMETER :: To_M_atm = 9.86923e-3_f8  ! mol/m3/Pa -> M/atm

    !=======================================================================
    ! Init_Species_Database begins here!
    !=======================================================================

    ! Number of species
    nSpecies = Input_Opt%N_TRACERS

    ! Initialize the species vector
    CALL SpcData_Init( am_I_Root, nSpecies, SpcData, RC )
    IF ( RC /= GIGC_SUCCESS ) THEN
       PRINT*, '### Could not initialize species vector!'
       CALL EXIT( -999 )
    ENDIF
    
    ! Loop over tracers
    DO N = 1, nSpecies

       ! Translate species name to uppercase
       NameAllCaps = TRIM( Input_Opt%TRACER_NAME(N) )
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
          !
          ! NOTE: Rainout efficiency is a 3-element vector:
          ! Element 1: Efficiency when T < 237 K.
          ! Element 2: Efficiency when 237 K <= T < 258 K
          ! Element 3: Efficiency when T > 258 K.
          !==================================================================

          CASE( 'ACET' )
             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              Name          = NameAllCaps,                  &
                              FullName      = 'Acetone',                    &
                              MW_g          = 58.08_fp,                     &
                              EmMW_g        = 12.00_fp,                     &
                              MolecRatio    = 3.0_fp,                       &
                              Is_Advected   = T,                            &
                              Is_Gas        = T,                            &
                              Is_Drydep     = T,                            &
                              Is_Wetdep     = F,                            &
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
                              Name          = NameAllCaps,                  &
                              FullName      = 'Acetaldehyde',               &
                              MW_g          = 44.05_fp,                     &
                              EmMW_g        = 12.0_fp,                      &
                              MolecRatio    = 2.0_fp,                       &
                              Is_Advected   = T,                            &
                              Is_Gas        = T,                            &
                              Is_Drydep     = T,                            &
                              Is_Wetdep     = F,                            &
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
                              Name          = NameAllCaps,                  &
                              FullName      = 'Lumped >= C4 Alkanes',       &
                              MW_g          = 58.12_fp,                     &
                              EmMW_g        = 12.0_fp,                      &
                              MolecRatio    = 4.0_fp,                       &
                              Is_Advected   = T,                            &
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
             
             ! Rainout efficiency is 0.8 because this is an SOA tracer.
             ! Turn off rainout when 237 K <= T < 258 K.
             RainEff = (/ 0.8_fp, 0.0_fp, 0.8_fp /)

             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              Name          = NameAllCaps,                  &
                              FullName      = FullName,                     &
                              MW_g          = 150.0_fp,                     &
                              Is_Advected   = T,                            &
                              Is_Gas        = F,                            &
                              Is_Drydep     = T,                            &
                              Is_Wetdep     = T,                            &
                              DD_F0         = 0.0_fp,                       &
                              DD_Hstar_old  = 0.0_fp,                       &
                              WD_AerScavEff = 0.8_fp,                       &
                              WD_RainoutEff = RainEff,                      &
                              RC            = RC )

          CASE( 'ASOG1', 'ASOG2', 'ASOG3' )
             FullName = 'Lumped non-volatile gas products of light aromatics + IVOCs'
             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              Name          = NameAllCaps,                  &
                              FullName      = FullName,                     &
                              MW_g          = 150.0_fp,                     &
                              Is_Advected   = T,                            &
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
             
             ! These have identical properties except for 
             ! the names and rainout efficiencies
             !
             ! NOTE: Hydrophobic BC, which normally has a rainout efficiency
             ! of zero, is considered to be IN (ice nuclei) and therefore is
             ! allowed to be scavenged at temperatures below 258 K.
             SELECT CASE( NameAllCaps )
                CASE( 'BCPI' )
                   FullName = 'Hydrophilic black carbon aerosol'
                   RainEff  = (/ 1.0_fp, 0.0_fp, 1.0_fp /)
                CASE( 'BCPO' )
                   Fullname = 'Hydrophobic black carbon aerosol'
                   RainEff  = (/ 1.0_fp, 1.0_fp, 0.0_fp /)
             END SELECT

             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              Name          = NameAllCaps,                  &
                              FullName      = FullName,                     &
                              MW_g          = 12.01_fp,                     &
                              EmMW_g        = 12.0_fp,                      &
                              Is_Advected   = T,                            &
                              Is_Gas        = F,                            &
                              Is_Drydep     = T,                            &
                              Is_Wetdep     = T,                            &
                              DD_F0         = 0.0_fp,                       &
                              DD_Hstar_Old  = 0.0_fp,                       &
                              WD_AerScavEff = 1.0_fp,                       &
                              WD_RainoutEff = RainEff,                      &
                              RC            = RC )

          CASE( 'BENZ' )
             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              Name          = NameAllCaps,                  &
                              FullName      = 'Benzene',                    &
                              MW_g          = 78.11_fp,                     &
                              EmMW_g        = 12.0_fp,                      &
                              MolecRatio    = 6.0_fp,                       &
                              Is_Advected   = T,                            &
                              Is_Gas        = T,                            &
                              Is_Drydep     = F,                            &
                              Is_Wetdep     = F,                            &
                              RC            = RC )

          CASE( 'BR' )
             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              Name          = 'Br',                         &
                              FullName      = 'Atomic bromine',             &
                              MW_g          = 80.0_fp,                      &
                              Is_Advected   = T,                            &
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
                              Name          = 'Br2',                        &
                              FullName      = 'Molecular Bromine',          &
                              MW_g          = 160.0_fp,                     &
                              Is_Advected   = T,                            &
                              Is_Gas        = T,                            &
                              Is_Drydep     = T,                            &
                              Is_Wetdep     = T,                            &
                              DD_F0         = 0.0_fp,                       &
#if defined( NEW_HENRY_CONSTANTS )					    
                              Henry_K0      = 7.20e-3_f8 * To_M_atm,        &
                              Henry_CR      = 4400.0_f8,                    &
#else									    
                              DD_Hstar_old  = 7.60e-1_fp,                   &
                              Henry_K0      = 7.60e01_f8,                   &
                              Henry_CR      = 3720.0_f8,                    &
#endif									    
                              WD_RetFactor  = 0.0_fp,                       &
                              RC            = RC )

          CASE( 'BRCL' )
             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              Name          = 'BrCl',                       &
                              FullName      = 'Bromine chloride',           &
                              MW_g          = 115.0_fp,                     &
                              Is_Advected   = T,                            &
                              Is_Gas        = T,                            &
                              Is_Drydep     = F,                            &
                              Is_Wetdep     = F,                            &
                              RC            = RC )

          CASE( 'BRNO2' )
             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,               &
                              ModelID       = N,                            &
                              Name          = 'BrNO2',                      &
                              FullName      = 'Nitryl bromide',             &
                              MW_g          = 126.0_fp,                     &
                              Is_Advected   = T,                            &
                              Is_Gas        = T,                            &
                              Is_Drydep     = F,                            &
                              Is_Wetdep     = F,                            &
                              RC            = RC )

          CASE( 'BRNO3' )
             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              Name          = 'BrNO3',                      &
                              FullName      = 'Bromine nitrate',            &
                              MW_g          = 142.0_fp,                     &
                              Is_Advected   = T,                            &
                              Is_Gas        = T,                            &
                              Is_Drydep     = T,                            &
                              Is_Wetdep     = F,                            &
                              DD_F0         = 0.0_fp,                       &
                              DD_Hstar_old  = 1.00e+20_fp,                  &
                              RC            = RC )

          CASE( 'BRO' )
             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              Name          = 'BrO',                        &
                              FullName      = 'Bromine monoxide',           &
                              MW_g          = 96.0_fp,                      &
                              Is_Advected   = T,                            &
                              Is_Gas        = T,                            &
                              Is_Drydep     = F,                            &
                              Is_Wetdep     = F,                            &
                              RC            = RC )


          CASE( 'C2H6' )
             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              Name          = NameAllCaps,                  &
                              FullName      = 'Ethane',                     &
                              MW_g          = 30.07_fp,                     &
                              EmMW_g        = 12.0_fp,                      &
                              MolecRatio    = 2.0_fp,                       &
                              Is_Advected   = T,                            &
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
                              Name          = NameAllCaps,                  &
                              FullName      = 'Propane',                    &
                              MW_g          = 44.1_fp,                      &
                              EmMW_g        = 12.0_fp,                      &
                              MolecRatio    = 3.0_fp,                       &
                              Is_Advected   = T,                            &
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
                              Name          = 'CCl4',                       &
                              FullName      = 'Carbon tetrachloride',       &
                              MW_g          = 152.0_fp,                     &
                              Is_Advected   = T,                            &
                              Is_Gas        = T,                            &
                              Is_Drydep     = F,                            &
                              Is_Wetdep     = F,                            &
                              RC            = RC )
             
          CASE( 'CFC11' )
             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              Name          = NameAllCaps,                  &
                              FullName      = 'CFC-11',                     &
                              MW_g          = 137.0_fp,                     &
                              Is_Advected   = T,                            &
                              Is_Gas        = T,                            &
                              Is_Drydep     = F,                            &
                              Is_Wetdep     = F,                            &
                              RC            = RC )

          CASE( 'CFC12' )
             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              Name          = NameAllCaps,                  &
                              FullName      = 'CFC-12',                     &
                              MW_g          = 121.0_fp,                     &
                              Is_Advected   = T,                            &
                              Is_Gas        = T,                            &
                              Is_Drydep     = F,                            &
                              Is_Wetdep     = F,                            &
                              RC            = RC )

          CASE( 'CFCX' )
             FullName = 'CFC-113, CFC-114, CFC-115'
             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              Name          = NameAllCaps,                  &
                              FullName      = FullName,                     &
                              MW_g          = 187.0_fp,                     &
                              Is_Advected   = T,                            &
                              Is_Gas        = T,                            &
                              Is_Drydep     = F,                            &
                              Is_Wetdep     = F,                            &
                              RC            = RC )

          CASE( 'CH2BR2' )
             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              Name          = 'CH2Br2',                     &
                              FullName      = 'Dibromomethane',             &
                              MW_g          = 174.0_fp,                     &
                              Is_Advected   = T,                            &
                              Is_Gas        = T,                            &
                              Is_Drydep     = F,                            &
                              Is_Wetdep     = F,                            &
                              RC            = RC )

          CASE( 'CH2O', 'HCHO' )
             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              Name          = 'CH2O',                       &
                              FullName      = 'Formaldehyde',               &
                              MW_g          = 30.0_fp,                      &
                              Is_Advected   = T,                            &
                              Is_Gas        = T,                            &
                              Is_Drydep     = T,                            &
                              Is_Wetdep     = T,                            &
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
                              Name          = 'CH3Br',                      &
                              FullName      = 'Methyl bromide',             &
                              MW_g          = 95.0_fp,                      &
                              Is_Advected   = T,                            &
                              Is_Gas        = T,                            &
                              Is_Drydep     = F,                            &
                              Is_Wetdep     = F,                            &
                              RC            = RC )

          CASE( 'CH3CL' )
             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              Name          = 'CH3Cl',                      &
                              FullName      = 'Chloromethane',              &
                              MW_g          = 50.0_fp,                      &
                              Is_Advected   = T,                            &
                              Is_Gas        = T,                            &
                              Is_Drydep     = F,                            &
                              Is_Wetdep     = F,                            &
                              RC            = RC )

          CASE( 'CH3CCL3' )
             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              Name          = 'CH3CCl3',                    &
                              FullName      = 'Methyl chloroform',          &
                              MW_g          = 133.0_fp,                     &
                              Is_Advected   = T,                            &
                              Is_Gas        = T,                            &
                              Is_Drydep     = F,                            &
                              Is_Wetdep     = F,                            &
                              RC            = RC )

          CASE( 'CH4' )
             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              Name          = NameAllCaps,                  &
                              FullName      = 'Methane',                    &
                              MW_g          = 16.0_fp,                      &
                              Is_Advected   = T,                            &
                              Is_Gas        = T,                            &
                              Is_Drydep     = F,                            &
                              Is_Wetdep     = F,                            &
                              RC            = RC )

          CASE( 'CHBR3' )
             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              Name          = 'CHBr3',                      &
                              FullName      = 'Bromoform',                  &
                              MW_g          = 253.0_fp,                     &
                              Is_Advected   = T,                            &
                              Is_Gas        = T,                            &
                              Is_Drydep     = F,                            &
                              Is_Wetdep     = F,                            &
                              RC            = RC )

          CASE( 'CL' )
             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              Name          = 'Cl',                         &
                              FullName      = 'Atomic chlorine',            &
                              MW_g          = 35.0_fp,                      &
                              Is_Advected   = T,                            &
                              Is_Gas        = T,                            &
                              Is_Drydep     = F,                            &
                              Is_Wetdep     = F,                            &
                              RC            = RC )

          CASE( 'CL2' )
             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              Name          = 'Cl2',                        &
                              FullName      = 'Molecular chlorine',         &
                              MW_g          = 71.0_fp,                      &
                              Is_Advected   = T,                            &
                              Is_Gas        = T,                            &
                              Is_Drydep     = F,                            &
                              Is_Wetdep     = F,                            &
                              RC            = RC )

          CASE( 'CL2O2' )
             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              Name          = 'Cl2O2',                      &
                              FullName      = 'Dichlorine dioxide',         &
                              MW_g          = 103.0_fp,                     &
                              Is_Advected   = T,                            &
                              Is_Gas        = T,                            &
                              Is_Drydep     = F,                            &
                              Is_Wetdep     = F,                            &
                              RC            = RC )

          CASE( 'CLNO2' )
             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              Name          = 'ClNO2',                      &
                              FullName      = 'Nitryl chloride',            &
                              MW_g          = 81.0_fp,                      &
                              Is_Advected   = T,                            &
                              Is_Gas        = T,                            &
                              Is_Drydep     = F,                            &
                              Is_Wetdep     = F,                            &
                              RC            = RC )
          CASE( 'CLNO3' )
             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              Name          = 'ClNO3',                      &
                              FullName      = 'Chlorine nitrate',           &
                              MW_g          = 97.0_fp,                      &
                              Is_Advected   = T,                            &
                              Is_Gas        = T,                            &
                              Is_Drydep     = F,                            &
                              Is_Wetdep     = F,                            &
                              RC            = RC )

          CASE( 'CLO' )
             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              Name          = 'ClO',                        &
                              FullName      = 'Chlorine monoxide',          &
                              MW_g          = 51.0_fp,                      &
                              Is_Advected   = T,                            &
                              Is_Gas        = T,                            &
                              Is_Drydep     = F,                            &
                              Is_Wetdep     = F,                            &
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
                              Name          = Name,                         &
                              FullName      = 'Chlorine dioxide',           &
                              MW_g          = 67.0_fp,                      &
                              Is_Advected   = T,                            &
                              Is_Gas        = T,                            &
                              Is_Drydep     = F,                            &
                              Is_Wetdep     = F,                            &
                              RC            = RC )

          CASE( 'CO' )
             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              Name          = NameAllCaps,                  &
                              FullName      = 'Carbon monoxide',            &
                              MW_g          = 28.0_fp,                      &
                              Is_Advected   = T,                            &
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
                              Name          = NameAllCaps,                  &
                              FullName      = 'Dimethyl sulfide',           &
                              MW_g          = 62.0_fp,                      &
                              Is_Advected   = T,                            &
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

             ! Dust species are considered to be IN (ice nuclei),
             ! so we allow rainout at all temperatures
             RainEff = (/ 1.0_fp, 1.0_fp, 1.0_fp /)

             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              Name          = NameAllCaps,                  &
                              FullName      = FullName,                     &
                              MW_g          = 29.0_fp,                      &
                              Is_Advected   = T,                            &
                              Is_Gas        = F,                            &
                              Is_Drydep     = T,                            &
                              Is_Wetdep     = T,                            &
                              DD_A_Density  = 2500.0_fp,                    &
                              DD_A_Radius   = 7.3e-7_fp,                    &
                              DD_F0         = 0.0_fp,                       &
                              DD_Hstar_Old  = 0.0_fp,                       &
                              WD_AerScavEff = 1.0_fp,                       &
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

             ! Dust species are considered to be IN (ice nuclei),
             ! so we allow rainout at all temperatures
             RainEff = (/ 1.0_fp, 1.0_fp, 1.0_fp /)

             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              Name          = NameAllCaps,                  &
                              FullName      = FullName,                     &
                              MW_g          = 29.0_fp,                      &
                              Is_Advected   = T,                            &
                              Is_Gas        = F,                            &
                              Is_Drydep     = T,                            &
                              Is_Wetdep     = T,                            &
                              DD_A_Density  = 2650.0_fp,                    &
                              DD_A_Radius   = 1.4e-6_fp,                    &
                              DD_F0         = 0.0_fp,                       &
                              DD_Hstar_Old  = 0.0_fp,                       &
                              WD_AerScavEff = 1.0_fp,                       &
                              WD_CoarseAer  = T,                            &
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

             ! Dust species are considered to be IN (ice nuclei),
             ! so we allow rainout at all temperatures
             RainEff = (/ 1.0_fp, 1.0_fp, 1.0_fp /)

             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              Name          = NameAllCaps,                  &
                              FullName      = FullName,                     &
                              MW_g          = 29.0_fp,                      &
                              Is_Advected   = T,                            &
                              Is_Gas        = F,                            &
                              Is_Drydep     = T,                            &
                              Is_Wetdep     = T,                            &
                              DD_A_Density  = 2650.0_fp,                    &
                              DD_A_Radius   = 2.4e-6_fp,                    &
                              DD_F0         = 0.0_fp,                       &
                              DD_Hstar_Old  = 0.0_fp,                       &
                              WD_AerScavEff = 1.0_fp,                       &
                              WD_CoarseAer  = T,                            &
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

             ! Dust species are considered to be IN (ice nuclei),
             ! so we allow rainout at all temperatures
             RainEff = (/ 1.0_fp, 1.0_fp, 1.0_fp /)

             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              Name          = NameAllCaps,                  &
                              FullName      = FullName,                     &
                              MW_g          = 29.0_fp,                      &
                              Is_Advected   = T,                            &
                              Is_Gas        = F,                            &
                              Is_Drydep     = T,                            &
                              Is_Wetdep     = T,                            &
                              DD_A_Density  = 2650.0_fp,                    &
                              DD_A_Radius   = 4.50e-6_fp,                   &
                              DD_F0         = 0.0_fp,                       &
                              DD_Hstar_Old  = 0.0_fp,                       &
                              WD_AerScavEff = 1.0_fp,                       &
                              WD_CoarseAer  = T,                            &
                              WD_RainoutEff = RainEff,                      &
                              RC            = RC )

          CASE( 'GLYC' )
             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              Name          = NameAllCaps,                  &
                              FullName      = 'Glycoaldehyde',              &
                              MW_g          = 60.0_fp,                      &
                              Is_Advected   = T,                            &
                              Is_Gas        = T,                            &
                              Is_Drydep     = T,                            &
                              Is_Wetdep     = T,                            &
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
                              Name          = NameAllCaps,                  &
                              FullName      = 'Water vapor',                &
                              MW_g          = 18.0_fp,                      &
                              Is_Advected   = T,                            &
                              Is_Gas        = T,                            &
                              Is_Drydep     = F,                            &
                              Is_Wetdep     = F,                            &
                              RC            = RC )

          CASE( 'H2O2' )
             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              Name          = NameAllCaps,                  &
                              FullName      = 'Hydrogen peroxide',          &
                              MW_g          = 34.0_fp,                      &
                              Is_Advected   = T,                            &
                              Is_Gas        = T,                            &
                              Is_Drydep     = T,                            &
                              Is_Wetdep     = T,                            &
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
                              Name          = NameAllCaps,                  &
                              FullName      = 'Hydroxyacetone',             &
                              MW_g          = 74.0_fp,                      &
                              Is_Advected   = T,                            &
                              Is_Gas        = T,                            &
                              Is_Drydep     = T,                            &
                              Is_Wetdep     = F,                            &
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
                              Name          = NameAllCaps,                  &
                              FullName      = 'H-1211',                     &
                              MW_g          = 165.0_fp,                     &
                              Is_Advected   = T,                            &
                              Is_Gas        = T,                            &
                              Is_Drydep     = F,                            &
                              Is_Wetdep     = F,                            &
                              RC            = RC )

          CASE( 'H1301' )
             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              Name          = NameAllCaps,                  &
                              FullName      = 'H-1301',                     &
                              MW_g          = 149.0_fp,                     &
                              Is_Advected   = T,                            &
                              Is_Gas        = T,                            &
                              Is_Drydep     = F,                            &
                              Is_Wetdep     = F,                            &
                              RC            = RC )

          CASE( 'H2402' )
             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              Name          = NameAllCaps,                  &
                              FullName      = 'H-2402',                     &
                              MW_g          = 260.0_fp,                     &
                              Is_Advected   = T,                            &
                              Is_Gas        = T,                            &
                              Is_Drydep     = F,                            &
                              Is_Wetdep     = F,                            &
                              RC            = RC )

          CASE( 'HCFC22' )
             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              Name          = NameAllCaps,                  &
                              FullName      = 'HCFC-22',                    &
                              MW_g          = 86.0_fp,                      &
                              Is_Advected   = T,                            &
                              Is_Gas        = T,                            &
                              Is_Drydep     = F,                            &
                              Is_Wetdep     = F,                            &
                              RC            = RC )

          CASE( 'HCFCX' )
             FullName = 'HCFC-123, HCFC-141b, HCFC-142b'

             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              Name          = NameAllCaps,                  &
                              FullName      = FullName,                     &
                              MW_g          = 117.0_fp,                     &
                              Is_Advected   = T,                            &
                              Is_Gas        = T,                            &
                              Is_Drydep     = F,                            &
                              Is_Wetdep     = F,                            &
                              RC            = RC )

          CASE( 'HCL' )
             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              Name          = 'HCl',                        &
                              FullName      = 'Hydrochloric acid',          &
                              MW_g          = 36.0_fp,                      &
                              Is_Advected   = T,                            &
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
                              Name          = NameAllCaps,                  &
                              FullName      = 'Nitrous acid',               &
                              MW_g          = 47.0_fp,                      &
                              Is_Advected   = T,                            &
                              Is_Gas        = T,                            &
                              Is_Drydep     = F,                            &
                              Is_Wetdep     = F,                            &
                              RC            = RC )

          CASE( 'HNO3' )
             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              Name          = NameAllCaps,                  &
                              FullName      = 'Nitric acid',                &
                              MW_g          = 63.0_fp,                      &
                              Is_Advected   = T,                            &
                              Is_Gas        = T,                            &
                              Is_Drydep     = T,                            &
                              Is_Wetdep     = T,                            &
                              DD_F0         = 0.0_fp,                       &
#if defined( NEW_HENRY_CONSTANTS )                                          
                              Henry_K0      = 2.10e+3_f8 * To_M_atm,        &
                              Henry_CR      = 8700.0_f8,                    &
#else                                                                       
                              DD_Hstar_old  = 1.0e+14_fp,                   &
                              Henry_K0      = 8.3e+4_f8,                    &
                              Henry_CR      = 7400.0_f8,                    &
#endif                                      
                              RC            = RC )

                              ! Save HNO3 drydep index for use further down
                              DryDepID_HNO3 = SpcData(N)%Info%DryDepID
            
          CASE( 'HNO4' )
             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              Name          = NameAllCaps,                  &
                              FullName      = 'Pernitric acid',             &
                              MW_g          = 79.0_fp,                      &
                              Is_Advected   = T,                            &
                              Is_Gas        = T,                            &
                              Is_Drydep     = F,                            &
                              Is_Wetdep     = F,                            &
#if defined( NEW_HENRY_CONSTANTS )					    
                              Henry_K0      = 3.90e1X_f8 * To_M_atm,        &
                              Henry_CR      = 8400.0_f8,                    &
#endif									    
                              RC            = RC )

          CASE( 'HBR' )
             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              Name          = 'HBr',                        &
                              FullName      = 'Hypobromic acid',            &
                              MW_g          = 81.0_fp,                      &
                              Is_Advected   = T,                            &
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
                              Name          = 'HOBr',                       &
                              FullName      = 'Hypobromous acid',           &
                              MW_g          = 97.0_fp,                      &
                              Is_Advected   = T,                            &
                              Is_Gas        = T,                            &
                              Is_Drydep     = T,                            &
                              Is_Wetdep     = T,                            &
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
                              Name          = 'HOCl',                       &
                              FullName      = 'Hypochlorous acid',          &
                              MW_g          = 52.0_fp,                      &
                              Is_Advected   = T,                            &
                              Is_Gas        = T,                            &
                              Is_Drydep     = F,                            &
                              Is_Wetdep     = F,                            &
                              RC            = RC )

          CASE( 'ISOP' )
             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              Name          = NameAllCaps,                  &
                              FullName      = 'Isoprene',                   &
                              MW_g          = 68.12_fp,                     &
                              EmMW_g        = 12.0_fp,                      &
                              MolecRatio    = 5e+0_fp,                      &
                              Is_Advected   = T,                            &
                              Is_Gas        = T,                            &
                              Is_Drydep     = F,                            &
                              Is_Wetdep     = F,                            &
#if defined( NEW_HENRY_CONSTANTS )                                          
                              Henry_K0      = 3.40e-4_f8 * To_M_atm,        &
                              Henry_CR      = 4400.0_f8,                    &
#endif                                      
                              RC            = RC )

          CASE( 'IEPOX' )
             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              Name          = NameAllCaps,                  &
                              FullName      = 'Isoprene epoxide',           &
                              MW_g          = 118.0_fp,                     &
                              Is_Advected   = T,                            &
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

!          CASE( 'ISN1' )
!             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
!                              ThisSpc       = SpcData(N)%Info,              &
!                              ModelID       = N,                            &
!                              DryDepID      = DryDepID_HNO3
!                              Name          = NameAllCaps,                  &
!                              FullName      = '',                           &
!                              MW_g          = __.0_fp,                      &
!                              MolecRatio    = 1.0_fp,                       &
!                              Is_Advected   = T,                            &
!                              Is_Gas        = T,                            &
!                              Is_Drydep     = T,                            &
!                              Is_Wetdep     = T,                            &
!                              DD_F0         = 0.0_fp,                       &
!                              DD_Hstar_old  = 0.0_fp,                      &
!                              RC            = RC )

          CASE( 'ISOA1', 'ISOA2', 'ISOA3' )
             FullName = 'Lumped semivolatile gas products of isoprene oxidation'

             ! Rainout efficiency is 0.8 because this is a SOA species.
             ! Turn off rainout when 237 K < T < 258 K.
             RainEff = (/ 0.8_fp, 0.0_fp, 0.8_fp /)

             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              Name          = NameAllCaps,                  &
                              FullName      = FullName,                     &
                              MW_g          = 150.0_fp,                     &
                              Is_Advected   = T,                            &
                              Is_Gas        = F,                            &
                              Is_Drydep     = T,                            &
                              Is_Wetdep     = T,                            &
                              DD_F0         = 0.0_fp,                       &
                              DD_HStar_old  = 0.0_fp,                       &
                              WD_AerScavEff = 0.8_fp,                       &
                              RC            = RC )

          CASE( 'ISOG1', 'ISOG2', 'ISOG3' )
             FullName = 'Lumped semivolatile gas products of isoprene oxidation'

             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              Name          = NameAllCaps,                  &
                              FullName      = FullName,                     &
                              MW_g          = 150.0_fp,                     &
                              Is_Advected   = T,                            &
                              Is_Gas        = T,                            &
                              Is_Drydep     = T,                            &
                              Is_Wetdep     = T,                            &
                              DD_F0         = 0.0_fp,                       &
                              DD_Hstar_old  = 1.0e+5_fp,                    &
                              WD_RetFactor  = 2.0e-2_fp,                    &
                              RC            = RC )

          CASE( 'ISOPN' )

             ! First create isoprene nitrate species name
             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              Name          = NameAllCaps,                  &
                              FullName      = 'Isoprene nitrate',           &
                              MW_g          = 147.0_fp,                     &
                              Is_Advected   = T,                            &
                              Is_Gas        = T,                            &
                              Is_Drydep     = T,                            &
                              Is_Wetdep     = T,                            &
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
        
! Leave for future expansion 
!          CASE( 'ISOPNB', 'ISOPND'  )
!             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
!                              ThisSpc       = SpcData(N)%Info,              &
!                              ModelID       = N,                            &
!                              Name          = NameAllCaps,                  &
!                              FullName      = 'Isoprene nitrate',           &
!                              MW_g          = 147.0_fp,                     &
!                              Is_Advected   = F,                            &
!                              Is_Gas        = T,                            &
!                              Is_Drydep     = T,                            &
!                              Is_Wetdep     = F,                            &
!                              DD_F0         = 1.0_fp,                       &
!                              DD_Hstar_old  = 1.70e+4_fp,                   &
!                              RC            = RC )

          CASE( 'LIMO' )
             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              Name          = NameAllCaps,                  &
                              FullName      = 'Limonene',                   &
                              MW_g          = 136.23_fp,                    &
                              MolecRatio    = 1.0_fp,                       &
                              Is_Advected   = T,                            &
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
                              Name          = NameAllCaps,                  &
                              FullName      = 'Methacrolein',               &
                              MW_g          = 70.0_fp,                      &
                              Is_Advected   = T,                            &
                              Is_Gas        = T,                            &
                              Is_Drydep     = T,                            &
                              Is_Wetdep     = F,                            &
                              DD_F0         = 1.0_fp,                       &
#if defined( NEW_HENRY_CONSTANTS )                                          
                              Henry_K0      = 4.8e-2_f8 * To_M_atm,         &
                              Henry_CR      = 4300.0_f8,                    &
#else                                                                       
                              DD_Hstar_old  = 6.5e+0_fp,                    &
#endif
                              RC            = RC )

          CASE( 'MAP' )
             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              Name          = NameAllCaps,                  &
                              FullName      = 'Peroxyacetic acid',          &
                              MW_g          = 76.0_fp,                      &
                              Is_Advected   = T,                            &
                              Is_Gas        = T,                            &
                              Is_Drydep     = T,                            &
                              Is_Wetdep     = T,                            &
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
                              Name          = NameAllCaps,                  &
                              FullName      = 'Peroxyacetyl nitrate',       &
                              MW_g          = 72.11_fp,                     &
                              EmMW_g        = 12.0_fp,                      &
                              MolecRatio    = 4.0_fp,                       &
                              Is_Advected   = T,                            &
                              Is_Gas        = T,                            &
                              Is_Drydep     = F,                            &
                              Is_Wetdep     = F,                            &
#if defined( NEW_HENRY_CONSTANTS )                                          
                              Henry_K0      = 2.90e+02_f8 * To_M_atm,       &
                              Henry_CR      = 5700.0_f8,                    &
#endif                                      
                              RC            = RC )

          CASE( 'MMN' )
             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              Name          = NameAllCaps,                  &
                              FullName      = 'Nitrate from MACR + MVK',    &
                              MW_g          = 149.0_fp,                     &
                              Is_Advected   = T,                            &
                              Is_Gas        = T,                            &
                              Is_Drydep     = T,                            &
                              Is_Wetdep     = T,                            &
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

          CASE( 'MOBA' )
             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              Name          = NameAllCaps,                  &
                              FullName      = '5C acid from isoprene',      &
                              MW_g          = 114.0_fp,                     &
                              Is_Advected   = T,                            &
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
             
             ! Turn off rainout when 237 K <= T < 258 K
             RainEff = (/ 1.0_fp, 0.0_fp, 1.0_fp /)

             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              Name          = NameAllCaps,                  &
                              FullName      = 'Hydrophilic marine OC',      &
                              MW_g          = 12.01_fp,                     &
                              EmMW_g        = 12.0_fp,                      &
                              MolecRatio    = 1.0_fp,                       &
                              Is_Advected   = T,                            &
                              Is_Gas        = T,                            &
                              Is_Drydep     = F,                            &
                              Is_Wetdep     = F,                            &
                              DD_Hstar_old  = 0.0_fp,                       &
                              WD_AerScavEff = 1.0_fp,                       &
                              WD_RainoutEff = RainEff,                      &
                              RC            = RC )

          CASE( 'MOPO' )

             ! Turn off rainout because MOPO is hydrophobic
             RainEff = (/ 0.0_fp, 0.0_fp, 0.0_fp /)

             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              Name          = NameAllCaps,                  &
                              FullName      = 'Hydrophobic marine OC',      &
                              MW_g          = 12.01_fp,                     &
                              EmMW_g        = 12.0_fp,                      &
                              MolecRatio    = 1.0_fp,                       &
                              Is_Advected   = T,                            &
                              Is_Gas        = T,                            &
                              Is_Drydep     = F,                            &
                              Is_Wetdep     = F,                            &
                              DD_Hstar_old  = 0.0_fp,                       &
                              WD_AerScavEff = 0.0_fp,                       &
                              WD_RainOutEff = RainEff,                      &
                              RC            = RC )

          CASE( 'MP', 'CH3OOH' )
             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              Name          = 'MP',                         &
                              FullName      = 'Methyl hydro peroxide',      &
                              MW_g          = 48.0_fp,                      &
                              Is_Advected   = T,                            &
                              Is_Gas        = T,                            &
                              Is_Drydep     = F,                            &
                              Is_Wetdep     = T,                            &
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
                              Name          = NameAllCaps,                  &
                              FullName      = 'Methyl peroxy nitrate',      &
                              MW_g          = 93.0_fp,                      &
                              Is_Advected   = T,                            &
                              Is_Gas        = T,                            &
                              Is_Drydep     = F,                            &
                              Is_Wetdep     = F,                            &
                              RC            = RC )

          CASE( 'MSA' )

             ! Turn off rainout when 237 K <= T < 258 K
             RainEff = (/ 1.0_fp, 0.0_fp, 1.0_fp /)

             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              Name          = NameAllCaps,                  &
                              FullName      = 'Methyl sulfonic acid',       &
                              MW_g          = 96.0_fp,                      &
                              Is_Advected   = T,                            &
                              Is_Gas        = F,                            &
                              Is_Drydep     = T,                            &
                              Is_Wetdep     = T,                            &
                              DD_F0         = 0.0_fp,                       &
                              DD_Hstar_old  = 0.0_fp,                       &
                              WD_AerScavEff = 1.0_fp,                       &
                              WD_RainoutEff = RainEff,                      &
                              RC            = RC )

          CASE( 'MTPA' )
             FullName = 'a-pinene, b-pinene, sabinene, carene'
             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              Name          = NameAllCaps,                  &
                              FullName      = FullName,                     &
                              MW_g          = 136.23_fp,                    &
                              Is_Advected   = T,                            &
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
                              Name          = NameAllCaps,                  &
                              FullName      = FullName,                     &
                              MW_g          = 136.23_fp,                    &
                              MolecRatio    = 1.0_fp,                       &
                              Is_Advected   = T,                            &
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
                              Name          = NameAllCaps,                  &
                              FullName      = 'Methyl vinyl ketone',        &
                              MW_g          = 70.0_fp,                      &
                              Is_Advected   = T,                            &
                              Is_Gas        = T,                            &
                              Is_Drydep     = T,                            &
                              Is_Wetdep     = F,                            &
                              DD_F0         = 1.0_fp,                       &
#if defined( NEW_HENRY_CONSTANTS )                                          
                              Henry_K0      = 2.6e-1_f8 * To_M_atm,         &
                              Henry_CR      = 4800.0_f8,                    &
#else                                                                       
                              DD_Hstar_old  = 4.4e+1_fp,                    &
#endif                                      
                              RC            = RC )

          CASE( 'NAP' )
             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              Name          = NameAllCaps,                  &
                              FullName      = 'Naphtalene/IVOC surrogate',  &
                              MW_g          = 128.27_fp,                    &
                              EmMw_g        = 12.0_fp,                      &
                              MolecRatio    = 10.0_fp,                      &
                              Is_Advected   = T,                            &
                              Is_Gas        = T,                            &
                              Is_Drydep     = F,                            &
                              Is_Wetdep     = F,                            &
                              RC            = RC )

          CASE( 'N2O' )
             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              Name          = NameAllCaps,                  &
                              FullName      = 'Nitrous oxide',              &
                              MW_g          = 44.0_fp,                      &
                              Is_Advected   = T,                            &
                              Is_Gas        = T,                            &
                              Is_Drydep     = F,                            &
                              Is_Wetdep     = F,                            &
                              RC            = RC )

          CASE( 'N2O5' )
             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              DryDepID      = DryDepID_HNO3,                &
                              Name          = NameAllCaps,                  &
                              FullName      = 'Dinitrogen pentoxide',       &
                              MW_g          = 105.0_fp,                     &
                              Is_Advected   = T,                            &
                              Is_Gas        = T,                            &
                              Is_Drydep     = T,                            &
                              Is_Wetdep     = F,                            &
                              DD_F0         = 0.0_fp,                       &
#if defined( NEW_HENRY_CONSTANTS )					    
                              Henry_K0      = 2.10e-2_f8 * To_M_atm,        &
                              Henry_CR      = 3400.0_f8,                    &
#else									    
                              DD_Hstar_old  = 0.0_fp,                       &
#endif									    
                              RC            = RC )

          CASE( 'NH3' )
             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              Name          = NameAllCaps,                  &
                              FullName      = 'Ammonia',                    &
                              MW_g          = 17.0_fp,                      &
                              Is_Advected   = T,                            &
                              Is_Gas        = T,                            &
                              Is_Drydep     = T,                            &
                              Is_Wetdep     = T,                            &
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

             ! Turn off rainout when 237 K <= T < 258 K
             RainEff = (/ 1.0_fp, 0.0_fp, 1.0_fp /)

             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              Name          = NameAllCaps,                  &
                              FullName      = 'Ammonium',                   &
                              MW_g          = 18.0_fp,                      &
                              Is_Advected   = T,                            &
                              Is_Gas        = F,                            &
                              Is_Drydep     = T,                            &
                              Is_Wetdep     = T,                            &
                              DD_F0         = 0.0_fp,                       &
                              DD_Hstar_old  = 0.0_fp,                       &
                              WD_AerScavEff = 1.0_fp,                       &
                              WD_RainoutEff = RainEff,                      &
                              RC            = RC )

          CASE( 'NIT' )

             ! Turn off rainout when 237 K <= T < 258 K
             RainEff = (/ 1.0_fp, 0.0_fp, 1.0_fp /)

             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              Name          = NameAllCaps,                  &
                              FullName      = 'Inorganic nitrates',         &
                              MW_g          = 62.0_fp,                      &
                              Is_Advected   = T,                            &
                              Is_Gas        = F,                            &
                              Is_Drydep     = T,                            &
                              Is_Wetdep     = T,                            &
                              DD_F0         = 0.0_fp,                       &
                              DD_Hstar_Old  = 0.0_fp,                       &
                              WD_AerScavEff = 1.0_fp,                       &
                              WD_RainoutEff = RainEff,                      &
                              RC            = RC )

          CASE( 'NITS' )
             
             A_Radius = ( Input_Opt%SALC_REDGE_um(1) +                      &
                          Input_Opt%SALC_REDGE_um(2)  ) * 0.5e-6_fp
             Fullname = 'Inorganic nitrates on surface of seasalt aerosol'

             ! Turn off rainout when 237 K <= T < 258 K
             RainEff  = (/ 1.0_fp, 0.0_fp, 1.0_fp /)

             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              Name          = 'NITs',                       &
                              FullName      = FullName,                     &
                              MW_g          = 62.0_fp,                      &
                              Is_Advected   = T,                            &
                              Is_Gas        = F,                            &
                              Is_Drydep     = T,                            &
                              Is_Wetdep     = T,                            &
                              DD_A_Density  = 2200.0_fp,                    &
                              DD_A_Radius   = A_Radius,                     &
                              DD_F0         = 0.0_fp,                       &
                              DD_Hstar_Old  = 0.0_fp,                       &
                              WD_AerScavEff = 1.0_fp,                       &
                              WD_RainoutEff = RainEff,                      &
                              RC            = RC )

          CASE( 'NO' )
             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              Name          = NameAllCaps,                  &
                              FullName      = 'Nitrogen oxide',             &
                              MW_g          = 30.0_fp,                      &
                              Is_Advected   = T,                            &
                              Is_Gas        = T,                            &
                              Is_Drydep     = F,                            &
                              Is_Wetdep     = F,                            &
#if defined( NEW_HENRY_CONSTANTS )                                          
                              Henry_K0      = 1.90e-5_f8 * To_M_atm,        &
                              Henry_CR      = 1600.0_f8,                    &
#endif                                      
                              RC            = RC )

          CASE( 'NO2' )
             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              Name          = NameAllCaps,                  &
                              FullName      = 'Nitrogen dioxide',           &
                              MW_g          = 46.0_fp,                      &
                              Is_Advected   = T,                            &
                              Is_Gas        = T,                            &
                              Is_Drydep     = T,                            &
                              Is_Wetdep     = F,                            &
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
                              Name          = NameAllCaps,                  &
                              FullName      = 'Nitrate radical',            &
                              MW_g          = 62.0_fp,                      &
                              MolecRatio    = 1.0_fp,                       &
                              Is_Advected   = T,                            &
                              Is_Gas        = T,                            &
                              Is_Drydep     = F,                            &
                              Is_Wetdep     = F,                            &
#if defined( NEW_HENRY_CONSTANTS )					    
                              Henry_K0      = 3.80e-4_f8 * To_M_atm,        &
                              Henry_CR      = 1900.0_f8,                    &
#endif									    
                              RC            = RC )


          CASE( 'O3' )
             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              Name          = NameAllCaps,                  &
                              FullName      = 'Ozone',                      &
                              MW_g          = 48.0_fp,                      &
                              Is_Advected   = T,                            &
                              Is_Gas        = T,                            &
                              Is_Drydep     = T,                            &
                              Is_Wetdep     = F,                            &
                              DD_F0         = 1.0_fp,                       &
#if defined( NEW_HENRY_CONSTANTS )                                          
                              Henry_K0      = 1.90e-05_f8 * To_M_atm,       &
                              Henry_CR      = 1600.0_f8,                    &
#else                                                                       
                              DD_Hstar_old  = 1.0e-2_fp,                    &
#endif                                      
                              RC            = RC )

          CASE( 'OCPI', 'OCPO' )
             
             ! These have identical properties except for 
             ! the names and rainout efficiencies
             !
             ! Turn off rainout for hydrophilic OC when 237 K <= T < 258 K.
             ! Turn off rainout for hydrophobic OC, for all temperatures.
             SELECT CASE( NameAllCaps )
                CASE( 'OCPI' )
                   FullName = 'Hydrophilic organic carbon aerosol'
                   RainEff  = (/ 1.0_fp, 0.0_fp, 1.0_fp /)
                CASE( 'OCPO' )
                   Fullname = 'Hydrophobic organic carbon aerosol'
                   RainEff  = (/ 0.0_fp, 0.0_fp, 0.0_fp /)
             END SELECT

             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              Name          = NameAllCaps,                  &
                              FullName      = FullName,                     &
                              MW_g          = 12.01_fp,                     &
                              EmMW_g        = 12.0_fp,                      &
                              Is_Advected   = T,                            &
                              Is_Gas        = F,                            &
                              Is_Drydep     = T,                            &
                              Is_Wetdep     = T,                            &
                              DD_F0         = 0.0_fp,                       &
                              DD_Hstar_Old  = 0.0_fp,                       &
                              WD_AerScavEff = 1.0_fp,                       &
                              WD_RainoutEff = RainEff,                      &
                              RC            = RC )

          CASE( 'OPOA1', 'OPOA2' )
             FullName = 'Lumped aerosol product of SVOC oxidation'

             ! Rainout efficiency is 0.8 because these are SOA species.
             ! Turn off rainout when when 237 K <= T < 258 K.
             RainEff  = (/ 0.8_fp, 0.0_fp, 0.8_fp /) 

             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              Name          = NameAllCaps,                  &
                              FullName      = FullName,                     &
                              MW_g          = 12.01_fp,                     &
                              EmMW_g        = 12.0_fp,                      &
                              MolecRatio    = 1.0_fp,                       &
                              Is_Advected   = T,                            &
                              Is_Gas        = F,                            &
                              Is_Drydep     = T,                            &
                              Is_Wetdep     = T,                            &
                              DD_F0         = 0.0_fp,                       &
                              DD_Hstar_old  = 0.0_fp,                       &
                              WD_AerScavEff = 0.8_fp,                       &
                              WD_RainoutEff = RainEff,                      &
                              RC            = RC )

          CASE( 'OPOG1', 'OPOG2' )
             FullName = 'Lumped gas product of SVOC oxidation'
             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              Name          = NameAllCaps,                  &
                              FullName      = FullName,                     &
                              MW_g          = 12.01_fp,                     &
                              EmMW_g        = 12.0_fp,                      &
                              MolecRatio    = 1.0_fp,                       &
                              Is_Advected   = T,                            &
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
                              Name          = NameAllCaps,                  &
                              FullName      = 'Carbonyl sulfide',           &
                              MW_g          = 60.0_fp,                      &
                              Is_Advected   = T,                            &
                              Is_Gas        = T,                            &
                              Is_Drydep     = F,                            &
                              Is_Wetdep     = F,                            &
                              RC            = RC )

          CASE( 'PAN' )
             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              Name          = NameAllCaps,                  &
                              FullName      = 'Peroxyacetyl nitrate',       &
                              MW_g          = 121.0_fp,                     &
                              Is_Advected   = T,                            &
                              Is_Gas        = T,                            &
                              Is_Drydep     = T,                            &
                              Is_Wetdep     = F,                            &
                              DD_F0         = 1.0e+0_fp,                    &
#if defined( NEW_HENRY_CONSTANTS )                                          
                              Henry_K0      = 2.90e+02_f8 * To_M_atm,       &
                              Henry_CR      = 5700.0_f8,                    &
#else                                                                       
                              DD_Hstar_old  = 3.60e+0_fp,                   &
#endif                                      
                              RC            = RC )

                              ! Save PAN drydep index for use further down
                              DryDepID_PAN = SpcData(N)%Info%DryDepID

          CASE( 'PMN' )
             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              Name          = NameAllCaps,                  &
                              FullName      = 'Peroxymethacroyl nitrate',   &
                              MW_g          = 147.0_fp,                     &
                              Is_Advected   = T,                            &
                              Is_Gas        = T,                            &
                              Is_Drydep     = T,                            &
                              Is_Wetdep     = F,                            &
                              DryDepId      = DryDepID_PAN,                 &
                              DD_F0         = 0.0_fp,                       &
#if defined( NEW_HENRY_CONSTANTS )                                          
                              Henry_K0      = 1.70e-2_f8 * To_M_atm,        &
#else                                                                       
                              DD_Hstar_old  = 0.0_fp,                       &
#endif
                              RC            = RC )

          CASE( 'PPN' )
             FullName = 'Lumped peroxypropionyl nitrate'
             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              DryDepID      = DryDepID_PAN,                 &
                              Name          = NameAllCaps,                  &
                              FullName      = FULLNAME,                     &
                              MW_g          = 135.0_fp,                     &
                              Is_Advected   = T,                            &
                              Is_Gas        = T,                            &
                              Is_Drydep     = T,                            &
                              Is_Wetdep     = F,                            &
                              DD_F0         = 0.0_fp,                       &
#if defined( NEW_HENRY_CONSTANTS )                                          
                              Henry_K0      = 2.9e-2_f8 * To_M_atm,         &
                              Henry_CR      = _f8,                          &
#else                                                                       
                              DD_Hstar_old  = 0.0_fp,                       &
#endif                                                                      
                              RC            = RC )

          CASE( 'POA1', 'POA2' )

             ! Turn off rainout because these are hydrophobic species.
             RainEff = (/ 0.0_fp, 0.0_fp, 0.0_fp /)

             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              Name          = NameAllCaps,                  &
                              FullName      ='Lumped aerosol primary SVOCs',&
                              MW_g          = 12.01_fp,                     &
                              EmMW_g        = 12.0_fp,                      &
                              MolecRatio    = 1.0_fp,                       &
                              Is_Advected   = T,                            &
                              Is_Gas        = F,                            &
                              Is_Drydep     = T,                            &
                              Is_Wetdep     = T,                            &
                              DD_F0         = 0.0_fp,                       &
                              DD_Hstar_old  = 0.0_fp,                       &
                              WD_AerScavEff = 1.0_fp,                       &
                              WD_RainoutEff = RainEff,                      &
                              RC            = RC )

          CASE( 'POG1', 'POG2' )
             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              Name          = NameAllCaps,                  &
                              FullName      = 'Lumped gas primary SVOCs',   &
                              MW_g          = 12.01_fp,                     &
                              EmMW_g        = 12.0_fp,                      &
                              MolecRatio    = 1.0_fp,                       &
                              Is_Advected   = T,                            &
                              Is_Gas        = T,                            &
                              Is_Drydep     = T,                            &
                              Is_Wetdep     = T,                            &
                              DD_F0         = 0.0_fp,                       &
                              DD_Hstar_old  = 9.50e+0_f8,                   &
                              Henry_K0      = 9.50e+0_f8,                   &
                              Henry_CR      = 4700.0_f8,                    &
                              WD_RetFactor  = 2.0e-2_fp,                    &
                              RC            = RC )

          CASE( 'PRPE' )
             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              Name          = NameAllCaps,                  &
                              FullName      = 'Lumped >= C3 alkenes',       &
                              MW_g          = 42.08_fp,                     &
                              EmMW_g        = 12.0_fp,                      &
                              MolecRatio    = 3.0_fp,                       &
                              Is_Advected   = T,                            &
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
                              Name          = NameAllCaps,                  &
                              FullName      = 'Propanone nitrate',          &
                              MW_g          = 119.0_fp,                     &
                              Is_Advected   = T,                            &
                              Is_Gas        = T,                            &
                              Is_Drydep     = T,                            &
                              Is_Wetdep     = T,                            &
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
             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              DryDepID      = DryDepId_PAN,                 &
                              Name          = NameAllCaps,                  &
                              FullName      = 'Lumped alkyl nitrate',       &
                              MW_g          = 119.0_fp,                     &
                              Is_Advected   = T,                            &
                              Is_Gas        = T,                            &
                              Is_Drydep     = T,                            &
                              Is_Wetdep     = F,                            &
                              DD_F0         = 0.0_fp,                       &
#if defined( NEW_HENRY_CONSTANTS )                                          
                              Henry_K0      = 1.0e-2_f8 * To_M_atm,         &
                              Henry_CR      = 5800.0_f8,                    &
#else                                                                       
                              DD_Hstar_old  = 0.0_fp,                       &
#endif
                              RC            = RC )


          CASE( 'RCHO' )
             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              Name          = NameAllCaps,                  &
                              FullName      = 'Lumped aldehyde >= C3',      &
                              MW_g          = 58.0_fp,                      &
                              Is_Advected   = T,                            &
                              Is_Gas        = T,                            &
                              Is_Drydep     = F,                            &
                              Is_Wetdep     = F,                            &
#if defined( NEW_HENRY_CONSTANTS )                                          
                              Henry_K0      = 9.5e-2_f8 * To_M_atm,         &
                              Henry_CR      = 6200.0_f8,                    &
#endif
                              RC            = RC )

          CASE( 'RIP' )
             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              Name          = NameAllCaps,                  &
                              FullName      = 'Peroxide from RIO2',         &
                              MW_g          = 118.0_fp,                     &
                              Is_Advected   = T,                            &
                              Is_Gas        = T,                            &
                              Is_Drydep     = T,                            &
                              Is_Wetdep     = T,                            &
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
             
             A_Radius = ( Input_Opt%SALA_REDGE_um(1) +                      &
                          Input_Opt%SALA_REDGE_um(2)  ) * 0.5e-6_fp
             FullName = 'Accumulation mode sea salt aerosol'

             ! Turn off rainout when 237 K <= T < 258 K
             RainEff = (/ 1.0_fp, 0.0_fp, 1.0_fp /)   

             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              Name          = NameAllCaps,                  &
                              FullName      = Fullname,                     &
                              MW_g          = 31.4_fp,                      &
                              Is_Advected   = T,                            &
                              Is_Gas        = F,                            &
                              Is_Drydep     = T,                            &
                              Is_Wetdep     = T,                            &
                              DD_A_Density  = 2200.0_fp,                    &
                              DD_A_Radius   = A_Radius,                     &
                              DD_F0         = 0.0_fp,                       &
                              DD_Hstar_old  = 0.0_fp,                       &
                              WD_AerScavEff = 1.0_fp,                       &
                              WD_RainoutEff = RainEff,                      &
                              RC            = RC )

          CASE( 'SALC' )
             
             A_Radius = ( Input_Opt%SALC_REDGE_um(1) +                      &
                          Input_Opt%SALC_REDGE_um(2)  ) * 0.5e-6_fp
             FullName = 'Coarse mode sea salt aerosol'
 
             ! Turn off rainout when 237 K <= T < 258 K
             RainEff = (/ 1.0_fp, 0.0_fp, 1.0_fp /) 

             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              Name          = NameAllCaps,                  &
                              FullName      = Fullname,                     &
                              MW_g          = 31.4_fp,                      &
                              Is_Advected   = T,                            &
                              Is_Gas        = F,                            &
                              Is_Drydep     = T,                            &
                              Is_Wetdep     = T,                            &
                              DD_A_Density  = 2200.0_fp,                    &
                              DD_A_Radius   = A_Radius,                     &
                              DD_F0         = 0.0_fp,                       &
                              DD_Hstar_old  = 0.0_fp,                       &
                              WD_AerScavEff = 1.0_fp,                       &
                              WD_CoarseAer  = T,                            &
                              WD_RainoutEff = RainEff,                      &
                              RC            = RC )

          CASE( 'SESQ' )
             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              Name          = NameAllCaps,                  &
                              FullName      = 'Sesquiterpene',              &
                              MW_g          = 150.0_fp,                     &
                              Is_Advected   = T,                            &
                              Is_Gas        = T,                            &
                              Is_Drydep     = F,                            &
                              Is_Wetdep     = F,                            &
                              RC            = RC )

          CASE( 'SO2' )

             ! SO2 wet-deposits like an aerosol. 
             ! Turn off rainout when 237 K < T < 258K.
             RainEff = (/ 1.0_fp, 0.0_fp, 1.0_fp /)

             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              Name          = NameAllCaps,                  &
                              FullName      = 'Sulfur dioxide',             &
                              MW_g          = 64.0_fp,                      &
                              Is_Advected   = T,                            &
                              Is_Gas        = T,                            &
                              Is_Drydep     = T,                            &
                              Is_Wetdep     = T,                            &
                              DD_F0         = 0.0_fp,                       &
#if defined( NEW_HENRY_CONSTANTS )					    
                              Henry_K0      = 1.30e-2_f8 * To_M_atm,        &
                              Henry_CR      = 2900.0_f8,                    &
#else									    
                              DD_Hstar_old  = 1.00e+5_fp,                   &
#endif									    
                              WD_RainoutEff = RainEff,                      &
                              RC            = RC )

          CASE( 'SO4' )

             ! Turn off rainout when 237 K <= T < 258 K
             RainEff = (/ 1.0_fp, 0.0_fp, 1.0_fp /)   

             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              Name          = 'SO4',                        &
                              FullName      = 'Sulfate',                    &
                              MW_g          = 96.0_fp,                      &
                              Is_Advected   = T,                            &
                              Is_Gas        = F,                            &
                              Is_Drydep     = T,                            &
                              Is_Wetdep     = T,                            &
                              DD_F0         = 0.0_fp,                       &
                              DD_Hstar_Old  = 0.0_fp,                       &
                              WD_AerScavEff = 1.0_fp,                       &
                              WD_RainoutEff = RainEff,                      &
                              RC            = RC )

          CASE( 'SO4S' )
             
             A_Radius = ( Input_Opt%SALC_REDGE_um(1) +                      &
                          Input_Opt%SALC_REDGE_um(2)  ) * 0.5e-6_fp
             Fullname = 'Sulfate on surface of seasalt aerosol'

             ! Turn off rainout when 237 K <= T < 258 K
             RainEff = (/ 1.0_fp, 0.0_fp, 1.0_fp /)   

             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              Name          = 'SO4s',                       &
                              FullName      = FullName,                     &
                              MW_g          = 96.0_fp,                      &
                              Is_Advected   = T,                            &
                              Is_Gas        = F,                            &
                              Is_Drydep     = T,                            &
                              Is_Wetdep     = T,                            &
                              DD_A_Density  = 2200.0_fp,                    &
                              DD_A_Radius   = A_Radius,                     &
                              DD_F0         = 0.0_fp,                       &
                              DD_Hstar_Old  = 0.0_fp,                       &
                              WD_AerScavEff = 1.0_fp,                       &
                              WD_RainoutEff = RainEff,                      &
                              RC            = RC )

          CASE( 'TOLU' )
             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              Name          = NameAllCaps,                  &
                              FullName      = 'Toluene',                    &
                              MW_g          = 92.14_fp,                     &
                              EmMW_g        = 12.0_fp,                      &
                              MolecRatio    = 7.0_fp,                       &
                              Is_Advected   = T,                            &
                              Is_Gas        = T,                            &
                              Is_Drydep     = F,                            &
                              Is_Wetdep     = F,                            &
                              RC            = RC )

          CASE( 'TSOA0', 'TSOA1', 'TSOA2', 'TSOA3' )
             FullName = 'Lumped semivolatile aerosol products of monoterpene + sesquiterpene oxidation'

             ! Rainout efficiency is 0.8 because these are SOA species.
             ! Turn off rainout when 237 K <= T < 258 K.
             RainEff = (/ 0.8_fp, 0.0_fp, 0.8_fp /)   

             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              Name          = NameAllCaps,                  &
                              FullName      = FullName,                     &
                              MW_g          = 150.0_fp,                     &
                              Is_Advected   = T,                            &
                              Is_Gas        = F,                            &
                              Is_Drydep     = T,                            &
                              Is_Wetdep     = T,                            &
                              DD_F0         = 0.0_fp,                       &
                              DD_Hstar_old  = 0.0_fp,                       &
                              WD_AerScavEff = 0.8_fp,                       &
                              WD_RainoutEff = RainEff,                      &
                              RC            = RC )

          CASE( 'TSOG0', 'TSOG1', 'TSOG2', 'TSOG3' )
             FullName = 'Lumped semivolatile gas products of monoterpene + sesquiterpene oxidation'
             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              Name          = NameAllCaps,                  &
                              FullName      = FullName,                     &
                              MW_g          = 150.0_fp,                     &
                              Is_Advected   = T,                            &
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
                              Name          = NameAllCaps,                  &
                              FullName      = 'Xylene',                     &
                              MW_g          = 106.16_fp,                    &
                              EmMW_g        = 12.0_fp,                      &
                              MolecRatio    = 8.0_fp,                       &
                              Is_Advected   = T,                            &
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
                              Is_Advected   = T,                            &
                              Is_Gas        = F,                            &
                              Is_Drydep     = F,                            &
                              Is_Wetdep     = F,                            &
                              RC            = RC )

          CASE( 'PB', '210PB', 'PB210' )

             ! Turn off rainout when 237 K <= T < 258 K
             RainEff = (/ 1.0_fp, 0.0_fp, 1.0_fp /)   

             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              Name          = 'Pb',                         &
                              FullName      = 'Lead-210 isotope',           &
                              MW_g          = 210.0_fp,                     &
                              Is_Advected   = T,                            &
                              Is_Gas        = F,                            &
                              Is_Drydep     = T,                            &
                              Is_Wetdep     = T,                            &
                              DD_F0         = 0.0_fp,                       &
                              DD_HStar_old  = 0.0_fp,                       &
                              WD_AerScavEff = 1.0_fp,                       &
                              WD_RainoutEff = RainEff,                      &
                              RC            = RC )

          CASE( 'BE', '7BE', 'BE7' )

             ! Turn off rainout when 237 K <= T < 258 K
             RainEff = (/ 1.0_fp, 0.0_fp, 1.0_fp /)   

             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              Name          = 'Be7',                        &
                              FullName      = 'Beryllium-7 isotope',        &
                              MW_g          = 7.0_fp,                       &
                              Is_Advected   = T,                            &
                              Is_Gas        = F,                            &
                              Is_Drydep     = T,                            &
                              Is_Wetdep     = T,                            &
                              DD_F0         = 0.0_fp,                       &
                              DD_HStar_old  = 0.0_fp,                       &
                              WD_AerScavEff = 1.0_fp,                       &
                              WD_RainoutEff = RainEff,                      &
                              RC            = RC )

          !==================================================================
          ! Species for the Hg specialty simulation
          !==================================================================

          CASE( 'HG0' )
             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              Name          = 'Hg0',                        &
                              FullName      = 'Elemental mercury',          &
                              MW_g          = 201.0_fp,                     &
                              Is_Advected   = T,                            &
                              Is_Gas        = T,                            &
                              Is_Drydep     = T,                            &
                              Is_Wetdep     = F,                            &
                              DD_F0         = 1.0e-5_fp,                    &
                              DD_Hstar_old  = 1.10e-1_fp,                   &
                              RC            = RC )

          CASE( 'HG2' )
             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              Name          = 'Hg2',                        &
                              FullName      = 'Divalent mercury',           &
                              MW_g          = 201.0_fp,                     &
                              Is_Advected   = T,                            &
                              Is_Gas        = T,                            &
                              Is_Drydep     = T,                            &
                              Is_Wetdep     = T,                            &
                              DD_F0         = 1.0_fp,                       &
#if defined( NEW_HENRY_CONSTANTS )					    
                              Henry_K0      = 1.40e+4_f8 * To_M_atm,        &
                              Henry_CR      = 5300.0_f8,                    &
#else									    
                              DD_Hstar_old  = 1.00e+14_fp,                  &
                              Henry_K0      = 1.40e+6_f8,                   &
                              Henry_CR      = 8400.0_f8,                    &
#endif									    
                              WD_RetFactor  = 0.0_fp,                       &
                              RC            = RC )

          CASE( 'HGP' )

             ! Turn off rainout when 237 K <= T < 258 K
             RainEff = (/ 1.0_fp, 0.0_fp, 1.0_fp /)   

             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              Name          = 'HgP',                        &
                              FullName      = 'Particulate mercury',        &
                              MW_g          = 201.0_fp,                     &
                              Is_Advected   = T,                            &
                              Is_Gas        = F,                            &
                              Is_Drydep     = T,                            &
                              Is_Wetdep     = T,                            &
                              DD_F0         = 0.0_fp,                       &
                              DD_Hstar_old  = 0.0_fp,                       &
                              WD_AerScavEff = 1.0_fp,                       &
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
                              Is_Advected   = T,                            &
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
             !
             ! POPS on hydrophobic BC, which normally has a rainout 
             ! efficiency of zero, is considered to be IN (ice nuclei) and 
             ! therefore is allowed to be scavenged when T < 258 K.
             !
             ! POPS on hydrophobic OC does not rainout.
             SELECT CASE( TRIM( NameAllCaps ) ) 
                CASE( 'POPPBCPO' ) 
                   FullName = TRIM( FullName ) // ' hydrophobic black carbon'
                   RainEff  = (/ 1.0_fp, 1.0_fp, 0.0_fp /)
                CASE( 'POPPOCPO' )
                   FullName = TRIM( FullName ) // ' hydrophobic organic carbon'
                   RainEff  = (/ 0.0_fp, 0.0_fp, 0.0_fp /)
             END SELECT

             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              Name          = NameAllCaps,                  &
                              FullName      = FullName,                     &
                              MW_g          = MW_g,                         &
                              Is_Advected   = T,                            &
                              Is_Gas        = F,                            &
                              Is_Drydep     = T,                            &
                              Is_Wetdep     = T,                            &
                              DD_F0         = 0.0_fp,                       &
                              DD_Hstar_old  = 0.0_fp,                       &
                              WD_AerScavEff = 1.0_fp,                       &
                              WD_RainoutEff = RainEff,                      &
                              RC            = RC )

          CASE( 'POPPBCPI', 'POPPOCPI' )

             ! These have identica properties except for the mol weights ...
             ! Get info from the POP menu
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
             !
             ! Turn off rainout when 237 K <= T < 258 K
             SELECT CASE( TRIM( NameAllCaps ) ) 
                CASE( 'POPPBCPI' ) 
                   FullName = TRIM( FullName ) // ' hydrophilic black carbon'
                   RainEff  = (/ 1.0_fp, 0.0_fp, 1.0_fp /)   
                CASE( 'POPPOCPI' )
                   FullName = TRIM( FullName ) // ' hydrophilic organic carbon'
                   RainEff  = (/ 1.0_fp, 0.0_fp, 1.0_fp /)   
             END SELECT
             
             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              Name          = NameAllCaps,                  &
                              FullName      = FullName,                     &
                              MW_g          = MW_g,                         &
                              Is_Advected   = T,                            &
                              Is_Gas        = F,                            &
                              Is_Drydep     = T,                            &
                              Is_Wetdep     = T,                            &
                              DD_F0         = 0.0_fp,                       &
                              DD_Hstar_old  = 0.0_fp,                       &
                              WD_AerScavEff = 1.0_fp,                       &
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

             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              Name          = Name,                         &
                              FullName      = FullName,                     &
                              MW_g          = 44.0_fp,                      &
                              Is_Advected   = T,                            &
                              Is_Gas        = T,                            &
                              Is_Drydep     = F,                            &
                              Is_Wetdep     = F,                            &
                              RC            = RC )

#if defined( TOMAS )
          !==================================================================
          ! Species for the TOMAS microphysics simulations
          !==================================================================
             
          CASE( 'H2SO4' )
             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              Name          = NameAllCaps,                  &
                              FullName      = 'Sulfuric acid',              &
                              MW_g          = 98.0_fp,                      &
                              MolecRatio    = 1.0_fp,                       &
                              Is_Advected   = T,                            &
                              Is_Gas        = T,                            &
                              Is_Drydep     = T,                            &
                              Is_Wetdep     = F,                            &
                              DD_F0         = 0.0_fp,                       &
                              DD_Hstar_old  = 1.00e+5_fp,                   &
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

             ! Turn off rainout when 237 K <= T < 258 K
             RainEff = (/ 1.0_fp, 0.0_fp, 1.0_fp /)   

             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              Name          = NameAllCaps,                  &
                              FullName      = FullName,                     &
                              MW_g          = 1.0_fp,                       &
                              Is_Advected   = T,                            &
                              Is_Gas        = F,                            &
                              Is_Drydep     = T,                            &
                              Is_Wetdep     = T,                            &
                              DD_F0         = 0.0_fp,                       &
                              DD_Hstar_old  = 0.0_fp,                       &
                              WD_AerScavEff = 1.0_fp,                       &
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

             ! Turn off rainout when 237 K <= T < 258 K
             RainEff = (/ 1.0_fp, 0.0_fp, 1.0_fp /)   

             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              Name          = NameAllCaps,                  &
                              FullName      = FullName,                     &
                              MW_g          = 96.0_fp,                      &
                              Is_Advected   = T,                            &
                              Is_Gas        = F,                            &
                              Is_Drydep     = T,                            &
                              Is_Wetdep     = T,                            &
                              DD_F0         = 0.0_fp,                       &
                              DD_Hstar_old  = 0.0_fp,                       &
                              WD_AerScavEff = 1.0_fp,                       &
                              WD_SizeResAer = T,                            &
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

             ! Turn off rainout when 237 K <= T < 258 K
             RainEff = (/ 1.0_fp, 0.0_fp, 1.0_fp /)   

             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              Name          = NameAllCaps,                  &
                              FullName      = FullName,                     &
                              MW_g          = 58.5_fp,                      &
                              Is_Advected   = T,                            &
                              Is_Gas        = F,                            &
                              Is_Drydep     = T,                            &
                              Is_Wetdep     = T,                            &
                              DD_F0         = 0.0_fp,                       &
                              DD_Hstar_old  = 0.0_fp,                       &
                              WD_AerScavEff = 1.0_fp,                       &
                              WD_SizeResAer = T,                            &
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
             FullName = 'Hydrophobic elemental carbon, size bin ='
             C        = LEN_TRIM( NameAllCaps )
             IF ( C == 5 ) THEN 
                FullName = TRIM( FullName ) // ' ' // NameAllCaps(C:C)
             ELSE
                FullName = TRIM( FullName ) // ' ' // NameAllCaps(C-1:C)
             ENDIF

             ! Turn off rainout when 237 K <= T < 258 K
             RainEff = (/ 1.0_fp, 0.0_fp, 1.0_fp /)   

             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              Name          = NameAllCaps,                  &
                              FullName      = FullName,                     &
                              MW_g          = 12.0_fp,                      &
                              Is_Advected   = T,                            &
                              Is_Gas        = F,                            &
                              Is_Drydep     = T,                            &
                              Is_Wetdep     = T,                            &
                              DD_F0         = 0.0_fp,                       &
                              DD_Hstar_old  = 0.0_fp,                       &
                              WD_AerScavEff = 1.0_fp,                       &
                              WD_SizeResAer = T,                            &
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

             ! Turn off rainout when 237 K <= T < 258 K
             RainEff = (/ 1.0_fp, 0.0_fp, 1.0_fp /)   

             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              Name          = NameAllCaps,                  &
                              FullName      = FullName,                     &
                              MW_g          = 12.0_fp,                      &
                              Is_Advected   = T,                            &
                              Is_Gas        = F,                            &
                              Is_Drydep     = T,                            &
                              Is_Wetdep     = T,                            &
                              DD_F0         = 0.0_fp,                       &
                              DD_Hstar_old  = 0.0_fp,                       &
                              WD_AerScavEff = 1.0_fp,                       &
                              WD_SizeResAer = T,                            &
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
             FullName = 'Hydrophilic organic carbon, size bin ='
             C        = LEN_TRIM( NameAllCaps )
             IF ( C == 5 ) THEN 
                FullName = TRIM( FullName ) // ' ' // NameAllCaps(C:C)
             ELSE
                FullName = TRIM( FullName ) // ' ' // NameAllCaps(C-1:C)
             ENDIF

             ! Turn off rainout when 237 K <= T < 258 K
             !%%% NOTE: Should these be considered "IN" as well?      %%%
             !%%% Ask the TOMAS team for clarification (bmy, 9/24/15) %%%
             RainEff = (/ 1.0_fp, 0.0_fp, 1.0_fp /)   

             CALL Spc_Create( am_I_Root     = am_I_Root,                    &
                              ThisSpc       = SpcData(N)%Info,              &
                              ModelID       = N,                            &
                              Name          = NameAllCaps,                  &
                              FullName      = FullName,                     &
                              MW_g          = 100.0_fp,                     &
                              Is_Advected   = T,                            &
                              Is_Gas        = F,                            &
                              Is_Drydep     = T,                            &
                              Is_Wetdep     = T,                            &
                              DD_F0         = 0.0_fp,                       &
                              DD_Hstar_old  = 0.0_fp,                       &
                              WD_AerScavEff = 1.0_fp,                       &
                              WD_SizeResAer = T,                            &
                              WD_RainoutEff = RainEff,                      &
                              RC            = RC )

#endif

          !==================================================================
          ! Species not found!  Stop with error message
          !==================================================================
          CASE DEFAULT
             WRITE( 6, '(a)' ) REPEAT( '=', 79 )
             WRITE( 6, 100 ) TRIM( NameAllCaps )
             WRITE( 6, 110 )
             WRITE( 6, 120 )
 100         FORMAT( 'Species ', a, ' not found in the data base!'      )
 110         FORMAT( 'Please add the information to the CASE statement' )
 120         FORMAT( 'in module species_database_mod.F90!'              )
             RC = -1
             RETURN

       END SELECT

       ! Error
       IF ( RC /= GIGC_SUCCESS ) THEN
          PRINT*, '### Could not initialize species vector!'
          CALL EXIT( -999 ) 
       ENDIF

!       ! Print info about each species
!       CALL Spc_Print( am_I_Root, SpcData(N)%Info, RC )

    ENDDO

  END SUBROUTINE Init_Species_Database
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: cleanup_define_species
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
    USE GIGC_ErrCode_Mod
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
    RC = GIGC_SUCCESS

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
!  with gigc_state_chm_mod.F90 and species_mod.F90.
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

  END SUBROUTINE TRANUC
!EOC
END MODULE Species_Database_Mod
