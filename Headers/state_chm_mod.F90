!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: state_chm_mod.F90
!
! !DESCRIPTION: Module STATE\_CHM\_MOD contains the derived type
!  used to define the Chemistry State object for GEOS-Chem.
!\\
!\\
!  This module also contains the routines that allocate and deallocate memory
!  to the Chemistry State object.  The chemistry state object is not defined
!  in this module.  It must be be declared as variable in the top-level
!  driver routine, and then passed to lower-level routines as an argument.
!\\
!\\
! !INTERFACE:
!
MODULE State_Chm_Mod
!
! USES:
!
  USE Dictionary_M, ONLY : dictionary_t  ! Fortran hash table type
  USE ErrCode_Mod                        ! Error handling
  USE PhysConstants                      ! Physical constants
  USE Precision_Mod                      ! GEOS-Chem precision types
  USE Registry_Mod                       ! Registry module
  USE Species_Mod                        ! For species database object

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC :: Init_State_Chm
  PUBLIC :: Cleanup_State_Chm
  PUBLIC :: Get_Metadata_State_Chm
  PUBLIC :: Ind_
!
! !PRIVATE MEMBER FUNCTIONS
!
  PRIVATE :: Register_ChmField
!
! !PRIVATE DATA MEMBERS:
!
  TYPE(SpcPtr), PRIVATE, POINTER :: SpcDataLocal(:)  ! Local pointer to
                                                     ! StateChm%SpcData for
                                                     ! availability to IND_

  TYPE(dictionary_t), PRIVATE    :: SpcDictLocal     ! Private copy of the
                                                     ! Fortran Hash table for
                                                     ! availability to IND_


  INTEGER, PRIVATE               :: nChmState = 0    ! # chemistry states,
                                                     ! this CPU
!
! !PUBLIC DATA MEMBERS:
!
  !=========================================================================
  ! Derived type for Chemistry State
  !=========================================================================
  TYPE, PUBLIC :: ChmState

     !----------------------------------------------------------------------
     ! Count of each type of species
     !----------------------------------------------------------------------
     INTEGER                    :: nSpecies             ! # species (all)
     INTEGER                    :: nAdvect              ! # advected species
     INTEGER                    :: nAeroSpc             ! # of Aerosol Species
     INTEGER                    :: nAeroType            ! # of Aerosol Types
     INTEGER                    :: nDryAlt              ! # dryalt species
     INTEGER                    :: nDryDep              ! # drydep species
     INTEGER                    :: nGasSpc              ! # gas phase species
     INTEGER                    :: nHygGrth             ! # hygroscopic growth
     INTEGER                    :: nKppVar              ! # KPP variable species
     INTEGER                    :: nKppFix              ! # KPP fixed species
     INTEGER                    :: nKppSpc              ! # KPP chem species
     INTEGER                    :: nLoss                ! # of loss species
     INTEGER                    :: nPhotol              ! # photolysis species
     INTEGER                    :: nProd                ! # of prod species
     INTEGER                    :: nWetDep              ! # wetdep species

     !----------------------------------------------------------------------
     ! Mapping vectors to subset types of species
     !----------------------------------------------------------------------
     INTEGER,           POINTER :: Map_Advect (:      ) ! Advected species IDs
     INTEGER,           POINTER :: Map_Aero   (:      ) ! Aerosol species IDs
     INTEGER,           POINTER :: Map_DryAlt (:      ) ! Dryalt species IDs
     INTEGER,           POINTER :: Map_DryDep (:      ) ! Drydep species IDs
     INTEGER,           POINTER :: Map_GasSpc (:      ) ! Gas species IDs
     INTEGER,           POINTER :: Map_HygGrth(:      ) ! HygGrth species IDs
     INTEGER,           POINTER :: Map_KppVar (:      ) ! Kpp variable spc IDs
     INTEGER,           POINTER :: Map_KppFix (:      ) ! KPP fixed species IDs
     INTEGER,           POINTER :: Map_KppSpc (:      ) ! KPP chem species IDs
     INTEGER,           POINTER :: Map_Loss   (:      ) ! Loss diag species
     CHARACTER(LEN=36), POINTER :: Name_Loss  (:      ) !  ID's and names
     INTEGER,           POINTER :: Map_Photol (:      ) ! Photolysis species IDs
     INTEGER,           POINTER :: Map_Prod   (:      ) ! Prod diag species
     CHARACTER(LEN=36), POINTER :: Name_Prod  (:      ) !  ID and names
     INTEGER,           POINTER :: Map_WetDep (:      ) ! Wetdep species IDs
     INTEGER,           POINTER :: Map_WL     (:      ) ! Wavelength bins in fjx

#if defined( MODEL_GEOS )
     ! For drydep
     REAL(fp),          POINTER :: DryDepRa2m (:,:    ) ! 2m  aerodynamic resistance
     REAL(fp),          POINTER :: DryDepRa10m(:,:    ) ! 10m aerodynamic resistance
#endif

     !----------------------------------------------------------------------
     ! Physical properties & indices for each species
     !----------------------------------------------------------------------
     TYPE(SpcPtr),      POINTER :: SpcData    (:      ) ! GC Species database
     TYPE(dictionary_t)         :: SpcDict              ! Species dictionary

     !----------------------------------------------------------------------
     ! Chemical species
     !----------------------------------------------------------------------
     REAL(fp),          POINTER :: Species    (:,:,:,:) ! Species concentration
                                                        !  [kg/kg dry air]
     CHARACTER(LEN=20)          :: Spc_Units            ! Species units

     !----------------------------------------------------------------------
     ! Boundary conditions
     !----------------------------------------------------------------------
     REAL(fp),          POINTER :: BoundaryCond(:,:,:,:)! Boundary conditions
                                                        !  [kg/kg dry air]

     !----------------------------------------------------------------------
     ! Aerosol quantities
     !----------------------------------------------------------------------
     REAL(fp),          POINTER :: AeroArea   (:,:,:,:) ! Aerosol Area [cm2/cm3]
     REAL(fp),          POINTER :: AeroRadi   (:,:,:,:) ! Aerosol Radius [cm]
     REAL(fp),          POINTER :: WetAeroArea(:,:,:,:) ! Aerosol Area [cm2/cm3]
     REAL(fp),          POINTER :: WetAeroRadi(:,:,:,:) ! Aerosol Radius [cm]
     REAL(fp),          POINTER :: AeroH2O    (:,:,:,:) ! Aerosol water [cm3/cm3]
     REAL(fp),          POINTER :: GammaN2O5  (:,:,:,:) ! N2O5 aerosol uptake [unitless]
     REAL(fp),          POINTER :: SSAlk      (:,:,:,:) ! Sea-salt alkalinity[-]
     REAL(fp),          POINTER :: H2O2AfterChem(:,:,:) ! H2O2, SO2 [v/v]
     REAL(fp),          POINTER :: SO2AfterChem (:,:,:) !  after sulfate chem
     REAL(fp),          POINTER :: OMOC_POA       (:,:) ! OM:OC Ratio (OCFPOA) [unitless]
     REAL(fp),          POINTER :: OMOC_OPOA      (:,:) ! OM:OC Ratio (OCFOPOA) [unitless]

     !----------------------------------------------------------------------
     ! Fields for nitrogen deposition
     !----------------------------------------------------------------------
     REAL(fp),          POINTER :: DryDepNitrogen (:,:) ! Dry deposited N
     REAL(fp),          POINTER :: WetDepNitrogen (:,:) ! Wet deposited N

     !----------------------------------------------------------------------
     ! Cloud quantities
     !----------------------------------------------------------------------
     REAL(fp),          POINTER :: pHCloud    (:,:,:  ) ! Cloud pH [-]
     REAL(fp),          POINTER :: isCloud    (:,:,:  ) ! Cloud presence [-]

     !----------------------------------------------------------------------
     ! Fields for KPP solver
     !----------------------------------------------------------------------
     REAL(fp),          POINTER :: KPPHvalue  (:,:,:  ) ! H-value for Rosenbrock
                                                        !  solver
     !----------------------------------------------------------------------
     ! Fields for UCX mechanism
     !----------------------------------------------------------------------
     REAL(f4),          POINTER :: STATE_PSC  (:,:,:  ) ! PSC type (see Kirner
                                                        !  et al. 2011, GMD)
     REAL(fp),          POINTER :: KHETI_SLA  (:,:,:,:) ! Strat. liquid aerosol
                                                        !  reaction cofactors

     !----------------------------------------------------------------------
     ! For isoprene SOA via ISORROPIA
     !----------------------------------------------------------------------
     REAL(fp),          POINTER :: pHSav      (:,:,:  ) ! ISORROPIA aerosol pH
     REAL(fp),          POINTER :: HplusSav   (:,:,:  ) ! H+ concentration [M]
     REAL(fp),          POINTER :: WaterSav   (:,:,:  ) ! ISORROPIA aerosol H2O
     REAL(fp),          POINTER :: SulRatSav  (:,:,:  ) ! Sulfate conc [M]
     REAL(fp),          POINTER :: NaRatSav   (:,:,:  ) ! Nitrate conc [M]
     REAL(fp),          POINTER :: AcidPurSav (:,:,:  ) !
     REAL(fp),          POINTER :: BiSulSav   (:,:,:  ) ! Bisulfate conc [M]

     !----------------------------------------------------------------------
     ! For the tagged Hg simulation
     !----------------------------------------------------------------------
     INTEGER                    :: N_HG_CATS            ! # of Hg categories
     INTEGER,           POINTER :: Hg0_Id_List(:      ) ! Hg0 cat <-> tracer #
     INTEGER,           POINTER :: Hg2_Id_List(:      ) ! Hg2 cat <-> tracer #
     INTEGER,           POINTER :: HgP_Id_List(:      ) ! HgP cat <-> tracer #
     CHARACTER(LEN=4),  POINTER :: Hg_Cat_Name(:      ) ! Category names

     REAL(fp),          POINTER :: OceanHg0(:,:,:)      ! Hg(0)  ocean mass [kg]
     REAL(fp),          POINTER :: OceanHg2(:,:,:)      ! Hg(II) ocean mass [kg]
     REAL(fp),          POINTER :: OceanHgP(:,:,:)      ! HgP    ocean mass [kg]
     REAL(fp),          POINTER :: SnowHgOcean(:,:,:)   ! Reducible Hg snowpack
                                                        !  on ocean [kg]
     REAL(fp),          POINTER :: SnowHgLand(:,:,:)    ! Reducible Hg snowpack
                                                        !  on land [kg]
     REAL(fp),          POINTER :: SnowHgOceanStored(:,:,:) ! Non-reducible Hg
                                                            !  snowpack on ocean
     REAL(fp),          POINTER :: SnowHgLandStored(:,:,:)  ! Non-reducible Hg
                                                            !  snowpack on land

     !----------------------------------------------------------------------
     ! For HOBr + S(IV) heterogeneous chemistry
     !----------------------------------------------------------------------
     REAL(fp),          POINTER :: HSO3_AQ    (:,:,:  ) ! Cloud bisulfite[mol/l]
     REAL(fp),          POINTER :: SO3_AQ     (:,:,:  ) ! Cloud sulfite  [mol/l]
     REAL(fp),          POINTER :: fupdateHOBr(:,:,:  ) ! Correction factor for
                                                        ! HOBr removal by SO2
                                                        ! [unitless]

     !----------------------------------------------------------------------
     ! Fields for dry deposition
     !----------------------------------------------------------------------
     REAL(fp),          POINTER :: DryDepSav  (:,:,:  ) ! Drydep freq [s-1]

     !----------------------------------------------------------------------
     ! Fields for Linoz stratospheric ozone algorithm
     !----------------------------------------------------------------------
     REAL(fp),          POINTER :: TLSTT      (:,:,:,:) ! TLSTT (I,J,L,LINOZ_NFIELDS)

     !----------------------------------------------------------------------
     ! Fields for Gan Luo et al Wetdep scheme (GMD-12-3439-2019)
     !----------------------------------------------------------------------
     REAL(fp),          POINTER :: PSO4s      (:,:,:  )
     REAL(fp),          POINTER :: QQ3D       (:,:,:  )

     !----------------------------------------------------------------------
     ! Fields for non-local PBL mixing
     !----------------------------------------------------------------------
     REAL(fp),          POINTER :: SurfaceFlux(:,:,:  )

     !----------------------------------------------------------------------
     ! Registry of variables contained within State_Chm
     !----------------------------------------------------------------------
     CHARACTER(LEN=4)           :: State     = 'CHEM'   ! Name of this state
     TYPE(MetaRegItem), POINTER :: Registry  => NULL()  ! Registry object
     TYPE(dictionary_t)         :: RegDict              ! Registry lookup table

  END TYPE ChmState
!
! !REMARKS:
!
! !REVISION HISTORY:
!  19 Oct 2012 - R. Yantosca - Initial version, based on "gc_type2_mod.F90"
!  See https://github.com/geoschem/geos-chem for complete history
!------------------------------------------------------------------------------
!BOC
!
! !MODULE INTERFACES:
!
  INTERFACE Register_ChmField
     MODULE PROCEDURE Register_ChmField_R4_3D
     MODULE PROCEDURE Register_ChmField_Rfp_2D
     MODULE PROCEDURE Register_ChmField_Rfp_3D
     MODULE PROCEDURE Register_ChmField_Rfp_4D
  END INTERFACE Register_ChmField

CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Init_State_Chm
!
! !DESCRIPTION: Routine INIT\_STATE\_CHM allocates and initializes the
!  pointer fields of the chemistry state object.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Init_State_Chm( Input_Opt, State_Chm, State_Grid, RC )
!
! !USES:
!
    USE CharPak_Mod,          ONLY : To_UpperCase
    USE CMN_Size_Mod,         ONLY : NDUST, NAER
    USE GCKPP_Parameters,     ONLY : NSPEC
    USE Input_Opt_Mod,        ONLY : OptInput
    USE Species_Database_Mod, ONLY : Init_Species_Database
    USE State_Grid_Mod,       ONLY : GrdState
    USE CMN_FJX_MOD,          ONLY : W_         ! For UVFlx diagnostic
!
! !INPUT PARAMETERS:
!
    TYPE(GrdState), INTENT(IN)    :: State_Grid  ! Grid State object
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(INOUT) :: Input_Opt   ! Input Options object
    TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Return code
!
! !REMARKS:
!  In the near future we will put some error trapping on the allocations
!  so that we can stop the simulation if the allocations cannot be made.
!
! !REVISION HISTORY:
!  19 Oct 2012 - R. Yantosca - Renamed from gc_type2_mod.F90
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER                :: N, C, IM, JM, LM
    INTEGER                :: N_Hg0_CATS, N_Hg2_CATS, N_HgP_CATS
    INTEGER                :: nKHLSA, nAerosol, nMatches

    ! Strings
    CHARACTER(LEN=255)     :: ErrMsg, ThisLoc, ChmID

    ! Pointers
    TYPE(Species), POINTER :: ThisSpc
    INTEGER,       POINTER :: CheckIds(:)
    REAL(fp),      POINTER :: Ptr2data(:,:,:)

    ! Error handling
    RC      = GC_SUCCESS
    ErrMsg  = ''
    ThisLoc = ' -> at Init_State_Chm (in Headers/state_chm_mod.F90)'

    !=======================================================================
    ! Initialization
    !=======================================================================

    ! Count the # of chemistry states we have initialized, so SpcData(Local)
    ! is not deallocated until the last ChmState is cleaned up.
    ! This avoids dangling pointers with detrimental effects. (hplin, 8/3/18)
    nChmState                   =  nChmState + 1

    ! Shorten grid parameters for readability
    IM                          =  State_Grid%NX ! # latitudes
    JM                          =  State_Grid%NY ! # longitudes
    LM                          =  State_Grid%NZ ! # levels

    ! Number of aerosols
    nAerosol                    =  NDUST + NAER

    ! Number of each type of species
    State_Chm%nSpecies          =  0
    State_Chm%nAdvect           =  0
    State_Chm%nAeroSpc          =  0
    State_Chm%nAeroType         =  0
    State_Chm%nDryAlt           =  0
    State_Chm%nDryDep           =  0
    State_Chm%nGasSpc           =  0
    State_Chm%nHygGrth          =  0
    State_Chm%nKppVar           =  0
    State_Chm%nKppFix           =  0
    State_Chm%nKppSpc           =  0
    State_Chm%nLoss             =  0
    State_Chm%nPhotol           =  0
    State_Chm%nProd             =  0
    State_Chm%nWetDep           =  0


    ! Mapping vectors for subsetting each type of species
    State_Chm%Map_Advect        => NULL()
    State_Chm%Map_Aero          => NULL()
    State_Chm%Map_DryAlt        => NULL()
    State_Chm%Map_DryDep        => NULL()
    State_Chm%Map_GasSpc        => NULL()
    State_Chm%Map_HygGrth       => NULL()
    State_Chm%Map_KppVar        => NULL()
    State_Chm%Map_KppFix        => NULL()
    State_Chm%Map_KppSpc        => NULL()
    State_Chm%Name_Loss         => NULL()
    State_Chm%Map_Loss          => NULL()
    State_Chm%Map_Photol        => NULL()
    State_Chm%Name_Prod         => NULL()
    State_Chm%Map_Prod          => NULL()
    State_Chm%Map_WetDep        => NULL()
    State_Chm%Map_WL            => NULL()

    ! Chemical species
    State_Chm%Species           => NULL()
    State_Chm%Spc_Units         = ''

    ! Boundary conditions
    State_Chm%BoundaryCond      => NULL()

    ! Species database
    State_Chm%SpcData           => NULL()
    ThisSpc                     => NULL()

    ! Aerosol parameters
    State_Chm%AeroArea          => NULL()
    State_Chm%AeroRadi          => NULL()
    State_Chm%WetAeroArea       => NULL()
    State_Chm%WetAeroRadi       => NULL()
    State_Chm%AeroH2O           => NULL()
    State_Chm%GammaN2O5         => NULL()
    State_Chm%OMOC_POA          => NULL()
    State_Chm%OMOC_OPOA         => NULL()

    ! Isoprene SOA
    State_Chm%pHSav             => NULL()
    State_Chm%HplusSav          => NULL()
    State_Chm%WaterSav          => NULL()
    State_Chm%SulRatSav         => NULL()
    State_Chm%NaRatSav          => NULL()
    State_Chm%AcidPurSav        => NULL()
    State_Chm%BisulSav          => NULL()

    ! Fields for KPP solver
    State_Chm%KPPHvalue         => NULL()

    ! Fields for UCX mechanism
    State_Chm%STATE_PSC         => NULL()
    State_Chm%KHETI_SLA         => NULL()

    ! pH/alkalinity
    State_Chm%pHCloud           => NULL()
    State_Chm%isCloud           => NULL()
    State_Chm%SSAlk             => NULL()

    ! Fields for sulfate chemistry
    State_Chm%H2O2AfterChem     => NULL()
    State_Chm%SO2AfterChem      => NULL()

    ! Fields for nitrogen deposition
    State_Chm%DryDepNitrogen    => NULL()
    State_Chm%WetDepNitrogen    => NULL()

    ! Hg species indexing
    N_Hg0_CATS                  =  0
    N_Hg2_CATS                  =  0
    N_HgP_CATS                  =  0
    State_Chm%N_Hg_CATS         =  0
    State_Chm%Hg_Cat_Name       => NULL()
    State_Chm%Hg0_Id_List       => NULL()
    State_Chm%Hg2_Id_List       => NULL()
    State_Chm%HgP_Id_List       => NULL()
    State_Chm%OceanHg0          => NULL()
    State_Chm%OceanHg2          => NULL()
    State_Chm%OceanHgP          => NULL()
    State_Chm%SnowHgOcean       => NULL()
    State_Chm%SnowHgLand        => NULL()
    State_Chm%SnowHgOceanStored => NULL()
    State_Chm%SnowHgLandStored  => NULL()

    ! For HOBr + S(IV) chemistry
    State_Chm%HSO3_AQ           => NULL()
    State_Chm%SO3_AQ            => NULL()
    State_Chm%fupdateHOBr       => NULL()

    ! For Luo et al wetdep
    State_Chm%PSO4s             => NULL()
    State_Chm%QQ3D              => NULL()

    ! For LINOZ
    State_Chm%TLSTT             => NULL()

    State_Chm%DryDepSav         => NULL()
    State_Chm%SurfaceFlux       => NULL()

    ! Local variables
    Ptr2data                    => NULL()

    !=======================================================================
    ! Populate the species database object field
    ! (assumes Input_Opt has already been initialized)
    !=======================================================================

    ! If the species database has already been initialized in this CPU,
    ! SpcDataLocal in State_Chm_Mod already contains a copy of the species data.
    ! It can be directly associated to this new chemistry state.
    ! (assumes one CPU will run one copy of G-C with the same species DB)
    IF ( ASSOCIATED( SpcDataLocal ) ) THEN
        State_Chm%SpcData => SpcDataLocal
    ELSE
        CALL Init_Species_Database( Input_Opt = Input_Opt,                   &
                                    SpcData   = State_Chm%SpcData,           &
                                    RC        = RC                           )

        ! Point to a private module copy of the species database
        ! which will be used by the Ind_ indexing function
        SpcDataLocal => State_Chm%SpcData
    ENDIF

    !=======================================================================
    ! Before proceeding, make sure none of the species has a blank name,
    ! because this has the potential to halt the run inadvertently.
    !=======================================================================

    ! The total number of species is the size of SpcData
    State_Chm%nSpecies = SIZE( State_Chm%SpcData )

    ! Exit if any species name is blank
    DO N = 1, State_Chm%nSpecies
       IF ( LEN_TRIM(  State_Chm%SpcData(N)%Info%Name ) == 0 ) THEN
          WRITE( ErrMsg, '("Species number ", i6, " has a blank name!")' ) N
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
    ENDDO

    !=======================================================================
    ! Determine the number of advected, drydep, wetdep, and total species
    !=======================================================================

    ! Get the number of advected, dry-deposited, KPP chemical species,
    ! and and wet-deposited species.  Also return the # of Hg0, Hg2, and
    ! HgP species (these are zero unless the Hg simulation is used).
    CALL Spc_GetNumSpecies( nAdvect  = State_Chm%nAdvect,                  &
                            nAeroSpc = State_Chm%nAeroSpc,                 &
                            nDryAlt  = State_Chm%nDryAlt,                  &
                            nDryDep  = State_Chm%nDryDep,                  &
                            nGasSpc  = State_Chm%nGasSpc,                  &
                            nHygGrth = State_Chm%nHygGrth,                 &
                            nKppVar  = State_Chm%nKppVar,                  &
                            nKppFix  = State_Chm%nKppFix,                  &
                            nKppSpc  = State_Chm%nKppSpc,                  &
                            nPhotol  = State_Chm%nPhotol,                  &
                            nWetDep  = State_Chm%nWetDep,                  &
                            nHg0Cats = N_Hg0_CATS,                         &
                            nHg2Cats = N_Hg2_CATS,                         &
                            nHgPCats = N_HgP_CATS                         )

    ! Also get the number of the prod/loss species.  For fullchem simulations,
    ! the prod/loss species are listed in FAM_NAMES in gckpp_Monitor.F90,
    ! but for certain other simulations (tagO3, tagCO), advected species
    ! can have prod and loss diagnostic entries.
    CALL GetNumProdLossSpecies( Input_Opt, State_Chm, RC )
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Error encountered in "GetNumProdLossSpecies"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    !########################################################################
    !### Save species database info to a HEMCO_sa_Spec.rc file for use with
    !### the HEMCO standalone simulation.  Uncomment this if you need it.
    !### (bmy, 9/26/18)
    !###
    !
    !! Open file
    !OPEN( 700, FILE = 'HEMCO_sa_Spec.rc', STATUS = 'UNKNOWN', IOSTAT=RC )
    !
    !! Write data
    !DO N = 1, State_Chm%nAdvect
    !   WRITE( 700, 700 ) N, State_Chm%SpcData(N)%Info%Name,                  &
    !                        State_Chm%SpcData(N)%Info%Mw_g,                  &
    !                        State_Chm%SpcData(N)%Info%EmMw_g,                &
    !                        State_Chm%SpcData(N)%Info%MolecRatio,            &
    !                        MAX(State_Chm%SpcData(N)%Info%Henry_K0, 0.0_fp), &
    !                        MAX(State_Chm%SpcData(N)%Info%Henry_CR, 0.0_fp), &
    !                        MAX(State_Chm%SpcData(N)%Info%Henry_pKa,0.0_fp)
    !
    !   700 FORMAT( i4, 1x, a10, 1x, 2f9.2, f5.1, 2x, es13.6, 2f10.2 )
    !ENDDO
    !
    !! Close file
    !CLOSE( 700 )
    !STOP
    !########################################################################

    !=======================================================================
    ! Populate the species lookup table, for quick index lookup via Ind_
    !=======================================================================

    ! Initialize the species lookup table
    CALL State_Chm%SpcDict%Init( State_Chm%nSpecies )

    ! Populate the species lookup table
    DO N = 1, State_Chm%nSpecies
       ThisSpc => SpcDataLocal(N)%Info
       CALL State_Chm%SpcDict%Set( To_UpperCase( TRIM( ThisSpc%Name ) ),     &
                                   ThisSpc%ModelId                          )
       ThisSpc => NULL()
    ENDDO

    ! Error check: make sure we have no hash collisions that would
    ! assign more than one species to the same ModelId value
    ALLOCATE( CheckIds( State_Chm%nSpecies ), STAT=RC )
    DO N = 1, State_Chm%nSpecies
       CheckIds(N) = SpcDataLocal(N)%Info%ModelId
    ENDDO
    DO N = 1, State_Chm%nSpecies
       nMatches = COUNT( CheckIds(N) == CheckIds )
       IF ( nMatches > 1 ) THEN
          ErrMsg = 'Species: ' // TRIM( SpcDataLocal(N)%Info%Name )       // &
                   'maps to more than one ModelID value!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          CheckIds => NULL()
          RETURN
       ENDIF
    ENDDO
    IF ( ASSOCIATED( CheckIds ) ) DEALLOCATE( CheckIds )

    ! If there are no hash collisions, then species lookup table
    ! to a local shadow variable for use with the Ind_ function.
    SpcDictLocal = State_Chm%SpcDict

    !### Debug: Show the values in the lookup table
    !###CALL State_Chm%SpcDict%Show()

    !=======================================================================
    ! Exit if this is a dry-run simulation
    !=======================================================================
    IF ( Input_Opt%DryRun ) THEN
       RC = GC_SUCCESS
       RETURN
    ENDIF

    !=======================================================================
    ! Allocate and initialize mapping vectors to subset species
    !=======================================================================

    IF ( State_Chm%nAdvect > 0 ) THEN
       ALLOCATE( State_Chm%Map_Advect( State_Chm%nAdvect ), STAT=RC )
       CALL GC_CheckVar( 'State_Chm%Map_Advect', 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Chm%Map_Advect = 0
    ELSE
       ErrMsg = 'No advected species specified!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    IF ( State_Chm%nAeroSpc > 0 ) THEN
       ALLOCATE( State_Chm%Map_Aero( State_Chm%nAeroSpc ), STAT=RC )
       CALL GC_CheckVar( 'State_Chm%Map_Aero', 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Chm%Map_Aero = 0
    ENDIF

    IF (  State_Chm%nDryAlt > 0 ) THEN
       ALLOCATE( State_Chm%Map_DryAlt( State_Chm%nDryAlt ), STAT=RC )
       CALL GC_CheckVar( 'State_Chm%Map_DryAlt', 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Chm%Map_DryAlt = 0
    ENDIF

    IF (  State_Chm%nDryDep > 0 ) THEN
       ALLOCATE( State_Chm%Map_Drydep( State_Chm%nDryDep ), STAT=RC )
       CALL GC_CheckVar( 'State_Chm%Map_Drydep', 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Chm%Map_DryDep = 0
    ENDIF

    IF ( State_Chm%nGasSpc > 0 ) THEN
       ALLOCATE( State_Chm%Map_GasSpc( State_Chm%nGasSpc ), STAT=RC )
       CALL GC_CheckVar( 'State_Chm%Map_GasSpc', 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Chm%Map_GasSpc = 0
    ENDIF

    IF ( State_Chm%nHygGrth > 0 ) THEN
       ALLOCATE( State_Chm%Map_HygGrth( State_Chm%nHygGrth ), STAT=RC )
       CALL GC_CheckVar( 'State_Chm%Map_HygGrth', 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Chm%Map_HygGrth = 0
    ENDIF

    IF ( State_Chm%nKppVar > 0 ) THEN
       ALLOCATE( State_Chm%Map_KppVar( State_Chm%nKppVar ), STAT=RC )
       CALL GC_CheckVar( 'State_Chm%Map_KppVar', 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Chm%Map_KppVar = 0
    ENDIF

    IF ( State_Chm%nKppFix > 0 ) THEN
       ALLOCATE( State_Chm%Map_KppFix( State_Chm%nKppFix ), STAT=RC )
       CALL GC_CheckVar( 'State_Chm%Map_KppFix', 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Chm%Map_KppFix = 0
    ENDIF

    IF ( NSPEC > 0 ) THEN
       ALLOCATE( State_Chm%Map_KppSpc( NSPEC ), STAT=RC )
       CALL GC_CheckVar( 'State_Chm%Map_KppSpc', 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Chm%Map_KppSpc = 0
    ENDIF

    IF ( State_Chm%nLoss > 0 ) THEN
       ALLOCATE( State_Chm%Name_Loss( State_Chm%nLoss ), STAT=RC )
       CALL GC_CheckVar( 'State_Chm%Name_Loss', 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Chm%Name_Loss = ''

       ALLOCATE( State_Chm%Map_Loss( State_Chm%nLoss ), STAT=RC )
       CALL GC_CheckVar( 'State_Chm%Map_Loss', 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Chm%Map_Loss = 0
    ENDIF

    IF ( State_Chm%nPhotol > 0 ) THEN
       ALLOCATE( State_Chm%Map_Photol( State_Chm%nPhotol ), STAT=RC )
       CALL GC_CheckVar( 'State_Chm%Map_Photol', 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Chm%Map_Photol = 0
    ENDIF

    IF ( State_Chm%nProd >0 ) THEN
       ALLOCATE( State_Chm%Name_Prod( State_Chm%nProd ), STAT=RC )
       CALL GC_CheckVar( 'State_Chm%Name_Prod', 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Chm%Name_Prod = ''

       ALLOCATE( State_Chm%Map_Prod( State_Chm%nProd ), STAT=RC )
       CALL GC_CheckVar( 'State_Chm%Map_Prod', 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Chm%Map_Prod = 0
    ENDIF

    IF ( State_Chm%nWetDep > 0 ) THEN
       ALLOCATE( State_Chm%Map_WetDep( State_Chm%nWetDep ), STAT=RC )
       CALL GC_CheckVar( 'State_Chm%Map_Wetdep', 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Chm%Map_WetDep = 0
    ENDIF

    IF ( W_ > 0 ) THEN
       ALLOCATE( State_Chm%Map_WL( W_ ), STAT=RC )
       CALL GC_CheckVar( 'State_Chm%Map_WL', 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Chm%Map_WL = 0
    ENDIF

    !=======================================================================
    ! Set up the species mapping vectors
    !=======================================================================
    IF ( Input_Opt%amIRoot ) THEN
       WRITE( 6,'(/,a)' ) 'ADVECTED SPECIES MENU'
       WRITE( 6,'(  a)' ) REPEAT( '-', 48 )
       WRITE( 6,'(  a)' ) '  #  Species Name'
    ENDIF

    ! Loop over all species
    DO N = 1, State_Chm%nSpecies

       ! GEOS-Chem Species Database entry for species # N
       ThisSpc => State_Chm%SpcData(N)%Info

       !--------------------------------------------------------------------
       ! Set up the mapping for ADVECTED SPECIES
       !--------------------------------------------------------------------
       IF ( ThisSpc%Is_Advected ) THEN

          ! Update the mapping vector of advected species
          C                       = ThisSpc%AdvectId
          State_Chm%Map_Advect(C) = ThisSpc%ModelId

          ! Print to screen
          IF ( Input_Opt%amIRoot ) THEN
             WRITE( 6, 100 ) ThisSpc%ModelId, ThisSpc%Name
          ENDIF

       ENDIF

       !--------------------------------------------------------------------
       ! Set up the mapping for AEROSOL SPECIES
       !--------------------------------------------------------------------
       IF ( ThisSpc%Is_Aero ) THEN
          C                     = ThisSpc%AeroId
          State_Chm%Map_Aero(C) = ThisSpc%ModelId
       ENDIF

       !--------------------------------------------------------------------
       ! Set up the mapping for DRYDEP SPECIES TO SAVE AT A GIVEN ALTITUDE
       !--------------------------------------------------------------------
       IF ( ThisSpc%Is_DryAlt ) THEN
          C                       = ThisSpc%DryAltId
          State_Chm%Map_DryAlt(C) = ThisSpc%ModelId
       ENDIF

       !--------------------------------------------------------------------
       ! Set up the mapping for DRYDEP SPECIES
       !--------------------------------------------------------------------
       IF ( ThisSpc%Is_DryDep ) THEN
          C                       = ThisSpc%DryDepId
          State_Chm%Map_Drydep(C) = ThisSpc%ModelId
       ENDIF

       !--------------------------------------------------------------------
       ! Set up the mapping for GAS SPECIES
       !--------------------------------------------------------------------
       IF ( ThisSpc%Is_Gas ) THEN
          C                       = ThisSpc%GasSpcId
          State_Chm%Map_GasSpc(C) = ThisSpc%ModelId
       ENDIF

       !--------------------------------------------------------------------
       ! Set up the mapping for HYGROSCOPIC GROWTH SPECIES
       !--------------------------------------------------------------------
       IF ( ThisSpc%Is_HygroGrowth ) THEN
          C                        = ThisSpc%HygGrthId
          State_Chm%Map_HygGrth(C) = ThisSpc%ModelId
       ENDIF

       !--------------------------------------------------------------------
       ! Set up the mapping for KPP ACTIVE (VARIABLE) SPECIES
       !--------------------------------------------------------------------
       IF ( ThisSpc%Is_ActiveChem ) THEN
          C                       = ThisSpc%KppVarId
          State_Chm%Map_KppVar(C) = ThisSpc%ModelId
       ENDIF

       !--------------------------------------------------------------------
       ! Set up the mapping for KPP FIXED SPECIES
       !--------------------------------------------------------------------
       IF ( ThisSpc%Is_FixedChem ) THEN
          C                       = ThisSpc%KppFixId
          State_Chm%Map_KppFix(C) = ThisSpc%ModelId
       ENDIF

       !--------------------------------------------------------------------
       ! Set up the mapping for SPECIES IN THE KPP MECHANISM
       !--------------------------------------------------------------------
       IF ( ThisSpc%Is_Kpp ) THEN
          C                       = ThisSpc%KppSpcId
          State_Chm%Map_KppSpc(C) = ThisSpc%ModelId
       ENDIF

       !--------------------------------------------------------------------
       ! Set up the mapping for PHOTOLYSIS SPECIES
       !--------------------------------------------------------------------
       IF ( ThisSpc%Is_Photolysis ) THEN
          C                       = ThisSpc%PhotolId
          State_Chm%Map_Photol(C) = ThisSpc%ModelId
       ENDIF

       !--------------------------------------------------------------------
       ! Set up the mapping for WETDEP SPECIES
       !--------------------------------------------------------------------
       IF ( ThisSpc%Is_WetDep ) THEN
          C                       = ThisSpc%WetDepId
          State_Chm%Map_WetDep(C) = ThisSpc%ModelId
       ENDIF

       ! Free pointer
       ThisSpc => NULL()

    ENDDO

    !-----------------------------------------------------------------------
    ! Set up the mapping for UVFlux Diagnostics
    ! placeholder for now since couldn't figure out how to read in WL from file
    !-----------------------------------------------------------------------
    IF ( W_ > 0 ) THEN

       DO N = 1, W_
          !
          ! Define identifying string
                State_Chm%Map_WL(N) = 0
       ENDDO
    ENDIF

    !-----------------------------------------------------------------------
    ! Set up the mapping for PRODUCTION AND LOSS DIAGNOSTIC SPECIES
    !-----------------------------------------------------------------------
    IF ( State_Chm%nProd > 0 .or. State_Chm%nLoss > 0 ) THEN
       CALL MapProdLossSpecies( Input_Opt, State_Chm, RC )
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Error encountered in "MapProdLossSpecies"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
    ENDIF

    !=======================================================================
    ! Allocate and initialize chemical species fields
    !=======================================================================
    chmID = 'Species'
    ALLOCATE( State_Chm%Species( IM, JM, LM, State_Chm%nSpecies ), STAT=RC )
    CALL GC_CheckVar( 'State_Chm%Species', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Chm%Species = 0.0_fp
    CALL Register_ChmField( Input_Opt, chmID, State_Chm%Species, State_Chm, RC )
    CALL GC_CheckVar( 'State_Chm%Species', 1, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    !=======================================================================
    ! Allocate and initialize boundary condition fields
    !=======================================================================
    chmID = 'BoundaryCond'
    ALLOCATE( State_Chm%BoundaryCond( IM, JM, LM, State_Chm%nSpecies ), STAT=RC)
    CALL GC_CheckVar( 'State_Chm%BoundaryCond', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Chm%BoundaryCond = 0.0_fp
    CALL Register_ChmField( Input_Opt, chmID, State_Chm%BoundaryCond, &
                            State_Chm, RC )
    CALL GC_CheckVar( 'State_Chm%BoundaryCond', 1, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

#if defined( MODEL_GEOS )
    !=======================================================================
    ! Allocate and initialize aerodynamic resistance fields
    !=======================================================================
    ALLOCATE( State_Chm%DryDepRa2m( IM, JM ), STAT=RC )
    CALL GC_CheckVar( 'State_Chm%DryDepRa2m', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Chm%DryDepRa2m = 0.0_fp

    ALLOCATE( State_Chm%DryDepRa10m( IM, JM ), STAT=RC )
    CALL GC_CheckVar( 'State_Chm%DryDepRa10m', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Chm%DryDepRa10m = 0.0_fp
#endif

    !=======================================================================
    ! Allocate and initialize quantities that are only relevant for the
    ! the various fullchem simulations or the aerosol-only simulation
    !=======================================================================
    IF ( Input_Opt%ITS_A_FULLCHEM_SIM .or. Input_Opt%ITS_AN_AEROSOL_SIM ) THEN

       ! Save nAerosol to State_Chm
       State_Chm%nAeroType = nAerosol

       !--------------------------------------------------------------------
       ! AeroArea
       !--------------------------------------------------------------------
       ALLOCATE( State_Chm%AeroArea(IM,JM,LM,State_Chm%nAeroType), STAT=RC )
       CALL GC_CheckVar( 'State_Chm%AeroArea', 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Chm%AeroArea = 0.0_fp

       ! Loop over all entries to register each category individually
       DO N = 1, State_Chm%nAeroType

          ! Define identifying string
          SELECT CASE( N )
             CASE( 1  )
                chmID = 'AeroAreaMDUST1'
             CASE( 2  )
                chmID = 'AeroAreaMDUST2'
             CASE( 3  )
                chmID = 'AeroAreaMDUST3'
             CASE( 4  )
                chmID = 'AeroAreaMDUST4'
             CASE( 5  )
                chmID = 'AeroAreaMDUST5'
             CASE( 6  )
                chmID = 'AeroAreaMDUST6'
             CASE( 7  )
                chmID = 'AeroAreaMDUST7'
             CASE( 8  )
                chmID = 'AeroAreaSULF'
             CASE( 9  )
                chmID = 'AeroAreaBC'
             CASE( 10 )
                chmID = 'AeroAreaOC'
             CASE( 11 )
                chmID = 'AeroAreaSSA'
             CASE( 12 )
                chmID = 'AeroAreaSSC'
             CASE( 13 )
                chmID = 'AeroAreaBGSULF'
             CASE( 14 )
                chmID  = 'AeroAreaICEI'
             CASE DEFAULT
                ErrMsg = 'State_Chm%nAeroType exceeds the number of defined' &
                         // ' dry aerosol area categories'
                CALL GC_Error( ErrMsg, RC, ThisLoc )
                RETURN
          END SELECT

          CALL Register_ChmField( Input_Opt, chmID, State_Chm%AeroArea,     &
                                  State_Chm, RC,    Ncat=N )
          CALL GC_CheckVar( 'State_Chm%AeroArea', 1, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDDO

       !--------------------------------------------------------------------
       ! AeroRadi
       !--------------------------------------------------------------------
       ALLOCATE( State_Chm%AeroRadi(IM,JM,LM,State_Chm%nAeroType), STAT=RC )
       CALL GC_CheckVar( 'State_Chm%AeroRadi', 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Chm%AeroRadi    = 0.0_fp

       ! Loop over all entries to register each category individually
       DO N = 1, State_Chm%nAeroType

          ! Define identifying string
          SELECT CASE( N )
             CASE( 1  )
                chmID = 'AeroRadiMDUST1'
             CASE( 2  )
                chmID = 'AeroRadiMDUST2'
             CASE( 3  )
                chmID = 'AeroRadiMDUST3'
             CASE( 4  )
                chmID = 'AeroRadiMDUST4'
             CASE( 5  )
                chmID = 'AeroRadiMDUST5'
             CASE( 6  )
                chmID = 'AeroRadiMDUST6'
             CASE( 7  )
                chmID = 'AeroRadiMDUST7'
             CASE( 8  )
                chmID = 'AeroRadiSULF'
             CASE( 9  )
                chmID = 'AeroRadiBC'
             CASE( 10 )
                chmID = 'AeroRadiOC'
             CASE( 11 )
                chmID = 'AeroRadiSSA'
             CASE( 12 )
                chmID = 'AeroRadiSSC'
             CASE( 13 )
                chmID = 'AeroRadiBGSULF'
             CASE( 14 )
                chmID = 'AeroRadiICEI'
             CASE DEFAULT
                ErrMsg = 'State_Chm%nAeroType exceeds the number of defined' &
                         // ' dry aerosol radius categories'
                CALL GC_Error( ErrMsg, RC, ThisLoc )
                RETURN
          END SELECT

          CALL Register_ChmField( Input_Opt, chmID, State_Chm%AeroRadi,      &
                                  State_Chm, RC,    Ncat=N                  )
          CALL GC_CheckVar( 'State_Chm%AeroRadi', 1, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDDO

       !--------------------------------------------------------------------
       ! WetAeroArea
       !--------------------------------------------------------------------
       ALLOCATE( State_Chm%WetAeroArea(IM,JM,LM,State_Chm%nAeroType), STAT=RC )
       CALL GC_CheckVar( 'State_Chm%WetAeroArea', 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Chm%WetAeroArea = 0.0_fp

       ! Loop over all entries to register each category individually
       DO N = 1, State_Chm%nAeroType

          ! Define identifying string
          SELECT CASE( N )
             CASE( 1  )
                chmID = 'WetAeroAreaMDUST1'
             CASE( 2  )
                chmID = 'WetAeroAreaMDUST2'
             CASE( 3  )
                chmID = 'WetAeroAreaMDUST3'
             CASE( 4  )
                chmID = 'WetAeroAreaMDUST4'
             CASE( 5  )
                chmID = 'WetAeroAreaMDUST5'
             CASE( 6  )
                chmID = 'WetAeroAreaMDUST6'
             CASE( 7  )
                chmID = 'WetAeroAreaMDUST7'
             CASE( 8  )
                chmID = 'WetAeroAreaSULF'
             CASE( 9  )
                chmID = 'WetAeroAreaBC'
             CASE( 10 )
                chmID = 'WetAeroAreaOC'
             CASE( 11 )
                chmID = 'WetAeroAreaSSA'
             CASE( 12 )
                chmID = 'WetAeroAreaSSC'
             CASE( 13 )
                chmID = 'WetAeroAreaBGSULF'
             CASE( 14 )
                chmID = 'WetAeroAreaICEI'
             CASE DEFAULT
                ErrMsg = 'State_Chm%nAeroType exceeds the number of defined' &
                         // ' wet aerosol area categories'
                CALL GC_Error( ErrMsg, RC, ThisLoc )
                RETURN
          END SELECT

          CALL Register_ChmField( Input_Opt, chmID, State_Chm%WetAeroArea,   &
                                  State_Chm, RC,    Ncat=N                  )
          CALL GC_CheckVar( 'State_Chm%WetAeroArea', 1, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDDO

       !--------------------------------------------------------------------
       ! WetAeroRadi
       !--------------------------------------------------------------------
       ALLOCATE( State_Chm%WetAeroRadi(IM,JM,LM,State_Chm%nAeroType), STAT=RC )
       CALL GC_CheckVar( 'State_Chm%WetAeroRadi', 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Chm%WetAeroRadi = 0.0_fp

       ! Loop over all entries to register each category individually
       DO N = 1, State_Chm%nAeroType

          ! Define identifying string
          SELECT CASE( N )
             CASE( 1  )
                chmID = 'WetAeroRadiMDUST1'
             CASE( 2  )
                chmID = 'WetAeroRadiMDUST2'
             CASE( 3  )
                chmID = 'WetAeroRadiMDUST3'
             CASE( 4  )
                chmID = 'WetAeroRadiMDUST4'
             CASE( 5  )
                chmID = 'WetAeroRadiMDUST5'
             CASE( 6  )
                chmID = 'WetAeroRadiMDUST6'
             CASE( 7  )
                chmID = 'WetAeroRadiMDUST7'
             CASE( 8  )
                chmID = 'WetAeroRadiSULF'
             CASE( 9  )
                chmID = 'WetAeroRadiBC'
             CASE( 10 )
                chmID = 'WetAeroRadiOC'
             CASE( 11 )
                chmID = 'WetAeroRadiSSA'
             CASE( 12 )
                chmID = 'WetAeroRadiSSC'
             CASE( 13 )
                chmID = 'WetAeroRadiBGSULF'
             CASE( 14 )
                chmID = 'WetAeroRadiICEI'
             CASE DEFAULT
                ErrMsg = 'State_Chm%nAeroType exceeds the number of defined' &
                         // ' wet aerosol radius categories'
                CALL GC_Error( ErrMsg, RC, ThisLoc )
                RETURN
          END SELECT

          CALL Register_ChmField( Input_Opt, chmID, State_Chm%WetAeroRadi,   &
                                  State_Chm, RC,    Ncat=N )
          CALL GC_CheckVar( 'State_Chm%WetAeroRadi', 1, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDDO

       !--------------------------------------------------------------------
       ! AeroH2O
       !--------------------------------------------------------------------
       ALLOCATE( State_Chm%AeroH2O( IM, JM, LM, State_Chm%nAeroType ), STAT=RC )
       CALL GC_CheckVar( 'State_Chm%AeroH2O', 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Chm%AeroH2O = 0.0_fp

       ! Loop over all entries to register each category individually
       DO N = 1, State_Chm%nAeroType

          ! Define identifying string
          SELECT CASE( N )
             CASE( 1  )
                chmID = 'AeroH2OMDUST1'
             CASE( 2  )
                chmID = 'AeroH2OMDUST2'
             CASE( 3  )
                chmID = 'AeroH2OMDUST3'
             CASE( 4  )
                chmID = 'AeroH2OMDUST4'
             CASE( 5  )
                chmID = 'AeroH2OMDUST5'
             CASE( 6  )
                chmID = 'AeroH2OMDUST6'
             CASE( 7  )
                chmID = 'AeroH2OMDUST7'
             CASE( 8  )
                chmID = 'AeroH2OSULF'
             CASE( 9  )
                chmID = 'AeroH2OBC'
             CASE( 10 )
                chmID = 'AeroH2OOC'
             CASE( 11 )
                chmID = 'AeroH2OSSA'
             CASE( 12 )
                chmID = 'AeroH2OSSC'
             CASE( 13 )
                chmID = 'AeroH2OBGSULF'
             CASE( 14 )
                chmID = 'AeroH2OICEI'
             CASE DEFAULT
                ErrMsg = 'State_Chm%nAeroType exceeds the number of defined' &
                         // ' aerosol H2O categories'
                CALL GC_Error( ErrMsg, RC, ThisLoc )
                RETURN
          END SELECT

          CALL Register_ChmField( Input_Opt, chmID, State_Chm%AeroH2O,   &
                                  State_Chm, RC,    Ncat=N )
          CALL GC_CheckVar( 'State_Chm%AeroH2O', 1, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDDO

       !--------------------------------------------------------------------
       ! GammaN2O5
       !--------------------------------------------------------------------
       ALLOCATE( State_Chm%GammaN2O5( IM, JM, LM, 4 ), STAT=RC )
       CALL GC_CheckVar( 'State_Chm%GammaN2O5', 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Chm%GammaN2O5 = 0.0_fp

       ! Loop over all entries to register each category individually
       DO N = 1, 4

          ! Define identifying string
          SELECT CASE( N )
             CASE( 1  )
                chmID = 'GammaN2O5H2O'
             CASE( 2  )
                chmID = 'GammaN2O5HCl'
             CASE( 3  )
                chmID = 'GammaN2O5SS'
             CASE( 4  )
                chmID = 'YieldClNO2'
             CASE DEFAULT
                ErrMsg = 'State_Chm%GammaN2O5 exceeds the number of defined' &
                         // ' N2O5 uptake categories'
                CALL GC_Error( ErrMsg, RC, ThisLoc )
                RETURN
          END SELECT

          CALL Register_ChmField( Input_Opt, chmID, State_Chm%GammaN2O5,     &
                                  State_Chm, RC,    Ncat=N )
          CALL GC_CheckVar( 'State_Chm%GammaN2O5', 1, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDDO

       !--------------------------------------------------------------------
       ! OM:OC Ratios
       !--------------------------------------------------------------------
       chmId = 'OMOCpoa'
       ALLOCATE( State_Chm%OMOC_POA( IM, JM ), STAT=RC )
       CALL GC_CheckVar( 'State_Chm%OMOC_POA', 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Chm%OMOC_POA = 0.0_fp
       CALL Register_ChmField( Input_Opt, chmID, State_Chm%OMOC_POA,         &
                               State_Chm, RC                                )
       CALL GC_CheckVar( 'State_Chm%OMOC_POA', 1, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       !--------------------------------------------------------------------
       chmId = 'OMOCopoa'
       ALLOCATE( State_Chm%OMOC_OPOA( IM, JM ), STAT=RC )
       CALL GC_CheckVar( 'State_Chm%OMOC_OPOA', 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Chm%OMOC_OPOA = 0.0_fp
       CALL Register_ChmField( Input_Opt, chmID, State_Chm%OMOC_OPOA,        &
                               State_Chm, RC                                )
       CALL GC_CheckVar( 'State_Chm%OMOC_OPOA', 1, RC )
       IF ( RC /= GC_SUCCESS ) RETURN

       !--------------------------------------------------------------------
       ! phSav
       !--------------------------------------------------------------------
       chmId = 'phSav'
       ALLOCATE( State_Chm%phSav( IM, JM, LM ), STAT=RC )
       CALL GC_CheckVar( 'State_Chm%phSav', 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Chm%phSav = 0.0_fp
       CALL Register_ChmField( Input_Opt, chmID, State_Chm%phSav,            &
                               State_Chm, RC                                )
       CALL GC_CheckVar( 'State_Chm%phSav', 1, RC )
       IF ( RC /= GC_SUCCESS ) RETURN

       !--------------------------------------------------------------------
       ! HplusSav
       !--------------------------------------------------------------------
       chmID  = 'HplusSav'
       ALLOCATE( State_Chm%HplusSav( IM, JM, LM ), STAT=RC )
       CALL GC_CheckVar( 'State_Chm%HplusSav', 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Chm%HplusSav = 0.0_fp
       CALL Register_ChmField( Input_Opt, chmID, State_Chm%HplusSav,         &
                               State_Chm, RC                                )
       CALL GC_CheckVar( 'State_Chm%HplusSav', 1, RC )
       IF ( RC /= GC_SUCCESS ) RETURN

       !--------------------------------------------------------------------
       ! WaterSav
       !--------------------------------------------------------------------
       chmID  = 'WaterSav'
       ALLOCATE( State_Chm%WaterSav( IM, JM, LM ), STAT=RC )
       CALL GC_CheckVar( 'State_Chm%WaterSav', 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Chm%WaterSav = 0.0_fp
       CALL Register_ChmField( Input_Opt, chmID, State_Chm%WaterSav,         &
                               State_Chm, RC                                )
       CALL GC_CheckVar( 'State_Chm%WaterSav', 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN

       !--------------------------------------------------------------------
       ! SulRatSav
       !--------------------------------------------------------------------
       chmID  = 'SulRatSav'
       ALLOCATE( State_Chm%SulRatSav( IM, JM, LM ), STAT=RC )
       CALL GC_CheckVar( 'State_Chm%SulRatSav', 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Chm%SulRatSav = 0.0_fp
       CALL Register_ChmField( Input_Opt, chmID, State_Chm%SulRatSav,        &
                               State_Chm, RC                                )
       CALL GC_CheckVar( 'State_Chm%SulRatSav', 1, RC )
       IF ( RC /= GC_SUCCESS ) RETURN

       !--------------------------------------------------------------------
       ! NaRatSav
       !--------------------------------------------------------------------
       chmID  = 'NaRatSav'
       ALLOCATE( State_Chm%NaRatSav( IM, JM, LM ), STAT=RC )
       CALL GC_CheckVar( 'State_Chm%NaRatSav', 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Chm%NaRatSav = 0.0_fp
       CALL Register_ChmField( Input_Opt, chmID, State_Chm%NaRatSav,         &
                               State_Chm, RC                                )
       CALL GC_CheckVar( 'State_Chm%NaRatSav', 1, RC )
       IF ( RC /= GC_SUCCESS ) RETURN

       !--------------------------------------------------------------------
       ! AcidPurSav
       !--------------------------------------------------------------------
       chmID  = 'AcidPurSav'
       ALLOCATE( State_Chm%AcidPurSav( IM, JM, LM ), STAT=RC )
       CALL GC_CheckVar( 'State_Chm%AcidPurSav', 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Chm%AcidPurSav = 0.0_fp
       CALL Register_ChmField( Input_Opt, chmID, State_Chm%AcidPurSav,       &
                               State_Chm, RC                                )
       CALL GC_CheckVar( 'State_Chm%AcidPurSav', 1, RC )
       IF ( RC /= GC_SUCCESS ) RETURN

       !--------------------------------------------------------------------
       ! BisulSav
       !--------------------------------------------------------------------
       chmID  = 'BisulSav'
       ALLOCATE( State_Chm%BisulSav( IM, JM, LM ), STAT=RC )
       CALL GC_CheckVar( 'State_Chm%BiSulSav', 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Chm%BisulSav = 0.0_fp
       CALL Register_ChmField( Input_Opt, chmID, State_Chm%BisulSav,         &
                               State_Chm, RC                                )
       CALL GC_CheckVar( 'State_Chm%BiSulSav', 1, RC )
       IF ( RC /= GC_SUCCESS ) RETURN

       !--------------------------------------------------------------------
       ! phCloud
       !--------------------------------------------------------------------
       chmId = 'pHCloud'
       ALLOCATE( State_Chm%pHCloud( IM, JM, LM ), STAT=RC )
       CALL GC_CheckVar( 'State_Chm%pHCloud', 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Chm%pHCloud = 0.0_fp
       CALL Register_ChmField( Input_Opt, chmID, State_Chm%pHCloud,          &
                               State_Chm, RC                                )
       CALL GC_CheckVar( 'State_Chm%pHCloud', 1, RC )
       IF ( RC /= GC_SUCCESS ) RETURN

       !--------------------------------------------------------------------
       ! isCloud
       ! jmm 3/1/19
       !--------------------------------------------------------------------
       chmId = 'isCloud'
       ALLOCATE( State_Chm%isCloud( IM, JM, LM ), STAT=RC )
       CALL GC_CheckVar( 'State_Chm%isCloud', 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Chm%isCloud = 0.0_fp
       CALL Register_ChmField( Input_Opt, chmID, State_Chm%isCloud,          &
                               State_Chm, RC                                )
       CALL GC_CheckVar( 'State_Chm%isCloud', 1, RC )
       IF ( RC /= GC_SUCCESS ) RETURN

       !--------------------------------------------------------------------
       ! SSAlk
       !--------------------------------------------------------------------
       ALLOCATE( State_Chm%SSAlk( IM, JM, LM, 2 ), STAT=RC )
       CALL GC_CheckVar( 'State_Chm%SSAlk', 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Chm%SSAlk = 0.0_fp

       ! Register accumulation mode as category 1
       chmId = 'SSAlkAccum'
       CALL Register_ChmField( Input_Opt, chmID, State_Chm%SSAlk,            &
                               State_Chm, RC,    nCat=1                     )
       CALL GC_CheckVar( 'State-Chm%SSAlk', 1, RC )
       IF ( RC /= GC_SUCCESS ) RETURN

       ! Register coarse mode as category 1
       chmId = 'SSAlkCoarse'
       CALL Register_ChmField( Input_Opt, chmID, State_Chm%SSAlk,            &
                               State_Chm, RC,    nCat=2                     )
       CALL GC_CheckVar( 'State_Chm%SSAlk', 1, RC )
       IF ( RC /= GC_SUCCESS ) RETURN

       !------------------------------------------------------------------
       ! HSO3_AQ
       !------------------------------------------------------------------
       chmID = 'HSO3AQ'
       ALLOCATE( State_Chm%HSO3_AQ( IM, JM, LM ) , STAT=RC )
       CALL GC_CheckVar( 'State_Chm%HSO3_AQ', 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Chm%HSO3_AQ = 0.0_fp
       CALL Register_ChmField( Input_Opt, chmID, State_Chm%HSO3_AQ,          &
                               State_Chm, RC                                )
       CALL GC_CheckVar( 'State_Chm%HSO3_AQ', 1, RC )
       IF ( RC /= GC_SUCCESS ) RETURN

       !------------------------------------------------------------------
       ! SO3_AQ
       !------------------------------------------------------------------
       chmID = 'SO3AQ'
       ALLOCATE( State_Chm%SO3_AQ( IM, JM, LM ) , STAT=RC )
       CALL GC_CheckVar( 'State_Chm%SO3_AQ', 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Chm%SO3_AQ = 0.0_fp
       CALL Register_ChmField( Input_Opt, chmID, State_Chm%SO3_AQ,           &
                               State_Chm, RC                                )
       CALL GC_CheckVar( 'State_Chm%SO3_AQ', 1, RC )
       IF ( RC /= GC_SUCCESS ) RETURN

       !------------------------------------------------------------------
       ! fupdateHOBr
       !------------------------------------------------------------------
       chmID = 'fupdateHOBr'
       ALLOCATE( State_Chm%fupdateHOBr( IM, JM, LM ) , STAT=RC )
       CALL GC_CheckVar( 'State_Chm%fupdateHOBr', 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Chm%fupdateHOBr = 0.0_fp
       CALL Register_ChmField( Input_Opt, chmID, State_Chm%fupdateHOBr,     &
                               State_Chm, RC                               )
       CALL GC_CheckVar( 'State_Chm%fupdateHOBr', 1, RC )
       IF ( RC /= GC_SUCCESS ) RETURN

       !------------------------------------------------------------------
       ! DryDepNitrogen
       !------------------------------------------------------------------
       chmID = 'DryDepNitrogen'
       ALLOCATE( State_Chm%DryDepNitrogen( IM, JM ) , STAT=RC )
       CALL GC_CheckVar( 'State_Chm%DryDepNitrogen', 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Chm%DryDepNitrogen = 0.0_fp
       CALL Register_ChmField( Input_Opt, chmID, State_Chm%DryDepNitrogen,   &
                               State_Chm, RC                                )
       CALL GC_CheckVar( 'State_Chm%DryDepNitrogen', 1, RC )
       IF ( RC /= GC_SUCCESS ) RETURN

       !------------------------------------------------------------------
       ! WetDepNitrogen
       !------------------------------------------------------------------
       chmID = 'WetDepNitrogen'
       ALLOCATE( State_Chm%WetDepNitrogen( IM, JM ) , STAT=RC )
       CALL GC_CheckVar( 'State_Chm%WetDepNitrogen', 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Chm%WetDepNitrogen = 0.0_fp
       CALL Register_ChmField( Input_Opt, chmID, State_Chm%WetDepNitrogen,   &
                               State_Chm, RC                                )
       CALL GC_CheckVar( 'State_Chm%WetDepNitrogen', 1, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    !=======================================================================
    ! Allocate and initialize quantities for wet deposition routines
    !=======================================================================

    !------------------------------------------------------------------
    ! H2O2AfterChem
    !------------------------------------------------------------------
    chmID = 'H2O2AfterChem'
    ALLOCATE( State_Chm%H2O2AfterChem( IM, JM, LM ) , STAT=RC )
    CALL GC_CheckVar( 'State_Chm%H2O2AfterChem', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Chm%H2O2AfterChem = 0.0_fp
    CALL Register_ChmField( Input_Opt, chmID, State_Chm%H2O2AfterChem,    &
                            State_Chm, RC                                )
    CALL GC_CheckVar( 'State_Chm%H2O2AfterChem', 1, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    !------------------------------------------------------------------
    ! SO2AfterChem
    !------------------------------------------------------------------
    chmID = 'SO2AfterChem'
    ALLOCATE( State_Chm%SO2AfterChem( IM, JM, LM ) , STAT=RC )
    CALL GC_CheckVar( 'State_Chm%SO2AfterChem', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Chm%SO2AfterChem = 0.0_fp
    CALL Register_ChmField( Input_Opt, chmID, State_Chm%SO2AfterChem,     &
                            State_Chm, RC                                )
    CALL GC_CheckVar( 'State_Chm%SO2AfterChem', 1, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    !=======================================================================
    ! Allocate and initialize fields for KPP solver
    !=======================================================================
    IF ( Input_Opt%ITS_A_FULLCHEM_SIM ) THEN

       !--------------------------------------------------------------------
       ! KPPHvalue
       !--------------------------------------------------------------------
       chmID = 'KPPHvalue'
       ALLOCATE( State_Chm%KPPHvalue( IM, JM, LM ), STAT=RC )
       CALL GC_CheckVar( 'State_Chm%KPPHvalue', 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Chm%KPPHvalue = 0.0_fp
       CALL Register_ChmField( Input_Opt, chmID, State_Chm%KPPHvalue,        &
                               State_Chm, RC                                )
       CALL GC_CheckVar( 'State_Chm%KPPHvalue', 1, RC )
       IF ( RC /= GC_SUCCESS ) RETURN

    ENDIF

    !=======================================================================
    ! Allocate and initialize fields for UCX mechamism
    !=======================================================================
    IF ( Input_Opt%ITS_A_FULLCHEM_SIM .and. Input_Opt%LUCX ) THEN

       !--------------------------------------------------------------------
       ! STATE_PSC
       !--------------------------------------------------------------------
       chmID = 'StatePSC'
       ALLOCATE( State_Chm%STATE_PSC( IM, JM, LM ), STAT=RC )
       CALL GC_CheckVar( 'State_Chm%STATE_PSC', 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Chm%STATE_PSC = 0.0_f4
       CALL Register_ChmField( Input_Opt, chmID, State_Chm%STATE_PSC,        &
                            State_Chm, RC )
       CALL GC_CheckVar( 'State_Chm%STATE_PSC', 1, RC )
       IF ( RC /= GC_SUCCESS ) RETURN

       !--------------------------------------------------------------------
       ! KHETI_SLA
       !-------------------------------------------------------------------
       nKHLSA = 11
       ALLOCATE( State_Chm%KHETI_SLA ( IM, JM, LM, nKHLSA ), STAT=RC )
       CALL GC_CheckVar( 'State_Chm%KHETISLA', 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Chm%KHETI_SLA = 0.0_fp

       ! Loop over all entries to register each category individually
       DO N = 1, nKHLSA

          ! Define identifying string
          SELECT CASE( N )
             CASE( 1  )
                chmID = 'KhetiSLAN2O5H2O'
             CASE( 2  )
                chmID = 'KhetiSLAN2O5HCl'
             CASE( 3  )
                chmID = 'KhetiSLAClNO3H2O'
             CASE( 4  )
                chmID = 'KhetiSLAClNO3HCl'
             CASE( 5  )
                chmID = 'KhetiSLAClNO3HBr'
             CASE( 6  )
                chmID = 'KhetiSLABrNO3H2O'
             CASE( 7  )
                chmID = 'KhetiSLABrNO3HCl'
             CASE( 8  )
                chmID = 'KhetiSLAHOClHCl'
             CASE( 9  )
                chmID = 'KhetiSLAHOClHBr'
             CASE( 10 )
                chmID = 'KhetiSLAHOBrHCl'
             CASE( 11 )
                chmID = 'KhetiSLAHOBrHBr'
             CASE DEFAULT
                ErrMsg = 'nKHLSA exceeds the number of defined' &
                       // ' KHETI_SLA categories'
                CALL GC_Error( ErrMsg, RC, ThisLoc )
                RETURN
          END SELECT

          CALL Register_ChmField( Input_Opt, chmID, State_Chm%KHETI_SLA, &
                                  State_Chm, RC,    Ncat=N )
          CALL GC_CheckVar( 'State_Chm%KHETISLA', 1, RC )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDDO
    ENDIF

    !=======================================================================
    ! Special handling for the Hg and tagHg simulations: get the # of Hg
    ! categories for total & tagged tracers from the species database
    !=======================================================================
    IF ( Input_Opt%ITS_A_MERCURY_SIM ) THEN

       ! Hg0, Hg2, HgP should all have the same number of categories as
       ! returned from the species database.  If not, there's an error.
       IF ( N_Hg0_CATS == N_Hg2_CATS .and. N_Hg0_CATS == N_HgP_CATS ) THEN
          State_Chm%N_Hg_CATS = N_Hg0_CATS
       ELSE
          ErrMsg = 'Inconsistent number of Hg categories!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

       ! Index array: Hg0 species # <--> Hg0 category #
       ALLOCATE( State_Chm%Hg0_Id_List( State_Chm%N_Hg_CATS ), STAT=RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Chm%Hg0_Id_List = 0

       ! Index array: Hg2 species # <--> Hg0 category #
       ALLOCATE( State_Chm%Hg2_Id_List( State_Chm%N_Hg_CATS ), STAT=RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Chm%Hg2_Id_List = 0

       ! Index array: HgP species # <--> Hg0 category #
       ALLOCATE( State_Chm%HgP_Id_List( State_Chm%N_Hg_CATS ), STAT=RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Chm%HgP_Id_List = 0

       ! Hg category names
       ALLOCATE( State_Chm%Hg_Cat_Name( State_Chm%N_Hg_CATS ), STAT=RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Chm%Hg_Cat_Name = ''

       ! Loop over all species
       DO N = 1, State_Chm%nSpecies

          ! Point to Species Database entry for Hg species N
          ThisSpc => State_Chm%SpcData(N)%Info

          ! Populate the Hg0 index array
          IF ( ThisSpc%Is_Hg0 ) THEN
             State_Chm%Hg0_Id_List(ThisSpc%Hg_Cat) = ThisSpc%ModelId
          ENDIF

          ! Populate the Hg2 index array
          IF ( ThisSpc%Is_Hg2 ) THEN
             State_Chm%Hg2_Id_List(ThisSpc%Hg_Cat) = ThisSpc%ModelId
          ENDIF

          ! Populate the HgP index array
          IF ( ThisSpc%Is_HgP ) THEN
             State_Chm%HgP_Id_List(ThisSpc%Hg_Cat) = ThisSpc%ModelId
          ENDIF

          ! Free pointer
          ThisSpc => NULL()
       ENDDO

       ! Loop over Hg categories (except the first
       DO C = 2, State_Chm%N_Hg_CATS

          ! Hg0 tracer number corresponding to this category
          N                        =  State_Chm%Hg0_Id_List(C)

          ! The category name (e.g. "_can") follows the "Hg0"
          ThisSpc                  => State_Chm%SpcData(N)%Info
          State_Chm%Hg_Cat_Name(C) =  ThisSpc%Name(4:7)
          ThisSpc                  => NULL()
       ENDDO

       !--------------------------------------------------------------------
       ! Hg(0) ocean mass
       !--------------------------------------------------------------------
       chmID = 'OceanHg0'
       ALLOCATE( State_Chm%OceanHg0( IM, JM, State_Chm%N_Hg_CATS ), STAT=RC )
       CALL GC_CheckVar( 'State_Chm%OceanHg0', 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Chm%OceanHg0 = 0.0_fp
       CALL Register_ChmField( Input_Opt, chmID, State_Chm%OceanHg0,      &
                               State_Chm, RC                             )
       CALL GC_CheckVar( 'State_Chm%OceanHg0', 1, RC )
       IF ( RC /= GC_SUCCESS ) RETURN

       !--------------------------------------------------------------------
       ! Hg(II) ocean mass
       !--------------------------------------------------------------------
       chmID = 'OceanHg2'
       ALLOCATE( State_Chm%OceanHg2( IM, JM, State_Chm%N_Hg_CATS ), STAT=RC )
       CALL GC_CheckVar( 'State_Chm%OceanHg2', 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Chm%OceanHg2 = 0.0_fp
       CALL Register_ChmField( Input_Opt, chmID, State_Chm%OceanHg2,      &
                               State_Chm, RC                             )
       CALL GC_CheckVar( 'State_Chm%OceanHg2', 1, RC )
       IF ( RC /= GC_SUCCESS ) RETURN

       !--------------------------------------------------------------------
       ! HgP ocean mass
       !--------------------------------------------------------------------
       chmID = 'OceanHgP'
       ALLOCATE( State_Chm%OceanHgP (IM, JM, State_Chm%N_Hg_CATS ), STAT=RC )
       CALL GC_CheckVar( 'State_Chm%OceanHgP', 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Chm%OceanHgP = 0.0_fp
       CALL Register_ChmField( Input_Opt, chmID, State_Chm%OceanHgP,      &
                               State_Chm, RC                             )
       CALL GC_CheckVar( 'State_Chm%OceanHgP', 1, RC )
       IF ( RC /= GC_SUCCESS ) RETURN

       !--------------------------------------------------------------------
       ! Reducible Hg snowpack on ocean
       !--------------------------------------------------------------------
       chmID = 'SnowHgOcean'
       ALLOCATE( State_Chm%SnowHgOcean( IM, JM, State_Chm%N_Hg_CATS ), STAT=RC )
       CALL GC_CheckVar( 'State_Chm%SnowHgOcean', 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Chm%SnowHgOcean = 0.0_fp
       CALL Register_ChmField( Input_Opt, chmID, State_Chm%SnowHgOcean,   &
                               State_Chm, RC                             )
       CALL GC_CheckVar( 'State_Chm%SnowHgOcean', 1, RC )
       IF ( RC /= GC_SUCCESS ) RETURN

       !--------------------------------------------------------------------
       ! Reducible Hg snowpack on land
       !--------------------------------------------------------------------
       chmID = 'SnowHgLand'
       ALLOCATE( State_Chm%SnowHgLand( IM, JM, State_Chm%N_Hg_CATS ), STAT=RC )
       CALL GC_CheckVar( 'State_Chm%SnowHgLand', 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Chm%SnowHgLand = 0.0_fp
       CALL Register_ChmField( Input_Opt, chmID, State_Chm%SnowHgLand,    &
                               State_Chm, RC                             )
       CALL GC_CheckVar( 'State_Chm%SnowHgLand', 1, RC )
       IF ( RC /= GC_SUCCESS ) RETURN

       !--------------------------------------------------------------------
       ! Non-reducible Hg snowpack on ocean
       !--------------------------------------------------------------------
       chmID = 'SnowHgOceanStored'
       ALLOCATE( State_Chm%SnowHgOceanStored(IM, JM, State_Chm%N_Hg_CATS ), &
                 STAT=RC )
       CALL GC_CheckVar( 'State_Chm%SnowHgOceanStored', 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Chm%SnowHgOceanStored = 0.0_fp
       CALL Register_ChmField( Input_Opt, chmID, State_Chm%SnowHgOceanStored, &
                               State_Chm, RC                             )
       CALL GC_CheckVar( 'State_Chm%SnowHgOceanStored', 1, RC )
       IF ( RC /= GC_SUCCESS ) RETURN

       !--------------------------------------------------------------------
       ! Non-reducible Hg snowpack on land
       !--------------------------------------------------------------------
       chmID = 'SnowHgLandStored'
       ALLOCATE( State_Chm%SnowHgLandStored(IM, JM, State_Chm%N_Hg_CATS ),   &
                 STAT=RC )
       CALL GC_CheckVar( 'State_Chm%SnowHgLandStored', 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Chm%SnowHgLandStored = 0.0_fp
       CALL Register_ChmField( Input_Opt, chmID, State_Chm%SnowHgLandStored, &
                               State_Chm, RC                                )
       CALL GC_CheckVar( 'State_Chm%SnowHgLandStored', 1, RC )
       IF ( RC /= GC_SUCCESS ) RETURN

    ENDIF


    !=======================================================================
    ! Allocate fields for various GeosCore modules
    !=======================================================================
    !------------------------------------------------------------------
    ! DryDepSav
    !------------------------------------------------------------------
    IF ( State_Chm%nDryDep > 0 ) THEN
        chmID = 'DryDepSav'
        ALLOCATE( State_Chm%DryDepSav( IM, JM, State_Chm%nDryDep ) , STAT=RC )
        CALL GC_CheckVar( 'State_Chm%DryDepSav', 0, RC )
        IF ( RC /= GC_SUCCESS ) RETURN
        State_Chm%DryDepSav = 0.0_fp
        CALL Register_ChmField( Input_Opt, chmID, State_Chm%DryDepSav,       &
                                State_Chm, RC                               )
        CALL GC_CheckVar( 'State_Chm%DryDepSav', 1, RC )
        IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    !------------------------------------------------------------------
    ! TLSTT (Linoz)
    !------------------------------------------------------------------
    IF ( Input_Opt%LLINOZ .AND. Input_Opt%LINOZ_NFIELDS > 0 ) THEN
        chmID = 'TLSTT'
        ALLOCATE( State_Chm%TLSTT( IM, JM, LM, Input_Opt%LINOZ_NFIELDS ),    &
                  STAT=RC )
        CALL GC_CheckVar( 'State_Chm%TLSTT', 0, RC )
        IF ( RC /= GC_SUCCESS ) RETURN
        State_Chm%TLSTT = 0.0_fp

        ! Do not register this field as it is internal
        ! to the linoz_mod module state. (hplin, 1/24/19)
        ! Note: We might want to implement support for implementing a 4th
        ! dimension later.
    ENDIF

    !------------------------------------------------------------------
    ! Gan Luo et al wetdep fields
    !------------------------------------------------------------------
    IF ( Input_Opt%LWETD .or. Input_Opt%LCONV ) THEN

        ! PSO4s
        chmID = 'PSO4s'
        ALLOCATE( State_Chm%PSO4s( IM, JM, LM ), STAT=RC )
        CALL GC_CheckVar( 'State_Chm%PSO4s', 0, RC )
        IF ( RC /= GC_SUCCESS ) RETURN
        State_Chm%PSO4s = 0.0_fp
        CALL Register_ChmField( Input_Opt, chmID, State_Chm%PSO4s,           &
                                State_Chm, RC                               )
        CALL GC_CheckVar( 'State_Chm%PSO4s', 1, RC )
        IF ( RC /= GC_SUCCESS ) RETURN

        ! QQ3D
        chmID = 'QQ3D'
        ALLOCATE( State_Chm%QQ3D( IM, JM, LM ), STAT=RC )
        CALL GC_CheckVar( 'State_Chm%QQ3D', 0, RC )
        IF ( RC /= GC_SUCCESS ) RETURN
        State_Chm%QQ3D = 0.0_fp
        CALL Register_ChmField( Input_Opt, chmID, State_Chm%QQ3D,            &
                                State_Chm, RC                               )
        CALL GC_CheckVar( 'State_Chm%QQ3D', 1, RC )
        IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    !------------------------------------------------------------------
    ! Surface flux for non-local PBL mixing
    !------------------------------------------------------------------
    IF ( Input_Opt%LTURB .and. Input_Opt%LNLPBL ) THEN
       chmID = 'SurfaceFlux'
       ALLOCATE( State_Chm%SurfaceFlux( IM, JM, State_Chm%nAdvect ), STAT=RC )
       CALL GC_CheckVar( 'State_Chm%SurfaceFlux', 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Chm%SurfaceFlux = 0.0_fp
       CALL Register_ChmField( Input_Opt, chmID, State_Chm%SurfaceFlux,      &
                               State_Chm, RC                                )
       CALL GC_CheckVar( 'State_Chm%SurfaceFlux', 1, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    !========================================================================
    ! Once we are done registering all fields, we need to define the
    ! registry lookup table.  This algorithm will avoid hash collisions.
    !========================================================================
    CALL Registry_Set_LookupTable( Registry  = State_Chm%Registry,           &
                                   RegDict   = State_Chm%RegDict,            &
                                   RC        = RC                           )

    ! Trap potential errors
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Error encountered in routine "Registry_Set_LookupTable"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    !========================================================================
    ! Print out the list of registered fields
    !========================================================================
    IF ( Input_Opt%amIRoot ) THEN
       WRITE( 6, 10 )
10     FORMAT( /, 'Registered variables contained within the State_Chm object:')
       WRITE( 6, '(a)' ) REPEAT( '=', 79 )
    ENDIF

    ! Print registered fields
    CALL Registry_Print( Input_Opt   = Input_Opt,             &
                         Registry    = State_Chm%Registry,    &
                         ShortFormat = .TRUE.,                &
                         RC          = RC                    )

    ! Trap potential errors
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Error encountered in routine "Registry Print"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    !=======================================================================
    ! Cleanup and quit
    !=======================================================================

    ! Free pointer for safety's sake
    ThisSpc => NULL()

    ! Echo output
    IF ( Input_Opt%amIRoot ) THEN
       print*, REPEAT( '#', 79 )
    ENDIF

    ! Format statement
100 FORMAT( I3, 2x, A31 )
110 FORMAT( 5x, '===> ', f4.1, 1x, A6  )
120 FORMAT( 5x, '---> ', f4.1, 1x, A4  )

  END SUBROUTINE Init_State_Chm
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Cleanup_State_Chm
!
! !DESCRIPTION: Routine CLEANUP\_STATE\_CHM deallocates all fields
!  of the chemistry state object.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Cleanup_State_Chm( State_Chm, RC )
!
! !USES:
!
    USE Species_Database_Mod, ONLY : Cleanup_Species_Database
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ChmState), INTENT(INOUT) :: State_Chm    ! Chemistry State object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC           ! Return code
!
! !REMARKS:
!
! !REVISION HISTORY:
!  15 Oct 2012 - R. Yantosca - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARAIBLES
!
    ! Strings
    CHARACTER(LEN=255) :: ErrMsg, ThisLoc

    !=======================================================================
    ! Initialize
    !=======================================================================
    RC      = GC_SUCCESS
    ErrMsg  = ''
    ThisLoc = ' -> at Cleanup_State_Chm (in Headers/state_chm_mod.F90)'

    !=======================================================================
    ! Deallocate and nullify pointer fields of State_Chm
    !=======================================================================
    IF ( ASSOCIATED( State_Chm%Map_Advect ) ) THEN
       DEALLOCATE( State_Chm%Map_Advect, STAT=RC )
       CALL GC_CheckVar( 'State_Chm%Map_Advect', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Chm%Map_Advect => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Chm%Map_Aero ) ) THEN
       DEALLOCATE( State_Chm%Map_Aero, STAT=RC )
       CALL GC_CheckVar( 'State_Chm%Map_Aero', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Chm%Map_Aero => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Chm%Map_DryDep ) ) THEN
       DEALLOCATE( State_Chm%Map_DryDep, STAT=RC )
       CALL GC_CheckVar( 'State_Chm%Map_Drydep', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Chm%Map_DryDep => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Chm%Map_GasSpc ) ) THEN
       DEALLOCATE( State_Chm%Map_GasSpc, STAT=RC )
       CALL GC_CheckVar( 'State_Chm%Map_GasSpc', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Chm%Map_GasSpc => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Chm%Map_HygGrth ) ) THEN
       DEALLOCATE( State_Chm%Map_HygGrth, STAT=RC )
       CALL GC_CheckVar( 'State_Chm%Map_HygGrth', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Chm%Map_HygGrth => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Chm%Map_KppVar ) ) THEN
       DEALLOCATE( State_Chm%Map_KppVar, STAT=RC )
       CALL GC_CheckVar( 'State_Chm%Map_KppVar', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Chm%Map_KppVar => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Chm%Map_KppFix ) ) THEN
       DEALLOCATE( State_Chm%Map_KppFix, STAT=RC )
       CALL GC_CheckVar( 'State_Chm%Map_KppFix', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Chm%Map_KppFix => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Chm%Map_KppSpc ) ) THEN
       DEALLOCATE( State_Chm%Map_KppSpc, STAT=RC )
       CALL GC_CheckVar( 'State_Chm%Map_KppSpc', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Chm%Map_KppSpc => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Chm%Name_Loss ) ) THEN
       DEALLOCATE( State_Chm%Name_Loss, STAT=RC )
       CALL GC_CheckVar( 'State_Chm%Name_Loss', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Chm%Name_Loss => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Chm%Map_Loss ) ) THEN
       DEALLOCATE( State_Chm%Map_Loss, STAT=RC )
       CALL GC_CheckVar( 'State_Chm%Map_Loss', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Chm%Map_Loss => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Chm%Map_Photol ) ) THEN
       DEALLOCATE( State_Chm%Map_Photol, STAT=RC )
       CALL GC_CheckVar( 'State_Chm%Map_Photol', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Chm%Map_Photol => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Chm%Name_Prod ) ) THEN
       DEALLOCATE( State_Chm%Name_Prod, STAT=RC )
       CALL GC_CheckVar( 'State_Chm%Name_Prod', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Chm%Name_Prod => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Chm%Map_Prod ) ) THEN
       DEALLOCATE( State_Chm%Map_Prod, STAT=RC )
       CALL GC_CheckVar( 'State_Chm%Map_Prod', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Chm%Map_Prod => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Chm%Map_WetDep ) ) THEN
       DEALLOCATE( State_Chm%Map_WetDep, STAT=RC )
       CALL GC_CheckVar( 'State_Chm%Map_WetDep', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Chm%Map_WetDep => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Chm%Map_WL ) ) THEN
       DEALLOCATE( State_Chm%Map_WL, STAT=RC )
       CALL GC_CheckVar( 'State_Chm%Map_WL', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Chm%Map_WL => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Chm%Species ) ) THEN
       DEALLOCATE( State_Chm%Species, STAT=RC )
       CALL GC_CheckVar( 'State_Chm%Species', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Chm%Species => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Chm%BoundaryCond ) ) THEN
       DEALLOCATE( State_Chm%BoundaryCond, STAT=RC )
       CALL GC_CheckVar( 'State_Chm%BoundaryCond', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Chm%BoundaryCond => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Chm%Hg_Cat_Name ) ) THEN
       DEALLOCATE( State_Chm%Hg_Cat_Name, STAT=RC )
       CALL GC_CheckVar( 'State_Chm%Hg_Cat_Name', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Chm%Hg_Cat_Name => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Chm%Hg0_Id_List ) ) THEN
       DEALLOCATE( State_Chm%Hg0_Id_List, STAT=RC )
       CALL GC_CheckVar( 'State_Chm%Hg0_Id_List', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Chm%Hg0_Id_List => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Chm%Hg2_Id_List ) ) THEN
       DEALLOCATE( State_Chm%Hg2_Id_List, STAT=RC )
       CALL GC_CheckVar( 'State_Chm%Hg2_Id_List', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Chm%Hg2_Id_List => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Chm%HgP_Id_List ) ) THEN
       DEALLOCATE( State_Chm%HgP_Id_List, STAT=RC )
       CALL GC_CheckVar( 'State_Chm%HgP_Id_List', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Chm%HgP_Id_List => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Chm%AeroArea ) ) THEN
       DEALLOCATE( State_Chm%AeroArea, STAT=RC )
       CALL GC_CheckVar( 'State_Chm%AeroArea', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Chm%AeroArea => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Chm%AeroRadi ) ) THEN
       DEALLOCATE( State_Chm%AeroRadi, STAT=RC )
       CALL GC_CheckVar( 'State_Chm%AeroRadi', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Chm%AeroRadi => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Chm%WetAeroArea ) ) THEN
       DEALLOCATE( State_Chm%WetAeroArea, STAT=RC )
       CALL GC_CheckVar( 'State_Chm%WetAeroArea', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Chm%WetAeroArea => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Chm%WetAeroRadi ) ) THEN
       DEALLOCATE( State_Chm%WetAeroRadi, STAT=RC )
       CALL GC_CheckVar( 'State_Chm%WetAeroRadi', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Chm%WetAeroRadi => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Chm%AeroH2O ) ) THEN
       DEALLOCATE( State_Chm%AeroH2O, STAT=RC )
       CALL GC_CheckVar( 'State_Chm%AeroH2O', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Chm%AeroH2O => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Chm%GammaN2O5 ) ) THEN
       DEALLOCATE( State_Chm%GammaN2O5, STAT=RC )
       CALL GC_CheckVar( 'State_Chm%GammaN2O5', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Chm%GammaN2O5 => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Chm%OMOC_POA ) ) THEN
       DEALLOCATE( State_Chm%OMOC_POA, STAT=RC )
       CALL GC_CheckVar( 'State_Chm%OMOC_POA', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Chm%OMOC_POA => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Chm%OMOC_OPOA ) ) THEN
       DEALLOCATE( State_Chm%OMOC_OPOA, STAT=RC )
       CALL GC_CheckVar( 'State_Chm%OMOC_OPOA', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Chm%OMOC_OPOA => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Chm%phSav ) ) THEN
       DEALLOCATE( State_Chm%phSav, STAT=RC  )
       CALL GC_CheckVar( 'State_Chm%phSav', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Chm%pHSav => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Chm%HplusSav ) ) THEN
       DEALLOCATE( State_Chm%HplusSav, STAT=RC )
       CALL GC_CheckVar( 'State_Chm%HplusSav', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Chm%HplusSav => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Chm%WaterSav ) ) THEN
       DEALLOCATE( State_Chm%WaterSav, STAT=RC )
       CALL GC_CheckVar( 'State_Chm%WaterSav', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Chm%WaterSav => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Chm%SulRatSav ) ) THEN
       DEALLOCATE( State_Chm%SulRatSav, STAT=RC )
       CALL GC_CheckVar( 'State_Chm%SulRatSav', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Chm%SulRatSav => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Chm%NaRatSav ) ) THEN
       DEALLOCATE( State_Chm%NaRatSav, STAT=RC )
       CALL GC_CheckVar( 'State_Chm%NaRatSav', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Chm%NaRatSav => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Chm%AcidPurSav ) ) THEN
       DEALLOCATE( State_Chm%AcidPurSav, STAT=RC )
       CALL GC_CheckVar( 'State_Chm%AcidPurSav', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Chm%AcidPurSav => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Chm%BisulSav ) ) THEN
       DEALLOCATE( State_Chm%BisulSav, STAT=RC )
       CALL GC_CheckVar( 'State_Chm%BiSulSav', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Chm%BisulSav => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Chm%pHCloud ) ) THEN
       DEALLOCATE( State_Chm%pHCloud, STAT=RC )
       CALL GC_CheckVar( 'State_Chm%pHCloud', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Chm%pHCloud => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Chm%isCloud ) ) THEN
       DEALLOCATE( State_Chm%isCloud, STAT=RC )
       CALL GC_CheckVar( 'State_Chm%isCloud', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Chm%isCloud => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Chm%SSAlk ) ) THEN
       DEALLOCATE( State_Chm%SSAlk, STAT=RC )
       CALL GC_CheckVar( 'State_Chm%SSAlk', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Chm%SSAlk => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Chm%H2O2AfterChem ) ) THEN
       DEALLOCATE( State_Chm%H2O2AfterChem, STAT=RC )
       CALL GC_CheckVar( 'State_Chm%H2O2AfterChem', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ASSOCIATED( State_Chm%SO2AfterChem ) ) THEN
       DEALLOCATE( State_Chm%SO2AfterChem, STAT=RC )
       CALL GC_CheckVar( 'State_Chm%SO2AfterChem', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ASSOCIATED( State_Chm%DryDepNitrogen ) ) THEN
       DEALLOCATE( State_Chm%DryDepNitrogen, STAT=RC )
       CALL GC_CheckVar( 'State_Chm%DryDepNitrogen', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ASSOCIATED( State_Chm%WetDepNitrogen ) ) THEN
       DEALLOCATE( State_Chm%WetDepNitrogen, STAT=RC )
       CALL GC_CheckVar( 'State_Chm%WetDepNitrogen', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ASSOCIATED( State_Chm%KPPHvalue ) ) THEN
       DEALLOCATE( State_Chm%KPPHvalue, STAT=RC  )
       CALL GC_CheckVar( 'State_Chm%KPPHvalue', 2, RC )
       RETURN
    ENDIF

    IF ( ASSOCIATED( State_Chm%STATE_PSC ) ) THEN
       DEALLOCATE( State_Chm%STATE_PSC, STAT=RC  )
       CALL GC_CheckVar( 'State_Chm%State_PSC', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Chm%STATE_PSC => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Chm%KHETI_SLA ) ) THEN
       DEALLOCATE( State_Chm%KHETI_SLA, STAT=RC  )
       CALL GC_CheckVar( 'State_Chm%KHETI_SLA', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Chm%KHETI_SLA => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Chm%HSO3_AQ ) ) THEN
       DEALLOCATE( State_Chm%HSO3_AQ, STAT=RC )
       CALL GC_CheckVar( 'State_Chm%HSO3_AQ', 3, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Chm%HSO3_AQ => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Chm%SO3_AQ ) ) THEN
       DEALLOCATE( State_Chm%SO3_AQ, STAT=RC )
       CALL GC_CheckVar( 'State_Chm%SO3_AQ', 3, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Chm%SO3_AQ => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Chm%fupdateHOBr ) ) THEN
       DEALLOCATE( State_Chm%fupdateHOBr, STAT=RC )
       CALL GC_CheckVar( 'State_Chm%fupdateHOBr', 3, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Chm%fupdateHOBr => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Chm%OceanHg0 ) ) THEN
       DEALLOCATE( State_Chm%OceanHg0, STAT=RC )
       CALL GC_CheckVar( 'State_Chm%OceanHg0', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ASSOCIATED( State_Chm%OceanHg2 ) ) THEN
       DEALLOCATE( State_Chm%OceanHg2, STAT=RC )
       CALL GC_CheckVar( 'State_Chm%OceanHg2', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ASSOCIATED( State_Chm%OceanHgP ) ) THEN
       DEALLOCATE( State_Chm%OceanHgP, STAT=RC )
       CALL GC_CheckVar( 'State_Chm%OceanHgP', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ASSOCIATED( State_Chm%SnowHgOcean ) ) THEN
       DEALLOCATE( State_Chm%SnowHgOcean, STAT=RC )
       CALL GC_CheckVar( 'State_Chm%SnowHgOcean', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ASSOCIATED( State_Chm%SnowHgLand ) ) THEN
       DEALLOCATE( State_Chm%SnowHgLand, STAT=RC )
       CALL GC_CheckVar( 'State_Chm%SnowHgLand', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ASSOCIATED( State_Chm%SnowHgOceanStored ) ) THEN
       DEALLOCATE( State_Chm%SnowHgOceanStored, STAT=RC )
       CALL GC_CheckVar( 'State_Chm%SnowHgOceanStored', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ASSOCIATED( State_Chm%SnowHgLandStored ) ) THEN
       DEALLOCATE( State_Chm%SnowHgLandStored, STAT=RC )
       CALL GC_CheckVar( 'State_Chm%SnowHgLandStored', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

#if defined( MODEL_GEOS )
    ! Aerodynamic resistance @ 2m
    IF ( ASSOCIATED( State_Chm%DryDepRa2m ) ) THEN
       DEALLOCATE( State_Chm%DryDepRa2m, STAT=RC )
       CALL GC_CheckVar( 'State_Chm%DryDepRa2m', 3, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Chm%DryDepRa2m => NULL()
    ENDIF

    ! Aerodynamic resistance @ 10m
    IF ( ASSOCIATED( State_Chm%DryDepRa10m ) ) THEN
       DEALLOCATE( State_Chm%DryDepRa10m, STAT=RC )
       CALL GC_CheckVar( 'State_Chm%DryDepRa10m', 3, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Chm%DryDepRa10m => NULL()
    ENDIF
#endif

    IF ( ASSOCIATED( State_Chm%DryDepSav ) ) THEN
       DEALLOCATE( State_Chm%DryDepSav, STAT=RC )
       CALL GC_CheckVar( 'State_Chm%DryDepSav', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Chm%DryDepSav => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Chm%TLSTT ) ) THEN
       DEALLOCATE( State_Chm%TLSTT, STAT=RC )
       CALL GC_CheckVar( 'State_Chm%TLSTT', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Chm%TLSTT => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Chm%PSO4s ) ) THEN
       DEALLOCATE( State_Chm%PSO4s, STAT=RC )
       CALL GC_CheckVar( 'State_Chm%PSO4s', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Chm%PSO4s => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Chm%QQ3D ) ) THEN
       DEALLOCATE( State_Chm%QQ3D, STAT=RC )
       CALL GC_CheckVar( 'State_Chm%QQ3D', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Chm%QQ3D => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Chm%SurfaceFlux ) ) THEN
       DEALLOCATE( State_Chm%SurfaceFlux, STAT=RC )
       CALL GC_CheckVar( 'State_Chm%SurfaceFlux', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Chm%SurfaceFlux => NULL()
    ENDIF

    !-----------------------------------------------------------------------
    ! Template for deallocating more arrays, replace xxx with field name
    !-----------------------------------------------------------------------
    !IF ( ASSOCIATED( State_Chm%xxx ) ) THEN
    !   DEALLOCATE( State_Chm%xxx, STAT=RC  )
    !   CALL GC_CheckVar( 'State_Chm%xxx', 2, RC )
    !   IF ( RC /= GC_SUCCESS ) RETURN
    !   State_Chm%xxx => NULL()
    !ENDIF

    !=======================================================================
    ! Deallocate the species database object field
    !=======================================================================

    ! This operation should ONLY be done if there are no remaining chemistry
    ! states in the system, as destroying this State_Chm%SpcData will destroy
    ! all %SpcDatas, incl. state_chm_mod.F90's SpcDataLocal, in this CPU.
    !
    ! The variable state_chm_mod.F90 nChmState keeps track of the # of chemistry
    ! states initialized in the system. (hplin, 8/3/18)
    IF ( nChmState == 1 ) THEN
       CALL Cleanup_Species_Database( State_Chm%SpcData, RC )
       CALL GC_CheckVar( 'State_Chm%SpcData', 3, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    ! Nullify the State_Chm%SpcData object
    State_Chm%SpcData => NULL()

    !=======================================================================
    ! Destroy the registry of fields for this module
    !=======================================================================
    CALL Registry_Destroy( State_Chm%Registry, State_Chm%RegDict, RC )
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Could not destroy registry object State_Chm%Registry!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    ! Nullify the registry object
    State_Chm%Registry => NULL()

    !=======================================================================
    ! Decrease the counter of chemistry states in this CPU
    !=======================================================================
    nChmState = nChmState - 1

  END SUBROUTINE Cleanup_State_Chm
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get_Metadata_State_Chm
!
! !DESCRIPTION: Subroutine GET\_METADATA\_STATE\_CHM retrieves basic
!  information about each State\_Chm field.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Get_Metadata_State_Chm( am_I_Root,  metadataID, Found,      &
                                     RC,         Desc,       Units,      &
                                     PerSpecies, Rank,       Type,       &
                                     VLoc )
!
! !USES:
!
    USE Charpak_Mod,        ONLY : To_UpperCase
    USE Registry_Params_Mod
!
! !INPUT PARAMETERS:
!
    LOGICAL,             INTENT(IN)  :: am_I_Root
    CHARACTER(LEN=*),    INTENT(IN)  :: metadataID  ! State_Chm field name
!
! !OUTPUT PARAMETERS:
!
    LOGICAL,             INTENT(OUT)           :: Found      ! Item found?
    INTEGER,             INTENT(OUT)           :: RC         ! Return code
    CHARACTER(LEN=255),  INTENT(OUT), OPTIONAL :: Desc       ! Long name string
    CHARACTER(LEN=255),  INTENT(OUT), OPTIONAL :: Units      ! Units string
    CHARACTER(LEN=255),  INTENT(OUT), OPTIONAL :: PerSpecies ! Max spc wildcard
    INTEGER,             INTENT(OUT), OPTIONAL :: Rank       ! # of dimensions
    INTEGER,             INTENT(OUT), OPTIONAL :: Type       ! Desc of data type
    INTEGER,             INTENT(OUT), OPTIONAL :: VLoc       ! Vert placement
!
! !REMARKS:
!
! !REVISION HISTORY:
!  02 Oct 2017 - E. Lundgren - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    CHARACTER(LEN=255) :: ErrMsg, ThisLoc, Name_AllCaps
    LOGICAL            :: isDesc, isUnits, isRank, isType, isVLoc, isSpecies

    !=======================================================================
    ! Initialize
    !=======================================================================

    ! Assume success
    RC    =  GC_SUCCESS
    ThisLoc = ' -> at Get_Metadata_State_Chm (in Headers/state_chm_mod.F90)'
    Found = .TRUE.

    ! Optional arguments present?
    isDesc    = PRESENT( Desc  )
    isUnits   = PRESENT( Units )
    isRank    = PRESENT( Rank  )
    isType    = PRESENT( Type  )
    isVLoc    = PRESENT( VLoc  )
    isSpecies = PRESENT( PerSpecies )

    ! Set defaults for optional arguments. Assume type and vertical
    ! location are real (flexible precision) and center unless specified
    ! otherwise
    IF ( isUnits ) Units = ''
    IF ( isDesc  ) Desc  = ''
    IF ( isRank  ) Rank  = -1              ! Initialize # dims as bad value
    IF ( isType  ) Type  = KINDVAL_FP      ! Assume real(fp) for State_Chm flds
    IF ( isVLoc  ) VLoc  = VLocationCenter ! Assume vertically centered
    IF ( isSpecies ) PerSpecies = ''       ! Assume not per species

    ! Convert name to uppercase
    Name_AllCaps = To_Uppercase( TRIM( metadataID ) )

    !=======================================================================
    ! Values for Retrieval (string comparison slow but happens only once)
    !=======================================================================
    SELECT CASE ( TRIM( Name_AllCaps ) )

       CASE ( 'SPECIES' )
          IF ( isDesc    ) Desc  = 'Concentration for species'
          IF ( isUnits   ) Units = 'varies'
          IF ( isRank    ) Rank  = 3
          IF ( isSpecies ) PerSpecies = 'ALL'

       CASE( 'BOUNDARYCOND' )
          IF ( isDesc    ) Desc  = 'Boundary conditions for species'
          IF ( isUnits   ) Units = 'v/v'
          IF ( isRank    ) Rank  = 3
          IF ( isSpecies ) PerSpecies = 'ALL'

       CASE ( 'AEROAREAMDUST1' )
          IF ( isDesc  ) Desc  = 'Dry aerosol area for mineral dust (0.15 um)'
          IF ( isUnits ) Units = 'cm2 cm-3'
          IF ( isRank  ) Rank  = 3

       CASE ( 'AEROAREAMDUST2' )
          IF ( isDesc  ) Desc  = 'Dry aerosol area for mineral dust (0.25 um)'
          IF ( isUnits ) Units = 'cm2 cm-3'
          IF ( isRank  ) Rank  = 3

       CASE ( 'AEROAREAMDUST3' )
          IF ( isDesc  ) Desc  = 'Dry aerosol area for mineral dust (0.4 um)'
          IF ( isUnits ) Units = 'cm2 cm-3'
          IF ( isRank  ) Rank  = 3

       CASE ( 'AEROAREAMDUST4' )
          IF ( isDesc  ) Desc  = 'Dry aerosol area for mineral dust (0.8 um)'
          IF ( isUnits ) Units = 'cm2 cm-3'
          IF ( isRank  ) Rank  = 3

       CASE ( 'AEROAREAMDUST5' )
          IF ( isDesc  ) Desc  = 'Dry aerosol area for mineral dust (1.5 um)'
          IF ( isUnits ) Units = 'cm2 cm-3'
          IF ( isRank  ) Rank  = 3

       CASE ( 'AEROAREAMDUST6' )
          IF ( isDesc  ) Desc  = 'Dry aerosol area for mineral dust (2.5 um)'
          IF ( isUnits ) Units = 'cm2 cm-3'
          IF ( isRank  ) Rank  = 3

       CASE ( 'AEROAREAMDUST7' )
          IF ( isDesc  ) Desc  = 'Dry aerosol area for mineral dust (4.0 um)'
          IF ( isUnits ) Units = 'cm2 cm-3'
          IF ( isRank  ) Rank  = 3

       CASE ( 'AEROAREASULF' )
          IF ( isDesc  ) Desc  = 'Dry aerosol area for black carbon'
          IF ( isUnits ) Units = 'cm2 cm-3'
          IF ( isRank  ) Rank  = 3

       CASE ( 'AEROAREABC' )
          IF ( isDesc  ) Desc  = 'Dry aerosol area for black carbon'
          IF ( isUnits ) Units = 'cm2 cm-3'
          IF ( isRank  ) Rank  = 3

       CASE ( 'AEROAREAOC' )
          IF ( isDesc  ) Desc  = 'Dry aerosol area for organic carbon'
          IF ( isUnits ) Units = 'cm2 cm-3'
          IF ( isRank  ) Rank  = 3

       CASE ( 'AEROAREASSA' )
          IF ( isDesc  ) Desc  = 'Dry aerosol area for sea salt,' &
                                 // ' accumulation mode'
          IF ( isUnits ) Units = 'cm2 cm-3'
          IF ( isRank  ) Rank  = 3

       CASE ( 'AEROAREASSC' )
          IF ( isDesc  ) Desc  = 'Dry aerosol area for sea salt, coarse mode'
          IF ( isUnits ) Units = 'cm2 cm-3'
          IF ( isRank  ) Rank  = 3

       CASE ( 'AEROAREABGSULF' )
          IF ( isDesc  ) Desc  = 'Dry aerosol area for background' &
                                 // ' stratospheric sulfate'
          IF ( isUnits ) Units = 'cm2 cm-3'
          IF ( isRank  ) Rank  = 3

       CASE ( 'AEROAREAICEI' )
          IF ( isDesc  ) Desc  = 'Dry aerosol area for irregular ice cloud' &
                                 // ' (Mischenko)'
          IF ( isUnits ) Units = 'cm2 cm-3'
          IF ( isRank  ) Rank  = 3

       CASE ( 'AERORADIMDUST1' )
          IF ( isDesc  ) Desc  = 'Dry aerosol radius for mineral dust (0.15 um)'
          IF ( isUnits ) Units = 'cm'
          IF ( isRank  ) Rank  = 3

       CASE ( 'AERORADIMDUST2' )
          IF ( isDesc  ) Desc  = 'Dry aerosol radius for mineral dust (0.25 um)'
          IF ( isUnits ) Units = 'cm'
          IF ( isRank  ) Rank  = 3

       CASE ( 'AERORADIMDUST3' )
          IF ( isDesc  ) Desc  = 'Dry aerosol radius for mineral dust (0.4 um)'
          IF ( isUnits ) Units = 'cm'
          IF ( isRank  ) Rank  = 3

       CASE ( 'AERORADIMDUST4' )
          IF ( isDesc  ) Desc  = 'Dry aerosol radius for mineral dust (0.8 um)'
          IF ( isUnits ) Units = 'cm'
          IF ( isRank  ) Rank  = 3

       CASE ( 'AERORADIMDUST5' )
          IF ( isDesc  ) Desc  = 'Dry aerosol radius for mineral dust (1.5 um)'
          IF ( isUnits ) Units = 'cm'
          IF ( isRank  ) Rank  = 3

       CASE ( 'AERORADIMDUST6' )
          IF ( isDesc  ) Desc  = 'Dry aerosol radius for mineral dust (2.5 um)'
          IF ( isUnits ) Units = 'cm'
          IF ( isRank  ) Rank  = 3

       CASE ( 'AERORADIMDUST7' )
          IF ( isDesc  ) Desc  = 'Dry aerosol radius for mineral dust (4.0 um)'
          IF ( isUnits ) Units = 'cm'
          IF ( isRank  ) Rank  = 3

       CASE ( 'AERORADISULF' )
          IF ( isDesc  ) Desc  = 'Dry aerosol radius for tropospheric sulfate'
          IF ( isUnits ) Units = 'cm'
          IF ( isRank  ) Rank  = 3

       CASE ( 'AERORADIBC' )
          IF ( isDesc  ) Desc  = 'Dry aerosol radius for black carbon'
          IF ( isUnits ) Units = 'cm'
          IF ( isRank  ) Rank  = 3

       CASE ( 'AERORADIOC' )
          IF ( isDesc  ) Desc  = 'Dry aerosol radius for organic carbon'
          IF ( isUnits ) Units = 'cm'
          IF ( isRank  ) Rank  = 3

       CASE ( 'AERORADISSA' )
          IF ( isDesc  ) Desc  = 'Dry aerosol radius for sea salt,' &
                                 // ' accumulation mode'
          IF ( isUnits ) Units = 'cm'
          IF ( isRank  ) Rank  = 3

       CASE ( 'AERORADISSC' )
          IF ( isDesc  ) Desc  = 'Dry aerosol radius for sea salt, coarse mode'
          IF ( isUnits ) Units = 'cm'
          IF ( isRank  ) Rank  = 3

       CASE ( 'AERORADIBGSULF' )
          IF ( isDesc  ) Desc  = 'Dry aerosol radius for background' &
                                 // ' stratospheric sulfate'
          IF ( isUnits ) Units = 'cm'
          IF ( isRank  ) Rank  = 3

       CASE ( 'AERORADIICEI' )
          IF ( isDesc  ) Desc  = 'Dry aerosol radius for irregular ice' &
                                 // ' cloud (Mischenko)'
          IF ( isUnits ) Units = 'cm'
          IF ( isRank  ) Rank  = 3

       CASE ( 'WETAEROAREAMDUST1' )
          IF ( isDesc  ) Desc  = 'Wet aerosol area for mineral dust (0.15 um)'
          IF ( isUnits ) Units = 'cm2 cm-3'
          IF ( isRank  ) Rank  = 3

       CASE ( 'WETAEROAREAMDUST2' )
          IF ( isDesc  ) Desc  = 'Wet aerosol area for mineral dust (0.25 um)'
          IF ( isUnits ) Units = 'cm2 cm-3'
          IF ( isRank  ) Rank  = 3

       CASE ( 'WETAEROAREAMDUST3' )
          IF ( isDesc  ) Desc  = 'Wet aerosol area for mineral dust (0.4 um)'
          IF ( isUnits ) Units = 'cm2 cm-3'
          IF ( isRank  ) Rank  = 3

       CASE ( 'WETAEROAREAMDUST4' )
          IF ( isDesc  ) Desc  = 'Wet aerosol area for mineral dust (0.8 um)'
          IF ( isUnits ) Units = 'cm2 cm-3'
          IF ( isRank  ) Rank  = 3

       CASE ( 'WETAEROAREAMDUST5' )
          IF ( isDesc  ) Desc  = 'Wet aerosol area for mineral dust (1.5 um)'
          IF ( isUnits ) Units = 'cm2 cm-3'
          IF ( isRank  ) Rank  = 3

       CASE ( 'WETAEROAREAMDUST6' )
          IF ( isDesc  ) Desc  = 'Wet aerosol area for mineral dust (2.5 um)'
          IF ( isUnits ) Units = 'cm2 cm-3'
          IF ( isRank  ) Rank  = 3

       CASE ( 'WETAEROAREAMDUST7' )
          IF ( isDesc  ) Desc  = 'Wet aerosol area for mineral dust (4.0 um)'
          IF ( isUnits ) Units = 'cm2 cm-3'
          IF ( isRank  ) Rank  = 3

       CASE ( 'WETAEROAREASULF' )
          IF ( isDesc  ) Desc  = 'Wet aerosol area for tropospheric sulfate'
          IF ( isUnits ) Units = 'cm2 cm-3'
          IF ( isRank  ) Rank  = 3

       CASE ( 'WETAEROAREABC' )
          IF ( isDesc  ) Desc  = 'Wet aerosol area for black carbon'
          IF ( isUnits ) Units = 'cm2 cm-3'
          IF ( isRank  ) Rank  = 3

       CASE ( 'WETAEROAREAOC' )
          IF ( isDesc  ) Desc  = 'Wet aerosol area for organic carbon'
          IF ( isUnits ) Units = 'cm2 cm-3'
          IF ( isRank  ) Rank  = 3

       CASE ( 'WETAEROAREASSA' )
          IF ( isDesc  ) Desc  = 'Wet aerosol area for sea salt,' &
                                 // ' accumulation mode'
          IF ( isUnits ) Units = 'cm2 cm-3'
          IF ( isRank  ) Rank  = 3

       CASE ( 'WETAEROAREASSC' )
          IF ( isDesc  ) Desc  = 'Wet aerosol area for sea salt, coarse mode'
          IF ( isUnits ) Units = 'cm2 cm-3'
          IF ( isRank  ) Rank  = 3

       CASE ( 'WETAEROAREABGSULF' )
          IF ( isDesc  ) Desc  = 'Wet aerosol area for background' &
                                 // ' stratospheric sulfate'
          IF ( isUnits ) Units = 'cm2 cm-3'
          IF ( isRank  ) Rank  = 3

       CASE ( 'WETAEROAREAICEI' )
          IF ( isDesc  ) Desc  = 'Wet aerosol area for irregular ice cloud' &
                                 // ' (Mischenko)'
          IF ( isUnits ) Units = 'cm2 cm-3'
          IF ( isRank  ) Rank  = 3

       CASE ( 'WETAERORADIMDUST1' )
          IF ( isDesc  ) Desc  = 'Wet aerosol radius for mineral dust (0.15 um)'
          IF ( isUnits ) Units = 'cm'
          IF ( isRank  ) Rank  = 3

       CASE ( 'WETAERORADIMDUST2' )
          IF ( isDesc  ) Desc  = 'Wet aerosol radius for mineral dust (0.25 um)'
          IF ( isUnits ) Units = 'cm'
          IF ( isRank  ) Rank  = 3

       CASE ( 'WETAERORADIMDUST3' )
          IF ( isDesc  ) Desc  = 'Wet aerosol radius for mineral dust (0.4 um)'
          IF ( isUnits ) Units = 'cm'
          IF ( isRank  ) Rank  = 3

       CASE ( 'WETAERORADIMDUST4' )
          IF ( isDesc  ) Desc  = 'Wet aerosol radius for mineral dust (0.8 um)'
          IF ( isUnits ) Units = 'cm'
          IF ( isRank  ) Rank  = 3

       CASE ( 'WETAERORADIMDUST5' )
          IF ( isDesc  ) Desc  = 'Wet aerosol radius for mineral dust (1.5 um)'
          IF ( isUnits ) Units = 'cm'
          IF ( isRank  ) Rank  = 3

       CASE ( 'WETAERORADIMDUST6' )
          IF ( isDesc  ) Desc  = 'Wet aerosol radius for mineral dust (2.5 um)'
          IF ( isUnits ) Units = 'cm'
          IF ( isRank  ) Rank  = 3

       CASE ( 'WETAERORADIMDUST7' )
          IF ( isDesc  ) Desc  = 'Wet aerosol radius for mineral dust (4.0 um)'
          IF ( isUnits ) Units = 'cm'
          IF ( isRank  ) Rank  = 3

       CASE ( 'WETAERORADISULF' )
          IF ( isDesc  ) Desc  = 'Wet aerosol radius for tropospheric sulfate'
          IF ( isUnits ) Units = 'cm'
          IF ( isRank  ) Rank  = 3

       CASE ( 'WETAERORADIBC' )
          IF ( isDesc  ) Desc  = 'Wet aerosol radius for black carbon'
          IF ( isUnits ) Units = 'cm'
          IF ( isRank  ) Rank  = 3

       CASE ( 'WETAERORADIOC' )
          IF ( isDesc  ) Desc  = 'Wet aerosol radius for organic carbon'
          IF ( isUnits ) Units = 'cm'
          IF ( isRank  ) Rank  = 3

       CASE ( 'WETAERORADISSA' )
          IF ( isDesc  ) Desc= 'Wet aerosol radius for sea salt,' &
                               // ' accumulation mode'
          IF ( isUnits ) Units = 'cm'
          IF ( isRank  ) Rank  = 3

       CASE ( 'WETAERORADISSC' )
          IF ( isDesc  ) Desc  = 'Wet aerosol radius for sea salt, coarse mode'
          IF ( isUnits ) Units = 'cm'
          IF ( isRank  ) Rank  = 3

       CASE ( 'WETAERORADIBGSULF' )
          IF ( isDesc  ) Desc  = 'Wet aerosol radius for background' &
                                // ' stratospheric sulfate'
          IF ( isUnits ) Units = 'cm'
          IF ( isRank  ) Rank  = 3

       CASE ( 'WETAERORADIICEI' )
          IF ( isDesc  ) Desc  = 'Wet aerosol radius for irregular ice cloud' &
                                // ' (Mischenko)'
          IF ( isUnits ) Units = 'cm'
          IF ( isRank  ) Rank  = 3

       CASE ( 'AEROH2OMDUST1' )
          IF ( isDesc  ) Desc  = 'Aerosol H2O content for mineral dust (0.15 um)'
          IF ( isUnits ) Units = 'cm3(H2O) cm-3(air)'
          IF ( isRank  ) Rank  = 3

       CASE ( 'AEROH2OMDUST2' )
          IF ( isDesc  ) Desc  = 'Aerosol H2O content for mineral dust (0.25 um)'
          IF ( isUnits ) Units = 'cm3(H2O) cm-3(air)'
          IF ( isRank  ) Rank  = 3

       CASE ( 'AEROH2OMDUST3' )
          IF ( isDesc  ) Desc  = 'Aerosol H2O content for mineral dust (0.4 um)'
          IF ( isUnits ) Units = 'cm3(H2O) cm-3(air)'
          IF ( isRank  ) Rank  = 3

       CASE ( 'AEROH2OMDUST4' )
          IF ( isDesc  ) Desc  = 'Aerosol H2O content for mineral dust (0.8 um)'
          IF ( isUnits ) Units = 'cm3(H2O) cm-3(air)'
          IF ( isRank  ) Rank  = 3

       CASE ( 'AEROH2OMDUST5' )
          IF ( isDesc  ) Desc  = 'Aerosol H2O content for mineral dust (1.5 um)'
          IF ( isUnits ) Units = 'cm3(H2O) cm-3(air)'
          IF ( isRank  ) Rank  = 3

       CASE ( 'AEROH2OMDUST6' )
          IF ( isDesc  ) Desc  = 'Aerosol H2O content for mineral dust (2.5 um)'
          IF ( isUnits ) Units = 'cm3(H2O) cm-3(air)'
          IF ( isRank  ) Rank  = 3

       CASE ( 'AEROH2OMDUST7' )
          IF ( isDesc  ) Desc  = 'Aerosol H2O content for mineral dust (4.0 um)'
          IF ( isUnits ) Units = 'cm3(H2O) cm-3(air)'
          IF ( isRank  ) Rank  = 3

       CASE ( 'AEROH2OSULF' )
          IF ( isDesc  ) Desc  = 'Aerosol H2O content for tropospheric sulfate'
          IF ( isUnits ) Units = 'cm3(H2O) cm-3(air)'
          IF ( isRank  ) Rank  = 3

       CASE ( 'AEROH2OBC' )
          IF ( isDesc  ) Desc  = 'Aerosol H2O content for black carbon'
          IF ( isUnits ) Units = 'cm3(H2O) cm-3(air)'
          IF ( isRank  ) Rank  = 3

       CASE ( 'AEROH2OOC' )
          IF ( isDesc  ) Desc  = 'Aerosol H2O content for organic carbon'
          IF ( isUnits ) Units = 'cm3(H2O) cm-3(air)'
          IF ( isRank  ) Rank  = 3

       CASE ( 'AEROH2OSSA' )
          IF ( isDesc  ) Desc= 'Aerosol H2O content for sea salt,' &
                               // ' accumulation mode'
          IF ( isUnits ) Units = 'cm3(H2O) cm-3(air)'
          IF ( isRank  ) Rank  = 3

       CASE ( 'AEROH2OSSC' )
          IF ( isDesc  ) Desc  = 'Aerosol H2O content for sea salt, coarse mode'
          IF ( isUnits ) Units = 'cm3(H2O) cm-3(air)'
          IF ( isRank  ) Rank  = 3

       CASE ( 'AEROH2OBGSULF' )
          IF ( isDesc  ) Desc  = 'Aerosol H2O content for background' &
                                // ' stratospheric sulfate'
          IF ( isUnits ) Units = 'cm3(H2O) cm-3(air)'
          IF ( isRank  ) Rank  = 3

       CASE ( 'AEROH2OICEI' )
          IF ( isDesc  ) Desc  = 'Aerosol H2O content for irregular ice cloud' &
                                // ' (Mischenko)'
          IF ( isUnits ) Units = 'cm3(H2O) cm-3(air)'
          IF ( isRank  ) Rank  = 3


       CASE ( 'GAMMAN2O5H2O' )
          IF ( isDesc  ) Desc  = 'Sticking coefficient for N2O5 + H2O reaction'
          IF ( isUnits ) Units = '1'
          IF ( isRank  ) Rank  = 3

       CASE ( 'GAMMAN2O5HCL' )
          IF ( isDesc  ) Desc  = 'Sticking coefficient for N2O5 + HCl reaction'
          IF ( isUnits ) Units = '1'
          IF ( isRank  ) Rank  = 3

       CASE ( 'GAMMAN2O5SS' )
          IF ( isDesc  ) Desc  = 'Sticking coefficient for N2O5 + SS reaction'
          IF ( isUnits ) Units = '1'
          IF ( isRank  ) Rank  = 3

       CASE ( 'YIELDCLNO2' )
          IF ( isDesc  ) Desc  = 'Production yield coefficient for ClNO2' &
                                 // ' from N2O5 aerosol uptake'
          IF ( isUnits ) Units = '1'
          IF ( isRank  ) Rank  = 3

       CASE ( 'KPPHVALUE' )
          IF ( isDesc  ) Desc  = 'H-value for Rosenbrock solver'
          IF ( isUnits ) Units = '1'
          IF ( isRank  ) Rank  = 3

       CASE ( 'OMOCPOA' )
          IF ( isDesc  ) Desc  = 'OM:OC ratio for POA (from /aerosol_mod.F90)'
          IF ( isUnits ) Units = '1'
          IF ( isRank  ) Rank  = 2

       CASE ( 'OMOCOPOA' )
          IF ( isDesc  ) Desc  = 'OM:OC ratio for OPOA (from /aerosol_mod.F90)'
          IF ( isUnits ) Units = '1'
          IF ( isRank  ) Rank  = 2

       CASE ( 'STATEPSC' )
          IF ( isDesc  ) Desc  = 'Polar stratospheric cloud type (cf Kirner' &
                                // ' et al 2011, GMD)'
          IF ( isUnits ) Units = 'count'
          IF ( isRank  ) Rank  = 3

       CASE ( 'KHETISLAN2O5H2O' )
          IF ( isDesc  ) Desc  = 'Sticking coefficient for N2O5 + H2O reaction'
          IF ( isUnits ) Units = '1'
          IF ( isRank  ) Rank  = 3

       CASE ( 'KHETISLAN2O5HCL' )
          IF ( isDesc  ) Desc  = 'Sticking coefficient for N2O5 + H2O reaction'
          IF ( isUnits ) Units = '1'
          IF ( isRank  ) Rank  = 3

       CASE ( 'KHETISLACLNO3H2O' )
          IF ( isDesc  ) Desc  = 'Sticking coefficient for ClNO3 + H2O reaction'
          IF ( isUnits ) Units = '1'
          IF ( isRank  ) Rank  = 3

       CASE ( 'KHETISLACLNO3HCL' )
          IF ( isDesc  ) Desc  = 'Sticking coefficient for ClNO3 + HCl reaction'
          IF ( isUnits ) Units = '1'
          IF ( isRank  ) Rank  = 3

       CASE ( 'KHETISLACLNO3HBR' )
          IF ( isDesc  ) Desc  = 'Sticking coefficient for ClNO3 + HBr reaction'
          IF ( isUnits ) Units = '1'
          IF ( isRank  ) Rank  = 3

       CASE ( 'KHETISLABRNO3H2O' )
          IF ( isDesc  ) Desc  = 'Sticking coefficient for BrNO3 + H2O reaction'
          IF ( isUnits ) Units = '1'
          IF ( isRank  ) Rank  = 3

       CASE ( 'KHETISLABRNO3HCL' )
          IF ( isDesc  ) Desc  = 'Sticking coefficient for BrNO3 + HCl reaction'
          IF ( isUnits ) Units = '1'
          IF ( isRank  ) Rank  = 3

       CASE ( 'KHETISLAHOCLHCL' )
          IF ( isDesc  ) Desc  = 'Sticking coefficient for HOCl + HCl reaction'
          IF ( isUnits ) Units = '1'
          IF ( isRank  ) Rank  = 3

       CASE ( 'KHETISLAHOCLHBR' )
          IF ( isDesc  ) Desc  = 'Sticking coefficient for HClr + HBr reaction'
          IF ( isUnits ) Units = '1'
          IF ( isRank  ) Rank  = 3

       CASE ( 'KHETISLAHOBRHCL' )
          IF ( isDesc  ) Desc  = 'Sticking coefficient for HOBr + HCl reaction'
          IF ( isUnits ) Units = '1'
          IF ( isRank  ) Rank  = 3

       CASE ( 'KHETISLAHOBRHBR' )
          IF ( isDesc  ) Desc  = 'Sticking coefficient for HOBr + HBr reaction'
          IF ( isUnits ) Units = '1'
          IF ( isRank  ) Rank  = 3

       CASE( 'PHSAV' )
          IF ( isDesc  ) Desc  = 'ISORROPIA aerosol pH'
          IF ( isUnits ) Units = '1'
          IF ( isRank  ) Rank  = 3

       CASE( 'HPLUSSAV' )
          IF ( isDesc  ) Desc  = 'ISORROPIA H+ concentration'
          IF ( isUnits ) Units = 'mol L-1'
          IF ( isRank  ) Rank  = 3

       CASE( 'WATERSAV' )
          IF ( isDesc  ) Desc  = 'ISORROPIA aerosol water concentration'
          IF ( isUnits ) Units = 'ug m-3'
          IF ( isRank  ) Rank  = 3

       CASE( 'SULRATSAV' )
          IF ( isDesc  ) Desc  = 'ISORROPIA sulfate concentration'
          IF ( isUnits ) Units = 'mol L-1'
          IF ( isRank  ) Rank  = 3

       CASE( 'NARATSAV' )
          IF ( isDesc  ) Desc  = 'ISORROPIA Na+ concentration'
          IF ( isUnits ) Units = 'mol L-1'
          IF ( isRank  ) Rank  = 3

       CASE( 'ACIDPURSAV' )
          IF ( isDesc  ) Desc  = 'ISORROPIA ACIDPUR'
          IF ( isUnits ) Units = 'mol L-1'
          IF ( isRank  ) Rank  = 3

       CASE( 'BISULSAV' )
          IF ( isDesc  ) Desc  = 'ISORROPIA Bisulfate (general acid)' &
                                 // ' concentration'
          IF ( isUnits ) Units = 'mol L-1'
          IF ( isRank  ) Rank  =  3

       CASE( 'PHCLOUD' )
          IF ( isDesc  ) Desc  = 'Cloud pH'
          IF ( isUnits ) Units = '1'
          IF ( isRank  ) Rank  =  3

       CASE( 'ISCLOUD' )
          IF ( isDesc  ) Desc  = 'Cloud presence'
          IF ( isUnits ) Units = '1'
          IF ( isRank  ) Rank  =  3

       CASE( 'SSALKACCUM' )
          IF ( isDesc  ) Desc  = 'Sea salt alkalinity, accumulation mode'
          IF ( isUnits ) Units = '1'
          IF ( isRank  ) Rank  =  3

       CASE( 'SSALKCOARSE' )
          IF ( isDesc  ) Desc  = 'Sea salt alkalinity, coarse mode'
          IF ( isUnits ) Units = '1'
          IF ( isRank  ) Rank  =  3

       CASE ( 'HSO3AQ' )
          IF ( isDesc  ) Desc  = 'Cloud bisulfite concentration'
          IF ( isUnits ) Units = 'mol L-1'
          IF ( isRank  ) Rank  =  3

       CASE ( 'SO3AQ' )
          IF ( isDesc  ) Desc  = 'Cloud sulfite concentration'
          IF ( isUnits ) Units = 'mol L-1'
          IF ( isRank  ) Rank  =  3

       CASE ( 'FUPDATEHOBR' )
          IF ( isDesc  ) Desc  = 'Correction factor for HOBr removal by SO2'
          IF ( isUnits ) Units = '1'
          IF ( isRank  ) Rank  =  3

       CASE ( 'H2O2AFTERCHEM' )
          IF ( isDesc  ) Desc  = 'H2O2 after sulfate chemistry'
          IF ( isUnits ) Units = 'mol mol-1'
          IF ( isRank  ) Rank  =  3

       CASE ( 'SO2AFTERCHEM' )
          IF ( isDesc  ) Desc  = 'SO2 after sulfate chemistry'
          IF ( isUnits ) Units = 'mol mol-1'
          IF ( isRank  ) Rank  =  3

       CASE ( 'DRYDEPNITROGEN' )
          IF ( isDesc  ) Desc  = 'Dry deposited nitrogen'
          IF ( isUnits ) Units = 'molec cm-2 s-1'
          IF ( isRank  ) Rank  =  2

       CASE ( 'WETDEPNITROGEN' )
          IF ( isDesc  ) Desc  = 'Wet deposited nitrogen'
          IF ( isUnits ) Units = 'molec cm-2 s-1'
          IF ( isRank  ) Rank  =  2

       CASE( 'OCEANHG0' )
          IF ( isDesc  ) Desc  = 'Hg(0) ocean mass'
          IF ( isUnits ) Units = 'kg'
          IF ( isRank  ) Rank  = 2
          IF ( isSpecies ) PerSpecies = 'HgCat'

       CASE( 'OCEANHG2' )
          IF ( isDesc  ) Desc  = 'Hg(II) ocean mass'
          IF ( isUnits ) Units = 'kg'
          IF ( isRank  ) Rank  = 2
          IF ( isSpecies ) PerSpecies = 'HgCat'

       CASE( 'OCEANHGP' )
          IF ( isDesc  ) Desc  = 'HgP ocean mass'
          IF ( isUnits ) Units = 'kg'
          IF ( isRank  ) Rank  = 2
          IF ( isSpecies ) PerSpecies = 'HgCat'

       CASE( 'SNOWHGOCEAN' )
          IF ( isDesc  ) Desc  = 'Reducible Hg snowpack on ocean'
          IF ( isUnits ) Units = 'kg'
          IF ( isRank  ) Rank  = 2
          IF ( isSpecies ) PerSpecies = 'HgCat'

       CASE( 'SNOWHGLAND' )
          IF ( isDesc  ) Desc  = 'Reducible Hg snowpack on land'
          IF ( isUnits ) Units = 'kg'
          IF ( isRank  ) Rank  = 2
          IF ( isSpecies ) PerSpecies = 'HgCat'

       CASE( 'SNOWHGOCEANSTORED' )
          IF ( isDesc  ) Desc  = 'Non-reducible Hg snowpack on ocean'
          IF ( isUnits ) Units = 'kg'
          IF ( isRank  ) Rank  = 2
          IF ( isSpecies ) PerSpecies = 'HgCat'

       CASE( 'SNOWHGLANDSTORED' )
          IF ( isDesc  ) Desc  = 'Non-reducible Hg snowpack on land'
          IF ( isUnits ) Units = 'kg'
          IF ( isRank  ) Rank  = 2
          IF ( isSpecies ) PerSpecies = 'HgCat'

       CASE( 'DRYDEPSAV' )
          IF ( isDesc  ) Desc  = 'Dry deposition frequencies'
          IF ( isUnits ) Units = 's-1'
          IF ( isRank  ) Rank  = 3

       CASE( 'TLSTT' )
          IF ( isDesc  ) Desc  = 'TLSTT'
          IF ( isUnits ) Units = ''
          IF ( isRank  ) Rank  = 4

       CASE( 'PSO4S' )
          IF ( isDesc  ) Desc  = 'PSO4s'
          IF ( isUnits ) Units = '1'
          IF ( isRank  ) Rank  = 3

       CASE( 'QQ3D' )
          IF ( isDesc  ) Desc  = 'Rate of new precipitation formation'
          IF ( isUnits ) Units = 'cm3 H2O cm-3 air'
          IF ( isRank  ) Rank  = 3

       CASE( 'SURFACEFLUX' )
          IF ( isDesc  ) Desc  = 'Surface flux (E-D) for non-local PBL mixing'
          IF ( isUnits ) Units = 'kg m-2 s-1'
          IF ( isRank  ) Rank  = 3

       CASE DEFAULT
          Found = .False.
          ErrMsg = 'Metadata not found for State_Chm field ' // &
                   TRIM( metadataID ) // ' when search for all caps name ' &
                   // TRIM( Name_AllCaps )
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          IF ( RC /= GC_SUCCESS ) RETURN

    END SELECT

   END SUBROUTINE Get_Metadata_State_Chm
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Register_ChmField_R4_3D
!
! !DESCRIPTION: Registers a 3-dimensional, 4-byte floating point array field
!  of the State\_Chm object.  This allows the diagnostic modules get a pointer
!  to the field by searching on the field name.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Register_ChmField_R4_3D( Input_Opt,  metadataID, Ptr2Data,  &
                                      State_Chm,  RC )
!
! !USES:
!
    USE Input_Opt_Mod,      ONLY : OptInput
    USE Registry_Params_Mod
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput),    INTENT(IN)    :: Input_Opt       ! Input Options object
    CHARACTER(LEN=*),  INTENT(IN)    :: metadataID      ! Name
    REAL(f4),          POINTER       :: Ptr2Data(:,:,:) ! pointer to data
    TYPE(ChmState),    INTENT(IN)    :: State_Chm       ! Obj for chem state
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,           INTENT(OUT)   :: RC              ! Success/failure
!
! !REMARKS:
!
! !REVISION HISTORY:
!  20 Sep 2017 - E. Lundgren - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    CHARACTER(LEN=512)     :: ErrMsg
    CHARACTER(LEN=255)     :: ErrMsg_reg,  ThisLoc
    CHARACTER(LEN=255)     :: desc, units, perSpecies
    CHARACTER(LEN=255)     :: thisSpcName, thisSpcDesc
    INTEGER                :: N, rank, type,  vloc
    LOGICAL                :: found, onEdges
    TYPE(Species), POINTER :: SpcInfo

    !-----------------------------------------------------------------------
    ! Initialize
    !-----------------------------------------------------------------------
    RC = GC_SUCCESS
    ThisLoc = ' -> at Register_ChmField_R4_3D (in Headers/state_chm_mod.F90)'
    ErrMsg  = ''
    ErrMsg_reg = 'Error encountered while registering State_Chm%'

    !-----------------------------------------------------------------------
    ! Get metadata
    !-----------------------------------------------------------------------
    CALL Get_Metadata_State_Chm( Input_Opt%amIRoot, metadataID,  Found,  RC, &
                                 desc=desc, units=units, rank=rank,          &
                                 type=type, vloc=vloc, perSpecies=perSpecies )

    ! Trap potential errors
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = TRIM( ErrMsg_reg ) // TRIM( MetadataID ) //                  &
                '; Abnormal exit from routine "Get_Metadata_State_Chm"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    !-----------------------------------------------------------------------
    ! If not tied to species then simply register the single field
    !-----------------------------------------------------------------------
    IF ( perSpecies == '' ) THEN

       ! Check that metadata consistent with data size
       IF ( rank /= 3 ) THEN
          ErrMsg = 'Data dims and metadata rank do not match for '           &
                   // TRIM(metadataID)
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

       ! Is the data placed on vertical edges?
       onEdges = ( vLoc == VLocationEdge )

       ! Add field to registry
       CALL Registry_AddField( Input_Opt    = Input_Opt,                     &
                               Registry     = State_Chm%Registry,            &
                               State        = State_Chm%State,               &
                               Variable     = metadataID,                    &
                               Description  = desc,                          &
                               Units        = units,                         &
                               Data3d_4     = Ptr2Data,                      &
                               OnLevelEdges = onEdges,                       &
                               RC           = RC                            )

       ! Trap potential errors
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = TRIM( ErrMsg_reg ) // TRIM( MetadataID ) //               &
                  '; Abnormal exit from routine "Registry_AddField"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

    !-----------------------------------------------------------------------
    ! Otherwise exit with error
    !-----------------------------------------------------------------------
    ELSE

       ! Error: cannot register field!
       ErrMsg = 'Handling of PerSpecies metadata ' // TRIM(perSpecies) //    &
                ' is not implemented for this combo of data type and size'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN

    ENDIF

  END SUBROUTINE Register_ChmField_R4_3D
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Register_ChmField_Rfp_2D
!
! !DESCRIPTION: Registers a 2-dimensional, flexible-precision array field
!  of the State\_Chm object.  This allows the diagnostic modules get a pointer
!  to the field by searching on the field name.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Register_ChmField_Rfp_2D( Input_Opt, metadataID, Ptr2Data,  &
                                       State_Chm, RC )
!
! !USES:
!
    USE Input_Opt_Mod,      ONLY : OptInput
    USE Registry_Params_Mod
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput),    INTENT(IN)    :: Input_Opt       ! Input Options object
    CHARACTER(LEN=*),  INTENT(IN)    :: metadataID      ! State_Chm field ID
    REAL(fp),          POINTER       :: Ptr2Data(:,:)   ! pointer to data
    TYPE(ChmState),    INTENT(IN)    :: State_Chm       ! Obj for chem state
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,           INTENT(OUT)   :: RC              ! Success/failure
!
! !REMARKS:
!
! !REVISION HISTORY:
!  20 Sep 2017 - E. Lundgren - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    CHARACTER(LEN=512)     :: ErrMsg
    CHARACTER(LEN=255)     :: ErrMsg_reg, ThisLoc
    CHARACTER(LEN=255)     :: desc, units, perSpecies
    CHARACTER(LEN=255)     :: thisSpcName, thisSpcDesc
    INTEGER                :: N, rank, type,  vloc
    LOGICAL                :: found, onEdges
    TYPE(Species), POINTER :: SpcInfo

    !-----------------------------------------------------------------------
    ! Initialize
    !-----------------------------------------------------------------------
    RC      = GC_SUCCESS
    ThisLoc = ' -> at Register_ChmField_Rfp_2D (in Headers/state_chm_mod.F90)'
    ErrMsg  = ''
    ErrMsg_reg = 'Error encountered while registering State_Chm%'

    !-----------------------------------------------------------------------
    ! Get metadata
    !-----------------------------------------------------------------------
    CALL Get_Metadata_State_Chm( Input_Opt%amIRoot, metadataID,  Found,  RC, &
                                 desc=desc, units=units, rank=rank,          &
                                 type=type, vloc=vloc, perSpecies=perSpecies )

    ! Trap potential errors
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = TRIM( ErrMsg_reg ) // TRIM( MetadataID ) //                  &
                '; Abnormal exit from routine "Get_Metadata_State_Chm"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    ! Is the data placed on vertical edges?
    onEdges = ( vLoc == VLocationEdge )

    !-----------------------------------------------------------------------
    ! If not tied to species then simply register the single field
    !-----------------------------------------------------------------------
    IF ( perSpecies == '' ) THEN

       ! Check that metadata consistent with data size
       IF ( rank /= 2 ) THEN
          ErrMsg = 'Data dims and metadata rank do not match for '           &
                   // TRIM(metadataID)
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

       ! Add field to registry
       CALL Registry_AddField( Input_Opt    = Input_Opt,                     &
                               Registry     = State_Chm%Registry,            &
                               State        = State_Chm%State,               &
                               Variable     = metadataID,                    &
                               Description  = desc,                          &
                               Units        = units,                         &
                               Data2d       = Ptr2Data,                      &
                               RC           = RC                            )

       ! Trap potential errors
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = TRIM( ErrMsg_reg ) // TRIM( MetadataID ) //               &
                  '; Abnormal exit from routine "Registry_AddField"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

    !-----------------------------------------------------------------------
    ! Otherwise exit with error
    !-----------------------------------------------------------------------
    ELSE

       ErrMsg = 'Handling of PerSpecies metadata ' // TRIM(perSpecies) //    &
                ' is not implemented for this combo of data type and size'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN

    ENDIF

  END SUBROUTINE Register_ChmField_Rfp_2D
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Register_ChmField_Rfp_3D
!
! !DESCRIPTION: Registers a 3-dimensional, flexible-precision array field
!  of the State\_Chm object.  This allows the diagnostic modules get a pointer
!  to the field by searching on the field name.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Register_ChmField_Rfp_3D( Input_Opt, metadataID, Ptr2Data,  &
                                       State_Chm, RC )
!
! !USES:
!
    USE Input_Opt_Mod,      ONLY : OptInput
    USE Registry_Params_Mod
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput),    INTENT(IN)    :: Input_Opt       ! Input Options object
    CHARACTER(LEN=*),  INTENT(IN)    :: metadataID      ! State_Chm field ID
    REAL(fp),          POINTER       :: Ptr2Data(:,:,:) ! pointer to data
    TYPE(ChmState),    INTENT(IN)    :: State_Chm       ! Obj for chem state
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,           INTENT(OUT)   :: RC              ! Success/failure
!
! !REMARKS:
!
! !REVISION HISTORY:
!  20 Sep 2017 - E. Lundgren - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    CHARACTER(LEN=512)     :: ErrMsg
    CHARACTER(LEN=255)     :: ErrMsg_reg, ThisLoc
    CHARACTER(LEN=255)     :: desc, units, perSpecies
    CHARACTER(LEN=255)     :: thisSpcName, thisSpcDesc
    INTEGER                :: N, rank, type,  vloc
    LOGICAL                :: found, onEdges
    TYPE(Species), POINTER :: SpcInfo

    !-----------------------------------------------------------------------
    ! Initialize
    !-----------------------------------------------------------------------
    RC      = GC_SUCCESS
    ThisLoc = ' -> at Register_ChmField_Rfp_3D (in Headers/state_chm_mod.F90)'
    ErrMsg  = ''
    ErrMsg_reg = 'Error encountered while registering State_Chm%'

    !-----------------------------------------------------------------------
    ! Get metadata
    !-----------------------------------------------------------------------
    CALL Get_Metadata_State_Chm( Input_Opt%amIRoot, metadataID,  Found,  RC, &
                                 desc=desc, units=units, rank=rank,          &
                                 type=type, vloc=vloc, perSpecies=perSpecies )

    ! Trap potential errors
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = TRIM( ErrMsg_reg ) // TRIM( MetadataID ) //                  &
                '; Abnormal exit from routine "Get_Metadata_State_Chm"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    ! Is the data placed on vertical edges?
    onEdges = ( vLoc == VLocationEdge )

    !-----------------------------------------------------------------------
    ! If not tied to species then simply register the single field
    !-----------------------------------------------------------------------
    IF ( perSpecies == '' ) THEN

       ! Check that metadata consistent with data size
       IF ( rank /= 3 ) THEN
          ErrMsg = 'Data dims and metadata rank do not match for '           &
                   // TRIM(metadataID)
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

       ! Add field to registry
       CALL Registry_AddField( Input_Opt    = Input_Opt,                     &
                               Registry     = State_Chm%Registry,            &
                               State        = State_Chm%State,               &
                               Variable     = metadataID,                    &
                               Description  = desc,                          &
                               Units        = units,                         &
                               OnLevelEdges = onEdges,                       &
                               Data3d       = Ptr2Data,                      &
                               RC           = RC                            )

       ! Trap potential errors
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = TRIM( ErrMsg_reg ) // TRIM( MetadataID ) //               &
                  '; Abnormal exit from routine "Registry_AddField"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

    !-----------------------------------------------------------------------
    ! If tied to Hg category then register each one
    !-----------------------------------------------------------------------
    ELSE IF ( perSpecies == 'HgCat' ) THEN

       ! Loop over all species
       DO N = 1, State_Chm%N_Hg_CATS

          ! Append Hg category to name and description for tagged Hg
          IF ( N == 1 ) THEN
             thisSpcName = TRIM( metadataID )
             thisSpcDesc = TRIM( Desc       )
          ELSE
             thisSpcName = TRIM( metadataID ) // '_' // &
                  TRIM(State_Chm%Hg_Cat_Name(N))
             thisSpcDesc = TRIM( Desc       ) // ' ' // &
                  TRIM(State_Chm%Hg_Cat_Name(N))
          ENDIF

          ! Add field to registry
          CALL Registry_AddField( Input_Opt    = Input_Opt,                  &
                                  Registry     = State_Chm%Registry ,        &
                                  State        = State_Chm%State,            &
                                  Variable     = thisSpcName,                &
                                  Description  = thisSpcDesc,                &
                                  Units        = units,                      &
                                  OnLevelEdges = onEdges,                    &
                                  Data2d       = Ptr2Data(:,:,N),            &
                                  RC           = RC                         )

          ! Free pointers
          SpcInfo => NULL()

          ! Trap potential errors
          IF ( RC /= GC_SUCCESS ) THEN
             ErrMsg = TRIM( ErrMsg_reg ) // TRIM( MetadataID ) //            &
                      '; Abnormal exit from routine "Registry_AddField"!'
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN
          ENDIF

       ENDDO

    !-----------------------------------------------------------------------
    ! Otherwise exit with error
    !-----------------------------------------------------------------------
    ELSE

       ErrMsg = 'Handling of PerSpecies metadata ' // TRIM(perSpecies) //    &
                ' is not implemented for this combo of data type and size'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN

    ENDIF

  END SUBROUTINE Register_ChmField_Rfp_3D
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Register_ChmField_Rfp_4D
!
! !DESCRIPTION: Registers a 4-dimensional, flexible-precision array field
!  of the State\_Chm object.  This allows the diagnostic modules get a pointer
!  to the field by searching on the field name.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Register_ChmField_Rfp_4D( Input_Opt,  metadataID, Ptr2Data,     &
                                       State_Chm,  RC,         Ncat         )
!
! !USES:
!
    USE Input_Opt_Mod,      ONLY : OptInput
    USE Registry_Params_Mod
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput),    INTENT(IN)    :: Input_Opt         ! Input Options object
    CHARACTER(LEN=*),  INTENT(IN)    :: metadataID        ! State_Chm field id
    REAL(fp),          POINTER       :: Ptr2Data(:,:,:,:) ! pointer to data
    TYPE(ChmState),    INTENT(IN)    :: State_Chm         ! Obj for chem state
    INTEGER,           OPTIONAL      :: Ncat              ! category index
!
! !INPUT/OUTPUT PARAMETERS:
!
!
! !OUTPUT PARAMETERS:
!
    INTEGER,           INTENT(OUT)   :: RC              ! Success/failure
!
! !REMARKS:
!
! !REVISION HISTORY:
!  20 Sep 2017 - E. Lundgren - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    CHARACTER(LEN=512)     :: ErrMsg
    CHARACTER(LEN=255)     :: ErrMsg_reg, ThisLoc
    CHARACTER(LEN=255)     :: desc, units, perSpecies
    CHARACTER(LEN=255)     :: thisSpcName, thisSpcDesc
    INTEGER                :: N, rank, type,  vloc
    LOGICAL                :: found, onEdges
    TYPE(Species), POINTER :: SpcInfo

    !-----------------------------------------------------------------------
    ! Initialize
    !-----------------------------------------------------------------------
    RC = GC_SUCCESS
    ThisLoc = ' -> at Register_ChmField_Rfp_4D (in Headers/state_chm_mod.F90)'
    ErrMsg_reg = 'Error encountered while registering State_Chm%'

    !-----------------------------------------------------------------------
    ! Initialize
    !-----------------------------------------------------------------------
    CALL Get_Metadata_State_Chm( Input_Opt%amIRoot, metadataID,  Found,  RC, &
                                 desc=desc, units=units, rank=rank,          &
                                 type=type, vloc=vloc, perSpecies=perSpecies )

    ! Trap potential errors
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = TRIM( ErrMsg_reg ) // TRIM( MetadataID ) //                  &
                '; Abnormal exit from routine "Get_Metadata_State_Chm"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    !-----------------------------------------------------------------------
    ! Check that metadata consistent with data size
    !-----------------------------------------------------------------------
    IF ( rank /= 3 ) THEN
       ErrMsg = 'Data dims and metadata rank do not match for ' &
                // TRIM(metadataID)
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    ! Is the data placed on level edges?
    onEdges = ( VLoc == VLocationEdge )

    !-----------------------------------------------------------------------
    ! If tied to all species then register each one
    !-----------------------------------------------------------------------
    IF ( perSpecies == 'ALL' ) THEN

       ! Loop over all species
       DO N = 1, State_Chm%nSpecies

          ! Get name from species database for name and description tags
          SpcInfo     => State_Chm%SpcData(N)%Info
          thisSpcName = TRIM( metadataID ) // '_' // TRIM( SpcInfo%Name )
          thisSpcDesc = TRIM( Desc       ) // ' ' // TRIM( SpcInfo%Name )

          ! Add field to registry
          CALL Registry_AddField( Input_Opt    = Input_Opt,                  &
                                  Registry     = State_Chm%Registry ,        &
                                  State        = State_Chm%State,            &
                                  Variable     = thisSpcName,                &
                                  Description  = thisSpcDesc,                &
                                  Units        = units,                      &
                                  OnLevelEdges = onEdges,                    &
                                  Data3d       = Ptr2Data(:,:,:,N),          &
                                  RC           = RC                         )

          ! Free pointers
          SpcInfo => NULL()

          ! Trap potential errors
          IF ( RC /= GC_SUCCESS ) THEN
             ErrMsg = TRIM( ErrMsg_reg ) // TRIM( MetadataID ) //            &
                      '; Abnormal exit from routine "Registry_AddField"!'
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN
          ENDIF

       ENDDO

    !-----------------------------------------------------------------------
    ! If tied to a given category, only registry that one
    !-----------------------------------------------------------------------
    ELSE IF ( PRESENT(Ncat) ) THEN

       ! Add field to registry
       CALL Registry_AddField( Input_Opt    = Input_Opt,                     &
                               Registry     = State_Chm%Registry ,           &
                               State        = State_Chm%State,               &
                               Variable     = metadataID ,                   &
                               Description  = desc,                          &
                               Units        = units,                         &
                               OnLevelEdges = onEdges,                       &
                               Data3d       = Ptr2Data(:,:,:,Ncat),          &
                               RC           = RC                            )

       ! Trap potential errors
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = TRIM( ErrMsg_reg ) // TRIM( MetadataID ) //               &
                   '; Abnormal exit from Registry_AddField!'
          CALL GC_Error( ErrMsg_reg, RC, ThisLoc )
          RETURN
       ENDIF

    !-----------------------------------------------------------------------
    ! Otherwise, exit with error
    !-----------------------------------------------------------------------
    ELSE
       ErrMsg = 'Handling of PerSpecies metadata ' // TRIM(perSpecies) // &
                ' is not implemented for this combo of data type and size!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

  END SUBROUTINE Register_ChmField_Rfp_4D
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Ind_
!
! !DESCRIPTION: Function IND\_ returns the index of an advected species or
!  chemical species contained in the chemistry state object by name.
!\\
!\\
! !INTERFACE:
!
  FUNCTION Ind_( name, flag ) RESULT( Indx )
!
! !USES:
!
    USE CharPak_Mod, ONLY : To_UpperCase

!
! !INPUT PARAMETERS:
!
    CHARACTER(LEN=*),           INTENT(IN) :: name  ! Species name
    CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: flag  ! Species type
!
! !RETURN VALUE:
!
    INTEGER                                :: Indx  ! Index of this species
!
! !REMARKS:
!   Values of FLAG (case-insensitive):
!   'A' : Returns advected species index
!   'D' : Returns dry-deposition species index
!   'F' : Returns KPP fixed species index
!   'G' : Returns gas-phase species index
!   'H' : Returns hygroscopic-growth species index
!   'K' : Returns KPP master species index
!   'P' : Returns photolysis species index
!   'S' : Returns master species index (aka "ModelId")
!   'V' : Returns KPP variable species index
!   'W' : Returns wet-deposition species index
!
! !REVISION HISTORY:
!  07 Oct 2016 - M. Long     - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER :: N

    !=====================================================================
    ! Ind_ begins here!
    !=====================================================================

    ! Get the ModelId value from the lookup table
    ! NOTE: -1 is used to denote missing species.
    N    = SpcDictLocal%Get( To_UpperCase( TRIM( Name ) ) )
    Indx = N

    ! If N is negative then return -1 to denote the species was not found.
    ! If FLAG is not passed, RETURN the ModelId (regardless of whether the
    ! species was found or missing).
    IF ( ( N < 0 ).or. ( .not. PRESENT( Flag ) ) ) RETURN

    ! For species that were found, return the index specified by FLAG
    SELECT CASE( Flag(1:1) )

       ! Advected species flag
       CASE( 'A', 'a' )
          Indx = SpcDataLocal(N)%Info%AdvectID
          RETURN

       ! Dry-deposited species ID
       CASE( 'D', 'd' )
          Indx = SpcDataLocal(N)%Info%DryDepId
          RETURN

       ! KPP fixed species ID
       CASE( 'F', 'f' )
          Indx = SpcDataLocal(N)%Info%KppFixId
          RETURN

       ! Gas-phase species ID
       CASE( 'G', 'g' )
          Indx = SpcDataLocal(N)%Info%GasSpcId
          RETURN

       ! Hygroscopic growth species ID
       CASE( 'H', 'h' )
          Indx = SpcDataLocal(N)%Info%HygGrthId
          RETURN

       ! KPP chemical species ID
       CASE( 'K', 'k' )
          Indx = SpcDataLocal(N)%Info%KppSpcId
          RETURN

       ! Photolysis species ID
       CASE( 'P', 'p' )
          Indx = SpcDataLocal(N)%Info%PhotolId
          RETURN

       ! Species/ModelID
       CASE ( 'S', 's' )
          Indx = SpcDataLocal(N)%Info%ModelID
          RETURN

       ! KPP variable species ID
       CASE( 'V', 'v' )
          Indx = SpcDataLocal(N)%Info%KppVarId
          RETURN

       ! WetDep ID
       CASE( 'W', 'w' )
          Indx = SpcDataLocal(N)%Info%WetDepId
          RETURN

       CASE DEFAULT
          ! Pass

     END SELECT

  END FUNCTION Ind_
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: GetNumProdLossSpecies
!
! !DESCRIPTION: Saves the number of production and loss diagnostic species
!  in the State\_Chm\%nProdLoss variable.  This will be used to set up the
!  State\_Chm\%Map\_ProdLoss species index vector.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE GetNumProdLossSpecies( Input_Opt, State_Chm, RC )
!
! !USES:
!
    USE GcKpp_Monitor,    ONLY : Fam_Names
    USE GcKpp_Parameters, ONLY : nFam
    USE Input_Opt_Mod,    ONLY : OptInput
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(INOUT) :: Input_Opt   ! Input Options object
    TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Return code
!
! !REMARKS:
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
    ! Scalars
    INTEGER            :: N

    ! Strings
    CHARACTER(LEN=255) :: ErrMsg, ThisLoc

    !=======================================================================
    ! GetProdLossSpecies begins here!
    !=======================================================================

    ! Initialize
    RC      = GC_SUCCESS
    ErrMsg  = ''
    ThisLoc = &
         ' -> at GetNumProdLossSpecies (in module Headers/state_chm_mod.F90)'

    !=======================================================================
    ! Get the number of prod and loss species depending on the simulation
    !=======================================================================
    IF ( Input_Opt%ITS_A_FULLCHEM_SIM ) THEN

       !------------------------------
       ! Full-chemistry simulations
       !------------------------------

       ! Get the # of prod/loss species by querying the first leter of
       ! the species in the Fam_Names array (in gckpp_Monitor.F90)
       DO N = 1, nFam
          IF ( Fam_Names(N)(1:1) == 'L' ) State_Chm%nLoss = State_Chm%nLoss + 1
          IF ( Fam_Names(N)(1:1) == 'P' ) State_Chm%nProd = State_Chm%nProd + 1
       ENDDO

    ELSE IF ( Input_Opt%ITS_A_TAGCO_SIM ) THEN

       !------------------------------
       ! Tagged CO simulation
       !------------------------------

       ! Each advected species can have a loss diagnostic attached ...
       State_Chm%nLoss = State_Chm%nAdvect

       ! ... but no prod diagnostics.  These will get archived by separate
       ! array fields of the State_Diag object (e.g. ProdCOfromISOP, etc.)
       State_Chm%nProd = 0

    ELSE IF ( Input_Opt%ITS_A_TAGO3_SIM ) THEN

       !------------------------------
       ! Tagged O3 simulation
       !-----------------------------

       ! Each advected species can have a prod and loss diagnostic attached
       State_Chm%nLoss = State_Chm%nAdvect
       State_Chm%nProd = State_Chm%nAdvect

    ELSE

       ! Other simulations do not have a prod/loss functionality
       ! but this can be added in if necessary
       State_Chm%nLoss = 0
       State_Chm%nProd = 0

    ENDIF

  END SUBROUTINE GetNumProdLossSpecies
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: MapProdLossSpecies
!
! !DESCRIPTION: Stores the ModelId (from the GEOS-Chem Species Database) of
!  each prod/loss diagnostic species in the State\_Chm\%Map\_ProdLoss vector.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE MapProdLossSpecies( Input_Opt, State_Chm, RC )
!
! !USES:
!
    USE GcKpp_Monitor,    ONLY : Fam_Names
    USE GcKpp_Parameters, ONLY : nFam
    USE Input_Opt_Mod,    ONLY : OptInput
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(INOUT) :: Input_Opt   ! Input Options object
    TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Return code
!
! !REMARKS:
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
    ! Scalars
    INTEGER            :: Id,     N
    INTEGER            :: P,      L

    ! Strings
    CHARACTER(LEN=36)  :: Name
    CHARACTER(LEN=255) :: ErrMsg, ThisLoc

    !=======================================================================
    ! GetProdLossSpecies begins here!
    !=======================================================================

    ! Initialize
    RC      = GC_SUCCESS
    P       = 0
    L       = 0
    ErrMsg  = ''
    ThisLoc = &
         ' -> at MapProdLossSpecies (in module Headers/state_chm_mod.F90)'

    !=======================================================================
    ! Get the number of prod and loss species depending on the simulation
    !=======================================================================
    IF ( Input_Opt%ITS_A_FULLCHEM_SIM ) THEN

       !--------------------------------------------------------------------
       ! Full-chemistry simulations
       !--------------------------------------------------------------------

       ! Loop over the number of prod/loss species
       DO N = 1, nFam

          ! Get the KPP prod/loss species from the FAM_NAMES
          ! array in the gckpp_Parameters.F90 module.
          ! NOTE: This is the KPP ID number (index of "VAR" array)
          ! and not the GEOS-Chem "master" species index!!!
          Id = Ind_( TRIM( Fam_Names(N) ), 'K' )

          ! Add the species
          IF ( Id > 0 ) THEN

             ! KPP prod/loss species name
             Name = TRIM( Fam_Names(N) )

             ! Fix the name so that it is of the form Prod_<spcname> or
             ! Loss_<spcname>.  This will facilitate the new diagnostics.
             IF ( Name(1:1) == 'L' ) THEN
                L                      = L + 1
                State_Chm%Map_Loss(L)  = Id
                State_Chm%Name_Loss(L) = 'Loss_' // TRIM( Name(2:) )
             ELSE IF ( Name(1:1) == 'P' ) THEN
                P                      = P + 1
                State_Chm%Map_Prod(P)  = Id
                State_Chm%Name_Prod(P) = 'Prod_' // TRIM( Name(2:) )
             ELSE
                ErrMsg = 'Invalid prod/loss species name!' //                &
                          TRIM( Fam_Names(N) )
                CALL GC_Error( ErrMsg, RC, ThisLoc )
                RETURN
             ENDIF

          ELSE

             ! Invalid species, exit with error!
             ErrMsg = 'Could not locate KPP prod/loss species: ' //          &
                      TRIM( Fam_Names(N) )                       // '!'
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN

          ENDIF

       ENDDO

    ELSE IF ( Input_Opt%ITS_A_TAGCO_SIM ) THEN

       !--------------------------------------------------------------------
       ! Tagged CO simulations
       !--------------------------------------------------------------------

       ! Each advected species can have an attached loss diagnostic ...
       DO Id = 1, State_Chm%nLoss
          Name = 'Loss_' // TRIM( State_Chm%SpcData(Id)%Info%Name )
          State_Chm%Name_Loss(Id) = TRIM( Name )
          State_Chm%Map_Loss(Id)  = Id
       ENDDO

       ! ... but not an attached prod diagnostic.  These will be
       ! archived by other fields of the State_Diag object.
       State_Chm%Name_Prod => NULL()
       State_Chm%Map_Prod  => NULL()

    ELSE IF ( Input_Opt%ITS_A_TAGO3_SIM ) THEN

       !--------------------------------------------------------------------
       ! Tagged O3 simulations
       !--------------------------------------------------------------------

       ! Each advected species can have an attached loss diagnostic ...
       DO Id = 1, State_Chm%nLoss
          Name = 'Loss_' // TRIM( State_Chm%SpcData(Id)%Info%Name )
          State_Chm%Name_Loss(Id) = TRIM( Name )
          State_Chm%Map_Loss(Id)  = Id
       ENDDO

       ! ... as well as an attached production diagnostic
       DO Id = 1, State_Chm%nProd
          Name = 'Prod_' // TRIM( State_Chm%SpcData(Id)%Info%Name )
          State_Chm%Name_Prod(Id) = TRIM( Name )
          State_Chm%Map_Prod(Id)  = Id
       ENDDO

    ELSE

       !--------------------------------------------------------------------
       ! Other simulations do not have prod/loss capability
       !--------------------------------------------------------------------
       State_Chm%Name_Prod => NULL()
       State_Chm%Map_Prod  => NULL()
    ENDIF

  END SUBROUTINE MapProdLossSpecies
!EOC
END MODULE State_Chm_Mod
