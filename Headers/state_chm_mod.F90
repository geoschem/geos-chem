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
  USE Species_Mod                        ! For species database and conc objects

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
  PRIVATE :: Init_and_Register
  PRIVATE :: Init_Mapping_Vectors
  PRIVATE :: Register_ChmField
  PRIVATE :: Zero_State_Chm
!
! !PRIVATE DATA MEMBERS:
!
  TYPE(SpcPtr), PRIVATE, POINTER :: SpcDataLocal(:)  ! Local pointer to
                                                     ! StateChm%SpcData for
                                                     ! availability to IND_
  TYPE(SpcIndCt)                 :: SpcCount

  TYPE(dictionary_t), PRIVATE    :: SpcDictLocal     ! Private copy of the
                                                     ! Fortran Hash table for
                                                     ! availability to IND_

  INTEGER, PRIVATE               :: nChmState = 0    ! # chemistry states,

  !==========================================================================
  ! Derived type for Chemistry State
  !==========================================================================
  TYPE, PUBLIC :: ChmState

     !-----------------------------------------------------------------------
     ! Count of each type of species
     !-----------------------------------------------------------------------
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
     INTEGER                    :: nOmitted             ! # of omitted species
     INTEGER                    :: nPhotol              ! # photolysis species
     INTEGER                    :: nProd                ! # of prod species
     INTEGER                    :: nRadNucl             ! # of radionuclides
     INTEGER                    :: nWetDep              ! # wetdep species

     !-----------------------------------------------------------------------
     ! Mapping vectors to subset types of species
     !-----------------------------------------------------------------------
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
     INTEGER,           POINTER :: Map_RadNucl(:      ) ! Radionuclide IDs
     INTEGER,           POINTER :: Map_WetDep (:      ) ! Wetdep species IDs
     INTEGER,           POINTER :: Map_WL     (:      ) ! Wavelength bins in fjx

     !-----------------------------------------------------------------------
     ! Physical properties & indices for each species
     !-----------------------------------------------------------------------
     TYPE(SpcPtr),      POINTER :: SpcData    (:      ) ! GC Species database
     TYPE(dictionary_t)         :: SpcDict              ! Species dictionary

     !-----------------------------------------------------------------------
     ! Chemical species
     !-----------------------------------------------------------------------
     TYPE(SpcConc),     POINTER :: Species (:      )    ! Vector for species
                                                        ! concentrations
                                                        !  [kg/kg dry air]
     CHARACTER(LEN=20)          :: Spc_Units            ! Species units

#ifdef ADJOINT
     REAL(fp),          POINTER :: SpeciesAdj (:,:,:,:) ! Species adjoint variables
     REAL(fp),          POINTER :: CostFuncMask(:,:,:)  ! cost function volume mask
#endif

     !----------------------------------------------------------------------
     ! Boundary conditions
     !-----------------------------------------------------------------------
     REAL(fp),          POINTER :: BoundaryCond(:,:,:,:)! Boundary conditions
                                                        !  [kg/kg dry air]

     !-----------------------------------------------------------------------
     ! RRTMG state variables
     !-----------------------------------------------------------------------
     INTEGER                    :: RRTMG_iSeed
     INTEGER                    :: RRTMG_iCld

     !-----------------------------------------------------------------------
     ! Aerosol quantities
     !-----------------------------------------------------------------------
     REAL(fp),          POINTER :: AeroArea   (:,:,:,:) ! Aerosol Area [cm2/cm3]
     REAL(fp),          POINTER :: AeroRadi   (:,:,:,:) ! Aerosol Radius [cm]
     REAL(fp),          POINTER :: WetAeroArea(:,:,:,:) ! Aerosol Area [cm2/cm3]
     REAL(fp),          POINTER :: WetAeroRadi(:,:,:,:) ! Aerosol Radius [cm]
     REAL(fp),          POINTER :: AeroH2O    (:,:,:,:) ! Aerosol water [cm3/cm3]
     REAL(fp),          POINTER :: GammaN2O5  (:,:,:,:) ! N2O5 aerosol uptake [unitless]
     REAL(fp),          POINTER :: SSAlk      (:,:,:,:) ! Sea-salt alkalinity[-]
     REAL(fp),          POINTER :: H2O2AfterChem(:,:,:) ! H2O2, SO2 [v/v]
     REAL(fp),          POINTER :: SO2AfterChem (:,:,:) !  after sulfate chem
     REAL(fp),          POINTER :: OMOC           (:,:) ! OM:OC Ratio [unitless]
     REAL(fp),          POINTER :: OMOC_POA       (:,:) ! OM:OC Ratio (OCFPOA) [unitless]
     REAL(fp),          POINTER :: OMOC_OPOA      (:,:) ! OM:OC Ratio (OCFOPOA) [unitless]
     REAL(fp),          POINTER :: ACLArea      (:,:,:) ! Fine Cl- Area [cm2/cm3]
     REAL(fp),          POINTER :: ACLRadi      (:,:,:) ! Fine Cl- Radius [cm]
     REAL(fp),          POINTER :: QLxpHCloud   (:,:,:) !
     REAL(fp),          POINTER :: SoilDust   (:,:,:,:) ! Soil dust [kg/m3]
     REAL(fp),          POINTER :: ORVCsesq     (:,:,:) ! Sesquiterpenes mass [kg/box]

     !-----------------------------------------------------------------------
     ! Fields for nitrogen deposition
     !-----------------------------------------------------------------------
     REAL(fp),          POINTER :: DryDepNitrogen (:,:) ! Dry deposited N
     REAL(fp),          POINTER :: WetDepNitrogen (:,:) ! Wet deposited N

     !-----------------------------------------------------------------------
     ! Cloud quantities
     !-----------------------------------------------------------------------
     REAL(fp),          POINTER :: pHCloud    (:,:,:  ) ! Cloud pH [-]
     REAL(fp),          POINTER :: isCloud    (:,:,:  ) ! Cloud presence [-]

     !-----------------------------------------------------------------------
     ! Fields for KPP solver
     !-----------------------------------------------------------------------
     REAL(fp),          POINTER :: KPPHvalue  (:,:,:  ) ! H-value for Rosenbrock
                                                        !  solver
     !-----------------------------------------------------------------------
     ! Fields for UCX mechanism
     !-----------------------------------------------------------------------
     REAL(f4),          POINTER :: STATE_PSC  (:,:,:  ) ! PSC type (see Kirner
                                                        !  et al. 2011, GMD)
     REAL(fp),          POINTER :: KHETI_SLA  (:,:,:,:) ! Strat. liquid aerosol
                                                        !  reaction cofactors

     !-----------------------------------------------------------------------
     ! For isoprene SOA via ISORROPIA
     !-----------------------------------------------------------------------
     REAL(fp),          POINTER :: IsorropAeropH  (:,:,:,:) ! ISORROPIA aero pH
     REAL(fp),          POINTER :: IsorropHplus   (:,:,:,:) ! H+ conc [M]
     REAL(fp),          POINTER :: IsorropAeroH2O (:,:,:,:) ! ISORROPIA aero H2O
     REAL(fp),          POINTER :: IsorropSulfate (:,:,:  ) ! Sulfate conc [M]
     REAL(fp),          POINTER :: IsorropNitrate (:,:,:,:) ! Nitrate conc [M]
     REAL(fp),          POINTER :: IsorropChloride(:,:,:,:) ! Chloride conc [M]
     REAL(fp),          POINTER :: IsorropBisulfate(:,:,:  )! Bisulfate conc [M]

     !-----------------------------------------------------------------------
     ! For the tagged Hg simulation
     !-----------------------------------------------------------------------
     REAL(fp),          POINTER :: OceanHg0(:,:)        ! Hg(0)  ocean mass [kg]
     REAL(fp),          POINTER :: OceanHg2(:,:)        ! Hg(II) ocean mass [kg]
     REAL(fp),          POINTER :: OceanHgP(:,:)        ! HgP    ocean mass [kg]
     REAL(fp),          POINTER :: SnowHgOcean(:,:)     ! Reducible Hg snowpack
                                                        !  on ocean [kg]
     REAL(fp),          POINTER :: SnowHgLand(:,:)      ! Reducible Hg snowpack
                                                        !  on land [kg]
     REAL(fp),          POINTER :: SnowHgOceanStored(:,:)   ! Non-reducible Hg
                                                            !  snowpack on ocean
     REAL(fp),          POINTER :: SnowHgLandStored(:,:)    ! Non-reducible Hg
                                                            !  snowpack on land

     !----------------------------------------------------------------------
     ! For HOBr + S(IV) heterogeneous chemistry
     !----------------------------------------------------------------------
     REAL(fp),          POINTER :: HSO3_AQ    (:,:,:  ) ! Cloud bisulfite/SO2 ratio
     REAL(fp),          POINTER :: SO3_AQ     (:,:,:  ) ! Cloud sulfite/SO2 ratio
     REAL(fp),          POINTER :: fupdateHOBr(:,:,:  ) ! Correction factor for
                                                        ! HOBr removal by SO2
                                                        ! [unitless]
     REAL(fp),          POINTER :: fupdateHOCl(:,:,:  ) ! Correction factor for
                                                        ! HOCl removal by SO2
                                                        ! [unitless]
     !-----------------------------------------------------------------------
     ! Fields for dry deposition
     !-----------------------------------------------------------------------
     REAL(fp),          POINTER :: Iodide       (:,:  ) ! Ocn surf iodide [nM]
     REAL(fp),          POINTER :: Salinity     (:,:  ) ! Ocn surf salinity [PSU]
     REAL(fp),          POINTER :: DryDepFreq (:,:,:  ) ! Drydep freq [s-1]
     REAL(f8),          POINTER :: DryDepVel  (:,:,:  ) ! Dry deposition velocities
                                                        ! [m/s] - use REAL8 in drydep
#if defined( MODEL_GEOS )
     REAL(fp),          POINTER :: DryDepRa2m (:,:    ) ! 2m  aerodynamic resistance
     REAL(fp),          POINTER :: DryDepRa10m(:,:    ) ! 10m aerodynamic resistance
#endif
     REAL(fp),          POINTER :: JOH        (:,:    ) ! OH J-value
     REAL(fp),          POINTER :: JNO2       (:,:    ) ! NO2 J-value

     !-----------------------------------------------------------------------
     ! Fields for non-local PBL mixing
     !-----------------------------------------------------------------------
     REAL(fp),          POINTER :: SurfaceFlux(:,:,:  )

     !-----------------------------------------------------------------------
     ! Fields for Linoz stratospheric ozone algorithm
     !----------------------------------------------------------------------
     REAL(fp),          POINTER :: TLSTT      (:,:,:,:) ! TLSTT (I,J,L,LINOZ_NFIELDS)

     !------------------------------------------------------------------------
     ! Fields for Gan Luo et al Wetdep scheme (GMD-12-3439-2019)
     !----------------------------------------------------------------------
     REAL(fp),          POINTER :: KRATE      (:,:,:  )
     REAL(fp),          POINTER :: QQ3D       (:,:,:  )
     REAL(fp),          POINTER :: pHRain     (:,:,:  ) ! Rain pH [-]
     REAL(fp),          POINTER :: QQpHRain   (:,:,:  ) ! Rain pH*QQ3D [-]
     REAL(fp),          POINTER :: QQRain     (:,:,:  ) ! Rain QQ3D [-]

     !-----------------------------------------------------------------------
     ! Fields for setting mean surface CH4 from HEMCO
     !-----------------------------------------------------------------------
     REAL(fp),          POINTER :: SFC_CH4    (:,:    )

     !-----------------------------------------------------------------------
     ! Fields for TOMS overhead ozone column data
     !-----------------------------------------------------------------------
     REAL(fp),          POINTER :: TO3_DAILY  (:,:    ) ! Daily overhead ozone
     REAL(fp),          POINTER :: TOMS1      (:,:    )
     REAL(fp),          POINTER :: TOMS2      (:,:    )

     !-----------------------------------------------------------------------
     ! Fields for UCX (moved from module)
     ! Many of these fields are not sized NX x NY, and thus cannot be allocated
     ! by Init_and_Register. They will be handled by Init_UCX as appropriate,
     ! and only stored in State_Chm so they can be separated per chemistry state.
     ! (hplin, 1/5/23)
     !-----------------------------------------------------------------------
     REAL(fp),          POINTER :: UCX_REGRID (:,:    )
     REAL(fp),          POINTER :: UCX_PLEVS  (:      ) ! Pressure levels of 2D data (hPa)
     REAL(fp),          POINTER :: UCX_LATS   (:      ) ! Latitude edges of 2D data (deg)
     REAL(fp),          POINTER :: RAD_AER    (:,:,:,:) ! Strat. aerosol radius (cm)
     REAL(fp),          POINTER :: KG_AER     (:,:,:,:) ! Aerosol mass (kg/box)
     REAL(fp),          POINTER :: SAD_AER    (:,:,:,:) ! Aerosol surface area density (cm2/cm3)
     REAL(fp),          POINTER :: NDENS_AER  (:,:,:,:) ! Aerosol number density (#/m3)
     REAL(fp),          POINTER :: RHO_AER    (:,:,:,:) ! Aerosol mass density (kg/m3 aerosol)
     REAL(fp),          POINTER :: AERFRAC    (:,:,:,:) ! Mass fraction of species in liquid aerosols
     INTEGER,           POINTER :: AERFRACIND (:      ) ! Indices of liquid aerosol species
     REAL(fp),          POINTER :: NOX_O      (:,:,:,:) ! Monthly mean noontime O3P/O1D for NOx calcs
     REAL(fp),          POINTER :: NOX_J      (:,:,:,:) ! Monthly mean noontime J-rates for NOx calcs
     REAL(fp),          POINTER :: SO4_TOPPHOT(:,:    ) ! Photolysis rate at the top of the chemgrid (1/s)

     !=================================================================
     ! Variables to use NOx coefficients in ESMF / grid-independent
     ! envionment. The NOx coefficients are climatological 2D
     ! (lat/lev/12 months) data that are currently available for
     ! horizontal (latitude) resolutions of 2 and 4 degrees. For other
     ! resolutions, the horizontal data becomes mapped onto the
     ! simulation grid (see GET_JJNOX).
     ! Similar to the surface mixing ratio boundary conditions, we now
     ! read all the NOx coefficients during initialization to avoid
     ! additional I/O calls during run time (ckeller, 05/12/2014).
     !=================================================================
     REAL(fp),          POINTER :: NOXCOEFF   (:,:,:,:)
     REAL(fp),          POINTER :: NOXLAT     (:      )
     INTEGER                    :: JJNOXCOEFF

     !-----------------------------------------------------------------------
     ! Switches to enable SO2 cloud chemistry and seasalt chemistry in
     ! sulfate_mod (TRUE) or in the KPP mechanism (FALSE).
     !-----------------------------------------------------------------------
     LOGICAL                    :: Do_SulfateMod_Cld
     LOGICAL                    :: Do_SulfateMod_SeaSalt

#if defined(MODEL_CESM)
     !-----------------------------------------------------------------------
     ! Fields for CESM interface to GEOS-Chem
     !-----------------------------------------------------------------------
     REAL(fp),          POINTER :: H2SO4_PRDR (:,:,:  ) ! H2SO4 prod rate [mol/mol]
#endif

     !-----------------------------------------------------------------------
     ! Fields for CH4 specialty simulation
     !-----------------------------------------------------------------------
     REAL(fp),          POINTER :: BOH        (:,:,:  ) ! OH values [molec/cm3]
     REAL(fp),          POINTER :: BCl        (:,:,:  ) ! Cl values [v/v]
     REAL(fp),          POINTER :: CH4_EMIS   (:,:,:  ) ! CH4 emissions [kg/m2/s].
                                                        ! third dim is cat, total 15

     !-----------------------------------------------------------------------
     ! Registry of variables contained within State_Chm
     !-----------------------------------------------------------------------
     CHARACTER(LEN=4)           :: State     = 'CHEM'   ! Name of this state
     TYPE(MetaRegItem), POINTER :: Registry  => NULL()  ! Registry object
     TYPE(dictionary_t)         :: RegDict              ! Registry lookup table

     !-----------------------------------------------------------------------
     ! Carbon stuff for GEOS 
     !-----------------------------------------------------------------------
#if defined( MODEL_GEOS )
     INTEGER            :: CO2fromGOCART
     CHARACTER(LEN=255) :: impCO2name
     INTEGER            :: numphoto
     INTEGER            :: nxdo
     INTEGER            :: nlam
     INTEGER            :: nsza
     INTEGER            :: numo3
     INTEGER            :: nts
     INTEGER            :: aqsize

     REAL, POINTER      :: sdat(:,:,:,:)
     REAL, POINTER      :: o2jdat(:,:,:)
     REAL, POINTER      :: sza_tab(:)
     REAL, POINTER      :: o3_tab(:,:)
     REAL, POINTER      :: xtab(:,:,:)
     REAL, POINTER      :: CH2O_aq(:)
     REAL, POINTER      :: rlam(:)
#endif

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
  INTERFACE Init_and_Register
     MODULE PROCEDURE Init_and_Register_R4_2D
     MODULE PROCEDURE Init_and_Register_R4_3D
     MODULE PROCEDURE Init_and_Register_R4_4D
     MODULE PROCEDURE Init_and_Register_R8_2D
     MODULE PROCEDURE Init_and_Register_R8_3D
     MODULE PROCEDURE Init_and_Register_R8_4D
  END INTERFACE Init_and_Register

  INTERFACE Register_ChmField
     MODULE PROCEDURE Register_ChmField_R4_2D
     MODULE PROCEDURE Register_ChmField_R4_3D
     MODULE PROCEDURE Register_ChmField_R4_4D
     MODULE PROCEDURE Register_ChmField_R8_2D
     MODULE PROCEDURE Register_ChmField_R8_3D
     MODULE PROCEDURE Register_ChmField_R8_4D
  END INTERFACE Register_ChmField

CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Zero_State_Chm
!
! !DESCRIPTION: Nullifies and/or zeroes all fields of State\_Chm.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Zero_State_Chm( State_Chm, RC )
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ChmState), INTENT(INOUT) :: State_Chm
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC
!
! !REVISION HISTORY:
!  21 Sep 2020 - R. Yantosca - Initial version
!  See the subsequent Git history with the gitk browser!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    !========================================================================
    ! Assume success
    !========================================================================
    RC = GC_SUCCESS

    !========================================================================
    ! Initialize or nullify each member of State_Chm
    ! This will prevent potential deallocation errors
    !========================================================================

    ! Counters
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
    State_Chm%nOmitted          =  0
    State_Chm%nSpecies          =  0
    State_Chm%nWetDep           =  0

    ! Mapping vectors
    State_Chm%Map_Advect        => NULL()
    State_Chm%Map_Aero          => NULL()
    State_Chm%Map_DryAlt        => NULL()
    State_Chm%Map_DryDep        => NULL()
    State_Chm%Map_GasSpc        => NULL()
    State_Chm%Map_HygGrth       => NULL()
    State_Chm%Map_KppVar        => NULL()
    State_Chm%Map_KppFix        => NULL()
    State_Chm%Map_KppSpc        => NULL()
    State_Chm%Map_Loss          => NULL()
    State_Chm%Name_Loss         => NULL()
    State_Chm%Map_Photol        => NULL()
    State_Chm%Map_Prod          => NULL()
    State_Chm%Name_Prod         => NULL()
    State_Chm%Map_RadNucl       => NULL()
    State_Chm%Map_WetDep        => NULL()
    State_Chm%Map_WL            => NULL()

    ! Species-based quantities
    State_Chm%SpcData           => NULL()
    State_Chm%Species           => NULL()
    State_Chm%Spc_Units         =  ''
    State_Chm%BoundaryCond      => NULL()

#ifdef ADJOINT
    ! Chemical species adjoint variables
    State_Chm%SpeciesAdj    => NULL()
    State_Chm%CostFuncMask  => NULL()
#endif

    ! RRTMG state
    State_Chm%RRTMG_iSeed       = 0
    State_Chm%RRTMG_iCld        = 0

    ! Aerosol and chemistry quantities
    State_Chm%AeroArea          => NULL()
    State_Chm%AeroRadi          => NULL()
    State_Chm%WetAeroArea       => NULL()
    State_Chm%WetAeroRadi       => NULL()
    State_Chm%AeroH2O           => NULL()
    State_Chm%GammaN2O5         => NULL()
    State_Chm%SSAlk             => NULL()
    State_Chm%OMOC              => NULL()
    State_Chm%OMOC_POA          => NULL()
    State_Chm%OMOC_OPOA         => NULL()
    State_Chm%DryDepNitrogen    => NULL()
    State_Chm%WetDepNitrogen    => NULL()
    State_Chm%pHCloud           => NULL()
    State_Chm%isCloud           => NULL()
    State_Chm%QLxpHCloud        => NULL()
    State_Chm%ORVCsesq          => NULL()
    State_Chm%KPPHvalue         => NULL()
    State_Chm%STATE_PSC         => NULL()
    State_Chm%KHETI_SLA         => NULL()
    State_Chm%ACLArea           => NULL()
    State_Chm%ACLRadi           => NULL()
    State_Chm%SoilDust          => NULL()
    State_Chm%HSO3_AQ           => NULL()
    State_Chm%SO3_AQ            => NULL()
    State_Chm%fupdateHOBr       => NULL()
    State_Chm%fupdateHOCl       => NULL()
    State_Chm%TLSTT             => NULL()
    State_Chm%TO3_DAILY         => NULL()
    State_Chm%TOMS1             => NULL()
    State_Chm%TOMS2             => NULL()
    State_Chm%BOH               => NULL()
    State_Chm%BCl               => NULL()
    State_Chm%CH4_EMIS          => NULL()
    State_Chm%SFC_CH4           => NULL()

    State_Chm%UCX_REGRID        => NULL()
    State_Chm%UCX_PLEVS         => NULL()
    State_Chm%UCX_LATS          => NULL()
    State_Chm%RAD_AER           => NULL()
    State_Chm%KG_AER            => NULL()
    State_Chm%SAD_AER           => NULL()
    State_Chm%NDENS_AER         => NULL()
    State_Chm%RHO_AER           => NULL()
    State_Chm%AERFRAC           => NULL()
    State_Chm%AERFRACIND        => NULL()
    State_Chm%NOX_O             => NULL()
    State_Chm%NOX_J             => NULL()
    State_Chm%NOXCOEFF          => NULL()
    State_Chm%NOXLAT            => NULL()

    ! Emissions and drydep quantities
    State_Chm%Iodide            => NULL()
    State_Chm%Salinity          => NULL()
    State_Chm%DryDepFreq        => NULL()
    State_Chm%DryDepVel         => NULL()
#ifdef MODEL_GEOS
    State_Chm%DryDepRa2m        => NULL()
    State_Chm%DryDepRa10m       => NULL()
#endif
    State_Chm%JOH               => NULL()
    State_Chm%JNO2              => NULL()

    ! Non-local PBL mixing quantities
    State_Chm%SurfaceFlux       => NULL()

    ! Wetdep quantities
    State_Chm%H2O2AfterChem     => NULL()
    State_Chm%SO2AfterChem      => NULL()
    State_Chm%QQ3D              => NULL()
    State_Chm%KRATE             => NULL()
    State_Chm%QQ3D              => NULL()
    State_Chm%pHRain            => NULL()
    State_Chm%QQpHRain          => NULL()
    State_Chm%QQRain            => NULL()

    ! Isoprene SOA
    State_Chm%IsorropAeropH     => NULL()
    State_Chm%IsorropHplus      => NULL()
    State_Chm%IsorropAeroH2O    => NULL()
    State_Chm%IsorropSulfate    => NULL()
    State_Chm%IsorropNitrate    => NULL()
    State_Chm%IsorropChloride   => NULL()
    State_Chm%IsorropBisulfate  => NULL()

    ! Hg simulation quantities
    State_Chm%OceanHg0          => NULL()
    State_Chm%OceanHg2          => NULL()
    State_Chm%OceanHgP          => NULL()
    State_Chm%SnowHgOcean       => NULL()
    State_Chm%SnowHgLand        => NULL()
    State_Chm%SnowHgOceanStored => NULL()
    State_Chm%SnowHgLandStored  => NULL()

    ! Flags to toggle sulfate-mod computations or KPP computations
    ! TRUE  = use sulfate_mod
    ! FALSE = use KPP computations
    State_Chm%Do_SulfateMod_Cld     = .FALSE.
    State_Chm%Do_SulfateMod_SeaSalt = .FALSE.

#ifdef MODEL_GEOS
    ! Add quantities for coupling to the NASA/GEOS ESM
    State_Chm%CO2fromGOCART     = .FALSE.
    State_Chm%impCO2name        = "unknown" 
    State_Chm%numphoto          = 0
    State_Chm%nxdo              = 0
    State_Chm%nlam              = 0
    State_Chm%nsza              = 0
    State_Chm%numo3             = 0
    State_Chm%nts               = 0
    State_Chm%aqsize            = 0
    State_Chm%sdat              => NULL()
    State_Chm%o2jdat            => NULL()
    State_Chm%sza_tab           => NULL()
    State_Chm%o3_tab            => NULL()
    State_Chm%xtab              => NULL()
    State_Chm%CH2O_aq           => NULL()
    State_Chm%rlam              => NULL()
#endif
#ifdef MODEL_CESM
    ! Add quantities for coupling to CESM
    State_Chm%H2SO4_PRDR        => NULL()
#endif

  END SUBROUTINE Zero_State_Chm
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
    USE ErrCode_Mod
    USE CharPak_Mod,          ONLY : To_UpperCase
    USE CMN_Size_Mod,         ONLY : NDUST, NAER
    USE GCKPP_Parameters,     ONLY : NSPEC
    USE Input_Opt_Mod,        ONLY : OptInput
    USE Species_Database_Mod, ONLY : Init_Species_Database
    USE State_Grid_Mod,       ONLY : GrdState
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
    INTEGER                 :: C,          N
    INTEGER                 :: nKHLSA,     nAerosol,   nMatches

    ! Strings
    CHARACTER(LEN=255)      :: errMsg_ir,  errMsg
    CHARACTER(LEN=255)      :: chmId,      thisLoc

    ! String arrays
    CHARACTER(LEN=31)       :: fieldId(14)

    ! Objects
    TYPE(Species),  POINTER :: ThisSpc
    INTEGER,        POINTER :: CheckIds(:)
    REAL(fp),       POINTER :: Ptr2data(:,:,:)

    !========================================================================
    ! Init_State_Chm begins here!
    !========================================================================

    ! Initialize
    RC         =  GC_SUCCESS
    nAerosol   =  NDUST + NAER
    Ptr2data   => NULL()
    ThisSpc    => NULL()
    errMsg     =  ''
    errMsg_ir  =  'Error encountered in "Init_and_Register", chmId = '
    thisLoc    =  &
         ' -> at Init_State_Chm (in module Headers/state_chm_mod.F90)'

    ! Nullify or zero all State_Chm variables
    CALL Zero_State_Chm( State_Chm, RC )

    ! Count the # of chemistry states we have initialized, so SpcData(Local)
    ! is not deallocated until the last ChmState is cleaned up.
    ! This avoids dangling pointers with detrimental effects. (hplin, 8/3/18)
    nChmState = nChmState + 1

    ! Nullify SpcDataLocal upon the first allocation of the Chemistry
    ! State object.  This will avoid undefined pointer errors.
    IF ( nChmState == 1 ) THEN
       SpcDataLocal => NULL()
    ENDIF

    !========================================================================
    ! Do sulfur sea-salt and in-cloud chemistry as part of the KPP-generated
    ! chemical mechanism for all full-chemistry simulations.  For aerosol-
    ! only simulations, do the sulfur chemistry rxns in sulfate_mod.
    !========================================================================
!    IF ( Input_Opt%ITS_A_FULLCHEM_SIM ) THEN
!       State_Chm%Do_SulfateMod_Seasalt = .FALSE.
!       State_Chm%Do_SulfateMod_Cld     = .FALSE.
!    ENDIF

    !========================================================================
    ! Populate the species database object field
    ! (assumes Input_Opt has already been initialized)
    !========================================================================
    IF ( ASSOCIATED( SpcDataLocal ) ) THEN

       ! If the species database has already been initialized on this core,
       ! State_Chm%SpcDataLocal in already contains a copy of the species
       ! metadata.  It can be directly associated to this new chemistry state.
       ! (assumes one core will run one copy of G-C with the same species DB)
       State_Chm%SpcData => SpcDataLocal

    ELSE

       ! Otherwise, initialize the species database by reading the YAML file.
       CALL Init_Species_Database( Input_Opt = Input_Opt,                    &
                                   SpcData   = State_Chm%SpcData,            &
                                   SpcCount  = SpcCount,                     &
                                   RC        = RC                           )

       ! Trap potential errors
       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = 'Error encountered in routine "Init_Species_Database"!'
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       ! Point to a private module copy of the species database
       ! which will be used by the Ind_ indexing function
       SpcDataLocal => State_Chm%SpcData

    ENDIF

    !========================================================================
    ! Before proceeding, make sure none of the species has a blank name,
    ! because this has the potential to halt the run inadvertently.
    !========================================================================

    ! Total number of "real" species (excluding "dummy" placeholder species)
    State_Chm%nSpecies =  SpcCount%nRealSpc

    ! Exit if any species name is blank
    DO N = 1, State_Chm%nSpecies
       IF ( LEN_TRIM(  State_Chm%SpcData(N)%Info%Name ) == 0 ) THEN
          WRITE( ErrMsg, '("Species number ", i6, " has a blank name!")' ) N
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
    ENDDO

    !========================================================================
    ! Determine the number of advected, drydep, wetdep, and total species
    !========================================================================

    ! Get the number of advected, dry-deposited, KPP chemical species,
    ! and and wet-deposited species.  Also return the # of Hg0, Hg2, and
    ! HgP species (these are zero unless the Hg simulation is used).
    State_Chm%nAdvect  = SpcCount%nAdvect
    State_Chm%nAeroSpc = SpcCount%nAeroSpc
    State_Chm%nDryAlt  = SpcCount%nDryAlt
    State_Chm%nDryDep  = SpcCount%nDryDep
    State_Chm%nGasSpc  = SpcCount%nGasSpc
    State_Chm%nHygGrth = SpcCount%nHygGrth
    State_Chm%nKppVar  = SpcCount%nKppVar
    State_Chm%nKppFix  = SpcCount%nKppFix
    State_Chm%nKppSpc  = SpcCount%nKppSpc
    State_Chm%nOmitted = SpcCount%nOmitted
    State_Chm%nPhotol  = SpcCount%nPhotol
    State_Chm%nRadNucl = SpcCount%nRadNucl
    State_Chm%nWetDep  = SpcCount%nWetDep

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

    !========================================================================
    ! Populate the species lookup table, for quick index lookup via Ind_
    !========================================================================

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
    !CALL State_Chm%SpcDict%Show()

    !========================================================================
    ! Exit if this is a dry-run simulation
    !========================================================================
    IF ( Input_Opt%DryRun ) THEN
       RC = GC_SUCCESS
       RETURN
    ENDIF

    !========================================================================
    ! Initialize the 1-D mapping vectors (e.g. State_Chm%Map_DryDep)
    !========================================================================
    CALL Init_Mapping_Vectors( Input_Opt, State_Chm, RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error encountered in "Init_Mapping_Vectors" routine!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF


    !========================================================================
    ! Allocate and initialize chemical species fields
    !========================================================================
    ALLOCATE( State_Chm%Species( State_Chm%nSpecies ), STAT=RC )
    DO N = 1, State_Chm%nSpecies
#if defined ( MODEL_GCHPCTM )
       ! Species concentration array pointers will be set to point to MAPL internal state
       ! every timestep when intstate level values are flipped to match GEOS-Chem standard
       State_Chm%Species(N)%Conc => NULL()
#else
       ALLOCATE( State_Chm%Species(N)%Conc( State_Grid%NX, &
                                            State_Grid%NY, &
                                            State_Grid%NZ ), STAT=RC )
       State_Chm%Species(N)%Conc = 0.0_f8
#endif
    ENDDO

#ifdef ADJOINT
    !========================================================================
    ! Allocate and initialize chemical species fields
    !========================================================================
    chmID = 'SpeciesAdj'
    CALL Init_and_Register(                                                  &
         Input_Opt  = Input_Opt,                                             &
         State_Chm  = State_Chm,                                             &
         State_Grid = State_Grid,                                            &
         chmId      = chmId,                                                 &
         Ptr2Data   = State_Chm%SpeciesAdj,                                  &
         nSlots     = State_Chm%nSpecies,                                    &
         RC         = RC                                                    )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( chmId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !========================================================================
    ! Allocate and initialize chemical species fields
    !========================================================================
    chmID = 'CostFuncMask'
    CALL Init_and_Register(                                                  &
         Input_Opt  = Input_Opt,                                             &
         State_Chm  = State_Chm,                                             &
         State_Grid = State_Grid,                                            &
         chmId      = chmId,                                                 &
         Ptr2Data   = State_Chm%CostFuncMask,                                &
         RC         = RC                                                    )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( chmId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
#endif

    !========================================================================
    ! Boundary conditions (only needed for nested grid simulations)
    !========================================================================
    IF ( State_Grid%NestedGrid ) THEN
       chmID = 'BoundaryCond'
       CALL Init_and_Register(                                               &
            Input_Opt  = Input_Opt,                                          &
            State_Chm  = State_Chm,                                          &
            State_Grid = State_Grid,                                         &
            chmId      = chmId,                                              &
            Ptr2Data   = State_Chm%BoundaryCond,                             &
            nSlots     = State_Chm%nAdvect,                                  &
            RC         = RC                                                 )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( chmId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF
    ENDIF

    !========================================================================
    ! Allocate and initialize quantities that are only relevant for the
    ! the various FULLCHEM simulations or the AEROSOL-ONLY simulation
    !========================================================================
    IF ( Input_Opt%ITS_A_FULLCHEM_SIM .or. Input_Opt%ITS_AN_AEROSOL_SIM ) THEN

       ! Save nAerosol to State_Chm
       State_Chm%nAeroType = nAerosol

       !---------------------------------------------------------------------
       ! AeroArea
       !---------------------------------------------------------------------
       fieldId = (/ 'AeroAreaMDUST1   ', 'AeroAreaMDUST2   ',                &
                    'AeroAreaMDUST3   ', 'AeroAreaMDUST4   ',                &
                    'AeroAreaMDUST5   ', 'AeroAreaMDUST6   ',                &
                    'AeroAreaMDUST7   ', 'AeroAreaSULF     ',                &
                    'AeroAreaBC       ', 'AeroAreaOC       ',                &
                    'AeroAreaSSA      ', 'AeroAreaSSC      ',                &
                    'AeroAreaBGSULF   ', 'AeroAreaICEI     '                /)

       ! Allocate and register each field individually
       DO N = 1, State_Chm%nAeroType
          CALL Init_and_Register(                                            &
               Input_Opt  = Input_Opt,                                       &
               State_Chm  = State_Chm,                                       &
               State_Grid = State_Grid,                                      &
               chmId      = TRIM( fieldId(N) ),                              &
               Ptr2Data   = State_Chm%AeroArea,                              &
               nSlots     = State_Chm%nAeroType,                             &
               nCat       = N,                                               &
               RC         = RC                                              )

          IF ( RC /= GC_SUCCESS ) THEN
             errMsg = TRIM( errMsg_ir ) // TRIM( chmId )
             CALL GC_Error( errMsg, RC, thisLoc )
             RETURN
          ENDIF
        ENDDO

       !---------------------------------------------------------------------
       ! AeroRadi
       !---------------------------------------------------------------------
       fieldId = (/ 'AeroRadiMDUST1   ', 'AeroRadiMDUST2   ',                &
                    'AeroRadiMDUST3   ', 'AeroRadiMDUST4   ',                &
                    'AeroRadiMDUST5   ', 'AeroRadiMDUST6   ',                &
                    'AeroRadiMDUST7   ', 'AeroRadiSULF     ',                &
                    'AeroRadiBC       ', 'AeroRadiOC       ',                &
                    'AeroRadiSSA      ', 'AeroRadiSSC      ',                &
                    'AeroRadiBGSULF   ', 'AeroRadiICEI     '               /)

       ! Allocate and register each field individually
       DO N = 1, State_Chm%nAeroType
          CALL Init_and_Register(                                            &
               Input_Opt  = Input_Opt,                                       &
               State_Chm  = State_Chm,                                       &
               State_Grid = State_Grid,                                      &
               chmId      = TRIM( fieldId(N) ),                              &
               Ptr2Data   = State_Chm%AeroRadi,                              &
               nSlots     = State_Chm%nAeroType,                             &
               nCat       = N,                                               &
               RC         = RC                                              )

          IF ( RC /= GC_SUCCESS ) THEN
             errMsg = TRIM( errMsg_ir ) // TRIM( chmId )
             CALL GC_Error( errMsg, RC, thisLoc )
             RETURN
          ENDIF
       ENDDO

       !---------------------------------------------------------------------
       ! WetAeroArea
       !---------------------------------------------------------------------
       fieldId = (/ 'WetAeroAreaMDUST1', 'WetAeroAreaMDUST2',                &
                    'WetAeroAreaMDUST3', 'WetAeroAreaMDUST4',                &
                    'WetAeroAreaMDUST5', 'WetAeroAreaMDUST6',                &
                    'WetAeroAreaMDUST7', 'WetAeroAreaSULF  ',                &
                    'WetAeroAreaBC    ', 'WetAeroAreaOC    ',                &
                    'WetAeroAreaSSA   ', 'WetAeroAreaSSC   ',                &
                    'WetAeroAreaBGSULF', 'WetAeroAreaICEI  '               /)

       ! Allocate and register each field individually
       DO N = 1, State_Chm%nAeroType
          CALL Init_and_Register(                                            &
               Input_Opt  = Input_Opt,                                       &
               State_Chm  = State_Chm,                                       &
               State_Grid = State_Grid,                                      &
               chmId      = TRIM( fieldId(N) ),                              &
               Ptr2Data   = State_Chm%WetAeroArea,                           &
               nSlots     = State_Chm%nAeroType,                             &
               nCat       = N,                                               &
               RC         = RC                                              )

          IF ( RC /= GC_SUCCESS ) THEN
             errMsg = TRIM( errMsg_ir ) // TRIM( chmId )
             CALL GC_Error( errMsg, RC, thisLoc )
             RETURN
          ENDIF
       ENDDO

       !---------------------------------------------------------------------
       ! WetAeroRadi
       !---------------------------------------------------------------------
       fieldId = (/ 'WetAeroRadiMDUST1', 'WetAeroRadiMDUST2',                &
                    'WetAeroRadiMDUST3', 'WetAeroRadiMDUST4',                &
                    'WetAeroRadiMDUST5', 'WetAeroRadiMDUST6',                &
                    'WetAeroRadiMDUST7', 'WetAeroRadiSULF  ',                &
                    'WetAeroRadiBC    ', 'WetAeroRadiOC    ',                &
                    'WetAeroRadiSSA   ', 'WetAeroRadiSSC   ',                &
                    'WetAeroRadiBGSULF', 'WetAeroRadiICEI  '               /)

       ! Allocate and register each field individually
       DO N = 1, State_Chm%nAeroType
          CALL Init_and_Register(                                            &
               Input_Opt  = Input_Opt,                                       &
               State_Chm  = State_Chm,                                       &
               State_Grid = State_Grid,                                      &
               chmId      = TRIM( fieldId(N) ),                              &
               Ptr2Data   = State_Chm%WetAeroRadi,                           &
               nSlots     = State_Chm%nAeroType,                             &
               nCat       = N,                                               &
               RC         = RC                                              )

          IF ( RC /= GC_SUCCESS ) THEN
             errMsg = TRIM( errMsg_ir ) // TRIM( chmId )
             CALL GC_Error( errMsg, RC, thisLoc )
             RETURN
          ENDIF
       ENDDO

       !---------------------------------------------------------------------
       ! AeroH2O
       !---------------------------------------------------------------------
       fieldId = (/ 'AeroH2OMDUST1    ', 'AeroH2OMDUST2    ',                &
                    'AeroH2OMDUST3    ', 'AeroH2OMDUST4    ',                &
                    'AeroH2OMDUST5    ', 'AeroH2OMDUST6    ',                &
                    'AeroH2OMDUST7    ', 'AeroH2OSNA       ',                &
                    'AeroH2OBC        ', 'AeroH2OOC        ',                &
                    'AeroH2OSSA       ', 'AeroH2OSSC       ',                &
                    'AeroH2OBGSULF    ', 'AeroH2OICEI      '               /)

       ! Allocate and register each field individually
       DO N = 1, State_Chm%nAeroType
          CALL Init_and_Register(                                            &
               Input_Opt  = Input_Opt,                                       &
               State_Chm  = State_Chm,                                       &
               State_Grid = State_Grid,                                      &
               chmId      = TRIM( fieldId(N) ),                              &
               Ptr2Data   = State_Chm%AeroH2O,                               &
               nSlots     = State_Chm%nAeroType,                             &
               nCat       = N,                                               &
               RC         = RC                                              )

          IF ( RC /= GC_SUCCESS ) THEN
             errMsg = TRIM( errMsg_ir ) // TRIM( chmId )
             CALL GC_Error( errMsg, RC, thisLoc )
             RETURN
          ENDIF
       ENDDO

       !---------------------------------------------------------------------
       ! SoilDust, tmf 10/26/21
       !---------------------------------------------------------------------
       fieldId(1) = 'SoilDUST1'
       fieldId(2) = 'SoilDUST2'
       fieldId(3) = 'SoilDUST3'
       fieldId(4) = 'SoilDUST4'
       fieldId(5) = 'SoilDUST5'
       fieldId(6) = 'SoilDUST6'
       fieldId(7) = 'SoilDUST7'

       ! Allocate and register each field individually
       DO N = 1, NDUST
          CALL Init_and_Register(                                               &
               Input_Opt  = Input_Opt,                                          &
               State_Chm  = State_Chm,                                          &
               State_Grid = State_Grid,                                         &
               chmId      = TRIM( fieldId(N) ),                                 &
               Ptr2Data   = State_Chm%SoilDust,                                 &
               nSlots     = NDUST,                                              &
               nCat       = N,                                                  &
               RC         = RC                                                 )

          IF ( RC /= GC_SUCCESS ) THEN
             errMsg = TRIM( errMsg_ir ) // TRIM( chmId )
             CALL GC_Error( errMsg, RC, thisLoc )
             RETURN
          ENDIF
       ENDDO

       !---------------------------------------------------------------------
       ! AClArea, xnw 1/20/18
       !---------------------------------------------------------------------
       chmID = 'AClArea'
       CALL Init_and_Register(                                               &
            Input_Opt  = Input_Opt,                                          &
            State_Chm  = State_Chm,                                          &
            State_Grid = State_Grid,                                         &
            chmId      = chmId,                                              &
            Ptr2Data   = State_Chm%AClArea,                                  &
            RC         = RC                                                 )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( chmId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       !---------------------------------------------------------------------
       ! AClRadi, xnw 1/20/18
       !---------------------------------------------------------------------
       chmID = 'AClRadi'
       CALL Init_and_Register(                                               &
            Input_Opt  = Input_Opt,                                          &
            State_Chm  = State_Chm,                                          &
            State_Grid = State_Grid,                                         &
            chmId      = chmId,                                              &
            Ptr2Data   = State_Chm%AClRadi,                                  &
            RC         = RC                                                 )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( chmId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       !---------------------------------------------------------------------
       ! GammaN2O5
       !---------------------------------------------------------------------
       fieldId(1) = 'GammaN2O5overall '
       fieldId(2) = 'GammaN2O5fine    '
       fieldId(3) = 'YieldClNO2fine   '

       ! Allocate and register each field individually
       DO N = 1, 3
          CALL Init_and_Register(                                            &
               Input_Opt  = Input_Opt,                                       &
               State_Chm  = State_Chm,                                       &
               State_Grid = State_Grid,                                      &
               chmId      = TRIM( fieldId(N) ),                              &
               Ptr2Data   = State_Chm%GammaN2O5,                             &
               nSlots     = 3,                                               &
               nCat       = N,                                               &
               RC         = RC                                              )

          IF ( RC /= GC_SUCCESS ) THEN
             errMsg = TRIM( errMsg_ir ) // TRIM( chmId )
             CALL GC_Error( errMsg, RC, thisLoc )
             RETURN
          ENDIF
       ENDDO

       !---------------------------------------------------------------------
       ! OM:OC Ratios
       !---------------------------------------------------------------------
       chmId = 'OMOC'
       CALL Init_and_Register(                                               &
            Input_Opt  = Input_Opt,                                          &
            State_Chm  = State_Chm,                                          &
            State_Grid = State_Grid,                                         &
            chmId      = chmId,                                              &
            Ptr2Data   = State_Chm%OMOC,                                     &
            RC         = RC                                                 )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( chmId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       chmId = 'OMOCpoa'
       CALL Init_and_Register(                                               &
            Input_Opt  = Input_Opt,                                          &
            State_Chm  = State_Chm,                                          &
            State_Grid = State_Grid,                                         &
            chmId      = chmId,                                              &
            Ptr2Data   = State_Chm%OMOC_POA,                                 &
            RC         = RC                                                 )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( chmId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       chmId = 'OMOCopoa'
       CALL Init_and_Register(                                               &
            Input_Opt  = Input_Opt,                                          &
            State_Chm  = State_Chm,                                          &
            State_Grid = State_Grid,                                         &
            chmId      = chmId,                                              &
            Ptr2Data   = State_Chm%OMOC_OPOA,                                &
            RC         = RC                                                 )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( chmId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       !---------------------------------------------------------------------
       ! IsorropAeropH
       !---------------------------------------------------------------------
       fieldId(1) = 'IsorropAeropHAccum'
       fieldId(2) = 'IsorropAeropHCoarse'

       ! Allocate and register each field individually
       DO N = 1, 2
          CALL Init_and_Register(                                            &
               Input_Opt  = Input_Opt,                                       &
               State_Chm  = State_Chm,                                       &
               State_Grid = State_Grid,                                      &
               chmId      = TRIM( fieldId(N) ),                              &
               Ptr2Data   = State_Chm%IsorropAeropH,                         &
               nSlots     = 2,                                               &
               nCat       = N,                                               &
               RC         = RC                                              )

          IF ( RC /= GC_SUCCESS ) THEN
             errMsg = TRIM( errMsg_ir ) // TRIM( chmId )
             CALL GC_Error( errMsg, RC, thisLoc )
             RETURN
          ENDIF
       ENDDO

       !---------------------------------------------------------------------
       ! IsorropHplus
       !---------------------------------------------------------------------
       fieldId(1) = 'IsorropHplusAccum'
       fieldId(2) = 'IsorropHplusCoarse'

       ! Allocate and register each field individually
       DO N = 1, 2
          CALL Init_and_Register(                                            &
               Input_Opt  = Input_Opt,                                       &
               State_Chm  = State_Chm,                                       &
               State_Grid = State_Grid,                                      &
               chmId      = TRIM( fieldId(N) ),                              &
               Ptr2Data   = State_Chm%IsorropHplus,                          &
               nSlots     = 2,                                               &
               nCat       = N,                                               &
               RC         = RC                                              )

          IF ( RC /= GC_SUCCESS ) THEN
             errMsg = TRIM( errMsg_ir ) // TRIM( chmId )
             CALL GC_Error( errMsg, RC, thisLoc )
             RETURN
          ENDIF
       ENDDO

       !---------------------------------------------------------------------
       ! IsorropAeroH2O
       !---------------------------------------------------------------------
       fieldId(1) = 'IsorropAeroH2OAccum'
       fieldId(2) = 'IsorropAeroH2OCoarse'

       ! Allocate and register each field individually
       DO N = 1, 2
          CALL Init_and_Register(                                            &
               Input_Opt  = Input_Opt,                                       &
               State_Chm  = State_Chm,                                       &
               State_Grid = State_Grid,                                      &
               chmId      = TRIM( fieldId(N) ),                              &
               Ptr2Data   = State_Chm%IsorropAeroH2O,                        &
               nSlots     = 2,                                               &
               nCat       = N,                                               &
               RC         = RC                                              )

          IF ( RC /= GC_SUCCESS ) THEN
             errMsg = TRIM( errMsg_ir ) // TRIM( chmId )
             CALL GC_Error( errMsg, RC, thisLoc )
             RETURN
          ENDIF
       ENDDO

       !---------------------------------------------------------------------
       ! IsorropSulfate
       !---------------------------------------------------------------------
       chmId = 'IsorropSulfate'
       CALL Init_and_Register(                                               &
            Input_Opt  = Input_Opt,                                          &
            State_Chm  = State_Chm,                                          &
            State_Grid = State_Grid,                                         &
            chmId      = chmId,                                              &
            Ptr2Data   = State_Chm%IsorropSulfate,                           &
            RC         = RC                                                 )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( chmId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF


       !---------------------------------------------------------------------
       ! IsorropNitrate
       !---------------------------------------------------------------------
       fieldId(1) = 'IsorropNitrateAccum'
       fieldId(2) = 'IsorropNitrateCoarse'

       ! Allocate and register each field individually
       DO N = 1, 2
          CALL Init_and_Register(                                            &
               Input_Opt  = Input_Opt,                                       &
               State_Chm  = State_Chm,                                       &
               State_Grid = State_Grid,                                      &
               chmId      = TRIM( fieldId(N) ),                              &
               Ptr2Data   = State_Chm%IsorropNitrate,                        &
               nSlots     = 2,                                               &
               nCat       = N,                                               &
               RC         = RC                                              )

          IF ( RC /= GC_SUCCESS ) THEN
             errMsg = TRIM( errMsg_ir ) // TRIM( chmId )
             CALL GC_Error( errMsg, RC, thisLoc )
             RETURN
          ENDIF
       ENDDO

       !---------------------------------------------------------------------
       ! IsorropChloride
       !---------------------------------------------------------------------
       fieldId(1) = 'IsorropChlorideAccum'
       fieldId(2) = 'IsorropChlorideCoarse'

       ! Allocate and register each field individually
       DO N = 1, 2
          CALL Init_and_Register(                                            &
               Input_Opt  = Input_Opt,                                       &
               State_Chm  = State_Chm,                                       &
               State_Grid = State_Grid,                                      &
               chmId      = TRIM( fieldId(N) ),                              &
               Ptr2Data   = State_Chm%IsorropChloride,                       &
               nSlots     = 2,                                               &
               nCat       = N,                                               &
               RC         = RC                                              )

          IF ( RC /= GC_SUCCESS ) THEN
             errMsg = TRIM( errMsg_ir ) // TRIM( chmId )
             CALL GC_Error( errMsg, RC, thisLoc )
             RETURN
          ENDIF
       ENDDO

       !---------------------------------------------------------------------
       ! IsorropBisulfate
       !---------------------------------------------------------------------
       chmId  = 'IsorropBisulfate'
       CALL Init_and_Register(                                               &
            Input_Opt  = Input_Opt,                                          &
            State_Chm  = State_Chm,                                          &
            State_Grid = State_Grid,                                         &
            chmId      = chmId,                                              &
            Ptr2Data   = State_Chm%IsorropBisulfate,                         &
            RC         = RC                                                 )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( chmId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       !---------------------------------------------------------------------
       ! QLxphCloud
       !---------------------------------------------------------------------
       chmId = 'QLxpHCloud'
       CALL Init_and_Register(                                               &
            Input_Opt  = Input_Opt,                                          &
            State_Chm  = State_Chm,                                          &
            State_Grid = State_Grid,                                         &
            chmId      = chmId,                                              &
            Ptr2Data   = State_Chm%QLxpHCloud,                               &
            RC         = RC                                                 )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( chmId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       !---------------------------------------------------------------------
       ! ORVCsesq
       !---------------------------------------------------------------------
       chmId = 'ORVCsesq'
       CALL Init_and_Register(                                            &
            Input_Opt  = Input_Opt,                                       &
            State_Chm  = State_Chm,                                       &
            State_Grid = State_Grid,                                      &
            chmId      = chmId,                                           &
            Ptr2Data   = State_Chm%ORVCsesq,                              &
            RC         = RC                                              )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( chmId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       !---------------------------------------------------------------------
       ! isCloud (jmm 3/1/19)
       !---------------------------------------------------------------------
       chmId = 'isCloud'
       CALL Init_and_Register(                                               &
            Input_Opt  = Input_Opt,                                          &
            State_Chm  = State_Chm,                                          &
            State_Grid = State_Grid,                                         &
            chmId      = chmId,                                              &
            Ptr2Data   = State_Chm%isCloud,                                  &
            RC         = RC                                                 )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( chmId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       !---------------------------------------------------------------------
       ! SSAlk
       !---------------------------------------------------------------------
       fieldId(1) = 'SSAlkAccumMode'
       fieldId(2) = 'SSAlkCoarseMode'

       ! Allocate and register each field individually
       DO N = 1, 2
          CALL Init_and_Register(                                            &
               Input_Opt  = Input_Opt,                                       &
               State_Chm  = State_Chm,                                       &
               State_Grid = State_Grid,                                      &
               chmId      = TRIM( fieldId(N) ),                              &
               Ptr2Data   = State_Chm%SSAlk,                                 &
               nSlots     = 2,                                               &
               nCat       = N,                                               &
               RC         = RC                                              )

          IF ( RC /= GC_SUCCESS ) THEN
             errMsg = TRIM( errMsg_ir ) // TRIM( chmId )
             CALL GC_Error( errMsg, RC, thisLoc )
             RETURN
          ENDIF
       ENDDO

       !---------------------------------------------------------------------
       ! HSO3_AQ
       !---------------------------------------------------------------------
       chmId = 'HSO3AQ'
       CALL Init_and_Register(                                               &
            Input_Opt  = Input_Opt,                                          &
            State_Chm  = State_Chm,                                          &
            State_Grid = State_Grid,                                         &
            chmId      = chmId,                                              &
            Ptr2Data   = State_Chm%HSO3_AQ,                                  &
            RC         = RC                                                 )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( chmId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       !---------------------------------------------------------------------
       ! SO3_AQ
       !---------------------------------------------------------------------
       chmId = 'SO3AQ'
       CALL Init_and_Register(                                               &
            Input_Opt  = Input_Opt,                                          &
            State_Chm  = State_Chm,                                          &
            State_Grid = State_Grid,                                         &
            chmId      = chmId,                                              &
            Ptr2Data   = State_Chm%SO3_AQ,                                   &
            RC         = RC                                                 )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( chmId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       !---------------------------------------------------------------------
       ! fupdateHOBr
       !---------------------------------------------------------------------
       chmId = 'fupdateHOBr'
       CALL Init_and_Register(                                               &
            Input_Opt  = Input_Opt,                                          &
            State_Chm  = State_Chm,                                          &
            State_Grid = State_Grid,                                         &
            chmId      = chmId,                                              &
            Ptr2Data   = State_Chm%fupdateHOBr,                              &
            RC         = RC                                                 )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( chmId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       !---------------------------------------------------------------------
       ! fupdateHOCl
       !---------------------------------------------------------------------
       chmId = 'fupdateHOCl'
       CALL Init_and_Register(                                               &
            Input_Opt  = Input_Opt,                                          &
            State_Chm  = State_Chm,                                          &
            State_Grid = State_Grid,                                         &
            chmId      = chmId,                                              &
            Ptr2Data   = State_Chm%fupdateHOCl,                              &
            RC         = RC                                                 )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( chmId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       !---------------------------------------------------------------------
       ! DryDepNitrogen
       !---------------------------------------------------------------------
       chmId = 'DryDepNitrogen'
       CALL Init_and_Register(                                               &
            Input_Opt  = Input_Opt,                                          &
            State_Chm  = State_Chm,                                          &
            State_Grid = State_Grid,                                         &
            chmId      = chmId,                                              &
            Ptr2Data   = State_Chm%DryDepNitrogen,                           &
            RC         = RC                                                 )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( chmId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       !---------------------------------------------------------------------
       ! WetDepNitrogen
       !---------------------------------------------------------------------
       chmId = 'WetDepNitrogen'
       CALL Init_and_Register(                                               &
            Input_Opt  = Input_Opt,                                          &
            State_Chm  = State_Chm,                                          &
            State_Grid = State_Grid,                                         &
            chmId      = chmId,                                              &
            Ptr2Data   = State_Chm%WetDepNitrogen,                           &
            RC         = RC                                                 )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( chmId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       !------------------------------------------------------------------------
       ! TOMS_MOD
       ! Not registered to the registry as these are fields internal to the
       ! toms_mod module state.
       !------------------------------------------------------------------------
       chmId = 'TO3_DAILY'
       CALL Init_and_Register(                                               &
            Input_Opt  = Input_Opt,                                          &
            State_Chm  = State_Chm,                                          &
            State_Grid = State_Grid,                                         &
            chmId      = chmId,                                              &
            Ptr2Data   = State_Chm%TO3_DAILY,                                &
            noRegister = .TRUE.,                                             &
            RC         = RC                                                 )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( chmId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       chmId = 'TOMS1'
       CALL Init_and_Register(                                               &
            Input_Opt  = Input_Opt,                                          &
            State_Chm  = State_Chm,                                          &
            State_Grid = State_Grid,                                         &
            chmId      = chmId,                                              &
            Ptr2Data   = State_Chm%TOMS1,                                    &
            noRegister = .TRUE.,                                             &
            RC         = RC                                                 )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( chmId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       chmId = 'TOMS2'
       CALL Init_and_Register(                                               &
            Input_Opt  = Input_Opt,                                          &
            State_Chm  = State_Chm,                                          &
            State_Grid = State_Grid,                                         &
            chmId      = chmId,                                              &
            Ptr2Data   = State_Chm%TOMS2,                                    &
            noRegister = .TRUE.,                                             &
            RC         = RC                                                 )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( chmId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

    ENDIF ! ITS_A_FULLCHEM_SUM or ITS_AN_AEROSOL_SIM

    !========================================================================
    ! Allocate and initialize KPPHvalue (used by KPP-based simulations)
    !========================================================================
    IF ( Input_Opt%ITS_A_FULLCHEM_SIM      .or.                              &
         Input_Opt%ITS_A_MERCURY_SIM     ) THEN

       chmId = 'KPPHvalue'
       CALL Init_and_Register(                                               &
            Input_Opt  = Input_Opt,                                          &
            State_Chm  = State_Chm,                                          &
            State_Grid = State_Grid,                                         &
            chmId      = chmId,                                              &
            Ptr2Data   = State_Chm%KPPHvalue,                                &
            RC         = RC                                                 )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( chmId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF
    ENDIF

    !========================================================================
    ! Allocate and initialize fields for FULLCHEM or MERCURY simulations
    !========================================================================
    IF ( Input_Opt%ITS_A_FULLCHEM_SIM .or. Input_Opt%ITS_A_MERCURY_SIM ) THEN

       !---------------------------------------------------------------------
       ! STATE_PSC (polar stratospheric clouds)
       !---------------------------------------------------------------------
       chmId = 'StatePSC'
       CALL Init_and_Register(                                            &
            Input_Opt  = Input_Opt,                                       &
            State_Chm  = State_Chm,                                       &
            State_Grid = State_Grid,                                      &
            chmId      = chmId,                                           &
            Ptr2Data   = State_Chm%STATE_PSC,                             &
            RC         = RC                                              )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( chmId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       !---------------------------------------------------------------------
       ! KHETI_SLA
       !---------------------------------------------------------------------
       fieldId = (/ 'KhetiSLAN2O5H2O  ', 'KhetiSLAN2O5HCl  ',             &
                    'KhetiSLAClNO3H2O ', 'KhetiSLAClNO3HCl ',             &
                    'KhetiSLAClNO3HBr ', 'KhetiSLABrNO3H2O ',             &
                    'KhetiSLABrNO3HCl ', 'KhetiSLAHOClHCl  ',             &
                    'KhetiSLAHOClHBr  ', 'KhetiSLAHOBrHCl  ',             &
                    'KhetiSLAHOBrHBr  ', '                 ',             &
                    '                 ', '                 '            /)

       ! Allocate and register each field individually
       nKHLSA = 11
       DO N = 1, nKHLSA
          CALL Init_and_Register(                                         &
               Input_Opt  = Input_Opt,                                    &
               State_Chm  = State_Chm,                                    &
               State_Grid = State_Grid,                                   &
               chmId      = TRIM( fieldId(N) ),                           &
               Ptr2Data   = State_Chm%KHETI_SLA,                          &
               nSlots     = nKHLSA,                                       &
               nCat       = N,                                            &
               RC         = RC                                           )

          IF ( RC /= GC_SUCCESS ) THEN
             errMsg = TRIM( errMsg_ir ) // TRIM( chmId )
             CALL GC_Error( errMsg, RC, thisLoc )
             RETURN
          ENDIF
       ENDDO

       !------------------------------------------------------------------------
       ! TOMS_MOD
       ! Not registered to the registry as these are fields internal to the
       ! toms_mod module state.
       !------------------------------------------------------------------------
       chmId = 'TO3_DAILY'
       CALL Init_and_Register(                                               &
            Input_Opt  = Input_Opt,                                          &
            State_Chm  = State_Chm,                                          &
            State_Grid = State_Grid,                                         &
            chmId      = chmId,                                              &
            Ptr2Data   = State_Chm%TO3_DAILY,                                &
            noRegister = .TRUE.,                                             &
            RC         = RC                                                 )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( chmId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       chmId = 'TOMS1'
       CALL Init_and_Register(                                               &
            Input_Opt  = Input_Opt,                                          &
            State_Chm  = State_Chm,                                          &
            State_Grid = State_Grid,                                         &
            chmId      = chmId,                                              &
            Ptr2Data   = State_Chm%TOMS1,                                    &
            noRegister = .TRUE.,                                             &
            RC         = RC                                                 )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( chmId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       chmId = 'TOMS2'
       CALL Init_and_Register(                                               &
            Input_Opt  = Input_Opt,                                          &
            State_Chm  = State_Chm,                                          &
            State_Grid = State_Grid,                                         &
            chmId      = chmId,                                              &
            Ptr2Data   = State_Chm%TOMS2,                                    &
            noRegister = .TRUE.,                                             &
            RC         = RC                                                 )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( chmId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF
    ENDIF

    !========================================================================
    ! Allocate and initialize quantities for wet deposition routines
    !========================================================================

    !------------------------------------------------------------------------
    ! H2O2AfterChem
    !------------------------------------------------------------------------
    chmId = 'H2O2AfterChem'
    CALL Init_and_Register(                                                  &
         Input_Opt  = Input_Opt,                                             &
         State_Chm  = State_Chm,                                             &
         State_Grid = State_Grid,                                            &
         chmId      = chmId,                                                 &
         Ptr2Data   = State_Chm%H2O2AfterChem,                               &
         RC         = RC                                                    )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( chmId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !------------------------------------------------------------------------
    ! SO2AfterChem
    !------------------------------------------------------------------------
    chmId = 'SO2AfterChem'
    CALL Init_and_Register(                                                  &
         Input_Opt  = Input_Opt,                                             &
         State_Chm  = State_Chm,                                             &
         State_Grid = State_Grid,                                            &
         chmId      = chmId,                                                 &
         Ptr2Data   = State_Chm%SO2AfterChem,                                &
         RC         = RC                                                    )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( chmId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !------------------------------------------------------------------------
    ! phCloud
    !------------------------------------------------------------------------
    chmId = 'pHCloud'
    CALL Init_and_Register(                                                  &
         Input_Opt  = Input_Opt,                                             &
         State_Chm  = State_Chm,                                             &
         State_Grid = State_Grid,                                            &
         chmId      = chmId,                                                 &
         Ptr2Data   = State_Chm%pHcloud,                                     &
         RC         = RC                                                    )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( chmId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    ! Set default pHcloud value to 5.6, which is typical values of cloud
    ! water pH in the atmosphere.  This pH value reflects dissolved CO2
    ! in cloud water.  See geoschem/geos-chem Pull Request #779.
    State_Chm%pHCloud = 5.6_fp

#ifdef LUO_WETDEP
    !------------------------------------------------------------------------
    ! Gan Luo et al 2020 wetdep fields
    !------------------------------------------------------------------------
    IF ( Input_Opt%LWETD .or. Input_Opt%LCONV ) THEN

       ! %%% QQ3D %%%
       chmId = 'QQ3D'
       CALL Init_and_Register(                                               &
            Input_Opt  = Input_Opt,                                          &
            State_Chm  = State_Chm,                                          &
            State_Grid = State_Grid,                                         &
            chmId      = chmId,                                              &
            Ptr2Data   = State_Chm%QQ3D,                                     &
            RC         = RC                                                 )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( chmId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       ! %%% KRATE %%%
       chmId = 'KRATE'
       CALL Init_and_Register(                                               &
            Input_Opt  = Input_Opt,                                          &
            State_Chm  = State_Chm,                                          &
            State_Grid = State_Grid,                                         &
            chmId      = chmId,                                              &
            Ptr2Data   = State_Chm%KRATE,                                    &
            RC         = RC                                                 )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( chmId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       ! %%% phRain %%%
       chmId = 'pHrain'
       CALL Init_and_Register(                                               &
            Input_Opt  = Input_Opt,                                          &
            State_Chm  = State_Chm,                                          &
            State_Grid = State_Grid,                                         &
            chmId      = chmId,                                              &
            Ptr2Data   = State_Chm%pHrain,                                   &
            RC         = RC                                                 )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( chmId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       ! %%% QQpHrain %%%
       chmId = 'QQpHrain'
       CALL Init_and_Register(                                               &
            Input_Opt  = Input_Opt,                                          &
            State_Chm  = State_Chm,                                          &
            State_Grid = State_Grid,                                         &
            chmId      = chmId,                                              &
            Ptr2Data   = State_Chm%QQpHrain,                                 &
            RC         = RC                                                 )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( chmId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       ! %%% QQrain %%%
       chmId = 'QQrain'
       CALL Init_and_Register(                                               &
            Input_Opt  = Input_Opt,                                          &
            State_Chm  = State_Chm,                                          &
            State_Grid = State_Grid,                                         &
            chmId      = chmId,                                              &
            Ptr2Data   = State_Chm%QQrain,                                   &
            RC         = RC                                                 )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( chmId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF
    ENDIF
#endif

    !========================================================================
    ! Allocate fields for various GeosCore modules
    !========================================================================

    !------------------------------------------------------------------------
    ! Ocean surface iodide
    !------------------------------------------------------------------------
    IF ( State_Chm%nDryDep > 0 ) THEN
        chmId = 'Iodide'
        CALL Init_and_Register(                                              &
            Input_Opt  = Input_Opt,                                          &
            State_Chm  = State_Chm,                                          &
            State_Grid = State_Grid,                                         &
            chmId      = chmId,                                              &
            Ptr2Data   = State_Chm%Iodide,                                   &
            RC         = RC                                                 )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( chmId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF
    ENDIF

    !------------------------------------------------------------------------
    ! Ocean surface salinity
    !------------------------------------------------------------------------
    IF ( State_Chm%nDryDep > 0 ) THEN
        chmId = 'Salinity'
        CALL Init_and_Register(                                              &
            Input_Opt  = Input_Opt,                                          &
            State_Chm  = State_Chm,                                          &
            State_Grid = State_Grid,                                         &
            chmId      = chmId,                                              &
            Ptr2Data   = State_Chm%Salinity,                                 &
            RC         = RC                                                 )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( chmId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF
    ENDIF

    !------------------------------------------------------------------------
    ! DryDepFreq
    !------------------------------------------------------------------------
    IF ( State_Chm%nDryDep > 0 ) THEN
        chmId = 'DryDepFreq'
        CALL Init_and_Register(                                              &
            Input_Opt  = Input_Opt,                                          &
            State_Chm  = State_Chm,                                          &
            State_Grid = State_Grid,                                         &
            chmId      = chmId,                                              &
            Ptr2Data   = State_Chm%DryDepFreq,                               &
            nSlots     = State_Chm%nDryDep,                                  &
            RC         = RC                                                 )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( chmId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF
    ENDIF

    !------------------------------------------------------------------
    ! DryDepVel
    !------------------------------------------------------------------
    chmID = 'DryDepVel'
    CALL Init_and_Register(                                                  &
         Input_Opt  = Input_Opt,                                             &
         State_Chm  = State_Chm,                                             &
         State_Grid = State_Grid,                                            &
         chmId      = chmId,                                                 &
         Ptr2Data   = State_Chm%DryDepVel,                                   &
         nSlots     = State_Chm%nDryDep,                                     &
         RC         = RC                                                    )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( chmId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

#ifdef MODEL_GEOS
    !========================================================================
    ! Allocate and initialize aerodynamic resistance fields (GEOS-5 only)
    !========================================================================
    chmID = 'DryDepRa2m'
    CALL Init_and_Register(                                                  &
         Input_Opt  = Input_Opt,                                             &
         State_Chm  = State_Chm,                                             &
         State_Grid = State_Grid,                                            &
         chmId      = chmId,                                                 &
         Ptr2Data   = State_Chm%DryDepRa2m,                                  &
         RC         = RC                                                    )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( chmId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    chmID = 'DryDepRa10m'
    CALL Init_and_Register(                                                  &
         Input_Opt  = Input_Opt,                                             &
         State_Chm  = State_Chm,                                             &
         State_Grid = State_Grid,                                            &
         chmId      = chmId,                                                 &
         Ptr2Data   = State_Chm%DryDepRa10m,                                 &
         RC         = RC                                                    )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( chmId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
#endif

    ! J(OH) and J(NO2) are only used in fullchem simulations
    IF ( Input_Opt%ITS_A_FULLCHEM_SIM ) THEN

       !---------------------------------------------------------------------
       ! J(OH); needed for restart file input to HEMCO PARANOx extension
       !---------------------------------------------------------------------
       chmId = 'JOH'
       CALL Init_and_Register(                                               &
            Input_Opt  = Input_Opt,                                          &
            State_Chm  = State_Chm,                                          &
            State_Grid = State_Grid,                                         &
            chmId      = chmId,                                              &
            Ptr2Data   = State_Chm%JOH,                                      &
            RC         = RC                                                 )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( chmId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       !---------------------------------------------------------------------
       ! J(NO2); needed for restart file input to HEMCO PARANOx extension
       !---------------------------------------------------------------------
       chmId = 'JNO2'
       CALL Init_and_Register(                                               &
            Input_Opt  = Input_Opt,                                          &
            State_Chm  = State_Chm,                                          &
            State_Grid = State_Grid,                                         &
            chmId      = chmId,                                              &
            Ptr2Data   = State_Chm%JNO2,                                     &
            RC         = RC                                                 )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( chmId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF
    ENDIF

    !------------------------------------------------------------------
    ! Surface flux for non-local PBL mixing
    !------------------------------------------------------------------
    IF ( Input_Opt%LTURB .and. Input_Opt%LNLPBL ) THEN
       chmId = 'SurfaceFlux'
       CALL Init_and_Register(                                               &
            Input_Opt  = Input_Opt,                                          &
            State_Chm  = State_Chm,                                          &
            State_Grid = State_Grid,                                         &
            chmId      = chmId,                                              &
            Ptr2Data   = State_Chm%SurfaceFlux,                              &
            nSlots     = State_Chm%nAdvect,                                  &
            RC         = RC                                                 )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( chmId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF
    ENDIF

    !------------------------------------------------------------------
    ! TLSTT (Linoz)
    !------------------------------------------------------------------
    IF ( Input_Opt%LINOZ_NFIELDS > 0 ) THEN
        chmId = 'TLSTT'
        CALL Init_and_Register(                                              &
            Input_Opt  = Input_Opt,                                          &
            State_Chm  = State_Chm,                                          &
            State_Grid = State_Grid,                                         &
            chmId      = chmId,                                              &
            Ptr2Data   = State_Chm%TLSTT,                                    &
            nSlots     = Input_Opt%LINOZ_NFIELDS,                            &
            noRegister = .TRUE.,                                             &
            RC         = RC                                                 )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( chmId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF
    ENDIF

    !------------------------------------------------------------------------
    ! SFC_CH4
    ! Not registered to the registry as these are fields internal to the
    ! set_global_ch4_mod module state.
    !------------------------------------------------------------------------
    chmId = 'SFC_CH4'
    CALL Init_and_Register(                                                  &
         Input_Opt  = Input_Opt,                                             &
         State_Chm  = State_Chm,                                             &
         State_Grid = State_Grid,                                            &
         chmId      = chmId,                                                 &
         Ptr2Data   = State_Chm%SFC_CH4,                                     &
         noRegister = .TRUE.,                                                &
         RC         = RC                                                    )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( chmId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

#if defined(MODEL_CESM)
    IF ( Input_Opt%ITS_A_FULLCHEM_SIM ) THEN
       !---------------------------------------------------------------------
       ! H2SO4_PRDR: H2SO4 production rate [mol/mol] for MAM4 interface
       !---------------------------------------------------------------------
       chmId = 'H2SO4_PRDR'
       CALL Init_and_Register(                                               &
            Input_Opt  = Input_Opt,                                          &
            State_Chm  = State_Chm,                                          &
            State_Grid = State_Grid,                                         &
            chmId      = chmId,                                              &
            Ptr2Data   = State_Chm%H2SO4_PRDR,                               &
            RC         = RC                                                 )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( chmId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF
    ENDIF
#endif

    !=======================================================================
    ! Initialize State_Chm quantities pertinent to Hg simulations
    !=======================================================================
    IF ( Input_Opt%ITS_A_MERCURY_SIM ) THEN
       CALL Init_Hg_Simulation_Fields( Input_Opt, State_Chm, State_Grid,     &
                                       SpcCount,  RC                        )
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Error encountered in routine "Init_Hg_Simulation_Fields"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
    ENDIF

    !=======================================================================
    ! Initialize State_Chm quantities pertinent to CH4 simulations
    !=======================================================================
    IF ( Input_Opt%ITS_A_CH4_SIM ) THEN
        ! CH4_EMIS
        chmId = 'CH4_EMIS'
        CALL Init_and_Register(                                              &
            Input_Opt  = Input_Opt,                                          &
            State_Chm  = State_Chm,                                          &
            State_Grid = State_Grid,                                         &
            chmId      = chmId,                                              &
            Ptr2Data   = State_Chm%CH4_EMIS,                                 &
            nSlots     = 15,                                                 &
            RC         = RC                                                 )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( chmId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       ! Global OH and Cl from HEMCO input
       chmId = 'BOH'
       CALL Init_and_Register(                                               &
            Input_Opt  = Input_Opt,                                          &
            State_Chm  = State_Chm,                                          &
            State_Grid = State_Grid,                                         &
            chmId      = chmId,                                              &
            Ptr2Data   = State_Chm%BOH,                                      &
            RC         = RC                                                 )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( chmId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       chmId = 'BCl'
       CALL Init_and_Register(                                               &
            Input_Opt  = Input_Opt,                                          &
            State_Chm  = State_Chm,                                          &
            State_Grid = State_Grid,                                         &
            chmId      = chmId,                                              &
            Ptr2Data   = State_Chm%BCl,                                      &
            RC         = RC                                                 )

       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = TRIM( errMsg_ir ) // TRIM( chmId )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF
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
 10    FORMAT( /, 'Registered variables contained within the State_Chm object:')
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
! !IROUTINE: Init_Mapping_Vectors
!
! !DESCRIPTION: Initializes the 1-D mapping vectors in the State_Chm object.
!  This was split off from Init_State_Chm.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Init_Mapping_Vectors( Input_Opt, State_Chm, RC )
!
! !USES:
!
    USE GCKPP_Parameters, ONLY : NSPEC
    USE Input_Opt_Mod,    ONLY : OptInput
    USE Cmn_Fjx_Mod,      ONLY : W_
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(INOUT) :: Input_Opt   ! Input Options object
    TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure
!
! !REVISION HISTORY:
!  06 Jan 2015 - R. Yantosca - Initial version
!  See the subsequent Git history with the gitk browser!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER                :: C,      N

    ! Strings
    CHARACTER(LEN=255)     :: errMsg, thisLoc

    ! Objects
    TYPE(Species), POINTER :: ThisSpc

    !========================================================================
    ! Init_Mapping_Vectors begins here!
    !========================================================================

    ! Initialize
    RC         =  GC_SUCCESS
    ThisSpc    => NULL()
    errMsg     =  ''
    thisLoc    =  &
       ' -> at Init_Mapping_Vectors (in module Headers/state_chm_mod.F90)'

    !========================================================================
    ! Allocate and initialize mapping vectors to subset species
    !========================================================================
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

    !---------------------------------------------------------------------------
    ! NOTE: Need to also leave room for the omitted "dummy" KPP species in
    ! the mapping arrays so that the rest of the KPP indices will line up!
    !  -- Bob Yantosca (04 Jun 2021)
    N = State_Chm%nKppVar + State_Chm%nOmitted
    IF ( N > 0 ) THEN
       ALLOCATE( State_Chm%Map_KppVar( N ), STAT=RC )
       CALL GC_CheckVar( 'State_Chm%Map_KppVar', 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Chm%Map_KppVar = 0
    ENDIF

    N = State_Chm%nKppFix + State_Chm%nOmitted
    IF ( N > 0 ) THEN
       ALLOCATE( State_Chm%Map_KppFix( N ), STAT=RC )
       CALL GC_CheckVar( 'State_Chm%Map_KppFix', 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Chm%Map_KppFix = 0
    ENDIF

    N = State_Chm%nKppSpc + State_Chm%nOmitted
    IF ( N > 0 ) THEN
       ALLOCATE( State_Chm%Map_KppSpc( N ), STAT=RC )
       CALL GC_CheckVar( 'State_Chm%Map_KppSpc', 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Chm%Map_KppSpc = 0
    ENDIF
    !---------------------------------------------------------------------------

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

    IF ( State_Chm%nRadNucl > 0 ) THEN
       ALLOCATE( State_Chm%Map_RadNucl( State_Chm%nRadNucl ), STAT=RC )
       CALL GC_CheckVar( 'State_Chm%Map_RadNucl', 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Chm%Map_RadNucl = 0
    ENDIF

    IF ( State_Chm%nWetDep > 0 ) THEN
       ALLOCATE( State_Chm%Map_WetDep( State_Chm%nWetDep ), STAT=RC )
       CALL GC_CheckVar( 'State_Chm%Map_WetDep', 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Chm%Map_WetDep = 0
    ENDIF

    IF ( W_ > 0 ) THEN
       ALLOCATE( State_Chm%Map_WL( W_ ), STAT=RC )
       CALL GC_CheckVar( 'State_Chm%Map_WL', 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Chm%Map_WL = 0
    ENDIF

    !========================================================================
    ! Set up the species mapping vectors
    !========================================================================
    IF ( Input_Opt%amIRoot ) THEN
       WRITE( 6,'(/,a)' ) 'ADVECTED SPECIES MENU'
       WRITE( 6,'(  a)' ) REPEAT( '-', 48 )
       WRITE( 6,'(  a)' ) '  #  Species Name'
    ENDIF

    ! Loop over all species
    DO N = 1, State_Chm%nSpecies

       ! GEOS-Chem Species Database entry for species # N
       ThisSpc => State_Chm%SpcData(N)%Info

       !---------------------------------------------------------------------
       ! Set up the mapping for ADVECTED SPECIES
       !---------------------------------------------------------------------
       IF ( ThisSpc%Is_Advected ) THEN

          ! Update the mapping vector of advected species
          C                       = ThisSpc%AdvectId
          State_Chm%Map_Advect(C) = ThisSpc%ModelId

          ! Print to screen
          IF ( Input_Opt%amIRoot ) THEN
             WRITE( 6, 100 ) ThisSpc%ModelId, ThisSpc%Name
 100         FORMAT( I3, 2x, A31 )
          ENDIF

       ENDIF

       !---------------------------------------------------------------------
       ! Set up the mapping for AEROSOL SPECIES
       !---------------------------------------------------------------------
       IF ( ThisSpc%Is_Aerosol ) THEN
          C                     = ThisSpc%AerosolId
          State_Chm%Map_Aero(C) = ThisSpc%ModelId
       ENDIF

       !---------------------------------------------------------------------
       ! Set up the mapping for DRYDEP SPECIES TO SAVE AT A GIVEN ALTITUDE
       !---------------------------------------------------------------------
       IF ( ThisSpc%Is_DryAlt ) THEN
          C                       = ThisSpc%DryAltId
          State_Chm%Map_DryAlt(C) = ThisSpc%ModelId
       ENDIF

       !---------------------------------------------------------------------
       ! Set up the mapping for DRYDEP SPECIES
       !---------------------------------------------------------------------
       IF ( ThisSpc%Is_DryDep ) THEN
          C                       = ThisSpc%DryDepId
          State_Chm%Map_Drydep(C) = ThisSpc%ModelId
       ENDIF

       !---------------------------------------------------------------------
       ! Set up the mapping for GAS SPECIES
       !---------------------------------------------------------------------
       IF ( ThisSpc%Is_Gas ) THEN
          C                       = ThisSpc%GasSpcId
          !###
          !### Uncomment for debug print if Map_GasSpc goes out-of-bounds
          !### print*, '===> ', ThisSpc%Name, C, ThisSpc%ModelId
          !###
          State_Chm%Map_GasSpc(C) = ThisSpc%ModelId
       ENDIF

       !---------------------------------------------------------------------
       ! Set up the mapping for HYGROSCOPIC GROWTH SPECIES
       !---------------------------------------------------------------------
       IF ( ThisSpc%Is_HygroGrowth ) THEN
          C                        = ThisSpc%HygGrthId
          State_Chm%Map_HygGrth(C) = ThisSpc%ModelId
       ENDIF

       !---------------------------------------------------------------------
       ! Set up the mapping for KPP ACTIVE (VARIABLE) SPECIES
       !---------------------------------------------------------------------
       IF ( ThisSpc%Is_ActiveChem ) THEN
          C                       = ThisSpc%KppVarId
          State_Chm%Map_KppVar(C) = ThisSpc%ModelId
       ENDIF

       !---------------------------------------------------------------------
       ! Set up the mapping for KPP FIXED SPECIES
       !---------------------------------------------------------------------
       IF ( ThisSpc%Is_FixedChem ) THEN
          C                       = ThisSpc%KppFixId
          State_Chm%Map_KppFix(C) = ThisSpc%ModelId
       ENDIF

       !---------------------------------------------------------------------
       ! Set up the mapping for SPECIES IN THE KPP MECHANISM
       !---------------------------------------------------------------------
       IF ( ThisSpc%Is_Kpp ) THEN
          C                       = ThisSpc%KppSpcId
          State_Chm%Map_KppSpc(C) = ThisSpc%ModelId
       ENDIF

       !---------------------------------------------------------------------
       ! Set up the mapping for PHOTOLYSIS SPECIES
       !---------------------------------------------------------------------
       IF ( ThisSpc%Is_Photolysis ) THEN
          C                       = ThisSpc%PhotolId
          State_Chm%Map_Photol(C) = ThisSpc%ModelId
       ENDIF

       !---------------------------------------------------------------------
       ! Set up the mapping for WETDEP SPECIES
       !---------------------------------------------------------------------
       IF ( ThisSpc%Is_RadioNuclide ) THEN
          C                        = ThisSpc%RadNuclId
          State_Chm%Map_RadNucl(C) = ThisSpc%ModelId
       ENDIF

       !---------------------------------------------------------------------
       ! Set up the mapping for WETDEP SPECIES
       !---------------------------------------------------------------------
       IF ( ThisSpc%Is_WetDep ) THEN
          C                       = ThisSpc%WetDepId
          State_Chm%Map_WetDep(C) = ThisSpc%ModelId
       ENDIF

       ! Free pointer
       ThisSpc => NULL()

    ENDDO

    !------------------------------------------------------------------------
    ! Set up the mapping for UVFlux Diagnostics
    ! placeholder for now since couldn't figure out how to read in WL from file
    !------------------------------------------------------------------------
    IF ( W_ > 0 ) THEN

       ! Define identifying string
       DO N = 1, W_
          State_Chm%Map_WL(N) = 0
       ENDDO
    ENDIF

    !------------------------------------------------------------------------
    ! Set up the mapping for PRODUCTION AND LOSS DIAGNOSTIC SPECIES
    !------------------------------------------------------------------------
    IF ( State_Chm%nProd > 0 .or. State_Chm%nLoss > 0 ) THEN
       CALL MapProdLossSpecies( Input_Opt, State_Chm, RC )
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Error encountered in routine "MapProdLossSpecies"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
    ENDIF

  END SUBROUTINE Init_Mapping_Vectors
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Init_Hg_Simulation_Fields
!
! !DESCRIPTION: Initializes State_Chm quantities that pertain to the
!  Hg or tagged Hg simulations.  This was split off from Init_State_Chm.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Init_Hg_Simulation_Fields( Input_Opt, State_Chm, State_Grid,    &
                                        SpcCount,  RC                       )
!
! !USES:
!
    USE Input_Opt_Mod,  ONLY : OptInput
    USE State_Grid_Mod, ONLY : GrdState
!
! INPUT PARAMETERS:
!
    TYPE(GrdState), INTENT(IN)    :: State_Grid
    TYPE(SpcIndCt), INTENT(IN)    :: SpcCount
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(INOUT) :: Input_Opt
    TYPE(ChmState), INTENT(INOUT) :: State_Chm
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC
!
! !REVISION HISTORY:
!  23 Jan 2015 - R. Yantosca - Initial version
!  See the subsequent Git history with the gitk browser!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER                :: C
    INTEGER                :: N

    ! Strings
    CHARACTER(LEN=255)     :: chmId
    CHARACTER(LEN=255)     :: errMsg
    CHARACTER(LEN=255)     :: errMsg_ir
    CHARACTER(LEN=255)     :: thisLoc

    ! Objects
    TYPE(Species), POINTER :: ThisSpc

    !========================================================================
    ! Init_Hg_Simulation_Fields begins here!
    !========================================================================

    ! Initialize
    RC         = GC_SUCCESS
    errMsg     = ''
    errMsg_ir  =  'Error encountered in "Init_and_Register", chmId = '
    thisLoc    = &
     ' -> at Init_Hg_Simulation_Fields (in module Headers/state_chm_mod.F90)'

    !------------------------------------------------------------------------
    ! Hg(0) ocean mass
    !------------------------------------------------------------------------
    chmId = 'OceanHg0'
    CALL Init_and_Register(                                                  &
         Input_Opt  = Input_Opt,                                             &
         State_Chm  = State_Chm,                                             &
         State_Grid = State_Grid,                                            &
         chmId      = chmId,                                                 &
         Ptr2Data   = State_Chm%OceanHg0,                                    &
         RC         = RC                                                    )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( chmId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !------------------------------------------------------------------------
    ! Hg(II) ocean mass
    !------------------------------------------------------------------------
    chmId = 'OceanHg2'
    CALL Init_and_Register(                                                  &
         Input_Opt  = Input_Opt,                                             &
         State_Chm  = State_Chm,                                             &
         State_Grid = State_Grid,                                            &
         chmId      = chmId,                                                 &
         Ptr2Data   = State_Chm%OceanHg2,                                    &
         RC         = RC                                                    )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( chmId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !------------------------------------------------------------------------
    ! HgP ocean mass
    !------------------------------------------------------------------------
    chmId = 'OceanHgP'
    CALL Init_and_Register(                                                  &
         Input_Opt  = Input_Opt,                                             &
         State_Chm  = State_Chm,                                             &
         State_Grid = State_Grid,                                            &
         chmId      = chmId,                                                 &
         Ptr2Data   = State_Chm%OceanHgP,                                    &
         RC         = RC                                                    )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( chmId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !------------------------------------------------------------------------
    ! Reducible Hg snowpack on ocean
    !------------------------------------------------------------------------
    chmId = 'SnowHgOcean'
    CALL Init_and_Register(                                                  &
         Input_Opt  = Input_Opt,                                             &
         State_Chm  = State_Chm,                                             &
         State_Grid = State_Grid,                                            &
         chmId      = chmId,                                                 &
         Ptr2Data   = State_Chm%SnowHgOcean,                                 &
         RC         = RC                                                    )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( chmId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !------------------------------------------------------------------------
    ! Reducible Hg snowpack on land
    !------------------------------------------------------------------------
    chmId = 'SnowHgLand'
    CALL Init_and_Register(                                                  &
         Input_Opt  = Input_Opt,                                             &
         State_Chm  = State_Chm,                                             &
         State_Grid = State_Grid,                                            &
         chmId      = chmId,                                                 &
         Ptr2Data   = State_Chm%SnowHgLand,                                  &
         RC         = RC                                                    )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( chmId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !------------------------------------------------------------------------
    ! Non-reducible Hg snowpack on ocean
    !------------------------------------------------------------------------
    chmId = 'SnowHgOceanStored'
    CALL Init_and_Register(                                                  &
         Input_Opt  = Input_Opt,                                             &
         State_Chm  = State_Chm,                                             &
         State_Grid = State_Grid,                                            &
         chmId      = chmId,                                                 &
         Ptr2Data   = State_Chm%SnowHgOceanStored,                           &
         RC         = RC                                                    )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( chmId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !------------------------------------------------------------------------
    ! Non-reducible Hg snowpack on land
    !------------------------------------------------------------------------
    chmId = 'SnowHgLandStored'
    CALL Init_and_Register(                                                  &
         Input_Opt  = Input_Opt,                                             &
         State_Chm  = State_Chm,                                             &
         State_Grid = State_Grid,                                            &
         chmId      = chmId,                                                 &
         Ptr2Data   = State_Chm%SnowHgLandStored,                            &
         RC         = RC                                                    )

    !------------------------------------------------------------------------
    ! TOMS_MOD
    ! Not registered to the registry as these are fields internal to the
    ! toms_mod module state.
    !------------------------------------------------------------------------
    chmId = 'TO3_DAILY'
    CALL Init_and_Register(                                               &
         Input_Opt  = Input_Opt,                                          &
         State_Chm  = State_Chm,                                          &
         State_Grid = State_Grid,                                         &
         chmId      = chmId,                                              &
         Ptr2Data   = State_Chm%TO3_DAILY,                                &
         noRegister = .TRUE.,                                             &
         RC         = RC                                                 )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( chmId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    chmId = 'TOMS1'
    CALL Init_and_Register(                                               &
         Input_Opt  = Input_Opt,                                          &
         State_Chm  = State_Chm,                                          &
         State_Grid = State_Grid,                                         &
         chmId      = chmId,                                              &
         Ptr2Data   = State_Chm%TOMS1,                                    &
         noRegister = .TRUE.,                                             &
         RC         = RC                                                 )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( chmId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    chmId = 'TOMS2'
    CALL Init_and_Register(                                               &
         Input_Opt  = Input_Opt,                                          &
         State_Chm  = State_Chm,                                          &
         State_Grid = State_Grid,                                         &
         chmId      = chmId,                                              &
         Ptr2Data   = State_Chm%TOMS2,                                    &
         noRegister = .TRUE.,                                             &
         RC         = RC                                                 )

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( chmId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = TRIM( errMsg_ir ) // TRIM( chmId )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

  END SUBROUTINE Init_Hg_Simulation_Fields
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
    INTEGER            :: N
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

    IF ( ASSOCIATED( State_Chm%Map_RadNucl ) ) THEN
       DEALLOCATE( State_Chm%Map_RadNucl, STAT=RC )
       CALL GC_CheckVar( 'State_Chm%Map_WetDep', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Chm%Map_RadNucl => NULL()
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

    IF ( ASSOCIATED ( State_Chm%Species ) ) THEN
       DO N = 1, State_Chm%nSpecies
          IF ( ASSOCIATED( State_Chm%Species(N)%Conc ) ) THEN
#if !defined( MODEL_GCHPCTM )
             DEALLOCATE( State_Chm%Species(N)%Conc, STAT=RC )
             IF ( RC /= GC_SUCCESS ) RETURN
#endif
             State_Chm%Species(N)%Conc => NULL()
          ENDIF
       ENDDO
       DEALLOCATE( State_Chm%Species )
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

    IF ( ASSOCIATED( State_Chm%AClArea ) ) THEN
       DEALLOCATE( State_Chm%AClArea, STAT=RC )
       CALL GC_CheckVar( 'State_Chm%AClArea', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Chm%AClArea => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Chm%AClRadi ) ) THEN
       DEALLOCATE( State_Chm%AClRadi, STAT=RC )
       CALL GC_CheckVar( 'State_Chm%AClRadi', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Chm%AClRadi => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Chm%SoilDust ) ) THEN
       DEALLOCATE( State_Chm%SoilDust, STAT=RC )
       CALL GC_CheckVar( 'State_Chm%SoilDust', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Chm%SoilDust => NULL()
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

    IF ( ASSOCIATED( State_Chm%OMOC ) ) THEN
       DEALLOCATE( State_Chm%OMOC, STAT=RC )
       CALL GC_CheckVar( 'State_Chm%OMOC', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Chm%OMOC => NULL()
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

    IF ( ASSOCIATED( State_Chm%IsorropAeropH ) ) THEN
       DEALLOCATE( State_Chm%IsorropAeropH, STAT=RC  )
       CALL GC_CheckVar( 'State_Chm%IsorropAeropH', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Chm%IsorropAeropH => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Chm%IsorropHplus ) ) THEN
       DEALLOCATE( State_Chm%IsorropHplus, STAT=RC )
       CALL GC_CheckVar( 'State_Chm%IsorropHplus', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Chm%IsorropHplus => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Chm%IsorropAeroH2O ) ) THEN
       DEALLOCATE( State_Chm%IsorropAeroH2O, STAT=RC )
       CALL GC_CheckVar( 'State_Chm%IsorropAeroH2O', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Chm%IsorropAeroH2O => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Chm%IsorropSulfate ) ) THEN
       DEALLOCATE( State_Chm%IsorropSulfate, STAT=RC )
       CALL GC_CheckVar( 'State_Chm%IsorropSulfate', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Chm%IsorropSulfate => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Chm%IsorropNitrate ) ) THEN
       DEALLOCATE( State_Chm%IsorropNitrate, STAT=RC )
       CALL GC_CheckVar( 'State_Chm%IsorropNitrate', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Chm%IsorropNitrate => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Chm%IsorropChloride ) ) THEN
       DEALLOCATE( State_Chm%IsorropChloride, STAT=RC )
       CALL GC_CheckVar( 'State_Chm%IsorropChloride', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Chm%IsorropChloride => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Chm%IsorropBisulfate ) ) THEN
       DEALLOCATE( State_Chm%IsorropBisulfate, STAT=RC )
       CALL GC_CheckVar( 'State_Chm%IsorropBisulfate', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Chm%IsorropBisulfate => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Chm%pHCloud ) ) THEN
       DEALLOCATE( State_Chm%pHCloud, STAT=RC )
       CALL GC_CheckVar( 'State_Chm%pHCloud', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Chm%pHCloud => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Chm%QLxpHCloud ) ) THEN
       DEALLOCATE( State_Chm%QLxpHCloud, STAT=RC )
       CALL GC_CheckVar( 'State_Chm%QLxpHCloud', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Chm%QLxpHCloud => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Chm%ORVCsesq ) ) THEN
       DEALLOCATE( State_Chm%ORVCsesq, STAT=RC )
       CALL GC_CheckVar( 'State_Chm%ORVCsesq', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Chm%ORVCsesq => NULL()
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
       State_Chm%WetDepNitrogen => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Chm%KPPHvalue ) ) THEN
       DEALLOCATE( State_Chm%KPPHvalue, STAT=RC  )
       CALL GC_CheckVar( 'State_Chm%KPPHvalue', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Chm%KPPHvalue => NULL()
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

    IF ( ASSOCIATED( State_Chm%fupdateHOCl ) ) THEN
       DEALLOCATE( State_Chm%fupdateHOCl, STAT=RC )
       CALL GC_CheckVar( 'State_Chm%fupdateHOCl', 3, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Chm%fupdateHOCl => NULL()
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

    IF ( ASSOCIATED( State_Chm%DryDepVel ) ) THEN
       DEALLOCATE( State_Chm%DryDepVel, STAT=RC )
       CALL GC_CheckVar( 'State_Chm%DryDepVel', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Chm%DryDepVel => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Chm%DryDepFreq ) ) THEN
       DEALLOCATE( State_Chm%DryDepFreq, STAT=RC )
       CALL GC_CheckVar( 'State_Chm%DryDepFreq', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Chm%DryDepFreq => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Chm%Iodide ) ) THEN
       DEALLOCATE( State_Chm%Iodide, STAT=RC )
       CALL GC_CheckVar( 'State_Chm%Iodide', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Chm%Iodide => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Chm%Salinity ) ) THEN
       DEALLOCATE( State_Chm%Salinity, STAT=RC )
       CALL GC_CheckVar( 'State_Chm%Salinity', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Chm%Salinity => NULL()
    ENDIF

#ifdef MODEL_GEOS
    IF ( ASSOCIATED( State_Chm%DryDepRa2m ) ) THEN
       DEALLOCATE( State_Chm%DryDepRa2m, STAT=RC )
       CALL GC_CheckVar( 'State_Chm%DryDepRa2m', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Chm%DryDepRa2m => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Chm%DryDepRa10m ) ) THEN
       DEALLOCATE( State_Chm%DryDepRa10m, STAT=RC )
       CALL GC_CheckVar( 'State_Chm%DryDepRa10m', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Chm%DryDepRa10m => NULL()
    ENDIF

    ! CO2 photolysis stuff
    IF ( ASSOCIATED( State_Chm%sdat ) ) THEN
       DEALLOCATE( State_Chm%sdat, STAT=RC )
       CALL GC_CheckVar( 'State_Chm%sdat', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Chm%sdat => NULL()
    ENDIF
    IF ( ASSOCIATED( State_Chm%o2jdat ) ) THEN
       DEALLOCATE( State_Chm%o2jdat, STAT=RC )
       CALL GC_CheckVar( 'State_Chm%o2jdat', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Chm%o2jdat => NULL()
    ENDIF
    IF ( ASSOCIATED( State_Chm%sza_tab ) ) THEN
       DEALLOCATE( State_Chm%sza_tab, STAT=RC )
       CALL GC_CheckVar( 'State_Chm%sza_tab', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Chm%sza_tab => NULL()
    ENDIF
    IF ( ASSOCIATED( State_Chm%o3_tab ) ) THEN
       DEALLOCATE( State_Chm%o3_tab, STAT=RC )
       CALL GC_CheckVar( 'State_Chm%o3_tab', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Chm%o3_tab => NULL()
    ENDIF
    IF ( ASSOCIATED( State_Chm%xtab ) ) THEN
       DEALLOCATE( State_Chm%xtab, STAT=RC )
       CALL GC_CheckVar( 'State_Chm%xtab', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Chm%xtab => NULL()
    ENDIF
    IF ( ASSOCIATED( State_Chm%CH2O_aq ) ) THEN
       DEALLOCATE( State_Chm%CH2O_aq, STAT=RC )
       CALL GC_CheckVar( 'State_Chm%CH2O_aq', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Chm%CH2O_aq => NULL()
    ENDIF
    IF ( ASSOCIATED( State_Chm%rlam ) ) THEN
       DEALLOCATE( State_Chm%rlam, STAT=RC )
       CALL GC_CheckVar( 'State_Chm%rlam', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Chm%rlam => NULL()
    ENDIF
#endif

    IF ( ASSOCIATED( State_Chm%JOH ) ) THEN
       DEALLOCATE( State_Chm%JOH, STAT=RC )
       CALL GC_CheckVar( 'State_Chm%JOH', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Chm%JOH => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Chm%JNO2 ) ) THEN
       DEALLOCATE( State_Chm%JNO2, STAT=RC )
       CALL GC_CheckVar( 'State_Chm%JNO2', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Chm%JNO2 => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Chm%SurfaceFlux ) ) THEN
       DEALLOCATE( State_Chm%SurfaceFlux, STAT=RC )
       CALL GC_CheckVar( 'State_Chm%SurfaceFlux', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Chm%SurfaceFlux => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Chm%TLSTT ) ) THEN
       DEALLOCATE( State_Chm%TLSTT, STAT=RC )
       CALL GC_CheckVar( 'State_Chm%TLSTT', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Chm%TLSTT => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Chm%BOH ) ) THEN
       DEALLOCATE( State_Chm%BOH, STAT=RC )
       CALL GC_CheckVar( 'State_Chm%BOH', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Chm%BOH => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Chm%BCl ) ) THEN
       DEALLOCATE( State_Chm%BCl, STAT=RC )
       CALL GC_CheckVar( 'State_Chm%BCl', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Chm%BCl => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Chm%CH4_EMIS ) ) THEN
       DEALLOCATE( State_Chm%CH4_EMIS, STAT=RC )
       CALL GC_CheckVar( 'State_Chm%CH4_EMIS', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Chm%CH4_EMIS => NULL()
    ENDIF

#ifdef LUO_WETDEP
    IF ( ASSOCIATED( State_Chm%QQ3D ) ) THEN
       DEALLOCATE( State_Chm%QQ3D, STAT=RC )
       CALL GC_CheckVar( 'State_Chm%QQ3D', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Chm%QQ3D => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Chm%KRATE ) ) THEN
       DEALLOCATE( State_Chm%KRATE, STAT=RC )
       CALL GC_CheckVar( 'State_Chm%KRATE', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Chm%KRATE => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Chm%pHrain ) ) THEN
       DEALLOCATE( State_Chm%pHrain, STAT=RC )
       CALL GC_CheckVar( 'State_Chm%pHrain', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Chm%pHrain => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Chm%QQpHrain ) ) THEN
       DEALLOCATE( State_Chm%QQpHrain, STAT=RC )
       CALL GC_CheckVar( 'State_Chm%QQpHrain', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Chm%QQpHrain => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Chm%QQrain ) ) THEN
       DEALLOCATE( State_Chm%QQrain, STAT=RC )
       CALL GC_CheckVar( 'State_Chm%QQrain', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Chm%QQrain => NULL()
    ENDIF
#endif

    IF ( ASSOCIATED( State_Chm%SFC_CH4 ) ) THEN
       DEALLOCATE( State_Chm%SFC_CH4, STAT=RC )
       CALL GC_CheckVar( 'State_Chm%SFC_CH4', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Chm%SFC_CH4 => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Chm%TO3_DAILY ) ) THEN
       DEALLOCATE( State_Chm%TO3_DAILY, STAT=RC )
       CALL GC_CheckVar( 'State_Chm%TO3_DAILY', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Chm%TO3_DAILY => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Chm%TOMS1 ) ) THEN
       DEALLOCATE( State_Chm%TOMS1, STAT=RC )
       CALL GC_CheckVar( 'State_Chm%TOMS1', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Chm%TOMS1 => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Chm%TOMS2 ) ) THEN
       DEALLOCATE( State_Chm%TOMS2, STAT=RC )
       CALL GC_CheckVar( 'State_Chm%TOMS2', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Chm%TOMS2 => NULL()
    ENDIF

#if defined(MODEL_CESM)
    IF ( ASSOCIATED( State_Chm%H2SO4_PRDR ) ) THEN
       DEALLOCATE( State_Chm%H2SO4_PRDR, STAT=RC )
       CALL GC_CheckVar( 'State_Chm%H2SO4_PRDR', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Chm%H2SO4_PRDR => NULL()
    ENDIF
#endif

    IF ( ASSOCIATED( State_Chm%RAD_AER ) ) THEN
       DEALLOCATE( State_Chm%RAD_AER, STAT=RC )
       CALL GC_CheckVar( 'State_Chm%RAD_AER', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Chm%RAD_AER => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Chm%SAD_AER ) ) THEN
       DEALLOCATE( State_Chm%SAD_AER, STAT=RC )
       CALL GC_CheckVar( 'State_Chm%SAD_AER', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Chm%SAD_AER => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Chm%KG_AER ) ) THEN
       DEALLOCATE( State_Chm%KG_AER, STAT=RC )
       CALL GC_CheckVar( 'State_Chm%KG_AER', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Chm%KG_AER => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Chm%RHO_AER ) ) THEN
       DEALLOCATE( State_Chm%RHO_AER, STAT=RC )
       CALL GC_CheckVar( 'State_Chm%RHO_AER', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Chm%RHO_AER => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Chm%NDENS_AER ) ) THEN
       DEALLOCATE( State_Chm%NDENS_AER, STAT=RC )
       CALL GC_CheckVar( 'State_Chm%NDENS_AER', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Chm%NDENS_AER => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Chm%AERFRAC ) ) THEN
       DEALLOCATE( State_Chm%AERFRAC, STAT=RC )
       CALL GC_CheckVar( 'State_Chm%AERFRAC', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Chm%AERFRAC => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Chm%AERFRACIND ) ) THEN
       DEALLOCATE( State_Chm%AERFRACIND, STAT=RC )
       CALL GC_CheckVar( 'State_Chm%AERFRACIND', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Chm%AERFRACIND => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Chm%UCX_REGRID ) ) THEN
       DEALLOCATE( State_Chm%UCX_REGRID, STAT=RC )
       CALL GC_CheckVar( 'State_Chm%UCX_REGRID', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Chm%UCX_REGRID => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Chm%UCX_PLEVS ) ) THEN
       DEALLOCATE( State_Chm%UCX_PLEVS, STAT=RC )
       CALL GC_CheckVar( 'State_Chm%UCX_PLEVS', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Chm%UCX_PLEVS => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Chm%UCX_LATS ) ) THEN
       DEALLOCATE( State_Chm%UCX_LATS, STAT=RC )
       CALL GC_CheckVar( 'State_Chm%UCX_LATS', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Chm%UCX_LATS => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Chm%NOX_O ) ) THEN
       DEALLOCATE( State_Chm%NOX_O, STAT=RC )
       CALL GC_CheckVar( 'State_Chm%NOX_O', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Chm%NOX_O => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Chm%NOX_J ) ) THEN
       DEALLOCATE( State_Chm%NOX_J, STAT=RC )
       CALL GC_CheckVar( 'State_Chm%NOX_J', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Chm%NOX_J => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Chm%SO4_TOPPHOT ) ) THEN
       DEALLOCATE( State_Chm%SO4_TOPPHOT, STAT=RC )
       CALL GC_CheckVar( 'State_Chm%SO4_TOPPHOT', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Chm%SO4_TOPPHOT => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Chm%NOXCOEFF ) ) THEN
       DEALLOCATE( State_Chm%NOXCOEFF, STAT=RC )
       CALL GC_CheckVar( 'State_Chm%NOXCOEFF', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Chm%NOXCOEFF => NULL()
    ENDIF

    IF ( ASSOCIATED( State_Chm%NOXLAT ) ) THEN
       DEALLOCATE( State_Chm%NOXLAT, STAT=RC )
       CALL GC_CheckVar( 'State_Chm%NOXLAT', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Chm%NOXLAT => NULL()
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
                                     PerSpc,     Rank,       Type,       &
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
    LOGICAL,             INTENT(OUT)           :: Found   ! Item found?
    INTEGER,             INTENT(OUT)           :: RC      ! Return code
    CHARACTER(LEN=255),  INTENT(OUT), OPTIONAL :: Desc    ! Long name string
    CHARACTER(LEN=255),  INTENT(OUT), OPTIONAL :: Units   ! Units string
    CHARACTER(LEN=255),  INTENT(OUT), OPTIONAL :: PerSpc  ! Max spc wildcard
    INTEGER,             INTENT(OUT), OPTIONAL :: Rank    ! # of dimensions
    INTEGER,             INTENT(OUT), OPTIONAL :: Type    ! Desc of data type
    INTEGER,             INTENT(OUT), OPTIONAL :: VLoc    ! Vert placement
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
    LOGICAL            :: isDesc, isUnits, isRank, isType, isVLoc, isSpc

    !=======================================================================
    ! Initialize
    !=======================================================================

    ! Assume success
    RC      =  GC_SUCCESS
    ThisLoc = ' -> at Get_Metadata_State_Chm (in Headers/state_chm_mod.F90)'
    Found   = .TRUE.

    ! Optional arguments present?
    isDesc  = PRESENT( Desc   )
    isUnits = PRESENT( Units  )
    isRank  = PRESENT( Rank   )
    isType  = PRESENT( Type   )
    isVLoc  = PRESENT( VLoc   )
    isSpc   = PRESENT( PerSpc )

    ! Set defaults for optional arguments. Assume type and vertical
    ! location are real (flexible precision) and center unless specified
    ! otherwise
    IF ( isUnits ) Units  = ''
    IF ( isDesc  ) Desc   = ''
    IF ( isRank  ) Rank   = -1              ! Init # dims as bad value
    IF ( isType  ) Type   = KINDVAL_FP      ! Assume flexible precision
    IF ( isVLoc  ) VLoc   = VLocationCenter ! Assume vertically centered
    IF ( isSpc   ) PerSpc = ''              ! Assume not per species

    ! Convert name to uppercase
    Name_AllCaps = To_Uppercase( TRIM( metadataID ) )

    !=======================================================================
    ! Values for Retrieval (string comparison slow but happens only once)
    !=======================================================================
    SELECT CASE ( TRIM( Name_AllCaps ) )

       CASE ( 'SPECIES' )
          IF ( isDesc  ) Desc   = 'Concentration for species'
          IF ( isUnits ) Units  = 'varies'
          IF ( isRank  ) Rank   = 3
          IF ( isSpc   ) PerSpc = 'ALL'

#ifdef ADJOINT
       CASE ( 'SPECIESADJ' )
          IF ( isDesc  ) Desc   = 'Adjoint variables for species'
          IF ( isUnits ) Units  = 'varies'
          IF ( isRank  ) Rank   = 3
          IF ( isSpc   ) PerSpc = 'ALL'
       CASE ( 'COSTFUNCMASK' )
          IF ( isDesc    ) Desc  = 'Cost function volume mask'
          IF ( isUnits   ) Units = 'none'
          IF ( isRank    ) Rank  = 3
#endif

       CASE( 'BOUNDARYCOND' )
          IF ( isDesc  ) Desc   = 'Transport boundary conditions for species'
          IF ( isUnits ) Units  = 'v/v'
          IF ( isRank  ) Rank   = 3
          IF ( isSpc   ) PerSpc = 'ADV'

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
          IF ( isDesc  ) Desc  = &
               'Dry aerosol radius for mineral dust (0.15 um)'
          IF ( isUnits ) Units = 'cm'
          IF ( isRank  ) Rank  = 3

       CASE ( 'AERORADIMDUST2' )
          IF ( isDesc  ) Desc  = &
               'Dry aerosol radius for mineral dust (0.25 um)'
          IF ( isUnits ) Units = 'cm'
          IF ( isRank  ) Rank  = 3

       CASE ( 'AERORADIMDUST3' )
          IF ( isDesc  ) Desc  = &
               'Dry aerosol radius for mineral dust (0.4 um)'
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

!       CASE ( 'WETAERORADINITS' )
!          IF ( isDesc  ) Desc  = 'Wet aerosol radius for inorganic nitrates on' &
!                                  // 'surface of seasalt aerosol'
!          IF ( isUnits ) Units = 'cm'
!          IF ( isRank  ) Rank  = 3
!
!       CASE ( 'WETAERORADISALACL' )
!          IF ( isDesc  ) Desc  = 'Wet aerosol radius for chloride in Accumulation' &
!                                  // 'mode seasalt aerosol'
!          IF ( isUnits ) Units = 'cm'
!          IF ( isRank  ) Rank  = 3
!
!       CASE ( 'WETAERORADISALCCL' )
!          IF ( isDesc  ) Desc  = 'Wet aerosol radius for chloride in coarse' &
!                                  // 'mode seasalt aerosol'
!          IF ( isUnits ) Units = 'cm'
!          IF ( isRank  ) Rank  = 3
!
!       CASE ( 'WETAERORADISO4S' )
!          IF ( isDesc  ) Desc  = 'Wet aerosol radius for sulfate  on' &
!                                  // 'surface of seasalt aerosol'
!          IF ( isUnits ) Units = 'cm'
!          IF ( isRank  ) Rank  = 3

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

       CASE ( 'AEROH2OSNA' )
          IF ( isDesc  ) Desc  = 'Sulfur-nitrogen-ammonia aerosol water content'
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

       CASE ( 'SOILDUST1' )
          IF ( isDesc  ) Desc  = 'Dust aerosol concentration in bin 1'
          IF ( isUnits ) Units = 'kg/m3'
          IF ( isRank  ) Rank  =  3

       CASE ( 'SOILDUST2' )
          IF ( isDesc  ) Desc  = 'Dust aerosol concentration in bin 2'
          IF ( isUnits ) Units = 'kg/m3'
          IF ( isRank  ) Rank  =  3

       CASE ( 'SOILDUST3' )
          IF ( isDesc  ) Desc  = 'Dust aerosol concentration in bin 3'
          IF ( isUnits ) Units = 'kg/m3'
          IF ( isRank  ) Rank  =  3

       CASE ( 'SOILDUST4' )
          IF ( isDesc  ) Desc  = 'Dust aerosol concentration in bin 4'
          IF ( isUnits ) Units = 'kg/m3'
          IF ( isRank  ) Rank  =  3

       CASE ( 'SOILDUST5' )
          IF ( isDesc  ) Desc  = 'Dust aerosol concentration in bin 5'
          IF ( isUnits ) Units = 'kg/m3'
          IF ( isRank  ) Rank  =  3

       CASE ( 'SOILDUST6' )
          IF ( isDesc  ) Desc  = 'Dust aerosol concentration in bin 6'
          IF ( isUnits ) Units = 'kg/m3'
          IF ( isRank  ) Rank  =  3

       CASE ( 'SOILDUST7' )
          IF ( isDesc  ) Desc  = 'Dust aerosol concentration in bin 7'
          IF ( isUnits ) Units = 'kg/m3'
          IF ( isRank  ) Rank  =  3

       CASE ( 'GAMMAN2O5OVERALL' )
          IF ( isDesc  ) Desc = 'Sticking coefficient for Gamma N2O5 overall'
          IF ( isUnits ) Units = 'l'
          IF ( isRank  ) Rank = 3

       CASE ( 'GAMMAN2O5FINE' )
          IF ( isDesc  ) Desc = 'Sticking coefficient for Gamma N2O5 and fine aerosol'
          IF ( isUnits ) Units = 'l'
          IF ( isRank  ) Rank = 3

       CASE ( 'YIELDCLNO2FINE' )
          IF ( isDesc  ) Desc = 'Production yield coefficient for ClNO2 ' &
                               // ' from N2O5 fine aerosol uptake'
          IF ( isUnits ) Units = 'l'
          IF ( isRank  ) Rank = 3

       CASE ( 'KPPHVALUE' )
          IF ( isDesc  ) Desc  = 'H-value for Rosenbrock solver'
          IF ( isUnits ) Units = '1'
          IF ( isRank  ) Rank  = 3

       CASE ( 'OMOC' )
          IF ( isDesc  ) Desc  = 'OM:OC ratio as read by HEMCO (from /aerosol_mod.F90)'
          IF ( isUnits ) Units = '1'
          IF ( isRank  ) Rank  = 2

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

       CASE( 'ISORROPAEROPHACCUM' )
          IF ( isDesc  ) Desc  = 'ISORROPIA aerosol pH, accumulation mode'
          IF ( isUnits ) Units = '1'
          IF ( isRank  ) Rank  = 3

       CASE( 'ISORROPAEROPHCOARSE' )
          IF ( isDesc  ) Desc  = 'ISORROPIA aerosol pH, accumulation mode'
          IF ( isUnits ) Units = '1'
          IF ( isRank  ) Rank  = 3

       CASE( 'ISORROPHPLUSACCUM' )
          IF ( isDesc  ) Desc  = &
             'ISORROPIA H+ concentration, accumulation mode'
          IF ( isUnits ) Units = 'mol L-1'
          IF ( isRank  ) Rank  = 3

       CASE( 'ISORROPHPLUSCOARSE' )
          IF ( isDesc  ) Desc  = 'ISORROPIA H+ concentration, coarse mode'
          IF ( isUnits ) Units = 'mol L-1'
          IF ( isRank  ) Rank  = 3

       CASE( 'ISORROPAEROH2OACCUM' )
          IF ( isDesc  ) Desc  = &
             'ISORROPIA aerosol water concentration, accumulation mode'
          IF ( isUnits ) Units = 'ug m-3'
          IF ( isRank  ) Rank  = 3

       CASE( 'ISORROPAEROH2OCOARSE' )
          IF ( isDesc  ) Desc  = &
             'ISORROPIA aerosol water concentration, coarse mode'
          IF ( isUnits ) Units = 'ug m-3'
          IF ( isRank  ) Rank  = 3

       CASE( 'ISORROPSULFATE' )
          IF ( isDesc  ) Desc  = 'ISORROPIA sulfate concentration'
          IF ( isUnits ) Units = 'mol L-1'
          IF ( isRank  ) Rank  = 3

       CASE( 'ISORROPNITRATEACCUM' )
          IF ( isDesc  ) Desc  = &
             'ISORROPIA nitrate concentration, accumulation mode'
          IF ( isUnits ) Units = 'mol L-1'
          IF ( isRank  ) Rank  = 3

       CASE( 'ISORROPNITRATECOARSE' )
          IF ( isDesc  ) Desc  = &
             'ISORROPIA nitrate concentration, coarse mode'
          IF ( isUnits ) Units = 'mol L-1'
          IF ( isRank  ) Rank  = 3

       CASE( 'ISORROPCHLORIDEACCUM' )
          IF ( isDesc  ) Desc  = &
             'ISORROPIA chloride concentration, accumulation mode'
          IF ( isUnits ) Units = 'mol/L'
          IF ( isRank  ) Rank  = 3

       CASE( 'ISORROPCHLORIDECOARSE' )
          IF ( isDesc  ) Desc  = &
             'ISORROPIA chloride concentration, coarse mode'
          IF ( isUnits ) Units = 'mol/L'
          IF ( isRank  ) Rank  = 3

       CASE( 'ISORROPBISULFATE' )
          IF ( isDesc  ) Desc  = 'ISORROPIA Bisulfate (general acid)' &
                                 // ' concentration'
          IF ( isUnits ) Units = 'mol L-1'
          IF ( isRank  ) Rank  =  3

       CASE( 'PHCLOUD' )
          IF ( isDesc  ) Desc  = 'Cloud pH'
          IF ( isUnits ) Units = '1'
          IF ( isRank  ) Rank  =  3

       CASE( 'QLXPHCLOUD' )
          IF ( isDesc  ) Desc  = 'Cloud pH * Met_QL'
          IF ( isUnits ) Units = '1'
          IF ( isRank  ) Rank  =  3

       CASE( 'ORVCSESQ' )
          IF ( isDesc  ) Desc  = 'Sesquiterpenes mass'
          IF ( isUnits ) Units = 'kg'
          IF ( isRank  ) Rank  =  3

       CASE( 'ISCLOUD' )
          IF ( isDesc  ) Desc  = 'Cloud presence'
          IF ( isUnits ) Units = '1'
          IF ( isRank  ) Rank  =  3

       CASE( 'SSALKACCUMMODE' )
          IF ( isDesc  ) Desc  = 'Sea salt alkalinity, accumulation mode'
          IF ( isUnits ) Units = '1'
          IF ( isRank  ) Rank  =  3

       CASE( 'SSALKCOARSEMODE' )
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

       CASE ( 'FUPDATEHOCL' )
          IF ( isDesc  ) Desc  = 'Correction factor for HOCl removal by SO2'
          IF ( isUnits ) Units = '1'
          IF ( isRank  ) Rank  =  3

       CASE ( 'ACLAREA' )
          IF ( isDesc  ) Desc  = 'Dry aerosol area for fine mode Cl-'
          IF ( isUnits ) Units = 'cm2 cm-3'
          IF ( isRank  ) Rank  =  3

       CASE ( 'ACLRADI' )
          IF ( isDesc  ) Desc  = 'Dry aerosol radius for fine mode Cl-'
          IF ( isUnits ) Units = 'cm'
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
          IF ( isDesc  ) Desc   = 'Hg(0) ocean mass'
          IF ( isUnits ) Units  = 'kg'
          IF ( isRank  ) Rank   = 2

       CASE( 'OCEANHG2' )
          IF ( isDesc  ) Desc   = 'Hg(II) ocean mass'
          IF ( isUnits ) Units  = 'kg'
          IF ( isRank  ) Rank   = 2

       CASE( 'OCEANHGP' )
          IF ( isDesc  ) Desc   = 'HgP ocean mass'
          IF ( isUnits ) Units  = 'kg'
          IF ( isRank  ) Rank   = 2

       CASE( 'SNOWHGOCEAN' )
          IF ( isDesc  ) Desc   = 'Reducible Hg snowpack on ocean'
          IF ( isUnits ) Units  = 'kg'
          IF ( isRank  ) Rank   = 2

       CASE( 'SNOWHGLAND' )
          IF ( isDesc  ) Desc   = 'Reducible Hg snowpack on land'
          IF ( isUnits ) Units  = 'kg'
          IF ( isRank  ) Rank   = 2

       CASE( 'SNOWHGOCEANSTORED' )
          IF ( isDesc  ) Desc   = 'Non-reducible Hg snowpack on ocean'
          IF ( isUnits ) Units  = 'kg'
          IF ( isRank  ) Rank   = 2

       CASE( 'SNOWHGLANDSTORED' )
          IF ( isDesc  ) Desc   = 'Non-reducible Hg snowpack on land'
          IF ( isUnits ) Units  = 'kg'
          IF ( isRank  ) Rank   = 2

       CASE ( 'IODIDE' )
          IF ( isDesc  ) Desc  = 'Surface iodide concentration'
          IF ( isUnits ) Units = 'nM'
          IF ( isRank  ) Rank  = 2

       CASE ( 'SALINITY' )
          IF ( isDesc  ) Desc  = 'Salinity'
          IF ( isUnits ) Units = 'PSU'
          IF ( isRank  ) Rank  = 2

       CASE( 'DRYDEPFREQ' )
          IF ( isDesc  ) Desc   = 'Dry deposition frequencies'
          IF ( isUnits ) Units  = 's-1'
          IF ( isRank  ) Rank   = 2
          IF ( isSpc   ) perSpc = 'DRY'

       CASE( 'DRYDEPVEL' )
          IF ( isDesc    ) Desc   = 'Dry deposition velocities'
          IF ( isUnits   ) Units  = 'm s-1'
          IF ( isRank    ) Rank   = 2
          IF ( isSpc     ) perSpc = 'DRY'

#ifdef MODEL_GEOS
       CASE( 'DRYDEPRA2M' )
          IF ( isDesc    ) Desc  = '2 meter aerodynamic resistance'
          IF ( isUnits   ) Units = 's cm-1'
          IF ( isRank    ) Rank  = 2

       CASE( 'DRYDEPRA10M' )
          IF ( isDesc    ) Desc  = '10 meter aerodynamic resistance'
          IF ( isUnits   ) Units = 's cm-1'
          IF ( isRank    ) Rank  = 2

#endif
       CASE( 'JOH' )
          IF ( isDesc    ) Desc  = 'Surface J-values for reaction O3 + hv --> O2 + O'
          IF ( isUnits   ) Units = '1'
          IF ( isRank    ) Rank  = 2

       CASE( 'JNO2' )
          IF ( isDesc    ) Desc  = 'Surface J-values for reaction NO2 + hv --> NO + O'
          IF ( isUnits   ) Units = '1'
          IF ( isRank    ) Rank  = 2

       CASE( 'SURFACEFLUX' )
          IF ( isDesc  ) Desc   = 'Surface flux (E-D) for non-local PBL mixing'
          IF ( isUnits ) Units  = 'kg m-2 s-1'
          IF ( isRank  ) Rank   = 2
          IF ( isSpc   ) perSpc = 'ADV'

       CASE( 'TLSTT' )
          IF ( isDesc  ) Desc  = 'TLSTT'
          IF ( isUnits ) Units = ''
          IF ( isRank  ) Rank  = 4

       CASE( 'CH4_EMIS' )
          IF ( isDesc  ) Desc  = 'CH4 emissions by sector, CH4 specialty simulation only'
          IF ( isUnits ) Units = 'kg/m2/s'
          IF ( isRank  ) Rank  = 3

       CASE( 'BOH' )
          IF ( isDesc  ) Desc  = 'OH values, CH4 specialty simulation only'
          IF ( isUnits ) Units = 'molec/cm3'
          IF ( isRank  ) Rank  = 3

       CASE( 'BCL' )
          IF ( isDesc  ) Desc  = 'Cl values, CH4 specialty simulation only'
          IF ( isUnits ) Units = 'v/v'
          IF ( isRank  ) Rank  = 3

       CASE( 'QQ3D' )
          IF ( isDesc  ) Desc  = 'Rate of new precipitation formation'
          IF ( isUnits ) Units = 'cm3 H2O cm-3 air'
          IF ( isRank  ) Rank  = 3

       CASE( 'KRATE' )
          IF ( isDesc  ) Desc  = 'KRATE'
          IF ( isUnits ) Units = '1'
          IF ( isRank  ) Rank  = 3

       CASE( 'PHRAIN' )
          IF ( isDesc  ) Desc  = 'Rain pH'
          IF ( isUnits ) Units = '1'
          IF ( isRank  ) Rank  = 3

       CASE( 'QQPHRAIN' )
          IF ( isDesc  ) Desc  = 'QQRain pH'
          IF ( isUnits ) Units = '1'
          IF ( isRank  ) Rank  = 3

       CASE( 'QQRAIN' )
          IF ( isDesc  ) Desc  = 'QQRain'
          IF ( isUnits ) Units = '1'
          IF ( isRank  ) Rank  = 3

#if defined(MODEL_CESM)
       CASE( 'H2SO4_PRDR' )
          IF ( isDesc  ) Desc  = 'H2SO4 production rate in timestep'
          IF ( isUnits ) Units = 'mol mol-1'
          IF ( isRank  ) Rank  = 3
#endif

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
! !IROUTINE: Init_and_Register_R4_2D
!
! !DESCRIPTION: Allocates the data array for a State_Chm field,
!  and also adds the field into the State_Chm registry.
!  This particular routine is for 4-byte, 2-dimensional array fields.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Init_and_Register_R4_2D( Input_Opt, State_Chm, State_Grid,      &
                                      Ptr2Data,  chmId,     RC,              &
                                      noRegister                            )
!
! !USES:
!
    USE Input_Opt_Mod,  ONLY : OptInput
    USE State_Grid_Mod, ONLY : GrdState
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput),   INTENT(IN)  :: Input_Opt           ! Input Options
    TYPE(ChmState),   INTENT(IN)  :: State_Chm           ! Chemistry State
    TYPE(GrdState),   INTENT(IN)  :: State_Grid          ! Grid State
    CHARACTER(LEN=*), INTENT(IN)  :: chmId               ! Field name
    LOGICAL,          OPTIONAL    :: noRegister          ! Exit after init
!
! !INPUT/OUTPUT PARAMETERS:
!
    REAL(f4),         POINTER     :: Ptr2Data(:,:)       ! Pointer to data
!
! !OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(OUT) :: RC                  ! Success or failure?
!
! !REVISION HISTORY:
!  21 Sep 2020 - R. Yantosca - Initial version
!  See the subsequent Git history with the gitk browser!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER            :: NX, NY
    LOGICAL            :: doRegister

    ! Strings
    CHARACTER(LEN=255) :: arrayId

    !========================================================================
    ! Init_and_Register_R4_2D begins here!
    !========================================================================

    ! Initialize
    RC      = GC_SUCCESS
    arrayId = 'State_Chm%' // TRIM( chmId )

    IF ( PRESENT( noRegister ) ) THEN
       doRegister = ( .not. noRegister )
    ELSE
       doRegister = .TRUE.
    ENDIF

    !========================================================================
    ! Allocate the field array (if it hasn't already been allocated)
    !========================================================================
    IF ( .not. ASSOCIATED( Ptr2Data ) ) THEN

       ! Get array dimensions
       NX = State_Grid%NX
       NY = State_Grid%NY

       ! Allocate the array
       ALLOCATE( Ptr2Data( NX, NY ), STAT=RC )
       CALL GC_CheckVar( arrayId, 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       Ptr2Data = 0.0_f4

    ENDIF

    !========================================================================
    ! Register the field (unless we explicitly say not to)
    !========================================================================
    IF ( doRegister ) THEN
       CALL Register_ChmField( Input_Opt, chmId, Ptr2Data, State_Chm, RC )
       CALL GC_CheckVar( arrayId, 1, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

  END SUBROUTINE Init_and_Register_R4_2D
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Init_and_Register_R4_3D
!
! !DESCRIPTION: Allocates the data array for a State_Chm field,
!  and also adds the field into the State_Chm registry.
!  This particular routine is for 4-byte, 3-dimensional array fields.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Init_and_Register_R4_3D( Input_Opt, State_Chm, State_Grid,      &
                                      Ptr2Data,  chmId,     RC,              &
                                      nSlots,    nCat,      noRegister      )
!
! !USES:
!
    USE Input_Opt_Mod,  ONLY : OptInput
    USE State_Grid_Mod, ONLY : GrdState
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput),   INTENT(IN)  :: Input_Opt           ! Input Options
    TYPE(ChmState),   INTENT(IN)  :: State_Chm           ! Chemistry State
    TYPE(GrdState),   INTENT(IN)  :: State_Grid          ! Grid State
    CHARACTER(LEN=*), INTENT(IN)  :: chmId               ! Field name
    INTEGER,          OPTIONAL    :: nSlots              ! # slots, 3rd dim
    INTEGER,          OPTIONAL    :: nCat                ! Category index
    LOGICAL,          OPTIONAL    :: noRegister          ! Exit after init
!
! !INPUT/OUTPUT PARAMETERS:
!
    REAL(f4),         POINTER     :: Ptr2Data(:,:,:)     ! Pointer to data
!
! !OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(OUT) :: RC                  ! Success or failure?
!
! !REVISION HISTORY:
!  21 Sep 2020 - R. Yantosca - Initial version
!  See the subsequent Git history with the gitk browser!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER            :: NX, NY, NZ, NW
    LOGICAL            :: doRegister

    ! Strings
    CHARACTER(LEN=255) :: arrayId

    !========================================================================
    ! Init_and_Register_R4_3D begins here!
    !========================================================================

    ! Initialize
    RC      = GC_SUCCESS
    arrayId = 'State_Chm%' // TRIM( chmId )

    IF ( PRESENT( noRegister ) ) THEN
       doRegister = ( .not. noRegister )
    ELSE
       doRegister = .TRUE.
    ENDIF

    !========================================================================
    ! Allocate the field array (if it hasn't already been allocated)
    !========================================================================
    IF ( .not. ASSOCIATED( Ptr2Data ) ) THEN

       ! Get array ID and dimensions
       ! If optional nSlots is passed, use it for the 3rd dimension
       NX = State_Grid%NX
       NY = State_Grid%NY
       IF ( PRESENT( nSlots ) ) THEN
          NW = nSlots
       ELSE
          NW = State_Grid%NZ
       ENDIF

       ! Allocate the array
       ALLOCATE( Ptr2Data( NX, NY, NW ), STAT=RC )
       CALL GC_CheckVar( arrayId, 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       Ptr2Data = 0.0_f4

    ENDIF

    !========================================================================
    ! Register the field (unless we explicitly say not to)
    !========================================================================
    IF ( doRegister ) THEN
       CALL Register_ChmField( Input_Opt, chmId, Ptr2Data,                   &
                               State_Chm, RC,    nCat=nCat                  )
       CALL GC_CheckVar( arrayId, 1, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

  END SUBROUTINE Init_and_Register_R4_3D
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Init_and_Register_R4_4D
!
! !DESCRIPTION: Allocates the data array for a State_Chm field,
!  and also adds the field into the State_Chm registry.
!  This particular routine is for 4-byte, 4-dimensional arrays.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Init_and_Register_R4_4D( Input_Opt, State_Chm, State_Grid,      &
                                      Ptr2Data,  chmId,     nSlots,          &
                                      RC,        nCat,      noRegister      )
!
! !USES:
!
    USE Input_Opt_Mod,  ONLY : OptInput
    USE State_Grid_Mod, ONLY : GrdState
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput),   INTENT(IN)  :: Input_Opt           ! Input Options
    TYPE(ChmState),   INTENT(IN)  :: State_Chm           ! Chemistry State
    TYPE(GrdState),   INTENT(IN)  :: State_Grid          ! Grid State
    CHARACTER(LEN=*), INTENT(IN)  :: chmId               ! Field name

    INTEGER,          INTENT(IN)  :: nSlots              ! # of slots, 4th dim
    INTEGER,          OPTIONAL    :: nCat                ! Optional category
    LOGICAL,          OPTIONAL    :: noRegister          ! Exit after init
!
! !INPUT/OUTPUT PARAMETERS:
!
    REAL(f4),         POINTER     :: Ptr2Data(:,:,:,:)   ! Pointer to data
!
! !OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(OUT) :: RC                  ! Success or failure?
!
! !REVISION HISTORY:
!  21 Sep 2020 - R. Yantosca - Initial version
!  See the subsequent Git history with the gitk browser!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER            :: NX, NY, NZ
    LOGICAL            :: doRegister

    ! Strings
    CHARACTER(LEN=255) :: arrayId

    !========================================================================
    ! Init_and_Register_R4_4D begins here!
    !========================================================================

    ! Initialize
    RC      = GC_SUCCESS
    arrayId = 'State_Chm%' // TRIM( chmId )

    IF ( PRESENT( noRegister ) ) THEN
       doRegister = ( .not. noRegister )
    ELSE
       doRegister = .TRUE.
    ENDIF

    !========================================================================
    ! Allocate the field array
    !========================================================================
    IF ( .not. ASSOCIATED( Ptr2Data ) ) THEN

       ! Get array dimensions
       NX = State_Grid%NX
       NY = State_Grid%NY
       NZ = State_Grid%NZ

       ALLOCATE( Ptr2Data( NX, NY, NZ, nSlots ), STAT=RC )
       CALL GC_CheckVar( arrayId, 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       Ptr2Data = 0.0_f4

    ENDIF

    !========================================================================
    ! Register the field
    !========================================================================
    IF ( doRegister ) THEN
       CALL Register_ChmField( Input_Opt, chmId, Ptr2Data,                   &
                               State_Chm, RC,    nCat=nCat                  )
       CALL GC_CheckVar( arrayId, 1, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

  END SUBROUTINE Init_and_Register_R4_4D
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Init_and_Register_R8_2D
!
! !DESCRIPTION: Allocates the data array for a State_Chm field,
!  and also adds the field into the State_Chm registry.
!  This particular routine is for 8-byte, 2-dimensional fields.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Init_and_Register_R8_2D( Input_Opt, State_Chm, State_Grid,      &
                                      Ptr2Data,  chmId,     RC,              &
                                      noRegister                            )
!
! !USES:
!
    USE Input_Opt_Mod,  ONLY : OptInput
    USE State_Grid_Mod, ONLY : GrdState
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput),   INTENT(IN)  :: Input_Opt           ! Input Options
    TYPE(ChmState),   INTENT(IN)  :: State_Chm           ! Chemistry State
    TYPE(GrdState),   INTENT(IN)  :: State_Grid          ! Grid State
    CHARACTER(LEN=*), INTENT(IN)  :: chmId               ! Field name
    LOGICAL,          OPTIONAL    :: noRegister          ! Exit after init
!
! !INPUT/OUTPUT PARAMETERS:
!
    REAL(f8),         POINTER     :: Ptr2Data(:,:)       ! Pointer to data
!
! !OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(OUT) :: RC                  ! Success or failure?
!
! !REVISION HISTORY:
!  21 Sep 2020 - R. Yantosca - Initial version
!  See the subsequent Git history with the gitk browser!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER            :: NX, NY
    LOGICAL            :: doRegister

    ! Strings
    CHARACTER(LEN=255) :: arrayId

    !=======================================================================
    ! Init_and_Register_R8_2D begins here!
    !=======================================================================

    ! Initialize
    RC      = GC_SUCCESS
    arrayId = 'State_Chm%' // TRIM( chmId )

    IF ( PRESENT( noRegister ) ) THEN
       doRegister = ( .not. noRegister )
    ELSE
       doRegister = .TRUE.
    ENDIF

    !=======================================================================
    ! Allocate the field array (if it hasn't already been allocated)
    !=======================================================================
    IF ( .not. ASSOCIATED( Ptr2Data ) ) THEN

       ! Get array dimensions
       NX = State_Grid%NX
       NY = State_Grid%NY

       ! Allocate the data
       ALLOCATE( Ptr2Data( NX, NY ), STAT=RC )
       CALL GC_CheckVar( arrayId, 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       Ptr2Data = 0.0_f8

    ENDIF

    !=======================================================================
    ! Register the field
    !=======================================================================
    IF ( doRegister ) THEN
       CALL Register_ChmField( Input_Opt, chmId, Ptr2Data, State_Chm, RC )
       CALL GC_CheckVar( arrayId, 1, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

  END SUBROUTINE Init_and_Register_R8_2D
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Init_and_Register_R8_3D
!
! !DESCRIPTION: Allocates the data array for a State_Chm field,
!  and also adds the field into the State_Chm registry.
!  This particular routine is for 4-byte, 2-dimensional arrays.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Init_and_Register_R8_3D( Input_Opt, State_Chm, State_Grid,      &
                                      Ptr2Data,  chmId,     RC,              &
                                      nSlots,    nCat,      noRegister      )
!
! !USES:
!
    USE Input_Opt_Mod,  ONLY : OptInput
    USE State_Grid_Mod, ONLY : GrdState
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput),   INTENT(IN)  :: Input_Opt           ! Input Options
    TYPE(ChmState),   INTENT(IN)  :: State_Chm           ! Chemistry State
    TYPE(GrdState),   INTENT(IN)  :: State_Grid          ! Grid State
    CHARACTER(LEN=*), INTENT(IN)  :: chmId               ! Field name
    INTEGER,          OPTIONAL    :: nSlots              ! # slots, 3rd dim
    INTEGER,          OPTIONAL    :: nCat                ! Category index
    LOGICAL,          OPTIONAL    :: noRegister          ! Exit after init
!
! !INPUT/OUTPUT PARAMETERS:
!
    REAL(f8),         POINTER     :: Ptr2Data(:,:,:)     ! Pointer to data
!
! !OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(OUT) :: RC                  ! Success or failure?
!
! !REVISION HISTORY:
!  21 Sep 2020 - R. Yantosca - Initial version
!  See the subsequent Git history with the gitk browser!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER            :: NX, NY, NZ, NW
    LOGICAL            :: doRegister

    ! Strings
    CHARACTER(LEN=255) :: arrayId

    !========================================================================
    ! Init_and_Register_R8_3D begins here!
    !========================================================================

    ! Initialize
    RC      = GC_SUCCESS
    arrayId = 'State_Chm%' // TRIM( chmId )

    IF ( PRESENT( noRegister ) ) THEN
       doRegister = ( .not. noRegister )
    ELSE
       doRegister = .TRUE.
    ENDIF

    !========================================================================
    ! Allocate the field array (if it hasn't already been allocated)
    !========================================================================
    IF ( .not. ASSOCIATED( Ptr2Data ) ) THEN

       ! Get array dimensions
       NX = State_Grid%NX
       NY = State_Grid%NY
       IF ( PRESENT( nSlots ) ) THEN
          NW = nSlots
       ELSE
          NW = State_Grid%NZ
       ENDIF

       ! Allocate the array
       ALLOCATE( Ptr2Data( NX, NY, NW ), STAT=RC )
       CALL GC_CheckVar( arrayId, 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       Ptr2Data = 0.0_f8

    ENDIF

    !========================================================================
    ! Register the field (unless we explicitly say not to)
    !========================================================================
    IF ( doRegister ) THEN
       CALL Register_ChmField( Input_Opt, chmId, Ptr2Data,                   &
                               State_Chm, RC,    nCat=nCat                  )
       CALL GC_CheckVar( arrayId, 1, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

  END SUBROUTINE Init_and_Register_R8_3D
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Init_and_Register_R8_4D
!
! !DESCRIPTION: Allocates the data array for a State_Chm field,
!  and also adds the field into the State_Chm registry.
!  This particular routine is for 4-byte, 2-dimensional arrays.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Init_and_Register_R8_4D( Input_Opt, State_Chm, State_Grid,      &
                                      Ptr2Data,  chmId,     nSlots,          &
                                      RC,        nCat,      noRegister      )
!
! !USES:
!
    USE Input_Opt_Mod,  ONLY : OptInput
    USE State_Grid_Mod, ONLY : GrdState
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput),   INTENT(IN)  :: Input_Opt           ! Input Options
    TYPE(ChmState),   INTENT(IN)  :: State_Chm           ! Chemistry State
    TYPE(GrdState),   INTENT(IN)  :: State_Grid          ! Grid State
    CHARACTER(LEN=*), INTENT(IN)  :: chmId               ! Field name
    INTEGER,          INTENT(IN)  :: nSlots              ! # of slots, 4th dim
    INTEGER,          OPTIONAL    :: nCat                ! Optional category
    LOGICAL,          OPTIONAL    :: noRegister          ! Exit after init
!
! !INPUT/OUTPUT PARAMETERS:
!
    REAL(f8),         POINTER     :: Ptr2Data(:,:,:,:)   ! Pointer to data
!
! !OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(OUT) :: RC                  ! Success/failure!
!
! !REVISION HISTORY:
!  21 Sep 2020 - R. Yantosca - Initial version
!  See the subsequent Git history with the gitk browser!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER            :: NX, NY, NZ
    LOGICAL            :: doRegister

    ! Strings
    CHARACTER(LEN=255) :: arrayId

    !========================================================================
    ! Init_and_Register_R8_4D begins here!
    !========================================================================

    ! Initialize
    RC      = GC_SUCCESS
    arrayId = 'State_Chm%' // TRIM( chmId )

    IF ( PRESENT( noRegister ) ) THEN
       doRegister = ( .not. noRegister )
    ELSE
       doRegister = .TRUE.
    ENDIF

    !========================================================================
    ! Allocate the field array (if it hasn't already been allocated)
    !========================================================================
    IF ( .not. ASSOCIATED( Ptr2Data ) ) THEN

       ! Get array dimensions
       NX = State_Grid%NX
       NY = State_Grid%NY
       NZ = State_Grid%NZ

       ! Allocate the array
       ALLOCATE( Ptr2Data( NX, NY, NZ, nSlots ), STAT=RC )
       CALL GC_CheckVar( arrayId, 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       Ptr2Data = 0.0_f8

    ENDIF

    !========================================================================
    ! Register the field (unless we explicitly say not to)
    !========================================================================
    IF ( doRegister ) THEN
       CALL Register_ChmField( Input_Opt, chmId, Ptr2Data,                   &
                               State_Chm, RC,    nCat=nCat                  )
       CALL GC_CheckVar( arrayId, 1, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

  END SUBROUTINE Init_and_Register_R8_4D
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Test_for_Species_Dim
!
! !DESCRIPTION: Returns true if a State_Chm quantity has a species dimension.
!\\
!\\
! !INTERFACE:
!
  FUNCTION Test_for_Species_Dim( perSpc ) RESULT( returnCode )
!
! !INPUT PARAMETERS:
!
    CHARACTER(LEN=*), INTENT(IN) :: perSpc      ! PerSpc value from metadata
!
! !RETURN VALUE:
!
    INTEGER                      :: returnCode  !  1  = has species dimension
                                                !  0  = no species dimension
                                                ! -1  = unknown perSpc value
!
! !REVISION HISTORY:
!  06 Jan 2015 - R. Yantosca - Initial version
!  See the subsequent Git history with the gitk browser!
!EOP
!------------------------------------------------------------------------------
!BOC
    SELECT CASE( TRIM( perSpc ) )
       CASE( 'ADV', 'ALL', 'DRY', 'WET' )
          returnCode = 1
       CASE( '' )
          returnCode = 0
       CASE DEFAULT
          returnCode = -1
    END SELECT

  END FUNCTION Test_For_Species_Dim
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get_NumSlots
!
! !DESCRIPTION: Returns the number of slots with which to define a
!  species-based array of State_Chm.
!\\
!\\
! !INTERFACE:
!
  FUNCTION Get_NumSlots( perSpc, State_Chm ) RESULT( nSlots )
!
! !INPUT PARAMETERS:
!
    CHARACTER(LEN=*), INTENT(IN) :: perSpc      ! PerSpc value from metadata
    TYPE(ChmState),   INTENT(IN) :: State_Chm   ! Chemistry State object
!
! !RETURN VALUE:
!
    INTEGER                      :: nSlots      ! Number of slots

! !REVISION HISTORY:
!  23 Sep 2020 - R. Yantosca - Initial version
!  See the subsequent Git history with the gitk browser!
!EOP
!------------------------------------------------------------------------------
!BOC

    SELECT CASE( TRIM( perSpc ) )
       CASE( 'ADV'   )
          nSlots = State_Chm%nAdvect
       CASE( 'ALL'   )
          nSlots = State_Chm%nSpecies
       CASE( 'DRY'   )
          nSlots = State_Chm%nDryDep
       CASE( 'WET'   )
          nSlots = State_Chm%nWetDep
       CASE DEFAULT
          nSlots = -1
    END SELECT

  END FUNCTION Get_NumSlots
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get_Diagnostic_Name
!
! !DESCRIPTION: Returns the diagnostic name and description of a species-based
!  quantity (appending the species to the base name if necessary).
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Get_Diagnostic_Name( State_Chm, perSpc,    N,        name,      &
                                  desc,      diagName,  diagDesc            )
!
! !INPUT PARAMETERS:
!
    TYPE(ChmState),     INTENT(IN)  :: State_Chm  ! Chemistry State
    CHARACTER(LEN=*),   INTENT(IN)  :: perSpc     ! PerSpc value from metadata
    INTEGER,            INTENT(IN)  :: N          ! Diagnostic index
    CHARACTER(LEN=*),   INTENT(IN)  :: name       ! Name from metadata
    CHARACTER(LEN=*),   INTENT(IN)  :: desc       ! Description from metadata
!
! !OUTPUT PARAMETERS:
!
    CHARACTER(LEN=255), INTENT(OUT) :: diagName   ! Name        + species name
    CHARACTER(LEN=255), INTENT(OUT) :: diagDesc   ! Description + species name
!
! !REVISION HISTORY:
!  20 Jan 2015 - R. Yantosca - Initial version
!  See the subsequent Git history with the gitk browser!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER                :: modelId

    ! Objects
    TYPE(Species), POINTER :: ThisSpc

    !---------------------------------------------------------------------
    ! All other species-bound quantities
    !---------------------------------------------------------------------

    ! Get the species index from the diagnostic index
    ! depending on the value of PerSpc (bmy, 05 Oct 2021)
    modelId = N
    IF ( PerSpc == 'DRY' ) modelId = State_Chm%Map_DryDep(N)
    IF ( PerSpc == 'WET' ) modelId = State_Chm%Map_WetDep(N)
    
    ! Point to the proper species, by modelId
    ThisSpc => State_Chm%SpcData(modelId)%Info

    ! Append the species name to the diagnostic name with an underscore
    diagName = TRIM( name ) // '_' // TRIM( ThisSpc%Name )

    ! Append the species name to the diagnostic description
    diagDesc = TRIM( desc ) // ' ' // TRIM( ThisSpc%Name )

    ! Free pointer
    ThisSpc => NULL()

  END SUBROUTINE Get_Diagnostic_Name
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Register_ChmField_R4_2D
!
! !DESCRIPTION: Registers a 2-dimensional, 4-byte real array field
!  of the State\_Chm object.  This allows the diagnostic modules get
!  a pointer to the field by searching on the field name.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Register_ChmField_R4_2D( Input_Opt, metadataID, Ptr2Data,       &
                                      State_Chm, RC                         )
!
! !USES:
!
    USE Input_Opt_Mod,      ONLY : OptInput
    USE Registry_Params_Mod
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput),   INTENT(IN)  :: Input_Opt         ! Input Options
    CHARACTER(LEN=*), INTENT(IN)  :: metadataID        ! State_Chm field ID
    REAL(f4),         POINTER     :: Ptr2Data(:,:)     ! Pointer to data
    TYPE(ChmState),   INTENT(IN)  :: State_Chm         ! Chemistry State
!
! !OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(OUT) :: RC                ! Success or failure?
!
! !REMARKS:
!
! !REVISION HISTORY:
!  21 Sep 2020 - E. Lundgren - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    LOGICAL            :: found,       onEdges
    INTEGER            :: N,           rank,        hasSpeciesDim
    INTEGER            :: type,        vloc

    ! Strings
    CHARACTER(LEN=255) :: errMsg_reg,  thisLoc,     desc
    CHARACTER(LEN=255) :: thisSpcName, thisSpcDesc, perSpc
    CHARACTER(LEN=255) :: units
    CHARACTER(LEN=512) :: errMsg

    !========================================================================
    ! Initialize
    !========================================================================
    RC         = GC_SUCCESS
    errMsg     = ''
    errMsg_reg = 'Error encountered while registering State_Chm%'
    thisLoc    = &
       ' -> at Register_ChmField_R4_2D (in Headers/state_chm_mod.F90)'

    !========================================================================
    ! Get metadata
    !========================================================================
    CALL Get_MetaData_State_Chm(                                             &
         am_I_Root  = Input_Opt%amIRoot,                                     &
         metadataId = metadataId,                                            &
         found      = found,                                                 &
         desc       = desc,                                                  &
         units      = units,                                                 &
         rank       = rank,                                                  &
         type       = type,                                                  &
         vloc       = vloc,                                                  &
         perSpc     = perSpc,                                                &
         RC         = RC                                                    )

    ! Trap potential errors
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = TRIM( ErrMsg_reg ) // TRIM( MetadataID ) //                  &
                '; Abnormal exit from routine "Get_Metadata_State_Chm"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    ! Is the data placed on vertical edges?
    onEdges = ( vLoc == VLocationEdge )

    ! Test if the data has a species dimension
    hasSpeciesDim = Test_For_Species_Dim( perSpc )

    !------------------------------------------------------------------------
    ! If not tied to species then simply register the single field
    !------------------------------------------------------------------------
    IF ( hasSpeciesDim == 0 ) THEN

       ! Check that metadata consistent with data size
       IF ( rank /= 2 ) THEN
          ErrMsg = 'Data dims and metadata rank do not match for '           &
                   // TRIM(metadataID)
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

       ! Add field to registry
       CALL Registry_AddField(                                               &
            Input_Opt    = Input_Opt,                                        &
            Registry     = State_Chm%Registry,                               &
            State        = State_Chm%State,                                  &
            Variable     = TRIM( metadataID ),                               &
            Description  = TRIM( desc       ),                               &
            Units        = TRIM( units      ),                               &
            Data2d_4     = Ptr2Data,                                         &
            RC           = RC                                               )

       ! Trap potential errors
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = TRIM( ErrMsg_reg ) // TRIM( MetadataID ) //               &
                  '; Abnormal exit from routine "Registry_AddField"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

    !------------------------------------------------------------------------
    ! Otherwise exit with error
    !------------------------------------------------------------------------
    ELSE

       ! Error msg
       ErrMsg = 'Handling of PerSpc metadata ' // TRIM(perSpc) //    &
                ' is not implemented for this combo of data type and size'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN

    ENDIF

  END SUBROUTINE Register_ChmField_R4_2D
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
  SUBROUTINE Register_ChmField_R4_3D( Input_Opt,  metadataID, Ptr2Data,      &
                                      State_Chm,  RC,         nCat          )
!
! !USES:
!
    USE Input_Opt_Mod,      ONLY : OptInput
    USE Registry_Params_Mod
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput),   INTENT(IN)  :: Input_Opt        ! Input Options
    CHARACTER(LEN=*), INTENT(IN)  :: metadataID       ! State_Chm field ID
    REAL(f4),         POINTER     :: Ptr2Data(:,:,:)  ! pointer to data
    TYPE(ChmState),   INTENT(IN)  :: State_Chm        ! Chemistry State
    INTEGER,          OPTIONAL    :: nCat             ! Category index
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(OUT) :: RC               ! Success or failure?
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
    ! Scalars
    LOGICAL            :: found,       onEdges
    INTEGER            :: N,           rank,        type
    INTEGER            :: vloc,        nSlots,      hasSpeciesDim

    ! Strings
    CHARACTER(LEN=255) :: errMsg_reg,  thisLoc,     desc
    CHARACTER(LEN=255) :: thisSpcName, thisSpcDesc, perSpc
    CHARACTER(LEN=255) :: units
    CHARACTER(LEN=512) :: errMsg

    !=========================================================================
    ! Initialize
    !=========================================================================
    RC         = GC_SUCCESS
    errMsg     = ''
    errMsg_reg = 'Error encountered while registering State_Chm%'
    thisLoc    = &
       ' -> at Register_ChmField_R4_3D (in Headers/state_chm_mod.F90)'

    !=========================================================================
    ! Get metadata
    !=========================================================================
    CALL Get_MetaData_State_Chm(                                              &
         am_I_Root  = Input_Opt%amIRoot,                                      &
         metadataId = metadataId,                                             &
         found      = found,                                                  &
         desc       = desc,                                                   &
         units      = units,                                                  &
         rank       = rank,                                                   &
         type       = type,                                                   &
         vloc       = vloc,                                                   &
         perSpc     = perSpc,                                                 &
         RC         = RC                                                     )

    ! Trap potential errors
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = TRIM( ErrMsg_reg ) // TRIM( MetadataID ) //                  &
                '; Abnormal exit from routine "Get_Metadata_State_Chm"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    ! Test if the quantity has a species dimension
    hasSpeciesDim = Test_for_Species_Dim( perSpc )

    !------------------------------------------------------------------------
    ! If tied to a given category, only register that one
    !------------------------------------------------------------------------
    IF ( PRESENT( nCat ) ) THEN

       ! Check that metadata consistent with data size
       IF ( rank /= 2 ) THEN
          ErrMsg = 'Data dims and metadata rank do not match for '           &
                   // TRIM(metadataID)
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

       ! Add field to registry
       CALL Registry_AddField(                                               &
            Input_Opt    = Input_Opt,                                        &
            Registry     = State_Chm%Registry,                               &
            State        = State_Chm%State,                                  &
            Variable     = TRIM( metadataID ),                               &
            Description  = TRIM( desc       ),                               &
            Units        = TRIM( units      ),                               &
            OnLevelEdges = onEdges,                                          &
            Data2d_4     = Ptr2Data(:,:,nCat),                               &
            RC           = RC                                               )

       ! Trap potential errors
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = TRIM( ErrMsg_reg ) // TRIM( MetadataID ) //               &
                   '; Abnormal exit from Registry_AddField!'
          CALL GC_Error( ErrMsg_reg, RC, ThisLoc )
          RETURN
      ENDIF

    !------------------------------------------------------------------------
    ! If tied to a particular species, register each species individually
    !------------------------------------------------------------------------
    ELSE IF ( hasSpeciesDim == 1 ) THEN

       ! Check that metadata consistent with data size
       IF ( rank /= 2 ) THEN
          ErrMsg = 'Data dims and metadata rank do not match for '           &
                   // TRIM(metadataID)
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

       ! Get the length of the species-based dimension
       nSlots = Get_NumSlots( perSpc, State_Chm )

       ! Loop over all species
       DO N = 1, nSlots

          ! Append the species name to the diagnostic name & description
          CALL Get_Diagnostic_Name(                                          &
               State_Chm    = State_Chm,                                     &
               perSpc       = perSpc,                                        &
               N            = N,                                             &
               name         = metaDataId,                                    &
               desc         = desc,                                          &
               diagName     = thisSpcName,                                   &
               diagDesc     = thisSpcDesc                                   )

          ! Add field to registry
          CALL Registry_AddField( &
               Input_Opt    = Input_Opt,                                     &
               Registry     = State_Chm%Registry,                            &
               State        = State_Chm%State,                               &
               Variable     = TRIM( thisSpcName ),                           &
               Description  = TRIM( thisSpcDesc ),                           &
               Units        = TRIM( units       ),                           &
               OnLevelEdges = onEdges,                                       &
               Data2d_4     = Ptr2Data(:,:,N),                               &
               RC           = RC                                            )

          ! Trap potential errors
          IF ( RC /= GC_SUCCESS ) THEN
             ErrMsg = TRIM( ErrMsg_reg ) // TRIM( MetadataID ) //            &
                      '; Abnormal exit from routine "Registry_AddField"!'
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN
          ENDIF

       ENDDO

    !------------------------------------------------------------------------
    ! If not tied to species then simply register the single field
    !------------------------------------------------------------------------
    ELSE IF ( hasSpeciesDim == 0 ) THEN

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
       CALL Registry_AddField(                                               &
            Input_Opt    = Input_Opt,                                        &
            Registry     = State_Chm%Registry,                               &
            State        = State_Chm%State,                                  &
            Variable     = TRIM( metadataID ),                               &
            Description  = TRIM( desc       ),                               &
            Units        = TRIM( units      ),                               &
            Data3d_4     = Ptr2Data,                                         &
            OnLevelEdges = onEdges,                                          &
            RC           = RC                                               )

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

       ! Error msg
       ErrMsg = 'Handling of PerSpc metadata ' // TRIM(perSpc) //    &
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
! !IROUTINE: Register_ChmField_R4_4D
!
! !DESCRIPTION: Registers a 4-dimensional, 4-byte real array field
!  of the State\_Chm object.  This allows the diagnostic modules get
!  a pointer to the field by searching on the field name.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Register_ChmField_R4_4D( Input_Opt,  metadataID, Ptr2Data,     &
                                      State_Chm,  RC,         Ncat         )
!
! !USES:
!
    USE Input_Opt_Mod,      ONLY : OptInput
    USE Registry_Params_Mod
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput),   INTENT(IN)  :: Input_Opt         ! Input Options
    CHARACTER(LEN=*), INTENT(IN)  :: metadataID        ! State_Chm field Id
    REAL(f4),         POINTER     :: Ptr2Data(:,:,:,:) ! Pointer to data
    TYPE(ChmState),   INTENT(IN)  :: State_Chm         ! Chemistry State
    INTEGER,          OPTIONAL    :: Ncat              ! Category index
!
! !OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(OUT) :: RC                ! Success or failure?
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
    ! Scalars
    LOGICAL            :: found,       onEdges
    INTEGER            :: N,           rank,        hasSpeciesDim
    INTEGER            :: type,        vloc,        nSlots

    ! Strings
    CHARACTER(LEN=255) :: errMsg_reg,  thisLoc,     desc
    CHARACTER(LEN=255) :: thisSpcName, thisSpcDesc, perSpc
    CHARACTER(LEN=255) :: units
    CHARACTER(LEN=512) :: errMsg

    !========================================================================
    ! Initialize
    !========================================================================
    RC         = GC_SUCCESS
    errMsg     = ''
    errMsg_reg = 'Error encountered while registering State_Chm%'
    thisLoc    = &
       ' -> at Register_ChmField_R4_4D (in Headers/state_chm_mod.F90)'

    !========================================================================
    ! Get metadata
    !========================================================================
    CALL Get_MetaData_State_Chm(                                             &
         am_I_Root  = Input_Opt%amIRoot,                                     &
         metadataId = metadataId,                                            &
         found      = found,                                                 &
         desc       = desc,                                                  &
         units      = units,                                                 &
         rank       = rank,                                                  &
         type       = type,                                                  &
         vloc       = vloc,                                                  &
         perSpc     = perSpc,                                                &
         RC         = RC                                                    )

    ! Trap potential errors
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = TRIM( ErrMsg_reg ) // TRIM( MetadataID ) //                  &
                '; Abnormal exit from routine "Get_Metadata_State_Chm"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    ! Test if the quantity has a species dimension
    hasSpeciesDim = Test_For_Species_Dim( perSpc )

    !------------------------------------------------------------------------
    ! Check that metadata consistent with data size
    !------------------------------------------------------------------------
    IF ( rank /= 3 ) THEN
       ErrMsg = 'Data dims and metadata rank do not match for ' &
                // TRIM(metadataID)
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    ! Is the data placed on level edges?
    onEdges = ( VLoc == VLocationEdge )

    !------------------------------------------------------------------------
    ! If tied to a given category, only registry that one
    !------------------------------------------------------------------------
    IF ( PRESENT( Ncat ) ) THEN

       ! Add field to registry
       CALL Registry_AddField(                                               &
            Input_Opt    = Input_Opt,                                        &
            Registry     = State_Chm%Registry,                               &
            State        = State_Chm%State,                                  &
            Variable     = TRIM( metadataID ),                               &
            Description  = TRIM( desc       ),                               &
            Units        = TRIM( units      ),                               &
            OnLevelEdges = onEdges,                                          &
            Data3d_4     = Ptr2Data(:,:,:,Ncat),                             &
            RC           = RC                                               )

       ! Trap potential errors
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = TRIM( ErrMsg_reg ) // TRIM( MetadataID ) //               &
                   '; Abnormal exit from Registry_AddField!'
          CALL GC_Error( ErrMsg_reg, RC, ThisLoc )
          RETURN
       ENDIF

    !------------------------------------------------------------------------
    ! If tied to a particular species, register each species individually
    !------------------------------------------------------------------------
    ELSE IF ( hasSpeciesDim == 1 ) THEN

       ! Check that metadata consistent with data size
       IF ( rank /= 3 ) THEN
          ErrMsg = 'Data dims and metadata rank do not match for '           &
                   // TRIM(metadataID)
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

       ! Get the length of the species-based dimension
       nSlots = Get_NumSlots( perSpc, State_Chm )

       ! Loop over all species
       DO N = 1, nSlots

          ! Append the species name to the diagnostic name & description
          CALL Get_Diagnostic_Name(                                          &
               State_Chm    = State_Chm,                                     &
               perSpc       = perSpc,                                        &
               N            = N,                                             &
               name         = TRIM( metaDataId ),                            &
               desc         = TRIM( desc       ),                            &
               diagName     = thisSpcName,                                   &
               diagDesc     = thisSpcDesc                                   )

          ! Add field to registry
          CALL Registry_AddField( &
               Input_Opt    = Input_Opt,                                     &
               Registry     = State_Chm%Registry,                            &
               State        = State_Chm%State,                               &
               Variable     = TRIM( thisSpcName ),                           &
               Description  = TRIM( thisSpcDesc ),                           &
               Units        = TRIM( units       ),                           &
               OnLevelEdges = onEdges,                                       &
               Data3d_4     = Ptr2Data(:,:,:,N),                             &
               RC           = RC                                            )

          ! Trap potential errors
          IF ( RC /= GC_SUCCESS ) THEN
             ErrMsg = TRIM( ErrMsg_reg ) // TRIM( MetadataID ) //            &
                      '; Abnormal exit from routine "Registry_AddField"!'
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN
          ENDIF

       ENDDO

    !------------------------------------------------------------------------
    ! Otherwise, exit with error
    !------------------------------------------------------------------------
    ELSE

       ! Error msg
       ErrMsg = 'Handling of PerSpc metadata ' // TRIM(perSpc) // &
                ' is not implemented for this combo of data type and size!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN

    ENDIF

  END SUBROUTINE Register_ChmField_R4_4D
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Register_ChmField_R8_2D
!
! !DESCRIPTION: Registers a 2-dimensional, 8-byte real array field
!  of the State\_Chm object.  This allows the diagnostic modules get
!  a pointer to the field by searching on the field name.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Register_ChmField_R8_2D( Input_Opt, metadataID, Ptr2Data,       &
                                      State_Chm, RC                         )
!
! !USES:
!
    USE Input_Opt_Mod,      ONLY : OptInput
    USE Registry_Params_Mod
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput),   INTENT(IN)  :: Input_Opt       ! Input Options
    CHARACTER(LEN=*), INTENT(IN)  :: metadataID      ! State_Chm field ID
    REAL(f8),         POINTER     :: Ptr2Data(:,:)   ! Pointer to data
    TYPE(ChmState),   INTENT(IN)  :: State_Chm       ! Chemistry State
!
! !OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(OUT) :: RC              ! Success or failure?
!
! !REMARKS:
!
! !REVISION HISTORY:
!  21 Sep 2020 - E. Lundgren - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    LOGICAL            :: found,       onEdges
    INTEGER            :: N,           rank,        hasSpeciesDim
    INTEGER            :: type,        vloc

    ! Strings
    CHARACTER(LEN=255) :: errMsg_reg,  thisLoc,     desc
    CHARACTER(LEN=255) :: thisSpcName, thisSpcDesc, perSpc
    CHARACTER(LEN=255) :: units
    CHARACTER(LEN=512) :: errMsg

    !========================================================================
    ! Initialize
    !========================================================================
    RC         = GC_SUCCESS
    errMsg     = ''
    errMsg_reg = 'Error encountered while registering State_Chm%'
    thisLoc    = &
       ' -> at Register_ChmField_R8_2D (in Headers/state_chm_mod.F90)'

    !========================================================================
    ! Get metadata
    !========================================================================
    CALL Get_MetaData_State_Chm(                                             &
         am_I_Root  = Input_Opt%amIRoot,                                     &
         metadataId = metadataId,                                            &
         found      = found,                                                 &
         desc       = desc,                                                  &
         units      = units,                                                 &
         rank       = rank,                                                  &
         type       = type,                                                  &
         vloc       = vloc,                                                  &
         perSpc     = perSpc,                                                &
         RC         = RC                                                    )

    ! Trap potential errors
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = TRIM( ErrMsg_reg ) // TRIM( MetadataID ) //                  &
                '; Abnormal exit from routine "Get_Metadata_State_Chm"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    ! Is the data placed on vertical edges?
    onEdges = ( vLoc == VLocationEdge )

    ! Test if the quantity has a species dimension
    hasSpeciesDim = Test_for_Species_Dim( perSpc )

    !------------------------------------------------------------------------
    ! If not tied to species then simply register the single field
    !------------------------------------------------------------------------
    IF ( hasSpeciesDim == 0 ) THEN

       ! Check that metadata consistent with data size
       IF ( rank /= 2 ) THEN
          ErrMsg = 'Data dims and metadata rank do not match for '           &
                   // TRIM(metadataID)
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

       ! Add field to registry
       CALL Registry_AddField(                                               &
            Input_Opt    = Input_Opt,                                        &
            Registry     = State_Chm%Registry,                               &
            State        = State_Chm%State,                                  &
            Variable     = TRIM( metadataID ),                               &
            Description  = TRIM( desc       ),                               &
            Units        = TRIM( units      ),                               &
            Data2d_8     = Ptr2Data,                                         &
            RC           = RC                                               )

       ! Trap potential errors
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = TRIM( ErrMsg_reg ) // TRIM( MetadataID ) //               &
                  '; Abnormal exit from routine "Registry_AddField"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

    !------------------------------------------------------------------------
    ! Otherwise exit with error
    !------------------------------------------------------------------------
    ELSE

       ! Error msg
       ErrMsg = 'Handling of PerSpc metadata ' // TRIM(perSpc) //    &
                ' is not implemented for this combo of data type and size'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN

    ENDIF

  END SUBROUTINE Register_ChmField_R8_2D
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Register_ChmField_R8_3D
!
! !DESCRIPTION: Registers a 3-dimensional, 8-byte floating point array field
!  of the State\_Chm object.  This allows the diagnostic modules get a pointer
!  to the field by searching on the field name.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Register_ChmField_R8_3D( Input_Opt, metadataID, Ptr2Data,       &
                                      State_Chm, RC,         nCat           )
!
! !USES:
!
    USE Input_Opt_Mod,      ONLY : OptInput
    USE Registry_Params_Mod
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput),    INTENT(IN)  :: Input_Opt       ! Input Options object
    CHARACTER(LEN=*),  INTENT(IN)  :: metadataID      ! State_Chm field ID
    REAL(f8),          POINTER     :: Ptr2Data(:,:,:) ! Pointer to data
    TYPE(ChmState),    INTENT(IN)  :: State_Chm       ! Chemistry State
    INTEGER,           OPTIONAL    :: nCat            ! Category index
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,           INTENT(OUT) :: RC              ! Success or failure?
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
    ! Scalars
    LOGICAL            :: found,       onEdges
    INTEGER            :: N,           rank,        hasSpeciesDim
    INTEGER            :: type,        vloc,        nSlots

    ! Strings
    CHARACTER(LEN=255) :: errMsg_reg,  thisLoc,     desc
    CHARACTER(LEN=255) :: thisSpcName, thisSpcDesc, perSpc
    CHARACTER(LEN=255) :: units
    CHARACTER(LEN=512) :: errMsg

    !========================================================================
    ! Initialize
    !========================================================================
    RC         = GC_SUCCESS
    ErrMsg     = ''
    ErrMsg_reg = 'Error encountered while registering State_Chm%'
    ThisLoc    = &
       ' -> at Register_ChmField_R8_3D (in Headers/state_chm_mod.F90)'

    !========================================================================
    ! Get metadata
    !========================================================================
    CALL Get_MetaData_State_Chm(                                             &
         am_I_Root  = Input_Opt%amIRoot,                                     &
         metadataId = metadataId,                                            &
         found      = found,                                                 &
         desc       = desc,                                                  &
         units      = units,                                                 &
         rank       = rank,                                                  &
         type       = type,                                                  &
         vloc       = vloc,                                                  &
         perSpc     = perSpc,                                                &
         RC         = RC                                                    )

    ! Trap potential errors
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = TRIM( ErrMsg_reg ) // TRIM( MetadataID ) //                  &
                '; Abnormal exit from routine "Get_Metadata_State_Chm"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    ! Is the data placed on vertical edges?
    onEdges = ( vLoc == VLocationEdge )

    ! Test if the data has a species dimension
    hasSpeciesDim = Test_for_Species_Dim( perSpc )

    !------------------------------------------------------------------------
    ! If tied to a given category, only registry that one
    !------------------------------------------------------------------------
    IF ( PRESENT( nCat ) ) THEN

       ! Check that metadata consistent with data size
       IF ( rank /= 2 ) THEN
          ErrMsg = 'Data dims and metadata rank do not match for '           &
                   // TRIM(metadataID)
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

       ! Add field to registry
       CALL Registry_AddField(                                               &
            Input_Opt    = Input_Opt,                                        &
            Registry     = State_Chm%Registry,                               &
            State        = State_Chm%State,                                  &
            Variable     = TRIM( metadataID ),                               &
            Description  = TRIM( desc       ),                               &
            Units        = TRIM( units      ),                               &
            OnLevelEdges = onEdges,                                          &
            Data2d_8     = Ptr2Data(:,:,nCat),                               &
            RC           = RC                                               )

       ! Trap potential errors
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = TRIM( ErrMsg_reg ) // TRIM( MetadataID ) //               &
                   '; Abnormal exit from Registry_AddField!'
          CALL GC_Error( ErrMsg_reg, RC, ThisLoc )
          RETURN
       ENDIF

    !------------------------------------------------------------------------
    ! If tied to a particular species, register each species individually
    !------------------------------------------------------------------------
    ELSE IF ( hasSpeciesDim == 1 ) THEN

       ! Check that metadata consistent with data size
       IF ( rank /= 2 ) THEN
          ErrMsg = 'Data dims and metadata rank do not match for '           &
                   // TRIM(metadataID)
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

       ! Get the length of the species-based dimension
       nSlots = Get_NumSlots( perSpc, State_Chm )

       ! Loop over all species
       DO N = 1, nSlots

          ! Append the species name to the diagnostic name & description
          CALL Get_Diagnostic_Name(                                          &
               State_Chm    = State_Chm,                                     &
               perSpc       = perSpc,                                        &
               N            = N,                                             &
               name         = TRIM( metaDataId ),                            &
               desc         = TRIM( desc       ),                            &
               diagName     = thisSpcName,                                   &
               diagDesc     = thisSpcDesc                                   )

          ! Add field to registry
          CALL Registry_AddField( &
               Input_Opt    = Input_Opt,                                     &
               Registry     = State_Chm%Registry,                            &
               State        = State_Chm%State,                               &
               Variable     = TRIM( thisSpcName ),                           &
               Description  = TRIM( thisSpcDesc ),                           &
               Units        = TRIM( units       ),                           &
               OnLevelEdges = onEdges,                                       &
               Data2d_8     = Ptr2Data(:,:,N),                               &
               RC           = RC                                            )

          ! Trap potential errors
          IF ( RC /= GC_SUCCESS ) THEN
             ErrMsg = TRIM( ErrMsg_reg ) // TRIM( MetadataID ) //            &
                      '; Abnormal exit from routine "Registry_AddField"!'
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN
          ENDIF

       ENDDO

    !------------------------------------------------------------------------
    ! If not tied to species then simply register the single field
    !------------------------------------------------------------------------
    ELSE IF ( hasSpeciesDim == 0 ) THEN

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
       CALL Registry_AddField(                                               &
            Input_Opt    = Input_Opt,                                        &
            Registry     = State_Chm%Registry,                               &
            State        = State_Chm%State,                                  &
            Variable     = TRIM( metadataID ),                               &
            Description  = TRIM( desc       ),                               &
            Units        = TRIM( units      ),                               &
            Data3d_8     = Ptr2Data,                                         &
            OnLevelEdges = onEdges,                                          &
            RC           = RC                                               )

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

       ! Error msg
       ErrMsg = 'Handling of PerSpc metadata ' // TRIM(perSpc) //    &
                ' is not implemented for this combo of data type and size'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN

    ENDIF

  END SUBROUTINE Register_ChmField_R8_3D
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Register_ChmField_R8_4D
!
! !DESCRIPTION: Registers a 4-dimensional, 8-byte real array field
!  of the State\_Chm object.  This allows the diagnostic modules get
!  a pointer to the field by searching on the field name.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Register_ChmField_R8_4D( Input_Opt,  metadataID, Ptr2Data,      &
                                      State_Chm,  RC,         Ncat          )
!
! !USES:
!
    USE Input_Opt_Mod,      ONLY : OptInput
    USE Registry_Params_Mod
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput),   INTENT(IN)    :: Input_Opt         ! Input Options
    CHARACTER(LEN=*), INTENT(IN)    :: metadataID        ! State_Chm field ID
    REAL(f8),         POINTER       :: Ptr2Data(:,:,:,:) ! Pointer to data
    TYPE(ChmState),   INTENT(IN)    :: State_Chm         ! Chemistry State
    INTEGER,          OPTIONAL      :: Ncat              ! Category index
!
! !OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(OUT)   :: RC                ! Success or failure?
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
    ! Scalars
    LOGICAL            :: found,       onEdges
    INTEGER            :: N,           rank,        hasSpeciesDim
    INTEGER            :: type,        vloc,        nSlots

    ! Strings
    CHARACTER(LEN=255) :: errMsg_reg,  thisLoc,     desc
    CHARACTER(LEN=255) :: thisSpcName, thisSpcDesc, perSpc
    CHARACTER(LEN=255) :: units
    CHARACTER(LEN=512) :: errMsg

    !========================================================================
    ! Initialize
    !========================================================================
    RC         = GC_SUCCESS
    errMsg     = ''
    errMsg_reg = 'Error encountered while registering State_Chm%'
    thisLoc    = &
       ' -> at Register_ChmField_R8_4D (in Headers/state_chm_mod.F90)'

    !========================================================================
    ! Get metadata
    !========================================================================
    CALL Get_MetaData_State_Chm(                                             &
         am_I_Root  = Input_Opt%amIRoot,                                     &
         metadataId = metadataId,                                            &
         found      = found,                                                 &
         desc       = desc,                                                  &
         units      = units,                                                 &
         rank       = rank,                                                  &
         type       = type,                                                  &
         vloc       = vloc,                                                  &
         perSpc     = perSpc,                                                &
         RC         = RC                                                    )

    ! Trap potential errors
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = TRIM( ErrMsg_reg ) // TRIM( MetadataID ) //                  &
                '; Abnormal exit from routine "Get_Metadata_State_Chm"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    !------------------------------------------------------------------------
    ! Check that metadata consistent with data size
    !------------------------------------------------------------------------
    IF ( rank /= 3 ) THEN
       ErrMsg = 'Data dims and metadata rank do not match for ' &
                // TRIM(metadataID)
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    ! Is the data placed on level edges?
    onEdges = ( VLoc == VLocationEdge )

    ! Test if the quantity has a species dimension
    hasSpeciesDim = Test_for_Species_Dim( perSpc )

    !------------------------------------------------------------------------
    ! If tied to a given category, only registry that one
    !------------------------------------------------------------------------
    IF ( PRESENT( Ncat ) ) THEN

       ! Add field to registry
       CALL Registry_AddField(                                               &
            Input_Opt    = Input_Opt,                                        &
            Registry     = State_Chm%Registry,                               &
            State        = State_Chm%State,                                  &
            Variable     = TRIM( metadataID ),                               &
            Description  = TRIM( desc       ),                               &
            Units        = TRIM( units      ),                               &
            OnLevelEdges = onEdges,                                          &
            Data3d_8     = Ptr2Data(:,:,:,Ncat),                             &
            RC           = RC                                               )

       ! Trap potential errors
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = TRIM( ErrMsg_reg ) // TRIM( MetadataID ) //               &
                   '; Abnormal exit from Registry_AddField!'
          CALL GC_Error( ErrMsg_reg, RC, ThisLoc )
          RETURN
       ENDIF

    !------------------------------------------------------------------------
    ! If tied to a particular species, register each species individually
    !------------------------------------------------------------------------
    ELSE IF ( hasSpeciesDim == 1 ) THEN

       ! Get the length of the species-based dimension
       nSlots = Get_NumSlots( perSpc, State_Chm )

       ! Loop over all species
       DO N = 1, nSlots

          ! Append the species name to the diagnostic name & description
          CALL Get_Diagnostic_Name(                                          &
               State_Chm    = State_Chm,                                     &
               perSpc       = perSpc,                                        &
               N            = N,                                             &
               name         = TRIM( metaDataId ),                            &
               desc         = TRIM( desc       ),                            &
               diagName     = thisSpcName,                                   &
               diagDesc     = thisSpcDesc                                   )

          ! Add field to registry
          CALL Registry_AddField( &
               Input_Opt    = Input_Opt,                                     &
               Registry     = State_Chm%Registry,                            &
               State        = State_Chm%State,                               &
               Variable     = TRIM( thisSpcName ),                           &
               Description  = TRIM( thisSpcDesc ),                           &
               Units        = TRIM( units       ),                           &
               OnLevelEdges = onEdges,                                       &
               Data3d_8     = Ptr2Data(:,:,:,N),                             &
               RC           = RC                                            )

          ! Trap potential errors
          IF ( RC /= GC_SUCCESS ) THEN
             ErrMsg = TRIM( ErrMsg_reg ) // TRIM( MetadataID ) //            &
                      '; Abnormal exit from routine "Registry_AddField"!'
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN
          ENDIF

       ENDDO

    !------------------------------------------------------------------------
    ! Otherwise, exit with error
    !------------------------------------------------------------------------
    ELSE

       ! Error msg
       ErrMsg = 'Handling of PerSpc metadata ' // TRIM(perSpc ) // &
                ' is not implemented for this combo of data type and size!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN

    ENDIF

  END SUBROUTINE Register_ChmField_R8_4D
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
!   'A' or 'a' : Returns advected species index
!   'D' or 'd' : Returns dry-deposition species index
!   'F' or 'f' : Returns KPP fixed species index
!   'G' or 'g' : Returns gas-phase species index
!   'H' or 'h' : Returns hygroscopic-growth species index
!   'K' or 'k' : Returns KPP main species index
!   'N' or 'n' : Returns radionuclide species index
!   'P' or 'p' : Returns photolysis species index
!   'S' or 's' : Returns main species index (aka "ModelId")
!   'V' or 'v' : Returns KPP variable species index
!   'W' or 'w' : Returns wet-deposition species index
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

       ! Radionuclide chemical species ID
       CASE( 'N', 'n' )
          Indx = SpcDataLocal(N)%Info%RadNuclId
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
    IF ( Input_Opt%ITS_A_FULLCHEM_SIM .or. Input_Opt%ITS_A_MERCURY_SIM ) THEN

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

    ELSE IF ( Input_Opt%ITS_A_TAGO3_SIM         .or.                         &
              Input_Opt%ITS_A_CARBON_SIM      ) THEN

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
    IF ( Input_Opt%ITS_A_FULLCHEM_SIM .or. Input_Opt%ITS_A_MERCURY_SIM ) THEN

       !--------------------------------------------------------------------
       ! Full-chemistry simulations
       !--------------------------------------------------------------------

       ! Loop over the number of prod/loss species
       DO N = 1, nFam

          ! Get the KPP prod/loss species from the FAM_NAMES
          ! array in the gckpp_Parameters.F90 module.
          ! NOTE: This is the KPP ID number (index of "VAR" array)
          ! and not the GEOS-Chem "main" species index!!!
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

    ELSE IF ( Input_Opt%ITS_A_TAGO3_SIM         .or.                         &
              Input_Opt%ITS_A_CARBON_SIM      ) THEN

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
