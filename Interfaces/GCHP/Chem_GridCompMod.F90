#include "MAPL_Generic.h"

!------------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 910.1 and      !
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: GEOSCHEMchem_GridCompMod
!
! !DESCRIPTION: GEOSCHEMchem_GridComp is an ESMF5 gridded component 
! implementing the GEOS-Chem chemistry and related processes, including
! dry deposition, emissions, and wet deposition. In addition, the
! parameterizations for PBL mixing and convection as used in GEOS-Chem
! can be invoked by enabling the corresponding option in the GEOS-Chem
! input file (input.geos). In this case, the corresponding GEOS-5
! process must NOT be applied to the GC tracers, i.e. the tracers must 
! not be friendly to turbulence (if PBL mixing is used) and/or moist
! (for convection).
!\\
!\\
! This gridded component contains three run phases: 
!
!  -1: Phase -1 is the standard setting in GCHP. It executes all components. 
!      Phase is -1 if number of phases is set to 1 in config file GCHP.rc.
!
!   1: Phase 1 is used in GEOS-5. It executes convection, dry deposition,
!      and emissions and should be called before surface processes/turbulence. 
! 
!   2: Phase 2 is used in GEOS-5. It performs chemistry, and wet deposition, 
!      and should be called after turbulence.
!
! GEOS-5 only:
! All GEOS-Chem species are stored in the GEOSCHEMchem internal state object
! in units of kg/kg total. Since molecular weights are missing for 
! non-transported species (SPC_<XXX>), a value of 1g/mol is used when 
! converting molec/cm3 to kg/kg total. Hence, the concentrations of all
! non-transported species are in units of 'kg of species with molecular
! weight of 1g/mol per kg air'.
!\\
!\\
! !INTERFACE:
!
#ifdef MODEL_GEOS
MODULE GEOSCHEMchem_GridCompMod
#else
MODULE Chem_GridCompMod
#endif
!
! !USES:
!
  USE CMN_Size_Mod
  USE ESMF                                           ! ESMF library
  USE MAPL_Mod                                       ! MAPL library
  USE Charpak_Mod                                    ! String functions
  USE Hco_Types_Mod, ONLY : ConfigObj
  USE Input_Opt_Mod                                  ! Input Options obj
  USE GIGC_Chunk_Mod                                 ! GIGC IRF methods
  USE GIGC_HistoryExports_Mod
#if !defined( MODEL_GEOS )
  USE GIGC_ProviderServices_Mod
#endif
  USE ErrCode_Mod                                    ! Error numbers
  USE State_Chm_Mod                                  ! Chemistry State obj
  USE State_Diag_Mod                                 ! Diagnostics State obj
  USE State_Grid_Mod                                 ! Grid State obj
  USE State_Met_Mod                                  ! Meteorology State obj
  USE Species_Mod,   ONLY : Species

#if defined( MODEL_GEOS )
  USE MAPL_ConstantsMod                   ! Doesn't seem to be used. Needed?
  USE Chem_Mod                            ! Chemistry Base Class (chem_mie?)
  USE Chem_GroupMod                       ! For family transport
  USE PHYSCONSTANTS
#endif

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC   :: SetServices    ! Sets ESMF entry points
!
! !PRIVATE MEMBER FUNCTIONS:
!
  PRIVATE  :: Initialize_    ! Init method
  PRIVATE  :: Run1           ! Run wrapper phase 1
  PRIVATE  :: Run2           ! Run wrapper phase 2
  PRIVATE  :: Run_           ! Run method
  PRIVATE  :: Finalize_      ! Finalize method
  PRIVATE  :: Extract_       ! Get values from ESMF
#if defined( MODEL_GEOS )
  PRIVATE  :: Roundoff       ! Truncates a number
  PRIVATE  :: Print_Mean_OH  ! Mean OH lifetime
  PRIVATE  :: GlobalSum      ! Sums across PETs
#endif
!
! !PRIVATE TYPES:
!
  ! Legacy state
  TYPE GEOSCHEM_State
     PRIVATE
     TYPE(ESMF_Config)             :: myCF           ! Private ESMF Config obj
  END TYPE GEOSCHEM_State

  ! Hook for the ESMF
  TYPE GEOSCHEM_Wrap
     TYPE(GEOSCHEM_State), POINTER :: PTR => null()  ! Ptr to GEOSCHEM_State
  END TYPE GEOSCHEM_Wrap

  ! For passing from internal state to Chm_State and vice versa
#if defined( MODEL_GEOS )
  TYPE Int2SpcMap
     CHARACTER(LEN=255)            :: Name
     INTEGER                       :: ID
     REAL, POINTER                 :: Internal(:,:,:) => NULL()
  END TYPE Int2SpcMap
#else
  ! SDE 2016-03-28: Assume that Int2Spc is held as v/v dry
  TYPE Int2SpcMap
     CHARACTER(LEN=255)            :: TrcName
     INTEGER                       :: TrcID
     REAL                          :: TCVV
     REAL(ESMF_KIND_R8), POINTER   :: Internal(:,:,:) => NULL()
  END TYPE Int2SpcMap
#endif

  ! For mapping State_Chm%Tracers/Species arrays onto the internal state.
  TYPE(Int2SpcMap), POINTER        :: Int2Spc(:) => NULL()

  ! Objects for GEOS-Chem
  TYPE(OptInput)                   :: Input_Opt      ! Input Options
  TYPE(ChmState)                   :: State_Chm      ! Chemistry state
  TYPE(DgnState)                   :: State_Diag     ! Diagnostics state 
  TYPE(GrdState)                   :: State_Grid     ! Grid state 
  TYPE(MetState)                   :: State_Met      ! Meteorology state
  TYPE(Species),          POINTER  :: ThisSpc => NULL()
  TYPE(HistoryConfigObj), POINTER  :: HistoryConfig
  TYPE(ConfigObj),        POINTER  :: HcoConfig

#if defined( MODEL_GEOS )
  ! Is GEOS-Chem the provider for AERO, RATS, and/or Analysis OX? 
  LOGICAL                          :: DoAERO
  LOGICAL                          :: DoRATS
  LOGICAL                          :: DoANOX

  ! To set ozone from PCHEM 
  LOGICAL                          :: LANAO3
  LOGICAL                          :: LPCHEMO3
  INTEGER                          :: ANAO3L1
  INTEGER                          :: ANAO3L2
  INTEGER                          :: ANAO3L3
  INTEGER                          :: ANAO3L4
  REAL                             :: ANAO3FR 
  CHARACTER(LEN=ESMF_MAXSTR)       :: ANAO3FILE
#else
  LOGICAL                          :: isProvider ! provider to AERO, RATS, ANOX?
  LOGICAL                          :: calcOzone  ! if PTR_GCCTO3 is associated
#endif

  ! Number of run phases, 1 or 2. Set in the rc file; else default is 2.
  INTEGER                          :: NPHASE

  ! Is this being run as a CTM?
  INTEGER                          :: IsCTM

  ! Memory debug level
  INTEGER                          :: MemDebugLevel

#if defined( MODEL_GEOS )
  ! GEOS-5 only
  ! Flag to initialize species concentrations from external fields. Read 
  ! through GEOSCHEMchem_GridComp.rc. If true, initial species concentrations
  ! are read from an external file instead of taken from the internal state. 
  ! The field names in the external file are expected to be 'SPC_<XXX>'. 
  ! This option can be used to initialize a simulation using a restart file
  ! from a 'GEOS-Chem classic' CTM simulation.
  LOGICAL                          :: InitFromFile
#endif

#if defined( MODEL_GEOS )
  ! GEOS-5 only (also in gigc_providerservices_mod but don't use that yet):
  ! List here GEOS-Chem tracer names and corresponding names to be assigned
  ! to the AERO bundle (if GC is the AERO provider). The names in the AERO
  ! bundle must be the names that are expected by the irradiation component:
  ! - OCphobic, OCphilic, BCphobic, and BCphilic for hydrophobic and hydrophilic
  !   organic and black carbon, respectively
  ! - SO4 for SO4
  ! - du001 - du005 for the following five dust bins (see DU_GridComp.rc in
  !   GOCART):
  !   radius_lower: 0.1 1.0 1.8 3.0 6.0
  !   radius_upper: 1.0 1.8 3.0 6.0 10.0
  !
  !   The GEOS-Chem dust bins are: 
  !   Reff: 0.7 1.4 2.4 4.5
  !   Those become simply mapped onto the GOCART dust bins 1-4 
  !   (du001 ... du004).
  !
  ! - ss001-ss005 for the following five sea salt aerosol bins 
  !   (see SS_GridComp.rc
  !   in GOCART):
  !   radius_lower: 0.03 0.1 0.5 1.5 5.0
  !   radius_upper: 0.1  0.5 1.5 5.0 10.0
  !
  !   The GEOS-Chem sea salt aerosols are (SALA and SALC):
  !   radius_lower: 0.01 0.5
  !   radius_upper: 0.5  8.0
  !   SALA becomes mapped onto ss001 and ss002, and SALC onto ss003, ss004, 
  !   ss005. For now, we assume uniform size distribution within the 
  !   GEOS-Chem bins, i.e. the GEOS-Chem size bins are evenly split into the 
  !   GOCART bins. The fractions can be specified below.
  !   At some point, we may revisit these fractions (at least take into 
  !   account the log-normal behavior of the aerosol distribution)
  INTEGER, PARAMETER           :: NumAERO = 11

  CHARACTER(LEN=ESMF_MAXSTR)   :: GcNames(NumAero) = &
                                  (/ 'DST1',  'DST2',  'DST3',  'DST4',     &
                                     'SALA',  'SALC',  'BCPO',  'BCPI',     &
                                     'OCPO',  'OCPI',  'SO4 '                /)

  CHARACTER(LEN=ESMF_MAXSTR)   :: AeroNames(NumAero) = &
                         (/ 'du001   ', 'du002   ', 'du003   ', 'du004   ', &
                            'ss001   ', 'ss003   ', 'BCphobic', 'BCphilic', &
                            'OCphobic', 'OCphilic', 'SO4     '               /)

  ! Fraction of SALA in ss001 and ss002, respectively
  CHARACTER(LEN=ESMF_MAXSTR)   :: SALAnames(2) = (/ 'ss001', 'ss002' /)
  REAL, PARAMETER              :: SALAsplit(2) = (/  0.2,     0.8    /)

  ! Fraction of SALC in ss003, ss004, and ss005.
  CHARACTER(LEN=ESMF_MAXSTR) :: SALCnames(3) = (/ 'ss003', 'ss004' , 'ss005' /)
  REAL, PARAMETER            :: SALCsplit(3) = (/  0.13,    0.47,     0.4    /) 

  CHARACTER(LEN=ESMF_MAXSTR)   :: DST4names(2) = (/ 'du004', 'du005' /)
  REAL, PARAMETER              :: DST4split(2) = (/  1.00,    0.0    /) 

  ! Molecular weights (g/mol) used by GOCART
  REAL                         :: GocartMW(NumAero) = &
                                  (/ 100.0,  100.0,  100.0,  100.0 ,     &
                                      58.0,   58.0,  180.0,  180.0 ,     &
                                     180.0,  180.0,  132.0                /)

  CHARACTER(LEN=15), PARAMETER :: COLLIST(8) = (/ 'NO2', 'O3',   'CH4', 'CO', &
                                                  'BrO', 'CH2O', 'SO2', 'IO'  /)
#endif
 
  ! Pointers to import, export and internal state data. Declare them as 
  ! module variables so that we have to assign them only on first call.
#if !defined( MODEL_GEOS )
  ! NOTE: Any provider-related exports (e.g. H2O_TEND) are now handled within
  ! gigc_providerservices_mod.F90. Pointers are manually declared there and
  ! those declared in the .h file included below are not used. (ewl, 11/3/2017)
#endif

#if defined( MODEL_GEOS )
# include "GEOSCHEMCHEM_DeclarePointer___.h"
#else
# include "GIGCchem_DeclarePointer___.h"
#endif

#if !defined( MODEL_GEOS )
  ! GCHP only
  ! Use archived convection fields?
  ! If the attribute 'ARCHIVED_CONV' in the GEOS-Chem configuration file is set
  ! to '1', GEOS-Chem will use archived convection fields, imported through
  ! ExtData (ExtData must contain an entry for each of the pointers defined
  ! below). These data fields are then passed to the GEOS-Chem meteorlogical 
  ! state instead of the (instantaneous) fields imported from MOIST. The 
  ! fields imported through ExtData must be named 'ARCHIVED_PFI_CN', 
  ! 'ARCHIVED_PFL_CN', etc.
  LOGICAL           :: ArchivedConv
  REAL, POINTER     :: PTR_ARCHIVED_PFI_CN (:,:,:) => NULL()
  REAL, POINTER     :: PTR_ARCHIVED_PFL_CN (:,:,:) => NULL()
  REAL, POINTER     :: PTR_ARCHIVED_CNV_MFC(:,:,:) => NULL()
  REAL, POINTER     :: PTR_ARCHIVED_CNV_MFD(:,:,:) => NULL()
  REAL, POINTER     :: PTR_ARCHIVED_CNV_CVW(:,:,:) => NULL()
  REAL, POINTER     :: PTR_ARCHIVED_DQRC   (:,:,:) => NULL()
  REAL, POINTER     :: PTR_ARCHIVED_REV_CN (:,:,:) => NULL()
  REAL, POINTER     :: PTR_ARCHIVED_T      (:,:,:) => NULL()
#endif

#if defined( MODEL_GEOS )
  !! GEOS-5 only (also in gigc_providerservices but don't use yet):
  !! -RATS:
  !REAL, POINTER     :: CH4     (:,:,:) => NULL()
  !REAL, POINTER     :: N2O     (:,:,:) => NULL()
  !REAL, POINTER     :: CFC11   (:,:,:) => NULL()
  !REAL, POINTER     :: CFC12   (:,:,:) => NULL()
  !REAL, POINTER     :: HCFC22  (:,:,:) => NULL()
  !
  !! -Corresponding pointers to internal state. We now use these variables 
  !!  instead of the auto-generated pointers (GEOSCHEMCHEM_DeclarePointer___.h) 
  !!  to avoid compilation errors if these species are not defined in 
  !!  GEOS-Chem (e.g. for specialty sims). 
  !REAL, POINTER     :: PTR_O3      (:,:,:) => NULL()
  !REAL, POINTER     :: PTR_CH4     (:,:,:) => NULL()
  !REAL, POINTER     :: PTR_N2O     (:,:,:) => NULL()
  !REAL, POINTER     :: PTR_CFC11   (:,:,:) => NULL()
  !REAL, POINTER     :: PTR_CFC12   (:,:,:) => NULL()
  !REAL, POINTER     :: PTR_HCFC22  (:,:,:) => NULL()
  !REAL, POINTER     :: PTR_H2O     (:,:,:) => NULL()
  !
  !! GCCTO3 and GCCTTO3 are the pointers to the corresponding export state fields
  !REAL, POINTER     :: PTR_GCCTO3 (:,:) => NULL()
  !REAL, POINTER     :: PTR_GCCTTO3(:,:) => NULL()

  ! Perturb ozone?
  LOGICAL           :: PerturbO3
  LOGICAL           :: PerturbCO
  REAL              :: FIXPERT       ! Fixed perturbation 
  REAL, PARAMETER   :: MAXPERT = 0.1 ! Maximum perturbation

  ! For NOx diagnostics
  INTEGER           :: id_rc_no  = -1
  INTEGER           :: id_rc_no2 = -1
  INTEGER           :: id_jno2   = -1

  ! Mie table
  TYPE(Chem_Mie)     :: geoschemMieTable(2)
  INTEGER, PARAMETER :: instanceComputational = 1
  INTEGER, PARAMETER :: instanceData          = 2
#endif
!
! !REMARKS:
!  Developed for GEOS-5 release Fortuna 2.0 and later.
!                                                                             .
!  NOTES: 
!  - The abbreviation "PET" stands for "Persistent Execution Thread". 
!    It is a synomym for CPU.
#if defined( MODEL_GEOS )
!  - The internal state now holds all species in v/v
!    dry air. For backwards compatibility, the internal
!    state names can still use the tracer prefix, e.g.
!    'TRC_NO2' instead of 'SPC_NO2'.
!  - The current implementation is not the most 
!    efficient version of GEOS-Chem. In particular, the 
!    chemical solver does not use the archived 
!    integration time steps. This makes the code more
!    stable, but also slower. To remove this safety
!    net, activate 'RCNTRL(3) = HSAVE_KPP(I,J,L)'
!    around line 765 in gc_column/GeosCore/flexchem_mod.F90.
!
!  TODO:
!  - Spin off HEMCO code (use HEMCO_GridComp instead)
!  - Activate RCNTRL(3) in flexchem_mod
#endif
!
! !REVISION HISTORY:
!  06 Dec 2009 - A. da Silva - Initial version
!  10 Oct 2012 - R. Yantosca - Now references GC_Utils.F90
!  10 Oct 2012 - R. Yantosca - Updated for GEOS-Chem v9-01-03
!  16 Oct 2012 - R. Yantosca - Rename GC_MET object to State_Met
!  16 Oct 2012 - R. Yantosca - Rename GC_STATE object to State_Chm
!  17 Oct 2012 - R. Yantosca - Removed some old "column code" stuff
!  22 Oct 2012 - R. Yantosca - Now references renamed gigc_* modules
!  01 Nov 2012 - R. Yantosca - Now references gigc_input_opt_mod.F90
!  07 Nov 2012 - R. Yantosca - Removed Setup_GeoLoc_ routine
!  07 Nov 2012 - R. Yantosca - Now read placeholder values for input.geos
!  08 Nov 2012 - R. Yantosca - Now initialize Input_Opt%MAX_DIAG field
!  15 Mar 2013 - R. Yantosca - Remove IDENT object and Error_Trap routine
!  25 Mar 2014 - E. Nielsen  - ESMF-5
!  22 Sep 2014 - C. Keller   - Added two run phases
!  17 Oct 2014 - C. Keller   - Various updates to fill provider fields.
!  26 Nov 2014 - C. Keller   - Added H2O_HIST and O3_HIST. 
!  22 Feb 2015 - C. Keller   - Now check if geoschemchem_import_rst exist
#if defined( MODEL_GEOS )
!  08 May 2015 - C. Keller   - Update on Int2Trc. Also added Int2Spc to make
!                              sure that the internal state always contains the
!                              most current species values. This is critical 
!                              for the checkpoint files.
!  16 Sep 2015 - C. Keller   - Added H2O2s, SO2s, DRY_TOTN and WET_TOTN to 
!                              internal state. Values are passed from/to internal
!                              state to GEOS-Chem arrays in Include_Before_Run.H
!                              and Include_After_Run.H.
!  27 Feb 2017 - C. Keller   - Update to GEOS-Chem v11. GEOS-Chem does not
!                              distinguish between tracers and species any more.
!  07 Mar 2017 - C. Keller   - Species unit in internal state is now kg/kg total.
!  21 Apr 2017 - C. Keller   - Update to v11-02.
#else
!  06 Jun 2016 - M. Yannetti - Added Get_Transport.
!  19 Dec 2016 - M. Long     - Update for v11-01k
!  01 Sep 2017 - E. Lundgren - Enable automation of GCHP diagnostics
!  19 Sep 2017 - E. Lundgren - Remove Get_Transport
!  02 Nov 2017 - E. Lundgren - Remove unused private functions roundoff, 
!                              globalsum, and print_mean_oh
!  06 Nov 2017 - E. Lundgren - Abstract provider services to new module
!                              gigc_providerservices_mod.F90
#endif
!  22 May 2019 - M. Sulprizio- Added State_Grid object; Remove unused variables
!                              I_LO, J_LO, I_HI, and J_HI  
!EOP
!------------------------------------------------------------------------------
!BOC
CONTAINS
!EOC
!------------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 910.1 and      !
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: SetServices
!
! !DESCRIPTION: The SetServices routine does the following:
!
! \begin{itemize}
! \item Defines the Initialize method for the GEOSCHEMchem gridded component
! \item Defines the Run methods for the GEOSCHEMchem gridded component
!       (phase 1 and phase 2).
! \item Defines the Finalize method for the GEOSCHEMchem gridded component
! \item Attaches an internal state (which holds a private ESMF Config object)
!       to the GEOSCHEMchem gridded component.
! \end{itemize}
!
! !INTERFACE:
!
  SUBROUTINE SetServices( GC, RC )
!
! !USES:
!
    USE HCOI_ESMF_MOD,        ONLY : HCO_SetServices
    USE GCKPP_Model
    USE CHARPAK_MOD,          ONLY : STRSPLIT, CSTRIP
    USE inquireMod,           ONLY : findFreeLUN
    USE FILE_MOD,             ONLY : IOERROR
#if defined( MODEL_GEOS )
    USE SPECIES_DATABASE_MOD, ONLY : Spc_Info
    USE CMN_FJX_MOD
    USE GCKPP_Monitor
    USE GCKPP_Parameters
    USE Precision_Mod
#endif
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ESMF_GridComp), INTENT(INOUT) :: GC       ! Ref to this GridComp
!
! !OUTPUT PARAMETERS:
!
    INTEGER,             INTENT(OUT)   :: RC       ! Success or failure
!
! !REMARKS:
!  ESMF can only attach one Config object per Gridded Component.  The 
!  Config object that is defined from the "MAPL.rc" resource file is 
!  directly attached to the GEOSCHEMchem gridded component.  
!                                                                             .
!  To attach the Config object defined from the "GEOSCHEMchem_GridComp.rc" 
!  resource file, we must first create a derived type with a pointer to
!  the Config object, then attach that to the gridded component as an
!  "internal state" (also called "legacy state").
!
! !REVISION HISTORY:
!  06 Dec 2009 - A. da Silva - Initial version
!  07 Apr 2010 - R. Yantosca - Updated comments, cosmetic changes 
!  22 Sep 2014 - C. Keller   - Added two run phases
!  07 Aug 2017 - E. Lundgren - Add Olson and CHRL imports
!  14 Jul 2017 - E. Lundgren - Read simulation type to determine whether to
!                              add KPP species to the internal state
!  01 Sep 2017 - E. Lundgren - Call new subroutine HistoryExports_SetServices
!                              for GEOS-Chem state object diagnostics
!  12 Sep 2017 - E. Lundgren - Use species prefix "SPFX" from gigc_types_mod.F90
!  06 Nov 2017 - E. Lundgren - Abstract provider services to new module
!                              gigc_providerservices_mod.F90
!  08 Mar 2018 - E. Lundgren - "SPFX" now retrieved from gigc_historyexports_mod
!EOP
!------------------------------------------------------------------------------
!BOC
! 
! !LOCAL VARIABLES:
!
    TYPE(ESMF_CONFIG)             :: CF
    TYPE(GEOSCHEM_State), POINTER :: myState       ! Legacy state
    TYPE(GEOSCHEM_Wrap)           :: wrap          ! Wrapper for myState
    CHARACTER(LEN=ESMF_MAXSTR)    :: compName      ! Gridded Component name
    CHARACTER(LEN=ESMF_MAXSTR)    :: COMP_NAME     ! This syntax for mapl_acg.pl
    CHARACTER(LEN=ESMF_MAXSTR)    :: HcoConfigFile ! HEMCO configuration file
    CHARACTER(LEN=ESMF_MAXSTR)    :: SpcName       ! Registered species name
    CHARACTER(LEN=40)             :: AdvSpc(500)
    CHARACTER(LEN=255)            :: LINE, MSG, SUBSTRS(500)
    INTEGER                       :: N, I, J, IU_GEOS, IOS
    INTEGER                       :: Nadv, landTypeInt
    LOGICAL                       :: FOUND 
    LOGICAL                       :: EOF
    CHARACTER(LEN=60)             :: rstFile, landTypeStr, importName
    INTEGER                       :: restartAttr
    CHARACTER(LEN=ESMF_MAXSTR)    :: HistoryConfigFile ! HISTORY config file
    INTEGER                       :: T

#if !defined( MODEL_GEOS )
    TYPE(MAPL_MetaComp),  POINTER :: STATE => NULL()
#endif

#if defined( MODEL_GEOS )
    CHARACTER(LEN=ESMF_MAXSTR)    :: ProviderName  ! Provider name
    CHARACTER(LEN=ESMF_MAXSTR)    :: LongName      ! Long name for diagnostics
    CHARACTER(LEN=ESMF_MAXSTR)    :: ShortName
    CHARACTER(LEN=3)              :: III
    CHARACTER(LEN=6)              :: T_FJX
    CHARACTER(LEN=50)             :: T_REACT
    CHARACTER(LEN=120)            :: CLINE
    REAL(fp)                      :: F_FJX
    INTEGER                       :: JJ, NUNIT
    CHARACTER(LEN=255)            :: MYFRIENDLIES
    CHARACTER(LEN=31)             :: iName
    CHARACTER(LEN=127)            :: FullName, Formula 
    LOGICAL                       :: FriendDyn, FriendTurb
#endif

    __Iam__('SetServices')

    !=======================================================================
    ! Set services begins here 
    !=======================================================================

    ! Set up traceback info
    CALL ESMF_GridCompGet( GC, name=compName, __RC__ )

    ! NOTE: We need to use COMP_NAME for mapl_acg.pl script
    COMP_NAME = TRIM( compName )

    ! Identify this routine to MAPL
    Iam = TRIM(compName)//'::SetServices'
    
    !=======================================================================
    ! Wrap internal state for storing in this gridded component
    ! Rename this to a "legacy state"
    !=======================================================================
    ALLOCATE( myState, stat=STATUS )
    _VERIFY(STATUS)
    wrap%ptr => myState

    !=======================================================================
    ! Define an ESMF Config object from the Resource file and set it 
    ! as an "internal state" of the GEOSCHEMchem gridded component
    !=======================================================================
    myState%myCF = ESMF_ConfigCreate(__RC__)

#if defined( MODEL_GEOS )
    call ESMF_ConfigLoadFile( myState%myCF, 'GEOSCHEMchem_GridComp.rc', __RC__)
#else
    call ESMF_ConfigLoadFile( myState%myCF, 'GCHP.rc', __RC__)
#endif

#if !defined( MODEL_GEOS )
    ! Get generic state object
    CALL MAPL_GetObjectFromGC( GC, STATE, __RC__ )
    call MAPL_GetResource( STATE, IsCTM, label='GEOSChem_CTM:', & 
                           default=1, rc=status )
    _VERIFY(STATUS)
#endif

    !=======================================================================
    !                 %%% ESMF Functional Services %%%
    !=======================================================================

    ! Set the Initialize, Run, Finalize entry points
    CALL MAPL_GridCompSetEntryPoint( GC, ESMF_METHOD_INITIALIZE,  &
                                     Initialize_, __RC__ )
#if defined( MODEL_GEOS )
    CALL MAPL_GridCompSetEntryPoint( GC, ESMF_METHOD_RUN, Run1, __RC__ )
#endif
    CALL MAPL_GridCompSetEntryPoint( GC, ESMF_METHOD_RUN, Run2, __RC__ )
    CALL MAPL_GridCompSetEntryPoint( GC, ESMF_METHOD_FINALIZE,  &
                                     Finalize_, __RC__ )
        
    ! Store internal state with Config object in the gridded component
    CALL ESMF_UserCompSetInternalState( GC, 'GEOSCHEM_State', wrap, STATUS )
    _VERIFY(STATUS)
  
#if defined(MODEL_GEOS)
    ! new after abstracting to gigc_providerservices_mod, but do not use yet:
    !    CALL Provider_SetServices( MAPL_am_I_Root(), GC, isProvider, __RC__ )
    ! GEOS-5 (also in gigc_providerservices but do not use yet):

    ! Check if GEOS-Chem is set as the AERO and/or RATS provider
    ! ----------------------------------------------------------

    ! Get configuration
    CALL ESMF_GridCompGet( GC, CONFIG = CF, __RC__ )

    ! See if GC is the AERO provider
    DoAERO = .FALSE.
    CALL ESMF_ConfigGetAttribute( CF, ProviderName,       &
                                  Label="AERO_PROVIDER:", &
                                  Default="PCHEM",        &
                                  __RC__                   )
    IF ( ProviderName == "GEOSCHEMCHEM" ) DoAERO = .TRUE.

    ! See if GC is the RATS provider
    DoRATS = .FALSE.
    CALL ESMF_ConfigGetAttribute( CF, ProviderName,       &
                                  Label="RATS_PROVIDER:", &
                                  Default="PCHEM",        &
                                  __RC__                   )
    IF ( ProviderName == "GEOSCHEMCHEM" ) DoRATS = .TRUE.

    ! See if GC is the Analysis OX provider
    DoANOX = .FALSE.
    CALL ESMF_ConfigGetAttribute( CF, ProviderName,              &
                                  Label="ANALYSIS_OX_PROVIDER:", &
                                  Default="PCHEM",               &
                                  __RC__                          )
    IF ( ProviderName == "GEOSCHEMCHEM" ) DoANOX = .TRUE.
#endif

    !=======================================================================
    !                    %%% MAPL Data Services %%%
    !=======================================================================
!EOC
!BOP
!
! !IMPORT STATE:
!
#if defined( MODEL_GEOS )
#   include "GEOSCHEMCHEM_ImportSpec___.h"
#else
#   include "GIGCchem_ImportSpec___.h"
#endif

#if !defined( MODEL_GEOS )
    call MAPL_AddImportSpec(GC, &
       SHORT_NAME         = 'PLE',  &
       LONG_NAME          = 'pressure_level_edges',  &
       UNITS              = 'Pa', &
       PRECISION          = ESMF_KIND_R8, &
       DIMS               = MAPL_DimsHorzVert,    &
       VLOCATION          = MAPL_VLocationEdge,    &
                                                      RC=STATUS  )
    _VERIFY(STATUS)

    call MAPL_AddImportSpec(GC, &
       SHORT_NAME         = 'DryPLE',  &
       LONG_NAME          = 'dry_pressure_level_edges',  &
       UNITS              = 'Pa', &
       PRECISION          = ESMF_KIND_R8, &
       DIMS               = MAPL_DimsHorzVert,    &
       VLOCATION          = MAPL_VLocationEdge,    &
                                                      RC=STATUS  )
    _VERIFY(STATUS)
#endif

!
! !INTERNAL STATE:
!
#if defined( MODEL_GEOS )
#   include "GEOSCHEMCHEM_InternalSpec___.h"
#else
#   include "GIGCchem_InternalSpec___.h"
#endif

#if !defined( MODEL_GEOS )
    ! Determine if using a restart file for the internal state. Setting
    ! the GIGCchem_INTERNAL_RESTART_FILE to +none in GCHP.rc indicates
    ! skipping the restart file. Species concentrations will be retrieved
    ! from the species database, overwriting MAPL-assigned default values.
    CALL ESMF_ConfigGetAttribute( myState%myCF, rstFile, &
                                  Label = "GIGCchem_INTERNAL_RESTART_FILE:",&
                                  __RC__ ) 
    IF ( TRIM(rstFile) == '+none' ) THEN
       restartAttr = MAPL_RestartSkipInitial ! file does not exist;
                                             ! use background values
    ELSE
       restartAttr = MAPL_RestartOptional    ! try to read species from file;
                                             ! use background vals if not found
    ENDIF
#endif

!-- Read in species from input.geos and set FRIENDLYTO
    ! ewl TODO: This works but is not ideal. Look into how to remove it.

    ! Open input.geos and read a lines until hit advected species menu
    IU_GEOS = findFreeLun()
    OPEN( IU_GEOS, FILE='input.geos', STATUS='OLD', IOSTAT=IOS )
    IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_GEOS, 'READ_SPECIES_FROM_FILE:1' )
    DO
      READ( IU_GEOS, '(a)', IOSTAT=IOS ) LINE
      IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_GEOS, 'READ_SPECIES_FROM_FILE:2' )
      IF ( INDEX( LINE, 'ADVECTED SPECIES MENU' ) > 0 ) EXIT
    ENDDO

    ! Read in all advected species names and add them to internal state
    NADV=0
    DO WHILE( INDEX( LINE, 'TRANSPORT MENU' ) .le. 0) 
       READ( IU_GEOS, '(a)', IOSTAT=IOS ) LINE
       EOF = IOS < 0
       IF ( EOF ) RETURN
       CALL STRSPLIT( LINE(26:), ' ', SUBSTRS, N )
       IF ( INDEX( LINE, 'Species name' ) > 0 ) THEN

#if defined( MODEL_GEOS )
          ! Define friendliness to dynamics / turbulence 
          FriendDyn  = .TRUE.
          FriendTurb = .TRUE.
          IF ( N > 1 ) THEN
             IF ( TRIM(SUBSTRS(2))=="N" .OR. TRIM(SUBSTRS(2))=="NO" ) THEN
                FriendDyn = .FALSE.
             ENDIF
             IF ( N > 2 ) THEN
                IF ( TRIM(SUBSTRS(3))=="N" .OR. TRIM(SUBSTRS(3))=="NO" ) THEN
                   FriendTurb = .FALSE.
                ENDIF
             ELSE
                FriendTurb = FriendDyn
             ENDIF
          ENDIF
          IF (       FriendDyn .AND.       FriendTurb ) &
                                     MYFRIENDLIES = 'DYNAMICS:TURBULENCE'
          IF (       FriendDyn .AND. .NOT. FriendTurb ) &
                                     MYFRIENDLIES = 'DYNAMICS'
          IF ( .NOT. FriendDyn .AND.       FriendTurb ) &
                                     MYFRIENDLIES = 'TURBULENCE'
          IF ( .NOT. FriendDyn .AND. .NOT. FriendTurb ) &
                                     MYFRIENDLIES = 'GEOSCHEMCHEM'

          ! Get long name
          iName = TRIM(SUBSTRS(1))
          CALL Spc_Info ( am_I_Root = MAPL_am_I_Root(), &
                          iName=iName,                  &
                          KppSpcID=-1,                  &
                          oDiagName = FullName,         &
                          oFormula = Formula,           &
                          Found=Found,                  &
                          Underscores = .TRUE.,         &
                          RC = RC )
          IF ( .NOT. FOUND ) FullName = TRIM(SUBSTRS(1))

          call MAPL_AddInternalSpec(GC, &
               SHORT_NAME         = TRIM(TPFX)//TRIM(SUBSTRS(1)), &
               LONG_NAME          = TRIM(FullName)//                &
                                    ' mass mixing ratio total air', &
               UNITS              = 'kg kg-1',                &
               DIMS               = MAPL_DimsHorzVert,        &
               VLOCATION          = MAPL_VLocationCenter,     &
               !!!PRECISION          = ESMF_KIND_R8,             &
               FRIENDLYTO         = TRIM(MYFRIENDLIES),       &
               RC                 = RC  )
          Nadv = Nadv+1
          AdvSpc(Nadv) = TRIM(SUBSTRS(1))
          
          ! verbose
          if(MAPL_am_I_Root()) write(*,*) &
                   'GCC added to internal: TRC_'//TRIM(SUBSTRS(1)), &
                   '; Friends: ', TRIM(MYFRIENDLIES)
#else
          call MAPL_AddInternalSpec(GC, &
              SHORT_NAME         = TRIM(SPFX) // TRIM(SUBSTRS(1)),  &
              LONG_NAME          = TRIM(SUBSTRS(1)),  &
              UNITS              = 'mol mol-1', &
              DIMS               = MAPL_DimsHorzVert,    &
              VLOCATION          = MAPL_VLocationCenter,    &
              PRECISION          = ESMF_KIND_R8, &
              FRIENDLYTO         = 'DYNAMICS:TURBULENCE:MOIST',  &
              RESTART            = restartAttr, &
              RC                 = RC  )
         NADV = NADV+1
         AdvSpc(NADV) = TRIM(SUBSTRS(1))
#endif

       ENDIF
    ENDDO
    CLOSE( IU_GEOS )
!
!    CALL READ_SPECIES_FROM_FILE( GC, MAPL_am_I_Root(), restartAttr, &
!                                 AdvSpc, Nadv, SimType, RC )

!-- Add all additional species in KPP (careful not to add dummy species)
    IF ( Nadv > 50 ) THEN ! Exclude specialty sims
       DO I=1,NSPEC
          FOUND = .false.
       
          ! Skip dummy RR species for prod/loss diagnostic (mps, 8/23/16)
          SpcName = ADJUSTL( Spc_Names(I) )
          IF ( SpcName(1:2) == 'RR' ) CYCLE
       
          DO J=1,Nadv !Size of AdvSpc
#if defined( MODEL_GEOS )
             IF (trim(AdvSpc(J)) .eq. trim(SpcName)) THEN
                FOUND = .true.
                EXIT
             ENDIF
#else
             IF (trim(AdvSpc(J)) .eq. trim(SpcName)) FOUND = .true.
#endif
          END DO
       
#if defined( MODEL_GEOS )
          ! Add non-advected species to internal state 
          IF ( .NOT. Found ) THEN 

             ! Get long name
             iName = TRIM(SpcName)
             CALL Spc_Info ( am_I_Root = MAPL_am_I_Root(), iName=iName, &
                             KppSpcID=-1, oDiagName = FullName, &
                             oFormula = Formula, &
                             Found=Found, Underscores = .TRUE., RC = RC )
             IF ( .NOT. FOUND ) FullName = TRIM(SpcName)

             ! Error trap for POx and LOx. Their species names in the internal
             ! state must be all caps
             ! (ckeller, 3/11/19)
             IF ( TRIM(SpcName) == 'POx' ) SpcName = 'POX'
             IF ( TRIM(SpcName) == 'LOx' ) SpcName = 'LOX'

             ! Set some long names manually ...
             SELECT CASE ( TRIM(SpcName) )
                CASE ('OH')
                   FullName = 'Hydroxyl radical (OH, MW = 17.01 g mol-1)'
                CASE ('HO2')
                   FullName = 'Hydroperoxyl radical (HO2, MW = 33.01 g mol-1)'
                CASE ('O')
                   FullName = 'Molecular oxygen (O, MW = 16.01 g mol-1)'
             END SELECT

             CALL MAPL_AddInternalSpec(GC,                 &
                SHORT_NAME         = TRIM(SPFX)//TRIM(SpcName),      &
                LONG_NAME          = TRIM(FullName)//                &
                                     ' mass mixing ratio total air', &
                UNITS              = 'kg kg-1',            &
                !!!PRECISION          = ESMF_KIND_R8,          &
                DIMS               = MAPL_DimsHorzVert,    &
                FRIENDLYTO         = COMP_NAME,            &
                VLOCATION          = MAPL_VLocationCenter, &
                                                     __RC__ )
             ! verbose
             if(MAPL_am_I_Root()) write(*,*)  &
                       'GCC added to internal: '//TRIM(SPFX)//TRIM(SpcName)

             ! Also add an export as dry air
             CALL MAPL_AddExportSpec(GC,                               &
                SHORT_NAME         = TRIM(GPFX)//TRIM(SpcName)//'dry', &
                LONG_NAME          = TRIM(FullName)//                  &
                                     ' volume mixing ratio dry air',   &
                UNITS              = 'mol mol-1',                      &
                DIMS               = MAPL_DimsHorzVert,                &
                VLOCATION          = MAPL_VLocationCenter,             &
                                                                 __RC__ )
#else
          IF (Found .neqv. .true.) Then
          call MAPL_AddInternalSpec(GC, &
               SHORT_NAME         = TRIM(SPFX) // SpcName,  &
               LONG_NAME          = SpcName,  &
               UNITS              = 'mol mol-1', &
               PRECISION          = ESMF_KIND_R8, &
               DIMS               = MAPL_DimsHorzVert,    &
               VLOCATION          = MAPL_VLocationCenter,    &
               RESTART            = restartAttr,    &
               RC                 = STATUS  )
#endif
          Endif
       ENDDO
    ENDIF

#if !defined( MODEL_GEOS )
    ! Add other internal state variables as real8
    call MAPL_AddInternalSpec(GC, &
       SHORT_NAME         = 'DryDepNitrogen',  &
       LONG_NAME          = 'Dry deposited nitrogen',  &
       UNITS              = 'cm-2s-1', &
       DIMS               = MAPL_DimsHorzOnly,    &
       VLOCATION          = MAPL_VLocationCenter,    &
       PRECISION          = ESMF_KIND_R8, &
       FRIENDLYTO         = trim(COMP_NAME),    &
                                                      RC=STATUS  )
    _VERIFY(STATUS)

    call MAPL_AddInternalSpec(GC, &
       SHORT_NAME         = 'WetDepNitrogen',  &
       LONG_NAME          = 'Wet deposited nitrogen',  &
       UNITS              = 'cm-2s-1', &
       DIMS               = MAPL_DimsHorzOnly,    &
       VLOCATION          = MAPL_VLocationCenter,    &
       PRECISION          = ESMF_KIND_R8, &
       FRIENDLYTO         = trim(COMP_NAME),    &
                                                      RC=STATUS  )
    _VERIFY(STATUS)

    call MAPL_AddInternalSpec(GC, &
       SHORT_NAME         = 'H2O2AfterChem',  &
       LONG_NAME          = 'Soluble fraction H2O2',  &
       UNITS              = 'vv-1', &
       DIMS               = MAPL_DimsHorzVert,    &
       VLOCATION          = MAPL_VLocationCenter,    &
       PRECISION          = ESMF_KIND_R8, &
       FRIENDLYTO         = trim(COMP_NAME),    &
                                                      RC=STATUS  )
    _VERIFY(STATUS)

    call MAPL_AddInternalSpec(GC, &
       SHORT_NAME         = 'SO2AfterChem',  &
       LONG_NAME          = 'Soluble fraction SO2',  &
       UNITS              = 'vv-1', &
       DIMS               = MAPL_DimsHorzVert,    &
       VLOCATION          = MAPL_VLocationCenter,    &
       PRECISION          = ESMF_KIND_R8, &
       FRIENDLYTO         = trim(COMP_NAME),    &
                                                      RC=STATUS  )
    _VERIFY(STATUS)

    call MAPL_AddInternalSpec(GC, &
       SHORT_NAME         = 'KPPHvalue',  &
       LONG_NAME          = 'HSAVE for KPP',  &
       UNITS              = '1', &
       DIMS               = MAPL_DimsHorzVert,    &
       VLOCATION          = MAPL_VLocationCenter,    &
       PRECISION          = ESMF_KIND_R8, &
       FRIENDLYTO         = trim(COMP_NAME),    &
                                                      RC=STATUS  )
    _VERIFY(STATUS)

    call MAPL_AddInternalSpec(GC, &
       SHORT_NAME         = 'DELP_DRY',  &
       LONG_NAME          = 'Delta dry pressure across box',  &
       UNITS              = 'hPa', &
       DIMS               = MAPL_DimsHorzVert,    &
       VLOCATION          = MAPL_VLocationCenter,    &
       PRECISION          = ESMF_KIND_R8, &
       FRIENDLYTO         = trim(COMP_NAME),    &
                                                      RC=STATUS  )
    _VERIFY(STATUS)
#endif

#if defined( MODEL_GEOS )
!-- Add two extra advected species for use in family transport  (Manyin)

          CALL MAPL_AddInternalSpec(GC,                                    &
             SHORT_NAME         = 'TRC_Bry',                               &
             LONG_NAME          = 'Bromine group for use in transport',    &
             UNITS              = 'kg kg-1',                               &
!!!          PRECISION          = ESMF_KIND_R8,                            &
             DIMS               = MAPL_DimsHorzVert,                       &
             FRIENDLYTO         = 'DYNAMICS',                              &
             RESTART            = MAPL_RestartSkip,                        &
             VLOCATION          = MAPL_VLocationCenter,                    &
                                                  __RC__ )
          if(MAPL_am_I_Root()) write(*,*) 'GCC added to internal: TRC_Bry; Friendly to: DYNAMICS'

          CALL MAPL_AddInternalSpec(GC,                                    &
             SHORT_NAME         = 'TRC_Cly',                               &
             LONG_NAME          = 'Chlorine group for use in transport',   &
             UNITS              = 'kg kg-1',                               &
!!!          PRECISION          = ESMF_KIND_R8,                            &
             DIMS               = MAPL_DimsHorzVert,                       &
             FRIENDLYTO         = 'DYNAMICS',                              &
             RESTART            = MAPL_RestartSkip,                        &
             VLOCATION          = MAPL_VLocationCenter,                    &
                                                  __RC__ )
          if(MAPL_am_I_Root()) write(*,*) 'GCC added to internal: TRC_Cly; Friendly to: DYNAMICS'
!
!-- Add OX to the internal state if GEOS-Chem is the analysis OX provider
!   Make sure it is friendly to ANALYSIS! In GEOS-Chem, OX is diagnosed 
!   from O3, O3P, and O1D.
     IF ( DoANOX ) THEN
        CALL MAPL_AddInternalSpec(GC,                                    &
           SHORT_NAME         = 'OX',                                    &
           LONG_NAME          = 'odd_oxygen_volume_mixing_ratio',        &
           UNITS              = 'mol mol-1',                             &
           DIMS               = MAPL_DimsHorzVert,                       &
           FRIENDLYTO         = 'ANALYSIS:DYNAMICS:TURBULENCE:MOIST',    &
           RESTART            = MAPL_RestartSkip,                        &
           VLOCATION          = MAPL_VLocationCenter,                    &
                                                  __RC__ )
        if(MAPL_am_I_Root()) write(*,*) 'OX added to internal: Friendly to: ANALYSIS, DYNAMICS, TURBULENCE'

        ! Add additional RATs/ANOX exports
        call MAPL_AddExportSpec(GC,                                  &
           SHORT_NAME         = 'OX_TEND',                           &
           LONG_NAME          = 'tendency_of_odd_oxygen_mixing_ratio_due_to_chemistry', &
           UNITS              = 'mol mol-1 s-1',                     &
           DIMS               = MAPL_DimsHorzVert,                   &
           VLOCATION          = MAPL_VLocationCenter,                &
                                                     __RC__ )

        call MAPL_AddExportSpec(GC,                                  &
           SHORT_NAME         = 'O3',                                &
           LONG_NAME          = 'ozone_mass_mixing_ratio',           &
           UNITS              = 'kg kg-1',                           &
           DIMS               = MAPL_DimsHorzVert,                   &
           VLOCATION          = MAPL_VLocationCenter,                &
                                                     __RC__ )

        call MAPL_AddExportSpec(GC,                                  &
           SHORT_NAME         = 'O3PPMV',                            &
           LONG_NAME          = 'ozone_volume_mixing_ratio',         &
           UNITS              = 'ppmv',                              &
           DIMS               = MAPL_DimsHorzVert,                   &
           VLOCATION          = MAPL_VLocationCenter,                &
                                                 __RC__  )
     ENDIF
#endif

!
! !EXTERNAL STATE:
!
#if defined( MODEL_GEOS )
#   include "GEOSCHEMCHEM_ExportSpec___.h"
#else
#   include "GIGCchem_ExportSpec___.h"
#endif

    ! Read HISTORY config file and add exports for unique items
    CALL ESMF_ConfigGetAttribute( myState%myCF, HistoryConfigFile, &
                                  Label="HISTORY_CONFIG:",         &
                                  Default="HISTORY.rc", __RC__ )
    CALL HistoryExports_SetServices( MAPL_am_I_Root(), HistoryConfigFile, &
                                     GC, HistoryConfig, RC=STATUS )
    _VERIFY(STATUS)

!EOP
!BOC

#if defined( MODEL_GEOS )
    ! GEOS-5 only (should handle diags elsewhere, use State vars now):
    !-- Exports
    DO I=1,Nadv
       iName = TRIM(AdvSpc(I))
       CALL Spc_Info ( am_I_Root = MAPL_am_I_Root(), iName=iName, &
                       KppSpcID=-1, oDiagName = FullName, Found=Found, &
                       Underscores = .FALSE., RC = RC )
       IF ( .NOT. FOUND ) FullName = TRIM(SpcName)

       ! For all advected species, create placeholder export for deposition 
       ! diagnostics. These may not be defined for all species (some species 
       ! don't have dry or wet dep), but at this point the species metadata 
       ! are not yet fully defined. Convective wet deposition flux
       CALL MAPL_AddExportSpec(GC,                                            &
          SHORT_NAME         = 'WetLossConv_'//TRIM(AdvSpc(I)),               & 
          LONG_NAME          = TRIM(FullName)//' vertical integrated loss'//  &
                               ' in convective updrafts',                     &
          UNITS              = 'kg m-2 s-1',                                  &
          DIMS               = MAPL_DimsHorzOnly,                             &
          VLOCATION          = MAPL_VLocationNone,                            &
                                                                     __RC__ )
       ! 3D
       CALL MAPL_AddExportSpec(GC,                                            &
          SHORT_NAME         = 'WetLossConv3D_'//TRIM(AdvSpc(I)),             &
          LONG_NAME          = TRIM(FullName)//' loss in convective updrafts',&
          UNITS              = 'kg m-2 s-1',                                  &
          DIMS               = MAPL_DimsHorzVert,                             &
          VLOCATION          = MAPL_VLocationCenter,                          &
                                                                     __RC__ )
       ! Large scale wet deposition flux
       CALL MAPL_AddExportSpec(GC,                                            &
          SHORT_NAME         = 'WetLossLS_'//TRIM(AdvSpc(I)),                 & 
          LONG_NAME          = TRIM(FullName)//' vertical integrated loss'//  &
                               ' in large scale precipitation',               &
          UNITS              = 'kg m-2 s-1',                                  &
          DIMS               = MAPL_DimsHorzOnly,                             &
          VLOCATION          = MAPL_VLocationNone,                            &
                                                                     __RC__ )
       CALL MAPL_AddExportSpec(GC,                                            &
          SHORT_NAME         = 'WetLossLS3D_'//TRIM(AdvSpc(I)),               & 
          LONG_NAME          = TRIM(FullName)//                               &
                               ' loss in large scale precipitation',          &
          UNITS              = 'kg m-2 s-1',                                  &
          DIMS               = MAPL_DimsHorzVert,                             &
          VLOCATION          = MAPL_VLocationCenter,                          &
                                                                     __RC__ )
       ! Total wet deposition flux
       CALL MAPL_AddExportSpec(GC,                                            &
          SHORT_NAME         = 'WetLossTot_'//TRIM(AdvSpc(I)),                & 
          LONG_NAME          = TRIM(FullName)//' vertical integrated loss'//  &
                               '_due_to_wet_scavenging',                      &
          UNITS              = 'kg m-2 s-1',                                  &
          DIMS               = MAPL_DimsHorzOnly,                             &
          VLOCATION          = MAPL_VLocationNone,                            &
                                                                     __RC__ )
       CALL MAPL_AddExportSpec(GC,                                            &
          SHORT_NAME         = 'WetLossTot3D_'//TRIM(AdvSpc(I)),              & 
          LONG_NAME          = TRIM(FullName)//' loss due to wet scavenging', &

          UNITS              = 'kg m-2 s-1',                                  &
          DIMS               = MAPL_DimsHorzVert,                             &
          VLOCATION          = MAPL_VLocationCenter,                          &
                                                                     __RC__ )
       ! Dry deposition flux
       CALL MAPL_AddExportSpec(GC,                                            &
          SHORT_NAME         = 'DryDep_'//TRIM(AdvSpc(I)),                    & 
          LONG_NAME          = TRIM(FullName)//' dry deposition flux',        &
          UNITS              = 'molec cm-2 s-1',                              &
          DIMS               = MAPL_DimsHorzOnly,                             &
          VLOCATION          = MAPL_VLocationNone,                            &
                                                                     __RC__ )
       ! Dry deposition velocity
       CALL MAPL_AddExportSpec(GC,                                            &
          SHORT_NAME         = 'DryDepVel_'//TRIM(AdvSpc(I)),                 & 
          LONG_NAME          = TRIM(FullName)//' dry deposition velocity',    &
          UNITS              = 'cm s-1',                                      &
          DIMS               = MAPL_DimsHorzOnly,                             &
          VLOCATION          = MAPL_VLocationNone,                            &
                                                                     __RC__ )

       ! Also create export field in v/v dry air
       CALL MAPL_AddExportSpec(GC,                                            &
          SHORT_NAME         = TRIM(GPFX)//TRIM(AdvSpc(I))//'dry',            &
          LONG_NAME          = TRIM(FullName)//                               &
                               ' volume mixing ratio dry air',                &
          UNITS              = 'mol mol-1',                                   &
          DIMS               = MAPL_DimsHorzVert,                             &
          VLOCATION          = MAPL_VLocationCenter,                          &
                                                                     __RC__ )
       CALL MAPL_AddExportSpec(GC,                                            &
          SHORT_NAME         = TRIM(AdvSpc(I))//'vv_2m',                      & 
          LONG_NAME          = TRIM(FullName)//                               &
                               ' volume mixing ratio dry air at 2m',          &
          UNITS              = 'mol mol-1',                                   &
          DIMS               = MAPL_DimsHorzOnly,                             &
          VLOCATION          = MAPL_VLocationNone,                            &
                                                                     __RC__ )
       CALL MAPL_AddExportSpec(GC,                                            &
          SHORT_NAME         = TRIM(AdvSpc(I))//'vv_10m',                     & 
          LONG_NAME          = TRIM(FullName)//                               &
                               ' volume mixing ratio dry air at 10m',         &
          UNITS              = 'mol mol-1',                                   &
          DIMS               = MAPL_DimsHorzOnly,                             &
          VLOCATION          = MAPL_VLocationNone,                            &
                                                                     __RC__ )
    ENDDO

    ! Add exports for tropospheric and total column densities
    DO I=1,SIZE(COLLIST,1)
       SpcName = COLLIST(I)
       iName = TRIM(SpcName)
       CALL Spc_Info ( am_I_Root = MAPL_am_I_Root(), iName=iName, &
                       KppSpcID=-1, oDiagName = FullName,         &
                       Found = Found, Underscores = .TRUE., RC = RC )
       IF ( .NOT. Found ) FullName = TRIM(SpcName)
       CALL MAPL_AddExportSpec(GC,                                            &
          SHORT_NAME         = 'TOTCOL_'//TRIM(SpcName),                      & 
          LONG_NAME          = TRIM(FullName)//' total column density',        &
          UNITS              = '1.0e15 molec cm-2',                           &
          DIMS               = MAPL_DimsHorzOnly,                             &
          VLOCATION          = MAPL_VLocationNone,                            &
                                                                     __RC__ )
       CALL MAPL_AddExportSpec(GC,                                            &
          SHORT_NAME         = 'TROPCOL_'//TRIM(SpcName),                     & 
          LONG_NAME          = TRIM(FullName)//' tropospheric column density', &
          UNITS              = '1.0e15 molec cm-2',                           &
          DIMS               = MAPL_DimsHorzOnly,                             &
          VLOCATION          = MAPL_VLocationNone,                            &
                                                                     __RC__ )
    ENDDO
#endif

    !=======================================================================
    ! Add provider services, if any (AERO, RATS, Analysis Ox)
    !=======================================================================
#if defined( MODEL_GEOS )

    ! Add AERO and AERO_DP bundles to export state if GEOS-Chem is the 
    ! AERO provider
    ! ----------------------------------------------------------------
    IF ( DoAERO ) THEN
      
       ! The AERO bundle contains DUST, SALT, SO4, BC, and OC.
       ! These quantities will be obtained from the respective
       ! GEOS-Chem internal state quantities. 
       ! Fields are added to bundle in the initialize routine.
       call MAPL_AddExportSpec(GC,                                  &
          SHORT_NAME         = 'AERO',                              &
          LONG_NAME          = 'aerosol_mass_mixing_ratios',        &
          UNITS              = 'kg kg-1',                           &
          DIMS               = MAPL_DimsHorzVert,                   &
          VLOCATION          = MAPL_VLocationCenter,                &
          DATATYPE           = MAPL_StateItem,                      &
                                                            __RC__ )
       ! This bundle is needed by surface for snow albedo modification.
       ! At the moment, it is not filled by GEOS-Chem.
       call MAPL_AddExportSpec(GC,                                  &
          SHORT_NAME         = 'AERO_DP',                           &
          LONG_NAME          = 'aerosol_deposition',                &
          UNITS              = 'kg m-2 s-1',                        &
          DIMS               = MAPL_DimsHorzOnly,                   &
          DATATYPE           = MAPL_BundleItem,                     &
                                                            __RC__ )
      
       ! Fields of AERO_DP bundle:
       call MAPL_AddExportSpec(GC,                         &
          SHORT_NAME         = 'DUDP_DST1',                &
          LONG_NAME          = 'dust1_dry_depostion',      &
          UNITS              = 'kg m-2 s-1',               &
          DIMS               = MAPL_DimsHorzOnly,          &
          VLOCATION          = MAPL_VLocationNone,         &
                                                   __RC__ )
 
       call MAPL_AddExportSpec(GC,                         &
          SHORT_NAME         = 'DUDP_DST2',                &
          LONG_NAME          = 'dust2_dry_depostion',      &
          UNITS              = 'kg m-2 s-1',               &
          DIMS               = MAPL_DimsHorzOnly,          &
          VLOCATION          = MAPL_VLocationNone,         &
                                                   __RC__ )
  
       call MAPL_AddExportSpec(GC,                         &
          SHORT_NAME         = 'DUDP_DST3',                &
          LONG_NAME          = 'dust3_dry_depostion',      &
          UNITS              = 'kg m-2 s-1',               &
          DIMS               = MAPL_DimsHorzOnly,          &
          VLOCATION          = MAPL_VLocationNone,         &
                                                   __RC__ )
  
       call MAPL_AddExportSpec(GC,                         &
          SHORT_NAME         = 'DUDP_DST4',                &
          LONG_NAME          = 'dust4_dry_depostion',      &
          UNITS              = 'kg m-2 s-1',               &
          DIMS               = MAPL_DimsHorzOnly,          &
          VLOCATION          = MAPL_VLocationNone,         &
                                                   __RC__ )
 
       call MAPL_AddExportSpec(GC,                         &
          SHORT_NAME         = 'DUWT_DST1',                &
          LONG_NAME          = 'dust1_wet_depostion',      &
          UNITS              = 'kg m-2 s-1',               &
          DIMS               = MAPL_DimsHorzOnly,          &
          VLOCATION          = MAPL_VLocationNone,         &
                                                   __RC__ )
 
       call MAPL_AddExportSpec(GC,                         &
          SHORT_NAME         = 'DUWT_DST2',                &
          LONG_NAME          = 'dust2_wet_depostion',      &
          UNITS              = 'kg m-2 s-1',               &
          DIMS               = MAPL_DimsHorzOnly,          &
          VLOCATION          = MAPL_VLocationNone,         &
                                                   __RC__ )
  
       call MAPL_AddExportSpec(GC,                         &
          SHORT_NAME         = 'DUWT_DST3',                &
          LONG_NAME          = 'dust3_wet_depostion',      &
          UNITS              = 'kg m-2 s-1',               &
          DIMS               = MAPL_DimsHorzOnly,          &
          VLOCATION          = MAPL_VLocationNone,         &
                                                   __RC__ )
  
       call MAPL_AddExportSpec(GC,                         &
          SHORT_NAME         = 'DUWT_DST4',                &
          LONG_NAME          = 'dust4_wet_depostion',      &
          UNITS              = 'kg m-2 s-1',               &
          DIMS               = MAPL_DimsHorzOnly,          &
          VLOCATION          = MAPL_VLocationNone,         &
                                                   __RC__ )
 
       call MAPL_AddExportSpec(GC,                         &
          SHORT_NAME         = 'BCDP_BCPI',                &
          LONG_NAME          = 'BCPI_dry_depostion',       &
          UNITS              = 'kg m-2 s-1',               &
          DIMS               = MAPL_DimsHorzOnly,          &
          VLOCATION          = MAPL_VLocationNone,         &
                                                   __RC__ )
  
       call MAPL_AddExportSpec(GC,                         &
          SHORT_NAME         = 'BCDP_BCPO',                &
          LONG_NAME          = 'BCPO_dry_depostion',       &
          UNITS              = 'kg m-2 s-1',               &
          DIMS               = MAPL_DimsHorzOnly,          &
          VLOCATION          = MAPL_VLocationNone,         &
                                                   __RC__ )
  
       call MAPL_AddExportSpec(GC,                         &
          SHORT_NAME         = 'BCWT_BCPI',                &
          LONG_NAME          = 'BCPI_wet_depostion',       &
          UNITS              = 'kg m-2 s-1',               &
          DIMS               = MAPL_DimsHorzOnly,          &
          VLOCATION          = MAPL_VLocationNone,         &
                                                   __RC__ )
 
       call MAPL_AddExportSpec(GC,                         &
          SHORT_NAME         = 'BCWT_BCPO',                &
          LONG_NAME          = 'BCPO_wet_depostion',       &
          UNITS              = 'kg m-2 s-1',               &
          DIMS               = MAPL_DimsHorzOnly,          &
          VLOCATION          = MAPL_VLocationNone,         &
                                                   __RC__ )
 
       call MAPL_AddExportSpec(GC,                         &
          SHORT_NAME         = 'OCDP_OCPI',                &
          LONG_NAME          = 'OCPI_dry_depostion',       &
          UNITS              = 'kg m-2 s-1',               &
          DIMS               = MAPL_DimsHorzOnly,          &
          VLOCATION          = MAPL_VLocationNone,         &
                                                   __RC__ )
  
       call MAPL_AddExportSpec(GC,                         &
          SHORT_NAME         = 'OCDP_OCPO',                &
          LONG_NAME          = 'OCPO_dry_depostion',       &
          UNITS              = 'kg m-2 s-1',               &
          DIMS               = MAPL_DimsHorzOnly,          &
          VLOCATION          = MAPL_VLocationNone,         &
                                                   __RC__ )
  
       call MAPL_AddExportSpec(GC,                         &
          SHORT_NAME         = 'OCWT_OCPI',                &
          LONG_NAME          = 'OCPI_wet_depostion',       &
          UNITS              = 'kg m-2 s-1',               &
          DIMS               = MAPL_DimsHorzOnly,          &
          VLOCATION          = MAPL_VLocationNone,         &
                                                   __RC__ )

       call MAPL_AddExportSpec(GC,                         &
          SHORT_NAME         = 'OCWT_OCPO',                &
          LONG_NAME          = 'OCPO_wet_depostion',       &
          UNITS              = 'kg m-2 s-1',               &
          DIMS               = MAPL_DimsHorzOnly,          &
          VLOCATION          = MAPL_VLocationNone,         &
                                                   __RC__ )

       !!! to diagnose fields in AERO bundle
       call MAPL_AddExportSpec(GC,                         &
          SHORT_NAME         = 'AERO_OCphobic',            &
          LONG_NAME          = 'AERO_OCphobic',            &
          UNITS              = 'kg kg-1',                  &
          DIMS               = MAPL_DimsHorzVert,          &
          VLOCATION          = MAPL_VLocationCenter,       &
                                                   __RC__ )
  
       call MAPL_AddExportSpec(GC,                         &
          SHORT_NAME         = 'AERO_OCphilic',            &
          LONG_NAME          = 'AERO_OCphilic',            &
          UNITS              = 'kg kg-1',                  &
          DIMS               = MAPL_DimsHorzVert,          &
          VLOCATION          = MAPL_VLocationCenter,       &
                                                   __RC__ )
  
       call MAPL_AddExportSpec(GC,                         &
          SHORT_NAME         = 'AERO_BCphobic',            &
          LONG_NAME          = 'AERO_BCphobic',            &
          UNITS              = 'kg kg-1',                  &
          DIMS               = MAPL_DimsHorzVert,          &
          VLOCATION          = MAPL_VLocationCenter,       &
                                                   __RC__ )
  
       call MAPL_AddExportSpec(GC,                         &
          SHORT_NAME         = 'AERO_BCphilic',            &
          LONG_NAME          = 'AERO_BCphilic',            &
          UNITS              = 'kg kg-1',                  &
          DIMS               = MAPL_DimsHorzVert,          &
          VLOCATION          = MAPL_VLocationCenter,       &
                                                   __RC__ )
  
       call MAPL_AddExportSpec(GC,                         &
          SHORT_NAME         = 'AERO_SO4',                 &
          LONG_NAME          = 'AERO_SO4',                 &
          UNITS              = 'kg kg-1',                  &
          DIMS               = MAPL_DimsHorzVert,          &
          VLOCATION          = MAPL_VLocationCenter,       &
                                                   __RC__ )
  
       call MAPL_AddExportSpec(GC,                         &
          SHORT_NAME         = 'AERO_du001',               &
          LONG_NAME          = 'AERO_du001',               &
          UNITS              = 'kg kg-1',                  &
          DIMS               = MAPL_DimsHorzVert,          &
          VLOCATION          = MAPL_VLocationCenter,       &
                                                   __RC__ )
  
       call MAPL_AddExportSpec(GC,                         &
          SHORT_NAME         = 'AERO_du002',               &
          LONG_NAME          = 'AERO_du002',               &
          UNITS              = 'kg kg-1',                  &
          DIMS               = MAPL_DimsHorzVert,          &
          VLOCATION          = MAPL_VLocationCenter,       &
                                                   __RC__ )
  
       call MAPL_AddExportSpec(GC,                         &
          SHORT_NAME         = 'AERO_du003',               &
          LONG_NAME          = 'AERO_du003',               &
          UNITS              = 'kg kg-1',                  &
          DIMS               = MAPL_DimsHorzVert,          &
          VLOCATION          = MAPL_VLocationCenter,       &
                                                   __RC__ )
  
       call MAPL_AddExportSpec(GC,                         &
          SHORT_NAME         = 'AERO_du004',               &
          LONG_NAME          = 'AERO_du004',               &
          UNITS              = 'kg kg-1',                  &
          DIMS               = MAPL_DimsHorzVert,          &
          VLOCATION          = MAPL_VLocationCenter,       &
                                                   __RC__ )
        
       call MAPL_AddExportSpec(GC,                         &
          SHORT_NAME         = 'AERO_du005',               &
          LONG_NAME          = 'AERO_du005',               &
          UNITS              = 'kg kg-1',                  &
          DIMS               = MAPL_DimsHorzVert,          &
          VLOCATION          = MAPL_VLocationCenter,       &
                                                   __RC__ )
  
       call MAPL_AddExportSpec(GC,                         &
          SHORT_NAME         = 'AERO_ss001',               &
          LONG_NAME          = 'AERO_ss001',               &
          UNITS              = 'kg kg-1',                  &
          DIMS               = MAPL_DimsHorzVert,          &
          VLOCATION          = MAPL_VLocationCenter,       &
                                                   __RC__ )
  
       call MAPL_AddExportSpec(GC,                         &
          SHORT_NAME         = 'AERO_ss002',               &
          LONG_NAME          = 'AERO_ss002',               &
          UNITS              = 'kg kg-1',                  &
          DIMS               = MAPL_DimsHorzVert,          &
          VLOCATION          = MAPL_VLocationCenter,       &
                                                   __RC__ )
  
       call MAPL_AddExportSpec(GC,                         &
          SHORT_NAME         = 'AERO_ss003',               &
          LONG_NAME          = 'AERO_ss003',               &
          UNITS              = 'kg kg-1',                  &
          DIMS               = MAPL_DimsHorzVert,          &
          VLOCATION          = MAPL_VLocationCenter,       &
                                                   __RC__ )
  
       call MAPL_AddExportSpec(GC,                         &
          SHORT_NAME         = 'AERO_ss004',               &
          LONG_NAME          = 'AERO_ss004',               &
          UNITS              = 'kg kg-1',                  &
          DIMS               = MAPL_DimsHorzVert,          &
          VLOCATION          = MAPL_VLocationCenter,       &
                                                   __RC__ )
  
       call MAPL_AddExportSpec(GC,                         &
          SHORT_NAME         = 'AERO_ss005',               &
          LONG_NAME          = 'AERO_ss005',               &
          UNITS              = 'kg kg-1',                  &
          DIMS               = MAPL_DimsHorzVert,          &
          VLOCATION          = MAPL_VLocationCenter,       &
                                                   __RC__ )

    ENDIF ! DoAERO

    ! If GEOS-Chem is the RATS provider, we need to make sure that all 
    ! RATS quantities are available to irradiation. We will get these 
    ! quantities directly from the GEOS-Chem internal state, except for 
    ! H2O_TEND that is calculated explicitly.
    ! Since those fields are just copies of the GEOS-Chem internal
    ! species, we add them as export specs, i.e. no physics is applied
    ! to those fields.
    ! ----------------------------------------------------------------
    IF ( DoRATS ) THEN

       call MAPL_AddExportSpec(GC,                                &
          SHORT_NAME         = 'N2O',                               &
          LONG_NAME          = 'nitrous_oxide_volume_mixing_ratio', &
          UNITS              = 'mol mol-1',                         &
          DIMS               = MAPL_DimsHorzVert,                   &
          VLOCATION          = MAPL_VLocationCenter,                &
                                                            __RC__ )
  
       call MAPL_AddExportSpec(GC,                                &
          SHORT_NAME         = 'CFC11',                             &
          LONG_NAME          = 'CFC11_(CCl3F)_volume_mixing_ratio', &
          UNITS              = 'mol mol-1',                         &
          DIMS               = MAPL_DimsHorzVert,                   &
          VLOCATION          = MAPL_VLocationCenter,                &
                                                            __RC__ )
       
       call MAPL_AddExportSpec(GC,                                &
          SHORT_NAME         = 'CFC12',                             &
          LONG_NAME          = 'CFC12_(CCl2F2)_volume_mixing_ratio',&
          UNITS              = 'mol mol-1',                         &
          DIMS               = MAPL_DimsHorzVert,                   &
          VLOCATION          = MAPL_VLocationCenter,                &
                                                            __RC__ )
  
       call MAPL_AddExportSpec(GC,                                &
          SHORT_NAME         = 'HCFC22',                            &
          LONG_NAME          = 'HCFC22_(CHClF2)_volume_mixing_ratio', &
          UNITS              = 'mol mol-1',                         &
          DIMS               = MAPL_DimsHorzVert,                   &
          VLOCATION          = MAPL_VLocationCenter,                &
                                                            __RC__ )
    
       call MAPL_AddExportSpec(GC,                                &
          SHORT_NAME         = 'CH4',                               &
          LONG_NAME          = 'methane_volume_mixing_ratio',       &
          UNITS              = 'mol mol-1',                         &
          DIMS               = MAPL_DimsHorzVert,                   &
          VLOCATION          = MAPL_VLocationCenter,                &
                                                            __RC__ )
  
    ENDIF ! DoRATS

    ! Analysis ozone import
    CALL ESMF_ConfigGetAttribute( myState%myCF, I, &
             Label = "Use_ANA_O3:", Default = 0, __RC__ ) 
    IF ( I == 1 ) THEN
       CALL ESMF_ConfigGetAttribute( myState%myCF, J, &
                Label = "Use_PCHEM_O3:", Default = 0, __RC__ ) 
       IF ( J == 1 ) THEN
          call MAPL_AddImportSpec(GC,                   &
             SHORT_NAME         = 'PCHEM_O3',           &
             LONG_NAME          = 'PCHEM_ozone',        &
             UNITS              = 'kg/kg',              &
             DIMS               = MAPL_DimsHorzVert,    &
             VLOCATION          = MAPL_VLocationCenter, &
             RESTART            = MAPL_RestartSkip,     &
                                                  __RC__ )
       ELSE
          !call MAPL_AddImportSpec(GC,                   &
          !   SHORT_NAME         = 'ANA_O3',             &
          !   LONG_NAME          = 'Analysis_ozone',     &
          !   UNITS              = 'kg/kg',              &
          !   DIMS               = MAPL_DimsHorzVert,    &
          !   VLOCATION          = MAPL_VLocationCenter, &
          !   RESTART            = MAPL_RestartSkip,     &
          !                                        __RC__ )
       ENDIF
    ENDIF

    ! Diagnostics for applied ozone increment
    call MAPL_AddExportSpec(GC,                        &
       SHORT_NAME         = 'ANA_O3_INC',            &
       LONG_NAME          = 'Applied_ozone_increment', &
       UNITS              = 'kg/kg',                   &
       DIMS               = MAPL_DimsHorzVert,         &
       VLOCATION          = MAPL_VLocationCenter,      &
                                                 __RC__ )
#else
    CALL Provider_SetServices( MAPL_am_I_Root(), GC, isProvider, __RC__ )
#endif

#if defined( MODEL_GEOS )
    ! GEOS-5 only (should handle diags elsewhere):
    ! Reaction coefficients & rates 
    ! -----------------------------
    DO I = 1, NREACT
       WRITE(III,'(i3.3)') I
       ShortName = 'GCC_RC_'//TRIM(III) 
       LongName  = TRIM( ADJUSTL(EQN_NAMES(I)) )
       call MAPL_AddExportSpec(GC,                                  &
          SHORT_NAME         = TRIM(ShortName),                     &
          LONG_NAME          = TRIM(LongName),                      &
          UNITS              = 'cm3 s-1',                           &
          DIMS               = MAPL_DimsHorzVert,                   &
          VLOCATION          = MAPL_VLocationCenter,                &
                                                              __RC__ )

       ! Archive reaction coefficients for NOx diagnostics
       ! This is super hacky but not sure how else to do it without
       ! adding really elaborated syntax... 
       IF ( INDEX( EQN_NAMES(I), 'O3 + NO --> NO2 + O2' ) > 0 ) THEN
          id_rc_no  = I
       ELSEIF ( INDEX( EQN_NAMES(I), 'NO + O3 --> NO2 + O2' ) > 0 ) THEN
          id_rc_no  = I
       ELSEIF ( INDEX( EQN_NAMES(I), 'OH + NO2 --> HNO3' ) > 0 ) THEN
          id_rc_no2 = I
       ELSEIF ( INDEX( EQN_NAMES(I), 'NO2 + OH --> HNO3' ) > 0 ) THEN
          id_rc_no2 = I
       ENDIF

       ShortName = 'GCC_RR_'//TRIM(III) 
       LongName  = TRIM( ADJUSTL(EQN_NAMES(I)) )
       call MAPL_AddExportSpec(GC,                                  &
          SHORT_NAME         = TRIM(ShortName),                     &
          LONG_NAME          = TRIM(LongName),                      &
          UNITS              = 's-1',                               &
          DIMS               = MAPL_DimsHorzVert,                   &
          VLOCATION          = MAPL_VLocationCenter,                &
                                                              __RC__ )
    ENDDO

    ! Open FJX_j2j.dat to get reaction names. To be used as long names
    ! in the exports below.
    NUNIT = findFreeLUN()
    open (NUNIT,file='FJX_j2j.dat',status='old',form='formatted')
    read (NUNIT,'(a)') CLINE
    JJ = 0

    DO I = 1, JVN_
       IF ( JJ <= JVN_ ) THEN
          read (NUNIT,'(i4,1x,a50,4x,f5.3,2x,a6)') JJ,T_REACT,F_FJX,T_FJX
       ELSE
          JJ = JVN_ + 1
       ENDIF
       IF ( JJ > JVN_ ) THEN
          LongName = 'Unknown'
       ELSE
          LongName = T_REACT
       ENDIF
       CALL CSTRIP(LongName)

       WRITE(III,'(i3.3)') I
       ShortName = 'GCC_JVAL_'//TRIM(III) 
       call MAPL_AddExportSpec(GC,                                  &
          SHORT_NAME         = TRIM(ShortName),                     &
          LONG_NAME          = TRIM(LongName),                      &
          UNITS              = 's-1',                               &
          DIMS               = MAPL_DimsHorzVert,                   &
          VLOCATION          = MAPL_VLocationCenter,                &
                                                              __RC__ )
       ! Archive NO2 photolysis index for use in NOx diagnostics
       IF ( TRIM(LongName) == 'NO2PHOTONNOO' ) id_jno2 = I
    ENDDO
    close(NUNIT)
#endif

#if !defined( MODEL_GEOS )
    !=======================================================================
    !              %%% Test for archived convection fields %%%
    !=======================================================================
    CALL ESMF_ConfigGetAttribute( myState%myCF, I, &
            Label="ARCHIVED_CONV:", Default=0, __RC__ )
    ArchivedConv = ( I == 1 )

    ! Need to add archived convection fields to import state
    IF ( ArchivedConv ) THEN
       call MAPL_AddImportSpec(GC,                                  &
          SHORT_NAME         = 'ARCHIVED_PFI_CN',                   &
          LONG_NAME          = 'archived_PFI_CN',                   &
          UNITS              = 'kg m-2 s-1',                        &
          DIMS               = MAPL_DimsHorzVert,                   &
          VLOCATION          = MAPL_VLocationEdge,                  &
                                                            __RC__ )

       call MAPL_AddImportSpec(GC,                                  &
          SHORT_NAME         = 'ARCHIVED_PFL_CN',                   &
          LONG_NAME          = 'archived_PFL_CN',                   &
          UNITS              = 'kg m-2 s-1',                        &
          DIMS               = MAPL_DimsHorzVert,                   &
          VLOCATION          = MAPL_VLocationEdge,                  &
                                                            __RC__ )

       call MAPL_AddImportSpec(GC,                                  &
          SHORT_NAME         = 'ARCHIVED_CNV_MFC',                  &
          LONG_NAME          = 'archived_CNV_MFC',                  &
          UNITS              = 'kg m-2 s-1',                        &
          DIMS               = MAPL_DimsHorzVert,                   &
          VLOCATION          = MAPL_VLocationEdge,                  &
                                                            __RC__ )

       call MAPL_AddImportSpec(GC,                                  &
          SHORT_NAME         = 'ARCHIVED_CNV_MFD',                  &
          LONG_NAME          = 'archived_CNV_MFD',                  &
          UNITS              = 'kg m-2 s-1',                        &
          DIMS               = MAPL_DimsHorzVert,                   &
          VLOCATION          = MAPL_VLocationCenter,                &
                                                            __RC__ )

       call MAPL_AddImportSpec(GC,                                  &
          SHORT_NAME         = 'ARCHIVED_CNV_CVW',                  &
          LONG_NAME          = 'archived_CNV_CVW',                  &
          UNITS              = 'hPa s-1',                           &
          DIMS               = MAPL_DimsHorzVert,                   &
          VLOCATION          = MAPL_VLocationCenter,                &
                                                            __RC__ )

       call MAPL_AddImportSpec(GC,                                  &
          SHORT_NAME         = 'ARCHIVED_DQRC',                     &
          LONG_NAME          = 'archived_DQRC',                     &
          UNITS              = 'kg kg-1 s-1',                       &
          DIMS               = MAPL_DimsHorzVert,                   &
          VLOCATION          = MAPL_VLocationCenter,                &
                                                            __RC__ )

       call MAPL_AddImportSpec(GC,                                  &
          SHORT_NAME         = 'ARCHIVED_REV_CN',                   &
          LONG_NAME          = 'archived_REV_CN',                   &
          UNITS              = 'kg kg-1 s-1',                       &
          DIMS               = MAPL_DimsHorzVert,                   &
          VLOCATION          = MAPL_VLocationCenter,                &
                                                            __RC__ )

       call MAPL_AddImportSpec(GC,                                  &
          SHORT_NAME         = 'ARCHIVED_T',                        &
          LONG_NAME          = 'archived_T',                        &
          UNITS              = 'K',                                 &
          DIMS               = MAPL_DimsHorzVert,                   &
          VLOCATION          = MAPL_VLocationCenter,                &
                                                            __RC__ )
    ENDIF ! ArchivedConv 

#endif

    ! OLSON
    DO T = 1, NSURFTYPE
       landTypeInt = T-1
       WRITE ( landTypeStr, '(I2.2)' ) landTypeInt
       importName = 'OLSON' // TRIM(landTypeStr)
       CALL MAPL_AddImportSpec(GC,                                  &
          SHORT_NAME         = importName,                          &
          LONG_NAME          = 'OLSON_land_by_type',                &
          UNITS              = 'unitless',                          &
          DIMS               = MAPL_DimsHorzOnly,                   &
                                                            __RC__ )
    ENDDO

    ! Set HEMCO services
    ! --------------------
    CALL ESMF_ConfigGetAttribute( myState%myCF, HcoConfigFile, &
                                  Label="HEMCO_CONFIG:", &
                                  Default="HEMCO_Config.rc", __RC__ )
    CALL HCO_SetServices( MAPL_am_I_Root(), GC, HcoConfig,  &
                          TRIM(HcoConfigFile), __RC__ )    

    ! Set the Profiling timers
    ! ------------------------
    CALL MAPL_TimerAdd(GC, NAME="INITIALIZE", RC=status)
    _VERIFY(status)
    CALL MAPL_TimerAdd(GC, NAME="RUN", RC=status)
    _VERIFY(status)
    CALL MAPL_TimerAdd(GC, NAME="FINALIZE", RC=status)
    _VERIFY(status)

    CALL MAPL_TimerAdd(GC, NAME="DO_CHEM", RC=status)
    _VERIFY(status)
    CALL MAPL_TimerAdd(GC, NAME="CP_BFRE", RC=status)
    _VERIFY(status)
    CALL MAPL_TimerAdd(GC, NAME="CP_AFTR", RC=status)
    _VERIFY(status)

    ! More timers to be called in gigc_chunk_run 
    CALL MAPL_TimerAdd(GC, NAME="GC_CONV"  , __RC__)
    CALL MAPL_TimerAdd(GC, NAME="GC_EMIS"  , __RC__)
    CALL MAPL_TimerAdd(GC, NAME="GC_DRYDEP", __RC__)
    CALL MAPL_TimerAdd(GC, NAME="GC_FLUXES", __RC__)
    CALL MAPL_TimerAdd(GC, NAME="GC_TURB"  , __RC__)
    CALL MAPL_TimerAdd(GC, NAME="GC_CHEM"  , __RC__)
    CALL MAPL_TimerAdd(GC, NAME="GC_WETDEP", __RC__)
    CALL MAPL_TimerAdd(GC, NAME="GC_DIAGN" , __RC__)

    ! Generic Set Services
    ! --------------------
    CALL MAPL_GenericSetServices( GC, RC=status )
    _VERIFY(status)

    !=======================================================================
    ! All done
    !=======================================================================
    _RETURN(ESMF_SUCCESS)

  END SUBROUTINE SetServices
!EOC
!------------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 910.1 and      !
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Initialize_ 
!
! !DESCRIPTION: Initialize_ is the initialize method of the GEOSCHEMchem
!  gridded component.  This is a simple ESMF/MAPL wrapper which calls down
!  to the Initialize method of the GEOS-Chem column chemistry code.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Initialize_( GC, Import, Export, Clock, RC )
!
! !USES:
!
    USE TIME_MOD,  ONLY : GET_TS_CHEM, GET_TS_EMIS
    USE TIME_MOD,  ONLY : GET_TS_DYN,  GET_TS_CONV
#if defined( MODEL_GEOS )
    USE TENDENCIES_MOD, ONLY : Tend_CreateClass
    USE TENDENCIES_MOD, ONLY : Tend_Add
#endif
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ESMF_GridComp), INTENT(INOUT)         :: GC       ! Ref to GridComp
    TYPE(ESMF_State),    INTENT(INOUT), TARGET :: Import   ! Import State object
    TYPE(ESMF_State),    INTENT(INOUT), TARGET :: Export   ! Export State object
    TYPE(ESMF_Clock),    INTENT(INOUT)         :: Clock    ! ESMF clock object
!                                                  
! !OUTPUT PARAMETERS:                              
!                                                  
    INTEGER,             INTENT(OUT)           :: RC       ! Success or failure?
!
! !REMARKS:
!  We call routine Extract_ to return various values (i.e. grid parameters,
!  start & end dates, PET information, etc.) from the ESMF/MAPL environment.  
!  We then pass those to GEOS-Chem via routine GIGC_CHUNK_INIT, which is
!  located in GEOS-Chem module ./GEOS-Chem/ESMF/gigc_chunk_mod.F90.
!
! !REVISION HISTORY:
!  06 Dec 2009 - A. da Silva - Initial version
!  08 Apr 2010 - R. Yantosca - Now uses the updated Extract_ method.
!  14 Apr 2010 - R. Yantosca - Activated call to GC_CHUNK_INIT
!  15 Apr 2010 - R. Yantosca - Add extra error checks for dimensions
!  23 Apr 2010 - R. Yantosca - Now pass IDENT obj to GC_CHUNK_INIT routine
!  30 Apr 2010 - R. Yantosca - Now use 5 digits for PET
!  02 Jun 2010 - R. Yantosca - Now set Ident%VERBOSE to FALSE
!  09 Oct 2012 - R. Yantosca - Now call MAPL_Am_I_Root to test for root PET
!  08 Nov 2012 - R. Yantosca - Now pass options to G-C via Input_Opt object
!  04 Dec 2012 - R. Yantosca - Now get local PET grid indices as well as the
!                              global grid indices from the Extract_ function
!  26 Feb 2013 - R. Yantosca - Now read Input_Opt%MAX_DEP from the rc file
!  08 Mar 2013 - R. Yantosca - Now save the PET # (aka PET #) in Input_Opt
!  15 Mar 2013 - R. Yantosca - Remove IDENT object, which was a holdover from
!                              the GEOS-Chem column code
!  13 Oct 2014 - C. Keller   - Updated for HEMCO
!  24 Oct 2014 - C. Keller   - Updated for RATS/AERO/Analysis OX provider
!  23 Feb 2015 - C. Keller   - Now use local variable haveImpRst
!  06 Nov 2017 - E. Lundgren - Abstract provider services to new module
!                              gigc_providerservices_mod.F90
!  06 Mar 2018 - E. Lundgren - Remove obsolete variables; update usage of GC
!                              timesteps since now in seconds not minutes
!EOP
!------------------------------------------------------------------------------
!BOC
!
! LOCAL VARIABLES:
!
    ! Objects
    TYPE(ESMF_Grid)             :: Grid        ! ESMF Grid object
    TYPE(ESMF_Config)           :: MaplCF      ! ESMF Config obj (MAPL.rc)
    TYPE(ESMF_Config)           :: GeosCF      ! ESMF Config obj (GEOSCHEM*.rc) 
                                       
    ! Scalars                                   
    LOGICAL                     :: am_I_Root   ! Are we on the root PET?
    INTEGER                     :: myPet       ! # of the PET we are on 
    INTEGER                     :: NPES        ! # of total PETs in MPI world 
    INTEGER                     :: nymdB       ! GMT date @ start of simulation
    INTEGER                     :: nymdE       ! GMT date @ end of simulation
    INTEGER                     :: nymd        ! GMT date @ current time
    INTEGER                     :: nhmsB       ! GMT time @ start of simulation
    INTEGER                     :: nhmsE       ! GMT time @ end of simulation
    INTEGER                     :: nhms        ! GMT time @ current time
    INTEGER                     :: IM          ! # of longitudes on this PET
    INTEGER                     :: JM          ! # of latitudes  on this PET
    INTEGER                     :: LM          ! # of levels     on this PET
    INTEGER                     :: value_LLSTRAT ! # of strat. levels  
    INTEGER                     :: IM_WORLD    ! # of longitudes in global grid
    INTEGER                     :: JM_WORLD    ! # of latitudes  in global grid
    INTEGER                     :: LM_WORLD    ! # of levels     in global grid
    REAL                        :: tsChem      ! Chemistry timestep [s]
    REAL                        :: tsDyn       ! Dynamic timestep [s]
    CHARACTER(LEN=5)            :: petStr      ! String for PET #
    CHARACTER(LEN=ESMF_MAXSTR)  :: compName    ! Name of gridded component
    
    ! time step error checks 
    REAL                         :: ChemTS, EmisTS

    ! Pointer arrays
    REAL(ESMF_KIND_R4),  POINTER :: lonCtr(:,:) ! Lon centers on this PET [rad]
    REAL(ESMF_KIND_R4),  POINTER :: latCtr(:,:) ! Lat centers on this PET [rad]

    INTEGER                      :: I, nFlds, mpiComm
    TYPE(ESMF_STATE)             :: INTSTATE 
    TYPE(ESMF_Field)             :: GcFld

#if defined( MODEL_GEOS )
    ! Does GEOS-Chem restart file exist?
    ! Before broadcasting, we check if there is an import restart file for
    ! GEOS-Chem. This variable is then passed to Input_Opt after 
    ! initialization (and CPU broadcasting) of all GEOS-Chem variables.
    LOGICAL                      :: haveImpRst

    TYPE(MAPL_MetaComp), POINTER :: STATE

    ! For AERO
    INTEGER                      :: J, GCID
    LOGICAL                      :: FRIENDLY
    REAL                         :: GCMW, FRAC
    REAL, POINTER                :: Ptr3D(:,:,:) => NULL()
    TYPE(ESMF_STATE)             :: Aero
    TYPE(ESMF_FieldBundle)       :: AeroBdl 
    TYPE(ESMF_Field)             :: AeroFld
    CHARACTER(LEN=ESMF_MAXSTR)   :: GCName, AeroName, fieldName

    ! v11: species
    TYPE(Species), POINTER       :: SpcInfo

    ! To read various options 
    INTEGER                      :: DoIt 
    REAL                         :: Val, OzPause

    ! Mie table updates
    INTEGER                      :: instance
#else
    INTEGER                      :: N, trcID
    TYPE(ESMF_Time)              :: CurrTime    ! Current time of the ESMF clock
    TYPE(MAPL_MetaComp), POINTER :: STATE => NULL()
    REAL(ESMF_KIND_R8), POINTER  :: Ptr3D_R8(:,:,:) => NULL()
#endif

    __Iam__('Initialize_')

    !=======================================================================
    ! Initialization
    !=======================================================================

    ! Get my name and set-up traceback handle
    CALL ESMF_GridCompGet( GC, name=compName, __RC__ )

    ! Identify this routine to MAPL
    Iam = TRIM(compName)//'::Initialize_'
    
    ! Get my MAPL_Generic state
    ! -------------------------
    CALL MAPL_GetObjectFromGC(GC, STATE, RC=STATUS)
    _VERIFY(STATUS)

    !  Start timers
    !  ------------
    CALL MAPL_TimerOn( STATE, "TOTAL")
    CALL MAPL_TimerOn( STATE, "INITIALIZE")

    ! Initialize MAPL Generic
    CALL MAPL_GenericInitialize( GC, Import, Export, Clock, __RC__ )

    ! Get Internal state.
    CALL MAPL_Get ( STATE, INTERNAL_ESMF_STATE=INTSTATE, __RC__ ) 

    ! Initialize GEOS-Chem Input_Opt fields to zeros or equivalent
    CALL Set_Input_Opt( MAPL_am_I_Root(), Input_Opt, RC )
    _ASSERT(RC==GC_SUCCESS, 'Error calling Set_Input_Opt')

    ! Get various parameters from the ESMF/MAPL framework
#if defined( MODEL_GEOS )
    ! Note: variable haveImpRst must not yet be written into Input_Opt yet 
    ! since the variables of Input_Opt may be 'erased' during initialization
    ! of GEOS-Chem.
#endif
    CALL Extract_( GC,                        &  ! Ref to this Gridded Comp 
                   Clock,                     &  ! ESMF Clock object
                   Grid        = Grid,        &  ! ESMF Grid object
                   MaplCF      = MaplCF,      &  ! AGCM.rc/GCHP.rc config object
                   GeosCF      = GeosCF,      &  ! GEOSCHEM*.rc Config object
                   IM          = IM,          &  ! # of longitudes on this PET
                   JM          = JM,          &  ! # of latitudes  on this PET
                   LM          = LM,          &  ! # of levels     on this PET
                   IM_WORLD    = IM_WORLD,    &  ! # of lons in global grid
                   JM_WORLD    = JM_WORLD,    &  ! # of lats  in global grid
                   LM_WORLD    = LM_WORLD,    &  ! # of levels in global grid
                   nymdB       = nymdB,       &  ! YYYYMMDD @ start of sim
                   nhmsB       = nhmsB,       &  ! hhmmss   @ end   of sim
                   nymdE       = nymdE,       &  ! YYYMMDD  @ start of sim
                   nhmsE       = nhmsE,       &  ! hhmmss   @ end   of sim
                   tsChem      = tsChem,      &  ! Chemistry timestep [seconds]
                   tsDyn       = tsDyn,       &  ! Dynamics timestep  [seconds]
                   localPet    = myPet,       &  ! PET # that we are on now
                   petCount    = NPES,        &  ! Number of PETs in MPI World
                   mpiComm     = mpiComm,     &  ! MPI Communicator Handle
                   lonCtr      = lonCtr,      &  ! This PET's lon ctrs [radians]
                   latCtr      = latCtr,      &  ! This PET's lat ctrs [radians]
#if defined( MODEL_GEOS )
		   haveImpRst  = haveImpRst,  &  ! Does import restart exist? 
#endif
                   __RC__                      )

    ! Set MPI values in Input_Opt
    Input_Opt%thisCPU = myPet
    Input_Opt%MPIComm = mpiComm
    Input_Opt%numCPUs = NPES
    Input_Opt%isMPI   = .true.
    if ( MAPL_am_I_Root() ) Input_Opt%amIRoot = .true.

#if defined( MODEL_GEOS )
    Input_Opt%haveImpRst = haveImpRst
#endif

    ! MSL - shift from 0 - 360 to -180 - 180 degree grid
    where (lonCtr .gt. MAPL_PI ) lonCtr = lonCtr - 2*MAPL_PI

#if !defined( MODEL_GEOS )
    ! Get the memory debug level
    call ESMF_ConfigGetAttribute(GeosCF, MemDebugLevel, &
                                 Label="MEMORY_DEBUG_LEVEL:" , RC=STATUS)
    _VERIFY(STATUS)
#endif

    !=======================================================================
    ! Save values from the resource file (GCHP.rc for GCHP)
    !=======================================================================

    ! # of run phases
    CALL ESMF_ConfigGetAttribute( GeosCF, NPHASE,                   & 
                                  Default = 2,                      &
                                  Label   = "RUN_PHASES:",          &
                                  __RC__                           )
    _ASSERT(NPHASE==1.OR.NPHASE==2,'Error calling ESMF_ConfigGetAttribute on RUN_PHASES') 

#if defined( MODEL_GEOS )
    ! Top stratospheric level
    CALL ESMF_ConfigGetAttribute( GeosCF, value_LLSTRAT,            & 
                                  Default = LM,                     &
                                  Label   = "LLSTRAT:",             &
                                  __RC__                           )
    IF ( Input_Opt%AmIRoot ) THEN
       WRITE(*,*) 'GCC: top strat. level (LLSTRAT) is set to ',     &
                  value_LLSTRAT
    ENDIF

    ! FAST-JX settings: number of levels, number of EXTRAL iterations,
    ! print error if EXTRAL fails? 
    ! LLFASTJX: default is 1201 for LM=132, 601 otherwise
    IF ( LM == 132 ) THEN
       I = 1201
    ELSE
       I = 601
    ENDIF 
    CALL ESMF_ConfigGetAttribute( GeosCF, Input_Opt%LLFASTJX,     & 
                                  Default = I,                    &
                                  Label   = "LLFASTJX:",          &
                                  __RC__                          )

    ! FJX_EXTRAL_ITERMAX: default is 5 for LM=132, 1 otherwise
    IF ( LM == 132 ) THEN
       I = 5 
    ELSE
       I = 1 
    ENDIF 
    CALL ESMF_ConfigGetAttribute( GeosCF, Input_Opt%FJX_EXTRAL_ITERMAX, & 
                                  Default = I,                          &
                                  Label   = "FJX_EXTRAL_ITERMAX:",      &
                                  __RC__                                )

    ! FJX_EXTRAL_ERR: default is 1
    I = 1
    CALL ESMF_ConfigGetAttribute( GeosCF, DoIt,   & 
                                  Default = I,                        &
                                  Label   = "FJX_EXTRAL_ERR:",        &
                                  __RC__                              )
    Input_Opt%FJX_EXTRAL_ERR = ( DoIt == 1 )
    IF ( Input_Opt%AmIRoot ) THEN
       WRITE(*,*) 'GCC: Fast-JX settings:'
       WRITE(*,*) 'Number of FAST-JX levels      : ',Input_Opt%LLFASTJX
       WRITE(*,*) 'Max. no. of EXTRAL iterations : ',Input_Opt%FJX_EXTRAL_ITERMAX
       WRITE(*,*) 'Show EXTRAL overflow error    : ',Input_Opt%FJX_EXTRAL_ERR
    ENDIF

    ! Stop KPP if integration fails twice
    CALL ESMF_ConfigGetAttribute( GeosCF, DoIt,   & 
                                  Default = 1,                       &
                                  Label   = "KPP_STOP_IF_FAIL:",     &
                                  __RC__                              )
    Input_Opt%KppStop = ( DoIt == 1 )
    IF ( Input_Opt%AmIRoot ) THEN
       WRITE(*,*) 'Stop KPP if integration fails: ',Input_Opt%KppStop
    ENDIF

    ! Turn off three heterogenous reactions in stratosphere 
    CALL ESMF_ConfigGetAttribute( GeosCF, DoIt, Default = 0, &
                                  Label = "TurnOffHetRates:", __RC__ )
    Input_Opt%TurnOffHetRates = ( DoIt == 1 )
    IF ( Input_Opt%AmIRoot ) THEN
       WRITE(*,*) 'Disable selected het. reactions in stratosphere: ', &
                  Input_Opt%TurnOffHetRates
    ENDIF
#endif

    !=======================================================================
    ! Initialize GEOS-Chem (will also initialize HEMCO)
    !=======================================================================

    ! Initialize fields of the Grid State object
    CALL Init_State_Grid( Input_Opt, State_Grid, RC )
    _ASSERT(RC==GC_SUCCESS,'Error calling Init_State_Grid')
  
    ! Pass grid information obtained from Extract_ to State_Grid
    State_Grid%NX          = IM            ! # lons   on this PET
    State_Grid%NY          = JM            ! # lats   on this PET
    State_Grid%NZ          = LM            ! # levels on this PET
    State_Grid%GlobalNX    = IM_WORLD      ! # lons   in global grid
    State_Grid%GlobalNY    = JM_WORLD      ! # lats   in global grid
    State_Grid%NativeNZ    = LM_WORLD      ! # levels in global grid
    State_Grid%XMinOffset  = 1             ! X offset from global grid
    State_Grid%XMaxOffset  = State_Grid%NX ! X offset from global grid
    State_Grid%YMinOffset  = 1             ! Y offset from global grid
    State_Grid%YMaxOffset  = State_Grid%NY ! Y offset from global grid
    State_Grid%MaxTropLev  = 40            ! # trop. levels
#if defined( MODEL_GEOS )
    State_Grid%MaxStratLev = value_LLSTRAT ! # strat. levels
#else
    State_Grid%MaxStratLev = 59            ! # strat. levels
#endif

    ! Call the GIGC initialize routine
    CALL GIGC_Chunk_Init( nymdB     = nymdB,      & ! YYYYMMDD @ start of run
                          nhmsB     = nhmsB,      & ! hhmmss   @ start of run
                          nymdE     = nymdE,      & ! YYYYMMDD @ end of run
                          nhmsE     = nhmsE,      & ! hhmmss   @ end of run
                          tsChem    = tsChem,     & ! Chemical timestep [s]
                          tsDyn     = tsDyn,      & ! Dynamic  timestep [s]
                          lonCtr    = lonCtr,     & ! Lon centers [radians]
                          latCtr    = latCtr,     & ! Lat centers [radians]
#if !defined( MODEL_GEOS )
                          GC        = GC,         & ! Ref to this gridded comp
                          EXPORT    = EXPORT,     & ! Export state object
#endif
                          Input_Opt = Input_Opt,  & ! Input Options obj
                          State_Chm = State_Chm,  & ! Chemistry State obj
                          State_Diag= State_Diag, & ! Diagnostics State obj
                          State_Grid= State_Grid, & ! Grid State obj
                          State_Met = State_Met,  & ! Meteorology State obj
                          HcoConfig = HcoConfig,  & ! HEMCO config obj 
                          HistoryConfig = HistoryConfig, & ! History Config Obj
                          __RC__                 )

#if defined( MODEL_GEOS )
    !=======================================================================
    ! If GEOS-Chem is the AERO provider, initialize the AERO bundle here.
    ! All GEOS-Chem tracers possibly being added to the AERO bundle are
    ! listed at the beginning of the module. Here, we see which ones of 
    ! those are effectively defined and create a field in the bundle for
    ! them. The AERO names are given the names listed at the beginning of
    ! the module.
    ! GEOS-Chem tracers are in mol/mol, whereas the AERO bundle holds
    ! data in kg/kg. We therefore need to copy the data so that we can 
    ! change units independently.
    !=======================================================================
    IF ( DoAERO ) THEN

       ! Get AERO bundle
       CALL ESMF_StateGet( EXPORT, 'AERO', Aero, __RC__ )

       ! This attribute indicates if the aerosol optics method is implemented 
       ! or not. Radiation will not call the aerosol optics method unless this 
       ! attribute is explicitly set to true.
       call ESMF_AttributeSet(aero, name='implements_aerosol_optics_method', &
                              value=.true., __RC__)
       aeroBdl = ESMF_FieldBundleCreate(name='AEROSOLS', __RC__)
       call MAPL_StateAdd(Aero, aeroBdl, __RC__)

       ! Loop over all GC tracers that we may want to add to the AERO
       ! bundle
       DO I = 1, NumAERO

          ! Get GEOS-Chem tracer ID
          GCID = Ind_( TRIM(GcNames(I)) )

          ! If species is defined, copy field and add to AERO bundle
          IF ( GCID > 0 ) THEN

             ! This is the name in the internal state
             GCName = TRIM(SPFX) // TRIM(GcNames(I))

             ! Get field from internal state
             CALL ESMF_StateGet( INTSTATE, TRIM(GCName), GcFld, RC=RC )

             ! Try TRC_<NAME> if SPC_<NAME> not found
             IF ( RC /= ESMF_SUCCESS ) THEN
                GCName = 'TRC_'//TRIM(GcNames(I))
                CALL ESMF_StateGet( INTSTATE, TRIM(GCName), GcFld, RC=RC )
             ENDIF

             ! Error if none of the above found
             IF ( RC /= ESMF_SUCCESS ) THEN
                WRITE(*,*) 'Cannot fill AERO bundle - field not found in ' // &
                           'internal state: ' // TRIM(GCName)
                _ASSERT(.FALSE.,'Error filling AERO bundle')
             ENDIF
  
             ! Set number of fields to be created. This is only different from
             ! 1 for sea salt aerosols, which are mapped onto multiple AERO
             ! fields.
             NFLDS = 1
             IF ( TRIM(GcNames(I)) == 'SALA' ) NFLDS = 2
             IF ( TRIM(GcNames(I)) == 'SALC' ) NFLDS = 3
             IF ( TRIM(GcNames(I)) == 'DST4' ) NFLDS = 2

             ! Now create all fields
             DO J = 1, NFLDS
 
                ! AERO field name
                AeroName = TRIM(AeroNames(I))
                IF ( TRIM(GcNames(I)) == 'SALA' ) AeroName = SALAnames(J)
                IF ( TRIM(GcNames(I)) == 'SALC' ) AeroName = SALCnames(J)
                IF ( TRIM(GcNames(I)) == 'DST4' ) AeroName = DST4names(J)
 
                ! Create new field
                AeroFld = MAPL_FieldCreate( GcFld, name=AeroName, &
                                            DoCopy=.TRUE., __RC__  )
      
                ! Get molecular weight (g/mol)
                !GCMW = Input_Opt%Tracer_MW_G(GCID)
                GCMW = GocartMW(I)      

                ! Fraction of the GC field to be used in the AERO field
                FRAC = 1.0
                IF ( TRIM(GcNames(I)) == 'SALA' ) FRAC = SALAsplit(J)
                IF ( TRIM(GcNames(I)) == 'SALC' ) FRAC = SALCsplit(J)
                IF ( TRIM(GcNames(I)) == 'DST4' ) FRAC = DST4split(J)

                ! Pass GEOS-Chem field name, molecular weight and fraction 
                ! to be used to bundle for easier handling lateron
                CALL ESMF_AttributeSet ( AeroFld, NAME='GCNAME', &
                                         VALUE=GCName, __RC__ ) 
                CALL ESMF_AttributeSet ( AeroFld, NAME='GCMW',   &
                                         VALUE=GCMW,   __RC__ ) 
                CALL ESMF_AttributeSet ( AeroFld, NAME='FRAC',   &
                                         VALUE=FRAC,   __RC__ ) 
      
                ! Before adding to the bundle, convert data from mol/mol to 
                ! kg/kg. Data is now stored in kg/kg total. (ckeller, 3/7/17)
                CALL ESMF_FieldGet( AeroFld, farrayPtr=Ptr3D, __RC__ )
                !Ptr3D = Ptr3D * GCMW / MAPL_AIRMW * FRAC
                Ptr3D = Ptr3D * FRAC
                Ptr3D => NULL()
   
                ! Add to bundle
                CALL MAPL_FieldBundleAdd ( AeroBdl, AeroFld, __RC__ )
             ENDDO !J
          ENDIF
       ENDDO

       ! Mie table
       instance = instanceComputational
       geoschemMieTable(instance) = Chem_MieCreate(MaplCF, __RC__)
       call ESMF_AttributeSet(aero, name='mie_table_instance', &
                              value=instance, __RC__)

       ! state of the atmosphere
       call ESMF_AttributeSet(aero, name='air_pressure_for_aerosol_optics',             value='PLE', __RC__)
       call ESMF_AttributeSet(aero,   &
                              name='relative_humidity_for_aerosol_optics',  &
                              value='RH',  __RC__)
       ! 'cloud_area_fraction_in_atmosphere_layer_for_aerosol_optics'
       call ESMF_AttributeSet(aero,   &
                              name='cloud_area_fraction_for_aerosol_optics', &
                              value='',    __RC__) 

       ! aerosol optics
       call ESMF_AttributeSet(aero, name='band_for_aerosol_optics',                     value=0,     __RC__)
       call ESMF_AttributeSet(aero,  &
                     name='extinction_in_air_due_to_ambient_aerosol', &
                     value='EXT', __RC__)
       call ESMF_AttributeSet(aero,  &
                     name='single_scattering_albedo_of_ambient_aerosol', &
                     value='SSA', __RC__)
       call ESMF_AttributeSet(aero,  &
                     name='asymmetry_parameter_of_ambient_aerosol', &
                     value='ASY', __RC__)

       ! add PLE to Aero state
       call ESMF_AttributeGet(aero,  &
                     name='air_pressure_for_aerosol_optics',  &
                     value=fieldName, __RC__)
       if (fieldName /= '') then
          aeroFld = MAPL_FieldCreateEmpty(trim(fieldName), grid, __RC__)

          call MAPL_FieldAllocCommit(aeroFld, dims=MAPL_DimsHorzVert, &
                                     location=MAPL_VLocationEdge,     &
                                     typekind=MAPL_R4, hw=0, __RC__)
          call MAPL_StateAdd(aero, aeroFld, __RC__)
       end if

       ! add RH to Aero state
       call ESMF_AttributeGet(aero,  &
                      name='relative_humidity_for_aerosol_optics', &
                      value=fieldName, __RC__)
       if (fieldName /= '') then
          aeroFld = MAPL_FieldCreateEmpty(trim(fieldName), grid, __RC__)

          call MAPL_FieldAllocCommit(aeroFld, dims=MAPL_DimsHorzVert, &
                                     location=MAPL_VLocationCenter,   &
                                     typekind=MAPL_R4, hw=0, __RC__)
          call MAPL_StateAdd(aero, aeroFld, __RC__)
       end if

       ! add EXT to Aero state
       call ESMF_AttributeGet(aero,  &
                            name='extinction_in_air_due_to_ambient_aerosol', &
                            value=fieldName, __RC__)
       if (fieldName /= '') then
          aeroFld = MAPL_FieldCreateEmpty(trim(fieldName), grid, __RC__)

          call MAPL_FieldAllocCommit(aeroFld, dims=MAPL_DimsHorzVert, &
                                     location=MAPL_VLocationCenter,   &
                                     typekind=MAPL_R4, hw=0, __RC__)
          call MAPL_StateAdd(aero, aeroFld, __RC__)
       end if

       ! add SSA to aero state
       call ESMF_AttributeGet(aero,  &
                         name='single_scattering_albedo_of_ambient_aerosol',  &
                         value=fieldName, __RC__)
       if (fieldName /= '') then
          aeroFld = MAPL_FieldCreateEmpty(trim(fieldName), grid, __RC__)

          call MAPL_FieldAllocCommit(aeroFld, dims=MAPL_DimsHorzVert, &
                                     location=MAPL_VLocationCenter,   &
                                     typekind=MAPL_R4, hw=0, __RC__)
          call MAPL_StateAdd(aero, aeroFld, __RC__)
       end if

       ! add ASY to aero state
       call ESMF_AttributeGet(aero,   &
                              name='asymmetry_parameter_of_ambient_aerosol', &
                              value=fieldName, __RC__)
       if (fieldName /= '') then
          aeroFld = MAPL_FieldCreateEmpty(trim(fieldName), grid, __RC__)

          call MAPL_FieldAllocCommit(aeroFld, dims=MAPL_DimsHorzVert, &
                                     location=MAPL_VLocationCenter,   &
                                     typekind=MAPL_R4, hw=0, __RC__)
          call MAPL_StateAdd(aero, aeroFld, __RC__)
       end if

       ! attach the aerosol optics method
       call ESMF_MethodAdd(aero, label='aerosol_optics', &
                           userRoutine=aerosol_optics, __RC__)

       ! ---------------------------------------------------------------------
       ! Initialize the AERO_DP bundle
       ! ---------------------------------------------------------------------
       CALL ESMF_StateGet( EXPORT, 'AERO_DP', AeroBdl, __RC__ )

       ! Dust dry and wet deposition 
       CALL ESMF_StateGet( EXPORT, 'DUDP_DST1', AeroFld, __RC__ )
       CALL MAPL_FieldBundleAdd( AeroBdl, AeroFld, __RC__ )

       CALL ESMF_StateGet( EXPORT, 'DUDP_DST2', AeroFld, __RC__ )
       CALL MAPL_FieldBundleAdd( AeroBdl, AeroFld, __RC__ )

       CALL ESMF_StateGet( EXPORT, 'DUDP_DST3', AeroFld, __RC__ )
       CALL MAPL_FieldBundleAdd( AeroBdl, AeroFld, __RC__ )

       CALL ESMF_StateGet( EXPORT, 'DUDP_DST4', AeroFld, __RC__ )
       CALL MAPL_FieldBundleAdd( AeroBdl, AeroFld, __RC__ )
      
       CALL ESMF_StateGet( EXPORT, 'DUWT_DST1', AeroFld, __RC__ )
       CALL MAPL_FieldBundleAdd( AeroBdl, AeroFld, __RC__ )

       CALL ESMF_StateGet( EXPORT, 'DUWT_DST2', AeroFld, __RC__ )
       CALL MAPL_FieldBundleAdd( AeroBdl, AeroFld, __RC__ )

       CALL ESMF_StateGet( EXPORT, 'DUWT_DST3', AeroFld, __RC__ )
       CALL MAPL_FieldBundleAdd( AeroBdl, AeroFld, __RC__ )

       CALL ESMF_StateGet( EXPORT, 'DUWT_DST4', AeroFld, __RC__ )
       CALL MAPL_FieldBundleAdd( AeroBdl, AeroFld, __RC__ )

       ! Black carbon dry and wet depostion 
       CALL ESMF_StateGet( EXPORT, 'BCDP_BCPI', AeroFld, __RC__ )
       CALL MAPL_FieldBundleAdd( AeroBdl, AeroFld, __RC__ )
      
       CALL ESMF_StateGet( EXPORT, 'BCDP_BCPO', AeroFld, __RC__ )
       CALL MAPL_FieldBundleAdd( AeroBdl, AeroFld, __RC__ )
      
       CALL ESMF_StateGet( EXPORT, 'BCWT_BCPI', AeroFld, __RC__ )
       CALL MAPL_FieldBundleAdd( AeroBdl, AeroFld, __RC__ )

       CALL ESMF_StateGet( EXPORT, 'BCWT_BCPO', AeroFld, __RC__ )
       CALL MAPL_FieldBundleAdd( AeroBdl, AeroFld, __RC__ )

       ! Organic carbon dry and wet depostion 
       CALL ESMF_StateGet( EXPORT, 'OCDP_OCPI', AeroFld, __RC__ )
       CALL MAPL_FieldBundleAdd( AeroBdl, AeroFld, __RC__ )
      
       CALL ESMF_StateGet( EXPORT, 'OCDP_OCPO', AeroFld, __RC__ )
       CALL MAPL_FieldBundleAdd( AeroBdl, AeroFld, __RC__ )
      
       CALL ESMF_StateGet( EXPORT, 'OCWT_OCPI', AeroFld, __RC__ )
       CALL MAPL_FieldBundleAdd( AeroBdl, AeroFld, __RC__ )

       CALL ESMF_StateGet( EXPORT, 'OCWT_OCPO', AeroFld, __RC__ )
       CALL MAPL_FieldBundleAdd( AeroBdl, AeroFld, __RC__ )

    ENDIF ! DoAERO

    ! ckeller, 6/12/19: now do during run stage
    !CALL MAPL_GetPointer( INTSTATE, PTR_O3, 'TRC_O3', &
    !                      NotFoundOk=.TRUE., __RC__ )
    !
    !IF ( DoRATS ) THEN
    !   CALL MAPL_GetPointer( INTSTATE,    PTR_CH4, 'TRC_CH4',    __RC__ )
    !   CALL MAPL_GetPointer( INTSTATE,    PTR_N2O, 'TRC_N2O',    __RC__ )
    !   CALL MAPL_GetPointer( INTSTATE,  PTR_CFC11, 'TRC_CFC11',  __RC__ )
    !   CALL MAPL_GetPointer( INTSTATE,  PTR_CFC12, 'TRC_CFC12',  __RC__ )
    !   CALL MAPL_GetPointer( INTSTATE, PTR_HCFC22, 'TRC_HCFC22', __RC__ )
    !   CALL MAPL_GetPointer( INTSTATE,    PTR_H2O, 'TRC_H2O',    __RC__ )
    !ENDIF
#else
    IF ( isProvider ) THEN
       CALL Provider_Initialize( am_I_Root, State_Chm, State_Grid, &
                                 INTSTATE,  EXPORT,    __RC__ )
    ENDIF
#endif

    !=======================================================================
    ! Initialize the Int2Spc object. This is used to copy the tracer arrays
    ! from the internal state to State_Chm%Tracers, and vice versa.
    ! In this step, we also check for the friendlieness of the tracers. If
    ! the GEOS-Chem internal convection/turbulence schemes shall be used
    ! (as specified in input.geos), the tracers must not be friendly to
    ! the GEOS-5 moist / turbulence components!
    !=======================================================================
#if defined( MODEL_GEOS )
    nFlds = State_Chm%nSpecies
#else
    ! GCHP uses nAdvect (bug?):
    nFlds = State_Chm%nAdvect
#endif
    ALLOCATE( Int2Spc(nFlds), STAT=STATUS )
    _ASSERT(STATUS==0,'Int2Spc could not be allocated')

    ! Do for every tracer in State_Chm
    DO I = 1, nFlds

#if defined( MODEL_GEOS )
       SpcInfo => State_Chm%SpcData(I)%Info
#else
       ! Get info about this species from the species database
       N = State_Chm%Map_Advect(I)
       ThisSpc => State_Chm%SpcData(N)%Info
#endif

       ! Pass tracer name
#if defined( MODEL_GEOS )
       Int2Spc(I)%Name = TRIM(SpcInfo%Name)
#else
       Int2Spc(I)%TrcName = TRIM(ThisSpc%Name)
#endif

       ! Get tracer ID
#if defined( MODEL_GEOS )
       ! Do not use Ind_() to get tracer ID. It may return the wrong tracer ID
       ! for species with the same hash (ckeller, 10/6/17)
       !!!Int2Spc(I)%ID = Ind_( TRIM(Int2Spc(I)%Name) )
       Int2Spc(I)%ID = SpcInfo%ModelID
#else
       Int2Spc(I)%TrcID = IND_( TRIM(Int2Spc(I)%TrcName) )
#endif

       ! If tracer ID is not valid, make sure all vars are at least defined.
#if defined( MODEL_GEOS )
       IF ( Int2Spc(I)%ID <= 0 ) THEN
          Int2Spc(I)%Internal => NULL()
          CYCLE
       ENDIF
#else
       IF ( Int2Spc(I)%TrcID <= 0 ) THEN
          Int2Spc(I)%Internal => NULL()
          CYCLE
       ENDIF
#endif

       ! Get internal state field
#if defined( MODEL_GEOS )
       fieldName = TRIM(SPFX)//TRIM(Int2Spc(I)%Name)
       CALL ESMF_StateGet( INTSTATE, TRIM(fieldName), GcFld, RC=STATUS )
       IF ( STATUS /= ESMF_SUCCESS ) THEN
          fieldName = TRIM(TPFX)//TRIM(Int2Spc(I)%Name)
          CALL ESMF_StateGet( INTSTATE, TRIM(fieldName), GcFld, RC=STATUS )
       ENDIF 
#else
       CALL ESMF_StateGet( INTSTATE, TRIM(SPFX) // TRIM(Int2Spc(I)%TrcName), &
                           GcFld, RC=STATUS )
#endif

       ! This is mostly for testing 
       IF ( STATUS /= ESMF_SUCCESS ) THEN
          IF( am_I_Root ) THEN
             WRITE(*,*) 'Cannot find in internal state: ', TRIM(SPFX) &
#if defined( MODEL_GEOS )
                        //TRIM(Int2Spc(I)%Name),I
          ENDIF
          Int2Spc(I)%Internal => NULL()
          CYCLE
#else
                        //TRIM(Int2Spc(I)%TrcName),I
          ENDIF
          _ASSERT(.FALSE.,'Error finding internal state variable')
#endif
       ENDIF

#if defined( MODEL_GEOS )
       ! Check friendliness of field: the field must not be friendly to 
       ! moist and/or turbulence if the corresponding GEOS-Chem switches 
       ! are turned on!
       ! Check for friendliness to convection: only if GEOS-Chem convection
       ! is enabled
       IF ( Input_Opt%LCONV ) THEN
          FRIENDLY=.FALSE.
          CALL ESMF_AttributeGet( GcFld, NAME="FriendlyToMOIST", &
                                  VALUE=FRIENDLY, RC=STATUS )
          IF ( STATUS==ESMF_SUCCESS .AND. FRIENDLY ) THEN
             IF ( am_I_Root ) THEN
                WRITE(*,*) ' '
                WRITE(*,*) 'GEOS-Chem convection is turned on, but ' &
                           // 'tracer is also'
                WRITE(*,*) 'friendly to MOIST. Cannot do both: ',    &
                           TRIM(Int2Spc(I)%Name)
                WRITE(*,*) ' '
             ENDIF
             _ASSERT(.FALSE.,'Error in Friendly settings')
          ENDIF
       ENDIF

       ! Check for friendliness to turbulence: only if GEOS-Chem turbulence
       ! is enabled
       IF ( Input_Opt%LTURB ) THEN
          FRIENDLY=.FALSE.
          CALL ESMF_AttributeGet( GcFld, NAME="FriendlyToTURBULENCE", &
                                  VALUE=FRIENDLY, RC=STATUS )
          IF ( STATUS==ESMF_SUCCESS .AND. FRIENDLY ) THEN
             IF ( am_I_Root ) THEN
                WRITE(*,*) ' '
                WRITE(*,*)    &
                   'GEOS-Chem turbulence is turned on, but tracer is also'
                WRITE(*,*) 'friendly to TURBULENCE. Cannot do both: ', &
                    TRIM(Int2Spc(I)%Name)
                WRITE(*,*) ' '
             ENDIF
             _ASSERT(.FALSE.,'Error in Friendly settings')
          ENDIF
       ENDIF

       ! Get pointer to field
       CALL ESMF_FieldGet( GcFld, 0, Ptr3D, __RC__ )
       Int2Spc(I)%Internal => Ptr3D

       ! Verbose
       !IF ( am_I_Root ) THEN
       !    WRITE(*,*) 'Connected ',TRIM(SpcInfo%Name), &
       !    ' <--> ',TRIM(fieldName),'; MW: ',SpcInfo%EmMW_g
       !ENDIF

       ! Free pointers
       Ptr3D   => NULL()
       SpcInfo => NULL()
#else
       ! Get pointer to field
       CALL ESMF_FieldGet( GcFld, 0, Ptr3D_R8, __RC__ )
       Int2Spc(I)%Internal => Ptr3D_R8

       ! Free pointers
       Ptr3D_R8 => NULL()
       ThisSpc  => NULL()
#endif

    ENDDO

#if defined( MODEL_GEOS )
    !=======================================================================
    ! Make sure that GEOS-Chem calculates chemistry tendencies of O3 and H2O
    ! if it is the analysis OX and RATS provider, respectively. 
    !=======================================================================
    IF ( DoAnox ) THEN
       GCID = ind_('O3')
       IF ( GCID > 0 ) THEN
          CALL Tend_CreateClass( am_I_Root, Input_Opt, State_Chm,             &
                                 'CHEM',    __RC__ )
          CALL Tend_Add        ( am_I_Root, Input_Opt, State_Chm, State_Grid, &
                                 'CHEM',    GCID,      __RC__ )
       ENDIF
    ENDIF

    IF ( DoRATS ) THEN
       GCID = ind_('H2O')
       IF ( GCID > 0 ) THEN
          CALL Tend_CreateClass( am_I_Root, Input_Opt, State_Chm,             &
                                 'CHEM',    __RC__ )
          CALL Tend_Add        ( am_I_Root, Input_Opt, State_Chm, State_Grid, &
                                 'CHEM',    GCID,      __RC__ )
       ENDIF
    ENDIF
#endif

    !=======================================================================
    ! Error trap: make sure that chemistry / emission time step are same and
    ! correspond to the chemistry step set in GEOSCHEMchem_GridComp.rc.
    !=======================================================================
    ChemTS = GET_TS_CHEM() 
    EmisTS = GET_TS_EMIS()
    IF ( ChemTS /= tsChem .OR. EmisTS /= tsChem ) THEN
       WRITE(*,*) 'GEOS-Chem chemistry and/or emission time step do not'
       WRITE(*,*) 'agree with time step set in GEOSCHEMchem_GridComp.rc'
       WRITE(*,*) 'GEOS-Chem chemistry time step                 : ', ChemTS
       WRITE(*,*) 'GEOS-Chem emission  time step                 : ', EmisTS
       WRITE(*,*) 'CHEMISTRY_TIMESTEP in GCHP.rc                 : ', tsChem
       _ASSERT(.FALSE.,'Error in timesteps')
    ENDIF

    ! Also check for convection and dynamics time step.
    ChemTS = GET_TS_CONV()
    EmisTS = GET_TS_DYN()
    IF ( ChemTS /= tsDyn .OR. EmisTS /= tsDyn ) THEN
       WRITE(*,*) 'GEOS-Chem transport and/or convection time step do not'
       WRITE(*,*) 'agree with time step set in GEOSCHEMchem_GridComp.rc'
       WRITE(*,*) 'GEOS-Chem convection time step                : ', ChemTS
       WRITE(*,*) 'GEOS-Chem dynamics   time step                : ', EmisTS
       WRITE(*,*) 'RUN_DT in CAP.rc                              : ', tsDyn
       _ASSERT(.FALSE.,'Error in timesteps')
    ENDIF

#if !defined( MODEL_GEOS )
    IF ( ArchivedConv .AND. am_I_Root ) THEN
       WRITE(*,*) ' '
       WRITE(*,*) ' --------------------------------------------------- '
       WRITE(*,*) ' GEOS-Chem will be using archived convection fields! '
       WRITE(*,*) ' --------------------------------------------------- '
       WRITE(*,*) ' '
    ENDIF
#endif

#if defined( MODEL_GEOS )
    !=======================================================================
    ! Read GEOSCHEMchem settings 
    !=======================================================================
    IF ( am_I_Root ) THEN
       WRITE(*,*) TRIM(Iam), ': options from GEOSCHEMchem_GridComp.rc:'
    ENDIF

    ! Apply correction term to large-scale precipitation that should be
    ! convective?
    CALL ESMF_ConfigGetAttribute( GeosCF, Input_Opt%WETD_CONV_SCAL,          &
                                  Label   = "Convective_precip_correction:", &
                                  Default = 1.0d0,                           &
                                  __RC__                                      )
    IF ( Input_Opt%WETD_CONV_SCAL < 0d0 ) Input_Opt%WETD_CONV_SCAL = 1.0d0
    Input_Opt%WETD_CONV_SCAL = MAX(MIN(1.0d0,Input_Opt%WETD_CONV_SCAL),0.0d0)
    IF ( am_I_Root ) THEN
       WRITE(*,*)   &
       '- Convective precip correction (0=no washout, 1=no correction): ', &
       Input_Opt%WETD_CONV_SCAL
    ENDIF

    ! Use GMI O3 P/L 
    CALL ESMF_ConfigGetAttribute( GeosCF, DoIt, Label = "Use_GMI_O3_PL:", &
                                  Default = 0, __RC__ ) 
    Input_Opt%LSYNOZ = .FALSE.
    IF ( DoIt == 1 ) THEN
       Input_Opt%LGMIOZ = .TRUE.
       Input_Opt%LLINOZ = .FALSE.
    ELSE
       Input_Opt%LGMIOZ = .FALSE.
       !!!Input_Opt%LSYNOZ = .NOT. Input_Opt%LLINOZ
    ENDIF
    IF ( am_I_Root ) THEN
       WRITE(*,*) '- Use GMIOZ: ', Input_Opt%LGMIOZ
       WRITE(*,*) '- Use LINOZ: ', Input_Opt%LLINOZ 
       WRITE(*,*) '- Use SYNOZ: ', Input_Opt%LSYNOZ
    ENDIF

    ! Overwrite strat. O3 with ANA_OZ 
    ! Default settings
    LANAO3  = .FALSE.
    ANAO3L1 = value_LLSTRAT 
    ANAO3L2 = value_LLSTRAT
    ANAO3L3 = LM
    ANAO3L4 = LM
    ANAO3FR = 1.0
    CALL ESMF_ConfigGetAttribute( GeosCF, DoIt, Label = "Use_ANA_O3:", &
                                  Default = 0, __RC__ ) 
    IF ( DoIt == 1 ) THEN
       LANAO3 = .TRUE.
       CALL ESMF_ConfigGetAttribute( GeosCF, ANAO3L1, Label = "ANAO3L1:", &
                                     Default = value_LLSTRAT, __RC__ ) 
       CALL ESMF_ConfigGetAttribute( GeosCF, ANAO3L2, Label = "ANAO3L2:", &
                                     Default = value_LLSTRAT, __RC__ ) 
       CALL ESMF_ConfigGetAttribute( GeosCF, ANAO3L3, Label = "ANAO3L3:", &
                                     Default = LM, __RC__ ) 
       CALL ESMF_ConfigGetAttribute( GeosCF, ANAO3L4, Label = "ANAO3L4:", &
                                     Default = LM, __RC__ ) 
       CALL ESMF_ConfigGetAttribute( GeosCF, ANAO3FR, Label = "ANAO3FR:", &
                                     Default = 1.0, __RC__ ) 
       CALL ESMF_ConfigGetAttribute( GeosCF, ANAO3FILE, Label = "ANAO3TMPL:", &
                                     Default = '/dev/null', __RC__ )
       CALL ESMF_ConfigGetAttribute( GeosCF, I, Label = "Use_PCHEM_O3:", &
                                     Default = 0, __RC__ )
       LPCHEMO3 = ( I == 1 )
       ASSERT_(ANAO3L1 >  0)
       ASSERT_(ANAO3L1 <= LM)
       ASSERT_(ANAO3L2 >= ANAO3L1)
       ASSERT_(ANAO3L2 <= LM)
       ASSERT_(ANAO3L4 >= ANAO3L3)
       ASSERT_(ANAO3L4 <= LM)
       ASSERT_(ANAO3L3 >  ANAO3L2)
       ANAO3FR = MAX(0.0,MIN(1.0,ANAO3FR))
    ENDIF
    IF ( am_I_Root ) THEN
       WRITE(*,*) 'Overwrite with analysis ozone?  ',LANAO3
       IF ( LANAO3 ) THEN 
          WRITE(*,*) '-> Use internal PCHEM field?     : ',LPCHEMO3
          WRITE(*,*) '-> Analysis ozone bottom level 1 : ',ANAO3L1
          WRITE(*,*) '-> Analysis ozone bottom level 2 : ',ANAO3L2
          WRITE(*,*) '-> Analysis ozone top level 3    : ',ANAO3L3
          WRITE(*,*) '-> Analysis ozone top level 4    : ',ANAO3L4
          WRITE(*,*) '-> Analysis ozone blend factor   : ',ANAO3FR
       ENDIF
    ENDIF

    ! Tropopause options 
    CALL ESMF_ConfigGetAttribute( GeosCF, DoIt,                    &
                                  Label = "Cap_polar_tropopause:", &
                                  Default = 1, __RC__ ) 
    Input_Opt%LCAPTROP = ( DoIt == 1 )
    IF ( am_I_Root ) THEN
       WRITE(*,*) '- Cap polar tropopause: ', Input_Opt%LCAPTROP
    ENDIF

    ! Check for pertubation of O3 field 
    CALL ESMF_ConfigGetAttribute( GeosCF, DoIt, &
            Label="PERTURB_O3:", Default=0, __RC__ )
    PerturbO3 = ( DoIt == 1 )
    CALL ESMF_ConfigGetAttribute( GeosCF, DoIt, &
            Label="PERTURB_CO:", Default=0, __RC__ )
    PerturbCO = ( DoIt == 1 )
    CALL ESMF_ConfigGetAttribute( GeosCF, Val, &
            Label="FIX_PERT:", Default=-999.0, __RC__ )
    FIXPERT = Val
    IF ( am_I_Root ) THEN
       WRITE(*,*) '- Perturb O3: ', PerturbO3
       WRITE(*,*) '- Perturb CO: ', PerturbCO
       WRITE(*,*) '- Fix pert  : ', FIXPERT
    ENDIF

    !CALL ESMF_ConfigGetAttribute( GeosCF, Input_Opt%NOx_sensitivity,   &
    !                              Label   = "NOx_sensitivity_factor:", &
    !                              Default = -999.0d0,                  &
    !                              __RC__                                )
    !IF ( am_I_Root ) THEN
    !   WRITE(*,*) '- Use NOx sensitivity factor: ', Input_Opt%NOx_sensitivity
    !ENDIF

    ! Get internal state from external data 
    CALL ESMF_ConfigGetAttribute( GeosCF, DoIt, Label = "INIT_SPC_FROM_FILE:", &
                                  Default = 0, __RC__ ) 
    InitFromFile = ( DoIt == 1 )

    ! Always set stratospheric H2O 
    CALL ESMF_ConfigGetAttribute( GeosCF, DoIt, & 
          Label="Prescribe_strat_H2O:", Default=0, __RC__ )
    Input_Opt%AlwaysSetH2O = ( DoIt == 1 )
    IF ( am_I_Root ) THEN
       WRITE(*,*) '- Prescribe H2O in stratosphere: ', Input_Opt%AlwaysSetH2O
    ENDIF

    ! Compute vertical updraft velocity from online values 
    CALL ESMF_ConfigGetAttribute( GeosCF, DoIt, & 
          Label="Calc_VUD_online:", Default=0, __RC__ )
    Input_Opt%UseOnlineVUD = ( DoIt == 1 )
    IF ( am_I_Root ) THEN
       WRITE(*,*) '- Compute VUD online: ', Input_Opt%UseOnlineVUD
    ENDIF

    ! Turn on Family Transport
    CALL Init_GCC_Chem_Groups()

    !=======================================================================
    ! CH4 error checks 
    !=======================================================================

    ! CH4 surface boundary conditions have to be turned off in GEOS-5
    IF ( Input_Opt%LCH4SBC ) THEN
       IF ( am_I_Root ) THEN
          WRITE(*,*) 'Please disable CH4 boundary conditions in input.geos.rc'
          WRITE(*,*) 'CH4 boundary conditions will automatically be applied'
          WRITE(*,*) 'if the CH4 emissions flag in input.geos.rc is turned off.'
          WRITE(*,*) '=> Turn on surface BCs :---'
          WRITE(*,*) '   => CH4?             : F'
       ENDIF
       ASSERT_(.FALSE.)
    ENDIF

    ! If CH4 emissions are turned on, boundary conditions will not be set
    ! and CH4 emissions are expected to be provided via HEMCO 
    IF ( am_I_Root ) THEN
       IF ( Input_Opt%LCH4EMIS ) THEN
          WRITE(*,*) 'CH4 emissions are turned on - no CH4 boundary conditions will be applied'
          WRITE(*,*) 'and CH4 emissions are taken from HEMCO'
       ELSE
          WRITE(*,*) 'CH4 emissions are turned off - CH4 boundary conditions will be applied'
       ENDIF
    ENDIF

    !=======================================================================
    ! All done
    !=======================================================================
#endif

    ! Stop timers
    ! -----------
    CALL MAPL_TimerOff( STATE, "INITIALIZE")

    CALL MAPL_TimerOff( STATE, "TOTAL")

    ! Successful return
    _RETURN(ESMF_SUCCESS)

    ! Formats
100 FORMAT( '### ',                                           / &
            '### ', a ,                                       / &
            '### ', a, '  |  Initialization on PET # ', i5.5, / &
            '### ' )
200 FORMAT( '### ',                                           / &
            '### ', a, '  |  Execution on PET # ',      i5.5, / &
            '###' )

  END SUBROUTINE Initialize_

!EOC
!------------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 910.1 and      !
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Run1
!
! !DESCRIPTION: Run1 is a wrapper method for the phase 1 run phase of the 
!  GEOSCHEMchem gridded component. It calls down to the Run method of the 
!  GEOS-Chem column chemistry code.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Run1 ( GC, Import, Export, Clock, RC )
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ESMF_GridComp), INTENT(INOUT) :: GC       ! Ref to this GridComp
    TYPE(ESMF_State),    INTENT(INOUT) :: Import   ! Import State
    TYPE(ESMF_State),    INTENT(INOUT) :: Export   ! Export State
    TYPE(ESMF_Clock),    INTENT(INOUT) :: Clock    ! ESMF Clock object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,             INTENT(OUT)   :: RC       ! Error return code
!
! !REMARKS:
!
! !REVISION HISTORY:
!  22 Sep 2014 - C. Keller   - Initial version.
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    CHARACTER(LEN=ESMF_MAXSTR)  :: compName    ! Name of gridded component
    CHARACTER(LEN=ESMF_MAXSTR)  :: Iam
    INTEGER                     :: STATUS
    INTEGER                     :: PHASE

    !=======================================================================
    ! Run1 starts here 
    !=======================================================================

    ! Set up traceback info
    CALL ESMF_GridCompGet( GC, name=compName, __RC__ )

    ! Identify this routine to MAPL
    Iam = TRIM(compName)//'::Run1'

    ! Call run routine stage 1 if more than one phase. If not 2 phases, 
    ! such as in GCHP, then we do all chemistry related processes from 
    ! Run2 instead.
    IF ( NPHASE == 2 ) THEN
       PHASE = 1
       CALL Run_ ( GC, IMPORT, EXPORT, CLOCK, PHASE, __RC__ )
    ENDIF

    ! Return w/ success
    _RETURN(ESMF_SUCCESS)

  END SUBROUTINE Run1 
!EOC
!------------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 910.1 and      !
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Run2
!
! !DESCRIPTION: Run2 is a wrapper method for the phase 2 run phase of the 
!  GEOSCHEMchem gridded component. It calls down to the Run method of the 
!  GEOS-Chem column chemistry code.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Run2 ( GC, Import, Export, Clock, RC )
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ESMF_GridComp), INTENT(INOUT) :: GC       ! Ref to this GridComp
    TYPE(ESMF_State),    INTENT(INOUT) :: Import   ! Import State
    TYPE(ESMF_State),    INTENT(INOUT) :: Export   ! Export State
    TYPE(ESMF_Clock),    INTENT(INOUT) :: Clock    ! ESMF Clock object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,             INTENT(OUT)   :: RC       ! Error return code
!
! !REMARKS:
!
! !REVISION HISTORY:
!  22 Sep 2014 - C. Keller   - Initial version.
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    CHARACTER(LEN=ESMF_MAXSTR)  :: compName    ! Name of gridded component
    CHARACTER(LEN=ESMF_MAXSTR)  :: Iam
    INTEGER                     :: PHASE
    INTEGER                     :: STATUS

    !=======================================================================
    ! Run2 starts here 
    !=======================================================================

    ! Set up traceback info
    CALL ESMF_GridCompGet( GC, name=compName, __RC__ )

    ! Identify this routine to MAPL
    Iam = TRIM(compName)//'::Run2'

    ! Set phase number: this is 2 for multi-phase runs (e.g. GEOS-5), and
    ! is -1 for single-phase runs (e.g. GCHP). If set to -1, all processes 
    ! are called (drydep, emissions, chemistry, etc.)
    IF ( NPHASE == 1 ) THEN
       PHASE = -1
    ELSE
       PHASE = 2
    ENDIF

    ! Call run routine stage 2
    CALL Run_ ( GC, IMPORT, EXPORT, CLOCK, PHASE, __RC__ )

    ! Return w/ success
    _RETURN(ESMF_SUCCESS)

  END SUBROUTINE Run2 
!EOC
!------------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 910.1 and      !
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Run_
!
! !DESCRIPTION: Run_ is the run method of the GEOSCHEMchem gridded component.  
!  GC is a simple ESMF/MAPL wrapper which calls down to the Run method of 
!  the GEOS-Chem column chemistry code.
!  Note: this routine currently skips the call down to GEOS-Chem on the very 
!  first time it is invoked. The reason is that a number of met-variables seem 
!  to be undefined still (e.g. BXHEIGHT, T, etc), yielding to seg-faults and/or 
!  crazy results. 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Run_( GC, Import, Export, Clock, Phase, RC )
!
! !USES:
!
    USE HCO_INTERFACE_MOD,       ONLY : HcoState
    USE MAPL_MemUtilsMod
    USE Olson_Landmap_Mod,       ONLY : Compute_Olson_Landmap
    USE Precision_Mod

#if defined( MODEL_GEOS )
    ! To store archived variables in internal state (ckeller, 9/16/15)
    USE CARBON_MOD,              ONLY : ORVC_SESQ

    ! To archive selected reaction rates
    USE GCKPP_Parameters
    USE GCKPP_Monitor
    USE CMN_FJX_MOD,             ONLY : JVN_

    ! Fast-JX diagnostics 
    USE FAST_JX_MOD,             ONLY : EXTRAL_NLEVS, EXTRAL_NITER
#endif

!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ESMF_GridComp), INTENT(INOUT), TARGET :: GC     ! Ref to this GridComp
    TYPE(ESMF_State),    INTENT(INOUT), TARGET :: Import ! Import State
    TYPE(ESMF_State),    INTENT(INOUT), TARGET :: Export ! Export State
    TYPE(ESMF_Clock),    INTENT(INOUT)         :: Clock  ! ESMF Clock object
    INTEGER,             INTENT(IN   )         :: Phase  ! Run phase (-1/1/2)
!
! !OUTPUT PARAMETERS:
!
    INTEGER,             INTENT(  OUT)         :: RC     ! Error return code
!
! !REMARKS:
!  We call routine Extract_ to return various values (i.e. grid parameters,
!  start & end dates, PET information, etc.) from the ESMF/MAPL environment.  
!  We then pass those to GEOS-Chem via routine GIGC_CHUNK_RUN, which is
!  located in GEOS-Chem module ./GEOS-Chem/ESMF/gigc_chunk_mod.F90.

! !REVISION HISTORY:
!  06 Dec 2009 - A. da Silva - Initial version
!  08 Apr 2010 - R. Yantosca - Now uses the updated Extract_ method
!  09 Apr 2010 - R. Yantosca - Initialize Timing, GeoLoc objects
!  16 Apr 2010 - R. Yantosca - Now move the array assignments before & after
!                              the call to GC_CHUNK_RUN into separate
!                              include files, for clarity
!  30 Apr 2010 - R. Yantosca - Now use 5 digits for PET
!  02 Jun 2010 - R. Yantosca - Now use IDENT%VERBOSE to trigger debug output
!  09 Oct 2012 - R. Yantosca - Now call MAPL_Am_I_Root to test for root PET
!  16 Oct 2012 - R. Yantosca - Now Include freeform files Includes_Before_Run.H
!                              and Includes_After_Run.H
!  13 Feb 2013 - R. Yantosca - Now call MAPL_Get_SunInsolation to return
!                              solar zenith angle and related properties
!  08 Mar 2013 - R. Yantosca - Now save the PET # (aka PET #) in Input_Opt 
!  15 Mar 2013 - R. Yantosca - Remove IDENT object, which was a holdover from
!                              the GEOS-Chem column code
!  22 Sep 2014 - C. Keller   - Added Phase argument
!  24 Oct 2014 - C. Keller   - Now derive all O3 export quantities from Tracers
!                              instead of Species (Species are zero in the 
!                              stratosphere!). Removed species from internal
!                              state as no physics was applied to them anyways.
#if !defined( MODEL_GEOS )
!  29 Nov 2016 - E. Lundgren - Initialize Olson fractional land type, MODIS 
!                              LAI, and MODIS CHLR from ExtData imports
!  01 Sep 2017 - E. Lundgren - Enable automation of GCHP diagnostics by 
!                              setting data pointers and copying GC state values
!                              using the new HistoryConfig object
!  26 Jul 2017 - S. Eastham  - Read LAI from a single variable in file
!  06 Nov 2017 - E. Lundgren - Abstract provider services to new module
!                              gigc_providerservices_mod.F90
#endif
!EOP
!------------------------------------------------------------------------------
!BOC

!
! LOCAL VARIABLES:
!  
    ! Objects
    TYPE(ESMF_Grid)              :: Grid          ! ESMF Grid object
    TYPE(ESMF_Config)            :: MaplCF        ! Config (MAPL.rc)
    TYPE(ESMF_Config)            :: GeosCF        ! Config (GEOSCHEM*.rc)
    TYPE(ESMF_Alarm)             :: ALARM
    TYPE(ESMF_VM)                :: VM            ! ESMF VM object
                                                  
    ! Scalars                                     
    LOGICAL                      :: am_I_Root     ! Are we on the root PET?
    LOGICAL                      :: IsChemTime    ! Chemistry alarm proxy
    LOGICAL                      :: IsRunTime     ! Time to call GEOS-Chem
    LOGICAL                      :: IsTendTime    ! Time to calculate tendencies
    INTEGER                      :: IND           ! Species or tracer index
    INTEGER                      :: error         ! G-C error return code
    INTEGER(ESMF_KIND_I8)        :: advCount      ! # of clock advances
    INTEGER                      :: nymd          ! YYYY/MM/DD date
    INTEGER                      :: nhms          ! hh:mm:ss time
    INTEGER                      :: myPet         ! PET # we are on now
    INTEGER                      :: nPets         ! Total # of PETs
    INTEGER                      :: I, J, L       ! Loop indices
    INTEGER                      :: IM,JM,LM      ! Grid dimensions
    INTEGER                      :: LR, N         ! Loop indices
    INTEGER                      :: N_TRC         ! Shadow var: # of tracers
    INTEGER                      :: year          ! Current year    
    INTEGER                      :: month         ! Current month
    INTEGER                      :: day           ! Current day
    INTEGER                      :: dayOfYr       ! Current day of year
    INTEGER                      :: hour          ! Current hour
    INTEGER                      :: minute        ! Current minute
    INTEGER                      :: second        ! Current second
    REAL                         :: UTC           ! Universal time
    REAL                         :: tsChem        ! Chem timestep [sec]
    REAL                         :: tsDyn         ! Dynamic timestep [sec]
    REAL                         :: hElapsed      ! Elapsed time [hours]
    REAL*8                       :: lonDeg        ! Longitude [degrees]
    REAL*8                       :: latDeg        ! Latitude [degrees]
    REAL*8                       :: P1, P2        ! Pressure variables
    CHARACTER(LEN=4)             :: petStr        ! String for PET #
    CHARACTER(LEN=ESMF_MAXSTR)   :: compName      ! Gridded Component name
    
    ! Allocatable local arrays
    REAL,  ALLOCATABLE, TARGET   :: zenith(:,:)   ! Solar zenith angle
    REAL,  ALLOCATABLE, TARGET   :: solar(:,:)    ! Solar insolation

    ! Pointer arrays needed to initialize from imports
    REAL, POINTER                :: Ptr2d   (:,:)   => NULL()
    REAL, POINTER                :: Ptr3d   (:,:,:) => NULL()
    REAL(ESMF_KIND_R8), POINTER  :: Ptr2d_R8(:,:)   => NULL()
    REAL(ESMF_KIND_R8), POINTER  :: Ptr3d_R8(:,:,:) => NULL()

    ! Other pointer arrays
    REAL(ESMF_KIND_R4),  POINTER :: lonCtr  (:,:) ! Lon centers, this PET [rad]
    REAL(ESMF_KIND_R4),  POINTER :: latCtr  (:,:) ! Lat centers, this PET [rad
    TYPE(MAPL_MetaComp), POINTER :: STATE

    ! For CTM Mode
    ! ckeller, 8/22/19: In GEOS, PLE and AIRDENS are from the IMPORT state
#if !defined( MODEL_GEOS )
    REAL(ESMF_KIND_R8),  POINTER :: PLE(:,:,:)     => NULL() ! INTERNAL: PEDGE
    REAL,                POINTER :: AIRDENS(:,:,:) => NULL() ! INTERNAL: PEDGE
#endif

    ! For aero bundle
    INTEGER                      :: nAero
    TYPE(ESMF_STATE)             :: INTSTATE
    TYPE(ESMF_FieldBundle)       :: AeroBdl 
    TYPE(ESMF_Field)             :: AeroFld
    REAL, POINTER                :: GcPtr3d  (:,:,:) => NULL()
    REAL, POINTER                :: AeroPtr3d(:,:,:) => NULL()
    CHARACTER(LEN=ESMF_MAXSTR)   :: GcName
    REAL                         :: GCMW, FRAC

    ! Initialize variables used for reading Olson and MODIS LAI imports
    INTEGER            :: TT, VV, landTypeInt
    CHARACTER(len=64)  :: landTypeStr, varName, importName

#if defined( MODEL_GEOS )
    ! GEOS-5 only local variables
    INTEGER                      :: LB            ! Loop indices
    REAL                         :: lp1, lp2      ! lightning potentials
    REAL                         :: DFPAR_MAX     ! Global max of DFPAR
    REAL(fp)                     :: TS_TEND       ! tendency time step

    ! Other pointer arrays needed by Include_Before_Run.H
    REAL, POINTER                :: PtrTmp(:,:,:) => NULL()
    REAL, POINTER                :: PtrEmis(:,:)  => NULL() 

    TYPE(ESMF_STATE)             :: Aero

    ! Initialize everything to zero (from registry file)?
    INTEGER                      :: InitZero 

    ! Initialize species to values given in globchem.dat (from rc file)?
    INTEGER                      :: InitSpecs

    ! Some checks for replay runs
    LOGICAL                      :: FIRSTREWIND
    LOGICAL                      :: AFTERREWIND
    LOGICAL                      :: IsFirst
    INTEGER, SAVE                :: pymd = 0         ! previous date
    INTEGER, SAVE                :: phms = 0         ! previous time
    INTEGER, SAVE                :: nnRewind = 0

    ! NO2 columns
    REAL, POINTER                :: PTR_NO2(:,:,:)
    REAL, POINTER                :: TNO2(:,:)
    REAL, POINTER                :: SNO2(:,:)
    REAL, POINTER                :: FNO2(:,:)

    ! Perturb O3
    REAL                         :: Rnd(2)
    REAL, POINTER                :: Ptr3dA(:,:,:) => NULL() 
    REAL, POINTER                :: Ptr3dB(:,:,:) => NULL()

    ! NO2 export
    REAL, POINTER                :: TROPN(:,:), STRATN(:,:)

    ! Diagnose reaction rates
    INTEGER                      :: NDIAG
    INTEGER                      :: RRids(NREACT)
    INTEGER                      :: JVids(JVN_)
    CHARACTER(LEN=3)             :: III
    CHARACTER(LEN=ESMF_MAXSTR)   :: RRName      ! Name of reaction name 

    ! OX
    REAL, POINTER                :: OX(:,:,:)
    REAL, POINTER                :: O3(:,:,:)
    REAL, POINTER                :: O3PPMV(:,:,:)
    REAL, POINTER                :: OX_TEND(:,:,:)
    REAL, POINTER                :: PTR_O3P(:,:,:)
    REAL, POINTER                :: PTR_O1D(:,:,:)
    REAL, PARAMETER              :: OMW = 16.0

    ! GEOS-5 only (also in gigc_providerservices but don't use yet):
    ! -RATS:
    REAL, POINTER     :: CH4     (:,:,:) => NULL()
    REAL, POINTER     :: N2O     (:,:,:) => NULL()
    REAL, POINTER     :: CFC11   (:,:,:) => NULL()
    REAL, POINTER     :: CFC12   (:,:,:) => NULL()
    REAL, POINTER     :: HCFC22  (:,:,:) => NULL()

    REAL, POINTER     :: PTR_O3      (:,:,:) => NULL()
    REAL, POINTER     :: PTR_CH4     (:,:,:) => NULL()
    REAL, POINTER     :: PTR_N2O     (:,:,:) => NULL()
    REAL, POINTER     :: PTR_CFC11   (:,:,:) => NULL()
    REAL, POINTER     :: PTR_CFC12   (:,:,:) => NULL()
    REAL, POINTER     :: PTR_HCFC22  (:,:,:) => NULL()
    REAL, POINTER     :: PTR_H2O     (:,:,:) => NULL()

    ! GCCTO3 and GCCTTO3 are the pointers to the corresponding export state fields
    REAL, POINTER     :: PTR_GCCTO3 (:,:) => NULL()
    REAL, POINTER     :: PTR_GCCTTO3(:,:) => NULL()

#else
    ! GCHP only local variables
    INTEGER                      :: trcID, RST
    REAL                         :: COEFF
    CHARACTER(LEN=ESMF_MAXSTR)   :: trcNAME,hcoNAME
    TYPE(ESMF_Field      )       :: trcFIELD
    TYPE(ESMF_FieldBundle)       :: trcBUNDLE
    REAL              , POINTER  :: fPtrArray(:,:,:)
    REAL(ESMF_KIND_R8), POINTER  :: fPtrVal, fPtr1D(:)
    INTEGER                      :: IMAXLOC(1)

#endif

    ! First call?
    LOGICAL, SAVE                :: FIRST = .TRUE.

    __Iam__('Run_')

    !=======================================================================
    ! Run starts here
    !=======================================================================

    ! Are we on the root PET?
    am_I_Root = MAPL_Am_I_Root()

    ! Set up traceback info
    CALL ESMF_GridCompGet( GC, name=compName, __RC__ )

    ! Identify this routine to MAPL
    Iam = TRIM(compName)//'::Run_'

#if !defined( MODEL_GEOS )
    ! Get the VM for optional memory prints (level >= 2)
    !-----------------------------------
    if ( MemDebugLevel > 0 ) THEN
       call ESMF_VmGetCurrent(VM, RC=STATUS)
       _VERIFY(STATUS)
    endif
#endif

    ! Get my MAPL_Generic state
    ! -------------------------
    CALL MAPL_GetObjectFromGC(GC, STATE, __RC__)

    ! Query the chemistry alarm.
    ! This checks if it's time to do chemistry, based on the time step
    ! set in AGCM.rc (GEOSCHEMCHEM_DT:). If the GEOS-Chem time step is not
    ! specified in AGCM.rc, the heartbeat will be taken (set in MAPL.rc).
    ! ----------------------------------------------------------------------
    CALL MAPL_Get(STATE, RUNALARM=ALARM, __RC__)
    IsChemTime = ESMF_AlarmIsRinging(ALARM, __RC__)

    ! Turn off alarm: only if it was on and this is phase 2 (don't turn off
    ! after phase 1 since this would prevent phase 2 from being executed).
    IF ( IsChemTime .AND. PHASE /= 1 ) THEN
       CALL ESMF_AlarmRingerOff(ALARM, __RC__ )
    ENDIF

    ! Get Internal state
    CALL MAPL_Get ( STATE, INTERNAL_ESMF_STATE=INTSTATE, __RC__ )

    ! ----------------------------------------------------------------------
    ! Check if we need to call the GEOS-Chem driver. The GEOS-Chem driver 
    ! contains routines for the following processes:
    !
    ! Phase 1:
    ! (1) Convection:     --> Dynamics time step  (optional)
    ! (2) Dry deposition  --> Chemistry time step
    ! (3) Emissions       --> Chemistry time step
    !
    ! Phase 2:
    ! (4) Turbulence      --> Dynamics time step  (optional)
    ! (5) Chemistry       --> Chemistry time step
    ! (6) Wet deposition  --> Dynamics time step
    ! 
    ! Phase -1:
    ! Includes all of the above
    !
    ! Convection and turbulence are only called if the corresponding 
    ! switches are turned on in the GEOS-Chem input file (input.geos).
    !
    ! To avoid unnecessary calls to the GEOS-Chem driver routine, we 
    ! check here if it's time for any of the processes listed above.
    ! The IsChemTime variable will be passed down to the GEOS-Chem driver
    ! to ensure that chemistry is only executed if it's time to do so.
    !
    ! The O3 and H2O tendencies will only be calculated when doing chemistry 
    ! (set to zero otherwise). All other export variables are updated every 
    ! time GEOS-Chem is called.
    ! ----------------------------------------------------------------------
    IsRunTime = IsChemTime
    IF ( Input_Opt%LCONV .AND. Phase /= 2 ) IsRunTime = .TRUE.
    IF ( Input_Opt%LTURB .AND. Phase /= 1 ) IsRunTime = .TRUE.
    IF ( Input_Opt%LWETD .AND. Phase /= 1 ) IsRunTime = .TRUE.

#if defined( MODEL_GEOS )
    !!! always run
    !IsRunTime = .TRUE.
#endif

    ! Is it time to update tendencies?
    ! Tendencies shall only be updated when chemistry is done, which is 
    ! Phase -1 or 2.
    IsTendTime = ( IsChemTime .AND. Phase /= 1 )

    ! Start timers
    ! ------------
    CALL MAPL_TimerOn(STATE, "TOTAL")

    ! Get pointers to fields in import, internal, and export states defined
    ! in the registry file. This has to be done on the first call only.
    IF ( FIRST ) THEN
#if defined( MODEL_GEOS )
#      include "GEOSCHEMCHEM_GetPointer___.h"
#else
#      include "GIGCchem_GetPointer___.h"

       !IF ( IsCTM ) THEN
       call MAPL_GetPointer ( IMPORT, PLE,      'PLE',     __RC__ )
       call MAPL_GetPointer ( IMPORT, AIRDENS,  'AIRDENS', __RC__ )
       !ENDIF

       ! Set up pointers if GEOS-Chem is a provider
       !IF ( isProvider ) THEN
       CALL Provider_SetPointers( am_I_Root, EXPORT, calcOzone, __RC__ )
       !ENDIF

       ! Pass IMPORT/EXPORT object to HEMCO state object
       !CALL GetHcoState( HcoState )
       _ASSERT(ASSOCIATED(HcoState),'HcoState is not associated')
       HcoState%GRIDCOMP => GC
       HcoState%IMPORT   => IMPORT
       HcoState%EXPORT   => EXPORT
       !HcoState => NULL()

       ! To use archived convection fields
       IF ( ArchivedConv ) THEN
          CALL MAPL_GetPointer ( IMPORT, PTR_ARCHIVED_PFI_CN ,           &
                                 'ARCHIVED_PFI_CN'  , notFoundOK=.TRUE., &
                                 __RC__ )
          CALL MAPL_GetPointer ( IMPORT, PTR_ARCHIVED_PFL_CN ,           &
                                 'ARCHIVED_PFL_CN'  , notFoundOK=.TRUE., &
                                  __RC__ )
          CALL MAPL_GetPointer ( IMPORT, PTR_ARCHIVED_CNV_MFC,           &
                                 'ARCHIVED_CNV_MFC' , notFoundOK=.TRUE., &
                                 __RC__ )
          CALL MAPL_GetPointer ( IMPORT, PTR_ARCHIVED_CNV_MFD,           &
                                 'ARCHIVED_CNV_MFD' , notFoundOK=.TRUE., &
                                 __RC__ )
          CALL MAPL_GetPointer ( IMPORT, PTR_ARCHIVED_CNV_CVW,           &
                                 'ARCHIVED_CNV_CVW' , notFoundOK=.TRUE., &
                                 __RC__ )
          CALL MAPL_GetPointer ( IMPORT, PTR_ARCHIVED_DQRC   ,           &
                                 'ARCHIVED_DQRC'    , notFoundOK=.TRUE., &
                                 __RC__ )
          CALL MAPL_GetPointer ( IMPORT, PTR_ARCHIVED_REV_CN ,           &
                                 'ARCHIVED_PFI_CN'  , notFoundOK=.TRUE., &
                                 __RC__ )
          CALL MAPL_GetPointer ( IMPORT, PTR_ARCHIVED_T      ,           &  
                                 'ARCHIVED_T'       , notFoundOK=.TRUE., &
                                  __RC__ )
       ENDIF
#endif
    ENDIF

#if defined( MODEL_GEOS )
    ! Eventually get pointers to GCCTO3 and GCCTTO3. Those fields are optional
    ! and are only filled if defined and required.
    CALL MAPL_GetPointer ( EXPORT, PTR_GCCTO3,   'GCCTO3', notFoundOK=.TRUE., __RC__ )
    CALL MAPL_GetPointer ( EXPORT, PTR_GCCTTO3, 'GCCTTO3', notFoundOK=.TRUE., __RC__ )
!---

! GCHP ends the (if FIRST) block and then links HEMCO state to GC objects:
!    ENDIF

    CALL MAPL_GetPointer( INTSTATE, PTR_O3, 'TRC_O3', NotFoundOk=.TRUE., __RC__ )

    ! Get pointers to analysis OX exports
    IF ( DoANOX ) THEN
       CALL MAPL_GetPointer ( EXPORT, OX_TEND, 'OX_TEND' , __RC__ )
       CALL MAPL_GetPointer ( EXPORT,      OX, 'OX'      , __RC__ )
       CALL MAPL_GetPointer ( EXPORT,      O3, 'O3'      , __RC__ )
       CALL MAPL_GetPointer ( EXPORT,  O3PPMV, 'O3PPMV'  , __RC__ )
    ENDIF

    ! Get pointers to RATS exports
    IF ( DoRATS) THEN
       CALL MAPL_GetPointer ( EXPORT, H2O_TEND, 'H2O_TEND' , __RC__ )
       CALL MAPL_GetPointer ( EXPORT,      CH4, 'CH4'      , __RC__ )
       CALL MAPL_GetPointer ( EXPORT,      N2O, 'N2O'      , __RC__ )
       CALL MAPL_GetPointer ( EXPORT,    CFC11, 'CFC11'    , __RC__ )
       CALL MAPL_GetPointer ( EXPORT,    CFC12, 'CFC12'    , __RC__ )
       CALL MAPL_GetPointer ( EXPORT,   HCFC22, 'HCFC22'   , __RC__ )
       CALL MAPL_GetPointer( INTSTATE,    PTR_CH4, 'TRC_CH4',    __RC__ )
       CALL MAPL_GetPointer( INTSTATE,    PTR_N2O, 'TRC_N2O',    __RC__ )
       CALL MAPL_GetPointer( INTSTATE,  PTR_CFC11, 'TRC_CFC11',  __RC__ )
       CALL MAPL_GetPointer( INTSTATE,  PTR_CFC12, 'TRC_CFC12',  __RC__ )
       CALL MAPL_GetPointer( INTSTATE, PTR_HCFC22, 'TRC_HCFC22', __RC__ )
       CALL MAPL_GetPointer( INTSTATE,    PTR_H2O, 'TRC_H2O',    __RC__ )
   ENDIF

   ! Link HEMCO state to gridcomp objects
   ASSERT_(ASSOCIATED(HcoState))
   HcoState%GRIDCOMP => GC
   HcoState%IMPORT   => IMPORT
   HcoState%EXPORT   => EXPORT
#endif

    ! Run when it's time to do so
    ! Always run on first call to make sure that all variables become
    ! properly specified and initialized.
    ! ------------------------------------------------------------------
    RunningGEOSChem: IF(IsRunTime .OR. FIRST) THEN

       CALL MAPL_TimerOn(STATE, "RUN"  )

       ! Get various parameters from the ESMF/MAPL framework
       CALL Extract_( GC,                   &  ! Ref to this Gridded Component
                      Clock,                &  ! ESMF Clock object
                      Grid      = Grid,     &  ! ESMF Grid object
                      MaplCf    = MaplCF,   &  ! ESMF Config obj (MAPL*.rc) 
                      GeosCf    = GeosCF,   &  ! ESMF Config obj (GEOSCHEM*.rc)
                      tsChem    = tsChem,   &  ! Chemistry timestep [min]
                      tsDyn     = tsDyn,    &  ! Dynamic timestep [min]
                      nymd      = nymd,     &  ! Current YYYY/MM/DD date
                      nhms      = nhms,     &  ! Current hh:mm:ss time
                      year      = year,     &  ! Current year
                      month     = month,    &  ! Current month
                      day       = day,      &  ! Current day
                      dayOfYr   = dayOfYr,  &  ! Current day of year
                      hour      = hour,     &  ! Current hour
                      minute    = minute,   &  ! Current minute
                      helapsed  = hElapsed, &  ! Elapsed hours
                      advCount  = advCount, &  ! # of times clock has advanced
                      utc       = utc,      &  ! Universal time [hours]
                      localpet  = myPet,    &  ! # of the PET we are on now
                      petCount  = nPets,    &  ! Total # of PETs
                      __RC__ )

       ! For convenience, set grid dimension variables. These are being
       ! used in Includes_Before_Run.H (ckeller, 8/22/19)
       IM = State_Grid%NX
       JM = State_Grid%NY
       LM = State_Grid%NZ

       ! Allocate GMAO_ZTH (declared at top of module)
       IF ( .not. ALLOCATED( zenith ) ) THEN
          ALLOCATE( zenith(State_Grid%NX,State_Grid%NY), STAT=STATUS)
          _VERIFY(STATUS)
       ENDIF
       
       ! Allocate GMAO_SLR (declared @ top of module)
       IF ( .not. ALLOCATED( solar ) ) THEN
          ALLOCATE( solar(State_Grid%NX,State_Grid%NY), STAT=STATUS)
          _VERIFY(STATUS)
       ENDIF
       
       ! Call EXTRACT a second time to get the solar zenith
       ! angle and solar insolation fields
       CALL Extract_( GC,                   &  ! Ref to this Gridded Component
                      Clock,                &  ! ESMF Clock object
                      Grid      = Grid,     &  ! ESMF Grid object
                      MaplCf    = MaplCF,   &  ! ESMF Config obj (MAPL*.rc) 
                      GeosCf    = GeosCF,   &  ! ESMF Config obj (GEOSCHEM*.rc)
                      lonCtr    = lonCtr,   &  ! Lon centers on this PET [rad]
                      latCtr    = latCtr,   &  ! Lat centers on this PET [rad]
                      ZTH       = zenith,   &  ! Solar zenith angle
                      SLR       = solar,    &  ! Solar insolation
                      __RC__ )

#if !defined( MODEL_GEOS )
       ! MSL - shift from 0 - 360 to -180 - 180 degree grid
       where (lonCtr .gt. MAPL_PI ) lonCtr = lonCtr - 2*MAPL_PI
#endif

       ! Pass grid area [m2] obtained from dynamics component to State_Grid
       State_Grid%Area_M2 = AREA
       !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       ! KLUDGE (mps, 5/23/19):
       ! Copy to State_Met%AREA_M2 to avoid breaking GCHP benchmarks, which
       ! require the AREA_M2 field saved out to the StateMet diagnostic
       ! collection for things like computing emission totals.
       !
       State_Met%Area_M2 = State_Grid%Area_M2
       !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if defined( MODEL_GEOS )
       ! Check if this time is before the datetime of the prev timestep, e.g.
       ! if this is after a clock rewind
       AFTERREWIND = .FALSE.
       FIRSTREWIND = .FALSE.
       IF ( nymd < pymd ) THEN
          AFTERREWIND = .TRUE.
       ELSEIF ( (nymd == pymd) .AND. (nhms < phms) ) THEN
          AFTERREWIND = .TRUE.
       ENDIF
       
       ! If this is after a rewind, check if it's the first rewind. In this
       ! case, we need to re-do some first-time assignments to make sure that
       ! we reset all variables to the initial state!
       IF ( AFTERREWIND ) THEN
          nnRewind = nnRewind + 1
          IF ( nnRewind == 1 ) FIRSTREWIND = .TRUE.
       ENDIF
       !=======================================================================
       ! # of reaction rates to be diagnosed
       !=======================================================================
       IF ( Input_Opt%NN_RxnRates < 0 ) THEN
          NDIAG    = 0
          RRids(:) = -999
          DO I = 1, NREACT
             WRITE(III,'(i3.3)') I
             RRName = 'GCC_RR_'//TRIM(III) 
             CALL MAPL_GetPointer ( EXPORT, Ptr3D, TRIM(RRName), &
                                    NotFoundOk=.TRUE., __RC__ )
             IF ( ASSOCIATED(Ptr3D) ) THEN
                NDIAG        = NDIAG + 1
                RRids(NDIAG) = I
                IF ( am_I_Root ) THEN
                   WRITE(*,*) 'Will diagnose reaction rate ',TRIM(III),': ', &
                              TRIM(ADJUSTL(EQN_NAMES(I)))
                ENDIF
             ENDIF
          ENDDO
          Input_Opt%NN_RxnRates = NDIAG
          IF ( NDIAG > 0 ) THEN
             ALLOCATE(Input_Opt%RxnRates_IDs(NDIAG))
             Input_Opt%RxnRates_IDs(1:NDIAG) = RRids(1:NDIAG)
             ALLOCATE(State_Diag%RxnRate(State_Grid%NX,State_Grid%NY,&
                                          State_Grid%NZ,NDIAG))
             State_Diag%RxnRate = 0.0
          ENDIF
       ENDIF
       
       IF ( Input_Opt%NN_RxnRconst < 0 ) THEN
          NDIAG    = 0
          RRids(:) = -999
          DO I = 1, NREACT
             ! Always archive reaction coeffs needed for diagnostics
             IF ( I == id_rc_no .OR. I == id_rc_no2 ) THEN
                NDIAG        = NDIAG + 1
                RRids(NDIAG) = I
             ELSE
                WRITE(III,'(i3.3)') I
                RRName = 'GCC_RC_'//TRIM(III) 
                CALL MAPL_GetPointer ( EXPORT, Ptr3D, TRIM(RRName),  &
                                       NotFoundOk=.TRUE., __RC__ )
                IF ( ASSOCIATED(Ptr3D) ) THEN
                   NDIAG        = NDIAG + 1
                   RRids(NDIAG) = I
                   !IF ( am_I_Root ) THEN
                   !   WRITE(*,*) 'Will diagnose reaction rate constant ', &
                   !               TRIM(III),': ',TRIM(ADJUSTL(EQN_NAMES(I)))
                   !ENDIF
                ENDIF
             ENDIF
          ENDDO
          Input_Opt%NN_RxnRconst = NDIAG
          IF ( NDIAG > 0 ) THEN
             ALLOCATE(Input_Opt%RxnRconst_IDs(NDIAG))
             Input_Opt%RxnRconst_IDs(1:NDIAG) = RRids(1:NDIAG)
             ALLOCATE(State_Diag%RxnRconst(State_Grid%NX,State_Grid%NY, &
                                           State_Grid%NZ,NDIAG))
             State_Diag%RxnRconst = 0.0
          ENDIF
       ENDIF
       
       IF ( Input_Opt%NN_Jvals < 0 ) THEN
          NDIAG    = 0
          JVids(:) = -999
          DO I = 1, JVN_
             ! Always archive Jvalues needed for diagnostics
             IF ( I == id_jno2 ) THEN
                NDIAG        = NDIAG + 1
                JVids(NDIAG) = I
                IF ( am_I_Root ) THEN
                   WRITE(*,*) 'Will diagnose J-value',TRIM(III),NDIAG
                ENDIF
             ELSE
                WRITE(III,'(i3.3)') I
                RRName = 'GCC_JVAL_'//TRIM(III) 
                CALL MAPL_GetPointer ( EXPORT, Ptr3D, TRIM(RRName), &
                                       NotFoundOk=.TRUE., __RC__ )
                IF ( ASSOCIATED(Ptr3D) ) THEN
                   NDIAG        = NDIAG + 1
                   JVids(NDIAG) = I
                   IF ( am_I_Root ) THEN
                      WRITE(*,*) 'Will diagnose J-value',TRIM(III),NDIAG
                   ENDIF
                ENDIF
             ENDIF
          ENDDO
          Input_Opt%NN_JVals = NDIAG
          IF ( NDIAG > 0 ) THEN
             ALLOCATE(Input_Opt%JVal_IDs(NDIAG))
             Input_Opt%JVal_IDs(1:NDIAG) = JVids(1:NDIAG)
             ALLOCATE(State_Diag%JValIndiv(State_Grid%NX,State_Grid%NY, &
                                           State_Grid%NZ,NDIAG))
             State_Diag%JValIndiv = 0.0
          ENDIF
       ENDIF
#endif
       
       !=======================================================================
       ! Prevent use of (occasional) MAPL_UNDEF tropopause pressures
       !=======================================================================
       
       ! GCCTROPP contains the last valid tropopause pressure
       WHERE ( TROPP /= MAPL_UNDEF ) GCCTROPP = TROPP
       
       ! If any values in GCCTROPP are undefined, stop the run
       IF ( ANY( GCCTROPP == MAPL_UNDEF ) ) THEN
          PRINT *,TRIM(Iam)//": At least one invalid tropopause pressure."
          STATUS = GC_FAILURE
          _VERIFY(STATUS)
       ENDIF

!       !=======================================================================
!       ! pre-Run method array assignments. This passes the tracer arrays from
!       ! the internal state to State_Chm. On the first call, it also fills the
!       ! internal species arrays in State_Chm with the values read from the
!       ! restart file (and stored in the internal state).
!       !=======================================================================
       CALL MAPL_TimerOn(STATE, "CP_BFRE")
#include "Includes_Before_Run.H"
       CALL MAPL_TimerOff(STATE, "CP_BFRE")

#if !defined( MODEL_GEOS )
      !=======================================================================
      ! Pass advected tracers from internal state to GEOS-Chem tracers array
      !=======================================================================
      DO I = 1, SIZE(Int2Spc,1)
         IF ( Int2Spc(I)%TrcID <= 0 ) CYCLE
         State_Chm%Species(:,:,:,Int2Spc(I)%TrcID) = Int2Spc(I)%Internal
      ENDDO
      
      ! Flip in the vertical
      State_Chm%Species   = State_Chm%Species( :, :, State_Grid%NZ:1:-1, : )

       !=======================================================================
       ! On first call, also need to initialize the species from restart file.
       ! Only need to do this for species that are not advected, i.e. species
       ! that are not tracers (all other species arrays will be filled with
       ! tracer values anyways!).
       ! We only need to do this on the first call because afterwards, species
       ! array already contains values from previous chemistry time step
       ! (advected species will be updated with tracers)
       ! ckeller, 10/27/2014
       !=======================================================================
       IF ( FIRST ) THEN
       
          ! Get Generic State
          call MAPL_GetObjectFromGC ( GC, STATE, RC=STATUS)
          _VERIFY(STATUS)
          ! Get Internal state
          CALL MAPL_Get ( STATE, INTERNAL_ESMF_STATE=INTERNAL, __RC__ ) 
       
          ! Loop over all species and get info from spc db
          DO N = 1, State_Chm%nSpecies
             ThisSpc => State_Chm%SpcData(N)%Info
             IF ( TRIM(ThisSpc%Name) == '' ) CYCLE
             IND = IND_( TRIM(ThisSpc%Name ) )
             IF ( IND < 0 ) CYCLE
       
             ! Get data from internal state and copy to species array
             CALL MAPL_GetPointer( INTERNAL, Ptr3D_R8, TRIM(SPFX) //          &
                                   TRIM(ThisSpc%Name), notFoundOK=.TRUE.,     &
                                   __RC__ )
             IF ( .NOT. ASSOCIATED(Ptr3D_R8) ) THEN
                IF ( MAPL_am_I_Root()) WRITE(*,*)                             &
                   'Could not find species in INTERNAL state - will be ' //   &
                   'initialized to zero: ', TRIM(SPFX), TRIM(ThisSpc%Name)
                State_Chm%Species(:,:,:,IND) = 1d-26
                CYCLE
             ENDIF
             State_Chm%Species(:,:,:,IND) = Ptr3D_R8(:,:,State_Grid%NZ:1:-1)
             if ( MAPL_am_I_Root()) WRITE(*,*)                                &
             'Initialized species from INTERNAL state: ', TRIM(ThisSpc%Name)

             ! Determine if species in restart file
             CALL ESMF_StateGet( INTERNAL, TRIM(SPFX) // TRIM(ThisSpc%Name),  &
                  trcFIELD, RC=RC )
             CALL ESMF_AttributeGet( trcFIELD, NAME="RESTART",                &
                  VALUE=RST, RC=STATUS )
       
             ! Set spc conc to background value if rst skipped or var not there
             IF ( RC /= ESMF_SUCCESS .OR. RST == MAPL_RestartBootstrap .OR.   &
                      RST == MAPL_RestartSkipInitial ) THEN
                DO L = 1, State_Grid%NZ
                DO J = 1, State_Grid%NY
                DO I = 1, State_Grid%NX
                   IF ( L > State_Grid%MaxChemLev .AND. &
                            ( .NOT. ThisSpc%Is_Advected ) ) THEN
                      ! For non-advected spc at L > MaxChemLev, use small number
                      State_Chm%Species(I,J,L,IND) = 1.0E-30_FP           
                   ELSE
                      ! For all other cases, use the background value in spc db
                      State_Chm%Species(I,J,L,IND) = ThisSpc%BackgroundVV 
                   ENDIF
                ENDDO
                ENDDO
                ENDDO
                Ptr3D_R8(:,:,:) = State_Chm%Species(:,:,State_Grid%NZ:1:-1,IND)
                IF ( MAPL_am_I_Root()) THEN
                   WRITE(*,*)  &
                   '   WARNING: using background values from species database'
                ENDIF
             ENDIF
             ThisSpc => NULL()
          ENDDO
       ENDIF

       !=======================================================================
       ! On first call, initialize certain State_Chm and State_Met arrays from
       ! imports if they are found (ewl, 12/13/18)
       !=======================================================================
       IF ( FIRST ) THEN
          CALL MAPL_GetPointer( INTSTATE, Ptr3d_R8, 'H2O2AfterChem',  &
                                notFoundOK=.TRUE., __RC__ )
          IF ( ASSOCIATED(Ptr3d_R8) .AND. &
               ASSOCIATED(State_Chm%H2O2AfterChem) ) THEN
             State_Chm%H2O2AfterChem = Ptr3d_R8(:,:,State_Grid%NZ:1:-1)
          ENDIF
          Ptr3d_R8 => NULL()
          
          CALL MAPL_GetPointer( INTSTATE, Ptr3d_R8, 'SO2AfterChem',   &
                                notFoundOK=.TRUE., __RC__ )
          IF ( ASSOCIATED(Ptr3d_R8) .AND. &
               ASSOCIATED(State_Chm%SO2AfterChem) ) THEN
             State_Chm%SO2AfterChem = Ptr3d_R8(:,:,State_Grid%NZ:1:-1)
          ENDIF
          Ptr3d_R8 => NULL()
          
          CALL MAPL_GetPointer( INTSTATE, Ptr2d_R8, 'DryDepNitrogen', &
                                notFoundOK=.TRUE., __RC__ )
          IF ( ASSOCIATED(Ptr2d_R8) .AND. &
               ASSOCIATED(State_Chm%DryDepNitrogen) ) THEN
             State_Chm%DryDepNitrogen = Ptr2d_R8
          ENDIF
          Ptr2d_R8 => NULL()
          
          CALL MAPL_GetPointer( INTSTATE, Ptr2d_R8, 'WetDepNitrogen', &
                                notFoundOK=.TRUE., __RC__ )
          IF ( ASSOCIATED(Ptr2d_R8) .AND. &
               ASSOCIATED(State_Chm%WetDepNitrogen) ) THEN
             State_Chm%WetDepNitrogen = Ptr2d_R8
          ENDIF
          Ptr2d_R8 => NULL()
          
          CALL MAPL_GetPointer( INTSTATE, Ptr3d_R8, 'KPPHvalue' ,     &
                                notFoundOK=.TRUE., __RC__ )
          IF ( ASSOCIATED(Ptr3d_R8) .AND. &
               ASSOCIATED(State_Chm%KPPHvalue) ) THEN
             State_Chm%KPPHvalue(:,:,1:State_Grid%MaxChemLev) =       &
            Ptr3d_R8(:,:,State_Grid%NZ:State_Grid%NZ-State_Grid%MaxChemLev+1:-1)
          ENDIF
          Ptr3d_R8 => NULL()

          CALL MAPL_GetPointer( INTSTATE, Ptr3d_R8, 'DELP_DRY' ,     &
                                notFoundOK=.TRUE., __RC__ )
          IF ( ASSOCIATED(Ptr3d_R8) .AND. &
               ASSOCIATED(State_Met%DELP_DRY) ) THEN
             State_Met%DELP_DRY(:,:,1:State_Grid%NZ) =       &
                                  Ptr3d_R8(:,:,State_Grid%NZ:1:-1)
          ENDIF
          Ptr3d_R8 => NULL()
       ENDIF
#endif

#if defined( MODEL_GEOS )
       ! Eventually initialize species concentrations from external field. 
       IF ( InitFromFile ) THEN 
          IsFirst = ( FIRST .OR. FIRSTREWIND )
          CALL InitFromFile_( GC, Import, INTSTATE, Export, Clock, &
                              Input_Opt, State_Met, State_Chm, Q,  &
                              PLE, GCCTROPP, IsFirst, __RC__ ) 
       ENDIF

       !=======================================================================
       ! Error trap: make sure that PBL height is defined.
       ! Some fields needed by GEOS-Chem are only filled after the first 
       ! GEOS-Chem run. We need to avoid that GEOS-Chem is called in these 
       ! cases, since those fields are still zero and would cause seg-faults.
       ! The PBL height is a good proxy variable, and if those values are ok
       ! all others seem to be fine as well. 
       ! We do this error check on every time step (not only on the first
       ! one) to also catch the case where the time is reset to the initial
       ! conditions (replay mode).
       ! (ckeller, 4/24/2015). 
       !=======================================================================
       IF ( ANY(State_Met%PBLH <= 0.0_fp) ) THEN
          Input_Opt%haveImpRst = .FALSE. 

          ! Warning message
          IF ( am_I_Root ) THEN
             write(*,*) ' '
             write(*,*)      &
                'At least one PBLH value in GEOS-Chem is zero - skip time step'
             write(*,*) ' '
          ENDIF
       ENDIF

#if defined( MODEL_GEOS )
       !=======================================================================
       ! Also make sure that radiation fields are available. State_Met%OPTD is
       ! a good proxy since it is composed of three imports from RADIATION.
       ! (ckeller, 11/25/2015)
       ! GlobalSum does not seem to work properly anymore, skip this error
       ! trap (ckeller, 8/27/2019).
       !=======================================================================
       !DFPAR_MAX = GlobalSum( GC, DataPtr2D=DFPAR, maximum=.TRUE., __RC__ )
       !IF ( DFPAR_MAX == 0.0 ) THEN
       !   Input_Opt%haveImpRst = .FALSE. 
       !
       !   ! Warning message
       !   IF ( am_I_Root ) THEN
       !      write(*,*) ' '
       !      write(*,*)    &
       !            'All GEOS-Chem radiation imports are zero - skip time step'
       !      write(*,*) ' '
       !   ENDIF
       !ENDIF
#endif

       !=======================================================================
       ! Handling of species/tracer initialization. Default practice is to take
       ! whatever values are in the restarts. However, it is possible to
       ! initialize everything to zero and/or to set species' concentration to
       ! the values set in globchem.dat.rc. These options can be set in the
       ! GEOSCHEMchem GridComp registry. (ckeller, 2/4/16)
       !=======================================================================
       IF ( FIRST .OR. FIRSTREWIND ) THEN
          ! Check if zero initialization option is selected. If so, make sure 
          ! all concentrations are initialized to zero! 
          CALL ESMF_ConfigGetAttribute( GeosCF, InitZero, Default=0, &
                                        Label = "INIT_ZERO:", __RC__ ) 
          IF ( InitZero == 1 ) THEN
             State_Chm%Species = 0.0d0
             IF ( am_I_Root ) THEN
                write(*,*) ' '
                write(*,*) ' '
                write(*,*)     &
                 '### ALL GEOS-CHEM CONCENTRATIONS INITIALIZED TO ZERO !!! ###'
                write(*,*) ' '
                write(*,*) ' '
             ENDIF
          ENDIF
       
          ! Check if species shall be initialized to values set in globchem.dat 
          CALL ESMF_ConfigGetAttribute( GeosCF, InitSpecs, Default=0, &
                                        Label = "INIT_SPECS:", __RC__ )
          IF ( InitSpecs == 1 ) Input_Opt%LINITSPEC = .TRUE.
       ENDIF
#endif

       !=======================================================================
       ! Set Olson land map types from import of Olson file. 
       !=======================================================================
       ! We are currently using land type fractions derived from the 2001
       ! Olson land map instead of GEOS5 vegetation type fractions. Fractions
       ! are calculated by ExtData using conservate fractional regridding of
       ! the native 0.25x0.25 resolution file. (ewl, 11/29/16)
       !
       ! Previous:
       ! Set land types in State_Met from GEOS5 vegetation type fractions or
       ! OLSON land type fractions. For now, the land types are treated as 
       ! static and obtained from offline fields. The routine below thus needs
       ! to be called only once.
       ! Once the GEOS-5 land types are dynamic, we should import those from
       ! the surface component (field ITY, or better: vegetation type fractions
       ! per grid box).                                   (ckeller, 01/06/2015)
       !
       !=======================================================================
#if defined( MODEL_GEOS )
       IF ( FIRST .OR. FIRSTREWIND ) THEN
#else
       IF ( FIRST ) THEN
#endif

          ! Set Olson fractional land type from import (ewl)
          If (am_I_Root) Write(*,'(a)') 'Initializing land type ' // &
                           'fractions from Olson imports'
          Ptr2d => NULL()
          DO TT = 1, NSURFTYPE
          
             ! Create two-char string for land type
             landTypeInt = TT-1
             WRITE ( landTypeStr, '(I2.2)' ) landTypeInt
             importName = 'OLSON' // TRIM(landTypeStr)
          
             ! Get pointer and populate State_Met variable
             CALL MAPL_GetPointer ( IMPORT, Ptr2D, TRIM(importName),  &
                                    notFoundOK=.TRUE., __RC__ )
             If ( Associated(Ptr2D) ) Then
                State_Met%LandTypeFrac(:,:,TT) = Ptr2D(:,:)
             ELSE
                WRITE(6,*) TRIM(importName) // ' pointer is not associated'
             ENDIF
             Ptr2D => NULL()
          ENDDO

          ! Compute State_Met variables IREG, ILAND, IUSE, and FRCLND
          CALL Compute_Olson_Landmap( Input_Opt, State_Grid, State_Met, RC )
          _ASSERT(RC==GC_SUCCESS,'Error calling Compute_Olson_Landmap')
       ENDIF

#if defined( MODEL_GEOS )
       !=======================================================================
       ! Set ozone to values from PCHEM if this option is selected 
       !=======================================================================
       IF ( PHASE /= 2 .AND. LANAO3 ) THEN
          CALL SetAnaO3_( GC, Import, INTSTATE, Export, Clock, &
                          Input_Opt,  State_Met, State_Chm, Q, __RC__ ) 
       ENDIF
#endif

       !=======================================================================
       ! Get total ozone column from GEOS-Chem export variable.
       ! Need to calculate from restart variables on first call!
       !=======================================================================
#if defined( MODEL_GEOS )
       IF ( PHASE /= 1 ) THEN
          CALL CalcTotOzone_( am_I_Root, Input_Opt, PTR_O3, PLE, GCCTROPP, &
                              State_Met%TO3, PTR_GCCTTO3, __RC__ )
       ENDIF
#else
       IF ( calcOzone ) THEN
          IF ( FIRST ) THEN
             CALL CalcTotalOzone( am_I_Root, PLE, GCCTROPP, __RC__ )
          ENDIF
          CALL SetStateMetTO3( am_I_Root, State_Met, __RC__ )
       ENDIF
#endif

       !=======================================================================
       ! Execute GEOS-Chem on multiple PETs
       !=======================================================================
       
       ! Fix negatives!
       ! These can be brought in as an artifact of convection.
       WHERE ( State_Chm%Species < 0.0e0 )
          State_Chm%Species = 1.0e-36
       END WHERE 
       
       ! Execute GEOS-Chem if it's time to run it
       IF ( IsRunTime ) THEN
       
          ! This is mostly for testing
#if defined( MODEL_GEOS )
          IF ( FIRST .AND. Input_Opt%haveImpRst ) THEN
#else
          IF ( FIRST ) THEN
#endif
             IF ( am_I_Root ) THEN
                WRITE(*,*) ''
                WRITE(*,*) 'Doing warm GEOS-Chem restart'
                WRITE(*,*) ''
             ENDIF
          ENDIF
       
#if defined( MODEL_GEOS )
          ! Only if restart file exists...
          IF ( Input_Opt%haveImpRst ) THEN
#endif

#if !defined( MDOEL_GEOS )
             ! Optional memory prints (level >= 2)
             if ( MemDebugLevel > 0 ) THEN
                call ESMF_VMBarrier(vm, RC=STATUS)
                _VERIFY(STATUS)
                call MAPL_MemUtilsWrite(VM, &
                  'Chem_GridCompMod, before chunk_run', RC=STATUS )
                _VERIFY(STATUS)
             endif
#endif
       
             CALL MAPL_TimerOn(STATE, "DO_CHEM")
       
#if !defined( MODEL_GEOS )
             ! NOTE: Second was not extracted previously; set to 0 for now
             second = 0
#endif

             ! Run the GEOS-Chem column chemistry code for the given phase
             CALL GIGC_Chunk_Run( GC         = GC,         & ! Grid comp ref. 
                                  nymd       = nymd,       & ! Current YYYYMMDD
                                  nhms       = nhms,       & ! Current hhmmss
                                  year       = year,       & ! Current year
                                  month      = month,      & ! Current month
                                  day        = day,        & ! Current day
                                  dayOfYr    = dayOfYr,    & ! Current doy
                                  hour       = hour,       & ! Current hour
                                  minute     = minute,     & ! Current minute
                                  second     = second,     & ! Current second
                                  utc        = utc,        & ! Current UTC [hrs]
                                  hElapsed   = hElapsed,   & ! Elapsed hours
                                  Input_Opt  = Input_Opt,  & ! Input Options
                                  State_Chm  = State_Chm,  & ! Chemistry State
                                  State_Diag = State_Diag, & ! Diagnostics State
                                  State_Grid = State_Grid, & ! Grid State
                                  State_Met  = State_Met,  & ! Meteorology State
                                  Phase      = Phase,      & ! Run phase
                                  IsChemTime = IsChemTime, & ! Time for chem?
#if defined( MODEL_GEOS )
                                  FrstRewind = FirstRewind,& ! First rewind?
#endif
                                  __RC__                  )  ! Success or fail?
       
             CALL MAPL_TimerOff(STATE, "DO_CHEM")
       
#if !defined( MODEL_GEOS )
             ! Optional memory prints (level >= 2)
             if ( MemDebugLevel > 0 ) THEN
                call ESMF_VMBarrier(vm, RC=STATUS)
                _VERIFY(STATUS)
                call MAPL_MemUtilsWrite(VM, &
                  'Chem_GridCompMod, after  chunk_run', RC=STATUS )
                _VERIFY(STATUS)
             endif
#endif

#if !defined( MODEL_GEOS )
             where( State_Met%HFLUX .eq. 0.) State_Met%HFLUX = 1e-5
#endif

#if defined( MODEL_GEOS )
          ! Restart file does not exist:
          ELSE
             IF ( am_I_Root ) THEN
                WRITE(*,*) ''
                WRITE(*,*) ' SKIP GEOS-CHEM PHASE ', Phase, &
                           ' BECAUSE IMPORT RESTART FILE IS MISSING'
                WRITE(*,*) ''
             ENDIF
          ENDIF 
#endif
       
       ENDIF !IsRunTime

       !=======================================================================
       ! post-Run method array assignments. This copies the values back from
       ! the State_Chm tracer arrays to the internal state, so that they can
       ! be seen by other components (moist, turbulence, ...)
       !=======================================================================
       
       CALL MAPL_TimerOn(STATE, "CP_AFTR")
#if defined( MODEL_GEOS )
#      include "Includes_After_Run.H"
#endif

#if !defined( MODEL_GEOS )
       State_Chm%Species = State_Chm%Species(:,:,State_Grid%NZ:1:-1,:)
       
       DO I = 1, SIZE(Int2Spc,1)
          IF ( Int2Spc(I)%TrcID <= 0 ) CYCLE
          Int2Spc(I)%Internal = State_Chm%Species(:,:,:,Int2Spc(I)%TrcID)
       ENDDO
#endif

       CALL MAPL_TimerOff(STATE, "CP_AFTR")
       
       ! Stop timer
       ! ----------
       CALL MAPL_TimerOff(STATE, "RUN"  )

#if !defined( MODEL_GEOS )
       ! Fill bundles only on chemistry time steps and after phase 2 
       ! -----------------------------------------------------------
       IF ( IsTendTime ) THEN

          IF ( isProvider ) THEN
             CALL Provider_FillBundles( am_I_Root, tsChem,    PLE, GCCTROPP, &
                                        STATE,     Input_Opt, GC,  EXPORT,   &
                                        __RC__ )
          ENDIF

          IF ( calcOzone ) THEN
             !================================================================
             ! Total ozone and total tropospheric ozone for export [dobsons].
             ! 2.69E+20 per dobson.
             !================================================================
             CALL CalcTotalOzone( am_I_Root, PLE, GCCTROPP, __RC__ )
          ENDIF

       ENDIF ! IsTendTime
#endif

#if defined( MODEL_GEOS )
       !=======================================================================
       ! Perturb O3 by random amount if specified so
       !=======================================================================
       IF ( PHASE /= 1 .AND. ( PerturbO3 .OR. PerturbCO ) ) THEN 
          CALL MAPL_GetPointer( INTSTATE, Ptr3DA, 'TRC_O3', NotFoundOk=.TRUE.,&
                                 __RC__ )
          CALL MAPL_GetPointer( INTSTATE, Ptr3DB, 'TRC_CO', NotFoundOk=.TRUE.,&
                                 __RC__ )
          IF ( ASSOCIATED(Ptr3DA) .OR. ASSOCIATED(Ptr3DB) ) THEN
             DO L=1,State_Grid%NZ
             DO J=1,State_Grid%NY
             DO I=1,State_Grid%NX
                IF ( FIXPERT >= 0.0 ) Rnd(1) = FIXPERT
       
                ! O3
                IF ( PerturbO3 ) THEN
                   IF ( FIXPERT < 0.0 ) THEN
                      CALL RANDOM_NUMBER(Harvest=Rnd)
                      IF ( Rnd(2) >= 0.5 ) Rnd(1) = Rnd(1) * -1.0
                      Rnd(1) = 1.0 + ( Rnd(1) * MAXPERT )
                   ENDIF
                   Ptr3DA(I,J,L) = Ptr3DA(I,J,L) * Rnd(1)
                ENDIF
           
                ! CO
                IF ( PerturbCO ) THEN
                   IF ( FIXPERT < 0.0 ) THEN
                      CALL RANDOM_NUMBER(Harvest=Rnd)
                      IF ( Rnd(2) >= 0.5 ) Rnd(1) = Rnd(1) * -1.0
                      Rnd(1) = 1.0 + ( Rnd(1) * MAXPERT )
                   ENDIF
                   Ptr3DB(I,J,L) = Ptr3DB(I,J,L) * Rnd(1)
                ENDIF
       
             ENDDO
             ENDDO
             ENDDO
          ENDIF
       
          IF ( am_I_Root ) THEN
             IF ( FIXPERT > 0.0 ) THEN
                IF ( PerturbO3 ) write(*,*) &
                     'GCC O3 concentrations perturbed by factor of ', Rnd(1)
                IF ( PerturbCO ) write(*,*) &
                     'GCC CO concentrations perturbed by factor of ', Rnd(1)
             ELSE
                IF ( PerturbO3 ) write(*,*) &
                     'GCC O3 concentrations randomly perturbed'
                IF ( PerturbCO ) write(*,*) &
                     'GCC CO concentrations randomly perturbed'
             ENDIF
          ENDIF
       ENDIF

!      !=======================================================================
!      ! Tendencies (needs more testing - ignore for now)
!      !=======================================================================
!      IF ( IsRunTime ) THEN
!         CALL CalcTendencies_( am_I_Root, Input_Opt, State_Met, State_Chm,   &
!                               IsChemTime, Phase, EXPORT, OX_TEND, H2O_TEND, &
!                               State_Grid%NZ, __RC__ )
!      ENDIF

       ! Archive last active time steps
       pymd = nymd
       phms = nhms
       
#endif

    ENDIF RunningGEOSChem

    !=======================================================================
    ! Diagnostics 
    !=======================================================================

#if defined( MODEL_GEOS )

    !=======================================================================
    ! NOx diagnostics 
    !=======================================================================
    IF ( IsTendTime .AND. Input_Opt%haveImpRst ) THEN 
       CALL NOxDiagnostics_ ( am_I_Root, Input_Opt, State_Met, State_Chm, &
                              State_Diag, IMPORT, EXPORT, IntState, __RC__ )
    ENDIF

    !=======================================================================
    ! Dry volume mixing ratios and PM2.5 diagnostics 
    !=======================================================================
    CALL CalcSpeciesDiagnostics_ ( am_I_Root, Input_Opt, State_Met, State_Chm, &
                                   State_Diag, IMPORT, EXPORT, IntState, &
                                   Q, __RC__ )


    !=======================================================================
    ! Mass-weighted OH
    !=======================================================================
    CALL MassWeightedOH_ ( am_I_Root, Input_Opt, State_Met, State_Chm, &
                           State_Diag, IMPORT, EXPORT, IntState, Q, PLE, TROPP, __RC__ )

    !=======================================================================
    ! Fill ozone export states if GC is the analysis OX provider:
    !      OX: volume mixing ratio
    !      O3: mass mixing ratio
    !  O3PPMV: volume mixing ratio in ppm
    ! OX_TEND: mol mol-1 s-1
    !
    ! GEOS-Chem tracer:
    ! PTR_O3: kg kg-1 total air
    !=======================================================================
     IF ( DoANOX ) THEN
       ASSERT_(ASSOCIATED(PTR_O3))
       IF ( ASSOCIATED(O3     ) ) O3      = PTR_O3
       IF ( ASSOCIATED(O3PPMV ) ) O3PPMV  = PTR_O3  * MAPL_AIRMW / MAPL_O3MW  &
                                          * 1.00E+06
       IF ( ASSOCIATED(OX) ) THEN
          CALL MAPL_GetPointer( INTSTATE, PTR_O3P, 'SPC_O'  , __RC__ )
          CALL MAPL_GetPointer( INTSTATE, PTR_O1D, 'SPC_O1D', __RC__ )
          OX =        PTR_O3  * MAPL_AIRMW / MAPL_O3MW
          OX = OX + ( PTR_O3P * MAPL_AIRMW / OMW )
          OX = OX + ( PTR_O1D * MAPL_AIRMW / OMW )
       ENDIF
    ENDIF

! GEOS-5 (also in gigc_providerservices_mod routine Provider_FillBundles, but 
! might not be in the right place anymore):
    !=======================================================================
    ! Fill RATS export states if GC is the RATS provider
    ! The tracer concentrations of the RATS export states are in mol mol-1.
    !=======================================================================
    IF ( DoRATS) THEN
       IF ( ASSOCIATED(CH4   ) )    CH4 = PTR_CH4    * MAPL_AIRMW /  16.00
       IF ( ASSOCIATED(N2O   ) )    N2O = PTR_N2O    * MAPL_AIRMW /  44.00
       IF ( ASSOCIATED(CFC11 ) )  CFC11 = PTR_CFC11  * MAPL_AIRMW / 137.37
       IF ( ASSOCIATED(CFC12 ) )  CFC12 = PTR_CFC12  * MAPL_AIRMW / 120.91
       IF ( ASSOCIATED(HCFC22) ) HCFC22 = PTR_HCFC22 * MAPL_AIRMW /  86.47
    ENDIF

    !=======================================================================
    ! Total and tropospheric columns 
    !=======================================================================
    CALL CalcColumns_( am_I_Root, Input_Opt, State_Chm, EXPORT,     &
                       INTSTATE, PLE, GCCTROPP, __RC__ )

    !=======================================================================
    ! Total ozone and total tropospheric ozone for export [dobsons].
    ! 2.69E+20 per dobson.
    !=======================================================================
    CALL CalcTotOzone_( am_I_Root, Input_Opt, PTR_O3, PLE, GCCTROPP, &
                        State_Met%TO3, PTR_GCCTTO3, __RC__ )
    IF ( ASSOCIATED(PTR_GCCTO3) ) PTR_GCCTO3 = State_Met%TO3

    ! O3 mass in kg/m2
    CALL MAPL_GetPointer( EXPORT, Ptr3D, 'O3_MASS', NotFoundOk=.TRUE., __RC__ )
    IF ( ASSOCIATED(Ptr3D) ) THEN
       LB = LBOUND(PLE,3)
       DO L=1,State_Grid%NZ
          Ptr3D(:,:,L) = PTR_O3(:,:,L) * ( g0_100 * ( PLE(:,:,L+LB) -  &
                                                      PLE(:,:,L+LB-1) ) )
       ENDDO
    ENDIF
    Ptr3D => NULL()

    !=======================================================================
    ! NO2 columns
    !=======================================================================
    IF ( PHASE /= 1 ) THEN
       CALL MAPL_GetPointer( INTSTATE, PTR_NO2, 'TRC_NO2',      &
                             NotFoundOk=.TRUE., __RC__ )
       CALL MAPL_GetPointer( EXPORT, TNO2, 'NO2_TROPCOLUMN',    &
                             NotFoundOk=.TRUE., __RC__ )
       CALL MAPL_GetPointer( EXPORT, SNO2, 'NO2_TOTCOLUMN',     &
                             NotFoundOk=.TRUE., __RC__ )
       CALL MAPL_GetPointer( EXPORT, FNO2, 'NO2_PBL_FRAC',      &
                             NotFoundOk=.TRUE., __RC__ )

       CALL CalcNO2Column_( am_I_Root, Input_Opt, PTR_NO2, PLE, &
                            PPBL, GCCTROPP, TNO2, SNO2, FNO2, RC )
       _ASSERT(RC==ESMF_SUCCESS,'Error calling CalcNO2Column_')
       PTR_NO2 => NULL()
       TNO2    => NULL()
       SNO2    => NULL()
    ENDIF

    !=======================================================================
    ! Fill AERO bundle if GEOS-Chem is the AERO provider.
    ! For every field of the AERO bundle, we will copy the corresponding
    ! GEOS-Chem tracer field, converting units from mol mol-1 to kg kg-1.
    !=======================================================================
    IF ( DoAERO ) THEN

       ! Get Internal state
       CALL MAPL_Get ( STATE, INTERNAL_ESMF_STATE=INTSTATE, __RC__ ) 

       ! Get AERO bundle
       CALL ESMF_StateGet( EXPORT, 'AERO',     Aero,    __RC__ )
       CALL ESMF_StateGet( Aero,   'AEROSOLS', AeroBdl, __RC__ )

       ! Number of fields in the AERO Bundle
       CALL ESMF_FieldBundleGet ( AeroBdl, FieldCount=nAero, __RC__ )

       ! Update every field
       DO N = 1, nAero

          ! Get field
          CALL ESMF_FieldBundleGet( AeroBdl, N, AeroFld, __RC__ )

          ! Extract GC tracer name, molecular weight and fraction to be used
          CALL ESMF_AttributeGet( AeroFld, NAME='GCNAME', VALUE=GcName, __RC__ )
          CALL ESMF_AttributeGet( AeroFld, NAME='GCMW'  , VALUE=GCMW,   __RC__ )
          CALL ESMF_AttributeGet( AeroFld, NAME='FRAC',   VALUE=FRAC,  __RC__ ) 

          ! Get pointer to Aero data
          CALL ESMF_FieldGet( AeroFld, farrayPtr=AeroPtr3D, __RC__ )

          ! Get pointer to GC data
          CALL MAPL_GetPointer ( INTSTATE, GcPtr3D, TRIM(GcName), __RC__ )

          ! Pass GC to AERO. Convert from mol/mol to kg/kg. Only use the 
          ! fraction specified during initialization (different from 1 for
          ! sea salt aerosols only)
          !AeroPtr3D = GcPtr3D * FRAC * GCMW / MAPL_AIRMW
          AeroPtr3D = GcPtr3D * FRAC

          !!! writing to diagnostics
          GcPtr3D   => NULL()
          CALL ESMF_FieldGet( AeroFld, NAME=GcName, __RC__ )
          CALL MAPL_GetPointer ( EXPORT, GcPtr3D, 'AERO_'//TRIM(GcName), &
                                 NotFoundOk=.TRUE., __RC__ )
          IF ( ASSOCIATED(GcPtr3D) ) GcPtr3D = AeroPtr3D
          !!!

          ! Free pointers
          GcPtr3D   => NULL()
          AeroPtr3D => NULL()
       ENDDO
  
       ! Fill AERO_DP bundle
       CALL FillAeroDP ( am_I_Root, GC, EXPORT, __RC__ )
    ENDIF ! DoAero

    !=======================================================================
    ! Write chemistry top level into CHEMTOP diagnostics array
    !=======================================================================
    IF ( Phase /= 1 ) THEN
       CALL MAPL_GetPointer( EXPORT, Ptr2D, 'CHEMTOP',               &
                             NotFoundOk=.TRUE., __RC__ )
       IF ( ASSOCIATED(Ptr2D) ) THEN
          DO J = 1, State_Grid%NY
          DO I = 1, State_Grid%NX
             Ptr2D(I,J) = State_Grid%NZ - State_Met%ChemGridLev(I,J) + 1
          ENDDO
          ENDDO
       ENDIF
       Ptr2D => NULL()

       ! chemistry tropopause
       CALL MAPL_GetPointer( EXPORT, Ptr2D, 'CHEM_TROPP',            &
                             NotFoundOk=.TRUE., __RC__ )
       IF ( ASSOCIATED(Ptr2D) ) THEN
          Ptr2D(:,:) = State_Met%TROPP(:,:) * 100.0 ! hPa -> Pa 
       ENDIF
       Ptr2D => NULL()
    ENDIF

    IF ( Phase /= 2 ) THEN
       ! convective cloud top height
       CALL MAPL_GetPointer( EXPORT, Ptr2D, 'CONV_CLDTOP',           &
                             NotFoundOk=.TRUE., __RC__ )
       IF ( ASSOCIATED(Ptr2D) ) THEN
          DO J = 1, State_Grid%NY
          DO I = 1, State_Grid%NX
             DO L = State_Grid%NZ,1,-1
                IF ( State_Met%CMFMC(I,J,L) > 0.0d0 ) THEN
                   Ptr2D(I,J) = State_Grid%NZ - L + 1
                   EXIT
                ENDIF
             ENDDO
          ENDDO
          ENDDO
       ENDIF
       Ptr2D => NULL()
    ENDIF

    !=======================================================================
    ! Lightning potential (from GEOS lightning flash rates and convective 
    ! fraction) 
    !=======================================================================
    if ( phase /= 1 ) then
       CALL MAPL_GetPointer( EXPORT, Ptr2D, 'FJX_EXTRAL_NLEVS', &
                             NotFoundOk=.TRUE., __RC__ )
       IF ( ASSOCIATED(Ptr2d) ) Ptr2d = EXTRAL_NLEVS
       CALL MAPL_GetPointer( EXPORT, Ptr2D, 'FJX_EXTRAL_NITER', &
                             NotFoundOk=.TRUE., __RC__ )
       IF ( ASSOCIATED(Ptr2d) ) Ptr2d = EXTRAL_NITER
    endif

    !=======================================================================
    ! Lightning potential (from GEOS lightning flash rates and convective 
    ! fraction) 
    !=======================================================================
    IF ( Phase /= 2 ) THEN
       ! convective cloud top height
       CALL MAPL_GetPointer( EXPORT, Ptr2D, 'LightningPotential', &
                             NotFoundOk=.TRUE., __RC__ )
       IF ( ASSOCIATED(Ptr2D) ) THEN
          _ASSERT(ASSOCIATED(LWI),'LWI is not associated') ! Land-Water-Ice flag
          _ASSERT(ASSOCIATED(LFR),'LFR is not associated') ! GEOS lightning flash rate
          _ASSERT(ASSOCIATED(CNV_FRC),'CNV_FRC is not associated') ! Convective fraction
          CALL MAPL_GetPointer( EXPORT, PtrEmis, 'EMIS_NO_LGHT',  &
                                NotFoundOk=.TRUE., __RC__ )
          Ptr2D = 0.0
          DO J = 1, State_Grid%NY
          DO I = 1, State_Grid%NX
             lp1 = 0.0
             lp2 = 0.0

             ! If there are HEMCO lightning emissions in current grid box set 
             ! lightning potential accordingly 
             IF ( ASSOCIATED(PtrEmis) ) THEN
                IF ( LWI(I,J) == 1 ) THEN
                   lp1 = PtrEmis(I,J) / 1.0e-11 ! Land
                ELSE
                   lp1 = PtrEmis(I,J) / 1.0e-13 ! Water/Ice
                ENDIF
                lp1 = MIN(MAX(0.25,lp1),1.00)
             ENDIF
             
             ! Lightning flash rate
             IF ( LFR(I,J) > 0.0 ) THEN
                IF ( LWI(I,J) == 1 ) THEN
                   lp2 = LFR(I,J) / 5.0e-07 ! Land
                ELSE
                   lp2 = LFR(I,J) / 1.0e-08 ! Water/Ice
                ENDIF
                lp2 = MIN(MAX(0.25,lp2),1.00)

             ! Convective fraction
             ELSE
                lp2 = CNV_FRC(I,J)
             ENDIF

             ! Take highest value
             Ptr2D(I,J) = MAX(lp1,lp2)
          ENDDO
          ENDDO
       ENDIF
       Ptr2D   => NULL()
       PtrEmis => NULL()
    ENDIF 

#else
    !=======================================================================
    ! If we were not doing chemistry, make sure that all tendencies are
    ! zero. We ignore the tendencies that may arise due to physical
    ! processes covered by GEOS-Chem (e.g. convection).
    !=======================================================================
    IF ( .NOT. IsTendTime ) THEN
       CALL Provider_ZeroTendencies( am_I_Root, __RC__ )
    ENDIF
#endif

    !=======================================================================
    ! Copy HISTORY.rc diagnostic data to exports. Includes HEMCO emissions 
    ! diagnostics but excludes internal state and exports created explicitly
    ! in Chem_GridCompMod. NOTE: Exports created explicitly in Chem_GridCompMod
    ! will eventually be moved elsewhere as diagnostics for use with GEOS-5.
    ! (ewl, 11/2/17)
    !=======================================================================
    IF ( FIRST ) THEN
       CALL HistoryExports_SetDataPointers( am_I_Root,     EXPORT,    &
                                            HistoryConfig, State_Chm, &
                                            State_Diag,    State_Met, &
                                            STATUS )
       _VERIFY(STATUS)
    ENDIF
    CALL CopyGCStates2Exports( am_I_Root, Input_Opt, HistoryConfig, STATUS )
    _VERIFY(STATUS)

    !=======================================================================
    ! All done
    !=======================================================================

    IF ( ALLOCATED( zenith     ) ) DEALLOCATE( zenith     )
    IF ( ALLOCATED( solar      ) ) DEALLOCATE( solar      )

    ! Stop timer
    ! ----------
    CALL MAPL_TimerOff(STATE, "TOTAL")
  
#if defined( MODEL_GEOS )
    ! The restart file should exist at least after the first full cycle,
    ! e.g. after phase 1 and phase 2 has been called once.
    IF ( FIRST .AND. Phase == 1 ) THEN
       Input_Opt%haveImpRst = Input_Opt%haveImpRst
    ELSE
       Input_Opt%haveImpRst = .TRUE.
    ENDIF
#endif

    ! Update first flags
    FIRST = .FALSE.

#if defined( MODEL_GEOS )
    ! Unlink HEMCO state from gridcomp objects
    HcoState%GRIDCOMP => NULL()
    HcoState%IMPORT   => NULL() 
    HcoState%EXPORT   => NULL()
#endif

    ! Successful return
    _RETURN(ESMF_SUCCESS)

    ! Formats
100 FORMAT( '---> DATE: ', i4.4, '/', i2.2, '/', i2.2,      &
            '  GMT: ', i2.2, ':', i2.2, '  X-HRS: ', f11.3 )
110 FORMAT( 'Box (',i3,',',i3,') on PET ', i3, ' has coords: ', 2f7.2, &
               ' LocT = ', f9.4 )
200 FORMAT( '### ',                                           / &
            '### ', a, '  |  Execution on PET # ',      i5.5, / &
            '###' )

  END SUBROUTINE Run_
!EOC
!------------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 910.1 and      !
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Finalize_
!
! !DESCRIPTION: Finalize_ is the finalize method of the GEOSCHEMchem gridded 
!  component.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Finalize_( GC, Import, Export, Clock, RC )
!
! !USES:
!
    USE Input_Opt_Mod,         ONLY : OptInput
#if !defined( MODEL_GEOS )
    USE Input_Opt_Mod,         ONLY : Cleanup_Input_Opt
#endif
    USE State_Chm_Mod,         ONLY : ChmState, Cleanup_State_Chm
    USE State_Diag_Mod,        ONLY : DgnState, Cleanup_State_Diag
    USE State_Grid_Mod,        ONLY : GrdState, Cleanup_State_Grid
    USE State_Met_Mod,         ONLY : MetState, Cleanup_State_Met
    USE HCOI_GC_MAIN_MOD,      ONLY : HCOI_GC_FINAL
#if defined( MODEL_GEOS )
    USE HCO_INTERFACE_MOD,     ONLY : HcoState
#endif
!
! !INPUT/OUTPUT PARAMETERS:
!
#if defined( MODEL_GEOS )
    TYPE(ESMF_GridComp), INTENT(INOUT), TARGET :: GC       ! Ref. to this GC
    TYPE(ESMF_State),    INTENT(INOUT), TARGET :: Import   ! Import State
    TYPE(ESMF_State),    INTENT(INOUT), TARGET :: Export   ! Export State
#else
    TYPE(ESMF_GridComp), INTENT(INOUT)         :: GC       ! Ref. to this GC
    TYPE(ESMF_State),    INTENT(INOUT)         :: Import   ! Import State
    TYPE(ESMF_State),    INTENT(INOUT)         :: Export   ! Export State
#endif
    TYPE(ESMF_Clock),    INTENT(INOUT)         :: Clock    ! ESMF Clock object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,             INTENT(OUT)           :: RC       ! Success or failure?
!
! !REMARKS:
!  We call routine Extract_ to return various values (i.e. grid parameters,
!  start & end dates, PET information, etc.) from the ESMF/MAPL environment.  
!  We then pass those to GEOS-Chem via routine GIGC_CHUNK_FINAL, which is
!  located in GEOS-Chem module ./GEOS-Chem/ESMF/gigc_chunk_mod.F90.
!
! !REVISION HISTORY:
!  01 Dec 2009 - A. Da Silva - Initial version
!  08 Apr 2010 - R. Yantosca - Now finalize myState%CF and myState
!  15 Apr 2010 - R. Yantosca - Activate call to GC_CHUNK_FINAL
!  30 Apr 2010 - R. Yantosca - Now use 5 digits for PET
!  02 Jun 2010 - R. Yantosca - Now set Ident%VERBOSE to FALSE
!  09 Oct 2012 - R. Yantosca - Now call MAPL_Am_I_Root to test for root PET
!  08 Mar 2013 - R. Yantosca - Now save the PET # (aka PET #) in Input_Opt
!  15 Mar 2013 - R. Yantosca - Remove IDENT object, which was a holdover from
!                              the GEOS-Chem column code
!  27 Oct 2014 - C. Keller   - Now save species that are not advected into
!                              internal state to ensure they are written into
!                              the restart file.
!  08 May 2015 - C. Keller   - Removed species --> internal copying because
!                              this is now done on every run-time step (GEOS-5)
!  07 Aug 2017 - E. Lundgren - Use species database instead of State_Chm vars 
!                              Spec_Name and Spec_ID
!EOP
!------------------------------------------------------------------------------
!BOC
!
! LOCAL VARIABLES:
!
    ! Objects
    TYPE(ESMF_Grid)            :: Grid        ! ESMF Grid object
    TYPE(ESMF_Config)          :: MaplCF      ! Config (MAPL.rc)
    TYPE(ESMF_Config)          :: GeosCF      ! Config (GEOSCHEM*.rc)
 
    ! Scalars
    LOGICAL                    :: am_I_Root   ! Are we on the root PET?
    CHARACTER(LEN=ESMF_MAXSTR) :: compName    ! Gridded component name
    INTEGER                    :: error       ! GEOS-Chem error code
    INTEGER                    :: myPet       ! # of PET we are on now
    INTEGER                    :: I,  J,  L   ! Loop indices
    REAL                       :: UTC         ! UTC time [hours]
    
    ! Pointers
    TYPE(MAPL_MetaComp), POINTER :: STATE
#if !defined( MODEL_GEOS )
    TYPE(Species),       POINTER :: ThisSpc

    ! For species copying
    INTEGER                     :: IND
    TYPE(ESMF_STATE)            :: INTSTATE
    REAL, POINTER               :: Ptr2D(:,:)      => NULL()
    REAL, POINTER               :: Ptr3D(:,:,:)    => NULL()
    REAL(ESMF_KIND_R8), POINTER :: Ptr2D_R8(:,:)   => NULL()
    REAL(ESMF_KIND_R8), POINTER :: Ptr3D_R8(:,:,:) => NULL()
#endif

    __Iam__('Finalize_')

    !=======================================================================
    ! Initialization
    !=======================================================================

    ! Are we on the root PET
    am_I_Root = MAPL_Am_I_Root()

    ! Set up traceback info
    CALL ESMF_GridCompGet( GC, name=compName, __RC__ )

    ! Identify this routine to MAPL
    Iam = TRIM(compName)//'::Finalize_'

    ! Get my MAPL_Generic state
    ! -------------------------
    CALL MAPL_GetObjectFromGC(GC, STATE, RC=STATUS)
    _VERIFY(STATUS)

    ! Start timers
    ! ------------
!    CALL MAPL_TimerOn(STATE, "TOTAL")
!    CALL MAPL_TimerOn(STATE, "FINALIZE")

    ! Get various parameters from the ESMF/MAPL framework
    CALL Extract_( GC,                 &    ! Ref to this Gridded Component
                   Clock,              &    ! ESMF Clock object
                   Grid     = Grid,    &    ! ESMF Grid object
                   MaplCF   = MaplCF,  &    ! ESMF Config obj (MAPL.rc)
                   GeosCF   = GeosCF,  &    ! ESMF Config obj (GEOSCHEM*.rc)
                   utc      = utc,     &    ! Universal time [hours]
                   localPET = myPet,   &    ! PET # we are on now 
                   __RC__ )

#if !defined( MODEL_GEOS )
    !=========================================================================
    ! Archive species in internal state. Also include certain State_Chm
    ! and State_Met arrays needed for continuity across runs (ewl, 12/13/18)
    !=========================================================================

    ! Get Internal state
    CALL MAPL_Get ( STATE, INTERNAL_ESMF_STATE=INTSTATE, __RC__ ) 

    ! Loop over all species
    DO I = 1, State_Chm%nSpecies
 
       ! Get info about this species from the species database
       ThisSpc => State_Chm%SpcData(I)%Info

       ! Skip if empty
       IF ( TRIM(ThisSpc%Name) == '' ) CYCLE

       ! Is this a tracer?
       IND = IND_( TRIM(ThisSpc%Name) )
       IF ( IND >= 0 ) CYCLE

       ! Get data from internal state and copy to species array
       CALL MAPL_GetPointer( INTSTATE, Ptr3D_R8, TRIM(ThisSpc%Name), &
                             notFoundOK=.TRUE., __RC__ )
       IF ( .NOT. ASSOCIATED(Ptr3D_R8) ) CYCLE
       Ptr3D_R8 = State_Chm%Species(:,:,State_Grid%NZ:1:-1,IND)
       Ptr3D_R8 => NULL()

       ! Verbose 
       if ( MAPL_am_I_Root()) write(*,*)                &
                'Species written to INTERNAL state: ',  &
                TRIM(ThisSpc%Name)
    ENDDO

    CALL MAPL_GetPointer( INTSTATE, Ptr3d_R8, 'H2O2AfterChem',  &
                          notFoundOK=.TRUE., __RC__ ) 
    IF ( ASSOCIATED(Ptr3d_R8) .AND. &
         ASSOCIATED(State_Chm%H2O2AfterChem) ) THEN
       Ptr3d_R8(:,:,State_Grid%NZ:1:-1) = State_Chm%H2O2AfterChem
    ENDIF
    Ptr3d_R8 => NULL()
    
    CALL MAPL_GetPointer( INTSTATE, Ptr3d_R8, 'SO2AfterChem',   &
                          notFoundOK=.TRUE., __RC__ ) 
    IF ( ASSOCIATED(Ptr3d_R8) .AND. &
         ASSOCIATED(State_Chm%SO2AfterChem) ) THEN
       Ptr3d_R8(:,:,State_Grid%NZ:1:-1) = State_Chm%SO2AfterChem
    ENDIF
    Ptr3d_R8 => NULL()
    
    CALL MAPL_GetPointer( INTSTATE, Ptr2d_R8, 'DryDepNitrogen', &
                          notFoundOK=.TRUE., __RC__ ) 
    IF ( ASSOCIATED(Ptr2d_R8) .AND. &
         ASSOCIATED(State_Chm%DryDepNitrogen) ) THEN
       Ptr2d_R8 = State_Chm%DryDepNitrogen
    ENDIF
    Ptr2d_R8 => NULL()
    
    CALL MAPL_GetPointer( INTSTATE, Ptr2d_R8, 'WetDepNitrogen', &
                          notFoundOK=.TRUE., __RC__ ) 
    IF ( ASSOCIATED(Ptr2d_R8) .AND. &
         ASSOCIATED(State_Chm%WetDepNitrogen) ) THEN
       Ptr2d_R8 = State_Chm%WetDepNitrogen
    ENDIF
    Ptr2d_R8 => NULL()
    
    CALL MAPL_GetPointer( INTSTATE, Ptr3d_R8, 'KPPHvalue' ,     &
                          notFoundOK=.TRUE., __RC__ ) 
    IF ( ASSOCIATED(Ptr3d_R8) .AND. &
         ASSOCIATED(State_Chm%KPPHvalue) ) THEN
       Ptr3d_R8(:,:,1:State_Grid%NZ-State_Grid%MaxChemLev) = 0.0
       Ptr3d_R8(:,:,State_Grid%NZ:State_Grid%NZ-State_Grid%MaxChemLev+1:-1) = &
          State_Chm%KPPHvalue(:,:,1:State_Grid%MaxChemLev)
    ENDIF
    Ptr3d_R8 => NULL()

    CALL MAPL_GetPointer( INTSTATE, Ptr3d_R8, 'DELP_DRY' ,     &
                          notFoundOK=.TRUE., __RC__ ) 
    IF ( ASSOCIATED(Ptr3d_R8) .AND. &
         ASSOCIATED(State_Met%DELP_DRY) ) THEN
       Ptr3d_R8(:,:,State_Grid%NZ:1:-1) =  &
                 State_Met%DELP_DRY(:,:,1:State_Grid%NZ)
    ENDIF
    Ptr3d_R8 => NULL()
#endif

#if defined( MODEL_GEOS )
    !=======================================================================
    ! Print end-of-simulation output
    !=======================================================================

    ! GlobalSum used within Print_Mean_OH does not work anymore - skip mean
    ! OH printout for now (ckeller, 8/27/2019).
    !! Print mean OH value to GEOS-Chem log-file
    !IF ( Input_Opt%ITS_A_FULLCHEM_SIM ) THEN
    !   CALL Print_Mean_OH( GC, State_Grid, logLun, __RC__ )
    !ENDIF
#endif

    !=======================================================================
    ! Finalize the Gridded Component
    !=======================================================================

#if defined( MODEL_GEOS )
    ! Link HEMCO state to gridcomp objects
    _ASSERT(ASSOCIATED(HcoState),'HcoState is not associated')
    HcoState%GRIDCOMP => GC
    HcoState%IMPORT   => IMPORT
    HcoState%EXPORT   => EXPORT
#endif

    ! Finalize HEMCO
    CALL HCOI_GC_FINAL( .FALSE., RC )
    IF ( Input_Opt%AmIRoot ) THEN
       IF ( RC == GC_SUCCESS ) THEN
          write(*,'(a)') 'HEMCO::Finalize... OK.'
       ELSE
          write(*,'(a)') 'HEMCO::Finalize... FAILURE.'
       ENDIF
    ENDIF

    ! Deallocate fields of the Chemistry State object
    CALL Cleanup_State_Chm( State_Chm, RC )
    IF ( Input_Opt%AmIRoot ) THEN
       IF ( RC == GC_SUCCESS ) THEN
          write(*,'(a)') 'Chem::State_Chm Finalize... OK.'
       ELSE
          write(*,'(a)') 'Chem::State_Chm Finalize... FAILURE.'
       ENDIF
    ENDIF

    ! Deallocate fields of the Diagnostics State object
    CALL Cleanup_State_Diag( State_Diag, RC )
    IF ( Input_Opt%AmIRoot ) THEN
       IF ( RC == GC_SUCCESS ) THEN
          write(*,'(a)') 'Chem::State_Diag Finalize... OK.'
       ELSE
          write(*,'(a)') 'Chem::State_Diag Finalize... FAILURE.'
       ENDIF
    ENDIF

    ! Deallocate fields of the Grid State object
    CALL Cleanup_State_Grid( State_Grid, RC )
    IF ( Input_Opt%AmIRoot ) THEN
       IF ( RC == GC_SUCCESS ) THEN
          write(*,'(a)') 'Chem::State_Grid Finalize... OK.'
       ELSE
          write(*,'(a)') 'Chem::State_Grid Finalize... FAILURE.'
       ENDIF
    ENDIF

    ! Deallocate fields of the Meteorology State object
    CALL Cleanup_State_Met( State_Met, RC )
    IF ( Input_Opt%AmIRoot ) THEN
       IF ( RC == GC_SUCCESS ) THEN
          write(*,'(a)') 'Chem::State_Met Finalize... OK.'
       ELSE
          write(*,'(a)') 'Chem::State_Met Finalize... FAILURE.'
       ENDIF
    ENDIF

#if !defined( MODEL_GEOS )
    ! Deallocate fields of the Input Options object
    ! The call to Cleanup_Input_Opt causes a memory leak error. Comment
    ! for now (ckeller, 11/29/16).
    ! Does this still cause a memory leak? (ewl, 12/14/18)
     CALL Cleanup_Input_Opt( Input_Opt, RC )
    IF ( Input_Opt%AmIRoot ) THEN
       IF ( RC == GC_SUCCESS ) THEN
          write(*,'(a)') 'Chem::Input_Opt Finalize... OK.'
       ELSE
          write(*,'(a)') 'Chem::Input_Opt Finalize... FAILURE.'
       ENDIF
    ENDIF
#endif

    ! Free Int2Spc pointer
    IF ( ASSOCIATED(Int2Spc) ) THEN
       DO I=1,SIZE(Int2Spc,1)
          Int2Spc(I)%Internal => NULL()
       ENDDO
       DEALLOCATE(Int2Spc)
    ENDIF

    ! Deallocate the history interface between GC States and ESMF Exports
    CALL Destroy_HistoryConfig( am_I_Root, HistoryConfig, RC )

#if defined( MODEL_GEOS )
    ! Free local pointers
    !O3               => NULL()
    !O3PPMV           => NULL()
    !OX               => NULL()
    !OX_TEND          => NULL()
    !CH4              => NULL()
    !N2O              => NULL()
    !CFC11            => NULL()
    !CFC12            => NULL()
    !HCFC22           => NULL()
    !H2O_TEND         => NULL()
    !PTR_O3           => NULL()
    !PTR_CH4          => NULL()
    !PTR_N2O          => NULL()
    !PTR_CFC11        => NULL()
    !PTR_CFC12        => NULL()
    !PTR_HCFC22       => NULL()
    !PTR_H2O          => NULL()
    !PTR_GCCTO3       => NULL()
    !PTR_GCCTTO3      => NULL()
#else
    ! Deallocate provide pointers and arrays
    CALL Provider_Finalize( am_I_Root, __RC__ )
#endif

    ! Finalize MAPL Generic
    CALL MAPL_GenericFinalize( GC, Import, Export, Clock, __RC__ )

    !=======================================================================
    ! All done
    !=======================================================================

    ! Stop timers
    ! -----------
!    CALL MAPL_TimerOff(STATE, "FINALIZE")
!    CALL MAPL_TimerOff(STATE, "TOTAL")

    ! Successful return
    _RETURN(ESMF_SUCCESS)

  END SUBROUTINE Finalize_
!EOC
!------------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 910.1 and      !
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Extract_
!
! !DESCRIPTION: GC routine extracts several common quantities from the 
!  ESMF/MAPL environment so that they can be later passed down to the 
!  grid-independent GEOS-Chem code.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Extract_( GC,         Clock,    Grid,    MaplCF, GeosCF,    &
                       localPet,   petCount,                             &
                       IM,         JM,       LM,                         &     
                       IM_WORLD,   JM_WORLD, LM_WORLD,                   &
                       lonCtr,     latCtr,   advCount,                   &
                       nymdB,      nymdE,    nymd,    nhmsB,  nhmsE,     &
                       nhms,       year,     month,   day,    dayOfYr,   &
                       hour,       minute,   second,  utc,    hElapsed,  &
                       tsChem,     tsDyn,    mpiComm, ZTH,   SLR,        &
#if defined( MODEL_GEOS )
                       haveImpRst,                                       &
#endif
                       RC )
!
! !INPUT PARAMETERS:
!
    TYPE(ESMF_Clock),    INTENT(IN)            :: Clock       ! ESMF clock obj 
!                                                             
! !INPUT/OUTPUT PARAMETERS:                                   
!                                                             
    TYPE(ESMF_GridComp), INTENT(INOUT)         :: GC          ! GC grid comp
!                                                             
! !OUTPUT PARAMETERS:                                         
!              
    !----------------------------------
    ! ESMF and/or MAPL quantities
    !----------------------------------
    TYPE(ESMF_Grid),     INTENT(OUT), OPTIONAL :: Grid        ! ESMF Grid obj
    TYPE(ESMF_Config),   INTENT(OUT), OPTIONAL :: MaplCF      ! AGCM.rc
    TYPE(ESMF_Config),   INTENT(OUT), OPTIONAL :: GeosCF      ! GEOSCHEM*.rc
    INTEGER,             INTENT(OUT), OPTIONAL :: localPet    ! This PET
    INTEGER,             INTENT(OUT), OPTIONAL :: petCount    ! Total # of PETs
    INTEGER,             INTENT(OUT), OPTIONAL :: mpiComm     ! MPI Comm Handle

    !----------------------------------
    ! Local grid coordinates 
    ! (defined on the current CPU)
    !----------------------------------
    INTEGER,             INTENT(OUT), OPTIONAL :: IM          ! Total # lons
    INTEGER,             INTENT(OUT), OPTIONAL :: JM          ! Total # lats
    INTEGER,             INTENT(OUT), OPTIONAL :: LM          ! Total # levs
    REAL(ESMF_KIND_R4),  POINTER,     OPTIONAL :: lonCtr(:,:) ! Lon ctrs [rad]
    REAL(ESMF_KIND_R4),  POINTER,     OPTIONAL :: latCtr(:,:) ! Lat ctrs [rad]

    !----------------------------------
    ! Global grid coordinates
    !----------------------------------
    INTEGER,             INTENT(OUT), OPTIONAL :: IM_WORLD    ! Global # lons
    INTEGER,             INTENT(OUT), OPTIONAL :: JM_WORLD    ! Global # lats
    INTEGER,             INTENT(OUT), OPTIONAL :: LM_WORLD    ! Global # levs
                                                              
    !----------------------------------
    ! Date and time variables
    !----------------------------------
   INTEGER(ESMF_KIND_I8),INTENT(OUT), OPTIONAL :: advCount    ! # of clock advs
    INTEGER,             INTENT(OUT), OPTIONAL :: nymdB       ! YYYYMMDD @ start
    INTEGER,             INTENT(OUT), OPTIONAL :: nymdE       ! YYYYMMDD @ end
    INTEGER,             INTENT(OUT), OPTIONAL :: nymd        ! YYYYMMDD now
    INTEGER,             INTENT(OUT), OPTIONAL :: nhmsB       ! hhmmss @ start
    INTEGER,             INTENT(OUT), OPTIONAL :: nhmsE       ! hhmmss @ end
    INTEGER,             INTENT(OUT), OPTIONAL :: nhms        ! hhmmss now
    INTEGER,             INTENT(OUT), OPTIONAL :: year        ! UTC year 
    INTEGER,             INTENT(OUT), OPTIONAL :: month       ! UTC month
    INTEGER,             INTENT(OUT), OPTIONAL :: day         ! UTC day
    INTEGER,             INTENT(OUT), OPTIONAL :: dayOfYr     ! UTC day of year
    INTEGER,             INTENT(OUT), OPTIONAL :: hour        ! UTC hour
    INTEGER,             INTENT(OUT), OPTIONAL :: minute      ! UTC minute
    INTEGER,             INTENT(OUT), OPTIONAL :: second      ! UTC second
    REAL,                INTENT(OUT), OPTIONAL :: utc         ! UTC time [hrs]
    REAL,                INTENT(OUT), OPTIONAL :: hElapsed    ! Elapsed hours

    !-----------------------------------                     
    ! Timestep variables [seconds]          
    !-----------------------------------                     
    REAL,                INTENT(OUT), OPTIONAL :: tsChem      ! Chemistry
    REAL,                INTENT(OUT), OPTIONAL :: tsDyn       ! Dynamics

    !-----------------------------------                     
    ! Solar parameters
    !-----------------------------------                     
    REAL,                INTENT(OUT), OPTIONAL :: ZTH(:,:)    ! Solar zth angle
    REAL,                INTENT(OUT), OPTIONAL :: SLR(:,:)    ! Insolation

#if defined( MODEL_GEOS )
    !-----------------------------------------------
    ! Optional import restart file existence inquiry
    !-----------------------------------------------
    LOGICAL,             INTENT(OUT), OPTIONAL :: haveImpRst ! Import rst exist?
#endif

    !-----------------------------------                        
    ! Return code 
    !-----------------------------------                     
    INTEGER,             INTENT(OUT), OPTIONAL :: RC          ! 0 = all is well
!
! !REMARKS:
!  If you need to obtain a quantity not returned by this routine, you can
!  manually extract it from the MaplCF or GeosCF configuration objects.
!
! !REVISION HISTORY:
!  01 Dec 2009 - A. Da Silva - Initial version
!  07 Apr 2010 - R. Yantosca - Added ProTeX headers
!  08 Apr 2010 - R. Yantosca - Make all outputs optional
!  08 Apr 2010 - R. Yantosca - Added outputs for localPet, petCount
!  08 Apr 2010 - R. Yantosca - Added outputs for individual time values
!                              as well as elapsed time (hours)
!  13 Apr 2010 - R. Yantosca - Now take tsDyn from the MAPL "RUN_DT:" setting
!  30 Nov 2012 - R. Yantosca - Now return IM_WORLD, JM_WORLD, LM_WORLD
!  30 Nov 2012 - R. Yantosca - Now return local indices I_LO, J_LO, I_HI, J_HI
!  05 Dec 2012 - R. Yantosca - Removed latEdg argument; cosmetic changes
!  13 Feb 2013 - E. Nielsen  - Restart file inquiry for GEOS-5
!  05 Jan 2016 - S. D. Eastham - Fixed order of time calls
!  02 Nov 2017 - E. Lundgren - Replace use of local GridGetInterior with 
!                              call to MAPL_GridGetInterior (now public)
!EOP
!------------------------------------------------------------------------------
!BOC
!
! LOCAL VARIABLES:
! 
    ! Objects
    TYPE(ESMF_Time)               :: startTime      ! ESMF start time obj
    TYPE(ESMF_Time)               :: stopTime       ! ESMF stop time obj
    TYPE(ESMF_Time)               :: currTime       ! ESMF current time obj
    TYPE(ESMF_TimeInterval)       :: elapsedTime    ! ESMF elapsed time obj
    TYPE(ESMF_TimeInterval)       :: chemInterval   ! chemistry interval
    TYPE(ESMF_ALARM)              :: ALARM          ! Run alarm 
    TYPE(ESMF_VM)                 :: VM             ! ESMF VM object
    TYPE(GEOSCHEM_State), POINTER :: myState        ! Legacy state
    TYPE(GEOSCHEM_Wrap)           :: wrap           ! Wrapper for myState
    TYPE(MAPL_MetaComp),  POINTER :: STATE          ! MAPL MetaComp object
    TYPE(MAPL_SunOrbit)           :: sunOrbit

    ! Scalars
    CHARACTER(len=ESMF_MAXSTR)    :: compName       ! Gridded component name
    CHARACTER(len=ESMF_MAXSTR)    :: importRstFN    ! Import restart file name
    INTEGER(ESMF_KIND_I8)         :: count          ! # of clock advances
    INTEGER                       :: locDims(3)     ! Array for local dims
    INTEGER                       :: globDims(3)    ! Array for global dims
    INTEGER                       :: doy            ! Day of year (0-365/366)
    INTEGER                       :: yyyy, mm, dd   ! Year, month, day
    INTEGER                       :: h,    m,  s    ! Hour, minute, seconds
    INTEGER                       :: IL,   IU       ! Min/max local lon indices
    INTEGER                       :: JL,   JU       ! Min/max local lat indices
    REAL                          :: elapsedHours   ! Elapsed hours of run
    REAL(ESMF_KIND_R8)            :: dt_r8          ! chemistry timestep

    __Iam__('Extract_')

    !=======================================================================
    ! Initialization
    !=======================================================================

    ! Get my name and set-up traceback handle
    CALL ESMF_GridCompGet( GC, name=compName, vm=VM, __RC__ )
    Iam = TRIM(compName)//'::Extract_'

    ! Get the internal state which holds the private Config object
    CALL ESMF_UserCompGetInternalState( GC, 'GEOSCHEM_State', wrap, STATUS )
    _VERIFY(STATUS)
    myState => wrap%ptr

    ! Get generic state object
    CALL MAPL_GetObjectFromGC( GC, STATE, __RC__ )

    ! Assume successful return
    IF ( PRESENT( RC ) ) RC = ESMF_SUCCESS

    ! Zero variables
    locDims  = 0
    globDims = 0
    IL       = 0
    JL       = 0
    IU       = 0
    JU       = 0

    !=======================================================================
    ! Extract information from ESMF VM object
    !=======================================================================


    ! Index of the PET we are on now
    IF ( PRESENT( localPet ) ) THEN
       CALL ESMF_VmGet( VM, localPet=localPet, __RC__ )
    ENDIF

    ! Total # of PETs used by this gridded component
    IF ( PRESENT( petCount ) ) THEN
       CALL ESMF_VmGet( VM, petCount=petCount, __RC__ )
    ENDIF

    ! Global MPI Communicator Handle
    IF ( PRESENT( mpiComm ) ) THEN
       CALL ESMF_VmGet( VM, mpicommunicator=mpiComm, __RC__ )
    ENDIF

    !=======================================================================
    ! Extract information from ESMF Config objects
    !=======================================================================

    ! Get the Config object 
    CALL ESMF_GridCompGet( GC, Config=MaplCF, __RC__ )
    
    ! Get the Config object based on "GEOSCHEMchem_GridComp.rc"
    GeosCF = myState%myCF

    ! Dynamic timestep (in seconds)
    IF ( PRESENT( tsDyn ) ) THEN
       CALL ESMF_ConfigGetAttribute( MaplCF, tsDyn, Default=1800.,        &
                                     Label="RUN_DT:",             __RC__ )
    ENDIF

    ! Chemistry timestep (in seconds)
    IF ( PRESENT( tsChem ) ) THEN
        CALL MAPL_Get( STATE, RUNALARM=ALARM, __RC__ ) 
        CALL ESMF_AlarmGet( ALARM, RingInterval=chemInterval, __RC__ )
        CALL ESMF_TimeIntervalGet( chemInterval, s_r8=dt_r8, __RC__ )
        tsChem = real(dt_r8)

        IF(tsChem < tsDyn) THEN
           IF( MAPL_AM_I_ROOT() ) THEN
#if defined( MODEL_GEOS )
              WRITE(6,*) 'GEOSCHEMCHEM_DT cannot be less than RUN_DT'
#else
              WRITE(6,*) 'Chem_DT cannot be less than RUN_DT'
#endif
           ENDIF
           STATUS = 1
           _VERIFY(STATUS)
        ENDIF
    ENDIF

#if defined( MODEL_GEOS )
    ! Simulation dates. Legacy stuff, not used. 
    ! Set to dummy values (ckeller, 1/18/18)
    !IF ( PRESENT( nymdB ) ) nymdB = 20130701 
    !IF ( PRESENT( nhmsB ) ) nhmsB = 000000 
    !IF ( PRESENT( nymdE ) ) nymdE = 20130701 
    !IF ( PRESENT( nhmsE ) ) nhmsE = 000000 

    !=======================================================================
    ! Does the import restart file exist?
    !=======================================================================
    
    ! Import restart file name
    IF ( PRESENT ( haveImpRst ) ) THEN
       CALL ESMF_ConfigGetAttribute( MaplCF, importRstFN,                      &
                                  DEFAULT = "geoschemchem_import_rst.nc4",     &
                                  LABEL = "GEOSCHEMCHEM_IMPORT_RESTART_FILE:", &
                                  __RC__ )

       ! remove bootstrap parameter
       IF( importRstFN(1:1) == '-' .OR. &  
           importRstFN(1:1) == '+'       ) THEN
          importRstFN = importRstFN(2:LEN(TRIM(importRstFN)))
       ENDIF
       INQUIRE( FILE=TRIM( importRstFN ), EXIST=haveImpRst )
       IF( MAPL_AM_I_ROOT() ) THEN
          PRINT *," ",TRIM( importRstFN )," exists: ", haveImpRst
          PRINT *," "
       END IF
    ENDIF
#endif

    !=======================================================================
    ! Extract time/date information
    !=======================================================================
    
    ! Get the ESMF time object
    CALL ESMF_ClockGet( Clock,                    &
                        startTime    = startTime, &
                        stopTime     = stopTime,  &
                        currTime     = currTime,  &
                        advanceCount = count,     &
                         __RC__ )

    ! Get starting-time fields from the time object
    CALL ESMF_TimeGet( startTime, yy=yyyy, mm=mm, dd=dd, dayOfYear=doy, &
                                  h=h,     m=m,   s=s,   __RC__ )

    ! Save fields for return
    IF ( PRESENT( nymd     ) ) CALL MAPL_PackTime( nymd, yyyy, mm, dd )
    IF ( PRESENT( nhms     ) ) CALL MAPL_PackTime( nhms, h,    m,  s  )

    ! Get ending-time fields from the time object
    CALL ESMF_TimeGet( stopTime, yy=yyyy, mm=mm, dd=dd, dayOfYear=doy, &
                                 h=h,     m=m,   s=s,   __RC__ )

    ! Save packed fields for return
    IF ( PRESENT( nymdE    ) ) CALL MAPL_PackTime( nymdE, yyyy, mm, dd )
    IF ( PRESENT( nhmsE    ) ) CALL MAPL_PackTime( nhmsE, h,    m,  s  )

    IF ( PRESENT( advCount ) ) advCount = count
 
    !=======================================================================
    ! SDE 2017-01-05: The following calls must be kept as a single block,
    ! or the wrong date/time elements will be returned (the yyyy/mm/dd    
    ! etc variables are re-used). Specifically, the output variables must
    ! be set now, before the variables are re-used.
    !=======================================================================
    ! Start of current-time block
    !=======================================================================
    ! Get current-time fields from the time object
    CALL ESMF_TimeGet( currTime, yy=yyyy, mm=mm, dd=dd, dayOfYear=doy, &
                                 h=h,     m=m,   s=s,   __RC__ )
 
    ! Save packed fields for return
    IF ( PRESENT( nymd     ) ) CALL MAPL_PackTime( nymd, yyyy, mm, dd )
    IF ( PRESENT( nhms     ) ) CALL MAPL_PackTime( nhms, h,    m,  s  )

    ! Save the various extacted current-time fields for return
    IF ( PRESENT( year     ) ) year     = yyyy
    IF ( PRESENT( month    ) ) month    = mm
    IF ( PRESENT( day      ) ) day      = dd
    IF ( PRESENT( dayOfYr  ) ) dayOfYr  = doy
    IF ( PRESENT( hour     ) ) hour     = h
    IF ( PRESENT( minute   ) ) minute   = m
    IF ( PRESENT( second   ) ) second   = s
    IF ( PRESENT( utc      ) ) utc      = ( DBLE( h )        ) + & 
                                          ( DBLE( m )/60d0   ) + &
                                          ( DBLE( s )/3600d0 )

    !=======================================================================
    ! End of current-time block
    !=======================================================================

    CALL ESMF_TimeGet( startTime, yy=yyyy, mm=mm, dd=dd, dayOfYear=doy, &
                                 h=h,     m=m,   s=s,   __RC__ )

    ! Save fields for return
    IF ( PRESENT( nymdB    ) ) CALL MAPL_PackTime( nymdB, yyyy, mm, dd )
    IF ( PRESENT( nhmsB    ) ) CALL MAPL_PackTime( nhmsB, h,    m,  s  )

    CALL ESMF_TimeGet( stopTime, yy=yyyy, mm=mm, dd=dd, dayOfYear=doy, &
                                 h=h,     m=m,   s=s,   __RC__ )

    ! Save fields for return
    IF ( PRESENT( nymdE    ) ) CALL MAPL_PackTime( nymdE, yyyy, mm, dd )
    IF ( PRESENT( nhmsE    ) ) CALL MAPL_PackTime( nhmsE, h,    m,  s  )

    IF ( PRESENT( advCount ) ) advCount = count

    ! Compute elapsed time since start of simulation
    elapsedTime = currTime - startTime

    ! Get time fields from the elapsedTime object
    CALL ESMF_TimeIntervalGet( elapsedTime, h=h, m=m, s=s, __RC__ )

    ! Convert to decimal hours
    elapsedHours = DBLE( h ) + ( DBLE( m )/60d0 ) + ( DBLE( s )/3600d0 )
    
    ! Save fields for return
    IF ( PRESENT( hElapsed ) ) hElapsed = elapsedHours

    !=======================================================================
    ! Extract grid information
    !=======================================================================
    IF ( PRESENT( Grid ) ) THEN
    
       ! Get the ESMF grid attached to this gridded component
       CALL ESMF_GridCompGet( GC, grid=Grid, __RC__ )

       ! Get # of dimensions on this pet, and globally
       CALL MAPL_GridGet( Grid,                                        &
                          localCellCountPerDim  = locDims,             &
                          globalCellCountPerDim = globDims,            &
                          __RC__ )
          
#if defined( MODEL_GEOS )
       ! Get the upper and lower bounds of on each PET
       CALL GridGetInterior( Grid, IL, IU, JL, JU, __RC__  )
#else
       ! Get the upper and lower bounds of on each PET using MAPL
       CALL MAPL_GridGetInterior( Grid, IL, IU, JL, JU )
#endif

    ENDIF

    ! Save fields for return
    IF ( PRESENT( IM       ) ) IM       = locDims(1)
    IF ( PRESENT( JM       ) ) JM       = locDims(2)
    IF ( PRESENT( LM       ) ) LM       = locDims(3)
    IF ( PRESENT( IM_WORLD ) ) IM_WORLD = globDims(1)
    IF ( PRESENT( JM_WORLD ) ) JM_WORLD = globDims(2)
    IF ( PRESENT( LM_WORLD ) ) LM_WORLD = globDims(3)

    ! Longitude values on this PET
    IF ( PRESENT( lonCtr ) ) THEN
       CALL MAPL_Get( STATE, lons=lonCtr, __RC__ )
    ENDIF

    ! Latitude values on this PET
    IF ( PRESENT( latCtr ) ) THEN
       CALL MAPL_Get( STATE, lats=latCtr, __RC__ )
    ENDIF

    !=======================================================================
    ! Get solar zenith angle enformation
    !=======================================================================
    IF ( PRESENT( ZTH    ) .and. PRESENT( SLR    )  .and. &
         PRESENT( lonCtr ) .and. PRESENT( latCtr ) ) THEN
         
       ! Get the Orbit object (of type MAPL_SunOrbit),
       ! which is used in the call to MAPL_SunGetInsolation
       CALL MAPL_Get( STATE,                       &
                      LONS      = lonCtr,             &
                      LATS      = latCtr,             &
                      ORBIT     = sunOrbit,           &
                      __RC__                         )

       ! Get the solar zenith angle and solar insolation
       ! NOTE: ZTH, SLR are allocated outside of this routine
       CALL MAPL_SunGetInsolation( LONS  = lonCtr,    &
                                   LATS  = latCtr,    &
                                   ORBIT = sunOrbit,  &
                                   ZTH   = ZTH,       &
                                   SLR   = SLR,       &
                                   CLOCK = Clock,     &
                                   __RC__            )

    ENDIF

    !=======================================================================
    ! All done
    !=======================================================================
    _RETURN(ESMF_SUCCESS)

  END SUBROUTINE Extract_
!EOC
#if defined( MODEL_GEOS )
!------------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 910.1 and      !
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: CalcTotOzone_ 
!
! !DESCRIPTION: CalcTotOzone_ calculates total ozone for the entire
!  atmosphere and troposphere only (in dobsons) and writes them into
!  the export variables GCCTO3 and GCCTTO3, respectively. Expects O3 input
!  in kg/kg total.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE CalcTotOzone_ ( am_I_Root, Input_Opt, O3, PLE, TROPP, TO3, &
                             TTO3, RC )
!
! !USES:
!
    USE Precision_Mod
!
! !INPUT PARAMETERS:
!
    LOGICAL                          :: am_I_Root
     TYPE(OptInput)                  :: Input_Opt 
    REAL,     POINTER                :: O3   (:,:,:)
    REAL,     POINTER                :: PLE  (:,:,:)
    REAL,     POINTER                :: TROPP(:,:  )
!                                                             
! !INPUT/OUTPUT PARAMETERS:                                         
!              
    REAL(fp), POINTER                :: TO3 (:,:)
    REAL,     POINTER                :: TTO3(:,:)
!                                                             
! !OUTPUT PARAMETERS:                                         
!              
    INTEGER, INTENT(OUT), OPTIONAL   :: RC
!
! !REVISION HISTORY:
!  25 Oct 2014 - C. Keller   - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! LOCAL VARIABLES:
! 
    REAL,  ALLOCATABLE           :: DUsLayerL(:,:)! Dobsons in a layer, 
                                                  !  for total ozone
    REAL,  ALLOCATABLE           :: wgt(:,:)      ! Layer thickness weighting
                                                  !  for total ozone
    REAL                         :: const
    INTEGER                      :: IM, JM, LM, LB, L, STATUS
    CHARACTER(LEN=ESMF_MAXSTR)   :: Iam

    !=======================================================================
    ! CalcTotOzone_ begins here
    !=======================================================================

    ! Traceback handle
    Iam = 'CalcTotOzone_'

    ! Nothing to do if it's a Rn simulation
    IF ( Input_Opt%ITS_A_RnPbBe_SIM ) THEN
       RC = ESMF_SUCCESS
       RETURN
    ENDIF

    ! Nothing to do if neither of the arrays is associated
    IF ( .NOT. ASSOCIATED(TO3) .AND. .NOT. ASSOCIATED(TTO3) ) THEN
       RC = ESMF_SUCCESS
       RETURN
    ENDIF

    ! Make sure O3 is defined
    _ASSERT(ASSOCIATED(O3),'O3 is not associated')

    ! Grid size
    IM = SIZE(O3,1)
    JM = SIZE(O3,2)
    LM = SIZE(O3,3)

    ! Lower bound of PLE 3rd dim
    LB = LBOUND(PLE,3)

    ! Reset values 
    IF ( ASSOCIATED(TO3 ) ) TO3  = 0.0
    IF ( ASSOCIATED(TTO3) ) TTO3 = 0.0

    ! Allocate local variables
    ALLOCATE(DUsLayerL(IM,JM), STAT=STATUS)
    _VERIFY(STATUS)
    ALLOCATE(wgt(IM,JM), STAT=STATUS)
    _VERIFY(STATUS)

    ! constant 
    const = 0.01 * MAPL_AVOGAD / ( MAPL_GRAV * (MAPL_AIRMW/1000.0) )
    const = const * MAPL_AIRMW / MAPL_O3MW ! convert kg/kg total to v/v total
 
    ! Calculate total ozone
    DO L = 1,LM

!       DUsLayerL(:,:) = AIRDENS(:,:,L) / (MAPL_AIRMW/1000.0) * MAPL_AVOGAD &
!                        * O3(:,:,L) / 2.69e20

       DUsLayerL(:,:) = O3(:,:,L) * ((PLE(:,:,L+LB)-PLE(:,:,L+LB-1))/100.0) &
                        * const / 2.69e16 / 1000.0
 
!       DUsLayerL(:,:) = O3(:,:,L)*(PLE(:,:,L)-PLE(:,:,L-1))               &
!                        *(MAPL_AVOGAD/2.69E+20)/(MAPL_AIRMW*MAPL_GRAV)

       IF ( ASSOCIATED(TO3) ) TO3 = TO3+DUsLayerL
       IF ( ASSOCIATED(TTO3) ) THEN
          wgt  = MAX(0.0,MIN(1.0,(PLE(:,:,L+LB)-TROPP(:,:)) &
                 /(PLE(:,:,L+LB)-PLE(:,:,L+LB-1))))
          TTO3 = TTO3+DUsLayerL*wgt
       END IF
    END DO
 
    ! Cleanup
    DEALLOCATE(DUsLayerL, STAT=STATUS)
    _VERIFY(STATUS)
    DEALLOCATE(wgt, STAT=STATUS)
    _VERIFY(STATUS)

    ! Successful return
    RC = ESMF_SUCCESS

  END SUBROUTINE CalcTotOzone_
!EOC
!------------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 910.1 and      !
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: CalcNO2Column_ 
!
! !DESCRIPTION: CalcNO2Column_ calculates total NO2 column (troposphere and
!  stratosphere. Expects NO2 input in kg/kg total air. 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE CalcNO2Column_ ( am_I_Root, Input_Opt, NO2, PLE, PPBL,  &
                              TROPP, TNO2, SNO2, FNO2, RC, omi, LAYERS )
!
! !INPUT PARAMETERS:
!
    LOGICAL                          :: am_I_Root
     TYPE(OptInput)                  :: Input_Opt
    REAL,     POINTER                :: NO2  (:,:,:)
    REAL,     POINTER                :: PLE  (:,:,:)
    REAL,     POINTER                :: PPBL (:,:)
    REAL,     POINTER                :: TROPP(:,:)
    LOGICAL,  INTENT(IN), OPTIONAL   :: omi
!                                                             
! !INPUT/OUTPUT PARAMETERS:                                         
!              
    REAL,     POINTER                :: TNO2(:,:)
    REAL,     POINTER                :: SNO2(:,:)
    REAL,     POINTER                :: FNO2(:,:)
!                                                             
! !OUTPUT PARAMETERS:                                         
!              
    REAL,    POINTER,     OPTIONAL   :: LAYERS(:,:,:)
    INTEGER, INTENT(OUT), OPTIONAL   :: RC
!
! !REVISION HISTORY:
!  25 Jul 2016 - C. Keller   - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! LOCAL VARIABLES:
! 
    REAL,  ALLOCATABLE           :: DUsLayerL(:,:)! NO2 in a layer
    REAL,  ALLOCATABLE           :: wgt(:,:)      ! Layer thickness weighting
    REAL,  ALLOCATABLE           :: PBLTNO2(:,:)  ! NO2 PBL column
    REAL                         :: DP, const, iwgt
    INTEGER                      :: IM, JM, LM, LB, I, J, L, STATUS
    CHARACTER(LEN=ESMF_MAXSTR)   :: Iam

    LOGICAL                      :: domi

    !=======================================================================
    ! CalcNO2Column_ begins here
    !=======================================================================

    ! Traceback handle
    Iam = 'CalcNO2Column_'

    domi = .FALSE.
    IF ( PRESENT(omi) ) domi=omi

    ! Nothing to do if it's a Rn simulation
    IF ( Input_Opt%ITS_A_RnPbBe_SIM ) THEN
       RC = ESMF_SUCCESS
       RETURN
    ENDIF

    ! Nothing to do if neither of the arrays is associated
    IF ( .NOT. ASSOCIATED(TNO2) .AND. .NOT. ASSOCIATED(SNO2) &
       .AND. .NOT. ASSOCIATED(FNO2) ) THEN
       RC = ESMF_SUCCESS
       RETURN
    ENDIF

    ! Make sure O3 is defined
    _ASSERT(ASSOCIATED(NO2),'NO2 is not associated')

    ! Grid size
    IM = SIZE(NO2,1)
    JM = SIZE(NO2,2)
    LM = SIZE(NO2,3)

    ! Lower bound of PLE 3rd dim
    LB = LBOUND(PLE,3)

    ! Reset values 
    IF ( ASSOCIATED(TNO2) ) TNO2 = 0.0
    IF ( ASSOCIATED(SNO2) ) SNO2 = 0.0
    IF ( ASSOCIATED(FNO2) ) FNO2 = 0.0
    IF ( PRESENT(LAYERS)  ) LAYERS = 0.0

    ! Allocate local variables
    ALLOCATE(DUsLayerL(IM,JM), STAT=STATUS)
    _VERIFY(STATUS)
    ALLOCATE(wgt(IM,JM), STAT=STATUS)
    _VERIFY(STATUS)
    IF ( ASSOCIATED(FNO2) ) THEN
       _ASSERT(ASSOCIATED(PPBL),'PPBL is not associated')
       ALLOCATE(PBLTNO2(IM,JM), STAT=STATUS)
       _VERIFY(STATUS)
       PBLTNO2(:,:) = 0.0
    ENDIF

    ! constant. Note: MAPL_AVOGAD is kmol-1, MAPL_AIRMW is kg kmol-1, 
    ! MAPL_GRAV is m s-2.
    const = MAPL_AVOGAD / 1000.0 / ( MAPL_GRAV * MAPL_AIRMW / 1000.0 )
    const = const * MAPL_AIRMW / 46.0 ! convert kg/kg total to v/v total

    ! Calculate total NOx
    DO L = 1,LM
    DO J = 1,JM
    DO I = 1,IM

       DUsLayerL(I,J) = NO2(I,J,L) * ( (PLE(I,J,L+LB)-PLE(I,J,L+LB-1)) ) &
                        * const

       ! rescale: molec/m2 --> molec/cm2
       DUsLayerL(I,J) = DUsLayerL(I,J) / 1.0e4

       ! compute 1e15 molec/cm2 so that numbers are in better range
       DUsLayerL(I,J) = DUsLayerL(I,J) / 1.0e15

       IF ( ASSOCIATED(SNO2) ) THEN
          SNO2(I,J) = SNO2(I,J)+DUsLayerL(I,J)  !*wgt
       END IF

       wgt(I,J)  = MAX(0.0,MIN(1.0,(PLE(I,J,L)-TROPP(I,J))/(PLE(I,J,L) &
                   -PLE(I,J,L-1))))

       IF ( ASSOCIATED(TNO2) ) THEN
          TNO2(I,J) = TNO2(I,J)+DUsLayerL(I,J)*wgt(I,J)
       END IF

       IF ( PRESENT(LAYERS) ) THEN
          LAYERS(I,J,L) = DUsLayerL(I,J)*1.0e15*wgt(I,J)
       END IF

       ! Compute fraction of surface layer relative to total trop. layer
       IF ( ASSOCIATED(FNO2) ) THEN
          FNO2(I,J)    = FNO2(I,J) + DUsLayerL(I,J) * wgt(I,J)
          iwgt         = MAX(0.0,MIN(1.0,(PLE(I,J,L)-PPBL(I,J)) &
                         /(PLE(I,J,L)-PLE(I,J,L-1))))
          PBLTNO2(I,J) = PBLTNO2(I,J) + DUsLayerL(I,J) * iwgt
          IF ( L == LM ) THEN
             FNO2(I,J) = PBLTNO2(I,J) / FNO2(I,J)
             !FNO2(I,J) = ( DUsLayerL(I,J)*wgt(I,J) ) / FNO2(I,J)
          END IF
       END IF

    END DO
    END DO
    END DO

    ! Cleanup
    DEALLOCATE(DUsLayerL, STAT=STATUS)
    _VERIFY(STATUS)
    DEALLOCATE(wgt, STAT=STATUS)
    _VERIFY(STATUS)
    IF ( ALLOCATED(PBLTNO2) ) DEALLOCATE(PBLTNO2)
    ! Successful return
    RC = ESMF_SUCCESS

  END SUBROUTINE CalcNO2Column_
!EOC
!------------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 910.1 and      !
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: CalcColumns_ 
!
! !DESCRIPTION: CalcColumns_ calculates total and tropospheric columns for a
!  number of species. 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE CalcColumns_ ( am_I_Root, Input_Opt, State_Chm, EXPORT, &
                            INTSTATE, PLE, TROPP, RC )
!
! !INPUT PARAMETERS:
!
    LOGICAL                                 :: am_I_Root
    TYPE(OptInput)                          :: Input_Opt 
    TYPE(ChmState)                          :: State_Chm 
    TYPE(ESMF_State), INTENT(INOUT), TARGET :: Export
    TYPE(ESMF_STATE)                        :: INTSTATE 
    REAL,     POINTER                       :: PLE  (:,:,:)
    REAL,     POINTER                       :: TROPP(:,:  )
!                                                             
! !INPUT/OUTPUT PARAMETERS:                                         
!              
!                                                             
! !OUTPUT PARAMETERS:                                         
!              
    INTEGER, INTENT(OUT), OPTIONAL   :: RC
!
! !REVISION HISTORY:
!  25 Oct 2014 - C. Keller   - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! LOCAL VARIABLES:
! 
    REAL,  POINTER               :: ExpTOTCOL(:,:)
    REAL,  POINTER               :: ExpTRPCOL(:,:)
    REAL,  POINTER               :: IntSpc   (:,:,:)
    REAL,  ALLOCATABLE           :: DUsLayerL(:,:)! Dobsons in a layer, 
                                                  !  for total ozone
    REAL,  ALLOCATABLE           :: wgt(:,:)      ! Layer thickness weighting
                                                  !  for total ozone
    REAL                         :: MW, const
    INTEGER                      :: I, ID, IM, JM, LM, LB, L, STATUS
    CHARACTER(LEN=ESMF_MAXSTR)   :: Iam
    CHARACTER(LEN=15)            :: ISPEC

    !=======================================================================
    ! CalcColumns_ begins here
    !=======================================================================

    ! Traceback handle
    Iam = 'CalcColumns_'

    ! Grid size
    IM = SIZE(PLE,1)
    JM = SIZE(PLE,2)
    LM = SIZE(PLE,3)-1
   
    ! Allocate local variables
    ALLOCATE(DUsLayerL(IM,JM), STAT=STATUS)
    _VERIFY(STATUS)
    ALLOCATE(wgt(IM,JM), STAT=STATUS)
    _VERIFY(STATUS)
   
    ! Do for all species
    DO I=1,SIZE(COLLIST,1)

       ! For convenience
       ISPEC = COLLIST(I)

       ! Get diagnostics from export state
       CALL MAPL_GetPointer ( EXPORT, ExpTOTCOL,  'TOTCOL_'//TRIM(ISPEC), &
                              NotFoundOK=.TRUE., __RC__ )
       CALL MAPL_GetPointer ( EXPORT, ExpTRPCOL, 'TROPCOL_'//TRIM(ISPEC), &
                              NotFoundOK=.TRUE., __RC__ )
   
       ! Nothing to do if neither of the arrays is associated
       IF ( .NOT. ASSOCIATED(ExpTOTCOL) .AND. .NOT. ASSOCIATED(ExpTRPCOL) ) &
           CYCLE 
   
       ! Get species from internal state
       CALL MAPL_GetPointer ( INTSTATE, IntSpc, 'TRC_'//TRIM(ISPEC), __RC__ )
   
       ! Lower bound of PLE 3rd dim
       LB = LBOUND(PLE,3)
   
       ! Reset values 
       IF ( ASSOCIATED(ExpTOTCOL) ) ExpTOTCOL = 0.0
       IF ( ASSOCIATED(ExpTRPCOL) ) ExpTRPCOL = 0.0
   
       ! Get molecular weight of species
       ID = IND_(TRIM(ISPEC))
       MW = State_Chm%SpcData(ID)%Info%EmMW_g 
   
       ! constant 
       !const = MAPL_AVOGAD / 1000.0 / ( MAPL_GRAV * ( MAPL_AIRMW / 1000.0 ) )
       !const = const * MAPL_AIRMW  / MW 
       const = MAPL_AVOGAD / ( MAPL_GRAV * MW )
    
       ! Calculate total and trop. column
       DO L = 1,LM
          DUsLayerL(:,:) = IntSpc(:,:,L) * ( PLE(:,:,L+LB) &
                           - PLE(:,:,L+LB-1) ) * const 
          ! rescale: molec/m2 --> molec/cm2
          DUsLayerL(:,:) = DUsLayerL(:,:) / 1.0e4 / 1.0e15
          IF ( ASSOCIATED(ExpTOTCOL) ) ExpTOTCOL = ExpTOTCOL+DUsLayerL
          IF ( ASSOCIATED(ExpTRPCOL) ) THEN
             wgt  = MAX(0.0,MIN(1.0,(PLE(:,:,L+LB)-TROPP(:,:)) &
                    /(PLE(:,:,L+LB)-PLE(:,:,L+LB-1))))
             ExpTRPCOL = ExpTRPCOL+DUsLayerL*wgt
          END IF
       END DO
    ENDDO   
 
    ! Cleanup
    DEALLOCATE(DUsLayerL, STAT=STATUS)
    _VERIFY(STATUS)
    DEALLOCATE(wgt, STAT=STATUS)
    _VERIFY(STATUS)

    ! Successful return
    RC = ESMF_SUCCESS

  END SUBROUTINE CalcColumns_
!EOC
!------------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 910.1 and      !
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: CalcTendencies_ 
!
! !DESCRIPTION: CalcTendencies_ computes tendencies. 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE CalcTendencies_( am_I_Root, Input_Opt, State_Met, State_Chm, &
                              IsChemTime, Phase, EXPORT, OX_TEND, H2O_TEND, &
                              LM, RC )
!
! !USES:
!
    USE TENDENCIES_MOD,          ONLY : Tend_Get
!
! !INPUT/OUTPUT PARAMETERS:
!
    LOGICAL,             INTENT(IN)            :: am_I_Root
    TYPE(OptInput),      INTENT(INOUT)         :: Input_Opt 
    TYPE(MetState),      INTENT(INOUT)         :: State_Met 
    TYPE(ChmState),      INTENT(INOUT)         :: State_Chm 
    LOGICAL,             INTENT(IN)            :: IsChemTime
    INTEGER,             INTENT(IN)            :: Phase    
    TYPE(ESMF_State),    INTENT(INOUT)         :: Export   ! Export State
    REAL,                POINTER               :: OX_TEND(:,:,:)
    REAL,                POINTER               :: H2O_TEND(:,:,:)
    INTEGER,             INTENT(IN)            :: LM
!
! !OUTPUT PARAMETERS:
!
    INTEGER,             INTENT(OUT)           :: RC       ! Success or failure?
!
! !REMARKS:
!
! !REVISION HISTORY:
!  05 Dec 2017 - C. Keller   - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! LOCAL VARIABLES:
!
    ! Objects
 
    ! Scalars
    INTEGER                    :: STATUS
    INTEGER                    :: I, N, IND
    INTEGER                    :: Stage
    CHARACTER(LEN=ESMF_MAXSTR) :: Iam           ! Gridded component name
    CHARACTER(LEN=ESMF_MAXSTR) :: FldName
    REAL, POINTER              :: Ptr3D(:,:,:), Ptr2D(:,:)
    REAL(f4), POINTER          :: Tend(:,:,:)

    
    !=======================================================================
    ! Initialization
    !=======================================================================

    ! Identify this routine to MAPL
    Iam = 'GCC::CalcTendencies_'

    ! O3 tendencies due to chemistry.
    IF ( IsChemTime .AND. Phase==2 ) THEN
       CALL MAPL_GetPointer( EXPORT, Ptr3D, 'MTEND_CHEM_O3', &
                             NotFoundOk=.TRUE., __RC__ )
       IF ( ASSOCIATED(OX_TEND) .OR. ASSOCIATED(Ptr3D) ) THEN

          ! Get O3 tracer ID
          IND = Ind_('O3')
          _ASSERT(IND>0,'O3 tracer ID is negative')

          ! Get tracer tendency in kg/kg dry air/s 
          CALL Tend_Get( am_I_Root, Input_Opt, 'CHEM', IND, Stage, &
                         Tend, __RC__ )
          ! Fill OX_TEND (in kg/kg/s)
          IF ( ASSOCIATED(Tend) .AND. Stage==2 ) THEN
             OX_TEND = Tend(:,:,LM:1:-1)
          ELSE
             OX_TEND = 0.0
          ENDIF

          ! export as kg/m2/s:
          IF ( ASSOCIATED(Ptr3D) .AND. ASSOCIATED(Tend) .AND. Stage==2 ) THEN
             Ptr3D(:,:,LM:1:-1) = Tend(:,:,1:LM) * ( g0_100 &
                                  * State_Met%DELP_DRY(:,:,LM:1:-1) )
          ENDIF

          ! Cleanup
          Ptr3D => NULL()
          Tend  => NULL()
       ENDIF
    ENDIF

    ! H2O tendencies due to chemistry.
    IF ( ASSOCIATED(H2O_TEND) .AND. IsChemTime .AND. Phase==2 ) THEN
       ! Initialize
       H2O_TEND = 0.0
       ! Get tracer ID
       IND = Ind_('H2O')
       IF ( IND > 0 ) THEN
          ! Get tendency due to chemistry in kg/kg/s
          CALL Tend_Get( am_I_Root, Input_Opt, 'CHEM', IND, Stage, &
                         Tend, __RC__ )
          ! Fill H2O_TEND. Convert to kg/kg/s
          IF ( ASSOCIATED(Tend) .AND. Stage==2 ) THEN
             H2O_TEND = Tend(:,:,LM:1:-1)
          ENDIF

!          ! It looks like H2O_TEND is not preserved in the MAPL restart files, 
!          ! which causes the time regression to fail. Force GCC H2O tendencies
!          ! to zero until this problem is resolved! (ckeller, 11/02/2015) 
!          H2O_TEND = 0.0
!          if(am_I_Root) write(*,*) 'GCC: H2O_TEND set to zero!'
       ENDIF
       ! Cleanup
       Tend => NULL()
    ENDIF

    ! Drydep mass tendencies
    ! MTEND_FLUX is a 2D diagnostics and reflects the mass flux due to dry 
    ! deposition
    IF ( PHASE == 1 ) THEN 
       DO I = 1, State_Chm%nDryDep
          N = State_Chm%Map_DryDep(I)
          FldName = 'MTEND_FLUX_'//TRIM(State_Chm%SpcData(N)%Info%Name)
          CALL MAPL_GetPointer( Export, Ptr2D, TRIM(FldName), &
                                NotFoundOk=.TRUE., __RC__ )
          IF ( ASSOCIATED(Ptr2D) ) THEN
             ! Tracer index 
             IND = Ind_(TRIM(State_Chm%SpcData(N)%Info%Name))
             _ASSERT(IND>0,'Tracer ID is negative')
             ! Get tracer tendency in kg/kg dry air/s 
             CALL Tend_Get( am_I_Root, Input_Opt, 'FLUX', IND, Stage, &
                            Tend, __RC__ )
             ! Export as kg/m2/s:
             IF ( ASSOCIATED(Tend) .AND. Stage==2 ) THEN
                Tend(:,:,1:LM) = Tend(:,:,1:LM) * ( g0_100  &
                                 * State_Met%DELP_DRY(:,:,LM:1:-1) )
                Ptr2D(:,:) = SUM(Tend,DIM=3) 
             ENDIF
             Ptr2D => NULL()
          END IF
       ENDDO
    ENDIF

    ! Wetdep mass tendencies
    ! MTEND_WETD is a 2D diagnostics and reflects the mass flux due to 
    ! large-scale wet deposition
    IF ( PHASE == 2 ) THEN 
       DO I = 1, State_Chm%nWetDep
          N = State_Chm%Map_WetDep(I)
          FldName = 'MTEND_WETD_'//TRIM(State_Chm%SpcData(N)%Info%Name)
          CALL MAPL_GetPointer( Export, Ptr2D, TRIM(FldName), &
                                NotFoundOk=.TRUE., __RC__ )
          IF ( ASSOCIATED(Ptr2D) ) THEN
             ! Tracer index 
             IND = Ind_(TRIM(State_Chm%SpcData(N)%Info%Name))
             _ASSERT(IND>0,'Tracer ID is negative')
             ! Get tracer tendency in kg/kg dry air/s 
             CALL Tend_Get( am_I_Root, Input_Opt, 'WETD', IND, Stage, &
                            Tend, __RC__ )
             ! Export as kg/m2/s:
             IF ( ASSOCIATED(Tend) .AND. Stage==2 ) THEN
                Tend(:,:,1:LM) = Tend(:,:,1:LM) * ( g0_100 &
                                 * State_Met%DELP_DRY(:,:,LM:1:-1) )
                Ptr2D(:,:) = SUM(Tend,DIM=3) 
             ENDIF
             Ptr2D => NULL()
          END IF
       ENDDO
    ENDIF

    ! Convection mass tendencies
    ! MTEND_CONV is a 2D diagnostics and reflects the mass flux due to 
    ! convective wet deposition
    IF ( PHASE == 1 ) THEN 
       DO I = 1, State_Chm%nWetDep
          N = State_Chm%Map_WetDep(I)
          FldName = 'MTEND_CONV_'//TRIM(State_Chm%SpcData(N)%Info%Name)
          CALL MAPL_GetPointer( Export, Ptr2D, TRIM(FldName),  &
                                NotFoundOk=.TRUE., __RC__ )
          IF ( ASSOCIATED(Ptr2D) ) THEN
             ! Tracer index 
             IND = Ind_(TRIM(State_Chm%SpcData(N)%Info%Name))
             _ASSERT(IND>0,'Tracer ID is negative')
             ! Get tracer tendency in kg/kg dry air/s 
             CALL Tend_Get( am_I_Root, Input_Opt, 'CONV', IND, Stage, &
                            Tend, __RC__ )
             ! Export as kg/m2/s:
             IF ( ASSOCIATED(Tend) .AND. Stage==2 ) THEN
                Tend(:,:,1:LM) = Tend(:,:,1:LM) * ( g0_100  &
                                 * State_Met%DELP_DRY(:,:,LM:1:-1) )
                Ptr2D(:,:) = SUM(Tend,DIM=3) 
             ENDIF
             Ptr2D => NULL()
          END IF
       ENDDO
    ENDIF

    ! Successful return
    _RETURN(ESMF_SUCCESS)

  END SUBROUTINE CalcTendencies_ 
!EOC
!------------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 910.1 and      !
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: NOxDiagnostics_
!
! !DESCRIPTION: NOxDiagnostics_ computes some NOx diagnostics 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE NOxDiagnostics_( am_I_Root, Input_Opt, State_Met, State_Chm, &
                              State_Diag, IMPORT, EXPORT, INTSTATE, RC )
!
! !USES:
!
!
! !INPUT/OUTPUT PARAMETERS:
!
    LOGICAL,             INTENT(IN)            :: am_I_Root
    TYPE(OptInput),      INTENT(INOUT)         :: Input_Opt 
    TYPE(MetState),      INTENT(INOUT)         :: State_Met 
    TYPE(ChmState),      INTENT(INOUT)         :: State_Chm 
    TYPE(DgnState),      INTENT(INOUT)         :: State_Diag
    TYPE(ESMF_State),    INTENT(INOUT)         :: Import   ! Import State
    TYPE(ESMF_State),    INTENT(INOUT)         :: Export   ! Export State
    TYPE(ESMF_STATE),    INTENT(INOUT)         :: INTSTATE
!
! !OUTPUT PARAMETERS:
!
    INTEGER,             INTENT(INOUT)         :: RC       ! Success or failure?
!
! !REMARKS:
!
! !REVISION HISTORY:
!  05 Dec 2017 - C. Keller   - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! LOCAL VARIABLES:
!
    ! Objects
 
    ! Scalars
    INTEGER                    :: STATUS
    INTEGER                    :: I, IM, JM, LM
    CHARACTER(LEN=ESMF_MAXSTR) :: Iam           ! Gridded component name

    ! NO2 to NOx ratio & NOx lifetime
    REAL, POINTER              :: NO2toNOx(:,:,:), NOxtau(:,:,:) 
    REAL, ALLOCATABLE          :: JNO2(:,:,:)
    REAL, ALLOCATABLE          :: O3eff(:,:,:)
    REAL, ALLOCATABLE          :: k1O3(:,:,:), k3OH(:,:,:)
    REAL, ALLOCATABLE          :: tmp3d(:,:,:)
    INTEGER                    :: idJNO2, idk1, idk3

    !=======================================================================
    ! Routine starts here 
    !=======================================================================

    ! Identify this routine to MAPL
    Iam = 'GCC::NOxDiagnostics_'

    !=======================================================================
    ! NOx diagnostics 
    !=======================================================================
    CALL MAPL_GetPointer( EXPORT, NO2toNOx, 'NO2toNOx', &
                          NotFoundOk=.TRUE., __RC__ )
    CALL MAPL_GetPointer( EXPORT, NOxtau,   'NOx_tau' , &
                          NotFoundOk=.TRUE., __RC__ )
    IF ( ASSOCIATED(NO2toNOx) .OR. ASSOCIATED(NOxtau) ) THEN
       ! Get j-value index
       idJNO2 = -1
       DO I=1,Input_Opt%NN_JVals
          IF ( Input_Opt%JVal_IDs(I) == id_jno2 ) THEN
             idJNO2 = I
             EXIT
          ENDIF 
       ENDDO
       ! Get reaction coeffs indeces 
       idk1 = -1
       idk3 = -1
       DO I=1,Input_Opt%NN_RxnRconst
          IF ( Input_Opt%RxnRconst_IDs(I) == id_rc_no ) THEN
             idk1 = I
          ENDIF 
          IF ( Input_Opt%RxnRconst_IDs(I) == id_rc_no2 ) THEN
             idk3 = I
          ENDIF 
       ENDDO
       _ASSERT( idJNO2 > 0,'idJNO2 is negative' )
       _ASSERT( idk1   > 0,'idk1 is negative' )
       _ASSERT( ASSOCIATED(State_Diag%O3concAfterChem),'O3concAfterChem is not associated'  )
       _ASSERT( ASSOCIATED(State_Diag%RO2concAfterChem),'RO2concAfterChem is not associated' )
       ! Allocate local arrays
       IM = SIZE(State_Diag%O3concAfterChem,1)
       JM = SIZE(State_Diag%O3concAfterChem,2)
       LM = SIZE(State_Diag%O3concAfterChem,3)
       ALLOCATE(JNO2(IM,JM,LM))
       ALLOCATE(O3eff(IM,JM,LM)) 
       ALLOCATE(k1O3(IM,JM,LM))
       ALLOCATE(k3OH(IM,JM,LM))
       ALLOCATE(tmp3d(IM,JM,LM))
       JNO2  = 0.0
       O3eff = 0.0
       k1O3  = 0.0
       k3OH  = 0.0
       tmp3d = 0.0

       jNO2(:,:,:)  = State_Diag%JValIndiv(:,:,LM:1:-1,idJNO2)
       O3eff(:,:,:) = State_Diag%O3concAfterChem (:,:,LM:1:-1) + &
                      State_Diag%RO2concAfterChem(:,:,LM:1:-1)
       k1O3         = State_Diag%RxnRconst(:,:,LM:1:-1,idk1) &
                    * O3eff

       ! Compute NO2 to NOx ratio
       IF ( ASSOCIATED(NO2toNOx) ) THEN
          tmp3d    = jNO2 + k1O3
          where ( tmp3d /= 0.0 ) 
             NO2toNOx = k1O3 / tmp3d
          else where
             NO2toNOx = 0.0
          end where 
       ENDIF
       ! Compute NOx chemical lifetime
       IF ( ASSOCIATED(NOxtau) ) THEN
          _ASSERT( idk3 > 0,'idk3 is negative' )
          _ASSERT( ASSOCIATED(State_Diag%OHconcAfterChem),'OHconcAfterChem is negative' )
          k3OH      = State_Diag%RxnRconst(:,:,LM:1:-1,idk3) &
                    * State_Diag%OHconcAfterChem(:,:,LM:1:-1)
          tmp3d = 0.0
          where ( k1O3 /= 0.0 ) 
             tmp3d = ( jNO2 + k1O3 ) / k1O3 
          else where
             tmp3d = 0.0
          end where
          where ( k3OH /= 0.0 )
             NOxtau = tmp3d / k3OH
          else where
             NOxtau = 0.0
          end where 
       ENDIF

       ! Cleanup
       DEALLOCATE( jNO2, O3eff, k1O3, k3OH, tmp3d )
    ENDIF 

    !=======================================================================
    ! All done 
    !=======================================================================

    _RETURN(ESMF_SUCCESS)

  END SUBROUTINE NOxDiagnostics_
!EOC
!------------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 910.1 and      !
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: CalcSpeciesDiagnostics_
!
! !DESCRIPTION: CalcSpeciesDiagnostics_ computes species' diagnostics 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE CalcSpeciesDiagnostics_( am_I_Root, Input_Opt, State_Met, &
                                      State_Chm, State_Diag, IMPORT, EXPORT, &
                                      INTSTATE, Q, RC )
!
! !USES:
!
    USE TENDENCIES_MOD,          ONLY : Tend_Get
!
! !INPUT/OUTPUT PARAMETERS:
!
    LOGICAL,             INTENT(IN)            :: am_I_Root
    TYPE(OptInput),      INTENT(INOUT)         :: Input_Opt 
    TYPE(MetState),      INTENT(INOUT)         :: State_Met 
    TYPE(ChmState),      INTENT(INOUT)         :: State_Chm 
    TYPE(DgnState),      INTENT(INOUT)         :: State_Diag
    TYPE(ESMF_State),    INTENT(INOUT)         :: Import   ! Import State
    TYPE(ESMF_State),    INTENT(INOUT)         :: Export   ! Export State
    TYPE(ESMF_STATE),    INTENT(INOUT)         :: INTSTATE
    REAL,                POINTER               :: Q(:,:,:)
!
! !OUTPUT PARAMETERS:
!
    INTEGER,             INTENT(INOUT)         :: RC       ! Success or failure?
!
! !REMARKS:
!
! !REVISION HISTORY:
!  05 Dec 2017 - C. Keller   - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! LOCAL VARIABLES:
!
    ! Objects
 
    ! Scalars
    INTEGER                    :: STATUS
    INTEGER                    :: I, J, N, IM, JM, LM, DryID
    LOGICAL                    :: IsBry, IsNOy,  IsCly, IsOrgCl, DoPM25
    LOGICAL                    :: RunMe, IsPM25, IsSOA, IsNi, IsSu, IsOC, IsBC, IsDu, IsSS
    CHARACTER(LEN=ESMF_MAXSTR) :: Iam           ! Gridded component name
    CHARACTER(LEN=ESMF_MAXSTR) :: FieldName, SpcName
    REAL                       :: MW
    REAL                       :: Ra_z0, Ra_2m, Ra_10m
    REAL                       :: t_, rho_, fhs_, ustar_, dz_, z0h_
    REAL                       :: BrCoeff, ClCoeff, OrgClCoeff
!    REAL, ALLOCATABLE          :: Ra2m(:,:), Ra10m(:,:)
!    REAL, ALLOCATABLE          :: Lz0(:,:),  L2m(:,:), L10m(:,:)
    REAL, POINTER              :: Ptr3D(:,:,:), PtrTmp(:,:,:), Ptr2D(:,:)
    REAL, POINTER              :: Bry(:,:,:), NOy(:,:,:), Ptr2m(:,:), Ptr10m(:,:)
    REAL, POINTER              :: Cly(:,:,:), OrgCl(:,:,:) 
    REAL, ALLOCATABLE          :: CONV(:,:), MyPM25(:,:,:)
    TYPE(Species), POINTER     :: SpcInfo

    ! For PM25 diagnostics
    REAL, POINTER, DIMENSION(:,:,:)  :: PtrPM25
    REAL, POINTER, DIMENSION(:,:)    :: PtrPM25_2m
    REAL, POINTER, DIMENSION(:,:)    :: PtrPM25_10m
    REAL, POINTER, DIMENSION(:,:,:)  :: PtrPM25_SOA
    REAL, POINTER, DIMENSION(:,:,:)  :: PtrPM25_nitrates
    REAL, POINTER, DIMENSION(:,:,:)  :: PtrPM25_sulfates
    REAL, POINTER, DIMENSION(:,:,:)  :: PtrPM25_OC
    REAL, POINTER, DIMENSION(:,:,:)  :: PtrPM25_BC
    REAL, POINTER, DIMENSION(:,:,:)  :: PtrPM25_dust
    REAL, POINTER, DIMENSION(:,:,:)  :: PtrPM25_seasalt
    REAL, POINTER, DIMENSION(:,:,:)  :: PtrNH4
    REAL, POINTER, DIMENSION(:,:,:)  :: PtrNIT
    REAL, POINTER, DIMENSION(:,:,:)  :: PtrSO4
    REAL, POINTER, DIMENSION(:,:,:)  :: PtrBCPI
    REAL, POINTER, DIMENSION(:,:,:)  :: PtrBCPO
    REAL, POINTER, DIMENSION(:,:,:)  :: PtrOCPI
    REAL, POINTER, DIMENSION(:,:,:)  :: PtrOCPO
    REAL, POINTER, DIMENSION(:,:,:)  :: PtrDST1
    REAL, POINTER, DIMENSION(:,:,:)  :: PtrDST2
    REAL, POINTER, DIMENSION(:,:,:)  :: PtrSALA
    REAL, POINTER, DIMENSION(:,:,:)  :: PtrTSOA0
    REAL, POINTER, DIMENSION(:,:,:)  :: PtrTSOA1
    REAL, POINTER, DIMENSION(:,:,:)  :: PtrTSOA2
    REAL, POINTER, DIMENSION(:,:,:)  :: PtrTSOA3
    REAL, POINTER, DIMENSION(:,:,:)  :: PtrISOA1
    REAL, POINTER, DIMENSION(:,:,:)  :: PtrISOA2
    REAL, POINTER, DIMENSION(:,:,:)  :: PtrISOA3
    REAL, POINTER, DIMENSION(:,:,:)  :: PtrASOAN
    REAL, POINTER, DIMENSION(:,:,:)  :: PtrASOA1
    REAL, POINTER, DIMENSION(:,:,:)  :: PtrASOA2
    REAL, POINTER, DIMENSION(:,:,:)  :: PtrASOA3


    REAL, PARAMETER                  :: Lv = 2257000.0 !J kg-1
    REAL, PARAMETER                  :: Cv = 1004.0 !J K-1 kg-1 

    REAL, PARAMETER                  :: kg2ug    = 1.0e9

    LOGICAL, SAVE                    :: FIRST = .TRUE.
 
    !=======================================================================
    ! Routine starts here 
    !=======================================================================

    ! Identify this routine to MAPL
    Iam = 'GCC::CalcSpeciesDiagnostics_'

    ! Grid size
    IM = SIZE(Q,1)
    JM = SIZE(Q,2)
    LM = SIZE(Q,3)

    !=======================================================================
    ! Aerodynamic resistance
    !=======================================================================
!    ALLOCATE(Ra2m(IM,JM),Ra10m(IM,JM))
!    ALLOCATE(Lz0(IM,JM))
!    ! Compute aerodynamic resistance between 2m/10m and surface mid layer.
!    DO J=1,JM
!    DO I=1,IM
!       ! Ra from surface to surface mid-layer
!       ! Use t2m for all calculations. Only used for Monin-Obukhov.
!       t_         = T2M(I,J)
!       rho_       = AIRDENS(I,J,LM)
!       fhs_       = SH(I,J) 
!       ustar_     = USTAR(I,J)
!       dz_        = DELP(I,J,LM) / ( MAPL_GRAV * rho_ )
!       z0h_       = Z0H(I,J)
!       Ra_z0      = aerodynamic_resistance(t_, rho_, fhs_, ustar_, dz_, z0h_ )
!       Lz0(I,J)   = monin_obukhov_length(t_, rho_, fhs_, ustar_ ) 
!       ! Ra from surface to 2m 
!       !t_         = T2M(I,J)
!       dz_        = 2. * 2.
!       Ra_2m      = aerodynamic_resistance(t_, rho_, fhs_, ustar_, dz_, z0h_ )
!       Ra2m(I,J)  = 0.01 * ( Ra_z0 - Ra_2m )  ! convert to s cm-1
!       ! Ra from surface to 10m 
!       !t_         = T10M(I,J)
!       dz_        = 2. * 10.
!       Ra_10m     = aerodynamic_resistance(t_, rho_, fhs_, ustar_, dz_, z0h_ )
!       Ra10m(I,J) = 0.01 * ( Ra_z0 - Ra_10m ) ! convert to s cm-1
!    ENDDO
!    ENDDO

    CALL MAPL_GetPointer( EXPORT, Ptr2d,  'DryDepRa2m', &
                          NotFoundOk=.TRUE., __RC__ )
    IF ( ASSOCIATED(Ptr2d) ) Ptr2d = State_Chm%DryDepRa2m
    CALL MAPL_GetPointer( EXPORT, Ptr2d, 'DryDepRa10m', &
                          NotFoundOk=.TRUE., __RC__ )
    IF ( ASSOCIATED(Ptr2d) ) Ptr2d = State_Chm%DryDepRa10m
!    CALL MAPL_GetPointer( EXPORT, Ptr2d, 'MoninObukhov', &
!                          NotFoundOk=.TRUE., __RC__ )
!    IF ( ASSOCIATED(Ptr2d) ) Ptr2d = Lz0
!    DEALLOCATE(Lz0)

    !=======================================================================
    ! Exports in dry vol mixing ratio (v/v dry). Includes NOy. Convert from
    ! kg/kg total. 
    !=======================================================================
    CALL MAPL_GetPointer( EXPORT, NOy, 'NOy', NotFoundOk=.TRUE., __RC__ )
    IF ( ASSOCIATED(NOy) ) NOy = 0.0
    CALL MAPL_GetPointer( EXPORT, Bry, 'Bry', NotFoundOk=.TRUE., __RC__ )
    IF ( ASSOCIATED(Bry) ) Bry = 0.0
    CALL MAPL_GetPointer( EXPORT, Cly, 'Cly', NotFoundOk=.TRUE., __RC__ )
    IF ( ASSOCIATED(Cly) ) Cly = 0.0
    CALL MAPL_GetPointer( EXPORT, OrgCl, 'OrganicCl', NotFoundOk=.TRUE., __RC__ )
    IF ( ASSOCIATED(OrgCl) ) OrgCl = 0.0

    ! Check for PM25 diagnostics
    CALL MAPL_GetPointer( EXPORT, PtrPM25         , 'myPM25'         , &
                          NotFoundOk=.TRUE., __RC__ )
    CALL MAPL_GetPointer( EXPORT, PtrPM25_2m      , 'myPM25_2m'      , &
                          NotFoundOk=.TRUE., __RC__ )
    CALL MAPL_GetPointer( EXPORT, PtrPM25_10m     , 'myPM25_10m'     , & 
                          NotFoundOk=.TRUE., __RC__ )
    CALL MAPL_GetPointer( EXPORT, PtrPM25_SOA     , 'myPM25_SOA'     , &
                          NotFoundOk=.TRUE., __RC__ )
    CALL MAPL_GetPointer( EXPORT, PtrPM25_nitrates, 'myPM25_nitrates', &
                          NotFoundOk=.TRUE., __RC__ )
    CALL MAPL_GetPointer( EXPORT, PtrPM25_sulfates, 'myPM25_sulfates', &
                          NotFoundOk=.TRUE., __RC__ )
    CALL MAPL_GetPointer( EXPORT, PtrPM25_OC      , 'myPM25_OC'      , &
                          NotFoundOk=.TRUE., __RC__ )
    CALL MAPL_GetPointer( EXPORT, PtrPM25_BC      , 'myPM25_BC'      , &
                          NotFoundOk=.TRUE., __RC__ )
    CALL MAPL_GetPointer( EXPORT, PtrPM25_dust    , 'myPM25_dust'    , &
                          NotFoundOk=.TRUE., __RC__ )
    CALL MAPL_GetPointer( EXPORT, PtrPM25_seasalt , 'myPM25_seasalt' , &
                          NotFoundOk=.TRUE., __RC__ )
    DoPM25 = ( ASSOCIATED(PtrPM25)          .OR. ASSOCIATED(PtrPM25_2m)   .OR. &
               ASSOCIATED(PtrPM25_10m)      .OR. ASSOCIATED(PtrPM25_SOA)  .OR. &
               ASSOCIATED(PtrPM25_nitrates) .OR. ASSOCIATED(PtrPM25_sulfates) .OR. &
               ASSOCIATED(PtrPM25_OC)       .OR. ASSOCIATED(PtrPM25_BC)   .OR. &
               ASSOCIATED(PtrPM25_dust)     .OR. ASSOCIATED(PtrPM25_seasalt)        ) 
    IF ( ASSOCIATED(PtrPM25         ) ) PtrPM25          = 0.0
    IF ( ASSOCIATED(PtrPM25_2m      ) ) PtrPM25_2m       = 0.0
    IF ( ASSOCIATED(PtrPM25_10m     ) ) PtrPM25_10m      = 0.0
    IF ( ASSOCIATED(PtrPM25_nitrates) ) PtrPM25_nitrates = 0.0
    IF ( ASSOCIATED(PtrPM25_sulfates) ) PtrPM25_sulfates = 0.0
    IF ( ASSOCIATED(PtrPM25_SOA     ) ) PtrPM25_SOA      = 0.0
    IF ( ASSOCIATED(PtrPM25_OC      ) ) PtrPM25_OC       = 0.0
    IF ( ASSOCIATED(PtrPM25_BC      ) ) PtrPM25_BC       = 0.0
    IF ( ASSOCIATED(PtrPM25_dust    ) ) PtrPM25_dust     = 0.0
    IF ( ASSOCIATED(PtrPM25_seasalt ) ) PtrPM25_seasalt  = 0.0

    DO N=1,State_Chm%nSpecies
       SpcInfo   => State_Chm%SpcData(N)%Info ! Species database
       SpcName   =  TRIM(SpcInfo%Name)
!       ! Ignore family species
!       IF ( LEN(SpcName) > 1 ) THEN
!          IF ( SpcName(1:2) == 'T0' .OR.  &
!               SpcName(1:2) == 'T1' .OR.  &
!               SpcName(1:2) == 'T2'        ) CYCLE
!       ENDIF
!       IF ( LEN(SpcName) > 2 ) THEN
!          IF ( SpcName(1:3) == 'PT0' .OR.  &
!               SpcName(1:3) == 'PT1' .OR.  &
!               SpcName(1:3) == 'PT2'        ) CYCLE
!       ENDIF
       FieldName = TRIM(GPFX)//TRIM(SpcName)//'dry' ! FieldName

       ! See if field exists in export
       CALL MAPL_GetPointer( EXPORT, Ptr3D, FieldName, &
                             NotFoundOk=.TRUE., __RC__ )

       ! Need to fill at least one export?
       RunMe = .FALSE.
       IsNi  = .FALSE. 
       IsSu  = .FALSE.
       IsBC  = .FALSE.
       IsOC  = .FALSE.
       IsDu  = .FALSE.
       IsSS  = .FALSE.
       IsSOA = .FALSE.

       ! Is this a NOy species?
       IF ( ASSOCIATED(NOy) ) THEN
          SELECT CASE ( TRIM(SpcName) )
             CASE ( 'BrNO3', 'ClNO3', 'DHDN', 'ETHLN', 'HNO2', &
                    'HNO3',  'HNO4',  'HONIT'  )
                IsNOy = .TRUE.
             CASE ( 'IONITA', 'IPMN', 'ISN1', 'ISNIOA', 'ISNIOG' )
                IsNOy = .TRUE.
             CASE ( 'ISOPNB', 'ISOPND', 'MACRN', 'MPN', 'MVKN', &
                    'N2O5',   'NIT',    'NO',    'NO2', 'NO3' )
                IsNOy = .TRUE.
             CASE ( 'NPMN', 'ONIT', 'PAN', 'PROPNN', 'R4N2' )
                IsNOy = .TRUE.
             CASE DEFAULT
                IsNOy = .FALSE.
          END SELECT
       ELSE
          IsNOy = .FALSE.
       ENDIF
       IF ( IsNOy ) RunMe = .TRUE.

       ! Is this a Bry species?
       BrCoeff = 0.0
       IF ( ASSOCIATED(Bry) ) THEN
          SELECT CASE ( TRIM(SpcName) )
             CASE ( 'Br', 'BrO', 'HOBr', 'HBr', 'BrNO2', 'BrNO3', 'BrCl', 'IBr' )
                BrCoeff = 1.0
                IsBry   = .TRUE.
             CASE ( 'Br2' )
                BrCoeff = 2.0
                IsBry   = .TRUE.
             CASE DEFAULT
                IsBry = .FALSE.
          END SELECT
       ELSE
          IsBry = .FALSE.
       ENDIF
       IF ( IsBry ) RunMe = .TRUE.

       ! Is this a Cly species?
       ClCoeff = 0.0
       IF ( ASSOCIATED(Cly) ) THEN
          SELECT CASE ( TRIM(SpcName) )
             CASE ( 'Cl', 'ClO', 'OClO', 'ClOO', 'HOCl', 'HCl', 'ClNO2', 'ClNO3', 'BrCl', 'ICl' )
                ClCoeff = 1.0
                IsCly   = .TRUE.
             CASE ( 'Cl2', 'Cl2O2' )
                ClCoeff = 2.0
                IsCly   = .TRUE.
             CASE DEFAULT
                IsCly = .FALSE.
          END SELECT
       ELSE
          IsCly = .FALSE.
       ENDIF
       IF ( IsCly ) RunMe = .TRUE.

       ! Is this an OrgCl species?
       OrgClCoeff = 0.0
       IF ( ASSOCIATED(Cly) ) THEN
          SELECT CASE ( TRIM(SpcName) )
             CASE ( 'H1211', 'CFC115', 'CH3Cl', 'HCFC142b', 'HCFC22', 'CH2ICl' )
                OrgClCoeff = 1.0
                IsOrgCl    = .TRUE.
             CASE ( 'CFC114', 'CFC12', 'HCFC141b', 'HCFC123', 'CH2Cl2' )
                OrgClCoeff = 2.0
                IsOrgCl    = .TRUE.
             CASE ( 'CFC11', 'CFC113', 'CH3CCl3', 'CHCl3' ) 
                OrgClCoeff = 3.0
                IsOrgCl    = .TRUE.
             CASE ( 'CCl4' ) 
                OrgClCoeff = 4.0
                IsOrgCl    = .TRUE.
             CASE DEFAULT
                IsOrgCl = .FALSE.
          END SELECT
       ELSE
          IsOrgCl = .FALSE.
       ENDIF
       IF ( IsOrgCl ) RunMe = .TRUE.

       ! Is this a PM25 component?
       IF ( DoPM25 ) THEN
          SELECT CASE ( TRIM(SpcName) )
             CASE ( 'NH4', 'NIT' )
                IsNi  = .TRUE.
                RunMe = .TRUE.
             CASE ( 'SO4' )
                IsSu  = .TRUE.
                RunMe = .TRUE.
             CASE ( 'BCPI', 'BCPO' ) 
                IsBC  = .TRUE.
                RunMe = .TRUE.
             CASE ( 'OCPI', 'OCPO' ) 
                IsOC  = .TRUE.
                RunMe = .TRUE.
             CASE ( 'DST1', 'DST2' ) 
                IsDu  = .TRUE.
                RunMe = .TRUE.
             CASE ( 'SALA' )
                IsSS  = .TRUE.
                RunMe = .TRUE.
             CASE ( 'TSOA0', 'TSOA1', 'TSOA2', 'TSOA3', 'ISOA1', &
                    'ISOA2', 'ISOA3', 'ASOAN', 'ASOA1', 'ASOA2', 'ASOA3' )
                IsSOA = .TRUE.
                RunMe = .TRUE.
          END SELECT
       ENDIF
       IsPM25 = ( IsSOA .OR. IsNi .OR. IsSu .OR. IsOC .OR. IsBC .OR. &
                  IsDu .OR. IsSS )

       ! Check if we want concentration at 2M/10m
       CALL MAPL_GetPointer( EXPORT, Ptr2m , TRIM(SpcName)//'vv_2m' , &
                             NotFoundOk=.TRUE., __RC__ )
       CALL MAPL_GetPointer( EXPORT, Ptr10m, TRIM(SpcName)//'vv_10m', &
                             NotFoundOk=.TRUE., __RC__ )

       ! Fill exports
       IF ( ASSOCIATED(Ptr3D)  .OR. IsNOy .OR. ASSOCIATED(Ptr2m) .OR. &
            ASSOCIATED(Ptr10m) ) RunMe = .TRUE.
       IF ( RunMe ) THEN

          MW = SpcInfo%EmMW_g
          IF ( MW > 0.0 ) THEN
             FieldName = 'TRC_'//TRIM(SpcName)
          ELSE
             ! Get species and set MW to 1.0. This is ok because the internal
             ! state uses a MW of 1.0 for all species
             FieldName = 'SPC_'//TRIM(SpcName)
             MW = 1.0
             ! Cannot add to NOy if MW is unknown because it would screw up 
             ! unit conversion
             IF ( IsNOy ) THEN
                IsNOy = .FALSE.
                IF ( am_I_Root .AND. FIRST ) THEN
                   write(*,*) 'WARNING: Ignore species for NOy computation' //&
                              '  because MW is unknown: ', TRIM(SpcName)
                ENDIF
             ENDIF
             !IF ( ASSOCIATED(Ptr3D) .AND. am_I_Root .AND. FIRST ) THEN
             !   write(*,*)    &
             !      'WARNING: Attempt conversion of kg/kg total'// &
             !      ' to v/v dry but MW is unknown: ', TRIM(SpcName)
             !ENDIF
          ENDIF
          CALL MAPL_GetPointer( INTSTATE, PtrTmp, FieldName, NotFoundOK=.TRUE., RC=STATUS )
          ! On first fail try field with other prefix
          !IF ( STATUS /= ESMF_SUCCESS ) THEN
          IF ( .NOT. ASSOCIATED(PtrTmp) ) THEN
             IF ( FieldName(1:4)=='TRC_' ) THEN
                FieldName = 'SPC_'//TRIM(SpcName)
             ELSE
                FieldName = 'TRC_'//TRIM(SpcName)
             ENDIF
             CALL MAPL_GetPointer( INTSTATE, PtrTmp, FieldName, RC=STATUS )
             IF ( STATUS /= ESMF_SUCCESS ) THEN
                WRITE(*,*) 'Error reading ',TRIM(SpcName)
                VERIFY_(STATUS)
             ENDIF
          ENDIF

          !====================================================================
          ! Export in v/v
          !====================================================================
          IF ( ASSOCIATED(Ptr3D) ) Ptr3D = PtrTmp * ( MAPL_AIRMW / MW ) / &
                                           ( 1.0 - Q )

          !====================================================================
          ! NOy concentration
          !====================================================================
          IF ( IsNOy ) NOy = NOy + PtrTmp * ( MAPL_AIRMW / MW ) / ( 1.0 - Q )

          !====================================================================
          ! Bry concentration
          IF ( IsBry ) Bry = Bry + BrCoeff * PtrTmp * ( MAPL_AIRMW / MW ) / ( 1.0 - Q )
          !====================================================================

          !====================================================================
          ! Cly concentration
          !====================================================================
          IF ( IsCly ) Cly = Cly + ClCoeff * PtrTmp * ( MAPL_AIRMW / MW ) / ( 1.0 - Q )

          !====================================================================
          ! OrgCl concentration
          !====================================================================
          IF ( IsOrgCl ) OrgCl = OrgCl + OrgClCoeff * PtrTmp * ( MAPL_AIRMW / MW ) / ( 1.0 - Q )

          !==================================================================
          ! PM2.5 diagnostics. Calculate PM25 according to 
          ! http://wiki.seas.harvard.edu/geos-chem/index.php/Particulate_matter_in_GEOS-Chem
          !===================================================================
          IF ( IsPM25 ) THEN 
             ALLOCATE(MyPM25(IM,JM,LM))
             IF ( IsNi  ) THEN
                MyPM25 = 1.33 * PtrTmp * AIRDENS / ( 1.0-Q ) * kg2ug
                IF ( ASSOCIATED(PtrPM25_nitrates) ) &
                   PtrPM25_nitrates = PtrPM25_nitrates + MyPM25
             ENDIF
             IF ( IsSu  ) THEN
                MyPM25 = 1.33 * PtrTmp * AIRDENS / ( 1.0-Q ) * kg2ug
                IF ( ASSOCIATED(PtrPM25_sulfates) ) &
                   PtrPM25_sulfates = PtrPM25_sulfates + MyPM25
             ENDIF
             IF ( IsBC  ) THEN
                MyPM25 = PtrTmp * AIRDENS / ( 1.0-Q ) * kg2ug
                IF ( ASSOCIATED(PtrPM25_BC) ) &
                   PtrPM25_BC = PtrPM25_BC + MyPM25
             ENDIF
             IF ( IsSOA ) THEN
                MyPM25 = 1.16 * PtrTmp * AIRDENS / ( 1.0-Q ) * kg2ug
                IF ( ASSOCIATED(PtrPM25_SOA) ) &
                   PtrPM25_SOA = PtrPM25_SOA + MyPM25
             ENDIF
             IF ( IsSS  ) THEN
                MyPM25 = 1.86 * PtrTmp * AIRDENS / ( 1.0-Q ) * kg2ug
                IF ( ASSOCIATED(PtrPM25_seasalt) ) &
                   PtrPM25_seasalt = PtrPM25_seasalt + MyPM25
             ENDIF                                 
             IF ( IsOC  ) THEN
                MyPM25 = 2.1  * PtrTmp * AIRDENS / ( 1.0-Q ) * kg2ug
                IF ( TRIM(SpcName) == 'OCPI' ) MyPM25 = 1.16 * MyPM25 
                IF ( ASSOCIATED(PtrPM25_OC) ) &
                   PtrPM25_OC = PtrPM25_OC + MyPM25
             ENDIF
             IF ( IsDu ) THEN 
                MyPM25 = PtrTmp * AIRDENS / ( 1.0-Q ) * kg2ug
                IF ( TRIM(SpcName) == 'DST2' ) MyPM25 = 0.38 * MyPM25 
                IF ( ASSOCIATED(PtrPM25_dust) ) &
                   PtrPM25_dust = PtrPM25_dust + MyPM25
             ENDIF
 
             ! Total PM25
             IF ( ASSOCIATED(PtrPM25) ) PtrPM25 = PtrPM25 + MyPM25
          ENDIF

          !====================================================================
          ! 2m and 10m concentration, corrected by aerodynamic resistance 
          ! (Zhang et al. 2012)
          !====================================================================
          IF ( ASSOCIATED(Ptr2m) .OR. &
             (IsPM25.AND.ASSOCIATED(PtrPM25_2m)) ) THEN
             ALLOCATE(CONV(SIZE(Ptr2m,1),SIZE(Ptr2m,2)))
             DryID = SpcInfo%DryDepID
             IF ( DryID > 0 ) THEN
                IF ( .NOT. ASSOCIATED(State_Diag%DryDepVel ) ) THEN 
                   WRITE(*,*)     &
                       'Cannot compute 2M concentration - need diagnostics '// &
                       'of drydep velocity - please activate - ' //TRIM(SpcName)
                   ASSERT_(.FALSE.)
                ENDIF      
                CONV(:,:) = 1.0 - ( State_Chm%DryDepRa2m(:,:) *  &
                            State_Diag%DryDepVel(:,:,DryID) )
             ELSE
                CONV(:,:) = 1.0
             ENDIF
             ! Don't allow negative adjustment factor. 
             ! Arbitrarily restrict correction to 0.25 
             WHERE ( CONV < 0.0 )
                CONV = 0.25
             ENDWHERE
             IF ( ASSOCIATED(Ptr2m) ) Ptr2m = CONV * PtrTmp(:,:,LM) * &
                     ( MAPL_AIRMW / MW ) / ( 1.0 - Q(:,:,LM) )
             IF ( IsPM25 .AND. ASSOCIATED(PtrPM25_2m) )  &
                     PtrPM25_2m = PtrPM25_2m + CONV * MyPM25(:,:,LM)
             DEALLOCATE(CONV) 
          ENDIF
          ! 10m concentration, corrected by aerodynamic resistance 
          ! (Zhang et al. 2012)
          IF ( ASSOCIATED(Ptr10m) .OR.  &
           (IsPM25.AND.ASSOCIATED(PtrPM25_10m)) )  THEN
             LM = SIZE(PtrTmp,3)
             ALLOCATE(CONV(SIZE(Ptr10m,1),SIZE(Ptr10m,2)))
             DryID = SpcInfo%DryDepID
             IF ( DryID > 0 ) THEN
                IF ( .NOT. ASSOCIATED(State_Diag%DryDepVel ) ) THEN
                   WRITE(*,*)      &
                     'Cannot compute 10M concentration - need diagnostics '// &
                     'of drydep velocity - please activate - ' //TRIM(SpcName)
                   ASSERT_(.FALSE.)
                ENDIF      
                CONV(:,:) = 1.0 - ( State_Chm%DryDepRa10m(:,:) * &
                            State_Diag%DryDepVel(:,:,DryID) )
             ELSE
                CONV(:,:) = 1.0
             ENDIF  
             ! Don't allow negative adjustment factor. Arbitrarily restrict 
             ! correction to 0.25 
             WHERE ( CONV < 0.0 )
                CONV = 0.25
             ENDWHERE
             IF ( ASSOCIATED(Ptr10m) ) Ptr10m = CONV * PtrTmp(:,:,LM) * &
                                 ( MAPL_AIRMW / MW ) / ( 1.0 - Q(:,:,LM) )
             IF ( IsPM25 .AND. ASSOCIATED(PtrPM25_10m) ) PtrPM25_10m =  &
                                 PtrPM25_10m + CONV * MyPM25(:,:,LM)
             DEALLOCATE(CONV) 
          ENDIF
          ! Cleanup
          IF (ALLOCATED(MyPM25)) DEALLOCATE(MyPM25)
       ENDIF
    ENDDO

    !=======================================================================
    ! All done 
    !=======================================================================

    ! Cleanup
    !IF ( ALLOCATED(Ra2m ) ) DEALLOCATE(Ra2m )
    !IF ( ALLOCATED(Ra10m) ) DEALLOCATE(Ra10m)

    ! Successful return
    FIRST = .FALSE.
    RETURN_(ESMF_SUCCESS)

  END SUBROUTINE CalcSpeciesDiagnostics_
!EOC
!------------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 910.1 and      !
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: MassWeightedOH_ 
!
! !DESCRIPTION: MassWeightedOH_ computes vertically integrated OH weighted 
!  by mass (troposphere only). 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE MassWeightedOH_ ( am_I_Root, Input_Opt, State_Met, State_Chm, &
                               State_Diag, IMPORT, EXPORT, INTSTATE, Q, PLE, TROPP, RC )
!
! !USES:
!
!
! !INPUT/OUTPUT PARAMETERS:
!
    LOGICAL,             INTENT(IN)            :: am_I_Root
    TYPE(OptInput),      INTENT(INOUT)         :: Input_Opt
    TYPE(MetState),      INTENT(INOUT)         :: State_Met
    TYPE(ChmState),      INTENT(INOUT)         :: State_Chm
    TYPE(DgnState),      INTENT(INOUT)         :: State_Diag
    TYPE(ESMF_State),    INTENT(INOUT)         :: Import   ! Import State
    TYPE(ESMF_State),    INTENT(INOUT)         :: Export   ! Export State
    TYPE(ESMF_STATE),    INTENT(INOUT)         :: INTSTATE
    REAL,                POINTER               :: Q(:,:,:)
    REAL,                POINTER               :: PLE  (:,:,:)
    REAL,                POINTER               :: TROPP(:,:  )
!
! !OUTPUT PARAMETERS:
!
    INTEGER,             INTENT(INOUT)         :: RC       ! Success or failure?
!
! !REMARKS:
!
! !REVISION HISTORY:
!  01 Feb 2019 - C. Keller   - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! LOCAL VARIABLES:
!
    ! Objects

    ! Scalars
    INTEGER                    :: STATUS
    INTEGER                    :: I, J, L, N, IM, JM, LM, LB, id_OH
    CHARACTER(LEN=ESMF_MAXSTR) :: Iam           ! Gridded component name
    CHARACTER(LEN=ESMF_MAXSTR) :: SpcName
    REAL                       :: MW_kg
    REAL, POINTER              :: MeanOH(:,:)
    REAL*8                     :: air_mass, oh_mass, oh_molec, wgt, mwgt, molair
    TYPE(Species), POINTER     :: SpcInfo

    !=======================================================================
    ! Routine starts here 
    !=======================================================================

    ! Identify this routine to MAPL
    Iam = 'GCC::MassWeightedOH_'

    ! Check for diagnostics
    CALL MAPL_GetPointer( EXPORT, MeanOH,  'OH_AirMassWgt', NotFoundOk=.TRUE., __RC__ )
    IF ( ASSOCIATED(MeanOH) ) THEN
       ! Reset
       MeanOH = 0.0

       ! Get OH species index 
       id_OH = -1
       DO N=1,State_Chm%nSpecies
          SpcInfo   => State_Chm%SpcData(N)%Info ! Species database
          IF ( TRIM(SpcInfo%Name) == "OH" ) THEN
             id_OH = N
             IF ( SpcInfo%emMW_g > 0.0 ) THEN
                MW_kg = SpcInfo%emMW_g * 1.e-3
             ELSE
                MW_kg = 1.0 * 1.e-3
             ENDIF
             EXIT
          ENDIF
       ENDDO
       ASSERT_( id_OH > 0 )

       ! Grid size
       IM = SIZE(Q,1)
       JM = SIZE(Q,2)
       LM = SIZE(Q,3)

       ! Lower bound of PLE 3rd dim
       LB = LBOUND(PLE,3)

       ! Compute mass weighted OH per grid box, troposphere only
       DO J = 1,JM
       DO I = 1,IM
          air_mass = 0.0d0
          oh_mass  = 0.0d0
          DO L = 1,LM
             wgt  = MAX(0.0,MIN(1.0,(PLE(I,J,L+LB)-TROPP(I,J)) &
                  /(PLE(I,J,L+LB)-PLE(I,J,L+LB-1))))
             IF ( wgt > 0.0d0 ) THEN
                molair   = State_Met%AIRNUMDEN(I,J,LM-L+1) / AVO 
                mwgt     = wgt * molair * State_Met%AIRVOL(I,J,LM-L+1) * 1d-9
                air_mass = air_mass + mwgt 
                oh_molec = State_Chm%Species(I,J,LM-L+1,id_OH) / ( 1 - Q(I,J,L) ) &
                         * State_Met%AIRDEN(I,J,LM-L+1) * ( AVO / MW_kg ) / 1d+9
                oh_mass  = oh_mass +  oh_molec * 1d-5 * mwgt
             ENDIF
          ENDDO
          IF ( air_mass > 0.0d0 ) MeanOH(I,J) = oh_mass / air_mass * 1e5 *1e3
       ENDDO
       ENDDO
    ENDIF
 
    RETURN_(ESMF_SUCCESS)

  END SUBROUTINE MassWeightedOH_
!------------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 910.1 and      !
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: InitFromFile_
!
! !DESCRIPTION: InitFromFile_ initializes the GEOS-Chem species values from
!  external data. 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE InitFromFile_( GC, Import, Internal, Export, Clock, Input_Opt, &
                            State_Met, State_Chm, Q, PLE, TROPP, First, RC )
!
! !USES:
!
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ESMF_GridComp), INTENT(INOUT)         :: GC       ! Ref. to this GridComp
    TYPE(ESMF_State),    INTENT(INOUT)         :: Import   ! Import State
    TYPE(ESMF_STATE),    INTENT(INOUT)         :: Internal ! Internal state
    TYPE(ESMF_State),    INTENT(INOUT)         :: Export   ! Export State
    TYPE(ESMF_Clock),    INTENT(INOUT)         :: Clock    ! ESMF Clock object
    TYPE(OptInput)                             :: Input_Opt 
    TYPE(MetState)                             :: State_Met 
    TYPE(ChmState)                             :: State_Chm 
    TYPE(ESMF_Time)                            :: currTime
    REAL,                INTENT(IN)            :: Q(:,:,:)
    REAL,                POINTER               :: PLE(:,:,:)
    REAL,                POINTER               :: TROPP(:,:)
    LOGICAL,             INTENT(IN)            :: First
!
! !OUTPUT PARAMETERS:
!
    INTEGER,             INTENT(OUT)           :: RC       ! Success or failure?
!
! !REMARKS:
!
! !REVISION HISTORY:
!  18 Mar 2017 - C. Keller   - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! LOCAL VARIABLES:
!
    ! Objects
 
    ! Scalars
    LOGICAL                    :: am_I_Root     ! Are we on the root PET?
    TYPE(ESMF_Config)          :: GeosCF      ! ESMF Config obj (GEOSCHEM*.rc) 
    CHARACTER(LEN=ESMF_MAXSTR) :: Iam, compName ! Gridded component name
    CHARACTER(LEN=ESMF_MAXSTR) :: FldName
    CHARACTER(LEN=ESMF_MAXSTR) :: SpcName
    CHARACTER(LEN=255)         :: ifile
    TYPE(MAPL_SimpleBundle)    :: VarBundle
    TYPE(ESMF_Grid)            :: grid
    TYPE(ESMF_TIME)            :: time
    TYPE(ESMF_Field)           :: iFld 
    REAL, POINTER              :: Ptr3D(:,:,:) => NULL()
    REAL, ALLOCATABLE          :: Scal(:,:), Temp(:,:), wgt1(:,:), wgt2(:,:)
    INTEGER                    :: varID, fid, N, L, LM2
    INTEGER                    :: L1, L2, INC
    INTEGER                    :: IM, JM, LM, LB
    INTEGER                    :: STATUS
    INTEGER                    :: nymd, nhms, yy, mm, dd, h, m, s, incSecs
    REAL                       :: MW
    LOGICAL                    :: FileExists
    INTEGER                    :: DoIt, idx, x1, x2
    LOGICAL                    :: ReadGMI 
    LOGICAL                    :: OnGeosLev 
    LOGICAL                    :: AboveTroppOnly
    LOGICAL                    :: IsInPPBV
    LOGICAL                    :: DoUpdate
    INTEGER                    :: TopLev
    CHARACTER(LEN=ESMF_MAXSTR) :: GmiTmpl 

    ! Read GMI file
    CHARACTER(LEN=ESMF_MAXSTR) :: GmiFldName
    CHARACTER(LEN=255)         :: Gmiifile
    TYPE(MAPL_SimpleBundle)    :: GmiVarBundle
    TYPE(ESMF_TIME)            :: Gmitime
    LOGICAL                    :: GmiFileExists
    INTEGER, SAVE              :: OnlyOnFirst = -999 

    !=======================================================================
    ! Initialization
    !=======================================================================

    ! Are we on the root PET
    am_I_Root = MAPL_Am_I_Root()

    ! Set up traceback info
    CALL ESMF_GridCompGet( GC, name=compName, grid=grid, __RC__ )

    ! Identify this routine to MAPL
    Iam = TRIM(compName)//'::InitFromFile_'

    ! Get GEOS-Chem resource file
    CALL Extract_( GC, Clock, GeosCF=GeosCF, __RC__ )

    ! Check if we need to do update
    IF ( OnlyOnFirst < 0 ) THEN 
        CALL ESMF_ConfigGetAttribute( GeosCF, DoIt, Label = 'ONLY_ON_FIRST_STEP:', Default=1, __RC__ )
        IF ( DoIt == 1 ) THEN
           OnlyOnFirst = 1
        ELSE
           OnlyOnFirst = 0
        ENDIF
    ENDIF
    DoUpdate = .FALSE.
    IF ( OnlyOnFirst == 0             ) DoUpdate = .TRUE.
    IF ( OnlyOnFirst == 1 .AND. First ) DoUpdate = .TRUE.

    ! Do the following only if we need to...
    IF ( DoUpdate ) THEN

    ! Array size 
    IM = SIZE(Q,1)
    JM = SIZE(Q,2)
    LM = SIZE(Q,3)

    ! Lower bound of PLE 3rd dim
    LB = LBOUND(PLE,3)

    ! Name of file to read internal state fields from
    CALL ESMF_ConfigGetAttribute( GeosCF, ifile, Label = "INIT_SPC_FILE:", __RC__ ) 
    IF ( am_I_Root ) WRITE(*,*) TRIM(Iam)//': reading species from '//TRIM(ifile)

    ! Check if file exists
    INQUIRE( FILE=TRIM(ifile), EXIST=FileExists )
    IF ( .NOT. FileExists ) THEN
       IF ( am_I_Root ) WRITE(*,*) 'File does not exist: ',TRIM(ifile)
       ASSERT_(.FALSE.)
    ENDIF

    ! Check for other flags
    CALL ESMF_ConfigGetAttribute( GeosCF, DoIt, Label = 'DATA_ON_GEOS_LEVELS:', Default=0, __RC__ )
    OnGeosLev = ( DoIt == 1 )
    CALL ESMF_ConfigGetAttribute( GeosCF, DoIt, Label = 'ONLY_ABOVE_TROPOPAUSE:', Default=0, __RC__ )
    AboveTroppOnly = ( DoIt == 1 )
    CALL ESMF_ConfigGetAttribute( GeosCF, DoIt, Label = 'DATA_IS_IN_PPBV:', Default=1, __RC__ )
    IsInPPBV = ( DoIt == 1 )
    CALL ESMF_ConfigGetAttribute( GeosCF, TopLev, Label = 'DO_NOT_OVERWRITE_ABOVE_LEVEL:', Default=LM, __RC__ )
    IF ( TopLev < 0 ) TopLev = LM

    ! Verbose
    IF ( am_I_Root ) THEN
       WRITE(*,*) 'Will use the following settings to overwrite restart variables:'
       WRITE(*,*) 'External data is in ppbv: ',IsInPPBV
       WRITE(*,*) 'External data is on GEOS levels: ',OnGeosLev
       WRITE(*,*) 'Only overwrite above tropopause: ',AboveTroppOnly
       WRITE(*,*) 'Maximum valid level (will be used above that level): ',TopLev
    ENDIF

    ! Check for GMI flags
    CALL ESMF_ConfigGetAttribute( GeosCF, DoIt, Label = "USE_GMI_MESO:", Default = 0, __RC__ )
    ReadGMI = ( DoIt == 1 )
    IF ( ReadGMI ) THEN
       CALL ESMF_ConfigGetAttribute( GeosCF, GmiTmpl, Label = "GMI_TEMPLATE:", __RC__ )
    ENDIF

    ! Get time stamp on file
    call GFIO_Open( TRIM(ifile), 1, fid, STATUS )
    ASSERT_(STATUS==0)
    call GetBegDateTime ( fid, nymd, nhms, incSecs, STATUS )
    ASSERT_(STATUS==0)
    caLL GFIO_Close( fid, STATUS )
    ASSERT_(STATUS==0)
    yy = nymd/10000
    mm = (nymd-yy*10000) / 100
    dd = nymd - (10000*yy + mm*100)
    h  = nhms/10000
    m  = (nhms- h*10000) / 100
    s  = nhms - (10000*h  +  m*100)
    call ESMF_TimeSet(time, yy=yy, mm=mm, dd=dd, h=h, m=m, s=s)

    ! Read file
    VarBundle = MAPL_SimpleBundleRead ( TRIM(iFile), 'GCCinit', grid, time, __RC__ )
          
    ! Scal is the array with scale factors 
    ALLOCATE(Scal(IM,JM),Temp(IM,JM),wgt1(IM,JM),wgt2(IM,JM))
    Scal(:,:) = 1.0    
    Temp(:,:) = 0.0
    wgt1(:,:) = 0.0
    wgt2(:,:) = 1.0

    ! Loop over all species
    DO N = 1, State_Chm%nSpecies

       ! Molecular weight. Note: -1.0 for non-advected species
       MW = State_Chm%SpcData(N)%Info%emMW_g
       IF ( MW < 0.0 ) MW = 1.0

       ! Construct field name
       SpcName = TRIM(State_Chm%SpcData(N)%Info%Name)

       ! Check if variable is in file
       FldName = 'SPC_'//TRIM(SpcName)
       VarID = MAPL_SimpleBundleGetIndex ( VarBundle, trim(FldName), 3, RC=STATUS, QUIET=.TRUE. )
   
       ! Check other fieldname if default one is not found
       IF ( VarID <= 0 ) THEN
          FldName = 'TRC_'//TRIM(SpcName)
          VarID = MAPL_SimpleBundleGetIndex ( VarBundle, trim(FldName), 3, RC=STATUS, QUIET=.TRUE. )
       ENDIF
       IF ( VarID <= 0 ) THEN
          FldName = TRIM(SpcName)
          VarID = MAPL_SimpleBundleGetIndex ( VarBundle, trim(FldName), 3, RC=STATUS, QUIET=.TRUE. )
       ENDIF
       IF ( VarID > 0 ) THEN
          ! Make sure vertical dimensions match
          LM2 = SIZE(VarBundle%r3(VarID)%q,3)

          ! Error if vertical dimensions do not agree
          IF ( LM2 /= LM ) THEN
             IF ( am_I_Root ) THEN
                WRITE(*,*) 'Wrong # of vert. levels for variable ',TRIM(FldName), ' ',LM2,' vs. ',LM
             ENDIF
             ASSERT_( LM==LM2 )
          ENDIF

          ! Loop over all vertical levels
          DO L = 1, LM
             ! Scale factor for unit conversion
             IF ( IsInPPBV ) THEN
                Scal(:,:) =  MW / MAPL_AIRMW * ( 1 - Q(:,:,L) )
                IF(L==1 .and. am_I_Root ) WRITE(*,*) 'Convert units from ppbv to kg/kg: ',TRIM(FldName), MW
             ENDIF

             ! Pass to temporary array
             IF ( OnGeosLev ) THEN
                Temp(:,:) = VarBundle%r3(VarID)%q(:,:,L) * Scal 
             ELSE
                Temp(:,:) = VarBundle%r3(VarID)%q(:,:,LM-L+1) * Scal 
             ENDIF

             ! Flag for stratosphere only
             IF ( AboveTroppOnly ) THEN
                wgt1 = MAX(0.0,MIN(1.0,(PLE(:,:,L+LB)-TROPP(:,:))/(PLE(:,:,L+LB)-PLE(:,:,L+LB-1))))
                wgt2 = 1.0 - wgt1 
             ENDIF

             ! Pass to State_Chm
             State_Chm%Species(:,:,LM-L+1,N) = State_Chm%Species(:,:,LM-L+1,N)*wgt1 + Temp(:,:)*wgt2
          ENDDO
 
          ! Check for cap at given level 
          IF ( TopLev < LM ) THEN
             DO L = TopLev+1,LM
                State_Chm%Species(:,:,L,N) = State_Chm%Species(:,:,TopLev,N)
             ENDDO
             IF ( am_I_Root ) WRITE(*,*) 'Extend values from level ',TopLev,' to top of atmosphere: ',TRIM(FldName)
          ENDIF

          ! Verbose
          IF ( am_I_Root ) WRITE(*,*) 'Species initialized from external field: ',TRIM(FldName)

       ELSE
          IF ( am_I_Root ) WRITE(*,*) 'Species unchanged, field not found for species ',TRIM(SpcName)
       ENDIF

       ! ---------------------------
       ! Try to read GMI data
       ! ---------------------------
       IF ( ReadGMI ) THEN
          ! Get file name
          Gmiifile = GmiTmpl
          idx = INDEX(Gmiifile,'%spc')
          IF ( idx > 0 ) THEN
             x1 = idx + 4
             x2 = LEN(TRIM(Gmiifile))
             Gmiifile = TRIM(Gmiifile(1:idx-1))//TRIM(SpcName)//TRIM(Gmiifile(x1:x2))
          ENDIF
          INQUIRE( FILE=TRIM(Gmiifile), EXIST=GmiFileExists )
   
          IF ( GmiFileExists ) THEN 
   
             ! Get time stamp on file
             call GFIO_Open( Gmiifile, 1, fid, STATUS )
             ASSERT_(STATUS==0)
             call GetBegDateTime ( fid, nymd, nhms, incSecs, STATUS )
             ASSERT_(STATUS==0)
             caLL GFIO_Close( fid, STATUS )
             ASSERT_(STATUS==0)
             yy = nymd/10000
             mm = (nymd-yy*10000) / 100
             dd = nymd - (10000*yy + mm*100)
             h  = nhms/10000
             m  = (nhms- h*10000) / 100
             s  = nhms - (10000*h  +  m*100)
             call ESMF_TimeSet(Gmitime, yy=yy, mm=7, dd=6, h=h, m=m, s=s)
  
             ! Read data 
             GmiVarBundle = MAPL_SimpleBundleRead ( TRIM(GmiiFile), 'GCCinitGMI', grid, Gmitime, __RC__ )
      
             ! Check if variable is in file
             VarID = MAPL_SimpleBundleGetIndex ( GmiVarBundle, 'species', 3, RC=STATUS, QUIET=.TRUE. )
             IF ( VarID > 0 ) THEN
                ! Pass to State_Chm, convert v/v to kg/kg.
                State_Chm%Species(:,:,60:72,N) = VarBundle%r3(VarID)%q(:,:,13:1:-1) * MW / MAPL_AIRMW * ( 1 - Q(:,:,13:1:-1) )
                IF ( am_I_Root ) WRITE(*,*) 'Use GMI concentrations in mesosphere: ',TRIM(SpcName)
             ENDIF
   
          ELSE 
             IF ( am_I_Root ) WRITE(*,*) 'No GMI file found: ',TRIM(Gmiifile)
          ENDIF
       ENDIF

    ENDDO

    ! Deallocate helper array
    IF ( ALLOCATED(Scal) ) DEALLOCATE(Scal) 
    IF ( ALLOCATED(Temp) ) DEALLOCATE(Temp) 
    IF ( ALLOCATED(wgt1) ) DEALLOCATE(wgt1)
    IF ( ALLOCATED(wgt2) ) DEALLOCATE(wgt2)

    ! All done
    CALL MAPL_SimpleBundleDestroy ( VarBundle, __RC__ )

    ENDIF ! DoUpdate 

    ! Return
    RETURN_(ESMF_SUCCESS)

  END SUBROUTINE InitFromFile_ 
!EOC
!------------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 910.1 and      !
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: SetAnaO3_ 
!
! !DESCRIPTION: 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE SetAnaO3_( GC, Import, Internal, Export, Clock, &
                        Input_Opt,  State_Met, State_Chm, Q, RC )
!
! !USES:
!
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ESMF_GridComp), INTENT(INOUT)         :: GC       ! Ref. to this GridComp
    TYPE(ESMF_State),    INTENT(INOUT)         :: Import   ! Import State
    TYPE(ESMF_STATE),    INTENT(INOUT)         :: Internal ! Internal state
    TYPE(ESMF_State),    INTENT(INOUT)         :: Export   ! Export State
    TYPE(ESMF_Clock),    INTENT(INOUT)         :: Clock    ! ESMF Clock object
    TYPE(OptInput)                             :: Input_Opt 
    TYPE(MetState)                             :: State_Met 
    TYPE(ChmState)                             :: State_Chm 
    REAL,                INTENT(IN)            :: Q(:,:,:)
!
! !OUTPUT PARAMETERS:
!
    INTEGER,             INTENT(OUT)           :: RC       ! Success or failure?
!
! !REMARKS:
!
! !REVISION HISTORY:
!  18 Mar 2017 - C. Keller   - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! LOCAL VARIABLES:
!
    ! Objects
 
    ! Scalars
    LOGICAL                    :: am_I_Root     ! Are we on the root PET?
    LOGICAL                    :: TimeForAna
    REAL, POINTER              :: ANAO3(:,:,:)
    REAL, POINTER              :: O3INC(:,:,:)
    CHARACTER(LEN=ESMF_MAXSTR) :: Iam, compName ! Gridded component name
    CHARACTER(LEN=ESMF_MAXSTR) :: SpcName
    INTEGER                    :: I, J, L, N, LR, IM, JM, LM 
    INTEGER                    :: STATUS
    REAL                       :: O3new, MW, ifrac
    INTEGER                    :: idx
    CHARACTER(LEN=2)           :: SL1, SL2 

    ! To read from file
    LOGICAL                    :: HasFile 
    CHARACTER(LEN=ESMF_MAXSTR) :: ifile
    CHARACTER(LEN=ESMF_MAXSTR) :: VarName
    TYPE(MAPL_SimpleBundle)    :: VarBundle
    TYPE(ESMF_Grid)            :: grid
    TYPE(ESMF_TIME)            :: currTime 
    TYPE(ESMF_TIME)            :: fileTime 
    TYPE(ESMF_Field)           :: iFld
    INTEGER                    :: VarID, fid 
    INTEGER                    :: nymd, nhms, yy, mm, dd, h, m, s, incSecs    
    CHARACTER(LEN=4)           :: syy
    CHARACTER(LEN=2)           :: smm, sdd, sh, sm
 
    !=======================================================================
    ! Initialization
    !=======================================================================

    ! Are we on the root PET
    am_I_Root = MAPL_Am_I_Root()

    ! Set up traceback info
    CALL ESMF_GridCompGet( GC, name=compName, grid=grid,  __RC__ )

    ! Identify this routine to MAPL
    Iam = TRIM(compName)//'::SetAnaO3_'

    ! Is it time to do the analysis? For now, do analysis every three hours
    ! (middle of 3-hour time step)
    ! Get current time
    CALL ESMF_ClockGet( Clock, currTime = currTime, __RC__ )
    CALL ESMF_TimeGet( currTime, yy=yy, mm=mm, dd=dd, h=h, m=m, s=s, __RC__ )
    IF ( m==30 .AND. MOD(h,3)==1 ) THEN
       TimeForAna = .TRUE.
    ELSE
       TimeForAna = .FALSE.
    ENDIF

    ! Diagnostics
    CALL MAPL_GetPointer ( Export, O3INC, 'ANA_O3_INC', NotFoundOk=.TRUE., __RC__ )
    IF ( ASSOCIATED(O3INC) ) O3INC = 0.0
   
    ! Get target ozone field 
    IF ( TimeForAna ) THEN

       ! Verbose
       WRITE(SL1,'(I2)') ANAO3L1
       WRITE(SL2,'(I2)') ANAO3L4
       IF ( am_I_Root ) WRITE(*,*) 'GEOS-Chem: nudging ozone between level '//SL1//' and '//SL2

       ! Get analysis O3 field from import in kg/kg
       IF ( LPCHEMO3 ) THEN
          CALL MAPL_GetPointer ( IMPORT, ANAO3, 'PCHEM_O3', __RC__ )
   
       ! If not PCHEM, try to read it from the specified file
       ELSE
          !CALL MAPL_GetPointer ( IMPORT, ANAO3, 'ANA_O3', __RC__ )
          ! Parse file name
          ifile = ANAO3FILE

          ! testing only
          IF ( am_I_Root ) WRITE(*,*) 'parsing '//TRIM(ifile)

          write(syy,'(I4.4)') yy
          CALL ReplaceChar ( ifile, '%y4', syy )
          write(smm,'(I2.2)') mm
          CALL ReplaceChar ( ifile, '%m2', smm )
          write(sdd,'(I2.2)') dd
          CALL ReplaceChar ( ifile, '%d2', sdd )
          write(sh,'(I2.2)') h
          CALL ReplaceChar ( ifile, '%h2', sh  )
          write(sm,'(I2.2)') m
          CALL ReplaceChar ( ifile, '%n2', sm  )
  
          ! testing only
          IF ( am_I_Root ) WRITE(*,*) 'parsed file: '//TRIM(ifile)

          ! Check if file exists
          INQUIRE( FILE=ifile, EXIST=HasFile ) 
          IF ( HasFile ) THEN
 
             ! Get time stamp on file
             call GFIO_Open( ifile, 1, fid, STATUS )
             ASSERT_(STATUS==0)
             call GetBegDateTime ( fid, nymd, nhms, incSecs, STATUS )
             ASSERT_(STATUS==0)
             caLL GFIO_Close( fid, STATUS )
             ASSERT_(STATUS==0)
             yy = nymd/10000
             mm = (nymd-yy*10000) / 100
             dd = nymd - (10000*yy + mm*100)
             h  = nhms/10000
             m  = (nhms- h*10000) / 100
             s  = nhms - (10000*h  +  m*100)
             call ESMF_TimeSet(fileTime, yy=yy, mm=mm, dd=dd, h=h, m=m, s=s)
   
             ! Read file
             VarBundle =  MAPL_SimpleBundleRead ( TRIM(ifile), 'GCCAnaO3', grid, fileTime, __RC__ )
             VarName   =  'O3'
             VarID     =  MAPL_SimpleBundleGetIndex ( VarBundle, trim(VarName), 3, RC=STATUS, QUIET=.TRUE. )
             ANAO3     => VarBundle%r3(VarID)%q
             IF ( am_I_Root ) WRITE(*,*) 'Use analysis ozone from '//TRIM(ifile)
          ELSE
             TimeForAna = .FALSE.
             IF ( am_I_Root ) WRITE(*,*) 'File not found - skip ozone nudging: '//TRIM(ifile)
          ENDIF   
       ENDIF
    ENDIF

    ! Apply ozone analysis increment if it's time to do so
    IF ( TimeForAna ) THEN  
 
       ! Select ozone index
       ! Loop over all species
       N = -1
       DO I = 1, State_Chm%nSpecies
          IF ( TRIM(State_Chm%SpcData(I)%Info%Name) == 'O3' ) THEN
             N = I 
          ENDIF
       ENDDO
       ASSERT_(N > 0)

       ! Molecular weight. Note: -1.0 for non-advected species
   !    MW = State_Chm%SpcData(I)%Info%emMW_g
   !    IF ( MW < 0.0 ) MW = 48.0 
   
       ! # of vertical levels of Q
       IM = SIZE(ANAO3,1)
       JM = SIZE(ANAO3,2)
       LM = SIZE(ANAO3,3)
   
       ! Loop over all relevant levels
       DO L = ANAO3L1, ANAO3L4
   
          ! LR is the reverse of L
          LR = LM - L + 1
   
          ! Get fraction of analysis field to be used. 
          ! This goes gradually from 0.0 to ANAO3FR between ANAO3L1 to ANAO3L2, 
          ! and is ANAO3FR above
          !IF ( ( L >= ANAO3L2 ) .OR. ( ANAO3L1 == ANAO3L2 ) ) THEN
          ifrac = 0.0
          IF ( L >= ANAO3L2 .AND. L <=  ANAO3L3 ) THEN
             ifrac = ANAO3FR 
          ELSEIF ( L < ANAO3L2 ) THEN
             ifrac = ANAO3FR * ( (L-ANAO3L1) / (ANAO3L2-ANAO3L1) )
          ELSEIF ( L > ANAO3L3 ) THEN
             ifrac = ANAO3FR * ( (ANAO3L4-L) / (ANAO3L4-ANAO3L3) )
          ENDIF
          ifrac = max(0.0,min(1.0,ifrac))
   
          ! Pass to State_Chm species array. PCHEM ozone should never be zero or smaller!
          DO J = 1,JM
          DO I = 1,IM
             IF ( ANAO3(I,J,LR) > 0.0 ) THEN
                O3new = ( (1.0-ifrac) * State_Chm%Species(I,J,L,N) ) &
                      + (      ifrac  * ANAO3(I,J,LR)              )
                IF ( ASSOCIATED(O3INC) ) O3INC(I,J,LR) = O3new - State_Chm%Species(I,J,L,N) 
                State_Chm%Species(I,J,L,N) = O3new 
             ENDIF
          ENDDO
          ENDDO
       ENDDO
   
       ! Clean up
       IF ( ASSOCIATED(ANAO3 ) ) ANAO3 => NULL()
       IF ( .NOT. LPCHEMO3 ) CALL MAPL_SimpleBundleDestroy ( VarBundle, __RC__ )
    ENDIF

    ! All done
    RETURN_(ESMF_SUCCESS)

  END SUBROUTINE SetAnaO3_ 
!EOC

!------------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 910.1 and      !
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ReplaceChar 
!
! !DESCRIPTION: Replaces all characters in a string 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE ReplaceChar ( str, pattern, replace )
!
! !INPUT PARAMETERS:
!
    CHARACTER(LEN=*), INTENT(IN)    :: pattern 
    CHARACTER(LEN=*), INTENT(IN)    :: replace
!                                                             
! !INPUT/OUTPUT PARAMETERS:                                         
!              
    CHARACTER(LEN=*), INTENT(INOUT) :: str 
!
! !REVISION HISTORY:
!  20 Dec 2018 - C. Keller   - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! LOCAL VARIABLES:
!
    INTEGER  :: I, LP, LS

    LP = LEN(TRIM(pattern))
    LS = LEN(TRIM(str))+1
    I  = INDEX(TRIM(str),TRIM(pattern))
    DO WHILE ( I > 0 ) 
       str = TRIM(str(1:I-1))//TRIM(replace)//TRIM(str(I+LP:LS))
       I = INDEX(TRIM(str),TRIM(pattern))
    ENDDO

  END SUBROUTINE ReplaceChar
!EOC
! GEOS-5 routine moved to gigc_providerservices_mod but needs updating so
! use this for now:
!------------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 910.1 and      !
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: FillAeroDP
!
! !DESCRIPTION: FillAeroDP fills the AERO_DP bundle
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE FillAeroDP ( am_I_Root, GC, EXPORT, RC ) 
!
! !USES:
!
    USE HCO_ERROR_MOD
    USE HCO_TYPES_MOD,     ONLY : DiagnCont
    USE HCO_DIAGN_MOD,     ONLY : Diagn_Get
    USE HCO_INTERFACE_MOD, ONLY : HcoState
!
! !INPUT PARAMETERS:
!
    LOGICAL                            :: am_I_Root
!                                                             
! !INPUT/OUTPUT PARAMETERS:                                         
!              
    TYPE(ESMF_GridComp), INTENT(INOUT) :: GC       ! Ref to this GridComp
    TYPE(ESMF_State),    INTENT(INOUT) :: Export   ! Export State
!                                                             
! !OUTPUT PARAMETERS:                                         
!              
    INTEGER, INTENT(OUT), OPTIONAL     :: RC
!
! !REVISION HISTORY:
!  30 Mar 2015 - C. Keller   - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! LOCAL VARIABLES:
!

    REAL, POINTER                :: Ptr2d(:,:) => NULL()
    INTEGER                      :: I, J, N, TrcID
    CHARACTER(LEN= 2)            :: Prfx 
    CHARACTER(LEN=15)            :: TrcName 
    CHARACTER(LEN=ESMF_MAXSTR)   :: ExpName 

    ! Hemco diagnostics
    INTEGER                      :: DgnID
    INTEGER                      :: FLAG, ERR
    TYPE(DiagnCont), POINTER     :: DgnCont => NULL()

    ! Error handling
    INTEGER                      :: STATUS
    CHARACTER(LEN=ESMF_MAXSTR)   :: Iam

    !=======================================================================
    ! FillAeroDP begins here
    !=======================================================================

    ! Traceback handle
    Iam = 'FillAeroDP'

    ! There are 8 species in total
    DO N = 1, 8
 
       ! Get species ID
       SELECT CASE ( N )
          CASE ( 1 )
             TrcName = 'DST1'
             Prfx    = 'DU'
          CASE ( 2 )
             TrcName = 'DST2'
             Prfx    = 'DU'
          CASE ( 3 )
             TrcName = 'DST3'
             Prfx    = 'DU'
          CASE ( 4 )
             TrcName = 'DST4'
             Prfx    = 'DU'
          CASE ( 5 )
             TrcName = 'BCPI'
             Prfx    = 'BC'
          CASE ( 6 )
             TrcName = 'BCPO'
             Prfx    = 'BC'
          CASE ( 7 )
             TrcName = 'OCPI'
             Prfx    = 'OC'
          CASE ( 8 )
             TrcName = 'OCPO'
             Prfx    = 'OC'
          CASE DEFAULT
             TrcName = 'YeahYeahYeah'
       END SELECT

       ! Get GEOS-Chem tracer ID
       TrcID = Ind_( TRIM(TrcName) )

       ! Only if tracer is defined...
       IF ( TrcID <= 0 ) CYCLE 

       ! Dry dep and wet dep
       DO I = 1, 2
   
          IF ( I == 1 ) THEN
             ExpName = TRIM(Prfx)//'DP_'//TRIM(TrcName)
          ELSEIF ( I == 2 ) THEN 
             ExpName = TRIM(Prfx)//'WT_'//TRIM(TrcName)
          ENDIF

          ! Get pointer
          CALL MAPL_GetPointer( EXPORT, Ptr2D, TRIM(ExpName),   &
                                notFoundOk=.TRUE., __RC__ )

          ! Skip if not defined
          IF ( .NOT. ASSOCIATED(Ptr2D) ) CYCLE             

          ! Reset
          Ptr2D = 0.0
      
          ! ------------------
          ! Dry deposition
          ! ------------------
          IF ( I == 1 ) THEN
            
             ! Get diagnostics 
             DgnID = 44500 + TrcID
             CALL Diagn_Get( am_I_Root, HcoState, .FALSE., DgnCont,  &
                             FLAG, ERR, cID=DgnID, AutoFill=-1,      &
                             COL=Input_Opt%DIAG_COLLECTION ) 

             ! Error check 
             _ASSERT( ERR == HCO_SUCCESS,'Error calling Diagn_Get' )

             ! Add to array if diagnostics is defined
             ! GEOS-Chem diagnostics is in kg m-2 s-1.
             IF ( FLAG == HCO_SUCCESS ) THEN
                IF ( ASSOCIATED(DgnCont%Arr2D%Val) ) THEN
                   Ptr2D = Ptr2D + DgnCont%Arr2D%Val
                ENDIF
             ENDIF
        
          ! ------------------
          ! Wet depostion 
          ! ------------------
          ELSEIF ( I == 2 ) THEN

             ! Convective and wet scavenging
             DO J = 1, 2

                SELECT CASE ( J ) 
                   ! Convection:
                   CASE ( 1 ) 
                      DgnID = 38000 + TrcID
                   ! Wet deposition
                   CASE ( 2 ) 
                      DgnID = 39000 + TrcID
                   CASE DEFAULT
                      DgnID = -1
                END SELECT

                ! Get diagnostics 
                CALL Diagn_Get( am_I_Root, HcoState, .FALSE., DgnCont,  &
                                FLAG, ERR, cID=DgnID, AutoFill=-1,      &
                                COL=Input_Opt%DIAG_COLLECTION ) 

                ! Error check 
                _ASSERT( ERR == HCO_SUCCESS,'Error calling Diagn_Get' )

                ! Add to array if diagnostics is defined. GEOS-Chem 
                ! diagnostics is already in kg m-2 s-1.
                IF ( FLAG == HCO_SUCCESS ) THEN
                   IF ( ASSOCIATED(DgnCont%Arr2D%Val) ) THEN
                      Ptr2D = Ptr2D + DgnCont%Arr2D%Val
                   ELSEIF ( ASSOCIATED(DgnCont%Arr3D%Val) ) THEN
                      Ptr2D = Ptr2D + SUM(DgnCont%Arr3D%Val,DIM=3)
                   ENDIF
                ENDIF
             ENDDO !J 
          ENDIF

       ENDDO !I
    ENDDO !N

    ! Successful return
    RC = ESMF_SUCCESS

  END SUBROUTINE FillAeroDP 
!EOC
!------------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 910.1 and      !
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: RoundOff
!
! !DESCRIPTION: Rounds a number X to N decimal places of precision.
!\\
!\\
! !INTERFACE:
!
  FUNCTION RoundOff( X, N ) RESULT( Y )
!
! !INPUT PARAMETERS:
! 
    REAL*8,  INTENT(IN) :: X   ! Number to be rounded
    INTEGER, INTENT(IN) :: N   ! Number of decimal places to keep
!
! !RETURN VALUE:
!
    REAL*8              :: Y   ! Number rounded to N decimal places
!
! !REMARKS:
!  The algorithm to round X to N decimal places is as follows:
!  (1) Multiply X by 10**(N+1)
!  (2) If X < 0, then add -5 to X; otherwise add 5 to X
!  (3) Round X to nearest integer
!  (4) Divide X by 10**(N+1)
!  (5) Truncate X to N decimal places: INT( X * 10**N ) / 10**N
!                                                                             .
!  Rounding algorithm from: Hultquist, P.F, "Numerical Methods for Engineers 
!   and Computer Scientists", Benjamin/Cummings, Menlo Park CA, 1988, p. 20.
!                                                                             .
!  Truncation algorithm from: http://en.wikipedia.org/wiki/Truncation
!                                                                             .
!  The two algorithms have been merged together for efficiency.
!
! !REVISION HISTORY:
!  14 Jul 2010 - R. Yantosca - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES
!
    ! Round and truncate X to N decimal places
    Y = INT( NINT( X*(10d0**(N+1)) + SIGN( 5d0, X ) ) / 10d0 ) / (10d0**N)

  END FUNCTION RoundOff
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: print_mean_oh
!
! !DESCRIPTION: Subroutine Print\_Mean\_OH prints the average mass-weighted OH 
!  concentration at the end of a simulation.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Print_Mean_OH( GC, State_Grid, RC )
!
! !USES:
!
    USE DIAG_OH_MOD,    ONLY : OH_MASS, AIR_MASS
    USE State_Grid_Mod, ONLY : GrdState
!
! !INPUT PARAMETERS:
!
    TYPE(ESMF_GridComp), INTENT(INOUT) :: GC             ! Ref to this GridComp
    TYPE(GrdState),      INTENT(IN)    :: State_Grid     ! Grid State obj
!
! !INPUT PARAMETERS:
!
    INTEGER,             INTENT(OUT)   :: RC             ! Return code
!
! !REMARKS:
!  The AIRMASS and OHMASS variables are used to pass the data values from the 
!  ESMF state, which is REAL*4.  We need to make sure that we do not store into
!  the ESMF state any values that exceed 1e38, which is the maximum allowable 
!  value.  Therefore, AIRMASS and OHMASS will store the values divided by the
!  scale factor OH_SCALE = 1d20 in order to prevent this overflow situation.
!
! !REVISION HISTORY: 
!  01 Jul 2010 - R. Yantosca - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES
!
    TYPE(MAPL_MetaComp), POINTER :: STATE          ! MAPL MetaComp object
    INTEGER                      :: LMAX
    INTEGER                      :: I,J
    CHARACTER(LEN=ESMF_MAXSTR)   :: compName      ! Name of gridded component
    REAL*8                       :: SUM_OHMASS   
    REAL*8                       :: SUM_MASS   
    REAL*8                       :: OHCONC
    REAL*8                       :: SUM_OHMASS_NH
    REAL*8                       :: SUM_MASS_NH
    REAL*8                       :: OHCONC_NH
    REAL*8                       :: SUM_OHMASS_SH
    REAL*8                       :: SUM_MASS_SH
    REAL*8                       :: OHCONC_SH
    REAL(ESMF_KIND_R4),  POINTER :: lonCtr(:,:) ! Lon centers on this PET [rad]
    REAL(ESMF_KIND_R4),  POINTER :: latCtr(:,:) ! Lat centers on this PET [rad]
    REAL,     POINTER            :: airMass(:,:,:)    ! Air mass [molec air]
    REAL,     POINTER            :: ohMass (:,:,:)    ! Mass-weighted OH
    REAL,     POINTER            :: airMass_nh(:,:,:) ! Air mass [molec air]
    REAL,     POINTER            :: ohMass_nh (:,:,:) ! Mass-weighted OH
    REAL,     POINTER            :: airMass_sh(:,:,:) ! Air mass [molec air]
    REAL,     POINTER            :: ohMass_sh (:,:,:) ! Mass-weighted OH
                                                    !  [molec OH/cm3 * 
                                                    !   molec air]
 
    __Iam__('Print_Mean_OH')

    ! Traceback info
    CALL ESMF_GridCompGet( GC, name=compName, __RC__ )
    Iam = TRIM(compName)//'::Print_Mean_OH'

    ! Get generic state object
    CALL MAPL_GetObjectFromGC( GC, STATE, __RC__ )

    ! Make sure diagnostics are allocated...
    IF ( ALLOCATED(AIR_MASS) .AND. ALLOCATED(OH_MASS) ) THEN

       ! Pointers to accumulated statistics. Seems like I have to cast from
       ! R8 to R here...
       ALLOCATE(airMass   (State_Grid%NX,State_Grid%NY,State_Grid%NZ))
       ALLOCATE(ohMass    (State_Grid%NX,State_Grid%NY,State_Grid%NZ)) 
       ALLOCATE(airMass_nh(State_Grid%NX,State_Grid%NY,State_Grid%NZ))
       ALLOCATE(ohMass_nh (State_Grid%NX,State_Grid%NY,State_Grid%NZ))
       ALLOCATE(airMass_sh(State_Grid%NX,State_Grid%NY,State_Grid%NZ))
       ALLOCATE(ohMass_sh (State_Grid%NX,State_Grid%NY,State_Grid%NZ))
       airMass    = 0.0
       ohMass     = 0.0
       airMass_nh = 0.0
       ohMass_nh  = 0.0
       airMass_sh = 0.0
       ohMass_sh  = 0.0
   
       ! OH diagnostics can have reduced vertical axis...
       LMAX = SIZE(OH_MASS,3)

       ! Pass to pointers.
       airMass(:,:,1:LMAX) = AIR_MASS(:,:,1:LMAX) / MAPL_AVOGAD
       ohMass (:,:,1:LMAX) = OH_MASS(:,:,1:LMAX)  / MAPL_AVOGAD
 
       ! Total Mass-weighted OH [molec OH/cm3] * [molec air]
       SUM_OHMASS = GlobalSum( GC, DataPtr3D=ohMass,  __RC__ )
   
       ! Atmospheric air mass [molec air]
       SUM_MASS   = GlobalSum( GC, DataPtr3D=airMass, __RC__ )

       ! Get the Orbit object (of type MAPL_SunOrbit),
       ! which is used in the call to MAPL_SunGetInsolation
       CALL MAPL_Get( STATE,                       &
                      LONS      = lonCtr,             &
                      LATS      = latCtr,             &
                      __RC__                         )

       ! SH/NH
       DO J=1,State_Grid%NY
       DO I=1,State_Grid%NX
          IF ( latCtr(I,J)/PI_180 >= 0.0 ) THEN
             airMass_nh(I,J,1:LMAX) = AIR_MASS(I,J,1:LMAX) / MAPL_AVOGAD
             ohMass_nh(I,J,1:LMAX)  = OH_MASS(I,J,1:LMAX) / MAPL_AVOGAD
          ELSE
             airMass_sh(I,J,1:LMAX) = AIR_MASS(I,J,1:LMAX) / MAPL_AVOGAD
             ohMass_sh(I,J,1:LMAX)  = OH_MASS(I,J,1:LMAX) / MAPL_AVOGAD
          ENDIF
       ENDDO
       ENDDO
 
       ! Atmospheric air mass
       SUM_OHMASS_NH = GlobalSum( GC, DataPtr3D=ohMass_nh, __RC__ )
       SUM_MASS_NH   = GlobalSum( GC, DataPtr3D=airMass_nh, __RC__ )

       SUM_OHMASS_SH = GlobalSum( GC, DataPtr3D=ohMass_sh, __RC__ )
       SUM_MASS_SH   = GlobalSum( GC, DataPtr3D=airMass_sh, __RC__ )
 
       ! Cleanup
       DEALLOCATE(airMass,ohMass)
       DEALLOCATE(airMass_nh,ohMass_nh)
       DEALLOCATE(airMass_sh,ohMass_sh)

    ! Will cause error msg below
    ELSE
       SUM_MASS = 0d0
    ENDIF   

    ! Avoid divide-by-zero errors 
    IF ( MAPL_am_I_Root() ) THEN
       IF ( SUM_MASS > 0d0 ) THEN 
               
          ! Divide OH by [molec air] and report as [1e5 molec/cm3]
          OHCONC    = ( SUM_OHMASS    / SUM_MASS    ) / 1d5
          OHCONC_NH = ( SUM_OHMASS_NH / SUM_MASS_NH ) / 1d5
          OHCONC_SH = ( SUM_OHMASS_SH / SUM_MASS_SH ) / 1d5
            
          ! Write value to log file
          WRITE( *, '(/,a)' ) REPEAT( '=', 79 ) 
          WRITE( *, *       ) 'GEOSCHEMchem mass-weighted OH concentration'
          WRITE( *, *       ) 'Mean OH = ', OHCONC, ' [1e5 molec/cm3]' 
          WRITE( *, *       ) 'OH NH/SH ratio: ', OHCONC_NH/OHCONC_SH
          WRITE( *, '(  a)' ) REPEAT( '=', 79 ) 
   
       ELSE
   
          ! Write error msg if SUM_MASS is zero
          WRITE( *, '(/,a)' ) REPEAT( '=', 79 ) 
          WRITE( *, '(  a)' )  &
                  'Could not print mass-weighted OH of GEOSCHEMchem!'
          WRITE( *, '(  a)' ) 'Atmospheric air mass is zero!'
          WRITE( *, '(  a)' ) REPEAT( '=', 79 ) 
          
       ENDIF
    ENDIF

    !=======================================================================
    ! All done
    !=======================================================================
    _RETURN(ESMF_SUCCESS)

  END SUBROUTINE Print_Mean_OH
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: GlobalSum
!
! !DESCRIPTION: Subroutine GlobalSum prints the sum (or minimum, or maximum)
!  of an array across all PETs.  Calls the ESMF_VMAllFullReduce function
!  to do the array reduction.  The default is to compute the array sum
!  unless MINIMUM=.TRUE. or MAXIMUM=.TRUE.
!\\
!\\
! !INTERFACE:
!
  FUNCTION GlobalSum( GC, dataPtr3D, dataPtr2D, minimum, maximum, RC ) &
                      RESULT( value )
!
! !INPUT PARAMETERS:
!
    TYPE(ESMF_GridComp),         INTENT(INOUT) :: GC               ! GC name
    REAL,    POINTER, OPTIONAL,  INTENT(IN)    :: dataPtr3D(:,:,:) ! Data array
    REAL,    POINTER, OPTIONAL,  INTENT(IN)    :: dataPtr2D(:,:)   ! Data array
    LOGICAL, OPTIONAL,           INTENT(IN)    :: minimum          ! Get min?
    LOGICAL, OPTIONAL,           INTENT(IN)    :: maximum          ! Get max?
!
! !OUTPUT PARAMETERS:
!
    INTEGER,                     INTENT(OUT)   :: RC         ! Return code
!
! !RETURN VALUE:
!
    REAL                                       :: value      ! Sum, max, or min
! 
! !REVISION HISTORY: 
!  01 Jul 2010 - R. Yantosca - Initial version
!  30 Nov 2015 - C. Keller   - Added 3D/2D option
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES
!
    ! Objects
    TYPE(ESMF_VM)              :: VM            ! ESMF VM object

    ! Scalars
    INTEGER                    :: nSize         ! Size of 3-D array
    CHARACTER(LEN=ESMF_MAXSTR) :: compName      ! Gridded component name
    __Iam__('GlobalSum')

    ! Arrays 
    REAL, ALLOCATABLE          :: data1d(:)     ! 1-D array for reduction

    !=======================================================================
    ! Initialization
    !=======================================================================

    ! Traceback info
    !CALL ESMF_GridCompGet( GC, name=compName, vm=VM, __RC__ )
    CALL ESMF_GridCompGet( GC, name=compName, __RC__ )
    CALL ESMF_VMGetCurrent(VM, __RC__ )
    Iam = TRIM(compName)//'::GlobalSum'

    !=======================================================================
    ! Do the reduction operation across all CPU's: sum, max, or min
    !=======================================================================

    ! Create a 1-D vector
    IF ( PRESENT( dataPtr3D ) ) THEN
       nSize = SIZE( dataPtr3D )
    ELSEIF ( PRESENT( dataPtr2D ) ) THEN
       nSize = SIZE( dataPtr2D )
    ENDIF
    ALLOCATE( data1d( nSize ), STAT=STATUS )
    _VERIFY(STATUS)
    data1d(:) = 0.0

    ! Rearrange into a 1-D array
    IF ( PRESENT( dataPtr3D ) ) THEN
       data1d = RESHAPE( dataPtr3D, (/1/) )
    ELSEIF ( PRESENT( dataPtr2D ) ) THEN
       data1d = RESHAPE( dataPtr2D, (/1/) )
    ENDIF
    
    ! Compute the sum over all PETS
    IF ( PRESENT( maximum ) ) THEN
       CALL ESMF_VMAllFullReduce( VM, data1d, value, nSize,  &
                                  ESMF_REDUCE_MAX, __RC__ )
    ELSE IF ( PRESENT( minimum ) ) THEN
       CALL ESMF_VMAllFullReduce( VM, data1d, value, nSize,  &
                                  ESMF_REDUCE_MIN, __RC__ )
    ELSE 
       CALL ESMF_VMAllFullReduce( VM, data1d, value, nSize,  &
                                  ESMF_REDUCE_SUM, __RC__ )
    ENDIF

    ! Deallocate temporary array
    IF( ALLOCATED( data1D ) ) DEALLOCATE( data1d )

  END FUNCTION GlobalSum
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: gridGetInterior
!
! !DESCRIPTION: Given an ESMF grid, returns the lower and upper longitude
!  and latitude indices on a given PET.
!\\
!\\
! !INTERFACE:
!
    SUBROUTINE GridGetInterior( Grid, I1, IN, J1, JN, RC )
!
! !INPUT PARAMETERS: 
!
    TYPE(ESMF_Grid), INTENT(IN)  :: Grid   ! ESMF Grid object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,         INTENT(OUT) :: I1     ! Lower lon index on this PET
    INTEGER,         INTENT(OUT) :: IN     ! Upper lon index on this PET
    INTEGER,         INTENT(OUT) :: J1     ! Lower lat index on this PET
    INTEGER,         INTENT(OUT) :: JN     ! Upper lat index on this PET
    INTEGER,         INTENT(OUT) :: RC     ! Success/failure
!
! !REMARKS:
!  This was a PRIVATE routine named MAPL_GridGetInterior within MAPL_Base.F90.
!  I have pulled the source code from there.
! 
! !REVISION HISTORY: 
!  30 Nov 2012 - R. Yantosca - Initial version, based on MAPL_Base
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    TYPE(ESMF_DistGrid)             :: distGrid
    TYPE(ESMF_DELayout)             :: LAYOUT
    INTEGER,            ALLOCATABLE :: AL(:,:)
    INTEGER,            ALLOCATABLE :: AU(:,:)
    INTEGER                         :: nDEs
    INTEGER                         :: deId
    INTEGER                         :: gridRank
    INTEGER                         :: deList(1)

    ! Identify this routine to MAPL
    __Iam__('GridGetInterior')

    ! Get ESMF DistGrid object
    CALL ESMF_GridGet    ( GRID,                           &
                           dimCount        = gridRank,     &
                           distGrid        = distGrid,     &
                           __RC__ )                    
 
    ! Get ESMF DELayout object
    CALL ESMF_DistGridGet( distGrid,                       &
                           delayout        = layout,       &
                           __RC__ )                    
                                                       
                
    ! Get the # of DE's and the list of DE's
    CALL ESMF_DELayoutGet( layout,                         &
                           deCount         = nDEs,         &
                           localDeList     = deList,       &
                           __RC__ )

    deId = deList(1)

    ! Allocate memory
    ALLOCATE( AL( gridRank, 0:nDEs-1 ), stat=status )
    ALLOCATE( AU( gridRank, 0:nDEs-1 ), stat=status )

    ! Get the min/max lon/lat values on each PET
    CALL ESMF_DistGridGet( distGrid,                       &
                           minIndexPDe = AL,               &
                           maxIndexPDe = AU,               &
                           __RC__ )

    ! Local Lon indices
    I1 = AL( 1, deId )   ! Lower
    IN = AU( 1, deId )   ! Upper

    ! Local lat indices
    J1 = AL( 2, deId )   ! Lower
    JN = AU( 2, deId )   ! Upper
 
    ! Free memory
    DEALLOCATE(AU, AL)
   
    ! Return successfully
    RC = GC_SUCCESS

  END SUBROUTINE GridGetInterior

  ! Adapted from the GOCART interface
  subroutine aerosol_optics(state, rc)

    implicit none

    ! Arguments
    ! ---------
    type(ESMF_State)     :: state
    integer, intent(out) :: rc


    ! Local
    ! ---------
    integer                                 :: n_aerosols
    character(len=ESMF_MAXSTR), allocatable :: aerosol_names(:)
    type(ESMF_FieldBundle)                  :: aerosols

    real, dimension(:,:,:), pointer         :: ple
    real, dimension(:,:,:), pointer         :: rh
    real, dimension(:,:,:), pointer         :: var
    real, dimension(:,:,:), pointer         :: q
    real, dimension(:,:,:,:), pointer       :: q_4d

    real, dimension(:,:,:), allocatable     :: dp, f_p

    character(len=ESMF_MAXSTR)              :: fld_name
    type(ESMF_Field)                        :: fld

    real, dimension(:,:,:,:), allocatable   :: ext, ssa, asy ! (lon:,lat:,lev:,band:)

    integer                                 :: n
    integer                                 :: i1, j1, i2, j2, km

    integer                                 :: band, offset

    integer                                 :: instance

    integer                                 :: STATUS
    character(len=ESMF_MAXSTR)              :: Iam

    integer, parameter                      :: n_bands = 1

    real    :: x
    integer :: i, j, k

    Iam = 'GEOSCHEMCHEM::aerosol_optics()'

    ! Mie Table instance/index
    ! ------------------------
    call ESMF_AttributeGet(state, name='mie_table_instance',  &
                           value=instance, __RC__)

    ! Radiation band
    ! --------------
    band = 0
    call ESMF_AttributeGet(state, name='band_for_aerosol_optics',  &
                           value=band, __RC__)
    offset = band - n_bands

    ! Pressure at layer edges 
    ! ------------------------
    call ESMF_AttributeGet(state, name='air_pressure_for_aerosol_optics', &
                           value=fld_name, __RC__)
    call MAPL_GetPointer(state, ple, trim(fld_name), __RC__)

    i1 = lbound(ple, 1); i2 = ubound(ple, 1)
    j1 = lbound(ple, 2); j2 = ubound(ple, 2)
    km = ubound(ple, 3)

    ! Relative humidity
    ! -----------------
    call ESMF_AttributeGet(state, name='relative_humidity_for_aerosol_optics', &
                           value=fld_name, __RC__)
    call MAPL_GetPointer(state, rh, trim(fld_name), __RC__)

    i1 = lbound(rh, 1); i2 = ubound(rh, 1)
    j1 = lbound(rh, 2); j2 = ubound(rh, 2)
    km = ubound(rh, 3)

    call ESMF_StateGet(state, 'AEROSOLS', aerosols, __RC__)
    call ESMF_FieldBundleGet(aerosols, fieldCount=n_aerosols, __RC__)

    allocate(aerosol_names(n_aerosols), __STAT__)

    call ESMF_FieldBundleGet(aerosols, FieldNameList=aerosol_names, __RC__)

    allocate(ext(i1:i2,j1:j2,km,n_bands), &
         ssa(i1:i2,j1:j2,km,n_bands), &
         asy(i1:i2,j1:j2,km,n_bands), __STAT__)

    allocate(q_4d(i1:i2,j1:j2,km,n_aerosols), __STAT__)

#if (0)
    allocate(dp(i1:i2,j1:j2,km), f_p(i1:i2,j1:j2,km), __STAT__)

    dp  = ple(:,:,1:km) - ple(:,:,0:km-1)
    f_p = dp / MAPL_GRAV

    do n = 1, n_aerosols
       call ESMF_FieldBundleGet(aerosols, trim(aerosol_names(n)),  &
                                field=fld, __RC__)
       call ESMF_FieldGet(fld, farrayPtr=q, __RC__)

       q_4d(:,:,:,n) = f_p * q
    end do

    call ESMF_AttributeGet(state, name='mie_table_instance',  &
                           value=instance, __RC__)
    call mie_(geoschemMieTable(instance), aerosol_names, n_bands, &
              offset, q_4d, rh, ext, ssa, asy, __RC__)

    deallocate(dp, f_p, __STAT__)
#else
    do n = 1, n_aerosols
       call ESMF_FieldBundleGet(aerosols, trim(aerosol_names(n)), &
                                field=fld, __RC__)
       call ESMF_FieldGet(fld, farrayPtr=q, __RC__)

       do k = 1, km
          do j = j1, j2
             do i = i1, i2
                x = ((PLE(i,j,k) - PLE(i,j,k-1))*0.01)*(100./MAPL_GRAV)
                q_4d(i,j,k,n) = x * q(i,j,k)
             end do
          end do
       end do
    end do

    call mie_(geoschemMieTable(instance), aerosol_names, n_bands,  &
              offset, q_4d, rh, ext, ssa, asy, __RC__)
#endif

    call ESMF_AttributeGet(state,                                            &
                           name='extinction_in_air_due_to_ambient_aerosol',  &
                           value=fld_name, __RC__)
    if (fld_name /= '') then 
       call MAPL_GetPointer(state, var, trim(fld_name), __RC__)
       var = ext(:,:,:,1)
    end if

    call ESMF_AttributeGet(state,                                             &
                           name='single_scattering_albedo_of_ambient_aerosol',&
                           value=fld_name, __RC__)
    if (fld_name /= '') then 
       call MAPL_GetPointer(state, var, trim(fld_name), __RC__)
       var = ssa(:,:,:,1)
    end if

    call ESMF_AttributeGet(state,                                         &
                           name='asymmetry_parameter_of_ambient_aerosol', &
                           value=fld_name, __RC__)
    if (fld_name /= '') then 
       call MAPL_GetPointer(state, var, trim(fld_name), __RC__)
       var = asy(:,:,:,1)
    end if

    deallocate(aerosol_names, ext, ssa, asy, q_4d, __STAT__)

    _RETURN(ESMF_SUCCESS)

  contains 

    subroutine mie_(mie_table, aerosol, nb, offset, q, rh, ext, ssa, asy, rc)

      implicit none

      type(Chem_Mie),    intent(inout):: mie_table    ! mie table
      character(len=*),  intent(in )  :: aerosol(:)   ! list of aerosols
      integer,           intent(in )  :: nb           ! number of bands
      integer,           intent(in )  :: offset       ! bands offset 
      real,              intent(in )  :: q(:,:,:,:)   ! aerosol mass mixing 
                                                      ! ratio, kg kg-1
      real,              intent(in )  :: rh(:,:,:)    ! relative humidity

      real,              intent(out)  :: ext(:,:,:,:) ! extinction
      real,              intent(out)  :: ssa(:,:,:,:) ! SSA
      real,              intent(out)  :: asy(:,:,:,:) ! asymmetry parameter

      integer,           intent(out)  :: rc

      ! local
      integer :: STATUS
      character(len=ESMF_MAXSTR) :: Iam='aerosol_optics::mie_' 

      integer :: l, idx, na

      real(kind=8) :: ext_(size(ext,1),size(ext,2),size(ext,3),size(ext,4))
      real(kind=8) :: ssa_(size(ext,1),size(ext,2),size(ext,3),size(ext,4))
      real(kind=8) :: asy_(size(ext,1),size(ext,2),size(ext,3),size(ext,4))

      na = size(aerosol)

      _ASSERT (na == size(q,4),'Error in number of aerosols')

      ext_ = 0.0d0
      ssa_ = 0.0d0
      asy_ = 0.0d0

      do l = 1, na
         idx = Chem_MieQueryIdx(mie_table, trim(aerosol(l)), __RC__)

         call Chem_MieQueryAllBand4D(mie_table, idx, nb, offset, &
                                     q(:,:,:,l), rh, ext, ssa, asy, __RC__)

         ext_ = ext_ +          ext     ! total extinction
         ssa_ = ssa_ +     (ssa*ext)    ! total scattering
         asy_ = asy_ + asy*(ssa*ext)    ! sum of (asy * sca)
      end do

      ext = ext_
      ssa = ssa_
      asy = asy_

      _RETURN(ESMF_SUCCESS)

    end subroutine mie_

  end subroutine aerosol_optics
!EOC
!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: aerodynamic_resistance --- calculates the aerodynamic resistance
!
! !INTERFACE:

   function aerodynamic_resistance(temperature, density_air, flux_sh, &
                                   friction_velocity, dz, z0h) result (r_a)
! !USES:

   implicit None

   real :: r_a                           ! Aerodynamic resistance, 'm-1 s'

! !INPUT/OUTPUT PARAMETERS:

! !INPUT PARAMETERS:

   real, intent(in) :: temperature       ! temperature,       'K'
   real, intent(in) :: density_air       ! density of air,    'kg m-3'
   real, intent(in) :: flux_sh           ! sensible heat flux at the surface, 
                                         ! 'W m-2'
   real, intent(in) :: friction_velocity ! friction velocity, 'm s-1'
   real, intent(in) :: dz                ! depth of the surface layer, 'm'
   real, intent(in) :: z0h               ! roughness height for sensible heat,
                                         ! 'm'

! !OUTPUT PARAMETERS:


! !DESCRIPTION: calculates the aerodynamic resistance (see Eq. XX, Seinfeld
!               and Pandis)
!
! !REVISION HISTORY:
!
!  15Dec2011  A. Darmenov
!
!EOP
!-------------------------------------------------------------------------

                    __Iam__('aerodynamic_resistance')

   ! local
   real, parameter :: k = MAPL_KARMAN 

   real :: z_ref
   real :: f
   real :: psi_h
   real :: log_f
   real :: z0h_
   real :: eps
   real :: L

   z_ref = 0.5 * dz

   L = monin_obukhov_length(temperature, density_air, flux_sh, &
                            friction_velocity)

   f = z_ref / L

   if(f > 1.0) then
       f = 1.0
   end if

   if ( (f > 0.0) .and. (f <= 1.0)) then
       psi_h = -5.0*f
   else if (f < 0.0) then
       eps = min(1.0, -f)
       log_f = log(eps)
       psi_h = exp(0.598 + 0.39*log_f - 0.09*(log_f**2))
   endif

   z0h_ = max(z0h, 1e2 * tiny(1.0))

   if ( friction_velocity > tiny(1.0) ) then 
      r_a = (log(z_ref / z0h_) - psi_h) / (k * friction_velocity)
   else
      r_a = 0.0
   endif

   end function aerodynamic_resistance
!EOC
!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: monin_obukhov_length --- calculates the Monin-Obukhov length
!
! !INTERFACE:

   function monin_obukhov_length(temperature, density_air, flux_sh, &
                                 friction_velocity) result (L)
! !USES:

   implicit None

   real :: L                             !  Monin-Obukhov length, 'm'

! !INPUT/OUTPUT PARAMETERS:

! !INPUT PARAMETERS:

   real, intent(in) :: temperature ! temperature,       'K'
   real, intent(in) :: density_air ! density of air,    'kg m-3'
   real, intent(in) :: flux_sh     ! sensible heat flux at the surface, 'W m-2'
   real, intent(in) :: friction_velocity ! friction velocity, 'm s-1'

! !OUTPUT PARAMETERS:


! !DESCRIPTION: calculates the calculates the Monin-Obukhov length
!
! !REVISION HISTORY:
!
!  15Dec2011  A. Darmenov
!
!EOP
!-------------------------------------------------------------------------

                    __Iam__('monin_obukhov_length')

   ! local
   real, parameter :: k   = MAPL_KARMAN 
   real, parameter :: g   = MAPL_GRAV
   real, parameter :: c_p = MAPL_CP

   if (abs(flux_sh) > 1e3*epsilon(0.0)) then
       L = - density_air * c_p * temperature * friction_velocity**3 / &
           (k * g * flux_sh)
   else
       L = 1/(1e3*epsilon(0.0))
   end if

   end function monin_obukhov_length

!EOC
#endif

#ifdef MODEL_GEOS
 END MODULE GEOSCHEMchem_GridCompMod
#else
 END MODULE Chem_GridCompMod
#endif
