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
     INTEGER                    :: nSpecies             ! # of species
     INTEGER                    :: nAdvect              ! # of advected species
     INTEGER                    :: nDryDep              ! # of drydep species
     INTEGER                    :: nKppSpc              ! # of KPP chem species
     INTEGER                    :: nWetDep              ! # of wetdep species

     !----------------------------------------------------------------------
     ! Mapping vectors to subset types of species
     !----------------------------------------------------------------------
     INTEGER,           POINTER :: Map_Advect (:      ) ! Advected species ID's
     INTEGER,           POINTER :: Map_DryDep (:      ) ! Drydep species ID's
     INTEGER,           POINTER :: Map_KppSpc (:      ) ! KPP chem species ID's
     INTEGER,           POINTER :: Map_WetDep (:      ) ! Wetdep species IDs'

     !----------------------------------------------------------------------
     ! Physical properties & indices for each species
     !----------------------------------------------------------------------
     TYPE(SpcPtr),      POINTER :: SpcData    (:      ) ! GC Species database

     !----------------------------------------------------------------------
     ! Chemical species
     !----------------------------------------------------------------------
     REAL(fp),          POINTER :: Species    (:,:,:,:) ! Species [molec/cm3]
     CHARACTER(LEN=20)          :: Spc_Units            ! Species units

     !----------------------------------------------------------------------
     ! Aerosol quantities
     !----------------------------------------------------------------------
     INTEGER                    :: nAero                ! # of Aerosol Types
     REAL(fp),          POINTER :: AeroArea   (:,:,:,:) ! Aerosol Area [cm2/cm3]
     REAL(fp),          POINTER :: AeroRadi   (:,:,:,:) ! Aerosol Radius [cm]
     REAL(fp),          POINTER :: WetAeroArea(:,:,:,:) ! Aerosol Area [cm2/cm3]
     REAL(fp),          POINTER :: WetAeroRadi(:,:,:,:) ! Aerosol Radius [cm]
     REAL(fp),          POINTER :: pHCloud    (:,:,:  ) ! Cloud pH [-]
     REAL(fp),          POINTER :: SSAlk      (:,:,:,:) ! Sea-salt alkalinity[-]

     !----------------------------------------------------------------------
     ! Fields for UCX mechanism
     !----------------------------------------------------------------------
     REAL(f4),          POINTER :: STATE_PSC  (:,:,:  ) ! PSC type (see Kirner
                                                        !  et al. 2011, GMD)
     REAL(fp),          POINTER :: KHETI_SLA  (:,:,:,:) ! Strat. liquid aerosol
                                                        !  reaction cofactors

     !----------------------------------------------------------------------
     ! For the tagged Hg simulation
     !----------------------------------------------------------------------
     INTEGER                    :: N_HG_CATS            ! # of Hg categories
     INTEGER,           POINTER :: Hg0_Id_List(:      ) ! Hg0 cat <-> tracer #
     INTEGER,           POINTER :: Hg2_Id_List(:      ) ! Hg2 cat <-> tracer #
     INTEGER,           POINTER :: HgP_Id_List(:      ) ! HgP cat <-> tracer #
     CHARACTER(LEN=4),  POINTER :: Hg_Cat_Name(:      ) ! Category names

     !----------------------------------------------------------------------
     ! For isoprene SOA
     !----------------------------------------------------------------------
     REAL(fp),          POINTER :: PH_SAV     (:,:,:  ) ! ISORROPIA aerosol pH
     REAL(fp),          POINTER :: HPLUS_SAV  (:,:,:  ) ! H+ concentration [M]
     REAL(fp),          POINTER :: WATER_SAV  (:,:,:  ) ! ISORROPIA aerosol H2O
     REAL(fp),          POINTER :: SULRAT_SAV (:,:,:  ) ! Sulfate conc [M]
     REAL(fp),          POINTER :: NARAT_SAV  (:,:,:  ) ! Nitrate conc [M]
     REAL(fp),          POINTER :: ACIDPUR_SAV(:,:,:  ) !
     REAL(fp),          POINTER :: BISUL_SAV  (:,:,:  ) ! Bisulfate conc [M]

     !----------------------------------------------------------------------
     ! For HOBr + S(IV) heterogeneous chemistry
     !----------------------------------------------------------------------
     REAL(fp),          POINTER :: HSO3_AQ    (:,:,:  ) ! Cloud bisulfite[mol/l]
     REAL(fp),          POINTER :: SO3_AQ     (:,:,:  ) ! Cloud sulfite  [mol/l]
     REAL(fp),          POINTER :: fupdateHOBr(:,:,:  ) ! Correction factor for
                                                        ! HOBr removal by SO2
                                                        ! [unitless]

     !----------------------------------------------------------------------
     ! Registry of variables contained within State_Chm
     !----------------------------------------------------------------------
     CHARACTER(LEN=4)           :: State     = 'CHEM'   ! Name of this state
     TYPE(MetaRegItem), POINTER :: Registry  => NULL()  ! Registry object

  END TYPE ChmState
!
! !REMARKS:
!                                                                             
! !REVISION HISTORY:
!  19 Oct 2012 - R. Yantosca - Initial version, based on "gc_type2_mod.F90"
!  26 Oct 2012 - R. Yantosca - Add fields for stratospheric chemistry
!  26 Feb 2013 - M. Long     - Add DEPSAV to derived type ChmState
!  07 Mar 2013 - R. Yantosca - Add Register_Tracer subroutine
!  07 Mar 2013 - R. Yantosca - Now make POSITION a locally SAVEd variable
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
!  19 May 2014 - C. Keller   - Removed Trac_Btend. DepSav array covers now
!                              all species.
!  03 Dec 2014 - M. Yannetti - Added PRECISION_MOD
!  11 Dec 2014 - R. Yantosca - Keep JLOP and JLOP_PREV for ESMF runs only
!  17 Feb 2015 - E. Lundgren - New tracer units kg/kg dry air (previously kg)
!  13 Aug 2015 - E. Lundgren - Add tracer units string to ChmState derived type 
!  28 Aug 2015 - R. Yantosca - Remove strat chemistry fields, these are now
!                              handled by the HEMCO component
!  05 Jan 2016 - E. Lundgren - Use global physical constants
!  28 Jan 2016 - M. Sulprizio- Add STATE_PSC, KHETI_SLA. These were previously
!                              local arrays in ucx_mod.F, but now need to be
!                              accessed in gckpp_HetRates.F90.
!  12 May 2016 - M. Sulprizio- Add WetAeroArea, WetAeroRadi to replace 1D arrays
!                              WTARE, WERADIUS previously in comode_mod.F
!  18 May 2016 - R. Yantosca - Add mapping vectors for subsetting species
!  07 Jun 2016 - M. Sulprizio- Remove routines Get_Indx, Register_Species, and
!                              Register_Tracer made obsolete by the species
!                              database.
!  22 Jun 2016 - R. Yantosca - Rename Id_Hg0 to Hg0_Id_List, Id_Hg2 to
!                              Hg2_Id_List, and Id_HgP to HgP_Id_List
!  16 Aug 2016 - M. Sulprizio- Rename from gigc_state_chm_mod.F90 to
!                              state_chm_mod.F90. The "gigc" nomenclature is
!                              no longer used.
!  23 Aug 2016 - M. Sulprizio- Remove tracer fields from State_Chm. These are
!                              now entirely replaced with the species fields.
!  08 Jun 2017 - M. Sulprizio- Add fields for isoprene SOA updates from E.Marais
!  29 Jun 2017 - R. Yantosca - Add fields of State_Chm to the registry
!  29 Jun 2017 - R. Yantosca - Remove Spec_Id, it's no longer used
!  30 Jun 2017 - R. Yantosca - Now register variables of State_Chm
!  31 Jul 2017 - R. Yantosca - Add fixes in registering ISORROPIA *_SAV fields
!  26 Sep 2017 - E. Lundgren - Remove Lookup_State_Chm and Print_State_Chm
!  02 Oct 2017 - E. Lundgren - Abstract metadata and routine to add to Registry
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !MODULE INTERFACES:
!
  INTERFACE Register_ChmField
     MODULE PROCEDURE Register_ChmField_R4_3D
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
  SUBROUTINE Init_State_Chm( am_I_Root, IM,        JM,         &   
                             LM,        Input_Opt, State_Chm,  &
                             nAerosol,  RC                    )
!
! !USES:
!
    USE GCKPP_Parameters,     ONLY : NSPEC
    USE Input_Opt_Mod,        ONLY : OptInput
    USE Species_Database_Mod, ONLY : Init_Species_Database
!
! !INPUT PARAMETERS:
! 
    LOGICAL,        INTENT(IN)    :: am_I_Root   ! Is this the root CPU?
    INTEGER,        INTENT(IN)    :: IM          ! # longitudes on this PET
    INTEGER,        INTENT(IN)    :: JM          ! # longitudes on this PET
    INTEGER,        INTENT(IN)    :: LM          ! # longitudes on this PET
    INTEGER,        INTENT(IN)    :: nAerosol    ! # aerosol species
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
!  19 Oct 2012 - R. Yantosca - Now pass all dimensions as arguments
!  26 Oct 2012 - R. Yantosca - Now allocate Strat_P, Strat_k fields
!  26 Oct 2012 - R. Yantosca - Add nSchem, nSchemBry as arguments
!  01 Nov 2012 - R. Yantosca - Don't allocate strat chem fields if nSchm=0
!                              and nSchmBry=0 (i.e. strat chem is turned off)
!  26 Feb 2013 - M. Long     - Now pass Input_Opt via the argument list
!  26 Feb 2013 - M. Long     - Now allocate the State_Chm%DEPSAV field
!  11 Dec 2014 - R. Yantosca - Remove TRAC_TEND and DEPSAV fields
!  13 Aug 2015 - E. Lundgren - Initialize trac_units to ''
!  28 Aug 2015 - R. Yantosca - Remove stratospheric chemistry fields; 
!                              these are all now read in via HEMCO
!  28 Aug 2015 - R. Yantosca - Also initialize the species database object
!  09 Oct 2015 - R. Yantosca - Bug fix: set State_Chm%SpcData to NULL
!  16 Dec 2015 - R. Yantosca - Now overwrite the Input_Opt%TRACER_MW_G and
!                              related fields w/ info from species database
!  29 Apr 2016 - R. Yantosca - Don't initialize pointers in declaration stmts
!  02 May 2016 - R. Yantosca - Nullify Hg index fields for safety's sake
!  18 May 2016 - R. Yantosca - Now determine the # of each species first,
!                              then allocate fields of State_Chm
!  18 May 2016 - R. Yantosca - Now populate the species mapping vectors
!  30 Jun 2016 - M. Sulprizio- Remove nSpecies as an input argument. This is now
!                              initialized as the size of SpcData.
!  22 Jul 2016 - E. Lundgren - Initialize spc_units to ''
!  28 Nov 2016 - R. Yantosca - Only allocate STATE_PSC and KHETI_SLA for UCX
!                              simulations; set to NULL otherwise
!  28 Nov 2016 - R. Yantosca - Only allocate State_Chm%*Aero* fields for
!                              fullchem and/or aerosol-only simulations
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER                :: N,          C
    INTEGER                :: N_Hg0_CATS, N_Hg2_CATS, N_HgP_CATS
    INTEGER                :: nKHLSA

    ! Strings
    CHARACTER(LEN=255)     :: ErrMsg, ThisLoc, ChmID

    ! Pointers
    TYPE(Species), POINTER :: ThisSpc
    REAL(fp),      POINTER :: Ptr2data(:,:,:)

    ! Error handling
    RC = GC_SUCCESS
    ErrMsg = ''
    ThisLoc = ' -> at Init_State_Chm (in Headers/state_chm_mod.F90)'

    !=====================================================================
    ! Initialization
    !=====================================================================

    ! Number of each type of species
    State_Chm%nSpecies    =  0
    State_Chm%nAdvect     =  0
    State_Chm%nDryDep     =  0
    State_Chm%nKppSpc     =  0
    State_Chm%nWetDep     =  0

    ! Mapping vectors for subsetting each type of species
    State_Chm%Map_Advect  => NULL()
    State_Chm%Map_DryDep  => NULL()
    State_Chm%Map_KppSpc  => NULL()
    State_Chm%Map_WetDep  => NULL()

    ! Chemical species
    State_Chm%Species     => NULL()
    State_Chm%Spc_Units   = ''

    ! Species database
    State_Chm%SpcData     => NULL()
    ThisSpc               => NULL()

    ! Aerosol parameters
    State_Chm%nAero       = 0
    State_Chm%AeroArea    => NULL()
    State_Chm%AeroRadi    => NULL()
    State_Chm%WetAeroArea => NULL()
    State_Chm%WetAeroRadi => NULL()

    ! Fields for UCX mechanism
    State_Chm%STATE_PSC   => NULL()
    State_Chm%KHETI_SLA   => NULL()   

    ! pH/alkalinity
    State_Chm%pHCloud     => NULL()
    State_Chm%SSAlk       => NULL()

    ! Hg species indexing
    N_Hg0_CATS            =  0
    N_Hg2_CATS            =  0
    N_HgP_CATS            =  0
    State_Chm%N_Hg_CATS   =  0
    State_Chm%Hg_Cat_Name => NULL()
    State_Chm%Hg0_Id_List => NULL()
    State_Chm%Hg2_Id_List => NULL()
    State_Chm%HgP_Id_List => NULL()

    ! For isoprene SOA
    State_Chm%PH_SAV      => NULL()
    State_Chm%HPLUS_SAV   => NULL()
    State_Chm%WATER_SAV   => NULL()
    State_Chm%SULRAT_SAV  => NULL()
    State_Chm%NARAT_SAV   => NULL()
    State_Chm%ACIDPUR_SAV => NULL()
    State_Chm%BISUL_SAV   => NULL()

    ! For HOBr + S(IV) chemistry
    State_Chm%HSO3_AQ     => NULL()
    State_Chm%SO3_AQ      => NULL()
    State_Chm%fupdateHOBr => NULL()

    ! Local variables
    Ptr2data => NULL()

    !=====================================================================
    ! Populate the species database object field
    ! (assumes Input_Opt has already been initialized)
    !=====================================================================
    CALL Init_Species_Database( am_I_Root = am_I_Root,          &
                                Input_Opt = Input_Opt,          &
                                SpcData   = State_Chm%SpcData,  &
                                RC        = RC                 )
    SpcDataLocal => State_Chm%SpcData

    !=====================================================================
    ! Determine the number of advected, drydep, wetdep, and total species
    !=====================================================================

    ! The total number of species is the size of SpcData
    State_Chm%nSpecies = SIZE( State_Chm%SpcData )

    ! Get the number of advected, dry-deposited, KPP chemical species,
    ! and and wet-deposited species.  Also return the # of Hg0, Hg2, and 
    ! HgP species (these are zero unless the Hg simulation is used).
    CALL Spc_GetNumSpecies( State_Chm%nAdvect,  &
                            State_Chm%nDryDep,  &
                            State_Chm%nKppSpc,  &
                            State_Chm%nWetDep,  &
                            N_Hg0_CATS,         &
                            N_Hg2_CATS,         &
                            N_HgP_CATS         )

    !=====================================================================
    ! Allocate and initialize mapping vectors to subset species
    !=====================================================================

    ALLOCATE( State_Chm%Map_Advect( State_Chm%nAdvect  ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%Map_Advect', 0, RC )  
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Chm%Map_Advect = 0

    ALLOCATE( State_Chm%Map_Drydep( State_Chm%nDryDep  ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%Map_Drydep', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Chm%Map_DryDep = 0

    ALLOCATE( State_Chm%Map_KppSpc( NSPEC ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%Map_KppSpc', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Chm%Map_KppSpc = 0

    ALLOCATE( State_Chm%Map_WetDep( State_Chm%nWetDep  ), STAT=RC )
    CALL GC_CheckVar( 'State_Met%Map_Wetdep', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Chm%Map_WetDep = 0

    !=======================================================================
    ! Set up the species mapping vectors
    !=======================================================================
    IF ( am_I_Root ) THEN
       WRITE( 6,'(/,a)' ) 'ADVECTED SPECIES MENU'
       WRITE( 6,'(  a)' ) REPEAT( '-', 48 )
       WRITE( 6,'(  a)' ) '  # Species Name  g/mole'
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
          C                         = ThisSpc%AdvectId
          State_Chm%Map_Advect(C)   = ThisSpc%ModelId
          
          ! Print to screen
          IF ( am_I_Root ) THEN
             WRITE( 6, 100 ) ThisSpc%ModelId, ThisSpc%Name, ThisSpc%MW_g
          ENDIF

       ENDIF

       !--------------------------------------------------------------------
       ! Set up the mapping for DRYDEP SPECIES
       !--------------------------------------------------------------------
       IF ( ThisSpc%Is_DryDep ) THEN
          C                         = ThisSpc%DryDepId
          State_Chm%Map_Drydep(C)   = ThisSpc%ModelId
       ENDIF

       !--------------------------------------------------------------------
       ! Set up the mapping for SPECIES IN THE KPP MECHANISM
       !--------------------------------------------------------------------
       IF ( ThisSpc%Is_Kpp ) THEN
          C                         = ThisSpc%KppSpcId
          State_Chm%Map_KppSpc(C)   = ThisSpc%ModelId
       ENDIF

       !--------------------------------------------------------------------
       ! Set up the mapping for SPECIES IN THE KPP MECHANISM
       !--------------------------------------------------------------------
       IF ( ThisSpc%Is_WetDep ) THEN
          C                         = ThisSpc%WetDepId
          State_Chm%Map_WetDep(C)   = ThisSpc%ModelId
       ENDIF

       ! Free pointer
       ThisSpc => NULL()

    ENDDO

    !=====================================================================
    ! Allocate and initialize chemical species fields
    !=====================================================================   
    chmID = 'SPC'
    ALLOCATE( State_Chm%Species( IM, JM, LM, State_Chm%nSpecies ), STAT=RC )
    CALL GC_CheckVar( 'State_Chm%Species', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Chm%Species = 0.0_fp
    CALL Register_ChmField( am_I_Root, chmID, State_Chm%Species, State_Chm, RC )

    !=====================================================================
    ! Allocate and initialize aerosol area and radius fields
    ! These are only relevant for fullchem or aerosol-only simulations
    !=====================================================================
    IF ( Input_Opt%ITS_A_FULLCHEM_SIM .or. Input_Opt%ITS_AN_AEROSOL_SIM ) THEN

       ! Save nAerosol to State_Chm
       State_Chm%nAero = nAerosol

       !------------------------------------------------------------------
       ! AEROAREA
       !------------------------------------------------------------------
       ALLOCATE( State_Chm%AeroArea( IM, JM, LM, State_Chm%nAero ), STAT=RC )
       CALL GC_CheckVar( 'State_Chm%AeroArea', 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Chm%AeroArea = 0.0_fp

       ! Loop over all entries to register each category individually
       DO N = 1, State_Chm%nAero

          ! Define identifying string
          SELECT CASE( N )
             CASE( 1  )
                chmID = 'AEROAREA_MDUST1'
             CASE( 2  )
                chmID = 'AEROAREA_MDUST2'
             CASE( 3  )
                chmID = 'AEROAREA_MDUST3'
             CASE( 4  )
                chmID = 'AEROAREA_MDUST4'
             CASE( 5  )
                chmID = 'AEROAREA_MDUST5'
             CASE( 6  )
                chmID = 'AEROAREA_MDUST6'
             CASE( 7  )
                chmID = 'AEROAREA_MDUST7'
             CASE( 8  )
                chmID = 'AEROAREA_SULF'
             CASE( 9  )
                chmID = 'AEROAREA_BC'
             CASE( 10 )
                chmID = 'AEROAREA_OC'
             CASE( 11 )
                chmID = 'AEROAREA_SSA'
             CASE( 12 )
                chmID = 'AEROAREA_SSC'
             CASE( 13 )
                chmID = 'AEROAREA_BGSULF'
             CASE( 14 )
                chmID  = 'AEROAREA_ICEI'
             CASE DEFAULT
                ErrMsg = 'State_Chm%nAero exceeds the number of defined' &
                         // ' dry aerosol area categories'
                CALL GC_Error( ErrMsg, RC, ThisLoc )
                RETURN
          END SELECT

          CALL Register_ChmField( am_I_Root, chmID, State_Chm%AeroArea, &
                                  State_Chm, RC,    Ncat=N )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDDO

       !------------------------------------------------------------------
       ! AERORADI
       !------------------------------------------------------------------
       ALLOCATE( State_Chm%AeroRadi( IM, JM, LM, State_Chm%nAero ), STAT=RC )
       CALL GC_CheckVar( 'State_Chm%AeroRadi', 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Chm%AeroRadi    = 0.0_fp

       ! Loop over all entries to register each category individually
       DO N = 1, State_Chm%nAero

          ! Define identifying string
          SELECT CASE( N )
             CASE( 1  )
                chmID = 'AERORADI_MDUST1'
             CASE( 2  )
                chmID = 'AERORADI_MDUST2'
             CASE( 3  )
                chmID = 'AERORADI_MDUST3'
             CASE( 4  )
                chmID = 'AERORADI_MDUST4'
             CASE( 5  )
                chmID = 'AERORADI_MDUST5'
             CASE( 6  )
                chmID = 'AERORADI_MDUST6'
             CASE( 7  )
                chmID = 'AERORADI_MDUST7'
             CASE( 8  )
                chmID = 'AERORADI_SULF'
             CASE( 9  )
                chmID = 'AERORADI_BC'
             CASE( 10 )
                chmID = 'AERORADI_OC'
             CASE( 11 )
                chmID = 'AERORADI_SSA'
             CASE( 12 )
                chmID = 'AERORADI_SSC'
             CASE( 13 )
                chmID = 'AERORADI_BGSULF'
             CASE( 14 )
                chmID = 'AERORADI_ICEI'
             CASE DEFAULT
                ErrMsg = 'State_Chm%nAero exceeds the number of defined' &
                         // ' dry aerosol radius categories'
                CALL GC_Error( ErrMsg, RC, ThisLoc )
                RETURN
          END SELECT

          CALL Register_ChmField( am_I_Root, chmID, State_Chm%AeroRadi, &
                                  State_Chm, RC,    Ncat=N )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDDO

       !------------------------------------------------------------------
       ! WETAEROAREA
       !------------------------------------------------------------------
       ALLOCATE( State_Chm%WetAeroArea( IM, JM, LM, State_Chm%nAero ), STAT=RC )
       CALL GC_CheckVar( 'State_Chm%WetAeroArea', 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Chm%WetAeroArea = 0.0_fp

       ! Loop over all entries to register each category individually
       DO N = 1, State_Chm%nAero

          ! Define identifying string
          SELECT CASE( N )
             CASE( 1  )
                chmID = 'WETAEROAREA_MDUST1'
             CASE( 2  )
                chmID = 'WETAEROAREA_MDUST2'
             CASE( 3  )
                chmID = 'WETAEROAREA_MDUST3'
             CASE( 4  )
                chmID = 'WETAEROAREA_MDUST4'
             CASE( 5  )
                chmID = 'WETAEROAREA_MDUST5'
             CASE( 6  )
                chmID = 'WETAEROAREA_MDUST6'
             CASE( 7  )
                chmID = 'WETAEROAREA_MDUST7'
             CASE( 8  )
                chmID = 'WETAEROAREA_SULF'
             CASE( 9  )
                chmID = 'WETAEROAREA_BC'
             CASE( 10 )
                chmID = 'WETAEROAREA_OC'
             CASE( 11 )
                chmID = 'WETAEROAREA_SSA'
             CASE( 12 )
                chmID = 'WETAEROAREA_SSC'
             CASE( 13 )
                chmID = 'WETAEROAREA_BGSULF'
             CASE( 14 )
                chmID = 'WETAEROAREA_ICEI'
             CASE DEFAULT
                ErrMsg = 'State_Chm%nAero exceeds the number of defined' &
                         // ' wet aerosol area categories'
                CALL GC_Error( ErrMsg, RC, ThisLoc )
                RETURN
          END SELECT

          CALL Register_ChmField( am_I_Root, chmID, State_Chm%WetAeroArea, &
                                  State_Chm, RC,    Ncat=N )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDDO

       !------------------------------------------------------------------
       ! WETAERORADI
       !------------------------------------------------------------------
       ALLOCATE( State_Chm%WetAeroRadi( IM, JM, LM, State_Chm%nAero ), STAT=RC )
       CALL GC_CheckVar( 'State_Chm%WetAeroRadi', 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Chm%WetAeroRadi = 0.0_fp

       ! Loop over all entries to register each category individually
       DO N = 1, State_Chm%nAero

          ! Define identifying string
          SELECT CASE( N )
             CASE( 1  )
                chmID = 'WETAERORADI_MDUST1'
             CASE( 2  )
                chmID = 'WETAERORADI_MDUST2'
             CASE( 3  )
                chmID = 'WETAERORADI_MDUST3'
             CASE( 4  )
                chmID = 'WETAERORADI_MDUST4'
             CASE( 5  )
                chmID = 'WETAERORADI_MDUST5'
             CASE( 6  )
                chmID = 'WETAERORADI_MDUST6'
             CASE( 7  )
                chmID = 'WETAERORADI_MDUST7'
             CASE( 8  )
                chmID = 'WETAERORADI_SULF'
             CASE( 9  )
                chmID = 'WETAERORADI_BC'
             CASE( 10 )
                chmID = 'WETAERORADI_OC'
             CASE( 11 )
                chmID = 'WETAERORADI_SSA'
             CASE( 12 )
                chmID = 'WETAERORADI_SSC'
             CASE( 13 )
                chmID = 'WETAERORADI_BGSULF'
             CASE( 14 )
                chmID = 'WETAERORADI_ICEI'
             CASE DEFAULT
                ErrMsg = 'State_Chm%nAero exceeds the number of defined' &
                         // ' wet aerosol radius categories'
                CALL GC_Error( ErrMsg, RC, ThisLoc )
                RETURN
          END SELECT          

          CALL Register_ChmField( am_I_Root, chmID, State_Chm%WetAeroRadi, &
                                  State_Chm, RC,    Ncat=N )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDDO

    ENDIF

    !=====================================================================
    ! Allocate and initialize fields for halogen chemistry
    !=====================================================================
    
    ALLOCATE( State_Chm%pHCloud    ( IM, JM, LM                    ), STAT=RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Chm%pHCloud = 0e+0_fp

    ALLOCATE( State_Chm%SSAlk      ( IM, JM, LM, 2                 ), STAT=RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Chm%SSAlk = 0e+0_fp

    !=====================================================================
    ! Allocate and initialize fields for UCX mechamism
    !=====================================================================
#if defined( UCX )

    !---------------------------------------------------------------------
    ! STATE_PSC
    !---------------------------------------------------------------------
    chmID = 'STATE_PSC'
    ALLOCATE( State_Chm%STATE_PSC( IM, JM, LM ), STAT=RC )
    CALL GC_CheckVar( 'State_Chm%STATE_PSC', 0, RC )    
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Chm%STATE_PSC = 0.0_f4
    CALL Register_ChmField( am_I_Root, chmID, State_Chm%STATE_PSC, &
                            State_Chm, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    !---------------------------------------------------------------------
    ! KHETI_SLA
    !---------------------------------------------------------------------
    nKHLSA = 11 ! TODO: should this be in CMN_SIZE_Mod?
    ALLOCATE( State_Chm%KHETI_SLA ( IM, JM, LM, nKHLSA ), STAT=RC )
    CALL GC_CheckVar( 'State_Chm%KHETI_SLA', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    State_Chm%KHETI_SLA = 0.0_fp

    ! Loop over all entries to register each category individually
    DO N = 1, nKHLSA

       ! Define identifying string
       SELECT CASE( N )
          CASE( 1  ) 
             chmID = 'KHSLA_N2O5+H2O'
          CASE( 2  ) 
             chmID = 'KHSLA_N2O5+HCl'
          CASE( 3  ) 
             chmID = 'KHSLA_ClNO3+H2O'
          CASE( 4  ) 
             chmID = 'KHSLA_ClNO3+HCl'
          CASE( 5  ) 
             chmID = 'KHSLA_ClNO3+HBr'
          CASE( 6  ) 
             chmID = 'KHSLA_BrNO3+H2O'
          CASE( 7  ) 
             chmID = 'KHSLA_BrNO3+HCl'
          CASE( 8  ) 
             chmID = 'KHSLA_HOCl+HCl'
          CASE( 9  ) 
             chmID = 'KHSLA_HOCl+HBr'
          CASE( 10 ) 
             chmID = 'KHSLA_HOBr+HCl'
          CASE( 11 ) 
             chmID = 'KHSLA_HOBr+HBr'
          CASE DEFAULT
             ErrMsg = 'nKHLSA exceeds the number of defined' &
                      // ' KHETI_SLA categories'
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN
       END SELECT

       CALL Register_ChmField( am_I_Root, chmID, State_Chm%KHETI_SLA, &
                               State_Chm, RC,    Ncat=N )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDDO
#endif

    !=====================================================================
    ! Allocate and initialize isoprene SOA fields
    ! These are only relevant for fullchem or aerosol-only simulations
    !=====================================================================
    IF ( Input_Opt%ITS_A_FULLCHEM_SIM .or. Input_Opt%ITS_AN_AEROSOL_SIM ) THEN

       !------------------------------------------------------------------
       ! PH_SAV
       !------------------------------------------------------------------
       chmID = 'PH_SAV'
       ALLOCATE( State_Chm%PH_SAV( IM, JM, LM ) , STAT=RC )
       CALL GC_CheckVar( 'State_Chm%PH_SAV', 0, RC )    
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Chm%PH_SAV = 0.0_fp
       CALL Register_ChmField( am_I_Root, chmID, State_Chm%PH_SAV, &
                               State_Chm, RC )
       IF ( RC /= GC_SUCCESS ) RETURN

       !------------------------------------------------------------------
       ! HPLUS_SAV
       !------------------------------------------------------------------
       chmID = 'HPLUS_SAV'
       ALLOCATE( State_Chm%HPLUS_SAV( IM, JM, LM ) , STAT=RC )
       CALL GC_CheckVar( 'State_Chm%HPLUS_SAV', 0, RC )    
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Chm%HPLUS_SAV = 0.0_fp
       CALL Register_ChmField( am_I_Root, chmID, State_Chm%HPLUS_SAV, &
                               State_Chm, RC )
       IF ( RC /= GC_SUCCESS ) RETURN

       !------------------------------------------------------------------
       ! WATER_SAV
       !------------------------------------------------------------------
       chmID = 'WATER_SAV'
       ALLOCATE( State_Chm%WATER_SAV( IM, JM, LM ) , STAT=RC )
       CALL GC_CheckVar( 'State_Chm%WATER_SAV', 0, RC )    
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Chm%WATER_SAV = 0.0_fp
       CALL Register_ChmField( am_I_Root, chmID, State_Chm%WATER_SAV, &
                               State_Chm, RC )
       IF ( RC /= GC_SUCCESS ) RETURN

       !------------------------------------------------------------------
       ! SULRAT_SAV
       !------------------------------------------------------------------
       chmID = 'SULRAT_SAV'
       ALLOCATE( State_Chm%SULRAT_SAV( IM, JM, LM ) , STAT=RC )
       CALL GC_CheckVar( 'State_Chm%SULRAT_SAV', 0, RC )    
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Chm%SULRAT_SAV = 0.0_fp
       CALL Register_ChmField( am_I_Root, chmID, State_Chm%SULRAT_SAV, &
                               State_Chm, RC )
       IF ( RC /= GC_SUCCESS ) RETURN

       !------------------------------------------------------------------
       ! NARAT_SAV
       !------------------------------------------------------------------
       chmID = 'NARAT_SAV'
       ALLOCATE( State_Chm%NARAT_SAV( IM, JM, LM ) , STAT=RC )
       CALL GC_CheckVar( 'State_Chm%NARAT_SAV', 0, RC )    
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Chm%NARAT_SAV = 0.0_fp
       CALL Register_ChmField( am_I_Root, chmID, State_Chm%NARAT_SAV, &
                               State_Chm, RC )
       IF ( RC /= GC_SUCCESS ) RETURN

       !------------------------------------------------------------------
       ! ACIDPUR_SAV
       !------------------------------------------------------------------
       chmID = 'ACIDPUR_SAV'
       ALLOCATE( State_Chm%ACIDPUR_SAV( IM, JM, LM ) , STAT=RC )
       CALL GC_CheckVar( 'State_Chm%ACIDPUR_SAV', 0, RC )    
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Chm%ACIDPUR_SAV = 0.0_fp
       CALL Register_ChmField( am_I_Root, chmID, State_Chm%ACIDPUR_SAV, &
                               State_Chm, RC )
       IF ( RC /= GC_SUCCESS ) RETURN

       !------------------------------------------------------------------
       ! BISUL_SAV
       !------------------------------------------------------------------
       chmID = 'BISUL_SAV'
       ALLOCATE( State_Chm%BISUL_SAV( IM, JM, LM ) , STAT=RC )
       CALL GC_CheckVar( 'State_Chm%BISUL_SAV', 0, RC )    
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Chm%BISUL_SAV = 0.0_fp
       CALL Register_ChmField( am_I_Root, chmID, State_Chm%BISUL_SAV, &
                               State_Chm, RC )
       IF ( RC /= GC_SUCCESS ) RETURN

       !------------------------------------------------------------------
       ! HSO3_AQ
       !------------------------------------------------------------------
       chmID = 'HSO3_AQ'
       ALLOCATE( State_Chm%HSO3_AQ( IM, JM, LM ) , STAT=RC )
       CALL GC_CheckVar( 'State_Chm%HSO3_AQ', 0, RC )    
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Chm%HSO3_AQ = 0.0_fp
       CALL Register_ChmField( am_I_Root, chmID, State_Chm%HSO3_AQ, &
                               State_Chm, RC )
       IF ( RC /= GC_SUCCESS ) RETURN

       !------------------------------------------------------------------
       ! SO3_AQ
       !------------------------------------------------------------------
       chmID = 'SO3_AQ'
       ALLOCATE( State_Chm%SO3_AQ( IM, JM, LM ) , STAT=RC )
       CALL GC_CheckVar( 'State_Chm%SO3_AQ', 0, RC )    
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Chm%SO3_AQ = 0.0_fp
       CALL Register_ChmField( am_I_Root, chmID, State_Chm%SO3_AQ, &
                               State_Chm, RC )
       IF ( RC /= GC_SUCCESS ) RETURN

       !------------------------------------------------------------------
       ! fupdateHOBr
       !------------------------------------------------------------------
       chmID = 'fupdateHOBr'
       ALLOCATE( State_Chm%fupdateHOBr( IM, JM, LM ) , STAT=RC )
       CALL GC_CheckVar( 'State_Chm%fupdateHOBr', 0, RC )    
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Chm%fupdateHOBr = 0.0_fp
       CALL Register_ChmField( am_I_Root, chmID, State_Chm%fupdateHOBr, &
                               State_Chm, RC )
       IF ( RC /= GC_SUCCESS ) RETURN

    ENDIF

    !=======================================================================
    ! Print out the list of registered fields
    !=======================================================================
    IF ( am_I_Root ) THEN
       WRITE( 6, 10 )
10     FORMAT( /, 'Registered variables contained within the State_Chm object:' )
       WRITE( 6, '(a)' ) REPEAT( '=', 79 )
    ENDIF
    CALL Registry_Print( am_I_Root   = am_I_Root,             &
                         Registry    = State_Chm%Registry,    &
                         ShortFormat = .TRUE.,                &
                         RC          = RC                    )

    !=======================================================================
    ! Special handling for the Hg and tagHg  simulations: get the # of Hg
    ! categories for total & tagged tracers from the species database
    !=======================================================================
    IF ( Input_Opt%ITS_A_MERCURY_SIM ) THEN

       ! Hg0, Hg2, HgP should all have the same number of categories as
       ! returned from the species database.  If not, there's an error.
       IF ( N_Hg0_CATS == N_Hg2_CATS .and. N_Hg0_CATS == N_HgP_CATS ) THEN
          State_Chm%N_Hg_CATS = N_Hg0_CATS
       ELSE
          RC = GC_FAILURE
          PRINT*, '### Inconsistent number of Hg categories!'
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
       
    ENDIF
   
    !=======================================================================
    ! Cleanup and quit
    !=======================================================================

    ! Trap error
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Error encountered!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    ! Free pointer for safety's sake
    ThisSpc => NULL()

    ! Echo output
    IF ( am_I_Root ) THEN
       print*, REPEAT( '#', 79 )
    ENDIF 

    ! Format statement
100 FORMAT( I3, 1x, A10, 3x, F7.2 )
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
  SUBROUTINE Cleanup_State_Chm( am_I_Root, State_Chm, RC )
!
! !USES:
!
    USE Species_Database_Mod, ONLY : Cleanup_Species_Database
!
! !INPUT PARAMETERS:
! 
    LOGICAL,        INTENT(IN)    :: am_I_Root    ! Is this the root CPU?
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
!  For now the am_I_Root and RC arguments are not used.  We include these
!  for consistency and also to facilitate future expansion. (bmy, 10/16/12)
!
! !REVISION HISTORY: 
!  15 Oct 2012 - R. Yantosca - Initial version
!  26 Oct 2012 - R. Yantosca - Now deallocate Strat_P, Strat_k fields
!  26 Feb 2013 - M. Long     - Now deallocate State_Chm%DEPSAV
!  11 Dec 2014 - R. Yantosca - Remove TRAC_TEND and DEPSAV fields
!  28 Aug 2015 - R. Yantosca - Remove stratospheric chemistry fields; 
!                              these are all now read in via HEMCO
!  28 Aug 2015 - R. Yantosca - Also initialize the species database object
!  29 Jun 2017 - R. Yantosca - Add error checks for deallocations.  Also
!                              destroy the registry of State_Chm fields.
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
    ! Deallocate fields
    !=======================================================================
    IF ( ASSOCIATED( State_Chm%Map_Advect) ) THEN
       DEALLOCATE( State_Chm%Map_Advect, STAT=RC )
       CALL GC_CheckVar( 'State_Chm%Map_Advect', 3, RC )
       RETURN
    ENDIF

    IF ( ASSOCIATED( State_Chm%Map_DryDep ) ) THEN
       DEALLOCATE( State_Chm%Map_DryDep, STAT=RC )
       CALL GC_CheckVar( 'State_Chm%Map_Drydep', 3, RC )
       RETURN
    ENDIF

    IF ( ASSOCIATED( State_Chm%Map_KppSpc ) ) THEN
       DEALLOCATE( State_Chm%Map_KppSpc, STAT=RC )
       CALL GC_CheckVar( 'State_Chm%Map_KppSpc', 3, RC )
       RETURN
    ENDIF

    IF ( ASSOCIATED( State_Chm%Map_WetDep ) ) THEN
       DEALLOCATE( State_Chm%Map_WetDep, STAT=RC )
       CALL GC_CheckVar( 'State_Chm%Map_WetDep', 3, RC )
       RETURN
    ENDIF

    IF ( ASSOCIATED( State_Chm%Species ) ) THEN
       DEALLOCATE( State_Chm%Species, STAT=RC )
       CALL GC_CheckVar( 'State_Chm%Map_Species', 3, RC )
       RETURN
    ENDIF

    IF ( ASSOCIATED( State_Chm%Hg_Cat_Name ) ) THEN
       DEALLOCATE( State_Chm%Hg_Cat_Name, STAT=RC )
       CALL GC_CheckVar( 'State_Chm%Hg_Cat_Name', 3, RC )
       RETURN
    ENDIF

    IF ( ASSOCIATED( State_Chm%Hg0_Id_List ) ) THEN
       DEALLOCATE( State_Chm%Hg0_Id_List, STAT=RC ) 
       CALL GC_CheckVar( 'State_Chm%Hg0_Id_List', 3, RC )
       RETURN
    ENDIF

    IF ( ASSOCIATED( State_Chm%Hg2_Id_List ) ) THEN
       DEALLOCATE( State_Chm%Hg2_Id_List, STAT=RC )
       CALL GC_CheckVar( 'State_Chm%Hg2_Id_List', 3, RC )
       RETURN
    ENDIF

    IF ( ASSOCIATED( State_Chm%HgP_Id_List ) ) THEN
       DEALLOCATE( State_Chm%HgP_Id_List, STAT=RC )
       CALL GC_CheckVar( 'State_Chm%HgP_Id_List', 3, RC )
       RETURN
    ENDIF

    IF ( ASSOCIATED( State_Chm%AeroArea ) ) THEN
       DEALLOCATE( State_Chm%AeroArea, STAT=RC )
       CALL GC_CheckVar( 'State_Chm%AeroArea', 3, RC )
       RETURN
    ENDIF

    IF ( ASSOCIATED( State_Chm%AeroRadi ) ) THEN
       DEALLOCATE( State_Chm%AeroRadi, STAT=RC )
       CALL GC_CheckVar( 'State_Chm%AeroRadi', 3, RC )
       RETURN
    ENDIF

    IF ( ASSOCIATED( State_Chm%WetAeroArea ) ) THEN
       DEALLOCATE( State_Chm%WetAeroArea, STAT=RC )
       CALL GC_CheckVar( 'State_Chm%WetAeroArea', 3, RC )
       RETURN
    ENDIF

    IF ( ASSOCIATED( State_Chm%WetAeroRadi ) ) THEN
       DEALLOCATE( State_Chm%WetAeroRadi, STAT=RC )
       CALL GC_CheckVar( 'State_Chm%WetAeroRadi', 3, RC )
       RETURN
    ENDIF

    IF ( ASSOCIATED(State_Chm%pHCloud) ) THEN
       DEALLOCATE(State_Chm%pHCloud)
    ENDIF

    IF ( ASSOCIATED(State_Chm%SSAlk) ) THEN
       DEALLOCATE(State_Chm%SSAlk)
    ENDIF

    IF ( ASSOCIATED( State_Chm%STATE_PSC ) ) THEN
       DEALLOCATE( State_Chm%STATE_PSC, STAT=RC  )
       CALL GC_CheckVar( 'State_Chm%State_PSC', 3, RC )
       RETURN
    ENDIF
       
    IF ( ASSOCIATED( State_Chm%KHETI_SLA ) ) THEN
       DEALLOCATE( State_Chm%KHETI_SLA, STAT=RC  )
       CALL GC_CheckVar( 'State_Chm%KHETI_SLA', 3, RC )
       RETURN
    ENDIF

    IF ( ASSOCIATED( State_Chm%PH_SAV ) ) THEN
       DEALLOCATE( State_Chm%PH_SAV, STAT=RC )
       CALL GC_CheckVar( 'State_Chm%PH_SAV', 3, RC )
       RETURN
    ENDIF

    IF ( ASSOCIATED( State_Chm%HPLUS_SAV ) ) THEN
       DEALLOCATE( State_Chm%HPLUS_SAV, STAT=RC )
       CALL GC_CheckVar( 'State_Chm%HPLUS_SAV', 3, RC )
       RETURN
    ENDIF

    IF ( ASSOCIATED( State_Chm%WATER_SAV ) ) THEN
       DEALLOCATE( State_Chm%WATER_SAV, STAT=RC )
       CALL GC_CheckVar( 'State_Chm%WATER_SAV', 3, RC )
       RETURN
    ENDIF

    IF ( ASSOCIATED( State_Chm%SULRAT_SAV ) ) THEN
       DEALLOCATE( State_Chm%SULRAT_SAV, STAT=RC )
       CALL GC_CheckVar( 'State_Chm%SULRAT_SAV', 3, RC )
       RETURN
    ENDIF

    IF ( ASSOCIATED( State_Chm%NARAT_SAV ) ) THEN
       DEALLOCATE( State_Chm%NARAT_SAV, STAT=RC )
       CALL GC_CheckVar( 'State_Chm%NARAT_SAV', 3, RC )
       RETURN
    ENDIF

    IF ( ASSOCIATED( State_Chm%ACIDPUR_SAV ) ) THEN
       DEALLOCATE( State_Chm%ACIDPUR_SAV, STAT=RC )
       CALL GC_CheckVar( 'State_Chm%ACIDPUR_SAV', 3, RC )
       RETURN
    ENDIF

    IF ( ASSOCIATED( State_Chm%BISUL_SAV ) ) THEN
       DEALLOCATE( State_Chm%BISUL_SAV, STAT=RC )
       CALL GC_CheckVar( 'State_Chm%BISUL_SAV', 3, RC )
       RETURN
    ENDIF

    IF ( ASSOCIATED( State_Chm%HSO3_AQ ) ) THEN
       DEALLOCATE( State_Chm%HSO3_AQ, STAT=RC )
       CALL GC_CheckVar( 'State_Chm%HSO3_AQ', 3, RC )
       RETURN
    ENDIF

    IF ( ASSOCIATED( State_Chm%SO3_AQ ) ) THEN
       DEALLOCATE( State_Chm%SO3_AQ, STAT=RC )
       CALL GC_CheckVar( 'State_Chm%SO3_AQ', 3, RC )
       RETURN
    ENDIF

    IF ( ASSOCIATED( State_Chm%fupdateHOBr ) ) THEN
       DEALLOCATE( State_Chm%fupdateHOBr, STAT=RC )
       CALL GC_CheckVar( 'State_Chm%fupdateHOBr', 3, RC )
       RETURN
    ENDIF

    !=======================================================================
    ! Deallocate the species database object field
    !=======================================================================
    CALL Cleanup_Species_Database( am_I_Root, State_Chm%SpcData, RC )

    !=======================================================================
    ! Destroy the registry of fields for this module
    !=======================================================================
    CALL Registry_Destroy( am_I_Root, State_Chm%Registry, RC )
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Could not destroy registry object State_Chm%Registry!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

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
!  information about each State_Chm field.
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
    USE Registry_Params_Mod
!
! !INPUT PARAMETERS:
! 
    LOGICAL,             INTENT(IN)  :: am_I_Root   ! Is this the root CPU?
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

       CASE ( 'SPC' )
          IF ( isDesc  ) Desc  = 'Concentration for species'
          IF ( isUnits ) Units = 'varies'
          IF ( isRank  ) Rank  = 3
          IF ( isSpecies ) PerSpecies = 'ALL'

       CASE ( 'AEROAREA_MDUST1' )
          IF ( isDesc  ) Desc  = 'Dry aerosol area for mineral dust (0.15 um)'
          IF ( isUnits ) Units = 'cm2 cm-3'
          IF ( isRank  ) Rank  = 3

       CASE ( 'AEROAREA_MDUST2' )
          IF ( isDesc  ) Desc  = 'Dry aerosol area for mineral dust (0.25 um)'
          IF ( isUnits ) Units = 'cm2 cm-3'
          IF ( isRank  ) Rank  = 3

       CASE ( 'AEROAREA_MDUST3' )
          IF ( isDesc  ) Desc  = 'Dry aerosol area for mineral dust (0.4 um)'
          IF ( isUnits ) Units = 'cm2 cm-3'
          IF ( isRank  ) Rank  = 3

       CASE ( 'AEROAREA_MDUST4' )
          IF ( isDesc  ) Desc  = 'Dry aerosol area for mineral dust (0.8 um)'
          IF ( isUnits ) Units = 'cm2 cm-3'
          IF ( isRank  ) Rank  = 3

       CASE ( 'AEROAREA_MDUST5' )
          IF ( isDesc  ) Desc  = 'Dry aerosol area for mineral dust (1.5 um)'
          IF ( isUnits ) Units = 'cm2 cm-3'
          IF ( isRank  ) Rank  = 3

       CASE ( 'AEROAREA_MDUST6' )
          IF ( isDesc  ) Desc  = 'Dry aerosol area for mineral dust (2.5 um)'
          IF ( isUnits ) Units = 'cm2 cm-3'
          IF ( isRank  ) Rank  = 3

       CASE ( 'AEROAREA_MDUST7' )
          IF ( isDesc  ) Desc  = 'Dry aerosol area for mineral dust (4.0 um)'
          IF ( isUnits ) Units = 'cm2 cm-3'
          IF ( isRank  ) Rank  = 3

       CASE ( 'AEROAREA_SULF' )
          IF ( isDesc  ) Desc  = 'Dry aerosol area for black carbon'
          IF ( isUnits ) Units = 'cm2 cm-3'
          IF ( isRank  ) Rank  = 3

       CASE ( 'AEROAREA_BC' )
          IF ( isDesc  ) Desc  = 'Dry aerosol area for black carbon'
          IF ( isUnits ) Units = 'cm2 cm-3'
          IF ( isRank  ) Rank  = 3

       CASE ( 'AEROAREA_OC' )
          IF ( isDesc  ) Desc  = 'Dry aerosol area for organic carbon'
          IF ( isUnits ) Units = 'cm2 cm-3'
          IF ( isRank  ) Rank  = 3

       CASE ( 'AEROAREA_SSA' )
          IF ( isDesc  ) Desc  = 'Dry aerosol area for sea salt,' &
                                 // ' accumulation mode'
          IF ( isUnits ) Units = 'cm2 cm-3'
          IF ( isRank  ) Rank  = 3

       CASE ( 'AEROAREA_SSC' )
          IF ( isDesc  ) Desc  = 'Dry aerosol area for sea salt, coarse mode'
          IF ( isUnits ) Units = 'cm2 cm-3'
          IF ( isRank  ) Rank  = 3

       CASE ( 'AEROAREA_BGSULF' )
          IF ( isDesc  ) Desc  = 'Dry aerosol area for background' &
                                 // ' stratospheric sulfate'
          IF ( isUnits ) Units = 'cm2 cm-3'
          IF ( isRank  ) Rank  = 3

       CASE ( 'AEROAREA_ICEI' )
          IF ( isDesc  ) Desc  = 'Dry aerosol area for irregular ice cloud' &
                                 // ' (Mischenko)'
          IF ( isUnits ) Units = 'cm2 cm-3'
          IF ( isRank  ) Rank  = 3

       CASE ( 'AERORADI_MDUST1' )
          IF ( isDesc  ) Desc  = 'Dry aerosol radius for mineral dust (0.15 um)'
          IF ( isUnits ) Units = 'cm'
          IF ( isRank  ) Rank  = 3

       CASE ( 'AERORADI_MDUST2' )
          IF ( isDesc  ) Desc  = 'Dry aerosol radius for mineral dust (0.25 um)'
          IF ( isUnits ) Units = 'cm'
          IF ( isRank  ) Rank  = 3

       CASE ( 'AERORADI_MDUST3' )
          IF ( isDesc  ) Desc  = 'Dry aerosol radius for mineral dust (0.4 um)'
          IF ( isUnits ) Units = 'cm'
          IF ( isRank  ) Rank  = 3

       CASE ( 'AERORADI_MDUST4' )
          IF ( isDesc  ) Desc  = 'Dry aerosol radius for mineral dust (0.8 um)'
          IF ( isUnits ) Units = 'cm'
          IF ( isRank  ) Rank  = 3

       CASE ( 'AERORADI_MDUST5' )
          IF ( isDesc  ) Desc  = 'Dry aerosol radius for mineral dust (1.5 um)'
          IF ( isUnits ) Units = 'cm'
          IF ( isRank  ) Rank  = 3

       CASE ( 'AERORADI_MDUST6' )
          IF ( isDesc  ) Desc  = 'Dry aerosol radius for mineral dust (2.5 um)'
          IF ( isUnits ) Units = 'cm'
          IF ( isRank  ) Rank  = 3

       CASE ( 'AERORADI_MDUST7' )
          IF ( isDesc  ) Desc  = 'Dry aerosol radius for mineral dust (4.0 um)'
          IF ( isUnits ) Units = 'cm'
          IF ( isRank  ) Rank  = 3

       CASE ( 'AERORADI_SULF' )
          IF ( isDesc  ) Desc  = 'Dry aerosol radius for tropospheric sulfate'
          IF ( isUnits ) Units = 'cm'
          IF ( isRank  ) Rank  = 3

       CASE ( 'AERORADI_BC' )
          IF ( isDesc  ) Desc  = 'Dry aerosol radius for black carbon'
          IF ( isUnits ) Units = 'cm'
          IF ( isRank  ) Rank  = 3


       CASE ( 'AERORADI_OC' )
          IF ( isDesc  ) Desc  = 'Dry aerosol radius for organic carbon'
          IF ( isUnits ) Units = 'cm'
          IF ( isRank  ) Rank  = 3

       CASE ( 'AERORADI_SSA' )
          IF ( isDesc  ) Desc  = 'Dry aerosol radius for sea salt,' & 
                                 // ' accumulation mode'
          IF ( isUnits ) Units = 'cm'
          IF ( isRank  ) Rank  = 3

       CASE ( 'AERORADI_SSC' )
          IF ( isDesc  ) Desc  = 'Dry aerosol radius for sea salt, coarse mode'
          IF ( isUnits ) Units = 'cm'
          IF ( isRank  ) Rank  = 3

       CASE ( 'AERORADI_BGSULF' )
          IF ( isDesc  ) Desc  = 'Dry aerosol radius for background' &
                                 // ' stratospheric sulfate'
          IF ( isUnits ) Units = 'cm'
          IF ( isRank  ) Rank  = 3

       CASE ( 'AERORADI_ICEI' )
          IF ( isDesc  ) Desc  = 'Dry aerosol radius for irregular ice' &
                                 // ' cloud (Mischenko)'
          IF ( isUnits ) Units = 'cm'
          IF ( isRank  ) Rank  = 3

       CASE ( 'WETAEROAREA_MDUST1' )
          IF ( isDesc  ) Desc  = 'Wet aerosol area for mineral dust (0.15 um)'
          IF ( isUnits ) Units = 'cm2 cm-3'
          IF ( isRank  ) Rank  = 3

       CASE ( 'WETAEROAREA_MDUST2' )
          IF ( isDesc  ) Desc  = 'Wet aerosol area for mineral dust (0.25 um)'
          IF ( isUnits ) Units = 'cm2 cm-3'
          IF ( isRank  ) Rank  = 3

       CASE ( 'WETAEROAREA_MDUST3' )
          IF ( isDesc  ) Desc  = 'Wet aerosol area for mineral dust (0.4 um)'
          IF ( isUnits ) Units = 'cm2 cm-3'
          IF ( isRank  ) Rank  = 3

       CASE ( 'WETAEROAREA_MDUST4' )
          IF ( isDesc  ) Desc  = 'Wet aerosol area for mineral dust (0.8 um)'
          IF ( isUnits ) Units = 'cm2 cm-3'
          IF ( isRank  ) Rank  = 3

       CASE ( 'WETAEROAREA_MDUST5' )
          IF ( isDesc  ) Desc  = 'Wet aerosol area for mineral dust (1.5 um)'
          IF ( isUnits ) Units = 'cm2 cm-3'
          IF ( isRank  ) Rank  = 3

       CASE ( 'WETAEROAREA_MDUST6' )
          IF ( isDesc  ) Desc  = 'Wet aerosol area for mineral dust (2.5 um)'
          IF ( isUnits ) Units = 'cm2 cm-3'
          IF ( isRank  ) Rank  = 3

       CASE ( 'WETAEROAREA_MDUST7' )
          IF ( isDesc  ) Desc  = 'Wet aerosol area for mineral dust (4.0 um)'
          IF ( isUnits ) Units = 'cm2 cm-3'
          IF ( isRank  ) Rank  = 3


       CASE ( 'WETAEROAREA_SULF' )
          IF ( isDesc  ) Desc  = 'Wet aerosol area for tropospheric sulfate'
          IF ( isUnits ) Units = 'cm2 cm-3'
          IF ( isRank  ) Rank  = 3

       CASE ( 'WETAEROAREA_BC' )
          IF ( isDesc  ) Desc  = 'Wet aerosol area for black carbon'
          IF ( isUnits ) Units = 'cm2 cm-3'
          IF ( isRank  ) Rank  = 3

       CASE ( 'WETAEROAREA_OC' )
          IF ( isDesc  ) Desc  = 'Wet aerosol area for organic carbon'
          IF ( isUnits ) Units = 'cm2 cm-3'
          IF ( isRank  ) Rank  = 3

       CASE ( 'WETAEROAREA_SSA' )
          IF ( isDesc  ) Desc  = 'Wet aerosol area for sea salt,' &
                                 // ' accumulation mode'
          IF ( isUnits ) Units = 'cm2 cm-3'
          IF ( isRank  ) Rank  = 3

       CASE ( 'WETAEROAREA_SSC' )
          IF ( isDesc  ) Desc  = 'Wet aerosol area for sea salt, coarse mode'
          IF ( isUnits ) Units = 'cm2 cm-3'
          IF ( isRank  ) Rank  = 3

       CASE ( 'WETAEROAREA_BGSULF' )
          IF ( isDesc  ) Desc  = 'Wet aerosol area for background' &
                                 // ' stratospheric sulfate'
          IF ( isUnits ) Units = 'cm2 cm-3'
          IF ( isRank  ) Rank  = 3

       CASE ( 'WETAEROAREA_ICEI' )
          IF ( isDesc  ) Desc  = 'Wet aerosol area for irregular ice cloud' &
                                 // ' (Mischenko)'
          IF ( isUnits ) Units = 'cm2 cm-3'
          IF ( isRank  ) Rank  = 3


       CASE ( 'WETAERORADI_MDUST1' )
          IF ( isDesc  ) Desc  = 'Wet aerosol radius for mineral dust (0.15 um)'
          IF ( isUnits ) Units = 'cm'
          IF ( isRank  ) Rank  = 3

       CASE ( 'WETAERORADI_MDUST2' )
          IF ( isDesc  ) Desc  = 'Wet aerosol radius for mineral dust (0.25 um)'
          IF ( isUnits ) Units = 'cm'
          IF ( isRank  ) Rank  = 3

       CASE ( 'WETAERORADI_MDUST3' )
          IF ( isDesc  ) Desc  = 'Wet aerosol radius for mineral dust (0.4 um)'
          IF ( isUnits ) Units = 'cm'
          IF ( isRank  ) Rank  = 3

       CASE ( 'WETAERORADI_MDUST4' )
          IF ( isDesc  ) Desc  = 'Wet aerosol radius for mineral dust (0.8 um)'
          IF ( isUnits ) Units = 'cm'
          IF ( isRank  ) Rank  = 3

       CASE ( 'WETAERORADI_MDUST5' )
          IF ( isDesc  ) Desc  = 'Wet aerosol radius for mineral dust (1.5 um)'
          IF ( isUnits ) Units = 'cm'
          IF ( isRank  ) Rank  = 3

       CASE ( 'WETAERORADI_MDUST6' )
          IF ( isDesc  ) Desc  = 'Wet aerosol radius for mineral dust (2.5 um)'
          IF ( isUnits ) Units = 'cm'
          IF ( isRank  ) Rank  = 3

       CASE ( 'WETAERORADI_MDUST7' )
          IF ( isDesc  ) Desc  = 'Wet aerosol radius for mineral dust (4.0 um)'
          IF ( isUnits ) Units = 'cm'
          IF ( isRank  ) Rank  = 3

       CASE ( 'WETAERORADI_SULF' )
          IF ( isDesc  ) Desc  = 'Wet aerosol radius for tropospheric sulfate'
          IF ( isUnits ) Units = 'cm'
          IF ( isRank  ) Rank  = 3

       CASE ( 'WETAERORADI_BC' )
          IF ( isDesc  ) Desc  = 'Wet aerosol radius for black carbon'
          IF ( isUnits ) Units = 'cm'
          IF ( isRank  ) Rank  = 3

       CASE ( 'WETAERORADI_OC' )
          IF ( isDesc  ) Desc  = 'Wet aerosol radius for organic carbon'
          IF ( isUnits ) Units = 'cm'
          IF ( isRank  ) Rank  = 3

       CASE ( 'WETAERORADI_SSA' )
          IF ( isDesc  ) Desc= 'Wet aerosol radius for sea salt,' &
                               // ' accumulation mode'
          IF ( isUnits ) Units = 'cm'
          IF ( isRank  ) Rank  = 3

       CASE ( 'WETAERORADI_SSC' )
          IF ( isDesc  ) Desc  = 'Wet aerosol radius for sea salt, coarse mode'
          IF ( isUnits ) Units = 'cm'
          IF ( isRank  ) Rank  = 3

       CASE ( 'WETAERORADI_BGSULF' )
          IF ( isDesc  ) Desc  = 'Wet aerosol radius for background' &
                                // ' stratospheric sulfate'
          IF ( isUnits ) Units = 'cm'
          IF ( isRank  ) Rank  = 3

       CASE ( 'WETAERORADI_ICEI' )
          IF ( isDesc  ) Desc  = 'Wet aerosol radius for irregular ice cloud' &
                                // ' (Mischenko)'
          IF ( isUnits ) Units = 'cm'
          IF ( isRank  ) Rank  = 3

       CASE ( 'STATE_PSC' )
          IF ( isDesc  ) Desc  = 'Polar stratospheric cloud type (cf Kirner' &
                                // ' et al 2011, GMD)'
          IF ( isUnits ) Units = 'count'
          IF ( isRank  ) Rank  = 3

       CASE ( 'KHSLA_N2O5+H2O' )
          IF ( isDesc  ) Desc  = 'Sticking coeeficient for N2O5 + H2O reaction'
          IF ( isUnits ) Units = '1'
          IF ( isRank  ) Rank  = 3

       CASE ( 'KHSLA_N2O5+HCL' )
          IF ( isDesc  ) Desc  = 'Sticking coeeficient for N2O5 + H2O reaction'
          IF ( isUnits ) Units = '1'
          IF ( isRank  ) Rank  = 3

       CASE ( 'KHSLA_CLNO3+H2O' )
          IF ( isDesc  ) Desc  = 'Sticking coeeficient for ClNO3 + H2O reaction'
          IF ( isUnits ) Units = '1'
          IF ( isRank  ) Rank  = 3

       CASE ( 'KHSLA_CLNO3+HCL' )
          IF ( isDesc  ) Desc  = 'Sticking coeeficient for ClNO3 + HCl reaction'
          IF ( isUnits ) Units = '1'
          IF ( isRank  ) Rank  = 3

       CASE ( 'KHSLA_CLNO3+HBR' )
          IF ( isDesc  ) Desc  = 'Sticking coeeficient for ClNO3 + HBr reaction'
          IF ( isUnits ) Units = '1'
          IF ( isRank  ) Rank  = 3

       CASE ( 'KHSLA_BRNO3+H2O' )
          IF ( isDesc  ) Desc  = 'Sticking coeeficient for BrNO3 + H2O reaction'
          IF ( isUnits ) Units = '1'
          IF ( isRank  ) Rank  = 3

       CASE ( 'KHSLA_BRNO3+HCL' )
          IF ( isDesc  ) Desc  = 'Sticking coeeficient for BrNO3 + HCl reaction'
          IF ( isUnits ) Units = '1'
          IF ( isRank  ) Rank  = 3

       CASE ( 'KHSLA_HOCL+HCL' )
          IF ( isDesc  ) Desc  = 'Sticking coeeficient for HOCl + HCl reaction'
          IF ( isUnits ) Units = '1'
          IF ( isRank  ) Rank  = 3

       CASE ( 'KHSLA_HOCL+HBR' )
          IF ( isDesc  ) Desc  = 'Sticking coeeficient for HClr + HBr reaction'
          IF ( isUnits ) Units = '1'
          IF ( isRank  ) Rank  = 3

       CASE ( 'KHSLA_HOBR+HCL' )
          IF ( isDesc  ) Desc  = 'Sticking coeeficient for HOBr + HCl reaction'
          IF ( isUnits ) Units = '1'
          IF ( isRank  ) Rank  = 3

       CASE ( 'KHSLA_HOBR+HBR' )
          IF ( isDesc  ) Desc  = 'Sticking coeeficient for HOBr + HBr reaction'
          IF ( isUnits ) Units = '1'
          IF ( isRank  ) Rank  = 3

       ! NOTE: if the rest of these are only used for diagnostics, 
       !       consider moving to State_Diag and out of State_Chm.
       CASE ( 'PH_SAV' )
          IF ( isDesc  ) Desc  = 'ISORROPIA aerosol pH'
          IF ( isUnits ) Units = '1'
          IF ( isRank  ) Rank  = 3

       CASE ( 'HPLUS_SAV' )
          IF ( isDesc  ) Desc  = 'ISORROPIA H+ concentration'
          IF ( isUnits ) Units = 'mol L-1'
          IF ( isRank  ) Rank  = 3

       CASE ( 'WATER_SAV' )
          IF ( isDesc  ) Desc  = 'ISORROPIA aerosol water concentration'
          IF ( isUnits ) Units = 'ug m-3'
          IF ( isRank  ) Rank  = 3

       CASE ( 'SULRAT_SAV' )
          IF ( isDesc  ) Desc  = 'ISORROPIA sulfate concentration'
          IF ( isUnits ) Units = 'M'
          IF ( isRank  ) Rank  = 3

       CASE ( 'NARAT_SAV' )
          IF ( isDesc  ) Desc  = 'ISORROPIA sulfate concentration'
          IF ( isUnits ) Units = 'M'
          IF ( isRank  ) Rank  = 3

       CASE ( 'ACIDPUR_SAV' )
          IF ( isDesc  ) Desc  = 'ISORROPIA ACIDPUR'
          IF ( isUnits ) Units = 'M'
          IF ( isRank  ) Rank  = 3

       CASE ( 'BISUL_SAV' )
          IF ( isDesc  ) Desc  = 'ISORROPIA Bisulfate (general acid)' &
                                 // ' concentration'
          IF ( isUnits ) Units = 'M'
          IF ( isRank  ) Rank  =  3

       CASE ( 'HSO3_AQ' )
          IF ( isDesc  ) Desc  = 'Cloud bisulfite concentration'
          IF ( isUnits ) Units = 'mol/L'
          IF ( isRank  ) Rank  =  3

       CASE ( 'SO3_AQ' )
          IF ( isDesc  ) Desc  = 'Cloud sulfite concentration'
          IF ( isUnits ) Units = 'mol/L'
          IF ( isRank  ) Rank  =  3

       CASE ( 'fupdateHOBr' )
          IF ( isDesc  ) Desc  = 'Correction factor for HOBr removal by SO2'
          IF ( isUnits ) Units = '1'
          IF ( isRank  ) Rank  =  3


       CASE DEFAULT
          Found = .False.
          ErrMsg = 'Metadata not found for State_Chm field ' // &
                   TRIM( metadataID ) // ' when search for all caps name ' &
                   // TRIM( Name_AllCaps )
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN

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
! !DESCRIPTION:
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Register_ChmField_R4_3D( am_I_Root,  metadataID, Ptr2Data,  &
                                      State_Chm,  RC )
!
! !USES:
!
!
! !INPUT PARAMETERS:
!
    LOGICAL,           INTENT(IN)    :: am_I_Root       ! Root CPU?
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
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!   
    CHARACTER(LEN=255)     :: ErrMsg, ErrMsg_reg, ThisLoc
    CHARACTER(LEN=255)     :: desc, units, perSpecies
    CHARACTER(LEN=255)     :: thisSpcName, thisSpcDesc
    INTEGER                :: N, rank, type,  vloc
    LOGICAL                :: found
    TYPE(Species), POINTER :: SpcInfo

    ! Initialize
    RC = GC_SUCCESS
    ThisLoc = ' -> at Register_ChmField_R4_3D (in Headers/state_chm_mod.F90)'
    ErrMsg_reg = 'Error encountered while registering State_Chm field'

    CALL Get_Metadata_State_Chm( am_I_Root, metadataID,  Found,  RC,   &
                                 desc=desc, units=units, rank=rank,    &
                                 type=type, vloc=vloc, perSpecies=perSpecies )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL GC_Error( ErrMsg_reg, RC, ThisLoc )
       RETURN
    ENDIF
    
    ! If not tied to species then simply register the single field
    IF ( perSpecies == '' ) THEN
       
       ! Check that metadata consistent with data size
       IF ( rank /= 3 ) THEN
          ErrMsg = 'Data dims and metadata rank do not match for ' &
                   // TRIM(metadataID)
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

       CALL Registry_AddField( am_I_Root    = am_I_Root,            &
                               Registry     = State_Chm%Registry,   &
                               State        = State_Chm%State,      &
                               Variable     = metadataID,           &
                               Description  = desc,                 &
                               Units        = units,                &
                               Data3d_4     = Ptr2Data,             &
                               RC           = RC                   )
       IF ( RC /= GC_SUCCESS ) THEN
          CALL GC_Error( ErrMsg_reg, RC, ThisLoc )
          RETURN
       ENDIF

    ELSE
       ErrMsg = 'Handling of PerSpecies metadata ' // TRIM(perSpecies) // &
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
! !IROUTINE: Register_ChmField_Rfp_3D
!
! !DESCRIPTION:
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Register_ChmField_Rfp_3D( am_I_Root, metadataID, Ptr2Data,  &
                                       State_Chm, RC )
!
! !USES:
!
!
! !INPUT PARAMETERS:
!
    LOGICAL,           INTENT(IN)    :: am_I_Root       ! Root CPU?
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
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!   
    CHARACTER(LEN=255)     :: ErrMsg, ErrMsg_reg, ThisLoc
    CHARACTER(LEN=255)     :: desc, units, perSpecies
    CHARACTER(LEN=255)     :: thisSpcName, thisSpcDesc
    INTEGER                :: N, rank, type,  vloc
    LOGICAL                :: found
    TYPE(Species), POINTER :: SpcInfo

    ! Initialize
    RC = GC_SUCCESS
    ThisLoc = ' -> at Register_ChmField_Rfp_3D (in Headers/state_chm_mod.F90)'
    ErrMsg_reg = 'Error encountered while registering State_Chm field'

    CALL Get_Metadata_State_Chm( am_I_Root, metadataID,  Found,  RC,   &
                                 desc=desc, units=units, rank=rank,    &
                                 type=type, vloc=vloc, perSpecies=perSpecies )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL GC_Error( ErrMsg_reg, RC, ThisLoc )
       RETURN
    ENDIF
    
    ! If not tied to species then simply register the single field
    IF ( perSpecies == '' ) THEN
       
       ! Check that metadata consistent with data size
       IF ( rank /= 3 ) THEN
          ErrMsg = 'Data dims and metadata rank do not match for ' &
                   // TRIM(metadataID)
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
       CALL Registry_AddField( am_I_Root    = am_I_Root,            &
                               Registry     = State_Chm%Registry,   &
                               State        = State_Chm%State,      &
                               Variable     = metadataID,           &
                               Description  = desc,                 &
                               Units        = units,                &
                               Data3d       = Ptr2Data,             &
                               RC           = RC                   )
       IF ( RC /= GC_SUCCESS ) THEN
          CALL GC_Error( ErrMsg_reg, RC, ThisLoc )
          RETURN
       ENDIF

    ELSE
       ErrMsg = 'Handling of PerSpecies metadata ' // TRIM(perSpecies) // &
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
! !DESCRIPTION:
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Register_ChmField_Rfp_4D( am_I_Root,  metadataID, Ptr2Data,  &
                                       State_Chm,  RC,         Ncat )
!
! !USES:
!
!
! !INPUT PARAMETERS:
!
    LOGICAL,           INTENT(IN)    :: am_I_Root         ! Root CPU?
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
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!   
    CHARACTER(LEN=255)     :: ErrMsg, ErrMsg_reg, ThisLoc
    CHARACTER(LEN=255)     :: desc, units, perSpecies
    CHARACTER(LEN=255)     :: thisSpcName, thisSpcDesc
    INTEGER                :: N, rank, type,  vloc
    LOGICAL                :: found
    TYPE(Species), POINTER :: SpcInfo

    ! Initialize
    RC = GC_SUCCESS
    ThisLoc = ' -> at Register_ChmField_Rfp_4D (in Headers/state_chm_mod.F90)'
    ErrMsg_reg = 'Error encountered while registering State_Chm field'

    CALL Get_Metadata_State_Chm( am_I_Root, metadataID,  Found,  RC,   &
                                 desc=desc, units=units, rank=rank,    &
                                 type=type, vloc=vloc, perSpecies=perSpecies )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL GC_Error( ErrMsg_reg, RC, ThisLoc )
       RETURN
    ENDIF

    ! Check that metadata consistent with data size
    IF ( rank /= 3 ) THEN
       ErrMsg = 'Data dims and metadata rank do not match for ' &
                // TRIM(metadataID)
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF
    
    ! If tied to all species then register each one
    IF ( perSpecies == 'ALL' ) THEN       
       DO N = 1, State_Chm%nSpecies
          SpcInfo  => State_Chm%SpcData(N)%Info
          thisSpcName = TRIM( metadataID ) // '__' // TRIM( SpcInfo%Name )
          thisSpcDesc = TRIM( Desc ) // ' ' // TRIM( SpcInfo%Name )
          CALL Registry_AddField( am_I_Root    = am_I_Root,            &
                                  Registry     = State_Chm%Registry ,  &
                                  State        = State_Chm%State,      &
                                  Variable     = thisSpcName,          &
                                  Description  = thisSpcDesc,          &
                                  Units        = units,                &
                                  Data3d       = Ptr2Data(:,:,:,N),    &
                                  RC           = RC                   )
          SpcInfo => NULL()
          IF ( RC /= GC_SUCCESS ) THEN
             CALL GC_Error( ErrMsg_reg, RC, ThisLoc )
             RETURN
          ENDIF
       ENDDO
    ! If tied to a given category, only registry that one
    ELSE IF ( PRESENT(Ncat) ) THEN
       CALL Registry_AddField( am_I_Root    = am_I_Root,            &
                               Registry     = State_Chm%Registry ,  &
                               State        = State_Chm%State,      &
                               Variable     = metadataID ,          &
                               Description  = desc,                 &
                               Units        = units,                &
                               Data3d       = Ptr2Data(:,:,:,Ncat), &
                               RC           = RC                   )
    ELSE
       ErrMsg = 'Handling of PerSpecies metadata ' // TRIM(perSpecies) // &
                ' is not implemented for this combo of data type and size'
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
    USE CHARPAK_MOD, ONLY : TRANUC
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
! !REMARKS
!
! !REVISION HISTORY: 
!  07 Oct 2016 - M. Long     - Initial version
!  15 Jun 2016 - M. Sulprizio- Make species name uppercase before computing hash
!  17 Aug 2016 - M. Sulprizio- Tracer flag 'T' is now advected species flag 'A'
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER           :: N, Hash
    CHARACTER(LEN=14) :: Name14

    !=====================================================================
    ! Ind_ begins here!
    !=====================================================================

    ! Initialize the output value
    Indx   = -1

    ! Make species name uppercase for hash algorithm
    Name14 = Name
    CALL TRANUC( Name14 )

    ! Compute the hash corresponding to the given species name
    Hash   = Str2Hash( Name14 )

    ! Loop over all entries in the Species Database object
    DO N = 1, SIZE( SpcDataLocal )

       ! Compare the hash we just created against the list of
       ! species name hashes stored in the species database
       IF( Hash == SpcDataLocal(N)%Info%NameHash  ) THEN

          IF (.not. PRESENT(Flag)) THEN

             ! Default to Species/ModelID
             Indx = SpcDataLocal(N)%Info%ModelID
             RETURN

          ELSE

             ! Only need first character of the flag for this.
             IF     (flag(1:1) .eq. 'S' .or. flag(1:1) .eq. 's') THEN

                ! Species/ModelID
                Indx = SpcDataLocal(N)%Info%ModelID
                RETURN

             ELSEIF (flag(1:1) .eq. 'A' .or. flag(1:1) .eq. 'a') THEN

                ! Advected species flag
                Indx = SpcDataLocal(N)%Info%AdvectID
                RETURN

             ELSEIF (flag(1:1) .eq. 'K' .or. flag(1:1) .eq. 'k') THEN

                ! KPP species ID
                Indx = SpcDataLocal(N)%Info%KppSpcId
                RETURN

             ELSEIF (flag(1:1) .eq. 'V' .or. flag(1:1) .eq. 'v') THEN

                ! KPP VAR ID
                Indx = SpcDataLocal(N)%Info%KppVarId
                RETURN

             ELSEIF (flag(1:1) .eq. 'F' .or. flag(1:1) .eq. 'f') THEN

                ! KPP FIX ID
                Indx = SpcDataLocal(N)%Info%KppFixId
                RETURN

             ELSEIF (flag(1:1) .eq. 'W' .or. flag(1:1) .eq. 'w') THEN

                ! WetDep ID
                Indx = SpcDataLocal(N)%Info%WetDepId
                RETURN

             ELSEIF (flag(1:1) .eq. 'D' .or. flag(1:1) .eq. 'd') THEN

                ! DryDep ID
                Indx = SpcDataLocal(N)%Info%DryDepId
                RETURN

             ENDIF

          ENDIF
          EXIT

       ENDIF

    ENDDO

    RETURN

  END FUNCTION Ind_
!EOC
END MODULE State_Chm_Mod
