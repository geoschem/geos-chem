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
#if defined( MODEL_GCHPCTM)
  USE ESMF
#endif
  USE Precision_Mod

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC :: SpcData_Init
  PUBLIC :: SpcData_Cleanup
  PUBLIC :: Spc_Print
!
! !PUBLIC TYPES:
!
  !=========================================================================
  ! Type for species index counters
  !=========================================================================
  TYPE, PUBLIC :: SpcIndCt
     INTEGER :: nAdvect  ! # of advected species
     INTEGER :: nAeroSpc ! # of aerosol species
     INTEGER :: nDryAlt  ! # of dry-dep species to save @ user-defined altitude
     INTEGER :: nDryDep  ! # of dry-deposited species
     INTEGER :: nGasSpc  ! # of gas-phase species
     INTEGER :: nHygGrth ! # of hygroscopic growth spc
     INTEGER :: nKppVar  ! # of variable KPP species
     INTEGER :: nKppFix  ! # of fixed KPP species
     INTEGER :: nKppSpc  ! # of species in KPP matrix
     INTEGER :: nPhotol  ! # of photolysis species
     INTEGER :: nOmitted ! # of omitted species
     INTEGER :: nRadNucl ! # of radionuclide species
     INTEGER :: nRealSpc ! # of total species (w/ omitted species removed)
     INTEGER :: nWetDep  ! # of wet-deposited species
     INTEGER :: nHg0     ! # of Hg0 tracers
     INTEGER :: nHg2     ! # of Hg2 tracers
     INTEGER :: nHgP     ! # of HgP tracers
  END TYPE SpcIndCt

  !=========================================================================
  ! Type for the Species Database object (vector of type Species)
  !=========================================================================
  TYPE, PUBLIC :: SpcPtr
     TYPE(Species), POINTER :: Info         ! Single entry of Species Database
  END TYPE SpcPtr

  !=========================================================================
  ! Type for single species concentrations
  !=========================================================================
  TYPE, PUBLIC :: SpcConc
#if defined( MODEL_GCHPCTM )
     REAL(ESMF_KIND_R8), POINTER :: Conc(:,:,:)
#else
     REAL(fp), POINTER :: Conc(:,:,:)
#endif
  END TYPE SpcConc

  !=========================================================================
  ! Type for individual species information
  ! (i.e. this is a single entry in the Species Database)
  !=========================================================================
  TYPE, PUBLIC :: Species

     ! Indices
     INTEGER            :: ModelId          ! Model species Id
     INTEGER            :: AdvectId         ! Advection index
     INTEGER            :: AerosolId        ! Aerosol species index
     INTEGER            :: DryAltId         ! Dry dep species at altitude Id
     INTEGER            :: DryDepId         ! Dry deposition index
     INTEGER            :: GasSpcId         ! Gas-phase species index
     INTEGER            :: HygGrthId        ! Hygroscopic growth species index
     INTEGER            :: KppVarId         ! KPP variable species index
     INTEGER            :: KppFixId         ! KPP fixed spcecies index
     INTEGER            :: KppSpcId         ! KPP species index
     INTEGER            :: OmittedId        ! Omitted species index
     INTEGER            :: PhotolId         ! Photolysis index
     INTEGER            :: RadNuclId        ! Radionuclide index
     INTEGER            :: WetDepId         ! Wet deposition index

     ! Names
     CHARACTER(LEN=31)  :: Name             ! Short name
     CHARACTER(LEN=80)  :: FullName         ! Long name
     CHARACTER(LEN=80)  :: Formula          ! Chemical formula

     ! Logical switches
     LOGICAL            :: Is_Advected      ! Is it advected?
     LOGICAL            :: Is_Aerosol       ! Is it an aerosol species?
     LOGICAL            :: Is_DryAlt        ! Is it a dry-dep species that we
                                            !  want to save at a given altitude?
     LOGICAL            :: Is_DryDep        ! Is it dry-deposited?
     LOGICAL            :: Is_Gas           ! Is it a gas?  If not, aerosol.
     LOGICAL            :: Is_HygroGrowth   ! Does it have hygroscropic growth?
     LOGICAL            :: Is_ActiveChem    ! Is it an active chemical species?
     LOGICAL            :: Is_FixedChem     ! Is it a fixed chemical species?
     LOGICAL            :: Is_Kpp           ! Is it in the KPP mechanism?
     LOGICAL            :: Is_Omitted       ! Is it omitted from the database?
     LOGICAL            :: Is_Photolysis    ! Is it an photolysis species?
     LOGICAL            :: Is_RadioNuclide  ! Is it a radionuclide species?
     LOGICAL            :: Is_WetDep        ! Is it wet-deposited?
     LOGICAL            :: Is_InRestart     ! Is it in the restart file?

     ! Molecular weights
     REAL(fp)           :: MW_g             ! Species molecular weight [g/mol]

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
     REAL(fp)           :: DD_Hstar         ! HSTAR value in drydep_mod [M/atm]

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
     LOGICAL            :: Is_Hg0           ! Is a Hg0 species?
     LOGICAL            :: Is_Hg2           ! Is a Hg2 species?
     LOGICAL            :: Is_HgP           ! Is a HgP species?

  END TYPE Species
!
! !DEFINED PARAMETERS:
!
  !=========================================================================
  ! Missing species concentration value if not in restart file and special
  ! background value not defined
  !=========================================================================
  REAL(fp), PARAMETER, PUBLIC :: MISSING_VV  = 1.0e-20_fp ! Missing spc conc
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
  SUBROUTINE SpcData_Init( Input_Opt, nSpecies, SpcData, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE Input_Opt_Mod, ONLY : OptInput
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput),       INTENT(IN)    :: Input_Opt    ! Input Options object
    INTEGER,              INTENT(IN)    :: nSpecies     ! # of species
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(SpcPtr),         POINTER       :: SpcData(:)   ! Species database
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
    ! Scalars
    INTEGER            :: N

    ! strings
    CHARACTER(LEN=255) :: varId

    !=====================================================================
    ! SpcData_Init begins here!
    !=====================================================================

    ! Assume success
    RC = GC_SUCCESS

    ! Check if already allocated
    IF ( ASSOCIATED( SpcData ) ) THEN
       CALL SpcData_Cleanup( SpcData )
    ENDIF

    ! Allocate the species database object
    varId = "State_Chm%SpcData"
    ALLOCATE( SpcData( nSpecies ), STAT=RC )
    CALL GC_CheckVar( varId, 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    ! Initialize each entry in the species database object
    DO N = 1, nSpecies

       ! Allocate
       WRITE( varId, 100 ) N
 100   FORMAT( 'State_Chm%SpcData(', i6, ')%Info' )
       ALLOCATE( SpcData(N)%Info, STAT=RC )
       CALL GC_CheckVar( varId, 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN

       ! Set all fields to missing values
       CALL Spc_Zero( SpcData(N)%Info )
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
  SUBROUTINE SpcData_Cleanup( SpcData )
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(SpcPtr), POINTER :: SpcData(:)  ! Species database object
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
    IF ( ASSOCIATED( SpcData ) ) THEN

       ! First get the size of the SpecDb object
       nSpecies = SIZE( SpcData )

       ! If there are more than 0 elements ...
       IF ( nSpecies > 0 ) THEN

          ! Nullify each entry in the species database
          DO N = 1, nSpecies
             IF( ASSOCIATED( SpcData(N)%Info ) ) THEN
                DEALLOCATE( SpcData(N)%Info )
             ENDIF
          ENDDO

          ! And free the object's memory
          DEALLOCATE( SpcData )
       ENDIF
    ENDIF

  END SUBROUTINE SpcData_Cleanup
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Spc_Zero
!
! !DESCRIPTION: Sets all fields of an object of type Species
!  to missing values.  Called at initialization.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Spc_Zero( Spc )
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(Species), INTENT(INOUT) :: Spc
!
! !REVISION HISTORY:
!  28 Apr 2020 - R. Yantosca - Initial version
!  See the subsequent Git history with the gitk browser!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Boolean/Logical
    Spc%DD_AeroDryDep   = MISSING_BOOL
    Spc%DD_DustDryDep   = MISSING_BOOL
    Spc%Is_ActiveChem   = MISSING_BOOL
    Spc%Is_Advected     = MISSING_BOOL
    Spc%Is_Aerosol      = MISSING_BOOL
    Spc%Is_DryAlt       = MISSING_BOOL
    Spc%Is_DryDep       = MISSING_BOOL
    Spc%Is_FixedChem    = MISSING_BOOL
    Spc%Is_Gas          = MISSING_BOOL
    Spc%Is_Hg0          = MISSING_BOOL
    Spc%Is_Hg2          = MISSING_BOOL
    Spc%Is_HgP          = MISSING_BOOL
    Spc%Is_HygroGrowth  = MISSING_BOOL
    Spc%Is_InRestart    = MISSING_BOOL
    Spc%Is_Kpp          = MISSING_BOOL
    Spc%Is_Omitted      = MISSING_BOOL
    Spc%Is_Photolysis   = MISSING_BOOL
    Spc%Is_RadioNuclide = MISSING_BOOL
    Spc%Is_WetDep       = MISSING_BOOL
    Spc%MP_SizeResAer   = MISSING_BOOL
    Spc%MP_SizeResNum   = MISSING_BOOL
    Spc%WD_CoarseAer    = MISSING_BOOL
    Spc%WD_Is_H2SO4     = MISSING_BOOL
    Spc%WD_Is_HNO3      = MISSING_BOOL
    Spc%WD_Is_SO2       = MISSING_BOOL
    Spc%WD_LiqAndGas    = MISSING_BOOL

    ! Integers
    Spc%AdvectId        = MISSING_INT
    Spc%AerosolId       = MISSING_INT
    Spc%DryAltId        = MISSING_INT
    Spc%DryDepId        = MISSING_INT
    Spc%GasSpcId        = MISSING_INT
    Spc%HygGrthId       = MISSING_INT
    Spc%KppFixId        = MISSING_INT
    Spc%KppSpcId        = MISSING_INT
    Spc%KppVarId        = MISSING_INT
    Spc%ModelId         = MISSING_INT
    Spc%OmittedId       = MISSING_INT
    Spc%PhotolId        = MISSING_INT
    Spc%RadNuclId       = MISSING_INT
    Spc%WetDepId        = MISSING_INT

    ! Reals (floating precision)
    Spc%BackgroundVV    = MISSING
    Spc%DD_DvzAerSnow   = MISSING
    Spc%DD_DvzMinVal    = MISSING
    Spc%DD_F0           = MISSING
    Spc%DD_KOA          = MISSING
    Spc%DD_Hstar        = MISSING
    Spc%Density         = MISSING
    Spc%MW_g            = MISSING
    Spc%Radius          = MISSING
    Spc%WD_AerScavEff   = MISSING
    Spc%WD_ConvFacI2G   = MISSING
    Spc%WD_KcScaleFac   = MISSING
    Spc%WD_RainoutEff   = MISSING
    Spc%WD_RetFactor    = MISSING

    ! Reals (8-byte precision)
    Spc%Henry_CR        = MISSING_DBLE
    Spc%Henry_K0        = MISSING_DBLE
    Spc%Henry_PKA       = MISSING_DBLE

    ! Strings
    Spc%Formula         = MISSING_STR
    Spc%FullName        = MISSING_STR
    Spc%Name            = MISSING_STR

   END SUBROUTINE Spc_Zero
!BOC
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

    !========================================================================
    ! Spc_Create begins here!
    !========================================================================
    IF ( Input_Opt%amIRoot .and. ( .not. ThisSpc%Is_Omitted ) ) THEN

       !---------------------------------------------------------------------
       ! Print general species info
       !---------------------------------------------------------------------
       WRITE( 6, "(a)" )     REPEAT( "=", 79 )
       WRITE( 6, 100 )       "ModelId        ",  ThisSpc%ModelId
       WRITE( 6, 110 )       "Name           ",  TRIM( ThisSpc%Name     )
       WRITE( 6, 110 )       "FullName       ",  TRIM( ThisSpc%FullName )
       WRITE( 6, 110 )       "Formula        ",  TRIM( ThisSpc%Formula  )
       WRITE( 6, 121 )       "MW_g           ",  ThisSpc%MW_g
       IF ( ThisSpc%Is_Gas ) THEN
          WRITE( 6, "(a)" )  "Gas or aerosol  : GAS"
       ELSE IF ( ThisSpc%Is_Aerosol ) THEN
          WRITE( 6, "(a)" )  "Gas or aerosol  : AEROSOL"
       ENDIF
       IF ( ThisSpc%Is_RadioNuclide ) THEN
          WRITE( 6, "(a)" )  "Radionuclide?   : YES"
       ENDIF

       !--------------------------------------------------------------------
       ! Print Henry"s Law info (only applicable to gas-phase species)
       !--------------------------------------------------------------------
       IF ( ThisSpc%Is_Gas ) THEN
          IF ( ThisSpc%Henry_K0 > ZERO_DBLE ) THEN
             WRITE( 6, 120 ) "Henry_K0       ", ThisSpc%Henry_K0
          ENDIF

          IF ( ThisSpc%Henry_CR > ZERO_DBLE ) THEN
             WRITE( 6, 120 ) "Henry_CR       ", ThisSpc%Henry_CR
          ENDIF
       ENDIF

       !--------------------------------------------------------------------
       ! Print aerosol-specific properties
       !--------------------------------------------------------------------
       IF ( ThisSpc%Is_Aerosol ) THEN
          IF ( ThisSpc%Density > ZERO ) THEN
             WRITE( 6, 121 ) "Density        ", ThisSpc%Density
          ENDIF

          IF ( ThisSpc%Radius > ZERO ) THEN
             WRITE( 6, 120 ) "Radius         ", ThisSpc%Radius
          ENDIF

          IF ( ThisSpc%Is_HygroGrowth ) THEN
             WRITE( 6, 130 ) "Is_HygroGrowth ", ThisSpc%Is_HygroGrowth
             WRITE( 6, 100 ) "HygGrthId      ", ThisSpc%HygGrthId
          ENDIF

          ! Microphysics properties
          IF ( ThisSpc%MP_SizeResAer ) THEN
             WRITE( 6, 130 ) "MP_SizeResAer  ", ThisSpc%MP_SizeResAer
          ENDIF
          IF ( ThisSpc%MP_SizeResNum ) THEN
             WRITE( 6, 130 ) "MP_SizeResNum  ", ThisSpc%MP_SizeResNum
          ENDIF
       ENDIF

       !--------------------------------------------------------------------
       ! Is the species advected?
       !--------------------------------------------------------------------
       IF ( ThisSpc%Is_Advected ) THEN
          WRITE( 6, 130 )    "Is_Advected    ", ThisSpc%Is_Advected
          WRITE( 6, 100 )    "AdvectId       ", ThisSpc%AdvectId
       ENDIF

       !--------------------------------------------------------------------
       ! Is the species in the KPP mechanism and is it photolyzed?
       !--------------------------------------------------------------------
       IF ( ThisSpc%Is_Kpp ) THEN
          WRITE( 6, 130 )    "Is_Kpp         ", ThisSpc%Is_Kpp
          WRITE( 6, 100 )    "KppSpcId       ", ThisSpc%KppSpcId

          IF ( ThisSpc%Is_ActiveChem ) THEN
             WRITE( 6, 130 ) "Is_ActiveChem  ", ThisSpc%Is_ActiveChem
             WRITE( 6, 100 ) "KppVarId       ", ThisSpc%KppVarId
          ENDIF

          IF ( ThisSpc%Is_FixedChem ) THEN
             WRITE( 6, 130 ) "Is_FixedChem   ", ThisSpc%Is_FixedChem
             WRITE( 6, 100 ) "KppFixId       ", ThisSpc%KppFixId
          ENDIF
       ENDIF

       !--------------------------------------------------------------------
       ! Is the species photolyzed
       !--------------------------------------------------------------------
       IF ( ThisSpc%Is_Photolysis ) THEN
          WRITE( 6, 130 ) "Is_Photolysis  ", ThisSpc%Is_Photolysis
          WRITE( 6, 100 ) "PhotolId       ", ThisSpc%PhotolId
       ENDIF

       !--------------------------------------------------------------------
       ! Is the species dry-deposited?
       !--------------------------------------------------------------------
       IF ( ThisSpc%Is_DryDep ) THEN
          WRITE( 6, 130 ) "Is_DryDep      ", ThisSpc%Is_DryDep
          WRITE( 6, 100 ) "DryDepID       ", ThisSpc%DryDepId

          IF ( ThisSpc%DD_AeroDryDep ) THEN
             WRITE( 6, 130 ) "DD_AeroDryDep  ", ThisSpc%DD_AeroDryDep
          ENDIF

          IF ( ThisSpc%DD_DustDryDep ) THEN
             WRITE( 6, 130 ) "DD_DustDryDep  ", ThisSpc%DD_DustDryDep
          ENDIF

          IF ( ThisSpc%DD_DvzAerSnow > ZERO ) THEN
             WRITE( 6, 121 ) "DD_DvzAerSnow  ", ThisSpc%DD_DvzAerSnow
          ENDIF

          IF ( SUM( ThisSpc%DD_DvzMinVal ) > ZERO ) THEN
             WRITE( 6, 140 ) "DD_DvzMinVal   ", ThisSpc%DD_DvzMinVal
          ENDIF

          IF ( ThisSpc%DD_F0 > ZERO ) THEN
             WRITE( 6, 120 ) "DD_F0          ", ThisSpc%DD_F0
          ENDIF

          IF ( ThisSpc%DD_KOA > ZERO ) THEN
             WRITE( 6, 120 ) "DD_KOA         ", ThisSpc%DD_KOA
          ENDIF

          IF ( ThisSpc%DD_Hstar > ZERO ) THEN
             WRITE( 6, 120 ) "DD_Hstar       ", ThisSpc%DD_Hstar
          ENDIF

          IF ( ThisSpc%Is_DryAlt ) THEN
             WRITE( 6, 130 ) "Is_DryAlt      ", ThisSpc%Is_DryAlt
          ENDIF
       ENDIF

       !--------------------------------------------------------------------
       ! Is the species wet-deposited?
       !--------------------------------------------------------------------
       IF ( ThisSpc%Is_WetDep ) THEN
          WRITE( 6, 130 )    "Is_WetDep      ", ThisSpc%Is_WetDep
          WRITE( 6, 100 )    "WetDepID       ", ThisSpc%WetDepId

          IF ( ThisSpc%WD_LiqAndGas ) THEN
             WRITE( 6, 130 ) "WD_LiqAndGas   ", ThisSpc%WD_LiqAndGas
             WRITE( 6, 120 ) "WD_ConvFacI2G  ", ThisSpc%WD_ConvFacI2G
          ENDIF

          IF ( ThisSpc%WD_CoarseAer ) THEN
             WRITE( 6, 130 ) "WD_CoarseAer   ", ThisSpc%WD_CoarseAer
          ENDIF

          IF ( ThisSpc%WD_AerScavEff > ZERO ) THEN
             WRITE( 6, 120 ) "WD_AerScavEff  ", ThisSpc%WD_AerScavEff
          ENDIF

          IF ( SUM( ThisSpc%WD_KcScaleFac ) > ZERO ) THEN
             WRITE( 6, 140 ) "WD_KcScaleFac  ", ThisSpc%WD_KcScaleFac
          ENDIF

          IF ( SUM( ThisSpc%WD_RainoutEff ) > ZERO ) THEN
             WRITE( 6, 140 ) "WD_RainoutEff  ", ThisSpc%WD_RainoutEff
          ENDIF

          IF ( ThisSpc%WD_RetFactor > ZERO ) THEN
             WRITE( 6, 121 ) "WD_RetFactor   ", ThisSpc%WD_RetFactor
          ENDIF

          IF ( ThisSpc%WD_Is_H2SO4 ) THEN
             WRITE( 6, 130 ) "WD_Is_H2SO4    ", ThisSpc%WD_Is_H2SO4
          ENDIF

          IF ( ThisSpc%WD_Is_HNO3 ) THEN
             WRITE( 6, 130 ) "WD_Is_HNO3     ",  ThisSpc%WD_Is_HNO3
          ENDIF

          IF ( ThisSpc%WD_Is_SO2 ) THEN
             WRITE( 6, 130 ) "WD_Is_SO2      ",  ThisSpc%WD_Is_SO2
          ENDIF
       ENDIF

       !--------------------------------------------------------------------
       ! Is the species a mercury species?
       !--------------------------------------------------------------------
       IF ( ThisSpc%Is_Hg0 ) THEN
          WRITE( 6, 130 )    "Is_Hg0         ",  ThisSpc%Is_Hg0
       ENDIF

       IF ( ThisSpc%Is_Hg2 ) THEN
          WRITE( 6, 130 )    "Is_Hg2         ",  ThisSpc%Is_Hg2
       ENDIF

       IF ( ThisSpc%Is_HgP ) THEN
          WRITE( 6, 130 )    "Is_HgP         ",  ThisSpc%Is_HgP
       ENDIF

       !--------------------------------------------------------------------
       ! Print default background concentration
       !--------------------------------------------------------------------
       IF ( ThisSpc%BackgroundVV > ZERO ) THEN
          WRITE( 6, 120 )    "BackgroundVV   ", ThisSpc%BackgroundVV
       ENDIF

       !--------------------------------------------------------------------
       ! Format statements
       !--------------------------------------------------------------------
 100   FORMAT( a, " : ", i8          )
 110   FORMAT( a, " : ", a           )
 120   FORMAT( a, " : ", es13.6      )
 121   FORMAT( a, " : ", f8.2        )
 130   FORMAT( a, " : ", L1          )
 140   FORMAT( a, " : ", 3(f8.2, 1x) )
    ENDIF

  END SUBROUTINE Spc_Print
!EOC
END MODULE Species_Mod
