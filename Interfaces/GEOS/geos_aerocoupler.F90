#include "MAPL_Generic.h"
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Model                            !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: geos_aerocupler
!
! !DESCRIPTION: Module with routines and variables to couple GEOS-Chem aerosols 
!  with GEOS 
!\\
!\\
! !INTERFACE:
!
MODULE GEOS_AeroCoupler
!
! !USES:
!
  ! MAPL/ESMF
  USE ESMF     
  USE MAPL_Mod 
  USE Chem_Mod
  ! GEOS-Chem
  USE Precision_Mod
  USE ErrCode_Mod                                    ! Error numbers
  USE Input_Opt_Mod,         ONLY : OptInput
  USE State_Chm_Mod,         ONLY : ChmState, Ind_   ! Chemistry State obj
  USE State_Met_Mod,         ONLY : MetState         ! Meteorology State obj
  USE State_Diag_Mod,        ONLY : DgnState         ! Diagnostics State obj
  USE State_Grid_Mod,        ONLY : GrdState         ! Grid State obj

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC   :: GEOS_FillAeroBundle 
  PUBLIC   :: GEOS_AeroSetServices
  PUBLIC   :: GEOS_AeroInit
  PUBLIC   :: GEOS_AerosolOptics
!
! !PRIVATE MEMBER FUNCTIONS:
!
  PRIVATE  :: FillAeroDP_
!
! !PUBLIC TYPES
!
!
! !PRIVATE TYPES:
!
  ! Mie table
  TYPE(Chem_Mie)             :: geoschemMieTable(2)
  INTEGER, PARAMETER         :: instanceComputational = 1

  ! GEOS-5 only (also in gchp_providerservices_mod but don't use that yet):
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

! !REVISION HISTORY:
!  05 Jul 2022 - C. Keller - initial version (refactored Chem_GridCompMod)
!  See https://github.com/geoschem/geos-chem for full history
!EOP
!------------------------------------------------------------------------------
!BOC
CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Model                            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: GEOS_AeroSetServices
!
! !DESCRIPTION: Set services for aerosol coupling (if requested) 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE GEOS_AeroSetServices( GC, DoAERO, RC ) 
!
! !USES:
!
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ESMF_GridComp), INTENT(INOUT), TARGET :: GC       ! Ref to this GridComp
    LOGICAL,             INTENT(INOUT)         :: DoAERO   ! Is GCC the aerosol provider?
    INTEGER,             INTENT(INOUT)         :: RC       ! Success or failure?
!
! !REMARKS:
!
!
! !REMARKS:
!
! !REVISION HISTORY:
!  12 Jul 2022 - C. Keller   - Initial version (refactored Chem_GridCompMod)
!  See https://github.com/geoschem/geos-chem for history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! LOCAL VARIABLES:
!
    TYPE(ESMF_CONFIG)            :: CF
    CHARACTER(LEN=ESMF_MAXSTR)   :: ProviderName 

    __Iam__('GEOS_AeroSetServices')

    !=======================================================================
    ! GEOS_AeroSetServices starts here
    !=======================================================================

    ! Get configuration
    CALL ESMF_GridCompGet( GC, CONFIG = CF, __RC__ )

    ! See if GC is the AERO provider
    DoAERO = .FALSE.
    CALL ESMF_ConfigGetAttribute( CF, ProviderName,       &
                                  Label="AERO_PROVIDER:", &
                                  Default="PCHEM",        &
                                  __RC__                   )
    IF ( ProviderName == "GEOSCHEMCHEM" ) DoAERO = .TRUE.

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

    _RETURN(ESMF_SUCCESS)

    END SUBROUTINE GEOS_AeroSetServices
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Model                            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: GEOS_AeroInit
!
! !DESCRIPTION: Routine to initialize the aerosol bundle. 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE GEOS_AeroInit( GC, MaplCF, INTERNAL, EXPORT, Grid, RC ) 
!
! !USES:
!
    USE DiagList_Mod,       ONLY : SPFX
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ESMF_GridComp), INTENT(INOUT), TARGET :: GC       ! Ref to this GridComp
    TYPE(ESMF_Config),   INTENT(INOUT)         :: MaplCF   ! MAPL.rc
    TYPE(ESMF_State),    INTENT(INOUT)         :: INTERNAL ! Internal State 
    TYPE(ESMF_State),    INTENT(INOUT)         :: Export   ! Export State
    TYPE(ESMF_Grid),     INTENT(INOUT)         :: Grid     ! Grid object
    INTEGER,             INTENT(INOUT)         :: RC       ! Success or failure?
!
! !REMARKS:
!
!
! !REMARKS:
!
! !REVISION HISTORY:
!  12 Jul 2022 - C. Keller   - Initial version (refactored Chem_GridCompMod)
!  See https://github.com/geoschem/geos-chem for history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! LOCAL VARIABLES:
!
    ! Aerosol bundle
    REAL                         :: GCMW, FRAC
    INTEGER                      :: instance
    INTEGER                      :: I, J, NFLDS
    INTEGER                      :: GCID
    TYPE(ESMF_STATE)             :: Aero
    TYPE(ESMF_FieldBundle)       :: AeroBdl
    TYPE(ESMF_Field)             :: AeroFld, GcFld
    LOGICAL                      :: DynFriend, FRIENDLY
    LOGICAL                      :: IsPresent
    CHARACTER(LEN=ESMF_MAXSTR)   :: GCName, AeroName, fieldName
    REAL, POINTER                :: Ptr3D(:,:,:) => NULL()

    __Iam__('GEOS_AeroInit')

    !=======================================================================
    ! GEOS_AeroInit starts here
    !=======================================================================

       ! Get AERO bundle
       CALL ESMF_StateGet( EXPORT, 'AERO', Aero, __RC__ )

       ! This attribute indicates if the aerosol optics method is implemented
       ! or not. Radiation will not call the aerosol optics method unless this
       ! attribute is explicitly set to true.
       call ESMF_AttributeSet(Aero, name='implements_aerosol_optics_method', &
                              value=.true., __RC__)
       AeroBdl = ESMF_FieldBundleCreate(name='AEROSOLS', __RC__)
       call MAPL_StateAdd(Aero, AeroBdl, __RC__)

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
             CALL ESMF_StateGet( INTERNAL, TRIM(GCName), GcFld, RC=RC )

             ! Try TRC_<NAME> if SPC_<NAME> not found
             IF ( RC /= ESMF_SUCCESS ) THEN
                GCName = 'TRC_'//TRIM(GcNames(I))
                CALL ESMF_StateGet( INTERNAL, TRIM(GCName), GcFld, RC=RC )
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
                           userRoutine=GEOS_AerosolOptics, __RC__)

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

    ! Return
    _RETURN(ESMF_SUCCESS)

    END SUBROUTINE GEOS_AeroInit
!EOC
!------------------------------------------------------------------------------








!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Model                            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: GEOS_FillAeroBundle
!
! !DESCRIPTION: Routine to fill the aerosol bundle. 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE GEOS_FillAeroBundle( GC, EXPORT, State_Chm, State_Grid, Input_Opt, RC ) 
!
! !USES:
!
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ESMF_GridComp), INTENT(INOUT), TARGET :: GC     ! Ref to this GridComp
    TYPE(ESMF_State),    INTENT(INOUT)         :: Export   ! Export State
    TYPE(ChmState),      INTENT(INOUT)         :: State_Chm
    TYPE(GrdState),      INTENT(INOUT)         :: State_Grid
    TYPE(OptInput),      INTENT(INOUT)         :: Input_Opt
!
! !OUTPUT PARAMETERS:
!
    INTEGER,             INTENT(INOUT)         :: RC       ! Success or failure?
!
! !REMARKS:
!
! !REVISION HISTORY:
!  05 Jul 2022 - C. Keller   - Initial version (refactored Chem_GridCompMod)
!  See https://github.com/geoschem/geos-chem for history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! LOCAL VARIABLES:
!
    LOGICAL                      :: am_I_Root
    REAL                         :: GCMW, FRAC
    INTEGER                      :: N, IndSpc
    INTEGER                      :: nAero, nLen, GCID
    TYPE(ESMF_STATE)             :: Aero
    TYPE(ESMF_FieldBundle)       :: AeroBdl
    TYPE(ESMF_Field)             :: AeroFld
    CHARACTER(LEN=ESMF_MAXSTR)   :: GCName, AeroName
    REAL(fp), POINTER            :: GcPtr3d  (:,:,:) => NULL()
    REAL, POINTER                :: AeroPtr3d(:,:,:) => NULL()

    INTEGER, PARAMETER           :: NRATS = 5
    CHARACTER(LEN=15), PARAMETER :: RatsNames(NRATS) = (/ 'CH4', 'N2O', 'CFC11', 'CFC12', 'HCFC22' /)

    __Iam__('GEOS_FillAeroBundle')

    !=======================================================================
    ! GEOS_FillAeroBundle starts here
    !=======================================================================

    ! Are we on the root PET?
    am_I_Root = MAPL_Am_I_Root()

    ! For every field of the AERO bundle, we will copy the corresponding
    ! GEOS-Chem tracer field, converting units from mol mol-1 to kg kg-1.

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
       nlen = LEN(TRIM(GcName))
       IndSpc = Ind_(TRIM(GcName(5:nlen)))
       ASSERT_(IndSpc>0)
       GcPtr3D => State_Chm%Species(IndSpc)%Conc(:,:,State_Grid%NZ:1:-1)
       !CALL MAPL_GetPointer ( INTSTATE, GcPtr3D, TRIM(GcName), __RC__ )

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

       ! Free pointers
       GcPtr3D   => NULL()
       AeroPtr3D => NULL()
    ENDDO

    ! Fill AERO_DP bundle
    CALL FillAeroDP_ ( am_I_Root, GC, EXPORT, Input_Opt, __RC__ )

    _RETURN(ESMF_SUCCESS)

    END SUBROUTINE GEOS_FillAeroBundle
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Model                            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: GEOS_AerosolOptics 
!
! !DESCRIPTION: Aerosol optics routine, adapted from GOCART 
!\\
!\\
! !INTERFACE:
!
  ! Adapted from the GOCART interface
  subroutine GEOS_AerosolOptics(state, rc)
!
! !USES:
!
!
! !PARAMETERS:
!
    type(ESMF_State)     :: state
    integer, intent(out) :: rc
!
! !REVISION HISTORY:
!  06 Jul 2022 - C. Keller   - Initial version (from Chem_GridCompMod)
!  See https://github.com/geoschem/geos-chem for history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! LOCAL VARIABLES:
!
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

    Iam = 'GEOSCHEMCHEM::GEOS_AerosolOptics()'

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

  end subroutine GEOS_AerosolOptics 
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Model                            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: FillAeroDP_
!
! !DESCRIPTION: FillAeroDP_ fills the AERO_DP bundle
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE FillAeroDP_ ( am_I_Root, GC, EXPORT, Input_Opt, RC )
!
! !USES:
!
    USE HCO_ERROR_MOD
    USE HCO_TYPES_MOD,     ONLY : DiagnCont
    USE HCO_DIAGN_MOD,     ONLY : Diagn_Get
    USE HCO_State_GC_Mod,  ONLY : HcoState
!
! !INPUT PARAMETERS:
!
    LOGICAL                            :: am_I_Root
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ESMF_GridComp), INTENT(INOUT) :: GC       ! Ref to this GridComp
    TYPE(ESMF_State),    INTENT(INOUT) :: Export   ! Export State
    TYPE(OptInput),      INTENT(INOUT) :: Input_Opt
!
! !OUTPUT PARAMETERS:
!
    INTEGER, INTENT(OUT), OPTIONAL     :: RC
!
! !REVISION HISTORY:
!  30 Mar 2015 - C. Keller   - Initial version
!  See https://github.com/geoschem/geos-chem for history
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
    ! FillAeroDP_ begins here
    !=======================================================================

    ! Traceback handle
    Iam = 'FillAeroDP_'

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
             CALL Diagn_Get( HcoState, .FALSE., DgnCont,  &
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
                CALL Diagn_Get( HcoState, .FALSE., DgnCont,  &
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

  END SUBROUTINE FillAeroDP_
!EOC
END MODULE GEOS_AeroCoupler
