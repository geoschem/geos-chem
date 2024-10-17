!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: aermass_container_mod.F90
!
! !DESCRIPTION: Module AERMASS\_CONTAINER\_MOD contains the derived type used
!  to store aerosol data used for computing optical depths and diagnostics.
!\\
!\\
! !INTERFACE:
!
MODULE AerMass_Container_Mod
!
! USES:
!
  USE ErrCode_Mod
  USE Precision_Mod

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC :: Init_AerMass_Container
  PUBLIC :: Cleanup_AerMass_Container
!
! !PUBLIC DATA MEMBERS:
!
  !=========================================================================
  ! Derived type for aerosol data container
  !=========================================================================
  TYPE, PUBLIC :: AerMassContainer

     ! Logical flags
     LOGICAL  :: Is_OCPI
     LOGICAL  :: Is_OCPO
     LOGICAL  :: Is_BC
     LOGICAL  :: Is_SO4
     LOGICAL  :: Is_HMS
     LOGICAL  :: Is_NH4
     LOGICAL  :: Is_NIT
     LOGICAL  :: Is_DST
     LOGICAL  :: Is_SAL
     LOGICAL  :: Is_POA
     LOGICAL  :: Is_OPOA
     LOGICAL  :: Is_TSOA
     LOGICAL  :: Is_ASOA
     LOGICAL  :: Is_SOAGX
     LOGICAL  :: Is_SimpleSOA
     LOGICAL  :: Is_ComplexSOA

     ! Arrays
     !========================================================================
     ! OCFPOA      : OM/OC for POA                      [unitless] - used in carbon_mod
     ! OCFOPOA     : OM/OC for OPOA, OCPI, OCPO         [unitless] - used in carbon_mod
     ! BCPI        : Hydrophilic black carbon aerosol   [kg/m3]
     ! BCPO        : Hydrophobic black carbon aerosol   [kg/m3]
     ! OCPI        : Hydrophilic organic carbon aerosol [kg/m3]
     ! OCPO        : Hydrophobic organic carbon aerosol [kg/m3]
     ! OCPISOA     : Hydrophilic OC + SOA aerosol       [kg/m3]
     ! SALA        : Accumulation mode seasalt aerosol  [kg/m3]
     ! ACL         : Accumulation mode Cl aerosol       [kg/m3]
     ! SALC        : Coarse mode seasalt aerosol        [kg/m3]
     ! SO4_NH4_NIT : Lumped SO4-NH4-NIT aerosol         [kg/m3]
     ! SO4         : Sulfate aerosol                    [kg/m3]
     ! HMS         : Hydroxymethane sulfonate aerosol   [kg/m3]
     ! NH4         : Ammonium aerosol                   [kg/m3]
     ! NIT         : Inorganic nitrate aerosol          [kg/m3]
     ! SLA         : Stratospheric liquid aerosol       [kg/m3]
     ! SPA         : Stratospheric particulate aerosol  [kg/m3]
     ! TSOA        : Terpene SOA                        [kg/m3]
     ! ASOA        : Aromatic + IVOC SOA                [kg/m3]
     ! OPOA        : Aerosol product of SVOC oxidation  [kg/m3]
     ! SOAGX       : SOA product of GLYX                [kg/m3]
     ! SOAIE       : SOA product of IEPOX & HMML        [kg/m3]
     ! PM25        : Particulate matter < 2.5 um        [kg/m3]
     ! PM10        : Particulate matter < 10 um        [kg/m3]
     ! PDER        : Parameterized effective Radius for SNA and OM [nm] - used for AOD calcualtion (H. Zhu)
     ! ISOAAQ      : Isoprene SOA (aqueous formation)   [kg/m3]
     ! SOAS        : Simple SOA                         [kg/m3]
     ! FRAC_SNA    :
     ! DAERSL      : Mass of hydrophobic aerosol (Mian Chin)
     ! WAERSL      : Mass of hydrophilic aerosol (Mian Chin)
     !========================================================================
     REAL(fp), POINTER :: OCFPOA     (:,:)
     REAL(fp), POINTER :: OCFOPOA    (:,:)
     REAL(fp), POINTER :: BCPI       (:,:,:)
     REAL(fp), POINTER :: BCPO       (:,:,:)
     REAL(fp), POINTER :: OCPI       (:,:,:)
     REAL(fp), POINTER :: OCPO       (:,:,:)
     REAL(fp), POINTER :: OCPISOA    (:,:,:)
     REAL(fp), POINTER :: SALA       (:,:,:)
     REAL(fp), POINTER :: ACL        (:,:,:)
     REAL(fp), POINTER :: SALC       (:,:,:)
     REAL(fp), POINTER :: SO4_NH4_NIT(:,:,:)
     REAL(fp), POINTER :: SO4        (:,:,:)
     REAL(fp), POINTER :: HMS        (:,:,:)
     REAL(fp), POINTER :: NH4        (:,:,:)
     REAL(fp), POINTER :: NIT        (:,:,:)
     REAL(fp), POINTER :: SLA        (:,:,:)
     REAL(fp), POINTER :: SPA        (:,:,:)
     REAL(fp), POINTER :: TSOA       (:,:,:)
     REAL(fp), POINTER :: ASOA       (:,:,:)
     REAL(fp), POINTER :: OPOA       (:,:,:)
     REAL(fp), POINTER :: SOAGX      (:,:,:)
     REAL(fp), POINTER :: PM25       (:,:,:)
     REAL(fp), POINTER :: PM10       (:,:,:)
     REAL(fp), POINTER :: PDER       (:,:,:) !H. Zhu
     REAL(fp), POINTER :: SNAOM      (:,:,:) !H. Zhu
     REAL(fp), POINTER :: R_OMSNA    (:,:,:) !H. Zhu
     REAL(fp), POINTER :: ISOAAQ     (:,:,:)
     REAL(fp), POINTER :: SOAS       (:,:,:)
     REAL(fp), POINTER :: FRAC_SNA   (:,:,:,:)
     REAL(fp), POINTER :: DAERSL     (:,:,:,:)
     REAL(fp), POINTER :: WAERSL     (:,:,:,:)

  END TYPE AerMassContainer
!
! !REMARKS:
! 
! !REVISION HISTORY:
!  28 Mar 2023 - E. Lundgren- Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Init_AerMass_Container
!
! !DESCRIPTION: Subroutine INIT\_AER\_Container allocates and initializes
! the Aer container object.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Init_AerMass_Container( Input_Opt, State_Grid, Aer, RC )
!
! !USES:
!
    USE CMN_Size_Mod,   ONLY : NAER
    USE Input_Opt_Mod,  ONLY : OptInput
    USE State_Grid_Mod, ONLY : GrdState
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput),      INTENT(IN)  :: Input_Opt  ! Input Options object
    TYPE(GrdState),      INTENT(IN)  :: State_Grid ! Grid object
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(AerMassContainer),  POINTER :: Aer        ! Aerosol data container
!
! !OUTPUT PARAMETERS:
!
    INTEGER,             INTENT(OUT) :: RC         ! Success or failure?
!
! !REVISION HISTORY:
!  28 Mar 2023 - E. Lundgren- Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    CHARACTER(LEN=255) :: errMsg, thisLoc
    INTEGER            :: NX, NY, NZ

    !======================================================================
    ! Init_AerMass_Container starts here
    !======================================================================

    ! Initialize local variables
    RC = GC_SUCCESS
    thisLoc = ' -> at Init_AerMass_Container (in module Headers/aermass_container_mod.F90)'
    NX = State_Grid%NX
    NY = State_Grid%NY
    NZ = State_Grid%NZ

    ! Exit immediately if this is a dry-run
    IF ( Input_Opt%DryRun ) RETURN

    !======================================================================
    ! Initialize logical flags (will be updated during init_aerosol
    !======================================================================
    Aer%Is_OCPI       = .FALSE.
    Aer%Is_OCPO       = .FALSE.
    Aer%Is_BC         = .FALSE.
    Aer%Is_SO4        = .FALSE.
    Aer%Is_HMS        = .FALSE.
    Aer%Is_NH4        = .FALSE.
    Aer%Is_NIT        = .FALSE.
    Aer%Is_DST        = .FALSE.
    Aer%Is_SAL        = .FALSE.
    Aer%Is_POA        = .FALSE.
    Aer%Is_OPOA       = .FALSE.
    Aer%Is_TSOA       = .FALSE.
    Aer%Is_ASOA       = .FALSE.
    Aer%Is_SOAGX      = .FALSE.
    Aer%Is_SimpleSOA  = .FALSE.
    Aer%Is_ComplexSOA = .FALSE.

    !======================================================================
    ! Initialize arrays
    !======================================================================

    ! OM/OC for POA
    ALLOCATE( Aer%OCFPOA( NX, NY ), STAT=RC )
    CALL GC_CheckVar( 'OCFPOA', 0, RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error allocating array OCFPOA!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    Aer%OCFPOA = 0.0_fp

    ! OM/OC for OPOA, OCPI, OCPO
    ALLOCATE( Aer%OCFOPOA( NX, NY ), STAT=RC )
    CALL GC_CheckVar( 'OCFOPOA', 0, RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error allocating array OCFOPOA!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    Aer%OCFOPOA = 0.0_fp

    ALLOCATE( Aer%BCPI( NX, NY, NZ ), STAT=RC )
    CALL GC_CheckVar( 'BCPI', 0, RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error allocating array BCPI!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    Aer%BCPI = 0.0_fp

    ALLOCATE( Aer%BCPO( NX, NY, NZ ), STAT=RC )
    CALL GC_CheckVar( 'BCPO', 0, RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error allocating array BCPO!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    Aer%BCPO = 0.0_fp

    ALLOCATE( Aer%OCPI( NX, NY, NZ ), STAT=RC )
    CALL GC_CheckVar( 'OCPI', 0, RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error allocating array OCPI!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    Aer%OCPI = 0.0_fp

    ALLOCATE( Aer%OCPO( NX, NY, NZ ), STAT=RC )
    CALL GC_CheckVar( 'OCPO', 0, RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error allocating array OCPO!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    Aer%OCPO = 0.0_fp

    ALLOCATE( Aer%OCPISOA( NX, NY, NZ ), STAT=RC )
    CALL GC_CheckVar( 'OCPISOA', 0, RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error allocating array OCPISOA!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    Aer%OCPISOA = 0.0_fp

    ALLOCATE( Aer%SALA( NX, NY, NZ ), STAT=RC )
    CALL GC_CheckVar( 'SALA', 0, RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error allocating array SALA!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    Aer%SALA = 0.0_fp

    ALLOCATE( Aer%SALC( NX, NY, NZ ), STAT=RC )
    CALL GC_CheckVar( 'SALC', 0, RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error allocating array SALC!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    Aer%SALC = 0.0_fp

    ALLOCATE( Aer%ACL( NX, NY, NZ ),  STAT=RC )
    CALL GC_CheckVar( 'ACL', 0, RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error allocating array ACL!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    Aer%ACL = 0.0_fp

    ALLOCATE( Aer%SO4_NH4_NIT( NX, NY, NZ ), STAT=RC )
    CALL GC_CheckVar( 'SO4_NH4_NIT', 0, RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error allocating array SO4_NH4_NIT!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    Aer%SO4_NH4_NIT = 0.0_fp

    ALLOCATE( Aer%SO4( NX, NY, NZ ), STAT=RC )
    CALL GC_CheckVar( 'SO4', 0, RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error allocating array SO4!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    Aer%SO4 = 0.0_fp

    ! (jmm, 06/30/18)
    ALLOCATE( Aer%HMS( NX, NY, NZ ), STAT=RC )
    CALL GC_CheckVar( 'HMS', 0, RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error allocating array HMS!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    Aer%HMS = 0.0_fp

    ALLOCATE( Aer%NH4( NX, NY, NZ ), STAT=RC )
    CALL GC_CheckVar( 'NH4', 0, RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error allocating array NH4!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    Aer%NH4 = 0.0_fp

    ALLOCATE( Aer%NIT( NX, NY, NZ ), STAT=RC )
    CALL GC_CheckVar( 'NIT', 0, RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error allocating array NIT!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    Aer%NIT = 0.0_fp

    ALLOCATE( Aer%SLA( NX, NY, NZ ), STAT=RC )
    CALL GC_CheckVar( 'SLA', 0, RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error allocating array SLA!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    Aer%SLA = 0.0_fp

    ALLOCATE( Aer%SPA( NX, NY, NZ ), STAT=RC )
    CALL GC_CheckVar( 'SPA', 0, RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error allocating array SPA!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    Aer%SPA   = 0.0_fp

    ALLOCATE( Aer%TSOA( NX, NY, NZ ), STAT=RC )
    CALL GC_CheckVar( 'TSOA', 0, RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error allocating array TSOA!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    Aer%TSOA = 0.0_fp

    ALLOCATE( Aer%ASOA( NX, NY, NZ ), STAT=RC )
    CALL GC_CheckVar( 'ASOA', 0, RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error allocating array ASOA!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    Aer%ASOA = 0.0_fp

    ALLOCATE( Aer%OPOA( NX, NY, NZ ), STAT=RC )
    CALL GC_CheckVar( 'OPOA', 0, RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error allocating array OPOA!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    Aer%OPOA = 0.0_fp

    ALLOCATE( Aer%PM25( NX, NY, NZ ), STAT=RC )
    CALL GC_CheckVar( 'PM25', 0, RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error allocating array PM25!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    Aer%PM25 = 0.0_fp

    !zhaisx
    ALLOCATE( Aer%PM10( NX, NY, NZ ), STAT=RC )
    CALL GC_CheckVar( 'PM10', 0, RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error allocating array PM10!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    Aer%PM10 = 0.0_fp
    ! H. Zhu
    ALLOCATE( Aer%PDER( NX, NY, NZ ), STAT=RC )
    CALL GC_CheckVar( 'PDER', 0, RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error allocating array PDER!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    Aer%PDER = 0.0_fp
    ! H. Zhu
    ALLOCATE( Aer%SNAOM( NX, NY, NZ ), STAT=RC )
    CALL GC_CheckVar( 'SNAOM', 0, RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error allocating array SNAOM!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    Aer%SNAOM = 0.0_fp    
    ! H. Zhu
    ALLOCATE( Aer%R_OMSNA( NX, NY, NZ ), STAT=RC )
    CALL GC_CheckVar( 'R_OMSNA', 0, RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error allocating array R_OMSNA!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    Aer%R_OMSNA = 0.0_fp


    ALLOCATE( Aer%SOAGX( NX, NY, NZ ), STAT=RC )
    CALL GC_CheckVar( 'SOAGX', 0, RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error allocating array SOAGX!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    Aer%SOAGX = 0.0_fp

    ! Mechanistic isoprene SOA (eam, 2014):
    ALLOCATE( Aer%ISOAAQ( NX, NY, NZ ), STAT=RC )
    CALL GC_CheckVar( 'ISOAAQ', 0, RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error allocating array ISOAAQ!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    Aer%ISOAAQ = 0.0_fp

    ! Simple SOA
    ALLOCATE( Aer%SOAS( NX, NY, NZ ), STAT=RC )
    CALL GC_CheckVar( 'SOAS', 0, RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error allocating array SOAS!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    Aer%SOAS = 0.0_fp

    ALLOCATE( Aer%FRAC_SNA( NX, NY, NZ, 3 ), STAT=RC )
    CALL GC_CheckVar( 'FRAC_SNA', 0, RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error allocating array FRAC_SNA!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    Aer%FRAC_SNA = 0.0_fp

    ! Mass of hydrophobic aerosol from Mian Chin
    ALLOCATE( Aer%DAERSL( NX, NY, NZ, 2 ), STAT=RC )
    CALL GC_CheckVar( 'DAERSL', 0, RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error allocating array DAERSL!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    Aer%DAERSL = 0.0_fp

    ! Mass of hydrophilic aerosol from Mian Chin
    ALLOCATE( Aer%WAERSL( NX, NY, NZ, NAER ), STAT=RC )
    CALL GC_CheckVar( 'WAERSL', 0, RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error allocating array WAERSL!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    Aer%WAERSL = 0.0_fp

  END SUBROUTINE Init_AerMass_Container
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Cleanup_AerMass_Container
!
! !DESCRIPTION: Subroutine CLEANUP\_AER\_CONTAINER deallocates all fields
!  of the aer container object.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Cleanup_AerMass_Container( Aer, RC )
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(AerMassContainer), POINTER :: Aer   ! Aer data container
!
! !OUTPUT PARAMETERS:
!
    INTEGER,            INTENT(OUT) :: RC    ! Return code
!
! !REVISION HISTORY:
!  28 Mar 2023 - E. Lundgren- Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC

    !======================================================================
    ! Cleanup_AerMass_Container starts here
    !======================================================================

    ! Assume success
    RC = GC_SUCCESS

    ! Deallocate arrays and nullify pointer
    IF ( ASSOCIATED( Aer%OCFPOA ) ) THEN
       DEALLOCATE( Aer%OCFPOA, STAT=RC )
       CALL GC_CheckVar( 'Aer%OCFPOA', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       Aer%OCFPOA => NULL()
    ENDIF

    IF ( ASSOCIATED( Aer%OCFOPOA ) ) THEN
       DEALLOCATE( Aer%OCFOPOA, STAT=RC )
       CALL GC_CheckVar( 'Aer%OCFOPOA', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       Aer%OCFOPOA => NULL()
    ENDIF

    IF ( ASSOCIATED( Aer%BCPI ) ) THEN
       DEALLOCATE( Aer%BCPI, STAT=RC )
       CALL GC_CheckVar( 'Aer%BCPI', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       Aer%BCPI => NULL()
    ENDIF

    IF ( ASSOCIATED( Aer%BCPO ) ) THEN
       DEALLOCATE( Aer%BCPO, STAT=RC )
       CALL GC_CheckVar( 'Aer%BCPO', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       Aer%BCPO => NULL()
    ENDIF

    IF ( ASSOCIATED( Aer%OCPI ) ) THEN
       DEALLOCATE( Aer%OCPI, STAT=RC )
       CALL GC_CheckVar( 'Aer%OCPI', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       Aer%OCPI => NULL()
    ENDIF

    IF ( ASSOCIATED( Aer%OCPO ) ) THEN
       DEALLOCATE( Aer%OCPO, STAT=RC )
       CALL GC_CheckVar( 'Aer%OCPO', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       Aer%OCPO => NULL()
    ENDIF

    IF ( ASSOCIATED( Aer%OCPISOA ) ) THEN
       DEALLOCATE( Aer%OCPISOA, STAT=RC )
       CALL GC_CheckVar( 'Aer%OCPISOA', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       Aer%OCPISOA => NULL()
    ENDIF

    IF ( ASSOCIATED( Aer%SALA ) ) THEN
       DEALLOCATE( Aer%SALA, STAT=RC )
       CALL GC_CheckVar( 'Aer%SALA', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       Aer%SALA => NULL()
    ENDIF

    IF ( ASSOCIATED( Aer%ACL ) ) THEN
       DEALLOCATE( Aer%ACL, STAT=RC )
       CALL GC_CheckVar( 'Aer%ACL', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       Aer%ACL => NULL()
    ENDIF

    IF ( ASSOCIATED( Aer%SALC ) ) THEN
       DEALLOCATE( Aer%SALC, STAT=RC )
       CALL GC_CheckVar( 'Aer%SALC', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       Aer%SALC => NULL()
    ENDIF

    IF ( ASSOCIATED( Aer%SO4_NH4_NIT ) ) THEN
       DEALLOCATE( Aer%SO4_NH4_NIT, STAT=RC )
       CALL GC_CheckVar( 'Aer%SO4_NH4_NIT', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       Aer%SO4_NH4_NIT => NULL()
    ENDIF

    IF ( ASSOCIATED( Aer%SO4 ) ) THEN
       DEALLOCATE( Aer%SO4, STAT=RC )
       CALL GC_CheckVar( 'Aer%SO4', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       Aer%SO4 => NULL()
    ENDIF

    IF ( ASSOCIATED( Aer%HMS ) ) THEN
       DEALLOCATE( Aer%HMS, STAT=RC )
       CALL GC_CheckVar( 'Aer%HMS', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       Aer%HMS => NULL()
    ENDIF

    IF ( ASSOCIATED( Aer%NH4 ) ) THEN
       DEALLOCATE( Aer%NH4, STAT=RC )
       CALL GC_CheckVar( 'Aer%NH4', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       Aer%NH4 => NULL()
    ENDIF

    IF ( ASSOCIATED( Aer%NIT ) ) THEN
       DEALLOCATE( Aer%NIT, STAT=RC )
       CALL GC_CheckVar( 'Aer%NIT', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       Aer%NIT => NULL()
    ENDIF

    IF ( ASSOCIATED( Aer%SLA ) ) THEN
       DEALLOCATE( Aer%SLA, STAT=RC )
       CALL GC_CheckVar( 'Aer%SLA', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       Aer%SLA => NULL()
    ENDIF

    IF ( ASSOCIATED( Aer%SPA ) ) THEN
       DEALLOCATE( Aer%SPA, STAT=RC )
       CALL GC_CheckVar( 'Aer%SPA', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       Aer%SPA => NULL()
    ENDIF

    IF ( ASSOCIATED( Aer%TSOA ) ) THEN
       DEALLOCATE( Aer%TSOA, STAT=RC )
       CALL GC_CheckVar( 'Aer%TSOA', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       Aer%TSOA => NULL()
    ENDIF

    IF ( ASSOCIATED( Aer%ASOA ) ) THEN
       DEALLOCATE( Aer%ASOA, STAT=RC )
       CALL GC_CheckVar( 'Aer%ASOA', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       Aer%ASOA => NULL()
    ENDIF

    IF ( ASSOCIATED( Aer%OPOA ) ) THEN
       DEALLOCATE( Aer%OPOA, STAT=RC )
       CALL GC_CheckVar( 'Aer%OPOA', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       Aer%OPOA => NULL()
    ENDIF

    IF ( ASSOCIATED( Aer%SOAGX ) ) THEN
       DEALLOCATE( Aer%SOAGX, STAT=RC )
       CALL GC_CheckVar( 'Aer%SOAGX', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       Aer%SOAGX => NULL()
    ENDIF

    IF ( ASSOCIATED( Aer%PM25 ) ) THEN
       DEALLOCATE( Aer%PM25, STAT=RC )
       CALL GC_CheckVar( 'Aer%PM25', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       Aer%PM25 => NULL()
    ENDIF

    IF ( ASSOCIATED( Aer%PM10 ) ) THEN
       DEALLOCATE( Aer%PM10, STAT=RC )
       CALL GC_CheckVar( 'Aer%PM10', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       Aer%PM10 => NULL()
    ENDIF
    ! H. Zhu
    IF ( ASSOCIATED( Aer%PDER ) ) THEN
       DEALLOCATE( Aer%PDER, STAT=RC )
       CALL GC_CheckVar( 'Aer%PDER', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       Aer%PDER => NULL()
    ENDIF
    ! H. Zhu
    IF ( ASSOCIATED( Aer%SNAOM ) ) THEN
       DEALLOCATE( Aer%SNAOM, STAT=RC )
       CALL GC_CheckVar( 'Aer%SNAOM', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       Aer%SNAOM => NULL()
    ENDIF    
    ! H. Zhu
    IF ( ASSOCIATED( Aer%R_OMSNA ) ) THEN
       DEALLOCATE( Aer%R_OMSNA, STAT=RC )
       CALL GC_CheckVar( 'Aer%R_OMSNA', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       Aer%R_OMSNA => NULL()
    ENDIF

    IF ( ASSOCIATED( Aer%ISOAAQ ) ) THEN
       DEALLOCATE( Aer%ISOAAQ, STAT=RC )
       CALL GC_CheckVar( 'Aer%ISOAAQ', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       Aer%ISOAAQ => NULL()
    ENDIF

    IF ( ASSOCIATED( Aer%SOAS ) ) THEN
       DEALLOCATE( Aer%SOAS, STAT=RC )
       CALL GC_CheckVar( 'Aer%SOAS', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       Aer%SOAS => NULL()
    ENDIF

    IF ( ASSOCIATED( Aer%FRAC_SNA ) ) THEN
       DEALLOCATE( Aer%FRAC_SNA, STAT=RC )
       CALL GC_CheckVar( 'Aer%FRAC_SNA', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       Aer%FRAC_SNA => NULL()
    ENDIF

    IF ( ASSOCIATED( Aer%DAERSL ) ) THEN
       DEALLOCATE( Aer%DAERSL, STAT=RC )
       CALL GC_CheckVar( 'Aer%DAERSL', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       Aer%DAERSL => NULL()
    ENDIF

    IF ( ASSOCIATED( Aer%WAERSL ) ) THEN
       DEALLOCATE( Aer%WAERSL, STAT=RC )
       CALL GC_CheckVar( 'Aer%WAERSL', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       Aer%WAERSL => NULL()
    ENDIF

  END SUBROUTINE Cleanup_AerMass_Container
!EOC

END MODULE AerMass_Container_Mod
