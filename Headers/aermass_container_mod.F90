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
     ! ISOAAQ      : Isoprene SOA (aqueous formation)   [kg/m3]
     ! SOAS        : Simple SOA                         [kg/m3]
     ! FRAC_SNA    :
     ! DAERSL      : Mass of hydrophobic aerosol (Mian Chin)
     ! WAERSL      : Mass of hydrophilic aerosol (Mian Chin)
     !========================================================================
     REAL(fp), ALLOCATABLE :: OCFPOA     (:,:)
     REAL(fp), ALLOCATABLE :: OCFOPOA    (:,:)
     REAL(fp), ALLOCATABLE :: BCPI       (:,:,:)
     REAL(fp), ALLOCATABLE :: BCPO       (:,:,:)
     REAL(fp), ALLOCATABLE :: OCPI       (:,:,:)
     REAL(fp), ALLOCATABLE :: OCPO       (:,:,:)
     REAL(fp), ALLOCATABLE :: OCPISOA    (:,:,:)
     REAL(fp), ALLOCATABLE :: SALA       (:,:,:)
     REAL(fp), ALLOCATABLE :: ACL        (:,:,:)
     REAL(fp), ALLOCATABLE :: SALC       (:,:,:)
     REAL(fp), ALLOCATABLE :: SO4_NH4_NIT(:,:,:)
     REAL(fp), ALLOCATABLE :: SO4        (:,:,:)
     REAL(fp), ALLOCATABLE :: HMS        (:,:,:)
     REAL(fp), ALLOCATABLE :: NH4        (:,:,:)
     REAL(fp), ALLOCATABLE :: NIT        (:,:,:)
     REAL(fp), ALLOCATABLE :: SLA        (:,:,:)
     REAL(fp), ALLOCATABLE :: SPA        (:,:,:)
     REAL(fp), ALLOCATABLE :: TSOA       (:,:,:)
     REAL(fp), ALLOCATABLE :: ASOA       (:,:,:)
     REAL(fp), ALLOCATABLE :: OPOA       (:,:,:)
     REAL(fp), ALLOCATABLE :: SOAGX      (:,:,:)
     REAL(fp), ALLOCATABLE :: PM25       (:,:,:)
     REAL(fp), ALLOCATABLE :: PM10       (:,:,:)
     REAL(fp), ALLOCATABLE :: ISOAAQ     (:,:,:)
     REAL(fp), ALLOCATABLE :: SOAS       (:,:,:)
     REAL(fp), ALLOCATABLE :: FRAC_SNA   (:,:,:,:)
     REAL(fp), ALLOCATABLE :: DAERSL     (:,:,:,:)
     REAL(fp), ALLOCATABLE :: WAERSL     (:,:,:,:)

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
    ! Initialize arrays
    !======================================================================

    ! OM/OC for POA
    ALLOCATE( Aer%OCFPOA( NX, NY ), STAT=RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error allocating array OCFPOA!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    Aer%OCFPOA = 0.0_fp

    ! OM/OC for OPOA, OCPI, OCPO
    ALLOCATE( Aer%OCFOPOA( NX, NY ), STAT=RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error allocating array OCFOPOA!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    Aer%OCFOPOA = 0.0_fp

    ALLOCATE( Aer%BCPI( NX, NY, NZ ), STAT=RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error allocating array BCPI!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    Aer%BCPI = 0.0_fp

    ALLOCATE( Aer%BCPO( NX, NY, NZ ), STAT=RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error allocating array BCPO!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    Aer%BCPO = 0.0_fp

    ALLOCATE( Aer%OCPI( NX, NY, NZ ), STAT=RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error allocating array OCPI!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    Aer%OCPI = 0.0_fp

    ALLOCATE( Aer%OCPO( NX, NY, NZ ), STAT=RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error allocating array OCPO!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    Aer%OCPO = 0.0_fp

    ALLOCATE( Aer%OCPISOA( NX, NY, NZ ), STAT=RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error allocating array OCPISOA!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    Aer%OCPISOA = 0.0_fp

    ALLOCATE( Aer%SALA( NX, NY, NZ ), STAT=RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error allocating array SALA!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    Aer%SALA = 0.0_fp

    ALLOCATE( Aer%SALC( NX, NY, NZ ), STAT=RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error allocating array SALC!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    Aer%SALC = 0.0_fp

    ALLOCATE( Aer%ACL( NX, NY, NZ ),  STAT=RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error allocating array ACL!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    Aer%ACL = 0.0_fp

    ALLOCATE( Aer%SO4_NH4_NIT( NX, NY, NZ ), STAT=RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error allocating array SO4_NH4_NIT!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    Aer%SO4_NH4_NIT = 0.0_fp

    ALLOCATE( Aer%SO4( NX, NY, NZ ), STAT=RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error allocating array SO4!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    Aer%SO4 = 0.0_fp

    ! (jmm, 06/30/18)
    ALLOCATE( Aer%HMS( NX, NY, NZ ), STAT=RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error allocating array HMS!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    Aer%HMS = 0.0_fp

    ALLOCATE( Aer%NH4( NX, NY, NZ ), STAT=RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error allocating array NH4!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    Aer%NH4 = 0.0_fp

    ALLOCATE( Aer%NIT( NX, NY, NZ ), STAT=RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error allocating array NIT!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    Aer%NIT = 0.0_fp

    ALLOCATE( Aer%SLA( NX, NY, NZ ), STAT=RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error allocating array SLA!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    Aer%SLA = 0.0_fp

    ALLOCATE( Aer%SPA( NX, NY, NZ ), STAT=RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error allocating array SPA!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    Aer%SPA   = 0.0_fp

    ALLOCATE( Aer%TSOA( NX, NY, NZ ), STAT=RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error allocating array TSOA!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    Aer%TSOA = 0.0_fp

    ALLOCATE( Aer%ASOA( NX, NY, NZ ), STAT=RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error allocating array ASOA!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    Aer%ASOA = 0.0_fp

    ALLOCATE( Aer%OPOA( NX, NY, NZ ), STAT=RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error allocating array OPOA!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    Aer%OPOA = 0.0_fp

    ALLOCATE( Aer%PM25( NX, NY, NZ ), STAT=RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error allocating array PM25!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    Aer%PM25 = 0.0_fp

    !zhaisx
    ALLOCATE( Aer%PM10( NX, NY, NZ ), STAT=RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error allocating array PM10!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    Aer%PM10 = 0.0_fp

    ALLOCATE( Aer%SOAGX( NX, NY, NZ ), STAT=RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error allocating array SOAGX!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    Aer%SOAGX = 0.0_fp

    ! Mechanistic isoprene SOA (eam, 2014):
    ALLOCATE( Aer%ISOAAQ( NX, NY, NZ ), STAT=RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error allocating array ISOAAQ!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    Aer%ISOAAQ = 0.0_fp

    ! Simple SOA
    ALLOCATE( Aer%SOAS( NX, NY, NZ ), STAT=RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error allocating array SOAS!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    Aer%SOAS = 0.0_fp

    ALLOCATE( Aer%FRAC_SNA( NX, NY, NZ, 3 ), &
              STAT=RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error allocating array FRAC_SNA!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    Aer%FRAC_SNA = 0.0_fp

    ! Mass of hydrophobic aerosol from Mian Chin
    ALLOCATE( Aer%DAERSL( NX, NY, NZ, 2 ), STAT=RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error allocating array DAERSL!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    Aer%DAERSL = 0.0_fp

    ! Mass of hydrophilic aerosol from Mian Chin
    ALLOCATE( Aer%WAERSL( NX, NY, NZ, NAER ), STAT=RC )
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
    IF ( ASSOCIATED( Aer ) ) THEN
       IF (ALLOCATED(Aer%OCFPOA       )) DEALLOCATE(Aer%OCFPOA       )
       IF (ALLOCATED(Aer%OCFOPOA      )) DEALLOCATE(Aer%OCFOPOA      )
       IF (ALLOCATED(Aer%BCPI         )) DEALLOCATE(Aer%BCPI         )
       IF (ALLOCATED(Aer%BCPO         )) DEALLOCATE(Aer%BCPO         )
       IF (ALLOCATED(Aer%OCPI         )) DEALLOCATE(Aer%OCPI         )
       IF (ALLOCATED(Aer%OCPO         )) DEALLOCATE(Aer%OCPO         )
       IF (ALLOCATED(Aer%OCPISOA      )) DEALLOCATE(Aer%OCPISOA      )
       IF (ALLOCATED(Aer%SALA         )) DEALLOCATE(Aer%SALA         )
       IF (ALLOCATED(Aer%ACL          )) DEALLOCATE(Aer%ACL          )
       IF (ALLOCATED(Aer%SALC         )) DEALLOCATE(Aer%SALC         )
       IF (ALLOCATED(Aer%SO4_NH4_NIT  )) DEALLOCATE(Aer%SO4_NH4_NIT  )
       IF (ALLOCATED(Aer%SO4          )) DEALLOCATE(Aer%SO4          )
       IF (ALLOCATED(Aer%HMS          )) DEALLOCATE(Aer%HMS          )
       IF (ALLOCATED(Aer%NH4          )) DEALLOCATE(Aer%NH4          )
       IF (ALLOCATED(Aer%NIT          )) DEALLOCATE(Aer%NIT          )
       IF (ALLOCATED(Aer%SLA          )) DEALLOCATE(Aer%SLA          )
       IF (ALLOCATED(Aer%SPA          )) DEALLOCATE(Aer%SPA          )
       IF (ALLOCATED(Aer%TSOA         )) DEALLOCATE(Aer%TSOA         )
       IF (ALLOCATED(Aer%ASOA         )) DEALLOCATE(Aer%ASOA         )
       IF (ALLOCATED(Aer%OPOA         )) DEALLOCATE(Aer%OPOA         )
       IF (ALLOCATED(Aer%SOAGX        )) DEALLOCATE(Aer%SOAGX        )
       IF (ALLOCATED(Aer%PM25         )) DEALLOCATE(Aer%PM25         )
       IF (ALLOCATED(Aer%PM10         )) DEALLOCATE(Aer%PM10         )
       IF (ALLOCATED(Aer%ISOAAQ       )) DEALLOCATE(Aer%ISOAAQ       )
       IF (ALLOCATED(Aer%SOAS         )) DEALLOCATE(Aer%SOAS         )
       IF (ALLOCATED(Aer%FRAC_SNA     )) DEALLOCATE(Aer%FRAC_SNA     )
       IF (ALLOCATED(Aer%DAERSL       )) DEALLOCATE(Aer%DAERSL       )
       IF (ALLOCATED(Aer%WAERSL       )) DEALLOCATE(Aer%WAERSL       )
       Aer => NULL()
    ENDIF

  END SUBROUTINE Cleanup_AerMass_Container
!EOC

END MODULE AerMass_Container_Mod
