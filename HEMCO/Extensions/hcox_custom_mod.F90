!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: hcox_custom_mod.F90
!
! !DESCRIPTION: Customizable HEMCO emission extension. 
!\\
!\\
! !INTERFACE:
!
MODULE HCOX_Custom_Mod 
!
! !USES:
!
  USE HCO_Error_MOD
  USE HCO_Diagn_MOD
  USE HCOX_State_MOD, ONLY : Ext_State
  USE HCO_State_MOD,  ONLY : HCO_State 

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC :: HCOX_Custom_Run
  PUBLIC :: HCOX_Custom_Init
  PUBLIC :: HCOX_Custom_Final
!
! !REVISION HISTORY:
!  13 Dec 2013 - C. Keller   - Initial version 
!  06 Jun 2014 - R. Yantosca - Cosmetic changes in ProTeX headers
!  06 Jun 2014 - R. Yantosca - Now indented with F90 free-format
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !MODULE VARIABLES:
!
  INTEGER                         :: ExtNr   = -1
  INTEGER                         :: nOcWind = -1
  INTEGER                         :: nIceSrc = -1
  INTEGER,  ALLOCATABLE           :: OcWindIDs(:)
  INTEGER,  ALLOCATABLE           :: IceSrcIDs(:)
  REAL(hp), ALLOCATABLE, TARGET   :: FLUXICE(:,:)
  REAL(hp), ALLOCATABLE, TARGET   :: FLUXWIND(:,:)

CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCOX_Custom_Run
!
! !DESCRIPTION: Subroutine HCOX\_Custom\_Run is the driver routine 
! for the customizable HEMCO extension. 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCOX_Custom_Run( am_I_Root, ExtState, HcoState, RC )
!
! !USES:
!
    USE HCO_FluxArr_Mod,  ONLY : HCO_EmisAdd
    USE HCO_GeoTools_Mod, ONLY : HCO_LANDTYPE
!
! !INPUT PARAMETERS:
!
    LOGICAL,         INTENT(IN   ) :: am_I_Root   ! Are we on the root CPU?
    TYPE(Ext_State), POINTER       :: ExtState    ! Module options
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(HCO_State), POINTER       :: HcoState    ! Hemco state 
    INTEGER,         INTENT(INOUT) :: RC          ! Success or failure
!
! !REMARKS:
!  
!
! !REVISION HISTORY:
!  13 Dec 2013 - C. Keller   - Initial version 
!  05 Jun 2014 - R. Yantosca - Now store the results of HCO_LANDTYPE
!                              in a PRIVATE variable for the !OMP loop
!  05 Jun 2014 - R. Yantosca - Cosmetic changes
!  06 Jun 2014 - R. Yantosca - Now indented with F90 free-format
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER             :: I, J, N, AS, LANDTYPE
    INTEGER             :: tmpID
    REAL*8              :: W10M
    LOGICAL             :: ERR
!
! !DEFINED PARAMETERS:
!
    REAL*8,   PARAMETER :: SCALICE  = 1.0d-14
    REAL*8,   PARAMETER :: SCALWIND = 1.0d-14

    !=================================================================
    ! HCOX_CUSTOM_RUN begins here!
    !=================================================================

    ! Enter
    CALL HCO_ENTER( HcoState%Config%Err, 'HCOX_Custom_Run (hcox_custom_mod.F90)', RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Set error flag
    ERR = .FALSE.

    ! Sanity check: return if extension not turned on
    IF ( .NOT. ExtState%Custom ) RETURN

    ! Initialize flux arrays
    ALLOCATE ( FLUXICE( HcoState%NX,HcoState%NY),        &
               FLUXWIND(HcoState%NX,HcoState%NY), STAT=AS )
    IF ( AS/= 0 ) THEN
       CALL HCO_ERROR( HcoState%Config%Err, 'ALLOCATION ERROR', RC )
       RETURN
    ENDIF
    FLUXICE  = 0.0_hp
    FLUXWIND = 0.0_hp

!$OMP PARALLEL DO                                            &
!$OMP DEFAULT( SHARED )                                      &
!$OMP PRIVATE( I, J, W10M, LANDTYPE                        ) &
!$OMP SCHEDULE( DYNAMIC )
    ! Loop over surface grid boxes
    DO J = 1, HcoState%NY
    DO I = 1, HcoState%NX

       ! Get the land type for grid box (I,J)
       LANDTYPE = HCO_LANDTYPE( ExtState%WLI%Arr%Val(I,J),  &
                                ExtState%ALBD%Arr%Val(I,J) )

       ! Check surface type
       ! Ocean:
       IF ( LANDTYPE == 0 ) THEN

          ! 10m wind speed [m/s]
          W10M = ExtState%U10M%Arr%Val(I,J)**2 + &
                 ExtState%V10M%Arr%Val(I,J)**2
          W10M = SQRT(W10M)

          ! Set flux to wind speed
          FLUXWIND(I,J) = W10M * SCALWIND

         ! Ice:         
       ELSE IF ( LANDTYPE == 2 ) THEN

          ! Set uniform flux
          FLUXICE(I,J) = SCALICE
       ENDIF

    ENDDO !I
    ENDDO !J
!$OMP END PARALLEL DO

    ! Check exit status
    IF ( ERR ) THEN
       RC = HCO_FAIL
       RETURN 
    ENDIF

    ! Add wind fluxes to emission arrays & diagnostics 
    DO N = 1, nOcWind

       ! Emissions array
       CALL HCO_EmisAdd( am_I_Root, HcoState, FLUXWIND, OcWindIDs(N), &
                         RC,        ExtNr=ExtNr )
       IF ( RC /= HCO_SUCCESS ) RETURN
    ENDDO !N

    ! Add ice fluxes to emission arrays & diagnostics 
    DO N = 1, nIceSrc

       ! Emissions array
       CALL HCO_EmisAdd( am_I_Root, HcoState, FLUXICE, IceSrcIDs(N), &
                         RC,        ExtNr=ExtNr )
       IF ( RC /= HCO_SUCCESS ) RETURN
    ENDDO !N

    ! Return w/ success
    CALL HCO_LEAVE( HcoState%Config%Err,RC )

  END SUBROUTINE HCOX_Custom_Run
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCOX_Custom_Init
!
! !DESCRIPTION: Subroutine HCOX\_Custom\_Init initializes the HEMCO
! CUSTOM extension.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCOX_Custom_Init( am_I_Root, HcoState, ExtName, &
                               ExtState,  RC                  )
!
! !USES:
!
    USE HCO_ExtList_Mod,    ONLY : GetExtNr
    USE HCO_STATE_MOD,      ONLY : HCO_GetExtHcoID
!
! !INPUT PARAMETERS:
!
    LOGICAL,          INTENT(IN   ) :: am_I_Root
    CHARACTER(LEN=*), INTENT(IN   ) :: ExtName    ! Extension name
    TYPE(Ext_State),  POINTER       :: ExtState   ! Module options      
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(HCO_State),  POINTER       :: HcoState   ! Hemco state 
    INTEGER,          INTENT(INOUT) :: RC 

! !REVISION HISTORY:
!  13 Dec 2013 - C. Keller   - Now a HEMCO extension
!  06 Jun 2014 - R. Yantosca - Cosmetic changes in ProTeX headers
!  06 Jun 2014 - R. Yantosca - Now indented using F90 free-format
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER                        :: N, nSpc, AS
    INTEGER,           ALLOCATABLE :: HcoIDs(:)
    LOGICAL                        :: verb
    CHARACTER(LEN=31), ALLOCATABLE :: SpcNames(:)
    CHARACTER(LEN=255)             :: MSG 

    !=================================================================
    ! HCOX_CUSTOM_INIT begins here!
    !=================================================================

    ! Extension Nr.
    ExtNr = GetExtNr( HcoState%Config%ExtList, TRIM(ExtName) )
    IF ( ExtNr <= 0 ) RETURN

    ! Enter
    CALL HCO_ENTER( HcoState%Config%Err, 'HCOX_Custom_Init (hcox_custom_mod.F90)', RC )
    IF ( RC /= HCO_SUCCESS ) RETURN
    verb = HCO_IsVerb(HcoState%Config%Err,1) 

    ! Set species IDs      
    CALL HCO_GetExtHcoID( HcoState, ExtNr, HcoIDs, SpcNames, nSpc, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Assume first half are 'wind species', second half are ice.
    IF ( MOD(nSpc,2) /= 0 ) THEN
       MSG = 'Cannot set species IDs for custom emission module!'
       CALL HCO_ERROR(HcoState%Config%Err,MSG, RC )
       RETURN
    ENDIF
 
    ! Pass # of sources 
    nOcWind = nSpc / 2
    nIceSrc = nSpc / 2

    ! Allocate vector w/ the species IDs 
    ALLOCATE ( OcWindIDs(nOcWind) ) 
    ALLOCATE ( IceSrcIDs(nIceSrc) ) 
    OcWindIDs(:) = HcoIDs(1:nOcWind)
    N = nOcWind + 1
    IceSrcIDs(:) = HcoIDs(N:nSpc)

    ! Allocate flux arrays (to be filled every time step)
    ALLOCATE ( FLUXWIND(HcoState%NX, HcoState%NY), &
               FLUXICE (HcoState%NX, HcoState%NY), STAT=AS )
    IF ( AS/=0 ) THEN
       MSG = 'Allocation error: FLUXWIND, FLUXICE'
       CALL HCO_ERROR(HcoState%Config%Err,MSG, RC )
       RETURN
    ENDIF

    ! Verbose mode
    IF ( verb ) THEN
       MSG = 'Use custom emissions module (extension module)'
       CALL HCO_MSG(HcoState%Config%Err,MSG )

       MSG = 'Use the following species (Name: HcoID):'
       CALL HCO_MSG(HcoState%Config%Err,MSG)
       DO N = 1, nSpc
          WRITE(MSG,*) TRIM(SpcNames(N)), ':', HcoIDs(N)
          CALL HCO_MSG(HcoState%Config%Err,MSG)
       ENDDO
    ENDIF

    ! Activate met fields required by this extension
    ExtState%U10M%DoUse = .TRUE. 
    ExtState%V10M%DoUse = .TRUE. 
    ExtState%ALBD%DoUse = .TRUE. 
    ExtState%WLI%DoUse  = .TRUE. 

    ! Activate this extension
    ExtState%Custom = .TRUE.

    ! Leave w/ success
    IF ( ALLOCATED(HcoIDs  ) ) DEALLOCATE(HcoIDs  )
    IF ( ALLOCATED(SpcNames) ) DEALLOCATE(SpcNames)

    CALL HCO_LEAVE( HcoState%Config%Err,RC ) 

  END SUBROUTINE HCOX_Custom_Init
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCOX_Custom_Final
!
! !DESCRIPTION: Subroutine HCOX\_Custom\_Final finalizes the HEMCO
! CUSTOM extension.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCOX_Custom_Final
!
! !REVISION HISTORY:
!  13 Dec 2013 - C. Keller   - Now a HEMCO extension
!  06 Jun 2014 - R. Yantosca - Cosmetic changes in ProTeX headers
!  06 Jun 2014 - R. Yantosca - Now indended with F90 free-format
!EOP
!------------------------------------------------------------------------------
!BOC

    !=================================================================
    ! HCOX_CUSTOM_FINAL begins here!
    !=================================================================

    IF ( ALLOCATED(OcWindIDs) ) DEALLOCATE ( OcWindIDs )
    IF ( ALLOCATED(IceSrcIDs) ) DEALLOCATE ( IceSrcIDs )
    IF ( ALLOCATED(FLUXICE  ) ) DEALLOCATE ( FLUXICE   )
    IF ( ALLOCATED(FLUXWIND ) ) DEALLOCATE ( FLUXWIND  )

  END SUBROUTINE HCOX_Custom_Final
!EOC
END MODULE HCOX_Custom_Mod
