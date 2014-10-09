#if defined( DEVEL ) || defined( EXTERNAL_GRID ) || defined( EXTERNAL_FORCING )
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: gigc_state_phy_mod
!
! !DESCRIPTION: Module GIGC\_STATE\_PHY\_MOD contains the derived type
!  used to define the Physics State object for the Grid-Independent 
!  GEOS-Chem implementation (abbreviated "GIGC").
!\\
!\\
!  This module also contains the routines that allocate and deallocate memory 
!  to the Physics State object.  The chemistry state object is not defined
!  in this module.  It must be be declared as variable in the top-level 
!  driver routine, and then passed to lower-level routines as an argument.
!\\
!\\
! !INTERFACE: 
!
MODULE GIGC_State_Phy_Mod
!
! USES:
!
  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC :: Init_GIGC_State_Phy
  PUBLIC :: Cleanup_GIGC_State_Phy
!
! !PUBLIC DATA MEMBERS:
!
  !=========================================================================
  ! Derived type for physics state
  !=========================================================================
  TYPE, PUBLIC :: PhyState

     ! Scalars
     INTEGER         :: BEGIN_I
     INTEGER         :: END_I
     INTEGER         :: BEGIN_J
     INTEGER         :: END_J
     INTEGER         :: NINDX
     INTEGER         :: NJNDX

     ! 1-D arrays
     REAL*8, POINTER :: LAT     (:  )  ! LATITUDE (DEG)
     REAL*8, POINTER :: LON     (:  )  ! LONGITUDE (DEG)
     REAL*8, POINTER :: PS      (:  )  ! SURFACE PRESSURE
     REAL*8, POINTER :: PSDRY   (:  )  ! DRY SURFACE PRESSURE
     REAL*8, POINTER :: PHIS    (:  )  ! SURFACE GEOPOTENTIAL
     REAL*8, POINTER :: ULAT    (:  )  ! UNIQUE LATITUDES  (DEG)
     REAL*8, POINTER :: ULON    (:  )  ! UNIQUE LONGITUDES (DEG)

     ! 2-D arrays
     REAL*8, POINTER :: T       (:,:)   ! TEMPERATURE (K)
     REAL*8, POINTER :: U       (:,:)   ! ZONAL WIND (M/S)
     REAL*8, POINTER :: V       (:,:)   ! MERIDIONAL WIND (M/S)
     REAL*8, POINTER :: OMEGA   (:,:)   ! VERTICAL PRESSURE VELOCITY (PA/S) 
     REAL*8, POINTER :: PMID    (:,:)   ! MIDPOINT PRESSURE (PA) 
     REAL*8, POINTER :: PMIDDRY (:,:)   ! MIDPOINT PRESSURE DRY (PA) 
     REAL*8, POINTER :: PDEL    (:,:)   ! LAYER THICKNESS (PA)
     REAL*8, POINTER :: PDELDRY (:,:)   ! LAYER THICKNESS DRY (PA)
     REAL*8, POINTER :: RPDEL   (:,:)   ! 1/LAYER THICKNESS (PA)
     REAL*8, POINTER :: RPDELDRY(:,:)   ! 1/LAYER THICKNESS DRY (PA)
     REAL*8, POINTER :: UZM     (:,:)   ! ZONAL WIND FOR QBO (M/S)
     REAL*8, POINTER :: ZM      (:,:)   ! GEOPOTENTIAL HEIGHT ABOVE 
                                            !  SURFACE AT MIDPOINTS (M)

! Leave commented out for now         
!         REAL*8, DIMENSION(PCOLS,PVER+1)           :: &
!              PINT,    &  ! INTERFACE PRESSURE (PA)
!              PINTDRY, &  ! INTERFACE PRESSURE DRY (PA) 
!              LNPINT,  &  ! LN(PINT)
!              LNPINTDRY,& ! LOG INTERFACE PRESSURE DRY (PA) 
!              ZI          ! GEOPOTENTIAL HEIGHT ABOVE SURFACE AT INTERFACES (M)
!         
!         INTEGER, DIMENSION(PCOLS) :: &
!              LATMAPBACK, &! MAP FROM COLUMN TO UNIQUE LAT FOR THAT COLUMN
!              LONMAPBACK, &! MAP FROM COLUMN TO UNIQUE LON FOR THAT COLUMN
!              CID         ! UNIQUE COLUMN ID
!         INTEGER :: ULATCNT, & ! NUMBER OF UNIQUE LATS IN CHUNK
!              ULONCNT     ! NUMBER OF UNIQUE LONS IN CHUNK

  END TYPE PhyState

  !=========================================================================
  ! Other variables
  !=========================================================================

  ! Position value used for registering CSPEC parameters in the chemical state
  INTEGER,  SAVE :: POSITION = 1 
!
! !REMARKS:
!  This may be just for the BCC interface.
!                                                                             
! !REVISION HISTORY:
!  15 Oct 2012 - M. Long     - Initial version, based on gc_type_mod.F
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
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
! !IROUTINE: init_gigc_state_phy
!
! !DESCRIPTION: Subroutine INIT\_GIGC\_STATE\_PHY allocates all fields of 
!  the Grid-Indpendent GEOS-Chem (aka "GIGC") Physics State object.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Init_GIGC_State_Phy( am_I_Root, IM, JM, LM, State_Phy, RC )
!
! !USES:
!
    USE GIGC_ErrCode_Mod                         ! Error codes
!
! !INPUT PARAMETERS:
! 
    LOGICAL,        INTENT(IN)    :: am_I_Root   ! Is this the root CPU?
    INTEGER,        INTENT(IN)    :: IM          ! # longitudes on this PET
    INTEGER,        INTENT(IN)    :: JM          ! # longitudes on this PET
    INTEGER,        INTENT(IN)    :: LM          ! # longitudes on this PET
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(PhyState), INTENT(INOUT) :: State_Phy   ! Physics State object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Return code
!
! !REVISION HISTORY: 
!  19 Oct 2012 - R. Yantosca - Initial version, based on gc_environment_mod.F90
!  19 Oct 2012 - R. Yantosca - Now pass all dimensions as arguments
!EOP
!------------------------------------------------------------------------------
!BOC

    ! Assume success
    RC = GIGC_SUCCESS
!
!    ! 1-D arrays
!   ALLOCATE( State_Met%LAT     (JM  )  
!   IF ( RC /= GIGC_SUCCESS ) RETURN
!
!   ALLOCATE( State_Met%LON     (IM  )  
!   IF ( RC /= GIGC_SUCCESS ) RETURN
!
!   ALLOCATE( State_Met%PS      (:  )  
!   IF ( RC /= GIGC_SUCCESS ) RETURN
!
!   ALLOCATE( State_Met%PSDRY   (:  )  
!   IF ( RC /= GIGC_SUCCESS ) RETURN
!
!   ALLOCATE( State_Met%PHIS    (:  )  
!   IF ( RC /= GIGC_SUCCESS ) RETURN
!
!   ALLOCATE( State_Met%ULAT    (:  )  
!   IF ( RC /= GIGC_SUCCESS ) RETURN
!
!   ALLOCATE( State_Met%ULON    (:  )  
!   IF ( RC /= GIGC_SUCCESS ) RETURN
!
!    ! 2-D arrays
!   ALLOCATE( State_Met%T       (:,:)   
!   IF ( RC /= GIGC_SUCCESS ) RETURN
!
!   ALLOCATE( State_Met%U       (:,:)   
!   IF ( RC /= GIGC_SUCCESS ) RETURN
!
!   ALLOCATE( State_Met%V       (:,:)   
!   IF ( RC /= GIGC_SUCCESS ) RETURN
!
!   ALLOCATE( State_Met%OMEGA   (:,:)   
!   IF ( RC /= GIGC_SUCCESS ) RETURN
!
!   ALLOCATE( State_Met%PMID    (:,:)   
!   IF ( RC /= GIGC_SUCCESS ) RETURN
!
!   ALLOCATE( State_Met%PMIDDRY (:,:)   
!   IF ( RC /= GIGC_SUCCESS ) RETURN
!
!   ALLOCATE( State_Met%PDEL    (:,:)   
!   IF ( RC /= GIGC_SUCCESS ) RETURN
!
!   ALLOCATE( State_Met%PDELDRY (:,:)   
!   IF ( RC /= GIGC_SUCCESS ) RETURN
!
!   ALLOCATE( State_Met%RPDEL   (:,:)   
!   IF ( RC /= GIGC_SUCCESS ) RETURN
!
!   ALLOCATE( State_Met%RPDELDRY(:,:)   
!   IF ( RC /= GIGC_SUCCESS ) RETURN
!
!   ALLOCATE( State_Met%UZM     (:,:)   
!   IF ( RC /= GIGC_SUCCESS ) RETURN
!
!   ALLOCATE( State_Met%ZM      (:,:)   
!   IF ( RC /= GIGC_SUCCESS ) RETURN
                                      
  END SUBROUTINE Init_GIGC_State_Phy
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: cleanup_gigc_state_phy
!
! !DESCRIPTION: Subroutine CLEANUP\_GIGC\_STATE\_PHY allocates all fields 
!  of the Grid-Independent GEOS-Chem (aka "GIGC") Physics State object.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Cleanup_GIGC_State_Phy( am_I_Root, State_Phy, RC )
!
! !USES:
!
    USE GIGC_ErrCode_Mod                         ! Error codes
!
! !INPUT PARAMETERS:
! 
    LOGICAL,        INTENT(IN)    :: am_I_Root   ! Is this the root CPU?
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(PhyState), INTENT(INOUT) :: State_Phy   ! Physics State object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Return code
!
! !REMARKS:
!  For consistency, maybe this should be moved to a different module.
!
! !REVISION HISTORY: 
!  19 Oct 2012 - R. Yantosca - Initial version, based on gc_environment_mod.F90
!EOP
!------------------------------------------------------------------------------
!BOC
    ! Assume success
    RC = GIGC_SUCCESS

    ! Deallocate fields
    IF ( ASSOCIATED( State_Phy%LAT      )  ) DEALLOCATE( State_Phy%LAT      ) 
    IF ( ASSOCIATED( State_Phy%LON      )  ) DEALLOCATE( State_Phy%LON      ) 
    IF ( ASSOCIATED( State_Phy%PS       )  ) DEALLOCATE( State_Phy%PS       )
    IF ( ASSOCIATED( State_Phy%PSDRY    )  ) DEALLOCATE( State_Phy%PSDRY    )
    IF ( ASSOCIATED( State_Phy%PHIS     )  ) DEALLOCATE( State_Phy%PHIS     )
    IF ( ASSOCIATED( State_Phy%ULAT     )  ) DEALLOCATE( State_Phy%ULAT     )
    IF ( ASSOCIATED( State_Phy%ULON     )  ) DEALLOCATE( State_Phy%ULON     )
    IF ( ASSOCIATED( State_Phy%T        )  ) DEALLOCATE( State_Phy%T        )
    IF ( ASSOCIATED( State_Phy%U        )  ) DEALLOCATE( State_Phy%U        )
    IF ( ASSOCIATED( State_Phy%V        )  ) DEALLOCATE( State_Phy%V        )
    IF ( ASSOCIATED( State_Phy%OMEGA    )  ) DEALLOCATE( State_Phy%OMEGA    )
    IF ( ASSOCIATED( State_Phy%PMID     )  ) DEALLOCATE( State_Phy%PMID     )
    IF ( ASSOCIATED( State_Phy%PMIDDRY  )  ) DEALLOCATE( State_Phy%PMIDDRY  )
    IF ( ASSOCIATED( State_Phy%RPDEL    )  ) DEALLOCATE( State_Phy%RPDEL    )
    IF ( ASSOCIATED( State_Phy%RPDELDRY )  ) DEALLOCATE( State_Phy%RPDELDRY )
    IF ( ASSOCIATED( State_Phy%UZM      )  ) DEALLOCATE( State_Phy%UZM      )
    IF ( ASSOCIATED( State_Phy%ZM       )  ) DEALLOCATE( State_Phy%ZM       )

  END SUBROUTINE Cleanup_GIGC_State_Phy
!EOC
END MODULE GIGC_State_Phy_Mod
#endif
