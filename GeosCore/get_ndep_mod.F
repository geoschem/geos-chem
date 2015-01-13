!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!     
! !MODULE: get_ndep_mod
!     
! !DESCRIPTION: Module GET\_NDEP\_MOD contains routines for computing the
! accumulated nitrogen dry and wet deposition between emission time steps. 
! These variables are needed for soil NOx emission calculations.
!\\
!\\
! This module is basically a simple wrapper module to save out the nitrogen
! dry and wet deposition rates and pass them to HEMCO for soil NOx emission 
! calculation (via hcoi\_gc\_main\_mod.F90).
!\\
!\\
! IMPORTANT: Routine RESET\_DEP\_N resets the deposition arrays to zero. It
! is called in hcoi\_gc\_main\_mod.F90 after the emission calculations.
!\\
!\\ 
! !INTERFACE: 
!     
      MODULE GET_NDEP_MOD
!
! !USES:
!
      USE PRECISION_MOD    ! For GEOS-Chem Precision (fp)

      IMPLICIT NONE
      PRIVATE
!
! !PUBLIC DATA MEMBERS:
!
      ! Soil NOx deposited N arrays
      REAL(fp), PUBLIC, ALLOCATABLE, TARGET :: DRY_TOTN     (:,:) ! Drydep'd N
      REAL(fp), PUBLIC, ALLOCATABLE, TARGET :: WET_TOTN     (:,:) ! Wetdep'd N
! 
! !PUBLIC MEMBER FUNCTIONS:
!      
      PUBLIC  :: SOIL_DRYDEP
      PUBLIC  :: SOIL_WETDEP
      PUBLIC  :: RESET_DEP_N
      PUBLIC  :: Init_Get_Ndep
      PUBLIC  :: Cleanup_Get_Ndep
!
! !REVISION HISTORY:
!  23 Oct 2012 - M. Payer    - Added ProTeX headers
!  25 Jul 2014 - R. Yantosca - Moved module variables DEP_RESERVOIR, DRY_TOTN,
!                              and WET_TOTN here from Headers/commsoil_mod.F90
!  25 Jul 2014 - R. Yantosca - Add routine INIT_GET_NDEP
!  25 Jul 2014 - R. Yantosca - Add routine CLEANUP_GET_NDEP
!  09 Oct 2014 - C. Keller   - Removed obsolete routines DEP_RESERVOIR, 
!                              GET_DEP_N, SOURCE_DRYN, and SOURCE_WETN. These
!                              are now handled in hcox_soilnox_mod.F90.
!  13 Nov 2014 - M. Yannetti - Added PRECISION_MOD
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!  
      REAL(fp), PARAMETER :: CONVHNO3 = 0.222 ! MWN/MWHNO3
      REAL(fp), PARAMETER :: CONVNH4  = 0.777 ! MWN/MWNH4
      REAL(fp), PARAMETER :: CONVNH3  = 0.823 ! MWN/MWNH3
      REAL(fp), PARAMETER :: CONVNIT  = 0.226 ! MWN/MWNIT

      CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: soil_drydep
!
! !DESCRIPTION: Subroutine SOIL\_DRYDEP holds dry deposited species
!               [molec/cm2/s]. This is called from dry\_dep\_mod.F.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE SOIL_DRYDEP( I, J, L, NN, TDRYFX )
!
! !USES:
!
      USE TRACERID_MOD, ONLY : IDTNH4,  IDTNIT, IDTNH4aq 
      USE TRACERID_MOD, ONLY : IDTHNO3, IDTNH3, IDTNITs
      USE TRACERID_MOD, ONLY : IDTNO2,  IDTPAN
!
! !INPUT PARAMETERS: 
!
      INTEGER, INTENT(IN)   :: I          ! I
      INTEGER, INTENT(IN)   :: J          ! J
      INTEGER, INTENT(IN)   :: L          ! Level
      INTEGER, INTENT(IN)   :: NN         ! Dry Dep Tracer # 
      REAL(fp), INTENT(IN)   :: TDRYFX     ! Dry dep flux [molec/cm2/s] 
!
! !REVISION HISTORY:
!  23 Oct 2012 - M. Payer    - Added ProTeX headers
!  14 Mar 2013 - M. Payer    - Replace NOx with NO2 as part of removal of
!                              NOx-Ox partitioning
!EOP
!------------------------------------------------------------------------------
!BOC
!

      ! ckeller (14/04/02)
      IF ( NN == IDTNO2 .OR. NN == IDTPAN .OR. NN == IDTHNO3  .OR.
     &     NN == IDTNH3 .OR. NN == IDTNH4 .OR. NN == IDTNH4aq .OR. 
     &     NN == IDTNIT .OR. NN == IDTNITs                        ) THEN
         DRY_TOTN(I,J) = DRY_TOTN(I,J) + TDRYFX
      ENDIF

      END SUBROUTINE SOIL_DRYDEP
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: soil_wetdep
!
! !DESCRIPTION: Subroutine SOIL\_WETDEP holds wet deposited species
!               [molec/cm2/s]. This is called from wetscav\_mod.F.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE SOIL_WETDEP( I, J, L, NN, TWETFX )
!
! !USES:
!
      USE TRACERID_MOD, ONLY : IDTNH4,    IDTNIT, IDTNH4aq 
      USE TRACERID_MOD, ONLY : IDTHNO3,   IDTNH3, IDTNITs
!
! !INPUT PARAMETERS: 
!
      INTEGER, INTENT(IN)   :: I          ! I
      INTEGER, INTENT(IN)   :: J          ! J
      INTEGER, INTENT(IN)   :: L          ! Level
      INTEGER, INTENT(IN)   :: NN         ! Wet Dep Tracer # 
      REAL(fp) , INTENT(IN)   :: TWETFX     ! Wet dep flux [kg/s] 
!
! !REVISION HISTORY:
!  23 Oct 2012 - M. Payer    - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
      ! ckeller (14/04/02)
      REAL(fp)        :: SCAL

      ! ckeller (14/04/02)
      SCAL = 0e+0_fp
      IF ( NN == IDTHNO3  ) SCAL = CONVHNO3 
      IF ( NN == IDTNH3   ) SCAL = CONVNH3
      IF ( NN == IDTNH4   ) SCAL = CONVNH4
      IF ( NN == IDTNH4aq ) SCAL = CONVNH4
      IF ( NN == IDTNIT   ) SCAL = CONVNIT
      IF ( NN == IDTNITs  ) SCAL = CONVNIT

      IF ( SCAL > 0e+0_fp ) THEN
         WET_TOTN(I,J) = WET_TOTN(I,J) + ( TWETFX * SCAL )
      ENDIF

      END SUBROUTINE SOIL_WETDEP
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: reset_dep_N
!
! !DESCRIPTION: Subroutine RESET\_DEP\_N resets the dry and wet deposition 
!               arrays and variables so that they can be refilled. 
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE RESET_DEP_N()
!
! !REVISION HISTORY:
!  03 Apr 2014 - C. Keller   - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC

      ! Reset all variables
      DRY_TOTN = 0e+0_fp
      WET_TOTN = 0e+0_fp

      END SUBROUTINE RESET_DEP_N
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Init_Get_Ndep
!
! !DESCRIPTION: Routine INIT\_GET\_NDEP allocates all module arrays.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE Init_Get_Ndep( am_I_Root, RC )
!
! !USES:
!
      USE CMN_SIZE_Mod
      USE GIGC_ErrCode_Mod
!
! !INPUT PARAMETERS: 
!
      LOGICAL, INTENT(IN)  :: am_I_Root   ! Are we on the root CPU?
!
! !OUTPUT PARAMETERS:
!
      INTEGER, INTENT(OUT) :: RC          ! Success or failure
!
! !REMARKS:
!  We now allocate the DEP_RESERVOIR, DRY_TOTN, and WET_TOTN arrays within
!  get_ndep_mod.F.  These were formerly contained in Headers/commsoil_mod.F90,
!  which was rendered mostly obsolete by HEMCO.
! 
! !REVISION HISTORY: 
!  25 Jul 2014 - R. Yantosca - Initial version
!  09 Oct 2014 - C. Keller   - Removed obsolete variable DEP_RESERVOIR
!EOP
!------------------------------------------------------------------------------
!BOC

      ! Allocate arrays
      ALLOCATE( WET_TOTN     ( IIPAR, JJPAR ), STAT=RC )
      IF ( RC /= GIGC_SUCCESS ) RETURN

      ALLOCATE( DRY_TOTN     ( IIPAR, JJPAR ), STAT=RC )
      IF ( RC /= GIGC_SUCCESS ) RETURN

      ! Zero arrays
      WET_TOTN      = 0e+0_fp
      DRY_TOTN      = 0e+0_fp

      END SUBROUTINE Init_Get_Ndep
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: cleanup_commsoil
!
! !DESCRIPTION: Subroutine CLEANUP\_COMMSOIL deallocates all module arrays.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE Cleanup_Get_NDep( am_I_Root, RC )
!
! !USES:
!
      USE GIGC_ErrCode_Mod
!
! !INPUT PARAMETERS: 
!     
      LOGICAL, INTENT(IN)  :: am_I_Root   ! Are we on the root CPU?
!
! !OUTPUT PARAMETERS:
!
      INTEGER, INTENT(OUT) :: RC          ! Return code
! 
! !REVISION HISTORY: 
!  25 Jul 2014 - R. Yantosca - Initial version
!  09 Oct 2014 - C. Keller   - Removed obsolete variable DEP_RESERVOIR
!EOP
!------------------------------------------------------------------------------
!BOC
      ! Assume success
      RC = GIGC_SUCCESS
      
      ! Deallocate arrays
      IF ( ALLOCATED( DRY_TOTN      ) ) DEALLOCATE( DRY_TOTN      )
      IF ( ALLOCATED( WET_TOTN      ) ) DEALLOCATE( WET_TOTN      )

      END SUBROUTINE Cleanup_Get_Ndep
!EOC
      END MODULE GET_NDEP_MOD
