!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!     
! !MODULE: get_ndep_mod.F
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
      REAL(fp), PUBLIC, POINTER   :: DRY_TOTN (:,:) => NULL() ! Drydep'd N
      REAL(fp), PUBLIC, POINTER   :: WET_TOTN (:,:) => NULL() ! Wetdep'd N
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
!  16 Jun 2016 - C. Miller   - Now define species ID flags with Ind_
!  17 Jun 2016 - R. Yantosca - Add species ID flags as module variables
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
!
! !LOCAL VARAIABLES:
!     
      ! Species ID flags (formerly in tracerid_mod.F)
      INTEGER             :: id_HNO3
      INTEGER             :: id_NH3
      INTEGER             :: id_NH4
      INTEGER             :: id_NH4aq
      INTEGER             :: id_NIT
      INTEGER             :: id_NITs
      INTEGER             :: id_NO2
      INTEGER             :: id_PAN

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
      USE State_Chm_Mod,      only : Ind_
!
! !INPUT PARAMETERS: 
!
      INTEGER, INTENT(IN)   :: I          ! I
      INTEGER, INTENT(IN)   :: J          ! J
      INTEGER, INTENT(IN)   :: L          ! Level
      INTEGER, INTENT(IN)   :: NN         ! Dry Dep Tracer # 
      REAL(fp), INTENT(IN)  :: TDRYFX     ! Dry dep flux [molec/cm2/s] 
!
! !REVISION HISTORY:
!  23 Oct 2012 - M. Payer    - Added ProTeX headers
!  14 Mar 2013 - M. Payer    - Replace NOx with NO2 as part of removal of
!                              NOx-Ox partitioning
!  16 Jun 2016 - C. Miller   - Now use 

!EOP
!------------------------------------------------------------------------------
!BOC
!
      ! Update the reservoir if it's a nitrogen species
      IF ( NN == id_NO2  .OR. NN == id_PAN     .OR. 
     &     NN == id_HNO3 .OR. NN == id_NH3     .OR.
     &     NN == id_NH4  .OR. NN == id_NH4aq   .OR.
     &     NN == id_NIT  .OR. NN == id_NITs  ) THEN
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
      USE State_Chm_Mod,      Only : Ind_
!
! !INPUT PARAMETERS: 
!
      INTEGER,  INTENT(IN) :: I          ! I
      INTEGER,  INTENT(IN) :: J          ! J
      INTEGER,  INTENT(IN) :: L          ! Level
      INTEGER,  INTENT(IN) :: NN         ! Wet Dep Tracer # 
      REAL(fp), INTENT(IN) :: TWETFX     ! Wet dep flux [kg/s] 
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
      IF ( NN == id_HNO3  ) SCAL = CONVHNO3 
      IF ( NN == id_NH3   ) SCAL = CONVNH3
      IF ( NN == id_NH4   ) SCAL = CONVNH4
      IF ( NN == id_NH4aq ) SCAL = CONVNH4
      IF ( NN == id_NIT   ) SCAL = CONVNIT
      IF ( NN == id_NITs  ) SCAL = CONVNIT

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
      USE ErrCode_Mod
      USE State_Chm_Mod,      ONLY : Ind_
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

      ! Define species ID flags
      id_HNO3  = Ind_('HNO3' )
      id_NH3   = Ind_('NH3'  )
      id_NH4   = Ind_('NH4'  )
      id_NH4aq = Ind_('NH4aq')
      id_NIT   = Ind_('NIT'  )
      id_NITs  = Ind_('NITs' )
      id_NO2   = Ind_('NO2'  )
      id_PAN   = Ind_('PAN'  )

      ! Allocate arrays
      ALLOCATE( WET_TOTN     ( IIPAR, JJPAR ), STAT=RC )
      IF ( RC /= GC_SUCCESS ) RETURN

      ALLOCATE( DRY_TOTN     ( IIPAR, JJPAR ), STAT=RC )
      IF ( RC /= GC_SUCCESS ) RETURN

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
      USE ErrCode_Mod
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
      RC = GC_SUCCESS
      
      ! Deallocate arrays
      IF ( ASSOCIATED( DRY_TOTN      ) ) DEALLOCATE( DRY_TOTN      )
      IF ( ASSOCIATED( WET_TOTN      ) ) DEALLOCATE( WET_TOTN      )

      END SUBROUTINE Cleanup_Get_Ndep
!EOC
      END MODULE GET_NDEP_MOD
