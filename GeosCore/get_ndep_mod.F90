!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: get_ndep_mod.F90
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
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC  :: SOIL_DRYDEP
  PUBLIC  :: SOIL_WETDEP
  PUBLIC  :: RESET_DEP_N
  PUBLIC  :: Init_Get_Ndep
!
! !REVISION HISTORY:
!  See https://github.com/geoschem/geos-chem for complete history
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
  SUBROUTINE SOIL_DRYDEP( I, J, L, NN, TDRYFX, State_Chm )
!
! !USES:
!
    USE State_Chm_Mod, ONLY : ChmState
!
! !INPUT PARAMETERS:
!
    INTEGER,  INTENT(IN)  :: I          ! I
    INTEGER,  INTENT(IN)  :: J          ! J
    INTEGER,  INTENT(IN)  :: L          ! Level
    INTEGER,  INTENT(IN)  :: NN         ! Dry Dep Tracer #
    REAL(fp), INTENT(IN)  :: TDRYFX     ! Dry dep flux [molec/cm2/s]
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
!
! !REVISION HISTORY:
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
    ! Update the reservoir if it's a nitrogen species
    IF ( NN == id_NO2  .OR. NN == id_PAN     .OR. &
         NN == id_HNO3 .OR. NN == id_NH3     .OR. &
         NN == id_NH4  .OR. NN == id_NH4aq   .OR. &
         NN == id_NIT  .OR. NN == id_NITs  ) THEN
       State_Chm%DryDepNitrogen(I,J) = State_Chm%DryDepNitrogen(I,J) + TDRYFX
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
  SUBROUTINE SOIL_WETDEP( I, J, L, NN, TWETFX, State_Chm )
!
! !USES:
!
    USE State_Chm_Mod, ONLY : ChmState
!
! !INPUT PARAMETERS:
!
    INTEGER,  INTENT(IN) :: I          ! I
    INTEGER,  INTENT(IN) :: J          ! J
    INTEGER,  INTENT(IN) :: L          ! Level
    INTEGER,  INTENT(IN) :: NN         ! Wet Dep Tracer #
    REAL(fp), INTENT(IN) :: TWETFX     ! Wet dep flux [kg/s]
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
!
! !REVISION HISTORY:
!  See https://github.com/geoschem/geos-chem for complete history
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
       State_Chm%WetDepNitrogen(I,J) = State_Chm%WetDepNitrogen(I,J) + &
                                       ( TWETFX * SCAL )
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
  SUBROUTINE RESET_DEP_N( State_Chm )
!
! !USES:
!
    USE State_Chm_Mod,      ONLY : ChmState
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
!
! !REVISION HISTORY:
!  03 Apr 2014 - C. Keller   - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC

    ! Reset all variables
    State_Chm%DryDepNitrogen = 0e+0_fp
    State_Chm%WetDepNitrogen = 0e+0_fp

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
  SUBROUTINE Init_Get_Ndep( Input_Opt, State_Chm, State_Diag, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE State_Chm_Mod,      ONLY : Ind_
    USE Input_Opt_Mod,      ONLY : OptInput
    USE State_Chm_Mod,      ONLY : ChmState
    USE State_Diag_Mod,     ONLY : DgnState
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
    TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
    TYPE(DgnState), INTENT(INOUT) :: State_Diag  ! Diagnostics State object
!
! !OUTPUT PARAMETERS:
!
    INTEGER, INTENT(OUT) :: RC          ! Success or failure
!
! !REVISION HISTORY:
!  25 Jul 2014 - R. Yantosca - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
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

  END SUBROUTINE Init_Get_Ndep
!EOC
END MODULE GET_NDEP_MOD
