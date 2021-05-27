!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: gckpp_HetRates
!
! !DESCRIPTION: FlexChem module for aerosol and aqueous chemistry, via KPP.
!\\
!\\
! !INTERFACE:

MODULE GCKPP_AQRATES
!
! !USES:
!
  USE CMN_SIZE_Mod,     ONLY : NDSTBIN
  USE ErrCode_Mod
  USE Error_Mod,        ONLY : IS_SAFE_DIV, SAFE_DIV
  USE GcKpp_Global
  USE GcKpp_Precision
  USE GcKpp_Parameters
  USE GCKPP_HetRates,   ONLY : ARSL1K, KIIR1Ltd, KIIR1R2L
  USE State_Chm_Mod,    ONLY : ChmState
  USE State_Chm_Mod,    ONLY : Ind_
  USE State_Met_Mod,    ONLY : MetState
  USE State_Grid_Mod,   ONLY : GrdState
  USE Input_Opt_Mod,    ONLY : OptInput
  USE PhysConstants
  USE Precision_Mod
  Use Pressure_Mod,     ONLY : Get_PCenter

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC  :: SET_SEASALT_S
  PUBLIC  :: INIT_AQRATES
!
! !PRIVATE MEMBER FUNCTIONS:
!
  ! These functions are used for all mechanisms
!! PRIVATE
!!$OMP THREADPRIVATE( 
!!$OMP THREADPRIVATE( )
! !PRIVATE DATA MEMBERS:
!
  ! Scalars
  INTEGER, SAVE :: id_SALAAL, id_SALCAL, id_SO2, id_O3, id_HCl, id_HNO3
! !REMARKS:
!
! !REVISION HISTORY:
!  24 Mar 2021 - M. Long - Initial Version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
  CONTAINS

    SUBROUTINE INIT_AQRATES( RC )

      CHARACTER(LEN=255)            :: ErrMsg,   ThisLoc
      INTEGER,        INTENT(OUT)   :: RC         ! Success or failure

      ! Initialize
      RC      = GC_SUCCESS
      errMsg  = ''
      thisLoc = ' -> at INIT_AQRATES (in module KPP/fullchem/gckpp_AqRates.F90)'
      
      id_SALAAL = IND_('SALAAL')
      id_SALCAL = IND_('SALCAL')

       IF ( id_SALAAL < 0 ) THEN
          errMsg = 'SALAAL is not a defined species in this simulation!!!'
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF
       IF ( id_SALCAL < 0 ) THEN
          errMsg = 'SALCAL is not a defined species in this simulation!!!'
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

      id_SO2  = IND_('SO2')
      id_O3   = IND_('O3')
      id_HCl  = IND_('HCl')
      id_HNO3 = IND_('HNO3')

       IF ( id_SALAAL < 0 ) THEN
          errMsg = 'SALAAL is not a defined species in this simulation!!!'
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

    END SUBROUTINE INIT_AQRATES

!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Set_Aq
!
! !DESCRIPTION: Main aqueous/aerosol chemistry driver routine.  Sets up the
!  vector of aqueous chemistry rates for the KPP chemistry solver.
!\\
!\\
! !INTERFACE:
!
    SUBROUTINE SET_SEASALT_S( I, J, L, Input_Opt, State_Chm, State_Grid, &
         State_Met, RC )
!
! !USES:
!
      USE GcKpp_Global
!
! !INPUT PARAMETERS:
!
      INTEGER,        INTENT(IN)    :: I, J, L    ! Lon, lat, level indices
      TYPE(MetState), INTENT(IN)    :: State_Met  ! Meteorology State object
      TYPE(GrdState), INTENT(IN)    :: State_Grid  ! Grid State object
      TYPE(OptInput), INTENT(IN)    :: Input_Opt  ! Input Options object
!
! !INPUT/OUTPUT PARAMETERS:
!
      TYPE(ChmState), INTENT(INOUT) :: State_Chm  ! Chemistry State object
!
! !REMARKS:
!
! ! Reaction List (by K_MT() index)
! 1) SO2 + O3 + 2SALAAL --> SO4mm + O2 : From Sulfate_mod - 24 Mar 2021
!
! !REVISION HISTORY:
!  24 Mar 2021 - M. Long - Initial Version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      REAL(fp)                      :: ALK1, ALK2, KT1, KT2, KT1N, KT2N, KT1L, KT2L ! GET_ALK()
      REAL(fp)                      :: gamma, cvd
      REAL(fp)                      :: k_ex, k_ex2 ! Gas exchange rate. 1-forward. 2-backward

      REAL(fp)              :: ALK_d   (NDSTBIN)
      REAL(fp)              :: ALKA_d  (NDSTBIN)
      REAL(fp)              :: PSO4_d  (NDSTBIN)
      REAL(fp)              :: PNIT_d  (NDSTBIN)
      REAL(fp)              :: PH2SO4_d(NDSTBIN)
      REAL(fp)              :: KTN(NDSTBIN)
      REAL(fp)              :: KTS(NDSTBIN)
      REAL(fp)              :: KTH(NDSTBIN)

      CHARACTER(LEN=255)            :: ErrMsg,   ThisLoc
      INTEGER                       :: N
      INTEGER,        INTENT(OUT)   :: RC         ! Success or failure

      ! Initialize
      RC      = GC_SUCCESS
      errMsg  = ''
      thisLoc = ' -> at SET_SEASALT_SRATES (in module KPP/fullchem/gckpp_AqRates.F90)'

      !====================================================================
      ! SET_SEASALT_S begins here!
      !====================================================================
      K_MT = 0._dp

      IF (State_Chm%Species(I,J,L,id_SALAAL) .gt. 1.e-1_fp) then
      ! SO2 + O3 ...
         k_ex  = ARSL1K( State_Chm%WetAeroArea(I,J,L,11), &
              STATE_HET%XRADI(11), NUMDEN, 0.11_dp, SR_TEMP, SR_MW(ind_SO2))
         K_MT(1) = kIIR1Ltd(C(ind_SO2), C(ind_O3), k_ex)/State_Chm%Species(I,J,L,id_SALAAL)!**2
      endif

      ! HCl
      k_ex = ARSL1K( State_Chm%WetAeroArea(I,J,L,11), &
           STATE_HET%XRADI(11), NUMDEN, 0.074_dp, SR_TEMP, SR_MW(ind_HCl))

      IF (State_Chm%Species(I,J,L,id_SALAAL) .gt. 1.e-1_fp) &
         K_MT(2) = k_ex!/State_Chm%Species(I,J,L,id_SALAAL)

      ! HNO3
      k_ex = ARSL1K( State_Chm%WetAeroArea(I,J,L,11), &
           STATE_HET%XRADI(11), NUMDEN, 0.5_dp, SR_TEMP, SR_MW(ind_HNO3))

      IF (State_Chm%Species(I,J,L,id_SALAAL) .gt. 1.e-1_fp) &
         K_MT(3) = k_ex!/State_Chm%Species(I,J,L,id_SALAAL)

      ! End fine seasalt

      ! Coarse seasalt
      IF (State_Chm%Species(I,J,L,id_SALCAL) .gt. 1.e-1_fp) then
         ! SO2 + O3 ...
         k_ex  = ARSL1K( State_Chm%WetAeroArea(I,J,L,12), &
              STATE_HET%XRADI(12), NUMDEN, 0.11_dp, SR_TEMP, SR_MW(ind_SO2))
         K_MT(4) = kIIR1Ltd(C(ind_SO2), C(ind_O3), k_ex)/State_Chm%Species(I,J,L,id_SALCAL)!**2
      endif
      
      ! HCl
      k_ex = ARSL1K( State_Chm%WetAeroArea(I,J,L,12), &
           STATE_HET%XRADI(12), NUMDEN, 0.074_dp, SR_TEMP, SR_MW(ind_HCl))

      IF (State_Chm%Species(I,J,L,id_SALCAL) .gt. 1.e-1_fp) &
         K_MT(5) = k_ex!/State_Chm%Species(I,J,L,id_SALCAL)
      
      ! HNO3
      k_ex = ARSL1K( State_Chm%WetAeroArea(I,J,L,12), &
           STATE_HET%XRADI(12), NUMDEN, 0.5_dp, SR_TEMP, SR_MW(ind_HNO3))

      IF (State_Chm%Species(I,J,L,id_SALCAL) .gt. 1.e-1_fp) &
         K_MT(6) = k_ex!/State_Chm%Species(I,J,L,id_SALCAL)

      ! End coarse seasalt
      
    END SUBROUTINE SET_SEASALT_S
!EOC 
END MODULE GCKPP_AQRATES
