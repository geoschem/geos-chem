!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: sfcvmr_mod.F90
!
! !DESCRIPTION: Module sfcvmr\_mod.F90 is a simple module which forces 
!  surface concentrations of relevant species to pre-determined values.
!\\
!\\
! !INTERFACE:
!
MODULE SFCVMR_MOD
!
! !USES:
!
  USE PRECISION_MOD

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC :: fixSfcVMR
!
! !REVISION HISTORY:
!  24 Dec 2016 - S. D. Eastham - Initial version. 
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !PRIVATE TYPES: 
!

CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: FIXSFCVMR
!
! !DESCRIPTION: Subroutine FIXSFCVMR fixes the VMR of selected species
! throughout the PBL to observed values.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE fixSfcVMR( am_I_Root, Input_Opt, State_Met, State_Chm, RC ) 
!
! !USES:
!
    USE ErrCode_Mod
    USE CMN_SIZE_MOD
    USE Input_Opt_Mod,      ONLY : OptInput
    USE State_Met_Mod,      ONLY : MetState
    USE State_Chm_Mod,      ONLY : ChmState
    USE State_Chm_Mod,      ONLY : Ind_

    ! Needed for the new CHxCly boundary condition
    USE PBL_MIX_MOD,        ONLY : GET_FRAC_UNDER_PBLTOP
    Use PhysConstants,      Only : AirMW
!
! !INPUT PARAMETERS:
!
    LOGICAL,          INTENT(IN   )  :: am_I_Root  ! root CPU?
    TYPE(MetState),   INTENT(IN   )  :: State_Met  ! Met state
    TYPE(ChmState),   INTENT(IN   )  :: State_Chm  ! Chemistry state 
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(OptInput),   INTENT(INOUT)  :: Input_Opt  ! Input opts
    INTEGER,          INTENT(INOUT)  :: RC         ! Failure or success
!
! !REVISION HISTORY: 
!  27 Aug 2014 - C. Keller   - Initial version 
!  16 Jun 2016 - J. Sheng    - Added tracer index retriever
!  20 Jun 2016 - R. Yantosca - Now define species IDs only in the INIT phase
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Loop indices
    Integer           :: I, J, L
    ! Species index
    Integer           :: id_Spc
    ! Pointer to the species array
    Real(fp), Pointer :: Spc(:,:,:,:)

    !=================================================================
    ! FIXSFCVMR begins here!
    !=================================================================

    ! Assume success
    RC        = GC_SUCCESS

    ! Get a pointer to the species array
    Spc => State_Chm%Species

    ! ---------------------------------------------------
    ! JAS, 9/17/15: 
    ! Set mixing ratio of CH3Cl, CH2Cl2, and CHCl3 in PBL
    ! SDE 2016-12-14: This now replaces the UCX routine
    ! which previously gave only CH3Cl.
    ! ---------------------------------------------------
    ! Set CH3Cl mixing ratio in PBL
    id_Spc = Ind_('CH3Cl')
    IF ( id_Spc > 0 ) THEN
!$OMP PARALLEL DO                                                 &
!$OMP DEFAULT( SHARED )                                           &
!$OMP PRIVATE( I, J, L )
       DO L = 1, LLPAR
       DO J = 1, JJPAR
       DO I = 1, IIPAR
          IF (GET_FRAC_UNDER_PBLTOP(I,J,L)>0e+0_fp) THEN
             Spc(I,J,L,id_Spc) = 550e-12_fp / ( AIRMW / &
                State_Chm%SpcData(id_Spc)%Info%emMW_g )
          ENDIF  ! end selection of PBL boxes
       ENDDO
       ENDDO
       ENDDO
!$OMP END PARALLEL DO
    ENDIF

    ! Set CH2Cl2 mixing ratio in PBL
    id_Spc = Ind_('CH2Cl2')
    IF ( id_Spc > 0 ) THEN
!$OMP PARALLEL DO                                                 &
!$OMP DEFAULT( SHARED )                                           &
!$OMP PRIVATE( I, J, L )
       DO L = 1, LLPAR
       DO J = 1, JJPAR
       DO I = 1, IIPAR
          IF (GET_FRAC_UNDER_PBLTOP(I,J,L)>0e+0_fp) THEN
             Spc(I,J,L,id_Spc) = 20e-12_fp / ( AIRMW / &
                State_Chm%SpcData(id_Spc)%Info%emMW_g )
          ENDIF  ! end selection of PBL boxes
       ENDDO
       ENDDO
       ENDDO
!$OMP END PARALLEL DO
    ENDIF

    ! Set CHCl3 mixing ratio in PBL
    id_Spc = Ind_('CHCl3')
    IF ( id_Spc > 0 ) THEN
!$OMP PARALLEL DO                                                 &
!$OMP DEFAULT( SHARED )                                           &
!$OMP PRIVATE( I, J, L )
       DO L = 1, LLPAR
       DO J = 1, JJPAR
       DO I = 1, IIPAR
          IF (GET_FRAC_UNDER_PBLTOP(I,J,L)>0e+0_fp) THEN
             Spc(I,J,L,id_Spc) = 7e-12_fp / ( AIRMW / &
                State_Chm%SpcData(id_Spc)%Info%emMW_g )
          ENDIF  ! end selection of PBL boxes
       ENDDO
       ENDDO
       ENDDO
!$OMP END PARALLEL DO
    ENDIF

    ! Nullify pointer
    Nullify(Spc)

    IF ( RC/=GC_SUCCESS ) RETURN 

  END SUBROUTINE fixSfcVMR
!EOC
END MODULE SFCVMR_MOD
