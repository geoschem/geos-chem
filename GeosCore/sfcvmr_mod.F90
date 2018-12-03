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
    USE ERROR_MOD,          ONLY : ERROR_STOP
    USE Input_Opt_Mod,      ONLY : OptInput
    USE State_Met_Mod,      ONLY : MetState
    USE State_Chm_Mod,      ONLY : ChmState
    USE State_Chm_Mod,      ONLY : Ind_
    USE GC_GRID_MOD,        ONLY : GET_YMID
    USE TIME_MOD,           ONLY : Get_Month

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
    Integer           :: I, J, L, MONTH
    ! Species index
    Integer           :: id_Spc
    ! Pointer to the species array
    Real(fp), Pointer :: Spc(:,:,:,:)
    
    Real(fp)          :: A3090S (12), A0030S (12), A0030N(12)
    Real(fp)          :: A3060N (12), A6090N (12) 
    Real(fp)          :: YLAT, PPT

    !=================================================================
    ! FIXSFCVMR begins here!
    !=================================================================

    ! Assume success
    RC        = GC_SUCCESS

    ! Get a pointer to the species array
    Spc => State_Chm%Species

    ! Current month
    MONTH = Get_Month()

    YLAT = 0.0
    PPT = 0.0_fp

    ! ---------------------------------------------------
    ! JAS, 9/17/15: 
    ! Set mixing ratio of CH3Cl, CH2Cl2, and CHCl3 in PBL
    ! SDE 2016-12-14: This now replaces the UCX routine
    ! which previously gave only CH3Cl.
    ! ---------------------------------------------------
    ! Set CH3Cl mixing ratio in PBL
    id_Spc = Ind_('CH3Cl')

    A3090S = (/514,512,512,509,512,521,530,533,536,534,531,527/)
    A0030S = (/549,549,553,564,574,581,587,575,564,556,549,544/)
    A0030N = (/605,610,609,613,607,580,554,537,539,557,570,579/)
    A3060N = (/573,572,592,598,590,572,545,527,512,532,547,554/)
    A6090N = (/532,547,574,582,572,539,495,464,466,488,503,515/)

    IF ( id_Spc > 0 ) THEN
!$OMP PARALLEL DO                                                 &
!$OMP DEFAULT( SHARED )                                           &
!$OMP PRIVATE( I, J, L, YLAT, PPT)
       DO L = 1, LLPAR
       DO J = 1, JJPAR
       DO I = 1, IIPAR
          YLAT = GET_YMID( I, J, L )
          IF ( YLAT < -30.0_fp ) THEN
             PPT = A3090S(MONTH)
          ELSE IF ( YLAT >= -30.0_fp .and. YLAT < 0.0_fp  ) THEN
             PPT = A0030S(MONTH)
          ELSE IF ( YLAT >=   0.0_fp .and. YLAT < 30.0_fp ) THEN
             PPT = A0030N(MONTH)
          ELSE IF ( YLAT >=  30.0_fp .and. YLAT < 60.0_fp ) THEN
             PPT = A3060N(MONTH)
          ELSE
             PPT = A6090N(MONTH)
          ENDIF
          IF (GET_FRAC_UNDER_PBLTOP(I,J,L)>0e+0_fp) THEN
             Spc(I,J,L,id_Spc) = PPT*1e-12_fp / ( AIRMW / &
                State_Chm%SpcData(id_Spc)%Info%emMW_g )
          ENDIF  ! end selection of PBL boxes
       ENDDO
       ENDDO
       ENDDO
!$OMP END PARALLEL DO
    ENDIF

    ! Set CH2Cl2 mixing ratio in PBL
    id_Spc = Ind_('CH2Cl2')

    A3090S = (/10,10,10,11,12,12,13,14,13,13,12,11/)
    A0030S = (/15,17,17,18,18,18,19,20,19,18,17,18/)
    A0030N = (/45,45,44,39,45,42,39,34,31,29,35,46/)
    A3060N = (/54,56,56,58,59,52,46,41,41,45,50,54/)
    A6090N = (/65,66,63,63,61,56,51,47,44,47,54,61/)

    IF ( id_Spc > 0 ) THEN
!$OMP PARALLEL DO                                                 &
!$OMP DEFAULT( SHARED )                                           &
!$OMP PRIVATE( I, J, L, YLAT, PPT)
       DO L = 1, LLPAR
       DO J = 1, JJPAR
       DO I = 1, IIPAR
          YLAT = GET_YMID( I, J, L )
          IF ( YLAT < -30.0_fp ) THEN
             PPT = A3090S(MONTH)
          ELSE IF ( YLAT >= -30.0_fp .and. YLAT < 0.0_fp  ) THEN
             PPT = A0030S(MONTH)
          ELSE IF ( YLAT >=   0.0_fp .and. YLAT < 30.0_fp ) THEN
             PPT = A0030N(MONTH)
          ELSE IF ( YLAT >=  30.0_fp .and. YLAT < 60.0_fp ) THEN
             PPT = A3060N(MONTH)
          ELSE
             PPT = A6090N(MONTH)
          ENDIF
          IF (GET_FRAC_UNDER_PBLTOP(I,J,L)>0e+0_fp) THEN
             Spc(I,J,L,id_Spc) = PPT*1e-12_fp / ( AIRMW / &
                State_Chm%SpcData(id_Spc)%Info%emMW_g )
          ENDIF  ! end selection of PBL boxes
       ENDDO
       ENDDO
       ENDDO
!$OMP END PARALLEL DO
    ENDIF

    ! Set CHCl3 mixing ratio in PBL
    id_Spc = Ind_('CHCl3')

    A3090S = (/5, 5, 5, 5, 6, 6, 7, 7, 7, 6, 6, 6 /)
    A0030S = (/5, 5, 5, 6, 6, 6, 6, 6, 6, 5, 5, 5 /)
    A0030N = (/10,10,10,8, 9, 9, 9, 8, 8, 7, 8, 10/)
    A3060N = (/14,14,14,14,14,12,11,12,14,14,15,14/)
    A6090N = (/16,15,15,14,14,13,13,13,13,14,15,16/)

    IF ( id_Spc > 0 ) THEN
!$OMP PARALLEL DO                                                 &
!$OMP DEFAULT( SHARED )                                           &
!$OMP PRIVATE( I, J, L, YLAT )
       DO L = 1, LLPAR
       DO J = 1, JJPAR
       DO I = 1, IIPAR
          YLAT = GET_YMID( I, J, L )
          IF ( YLAT < -30.0_fp ) THEN
             PPT = A3090S(MONTH)
          ELSE IF ( YLAT >= -30.0_fp .and. YLAT < 0.0_fp  ) THEN
             PPT = A0030S(MONTH)
          ELSE IF ( YLAT >=   0.0_fp .and. YLAT < 30.0_fp ) THEN
             PPT = A0030N(MONTH)
          ELSE IF ( YLAT >=  30.0_fp .and. YLAT < 60.0_fp ) THEN
             PPT = A3060N(MONTH)
          ELSE
             PPT = A6090N(MONTH)
          ENDIF
          IF (GET_FRAC_UNDER_PBLTOP(I,J,L)>0e+0_fp) THEN
             Spc(I,J,L,id_Spc) = PPT*1e-12_fp / ( AIRMW / &
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
