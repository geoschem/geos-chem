#ifdef BPCH_DIAG
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: diag1.F90
!
! !DESCRIPTION: Subroutine DIAG1 accumulates diagnostic quantities on every
!  dynamic timestep.
!\\
!\\
! !INTERFACE:
!
SUBROUTINE DIAG1( Input_Opt, State_Chm, State_Grid, State_Met, RC )
!
! !USES:
!
  USE DIAG03_MOD,         ONLY : AD03_RGM, AD03_PBM, ND03
  USE ErrCode_Mod
  USE HCO_DIAGN_MOD,      ONLY : Diagn_Update
  USE HCO_ERROR_MOD
  USE Input_Opt_Mod,      ONLY : OptInput
  USE PhysConstants
  USE PRECISION_MOD
  USE Species_Mod,        ONLY : Species
  USE State_Chm_Mod,      ONLY : ChmState
  USE State_Chm_Mod,      ONLY : Ind_
  USE State_Grid_Mod,     ONLY : GrdState
  USE State_Met_Mod,      ONLY : MetState
  USE HCO_INTERFACE_MOD,  ONLY : HcoState

  IMPLICIT NONE
!
! !INPUT PARAMETERS:
!
  TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
  TYPE(GrdState), INTENT(IN)    :: State_Grid  ! Grid State object
  TYPE(MetState), INTENT(IN)    :: State_Met   ! Meteorology State object
!
! !INPUT/OUTPUT PARAMETERS:
!
  TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
!
! !OUTPUT PARAMETERS:
!
  INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
!
! !REMARKS:
!  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!  %%%  NOTE: THIS MODULE WILL BE A STUB UNLESS GEOS-Chem IS COMPILED    %%%
!  %%%  WITH THE BPCH_DIAG=y OPTION. (bmy, 10/4/19)                      %%%
!  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! !REVISION HISTORY:
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
  ! Scalars
  INTEGER       :: I,          J,          K
  INTEGER       :: L,          N,          NA
  REAL(fp)      :: EmMW_g,     P0,         Spc_VV

  ! SAVEd scalars
  LOGICAL, SAVE :: FIRST   = .TRUE.
  LOGICAL, SAVE :: Do_ND03 = .FALSE.
  INTEGER, SAVE :: id_Hg2  = -1
  INTEGER, SAVE :: id_HgP  = -1
  REAL(fp),SAVE :: EmMW_g_Hg2 = -1.0_fp
  REAL(fp),SAVE :: EmMW_g_HgP = -1.0_fp

  !=================================================================
  ! DIAG1 begins here!
  !=================================================================

  ! Initialize
  RC      =  GC_SUCCESS

  ! Are certain diagnostics turned on?
  IF ( .not. Input_Opt%ND03 > 0 ) RETURN

  !=================================================================
  ! First-time setup
  !=================================================================
  IF ( FIRST ) THEN

     ! Define species ID's  on the first call
     id_Hg2  = Ind_( 'Hg2' )
     id_HgP  = Ind_( 'HgP' )

     ! Pre-save the molecular weight of Hg2 [g]
     IF ( id_Hg2 > 0 ) THEN
        EmMW_g_Hg2 = State_Chm%SpcData(id_Hg2)%Info%EmMW_g
     ENDIF

     ! Pre-save the molecuar weight of HgP [s]
     IF ( id_HgP > 0 ) THEN
        EmMW_g_HgP = State_Chm%SpcData(id_HgP)%Info%EmMW_g
     ENDIF

     ! Reset first-time flag
     FIRST = .FALSE.
  ENDIF

  !=================================================================
  ! Archive diagnostics.  For better efficiency, place everything
  ! within a single OpenMP parallel loop, so that we can enclose
  ! more work within the parallel region than the prior code.
  !=================================================================
  IF ( Input_Opt%ND03 > 0 ) THEN

     ! Loop over advected species
     !$OMP PARALLEL DO &
     !$OMP DEFAULT( SHARED                             ) &
     !$OMP PRIVATE( I, J, L, N, NA, EmMw_g, P0, Spc_VV ) &
     !$OMP SCHEDULE( DYNAMIC, 4                        )
     DO NA = 1, State_Chm%nAdvect

        ! Initialize
        N      = State_Chm%Map_Advect(NA)          ! Species ID #
        EmMW_g = State_Chm%SpcData(N)%Info%EmMW_g  ! Emitted MW (g)

        ! Loop over grid boxes
        DO L = 1, State_Grid%NZ
        DO J = 1, State_Grid%NY
        DO I = 1, State_Grid%NX

           ! Initialize
           P0      = 0.0_fp
           Spc_VV  = 0.0_fp

           !===========================================================
           ! ND03: Hg speciality simulation diagnostics
           !===========================================================
           IF ( NA == 1 ) THEN

              !--------------------------------------------------------
              ! Reactive gaseous mercury [pptv]
              !--------------------------------------------------------
              IF ( id_Hg2 > 0 ) THEN
                 AD03_RGM(I,J,L) = AD03_RGM(I,J,L) &
                                 + State_Chm%Species(I,J,L,id_Hg2) &
                                 * ( AIRMW / EmMw_g_Hg2 * 1e+12_fp )
              ENDIF

              !--------------------------------------------------------
              ! Particulate bound mercury [pptv]
              !--------------------------------------------------------
              IF ( id_HgP > 0 ) THEN
                 AD03_PBM(I,J,L) = AD03_PBM(I,J,L) &
                                 + State_Chm%Species(I,J,L,id_HgP) &
                                 * ( AIRMW / EmMW_g_HgP * 1e+12_fp )
              ENDIF
           ENDIF
        ENDDO
        ENDDO
        ENDDO
     ENDDO
     !$OMP END PARALLEL DO
  ENDIF
END SUBROUTINE DIAG1
!EOC
#endif
