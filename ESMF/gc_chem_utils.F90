MODULE GC_CHEM_UTILS

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: GET_SPC_INDX, GET_PBLLEV

  CONTAINS

  INTEGER FUNCTION GET_SPC_INDX( SPC_NAME, GC_SPC_IDS, GC_SPC_NAMES )
    !-----------------------------------------------------------------------
    !     ... RETURN OVERALL SPECIES INDEX ASSOCIATED WITH SPC_NAME
    !-----------------------------------------------------------------------

    USE COMODE_LOOP_MOD, ONLY : IGAS

    IMPLICIT NONE

    !-----------------------------------------------------------------------
    !     ... DUMMY ARGUMENTS
    !-----------------------------------------------------------------------
    CHARACTER(LEN=*), INTENT(IN) :: SPC_NAME
    CHARACTER(LEN=*), INTENT(IN) :: GC_SPC_NAMES(:)
    INTEGER, INTENT(IN)          :: GC_SPC_IDS(:)

    !-----------------------------------------------------------------------
    !     ... LOCAL VARIABLES
    !-----------------------------------------------------------------------
    INTEGER :: M

    GET_SPC_INDX = -1
    DO M = 1, SIZE(GC_SPC_NAMES)
       IF( TRIM( SPC_NAME ) == TRIM( GC_SPC_NAMES(M) ) ) THEN
          GET_SPC_INDX = GC_SPC_IDS(M)
          EXIT
       END IF
    END DO

  END FUNCTION GET_SPC_INDX

  SUBROUTINE get_pbllev!(pblht,zmid,ncol,pbllev)
!    USE ppgrid

!    integer,  intent(in)   :: ncol ! Columns per chunk.                                                                                    
!
!    ! Calculate pbllev_cam across ncol                                                                                                     
!    real, intent(in)   :: pblht(pcols)
!    real, intent(in)   :: zmid(pcols, pver)
!    real, intent(out)  :: pbllev(pcols)
!    integer                :: i, z
!
!    pbllev = 0
!
!    DO i=1, ncol
!       DO z=1, pver-1
!          IF (pblht(i) <= zmid(i,pver-(z))) THEN
!             IF (pblht(i) > zmid(i,pver-(z-1))) THEN
!                pbllev(i) = z
!             ENDIF
!          ENDIF
!       ENDDO
!    ENDDO
  END SUBROUTINE get_pbllev

end module gc_chem_utils
