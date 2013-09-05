!------------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 910.1 and      !
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: GIGC_Chem_Utils
!
! !DESCRIPTION: Utility module for the ESMF interface to GEOS-Chem.
!\\
!\\
! !INTERFACE:
!
MODULE GIGC_Chem_Utils
!
! !USES:
!
  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC :: GET_SPC_INDX
  PUBLIC :: GET_PBLLEV
!
! !REVISION HISTORY:
!  09 Oct 2012 - M. Long     - Initial version
!  09 Oct 2012 - R. Yantosca - Added ProTeX headers
!  09 Oct 2012 - R. Yantosca - Use F90 free-format indenting (Emacs F90 mode)
!  22 Oct 2012 - R. Yantosca - Renamed to gigc_chem_utils.F90

!EOP
!------------------------------------------------------------------------------
!BOC
  CONTAINS
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_spc_indx
!
! !DESCRIPTION: Returns the index of a chemical species given its name.
!\\
!\\
! !INTERFACE:
!
  FUNCTION Get_Spc_Indx( SPC_NAME, GC_SPC_IDS, GC_SPC_NAMES ) RESULT( SPC_INDX )
!
! !INPUT PARAMETERS:
!
    CHARACTER(LEN=*), INTENT(IN) :: SPC_NAME          ! Species name
    CHARACTER(LEN=*), INTENT(IN) :: GC_SPC_NAMES(:)   ! Names of all species
    INTEGER,          INTENT(IN) :: GC_SPC_IDS(:)     ! ID's  of all species
!
! !RETURN VALUE:
!
    INTEGER                      :: SPC_INDX          ! Index of this species 
!
! !REVISION HISTORY: 
!  09 Oct 2012 - M. Long     - Initial version
!  09 Oct 2012 - R. Yantosca - Added ProTeX headers
!  09 Oct 2012 - R. Yantosca - Declare INDX as integer return value
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER :: M

    ! Initialize
    SPC_INDX = -1

    ! Loop over all species names
    DO M = 1, SIZE( GC_SPC_NAMES )

       ! Return the index of the sought-for species
       IF( TRIM( SPC_NAME ) == TRIM( GC_SPC_NAMES(M) ) ) THEN
          SPC_INDX = GC_SPC_IDS(M)
          EXIT
       ENDIF

    ENDDO

  END FUNCTION Get_Spc_Indx
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_pbllev
!
! !DESCRIPTION: Get the level of the PBL
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE get_pbllev!(pblht,zmid,ncol,pbllev)
!
! !USES:
!
    !USE ppgrid
!
! !INPUT PARAMETERS:
!
    ! Comment these out for now
    !integer,  intent(in)   :: ncol ! Columns per chunk.
    !
    ! Calculate pbllev_cam across ncol
    !real, intent(in)   :: pblht(pcols)
    !real, intent(in)   :: zmid(pcols, pver)
    !real, intent(out)  :: pbllev(pcols)

!
! !REMARKS:
!  Comment out until further notice.
!
! !REVISION HISTORY: 
!  09 Oct 2012 - M. Long     - Initial version
!  09 Oct 2012 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!

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
!EOC

END MODULE GIGC_Chem_Utils
