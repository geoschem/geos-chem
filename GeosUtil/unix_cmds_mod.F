!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: unix_cmds_mod
!
! !DESCRIPTION: Module UNIX\_CMDS\_MOD contains variables which contain file 
!  suffixes and various Unix command strings.
!\\
!\\
! !INTERFACE: 
!
      MODULE UNIX_CMDS_MOD
! 
! !USES:
!
      IMPLICIT NONE
      PUBLIC
!
! !PUBLIC DATA MEMBERS:
!
      ! Unix cmd and file suffix strings for ...
      CHARACTER(LEN=255) :: BACKGROUND   ! Background operator  ( ' &'    ) 
      CHARACTER(LEN=255) :: REDIRECT     ! Redirection operator ( ' >'    )
      CHARACTER(LEN=255) :: REMOVE_CMD   ! File/dir remove cmd  ( 'rm'    )
      CHARACTER(LEN=255) :: SEPARATOR    ! Dir path separator   ( '/'     )
      CHARACTER(LEN=255) :: SPACE        ! Blank space          ( ' '     )
      CHARACTER(LEN=255) :: UNZIP_CMD    ! Unzip command        ( 'gzcat' )
      CHARACTER(LEN=255) :: WILD_CARD    ! Wild card operator   ( '*'     )
      CHARACTER(LEN=255) :: A3_SUFFIX    ! !%%% OBSOLETE %%%
      CHARACTER(LEN=255) :: A6_SUFFIX    ! !%%% OBSOLETE %%%
      CHARACTER(LEN=255) :: I6_SUFFIX    ! !%%% OBSOLETE %%%
      CHARACTER(LEN=255) :: PH_SUFFIX    ! !%%% OBSOLETE %%%
      CHARACTER(LEN=255) :: KZZ_SUFFIX   ! !%%% OBSOLETE %%%
      CHARACTER(LEN=255) :: GRID_SUFFIX  ! !%%% OBSOLETE %%%
      CHARACTER(LEN=255) :: ZIP_SUFFIX   ! Zipped file suffix   ( '.gz'   )
!
!
! !REMARKS:
!  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!  %%%  NOTE: THIS MODULE IS NOW OBSOLETE AND IS SLATED TO BE REMOVED!   %%%
!  %%%  WE NOW USE THE UNIX CMDS FROM Input_Opt INSTEAD OF FROM HERE.    %%%
!  %%%     -- Bob Yantosca, 20 Aug 2013                                  %%%
!  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! !REVISION HISTORY:
!  09 Jul 2004 - R. Yantosca - Initial version
!  20 Nov 2009 - R. Yantosca - Added ProTeX header
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
!EOP
!------------------------------------------------------------------------------
!BOC
      END MODULE UNIX_CMDS_MOD
!EOC
