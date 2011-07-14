! $Id: unix_cmds_mod.f,v 1.1 2009/11/20 21:43:01 bmy Exp $
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: unix_cmds_mod.f
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
! !REVISION HISTORY:
!  09 Jul 2004 - R. Yantosca - Initial version
!  20 Nov 2009 - R. Yantosca - Added ProTeX header
!EOP
!------------------------------------------------------------------------------
!BOC
      END MODULE UNIX_CMDS_MOD
!EOC
