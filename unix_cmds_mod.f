! $Id: unix_cmds_mod.f,v 1.1 2004/09/21 18:04:20 bmy Exp $
      MODULE UNIX_CMDS_MOD
!
!******************************************************************************
!  Module UNIX_CMDS_MOD contains variables which contain file suffixes and
!  Unix command strings that are used to unzip met field data. (bmy, 7/9/04)
!
!  Module Variables:
!  ============================================================================
!  (1 ) BACKGROUND  : String for background operator  (' &'    in Unix) 
!  (2 ) REDIRECT    : String for redirection operator (' >'    in Unix)
!  (3 ) REMOVE_CMD  : String for remove command       ('rm'    in Unix)
!  (4 ) SEPARATOR   : String for dir path separator   ('/'     in Unix)
!  (5 ) SPACE       : String for blank spaces         (' '     in Unix)    
!  (6 ) STAR        : String for wild card operator   ('*'     in Unix)
!  (7 ) UNZIP_CMD   : String for unzip command        ('gzcat' in Unix) 
!  (8 ) A3_SUFFIX   : Suffix for DAO A-3  (Average 3h      ) met fields 
!  (9 ) A6_SUFFIX   : Suffix for DAO A-6  (Average 6h      ) met fields 
!  (10) I6_SUFFIX   : Suffix for DAO I-6  (Instantaneous 6h) met fields 
!  (11) PH_SUFFIX   : Suffix for DAO PHIS (geopotential hts) met fields 
!  (12) KZZ_SUFFIX  : Suffix for DAO KZZ  (Average 3h      ) met fields 
!  (13) GRID_SUFFIX : Suffix for grid resolution
!  (14) ZIP_SUFFIX  : Suffix for denoting compressed files
!
!   NOTES:
!******************************************************************************
!
      IMPLICIT NONE
      
      !=================================================================
      ! MODULE VARIABLES
      !=================================================================
      CHARACTER(LEN=255) :: BACKGROUND
      CHARACTER(LEN=255) :: REDIRECT 
      CHARACTER(LEN=255) :: REMOVE_CMD
      CHARACTER(LEN=255) :: SEPARATOR
      CHARACTER(LEN=255) :: SPACE
      CHARACTER(LEN=255) :: UNZIP_CMD
      CHARACTER(LEN=255) :: WILD_CARD  
      CHARACTER(LEN=255) :: A3_SUFFIX
      CHARACTER(LEN=255) :: A6_SUFFIX
      CHARACTER(LEN=255) :: I6_SUFFIX
      CHARACTER(LEN=255) :: PH_SUFFIX
      CHARACTER(LEN=255) :: KZZ_SUFFIX
      CHARACTER(LEN=255) :: GRID_SUFFIX
      CHARACTER(LEN=255) :: ZIP_SUFFIX

      ! End of module
      END MODULE UNIX_CMDS_MOD
