!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: hco_exts_mod
!
! !DESCRIPTION: Module HCO\_EXTS\_MOD contains routines to
! organize HEMCO extensions. 
! 
! \\
! !INTERFACE: 
!
      MODULE HCO_EXTS_MOD
!
! !USES:
!
      IMPLICIT NONE

      PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
      PUBLIC :: AddExt
      PUBLIC :: GetExtNr 
      PUBLIC :: ExtNrInUse
      PUBLIC :: ExtFinal
!
! !PRIVATE MEMBER FUNCTIONS:
!
!
! !REVISION HISTORY:
!  02 Oct 2013 - C. Keller: Initial version
!
!EOP
!-----------------------------------------------------------------------------
!BOC
!
! !MODULE TYPES/VARIABLES:
!
      TYPE Ext 
         CHARACTER(LEN=255)       :: ExtName
         INTEGER                  :: ExtNr
         TYPE(Ext), POINTER       :: NextExt
      END TYPE Ext

      ! Internally used ExtList 
      TYPE(Ext), POINTER          :: ExtList => NULL()
!
! !MODULE INTERFACES: 
!
      CONTAINS
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: AddExt 
!
! !DESCRIPTION: Subroutine AddExt adds the passed extension name to the
! list of used extensions. 
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE AddExt( ExtName, ExtNr ) 
!
! !USES:
!
      USE CHARPAK_MOD,  ONLY : TRANLC
!
! !ARGUMENTS:
!
      CHARACTER(LEN=*), INTENT(IN   )      :: ExtName 
      INTEGER,          INTENT(IN   )      :: ExtNr
!
! !REVISION HISTORY:
!  03 Oct 2013 - C. Keller - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !INTERNAL VARIABLES:
!
      TYPE(Ext), POINTER  :: NewExt => NULL()
      CHARACTER(LEN=255)  :: lcName

      !======================================================================
      ! AddExt 
      !======================================================================

      ! Allocate type 
      ALLOCATE ( NewExt )
      NewExt%NextExt => NULL()

      ! Set extension name and nr
      lcName         = TRIM(ExtName)
      CALL TRANLC( lcName )
      NewExt%ExtName = lcName
      NewExt%ExtNr   = ExtNr

      ! Place at beginning of extension list
      NewExt%NextExt => ExtList
      ExtList        => NewExt

      END SUBROUTINE AddExt 
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: GetExtNr 
!
! !DESCRIPTION: Function GetExtNr returns the extension number of
! extension ExtName. Returns -1 if no extension with the given name is
! found. 
!\\
!\\
! !INTERFACE:
!
      FUNCTION GetExtNr( ExtName ) Result ( ExtNr )
!
! !USES:
!
      USE CHARPAK_MOD,  ONLY : TRANLC
!
! !INPUT ARGUMENTS:
!
      CHARACTER(LEN=*), INTENT(IN   )      :: ExtName 
!
! !OUTPUT ARGUMENT:
!
      INTEGER                              :: ExtNr
!
! !REVISION HISTORY:
!  03 Oct 2013 - C. Keller - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !INTERNAL VARIABLES:
!
      TYPE(Ext), POINTER  :: ThisExt => NULL()
      CHARACTER(LEN=255)  :: LCname

      !======================================================================
      ! GetExtNr begins here
      !======================================================================

      ! Pass name to module 
      LCname = TRIM(ExtName)

      ! Init output
      ExtNr = -1 

      ! Point to header of extensions list
      ThisExt => ExtList

      ! Set to lower case
      CALL TRANLC( LCname )

      ! Loop over all used extensions and check if any of them matches
      ! ExtName.
      DO WHILE ( ASSOCIATED ( ThisExt ) ) 

         ! Compare extension names
         IF ( TRIM(ThisExt%ExtName) == TRIM(LCname) ) THEN
            ExtNr = ThisExt%ExtNr
            EXIT
         ENDIF

         ! Advance to next extension
         ThisExt => ThisExt%NextExt
      ENDDO

      ! Cleanup
      ThisExt => NULL()

      END FUNCTION GetExtNr 
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: ExtNrInUse
!
! !DESCRIPTION: Function ExtNrInUse checks if extension number ExtNr is
! in the list of used extensions or not. 
!\\
!\\
! !INTERFACE:
!
      FUNCTION ExtNrInUse( ExtNr ) Result ( InUse )
!
! !INPUT ARGUMENTS:
!
      INTEGER,          INTENT(IN   )      :: ExtNr 
!
! !OUTPUT ARGUMENT:
!
      LOGICAL                              :: InUse
!
! !REVISION HISTORY:
!  03 Oct 2013 - C. Keller - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !INTERNAL VARIABLES:
!
      TYPE(Ext), POINTER  :: ThisExt => NULL()
      CHARACTER(LEN=255)  :: LCname

      !======================================================================
      ! ExtNrInUse begins here
      !======================================================================

      ! Init output
      InUse = .FALSE.

      ! Point to header of extensions list
      ThisExt => ExtList

      ! Loop over all used extensions and check if any of them matches
      ! ExtName.
      DO WHILE ( ASSOCIATED ( ThisExt ) ) 

         ! Compare extension names
         IF ( ThisExt%ExtNr == ExtNr ) THEN
            InUse = .TRUE.
            EXIT
         ENDIF

         ! Advance to next extension
         ThisExt => ThisExt%NextExt
      ENDDO

      ! Cleanup
      ThisExt => NULL()

      END FUNCTION ExtNrInUse
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: ExtFinal
!
! !DESCRIPTION: Function ExtFinal finalizes this module. 
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE ExtFinal 
!
! !INPUT ARGUMENTS:
!
!
! !REVISION HISTORY:
!  03 Oct 2013 - C. Keller - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !INTERNAL VARIABLES:
!
      TYPE(Ext), POINTER  :: ThisExt => NULL()
      TYPE(Ext), POINTER  :: NextExt => NULL()

      !======================================================================
      ! ExtFinal begins here
      !======================================================================

      ! Point to header of extensions list
      ThisExt => ExtList

      ! Loop over all extensions and deallocate the types 
      DO WHILE ( ASSOCIATED ( ThisExt ) ) 

         ! First set pointer to next entry
         NextExt => ThisExt%NextExt

         ! Now clean up this entry
         ThisExt%NextExt => NULL()
         DEALLOCATE ( ThisExt )

         ! Advance to next extension
         ThisExt => NextExt
      ENDDO

      ! Cleanup
      ThisExt => NULL()

      END SUBROUTINE ExtFinal 
!EOC
      END MODULE HCO_EXTS_MOD 
!EOM
