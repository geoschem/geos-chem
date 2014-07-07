!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: hcox_extlist_mod
!
! !DESCRIPTION: Module HCOX\_EXTLIST\_MOD contains routines and
! variables to organize HEMCO extensions and the corresponding
! settings.
! This is done through the ExtList object, which is a simple list 
! containing all enabled HEMCO extensions (name and ext. ID) and
! the corresponding options as defined in the HEMCO configuration
! file. 
! \\
! !INTERFACE: 
!
      MODULE HCOX_EXTLIST_MOD
!
! !USES:
!
      USE HCO_ERROR_MOD

      IMPLICIT NONE
      PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
      PUBLIC :: AddExt
      PUBLIC :: AddExtOpt
      PUBLIC :: GetExtOpt
      PUBLIC :: GetExtHcoID
      PUBLIC :: GetExtNr 
      PUBLIC :: ExtNrInUse
      PUBLIC :: ExtFinal
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
      ! Type holding name, species, options and extension number of
      ! an extensions (as defined in the HEMCO configuration file)
      TYPE Ext 
         CHARACTER(LEN=255)       :: ExtName  ! Name
         CHARACTER(LEN=255)       :: Spcs     ! Species
         CHARACTER(LEN=1023)      :: Opts     ! Options
         INTEGER                  :: ExtNr    ! Ext. number
         TYPE(Ext), POINTER       :: NextExt  ! next list item
      END TYPE Ext

      ! Private linked list carrying information of all enabled extensions 
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
! !DESCRIPTION: Subroutine AddExt adds a new extension to the extensions
! list. The extension name, number and species (multiple species separated
! by the HEMCO separator sign) need to be provided. Extension options are 
! left blank but can be added lateron using AddExtOpt. 
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE AddExt( ExtName, ExtNr, Spcs ) 
!
! !USES:
!
      USE CHARPAK_MOD,  ONLY : TRANLC
!
! !ARGUMENTS:
!
      CHARACTER(LEN=*), INTENT(IN   )      :: ExtName 
      INTEGER,          INTENT(IN   )      :: ExtNr
      CHARACTER(LEN=*), INTENT(IN   )      :: Spcs
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

      ! Extension species 
      NewExt%Spcs    = Spcs

      ! Initialize extension options. These will be filled lateron
      NewExt%Opts = ''

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
! !ROUTINE: AddExtOpt 
!
! !DESCRIPTION: Function AddExtOpt appends the given string to the options
! character of the desired extension (identified by its extension number.
! The options string is expected to contain an option name and value,
! separated by a colon (:).
! Function GetExtOpt can be used to extract the option value at a later
! point. 
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE AddExtOpt ( Opt, ExtNr, RC ) 
!
! !INPUT ARGUMENTS:
!
      CHARACTER(LEN=*), INTENT(IN   )      :: Opt   ! Option name & value 
      INTEGER,          INTENT(IN   )      :: ExtNr ! Add to this extension
      INTEGER,          INTENT(INOUT)      :: RC
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
      CHARACTER(LEN=255)  :: LCname, MSG
      CHARACTER(LEN=1023) :: NewOpt

      !======================================================================
      ! AddExtOpt begins here
      !======================================================================

      ! Find extension of interest 
      ThisExt => ExtList
      DO WHILE ( ASSOCIATED ( ThisExt ) ) 
         IF ( ThisExt%ExtNr == ExtNr ) EXIT 
         ThisExt => ThisExt%NextExt
      ENDDO

      IF ( .NOT. ASSOCIATED( ThisExt ) ) THEN
         WRITE(MSG,*) 'Cannot find extension Nr. ', ExtNr
         CALL HCO_ERROR(MSG,RC,THISLOC='AddExtOpt (hco_config_mod)')
         RETURN
      ENDIF

      ! Create new options string. Eventually add to existing options, use
      ! ':::' to separate them. 
      NewOpt = TRIM(ThisExt%Opts)
      IF ( TRIM(NewOpt) == '' ) THEN
         NewOpt = TRIM(Opt)
      ELSE
         NewOpt = TRIM(NewOpt) // ':::' // TRIM(Opt)
      ENDIF

      ! Adjust options object
      ThisExt%Opts = ''
      ThisExt%Opts = NewOpt

      ! Cleanup and leave
      ThisExt => NULL()
      RC = HCO_SUCCESS

      END SUBROUTINE AddExtOpt 
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: GetExtOpt
!
! !DESCRIPTION: Function GetExtOpt returns the option value for a given 
! extension and option name. The type of the return value depends on the
! provided argument (real, boolean, character). 
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE GetExtOpt ( ExtNr,      OptName,   OptValSp,   &
                             OptValDp,   OptValInt, OptValBool, &
                             OptValChar, RC                      )
!
! !USES:
!
      USE CHARPAK_MOD,       ONLY : STRSPLIT, TRANLC
!
! !INPUT ARGUMENTS:
!
      INTEGER,          INTENT(IN   )            :: ExtNr 
      CHARACTER(LEN=*), INTENT(IN   )            :: OptName 
      REAL(sp),         INTENT(  OUT), OPTIONAL  :: OptValSp
      REAL(dp),         INTENT(  OUT), OPTIONAL  :: OptValDp
      INTEGER,          INTENT(  OUT), OPTIONAL  :: OptValInt
      LOGICAL,          INTENT(  OUT), OPTIONAL  :: OptValBool
      CHARACTER(LEN=*), INTENT(  OUT), OPTIONAL  :: OptValChar
      INTEGER,          INTENT(INOUT)            :: RC
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
      CHARACTER(LEN=255)  :: MSG, LOC, SUBSTR(255)
      CHARACTER(LEN=1023) :: STR
      INTEGER             :: N, cnt
      INTEGER             :: BGN, FIN, STRLEN, I

      !======================================================================
      ! GetExtOpt begins here
      !======================================================================

      ! Init
      LOC = 'GetExtOpt (hco_config_mod)'

      ! Find extension of interest 
      ThisExt => ExtList
      DO WHILE ( ASSOCIATED ( ThisExt ) ) 
         IF ( ThisExt%ExtNr == ExtNr ) EXIT 
         ThisExt => ThisExt%NextExt
      ENDDO

      ! Error if extension not found
      IF ( .NOT. ASSOCIATED( ThisExt ) ) THEN
         WRITE(MSG,*) 'Cannot find extension Nr. ', ExtNr
         CALL HCO_ERROR( MSG, RC, THISLOC=LOC )
         RETURN
      ENDIF

      ! Opts contains all options name-value pairs for the given
      ! extension, all separated by ':::'. Thus only check one
      ! option per time, i.e. walk through the Opts string until
      ! the option name of interest is found.
      ! BGN and FIN denote the start and end of the currently checked 
      ! Opts chunk. 
      cnt    = 0
      STRLEN = LEN(TRIM(ThisExt%Opts))
      BGN    = 1
      FIN    = STRLEN

      ! Do until option name of interest is found
      DO
         ! Reduce to valid range
         IF ( BGN >= STRLEN ) THEN
            MSG = '(A) Cannot find option ' // TRIM(OptName) // &
                  ' in ' // TRIM(ThisExt%Opts)
            CALL HCO_ERROR(MSG,RC,THISLOC=LOC )
            RETURN
         ENDIF

         ! I denotes the location of next option delimiter.
         I = INDEX( TRIM(ThisExt%Opts(BGN:FIN)), ':::' )
         IF ( I > 0 ) THEN
            FIN = BGN + I - 1
         ENDIF
 
         ! Extract value if this is the option of interest
         ! --> Assume here that options are separated by colon sign (:)
         !     and that the option value is in the second column.
         IF ( INDEX( ThisExt%Opts(BGN:FIN), TRIM(OptName) ) > 0 ) THEN
            CALL STRSPLIT( ThisExt%Opts(BGN:FIN), ':', SUBSTR, N )
            IF ( N < 2 ) THEN
               MSG = 'Option has too few elements: ' // &
                  TRIM(ThisExt%Opts(BGN:FIN)) 
               CALL HCO_ERROR(MSG,RC,THISLOC=LOC)
               RETURN
            ENDIF

            ! Pass option value to output
            IF ( PRESENT(OptValSp) ) THEN
               READ( SUBSTR(2), * ) OptValSp
            ELSEIF ( PRESENT(OptValDp) ) THEN
               READ( SUBSTR(2), * ) OptValDp
            ELSEIF ( PRESENT(OptValInt) ) THEN
               READ( SUBSTR(2), * ) OptValInt
            ELSEIF ( PRESENT(OptValBool) ) THEN
               CALL TRANLC( TRIM(SUBSTR(2)) )
               IF ( INDEX( TRIM(SUBSTR(2)), 'true' ) > 0 ) THEN
                  OptValBool = .TRUE.
               ELSE
                  OptValBool = .FALSE.
               ENDIF
            ELSEIF ( PRESENT(OptValChar) ) THEN
               OptValChar = ADJUSTL( TRIM(SUBSTR(2)) )
            ELSE
               MSG = 'Invalid option output element: ' // TRIM(OptName)
               CALL HCO_ERROR( MSG, RC, THISLOC=LOC )
               RETURN
            ENDIF

            ! Leave loop here
            EXIT
         ENDIF

         ! Update valid string range. The new chunk now starts at end+3
         ! (skip delimiter characters). 
         BGN = FIN + 3 
         FIN = STRLEN  ! End of string

         ! Error trap: don't allow more than 100 iterations!
         cnt = cnt + 1
         IF ( cnt > 100 ) THEN
            CALL HCO_ERROR('CNT>100',RC,THISLOC=LOC)
            RETURN
         ENDIF
      ENDDO 

      ! Cleanup and leave
      ThisExt => NULL()
      RC = HCO_SUCCESS

      END SUBROUTINE GetExtOpt
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: GetExtHcoID 
!
! !DESCRIPTION: Subroutine GetExtHcoID returns the HEMCO species IDs and
! names for all species assigned to the given extension (identified by its
! extension number).
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE GetExtHcoID( HcoState, ExtNr, HcoIDs, &
                              SpcNames, nSpc,  RC       ) 
!
! !USES:
!
      USE HCO_STATE_MOD,       ONLY : HCO_State, HCO_GetHcoID
      USE CHARPAK_MOD,         ONLY : STRSPLIT
!
! !ARGUMENTS:
!
      TYPE(HCO_State),               POINTER          :: HcoState 
      INTEGER,                       INTENT(IN   )    :: ExtNr       ! Extension Nr. 
      INTEGER,          ALLOCATABLE, INTENT(  OUT)    :: HcoIDs(:)   ! Species IDs
      CHARACTER(LEN=*), ALLOCATABLE, INTENT(  OUT)    :: SpcNames(:) ! Species names
      INTEGER,                       INTENT(INOUT)    :: nSpc        ! # of species
      INTEGER,                       INTENT(INOUT)    :: RC 
!
! !REVISION HISTORY:
!  10 Jan 2014 - C. Keller: Initialization (update)
!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL ARGUMENTS:
!
      INTEGER                   :: I, AS
      CHARACTER(LEN=255)        :: MSG, LOC, SpcStr, SUBSTR(255)
      TYPE(Ext), POINTER        :: ThisExt => NULL()

      !======================================================================
      ! GetExtHcoID begins here
      !======================================================================

      ! Enter
      LOC = 'GetExtHcoID (HCO_CONFIG_MOD.F90)'

      ! Find extension of interest 
      ThisExt => ExtList
      DO WHILE ( ASSOCIATED ( ThisExt ) ) 
         IF ( ThisExt%ExtNr == ExtNr ) EXIT 
         ThisExt => ThisExt%NextExt
      ENDDO

      IF ( .NOT. ASSOCIATED( ThisExt ) ) THEN
         WRITE(MSG,*) 'Cannot find extension Nr. ', ExtNr
         CALL HCO_ERROR( MSG, RC, THISLOC=LOC )
         RETURN
      ENDIF

      ! Get species string
      SpcStr  = TRIM(ThisExt%Spcs)
      ThisExt => NULL()

      ! Split character
      CALL STRSPLIT( SpcStr, HCO_SEP(), SUBSTR, nSpc )
      IF ( nSpc == 0 ) RETURN 

      ! Allocate arrays 
      IF ( ALLOCATED(HcoIDs  ) ) DEALLOCATE(HcoIDs  ) 
      IF ( ALLOCATED(SpcNames) ) DEALLOCATE(SpcNames) 
      ALLOCATE(HcoIDs(nSpc), SpcNames(nSpc), STAT=AS)
      IF ( AS/=0 ) THEN
         CALL HCO_ERROR('HcoIDs Alloc. error', RC, THISLOC=LOC)
         RETURN
      ENDIF

      ! Extract species information
      DO I = 1, nSpc
         HcoIDs(I)   = HCO_GetHcoID( TRIM(SUBSTR(I)), HcoState )
         SpcNames(I) = TRIM(SUBSTR(I))
      ENDDO

      ! Return w/ success
      RC = HCO_SUCCESS 

      END SUBROUTINE GetExtHcoID 
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
! !DESCRIPTION: Function ExtFinal finalizes the extensions list. 
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE ExtFinal 
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
      END MODULE HCOX_EXTLIST_MOD 
!EOM
