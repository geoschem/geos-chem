!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: hco_extlist_mod
!
! !DESCRIPTION: Module HCO\_EXTLIST\_MOD contains routines and
! variables to organize HEMCO extensions and the corresponding
! settings (options). This is done through the ExtList object,
! which is a simple list containing all enabled HEMCO extensions
! (name and ext. ID) and the corresponding options, as defined in
! the HEMCO configuration file. The general HEMCO settings are
! stored as options of the HEMCO core extension (Extension number
! = 0). The CORE extension is activated in every HEMCO run, while
! all other extensions are only activated if enabled in the
! configuration file.
!\\
!\\
! Extension number -999 is used as 'wildcard' value, e.g. data
! containers with extension number -999 will always be read by
! HEMCO but will be ignored for emission calculation. This is
! particularly useful for data fields that shall be used outside
! of HEMCO, e.g. stratospheric chemistry prod/loss rates, etc.
!\\
!\\
! Extension options are 'flexible' in a sense that any option
! name/value pair can be assigned to an extension. The value of
! any of these options can be queried using subroutine GetExtOpt
! or function HCO\_GetOpt. In fact, the HEMCO filename parser
! (in hco\_chartools\_mod.F90) will attempt to find an option
! value for any HEMCO 'token' (a character starting with the
! HEMCO token sign (which is, the dollar sign '\$'). This allows
! the user to specify as many individual tokens as HEMCO
! settings as needed.
!\\
!\\
! !INTERFACE:
!
MODULE HCO_ExtList_Mod
!
! !USES:
!
  USE HCO_Error_Mod
  USE HCO_Types_Mod

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC  :: AddExt
  PUBLIC  :: AddExtOpt
  PUBLIC  :: GetExtOpt
  PUBLIC  :: GetExtNr
  PUBLIC  :: GetExtSpcStr
  PUBLIC  :: GetExtSpcVal
  PUBLIC  :: SetExtNr
  PUBLIC  :: ExtNrInUse
  PUBLIC  :: ExtFinal
  PUBLIC  :: HCO_GetOpt
  PUBLIC  :: HCO_SetDefaultToken
  PUBLIC  :: HCO_ROOT

  PRIVATE :: HCO_AddOpt
  PRIVATE :: HCO_CleanupOpt

  ! Core extension number
  INTEGER, PARAMETER, PUBLIC  :: CoreNr = -1
!
! !REVISION HISTORY:
!  02 Oct 2013 - C. Keller   - Initial version
!  01 Jul 2014 - R. Yantosca - Now use F90 free-format indentation
!  01 Jul 2014 - R. Yantosca - Cosmetic changes in ProTeX headers
!  30 Sep 2014 - R. Yantosca - ThisExt%Spcs now has 2047 chars for extensions
!                              having many individual species
!  20 Sep 2015 - C. Keller   - Reorganize options in linked lists. Tokens are
!                              now the same as options and can be flexibly
!                              set by the user.
!  24 Aug 2017 - M. Sulprizio- Remove support for GCAP, GEOS-4, GEOS-5 and MERRA
!  29 Aug 2018 - M. Sulprizio- Users can now specify $MET or $met to use
!                              uppercase or lowercase strings for met field
!EOP
!-----------------------------------------------------------------------------
!BOC
!
! !PRIVATE TYPES:
!
  ! Lenght of maximum token character length
  CHARACTER(LEN=OPTLEN), PARAMETER  :: EMPTYOPT = '---'

  !---------------------------------------------------------------------------
  ! Default tokens
  ! HEMCO has three tokens that can be specified in the HEMCO configuration
  ! file: ROOT (root directory), MET/met (met fields), and RES (horizontal
  ! resolution). These tokens can be used in file names to be dynamically
  ! replaced, e.g. file.$MET.$RES.nc becomes file.GEOSFP.4x5.nc if MET is set
  ! to 'GEOSFP' and RES to '4x5'. Opts are also allowed for dates ($YYYY,
  ! $MM, $DD, $HH, see routine HCO_CharParse).
  ! The default tokens below will be used if by default, i.e. if the
  ! corresponding token is not specified in the HEMCO configuration file.
  !---------------------------------------------------------------------------

  ! Default root directory
  CHARACTER(LEN=OPTLEN), PARAMETER :: DEF_ROOT = '/please/provide/root/path'

  ! Default values for characters that can be changed
  ! through the configuration file
  CHARACTER(LEN=1),    PARAMETER :: DEF_COLON     = ':'
  CHARACTER(LEN=1),    PARAMETER :: DEF_SEPARATOR = '/'
  CHARACTER(LEN=1),    PARAMETER :: DEF_WILDCARD  = '*'

  ! Met field and grid tokens
  CHARACTER(LEN=15)              :: DEF_MET_UC = 'UNKNOWN_MET'
  CHARACTER(LEN=15)              :: DEF_MET_LC = 'unknown_met'
  CHARACTER(LEN=15)              :: DEF_MET_EXT= 'unknown_extension'
  CHARACTER(LEN=15)              :: DEF_CN_YR  = 'unknown_year'
  CHARACTER(LEN=15)              :: DEF_RES    = 'unknown_res'
  CHARACTER(LEN=15)              :: DEF_NC_VER = 'nc'

  INTERFACE GetExtSpcVal
     MODULE PROCEDURE GetExtSpcVal_Char
     MODULE PROCEDURE GetExtSpcVal_Int
     MODULE PROCEDURE GetExtSpcVal_Sp
  END INTERFACE GetExtSpcVal

CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: AddExt
!
! !DESCRIPTION: Subroutine AddExt adds a new extension to the extensions
! list. The extension name, number and species (multiple species separated
! by the HEMCO separator sign) need to be provided. Extension options are
! left blank but can be added lateron using AddExtOpt.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE AddExt( HcoConfig, ExtName, ExtNr, InUse, Spcs, RC )
!
! !USES:
!
    USE CHARPAK_MOD,  ONLY : TRANLC
!
! !INPUT PARAMETERS::
!
    TYPE(ConfigObj),  POINTER       :: HcoConfig
    CHARACTER(LEN=*), INTENT(IN   ) :: ExtName
    INTEGER,          INTENT(IN   ) :: ExtNr
    LOGICAL,          INTENT(IN   ) :: InUse
    CHARACTER(LEN=*), INTENT(IN   ) :: Spcs
!
! !INPUT/OUTPUT PARAMETERS::
!
    INTEGER,          INTENT(INOUT) :: RC
!
! !REVISION HISTORY:
!  03 Oct 2013 - C. Keller - Initial version
!  20 Sep 2015 - C. Keller - Options are now linked list
!  12 Dec 2015 - C. Keller - Added argument InUse
!  26 Oct 2016 - R. Yantosca - Don't nullify local ptrs in declaration stmts
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !INTERNAL VARIABLES:
!
    INTEGER             :: OrigExtNr
    CHARACTER(LEN=255)  :: MSG, lcName
    TYPE(Ext), POINTER  :: NewExt
    TYPE(Ext), POINTER  :: ThisExt

    !======================================================================
    ! AddExt
    !======================================================================

    ! Init
    NewExt  => NULL()
    ThisExt => NULL()

    ! All extension names are lower case
    lcName = TRIM(ExtName)
    CALL TRANLC( lcName )

    ! Check if extension already exists
    OrigExtNr = GetExtNr( HcoConfig%ExtList, TRIM(lcName) )
    IF ( OrigExtNr /= -999 ) THEN

       ! Return w/ error if extension numbers do not match
       IF ( OrigExtNr /= ExtNr ) THEN
          WRITE(MSG,*) 'Cannot create extension - extension already exists', &
                       TRIM(lcName), ExtNr, OrigExtNr
          CALL HCO_ERROR(HcoConfig%Err,MSG,RC,THISLOC='AddExt (hco_extlist_mod.F90)')
          RETURN

       ! Nothing to do otherwise
       ELSE
          RC = HCO_SUCCESS
          RETURN
       ENDIF
    ENDIF

    ! Allocate type
    ALLOCATE ( NewExt )
    NewExt%NextExt => NULL()

    ! Set extension name
    NewExt%ExtName = lcName

    ! Set extension number and species. Set to invalid values if not used.
    IF ( InUse ) THEN
       NewExt%ExtNr = ExtNr
       NewExt%Spcs  = Spcs
    ELSE
       NewExt%ExtNr = -1
       NewExt%Spcs  = 'None'
    ENDIF

    ! Initialize extension options. These will be filled lateron
    NewExt%Opts    => NULL()

    ! Place at end of extension list
    ! The rational for this is simply that the more often used extension options
    ! (of the HEMCO 'core' and 'base' extensions) are set first and we better have
    ! them at the beginning of ExtList for better efficiency.

    ! If this is the first entry...
    IF ( .NOT. ASSOCIATED(HcoConfig%ExtList) ) THEN
       HcoConfig%ExtList => NewExt

    ! Otherwise, scan list until last extension is encountered
    ELSE
       ThisExt => HcoConfig%ExtList
       DO WHILE(ASSOCIATED(ThisExt))
          IF ( .NOT. ASSOCIATED(ThisExt%NextExt) ) EXIT
          ThisExt => ThisExt%NextExt
       END DO

       ! Append new extension to the list
       ThisExt%NextExt => NewExt
    ENDIF

    ! Verbose
    IF ( HcoConfig%amIRoot .AND. HCO_IsVerb(HcoConfig%Err,2) .AND. InUse ) THEN
       WRITE(MSG,*) 'Added HEMCO extension: ', TRIM(ExtName), ExtNr
       CALL HCO_MSG(HcoConfig%Err,MSG)
    ENDIF

    ! Cleanup
    ThisExt => NULL()

    ! Return w/ success
    RC = HCO_SUCCESS

  END SUBROUTINE AddExt
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: AddExtOpt
!
! !DESCRIPTION: Function AddExtOpt appends the given string to the options
! character of the desired extension (identified by its extension number).
! The options string is expected to contain an option name and value,
! separated by a colon (:).
! Function GetExtOpt can be used to extract the option value at a later
! point.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE AddExtOpt( HcoConfig, Opt, ExtNr, RC, IgnoreIfExist )
!
! !USES:
!
  USE CHARPAK_MOD,       ONLY : STRSPLIT, TRANLC
!
! !INPUT PARAMETERS:
!
    TYPE(ConfigObj),  POINTER                 :: HcoConfig      ! Configuration object
    CHARACTER(LEN=*), INTENT(IN   )           :: Opt            ! Option name & value
    INTEGER,          INTENT(IN   )           :: ExtNr          ! Add to this extension
    LOGICAL,          INTENT(IN   ), OPTIONAL :: IgnoreIfExist  ! Ignore this entry if it exists already?
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(INOUT)           :: RC
!
! !REVISION HISTORY:
!  03 Oct 2013 - C. Keller - Initial version
!  20 Sep 2015 - C. Keller - Options are now linked list
!  12 Dec 2015 - C. Keller - Added argument IgnoreIfExist
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !INTERNAL VARIABLES:
!
    INTEGER               :: IDX
    CHARACTER(LEN=255)    :: MSG
    CHARACTER(LEN=OPTLEN) :: TmpStr, OptName, OptValue

    !======================================================================
    ! AddExtOpt begins here
    !======================================================================

    ! Parse option name and option value. These must be separated by colon.
    IDX = INDEX( TRIM(Opt), ':' )

    ! Error check
    IF ( IDX <= 0 ) THEN
       MSG = 'Cannot extract option name/value pair - these must be ' // &
             'separated by a colon (:) character: ' // TRIM(Opt)
       CALL HCO_ERROR(HcoConfig%Err,MSG,RC,THISLOC='AddExtOpt (hco_extlist_mod)')
       RETURN
    ENDIF

    ! Now split option name / value pair
    OptName  = Opt(1:(IDX-1))
    OptValue = Opt((IDX+1):LEN(Opt))

    ! Also check for '-->' option indicatior. This string must be stripped
    ! off the option name!
    IDX = INDEX( TRIM(OptName), '-->' )
    IF ( IDX > 0 ) THEN
       TmpStr  = OptName( (IDX+3) : LEN(TRIM(OptName)) )
       OptName = TmpStr
    ENDIF

    ! Pass to options
    CALL HCO_AddOpt( HcoConfig, OptName, OptValue, ExtNr, RC, &
                     IgnoreIfExist=IgnoreIfExist )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Cleanup and leave
    RC = HCO_SUCCESS

  END SUBROUTINE AddExtOpt
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: GetExtOpt
!
! !DESCRIPTION: Function GetExtOpt returns the option value for a given
! extension and option name. The type of the return value depends on the
! provided argument (real, boolean, character). The optional output
! argument FOUND returns TRUE if the given option name was found, and
! FALSE otherwise. If the FOUND argument is provided, no error is
! returned if the option name is not found!
! If the ExtNr is set to -999, the settings of all extensions are searched.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE GetExtOpt ( HcoConfig,  ExtNr,      OptName,   OptValHp, &
                         OptValSp,   OptValDp,   OptValInt,           &
                         OptValBool, OptValChar, Found,     RC )
!
! !USES:
!
  USE CHARPAK_MOD,       ONLY : STRSPLIT, TRANLC
!
! !INPUT PARAMETERS:
!
    TYPE(ConfigObj),  POINTER                  :: HcoConfig
    INTEGER,          INTENT(IN   )            :: ExtNr
    CHARACTER(LEN=*), INTENT(IN   )            :: OptName
!
! !OUTPUT PARAMETERS:
!
    REAL(hp),         INTENT(  OUT), OPTIONAL  :: OptValHp
    REAL(sp),         INTENT(  OUT), OPTIONAL  :: OptValSp
    REAL(dp),         INTENT(  OUT), OPTIONAL  :: OptValDp
    INTEGER,          INTENT(  OUT), OPTIONAL  :: OptValInt
    LOGICAL,          INTENT(  OUT), OPTIONAL  :: OptValBool
    CHARACTER(LEN=*), INTENT(  OUT), OPTIONAL  :: OptValChar
    LOGICAL,          INTENT(  OUT), OPTIONAL  :: Found
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(INOUT)            :: RC
!
! !REVISION HISTORY:
!  03 Oct 2013 - C. Keller   - Initial version
!  13 Jan 2015 - R. Yantosca - Add optional variable of flex precision (hp)
!  14 Feb 2015 - C. Keller   - Add option to search all extensions (ExtNr=-999).
!  17 Apr 2015 - C. Keller   - Passed option OptName must now exactly match the
!                              stored option name to avoid ambiguity.
!  20 Sep 2015 - C. Keller   - Options are now linked list.
!  20 Jan 2016 - C. Keller   - Bug fix: boolean options are now case insensitive.
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !INTERNAL VARIABLES:
!
    CHARACTER(LEN=OPTLEN) :: OptValue
    LOGICAL               :: OptFound
    CHARACTER(LEN=255)    :: MSG, LOC

    !======================================================================
    ! GetExtOpt begins here
    !======================================================================

    ! Init
    LOC = 'GetExtOpt (hco_extlist_mod)'

    ! Get option
    OptValue = HCO_GetOpt( HcoConfig%ExtList, OptName, ExtNr=ExtNr )
    IF ( TRIM(OptValue) == TRIM(EMPTYOPT) ) THEN
       OptFound = .FALSE.
    ELSE
       OptFound = .TRUE.
    ENDIF

    ! Check if option was found. Handling depends on presence
    ! of argument 'FOUND'. If FOUND is not present and option
    ! was not found, return with error.
    IF ( PRESENT(FOUND) ) THEN
       FOUND = OptFound
    ELSEIF ( .NOT. OptFound ) THEN
       WRITE(MSG,*) '(A) Cannot find option ', TRIM(OptName),  &
       ' in extension ', ExtNr
       CALL HCO_ERROR(HcoConfig%Err,MSG,RC,THISLOC=LOC )
       RETURN
    ENDIF

    ! Pass option value to output
    IF ( PRESENT(OptValSp) ) THEN
       IF ( OptFound ) THEN
          READ( OptValue, * ) OptValSp
       ELSE
          OptValSp = -999.0_sp
       ENDIF
    ELSEIF ( PRESENT(OptValDp) ) THEN
       IF ( OptFound ) THEN
          READ( OptValue, * ) OptValDp
       ELSE
          OptValDp = -999.0_dp
       ENDIF
    ELSEIF ( PRESENT(OptValHp) ) THEN
       IF ( OptFound ) THEN
          READ( OptValue, * ) OptValHp
       ELSE
          OptValHp = -999.0_hp
       ENDIF
    ELSEIF ( PRESENT(OptValInt) ) THEN
       IF ( OptFound ) THEN
          READ( OptValue, * ) OptValInt
       ELSE
          OptValInt = -999
       ENDIF
    ELSEIF ( PRESENT(OptValBool) ) THEN
       IF ( OptFound ) THEN
          CALL TRANLC( OptValue )
          IF ( INDEX( TRIM(OptValue), 'true' ) > 0 ) THEN
             OptValBool = .TRUE.
          ELSE
             OptValBool = .FALSE.
          ENDIF
       ELSE
          OptValBool = .FALSE.
       ENDIF
    ELSEIF ( PRESENT(OptValChar) ) THEN
       IF ( OptFound ) THEN
          OptValChar = ADJUSTL( TRIM(OptValue) )
       ELSE
          OptValChar = ''
       ENDIF
    ENDIF

    ! Cleanup and leave
    RC = HCO_SUCCESS

  END SUBROUTINE GetExtOpt
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: GetExtNr
!
! !DESCRIPTION: Function GetExtNr returns the extension number of
! extension ExtName. Returns -999 if no extension with the given name is
! found.
!\\
!\\
! !INTERFACE:
!
  FUNCTION GetExtNr( ExtList, ExtName ) Result ( ExtNr )
!
! !USES:
!
    USE CHARPAK_MOD,  ONLY : TRANLC
!
! !INPUT PARAMETERS:
!
    TYPE(Ext),        POINTER       :: ExtList
    CHARACTER(LEN=*), INTENT(IN   ) :: ExtName
!
! !RETURN VALUE:
!
    INTEGER                         :: ExtNr
!
! !REVISION HISTORY:
!  03 Oct 2013 - C. Keller - Initial version
!  26 Oct 2016 - R. Yantosca - Don't nullify local ptrs in declaration stmts
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !INTERNAL VARIABLES:
!
    TYPE(Ext), POINTER  :: ThisExt
    CHARACTER(LEN=255)  :: LCname

    !======================================================================
    ! GetExtNr begins here
    !======================================================================

    ! Init
    ThisExt => NULL()

    ! Pass name to module
    LCname = TRIM(ExtName)

    ! Set to lower case
    CALL TRANLC( LCname )

    ! Init output
    ExtNr = -999

    ! Point to header of extensions list
    ThisExt => ExtList

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
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: GetExtSpcStr
!
! !DESCRIPTION: Subroutine GetExtSpcStr returns the HEMCO species names
! string of all species assigned to the given extension (identified by its
! extension number).
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE GetExtSpcStr( HcoConfig, ExtNr, SpcStr, RC )
!
! !INPUT PARAMETERS:
!
    TYPE(ConfigObj),  POINTER       :: HcoConfig
    INTEGER,          INTENT(IN   ) :: ExtNr    ! Extension Nr.
!
! !OUTPUT PARAMETERS:
!
    CHARACTER(LEN=*), INTENT(  OUT) :: SpcStr   ! Species string
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(INOUT) :: RC       ! Success or failure?
!
! !REVISION HISTORY:
!  10 Jan 2014 - C. Keller: Initialization (update)
!  26 Oct 2016 - R. Yantosca - Don't nullify local ptrs in declaration stmts
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL ARGUMENTS:
!
    CHARACTER(LEN=255)        :: MSG, LOC
    TYPE(Ext), POINTER        :: ThisExt

    !======================================================================
    ! GetExtSpcStr begins here
    !======================================================================

    ! Enter
    LOC = 'GetExtSpcStr (hco_extlist_mod.F90)'
    RC  = HCO_SUCCESS
    ThisExt => NULL()

    ! Find extension of interest
    ThisExt => HcoConfig%ExtList
    DO WHILE ( ASSOCIATED ( ThisExt ) )
       IF ( ThisExt%ExtNr == ExtNr ) EXIT
       ThisExt => ThisExt%NextExt
    ENDDO

    IF ( .NOT. ASSOCIATED( ThisExt ) ) THEN
       WRITE(MSG,*) 'Cannot find extension Nr. ', ExtNr
       CALL HCO_ERROR( HcoConfig%Err, MSG, RC, THISLOC=LOC )
       RETURN
    ENDIF

    ! Get species string
    SpcStr  = TRIM(ThisExt%Spcs)
    ThisExt => NULL()

    ! Return w/ success
    RC = HCO_SUCCESS

    END SUBROUTINE GetExtSpcStr
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: GetExtSpcVal_Sp
!
! !DESCRIPTION: Subroutine GetExtSpcVal\_Sp returns single precision values
! associated with the species for a given extension. Specifically, this routine
! searches for extension setting '<Prefix>\_SpecName' for every species passed
! through input argument SpcNames and writes those into output argument SpcScal.
! The default value DefValue is assigned to all elements of SpcScal with no
! corresponding extension setting.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE GetExtSpcVal_Sp( HcoConfig, ExtNr,  NSPC,    SpcNames, &
                              Prefix,  DefValue, SpcScal, RC         )
!
! !INPUT PARAMETERS:
!
    TYPE(ConfigObj),       POINTER       :: HcoConfig
    INTEGER,               INTENT(IN   ) :: ExtNr          ! Extension Nr.
    INTEGER,               INTENT(IN   ) :: NSPC           ! # of species
    CHARACTER(LEN=*),      INTENT(IN   ) :: SpcNames(NSPC) ! Species string
    CHARACTER(LEN=*),      INTENT(IN   ) :: Prefix         ! search prefix
    REAL(sp),              INTENT(IN   ) :: DefValue       ! default value
!
! !INPUT/OUTPUT PARAMETERS:
!
    REAL(sp), ALLOCATABLE, INTENT(INOUT) :: SpcScal(:)     ! Species scale factors
    INTEGER,               INTENT(INOUT) :: RC             ! Success or failure?
!
! !REVISION HISTORY:
!  10 Jun 2015 - C. Keller - Initial version
!  20 Sep 2015 - C. Keller - Now allocate output array in this routine.
!EOP
!------------------------------------------------------------------------------
!BOC

    !======================================================================
    ! GetExtSpcVal_Sp begins here
    !======================================================================

    ! Make sure output is properly allocated
    IF ( ALLOCATED(SpcScal) ) DEALLOCATE(SpcScal)
    ALLOCATE(SpcScal(NSPC))
    SpcScal=DefValue

    CALL GetExtSpcVal_Dr ( HcoConfig, ExtNr, NSPC, SpcNames, Prefix, RC, &
                           DefVal_SP=DefValue, SpcScal_SP=SpcScal )

    END SUBROUTINE GetExtSpcVal_sp
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: GetExtSpcVal_Int
!
! !DESCRIPTION: Subroutine GetExtSpcVal\_Int returns integer values
! associated with the species for a given extension. Specifically, this routine
! searches for extension setting '<Prefix>\_SpecName' for every species passed
! through input argument SpcNames and writes those into output argument SpcScal.
! The default value DefValue is assigned to all elements of SpcScal with no
! corresponding extension setting.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE GetExtSpcVal_Int( HcoConfig, ExtNr,    NSPC,    SpcNames, &
                               Prefix,    DefValue, SpcScal, RC         )
!
! !INPUT PARAMETERS:
!
    TYPE(ConfigObj),       POINTER       :: HcoConfig
    INTEGER,               INTENT(IN   ) :: ExtNr          ! Extension Nr.
    INTEGER,               INTENT(IN   ) :: NSPC           ! # of species
    CHARACTER(LEN=*),      INTENT(IN   ) :: SpcNames(NSPC) ! Species string
    CHARACTER(LEN=*),      INTENT(IN   ) :: Prefix         ! search prefix
    INTEGER,               INTENT(IN   ) :: DefValue       ! default value
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,  ALLOCATABLE, INTENT(INOUT) :: SpcScal(:)     ! Species scale factors
    INTEGER,               INTENT(INOUT) :: RC             ! Success or failure?
!
! !REVISION HISTORY:
!  10 Jun 2015 - C. Keller - Initial version
!  20 Sep 2015 - C. Keller - Now allocate output array in this routine.
!EOP
!------------------------------------------------------------------------------
!BOC

    !======================================================================
    ! GetExtSpcVal_Int begins here
    !======================================================================

    ! Make sure output is properly allocated
    IF ( ALLOCATED(SpcScal) ) DEALLOCATE(SpcScal)
    ALLOCATE(SpcScal(NSPC))
    SpcScal=DefValue

    CALL GetExtSpcVal_Dr ( HcoConfig, ExtNr, NSPC, SpcNames, Prefix, RC, &
                           DefVal_IN=DefValue, SpcScal_IN=SpcScal )

    END SUBROUTINE GetExtSpcVal_Int
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: GetExtSpcVal_Char
!
! !DESCRIPTION: Subroutine GetExtSpcVal\_Char returns character values
! associated with the species for a given extension. Specifically, this routine
! searches for extension setting '<Prefix>\_SpecName' for every species passed
! through input argument SpcNames and writes those into output argument SpcScal.
! The default value DefValue is assigned to all elements of SpcScal with no
! corresponding extension setting.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE GetExtSpcVal_Char( HcoConfig, ExtNr,    NSPC,    SpcNames, &
                                Prefix,    DefValue, SpcScal, RC )
!
! !INPUT PARAMETERS:
!
    TYPE(ConfigObj),               POINTER       :: HcoConfig
    INTEGER,                       INTENT(IN   ) :: ExtNr          ! Extension Nr.
    INTEGER,                       INTENT(IN   ) :: NSPC           ! # of species
    CHARACTER(LEN=*),              INTENT(IN   ) :: SpcNames(NSPC) ! Species string
    CHARACTER(LEN=*),              INTENT(IN   ) :: Prefix         ! search prefix
    CHARACTER(LEN=*),              INTENT(IN   ) :: DefValue       ! default value
!
! !INPUT/OUTPUT PARAMETERS:
!
    CHARACTER(LEN=*), ALLOCATABLE, INTENT(INOUT) :: SpcScal(:)     ! Species scale factors
    INTEGER,                       INTENT(INOUT) :: RC             ! Success or failure?
!
! !REVISION HISTORY:
!  10 Jun 2015 - C. Keller - Initial version
!  20 Sep 2015 - C. Keller - Now allocate output array in this routine.
!EOP
!------------------------------------------------------------------------------
!BOC

    !======================================================================
    ! GetExtSpcVal_Char begins here
    !======================================================================

    ! Make sure output is properly allocated
    IF ( ALLOCATED(SpcScal) ) DEALLOCATE(SpcScal)
    ALLOCATE(SpcScal(NSPC))
    SpcScal=DefValue

    CALL GetExtSpcVal_Dr ( HcoConfig, ExtNr, NSPC, SpcNames, Prefix, RC, &
                           DefVal_Char=DefValue, SpcScal_Char=SpcScal )

    END SUBROUTINE GetExtSpcVal_char
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: GetExtSpcVal_Dr
!
! !DESCRIPTION: Subroutine GetExtSpcVal\_Dr is the GetExtSpcVal driver routine.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE GetExtSpcVal_Dr( HcoConfig, ExtNr, NSPC,    &
                              SpcNames,  Prefix, RC,     &
                              DefVal_SP, SpcScal_SP,     &
                              DefVal_Char, SpcScal_Char, &
                              DefVal_IN, SpcScal_IN       )
!
! !INPUT PARAMETERS:
!
    TYPE(ConfigObj),               POINTER                 :: HcoConfig
    INTEGER,                       INTENT(IN   )           :: ExtNr          ! Extension Nr.
    INTEGER,                       INTENT(IN   )           :: NSPC           ! # of species
    CHARACTER(LEN=*),              INTENT(IN   )           :: SpcNames(NSPC) ! Species string
    CHARACTER(LEN=*),              INTENT(IN   )           :: Prefix         ! search prefix
    REAL(sp),                      INTENT(IN   ), OPTIONAL :: DefVal_SP      ! default value
    INTEGER,                       INTENT(IN   ), OPTIONAL :: DefVal_IN      ! default value
    CHARACTER(LEN=*),              INTENT(IN   ), OPTIONAL :: DefVal_Char    ! default value
!
! !OUTPUT PARAMETERS:
!
    REAL(sp),                      INTENT(  OUT), OPTIONAL :: SpcScal_SP(NSPC)   ! Species values
    INTEGER,                       INTENT(  OUT), OPTIONAL :: SpcScal_IN(NSPC)   ! Species values
    CHARACTER(LEN=*),              INTENT(  OUT), OPTIONAL :: SpcScal_Char(NSPC) ! Species values
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,               INTENT(INOUT)           :: RC       ! Success or failure?
!
! !REVISION HISTORY:
!  10 Jun 2015 - C. Keller   - Initial version
!  20 Sep 2015 - C. Keller   - Options are now linked list.
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL ARGUMENTS:
!
    INTEGER              :: I
    LOGICAL              :: FND
    REAL(sp)             :: iScal_sp
    INTEGER              :: iScal_in
    CHARACTER(LEN=255)   :: iScal_char
    CHARACTER(LEN= 61)   :: IOptName
    CHARACTER(LEN=255)   :: MSG
    CHARACTER(LEN=255)   :: LOC = 'GetExtSpcVal_Dr (hco_extlist_mod.F90)'

    !======================================================================
    ! GetExtSpcVal_Dr begins here
    !======================================================================

    ! Do for every species
    DO I = 1, NSPC
       IOptName = TRIM(Prefix)//'_'//TRIM(SpcNames(I))

       IF ( PRESENT(SpcScal_sp) ) THEN
          CALL GetExtOpt ( HcoConfig, ExtNr, IOptName, OptValSp=iScal_sp, FOUND=FND, RC=RC )
          IF ( RC /= HCO_SUCCESS ) RETURN
          IF ( FND ) SpcScal_sp(I) = iScal_sp
       ENDIF
       IF ( PRESENT(SpcScal_in) ) THEN
          CALL GetExtOpt ( HcoConfig, ExtNr, IOptName, OptValInt=iScal_in, FOUND=FND, RC=RC )
          IF ( RC /= HCO_SUCCESS ) RETURN
          IF ( FND ) SpcScal_in(I) = iScal_in
       ENDIF
       IF ( PRESENT(SpcScal_char) ) THEN
          CALL GetExtOpt ( HcoConfig, ExtNr, IOptName, OptValChar=iScal_char, FOUND=FND, RC=RC )
          IF ( RC /= HCO_SUCCESS ) RETURN
          IF ( FND ) SpcScal_char(I) = iScal_char
       ENDIF
    ENDDO

    ! Return w/ success
    RC = HCO_SUCCESS

    END SUBROUTINE GetExtSpcVal_Dr
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: SetExtNr
!
! !DESCRIPTION: Subroutine SetExtNr overwrites the extension number of a
! given extension. The extension of interest is provided in argument
! ExtName. If this argument is omitted, the extension numbers of all
! extensions currently listed in ExtList will be set to the provided
! number. This is useful to disable all extensions by setting the ExtNr
! to a negative value.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE SetExtNr( HcoConfig, ExtNr, ExtName, RC )
!
! !USES:
!
    USE CHARPAK_MOD,  ONLY : TRANLC
!
! !INPUT PARAMETERS:
!
    TYPE(ConfigObj),  POINTER                 :: HcoConfig
    INTEGER,          INTENT(IN   )           :: ExtNr
    CHARACTER(LEN=*), INTENT(IN   ), OPTIONAL :: ExtName
!
! !INPUT/OUTPUT PARAMETER:
!
    INTEGER,          INTENT(INOUT)           :: RC
!
! !REVISION HISTORY:
!  12 Jan 2015 - C. Keller - Initial version
!  26 Oct 2016 - R. Yantosca - Don't nullify local ptrs in declaration stmts
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !INTERNAL VARIABLES:
!
    TYPE(Ext), POINTER  :: ThisExt
    CHARACTER(LEN=255)  :: LCname
    CHARACTER(LEN=255)  :: MSG
    LOGICAL             :: verb, overwrite

    !======================================================================
    ! SetExtNr begins here
    !======================================================================

    ! verbose?
    verb = HCO_IsVerb( HcoConfig%Err, 1 )

    ! Pass name to module and set to lower case
    IF ( PRESENT(ExtName) ) THEN
       LCname = TRIM(ExtName)
       CALL TRANLC( LCname )  ! lower case
    ELSE
       LCname = ''
    ENDIF

    ! Point to header of extensions list
    ThisExt => HcoConfig%ExtList

    ! Loop over all used extensions and check if any of them matches
    ! ExtName.
    DO WHILE ( ASSOCIATED ( ThisExt ) )

       ! Overwrite this ExtNr?
       overwrite = .FALSE.

       ! If argument ExtName is given, only overwrite extension number
       ! of that particular extension.
       IF ( PRESENT(ExtName) ) THEN
          IF ( TRIM(ThisExt%ExtName) == TRIM(LCname) ) overwrite = .TRUE.

       ! If argument is not given, overwrite all extensions except for
       ! HEMCO core
       ELSEIF ( ThisExt%ExtNr /= CoreNr ) THEN
          overwrite = .TRUE.
       ENDIF

       ! Overwrite extension number if needed
       IF ( overwrite ) THEN
          ThisExt%ExtNr = ExtNr
          IF ( verb ) THEN
             WRITE(MSG,*) 'Force ExtNr of extension ', TRIM(ThisExt%ExtName), &
                          ' to ', ExtNr
             CALL HCO_MSG(HcoConfig%Err,MSG)
          ENDIF
       ENDIF

       ! Advance to next extension
       ThisExt => ThisExt%NextExt
    ENDDO

    ! Cleanup
    ThisExt => NULL()

    ! Return w/ success
    RC = HCO_SUCCESS

  END SUBROUTINE SetExtNr
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ExtNrInUse
!
! !DESCRIPTION: Function ExtNrInUse checks if extension number ExtNr is
! in the list of used extensions or not.
!\\
!\\
! !INTERFACE:
!
  FUNCTION ExtNrInUse( ExtList, ExtNr ) Result ( InUse )
!
! !INPUT PARAMETERS:
!
    TYPE(Ext), POINTER       :: ExtList
    INTEGER,   INTENT(IN   ) :: ExtNr
!
! !RETURN VALUE::
!
    LOGICAL                  :: InUse
!
! !REVISION HISTORY:
!  03 Oct 2013 - C. Keller - Initial version
!  26 Oct 2016 - R. Yantosca - Don't nullify local ptrs in declaration stmts
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !INTERNAL VARIABLES:
!
    TYPE(Ext), POINTER  :: ThisExt
    CHARACTER(LEN=255)  :: LCname

    !======================================================================
    ! ExtNrInUse begins here
    !======================================================================

    ! Use number -999 for wildcard values
    IF ( ExtNr == -999 ) THEN
       InUse = .TRUE.
       RETURN
    ENDIF

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
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ExtFinal
!
! !DESCRIPTION: Function ExtFinal finalizes the extensions list.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE ExtFinal( ExtList )
!
! !INPUT/OUTPUT ARGUMENT:
!
    TYPE(Ext), POINTER       :: ExtList
!
! !REVISION HISTORY:
!  03 Oct 2013 - C. Keller - Initial version
!  20 Sep 2015 - C. Keller - Options are now linked list.
!  26 Oct 2016 - R. Yantosca - Don't nullify local ptrs in declaration stmts
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !INTERNAL VARIABLES:
!
    TYPE(Ext), POINTER  :: ThisExt
    TYPE(Ext), POINTER  :: NextExt

    !======================================================================
    ! ExtFinal begins here
    !======================================================================

    ! Point to header of extensions list
    ThisExt => ExtList
    NextExt => NULL()

    ! Loop over all extensions and deallocate the types
    DO WHILE ( ASSOCIATED ( ThisExt ) )

       ! First set pointer to next entry
       NextExt => ThisExt%NextExt

       ! Now clean up this entry
       ThisExt%NextExt => NULL()
       CALL HCO_CleanupOpt( ThisExt%Opts )
       DEALLOCATE ( ThisExt )

       ! Advance to next extension
       ThisExt => NextExt
    ENDDO

    ! Cleanup
    ThisExt => NULL()
    ExtList => NULL()

  END SUBROUTINE ExtFinal
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCO_AddOpt
!
! !DESCRIPTION: Subroutine HCO\_AddOpt adds a option name/value pair to the
! list of options.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCO_AddOpt ( HcoConfig, OptName, OptValue, ExtNr, RC, &
                          VERB,      IgnoreIfExist )
!
! !INPUT PARAMETERS:
!
    TYPE(ConfigObj),  POINTER                 :: HcoConfig      ! HEMCO config obj
    CHARACTER(LEN=*), INTENT(IN   )           :: OptName        ! OptName
    CHARACTER(LEN=*), INTENT(IN   )           :: OptValue       ! OptValue
    INTEGER,          INTENT(IN   )           :: ExtNr          ! Extension Nr.
    LOGICAL,          INTENT(IN   ), OPTIONAL :: VERB           ! Verbose on
    LOGICAL,          INTENT(IN   ), OPTIONAL :: IgnoreIfExist  ! Ignore if already exists
!
! !OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(INOUT)           :: RC             ! Return code
!
! !REVISION HISTORY:
!  18 Sep 2015 - C. Keller - Initial version
!  12 Dec 2015 - C. Keller - Added argument IgnoreIfExist
!  26 Oct 2016 - R. Yantosca - Don't nullify local ptrs in declaration stmts
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    TYPE(Ext), POINTER      :: ThisExt
    TYPE(Opt), POINTER      :: NewOpt
    CHARACTER(LEN=OPTLEN)   :: DUM
    LOGICAL                 :: Exists
    LOGICAL                 :: VRB
    LOGICAL                 :: Ignore
    CHARACTER(LEN=255)      :: MSG
    CHARACTER(LEN=255)      :: LOC = 'HCO_AddOpt (hco_chartools_mod.F90)'

    !=================================================================
    ! HCO_AddOpt begins here!
    !=================================================================

    ! Nullify
    ThisExt => NULL()
    NewOpt  => NULL()

    ! Init optional variables
    VRB    = .TRUE.
    Ignore = .FALSE.
    IF ( PRESENT(VERB)          ) VRB    = VERB
    IF ( PRESENT(IgnoreIfExist) ) Ignore = IgnoreIfExist

    ! Check if this option already exists
    DUM = HCO_GetOpt( HcoConfig%ExtList, OptName, ExtNr=ExtNr )

    ! If option already exists...
    IF ( TRIM(DUM) /= TRIM(EMPTYOPT) ) THEN

       ! Can leave here if we shall ignore the option if it already exists
       IF ( Ignore ) THEN

          ! If option exists and is the same, nothing to do
          IF ( TRIM(DUM) /= ADJUSTL(TRIM(OptValue)) ) THEN
             WRITE(*,*) 'Option is already defined - use original value of ', &
                        TRIM(DUM), ' and ignore the following value: ',       &
                        TRIM(OptName), ': ', TRIM(OptValue)
          ENDIF
          RC = HCO_SUCCESS
          RETURN

       ! If ignore flag is false:
       ELSE
          ! Error if values are not the same
          IF ( TRIM(DUM) /= ADJUSTL(TRIM(OptValue)) ) THEN
             MSG = 'Cannot add option pair: '//TRIM(OptName)//': '//TRIM(OptValue) &
                // ' - option already exists: '//TRIM(OptName)//': '//TRIM(DUM)
             CALL HCO_ERROR ( HcoConfig%Err, MSG, RC, THISLOC=LOC )
             RETURN
          ! Return with no error if values are the same
          ELSE
             RC = HCO_SUCCESS
             RETURN
          ENDIF
       ENDIF
    ENDIF

    ! Find extension of interest
    ThisExt => HcoConfig%ExtList
    DO WHILE ( ASSOCIATED ( ThisExt ) )
       IF ( ThisExt%ExtNr == ExtNr ) EXIT
       ThisExt => ThisExt%NextExt
    ENDDO

    IF ( .NOT. ASSOCIATED( ThisExt ) ) THEN
       WRITE(MSG,*) 'Cannot add option to extension Nr. ', ExtNr
       MSG = TRIM(MSG) // '. Make sure this extension is activated!'
       CALL HCO_ERROR(HcoConfig%Err,MSG,RC,THISLOC='AddOpt (hco_extlist_mod)')
       RETURN
    ENDIF

    ! Create new option
    ALLOCATE(NewOpt)
    NewOpt%OptName  = ADJUSTL( TRIM(OptName ) )
    NewOpt%OptValue = ADJUSTL( TRIM(OptValue) )

    ! Add to option linked list
    IF ( ASSOCIATED(ThisExt%Opts) ) THEN
       NewOpt%NextOpt => ThisExt%Opts
    ELSE
       NewOpt%NextOpt => NULL()
    ENDIF
    ThisExt%Opts => NewOpt

    ! Verbose
    IF ( VRB .AND. HcoConfig%amIRoot .AND. HCO_IsVerb(HcoConfig%Err,2) ) THEN
       MSG = 'Added the following option: ' // TRIM(OptName)//': '//TRIM(OptValue)
       CALL HCO_MSG(HcoConfig%Err,MSG)
    ENDIF

    ! Cleanup and return w/ success
    ThisExt => NULL()
    RC      =  HCO_SUCCESS

  END SUBROUTINE HCO_AddOpt
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCO_GetOpt
!
! !DESCRIPTION: Subroutine HCO\_GetOpt returns a option value for the given
! option name.
!\\
!\\
! !INTERFACE:
!
  FUNCTION HCO_GetOpt ( ExtList, OptName, ExtNr ) RESULT ( OptValue )
!
! !INPUT PARAMETERS:
!
    TYPE(Ext),        POINTER                 :: ExtList        ! Extension list
    CHARACTER(LEN=*), INTENT(IN   )           :: OptName  ! OptName
    INTEGER,          INTENT(IN   ), OPTIONAL :: ExtNr      ! Extension Nr.
!
! !OUTPUT PARAMETERS:
!
    CHARACTER(LEN=OPTLEN)                     :: OptValue ! OptValue
!
! !REVISION HISTORY:
!  18 Sep 2015 - C. Keller - Initial version
!  26 Oct 2016 - R. Yantosca - Don't nullify local ptrs in declaration stmts
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER               :: ThisExtNr
    LOGICAL               :: OptFound
    TYPE(Opt), POINTER    :: ThisOpt
    TYPE(Ext), POINTER    :: ThisExt

    !=================================================================
    ! HCO_GetOpt begins here!
    !=================================================================

    ! Init
    OptValue = EMPTYOPT
    OptFound = .FALSE.
    ThisOpt  => NULL()
    ThisExt  => NULL()

    ! Extension number to search for. If not explicitly set through the
    ! input argument, set to -999 to search all extensions.
    IF ( PRESENT(ExtNr) ) THEN
       ThisExtNr = ExtNr
    ELSE
       ThisExtNr = -999
    ENDIF

    ! Find extension of interest
    ThisExt => ExtList
    DO WHILE ( ASSOCIATED ( ThisExt ) )

       ! Check if this is the extension of interest. If extension number
       ! is set to -999, scan through all extensions.
       IF ( ThisExtNr /= -999 .AND. ThisExt%ExtNr /= ThisExtNr ) THEN
          ThisExt => ThisExt%NextExt
          CYCLE
       ENDIF

       ! Walk through token list until we find the given value
       ThisOpt => ThisExt%Opts
       DO WHILE ( ASSOCIATED(ThisOpt) )

          ! Check if this is the token of interest
          IF ( TRIM(ThisOpt%OptName) == ADJUSTL(TRIM(OptName)) ) THEN
             OptValue = ADJUSTL( TRIM(ThisOpt%OptValue) )
             OptFound = .TRUE.
             EXIT
          ENDIF

          ! Advance in list
          ThisOpt => ThisOpt%NextOpt
       END DO

       ! Advance to next extension
       IF ( OptFound ) THEN
          ThisExt => NULL()
       ELSE
          ThisExt => ThisExt%NextExt
       ENDIF
    ENDDO

    ! Free pointer
    ThisOpt => NULL()
    ThisExt => NULL()

  END FUNCTION HCO_GetOpt
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCO_ROOT
!
! !DESCRIPTION: Function HCO\_ROOT returns the root character string. This is
! a wrapper routine equivalent to HCO\_GetOpt('ROOT'). Since the ROOT character
! is called very frequently, it is recommended to use this routine instead.
!\\
!\\
! !INTERFACE:
!
  FUNCTION HCO_ROOT ( HcoConfig ) RESULT ( OutRoot )
!
! !INPUT PARAMETERS:
!
!
! !OUTPUT PARAMETERS:
!
    TYPE(ConfigObj), POINTER                  :: HcoConfig
    CHARACTER(LEN=OPTLEN)                     :: OutRoot ! Root output
!
! !REVISION HISTORY:
!  18 Sep 2015 - C. Keller - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC

    OutRoot = HcoConfig%ROOT

  END FUNCTION HCO_ROOT
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCO_CleanupOpt
!
! !DESCRIPTION: Subroutine HCO\_CleanupOpt cleans up the given options linked
!  list.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCO_CleanupOpt ( OptList )
!
! !INPUT PARAMETERS:
!
!
! !OUTPUT PARAMETERS:
!
    TYPE(Opt), POINTER        :: OptList
!
! !REVISION HISTORY:
!  18 Sep 2015 - C. Keller - Initial version
!  26 Oct 2016 - R. Yantosca - Don't nullify local ptrs in declaration stmts
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    TYPE(Opt), POINTER    :: ThisOpt
    TYPE(Opt), POINTER    :: NextOpt

    !=================================================================
    ! HCO_CleanupOpt begins here!
    !=================================================================

    ! Walk through option list until we find the given value
    NextOpt => NULL()
    ThisOpt => OptList
    DO WHILE ( ASSOCIATED(ThisOpt) )

       ! Archive next option in list
       NextOpt => ThisOpt%NextOpt

       ! Cleanup option
       ThisOpt%NextOpt => NULL()
       NULLIFY(ThisOpt)

       ! Go to next option in list (previously archived)
       ThisOpt => NextOpt
    END DO

    ! Free pointer
    ThisOpt => NULL()
    NextOpt => NULL()

  END SUBROUTINE HCO_CleanupOpt
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCO_SetDefaultToken
!
! !DESCRIPTION: Subroutine HCO\_SetDefaultToken is a wrapper routine to
! initialize the default set of HEMCO tokens. These can be obtained at any
! place in the HEMCO code via subroutine HCO\_GetOpt, e.g. HCO\_GetOpt('RES')
! will return the 'RES' token.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCO_SetDefaultToken( CF, RC )
!
! !USES:
!
!
! !INPUT PARAMETERS:
!
    TYPE(ConfigObj),  POINTER                 :: CF         ! Config object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(INOUT)           :: RC         ! Return code
!
! !REVISION HISTORY:
!  18 Sep 2015 - C. Keller - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    CHARACTER(LEN=OPTLEN)     :: DUM
    LOGICAL                   :: FOUND

    !=================================================================
    ! HCO_SetDefaultToken begins here!
    !=================================================================

    IF ( Trim(CF%MetField) == 'GEOSFP' ) THEN
       DEF_MET_UC = 'GEOSFP'
       DEF_MET_LC = 'geosfp'
       DEF_CN_YR  = '2011'  ! Constant met fld year
       DEF_NC_VER = 'nc'    ! NetCDF extension
    ELSE IF ( TRIM(CF%MetField) == 'MERRA2' ) THEN
       DEF_MET_UC = 'MERRA2'
       DEF_MET_LC = 'merra2'
       DEF_CN_YR  = '2015'  ! Constant met fld year
       DEF_NC_VER = 'nc4'   ! NetCDF extension
    ENDIF

    IF ( TRIM(CF%GridRes) == '4.0x5.0' ) THEN
       DEF_RES = '4x5'
    ELSE IF ( TRIM(CF%GridRes) == '2.0x2.5' ) THEN
       DEF_RES = '2x25'
    ELSE IF ( TRIM(CF%GridRes) == '0.5x0.625' ) THEN
       DEF_RES = '05x0625'
    ELSE IF ( TRIM(CF%GridRes) == '0.25x0.3125' ) THEN
       DEF_RES = '025x03125'
    ENDIF

    ! Wildcard character
    CALL GetExtOpt( CF, CoreNr, 'Wildcard', OptValChar=DUM, Found=FOUND, RC=RC )
    IF ( RC /= HCO_SUCCESS ) RETURN
    IF ( .NOT. FOUND) DUM = DEF_WILDCARD
    CALL HCO_AddOpt( CF, 'Wildcard', DUM, CoreNr, RC, VERB=.FALSE. )

    ! Separator
    CALL GetExtOpt( CF, CoreNr, 'Separator', OptValChar=DUM, Found=FOUND, RC=RC )
    IF ( RC /= HCO_SUCCESS ) RETURN
    IF ( .NOT. FOUND) DUM = DEF_SEPARATOR
    CALL HCO_AddOpt( CF, 'Separator', DUM, CoreNr, RC, VERB=.FALSE. )

    ! Colon
    CALL GetExtOpt( CF, CoreNr, 'Colon', OptValChar=DUM, Found=FOUND, RC=RC )
    IF ( RC /= HCO_SUCCESS ) RETURN
    IF ( .NOT. FOUND) DUM = DEF_COLON
    CALL HCO_AddOpt( CF, 'Colon', DUM, CoreNr, RC, VERB=.FALSE. )

    ! Root directory
    CALL GetExtOpt( CF, CoreNr, 'ROOT', OptValChar=DUM, Found=FOUND, RC=RC )
    IF ( RC /= HCO_SUCCESS ) RETURN
    IF ( .NOT. FOUND) DUM = DEF_ROOT
    CALL HCO_AddOpt( CF, 'ROOT', DUM, CoreNr, RC, VERB=.FALSE. )

    ! Also save in local variable (for fast access via HCO_ROOT)
    CF%ROOT = ADJUSTL( TRIM(DUM) )

    ! Meteorology token (uppercase)
    CALL GetExtOpt( CF, CoreNr, 'MET', OptValChar=DUM, Found=FOUND, RC=RC )
    IF ( RC /= HCO_SUCCESS ) RETURN
    IF ( .NOT. FOUND) DUM = DEF_MET_UC
    CALL HCO_AddOpt( CF, 'MET', DUM, CoreNr, RC, VERB=.FALSE. )

    ! Meteorology token (lowercase)
    CALL GetExtOpt( CF, CoreNr, 'met', OptValChar=DUM, Found=FOUND, RC=RC )
    IF ( RC /= HCO_SUCCESS ) RETURN
    IF ( .NOT. FOUND) DUM = DEF_MET_LC
    CALL HCO_AddOpt( CF, 'met', DUM, CoreNr, RC, VERB=.FALSE. )

    ! Year for constant met fields
    CALL GetExtOpt( CF, CoreNr, 'CNYR', OptValChar=DUM, Found=FOUND, RC=RC )
    IF ( RC /= HCO_SUCCESS ) RETURN
    IF ( .NOT. FOUND) DUM = DEF_CN_YR
    CALL HCO_AddOpt( CF, 'CNYR', DUM, CoreNr, RC, VERB=.FALSE. )

    ! NetCDF version extension
    CALL GetExtOpt( CF, CoreNr, 'NC', OptValChar=DUM, Found=FOUND, RC=RC )
    IF ( RC /= HCO_SUCCESS ) RETURN
    IF ( .NOT. FOUND) DUM = DEF_NC_VER
    CALL HCO_AddOpt( CF, 'NC', DUM, CoreNr, RC, VERB=.FALSE. )

    ! Resolution token
    CALL GetExtOpt( CF, CoreNr, 'RES', OptValChar=DUM, Found=FOUND, RC=RC )
    IF ( RC /= HCO_SUCCESS ) RETURN
    IF ( .NOT. FOUND ) DUM = DEF_RES
    CALL HCO_AddOpt( CF, 'RES', DUM, CoreNr, RC, VERB=.FALSE. )

    ! Return w/ success
    RC =  HCO_SUCCESS

  END SUBROUTINE HCO_SetDefaultToken
!EOC
END MODULE HCO_ExtList_Mod
