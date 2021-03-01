!------------------------------------------------------------------------------
! QFYAML: Bob Yantosca | yantosca@seas.harvard.edu | Apr 2020
! Based on existing package https://github.com/jannisteunissen/config_fortran
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: yaml_mod.F90
!
! !DESCRIPTION: Contains routines for reading a YAML file into Fortran,
!  based off the "config_fortran" package of H. J. Teunissen.
!\\
!\\
! !INTERFACE:
!
MODULE QFYAML_Mod
!
! !USES:
!
  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC TYPES:
!
  PUBLIC :: yp
  PUBLIC :: QFYAML_t
!
! !PUBLIC DATA MEMBERS:
!
  ! Constants
  PUBLIC :: QFYAML_Success
  PUBLIC :: QFYAML_Failure
  PUBLIC :: QFYAML_NamLen
  PUBLIC :: QFYAML_StrLen
  PUBLIC :: QFYAML_MaxArr

  ! Public methods
  PUBLIC :: QFYAML_Add
  PUBLIC :: QFYAML_Add_Get
  PUBLIC :: QFYAML_CleanUp
  PUBLIC :: QFYAML_Get
  PUBLIC :: QFYAML_Check
  PUBLIC :: QFYAML_Init
  PUBLIC :: QFYAML_Merge
  PUBLIC :: QFYAML_Update
!
! !REMARKS:
!  QFYAML -- The Quick Fortran YAML parser!
!
!  I developed this package because I needed a quick-and-dirty YAML parser
!  in Fortran for reading in certain YAML files (e.g. species database)
!  into the GEOS-Chem model. I found that certain Fortran YAML parsers
!  either did not support mapping, or required the bleeding edge versions
!  of Fortran compilers.
!
!  The back end is code that was taken from the "config_fortran" package
!  (https://github.com/jannisteunissen/config_fortran) by H. J. Teunissen.
!  and subsequently modified by myse;lf
!
!  The front end (parser) has been modified to accept YAML format instead of
!  configuration file format.  Not all features of YAML have been implemented.
!  At present, I have only tested with YAML mappings but as time allows I
!  can try to add other YAML features.
!
!  At present, nested levels of variables are not supported, but
!  might bein the future.
!
!  I have removed some routines that are not as pertinent to YAML input
!  from the original config-fortran code.
!
!      -- Bob Yantosca (15 Apr 2020), yantosca@seas.harvard.edu
!
! !REVISION HISTORY:
!  15 Apr 2020 - R. Yantosca - Initial version
!  See the subsequent Git history with the gitk browser!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
  ! The precision kind-parameter (4-byte)
  INTEGER, PARAMETER :: yp                    = KIND( 0.0e0 )

  ! Success/failure return codes
  INTEGER, PARAMETER :: QFYAML_Success        =  0
  INTEGER, PARAMETER :: QFYAML_Failure        = -1

  ! Numeric type constants
  INTEGER, PARAMETER :: QFYAML_num_types      =  4
  INTEGER, PARAMETER :: QFYAML_integer_type   =  1
  INTEGER, PARAMETER :: QFYAML_real_type      =  2
  INTEGER, PARAMETER :: QFYAML_string_type    =  3
  INTEGER, PARAMETER :: QFYAML_bool_type      =  4
  INTEGER, PARAMETER :: QFYAML_unknown_type   =  0

  ! How was the YAML file set?
  INTEGER, PARAMETER :: QFYAML_set_by_default =  1
  INTEGER, PARAMETER :: QFYAML_set_by_file    =  3

  ! Maxima
  INTEGER, PARAMETER :: QFYAML_NamLen         =  80    ! Max len for names
  INTEGER, PARAMETER :: QFYAML_StrLen         =  255   ! Max len for strings
  INTEGER, PARAMETER :: QFYAML_MaxArr         =  1000  ! Max entries per array

  CHARACTER(LEN=7), PARAMETER :: QFYAML_type_names(0:QFYAML_num_types) = &
      (/ 'storage', 'integer', 'real   ', 'string ', 'bool   ' /)

  ! The separator(s) for array-like variables (space, comma, ', ", and tab)
  CHARACTER,         PARAMETER :: tab_char = char(9)
  CHARACTER(LEN=*),  PARAMETER :: QFYAML_separators = " ,'"""//tab_char

  ! Bracket characters
  CHARACTER(LEN=4),  PARAMETER :: QFYAML_brackets = "{}[]"

  ! The separator for categories (stored in var_name)
  CHARACTER(LEN=1),  PARAMETER :: QFYAML_category_separator = "%"

  ! The default string for data that is not yet stored
  CHARACTER(LEN=21), PARAMETER :: unstored_data_string="__UNSTORED_DATA_STRING"

  ! Type for a single variable
  TYPE, PRIVATE :: QFYAML_var_t
     PRIVATE
     CHARACTER(LEN=QFYAML_namlen)              :: category
     CHARACTER(LEN=QFYAML_namlen)              :: var_name
     CHARACTER(LEN=QFYAML_strlen)              :: description
     INTEGER                                   :: var_type
     INTEGER                                   :: var_size
     LOGICAL                                   :: dynamic_size
     LOGICAL                                   :: used
     INTEGER                                   :: set_by=QFYAML_set_by_default
     CHARACTER(LEN=QFYAML_strlen)              :: stored_data
     CHARACTER(LEN=QFYAML_namlen)              :: anchor_ptr
     CHARACTER(LEN=QFYAML_namlen)              :: anchor_tgt
     REAL(yp),                     ALLOCATABLE :: real_data(:)
     INTEGER,                      ALLOCATABLE :: int_data(:)
     CHARACTER(LEN=QFYAML_strlen), ALLOCATABLE :: char_data(:)
     LOGICAL,                      ALLOCATABLE :: bool_data(:)
  END TYPE QFYAML_var_t

  ! Type for the list of variables
  TYPE :: QFYAML_t
     LOGICAL                                   :: sorted = .false.
     INTEGER                                   :: num_vars = 0
     TYPE(QFYAML_var_t),           ALLOCATABLE :: vars(:)
  END TYPE QFYAML_t

  ! Interface to add variables to the configuration
  INTERFACE QFYAML_Add
     MODULE PROCEDURE  Add_Real,       Add_Real_Array
     MODULE PROCEDURE  Add_Int,        Add_Int_Array
     MODULE PROCEDURE  Add_String,     Add_String_Array
     MODULE PROCEDURE  Add_Bool,       Add_Bool_Array
  END INTERFACE QFYAML_Add

  ! Interface to get variables from the configuration
  INTERFACE QFYAML_Get
     MODULE PROCEDURE  Get_Real,       Get_Real_Array
     MODULE PROCEDURE  Get_Int,        Get_Int_Array
     MODULE PROCEDURE  Get_Bool,       Get_Bool_Array
     MODULE PROCEDURE  Get_String,     Get_String_Array
  END INTERFACE QFYAML_Get

  ! Interface to get variables from the configuration
  INTERFACE QFYAML_Add_Get
     MODULE PROCEDURE  Add_Get_Real,   Add_Get_Real_Array
     MODULE PROCEDURE  Add_Get_Int,    Add_Get_Int_Array
     MODULE PROCEDURE  Add_Get_Bool,   Add_Get_Bool_Array
     MODULE PROCEDURE  Add_Get_String, Add_Get_String_Array
  END INTERFACE QFYAML_Add_Get

  ! Interface to get variables from the configuration
  INTERFACE QFYAML_Update
     MODULE PROCEDURE  Update_Real,    Update_Real_Array
     MODULE PROCEDURE  Update_Int,     Update_Int_Array
     MODULE PROCEDURE  Update_Bool,    Update_Bool_Array
     MODULE PROCEDURE  Update_String,  Update_String_Array
  END INTERFACE QFYAML_Update

CONTAINS
!EOC
!------------------------------------------------------------------------------
! QFYAML: Bob Yantosca | yantosca@seas.harvard.edu | Apr 2020
! Based on existing package https://github.com/jannisteunissen/config_fortran
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Handle_Error
!
! !DESCRIPTION: This routine will be called if an error occurs in one of
!  the subroutines of this module.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Handle_Error( errMsg, RC, thisLoc )
!
! !INPUT PARAMETERS:
!
    CHARACTER(LEN=*), INTENT(IN)  :: errMsg    ! Error message
    CHARACTER(LEN=*), OPTIONAL    :: thisLoc
!
! INPUT/OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(OUT) :: RC        ! Return code
!
! !REVISION HISTORY:
!  15 Apr 2020 - R. Yantosca - Initial version
!  See the subsequent Git history with the gitk browser!
!EOP
!------------------------------------------------------------------------------
!BOC
    !=======================================================================
    ! Handle_Error begins here!
    !=======================================================================
    WRITE( 6, "(a)" ) REPEAT( "=", 79 )
    WRITE( 6, "(a)" ) "QFYAML ERROR: " // TRIM( errMsg )
    IF ( PRESENT( thisLoc ) ) WRITE( 6, '(a)' ) TRIM( thisLoc )
    WRITE( 6, "(a)" ) REPEAT( "=", 79 )
    WRITE( 6, "(a)" )

    ! Return failure
    RC = QFYAML_FAILURE

  END SUBROUTINE Handle_Error
!EOC
!------------------------------------------------------------------------------
! QFYAML: Bob Yantosca | yantosca@seas.harvard.edu | Apr 2020
! Based on existing package https://github.com/jannisteunissen/config_fortran
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get_Var_Index
!
! !DESCRIPTION:  Return the index of the variable with name 'var_name',
!  or -1 if not found.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Get_Var_Index( yml, var_name, ix )
!
! !INPUT PARAMETERS:
!
    TYPE(QFYAML_t),   INTENT(IN)  :: yml
    CHARACTER(LEN=*), INTENT(IN)  :: var_name
!
! !OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(OUT) :: ix
!
! !REVISION HISTORY:
!  15 Apr 2020 - R. Yantosca - Initial version
!  See the subsequent Git history with the gitk browser!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER :: i

    !=======================================================================
    ! Get_Var_Index begins here!
    !=======================================================================

    ! Initialize
    ix = -1

    IF ( yml%sorted ) THEN

       ! If the variable names have been sorted, use a binary search
       CALL Binary_Search_Variable( yml, var_name, ix )

    ELSE

       ! Otherwise use a linear search
       DO i = 1, yml%num_vars
          IF ( TRIM( yml%vars(i)%var_name ) == TRIM( var_name ) ) THEN
             ix = i
             EXIT
          ENDIF
       ENDDO

    ENDIF

  END SUBROUTINE Get_Var_Index
!EOC
!------------------------------------------------------------------------------
! QFYAML: Bob Yantosca | yantosca@seas.harvard.edu | Apr 2020
! Based on existing package https://github.com/jannisteunissen/config_fortran
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get_Anchor_Info
!
! !DESCRIPTION: Returns information about a variable containing an anchor
!  target field: the index, the category name, and the variable name (minus
!  the category).  Missing values are returned if not found.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Get_Anchor_Info( yml, anchor_ptr, begin_ix, end_ix, anchor_cat )
!
! !INPUT PARAMETERS:
!
    TYPE(QFYAML_t),               INTENT(IN)  :: yml        ! Config object
    CHARACTER(LEN=*),             INTENT(IN)  :: anchor_ptr ! Anchor to match
!
! !OUTPUT PARAMETERS:
!
    INTEGER,                      INTENT(OUT) :: begin_ix   ! 1st var w/ anchor
    INTEGER,                      INTENT(OUT) :: end_ix     ! last var w/ anchor
    CHARACTER(LEN=QFYAML_namlen), INTENT(OUT) :: anchor_cat ! Anchor category
!
! !REVISION HISTORY:
!  15 Apr 2020 - R. Yantosca - Initial version
!  See the subsequent Git history with the gitk browser!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER :: i

    !=======================================================================
    ! Get_Anchor_Info begins here!
    !=======================================================================

    ! Initialize
    begin_ix   = 0
    end_ix     = 0
    anchor_cat = "UNKNOWN"

    ! Linear search
    DO i = 1, yml%num_vars
       IF ( TRIM( yml%vars(i)%anchor_tgt ) == TRIM( anchor_ptr ) ) THEN
          IF ( begin_ix == 0 ) begin_ix = i
          end_ix = i
       ENDIF
    ENDDO

    ! Also return the category for this anchor
    IF ( begin_ix > 0 ) THEN
       anchor_cat = yml%vars(begin_ix)%category
    ENDIF

  END SUBROUTINE Get_Anchor_Info
!EOC
!------------------------------------------------------------------------------
! QFYAML: Bob Yantosca | yantosca@seas.harvard.edu | Apr 2020
! Based on existing package https://github.com/jannisteunissen/config_fortran
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: QFYAML_Init
!
! !DESCRIPTION: Initializes a QFYAML_t object from a YAML file.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE QFYAML_Init( fileName, yml, yml_anchored, RC )
!
! !INPUT PARAMETERS:
!
    CHARACTER(LEN=*), INTENT(IN)    :: fileName
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(QFYAML_t),   INTENT(INOUT) :: yml
    TYPE(QFYAML_t),   INTENT(INOUT) :: yml_anchored
!
! !OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(OUT)   :: RC
!
! !REVISION HISTORY:
!  15 Apr 2020 - R. Yantosca - Initial version
!  See the subsequent Git history with the gitk browser!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Strings
    CHARACTER(LEN=QFYAML_strlen) :: errMsg, thisLoc

    !=======================================================================
    ! QFYAML_Init begins here!
    !=======================================================================

    ! Initialize
    RC      = QFYAML_success
    errMsg  = ''
    thisLoc = ' -> at QFYAML_Init (in module qfyaml_mod.F90)'

    ! Read the YML file
    CALL QFYAML_Read_File( yml, fileName, yml_anchored, RC )

    ! Trap potential errors
    IF ( RC /= QFYAML_Success ) THEN
       errMsg = 'Error encountered in "QFYAML_Read_File"!'
       CALL Handle_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    ! Sort the variable names in the yml object for faster search
    CALL QFYAML_Sort( yml )

  END SUBROUTINE QFYAML_Init
!EOC
!------------------------------------------------------------------------------
! QFYAML: Bob Yantosca | yantosca@seas.harvard.edu | Apr 2020
! Based on existing package https://github.com/jannisteunissen/config_fortran
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: QFYAML_Merge
!
! !DESCRIPTION: Concatetenates two QFYAML_t objects together.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE QFYAML_Merge( yml1, yml2, yml, RC )
!
! !INPUT PARAMETERS:
!
    TYPE(QFYAML_t),   INTENT(IN)  :: yml1
    TYPE(QFYAML_t),   INTENT(IN)  :: yml2
!
! !OUTPUT PARAMETERS:
!
    TYPE(QFYAML_t),   INTENT(OUT) :: yml
    INTEGER,          INTENT(OUT) :: RC
!
! !REVISION HISTORY:
!  15 Apr 2020 - R. Yantosca - Initial version
!  See the subsequent Git history with the gitk browser!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER                      :: N

    ! Strings
    CHARACTER(LEN=QFYAML_strlen) :: errMsg
    CHARACTER(LEN=QFYAML_strlen) :: thisLoc

    !=======================================================================
    ! QFYAML_Init begins here!
    !=======================================================================

    ! Initialize
    RC      = QFYAML_success
    N       = 0
    errMsg  = ''
    thisLoc = ' -> at QFYAML_Merge (in module qfyaml_mod.F90)'

    ! Total number of variables
    yml%num_vars = yml1%num_vars + yml2%num_vars

    ! Allocate yml%vars
    ALLOCATE( yml%vars( yml%num_vars ), STAT=RC )
    IF ( RC /= QFYAML_Success ) THEN
       errMsg = 'Could not allocate the yml%vars object!"!'
       CALL Handle_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    ! Add the variables from the first object
    DO N = 1, yml1%num_vars
       yml%vars(N) = yml1%vars(N)
    ENDDO

    ! Add the variables from the second object
    DO N = 1, yml2%num_vars
       yml%vars(N + yml1%num_vars) = yml2%vars(N)
    ENDDO

    ! Sort the variable names in the yml object for faster search
    CALL QFYAML_Sort( yml )

  END SUBROUTINE QFYAML_Merge
!EOC
!------------------------------------------------------------------------------
! QFYAML: Bob Yantosca | yantosca@seas.harvard.edu | Apr 2020
! Based on existing package https://github.com/jannisteunissen/config_fortran
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: QFYAML_Read_File
!
! !DESCRIPTION: Read variables from a YAML file into a configuration object.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE QFYAML_Read_File( yml, fileName, yml_anchored, RC )
!
! !INPUT PARAMETERS:
!
    CHARACTER(LEN=*), INTENT(IN)    :: fileName      ! YAML file to read
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(QFYAML_t),   INTENT(INOUT) :: yml           ! Configuration object
    TYPE(QFYAML_t),   INTENT(INOUT) :: yml_anchored  ! "" for anchored vars
!
! !OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(OUT  ) :: RC
!
! !REVISION HISTORY:
!  15 Apr 2020 - R. Yantosca - Initial version
!  See the subsequent Git history with the gitk browser!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    LOGICAL                      :: valid_syntax
    INTEGER                      :: anchor_ix
    INTEGER                      :: begin_ix
    INTEGER                      :: end_ix
    INTEGER                      :: I
    INTEGER                      :: io_state
    INTEGER                      :: N
    INTEGER                      :: line_number
    INTEGER                      :: my_unit

    ! Strings
    CHARACTER(LEN=QFYAML_namlen) :: line_fmt
    CHARACTER(LEN=QFYAML_namlen) :: category
    CHARACTER(LEN=QFYAML_namlen) :: anchor_cat
    CHARACTER(LEN=QFYAML_namlen) :: anchor_ptr
    CHARACTER(LEN=QFYAML_namlen) :: anchor_tgt
    CHARACTER(LEN=QFYAML_namlen) :: var_pt_to_anchor
    CHARACTER(LEN=QFYAML_namlen) :: var_w_anchor
    CHARACTER(LEN=QFYAML_namlen) :: var_name
    CHARACTER(LEN=QFYAML_strlen) :: errMsg
    CHARACTER(LEN=QFYAML_strlen) :: line
    CHARACTER(LEN=QFYAML_strlen) :: thisLoc

    !=======================================================================
    ! QFYAML_READ_FILE begins here!
    !=======================================================================

    ! Initialize
    RC           = QFYAML_Success
    anchor_ptr   = ""
    anchor_tgt   = ""
    category     = ""
    errMsg       = ""
    thisLoc      = " -> at QFYAML_Read_File (in module qfyaml_mod.F90)"
    line_number  = 0
    my_unit      = 777

    !=======================================================================
    ! First pass: read the file
    !=======================================================================

    ! Open the file
    OPEN( my_unit, FILE=TRIM(filename), STATUS="old", ACTION="read", IOSTAT=RC)
    IF ( RC /= QFYAML_SUCCESS ) THEN
       errMsg = 'Could not open file: ' // TRIM( fileName )
       CALL Handle_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    ! Create format line
    WRITE( line_fmt, "(a,i0,a)") "(a", QFYAML_strlen, ")"

    ! Start looping
    DO

       ! Read each line and increment the count
       READ( my_unit, FMT=TRIM(line_fmt), ERR=998, end=999) line
       line_number = line_number + 1

       ! Parse each line for information.  This will also add
       ! each found variable to the "yml" configuration object.
       ! YAML anchors will also be referenced.
       CALL Parse_Line( yml          = yml,                                  &
                        yml_anchored = yml_anchored,                         &
                        set_by       = QFYAML_set_by_file,                   &
                        line_arg     = line,                                 &
                        valid_syntax = valid_syntax,                         &
                        category     = category,                             &
                        anchor_tgt   = anchor_tgt,                           &
                        anchor_ptr   = anchor_ptr,                           &
                        RC           = RC                                   )

       ! Trap potential errors
       IF ( .not. valid_syntax ) THEN
          WRITE( errMsg, *) "Cannot read line ", line_number, &
               " from ", TRIM(filename)
          CALL Handle_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF
    ENDDO

    ! Error handling
998 WRITE( errMsg, "(a,i0,a,i0)" ) " IOSTAT = ", io_state, &
         " while reading from " // trim(filename) // " at line ", &
         line_number
    CALL Handle_Error( errMsg, RC, thisLoc )
    RETURN

    ! Routine ends here if the end of "filename" is reached
999 CLOSE( my_unit, iostat=io_state )

    !=======================================================================
    ! Second pass: Create variables that point to YAML anchors
    ! with all of the corresponding properties.
    !=======================================================================
    DO N = 1, yml_anchored%num_vars

       ! Get proprties of each variable that points to an anchor
       anchor_ptr = yml_anchored%vars(N)%anchor_ptr
       category   = yml_anchored%vars(N)%category
       var_name   = yml_anchored%vars(N)%var_name

       ! Find all target variables with the given value of anchor_ptr,
       ! and return the start and ending indices in the yml config object.
       CALL Get_Anchor_Info( yml        = yml,                                &
                             anchor_ptr = anchor_ptr,                         &
                             begin_ix   = begin_ix,                           &
                             end_ix     = end_ix,                             &
                             anchor_cat = anchor_cat )

       ! Loop over all target variables containing the value of anchor_ptr
       DO anchor_ix = begin_ix, end_ix

          ! Variable with the anchor
          var_w_anchor = yml%vars(anchor_ix)%var_name

          ! Variable that we want to point to the anchor
          I = INDEX( var_w_anchor, QFYAML_category_separator )
          var_pt_to_anchor = TRIM( category )                             // &
                             QFYAML_category_separator                    // &
                             var_w_anchor(I+1:)

          ! Create a new variable for this category,
          ! copying the fields of the variable with the anchor.
          CALL Copy_Anchor_Variable( yml              = yml,                 &
                                     anchor_ix        = anchor_ix,           &
                                     var_w_anchor     = var_w_anchor,        &
                                     var_pt_to_anchor = var_pt_to_anchor,    &
                                     RC               = RC                  )

          ! Trap potential errors
          IF ( RC /= QFYAML_Success ) THEN
             errMsg = 'Error encountered in "Copy_Anchor_Variable"!'
             CALL Handle_Error( errMsg, RC, thisLoc )
             RETURN
          ENDIF
       ENDDO
    ENDDO

  END SUBROUTINE QFYAML_Read_File
!EOC
!------------------------------------------------------------------------------
! QFYAML: Bob Yantosca | yantosca@seas.harvard.edu | Apr 2020
! Based on existing package https://github.com/jannisteunissen/config_fortran
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Parse_Line
!
! !DESCRIPTION: Parses a single line of a YAML file and adds the relevant
!  variables to the yml configuration object.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Parse_Line( yml,        yml_anchored, set_by,                   &
                         line_arg,   valid_syntax, category,                 &
                         anchor_ptr, anchor_tgt,   RC                       )
!
! !INPUT PARAMETERS:
!
    INTEGER,                      INTENT(IN)    :: set_by
    CHARACTER(LEN=*),             INTENT(IN)    :: line_arg
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(QFYAML_t),               INTENT(INOUT) :: yml
    TYPE(QFYAML_t),               INTENT(INOUT) :: yml_anchored
!
! !OUTPUT PARAMETERS:
!
    LOGICAL,                      INTENT(OUT)   :: valid_syntax
    CHARACTER(LEN=QFYAML_namlen), INTENT(OUT)   :: category
    CHARACTER(LEN=QFYAML_namlen), INTENT(OUT)   :: anchor_ptr
    CHARACTER(LEN=QFYAML_namlen), INTENT(OUT)   :: anchor_tgt
    INTEGER,                      INTENT(OUT)   :: RC
!
! !REVISION HISTORY:
!  15 Apr 2020 - R. Yantosca - Initial version
!  See the subsequent Git history with the gitk browser!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    LOGICAL                      :: append
    INTEGER                      :: ix
    INTEGER                      :: ampsnd_ix
    INTEGER                      :: anchor_ix
    INTEGER                      :: colon_ix
    INTEGER                      :: star_ix
    INTEGER                      :: trim_len

    ! Strings
    CHARACTER(LEN=QFYAML_namlen) :: var_name
    CHARACTER(LEN=QFYAML_strlen) :: errMsg
    CHARACTER(LEN=QFYAML_strlen) :: line
    CHARACTER(LEN=QFYAML_strlen) :: thisLoc

    !=======================================================================
    ! Parse_Line begins here!
    !=======================================================================

    ! Initialize
    RC           = QFYAML_Success
    colon_ix     = 0
    ampsnd_ix    = 0
    star_ix      = 0
    line         = line_arg
    valid_syntax = .true.
    errMsg       = ''
    thisLoc      = '-> at Parse_Line (in module qfyaml_mod.F90)'

    ! Strip all comments from the line
    CALL Trim_Comment(line, '#;')

    ! Skip empty lines
    IF ( line == "" ) RETURN

    ! Get the length of the line (excluding trailing whitespace)
    trim_len  = LEN_TRIM( line )

    ! Look for the positions of certain characters
    colon_ix  = INDEX( line, ":" )
    ampsnd_ix = INDEX( line, "&" )
    star_ix   = INDEX( line, "*" )

    ! If the text is flush with the first column
    ! and has a colon in the line, then it's a category
    ! FUTURE TODO: Check the indentation level to add nested categories
    IF ( line(1:1) /= "" .and. colon_ix > 0 ) THEN

       ! If there is nothing after the colon ...
       IF ( colon_ix == trim_len ) THEN

          ! Then this indicates a category name, so return
          ! and there is no anchor present
          anchor_tgt = ""
          category   = line(1:colon_ix-1)
          IF ( category == "'NO'" ) category = "NO"   ! Avoid clash w/ FALSE
          RETURN

       ! If there is an ampersand following the colon
       ! then this denotes a YAML anchor
       ELSE IF ( colon_ix > 0 .and. ampsnd_ix > 0 ) THEN

          ! Return anchor target and category name
          anchor_tgt = line(ampsnd_ix+1:trim_len)
          category   = line(1:colon_ix-1)
          IF ( category == "'NO'" ) category = "NO"   ! Avoid clash w/ FALSE
          RETURN
       ENDIF
    ENDIF

    ! ... otherwise we'll use a category of "" for the variable name

    ! define the variable name
    append = .FALSE.
    var_name = line(1:colon_ix-1)

    ! If there are less than two spaces or a tab, reset to no category
    IF (var_name(1:2) /= " " .and. var_name(1:1) /= tab_char) THEN
       category = ""
    ENDIF

    ! Replace leading tabs by spaces
    ix = VERIFY(var_name, tab_char) ! Find first non-tab character
    var_name(1:ix-1) = ""

    ! Remove leading blanks
    var_name = ADJUSTL(var_name)

    ! Test if the variable is a YAML anchor
    IF ( var_name == "<<" ) THEN

       !--------------------------------------------------------------------
       ! Variable points to a YAML anchor
       !--------------------------------------------------------------------

       ! Add category if it is defined
       IF ( category /= "" ) THEN
          var_name = TRIM( category ) // QFYAML_category_separator // var_name
       ENDIF

       ! Get the name of the anchor we want to point to
       anchor_ptr = line(star_ix+1:)

       ! Add the variable to the extra configuration object
       ! which we will pass back to routine QFYAML_Read_File.
       ! There we will create a new variable with all of the
       ! properties of the anchor target.
       CALL Add_Variable( yml            = yml_anchored,                     &
                          append         = append,                           &
                          set_by         = set_by,                           &
                          line_arg       = "",                               &
                          anchor_ptr_arg = anchor_ptr,                       &
                          anchor_tgt_arg = anchor_tgt,                       &
                          category_arg   = category,                         &
                          var_name_arg   = var_name,                         &
                          RC             = RC                               )

       ! Trap potential errors
       IF ( RC /= QFYAML_Success ) THEN
          errMsg = 'Error encountered in "Add_Variable" (points to anchor)!'
          CALL Handle_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

    ELSE

       !--------------------------------------------------------------------
       ! Variable does NOT point to a YAML anchor
       !--------------------------------------------------------------------

       ! Add category if it is defined
       IF ( category /= "" ) THEN
          var_name = TRIM( category ) // QFYAML_category_separator // var_name
       ENDIF

       ! Set line to the values behind the '=' sign
       line = line(colon_ix+1:)

       ! Add the variable to the config object
       CALL Add_Variable( yml            = yml,                              &
                          append         = append,                           &
                          set_by         = set_by,                           &
                          line_arg       = line,                             &
                          anchor_ptr_arg = anchor_ptr,                       &
                          anchor_tgt_arg = anchor_tgt,                       &
                          category_arg   = category,                         &
                          var_name_arg   = var_name,                         &
                          RC             = RC                               )

       ! Trap potential errors
       IF ( RC /= QFYAML_Success ) THEN
          errMsg = 'Error encountered in "Add_Variable"!'
          CALL Handle_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

    ENDIF

  END SUBROUTINE Parse_Line
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Add_Variable
!
! !DESCRIPTION: Adds a new variable to the config object.  Either it creates
!  the variable as "stored" (aka deferred) or actual.  This was split off
!  from the routine Parse\_Line above.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Add_Variable( yml,            line_arg,     anchor_ptr_arg,     &
                           anchor_tgt_arg, category_arg, var_name_arg,       &
                           set_by,         append,       RC                 )
!
! !INPUT PARAMETERS:
!
    LOGICAL,                      INTENT(IN)    :: append
    INTEGER,                      INTENT(IN)    :: set_by
    CHARACTER(LEN=*),             INTENT(IN)    :: line_arg
    CHARACTER(LEN=QFYAML_namlen), INTENT(IN)    :: anchor_ptr_arg
    CHARACTER(LEN=QFYAML_namlen), INTENT(IN)    :: anchor_tgt_arg
    CHARACTER(LEN=QFYAML_namlen), INTENT(IN)    :: category_arg
    CHARACTER(LEN=QFYAML_namlen), INTENT(IN)    :: var_name_arg
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(QFYAML_t),               INTENT(INOUT) :: yml
!
! !OUTPUT PARAMETERS:
!
    INTEGER,                      INTENT(OUT)   :: RC
!
! !REVISION HISTORY:
!  15 Apr 2020 - R. Yantosca - Initial version
!  See the subsequent Git history with the gitk browser!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER            :: ix

    ! Strings
    CHARACTER(LEN=255) :: errMsg, thisLoc

    !=======================================================================
    ! Add_Variable_Var begins here!
    !=======================================================================

    ! Initialize
    RC      = QFYAML_success
    errMsg  = ''
    thisLoc = ' -> at Add_Variable (in module qfyaml_mod.F90)'

    ! Find variable corresponding to name in file
    CALL Get_Var_Index( yml, var_name_arg, ix )

    IF ( ix <= 0 ) then

       !--------------------------------------------------------------------
       ! Variable is not already present in the yml object
       !--------------------------------------------------------------------

       ! Prepare to store the data as a string
       CALL Prepare_Store_Var( yml          = yml,                           &
                               var_name     = TRIM( var_name_arg ),          &
                               var_type     = QFYAML_unknown_type,           &
                               var_size     = 1,                             &
                               description  = "Not yet created",             &
                               ix           = ix,                            &
                               dynamic_size = .FALSE.,                       &
                               RC           = RC                            )

       IF ( RC /= QFYAML_Success ) THEN
          errMsg = 'Error encountered in "Prepare_Store_Var"!'
          CALL Handle_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       ! Store the value of the mapping in the "stored_data" field
       yml%vars(ix)%stored_data = TRIM( line_arg )

    ELSE

       !--------------------------------------------------------------------
       ! Variable is already present in the yml object
       !--------------------------------------------------------------------
       IF ( append ) THEN

          ! Append data to data that is already present
          yml%vars(ix)%stored_data = TRIM( yml%vars(ix)%stored_data )     // &
                                     TRIM( line_arg                 )

       ELSE

          ! Or store overwrite existing data
          yml%vars(ix)%stored_data = line_arg

       ENDIF

       ! If type is known, read in values
       IF ( yml%vars(ix)%var_type /= QFYAML_unknown_type ) THEN
          CALL Read_Variable( yml%vars(ix), RC )
          IF ( RC /= QFYAML_Success ) THEN
             errMsg = 'Error encountered in "Read_Variable"!'
             CALL Handle_Error( errMsg, RC, thisLoc )
             RETURN
          ENDIF
       ENDIF

    ENDIF

    ! Store other fields of this variable
    yml%vars(ix)%anchor_tgt = anchor_tgt_arg
    yml%vars(ix)%anchor_ptr = anchor_ptr_arg
    yml%vars(ix)%category   = category_arg
    yml%vars(ix)%set_by     = set_by

  END SUBROUTINE Add_Variable
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Copy_Anchor_Variable
!
! !DESCRIPTION: Adds a new variable that is a copy of a variable with a
!  YAML anchor.  The new variable will contain all of the field values of
!  the old variable.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Copy_Anchor_Variable( yml,          anchor_ix,                  &
                                   var_w_anchor, var_pt_to_anchor,           &
                                   RC                                       )
!
! !INPUT PARAMETERS:
!
    INTEGER,                      INTENT(IN)    :: anchor_ix
    CHARACTER(LEN=QFYAML_namlen), INTENT(IN)    :: var_w_anchor
    CHARACTER(LEN=QFYAML_namlen), INTENT(IN)    :: var_pt_to_anchor
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(QFYAML_t),               INTENT(INOUT) :: yml
!
! !OUTPUT PARAMETERS:
!
    INTEGER,                      INTENT(OUT)   :: RC
!
! !REVISION HISTORY:
!  15 Apr 2020 - R. Yantosca - Initial version
!  See the subsequent Git history with the gitk browser!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER :: ix

    ! Strings
    CHARACTER(LEN=QFYAML_strlen) :: errMsg
    CHARACTER(LEN=QFYAML_strlen) :: thisLoc

    !=======================================================================
    ! Copy_Anchor_Variable begins here!
    !=======================================================================

    ! Initialize
    RC           = QFYAML_Success
    errMsg       = ''
    thisLoc      = ' -> at Copy_Anchor_Variable (in module qfyaml_mod.F90)'

    ! Prepare to store the data as a string
    CALL Prepare_Store_Var( yml          = yml,                              &
                            var_name     = TRIM( var_pt_to_anchor ),         &
                            var_type     = QFYAML_unknown_type,              &
                            var_size     = 1,                                &
                            description  = "Not yet created",                &
                            ix           = ix,                               &
                            dynamic_size = .FALSE.,                          &
                            RC           = RC                               )

    ! Trap potential errors
    IF ( RC /= QFYAML_Success ) THEN
       errMsg = 'Error encountered at "Prepare_Store_Var"!'
       CALL Handle_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    ! Copy each field of the "anchor" variable to the new variable
    yml%vars(ix)%category    = yml%vars(anchor_ix)%category
    yml%vars(ix)%description = yml%vars(anchor_ix)%description
    yml%vars(ix)%set_by      = yml%vars(anchor_ix)%set_by
    yml%vars(ix)%stored_data = yml%vars(anchor_ix)%stored_data
    yml%vars(ix)%anchor_ptr  = yml%vars(anchor_ix)%anchor_tgt
    yml%vars(ix)%anchor_tgt  = ""

  END SUBROUTINE Copy_Anchor_Variable
!EOC
!------------------------------------------------------------------------------
! QFYAML: Bob Yantosca | yantosca@seas.harvard.edu | Apr 2020
! Based on existing package https://github.com/jannisteunissen/config_fortran
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Read_Variable
!
! !DESCRIPTION: Get the start and end positions of the line content,
!  and the number of entries.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Read_Variable( var, RC )
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(QFYAML_var_t), INTENT(INOUT) :: var
!
! !OUTPUT PARAMETERS:
!
    INTEGER,            INTENT(OUT)   :: RC
!
! !REVISION HISTORY:
!  15 Apr 2020 - R. Yantosca - Initial version
!  See the subsequent Git history wi th the gitk browser!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER                      :: n
    INTEGER                      :: n_entries
    INTEGER                      :: stat

    ! Arrays
    INTEGER                      :: ix_start(QFYAML_MaxArr)
    INTEGER                      :: ix_end(QFYAML_MaxArr)

    ! Strings
    CHARACTER(LEN=QFYAML_strlen) :: errMsg
    CHARACTER(LEN=QFYAML_strlen) :: thisLoc
    CHARACTER(LEN=QFYAML_strlen) :: s1
    CHARACTER(LEN=QFYAML_strlen) :: s2
    CHARACTER(LEN=QFYAML_strlen) :: s3
    CHARACTER(LEN=QFYAML_strlen) :: s4

    !=======================================================================
    ! Read_Variable begins here!
    !=======================================================================

    ! Initialize
    RC        = QFYAML_Success
    ix_start  = 0
    ix_end    = 0
    n         = 0
    n_entries = 0
    errMsg    = ''
    thisLoc   = ' -> at Read_Variable (in module qfyaml_mod.F90)'

    ! Read the stored data field
    CALL Get_Fields_String( var%stored_data, QFYAML_separators,              &
                            QFYAML_brackets, QFYAML_MaxArr,                  &
                            n_entries,       ix_start,                       &
                            ix_end                                          )

    IF ( var%var_size /= n_entries ) THEN

       IF ( .not. var%dynamic_size ) THEN

          ! Allow strings of length 1 to be automatically concatenated
          IF ( var%var_type == QFYAML_string_type .and.                      &
               var%var_size == 1                        ) THEN
             var%char_data(1) =                                              &
                  TRIM(var%stored_data(ix_start(1):ix_end(n_entries)))
             RETURN

          ELSE

             ! Return with error
             errMsg = "Variable " // TRIM( var%var_name )                 // &
                      " has the wrong size"
             CALL Handle_Error( errMsg, RC, thisLoc )
             RETURN
          ENDIF
       ELSE
          var%var_size = n_entries
          CALL Resize_Storage(var)
       ENDIF
    ENDIF

    DO n = 1, n_entries
       stat = 0
       SELECT CASE ( var%var_type )

          CASE( QFYAML_INTEGER_type )
             READ( var%stored_data(ix_start(n):ix_end(n)), *, iostat=stat )  &
                  var%int_data(n)

          CASE( QFYAML_real_type )
             READ( var%stored_data(ix_start(n):ix_end(n)), *, iostat=stat )  &
                  var%real_data(n)

          CASE( QFYAML_string_type )
             var%char_data(n) = TRIM( var%stored_data(ix_start(n):ix_end(n)) )

          CASE( QFYAML_bool_type )
             READ( var%stored_data(ix_start(n):ix_end(n)), *, iostat=stat )  &
                  var%bool_data(n)
       END SELECT

       ! Exit with error if the variable can't be read properly
       IF ( stat /= 0 ) THEN
          s1 = "Error parsing YAML file!"
          s2 = "Reading variable : " // var%var_name
          s3 = "Variable type    : " // QFYAML_type_names(var%var_type)
          s4 = "Parsing value    : " // var%stored_data(ix_start(n):ix_end(n))
          errMsg = TRIM( s1 ) // NEW_LINE( 'a' ) //                          &
                   TRIM( s2 ) // NEW_LINE( 'a' ) //                          &
                   TRIM( s3 ) // NEW_LINE( 'a' ) //                          &
                   TRIM( s4 )
          CALL Handle_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF
    ENDDO

  END SUBROUTINE Read_Variable
!EOC
!------------------------------------------------------------------------------
! QFYAML: Bob Yantosca | yantosca@seas.harvard.edu | Apr 2020
! Based on existing package https://github.com/jannisteunissen/config_fortran
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Trim_Comment
!
! !DESCRIPTION:  Strip comments, but only outside quoted strings
!  (so that var = '#yolo' is valid when # is a comment char)
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Trim_Comment(line, comment_chars)
!
! !INPUT PARAMETERS:
!
    CHARACTER(LEN=*), INTENT(IN)    :: comment_chars

!
! !INPUT/OUTPUT PARAMETERS:
!
    CHARACTER(LEN=*), INTENT(INOUT) :: line
!
! !REVISION HISTORY:
!  15 Apr 2020 - R. Yantosca - Initial version
!  See the subsequent Git history with the gitk browser!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER          :: n

    ! Strings
    CHARACTER(LEN=1) :: current_char
    CHARACTER(LEN=1) :: need_char

    !=======================================================================
    ! Trim_Comment begins here!
    !=======================================================================
    need_char = ""

    DO n = 1, LEN(line)
       current_char = line(n:n)

       IF (need_char == "") THEN
          IF (current_char == "'") THEN
             need_char = "'"    ! Open string
          ELSE IF (current_char == '"') THEN
             need_char = '"'    ! Open string
          ELSE IF (INDEX(comment_chars, current_char) /= 0) THEN
             line = line(1:n-1) ! Trim line up to comment character
             EXIT
          ENDIF
       ELSE IF (current_char == need_char) THEN
          need_char = ""        ! Close string
       ENDIF

    ENDDO

  END SUBROUTINE Trim_Comment
!EOC
!------------------------------------------------------------------------------
! QFYAML: Bob Yantosca | yantosca@seas.harvard.edu | Apr 2020
! Based on existing package https://github.com/jannisteunissen/config_fortran
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: QFYAML_Check
!
! !DESCRIPTION: Checks a QFYAML configuration variable.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE QFYAML_Check( yml, RC )
!
! !INPUT PARAMETERS:
!
    TYPE(QFYAML_t), INTENT(IN)  :: yml
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT) :: RC
!
! !REVISION HISTORY:
!  15 Apr 2020 - R. Yantosca - Initial version
!  See the subsequent Git history with the gitk browser!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER                      :: n

    ! Strings
    CHARACTER(LEN=QFYAML_strlen) :: errMsg
    CHARACTER(LEN=QFYAML_strlen) :: thisLoc

    !=======================================================================
    ! QFYAML_Check begins here!
    !=======================================================================

    ! Initialize
    RC      = QFYAML_success
    errMsg  = ''
    thisLoc = ' -> at QFYAML_Check (in module qfyaml_mod.F90)'

    DO n = 1, yml%num_vars

       IF ( yml%vars(n)%var_type == QFYAML_unknown_type) THEN
          errMsg = "Unknown variable " // TRIM( yml%vars(n)%var_name )    // &
                   " specified"
          CALL Handle_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

    ENDDO

  END SUBROUTINE QFYAML_check
!EOC
!------------------------------------------------------------------------------
! QFYAML: Bob Yantosca | yantosca@seas.harvard.edu | Apr 2020
! Based on existing package https://github.com/jannisteunissen/config_fortran
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Split_Category
!
! !DESCRIPTION: splits the category and the var name
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Split_Category( variable, category, var_name )
!
! !INPUT PARAMETERS:
!
    TYPE(QFYAML_var_t),       INTENT(IN)  :: variable

!
! !INPUT/OUTPUT PARAMETERS:
!
    CHARACTER(QFYAML_namlen), INTENT(OUT) :: category
    CHARACTER(QFYAML_namlen), INTENT(OUT) :: var_name
!
! !REMARKS:
!  TO DO: Support nested categories
!
! !REVISION HISTORY:
!  15 Apr 2020 - R. Yantosca - Initial version
!  See the subsequent Git history with the gitk browser!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER :: ix

    !=======================================================================
    ! Split_Category begins here!
    !=======================================================================
    ix = INDEX( variable%var_name, QFYAML_category_separator )

    IF ( ix == 0 ) THEN
       category = ""
       var_name = variable%var_name
    ELSE
       category = variable%var_name(1:ix-1)
       var_name = variable%var_name(ix+1:)
    ENDIF

  END SUBROUTINE split_category
!EOC
!------------------------------------------------------------------------------
! QFYAML: Bob Yantosca | yantosca@seas.harvard.edu | Apr 2020
! Based on existing package https://github.com/jannisteunissen/config_fortran
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Resize_Storage
!
! !DESCRIPTION: Resize the storage size of variable, which can be of type
!  INTEGER, LOGICAL, REAL, or CHARACTER
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Resize_Storage( variable )
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(QFYAML_var_t), INTENT(INOUT) :: variable
!
! !REVISION HISTORY:
!  15 Apr 2020 - R. Yantosca - Initial version
!  See the subsequent Git history with the gitk browser!
!EOP
!------------------------------------------------------------------------------
!BOC
    !=======================================================================
    ! Resize_Storage begins here!
    !=======================================================================
    SELECT CASE ( variable%var_type )

       CASE( QFYAML_INTEGER_type )
          DEALLOCATE( variable%int_data )
          ALLOCATE( variable%int_data(variable%var_size) )

       CASE( QFYAML_bool_type )
          DEALLOCATE( variable%bool_data )
          ALLOCATE( variable%bool_data(variable%var_size) )

       CASE( QFYAML_real_type )
          DEALLOCATE( variable%real_data )
          ALLOCATE( variable%real_data(variable%var_size) )

       CASE( QFYAML_string_type )
          DEALLOCATE( variable%char_data )
          ALLOCATE( variable%char_data(variable%var_size) )

    END SELECT

  END SUBROUTINE resize_storage
!EOC
!------------------------------------------------------------------------------
! QFYAML: Bob Yantosca | yantosca@seas.harvard.edu | APR 2020
! Based on existing package https://github.com/jannisteunissen/config_fortran
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Prepare_Store_Var
!
! !DESCRIPTION: Helper routine to store variables. This is useful because
!  a lot of the same code is executed for the different types of variables.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Prepare_Store_Var( yml,      var_name,    var_type,             &
                                var_size, description, ix,                   &
                                RC,       dynamic_size                      )
!
! !INPUT PARAMETERS:
!
    CHARACTER(LEN=*), INTENT(IN)   :: var_name
    INTEGER,          INTENT(IN)   :: var_type
    INTEGER,          INTENT(IN)   :: var_size
    CHARACTER(LEN=*), INTENT(IN)   :: description
    LOGICAL,          OPTIONAL     :: dynamic_size
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(QFYAML_t),  INTENT(INOUT) :: yml

!
! !OUTPUT PARAMETERS:
!
    INTEGER,         INTENT(OUT)   :: ix          ! Index of variable
    INTEGER,         INTENT(OUT)   :: RC
!
! !REVISION HISTORY:
!  15 Apr 2020 - R. Yantosca - Initial version
!  See the subsequent Git history with the gitk browser!
!EOP
!------------------------------------------------------------------------------
!BOC
    ! Strings
    CHARACTER(LEN=QFYAML_strlen) :: errMsg
    CHARACTER(LEN=QFYAML_strlen) :: thisLoc

    !=======================================================================
    ! Prepare_Store_Var begins here!
    !=======================================================================

    ! Initialize
    RC      = QFYAML_Success
    errMsg  = ''
    thisLoc = ' -> at Prepare_Store_Var (in module qfyaml_mod.F90)'

    ! Check if variable already exists
    CALL get_var_index(yml, var_name, ix)

    IF (ix == -1) THEN ! Create a new variable
       CALL ensure_free_storage(yml)
       yml%sorted               = .false.
       ix                       = yml%num_vars + 1
       yml%num_vars             = yml%num_vars + 1
       yml%vars(ix)%used        = .false.
       yml%vars(ix)%stored_data = unstored_data_string
    ELSE
       ! Only allowed when the variable is not yet created
       IF (yml%vars(ix)%var_type /= QFYAML_unknown_type) then
          errMsg = "Variable " // trim(var_name) // " already exists"
          CALL Handle_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF
    ENDIF

    yml%vars(ix)%var_name    = var_name
    yml%vars(ix)%description = description
    yml%vars(ix)%var_type    = var_type
    yml%vars(ix)%var_size    = var_size

    IF ( PRESENT( dynamic_size ) ) THEN
       yml%vars(ix)%dynamic_size = dynamic_size
    ELSE
       yml%vars(ix)%dynamic_size = .false.
    ENDIF

    SELECT CASE ( var_type )

       CASE( QFYAML_INTEGER_type )
          ALLOCATE( yml%vars(ix)%int_data(var_size) )

       CASE( QFYAML_real_type )
          ALLOCATE( yml%vars(ix)%real_data(var_size) )

       CASE( QFYAML_string_type )
          ALLOCATE( yml%vars(ix)%char_data(var_size) )

       CASE( QFYAML_bool_type )
          ALLOCATE( yml%vars(ix)%bool_data(var_size) )

    END SELECT

  END SUBROUTINE Prepare_Store_Var
!EOC
!------------------------------------------------------------------------------
! QFYAML: Bob Yantosca | yantosca@seas.harvard.edu | Apr 2020
! Based on existing package https://github.com/jannisteunissen/config_fortran
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Prepare_Get_Var
!
! !DESCRIPTION: Helper routine to get variables. This is useful because a
!  lot of the same code is executed for the different types of variables.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Prepare_Get_Var( yml, var_name, var_type, var_size, ix, RC )
!
! !INPUT PARAMETERS:
!
    CHARACTER(LEN=*), INTENT(IN)    :: var_name
    INTEGER,          INTENT(IN)    :: var_type
    INTEGER,          INTENT(IN)    :: var_size
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(QFYAML_t),   INTENT(INOUT) :: yml
!
! !OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(OUT)   :: ix
    INTEGER,          INTENT(OUT)   :: RC
!
! !REVISION HISTORY:
!  15 Apr 2020 - R. Yantosca - Initial version
!  See the subsequent Git history with the gitk browser!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    CHARACTER(LEN=QFYAML_strlen) :: errMsg
    CHARACTER(LEN=QFYAML_strlen) :: thisLoc

    !=======================================================================
    ! Prepare_Get_Var begins here!
    !=======================================================================

    ! Initialize
    RC      = QFYAML_success
    errMsg  = ""
    thisLoc = " -> at Prepare_Get_Var (in module qfyaml_mod.F90)"

    CALL Get_Var_Index( yml, var_name, ix )

    IF ( ix == QFYAML_Failure ) THEN

       ! Couldn't find variable, exit with error
       errMsg = "QFYAML_get: variable " // TRIM( var_name) // " not found!"
       CALL Handle_Error( errMsg, RC, thisLoc )
       RETURN

    ELSE IF ( yml%vars(ix)%var_type /= var_type ) THEN

       ! Variable is different type than expected: exit wit error
       WRITE( errMsg, "(a)" )                                                &
            "Variable " // TRIM( var_name ) // " has different type ("    // &
            TRIM( QFYAML_type_names( yml%vars(ix)%var_type ) )            // &
            ") than requested ("                                          // &
            TRIM( QFYAML_type_names(var_type ) ) // ")"
       CALL Handle_Error( errMsg, RC, thisLoc )
       RETURN

    ELSE IF ( yml%vars(ix)%var_size /= var_size ) THEN

       ! Variable has different size than requested: exit w/ error
       WRITE( errMsg, "(a,i0,a,i0,a)" )                                      &
            "Variable " // TRIM( var_name ) // " has different size (",      &
            yml%vars(ix)%var_size, ") than requested (", var_size, ")"
       CALL Handle_Error( errMsg, RC, thisLoc )

    ELSE

       ! All good, variable will be used
       yml%vars(ix)%used = .true.

    ENDIF

  END SUBROUTINE Prepare_Get_Var
!EOC
!------------------------------------------------------------------------------
! QFYAML: Bob Yantosca | yantosca@seas.harvard.edu | Apr 2020
! Based on existing package https://github.com/jannisteunissen/config_fortran
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ensure_free_storage
!
! !DESCRIPTION: Routine to ensure that enough storage is allocated for the
!  configuration type. If not the new size will be twice as much as the
!  current size. If no storage is allocated yet a minumum amount of storage
!  is allocated.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Ensure_Free_Storage( yml )

!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(QFYAML_t), INTENT(INOUT)   :: yml
!
! !REVISION HISTORY:
!  15 Apr 2020 - R. Yantosca - Initial version
!  See the subsequent Git history with the gitk browser!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    TYPE(QFYAML_var_t), ALLOCATABLE :: yml_copy(:)
    INTEGER                         :: cur_size, new_size
!
! !DEFINED PARAMETERS:
!
    INTEGER,            PARAMETER   :: min_dyn_size = 100

    !=======================================================================
    ! Ensure_Free_storage begins here!
    !=======================================================================
    IF ( ALLOCATED( yml%vars ) ) THEN
       cur_size = SIZE( yml%vars )

       IF ( cur_size < yml%num_vars + 1 ) THEN
          new_size = 2 * cur_size
          ALLOCATE( yml_copy( cur_size ) )
          yml_copy = yml%vars
          DEALLOCATE( yml%vars )
          ALLOCATE( yml%vars( new_size ) )
          yml%vars(1:cur_size) = yml_copy
       ENDIF
    ELSE
       ALLOCATE( yml%vars( min_dyn_size ) )
    ENDIF

  END SUBROUTINE ensure_free_storage
!EOC
!------------------------------------------------------------------------------
! QFYAML: Bob Yantosca | yantosca@seas.harvard.edu | Apr 2020
! Based on existing package https://github.com/jannisteunissen/config_fortran
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get_Fields_String
!
! !DESCRIPTION: Routine to find the indices of entries in a string.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Get_Fields_String( line_arg, delims,  brackets,                 &
                                n_max,    n_found, ixs_start, ixs_end       )
!
! !INPUT PARAMETERS:
!
    CHARACTER(LEN=*), INTENT(IN)   :: line_arg         ! Line to read
    CHARACTER(LEN=*), INTENT(IN)   :: delims           ! Accepted delimiters
    CHARACTER(LEN=*), INTENT(IN)   :: brackets         ! brackets
    INTEGER,          INTENT(IN)   :: n_max            ! Max entries to read
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,         INTENT(INOUT) :: n_found          ! # of entries found
    INTEGER,         INTENT(INOUT) :: ixs_start(n_max) ! start pt. of ith entry
    INTEGER,         INTENT(INOUT) :: ixs_end(n_max)   ! end pt.   of ith entry
!
! !REVISION HISTORY:
!  15 Apr 2020 - R. Yantosca - Initial version
!  See the subsequent Git history with the gitk browser!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER                      :: B, ix, ix_prev

    ! Strings
    CHARACTER(LEN=1            ) :: bkt
    CHARACTER(LEN=QFYAML_strlen) :: line

    !=======================================================================
    ! Get_Fields_String begins here!
    !=======================================================================

    ! Initialize
    ix      = 0
    ix_prev = 0
    n_found = 0
    line    = line_arg

    !=======================================================================
    ! Get_Fields_String begins here!
    !=======================================================================

    ! Strip out brackets from the line
    DO B = 1, LEN( QFYAML_brackets )
       bkt = QFYAML_brackets(B:B)
       ix  = INDEX( line, bkt )
       IF ( ix > 0 ) line(ix:ix) = " "
    ENDDO

    ! Parse the values
    ix = 0
    DO WHILE (n_found < n_max)

       ! Find the starting point of the next entry (a non-delimiter value)
       ix = VERIFY(line(ix_prev+1:), delims)
       IF (ix == 0) EXIT

       n_found            = n_found + 1
       ixs_start(n_found) = ix_prev + ix ! the absolute position in 'line2'

       ! Get the end point of the current entry (next delimiter index minus one)
       ix = SCAN( line(ixs_start(n_found)+1:), delims) - 1

       IF (ix == -1) THEN              ! If there is no last delimiter,
          ixs_end(n_found) = LEN(line) ! the end of the line is the endpoint
       ELSE
          ixs_end(n_found) = ixs_start(n_found) + ix
       ENDIF

       ix_prev = ixs_end(n_found) ! We continue to search from here
    ENDDO

  END SUBROUTINE get_fields_string
!EOC
!------------------------------------------------------------------------------
! QFYAML: Bob Yantosca | yantosca@seas.harvard.edu | Apr 2020
! Based on existing package https://github.com/jannisteunissen/config_fortran
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Binary_Search_Variable
!
! !DESCRIPTION: Performs a binary search for a given variable.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Binary_Search_Variable( yml, var_name, ix )
!
! !INPUT PARAMETERS:
!
    TYPE(QFYAML_t),   INTENT(IN)  :: yml
    CHARACTER(LEN=*), INTENT(IN)  :: var_name
!
! !OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(OUT) :: ix
!
! !REVISION HISTORY:
!  15 Apr 2020 - R. Yantosca - Initial version
!  See the subsequent Git history with the gitk browser!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER :: i_min, i_max, i_mid

    i_min = 1
    i_max = yml%num_vars
    ix    = - 1

    !=======================================================================
    ! Binary_Search_Variable begins here!
    !=======================================================================

    DO WHILE (i_min < i_max)
       i_mid = i_min + (i_max - i_min) / 2
       IF ( LLT( yml%vars(i_mid)%var_name, var_name ) ) THEN
          i_min = i_mid + 1
       ELSE
          i_max = i_mid
       ENDIF
    ENDDO

    ! If not found, binary_search_variable is not set here, and stays -1
    IF ( i_max == i_min .and. yml%vars(i_min)%var_name == var_name ) THEN
       ix = i_min
    ELSE
       ix = -1
    ENDIF

  END SUBROUTINE Binary_Search_Variable
!EOC
!------------------------------------------------------------------------------
! QFYAML: Bob Yantosca | yantosca@seas.harvard.edu | Apr 2020
! Based on existing package https://github.com/jannisteunissen/config_fortran
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: QFYAML_Sort
!
! !DESCRIPTION: Sorts the variables list for faster lookup
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE QFYAML_Sort( yml )
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(QFYAML_t), INTENT(INOUT) :: yml
!
! !REVISION HISTORY:
!  15 Apr 2020 - R. Yantosca - Initial version
 !  See the subsequent Git history with the gitk browser!
!EOP
!------------------------------------------------------------------------------
!BOC
    !=======================================================================
    ! QFYAML_Sort begins here!
    !=======================================================================

    ! Sort the list
    CALL Qsort( yml%vars(1:yml%num_vars) )

    ! Indicate that we have sorted
    yml%sorted = .TRUE.

  END SUBROUTINE QFYAML_sort

!EOC
!------------------------------------------------------------------------------
! QFYAML: Bob Yantosca | yantosca@seas.harvard.edu
! Based on existing package https://github.com/jannisteunissen/config_fortran
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Qsort
!
! !DESCRIPTION: Simple implementation of quicksort algorithm to sort
!  the variable list alphabetically.
!\\
!\\
! !INTERFACE:
!
  RECURSIVE SUBROUTINE Qsort( list )
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(QFYAML_var_t), INTENT(INOUT) :: list(:)
!
! !REVISION HISTORY:
!  15 Apr 2020 - R. Yantosca - Initial version
!  See the subsequent Git history with the gitk browser!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER :: split_pos

    !=======================================================================
    ! Qsort begins here!
    !=======================================================================
    IF ( SIZE(list) > 1 ) then
       CALL Partition_Var_List( list, split_pos )
       CALL Qsort( list(:split_pos-1) )
       CALL Qsort( list(split_pos:) )
    ENDIF

  END SUBROUTINE qsort
!EOC
!------------------------------------------------------------------------------
! QFYAML: Bob Yantosca | yantosca@seas.harvard.edu | Apr 2020
! Based on existing package https://github.com/jannisteunissen/config_fortran
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Partition_Var_List
!
! !DESCRIPTION: Helper routine for quicksort, to perform partitioning
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Partition_Var_List(list, marker)
!
! !INPUT PARAMETERS:
!
    TYPE(QFYAML_var_t), INTENT(INOUT) :: list(:)
!
! !OUTPUT PARAMETERS:
!
    INTEGER,            INTENT(OUT)   :: marker
!
! !REVISION HISTORY:
!  15 Apr2015 - R. Yantosca - Initial version
!  See the subsequent Git history with the gitk browser!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER                      :: left
    INTEGER                      :: right
    INTEGER                      :: pivot_ix

    ! Strings
    CHARACTER(LEN=QFYAML_namlen) :: pivot_value

    ! Objects
    TYPE(QFYAML_var_t)           :: temp

    !=======================================================================
    ! Partition_Var_List begins here!
    !=======================================================================

    left        = 0
    right       = SIZE( list ) + 1

    ! Take the middle element as pivot
    pivot_ix    = SIZE( list ) / 2
    pivot_value = list(pivot_ix)%var_name

    DO WHILE ( left < right )

       right = right - 1
       DO WHILE ( LGT( list(right)%var_name, pivot_value ) )
          right = right - 1
       ENDDO

       left = left + 1
       DO WHILE (LGT( pivot_value, list(left)%var_name ) )
          left = left + 1
       ENDDO

       IF ( left < right ) THEN
          temp = list(left)
          list(left) = list(right)
          list(right) = temp
       ENDIF
    ENDDO

    IF ( left == right ) THEN
       marker = left + 1
    ELSE
       marker = left
    ENDIF

  END SUBROUTINE partition_var_list
!EOC
!------------------------------------------------------------------------------
! QFYAML: Bob Yantosca | yantosca@seas.harvard.edu | Apr 2020
! Based on existing package https://github.com/jannisteunissen/config_fortran
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: QFYAML_CleanUp
!
! !DESCRIPTION: Clear all data from a QFYAML_t object, so that it can be reused.
!  Note that this also happens automatically when such an object goes out
!  of scope.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE QFYAML_CleanUp( yml )
!
! !USES:
!
    IMPLICIT NONE
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(QFYAML_t), INTENT(INOUT) :: yml
!
! !REVISION HISTORY:
!  15 Apr 2020
!  See the subsequent Git history with the gitk browser!
!EOP
!------------------------------------------------------------------------------
!BOC
    !=======================================================================
    ! QFYAML_CleanUp begins here!
    !=======================================================================

    ! Reset scalars
    yml%sorted   = .false.
    yml%num_vars = 0

    ! Deallocate variables array
    IF ( ALLOCATED( yml%vars ) ) DEALLOCATE( yml%vars )

  END SUBROUTINE QFYAML_CleanUp
!EOC
!------------------------------------------------------------------------------
! QFYAML: Bob Yantosca | yantosca@seas.harvard.edu | Apr 2020
! Based on existing package https://github.com/jannisteunissen/config_fortran
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: QFYAML_Get_Size
!
! !DESCRIPTION: Get the size of a variable.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE QFYAML_Get_Size( yml, var_name, res, RC )
!
! !INPUT PARAMETERS:
!
    TYPE(QFYAML_t),   INTENT(IN)  :: yml
    CHARACTER(LEN=*), INTENT(IN)  :: var_name
!
! !OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(OUT) :: res
    INTEGER,          INTENT(OUT) :: RC
!
! !REVISION HISTORY:
!  15 Apr 2020 - R. Yantosca - Initial version
!  See the subsequent Git history with the gitk browser!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER                      :: ix

    ! Strings
    CHARACTER(LEN=QFYAML_strlen) :: errMsg
    CHARACTER(LEN=QFYAML_strlen) :: thisLoc

    !=======================================================================
    ! Prepare_Store_Var begins here!
    !=======================================================================

    ! Initialize
    RC      = QFYAML_success
    errMsg  = ''
    thisLoc = ' -> at QFYAML_Get_Size (in module qfyaml_mod.F90)'

    ! Get variable index
    CALL Get_Var_Index( yml, var_name, ix )

    ! Return var_size or exit w/ error
    IF ( ix /= QFYAML_Failure ) THEN
       res = yml%vars(ix)%var_size
    ELSE
       res = QFYAML_Failure
       errMsg = "Variable [" // TRIM( var_name ) // "] not found"
       CALL Handle_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

  END SUBROUTINE QFYAML_Get_Size
!EOC
!------------------------------------------------------------------------------
! QFYAML: Bob Yantosca | yantosca@seas.harvard.edu | Apr 2020
! Based on existing package https://github.com/jannisteunissen/config_fortran
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: QFYAML_get_type
!
! !DESCRIPTION: Get the type of a given variable of a configuration type
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE QFYAML_Get_Type( yml, var_name, res, RC )
!
! !INPUT PARAMETERS:
!
    TYPE(QFYAML_t),   INTENT(IN)  :: yml
    CHARACTER(LEN=*), INTENT(IN)  :: var_name
!
! !OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(OUT) :: res
    INTEGER,          INTENT(OUT) :: RC
!
! !REVISION HISTORY:
!  15 Apr 2020 - R. Yantosca - Initial version
!  See the subsequent Git history with the gitk browser!
!EOP
!------------------------------------------------------------------------------
!
    ! Scalars
    INTEGER                      :: ix

    ! Strings
    CHARACTER(LEN=QFYAML_strlen) :: errMsg
    CHARACTER(LEN=QFYAML_strlen) :: thisLoc

    !=======================================================================
    ! Prepare_Store_Var begins here!
    !=======================================================================

    ! Initialize
    RC      = QFYAML_success
    errMsg  = ''
    thisLoc = ' -> QFYAML_Get_Type (in module qfyaml_mod.F90)'

    ! Find the variable index
    CALL Get_Var_Index( yml, var_name, ix )

    ! Retrurn var_type or exit w/ error
    IF ( ix /= QFYAML_Failure ) THEN
       res = yml%vars(ix)%var_type
    ELSE
       res = QFYAML_Failure
       errMsg = "Variable [" // TRIM( var_name ) // "] not found"
       CALL Handle_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

  END SUBROUTINE QFYAML_get_type
!EOC
!############################################################################
!### HERE FOLLOWS OVERLOADED MODULE PROCEDURES.  THESE ARE SIMPLE SO
!### WE WILL OMIT ADDING SUBROUTINE HEADERS
!############################################################################

  SUBROUTINE Add_Real( yml, var_name, real_data, comment, RC )
    !
    ! Add a YAML variable with a REAL value
    !
    TYPE(QFYAML_t),   INTENT(INOUT) :: yml
    CHARACTER(LEN=*), INTENT(IN   ) :: var_name
    REAL(yp),         INTENT(IN   ) :: real_data
    CHARACTER(LEN=*), INTENT(IN   ) :: comment
    INTEGER,          INTENT(OUT  ) :: RC

    INTEGER                         :: ix
    CHARACTER(LEN=QFYAML_strlen)    :: errMsg
    CHARACTER(LEN=QFYAML_strlen)    :: thisLoc

    ! Initialize
    RC      = QFYAML_success
    errMsg  = ''
    thisLoc = '-> at Add_Real (in module qfyaml_mod.F90)'

    CALL Prepare_Store_Var( yml, var_name, QFYAML_real_type,                 &
                            1,   comment,  ix,               RC )
    IF ( RC /= QFYAML_Success ) THEN
       errMsg = 'Error encountered in "Prepare_Store_Var"!'
       CALL Handle_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    IF ( yml%vars(ix)%stored_data /= unstored_data_string ) THEN
       CALL Read_Variable( yml%vars(ix), RC )
       IF ( RC /= QFYAML_Success ) THEN
          errMsg = 'Error encountered in "Read_Variable"!'
          CALL Handle_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF
    ELSE
       yml%vars(ix)%real_data(1) = real_data
    ENDIF

  END SUBROUTINE Add_Real

  SUBROUTINE Add_Real_Array( yml,     var_name,  real_data,                  &
                             comment, RC,        dynamic_size               )
    !
    ! Add a YAML variable with an array of type REAL
    !
    TYPE(QFYAML_t),   INTENT(INOUT) :: yml
    CHARACTER(LEN=*), INTENT(IN   ) :: var_name
    REAL(yp),         INTENT(IN   ) :: real_data(:)
    CHARACTER(LEN=*), INTENT(IN   ) :: comment
    INTEGER,          INTENT(OUT  ) :: RC
    LOGICAL,          OPTIONAL      :: dynamic_size

    INTEGER                         :: ix
    CHARACTER(LEN=QFYAML_strlen)    :: errMsg
    CHARACTER(LEN=QFYAML_strlen)    :: thisLoc

    ! Initialize
    RC      = QFYAML_success
    errMsg  = ''
    thisLoc = ' -> at Add_Real_Array (in module qfyaml_mod.F90)'

    CALL Prepare_Store_Var( yml,             var_name,    QFYAML_real_type,  &
                            SIZE(real_data), comment,     ix,                &
                            RC,              dynamic_size                   )
    IF ( RC /= QFYAML_Success ) THEN
       errMsg = 'Error encountered in "Prepare_Store_Var"!'
       CALL Handle_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    IF ( yml%vars(ix)%stored_data /= unstored_data_string ) THEN
       CALL Read_Variable( yml%vars(ix), RC )
       IF ( RC /= QFYAML_Success ) THEN
          errMsg = 'Error encountered in "Read_Variable"!'
          CALL Handle_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF
    ELSE
       yml%vars(ix)%real_data = real_data
    ENDIF

  END SUBROUTINE Add_Real_Array

  SUBROUTINE Add_Int( yml, var_name, int_data, comment, RC )
    !
    ! Add a YAML variable with an INTEGER value
    !
    TYPE(QFYAML_t),   INTENT(INOUT) :: yml
    CHARACTER(LEN=*), INTENT(IN   ) :: var_name
    INTEGER,          INTENT(IN   ) :: int_data
    CHARACTER(LEN=*), INTENT(IN   ) :: comment
    INTEGER,          INTENT(OUT  ) :: RC

    INTEGER                         :: ix
    CHARACTER(LEN=QFYAML_strlen)    :: errMsg
    CHARACTER(LEN=QFYAML_strlen)    :: thisLoc

    ! Initialize
    RC      = QFYAML_success
    errMsg  = ''
    thisLoc = ' -> at Add_Int (in module qfyaml_mod.F90)'

    CALL Prepare_Store_Var( yml, var_name, QFYAML_INTEGER_type,              &
                            1,   comment,  ix,                  RC          )
    IF ( RC /= QFYAML_Success ) THEN
       errMsg = 'Error encountered in "Prepare_Store_Var"!'
       CALL Handle_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    IF ( yml%vars(ix)%stored_data /= unstored_data_string ) THEN
       CALL Read_Variable( yml%vars(ix), RC )
       IF ( RC /= QFYAML_Success ) THEN
          errMsg = 'Error encountered in "Read_Variable"!'
          CALL Handle_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF
    ELSE
       yml%vars(ix)%int_data(1) = int_data
    ENDIF

  END SUBROUTINE Add_Int

  SUBROUTINE Add_Int_Array( yml,     var_name, int_data,                     &
                            comment, RC,       dynamic_size                 )
    !
    ! Add a YAML variable with an array of type INTEGER
    !
    TYPE(QFYAML_t),   INTENT(INOUT) :: yml
    CHARACTER(LEN=*), INTENT(IN   ) :: var_name
    INTEGER,          INTENT(IN   ) :: int_data(:)
    CHARACTER(LEN=*), INTENT(IN   ) :: comment
    LOGICAL,          OPTIONAL      :: dynamic_size
    INTEGER,          INTENT(OUT  ) :: RC

    INTEGER                         :: ix
    CHARACTER(LEN=QFYAML_strlen)    :: errMsg
    CHARACTER(LEN=QFYAML_strlen)    :: thisLoc

    ! Initialize
    RC      = QFYAML_success
    errMsg  = ''
    thisLoc = ' -> at Add_Int_Array (in module qfyaml_mod.F90)'

    CALL Prepare_Store_Var( yml,            var_name, QFYAML_INTEGER_type,   &
                            SIZE(int_data), comment,  ix,                    &
                            RC,             dynamic_size                    )
    IF ( RC /= QFYAML_Success ) THEN
       errMsg = 'Error encountered in "Prepare_Store_Var"!'
       CALL Handle_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    IF ( yml%vars(ix)%stored_data /= unstored_data_string ) THEN
       CALL Read_Variable( yml%vars(ix), RC )
       IF ( RC /= QFYAML_Success ) THEN
          errMsg = 'Error encountered in "Read_Variable"!'
          CALL Handle_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF
    ELSE
       yml%vars(ix)%int_data = int_data
    ENDIF
  END SUBROUTINE Add_Int_Array

  SUBROUTINE Add_String( yml, var_name, char_data, comment, RC )
    !
    ! Add a YAML variable with an CHARACTER value
    !
    TYPE(QFYAML_t),   INTENT(INOUT) :: yml
    CHARACTER(LEN=*), INTENT(IN   ) :: var_name
    CHARACTER(LEN=*), INTENT(IN   ) :: comment
    CHARACTER(LEN=*), INTENT(IN   ) :: char_data
    INTEGER,          INTENT(OUT  ) :: RC

    INTEGER                         :: ix
    CHARACTER(LEN=QFYAML_strlen)    :: errMsg
    CHARACTER(LEN=QFYAML_strlen)    :: thisLoc

    ! Initialize
    RC      = QFYAML_success
    errMsg  = ''
    thisLoc = ' -> at Add_String (in module qfyaml_mod.F90)'

    CALL Prepare_Store_Var( yml, var_name, QFYAML_string_type,               &
                            1,   comment,  ix,                  RC          )

    IF ( RC /= QFYAML_Success ) THEN
       errMsg = 'Error encountered in "Prepare_Store_Var"!'
       CALL Handle_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    IF ( yml%vars(ix)%stored_data /= unstored_data_string ) THEN
       CALL Read_Variable( yml%vars(ix), RC )
       IF ( RC /= QFYAML_Success ) THEN
          errMsg = 'Error encountered in "Read_Variable"!'
          CALL Handle_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF
    ELSE
       yml%vars(ix)%char_data(1) = char_data
    ENDIF

  END SUBROUTINE Add_String

  SUBROUTINE Add_String_Array( yml,     var_name,  char_data,                &
                               comment, RC,        dynamic_size             )
    !
    ! Add a YAML variable with an array of type character
    !
    TYPE(QFYAML_t),   INTENT(INOUT) :: yml
    CHARACTER(LEN=*), INTENT(IN   ) :: var_name
    CHARACTER(LEN=*), INTENT(IN   ) :: comment
    CHARACTER(LEN=*), INTENT(IN   ) :: char_data(:)
    INTEGER,          INTENT(OUT  ) :: RC
    LOGICAL,          OPTIONAL      :: dynamic_size

    INTEGER                         :: ix
    CHARACTER(LEN=QFYAML_strlen)    :: errMsg
    CHARACTER(LEN=QFYAML_strlen)    :: thisLoc

    ! Initialize
    RC      = QFYAML_success
    errMsg  = ''
    thisLoc = ' -> at Add_String_Array (in module qfyaml_mod.F90)'

    CALL Prepare_Store_Var( yml,             var_name,   QFYAML_string_type, &
                            SIZE(char_data), comment,    ix,                 &
                            RC,              dynamic_size                   )

    IF ( yml%vars(ix)%stored_data /= unstored_data_string ) THEN
       CALL Read_Variable( yml%vars(ix), RC )
       IF ( RC /= QFYAML_Success ) THEN
          errMsg = 'Error encountered in "Read_Variable"!'
          CALL Handle_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF
    ELSE
       yml%vars(ix)%char_data = char_data
    ENDIF

  END SUBROUTINE Add_String_Array

  SUBROUTINE Add_Bool( yml, var_name, bool_data, comment, RC )
    !
    ! Add a YAML variable with an logical value
    !
    TYPE(QFYAML_t),   INTENT(INOUT) :: yml
    CHARACTER(LEN=*), INTENT(IN   ) :: var_name
    LOGICAL,          INTENT(IN   ) :: bool_data
    CHARACTER(LEN=*), INTENT(IN   ) :: comment
    INTEGER,          INTENT(OUT  ) :: RC

    INTEGER                         :: ix
    CHARACTER(LEN=QFYAML_strlen)    :: errMsg
    CHARACTER(LEN=QFYAML_strlen)    :: thisLoc

    ! Initialize
    RC      = QFYAML_success
    errMsg  = ''
    thisLoc = ' -> at Add_Bool (in module qfyaml_mod.F90)'

    CALL Prepare_Store_Var( yml, var_name, QFYAML_bool_type,                 &
                            1,   comment,  ix,                RC            )

    IF ( RC /= QFYAML_Success ) THEN
       errMsg = 'Error encountered in "Prepare_Store_Var"!'
       CALL Handle_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    IF ( yml%vars(ix)%stored_data /= unstored_data_string ) THEN
       CALL Read_Variable( yml%vars(ix), RC )
       IF ( RC /= QFYAML_Success ) THEN
          errMsg = 'Error encountered in "Read_Variable"!'
          CALL Handle_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF
    ELSE
       yml%vars(ix)%bool_data(1) = bool_data
    ENDIF

  END SUBROUTINE Add_Bool

  SUBROUTINE Add_Bool_Array(yml,     var_name, bool_data,                    &
                            comment, RC,       dynamic_size                 )
    !
    ! Add a YAML variable with an array of type LOGICAL
    !
    TYPE(QFYAML_t),   INTENT(INOUT) :: yml
    CHARACTER(LEN=*), INTENT(IN   ) :: var_name
    LOGICAL,          INTENT(IN   ) :: bool_data(:)
    CHARACTER(LEN=*), INTENT(IN   ) :: comment
    INTEGER,          INTENT(OUT  ) :: RC
    LOGICAL,          OPTIONAL      :: dynamic_size

    INTEGER                         :: ix
    CHARACTER(LEN=QFYAML_strlen)    :: errMsg
    CHARACTER(LEN=QFYAML_strlen)    :: thisLoc

    ! Initialize
    RC      = QFYAML_success
    errMsg  = ''
    thisLoc = ' -> at Add_Bool_Array (in module qfyaml_mod.F90)'

    CALL Prepare_Store_Var( yml,             var_name,    QFYAML_bool_type,  &
                            SIZE(bool_data), comment,     ix,                &
                            RC,              dynamic_size                   )
    IF ( RC /= QFYAML_Success ) THEN
       errMsg = 'Error encountered in "Prepare_Store_Var"!'
       CALL Handle_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    IF ( yml%vars(ix)%stored_data /= unstored_data_string ) then
       CALL Read_Variable( yml%vars(ix), RC )
       IF ( RC /= QFYAML_Success ) THEN
          errMsg = 'Error encountered in "Read_Variable"!'
          CALL Handle_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF
    ELSE
       yml%vars(ix)%bool_data = bool_data
    ENDIF

  END SUBROUTINE Add_Bool_Array

  SUBROUTINE Get_Real_Array( yml, var_name, real_data, RC )
    !
    ! Get a real array of a given name
    !
    TYPE(QFYAML_t),   INTENT(INOUT) :: yml
    CHARACTER(LEN=*), INTENT(IN   ) :: var_name
    REAL(yp),         INTENT(INOUT) :: real_data(:)
    INTEGER,          INTENT(OUT  ) :: RC

    INTEGER                         :: ix
    CHARACTER(LEN=QFYAML_strlen)    :: errMsg
    CHARACTER(LEN=QFYAML_strlen)    :: thisLoc

    ! Initialize
    RC      = QFYAML_success
    errMsg  = ''
    thisLoc = ' -> at Get_Real_Array (in module qfyaml_mod.F90)'

    CALL Prepare_Get_Var( yml,             var_name, QFYAML_real_type,       &
                          SIZE(real_data), ix,       RC                     )
    IF ( RC /= QFYAML_Success ) THEN
       errMsg = 'Error encountered in "Prepare_Get_Var"!'
       CALL Handle_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    real_data = yml%vars(ix)%real_data

  END SUBROUTINE Get_Real_Array

  SUBROUTINE Get_Int_Array( yml, var_name, int_data, RC )
    !
    ! Get a INTEGER array of a given name
    !

    TYPE(QFYAML_t),   INTENT(INOUT) :: yml
    CHARACTER(LEN=*), INTENT(IN   ) :: var_name
    INTEGER,          INTENT(INOUT) :: int_data(:)
    INTEGER,          INTENT(OUT  ) :: RC

    INTEGER                         :: ix
    CHARACTER(LEN=QFYAML_strlen)    :: errMsg
    CHARACTER(LEN=QFYAML_strlen)    :: thisLoc

    ! Initialize
    RC      = QFYAML_success
    errMsg  = ''
    thisLoc = ' -> at Get_Int_Array (in module qfyaml_mod.F90)'

    CALL Prepare_Get_Var( yml,            var_name, QFYAML_integer_type,     &
                          SIZE(int_data), ix,       RC                      )
    IF ( RC /= QFYAML_Success ) THEN
       errMsg = 'Error encountered in "Prepare_Get_Var"!'
       CALL Handle_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    int_data    = yml%vars(ix)%int_data

  END SUBROUTINE Get_Int_Array

  SUBROUTINE Get_String_Array( yml, var_name, char_data, RC )
    !
    ! Get a character array of a given name
    !
    TYPE(QFYAML_t),   INTENT(INOUT) :: yml
    CHARACTER(LEN=*), INTENT(IN   )    :: var_name
    CHARACTER(LEN=*), INTENT(INOUT) :: char_data(:)
    INTEGER,          INTENT(OUT  ) :: RC

    INTEGER                         :: ix
    CHARACTER(LEN=QFYAML_strlen)    :: errMsg
    CHARACTER(LEN=QFYAML_strlen)    :: thisLoc

    ! Initialize
    RC      = QFYAML_success
    errMsg  = ''
    thisLoc = ' -> at Get_String_Array (in module qfyaml_mod.F90)'

    CALL Prepare_Get_Var( yml,             var_name, QFYAML_string_type,     &
                          SIZE(char_data), ix,       RC                     )
    IF ( RC /= QFYAML_Success ) THEN
       errMsg = 'Error encountered in "Prepare_Get_Var"!'
       CALL Handle_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    char_data = yml%vars(ix)%char_data

  END SUBROUTINE Get_String_Array

  SUBROUTINE Get_Bool_Array( yml, var_name, bool_data, RC )
    !
    ! Get a LOGICAL array of a given name
    !
    TYPE(QFYAML_t),   INTENT(INOUT) :: yml
    CHARACTER(LEN=*), INTENT(IN   ) :: var_name
    LOGICAL,          INTENT(INOUT) :: bool_data(:)
    INTEGER,          INTENT(OUT  ) :: RC

    INTEGER                         :: ix
    CHARACTER(LEN=QFYAML_strlen)    :: errMsg
    CHARACTER(LEN=QFYAML_strlen)    :: thisLoc

    ! Initialize
    RC      = QFYAML_success
    errMsg  = ''
    thisLoc = ' -> at Get_Bool_Array (in module qfyaml_mod.F90)'

    CALL Prepare_Get_Var( yml,             var_name, QFYAML_bool_type,       &
                          SIZE(bool_data), ix,       RC                     )
    IF ( RC /= QFYAML_Success ) THEN
       errMsg = 'Error encountered in "Prepare_Get_Var"!'
       CALL Handle_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    bool_data = yml%vars(ix)%bool_data

  END SUBROUTINE Get_Bool_Array

  SUBROUTINE Get_Real( yml, var_name, res, RC )
    !
    ! Get a real value of a given name
    !
    TYPE(QFYAML_t),   INTENT(INOUT) :: yml
    CHARACTER(LEN=*), INTENT(IN   ) :: var_name
    REAL(yp),         INTENT(OUT  ) :: res
    INTEGER,          INTENT(OUT  ) :: RC

    INTEGER                         :: ix
    CHARACTER(LEN=QFYAML_strlen)    :: errMsg
    CHARACTER(LEN=QFYAML_strlen)    :: thisLoc

    ! Initialize
    RC      = QFYAML_success
    errMsg  = ''
    thisLoc = ' -> at Get_Real (in module qfyaml_mod.F90)'

    CALL Prepare_Get_Var( yml, var_name, QFYAML_real_type, 1, ix, RC )
    IF ( RC /= QFYAML_Success ) THEN
       errMsg = 'Error encountered in "Prepare_Get_Var"!'
       CALL Handle_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    res = yml%vars(ix)%real_data(1)

  END SUBROUTINE Get_Real

  SUBROUTINE Get_Int( yml, var_name, res, RC )
    !
    ! Get a INTEGER value of a given name
    !
    TYPE(QFYAML_t),   INTENT(INOUT) :: yml
    CHARACTER(LEN=*), INTENT(IN   ) :: var_name
    INTEGER,          INTENT(INOUT) :: res
    INTEGER,          INTENT(OUT  ) :: RC

    INTEGER                         :: ix
    CHARACTER(LEN=QFYAML_strlen)    :: errMsg
    CHARACTER(LEN=QFYAML_strlen)    :: thisLoc

    ! Initialize
    RC      = QFYAML_success
    errMsg  = ''
    thisLoc = ' -> at Get_Int (in module qfyaml_mod.F90)'

    CALL Prepare_Get_Var( yml, var_name, QFYAML_integer_type, 1, ix, RC )
    IF ( RC /= QFYAML_Success ) THEN
       errMsg = 'Error encountered in "Prepare_Get_Var"!'
       CALL Handle_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    res = yml%vars(ix)%int_data(1)

  END SUBROUTINE Get_Int

  SUBROUTINE Get_Bool( yml, var_name, res, RC )
    !
    ! Get a LOGICAL value of a given name
    !
    TYPE(QFYAML_t),   INTENT(INOUT) :: yml
    CHARACTER(LEN=*), INTENT(IN   ) :: var_name
    LOGICAL,          INTENT(OUT  ) :: res
    INTEGER,          INTENT(OUT  ) :: RC

    INTEGER                         :: ix
    CHARACTER(LEN=QFYAML_strlen)    :: errMsg
    CHARACTER(LEN=QFYAML_strlen)    :: thisLoc

    ! Initialize
    RC      = QFYAML_success
    errMsg  = ''
    thisLoc = ' -> at Get_Bool (in module qfyaml_mod.F90)'

    CALL Prepare_Get_Var( yml, var_name, QFYAML_bool_type, 1, ix, RC )
    IF ( RC /= QFYAML_Success ) THEN
       errMsg = 'Error encountered in "Prepare_Get_Var"!'
       CALL Handle_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    res = yml%vars(ix)%bool_data(1)

  END SUBROUTINE Get_Bool

  SUBROUTINE Get_String( yml, var_name, res, RC )
    !
    ! Get a character value of a given name
    !
    TYPE(QFYAML_t),   INTENT(INOUT) :: yml
    CHARACTER(LEN=*), INTENT(IN   ) :: var_name
    CHARACTER(LEN=*), INTENT(OUT  ) :: res
    INTEGER,          INTENT(OUT  ) :: RC

    INTEGER                         :: ix
    CHARACTER(LEN=QFYAML_strlen)    :: errMsg
    CHARACTER(LEN=QFYAML_strlen)    :: thisLoc

    ! Initialize
    RC      = QFYAML_success
    errMsg  = ''
    thisLoc = ' -> at Get_String (in module qfyaml_mod.F90)'

    CALL Prepare_Get_Var( yml, var_name, QFYAML_string_type, 1, ix, RC )
    IF ( RC /= QFYAML_Success ) THEN
       errMsg = 'Error encountered in "Prepare_Get_Var"!'
       CALL Handle_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    res = yml%vars(ix)%char_data(1)

  END SUBROUTINE Get_String

  SUBROUTINE Add_Get_Real_Array( yml,     var_name,  real_data,              &
                                 comment, RC,        dynamic_size           )
    !
    ! Get or add a real array of a given name
    !
    TYPE(QFYAML_t),   INTENT(INOUT) :: yml
    CHARACTER(LEN=*), INTENT(IN   ) :: var_name
    REAL(yp),         INTENT(INOUT) :: real_data(:)
    CHARACTER(LEN=*), INTENT(IN   ) :: comment
    INTEGER,          INTENT(OUT  ) :: RC
    LOGICAL,          OPTIONAL      :: dynamic_size

    CHARACTER(LEN=QFYAML_strlen)    :: errMsg
    CHARACTER(LEN=QFYAML_strlen)    :: thisLoc

    ! Initialize
    RC      = QFYAML_success
    errMsg  = ''
    thisLoc = ' -> at Add_Get_Real_Array (in module qfyaml_mod.F90)'

    CALL Add_Real_Array( yml, var_name, real_data, comment, RC, dynamic_size )
    IF ( RC /= QFYAML_Success ) THEN
       errMsg = 'Error encountered in "Add_Real_Array"!'
       CALL Handle_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    CALL Get_Real_Array( yml, var_name, real_data, RC )
    IF ( RC /= QFYAML_Success ) THEN
       errMsg = 'Error encountered in "Add_Real_Array"!'
       CALL Handle_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

  END SUBROUTINE Add_Get_Real_Array

  SUBROUTINE Add_Get_Int_Array( yml,     var_name,  int_data,                &
                                comment, RC,        dynamic_size            )
    !
    ! Get or add a INTEGER array of a given name
    !
    TYPE(QFYAML_t),   INTENT(INOUT) :: yml
    CHARACTER(LEN=*), INTENT(IN   ) :: var_name
    INTEGER,          INTENT(INOUT) :: int_data(:)
    CHARACTER(LEN=*), INTENT(IN   ) :: comment
    LOGICAL,          OPTIONAL      :: dynamic_size
    INTEGER,          INTENT(OUT  ) :: RC

    CHARACTER(LEN=QFYAML_strlen)    :: errMsg
    CHARACTER(LEN=QFYAML_strlen)    :: thisLoc

    ! Initialize
    RC      = QFYAML_success
    errMsg  = ''
    thisLoc = ' -> at Add_Get_Int_Array (in module qfyaml_mod.F90)'

    CALL Add_Int_Array( yml, var_name, int_data, comment, RC, dynamic_size )
    IF ( RC /= QFYAML_Success ) THEN
       errMsg = 'Error encountered in "Add_Int_Array"!'
       CALL Handle_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    CALL Get_Int_Array( yml, var_name, int_data, RC )
    IF ( RC /= QFYAML_Success ) THEN
       errMsg = 'Error encountered in "Get_Int_Array"!'
       CALL Handle_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

  END SUBROUTINE Add_Get_Int_Array

  SUBROUTINE Add_Get_String_Array( yml,     var_name,     char_data,         &
                                   comment, dynamic_size, RC                )
    !
    ! Get or add a character array of a given name
    !
    TYPE(QFYAML_t),   INTENT(INOUT) :: yml
    CHARACTER(LEN=*), INTENT(IN   ) :: var_name
    CHARACTER(LEN=*), INTENT(INOUT) :: char_data(:)
    CHARACTER(LEN=*), INTENT(IN   ) :: comment
    LOGICAL,          optional      :: dynamic_size
    INTEGER,          INTENT(OUT  ) :: RC

    CHARACTER(LEN=QFYAML_strlen)    :: errMsg
    CHARACTER(LEN=QFYAML_strlen)    :: thisLoc

    ! Initialize
    RC      = QFYAML_success
    errMsg  = ''
    thisLoc = ' -> at Add_Get_String_Array (in module qfyaml_mod.F90)'

    CALL Add_String_Array(yml, var_name, char_data, comment, RC, dynamic_size )
    IF ( RC /= QFYAML_Success ) THEN
       errMsg = 'Error encountered in "Add_String_Array"!'
       CALL Handle_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    CALL Get_String_Array( yml, var_name, char_data, RC )
    IF ( RC /= QFYAML_Success ) THEN
       errMsg = 'Error encountered in "Get_String_Array"!'
       CALL Handle_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

  END SUBROUTINE Add_Get_String_Array

  SUBROUTINE Add_Get_Bool_Array(yml,     var_name,     bool_data,            &
                                comment, dynamic_size, RC                   )
    !
    ! Get or add a LOGICAL array of a given name
    !
    TYPE(QFYAML_t),   INTENT(INOUT) :: yml
    CHARACTER(LEN=*), INTENT(IN   ) :: var_name
    LOGICAL,          INTENT(INOUT) :: bool_data(:)
    CHARACTER(LEN=*), INTENT(IN   ) :: comment
    LOGICAL,          OPTIONAL      :: dynamic_size
    INTEGER,          INTENT(OUT  ) :: RC

    CHARACTER(LEN=QFYAML_strlen)    :: errMsg
    CHARACTER(LEN=QFYAML_strlen)    :: thisLoc

    ! Initialize
    RC      = QFYAML_success
    errMsg  = ''
    thisLoc = ' -> at Add_Get_Bool_Array (in module qfyaml_mod.F90)'

    CALL Add_Bool_Array( yml, var_name, bool_data, comment, RC, dynamic_size )
    IF ( RC /= QFYAML_Success ) THEN
       errMsg = 'Error encountered in "Add_Bool_Array"!'
       CALL Handle_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    CALL Get_Bool_Array( yml, var_name, bool_data, RC )
    IF ( RC /= QFYAML_Success ) THEN
       errMsg = 'Error encountered in "Get_Bool_Array"!'
       CALL Handle_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

  END SUBROUTINE Add_Get_Bool_Array

  SUBROUTINE Add_Get_Real( yml, var_name, real_data, comment, RC )
    !
    ! Get or add a real value of a given name
    !
    TYPE(QFYAML_t),   INTENT(INOUT) :: yml
    CHARACTER(LEN=*), INTENT(IN   ) :: var_name
    REAL(yp),         INTENT(INOUT) :: real_data
    CHARACTER(LEN=*), INTENT(IN   ) :: comment
    INTEGER,          INTENT(OUT  ) :: RC

    CHARACTER(LEN=QFYAML_strlen)    :: errMsg
    CHARACTER(LEN=QFYAML_strlen)    :: thisLoc

    ! Initialize
    RC      = QFYAML_success
    errMsg  = ''
    thisLoc = ' -> at Add_Get_Real (in module qfyaml_mod.F90)'

    CALL Add_Real( yml, var_name, real_data, comment, RC )
    IF ( RC /= QFYAML_Success ) THEN
       errMsg = 'Error encountered in "Add_Real"!'
       CALL Handle_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    CALL Get_Real( yml, var_name, real_data, RC )
    IF ( RC /= QFYAML_Success ) THEN
       errMsg = 'Error encountered in "Get_Real"!'
       CALL Handle_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

  END SUBROUTINE Add_Get_Real

  SUBROUTINE Add_Get_Int( yml, var_name, int_data, comment, RC )
    !
    ! Get or add a INTEGER value of a given name
    !
    TYPE(QFYAML_t),   INTENT(INOUT) :: yml
    CHARACTER(LEN=*), INTENT(IN   ) :: var_name
    INTEGER,          INTENT(INOUT) :: int_data
    CHARACTER(LEN=*), INTENT(IN   ) :: comment
    INTEGER,          INTENT(OUT  ) :: RC

    CHARACTER(LEN=QFYAML_strlen)    :: errMsg
    CHARACTER(LEN=QFYAML_strlen)    :: thisLoc

    ! Initialize
    RC      = QFYAML_success
    errMsg  = ''
    thisLoc = ' -> at Add_Get_Int (in module qfyaml_mod.F90)'

    CALL Add_Int( yml, var_name, int_data, comment, RC )
    IF ( RC /= QFYAML_Success ) THEN
       errMsg = 'Error encountered in "Add_Int"!'
       CALL Handle_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    CALL Get_Int( yml, var_name, int_data, RC )
    IF ( RC /= QFYAML_Success ) THEN
       errMsg = 'Error encountered in "Get_Int"!'
       CALL Handle_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

  END SUBROUTINE Add_Get_Int

  SUBROUTINE Add_Get_Bool( yml, var_name, bool_data, comment, RC)
    !
    ! Get or add a LOGICAL value of a given name
    !
    TYPE(QFYAML_t),   INTENT(INOUT) :: yml
    CHARACTER(LEN=*), INTENT(IN   ) :: var_name
    LOGICAL,          INTENT(INOUT) :: bool_data
    CHARACTER(LEN=*), INTENT(IN   ) :: comment
    INTEGER,          INTENT(OUT  ) :: RC

    CHARACTER(LEN=QFYAML_strlen)    :: errMsg
    CHARACTER(LEN=QFYAML_strlen)    :: thisLoc

    ! Initialize
    RC      = QFYAML_success
    errMsg  = ''
    thisLoc = ' -> at Add_Get_Bool (in module qfyaml_mod.F90)'

    CALL Add_Bool( yml, var_name, bool_data, comment, RC )
    IF ( RC /= QFYAML_Success ) THEN
       errMsg = 'Error encountered in "Add_"!'
       CALL Handle_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    CALL Get_Bool( yml, var_name, bool_data, RC )
    IF ( RC /= QFYAML_Success ) THEN
       errMsg = 'Error encountered in "Add_"!'
       CALL Handle_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

  END SUBROUTINE Add_Get_Bool

  SUBROUTINE Add_Get_String( yml, var_name, string_data, comment, RC )

    ! Get a character value of a given name

    TYPE(QFYAML_t),   INTENT(INOUT) :: yml
    CHARACTER(LEN=*), INTENT(IN   ) :: var_name
    CHARACTER(LEN=*), INTENT(INOUT) :: string_data
    CHARACTER(LEN=*), INTENT(IN   ) :: comment
    INTEGER,          INTENT(OUT  ) :: RC

    CHARACTER(LEN=QFYAML_strlen)    :: errMsg
    CHARACTER(LEN=QFYAML_strlen)    :: thisLoc

    ! Initialize
    RC      = QFYAML_success
    errMsg  = ''
    thisLoc = ' -> at Add_Int (in module qfyaml_mod.F90)'

    CALL Add_string( yml, var_name, string_data, comment, RC )
    IF ( RC /= QFYAML_Success ) THEN
       errMsg = 'Error encountered in "Add_String"!'
       CALL Handle_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    CALL Get_String( yml, var_name, string_data, RC )
    IF ( RC /= QFYAML_Success ) THEN
       errMsg = 'Error encountered in "Get_String"!'
       CALL Handle_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

  END SUBROUTINE Add_Get_String

  SUBROUTINE Update_Real_Array( yml, var_name, real_data )
    !
    ! Get or add a real array of a given name
    !
    TYPE(QFYAML_t),   INTENT(INOUT) :: yml
    CHARACTER(LEN=*), INTENT(IN   ) :: var_name
    REAL(yp),         INTENT(INOUT) :: real_data(:)
    INTEGER                         :: ix

    CALL Get_Var_Index( yml, var_name, ix )
    IF ( ix > 0 ) THEN
       yml%vars(ix)%real_data = real_data
       real_data = yml%vars(ix)%real_data
    ENDIF

  END SUBROUTINE Update_Real_Array

  SUBROUTINE Update_Int_Array(yml, var_name, int_data)
    !
    ! Get or add a INTEGER array of a given name
    !
    TYPE(QFYAML_t),   INTENT(INOUT) :: yml
    CHARACTER(LEN=*), INTENT(IN   ) :: var_name
    INTEGER,          INTENT(INOUT) :: int_data(:)
    INTEGER                         :: ix

    CALL Get_Var_Index(yml, var_name, ix)
    IF ( ix > 0 ) THEN
       yml%vars(ix)%int_data = int_data
       int_data = yml%vars(ix)%int_data
    ENDIF

  END SUBROUTINE Update_Int_Array

  SUBROUTINE Update_String_Array( yml, var_name, char_data )
    !
    ! Get or add a character array of a given name
    !
    TYPE(QFYAML_t),   INTENT(INOUT) :: yml
    CHARACTER(LEN=*), INTENT(IN   ) :: var_name
    CHARACTER(LEN=*), INTENT(INOUT) :: char_data(:)
    INTEGER                         :: ix

    CALL Get_Var_Index( yml, var_name, ix )
    IF ( ix > 0 ) THEN
       yml%vars(ix)%char_data = char_data
       char_data = yml%vars(ix)%char_data
    ENDIF

  END SUBROUTINE Update_String_Array

  SUBROUTINE Update_Bool_Array( yml, var_name, bool_data )
    !
    ! Get or add a LOGICAL array of a given name
    !
    TYPE(QFYAML_t),   INTENT(INOUT) :: yml
    CHARACTER(LEN=*), INTENT(IN   ) :: var_name
    LOGICAL,          INTENT(INOUT) :: bool_data(:)
    INTEGER                         :: ix

    CALL Get_Var_Index( yml, var_name, ix )
    IF ( ix > 0 ) THEN
       yml%vars(ix)%bool_data = bool_data
       bool_data = yml%vars(ix)%bool_data
    ENDIF

  END SUBROUTINE Update_Bool_Array

  SUBROUTINE Update_Real( yml, var_name, real_data )
    !
    ! Get or add a real value of a given name
    !
    TYPE(QFYAML_t),   INTENT(INOUT) :: yml
    CHARACTER(LEN=*), INTENT(IN   ) :: var_name
    REAL(yp),         INTENT(INOUT) :: real_data
    INTEGER                         :: ix

    CALL Get_Var_Index(yml, var_name, ix)
    IF ( ix > 0 ) THEN
       yml%vars(ix)%real_data(1) = real_data
       real_data = yml%vars(ix)%real_data(1)
    ENDIF

  END SUBROUTINE Update_Real

  SUBROUTINE Update_Int(yml, var_name, int_data)
    !
    ! Get or add a INTEGER value of a given name
    !
    TYPE(QFYAML_t),   INTENT(INOUT) :: yml
    CHARACTER(LEN=*), INTENT(IN   ) :: var_name
    INTEGER,          INTENT(INOUT) :: int_data
    INTEGER                         :: ix

    CALL Get_Var_Index( yml, var_name, ix )
    IF ( ix > 0 ) THEN
       yml%vars(ix)%int_data(1) = int_data
       int_data = yml%vars(ix)%int_data(1)
    ENDIF

  END SUBROUTINE Update_Int

  SUBROUTINE Update_Bool( yml, var_name, bool_data )
    !
    ! Get or add a LOGICAL value of a given name
    !
    TYPE(QFYAML_t),   INTENT(INOUT) :: yml
    CHARACTER(LEN=*), INTENT(IN   ) :: var_name
    LOGICAL,          INTENT(INOUT) :: bool_data
    INTEGER                         :: ix

    CALL Get_Var_Index( yml, var_name, ix )
    IF ( ix > 0 ) THEN
       yml%vars(ix)%bool_data(1) = bool_data
       bool_data = yml%vars(ix)%bool_data(1)
    ENDIF

  END SUBROUTINE Update_Bool

  SUBROUTINE Update_String(yml, var_name, string_data)
    !
    ! Get a character value of a given name
    !
    TYPE(QFYAML_t),   INTENT(INOUT) :: yml
    CHARACTER(LEN=*), INTENT(IN   ) :: var_name
    CHARACTER(LEN=*), INTENT(INOUT) :: string_data
    INTEGER                         :: ix

    CALL Get_Var_Index( yml, var_name, ix )
    IF ( ix > 0 ) THEN
       yml%vars(ix)%char_data(1) = string_data
       string_data = yml%vars(ix)%char_data(1)
    ENDIF

  END SUBROUTINE Update_String

END MODULE QFYAML_Mod
