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
  ! The precision kind-parameter
  INTEGER, PARAMETER :: yp                    = kind(0.0e0)

  ! Success return code
  INTEGER, PARAMETER :: QFYAML_Success        = 0

  ! Numeric type constants
  INTEGER, PARAMETER :: QFYAML_num_types      = 4
  INTEGER, PARAMETER :: QFYAML_integer_type   = 1
  INTEGER, PARAMETER :: QFYAML_real_type      = 2
  INTEGER, PARAMETER :: QFYAML_string_type    = 3
  INTEGER, PARAMETER :: QFYAML_bool_type      = 4
  INTEGER, PARAMETER :: QFYAML_unknown_type   = 0

  ! How was the YAML file set?
  INTEGER, PARAMETER :: QFYAML_set_by_default = 1
  INTEGER, PARAMETER :: QFYAML_set_by_file    = 3

  ! Maxima
  INTEGER, PARAMETER :: QFYAML_NamLen         = 80    ! Max len for names
  INTEGER, PARAMETER :: QFYAML_StrLen         = 255   ! Max len for strings
  INTEGER, PARAMETER :: QFYAML_MaxArr         = 1000  ! Max entries per array

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
  TYPE, PUBLIC :: QFYAML_t
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
     MODULE PROCEDURE  Update_Real,   Update_Real_Array
     MODULE PROCEDURE  Update_Int,    Update_Int_Array
     MODULE PROCEDURE  Update_Bool,   Update_Bool_Array
     MODULE PROCEDURE  Update_String, Update_String_Array
  END INTERFACE QFYAML_Update

CONTAINS
!EOC
!------------------------------------------------------------------------------
! QFYAML: Bob Yantosca | yantosca@seas.harvard.edu | Apr 2020
! Based on existing package https://github.com/jannisteunissen/config_fortran
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: handle_error
!
! !DESCRIPTION: This routine will be called if an error occurs in one of
!  the subroutines of this module.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Handle_Error( err_string )
!
! !INPUT PARAMETERS:
!
    CHARACTER(LEN=*), INTENT(IN) :: err_string
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
    print *, "The following error occured in m_config:"
    print *, trim(err_string)

    ! It is usually best to quit after an error, to make sure the error message
    ! is not overlooked in the program's output
    error stop
  END SUBROUTINE Handle_Error

!EOC
!------------------------------------------------------------------------------
! QFYAML: Bob Yantosca | yantosca@seas.harvard.edu | Apr 2020
! Based on existing package https://github.com/jannisteunissen/config_fortran
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_var_index
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

    ! Initialize
    ix = -1

    IF ( yml%sorted ) THEN

       ! If the variable names have been sorted, use a binary search
       CALL Binary_Search_Variable( yml, var_name, ix )

    ELSE

       ! Otherwise use a linear search
       DO i = 1, yml%num_vars
          IF ( TRIM(yml%vars(i)%var_name) == TRIM(var_name) ) THEN
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
! !IROUTINE: get_anchor_info
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
    ! Read the YML file
    CALL QFYAML_Read_File( yml, fileName, yml_anchored )

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
! !IROUTINE: QFYAML_Read_File
!
! !DESCRIPTION: Read variables from a YAML file into a configuration object.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE QFYAML_Read_File( yml, fileName, yml_anchored )
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
    CHARACTER(LEN=QFYAML_strlen) :: err_string
    CHARACTER(LEN=QFYAML_strlen) :: line

    !=======================================================================
    ! QFYAML_READ_FILE begins here!
    !=======================================================================

    ! Initialize
    anchor_ptr   = ""
    anchor_tgt   = ""
    category     = ""
    line_number  = 0
    my_unit      = 777

    !=======================================================================
    ! First pass: read the file
    !=======================================================================

    ! Open the file
    OPEN( my_unit, FILE=TRIM(filename), STATUS="old", ACTION="read")
    WRITE( line_fmt, "(A,I0,A)") "(A", QFYAML_strlen, ")"

    ! Start looping
    DO

       ! Read each line and increment the count
       READ( my_unit, FMT=trim(line_fmt), ERR=998, end=999) line
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
                        anchor_ptr   = anchor_ptr                           )

       ! Throw an error if the line is not in proper YAML format
       IF ( .not. valid_syntax ) THEN
          WRITE(err_string, *) "Cannot read line ", line_number, &
               " from ", TRIM(filename)
          CALL Handle_Error(err_string)
       ENDIF
    ENDDO

    ! Error handling
998 WRITE(err_string, "(A,I0,A,I0)") " IOSTAT = ", io_state, &
         " while reading from " // trim(filename) // " at line ", &
         line_number
    CALL Handle_Error("QFYAML_read_file:" // err_string)

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
                                     var_pt_to_anchor = var_pt_to_anchor    )

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
! !IROUTINE: parse line
!
! !DESCRIPTION: Parses a single line of a YAML file and adds the relevant
!  variables to the yml configuration object.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Parse_Line( yml,          yml_anchored, set_by,     line_arg,   &
                         valid_syntax, category,     anchor_ptr, anchor_tgt )
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
    CHARACTER(LEN=QFYAML_strlen) :: line

    !=======================================================================
    ! PARSE_LINE begins here!
    !=======================================================================

    ! Initialize
    valid_syntax = .true.
    colon_ix     = 0
    ampsnd_ix    = 0
    star_ix      = 0
    line         = line_arg

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
          RETURN

       ! If there is an ampersand following the colon
       ! then this denotes a YAML anchor
       ELSE IF ( colon_ix > 0 .and. ampsnd_ix > 0 ) THEN

          ! Return anchor target and category name
          anchor_tgt = line(ampsnd_ix+1:trim_len)
          category   = line(1:colon_ix-1)
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
    IF ( var_name == ">>" ) THEN

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
                          var_name_arg   = var_name                         )

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
                          var_name_arg   = var_name                         )

    ENDIF

  END SUBROUTINE Parse_Line
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: add_variable
!
! !DESCRIPTION: Adds a new variable to the config object.  Either it creates
!  the variable as "stored" (aka deferred) or actual.  This was split off
!  from the routine Parse\_Line above.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Add_Variable( yml,            line_arg,     anchor_ptr_arg,      &
                           anchor_tgt_arg, category_arg, var_name_arg,        &
                           set_by,         append                            )
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
                               dynamic_size = .FALSE.                       )

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
          CALL Read_Variable( yml%vars(ix) )
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
! !IROUTINE: copy_anchor_variable
!
! !DESCRIPTION: Adds a new variable that is a copy of a variable with a
!  YAML anchor.  The new variable will contain all of the field values of
!  the old variable.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Copy_Anchor_Variable( yml,          anchor_ix,                  &
                                   var_w_anchor, var_pt_to_anchor           )
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

    ! Prepare to store the data as a string
    CALL Prepare_Store_Var( yml          = yml,                              &
                            var_name     = TRIM( var_pt_to_anchor ),         &
                            var_type     = QFYAML_unknown_type,              &
                            var_size     = 1,                                &
                            description  = "Not yet created",                &
                            ix           = ix,                               &
                            dynamic_size = .FALSE.                          )

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
! !IROUTINE: read_variable
!
! !DESCRIPTION: Get the start and end positions of the line content,
!  and the number of entries.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE read_variable(var)
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(QFYAML_var_t), INTENT(INOUT) :: var
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
    INTEGER :: n, n_entries
    INTEGER :: ix_start(QFYAML_MaxArr)
    INTEGER :: ix_end(QFYAML_MaxArr), stat

    ! Initialize
    ix_start  = 0
    ix_end    = 0
    n         = 0
    n_entries = 0

    CALL get_fields_string(var%stored_data, QFYAML_separators,  &
                           QFYAML_brackets, QFYAML_MaxArr,      &
                           n_entries,       ix_start,           &
                           ix_end                              )

    if (var%var_size /= n_entries) then

       if (.not. var%dynamic_size) then
          ! Allow strings of length 1 to be automatically concatenated
          if (var%var_type == QFYAML_string_type .and. var%var_size == 1) then
             var%char_data(1) = trim(var%stored_data(ix_start(1):ix_end(n_entries)))
             return ! Leave routine
          else
             CALL handle_error("read_variable: variable " // &
                  & trim(var%var_name) // " has the wrong size")
          end if
       else
          var%var_size = n_entries
          CALL resize_storage(var)
       end if
    end if

    do n = 1, n_entries
       stat = 0
       select case (var%var_type)
       case (QFYAML_INTEGER_type)
          read(var%stored_data(ix_start(n):ix_end(n)), *, iostat=stat) var%int_data(n)
       case (QFYAML_real_type)
          read(var%stored_data(ix_start(n):ix_end(n)), *, iostat=stat) var%real_data(n)
       case (QFYAML_string_type)
          var%char_data(n) = trim(var%stored_data(ix_start(n):ix_end(n)))
       case (QFYAML_bool_type)
          read(var%stored_data(ix_start(n):ix_end(n)), *, iostat=stat) var%bool_data(n)
       end select

       if(stat /= 0) then
          write (*, *) "** m_config error **"
          write (*, *) "reading variable: ", trim(var%var_name)
          write (*, *) "variable type:    ", trim(QFYAML_type_names(var%var_type))
          write (*, *) "parsing value:    ", var%stored_data(ix_start(n):ix_end(n))
          write (*, "(A,I0)") " iostat value:     ", stat
          stop
       endif
    end do
  END SUBROUTINE read_variable
!EOC
!------------------------------------------------------------------------------
! QFYAML: Bob Yantosca | yantosca@seas.harvard.edu | Apr 2020
! Based on existing package https://github.com/jannisteunissen/config_fortran
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: trim_comment
!
! !DESCRIPTION:  Strip comments, but only outside quoted strings
!  (so that var = '#yolo' is valid when # is a comment char)
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE trim_comment(line, comment_chars)
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

  END SUBROUTINE trim_comment
!EOC
!------------------------------------------------------------------------------
! QFYAML: Bob Yantosca | yantosca@seas.harvard.edu | Apr 2020
! Based on existing package https://github.com/jannisteunissen/config_fortran
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE:
!
! !DESCRIPTION:
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE QFYAML_check(yml)
!
! !INPUT PARAMETERS:
!
    TYPE(QFYAML_t), INTENT(IN)       :: yml
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
    CHARACTER(LEN=QFYAML_strlen) :: err_string

    do n = 1, yml%num_vars
       if (yml%vars(n)%var_type == QFYAML_unknown_type) then
          write(err_string, *) "QFYAML_check: unknown variable ", &
               trim(yml%vars(n)%var_name), " specified"
          CALL handle_error(err_string)
       end if
    end do
  END SUBROUTINE QFYAML_check
!EOC
!------------------------------------------------------------------------------
! QFYAML: Bob Yantosca | yantosca@seas.harvard.edu | Apr 2020
! Based on existing package https://github.com/jannisteunissen/config_fortran
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: split_category
!
! !DESCRIPTION: splits the category and the var name
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE split_category(variable, category, var_name)
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

    ix = index(variable%var_name, QFYAML_category_separator)

    if (ix == 0) then
       category = ""
       var_name = variable%var_name
    else
       category = variable%var_name(1:ix-1)
       var_name = variable%var_name(ix+1:)
    end if

  END SUBROUTINE split_category
!EOC
!------------------------------------------------------------------------------
! QFYAML: Bob Yantosca | yantosca@seas.harvard.edu | Apr 2020
! Based on existing package https://github.com/jannisteunissen/config_fortran
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: resize_storage
!
! !DESCRIPTION: Resize the storage size of variable, which can be of type
!  INTEGER, LOGICAL, REAL, or CHARACTER
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE resize_storage(variable)
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
!
! !LOCAL VARIABLES:
!
    SELECT CASE (variable%var_type)
       CASE (QFYAML_INTEGER_type)
          DEALLOCATE( variable%int_data )
          ALLOCATE( variable%int_data(variable%var_size) )
       CASE (QFYAML_bool_type)
          DEALLOCATE( variable%bool_data )
          ALLOCATE( variable%bool_data(variable%var_size) )
       CASE (QFYAML_real_type)
          DEALLOCATE( variable%real_data )
          ALLOCATE( variable%real_data(variable%var_size) )
       CASE (QFYAML_string_type)
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
! !IROUTINE: prepare_store_var
!
! !DESCRIPTION: Helper routine to store variables. This is useful because
!  a lot of the same code is executed for the different types of variables.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE prepare_store_var(yml, var_name, var_type, var_size, &
                               description, ix, dynamic_size)

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
!
! !REVISION HISTORY:
!  15 Apr 2020 - R. Yantosca - Initial version
!  See the subsequent Git history with the gitk browser!
!EOP
!------------------------------------------------------------------------------
!BOC

    ! Check if variable already exists
    CALL get_var_index(yml, var_name, ix)

    if (ix == -1) then ! Create a new variable
       CALL ensure_free_storage(yml)
       yml%sorted               = .false.
       ix                       = yml%num_vars + 1
       yml%num_vars             = yml%num_vars + 1
       yml%vars(ix)%used        = .false.
       yml%vars(ix)%stored_data = unstored_data_string
    else
       ! Only allowed when the variable is not yet created
       if (yml%vars(ix)%var_type /= QFYAML_unknown_type) then
          CALL handle_error("prepare_store_var: variable " // &
               & trim(var_name) // " already exists")
       end if
    end if

    yml%vars(ix)%var_name    = var_name
    yml%vars(ix)%description = description
    yml%vars(ix)%var_type    = var_type
    yml%vars(ix)%var_size    = var_size

    if (present(dynamic_size)) then
       yml%vars(ix)%dynamic_size = dynamic_size
    else
       yml%vars(ix)%dynamic_size = .false.
    end if

    select case (var_type)
    case (QFYAML_INTEGER_type)
       allocate( yml%vars(ix)%int_data(var_size) )
    case (QFYAML_real_type)
       allocate( yml%vars(ix)%real_data(var_size) )
    case (QFYAML_string_type)
       allocate( yml%vars(ix)%char_data(var_size) )
    case (QFYAML_bool_type)
       allocate( yml%vars(ix)%bool_data(var_size) )
    end select

  END SUBROUTINE prepare_store_var
!EOC
!------------------------------------------------------------------------------
! QFYAML: Bob Yantosca | yantosca@seas.harvard.edu | Apr 2020
! Based on existing package https://github.com/jannisteunissen/config_fortran
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: prepare_get_var
!
! !DESCRIPTION: Helper routine to get variables. This is useful because a
!  lot of the same code is executed for the different types of variables.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Prepare_Get_Var(yml, var_name, var_type, var_size, ix)
!
! !INPUT PARAMETERS:
!
    CHARACTER(LEN=*), INTENT(IN)    :: var_name
    INTEGER,          INTENT(IN)    :: var_type, var_size
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(QFYAML_t),   INTENT(INOUT) :: yml
!
! !OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(OUT)   :: ix
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
    CHARACTER(LEN=QFYAML_strlen) :: err_string

    CALL get_var_index(yml, var_name, ix)

    IF (ix == -1) THEN
       CALL handle_error("QFYAML_get: variable "//var_name//" not found")

    ELSE IF (yml%vars(ix)%var_type /= var_type) THEN
       WRITE(err_string, fmt="(A)") "QFYAML_get: variable " &
            // var_name // " has different type (" // &
            TRIM(QFYAML_type_names(yml%vars(ix)%var_type)) // &
            ") than requested (" // TRIM(QFYAML_type_names(var_type)) // ")"
       CALL handle_error(err_string)

    ELSE IF (yml%vars(ix)%var_size /= var_size) THEN
       WRITE(err_string, fmt="(A,I0,A,I0,A)") "QFYAML_get: variable " &
            // var_name // " has different size (", yml%vars(ix)%var_size, &
            ") than requested (", var_size, ")"
       CALL handle_error(err_string)

    ELSE                        ! All good, variable will be used
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
  SUBROUTINE ensure_free_storage(yml)

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
    TYPE(QFYAML_var_t), allocatable :: yml_copy(:)
    INTEGER, parameter           :: min_dyn_size = 100
    INTEGER                      :: cur_size, new_size

    if (allocated(yml%vars)) then
       cur_size = size(yml%vars)

       if (cur_size < yml%num_vars + 1) then
          new_size = 2 * cur_size
          allocate(yml_copy(cur_size))
          yml_copy = yml%vars
          deallocate(yml%vars)
          allocate(yml%vars(new_size))
          yml%vars(1:cur_size) = yml_copy
       end if
    else
       allocate(yml%vars(min_dyn_size))
    end if

  END SUBROUTINE ensure_free_storage
!EOC
!------------------------------------------------------------------------------
! QFYAML: Bob Yantosca | yantosca@seas.harvard.edu | Apr 2020
! Based on existing package https://github.com/jannisteunissen/config_fortran
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE:
!
! !DESCRIPTION: Routine to find the indices of entries in a string
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE get_fields_string(line_arg, delims,  brackets,          &
                               n_max,    n_found, ixs_start, ixs_end)
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

    ! Initialize
    ix      = 0
    ix_prev = 0
    n_found = 0
    line   = line_arg

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
! !IROUTINE: binary_search_variable
!
! !DESCRIPTION: Performs a binary search for the variable 'var_name'
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE binary_search_variable(yml, var_name, ix)
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

    do while (i_min < i_max)
       i_mid = i_min + (i_max - i_min) / 2
       if ( llt(yml%vars(i_mid)%var_name, var_name) ) then
          i_min = i_mid + 1
       else
          i_max = i_mid
       end if
    end do

    ! If not found, binary_search_variable is not set here, and stays -1
    if (i_max == i_min .and. yml%vars(i_min)%var_name == var_name) then
       ix = i_min
    else
       ix = -1
    end if

  END SUBROUTINE binary_search_variable
!EOC
!------------------------------------------------------------------------------
! QFYAML: Bob Yantosca | yantosca@seas.harvard.edu | Apr 2020
! Based on existing package https://github.com/jannisteunissen/config_fortran
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: QFYAML_sort
!
! !DESCRIPTION: Sorts the variables list for faster lookup
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE QFYAML_sort( yml )
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
!
! !LOCAL VARIABLES:
!

    ! Sort the list
    CALL qsort( yml%vars(1:yml%num_vars) )

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
! !IROUTINE: qsort
!
! !DESCRIPTION: Simple implementation of quicksort algorithm to sort
!  the variable list alphabetically.
!\\
!\\
! !INTERFACE:
!
  RECURSIVE SUBROUTINE qsort( list )
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

    IF ( SIZE(list) > 1 ) then
       CALL partition_var_list(list, split_pos)
       CALL qsort( list(:split_pos-1) )
       CALL qsort( list(split_pos:) )
    ENDIF

  END SUBROUTINE qsort
!EOC
!------------------------------------------------------------------------------
! QFYAML: Bob Yantosca | yantosca@seas.harvard.edu | Apr 2020
! Based on existing package https://github.com/jannisteunissen/config_fortran
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: partition_var_list
!
! !DESCRIPTION: Helper routine for quicksort, to perform partitioning
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE partition_var_list(list, marker)
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

    INTEGER                        :: left, right, pivot_ix
    TYPE(QFYAML_var_t)                :: temp
    CHARACTER(LEN=QFYAML_namlen)    :: pivot_value

    left        = 0
    right       = size(list) + 1

    ! Take the middle element as pivot
    pivot_ix    = size(list) / 2
    pivot_value = list(pivot_ix)%var_name

    do while (left < right)

       right = right - 1
       do while (lgt(list(right)%var_name, pivot_value))
          right = right - 1
       end do

       left = left + 1
       do while (lgt(pivot_value, list(left)%var_name))
          left = left + 1
       end do

       if (left < right) then
          temp = list(left)
          list(left) = list(right)
          list(right) = temp
       end if
    end do

    if (left == right) then
       marker = left + 1
    else
       marker = left
    end if
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
  SUBROUTINE QFYAML_CleanUp(yml)
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
    yml%sorted   = .false.
    yml%num_vars = 0

    IF ( ALLOCATED( yml%vars ) ) DEALLOCATE( yml%vars )

  END SUBROUTINE QFYAML_CleanUp
!EOC
!------------------------------------------------------------------------------
! QFYAML: Bob Yantosca | yantosca@seas.harvard.edu | Apr 2020
! Based on existing package https://github.com/jannisteunissen/config_fortran
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE:
!
! !DESCRIPTION: Get the size of a variable

!\\
!\\
! !INTERFACE:
!
  SUBROUTINE QFYAML_Get_Size(yml, var_name, res)
!
! !INPUT PARAMETERS:
!
    TYPE(QFYAML_t),   INTENT(IN)  :: yml
    CHARACTER(LEN=*), INTENT(IN)  :: var_name
!
! !OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(OUT) :: res
!
! !REVISION HISTORY:
!  15 Apr 2020 - R. Yantosca - Initial version
!  See the subsequent Git history with the gitk browser!
!EOP
!------------------------------------------------------------------------------
!
    INTEGER :: ix

    CALL get_var_index(yml, var_name, ix)
    IF (ix /= -1) THEN
       res = yml%vars(ix)%var_size
    ELSE
       res = -1
       CALL handle_error("QFYAML_get_size: variable ["//var_name//"] not found")
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
  SUBROUTINE QFYAML_Get_Type(yml, var_name, res)
!
! !INPUT PARAMETERS:
!
    TYPE(QFYAML_t),   INTENT(IN)  :: yml
    CHARACTER(LEN=*), INTENT(IN)  :: var_name
!
! !OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(OUT) :: res
!
! !REVISION HISTORY:
!  15 Apr 2020 - R. Yantosca - Initial version
!  See the subsequent Git history with the gitk browser!
!EOP
!------------------------------------------------------------------------------
!
    INTEGER :: ix

    CALL get_var_index(yml, var_name, ix)

    IF (ix /= -1) THEN
       res = yml%vars(ix)%var_type
    ELSE
       res = -1
       CALL handle_error("QFYAML_get_type: variable ["//var_name//"] not found")
    ENDIF

  END SUBROUTINE QFYAML_get_type
!EOC
!############################################################################
!### HERE FOLLOWS OVERLOADED MODULE PROCEDURES.  THESE ARE SIMPLE SO
!### WE WILL OMIT ADDING SUBROUTINE HEADERS
!############################################################################

  ! Add a configuration variable with a real value
  SUBROUTINE Add_Real(yml, var_name, real_data, comment)
    TYPE(QFYAML_t),   INTENT(INOUT) :: yml
    CHARACTER(LEN=*), INTENT(IN   ) :: var_name
    REAL(yp),         INTENT(IN   ) :: real_data
    CHARACTER(LEN=*), INTENT(IN   ) :: comment
    INTEGER                         :: ix

    CALL prepare_store_var(yml, var_name, QFYAML_real_type, 1, comment, ix)

    if (yml%vars(ix)%stored_data /= unstored_data_string) then
       CALL read_variable(yml%vars(ix))
    else
       yml%vars(ix)%real_data(1) = real_data
    end if
  END SUBROUTINE Add_Real

  ! Add a configuration variable with an array of type REAL
  SUBROUTINE Add_Real_Array(yml, var_name, real_data, comment, dynamic_size)
    TYPE(QFYAML_t),   INTENT(INOUT) :: yml
    CHARACTER(LEN=*), INTENT(IN   ) :: var_name
    REAL(yp),         INTENT(IN   ) :: real_data(:)
    CHARACTER(LEN=*), INTENT(IN   ) :: comment
    LOGICAL,          OPTIONAL      :: dynamic_size
    INTEGER                         :: ix

    CALL prepare_store_var(yml, var_name, QFYAML_real_type, &
         size(real_data), comment, ix, dynamic_size)

    if (yml%vars(ix)%stored_data /= unstored_data_string) then
       CALL read_variable(yml%vars(ix))
    else
       yml%vars(ix)%real_data = real_data
    end if
  END SUBROUTINE Add_Real_Array

  ! Add a configuration variable with an INTEGER value
  SUBROUTINE Add_Int(yml, var_name, int_data, comment)
    TYPE(QFYAML_t),   INTENT(INOUT) :: yml
    CHARACTER(LEN=*), INTENT(IN   ) :: var_name
    INTEGER,          INTENT(IN   ) :: int_data
    CHARACTER(LEN=*), INTENT(IN   ) :: comment
    INTEGER                         :: ix

    CALL prepare_store_var(yml, var_name, QFYAML_INTEGER_type, 1, comment, ix)

    if (yml%vars(ix)%stored_data /= unstored_data_string) then
       CALL read_variable(yml%vars(ix))
    else
       yml%vars(ix)%int_data(1) = int_data
    end if
  END SUBROUTINE Add_Int

  ! Add a configuration variable with an array of type INTEGER
  SUBROUTINE Add_Int_Array(yml, var_name, int_data, comment, dynamic_size)
    TYPE(QFYAML_t),   INTENT(INOUT) :: yml
    CHARACTER(LEN=*), INTENT(IN   ) :: var_name
    INTEGER,          INTENT(IN   ) :: int_data(:)
    CHARACTER(LEN=*), INTENT(IN   ) :: comment
    LOGICAL,          OPTIONAL      :: dynamic_size
    INTEGER                         :: ix

    CALL prepare_store_var(yml, var_name, QFYAML_INTEGER_type, &
         size(int_data), comment, ix, dynamic_size)

    if (yml%vars(ix)%stored_data /= unstored_data_string) then
       CALL read_variable(yml%vars(ix))
    else
       yml%vars(ix)%int_data = int_data
    end if
  END SUBROUTINE Add_Int_Array

  ! Add a configuration variable with an character value
  SUBROUTINE Add_String(yml, var_name, char_data, comment)
    TYPE(QFYAML_t),   INTENT(INOUT) :: yml
    CHARACTER(LEN=*), INTENT(IN   ) :: var_name
    CHARACTER(LEN=*), INTENT(IN   ) :: comment
    CHARACTER(LEN=*), INTENT(IN   ) :: char_data
    INTEGER                         :: ix

    CALL prepare_store_var(yml, var_name, QFYAML_string_type, 1, comment, ix)
    if (yml%vars(ix)%stored_data /= unstored_data_string) then
       CALL read_variable(yml%vars(ix))
    else
       yml%vars(ix)%char_data(1) = char_data
    end if
  END SUBROUTINE Add_String

  ! Add a configuration variable with an array of type character
  SUBROUTINE Add_String_Array(yml, var_name, char_data, &
                              comment, dynamic_size)
    TYPE(QFYAML_t),   INTENT(INOUT) :: yml
    CHARACTER(LEN=*), INTENT(IN   ) :: var_name
    CHARACTER(LEN=*), INTENT(IN   ) :: comment
    CHARACTER(LEN=*), INTENT(IN   ) :: char_data(:)
    LOGICAL,          OPTIONAL      :: dynamic_size
    INTEGER                         :: ix

    CALL prepare_store_var(yml, var_name, QFYAML_string_type, &
         size(char_data), comment, ix, dynamic_size)

    if (yml%vars(ix)%stored_data /= unstored_data_string) then
       CALL read_variable(yml%vars(ix))
    else
       yml%vars(ix)%char_data = char_data
    end if
  END SUBROUTINE Add_String_Array

  ! Add a configuration variable with an logical value
  SUBROUTINE Add_Bool(yml, var_name, bool_data, comment)
    TYPE(QFYAML_t),   INTENT(INOUT) :: yml
    CHARACTER(LEN=*), INTENT(IN   ) :: var_name
    LOGICAL,          INTENT(IN   ) :: bool_data
    CHARACTER(LEN=*), INTENT(IN   ) :: comment
    INTEGER                         :: ix

    CALL prepare_store_var(yml, var_name, QFYAML_bool_type, 1, comment, ix)

    if (yml%vars(ix)%stored_data /= unstored_data_string) then
       CALL read_variable(yml%vars(ix))
    else
       yml%vars(ix)%bool_data(1) = bool_data
    end if
  END SUBROUTINE Add_Bool

  ! Add a configuration variable with an array of type LOGICAL
  SUBROUTINE Add_Bool_Array(yml, var_name, bool_data, &
                            comment, dynamic_size)
    TYPE(QFYAML_t),   INTENT(INOUT) :: yml
    CHARACTER(LEN=*), INTENT(IN   ) :: var_name
    LOGICAL,          INTENT(IN   ) :: bool_data(:)
    CHARACTER(LEN=*), INTENT(IN   ) :: comment
    LOGICAL,          OPTIONAL      :: dynamic_size
    INTEGER                         :: ix

    CALL prepare_store_var(yml, var_name, QFYAML_bool_type, &
         size(bool_data), comment, ix, dynamic_size)

    if (yml%vars(ix)%stored_data /= unstored_data_string) then
       CALL read_variable(yml%vars(ix))
    else
       yml%vars(ix)%bool_data = bool_data
    end if
  END SUBROUTINE Add_Bool_Array

  ! Get a real array of a given name
  SUBROUTINE Get_Real_Array(yml, var_name, real_data)
    TYPE(QFYAML_t),   INTENT(INOUT) :: yml
    CHARACTER(LEN=*), INTENT(IN   ) :: var_name
    REAL(yp),         INTENT(INOUT) :: real_data(:)
    INTEGER                         :: ix

    CALL prepare_get_var(yml, var_name, QFYAML_real_type, &
         size(real_data), ix)
    real_data = yml%vars(ix)%real_data
  END SUBROUTINE Get_Real_Array

  ! Get a INTEGER array of a given name
  SUBROUTINE Get_Int_Array(yml, var_name, int_data)
    TYPE(QFYAML_t),   INTENT(INOUT) :: yml
    CHARACTER(LEN=*), INTENT(IN   ) :: var_name
    INTEGER,          INTENT(INOUT) :: int_data(:)
    INTEGER                         :: ix

    CALL prepare_get_var(yml, var_name, QFYAML_integer_type, &
         size(int_data), ix)
    int_data    = yml%vars(ix)%int_data
  END SUBROUTINE Get_Int_Array

  ! Get a character array of a given name
  SUBROUTINE Get_String_Array(yml, var_name, char_data)
    TYPE(QFYAML_t),   INTENT(INOUT) :: yml
    CHARACTER(LEN=*), INTENT(IN   )    :: var_name
    CHARACTER(LEN=*), INTENT(INOUT) :: char_data(:)
    INTEGER                         :: ix

    CALL prepare_get_var(yml, var_name, QFYAML_string_type, &
         size(char_data), ix)
    char_data = yml%vars(ix)%char_data
  END SUBROUTINE Get_String_Array

  ! Get a LOGICAL array of a given name
  SUBROUTINE Get_Bool_Array(yml, var_name, bool_data)
    TYPE(QFYAML_t),   INTENT(INOUT) :: yml
    CHARACTER(LEN=*), INTENT(IN   ) :: var_name
    LOGICAL,          INTENT(INOUT) :: bool_data(:)
    INTEGER                         :: ix

    CALL prepare_get_var(yml, var_name, QFYAML_bool_type, &
         size(bool_data), ix)
    bool_data = yml%vars(ix)%bool_data
  END SUBROUTINE Get_Bool_Array

  ! Get a real value of a given name
  SUBROUTINE Get_Real(yml, var_name, res)
    TYPE(QFYAML_t),   INTENT(INOUT) :: yml
    CHARACTER(LEN=*), INTENT(IN   ) :: var_name
    REAL(yp),         INTENT(OUT  ) :: res
    INTEGER                         :: ix

    CALL prepare_get_var(yml, var_name, QFYAML_real_type, 1, ix)
    res = yml%vars(ix)%real_data(1)
  END SUBROUTINE Get_Real

  ! Get a INTEGER value of a given name
  SUBROUTINE Get_Int(yml, var_name, res)
    TYPE(QFYAML_t),   INTENT(INOUT) :: yml
    CHARACTER(LEN=*), INTENT(IN   ) :: var_name
    INTEGER,          INTENT(INOUT) :: res
    INTEGER                         :: ix

    CALL prepare_get_var(yml, var_name, QFYAML_integer_type, 1, ix)
    res = yml%vars(ix)%int_data(1)
  END SUBROUTINE Get_Int

  ! Get a LOGICAL value of a given name
  SUBROUTINE Get_Bool(yml, var_name, res)
    TYPE(QFYAML_t),   INTENT(INOUT) :: yml
    CHARACTER(LEN=*), INTENT(IN   ) :: var_name
    LOGICAL,          INTENT(OUT  ) :: res
    INTEGER                         :: ix

    CALL prepare_get_var(yml, var_name, QFYAML_bool_type, 1, ix)
    res = yml%vars(ix)%bool_data(1)
  END SUBROUTINE Get_Bool

  ! Get a character value of a given name
  SUBROUTINE Get_String(yml, var_name, res)
    TYPE(QFYAML_t),   INTENT(INOUT) :: yml
    CHARACTER(LEN=*), INTENT(IN   ) :: var_name
    CHARACTER(LEN=*), INTENT(OUT  ) :: res
    INTEGER                         :: ix

    CALL prepare_get_var(yml, var_name, QFYAML_string_type, 1, ix)
    res = yml%vars(ix)%char_data(1)
  END SUBROUTINE Get_String

  ! Get or add a real array of a given name
  SUBROUTINE Add_Get_Real_Array(yml, var_name, real_data, &
                                comment, dynamic_size)
    TYPE(QFYAML_t),   INTENT(INOUT) :: yml
    CHARACTER(LEN=*), INTENT(IN   ) :: var_name
    REAL(yp),         INTENT(INOUT) :: real_data(:)
    CHARACTER(LEN=*), INTENT(IN   ) :: comment
    LOGICAL,          OPTIONAL      :: dynamic_size

    CALL Add_Real_Array(yml, var_name, real_data, comment, dynamic_size)
    CALL Get_Real_Array(yml, var_name, real_data)
  END SUBROUTINE Add_Get_Real_Array

  ! Get or add a INTEGER array of a given name
  SUBROUTINE Add_Get_Int_Array(yml, var_name, int_data, &
                               comment, dynamic_size)
    TYPE(QFYAML_t),   INTENT(INOUT) :: yml
    CHARACTER(LEN=*), INTENT(IN   ) :: var_name
    INTEGER,          INTENT(INOUT) :: int_data(:)
    CHARACTER(LEN=*), INTENT(IN   ) :: comment
    LOGICAL,          OPTIONAL      :: dynamic_size

    CALL Add_Int_Array(yml, var_name, int_data, comment, dynamic_size)
    CALL Get_Int_Array(yml, var_name, int_data)
  END SUBROUTINE Add_Get_Int_Array

  ! Get or add a character array of a given name
  SUBROUTINE Add_Get_String_Array(yml, var_name, char_data, &
                                  comment, dynamic_size)
    TYPE(QFYAML_t),   INTENT(INOUT) :: yml
    CHARACTER(LEN=*), INTENT(IN   ) :: var_name
    CHARACTER(LEN=*), INTENT(INOUT) :: char_data(:)
    CHARACTER(LEN=*), INTENT(IN   ) :: comment
    LOGICAL,          optional      :: dynamic_size

    CALL add_string_Array(yml, var_name, char_data, comment, dynamic_size)
    CALL Get_String_Array(yml, var_name, char_data)
  END SUBROUTINE Add_Get_String_Array

  ! Get or add a LOGICAL array of a given name
  SUBROUTINE Add_Get_Bool_Array(yml, var_name, bool_data, &
                                comment, dynamic_size)
    TYPE(QFYAML_t),   INTENT(INOUT) :: yml
    CHARACTER(LEN=*), INTENT(IN   ) :: var_name
    LOGICAL,          INTENT(INOUT) :: bool_data(:)
    CHARACTER(LEN=*), INTENT(IN   ) :: comment
    LOGICAL,          OPTIONAL      :: dynamic_size

    CALL Add_Bool_Array(yml, var_name, bool_data, comment, dynamic_size)
    CALL Get_Bool_Array(yml, var_name, bool_data)
  END SUBROUTINE Add_Get_Bool_Array

  ! Get or add a real value of a given name
  SUBROUTINE Add_Get_Real(yml, var_name, real_data, comment)
    TYPE(QFYAML_t),   INTENT(INOUT) :: yml
    CHARACTER(LEN=*), INTENT(IN   ) :: var_name
    REAL(yp),         INTENT(INOUT) :: real_data
    CHARACTER(LEN=*), INTENT(IN   ) :: comment

    CALL Add_Real(yml, var_name, real_data, comment)
    CALL Get_Real(yml, var_name, real_data)
  END SUBROUTINE Add_Get_Real

  ! Get or add a INTEGER value of a given name
  SUBROUTINE Add_Get_Int(yml, var_name, int_data, comment)
    TYPE(QFYAML_t),   INTENT(INOUT) :: yml
    CHARACTER(LEN=*), INTENT(IN   ) :: var_name
    INTEGER,          INTENT(INOUT) :: int_data
    CHARACTER(LEN=*), INTENT(IN   ) :: comment

    CALL Add_Int(yml, var_name, int_data, comment)
    CALL Get_Int(yml, var_name, int_data)
  END SUBROUTINE Add_Get_Int

  ! Get or add a LOGICAL value of a given name
  SUBROUTINE Add_Get_Bool(yml, var_name, bool_data, comment)
    TYPE(QFYAML_t),   INTENT(INOUT) :: yml
    CHARACTER(LEN=*), INTENT(IN   ) :: var_name
    LOGICAL,          INTENT(INOUT) :: bool_data
    CHARACTER(LEN=*), INTENT(IN   ) :: comment

    CALL Add_Bool(yml, var_name, bool_data, comment)
    CALL Get_Bool(yml, var_name, bool_data)
  END SUBROUTINE Add_Get_Bool

  ! Get a character value of a given name
  SUBROUTINE Add_Get_String(yml, var_name, string_data, comment)
    TYPE(QFYAML_t),   INTENT(INOUT) :: yml
    CHARACTER(LEN=*), INTENT(IN   ) :: var_name
    CHARACTER(LEN=*), INTENT(INOUT) :: string_data
    CHARACTER(LEN=*), INTENT(IN   ) :: comment

    CALL add_string(yml, var_name, string_data, comment)
    CALL Get_String(yml, var_name, string_data)
  END SUBROUTINE Add_Get_String

  ! Get or add a real array of a given name
  SUBROUTINE Update_Real_Array(yml, var_name, real_data)
    TYPE(QFYAML_t),   INTENT(INOUT) :: yml
    CHARACTER(LEN=*), INTENT(IN   ) :: var_name
    REAL(yp),         INTENT(INOUT) :: real_data(:)
    INTEGER                         :: ix

    CALL Get_Var_Index(yml, var_name, ix)
    IF ( ix > 0 ) THEN 
       yml%vars(ix)%real_data = real_data
       real_data = yml%vars(ix)%real_data
    ENDIF
  END SUBROUTINE Update_Real_Array

  ! Get or add a INTEGER array of a given name
  SUBROUTINE Update_Int_Array(yml, var_name, int_data)
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

  ! Get or add a character array of a given name
  SUBROUTINE Update_String_Array(yml, var_name, char_data)
    TYPE(QFYAML_t),   INTENT(INOUT) :: yml
    CHARACTER(LEN=*), INTENT(IN   ) :: var_name
    CHARACTER(LEN=*), INTENT(INOUT) :: char_data(:)
    INTEGER                         :: ix

    CALL Get_Var_Index(yml, var_name, ix)
    IF ( ix > 0 ) THEN 
       yml%vars(ix)%char_data = char_data
       char_data = yml%vars(ix)%char_data
    ENDIF
  END SUBROUTINE Update_String_Array

  ! Get or add a LOGICAL array of a given name
  SUBROUTINE Update_Bool_Array(yml, var_name, bool_data)
    TYPE(QFYAML_t),   INTENT(INOUT) :: yml
    CHARACTER(LEN=*), INTENT(IN   ) :: var_name
    LOGICAL,          INTENT(INOUT) :: bool_data(:)
    INTEGER                         :: ix

    CALL Get_Var_Index(yml, var_name, ix)
    IF ( ix > 0 ) THEN 
       yml%vars(ix)%bool_data = bool_data
       bool_data = yml%vars(ix)%bool_data
    ENDIF
  END SUBROUTINE Update_Bool_Array

  ! Get or add a real value of a given name
  SUBROUTINE Update_Real(yml, var_name, real_data)
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

  ! Get or add a INTEGER value of a given name
  SUBROUTINE Update_Int(yml, var_name, int_data)
    TYPE(QFYAML_t),   INTENT(INOUT) :: yml
    CHARACTER(LEN=*), INTENT(IN   ) :: var_name
    INTEGER,          INTENT(INOUT) :: int_data
    INTEGER                         :: ix

    CALL Get_Var_Index(yml, var_name, ix)
    IF ( ix > 0 ) THEN 
       yml%vars(ix)%int_data(1) = int_data
       int_data = yml%vars(ix)%int_data(1)
    ENDIF
  END SUBROUTINE Update_Int

  ! Get or add a LOGICAL value of a given name
  SUBROUTINE Update_Bool(yml, var_name, bool_data)
    TYPE(QFYAML_t),   INTENT(INOUT) :: yml
    CHARACTER(LEN=*), INTENT(IN   ) :: var_name
    LOGICAL,          INTENT(INOUT) :: bool_data
    INTEGER                         :: ix

    CALL Get_Var_Index(yml, var_name, ix)
    IF ( ix > 0 ) THEN 
       yml%vars(ix)%bool_data(1) = bool_data
       bool_data = yml%vars(ix)%bool_data(1)
    ENDIF
  END SUBROUTINE Update_Bool

  ! Get a character value of a given name
  SUBROUTINE Update_String(yml, var_name, string_data)
    TYPE(QFYAML_t),   INTENT(INOUT) :: yml
    CHARACTER(LEN=*), INTENT(IN   ) :: var_name
    CHARACTER(LEN=*), INTENT(INOUT) :: string_data
    INTEGER                         :: ix
    
    CALL Get_Var_Index(yml, var_name, ix)
    IF ( ix > 0 ) THEN 
       yml%vars(ix)%char_data(1) = string_data
       string_data = yml%vars(ix)%char_data(1)
    ENDIF
  END SUBROUTINE Update_String

END MODULE QFYAML_Mod
