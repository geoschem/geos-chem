!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: hco_emisll_mod 
!
! !DESCRIPTION: Module HCO\_EMISLL\_MOD contains routines for creating / 
!  editing emission fields. Each emission field contains all information
!  associated with the emission field, which includes direct field 
!  information (dimensions, temporal resolution, units, file source, etc.), 
!  as well as information related to the calculation of the effective
!  Geos-Chem emissions, such as the field priority and hierarchy, the
!  tracer ID, associated scale factors, and more.\\ 
!  This module only contains the basic routines such as adding, editing or
!  removing fields, but not processing routines such as calculating the
!  current emissions, setting the correct time indices or associating
!  scale factors to emission fields. See NG\_EMIS\_CALC\_MOD, 
!  NG\_EMIS\_TIME\_MOD and NG\_EMIS\_SCALE\_MOD for more details on 
!  processing routines.\\
! \\
! ============================================================================
! CONCEPTS
! ============================================================================
! \\
! (1) Emission fields
! All emission field information such as the field name, the emission
! data, associated scale factors, the tracer ID, the emission field
! category and priority, etc., are stored in an emission field structure
! structure (derived type EMISCONT in NG\_EMIS\_TYPE\_MOD.F).
! All containers are linked to each other through the EMISSIONS linked
! list.
! (2) EMISSIONS linked list
! The core of the module is the EMISSIONS linked list, consisting of
! emission fields of type 'EMISCONT':\\
! \\
! EMISSIONS   --------- NEXT    ---------  NEXT   ---------  NEXT
! ---------> | FIELD 3 | ----> | FIELD 2 | ----> | FIELD 1 | ----> NULL()
!             ---------         ---------         ---------
! \\
! The pointer EMISSIONS always points to the head of the linked list.
! The last container in the list ('FIELD 1') always points to a nullified 
! 'EMISCONT' structure.\\
! The elements of the linked list can represent 'real' emission data
! (base emissions, e.g. NOx emissions for year 2000), scale factors (e.g.
! time of day scale factors), or masks (e.g. NEI mask). The emission fields 
! are 4-dimensional, covering the 3 spatial dimensions and time.\\
! \\
! (3) Emission calculations
! The emissions linked list contains all information required to
! calculate the final emission array. Calculations are done through
! multiplication of base emissions and scale factors. The scale factors
! for each base field are specified in SclID.
! See NG\_EMIS\_CALC\_MOD.F90 for more details on the emission
! calculation.\\
! \\
! (4) 4th dimension
! The 4th dimension handles the temporal variability of an emission
! field. For instance, diurnal scale factors can all be stored in the 
! same emission array (so that the 4th dimension becomes 24), and the 
! currently valid 3D-array is then automatically recognized based upon 
! the specification of TempRes (being 'HOURLY' in this example).\\
! See NG\_EMIS\_TIME\_MOD.F90 for more details.\\
! \\
! (5) Vertical extension
! Emissions can extend over more than one level. However, the emission
! levels have to be the same as the GEOS-Chem levels since vertical
! interpolation is not yet supported.\\
! Also, the levels are always assumed to start at the surface, so if the
! emission array has four vertical levels, these emissions will be
! emitted to the four lowest GEOS-Chem levels.   
! \\
! (6) Emission field categories and hierarchies
! Emission field categories and hierarchies define the way emissions are 
! taken into account. Generally speaking, emissions of higher hierarchy 
! overwrite lower-hierarchy emission values within the same emission 
! category (e.g if RETRO emissions have hierarchy 10, and GEIA emissions 
! hierarchy 5, RETRO emissions will be used instead of the GEIA
! emissions wherever RETRO emissions are available).\\ 
! The so defined emissions per category are then added up, i.e. fields
! of different categories never overwrite each other.
! If emissions shall be added to the existing data (e.g. lightning NOx 
! emissions), the priority can be set to 0 which will always add the 
! emissions within the specified category. This is equivalent to
! defining an own category for this field.
! If two emission fields have the same hierarchy (e.g. RETRO_BUTANES 
! and RETRO_PENTANES), these emissions won't overwrite each other but 
! both will be added to the emission array (unless there exist another
! field with higher hierarchy).\\ 
! This concept requires that the fields are sorted according to 
! their categories and hierarchies, with increasing hierarchies 
! within the same category. Fields with hierarchy zero are added at the 
! end of each category. This is done by subroutine Add2EmisLL, which places
! a new field in the proper slot within the EMISSIONS linked list.\\
! \\
! (7) Negative values / zeros
! Negative emissions are not supported and always ignored, i.e. grid
! boxes with negative emissions are not considered.
! Zeros are supported.
! \\
! (8) Scalar fields
! Scalar fields are uniform fields, i.e. fields that contain only one
! single value. The flag IsScalar has to be set for these fields,
! otherwise the code will crash when attempting to regrid the field to
! proper grid dimensions.
! \\
! (9) Scale factors
! Base emission data can be linked to other emission fields (scale
! factors) through the Scl_ID vector. This vector contains just the
! field IDs of the scale fields of interest. The parameter SclMax (see
! ng_emis_type_mod.F90) determines the maximum number of scale fields. 
! \\
! (10) Field IDs
! Each field has a unique ID, which makes it possible to easily
! identify / search for a field. The field ID starts with ID = 1 
! and increases by 1 for every new field.\\
! \\
! (11) IDList
! Accessing the fields through Find_EmisCont can become inefficient as
! the search always starts from the first element of the linked list and
! walks through it until the desired field is found. To avoid the 
! usage of Find_EmisCont too often, a 'quick access list' (IDList) is 
! created once that all fields are defined. The IDList is a vector
! that contains pointers to all emission fields. The position of the
! pointers corresponds to the field ID, i.e. IDList(3) will point
! to the field with ID 3.\\
! The IDList is created in subroutine CALC_EMIS. It is redefined every
! time an emission field has been added/removed from the emissions
! linked list. The two variables nnEmisCont and nnIDList are used to
! determine if this is the case or not. Typically, the IDList is only
! created at the first call.
! \\
! (12) Masks
! Masks are treated somewhat different than other emission fields in a
! sense that masks can be added together and will only be applied after
! all other emission calculations. Mask fields are identified through
! the flag IsMask.\\
! Currently, mask values are assumed to be binary, i.e. they are either
! 0 or 1. All grid boxes where the MASK is zero will be ignored for the
! emission calculation. Multiple masks can be used, in which case all 
! masks are added up and all grid boxes where at least one mask applies 
! are considered.\\
! \\
! (13) Field updates
! The numeric vector FLDTIME contains the year, month, day and hour of
! the emission field (with undefined time stamps being -1). The vector 
! UPDATE determine if (and how often) the field is updated.
! UPDATE is also of length four and contains the annual, monthly, daily
! and hourly time increases before a field will be updated. Hence, if
! UPDATE(2) is set to 1, the emission field will be updated every month.
! By default, all UPDATE values are set to zero, i.e. emissions are
! never updated.
! See NG\_EMIS\_UPDATE\_MOD.F90 for more details on the automatic field 
! updating. 
! \\
! (14) Field operators
! The operator flag can be used to perform special operations on the 
! data values before they are used. The following operator flags are
! currently defined: (1) Takes the inverse of the value (1/X) if X is
! non-zero; (2) Squares the value before using it (X**2); (3) Mirrors
! the value (1-X). This is only supported for mask fields. 
!
! Module Routines:
! ============================================================================
!
! - PUBLIC:
!
!  ( 1) Add_EmisCont        : Creates a new emission field structure or 
!                             updates existing content.
!  ( 2) Find_EmisCont       : Finds a field and points to it. It's 
!                             possible to search either for the field
!                             name or the field ID.
!  ( 3) Cleanup_EmisLL      : To remove all emission fields
!  ( 4) DEFAULT_EmisContArgs     : Defines the arguments for subroutine
!                             Add_EmisCont, incl. defaults.
!  ( 5) CREATE_IDList       : Creates the IDList with pointers to all
!                             emission fields.
!  ( 6) Show_EmisLL         : Promts the container names of all emission
!                             containers in the EMISSIONS linked list.
!  ( 7) Add_SclID           : Assigns a scale factor field ID to a
!                             (base) emission field. 
!  ( 9) GET_nnEmisCont      : Returns the number of containers currently
!                             defined in the linked list
!  (10) GET_nnIDList        : Returns the number of elements of IDList
!
! - PRIVATE:
!
!  (12) Add2EmisLL          : Adds a field to the fields linked list
!  (13) REMOVE_FIELD        : Removes an emission field.
! \\
! !INTERFACE: 
!
      MODULE HCO_EMISLL_MOD 
!
! !USES:
!
      USE HCO_ERROR_MOD

      IMPLICIT NONE

      PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
      PUBLIC :: Add_EmisCont
      PUBLIC :: Find_EmisCont
      PUBLIC :: Cleanup_EmisLL
      PUBLIC :: Create_IDList 
      PUBLIC :: Show_EmisLL
      PUBLIC :: Add_SclID
      PUBLIC :: Get_nnEmisCont
      PUBLIC :: Get_nnIDList
      PUBLIC :: Pnt2EmisCont
 
      ! Public types/variables
      PUBLIC :: EmisCont
      PUBLIC :: SclMax
!
! !PRIVATE:
!
      ! Private integers
      ! nnEmisCont is the number of defined containers in the emissions
      ! linked list, and nnIDList is the number of entries in the IDList. 
      ! Write out both numbers so that we can figure out if we have to 
      ! recalculate the IDList.
      PRIVATE :: nnEmisCont
      PRIVATE :: nnIDList
      PRIVATE :: IDList

      ! Private functions/subroutines
      PRIVATE :: Find_EmisCont_ID
      PRIVATE :: Find_EmisCont_Name
      PRIVATE :: Add2EmisLL
      PRIVATE :: InitEmisCont
      PRIVATE :: CleanupEmisCont
!
! !REVISION HISTORY:
!  04 Dec 2012 - C. Keller: Initialization
!  22 Aug 2013 - C. Keller: Renamed from NG_EMIS_MOD to HCO_EMISLL_MOD
!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !MODULE TYPES:
!
      ! Emission container
      TYPE EmisCont
         INTEGER,           POINTER  :: DataType
         INTEGER,           POINTER  :: ExtNr
         CHARACTER(LEN=31), POINTER  :: cName
         CHARACTER(LEN=31), POINTER  :: TempRes
         REAL*8,            POINTER  :: Array(:,:,:,:)
         INTEGER,           POINTER  :: cID
         INTEGER,           POINTER  :: scalID
         INTEGER,           POINTER  :: TrcID
         INTEGER,           POINTER  :: Cat
         INTEGER,           POINTER  :: Hier
         INTEGER,           POINTER  :: TimeSlc(:)
         LOGICAL,           POINTER  :: IsScalar
         INTEGER,           POINTER  :: Oper
         INTEGER,           POINTER  :: Scal_cID(:)
         TYPE(EmisCont),    POINTER  :: NextCont
      END TYPE EmisCont

      ! For the IDList
      TYPE IDListPnt
         TYPE(EmisCont), POINTER  :: PNT ! Pointer to field
      ENDTYPE IDListPnt
!
! !MODULE ARGUMENTS:
!
      TYPE(IDListPnt), POINTER      :: IDList(:)  => NULL()
      INTEGER                       :: nnEmisCont = 0
      INTEGER                       :: nnIDList   = 0
!
! !MODULE PARAMETER:
!
      ! Maximum possible number of scale factors per field
      INTEGER, PARAMETER :: SclMax = 10
!
! !MODULE INTERFACES:
!
      INTERFACE Find_EmisCont
         MODULE PROCEDURE Find_EmisCont_Name
         MODULE PROCEDURE Find_EmisCont_ID
      END INTERFACE

      !----------------------------------------------------------------
      ! MODULE ROUTINES follow below
      !----------------------------------------------------------------

      CONTAINS
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: Add_EmisCont
!
! !DESCRIPTION: Subroutine Add\_EmisCont creates a new emission container and 
!  sets the corresponding field parameter.
!  If the OVERWRITE argument is enabled, existing data will be
!  overwritten using the variables in ARGS. Note that in this case, previously 
!  set scale factor indices will be maintained!
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE Add_EmisCont( aIR,     HcoState, DataType, cName,   &
                               Array,   scalID,    TrcID,    Oper,    &
                               Cat,     Hier,      TempRes,  EmisLL,  &
                               cID,     ExtNr,     RC )
!
! !USES:
!
      USE HCO_TIME_MOD,     ONLY : Set_tSlc
      USE HCO_TYPE_MOD,     ONLY : HCO_State 
!
! !ARGUMENTS:
!
      ! Container arguments
      LOGICAL,           INTENT(IN)         :: aIR 
      TYPE(HCO_State),   POINTER            :: HcoState
      INTEGER,           INTENT(IN), TARGET :: DataType
      CHARACTER(LEN=31), INTENT(IN), TARGET :: cName
      REAL*8,            POINTER            :: Array(:,:,:,:)
      INTEGER,           INTENT(IN), TARGET :: scalID
      INTEGER,           INTENT(IN), TARGET :: TrcID
      INTEGER,           INTENT(IN), TARGET :: Oper
      INTEGER,           INTENT(IN), TARGET :: Cat
      INTEGER,           INTENT(IN), TARGET :: Hier
      INTEGER,           INTENT(IN), TARGET :: ExtNr
      CHARACTER(LEN=31), INTENT(IN)         :: TempRes

      ! Emissions linked list    
      TYPE(EmisCont), POINTER               :: EmisLL

      ! Field ID assigned to this field
      INTEGER,       INTENT(OUT)            :: cID

      ! Error status
      INTEGER,       INTENT(INOUT)          :: RC  
!
! !REVISION HISTORY:
!  04 Dec 2012 - C. Keller - Initialization
!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL ARGUMENTS:
!
      INTEGER                                 :: NLON, NLAT, NLEV, NTIME
      TYPE(EmisCont), POINTER                 :: NewCnt => NULL()
      LOGICAL                                 :: FOUND, NEW, VERBOSE 
      CHARACTER(LEN=255)                      :: MSG, LOC

      !======================================================================
      ! Add_EmisCont begins here!
      !======================================================================

      ! Enter
      LOC = 'Add_EmisCont (HCO_EMISLL_MOD.F90)'

      ! Set verbose flag
      VERBOSE = HcoState%verbose

      ! Check if field already exists. If this is the case and the
      ! OVERWRITE option has been disabled, return with status = 2.
      ! If OVERWRITE is on, set status to 1 and proceed normally.
      CALL FIND_EMISCONT ( TRIM(cName), EmisLL, FOUND, NewCnt )

      ! Field exists?
      IF ( FOUND ) THEN
         ! Write out FID in any case!
         cID = NewCnt%cID
      ENDIF

      ! For convenience, define variable NEW
      NEW = .NOT. FOUND

      ! Create the new emission structure if this is a new field
      ! (still blank)
      IF ( NEW ) CALL InitEmisCont ( NewCnt )

      ! Get field dimensions
      NLON  = SIZE( Array, 1 )
      NLAT  = SIZE( Array, 2 )
      NLEV  = SIZE( Array, 3 )
      NTIME = SIZE( Array, 4 )

      ! Set field type and extension number
      IF ( NEW ) THEN
         NewCnt%DataType => DataType
         NewCnt%ExtNr    => ExtNr
      ENDIF

      ! Set field data operator
      ! This attribute is not used in the emissions calculations as the
      ! values will be written into the container array with the data
      ! operator already being applied, but storing it along with the
      ! data makes it easier to follow the data flow. 
      IF ( NEW ) THEN
         NewCnt%OPER => Oper
      ENDIF

      ! Add the scalar flag. Fields with one spatial dimension (i.e.
      ! NLON = NLAT = NLEV = 1) are treated as single scalar, which is 
      ! uniformely applied. This helps to reduce memory amount since we
      ! don't have to populate each grid box with the same value.
      IF ( NEW ) ALLOCATE ( NewCnt%IsScalar )
      IF ( NLON == 1 .AND. NLAT == 1 .AND. NLEV == 1 ) THEN
         NewCnt%IsScalar = .TRUE.
      ELSE
         NewCnt%IsScalar = .FALSE.
      ENDIF

      ! Mask field cannot be a scalar field!
      IF ( (NewCnt%DataType==3) .AND. NewCnt%IsScalar ) THEN
         MSG = 'Mask field cannot be scalar: ' // TRIM(cName)
         CALL HCO_ERROR ( MSG, LOC, RC ); RETURN
      ENDIF

      ! Fill data array according to field operator
      IF ( ASSOCIATED(NewCnt%Array) ) DEALLOCATE(NewCnt%Array)
      ALLOCATE ( NewCnt%Array(NLON,NLAT,NLEV,NTIME) )

      ! Operator -1: 1/x
      IF ( NewCnt%Oper == -1 ) THEN
         WHERE ( Array > 0d0 ) 
            NewCnt%Array = 1d0 / Array
         ELSE WHERE
            NewCnt%Array = 0d0
         END WHERE

      ! Operator 1: 1*x
      ELSEIF ( NewCnt%Oper == 1 ) THEN
         NewCnt%Array(:,:,:,:) = Array(:,:,:,:)

      ! Operator 2: x*x
      ELSEIF ( NewCnt%Oper == 2 ) THEN
         NewCnt%Array(:,:,:,:) = Array(:,:,:,:) * Array(:,:,:,:)

      ! Operator 3: 1-x (masks only)
      ELSEIF ( NewCnt%Oper == 3 .AND. (NewCnt%DataType == 3) ) THEN
         NewCnt%Array(:,:,:,:) = 1d0 - Array(:,:,:,:)
         WHERE ( NewCnt%Array < 0d-03 ) 
            NewCnt%Array = 0d0
         END WHERE
      ELSE
         MSG = 'Invalid operator for field ' // TRIM(cName)
         CALL HCO_ERROR ( MSG, LOC, RC ); RETURN
      ENDIF

      ! Set all positive values to 1 if mask flag is set
      IF ( NewCnt%DataType == 3 ) THEN
         WHERE ( NewCnt%Array > 0d0 )
            NewCnt%Array = 1d0
         END WHERE
      ENDIF

      ! Establish the appropriate pointer for the
      ! 4th dimension (temporal resolution) of the field array
      IF ( NEW ) THEN
         ALLOCATE ( NewCnt%TempRes )
         NewCnt%TempRes = TempRes
         CALL Set_tSlc ( NewCnt%TimeSlc, NewCnt%TempRes, NLON )
      ENDIF      

      ! Set scale factor ID
      IF ( NEW ) THEN
         NewCnt%scalID => scalID
      ENDIF 

      ! Update number of emission fields and set container ID
      ! ==> Only if it's a new field. Otherwise, the cID is already
      ! defined.
      IF ( NEW ) THEN
         ALLOCATE ( NewCnt%cID )
         nnEmisCont = nnEmisCont + 1
         cID        = nnEmisCont
         NewCnt%cID = cID
      ENDIF

      ! Set emission field name
      IF ( NEW ) THEN
         NewCnt%cName => cName
      ENDIF

      ! Set emission tracer ID
      IF ( NEW ) THEN
         NewCnt%TrcID => TrcID
      ENDIF

      ! Set field category and hierarchy
      IF ( NEW ) THEN
         NewCnt%Cat  => Cat
         NewCnt%Hier => Hier
      ENDIF

      ! Initialize vector with scale factor IDs
      IF ( NEW ) THEN
         ALLOCATE ( NewCnt%Scal_cID(SclMax) )
         NewCnt%Scal_cID(:) = -1
      ENDIF

      ! For testing, promt the most important field properties
      IF ( aIR .and. verbose ) THEN
         write(6,*) '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
         IF ( NEW ) THEN
            write(6,*) 'New emission field created:'
         ELSE
            write(6,*) 'Emission field updated:'
         ENDIF
         write(6,*) 'Field name  : ', TRIM(NewCnt%cName)
         write(6,*) 'Extension Nr: ', NewCnt%ExtNr
         write(6,*) 'Data type   : ', NewCnt%DataType
         write(6,*) 'Container ID: ', NewCnt%cID
         write(6,*) 'Tracer ID   : ', NewCnt%TrcID
         write(6,*) 'Scal ID     : ', NewCnt%ScalID
         write(6,*) 'Tempres     : ', TRIM(NewCnt%TempRes)
         write(6,*) 'Array dim   : ', NLON, NLAT, NLEV, NTIME
         write(6,*) 'Category    : ', NewCnt%Cat
         write(6,*) 'Hierarchy   : ', NewCnt%Hier
         write(6,*) 'Scalar?       ', NewCnt%IsScalar
         write(6,*) 'Operator    : ', NewCnt%OPER
         write(6,*) 'Emis total1 : ', SUM(NewCnt%Array(:,:,1,:))
         write(6,*) 'Passed array: ', SUM(       Array(:,:,1,:))
         if ( size(NewCnt%Array,3)>1 ) then
            write(6,*) 'Emis total2: ', SUM(NewCnt%Array(:,:,2,:))
         endif
         write(6,*) '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
      ENDIF

      ! Add the new field to the Emissions linked list
      IF ( NEW ) THEN
         CALL Add2EmisLL ( NewCnt, EmisLL, RC )
         IF ( RC /= HCO_SUCCESS ) RETURN
      ENDIF

      ! Leave w/ success
      RC = HCO_SUCCESS

      END SUBROUTINE Add_EmisCont
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: Find_EmisCont_Name
!
! !DESCRIPTION: Subroutine Find\_EmisCont\_Nanem searches for a container 
! name and returns a pointer pointing to this container. 
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE Find_EmisCont_Name ( NME, EmisLL, FOUND, OutCnt )
!
! !ARGUMENTS:
!
      CHARACTER(LEN=*), INTENT(IN )           :: NME
      TYPE(EmisCont),   POINTER               :: EmisLL
      LOGICAL,          INTENT(OUT)           :: FOUND
      TYPE(EmisCont),   POINTER, OPTIONAL     :: OutCnt
!
! !REVISION HISTORY:
!  04 Dec 2012 - C. Keller: Initialization
!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL ARGUMENTS:
!
      TYPE(EMISCONT),   POINTER         :: CurrCnt => NULL() 

      !======================================================================
      ! Find_EmisCont_Name begins here!
      !======================================================================

      ! Initialize
      FOUND  = .FALSE.

      ! Check if an emission field has already been defined or not
      IF ( nnEmisCont == 0 ) RETURN

      ! Make CurrCnt point to first element of the EMISSIONS linked list
      CurrCnt => EmisLL

      ! Loop over EMISSIONS linked list
      DO WHILE ( ASSOCIATED ( CurrCnt ) )

         ! Check if current field is the wanted one
         IF ( TRIM(CurrCnt%cName) == TRIM(NME) ) THEN
            IF ( PRESENT(OutCnt) ) OutCnt => CurrCnt
            FOUND = .TRUE.
            RETURN 
         ENDIF

         ! Advance to next field otherwise
         CurrCnt => CurrCnt%NEXTCONT
      ENDDO

      ! Cleanup
      CurrCnt => NULL()

      END SUBROUTINE Find_EmisCont_Name
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: Find_EmisCont_cID
!
! !DESCRIPTION: Subroutine Find\_EmisCont\_cID searches for a container ID 
! and returns a pointer pointing to this container. 
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE Find_EmisCont_ID ( ID, EmisLL, IsScalID, FOUND, OutCnt)
!
! !ARGUMENTS:
!
      INTEGER,          INTENT(IN )           :: ID
      TYPE(EmisCont),   POINTER               :: EmisLL
      INTEGER,          INTENT(IN )           :: IsScalID
      LOGICAL,          INTENT(OUT)           :: FOUND
      TYPE(EmisCont),   POINTER, OPTIONAL     :: OutCnt
!
! !REVISION HISTORY:
!  04 Dec 2012 - C. Keller: Initialization
!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL ARGUMENTS:
!
      TYPE(EmisCont),   POINTER         :: CurrCnt => NULL() 
      INTEGER                           :: thisID

      !======================================================================
      ! Find_EmisCont_ID begins here!
      !======================================================================

      ! Initialize
      FOUND  = .FALSE.

      ! Check if an emission field has already been defined or not
      IF ( nnEmisCont == 0 ) RETURN

      ! Make CurrCnt point to first element of the EMISSIONS linked list
      CurrCnt => EmisLL

      ! Loop over EMISSIONS linked list
      DO WHILE ( ASSOCIATED ( CurrCnt ) )

         ! Get the current container or original ID
         IF ( IsScalID == 1 ) THEN
            thisID = CurrCnt%scalID
         ELSE
            thisID = CurrCnt%cID
         ENDIF

         ! Check if current field is the wanted one
         IF ( thisID == ID ) THEN
            IF ( PRESENT(OutCnt) ) OutCnt => CurrCnt
            FOUND = .TRUE.
            RETURN 
         ENDIF

         ! Advance to next field otherwise
         CurrCnt => CurrCnt%NEXTCONT
      ENDDO

      ! Cleanup
      CurrCnt => NULL()

      END SUBROUTINE Find_EmisCont_ID
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: Add2EmisLL
!
! !DESCRIPTION: Subroutine Add2EmisLL adds the specified emission
! field to the EMISSIONS linked list. The new field is placed
! according to its category and hierarchy.\\ 
! Base emission fields (hierarchy > 0) are sorted based on tracer ID, 
! category and hierarchy (for fields of same category). 
! Scale fields (= tracer ID < 1) are added at the end of EmisLL.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE Add2EmisLL ( NewCnt, EmisLL, RC )
!
! !ARGUMENTS
!
      TYPE(EmisCont), POINTER       :: NewCnt
      TYPE(EmisCont), POINTER       :: EmisLL
      INTEGER,        INTENT(INOUT) :: RC
!
! !REVISION HISTORY:
!  06 Dec 2012 - C. Keller - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER                   :: NEWCAT, NEWHIR, NEWTRC
      TYPE(EMISCONT), POINTER   :: tmpCnt => NULL()

      !======================================================================
      ! Add2EmisLL begins here!
      !======================================================================

      ! Assume success until otherwise
      RC = HCO_SUCCESS

      ! If this is the first container, we can simply place it at the
      ! beginning of the list.
      IF ( nnEmisCont == 1 ) THEN
         EmisLL => NewCnt
         RETURN ! Leave routine
      ENDIF

      ! Special case where the linked list consists of scale factors
      ! only: In this case, we can place the new container at the head 
      ! of the list no matter what content it has!
      IF ( EmisLL%TrcID <= 0 ) THEN

         ! Place at head and leave
         NewCnt%NEXTCONT => EmisLL 
         EmisLL          => NewCnt
         RETURN 
      ENDIF

      ! Get field tracer ID, category and priority of the new container 
      NEWTRC  = NewCnt%TrcID 
      NEWCAT  = NewCnt%Cat
      NEWHIR  = NewCnt%Hier

      ! Containers are listed with increasing tracer ID. If the current
      ! container has lower tracerID than the first container, just add
      ! it at the beginning of the list.     
      IF ( (NEWTRC > 0) .AND. (NEWTRC < EmisLL%TrcID) ) THEN
         ! Place at head and leave
         NewCnt%NEXTCONT => EmisLL
         EmisLL          => NewCnt
         RETURN
      ENDIF

      ! For the special case that the current field has the same tracer
      ! ID as the first container in the list, but lower category or
      ! same category and lower hierarchy so that the new container has 
      ! to be placed before the first container in the list:
      IF ( NEWTRC == EmisLL%TrcID ) THEN
         IF ( (EmisLL%Cat  >  NEWCAT) .OR.      &
              (EmisLL%Cat  == NEWCAT  .AND.     &
               EmisLL%Hier >  NEWHIR) ) THEN

            ! Place at head and leave
            NewCnt%NEXTCONT => EmisLL
            EmisLL          => NewCnt
            RETURN
         ENDIF
      ENDIF

      ! Set tmpCnt pointer to the beginning of the linked list
      tmpCnt => EmisLL

      ! If the new container contains base data (i.e. has a valid tracer 
      ! ID), we have to move the tmpCnt pointer to the position where the 
      ! next container is one of the following: (a) the first container 
      ! with the same tracer ID as the new container; (b) a container with 
      ! higher tracer ID; (c) scale factors. 
      ! From there, we can determine where to place the container exactly.
      IF ( NEWTRC > 0 ) THEN

         ! Loop over list
         DO WHILE ( ASSOCIATED ( tmpCnt%NextCont ) )
             
            ! Check if next container's tracer ID is higher or if it's a
            ! scale factor, in which case we have to exit.
            IF ( tmpCnt%NextCont%TrcID >  NEWTRC .OR. & 
                 tmpCnt%NextCont%DataType > 1 ) THEN
               EXIT
            ENDIF
 
            ! Check if next container has the same tracer ID but a
            ! higher category or the same category but higher hierarchy,
            ! in which case we have to exit.
            IF ( tmpCnt%NextCont%TrcID == NEWTRC ) THEN
               IF ( tmpCnt%NextCont%Cat  >  NEWCAT ) THEN
                  EXIT
               ENDIF
               IF ( tmpCnt%NextCont%Cat  == NEWCAT .AND. &
                    tmpCnt%NextCont%Hier >  NEWHIR ) THEN
                  EXIT
               ENDIF
            ENDIF

            ! Advance in list if none of the above checks was true.
            tmpCnt => tmpCnt%NextCont
         ENDDO

      ! Scale factors (tracer IDs less or equal than zero) are collected at
      ! the end of the list. Hence, make tmpCnt pointer point to the
      ! last container w/ base emissions (or the last container in
      ! the list).
      ELSE

         ! Loop over list 
         DO WHILE ( ASSOCIATED ( tmpCnt%NextCont ) )

            ! Check if next container is scale factor 
            IF ( tmpCnt%NextCont%DataType > 1 ) EXIT

            ! Advance in list
            tmpCnt => tmpCnt%NextCont
         ENDDO

      ENDIF

      ! Add this container AFTER current one
      NewCnt%NextCont => tmpCnt%NextCont
      tmpCnt%NextCont => NewCnt 

      ! Cleanup
      tmpCnt => NULL()

      END SUBROUTINE Add2EmisLL
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: InitEmisCont
!
! !DESCRIPTION: Subroutine InitEmisCont initializes an emission
! container. 
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE InitEmisCont( ThisECont )
!
! !ARGUMENTS:
!
      TYPE(EmisCont), POINTER    :: ThisECont
!
! !REVISION HISTORY:
!  04 Dec 2012 - C. Keller: Initialization
!
!EOP
!------------------------------------------------------------------------------
!BOC

      !======================================================================
      ! InitEmisCont begins here!
      !======================================================================

      ! Allocate the new container
      ALLOCATE( ThisECont )

      ! Nullify all pointers
      NULLIFY ( ThisECont%Array     )
      NULLIFY ( ThisECont%cID       )
      NULLIFY ( ThisECont%scalID    )
      NULLIFY ( ThisECont%DataType  )
      NULLIFY ( ThisECont%ExtNr     )
      NULLIFY ( ThisECont%isScalar  )
      NULLIFY ( ThisECont%TempRes   )
      NULLIFY ( ThisECont%TimeSlc   )
      NULLIFY ( ThisECont%cName     )
      NULLIFY ( ThisECont%TrcID     )
      NULLIFY ( ThisECont%Cat       )
      NULLIFY ( ThisECont%Hier      )
      NULLIFY ( ThisECont%OPER      )
      NULLIFY ( ThisECont%NEXTCONT  )
      NULLIFY ( ThisECont%Scal_cID  )

      END SUBROUTINE InitEmisCont
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: CleanupEmisCont
!
! !DESCRIPTION: Subroutine CleanupEmisCont removes an emission
! container.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE CleanupEmisCont ( ThisECont )
!
! !ARGUMENTS:
!
      TYPE(EmisCont), POINTER      :: ThisECont
!
! !REVISION HISTORY:
!  04 Dec 2012 - C. Keller: Initialization
!
!EOP
!------------------------------------------------------------------------------
!BOC

      !======================================================================
      ! CleanupEmisCont begins here!
      !======================================================================

      ! Remove the container
      IF ( ASSOCIATED( ThisECont ) ) THEN

         ! First nullify pointer to the next container 
         ThisECont%NEXTCONT => NULL()

         ! Deallocatable variables. These pointers can be
         ! entirely removed as they do not point to data
         ! that is the target of multiple pointers

         IF ( ASSOCIATED ( ThisECont%TempRes ) ) THEN
            DEALLOCATE ( ThisECont%TempRes  )
         ENDIF

         IF ( ASSOCIATED ( ThisECont%Array ) ) THEN
            DEALLOCATE ( ThisECont%Array      )
         ENDIF

         IF ( ASSOCIATED ( ThisECont%cID ) ) THEN
            DEALLOCATE ( ThisECont%cID    )
         ENDIF

         IF ( ASSOCIATED ( ThisECont%isScalar ) ) THEN
            DEALLOCATE ( ThisECont%isScalar)
         ENDIF

         IF ( ASSOCIATED ( ThisECont%Scal_cID ) ) THEN
            DEALLOCATE ( ThisECont%Scal_cID )
         ENDIF

         ! Nullify those pointers which point to content used
         ! by some other pointers too. We must not deallocate
         ! this content as this would leave some other pointers
         ! unspecified.
         ThisECont%TimeSlc  => NULL()
         ThisECont%cName    => NULL()
         ThisECont%ScalID   => NULL()
         ThisECont%TrcID    => NULL()
         ThisECont%Cat      => NULL()
         ThisECont%Hier     => NULL()
         ThisECont%DataType => NULL()
         ThisECont%ExtNr    => NULL()
         ThisECont%Oper     => NULL()

         ! Deallocate the emission container
         IF ( ASSOCIATED ( ThisECont ) ) DEALLOCATE ( ThisECont )

      ENDIF

      END SUBROUTINE CleanupEmisCont
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: Cleanup_EmisLL
!
! !DESCRIPTION: Subroutine Cleanup\_EmisLL clears the emissions linked list. 
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE Cleanup_EmisLL ( EmisLL )
!
! !ARGUMENTS
!
      TYPE(EmisCont), POINTER  :: EmisLL
!
! !REVISION HISTORY:
!  04 Dec 2012 - C. Keller: Initialization
!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL ARGUMENTS:
!
      TYPE(EmisCont), POINTER   :: CurrCnt => NULL()
      TYPE(EmisCont), POINTER   :: NextCnt => NULL()
      INTEGER                   :: II
     
      !======================================================================
      ! Cleanup_EmisLL begins here!
      !======================================================================

      ! First remove the list ID 
      IF ( ASSOCIATED ( IDList ) ) THEN
         DO II = 1, nnEmisCont
            NULLIFY ( IDList(II)%PNT )
         ENDDO
         DEALLOCATE ( IDList )
      ENDIF

      ! Set CurrCnt to first entry of the EMISSIONS list
      CurrCnt => EmisLL
 
      ! Loop over EMISSIONS linked list
      DO WHILE ( ASSOCIATED ( CurrCnt ) ) 

         ! Set next pointer
         NextCnt => CurrCnt%NEXTCONT

         ! Remove data
         CALL CleanupEmisCont ( CurrCnt )

         ! Update pointer
         CurrCnt => NextCnt

      ENDDO

      ! Reset field counter
      nnEmisCont = 0

      END SUBROUTINE Cleanup_EmisLL 
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: CREATE_IDList
!
! !DESCRIPTION: Subroutine CREATE_IDList creates a vector of pointers 
! (IDList) pointing to all available emission fields.
! The numbering of the emission fields corresponds to the field FIDs,
! i.e. the IDList pointer at position 3 will point to field with FID #3.
! This allows a quick access to all fields.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE CREATE_IDList ( am_I_Root, EmisLL, HcoState, RC )
!
! !USES:
!
      USE HCO_TYPE_MOD, ONLY : HCO_State 
!
! !ARGUMENTS:
!
      LOGICAL,         INTENT(IN)    :: am_I_Root
      TYPE(EmisCont),  POINTER       :: EmisLL
      TYPE(HCO_State), POINTER       :: HcoState
      INTEGER,         INTENT(INOUT) :: RC
!
! !REVISION HISTORY:
!  24 Aug 2012 - C. Keller - Initial Version
!EOP
!------------------------------------------------------------------------------
!BOC
!
!
! !LOCAL VARIABLES:
!
      INTEGER                   :: II
      TYPE(EMISCONT), POINTER   :: CurrCnt => NULL()
      LOGICAL                   :: verbose

      !======================================================================
      ! CREATE_IDList begins here
      !======================================================================

      ! Set verbose flag
      verbose = HcoState%verbose

      ! Prompt some output
      IF ( am_I_Root .and. verbose ) THEN
         WRITE(6,*) '=================================================='
         WRITE(6,*) 'Emission fields ID list is being created...'
      ENDIF

      ! Eventually cleanup the list
      IF ( ASSOCIATED ( IDList ) ) THEN
         DO II = 1, NNIDList
            NULLIFY ( IDList(II)%PNT )
         ENDDO
         DEALLOCATE ( IDList )
      ENDIF

      ! Leave if no emission fields defined 
      IF ( nnEmisCont == 0 ) THEN
         IF ( am_I_Root .and. verbose ) THEN
            WRITE(6,*) 'no emission fields defined!'
         ENDIF
         RC = HCO_SUCCESS
         RETURN
      ENDIF

      ! verbose 
      IF ( am_I_Root .and. verbose ) THEN
         WRITE(6,*) 'Number of available fields: ', nnEmisCont
      ENDIF

      ! Allocate IDList
      ALLOCATE ( IDList(nnEmisCont) )

      ! Nullify all pointers. The pointers will be
      ! connected to the emission fields / ids in the next step.
      DO II = 1, nnEmisCont
         NULLIFY( IDList(II)%PNT )
      ENDDO !II

      ! Now set the quicklist pointers 
      IILOOP: DO II = 1, nnEmisCont

         ! testing only:
         IF ( am_I_Root .and. verbose ) THEN
            WRITE(6,*) 'Trying to create IDList for field ', II
         ENDIF

         ! Set CurrCnt to head of emission fields linked list
         CurrCnt => EmisLL

         DO WHILE ( ASSOCIATED ( CurrCnt ) ) 

            ! Check if current field is the one with the correct FID
            IF ( CurrCnt%cID == II ) THEN

               ! Set pointer to emission field
               IDList(II)%PNT => CurrCnt
 
               ! testing only:
               IF ( am_I_Root .and. verbose ) THEN
                  WRITE(6,*) 'IDList entry created for field ', II
               ENDIF

               ! Advance in loop
               CYCLE IILOOP
            ENDIF

            ! If current field is not the right one, advance to next field
            CurrCnt => CurrCnt%NEXTCONT

         ENDDO

      ENDDO IILOOP

      ! Prompt some output
      IF ( am_I_Root .and. verbose ) THEN
         WRITE(6,*) 'IDList created!'
         WRITE(6,*) '=================================================='
      ENDIF

      ! Set NNIDList to total number of fields
      nnIDList = nnEmisCont

      ! Cleanup and leave
      CurrCnt => NULL()
      RC = HCO_SUCCESS

      END SUBROUTINE CREATE_IDList
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: Add_SclID
!
! !DESCRIPTION: Subroutine Add\_SclID assigns the given scale factor field to 
!  the list of scale factors of the given base field.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE Add_SclID( base_cID, scalID, EmisLL, RC )
!
! !USES:
!
!
! !ARGUMENTS
!
      INTEGER,          INTENT(IN)    :: base_cID
      INTEGER,          INTENT(IN)    :: scalID
      TYPE(EmisCont),   POINTER       :: EmisLL 
      INTEGER,          INTENT(INOUT) :: RC
!
! !REVISION HISTORY:
!  06 Dec 2012 - C. Keller - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      TYPE(EmisCont), POINTER          :: BaseCnt  => NULL()
      TYPE(EmisCont), POINTER          :: ScalCnt  => NULL()
      LOGICAL                          :: FOUND
      INTEGER                          :: scal_cID, II 
      CHARACTER(LEN=255)               :: MSG, LOC

      !======================================================================
      ! Add_SclID begins here!
      !======================================================================

      ! Enter
      LOC = 'Add_SclID'

      ! Point to base field by looking up the container ID
      CALL Find_EmisCont ( base_cID, EmisLL, 0, FOUND, BaseCnt )
      IF ( .NOT. FOUND ) THEN
         MSG = 'Base container not found: ' // TRIM(BaseCnt%cName)
         CALL HCO_ERROR ( MSG, LOC, RC); RETURN
      ENDIF

      ! Point to scale field by lookin up the original ID
      CALL Find_EmisCont ( scalID, EmisLL, 1, FOUND, ScalCnt )
      IF ( .NOT. FOUND ) THEN
         MSG = 'Scale container not found: ' // TRIM(ScalCnt%cName)
         CALL HCO_ERROR ( MSG, LOC, RC); RETURN
      ENDIF

      ! Get container ID of the desired scale factor
      scal_cID = ScalCnt%cID

      ! Loop over the vector with the scale factor field IDs to find a
      ! free slot
      DO II = 1, SclMax

         ! Check if this is an empty slot. In this case, assign ID to
         ! this slot and exit loop
         IF ( BaseCnt%Scal_cID(II) < 0 ) THEN
            BaseCnt%Scal_cID(II) = scal_cID
            EXIT

         ! Check if this field ID has already been assigned, in which
         ! case we leave!
         ELSEIF ( BaseCnt%Scal_cID(II) == scal_cID ) THEN
            EXIT
         ENDIF
           
         ! Error if all slots are already filled
         IF ( II == SclMax ) THEN
            MSG = 'Maximum number of scale factors reached!'
            CALL HCO_ERROR ( MSG, LOC, RC); RETURN
         ENDIF

      ENDDO !II

      ! Cleanup 
      BaseCnt => NULL()
      ScalCnt => NULL()

      ! Return w/ success
      RC = HCO_SUCCESS

      END SUBROUTINE Add_SclID
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: Show_EmisLL
!
! !DESCRIPTION: Subroutine Show\_EmisLL prompts the names of all defined
! emission containers.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE Show_EmisLL ( am_I_Root, EmisLL )
!
! ARGUMENTS
!
      LOGICAL,        INTENT(IN)  :: am_I_Root
      TYPE(EmisCont), POINTER     :: EmisLL
!
! !REVISION HISTORY:
!  06 Dec 2012 - C. Keller - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      TYPE(EmisCont), POINTER   :: CurrCnt => NULL()

      !======================================================================
      ! Show_EmisLL begins here!
      !======================================================================

      ! Leave here if not root
      IF ( .NOT. am_I_Root ) RETURN

      ! Check if an emission family has already been defined or not
      IF ( nnEmisCont == 0 ) THEN
         WRITE(6,*) 'No emission fields defined yet'
         RETURN
      ENDIF

      ! Set CurrCnt to first field
      CurrCnt => EmisLL

      ! Start writing
      WRITE(6,*) 'Currently defined emission fields:'
      DO WHILE ( ASSOCIATED( CurrCnt ) )

         ! Write field name, ID, category and hierarchy
         WRITE ( 6, 100 ) TRIM(CurrCnt%cName)
         WRITE ( 6, 110 ) CurrCnt%DataType, CurrCnt%ExtNr
         WRITE ( 6, 120 ) CurrCnt%cID,      CurrCnt%TrcID 
         WRITE ( 6, 130 ) CurrCnt%Cat,      CurrCnt%Hier
         
         ! Advance to next field
         CurrCnt => CurrCnt%NEXTCONT

      ENDDO

 100  FORMAT( ' -> ', a, ':' )
 110  FORMAT( '    DataType: ',i5,'; ExtNr : ',i5)
 120  FORMAT( '    cID     : ',i5,'; Tracer: ',i5)
 130  FORMAT( '    Cat     : ',i5,'; Hier  : ',i5) 

      ! Cleanup
      CurrCnt => NULL()

      END SUBROUTINE Show_EmisLL
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: Get_nnEmisCont
!
! !DESCRIPTION: Function Get\_nnEmisCont returns the number of containers
! currently defined in the linked list.
!\\
!\\
! !INTERFACE:
!
      FUNCTION Get_nnEmisCont RESULT ( NN )
!
! !RETURN VALUE:
!
      INTEGER :: NN
!
! !REVISION HISTORY:
!  11 Apr 2012 - C. Keller - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
      NN = nnEmisCont

      END FUNCTION Get_nnEmisCont
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: GET_nnIDList
!
! !DESCRIPTION: Function GET\_nnIDList returns the number of pointers
! currently defined in the IDList 
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_nnIDList RESULT ( NN )
!
! !RETURN VALUE:
!
      INTEGER :: NN
!
! !REVISION HISTORY:
!  11 Apr 2012 - C. Keller - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
      NN = nnIDList

      END FUNCTION GET_nnIDList
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: Pnt2EmisCont
!
! !DESCRIPTION: Subroutine Pnt2EmisCont returns a pointer pointing to
! the linked list entry with the given field ID.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE Pnt2EmisCont ( ID, ThisCont, RC ) 
!
! !ARGUMENTS:
!
      INTEGER,        INTENT(IN)     :: ID
      TYPE(EmisCont), POINTER        :: ThisCont
      INTEGER,        INTENT(INOUT)  :: RC
!
! !REVISION HISTORY:
!  11 Apr 2012 - C. Keller - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      CHARACTER(LEN=255)  :: MSG, LOC

      !======================================================================
      ! Pnt2EmisCont begins here!
      !======================================================================

      ! Enter
      LOC = 'Pnt2EmisCont'

      ! Check input 
      IF ( ID > nnIDList ) THEN
         MSG = 'ID higher than number of fields in IDList'
         CALL HCO_ERROR ( MSG, LOC, RC); RETURN
      ENDIF

      ! Point to field with the given field ID
      ThisCont => IDList(ID)%PNT
 
      ! Leave
      RC = HCO_SUCCESS

      END SUBROUTINE Pnt2EmisCont
!EOC
      END MODULE HCO_EMISLL_MOD
!EOF
