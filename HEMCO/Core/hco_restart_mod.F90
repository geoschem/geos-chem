!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: hco_restart_mod.F90
!
! !DESCRIPTION: Module HCO\_RESTART\_MOD contains wrapper routines to define,
! get and write restart fields.  
!\\
!\\
! Restart variables are required by some of the HEMCO extensions. The
! HEMCO restart variables are organized through the HEMCO diagnostics module. 
! . At the end of a simulation, a HEMCO restart file is written that contains 
! all diagnostics that have not been written out yet. In particular, all 
! diagnostic fields with assigned output frequency 'End' are added to the 
! HEMCO restart file. The fields in the HEMCO restart file can then be used for 
! a 'warm' restart by referencing to them in the HEMCO configuration file. 
! Within a module, the fields can be fetched on the first run time step via
! HCO\_GetPtr.
!\\
!\\
! In an ESMF/MAPL environment, it is also possible to store the restart 
! fields in the ESMF internal state object. In this case, the fields must
! be copied from the internal state object on the first run time step, and be
! passed back to the internal state upon finalization.
!\\
!\\
! This module contains wrapper routines to define, obtain and write restart
! fields for the two aforementioned restart types. In an ESMF environment, it
! is always checked first if a given variable exists in the internal state, in
! which case this field is obtained/written. If no internal state object exist,
! an attempt is made to obtain the field through the HEMCO data list, i.e. it
! is checked if the restart variable is specified in the HEMCO configuration
! file.
!\\
!\\
! !INTERFACE: 
!      
MODULE HCO_RESTART_MOD
!
! !USES:
!      
  USE HCO_ERROR_MOD

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!      
  ! defined in all environment
  PUBLIC :: HCO_RestartDefine
  PUBLIC :: HCO_RestartGet
  PUBLIC :: HCO_RestartWrite
!
! !PRIVATE MEMBER FUNCTIONS:
!      
#if defined(ESMF_)
  PRIVATE :: HCO_CopyFromIntnal_ESMF
#endif
!
! !REVISION HISTORY:
!  10 Mar 2015 - C. Keller   - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: HCO_RestartDefine
!
! !DESCRIPTION: Subroutine HCO\_RestartDefine defines a restart diagnostics.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCO_RestartDefine( am_I_Root, HcoState, Name, Arr2D, &
                                Unit,      RC                      )
!
! !USES:
!
    USE HCO_DIAGN_MOD,    ONLY : Diagn_Create
    USE HCO_STATE_MOD,    ONLY : HCO_State
!
! !INPUT ARGUMENTS:
!
    LOGICAL,             INTENT(IN   )         :: am_I_Root
    TYPE(HCO_State),     POINTER               :: HcoState
    CHARACTER(LEN=*),    INTENT(IN   )         :: Name 
    REAL(sp),            INTENT(IN   ), TARGET :: Arr2D(HcoState%NX,HcoState%NY)
    CHARACTER(LEN=*),    INTENT(IN   )         :: Unit 
!
! !INPUT/OUTPUT ARGUMENTS:
!
    INTEGER,             INTENT(INOUT)         :: RC
!
! !REVISION HISTORY:
!  29 Aug 2013 - C. Keller - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! ================================================================
    ! HCO_RestartDefine begins here
    ! ================================================================

    ! Define diagnostics array
    CALL Diagn_Create ( am_I_Root,               & 
                        HcoState   = HcoState,   &
                        cName      = TRIM(Name), &
                        ExtNr      = -1,         &
                        Cat        = -1,         &
                        Hier       = -1,         &
                        HcoID      = -1,         &
                        SpaceDim   =  2,         &
                        OutUnit    = TRIM(Unit), &
                        WriteFreq  = 'End',      &
                        AutoFill   = 0,          &
                        Trgt2D     = Arr2D,      &
                        RC         = RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Return w/ success
    RC = HCO_SUCCESS

  END SUBROUTINE HCO_RestartDefine
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: HCO_RestartGet
!
! !DESCRIPTION: Subroutine HCO\_RestartGet attempts to read a restart field.
! In an ESMF environment, it first checks if the given field (name) is included
! in the internal state object, in which case the data object is filled with
! these values. If not found or if not in an ESMF environment, the HEMCO data
! list (specified in the HEMCO configuration file) is searched. A default value
! can be specified in case that no field could be imported via ESMF and/or the
! HEMCO interface. 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCO_RestartGet( am_I_Root, HcoState, Name,   Arr2D, &
                             RC,        FOUND,    Def2D,  DefVal  )
!
! !USES:
!
    USE HCO_STATE_MOD,    ONLY : HCO_State
    USE HCO_EMISLIST_MOD, ONLY : HCO_GetPtr
!
! !INPUT ARGUMENTS:
!
    LOGICAL,             INTENT(IN   )           :: am_I_Root
    TYPE(HCO_State),     POINTER                 :: HcoState
    CHARACTER(LEN=*),    INTENT(IN   )           :: Name 
    REAL(sp),            INTENT(IN   ), OPTIONAL :: Def2D(HcoState%NX,HcoState%NY)
    REAL(sp),            INTENT(IN   ), OPTIONAL :: DefVal
!
! !OUTPUT ARGUMENTS:
!
    LOGICAL,             INTENT(  OUT), OPTIONAL :: FOUND
!
! !INPUT/OUTPUT ARGUMENTS:
!
    REAL(sp),            INTENT(INOUT)           :: Arr2D(HcoState%NX,HcoState%NY)
    INTEGER,             INTENT(INOUT)           :: RC
!
! !REVISION HISTORY:
!  29 Aug 2013 - C. Keller - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    REAL(sp), POINTER    :: Ptr2D(:,:) => NULL()
    LOGICAL              :: FILLED
    CHARACTER(LEN=255)   :: MSG

    ! ================================================================
    ! HCO_RestartGet begins here
    ! ================================================================

    ! Is the output array filled yet?
    FILLED = .FALSE.

    ! ------------------------------------------------------------------
    ! Try to get from ESMF internal state 
    ! ------------------------------------------------------------------
#if defined(ESMF_)
    CALL HCO_CopyFromIntnal_ESMF( am_I_Root, HcoState, TRIM(Name),   &
                                  Arr2D,     1,        FILLED,     RC )
    IF ( RC /= HCO_SUCCESS ) RETURN
 
    ! Log output 
    IF ( am_I_Root .AND. FILLED ) THEN
       MSG = 'Obtained restart variable from ESMF internal state: '//TRIM(Name)
       CALL HCO_MSG(MSG)
    ENDIF
#endif

    ! ------------------------------------------------------------------
    ! If not yet filled, try to get from HEMCO configuration
    ! ------------------------------------------------------------------
    IF ( .NOT. FILLED ) THEN

       ! Try to get pointer from HEMCO configuration
       CALL HCO_GetPtr( am_I_Root, TRIM(Name), Ptr2D, RC, FOUND=FILLED )
       IF ( RC /= HCO_SUCCESS ) RETURN
      
       ! Eventually pass data
       IF ( FILLED ) THEN
          Arr2D = Ptr2D

          ! Log output 
          IF ( am_I_Root ) THEN
             MSG = 'Obtained restart variable from HEMCO config: '//TRIM(Name)
             CALL HCO_MSG(MSG)
          ENDIF
       ENDIF

       ! Cleanup
       Ptr2D => NULL()
    ENDIF

    ! ------------------------------------------------------------------
    ! If still not filled, assign default values 
    ! ------------------------------------------------------------------
    IF ( .NOT. FILLED ) THEN
       IF ( PRESENT(Def2D) ) THEN
          Arr2D = Def2D
          IF ( am_I_Root ) THEN
             MSG = 'Filled restart variable with default 2D field: '//TRIM(Name)
             CALL HCO_MSG(MSG)
          ENDIF
       ELSEIF( PRESENT(DefVal) ) THEN
          Arr2D = DefVal
          IF ( am_I_Root ) THEN
             MSG = 'Filled restart variable with default scalar: '//TRIM(Name)
             CALL HCO_MSG(MSG)
          ENDIF
       ENDIF
    ENDIF

    ! ------------------------------------------------------------------
    ! Leave 
    ! ------------------------------------------------------------------
    IF ( PRESENT(FOUND) ) THEN
       FOUND = FILLED
    ENDIF

    ! Return w/ success
    RC = HCO_SUCCESS

  END SUBROUTINE HCO_RestartGet
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: HCO_RestartWrite
!
! !DESCRIPTION: Subroutine HCO\_RestartWrite writes a restart variable to the
! ESMF internal state. This is only of relevance in an ESMF environment. The
! 'regular' HEMCO diagnostics created in HCO\_RestartDefine becomes 
! automatically written to disk. 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCO_RestartWrite( am_I_Root, HcoState, Name, Arr2D, RC, FOUND ) 
!
! !USES:
!
    USE HCO_STATE_MOD,    ONLY : HCO_State
!
! !INPUT ARGUMENTS:
!
    LOGICAL,             INTENT(IN   )           :: am_I_Root
    TYPE(HCO_State),     POINTER                 :: HcoState
    CHARACTER(LEN=*),    INTENT(IN   )           :: Name 
!
! !OUTPUT ARGUMENTS:
!
    LOGICAL,             INTENT(  OUT), OPTIONAL :: FOUND
!
! !INPUT/OUTPUT ARGUMENTS:
!
    REAL(sp),            INTENT(INOUT)           :: Arr2D(HcoState%NX,HcoState%NY)
    INTEGER,             INTENT(INOUT)           :: RC
!
! !REVISION HISTORY:
!  29 Aug 2013 - C. Keller - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    LOGICAL   :: WRITTEN

    ! ================================================================
    ! HCO_RestartWrite begins here
    ! ================================================================

    ! Data written to internal state?
    WRITTEN = .FALSE.
    
#if defined(ESMF_)
    CALL HCO_CopyFromIntnal_ESMF( am_I_Root, HcoState, TRIM(Name), &
                                  Arr2D,     -1,       WRITTEN,    RC ) 
#endif

    ! Pass to output
    IF ( PRESENT(FOUND) ) THEN
       FOUND = WRITTEN
    ENDIF

    ! Return w/ success
    RC = HCO_SUCCESS

  END SUBROUTINE HCO_RestartWrite
!EOC
#if defined(ESMF_)
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: HCO_CopyFromIntnal_ESMF
!
! !DESCRIPTION: Subroutine HCO\_CopyFromIntnal\_ESMF attempts to transfer 
! data to and from the ESMF/MAPL internal state. 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCO_CopyFromIntnal_ESMF ( am_I_Root, HcoState,  Name,    &
                                    Arr2D,     Direction, Found, RC )
!
! !USES:
!
#include "MAPL_Generic.h"
    USE ESMF
    USE MAPL_Mod
    USE HCO_STATE_MOD,   ONLY : Hco_State
!
! !ARGUMENTS:
!
    LOGICAL,             INTENT(IN   )   :: am_I_Root
    TYPE(HCO_State),     POINTER         :: HcoState
    CHARACTER(LEN=*),    INTENT(IN   )   :: Name
    REAL(sp),            INTENT(INOUT)   :: Arr2D(HcoState%NX,HcoState%NY)
    INTEGER,             INTENT(IN   )   :: Direction    ! 1: internal to Arr2D; -1: Arr2D to internal
    LOGICAL,             INTENT(  OUT)   :: Found
    INTEGER,             INTENT(INOUT)   :: RC
!
! !REVISION HISTORY:
!  10 Mar 2015 - C. Keller - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER                      :: STAT
    TYPE(MAPL_MetaComp), POINTER :: STATE
    TYPE(ESMF_STATE)             :: INTERNAL
    REAL,                POINTER :: Ptr2D(:,:)   => NULL()

    ! ================================================================
    ! HCO_CopyFromIntnal_ESMF begins here
    ! ================================================================

    ! For MAPL/ESMF error handling (defines Iam and STATUS)
    __Iam__('HCO_CopyFromIntnal_ESMF (HCOI_ESMF_MOD.F90)')

    ! Get internal state
    CALL MAPL_GetObjectFromGC( HcoState%GridComp, STATE, __RC__ )
    CALL MAPL_Get ( STATE, INTERNAL_ESMF_STATE=INTERNAL, __RC__ )

    ! Try to import field
    CALL MAPL_GetPointer( INTERNAL, Ptr2D, TRIM(Name), &
                          NotFoundOk=.TRUE., __RC__ ) 
   
    ! Eventually copy data to or from output array 
    IF ( ASSOCIATED(Ptr2D) ) THEN
    
       ! Make sure we can copy the data
       ASSERT_(SIZE(Arr2D,1)==SIZE(Ptr2D,1))   
       ASSERT_(SIZE(Arr2D,2)==SIZE(Ptr2D,2))
    
       ! transfer direction must be 1 or -1
       ASSERT_(Direction==1 .OR. Direction==-1)
 
       ! Transfer data
       IF ( Direction == 1 ) THEN
          Arr2D = Ptr2D
       ELSEIF ( Direction == -1 ) THEN
          Ptr2D = Arr2D
       ENDIF
       Found = .TRUE.
    ELSE
       Found = .FALSE.
    ENDIF

    ! Cleanup
    Ptr2D => NULL()

    ! Return success
    RC = HCO_SUCCESS 

  END SUBROUTINE HCO_CopyFromIntnal_ESMF 
!EOC
#endif
END MODULE HCO_RESTART_MOD
