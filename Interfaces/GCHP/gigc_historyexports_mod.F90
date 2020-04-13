!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: gigc_historyexports_mod.F90
!
! !DESCRIPTION: Module GIGC\_HistoryExports\_Mod serves as the interface
! between the HISTORY configuration file, the GEOS-Chem State registry,
! and the ESMF Export State.
!\\
!\\
! !INTERFACE: 
!
MODULE GIGC_HistoryExports_Mod
!
! !USES:
!
#include "MAPL_Generic.h"
  USE DiagList_Mod
  USE ErrCode_Mod
  USE Precision_Mod
  USE MAPL_Mod

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC :: HistoryExports_SetServices
  PUBLIC :: HistoryExports_SetDataPointers
  PUBLIC :: CopyGCStates2Exports
  PUBLIC :: Destroy_HistoryConfig 
!
! !PRIVATE:
!
  PRIVATE :: Init_HistoryConfig
  PRIVATE :: Init_HistoryExport
  PRIVATE :: Init_HistoryExportsList
  PRIVATE :: Append_HistoryExportsList
  PRIVATE :: Check_HistoryExportsList
  PRIVATE :: Print_HistoryExportsList
!
! !PUBLIC TYPES
!
  ! History Configuration Object 
  TYPE, PUBLIC :: HistoryConfigObj

     CHARACTER(LEN=255)                   :: ROOT ! TODO: needed?
     CHARACTER(LEN=255)                   :: ConfigFileName
     LOGICAL                              :: ConfigFileRead
     TYPE(HistoryExportsListObj), POINTER :: HistoryExportsList
     TYPE(DgnList)                        :: DiagList
 
 END TYPE HistoryConfigObj
!
! !PRIVATE TYPES
!
  ! History Exports Linked List
  TYPE :: HistoryExportsListObj

     TYPE(HistoryExportObj), POINTER :: head
     INTEGER                         :: numExports

  END TYPE HistoryExportsListObj

  ! History Export Object
  TYPE :: HistoryExportObj

     CHARACTER(LEN=255)              :: name 
     CHARACTER(LEN=255)              :: metadataID
     CHARACTER(LEN=255)              :: registryID
     CHARACTER(LEN=255)              :: long_name  
     CHARACTER(LEN=255)              :: units       
     INTEGER                         :: vloc
     INTEGER                         :: rank        
     INTEGER                         :: type
     LOGICAL                         :: isMet
     LOGICAL                         :: isChem
     LOGICAL                         :: isDiag
     TYPE(HistoryExportObj), POINTER :: next

     ! Pointers to ESMF Export and GEOS-Chem State
     ! TODO: for now, include all possible data types in the registry. 
     REAL,     POINTER :: ExportData2d(:,:)
     REAL,     POINTER :: ExportData3d(:,:,:)
     REAL(fp), POINTER :: GCStateData0d
     REAL(fp), POINTER :: GCStateData1d(:)
     REAL(fp), POINTER :: GCStateData2d(:,:)
     REAL(fp), POINTER :: GCStateData3d(:,:,:)
     REAL(f4), POINTER :: GCStateData0d_4
     REAL(f4), POINTER :: GCStateData1d_4(:)
     REAL(f4), POINTER :: GCStateData2d_4(:,:)
     REAL(f4), POINTER :: GCStateData3d_4(:,:,:) 
     REAL(f8), POINTER :: GCStateData0d_8
     REAL(f8), POINTER :: GCStateData1d_8(:)
     REAL(f8), POINTER :: GCStateData2d_8(:,:)
     REAL(f8), POINTER :: GCStateData3d_8(:,:,:)
     INTEGER,  POINTER :: GCStateData0d_I
     INTEGER,  POINTER :: GCStateData1d_I(:)
     INTEGER,  POINTER :: GCStateData2d_I(:,:)
     INTEGER,  POINTER :: GCStateData3d_I(:,:,:)

  END TYPE HistoryExportObj
!
! !PUBLIC PARAMETERS
!
  ! Prefix of the species names in the internal state and HISTORY.rc
  CHARACTER(LEN=4), PUBLIC, PARAMETER  :: SPFX = 'SPC_'

#if defined( MODEL_GEOS )
  ! GEOS-Chem advected species in GEOS use prefix TRC_. Use GCD_ for GEOS-only.
  CHARACTER(LEN=4), PUBLIC, PARAMETER  :: TPFX = 'TRC_'
  CHARACTER(LEN=4), PUBLIC, PARAMETER  :: GPFX = 'GCD_'
#endif
!
! !REVISION HISTORY:
!  01 Sep 2017 - E. Lundgren - Initial version
!  08 Mar 2018 - E. Lundgren - Define the internal state prefix expected in 
!                              HISTORY.rc within this module
!EOP
!------------------------------------------------------------------------------
!BOC

CONTAINS
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Init_HistoryConfig 
!
! !DESCRIPTION:
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Init_HistoryConfig ( am_I_Root, HistoryConfig, configFile, RC )
!
! !USES:
!
!
! !INPUT PARAMETERS:
!
    LOGICAL,             INTENT(IN) :: am_I_Root
    CHARACTER(LEN=*),    INTENT(IN) :: configFile
!
! !OUTPUT PARAMETERS:
!
    TYPE(HistoryConfigObj), POINTER :: HistoryConfig
    INTEGER, INTENT(OUT)            :: RC 
!
! !REVISION HISTORY:
!  01 Sep 2017 - E. Lundgren - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
    __Iam__('Init_HistoryConfig (gigc_historyexports_mod.F90)')
    ALLOCATE(HistoryConfig)
    HistoryConfig%ROOT               =  ''
    HistoryConfig%ConfigFileName     =  TRIM(configFile)
    HistoryConfig%ConfigFileRead     =  .FALSE.
    CALL Init_DiagList( am_I_Root, configFile, HistoryConfig%DiagList, RC )
    IF ( RC == GC_FAILURE ) THEN
       _ASSERT(.FALSE., 'informative message here')
       RETURN
    ENDIF
    ! Optional debugging
    !CALL Print_DiagList( am_I_Root, HistoryConfig%DiagList, RC )

    CALL Init_HistoryExportsList( am_I_Root, HistoryConfig, RC )
    IF ( RC == GC_FAILURE ) THEN
       _ASSERT(.FALSE., 'informative message here')
       RETURN
    ENDIF
    ! Optional debugging
    !CALL Print_HistoryExportsList( am_I_Root, HistoryConfig, RC )

  END SUBROUTINE Init_HistoryConfig
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Init_HistoryExportsList 
!
! !DESCRIPTION: 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Init_HistoryExportsList ( am_I_Root, HistoryConfig, RC )
!
! !USES:
!
    USE State_Chm_Mod,    ONLY: Get_Metadata_State_Chm
    USE State_Diag_Mod,   ONLY: Get_Metadata_State_Diag 
    USE State_Met_Mod,    ONLY: Get_Metadata_State_Met
!
! !INPUT PARAMETERS:
!
    LOGICAL,          INTENT(IN)    :: am_I_Root
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(HistoryConfigObj), POINTER :: HistoryConfig
!
! !OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(OUT)   :: RC
!
! !REVISION HISTORY:
!  01 Sep 2017 - E. Lundgren - initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER               :: N, rank, vloc, type
    CHARACTER(LEN=255)    :: ErrMsg, desc, units, tag
    LOGICAL               :: isMet, isChem, isDiag, found
    TYPE(HistoryExportObj),  POINTER :: NewHistExp
    TYPE(DgnItem),           POINTER :: current

    ! ================================================================
    ! Init_HistoryExportsList begins here
    ! ================================================================
    __Iam__('Init_HistoryExportsList (gigc_historyexports_mod.F90)')

    ! Init
    NewHistExp => NULL()

    ! Create HistoryExportsList object
    ALLOCATE(HistoryConfig%HistoryExportsList)
    HistoryConfig%HistoryExportsList%numExports = 0
    HistoryConfig%HistoryExportsList%head => NULL()

    ! Loop over entries in DiagList
    current => HistoryConfig%DiagList%head
    DO WHILE ( ASSOCIATED( current ) )

#if defined ( MODEL_GEOS )
       ! Skip State_Chm%Species entries since in internal state.
       ! Also skip GEOS-only diagnostics.
       IF ( ( INDEX( current%name,  TRIM(TPFX) ) == 1 ) .OR.   &
            ( INDEX( current%name,  TRIM(SPFX) ) == 1 ) .OR.   &
            ( INDEX( current%name,  TRIM(GPFX) ) == 1 ) ) THEN
          current => current%next
          CYCLE
       ENDIF
#else
       ! Skip State_Chm%Species entries since in internal state.
       IF ( INDEX( current%name,  TRIM(SPFX) ) == 1 ) THEN
          current => current%next
          CYCLE
       ENDIF
#endif

       ! Skip emissions diagnostics since handled by HEMCO
       IF ( INDEX( current%name,  'EMIS' ) == 1 .or. &
            INDEX( current%name,  'INV'  ) == 1 ) THEN
          current => current%next
          CYCLE
       ENDIF

       ! Check history exports list to see if already added (unless wildcard)
       ! TODO: consider making the call a function that returns a logical
       IF ( .NOT. current%isWildcard ) THEN
          CALL Check_HistoryExportsList( am_I_Root, current%name,           &
                                         HistoryConfig%HistoryExportsList,  &
                                         found, RC                         )
          IF ( found ) THEN
             current => current%next
             CYCLE
          ENDIF
       ENDIF

       ! Get metadata using metadataID and state
       ! If isTagged, then append to description
       ! If isWildcard, shouldn't get here
       ! The name of the export is simply name
       Found = .TRUE.
       isMet  = .FALSE.
       isChem = .FALSE.
       isDiag = .FALSE.
       IF ( TRIM(current%state) == 'MET' ) THEN
          isMet = .TRUE.
          CALL Get_Metadata_State_Met( am_I_Root, current%metadataID,     &
                                       Found, RC, desc=desc, units=units, &
                                       rank=rank, type=type, vloc=vloc )
          ! TODO: need to add found to outputs of get_metadata_state_met
       ELSEIF ( TRIM(current%state) == 'CHEM' ) THEN
          isCHEM = .TRUE.
          CALL Get_Metadata_State_Chm( am_I_Root, current%metadataID,     &
                                       Found, RC, desc=desc, units=units, &
                                       rank=rank, type=type, vloc=vloc )
       ELSEIF ( TRIM(current%state) == 'DIAG' ) THEN
          isDIAG = .TRUE.
          CALL Get_Metadata_State_Diag( am_I_Root, current%metadataID,     &
                                        Found, Rc, desc=desc, units=units, &
                                        rank=rank, type=type, vloc=vloc )
       ELSEIF ( TRIM(current%state) == 'GEOS5' ) THEN
          ! Skip it
          current => current%next
          CYCLE
       ELSE
          ErrMsg = "Unknown state of item " // TRIM(current%name) // &
                   " in DiagList: " // TRIM(current%state)
          EXIT
       ENDIF
       IF ( Found .eqv. .FALSE. ) THEN
          ErrMsg = "Metadata not found for " // TRIM(current%name)
          EXIT
       ENDIF

       ! If wildcard is present
       IF ( current%isWildcard ) THEN
          ! Do nothing. This should never happen at this point since
          ! Init_DiagList will exit with an error if wildcard is
          ! encountered in HISTORY.rc while compiling with ESMF_.

          ! When it comes time to implement, create exports in a loop,
          ! either for all species or for advected species only. Include 
          ! a check that the export was not already created. Loop over 
          ! AdvNames if wildcard is ADV. Loop over SpecNames for all other 
          ! cases, passing not found = OK so that not all are necessarily 
          ! output. Later on, after species database is initialized, exports 
          ! for only species in the specific wildcard will be associated 
          ! with data and thus included in the output file.

          ! If the meantime, skip wildcards if it gets here.
          current => current%next
          CYCLE
       ENDIF

       ! If this item is for a specific tag, append description. 
       ! This will need revisiting since there may be tag-dependent
       ! strings to append to long names
       IF ( current%isTagged ) THEN
          desc = TRIM(desc) // " for " // TRIM(current%tag)
       ENDIF

       ! Create a new HistoryExportObj object
       CALL Init_HistoryExport( am_I_Root, NewHistExp,         &
                                name=current%name,             &
                                metadataID=current%metadataID, &
                                registryID=current%registryID, &
                                long_name=desc,                &
                                units=units,                   &
                                vloc=vloc,                     &
                                rank=rank,                     &
                                type=type,                     &
                                isMet=isMet,                   &
                                isChem=isChem,                 &
                                isDiag=isDiag,                 &
                                RC=RC )
       IF ( RC == GC_FAILURE ) THEN
          ErrMsg = "History export init fail for " // TRIM(current%name)
          EXIT
       ENDIF
       
       ! Add new HistoryExportObj to linked list
       CALL Append_HistoryExportsList( am_I_Root,     NewHistExp, &
                                       HistoryConfig, RC       )
       IF ( RC == GC_FAILURE ) THEN
          ErrMsg = "History export append fail for " // TRIM(current%name)
          EXIT
       ENDIF

       ! Set up for next item in DiagList
       current => current%next

    ENDDO
    current => NULL()

    IF ( RC == GC_SUCCESS ) THEN
       HistoryConfig%ConfigFileRead = .TRUE.
    ELSE
       CALL GC_ERROR( ErrMsg, RC, Iam )
       _ASSERT(.FALSE., 'informative message here')
       RETURN
    ENDIF

  END SUBROUTINE Init_HistoryExportsList
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Init_HistoryExport 
!
! !DESCRIPTION: 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Init_HistoryExport ( am_I_Root,  NewHistExp, name,         &
                                  metadataID, registryID, long_name,    &
                                  units,      vloc,       rank,         &
                                  type,       isMet,      isChem,       &
                                  isDiag,     RC  )
!
! !INPUT PARAMETERS:
!
    LOGICAL,             INTENT(IN)    :: am_I_Root
!
! !OUTPUT PARAMETERS:
!
    TYPE(HistoryExportObj), POINTER :: NewHistExp
    CHARACTER(LEN=*), OPTIONAL      :: name
    CHARACTER(LEN=*), OPTIONAL      :: metadataID
    CHARACTER(LEN=*), OPTIONAL      :: registryID
    CHARACTER(LEN=*), OPTIONAL      :: long_name
    CHARACTER(LEN=*), OPTIONAL      :: units
    INTEGER,          OPTIONAL      :: vloc 
    INTEGER,          OPTIONAL      :: rank
    INTEGER,          OPTIONAL      :: type
    LOGICAL,          OPTIONAL      :: isMet
    LOGICAL,          OPTIONAL      :: isChem
    LOGICAL,          OPTIONAL      :: isDiag
    INTEGER,          OPTIONAL      :: RC
!
! !REVISION HISTORY:
!  01 Sep 2017 - E. Lundgren - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
    __Iam__('Init_HistoryExport (gigc_historyexports_mod.F90)')
    ALLOCATE(NewHistExp)
    NewHistExp%name        = TRIM(name)
    NewHistExp%metadataID  = TRIM(metadataID)
    NewHistExp%registryID  = TRIM(registryID)
    NewHistExp%long_name   = TRIM(long_name)
    NewHistExp%units       = TRIM(units)
    NewHistExp%vloc        = vloc
    NewHistExp%rank        = rank
    NewHistExp%type        = type
    NewHistExp%isMet       = isMet
    NewHistExp%isChem      = isChem
    NewHistExp%isDiag      = isDiag
    NewHistExp%next          => NULL()
    NewHistExp%ExportData2d  => NULL()
    NewHistExp%ExportData3d  => NULL()
    NewHistExp%GCStateData0d => NULL()
    NewHistExp%GCStateData1d => NULL()
    NewHistExp%GCStateData2d => NULL()
    NewHistExp%GCStateData3d => NULL()
    NewHistExp%GCStateData0d_4 => NULL()
    NewHistExp%GCStateData1d_4 => NULL()
    NewHistExp%GCStateData2d_4 => NULL()
    NewHistExp%GCStateData3d_4 => NULL()
    NewHistExp%GCStateData0d_8 => NULL()
    NewHistExp%GCStateData1d_8 => NULL()
    NewHistExp%GCStateData2d_8 => NULL()
    NewHistExp%GCStateData3d_8 => NULL()
    NewHistExp%GCStateData0d_I => NULL()
    NewHistExp%GCStateData1d_I => NULL()
    NewHistExp%GCStateData2d_I => NULL()
    NewHistExp%GCStateData3d_I => NULL()

  END SUBROUTINE Init_HistoryExport
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Append_HistoryExportsList 
!
! !DESCRIPTION: 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Append_HistoryExportsList ( am_I_Root,     HistoryExport, &
                                         HistoryConfig, RC        )
!
! !USES:
!
!
! !INPUT PARAMETERS:
!
    LOGICAL,          INTENT(IN)    :: am_I_Root
    TYPE(HistoryExportObj), POINTER :: HistoryExport
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(HistoryConfigObj), POINTER :: HistoryConfig
!
! !OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(OUT)   :: RC
!
! !REVISION HISTORY:
!  01 Sep 2017 - E. Lundgren - initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    TYPE(HistoryExportObj),  POINTER :: NewHistExp

    ! ================================================================
    ! Append_HistoryExportsList begins here
    ! ================================================================
    __Iam__('Append_HistoryExportsList (gigc_historyexports_mod.F90)')

    ! Add new object to the beginning of the linked list
    HistoryExport%next => HistoryConfig%HistoryExportsList%head
    HistoryConfig%HistoryExportsList%head => HistoryExport

    ! Update # of list items
    HistoryConfig%HistoryExportsList%numExports = &
         HistoryConfig%HistoryExportsList%numExports + 1

  END SUBROUTINE Append_HistoryExportsList
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Check_HistoryExportsList
!
! !DESCRIPTION:
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Check_HistoryExportsList ( am_I_Root, name,  &
                                        ExportsList, found, RC )
!
! !INPUT PARAMETERS:
!
    LOGICAL,           INTENT(IN)        :: am_I_Root
    CHARACTER(LEN=*),  INTENT(IN)        :: name
    TYPE(HistoryExportsListObj), POINTER :: ExportsList
!
! !OUTPUT PARAMETERS:
!
    LOGICAL, INTENT(OUT)               :: found
    INTEGER, INTENT(OUT)               :: RC 
!
! !REVISION HISTORY:
!  12 Sep 2017 - E. Lundgren - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    TYPE(HistoryExportObj), POINTER :: current

    __Iam__('Check_HistoryExportsList (gigc_historyexports_mod.F90)')

    ! Assume not found
    found = .False.

    current => ExportsList%head
    DO WHILE ( ASSOCIATED( current ) )
       IF ( current%name == name ) THEN
          found = .TRUE.
          RETURN
       ENDIF
       current => current%next    
    ENDDO
    current => NULL()

  END SUBROUTINE Check_HistoryExportsList
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !!
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HistoryExports_SetServices
!
  ! !DESCRIPTION:
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HistoryExports_SetServices( am_I_Root, config_file, GC,  &
                                         HistoryConfig,  RC )
!
! !USES:
!
    USE ESMF, ONLY : ESMF_GridComp
    USE Registry_Params_Mod
!
! !INPUT PARAMETERS:
!
    LOGICAL,             INTENT(IN)    :: am_I_Root
    CHARACTER(LEN=*),    INTENT(IN)    :: config_file
!
! !INPUT AND OUTPUT PARAMETERS:
!
    TYPE(ESMF_GridComp), INTENT(INOUT) :: GC             ! Gridded Component
!
! !OUTPUT PARAMETERS:
!
    TYPE(HistoryConfigObj), POINTER    :: HistoryConfig  ! History config object
    INTEGER,             INTENT(OUT)   :: RC
!
! !REMARKS:
!  !
! !REVISION HISTORY: 
!  01 Sep 2017 - E. Lundgren - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    CHARACTER(LEN=255)                   :: ErrMsg
    TYPE(HistoryExportObj),      POINTER :: current

    ! ================================================================
    ! HistoryExports_SetServices begins here
    ! ================================================================

    ! For MAPL/ESMF error handling (defines Iam and STATUS)
    __Iam__('HistoryExports_SetServices (gigc_historyexports_mod.F90)')

    ! Create a config object if it does not already exist
    IF ( .NOT. ASSOCIATED(HistoryConfig) ) THEN
       CALL Init_HistoryConfig( am_I_Root, HistoryConfig, config_file, RC )
       IF ( RC == GC_FAILURE ) THEN
          _ASSERT(.FALSE., 'informative message here')
          RETURN
       ENDIF
    ENDIF

    ! Loop over the History Exports list to add one export per item
    IF ( am_I_Root ) THEN
       WRITE(6,*) " "
       WRITE(6,*) "Adding history variables to GCHP Export State:"
    ENDIF
    current => HistoryConfig%HistoryExportsList%head
    DO WHILE ( ASSOCIATED( current ) )
       IF ( am_I_Root ) PRINT *, "adding export: ", TRIM(current%name)       
       ! Create an export for this item
       IF ( current%rank == 3 ) THEN
          IF ( current%vloc == VLocationCenter ) THEN
             CALL MAPL_AddExportSpec(GC,                                     &
                                     SHORT_NAME = TRIM(current%name),        &
                                     LONG_NAME  = TRIM(current%long_name),   &
                                     UNITS      = TRIM(current%units),       &
                                     DIMS       = MAPL_DimsHorzVert,         &
                                     VLOCATION  = MAPL_VLocationCenter,      &
                                     RC         = RC                         )
          IF ( RC == GC_FAILURE ) THEN
             ErrMsg =  "Problem adding 3D export for " // TRIM(current%name)
             EXIT
          ENDIF
         ELSEIF ( current%vloc == VLocationEdge ) THEN
            CALL MAPL_AddExportSpec(GC,                                     &
                                    SHORT_NAME = TRIM(current%name), &
                                    LONG_NAME  = TRIM(current%long_name),   &
                                    UNITS      = TRIM(current%units),       &
                                    DIMS       = MAPL_DimsHorzVert,         &
                                    VLOCATION  = MAPL_VLocationEdge,        &
                                    RC         = STATUS                    )
         ELSE
            IF ( am_I_Root ) THEN
               PRINT *, "Unknown vertical location for ", &
                        TRIM(current%name)
            ENDIF
         ENDIF
       ELSEIF ( current%rank == 2 ) THEN
          CALL MAPL_AddExportSpec(GC,                                     &
                                  SHORT_NAME = TRIM(current%name), &
                                  LONG_NAME  = TRIM(current%long_name),   &
                                  UNITS      = TRIM(current%units),       &
                                  DIMS       = MAPL_DimsHorzOnly,         &
                                  RC         = RC                        )
          IF ( RC == GC_FAILURE ) THEN
             ErrMsg =  "Problem adding 2D export for " // TRIM(current%name)
             EXIT
          ENDIF
       ELSE
          ErrMsg = "Problem adding export for " // TRIM(current%name) // &
                   ". Rank is only implemented for 2 or 3!"
          RC = GC_FAILURE
          EXIT
       ENDIF

       current => current%next    
    ENDDO
    current => NULL()

    IF ( RC == GC_FAILURE ) THEN
       CALL GC_ERROR( ErrMsg, RC, Iam )
       _ASSERT(.FALSE., 'informative message here')
       RETURN
    ENDIF
    
  END SUBROUTINE HistoryExports_SetServices
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: CopyGCStates2Exports
!
! !DESCRIPTION:
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE CopyGCStates2Exports( am_I_Root, Input_Opt, HistoryConfig, RC )
!
! !USES:
!
  USE HCOI_GC_Main_Mod, ONLY : HCOI_GC_WriteDiagn
  USE Input_Opt_Mod,    ONLY : OptInput
!
! !INPUT PARAMETERS:
!
    LOGICAL,        INTENT(IN)    :: am_I_Root
    TYPE(OptInput), INTENT(IN)    :: Input_Opt
!
! !INPUT AND OUTPUT PARAMETERS:
!
    TYPE(HistoryConfigObj), POINTER :: HistoryConfig  ! History config object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,            INTENT(OUT) :: RC
!
! !REMARKS:
!  !
! !REVISION HISTORY: 
!  01 Sep 2017 - E. Lundgren - Initial version
!  02 Nov 2017 - E. Lundgren - Copy HEMCO data to emissions exports
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER                         :: LMAX
    CHARACTER(LEN=255)              :: ErrMsg
    TYPE(HistoryExportObj), POINTER :: current

    ! ================================================================
    ! CopyGCStates2Exports begins here
    ! ================================================================
    __Iam__('CopyGCStates2Exports (gigc_historyexports_mod.F90)')

    ! Loop over the History Exports list
    current => HistoryConfig%HistoryExportsList%head
    DO WHILE ( ASSOCIATED( current ) )

       IF ( current%rank == 2 ) THEN
          IF ( ASSOCIATED ( current%GCStateData2d ) ) THEN
             current%ExportData2d = current%GCStateData2d
          ELSEIF ( ASSOCIATED ( current%GCStateData2d_4 ) ) THEN
             current%ExportData2d = current%GCStateData2d_4
          ELSEIF ( ASSOCIATED ( current%GCStateData2d_8 ) ) THEN
             current%ExportData2d = current%GCStateData2d_8
          ELSEIF ( ASSOCIATED ( current%GCStateData2d_I ) ) THEN
             ! Convert integer to float (integers not allowed in MAPL exports)
             current%ExportData2d = FLOAT(current%GCStateData2d_I)
          ELSE
             RC = GC_FAILURE
             ErrMsg = "No GC 2D pointer found for " // TRIM(current%name)
             EXIT
          ENDIF
       ELSEIF ( current%rank == 3 ) THEN
          IF ( ASSOCIATED ( current%GCStateData3d ) ) THEN
             current%ExportData3d = current%GCStateData3d
          ELSEIF ( ASSOCIATED ( current%GCStateData3d_4 ) ) THEN
             current%ExportData3d = current%GCStateData3d_4
          ELSEIF ( ASSOCIATED ( current%GCStateData3d_8 ) ) THEN
             current%ExportData3d = current%GCStateData3d_8
          ELSEIF ( ASSOCIATED ( current%GCStateData3d_I ) ) THEN
             current%ExportData3d = FLOAT(current%GCStateData3d_I)
          ELSE
             RC = GC_FAILURE
             ErrMsg = "No GC 3D pointer found for " // TRIM(current%name)
             EXIT
          ENDIF
#if defined( MODEL_GEOS )
          ! If using GEOS-5, flip the data vertically to match model
          ! convention
          LMAX = SIZE(current%ExportData3d, 3)
          current%ExportData3d(:,:,1:LMAX) =  &
                              current%ExportData3d(:,:,LMAX:1:-1)
#endif
       ENDIF

       current => current%next    
    ENDDO
    current => NULL()

    ! Error handling
    IF ( RC == GC_FAILURE ) THEN
       CALL GC_ERROR( ErrMsg, RC, Iam )
       _ASSERT(.FALSE., 'informative message here')
       RETURN
    ENDIF

    ! Copy emissions data to MAPL exports via HEMCO
    CALL HCOI_GC_WriteDiagn( Input_Opt, .FALSE., RC )
    IF ( RC == GC_FAILURE ) THEN
       ErrMsg = "Error copying emissions data to MAPL via HEMCO"
       CALL GC_ERROR( ErrMsg, RC, Iam )
       _ASSERT(.FALSE., 'informative message here')
       RETURN
    ENDIF
    
  END SUBROUTINE CopyGCStates2Exports
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Print_HistoryExportsList
!
! !DESCRIPTION:
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Print_HistoryExportsList( am_I_Root, HistoryConfig, RC )
!
! !USES:
!
!
! !INPUT PARAMETERS:
!
    LOGICAL,             INTENT(IN) :: am_I_Root
!
! !INPUT AND OUTPUT PARAMETERS:
!
    TYPE(HistoryConfigObj), POINTER :: HistoryConfig  ! History config object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,            INTENT(OUT) :: RC
!
! !REMARKS:
!  !
! !REVISION HISTORY: 
!  01 Sep 2017 - E. Lundgren - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    TYPE(HistoryExportObj), POINTER :: current

    ! ================================================================
    ! Print_HistoryExportsList begins here
    ! ================================================================
    __Iam__('Print_HistoryExportsList (gigc_historyexports_mod.F90)')

    ! Loop over the History Exports list
    current => HistoryConfig%HistoryExportsList%head
    IF ( am_I_Root ) PRINT *, '==========================='
    IF ( am_I_Root ) PRINT *, 'History Exports List:'
    IF ( am_I_Root ) PRINT *, ' '
    DO WHILE ( ASSOCIATED( current ) )
       IF ( am_I_Root ) THEN
          PRINT *, "Name:        ",   TRIM(current%name) 
          PRINT *, " MetadataID: ",   TRIM(current%metadataID)
          PRINT *, " RegistryID: ",   TRIM(current%registryID)
          PRINT *, " Long name:  ",   TRIM(current%long_name)  
          PRINT *, " Units:      ",   TRIM(current%units)       
          PRINT *, " Vert loc:   ",   current%vloc
          PRINT *, " Rank:       ",   current%rank        
          PRINT *, " Type:       ",   current%type
          PRINT *, " isMet:      ",   current%isMet
          PRINT *, " isChem:     ",   current%isChem
          PRINT *, " isDiag:     ",   current%isDiag
          PRINT *, " "
       ENDIF
       current => current%next    
    ENDDO
    IF ( am_I_Root ) PRINT *, '==========================='
    current => NULL()
    
  END SUBROUTINE Print_HistoryExportsList
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HistoryExports_SetDataPointers
!
! !DESCRIPTION:
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HistoryExports_SetDataPointers( am_I_Root,     EXPORT,    &
                                             HistoryConfig, State_Chm, &
                                             State_Diag,    State_Met, &
                                             RC                       )
!
! !USES:
!
    USE ESMF,           ONLY : ESMF_State
    USE Registry_Mod,   ONLY : Registry_Lookup
    USE State_Chm_Mod,  ONLY : ChmState
    USE State_Diag_Mod, ONLY : DgnState
    USE State_Met_Mod,  ONLY : MetState
!
! !INPUT PARAMETERS:
!
    LOGICAL,         INTENT(IN)     :: am_I_Root
!
! !INPUT AND OUTPUT PARAMETERS:
!
    TYPE(ESMF_State), INTENT(INOUT), TARGET :: EXPORT ! Export state obj
    TYPE(HistoryConfigObj), POINTER :: HistoryConfig  ! History config obj
    TYPE(ChmState),   INTENT(INOUT) :: State_Chm      ! Chemistry State obj
    TYPE(MetState),   INTENT(INOUT) :: State_Met      ! Meteorology State obj
    TYPE(DgnState),   INTENT(INOUT) :: State_Diag     ! Diagnostics State obj
!
! !OUTPUT PARAMETERS:
!
    INTEGER,             INTENT(OUT)   :: RC
!
! !REMARKS:
!  !
! !REVISION HISTORY: 
!  01 Sep 2017 - E. Lundgren - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    CHARACTER(LEN=255)              :: ErrMsg
    TYPE(HistoryExportObj), POINTER :: current

    ! ================================================================
    ! HistoryExports_SetDataPointers begins here
    ! ================================================================
    __Iam__('HistoryExports_SetDataPointers')

    IF ( am_I_Root ) THEN
       WRITE(6,*) " "
       WRITE(6,*) "Setting history variable pointers to GC and Export States:"
    ENDIF

    ! Loop over the History Exports list
    current => HistoryConfig%HistoryExportsList%head
    DO WHILE ( ASSOCIATED( current ) )

       ! Get pointer to GC state data
       !IF ( am_I_Root ) WRITE(6,*) current%name
       IF ( current%isMET ) THEN
          CALL Registry_Lookup( am_I_Root = am_I_Root,               &
                                Registry  = State_Met%Registry,      &
                                State     = State_Met%State,         &
                                RegDict   = State_Met%RegDict,       &
                                Variable  = current%registryID,      &
                                Ptr2d     = current%GCStateData2d,   &
                                Ptr2d_4   = current%GCStateData2d_4, &
                                Ptr2d_8   = current%GCStateData2d_8, &
                                Ptr2d_I   = current%GCStateData2d_I, &
                                Ptr3d     = current%GCStateData3d,   &
                                Ptr3d_4   = current%GCStateData3d_4, &
                                Ptr3d_8   = current%GCStateData3d_8, &
                                Ptr3d_I   = current%GCStateData3d_I, &
                                RC        = RC                      )
       ELSEIF ( current%isChem ) THEN
          CALL Registry_Lookup( am_I_Root = am_I_Root,               &
                                Registry  = State_Chm%Registry,      &
                                RegDict   = State_Chm%RegDict,       &
                                State     = State_Chm%State,         &
                                Variable  = current%registryID,      &
                                Ptr2d     = current%GCStateData2d,   &
                                Ptr2d_4   = current%GCStateData2d_4, &
                                Ptr2d_8   = current%GCStateData2d,   &
                                Ptr2d_I   = current%GCStateData2d_I, &
                                Ptr3d     = current%GCStateData3d,   &
                                Ptr3d_4   = current%GCStateData3d_4, &
                                Ptr3d_8   = current%GCStateData3d,   &
                                Ptr3d_I   = current%GCStateData3d_I, &
                                RC        = RC                      )
       ELSEIF ( current%isDiag ) THEN
          CALL Registry_Lookup( am_I_Root = am_I_Root,               &
                                Registry  = State_Diag%Registry,     &
                                RegDict   = State_Diag%RegDict,      &
                                State     = State_Diag%State,        &
                                Variable  = current%registryID,      &
                                Ptr2d     = current%GCStateData2d,   &
                                Ptr2d_4   = current%GCStateData2d_4, &
                                Ptr2d_8   = current%GCStateData2d,   &
                                Ptr2d_I   = current%GCStateData2d_I, &
                                Ptr3d     = current%GCStateData3d,   &
                                Ptr3d_4   = current%GCStateData3d_4, &
                                Ptr3d_8   = current%GCStateData3d,   &
                                Ptr3d_I   = current%GCStateData3d_I, &
                                RC        = RC                      )
       ENDIF
       IF ( RC == GC_FAILURE ) THEN
          ErrMsg = "Registry pointer not found for " // TRIM(current%name) // &
                   ". Check that the tag (e.g. species) is valid "         // &
                   "for this diagnostic."
          EXIT
       ENDIF

       ! For MAPL export, need to pass a pointer of the right dimension
       IF ( current%rank == 2 ) THEN
          CALL MAPL_GetPointer ( EXPORT, current%ExportData2d, &
                                 current%name, __RC__ )
       ELSEIF ( current%rank == 3 ) THEN
          CALL MAPL_GetPointer ( EXPORT, current%ExportData3d, &
                                 current%name, __RC__ )
       ENDIF

       IF ( Am_I_Root) THEN
          WRITE(6,*) TRIM(current%name)
       ENDIF

       current => current%next    
    ENDDO
    current => NULL()

    IF ( RC == GC_FAILURE ) THEN
       CALL GC_ERROR( ErrMsg, RC, Iam )
       _ASSERT(.FALSE., 'informative message here')
       RETURN
    ENDIF
    
  END SUBROUTINE HistoryExports_SetDataPointers
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Destroy_HistoryConfig 
!
! !DESCRIPTION: Subroutine Destroy_HistoryConfig deallocates a HistoryConfig
!  object and all of its member objects including the linked list of
!  HistoryExport objects.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Destroy_HistoryConfig ( am_I_Root, HistoryConfig, RC )
!
! !INPUT PARAMETERS:
!
    LOGICAL,            INTENT(IN)    :: am_I_Root     ! root CPU?
    TYPE(HistoryConfigObj), POINTER   :: HistoryConfig
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,            INTENT(INOUT) :: RC            ! Success?
!
! !REVISION HISTORY:
!  01 Sep 2017 - E. Lundgren - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    TYPE(HistoryExportObj), POINTER :: current
    TYPE(HistoryExportObj), POINTER :: next

    ! ================================================================
    ! Destroy_HistoryConfig begins here
    ! ================================================================
    __Iam__('Destroy_HistoryConfig (gigc_historyexports_mod.F90)')

    current => NULL()
    next => NULL()

    ! Destroy each item in the linked list of HistoryExport objects
    current => HistoryConfig%HistoryExportsList%head
    IF ( ASSOCIATED( current ) ) next => current%next
    DO WHILE ( ASSOCIATED( current ) )
       DEALLOCATE( current, STAT=RC )
       _ASSERT( RC == GC_SUCCESS, 'informative message here' )
       IF ( .NOT. ASSOCIATED ( next ) ) EXIT
       current => next
       next => current%next
    ENDDO

    ! Deallocate the HistoryExportsList object
    DEALLOCATE( HistoryConfig%HistoryExportsList, STAT=RC )
    _ASSERT( RC == GC_SUCCESS, 'informative message here' )

    ! Deallocate the HistoryConfig object
    DEALLOCATE( HistoryConfig, STAT=RC )
    _ASSERT( RC == GC_SUCCESS, 'informative message here' )

    ! Final cleanup
    current => NULL()
    next    => NULL()

  END SUBROUTINE Destroy_HistoryConfig
!EOC
END MODULE GIGC_HistoryExports_Mod
