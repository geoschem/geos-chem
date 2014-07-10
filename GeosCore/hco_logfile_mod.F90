!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: hco_logfile_mod
!
! !DESCRIPTION: Module HCO\_LOGFILE\_MOD contains some routines to write
! data into the HEMCO logfile.
!\\
!\\
! !INTERFACE: 
!
MODULE HCO_LOGFILE_MOD 
! 
! !USES:
!
  USE HCO_ERROR_MOD

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC  :: HCO_SPEC2LOG
!
! !REVISION HISTORY:
!  27 May 2014 - C. Keller   - Initialization
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
! !IROUTINE: hco_spec2log
!
! !DESCRIPTION: Subroutine HCO\_SPEC2LOG writes information of a species
! to the logfile.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCO_SPEC2LOG( am_I_Root, HcoState, ID )
!
! !USES:
!
    USE HCO_STATE_MOD, ONLY : HCO_State
!
!
! !INPUT PARAMETER
!
    LOGICAL,          INTENT(IN)     :: am_I_Root  ! Root CPU
    TYPE(HCO_State),  POINTER        :: HcoState   ! HEMCO state object
    INTEGER,          INTENT(IN)     :: ID         ! HEMCO species ID
!
! !REVISION HISTORY: 
!  27 May 2014 - C. Keller   - Initialization
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    LOGICAL             :: verb
    CHARACTER(LEN=255)  :: MSG 

    !=================================================================
    ! HCO_SPEC2LOG begins here 
    !=================================================================

    ! Verbose mode?
    verb = HCO_VERBOSE_CHECK() .AND. am_I_Root

    MSG = 'Species ' // TRIM(HcoState%Spc(ID)%SpcName)
    CALL HCO_MSG(MSG)
    IF ( verb ) THEN
       write(MSG,*) '--> HcoID         : ', HcoState%Spc(ID)%HcoID
       CALL HCO_MSG(MSG)
       write(MSG,*) '--> ModID         : ', HcoState%Spc(ID)%ModID
       CALL HCO_MSG(MSG)
       write(MSG,*) '--> MW (g/mol)    : ', HcoState%Spc(ID)%MW_g
       CALL HCO_MSG(MSG)
       write(MSG,*) '--> emitted MW    : ', HcoState%Spc(ID)%EmMW_g
       CALL HCO_MSG(MSG)
       write(MSG,*) '--> Molecule ratio: ', HcoState%Spc(ID)%MolecRatio
       CALL HCO_MSG(MSG)
       write(MSG,*) '--> Henry constant: ', HcoState%Spc(ID)%HenryK0
       CALL HCO_MSG(MSG)
       write(MSG,*) '--> Henry temp.   : ', HcoState%Spc(ID)%HenryCR
       CALL HCO_MSG(MSG)
       write(MSG,*) '--> Henry pKA     : ', HcoState%Spc(ID)%HenryPKA
       CALL HCO_MSG(MSG)
    ENDIF    

  END SUBROUTINE HCO_SPEC2LOG
!EOC
END MODULE HCO_LOGFILE_MOD
