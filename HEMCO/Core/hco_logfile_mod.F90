!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: hco_logfile_mod
!
! !DESCRIPTION: Module HCO\_LOGFILE\_MOD contains some wrapper routines to
! write data into the HEMCO logfile.
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
  PUBLIC  :: HCO_Spec2Log
  PUBLIC  :: HCO_PrintList
  PUBLIC  :: HCO_PrintDataCont
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
! !DESCRIPTION: Subroutine HCO\_Spec2Log writes information of a species
! to the logfile.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCO_Spec2Log( HcoState, ID )
!
! !USES:
!
    USE HCO_STATE_MOD, ONLY : HCO_State
!
!
! !INPUT PARAMETER
!
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
    CHARACTER(LEN=255)  :: MSG

    !=================================================================
    ! HCO_Spec2Log begins here
    !=================================================================

    MSG = 'Species ' // TRIM(HcoState%Spc(ID)%SpcName)
    CALL HCO_MSG(HcoState%Config%Err,MSG)
    IF ( HCO_IsVerb(HcoState%Config%Err,3) ) THEN
       write(MSG,*) '--> HcoID         : ', HcoState%Spc(ID)%HcoID
       CALL HCO_MSG(HcoState%Config%Err,MSG)
       write(MSG,*) '--> ModID         : ', HcoState%Spc(ID)%ModID
       CALL HCO_MSG(HcoState%Config%Err,MSG)
       write(MSG,*) '--> MW (g/mol)    : ', HcoState%Spc(ID)%MW_g
       CALL HCO_MSG(HcoState%Config%Err,MSG)
       write(MSG,*) '--> emitted MW    : ', HcoState%Spc(ID)%EmMW_g
       CALL HCO_MSG(HcoState%Config%Err,MSG)
       write(MSG,*) '--> Molecule ratio: ', HcoState%Spc(ID)%MolecRatio
       CALL HCO_MSG(HcoState%Config%Err,MSG)
       write(MSG,*) '--> Henry constant: ', HcoState%Spc(ID)%HenryK0
       CALL HCO_MSG(HcoState%Config%Err,MSG)
       write(MSG,*) '--> Henry temp.   : ', HcoState%Spc(ID)%HenryCR
       CALL HCO_MSG(HcoState%Config%Err,MSG)
       write(MSG,*) '--> Henry pKA     : ', HcoState%Spc(ID)%HenryPKA
       CALL HCO_MSG(HcoState%Config%Err,MSG)
    ENDIF

  END SUBROUTINE HCO_Spec2Log
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: HCO_PrintList
!
! !DESCRIPTION: Subroutine HCO\_PrintList displays the content of List.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE HCO_PrintList ( HcoState, List, Verbose )
!
! !USES:
!
      USE HCO_STATE_MOD,     ONLY : HCO_State
      USE HCO_TYPES_MOD,     ONLY : ListCont
!
! !INPUT ARGUMENTS:
!
      TYPE(HCO_STATE),POINTER    :: HcoState
      TYPE(ListCont), POINTER    :: List
      INTEGER,        INTENT(IN) :: Verbose
!
! !REVISION HISTORY:
!  20 Apr 2013 - C. Keller - Initial version
!  26 Oct 2016 - R. Yantosca - Don't nullify local ptrs in declaration stmts
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !ARGUMENTS:
!
      TYPE(ListCont), POINTER   :: TmpLct
      CHARACTER(LEN=255)        :: MSG

      ! ================================================================
      ! HCO_PrintList begins here
      ! ================================================================

      ! Point to first element
      TmpLct => List
      DO WHILE ( ASSOCIATED(TmpLct) )
         IF ( ASSOCIATED(TmpLct%Dct) ) THEN
            CALL HCO_PrintDataCont(HcoState,TmpLct%Dct,Verbose)
         ENDIF
         TmpLct => TmpLct%NextCont
      ENDDO

      TmpLct => NULL()

      END SUBROUTINE HCO_PrintList
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: HCO_PrintDataCont
!
! !DESCRIPTION: Subroutine HCO\_PrintDataCont displays the content of the
! data container Dct.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE HCO_PrintDataCont ( HcoState, Dct, Verbose )
!
! !USES
!
      USE HCO_STATE_MOD,     ONLY : HCO_State
      USE HCO_TYPES_MOD,     ONLY : DataCont, HCO_DCTTYPE_BASE
!
! !INPUT ARGUMENTS:
!
      TYPE(HCO_STATE),POINTER    :: HcoState
      TYPE(DataCont), POINTER    :: Dct
      INTEGER,        INTENT(IN) :: Verbose
!
! !REVISION HISTORY:
!  20 Apr 2013 - C. Keller - Initial version
!  16 Mar 2015 - M. Sulprizio- Now print min and max values for debugging
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !ARGUMENTS:
!
      CHARACTER(LEN=255) :: MSG
      INTEGER            :: nx, ny, nz, nt
      REAL(sp)           :: sm, mn, mx

      ! ================================================================
      ! HCO_PrintDataCont begins here
      ! ================================================================

      sm = 0.0_sp
      mn = 0.0_sp
      mx = 0.0_sp
      sm = 0.0_sp
      nx = 0
      ny = 0
      nz = 0
      nt = Dct%Dta%nt
      IF ( nt > 0 ) THEN
         IF ( Dct%Dta%spaceDim<=2 ) THEN
            IF ( ASSOCIATED(Dct%Dta%V2) ) THEN
               nx = SIZE(Dct%Dta%V2(1)%Val,1)
               ny = SIZE(Dct%Dta%V2(1)%Val,2)
               sm = SUM(Dct%Dta%V2(1)%Val)
               mn = MINVAL(Dct%Dta%V2(1)%Val)
               mx = MAXVAL(Dct%Dta%V2(1)%Val)
            ENDIF
         ELSE
            IF ( ASSOCIATED(Dct%Dta%V3) ) THEN
               nx = SIZE(Dct%Dta%V3(1)%Val,1)
               ny = SIZE(Dct%Dta%V3(1)%Val,2)
               nz = SIZE(Dct%Dta%V3(1)%Val,3)
               sm = SUM(Dct%Dta%V3(1)%Val)
               mn = MINVAL(Dct%Dta%V3(1)%Val)
               mx = MAXVAL(Dct%Dta%V3(1)%Val)
            ENDIF
         ENDIF
      ENDIF

      ! Print name for verbose > 0
      IF ( HCO_IsVerb(HcoState%Config%Err,1) ) THEN
         MSG = 'Container ' // TRIM(Dct%cName)
         CALL HCO_MSG(HcoState%Config%Err,MSG)
      ENDIF

      ! Eventually add details
      IF ( HCO_IsVerb(HcoState%Config%Err,Verbose) ) THEN

         ! General information
         write(MSG,*) '   -->Data type       : ', Dct%DctType
         CALL HCO_MSG(HcoState%Config%Err,MSG)
         write(MSG,*) '   -->Container ID    : ', Dct%cID
         CALL HCO_MSG(HcoState%Config%Err,MSG)
         write(MSG,*) '   -->Target ID       : ', Dct%targetID
         CALL HCO_MSG(HcoState%Config%Err,MSG)
         write(MSG,*) '   -->File data home?   ', Dct%DtaHome
         CALL HCO_MSG(HcoState%Config%Err,MSG)
         write(MSG,*) '   -->Source file     : ', TRIM(Dct%Dta%ncFile)
         CALL HCO_MSG(HcoState%Config%Err,MSG)
         write(MSG,*) '   -->ncRead?           ', Dct%Dta%ncRead
         CALL HCO_MSG(HcoState%Config%Err,MSG)
         write(MSG,*) '   -->Shared data file? ', Dct%Dta%DoShare
         CALL HCO_MSG(HcoState%Config%Err,MSG)
         IF ( Dct%Dta%ncRead ) THEN
            write(MSG,*) '   -->Source parameter: ', TRIM(Dct%Dta%ncPara)
            CALL HCO_MSG(HcoState%Config%Err,MSG)
            write(MSG,*) '   -->Year range      : ', Dct%Dta%ncYrs
            CALL HCO_MSG(HcoState%Config%Err,MSG)
            write(MSG,*) '   -->Month range     : ', Dct%Dta%ncMts
            CALL HCO_MSG(HcoState%Config%Err,MSG)
            write(MSG,*) '   -->Day range       : ', Dct%Dta%ncDys
            CALL HCO_MSG(HcoState%Config%Err,MSG)
            write(MSG,*) '   -->Hour range      : ', Dct%Dta%ncHrs
            CALL HCO_MSG(HcoState%Config%Err,MSG)
            write(MSG,*) '   -->SpaceDim        : ', Dct%Dta%SpaceDim
            CALL HCO_MSG(HcoState%Config%Err,MSG)
         ENDIF
         IF ( NZ > 0 ) THEN
            write(MSG,*) '   -->Array dimension : ', nx,ny,nz
         ELSE
            write(MSG,*) '   -->Array dimension : ', nx,ny
         ENDIF
         CALL HCO_MSG(HcoState%Config%Err,MSG)
         write(MSG,*) '   -->Array sum       : ', sm
         CALL HCO_MSG(HcoState%Config%Err,MSG)
         write(MSG,*) '   -->Array min & max : ', mn,mx
         CALL HCO_MSG(HcoState%Config%Err,MSG)
         write(MSG,*) '   -->Time dimension  : ', nt
         CALL HCO_MSG(HcoState%Config%Err,MSG)
         write(MSG,*) '   -->Delta t[h]      : ', Dct%Dta%DeltaT
         CALL HCO_MSG(HcoState%Config%Err,MSG)
         write(MSG,*) '   -->Local time?       ', Dct%Dta%IsLocTime
         CALL HCO_MSG(HcoState%Config%Err,MSG)
         IF ( ASSOCIATED(Dct%Dta%tIDx) ) THEN
            write(MSG,*) '   -->Tempres         : ', &
               TRIM(Dct%Dta%tIDx%TempRes)
            CALL HCO_MSG(HcoState%Config%Err,MSG)
         ENDIF
         write(MSG,*) '   -->OrigUnit        : ',TRIM(Dct%Dta%OrigUnit)
         CALL HCO_MSG(HcoState%Config%Err,MSG)
         write(MSG,*) '   -->Concentration?    ', Dct%Dta%IsConc
         CALL HCO_MSG(HcoState%Config%Err,MSG)
         write(MSG,*) '   -->Coverage        : ', Dct%Dta%Cover
         CALL HCO_MSG(HcoState%Config%Err,MSG)

         ! For base emissions
         IF ( Dct%DctType==HCO_DCTTYPE_BASE ) THEN
            write(MSG,*) '   -->Extension Nr    : ', Dct%ExtNr
            CALL HCO_MSG(HcoState%Config%Err,MSG)
            write(MSG,*) '   -->Species name    : ',TRIM(Dct%SpcName)
            CALL HCO_MSG(HcoState%Config%Err,MSG)
            write(MSG,*) '   -->HEMCO species ID: ', Dct%HcoID
            CALL HCO_MSG(HcoState%Config%Err,MSG)
            write(MSG,*) '   -->Category        : ', Dct%Cat
            CALL HCO_MSG(HcoState%Config%Err,MSG)
            write(MSG,*) '   -->Hierarchy       : ', Dct%Hier
            CALL HCO_MSG(HcoState%Config%Err,MSG)
            IF ( nz == 0 ) THEN
               write(MSG,*) '   -->2D emitted into : ', Dct%Dta%EmisL1, &
                            ' and ', Dct%Dta%EmisL2
               CALL HCO_MSG(HcoState%Config%Err,MSG)
            ENDIF

         ! For scale factors
         ELSE
            write(MSG,*) '   -->Scal ID         : ', Dct%ScalID
            CALL HCO_MSG(HcoState%Config%Err,MSG)
            write(MSG,*) '   -->Operator        : ', Dct%Oper
            CALL HCO_MSG(HcoState%Config%Err,MSG)
         ENDIF
      ENDIF

      END SUBROUTINE HCO_PrintDataCont
!EOC
END MODULE HCO_LOGFILE_MOD
