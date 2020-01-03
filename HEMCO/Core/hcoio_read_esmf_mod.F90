!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: hcoio_read_esmf_mod.F90
!
! !DESCRIPTION: Module HCOIO\_Read\_ESMF\_mod is the HEMCO interface for
!  data reading within the ESMF framework.
!\\
!\\
! !INTERFACE:
!
MODULE HCOIO_Read_ESMF_mod
!
! !USES:
!
  USE HCO_Types_Mod
  USE HCO_Error_Mod
  USE HCO_State_Mod,       ONLY : Hco_State

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
#if defined(ESMF_)
  PUBLIC  :: HCOIO_Read_ESMF
!
! !REVISION HISTORY:
!  22 Aug 2013 - C. Keller   - Initial version
!  01 Jul 2014 - R. Yantosca - Now use F90 free-format indentation
!  01 Jul 2014 - R. Yantosca - Cosmetic changes in ProTeX headers
!  22 Feb 2016 - C. Keller   - Split off from hcoio_dataread_mod.F90
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
! !IROUTINE: HCOIO_DataRead (ESMF/MAPL version)
!
! !DESCRIPTION: Interface routine between ESMF and HEMCO to obtain
! the data array for a given HEMCO data container. The data is obtained
! through the ExtData interface. The HEMCO source file attribute is taken
! to identify the ExtData pointer name.
!\\
!\\
! !INTERFACE:
  !
  SUBROUTINE HCOIO_Read_ESMF( HcoState, Lct, RC )
!
! !USES:
!
    USE ESMF
    USE MAPL_mod
    USE HCO_FILEDATA_MOD, ONLY : FileData_ArrInit

# include "MAPL_Generic.h"
!
! !INPUT PARAMETERS:
!
    TYPE(HCO_State),  POINTER        :: HcoState
    TYPE(ListCont),   POINTER        :: Lct
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(INOUT)  :: RC
!
! !REVISION HISTORY:
!  28 Aug 2013 - C. Keller   - Initial version
!  27 Aug 2014 - R. Yantosca - Err msg now displays hcoio_dataread_mod.F90
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER                    :: II, JJ, LL, TT
    INTEGER                    :: I, J, L, T
    INTEGER                    :: STAT
    REAL,             POINTER  :: Ptr3D(:,:,:)
    REAL,             POINTER  :: Ptr2D(:,:)
    TYPE(ESMF_State), POINTER  :: IMPORT
    CHARACTER(LEN=255)         :: MSG
    CHARACTER(LEN=255), PARAMETER :: LOC = 'HCOIO_READ_ESMF (hcoi_read_esmf_mod.F90)'
    CHARACTER(LEN=ESMF_MAXSTR) :: Iam

    !=================================================================
    ! HCOIO_READ_ESMF begins here
    !=================================================================

    ! For error handling
    Iam = LOC
    CALL HCO_ENTER( HcoState%Config%Err,  LOC, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Point to ESMF IMPORT object
    IMPORT => HcoState%IMPORT
    ASSERT_(ASSOCIATED(IMPORT))

    ! Init pointers
    Ptr3D => NULL()
    Ptr2D => NULL()

    ! Verbose?
    IF ( HCO_IsVerb(HcoState%Config%Err,2) ) THEN
       MSG = 'Reading from ExtData: ' // TRIM(Lct%Dct%Dta%ncFile)
       CALL HCO_MSG(HcoState%Config%Err,MSG)
    ENDIF

    !-----------------------------------------------------------------
    ! Read 3D data from ESMF
    !-----------------------------------------------------------------
    IF ( Lct%Dct%Dta%SpaceDim == 3 ) THEN

       ! Get data
       CALL MAPL_GetPointer( IMPORT, Ptr3D, &
                             TRIM(Lct%Dct%Dta%ncFile), RC=STAT )

       ! Check for MAPL error
       IF( STAT /= ESMF_SUCCESS ) THEN
          MSG = 'Cannot get xyz pointer: ' // TRIM(Lct%Dct%Dta%ncFile)
          CALL HCO_ERROR( HcoState%Config%Err, MSG, RC )
          RETURN
       ENDIF

       ! Get array dimensions
       II = SIZE(Ptr3D,1)
       JJ = SIZE(Ptr3D,2)
       LL = SIZE(Ptr3D,3)
       TT = 1

       ! Define HEMCO array if not yet defined.
       IF ( .NOT. ASSOCIATED(Lct%Dct%Dta%V3) ) THEN

          ! Use pointer if types match
          CALL FileData_ArrInit( Lct%Dct%Dta, TT, 0, 0, 0, RC )
          !CALL FileData_ArrInit( Lct%Dct%Dta, TT, II, JJ, LL, RC )
          IF ( RC /= HCO_SUCCESS ) RETURN
       ENDIF

       ! Pointer to data. HEMCO expects data to have surface level at
       ! index 1 ('up').
       Lct%Dct%Dta%V3(1)%Val => Ptr3D(:,:,LL:1:-1)
       !Lct%Dct%Dta%V3(1)%Val(:,:,:) = Ptr3D(:,:,LL:1:-1)

       ! ewl debugging
       if ( mapl_am_i_root() ) then
          print *, "HEMCO: array pointer vertically flipped relative to MAPL Import ", trim(Lct%Dct%Dta%ncFile)
       endif

    !-----------------------------------------------------------------
    ! Read 2D data from ESMF
    !-----------------------------------------------------------------
    ELSEIF ( Lct%Dct%Dta%SpaceDim == 2 ) THEN

       ! Get data
       CALL MAPL_GetPointer( IMPORT, Ptr2D, &
                             TRIM(Lct%Dct%Dta%ncFile), RC=STAT )

       ! Check for MAPL error
       IF( STAT /= ESMF_SUCCESS ) THEN
          MSG = 'Cannot get xy pointer: ' // TRIM(Lct%Dct%Dta%ncFile)
          CALL HCO_ERROR( HcoState%Config%Err, MSG, RC )
          RETURN
       ENDIF

       ! Get array dimensions
       II = SIZE(Ptr2D,1)
       JJ = SIZE(Ptr2D,2)
       LL = 1
       TT = 1

       ! Define HEMCO array pointer if not yet defined
       IF ( .NOT. ASSOCIATED(Lct%Dct%Dta%V2) ) THEN
          CALL FileData_ArrInit( Lct%Dct%Dta, TT, 0, 0, RC )
          !CALL FileData_ArrInit( Lct%Dct%Dta, TT, II, JJ, RC )
          IF ( RC /= HCO_SUCCESS ) RETURN
       ENDIF

       ! Pointer to data
       Lct%Dct%Dta%V2(1)%Val => Ptr2D
       !Lct%Dct%Dta%V2(1)%Val = Ptr2D

    ENDIF

    !-----------------------------------------------------------------
    ! Cleanup and leave
    !-----------------------------------------------------------------
    Ptr3D  => NULL()
    Ptr2D  => NULL()
    IMPORT => NULL()

    ! Return w/ success
    CALL HCO_LEAVE ( HcoState%Config%Err,  RC )

  END SUBROUTINE HCOIO_Read_ESMF
!EOC
#endif
END MODULE HCOIO_Read_ESMF_mod
