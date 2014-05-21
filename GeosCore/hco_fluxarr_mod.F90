!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: hco_fluxarr_mod
!
! !DESCRIPTION: Module HCO\_FLUXARR\_MOD contains routines to handle 
! the HEMCO flux arrays. These are the emissions and deposition arrays
! listed in the HEMCO state object. 
!\\
!\\
! !INTERFACE: 
!
MODULE HCO_FLUXARR_MOD
!
! USES:
!
  USE HCO_ERROR_MOD
  USE HCO_ARR_MOD
  USE HCO_STATE_MOD, ONLY : HCO_State

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC :: HCO_EmisAdd
  PUBLIC :: HCO_DepvAdd
  PUBLIC :: HCO_FluxarrReset
!
! !REMARKS:
!                                                                             
! !REVISION HISTORY:
!  05 Jan 2014 - C. Keller - Initial version, adapted from hco_state_mod.F90 
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !MODULE INTERFACES:
!
  INTERFACE HCO_EmisAdd
     MODULE PROCEDURE HCO_EmisAdd_3D_DP 
     MODULE PROCEDURE HCO_EmisAdd_3D_SP 
     MODULE PROCEDURE HCO_EmisAdd_2D_SP
     MODULE PROCEDURE HCO_EmisAdd_2D_DP
     MODULE PROCEDURE HCO_EmisAdd_DP
     MODULE PROCEDURE HCO_EmisAdd_SP
  END INTERFACE

  INTERFACE HCO_DepvAdd
     MODULE PROCEDURE HCO_DepvAdd_2D_SP
     MODULE PROCEDURE HCO_DepvAdd_2D_DP
     MODULE PROCEDURE HCO_DepvAdd_DP
     MODULE PROCEDURE HCO_DepvAdd_SP
  END INTERFACE

!
! !PRIVATE MEMBER FUNCTIONS:
!
  CONTAINS
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCO_FluxarrReset 
!
! !DESCRIPTION: Routine HCO\_FluxarrReset (re)sets all defined 3D flux arrays 
! (Emis) and 2D deposition arrays (Depv) of the passed HEMCO object. The
! (optional) argument FlxDir indicates whether only emissions (+1) or 
! deposition (-1) arrays shall be reset, or both (0, default).
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCO_FluxarrReset ( HcoState, RC, FlxDir )
!
! !INPUT/OUTPUT PARAMETERS:
! 
    TYPE(HCO_State), INTENT(INOUT)           :: HcoState
    INTEGER,         INTENT(INOUT)           :: RC 
    INTEGER,         INTENT(IN   ), OPTIONAL :: FlxDir    
!
! !REMARKS:
!
! !REVISION HISTORY: 
!  01 May 2013 - C. Keller - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
    INTEGER :: N, locFlxDir 

    !=====================================================================
    ! HCO_FluxarrReset begins here!
    !=====================================================================

    ! Set local flux direction flag
    IF ( PRESENT(FlxDir) ) THEN
       locFlxDir = FlxDir
    ELSE
       locFlxDir = 0
    ENDIF

    ! Loop over all arrays. 
    DO N = 1, HcoState%nSpc

       ! 3D flux rates array
       IF ( locFlxDir >= 0 ) THEN
          IF ( ASSOCIATED(HcoState%Spc(N)%Emis) ) THEN
             IF ( ASSOCIATED(HcoState%Spc(N)%Emis%Val) ) THEN
                HcoState%Spc(N)%Emis%Val = 0.0_hp
             ENDIF
          ENDIF
       ENDIF

       ! 2D deposition velocity array
       IF ( locFlxDir <= 0 ) THEN
          IF ( ASSOCIATED(HcoState%Spc(N)%Depv) ) THEN
             IF ( ASSOCIATED(HcoState%Spc(N)%Depv%Val) ) THEN
                HcoState%Spc(N)%Depv%Val = 0.0_hp
             ENDIF
          ENDIF
       ENDIF
    ENDDO

    ! Return w/ success
    RC = HCO_SUCCESS

  END SUBROUTINE HCO_FluxarrReset 
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCO_EmisAdd_3D_DP
!
! !DESCRIPTION: Routine HCO\_EmisAdd\_3D adds the 3D-array Arr3D 
! to the emissions array of species HcoID in HEMCO object HcoState.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCO_EmisAdd_3D_DP ( HcoState, Arr3D, HcoID, RC )
!
! !INPUT/OUTPUT PARAMETERS:
! 
    TYPE(HCO_State), INTENT(INOUT)   :: HcoState
    REAL(dp),        INTENT(IN   )   :: Arr3D( HcoState%NX, &
                                               HcoState%NY, &
                                               HcoState%NZ )
    INTEGER,         INTENT(IN   )   :: HcoID 
    INTEGER,         INTENT(INOUT)   :: RC 
!
! !REMARKS:
!
! !REVISION HISTORY: 
!  01 May 2013 - C. Keller - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC

    !=====================================================================
    ! HCO_EmisAdd_3D_DP begins here!
    !=====================================================================

    ! Make sure target flux array in HcoState is allocated
    CALL HCO_ArrAssert ( HcoState%Spc(HcoID)%Emis, &
                         HcoState%NX, HcoState%NY, HcoState%NZ, RC )
    IF ( RC /=HCO_SUCCESS ) RETURN

    ! Add array
    HcoState%Spc(HcoID)%Emis%Val(:,:,:) = &
       HcoState%Spc(HcoID)%Emis%Val(:,:,:) + Arr3D(:,:,:)

    ! Return w/ success
    RC = HCO_SUCCESS 

  END SUBROUTINE HCO_EmisAdd_3D_DP
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCO_EmisAdd_3D_SP
!
! !DESCRIPTION: Routine HCO\_EmisAdd\_3D adds the 3D-array Arr3D to the 
! emissions array of species HcoID in HEMCO object HcoState.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCO_EmisAdd_3D_SP ( HcoState, Arr3D, HcoID, RC )
!
! !INPUT/OUTPUT PARAMETERS:
! 
    TYPE(HCO_State), INTENT(INOUT)   :: HcoState
    REAL(sp),        INTENT(IN   )   :: Arr3D( HcoState%NX, &
                                               HcoState%NY, &
                                               HcoState%NZ )
    INTEGER,         INTENT(IN   )   :: HcoID 
    INTEGER,         INTENT(INOUT)   :: RC 
!
! !REMARKS:
!
! !REVISION HISTORY: 
!  01 May 2013 - C. Keller - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC

    !=====================================================================
    ! HCO_EmisAdd_3D_SP begins here!
    !=====================================================================

    ! Make sure target flux array in HcoState is allocated
    CALL HCO_ArrAssert ( HcoState%Spc(HcoID)%Emis, &
                         HcoState%NX, HcoState%NY, HcoState%NZ, RC )
    IF ( RC /=HCO_SUCCESS ) RETURN

    ! Add array
    HcoState%Spc(HcoID)%Emis%Val(:,:,:) = &
       HcoState%Spc(HcoID)%Emis%Val(:,:,:) + Arr3D(:,:,:)

    ! Return w/ success
    RC = HCO_SUCCESS

  END SUBROUTINE HCO_EmisAdd_3D_SP
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCO_EmisAdd_2D_DP
!
! !DESCRIPTION: Routine HCO\_EmisAdd\_2D\_DP adds the real*8 2D-array Arr2D 
! to the emission array of species HcoID in HEMCO object HcoState.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCO_EmisAdd_2D_DP ( HcoState, Arr2D, HcoID, RC )
!
! !INPUT/OUTPUT PARAMETERS:
! 
    TYPE(HCO_State), INTENT(INOUT)   :: HcoState
    REAL(dp),        INTENT(IN   )   :: Arr2D(HcoState%NX,HcoState%NY)
    INTEGER,         INTENT(IN   )   :: HcoID 
    INTEGER,         INTENT(INOUT)   :: RC 
!
! !REMARKS:
!
! !REVISION HISTORY: 
!  01 May 2013 - C. Keller - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC

    !=====================================================================
    ! HCO_EmisAdd_2D_DP begins here!
    !=====================================================================

    ! Make sure target flux array in HcoState is allocated
    CALL HCO_ArrAssert ( HcoState%Spc(HcoID)%Emis, &
                         HcoState%NX, HcoState%NY, HcoState%NZ, RC )
    IF ( RC /=HCO_SUCCESS ) RETURN

    ! Add array
    HcoState%Spc(HcoID)%Emis%Val(:,:,1) = &
       HcoState%Spc(HcoID)%Emis%Val(:,:,1) + Arr2D(:,:)

    ! Return w/ success
    RC = HCO_SUCCESS 

  END SUBROUTINE HCO_EmisAdd_2D_DP
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCO_EmisAdd_2D_SP
!
! !DESCRIPTION: Routine HCO\_EmisAdd\_2D\_SP adds the real*4 2D-array Arr2D 
! to the emission array of species HcoID in HEMCO object HcoState.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCO_EmisAdd_2D_SP ( HcoState, Arr2D, HcoID, RC )
!
! !INPUT/OUTPUT PARAMETERS:
! 
    TYPE(HCO_State), INTENT(INOUT)   :: HcoState
    REAL(sp),        INTENT(IN   )   :: Arr2D(HcoState%NX,HcoState%NY)
    INTEGER,         INTENT(IN   )   :: HcoID 
    INTEGER,         INTENT(INOUT)   :: RC 
!
! !REMARKS:
!
! !REVISION HISTORY: 
!  01 May 2013 - C. Keller - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC

    !=====================================================================
    ! HCO_EmisAdd_2D_SP begins here!
    !=====================================================================

    ! Make sure target flux array in HcoState is allocated
    CALL HCO_ArrAssert ( HcoState%Spc(HcoID)%Emis, &
                         HcoState%NX, HcoState%NY, HcoState%NZ, RC )
    IF ( RC /=HCO_SUCCESS ) RETURN

    ! Add array
    HcoState%Spc(HcoID)%Emis%Val(:,:,1) = &
       HcoState%Spc(HcoID)%Emis%Val(:,:,1) + Arr2D(:,:)

    ! Return w/ success
    RC = HCO_SUCCESS

  END SUBROUTINE HCO_EmisAdd_2D_SP
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCO_EmisAdd_DP
!
! !DESCRIPTION: Routine HCO\_EmisAdd\_DP adds value iVal to the emission
! array of species HcoID in HEMCO object HcoState. The value is placed at 
! location I, J, L of the array.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCO_EmisAdd_DP ( HcoState, iVal, HcoID, I, J, L, RC )
!
! !INPUT/OUTPUT PARAMETERS:
! 
    TYPE(HCO_State), INTENT(INOUT)   :: HcoState
    REAL(dp),        INTENT(IN   )   :: iVal 
    INTEGER,         INTENT(IN   )   :: HcoID 
    INTEGER,         INTENT(IN   )   :: I 
    INTEGER,         INTENT(IN   )   :: J 
    INTEGER,         INTENT(IN   )   :: L 
    INTEGER,         INTENT(INOUT)   :: RC 
!
! !REMARKS:
!
! !REVISION HISTORY: 
!  01 May 2013 - C. Keller - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
    CHARACTER(LEN=255) :: MSG

    !=====================================================================
    ! HCO_EmisAdd_DP begins here!
    !=====================================================================

    ! Check size dimensions
    IF ( I > HcoState%NX ) THEN
       MSG = 'Cannot add DP - i too high!'
       CALL HCO_ERROR ( MSG, RC )
       RETURN
    ENDIF
    IF ( J > HcoState%NY ) THEN
       MSG = 'Cannot add DP - j too high!'
       CALL HCO_ERROR ( MSG, RC )
       RETURN
    ENDIF
    IF ( L > HcoState%NZ ) THEN
       MSG = 'Cannot add DP - l too high!'
       CALL HCO_ERROR ( MSG, RC )
       RETURN
    ENDIF

    ! Make sure target flux array in HcoState is allocated
    CALL HCO_ArrAssert ( HcoState%Spc(HcoID)%Emis, &
                         HcoState%NX, HcoState%NY, HcoState%NZ, RC )
    IF ( RC /=HCO_SUCCESS ) RETURN

    ! Add array
    HcoState%Spc(HcoID)%Emis%Val(I,J,L) = &
       HcoState%Spc(HcoID)%Emis%Val(I,J,L) + iVal

    ! Return w/ success
    RC = HCO_SUCCESS 

  END SUBROUTINE HCO_EmisAdd_DP
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCO_EmisAdd_SP
!
! !DESCRIPTION: Routine HCO\_EmisAdd\_SP adds value iVal to the emission
! array of species HcoID in HEMCO object HcoState. The value is placed 
! at location I, J, L of the array.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCO_EmisAdd_SP ( HcoState, iVal, HcoID, I, J, L, RC )
!
! !INPUT/OUTPUT PARAMETERS:
! 
    TYPE(HCO_State), INTENT(INOUT)   :: HcoState
    REAL(sp),        INTENT(IN   )   :: iVal 
    INTEGER,         INTENT(IN   )   :: HcoID 
    INTEGER,         INTENT(IN   )   :: I 
    INTEGER,         INTENT(IN   )   :: J 
    INTEGER,         INTENT(IN   )   :: L 
    INTEGER,         INTENT(INOUT)   :: RC 
!
! !REMARKS:
!
! !REVISION HISTORY: 
!  01 May 2013 - C. Keller - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
    CHARACTER(LEN=255) :: MSG

    !=====================================================================
    ! HCO_EmisAdd_SP begins here!
    !=====================================================================

    ! Check size dimensions
    IF ( I > HcoState%NX ) THEN
       MSG = 'Cannot add SP - i too high!'
       CALL HCO_ERROR ( MSG, RC )
       RETURN
    ENDIF
    IF ( J > HcoState%NY ) THEN
       MSG = 'Cannot add SP - j too high!'
       CALL HCO_ERROR ( MSG, RC )
       RETURN
    ENDIF
    IF ( L > HcoState%NZ ) THEN
       MSG = 'Cannot add SP - l too high!'
       CALL HCO_ERROR ( MSG, RC )
       RETURN
    ENDIF

    ! Make sure target flux array in HcoState is allocated
    CALL HCO_ArrAssert ( HcoState%Spc(HcoID)%Emis, &
                         HcoState%NX, HcoState%NY, HcoState%NZ, RC )
    IF ( RC /=HCO_SUCCESS ) RETURN

    ! Add array
    HcoState%Spc(HcoID)%Emis%Val(I,J,L) = &
       HcoState%Spc(HcoID)%Emis%Val(I,J,L) + iVal

    ! Return w/ success
    RC = HCO_SUCCESS

  END SUBROUTINE HCO_EmisAdd_SP
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCO_DepvAdd_2D_DP
!
! !DESCRIPTION: Routine HCO\_DepvAdd\_2D\_DP adds the real*8 2D-array Arr2D 
! to the depostion array of species HcoID in HEMCO object HcoState.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCO_DepvAdd_2D_DP ( HcoState, Arr2D, HcoID, RC )
!
! !INPUT/OUTPUT PARAMETERS:
! 
    TYPE(HCO_State), INTENT(INOUT)   :: HcoState
    REAL(dp),        INTENT(IN   )   :: Arr2D(HcoState%NX,HcoState%NY)
    INTEGER,         INTENT(IN   )   :: HcoID 
    INTEGER,         INTENT(INOUT)   :: RC 
!
! !REMARKS:
!
! !REVISION HISTORY: 
!  01 May 2013 - C. Keller - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
    CHARACTER(LEN=255) :: MSG

    !=====================================================================
    ! HCO_DepvAdd_2D_DP begins here!
    !=====================================================================

    ! Make sure target flux array in HcoState is allocated
    CALL HCO_ArrAssert ( HcoState%Spc(HcoID)%Depv, &
                         HcoState%NX, HcoState%NY, RC )
    IF ( RC /=HCO_SUCCESS ) RETURN

    ! Add array
    HcoState%Spc(HcoID)%Depv%Val(:,:) = &
       HcoState%Spc(HcoID)%Depv%Val(:,:) + Arr2D(:,:)

    ! Return w/ success
    RC = HCO_SUCCESS 

  END SUBROUTINE HCO_DepvAdd_2D_DP
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCO_DepvAdd_2D_SP
!
! !DESCRIPTION: Routine HCO\_DepvAdd\_2D\_SP adds the real*4 2D-array Arr2D 
! to the depostion array of species HcoID in HEMCO object HcoState.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCO_DepvAdd_2D_SP ( HcoState, Arr2D, HcoID, RC )
!
! !INPUT/OUTPUT PARAMETERS:
! 
    TYPE(HCO_State), INTENT(INOUT)   :: HcoState
    REAL(sp),        INTENT(IN   )   :: Arr2D(HcoState%NX,HcoState%NY)
    INTEGER,         INTENT(IN   )   :: HcoID 
    INTEGER,         INTENT(INOUT)   :: RC 
!
! !REMARKS:
!
! !REVISION HISTORY: 
!  01 May 2013 - C. Keller - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC

    !=====================================================================
    ! HCO_DepvAdd_2D_SP begins here!
    !=====================================================================

    ! Make sure target flux array in HcoState is allocated
    CALL HCO_ArrAssert ( HcoState%Spc(HcoID)%Depv, &
                         HcoState%NX, HcoState%NY, RC )
    IF ( RC /=HCO_SUCCESS ) RETURN

    ! Add array
    HcoState%Spc(HcoID)%Depv%Val(:,:) = &
       HcoState%Spc(HcoID)%Depv%Val(:,:) + Arr2D(:,:)

    ! Return w/ success
    RC = HCO_SUCCESS 

  END SUBROUTINE HCO_DepvAdd_2D_SP
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCO_DepvAdd_DP
!
! !DESCRIPTION: Routine HCO\_DepvAdd\_DP adds value iVal to the deposition 
! array of species HcoID in HEMCO object HcoState. The value is placed at 
! location I, J of the array.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCO_DepvAdd_DP ( HcoState, iVal, HcoID, I, J, RC )
!
! !INPUT/OUTPUT PARAMETERS:
! 
    TYPE(HCO_State), INTENT(INOUT)   :: HcoState
    REAL(dp),        INTENT(IN   )   :: iVal 
    INTEGER,         INTENT(IN   )   :: HcoID 
    INTEGER,         INTENT(IN   )   :: I 
    INTEGER,         INTENT(IN   )   :: J 
    INTEGER,         INTENT(INOUT)   :: RC 
!
! !REMARKS:
!
! !REVISION HISTORY: 
!  01 May 2013 - C. Keller - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
    CHARACTER(LEN=255) :: MSG

    !=====================================================================
    ! HCO_DepvAdd_DP begins here!
    !=====================================================================

    ! Check size dimensions
    IF ( I > HcoState%NX ) THEN
       MSG = 'Cannot add DP - i too high!'
       CALL HCO_ERROR ( MSG, RC )
       RETURN
    ENDIF
    IF ( J > HcoState%NY ) THEN
       MSG = 'Cannot add DP - j too high!'
       CALL HCO_ERROR ( MSG, RC )
       RETURN
    ENDIF

    ! Make sure target flux array in HcoState is allocated
    CALL HCO_ArrAssert ( HcoState%Spc(HcoID)%Depv, &
                         HcoState%NX, HcoState%NY, RC )
    IF ( RC /=HCO_SUCCESS ) RETURN

    ! Add array
    HcoState%Spc(HcoID)%Depv%Val(I,J) = &
       HcoState%Spc(HcoID)%Depv%Val(I,J) + iVal 

    ! Return w/ success
    RC = HCO_SUCCESS

  END SUBROUTINE HCO_DepvAdd_DP
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCO_DepvAdd_SP
!
! !DESCRIPTION: Routine HCO\_DepvAdd\_SP adds value iVal to the deposition 
! array of species HcoID in HEMCO object HcoState. The value is placed at 
! location I, J of the array.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCO_DepvAdd_SP ( HcoState, iVal, HcoID, I, J, RC )
!
! !INPUT/OUTPUT PARAMETERS:
! 
    TYPE(HCO_State), INTENT(INOUT)   :: HcoState
    REAL(sp),        INTENT(IN   )   :: iVal 
    INTEGER,         INTENT(IN   )   :: HcoID 
    INTEGER,         INTENT(IN   )   :: I 
    INTEGER,         INTENT(IN   )   :: J 
    INTEGER,         INTENT(INOUT)   :: RC 
!
! !REMARKS:
!
! !REVISION HISTORY: 
!  01 May 2013 - C. Keller - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
    CHARACTER(LEN=255) :: MSG

    !=====================================================================
    ! HCO_DepvAdd_SP begins here!
    !=====================================================================

    ! Check size dimensions
    IF ( I > HcoState%NX ) THEN
       MSG = 'Cannot add iVal - i too high!'
       CALL HCO_ERROR ( MSG, RC )
       RETURN
    ENDIF
    IF ( J > HcoState%NY ) THEN
       MSG = 'Cannot add iVal - j too high!'
       CALL HCO_ERROR ( MSG, RC )
       RETURN
    ENDIF

    ! Make sure target flux array in HcoState is allocated
    CALL HCO_ArrAssert ( HcoState%Spc(HcoID)%Depv, &
                         HcoState%NX, HcoState%NY, RC )
    IF ( RC /=HCO_SUCCESS ) RETURN

    ! Add array
    HcoState%Spc(HcoID)%Depv%Val(I,J) = &
       HcoState%Spc(HcoID)%Depv%Val(I,J) + iVal

    ! Return w/ success
    RC = HCO_SUCCESS 

  END SUBROUTINE HCO_DepvAdd_SP
!EOC
END MODULE HCO_FLUXARR_MOD
