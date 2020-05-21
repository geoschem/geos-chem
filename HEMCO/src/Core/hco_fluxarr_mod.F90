!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: hco_fluxarr_mod.F90
!
! !DESCRIPTION: Module HCO\_FluxArr\_Mod contains routines to handle
! the HEMCO flux arrays. These are the emissions and deposition arrays
! listed in the HEMCO state object.
!\\
!\\
! !INTERFACE:
!
MODULE HCO_FluxArr_Mod
!
! USES:
!
  USE HCO_Error_Mod
  USE HCO_Arr_Mod
  USE HCO_Scale_Mod
  USE HCO_State_Mod, ONLY : HCO_State

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC  :: HCO_EmisAdd
  PUBLIC  :: HCO_DepvAdd
  PUBLIC  :: HCO_FluxarrReset
!
! !PRIVATE MEMBER FUNCTIONS:
!
  PRIVATE :: DiagnCheck
!
! !REMARKS:
!
! !REVISION HISTORY:
!  05 Jan 2014 - C. Keller - Initial version, adapted from hco_state_mod.F90
!  21 Oct 2014 - C. Keller - Added error check for negative values to HCO_EmisAdd
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !PRIVATE TYPES:
!
  INTERFACE HCO_EmisAdd
     MODULE PROCEDURE HCO_EmisAdd_3D_Dp
     MODULE PROCEDURE HCO_EmisAdd_3D_Sp
     MODULE PROCEDURE HCO_EmisAdd_2D_Sp
     MODULE PROCEDURE HCO_EmisAdd_2D_Dp
     MODULE PROCEDURE HCO_EmisAdd_Dp
     MODULE PROCEDURE HCO_EmisAdd_Sp
  END INTERFACE

  INTERFACE HCO_DepvAdd
     MODULE PROCEDURE HCO_DepvAdd_2D_Sp
     MODULE PROCEDURE HCO_DepvAdd_2D_Dp
     MODULE PROCEDURE HCO_DepvAdd_Dp
     MODULE PROCEDURE HCO_DepvAdd_Sp
  END INTERFACE

CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCO_FluxarrReset
!
! !DESCRIPTION: Routine HCO\_FluxarrReset (re)sets all data arrays
! of the passed HEMCO state object. The (optional) argument Typ
! indicates whether only emissions (1), deposition (2), or concentration
! (3) arrays shall be reset. To reset all, set Typ to 0 (default).
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCO_FluxarrReset( HcoState, RC, Typ )
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(HCO_State), POINTER                 :: HcoState
    INTEGER,         INTENT(INOUT)           :: RC
!
! !INPUT PARAMETERS:
!
    INTEGER,         INTENT(IN   ), OPTIONAL :: Typ
!
! !REMARKS:
!
! !REVISION HISTORY:
!  01 May 2013 - C. Keller - Initial version
!  21 Aug 2014 - C. Keller - Added concentration
!EOP
!------------------------------------------------------------------------------
!BOC
    INTEGER :: N, thisTyp

    !=====================================================================
    ! HCO_FluxarrReset begins here!
    !=====================================================================

    ! Set local flux direction flag
    IF ( PRESENT(Typ) ) THEN
       thisTyp = Typ
    ELSE
       thisTyp = 0
    ENDIF

    ! Loop over all arrays.
    DO N = 1, HcoState%nSpc

       ! 3D flux rates array
       IF ( thisTyp == 0 .OR. thisTyp == 1 ) THEN
          IF ( ASSOCIATED(HcoState%Spc(N)%Emis) ) THEN
             IF ( ASSOCIATED(HcoState%Spc(N)%Emis%Val) ) THEN
                HcoState%Spc(N)%Emis%Val = 0.0_hp
             ENDIF
          ENDIF
       ENDIF

       ! 2D deposition velocity array
       IF ( thisTyp == 0 .OR. thisTyp == 2 ) THEN
          IF ( ASSOCIATED(HcoState%Spc(N)%Depv) ) THEN
             IF ( ASSOCIATED(HcoState%Spc(N)%Depv%Val) ) THEN
                HcoState%Spc(N)%Depv%Val = 0.0_hp
             ENDIF
          ENDIF
       ENDIF

       ! 3D concentrations
       IF ( thisTyp == 0 .OR. thisTyp == 3 ) THEN
          IF ( ASSOCIATED(HcoState%Spc(N)%Conc) ) THEN
             IF ( ASSOCIATED(HcoState%Spc(N)%Conc%Val) ) THEN
                HcoState%Spc(N)%Conc%Val = 0.0_hp
             ENDIF
          ENDIF
       ENDIF

    ENDDO

    ! Return w/ success
    RC = HCO_SUCCESS

  END SUBROUTINE HCO_FluxarrReset
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCO_EmisAdd_3D_Dp
!
! !DESCRIPTION: Routine HCO\_EmisAdd\_3D adds the 3D-array Arr3D
! to the emissions array of species HcoID in HEMCO object HcoState.
! This routine also updates all autofill diagnostics that are defined
! for the givne species, extension number, emission category and
! hierarchy.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCO_EmisAdd_3D_Dp( HcoState,  Arr3D,    HcoID, &
                                RC,        ExtNr,    Cat,   Hier,  &
                                MinDiagnLev )
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(HCO_State), POINTER                 :: HcoState
    REAL(dp),        INTENT(INOUT)           :: Arr3D( HcoState%NX, &
                                                       HcoState%NY, &
                                                       HcoState%NZ )
    INTEGER,         INTENT(INOUT)           :: RC
!
! !INPUT PARAMETERS:
!
    INTEGER,         INTENT(IN   )           :: HcoID
    INTEGER,         INTENT(IN   ), OPTIONAL :: ExtNr
    INTEGER,         INTENT(IN   ), OPTIONAL :: Cat
    INTEGER,         INTENT(IN   ), OPTIONAL :: Hier
    INTEGER,         INTENT(IN   ), OPTIONAL :: MinDiagnLev
!
! !REVISION HISTORY:
!  01 May 2013 - C. Keller - Initial version
!  20 Apr 2015 - C. Keller - Added DiagnCheck
!  12 May 2017 - C. Keller - Added option to use uniform scale factor
!EOP
!------------------------------------------------------------------------------
!BOC

    !=====================================================================
    ! HCO_EmisAdd_3D_Dp begins here!
    !=====================================================================

    ! Make sure target flux array in HcoState is allocated
    CALL HCO_ArrAssert ( HcoState%Spc(HcoID)%Emis, &
                         HcoState%NX, HcoState%NY, HcoState%NZ, RC )
    IF ( RC /=HCO_SUCCESS ) RETURN

    ! Check for negative values. NegFlag determines the behavior for
    ! negative values: 2 = allow; 1 = set to zero + warning; 0 = error.
    IF ( HcoState%Options%NegFlag /= 2 ) THEN
       IF ( ANY( Arr3D < 0.0_hp ) ) THEN

          ! Negative flag is 1: set to zero and prompt warning
          IF ( HcoState%Options%NegFlag == 1 ) THEN
             WHERE ( Arr3D < 0.0_hp ) Arr3D = 0.0_hp
             CALL HCO_WARNING ( HcoState%Config%Err, &
                'Negative values found - set to zero!', &
                RC, THISLOC = 'HCO_EmisAdd (HCO_FLUXARR_MOD.F90)' )

          ! Negative flag is 0: return w/ error
          ELSE
             CALL HCO_ERROR ( HcoState%Config%Err, &
                'Negative values found!', &
                RC, THISLOC = 'HCO_EmisAdd (HCO_FLUXARR_MOD.F90)' )
             RETURN
          ENDIF
       ENDIF
    ENDIF

    ! Eventually add universal scale factor
    CALL HCO_ScaleArr( HcoState, HcoID, Arr3D, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Add array
    HcoState%Spc(HcoID)%Emis%Val(:,:,:) = &
       HcoState%Spc(HcoID)%Emis%Val(:,:,:) + Arr3D(:,:,:)

    ! Check for diagnostics
    CALL DiagnCheck( HcoState,   ExtNr=ExtNr, Cat=Cat, &
                     Hier=Hier,  HcoID=HcoID, Arr3D=Arr3D, &
                     MinDiagnLev=MinDiagnLev, RC=RC     )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Return w/ success
    RC = HCO_SUCCESS

  END SUBROUTINE HCO_EmisAdd_3D_Dp
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCO_EmisAdd_3D_Sp
!
! !DESCRIPTION: Routine HCO\_EmisAdd\_3D adds the 3D-array Arr3D to the
! emissions array of species HcoID in HEMCO object HcoState.
! This routine also updates all autofill diagnostics that are defined
! for the givne species, extension number, emission category and
! hierarchy.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCO_EmisAdd_3D_Sp ( HcoState,  Arr3D,    HcoID, &
                                 RC,        ExtNr,    Cat,   Hier,  &
                                 MinDiagnLev )
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(HCO_State), POINTER                 :: HcoState
    REAL(sp),        INTENT(INOUT)           :: Arr3D( HcoState%NX, &
                                                       HcoState%NY, &
                                                       HcoState%NZ )
    INTEGER,         INTENT(INOUT)           :: RC
!
! !INPUT PARAMETERS:
!
    INTEGER,         INTENT(IN   )           :: HcoID
    INTEGER,         INTENT(IN   ), OPTIONAL :: ExtNr
    INTEGER,         INTENT(IN   ), OPTIONAL :: Cat
    INTEGER,         INTENT(IN   ), OPTIONAL :: Hier
    INTEGER,         INTENT(IN   ), OPTIONAL :: MinDiagnLev
!
! !REVISION HISTORY:
!  01 May 2013 - C. Keller - Initial version
!  20 Apr 2015 - C. Keller - Added DiagnCheck
!  12 May 2017 - C. Keller - Added option to use uniform scale factor
!EOP
!------------------------------------------------------------------------------
!BOC

    !=====================================================================
    ! HCO_EmisAdd_3D_Sp begins here!
    !=====================================================================

    ! Make sure target flux array in HcoState is allocated
    CALL HCO_ArrAssert ( HcoState%Spc(HcoID)%Emis, &
                         HcoState%NX, HcoState%NY, HcoState%NZ, RC )
    IF ( RC /=HCO_SUCCESS ) RETURN

    ! Check for negative values. NegFlag determines the behavior for
    ! negative values: 2 = allow; 1 = set to zero + warning; 0 = error.
    IF ( HcoState%Options%NegFlag /= 2 ) THEN
       IF ( ANY( Arr3D < 0.0_sp ) ) THEN

          ! Negative flag is 1: set to zero and prompt warning
          IF ( HcoState%Options%NegFlag == 1 ) THEN
             WHERE ( Arr3D < 0.0_sp ) Arr3D = 0.0_sp
             CALL HCO_WARNING ( HcoState%Config%Err, &
               'Negative values found - set to zero!', &
                RC, THISLOC = 'HCO_EmisAdd (HCO_FLUXARR_MOD.F90)' )

          ! Negative flag is 0: return w/ error
          ELSE
             CALL HCO_ERROR ( HcoState%Config%Err, &
               'Negative values found!', &
                RC, THISLOC = 'HCO_EmisAdd (HCO_FLUXARR_MOD.F90)' )
             RETURN
          ENDIF
       ENDIF
    ENDIF

    ! Eventually add universal scale factor
    CALL HCO_ScaleArr( HcoState, HcoID, Arr3D, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Add array
    HcoState%Spc(HcoID)%Emis%Val(:,:,:) = &
       HcoState%Spc(HcoID)%Emis%Val(:,:,:) + Arr3D(:,:,:)

    ! Check for diagnostics
    CALL DiagnCheck( HcoState,   ExtNr=ExtNr, Cat=Cat, &
                     Hier=Hier,  HcoID=HcoID, Arr3Dsp=Arr3D, &
                     MinDiagnLev=MinDiagnLev, RC=RC     )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Return w/ success
    RC = HCO_SUCCESS

  END SUBROUTINE HCO_EmisAdd_3D_Sp
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCO_EmisAdd_2D_Dp
!
! !DESCRIPTION: Routine HCO\_EmisAdd\_2D\_Dp adds the real*8 2D-array Arr2D
! to the emission array of species HcoID in HEMCO object HcoState.
! This routine also updates all autofill diagnostics that are defined
! for the givne species, extension number, emission category and
! hierarchy.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCO_EmisAdd_2D_Dp( HcoState,  Arr2D,    HcoID, &
                                RC,        ExtNr,    Cat,   Hier,  &
                                MinDiagnLev )
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(HCO_State), POINTER                 :: HcoState
    REAL(dp),        INTENT(INOUT)           :: Arr2D(HcoState%NX,HcoState%NY)
    INTEGER,         INTENT(INOUT)           :: RC
!
! !INPUT PARAMETERS:
!
    INTEGER,         INTENT(IN   )           :: HcoID
    INTEGER,         INTENT(IN   ), OPTIONAL :: ExtNr
    INTEGER,         INTENT(IN   ), OPTIONAL :: Cat
    INTEGER,         INTENT(IN   ), OPTIONAL :: Hier
    INTEGER,         INTENT(IN   ), OPTIONAL :: MinDiagnLev
!
! !REVISION HISTORY:
!  01 May 2013 - C. Keller - Initial version
!  20 Apr 2015 - C. Keller - Added DiagnCheck
!  12 May 2017 - C. Keller - Added option to use uniform scale factor
!EOP
!------------------------------------------------------------------------------
!BOC

    !=====================================================================
    ! HCO_EmisAdd_2D_Dp begins here!
    !=====================================================================

    ! Make sure target flux array in HcoState is allocated
    CALL HCO_ArrAssert ( HcoState%Spc(HcoID)%Emis, &
                         HcoState%NX, HcoState%NY, HcoState%NZ, RC )
    IF ( RC /=HCO_SUCCESS ) RETURN

    ! Check for negative values. NegFlag determines the behavior for
    ! negative values: 2 = allow; 1 = set to zero + warning; 0 = error.
    IF ( HcoState%Options%NegFlag /= 2 ) THEN
       IF ( ANY( Arr2D < 0.0_hp ) ) THEN

          ! Negative flag is 1: set to zero and prompt warning
          IF ( HcoState%Options%NegFlag == 1 ) THEN
             WHERE ( Arr2D < 0.0_hp ) Arr2D = 0.0_hp
             CALL HCO_WARNING ( HcoState%Config%Err, &
               'Negative values found - set to zero!', &
                RC, THISLOC = 'HCO_EmisAdd (HCO_FLUXARR_MOD.F90)' )

          ! Negative flag is 0: return w/ error
          ELSE
             CALL HCO_ERROR ( HcoState%Config%Err, &
               'Negative values found!', &
                RC, THISLOC = 'HCO_EmisAdd (HCO_FLUXARR_MOD.F90)' )
             RETURN
          ENDIF
       ENDIF
    ENDIF

    ! Eventually add universal scale factor
    CALL HCO_ScaleArr( HcoState, HcoID, Arr2D, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Add array
    HcoState%Spc(HcoID)%Emis%Val(:,:,1) = &
       HcoState%Spc(HcoID)%Emis%Val(:,:,1) + Arr2D(:,:)

    ! Check for diagnostics
    CALL DiagnCheck( HcoState,   ExtNr=ExtNr, Cat=Cat, &
                     Hier=Hier,  HcoID=HcoID, Arr2D=Arr2D, &
                     MinDiagnLev=MinDiagnLev, RC=RC     )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Return w/ success
    RC = HCO_SUCCESS

  END SUBROUTINE HCO_EmisAdd_2D_Dp
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCO_EmisAdd_2D_Sp
!
! !DESCRIPTION: Routine HCO\_EmisAdd\_2D\_Sp adds the real*4 2D-array Arr2D
! to the emission array of species HcoID in HEMCO object HcoState.
! This routine also updates all autofill diagnostics that are defined
! for the givne species, extension number, emission category and
! hierarchy.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCO_EmisAdd_2D_Sp( HcoState,  Arr2D,    HcoID, &
                                RC,        ExtNr,    Cat,   Hier,  &
                                MinDiagnLev )
!
! !INPUT/OUTPUT PARAMETERS:
!
!
    TYPE(HCO_State), POINTER                 :: HcoState
    REAL(sp),        INTENT(INOUT)           :: Arr2D(HcoState%NX,HcoState%NY)
    INTEGER,         INTENT(INOUT)           :: RC
!
! !INPUT PARAMETERS:
!
    INTEGER,         INTENT(IN   )           :: HcoID
    INTEGER,         INTENT(IN   ), OPTIONAL :: ExtNr
    INTEGER,         INTENT(IN   ), OPTIONAL :: Cat
    INTEGER,         INTENT(IN   ), OPTIONAL :: Hier
    INTEGER,         INTENT(IN   ), OPTIONAL :: MinDiagnLev
!
! !REVISION HISTORY:
!  01 May 2013 - C. Keller - Initial version
!  20 Apr 2015 - C. Keller - Added DiagnCheck
!  12 May 2017 - C. Keller - Added option to use uniform scale factor
!EOP
!------------------------------------------------------------------------------
!BOC

    !=====================================================================
    ! HCO_EmisAdd_2D_Sp begins here!
    !=====================================================================

    ! Make sure target flux array in HcoState is allocated
    CALL HCO_ArrAssert ( HcoState%Spc(HcoID)%Emis, &
                         HcoState%NX, HcoState%NY, HcoState%NZ, RC )
    IF ( RC /=HCO_SUCCESS ) RETURN

    ! Check for negative values. NegFlag determines the behavior for
    ! negative values: 2 = allow; 1 = set to zero + warning; 0 = error.
    IF ( HcoState%Options%NegFlag /= 2 ) THEN
       IF ( ANY( Arr2D < 0.0_sp ) ) THEN

          ! Negative flag is 1: set to zero and prompt warning
          IF ( HcoState%Options%NegFlag == 1 ) THEN
             WHERE ( Arr2D < 0.0_sp ) Arr2D = 0.0_sp
             CALL HCO_WARNING ( HcoState%Config%Err, &
               'Negative values found - set to zero!', &
                RC, THISLOC = 'HCO_EmisAdd (HCO_FLUXARR_MOD.F90)' )

          ! Negative flag is 0: return w/ error
          ELSE
             CALL HCO_ERROR ( HcoState%Config%Err, &
               'Negative values found!', &
                RC, THISLOC = 'HCO_EmisAdd (HCO_FLUXARR_MOD.F90)' )
             RETURN
          ENDIF
       ENDIF
    ENDIF

    ! Eventually add universal scale factor
    CALL HCO_ScaleArr( HcoState, HcoID, Arr2D, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Add array
    HcoState%Spc(HcoID)%Emis%Val(:,:,1) = &
       HcoState%Spc(HcoID)%Emis%Val(:,:,1) + Arr2D(:,:)

    ! Check for diagnostics
    CALL DiagnCheck( HcoState,   ExtNr=ExtNr, Cat=Cat, &
                     Hier=Hier,  HcoID=HcoID, Arr2Dsp=Arr2D, &
                     MinDiagnLev=MinDiagnLev, RC=RC     )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Return w/ success
    RC = HCO_SUCCESS

  END SUBROUTINE HCO_EmisAdd_2D_Sp
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCO_EmisAdd_Dp
!
! !DESCRIPTION: Routine HCO\_EmisAdd\_Dp adds value iVal to the emission
! array of species HcoID in HEMCO object HcoState. The value is placed at
! location I, J, L of the array.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCO_EmisAdd_Dp( HcoState, iVal, HcoID, I, J, L, RC )
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(HCO_State), POINTER       :: HcoState
    REAL(dp),        INTENT(INOUT) :: iVal
    INTEGER,         INTENT(INOUT) :: RC
!
! !INPUT PARAMETERS:
!
    INTEGER,         INTENT(IN   ) :: HcoID
    INTEGER,         INTENT(IN   ) :: I
    INTEGER,         INTENT(IN   ) :: J
    INTEGER,         INTENT(IN   ) :: L
!
! !REVISION HISTORY:
!  01 May 2013 - C. Keller - Initial version
!  12 May 2017 - C. Keller - Added option to use uniform scale factor
!EOP
!------------------------------------------------------------------------------
!BOC
    CHARACTER(LEN=255) :: MSG

    !=====================================================================
    ! HCO_EmisAdd_Dp begins here!
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

    ! Check for negative values. NegFlag determines the behavior for
    ! negative values: 2 = allow; 1 = set to zero + warning; 0 = error.
    IF ( HcoState%Options%NegFlag /= 2 ) THEN
       IF ( iVal < 0.0_hp ) THEN

          ! Negative flag is 1: set to zero and prompt warning
          IF ( HcoState%Options%NegFlag == 1 ) THEN
             iVal = 0.0_hp
             CALL HCO_WARNING ( HcoState%Config%Err, &
               'Negative values found - set to zero!', &
                RC, THISLOC = 'HCO_EmisAdd (HCO_FLUXARR_MOD.F90)' )

          ! Negative flag is 0: return w/ error
          ELSE
             CALL HCO_ERROR ( HcoState%Config%Err, &
               'Negative values found!', &
                RC, THISLOC = 'HCO_EmisAdd (HCO_FLUXARR_MOD.F90)' )
             RETURN
          ENDIF
       ENDIF
    ENDIF

    ! Eventually add universal scale factor
    CALL HCO_ScaleArr( HcoState, HcoID, iVal, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Add array
    HcoState%Spc(HcoID)%Emis%Val(I,J,L) = &
       HcoState%Spc(HcoID)%Emis%Val(I,J,L) + iVal

    ! Return w/ success
    RC = HCO_SUCCESS

  END SUBROUTINE HCO_EmisAdd_Dp
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCO_EmisAdd_Sp
!
! !DESCRIPTION: Routine HCO\_EmisAdd\_Sp adds value iVal to the emission
! array of species HcoID in HEMCO object HcoState. The value is placed
! at location I, J, L of the array.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCO_EmisAdd_Sp( HcoState, iVal, HcoID, I, J, L, RC )
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(HCO_State), POINTER       :: HcoState
    REAL(sp),        INTENT(INOUT) :: iVal
    INTEGER,         INTENT(INOUT) :: RC
!
! !INPUT PARAMETERS:
!
    INTEGER,         INTENT(IN   ) :: HcoID
    INTEGER,         INTENT(IN   ) :: I
    INTEGER,         INTENT(IN   ) :: J
    INTEGER,         INTENT(IN   ) :: L
!
! !REVISION HISTORY:
!  01 May 2013 - C. Keller - Initial version
!  12 May 2017 - C. Keller - Added option to use uniform scale factor
!EOP
!------------------------------------------------------------------------------
!BOC
    CHARACTER(LEN=255) :: MSG

    !=====================================================================
    ! HCO_EmisAdd_Sp begins here!
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

    ! Check for negative values. NegFlag determines the behavior for
    ! negative values: 2 = allow; 1 = set to zero + warning; 0 = error.
    IF ( HcoState%Options%NegFlag /= 2 ) THEN
       IF ( iVal < 0.0_sp ) THEN

          ! Negative flag is 1: set to zero and prompt warning
          IF ( HcoState%Options%NegFlag == 1 ) THEN
             iVal = 0.0_sp
             CALL HCO_WARNING ( HcoState%Config%Err, &
               'Negative values found - set to zero!', &
                RC, THISLOC = 'HCO_EmisAdd (HCO_FLUXARR_MOD.F90)' )

          ! Negative flag is 0: return w/ error
          ELSE
             CALL HCO_ERROR ( HcoState%Config%Err, &
               'Negative values found!', &
                RC, THISLOC = 'HCO_EmisAdd (HCO_FLUXARR_MOD.F90)' )
             RETURN
          ENDIF
       ENDIF
    ENDIF

    ! Eventually add universal scale factor
    CALL HCO_ScaleArr( HcoState, HcoID, iVal, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Add array
    HcoState%Spc(HcoID)%Emis%Val(I,J,L) = &
       HcoState%Spc(HcoID)%Emis%Val(I,J,L) + iVal

    ! Return w/ success
    RC = HCO_SUCCESS

  END SUBROUTINE HCO_EmisAdd_Sp
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCO_DepvAdd_2D_Dp
!
! !DESCRIPTION: Routine HCO\_DepvAdd\_2D\_Dp adds the real*8 2D-array Arr2D
! to the depostion array of species HcoID in HEMCO object HcoState.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCO_DepvAdd_2D_Dp( HcoState, Arr2D, HcoID, RC )
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(HCO_State), POINTER       :: HcoState
    INTEGER,         INTENT(INOUT) :: RC
!
! !INPUT PARAMETERS:
!
    REAL(dp),        INTENT(IN   ) :: Arr2D(HcoState%NX,HcoState%NY)
    INTEGER,         INTENT(IN   ) :: HcoID
!
! !REVISION HISTORY:
!  01 May 2013 - C. Keller - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
    CHARACTER(LEN=255) :: MSG

    !=====================================================================
    ! HCO_DepvAdd_2D_Dp begins here!
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

  END SUBROUTINE HCO_DepvAdd_2D_Dp
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCO_DepvAdd_2D_Sp
!
! !DESCRIPTION: Routine HCO\_DepvAdd\_2D\_Sp adds the real*4 2D-array Arr2D
! to the depostion array of species HcoID in HEMCO object HcoState.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCO_DepvAdd_2D_Sp( HcoState, Arr2D, HcoID, RC )
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(HCO_State), POINTER       :: HcoState
    INTEGER,         INTENT(INOUT) :: RC
!
! !INPUT PARAMETERS:
!
    REAL(sp),        INTENT(IN   ) :: Arr2D(HcoState%NX,HcoState%NY)
    INTEGER,         INTENT(IN   ) :: HcoID
!
! !REVISION HISTORY:
!  01 May 2013 - C. Keller - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC

    !=====================================================================
    ! HCO_DepvAdd_2D_Sp begins here!
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

  END SUBROUTINE HCO_DepvAdd_2D_Sp
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCO_DepvAdd_Dp
!
! !DESCRIPTION: Routine HCO\_DepvAdd\_Dp adds value iVal to the deposition
! array of species HcoID in HEMCO object HcoState. The value is placed at
! location I, J of the array.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCO_DepvAdd_Dp( HcoState, iVal, HcoID, I, J, RC )
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(HCO_State), POINTER       :: HcoState
    INTEGER,         INTENT(INOUT) :: RC
!
! !INPUT PARAMETERS:
!
    REAL(dp),        INTENT(IN   ) :: iVal
    INTEGER,         INTENT(IN   ) :: HcoID
    INTEGER,         INTENT(IN   ) :: I
    INTEGER,         INTENT(IN   ) :: J
!
! !REVISION HISTORY:
!  01 May 2013 - C. Keller - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
    CHARACTER(LEN=255) :: MSG

    !=====================================================================
    ! HCO_DepvAdd_Dp begins here!
    !=====================================================================

    ! Check size dimensions
    IF ( I > HcoState%NX ) THEN
       MSG = 'Cannot add DP - i too high!'
       CALL HCO_ERROR ( HcoState%Config%Err, MSG, RC )
       RETURN
    ENDIF
    IF ( J > HcoState%NY ) THEN
       MSG = 'Cannot add DP - j too high!'
       CALL HCO_ERROR ( HcoState%Config%Err, MSG, RC )
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

  END SUBROUTINE HCO_DepvAdd_Dp
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCO_DepvAdd_Sp
!
! !DESCRIPTION: Routine HCO\_DepvAdd\_Sp adds value iVal to the deposition
! array of species HcoID in HEMCO object HcoState. The value is placed at
! location I, J of the array.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCO_DepvAdd_Sp( HcoState, iVal, HcoID, I, J, RC )
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(HCO_State), POINTER       :: HcoState
    INTEGER,         INTENT(INOUT) :: RC
!
! !INPUT PARAMETERS:
!
    REAL(sp),        INTENT(IN   ) :: iVal
    INTEGER,         INTENT(IN   ) :: HcoID
    INTEGER,         INTENT(IN   ) :: I
    INTEGER,         INTENT(IN   ) :: J
!
! !REVISION HISTORY:
!  01 May 2013 - C. Keller - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
    CHARACTER(LEN=255) :: MSG

    !=====================================================================
    ! HCO_DepvAdd_Sp begins here!
    !=====================================================================

    ! Check size dimensions
    IF ( I > HcoState%NX ) THEN
       MSG = 'Cannot add iVal - i too high!'
       CALL HCO_ERROR ( HcoState%Config%Err, MSG, RC )
       RETURN
    ENDIF
    IF ( J > HcoState%NY ) THEN
       MSG = 'Cannot add iVal - j too high!'
       CALL HCO_ERROR ( HcoState%Config%Err, MSG, RC )
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

  END SUBROUTINE HCO_DepvAdd_Sp
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: DiagnCheck
!
! !DESCRIPTION: Subroutine DiagnCheck checks if the given emission array needs
!               to be added to any auto-fill diagnostics. The diagnostics to be
!               filled (if any) depend on the passed extension number, emission
!               category and hierarchy, and the HEMCO species ID.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE DiagnCheck( HcoState,  ExtNr,    Cat,     &
                         Hier,      HcoID,    Arr3D,   Arr3Dsp, &
                         Arr2D,     Arr2Dsp,  MinDiagnLev,  RC   )
!
! !USES:
!
    USE HCO_DIAGN_MOD
!
! !INPUT PARAMETERS:
!
    INTEGER,         INTENT(IN   )           :: HcoID
    INTEGER,         INTENT(IN   ), OPTIONAL :: ExtNr
    INTEGER,         INTENT(IN   ), OPTIONAL :: Cat
    INTEGER,         INTENT(IN   ), OPTIONAL :: Hier
    INTEGER,         INTENT(IN   ), OPTIONAL :: MinDiagnLev
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(HCO_State), POINTER                 :: HcoState
    REAL(dp),        INTENT(INOUT), OPTIONAL :: Arr3D(   HcoState%NX, &
                                                         HcoState%NY, &
                                                         HcoState%NZ )
    REAL(sp),        INTENT(INOUT), OPTIONAL :: Arr3Dsp( HcoState%NX, &
                                                         HcoState%NY, &
                                                         HcoState%NZ )
    REAL(dp),        INTENT(INOUT), OPTIONAL :: Arr2D(   HcoState%NX, &
                                                         HcoState%NY )
    REAL(sp),        INTENT(INOUT), OPTIONAL :: Arr2Dsp( HcoState%NX, &
                                                         HcoState%NY )
    INTEGER,         INTENT(INOUT)           :: RC
!
! !REVISION HISTORY:
!  20 Apr 2015 - C. Keller - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER :: AFL, XT, CT, HR

    !=====================================================================
    ! DiagnCheck begins here!
    !=====================================================================

    ! Initialize values

    ! Autofill level:
    ! 1=species level, 2=ExtNr level, 3=Cat level, 4=Hier level
    AFL =  1

    ! ExtNr, Cat & Hier
    XT  = -1
    CT  = -1
    HR  = -1

    ! Set extension number, category, and hierarchy
    IF ( PRESENT(ExtNr) ) THEN
       XT  = ExtNr
       AFL = 2

       ! Consider category only within extension ...
       IF ( PRESENT(Cat) ) THEN
          IF ( Cat > 0 ) THEN
             CT  = Cat
             AFL = 3
          ENDIF

          ! Consider hierarchy only within category ...
          IF ( AFL==3 .AND. PRESENT(Hier) ) THEN
             IF ( Hier > 0 ) THEN
                HR  = Hier
                AFL = 4
             ENDIF
          ENDIF
       ENDIF
    ENDIF

    IF ( PRESENT(MinDiagnLev) ) THEN
       AFL = MIN(AFL,MinDiagnLev)
    ENDIF

    ! Check if we need to call diagnostics
    IF ( Diagn_AutoFillLevelDefined(HcoState%Diagn,AFL) ) THEN

       ! 3D HP array
       IF ( PRESENT(Arr3D) ) THEN
          CALL Diagn_Update( HcoState, ExtNr=XT,    Cat=CT,     &
                             Hier=HR,  HcoID=HcoID, AutoFill=1, &
                             Array3D=Arr3D, MinDiagnLev=MinDiagnLev, RC=RC )
          IF ( RC /= HCO_SUCCESS ) RETURN
       ENDIF

       ! 3D SP array
       IF ( PRESENT(Arr3Dsp) ) THEN
          CALL Diagn_Update( HcoState, ExtNr=XT,    Cat=CT,     &
                             Hier=HR,  HcoID=HcoID, AutoFill=1, &
                             Array3D=Arr3Dsp, MinDiagnLev=MinDiagnLev, RC=RC )
          IF ( RC /= HCO_SUCCESS ) RETURN
       ENDIF

       ! 2D HP array
       IF ( PRESENT(Arr2D) ) THEN
          CALL Diagn_Update( HcoState, ExtNr=XT,    Cat=CT,     &
                             Hier=HR,  HcoID=HcoID, AutoFill=1, &
                             Array2D=Arr2D, MinDiagnLev=MinDiagnLev, RC=RC )
          IF ( RC /= HCO_SUCCESS ) RETURN
       ENDIF

       ! 2D SP array
       IF ( PRESENT(Arr2Dsp) ) THEN
          CALL Diagn_Update( HcoState, ExtNr=XT,    Cat=CT,     &
                             Hier=HR,  HcoID=HcoID, AutoFill=1, &
                             Array2D=Arr2Dsp, MinDiagnLev=MinDiagnLev, RC=RC )
          IF ( RC /= HCO_SUCCESS ) RETURN
       ENDIF

    ENDIF

    ! Return w/ success
    RC = HCO_SUCCESS

  END SUBROUTINE DiagnCheck
!EOC
END MODULE HCO_FluxArr_Mod
