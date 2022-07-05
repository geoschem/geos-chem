#include "MAPL_Generic.h"
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Model                            !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: geos_aerocupler
!
! !DESCRIPTION: Module with routines and variables to couple GEOS-Chem aerosols 
!  with GEOS 
!\\
!\\
! !INTERFACE:
!
MODULE GEOS_AeroCoupler
!
! !USES:
!
  ! MAPL/ESMF
  USE ESMF     
  USE MAPL_Mod 
  ! GEOS-Chem
  USE Precision_Mod
  USE ErrCode_Mod                                    ! Error numbers
  USE Input_Opt_Mod,         ONLY : OptInput
  USE State_Chm_Mod,         ONLY : ChmState, Ind_   ! Chemistry State obj
  USE State_Met_Mod,         ONLY : MetState         ! Meteorology State obj
  USE State_Diag_Mod,        ONLY : DgnState         ! Diagnostics State obj
  USE State_Grid_Mod,        ONLY : GrdState         ! Grid State obj

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC   :: GEOS_FillAeroBundle 
!
! !PRIVATE MEMBER FUNCTIONS:
!
  PRIVATE  :: FillAeroDP_
!
! !PRIVATE TYPES:
!
! !REVISION HISTORY:
!  05 Jul 2022 - C. Keller - initial version (refactored Chem_GridCompMod)
!  See https://github.com/geoschem/geos-chem for full history
!EOP
!------------------------------------------------------------------------------
!BOC
CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Model                            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: GEOS_FillAeroBundle
!
! !DESCRIPTION: Routine to fill the aerosol bundle. 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE GEOS_FillAeroBundle( GC, EXPORT, State_Chm, State_Grid, Input_Opt, RC ) 
!
! !USES:
!
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ESMF_GridComp), INTENT(INOUT), TARGET :: GC     ! Ref to this GridComp
    TYPE(ESMF_State),    INTENT(INOUT)         :: Export   ! Export State
    TYPE(ChmState),      INTENT(INOUT)         :: State_Chm
    TYPE(GrdState),      INTENT(INOUT)         :: State_Grid
    TYPE(OptInput),      INTENT(INOUT)         :: Input_Opt
!
! !OUTPUT PARAMETERS:
!
    INTEGER,             INTENT(INOUT)         :: RC       ! Success or failure?
!
! !REMARKS:
!
! !REVISION HISTORY:
!  05 Jul 2022 - C. Keller   - Initial version (refactored Chem_GridCompMod)
!  See https://github.com/geoschem/geos-chem for history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! LOCAL VARIABLES:
!
    LOGICAL                      :: am_I_Root
    REAL                         :: GCMW, FRAC
    INTEGER                      :: N, IndSpc
    INTEGER                      :: nAero, nLen, GCID
    TYPE(ESMF_STATE)             :: Aero
    TYPE(ESMF_FieldBundle)       :: AeroBdl
    TYPE(ESMF_Field)             :: AeroFld
    CHARACTER(LEN=ESMF_MAXSTR)   :: GCName, AeroName
    REAL(fp), POINTER            :: GcPtr3d  (:,:,:) => NULL()
    REAL, POINTER                :: AeroPtr3d(:,:,:) => NULL()

    INTEGER, PARAMETER           :: NRATS = 5
    CHARACTER(LEN=15), PARAMETER :: RatsNames(NRATS) = (/ 'CH4', 'N2O', 'CFC11', 'CFC12', 'HCFC22' /)

    __Iam__('GEOS_FillAeroBundle')

    !=======================================================================
    ! GEOS_FillAeroBundle starts here
    !=======================================================================

    ! Are we on the root PET?
    am_I_Root = MAPL_Am_I_Root()

    ! For every field of the AERO bundle, we will copy the corresponding
    ! GEOS-Chem tracer field, converting units from mol mol-1 to kg kg-1.

    ! Get AERO bundle
    CALL ESMF_StateGet( EXPORT, 'AERO',     Aero,    __RC__ )
    CALL ESMF_StateGet( Aero,   'AEROSOLS', AeroBdl, __RC__ )

    ! Number of fields in the AERO Bundle
    CALL ESMF_FieldBundleGet ( AeroBdl, FieldCount=nAero, __RC__ )

    ! Update every field
    DO N = 1, nAero

       ! Get field
       CALL ESMF_FieldBundleGet( AeroBdl, N, AeroFld, __RC__ )

       ! Extract GC tracer name, molecular weight and fraction to be used
       CALL ESMF_AttributeGet( AeroFld, NAME='GCNAME', VALUE=GcName, __RC__ )
       CALL ESMF_AttributeGet( AeroFld, NAME='GCMW'  , VALUE=GCMW,   __RC__ )
       CALL ESMF_AttributeGet( AeroFld, NAME='FRAC',   VALUE=FRAC,  __RC__ )

       ! Get pointer to Aero data
       CALL ESMF_FieldGet( AeroFld, farrayPtr=AeroPtr3D, __RC__ )

       ! Get pointer to GC data
       nlen = LEN(TRIM(GcName))
       IndSpc = Ind_(TRIM(GcName(5:nlen)))
       ASSERT_(IndSpc>0)
       GcPtr3D => State_Chm%Species(IndSpc)%Conc(:,:,State_Grid%NZ:1:-1)
       !CALL MAPL_GetPointer ( INTSTATE, GcPtr3D, TRIM(GcName), __RC__ )

       ! Pass GC to AERO. Convert from mol/mol to kg/kg. Only use the
       ! fraction specified during initialization (different from 1 for
       ! sea salt aerosols only)
       !AeroPtr3D = GcPtr3D * FRAC * GCMW / MAPL_AIRMW
       AeroPtr3D = GcPtr3D * FRAC

       !!! writing to diagnostics
       GcPtr3D   => NULL()
       CALL ESMF_FieldGet( AeroFld, NAME=GcName, __RC__ )
       CALL MAPL_GetPointer ( EXPORT, GcPtr3D, 'AERO_'//TRIM(GcName), &
                              NotFoundOk=.TRUE., __RC__ )
       IF ( ASSOCIATED(GcPtr3D) ) GcPtr3D = AeroPtr3D

       ! Free pointers
       GcPtr3D   => NULL()
       AeroPtr3D => NULL()
    ENDDO

    ! Fill AERO_DP bundle
    CALL FillAeroDP_ ( am_I_Root, GC, EXPORT, Input_Opt, __RC__ )

    _RETURN(ESMF_SUCCESS)

    END SUBROUTINE GEOS_FillAeroBundle
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Model                            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: FillAeroDP_
!
! !DESCRIPTION: FillAeroDP_ fills the AERO_DP bundle
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE FillAeroDP_ ( am_I_Root, GC, EXPORT, Input_Opt, RC )
!
! !USES:
!
    USE HCO_ERROR_MOD
    USE HCO_TYPES_MOD,     ONLY : DiagnCont
    USE HCO_DIAGN_MOD,     ONLY : Diagn_Get
    USE HCO_State_GC_Mod,  ONLY : HcoState
!
! !INPUT PARAMETERS:
!
    LOGICAL                            :: am_I_Root
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ESMF_GridComp), INTENT(INOUT) :: GC       ! Ref to this GridComp
    TYPE(ESMF_State),    INTENT(INOUT) :: Export   ! Export State
    TYPE(OptInput),      INTENT(INOUT) :: Input_Opt
!
! !OUTPUT PARAMETERS:
!
    INTEGER, INTENT(OUT), OPTIONAL     :: RC
!
! !REVISION HISTORY:
!  30 Mar 2015 - C. Keller   - Initial version
!  See https://github.com/geoschem/geos-chem for history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! LOCAL VARIABLES:
!

    REAL, POINTER                :: Ptr2d(:,:) => NULL()
    INTEGER                      :: I, J, N, TrcID
    CHARACTER(LEN= 2)            :: Prfx
    CHARACTER(LEN=15)            :: TrcName
    CHARACTER(LEN=ESMF_MAXSTR)   :: ExpName

    ! Hemco diagnostics
    INTEGER                      :: DgnID
    INTEGER                      :: FLAG, ERR
    TYPE(DiagnCont), POINTER     :: DgnCont => NULL()

    ! Error handling
    INTEGER                      :: STATUS
    CHARACTER(LEN=ESMF_MAXSTR)   :: Iam

    !=======================================================================
    ! FillAeroDP_ begins here
    !=======================================================================

    ! Traceback handle
    Iam = 'FillAeroDP_'

    ! There are 8 species in total
    DO N = 1, 8

       ! Get species ID
       SELECT CASE ( N )
          CASE ( 1 )
             TrcName = 'DST1'
             Prfx    = 'DU'
          CASE ( 2 )
             TrcName = 'DST2'
             Prfx    = 'DU'
          CASE ( 3 )
             TrcName = 'DST3'
             Prfx    = 'DU'
          CASE ( 4 )
             TrcName = 'DST4'
             Prfx    = 'DU'
          CASE ( 5 )
             TrcName = 'BCPI'
             Prfx    = 'BC'
          CASE ( 6 )
             TrcName = 'BCPO'
             Prfx    = 'BC'
          CASE ( 7 )
             TrcName = 'OCPI'
             Prfx    = 'OC'
          CASE ( 8 )
             TrcName = 'OCPO'
             Prfx    = 'OC'
          CASE DEFAULT
             TrcName = 'YeahYeahYeah'
       END SELECT

       ! Get GEOS-Chem tracer ID
       TrcID = Ind_( TRIM(TrcName) )

       ! Only if tracer is defined...
       IF ( TrcID <= 0 ) CYCLE

       ! Dry dep and wet dep
       DO I = 1, 2

          IF ( I == 1 ) THEN
             ExpName = TRIM(Prfx)//'DP_'//TRIM(TrcName)
          ELSEIF ( I == 2 ) THEN
             ExpName = TRIM(Prfx)//'WT_'//TRIM(TrcName)
          ENDIF

          ! Get pointer
          CALL MAPL_GetPointer( EXPORT, Ptr2D, TRIM(ExpName),   &
                                notFoundOk=.TRUE., __RC__ )

          ! Skip if not defined
          IF ( .NOT. ASSOCIATED(Ptr2D) ) CYCLE

          ! Reset
          Ptr2D = 0.0

          ! ------------------
          ! Dry deposition
          ! ------------------
          IF ( I == 1 ) THEN

             ! Get diagnostics
             DgnID = 44500 + TrcID
             CALL Diagn_Get( HcoState, .FALSE., DgnCont,  &
                             FLAG, ERR, cID=DgnID, AutoFill=-1,      &
                             COL=Input_Opt%DIAG_COLLECTION )

             ! Error check
             _ASSERT( ERR == HCO_SUCCESS,'Error calling Diagn_Get' )

             ! Add to array if diagnostics is defined
             ! GEOS-Chem diagnostics is in kg m-2 s-1.
             IF ( FLAG == HCO_SUCCESS ) THEN
                IF ( ASSOCIATED(DgnCont%Arr2D%Val) ) THEN
                   Ptr2D = Ptr2D + DgnCont%Arr2D%Val
                ENDIF
             ENDIF

          ! ------------------
          ! Wet depostion
          ! ------------------
          ELSEIF ( I == 2 ) THEN

             ! Convective and wet scavenging
             DO J = 1, 2

                SELECT CASE ( J )
                   ! Convection:
                   CASE ( 1 )
                      DgnID = 38000 + TrcID
                   ! Wet deposition
                   CASE ( 2 )
                      DgnID = 39000 + TrcID
                   CASE DEFAULT
                      DgnID = -1
                END SELECT

                ! Get diagnostics
                CALL Diagn_Get( HcoState, .FALSE., DgnCont,  &
                                FLAG, ERR, cID=DgnID, AutoFill=-1,      &
                                COL=Input_Opt%DIAG_COLLECTION )

                ! Error check
                _ASSERT( ERR == HCO_SUCCESS,'Error calling Diagn_Get' )

                ! Add to array if diagnostics is defined. GEOS-Chem
                ! diagnostics is already in kg m-2 s-1.
                IF ( FLAG == HCO_SUCCESS ) THEN
                   IF ( ASSOCIATED(DgnCont%Arr2D%Val) ) THEN
                      Ptr2D = Ptr2D + DgnCont%Arr2D%Val
                   ELSEIF ( ASSOCIATED(DgnCont%Arr3D%Val) ) THEN
                      Ptr2D = Ptr2D + SUM(DgnCont%Arr3D%Val,DIM=3)
                   ENDIF
                ENDIF
             ENDDO !J
          ENDIF

       ENDDO !I
    ENDDO !N

    ! Successful return
    RC = ESMF_SUCCESS

  END SUBROUTINE FillAeroDP_
!EOC
END MODULE GEOS_AeroCoupler
