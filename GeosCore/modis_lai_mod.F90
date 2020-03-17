!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: modis_lai_mod.F90
!
! !DESCRIPTION: Module MODIS\_LAI\_MOD reads the MODIS LAI data at
!  native resolution and then regrids them to the GEOS-Chem resolution on the
!  fly.
!
! !INTERFACE:
!
MODULE Modis_Lai_Mod
!
! !USES:
!
  USE Precision_Mod

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC :: Get_XlaiNative_from_Hemco
  PUBLIC :: Compute_Xlai
!
! !REMARKS:
!  (1) MODIS LAI data resolution is the same as the Olson land map. The Olson
!      2001 landmap is the default, therefore, you will use MODIS LAI data at
!      0.25 x 0.25 resolution.
!  (2) In HEMCO, MEGAN uses 'offline' MODIS LAI (State_Met%MODISLAI) which
!      is computed in this module.
!      in an ESMF environment, in which case State_Met%LAI is used instead.
!  (3) MODIS LAI arrays and where they are used in GEOS-Chem:
!       (a) State_Met%XLAI     --> dry deposition routine DEPVEL
!       (b) State_Met%MODISLAI --> MEGAN (if using standard GC); Hg(0)
!                                  emissions; several diagnostics
!
! !REVISION HISTORY:
!  03 Apr 2012 - R. Yantosca - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_xlainative_from_hemco
!
! !DESCRIPTION: Copies the MODIS XLAI data from HEMCO pointers into
!  the State\_Met object.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Get_XlaiNative_from_HEMCO( Input_Opt, State_Met, RC )
!
! !USES:
!
    USE CMN_SIZE_Mod,      ONLY : NSURFTYPE
    USE ErrCode_Mod
    USE HCO_EmisList_Mod,  ONLY : Hco_GetPtr
    USE HCO_Interface_Mod, ONLY : HcoState
    USE Input_Opt_Mod,     ONLY : OptInput
    USE State_Met_Mod,     ONLY : MetState
    USE Time_Mod,          ONLY : ITS_A_NEW_DAY
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(MetState), INTENT(INOUT) :: State_Met   ! Meteorology State object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure
!
! !REVISION HISTORY:
!  14 Feb 2019 - R. Yantosca - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER            :: T

    ! Strings
    CHARACTER(LEN=6)   :: Name
    CHARACTER(LEN=255) :: ErrMsg, ThisLoc

    ! Pointers
    REAL(f4), POINTER  :: Ptr2D(:,:)

    !========================================================================
    ! Get_XlaiNative_from_HEMCO begins here!
    !========================================================================

    ! Initialize
    RC      = GC_SUCCESS
    ErrMsg  = ''
    ThisLoc = &
      ' -> at Get_XLAI_From_HEMCO (in module GeosCore/modis_lai_mod.F90)'

    ! Free pointer
    Ptr2d => NULL()

    ! Loop over the # of Olson land types
    DO T = 1, NSURFTYPE

       ! Get the HEMCO pointer to the LAI for each Olson type
       ! (variable names are XLAI00, XLAI01, .. XLAI72)
       WRITE( Name, 100 ) T-1
 100   FORMAT( 'XLAI' , i2.2 )
       CALL HCO_GetPtr( HcoState, Name, Ptr2D, RC )

       ! Trap potential errors
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Could not get pointer to HEMCO field: ' // TRIM( Name )
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

       ! Copy into State_Met%XLAI_NATIVE
       State_Met%XLAI_NATIVE(:,:,T) = Ptr2D

       ! Free pointer
       Ptr2D => NULL()

    ENDDO

  END SUBROUTINE Get_XlaiNative_from_HEMCO
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: compute_xlai
!
! !DESCRIPTION: Subroutine COMPUTE\_XLAI computes MODIS-based leaf
!  area indices (LAI) per land type and grid cell. This computation uses
!  offline 0.25x0.25 MODIS LAI and Olson landmap data regridded to
!  the cubed sphere. Variables set include State\_Met%XLAI.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Compute_Xlai( Input_Opt, State_Grid, State_Met, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE Input_Opt_Mod,  ONLY : OptInput
    USE State_Grid_Mod, ONLY : GrdState
    USE State_Met_Mod,  ONLY : MetState
    USE Time_Mod,       ONLY : Its_A_New_Day
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
    TYPE(GrdState), INTENT(IN)    :: State_Grid  ! Grid State object
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(MetState), INTENT(INOUT) :: State_Met   ! Meteorology State object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
!
! !REMARKS:
!  State_Met%XLAI_NATIVE is the LAI data as it comes in from either HEMCO
!  or the MAPL import state via ExtData.
!                                                                             .
!  State_Met%XLAI is used for inputs into the GEOS-Chem dry deposition code.
!  It is the LAI binned into the 11 dry-deposition land types.
!                                                                             .
!  State_Met%MODISLAI is the average LAI per grid box, averaged over all
!  land types.  This is needed for the HEMCO soil NOx extension.
!
! !REVISION HISTORY:
!  18 Oct 2016 - E. Lundgren - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER  :: I, J, S, T
    REAL(fp) :: landFrac

    !======================================================================
    ! Initialize
    !======================================================================

    ! Initialize
    RC                 = GC_SUCCESS
    State_Met%XLAI     = 0.0_fp
    State_Met%MODISLAI = 0.0_fp

    ! Loop over all grid cells
    DO J = 1, State_Grid%NY
    DO I = 1, State_Grid%NX

       ! Loop over all surface types present in this grid cell
       DO S = 1, State_Met%IREG(I,J)

          ! Set current surface type index
          T = State_Met%ILAND(I,J,S) + 1

          ! Get fraction of cell with this surface type by retrieving it
          ! from the land type fraction calculated by ExtData
          landFrac = State_Met%LandTypeFrac(I,J,T)

          ! Set XLAI to average LAI for this surface type ( as calculated
          ! by ExtData using zeros for coverage by other surface types )
          ! divided by the fractional coverage of this surface type.
          ! The resultant XLAI is the average LAI for only the area
          ! with the current surface type, and therefore is larger than
          ! XLAI_NATIVE when other surface types exist within the cell.
          !
          ! NOTE: Unlike XLAI_NATIVE and LandTypeFrac, the 3rd
          ! dimension indexes of XLAI are NOT surface types 1-73! Instead,
          ! It is the surface indexes ILAND-1 and therefore contains
          ! zeros beyond the number of surface types present in the cell
          ! (IREG). This is for backwards compatibility with GC classic
          ! legacy drydep code.
          IF ( landFrac .gt. 1.e-9_fp ) THEN
             State_Met%XLAI(I,J,S) = State_Met%XLAI_NATIVE(I,J,T) / landFrac
          ENDIF

       ENDDO

       ! Calculate average LAI for this grid cell across all land types
       State_Met%MODISLAI(I,J) = SUM( State_Met%XLAI_NATIVE(I,J,:) )

    ENDDO
    ENDDO

  END SUBROUTINE Compute_XLAI
!EOC
END MODULE Modis_Lai_Mod
