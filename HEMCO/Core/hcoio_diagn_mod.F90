!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: hcoio_diagn_mod.F90
!
! !DESCRIPTION: Module HCOIO\_Diagn\_Mod.F90 is the data interface module
! for the HEMCO diagnostics. It contains routines to write out diagnostics
! into a netCDF file.
!\\
!\\
! In an ESMF/MAPL environment, the HEMCO diagnostics are not directly 
! written to disk but passed to the gridded component export state, where 
! they can be picked up by the MAPL HISTORY component.
!\\
!\\
! !INTERFACE:
!
MODULE HCOIO_DIAGN_MOD
!
! !USES:
!
  USE HCO_ERROR_MOD
  USE HCO_DIAGN_MOD

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC :: HcoDiagn_Write 
  PUBLIC :: HCOIO_DIAGN_WRITEOUT
!
! !REMARKS:
!  HEMCO diagnostics are still in testing mode. We will fully activate them
!  at a later time.  They will be turned on when debugging & unit testing.
!
! !REVISION HISTORY:
!  04 May 2014 - C. Keller   - Initial version. 
!  11 Jun 2014 - R. Yantosca - Cosmetic changes in ProTeX headers
!  11 Jun 2014 - R. Yantosca - Now use F90 freeform indentation
!  28 Jul 2014 - C. Keller   - Removed GC specific initialization calls and
!                              moved to HEMCO core.
!  05 Aug 2014 - C. Keller   - Added dummy interface for ESMF.
!  03 Apr 2015 - C. Keller   - Added HcoDiagn_Write
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
  ! Fill value used in HEMCO diagnostics netCDF files.
!  REAL(hp), PARAMETER :: FillValue = 1.e-31_hp
  REAL(sp), PARAMETER :: FillValue = HCO_MISSVAL 

CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HcoDiagn_Write
!
! !DESCRIPTION: Subroutine HcoDiagn_Write is the wrapper routine to write out
! the content of the built-in HEMCO diagnostics. If input argument Restart is
! set to TRUE, only the restart collection will be written out. Otherwise,
! the default collection
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HcoDiagn_Write( am_I_Root, HcoState, Restart, RC )
!
! !USES:
!
    USE HCO_State_Mod,       ONLY : HCO_State
!
! !INPUT/OUTPUT PARAMETERS:
!
    LOGICAL,         INTENT(IN   )    :: am_I_Root    ! Root CPU?
    TYPE(HCO_State), POINTER          :: HcoState     ! HEMCO state object 
    LOGICAL,         INTENT(IN   )    :: Restart      ! write restart (enforced)?
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,         INTENT(INOUT)    :: RC           ! Return code
!
! !REVISION HISTORY: 
!  03 Apr 2015 - C. Keller   - Initial version 
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER            :: I, COL
    CHARACTER(LEN=255) :: MSG, LOC

    !=================================================================
    ! HcoDiagn_Write begins here!
    !=================================================================

    ! Init
    LOC = 'HcoDiagn_Write (hcoi_diagn_mod.F90)'

    ! To write restart (enforced)
    IF ( RESTART ) THEN
       CALL HCOIO_DIAGN_WRITEOUT ( am_I_Root,                       &
                                   HcoState,                        &
                                   ForceWrite  = .TRUE.,            &
                                   UsePrevTime = .FALSE.,           &
                                   COL         = HcoDiagnIDRestart, &
                                   RC          = RC                  )
       IF( RC /= HCO_SUCCESS) RETURN 

    ! Write all HEMCO diagnostics
    ELSE

       ! Loop over all collections that shall be written out.
       ! HCOIO_DIAGN_WRITEOUT will determine whether it is time to
       ! write a collection or not.
       DO I = 1, 3 

          ! Define collection ID
          SELECT CASE ( I ) 
             CASE ( 1 ) 
                COL = HcoDiagnIDDefault
             CASE ( 2 ) 
                COL = HcoDiagnIDRestart
             CASE ( 3 ) 
                COL = HcoDiagnIDManual
          END SELECT
  
          ! If not ESMF environment, never write the manual diagnostics
          ! to disk. Instead, the content of the manual diagnostics needs
          ! to be fetched explicitly.
#if       !defined ( ESMF_ ) 
          IF ( I == 3 ) CYCLE
#endif
 
          ! Restart file 
          CALL HCOIO_DIAGN_WRITEOUT ( am_I_Root,                       &
                                      HcoState,                        &
                                      ForceWrite  = .FALSE.,           &
                                      UsePrevTime = .FALSE.,           &
                                      COL         = COL,               &
                                      RC          = RC                  )
          IF(RC /= HCO_SUCCESS) RETURN 
       ENDDO
    ENDIF

  END SUBROUTINE HcoDiagn_Write 
!EOC
# if !defined(ESMF_)
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCOIO_Diagn_WriteOut
!
! !DESCRIPTION: Subroutine HCOIO\_Diagn\_WriteOut writes diagnostics to 
! netCDF file. If the ForceWrite flag is set to TRUE, all diagnostics are
! written out except they have already been written out during this time
! step. This option is usually only used at the end of a simulation run.
! If ForceWrite is False, only the diagnostics that are at the end of their
! time averaging interval are written. For example, if the current month
! is different from the previous (emissions) month, all diagnostics with 
! hourly, daily and monthly time averaging intervals are written out.
! If the optional argument OnlyIfFirst is set to TRUE, diagnostics will
! only be written out if its nnGetCalls is 1. This can be used to avoid
! that diagnostics will be written out twice. The nnGetCalls is reset to
! zero the first time a diagnostics is updated. For diagnostics that
! point to data stored somewhere else (i.e. that simply contain a data
! pointer, nnGetCalls is never reset and keeps counting.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCOIO_Diagn_WriteOut( am_I_Root,   HcoState, ForceWrite,  &
                                   RC,          PREFIX,   UsePrevTime, &
                                   OnlyIfFirst, COL                     )
!
! !USES:
!
    USE m_netCDF_io_define
    USE Ncdf_Mod,            ONLY : NC_Create
    USE Ncdf_Mod,            ONLY : NC_Close
    USE Ncdf_Mod,            ONLY : NC_Var_Def
    USE Ncdf_Mod,            ONLY : NC_Var_Write
    USE HCO_State_Mod,       ONLY : HCO_State
    USE JulDay_Mod,          ONLY : JulDay
    USE HCO_Clock_Mod
!
! !INPUT PARAMETERS:
!
    LOGICAL,                    INTENT(IN   ) :: am_I_Root   ! root CPU?
    TYPE(HCO_State),  POINTER                 :: HcoState    ! HEMCO state object 
    LOGICAL,                    INTENT(IN   ) :: ForceWrite  ! Write all diagnostics? 
    CHARACTER(LEN=*), OPTIONAL, INTENT(IN   ) :: PREFIX      ! File prefix
    LOGICAL,          OPTIONAL, INTENT(IN   ) :: UsePrevTime ! Use previous time 
    LOGICAL,          OPTIONAL, INTENT(IN   ) :: OnlyIfFirst ! Only write if nnDiagn is 1
    INTEGER,          OPTIONAL, INTENT(IN   ) :: COL         ! Collection Nr. 
!
! !INPUT/OUTPUT PARAMETERS:
!

    INTEGER,          INTENT(INOUT) :: RC          ! Failure or success
!
! !REVISION HISTORY: 
!  12 Sep 2013 - C. Keller   - Initial version 
!  11 Jun 2014 - R. Yantosca - Cosmetic changes in ProTeX headers
!  11 Jun 2014 - R. Yantosca - Now use F90 freeform indentation
!  19 Feb 2015 - C. Keller   - Added optional argument OnlyIfFirst
!  23 Feb 2015 - R. Yantosca - Now make Arr1D REAL(sp) so that we can write
!                              out lon & lat as float instead of double
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER                   :: I, PS, CNT, levIdTmp
    REAL(dp)                  :: GMT, JD1, JD1985, JD_DELTA, THISDAY
    REAL(sp)                  :: TMP, JD_DELTA_RND
    INTEGER                   :: YYYY, MM, DD, h, m, s
    REAL(sp), POINTER         :: nctime(:)
    REAL(sp), POINTER         :: Arr1D(:) => NULL()
    INTEGER,  POINTER         :: Int1D(:) => NULL()
    REAL(sp), POINTER         :: Arr3D(:,:,:) => NULL()
    REAL(sp), POINTER         :: Arr4D(:,:,:,:) => NULL()
    TYPE(DiagnCont), POINTER  :: ThisDiagn => NULL()
    INTEGER                   :: FLAG
    CHARACTER(LEN=255)        :: ncFile
    CHARACTER(LEN=255)        :: Pfx 
    CHARACTER(LEN=255)        :: MSG 
    CHARACTER(LEN=4 )         :: Yrs
    CHARACTER(LEN=2 )         :: Mts, Dys, hrs, mns 
    CHARACTER(LEN=31)         :: timeunit, myName, myUnit, OutOper
    INTEGER                   :: fId, lonId, latId, levId, TimeId
    INTEGER                   :: VarCt
    INTEGER                   :: nLon, nLat, nLev, nTime 
    INTEGER                   :: Prc
    INTEGER                   :: ResetFlag, MinResetFlag
    LOGICAL                   :: EOI, PrevTime
 
    CHARACTER(LEN=255), PARAMETER :: LOC = 'HCOIO_DIAGN_WRITEOUT (hcoio_diagn_mod.F90)' 
    !=================================================================
    ! HCOIO_DIAGN_WRITEOUT begins here!
    !=================================================================
  
    ! Init
    RC   = HCO_SUCCESS
    CNT  = 0

    ! Collection number
    PS = HcoDiagnIDDefault
    IF ( PRESENT(COL) ) PS = COL

    ! Get prefix and reset frequency of this diagnostics collection
    CALL DiagnCollection_Get( PS, Prefix=PFX, ResetFlag=ResetFlag, RC=RC ) 
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Check if there is at least one diagnostics to write:
    ! If current time stamp is not at the end of an interval - or
    ! if there is no diagnostics container in the list with a reset
    ! flag smaller or equal to MinResetFlag - there will be no matching
    ! container whatsoever. Can leave right here.
    ! EOI is the end-of-interval flag that will be used by routine
    ! Diagn_Get. If set to true, only the containers at the end of
    ! their averaging interval are returned.
    IF ( ForceWrite ) THEN
       MinResetFlag = -1
       EOI = .FALSE.
    ELSE
       MinResetFlag = HcoClock_GetMinResetFlag()
       EOI = .TRUE.
    ENDIF
    IF ( MinResetFlag > ResetFlag ) RETURN

    ! Inherit precision from HEMCO 
    Prc = HP
 
    ! Get PrevTime flag from input argument or set to default (=> TRUE)
    IF ( PRESENT(UsePrevTime) ) THEN
       PrevTime = UsePrevTime
    ELSE
       PrevTime = .TRUE.
    ENDIF

    !-----------------------------------------------------------------
    ! Create output file
    !-----------------------------------------------------------------

    ! Use HEMCO clock to create timestamp used in filename. Use previous
    ! time step if this option is selected.
    IF ( .NOT. PrevTime ) THEN
       CALL HcoClock_Get(cYYYY=YYYY,cMM=MM,cDD=DD,cH=h,cM=m,cS=s,RC=RC)
       IF ( RC /= HCO_SUCCESS ) RETURN
    ELSE
       CALL HcoClock_Get(pYYYY=YYYY,pMM=MM,pDD=DD,pH=h,pM=m,pS=s,RC=RC)
       IF ( RC /= HCO_SUCCESS ) RETURN
    ENDIF

    ! Define grid dimensions
    nLon  = HcoState%NX
    nLat  = HcoState%NY
    nLev  = HcoState%NZ
    nTime = 1 

    ! Initialize mirror variables
    allocate(Arr4D(nlon,nlat,nlev,ntime))
    allocate(Arr3D(nlon,nlat,ntime))
    Arr3D = 0.0_sp
    Arr4D = 0.0_sp

    ! Construct filename: diagnostics will be written into file
    ! PREFIX.YYYYMMDDhm.nc, where PREFIX is the input argument or
    ! (if not present) obtained from the HEMCO configuration file.
  
    ! Write datetime
    WRITE( Yrs, '(i4.4)' ) YYYY
    WRITE( Mts, '(i2.2)' ) MM
    WRITE( Dys, '(i2.2)' ) DD
    WRITE( hrs, '(i2.2)' ) h
    WRITE( mns, '(i2.2)' ) m

    ! Get prefix
    IF ( PRESENT(PREFIX) ) THEN
       Pfx = PREFIX
    ENDIF
    ncFile = TRIM(Pfx)//'.'//Yrs//Mts//Dys//hrs//mns//'.nc'

    ! verbose
    IF ( HCO_IsVerb(2) .AND. PS==1 ) THEN
       MSG = 'Write diagnostics into file '//TRIM(ncFile)
       CALL HCO_MSG( MSG )
    ENDIF

    ! Create output file
    CALL NC_CREATE( ncFile, nLon,  nLat,  nLev,  nTime, &
                    fId,    lonId, latId, levId, timeId, VarCt ) 

    !-----------------------------------------------------------------
    ! Write grid dimensions (incl. time) 
    !-----------------------------------------------------------------

    ! Add longitude 
    CALL NC_VAR_DEF ( fId, lonId, -1, -1, -1, &
                      'lon', 'Longitude', 'degrees_east', Prc, VarCt )
    ALLOCATE( Arr1D( nLon ) )
    Arr1D = HcoState%Grid%XMID%Val(:,1)
    CALL NC_VAR_WRITE ( fId, 'lon', Arr1D=Arr1D )
    DEALLOCATE( Arr1D )
    
    ! Add latitude
    CALL NC_VAR_DEF ( fId, -1, latId, -1, -1, &
                      'lat', 'Latitude', 'degrees_north', Prc, VarCt )
    ALLOCATE( Arr1D( nLat ) )
    Arr1D = HcoState%Grid%YMID%Val(1,:)
    CALL NC_VAR_WRITE ( fId, 'lat', Arr1D=Arr1D )
    DEALLOCATE( Arr1D )

    ! Add level 
    CALL NC_VAR_DEF ( fId, -1, levId, -1, -1, &
                      'lev', 'GEOS-Chem level', 'unitless', 1, VarCt )
    allocate(Int1D(nLev))
    DO I = 1, nLev
       Int1D(I) = I
    ENDDO
    CALL NC_VAR_WRITE ( fId, 'lev', Arr1D=Int1D )
    deallocate(Int1D)

    ! Add time 
    timeunit = 'hours since 1985-01-01 00:00:00 GMT'
    GMT = REAL(h,dp) + (REAL(m,dp)/60.0_dp) + (REAL(s,dp)/3600.0_dp)
    THISDAY  = DD + ( GMT / 24.0_dp )
    JD1      = JULDAY ( YYYY, MM, THISDAY )
    JD1985   = JULDAY ( 1985, 1,  0.0_dp  ) + 1.0_dp
    JD_DELTA = (JD1 - JD1985 ) * 24.0_dp

    ! Round to 2 digits after comma 
    JD_DELTA_RND = REAL(JD_DELTA,sp) * 100.0_sp
    TMP          = ANINT( JD_DELTA_RND )
    JD_DELTA_RND = TMP / 100.0_sp

    allocate(nctime(1)) 
    nctime(1) = JD_DELTA_RND
    CALL NC_VAR_DEF ( fId, -1, -1, -1, timeId, &
                      'time', 'Time', TRIM(timeunit), 4, VarCt )
    CALL NC_VAR_WRITE ( fId, 'time', Arr1D=nctime )
    deallocate(nctime)

    !-----------------------------------------------------------------
    ! Write out grid box areas 
    !-----------------------------------------------------------------
    myName = 'AREA'
    myUnit = 'm2'
    CALL NC_VAR_DEF ( fId, lonId, latId, -1, -1, &
                      TRIM(myName), 'Grid box area', TRIM(myUnit), Prc, VarCt )
    CALL NC_VAR_WRITE ( fId, TRIM(myName), Arr2D=HcoState%Grid%Area_M2%Val )

    !-----------------------------------------------------------------
    ! Write diagnostics 
    !-----------------------------------------------------------------

    ! Loop over all diagnostics in diagnostics list 
    ThisDiagn => NULL()
    DO WHILE ( .TRUE. )

       ! Get next diagnostics in list. This will return the next 
       ! diagnostics container that contains content. 
       CALL Diagn_Get ( am_I_Root, EOI, ThisDiagn, FLAG, RC, COL=PS ) 
       IF ( RC /= HCO_SUCCESS ) RETURN 
       IF ( FLAG /= HCO_SUCCESS ) EXIT

       ! Only write diagnostics if this is the first Diagn_Get call for
       ! this container and time step.
       IF ( PRESENT( OnlyIfFirst ) ) THEN
          IF ( OnlyIfFirst .AND. ThisDiagn%nnGetCalls > 1 ) CYCLE
       ENDIF 

       ! Define variable
       myName = ThisDiagn%cName
       myUnit = ThisDiagn%OutUnit
       IF ( ThisDiagn%SpaceDim == 3 ) THEN
          levIdTmp = levId
       ELSE
          levIdTmp = -1
       ENDIF

       ! Write out in single precision
       CALL NC_VAR_DEF ( fId, lonId, latId, levIdTmp, timeId, &
            TRIM(myName), TRIM(myName), TRIM(myUnit), SP, VarCt )

       ! Additional tracer attributes: long_name and _FillValue
       CALL NcBegin_Def( fID ) ! Reopen netCDF define mode
       CALL NcDef_var_attributes( fID, VarCt, "long_name",        &
            TRIM(ThisDiagn%long_name) )
       CALL NcDef_var_attributes( fID, VarCt, "averaging_method", &
            TRIM(ThisDiagn%AvgName  ) )
       CALL NcDef_var_attributes( fID, VarCt, "_FillValue",       &
            FillValue )
       CALL NcEnd_Def( fID )   ! Close netCDF define mode

       ! Mirror data and write to file. The mirroring is required in
       ! order to add the time dimension. Otherwise, the data would
       ! have no time information!
       IF ( ThisDiagn%SpaceDim == 3 ) THEN
          IF ( ASSOCIATED(ThisDiagn%Arr3D) ) THEN
             Arr4D(:,:,:,1) = ThisDiagn%Arr3D%Val
          ENDIF
          CALL NC_VAR_WRITE ( fId, TRIM(myName), Arr4D=Arr4D )
       ELSE
          IF ( ASSOCIATED(ThisDiagn%Arr2D) ) THEN
             Arr3D(:,:,1) = ThisDiagn%Arr2D%Val 
          ENDIF
          CALL NC_VAR_WRITE ( fId, TRIM(myName), Arr3D=Arr3D )
       ENDIF

       ! verbose
       IF ( HCO_IsVerb(2) .AND. PS==1 ) THEN
          MSG = '--- Added diagnostics: '//TRIM(myName)
          CALL HCO_MSG(MSG)
       ENDIF

    ENDDO

    !-----------------------------------------------------------------
    ! Close file
    !-----------------------------------------------------------------
    CALL NC_CLOSE ( fId )

    ! Cleanup
    Deallocate(Arr3D,Arr4D)
    ThisDiagn => NULL()
    
    ! Return 
    RC = HCO_SUCCESS

  END SUBROUTINE HCOIO_DiagN_WriteOut
!EOC
#else
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCOIO_Diagn_WriteOut
!
! !DESCRIPTION: Subroutine HCOIO\_Diagn\_WriteOut is the interface routine to
! link the HEMCO diagnostics arrays to the corresponding data pointers of the
! MAPL/ESMF history component. 
!\\
!\\
! Since the history component internally organizes many diagnostics tasks such
! as output scheduling, file writing, and data averaging, all HEMCO diagnostics
! are made available to the history component on every time step, e.g. the 
! entire content of the HEMCO diagnostics list is 'flushed' every time this
! subroutine is called. 
!\\
!\\
! For now, all diagnostics data is copied to the corresponding MAPL data
! pointer so that this routine works for cases where the HEMCO precision is
! not equal to the ESMF precision.
!\\
!\\
! Once the HEMCO precision is pegged to the ESMF precision, we can just 
! establish pointers between the export arrays and the diagnostics the first
! time this routine is called.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCOIO_Diagn_WriteOut( am_I_Root,  HcoState,    ForceWrite,  &
                                   RC,         PREFIX,      UsePrevTime, &
                                   InclManual, OnlyIfFirst, COL           )
!
! !USES:
!
    USE ESMF
    USE MAPL_MOD
    USE HCO_State_Mod, ONLY : HCO_State
    USE HCO_Clock_Mod, ONLY : HcoClock_GetMinResetFlag

# include "MAPL_Generic.h"
!
! !INPUT PARAMETERS:
!
    LOGICAL,                    INTENT(IN   ) :: am_I_Root   ! root CPU?
    TYPE(HCO_State),  POINTER                 :: HcoState    ! HEMCO state object 
    LOGICAL,                    INTENT(IN   ) :: ForceWrite  ! Write all diagnostics? 
    CHARACTER(LEN=*), OPTIONAL, INTENT(IN   ) :: PREFIX      ! File prefix
    LOGICAL,          OPTIONAL, INTENT(IN   ) :: UsePrevTime ! Use previous time 
    LOGICAL,          OPTIONAL, INTENT(IN   ) :: InclManual  ! Get manual diagn. too? 
    LOGICAL,          OPTIONAL, INTENT(IN   ) :: OnlyIfFirst !  
    INTEGER,          OPTIONAL, INTENT(IN   ) :: COL         ! Collection Nr. 
!
! !INPUT/OUTPUT PARAMETERS:
!

    INTEGER,          INTENT(INOUT) :: RC          ! Failure or success
!
! !REVISION HISTORY: 
!  05 Aug 2014 - C. Keller    - Initial version 
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    TYPE(DiagnCont), POINTER  :: ThisDiagn => NULL()
    INTEGER                   :: PS, FLAG, STAT
    CHARACTER(LEN=255)        :: MSG
    INTEGER                   :: MinResetFlag
    LOGICAL                   :: EOI, PrevTime
    REAL, POINTER             :: Ptr2D(:,:)   => NULL()
    REAL, POINTER             :: Ptr3D(:,:,:) => NULL()
    LOGICAL, SAVE             :: FIRST = .TRUE.

    CHARACTER(LEN=255), PARAMETER :: LOC = 'HCOIO_DIAGN_WRITEOUT (hcoio_diagn_mod.F90)'

    !=================================================================
    ! HCOIO_DIAGN_WRITEOUT begins here!
    !=================================================================

    ! Assume success until otherwise 
    RC  = HCO_SUCCESS

    ! If pointers are used, we need to call this routine only once
!    IF ( .NOT. FIRST ) RETURN

    ! Collection number
    PS = HcoDiagnIDDefault 
    IF ( PRESENT(COL) ) PS = COL

    ! In an ESMF environment, always get all diagnostics since output
    ! is scheduled through MAPL History!
    EOI = .FALSE.

    ! Get PrevTime flag from input argument or set to default (=> TRUE)
    IF ( PRESENT(UsePrevTime) ) THEN
       PrevTime = UsePrevTime
    ELSE
       PrevTime = .TRUE.
    ENDIF

    !-----------------------------------------------------------------
    ! Connect diagnostics to export state.
    !-----------------------------------------------------------------

    ! Loop over all diagnostics in diagnostics list 
    ThisDiagn => NULL()
    DO WHILE ( .TRUE. )

       ! Get next diagnostics in list. This will return the next 
       ! diagnostics container that contains content to be written
       ! out on this time step.
       CALL Diagn_Get ( am_I_Root, EOI, ThisDiagn, FLAG, RC, COL=PS ) 
       IF ( RC /= HCO_SUCCESS ) RETURN
       IF ( FLAG /= HCO_SUCCESS ) EXIT

       ! Only write diagnostics if this is the first Diagn_Get call for
       ! this container and time step.
       IF ( PRESENT(OnlyIfFirst) ) THEN 
          IF ( OnlyIfFirst .AND. ThisDiagn%nnGetCalls > 1 ) CYCLE
       ENDIF

       ! Get pointer to ESMF EXPORT field and pass data to it (if found):

       ! 2D...
       IF ( ThisDiagn%SpaceDim == 2 ) THEN
          CALL MAPL_GetPointer ( HcoState%EXPORT, Ptr2D, &
             TRIM(ThisDiagn%cName), NotFoundOk=.TRUE., RC=STAT )
          IF ( ASSOCIATED(Ptr2D) ) THEN
             IF ( ASSOCIATED(ThisDiagn%Arr2D) ) THEN
                Ptr2D = ThisDiagn%Arr2D%Val
                !Ptr2D => ThisDiagn%Arr2D%Val
             ENDIF
          ENDIF

       ! ... or 3D
       ELSEIF ( ThisDiagn%SpaceDim == 3 ) THEN
          CALL MAPL_GetPointer ( HcoState%EXPORT, Ptr3D, &
             TRIM(ThisDiagn%cName), NotFoundOk=.TRUE., RC=STAT )
          IF ( ASSOCIATED(Ptr3D) ) THEN
             IF ( ASSOCIATED(ThisDiagn%Arr3D) ) THEN
                Ptr3D(:,:,:) = ThisDiagn%Arr3D%Val(:,:,HcoState%NZ:1:-1)
                !Ptr3D => ThisDiagn%Arr3D%Val
             ENDIF
          ENDIF
       ENDIF

       ! Free pointer
       Ptr2D => NULL()
       Ptr3D => NULL()
    ENDDO

    ! Cleanup
    ThisDiagn => NULL()
    FIRST     = .FALSE.

    ! Return 
    RC = HCO_SUCCESS

  END SUBROUTINE HCOIO_DiagN_WriteOut
!EOC
#endif
END MODULE HCOIO_Diagn_Mod

