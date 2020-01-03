!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: hcoio_write_std_mod.F90
!
! !DESCRIPTION: Module HCOIO\_write\_std\_mod.F90 is the HEMCO data output
! interface for the 'standard' model environment. It contains routines to
! write out diagnostics into a netCDF file.
!\\
!\\
! !INTERFACE:
!
MODULE HCOIO_WRITE_STD_MOD
!
! !USES:
!
  USE HCO_ERROR_MOD
  USE HCO_DIAGN_MOD

  IMPLICIT NONE
  PRIVATE
#if !defined(ESMF_)
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC :: HCOIO_WRITE_STD
!
! !PRIVATE MEMBER FUNCTIONS:
!
  PRIVATE :: ConstructTimeStamp
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
!  22 Feb 2016 - C. Keller   - Split off from hcoio_diagn_mod.F90
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
! !IROUTINE: HCOIO_write_std
!
! !DESCRIPTION: Subroutine HCOIO\_write\_std writes diagnostics to
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
  SUBROUTINE HCOIO_write_std( HcoState, ForceWrite,  &
                              RC,          PREFIX,   UsePrevTime, &
                              OnlyIfFirst, COL                     )
!
! !USES:
!
    USE m_netCDF_io_define
    USE m_netcdf_io_read
    USE m_netcdf_io_open
    USE Ncdf_Mod,            ONLY : NC_Open
    USE Ncdf_Mod,            ONLY : NC_Read_Time
    USE Ncdf_Mod,            ONLY : NC_Read_Arr
    USE Ncdf_Mod,            ONLY : NC_Create
    USE Ncdf_Mod,            ONLY : NC_Close
    USE Ncdf_Mod,            ONLY : NC_Var_Def
    USE Ncdf_Mod,            ONLY : NC_Var_Write
    USE Ncdf_Mod,            ONLY : NC_Get_RefDateTime
    USE CHARPAK_Mod,         ONLY : TRANLC
    USE HCO_Chartools_Mod,   ONLY : HCO_CharParse
    USE HCO_State_Mod,       ONLY : HCO_State
    USE JulDay_Mod,          ONLY : JulDay
    USE HCO_EXTLIST_MOD,     ONLY : GetExtOpt, CoreNr
    USE HCO_Types_Mod,       ONLY : DiagnCont
    USE HCO_Clock_Mod

    ! Parameters for netCDF routines
    include "netcdf.inc"
!
! !INPUT PARAMETERS:
!
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
!  06 Nov 2015 - C. Keller   - Output time stamp is now determined from
!                              variable OutTimeStamp.
!  14 Jan 2016 - E. Lundgren - Create netcdf title out of filename prefix
!  20 Jan 2016 - C. Keller   - Added options DiagnRefTime and DiagnNoLevDim.
!  03 Mar 2016 - M. Sulprizio- Change netCDF format to netCDF-4
!  26 Oct 2016 - R. Yantosca - Don't nullify local ptrs in declaration stmts
!  21 Jan 2017 - C. Holmes   - Write all variable metadata in define mode, then
!                              switch to data mode just once. Much faster
!                              writing.
!  17 Feb 2017 - C. Holmes   - Enable netCDF-4 compression
!  08 Mar 2017 - R. Yantosca - Use unlimited time dimensions for netCDF files
!  29 Dec 2017 - C. Keller   - Now accept writing multiple time slices into
!                              same file.
!  03 Jan 2018 - R. Yantosca - Added more metadata for COARDS compliance.
!                              Also make TIME a 8-byte var to avoid roundoffs
!  05 Jan 2018 - R. Yantosca - Now print out all index variables as REAL*8
!  19 Oct 2018 - E. Lundgren - Disable writing multiple time slices to file
!                              until move of restart write from HEMCO to HISTORY
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER                   :: I, PS, CNT, levIdTmp, indexL, indexR
    REAL(dp)                  :: GMT, JD1, JD1985, JD_DELTA, THISDAY, P0
    REAL(sp)                  :: TMP, JD_DELTA_RND
    INTEGER                   :: YYYY, MM, DD, h, m, s
    REAL(sp), POINTER         :: nctime(:)
    REAL(dp), POINTER         :: Arr1D(:)
    INTEGER,  POINTER         :: Int1D(:)
    REAL(sp), POINTER         :: Arr3D(:,:,:)
    REAL(sp), POINTER         :: Arr4D(:,:,:,:)
    REAL(sp), POINTER         :: Arr4DOld(:,:,:,:)
    REAL*8,   POINTER         :: timeVec(:)
    REAL(hp), POINTER         :: hyam(:)
    REAL(hp), POINTER         :: hybm(:)
    TYPE(DiagnCont), POINTER  :: ThisDiagn
    INTEGER                   :: FLAG
    CHARACTER(LEN=255)        :: ncFile
    CHARACTER(LEN=255)        :: Pfx, title, Reference, Contact
    CHARACTER(LEN=255)        :: myLName, mySName, myFterm
    CHARACTER(LEN=255)        :: MSG
    CHARACTER(LEN=255)        :: RefTime
    CHARACTER(LEN=4 )         :: Yrs
    CHARACTER(LEN=2 )         :: Mts, Dys, hrs, mns
    CHARACTER(LEN=31)         :: myName, myUnit, OutOper
    CHARACTER(LEN=63)         :: timeunit
    INTEGER                   :: fId, lonId, latId, levId, TimeId
    INTEGER                   :: VarCt
    INTEGER                   :: nLon, nLat, nLev, nLevTmp, nTime
    INTEGER                   :: Prc,  L
    INTEGER                   :: lymd, lhms
    INTEGER                   :: refYYYY, refMM, refDD, refh, refm, refs
    LOGICAL                   :: EOI, DoWrite, PrevTime, FOUND
    LOGICAL                   :: NoLevDim, DefMode
    LOGICAL                   :: IsOldFile

    CHARACTER(LEN=255), PARAMETER :: LOC = 'HCOIO_WRITE_STD (hcoio_write_std_mod.F90)'

    !=================================================================
    ! HCOIO_WRITE_STD begins here!
    !=================================================================

    ! Init
    RC        =  HCO_SUCCESS
    CNT       =  0
    Arr1D     => NULL()
    Int1D     => NULL()
    Arr3D     => NULL()
    Arr4D     => NULL()
    Arr4DOld  => NULL()
    timeVec   => NULL()
    nctime    => NULL()
    ThisDiagn => NULL()

    ! Collection number
    PS = HcoState%Diagn%HcoDiagnIDDefault
    IF ( PRESENT(COL) ) PS = COL

    ! Check if it's time to write out this collection. Also set the
    ! end-of-interval EOI flag accordingly. This will be used lateron
    ! when calling Diagn_Get. Since all diagnostic containers in a
    ! given collection have the same output frequency, this is somewhat
    ! redundant (because we already check here if it is time to write
    ! out this particular collection). Keep it here for backwards
    ! consistency (ckeller, 8/6/2015).
    IF ( ForceWrite ) THEN
       DoWrite = .TRUE.
       EOI     = .FALSE.
    ELSE
       DoWrite = DiagnCollection_IsTimeToWrite( HcoState, PS )
       EOI     = .TRUE.
    ENDIF

    ! Create current time stamps (to be used to archive time stamps)
    CALL HcoClock_Get( HcoState%Clock,sYYYY=YYYY,sMM=MM,&
                       sDD=DD,sH=h,sM=m,sS=s,RC=RC)
    IF ( RC /= HCO_SUCCESS ) RETURN
    lymd = YYYY*10000 + MM*100 + DD
    lhms = h   *10000 + m *100 + s

    ! Leave here if it's not time to write diagnostics. On the first
    ! time step, set lastYMD and LastHMS to current dates.
    IF ( .NOT. DoWrite ) THEN
       IF ( .NOT. DiagnCollection_LastTimesSet(HcoState%Diagn,PS) ) THEN
          CALL DiagnCollection_Set ( HcoState%Diagn, COL=PS, &
                                     LastYMD=lymd, LastHMS=lhms, RC=RC )
       ENDIF
       RETURN
    ENDIF

    ! Inherit precision from HEMCO
    Prc = HP

    ! Get PrevTime flag from input argument or set to default (=> TRUE)
    IF ( PRESENT(UsePrevTime) ) THEN
       PrevTime = UsePrevTime
    ELSE
       PrevTime = .TRUE.
    ENDIF

    !-----------------------------------------------------------------
    ! Don't define level dimension if there are no 3D fields to write
    ! This is an optional feature. By default, all diagnostics have
    ! the full dimension definitions (lon,lat,lev,time) even if all
    ! output fields are only 2D. If the flag DiagnNoLevDim is
    ! enabled, the lev dimension is not defined if there are no 3D
    ! fields on the file.
    !-----------------------------------------------------------------
    NoLevDim = .FALSE.
    CALL GetExtOpt ( HcoState%Config, CoreNr, 'DiagnNoLevDim', &
                     OptValBool=NoLevDim, Found=Found, RC=RC )
    IF ( RC /= HCO_SUCCESS ) RETURN
    IF ( Found ) THEN
       IF ( NoLevDim ) THEN

          ! Loop over all diagnostics to see if any is 3D
          ThisDiagn => NULL()
          DO WHILE ( .TRUE. )

             ! Get next diagnostics in list. This will return the next
             ! diagnostics container that contains content.
             CALL Diagn_Get ( HcoState, EOI, &
                              ThisDiagn, FLAG, RC, COL=PS )
             IF ( RC /= HCO_SUCCESS ) RETURN
             IF ( FLAG /= HCO_SUCCESS ) EXIT

             ! If this is a 3D diagnostics, we must write the level
             ! coordinate
             IF ( ThisDiagn%SpaceDim == 3 ) THEN
                NoLevDim = .FALSE.
                EXIT
             ENDIF
          ENDDO
       ENDIF
    ENDIF

    !-----------------------------------------------------------------
    ! Create output file
    !-----------------------------------------------------------------

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
    CALL ConstructTimeStamp ( HcoState, PS, PrevTime, &
                              YYYY, MM, DD, h, m, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Write datetime
    WRITE( Yrs, '(i4.4)' ) YYYY
    WRITE( Mts, '(i2.2)' ) MM
    WRITE( Dys, '(i2.2)' ) DD
    WRITE( hrs, '(i2.2)' ) h
    WRITE( mns, '(i2.2)' ) m

    ! Get prefix
    IF ( PRESENT(PREFIX) ) THEN
       Pfx = PREFIX
    ELSE
       CALL DiagnCollection_Get( HcoState%Diagn, PS, PREFIX=Pfx, RC=RC )
       IF ( RC /= HCO_SUCCESS ) RETURN
    ENDIF
    ncFile = TRIM(Pfx)//'.'//Yrs//Mts//Dys//hrs//mns//'.nc'

    ! Multiple time slice update. Comment out for now since it causes
    ! timestamping the filename twice (ewl, 10/19/18)
    ! Add default time stamp if no time tokens are in the file template.
    ! This also ensures backward compatibility.
    !IF ( INDEX(TRIM(ncFile),'$') <= 0 ) THEN
    !   ncFile = TRIM(ncFile)//'.$YYYY$MM$DD$HH$MN.nc'
    !ENDIF
    !CALL HCO_CharParse ( HcoState%Config, ncFile, YYYY, MM, DD, h, m, RC )
    !IF ( RC /= HCO_SUCCESS ) RETURN

    ! Use filename prefix for title, replacing '_' with spaces
    ! NOTE: Prefix can only contain up to two underscores
    indexL = SCAN( Pfx, '_', .FALSE. ) ! Return left-most position
    indexR = SCAN( Pfx, '_', .TRUE.  ) ! Return right-most position
    IF ( indexL > 0 .AND. indexR > 0 ) THEN
       title = Pfx(1:indexL-1)        // ' ' //  &
               Pfx(indexL+1:indexR-1) // ' ' //  &
               Pfx(indexR+1:)
    ELSE IF ( indexL > 0 .AND. indexR == 0 ) THEN
       title = Pfx(1:indexL-1) // ' ' // Pfx(indexL+1:)
    ELSE
       title = Pfx
    ENDIF

    ! verbose
    IF ( HCO_IsVerb(HcoState%Config%Err,2) .AND. PS==1 ) THEN
       MSG = 'Write diagnostics into file '//TRIM(ncFile)
       CALL HCO_MSG( HcoState%Config%Err, MSG )
    ENDIF
    IF ( HCO_IsVerb(HcoState%Config%Err,3) .AND. PS==1 ) THEN
       WRITE(MSG,*) '--> write level dimension: ', .NOT.NoLevDim
       CALL HCO_MSG( HcoState%Config%Err, MSG )
    ENDIF

    ! Check if file already exists. If so, add new diagnostics to this file
    ! (instead of creating a new one)
    INQUIRE( FILE=ncFile, EXIST=IsOldFile )

    ! Disable multiple time slice update since causes an issue writing
    ! restart files. Re-enable when restart files are written via HISTORY
    ! rather than HEMCO by deleting the forcing of IsOldFile below.
    ! (ewl, 10/19/18)
    IsOldFile = .FALSE.

    ! If file exists, open file and get time dimension
    IF ( IsOldFile ) THEN
       CALL Ncop_Wr( fID, ncFile )
       CALL NC_READ_TIME( fID, ntime, timeunit, timeVec, RC=RC )

       ! new file will have one more time dimension
       ntime = ntime + 1

    ! Create output file
    ELSE

       ! Define a variable for the number of levels, which will either be -1
       ! (if all 2D data) or the number of levels in the grid (for 3D data).
       IF ( NoLevDim ) THEN
          nLevTmp = -1
       ELSE
          nLevTmp = nLev
       ENDIF

       ! Define extra metadata for global attributes
       Reference = 'http://wiki.geos-chem.org/The_HEMCO_Users_Guide'
       Contact   = 'GEOS-Chem Support Team (geos-chem-support@as.harvard.edu)'

       ! Create output file
       ! Pass CREATE_NC4 to make file format netCDF-4 (mps, 3/3/16)
       ! Now create netCDF file with time dimension as UNLIMITED (bmy, 3/8/17)
       CALL NC_Create( NcFile       = NcFile,                            &
                       Title        = Title,                             &
                       Reference    = Reference,                         &
                       Contact      = Contact,                           &
                       nLon         = nLon,                              &
                       nLat         = nLat,                              &
                       nLev         = nLevTmp,                           &
                       nTime        = NF_UNLIMITED,                      &
                       fId          = fId,                               &
                       lonId        = lonId,                             &
                       latId        = latId,                             &
                       levId        = levId,                             &
                       timeId       = timeId,                            &
                       VarCt        = VarCt,                             &
                       CREATE_NC4   =.TRUE.                             )

    ENDIF

    !-----------------------------------------------------------------
    ! Write grid dimensions (incl. time)
    !-----------------------------------------------------------------
    IF ( .NOT. IsOldFile ) THEN

       ! Write longitude axis variable ("lon") to file
       CALL NC_Var_Def( fId         = fId,                                &
                        lonId       = lonId,                              &
                        latId       = -1,                                 &
                        levId       = -1,                                 &
                        timeId      = -1,                                 &
                        VarName     = 'lon',                              &
                        VarLongName = 'Longitude',                        &
                        VarUnit     = 'degrees_east',                     &
                        Axis        = 'X',                                &
                        DataType    = dp,                                 &
                        VarCt       = VarCt,                              &
                        Compress    = .TRUE.                             )
       ALLOCATE( Arr1D( nLon ) )
       Arr1D = HcoState%Grid%XMID%Val(:,1)
       CALL NC_Var_Write( fId, 'lon', Arr1D=Arr1D )
       DEALLOCATE( Arr1D )

       ! Write latitude axis variable ("lat") to file
       CALL NC_Var_Def( fId         = fId,                              &
                        lonId       = -1,                               &
                        latId       = latId,                            &
                        levId       = -1,                               &
                        timeId      = -1,                               &
                        VarName     = 'lat',                            &
                        VarLongName = 'Latitude',                       &
                        VarUnit     = 'degrees_north',                  &
                        Axis        = 'Y',                              &
                        DataType    = dp,                               &
                        VarCt       = VarCt,                            &
                        Compress    = .TRUE.                           )
       ALLOCATE( Arr1D( nLat ) )
       Arr1D = HcoState%Grid%YMID%Val(1,:)
       CALL NC_Var_Write( fId, 'lat', Arr1D=Arr1D )
       DEALLOCATE( Arr1D )

       ! Write vertical grid parameters to file (if necessary)
       IF ( .NOT. NoLevDim ) THEN

          ! Reference pressure [Pa]
          P0 = 1.0e+05_dp

          ! Allocate vertical coordinate arrays
          ALLOCATE( Arr1D( nLev ) )
          ALLOCATE( hyam ( nLev ) )
          ALLOCATE( hybm ( nLev ) )

          ! Construct vertical level coordinates
          DO L = 1, nLev

             ! A parameter at grid midpoints
             hyam(L)  = ( HcoState%Grid%zGrid%Ap(L)                         &
                      +   HcoState%Grid%zGrid%Ap(L+1) ) * 0.5_dp

             ! B parameter at grid midpoints
             hybm(L)  = ( HcoState%Grid%zGrid%Bp(L)                         &
                      +   HcoState%Grid%zGrid%Bp(L+1) ) * 0.5_dp

             ! Vertical level coordinate
             Arr1d(L) = ( hyam(L) / P0 ) + hybm(L)

          ENDDO

          ! Write level axis variable ("lev") to file
          ! Define extra metadata for calls to NC_Var_Def
          myLName = 'hybrid level at midpoints ((A/P0)+B)'
          mySName = 'atmosphere_hybrid_sigma_pressure_coordinate'
          myFTerm = 'a: hyai b: hybi p0: P0 ps: PS'
          CALL NC_Var_Def( fId          = fId,                            &
                           lonId        = -1,                             &
                           latId        = -1,                             &
                           levId        = levId,                          &
                           timeId       = -1,                             &
                           VarName      = 'lev',                          &
                           VarLongName  = MyLName,                        &
                           StandardName = MySName,                        &
                           FormulaTerms = myFTerm,                        &
                           VarUnit      = 'level',                        &
                           Axis         = 'Z',                            &
                           Positive     = 'up',                           &
                           DataType     = dp,                             &
                           VarCt        = VarCt,                          &
                           Compress     = .TRUE.                         )
          CALL NC_Var_Write( fId, 'lev', Arr1D=Arr1D )

          ! Write hybrid A coordinate ("hyam") to file
          ! Define extra metadata for calls to NC_Var_Def
          myLName = 'hybrid A coefficient at layer midpoints'
          CALL NC_Var_Def( fId          = fId,                            &
                           lonId        = -1,                             &
                           latId        = -1,                             &
                           levId        = levId,                          &
                           timeId       = -1,                             &
                           VarName      = 'hyam',                         &
                           VarLongName  = MyLName,                        &
                           VarUnit      = 'Pa',                           &
                           DataType     = dp,                             &
                           VarCt        = VarCt,                          &
                           Compress     = .TRUE.                         )
          CALL NC_Var_Write ( fId, 'hyam', Arr1D=hyam )

          ! Write hybrid B coordinate ("hybm") to file
          ! Define extra metadata for calls to NC_Var_Def
          myLName = 'hybrid B coefficient at layer midpoints'
          CALL NC_Var_Def( fId          = fId,                           &
                           lonId        = -1,                            &
                           latId        = -1,                            &
                           levId        = levId,                         &
                           timeId       = -1,                            &
                           VarName      = 'hybm',                        &
                           VarLongName  = MyLName,                       &
                           VarUnit      = '1',                           &
                           DataType     = dp,                            &
                           VarCt        = VarCt,                         &
                           Compress     = .TRUE.                        )
          CALL NC_Var_Write( fId, 'hybm', Arr1D=hybm )

          ! Write out reference pressure (P0) to file
          CALL NC_Var_Def( fId         = fId,                             &
                           lonId       = -1,                              &
                           latId       = -1,                              &
                           levId       = -1,                              &
                           timeId      = -1,                              &
                           VarName     = 'P0',                            &
                           VarLongName = 'Reference pressure',            &
                           VarUnit     = 'Pa',                            &
                           DataType    = dp,                              &
                           VarCt       = VarCt,                           &
                           Compress    = .TRUE.                          )
          CALL NC_Var_Write( fId, 'P0', P0 )

          ! Deallocate arrays
          DEALLOCATE( Arr1d )
          DEALLOCATE( hyam  )
          DEALLOCATE( hybm  )

       ENDIF
    ENDIF

    !------------------------------------------------------------------------
    ! Write time axis variable ("time") to file
    !------------------------------------------------------------------------

    ! JD1 is the julian day of the data slice
    GMT     = REAL(h,dp) + (REAL(m,dp)/60.0_dp) + (REAL(s,dp)/3600.0_dp)
    THISDAY = DD + ( GMT / 24.0_dp )
    JD1     = JULDAY ( YYYY, MM, THISDAY )

    ! Check if reference time is given in HEMCO configuration file
    CALL GetExtOpt ( HcoState%Config, CoreNr, 'DiagnRefTime', &
                     OptValChar=RefTime, Found=Found, RC=RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Use specified reference time (if available)
    IF ( Found ) THEN
       timeunit = ADJUSTL(TRIM(RefTime))
       CALL TRANLC( timeunit )
       CALL NC_GET_REFDATETIME( timeunit, refYYYY, refMM, refDD, refh, &
                                refm, refs, RC )
       refs = 0
       IF ( RC /= HCO_SUCCESS ) RETURN
       GMT      = REAL(MAX(refh,0),dp) + (REAL(MAX(refm,0),dp)/60.0_dp) + &
                  (REAL(MAX(refs,0),dp)/3600.0_dp)
       THISDAY  = refDD + ( GMT / 24.0_dp )
       JD1985   = JULDAY ( refYYYY, refMM, THISDAY )

    ! Use current time if not found
    ELSE
       WRITE(timeunit,100) YYYY,MM,DD,h,m,s
       JD1985 = JD1
    ENDIF
100 FORMAT ( 'hours since ',i4.4,'-',i2.2,'-',i2.2,' ',i2.2,':',i2.2,':',i2.2,' GMT' )

    ! Calculate time value
    JD_DELTA = (JD1 - JD1985 )

    ! Default is 'days since'. Adjust for 'hours since', 'minutes since',
    ! 'seconds since'.
    IF ( timeunit(1:4) == 'days' ) THEN
       ! all ok
    ELSEIF ( timeunit(1:5) == 'hours' ) THEN
       JD_DELTA = JD_DELTA * 24.0_dp
    ELSEIF ( timeunit(1:7) == 'minutes' ) THEN
       JD_DELTA = JD_DELTA * 24.0_dp * 60.0_dp
    ELSEIF ( timeunit(1:7) == 'seconds' ) THEN
       JD_DELTA = JD_DELTA * 24.0_dp * 3600.0_dp
    ELSE
       MSG = 'Unrecognized output reference time, will ' // &
             'assume `days since`: '//TRIM(timeunit)
       CALL HCO_WARNING( MSG, WARNLEV=2, THISLOC=LOC, RC=RC )
    ENDIF

    ! Special case where we have an old file but it has the same time stamp: in
    ! that case simply overwrite the current values
    ! Comment out code for single precision rounded time (ewl, 10/18/18)
    !IF ( IsOldFile .AND. ntime == 2 .AND. timeVec(1) == JD_DELTA_RND ) THEN
    IF ( IsOldFile .AND. ntime == 2 ) THEN
       IF ( timeVec(1) == JD_DELTA ) THEN
          ntime = 1
       ENDIF
    ENDIF
    ALLOCATE( nctime(ntime) )
    IF ( IsOldFile .AND. ntime > 1 ) THEN
       nctime(1:ntime-1) = timeVec(:)
    ENDIF
    nctime(ntime) = JD_DELTA

    IF ( .NOT. IsOldFile ) THEN
       CALL NC_Var_Def( fId         = fId,                                &
                        lonId       = -1,                                 &
                        latId       = -1,                                 &
                        levId       = -1,                                 &
                        timeId      = timeId,                             &
                        VarName     = 'time',                             &
                        VarLongName = 'Time',                             &
                        VarUnit     = TimeUnit,                           &
                        Axis        = 'T',                                &
                        Calendar    = 'gregorian',                        &
                        DataType    = 8,                                  &
                        VarCt       = VarCt,                              &
                        Compress    = .TRUE.                             )
    ENDIF
    CALL NC_VAR_WRITE( fId, 'time', Arr1D=nctime )
    DEALLOCATE( nctime )
    IF ( ASSOCIATED(timeVec) ) DEALLOCATE( timeVec )

    !-----------------------------------------------------------------
    ! Write out grid box areas
    !-----------------------------------------------------------------

    IF ( .NOT. IsOldFile ) THEN
       CALL NC_Var_Def( fId         = fId,                                &
                        lonId       = lonId,                              &
                        latId       = latId,                              &
                        levId       = -1,                                 &
                        timeId      = -1,                                 &
                        VarName     = 'AREA',                             &
                        VarLongName = 'Grid box area',                    &
                        VarUnit     = 'm2',                               &
                        DataType    = Prc,                                &
                        VarCt       = VarCt,                              &
                        Compress    = .TRUE.                             )
       CALL NC_Var_Write ( fId, 'AREA', Arr2D=HcoState%Grid%Area_M2%Val )
    ENDIF

    !-----------------------------------------------------------------
    ! Write diagnostics
    !-----------------------------------------------------------------

    ! Run this section twice, first in define mode for metadata, then in
    ! data mode to write variables
    DO I=1,2

    ! Skip definition mode for existing file
    IF ( I==1 .AND. IsOldFile ) CYCLE

    IF (I==1) THEN
       ! Open netCDF define mode
       CALL NcBegin_Def( fID )
       DefMode=.TRUE.
    ELSE
!       IF ( .NOT. IsOldFile ) THEN
          ! Close netCDF define mode
          CALL NcEnd_Def( fID )
!       ENDIF
       DefMode=.False.
    ENDIF

    ! Loop over all diagnostics in diagnostics list
    ThisDiagn => NULL()
    DO WHILE ( .TRUE. )

       ! Get next diagnostics in list. This will return the next
       ! diagnostics container that contains content.
       CALL Diagn_Get ( HcoState, EOI, ThisDiagn, FLAG, RC, COL=PS )
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

       ! Error check: this should never happen!
       IF ( levIdTmp > 0 .AND. NoLevDim ) THEN
          MSG = 'Level dimension undefined but 3D container found: ' &
                // TRIM(myName)
          CALL HCO_ERROR(MSG,RC,THISLOC=LOC)
          RETURN
       ENDIF

       IF (DefMode) THEN

          !------------------------------------
          ! Define variables in define mode
          !------------------------------------

          ! Define variable as single precision
          CALL NC_Var_Def( fId          = fId,                               &
                           lonId        = lonId,                             &
                           latId        = latId,                             &
                           levId        = levIdTmp,                          &
                           timeId       = timeId,                            &
                           VarName      = TRIM(myName),                      &
                           VarLongName  = ThisDiagn%long_name,               &
                           VarUnit      = TRIM(myUnit),                      &
                           AvgMethod    = ThisDiagn%AvgName,                 &
                           MissingValue = FillValue,                         &
                           DataType     = sp,                                &
                           VarCt        = VarCt,                             &
                           DefMode      = DefMode,                           &
                           Compress     = .True.                            )

       ELSE

          !------------------------------------
          ! Write variables in data mode
          !------------------------------------

          IF ( IsOldFile .AND. ntime > 1 ) THEN
             IF ( ThisDiagn%SpaceDim == 3 ) THEN
                CALL NC_READ_ARR( fID, TRIM(myName), 1, nlon, 1, nlat, &
                                  1, nlev, 1, ntime-1, ncArr=Arr4DOld, RC=RC )
                Arr4D(:,:,:,1:ntime-1) = Arr4DOld(:,:,:,:)
             ELSE
                CALL NC_READ_ARR( fID, TRIM(myName), 1, nlon, 1, nlat, &
                                  -1, -1, 1, ntime-1, ncArr=Arr4DOld, RC=RC )
                Arr3D(:,:,1:ntime-1) = Arr4DOld(:,:,1,:)
             ENDIF
             IF ( ASSOCIATED(Arr4DOld) ) DEALLOCATE(Arr4DOld)
          ENDIF

          ! Mirror data and write to file. The mirroring is required in
          ! order to add the time dimension. Otherwise, the data would
          ! have no time information!
          IF ( ThisDiagn%SpaceDim == 3 ) THEN
             IF ( ASSOCIATED(ThisDiagn%Arr3D) ) THEN
                Arr4D(:,:,:,ntime) = ThisDiagn%Arr3D%Val
                Arr4D(:,:,:,1) = ThisDiagn%Arr3D%Val
             ENDIF
             CALL NC_VAR_WRITE ( fId, TRIM(myName), Arr4D=Arr4D )
          ELSE
             IF ( ASSOCIATED(ThisDiagn%Arr2D) ) THEN
                Arr3D(:,:,ntime) = ThisDiagn%Arr2D%Val
                Arr3D(:,:,1) = ThisDiagn%Arr2D%Val
             ENDIF
             CALL NC_VAR_WRITE ( fId, TRIM(myName), Arr3D=Arr3D )
          ENDIF

          ! verbose
          IF ( HCO_IsVerb(HcoState%Config%Err,2) .AND. PS==1 ) THEN
             MSG = '--- Added diagnostics: '//TRIM(myName)
             CALL HCO_MSG(HcoState%Config%Err,MSG)
          ENDIF
       ENDIF
    ENDDO
    ENDDO

    !-----------------------------------------------------------------
    ! Cleanup
    !-----------------------------------------------------------------

    ! Close file
    CALL NC_CLOSE ( fId )

    ! Cleanup local variables
    Deallocate(Arr3D,Arr4D)
    ThisDiagn => NULL()

    ! Archive time stamp
    CALL DiagnCollection_Set ( HcoState%Diagn, COL=PS, &
                               LastYMD=lymd, LastHMS=lhms, RC=RC )

    ! Return
    RC = HCO_SUCCESS

  END SUBROUTINE HCOIO_write_std
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ConstructTimeStamp
!
! !DESCRIPTION: Subroutine ConstructTimeStamp is a helper routine to construct
! the time stamp of a given diagnostics collection.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE ConstructTimeStamp ( HcoState, PS, PrevTime, Yr, Mt, Dy, hr, mn, RC )
!
! !USES:
!
    USE HCO_State_Mod,       ONLY : HCO_State
    USE HCO_Clock_Mod
    USE JULDAY_MOD
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(HCO_State), POINTER          :: HcoState     ! HEMCO state obj
    INTEGER,         INTENT(IN   )    :: PS           ! collecion ID
    LOGICAL,         INTENT(IN   )    :: PrevTime     ! Use previous time?
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,         INTENT(INOUT)    :: RC           ! Return code
!
! !OUTPUT PARAMETERS:
!
    INTEGER,         INTENT(  OUT)    :: Yr
    INTEGER,         INTENT(  OUT)    :: Mt
    INTEGER,         INTENT(  OUT)    :: Dy
    INTEGER,         INTENT(  OUT)    :: hr
    INTEGER,         INTENT(  OUT)    :: mn
!
! !REVISION HISTORY:
!  06 Nov 2015 - C. Keller   - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER            :: Y2, M2, D2, h2, n2, s2
    INTEGER            :: Y1, M1, D1, h1, n1, s1
    INTEGER            :: LastYMD, LastHMS
    INTEGER            :: YYYYMMDD, HHMMSS
    INTEGER            :: OutTimeStamp
    REAL(dp)           :: DAY, UTC, JD1, JD2, JDMID
    CHARACTER(LEN=255) :: MSG
    CHARACTER(LEN=255) :: LOC = 'ConstuctTimeStamp (hcoi_diagn_mod.F90)'

    !=================================================================
    ! ConstructTimeStamp begins here!
    !=================================================================

    ! Use HEMCO clock to create timestamp used in filename. Use previous
    ! time step if this option is selected.
    IF ( .NOT. PrevTime ) THEN
       CALL HcoClock_Get(HcoState%Clock,sYYYY=Y2,sMM=M2,&
                         sDD=D2,sH=h2,sM=n2,sS=s2,RC=RC)
       IF ( RC /= HCO_SUCCESS ) RETURN
    ELSE
       CALL HcoClock_Get(HcoState%Clock,pYYYY=Y2,pMM=M2,&
                         pDD=D2,pH=h2,pM=n2,pS=s2,RC=RC)
       IF ( RC /= HCO_SUCCESS ) RETURN
    ENDIF

    ! Get timestamp location for this collection
    CALL DiagnCollection_Get( HcoState%Diagn, PS, OutTimeStamp=OutTimeStamp, &
                              LastYMD=LastYMD, LastHMS=LastHMS, RC=RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Determine dates to be used:

    ! To use start date
    IF ( OutTimeStamp == HcoDiagnStart ) THEN
       Yr = FLOOR( MOD(LastYMD*1.d0, 100000000.d0 ) / 1.0d4 )
       Mt = FLOOR( MOD(LastYMD*1.d0, 10000.d0     ) / 1.0d2 )
       Dy = FLOOR( MOD(LastYMD*1.d0, 100.d0       ) / 1.0d0 )
       Hr = FLOOR( MOD(LastHMS*1.d0, 1000000.d0   ) / 1.0d4 )
       Mn = FLOOR( MOD(LastHMS*1.d0, 10000.d0     ) / 1.0d2 )

    ! Use mid point
    ELSEIF ( OutTimeStamp == HcoDiagnMid ) THEN

       ! Julian day of start interval:
       Y1 = FLOOR( MOD(LastYMD*1.d0, 100000000.d0 ) / 1.0d4 )
       M1 = FLOOR( MOD(LastYMD*1.d0, 10000.d0     ) / 1.0d2 )
       D1 = FLOOR( MOD(LastYMD*1.d0, 100.d0       ) / 1.0d0 )
       h1 = FLOOR( MOD(LastHMS*1.d0, 1000000.d0   ) / 1.0d4 )
       n1 = FLOOR( MOD(LastHMS*1.d0, 10000.d0     ) / 1.0d2 )
       s1 = FLOOR( MOD(LastHMS*1.d0, 100.d0       ) / 1.0d0 )

       UTC = ( REAL(h1,dp) / 24.0_dp    ) + &
             ( REAL(n1,dp) / 1440.0_dp  ) + &
             ( REAL(s1,dp) / 86400.0_dp )
       DAY = REAL(D1,dp) + UTC
       JD1 = JULDAY( Y1, M1, DAY )

       ! Julian day of end interval:
       UTC = ( REAL(h2,dp) / 24.0_dp    ) + &
             ( REAL(n2,dp) / 1440.0_dp  ) + &
             ( REAL(s2,dp) / 86400.0_dp )
       DAY = REAL(D2,dp) + UTC
       JD2 = JULDAY( Y2, M2, DAY )

       ! Julian day in the middle
       JDMID = ( JD1 + JD2 ) / 2.0_dp

       ! Tranlate back into dates
       CALL CALDATE( JDMID, YYYYMMDD, HHMMSS )
       Yr = FLOOR ( MOD( YYYYMMDD, 100000000) / 1.0e4_dp )
       Mt = FLOOR ( MOD( YYYYMMDD, 10000    ) / 1.0e2_dp )
       Dy = FLOOR ( MOD( YYYYMMDD, 100      ) / 1.0e0_dp )
       Hr = FLOOR ( MOD(   HHMMSS, 1000000  ) / 1.0e4_dp )
       Mn = FLOOR ( MOD(   HHMMSS, 10000    ) / 1.0e2_dp )

    ! Otherwise, use end date
    ELSE
       Yr = Y2
       Mt = M2
       Dy = D2
       Hr = h2
       Mn = n2
    ENDIF

    ! Return w/ success
    RC = HCO_SUCCESS

  END SUBROUTINE ConstructTimeStamp
!EOC
#endif
END MODULE HCOIO_WRITE_STD_MOD

