!BOM
# if !defined(ESMF_)
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: hcoi_gc_diagn_mod.F90
!
! !DESCRIPTION: Module HCOI\_GC\_Diagn\_Mod.F90 is the HEMCO - GEOS-Chem
! diagnostics interface module.
!\\
!\\
! This module contains wrapper routines to initialize, execute and finalize
! HEMCO diagnostics from within GEOS-Chem. These routines are called from
! the HEMCO- GEOS-Chem interface module (hcoi\_gc\_main\_mod.F90). 
! \\
! !INTERFACE:
!
MODULE HCOI_GC_DIAGN_MOD
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
  PUBLIC :: HCOI_DIAGN_INIT
  PUBLIC :: HCOI_DIAGN_WRITEOUT
!
! !REMARKS:
!  HEMCO diagnostics are still in testing.  We will fully activate them
!  at a later time.  They will be turned on when debugging & unit testing.
!
! !REVISION HISTORY:
!  04 May 2014 - C. Keller   - Initial version. 
!  11 Jun 2014 - R. Yantosca - Cosmetic changes in ProTeX headers
!  11 Jun 2014 - R. Yantosca - Now use F90 freeform indentation
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCOI_DiagN_Init
!
! !DESCRIPTION: Subroutine HCOI\_Diagn\_Init initializes the HEMCO diagnostics
! in GEOS-Chem. 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCOI_Diagn_Init( am_I_Root, HcoState, RC ) 
!
! !USES:
!
    USE HCO_State_Mod, ONLY : HCO_GetHcoID, HCO_State
!
! !INPUT PARAMETERS:
!
    LOGICAL,          INTENT(IN   )  :: am_I_Root  ! root CPU?
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(HCO_State),  POINTER        :: HcoState   ! HEMCO state object 
    INTEGER,          INTENT(INOUT)  :: RC         ! Failure or success
!
! !REVISION HISTORY: 
!  12 Sep 2013 - C. Keller   - Initial version 
!  11 Jun 2014 - R. Yantosca - Cosmetic changes in ProTeX headers
!  11 Jun 2014 - R. Yantosca - Now use F90 freeform indentation
!  25 Jul 2014 - R. Yantosca - Now save diagnostics at the end of the run,
!                              which facilitates unit testing & diff testing
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER                         :: I, N, AS
    CHARACTER(LEN=255)              :: MSG, LOC

    !=================================================================
    ! HCOI_DIAGN_INIT begins here!
    !=================================================================

    ! Init 
    LOC = 'HCOI_DIAGN_INIT (hcoi_gc_diagn_mod.F90)'
    RC  = HCO_SUCCESS

       !=================================================================
       ! Define manual diagnostics (no AutoFill)
       !=================================================================

       !=================================================================
       ! Define automatic diagnostics (AutoFill)
       !=================================================================

       ! Total NO
       I = HCO_GetHcoID( 'NO', HcoState )
       CALL Diagn_Create ( am_I_Root, &
                           HcoState,  &
                           cName    = 'NOtotal', &
                           ExtNr    = -1, &
                           Cat      = -1, &
                           Hier     = -1, &
                           HcoID    = I, &
                           SpaceDim = 2, &
                           LevIDx   = -1, &
                           OutUnit  = 'kg/m2/s', &
!                           WriteFreq = 'Hourly',  &
                           WriteFreq = 'End',  &
                           AutoFill  = 1, &
                           cID       = N, & 
                           RC        = RC ) 
       IF ( RC /= HCO_SUCCESS ) RETURN 

       ! NO from anthropogenic emissions 
       CALL Diagn_Create ( am_I_Root, &
                           HcoState,  &
                           cName    = 'NO_anthr', &
                           ExtNr    = 0,  &
                           Cat      = 1, &
                           Hier     = -1, &
                           HcoID    = I, &
                           SpaceDim = 2, &
                           OutUnit  = 'kg/m2/s', &
!                           WriteFreq = 'Hourly',  &
                           WriteFreq = 'End',  &
                           AutoFill  = 1, &
                           cID       = N, & 
                           RC        = RC ) 
       IF ( RC /= HCO_SUCCESS ) RETURN 

       ! NO from aircrafts (AEIC) 
       CALL Diagn_Create ( am_I_Root, &
                           HcoState,  &
                           cName     = 'NO_aircraft', &
                           ExtNr     = 0,  &
                           Cat       = 20, &
                           Hier      = -1, &
                           HcoID     = I, &
                           SpaceDim  = 2, &
                           LevIdx    = -1, &
                           OutUnit   = 'kg/m2/s', &
!                           WriteFreq = 'Hourly',  &
                           WriteFreq = 'End',  &
                           AutoFill  = 1, &
                           cID       = N, & 
                           RC        = RC ) 
       IF ( RC /= HCO_SUCCESS ) RETURN 

       ! NO from biomass burning 
       CALL Diagn_Create ( am_I_Root, &
                           HcoState,  &
                           cName    = 'NO_biomass', &
                           ExtNr    = 111, &
                           Cat      = -1,  &
                           Hier     = -1,  &
                           HcoID    = I, &
                           SpaceDim = 2, &
                           OutUnit  = 'kg/m2/s', &
!                           WriteFreq = 'Hourly',  &
                           WriteFreq = 'End',  &
                           AutoFill  = 1, &
                           cID       = N, & 
                           RC        = RC ) 
       IF ( RC /= HCO_SUCCESS ) RETURN 

       ! NO from ships (ParaNox) 
       I = HCO_GetHcoID( 'NO', HcoState )
       CALL Diagn_Create ( am_I_Root, &
                           HcoState,  &
                           cName    = 'ParaNOx', &
                           ExtNr    = 102,  &
                           Cat      = -1, &
                           Hier     = -1, &
                           HcoID    = I, &
                           SpaceDim = 2, &
                           LevIDx   = -1, &
                           OutUnit  = 'kg/m2/s', &
!                           WriteFreq = 'Hourly',  &
                           WriteFreq = 'End',  &
                           AutoFill  = 1, &
                           cID       = N, & 
                           RC        = RC ) 
       IF ( RC /= HCO_SUCCESS ) RETURN 

       ! NO from soils 
       I = HCO_GetHcoID( 'NO', HcoState )
       CALL Diagn_Create ( am_I_Root, &
                           HcoState,  &
                           cName    = 'SoilNOx', &
                           ExtNr    = 104,  &
                           Cat      = -1, &
                           Hier     = -1, &
                           HcoID    = I, &
                           SpaceDim = 2, &
                           OutUnit  = 'kg/m2/s', &
!                           WriteFreq = 'Hourly',  &
                           WriteFreq = 'End',  &
                           AutoFill  = 1, &
                           cID       = N, & 
                           RC        = RC ) 
       IF ( RC /= HCO_SUCCESS ) RETURN 

       ! NO from lightning 
       I = HCO_GetHcoID( 'NO', HcoState )
       CALL Diagn_Create ( am_I_Root, &
                           HcoState,  &
                           cName    = 'LightNOx', &
                           ExtNr    = 103,  &
                           Cat      = -1, &
                           Hier     = -1, &
                           HcoID    = I, &
                           SpaceDim = 2, &
                           LevIDx   = -1, &
                           OutUnit  = 'kg/m2/s', &
!                           WriteFreq = 'Hourly',  &
                           WriteFreq = 'End',  &
                           AutoFill  = 1, &
                           cID       = N, & 
                           RC        = RC ) 
       IF ( RC /= HCO_SUCCESS ) RETURN 

!       ! Total SO2
!       I = HCO_GetHcoID( 'SO2', HcoState )
!       CALL Diagn_Create ( am_I_Root, &
!                           HcoState,  &
!                           cName    = 'SO2total', &
!                           ExtNr    = -1, &
!                           Cat      = -1, &
!                           Hier     = -1, &
!                           HcoID    = I, &
!                           SpaceDim = 2, &
!                           LevIDx   = -1, &
!                           OutUnit  = 'kg/m2/s', &
!                           WriteFreq = 'Hourly',  &
!                           AutoFill  = 1, &
!                           cID       = N, & 
!                           RC        = RC ) 
!       IF ( RC /= HCO_SUCCESS ) RETURN 

    ! Leave w/ success
    RC = HCO_SUCCESS 

  END SUBROUTINE HCOI_Diagn_Init
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCOI_Diagn_WriteOut
!
! !DESCRIPTION: Subroutine HCOI\_Diagn\_WriteOut writes diagnostics to 
! netCDF file. If the WriteAll flag is set to TRUE, all diagnostics are
! written. This option is usually only used at the end of a simulation run.
! If WriteAll is False, only the diagnostics that are at the end of their
! time averaging interval are written. For example, if the current month
! is different from the previous (emissions) month, all diagnostics with 
! hourly, daily and monthly time averaging intervals are written out.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCOI_Diagn_WriteOut( am_I_Root, HcoState, WriteAll, &
                                  RC,        PREFIX,   UsePrevTime)
!
! !USES:
!

    USE Ncdf_Mod,      ONLY : NC_Create
    USE Ncdf_Mod,      ONLY : NC_Close
    USE Ncdf_Mod,      ONLY : NC_Var_Def
    USE Ncdf_Mod,      ONLY : NC_Var_Write
    USE HCO_State_Mod, ONLY : HCO_State
    USE JulDay_Mod,    ONLY : JulDay
    USE HCO_Clock_Mod, ONLY : HcoClock_Get
    USE HCO_Clock_Mod, ONLY : HcoClock_GetMinResetFlag
!
! !INPUT PARAMETERS:
!
    LOGICAL,          INTENT(IN   ) :: am_I_Root   ! root CPU?
    TYPE(HCO_State),  POINTER       :: HcoState    ! HEMCO state object 
    LOGICAL,          INTENT(IN   ) :: WriteAll    ! Write all diagnostics? 
    CHARACTER(LEN=*), OPTIONAL      :: PREFIX      ! File prefix
    LOGICAL,          OPTIONAL      :: UsePrevTime ! Use previous time 
!                                                  !  in filename?
! !INPUT/OUTPUT PARAMETERS:
!

    INTEGER,          INTENT(INOUT) :: RC          ! Failure or success
!
! !REVISION HISTORY: 
!  12 Sep 2013 - C. Keller    - Initial version 
!  11 Jun 2014 - R. Yantosca - Cosmetic changes in ProTeX headers
!  11 Jun 2014 - R. Yantosca - Now use F90 freeform indentation
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER                   :: I, CNT, levIdTmp
    REAL(dp)                  :: GMT, JD1, JD1985, JD_DELTA, THISDAY
    REAL(sp)                  :: TMP, JD_DELTA_RND
    INTEGER                   :: YYYY, MM, DD, h, m, s
    REAL(sp), POINTER         :: nctime(:)
    REAL(df), POINTER         :: Arr1D(:) => NULL()
    INTEGER,  POINTER         :: Int1D(:) => NULL()
    REAL(hp), POINTER         :: Arr3D(:,:,:) => NULL()
    REAL(hp), POINTER         :: Arr4D(:,:,:,:) => NULL()
    TYPE(DiagnCont), POINTER  :: ThisDiagn => NULL()
    INTEGER                   :: FLAG
    CHARACTER(LEN=255)        :: ncFile
    CHARACTER(LEN=255)        :: Pfx 
    CHARACTER(LEN=4 )         :: Yrs
    CHARACTER(LEN=2 )         :: Mts, Dys, hrs, mns 
    CHARACTER(LEN=31)         :: timeunit, myName, myUnit
    INTEGER                   :: fId, lonId, latId, levId, TimeId
    INTEGER                   :: VarCt
    INTEGER                   :: nLon, nLat, nLev, nTime 
    INTEGER                   :: Prc
    INTEGER                   :: MinResetFlag, MaxResetFlag
    LOGICAL                   :: EOI, PrevTime
    
    !=================================================================
    ! HCOI_DIAGN_WRITEOUT begins here!
    !=================================================================
  
    ! Init
    RC  = HCO_SUCCESS
    CNT = 0

       ! Inherit data precision from HEMCO
       Prc = HP
 
       ! Check if there is at least one diagnostics to write:
       ! If current time stamp is not at the end of an interval - or
       ! if there is no diagnostics container in the list with a reset
       ! flag smaller or equal to MinResetFlag - there will be no matching
       ! container whatsoever. Can leave right here.
       IF ( WriteAll ) THEN
          MinResetFlag = -1
       ELSE
          MinResetFlag = HcoClock_GetMinResetFlag()
       ENDIF
       MaxResetFlag = Diagn_GetMaxResetFlag()
       IF ( MinResetFlag > MaxResetFlag ) RETURN

       ! Get PrevTime flag from input argument or set to default (=> TRUE)
       IF ( PRESENT(UsePrevTime) ) THEN
          PrevTime = UsePrevTime
       ELSE
          PrevTime = .TRUE.
       ENDIF

       !-----------------------------------------------------------------
       ! Create output file
       !-----------------------------------------------------------------

       ! Get HEMCO time: if this is the last call (WriteAll flag
       ! true), get current time. Use previous time otherwise.
       IF ( .NOT. PrevTime ) THEN
          CALL HcoClock_Get(cYYYY=YYYY,cMM=MM,cDD=DD,cH=h,cM=m,cS=s,RC=RC)
          IF ( RC /= HCO_SUCCESS ) RETURN
          EOI = .FALSE.
       ELSE
          CALL HcoClock_Get(pYYYY=YYYY,pMM=MM,pDD=DD,pH=h,pM=m,pS=s,RC=RC)
          IF ( RC /= HCO_SUCCESS ) RETURN
          EOI = .TRUE.
       ENDIF

       ! Define grid dimensions
       nLon  = HcoState%NX
       nLat  = HcoState%NY
       nLev  = HcoState%NZ
       nTime = 1 

       ! Initialize mirror variables
       allocate(Arr4D(nlon,nlat,nlev,ntime))
       allocate(Arr3D(nlon,nlat,ntime))
       Arr3D = 0.0_hp
       Arr4D = 0.0_hp

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
       ELSE
          CALL Diagn_GetDiagnPrefix( Pfx, RC )
          IF ( RC /= HCO_SUCCESS ) RETURN
       ENDIF
       ncFile = TRIM(Pfx)//'.'//Yrs//Mts//Dys//hrs//mns//'.nc'

       ! testing only
       write(*,*) 'write diagnostics to ', TRIM(ncFile)

       ! Create output file
       CALL NC_CREATE( ncFile, nLon,  nLat,  nLev,  nTime, &
                       fId,    lonId, latId, levId, timeId, VarCt ) 

       !-----------------------------------------------------------------
       ! Write grid dimensions (incl. time) 
       !-----------------------------------------------------------------

       ! Add longitude 
       CALL NC_VAR_DEF ( fId, lonId, -1, -1, -1, &
                         'lon', 'Longitude', 'degrees_east', Prc, VarCt )
       Arr1D => HcoState%Grid%xmid(:,1)
       CALL NC_VAR_WRITE ( fId, 'lon', Arr1D=Arr1D )
       Arr1D => NULL()
       
       ! Add latitude
       CALL NC_VAR_DEF ( fId, -1, latId, -1, -1, &
                         'lat', 'Latitude', 'degrees_north', Prc, VarCt )
       Arr1D => HcoState%Grid%ymid(1,:)
       CALL NC_VAR_WRITE ( fId, 'lat', Arr1D=Arr1D )
       Arr1D => NULL()

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
       ! Write diagnostics 
       !-----------------------------------------------------------------

       ! Loop over all diagnostics in diagnostics list 
       ThisDiagn => NULL()
       DO WHILE ( .TRUE. )

          ! Get next diagnostics in list. This will return the next 
          ! diagnostics container that contains content to be written
          ! out on this time step.
          CALL Diagn_Get ( am_I_Root, HcoState, EOI, ThisDiagn, FLAG, RC )
          IF ( RC /= HCO_SUCCESS ) RETURN 
          IF ( FLAG /= HCO_SUCCESS ) EXIT

          ! Define variable
          myName = ThisDiagn%cName
          myUnit = ThisDiagn%OutUnit
          IF ( ThisDiagn%SpaceDim == 3 ) THEN
             levIdTmp = levId
          ELSE
             levIdTmp = -1
          ENDIF

          CALL NC_VAR_DEF ( fId, lonId, latId, levIdTmp, timeId, &
               TRIM(myName), TRIM(myName), TRIM(myUnit), Prc, VarCt)

          ! Mirror data and write to file. The mirroring is required in
          ! order to add the time dimension. Otherwise, the data would
          ! have no time information!
          IF ( ThisDiagn%SpaceDim == 3 ) THEN
             Arr4D(:,:,:,1) = thisdiagn%Arr3D%val
             CALL NC_VAR_WRITE ( fId, TRIM(myName), Arr4D=Arr4D )
          ELSE
             Arr3D(:,:,1) = thisdiagn%Arr2D%val 
             CALL NC_VAR_WRITE ( fId, TRIM(myName), Arr3D=Arr3D )
          ENDIF
       ENDDO

       !-----------------------------------------------------------------
       ! Close file
       !-----------------------------------------------------------------
       CALL NC_CLOSE ( fId )

       ! Cleanup
       deallocate(Arr3D,Arr4D)
       ThisDiagn => NULL()

    ! Return 
    RC = HCO_SUCCESS

  END SUBROUTINE HCOI_DiagN_WriteOut
!EOC
END MODULE HCOI_GC_Diagn_Mod
#endif

