!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: hco_main_mod
!
! !DESCRIPTION: Module NG\_EMISSIONS\_MOD is the driver routine of the
! new generation emissions module.\\
! It calls all the subroutines to create the emissions linked list
! (EmisLL) and to perform the emissions calculations.\\
! EmisLL is the core of the emissions module. It contains all arrays
! required to perform the emissions calculations at any given time,
! including base emissions, scale factors and masks.\\
! Data is passed to the EmisLL from the temporary data pool TempPool.
! TempPool contains detailed information on all data that has to be read
! at the current time, including netCDF file name, parameter, time
! stamp, etc. Once all information is collected, the arrays are read and
! afterwards passed to EmisLL. All data is regridded onto the simulation
! grid at this point, and base emission units are converted to
! kg/m2/s. Existing data in EmisLL will be overwritten. 
! TempPool is cleared after this step, but can be refilled at a later 
! point.\\
! The content of TempPool is defined depending upon the current time step
! and the simulation settings. For example, EDGAR information will only
! be added to TempPool if EDGAR emissions are enabled and if it's time
! to read them (which is, on the first emissions time step). Similarly,
! EPA information is only set if EPA emissions are switched on, if the
! current node includes NA, and if it's the first emission time step or
! if it's a new month.\\
! The emissions are calculated using a modular approach, i.e through
! multiplication of base emissions and corresponding scale factors (and
! masks), all of which is stored in EmisLL. The emissions are written into 
! the State_Chm object and have units of molec/cm2/s. The emission
! calculation can be performed on a grid resolution that differs from
! the simulation grid. This is particularly interesting if some of the
! emission data has high spatial inhomogeneity that cannot be properly
! resolved on the simulation grid. However, this also means that an
! additional regridding calculation is required on every emissions time
! step.\\ 
! Fields with non-linear emission calculations and specific emissions
! parameterizations (i.e. oceanic emissions) are not calculated through
! EmisLL but in a separate module, which is called after the 'regular'
! emissions calculations. Such online emissions are always added to the
! existing fields.\\
! It should be noted that all emissions are released into the vertical
! leve defined in the original file, i.e. boundary layer mixing is not
! accounted for. PBL mixing is performed in a separated step outside of
! this module.
!\\
!\\
! !INTERFACE: 
!
      MODULE HCO_MAIN_MOD 
! 
! !USES:
!
      USE HCO_ERROR_MOD
      USE HCO_TYPE_MOD,        ONLY : HCO_State 
      USE HCO_EMISLL_MOD,      ONLY : EmisCont
      USE HCO_TIME_MOD,        ONLY : HcoClock

      IMPLICIT NONE

      PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
      PUBLIC :: HCO_INIT 
      PUBLIC :: HCO_RUN 
      PUBLIC :: HCO_GetArray 
      PUBLIC :: HCO_GetState 
      PUBLIC :: HCO_FINAL 
!
! !REVISION HISTORY:
!  27 May 2012 - C. Keller    - Initialization
!  20 Apr 2013 - C. Keller    - Major update. Now use temporary pool for
!                               data reading. 
!  22 Aug 2013 - C. Keller    - Stripped from ng_emissions_mod.F
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !PUBLIC TYPES:
!
      ! (Internal) HEMCO state object
      TYPE(HCO_State), PUBLIC, POINTER    :: HcoStInt  => NULL()

      ! Emissions linked list
      TYPE(EmisCont), PUBLIC, POINTER     :: EmisLL    => NULL()
     
      ! Clock
      TYPE(HcoClock), PUBLIC, POINTER     :: Clock     => NULL() 

      CONTAINS

!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: hco_run 
!
! !DESCRIPTION: Subroutine HCO\_RUN is the HEMCO run routine. 
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE HCO_RUN ( am_I_Root, HcoState, RC ) 
!
! !USES:
!
      USE HCO_TIME_MOD,        ONLY : Set_Clock
      USE HCO_CALC_MOD,        ONLY : Calc_Emis
      USE HCO_READLIST_MOD,    ONLY : ReadData 
      USE HCO_READLIST_MOD,    ONLY : Write_EmisLL 
!
! !INPUT/OUTPUT PARAMETERS:
!
      LOGICAL,         INTENT(IN   )  :: am_I_Root   ! root CPU?
      TYPE(HCO_State), POINTER        :: HcoState    ! External state object
      INTEGER,         INTENT(INOUT)  :: RC          ! Success or failure?
!
! !REVISION HISTORY: 
!  27 May 2012 - C. Keller    - Initialization
!  22 Aug 2013 - C. Keller    - Stripped from ng_emissions_mod.F
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER                    :: ScalIDs(50)

      !=================================================================
      ! HCO_RUN
      !=================================================================

      ! Initialize scale factor IDs 
      ScalIDs(:) = -1

      !=================================================================
      ! 1. Set HEMCO clock 
      !=================================================================

      CALL SET_CLOCK ( Clock, HcoState, RC ) 
      IF ( RC /= HCO_SUCCESS ) RETURN

      !=================================================================
      ! 2. Read all data as defined in ReadList and pass it to the 
      ! Emissions linked list EmisLL.
      ! All data is on the emission grid and in units of mass/area/time. 
      !=================================================================

      ! Read/update data (arrays are written into ReadList)
      CALL ReadData ( am_I_Root, HcoState, Clock, RC )
      IF ( RC /= HCO_SUCCESS ) RETURN

      ! Add/update new data to emissions linked list.
      CALL Write_EmisLL ( am_I_Root, HcoState, EmisLL, Clock, RC ) 
      IF ( RC /= HCO_SUCCESS ) RETURN

      !=================================================================
      ! 3. Calculate the emissions for current time step based on the
      ! content of EmisLL. Emissions are written into the 
      ! HcoState%Emsr3D(TrcID)%Arr3D vectors. 
      ! All emissions in molec/cm2/s.
      !=================================================================

      CALL Calc_Emis ( am_I_Root, HcoState, EmisLL, Clock, RC ) 
      IF ( RC /= HCO_SUCCESS ) RETURN

      ! Leave w/ success
      RC = HCO_SUCCESS

      END SUBROUTINE HCO_RUN
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: hco_getArray 
!
! !DESCRIPTION: Subroutine HCO\_getarray returns the specified array
! of the emissions linked list. 
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE HCO_GetArray ( am_I_Root, ContName, Array, RC )
!
! !USES:
!
      USE HCO_EMISLL_MOD,     ONLY : Find_EmisCont
!
! !ARGUMENTS:
!
      LOGICAL,          INTENT(IN   )  :: am_I_Root      ! root CPU?
      CHARACTER(LEN=*), INTENT(IN   )  :: ContName       ! container name
      REAL*8,           POINTER        :: Array(:,:,:,:) ! output 4D array
      INTEGER,          INTENT(INOUT)  :: RC             ! Failure or success
!
! !REVISION HISTORY: 
!  04 Sep 2013 - C. Keller    - Initialization
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      TYPE(EmisCont), POINTER    :: ThisCont => NULL()
      LOGICAL                    :: FOUND
      CHARACTER(LEN=255)         :: MSG, LOC

      !=================================================================
      ! HCO_GetArray BEGINS HERE
      !=================================================================

      ! Enter
      LOC = 'HCO_GetArray (HCO_MAIN.F90)'

      ! Search for container in emissions linked list
      CALL Find_EmisCont ( TRIM(ContName), EmisLL, FOUND, ThisCont )

      ! Error check
      IF ( .NOT. FOUND ) THEN
         MSG = 'Container not found: ' // TRIM(ContName)
         CALL HCO_ERROR ( MSG, LOC, RC )
         RETURN
      ENDIF

      ! Point to array
      Array => ThisCont%Array

      ! Leave w/ success
      RC = HCO_SUCCESS

      END SUBROUTINE HCO_GetArray 
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: hco_getState 
!
! !DESCRIPTION: Subroutine HCO\_GetState returns the HEMCO state.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE HCO_GetState ( HcoState )
!
! !ARGUMENTS:
!
      TYPE(HCO_State), POINTER  :: HcoState 
!
! !REVISION HISTORY: 
!  04 Sep 2013 - C. Keller    - Initialization
!EOP
!------------------------------------------------------------------------------
!BOC

      !=================================================================
      ! HCO_GETSTATE BEGINS HERE
      !=================================================================

      ! Point to HEMCO state
      HcoState => HcoStInt

      END SUBROUTINE HCO_GetState 
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCO_INIT 
!
! !DESCRIPTION: Subroutine HCO\_INIT initializes the HEMCO derived
! types and arrays. 
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE HCO_INIT( am_I_Root,  nSpecies,    SpecIDs,    &
                             SpecNames,  SpecMW,    EmSpecMW,   &
                             MolecRatio, nLon,      nLat, nLev, & 
                             XMID,       YMID,      XEDGE,      &
                             YEDGE,      YSIN,      AREA_M2,    & 
                             BXHEIGHT_M, OnSimGrid, ConfigFile, &
#if defined(ESMF_)
                             IMPORT,                            &
#endif
                             TS_EMIS,  TS_DYN, RC )
!
! !USES:
!
      USE GIGC_ErrCode_Mod

      USE HCO_TYPE_MOD,        ONLY : Init_HCO_State
      USE HCO_TIME_MOD,        ONLY : Init_tSlc
      USE HCO_TIME_MOD,        ONLY : Init_Clock
      USE HCO_READLIST_MOD,    ONLY : ReadList_Init
      USE HCO_CONFIG_MOD,      ONLY : SetReadList 
      USE HCO_GRID_MOD,        ONLY : HCO_SET_GRID

      ! ESMF only:
#if defined(ESMF_)
      USE ESMF
#endif 
!
! !INPUT/OUTPUT PARAMETERS:
!
      LOGICAL,          INTENT(IN   )  :: am_I_Root  ! root CPU?
      INTEGER,          INTENT(IN   )  :: nSpecies   ! # of tracers
      INTEGER,          INTENT(IN   )  :: SpecIDs(:) ! tracerIDs 
      CHARACTER(LEN=*), INTENT(IN   )  :: SpecNames(:) 
      REAL*8,           INTENT(IN   )  :: SpecMW(:)
      REAL*8,           INTENT(IN   )  :: EmSpecMW(:)
      REAL*8,           INTENT(IN   )  :: MolecRatio(:)
      INTEGER,          INTENT(IN   )  :: nLon       ! # lons on this CPU 
      INTEGER,          INTENT(IN   )  :: nLat       ! # lats on this CPU 
      INTEGER,          INTENT(IN   )  :: nLev       ! # levs on this CPU 
      REAL*8,           INTENT(IN   )  :: XMID(:,:,:)
      REAL*8,           INTENT(IN   )  :: YMID(:,:,:)
      REAL*8,           INTENT(IN   )  :: XEDGE(:,:,:)
      REAL*8,           INTENT(IN   )  :: YEDGE(:,:,:)
      REAL*8,           INTENT(IN   )  :: YSIN(:,:,:)
      REAL*8,           INTENT(IN   )  :: AREA_M2(:,:,:)
      REAL*8,           INTENT(IN   )  :: BXHEIGHT_M(:,:,:)
      LOGICAL,          INTENT(IN   )  :: OnSimGrid 
      CHARACTER(LEN=*), INTENT(IN   )  :: ConfigFile 
#if defined(ESMF_)
      TYPE(ESMF_State), INTENT(IN   )  :: IMPORT     ! ESMF Import obj.
#endif
      REAL*8,           INTENT(IN   )  :: TS_EMIS    ! Emission timestep [s] 
      REAL*8,           INTENT(IN   )  :: TS_DYN     ! Dynamics timestep [s] 
      INTEGER,          INTENT(INOUT)  :: RC         ! Failure or success
!
! !REVISION HISTORY: 
!  27 May 2012 - C. Keller    - Initialization
!  22 Aug 2013 - C. Keller    - Stripped from ng_emissions_mod.F
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      REAL*8   :: XRNG(2), YRNG(2)

      !=================================================================
      ! HCO_INIT begins here!
      !=================================================================

      ! Initialize internal HEMCO state
      CALL Init_HCO_State( am_I_Root  = am_I_Root,   &
                           nSpecies   = nSpecies,    &
                           SpecIDs    = SpecIDs,     &
                           SpecNames  = SpecNames,   &
                           SpecMW     = SpecMW,      &
                           EmSpecMW   = EmSpecMW,    &
                           MolecRatio = MolecRatio,  &
                           ConfigFile = ConfigFile,  &
#if defined(ESMF_)
                           IMPORT     = IMPORT,      &
#endif
                           TS_EMIS    = TS_EMIS,     &
                           TS_DYN     = TS_DYN,      &
                           HcoState   = HcoStInt,    &
                           RC         = RC            )
      IF ( RC /= HCO_SUCCESS ) RETURN

      ! Initialize the Clock
      CALL Init_Clock ( Clock, RC )
      IF ( RC /= HCO_SUCCESS ) RETURN

      ! Initialize reading list (before config file!!)
      CALL ReadList_Init 

      ! Read config file with information of inventories and pass
      ! information to ReadList.
      ! Pass grid boundaries of this CPU to exclude data outside of the
      ! grid dimension.
      XRNG(1) = MINVAL(XEDGE)
      XRNG(2) = MAXVAL(XEDGE)
      YRNG(1) = MINVAL(YEDGE)
      YRNG(2) = MAXVAL(YEDGE)
      CALL SetReadList ( am_I_Root, HcoStInt, XRNG, YRNG, RC )
      IF ( RC /= HCO_SUCCESS ) RETURN

      ! Initialize the HEMCO grid
      CALL HCO_SET_GRID ( HcoState   = HcoStInt,    & 
                            I          = nLon,        &
                            J          = nLat,        &
                            L          = nLev,        &
                            XMID       = XMID,        &
                            YMID       = YMID,        &
                            XEDGE      = XEDGE,       &
                            YSIN       = YSIN,        &
                            AREA_M2    = AREA_M2,     &
                            BXHEIGHT_M = BXHEIGHT_M,  &
                            OnSimGrid  = OnSimGrid,   &
                            RC         = RC            )
      IF ( RC /= HCO_SUCCESS ) RETURN

      ! Initialize time slice pointers 
      CALL Init_tSlc ( HcoStInt%ISIZE, RC )
      IF ( RC /= HCO_SUCCESS ) RETURN

      ! Leave w/ success
      RC = HCO_SUCCESS

      END SUBROUTINE HCO_INIT 
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: hco_final 
!
! !DESCRIPTION: Subroutine HCO\_FINAL finalizes HEMCO. 
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE HCO_FINAL 
!
! !REVISION HISTORY: 
!  27 May 2012 - C. Keller    - Initialization
!  22 Aug 2013 - C. Keller    - Stripped from ng_emissions_mod.F
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!

      !=================================================================
      ! HCO_FINAL 
      !=================================================================

      CALL HCO_CLEANUP 

      END SUBROUTINE HCO_FINAL 
!EOC
      END MODULE HCO_MAIN_MOD 
