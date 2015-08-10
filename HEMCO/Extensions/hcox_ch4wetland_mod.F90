!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: hcox_ch4wetland_mod.F90
!
! !DESCRIPTION: Module HCOX\_CH4Wetland\_Mod contains routines to 
! calculate methane emissions (including rice) from wetlands. This code
! is adapted from the GEOS-Chem CH4 offline simulation.
!\\
!\\ 
! This is a HEMCO extension module that uses many of the HEMCO core
! utilities.
!\\
!\\
! This code can be used to calculate emissions from wetlands, from rice, or
! both. Both sources can be enabled/disabled in the HEMCO configuration file.
!\\
!\\
! The accepted species are: CH4 or CH4\_tot for total CH4, CH4\_ri for rice
! only, and CH4\_we for wetland only. Any combination of these species IDs
! can be provided (but at least one valid species ID needs to be given).
!\\
!\\
! References:
! \begin{itemize}
! \item Pickett-Heaps CA, Jacob DJ, Wecht KJ, et al. Magnitude and seasonality 
! of wetland methane emissions from the Hudson Bay Lowlands (Canada). ACP, 11, 
! 3773-3779, 2011.
! \end{itemize}
! !INTERFACE: 
!
MODULE HCOX_CH4WETLAND_Mod
!
! !USES:
! 
  USE HCO_Error_MOD
  USE HCO_Diagn_MOD
  USE HCO_State_MOD,  ONLY : HCO_State
  USE HCOX_State_MOD, ONLY : Ext_State

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC  :: HCOX_CH4WETLAND_INIT
  PUBLIC  :: HCOX_CH4WETLAND_RUN
  PUBLIC  :: HCOX_CH4WETLAND_FINAL
!
! !PRIVATE MEMBER FUNCTIONS:
!
  PRIVATE :: WETLAND_EMIS 
  PRIVATE :: RICE_EMIS 
!
! !REVISION HISTORY:
!  11 Sep 2014 - C. Keller   - Initial version
!  01 Oct 2013 - C. Keller   - Now a HEMCO extension module 
!  11 Dec 2013 - C. Keller   - Now define container name during initialization
!  01 Jul 2014 - R. Yantosca - Now use F90 free-format indentation
!  01 Jul 2014 - R. Yantosca - Cosmetic changes in ProTeX headers
!EOP
!------------------------------------------------------------------------------
!
! !PRIVATE VARIABLES:
!
  ! Module variables related to extension settings
  INTEGER                     :: ExtNr
  INTEGER                     :: IDTtot
  LOGICAL                     :: DoWetland 
  LOGICAL                     :: DoRice

  ! Pointers to data read through configuration file
  REAL(sp), POINTER           :: RICE        (:,:) => NULL()
  REAL(sp), POINTER           :: GWET_ANNUAL (:,:) => NULL()
  REAL(sp), POINTER           :: GWET_MONTHLY(:,:) => NULL()
  REAL(sp), POINTER           :: WETFRAC     (:,:) => NULL()
  REAL(sp), POINTER           :: LITTER_C    (:,:) => NULL()
  REAL(sp), POINTER           :: SOIL_C      (:,:) => NULL()
  REAL(sp), POINTER           :: MEAN_T      (:,:) => NULL()

CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCOX_CH4WETLAND_Run
!
! !DESCRIPTION: Subroutine HcoX\_CH4WETLAND\_Run is the run routine to 
! calculate oceanic emissions for the current time step.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCOX_CH4WETLAND_Run( am_I_Root, ExtState, HcoState, RC )
!
! !USES:
!
    USE HCO_EMISLIST_MOD, ONLY : HCO_GetPtr 
    USE HCO_FLUXARR_MOD,  ONLY : HCO_EmisAdd
!
! !INPUT PARAMETERS:
!
    LOGICAL,         INTENT(IN   ) :: am_I_Root  ! root CPU?
    TYPE(HCO_State), POINTER       :: HcoState   ! Output obj
    TYPE(Ext_State), POINTER       :: ExtState  ! Module options  
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,         INTENT(INOUT) :: RC         ! Success or failure?
!
! !REVISION HISTORY: 
!  11 Sep 2014 - C. Keller - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    LOGICAL, SAVE            :: FIRST   = .TRUE.
    LOGICAL, SAVE            :: DoDiagn = .FALSE.

    ! Array holding the CH4 emissions
    REAL(hp), TARGET         :: CH4wtl(HcoState%NX,HcoState%NY)
    REAL(hp), TARGET         :: CH4rce(HcoState%NX,HcoState%NY)
    TYPE(DiagnCont), POINTER :: TmpCnt     => NULL()

    ! Name of the manual diagnostic container
    CHARACTER(LEN=31),  PARAMETER :: DiagnWtl = 'CH4_WETLAND'
    CHARACTER(LEN=31),  PARAMETER :: DiagnRce = 'CH4_RICE'

    CHARACTER(LEN=255), PARAMETER :: LOC = &
       'HCOX_CH4WETLAND_Run (hcox_ch4wetland_mod.F90)'

    !=================================================================
    ! HCOX_CH4WETLAND_Run begins here!
    !=================================================================

    ! Enter
    CALL HCO_ENTER( LOC, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Return if extension disabled 
    IF ( .NOT. ExtState%Wetland_CH4 ) RETURN

    ! ---------------------------------------------------------------
    ! On first call, get pointers to data and check if manual 
    ! diagnostics will be used 
    ! ---------------------------------------------------------------
    IF ( FIRST ) THEN
    
       IF ( DoWetland ) THEN
          CALL HCO_GetPtr( am_I_Root, 'CH4_WETFRAC',  WETFRAC,  RC )
          IF ( RC /= HCO_SUCCESS ) RETURN
          CALL HCO_GetPtr( am_I_Root, 'CH4_LITTER_C', LITTER_C, RC )
          IF ( RC /= HCO_SUCCESS ) RETURN
          CALL HCO_GetPtr( am_I_Root, 'CH4_SOIL_C',   SOIL_C,   RC )
          IF ( RC /= HCO_SUCCESS ) RETURN
          CALL HCO_GetPtr( am_I_Root, 'CH4_MEAN_T',   MEAN_T,   RC )
          IF ( RC /= HCO_SUCCESS ) RETURN
          IF ( .NOT. DoDiagn ) THEN
             CALL DiagnCont_Find ( -1,-1,-1,-1,-1, TRIM(DiagnWtl), 0, DoDiagn, TmpCnt )
             TmpCnt => NULL()
          ENDIF
       ENDIF
 
       ! Fields required by rice emissions
       IF ( DoRice ) THEN
          CALL HCO_GetPtr( am_I_Root, 'CH4_RICE',    RICE,         RC )
          IF ( RC /= HCO_SUCCESS ) RETURN
          CALL HCO_GetPtr( am_I_Root, 'CH4_GWET_YR', GWET_ANNUAL,  RC )
          IF ( RC /= HCO_SUCCESS ) RETURN
          CALL HCO_GetPtr( am_I_Root, 'CH4_GWET_MT', GWET_MONTHLY, RC )
          IF ( RC /= HCO_SUCCESS ) RETURN
          IF ( .NOT. DoDiagn ) THEN
             CALL DiagnCont_Find ( -1,-1,-1,-1,-1, TRIM(DiagnRce), 0, DoDiagn, TmpCnt )
             TmpCnt => NULL()
          ENDIF
       ENDIF
    ENDIF

    ! ---------------------------------------------------------------
    ! Calculate emissions
    ! ---------------------------------------------------------------
    CH4rce = 0.0_hp
    CH4wtl = 0.0_hp

    ! If turned on, calculate wetland emissions
    IF ( DoWetland ) THEN
       CALL WETLAND_EMIS ( am_I_Root, HcoState, ExtState, CH4wtl, RC )
    ENDIF

    ! If turned on, calculate rice emissions. Previously calculated 
    ! wetland emissions may be subtracted!
    IF ( DoRice ) THEN
       CALL RICE_EMIS ( am_I_Root, HcoState, ExtState, CH4rce, RC )
    ENDIF

    ! ---------------------------------------------------------------
    ! Pass emissions to HEMCO state and update diagnostics 
    ! ---------------------------------------------------------------

    ! Adjust for double counting
    WHERE ( CH4wtl > CH4rce )
       CH4wtl = CH4wtl - CH4rce
    ENDWHERE

    ! Eventually update manual diagnostics
    IF ( DoDiagn ) THEN
       IF ( DoWetland ) THEN
          CALL Diagn_Update( am_I_Root, ExtNr=ExtNr, &
                             cName=TRIM(DiagnWtl), Array2D=CH4wtl, RC=RC )
          IF ( RC /= HCO_SUCCESS ) RETURN
       ENDIF
       IF ( DoRice ) THEN
          CALL Diagn_Update( am_I_Root, ExtNr=ExtNr, &
                             cName=TRIM(DiagnRce), Array2D=CH4rce, RC=RC )
          IF ( RC /= HCO_SUCCESS ) RETURN
       ENDIF
    ENDIF

    ! Total CH4 emissions [kg/m2/s]
    CH4wtl = CH4wtl + CH4rce 

    ! Set flux in HEMCO object [kg/m2/s]
    CALL HCO_EmisAdd ( am_I_Root, HcoState, CH4wtl, IDTtot, RC, ExtNr=ExtNr )
    IF ( RC /= HCO_SUCCESS ) THEN
       CALL HCO_ERROR( 'HCO_EmisAdd error: CH4wtl', RC )
       RETURN 
    ENDIF
 
    ! Leave w/ success
    CALL HCO_LEAVE ( RC ) 

  END SUBROUTINE HCOX_CH4WETLAND_Run
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: WETLAND_EMIS 
!
! !DESCRIPTION: Subroutine WETLAND\_EMIS is the driver routine for the CH4
! wetland emissions. It calculates wetland emissions and writes them into 
! the passed array CH4wtl. 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE WETLAND_EMIS ( am_I_Root, HcoState, ExtState, CH4wtl, RC )
!
! !USES:
!
    USE HCO_CLOCK_MOD,      ONLY : HcoClock_Get
!
! !INPUT PARAMETERS:
!
    LOGICAL,         INTENT(IN   ) :: am_I_Root  ! root CPU?
    TYPE(HCO_State), POINTER       :: HcoState   ! Output obj
    TYPE(Ext_State), POINTER       :: ExtState  ! Module options  
!
! !INPUT/OUTPUT PARAMETERS:
!
    REAL(hp),        INTENT(INOUT) :: CH4wtl(HcoState%NX,HcoState%NY) ! CH4 emis
    INTEGER,         INTENT(INOUT) :: RC         ! Success or failure?
!
! !REVISION HISTORY: 
!  (1 ) Adapted by Jérôme Drevet (3/06) from the BIOME-TG Wetland-Methane
!       scheme provided by Jed O. Kaplan.
!  (2 ) CH4 Emissions from Wetland depend on:
!               a - Soil Carbon content.
!               b - Vegetation type
!               c - Wetland area (%)
!               d - Soil moisture.
!       a, b, c are taken from the LPJ, a vegetation model. Data are provided 
!       by J.O.Kaplan. Soil moisture is read from GEOS Met input files.  
!  (3 ) Corrected order of DO loops (bmy, 10/1/09)
!  08 Feb 2012 - R. Yantosca - Treat GEOS-5.7.x in the same way as MERRA
!  01 Mar 2012 - R. Yantosca - Now use GET_AREA_M2(I,J,L) from grid_mod.F90
!  07 Mar 2012 - M. Payer    - Added ProTeX headers
!  09 Nov 2012 - M. Payer    - Replaced all met field arrays with State_Met
!                              derived type object
!  26 Sep 2013 - R. Yantosca - Renamed GEOS_57 Cpp switch to GEOS_FP
!  23 Jan 2014 - M. Sulprizio- Now zero wetland emissions if snow covers the
!                              ground. Also updated MOIST_SCALE and EMIT_FACT.
!                              (K. Wecht, C. Pickett-Heaps)
!  12 Feb 2014 - K. Wecht    - Updated for 0.25 x 0.3125 NA grid
!  09 Apr 2014 - R. Yantosca - Bug fix, extend #ifdef for MERRA met fields
!  11 Sep 2014 - C. Keller   - Now a HEMCO extension
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!

    INTEGER             :: I, J, LMD
    REAL(dp)            :: SECINMONTH
    REAL(dp)            :: REALWET
    REAL(dp)            :: EFF_GWET
    REAL(dp)            :: litterfast
    REAL(dp)            :: litterslow
    REAL(dp)            :: soilfast
    REAL(dp)            :: soilslow
    REAL(dp)            :: HETEROR 
    REAL(dp)            :: F_TEMP
    REAL(dp)            :: TROPICNESS
    REAL(dp)            :: EMIT_TROPIC 
    REAL(dp)            :: EMIT_TEMPER

    ! -------------
    ! PARAMETER:

!    ! Molecular weight of carbon (g/mol)
    REAL(dp), PARAMETER           :: MWC = 12.0_dp
 
    ! Max allowable snowdepth for emissions to take place [meters]
    REAL(dp), PARAMETER           :: MAX_SNOWDP = 0.0001_dp

    ! Comments from Jed Kaplan and Jerome Drevet below:
    !    (moist_scale can be between 0.07 and 0.14)
    !    (emit_fact can be between 0.001 and 0.005)
    ! Jerome's paper uses MOIST_SCALE = 0.19 (kjw, 6/9/09)
    !
    ! Jed and Jerome tune parameter values to the wetland emission
    !    routine using GEOS-4 met variables. I have tuned the model
    !    using GEOS-5 met variables and found that EMIT_FACT=.01 
    !    and MOIST_SCALE=.1 give a reasonable balance between
    !    extratropical and tropical emissions, kjw (7/20/10).
    ! 
    ! UPDATE 12/13/2011
    !    Carbon litter and soil maps used prior to GC v9-01-03 were
    !    incorrect. Carbon litter had wrong spatial distribution, and
    !    carbon soil values were too high and may also have had the 
    !    wrong spatial distribution.
    !    Updated carbon litter, carbon soil, and wetland fraction maps
    !    are included in v9-01-03. MOIST_SCALE and EMIT_FACT have been 
    !    tuned to maintain annual global average methane emissions 
    !    ~ 166.8 Tg/y. Boreal wetland emissions are ~ 28.9 Tg/y, 
    !    consistent with estimates described in Pickett-Heaps 2011,
    !    3rd pararaph of introduction. New parameter values:
    !      MOIST_SCALE = 0.19
    !      EMIT_FACT   = 0.023
    !    (kjw, 12/13/2011)
    !MOIST_SCALE=0.19d0
    !EMIT_FACT = .023d0
    !
    ! UPDATE 1/16/2014
    !    Use of the new carbon litter and soil maps changed global
    !    and regional emission totals. Here, I change MOIST_SCALE
    !    and EMIT_FACT to 1) match HBL emissions constrained by 
    !    surface and aircraft data in Pickett-Heaps et al. 2011. 
    !    and 2) preserve global total emissions (~165 Tg/y) used
    !    in comparison to NOAA GMD, SCIAMACHY, and INTEX-A in 
    !    Wecht et al. (2014).
    REAL(dp), PARAMETER :: MOIST_SCALE = 0.205_dp   ! tropical emissions
    REAL(dp), PARAMETER :: EMIT_FACT   = 0.018_dp   ! extratropical emissions

    CHARACTER(LEN=255), PARAMETER :: LOC = &
       'WETLAND_EMIS (hcox_ch4wetland_mod.F90)'


    !=================================================================
    ! WETLAND_EMIS begins here!
    !=================================================================

    ! Enter
    CALL HCO_ENTER( LOC, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

!$OMP PARALLEL DO                                                         &
!$OMP DEFAULT( SHARED )                                                   &
!$OMP PRIVATE( I, J, REALWET, F_TEMP, HETEROR, litterfast, litterslow )   &
!$OMP PRIVATE( soilfast, soilslow, TROPICNESS, EFF_GWET               )   &
!$OMP PRIVATE( EMIT_TROPIC, EMIT_TEMPER                               )   &
!$OMP SCHEDULE( DYNAMIC )
    DO J = 1, HcoState%NY 
    DO I = 1, HcoState%NX

       !===================================================================
       ! Calculate inundated fraction
       !
       ! REALWET calculation is based on maximum inundatable area (WETFRAC)
       ! and top soil moisture information
       !
       ! NOTE: WLI (land/water/ice flag) definition has changed between
       !   GEOS4 and GEOS5.  This contributes to the variance between GEOS4
       !   and GEOS5 wetland emissions.  Below is Jerome Drevet's and Jed
       !   Kaplan's original calculation of REALWET using GEOS4 and a
       !   modified calculation using GEOS5.
       !                                                     (kjw, 6/10/09)
       !===================================================================

       ! Init
       REALWET = 0.0_hp
   
       ! GEOS5 Calculation of inundated fraction
       ! We don't want emissions in frozen or snow-covered regions
       IF ( ExtState%TSKIN%Arr%Val(I,J) > 273d0      ) THEN
          IF ( ExtState%SNODP%Arr%Val(I,J) < MAX_SNOWDP ) THEN !SNOW DEPTH

             ! GEOS4 calculation of inundated fraction
#if   defined( GEOS_4 )
             ! We want emissions from land boxes only
             IF ( ExtState%WLI%Arr%Val(I,J) == 1 ) THEN

                ! If wetness>0.1, the wetland fraction is equal
                ! to the maximal potential wetland fraction
                IF (ExtState%GWETTOP%Arr%Val(I,J) > 0.1d0) THEN
                   REALWET = WETFRAC(I,J) / 100.0_hp
                ELSE
                   REALWET = 0.d0
                ENDIF
             ENDIF

#elif defined( GEOS_5 ) || defined( MERRA ) || defined( GEOS_FP )
             ! We want emissions from any box that contains some land
             ! FRLAND is fraction of grid box that is land
             IF ( ExtState%FRLAND%Arr%Val(I,J) > 0) THEN 

                ! Actual wetness of land /= GWETTOP because GWETTOP includes
                ! wetness in lakes, ocean, and ice.  Below is a scheme to
                ! calculate effective GWETTOP of the land fraction
                EFF_GWET = ( ExtState%GWETTOP%Arr%Val(I,J)       -   & 
                             ( ExtState%FROCEAN%Arr%Val(I,J)     +   &
                               ExtState%FRLAKE%Arr%Val(I,J)      +   &
                               ExtState%FRLANDIC%Arr%Val(I,J)      ) &
                           ) / ExtState%FRLAND%Arr%Val(I,J)

                ! Catch for negative EFF_GWET
                IF ( EFF_GWET < 0d0 ) THEN
                   EFF_GWET = 0d0
                ENDIF

                ! If wetness>0.1, the wetland fraction is equal
                ! to the maximal potential wetland fraction
                IF (EFF_GWET > 0.1d0) THEN
                   REALWET = WETFRAC(I,J) / 100.0_hp
                ELSE
                   REALWET = 0.0_dp
                ENDIF
             ENDIF
#endif

          ENDIF
       ENDIF

! TODO: need to add diagnostics here
!-------------------------------------------------------------------------
!      ! Update Wetland Fraction Diagnostic
!      GM = GET_MONTH()
!      IF ( ND60 > 0 ) THEN
!	 AD60(:,:) = AD60(:,:) + REALWET(:,:)/(24d0*MONTHDATES(GM))
!      ENDIF
!-------------------------------------------------------------------------

       !===================================================================
       ! Calculate CH4 emissions!
       !===================================================================
 
       IF ( ExtState%TSKIN%Arr%Val(I,J) < 233.d0 ) THEN
          F_TEMP = 0d0
       ELSE
          F_TEMP = exp(308.56d0*(1.0d0/56.02d0-&
                   1.0d0/(ExtState%TSKIN%Arr%Val(I,J)-227.13d0))) !Lloyd & Taylor 1994
       ENDIF
 
       ! Calculate Heterotrophic respiration
       litterfast = 0.985d0 * LITTER_C(i,j)
       litterslow = 0.015d0 * LITTER_C(i,j)
       soilfast =  0.985d0 * SOIL_C(i,j)
       soilslow =  0.015d0 * SOIL_C(i,j)

       ! The division by 12 seems to convert kgC/m2/yr to kgC/m2/mt. Since
       ! HEMCO input data is in kgC/m2/s, no need to do this anymore!
       HETEROR = 1d3* F_TEMP * ( litterfast*0.3d0   &
                               + litterslow*0.05d0  &
                               + soilfast*0.03d0    &
                               + soilslow*0.001d0 ) * 0.34d0 ! / 12d0
 
       ! Calculate "tropicness" of each box
       TROPICNESS = exp((MEAN_T(I,J) - 303.15d0) / 8d0)
       IF ( TROPICNESS < 0d0 ) THEN
          TROPICNESS = 0d0
       ENDIF
       IF ( TROPICNESS > 1d0 ) THEN
          TROPICNESS = 1d0
       ENDIF
         
       EMIT_TROPIC = 0.0d0
       EMIT_TEMPER = 0.0d0

!         EMIT_TROPIC = EMIT_TROPIC + HETEROR * MOIST_SCALE 
!     &                * REALWET(I,J) 
!
!         EMIT_TEMPER = EMIT_TEMPER + HETEROR * EMIT_FACT 
!     &                * REALWET(I,J) 

       EMIT_TROPIC = HETEROR * MOIST_SCALE * REALWET
       EMIT_TEMPER = HETEROR * EMIT_FACT   * REALWET

       ! kg/m2/s
       ! Note: CH4 emissions are calculated by scaling the soil carbon data,
       ! which appears to be in kgC. So I would think that we need to adjust
       ! by 16/12, but the original code does not seem to be doing this...
       CH4wtl(I,J) = TROPICNESS * EMIT_TROPIC + (1-TROPICNESS) * EMIT_TEMPER
       IF (CH4wtl(I,J) < 0.0_hp ) CH4wtl(I,J)=0.0_hp

    ENDDO
    ENDDO
!$OMP END PARALLEL DO

    ! Convert kgC to kgCH4 (needed?! See note above)
!    CH4wtl = CH4wtl / MWC * HcoState%Spc(IDTtot)%MW_g

    ! Leave w/ success
    CALL HCO_LEAVE ( RC ) 

  END SUBROUTINE WETLAND_EMIS 
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: RICE_EMIS 
!
! !DESCRIPTION: Subroutine RICE_EMIS is the driver routine for the CH4
! rice emissions. It calculates rice emissions and writes them into the
! passed array CH4.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE RICE_EMIS ( am_I_Root, HcoState, ExtState, CH4rce, RC )
!
! !USES:
!
!
! !INPUT PARAMETERS:
!
    LOGICAL,         INTENT(IN   ) :: am_I_Root  ! root CPU?
    TYPE(HCO_State), POINTER       :: HcoState   ! Output obj
    TYPE(Ext_State), POINTER       :: ExtState  ! Module options  
!
! !INPUT/OUTPUT PARAMETERS:
!
    REAL(hp),        INTENT(INOUT) :: CH4rce(HcoState%NX,HcoState%NY) ! CH4 emis
    INTEGER,         INTENT(INOUT) :: RC         ! Success or failure?
!
! !REMARKS:
!  Rice Emissions are scaled to GEOS soil wetness.  Scaling sceme developed
!     and implemented by Jerome Drevet.
!  Wetland emissions are modified by the presence of rice emissions.  Sceme
!     developed by Jerome Drevet.
!
! !REVISION HISTORY:
!  (1 ) CH4 emissions from rice calculated with a routine created by Jerome
!       Drevet.  Adapted as its own subroutine by Kevin Wecht (6/03/09)
!  (2 ) Corrected ordering of DO loops (bmy, 10/1/09)
!  07 Mar 2012 - M. Payer    - Added ProTeX headers
!  25 Mar 2013 - R. Yantosca - Now accept am_I_Root, Input_Opt, State_Chm, RC
!  09 Apr 2014 - R. Yantosca - Bug fix, extend #ifdef for MERRA met fields
!  11 Sep 2014 - C. Keller   - Now a HEMCO extension
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER             :: I, J
    REAL(dp)            :: wet_ratio

    CHARACTER(LEN=255), PARAMETER :: LOC = &
       'RICE_EMIS (hcox_ch4wetland_mod.F90)'

    !=================================================================
    ! RICE_EMIS begins here!
    !=================================================================

    ! Enter
    CALL HCO_ENTER( LOC, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    !scale rice emissions (by Jerome Drevet)
!$OMP PARALLEL DO                 &
!$OMP DEFAULT( SHARED )           &
!$OMP PRIVATE( I, J, wet_ratio )  &
!$OMP SCHEDULE( DYNAMIC )
    DO J = 1, HcoState%NY 
    DO I = 1, HcoState%NX
       wet_ratio = GWET_MONTHLY(I,J)/GWET_ANNUAL(I,J)-1.0_hp
       wet_ratio = wet_ratio * 2.0_dp
       wet_ratio = wet_ratio +1.0_dp
       IF (wet_ratio < 0) wet_ratio = 0.0_dp
       CH4rce(I,J) = RICE(I,J) * wet_ratio
    ENDDO
    ENDDO
!$OMP END PARALLEL DO

    ! Leave w/ success
    CALL HCO_LEAVE ( RC ) 

  END SUBROUTINE RICE_EMIS 
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCOX_CH4WETLAND_INIT
!
! !DESCRIPTION: Subroutine HCOX\_CH4WETLAND\_INIT initializes all module
! variables. 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCOX_CH4WETLAND_INIT( am_I_Root, HcoState, ExtName, ExtState, RC )
!
! !USES:
!
    USE HCO_ExtList_Mod,        ONLY : GetExtNr
    USE HCO_ExtList_Mod,        ONLY : GetExtOpt
    USE HCO_STATE_MOD,          ONLY : HCO_GetExtHcoID
!
! !INPUT PARAMETERS:
!
    LOGICAL,          INTENT(IN   )  :: am_I_Root    ! root CPU?
    TYPE(HCO_State),  POINTER        :: HcoState     ! Hemco State obj.
    CHARACTER(LEN=*), INTENT(IN   )  :: ExtName      ! Extension name
    TYPE(Ext_State),  POINTER        :: ExtState     ! Ext. obj.
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(INOUT)  :: RC           ! Return status
!
! !REVISION HISTORY:
!  11 Sep 2014 - C. Keller - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES
!
    ! Scalars
    INTEGER                        :: I, J, nSpc
    CHARACTER(LEN=255)             :: NAME_OC, MSG, ERR

    ! Arrays
    INTEGER,           ALLOCATABLE :: HcoIDs(:)
    CHARACTER(LEN=31), ALLOCATABLE :: SpcNames(:)

    !=================================================================
    ! HCOX_CH4WETLAND_Init begins here!
    !=================================================================

    ! Extension Nr.
    ExtNr = GetExtNr( TRIM(ExtName) )
    IF ( ExtNr <= 0 ) RETURN
 
    ! Enter 
    CALL HCO_ENTER ( 'HCOX_CH4WETLAND_Init (hcox_wetlands_ch4_mod.F90)', RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! HEMCO species IDs of species names defined in config. file 
    CALL HCO_GetExtHcoID( HcoState, ExtNr, HcoIDs, SpcNames, nSpc, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Only one species name must be given!
    IF ( nSpc /= 1 ) THEN
       MSG = 'CH4 wetland extension accepts only one species!'
       CALL HCO_ERROR ( MSG, RC ) 
       RETURN
    ENDIF
    IDTtot = HcoIDs(1)

    ! Make sure species ID is valid 
    IF ( IDTtot < 1 ) THEN
       MSG = 'Invalid CH4 species: ' // TRIM(SpcNames(1))
       CALL HCO_ERROR ( MSG, RC )
       RETURN
    ENDIF

    ! Check if wetlands and/or rice emissions shall be used and save
    ! in module variables DoWetland and DoRice.
    CALL GetExtOpt ( ExtNr, 'Wetlands', OptValBool=DoWetland, RC=RC )
    IF ( RC /= HCO_SUCCESS ) RETURN
    CALL GetExtOpt ( ExtNr, 'Rice',     OptValBool=DoRice,    RC=RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Make sure at least one source is used.
    IF ( .NOT. DoWetland .AND. .NOT. DoRice ) THEN
       MSG = 'Wetlands and rice emissions are both turned off!'
       CALL HCO_ERROR ( MSG, RC )
       RETURN
    ENDIF

    ! Verbose mode
    IF ( am_I_Root ) THEN
       MSG = 'Use wetland flux emissions (extension module)'
       CALL HCO_MSG( MSG,SEP1='-' )
       WRITE(MSG,*) 'CH4 species ID       : ', IDTtot
       CALL HCO_MSG( MSG )
       WRITE(MSG,*) 'Use wetlands         : ', DoWetland
       CALL HCO_MSG( MSG )
       WRITE(MSG,*) 'Use rice             : ', DoRice
       CALL HCO_MSG( MSG )
    ENDIF

    ! Register all required met fields
    ExtState%TSKIN%DoUse    = .TRUE.
    ExtState%SNODP%DoUse    = .TRUE.
    ExtState%WLI%DoUse      = .TRUE.
    ExtState%GWETTOP%DoUse  = .TRUE.
    ExtState%FRLAND%DoUse   = .TRUE.
    ExtState%FROCEAN%DoUse  = .TRUE.
    ExtState%FRLAKE%DoUse   = .TRUE.
    ExtState%FRLANDIC%DoUse = .TRUE.

    ! Enable extensions
    ExtState%Wetland_CH4 = .TRUE.

    ! Return w/ success
    IF ( ALLOCATED(HcoIDs  ) ) DEALLOCATE(HcoIDs  )
    IF ( ALLOCATED(SpcNames) ) DEALLOCATE(SpcNames)
    CALL HCO_LEAVE ( RC )

  END SUBROUTINE HCOX_CH4WETLAND_INIT
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCOX_CH4WETLAND_Final
!
! !DESCRIPTION: Subroutine HCOX\_CH4WETLAND\_Final deallocates 
!  all module arrays.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCOX_CH4WETLAND_Final()
!
! !REVISION HISTORY:
!  11 Sep 2014 - C. Keller - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
    !=================================================================
    ! HCOX_CH4WETLAND_Final begins here!
    !=================================================================

  END SUBROUTINE HCOX_CH4WETLAND_Final
!EOC
END MODULE HCOX_CH4WETLAND_Mod
!EOM
