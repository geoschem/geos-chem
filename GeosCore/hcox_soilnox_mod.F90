!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: hcox_soilnox_mod
!
! !DESCRIPTION: Module HCOX\_SOILNOX\_MOD contains routines to 
!  compute soil NOx emissions. 
!
! This is a HEMCO extension module that uses many of the HEMCO core
! utilities.
!
!  Original codes from:
!    HARVARD ATMOSPHERIC CHEMISTRY MODELING GROUP 
!    MODULE FOR SOIL NOX EMISSIONS           
!    by Yuhang Wang, Gerry Gardner, and Prof. Daniel Jacob 
!  Updated model code:
!    by  Rynda Hudman, Neil Moore, Randall Martin, and Bram Maasakkers
!
! !REMARKS:
!  The soil NOx code has been updated from the original implementation
!  of Yienger & Levy [1995] from  Wang et al., [1998] as summarized below.
!  
!  Old:
!  ENOx   = f( T, biome, w/d)  x Pulse(precip) x canopy uptake + FERT
! 
!  New:
!  ENOx   = f( T, biome, WFPS, Fert)  x Pulse(dryspell) x canopy uptake 
! 
!  1 - Update moisture treatment: soil moisture as a continuous variable 
!  using WFPS rather than discrete wet/dry states and purely exponential 
!  T impact (impact = -1. Tg N/yr)
! 
!  2 - Update to Fertilizer:  new fertilizer maps including chemical and 
!  manure fertilizer from Potter et al., [2010] distributed using MODIS EVI 
!  seasonality, online-N deposition as a fertilizer source, and N-fertilizer
!  source subject to T, WFPS, and pulsing like other N (impact = +1.3 Tg N/yr)
! 
!  3- Update Pulsing Scheme: Yan et al., [2005] (shorter, stronger pulses) 
!  (impact = +1. Tg N/yr). Also added restart file containing dry spell 
!  information to properly account for dry spell length in continuing runs. 
! 
!  References:
!  ============================================================================
!  (1 ) Wang, Y., D.J. Jacob, and J.A. Logan,  Global simulation of 
!        tropospheric O3-NOx-hydrocarbon chemistry, 1. Model formulation, 
!        J. Geophys. Res., 103/D9, 10, 713-10,726, 1998.
!  (2 ) Yienger, J.J, and H. Levy, Empirical model of global soil-biogenic
!        NOx emissions, J. Geophys. Res., 100, D6, 11,447-11464, June 20, 1995.
!  (3 ) Yan, X., T. Ohara, and H. Akimoto, Statistical modeling of global
!        soil NOx emissions, Global Biogeochem. Cycles, 19, GB3019,
!        doi:10.1029/2004GB002276, 2005.
!  (4 ) Potter, P., Ramankutty, N., Bennett, E., and Donner, S.: 
!        Characterizing the Spatial Patterns of Global Fertilizer Application
!        and Manure Production, Earth Interactions, 14, 1-22, 
!        10.1175/2009EI288.1, 2010.
!  (5 ) Moore, N.E., Improving global bottom-up biogenic soil NOx inventories,
!        Master's Thesis, Dalhousie University, 2007.
!  (6 ) Hudman, R.C., N.E. Moore, A.K. Mebust, R.V. Martin, A.R. Russell, 
!        L.C. Valin, and R.C Cohen, Steps toward a mechanistic model of global
!        soil nitric oxide emissions: implementation and space
!        based-constraints, Atmos. Chem. Phys., 12, 7779-7795,
!        doi:10.5194/acp-12-7779-2012, 2012.
!
! !INTERFACE:
!
MODULE HCOX_SOILNOX_MOD 
!
! !USES:
!
  USE HCO_ERROR_MOD
  USE HCO_DIAGN_MOD
  USE HCOX_State_MOD,     ONLY : Ext_State
  USE HCO_STATE_MOD,      ONLY : HCO_State

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC :: HcoX_SoilNox_Run
  PUBLIC :: HcoX_SoilNox_Init
  PUBLIC :: HcoX_SoilNox_Final
!
! !REVISION HISTORY:
!
!  17 Aug 2009 - R. Yantosca     - Columnized and cleaned up
!  17 Aug 2009 - R. Yantosca     - Added ProTeX headers
!  31 Jan 2011 - R. Hudman       - Added new code12259.perceus-ucb0
!  31 Jan 2011 - R. Hudman       - Updated headers 
!  29 Aug 2012 - J.D. Maasakkers - Implemented Jacob and Bakwin CRF
!  29 Aug 2012 - J.D. Maasakkers - Adapted code to work with new (online
!                                  regridded) landfraction, climate and
!                                  fertilizer data
!  29 Aug 2012 - J.D. Maasakkers - Removed all unused Wang et al. code
!                                  (comments) 
!  04 Nov 2013 - C. Keller       - Moved all soil NOx routines into one
!                                  module. Now a HEMCO extension.
!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !MODULE VARIABLES:
!
  INTEGER                       :: ExtNr          ! Extension number
  INTEGER                       :: IDTNO          ! NO tracer ID
  LOGICAL                       :: LFERTILIZERNOX ! Use fertilizer NOx?
  REAL*8                        :: FERT_SCALE     ! fertilizer scale factor

  ! # of MODIS/Koppen biome types
  INTEGER, PARAMETER            :: NBIOM_HSN = 24 

  ! Dry period length (from restart)
  REAL(hp), ALLOCATABLE, TARGET :: DRYPERIOD_HSN    (:,:  )

  ! Pulse factors (from restart)
  REAL(hp), ALLOCATABLE, TARGET :: PFACTOR_HSN      (:,:  )
  REAL(hp), ALLOCATABLE, TARGET :: GWET_PREV_HSN    (:,:  )

  ! Deposition reservoir (from restart)
  REAL(hp), ALLOCATABLE, TARGET :: DEP_RESERVOIR_HSN(:,:  )

  ! Instantaneous soil NOx and fertilizer
  REAL*8,  ALLOCATABLE          :: INST_SOIL_HSN    (:,:  )
  REAL*8,  ALLOCATABLE          :: INST_FERT_HSN    (:,:  )

  ! NOx in the canopy
  REAL*8,  ALLOCATABLE          :: CANOPYNOX    (:,:,:)

  ! MODIS landtype
  TYPE MODL
     REAL(hp), POINTER          :: VAL(:,:,:)
  ENDTYPE MODL
  TYPE(MODL), POINTER           :: LANDTYPE(:) => NULL()

  ! Soil fertilizer 
  REAL(hp), POINTER             :: SOILFERT     (:,:) => NULL()

  ! Fraction of arid and non-arid land
  REAL(hp), POINTER             :: CLIMARID     (:,:) => NULL()
  REAL(hp), POINTER             :: CLIMNARID    (:,:) => NULL()
!
! !DEFINED PARAMETERS:
!
  ! Canopy wind extinction coefficients
  ! (cf. Yienger & Levy [1995], Sec 5), now a function of the
  ! MODIS/KOPPEN biometype (J.D. Maasakkers)
  REAL*8,  PARAMETER :: SOILEXC(NBIOM_HSN) = (/           & 
                            0.10, 0.50, 0.10, 0.10, 0.10, &
                            0.10, 0.10, 0.10, 0.10, 1.00, &
                            1.00, 1.00, 1.00, 2.00, 4.00, &
                            4.00, 4.00, 4.00, 4.00, 4.00, &
                            4.00, 2.00, 0.10, 2.00         /)

  ! Steinkamp and Lawrence, 2011 A values, wet biome coefficients
  ! for each of the 24 soil biomes [ng N/m2/s].
  REAL*8,  PARAMETER :: A_BIOME(NBIOM_HSN) = (/           &
                            0.00, 0.00, 0.00, 0.00, 0.00, &
                            0.06, 0.09, 0.09, 0.01, 0.84, &
                            0.84, 0.24, 0.42, 0.62, 0.03, &
                            0.36, 0.36, 0.35, 1.66, 0.08, &
                            0.44, 0.57, 0.57, 0.57         /)

  ! "A" coefficients for converting surface temp to soil temp
  ! for each of the 24 soil biomes
  REAL*8,  PARAMETER :: SOILTA(NBIOM_HSN)  = (/           &
                            0.00, 0.92, 0.00, 0.66, 0.66, &
                            0.66, 0.66, 0.66, 0.66, 0.66, &
                            0.66, 0.66, 0.66, 0.66, 0.84, &
                            0.84, 0.84, 0.84, 0.84, 0.84, &
                            0.84, 1.03, 1.03, 1.03         /)

  ! "B" coefficients for converting surface temp to soil temp
  ! for each of the 24 soil biomes
  REAL*8,  PARAMETER :: SOILTB(NBIOM_HSN)  = (/           &
                            0.00, 4.40, 0.00, 8.80, 8.80, &
                            8.80, 8.80, 8.80, 8.80, 8.80, &
                            8.80, 8.80, 8.80, 8.80, 3.60, &
                            3.60, 3.60, 3.60, 3.60, 3.60, &
                            3.60, 2.90, 2.90, 2.90         /)

  ! MODIS/Koppen resistance values
  INTEGER, PARAMETER :: SNIMODIS(NBIOM_HSN) = (/          &
                               1,    2,    3,    4,    5, &
                               6,    7,    8,    9,   10, &
                              11,   12,   13,   14,   15, &
                              16,   17,   18,   19,   20, &
                              21,   22,   23,   24         /)

  INTEGER, PARAMETER :: SNIRI   (NBIOM_HSN) = (/          &
                            9999,  200, 9999, 9999, 9999, &
                            9999,  200,  200,  200,  200, &
                             200,  200,  200,  200,  200, &
                             200,  200,  400,  400,  200, &
                             200,  200, 9999,  200         /)

  INTEGER, PARAMETER :: SNIRLU  (NBIOM_HSN) = (/          &
                            9999, 9000, 9999, 9999, 9999, &
                            9999, 9000, 9000, 9000, 9000, &
                            9000, 9000, 9000, 9000, 9000, &
                            1000, 9000, 9000, 9000, 9000, &
                            1000, 9000, 9999, 9000         /)

  INTEGER, PARAMETER :: SNIRAC  (NBIOM_HSN) = (/          &
                               0,  300,    0,    0,    0, &
                               0,  100,  100,  100,  100, &
                             100,  100,  100,  100, 2000, &
                            2000, 2000, 2000, 2000, 2000, &
                            2000,  200,  100,  200         /)

  INTEGER, PARAMETER :: SNIRGSS (NBIOM_HSN) = (/          &
                               0,    0,  100, 1000,  100, &
                            1000,  350,  350,  350,  350, &
                             350,  350,  350,  350,  500, &
                             200,  500,  500,  500,  500, &
                             200,  150,  400,  150         /)

  INTEGER, PARAMETER :: SNIRGSO (NBIOM_HSN) = (/          &
                            2000, 1000, 3500,  400, 3500, &
                             400,  200,  200,  200,  200, &
                             200,  200,  200,  200,  200, &
                             200,  200,  200,  200,  200, &
                             200,  150,  300,  150         /)

  INTEGER, PARAMETER :: SNIRCLS (NBIOM_HSN) = (/          &
                            9999, 2500, 9999, 9999, 9999, &
                            9999, 2000, 2000, 2000, 2000, &
                            2000, 2000, 2000, 2000, 2000, &
                            9999, 2000, 2000, 2000, 2000, &
                            9999, 2000, 9999, 2000         /)


  INTEGER, PARAMETER :: SNIRCLO (NBIOM_HSN) = (/          &
                            9999, 1000, 1000, 9999, 1000, &
                            9999, 1000, 1000, 1000, 1000, & 
                            1000, 1000, 1000, 1000, 1000, &
                            9999, 1000, 1000, 1000, 1000, &
                            9999, 1000, 9999, 1000         /)

  INTEGER, PARAMETER :: SNIVSMAX(NBIOM_HSN) = (/          &
                              10,  100,  100,   10,  100, &
                              10,  100,  100,  100,  100, &
                             100,  100,  100,  100,  100, &
                             100,  100,  100,  100,  100, &
                             100,  100,  100,  100         /)

  ! Conversion factor from kg NO to ng N
  REAL*8, PARAMETER   :: kgNO_to_ngN = 4.666d11 !(14/30 * 1e12)

CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCOX_SOILNOX_RUN 
!
! !DESCRIPTION: Subroutine HcoX\_SoilNox\_Run is the driver routine to 
! calculate ship NOx emissions for the current time step. Emissions in
! [kg/m2/s] are added to the emissions array of the passed  
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HcoX_SoilNox_Run ( am_I_Root, ExtState, HcoState, RC )
!
! !USES:
!
    USE HCO_FLUXARR_MOD,    ONLY : HCO_EmisAdd
    USE HCO_EMISLIST_MOD,   ONLY : EmisList_GetDataArr
!
! !INPUT PARAMETERS:
!
    LOGICAL,         INTENT(IN   )  :: am_I_Root
    TYPE(Ext_State), POINTER        :: ExtState    ! Module options

!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(HCO_State), POINTER        :: HcoState   ! Output obj
    INTEGER,         INTENT(INOUT)  :: RC 

! !REVISION HISTORY:
!  05 Nov 2013 - C. Keller - Initial Version
!
! !NOTES: 
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    LOGICAL                :: aIR
    INTEGER                :: I, J
    REAL(hp), TARGET       :: FLUX_2D(HcoState%NX,HcoState%NY)
    REAL*8                 :: FERTDIAG, DEP_FERT, SOILFRT
    REAL*8, PARAMETER      :: SEC_PER_YEAR = 3.1536d7
    REAL*4                 :: TSEMIS
    REAL*8                 :: UNITCONV, IJFLUX
    REAL(hp), POINTER      :: TmpArr(:,:) => NULL()
    REAL(hp), POINTER      :: Arr2D (:,:) => NULL()
    LOGICAL, SAVE          :: FIRST = .TRUE.
    CHARACTER(LEN=255)     :: MSG

    !=================================================================
    ! HCOX_SOILNOX_RUN begins here!
    !=================================================================

    ! Enter 
    CALL HCO_ENTER ( 'HcoX_SoilNox_Run (HcoX_SoilNox_Mod.F90)', RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Return if extension disabled 
    IF ( .NOT. ExtState%SoilNOx ) RETURN

    ! Pass arguments from options object
    aIR   = am_I_Root

    ! Conversion factor from ng N to kg NO
    UNITCONV = 1.0d-12 / 14.0d0 * HcoState%Spc(IDTNO)%EmMW_g

    ! Get all arrays needed for SoilNOx calculation
    CALL EmisList_GetDataArr ( aIR, 'SOILNOX_LANDK1',  LANDTYPE(1)%VAL,  RC )
    IF ( RC /= HCO_SUCCESS ) RETURN
    CALL EmisList_GetDataArr ( aIR, 'SOILNOX_LANDK2',  LANDTYPE(2)%VAL,  RC )
    IF ( RC /= HCO_SUCCESS ) RETURN
    CALL EmisList_GetDataArr ( aIR, 'SOILNOX_LANDK3',  LANDTYPE(3)%VAL,  RC )
    IF ( RC /= HCO_SUCCESS ) RETURN
    CALL EmisList_GetDataArr ( aIR, 'SOILNOX_LANDK4',  LANDTYPE(4)%VAL,  RC )
    IF ( RC /= HCO_SUCCESS ) RETURN
    CALL EmisList_GetDataArr ( aIR, 'SOILNOX_LANDK5',  LANDTYPE(5)%VAL,  RC )
    IF ( RC /= HCO_SUCCESS ) RETURN
    CALL EmisList_GetDataArr ( aIR, 'SOILNOX_LANDK6',  LANDTYPE(6)%VAL,  RC )
    IF ( RC /= HCO_SUCCESS ) RETURN
    CALL EmisList_GetDataArr ( aIR, 'SOILNOX_LANDK7',  LANDTYPE(7)%VAL,  RC )
    IF ( RC /= HCO_SUCCESS ) RETURN
    CALL EmisList_GetDataArr ( aIR, 'SOILNOX_LANDK8',  LANDTYPE(8)%VAL,  RC )
    IF ( RC /= HCO_SUCCESS ) RETURN
    CALL EmisList_GetDataArr ( aIR, 'SOILNOX_LANDK9',  LANDTYPE(9)%VAL,  RC )
    IF ( RC /= HCO_SUCCESS ) RETURN
    CALL EmisList_GetDataArr ( aIR, 'SOILNOX_LANDK10', LANDTYPE(10)%VAL, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN
    CALL EmisList_GetDataArr ( aIR, 'SOILNOX_LANDK11', LANDTYPE(11)%VAL, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN
    CALL EmisList_GetDataArr ( aIR, 'SOILNOX_LANDK12', LANDTYPE(12)%VAL, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN
    CALL EmisList_GetDataArr ( aIR, 'SOILNOX_LANDK13', LANDTYPE(13)%VAL, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN
    CALL EmisList_GetDataArr ( aIR, 'SOILNOX_LANDK14', LANDTYPE(14)%VAL, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN
    CALL EmisList_GetDataArr ( aIR, 'SOILNOX_LANDK15', LANDTYPE(15)%VAL, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN
    CALL EmisList_GetDataArr ( aIR, 'SOILNOX_LANDK16', LANDTYPE(16)%VAL, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN
    CALL EmisList_GetDataArr ( aIR, 'SOILNOX_LANDK17', LANDTYPE(17)%VAL, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN
    CALL EmisList_GetDataArr ( aIR, 'SOILNOX_LANDK18', LANDTYPE(18)%VAL, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN
    CALL EmisList_GetDataArr ( aIR, 'SOILNOX_LANDK19', LANDTYPE(19)%VAL, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN
    CALL EmisList_GetDataArr ( aIR, 'SOILNOX_LANDK20', LANDTYPE(20)%VAL, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN
    CALL EmisList_GetDataArr ( aIR, 'SOILNOX_LANDK21', LANDTYPE(21)%VAL, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN
    CALL EmisList_GetDataArr ( aIR, 'SOILNOX_LANDK22', LANDTYPE(22)%VAL, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN
    CALL EmisList_GetDataArr ( aIR, 'SOILNOX_LANDK23', LANDTYPE(23)%VAL, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN
    CALL EmisList_GetDataArr ( aIR, 'SOILNOX_LANDK24', LANDTYPE(24)%VAL, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN
    CALL EmisList_GetDataArr ( aIR, 'SOILNOX_FERT',    SOILFERT,         RC )
    IF ( RC /= HCO_SUCCESS ) RETURN
    CALL EmisList_GetDataArr ( aIR, 'SOILNOX_ARID',    CLIMARID,         RC )
    IF ( RC /= HCO_SUCCESS ) RETURN
    CALL EmisList_GetDataArr ( aIR, 'SOILNOX_NONARID', CLIMNARID,        RC )
    IF ( RC /= HCO_SUCCESS ) RETURN
    
    ! NOTE:
    ! On first call, also initialize variables obtained from the 
    ! restart file specified in the configuration file. Note that
    ! this may lead to errorenous results if the time slices in the
    ! restart file don't match the simulation start date!!
    ! In future, we may want to include these values into the GEOS-Chem
    ! restart file! 
    IF ( FIRST ) THEN

       ! DEP_RESERVOIR. HEMCO converts this to kgNO/m2. Convert here
       ! back to ngN/m2.
       CALL EmisList_GetDataArr ( aIR, 'SOILNOX_DEPRES', TmpArr, RC )
       IF ( RC /= HCO_SUCCESS ) RETURN
       DEP_RESERVOIR_HSN(:,:) = TmpArr(:,:) * kgNO_to_ngN
       TmpArr => NULL()

       ! GWET_PREV [unitless]
       CALL EmisList_GetDataArr ( aIR, 'SOILNOX_GWET', TmpArr, RC )
       IF ( RC /= HCO_SUCCESS ) RETURN
       GWET_PREV_HSN(:,:) = TmpArr(:,:)
       TmpArr => NULL()
          
       ! PFACTOR [unitless]
       CALL EmisList_GetDataArr ( aIR, 'SOILNOX_PFACT', TmpArr, RC )
       IF ( RC /= HCO_SUCCESS ) RETURN
       PFACTOR_HSN(:,:) = TmpArr(:,:)
       TmpArr => NULL()
          
       ! DRYPERIOD [unitless]
       CALL EmisList_GetDataArr ( aIR, 'SOILNOX_DRYPER', TmpArr, RC )
       IF ( RC /= HCO_SUCCESS ) RETURN
       DRYPERIOD_HSN(:,:) = TmpArr(:,:)
       TmpArr => NULL()

       ! Print a warning
       MSG = 'Restart variables obtained from restart file - ' // &
            'I hope the correct time stamp is used!'
       CALL HCO_WARNING( MSG, RC )

       FIRST = .FALSE.
    ENDIF

    ! Now need to call GET_CANOPY_NOX to break ugly dependency between
    ! drydep and soil NOx emissions. (bmy, 6/22/09) 
    ! Now a function of the new MODIS/Koppen biome map (J.D. Maasakkers)
    CALL GET_CANOPY_NOX( aIR, HcoState, ExtState )

    ! Init
    FLUX_2D = 0d0 
    TSEMIS  = HcoState%TS_EMIS

    ! Loop over each land grid-box, removed loop over landpoints
!$OMP PARALLEL DO                                               &
!$OMP DEFAULT( SHARED )                                         &
!$OMP PRIVATE( I, J, DEP_FERT, SOILFRT, FERTDIAG, IJFLUX )
    DO J = 1, HcoState%NY 
    DO I = 1, HcoState%NX

       ! Do not calculate soil NOx emissions if there is no soil in
       ! the gridbox
       IF ( LANDTYPE(1)%VAL(I,J,1) == 1.0_hp ) CYCLE

       ! Get Deposited Fertilizer DEP_FERT [kg NO/m2]
       CALL GET_DEP_N( I, J, ExtState, HcoState, DEP_FERT )
       DEP_FERT = DEP_FERT / sec_per_year
       CALL FLUSH(6)

       ! Get N fertilizer reservoir associated with chemical and
       ! manure fertilizer [kg NO/m2/s]
       IF ( LFERTILIZERNOX ) THEN
          SOILFRT =  SOILFERT( I, J )  
       ELSE
          SOILFRT = 0d0
       ENDIF

       ! Convert to kg NO/m2/s
       !SOILFRT  = SOILFRT / sec_per_year

       ! Put in constraint if dry period gt 1 yr, keep at 1yr to
       ! avoid unrealistic pulse
       IF ( DRYPERIOD_HSN(I,J) > 8760d0 ) DRYPERIOD_HSN(I,J) = 8760d0
 
       ! Return NO emissions from soils [kg NO/m2/s]
       CALL SOIL_NOX_EMISSION( ExtState,                                &
                               TSEMIS, I, J,                            &
                               SOILFRT,                                 &
                               GWET_PREV_HSN(I,J), DRYPERIOD_HSN(I,J),  &
                               PFACTOR_HSN(I,J),   IJFLUX,              &
                               DEP_FERT,           FERTDIAG,            &
                               UNITCONV,                                &
                               CANOPYNOX(I,J,:)                  )
     
       ! Write out
       FLUX_2D(I,J) = IJFLUX

    ENDDO !J 
    ENDDO !I
!$OMP END PARALLEL DO

    !-----------------------------------------------------------------
    ! PASS TO HEMCO STATE AND UPDATE DIAGNOSTICS 
    !-----------------------------------------------------------------

    ! Add flux to emission array
    CALL HCO_EmisAdd( HcoState, FLUX_2D, IDTNO, RC)
    IF ( RC /= HCO_SUCCESS ) RETURN 

    ! Eventually update diagnostics
    IF ( Diagn_AutoFillLevelDefined(2) ) THEN
       Arr2D => FLUX_2D
       CALL Diagn_Update( am_I_Root, HcoState, ExtNr=ExtNr, &
                          Cat=-1, Hier=-1, HcoID=IDTNO,     &
                          AutoFill=1, Array2D=Arr2D, RC=RC   )
       IF ( RC /= HCO_SUCCESS ) RETURN 
       Arr2D => NULL() 
    ENDIF

    ! Cleanup 
    DO I = 1, NBIOM_HSN
       LANDTYPE(I)%VAL => NULL() 
    ENDDO
    SOILFERT  => NULL()
    CLIMARID  => NULL()
    CLIMNARID => NULL()

    ! Leave w/ success
    CALL HCO_LEAVE ( RC ) 

  END SUBROUTINE HcoX_SoilNox_Run
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCOX_SOILNOX_INIT 
!
! !DESCRIPTION: Subroutine HcoX\_SoilNox\_Init initializes the HEMCO
! SOILNOX extension.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HcoX_SoilNox_Init ( am_I_Root, HcoState, ExtName, &
                                 ExtState,    RC                  )
!
! !USES:
!
    USE HCOX_ExtList_Mod, ONLY : GetExtNr, GetExtHcoID, GetExtOpt
!
! !ARGUMENTS:
!
    LOGICAL,          INTENT(IN   )  :: am_I_Root
    TYPE(HCO_State),  POINTER        :: HcoState   ! Output obj
    CHARACTER(LEN=*), INTENT(IN   )  :: ExtName    ! Extension name
    TYPE(Ext_State),  POINTER        :: ExtState     ! Module options
    INTEGER,          INTENT(INOUT)  :: RC 

! !REVISION HISTORY:
!  05 Nov 2013 - C. Keller - Initial Version
!
! !NOTES: 
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    CHARACTER(LEN=255)             :: MSG, LOC
    CHARACTER(LEN=31), ALLOCATABLE :: SpcNames(:)
    INTEGER, ALLOCATABLE           :: HcoIDs(:)
    INTEGER                        :: nSpc, I, J, II, AS

    !=================================================================
    ! HCOX_SOILNOX_INIT begins here!
    !=================================================================

    ! Extension Nr.
    ExtNr = GetExtNr( TRIM(ExtName) )
    IF ( ExtNr <= 0 ) RETURN

    ! Enter 
    CALL HCO_ENTER ( 'HcoX_SoilNox_Init (HcoX_SoilNox_Mod.F90)', RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! ---------------------------------------------------------------------- 
    ! Get species IDs and settings 
    ! ---------------------------------------------------------------------- 

    ! Read settings specified in configuration file
    ! Note: the specified strings have to match those in 
    !       the config. file!
    CALL GetExtOpt ( ExtNr, 'fertilizer NO', &
                     OptValBool=LFERTILIZERNOX, RC=RC )
    IF ( RC /= HCO_SUCCESS ) RETURN
 
    ! Get global scale factor
    FERT_SCALE = HcoX_SoilNOx_GetFertScale()
 
    ! Get HEMCO species IDs
    CALL GetExtHcoID( HcoState, ExtNr, HcoIDs, SpcNames, nSpc, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN
    IF ( nSpc /= 1 ) THEN
       MSG = 'Module soil NOx accepts only one species!'
       CALL HCO_ERROR ( MSG, RC )
       RETURN
    ENDIF
    IDTNO = HcoIDs(1)

    ! Verbose mode
    MSG = 'Use soil NOx emissions (extension module)'
    CALL HCO_MSG( MSG )

    WRITE(MSG,*) '   - NOx species:', TRIM(SpcNames(1)), IDTNO
    CALL HCO_MSG(MSG)
    WRITE(MSG,*) '   - Use fertilizer NOx: ', LFERTILIZERNOX
    CALL HCO_MSG(MSG)
    WRITE(MSG,*) '   - Global scale factor: ', FERT_SCALE 
    CALL HCO_MSG(MSG)
    MSG = '   --> Restart variables are taken from file specified in'
    CALL HCO_MSG(MSG)
    MSG = '       the config. file - I hope this file contains the'
    CALL HCO_MSG(MSG)
    MSG = '       simulation start date!!!' 
    CALL HCO_MSG(MSG,SEP2='-')

    ! ---------------------------------------------------------------------- 
    ! Set module variables
    ! ---------------------------------------------------------------------- 

    ! horizontal dimensions
    I = HcoState%NX
    J = HcoState%NY

    ALLOCATE( DRYPERIOD_HSN    ( I, J        ), STAT=AS )
    IF ( AS /= 0 ) THEN
       CALL HCO_ERROR('DRYPERIOD_HSN',     RC )
       RETURN
    ENDIF

    ALLOCATE( PFACTOR_HSN      ( I, J        ), STAT=AS )
    IF ( AS /= 0 ) THEN
       CALL HCO_ERROR('PFACTOR_HSN',       RC )
       RETURN
    ENDIF

    ALLOCATE( GWET_PREV_HSN    ( I, J        ), STAT=AS )
    IF ( AS /= 0 ) THEN
       CALL HCO_ERROR('GWET_PREV_HSN',     RC )
       RETURN
    ENDIF

    ALLOCATE( INST_SOIL_HSN    ( I, J        ), STAT=AS )
    IF ( AS /= 0 ) THEN
       CALL HCO_ERROR('INST_SOIL_HSN',     RC )
       RETURN
    ENDIF

    ALLOCATE( INST_FERT_HSN    ( I, J        ), STAT=AS )
    IF ( AS /= 0 ) THEN
       CALL HCO_ERROR('INST_FERT_HSN',     RC )
       RETURN
    ENDIF

    ALLOCATE( DEP_RESERVOIR_HSN( I, J        ), STAT=AS )
    IF ( AS /= 0 ) THEN
       CALL HCO_ERROR('DEP_RESERVOIR_HSN', RC )
       RETURN
    ENDIF

    ALLOCATE( CANOPYNOX( I, J, NBIOM_HSN     ), STAT=AS )
    IF ( AS /= 0 ) THEN
       CALL HCO_ERROR('CANOPYNOX',         RC )
       RETURN
    ENDIF

    ! Reserve 24 pointers for land fractions for each Koppen category
    ALLOCATE ( LANDTYPE(NBIOM_HSN), STAT=AS )
    IF ( AS /= 0 ) THEN
       CALL HCO_ERROR('LANDTYPE',           RC )
       RETURN
    ENDIF
    DO II = 1,NBIOM_HSN
       LANDTYPE(NBIOM_HSN)%VAL => NULL()
    ENDDO

    ! Zero arrays
    DRYPERIOD_HSN     = 0.0_hp
    PFACTOR_HSN       = 0.0_hp
    GWET_PREV_HSN     = 0.0_hp
    INST_SOIL_HSN     = 0d0
    INST_FERT_HSN     = 0d0
    CANOPYNOX         = 0d0
    DEP_RESERVOIR_HSN = 0.0_hp

    ! ---------------------------------------------------------------------- 
    ! Set diagnostics 
    ! ---------------------------------------------------------------------- 
    CALL Diagn_Create ( am_I_Root, HcoState,      &    
                        cName      = 'PFACTOR',   &
                        ExtNr      = ExtNr,       &
                        Cat        = -1,          &
                        Hier       = -1,          &
                        HcoID      = IDTNO,       &
                        SpaceDim   = 2,           &
                        OutUnit    = 'unitless',  &
                        WriteFreq  = 'End',       &
                        AutoFill   = 0,           &
                        Trgt2D     = PFACTOR_HSN, &
                        cID        = II,          &
                        RC         = RC            )
    IF ( RC /= HCO_SUCCESS ) RETURN

    CALL Diagn_Create ( am_I_Root, HcoState,        &    
                        cName      = 'DRYPERIOD',   &
                        ExtNr      = ExtNr,         &
                        Cat        = -1,            &
                        Hier       = -1,            &
                        HcoID      = IDTNO,         &
                        SpaceDim   = 2,             &
                        OutUnit    = 'unitless',    &
                        WriteFreq  = 'End',         &
                        AutoFill   = 0,             &
                        Trgt2D     = DRYPERIOD_HSN, &
                        cID        = II,            &
                        RC         = RC              )
    IF ( RC /= HCO_SUCCESS ) RETURN

    CALL Diagn_Create ( am_I_Root, HcoState,        &    
                        cName      = 'GWET_PREV',   &
                        ExtNr      = ExtNr,         &
                        Cat        = -1,            &
                        Hier       = -1,            &
                        HcoID      = IDTNO,         &
                        SpaceDim   = 2,             &
                        OutUnit    = 'unitless',    &
                        WriteFreq  = 'End',         &
                        AutoFill   = 0,             &
                        Trgt2D     = GWET_PREV_HSN, &
                        cID        = II,            &
                        RC         = RC              )
    IF ( RC /= HCO_SUCCESS ) RETURN

    CALL Diagn_Create ( am_I_Root, HcoState,            &    
                        cName      = 'DEP_RESERVOIR',   &
                        ExtNr      = ExtNr,             &
                        Cat        = -1,                &
                        Hier       = -1,                &
                        HcoID      = IDTNO,             &
                        SpaceDim   = 2,                 &
                        OutUnit    = 'ngN/m2',          &
                        WriteFreq  = 'End',             &
                        AutoFill   = 0,                 &
                        Trgt2D     = DEP_RESERVOIR_HSN, &
                        cID        = II,                &
                        RC         = RC                  )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! ---------------------------------------------------------------------- 
    ! Set HEMCO extensions variables 
    ! ---------------------------------------------------------------------- 

    ! Activate required met fields
    ExtState%TSURFK%DoUse    = .TRUE. 
    ExtState%GWETTOP%DoUse   = .TRUE. 
    ExtState%SUNCOSmid%DoUse = .TRUE. 
    ExtState%U10M%DoUse      = .TRUE. 
    ExtState%V10M%DoUse      = .TRUE. 
    ExtState%GC_LAI%DoUse    = .TRUE. 
    ExtState%ALBD%DoUse      = .TRUE. 
    ExtState%RADSWG%DoUse    = .TRUE. 
    ExtState%CLDFRC%DoUse    = .TRUE. 

    ! Activate required deposition parameter
    ExtState%DRY_TOTN%DoUse = .TRUE.
    ExtState%WET_TOTN%DoUse = .TRUE.

    ! Enable module
    ExtState%SoilNOx = .TRUE.

    ! Leave w/ success
    IF ( ALLOCATED(HcoIDs  ) ) DEALLOCATE(HcoIDs  )
    IF ( ALLOCATED(SpcNames) ) DEALLOCATE(SpcNames)
    CALL HCO_LEAVE ( RC ) 

  END SUBROUTINE HcoX_SoilNox_Init
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCOX_SOILNOX_FINAL
!
! !DESCRIPTION: Subroutine HcoX\_SoilNox\_Final finalizes the HEMCO
! SOILNOX extension.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HcoX_SoilNox_Final
!
! !REVISION HISTORY:
!  05 Nov 2013 - C. Keller - Initial Version
!
! !NOTES: 
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER  :: I

    !=================================================================
    ! HCOX_SOILNOX_FINAL begins here!
    !=================================================================

    ! Deallocate arrays
    IF ( ALLOCATED(DRYPERIOD_HSN) ) DEALLOCATE ( DRYPERIOD_HSN )
    IF ( ALLOCATED(PFACTOR_HSN)   ) DEALLOCATE ( PFACTOR_HSN )
    IF ( ALLOCATED(GWET_PREV_HSN) ) DEALLOCATE ( GWET_PREV_HSN )
    IF ( ALLOCATED(INST_SOIL_HSN) ) DEALLOCATE ( INST_SOIL_HSN )
    IF ( ALLOCATED(INST_FERT_HSN) ) DEALLOCATE ( INST_FERT_HSN )
    IF ( ALLOCATED(CANOPYNOX) ) DEALLOCATE ( CANOPYNOX )
    IF ( ALLOCATED(DEP_RESERVOIR_HSN) ) DEALLOCATE ( DEP_RESERVOIR_HSN )

    ! Deallocate LANDTYPE vector 
    IF ( ASSOCIATED(LANDTYPE) ) THEN
       DO I = 1,NBIOM_HSN
          LANDTYPE(NBIOM_HSN)%VAL => NULL()
       ENDDO
       DEALLOCATE ( LANDTYPE )
    ENDIF

    ! Free pointers 
    CLIMARID  => NULL()
    CLIMNARID => NULL()
    SOILFERT  => NULL()

  END SUBROUTINE HcoX_SoilNox_Final
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HcoX_SoilNox_GetFertScale
!
! !DESCRIPTION: Function HcoX_SoilNox_GetFertScale returns the scale factor
! applied to fertilizer NOx emissions. 
!\\
!\\
! !INTERFACE:
!
  FUNCTION HcoX_SoilNox_GetFertScale RESULT ( FERT_SCALE )
!
! !ARGUMENTS:
!
    REAL*8 :: FERT_SCALE
!
! !REVISION HISTORY:
!  11 Dec 2013 - C. Keller - Initial version 
!
! !NOTES: 
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!

    ! Scale factor so that fertilizer emission = 1.8 Tg N/yr
    ! (Stehfest and Bouwman, 2006)
    ! before canopy reduction
    FERT_SCALE = 0.0068 
    ! Value calculated by running the 2x2.5 model
    ! For now, use this value for all resolutions since regular soil NOx
    ! emissions change with resolution as well (J.D. Maasakkers)

  END FUNCTION HcoX_SoilNox_GetFertScale
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
!IROUTINE: soil_nox_emission
!
! !DESCRIPTION: Subroutine SOIL\_NOX\_EMISSION computes the emission of soil and
!  fertilizer NOx for the GEOS-Chem model.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE SOIL_NOX_EMISSION( ExtState,   TS_EMIS,   I, J, &
                                SOILFRT,   &
                                GWET_PREV_HSN, DRYPERIOD_HSN, &
                                PFACTOR_HSN,   SOILNOx,   &
                                DEPN,      FERTDIAG,  &
                                UNITCONV,  R_CANOPY )
!
! !INPUT PARAMETERS:
!
    TYPE(Ext_State), POINTER :: ExtState     ! Module options
    REAL*4,  INTENT(IN)  :: TS_EMIS          ! Emission timestep [s]
    INTEGER, INTENT(IN)  :: I                ! grid box lon index 
    INTEGER, INTENT(IN)  :: J                ! grid box lat index 
    REAL*8,  INTENT(IN)  :: DEPN             ! Dry Dep Fert term [kg/m2/s]
    REAL*8,  INTENT(IN)  :: SOILFRT          ! Fertilizer emissions [kg/m2/s]
    REAL*8,  INTENT(IN)  :: UNITCONV         ! ng N to kg NO 

    !Input parameters for the canopy reduction factor
    REAL*8,  INTENT(IN)  :: R_CANOPY(:)      ! Resist. of canopy to NOx [1/s]
!
! !OUTPUT PARAMETERS:
!
    REAL*8,   INTENT(OUT) :: SOILNOx         ! Soil NOx emissions [molec/cm2/s]
    REAL(hp), INTENT(OUT) :: GWET_PREV_HSN   ! Soil Moisture Prev timestep
    REAL(hp), INTENT(OUT) :: DRYPERIOD_HSN   ! Dry period length in hours
    REAL(hp), INTENT(OUT) :: PFACTOR_HSN     ! Pulsing Factor
    REAL*8,   INTENT(OUT) :: FERTDIAG        ! Fert emissions [molec/cm2/s
!
! !REMARKS:
!  R_CANOPY is computed in routine GET_CANOPY_NOX of "canopy_nox_mod.f". 
!  This was originally in the GEOS-Chem dry deposition code, but was split 
!  off in order to avoid an ugly code dependency between the dry deposition
!  and soil NOx codes.
! 
!  As of v9-02, this module uses the MODIS/Koppen biome types instead
!  of the Olson land type / biome type, making it different from the original 
!  dry deposition code (J.D. Maasakkers)
!
! !REVISION HISTORY: 
!  17 Aug 2009 - R. Yantosca - Columnized and cleaned up
!  17 Aug 2009 - R. Yantosca - Added ProTeX headers
!  31 Jan 2011 - R. Hudman   - New Model added
!  23 Oct 2012 - M. Payer    - Now reference Headers/gigc_errcode_mod.F90
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER :: K
    REAL*8  :: BASE_TERM, CRF_TERM,  PULSE
    REAL*8  :: TC,        TEMP_TERM, WINDSQR
    REAL*8  :: WET_TERM,  A_FERT,    A_BIOM 
    REAL*8  :: LAI,       SUNCOS,    GWET
    REAL*8  :: ARID,      NARID

    !=================================================================
    ! Initialize
    !=================================================================

    ! Initialize
    SOILNOX        = 0d0
    FERTDIAG       = 0d0

    ! Surface temperature [C]
    TC             = ExtState%TSURFK%Arr%Val(I,J) - 273.15d0

    ! Surface wind speed, squared
    WINDSQR        = ExtState%U10M%Arr%Val(I,J)**2 + &
                     ExtState%V10M%Arr%Val(I,J)**2

    ! Leaf area index
    LAI = ExtState%GC_LAI%Arr%Val(I,J)

    ! Cosine of Solar Zenit Angle
    SUNCOS = ExtState%SUNCOSmid%Arr%Val(I,J)

    ! Top soil wetness [unitless]
    GWET = ExtState%GWETTOP%Arr%Val(I,J)

    !=================================================================
    ! Compute soil NOx emissions
    !=================================================================

    ! Cumulative multiplication factor (over baseline emissions) 
    ! that accounts for soil pulsing
    PULSE = PULSING( GWET, TS_EMIS, GWET_PREV_HSN, PFACTOR_HSN, DRYPERIOD_HSN )

    ! ------Loop Over MODIS/Koppen  Landtypes
    DO K = 1, 24

       ! Temperature-dependent term of soil NOx emissions [unitless]
       ! Use GWET instead of climo wet/dry
       TEMP_TERM = SOILTEMP( K , TC, GWET)

       ! Soil moisture scaling of soil NOx emissions
       ARID     = CLIMARID(I,J)
       NARID    = CLIMNARID(I,J)
       WET_TERM = SOILWET( GWET , ARID, NARID )

       ! Fertilizer emission [kg/m2/s] 
       A_FERT = FERTADD( SOILFRT , DEPN)

       ! Scale fertilizer emissions as specified
       ! (scale needed to force fert emiss of 1.8 Tg N/yr w/o canopy uptake)
       A_FERT = A_FERT * FERT_SCALE 

       ! Canopy reduction factor
       CRF_TERM  = SOILCRF( K, LAI, R_CANOPY(K), WINDSQR, SUNCOS )

       ! Base emission. ng N/m2/s --> kg NO/m2/s
       A_BIOM = A_BIOME(K) * UNITCONV

       ! SOILNOX includes fertilizer
       SOILNOX   = (SOILNOX                            &
                 + ( A_BIOM + A_FERT )                 &
                 * ( TEMP_TERM * WET_TERM * PULSE )    &
                 * LANDTYPE(K)%VAL(I,J,1)              &
                 * ( 1.d0 - CRF_TERM  )                 )

       ! FERTDIAG, only used for the fertilizer diagnostic
       FERTDIAG  = (FERTDIAG                           &
                 + ( A_FERT )                          &
                 * ( TEMP_TERM * WET_TERM * PULSE )    &
                 * LANDTYPE(K)%VAL(I,J,1)              &
                 * ( 1.d0 - CRF_TERM  )                 )

    ENDDO

  END SUBROUTINE SOIL_NOX_EMISSION
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_canopy_nox
!
! !DESCRIPTION: Subroutine GET\_CANOPY\_NOX computes the bulk surface 
!  resistance of the canopy to NOx.  This computation was originally done 
!  within legacy routine DEPVEL (in "drydep\_mod.f").  Moving this computation 
!  to GET\_CANOPY\_NOX now allows for a totally clean separation between 
!  dry deposition routines and emissions routines in GEOS-Chem.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE GET_CANOPY_NOX( am_I_Root, HcoState, ExtState )
!
! !USES:
!
    USE DRYDEP_TOOLBOX_MOD, ONLY : BIOFIT
!
! !ARGUMENTS:
!
    LOGICAL,         INTENT(IN   )  :: am_I_Root
    TYPE(HCO_State), POINTER        :: HcoState   ! Output obj
    TYPE(Ext_State), POINTER        :: ExtState    ! Module options
!
! !REMARKS:
!  For backwards compatibility, the bulk surface resistance is stored
!  in common block array CANOPYNOX in "commsoil.h".  Leave it like this
!  for the time being...we'll clean it up when we fix all of the soil
!  NOx routines.
!
! !REVISION HISTORY:
!  22 Jun 2009 - R. Yantosca     - Split off from "drydep_mod.f"
!  14 Jun 2012 - J.D. Maasakkers - Rewritten as a function of the
!                                     MODIS/Koppen biometype
!  09 Nov 2012 - M. Payer        - Replaced all met field arrays with State_Met
!                                   derived type object
!  13 Dec 2012 - R. Yantosca     - Removed ref to obsolete CMN_DEP_mod.F
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
    ! Molecular weight of water [kg]
    REAL*8, PARAMETER :: XMWH2O = 18d-3

    ! Surface pressure??? [Pa]
    REAL*8, PARAMETER :: PRESS  = 1.5d5
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER          :: I,     J,     K,      KK
    REAL*8           :: F0,    HSTAR, XMW              
    REAL*8           :: DTMP1, DTMP2, DTMP3,  DTMP4, GFACT, GFACI
    REAL*8           :: RT,    RAD0,  RIX,    RIXX,  RDC,   RLUXX
    REAL*8           :: RGSX,  RCLX,  TEMPK,  TEMPC
    REAL*8           :: LAI,   SUNCOS, CLDFRC

    ! Arrays
    REAL*8           :: RI  (NBIOM_HSN) 
    REAL*8           :: RLU (NBIOM_HSN)      
    REAL*8           :: RAC (NBIOM_HSN)      
    REAL*8           :: RGSS(NBIOM_HSN)     
    REAL*8           :: RGSO(NBIOM_HSN)     
    REAL*8           :: RCLS(NBIOM_HSN)     
    REAL*8           :: RCLO(NBIOM_HSN)  

    !=================================================================
    ! GET_CANOPY_NOX begins here!
    !=================================================================

    ! Set physical parameters
    HSTAR = 0.01d0              ! Henry's law constant
    F0    = 0.1d0               ! Reactivity factor for biological oxidation 
    XMW   = 46d-3               ! Molecular wt of NO2 (kg)

    ! Loop over surface boxes
    DO J = 1, HcoState%NY
    DO I = 1, HcoState%NX

       ! Surface temperature [K] and [C]
       TEMPK = ExtState%TSURFK%Arr%Val(I,J)
       TEMPC = ExtState%TSURFK%Arr%Val(I,J) - 273.15d0

       ! Compute bulk surface resistance for gases.    
       !                                  
       !  Adjust external surface resistances for temperature; 
       !  from Wesely [1989], expression given in text on p. 1296.        
       RT = 1000.0D0 * EXP( -TEMPC - 4.0d0 )

       !--------------------------------------------------------------
       ! Get surface resistances - loop over biome types K
       !
       ! The land types within each grid square are defined using the 
       ! Olson land-type database.  Each of the Olson land types is 
       ! assigned a corresponding "deposition land type" with 
       ! characteristic values of surface resistance components.  
       ! There are 74 Olson land-types but only 11 deposition 
       ! land-types (i.e., many of the Olson land types share the 
       ! same deposition characteristics).  Surface resistance 
       ! components for the "deposition land types" are from Wesely 
       ! [1989] except for tropical forests [Jacob and Wofsy, 1990] 
       ! and for tundra [Jacob et al., 1992].  All surface resistance 
       ! components are normalized to a leaf area index of unity.
       !--------------------------------------------------------------
       !Loop over all biometypes
       DO K = 1, NBIOM_HSN

          ! Skip if not present
          IF ( LANDTYPE(K)%VAL(I,J,1) == 0.0_hp ) CYCLE

          ! Set second loop variable to K to allow snow/ice correction    
          KK = K

          ! If the surface is snow or ice, then set K=3
          IF ( ExtState%ALBD%Arr%Val(I,J) > 0.4 ) KK = 3

          ! USE new MODIS/KOPPEN Biometypes to read data

          ! Read the internal resistance RI (minimum stomatal resistance 
          ! for water vapor, per unit area of leaf) from the IRI array; 
          ! a '9999' value means no deposition to stomata so we impose a 
          ! very large value for RI.
          RI(K) = DBLE( SNIRI(KK) )
          IF ( RI(K) >= 9999.D0 ) RI(K)= 1.D12

          ! Cuticular resistances IRLU read in from 'drydep.table'
          ! are per unit area of leaf; divide them by the leaf area index 
          ! to get a cuticular resistance for the bulk canopy.  If IRLU is 
          !'9999' it means there are no cuticular surfaces on which to 
          ! deposit so we impose a very large value for RLU.
          IF ( SNIRLU(KK) >= 9999 .OR. &
               ExtState%GC_LAI%Arr%Val(I,J) <= 0d0 ) THEN
             RLU(K)  = 1.D6
          ELSE
             RLU(K)= DBLE(SNIRLU(KK)) / ExtState%GC_LAI%Arr%Val(I,J) + RT
          ENDIF

          ! The following are the remaining resistances for the Wesely
          ! resistance-in-series model for a surface canopy
          ! (see Atmos. Environ. paper, Fig.1).  
          RAC(K)  = MAX( DBLE( SNIRAC(KK)  ),      1d0 )
          RGSS(K) = MAX( DBLE( SNIRGSS(KK) ) + RT, 1d0 )
          RGSO(K) = MAX( DBLE( SNIRGSO(KK) ) + RT, 1d0 ) 
          RCLS(K) =      DBLE( SNIRCLS(KK) ) + RT           
          RCLO(K) =      DBLE( SNIRCLO(KK) ) + RT 

          IF (  RAC(K) >= 9999.D0 ) RAC(K)  = 1d12
          IF ( RGSS(K) >= 9999.D0 ) RGSS(K) = 1d12
          IF ( RGSO(K) >= 9999.D0 ) RGSO(K) = 1d12
          IF ( RCLS(K) >= 9999.D0 ) RCLS(K) = 1d12         
          IF ( RCLO(K) >= 9999.D0 ) RCLO(K) = 1d12

          !-------------------------------------------------------------
          ! Adjust stomatal resistances for insolation and temperature:  
          ! 
          ! Temperature adjustment is from Wesely [1989], equation (3).
          ! 
          ! Light adjustment by the function BIOFIT is described by Wang 
          ! [1996].  It combines:
          !
          ! - Local dependence of stomal resistance on the intensity I 
          !   of light impinging the leaf; this is expressed as a 
          !   multiplicative factor I/(I+b) to the stomatal resistance 
          !   where b = 50 W m-2
          !   (equation (7) of Baldocchi et al. [1987])
          ! - Radiative transfer of direct and diffuse radiation in the 
          !   canopy using equations (12)-(16) from Guenther et al. 
          !   [1995]
          ! - Separate accounting of sunlit and shaded leaves using
          !   equation (12) of Guenther et al. [1995]
          ! - Partitioning of the radiation at the top of the canopy 
          !   into direct and diffuse components using a 
          !   parameterization to results from an atmospheric radiative 
          !   transfer model [Wang, 1996]
          !
          ! The dependent variables of the function BIOFIT are the leaf 
          ! area index (XYLAI), the cosine of zenith angle (SUNCOS) and 
          ! the fractional cloud cover (CFRAC).  The factor GFACI 
          ! integrates the light dependence over the canopy depth; so
          ! be scaled by LAI to yield a bulk canopy value because that's 
          ! already done in the GFACI formulation.
          !-------------------------------------------------------------

          ! Radiation @ sfc [W/m2]
          RAD0 = ExtState%RADSWG%Arr%Val(I,J)
            
          ! Internal resistance
          RIX  = RI(K)

          ! Skip the following block if the resistance RIX is high
          IF ( RIX < 9999d0 ) THEN
             GFACT = 100.0D0

             IF ( TEMPC > 0.D0 .AND. TEMPC < 40.D0) THEN
                GFACT = 400.D0 / TEMPC / ( 40.0D0 - TEMPC )
             ENDIF

             GFACI = 100.D0

             IF ( RAD0 > 0d0 .AND. ExtState%GC_LAI%Arr%Val(I,J) > 0d0 ) THEN

                LAI    = ExtState%GC_LAI%Arr%Val(I,J)
                SUNCOS = ExtState%SUNCOSmid%Arr%Val(I,J)
                CLDFRC = ExtState%CLDFRC%Arr%Val(I,J)

                GFACI= 1d0 /                           & 
                       BIOFIT( ExtState%DRYCOEFF,        & 
                       LAI,                    &
                       SUNCOS,                 &
                       CLDFRC,                 &
                       SIZE(ExtState%DRYCOEFF) )
             ENDIF
            
             RIX = RIX * GFACT * GFACI
          ENDIF

          ! Compute aerodynamic resistance to lower elements in lower 
          ! part of the canopy or structure, assuming level terrain - 
          ! equation (5) of Wesely [1989].                     
          RDC = 100.D0*(1.0D0+1000.0D0/(RAD0 + 10.D0))

          ! Loop over species; species-dependent corrections to resistances
          ! are from equations (6)-(9) of Wesely [1989].
          !
          ! NOTE: here we only consider NO2 (bmy, 6/22/09)
          RIXX   = RIX * DIFFG( TEMPK, PRESS, XMWH2O ) / &
                         DIFFG( TEMPK, PRESS, XMW    )   &
                 + 1.D0 / ( HSTAR/3000.D0 + 100.D0*F0  )

          RLUXX  = 1.D12

          IF ( RLU(K) < 9999.D0 ) THEN
             RLUXX = RLU(K) / ( HSTAR / 1.0D+05 + F0 )
          ENDIF

          ! To prevent virtually zero resistance to species with huge HSTAR, 
          ! such as HNO3, a minimum value of RLUXX needs to be set. 
          ! The rationality of the existence of such a minimum is 
          ! demonstrated by the observed relationship between Vd(NOy-NOx) 
          ! and Ustar in Munger et al.[1996]; Vd(HNO3) never exceeds 2 cm/s 
          ! in observations. The corresponding minimum resistance is 50 s/m.
          ! was introduced by J.Y. Liang on 7/9/95.
          RGSX = 1d0 / ( HSTAR/1d5/RGSS(K) + F0/RGSO(K) )
          RCLX = 1d0 / ( HSTAR/1d5/RCLS(K) + F0/RCLO(K) )

          ! Get the bulk surface resistance of the canopy
          ! from the network of resistances in parallel and in series 
          ! (Fig. 1 of Wesely [1989])
          DTMP1 = 1.D0 / RIXX
          DTMP2 = 1.D0 / RLUXX
          DTMP3 = 1.D0 / ( RAC(K) + RGSX )
          DTMP4 = 1.D0 / ( RDC      + RCLX )

          ! Save the within canopy depvel of NOx, used in calculating 
          ! the canopy reduction factor for soil emissions [1/s]
          CANOPYNOX(I,J,K) = DTMP1 + DTMP2 + DTMP3 + DTMP4

       ENDDO !K
    ENDDO !I
    ENDDO !J
   ! Return to calling program
  END SUBROUTINE GET_CANOPY_NOx
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: diffg
!
! !DESCRIPTION: Function DIFFG calculates the molecular diffusivity [m2/s] in 
!  air for a gas X of molecular weight XM [kg] at temperature TK [K] and 
!  pressure PRESS [Pa].
!\\
!\\
! !INTERFACE:
!
  FUNCTION DIFFG( TK, PRESS, XM ) RESULT( DIFF_G )
!
! !INPUT PARAMETERS:
!
    REAL*8, INTENT(IN) :: TK      ! Temperature [K]
    REAL*8, INTENT(IN) :: PRESS   ! Pressure [hPa]
    REAL*8, INTENT(IN) :: XM      ! Molecular weight of gas [kg]
!
! !RETURN VALUE:
!
    REAL*8             :: DIFF_G  ! Molecular diffusivity [m2/s]
!
! !REMARKS:
!  We specify the molecular weight of air (XMAIR) and the hard-sphere molecular
!  radii of air (RADAIR) and of the diffusing gas (RADX).  The molecular
!  radius of air is given in a Table on p. 479 of Levine [1988].  The Table
!  also gives radii for some other molecules.  Rather than requesting the user
!  to supply a molecular radius we specify here a generic value of 2.E-10 m for
!  all molecules, which is good enough in terms of calculating the diffusivity
!  as long as molecule is not too big.
! 
! !REVISION HISTORY:
!     22 Jun 2009 - R. Yantosca - Copied from "drydep_mod.f"
!EOP
!------------------------------------------------------------------------------
!BOC      
!
! !DEFINED PARAMETERS:
!
    REAL*8, PARAMETER  :: XMAIR  = 28.8d-3 
    REAL*8, PARAMETER  :: RADAIR = 1.2d-10
    REAL*8, PARAMETER  :: PI     = 3.1415926535897932d0
    REAL*8, PARAMETER  :: RADX   = 1.5d-10
    REAL*8, PARAMETER  :: RGAS   = 8.32d0
    REAL*8, PARAMETER  :: AVOGAD = 6.023d23
!
! !LOCAL VARIABLES:
!
    REAL*8             :: AIRDEN, Z, DIAM, FRPATH, SPEED      

    !=================================================================
    ! DIFFG begins here!
    !=================================================================

    ! Air density
    AIRDEN = ( PRESS * AVOGAD ) / ( RGAS * TK )

    ! DIAM is the collision diameter for gas X with air.
    DIAM   = RADX + RADAIR

    ! Calculate the mean free path for gas X in air: 
    ! eq. 8.5 of Seinfeld [1986];
    Z      = XM  / XMAIR
    FRPATH = 1d0 /( PI * SQRT( 1d0 + Z ) * AIRDEN*( DIAM**2 ) )

    ! Calculate average speed of gas X; eq. 15.47 of Levine [1988]
    SPEED  = SQRT( 8d0 * RGAS * TK / ( PI * XM ) )

    ! Calculate diffusion coefficient of gas X in air; 
    ! eq. 8.9 of Seinfeld [1986]
    DIFF_G = ( 3d0 * PI / 32d0 ) * ( 1d0 + Z ) * FRPATH * SPEED

    ! Return to calling program
  END FUNCTION DIFFG
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: get_dep_N
!
! !DESCRIPTION: Subroutine GET\_DEP\_N sums dry and wet deposition since prev.
!               timestep and calculates contribution to fertilizer N source.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE GET_DEP_N ( I, J, ExtState, HcoState, DEP_FERT )
!
! !INPUT PARAMETERS: 
!
    INTEGER,  INTENT(IN)         :: I
    INTEGER,  INTENT(IN)         :: J
    TYPE(Ext_State), POINTER     :: ExtState    ! Module options
    TYPE(HCO_State), POINTER     :: HcoState   ! Output obj
!
! !INPUT/OUTPUT PARAMETERS: 
!
    ! Dep emitted as Fert [ng N/m2]
    REAL*8 ,  INTENT(INOUT) :: DEP_FERT  
!
! !REVISION HISTORY:
!  23 Oct 2012 - M. Payer    - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
    REAL*8,  PARAMETER :: TAU_MONTHS   = 6. ! Decay rate of dep. N [months]
    REAL*8,  PARAMETER :: SECPERDAY    = 86400.d0
    REAL*8,  PARAMETER :: DAYSPERMONTH = 30.
!
! !LOCAL VARIABLES:
!
    REAL*8             :: DRYN  ! Dry dep. N since prev timestep
                                  ! Units ng N/m2/s     
    REAL*8             :: WETN  ! Wet dep. N since prev timestep 
    REAL*8             :: DEPN  ! dep. N since prev timestep 

    REAL*8             :: C1
    REAL*8             :: C2 
    REAL*8             :: TAU_SEC
    REAL*8             :: TS_SEC  

    !Total all N species & convert molec/cm2/s --> kg NO/m2/s
    DRYN = SOURCE_DRYN( I, J, ExtState, HcoState )

    !Total all N species & convert kg/s --> kg NO/m2/s
    WETN = SOURCE_WETN( I, J, ExtState, HcoState )

    ! Sum wet and dry deposition [kg NO/m2/s]
    DEPN = DRYN + WETN 

    ! Emission Timestep in seconds
    TS_SEC = HcoState%TS_EMIS

    !Do mass balance (see Intro to Atm Chem Chap. 3)
    !m(t) = m(0) * exp(-t/tau) + Source * tau * (1 - exp(-t/tau))

    !convert months -->  seconds (assume 30 days months)
    TAU_SEC = TAU_MONTHS * DAYSPERMONTH * SECPERDAY

    C1 = EXP( - TS_SEC / TAU_SEC)
    C2 = 1.d0 - C1        

    ! ngN/m2
    DEP_RESERVOIR_HSN(I,J) = ( DEP_RESERVOIR_HSN (I,J) * C1 ) &
                           + DEPN * TAU_SEC * C2

    ! 40% runoff. Convert ngN/m2 to kgNO/m2
    DEP_FERT = DEP_RESERVOIR_HSN(I,J) * 0.6d0 / kgNO_to_ngN 

  END SUBROUTINE  GET_DEP_N
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: source_dryN
!
! !DESCRIPTION: Subroutine SOURCE\_DRYN gets dry deposited Nitrogen since
!               last emission time step, converts to ng N/m2/s.
!\\
!\\
! !INTERFACE:
!
  FUNCTION SOURCE_DRYN( I, J, ExtState, HcoState ) RESULT( DRYN )
!
! !INPUT PARAMETERS: 
!
    INTEGER, INTENT(IN)             :: I           
    INTEGER, INTENT(IN)             :: J           
    TYPE(Ext_State), POINTER        :: ExtState    ! Module options
    TYPE(HCO_State), POINTER        :: HcoState   ! Output obj
!
! !RETURN VALUE:
!
    REAL*8               :: DRYN         !Dry dep. N since prev timestep
!
! !REVISION HISTORY:
!  23 Oct 2012 - M. Payer    - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    REAL*8,  PARAMETER   :: CM2_PER_M2  = 1.d4
    REAL*8               :: NTS
 
    ! Divide through by number of chemistry timesteps
    ! because DRY_TOTN is summed over chemistry timesteps
    ! need to get average

    !Molecules/cm2/s --> kg NO/m2/s 
    NTS  = HcoState%TS_EMIS / HcoState%TS_CHEM
    DRYN = ExtState%DRY_TOTN%Arr%Val(I,J) * CM2_PER_M2 / NTS / &
           HcoState%Phys%Avgdr * HcoState%Spc(IDTNO)%MW_g / 1000.0d0

  END FUNCTION SOURCE_DRYN

!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: source_wetN
!
! !DESCRIPTION: Subroutine SOURCE\_WETN gets wet deposited Nitrogen since
!               last emission time step, converts to ng N/m2/s.
!\\
!\\
! !INTERFACE:
!
  FUNCTION SOURCE_WETN( I, J, ExtState, HcoState ) RESULT(WETN )
!
! !INPUT PARAMETERS: 
!
    INTEGER, INTENT(IN)             :: I           
    INTEGER, INTENT(IN)             :: J           
    TYPE(Ext_State), POINTER        :: ExtState    ! Module options
    TYPE(HCO_State), POINTER        :: HcoState   ! Output obj
!
! !RETURN VALUE:
!
    REAL*8               :: WETN         !Dry dep. N since prev timestep
!
! !REVISION HISTORY:
!  23 Oct 2012 - M. Payer    - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    REAL*8               :: NTS,    AREA_M2     
 
    ! Divide through by number of transport timesteps
    ! because WET_TOTN is summed over transport timesteps
    ! need to get average

    NTS     = HcoState%TS_EMIS / HcoState%TS_DYN 
    AREA_M2 = HcoState%Grid%AREA_M2(I,J) 

    ! Total N wet dep
    WETN = ExtState%WET_TOTN%Arr%Val(I,J) / AREA_M2 / NTS

  END FUNCTION SOURCE_WETN
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: soiltemp
!
! !DESCRIPTION: Function SOILTEMP computes the temperature-dependent term
!  of the soil NOx emissions in ng N/m2/s and converts to molec/cm2/s
!\\
!\\
! !INTERFACE:
!
  FUNCTION SOILTEMP( NN, TC, GWET ) RESULT( SOIL_TEMP )
!
! !INPUT PARAMETERS: 
!
    INTEGER, INTENT(IN) :: NN            ! Soil biome type 
    REAL*8,  INTENT(IN) :: TC            ! Surface air temperature [C]
    REAL*8,  INTENT(IN) :: GWET          ! Top soil moisture
!
! !RETURN VALUE:
!
    REAL*8              :: SOIL_TEMP     ! Temperature-dependent term of
                                         ! soil NOx emissions [unitless]
!
! !REMARKS:
!    Based on Ormeci et al., [1999] and Otter et al., [1999]     
!    there exists and entirely exponential relationship between 
!    temperature and soil NOx emissions at constant soil moisture
!    Therefore we use the following relationship based 
!    on Yienger and Levy et al., [1995] for temperatures 0-30C:        
!                     
!                           
!         f(T) =  exp( 0.103+/-0.04 * T ) 
!           in ng N/m2/s    
!                      
!                                                                             
!     where T is the temperature in degrees Celsius....Below 
!     0 C, we assume emissions are zero because they are insignificant 
!     for the purposes of this global source. ...
! 
!  References:
!  ============================================================================
!  (1 ) Ormeci, B., S. L. Sanin, and J. J. Pierce, Laboratory study of 
!        NO flux from agricultural soil: Effects of soil moisture, pH, 
!        and temperature, J. Geophys. Res., 104 ,16211629, 1999.
!  (2 ) Otter, L. B., W. X. Yang, M. C. Scholes, and F. X. Meixner, 
!        Nitric oxide emissions from a southern African savanna, J. 
!        Geophys. Res., 105 , 20,69720,706, 1999.
!  (3 ) Yienger, J.J, and H. Levy, Empirical model of global soil-biogenic
!        NOx emissions, J. Geophys. Res., 100, D6, 11,447-11464, June 20, 1995.
!
! !REVISION HISTORY: 
!  17 Aug 2009 - R. Yantosca - Initial Version
!  17 Aug 2009 - R. Yantosca - Added ProTeX headers
!  31 Jan 2011 - R. Hudman   - Added new soil T dependance 
!  31 Jan 2011 - R. Hudman   - Updated headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES
!     
    REAL*8  :: TMMP

    !==============================================================
    ! 1) Convert from Surface Temp  --> Soil Temp 
    !==============================================================

    ! Save surface air temp in shadow variable TMMP
    TMMP   = TC

    ! DRY
    IF ( GWET < 0.3d0 ) THEN 
     
       ! Convert surface air temperature to model temperature
       ! by adding 5 degrees C to model temperature
       TMMP = TMMP + 5d0

    ! WET
    ELSE

       TMMP = SOILTA(NN) * TMMP + SOILTB(NN)

    ENDIF

    !==============================================================
    ! 2) Compute Temperature Dependence
    !==============================================================

    ! Compute the soil temperature dependence term according 
    ! to equations 9b, 9a of Yienger & Levy [1995].
    ! We now assume that soil response is exponential 0-30C
    ! based on latest observations, caps at 30C

    IF ( TMMP <= 0d0 ) THEN

       ! No soil emissions if temp below freezing
       SOIL_TEMP = 0d0

    ELSE 

       ! Caps temperature response at 30C
       IF ( TMMP >= 30.d0 ) TMMP = 30.d0 
      
       SOIL_TEMP =  EXP( 0.103 * TMMP )

    ENDIF
 
  END FUNCTION SOILTEMP
!EOC
!----------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: soilwet
!
! !DESCRIPTION: Function SOILWET returns the soil moisture scaling 
!  of soil NOx emissions (values from 0-1). 
!\\
!\\
! !INTERFACE:
!
  FUNCTION SOILWET( GWET , ARID, NONARID ) RESULT( WETSCALE )
!
! !INPUT PARAMETERS: 
!
    ! Top soil wetness [unitless]
    REAL*8, INTENT(IN) :: GWET 

    ! Fraction of arid & non-arid soil in the gridbox          
    REAL*8, INTENT(IN) :: ARID       
    REAL*8, INTENT(IN) :: NONARID       
!
! !RETURN_VALUE:
! 
    ! A scaling term between 0-1 based on soil moisture
    REAL*8             :: WETSCALE      
!
! !REMARKS:
!  Soil moisture and temperature and now decoupled, the temperature
!  term is scaled with a value from 0-1 based on water filled pore space
!  WFPS in top-soil. 
! 
!  From N.E. Moore thesis:
!  The response of SNOx is not monotonic to WFPS. SNOx are low for the 
!  extreme values of WFPS (0 and 1). For low values, emissions are 
!  substrate-limited. For high values, emissions are trapped and cannot 
!  diffuse to the surface [Yan et al., 2005]. SNOx dependence on soil 
!  moisture is best described as a Poisson function [Parsons et al., 1996; 
!  Otter et al., 1999; Pierce and Aneja, 2000; Kirkman et al., 2001; 
!  van Dijk and Meixner, 2001; van Dijk et al., 2002]:
! 
!     scaling = a*x*exp(-b*x^2) 
! 
!  where the values of a and b are chosen such that the maximum value 
!  (unity) occurs for WFPS=0.3, which laboratory and field measurements have 
!  found to be the optimal value for emissions in most soils. The typical 
!  range of values are 0.2 (arid) up to 0.45 (floodplain) 
!  [Yang and Meixner, 1997; Ormeci et al., 1999].  
! 
!  Rice paddies no longer have to be scaled as in the Yienger & Levy model.
! 
!  References:
!  ============================================================================
!  (1 ) Galbally, I. E., and R. Roy, Loss of fixed nitrogen from soils 
!        by nitric oxide exhalation, Nature, 275 , 734735, 1978.
!  (2 ) Kirkman, G. A., W. X. Yang, and F. X. Meixner, Biogenic nitric 
!        oxide emissions upscaling: An approach for Zimbabwe, Global 
!        Biogeochemical Cycles, 15 ,1005 1020, 2001.
!  (3 ) Ormeci, B., S. L. Sanin, and J. J. Pierce, Laboratory study of NO 
!        flux from agricultural soil: Effects of soil moisture, pH, and 
!        temperature, J. Geophys. Res., 104 , 16211629, 1999.
!  (4 ) Otter, L. B., W. X. Yang, M. C. Scholes, and F. X. Meixner, 
!        Nitric oxide emissions from a southern African savanna, J. 
!        Geophys. Res., 105 , 20,69720,706, 1999.
!  (5 ) Parsons, D. A., M. C. Scholes, R. J. Scholes, and J. S. Levine, 
!        Biogenic NO emissions from savanna soils as a function of fire 
!        regime, soil type, soil nitrogen, and water status, J. Geophys. 
!        Res., 101 , 23,68323,688, 1996.
!  (6 ) Pierce, J. J., and V. P. Aneja, Nitric oxide emissions from 
!        engineered soil systems, Journal of Environmental Engineering, 
!        pp. 225232, 2000.
!  (7 ) van Dijk, S. M., and J. H. Duyzer, Nitric oxide emissions from 
!        forest soils, J. Geophys. Res., 104 , 15,95515,961, 1999.
!  (8 ) van Dijk, S. M., and F. X. Meixner, Production and consumption of 
!        NO in forest and pasture soils from the Amazon basin, Water, Air,
!        and Soil Pollution: Focus 1 , pp. 119130, 2001.
!  (9 ) Yang, W. X., and F. X. Meixner, Gaseous Nitrogen Emissions from 
!        Grasslands, CAB Int., Wallingford, UK, 1997, 67-71.
!
! !REVISION HISTORY: 
!  17 Aug 2009 - R. Yantosca - Columnized and cleaned up
!  17 Aug 2009 - R. Yantosca - Added ProTeX headers
!  31 Jan 2011 - R. Hudman   - Rewrote scaling scheme
!  31 Jan 2011 - R.Hudman    - Updated ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
    !Scale by soil moisture           
    IF ( ARID .GE. NONARID .AND. ARID .NE. 0) THEN 
       !Arid, max Poisson = 0.2
       WETSCALE = 8.24 * GWET * EXP(-12.5*GWET*GWET)
    ELSE
       !Non-arid, max Poisson = 0.3
       WETSCALE = 5.5 * GWET * EXP( -5.55 * GWET * GWET) 
    ENDIF

  END FUNCTION SOILWET
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: soilcrf
!
! !DESCRIPTION: Computes the canopy reduction factor for the soil NOx
!  emissions according to Jacob \% Bakwin [1991] (and as used in Wang 
!  et al [1998]).
!\\
!\\
! !INTERFACE:
!	
  FUNCTION SOILCRF( K, LAI, CPYNOX, WINDSQR, SUNCOS ) RESULT( SOIL_CRF )
!
! !INPUT PARAMETERS: 
!
    INTEGER, INTENT(IN) :: K          ! Soil biome type
    REAL*8,  INTENT(IN) :: LAI        ! Leaf area index [cm2/cm2]
    REAL*8,  INTENT(IN) :: CPYNOX     ! Bulk sfc resistance to NOx [1/s]
    REAL*8,  INTENT(IN) :: WINDSQR    ! Square of sfc wind speed [m2/s2]
    REAL*8,  INTENT(IN) :: SUNCOS     ! Cosine of solar zenith angle
!
! !RETURN_VALUE:
! 
    REAL*8              :: SOIL_CRF   ! Canopy reduction factor (see below)
!
! !REMARKS:
!  Also note, CANOPYNOX (the bulk surface resistance to NOx) is computed
!  in routine GET_CANOPY_NOx (in "canopy_nox_mod.f") and is passed here
!  as an argument.
!
! !REVISION HISTORY: 
!  17 Aug 2009 - R. Yantosca - Initial Version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
    ! Ventilation velocity for NOx, day & night values [m/s]
    REAL*8,  PARAMETER :: VFDAY   = 1.0d-2
    REAL*8,  PARAMETER :: VFNIGHT = 0.2d-2 
!
! !LOCAL VARIABLES:
!
    REAL*8 :: VFNEW

    ! Pick proper ventilation velocity for day or night
    IF ( SUNCOS > 0d0 ) THEN
       VFNEW = VFDAY              
    ELSE 
       VFNEW = VFNIGHT            
    ENDIF

    ! If the leaf area index and the bulk surface resistance
    ! of the canopy to NOx deposition are both nonzero ...
    IF ( LAI > 0d0 .and. CPYNOX > 0d0 ) THEN

       ! Adjust the ventilation velocity.  
       ! NOTE: SOILEXC(21) is the canopy wind extinction 
       ! coefficient for the tropical rainforest biome.
       VFNEW    = (VFNEW * SQRT( WINDSQR/9d0 * 7d0/LAI     ) * &
                  ( SOILEXC(21)  / SOILEXC(K) ))

       ! Soil canopy reduction factor
       SOIL_CRF = CPYNOX / ( CPYNOX + VFNEW )

    ELSE
     
       ! Otherwise set the soil canopy reduction factor to zero
       SOIL_CRF = 0d0

    ENDIF

  END FUNCTION SOILCRF
!EOC
!-----------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: fertadd
!
! !DESCRIPTION: Function FERTADD computes fertilizer emissions
!\\
!\\
! !INTERFACE:
!
  FUNCTION FERTADD( SOILFRT, DEPN) RESULT( FERT_ADD )
!
! !INPUT PARAMETERS: 
!
    REAL*8, INTENT(IN) :: DEPN      ! N emissions from deposition
    REAL*8, INTENT(IN) :: SOILFRT   ! N emissions from fertilizers
                                    !  read in from disk and passed
                                    !  here as an argument [ng N/m2/s]
!
! !RETURN_VALUE:
! 
    REAL*8            :: FERT_ADD   ! Total Fert emissions
!
! !REMARKS:
!  We use a new spatially explicit data set of chemical and manure fert
!  (native resolution 0.5\B0x0.5\B0) from Potter et al., [2010] 
!  distributed using MODIS EVI seasonality as described in 
!  N.E. Moore thesis, and Hudman et al., in prep.
!    
!  In previous model, fertilizer emissions were emitted instantaneously as 
!  2.5% of applied fertilizer, independent of soil moisture/soil temperature, 
!  so that they were constant over the growing season. 
! 
!  Similar to the YL  parameterization, we now treat fertilizer emissions 
!  as part of the Aw. If we treat the wet biome coefficient as a measure of 
!  available N multiplied by a mean emission rate, we can treat fertilizer 
!  N in the same manner. 
! 
!  AW = SOILAW(BinewsoilAWS_08112011_emissonlyome) + N available in soil 
!       x mean emission rate
! 
!  Instead of choosing an emission rate for each box equivalent to 2.5% 
!  of applied N yearly as done in the YL scheme, we chose the mean emission 
!  rate so that the total global above canopy SNOx due to fertilizer matches 
!  observed estimates of fertilizer emissions of 1.8 Tg N yr-1 from Stehfest 
!  and Bouman [2006].  This treatment allows for interannual and daily 
!  variability in the strength  of response to temperature and precipitation. 
!  Note: this scaling must be set for each resolution. 
! 
!  References:
!  ============================================================================
!  (1 ) Potter, P., Ramankutty, N., Bennett, E.,  and Donner, S.: 
!        Characterizing the Spatial Patterns of Global Fertilizer 
!        Application and Manure Production, Earth Interactions, 
!        in press, 2010.
!  (2 ) Stehfest, E. and L. Bouwman, N2O and NO emission from 
!        agricultural fields and soils under natural vegetation: 
!        summarizing available measurement data and modeling
!        of global annual emissions, Nutrient Cycling in Agroecosystems 
!        (2006), 74:207-228 DOI 10.1007/s10705-006-9000-7.
!
! !REVISION HISTORY: 
!  17 Aug 2009 - R. Yantosca - Columnized and cleaned up
!  17 Aug 2009 - R. Yantosca - Added ProTeX headers
!  31 Jan 2011 - R. Hudman   - Rewrote pulsing scheme
!  31 Jan 2011 - R. Hudman   - Updated ProTex headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
    ! Initialize
    FERT_ADD = 0d0

    ! Soil fert and dep [ kg/m2/s ], a measure of N avail. in soil
    FERT_ADD = SOILFRT  + DEPN

    ! Convert [ng N/m2] --> [kg/m2/s]
    ! (scale needed to force fert emiss of 1.8 Tg N/yr w/o canopy uptake)
!    FERT_ADD = FERT_ADD * FERT_SCALE   

  END FUNCTION FERTADD
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: pulsing
!
! !DESCRIPTION: Function PULSING calculates the increase (or "pulse") of 
!  soil NOx emission that happens after preciptiation falls on dry soil.  
!                                                                             .
!  According to  Yan et al., [2005] , this pulsing process is thought to  
!  be due to a release of inorganic nitrogen trapped on top of the dry soil 
!  and a subsequent reactivation of water-stressed bacteria, which then 
!  metabolize the excess nitrogen. This can happen in seasonally dry
!  grasslands and savannahs or over freshly fertilized fields.
!\\
!\\
! !INTERFACE:
!
  FUNCTION PULSING( GWET,          TS_EMIS,                     & 
                    GWET_PREV_HSN, PFACTOR_HSN, DRYPERIOD_HSN ) &
                    RESULT( THE_PULSING )
!
! !INPUT PARAMETERS: 
!
    REAL*8, INTENT(IN)    :: GWET        ! Soil Moisture 
    REAL*4, INTENT(IN)    :: TS_EMIS     ! Emissions timestep [s]

! !INPUT/OUTPUT PARAMETERS:
!
    REAL(hp), INTENT(INOUT) :: GWET_PREV_HSN   ! Soil Moisture Prev timestep
    REAL(hp), INTENT(INOUT) :: PFACTOR_HSN     ! Pulsing Factor
    REAL(hp), INTENT(INOUT) :: DRYPERIOD_HSN   ! Dry period length in hours
!  
! !RETURN VALUE:
!
    REAL*8                :: THE_PULSING ! Factor to multiply baseline 
                                         ! emissions by to account for
                                         ! soil pulsing of all types
!
! !REMARKS:
!  Soil NOx emissions consist of baseline emissions plus discrete "pulsing"
!  episodes.  We follow thw Yan et al., [2005] algorithm, where the pulse
!  (relative to the flux prewetting) is determined by the antecedent dry 
!  period, with a simple logarithmic relationship,
! 
!     PFACTOR_HSN = 13.01 ln ( DRYPERIOD_HSN ) -  53.6
! 
!  where PFACTOR_HSN is the magnitude of peak flux relative to prewetting flux, 
!  and DRYPERIOD_HSN  is the length of the antecedent dry period in hours.
! 
!  The pulse decays with 
! 
!     PFACTOR_HSN = PFACTOR_HSN * EXP( -0.068d0 * DTSRCE )       
! 
!  References:
!  ============================================================================
!  (1 ) Yan, X., T. Ohara, and H. Akimoto (2005), Statistical modeling of 
!        global soil NOx emissions, Global Biogeochem. Cycles, 19, GB3019, 
!        doi:10.1029/2004GB002276.Section 2.3.3
!
! !REVISION HISTORY: 
!  17 Aug 2009 - R. Yantosca - Columnized and cleaned up
!  17 Aug 2009 - R. Yantosca - Added ProTeX headers
!  31 Jan 2011 - R. Hudman   - Rewrote pulsing scheme
!  31 Jan 2011 - R. Hudman   - Updated ProTex header
!  29 May 2013 - R. Yantosca - Bug fix: prevent log(0) from happening
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    REAL*8  :: DTSRCE, GDIFF 

    !=================================================================
    ! PULSING begins here!
    !=================================================================

    ! Emission timestep [s --> hours]
    DTSRCE = TS_EMIS / 3600d0

    ! If soil moisture less than 0.3 and no pulse is taking place
    IF ( GWET < 0.3D0 .and. PFACTOR_HSN == 1.D0) THEN

       ! Get change in soil moisture since previous timestep
       GDIFF = ( GWET - GWET_PREV_HSN )

       ! If change in soil moisture is > 0.01 (rains)
       IF ( GDIFF > 0.01 ) THEN

          ! Initialize new pulse factor (dry period hours)
          IF ( DRYPERIOD_HSN > 0d0 ) THEN
             PFACTOR_HSN = 13.01d0 * LOG( DRYPERIOD_HSN ) - 53.6d0
          ELSE
             PFACTOR_HSN = -53.6d0
          ENDIF

          ! Reinitialize dry period
          DRYPERIOD_HSN = 0d0

       ! If no rain (i.e.,  change in soil moisture is < 0.01)
       ELSE

          ! Add one timestep to dry period
          DRYPERIOD_HSN = DRYPERIOD_HSN + DTSRCE

       ENDIF

    ! If box is already pulsing , then decay pulse one timestep
    ELSEIF ( PFACTOR_HSN /= 1.d0) THEN

       ! Decay pulse
       PFACTOR_HSN   = PFACTOR_HSN * EXP( -0.068d0 * DTSRCE )

       ! Update dry period
       IF ( GWET < 0.3D0 ) DRYPERIOD_HSN = DRYPERIOD_HSN + DTSRCE

       ! If end of pulse
       IF ( PFACTOR_HSN < 1.d0 ) PFACTOR_HSN = 1.d0
      
    ENDIF

    ! Update soil moisture holder for previous timestep
    GWET_PREV_HSN = GWET

    ! Return the pulsing factor
    THE_PULSING = PFACTOR_HSN

  END FUNCTION PULSING
!EOC
END MODULE HCOX_SOILNOX_MOD
!EOM
