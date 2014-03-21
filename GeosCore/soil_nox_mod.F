!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: soil_nox_mod
!
! !DESCRIPTION: Module containing GEOS-Chem soil NOx emissions routines.
!\\
!\\
! !INTERFACE: 
!
      MODULE SOIL_NOx_MOD
! 
! !USES:
!
      USE CMN_SIZE_MOD

      IMPLICIT NONE
      PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
      PUBLIC  :: SOIL_NOX_EMISSION
!
! !PRIVATE MEMBER FUNCTIONS:
!
      PRIVATE :: SOILTEMP
      PRIVATE :: SOILWET
      PRIVATE :: SOILCRF
      PRIVATE :: FERTADD
      PRIVATE :: PULSING  
!
! !AUTHOR:
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
! !REVISION HISTORY:
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
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
      ! Conversion factor from [ng N/m2/s] to [molec/cm2/s]
      REAL*8,  PARAMETER :: UNITCONV = 4.3d9

      ! New soil biomes based on Steinkamp et al., 2011
      INTEGER, PARAMETER :: NSOIL    = 24

      ! Canopy wind extinction coefficients
      ! (cf. Yienger & Levy [1995], Sec 5), now a function of the
      ! MODIS/KOPPEN biometype (J.D. Maasakkers)
      REAL*8,  PARAMETER :: SOILEXC(NSOIL) = (/ 
     &                      0.10, 0.50, 0.10, 0.10, 0.10,
     &                      0.10, 0.10, 0.10, 0.10, 1.00,
     &                      1.00, 1.00, 1.00, 2.00, 4.00,
     &                      4.00, 4.00, 4.00, 4.00, 4.00,
     &                      4.00, 2.00, 0.10, 2.00        /)

      ! Steinkamp and Lawrence, 2011 A values, wet biome coefficients
      ! for each of the 24 soil biomes
      REAL*8,  PARAMETER :: A_BIOME(NSOIL) = (/ 
     &                      0.00, 0.00, 0.00, 0.00, 0.00,
     &                      0.06, 0.09, 0.09, 0.01, 0.84,
     &                      0.84, 0.24, 0.42, 0.62, 0.03,
     &                      0.36, 0.36, 0.35, 1.66, 0.08,
     &                      0.44, 0.57, 0.57, 0.57        /)

      ! "A" coefficients for converting surface temp to soil temp
      ! for each of the 24 soil biomes
      REAL*8,  PARAMETER :: SOILTA(NSOIL)  = (/
     &                      0.00, 0.92, 0.00, 0.66, 0.66,
     &                      0.66, 0.66, 0.66, 0.66, 0.66,
     &                      0.66, 0.66, 0.66, 0.66, 0.84,
     &                      0.84, 0.84, 0.84, 0.84, 0.84,
     &                      0.84, 1.03, 1.03, 1.03        /)

      ! "B" coefficients for converting surface temp to soil temp
      ! for each of the 24 soil biomes
      REAL*8,  PARAMETER :: SOILTB(NSOIL)  = (/         
     &                      0.00, 4.40, 0.00, 8.80, 8.80,
     &                      8.80, 8.80, 8.80, 8.80, 8.80,
     &                      8.80, 8.80, 8.80, 8.80, 3.60,
     &                      3.60, 3.60, 3.60, 3.60, 3.60,
     &                      3.60, 2.90, 2.90, 2.90        /)

      ! Scale factor so that fertilizer emission = 1.8 Tg N/yr
      ! (Stehfest and Bouwman, 2006)
      ! before canopy reduction
      REAL*8, PARAMETER :: FERT_SCALE = 0.0068 
      ! Value calculated by running the 2x2.5 model
      ! For now, use this value for all resolutions since regular soil NOx
      ! emissions change with resolution as well (J.D. Maasakkers)

      CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: soil_nox_emission
!
! !DESCRIPTION: Subroutine SOIL\_NOX\_EMISSION computes the emission of soil and
!  fertilizer NOx for the GEOS-Chem model.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE SOIL_NOx_EMISSION( TS_EMIS,   TK,
     &                              GWET,      SOILFRT,
     &                              GWET_PREV, DRYPERIOD,
     &                              PFACTOR,   SOILNOx, 
     &                              DEPN,      FERTDIAG, 
     &                              CLIM,      LANDFRAC,
     &                              RC,        SUNCOS,
     &                              U10M,      V10M,
     &				    R_CANOPY,  LAI)
! 
! !USES:
!
      USE GIGC_ERRCODE_MOD, ONLY : GIGC_SUCCESS
!
! !INPUT PARAMETERS:
!
      REAL*8,  INTENT(IN)  :: TS_EMIS     ! Emission timestep [min]
      REAL*8,  INTENT(IN)  :: TK          ! surface temperature [K]
      REAL*8,  INTENT(IN)  :: GWET        ! Top soil wetness [unitless]
      REAL*8,  INTENT(IN)  :: DEPN        ! Dry Dep Fert term [ng N/m2/s]
      REAL*8,  INTENT(IN)  :: SOILFRT     ! Fertilizer emissions [ng N/m2/s]
      REAL*4,  INTENT(IN)  :: CLIM(:)     ! CLIM(1), arid fraction
                                          ! CLIM(2), non-arid fraction
      REAL*4,  INTENT(IN)  :: LANDFRAC(:) ! Fraction of landtypes

      !Input parameters for the canopy reduction factor
      REAL*8,  INTENT(IN)  :: SUNCOS      ! Cosine of solar zenith angle
      REAL*8,  INTENT(IN)  :: U10M        ! E/W wind speed @ 10m altitude [m]
      REAL*8,  INTENT(IN)  :: V10M        ! N/S wind speed @ 10m altitude [m]
      REAL*8,  INTENT(IN)  :: R_CANOPY(:) ! Resist. of canopy to NOx [1/s]
      REAL*8,  INTENT(IN)  :: LAI         ! Leaf area indices [cm2/cm2]
!
! !OUTPUT PARAMETERS:
!
      REAL*8,  INTENT(OUT) :: SOILNOx     ! Soil NOx emissions [molec/cm2/s]
      REAL*4,  INTENT(OUT) :: GWET_PREV   ! Soil Moisture Prev timestep
      REAL*4,  INTENT(OUT) :: DRYPERIOD   ! Dry period length in hours
      REAL*4,  INTENT(OUT) :: PFACTOR     ! Pulsing Factor
      REAL*8,  INTENT(OUT) :: FERTDIAG    ! Fert emissions [molec/cm2/s
      INTEGER, INTENT(OUT) :: RC          ! Return code
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
      REAL*8  :: WET_TERM,  A_FERT    

      !=================================================================
      ! Initialize
      !=================================================================

      ! Assume successful return
      RC             = GIGC_SUCCESS

      ! Initialize
      SOILNOX        = 0d0
      FERTDIAG       = 0d0
 
      ! Surface temperature [C]
      TC             = TK - 273.15d0

      ! Surface wind speed, squared
      WINDSQR        = U10M**2 + V10M**2

      !=================================================================
      ! Compute soil NOx emissions
      !=================================================================

      ! Cumulative multiplication factor (over baseline emissions) 
      ! that accounts for soil pulsing
      PULSE = PULSING( GWET, TS_EMIS, GWET_PREV, PFACTOR, DRYPERIOD )

      ! ------Loop Over MODIS/Koppen  Landtypes
      DO K = 1, 24

         ! Temperature-dependent term of soil NOx emissions [unitless]
         ! Use GWET instead of climo wet/dry
         TEMP_TERM = SOILTEMP( K , TC, GWET)

         ! Soil moisture scaling of soil NOx emissions 
         WET_TERM = SOILWET( GWET , CLIM(:))

         ! Fertilizer emission  
         A_FERT = FERTADD( SOILFRT , DEPN)

	 ! Canopy reduction factor
         CRF_TERM  = SOILCRF( K, LAI, 
     &                        R_CANOPY(K),
     & 			      WINDSQR, SUNCOS )
       
         ! SOILNOX includes fertilizer
         SOILNOX   = (SOILNOX 
     &             + ( A_BIOME(K) + A_FERT )
     &             * ( TEMP_TERM * WET_TERM * PULSE ) 
     &             * LANDFRAC(K)
     &             * ( 1.d0 - CRF_TERM  ) )

         ! FERTDIAG, only used for the fertilizer diagnostic
         FERTDIAG  = (FERTDIAG
     &             + ( A_FERT )
     &             * ( TEMP_TERM * WET_TERM * PULSE ) 
     &             * LANDFRAC(K)
     &             * ( 1.d0 - CRF_TERM  ) )

      ENDDO

      END SUBROUTINE SOIL_NOX_EMISSION

!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
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
 
      ! Convert soil NOx  emissions to [molec/cm2/s]
      SOIL_TEMP = SOIL_TEMP * UNITCONV 

      END FUNCTION SOILTEMP
!EOC
!----------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
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
      FUNCTION SOILWET( GWET , CLIM) RESULT( WETSCALE )
!
! !INPUT PARAMETERS: 
!
      ! Fraction of arid & non-arid soil in the gridbox          
      REAL*4, INTENT(IN) :: CLIM(:)       

      ! Top soil wetness [unitless]
      REAL*8, INTENT(IN) :: GWET 
    
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
      IF ( CLIM(2) .GE. CLIM(1) .AND. CLIM(2) .NE. 0) THEN 
         !Arid, max Poisson = 0.2
         WETSCALE = 8.24 * GWET * EXP(-12.5*GWET*GWET)
      ELSE
         !Non-arid, max Poisson = 0.3
         WETSCALE = 5.5 * GWET * EXP( -5.55 * GWET * GWET) 
      ENDIF        

      END FUNCTION SOILWET
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
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
      FUNCTION SOILCRF( K,        LAI, 
     &                  CANOPYNOX, WINDSQR, 
     &                  SUNCOS              ) RESULT( SOIL_CRF )
!
! !INPUT PARAMETERS: 
!
      INTEGER, INTENT(IN) :: K          ! Soil biome type
      REAL*8,  INTENT(IN) :: LAI        ! Leaf area index [cm2/cm2]
      REAL*8,  INTENT(IN) :: CANOPYNOX  ! Bulk sfc resistance to NOx [1/s]
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
      IF ( LAI > 0d0 .and. CANOPYNOX > 0d0 ) THEN

         ! Adjust the ventilation velocity.  
         ! NOTE: SOILEXC(21) is the canopy wind extinction 
         ! coefficient for the tropical rainforest biome.
         VFNEW    = (VFNEW * SQRT( WINDSQR/9d0 * 7d0/LAI     ) *
     &                          ( SOILEXC(21)  / SOILEXC(K) ))

         ! Soil canopy reduction factor
         SOIL_CRF = CANOPYNOX / ( CANOPYNOX + VFNEW )

      ELSE
     
         ! Otherwise set the soil canopy reduction factor to zero
         SOIL_CRF = 0d0

      ENDIF

      END FUNCTION SOILCRF
!EOC
!-----------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
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
      FUNCTION FERTADD( SOILFERT, DEPN) RESULT( FERT_ADD )
!
! !INPUT PARAMETERS: 
!
      REAL*8, INTENT(IN) :: DEPN      ! N emissions from deposition
      REAL*8, INTENT(IN) :: SOILFERT  ! N emissions from fertilizers
                                      !  read in from disk and passed
                                      !  here as an argument [ng N/m2/s]
!
! !RETURN_VALUE:
! 
      REAL*8            :: FERT_ADD   ! Total Fert emissions

      REAL*8, PARAMETER :: SECPERYEAR = 86400.d0 * 365.
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

      ! Soil fert and dep [ ng N /m2 ], a measure of N avail. in soil
      FERT_ADD = SOILFERT  + DEPN

      ! Convert [ng N/m2] --> [ng N /m2/s]
      ! (scale needed to force fert emiss of 1.8 Tg N/yr w/o canopy uptake)
      FERT_ADD = FERT_ADD / SECPERYEAR * FERT_SCALE   

      END FUNCTION FERTADD
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
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
      FUNCTION PULSING( GWET,      TS_EMIS, 
     &                  GWET_PREV, PFACTOR,
     &                  DRYPERIOD         ) RESULT( THE_PULSING )
!
! !INPUT PARAMETERS: 
!
      REAL*8, INTENT(IN)    :: GWET        ! Soil Moisture 
      REAL*8, INTENT(IN)    :: TS_EMIS     ! Emissions timestep [min]

! !INPUT/OUTPUT PARAMETERS:
!
      REAL*4, INTENT(INOUT) :: GWET_PREV   ! soil moisture from prev. timestep
      REAL*4, INTENT(INOUT) :: PFACTOR     ! pulsing factor
      REAL*4, INTENT(INOUT) :: DRYPERIOD   ! dry period in # timesteps
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
!     PFACTOR = 13.01 ln ( DRYPERIOD ) -  53.6
! 
!  where PFACTOR is the magnitude of peak flux relative to prewetting flux, 
!  and DRYPERIOD  is the length of the antecedent dry period in hours.
! 
!  The pulse decays with 
! 
!     PFACTOR = PFACTOR * EXP( -0.068d0 * DTSRCE )       
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
!  28 Oct 2013 - R. Yantosca - Bug fix: prevent log(0) from happening 
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

      ! Emission timestep [min --> hours]
      DTSRCE = TS_EMIS / 60d0

      ! If soil moisture less than 0.3 and no pulse is taking place
      IF ( GWET < 0.3D0 .and. PFACTOR == 1.D0) THEN

         ! Get change in soil moisture since previous timestep
         GDIFF = ( GWET - GWET_PREV )

         ! If change in soil moisture is > 0.01 (rains)
         IF ( GDIFF > 0.01 ) THEN

!-----------------------------------------------------------------------------
! Prior to 10/28/13:
! Now add a numerical trap to prevent log(0) from occurring, which would
! cause a floating-point error. (mpayer, bmy, 10/28/13)
!            !Initialize new pulse factor (dry period hours)
!            PFACTOR = 13.01 * LOG( DRYPERIOD ) - 53.6
!-----------------------------------------------------------------------------
            ! Initialize new pulse factor (dry period hours)
            IF ( DRYPERIOD > 0 ) THEN
               PFACTOR = 13.01 * LOG( DRYPERIOD ) - 53.6
            ELSE
               PFACTOR = -53.6
            ENDIF

            ! If dry period < ~3 days then no pulse
            IF ( PFACTOR < 1.0 ) PFACTOR = 1.0

            ! Reinitialize dry period
            DRYPERIOD = 0

         ! If no rain (i.e.,  change in soil moisture is < 0.01)
         ELSE

            ! Add one timestep to dry period
            DRYPERIOD = DRYPERIOD + DTSRCE

         ENDIF

      ! If box is already pulsing , then decay pulse one timestep
      ELSEIF ( PFACTOR /= 1.d0) THEN

         ! Decay pulse
         PFACTOR   = PFACTOR * EXP( -0.068d0 * DTSRCE )

         ! Update dry period
         IF ( GWET < 0.3D0 ) DRYPERIOD = DRYPERIOD + DTSRCE

         ! If end of pulse
         IF ( PFACTOR < 1.d0 ) PFACTOR = 1.d0
      
      ENDIF

      ! Update soil moisture holder for previous timestep
      GWET_PREV = GWET

      ! Return the pulsing factor
      THE_PULSING = PFACTOR

      END FUNCTION PULSING
!EOC
      END MODULE SOIL_NOx_MOD

