! $Id: mercury_mod.f,v 1.24 2009/09/01 19:21:18 cdh Exp $
      MODULE LAND_MERCURY_MOD

!
!******************************************************************************
!  Module MERCURY_MOD contains variables and routines for the GEOS-CHEM 
!  mercury simulation. (eck, bmy, 12/14/04, 4/6/06)
!
!  Module Variables:
!  ============================================================================
!  (1 ) TRANSP      (REAL*8 ) : Plant transpiration rate [m/s]
!  (2 ) Hg0dryGEOS  (REAL*8 ) : Dry dep. Hg0 from restart file [kg/s]
!  (3 ) HgIIdryGEOS (REAL*8 ) : Dry dep. HgII from restart file [kg/s]
!  (4 ) HgIIwetGEOS (REAL*8 ) : Wet dep. HgII from restart file [kg/s]
!
!  Module Routines:
!  ===========================================================================
!  (1 ) BIOMASSHG            : Wildfire Hg emissions
!  (2 ) VEGEMIS              : Transpiration Hg emissions 
!  (3 ) SOILEMIS             : Soil Hg emissions
!  (4 ) LAND_MERCURY_FLUX    : Land model for Hg
!  (5 ) READ_NASA_TRANSP     : Read transpiration rates from file
!  (6 ) GTMM_DR              : Driver routine for GTMM code.
!  (7 ) MAKE_GTMM_RESTART    : Write restart file for deposition for GTMM
!  (8 ) READ_GTMM_RESTART    : Read restart file for deposition for GTMM
!  (9 ) INIT_LAND_MERCURY    : Allocates and zeroes all module arrays
!  (10) CLEANUP_LAND_MERCURY : Deallocates all module arrays
!
!  GEOS-CHEM modules referenced by land_mercury_mod.f
!  ============================================================================
!  (1 ) biomass_mod.f      : Wrapper for biomass emissions
!  (2 ) bpch2_mod.f        : Module to read bpch files
!  (3 ) dao_mod.f          : Module w/ arrays for DAO met fields
!  (4 ) error_mod.f        : Module for catching errors
!  (5 ) file_mod.f         : Module w/ file unit numbers
!  (6 ) grid_mod.f         : Module w/ horizontal grid information
!  (7 ) lai_mod.f          : Module w/ routines to read and store AVHRR LAI
!  (8 ) logical_mod.f      : Module w/ GEOS-CHEM logical switches
!  (9 ) ocean_mercury_mod.f: Module w/ routines to compute oceanic Hg fluxes
!  (10) time_mod.f         : Module w/ routines to compute date & time
!  (11) tracer_mod.f       :
!  (12) tracerid_mod.f     : Module w/ pointers to tracers & emissions
!  (13) transfer_mod.f     : Module w/ routines to copy data from REAL*4 to
!                            REAL*8
!
!  Nomenclature: 
!  ============================================================================
!  (1 ) Hg(0)  a.k.a. Hg0  : Elemental   mercury
!  (2 ) Hg(II) a.k.a. Hg2  : Divalent    mercury
!  (3 ) HgP                : Particulate mercury
!
!  Mercury Tracers (1-3 are always defined; 4-21 are defined for tagged runs)
!  ============================================================================
!  (1 ) Hg(0)              : Hg(0)  - total tracer
!  (2 ) Hg(II)             : Hg(II) - total tracer 
!  (3 ) HgP                : HgP    - total tracer
!  ------------------------+---------------------------------------------------
!  (4 ) Hg(0)_an           : Hg(0)  - North American anthropogenic
!  (5 ) Hg(0)_ae           : Hg(0)  - European Anthropogenic
!  (6 ) Hg(0)_aa           : Hg(0)  - Asian anthropogenic
!  (7 ) Hg(0)_ar           : Hg(0)  - Rest of World Anthropogenic
!  (8 ) Hg(0)_oc           : Hg(0)  - Ocean emission
!  (9 ) Hg(0)_ln           : Hg(0)  - Land reemission
!  (10) Hg(0)_nt           : Hg(0)  - Land natural emission
!  ------------------------+---------------------------------------------------
!  (11) Hg(II)_an          : Hg(II) - North American anthropogenic
!  (12) Hg(II)_ae          : Hg(II) - European Anthropogenic
!  (13) Hg(II)_aa          : Hg(II) - Asian anthropogenic
!  (14) Hg(II)_ar          : Hg(II) - Rest of World Anthropogenic
!  (15) Hg(II)_oc          : Hg(II) - Ocean emission
!  (16) Hg(II)_ln          : Hg(II) - Land reemission
!  (17) Hg(II)_nt          : Hg(II) - Land natural emission
!  ------------------------+---------------------------------------------------
!  (18) HgP_an             : HgP    - North American anthropogenic
!  (19) HgP_ae             : HgP    - European anthropogenic
!  (20) HgP_aa             : HgP    - Asian anthropogenic
!  (21) HgP_ar             : HgP    - Rest of world anthropogenic
!  ------------------------+---------------------------------------------------
!  (22) HgP_oc             : HgP    - Ocean emission        (FOR FUTURE USE)
!  (23) HgP_ln             : HgP    - Land reemission       (FOR FUTURE USE)
!  (24) HgP_nt             : HgP    - Land natural emission (FOR FUTURE USE)
!
!  References:
!  ============================================================================
!  (3 ) Selin, N., et al. (2007). "Chemical cycling and deposition of 
!       atmospheric mercury: Global constraints from observations." 
!       J. Geophys. Res. 112.
!  (4 ) Selin, N., et al. (2008). "Global 3-D land-ocean-atmospehre model
!       for mercury: present-day versus preindustrial cycles and
!       anthropogenic enrichment factors for deposition." Global
!       Biogeochemical Cycles 22: GB2011.
!  (5 ) Allison, J.D. and T.L. Allison (2005) "Partition coefficients for
!       metals in surface water, soil and waste." Rep. EPA/600/R-05/074,
!       US EPA, Office of Research and Development, Washington, D.C.
!  (6 ) Mintz, Y and G.K. Walker (1993). "Global fields of soil moisture
!       and land surface evapotranspiration derived from observed
!       precipitation and surface air temperature." J. Appl. Meteorol. 32 (8), 
!       1305-1334.
!
!  NOTES:
!  (1 ) Move all land emissions routine for mercury to this new module.
!       (ccc 9/10/09)
!******************************************************************************
!
      IMPLICIT NONE

      !=================================================================
      ! MODULE PRIVATE DECLARATIONS -- keep certain internal variables 
      ! and routines from being seen outside "mercury_mod.f"
      !=================================================================

      ! Make everything PRIVATE ...
      PRIVATE

      ! ... except these routines
      PUBLIC :: BIOMASSHG
      PUBLIC :: VEGEMIS
      PUBLIC :: SOILEMIS
      PUBLIC :: LAND_MERCURY_FLUX
      PUBLIC :: SNOWPACK_MERCURY_FLUX
      PUBLIC :: INIT_LAND_MERCURY
      PUBLIC :: CLEANUP_LAND_MERCURY

      ! ... except these variables
      PRIVATE :: TRANSP

      !=================================================================
      ! MODULE VARIABLES
      !=================================================================
      REAL*8,  ALLOCATABLE :: TRANSP(:,:)
      
      !=================================================================
      ! MODULE ROUTINES -- follow below the "CONTAINS" statement
      !=================================================================
      CONTAINS

!----------------------------------------------------------------------------
      SUBROUTINE LAND_MERCURY_FLUX( LFLUX, LHGSNOW )
!
!******************************************************************************
!  Subroutine LAND_MERCURY_FLUX calculates emissions of Hg(0) from 
!  prompt recycling of previously deposited mercury to land, in [kg/s].  
!  (eck, cdh, eds, 7/30/08)
!  
!  Arguments as Output
!  ============================================================================
!  (1 ) LFLUX (REAL*8) : Flux of Hg(0) from the promptly recycled land emissions!                        [kg/s]
!
!  NOTES:
!  (1 ) Now uses SNOWMAS from DAO_MOD for compatibility with GEOS-5.
!       (eds 7/30/08)
!  (2 ) Now includes REEMFRAC in parallelization; previous versions may have
!       overwritten variable. (cdh, eds 7/30/08)
!  (3 ) Now also reemit Hg(0) from ice surfaces, including sea ice 
!       (cdh, 8/19/08)
!******************************************************************************
!
      USE TRACERID_MOD,  ONLY : ID_Hg0,          N_Hg_CATS
      USE LOGICAL_MOD,   ONLY : LSPLIT
      USE TIME_MOD,      ONLY : GET_TS_EMIS
      USE DAO_MOD,       ONLY : SNOW, SNOMAS 
!      USE OCEAN_MERCURY_MOD, ONLY : WD_HGP, WD_HG2, DD_HGP, DD_HG2
      USE DEPO_MERCURY_MOD, ONLY : WD_HGP, WD_HG2, DD_HGP, DD_HG2
      USE DAO_MOD,       ONLY : IS_ICE, IS_LAND


#     include "CMN_SIZE"      ! Size parameters

      ! Arguments 
      REAL*8,  INTENT(OUT)  :: LFLUX(IIPAR,JJPAR,N_Hg_CATS)
      LOGICAL, INTENT(IN)   :: LHGSNOW
       
  
      REAL*8                :: DTSRCE, REEMFRAC, SNOW_HT
      REAL*8, PARAMETER     :: SEC_PER_YR = 365.25d0 * 86400d0
      INTEGER               :: I,      J,      NN

      !=================================================================
      ! LAND_MERCURY_FLUX begins here!
      !=================================================================

      ! Emission timestep [s]
      DTSRCE = GET_TS_EMIS() * 60d0     

!$OMP PARALLEL DO 
!$OMP+DEFAULT( SHARED ) 
!$OMP+PRIVATE( I, J, NN )
!$OMP+PRIVATE( REEMFRAC, SNOW_HT )  
      DO J  = 1, JJPAR
      DO I  = 1, IIPAR
      DO NN = 1, N_Hg_CATS
    
#if defined( GEOS_5 )
         ! GEOS5 snow height (water equivalent) in mm. (Docs wrongly say m)
         SNOW_HT = SNOMAS(I,J)
#else
         ! GEOS1-4 snow heigt (water equivalent) in mm
         SNOW_HT = SNOW(I,J)
#endif 
        
         ! If snow > 1mm on the ground, reemission fraction is 0.6,
         ! otherwise 0.2
         IF ( (SNOW_HT > 1D0) .OR. (IS_ICE(I,J)) ) THEN
            ! If snowpack model is on, then we don't do rapid reemission
            IF (LHGSNOW) THEN
               REEMFRAC=0d0
            ELSE
               REEMFRAC=0.6d0
            ENDIF 
         ELSE
            REEMFRAC=0.2d0
         ENDIF
         

         IF ( IS_LAND(I,J) .OR. IS_ICE(I,J) ) THEN 
            
            ! Mass of emitted Hg(0), kg
            LFLUX(I,J,NN) =
     &           ( WD_HgP(I,J,NN)+
     &           WD_Hg2(I,J,NN)+
     &           DD_HgP(I,J,NN)+
     &           DD_Hg2(I,J,NN) ) * REEMFRAC
            
            ! Emission rate of Hg(0). Convert kg /timestep -> kg/s
            LFLUX(I,J,NN) = LFLUX(I,J,NN) / DTSRCE
             
         ELSE
         
            ! No flux from non-land surfaces (water, sea ice)
            LFLUX(I,J,NN) = 0D0
         
         ENDIF

      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO
     
      ! Return to calling program
      END SUBROUTINE LAND_MERCURY_FLUX


!-----------------------------------------------------------------------------

      SUBROUTINE BIOMASSHG( EHg0_bb )
!
!******************************************************************************
!  Subroutine BIOMASSHG is the subroutine for Hg(0) emissions from biomass
!  burning. These emissions are active only for present day simulations and
!  not for preindustrial simulations (eck, cdh, eds, 7/30/08)
!
!  Emissions are based on an inventory of CO emissions from biomass burning 
!  (Duncan et al. J Geophys Res 2003), multiplied by a Hg/CO ratio in BB plumes
!  from Franz Slemr (Poster, EGU 2006).
!
!  Slemr surveyed emission factors from measurements worldwide. Although his
!  best estimate was 1.5e-7 mol Hg/ mol CO, we chose the highest value
!  (2.1e-7 mol Hg/ mol CO) in the range because the simulations shown in
!  Selin et al. (GBC 2008) required large Hg(0) emissions to sustain
!  reasonable atmospheric Hg(0) concentrations. (eck, 11/13/2008)   
!
!******************************************************************************
!     

      ! References to F90 modules
! IDBCO moved from BIOMASS_MOD to TRACERID_MOD. (ccc, 5/6/10)
!      USE BIOMASS_MOD,    ONLY: BIOMASS, IDBCO
      USE BIOMASS_MOD,    ONLY: BIOMASS
      USE TRACERID_MOD,    ONLY: IDBCO
      USE LOGICAL_MOD,    ONLY: LBIOMASS, LPREINDHG
      USE TIME_MOD,       ONLY: GET_TS_EMIS
      USE GRID_MOD,       ONLY: GET_AREA_CM2

#     include "CMN_SIZE"     ! Size parameters
#     include "CMN_DIAG"     ! Diagnostic arrays & switches

      !Arguments
      REAL*8, DIMENSION(:,:),INTENT(OUT) :: EHg0_bb
      
      ! Local variables
      REAL*8                 :: DTSRCE, E_CO, AREA_CM2
      INTEGER                :: I, J

      ! Hg molar mass, kg Hg/ mole Hg
      REAL*8,  PARAMETER   :: FMOL_HG     = 200.59d-3

      ! Hg/CO molar ratio in BB emissions, mol/mol
      ! emission factor 1.5e-7 molHg/molCO (Slemr et al poster EGU 2006)
      ! change emission factor to 2.1
      REAL*8,  PARAMETER   :: BBRatio_Hg_CO = 2.1D-7

      ! External functions
      REAL*8,  EXTERNAL      :: BOXVL

      !=================================================================
      ! BIOMASSHG begins here!
      !=================================================================

      ! DTSRCE is the number of seconds per emission timestep
      DTSRCE = GET_TS_EMIS() * 60d0

      ! Do biomass Hg emissions if biomass burning is on and it is a 
      ! present-day simulation (i.e. not preindustrial)
      IF ( LBIOMASS .AND. ( .NOT. LPREINDHG ) ) THEN

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( E_CO, I, J, AREA_CM2 )
         DO J = 1, JJPAR
         DO I = 1, IIPAR
 
            ! Grid box surface area, cm2
            AREA_CM2 = GET_AREA_CM2( J )

            ! Convert molec CO /cm3 /s -> mol CO /gridbox /s 
            E_CO = ( BIOMASS(I,J,IDBCO) / 6.022D23 ) * AREA_CM2

            ! Convert mol CO /gridbox /s to kg Hg /gridbox /s
            EHg0_bb(I,J) = E_CO * BBRatio_Hg_CO * FMOL_HG 
         
         ENDDO
         ENDDO
!$OMP END PARALLEL DO

      ELSE

         ! No emissions for preindustrial period, or when BB is turned off.
         EHg0_bb = 0D0

      ENDIF


      END SUBROUTINE BIOMASSHG

!-----------------------------------------------------------------------------

      SUBROUTINE VEGEMIS( LGCAPEMIS, EHg0_dist, EHg0_vg )
!
!******************************************************************************
!  Subroutine VEGEMIS is the subroutine for Hg(0) emissions from vegetation 
!  by evapotranspiration.
!  (eck, cdh, eds, 7/30/08)
!
!  Vegetation emissions are proportional to the evapotranspiration rate and the
!  soil water mercury content. We assume a constant concentration of mercury
!  in soil matter, based on the preindustrial and present-day simulations
!  described in Selin et al. (GBC 2008) and in SOILEMIS subroutine. From the
!  soil matter Hg concentration, we calculate a soil water Hg concentration in 
!  equilibrium (Allison and Allison, 2005).
!  NASA provides a climatology of evapotranspiration based on a water budget
!  model (Mintz and Walker, 1993).
!
! Calculate vegetation emissions following Xu et al (1999)
!    Fc = Ec Cw
!
!    Fc is Hg0 flux (ng m-2 s-1)
!    Ec is canopy transpiration (m s-1)
!    Cw is conc of Hg0 in surface soil water (ng m-3)
!
! Calculate Cw from the Allison and Allison (2005) equilibrium formula
!    Cw = Cs / Kd
!
!    Cs is the concentration of Hg is surface soil solids, ng/g
!    Kd is the equilibrium constant = [sorbed]/[dissolved]
!       log Kd = 3.8 L/kg -> Kd = 6310 L /kg = 6.31D-3 m3/g
!
! We assume a global mean Cs = 45 ng/g for the preindustrial period. In
! iterative simulations we redistribute this according to the deposition
! pattern while maintining the global mean. The scaling factor, EHg0_dist,
! also accounts for the anthropogenic enhancement of soil Hg in the present 
! day. 
!
!******************************************************************************
!       

      ! References to F90 modules      
      USE DAO_MOD,        ONLY: RADSWG, IS_LAND
      USE TIME_MOD,       ONLY: GET_MONTH, ITS_A_NEW_MONTH
      USE TIME_MOD,       ONLY: GET_TS_EMIS
      USE GRID_MOD,       ONLY: GET_AREA_M2

#     include "CMN_SIZE"     ! Size parameters
#     include "CMN_DEP"      ! FRCLND

      ! Arguments
      LOGICAL,                INTENT(IN)  :: LGCAPEMIS
      REAL*8, DIMENSION(:,:), INTENT(IN)  :: EHg0_dist
      REAL*8, DIMENSION(:,:), INTENT(OUT) :: EHg0_vg

      ! Local Variables
      REAL*8             :: DRYSOIL_HG, SOILWATER_HG, AREA_M2, VEG_EMIS
      INTEGER            :: I, J

      ! Soil Hg sorption to dissolution ratio, m3/g
      REAL*8, PARAMETER  :: Kd = 6.31D-3

      ! Preindustrial global mean soil Hg concentration, ng Hg /g dry soil
      REAL*8, PARAMETER  :: DRYSOIL_PREIND_HG = 45D0

      !=================================================================
      ! VEGEMIS begins here!
      !=================================================================

      ! No emissions through transpiration if we use Bess' GCAP emissions
      IF (LGCAPEMIS) THEN

         EHg0_vg = 0D0

      ELSE

         ! read GISS TRANSP monthly average
         IF ( ITS_A_NEW_MONTH() ) THEN 
            CALL READ_NASA_TRANSP
         ENDIF 

         ! loop over I,J
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, SOILWATER_HG, DRYSOIL_HG, VEG_EMIS, AREA_M2 ) 
         DO J=1, JJPAR
         DO I=1, IIPAR
        
            IF (IS_LAND(I,J)) THEN  

               ! Dry soil Hg concentration, ng Hg /g soil
               DRYSOIL_HG = DRYSOIL_PREIND_HG * EHg0_dist(I,J)

               ! Hg concentration in soil water, ng /m3
               SOILWATER_HG =  DRYSOIL_HG / Kd

               ! Emission from vegetation, ng /m2
               VEG_EMIS = SOILWATER_HG * TRANSP(I,J)

               ! convert from ng /m2 /s -> kg/gridbox/s
               ! Grid box surface area [m2]
               AREA_M2      = GET_AREA_M2( J )
               EHg0_vg(I,J) = VEG_EMIS * AREA_M2 * 1D-12

            ELSE

               ! No emissions from water and ice
               EHg0_vg(I,J) = 0D0

            ENDIF

         ENDDO
         ENDDO
!$OMP END PARALLEL DO

      ENDIF

      END SUBROUTINE VEGEMIS

!------------------------------------------------------------------------------
 
      SUBROUTINE SOILEMIS( EHg0_dist, EHg0_so )
!
!******************************************************************************
!  Subroutine SIOLEMIS is the subroutine for Hg(0) emissions from soils.
!  (eck, eds, 7/30/08)
!  
!  Soil emissions are a function of solar radiation at ground level 
!  (accounting for attenuation by leaf canopy) and surface temperature. 
!  The radiation dependence from Zhang et al. (2000) is multiplied by the 
!  temperature dependence from Poissant and Casimir (1998). 
!  Finally, this emission factor is multiplied by the soil mercury
!  concentration and scaled to meet the global emission total.
!
!  Comments on soil Hg concentration:
!  We chose the preindustrial value of 45 ng Hg /g dry soil as the mean of
!  the range quoted in Selin et al. (GBC 2008): 20-70 ng/g (Andersson, 1967; 
!  Shacklette et al., 1971; Richardson et al., 2003; Frescholtz and Gustin,
!  2004). Present-day soil concentrations are thought to be 15% greater than
!  preindustrial (Mason and Sheu 2002), but such a difference is much less
!  than the range of concentrations found today, so not well constrained.
!  We calculate the present-day soil Hg distribution by adding a global mean
!  6.75 ng/g (=0.15 * 45 ng/g) according to present-day Hg deposition.
!  (eck, 11/13/08)
!
!  Notes
!  (1 ) Added comments. (cdh, eds, 7/30/08)
!  (2 ) Now include light attenuation by the canopy after sunset. Emissions
!       change by < 1% in high-emission areas  (cdh, 8/13/2008)
!  (3 ) Removed FRCLND for consistency with other Hg emissions (cdh, 8/19/08) 
!******************************************************************************
!

      ! References to F90 modules      
      USE LAI_MOD,        ONLY: ISOLAI, MISOLAI, PMISOLAI, DAYS_BTW_M
      USE DAO_MOD,        ONLY: RADSWG, SUNCOS, TS, IS_LAND
      USE TIME_MOD,       ONLY: GET_MONTH, ITS_A_NEW_MONTH
      USE TIME_MOD,       ONLY: GET_TS_EMIS
      USE GRID_MOD,       ONLY: GET_AREA_M2
      USE DAO_MOD,        ONLY: SNOW, SNOMAS

#     include "CMN_SIZE"     ! Size parameters
#     include "CMN_DEP"      ! FRCLND

      ! Arguments
      REAL*8, DIMENSION(:,:), INTENT(IN) :: EHg0_dist
      REAL*8, DIMENSION(:,:), INTENT(OUT):: EHg0_so

      ! Local variables
      REAL*8             :: SOIL_EMIS, DIMLIGHT, TAUZ, LIGHTFRAC
      REAL*8             :: AREA_M2, DRYSOIL_HG, SNOW_HT
      INTEGER            :: I, J, JLOOP

      ! Preindustrial global mean soil Hg concentration, ng Hg /g dry soil
      REAL*8, PARAMETER  :: DRYSOIL_PREIND_HG = 45D0

      ! Scaling factor for emissions, g soil /m2 /h
      ! (This parameter is beta in Eq 3 of Selin et al., GBC 2008.
      ! The value in paper is actually DRYSOIL_PREIND_HG * SOIL_EMIS_FAC 
      ! and the stated units are incorrect. The paper should have stated
      ! beta = 1.5D15 / 45D0 = 3.3D13)
      ! This parameter is tuned in the preindustrial simulation 
      ! so that total deposition to soil equals total emission from soil,
      ! while also requiring global mean soil Hg concentration of 45 ng/g 
!      REAL*8, PARAMETER  :: SOIL_EMIS_FAC = 3.3D13
      REAL*8, PARAMETER  :: SOIL_EMIS_FAC = 2.4D-2 ! for sunlight function

      REAL*8              :: SUNCOSVALUE
      !=================================================================
      ! SOILEMIS begins here!
      !=================================================================

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, SOIL_EMIS, JLOOP,      SNOW_HT              ) 
!$OMP+PRIVATE( DRYSOIL_HG, TAUZ, LIGHTFRAC, AREA_M2, SUNCOSVALUE ) 
      DO J=1, JJPAR
      DO I=1, IIPAR
         
#if defined( GEOS_5 )
         ! GEOS5 snow height (water equivalent) in mm. (Docs wrongly say m)
         SNOW_HT = SNOMAS(I,J)
#else
         ! GEOS1-4 snow heigt (water equivalent) in mm
         SNOW_HT = SNOW(I,J)
#endif          
         
         IF ( IS_LAND(I,J) .AND. (SNOW_HT < 1d0) ) THEN     

            ! 1-D grid box index for SUNCOS
            JLOOP = ( (J-1) * IIPAR ) + I
         
            ! attenuate solar radiation based on function of leaf area index
            ! Jacob and Wofsy 1990 equations 8 & 9
            TAUZ = ISOLAI(I,J) * 0.5D0

            ! For very low and below-horizon solar zenith angles, use
            ! same attenuation as for SZA=85 degrees
            SUNCOSVALUE = MAX( SUNCOS(JLOOP), 0.09D0 )

            ! fraction of light reaching the surface is
            ! attenuated based on LAI
            LIGHTFRAC = EXP( -TAUZ / SUNCOSVALUE )

!            ! if there is sunlight
!            IF (SUNCOS(JLOOP) > 0d0 .and. RADSWG(I,J) > 0d0 ) THEN
!
!               ! attenuate solar radiation based on function of leaf area index
!               ! Jacob and Wofsy 1990 equations 8 & 9
!               TAUZ = ISOLAI(I,J) * 0.5D0
!
!               ! fraction of light reaching the surface is
!               ! attenuated based on LAI
!               LIGHTFRAC = EXP( -TAUZ / SUNCOS(JLOOP) )
!
!            ELSE
!
!               ! If the sun has set, then set the canopy attenuation to
!               ! the same as for a high solar zenith angle, 80 deg
!               LIGHTFRAC = EXP( -TAUZ / 0.17D0 )
!
!            ENDIF

            ! Dry soil Hg concentration, ng Hg /g soil
            DRYSOIL_HG = DRYSOIL_PREIND_HG * EHg0_dist(I,J)           

            ! Soil emissions, ng /m2 /h
            ! includes temperature and solar radiation effects
!            SOIL_EMIS = EXP( 1000D0 / TS(I,J) * -10.548D0 ) * 
!     &           EXP( 0.0011 * RADSWG(I,J) * LIGHTFRAC ) *
!     &           DRYSOIL_HG * SOIL_EMIS_FAC

! CDH try formula with just light dependence 10/18/2009
            SOIL_EMIS = 
     &           EXP( 0.0011 * RADSWG(I,J) * LIGHTFRAC ) *
     &           DRYSOIL_HG * SOIL_EMIS_FAC
     
            ! Grid box surface area [m2]
            AREA_M2   = GET_AREA_M2( J )
 
            ! convert soilnat from ng /m2 /h -> kg /gridbox /s
            EHg0_so(I,J) = SOIL_EMIS * AREA_M2 * 1D-12 / ( 60D0 * 60D0 )

         ELSE

            ! no soil emissions from water and ice
            EHg0_so(I,J) = 0D0
        
         ENDIF
        
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      WRITE(6,'(G12.3)') SUM(EHG0_SO)
    
      END SUBROUTINE SOILEMIS

!-----------------------------------------------------------------------------

      SUBROUTINE READ_NASA_TRANSP
!
!******************************************************************************
!  Subroutine READ_NASA_TRANSP reads monthly average transpirtation from NASA
!  http://gcmd.nasa.gov/records/GCMD_MINTZ_WALKER_SOIL_AND_EVAPO.html
!  for input into the vegetation emissions. (eck, 9/15/06)
!
!       Mintz, Y and G.K. Walker (1993). "Global fields of soil moisture
!       and land surface evapotranspiration derived from observed
!       precipitation and surface air temperature." J. Appl. Meteorol. 32 (8), 
!       1305-1334.
! 
!  Arguments as Input/Output:
!  ============================================================================
!  (1 ) TRANSP  : Transpiration [m/s]
!
!******************************************************************************
!

      ! References to F90 modules     
      USE TIME_MOD,       ONLY : GET_MONTH,  ITS_A_NEW_MONTH
      USE BPCH2_MOD,      ONLY : GET_TAU0,   READ_BPCH2
      USE TRANSFER_MOD,   ONLY : TRANSFER_2D

#     include "CMN_SIZE"      ! Size parameters

      ! Local variables
      INTEGER             :: I, J, L, MONTH, N
      REAL*4              :: ARRAY(IGLOB,JGLOB,1)
      REAL*8              :: XTAU
      CHARACTER(LEN=255)  :: FILENAME
      CHARACTER(LEN=2)    :: CMONTH(12) = (/ '01','02','03','04',
     &                                       '05','06','07','08',
     &                                       '09','10','11','12'/)
  
      !=================================================================
      ! READ_NASA_TRANSP begins here!
      !=================================================================

      ! Get the current month
      MONTH = GET_MONTH() 
     
!      FILENAME='/as/home/eck/transp/nasatransp_4x5.'
!     &        //CMONTH(MONTH)//'.bpch'
      FILENAME='/home/eck/emissions/transp/nasatransp_4x5.'
     &        //CMONTH(MONTH)//'.bpch'

      XTAU     = GET_TAU0(MONTH, 1, 1995 )

      ! Echo info
      WRITE( 6, 100 ) TRIM( FILENAME )
 100  FORMAT( '     - TRANSP_NASA: Reading ', a )     
 
      CALL READ_BPCH2( FILENAME, 'TRANSP-$', 1, 
     &                 XTAU,      IGLOB,     JGLOB,    
     &                 1,         ARRAY,   QUIET=.TRUE. )

      CALL TRANSFER_2D( ARRAY(:,:,1), TRANSP )
      
      ! convert from mm/month to m/s
      TRANSP = TRANSP * 1D-3 * 12D0 / ( 365D0 * 24D0 * 60D0 * 60D0 )
      
      END SUBROUTINE READ_NASA_TRANSP

!-----------------------------------------------------------------------------


      SUBROUTINE SNOWPACK_MERCURY_FLUX( FLUX, LHGSNOW )
!
!******************************************************************************
!  Subroutine SNOWPACK_MERCURY_FLUX calculates emission of Hg(0) from snow and
!  ice. Emissions are a linear function of Hg mass stored in the snowpack. The
!  Hg lifetime in snow is assumed to be 180 d when T< 270K and 7 d when T>270K
!
!     E = k * SNOW_HG     : k = 6D-8 if T<270K, 1.6D-6 otherwise
!
!  These time constants reflect the time scales of emission observed in the 
!  Arctic and in field studies. Holmes et al 2010
!
!  Arguments as Output
!  ============================================================================
!  (1 ) FLUX (REAL*8) : Flux of Hg(0) [kg/s]
!
!  NOTES:
!******************************************************************************
!
      ! References to F90 modules
      USE TRACERID_MOD,      ONLY : N_Hg_CATS
      USE TIME_MOD,          ONLY : GET_TS_EMIS
      USE DAO_MOD,           ONLY : T, SUNCOS
!      USE OCEAN_MERCURY_MOD, ONLY : SNOW_HG
      USE DEPO_MERCURY_MOD, ONLY : SNOW_HG

#     include "CMN_SIZE"      ! Size parameters

      ! Arguments 
      REAL*8,  INTENT(OUT)  :: FLUX(IIPAR,JJPAR,N_Hg_CATS)
      LOGICAL, INTENT(IN)   :: LHGSNOW

      ! Local variables
      INTEGER               :: I, J, NN, JLOOP
      LOGICAL, SAVE         :: FIRST
      REAL*8                :: DTSRCE, SNOW_HG_NEW, K_EMIT

      !=================================================================
      ! SNOWPACK_MERCURY_FLUX begins here!
      !=================================================================

      ! Initialize
      FLUX = 0D0

      ! Return to calling program if snowpack model is disabled
      IF (.NOT. LHGSNOW) RETURN

      ! Emission timestep [s]
      DTSRCE = GET_TS_EMIS() * 60d0      

      ! Emit Hg(0) at a steady rate, based on 180 d residence
      ! time in snowpack, based on cycle observed at Alert 
      ! (e.g. Steffen et al. 2008)
      K_EMIT = 6D-8

!$OMP PARALLEL DO 
!$OMP+DEFAULT( SHARED ) 
!$OMP+PRIVATE( I, J, NN )
!$OMP+PRIVATE( SNOW_HG_NEW, JLOOP, K_EMIT )
      DO J  = 1, JJPAR
      DO I  = 1, IIPAR

         ! 1-D grid box index for SUNCOS
         JLOOP = ( (J-1) * IIPAR ) + I

         ! If the sun is set, then no emissions, go to next box
         IF (SUNCOS(JLOOP)<0D0) CYCLE 

         ! Decrease residence time to 1 week when T > -3C
         IF (T(I,J,1) > 270D0) THEN
            K_EMIT = 1.6D-6
         ELSE
            K_EMIT = 6D-8
         ENDIF

         DO NN = 1, N_Hg_CATS

            ! Check if there is Hg that could be emitted
            IF (SNOW_HG(I,J,NN)>0D0) THEN

               ! New mass of snow in Hg
               SNOW_HG_NEW = SNOW_HG(I,J,NN) * EXP( - K_EMIT * DTSRCE )

               FLUX(I,J,NN) = MAX( SNOW_HG(I,J,NN) - SNOW_HG_NEW, 0D0 )

               ! Convert mass -> flux
               FLUX(I,J,NN) = FLUX(I,J,NN) / DTSRCE

               SNOW_HG(I,J,NN) = SNOW_HG_NEW

            ENDIF

         ENDDO

      ENDDO
      ENDDO
!$OMP END PARALLEL DO
      
      

      ! Return to calling program
      END SUBROUTINE SNOWPACK_MERCURY_FLUX

!------------------------------------------------------------------------------

      SUBROUTINE INIT_LAND_MERCURY( )

!******************************************************************************
!  Subroutine INIT_LAND_MERCURY allocates and zeroes all module arrays.
!  (ccc, 9/14/09)
!  
!  NOTES:
!******************************************************************************
!
      ! References to F90 modules
      USE ERROR_MOD,    ONLY : ALLOC_ERR
      USE TRACERID_MOD, ONLY : N_Hg_CATS

#     include "CMN_SIZE"     ! Size parameters

      ! Local variables
      INTEGER                      :: AS
      LOGICAL, SAVE         :: IS_INIT = .FALSE. 

      !=================================================================
      ! INIT_MERCURY begins here!
      !=================================================================

      !=================================================================
      ! Allocate arrays
      !=================================================================
      ALLOCATE( TRANSP( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'TRANSP' )
      TRANSP = 0d0

      ! Return to calling program
      END SUBROUTINE INIT_LAND_MERCURY

!-----------------------------------------------------------------------------

      SUBROUTINE CLEANUP_LAND_MERCURY

!******************************************************************************
!  Subroutine CLEANUP_LAND_MERCURY deallocates all module arrays.
!  (ccc, 9/14/09)
!  
!  NOTES:
!******************************************************************************
!
      IF ( ALLOCATED( TRANSP      ) ) DEALLOCATE( TRANSP      )

      ! Return to calling program
      END SUBROUTINE CLEANUP_LAND_MERCURY

!-----------------------------------------------------------------------------

      END MODULE LAND_MERCURY_MOD
