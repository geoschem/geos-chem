!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: land_mercury_mod
!
! !DESCRIPTION: Module LAND\_MERCURY\_MOD contains variables and routines for
! the land emissions for the GEOS-Chem mercury simulation.
!\\
!\\
! !INTERFACE:
!
      MODULE LAND_MERCURY_MOD
!
! !USES:
!
      IMPLICIT NONE
      PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
      PUBLIC :: BIOMASSHG
      PUBLIC :: VEGEMIS
      PUBLIC :: SOILEMIS
      PUBLIC :: LAND_MERCURY_FLUX
      PUBLIC :: GTMM_DR
      PUBLIC :: SNOWPACK_MERCURY_FLUX
      PUBLIC :: INIT_LAND_MERCURY
      PUBLIC :: CLEANUP_LAND_MERCURY
!
! !REVISION HISTORY:
!  02 Jun 2010 - N. E. Selin, C. Carouge - Initial version
!  02 Jun 2010 - C. Carouge  - Group all land emissions routine for mercury 
!                              into this new module.
!  13 Aug 2010 - R. Yantosca - Added modifications for MERRA
!  25 Aug 2010 - R. Yantosca - Treat MERRA in same way as GEOS-5
!  30 Aug 2010 - R. Yantosca - Added more ProTeX headers
!  12 Apr 2011 - J. Fisher   - Add missing code from Holmes 2010
!  08 Feb 2012 - R. Yantosca - Treat GEOS-5.7.x in the same way as MERRA
!  01 Mar 2012 - R. Yantosca - Now reference new grid_mod.F90
!  11 Apr 2012 - R. Yantosca - Now reference new modis_lai_mod.F90
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !PRIVATE TYPES:
!
      ! Plant transpiration rate [m/s]
      REAL*8,  ALLOCATABLE :: TRANSP(:,:)

      CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: land_mercury_flux
!
! !DESCRIPTION: Subroutine LAND\_MERCURY\_FLUX calculates emissions of Hg(0) 
!  from prompt recycling of previously deposited mercury to land, in [kg/s].  
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE LAND_MERCURY_FLUX( LFLUX, LHGSNOW, State_Met )
!
! !USES:
!
      USE DAO_MOD,            ONLY : IS_ICE, IS_LAND
      USE DEPO_MERCURY_MOD,   ONLY : WD_HGP, WD_HG2, DD_HGP, DD_HG2
      USE GIGC_State_Met_Mod, ONLY : MetState
      USE LOGICAL_MOD,        ONLY : LSPLIT
      USE TIME_MOD,           ONLY : GET_TS_EMIS
      USE TRACERID_MOD,       ONLY : ID_Hg0,          N_Hg_CATS

      USE CMN_SIZE_MOD             ! Size parameters
!
! !INPUT PARAMETERS:
!
      LOGICAL,        INTENT(IN)  :: LHGSNOW     ! Use Hg0 from snow?
      TYPE(MetState), INTENT(IN)  :: State_Met   ! Meteorology State object
!
! !OUTPUT PARAMETERS:
!
      REAL*8,         INTENT(OUT) :: LFLUX(IIPAR,JJPAR,N_Hg_CATS) ! Hg0 flux
                                                                  !  [kg/s]
!
! !REVISION HISTORY:
!  30 Aug 2010 - N. E. Selin, C. Holmes, B. Corbitt - Initial version
!  (1 ) Now uses SNOWMAS from DAO_MOD for compatibility with GEOS-5.
!       (eds 7/30/08)
!  (2 ) Now includes REEMFRAC in parallelization; previous versions may have
!       overwritten variable. (cdh, eds 7/30/08)
!  (3 ) Now also reemit Hg(0) from ice surfaces, including sea ice 
!       (cdh, 8/19/08)
!  13 Aug 2010 - R. Yantosca - Add modifications for MERRA
!  25 Aug 2010 - R. Yantosca - Treat MERRA in same way as GEOS-5
!  26 Apr 2011 - J. Fisher   - Use MERRA land fraction information
!  12 Apr 2011 - J. Fisher   - Add missing code from Holmes 2010
!  08 Feb 2012 - R. Yantosca - Treat GEOS-5.7.x in the same way as MERRA
!  09 Nov 2012 - M. Payer    - Replaced all met field arrays with State_Met
!                              derived type object
!  26 Sep 2013 - R. Yantosca - Renamed GEOS_57 Cpp switch to GEOS_FP
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!  
      REAL*8                :: DTSRCE, REEMFRAC
      REAL*8                :: FRAC_SNOW_OR_ICE, FRAC_SNOWFREE_LAND 
      LOGICAL               :: IS_SNOW_OR_ICE, IS_LAND_OR_ICE 
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
!!$OMP+PRIVATE( REEMFRAC, SNOW_HT )  
!$OMP+PRIVATE( REEMFRAC,       FRAC_SNOW_OR_ICE, FRAC_SNOWFREE_LAND )  
!$OMP+PRIVATE( IS_SNOW_OR_ICE, IS_LAND_OR_ICE                       )
      DO J  = 1, JJPAR
      DO I  = 1, IIPAR
      DO NN = 1, N_Hg_CATS
    

         ! For versions prior to MERRA, have to set snow/land fractions
         ! to 1 depending which type of box we are in.
         ! Assume here default of land box with no snow. We will reset
         ! fractions for MERRA here and for prior versions snow/ice boxes
         ! below (jaf, 4/26/11)
         FRAC_SNOW_OR_ICE = 0d0
         FRAC_SNOWFREE_LAND  = 1d0

#if   defined( MERRA ) || defined( GEOS_FP )
         ! Distinguish between ice/snow and snow-free land
         FRAC_SNOW_OR_ICE    = MIN( State_Met%FRSNO(I,J)     +
     &                              State_Met%FRSEAICE(I,J)  + 
     &                              State_Met%FRLANDIC(I,J), 1d0)
         FRAC_SNOWFREE_LAND  = MAX( State_Met%FRLAND(I,J)    - 
     &                              State_Met%FRSNO(I,J),    0d0)

         ! Set this to FALSE in MERRA to avoid overwriting fractions
         IS_SNOW_OR_ICE      = .FALSE.
         IS_LAND_OR_ICE      = (( FRAC_SNOWFREE_LAND > 0d0 ) .OR.
     &                          ( FRAC_SNOW_OR_ICE   > 0d0 ))

#elif defined( GEOS_5 )
         ! GEOS5 snow height (water equivalent) in mm. (Docs incorrectly say m)
         IS_SNOW_OR_ICE = (( IS_ICE ( I, J, State_Met     ))  .OR. 
     &                     ( IS_LAND( I, J, State_Met     )   .AND.
     &                       State_Met%SNOMAS(I,J) > 10d0 ))
         IS_LAND_OR_ICE = (( IS_LAND( I, J, State_Met     ))  .OR.
     &                     ( IS_ICE ( I, J, State_Met     )))

#else
         ! GEOS1-4 snow heigt (water equivalent) in mm
         IS_SNOW_OR_ICE = (( IS_ICE ( I, J, State_Met     ))  .OR.
     &                     ( IS_LAND( I, J, State_Met     )   .AND. 
     &                       State_Met%SNOW(I,J)   > 10d0 ))
         IS_LAND_OR_ICE = (( IS_LAND( I, J, State_Met     ))  .OR.
     &                     ( IS_ICE ( I, J, State_Met     )))
#endif

         ! Reset land fractions if we are in a version prior to MERRA
         ! and also in an ice- or snow-covered box (jaf, 4/26/11)
         IF ( IS_SNOW_OR_ICE ) THEN
            FRAC_SNOW_OR_ICE = 1d0
            FRAC_SNOWFREE_LAND  = 0d0
         ENDIF

         ! If snow or ice on the ground, reemission fraction is 0.6,
         ! otherwise 0.2
         IF ( IS_LAND_OR_ICE ) THEN

            IF (LHGSNOW) THEN
               REEMFRAC = 0.2d0 * FRAC_SNOWFREE_LAND
            ELSE
               REEMFRAC = 0.2d0 * FRAC_SNOWFREE_LAND +
     &                    0.6d0 * FRAC_SNOW_OR_ICE
            ENDIF
            
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
     
      END SUBROUTINE LAND_MERCURY_FLUX
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: biomasshg
!
! !DESCRIPTION: Subroutine BIOMASSHG is the subroutine for Hg(0) emissions 
!  from biomass burning. These emissions are active only for present day 
!  simulations and not for preindustrial simulations.
!\\
!\\
! !INTERFACE:
! 
      SUBROUTINE BIOMASSHG( EHg0_bb )
!
! !USES:
!     
      USE BIOMASS_MOD,  ONLY : BIOMASS
      USE TRACERID_MOD, ONLY : IDBCO
      USE LOGICAL_MOD,  ONLY : LBIOMASS, LPREINDHG
      USE TIME_MOD,     ONLY : GET_TS_EMIS
      USE GRID_MOD,     ONLY : GET_AREA_CM2

      USE CMN_SIZE_MOD     ! Size parameters
      USE CMN_DIAG_MOD     ! Diagnostic arrays & switches
!
! !OUTPUT PARAMETERS:
!
      REAL*8, DIMENSION(:,:),INTENT(OUT) :: EHg0_bb
!
! !REMARKS:
!  Emissions are based on an inventory of CO emissions from biomass burning 
!  (Duncan et al. J Geophys Res 2003), multiplied by a Hg/CO ratio in BB 
!  plumes from Franz Slemr (Poster, EGU 2006).
!                                                                             .
!  Slemr surveyed emission factors from measurements worldwide. Although his
!  best estimate was 1.5e-7 mol Hg/ mol CO, we chose the highest value
!  (2.1e-7 mol Hg/ mol CO) in the range because the simulations shown in
!  Selin et al. (GBC 2008) required large Hg(0) emissions to sustain
!  reasonable atmospheric Hg(0) concentrations. (eck, 11/13/2008)
!
! !REVISION HISTORY:
!  30 Jul 2008 - N. E. Selin, C. Holmes, B. Corbitt - Initial version 
!  12 Apr 2011 - J. Fisher   - Add missing code from Holmes 2010
!  01 Mar 2012 - R. Yantosca - Now use GET_AREA_CM2(I,J,L) from grid_mod.F90
!EOP
!------------------------------------------------------------------------------
!BOC
! 
! !LOCAL VARIABLES:
!     
      REAL*8                 :: DTSRCE, E_CO, AREA_CM2
      INTEGER                :: I, J

      ! Hg molar mass, kg Hg/ mole Hg
      REAL*8,  PARAMETER   :: FMOL_HG     = 200.59d-3

      ! Hg/CO molar ratio in BB emissions, mol/mol
      ! emission factor 1.5e-7 molHg/molCO (Slemr et al poster EGU 2006)
      ! change emission factor to 2.1
      ! Update to Holmes et al. 2010 version (jaf, 4/11/11)
      !REAL*8,  PARAMETER   :: BBRatio_Hg_CO = 2.1D-7
      ! Literature range is 70-240e-9. ARCTAS gives 80-130e-9 
      ! Choose 100e-9 (Holmes et al. 2010 )
      REAL*8,  PARAMETER   :: BBRatio_Hg_CO = 1d-7
 
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
            AREA_CM2 = GET_AREA_CM2( I, J, 1 )

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
!EOC
!-----------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: vegemis
!
! !DESCRIPTION: Subroutine VEGEMIS is the subroutine for Hg(0) emissions from 
!  vegetation by evapotranspiration.
!\\
!\\
! !INTERFACE:
!
      ! Bug fix: VEGEMIS shouldn't be tied to GCAP emissions
      ! (jaf, eds, 4/1/11)
      SUBROUTINE VEGEMIS( LVEGEMIS, EHg0_dist, EHg0_vg, State_Met )
!
! !USES:
!
      USE DAO_MOD,            ONLY : IS_LAND
      USE GIGC_State_Met_Mod, ONLY : MetState
      USE GRID_MOD,           ONLY : GET_AREA_M2
      USE TIME_MOD,           ONLY : GET_MONTH, ITS_A_NEW_MONTH
      USE TIME_MOD,           ONLY : GET_TS_EMIS

      USE CMN_SIZE_MOD             ! Size parameters
!
! !INPUT PARAMETERS:
!
      LOGICAL,                INTENT(IN)  :: LVEGEMIS
      REAL*8, DIMENSION(:,:), INTENT(IN)  :: EHg0_dist
      TYPE(MetState),         INTENT(IN)  :: State_Met   ! Met State object
!
! !OUTPUT PARAMETERS:
!
      REAL*8, DIMENSION(:,:), INTENT(OUT) :: EHg0_vg
!
!
! !REMARKS:
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

! !REVISION HISTORY:
!  30 Aug 2010 - N. Eckley, C. Holmes, B. Corbitt - Initial version
!  01 Mar 2012 - R. Yantosca - Now use GET_AREA_M2(I,J,L) from grid_mod.F90
!EOP
!------------------------------------------------------------------------------
!BOC
! 
! !LOCAL VARIABLES:
!     
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
      ! Bug fix: VEGEMIS shouldn't be tied to GCAP emissions
      ! (jaf, eds, 4/1/11)
      !IF (LGCAPEMIS) THEN
      IF (.NOT. LVEGEMIS) THEN

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
        
            IF (IS_LAND( I, J, State_Met )) THEN

               ! Dry soil Hg concentration, ng Hg /g soil
               DRYSOIL_HG = DRYSOIL_PREIND_HG * EHg0_dist(I,J)

               ! Hg concentration in soil water, ng /m3
               SOILWATER_HG =  DRYSOIL_HG / Kd

               ! Emission from vegetation, ng /m2
               VEG_EMIS = SOILWATER_HG * TRANSP(I,J)

               ! convert from ng /m2 /s -> kg/gridbox/s
               ! Grid box surface area [m2]
               AREA_M2      = GET_AREA_M2( I, J, 1 )
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
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: soilemis
!
! !DESCRIPTION: Subroutine SOILEMIS is the subroutine for Hg(0) emissions 
!  from soils.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE SOILEMIS( EHg0_dist, EHg0_so, State_Met )
!
! !USES:
!
      USE DAO_MOD,            ONLY : IS_LAND
      USE GIGC_State_Met_Mod, ONLY : MetState
      USE GRID_MOD,           ONLY : GET_AREA_M2
      USE MODIS_LAI_MOD,      ONLY : ISOLAI  => GC_LAI
      USE TIME_MOD,           ONLY : GET_MONTH, ITS_A_NEW_MONTH
      USE TIME_MOD,           ONLY : GET_TS_EMIS

      USE CMN_SIZE_MOD             ! Size parameters
!
! !INPUT PARAMETERS:
!
      REAL*8, DIMENSION(:,:), INTENT(IN)  :: EHg0_dist
      TYPE(MetState),         INTENT(IN)  :: State_Met   ! Met State object
!
! !OUTPUT PARAMETERS:
!
      REAL*8, DIMENSION(:,:), INTENT(OUT) :: EHg0_so
!
! !REMARKS:
!  Soil emissions are a function of solar radiation at ground level 
!  (accounting for attenuation by leaf canopy) and surface temperature. 
!  The radiation dependence from Zhang et al. (2000) is multiplied by the 
!  temperature dependence from Poissant and Casimir (1998). 
!  Finally, this emission factor is multiplied by the soil mercury
!  concentration and scaled to meet the global emission total.
!
!  Comments on soil Hg concentration:
!  ----------------------------------
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
! !REVISION HISTORY:
!  30 Aug 2010 - N. Eckley, B. Corbitt - Initial version
!  (1 ) Added comments. (cdh, eds, 7/30/08)
!  (2 ) Now include light attenuation by the canopy after sunset. Emissions
!       change by < 1% in high-emission areas  (cdh, 8/13/2008)
!  (3 ) Removed FRCLND for consistency with other Hg emissions (cdh, 8/19/08)
!  2 June 2010 - C. Carouge  - Solve  
!  13 Aug 2010 - R. Yantosca - Added modifications for MERRA
!  25 Aug 2010 - R. Yantosca - Treat MERRA in same way as GEOS-5
!  26 Apr 2011 - J. Fisher   - Use MERRA land fraction information
!  12 Apr 2011 - J. Fisher   - Bug fixes, add missing code from Holmes 2010
!  08 Feb 2012 - R. Yantosca - Treat GEOS-5.7.x in the same way as MERRA
!  10 Feb 2012 - R. Yantosca - Extend #if statement for SOIL_EMIS_FAC in order
!                              to get the code to compile w/o error.
!  01 Mar 2012 - R. Yantosca - Now use GET_AREA_CM2(I,J,L) from grid_mod.F90
!  11 Apr 2012 - R. Yantosca - Replace lai_mod.F with modis_lai_mod.F90
!  09 Nov 2012 - M. Payer    - Replaced all met field arrays with State_Met
!                              derived type object
!  28 Nov 2012 - R. Yantosca - Replace SUNCOS with State_Met%SUNCOS
!  26 Sep 2013 - R. Yantosca - Renamed GEOS_57 Cpp switch to GEOS_FP
!EOP
!------------------------------------------------------------------------------
!BOC
! 
! !LOCAL VARIABLES:
!     
      REAL*8             :: SOIL_EMIS, DIMLIGHT, TAUZ, LIGHTFRAC
      REAL*8             :: AREA_M2, DRYSOIL_HG
      INTEGER            :: I, J, JLOOP
      REAL*8             :: SUNCOSVALUE
      REAL*8             :: FRAC_SNOWFREE_LAND 
      LOGICAL            :: IS_SNOWFREE_LAND 
!
! !DEFINED PARAMETERS:
!
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
#if defined (GRID4x5)
      REAL*8, PARAMETER  :: SOIL_EMIS_FAC = 1.6D-2 ! for sunlight function, yzh
#elif defined(GRID2x25)
      REAL*8, PARAMETER :: SOIL_EMIS_FAC=1.6D-2
#elif defined(GRID05x0666)
      !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      ! lhu, 1/5/2012, add pre-defined SOIL_EMIS_FAC for 05x0667 simulation
      ! This is a non-physical value and would need to be changed for an
      ! actual mercury simulation
      REAL*8, PARAMETER :: SOIL_EMIS_FAC = 1.6D-2  ! yzh
      !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#elif defined(GRID025x03125)
      !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      ! lhu, 1/5/2012, add pre-defined SOIL_EMIS_FAC for 05x0667 simulation
      ! This is a non-physical value and would need to be changed for an
      ! actual mercury simulation
      REAL*8, PARAMETER :: SOIL_EMIS_FAC = -9.0D99
      !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#else  ! This sets a default value for GCM & ESMF COUPLING
      REAL*8, PARAMETER ::SOIL_EMIS_FAC = 2.4D-2*.5742d0 ! for sunlight function
#endif

      !=================================================================
      ! SOILEMIS begins here!
      !=================================================================

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!!$OMP+PRIVATE( I, J, SOIL_EMIS, JLOOP,      SNOW_HT              ) 
!$OMP+PRIVATE( I,          J,    SOIL_EMIS, JLOOP                ) 
!$OMP+PRIVATE( DRYSOIL_HG, TAUZ, LIGHTFRAC, AREA_M2, SUNCOSVALUE ) 
!$OMP+PRIVATE( IS_SNOWFREE_LAND, FRAC_SNOWFREE_LAND              )
      DO J=1, JJPAR
      DO I=1, IIPAR
         
         ! Except in MERRA, assume whole box is land if conditions met
         ! (jaf, 4/26/11)
         FRAC_SNOWFREE_LAND = 1d0

#if   defined( MERRA ) || defined( GEOS_FP )
         FRAC_SNOWFREE_LAND = MAX( State_Met%FRLAND(I,J) - 
     &                             State_Met%FRSNO(I,J), 0d0 )
         IS_SNOWFREE_LAND   = ( FRAC_SNOWFREE_LAND > 0d0 )

#elif defined( GEOS_5 )
         ! GEOS5 snow height (water equivalent) in mm. (Docs wrongly say m)
         IS_SNOWFREE_LAND = ( IS_LAND( I, J, State_Met ) .AND.
     &                      ( State_Met%SNOMAS(I,J) < 10d0 ))

#else
         ! GEOS1-4 snow heigt (water equivalent) in mm
         IS_SNOWFREE_LAND = ( IS_LAND( I, J, State_Met ) .AND.
     &                      ( State_Met%SNOW(I,J)   < 10d0 ))
#endif

         IF ( IS_SNOWFREE_LAND ) THEN

            ! 1-D grid box index for SUNCOS
            JLOOP = ( (J-1) * IIPAR ) + I
         
            ! attenuate solar radiation based on function of leaf area index
            ! Jacob and Wofsy 1990 equations 8 & 9
            TAUZ = ISOLAI(I,J) * 0.5D0

            ! For very low and below-horizon solar zenith angles, use
            ! same attenuation as for SZA=85 degrees
            SUNCOSVALUE = MAX( State_Met%SUNCOS(I,J), 0.09D0 )

            ! fraction of light reaching the surface is
            ! attenuated based on LAI
            LIGHTFRAC = EXP( -TAUZ / SUNCOSVALUE )

            ! Dry soil Hg concentration, ng Hg /g soil
            DRYSOIL_HG = DRYSOIL_PREIND_HG * EHg0_dist(I,J)           

            ! Soil emissions, ng /m2 /h
            SOIL_EMIS = 
     &           EXP( 0.0011 * State_Met%RADSWG(I,J) * LIGHTFRAC ) *
     &           DRYSOIL_HG * SOIL_EMIS_FAC
     
            ! Grid box surface area [m2]
            AREA_M2   = GET_AREA_M2( I, J, 1 )
 
            ! convert soilnat from ng /m2 /h -> kg /gridbox /s
            EHg0_so(I,J) = SOIL_EMIS * AREA_M2 * 1D-12 / ( 60D0 * 60D0 )

            ! Multiply by fractional land area
            EHg0_so(I,J) = EHg0_so(I,J) * FRAC_SNOWFREE_LAND

         ELSE

            ! no soil emissions from water and ice
            EHg0_so(I,J) = 0D0
        
         ENDIF
        
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      END SUBROUTINE SOILEMIS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: read_nasa_transp
!
! !DESCRIPTION: Subroutine READ\_NASA\_TRANSP reads monthly average 
!  transpiration [m/s] from NASA:
!  for input into the vegetation emissions.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE READ_NASA_TRANSP
!
! !USES:
!
      USE BPCH2_MOD,      ONLY : GET_TAU0,   READ_BPCH2
      USE BPCH2_MOD,      ONLY : GET_RES_EXT
      USE DIRECTORY_MOD,  ONLY : DATA_DIR
      USE TIME_MOD,       ONLY : GET_MONTH,  ITS_A_NEW_MONTH
      USE TRANSFER_MOD,   ONLY : TRANSFER_2D

      USE CMN_SIZE_MOD       ! Size parameters
!
! !REMARKS:
!  Data source:
!   http://gcmd.nasa.gov/records/GCMD_MINTZ_WALKER_SOIL_AND_EVAPO.html
!
!  References:
!  -----------
!  Mintz, Y and G.K. Walker (1993). "Global fields of soil moisture
!   and land surface evapotranspiration derived from observed
!   precipitation and surface air temperature." J. Appl. Meteorol. 32 (8), 
!  1305-1334. 
! 
! !REVISION HISTORY: 
!  15 Sep 2006 - N. E. Selin - Initial version
!  30 Aug 2010 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER             :: I, J, L, MONTH, N
      REAL*4              :: ARRAY(IIPAR,JJPAR,1)
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
     
      FILENAME=TRIM( DATA_DIR ) // 'mercury_201007/nasatransp_'
     &        // GET_RES_EXT() // '.'// CMONTH(MONTH) // '.bpch'

      XTAU     = GET_TAU0(MONTH, 1, 1995 )

      ! Echo info
      WRITE( 6, 100 ) TRIM( FILENAME )
 100  FORMAT( '     - TRANSP_NASA: Reading ', a )     
 
      CALL READ_BPCH2( FILENAME, 'TRANSP-$', 1, 
     &                 XTAU,      IIPAR,     JJPAR,    
     &                 1,         ARRAY,   QUIET=.TRUE. )

      CALL TRANSFER_2D( ARRAY(:,:,1), TRANSP )
      
      ! convert from mm/month to m/s
      TRANSP = TRANSP * 1D-3 * 12D0 / ( 365D0 * 24D0 * 60D0 * 60D0 )
      
      END SUBROUTINE READ_NASA_TRANSP
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: snowpack_mercury_flux
!
! !DESCRIPTION: Subroutine SNOWPACK\_MERCURY\_FLUX calculates emission of 
!  Hg(0) from snow and ice. 
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE SNOWPACK_MERCURY_FLUX( FLUX, LHGSNOW, State_Met )
!
! !USES:
!
      USE DEPO_MERCURY_MOD,   ONLY : SNOW_HG
      USE GIGC_State_Met_Mod, ONLY : MetState
      USE TIME_MOD,           ONLY : GET_TS_EMIS
      USE TRACERID_MOD,       ONLY : N_Hg_CATS

      USE CMN_SIZE_MOD             ! Size parameters
!
! !INPUT PARAMETERS: 
!
      LOGICAL,        INTENT(IN)  :: LHGSNOW     ! Use Hg from snow?
      TYPE(MetState), INTENT(IN)  :: State_Met   ! Meteorology State object
!
! !OUTPUT PARAMETERS:
!
      REAL*8,         INTENT(OUT) :: FLUX(IIPAR,JJPAR,N_Hg_CATS)  ! Hg0 flux
                                                                  !  [kg/s]
!
! !REMARKS:
!  Emissions are a linear function of Hg mass stored in the snowpack. The
!  Hg lifetime in snow is assumed to be 180 d when T< 270K and 7 d when T>270K
!
!     E = k * SNOW_HG     : k = 6D-8 if T<270K, 1.6D-6 otherwise
!
!  These time constants reflect the time scales of emission observed in the 
!  Arctic and in field studies. Holmes et al 2010
! 
! !REVISION HISTORY: 
!  15 Sep 2009 - C. Holmes, S. Carouge - Initial version
!  30 Aug 2010 - R. Yantosca - Added ProTex headers
!  12 Apr 2011 - J. Fisher   - Add missing code from Holmes 2010
!  09 Nov 2012 - M. Payer    - Replaced all met field arrays with State_Met
!                              derived type object
!  28 Nov 2012 - R. Yantosca - Replace SUNCOS with State_Met%SUNCOS
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER       :: I, J, NN, JLOOP
      LOGICAL, SAVE :: FIRST
      REAL*8        :: DTSRCE, SNOW_HG_NEW, K_EMIT

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
         IF ( State_Met%SUNCOS(I,J) < 0D0 ) CYCLE 

            ! Decrease residence time to 1 week when T > -3C
            ! Update to Holmes et al. 2010 version (jaf, 4/11/11)
            ! Emit Hg(0) at a steady rate, based on 180 d residence
            ! time in snowpack, based on cycle observed at Alert 
            ! (e.g. Steffen et al. 2008)
            IF ( State_Met%T(I,J,1) > 270D0) THEN
               ! NOW USE 3 weeks for broader summer peak
               !K_EMIT = 1.6D-6
               K_EMIT = 5D-7
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
      
      END SUBROUTINE SNOWPACK_MERCURY_FLUX
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: gtmm_dr
!
! !DESCRIPTION: GTMM\_DR is a driver to call GTMM from GEOS-Chem.
!\\
!\\
! !INTERFACE: 
!
      SUBROUTINE GTMM_DR( Hg0gtm, State_Met )
! 
! !USES:
!
      USE BPCH2_MOD
      USE DAO_MOD,            ONLY : IS_LAND
      USE DEPO_MERCURY_MOD,   ONLY : CHECK_DIMENSIONS
      USE DEPO_MERCURY_MOD,   ONLY : WD_Hg2, WD_HgP, DD_HgP, DD_Hg2
      USE DEPO_MERCURY_MOD,   ONLY : READ_GTMM_RESTART
      USE DIRECTORY_MOD,      ONLY : DATA_DIR
      USE FILE_MOD,           ONLY : IOERROR
      USE GIGC_State_Met_Mod, ONLY : MetState
      USE inquireMod,         ONLY : findFreeLun
      USE TIME_MOD,           ONLY : EXPAND_DATE, YMD_EXTRACT
      USE TIME_MOD,           ONLY : GET_NYMD, GET_NHMS
    
      USE CMN_SIZE_MOD             ! Size parameters
!
! !INPUT PARAMETERS: 
!
      ! Meteorology State object
      TYPE(MetState), INTENT(IN)  :: State_Met
!
! !OUTPUT PARAMETERS: 
!
      ! Emission of Hg0 calculated by GTMM for the month [kg/s]
      REAL*8,         INTENT(OUT) :: Hg0gtm(IIPAR, JJPAR) 
! 
! !REVISION HISTORY:
!  15 Sep 2009 - C. Carouge  - Initial version
!  03 Aug 2012 - R. Yantosca - Move calls to findFreeLUN out of DEVEL block
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER   :: YEAR    ! Current year 
      INTEGER   :: MONTH   ! Current month
      INTEGER   :: DAY
      INTEGER   :: NYMD, NHMS

      REAL*8, DIMENSION(IIPAR, JJPAR)   :: TSURF    ! Ground temperature
      REAL*8, DIMENSION(IIPAR, JJPAR)   :: PRECIP   ! Total precipitation
                                                    ! for the month
      REAL*8, DIMENSION(IIPAR, JJPAR)   :: SOLAR_W  ! Solar radiation for 
                                                    ! the month

      REAL*4, DIMENSION(IIPAR, JJPAR)   :: TRACER   ! Temporary array

      ! Monthly average deposition arrays
      REAL*8, DIMENSION(IIPAR, JJPAR)   :: Hg0mth_dry
      REAL*8, DIMENSION(IIPAR, JJPAR)   :: Hg2mth_dry
      REAL*8, DIMENSION(IIPAR, JJPAR)   :: Hg2mth_wet

      INTEGER               :: IOS, I, J, L, IU_FILE

      CHARACTER(LEN=255)    :: FILENAME

      ! For binary punch file, version 2.0
      INTEGER               :: NI,        NJ,      NL
      INTEGER               :: IFIRST,    JFIRST,  LFIRST
      INTEGER               :: HALFPOLAR, CENTER180
      INTEGER               :: NTRACER,   NSKIP
      REAL*4                :: LONRES,    LATRES
      REAL*8                :: ZTAU0,     ZTAU1
      CHARACTER(LEN=20)     :: MODELNAME
      CHARACTER(LEN=40)     :: CATEGORY
      CHARACTER(LEN=40)     :: UNIT     
      CHARACTER(LEN=40)     :: RESERVED


      !=================================================================
      ! GTMM_DR begins here!
      !=================================================================
      ! Initialise arrays
      NYMD      = GET_NYMD()
      NHMS      = GET_NHMS()

      TSURF   = 0d0
      PRECIP  = 0d0
      SOLAR_W = 0d0
      Hg0gtm  = 0d0

      ! Reset deposition arrays.
      Hg0mth_dry = 0d0
      Hg2mth_dry = 0d0
      Hg2mth_wet = 0d0

      CALL YMD_EXTRACT( NYMD, YEAR, MONTH, DAY )

      !=================================================================
      ! Read monthly meteorology fields
      !=================================================================

      ! Find a free file LUN
      IU_FILE  = findFreeLun()

      ! File name
      FILENAME = TRIM( DATA_DIR ) // 'mercury_201007/' //
     &           'mean_metfields/' // GET_NAME_EXT() // 
     &           '/mean_YYYYMM.bpch'

      ! Replace YYYY, MM, DD, HH tokens in FILENAME w/ actual values
      CALL EXPAND_DATE( FILENAME, NYMD, NHMS )
    
      ! Echo some input to the screen
      WRITE( 6, '(a)' ) REPEAT( '=', 79 )
      WRITE( 6, 100   ) 
      WRITE( 6, 110   ) TRIM( FILENAME ), IU_FILE
 100  FORMAT( 'G T M M  H g   M E T   F I L E   I N P U T' )
 110  FORMAT( /, 'GTMM_DR: Reading ', a, ' on unit ', i6  )

      ! Open the binary punch file for input
      CALL OPEN_BPCH2_FOR_READ( IU_FILE, FILENAME )
    
      !-----------------------------------------------------------------
      ! Read concentrations -- store in the TRACER array
      !-----------------------------------------------------------------
      DO 
         READ( IU_FILE, IOSTAT=IOS )                              
     &        MODELNAME, LONRES, LATRES, HALFPOLAR, CENTER180
       
         ! IOS < 0 is end-of-file, so exit
         IF ( IOS < 0 ) EXIT
       
         ! IOS > 0 is a real I/O error -- print error message
         IF ( IOS > 0 ) CALL IOERROR( IOS, IU_FILE, 'rd_gtmm_dr:1' )
       
         READ( IU_FILE, IOSTAT=IOS )                               
     &        CATEGORY, NTRACER,  UNIT, ZTAU0,  ZTAU1,  RESERVED,  
     &        NI,       NJ,       NL,   IFIRST, JFIRST, LFIRST,    
     &        NSKIP
       
         IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_FILE, 'rd_gtmm_dr:2' )
       
         READ( IU_FILE, IOSTAT=IOS )                               
     &        ( ( TRACER(I,J), I=1,NI ), J=1,NJ )
       
         IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_FILE, 'rd_gtmm_dr:3' )
       
         !--------------------------------------------------------------
         ! Assign data from the TRACER array to the arrays.
         !--------------------------------------------------------------
       
         ! Process dry deposition data 
         IF ( CATEGORY(1:8) == 'GMAO-2D' ) THEN 
          
            ! Make sure array dimensions are of global size
            ! (NI=IIPAR; NJ=JJPAR, NL=1), or stop the run
            CALL CHECK_DIMENSIONS( NI, NJ, NL )
          
            ! Save into arrays
            IF ( NTRACER == 54 .OR. NTRACER == 59 ) THEN
             
               !----------
               ! Surface temperature
               !----------
             
               ! Store surface temperature in TSURF array
               TSURF(:,:)   = TRACER(:,JJPAR:1:-1)
             
            ELSE IF ( NTRACER == 26 .OR. NTRACER == 29 ) THEN
             
               !----------
               ! Total precipitation
               !----------
             
               ! Store precipitation in PRECIP array
               PRECIP(:,:)   = TRACER(:,JJPAR:1:-1)
             
            ELSE IF ( NTRACER == 37 .OR. NTRACER == 51 ) THEN

               !----------
               ! Solar radiation
               !----------

               ! Store solar radiation in SOLAR_W array
               SOLAR_W(:,:) = TRACER(:,JJPAR:1:-1)

            ENDIF
         ENDIF
       
      ENDDO
    
      ! Close file
      CLOSE( IU_FILE )      
    
      !=================================================================
      ! Read GTMM restart file to get data from previous month
      !=================================================================
      CALL READ_GTMM_RESTART(NYMD,      NHMS, Hg0mth_dry, Hg2mth_dry, 
     &                       Hg2mth_wet )

      !=================================================================
      ! Call GTMM model
      !=================================================================
#if defined( GTMM_Hg )
      CALL GTMM_coupled(YEAR                    , MONTH  , 
     &                  Hg0mth_dry(:,JJPAR:1:-1), 
     &                  Hg2mth_dry(:,JJPAR:1:-1), 
     &                  Hg2mth_wet(:,JJPAR:1:-1), TSURF  , PRECIP, 
     &                  SOLAR_W                  , Hg0gtm(:,JJPAR:1:-1)
     &                  )
#endif

      ! Use LAND/OCEAN mask on the land emissions
      DO J = 1, JJPAR
      DO I = 1, IIPAR
         IF ( .NOT.(IS_LAND(I, J, State_Met )) ) Hg0gtm(I,J) = 0d0
      ENDDO
      ENDDO

      END SUBROUTINE GTMM_DR      
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_land_mercury
!
! !DESCRIPTION: Subroutine INIT\_LAND\_MERCURY allocates and zeroes all 
!  module arrays.
!\\
!\\
! !INTERFACE: 
!
      SUBROUTINE INIT_LAND_MERCURY
! 
! !USES:
!
      USE ERROR_MOD,    ONLY : ALLOC_ERR
      USE TRACERID_MOD, ONLY : N_Hg_CATS

      USE CMN_SIZE_MOD     ! Size parameters
! 
! !REVISION HISTORY:
!  14 Sep 2009 - C. Carouge  - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER              :: AS
      LOGICAL, SAVE :: IS_INIT = .FALSE. 

      !=================================================================
      ! INIT_MERCURY begins here!
      !=================================================================

      ! Allocate arrays
      ALLOCATE( TRANSP( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'TRANSP' )
      TRANSP = 0d0

      END SUBROUTINE INIT_LAND_MERCURY
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: cleanup_land_mercury
!
! !DESCRIPTION: Subroutine CLEANUP\_LAND\_MERCURY deallocates all module 
!  arrays. 
!\\
!\\
! !INTERFACE: 
!
      SUBROUTINE CLEANUP_LAND_MERCURY
! 
! !REVISION HISTORY:
!  14 Sep 2009 - C. Carouge  - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC

      ! Deallocate arrays
      IF ( ALLOCATED( TRANSP      ) ) DEALLOCATE( TRANSP      )

      END SUBROUTINE CLEANUP_LAND_MERCURY
!EOC
      END MODULE LAND_MERCURY_MOD
