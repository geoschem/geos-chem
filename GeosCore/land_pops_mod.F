!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: land_pops_mod
!
! !DESCRIPTION: Module LAND\_POPS\_MOD contains variables and routines for
! the land emissions for the GEOS-Chem persistent organic pollutant simulation.
!\\
!\\
! !INTERFACE:
!
      MODULE LAND_POPS_MOD
!
! !USES:
!
      IMPLICIT NONE
      PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
      PUBLIC :: VEGEMISPOP
      PUBLIC :: LAKEEMISPOP
      PUBLIC :: SOILEMISPOP
!      PUBLIC :: SNOWEMISPOP
!      PUBLIC :: INIT_LAND_POPS
!      PUBLIC :: CLEANUP_LAND_POPS
!
! !REVISION HISTORY:
!  21 Aug 2012 - C.L. Friedman - Initial version based on LAND_MERCURY_MOD

!EOP
!------------------------------------------------------------------------------
!BOC
!
! !PRIVATE TYPES:
!
      ! Plant transpiration rate [m/s]
!      REAL*8,  ALLOCATABLE :: TRANSP(:,:)
      REAL*8,  PARAMETER   :: SMALLNUM = 1D-20


      CONTAINS
!EOC

!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: soilemispop
!
! !DESCRIPTION: Subroutine SOILEMISPOP is the subroutine for secondary
!  POP emissions from soils.
!\\
!\\
! !INTERFACE:
!
         SUBROUTINE SOILEMISPOP( POP_SURF, F_OC_SOIL, EPOP_SOIL )
!
! !USES:
!
      USE TIME_MOD, ONLY      : GET_TS_EMIS, GET_TS_CHEM
      USE DAO_MOD, ONLY       : IS_ICE, IS_LAND, IS_WATER
      USE DAO_MOD, ONLY       : FRLAKE, SNICE, SNOW
      USE MODIS_LAI_MOD, ONLY : GC_LAI, GC_LAI_CM, GC_LAI_PM !hma, 7Mar2013
!      USE LAI_MOD, ONLY       : ISOLAI, MISOLAI, PMISOLAI
      USE DAO_MOD, ONLY       : TSKIN, AIRVOL, T, TS
      USE DAO_MOD, ONLY       : IS_ICE, IS_LAND, SNOMAS 
      USE GRID_MOD, ONLY      : GET_AREA_CM2, GET_YMID
      USE GRID_MOD, ONLY      : GET_AREA_M2
      USE TRACER_MOD, ONLY    : STT
      USE TRACERID_MOD, ONLY  : IDTPOPG
      USE DIAG53_MOD, ONLY    : AD53 !ND53,  LD53,  AD53
      USE GET_POPSINFO_MOD, ONLY  : GET_POP_DEL_H, GET_POP_KOA
      USE GET_POPSINFO_MOD, ONLY  : GET_POP_XMW

      ! Replace include statements (H. Amos, 7Mar2013)
      USE CMN_SIZE_MOD   ! Size parameters
      USE CMN_DIAG_MOD   ! Diagnostic arrays & switches
!#     include "CMN_SIZE"     ! Size parameters
!#     include "CMN_DIAG"     ! Diagnostic arrays & switches
!
! !INPUT PARAMETERS:
!
!      LOGICAL,                INTENT(IN) :: LSECEMISPOP
      REAL*8, DIMENSION(:,:), INTENT(IN) :: POP_SURF
      REAL*8, DIMENSION(:,:), INTENT(IN) :: F_OC_SOIL
!
! !OUTPUT PARAMETERS:
!
      REAL*8, DIMENSION(:,:), INTENT(OUT):: EPOP_SOIL
!
! !REMARKS:
!
!  
!
! !REVISION HISTORY:
!  
!EOP
!------------------------------------------------------------------------------
!BOC
! 
! !LOCAL VARIABLES:
!     

      INTEGER                :: I, J, L
      REAL*8                 :: POPG, POPG_GL!, SUNCOSVALUE
      REAL*8                 :: FRAC_SNOWFREE_LAND, TK_SURF 
      REAL*8                 :: KSA_T, FLUX, F_OC, K_MT
      REAL*8                 :: SOIL_CONC, DTSRCE, KOA_T
      REAL*8                 :: DIFF, FSOIL, FAIR, DS
      REAL*8                 :: DSA, DAD, DWD, PL
      REAL*8                 :: ZSOIL, ZAIR, TK
      REAL*8                 :: DTCHEM, NEWSOIL, E_KDEG
      REAL*8                 :: FRAC_LEAF, FRAC_LAKE, FRAC_SOIL
      REAL*8                 :: FUG_R
      LOGICAL                :: IS_SNOW_OR_ICE, IS_LAND_OR_ICE
      LOGICAL                :: IS_SNOWFREE_LAND

! !External functions
      REAL*8, EXTERNAL   :: BOXVL ! Grid box vol in cm3
!
! !DEFINED PARAMETERS:
      ! Delta H for POP [kJ/mol]. Delta H is enthalpy of phase transfer
      ! from gas phase to OC. For now we use Delta H for phase transfer 
      ! from the gas phase to the pure liquid state. 
      ! For PHENANTHRENE: 
      ! this is taken as the negative of the Delta H for phase transfer
      ! from the pure liquid state to the gas phase (Schwarzenbach,
      !  Gschwend, Imboden, 2003, pg 200, Table 6.3), or -74000 [J/mol].
      ! For PYRENE:
      ! this is taken as the negative of the Delta H for phase transfer
      ! from the pure liquid state to the gas phase (Schwarzenbach,
      ! Gschwend, Imboden, 2003, pg 200, Table 6.3), or -87000 [J/mol]. 
      ! For BENZO[a]PYRENE:
      ! this is also taken as the negative of the Delta H for phase transfer
      ! from the pure liquid state to the gas phase (Schwarzenbach,
      ! Gschwend, Imboden, 2003, pg 452, Prob 11.1), or -110,000 [J/mol]
      REAL*8     :: DEL_H

      ! R = universal gas constant for adjusting KOA for temp: 8.3145 [J/mol/K OR m3*Pa/K/mol]
      REAL*8, PARAMETER     :: R          = 8.31d0 

      ! Molecular weight
      ! For phe, 0.17823 kg/mol
      ! For pyr, 0.20225 kg/mol
      ! For BaP, 0,25231 kg/mol
      REAL*8     :: MW
!
      ! For PHENANTHRENE:
      ! log KOA_298 = 7.64, or 4.37*10^7 [unitless]
      ! For PYRENE:
      ! log KOA_298 = 8.86, or 7.24*10^8 [unitless]
      ! For BENZO[a]PYRENE:
      ! log KOA_298 = 11.48, or 3.02*10^11 [unitless]
      ! (Ma et al., J. Chem. Eng. Data, 2010, 55:819-825).
      REAL*8     :: KOA_298

      ! Set transfer velocity and diffusion coefficient values
      REAL*8, PARAMETER     :: KSA        = 1d0     ![m/h]
      REAL*8, PARAMETER     :: BA         = 0.04d0  ![m2/h]
      REAL*8, PARAMETER     :: BW         = 4d-6    ![m2/h]

      ! Set soil degradation rate
      REAL*8, PARAMETER     :: DEGR       = 3.5d-5  ![/h]

      REAL*8                :: DUM

      !=================================================================
      ! SOILEMISPOP begins here!
      !=================================================================

      DUM = 1.0
      DEL_H = GET_POP_DEL_H(DUM)
      KOA_298 = GET_POP_KOA(DUM)
      MW = GET_POP_XMW(DUM)

!      IF (.NOT. LSECEMISPOP) THEN

!         EPOP_SOIL = 0D0

!      ELSE

      ! Emission timestep [s]
      DTSRCE  = GET_TS_EMIS() * 60d0

      ! Chemistry timestep [h]
      DTCHEM = GET_TS_CHEM()

      DO J=1, JJPAR
      DO I=1, IIPAR

         ! Set logicals
         ! Is grid box covered by land/ice or by water? (IS_LAND_OR_ICE)
         ! IS_LAND will return non-ocean boxes but may still contain lakes
         ! If land, is it covered by snow/ice? (IS_SNOW_OR_ICE)
         IS_LAND_OR_ICE = ( (IS_LAND(I,J)) .OR. (IS_ICE(I,J)) ) 
         IS_SNOW_OR_ICE = 
     &      ( (IS_ICE(I,J)) .OR. (IS_LAND(I,J) .AND. SNOMAS(I,J)>10d0) )
         
!THIS IF-STATEMENT ISN'T WORKING. NEVER TRUE. DEBUG! (H. Amos, 16 May 2013)
c$$$#if defined ( GEOS_5 )
c$$$         ! Set logicals
c$$$         ! Is grid box covered by land/ice or by water? (IS_LAND_OR_ICE)
c$$$         ! IS_LAND will return non-ocean boxes but may still contain lakes
c$$$         ! If land, is it covered by snow/ice? (IS_SNOW_OR_ICE)
c$$$         IS_LAND_OR_ICE = ( (IS_LAND(I,J)) .OR. (IS_ICE(I,J)) ) 
c$$$         IS_SNOW_OR_ICE = 
c$$$     &      ( (IS_ICE(I,J)) .OR. (IS_LAND(I,J) .AND. SNOMAS(I,J)>10d0) )
c$$$         IF (I .EQ. 10 .AND. J .EQ. 10 ) PRINT*,'ICELAND Success!' ! hma debug
c$$$         
c$$$#elif defined ( GCAP )
c$$$         ! GCAP uses different snow variable 
c$$$         IS_LAND_OR_ICE = ( (IS_LAND(I,J)) .OR. (IS_ICE(I,J)) )    
c$$$         IS_SNOW_OR_ICE = 
c$$$     &      ( (IS_ICE(I,J)) .OR. (IS_LAND(I,J)  .AND. SNOW(I,J)>10d0) )
c$$$#endif

         ! Do soils routine only if we are on land that is not covered with
         ! snow or ice
         IF ((IS_LAND_OR_ICE) .AND. .NOT. ( IS_SNOW_OR_ICE )) THEN
!         IF( IS_LAND(I,J) ) THEN !hma debugging   
!            PRINT*,'Made it inside loop' !hma debug

            ! Get fraction of grid box covered by leaf surface area
            ! Do not consider different vegetation types for now
!            FRAC_LEAF = ISOLAI(I,J)
            ! DON'T DO THIS NOW B/C FRAC_LEAF CAN BE GREATER THAN 1D0
            ! PRESUMABLY CONSIDERING SURFACE AREA IN 3 DIMENSIONS

            ! Get fraction of grid box covered by lake surface area
            FRAC_LAKE = FRLAKE(I,J)

!THIS IF-STATEMENT ISN'T WORKING. NEVER TRUE. DEBUG! (H. Amos, 16 May 2013)
c$$$
c$$$#if defined ( GEOS_5 )
c$$$            ! Get fraction of grid box covered by leaf surface area
c$$$            ! Do not consider different vegetation types for now
c$$$!            FRAC_LEAF = ISOLAI(I,J)
c$$$            ! DON'T DO THIS NOW B/C FRAC_LEAF CAN BE GREATER THAN 1D0
c$$$            ! PRESUMABLY CONSIDERING SURFACE AREA IN 3 DIMENSIONS
c$$$
c$$$            ! Get fraction of grid box covered by lake surface area
c$$$            FRAC_LAKE = FRLAKE(I,J)
c$$$#elif defined ( GCAP )
c$$$            ! Lake surface is not defined for GCAP
c$$$            FRAC_LAKE = 0d0
c$$$#endif
            ! Get fraction of land remaining
            ! Assume the remaining land is soil and get OC content. If remaining land is not
            ! soil (e.g., desert), there should be a characteristically low OC content
            ! that will have little capacity to store POPs
            ! ONLY SUBTRACT FRAC LAKE NOW
            FRAC_SOIL = MAX(1d0 - FRAC_LAKE, 0d0)

            ! Get surface skin temp [K]
            TK_SURF = TSKIN(I,J)

!THIS IF-STATEMENT ISN'T WORKING. NEVER TRUE. DEBUG! (H. Amos, 16 May 2013)
c$$$#if defined ( GEOS_5 )
c$$$            ! Get surface skin temp [K]
c$$$            TK_SURF = TSKIN(I,J)
c$$$#elif defined ( GCAP )
c$$$            TK_SURF = TS(I,J)
c$$$#endif

            ! Get air temp [K]
            TK = T(I,J,1)

            ! Get gas phase air POP concentration at surface in mol/m3
            POPG = MAX( STT(I,J,1,IDTPOPG), SMALLNUM ) ![kg in gridbox]
            ! kg / (0.178 kg/mol) /m3 in gridbox
            POPG = POPG / MW / AIRVOL(I,J,1) ! mol/m3
            !WRITE(6,*) 'POPG (mol/m3) =', POPG

            ! Get soil concentration in top 5 cm of soil (following Howsam et al 2000)
            ! From Howsam et al, soil burdens are equal to 2.6 years deposition for PHE,
            ! 10 years for PYR, and 9.4 years for BaP
            ! Convert to mol/m3      
            ! 2.6 yrs * kg deposited to soil in 1 yr * / 0.178 kg/mol
            ! / area grid box (m2) / 0.05 m 
!Bug fix (H. Amos, 8 Mar 2013)
!            SOIL_CONC = 10d0 * POP_SURF(I,J) / MW / GET_AREA_M2( J ) /
!     &                  5d-2! mol/m3   
            SOIL_CONC = 10d0 * POP_SURF(I,J) / MW / GET_AREA_M2(I,J,1) /
     &                  5d-2! mol/m3   

            !WRITE(6,*) 'SOIL_CONC (mol/m3) =', SOIL_CONC

            ! Get rid of mass due to degradation
            ! Use rate constant for BaP from Mackay and Paterson 1991: 3.5*10^-5 /h
            ! Calculate exponential factor
            E_KDEG = EXP (-DEGR * DTCHEM)
            ! Adjust conc
            NEWSOIL = SOIL_CONC * E_KDEG

            ! Get foc from GTMM saved files (done in pops_mod.f)
            F_OC = F_OC_SOIL(I,J) 

            ! Define temperature-dependent KOA:
            KOA_T = KOA_298 * EXP((-DEL_H/R) * ((1d0/TK_SURF) - 
     &                  (1d0/298d0))) 

            ! Dimensionless coefficient (mol/m3 soil / mol/m3 air) 
            ! KSA = 1.5 (fTOC)*Koa
            KSA_T = 1.5 * F_OC * KOA_T
            KSA_T = MAX( KSA_T, SMALLNUM )

            ! Calculate fugacities from concentrations by dividing by "Z" values, 
            ! or the fugacity capacity in mol/m3*Pa following Mackay and Paterson, 1991

            ! fsoil = Csoil * R * T / KSA [Pa]
            ! where Csoil is in mol/m3, R is in Pa * m3 / mol K
            ! T is in K and KSA is dimensionless
            FSOIL = NEWSOIL * R * TK_SURF / KSA_T

            ! fair = Cair * R * T [Pa]
            ! where Cair is in mol/m3, R is in Pa * m3 / mol K and T is in K
            FAIR = POPG * R * TK

            ! Calculate the fugacity gradient [Pa]
            ! If the gradient is negative, fair is larger and the POP will diffuse from air to soil
            ! If the gradient is positive, fsoil is larger and POP will diffuse from soil to air
            DIFF = FSOIL - FAIR
!            WRITE (6,*) 'DIFF (mol/m3)=', NEWSOIL/KSA_T - POPG
            FUG_R = FSOIL/FAIR

            ! Calculate "Z" values from fugacities. 
            ! Z is the fugacity capacity in mol/m3*Pa. C = Z*f, so Z = C/f
            ZAIR  = POPG / FAIR ! (mol/m3) / (Pa)
            ZSOIL = NEWSOIL / FSOIL ! (mol/m3) / (Pa)

            ! Calculate the "D" value, or the transfer coefficient that describes the movement of 
            ! POP between phases (Mackay and Paterson, 1991). [mol/h*Pa]
            ! The D value for soil-air diffusion is given by
            ! Ds = 1 / (1/Dsa + 1/(Dad + Dwd))
            ! where Dsa is the air-side boundary layer diffusion parameter [mol/h*Pa]
            ! Dad is the diffusion parameter between soil particles and "soil air" [mol/h*Pa]
            ! Dwd is the diffusion parameter between soil particles and porewater [mol/h*Pa]
            ! Dsa is in series with soil-air and soil-water diffusion, which are in parallel
  
            ! Need to define each D value
            ! DSA = kSA * Zair
            ! where kSA is a mass transfer coefficient [m/h],
            ! Zair is the air fugacity capacity [mol/m3*Pa]
            ! ***********
            ! DAD = BA * Zair
            ! where BA is the molecular diffusivity in air [m2/h]
            ! ***********
            ! DWD = BW * Zsoil
            ! where BW is the molecular diffusivity in water [m2/h]
            ! **** PL = the soil diffusion pathlength, set to half the soil depth (0.025 m)
            DSA = KSA * ZAIR  ! (m/h) * (mol/m3*Pa) = mol/m2*h*Pa
            DAD = BA* ZAIR      ! (m2/h) * (mol/m3*Pa) = mol/m*h*Pa
            DWD = BW * ZSOIL    ! (m2/h) * (mol/m3*Pa) = mol/m*h*Pa
            DS = 1d0 / ( 1d0/DSA + PL/(DAD+DWD) ) ! mol/(m2*h*Pa) [* m3*Pa/K/mol = m/h/K

            ! Calculate Flux in mol/m2/h
            FLUX = DS * DIFF

            ! Change to units of ng/m2/d for storage
            FLUX = FLUX * 24d0 * MW * 1d12 

            ! Kludge soil emissions from poles for now
            ! Bug somewhere that allows GCAP versions to think some high polar boxes during 
            ! some months are land rather than ice - results in extremely high fluxes
!Bug fix (H. Amos, 8 Mar 2013)
!            IF  (GET_YMID(J) >  60  .OR. GET_YMID(J) < -60 ) THEN
            IF (GET_YMID(I,J,1) >  60  .OR. GET_YMID(I,J,1) < -60 ) THEN
               FLUX = 0d0
               DIFF = 0d0
            ENDIF

            ! Convert to an emission rate in kg/s for returning to pops_mod.f
!Bug fix (H. Amos, 8 Mar 2013)
!            EPOP_SOIL(I,J) = MAX(FLUX * GET_AREA_M2( J ) / 24d0 / 3600d0 
!     &                         / 1d12, 0d0 )
            EPOP_SOIL(I,J) = MAX(FLUX * GET_AREA_M2(I,J,1) / 24d0  
     &                         / 3600d0 / 1d12, 0d0 )


            ! Multiply the mass emitted by the fraction of land that is soil
            EPOP_SOIL(I,J) = FRAC_SOIL * EPOP_SOIL(I,J)
         
            ! If the flux is positive, then the direction will be from the soil to the air. 
            ! Store this in a diagnostic array. 
            ! If the flux is zero or negative, store it in a separate array.

            IF (ND53 > 0) THEN
               IF (FLUX > 0d0) THEN 
                  ! Store total mass emitted from soil [kg] in ND53 diagnostic.
                  ! (We don't care about the mass in the other direction right now)
                  AD53(I,J,5) = AD53(I,J,5) + (EPOP_SOIL(I,J) * DTSRCE)  
                  ! Store positive flux
                  AD53(I,J,8) = AD53(I,J,8) + FLUX
                  ! Make sure negative flux diagnostic has nothing added to it
                  AD53(I,J,9) = AD53(I,J,9) + 0d0
                  ! Store the soil/air fugacity ratio
                  AD53(I,J,14) = AD53(I,J,14) + FUG_R
             
               ELSE IF (FLUX <= 0d0) THEN
                  ! Store the negative flux
                  AD53(I,J,9) = AD53(I,J,9) + FLUX
                  ! Add nothing to positive flux or mass diagnostics
                  AD53(I,J,5) = AD53(I,J,5) + 0d0
                  AD53(I,J,8) = AD53(I,J,8) + 0d0
                  ! Continue to store the fugacity ratio
                  AD53(I,J,14) = AD53(I,J,14) + FUG_R
               ENDIF
            ENDIF

            ELSE 
            ! We are not on land or the land is covered with ice or snow
               EPOP_SOIL(I,J) = 0d0
               FLUX = 0d0
               AD53(I,J,8) = 0d0 
               AD53(I,J,5) = 0d0
               AD53(I,J,9) = 0d0
               AD53(I,J,14) = 0d0
            ENDIF
        
      ENDDO
      ENDDO


      END SUBROUTINE SOILEMISPOP
!EOC

!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: soilemispop
!
! !DESCRIPTION: Subroutine LAKEEMISPOP is the subroutine for secondary
!  POP emissions from lakes.
!\\
!\\
! !INTERFACE:
!
         SUBROUTINE LAKEEMISPOP( POP_SURF, EPOP_LAKE )
!
! !USES:
!
      USE TIME_MOD, ONLY      : GET_TS_EMIS, GET_TS_CHEM
      USE DAO_MOD, ONLY       : IS_ICE, IS_LAND, IS_WATER
      USE DAO_MOD, ONLY       : FRLAKE
      USE MODIS_LAI_MOD, ONLY : GC_LAI, GC_LAI_CM, GC_LAI_PM ! hma, 7Mar2013
!      USE LAI_MOD, ONLY       : ISOLAI, MISOLAI, PMISOLAI
      USE DAO_MOD, ONLY       : TSKIN, AIRVOL, T
      USE DAO_MOD, ONLY       : IS_ICE, IS_LAND, SNOMAS
      USE DAO_MOD, ONLY       : PSC2
      USE GRID_MOD, ONLY      : GET_AREA_CM2
      USE GRID_MOD, ONLY      : GET_AREA_M2
      USE TRACER_MOD, ONLY    : STT
      USE TRACERID_MOD, ONLY  : IDTPOPG
      USE DIAG53_MOD, ONLY    : AD53 !ND53,  LD53,  AD53
      USE GET_POPSINFO_MOD, ONLY : GET_POP_HSTAR, GET_POP_DEL_Hw
      USE GET_POPSINFO_MOD, ONLY : GET_POP_XMW
      
      ! Replace include statements (H. Amos, 7Mar2013)
      USE CMN_SIZE_MOD   ! Size parameters
      USE CMN_DIAG_MOD   ! Diagnostic arrays & switches
!#     include "CMN_SIZE"     ! Size parameters
!#     include "CMN_DIAG"     ! Diagnostic arrays & switches

!
! !INPUT PARAMETERS:
!
!      LOGICAL,                INTENT(IN) :: LSECEMISPOP
      REAL*8, DIMENSION(:,:), INTENT(IN) :: POP_SURF

!
! !OUTPUT PARAMETERS:
!
      REAL*8, DIMENSION(:,:), INTENT(OUT):: EPOP_LAKE
!
! !REMARKS:
!
! !REVISION HISTORY:
!  
!EOP
!------------------------------------------------------------------------------
!BOC
! 
! !LOCAL VARIABLES:
!     

      INTEGER                :: I, J, L
      REAL*8                 :: POP_AQ_TOT, POP_DISS
      REAL*8                 :: FRAC_SNOWFREE_LAND, TK_SURF 
      REAL*8                 :: KAW_T, FLUX, KOL_T
      REAL*8                 :: DTSRCE
      REAL*8                 :: KA_H2O, KA_POP, KW_CO2, KW_POP
      REAL*8                 :: DA_H2O, DA_POP, DW_CO2, DW_POP
      REAL*8                 :: SCH_CO2, SCH_POP, KOL
      REAL*8                 :: TK, PRESS
      REAL*8                 :: DTCHEM, E_KDEG, C_DISS
      REAL*8                 :: POPG, U10M, ALPHA
      REAL*8                 :: FRAC_LAKE, VISC_H2O
      LOGICAL                :: IS_SNOW_OR_ICE, IS_LAND_OR_ICE
      LOGICAL                :: IS_SNOWFREE_LAND

! !External functions
      REAL*8, EXTERNAL   :: BOXVL ! Grid box vol in cm3
!
! !DEFINED PARAMETERS:

      ! Delta H for POP:   
      ! For PHENANTHRENE: 
      ! this is the Delta H for phase transfer
      ! from air to water (Schwarzenbach,
      !  Gschwend, Imboden, 2003, pg 200, Table 6.3), or 47 [kJ/mol].
      ! For PYRENE: 43000 [kJ/mol]
      REAL*8     :: DEL_HW

      ! R = universal gas constant for adjusting KOA for temp: 8.3145d-3 [kJ/mol/K]
      REAL*8, PARAMETER     :: R          = 8.31d-3 

      ! Molecular weight
      ! For phe, 0.17823 kg/mol
      REAL*8     :: MWPOP

      ! Molecular weight of air
      REAL*8, PARAMETER     :: MWAIR      = 28.97d0 ! g/mol

      ! Molecular weight of water
      REAL*8, PARAMETER     :: MWH2O      = 18.1d0 ! g/mol

      ! Molar volumes calculated following Abraham and McGowan 1987 as summarized 
      ! by Schwarzenbach et al. 2003. 
      ! Each element is assigned a characteristic atomic volume, and an atomic volume
      ! of 6.56 cm3/mol is subtracted for each bond, no matter whether single, double,
      ! or triple 
      ! C = 16.35,  H = 8.71,  O = 12.43, N = 14.39, P = 24.87, F = 10.48
      ! Br = 26.21, I = 34.53, S = 22.91, Si = 26.83
      
      ! Molar volume of water
      ! 2*(8.71) + 12.43 - 2*(6.56) = 16.73
      REAL*8, PARAMETER     :: V_H2O  = 16.73d0 ! cm3/mol

      ! Molar volume of CO2
      ! 16.35 + 2*(12.43) - 2*(6.56) = 28.1
      REAL*8, PARAMETER     :: V_CO2  = 28.1d0  ! cm3/mol
     
      ! Molar volume of POP
      ! For PHE (C16H10):
      ! 16*(16.35) + 10*(8.71) - 29*(6.56) = 
      REAL*8, PARAMETER     :: V_POP = 538.94d0 ! cm3/mol

      ! Molar volume of air - average of gases in air
      REAL*8, PARAMETER     :: V_AIR  = 20.1d0     ! cm3/mol
      
      ! For PHENANTHRENE:
      ! log KAW_298 = -2.76, or 1.74*10-3 [unitless]
      ! For PYRENE:
      ! log KAW_298 = -3.27, or 5.37*10-4 [unitless]
      REAL*8     :: KAW_298 

      ! Set the kinematic viscosity of freshwater at 20C
!      REAL*8, PARAMETER     :: VISC_H2O       ! = 1d0     ![cm2/s]

      ! Set aqueous degradation rate
!      REAL*8, PARAMETER     :: DEGR       != 3.5d-5  ![/h]

      ! External functions
      REAL*8,  EXTERNAL     :: SFCWINDSQR 

      REAL*8              :: DUM

      !=================================================================
      ! LAKEEMISPOP begins here!
      !=================================================================

      DUM = 1.0
      DEL_HW = GET_POP_DEL_HW(DUM)
      MWPOP = GET_POP_XMW(DUM)
      KAW_298 = GET_POP_HSTAR(DUM)


!      IF (.NOT. LSECEMISPOP) THEN

!         EPOP_SOIL = 0D0

!      ELSE

#if defined ( GCAP )
      ! No definition of lakes, consider emission from lakes to be zero
       EPOP_LAKE(I,J) = 0d0
       FLUX = 0d0
       AD53(I,J,10) = 0d0 
       AD53(I,J,6) = 0d0
       AD53(I,J,11) = 0d0
       AD53(I,J,15) = 0d0
              
      ! Return to calling program
      RETURN

#elif defined ( GEOS_5 )
      
      ! Do lake emissions routine:

      ! Emission timestep [s]
      DTSRCE  = GET_TS_EMIS() * 60d0

      ! Chemistry timestep [h]
      DTCHEM = GET_TS_CHEM()

      DO J=1, JJPAR
      DO I=1, IIPAR

         ! Set logicals
         ! Is grid box covered by land/ice or by water? (IS_LAND_OR_ICE)
         ! IS_LAND will return non-ocean boxes but may still contain lakes
         ! If land, is it covered by snow/ice? (IS_SNOW_OR_ICE)
         IS_LAND_OR_ICE = ( (IS_LAND(I,J)) .OR. (IS_ICE(I,J)) ) 
         IS_SNOW_OR_ICE = 
     &      ( (IS_ICE(I,J)) .OR. (IS_LAND(I,J) .AND. SNOMAS(I,J)>10d0) )

         ! Do soils routine only if we are on land that is not covered with
         ! snow or ice
         IF ((IS_LAND_OR_ICE) .AND. .NOT. ( IS_SNOW_OR_ICE )) THEN

            ! Get fraction of grid box covered by lake surface area
            FRAC_LAKE = FRLAKE(I,J)
!            WRITE(6,*) 'FRLAKE =', FRLAKE(I,J)

            IF ( FRAC_LAKE > 0d0 ) THEN

            ! Get surface skin temp [K]
            TK_SURF = TSKIN(I,J)

            ! Get air temp [K]
            TK = T(I,J,1)

            ! Get surface pressure at end of dynamic time step [hPa]
            PRESS = PSC2(I,J)

            ! Convert to units of atm 
            PRESS = PRESS / 1013.25d0

            ! Get gas phase air POP concentration at surface in mol/m3
            POPG = MAX( STT(I,J,1,IDTPOPG), SMALLNUM ) ![kg in gridbox]
            ! kg / (kg/mol) /m3 in gridbox
            POPG = POPG / MWPOP / AIRVOL(I,J,1) ! mol/m3

            ! Get the dissolved POP concentration at lake surface in mol/m3
            ! Distribute the total deposited mass to a volume that best matches
            ! observed dissolved concentrations
            ! Future versions should consider aqueous particle concentrations
            ! and sinking rates and photolytic/microbial degradation

            ! Start with 1 m - scale by 100
!Bug fix (H. Amos, 8 Mar 2013)
!            C_DISS = POP_SURF(I,J) / MWPOP / GET_AREA_M2( J ) / 100d0
            C_DISS = POP_SURF(I,J) / MWPOP / GET_AREA_M2(I,J,1) / 100d0
            !WRITE(6,*) 'CDISS (mol/m3) =', C_DISS

            ! Get the wind speed at 10 m squared
            U10M = SQRT( SFCWINDSQR(I,J) ) ! m/s

            ! Need to calculate water-side and air-side mass transfer coefficients
            ! Start with air-side
            ! Relate air-side MTC of POP to that of H2O
            ! First, calculate air-side MTC for water following Schwarzenbach, 
            ! Gschwend, Imboden 2003
            KA_H2O = 0.2d0 * U10M + 0.3d0   ! cm/s

            ! Relate air-side MTC for water to that of POP via diffusivities
            ! following Schwarzenbach et al 2003
            ! Calculate temperature-dependent diffusivities in air first
            ! folling Fuller et al. 1966 (summarized by Schwarzenbach et al. 2003)
            DA_POP = 1d-3 * TK**1.75d0 * ( (1d0 / MWAIR ) + 
     &               (1d0 / (MWPOP*1d3) ) )**0.5d0 / ( PRESS * 
     &               ( V_AIR**(1d0/3d0) + V_POP**(1d0/3d0))**2d0 )  ! cm2/s

            DA_H2O = 1d-3 * TK**1.75d0 * ( (1d0 / MWAIR ) + 
     &               (1d0 / MWH2O ) )**0.5d0 / ( PRESS * 
     &               ( V_AIR**(1d0/3d0) + V_H2O**(1d0/3d0))**2d0 )  ! cm2/s

            ! Relate POP and H2O air-side MTCs
            KA_POP = KA_H2O * ( DA_POP / DA_H2O )**( 0.67d0 ) ! cm/s
 
            ! Now calculate water-side MTCs
            ! Start with calculating the water side MTC of CO2
            ! This depends on wind speed (Schwarzenbach et al 2003)
            ! Three different scenarios for the water surface under different wind speeds 
            ! are considered:
            ! Smooth Surface Regime (SSR): u10 <= 4.2 m/s
            ! Rough Surface Regime (RSR):  4.2 m/s < u10 <= 13 m/s
            ! Breaking Wave Regime (BWR):  u10 > 13 m/s
            ! Alpha, the exponent in the relationship for the water side MTC,
            ! also depends on the wind speed. Set this as well

            IF (U10M <= 4.2d0) THEN
               KW_CO2 = 0.65d-3       ! cm/s
               ALPHA  = 0.67d0
            ELSE IF (U10M > 4.2d0 .AND. U10M <= 13d0) THEN
               KW_CO2 = ( 0.79d0 * U10M - 2.68d0 ) * 1d-3
               ALPHA  = 0.5d0
            ELSE IF (U10M > 13d0) THEN 
               KW_CO2 = ( 1.64d0 * U10M - 13.69d0 ) * 1d-3
               ALPHA  = 0.50d0
            ENDIF

            ! Get the temperature-dependent kinematic viscosity of water
            IF (TK_SURF <= 273.15 ) THEN
               VISC_H2O = 1.787d-2      ! [cm2/s]
            ELSE IF (TK_SURF > 273.15 .AND. TK_SURF <= 278.15 ) THEN
               VISC_H2O = 1.518d-2
            ELSE IF (TK_SURF > 258.15 .AND. TK_SURF <= 283.15 ) THEN
               VISC_H2O = 1.307d-2
            ELSE IF (TK_SURF > 283.15 .AND. TK_SURF <= 287.15 ) THEN
               VISC_H2O = 1.139d-2
            ELSE IF (TK_SURF > 287.15 .AND. TK_SURF <= 293.15 ) THEN
               VISC_H2O = 1.002d-2
            ELSE IF (TK_SURF > 293.15 .AND. TK_SURF <= 298.15 ) THEN
               VISC_H2O = 0.89d-2
            ELSE IF (TK_SURF > 298.15 ) THEN
               VISC_H2O = 0.797d-2 

            ! Calculate the diffusivites of CO2 and POP in water
            DW_CO2 = ( 13.26 * 1d-5 ) / (( VISC_H2O*1d2 )**1.14d0 *
     &               (V_CO2)**0.589d0 )   ! [cm2/s]

            DW_POP = ( 13.26 * 1d-5 ) / (( VISC_H2O*1d2 )**1.14d0 *
     &               (V_POP)**0.589d0 )   ! [cm2/s]
 
            ! Calculate the Schmidt numbers for CO2 and POP
            SCH_CO2 = VISC_H2O / DW_CO2    ! [unitless]
            SCH_POP = VISC_H2O / DW_POP    ! [unitless]

            ! Calculate the water-side MTC for POP
            KW_POP = KW_CO2 * ( SCH_POP / SCH_CO2 ) ** (-ALPHA) ! [cm/s]

            ! Calculate the temperature-dependent dimensionless Henry's Law constant
            KAW_T = KAW_298 * EXP((-DEL_HW/R) * ((1d0/TK_SURF) - 
     &                  (1d0/298d0)))       ! [unitless]

            ! Now calculate the overall air-water MTC
            KOL_T = 1d0 / ( 1d0/KW_POP + 1d0/(KA_POP*KAW_T) ) ! [cm/s]

            ! Calculate Flux in ng/m2/day ! 
            FLUX = KOL_T * 3600d0 * 24d0 * ( C_DISS - POPG/KAW_T ) *
     &             MWPOP * 1d12 / 100d0

            ! Convert to an emission rate in kg/s for returning to pops_mod.f
            ! Only return it if it's positive
!Bug fix (H. Amos, 8 Mar 2013)
!            EPOP_LAKE(I,J) = MAX(FLUX * GET_AREA_M2( J ) / 24d0 / 3600d0 
!     &                         / 1d12, 0d0 )
            EPOP_LAKE(I,J) = MAX(FLUX * GET_AREA_M2(I,J,1) / 24d0  
     &                         / 3600d0 / 1d12, 0d0 )

            ! Multiply the mass emitted by the fraction of land that is water
            EPOP_LAKE(I,J) = FRAC_LAKE * EPOP_LAKE(I,J)
         
            ! If the flux is positive, then the direction will be from the soil to the air. 
            ! Store this in a diagnostic array. 
            ! If the flux is zero or negative, store it in a separate array.

            IF (ND53 > 0) THEN
               IF (FLUX > 0d0) THEN 
                  ! Store total mass emitted from soil [kg] in ND53 diagnostic.
                  ! (We don't care about the mass in the other direction right now)
                  AD53(I,J,6) = AD53(I,J,6) + (EPOP_LAKE(I,J) * DTSRCE)  
                  ! Store positive flux
                  AD53(I,J,10) = AD53(I,J,10) + FLUX
                  ! Make sure negative flux diagnostic has nothing added to it
                  AD53(I,J,11) = AD53(I,J,11) + 0d0
                  ! Store the soil/air fugacity ratio
                  AD53(I,J,15) = AD53(I,J,15) + (C_DISS / (POPG/KAW_T))
             
               ELSE IF (FLUX <= 0d0) THEN
                  ! Store the negative flux
                  AD53(I,J,11) = AD53(I,J,11) + FLUX
                  ! Add nothing to positive flux or mass diagnostics
                  AD53(I,J,6) = AD53(I,J,6) + 0d0
                  AD53(I,J,10) = AD53(I,J,10) + 0d0
                  ! Continue to store the fugacity ratio
                  AD53(I,J,15) = AD53(I,J,15) + (C_DISS / (POPG/KAW_T))
               ENDIF
            ENDIF
            ENDIF

            ELSE 
            ! We are not on land or the land is covered with ice or snow
            ! or we are land but there is no water
               EPOP_LAKE(I,J) = 0d0
               FLUX = 0d0
               AD53(I,J,10) = 0d0 
               AD53(I,J,6) = 0d0
               AD53(I,J,11) = 0d0
               AD53(I,J,15) = 0d0
            ENDIF
         ENDIF
        
      ENDDO
      ENDDO

#endif

      END SUBROUTINE LAKEEMISPOP
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: vegemispop
!
! !DESCRIPTION: Subroutine VEGEMISPOP is the subroutine for secondary
!  POP emissions from soils.
!\\
!\\
! !INTERFACE:
!
         SUBROUTINE VEGEMISPOP( POP_SURF, EPOP_VEG )
!
! !USES:
!
      USE TIME_MOD, ONLY      : GET_TS_EMIS, GET_TS_CHEM
      USE DAO_MOD, ONLY       : IS_ICE, IS_LAND, IS_WATER
      USE DAO_MOD, ONLY       : FRLAKE
      USE MODIS_LAI_MOD, ONLY : GC_LAI, GC_LAI_CM, GC_LAI_PM ! hma, 7Mar2013
!      USE LAI_MOD, ONLY       : ISOLAI, MISOLAI, PMISOLAI
      USE DAO_MOD, ONLY       : TSKIN, AIRVOL, T, TS, SNOW
      USE DAO_MOD, ONLY       : IS_ICE, IS_LAND, SNOMAS
      USE GRID_MOD, ONLY      : GET_AREA_CM2
      USE GRID_MOD, ONLY      : GET_AREA_M2
      USE TRACER_MOD, ONLY    : STT
      USE TRACERID_MOD, ONLY  : IDTPOPG
      USE DIAG53_MOD, ONLY    : AD53 !ND53,  LD53,  AD53
      USE GET_POPSINFO_MOD, ONLY    : GET_POP_XMW, GET_POP_HSTAR
      USE GET_POPSINFO_MOD, ONLY    : GET_POP_DEL_Hw, GET_POP_DEL_H
      USE GET_POPSINFO_MOD, ONLY    : GET_POP_KOA

      ! Replace include statements (H. Amos, 7Mar2013)
      USE CMN_SIZE_MOD  ! Size parameters
      USE CMN_DIAG_MOD  ! Diagnostic arrays & switches
!#     include "CMN_SIZE"     ! Size parameters
!#     include "CMN_DIAG"     ! Diagnostic arrays & switches

!
! !INPUT PARAMETERS:
!
!      LOGICAL,                INTENT(IN) :: LSECEMISPOP
      REAL*8, DIMENSION(:,:), INTENT(IN) :: POP_SURF

!
! !OUTPUT PARAMETERS:
!
      REAL*8, DIMENSION(:,:), INTENT(OUT):: EPOP_VEG
!
! !REMARKS:
!
!  
!
! !REVISION HISTORY:
!  
!EOP
!------------------------------------------------------------------------------
!BOC
! 
! !LOCAL VARIABLES:
!     

      INTEGER                :: I, J, L
      REAL*8                 :: POPG, POPG_GL!, SUNCOSVALUE
      REAL*8                 :: FRAC_SNOWFREE_LAND, TK_SURF 
      REAL*8                 :: KLA_T, FLUX, K_MT
      REAL*8                 :: LEAF_CONC, DTSRCE, KOA_T
      REAL*8                 :: DIFF, FLEAF, FAIR, DS
      REAL*8                 :: DAB_F, DC, DAL, PC, UC
      REAL*8                 :: ZLEAF, ZAIR, TK
      REAL*8                 :: DTCHEM, NEWLEAF, E_KDEG
      REAL*8                 :: LAI, KAW_T, KOW_T
      REAL*8                 :: FUG_R, NEW_LEAF, DLA
      LOGICAL                :: IS_SNOW_OR_ICE, IS_LAND_OR_ICE
      LOGICAL                :: IS_SNOWFREE_LAND

! !External functions
      REAL*8, EXTERNAL   :: BOXVL ! Grid box vol in cm3
!
! !DEFINED PARAMETERS:

      ! Delta H for POP [kJ/mol]. Delta H is enthalpy of phase transfer
      ! from gas phase to OC. For now we use Delta H for phase transfer 
      ! from the gas phase to the pure liquid state. 
      ! For PHENANTHRENE: 
      ! this is taken as the negative of the Delta H for phase transfer
      ! from the pure liquid state to the gas phase (Schwarzenbach,
      !  Gschwend, Imboden, 2003, pg 200, Table 6.3), or -74000 [J/mol].
      ! For PYRENE:
      ! this is taken as the negative of the Delta H for phase transfer
      ! from the pure liquid state to the gas phase (Schwarzenbach,
      ! Gschwend, Imboden, 2003, pg 200, Table 6.3), or -87000 [J/mol]. 
      ! For BENZO[a]PYRENE:
      ! this is also taken as the negative of the Delta H for phase transfer
      ! from the pure liquid state to the gas phase (Schwarzenbach,
      ! Gschwend, Imboden, 2003, pg 452, Prob 11.1), or -110,000 [J/mol]
      REAL*8  :: DEL_H   

      ! Delta H for POP [kJ/mol].  
      ! For PHENANTHRENE: 
      ! this is the Delta H for phase transfer
      ! from air to water (Schwarzenbach,
      !  Gschwend, Imboden, 2003, pg 200, Table 6.3), or 47000 [J/mol].
      ! For PYRENE: 43000 [J/mol]
      REAL*8  :: DEL_HW     

      ! R = universal gas constant for adjusting KOA for temp: 8.3145 [J/mol/K OR m3*Pa/K/mol]
      REAL*8, PARAMETER     :: R          = 8.31d0 

      ! Molecular weight
      ! For phe, 0.17823 kg/mol
      REAL*8    :: MW   
!
      ! For PHENANTHRENE:
      ! log KOA_298 = 7.64, or 4.37*10^7 [unitless]
      ! For PYRENE:
      ! log KOA_298 = 8.86, or 7.24*10^8 [unitless]
      ! For BENZO[a]PYRENE:
      ! log KOA_298 = 11.48, or 3.02*10^11 [unitless]
      ! (Ma et al., J. Chem. Eng. Data, 2010, 55:819-825).
      REAL*8  :: KOA_298 

      ! For PHENANTHRENE:
      ! log KAW_298 = -2.76, or 1.74*10-3 [unitless]
      ! For PYRENE:
      ! log KAW_298 = -3.27, or 5.37*10-4 [unitless]
      REAL*8    :: KAW_298 

      ! Set volume fractions of octanol and water in surface and reservoir
      ! leaf compartments [unitless]
      REAL*8, PARAMETER  :: OCT_SURF = 0.8d0
      REAL*8, PARAMETER  :: OCT_RES  = 0.02d0
      REAL*8, PARAMETER  :: H2O_RES  = 0.7d0

      ! Set thickness of different leaf compartments. Volumes calculated by
      ! multiplying thicknesses by leaf area index
      REAL*8, PARAMETER  :: SURF_THICK = 2d-6  ! m
      REAL*8, PARAMETER  :: RES_THICK  = 250d-6 ! m

      ! Set transfer velocity and diffusion coefficient values
      REAL*8, PARAMETER     :: UAB_F        = 9d0     ![m/h]

      ! Set soil degradation rate
      REAL*8, PARAMETER     :: DEGR       = 3.5d-5  ![/h]

      REAL*8                :: DUM

      !=================================================================
      ! VEGEMISPOP begins here!
      !=================================================================

      DUM = 1.0
      DEL_H = GET_POP_DEL_H(DUM)
      KOA_298 = GET_POP_KOA(DUM)
      DEL_HW = GET_POP_DEL_Hw(DUM)
      KAW_298 = GET_POP_HSTAR(DUM)
      MW = GET_POP_XMW(DUM)

!      IF (.NOT. LSECEMISPOP) THEN

!         EPOP_VEG = 0D0

!      ELSE

      ! Emission timestep [s]
      DTSRCE  = GET_TS_EMIS() * 60d0

      ! Chemistry timestep [h]
      DTCHEM = GET_TS_CHEM()

      DO J=1, JJPAR
      DO I=1, IIPAR

#if defined ( GEOS_5 )
         ! Set logicals
         ! Is grid box covered by land/ice or by water? (IS_LAND_OR_ICE)
         ! IS_LAND will return non-ocean boxes but may still contain lakes
         ! If land, is it covered by snow/ice? (IS_SNOW_OR_ICE)
         IS_LAND_OR_ICE = ( (IS_LAND(I,J)) .OR. (IS_ICE(I,J)) ) 
         IS_SNOW_OR_ICE = 
     &      ( (IS_ICE(I,J)) .OR. (IS_LAND(I,J) .AND. SNOMAS(I,J)>10d0) )
#elif defined ( GCAP )
         ! GCAP uses different snow variable
         IS_LAND_OR_ICE = (  (IS_LAND(I,J)) .OR. (IS_ICE(I,J)) ) 
         IS_SNOW_OR_ICE = 
     &      ( (IS_ICE(I,J)) .OR. (IS_LAND(I,J) .AND. SNOW(I,J)>10d0) )
#endif

         ! Do soils routine only if we are on land that is not covered with
         ! snow or ice
         IF ((IS_LAND_OR_ICE) .AND. .NOT. ( IS_SNOW_OR_ICE )) THEN

            ! Get fraction of grid box covered by leaf surface area
            ! Do not consider different vegetation types for now
            LAI = GC_LAI(I,J) !hma, 7Mar2013
!            LAI = ISOLAI(I,J)
!            WRITE(6,*) 'ISOLAI=', ISOLAI(I,J)

            IF ( LAI > 0d0 ) THEN

#if defined ( GEOS_5 )
            ! Get surface skin temp [K]
            TK_SURF = TSKIN(I,J)
#elif defined ( GCAP )
            TK_SURF = TS(I,J)
#endif
            ! Get air temp [K]
            TK = T(I,J,1)

            ! Get gas phase air POP concentration at surface in mol/m3
            POPG = MAX( STT(I,J,1,IDTPOPG), SMALLNUM ) ![kg in gridbox]
            ! kg / (0.178 kg/mol) /m3 in gridbox
            POPG = POPG / MW / AIRVOL(I,J,1) ! mol/m3

            ! Only consider partitioning into leaf surface (and not reservoir)
            ! for now following Mackay et al 2006 Environ Sci & Pollut Res
            ! Include reservoir when land-atm models become dynamic

            ! Assume that all leaf surfaces contain an average lipid content 
            ! of 80% (Mackay et al 2006) 

            ! Get leaf concentration 
            ! Convert to mol/m3      
            ! kg deposited to leaf in 1 yr * / 0.178 kg/mol
            ! / area grid box (m2) / surface thickness m 
!Bug fix (H. Amos, 8 Mar 2013)
!            LEAF_CONC = POP_SURF(I,J) / MW / GET_AREA_M2( J ) /
!     &                  SURF_THICK! mol/m3   
            LEAF_CONC = POP_SURF(I,J) / MW / GET_AREA_M2(I,J,1) /
     &                  SURF_THICK! mol/m3   

            ! Check concentration in leaves by assuming a density similar to water
            ! (1 g/cm3)

            ! No degradation/metabolism for now. Just scale leaf concentrations
            ! to match flux observations (are there any??)
            NEWLEAF = LEAF_CONC/1d4  !SCALING FACTOR
            !WRITE (6,*) 'Leaf conc pg/g=', NEWLEAF*MW*1d15/100**3

            ! Define temperature-dependent KOA:
            KOA_T = KOA_298 * EXP((-DEL_H/R) * ((1d0/TK_SURF) - 
     &                  (1d0/298d0))) 

            ! Calculate the temperature-dependent dimensionless Henry's Law constant
            KAW_T = KAW_298 * EXP((-DEL_HW/R) * ((1d0/TK_SURF) - 
     &                  (1d0/298d0)))       ! [unitless]

            ! Estimate the temperature-dependent dimensionless octanol-water constant
            KOW_T = KOA_T * KAW_T

            ! Define dimensionless leaf surface-air partition coefficient 
            ! (mol/m3 leaf / mol/m3 air) 
            ! KLA = foct_surf * Koa
            KLA_T = OCT_SURF * KOA_T
!            KLA_T = MAX( KLA_T, SMALLNUM )

            ! Calculate fugacities from concentrations by dividing by "Z" values, 
            ! or the fugacity capacity in mol/m3*Pa following Mackay and Paterson, 1991

            ! fleaf = Cleaf * R * T / KSA [Pa]
            ! where Cleaf is in mol/m3, R is in Pa * m3 / mol K
            ! T is in K and KLA is dimensionless
            FLEAF = NEWLEAF * R * TK_SURF / KLA_T

            ! fair = Cair * R * T [Pa]
            ! where Cair is in mol/m3, R is in Pa * m3 / mol K and T is in K
            FAIR = POPG * R * TK

            ! Calculate the fugacity gradient [Pa]
            ! If the gradient is negative, fair is larger and the POP will diffuse from air to soil
            ! If the gradient is positive, fsoil is larger and POP will diffuse from soil to air
            DIFF = FLEAF - FAIR
            FUG_R = FLEAF/FAIR

            ! Calculate "Z" values from fugacities. 
            ! Z is the fugacity capacity in mol/m3*Pa. C = Z*f, so Z = C/f
            ZAIR  = POPG / FAIR ! (mol/m3) / (Pa)
            ZLEAF = NEWLEAF / FLEAF ! (mol/m3) / (Pa)

            ! Calculate the "D" value, or the transfer coefficient that describes the movement of 
            ! POP between phases (Mackay and Paterson, 1991, Cousins and Mackay 2000, internal report).
            ! [mol/h*Pa]
            ! The D value for leaf surface-air gas diffusion is given by
            ! Dla = 1 / (1/Dc + 1/(Dab-f)) [mol/(Pa*h)]
            ! where Dab-f is the boundary layer diffusion  [mol/h*Pa]
            ! given by Dab-f = As * L * Uab-f * Za
            ! where As is the area of the land surface [m2], L is the leaf area index [m2/m2],
            ! Uab-f is a mass transfer coefficient for surface-air boundary layer diffusion [m/h],
            ! and Za is the fugacity capacity of the air [mol/(m3*Pa)]
            ! Dc is the cuticle diffusion, given by 
            ! Dc = As * L * Uc * Zf
            ! where As and L are as above, Uc is the cuticle mass transfer coefficient [m/h],
            ! and Zf is the fugacity capacity of the leaf surface (mol/(m3*Pa))
            ! Uc is given by 
            ! Uc = 3600 * Pc * 1/Kaw
            ! where Pc is the cuticle permeance (m/s) and Kaw is the dimensionless air-water
            ! partition coefficient. 
            ! Pc is given by 
            ! Log Pc = ((0.704 * log Kow - 11.2) + (-3.47 - 2.79 * logMW + 0.970 log Kow)) / 2
            ! (an average of two equations)
  
            ! Need to define each D value
            ! DAB_F:
            !  m/h * mol/(m3*Pa)  =  (mol/h*Pa*m2)
            DAB_F = UAB_F * ZAIR  !  mol/(h*Pa*m2)

            ! Calculate PC and then Uc in order to calculate Dc
            ! PC, UC are calculated according to Cousins and Mackay, Chemosphere, 2001, Table 2
            PC = 10** (( 0.704d0 * LOG (KOW_T) - 11.2d0 ) + ( -3.47d0 
     &           -2.79d0* LOG(MW*1000d0) + 0.97d0 * LOG(KOW_T))/2d0) ![m/s]

            UC = 3600d0 * PC * 1d0/KAW_T ! [m/h]

            DC = UC * ZLEAF    ! mol/(h*Pa*m2)

            ! Now calculate overall transfer  ! mol/(h*Pa*m2)
            DLA = 1d0 / (1d0/DC + 1d0/DAB_F)

            ! Calculate Flux in mol/h/m2
            FLUX = DLA * DIFF

            ! Change to units of ng/m2/d for storage
            FLUX = FLUX * 24d0 * MW * 1d12 

            ! Convert to an emission rate in kg/s for returning to pops_mod.f
            ! Only want to add rates that are positive
!Bug fix (H. Amos, 8 Mar 2013)
!            EPOP_VEG(I,J) = MAX(FLUX * GET_AREA_M2( J ) * LAI / 24d0 / 
!     &                         3600d0 / 1d12, 0d0)
            EPOP_VEG(I,J) = MAX(FLUX * GET_AREA_M2(I,J,1) * LAI / 24d0 / 
     &                         3600d0 / 1d12, 0d0)

            ! If the flux is positive, then the direction will be from the soil to the air. 
            ! Store this in a diagnostic array. 
            ! If the flux is zero or negative, store it in a separate array.

            IF (ND53 > 0) THEN
               IF (FLUX > 0d0) THEN 
                  ! Store total mass emitted from soil [kg] in ND53 diagnostic.
                  ! (We don't care about the mass in the other direction right now)
                  AD53(I,J,7) = AD53(I,J,7) + (EPOP_VEG(I,J) * DTSRCE)  
                  ! Store positive flux
                  AD53(I,J,12) = AD53(I,J,12) + FLUX
                  ! Make sure negative flux diagnostic has nothing added to it
                  AD53(I,J,13) = AD53(I,J,13) + 0d0
                  ! Store the soil/air fugacity ratio
                  AD53(I,J,16) = AD53(I,J,16) + FUG_R
             
               ELSE IF (FLUX <= 0d0) THEN
                  ! Store the negative flux
                  AD53(I,J,13) = AD53(I,J,13) + FLUX
                  ! Add nothing to positive flux or mass diagnostics
                  AD53(I,J,7) = AD53(I,J,7) + 0d0
                  AD53(I,J,12) = AD53(I,J,12) + 0d0
                  ! Continue to store the fugacity ratio
                  AD53(I,J,16) = AD53(I,J,16) + FUG_R
               ENDIF
            ENDIF

            ELSE 
            ! We are not on land or the land is covered with ice or snow
               EPOP_VEG(I,J) = 0d0
               FLUX = 0d0
               AD53(I,J,7) = 0d0 
               AD53(I,J,12) = 0d0
               AD53(I,J,13) = 0d0
               AD53(I,J,16) = 0d0
            ENDIF
         ENDIF
        
      ENDDO
      ENDDO


      END SUBROUTINE VEGEMISPOP
!EOC

!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_land_pops
!
! !DESCRIPTION: Subroutine INIT\_LAND\_POPS allocates and zeroes all 
!  module arrays.
!\\
!\\
! !INTERFACE: 
!
!      SUBROUTINE INIT_LAND_POPS
! 
! !USES:
!
!      USE ERROR_MOD,    ONLY : ALLOC_ERR

!      USE CMN_SIZE_MOD     ! Size parameters
! 
! !REVISION HISTORY:
!  28 Aug 2012 - C. Friedman  - Initial version based on INIT_LAND_MERCURY
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
!      INTEGER              :: AS
!      LOGICAL, SAVE :: IS_INIT = .FALSE. 

      !=================================================================
      ! INIT_LAND_POPS begins here!
      !=================================================================

      ! Allocate arrays
c$$$      ALLOCATE( TRANSP( IIPAR, JJPAR ), STAT=AS )
c$$$      IF ( AS /= 0 ) CALL ALLOC_ERR( 'TRANSP' )
c$$$      TRANSP = 0d0

!      END SUBROUTINE INIT_LAND_POPS
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!

! !IROUTINE: cleanup_land_pops
!
! !DESCRIPTION: Subroutine CLEANUP\_LAND\_POPS deallocates all module 
!  arrays. 
!\\
!\\
! !INTERFACE: 
!
!      SUBROUTINE CLEANUP_LAND_POPS
! 
! !REVISION HISTORY:
!  28 Aug 2012 - C. Friedman  - Initial version based on INIT_LAND_MERCURY
!EOP
!------------------------------------------------------------------------------
!BOC

      ! Deallocate arrays
!      IF ( ALLOCATED( TRANSP      ) ) DEALLOCATE( TRANSP      )

!      END SUBROUTINE CLEANUP_LAND_POPS
!EOC
      END MODULE LAND_POPS_MOD
