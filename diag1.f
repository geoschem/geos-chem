! $Id: diag1.f,v 1.22 2009/08/19 17:05:47 ccarouge Exp $
      SUBROUTINE DIAG1 
!
!******************************************************************************
!  Subroutine DIAG1 accumulates diagnostic quantities every NDIAG minutes
!  (bmy, bey, 6/16/98, 1/28/04)
!
!  NOTES:
!  (1 ) This subroutine was reconstructed from gmg's version of (10/10/97)
!  (2 ) GISS-specific code has been eliminated (bmy, 3/15/99)
!  (3 ) UWND, VWND, WW no longer needs to be passed (bmy, 4/7/99)
!  (4 ) Use F90 syntax for declarations, etc (bmy, 4/7/99)
!  (5 ) Remove counter KWACC...this is now redundant (bmy, 11/5/99)
!  (6 ) ND31, ND33, ND35, ND67, and ND69 now use dynamically 
!        allocatable arrays declared in "diag_mod.f". (bmy, 3/9/00)
!  (7 ) LTOTH is now an allocatable array in "diag_mod.f". (bmy, 3/17/00)
!  (8 ) Add parallel loops over tracer where expedient (bmy, 5/4/00)
!  (9 ) Updated comments and diagnostics list.  Also add more parallel
!        loops for ND31 and ND68.  (bmy, 6/21/00)
!  (10) Use NTRACE to dimension STT_VV instead of NNPAR (bmy, 10/17/00)
!  (11) Removed obsolete code from 10/17/00 (bmy, 12/21/00)
!  (12) Updated diagnostic list & comments, cosmetic changes (bmy, 6/19/01)
!  (13) Updated diagnostic list & comments (bmy, 9/4/01)
!  (14) Now reference AVGW from "dao_mod.f", and make sure it is allocated
!        before we reference it in the ND68 diagnostic.  Also reference PBL, 
!        PS, AIRDEN from "dao_mod.f". (bmy, 9/25/01)
!  (15) Removed obsolete code from 9/01 (bmy, 10/23/01)
!  (16) Renamed ND33 to "ATMOSPHERIC COLUMN SUM OF TRACER", since this is
!        a sum over all levels and not just in the troposphere.  Also
!        removed more obsolete code from 9/01.  Now use P(I,J)+PTOP instead
!        of PS, since that is the way to ensure that we use will be used
!        consistently.  Remove reference to PS from "dao_mod.f"(bmy, 4/11/02)
!  (17) Replaced all instances of IM with IIPAR and JM with JJPAR, in order
!        to prevent namespace confusion for the new TPCORE.  Also removed
!        obsolete, commented-out code.  Also now replaced reference to
!        P(IREF,JREF) with P(I,J). (bmy, 6/25/02)
!  (18) Replaced references to P(I,J) with call to GET_PEDGE(I,J,1) from
!        "pressure_mod.f"  Eliminated obsolete commented-out code from
!        6/02. (dsa, bdf, bmy, 8/20/02)
!  (19) Now reference AD, and BXHEIGHT from "dao_mod.f".  Removed obsolete 
!        code.  Now refEerence IDTOX from "tracerid_mod.f". (bmy, 11/6/02)
!  (20) Now replace DXYP(J) with routine GET_AREA_M2 from "grid_mod.f"
!        (bmy, 2/4/03)
!  (21) Now compute PBL top for ND67 for GEOS-4/fvDAS.  Also now include
!        SCALE_HEIGHT from header file "CMN_GCTM". (bmy, 6/23/03)
!  (22) Now references N_TRACERS, STT, and ITS_A_FULLCHEM_SIM from
!        "tracer_mod.f" (bmy, 7/20/04)
!  (23) Fixed ND67 PS-PBL for GCAP and GEOS-5 met fields (swu, bmy, 6/9/05)
!  (24) Now archive ND30 diagnostic for land/water/ice flags (bmy, 8/18/05)
!  (25) Now reference XNUMOL from "tracer_mod.f" (bmy, 10/25/05)
!  (26) Remove support for GEOS-1 and GEOS-STRAT met fields (bmy, 8/4/06)
!  (27) Added count for time in the troposphere - array AD54 (phs, 9/22/06)
!  (28) Now only archive O3 in ND45 and ND47 at chem timsteps (phs, 1/24/07)
!  (29) Bug fix: Update ND30 for both GEOS-3 and otherwise.  Also now save
!        3-D pressure edges in ND31 instead of PS-PTOP.  Revert to the !
!        pre-near-land ND30 diagnostic algorithm. (bmy, 1/28/04)
!  (30) Use LTO3 for O3 in ND45. (ccc, 7/20/09)
!******************************************************************************
!  List of GEOS-CHEM Diagnostics (bmy, 10/25/05)
!
!  FLAG  DIM'S    QUANTITY                                     UNITS
!  ---- --------  --------                                     -----
!  ND01 (I,J,L,3) RADON 222 - LEAD 210 - BERYLLIUM 7 SOURCE    [kg/s]
!
!  ND02 (I,J,L,3) RADON 222 - LEAD 210 - BERYLLIUM 7 DECAY     [kg/s]
!
!  ND03   ----    Free Diagnostic    
!  ND04   ----    Free Diagnostic  
!
!  ND05           PROD/LOSS for SULFATE CHEMISTRY QUANTITIES
!       (I,J,L)    P(SO2) from DMS + OH                        [kg S]
!       (I,J,L)    P(SO2) from DMS + NO3                       [kg S]
!       (I,J,L)    Total P(SO2)                                [kg S]
!       (I,J,L)    P(MSA) from DMS                             [kg S]
!       (I,J,L)    P(SO4) gas phase                            [kg S]         
!       (I,J,L)    P(SO4) aqueous phase                        [kg S]
!       (I,J,L)    Total P(SO4)                                [kg S]
!       (I,J,L)    L(OH) by DMS                                [kg OH] 
!       (I,J,L)    L(NO3) by DMS                               [kg NO3]
!       (I,J,L)    L(H2O2)                                     [kg H2O2]
!
!  ND06  (I,J)    DESERT DUST EMISSIONS                        [kg]
!
!  ND07           SOURCES OF BLACK CARBON & ORGANIC CARBON    
!       (I,J)      BLACK CARBON from anthro sources            [kg]
!       (I,J)      BLACK CARBON from biomass burning           [kg]
!       (I,J)      BLACK CARBON from biofuels                  [kg]
!       (I,J,L)    Hydrophilic BC from Hydrophobic BC          [kg]
!       (I,J)      ORGANIC CARBON from anthro sources          [kg]
!       (I,J)      ORGANIC CARBON from biomass burning         [kg]
!       (I,J)      ORGANIC CARBON from biofuels                [kg]
!       (I,J)      ORGANIC CARBON from biogenic sources        [kg]
!       (I,J)      Hydrophilic OC from Hydrophobic OC          [kg]
!
!  ND08  (I,J)    SEA SALT EMISSIONS                           [kg] 
!
!  ND09   ----    Free Diagnostic 
!  ND10   ----    Free Diagnostic 
!
!  ND11           ACETONE SOURCES & SINKS
!       (I,J)      Acetone source from MONOTERPENES            atoms C/cm2/s
!       (I,J)      Acetone source from METHYL BUTENOL          atoms C/cm2/s
!       (I,J)      Acetone source from DIRECT EMISSION         atoms C/cm2/s
!       (I,J)      Acetone source from DRY LEAF MATTER         atoms C/cm2/s
!       (I,J)      Acetone source from GRASSLANDS              atoms C/cm2/s
!       (I,J)      Acetone source from OCEANS                  atoms C/cm2/s
!       (I,J)      Acetone sink   from OCEANS                  atoms C/cm2/s
!
!  ND12 (I,J)     FRACTION OF BOUNDARY LAYER OCCUPIED BY       unitless
!                 LEVEL L (for new emissions in "setemis.f")
!
!  ND13           TROPOSPHERIC SULFUR EMISSIONS                   
!       (I,J)      DMS                                         kg S
!       (I,J,L)    Aircraft SO2             ( 1 <= L <= ND13 ) kg S
!       (I,J,L)    Anthro SO2               ( 1 <= L <= 2    ) kg S
!       (I,J)      Biomass SO2                                 kg S
!       (I,J,L)    Non-eruptive volcano SO2 ( 1 <= L <= ND13 ) kg S       
!       (I,J,L)    Eruptive     volcano SO2 ( 1 <= L <= ND13 ) kg S       
!       (I,J,L)    Anthro SO4               ( 1 <= L <= 2    ) kg S
!       (I,J)      Anthro NH3                                  kg NH3
!       (I,J)      Biomass NH3                                 kg NH3
!       (I,J)      Biofuel NH3                                 kg NH3
!
!  ND14 (I,J,L,N) UPWARD MASS FLUX DUE TO WET CONVECTION       kg/s 
!                 ( 1 <= L <= ND14 )
!
!  ND15 (I,J,L,N) MASS CHANGE DUE TO BOUNDARY-LAYER MIXING     kg/s 
!                 ( 1 <= L <= ND15 )
!
!  ND16 (I,J,L)   AREAL FRACTION OF LARGE-SCALE PRECIP         unitless
!       (I,J,L)   AREAL FRACTION OF CONVECTIVE  PRECIP         unitless
!                 ( 1 <= L <= ND16 )
!
!  ND17 (I,J,L,4) RAINOUT FRACTION IN LARGE-SCALE PRECIP       unitless
!       (I,J,L,4) RAINOUT FRACTION IN CONVECTIVE  PRECIP       unitless
!                 ( 1 <= L <= ND17 )
!
!  ND18 (I,J,L,4) WASHOUT FRACTION IN LARGE-SCALE PRECIP       unitless
!       (I,J,L,4) WASHOUT FRACTION IN CONVECTIVE  PRECIP       unitless
!                 ( 1 <= L <= ND18 )
!
!  ND19   ----    Free Diagnostic 
!
!  ND20 (I,J,L,2) SAVE O3 PROD/LOSS RATES TO DISK              molec/cm3/s
!                 ( 1 <= L <= LLTROP )
!
!  ND21 (I,J,L,3) CLOUD OPTICAL DEPTHS AND CLOUD FRACTIONS     unitless
!                 ( 1 <= L <= ND21 )
!
!  ND22 (I,J,L,6) J-VALUES: NO2, HNO3, H2O2, CH2O, OH, O3      s^-1
!                 ( 1 <= L <= ND22 )
!
!  ND23   ----    METHYL CHLOROFORM (CH3CCl3) LIFETIME         years
!
!  ND24 (I,J,L,N) EASTWARD  MASS FLUX BY TRANSPORT             kg/s 
!                 ( 1 <= L <= ND24 )
!
!  ND25 (I,J,L,N) NORTHWARD MASS FLUX BY TRANSPORT             kg/s 
!                 ( 1 <= L <= ND25 )
!
!  ND26 (I,J,L,N) UPWARD    MASS FLUX BY TRANSPORT             kg/s 
!                 ( 1 <= L <= ND26 )
!
!  ND27 (I,J,3)   STRATOSPHERIC INFLUX of NOX, OX, HNO3        kg/s 
!
!  ND28 (I,J,10)  BIOMASS BURNING EMISSIONS:                   molec/cm2/s
!                  NOX,  CO,   ALK4, ACET, MEK,
!                  ALD2, PRPE, C3H8, CH2O, C2H6
!
!  ND29           TROPOSPHERIC CO EMISSIONS                    molec/cm2/s
!       (I,J)      Anthropogenic
!       (I,J)      Biomass 
!       (I,J)      Biofuel
!       (I,J)      CO produced from monoterpenes
!       (I,J)      CO produced from methanol
!
!  ND30 (I,J)     PLOT LAND MAP FOR GISS or GEOS models        integers
!
!  ND31 (I,J)     SURFACE PRESSURE - PTOP                      mb
!
!  ND32           TROPOSPHERIC NOx EMISSIONS                   
!       (I,J,L)    Aircraft       ( 1 <= L <= LAIREMS )        molec/cm2/s
!       (I,J,L)    Anthropogenic  ( 1 <= L <= NOXEXTENT )      molec/cm2/s
!       (I,J)      Biomass                                     molec/cm2/s
!       (I,J)      Fertilization                               molec/cm2/s
!       (I,J,L)    Lightning      ( 1 <= L <= LCONVM  )        molec/cm2/s
!       (I,J)      Soils                                       molec/cm2/s     
!       (I,J)      Upper Boundary                              molec/cm2/s
!
!  ND33 (I,J,N)   ATMOSPHERIC COLUMN SUM OF TRACER             kg         
!
!  ND34 (I,J,10)  BIOFUEL BURNING EMISSIONS:                   molec/cm2/s
!                  NOx,  CO,   ALK4, ACET, MEK, 
!                  ALD2, PRPE, C3H8, CH2O, C2H6
!      
!  ND35 (I,J,N)   TRACER AT 500 HPA ( L = 9 )                  v/v
!
!  ND36 (I,J,9)   ANTHROPOGENIC EMISSIONS (for NSRCX == 3)     molec/cm2/s 
!                  NOx, CO,   ALK4, ACET, 
!                  MEK, ALD2, PRPE, C3H8, C2H6 
!                  OR CH3I EMISSIONS (for NSRCX == 2)          ng/m2/s
!
!  ND37 (I,J,L,4) FRACTION OF TRACER SCAVENGED BY CLOUD
!                 UPDRAFTS IN MOIST CONVECTION                 unitless
!                 ( 1 <= L <= ND37 )
!
!  ND38 (I,J,L,N) LOSS OF TRACER IN MOIST CONVECTION           kg/s     
!                 ( 1 <= L <= ND38 )
!
!  ND39 (I,J,L,N) LOSS OF TRACER IN AEROSOL WET DEPOSITION     kg/s
!                 ( 1 <= L <= ND39 )
!
!  ND40   ----    Free diagnostic                              kg/m2
!
!  ND41 (I,J)     AFTERNOON PBL DEPTH (1200-1600 LT)           m
!
!  ND42   ----    Free Diagnostic
!
!  ND43 (I,J,L,2) OH CONCENTRATIONS (from SMVGEAR) and         molec/cm3/s
!                 NO CONCENTRATIONS (from SMVGEAR)             v/v
!                 ( 1 <= L <= ND43 )
! 
!  ND44 (I,J,M)   DRYDEP FLUXES                                molec/cm2/s 
!       (I,J,M)   DRYDEP VELOCITIES                            cm/s
!                 ( M = NUMDEP )
! 
!  ND45 (I,J,L,N) TRACER CONCENTRATION, AVERAGED BETWEEN       v/v
!                 OTH_HR1 and OTH_HR2 (from "input.geos") 
!                 ( 1 <= L <= ND45 )
!
!  ND46 (I,J,5)   BIOGENIC EMISSIONS of ISOP, PRPE, ACET,      molec/cm2/s
!                 MONOTERPENES, METHYL BUTENOL
!
!  ND47 (I,J,L,N) DAILY (24-h) AVERAGE TRACER CONCENTRATIONS   v/v
!                 ( 1 <= L <= ND47 )
!
!  ND48 (I,J,L,N) TIME SERIES AT N = NNSTA LOCATIONS 
!                  FOR MS = 0, store tracer concentrations     v/v
!                  FOR MS = 1, store the following
!                  N = 1, store O3   RURAL                     molec/cm3
!                  N = 2, store OH   RURAL                     molec/cm3
!                  N = 3, store NOy  RURAL                     molec/cm3
!                  N = 4, store Drydep Vel for NO2             cm/s
!                  N = 5, store Drydep Vel for O3              cm/s
!                  N = 6, store Drydep Vel for PAN             cm/s
!                  N = 7, store Drydep Vel for HNO3            cm/s
!                  N = 8, store Drydep Vel for H2O2            cm/s
!                  N = 9, store NO   RURAL                     molec/cm3
!                  (get I,J,L,MS,N from "inptr.ctm")
!
!  ND49 (I,J,L,N) 3-D INSTANTANEOUS TRACER TIMESERIES          v/v
!                 SAVED IN BINARY PUNCH FILE FORMAT
!                 (get I,J,L,N from "timeseries.dat")
!
!  ND50 (I,J,L,N) 3-D 24-h AVERAGE TRACER TIMESERIES           v/v
!                 SAVED IN BINARY PUNCH FILE FORMAT
!                 (get I,J,L,N from "timeseries.dat")
!
!  ND51 (I,J,L,N) 3-D 1-4pm AVERAGE TRACER TIMESERIES          v/v 
!                 SAVED IN BINARY PUNCH FILE FORMAT             
!                 (get I,J,L,N from "timeseries.dat")
!
!  ND52  (I,J,L,N) HO2 aerosol uptake coefficient (gamma)     unitless
!
!  ND53   ----    Free Diagnostic 
!
!  ND54 (I,J,L)   Time in troposphere (fraction of total time) unitless
!
!  ND55           TROPOPAUSE DIAGNOSTICS
!       (I,J)      Tropopause level number                     unitless
!       (I,J)      Tropopause height                           km
!       (I,J)      Tropopause pressure                         hPa
!
!  ND60   ----    Free Diagnostic
!  
!  ND61   ----    Free Diagnostic 
!
!  ND62 (I,J,N)   INSTANTANEOUS COLUMN MIXING RATIO            v/v
!
!  ND63   ----    Free Diagnostic
!
!  ND65 (I,J,L,N) PRODUCTION & LOSS OF SELECTED                molec/cm3/s
!                 CHEMCIAL SPECIES (see "prodloss.dat")  
!
!  ND66           DAO 3-D FIELDS ( 1 <= L <= ND66 )
!       (I,J,L)    UWND   : U-winds                            m/s
!       (I,J,L)    VWND   : V-winds                            m/s
!       (I,J,L)    TMPU   : Temperature                        K
!       (I,J,L)    SPHU   : Specific humidity                  g H20/kg air   
!       (I,J,L)    CLDMAS : Convective Mass Flux               kg/m2/s 
!
!  ND67           DAO SURFACE FIELDS                            
!       (I,J)      HFLUX   : sensible heat flux from surface   W/m2  
!       (I,J)      RADSWG  : solar radiation @ ground          W/m2  
!       (I,J)      PREACC  : total precip. @ ground            mm/day 
!       (I,J)      PRECON  : conv.  precip. @ ground           mm/day 
!       (I,J)      TS      : surface air temperature           K    
!       (I,J)      RADSWT  : solar radiation @ atm. top        W/m2  
!       (I,J)      USTAR   : friction velocity                 m/s   
!       (I,J)      Z0      : surface roughness height          m    
!       (I,J)      PBL     : planetary boundary layer depth    hPa   
!       (I,J)      CLDFRC  : column cloud fraction             0 - 1  
!       (I,J)      U10M    : U-winds @ 10 meters altitude      m/s   
!       (I,J)      V10M    : V-winds @ 10 meters altitude      m/s   
!       (I,J)      PS-PBL  : Boundary Layer Top Pressure       hPa
!       (I,J)      ALBD    : Surface Albedo                    unitless
!       (I,J)      PHIS    : Geopotential Heights              m         
!       (I,J)      CLTOP   : Cloud Top Height                  levels
!       (I,J)      TROPP   : Tropopause pressure               hPa 
!       (I,J)      SLP     : Sea Level pressure                hPa
!       (I,J)      TSKIN   : Ground / sea surface temperature  K
!       (I,J)      PARDF   : Photosyn active diffuse rad       W/m2
!       (I,J)      PARDR   : Photosyn active direct rad        hPa
!       (I,J)      GWETTOP : Top soil wetness                  unitless
!
!  ND68           GRID BOX QUANTITIES ( 1 <= L <= ND68 )
!       (I,J,L)    BXHEIGHT : (grid box heights)               m 
!       (I,J,L)    AD       :(air mass in grid box)            kg
!       (I,J,L)    AVGW     : mixing ratio of water vapor      v/v
!       (I,J,L)    N_AIR    : Air number density               molec air/m3
!  
!  ND69 (I,J)     DXYP     : grid box surface areas            m2
!
!  ND70  ----     DEBUG PRINTOUT to stdout file via routine DEBUG_MSG
!*****************************************************************************
!
      ! References to F90 modules
      USE DAO_MOD,        ONLY : AD,  AIRDEN, AVGW,     BXHEIGHT 
      USE DAO_MOD,        ONLY : PBL, IS_ICE, IS_WATER, IS_LAND, IS_NEAR
      USE DIAG_MOD,       ONLY : AD30, AD31, AD33, AD35, AD45, AD54 
      USE DIAG_MOD,       ONLY : AD47, AD67, AD68, AD69, LTOTH, LTO3
      USE GRID_MOD,       ONLY : GET_AREA_M2
      USE PRESSURE_MOD,   ONLY : GET_PEDGE
      USE TIME_MOD,       ONLY : ITS_TIME_FOR_CHEM
      USE TRACER_MOD,     ONLY : N_TRACERS, STT, TCVV
      USE TRACER_MOD,     ONLY : ITS_A_FULLCHEM_SIM
      USE TRACER_MOD,     ONLY : XNUMOLAIR
      USE TRACERID_MOD,   ONLY : IDTOX
      USE TROPOPAUSE_MOD, ONLY : ITS_IN_THE_TROP

      IMPLICIT NONE

#     include "CMN_SIZE"  ! Size parameters
#     include "CMN_DIAG"  ! Diagnostic arrays & parameters
#     include "CMN_O3"    ! FRACO3
#     include "CMN_GCTM"  ! Physical constants

      ! Local variables
      LOGICAL            :: AVGW_ALLOCATED, IS_FULLCHEM, IS_CHEM
      INTEGER            :: I, J, K, L, N, NN, IREF, JREF, LN45
      REAL*8             :: FDTT, XLOCTM, AREA_M2
      REAL*8             :: STT_VV(IIPAR,JJPAR,LLPAR,N_TRACERS)

      !=================================================================
      ! DIAG1 begins here!
      !=================================================================
      
      ! Is it a fullchem run?
      IS_FULLCHEM = ITS_A_FULLCHEM_SIM()
      IS_CHEM     = ITS_TIME_FOR_CHEM()
      
      ! Compute conc. in mixing ratio for ND35, ND45, ND47 diagnostics 
      IF ( ND35 > 0 .or. ND45 > 0 .or. ND47 > 0 ) THEN
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, N )
         DO N = 1, N_TRACERS
         DO L = 1, LLPAR
         DO J = 1, JJPAR
         DO I = 1, IIPAR
            STT_VV(I,J,L,N) = 
     &           MAX( STT(I,J,L,N) * TCVV(N) / AD(I,J,L), 0d0 )
         ENDDO
         ENDDO
         ENDDO
         ENDDO
!$OMP END PARALLEL DO
      ENDIF

      !================================================================= 
      ! ND30: Land/water/ice flags
      !=================================================================
      IF ( ND30 > 0 ) THEN
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J )
         DO J = 1, JJPAR
         DO I = 1, IIPAR
            IF ( IS_WATER( I, J ) ) AD30(I,J) = AD30(I,J) + 0e0
            IF ( IS_LAND ( I, J ) ) AD30(I,J) = AD30(I,J) + 1e0
            IF ( IS_ICE  ( I, J ) ) AD30(I,J) = AD30(I,J) + 2e0
         ENDDO
         ENDDO
!$OMP END PARALLEL DO
      ENDIF

      !================================================================= 
      ! ND31: Surface pressure diagnostic (PS - PTOP) in hPa
      !=================================================================
      IF ( ND31 > 0 ) THEN
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L )
         DO L = 1, LD31
         DO J = 1, JJPAR
         DO I = 1, IIPAR
            AD31(I,J,L) = AD31(I,J,L) + GET_PEDGE( I, J, L ) 
         ENDDO
         ENDDO
         ENDDO
!$OMP END PARALLEL DO
      ENDIF

      !================================================================= 
      ! ND33: Atmospheric column sum of tracer [kg]
      !=================================================================  
      IF ( ND33 > 0 ) THEN   
!$OMP PARALLEL DO 
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, N )
         DO N = 1, N_TRACERS
         DO L = 1, LLPAR
         DO J = 1, JJPAR
         DO I = 1, IIPAR
            AD33(I,J,N) = AD33(I,J,N) + STT(I,J,L,N)
         ENDDO
         ENDDO
         ENDDO
         ENDDO  
!$OMP END PARALLEL DO
      ENDIF

      !================================================================= 
      ! ND35: 500 HPa fields.  
      !
      ! NOTES: 
      ! (1 ) Use level 9 for both GEOS-1 and GEOS-STRAT.
      !       They are both close to 500 hPa (bmy, 4/7/99)  
      !================================================================= 
      IF ( ND35 > 0 ) THEN
         L = 9
!$OMP PARALLEL DO 
!$OMP+DEFAULT( SHARED ) 
!$OMP+PRIVATE( I, J, N )
         DO N = 1, N_TRACERS
         DO J = 1, JJPAR
         DO I = 1, IIPAR
            AD35(I,J,N) = AD35(I,J,N) + STT_VV(I,J,L,N)
         ENDDO
         ENDDO
         ENDDO
!$OMP END PARALLEL DO
      ENDIF

      !================================================================= 
      ! ND45: Tracer (V/V) at level 1 to level LD45, averaged over 
      !       the time period OTH_HR1 to OTH_HR2.  
      !
      ! Store pure O3 (as opposed to Ox) in the NTRACE+1 location of 
      ! AD45.  FRACO3(I,J,L) is the fraction of Ox that is actually O3, 
      ! and is calculated in subroutine OHSAVE.
      !
      !  NOTES: 
      !  (1 ) AD45 array replaces the AIJ array for this diagnostic 
      !        (bmy, 3/22/99)
      !  (2 ) Add parallel loop over tracers (bmy, 5/4/00)
      !  (3 ) Use LTO3 and not LTOTH for O3. (ccc, 7/20/09)
      !================================================================= 
      IF ( ND45 > 0 ) THEN
!$OMP PARALLEL DO 
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, N )
         DO N = 1, N_TRACERS
            DO L = 1, LD45  
            DO J = 1, JJPAR 
            DO I = 1, IIPAR
               AD45(I,J,L,N) = AD45(I,J,L,N) + 
     &              STT_VV(I,J,L,N) * LTOTH(I,J)
            ENDDO
            ENDDO
            ENDDO

            ! With TS_DIAG, the diags are automatically calculated at the
            ! right time (= same time done for all processes). (ccc, 5/21/09)
!            IF ( N == IDTOX .and. IS_FULLCHEM .and. IS_CHEM ) THEN
            IF ( N == IDTOX .and. IS_FULLCHEM ) THEN
               DO L = 1, LD45
               DO J = 1, JJPAR
               DO I = 1, IIPAR
!---- Prior to (ccc, 7/20/09) ----
!                  AD45(I,J,L,N_TRACERS+1) = AD45(I,J,L,N_TRACERS+1) +
!     &                 ( STT_VV(I,J,L,N) * FRACO3(I,J,L) * LTOTH(I,J) )
                  AD45(I,J,L,N_TRACERS+1) = AD45(I,J,L,N_TRACERS+1) +
     &                 ( STT_VV(I,J,L,N) * FRACO3(I,J,L) * LTO3(I,J) )
               ENDDO            
               ENDDO            
               ENDDO            
            ENDIF
         ENDDO   
!$OMP END PARALLEL DO 
      ENDIF

      !================================================================= 
      ! ND47: Tracer (V/V) at level 1 to level LD45, 
      ! averaged from 0-24 hours
      !
      ! Added parallel loop over tracers (bmy, 5/4/00)
      !================================================================= 
      IF ( ND47 > 0 ) THEN
!$OMP PARALLEL DO 
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, N )
         DO N = 1, N_TRACERS
            DO L = 1, LD47
            DO J = 1, JJPAR
            DO I = 1, IIPAR
               AD47(I,J,L,N) = AD47(I,J,L,N) + STT_VV(I,J,L,N) 
            ENDDO   
            ENDDO
            ENDDO
            
            ! With TS_DIAG, the diags are automatically calculated at the
            ! right time (= same time done for all processes). (ccc, 5/21/09)
!            IF ( N == IDTOX .and. IS_FULLCHEM .and. IS_CHEM ) THEN
            IF ( N == IDTOX .and. IS_FULLCHEM ) THEN
               DO L = 1, LD47
               DO J = 1, JJPAR
               DO I = 1, IIPAR
                  AD47(I,J,L,N_TRACERS+1) = AD47(I,J,L,N_TRACERS+1) +
     &                 ( STT_VV(I,J,L,N) * FRACO3(I,J,L) )
               ENDDO
               ENDDO
               ENDDO
            ENDIF
         ENDDO 
!$OMP END PARALLEL DO   
      ENDIF

      !================================================================= 
      ! ND54: Count time the box was tropospheric
      !================================================================= 
      IF ( ND54 > 0 ) THEN

            DO L = 1, LD54
            DO J = 1, JJPAR
            DO I = 1, IIPAR
               IF ( ITS_IN_THE_TROP(I,J,L) )
     &              AD54(I,J,L) = AD54(I,J,L) + 1.
            ENDDO   
            ENDDO
            ENDDO
            
      ENDIF
          
          
      !=================================================================  
      ! ND67: Store PBL top pressure [hPa]
      !=================================================================  
      IF ( ND67 > 0 ) THEN
         DO J = 1, JJPAR
         DO I = 1, IIPAR

#if   defined( GEOS_3 )

            ! PBL is in [hPa], subtract from PSurface
            AD67(I,J,13) = AD67(I,J,13) + GET_PEDGE(I,J,1) - PBL(I,J) 

#else
            
            ! PBL is in [m], use hydrostatic law to get [hPa]
            AD67(I,J,13) = AD67(I,J,13) + 
     &           ( GET_PEDGE(I,J,1) * EXP( -PBL(I,J) / SCALE_HEIGHT ) )

#endif

         ENDDO
         ENDDO
      ENDIF

      !================================================================= 
      ! ND68: Quantity 1: BXHEIGHT(I,J,L) in meters 
      !       Quantity 2: AD(I,J,L)       in kg
      !       Quantity 3: AVGW(I,J,L)     in v/v
      !       Quantity 4: N_AIR(I,J,L)    in molecules air / m3
      !
      ! NOTE: AVGW is now an allocatable array from "dao_mod.f"
      !================================================================= 
      IF ( ND68 > 0 ) THEN

         ! Set a flag for whether AVGW is allocated or not
         AVGW_ALLOCATED = ALLOCATED( AVGW )

!$OMP PARALLEL DO 
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L )
         DO L = 1, LD68
         DO J = 1, JJPAR
         DO I = 1, IIPAR
            AD68(I,J,L,1) = AD68(I,J,L,1) + BXHEIGHT(I,J,L)
            AD68(I,J,L,2) = AD68(I,J,L,2) + AD(I,J,L)
            AD68(I,J,L,4) = AD68(I,J,L,4) + AIRDEN(L,I,J) * XNUMOLAIR

            ! Make sure AVGW is now allocated (bmy, 9/25/01)
            IF ( AVGW_ALLOCATED ) THEN
              AD68(I,J,L,3) = AD68(I,J,L,3) + AVGW(I,J,L)
           ENDIF
         ENDDO
         ENDDO
         ENDDO
!$OMP END PARALLEL DO
      ENDIF

      !================================================================= 
      ! ND69: Grid box surface areas [m2]
      !
      ! NOTE: Only save areas on the first timestep since the
      !       grid box surface areas are a time-invariant field.
      !================================================================= 
      IF ( ND69 > 0 ) THEN
         DO J = 1, JJPAR
            AREA_M2 = GET_AREA_M2( J )

            DO I = 1, IIPAR
               AD69(I,J,1) = AREA_M2
            ENDDO
         ENDDO
      ENDIF
      
      ! Return to calling program
      END SUBROUTINE DIAG1

