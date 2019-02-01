!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: dao_mod.F
!
! !DESCRIPTION: Module DAO\_MOD contains both arrays that hold DAO met fields, 
!  as well as subroutines that compute, interpolate, or otherwise process 
!  DAO met field data. 
!\\
!\\
! !INTERFACE: 
!
      MODULE DAO_MOD
!
! !USES:
!
      USE CMN_SIZE_MOD           ! Size parameters
      USE PhysConstants          ! Physical constants
      USE PRECISION_MOD          ! For GEOS-Chem Precision (fp)

      IMPLICIT NONE
      PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
      PUBLIC  :: AVGPOLE
      PUBLIC  :: AIRQNT
      PUBLIC  :: GET_COSINE_SZA
      PUBLIC  :: GET_OBK
      PUBLIC  :: INTERP
      PUBLIC  :: IS_LAND
      PUBLIC  :: IS_WATER
      PUBLIC  :: IS_ICE
      PUBLIC  :: IS_NEAR      
      PUBLIC  :: SET_DRY_SURFACE_PRESSURE
      PUBLIC  :: Set_Met_AgeOfAir
#if defined( ESMF_ ) || defined( EXTERNAL_GRID )
      PUBLIC  :: GIGC_Cap_Tropopause_Prs
#endif
!
! !REVISION HISTORY:
!  26 Jun 2010 - R. Yantosca - Initial version
!  (1 ) Added sea level pressure (SLP) met field for GEOS-3 (bmy, 10/10/00)
!  (2 ) Moved MAKE_QQ to "wetscav_mod.f" (bmy, 10/12/00)
!  (3 ) Now get LWI from ALBEDO for GEOS-3 in routines IS_LAND and
!        IS_WATER (bmy, 4/4/01)
!  (4 ) Define OPTDEP allocatable array for GEOS-3 -- this is the grid 
!        box optical depth and is now stored as a met field (bmy, 8/15/01)
!  (5 ) Updated comments (bmy, 9/4/01)
!  (6 ) Now make AVGW an allocatable module array.  Also replace obsolete
!        parameters {IJL}GCMPAR with IIPAR,JJPAR,LLPAR. (bmy, 9/27/01)
!  (7 ) Remove arguments LMAKEPW, PW, and LM from AIRQNT (bmy, 10/3/01)
!  (8 ) Remove obsolete code from 9/01 (bmy, 10/23/01)
!  (9 ) Bug fixes in IS_LAND and IS_WATER.  Also cosmetic changes and 
!        updated some comments. (mje, bmy, 1/9/02)
!  (10) Now add additional array PSC2 in order to pass to TPCORE, which will
!        fix the mixing ratio bug.  Compute PSC2 in subroutine INTERP.
!        Now bundle "convert_units.f" into "dao_mod.f".  Updated comments.
!        (bmy, 3/27/02)
!  (11) Updated comments (bmy, 5/28/02)
!  (12) Replaced all instances of IM with IIPAR and JM with JJPAR, in order
!        to prevent namespace confusion for the new TPCORE (bmy, 6/25/02)
!  (13) Eliminated PS, PSC arrays.  Now reference "pressure_mod.f".  Also
!        updated AIRQNT for hybrid grid.  Added routine MAKE_RH to this
!        module. (dsa, bdf, bmy, 8/27/02)
!  (14) Added arrays AD, BXHEIGHT, and T to "dao_mod.f".  Also removed 
!        obsolete code from 8/02 from several module routines.  Now 
!        references "error_mod.f".  Remove all references to QQ, it is now
!        declared in "wetscav_mod.f".  (bmy, 11/8/02)
!  (15) Now references "grid_mod.f".  Also added PHIS field, which was
!        formerly stored as PALTD in "CMN".  Added bug fix in routine
!        AVGPOLE for 1x1 nested grid. (bmy, 3/11/03)
!  (16) Added SUNCOSB array for SMVGEAR II.  Also removed KZZ array, since
!        that is now obsolete. (bmy, 4/28/03)
!  (17) Now moved MAKE_CLDFRC into "a6_read_mod.f".  Added HKBETA, HKETA, 
!        TSKIN, GWETTOP, ZMEU, ZMMD, ZMMU, PARDF, PARDR fields for 
!        GEOS-4/fvDAS. (bmy, 6/25/03)
!  (18) Added CLDFRC, RADSWG, RADLWG, SNOW arrays (bmy, 12/9/03)
!  (19) Added routine COPY_I6_FIELDS w/ parallel DO-loops (bmy, 4/13/04)
!  (20) Now also allocate AVGW for offline aerosol simulation (bmy, 9/28/04)
!  (21) AVGPOLE now uses NESTED_CH and NESTED_NA cpp switches (bmy, 12/1/04)
!  (22) Now modified for GEOS-5 and GCAP met fields (swu, bmy, 5/25/05)
!  (23) Now allocate SNOW and GWET for GCAP (bmy, 8/17/05)
!  (24) Now also add TSKIN for GEOS-3 (tmf, bmy, 10/20/05)
!  (25) Modifications for near-land formulation (ltm, bmy, 5/16/06)
!  (26) Remove support for GEOS-1 and GEOS-STRAT met fields (bmy, 8/4/06)
!  (27) Modified for variable tropopause (phs, bdf, 9/14/06)
!  (28) Add in extra fields for GEOS-5.  Updated COSSZA.  Now cap var trop 
!        at 200hPa near poles in INTERP (bmy, phs, 9/18/07)
!  (29) Bug fix in INIT_DAO for CMFMC array (bmy, jaf, 6/11/08)
!  (30) Add heat flux EFLUX for GEOS5. (lin, ccc, 5/29/09)
!  (31) Add fractions of land and water, FRLAND, FROCEAN, FRLANDIC, FRLAKE 
!        for methane (kjw, 8/18/09)
!  (32) Bug fix in AVGPOLE (bmy, 12/18/09)
!  (33) Remove obsolete SUNCOSB array (bmy, 4/28/10)
!  16 Aug 2010 - R. Yantosca - Added ProTeX headers
!  18 Aug 2010 - R. Yantosca - Added modifications for MERRA data
!  18 Aug 2010 - R. Yantosca - Move CMN_SIZE, CMN_DIAG to top of module
!  25 Aug 2010 - R. Yantosca - Now read LWI (land/water/ice) for MERRA met
!  05 Oct 2011 - R. Yantosca - Add SUNCOS_30 array to hold the cos(SZA)
!                              computed @ 30 mins after each GMT hour.
!  07 Oct 2011 - R. Yantosca - Rename SUNCOS30 to SUNCOS_MID, which is the
!                              cos(SZA) at the midpt of the chemistry timestep
!  06 Feb 2012 - R. Yantosca - Add modifications for GEOS-5.7.x met fields
!  06 Feb 2012 - R. Yantosca - Split up INIT_DAO into several routines
!  07 Feb 2012 - M. Payer    - Add subroutine GET_COSINE_SZA to compute sun
!                              angles at the current time and 5 hours prior to
!                              the current time (for the PARANOX ship emissions
!                              plume model) (R. Yantosca)
!  28 Feb 2012 - R. Yantosca - Removed support for GEOS-3
!  01 Mar 2012 - R. Yantosca - Now references the new grid_mod.F90
!  06 Mar 2012 - R. Yantosca - Now allocate TO3 for all met fields
!  21 Nov 2012 - R. Yantosca - Removed met fields now contained in State_met
!  21 Nov 2012 - R. Yantosca - Remove functions INIT_DAO_GCAP, INIT_DAO_GEOS4,
!                              INIT_DAO_GEOS5, INIT_DAO_GEOS57, INIT_DAO_MERRA
!  27 Nov 2012 - R. Yantosca - Removed obsolete AIRQNT_FULLGRID routine and 
!                              obsolete arrays AIRDEN_FULLGRID, T_FULLGRID
!  28 Nov 2012 - R. Yantosca - Removed SUNCOS, SUNCOS_MID, SUNCOS_MID_5hr
!  28 Nov 2012 - R. Yantosca - Removed routines INIT_DAO, INIT_DAO_DERIVED, and
!                              CLEANUP_DAO; we have no more allocatable arrays
!  14 Mar 2013 - M. Payer    - Restored routines AIRQNT_FULLGRID, INIT_DAO,
!                              CLEANUP_DAO and arrays AIRDEN_FULLGRID and
!                              T_FULLGRID. They are required to correct vertical
!                              regridding of OH for offline simulations.
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
!  07 Nov 2014 - M. Yannetti - Added PRECISION_MOD
!  23 Feb 2015 - E. Lundgren - Change AIRQNT and formulas to take into
!                              account the presence of water vapor
!  26 Feb 2015 - E. Lundgren - Move calculations of RH and AVGW from MAKE_RH
!                              and MAKE_AVGW to AIRQNT. Update algorithms
!                              to account for the presence of water vapor
!                              and to use [g water / kg moist air] for the
!                              definition of specific humidity (SPHU).
!  16 Apr 2015 - R. Yantosca - Removed AIRQNT_FULLGRID, this was used to regrid
!                              OH fields, but that is now done by HEMCO
!  16 Apr 2015 - R. Yantosca - Also remove INIT_DAO and CLEANUP_DAO, since
!                              these allocated/deallocated the FULLGRID arrays
!  29 Nov 2016 - R. Yantosca - grid_mod.F90 is now gc_grid_mod.F90
!  24 Aug 2017 - M. Sulprizio- Remove support for GCAP, GEOS-4, GEOS-5 and MERRA
!  08 Mar 2018 - E. Lundgren - Add GIGC_Cap_Tropopause_Prs from GCHP gchp_util
!  21 Dec 2018 - M. Sulprizio- Add new routine Set_Met_AgeOfAir
!EOP
!------------------------------------------------------------------------------
!BOC
      CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: avgpole
!
! !DESCRIPTION: Subroutine AVGPOLE computes average quantity near polar caps, 
!  defined by (J = 1, 2) and (J = JJPAR-1, JJPAR).  
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE AVGPOLE( Z ) 
!
! !USES:
!
      USE GC_GRID_MOD, ONLY : GET_AREA_M2
!
! !INPUT/OUTPUT PARAMETERS: 
!
      REAL(fp), INTENT(INOUT) :: Z(IIPAR,JJPAR)   ! Quantity to be averaged 
                                                !  over the pole (usually PS)
! 
! !REVISION HISTORY: 
!  30 Jan 1998 - R. Yantosca - Initial version
!  (1 ) AVGPOLE is written in Fixed-Form Fortran 90.  Use F90 syntax
!        for declarations, etc (bmy, 4/14/99)
!  (2 ) MAIN now passes the Harvard CTM variable for surface area of
!        a gridbox, DXYP(JJPAR), to AVGPOLE.  Use window offset
!        J+J0 when accessing DXYP.  Add JJPAR to the parameter list.
!  (3 ) Added this routine to "dao_mod.f" (bmy, 6/27/00)
!  (4 ) Updated comments (bmy, 4/4/01)
!  (5 ) Now replaced DXYP(J) with routine GET_AREA_M2 of "grid_mod.f"
!        Now also return immediately if GRID1x1 is selected. (bmy, 3/11/03)
!  (6 ) Now use cpp switches NESTED_CH and NESTED_NA to denote nested
!        grids...GRID1x1 can now also denote a global grid (bmy, 12/1/04)
!  (7 ) Also need to RETURN for 0.5 x 0.666 nested grid simulations 
!        (mpb, bmy, 12/18/09)
!  16 Aug 2010 - R. Yantosca - Added ProTeX headers
!  01 Mar 2012 - R. Yantosca - Now use GET_AREA_M2(I,J,L) from grid_mod.F90
!  26 Sep 2013 - R. Yantosca - Remove SEAC4RS C-preprocessor switch
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER :: I, J
      REAL(fp)  :: TOTAL_Z1, TOTAL_Z2, TOTAL_Z3, TOTAL_Z4
      REAL(fp)  :: TOTAL_A1, TOTAL_A2, TOTAL_A3, TOTAL_A4

      !=================================================================
      ! AVGPOLE begins here!                                                  
      !=================================================================

#if defined( GRID05x0625 ) || defined( GRID025x03125 )
#if defined( NESTED_CH ) || defined( NESTED_NA   ) || defined( NESTED_EU   ) || defined( NESTED_AS )
      ! NOTE: Only do this for 1x1 nested grids (bmy, 12/1/04)
      ! 1x1 window grid does not extend to poles
      RETURN
#endif
#endif

      TOTAL_Z1 = 0.
      TOTAL_Z2 = 0.
      TOTAL_Z3 = 0.
      TOTAL_Z4 = 0.
      TOTAL_A1 = 0.
      TOTAL_A2 = 0.
      TOTAL_A3 = 0.
      TOTAL_A4 = 0.

      DO I = 1, IIPAR
         TOTAL_Z1 = TOTAL_Z1 
     &            + GET_AREA_M2( I, 1,       1 ) * Z(I,1      )

         TOTAL_Z2 = TOTAL_Z2 
     &            + GET_AREA_M2( I, 2,       1 ) * Z(I,2      )

         TOTAL_Z3 = TOTAL_Z3 
     &            + GET_AREA_M2( I, JJPAR-1, 1 ) * Z(I,JJPAR-1)

         TOTAL_Z4 = TOTAL_Z4 
     &            + GET_AREA_M2( I, JJPAR,   1 ) * Z(I,JJPAR  )

         TOTAL_A1 = TOTAL_A1 + GET_AREA_M2( I, 1,       1 ) 
         TOTAL_A2 = TOTAL_A2 + GET_AREA_M2( I, 2,       1 )
         TOTAL_A3 = TOTAL_A3 + GET_AREA_M2( I, JJPAR-1, 1 )
         TOTAL_A4 = TOTAL_A4 + GET_AREA_M2( I, JJPAR,   1 )
      ENDDO

      DO I = 1, IIPAR
         Z(I,      1) = (TOTAL_Z1 + TOTAL_Z2) / (TOTAL_A1 + TOTAL_A2)
         Z(I,      2) = Z(I,1)
         Z(I,JJPAR-1) = (TOTAL_Z3 + TOTAL_Z4) / (TOTAL_A3 + TOTAL_A4)
         Z(I,JJPAR  ) = Z(I,JJPAR-1)
      ENDDO

      END SUBROUTINE AVGPOLE
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: airqnt
!
! !DESCRIPTION: Subroutine AIRQNT sets several members of State\_Met, the
!  meteorology object of derived type MetState, and optionally updates
!  the tracer concentrations to conserve tracer mass when air quantities
!  change. 
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE AIRQNT( am_I_Root, Input_Opt, State_Met,
     &                   State_Chm, RC,        Update_Mixing_Ratio )
!
! !USES:
!
      USE ErrCode_Mod    
      USE GC_GRID_MOD,   ONLY : Get_Area_M2
      USE Input_Opt_Mod, ONLY : OptInput
      USE State_Met_Mod, ONLY : MetState
      USE State_Chm_Mod, ONLY : ChmState
      USE PhysConstants, ONLY : AIRMW, AVO
      USE Pressure_Mod
      USE Time_Mod,      ONLY : Get_LocalTime
!
! !INPUT PARAMETERS:
!
      LOGICAL,        INTENT(IN)           :: am_I_Root ! On the root CPU?
      TYPE(OptInput), INTENT(IN)           :: Input_Opt ! Input Options obj
      LOGICAL,        INTENT(IN), OPTIONAL :: update_mixing_ratio 
                                                        ! Update mixing ratio?
                                                        ! Default is yes
!
! !INPUT/OUTPUT PARAMETERS:
!
      TYPE(MetState), INTENT(INOUT)        :: State_Met ! Met State obj
      TYPE(ChmState), INTENT(INOUT)        :: State_Chm ! Chemistry State obj
!
! !OUTPUT PARAMETERS:
!
      INTEGER,        INTENT(OUT)          :: RC        ! Success or failure?
!
! !REMARKS:
!  DAO met fields updated by AIRQNT:
!  ========================================================================
!  (1)  PEDGE     (REAL(fp)) : Moist air pressure at grid box bottom      [hPa]
!  (2)  PEDGE_DRY (REAL(fp)) : Dry air partial pressure at box bottom     [hPa]
!  (3)  PMID      (REAL(fp)) : Moist air pressure at grid box centroid    [hPa]
!  (4)  PMID_DRY  (REAL(fp)) : Dry air partial pressure at box centroid   [hPa]
!  (5)  PMEAN     (REAL(fp)) : Altitude-weighted mean moist air pressure  [hPa]
!  (6)  PMEAN_DRY (REAL(fp)) : Alt-weighted mean dry air partial pressure [hPa]
!  (7)  DELP      (REAL(fp)) : Delta-P extent of grid box                 [hPa]
!                              (Same for both moist and dry air since we
!                              assume constant water vapor pressure 
!                              across box) 
!  (8)  AIRDEN    (REAL(fp)) : Mean grid box dry air density            [kg/m^3]
!                              (defined as total dry air mass/box vol)
!  (9)  MAIRDEN   (REAL(fp)) : Mean grid box moist air density          [kg/m^3]
!                              (defined as total moist air mass/box vol)
!  (10) AD        (REAL(fp)) : Total dry air mass in grid box             [kg]
!  (11) ADMOIST   (REAL(fp)) : Total moist air mass in grid box           [kg]
!  (12) BXHEIGHT  (REAL(fp)) : Vertical height of grid box                [m]   
!  (13) AIRVOL    (REAL(fp)) : Volume of grid box                         [m^3]
!  (14) MOISTMW   (REAL(fp)) : Molecular weight of moist air in box     [g/mol]
!
! !REVISION HISTORY: 
!  30 Jan 1998 - R. Yantosca - Initial version
!  (1 ) AIRQNT is written in Fixed-Form Fortran 90.  Use F90 syntax
!        for declarations etc. (bmy, 4/14/99)
!  (2 ) AIRQNT can now compute PW from PS (if LMAKEPW=T) or PS from PW.
!  (3 ) AIRQNT should also be called after TPCORE, since TPCORE changes
!        the PW values.  AIRQNT must then be called to compute the post-TPCORE
!        values of AD, BXHEIGHT, AIRVOL, and AIRDEN.
!  (4 ) The AIRDEN and DELP arrays are now dimensioned as (LLPAR,IIPAR,JJPAR) 
!        for better efficiency when processing a whole (I,J) column layer by 
!        layer.  In FORTRAN, the best efficiency is obtained when the leftmost 
!        array index corresponds to the innermost loop.
!  (5 ) Remove PTOP from the arg list.  PTOP is now a parameter in 
!      "CMN_SIZE".  Also updated comments. (bmy, 2/22/00)
!  (6 ) Replace IM, JM, LM with IIPAR, JJPAR, LLPAR as loop boundaries.
!        This ensures that all quantities get defined up to the top of
!        the atmosphere. (bmy, 6/15/00)
!  (7 ) Added to "dao_mod.f" (bmy, 6/26/00)
!  (8 ) Updated comments (bmy, 4/4/01)
!  (9 ) P(IREF,JREF) is now P(I,J).  T(IREF,JREF,L) is now T(I,J,L).  Also
!        removed LM from the arg list, it is obsolete.  Also updated
!        comments. (bmy, 9/26/01)
!  (10) Remove PW -- it is now obsolete.  Also make PW a local variable,
!        we need to preserve the way it computes P so as to avoid numerical
!        drift. (bmy, 10/4/01)
!  (11) Removed obsolete code from 9/01 and 10/01 (bmy, 10/23/01)
!  (12) Removed LMAKEPW from arg list.  Added parallel DO loops (bmy, 11/15/01)
!  (13) Removed obsolete code from 11/01 (bmy, 1/9/02)
!  (14) Now rename G_SIGE to SIGE, and dimension it (1:LLPAR+1).  Updated
!        comments, cosmetic changes. (bmy, 4/4/02)
!  (15) Removed obsolete, commented-out code (bmy, 6/25/02)
!  (16) Removed PS, P, SIGE from the arg list for hybrid grid.  Now reference
!        routines GET_PEDGE and GET_BP from "pressure_mod.f".  Removed 
!        obsolete, commented-out code. (dsa, bdf, bmy, 8/27/02)
!  (17) Now only pass DXYP via the arg list -- the other arguments are actually
!        are already contained within "dao_mod.f" (bmy, 11/15/02)
!  (18) Now replace DXYP(JREF) with routine GET_AREA_M2 of "grid_mod.f".
!        (bmy, 3/11/03)
!  (19) Now move computation of DELP into main loop.  Also remove P, LOGP,
!        JREF, DSIG variables -- these are obsolete for fvDAS.  (bmy, 6/19/03)
!  16 Aug 2010 - R. Yantosca - Added ProTeX headers
!  01 Mar 2012 - R. Yantosca - Now use GET_AREA_M2(I,J,L) from grid_mod.F90
!  22 Oct 2012 - R. Yantosca - Now reference gigc_state_met_mod.F90
!  22 Oct 2012 - R. Yantosca - Renamed LOCAL_MET argument to State_Met
!  09 Nov 2012 - M. Payer    - Copy met field arrays to the State_Met derived
!                              type object
!  06 Nov 2014 - R. Yantosca - Now use State_Met%AIRDEN(I,J,L)
!  06 Nov 2014 - R. Yantosca - Now use State_Met%DELP(I,J,L)
!  07 Jan 2015 - E. Lundgren - Add new Met field MAIRDEN (moist air density)
!  13 Feb 2015 - E. Lundgren - Incorporate water vapor in calculations of 
!                              densities. Calculate dry air pressures.
!  23 Feb 2015 - E. Lundgren - Change box height formula to include water vapor
!  24 Feb 2015 - E. Lundgren - Add calculation of RH to subroutine. Use new
!                              formulate for H20 saturation pressure.
!  03 Mar 2015 - E. Lundgren - Set State_Met%TV (virtual temperature)
!  16 Apr 2015 - E. Lundgren - Add new State_Met variables PMEAN, PMEAN_DRY, 
!                              MOISTMW, and ADMOIST
!  28 Oct 2015 - E. Lundgren - Pass update_mixing_ratio flag, Input_Opt, and
!                              State_Chm to optionally update tracer units
!  30 Jun 2016 - R. Yantosca - Remove instances of STT.  Now get the advected
!                              species ID from State_Chm%Map_Advect.
!  15 Aug 2016 - E. Lundgren - Update conc w/ RH for all species
!  08 Jan 2018 - R. Yantosca - Now define query fields from State_Met here
!  17 Jan 2018 - R. Yantosca - Now compute tropopause quantities here, which
!                              allows us to retire GeosCore/chemgrid_mod.F
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      ! Scalars
      INTEGER             :: I,         J,          L
      INTEGER             :: L_CG,      L_TP,       N
      REAL(fp)            :: PEdge_Top, Esat,       Ev_mid,    Ev_edge
      REAL(fp)            :: Ev_mean,   PMEAN,      PMEAN_DRY, EsatA
      REAL(fp)            :: EsatB,     EsatC,      EsatD,     AREA_M2
      REAL(fp)            :: SPHU_kgkg, AVGW_moist, H,         FRAC
      REAL(fp)            :: Pb,        Pt
      LOGICAL             :: UpdtMR

      ! Arrays
      LOGICAL             :: IsLocNoon(IIPAR)
      REAL(f8)            :: LocTime  (IIPAR)
      REAL(f8)            :: Nearest  (IIPAR)

      ! Strings
      CHARACTER(LEN=255)  :: ErrMsg, ThisLoc
!
! DEFINED PARAMETERS:
!
      ! Conversion factors
      REAL(fp), PARAMETER :: TOFFSET = -2.7315e+2_fp   ! K   -> deg C
      REAL(fp), PARAMETER :: PCONV   = 1.00e+2_fp      ! hPa -> Pa
      REAL(fp), PARAMETER :: RHCONV  = 1.00e-2_fp      ! %   -> fraction
      REAL(fp), PARAMETER :: MCONV   = 1.00e-3_fp      ! g   -> kg

      ! Empirical parameters for water vapor saturation pressure
      ! (Source: Nordquist, 1973. "Numerical Approximiations of
      !  Selected Meteorological Parameters Related to Cloud Physics"
      !  Text quality clarifications from Stipanuk, 1973. "Algorithms
      !  for Generating a Skew-T, Log P Diagram and Computing Selected
      !  Meteorological Quantities") 
      REAL(fp), PARAMETER :: ESATP1  =  2.3832241e+1_fp
      REAL(fp), PARAMETER :: ESATP2  = -5.02808e+0_fp
      REAL(fp), PARAMETER :: ESATP3  =  8.1328e-3_fp
      REAL(fp), PARAMETER :: ESATP4  =  3.49149e+0_fp
      REAL(fp), PARAMETER :: ESATP5  = -1.3028844e+3_fp
      REAL(fp), PARAMETER :: ESATP6  = -1.3816e-7_fp
      REAL(fp), PARAMETER :: ESATP7  =  1.1344e+1_fp
      REAL(fp), PARAMETER :: ESATP8  = -3.03998e-2_fp
      REAL(fp), PARAMETER :: ESATP9  = -2.949076e+3_fp

      !=================================================================
      ! AIRQNT begins here! 
      !=================================================================

      ! Initialize
      RC       = GC_SUCCESS
      ErrMsg   = ''
      ThisLoc  = ' -> at AIRQNT (in module GeosCore/dao_mod.F)'

      ! Shadow variable for mixing ratio update
      UpdtMR = .TRUE.
      IF ( PRESENT(update_mixing_ratio) ) UpdtMR = update_mixing_ratio

      ! Pre-compute local solar time = UTC + Lon/15
      ! NOTE: Local time only depends on longitude
      DO I = 1, IIPAR
         LocTime(I) = GET_LOCALTIME( I, 1, 1 )
         Nearest(I) = ABS( LocTime(I) - 12.0_fp )
      ENDDO

      ! Find the box that is nearest to local noon
      IsLocNoon                      = .FALSE.
      IsLocNoon( MINLOC( Nearest ) ) = .TRUE.

      !=============================================================
      ! Update air quantities
      !=============================================================  

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED                                 )
!$OMP+PRIVATE( I,       J,         L,       Pedge_Top )
!$OMP+PRIVATE( EsatA,   EsatB,     EsatC,   EsatD     )
!$OMP+PRIVATE( Esat,    SPHU_kgkg, Ev_mid,  Ev_edge   )
!$OMP+PRIVATE( Area_M2, PMean,     Ev_mean, PMean_Dry )
      DO L = 1, LLPAR
      DO J = 1, JJPAR
      DO I = 1, IIPAR
                                                                          
         !=============================================================   
         ! Set wet air pressures [hPa]                                    
         ! (lower edge, delta, centroid, and spatially-weighted mean)     
         !=============================================================   
                                                                          
         ! Pressure at bottom edge of grid box [hPa]                      
         State_Met%PEDGE(I,J,L)      = GET_PEDGE( I, J, L   )

         ! Pressure at top edge of grid box [hPa]                         
         PEdge_Top                   = GET_PEDGE( I, J, L+1 )                 
                                                                          
         ! Pressure at bottom edge of grid box [hPa] (level LLPAR+1 only) 
         IF ( L == LLPAR ) THEN
            State_Met%PEDGE(I,J,L+1) = PEdge_Top
         ENDIF

         ! Pressure difference between top & bottom edges [hPa]           
         State_Met%DELP(I,J,L)       = State_Met%PEDGE(I,J,L)
     &                               - PEdge_Top

         ! Arithmetic average pressure in grid box [hPa] defined as       
         ! ( PEDGE(L) + PEDGE(L+1) ) / 2. This represents the grid box    
         ! mass centroid pressure. Use in the ideal gas law yields        
         ! a local air density at the centroid.                           
         State_Met%PMID(I,J,L)       = GET_PCENTER( I, J, L )
                                                                          
         !=============================================================   
         ! Calculate water vapor saturation pressure [hPa]                
         ! from temperature                                               
         !=============================================================   

         ! Prepare the intermediate terms
         EsatA = ESATP1 + ESATP2 * log10( State_Met%T(I,J,L) )            
         EsatB = ESATP3 * 10**( ESATP4 + ESATP5 / State_Met%T(I,J,L) )    
         EsatC = ESATP6 * 10**( ESATP7 + ESATP8 * State_Met%T(I,J,L) )    
         EsatD = ESATP9 / State_Met%T(I,J,L)                              

         ! Saturation water vapor pressure [hPa]
         Esat  = 10**( EsatA + EsatB + EsatC + EsatD )

         !=============================================================   
         ! Set AVGW, the water vapor volume mixing ratio from             
         ! specific humidity [mol H2O / mol dry air]                      
         !=============================================================   
         !                                                                
         ! vol H2O       dry air    H2O         kg H2O      kg wet air    
         ! ----------  = molec wt / molec wt * ---------- * ---------     
         ! vol dry air   [g/mol]    [g/mol]    kg wet air   kg dry air    
         !                                                                
         !      thus AVGW = AIRMW * SPHU_kgkg / ( H2OMW * (1-SPHU_kgkg) ) 
         !                                                                
         ! where,                                                         
         !        SPHU_kgkg = specific humidity in [kg/kg]                
         !                                                                
                                                                          
         ! Convert specific humidity from [g/kg] to [kg/kg]               
         SPHU_kgkg             = State_Met%SPHU(I,J,L) * MCONV

         ! Water vapor volume mixing ratio [v/v dry air]                  
         State_Met%AVGW(I,J,L) = AIRMW * SPHU_kgkg                        
     &                         / ( H2OMW * (1.0e+0_fp - SPHU_kgkg ) )
                                                                          
         !=============================================================   
         ! Calculate water vapor partial pressures [hPa] from relative    
         ! humidity
         !=============================================================   

         ! At vertical midpoint of grid box
         Ev_mid  = State_Met%RH(I,J,L) * RHCONV * Esat

         ! At bottom edge of grid box
         Ev_edge = State_Met%PEDGE(I,J,L) * Ev_mid
     &           / State_Met%PMID(I,J,L)

         !=============================================================   
         ! Set grid box height [m]                                        
         !=============================================================   
         !                                                                
         ! BXHEIGHT is the height (Delta-Z) of grid box (I,J,L)           
         ! in meters. It is calculated using the hypsometric equation.    
         ! A virtual temperature is calculated to enable use of           
         ! of moist air pressure with dry air molecular weight.           
         !                                                                
         !              Gas constant   virtual                            
         !              for dry air  * grid box                           
         !              [J/K/mol]      temp [K]       bottom edge P       
         ! height [m] = ----------------------- * ln(---------------)     
         !                     g [m/s2]                top edge P         
         !                                                                
         !  where,                                                        
         !                                                                
         !                         Grid box temperature [K]               
         !    Virtual  = -------------------------------------------      
         !    Temp [K]       H20 partial pressure          MW_H2O         
         !               1 - -------------------- * ( 1 - --------- )     
         !                    moist air  pressure         MW_dryair       
         !                                                                
         ! Source: Wallace and Hobbes "Atmospheric Science: An            
         !         Introductory Survey"                                   
         !                                                                
         ! Assume constant temperature and moisture across grid box. 
         !                                                                

         ! Grid box potential temperature [K]
         ! NOTE: Due to the parallelization, we cannot assume that
         ! State_Met%PEDGE(I,J,1) has been defined.  So always call
         ! GET_PEDGE(I,J,1) to return the proper surface pressure
         ! (bmy, 2/23/18)
         State_Met%THETA(I,J,L)    = State_Met%T(I,J,L)
     &                             * ( GET_PEDGE( I, J, 1 )
     &                             /   State_Met%PMID(I,J,L) )**0.286

         ! Grid box virtual temperature [K]                               
         State_Met%TV(I,J,L)       = State_Met%T(I,J,L) / (1 -
     &                               Ev_edge / State_Met%PEDGE(I,J,L)
     &                             * ( 1 - H2OMW / AIRMW ) )

         ! Grid box box height [m]                                        
         State_Met%BXHEIGHT(I,J,L) = Rdg0 * State_Met%TV(I,J,L)
     &                             * LOG( State_Met%PEDGE(I,J,L)
     &                             /      PEdge_Top )
                                                                          
         !==============================================================  
         ! Set grid box surface area [m2]                                 
         ! NOTE: We intend to gradually remove area use by                
         ! changing tracer units to mixing ratio                          
         !==============================================================  
                                                                          
#if defined( ESMF_ )                                                      
            AREA_M2                = State_Met%AREA_M2(I,J,1)
#else
            AREA_M2                = GET_AREA_M2( I,J,L )
#endif                                                                    
                                                                          
         !==============================================================  
         ! Set grid box volume [m3]                                       
         !==============================================================  
         !                                                                
         ! AIRVOL is the volume of grid box (I,J,L) in meters^3           
         !                                                                
         ! AREA_M2 is the Delta-X * Delta-Y surface area of grid          
         ! boxes (I,J,L=1), that is, at the earth's surface.              
         !                                                                
         ! Since the thickness of the atmosphere is much smaller          
         ! than the radius of the earth, we can make the "thin            
         ! atmosphere" approximation, namely:                             
         !                                                                
         !               (Rearth + h) ~ Rearth                            
         !                                                                
         ! Therefore, the Delta-X * Delta-Y surface area of grid          
         ! boxes that are above the earth's surface will be               
         ! approx. the same as AREA_M2.  Thus we are justified            
         ! in using AREA_M2 for grid boxes (I, J, L > 1)                  
         !                                                                
         State_Met%AIRVOL(I,J,L) = State_Met%BXHEIGHT(I,J,L) * AREA_M2
                                                                          
         !==============================================================
         ! Set grid box dry air partial pressures at grid box edges
         ! and at grid box centroid (arithmetic mean pressure level)
         ! [hPa]. Assume constant humidity across grid box.
         !============================================================== 

         ! Partial pressure of dry air at lower edge of grid box [hPa]
         State_Met%PEDGE_DRY(I,J,L) = State_Met%PEDGE(I,J,L) - Ev_edge

         ! Set dry air partial pressure for level LLPAR+1 lower edge
         IF ( L == LLPAR ) THEN
            State_Met%PEDGE_DRY(I,J,L+1) = Pedge_Top - Ev_edge
         ENDIF

         ! Partial pressure of dry air at box centroid [hPa]
         State_Met%PMID_DRY(I,J,L) = State_Met%PMID(I,J,L) - Ev_mid

         ! Set previous dry P difference to current dry P difference
         ! prior to updating with new met values
         State_Met%DP_DRY_PREV(I,J,L) = State_Met%DELP_DRY(I,J,L)

         ! Update dry pressure difference as calculated from the
         ! dry surface pressure
         State_Met%DELP_DRY(I,J,L) = GET_DELP_DRY(I,J,L)

         !==============================================================
         ! Set mass of dry air in grid box [kg]
         !==============================================================
         
         ! AD = mass of dry air in grid box (I,J,L) 
         !
         ! given by:        
         ! 
         !        Dry air pressure   Grid box   100 [Pa]   1 [s2]
         ! Mass = difference       * surface  * -------- * -----
         !        across grid        area       1 [hPa]    g [m]
         !        box [hPa]          [m2]                        
         !
         ! NOTES: 
         ! Dry air pressure difference across grid box is calculate
         ! from the surface dry pressure and GMAO A and B parameters
         ! (see GeosUtil/pressure_mod.F). The dry surface pressure
         ! that is used is calculated from the GMAO wet surface pressure 
         ! (or time-interpolated value) and the A and B parameters 
         ! (see SET_DRY_SURFACE_PRESSURE). It is not derived from the
         ! wet edge pressures. This distinction is important for
         ! mass conservation.
         State_Met%AD(I,J,L) = ( State_Met%DELP_DRY(I,J,L) * G0_100 )
     &                       * AREA_M2

         !==============================================================
         ! Set grid box dry air partial pressures at grid box 
         ! altitude-weighted mean pressure [hPa]
         ! Assume constant humidity across grid box.
         !==============================================================

         ! Mean altitude-weighted pressures in grid box [hPa] defined as   
         ! average P(z) over z. Use in the ideal gas law yields a         
         ! mean density equivalent to total mass per volume in grid box.  
         PMEAN     = State_Met%DELP(I,J,L) /
     &                              log( State_Met%PEDGE(I,J,L) /
     &                                   PEdge_Top )
         Ev_mean   = PMEAN * Ev_mid / State_Met%PMID(I,J,L)
         PMEAN_DRY = PMEAN - Ev_mean

         ! NOTE: Try the below definition in the future to change the
         ! AIRDEN equation to use the ideal gas law, thereby removing 
         ! area-independence of AIRDEN
         !State_Met%PMEAN_DRY( I,J,L ) = State_Met%DELP_DRY(I,J,L) /
         !&                              log( State_Met%PEDGE(I,J,L) /         
         !&                                   PEdge_Top )     

         !==============================================================
         ! Set grid box densities
         !============================================================== 
         
         ! Set grid box dry air density [kg/m3]
         !
         ! NOTE: Air density is derived from dry surface pressure
         ! following implementation of the moisture fix, and therefore
         ! its equation must be dry air mass divided by volume. 
         ! This is because the level pressures derived from the dry 
         ! surface pressure may be used for extracting mass per level, 
         ! but not a representative pressure for that level. The ideal 
         ! gas law requires a representative level pressure. Eventually 
         ! a representative pressure must be derived that, when used in
         ! the ideal gas law, will replicate the air density defined as
         ! mass/volume. Moving to use of the ideal gas law is necessary 
         ! for grid-independence (ewl, 9/16/16)
         State_Met%AIRDEN(I,J,L)    = State_Met%AD(I,J,L)
     &                              / State_Met%AIRVOL(I,J,L)

         ! Set grid box dry air number density [molec/cm3]
         State_Met%AIRNUMDEN(I,J,L) = State_Met%AIRDEN(I,J,L)
     &                              * 1e-3_fp *  AVO / AIRMW

         ! Set grid box moist air density [kg/m3] using the ideal gas law
         ! and pressures derived from GMAO pressure
         !
         !  MAIRDEN = density of moist air [kg/m^3],
         !  given by:
         !           
         !            Partial      Molec     Partial       Molec
         !            pressure   * wt of   + pressure of * wt of 
         !            of dry air   air       water vapor   water     
         ! Moist      [hPa]        [g/mol]   [hPa]         [g/mol]  
         ! Air     =  ------------------------------------------------        
         ! Density    Universal gas constant * Temp * 1000  * 0.01
         !                   [J/K/mol]         [K]   [g/kg]   [hPa/Pa]
         !
         ! NOTE: MAIRDEN is used in wetscav_mod only
         State_Met%MAIRDEN(I,J,L) = ( PMEAN_DRY 
     &                                * AIRMW + Ev_mean * H2OMW ) 
     &                              * PCONV * MCONV
     &                              / RSTARG / State_Met%T(I,J,L)


         !==============================================================
         ! Define the various query fields of State_Met
         !
         ! NOTE: For convenience, we set State_Met%InPbl in routine 
         ! COMPUTE_PBL_HEIGHT (in module GeosCore/pbl_mix_mod.F).
         !
         ! NOTE: For certain queries we test against level numbers,
         ! (e.g. LLSTRAT, LLCHEM), but should really test level 
         ! pressure edges, so that this algorithm will be robust if 
         ! we switch to different met fields or interface with 
         ! different ESMs.  Add this at a later time. (bmy, 1/8/18)
         !============================================================== 
 
         ! Is this grid box within the troposphere?
         State_Met%InTroposphere(I,J,L) = 
     &      ( State_Met%PEDGE(I,J,L) > State_Met%TROPP(I,J) )

         ! Is this grid box within the stratosphere or mesosphere?
         State_Met%InStratMeso(I,J,L) = 
     &      ( .not. State_Met%InTroposphere(I,J,L) )

         ! Is this grid box within the stratosphere (but not mesosphere)?
         State_Met%InStratosphere(I,J,L) =
     &      ( L <= LLSTRAT .and. State_Met%InStratMeso(I,J,L) )

         ! Is grid box (I,J,L) within the chemistry grid?
         IF ( L > LLCHEM ) THEN

            ! Chemistry is not done higher than the mesopause
            State_Met%InChemGrid(I,J,L) = .FALSE.

         ELSE

            IF ( Input_Opt%LUCX ) THEN

               ! UCX mechanisms: Chemistry grid goes up to stratopause
               State_Met%InChemGrid(I,J,L) = ( L <= LLSTRAT )

            ELSE

               ! "tropchem" mechanisms: Chemistry grid goes up to tropopause
               State_Met%InChemGrid(I,J,L) =
     &              State_Met%InTroposphere(I,J,L)

            ENDIF

         ENDIF

      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      !=================================================================
      ! Compute more tropopause and chemistry grid quantities.  This
      ! has to be done after the State_Met%InChemGrid etc. fields have 
      ! been totally initialized.
      !
      ! Also compute if it is near local noon in a grid box.
      ! This will be used for the J-value diagnostics.
      !=================================================================
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED                            )
!$OMP+PRIVATE( I, J, L_CG, L_TP, H, Pb, Pt, FRAC )
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         !--------------------------------------------------------------
         ! Local solar time and box that is nearest to local noon
         !--------------------------------------------------------------
         State_Met%LocalSolarTime(I,J) = LocTime(I)
         State_Met%IsLocalNoon(I,J)    = IsLocNoon(I)

         !--------------------------------------------------------------
         ! Compute highest level of chemistry grid and the troposphere
         !
         ! NOTE: The COUNT intrinsic function returns the number of 
         ! locations where the column State_Met%InChemGrid(I,J,:) is 
         ! TRUE.  This is also equivalent to the maximum level of the 
         ! chemistry grid.  Ditto for State_Met%InTroposphere.
         !--------------------------------------------------------------

         ! Highest level in the column at (I,J)
         ! that is fully within the chemistry grid
         L_CG  = COUNT( State_Met%InChemGrid(I,J,:) )

         ! Highest level in the column at (I,J)
         ! that is fully within the tropopause
         L_TP  = COUNT( State_Met%InTroposphere(I,J,:) )

         !--------------------------------------------------------------
         ! Compute tropopause height
         !--------------------------------------------------------------

         ! Get height (from surface to top edge) of all boxes that lie
         ! totally w/in the troposphere.  NOTE: Grid box (I,J,L_TP-1)
         ! is the highest purely tropospheric grid box in the column.
         !
         ! BMY COMMENT: L_TP may be the highest purely tropospheric
         ! grid box in the column, not L_TP-1. That might have been
         ! due to historical baggage.  In any case, this is only used
         ! to compute State_Met%TropHt which is a purely diagnostic
         ! output, so it won't make a big difference even if it is
         ! incorrect.  Check into this later on. (bmy, 1/17/18)
         H     = SUM( State_Met%BXHEIGHT(I,J,1:L_TP-1) )

         ! Get the pressures [hPa] at the bottom and top edges
         ! of the grid box in which the tropopause occurs
         Pb    = State_Met%PEDGE(I,J,L_TP  )  
         Pt    = State_Met%PEDGE(I,J,L_TP+1)

         ! FRAC is the fraction of the grid box (I,J,L_TP) 
         ! that lies totally within the troposphere
         FRAC  = ( Pb - State_Met%TROPP(I,J) ) / ( Pb - Pt ) 

         ! Add to H the height [m] of the purely tropospheric 
         ! fraction of grid box (I,J,L_TP)
         H     = H + ( FRAC * State_Met%BXHEIGHT(I,J,L_TP) )

         !--------------------------------------------------------------
         ! Save in the relevant fields of State_Met
         !--------------------------------------------------------------
         State_Met%ChemGridLev(I,J) = L_CG           ! Max chemgrid level
         State_Met%TropLev    (I,J) = L_TP           ! Max tropopause level
         State_Met%TropHt     (I,J) = H / 1000.0_fp  ! Troposphere ht [km]

      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      !=================================================================
      ! Update species concentrations with updated air quantities
      ! if update mixing ratio flag indicates to do so
      !=================================================================

      ! Update the mixing ratio with new air quantities such that species
      ! mass is conserved if (1) update_mixing_ratio flag is not passed, 
      ! or (2) update_mixing_ratio flag is passed as true.
      ! NOTE: The only places where mixing ratio is not currently updated
      ! following air quantity change is during GEOS-Chem initialization and  
      ! in transport after the pressure fixer is applied 
      IF ( UpdtMR ) THEN
!      IF ( .not. PRESENT( update_mixing_ratio )  
!     &               .or. update_mixing_ratio ) THEN

         ! The concentration update formula works only for dry mixing ratios
         ! (kg/kg or v/v) so check if units are correct
         IF ( State_Chm%Spc_units == 'kg/kg dry' .or.
     &        State_Chm%Spc_units == 'v/v dry' ) THEN

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, N )
            DO N = 1, State_Chm%nSpecies

               DO L = 1, LLPAR
               DO J = 1, JJPAR
               DO I = 1, IIPAR

               State_Chm%Species(I,J,L,N) = State_Chm%Species(I,J,L,N)
     &                                    * State_Met%DP_DRY_PREV(I,J,L)
     &                                    / State_Met%DELP_DRY(I,J,L) 
               ENDDO
               ENDDO
               ENDDO
               ENDDO
!$OMP END PARALLEL DO

         ELSE
            ErrMsg = 'Incorrect species units: ' // 
     &                TRIM( State_Chm%Spc_Units )
            CALL GC_Error( ErrMsg, RC, ThisLoc )
            RETURN
         ENDIF
      ENDIF

      END SUBROUTINE AIRQNT
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: interp
!
! !DESCRIPTION: Subroutine INTERP linearly interpolates GEOS-Chem I6 and I3
!  fields (winds, surface pressure, temperature, surface albedo, specific 
!  humidity etc.)  to the current dynamic timestep.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE INTERP( NTIME0, NTIME1, NTDT, Input_Opt, State_Met )
!
! !USES:
!
      USE GC_GRID_MOD,        ONLY : GET_YEDGE
      USE Input_Opt_Mod,      ONLY : OptInput
      USE State_Met_Mod,      ONLY : MetState
!
! !INPUT PARAMETERS: 
!
      INTEGER,        INTENT(IN)    :: NTIME0     ! Elapsed time [s] at 
                                                  !  start of outer time step
      INTEGER,        INTENT(IN)    :: NTIME1     ! Elapsed time [s] at 
                                                  !  current time
      INTEGER,        INTENT(IN)    :: NTDT       ! Dynamic timestep [s]
      TYPE(OptInput), INTENT(IN)    :: Input_Opt  ! Input Options object
!
! !INPUT/OUTPUT PARAMETERS: 
!
      TYPE(MetState), INTENT(INOUT) :: State_Met  ! Meteorology State object
! 
! !REMARKS:
!  Different met fields are archived at I6 (instantaneous 6-hr) and
!  I3 (instantaneous 3-hr) time resolution depending on the specific product.  
!  For example, relative humidity is an instantaneous 6hr field in MERRA and 
!  a 6-hr time averaged field in GEOS-5.
!
! !REVISION HISTORY: 
!  30 Jan 1998 - R. Yantosca - Initial version
!  (1 ) INTERP is written in Fixed-Form Fortran 90.
!  (2 ) Subtract PINT from PSC since the only subroutine that uses PSC
!        is TPCORE.  This prevents having to subtract and add PINT to PSC
!        before and after each call of TPCORE.
!  (3 ) Pass the Harvard CTM temperature variable T(IGCMPAR,JGCMPAR,LGCMPAR)
!        to INTERP via the argument list (instead of including file CMN).
!        It is computationally inefficient to keep two large arrays for
!        the same quantity.  Use the proper window offsets with T.
!  (4 ) Added to "dao_mod.f" (bmy, 6/26/00)
!  (5 ) Updated comments (bmy, 4/4/01)
!  (6 ) Replaced {IJL}GCMPAR w/ IIPAR,JJPAR,LLPAR.  Also now use parallel
!        DO-loop for interpolation.  Updated comments. (bmy, 9/26/01)
!  (7 ) Removed obsolete code from 9/01 (bmy, 10/23/01)
!  (8 ) Add PSC2 as the surface pressure at the end of the dynamic timestep.
!        This needs to be passed to TPCORE and AIRQNT so that the mixing ratio
!        can be converted to mass properly.  Removed PINT from the arg list,
!        since we don't need it anymore.  Also updated comments and made
!        some cosmetic changes.  (bmy, 3/27/02)
!  (9 ) Removed obsolete, commented-out code (bmy, 6/25/02)
!  (10) Eliminated PS, PSC from the arg list, for floating-pressure fix.
!        (dsa, bdf, bmy, 8/27/02)
!  (11) Met field arrays are module variables, so we don't need to pass them
!        as arguments. (bmy, 11/20/02)
!  (12) Removed NDT from the arg list since that is always 21600.  For GEOS-4
!        met fields, only interpolate PSC2; the other fields are 6-h averages.
!        Eliminate TC variable, it's obsolete.  Now use double precision to
!        compute TM and TC2 values.  Renamed NTIME to NTIME1 and NTIME1 to
!        NTIME0.  Updated comments. (bmy, 6/19/03)
!  (13) Now modified for GEOS-5 and GCAP met fields. (swu, bmy, 5/25/05)
!  (14) Remove support for GEOS-1 and GEOS-STRAT met fields (bmy, 8/4/06)
!  (15) Now interpolate TROPP, only if variable tropopause is used 
!        (phs, 9/12/06)
!  (16) Don't interpolate TROPP for GEOS-5 (bmy, 1/17/07)
!  (17) Now limit tropopause pressure to 200 mbar at latitudes above 60deg
!        (phs, 9/18/07)
!  16 Aug 2010 - R. Yantosca - Added ProTeX headers
!  18 Aug 2010 - R. Yantosca - Rewrite #if block logic for clarity
!  06 Feb 2012 - R. Yantosca - Add modifications for GEOS-5.7.x met fields
!  28 Feb 2012 - R. Yantosca - Removed support for GEOS-3
!  01 Mar 2012 - R. Yantosca - Now use GET_YEDGE(I,J,L) from new grid_mod.F90
!  09 Nov 2012 - M. Payer    - Replaced all met field arrays with State_Met
!                              derived type object
!  26 Sep 2013 - R. Yantosca - Renamed GEOS_57 Cpp switch to GEOS_FP
!  29 Oct 2013 - R. Yantosca - Now interpolate T_FULLGRID field for GEOS-FP met
!  12 May 2015 - E. Lundgren - Fix I3 interpolation bug where time fraction
!                              spans from 1 to 2, not 0 to 1, for hrs 3 to 6
!  11 Aug 2015 - R. Yantosca - MERRA2 behaves the same way as GEOS-FP
!  28 Oct 2015 - E. Lundgren - Set previous SPHU to current SPHU before update
!  04 May 2016 - E. Lundgren - Now also interpolate dry pressure calculated
!                              in SET_DRY_SURFACE_PRESSURE using GMAO wet P
!  07 Nov 2017 - R. Yantosca - Dynamic tropopause is now always assumed
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER   :: I,        J,        L
      REAL(fp)  :: D_NTIME0, D_NTIME1, D_NDT
      REAL(fp)  :: D_NTDT,   TM,       TC2
      REAL(fp)  :: YSOUTH,   YNORTH

      !=================================================================
      ! Initialization
      !=================================================================

      ! Convert time variables from FLOAT to DBLE
      D_NTIME0 = DBLE( NTIME0 )
      D_NTIME1 = DBLE( NTIME1 )
      D_NTDT   = DBLE( NTDT   )

      D_NDT    = 10800e+0_fp              ! For 3-hr instantaneous fields

      ! Fraction of timestep elapsed at mid point of this dyn timestep
      TM  = ( D_NTIME1 + D_NTDT/2e+0_fp - D_NTIME0 ) / D_NDT
      
      ! Fraction of timestep elapsed at the end of this dyn timestep
      TC2 = ( D_NTIME1 + D_NTDT     - D_NTIME0 ) / D_NDT 

      ! For I3 fields, need to reset fraction after 3 hours (ewl, 5/12/2015)
      IF ( TM > 1.0e+0_fp ) THEN
         TM  = TM  - 1.0e+0_fp
         TC2 = TC2 - 1.0e+0_fp
      ENDIF

      !=================================================================
      ! GEOS-FP or MERRA2 met fields: PSC2_WET, PSC2_DRY, TROPP, T, SPHU 
      !=================================================================
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, YSOUTH, YNORTH )
      DO L = 1, LLPAR
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         !----------------------------------------------------
         ! Interpolate 2D variables
         !----------------------------------------------------
         IF ( L == 1 ) THEN

            ! North & south edges of box
            YSOUTH     = GET_YEDGE( I, J,   L )
            YNORTH     = GET_YEDGE( I, J+1, L )

            ! Interpolate wet pressure [hPa] to the end of the dynamic timestep
            State_Met%PSC2_WET(I,J)  = State_Met%PS1_WET(I,J) +
     &                               ( State_Met%PS2_WET(I,J) -
     &                                 State_Met%PS1_WET(I,J) ) * TC2 
   
            ! Do the same for dry pressure [hPa] (ewl, 5/4/16)
            State_Met%PSC2_DRY(I,J)  = State_Met%PS1_DRY(I,J) +
     &                               ( State_Met%PS2_DRY(I,J) -
     &                                 State_Met%PS1_DRY(I,J) ) * TC2 

            ! Even though TROPP is a 3-hour average field, we 
            ! we still need to make sure to cap TROPP in the
            ! polar regions (if the entire box is outside 60S-60N)
            ! so that we don't do chemistry at an abnormally high
            ! altitude.  Set TROPP in the polar regions to 200 hPa.
            ! (jal, phs, bmy, 9/18/07)
            IF ( YSOUTH >= 60.0_fp .OR. YNORTH <= -60.0_fp ) THEN
               State_Met%TROPP(I,J) = MAX( State_Met%TROPP(I,J),
     &                                     200e+0_fp )
            ENDIF
         ENDIF

         !----------------------------------------------------
         ! Interpolate 3D variables
         !----------------------------------------------------

         ! Interpolate temperature
         State_Met%T   (I,J,L) = State_Met%TMPU1(I,J,L)   +
     &                         ( State_Met%TMPU2(I,J,L)   -
     &                           State_Met%TMPU1(I,J,L) ) * TM

         ! Interpolate specific humidity
         State_Met%SPHU(I,J,L) = State_Met%SPHU1(I,J,L)   +
     &                         ( State_Met%SPHU2(I,J,L)   -
     &                           State_Met%SPHU1(I,J,L) ) * TM 

      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      END SUBROUTINE INTERP
!EOC
#if defined( ESMF_ ) || defined( EXTERNAL_GRID )
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: gigc_cap_tropopause_prs
!
! !DESCRIPTION: Subroutine GIGC\_CAP\_TROPOPAUSE\_PRS caps the tropopause
!  pressure in polar latitudes to 200 hPa, so that we don't end up doing
!  troposheric chemistry too high over the poles.  This is done in the
!  standalone GEOS-Chem, and we also need to apply this when running
!  GEOS-Chem within the GEOS-5 GCM.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE GIGC_Cap_Tropopause_Prs( am_I_Root, IM,        JM, 
     &                                    Input_Opt, State_Met, RC  )
!
! !USES:
!
      USE ErrCode_Mod
      USE Input_Opt_Mod, ONLY : OptInput
      USE State_Met_Mod, ONLY : MetState
      USE GC_Grid_Mod,   ONLY : Get_XEdge
      USE GC_Grid_Mod,   ONLY : Get_YEdge
!
! !INPUT PARAMETERS:
!
      LOGICAL,        INTENT(IN)    :: am_I_Root     ! Are we on the root CPU?
      INTEGER,        INTENT(IN)    :: IM            ! # of lons on this CPU
      INTEGER,        INTENT(IN)    :: JM            ! # of lats on this CPU
      TYPE(OptInput), INTENT(IN)    :: Input_Opt     ! Input Options object
!
! !INPUT/OUTPUT PARAMETERS:
!
      TYPE(MetState), INTENT(INOUT) :: State_Met     ! Meteorology State object
!
! !OUTPUT PARAMETERS:
!
      INTEGER,        INTENT(OUT)   :: RC            ! Success or failure
!
! !REMARKS:
!  Jennifer Logan (see correspondence below) suggested that we should cap the 
!  variable tropopause at 200hPa in near-polar regions (90-60S and 60-90N), 
!  to avoid the problem with anomalously high tropopause heights at high 
!  latitudes. This fix was standardized in GEOS-Chem v7-04-13.
!                                                                             .
!  Jennifer Logan wrote:
!     I think we should restrict the tropopause at latitudes > 60 deg. to 
!     pressures greater than 200 mb (about 11 km). From Fig. 3 in Seidel and 
!     Randel, there are tropopause (TP) heights as high as 13.5 km in the 
!     Antarctic (median height is ~9.8 km, 250 mb), but I don't think we want 
!     to be doing trop. chem there. The median TP pressure at ~80 N is ~300 mb, 
!     compared to ~250 mb at 70-85 S. The extratropical TP heights are higher
!     (lower pressure) in the SH than in the NH according to Fig. 3. 
!     This approach is also very easy to explain in a paper. 
! 
! !REVISION HISTORY: 
!  14 Mar 2013 - R. Yantosca - Initial version
!  08 Mar 2018 - E. Lundgren - Move this rtn from obsolete GCHP/gchp_utils.F90
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER :: I,       J
      REAL*8  :: YSOUTH,  YNORTH
      
      ! Assume success
      RC = GC_SUCCESS
      
      ! Return if option not set
      IF ( .NOT. Input_Opt%LCAPTROP ) RETURN
      
      ! Loop over grid boxes on this PET
      DO J = 1, JM
      DO I = 1, IM
      
         ! North & south edges of box
         YSOUTH = GET_YEDGE( I, J,   1 )
         YNORTH = GET_YEDGE( I, J+1, 1 )
      
         ! Cap tropopause height at 200 hPa polewards of 60N and 60S
         IF ( YSOUTH >= 60d0 .or. YNORTH <= -60d0 ) THEN
            State_Met%TROPP(I,J) = MAX( State_Met%TROPP(I,J), 200d0 )
         ENDIF
      
      ENDDO
      ENDDO

      END SUBROUTINE GIGC_Cap_Tropopause_Prs
!EOC
#endif
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: set_dry_surface_pressure
!
! !DESCRIPTION: Subroutine SET\_DRY\_SURFACE\_PRESSURE sets the dry
!  surface pressures PS1\_DRY or PS2\_DRY by removing the water vapor from
!  the column constructed with MET pressure fields PS1\_WET or PS2\_WET. 
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE SET_DRY_SURFACE_PRESSURE( State_Met, PS_ID )
!
! !USES:
!
      USE CMN_SIZE_MOD 
      USE ERROR_MOD,            ONLY : CHECK_VALUE
      USE State_Met_Mod,        ONLY : MetState
      USE PRESSURE_MOD,         ONLY : GET_AP, GET_BP
!
! !INPUT PARAMETERS: 
!
      INTEGER,        INTENT(IN)     :: PS_ID       ! 1 = PS1, 2 = PS2
!
! !INPUT/OUTPUT PARAMETERS: 
!
      TYPE(MetState), INTENT(INOUT)  :: State_Met   ! Object for met fields
! 
! !REMARKS:
!  This subroutine is an adaptation of the GEOS-Chem moisture fix 
!  implemented by Meemong Lee (JPL) in the adjoint model. Like PS1_WET
!  and PS2_WET, from which PS1_DRY and PS2_DRY are derived, these values
!  are interpolated within routine INTERP to derive instantaneous PSC2_DRY.
!  Note that PSC2_WET and PSC2_DRY are not used to fetch pressures and 
!  calculate air quantities until after advection. Also note that the
!  dry surface pressure calculated in this routine may be used to 
!  calculate delta dry pressures across a level by utilitizing GMAO 
!  parameters Ap and Bp but should not be used to extract dry pressure 
!  edge values as a height proxy.
!
! !REVISION HISTORY:
!  03 May 2002 - E. Lundgren - Initial version (M. Lee, JPL) 
!  18 Feb 2017 - S. D. Eastham - Updated to include TOA pressure
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
! 
      INTEGER            :: I, J, L
      REAL(fp)           :: PEDGE_BOT, PEDGE_TOP

      ! Pointers
      REAL(fp), POINTER  :: PS_WET(:,:) => NULL()
      REAL(fp), POINTER  :: PS_DRY(:,:) => NULL()
      REAL(fp), POINTER  :: SPHU(:,:,:) => NULL()

      !=================================================================
      ! SET_DRY_SURFACE_PRESSURE begins here!
      !=================================================================

      ! Set pointer to the appropriate humidity
      IF ( PS_ID == 1 ) THEN
         SPHU   => State_Met%SPHU1
      ELSE IF ( PS_ID == 2 ) THEN
         SPHU   => State_Met%SPHU2
      ENDIF

      ! Set pointers to State_Met pressure fields.
      IF ( PS_ID == 1 ) THEN
         PS_WET => State_Met%PS1_WET
         PS_DRY => State_Met%PS1_DRY
      ELSE IF ( PS_ID == 2 ) THEN
         PS_WET => State_Met%PS2_WET
         PS_DRY => State_Met%PS2_DRY
      ELSE 
         ! Throw an error (to be added)
      ENDIF

      ! Reset dry surface pressure to TOA value
      PS_DRY = GET_AP(LLPAR+1)

      ! Calculate dry surface pressure from GMAO wet pressure as the 
      ! column sum of wet delta pressures with humidity removed
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, PEDGE_BOT, PEDGE_TOP )
      DO J = 1, JJPAR
      DO I = 1, IIPAR
      DO L = 1, LLPAR
         PEDGE_BOT   = GET_AP(L) + GET_BP(L) * PS_WET(I,J) 
         PEDGE_TOP   = GET_AP(L+1) + GET_BP(L+1) * PS_WET(I,J)
         PS_DRY(I,J) = PS_DRY(I,J) + ( PEDGE_BOT - PEDGE_TOP ) 
     &                     * ( 1.e+0_fp - SPHU(I,J,L) * 1.0e-3_fp )
      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      ! If dry pressure is negative, set equal to moist pressure
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J )
      DO J = 1, JJPAR
      DO I = 1, IIPAR
         IF ( PS_DRY(I,J) < 0.e+0_fp) THEN
            PS_DRY(I,J) = PS_WET(I,J)
         ENDIF
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      ! Nullify pointers
      PS_WET  => NULL()
      PS_DRY  => NULL()
      SPHU    => NULL()

      END SUBROUTINE SET_DRY_SURFACE_PRESSURE
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: is_land
!
! !DESCRIPTION: Function IS\_LAND returns TRUE if surface grid box (I,J) is 
!  a land box.
!\\
!\\
! !INTERFACE:
!
      FUNCTION IS_LAND( I, J, State_Met ) RESULT ( LAND )
!
! !USES:
!
      USE State_Met_Mod,      ONLY : MetState
      USE TIME_MOD,           ONLY : GET_YEAR
!
! !INPUT PARAMETERS: 
!
      INTEGER,        INTENT(IN)  :: I           ! Longitude index of grid box
      INTEGER,        INTENT(IN)  :: J           ! Latitude  index of grid box
      TYPE(MetState), INTENT(IN)  :: State_Met   ! Meteorology State object
!
! !RETURN VALUE:
!
      LOGICAL                     :: LAND        ! =T if it is a land box
! 
! !REVISION HISTORY: 
!  26 Jun 2000 - R. Yantosca - Initial version
!  (1 ) Now use ALBEDO field to determine land or land ice boxes for GEOS-3.
!        (bmy, 4/4/01)
!  (2 ) For 4x5 data, regridded albedo field can cause small inaccuracies
!        near the poles (bmy, 4/4/01)
!  (3 ) Add references to CMN_SIZE and CMN, so that we can use the JYEAR
!        variable to get the current year.  Also, for 1998, we need to compute
!        if is a land box or not from the surface albedo, since for this
!        year the LWI/SURFTYPE field is not given.  For other years than 1998,
!        we use LWI(I,J) < 50 as our land box criterion.  Deleted obsolete
!        code and updated comments.(mje, bmy, 1/9/02)
!  (4 ) Deleted GEOS-2 #ifdef statement.  GEOS-2 met fields never really
!        materialized, we use GEOS-3 instead. (bmy, 9/18/02)
!  (5 ) Now uses function GET_YEAR from "time_mod.f".  Removed reference
!        to CMN header file. (bmy, 3/11/03)
!  (6 ) Added code to determine land boxes for GEOS-4 (bmy, 6/18/03)
!  (7 ) Now modified for GEOS-5 and GCAP met fields (swu, bmy, 5/25/05)
!  (8 ) Now return TRUE only for land boxes (w/ no ice) (bmy, 8/10/05)
!  (9 ) Now use NINT to round LWI for GEOS-4/GEOS-5 (ltm, bmy, 5/9/06)
!  (10) Remove support for GEOS-1 and GEOS-STRAT met fields (bmy, 8/4/06)
!  16 Aug 2010 - R. Yantosca - Added ProTeX headers
!  25 Aug 2010 - R. Yantosca - Treat MERRA in the same way as GEOS-5
!  06 Feb 2012 - R. Yantosca - Treat GEOS-5.7.x in the same way as MERRA/GEOS-5
!  28 Feb 2012 - R. Yantosca - Removed support for GEOS-3
!  09 Nov 2012 - M. Payer    - Replaced all met field arrays with State_Met
!                              derived type object
!EOP
!------------------------------------------------------------------------------
!BOC

      ! LWI=1 and ALBEDO less than 69.5% is a LAND box 
      LAND = ( NINT( State_Met%LWI(I,J) ) == 1       .and.
     &               State_Met%ALBD(I,J)  <  0.695e+0_fp )

      END FUNCTION IS_LAND
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: is_water 
!
! !DESCRIPTION: Function IS\_WATER returns TRUE if surface grid box (I,J) is 
!  an ocean or an ocean-ice box.  
!\\
!\\
! !INTERFACE:
!
      FUNCTION IS_WATER( I, J, State_Met ) RESULT ( WATER )
!
! !USES:
!
      USE State_Met_Mod,      ONLY : MetState
      USE TIME_MOD,           ONLY : GET_YEAR
!
! !INPUT PARAMETERS: 
!
      INTEGER,         INTENT(IN) :: I           ! Longitude index of grid box
      INTEGER,         INTENT(IN) :: J           ! Latitude  index of grid box
      TYPE(MetState),  INTENT(IN) :: State_Met   ! Meteorology State object
!
! !RETURN VALUE:
!
      LOGICAL                     :: WATER       ! =T if this is a water box
! 
! !REVISION HISTORY: 
!  30 Jan 1998 - R. Yantosca - Initial version
!  (1 ) Now use ALBEDO field to determine water or water ice boxes for GEOS-3.
!        (bmy, 4/4/01)
!  (2 ) For 4x5 data, regridded albedo field can cause small inaccuracies
!        near the poles (bmy, 4/4/01)
!  (3 ) Add references to CMN_SIZE and CMN, so that we can use the JYEAR
!        variable to get the current year.  Also, for 1998, we need to compute
!        if is an ocean box or not from the surface albedo, since for this
!        year the LWI/SURFTYPE field is not given.  For other years than 1998,
!        we use LWI(I,J) >= 50 as our ocean box criterion.  Deleted obsolete
!        code and updated comments. (mje, bmy, 1/9/02)
!  (4 ) Deleted GEOS-2 #ifdef statement.  GEOS-2 met fields never really
!        materialized, we use GEOS-3 instead. (bmy, 9/18/02)
!  (5 ) Now uses function GET_YEAR from "time_mod.f".  Removed reference
!        to CMN header file. (bmy, 3/11/03)
!  (6 ) Added code to determine water boxes for GEOS-4 (bmy, 6/18/03)
!  (7 ) Now modified for GEOS-5 and GCAP met fields (swu, bmy, 5/25/05)
!  (8 ) Now remove test for sea ice (bmy, 8/10/05)
!  (9 ) Now use NINT to round LWI for GEOS-4/GEOS-5 (ltm, bmy, 5/9/06)
!  (10) Remove support for GEOS-1 and GEOS-STRAT met fields (bmy, 8/4/06
!  16 Aug 2010 - R. Yantosca - Added ProTeX headers
!  25 Aug 2010 - R. Yantosca - Treat MERRA in the same way as GEOS-5
!  06 Feb 2012 - R. Yantosca - Treat GEOS-5.7.x in the same way as MERRA/GEOS-5
!  28 Feb 2012 - R. Yantosca - Removed support for GEOS-3
!  09 Nov 2012 - M. Payer    - Replaced all met field arrays with State_Met
!                              derived type object
!EOP
!------------------------------------------------------------------------------
!BOC

      ! LWI=0 and ALBEDO less than 69.5% is a water box 
      WATER = ( NINT( State_Met%LWI(I,J) ) == 0       .and.
     &                State_Met%ALBD(I,J)  <  0.695e+0_fp )

      END FUNCTION IS_WATER
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: is_ice
!
! !DESCRIPTION: Function IS\_ICE returns TRUE if surface grid box (I,J) 
!  contains either land-ice or sea-ice. 
!\\
!\\
! !INTERFACE:
!
      FUNCTION IS_ICE( I, J, State_Met ) RESULT ( ICE )
!
! !USES:
!
      USE State_Met_Mod,      ONLY : MetState
      USE TIME_MOD,           ONLY : GET_YEAR
!
! !INPUT PARAMETERS: 
!
      INTEGER,         INTENT(IN) :: I           ! Longitude index of grid box
      INTEGER,         INTENT(IN) :: J           ! Latitude  index of grid box
      TYPE(MetState),  INTENT(IN) :: State_Met   ! Meteorology State object
!
! !RETURN VALUE:
!
      LOGICAL                     :: ICE         ! =T if this is an ice box
!
! 
! !REVISION HISTORY: 
!  09 Aug 2005 - R. Yantosca - Initial version
!  (1 ) Remove support for GEOS-1 and GEOS-STRAT met fields (bmy, 8/4/06)
!  16 Aug 2010 - R. Yantosca - Added ProTeX headers
!  25 Aug 2010 - R. Yantosca - Treat MERRA in the same way as GEOS-5
!  06 Feb 2012 - R. Yantosca - Treat GEOS-5.7.x in the same way as MERRA/GEOS-5
!  28 Feb 2012 - R. Yantosca - Removed support for GEOS-3
!EOP
!------------------------------------------------------------------------------
!BOC

      ! LWI=2 or ALBEDO > 69.5% is ice
      ICE = ( NINT( State_Met%LWI(I,J) ) == 2       .or.
     &              State_Met%ALBD(I,J)  >= 0.695e+0_fp )

      END FUNCTION IS_ICE
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: is_near
!
! !DESCRIPTION: Function IS\_NEAR returns TRUE if surface grid box (I,J) 
!  contains any land above a certain threshold (THRESH) or any of the 
!  adjacent boxes up to NEIGHBOR boxes away contain land.  
!\\
!\\
! !INTERFACE:
!
      FUNCTION IS_NEAR( I, J, THRESH, NEIGHBOR, State_Met )
     &   RESULT ( NEAR )
!
! !USES:
!
      USE State_Met_Mod,      ONLY : MetState
!
! !INPUT PARAMETERS: 
!
      ! Arguments
      INTEGER,        INTENT(IN) :: I, J       ! Lon & lat grid box indices
      INTEGER,        INTENT(IN) :: NEIGHBOR   ! # of neighbor boxes to consider
      REAL(fp),         INTENT(IN) :: THRESH     ! LWI threshold for near-land 
      TYPE(MetState), INTENT(IN) :: State_Met  ! Meteorology State object
!
! !RETURN VALUE:
!
      LOGICAL                    :: NEAR       ! # of near land boxes
!
! !REMARKS:
!  Typical values for:
!     GCAP   : THRESH =  0.2, NEIGHBOR = 1
!     GEOS-3 : THRESH = 80.0, NEIGHBOR = 1
!     GEOS-4 : THRESH =  0.2, NEIGHBOR = 1
!     GEOS-5 : THRESH =  0.2, NEIGHBOR = 1
!                                                                             .
!  NOTE: This routine is mostly obsolete now.
! 
! !REVISION HISTORY: 
!  09 May 2006 - R. Yantosca - Initial version
!  (1 ) Modified for GCAP and GEOS-3 met fields (bmy, 5/16/06)
!  (2 ) Remove support for GEOS-1 and GEOS-STRAT met fields (bmy, 8/4/06)
!  16 Aug 2010 - R. Yantosca - Added ProTeX headers
!  19 Aug 2010 - R. Yantosca - Rewrote logic of #if block for clarity
!  25 Aug 2010 - R. Yantosca - Treat MERRA in same way as GEOS-5
!  06 Feb 2012 - R. Yantosca - Treat GEOS-5.7.x in the same way as MERRA/GEOS-5
!  28 Feb 2012 - R. Yantosca - Removed support for GEOS-3
!  09 Nov 2012 - M. Payer    - Replaced all met field arrays with State_Met
!                              derived type object
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER :: NS, EW, LONGI, LATJ

      !=================================================================
      ! IS_NEAR begins here!
      !=================================================================

      ! Initialize
      NEAR = .FALSE.

      ! Loop over neighbor lat positions
      DO NS = -NEIGHBOR, NEIGHBOR

         ! Lat index of neighbor box
         LATJ = J + NS

         ! Special handling near poles
         IF ( LATJ < 1 .or. LATJ > JJPAR ) CYCLE

         ! Loop over neighbor lon positions
         DO EW = -NEIGHBOR, NEIGHBOR

            ! Lon index of neighbor box
            LONGI = I + EW

            ! Special handling near date line
            IF ( LONGI < 1     ) LONGI = LONGI + IIPAR 
            IF ( LONGI > IIPAR ) LONGI = LONGI - IIPAR
            
            ! If it's an ice box, skip to next neighbor
            IF ( IS_ICE( LONGI, LATJ, State_Met ) ) CYCLE

            !---------------------------------------------------
            ! LWI = 0.0 is ocean
            ! LWI = 1.0 is land
            ! LWI = 2.0 is ice 
            !
            ! with fractional values at land-water, land-ice,
            ! and water-ice boundaries.
            !
            ! It's near-land if THRESH <= LWI <= 1.0 
            !---------------------------------------------------
            IF ( State_Met%LWI(LONGI,LATJ) >  THRESH  .and.
     &           State_Met%LWI(LONGI,LATJ) <= 1e+0_fp    ) THEN

               ! We are in a near-land box
               NEAR = .TRUE.

               ! Break out of loop
               GOTO 999
            ENDIF

         ENDDO
      ENDDO

      ! Exit
 999  CONTINUE

      END FUNCTION IS_NEAR
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_obk
!
! !DESCRIPTION: Function GET\_OBK returns the Monin-Obhukov length at a grid 
!  box (I,J).
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_OBK( I, J, State_Met ) RESULT( OBK )
!
! !USES:
!
      USE State_Met_Mod,      ONLY : MetState
!
! !INPUT PARAMETERS: 
!
      INTEGER,         INTENT(IN) :: I           ! Longitude index
      INTEGER,         INTENT(IN) :: J           ! Latitude  index
      TYPE(MetState),  INTENT(IN) :: State_Met   ! Meteorology State object
!
! !RETURN VALUE:
!
      REAL(fp)                      :: OBK         ! Monin-Obhukhov length
!
! !REMARKS:
! 
! 
! !REVISION HISTORY: 
!  25 May 2005 - R. Yantosca - Initial version
!  16 Aug 2010 - R. Yantosca - Added ProTeX headers
!  09 Nov 2012 - M. Payer    - Replaced all met field arrays with State_Met
!                              derived type object
!EOP
!------------------------------------------------------------------------------
!BOC
!
      !=================================================================
      ! For all GEOS met fields:
      !
      ! The direct computation of the Monin-Obhukov length is:
      !
      !            - Air density * Cp * T(surface air) * Ustar^3 
      !    OBK =  -----------------------------------------------
      !              Kappa       * g  * Sensible Heat flux
      !
      ! Cp    = 1000 J / kg / K = specific heat of air at constant P
      ! Kappa = 0.4             = Von Karman's constant
      !
      !
      !  Also test the denominator in order to prevent div by zero.
      !=================================================================

      ! Local variables
      REAL(fp)            :: NUM, DEN

      ! Parameters
      REAL(fp), PARAMETER :: KAPPA = 0.4e+0_fp 
      REAL(fp), PARAMETER :: CP    = 1000.0e+0_fp

      ! Numerator
      NUM = - State_Met%AIRDEN(I,J,1) * CP                   *
     &        State_Met%TS(I,J)       * State_Met%USTAR(I,J) *
     &        State_Met%USTAR(I,J)    * State_Met%USTAR(I,J)

      ! Denominator
      DEN =  KAPPA         * g0       * State_Met%HFLUX(I,J) 

      ! Prevent div by zero
      IF ( ABS( DEN ) > 0e+0_fp ) THEN
         OBK = NUM / DEN
      ELSE
         OBK = 1.0e+5_fp
      ENDIF

      END FUNCTION GET_OBK
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_cosine_sza
!
! !DESCRIPTION: Routine GET\_COSINE\_SZA is a driver for calling the COSSZA 
!  routine from dao\_mod.F.  This routine calls COSSZA twice.  The first call
!  computes the sun angles at the current time and midpoint of the current 
!  chemistry time step.  The second call computes the sun angles 5 hours 
!  prior to the current time (for the PARANOX ship emissions plume model).
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE GET_COSINE_SZA( am_I_Root, Input_Opt, State_Met, RC  )
!
! USES:
!
      USE ErrCode_Mod
      USE Input_Opt_Mod,      ONLY : OptInput
      USE JULDAY_MOD,         ONLY : JULDAY
      USE State_Met_Mod,      ONLY : MetState
      USE TIME_MOD,           ONLY : GET_DAY_OF_YEAR
      USE TIME_MOD,           ONLY : GET_DAY
      USE TIME_MOD,           ONLY : GET_GMT
      USE TIME_MOD,           ONLY : GET_HOUR 
      USE TIME_MOD,           ONLY : GET_MINUTE
      USE TIME_MOD,           ONLY : GET_SECOND
      USE TIME_MOD,           ONLY : GET_MONTH
      USE TIME_MOD,           ONLY : GET_YEAR
!
! !INPUT PARAMETERS:
!
      LOGICAL,        INTENT(IN)    :: am_I_Root   ! Are we on the root CPU?
      TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
!
! !INPUT/OUTPUT PARAMETERS:
!
      TYPE(MetState), INTENT(INOUT) :: State_Met   ! Meteorology State object
!
! !OUTPUT PARAMETERS:
!
      INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
!
! !REVISION HISTORY: 
!  07 Feb 2012 - R. Yantosca - Initial version
!  27 Nov 2012 - R. Yantosca - Add am_I_root, Input_Opt, State_Met, RC args
!  27 Nov 2012 - R. Yantosca - Now pass State_Met to COSSZA so that the
!                              SUNCOS fields may be updated
!  28 Nov 2012 - R. Yantosca - Removed references to 1-D SUNCOS arrays
!  04 Dec 2014 - M. Yannetti - REAL*8 needs to stay *8 due to JULDAY
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARAIBLES
!
      ! Scalars
      INTEGER :: DOY,  DOY5, YEAR, MONTH,  DAY, HOUR, HOUR5, MINUTE
      INTEGER :: SECOND
      REAL*8  :: DDAY, JD,   JD5,  JDJan0, GMT, GMT5

      !=================================================================
      ! Get cosine(SZA) at the current time (SUNCOS) and at the
      ! midpoint of the chemistry timestep (SUNCOS_MID)
      !=================================================================

      ! Assume success
      RC     = GC_SUCCESS

      ! Current time
      DOY    = GET_DAY_OF_YEAR()                       ! Current day of year
      YEAR   = GET_YEAR()                              ! Current year
      MONTH  = GET_MONTH()                             ! Current month
      DAY    = GET_DAY()                               ! Current day of month
      HOUR   = GET_HOUR()                              ! Current GMT hour
      MINUTE = GET_MINUTE()                            ! Current GMT minutes
      SECOND = GET_SECOND()                            ! Current GMT seconds
      GMT    = GET_GMT()                               ! Current GMT
      DDAY   = DAY + ( HOUR/24d0     ) + ( MINUTE/1440d0 ) 
     &             + ( SECOND/3600d0 )                 ! Current decimal day 
! Current decimal day
      JD     = JULDAY( YEAR, MONTH, DDAY )             ! Current Julian date

      ! Compute cosine(SZA) quantities for the current time
      CALL COSSZA( DOY, HOUR,  .FALSE., State_Met )

! The following is now obsolete. SZA 5 hours ago is calculated in HEMCO
! PARANOx extension (ckeller, 5/22/15).
!
!      !=================================================================
!      ! Get cosine (SZA) at 5 hours behind the current time (SUNCOS_5hr)
!      ! and at the midpt of the chemistry timestep 5h ago (SUNCOS_MID_5hr)
!      !=================================================================
!
!      ! Time 5h ago
!      JDJan0 = JULDAY( YEAR, 1, 0d0 )                  ! Julian date on Jan 0
!      JD5    = JD  - ( 5d0 / 24d0 )                    ! Julian date 5h ago
!      DOY5   = JD5 - JDJan0                            ! Day of year 5h ago
!      GMT5   = GMT - 5d0                               ! GMT 5h ago
!      HOUR5  = INT( GMT5 )
!
!      ! Compute cosine(SZA) quantities for 5h ago
!      CALL COSSZA( DOY5, HOUR5, .TRUE., State_Met )

      END SUBROUTINE GET_COSINE_SZA
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: cossza
!
! !DESCRIPTION: COSSZA computes the cosine of the solar zenith angle, given
!  the day of the year and GMT hour.  The cosine of the solar zenith
!  angle is returned at both the current time and at the midpoint of the
!  chemistry timestep (i.e. for the centralized chemistry timestep option).
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE COSSZA( DOY, GMT_HOUR, DO_5hr_AGO, State_Met )
!
! !USES:
!
      USE GC_GRID_MOD,        ONLY : GET_YMID_R
      USE State_Met_Mod,      ONLY : MetState
      USE TIME_MOD,           ONLY : GET_MINUTE, GET_SECOND
      USE TIME_MOD,           ONLY : GET_LOCALTIME
      USE TIME_MOD,           ONLY : GET_TS_CHEM
!
! !INPUT PARAMETERS: 
!
      INTEGER,        INTENT(IN)    :: DOY          ! Day of the year
      INTEGER,        INTENT(IN)    :: GMT_HOUR     ! Hour of day
      LOGICAL,        INTENT(IN)    :: DO_5hr_AGO   ! Compute 5h ago?
!
! !INPUT/OUTPUT PARAMETERS:
!
      TYPE(MetState), INTENT(INOUT) :: State_Met    ! Meteorology State
!
! !REMARKS:
!  Hour angle (AHR) is a function of longitude.  AHR is zero at solar noon, 
!  and increases by 15 deg for every hour before or after solar noon.  Hour 
!  angle can be thought of as the time in hours since the sun last passed
!  the meridian (i.e. the time since the last local noon).
!                                                                             .
!  The cosine of the solar zenith angle (SZA) is given by:
!                                                                             .
!     cos(SZA) = sin(LAT)*sin(DEC) + cos(LAT)*cos(DEC)*cos(AHR) 
!                                                                             .
!  where LAT = the latitude angle, 
!        DEC = the solar declination angle,  
!        AHR = the hour angle, all in radians. 
!                                                                             .
!  If SUNCOS < 0, then the sun is below the horizon, and therefore does not 
!  contribute to any solar heating.  
!
! !REVISION HISTORY: 
!  21 Jan 1998 - R. Yantosca - Initial version
!  (1 ) COSSZA is written in Fixed-Form Fortran 90.
!  (2 ) Use IMPLICIT NONE
!  (3 ) Use C-preprocessor #include statement to include CMN_SIZE, which 
!        has IIPAR, JJPAR, LLPAR, IIPAR, JJPAR, LGLOB. 
!  (4 ) Use IM and JM (in CMN_SIZE) as loop limits.
!  (5 ) Include Harvard CTM common blocks and rename variables where needed.  
!  (6 ) Use SUNCOS(MAXIJ) instead of a 2D array, in order for compatibility
!        with the Harvard CTM subroutines.  SUNCOS loops over J, then I.
!  (7 ) Added DO WHILE loops to reduce TIMLOC into the range 0h - 24h.
!  (8 ) Cosmetic changes.  Also use F90 declaration statements (bmy, 6/5/00)
!  (9 ) Added to "dao_mod.f".  Also updated comments. (bmy, 9/27/01)
!  (10) Replaced all instances of IM with IIPAR and JM with JJPAR, in order
!        to prevent namespace confusion for the new TPCORE (bmy, 6/25/02)
!  (11) Deleted obsolete code from 6/02 (bmy, 8/21/02)
!  (12) Removed RLAT and XLON from the arg list.  Now compute these using 
!        functions from "grid_mod.f" (bmy, 2/3/03)
!  (13) Now uses GET_LOCALTIME from "time_mod.f" to get the local time. 
!        Added parallel DO loop. Removed NHMSb, NSEC arguments. (bmy, 2/13/07)
!  (14) Now compute SUNCOS at the midpoint of the relevant time interval
!        (i.e. the chemistry timestep).   Also make the A and B coefficients
!        parameters instead of variables. (bmy, 4/27/10)
!  16 Aug 2010 - R. Yantosca - Added ProTeX headers
!  05 Oct 2011 - R. Yantosca - Now also return the cosine of the solar 
!                              zenith angle at 30m after the GMT hour.
!  07 Oct 2011 - R. Yantosca - Now return SUNCOS_MID, the cos(SZA) at the
!                              midpt of the chem step (not always at 00:30).
!  07 Feb 2012 - R. Yantosca - Now add GMT_HOUR as a new argument, which !
!                              will facilitate computing sun angles 5h ago
!  01 Mar 2012 - R. Yantosca - Now use GET_YMID_R(I,J,L) from grid_mod.F90
!  01 Mar 2012 - R. Yantosca - Now use GET_LOCALTIME(I,J,L) from time_mod.F90
!  27 Nov 2012 - R. Yantosca - Update SUNCOS fields of the State_Met object
!  10 Mar 2017 - E. Lundgren - Bug fix: include partial hour in current time
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!      
      INTEGER              :: I,        J
      INTEGER              :: SECOND,   MINUTE,   TS_SUN,   FACTOR
      REAL(fp)             :: TIMLOC,   DEC,      C_DEC
      REAL(fp)             :: S_DEC,    R,        YMID_R
      REAL(fp)             :: C_YMID_R, S_YMID_R
      REAL(fp)             :: AHR,      SUNCOS,   SUNCOS_MID
      REAL(f8)             :: GMT_CUR,  GMT_MID 
!
! !DEFINED PARAMETERS:
!
      ! Coefficients for solar declination angle
      REAL(fp),  PARAMETER :: A0 = 0.006918e+0_fp
      REAL(fp),  PARAMETER :: A1 = 0.399912e+0_fp
      REAL(fp),  PARAMETER :: A2 = 0.006758e+0_fp
      REAL(fp),  PARAMETER :: A3 = 0.002697e+0_fp
      REAL(fp),  PARAMETER :: B1 = 0.070257e+0_fp
      REAL(fp),  PARAMETER :: B2 = 0.000907e+0_fp
      REAL(fp),  PARAMETER :: B3 = 0.000148e+0_fp

      !=================================================================
      ! Initialization   
      !=================================================================

      ! Quantities for central chemistry timestep
      TS_SUN   = GET_TS_CHEM()                         ! Chemistry interval
      SECOND   = GET_SECOND()                          ! Current seconds
      MINUTE   = GET_MINUTE()                          ! Current minutes
      FACTOR   = ( MINUTE * 60 + SECOND ) / TS_SUN     ! Multiplying factor

      ! GMT at the current time
      GMT_CUR  = DBLE( GMT_HOUR          )
     &         + ( DBLE( TS_SUN * FACTOR ) / 3600e+0_fp ) 

      ! GMT at the midpoint of the chemistry time interval
      GMT_MID  = ( DBLE( GMT_HOUR        )        )  
     &         + ( DBLE( TS_SUN * FACTOR ) / 3600e+0_fp ) 
     &         + ( DBLE( TS_SUN / 2      ) / 3600e+0_fp ) 

      ! Path length of earth's orbit traversed since Jan 1 [radians]
      R        = ( 2e+0_fp * PI / 365e+0_fp ) * DBLE( DOY - 1 ) 

      ! Solar declination angle (low precision formula) [radians]
      DEC      = A0 - A1*COS(         R ) + B1*SIN(         R )
     &              - A2*COS( 2e+0_fp*R ) + B2*SIN( 2e+0_fp*R )
     &              - A3*COS( 3e+0_fp*R ) + B3*SIN( 3e+0_fp*R )

      ! Pre-compute sin & cos of DEC outside of DO loops (for efficiency)
      S_DEC    = SIN( DEC )
      C_DEC    = COS( DEC )

      !=================================================================
      ! Compute cosine of solar zenith angle
      !=================================================================
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I,      J,   YMID_R, S_YMID_R,  C_YMID_R )
!$OMP+PRIVATE( TIMLOC, AHR, SUNCOS, SUNCOS_MID          )
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! Latitude of grid box [radians]
         YMID_R     = GET_YMID_R( I, J, 1 )

         ! Pre-compute sin & cos of DEC outside of I loop (for efficiency)
         S_YMID_R   = SIN( YMID_R )
         C_YMID_R   = COS( YMID_R )

         !==============================================================
         ! Compute cosine of SZA at the current GMT time
         !==============================================================

         ! Local time at box (I,J) [hours]
         TIMLOC     = GET_LOCALTIME( I, J, 1, GMT=GMT_CUR )

         ! Hour angle at box (I,J) [radians]
         AHR        = ABS( TIMLOC - 12e+0_fp ) * 15e+0_fp * PI_180
         
         ! Cosine of solar zenith angle at box (I,J) [unitless]
         SUNCOS     = ( S_YMID_R * S_DEC              ) 
     &              + ( C_YMID_R * C_DEC * COS( AHR ) )

         !==============================================================
         ! Compute cosine of SZA at the midpoint of the chem timestep
         ! Required for photolysis, chemistry, emissions, drydep
         !==============================================================

         ! Local time [hours] at box (I,J) at the midpt of the chem timestep
         TIMLOC     = GET_LOCALTIME( I, J, 1, GMT=GMT_MID )

         ! Hour angle at box (I,J) [radians]
         AHR        = ABS( TIMLOC - 12e+0_fp ) * 15e+0_fp * PI_180
         
         ! Corresponding cosine( SZA ) at box (I,J) [unitless]
         SUNCOS_MID = ( S_YMID_R * S_DEC              )
     &              + ( C_YMID_R * C_DEC * COS( AHR ) )

         !==============================================================
         ! Copy data into fields of the Meteorology State object
         !==============================================================
         IF ( DO_5hr_AGO ) THEN

!            ! COS(SZA) at midpt of chemistry timestep that was 5h ago
!            ! This is needed for the PARANOX ship plume model
!            State_Met%SUNCOSmid5(I,J) = SUNCOS_MID

         ELSE

            ! COS(SZA) at the current time
            State_Met%SUNCOS    (I,J) = SUNCOS
                               
            ! COS(SZA) @ the midpoint of the current chemistry timestep
            State_Met%SUNCOSmid (I,J) = SUNCOS_MID
                                                           
         ENDIF                                             

      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      END SUBROUTINE COSSZA
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: set_met_ageofair
!
! !DESCRIPTION: Subroutine Set\_Met\_AgeOfAir adds the time step (in seconds)
!  to every grid box every time step with a total sink at the surface every 
!  time step to reproduce GMI tracer mechanism.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE Set_Met_AgeOfAir( State_Met )
!
! !USES:
!
      USE CMN_SIZE_MOD
      USE State_Met_Mod,      ONLY : MetState
      USE TIME_MOD,           ONLY : GET_TS_DYN
!
! !OUTPUT PARAMETERS: 
!
      TYPE(MetState), INTENT(INOUT) :: State_Met
!
! !REVISION HISTORY: 
!  21 Dec 2018 - M. Sulprizio- Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER        :: I, J, L
      REAL(fp), SAVE :: TimeStep

      !=================================================================
      ! SET_MET_AGEOFAIR begins here!
      !=================================================================

      ! Get timestep [s]
      TimeStep = GET_TS_DYN()

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L )
      DO L=1,LLPAR
      DO J=1,JJPAR
      DO I=1,IIPAR

         IF ( L == 1 ) THEN
            ! Set the surface to a sink
            State_Met%AgeOfAir(I,J,L) = 0
         ELSE
            ! Otherwise add time step [s]
            State_Met%AgeOfAir(I,J,L) = State_Met%AgeOfAir(I,J,L) + 
     &                                  TimeStep
         ENDIF

      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      END SUBROUTINE Set_Met_AgeOfAir
!EOC
      END MODULE DAO_MOD
