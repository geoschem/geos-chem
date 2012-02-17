!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: dao_mod
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
      USE CMN_GCTM_MOD           ! Physical constants

      IMPLICIT NONE
#     include "define.h"

      PRIVATE

!
! !PUBLIC MEMBER FUNCTIONS:
!
      PUBLIC  :: AVGPOLE
      PUBLIC  :: AIRQNT
      PUBLIC  :: AIRQNT_FULLGRID
      PUBLIC  :: CLEANUP_DAO
      PUBLIC  :: CONVERT_UNITS
      PUBLIC  :: COPY_I6_FIELDS
      PUBLIC  :: COSSZA
      PUBLIC  :: GET_OBK
      PUBLIC  :: INIT_DAO
      PUBLIC  :: INTERP
      PUBLIC  :: IS_LAND
      PUBLIC  :: IS_WATER
      PUBLIC  :: IS_ICE
      PUBLIC  :: IS_NEAR
      PUBLIC  :: MAKE_AVGW
      PUBLIC  :: MAKE_RH
!
! !PUBLIC DATA MEMBERS:
!
      ! 2-D data fields
      REAL*8,  ALLOCATABLE, PUBLIC :: ALBD1   (:,:)
      REAL*8,  ALLOCATABLE, PUBLIC :: ALBD2   (:,:)
      REAL*8,  ALLOCATABLE, PUBLIC :: ALBD    (:,:)
      INTEGER, ALLOCATABLE, PUBLIC :: CLDTOPS (:,:)
      REAL*8,  ALLOCATABLE, PUBLIC :: CLDFRC  (:,:)
      REAL*8,  ALLOCATABLE, PUBLIC :: EFLUX   (:,:)
      REAL*8,  ALLOCATABLE, PUBLIC :: EVAP    (:,:)
      REAL*8,  ALLOCATABLE, PUBLIC :: FRLAKE  (:,:)
      REAL*8,  ALLOCATABLE, PUBLIC :: FRLAND  (:,:)
      REAL*8,  ALLOCATABLE, PUBLIC :: FROCEAN (:,:)
      REAL*8,  ALLOCATABLE, PUBLIC :: FRLANDIC(:,:)
      REAL*8,  ALLOCATABLE, PUBLIC :: FRSEAICE(:,:)
      REAL*8,  ALLOCATABLE, PUBLIC :: FRSNO   (:,:)
      REAL*8,  ALLOCATABLE, PUBLIC :: GRN     (:,:)
      REAL*8,  ALLOCATABLE, PUBLIC :: GWETROOT(:,:)
      REAL*8,  ALLOCATABLE, PUBLIC :: GWETTOP (:,:)
      REAL*8,  ALLOCATABLE, PUBLIC :: HFLUX   (:,:)
      REAL*8,  ALLOCATABLE, PUBLIC :: LAI     (:,:)
      REAL*8,  ALLOCATABLE, PUBLIC :: LWI_GISS(:,:)
      REAL*8,  ALLOCATABLE, PUBLIC :: LWI     (:,:)
      REAL*8,  ALLOCATABLE, PUBLIC :: MOLENGTH(:,:)
      REAL*8,  ALLOCATABLE, PUBLIC :: OICE    (:,:)      
      REAL*8,  ALLOCATABLE, PUBLIC :: PARDF   (:,:)
      REAL*8,  ALLOCATABLE, PUBLIC :: PARDR   (:,:)
      REAL*8,  ALLOCATABLE, PUBLIC :: PBL     (:,:)
      REAL*8,  ALLOCATABLE, PUBLIC :: PHIS    (:,:)
      REAL*8,  ALLOCATABLE, PUBLIC :: PREACC  (:,:)
      REAL*8,  ALLOCATABLE, PUBLIC :: PREANV  (:,:)
      REAL*8,  ALLOCATABLE, PUBLIC :: PRECON  (:,:)
      REAL*8,  ALLOCATABLE, PUBLIC :: PRELSC  (:,:)
      REAL*8,  ALLOCATABLE, PUBLIC :: PRECSNO (:,:)
      REAL*8,  ALLOCATABLE, PUBLIC :: PS1     (:,:)
      REAL*8,  ALLOCATABLE, PUBLIC :: PS2     (:,:)
      REAL*8,  ALLOCATABLE, PUBLIC :: PSC2    (:,:)
      REAL*8,  ALLOCATABLE, PUBLIC :: RADLWG  (:,:)
      REAL*8,  ALLOCATABLE, PUBLIC :: RADSWG  (:,:)
      REAL*8,  ALLOCATABLE, PUBLIC :: SEAICE00(:,:)
      REAL*8,  ALLOCATABLE, PUBLIC :: SEAICE10(:,:)
      REAL*8,  ALLOCATABLE, PUBLIC :: SEAICE20(:,:)
      REAL*8,  ALLOCATABLE, PUBLIC :: SEAICE30(:,:)
      REAL*8,  ALLOCATABLE, PUBLIC :: SEAICE40(:,:)
      REAL*8,  ALLOCATABLE, PUBLIC :: SEAICE50(:,:)
      REAL*8,  ALLOCATABLE, PUBLIC :: SEAICE60(:,:)
      REAL*8,  ALLOCATABLE, PUBLIC :: SEAICE70(:,:)
      REAL*8,  ALLOCATABLE, PUBLIC :: SEAICE80(:,:)
      REAL*8,  ALLOCATABLE, PUBLIC :: SEAICE90(:,:)
      REAL*8,  ALLOCATABLE, PUBLIC :: SLP     (:,:)
      REAL*8,  ALLOCATABLE, PUBLIC :: SNICE   (:,:)
      REAL*8,  ALLOCATABLE, PUBLIC :: SNODP   (:,:)
      REAL*8,  ALLOCATABLE, PUBLIC :: SNOMAS  (:,:)
      REAL*8,  ALLOCATABLE, PUBLIC :: SNOW    (:,:)
      REAL*8,  ALLOCATABLE, PUBLIC :: SUNCOS  (:  )
      REAL*8,  ALLOCATABLE, PUBLIC :: SUNCOS_MID(:  )
      REAL*8,  ALLOCATABLE, PUBLIC :: TO31    (:,:)
      REAL*8,  ALLOCATABLE, PUBLIC :: TO32    (:,:)
      REAL*8,  ALLOCATABLE, PUBLIC :: TO3     (:,:)
      REAL*8,  ALLOCATABLE, PUBLIC :: TTO3    (:,:)
      REAL*8,  ALLOCATABLE, PUBLIC :: TROPP1  (:,:)
      REAL*8,  ALLOCATABLE, PUBLIC :: TROPP2  (:,:)
      REAL*8,  ALLOCATABLE, PUBLIC :: TROPP   (:,:)
      REAL*8,  ALLOCATABLE, PUBLIC :: TS      (:,:)
      REAL*8,  ALLOCATABLE, PUBLIC :: TSKIN   (:,:)
      REAL*8,  ALLOCATABLE, PUBLIC :: U10M    (:,:)
      REAL*8,  ALLOCATABLE, PUBLIC :: USTAR   (:,:)  
      REAL*8,  ALLOCATABLE, PUBLIC :: V10M    (:,:)
      REAL*8,  ALLOCATABLE, PUBLIC :: Z0      (:,:)

      ! 3-D data fields
      REAL*8,  ALLOCATABLE, PUBLIC :: AD      (:,:,:)
      REAL*8,  ALLOCATABLE, PUBLIC :: AIRDEN  (:,:,:)
      REAL*8,  ALLOCATABLE, PUBLIC :: AIRVOL  (:,:,:)
      REAL*8,  ALLOCATABLE, PUBLIC :: AVGW    (:,:,:)
      REAL*8,  ALLOCATABLE, PUBLIC :: BXHEIGHT(:,:,:)
      REAL*8,  ALLOCATABLE, PUBLIC :: DQRCU   (:,:,:)
      REAL*8,  ALLOCATABLE, PUBLIC :: DQRLSAN (:,:,:)
      REAL*8,  ALLOCATABLE, PUBLIC :: CLDF    (:,:,:)
      REAL*8,  ALLOCATABLE, PUBLIC :: CLDMAS  (:,:,:)
      REAL*8,  ALLOCATABLE, PUBLIC :: CMFMC   (:,:,:)
      REAL*8,  ALLOCATABLE, PUBLIC :: DELP    (:,:,:)
      REAL*8,  ALLOCATABLE, PUBLIC :: DETRAINE(:,:,:)
      REAL*8,  ALLOCATABLE, PUBLIC :: DETRAINN(:,:,:)
      REAL*8,  ALLOCATABLE, PUBLIC :: DNDE    (:,:,:)
      REAL*8,  ALLOCATABLE, PUBLIC :: DNDN    (:,:,:)
      REAL*8,  ALLOCATABLE, PUBLIC :: DQIDTMST(:,:,:)
      REAL*8,  ALLOCATABLE, PUBLIC :: DQLDTMST(:,:,:)
      REAL*8,  ALLOCATABLE, PUBLIC :: DQRCON  (:,:,:)
      REAL*8,  ALLOCATABLE, PUBLIC :: DQRLSC  (:,:,:)    
      REAL*8,  ALLOCATABLE, PUBLIC :: DQVDTMST(:,:,:)
      REAL*8,  ALLOCATABLE, PUBLIC :: DTRAIN  (:,:,:)
      REAL*8,  ALLOCATABLE, PUBLIC :: ENTRAIN (:,:,:)
      REAL*8,  ALLOCATABLE, PUBLIC :: HKBETA  (:,:,:)
      REAL*8,  ALLOCATABLE, PUBLIC :: HKETA   (:,:,:)
      REAL*8,  ALLOCATABLE, PUBLIC :: MFXC    (:,:,:)
      REAL*8,  ALLOCATABLE, PUBLIC :: MFYC    (:,:,:)
      REAL*8,  ALLOCATABLE, PUBLIC :: MFZ     (:,:,:)
      REAL*8,  ALLOCATABLE, PUBLIC :: MOISTQ  (:,:,:)
      REAL*8,  ALLOCATABLE, PUBLIC :: OPTDEP  (:,:,:)
      REAL*8,  ALLOCATABLE, PUBLIC :: OPTD    (:,:,:)
      REAL*8,  ALLOCATABLE, PUBLIC :: PFICU   (:,:,:)
      REAL*8,  ALLOCATABLE, PUBLIC :: PFILSAN (:,:,:)
      REAL*8,  ALLOCATABLE, PUBLIC :: PFLCU   (:,:,:)
      REAL*8,  ALLOCATABLE, PUBLIC :: PFLLSAN (:,:,:)
      REAL*8,  ALLOCATABLE, PUBLIC :: PV      (:,:,:)
      REAL*8,  ALLOCATABLE, PUBLIC :: QI      (:,:,:)
      REAL*8,  ALLOCATABLE, PUBLIC :: QL      (:,:,:)
      REAL*8,  ALLOCATABLE, PUBLIC :: REEVAPCN(:,:,:)
      REAL*8,  ALLOCATABLE, PUBLIC :: REEVAPLS(:,:,:)
      REAL*8,  ALLOCATABLE, PUBLIC :: RH1     (:,:,:)
      REAL*8,  ALLOCATABLE, PUBLIC :: RH2     (:,:,:)
      REAL*8,  ALLOCATABLE, PUBLIC :: RH      (:,:,:)
      REAL*8,  ALLOCATABLE, PUBLIC :: SPHU1   (:,:,:)
      REAL*8,  ALLOCATABLE, PUBLIC :: SPHU2   (:,:,:)
      REAL*8,  ALLOCATABLE, PUBLIC :: SPHU    (:,:,:)
      REAL*8,  ALLOCATABLE, PUBLIC :: T       (:,:,:)
      REAL*8,  ALLOCATABLE, PUBLIC :: TAUCLI  (:,:,:)
      REAL*8,  ALLOCATABLE, PUBLIC :: TAUCLW  (:,:,:)
      REAL*8,  ALLOCATABLE, PUBLIC :: TMPU1   (:,:,:)
      REAL*8,  ALLOCATABLE, PUBLIC :: TMPU2   (:,:,:)
      REAL*8,  ALLOCATABLE, PUBLIC :: UPDE    (:,:,:)
      REAL*8,  ALLOCATABLE, PUBLIC :: UPDN    (:,:,:)
      REAL*8,  ALLOCATABLE, PUBLIC :: UWND1   (:,:,:)
      REAL*8,  ALLOCATABLE, PUBLIC :: UWND2   (:,:,:)
      REAL*8,  ALLOCATABLE, PUBLIC :: UWND    (:,:,:)
      REAL*8,  ALLOCATABLE, PUBLIC :: VWND1   (:,:,:)
      REAL*8,  ALLOCATABLE, PUBLIC :: VWND2   (:,:,:)
      REAL*8,  ALLOCATABLE, PUBLIC :: VWND    (:,:,:)
      REAL*8,  ALLOCATABLE, PUBLIC :: ZMEU    (:,:,:)
      REAL*8,  ALLOCATABLE, PUBLIC :: ZMMD    (:,:,:)
      REAL*8,  ALLOCATABLE, PUBLIC :: ZMMU    (:,:,:)
      REAL*8,  ALLOCATABLE, PUBLIC :: AIRDEN_FULLGRID(:,:,:)
      REAL*8,  ALLOCATABLE, PUBLIC :: T_FULLGRID     (:,:,:)
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
!EOP
!------------------------------------------------------------------------------
!BOC
      CONTAINS
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
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
      USE GRID_MOD, ONLY : GET_AREA_M2
!
! !INPUT/OUTPUT PARAMETERS: 
!
      REAL*8, INTENT(INOUT) :: Z(IIPAR,JJPAR)   ! Quantity to be averaged 
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
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER :: I, J
      REAL*8  :: TOTAL_Z1, TOTAL_Z2, TOTAL_Z3, TOTAL_Z4
      REAL*8  :: TOTAL_A1, TOTAL_A2, TOTAL_A3, TOTAL_A4

      !=================================================================
      ! AVGPOLE begins here!                                                  
      !=================================================================

#if   defined( GRID1x1   ) || defined( GRID05x0666 ) 
#if   defined( NESTED_CH ) || defined( NESTED_NA   ) 
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
         TOTAL_Z1 = TOTAL_Z1 + GET_AREA_M2(       1 ) * Z(I,      1)
         TOTAL_Z2 = TOTAL_Z2 + GET_AREA_M2(       2 ) * Z(I,      2)
         TOTAL_Z3 = TOTAL_Z3 + GET_AREA_M2( JJPAR-1 ) * Z(I,JJPAR-1)
         TOTAL_Z4 = TOTAL_Z4 + GET_AREA_M2(   JJPAR ) * Z(I,JJPAR  )
         TOTAL_A1 = TOTAL_A1 + GET_AREA_M2(       1 ) 
         TOTAL_A2 = TOTAL_A2 + GET_AREA_M2(       2 )
         TOTAL_A3 = TOTAL_A3 + GET_AREA_M2( JJPAR-1 )
         TOTAL_A4 = TOTAL_A4 + GET_AREA_M2(   JJPAR )
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
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: airqnt
!
! !DESCRIPTION: Subroutine AIRQNT calculates the volume [m\^3 and cm\^3], mass 
!  [kg], density, [kg/m\^3], and pressure thickness [hPa] of air for each grid 
!  box (I,J,L).   The quantity (surface pressure - PTOP) [hPa] at each surface 
!  grid box (I,J) is also computed. 
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE AIRQNT
!
! !USES:
!
      USE GRID_MOD,     ONLY : GET_AREA_M2
      USE PRESSURE_MOD, ONLY : GET_BP, GET_PEDGE
!
! !REMARKS:
!  DAO met fields updated by AIRQNT:
!  ========================================================================
!  (1 ) BXHEIGHT (REAL*8 ) : Vertical extent of a grid box   [m       ]
!  (2 ) DELP     (REAL*8 ) : Delta-P extent  of a grid box   [mb      ]
!  (3 ) AIRVOL   (REAL*8 ) : Volume  of air  in a grid box   [m^3     ]
!  (4 ) AD       (REAL*8 ) : Mass    of air  in a grid box   [kg      ]
!  (5 ) AIRDEN   (REAL*8 ) : Density of air  in a grid box   [kg/m^3  ]
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
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER :: I,  J,  L
      REAL*8  :: P1, P2, AREA_M2

      !=================================================================
      ! AIRQNT begins here! 
      !=================================================================
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, AREA_M2, P1, P2 )
      DO L = 1, LLPAR
      DO J = 1, JJPAR

         ! Grid box surface area [m2]
         AREA_M2 = GET_AREA_M2( J )

         DO I = 1, IIPAR
               
            ! Pressure at bottom edge of grid box [hPa]
            P1          = GET_PEDGE(I,J,L)

            ! Pressure at top edge of grid box [hPa]
            P2          = GET_PEDGE(I,J,L+1)

            ! Pressure difference between top & bottom edges [hPa]
            DELP(L,I,J) = P1 - P2
            
            !===========================================================
            ! BXHEIGHT is the height (Delta-Z) of grid box (I,J,L) 
            ! in meters. 
            !
            ! The formula for BXHEIGHT is just the hydrostatic eqn.  
            ! Rd = 287 J/K/kg is the value for the ideal gas constant
            ! R for air (M.W = 0.02897 kg/mol),  or 
            ! Rd = 8.31 J/(mol*K) / 0.02897 kg/mol. 
            !===========================================================
            BXHEIGHT(I,J,L) = Rdg0 * T(I,J,L) * LOG( P1 / P2 )

            !===========================================================
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
            !===========================================================
            AIRVOL(I,J,L) = BXHEIGHT(I,J,L) * AREA_M2

            !===========================================================
            ! AD = (dry) mass of air in grid box (I,J,L) in kg, 
            ! given by:        
            !
            !  Mass    Pressure        100      1        Surface area 
            !        = difference   *  ---  *  ---   *   of grid box 
            !          in grid box      1       g          AREA_M2
            !
            !   kg         mb          Pa      s^2           m^2
            !  ----  =    ----      * ----  * -----  *      -----
            !    1          1          mb       m             1
            !===========================================================
            AD(I,J,L) = DELP(L,I,J) * G0_100 * AREA_M2

            !===========================================================
            ! AIRDEN = density of air (AD / AIRVOL) in kg / m^3 
            !===========================================================
            AIRDEN(L,I,J) = AD(I,J,L) / AIRVOL(I,J,L)
         ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      END SUBROUTINE AIRQNT
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: airqnt_fullgrid
!
! !DESCRIPTION: Subroutine AIRQNT\_FULLGRID calculates the same quantities as 
!  AIRQNT, but for the full, unlumped vertical grid of the GEOS GCM.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE AIRQNT_FULLGRID
!
! !USES:
!
      USE GRID_MOD,     ONLY : GET_AREA_M2
      USE PRESSURE_MOD, ONLY : GET_PEDGE_FULLGRID
!
! !REMARKS:
!  DAO met fields updated by AIRQNT_FULLGRID:
!  ========================================================================
!  (1 ) AIRDEN_FULLGRID (REAL*8 ) : Density of air  in a grid box   [kg/m^3  ]
!
!  NOTES:
!  (1 ) Modified from AIRQNT in DAO_MOD (cdh, 1/22/09)
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER             :: I,  J,  L
      REAL*8              :: P1, P2, AREA_M2, DELP_TMP
      REAL*8              :: BXHEIGHT_TMP, AIRVOL_TMP, AD_TMP

      !=================================================================
      ! AIRQNT_FULLGRID begins here! 
      !=================================================================
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, AREA_M2, P1, P2 )
!$OMP+PRIVATE( DELP_TMP, BXHEIGHT_TMP, AIRVOL_TMP, AD_TMP )
      DO L = 1, LGLOB
      DO J = 1, JJPAR

         ! Grid box surface area [m2]
         AREA_M2 = GET_AREA_M2( J )

         DO I = 1, IIPAR
               
            ! Pressure at bottom edge of grid box [hPa]
            P1          = GET_PEDGE_FULLGRID(I,J,L)

            ! Pressure at top edge of grid box [hPa]
            P2          = GET_PEDGE_FULLGRID(I,J,L+1)

            ! Pressure difference between top & bottom edges [hPa]
            DELP_TMP = P1 - P2
            
            !===========================================================
            ! BXHEIGHT is the height (Delta-Z) of grid box (I,J,L) 
            ! in meters. 
            !
            ! The formula for BXHEIGHT is just the hydrostatic eqn.  
            ! Rd = 287 J/K/kg is the value for the ideal gas constant
            ! R for air (M.W = 0.02897 kg/mol),  or 
            ! Rd = 8.31 J/(mol*K) / 0.02897 kg/mol. 
            !===========================================================
            BXHEIGHT_TMP = Rdg0 * T_FULLGRID(I,J,L) * LOG( P1 / P2 )

            !===========================================================
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
            !===========================================================
            AIRVOL_TMP = BXHEIGHT_TMP * AREA_M2

            !===========================================================
            ! AD = (dry) mass of air in grid box (I,J,L) in kg, 
            ! given by:        
            !
            !  Mass    Pressure        100      1        Surface area 
            !        = difference   *  ---  *  ---   *   of grid box 
            !          in grid box      1       g          AREA_M2
            !
            !   kg         mb          Pa      s^2           m^2
            !  ----  =    ----      * ----  * -----  *      -----
            !    1          1          mb       m             1
            !===========================================================
            AD_TMP = DELP_TMP * G0_100 * AREA_M2

            !===========================================================
            ! AIRDEN = density of air (AD / AIRVOL) in kg / m^3 
            !===========================================================
            AIRDEN_FULLGRID(L,I,J) = AD_TMP / AIRVOL_TMP
         ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      END SUBROUTINE AIRQNT_FULLGRID
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: interp
!
! !DESCRIPTION: Subroutine INTERP linearly interpolates GEOS-Chem I6 fields 
!  (winds, surface pressure, temperature, surface albedo, specific humidity
!   etc.)  to the current dynamic timestep.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE INTERP( NTIME0, NTIME1, NTDT )
!
! !USES:
!
      USE GRID_MOD,    ONLY : GET_YEDGE
      USE LOGICAL_MOD, ONLY : LVARTROP
!
! !INPUT PARAMETERS: 
!
      INTEGER, INTENT(IN)  :: NTIME0   ! Elapsed time [s] at start of 6hr step
      INTEGER, INTENT(IN)  :: NTIME1   ! Elapsed time [s] at current time
      INTEGER, INTENT(IN)  :: NTDT     ! Length of dynamic timestep [s]
! 
! !REMARKS:
!  Different met fields are archived at I6 (instantaneous 6-hr) time 
!  resolution depending on the specific product.  For example, relative 
!  humidity is an instantaneous 6hr field in MERRA and a 6-hr time averaged
!  field in GEOS-5.
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
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER :: I,        J,        L
      REAL*8  :: D_NTIME0, D_NTIME1, D_NDT
      REAL*8  :: D_NTDT,   TM,       TC2
      REAL*8  :: YSOUTH,   YNORTH

      !=================================================================
      ! Initialization
      !=================================================================

      ! Convert time variables from FLOAT to DBLE
      D_NTIME0 = DBLE( NTIME0 )
      D_NTIME1 = DBLE( NTIME1 )
      D_NTDT   = DBLE( NTDT   )
      D_NDT    = 21600d0

      ! Fraction of 6h timestep elapsed at mid point of this dyn timestep
      TM  = ( D_NTIME1 + D_NTDT/2d0 - D_NTIME0 ) / D_NDT
      
      ! Fraction of 6h timestep elapsed at the end of this dyn timestep
      TC2 = ( D_NTIME1 + D_NTDT     - D_NTIME0 ) / D_NDT 

#if defined( GEOS_3 )

      !=================================================================
      ! Do the interpolation for GEOS-3 met fields 
      !=================================================================
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L )
      DO L = 1, LLPAR
      DO J = 1, JJPAR
      DO I = 1, IIPAR
         
         !----------------------------------------------------
         ! Interpolate 2D variables
         !----------------------------------------------------
         IF ( L == 1 ) THEN
            
            ! Interpolate pressure [hPa] to the end of the dynamic timestep
            PSC2(I,J)  = PS1(I,J)  + ( PS2(I,J) - PS1(I,J) ) * TC2 
  
            ! Interpolate ALBEDO to the midpt of the dynamic timestep
            ALBD(I,J) = ALBD1(I,J) + ( ALBD2(I,J) - ALBD1(I,J) ) * TM

            ! Interpolate TROPP to the midpt of the dynamic timestep
            IF ( LVARTROP ) THEN
               TROPP(I,J) = TROPP1(I,J) 
     &                    + ( TROPP2(I,J) - TROPP1(I,J) ) * TM
            ENDIF
         ENDIF
         
         !-----------------------------------------------------
         ! Interpolate 3D Variables
         !-----------------------------------------------------

         ! Interpolate to the midpt of the dynamic timestep
         UWND(I,J,L) = UWND1(I,J,L) + (UWND2(I,J,L) - UWND1(I,J,L)) * TM
         VWND(I,J,L) = VWND1(I,J,L) + (VWND2(I,J,L) - VWND1(I,J,L)) * TM
         SPHU(I,J,L) = SPHU1(I,J,L) + (SPHU2(I,J,L) - SPHU1(I,J,L)) * TM
         T(I,J,L)    = TMPU1(I,J,L) + (TMPU2(I,J,L) - TMPU1(I,J,L)) * TM

      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

#elif defined( GEOS_4 ) || defined( GCAP )

      !=================================================================
      ! Do the interpolation for GEOS-4 or GCAP met fields 
      !=================================================================
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, YSOUTH, YNORTH )
      DO J = 1, JJPAR

         ! North & south edges of box
         YSOUTH = GET_YEDGE( J   )
         YNORTH = GET_YEDGE( J+1 )

         DO I = 1, IIPAR

            ! Interpolate pressure [hPa] to the end of the dynamic timestep
            PSC2(I,J) = PS1(I,J) + ( PS2(I,J) - PS1(I,J) ) * TC2 

            ! Test if we are using the variable tropopause
            IF ( LVARTROP ) THEN
 
               ! Interpolate TROPP [hPa] to the midpt of the dynamic timestep
               TROPP(I,J) = TROPP1(I,J) 
     &                    + ( TROPP2(I,J) - TROPP1(I,J) ) * TM


               ! However, we still need to make sure to cap TROPP in the
               ! polar regions (if the entire box is outside 60S-60N)
               ! so that we don't do chemistry at an abnormally high
               ! altitude.  Set TROPP in the polar regions to 200 hPa.
               ! (jal, phs, bmy, 9/18/07)
               IF ( YSOUTH >= 60d0 .or. YNORTH <= -60d0 ) THEN
                  TROPP(I,J) = MAX( TROPP(I,J), 200d0 )
               ENDIF
            ENDIF
         ENDDO
      ENDDO
!$OMP END PARALLEL DO

#elif defined( GEOS_5 )

      !=================================================================
      ! Do the interpolation for GEOS-5 met fields 
      !=================================================================
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, YSOUTH, YNORTH )
      DO J = 1, JJPAR

         ! North & south edges of box
         YSOUTH = GET_YEDGE( J   )
         YNORTH = GET_YEDGE( J+1 )

         DO I = 1, IIPAR

            ! Interpolate pressure [hPa] to the end of the dynamic timestep
            PSC2(I,J) = PS1(I,J)  + ( PS2(I,J)  - PS1(I,J)  ) * TC2 

            ! Interpolate the GEOS-5 total O3 column [DU] to current time
            TO3(I,J)  = TO31(I,J) + ( TO32(I,J) - TO31(I,J) ) * TM

            ! Even though TROPP is a 3-hour average field, we 
            ! we still need to make sure to cap TROPP in the
            ! polar regions (if the entire box is outside 60S-60N)
            ! so that we don't do chemistry at an abnormally high
            ! altitude.  Set TROPP in the polar regions to 200 hPa.
            ! (jal, phs, bmy, 9/18/07)
            IF ( LVARTROP ) THEN
               IF ( YSOUTH >= 60d0 .or. YNORTH <= -60d0 ) THEN
                  TROPP(I,J) = MAX( TROPP(I,J), 200d0 )
               ENDIF
            ENDIF
         ENDDO
      ENDDO
!$OMP END PARALLEL DO

#elif defined( MERRA )

      !=================================================================
      ! Do the interpolation for MERRA met fields 
      !=================================================================
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, YSOUTH, YNORTH )
      DO L = 1, LLPAR
      DO J = 1, JJPAR

         ! North & south edges of box
         YSOUTH = GET_YEDGE( J   )
         YNORTH = GET_YEDGE( J+1 )

      DO I = 1, IIPAR

         !----------------------------------------------------
         ! Interpolate 2D variables
         !----------------------------------------------------
         IF ( L == 1 ) THEN
            
            ! Interpolate pressure [hPa] to the end of the dynamic timestep
            PSC2(I,J)  = PS1(I,J) + ( PS2(I,J) - PS1(I,J) ) * TC2 

            ! Even though TROPP is a 3-hour average field, we 
            ! we still need to make sure to cap TROPP in the
            ! polar regions (if the entire box is outside 60S-60N)
            ! so that we don't do chemistry at an abnormally high
            ! altitude.  Set TROPP in the polar regions to 200 hPa.
            ! (jal, phs, bmy, 9/18/07)
            IF ( LVARTROP ) THEN
               IF ( YSOUTH >= 60d0 .or. YNORTH <= -60d0 ) THEN
                  TROPP(I,J) = MAX( TROPP(I,J), 200d0 )
               ENDIF
            ENDIF
         ENDIF

         !----------------------------------------------------
         ! Interpolate 3D variables
         !----------------------------------------------------
         RH(I,J,L) = RH1(I,J,L) + ( RH2(I,J,L) - RH1(I,J,L) ) * TM

      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

#endif

      END SUBROUTINE INTERP
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
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
      FUNCTION IS_LAND( I, J ) RESULT ( LAND )
!
! !USES:
!
      USE TIME_MOD, ONLY : GET_YEAR
!
! !INPUT PARAMETERS: 
!
      INTEGER, INTENT(IN)  :: I           ! Longitude index of grid box
      INTEGER, INTENT(IN)  :: J           ! Latitude  index of grid box
!
! !RETURN VALUE:
!
      LOGICAL              :: LAND        ! =T if it is a land box
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
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
#if   defined( GEOS_3 )

      !--------------------------
      ! GEOS-3
      !--------------------------
      IF ( GET_YEAR() == 1998 ) THEN

         ! Fields for 1998 don't have LWI/SURFTYPE flags, so use albedo 
         ! as a proxy for land coverage instead: 0.08 < ALBEDO < 0.55
         LAND = ( ALBD(I,J) > 0.08d0 .and. ALBD(I,J) < 0.55d0 )

      ELSE

         ! Otherwise LWI < 50 and ALBEDO less than 69.5% is a water box 
         LAND = ( LWI(I,J) < 50 .and. ALBD(I,J) < 0.695d0 )

      ENDIF

#elif defined( GEOS_4 ) || defined( GEOS_5 ) || defined( MERRA )

      !--------------------------
      ! GEOS-4 / GEOS-5 / MERRA
      !--------------------------

      ! LWI=1 and ALBEDO less than 69.5% is a LAND box 
      LAND = ( NINT( LWI(I,J) ) == 1 .and. ALBD(I,J) < 0.695d0 )

#elif defined( GCAP )

      !--------------------------
      ! GCAP
      !--------------------------

      ! It's a land box if 50% or more of the box is covered by 
      ! land and less than 50% of the box is covered by ice
      LAND = ( LWI_GISS(I,J) >= 0.5d0 .and. SNICE(I,J) < 0.5d0 )
    
#endif

      END FUNCTION IS_LAND
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
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
      FUNCTION IS_WATER( I, J ) RESULT ( WATER )
!
! !USES:
!
      USE TIME_MOD, ONLY : GET_YEAR
!
! !INPUT PARAMETERS: 
!
      INTEGER, INTENT(IN)  :: I           ! Longitude index of grid box
      INTEGER, INTENT(IN)  :: J           ! Latitude  index of grid box
!
! !RETURN VALUE:
!
      LOGICAL              :: WATER       ! =T if this is a water box
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
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
#if   defined( GEOS_3 )
      
      !--------------------------
      ! GEOS-3 
      !--------------------------
      IF ( GET_YEAR() == 1998 ) THEN

         ! 1998 fields don't have LWI/SURFTYPE flags, so use albedo as 
         ! a proxy for water coverage: 55%  < ALBEDO < 69.5%
         WATER = ( ALBD(I,J) > 0.55d0 .and. ALBD(I,J) < 0.695d0 )

      ELSE

         ! Otherwise LWI >= 50 and ALBEDO less than 69.5% is a water box
         WATER = ( LWI(I,J) >= 50 .and. ALBD(I,J) < 0.695d0 )
         
      ENDIF

#elif defined( GEOS_4 ) || defined( GEOS_5 ) || defined( MERRA )
      
      !---------------------------
      ! GEOS-4 / GEOS-5 / MERRA
      !---------------------------

      ! LWI=0 and ALBEDO less than 69.5% is a water box 
      WATER = ( NINT( LWI(I,J) ) == 0 .and. ALBD(I,J) < 0.695d0 )

#elif defined( GCAP )

      !--------------------------
      ! GCAP
      !--------------------------

      ! It's a water box if less than 50% of the box is
      ! covered by land and less than 50% is covered by ice
      WATER = ( LWI_GISS(I,J) < 0.5d0 .and. SNICE(I,J) < 0.5d0 )

#endif

      END FUNCTION IS_WATER
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
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
      FUNCTION IS_ICE( I, J ) RESULT ( ICE )
!
! !USES:
!
      USE TIME_MOD, ONLY : GET_YEAR
!
! !INPUT PARAMETERS: 
!
      INTEGER, INTENT(IN)  :: I           ! Longitude index of grid box
      INTEGER, INTENT(IN)  :: J           ! Latitude  index of grid box
!
! !RETURN VALUE:
!
      LOGICAL              :: ICE         ! =T if this is an ice box
!
! 
! !REVISION HISTORY: 
!  09 Aug 2005 - R. Yantosca - Initial version
!  (1 ) Remove support for GEOS-1 and GEOS-STRAT met fields (bmy, 8/4/06)
!  16 Aug 2010 - R. Yantosca - Added ProTeX headers
!  25 Aug 2010 - R. Yantosca - Treat MERRA in the same way as GEOS-5
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
#if   defined( GEOS_3 )

      !--------------------------
      ! GEOS-3
      !-------------------------- 

      ! Fields for 1998 don't have LWI/SURFTYPE flags, so use albedo 
      ! as a proxy for water coverage instead: ALBEDO > 0.695
      ICE = ( ALBD(I,J) >= 0.695d0 )

#elif defined( GEOS_4 ) || defined( GEOS_5 ) || defined( MERRA )

      !--------------------------
      ! GEOS-4 / GEOS-5 / MERRA
      !--------------------------

      ! LWI=2 or ALBEDO > 69.5% is ice
      ICE = ( NINT( LWI(I,J) ) == 2 .or. ALBD(I,J) >= 0.695d0 )

#elif defined( GCAP )

      !--------------------------
      ! GCAP
      !--------------------------

      ! It's an ice box if 50% or more of the box is covered by ice
      ICE = ( SNICE(I,J) >= 0.5d0 )

#endif

      END FUNCTION IS_ICE
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
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
      FUNCTION IS_NEAR( I, J, THRESH, NEIGHBOR ) RESULT ( NEAR )
!
! !INPUT PARAMETERS: 
!
      ! Arguments
      INTEGER, INTENT(IN) :: I, J        ! Lon & lat grid box indices
      INTEGER, INTENT(IN) :: NEIGHBOR    ! # of neighbor boxes to consider
      REAL*8,  INTENT(IN) :: THRESH      ! LWI threshold for near-land 
!
! !RETURN VALUE:
!
      LOGICAL             :: NEAR        ! # of near land boxes
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
            IF ( IS_ICE( LONGI, LATJ ) ) CYCLE

#if   defined( GCAP ) 

            !---------------------------------------------------
            ! GCAP met fields
            !
            ! LWI_GISS = 0.0 means that the box is all water
            ! LWI_GISS = 1.0 means that the box is all land
            !
            ! with fractional values at land-water boundaries
            !
            ! It's near-land if THRESH <= LWI_GISS <= 1.0 
            !---------------------------------------------------
            IF ( LWI_GISS(LONGI,LATJ) >  THRESH .and.
     &           LWI_GISS(LONGI,LATJ) <= 1.0d0 ) THEN

               ! We are in a near-land box
               NEAR = .TRUE.

               ! Break out of loop
               GOTO 999
            ENDIF

#elif defined( GEOS_3 )

            !---------------------------------------------------
            ! GEOS-3 met fields
            !
            ! LWI < 10 is land
            ! LWI = 101 is water
            !
            ! with fractional values at land-water boundaries
            !
            ! Therefore if you pick a threshold value such
            ! as 80, then everything with LWI < THRESH is 
            ! sure to be a land box.
            !
            ! It's near land if LWI < THRESH.
            !---------------------------------------------------
            IF ( LWI(LONGI,LATJ) < THRESH ) THEN 

               ! We are in a near-land box
               NEAR = .TRUE.

               ! Break out of loop
               GOTO 999
            ENDIF

#elif defined( GEOS_4 ) || defined( GEOS_5 ) || defined( MERRA )

            !---------------------------------------------------
            ! GEOS-4 / GEOS-5 / MERRA met fields
            !
            ! LWI = 0.0 is ocean
            ! LWI = 1.0 is land
            ! LWI = 2.0 is ice 
            !
            ! with fractional values at land-water, land-ice,
            ! and water-ice boundaries.
            !
            ! It's near-land if THRESH <= LWI <= 1.0 
            !---------------------------------------------------
            IF ( LWI(LONGI,LATJ) >  THRESH  .and.
     &           LWI(LONGI,LATJ) <= 1d0    ) THEN

               ! We are in a near-land box
               NEAR = .TRUE.

               ! Break out of loop
               GOTO 999
            ENDIF

#endif

         ENDDO
      ENDDO

      ! Exit
 999  CONTINUE

      END FUNCTION IS_NEAR
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: make_avgw
!
! !DESCRIPTION: Subroutine MAKE\_AVGW converts DAO specific humidity SPHU 
!  to AVGW, which is the mixing ratio of water vapor. 
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE MAKE_AVGW
! 
! !REVISION HISTORY: 
!  30 Jan 1998 - R. Yantosca - Initial version
!  (1 ) AVGW was originally indexed by (L,I,J).  Reorder the indexing to
!        (I,J,L) to take advantage of the way FORTRAN stores by columns.
!        An (L,I,J) ordering can lead to excessive disk swapping.
!  (2 ) Now dimension AVGW as (IIPAR,JJPAR,LLPAR).  Also use parallel
!        DO-loop to compute AVGW.  Updated comments. (bmy, 9/24/01)
!  (3 ) Removed obsolete code from 9/01 (bmy, 10/23/01)
!  (4 ) SPHU and AVGW are declared w/in "dao_mod.f", so we don't need to pass
!        these as arguments anymore (bmy, 11/15/02)
!  16 Aug 2010 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER             :: I, IREF, J, JREF, L      
!
! !DEFINED PARAMETERS:
!
      ! Conversion factor
      REAL*8, PARAMETER   :: HCONV = 28.97d-3 / 18.0d0 

      !=================================================================
      ! MAKE_AVGW begins here!
      !
      ! In the original Harvard/GISS/Irvine CTM subroutines, 
      !    AVGW = log10( mixing ratio of water vapor ).  
      !
      ! In order to avoid costly log and exponentiation operations, 
      ! redefine AVGW, so that AVGW is the actual mixing ratio of water 
      ! vapor, and not the log10 of the mixing ratio.
      !
      ! The conversion from SPHU [g H2O/kg air] to [v/v] mixing ratio is:
      !
      !   g H2O  | mol H2O  | 28.97e-3 kg air    mol H2O     vol H2O
      ! ---------+----------+---------------- = --------- = ---------
      !   kg air | 18 g H2O |    mol air         mol air     vol air
      !
      !      thus AVGW (V/V) = SPHU (g/kg) * HCONV, 
      !
      ! where HCONV = the conversion factor ( 28.97e-3 / 18.0 ), 
      ! which is defined as a parameter at the top of this routine.
      !=================================================================
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L )
      DO L = 1, LLPAR
      DO J = 1, JJPAR
      DO I = 1, IIPAR
         AVGW(I,J,L) = HCONV * SPHU(I,J,L) 
      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      END SUBROUTINE MAKE_AVGW
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: make_rh
!
! !DESCRIPTION: Subroutine MAKE\_RH computes relative humidity from specific 
!  humidity and temperature. 
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE MAKE_RH
!
! !USES:
!
      USE PRESSURE_MOD, ONLY : GET_PCENTER
!
! !REMARKS:
!  Module variables used:
!  ===========================================================================
!  (1 ) SPHU (REAL*8) : Array containing 3-D specific humidity [g H2O/kg air]
!  (2 ) TMPU (REAL*8) : Array containing 3-D temperature field [K]
!  (3 ) RH   (REAL*8) : Output array for relative humidity     [%]
! 
! !REVISION HISTORY: 
!  13 Oct 1999 - R. Yantosca - Initial version
!  (1 ) Use F90 syntax for declarations, etc. 
!  (2 ) Cosmetic changes (bmy, 10/12/99)
!  (3 ) Now use GET_PCENTER from "pressure_mod.f" to compute the pressure
!        at the midpoint of grid box (I,J,L).  Updated comments, cosmetic
!        changes.  Added parallel DO-loops.  Remove reference to "CMN" 
!        header file.  Added to "dao_mod.f" (dsa, bdf, bmy, 8/27/02)
!  (4 ) Removed obsolete code from 8/02 (bmy, 9/18/02)
!  (5 ) Now remove SPHU, TMPU, RH from the arg list, since these are now
!        all contained w/in this dao_mod.f as module variables. (bmy, 9/23/02)
!  16 Aug 2010 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      REAL*8              :: ESAT, SHMB, PRES, TEMP
      INTEGER             :: I, J, L
!
! !DEFINED PARAMETERS:
!
      REAL*8, PARAMETER   :: A =   23.5518d0
      REAL*8, PARAMETER   :: B = 2937.4d0
      REAL*8, PARAMETER   :: C =   -4.9283d0

      !=================================================================
      ! MAKE_RH begins here!
      !=================================================================
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, PRES, TEMP, ESAT, SHMB )
      DO L = 1, LLPAR
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! Pressure at midpoint of box (I,J,L)
         PRES = GET_PCENTER(I,J,L)

         ! Temperature at grid box (I,J,L)
         TEMP = T(I,J,L)

         ! Saturation water vapor pressure in mbar 
         ! (from NASA GTE PEM-Tropics handbook)
         ESAT = ( 10d0**( A - ( B / TEMP ) ) ) * ( TEMP**C )
            
         ! Specific humidity in mb
         SHMB = SPHU(I,J,L) * 1.6072d-3 * PRES
            
         ! Relative humidity as a percentage
         RH(I,J,L) = ( SHMB / ESAT ) * 100d0 

      ENDDO
      ENDDO
      ENDDO  
!$OMP END PARALLEL DO

      END SUBROUTINE MAKE_RH
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
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
      FUNCTION GET_OBK( I, J ) RESULT( OBK )
!
! !INPUT PARAMETERS: 
!
      INTEGER, INTENT(IN) :: I     ! Longitude index
      INTEGER, INTENT(IN) :: J     ! Latitude  index
!
! !RETURN VALUE:
!
      REAL*8              :: OBK   ! Monin-Obhukhov length
!
! !REMARKS:
! 
! 
! !REVISION HISTORY: 
!  25 May 2005 - R. Yantosca - Initial version
!  16 Aug 2010 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC

#if   defined( GCAP )

      !=================================================================
      ! For GCAP met fields (based on GISS model)
      !=================================================================

      ! Monin-Obhukov length is a GCAP met field
      OBK = MOLENGTH(I,J)

#else

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
      REAL*8            :: NUM, DEN

      ! Parameters
      REAL*8, PARAMETER :: KAPPA = 0.4d0 
      REAL*8, PARAMETER :: CP    = 1000.0d0

      ! Numerator
      NUM = -AIRDEN(1,I,J) *  CP            * TS(I,J) *
     &       USTAR(I,J)    *  USTAR(I,J)    * USTAR(I,J)

      ! Denominator
      DEN =  KAPPA * g0 * HFLUX(I,J) 

      ! Prevent div by zero
      IF ( ABS( DEN ) > 0d0 ) THEN
         OBK = NUM / DEN
      ELSE
         OBK = 1.0d5
      ENDIF

#endif

      END FUNCTION GET_OBK
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: cossza
!
! !DESCRIPTION: COSSZA computes the cosine of the solar zenith angle at
!  the current time, and at the midpoint of the chemistry timestep interval.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE COSSZA( JDAY, SUNCOS, SUNCOS_MID )
!
! !USES:
!
      USE GRID_MOD, ONLY : GET_YMID_R
      USE TIME_MOD, ONLY : GET_HOUR
      USE TIME_MOD, ONLY : GET_MINUTE
      USE TIME_MOD, ONLY : GET_LOCALTIME
      USE TIME_MOD, ONLY : GET_TS_CHEM
!
! !INPUT PARAMETERS: 
!
      INTEGER, INTENT(IN)  :: JDAY              ! Day of year (0-365)
!
! !OUTPUT PARAMETERS:
!
      REAL*8,  INTENT(OUT) :: SUNCOS(MAXIJ)     ! cos(SZA) @ current time
      REAL*8,  INTENT(OUT) :: SUNCOS_MID(MAXIJ) ! cos(SZA) @ midpt of chem step
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
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!      
      INTEGER            :: FACTOR, I,        IJLOOP,  J
      INTEGER            :: HOUR,   MINUTE,   TS_SUN
      REAL*8             :: AHR,    GMT_MID,  R,       TIMLOC
      REAL*8             :: DEC,    C_DEC,    S_DEC
      REAL*8             :: YMID_R, C_YMID_R, S_YMID_R
!
! !DEFINED PARAMETERS:
!
      ! Coefficients for solar declination angle
      REAL*8,  PARAMETER :: A0 = 0.006918d0
      REAL*8,  PARAMETER :: A1 = 0.399912d0
      REAL*8,  PARAMETER :: A2 = 0.006758d0
      REAL*8,  PARAMETER :: A3 = 0.002697d0
      REAL*8,  PARAMETER :: B1 = 0.070257d0
      REAL*8,  PARAMETER :: B2 = 0.000907d0
      REAL*8,  PARAMETER :: B3 = 0.000148d0

      !=================================================================
      ! Initialization   
      !=================================================================

      ! Quantities for central chemistry timestep
      TS_SUN   = GET_TS_CHEM()                         ! Chemistry interval
      HOUR     = GET_HOUR()                            ! Current hour
      MINUTE   = GET_MINUTE()                          ! Current minutes
      FACTOR   = MINUTE / TS_SUN                       ! Multiplying factor

      ! GMT at the midpoint of the chemistry timestep
      GMT_MID  = ( DBLE( GET_HOUR()      )        )    ! Start of hour
     &         + ( DBLE( TS_SUN * FACTOR ) / 60d0 )    ! + steps already done
     &         + ( DBLE( TS_SUN / 2      ) / 60d0 )    ! + 1/2 chem timestep

      ! Path length of earth's orbit traversed since Jan 1 [radians]
      R        = ( 2d0 * PI / 365d0 ) * DBLE( JDAY - 1 ) 

      ! Solar declination angle (low precision formula) [radians]
      DEC      = A0 - A1*COS(     R ) + B1*SIN(     R )
     &              - A2*COS( 2d0*R ) + B2*SIN( 2d0*R )
     &              - A3*COS( 3d0*R ) + B3*SIN( 3d0*R )

      ! Pre-compute sin & cos of DEC outside of DO loops (for efficiency)
      S_DEC    = SIN( DEC )
      C_DEC    = COS( DEC )

      !=================================================================
      ! Compute cosine of solar zenith angle
      !=================================================================

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, YMID_R, S_YMID_R, C_YMID_R, IJLOOP, TIMLOC, AHR )

      ! Loop over latitude
      DO J = 1, JJPAR

         ! Latitude of grid box [radians]
         YMID_R   = GET_YMID_R( J )

         ! Pre-compute sin & cos of DEC outside of I loop (for efficiency)
         S_YMID_R = SIN( YMID_R )
         C_YMID_R = COS( YMID_R )

         ! Loop over longitude
         DO I = 1, IIPAR

            ! 1-D grid box index
            IJLOOP             = ( (J-1) * IIPAR ) + I

            !===========================================================
            ! Compute cosine of SZA at the current GMT time
            !===========================================================

            ! Local time at box (I,J) [hours]
            TIMLOC             = GET_LOCALTIME( I )

            ! Hour angle at box (I,J) [radians]
            AHR                = ABS( TIMLOC - 12d0 ) * 15d0 * PI_180
            
            ! Cosine of solar zenith angle at box (I,J) [unitless]
            SUNCOS(IJLOOP)     = ( S_YMID_R * S_DEC              ) 
     &                         + ( C_YMID_R * C_DEC * COS( AHR ) )

            !===========================================================
            ! Compute cosine(SZA) at the midpoint of the chem timestep
            ! Required for photolysis, chemistry, emissions, drydep
            !===========================================================

            ! Local time [hours] at box (I,J) at 30 minutes after the GMT hour
            TIMLOC             = GET_LOCALTIME( I, GMT=GMT_MID )

            ! Hour angle at box (I,J) [radians]
            AHR                = ABS( TIMLOC - 12d0 ) * 15d0 * PI_180
            
            ! Corresponding cosine( SZA ) at box (I,J) [unitless]
            SUNCOS_MID(IJLOOP) = ( S_YMID_R * S_DEC              )
     &                         + ( C_YMID_R * C_DEC * COS( AHR ) )

         ENDDO
      ENDDO
!$OMP END PARALLEL DO

      END SUBROUTINE COSSZA
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: convert_units
!
! !DESCRIPTION: Subroutine CONVERT\_UNITS converts the units of STT from [kg] 
!  to [v/v] mixing ratio, or vice versa.  
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE CONVERT_UNITS( IFLAG, N_TRACERS, TCVV, AD, STT ) 
!
! !USES:
!
      USE ERROR_MOD,  ONLY : ERROR_STOP
!
! !INPUT PARAMETERS: 
!

      ! =1 then convert from [kg ] --> [v/v]
      ! =2 then convert from [v/v] --> [kg ]
      INTEGER, INTENT(IN)    :: IFLAG

      ! Number of tracers
      INTEGER, INTENT(IN)    :: N_TRACERS 

      ! Array containing [Air MW / Tracer MW] for tracers
      REAL*8,  INTENT(IN)    :: TCVV(N_TRACERS)

      ! Array containing grid box air masses
      REAL*8,  INTENT(IN)    :: AD(IIPAR,JJPAR,LLPAR)
!
! !OUTPUT PARAMETERS:
!
      !  Array containing tracer conc. [kg] or [v/v]
      REAL*8,  INTENT(INOUT) :: STT(IIPAR,JJPAR,LLPAR,N_TRACERS)
!
! !REVISION HISTORY: 
!  15 Jun 1998 - R. Yantosca - Initial version
!  (1 ) CONVERT_UNITS is written in Fixed-Form Fortran 90.
!  (2 ) Cosmetic changes, updated comments (bmy, 4/19/00)
!  (3 ) Now use SELECT CASE statement.  Also added parallel DO-loops
!        with the new Open-MP compiler directives. (bmy, 4/27/00)
!  (4 ) Bundled into "dao_mod.f".  Now pass NTRACE, TCVV, AD, STT as args.
!        Now use explicit DO-loops for I-J-L w/in parallel loops.  Updated
!        comments, cosmetic changes. (bmy, 3/29/02)
!  (5 ) Removed obsolete, commented-out code.  Also now use F90 intrinsic
!        REPEAT to write a line of "="'s to the screen. (bmy, 6/25/02)
!  (6 ) Updated comments.  Now reference ERROR_STOP from "error_mod.f" 
!        (bmy, 10/15/02)
!  (7 ) Renamed NTRACE to N_TRACERS for consistency (bmy, 7/19/04)
!  16 Aug 2010 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER :: I, J, L, N

      !=================================================================
      ! CONVERT_UNITS begins here!
      !
      ! Most of the GEOS-CHEM subroutines require the tracer array 
      ! STT to be in units of [kg].  However, the cloud convection 
      ! (NFCLDMX), boundaryk layer mixing (TURBDAY), diffusion (DIFFUSE), 
      ! and transport (TPCORE) routines require STT to be in volume 
      ! mixing ratio [v/v].  
      !
      ! Therefore, before calling NFCLDMX, TURBDAY, DIFFUSE, or TPCORE, 
      ! call subroutine CONVERT_UNITS to convert STT from [kg] to [v/v].  
      !
      ! Also call CONVERT_UNITS again after calling NFCLDMX, TURBDAY, 
      ! DIFFUSE, or TPCORE to convert back from [v/v] to [kg].  
      !=================================================================
      SELECT CASE ( IFLAG )

         !==============================================================
         ! IFLAG = 1: Convert from [kg] -> [v/v] 
         !
         !  The conversion is as follows:
         !
         !   kg tracer(N)       1        Air mol wt     
         !   -----------  * -------- *  -------------   
         !        1          kg air     tracer mol wt   
         !
         !       moles tracer     volume tracer
         !   =   ------------  =  -------------
         !        moles air        volume air
         !
         ! Since the volume of a gas depends on the number of moles.
         ! Therefore, with:
         !
         !  TCMASS(N) = mol. wt. of tracer (AMU)
         !  TCVV(N)   = 28.97 / TCMASS(N)
         !            = mol. wt. of air (AMU) / mol. wt. of tracer (AMU)
         !  AD(I,J,L) = mass of air (kg) in grid box (I,J,L)
         !     
         ! the conversion is:
         ! 
         !  STT(I,J,L,N) [kg] * TCVV(N) / AD(I,J,L) = STT(I,J,L,N) [v/v]
         !==============================================================
         CASE ( 1 )

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED ) 
!$OMP+PRIVATE( I, J, L, N ) 
            DO N = 1, N_TRACERS
            DO L = 1, LLPAR
            DO J = 1, JJPAR
            DO I = 1, IIPAR
               STT(I,J,L,N) = STT(I,J,L,N) * TCVV(N) / AD(I,J,L)
            ENDDO
            ENDDO
            ENDDO
            ENDDO
!$OMP END PARALLEL DO

         !==============================================================
         ! IFLAG = 2: Convert from [v/v] -> [kg] 
         !
         ! From the above discussion, the reverse unit conversion 
         ! is given by:
         !
         !  STT(I,J,L,N) [v/v] * AD(I,J,L) / TCVV(N) = STT(I,J,L,N) [kg]
         !==============================================================
         CASE ( 2 )

!$OMP PARALLEL DO 
!$OMP+DEFAULT( SHARED ) 
!$OMP+PRIVATE( I, J, L, N )
            DO N = 1, N_TRACERS
            DO L = 1, LLPAR
            DO J = 1, JJPAR
            DO I = 1, IIPAR
               STT(I,J,L,N) = STT(I,J,L,N) * AD(I,J,L) / TCVV(N)
            ENDDO     
            ENDDO
            ENDDO
            ENDDO
!$OMP END PARALLEL DO

         !==============================================================
         ! Otherwise halt with an error message
         !==============================================================
         CASE DEFAULT
            CALL ERROR_STOP( 'Invalid IFLAG value (must be 1 or 2)!', 
     &                       'CONVERT_UNITS (dao_mod.f)' )
      END SELECT

      ! Return to calling program
      END SUBROUTINE CONVERT_UNITS
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: copy_i6_fields
!
! !DESCRIPTION: Subroutine COPY\_I6\_FIELDS copies the I-6 fields at the end 
!  of a 6-hr timestep.  The I-6 fields at the end of a given 6-hr timestep 
!  become the fields at the beginning of the next 6-hr timestep. 
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE COPY_I6_FIELDS
! 
! !REVISION HISTORY: 
!  13 Apr 2004 - R. Yantosca - Initial version
!  (1 ) Added parallel DO-loops (bmy, 4/13/04)
!  (2 ) Remove support for GEOS-1 and GEOS-STRAT met fields (bmy, 8/4/06)
!  (3 ) Added TROPP (phs 11/10/06)
!  (4 ) Don't copy TROPP2 to TROPP1 for GEOS-5 (bmy, 1/17/07) 
!  16 Aug 2010 - R. Yantosca - Added ProTeX headers
!  20 Aug 2010 - R. Yantosca - Rewrite #if block for clarity
!  20 Aug 2010 - R. Yantosca - Added #if block for MERRA met fields
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER :: I, J, L

      !=================================================================
      ! COPY_I6_FIELDS begins here!
      !=================================================================

#if   defined( GCAP )

      !-----------------
      ! GCAP met
      !-----------------
      PS1    = PS2          ! I6 surface pressure    [hPa]
      TROPP1 = TROPP2       ! I6 tropopause pressure [hPa]

#elif defined( GEOS_3 )

      !-----------------
      ! GEOS-3 met
      !-----------------
      PS1    = PS2          ! I6 surface pressure    [hPa]
      TROPP1 = TROPP2       ! I6 tropopause pressure [hPa]
      ALBD1  = ALBD2        ! I6 surface albedo      [unitless]
      UWND1  = UWND2        ! I6 zonal wind          [m/s]
      VWND1  = VWND2        ! I6 meridional wind     [m/s]
      TMPU1  = TMPU2        ! I6 temperature         [K]
      SPHU1  = SPHU2        ! I6 specific humidity   [g/kg]

#elif defined( GEOS_4 )

      !-----------------
      ! GEOS-4 met
      !-----------------
      PS1    = PS2          ! I6 surface pressure    [hPa]
      TROPP1 = TROPP2       ! I6 tropopause pressure [hPa]

#elif defined( GEOS_5 )

      !-----------------
      ! GEOS-5 met
      !-----------------
      PS1    = PS2          ! I6 surface pressure    [hPa] 

#elif defined( MERRA )

      !-----------------
      ! MERRA met
      !-----------------
      PS1    = PS2          ! I6 surface pressure    [hPa]
      RH1    = RH2          ! I6 relative humidity   [%]

#endif

      END SUBROUTINE COPY_I6_FIELDS
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_dao
!
! !DESCRIPTION: Subroutine INIT\_DAO allocates memory for all allocatable 
!  module arrays. 
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE INIT_DAO
!
! !USES:
!
      USE ERROR_MOD,   ONLY : ALLOC_ERR
      USE LOGICAL_MOD, ONLY : LWETD, LDRYD, LCHEM
      USE TRACER_MOD,  ONLY : ITS_AN_AEROSOL_SIM, ITS_A_FULLCHEM_SIM
!
! !REVISION HISTORY: 
!  26 Jun 2000 - R. Yantosca - Initial version
!  (1 ) Now allocate AVGW for either NSRCX == 3 or NSRCX == 5 (bmy, 9/24/01)
!  (2 ) Removed obsolete code from 9/01 (bmy, 10/23/01)
!  (3 ) Add PSC2 array for TPCORE mixing ratio fix.   (bmy, 3/27/02)
!  (4 ) Elimintated PS, PSC arrays for floating-pressure fix.
!        (dsa, bdf, bmy, 8/20/02)
!  (5 ) Added AD, BXHEIGHT, T to "dao_mod.f" as allocatable arrays, to remove
!        historical baggage and centralize variables.  Also remove GEOS_2 
!        flag from C-preprocessor statements.  Also allocate RH array
!        but only if we are doing a sulfate simulation.  Now references
!        ALLOC_ERR from "error_mod.f" (bmy, 10/15/02)
!  (6 ) Now allocate PHIS array (bmy, 3/11/03)
!  (7 ) Now allocate SUNCOSB array for SMVGEAR II.  Also removed KZZ array,
!        that is now obsolete. (bdf, bmy, 4/28/03)
!  (8 ) Now order all arrays in alphabetical order.  Also added new fields
!        for GEOS-4/fvDAS: HKBETA, HKETA, ZMEU, ZMMD, ZMMU, TSKIN, PARDF,
!        and PARDR. (bmy, 6/25/03)
!  (9 ) Now allocate CLDFRC, RADLWG, RADSWG, SNOW arrays.  USTAR, CLDFRC,
!        and Z0 and RADSWG are now 2-D arrays. (bmy, 12/9/03)
!  (10) Allocate RADLWG and SNOW for both GEOS-3 & GEOS-4 (bmy, 4/2/04)
!  (11) Now reference inquiry functions from "tracer_mod.f".  Now reference
!        LWETD, LDRYD, LCHEM from "logical_mod.f".  Now allocate RH regardless
!        of simulation. (bmy, 7/20/04)
!  (12) Now also allocate AVGW for offline aerosol simulations (bmy, 9/27/04)
!  (13) Now modified for GCAP met fields.  Removed references to CO-OH param 
!        simulation.  Now allocate AVGW only for fullchem or offline aerosol
!        simulations. (bmy, 6/24/05)
!  (14) Now allocate SNOW and GWETTOP for GCAP (bmy, 8/17/05)
!  (15) Now also add TSKIN for GEOS-3 (bmy, 10/20/05)
!  (16) Remove support for GEOS-1 and GEOS-STRAT met fields (bmy, 8/4/06)
!  (17) Reorganized for GEOS-5 met fields (bmy, 1/17/07)
!  (18) Bug fix: should be CMFMC=0 after allocating CMFMC (jaf, bmy, 6/11/08)
!  (19) Remove obsolete SUNCOSB array (bmy, 4/28/10)
!  16 Aug 2010 - R. Yantosca - Added ProTeX headers
!  18 Aug 2010 - R. Yantosca - Now allocate met fields for MERRA
!  20 Aug 2010 - R. Yantosca - Bug fix, now allocate REEVAPCN
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER :: AS

      !=================================================================
      ! INIT_DAO begins here!
      !=================================================================

      !-----------------------------------------------------------------
      ! Allocate met field arrays that are used for all met fields
      !-----------------------------------------------------------------

      ! Air mass
      ALLOCATE( AD( IIPAR, JJPAR, LLPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'AD' )
      AD = 0d0

      ! Air density
      ALLOCATE( AIRDEN( LLPAR, IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'AIRDEN' )
      AIRDEN = 0d0

      ! Air volume
      ALLOCATE( AIRVOL( IIPAR, JJPAR, LLPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'AIRVOL' )
      AIRVOL = 0d0
      
      ! Surface albedo
      ALLOCATE( ALBD( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'ALBD' )
      ALBD = 0d0

      ! AVGW (mixing ratio of H2O) is only used for NOx-Ox-HC or aerosol sims
      IF ( ITS_A_FULLCHEM_SIM() .or. ITS_AN_AEROSOL_SIM() ) THEN 
         ALLOCATE( AVGW( IIPAR, JJPAR, LLPAR ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'AVGW' )
         AVGW = 0d0
      ENDIF

      ! Grid box height
      ALLOCATE( BXHEIGHT( IIPAR, JJPAR, LLPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'BXHEIGHT' )
      BXHEIGHT = 0d0

      ! 3-D Cloud fraction
      ALLOCATE( CLDF( LLPAR, IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'CLDF' )
      CLDF = 0d0

      ! 2-D column cloud fraction
      ALLOCATE( CLDFRC( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'CLDFRC' )
      CLDFRC = 0d0

      ! Cloud top level
      ALLOCATE( CLDTOPS( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'CLDTOPS' )
      CLDTOPS = 0

      ! Pressure difference between levels
      ALLOCATE( DELP( LLPAR, IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'DELP' )
      DELP = 0d0

      ! Top soil wetness
      ALLOCATE( GWETTOP( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'GWETTOP' )
      GWETTOP = 0d0

      ! Sensible heat flux
      ALLOCATE( HFLUX( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'HFLUX' )
      HFLUX = 0d0

      ! Tendency of specific humidity
      ALLOCATE( MOISTQ( LLPAR, IIPAR, JJPAR ), STAT=AS ) 
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'MOISTQ' )
      MOISTQ = 0d0

      ! Optical depth
      ALLOCATE( OPTDEP( LLPAR, IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'OPTDEP' )
      OPTDEP = 0d0

      ! Optical depth
      ALLOCATE( OPTD( LLPAR, IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'OPTD' )
      OPTD = 0d0

      ! Diffuse PAR
      ALLOCATE( PARDF( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'PARDF' )
      PARDF = 0d0

      ! Direct PAR
      ALLOCATE( PARDR( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'PARDR' )
      PARDR = 0d0

      ! Mixed layer depth
      ALLOCATE( PBL( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'PBL' )
      PBL = 0d0

      ! Surface geopotential height
      ALLOCATE( PHIS( IIPAR, JJPAR ), STAT=AS )  
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'PHIS' )
      PHIS = 0d0

      ! Total precip at ground
      ALLOCATE( PREACC( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'PREACC' )
      PREACC = 0d0

      ! Convective precip at ground
      ALLOCATE( PRECON( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'PRECON' )
      PRECON = 0d0
   
      ! Pressure at beginning of 6hr timestep
      ALLOCATE( PS1( IIPAR, JJPAR ), STAT=AS )  
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'PS1' )
      PS1 = 0d0

      ! Pressure at end of 6hr timestep
      ALLOCATE( PS2( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'PS2' )
      PS2 = 0d0

      ! Pressure at end of dynamic timestep
      ALLOCATE( PSC2( IIPAR, JJPAR ), STAT=AS ) 
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'PSC2' )
      PSC2 = 0d0

      ! Longwave rad at ground
      ! NOTE: this is a net radiation for GEOS-5 (LWGNET)
      ALLOCATE( RADLWG( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'RADLWG' )
      RADLWG = 0d0

      ! Shortwave rad at ground
      ! NOTE: this is a net radiation for GEOS-5 (SWGNET)
      ALLOCATE( RADSWG( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'RADSWG' )
      RADSWG = 0d0

      ! Relative humidity
      ALLOCATE( RH( IIPAR, JJPAR, LLPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'RH' )
      RH = 0d0

      ! Sea level pressure
      ALLOCATE( SLP( IIPAR, JJPAR ), STAT=AS ) 
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'SLP' )
      SLP = 0d0      

      ! Specific humidity
      ALLOCATE( SPHU( IIPAR, JJPAR, LLPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'SPHU' )
      SPHU = 0d0

      ! Cosine of solar zenith angle
      ALLOCATE( SUNCOS( MAXIJ ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'SUNCOS' )
      SUNCOS = 0d0

      ! Cosine of solar zenith angle (computed @ midpt of chem timestep)
      ALLOCATE( SUNCOS_MID( MAXIJ ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'SUNCOS_MID' )
      SUNCOS_MID = 0d0

      ! Temperature
      ALLOCATE( T( IIPAR, JJPAR, LLPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'T' )
      T = 0d0

      ! Temperature on full vertical grid 
      ALLOCATE( T_FULLGRID( IIPAR, JJPAR, LGLOB ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'T_FULLGRID' )
      T_FULLGRID = 0d0

      ! Number density of air on full vertical grid
      ALLOCATE( AIRDEN_FULLGRID( LGLOB, IIPAR, JJPAR ), STAT=AS ) 
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'AIRDEN_FULLGRID' )
      AIRDEN_FULLGRID = 0d0

      ! Tropopause pressure
      ALLOCATE( TROPP( IIPAR, JJPAR ), STAT=AS ) 
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'TROPP' )
      TROPP = 0d0

      ! Surface temperature
      ALLOCATE( TS( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'TS' )
      TS = 0d0

      ! Skin (aka ground) temperature
      ALLOCATE( TSKIN( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'TSKIN' )
      TSKIN = 0d0

      ! 10m U-winds
      ALLOCATE( U10M( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'U10M' )
      U10M = 0d0

      ! Friction velocity
      ALLOCATE( USTAR( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'USTAR' )
      USTAR = 0d0

      ! U-wind
      ALLOCATE( UWND( IIPAR, JJPAR, LLPAR), STAT=AS ) 
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'UWND' )
      UWND = 0d0

      ! 10m V-wind
      ALLOCATE( V10M( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'V10M' )
      V10M = 0d0

      ! V-wind
      ALLOCATE( VWND( IIPAR, JJPAR, LLPAR ), STAT=AS ) 
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'VWND' )
      VWND = 0d0

      ! Roughness height
      ALLOCATE( Z0( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'Z0' )
      Z0 = 0d0

      !-----------------------------------------------------------------
      ! Allocate proper array for land/water/ice flags 
      !-----------------------------------------------------------------

#if   defined( GCAP )

      ! Land/water flags have to be REAL*8 for GCAP
      ALLOCATE( LWI_GISS( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'LWI_GISS' )
      LWI_GISS = 0d0

#else

      ! Land/water flags have to be INTEGER for GEOS
      ALLOCATE( LWI( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'LWI' )
      LWI = 0

#endif

#if   !defined( GEOS_5 )

      !-----------------------------------------------------------------
      ! Allocate met field arrays for everything EXCEPT GEOS-5
      !-----------------------------------------------------------------

      ! Snow depth for all other met fields
      ALLOCATE( SNOW( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'SNOW' )
      SNOW = 0d0

      ! TROPP1 is only defined for GEOS-3, GEOS-4, or GCAP
      ALLOCATE( TROPP1( IIPAR, JJPAR ), STAT=AS ) 
      IF ( AS /= 0 ) C ALL ALLOC_ERR( 'TROPP1' )
      TROPP1 = 0d0

      ! TROPP2 is only defined for GEOS-3, GEOS-4, or GCAP
      ALLOCATE( TROPP2( IIPAR, JJPAR ), STAT=AS ) 
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'TROPP2' )
      TROPP2 = 0d0

#endif

#if   defined( GEOS_3 )

      !-----------------------------------------------------------------
      ! Allocate met field arrays that are only used for GEOS-3
      !-----------------------------------------------------------------

      ! Albedo at start of 6-hr interval
      ALLOCATE( ALBD1( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'ALBD1' )
      ALBD1 = 0d0

      ! Albedo at end of 6-hr interval
      ALLOCATE( ALBD2( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'ALBD2' )
      ALBD2 = 0d0

      ! GEOS-3 cloud mass flux
      ALLOCATE( CLDMAS( IIPAR, JJPAR, LLPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'CLDMAS' )
      CLDMAS = 0d0

      ! GEOS-3 detrainment 
      ALLOCATE( DTRAIN( IIPAR, JJPAR, LLPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'DTRAIN' )
      DTRAIN = 0d0

      ! Specific humidity at start of 6-hr interval
      ALLOCATE( SPHU1( IIPAR, JJPAR, LLPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'SPHU1' )
      SPHU1 = 0d0

      ! Specific humidity at end of 6-hr interval
      ALLOCATE( SPHU2( IIPAR, JJPAR, LLPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'SPHU2' )
      SPHU2 = 0d0

      ! Temperature at start of 6-hr interval
      ALLOCATE( TMPU1( IIPAR, JJPAR, LLPAR ), STAT=AS ) 
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'TMPU1' )
      TMPU1 = 0d0

      ! Temperature at end of 6-hr interval
      ALLOCATE( TMPU2( IIPAR, JJPAR, LLPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'TMPU2' )
      TMPU2 = 0d0

      ! U-wind at start of 6-hr interval
      ALLOCATE( UWND1( IIPAR, JJPAR, LLPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'UWND1' )
      UWND1 = 0d0

      ! U-wind at end of 6-hr interval
      ALLOCATE( UWND2( IIPAR, JJPAR, LLPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'UWND2' )
      UWND2 = 0d0

      ! V-wind at start of 6-hr interval
      ALLOCATE( VWND1( IIPAR, JJPAR, LLPAR ), STAT=AS ) 
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'VWND1' )
      VWND1 = 0d0

      ! V-wind at end of 6-hr interval
      ALLOCATE( VWND2( IIPAR, JJPAR, LLPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'VWND2' )
      VWND2 = 0d0

#elif defined( GEOS_4 )

      !-----------------------------------------------------------------
      ! Allocate met field arrays that are only used for GEOS-4
      !-----------------------------------------------------------------

      ! Hack convection overshoot parameter
      ALLOCATE( HKBETA( IIPAR, JJPAR, LLPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'HKBETA' )
      HKBETA = 0d0

      ! Hack convection mass flux
      ALLOCATE( HKETA( IIPAR, JJPAR, LLPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'HKETA' )
      HKETA = 0d0

      ! Z&M updraft entrainment flux
      ALLOCATE( ZMEU( IIPAR, JJPAR, LLPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'ZMEU' )
      ZMEU = 0d0

      ! Z&M downdraft mass flux
      ALLOCATE( ZMMD( IIPAR, JJPAR, LLPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'ZMMD' )
      ZMMD = 0d0

      ! Z&M updraft mass flux
      ALLOCATE( ZMMU( IIPAR, JJPAR, LLPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'ZMMU' )
      ZMMU = 0d0

#elif defined( GEOS_5 ) || defined( MERRA )

      !-----------------------------------------------------------------
      ! Allocate met field arrays that are used for GEOS-5 or MERRA
      !-----------------------------------------------------------------

      ! GEOS-5 cloud mass flux 
      ALLOCATE( CMFMC( IIPAR, JJPAR, LLPAR+1 ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'CMFMC' )
      CMFMC = 0d0

      ! GEOS-5 tendency of ice in moist processes
      ALLOCATE( DQIDTMST( IIPAR, JJPAR, LLPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'DQIDTMST' )
      DQIDTMST = 0d0

      ! GEOS-5 tendency of ice in moist processes
      ALLOCATE( DQLDTMST( IIPAR, JJPAR, LLPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'DQLDTMST' )
      DQLDTMST = 0d0

      ! GEOS-5 tendency of ice in moist processes
      ALLOCATE( DQVDTMST( IIPAR, JJPAR, LLPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'DQVDTMST' )
      DQVDTMST = 0d0

      ! GEOS-5 detrainment 
      ALLOCATE( DTRAIN( IIPAR, JJPAR, LLPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'DTRAIN' )
      DTRAIN = 0d0

      ! Latent heat flux
      ALLOCATE( EFLUX( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'EFLUX' )
      EFLUX = 0d0

      ! GEOS-5 evapotranspiration flux
      ALLOCATE( EVAP( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'EVAP' )
      EVAP = 0d0

      ! Fraction of grid box that is land
      ALLOCATE( FRLAND( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'FRLAND' )
      FRLAND = 0d0

      ! Fraction of grid box that is lakes
      ALLOCATE( FRLAKE( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'FRLAKE' )
      FRLAKE = 0d0

      ! Fraction of grid box that is ocean
      ALLOCATE( FROCEAN( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'FROCEAN' )
      FROCEAN = 0d0

      ! Fraction of grid box that is land ice
      ALLOCATE( FRLANDIC( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'FRLANDIC' )
      FRLANDIC = 0d0

      ! GEOS-5 greenness index
      ALLOCATE( GRN( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'GRN' )
      GRN = 0d0

      ! GEOS-5 root soil moisture
      ALLOCATE( GWETROOT( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'GWETROOT' )
      GWETROOT = 0d0

      ! GEOS-5 root soil moisture
      ALLOCATE( LAI( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'LAI' )
      LAI = 0d0

      ! GEOS-5 "snow" (i.e. frozen) precipitation
      ALLOCATE( PRECSNO( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'PRECSNO' )
      PRECSNO = 0d0

      ! GEOS-5 potential vorticity
      ALLOCATE( PV( IIPAR, JJPAR, LLPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'PV' )
      PV = 0d0

      ! GEOS-5 ice mixing ratio
      ALLOCATE( QI( IIPAR, JJPAR, LLPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'QI' )
      QI = 0d0

      ! GEOS-5 ice mixing ratio
      ALLOCATE( QL( IIPAR, JJPAR, LLPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'QL' )
      QL = 0d0

      ! GEOS-5 snow depth (H2O equiv)
      ALLOCATE( SNOMAS( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'SNOMAS' )
      SNOMAS = 0d0

      ! GEOS-5 snow depth (geometric)
      ALLOCATE( SNODP( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'SNODP' )
      SNODP = 0d0

      ! GEOS-5 ice path optical depth
      ALLOCATE( TAUCLI( IIPAR, JJPAR, LLPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'TAUCLI' )
      TAUCLI = 0d0

      ! GEOS-5 water path optical depth
      ALLOCATE( TAUCLW( IIPAR, JJPAR, LLPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'TAUCLW' )
      TAUCLW = 0d0

      ! GEOS-5 total column ozone at beginning of 6hr timestep
      ALLOCATE( TO31( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'TO31' )
      TO31 = 0d0

      ! GEOS-5 total column ozone at end of 6hr timestep
      ALLOCATE( TO32( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'TO32' )
      TO32 = 0d0

      ! GEOS-5 total column ozone
      ALLOCATE( TO3( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'TO3' )
      TO3 = 0d0

      ! GEOS-5 total trop column ozone
      ALLOCATE( TTO3( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'TTO3' )
      TTO3 = 0d0

#if   defined( GEOS_5 )

      !-----------------------------------------------------------------
      ! Allocate met field arrays that are used ONLY for GEOS-5
      !-----------------------------------------------------------------

      ! GEOS-5 convective rain production
      ALLOCATE( DQRCON( IIPAR, JJPAR, LLPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'DQRCON' )
      DQRCON = 0d0

      ! GEOS-5 convective rain production
      ALLOCATE( DQRLSC( IIPAR, JJPAR, LLPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'DQRLSC' )
      DQRLSC = 0d0

#elif defined( MERRA )

      !-----------------------------------------------------------------
      ! Allocate met field arrays that are used ONLY for MERRA
      !-----------------------------------------------------------------

      ! Precipitation production rate -- convective
      ALLOCATE( DQRCU( IIPAR, JJPAR, LLPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'DQRCU' )
      DQRCU = 0d0

      ! Precipitation production rate -- LS + anvil
      ALLOCATE( DQRLSAN( IIPAR, JJPAR, LLPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'DQRLSAN' )
      DQRLSAN = 0d0

      ! Fraction of grid box that is land ice
      ALLOCATE( FRSEAICE( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'FRSEAICE' )
      FRSEAICE = 0d0

      ! Fraction of grid box that is snow
      ALLOCATE( FRSNO( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'FRSNO' )
      FRSNO = 0d0

      ! Downward flux of ice precipitation -- convective 
      ALLOCATE( PFICU( IIPAR, JJPAR, LLPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'PFICU' )
      PFICU = 0d0

      ! Downward flux of ice precipitation -- LS + anvil
      ALLOCATE( PFILSAN( IIPAR, JJPAR, LLPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'PFILSAN' )
      PFILSAN = 0d0

      ! Downward flux of liquid precipitation -- convective 
      ALLOCATE( PFLCU( IIPAR, JJPAR, LLPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'PFLCU' )
      PFLCU = 0d0

      ! Downward flux of liquid precipitation -- LS + anvil
      ALLOCATE( PFLLSAN( IIPAR, JJPAR, LLPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'PFILSAN' )
      PFLLSAN = 0d0

      ! Surface precipitation flux from anvils
      ALLOCATE( PREANV( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'PREANV' )
      PREANV = 0d0

      ! Surface precipitation flux from large-scale 
      ALLOCATE( PRELSC( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'PRELSC' )
      PRELSC = 0d0

      ! RH at current 3-hr timestep
      ALLOCATE( RH1( IIPAR, JJPAR, LLPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'RH1' )
      RH1 = 0d0

      ! RH at next 3-hr timestep
      ALLOCATE( RH2( IIPAR, JJPAR, LLPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'RH2' )
      RH2 = 0d0

      ! Evaporation of precipitating convective condensate 	
      ALLOCATE( REEVAPCN( IIPAR, JJPAR, LLPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'REEVAPCN' )
      REEVAPCN = 0d0

      ! Evaporation of precipitating LS & anvil condensate 	
      ALLOCATE( REEVAPLS( IIPAR, JJPAR, LLPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'REEVAPLS' )
      REEVAPLS = 0d0

      ! 0-10% sea ice coverage
      ALLOCATE( SEAICE00( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'SEAICE00' )
      SEAICE00 = 0d0

      ! 10-20% sea ice coverage
      ALLOCATE( SEAICE10( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'SEAICE10' )
      SEAICE10 = 0d0

      ! 20-30% sea ice coverage
      ALLOCATE( SEAICE20( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'SEAICE20' )
      SEAICE20 = 0d0

      ! 30-40% sea ice coverage
      ALLOCATE( SEAICE30( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'SEAICE30' )
      SEAICE30 = 0d0

      ! 40-50% sea ice coverage
      ALLOCATE( SEAICE40( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'SEAICE40' )
      SEAICE40 = 0d0

      ! 50-60% sea ice coverage
      ALLOCATE( SEAICE50( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'SEAICE50' )
      SEAICE50 = 0d0

      ! 60-70% sea ice coverage
      ALLOCATE( SEAICE60( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'SEAICE60' )
      SEAICE60 = 0d0

      ! 70-80% sea ice coverage
      ALLOCATE( SEAICE70( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'SEAICE70' )
      SEAICE70 = 0d0

      ! 80-90 sea ice coverage
      ALLOCATE( SEAICE80( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'SEAICE80' )
      SEAICE80 = 0d0

      ! 90-100% sea ice coverage
      ALLOCATE( SEAICE90( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'SEAICE90' )
      SEAICE90 = 0d0

#endif 

#elif defined( GCAP )

      !-----------------------------------------------------------------
      ! Allocate met field arrays that are only used for GCAP
      !-----------------------------------------------------------------

      ! DTRAINE is only defined for GCAP
      ALLOCATE( DETRAINE( IIPAR, JJPAR, LLPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'DETRAINE' )
      DETRAINE = 0d0  

      ! DETRAINN is only defined for GCAP
      ALLOCATE( DETRAINN( IIPAR, JJPAR, LLPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'DETRAINN' )
      DETRAINN = 0d0  

      ! DNDE is only defined for GCAP
      ALLOCATE( DNDE( IIPAR, JJPAR, LLPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'DNDE' )
      DNDE = 0d0  

      ! DNDN is only defined for GCAP
      ALLOCATE( DNDN( IIPAR, JJPAR, LLPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'DNDN' )
      DNDN = 0d0  

      ! ENTRAIN is only defined for GCAP
      ALLOCATE( ENTRAIN( IIPAR, JJPAR, LLPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'ENTRAIN' )
      ENTRAIN = 0d0  

      ! MOLENGTH is only defined for GCAP
      ALLOCATE( MOLENGTH( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'MOLENGTH' )
      MOLENGTH = 0d0

      ! OICE is only defined for GCAP
      ALLOCATE( OICE( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'OICE' )
      OICE = 0d0

      ! SNICE is only defined for GCAP
      ALLOCATE( SNICE( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'SNICE' )
      SNICE = 0d0

      ! UPDE is only defined for GCAP
      ALLOCATE( UPDE( IIPAR, JJPAR, LLPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'UPDE' )
      UPDE = 0d0  

      ! UPDN is only defined for GCAP
      ALLOCATE( UPDN( IIPAR, JJPAR, LLPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'UPDN' )
      UPDN = 0d0  

#endif

      ! Return to calling program
      END SUBROUTINE INIT_DAO
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: cleanup_dao
!
! !DESCRIPTION: Subroutine CLEANUP\_DAO deallocates all met field arrays. 
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE CLEANUP_DAO
! 
! !REVISION HISTORY: 
!  26 Jun 2000 - R. Yantosca - Initial version
!  (1 ) Now deallocate SLP met field for GEOS-3 (bmy, 10/10/00)
!  (2 ) Now deallocate OPTDEP met field for GEOS-3 (bmy, 8/15/01)
!  (3 ) Now deallocate AVGW (bmy, 9/24/01)
!  (4 ) Remove TAUCLD deallocation -- it's obsolete (bmy, 10/23/01)
!  (5 ) Add call to deallocate PSC2 array (bmy, 3/27/02)
!  (6 ) Elimintated PS, PSC arrays for floating-pressure fix.
!        (dsa, bdf, bmy, 8/20/02)
!  (7 ) Now deallocate AD, BXHEIGHT, and T arrays (bmy, 9/18/02)
!  (8 ) Now deallocate PHIS array (bmy, 3/11/03)
!  (9 ) Now deallocate SUNCOSB array.  Remove reference to KZZ, since
!        that is now obsolete. (bmy, 4/28/03)
!  (10) Now list all arrays in order.  Now also deallocate new arrays
!        for GEOS-4/fvDAS. (bmy, 6/25/03)
!  (11) Now deallocate CLDFRC, RADLWG, RADSWG, SNOW arrays (bmy, 12/9/03)
!  (12) Now deallocate GCAP met fields (bmy, 5/25/05)
!  (13) Remove support for GEOS-1 and GEOS-STRAT met fields (bmy, 8/4/06)
!  (14) Deallocate additional arrays for GEOS-5 (bmy, 1/17/07)
!  (15) Remove obsolete SUNCOSB (bmy, 4/28/10)
!  16 Aug 2010 - R. Yantosca - Added ProTeX headers
!  18 Aug 2010 - R. Yantosca - Now deallocate MERRA met field arrays
!  05 Oct 2011 - R. Yantosca - Now deallocate SUNCOS_MID
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      !=================================================================
      ! CLEANUP_DAO begins here!
      !=================================================================
      IF ( ALLOCATED( AD       ) ) DEALLOCATE( AD       ) 
      IF ( ALLOCATED( AIRDEN   ) ) DEALLOCATE( AIRDEN   )
      IF ( ALLOCATED( AIRVOL   ) ) DEALLOCATE( AIRVOL   )
      IF ( ALLOCATED( ALBD1    ) ) DEALLOCATE( ALBD1    )
      IF ( ALLOCATED( ALBD2    ) ) DEALLOCATE( ALBD2    )
      IF ( ALLOCATED( ALBD     ) ) DEALLOCATE( ALBD     )
      IF ( ALLOCATED( AVGW     ) ) DEALLOCATE( AVGW     )
      IF ( ALLOCATED( BXHEIGHT ) ) DEALLOCATE( BXHEIGHT )
      IF ( ALLOCATED( CLDF     ) ) DEALLOCATE( CLDF     )
      IF ( ALLOCATED( CLDFRC   ) ) DEALLOCATE( CLDFRC   )
      IF ( ALLOCATED( CLDMAS   ) ) DEALLOCATE( CLDMAS   )
      IF ( ALLOCATED( CLDTOPS  ) ) DEALLOCATE( CLDTOPS  )
      IF ( ALLOCATED( CMFMC    ) ) DEALLOCATE( CMFMC    )
      IF ( ALLOCATED( DELP     ) ) DEALLOCATE( DELP     )
      IF ( ALLOCATED( DETRAINE ) ) DEALLOCATE( DETRAINE )
      IF ( ALLOCATED( DETRAINN ) ) DEALLOCATE( DETRAINN ) 
      IF ( ALLOCATED( DNDE     ) ) DEALLOCATE( DNDE     ) 
      IF ( ALLOCATED( DNDN     ) ) DEALLOCATE( DNDN     )
      IF ( ALLOCATED( DQIDTMST ) ) DEALLOCATE( DQIDTMST )
      IF ( ALLOCATED( DQLDTMST ) ) DEALLOCATE( DQLDTMST )
      IF ( ALLOCATED( DQRCON   ) ) DEALLOCATE( DQRCON   )
      IF ( ALLOCATED( DQRCU    ) ) DEALLOCATE( DQRCU    )
      IF ( ALLOCATED( DQRLSC   ) ) DEALLOCATE( DQRLSC   )
      IF ( ALLOCATED( DQRLSAN  ) ) DEALLOCATE( DQRLSAN  )
      IF ( ALLOCATED( DQVDTMST ) ) DEALLOCATE( DQVDTMST )
      IF ( ALLOCATED( DTRAIN   ) ) DEALLOCATE( DTRAIN   )
      IF ( ALLOCATED( ENTRAIN  ) ) DEALLOCATE( ENTRAIN  ) 
      IF ( ALLOCATED( EVAP     ) ) DEALLOCATE( EVAP     ) 
      IF ( ALLOCATED( FRLAND   ) ) DEALLOCATE( FRLAND   )
      IF ( ALLOCATED( FRLAKE   ) ) DEALLOCATE( FRLAKE   )
      IF ( ALLOCATED( FROCEAN  ) ) DEALLOCATE( FROCEAN  )
      IF ( ALLOCATED( FRLANDIC ) ) DEALLOCATE( FRLANDIC )
      IF ( ALLOCATED( FRSEAICE ) ) DEALLOCATE( FRSEAICE )
      IF ( ALLOCATED( FRSNO    ) ) DEALLOCATE( FRSNO    )
      IF ( ALLOCATED( GRN      ) ) DEALLOCATE( GRN      ) 
      IF ( ALLOCATED( GWETROOT ) ) DEALLOCATE( GWETROOT ) 
      IF ( ALLOCATED( GWETTOP  ) ) DEALLOCATE( GWETTOP  )
      IF ( ALLOCATED( HFLUX    ) ) DEALLOCATE( HFLUX    )
      IF ( ALLOCATED( HKBETA   ) ) DEALLOCATE( HKBETA   )
      IF ( ALLOCATED( HKETA    ) ) DEALLOCATE( HKETA    )
      IF ( ALLOCATED( LAI      ) ) DEALLOCATE( LAI      )      
      IF ( ALLOCATED( LWI      ) ) DEALLOCATE( LWI      )
      IF ( ALLOCATED( LWI_GISS ) ) DEALLOCATE( LWI_GISS ) 
      IF ( ALLOCATED( MFXC     ) ) DEALLOCATE( MFXC     ) 
      IF ( ALLOCATED( MFYC     ) ) DEALLOCATE( MFYC     ) 
      IF ( ALLOCATED( MFZ      ) ) DEALLOCATE( MFZ      ) 
      IF ( ALLOCATED( MOLENGTH ) ) DEALLOCATE( MOLENGTH ) 
      IF ( ALLOCATED( MOISTQ   ) ) DEALLOCATE( MOISTQ   )
      IF ( ALLOCATED( OICE     ) ) DEALLOCATE( OICE     )  
      IF ( ALLOCATED( OPTD     ) ) DEALLOCATE( OPTD     )
      IF ( ALLOCATED( OPTDEP   ) ) DEALLOCATE( OPTDEP   )
      IF ( ALLOCATED( PARDF    ) ) DEALLOCATE( PARDF    )
      IF ( ALLOCATED( PARDR    ) ) DEALLOCATE( PARDR    )
      IF ( ALLOCATED( PBL      ) ) DEALLOCATE( PBL      )
      IF ( ALLOCATED( PHIS     ) ) DEALLOCATE( PHIS     )
      IF ( ALLOCATED( PFICU    ) ) DEALLOCATE( PFICU    )
      IF ( ALLOCATED( PFILSAN  ) ) DEALLOCATE( PFILSAN  )
      IF ( ALLOCATED( PFLCU    ) ) DEALLOCATE( PFLCU    )
      IF ( ALLOCATED( PFLLSAN  ) ) DEALLOCATE( PFLLSAN  )
      IF ( ALLOCATED( PREACC   ) ) DEALLOCATE( PREACC   )
      IF ( ALLOCATED( PREANV   ) ) DEALLOCATE( PREANV   )
      IF ( ALLOCATED( PRECON   ) ) DEALLOCATE( PRECON   )
      IF ( ALLOCATED( PRELSC   ) ) DEALLOCATE( PRELSC   )
      IF ( ALLOCATED( PRECSNO  ) ) DEALLOCATE( PRECSNO  )
      IF ( ALLOCATED( PS1      ) ) DEALLOCATE( PS1      )
      IF ( ALLOCATED( PS2      ) ) DEALLOCATE( PS2      )
      IF ( ALLOCATED( PSC2     ) ) DEALLOCATE( PSC2     )
      IF ( ALLOCATED( PV       ) ) DEALLOCATE( PV       )
      IF ( ALLOCATED( QI       ) ) DEALLOCATE( QI       )
      IF ( ALLOCATED( QL       ) ) DEALLOCATE( QL       )
      IF ( ALLOCATED( RADLWG   ) ) DEALLOCATE( RADLWG   )
      IF ( ALLOCATED( RADSWG   ) ) DEALLOCATE( RADSWG   )
      IF ( ALLOCATED( REEVAPCN ) ) DEALLOCATE( REEVAPCN )
      IF ( ALLOCATED( REEVAPLS ) ) DEALLOCATE( REEVAPLS )
      IF ( ALLOCATED( RH1      ) ) DEALLOCATE( RH1      )
      IF ( ALLOCATED( RH2      ) ) DEALLOCATE( RH2      )
      IF ( ALLOCATED( RH       ) ) DEALLOCATE( RH       )
      IF ( ALLOCATED( SEAICE00 ) ) DEALLOCATE( SEAICE00 )
      IF ( ALLOCATED( SEAICE10 ) ) DEALLOCATE( SEAICE10 )
      IF ( ALLOCATED( SEAICE20 ) ) DEALLOCATE( SEAICE20 )
      IF ( ALLOCATED( SEAICE30 ) ) DEALLOCATE( SEAICE30 )
      IF ( ALLOCATED( SEAICE40 ) ) DEALLOCATE( SEAICE40 )
      IF ( ALLOCATED( SEAICE50 ) ) DEALLOCATE( SEAICE50 )
      IF ( ALLOCATED( SEAICE60 ) ) DEALLOCATE( SEAICE60 )
      IF ( ALLOCATED( SEAICE70 ) ) DEALLOCATE( SEAICE70 )
      IF ( ALLOCATED( SEAICE80 ) ) DEALLOCATE( SEAICE80 )
      IF ( ALLOCATED( SEAICE90 ) ) DEALLOCATE( SEAICE90 )
      IF ( ALLOCATED( SLP      ) ) DEALLOCATE( SLP      )
      IF ( ALLOCATED( SNICE    ) ) DEALLOCATE( SNICE    )
      IF ( ALLOCATED( SNODP    ) ) DEALLOCATE( SNODP    )
      IF ( ALLOCATED( SNOMAS   ) ) DEALLOCATE( SNOMAS   )
      IF ( ALLOCATED( SNOW     ) ) DEALLOCATE( SNOW     )
      IF ( ALLOCATED( SPHU1    ) ) DEALLOCATE( SPHU1    )
      IF ( ALLOCATED( SPHU2    ) ) DEALLOCATE( SPHU2    )
      IF ( ALLOCATED( SPHU     ) ) DEALLOCATE( SPHU     )
      IF ( ALLOCATED( SUNCOS   ) ) DEALLOCATE( SUNCOS   )
      IF ( ALLOCATED( SUNCOS_MID ) ) DEALLOCATE( SUNCOS_MID )
      IF ( ALLOCATED( T        ) ) DEALLOCATE( T        )
      IF ( ALLOCATED( TAUCLI   ) ) DEALLOCATE( TAUCLI   )
      IF ( ALLOCATED( TAUCLW   ) ) DEALLOCATE( TAUCLW   )
      IF ( ALLOCATED( TO31     ) ) DEALLOCATE( TO31     )
      IF ( ALLOCATED( TO32     ) ) DEALLOCATE( TO32     )
      IF ( ALLOCATED( TO3      ) ) DEALLOCATE( TO3      )
      IF ( ALLOCATED( TTO3     ) ) DEALLOCATE( TTO3     )
      IF ( ALLOCATED( TMPU1    ) ) DEALLOCATE( TMPU1    )
      IF ( ALLOCATED( TMPU2    ) ) DEALLOCATE( TMPU2    )
      IF ( ALLOCATED( TROPP    ) ) DEALLOCATE( TROPP    )
      IF ( ALLOCATED( TROPP1   ) ) DEALLOCATE( TROPP1   )
      IF ( ALLOCATED( TROPP2   ) ) DEALLOCATE( TROPP2   )
      IF ( ALLOCATED( TS       ) ) DEALLOCATE( TS       )
      IF ( ALLOCATED( TSKIN    ) ) DEALLOCATE( TSKIN    )
      IF ( ALLOCATED( U10M     ) ) DEALLOCATE( U10M     )
      IF ( ALLOCATED( UPDE     ) ) DEALLOCATE( UPDE     ) 
      IF ( ALLOCATED( UPDN     ) ) DEALLOCATE( UPDN     ) 
      IF ( ALLOCATED( USTAR    ) ) DEALLOCATE( USTAR    )
      IF ( ALLOCATED( UWND     ) ) DEALLOCATE( UWND     )
      IF ( ALLOCATED( UWND1    ) ) DEALLOCATE( UWND1    )
      IF ( ALLOCATED( UWND2    ) ) DEALLOCATE( UWND2    )
      IF ( ALLOCATED( V10M     ) ) DEALLOCATE( V10M     ) 
      IF ( ALLOCATED( VWND     ) ) DEALLOCATE( VWND     )
      IF ( ALLOCATED( VWND1    ) ) DEALLOCATE( VWND1    )
      IF ( ALLOCATED( VWND2    ) ) DEALLOCATE( VWND2    )
      IF ( ALLOCATED( Z0       ) ) DEALLOCATE( Z0       )
      IF ( ALLOCATED( ZMEU     ) ) DEALLOCATE( ZMEU     )
      IF ( ALLOCATED( ZMMD     ) ) DEALLOCATE( ZMMD     )
      IF ( ALLOCATED( ZMMU     ) ) DEALLOCATE( ZMMU     )
      IF ( ALLOCATED( EFLUX    ) ) DEALLOCATE( EFLUX    )
      IF ( ALLOCATED( AIRDEN_FULLGRID ) ) DEALLOCATE( AIRDEN_FULLGRID )
      IF ( ALLOCATED( T_FULLGRID      ) ) DEALLOCATE( T_FULLGRID      )

      END SUBROUTINE CLEANUP_DAO
!EOC
      END MODULE DAO_MOD
