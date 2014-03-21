!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: Pjc_Pfix_Mod
!
! !DESCRIPTION: Module Pjc\_Pfix\_Mod contains routines which implements the 
!  Philip Cameron-Smith pressure fixer for the new fvDAS transport 
!  scheme. (bdf, bmy, 5/8/03, 10/27/03)
!\\
!\\
! !INTERFACE:
!
      MODULE Pjc_Pfix_Mod
!
! !USES:
!
      IMPLICIT NONE
!
! !PUBLIC MEMBER FUNCTIONS:
!
      PUBLIC  :: Do_Pjc_Pfix
      PUBLIC  :: Cleanup_Pjc_Pfix
!
! !PRIVATE MEMBER FUNCTIONS:
!
      PRIVATE :: Calc_Pressure
      PRIVATE :: Calc_Advection_Factors
      PRIVATE :: Adjust_Press
      PRIVATE :: Init_Press_Fix
      PRIVATE :: Do_Press_Fix_LLNL
      PRIVATE :: Average_Press_Poles
      PRIVATE :: Convert_Winds
      PRIVATE :: Calc_Horiz_Mass_Flux
      PRIVATE :: Calc_Divergence
      PRIVATE :: Set_Press_Terms
      PRIVATE :: Do_Divergence_Pole_Sum
      PRIVATE :: Xpavg
      PRIVATE :: Init_Pjc_Pfix
!
! !AUTHOR:
!  Philip Cameron-Smith and John Tannahill, GMI project @ LLNL (2003)
!  Brendan Field and Bob Yantosca (5/8/03) 
!  Modified for new GMI TPCORE by Claire Carouge (ccarouge@seas.harvard.edu)
! 
! !REVISION HISTORY:
!  (1 ) Bug fix for Linux/PGI compiler in routines ADJUST_PRESS and 
!        INIT_PRESS_FIX. (bmy, 6/23/03)
!  (2 ) Now make P1, P2 true surface pressure in DO_PJC_PFIX (bmy, 10/27/03)
!  01 Mar 2012 - R. Yantosca - Now reference new grid_mod.F90
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !PRIVATE DATA MEMBERS:
! 
      PRIVATE :: AI,        BI,      DAP,      DBK,     COSE_FV
      PRIVATE :: COSP_FV,   CLAT_FV, DLAT_FV,  ELAT_FV, GEOFAC
      PRIVATE :: GW_FV,     MCOR,    REL_AREA, RGW_FV,  SINE_FV
      PRIVATE :: GEOFAC_PC, DLON_FV, LOC_PROC, PR_DIAG, IMP_NBORDER 
      PRIVATE :: I1_GL,     I2_GL,   JU1_GL,   JV1_GL,  J2_GL
      PRIVATE :: K1_GL,     K2_GL,   ILO_GL,   IHI_GL,  JULO_GL
      PRIVATE :: JVLO_GL,   JHI_GL,  I1,       I2,      JU1
      PRIVATE :: JV1,       J2,      K1,       K2,      ILO
      PRIVATE :: IHI,       JULO,    JVLO,     JHI,     ILAT
      PRIVATE :: ILONG,     IVERT,   J1P,      J2P       

!  ============================================================================
!  Module Variables:
!  ============================================================================
!  (1 ) AI          (REAL*8 )  : Vertical coord "A" for hybrid grid [hPa]
!  (2 ) BI          (REAL*8 )  : Vertical coord "B" for hybrid grid [unitless]
!  (3 ) CLAT_FV     (REAL*8 )  : Grid box center latitude [radians] 
!  (4 ) COSE_FV     (REAL*8 )  : COSINE of grid box edge latitudes [radians]
!  (5 ) COSP_FV     (REAL*8 )  : COSINE of grid box center latitudes [radians]
!  (6 ) DAP         (REAL*8 )  : Delta-A vertical coordinate [hPa]
!  (7 ) DBK         (REAL*8 )  : Delta-B vertical coordinate [unitless]
!  (8 ) DLAT_FV     (REAL*8 )  : Latitude extent of grid boxes [radians]
!  (9 ) ELAT_FV     (REAL*8 )  : Grid box edge latitudes [radians]
!  (10) GEOFAC      (REAL*8 )  : Geometric factor for N-S advection
!  (11) GW_FV       (REAL*8 )  : Diff of SINE btw grid box lat edges [unitless]
!  (12) MCOR        (REAL*8 )  : Grid box surface areas [m2]
!  (13) REL_AREA    (REAL*8 )  : Relative surface area of grid box [fraction]
!  (14) RGW_FV      (REAL*8 )  : Reciprocal of GW_FV [radians
!  (15) SINE_FV     (REAL*8 )  : SINE of lat at grid box edges [unitless]
!  (16) GEOFAC_PC   (REAL*8 )  : Geometric factor for N-S advection @ poles
!  (17) DLON_FV     (REAL*8 )  : Longitude extent of a grid box [radians]
!  (18) LOC_PROC    (REAL*8 )  : Local processor number
!  (19) PR_DIAG     (LOGICAL)  : Flag for printing diagnostic message
!  (20) IMP_NBORDER (INTEGER)  : Used for ghost zones for MPI ??? 
!  (21) I1_GL       (INTEGER)  : ind of 1st  global lon       (no ghost zones)
!  (22) I2_GL       (INTEGER)  : ind of last global lon       (no ghost zones)
!  (23) JU1_GL      (INTEGER)  : ind of 1st  global "u" lat   (no ghost zones)
!  (24) JV1_GL      (INTEGER)  : ind of 1st  global "v" lat   (no ghost zones)
!  (25) J2_GL       (INTEGER)  : ind of last global "u&v" lat (no ghost zones)
!  (26) K1_GL       (INTEGER)  : ind of 1st  global alt       (no ghost zones)
!  (27) K2_GL       (INTEGER)  : ind of last global alt       (no ghost zones)
!  (28) ILO_GL      (INTEGER)  : I1_GL  - IMP_NBORDER        (has ghost zones)
!  (29) IHI_GL      (INTEGER)  : I2_GL  + IMP_NBORDER        (has ghost zones)
!  (30) JULO_GL     (INTEGER)  : JU1_GL - IMP_NBORDER        (has ghost zones)
!  (31) JVLO_GL     (INTEGER)  : JV1_GL - IMP_NBORDER        (has ghost zones)
!  (32) JHI_GL      (INTEGER)  : J2_GL  + IMP_NBORDER        (has ghost zones)
!  (33) I1          (INTEGER)  : ind of first local lon       (no ghost zones)
!  (34) I2          (INTEGER)  : ind of last  local lon       (no ghost zones)
!  (35) JU1         (INTEGER)  : ind of first local "u" lat   (no ghost zones)
!  (36) JV1         (INTEGER)  : ind of first local "v" lat   (no ghost zones)
!  (37) J2          (INTEGER)  : ind of last  local "u&v" lat (no ghost zones)
!  (38) K1          (INTEGER)  : index of first local alt     (no ghost zones)
!  (39) K2          (INTEGER)  : index of last  local alt     (no ghost zones)
!  (40) ILO         (INTEGER)  : I1  - IMP_NBORDER           (has ghost zones)
!  (41) IHI         (INTEGER)  : I2  + IMP_NBORDER           (has ghost zones)
!  (42) JULO        (INTEGER)  : JU1 - IMP_NBORDER           (has ghost zones)
!  (43) JVLO        (INTEGER)  : JV1 - IMP_NBORDER           (has ghost zones)
!  (44) JHI         (INTEGER)  : J2  + IMP_NBORDER           (has ghost zones)
!
      ! Allocatable arrays
      REAL*8, ALLOCATABLE :: AI(:)
      REAL*8, ALLOCATABLE :: BI(:)
      REAL*8, ALLOCATABLE :: CLAT_FV(:)
      REAL*8, ALLOCATABLE :: COSE_FV(:)
      REAL*8, ALLOCATABLE :: COSP_FV(:)
      REAL*8, ALLOCATABLE :: DAP(:)
      REAL*8, ALLOCATABLE :: DBK(:)
      REAL*8, ALLOCATABLE :: DLAT_FV(:)
      REAL*8, ALLOCATABLE :: ELAT_FV(:)
      REAL*8, ALLOCATABLE :: GEOFAC(:)
      REAL*8, ALLOCATABLE :: GW_FV(:)
      REAL*8, ALLOCATABLE :: MCOR(:,:)
      REAL*8, ALLOCATABLE :: REL_AREA(:,:)
      REAL*8, ALLOCATABLE :: RGW_FV(:)
      REAL*8, ALLOCATABLE :: SINE_FV(:)

      ! Scalar variables
      LOGICAL             :: PR_DIAG
      INTEGER             :: LOC_PROC
      REAL*8              :: GEOFAC_PC
      REAL*8              :: DLON_FV

      ! Dimensions for GMI code (from "imp_dims")
      INTEGER             :: IMP_NBORDER
      INTEGER             :: I1_GL,  I2_GL,   JU1_GL,  JV1_GL            
      INTEGER             :: J2_GL,  K1_GL,   K2_GL,   ILO_GL    
      INTEGER             :: IHI_GL, JULO_GL, JVLO_GL, JHI_GL    
      INTEGER             :: I1,     I2,      JU1,     JV1       
      INTEGER             :: J2,     K1,      K2,      ILO       
      INTEGER             :: IHI,    JULO,    JVLO,    JHI       
      INTEGER             :: ILAT,   ILONG,   IVERT,   J1P       
      INTEGER             :: J2P       

      !=================================================================
      ! MODULE ROUTINES -- follow below the "CONTAINS" statement
      !=================================================================
      CONTAINS

!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Do_Pjc_Pfix
!
! !DESCRIPTION: Subroutine Do\_Pjc\_Pfix is the driver routine for the Philip 
!  Cameron-Smith pressure fixer for the fvDAS transport scheme. 
!  (bdf, bmy, 5/8/03, 3/5/07)
!\\
!\\
!  We assume that the winds are on the A-GRID, since this is the input that 
!  the fvDAS transport scheme takes. (bdf, bmy, 5/8/03)
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE Do_Pjc_Pfix( D_DYN, P1, P2, UWND, VWND, XMASS, YMASS )
!
! !USES:
!
      USE CMN_SIZE_MOD    ! Size parameters
      USE CMN_GCTM_MOD    ! Physical constants
!
! !INPUT PARAMETERS: 
!
      ! Dynamic timestep [s]
      REAL*8,  INTENT(IN)  :: D_DYN

!------------------------------------------------------------------------------
! Prior to 11/4/13:
! Eliminate array temporaries by accepting assumed-shape arguments into this
! routine. (mpayer, 11/4/13)
!      ! True PSurface at middle of dynamic timestep [hPa]
!      REAL*8,  INTENT(IN)  :: P1(IIPAR,JJPAR)
!
!      ! True PSurface at end    of dynamic timestep [hPa]
!      REAL*8,  INTENT(IN)  :: P2(IIPAR,JJPAR)
!------------------------------------------------------------------------------
      ! True PSurface at middle of dynamic timestep [hPa]
      REAL*8,  INTENT(IN)  :: P1(:,:)

      ! True PSurface at end    of dynamic timestep [hPa]
      REAL*8,  INTENT(IN)  :: P2(:,:)

      ! Zonal (E-W) wind [m/s]
      REAL*8,  INTENT(IN)  :: UWND(IIPAR,JJPAR,LLPAR)

      ! Meridional (N-S) wind [m/s]
      REAL*8,  INTENT(IN)  :: VWND(IIPAR,JJPAR,LLPAR)
!
! !OUTPUT PARAMETERS:
!
      ! E-W mass fluxes [mixing ratio]
      REAL*8,  INTENT(OUT) :: XMASS(IIPAR,JJPAR,LLPAR)

      ! N-S mass fluxes [mixing ratio]
      REAL*8,  INTENT(OUT) :: YMASS(IIPAR,JJPAR,LLPAR)
!
! !AUTHOR:
!  Brendan Field and Bob Yantosca (5/8/03) 
!
! !REMARKS:
!  (1 ) Now P1 and P2 are "true" surface pressures, and not PS-PTOP.  If using
!        this P-fixer w/ GEOS-3 winds, pass true surface pressure to this
!        routine. (bmy, 10/27/03)
!  (2 ) Now define P2_TMP array for passing to ADJUST_PRESS (yxw, bmy, 3/5/07)
!
! !REVISION HISTORY: 
!  02 Dec 2008 - R. Yantosca - Updated documentation and added ProTeX headers.
!                              Declare all REAL variables as REAL*8.
!  04 Nov 2013 - M. Sulprizio- Eliminate array temporaries by accepting assumed-
!                              shape dummy arguments
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
! 
      LOGICAL, SAVE        :: FIRST = .TRUE.
      INTEGER              :: I, J, K
      REAL*8               :: P2_TMP(IIPAR,JJPAR)
!
! !DEFINED PARAMETERS:
! 
      LOGICAL, PARAMETER   :: INTERP_WINDS     = .TRUE.  ! winds are interp'd 
      INTEGER, PARAMETER   :: MET_GRID_TYPE    = 0       ! A-GRID
      INTEGER, PARAMETER   :: ADVEC_CONSRV_OPT = 0       ! 2=floating pressure 
      INTEGER, PARAMETER   :: PMET2_OPT        = 1       ! leave at 1 
      INTEGER, PARAMETER   :: PRESS_FIX_OPT    = 1       ! Turn on P-Fixer

      !=================================================================
      ! DO_PJC_PFIX begins here!
      !=================================================================

      ! Initialize on first call
      IF ( FIRST ) THEN

         ! Initialize/allocate module variables 
         CALL INIT_PJC_PFIX

         ! Calculate advection surface-area factors
         CALL CALC_ADVECTION_FACTORS( MCOR, REL_AREA, GEOFAC, GEOFAC_PC)

         ! Reset first-time flag
         FIRST = .FALSE.
      ENDIF


      ! Copy P2 into P2_TMP (yxw, bmy, 3/5/07)
      P2_TMP = P2
 

      ! Call PJC pressure fixer w/ the proper arguments
      ! NOTE: P1 and P2 are now "true" surface pressure, not PS-PTOP!!!
      CALL ADJUST_PRESS( 'GEOS-CHEM',        INTERP_WINDS,  
     &                   .TRUE.,             MET_GRID_TYPE, 
     &                   ADVEC_CONSRV_OPT,   PMET2_OPT, 
     &                   PRESS_FIX_OPT,      D_DYN, 
     &                   GEOFAC_PC,          GEOFAC, 
     &                   COSE_FV,            COSP_FV, 
     &                   REL_AREA,           DAP, 
     &                   DBK,                P1, 
     &                   P2_TMP,             P2_TMP,            
     &                   UWND,               VWND, 
     &                   XMASS,              YMASS )

      
      ! Return to calling program
      END SUBROUTINE Do_Pjc_Pfix
!EOC 
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Calc_Pressure
!
! !DESCRIPTION: Subroutine Calc\_Pressure recalculates the new surface 
!  pressure from the adjusted air masses XMASS and YMASS.  This is useful 
!  for debugging purposes. (bdf, bmy, 5/8/03)
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE Calc_Pressure( XMASS, YMASS, RGW_FV, PS_NOW, PS_AFTER )
!
! !USES:
!
      USE CMN_SIZE_MOD  ! Size parameters
      USE CMN_MOD       ! NTRACE, LPRT, LWINDO
!
! !INPUT PARAMETERS: 
!
      ! E-W mass flux from pressure fixer
      REAL*8, INTENT(IN)  :: XMASS(IIPAR,JJPAR,LLPAR)

      ! N-S mass flux from pressure fixer
      REAL*8, INTENT(IN)  :: YMASS(IIPAR,JJPAR,LLPAR)

      ! Surface pressure - PTOP at current time
      REAL*8, INTENT(IN)  :: PS_NOW(IIPAR,JJPAR)

      ! 1 / ( SINE(J+1) - SINE(J) ) -- latitude factor
      REAL*8, INTENT(IN)  :: RGW_FV(JJPAR)
!
! !OUTPUT PARAMETERS:
!
      ! Surface pressure - PTOP adjusted by P-fixer
      REAL*8, INTENT(OUT) :: PS_AFTER(IIPAR,JJPAR)
!
! !AUTHOR:
!   Brendan Field and Bob Yantosca (5/8/03)
! 
! !REVISION HISTORY: 
!   02 Dec 2008 - R. Yantosca - Updated documentation and added ProTeX headers.
!                               Declare all REAL variables as REAL*8.
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
! 
      INTEGER             :: I, J, L
      REAL*8              :: DELP(IIPAR,JJPAR,LLPAR)
      REAL*8              :: DELP1(IIPAR,JJPAR,LLPAR)
      REAL*8              :: PE(IIPAR,LLPAR+1,JJPAR)

      !=================================================================
      ! CALC_PRESSURE begins here!
      !=================================================================
      DO L = 1, LLPAR
      DO J = 1, JJPAR
      DO I = 1, IIPAR
         DELP1(I,J,L) = DAP(L) + ( DBK(L) * PS_NOW(I,J) )
      ENDDO
      ENDDO
      ENDDO

      DO L = 1, LLPAR
         DO J = 2, JJPAR-1

            DO I =1, IIPAR-1
               DELP(I,J,L) = DELP1(I,J,L) + 
     &                       XMASS(I,J,L) - XMASS(I+1,J,L) + 
     &                     ( YMASS(I,J,L) - YMASS(I,J+1,L) ) * RGW_FV(J)
            ENDDO
            
            DELP(IIPAR,J,L) = 
     &          DELP1(IIPAR,J,L) +
     &          XMASS(IIPAR,J,L) - XMASS(1,J,L) +
     &        ( YMASS(IIPAR,J,L) - YMASS(IIPAR,J+1,L) ) * RGW_FV(J)
         ENDDO

         DO I = 1, IIPAR
            DELP(I,1,L) = DELP1(I,1,L) - YMASS(I,2,L) * RGW_FV(1)
         ENDDO

         ! Compute average
         CALL XPAVG( DELP(1,1,L), IIPAR )

         DO I = 1, IIPAR
            DELP(I,JJPAR,L) = DELP1(I,JJPAR,L) + 
     &                        YMASS(I,JJPAR,L) * RGW_FV(JJPAR)
         ENDDO

         ! Compute average
         CALL XPAVG( DELP(1,JJPAR,L), IIPAR )
      ENDDO

      !=================================================================
      ! Make the pressure
      !=================================================================
      DO J = 1, JJPAR
         DO I = 1, IIPAR
            PE(I,1,J) = PTOP
         ENDDO

         DO L = 1,LLPAR
            DO I = 1,IIPAR
               PE(I,L+1,J) = PE(I,L,J) + DELP(I,J,L)
            ENDDO
         ENDDO

         DO I = 1,IIPAR
            PS_AFTER(I,J) = PE(I,LLPAR+1,J)
         ENDDO
      ENDDO

      ! Return to calling program
      END SUBROUTINE Calc_Pressure
!EOC 
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Calc_Advection_Factors
!
! !DESCRIPTION: Subroutine Calc\_Advection\_Factors calculates the relative 
!   area of each grid box, and the geometrical factors used by this modified 
!   version of TPCORE.  These geomoetrical DO assume that the space is 
!   regularly gridded, but do not assume any link between the surface area 
!   and the linear dimensions.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE Calc_Advection_Factors
     &  (mcor, rel_area, geofac, geofac_pc)
!
! !USES:
!
      USE CMN_SIZE_MOD   ! Size parameters
      USE CMN_GCTM_MOD   ! Physical constants
!
! !INPUT PARAMETERS: 
!
      ! Area of grid box (m^2)
      REAL*8, INTENT(IN)  :: mcor(i1_gl :i2_gl, ju1_gl:j2_gl)
!
! !OUTPUT PARAMETERS:
!
      ! relative surface area of grid box (fraction)
      REAL*8, INTENT(OUT) :: rel_area(i1_gl :i2_gl, ju1_gl:j2_gl)

      ! Geometrical factor for meridional advection; geofac uses
      ! correct spherical geometry, and replaces acosp as the
      ! meridional geometrical factor in tpcore
      REAL*8, INTENT(OUT) :: geofac(ju1_gl:j2_gl)

      ! Special geometrical factor (geofac) for Polar cap
      REAL*8, INTENT(OUT) :: geofac_pc
!
! !AUTHOR:
!   Philip Cameron-Smith and John Tannahill, GMI project @ LLNL (2003)
! 
! !REMARKS:
!   Now reference PI from "CMN_GCTM" for consistency.  Also force
!   double-precision with the "D" exponent. (bmy, 5/6/03)
!
! !REVISION HISTORY: 
!   02 Dec 2008 - R. Yantosca - Updated documentation and added ProTeX headers.
!                               Declare all REAL variables as REAL*8.
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
! 
      INTEGER :: ij
 
      REAL*8  :: dp           ! spacing in latitude (rad)
      REAL*8  :: ri2_gl
      REAL*8  :: rj2m1
      REAL*8  :: total_area

      !----------------
      !Begin execution.
      !----------------

      ri2_gl = i2_gl

      !---------------------------------
      !Set the relative area (rel_area).
      !---------------------------------

      total_area = Sum (mcor(:,:))

      rel_area(:,:) = mcor(:,:) / total_area


      !---------------------------------------------------------
      !Calculate geometrical factor for meridional advection.
      !Note that it is assumed that all grid boxes in a latitude
      !band are the same.
      !---------------------------------------------------------

      rj2m1 = j2_gl - 1
      dp    = PI / rj2m1

      do ij = ju1_gl, j2_gl
        geofac(ij) = dp / (2.0d0 * rel_area(1,ij) * ri2_gl)
      end do

      geofac_pc =
     &  dp / (2.0d0 * Sum (rel_area(1,ju1_gl:ju1_gl+1)) * ri2_gl)

      ! Return to calling program
      END SUBROUTINE Calc_Advection_Factors
!EOC 
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Adjust_Press
!
! !DESCRIPTION: Subroutine Adjust\_Press initializes and calls the 
!  pressure fixer code.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE Adjust_Press
     &  (metdata_name_org, do_timinterp_winds, new_met_rec,
     &   met_grid_type, advec_consrv_opt, pmet2_opt, press_fix_opt,
     &   tdt, geofac_pc, geofac, cose, cosp, rel_area, dap, dbk,
     &   pctm1, pctm2, pmet2, uu, vv, xmass, ymass)
!
! !INPUT PARAMETERS: 
!
      ! First  part of metdata_name, e.g., "NCAR"
      CHARACTER(LEN=*) :: metdata_name_org

      ! Time interpolate wind fields?
      LOGICAL :: do_timinterp_winds

      ! New met record?
      LOGICAL :: new_met_rec

      ! Met grid type, A or C
      INTEGER :: met_grid_type

      ! Advection_conserve option
      INTEGER :: advec_consrv_opt

      ! pmet2 option
      INTEGER :: pmet2_opt

      ! pressure fixer option
      INTEGER :: press_fix_opt

      ! Model time step [s]
      REAL*8  :: tdt
      
      ! Special geometrical factor (geofac) for Polar cap
      REAL*8  :: geofac_pc

      ! Geometrical factor for meridional advection; geofac uses
      ! correct spherical geometry, and replaces acosp as the
      ! meridional geometrical factor in tpcore
      REAL*8  :: geofac  (ju1_gl:j2_gl)

      ! Cosines of grid box edges and centers
      REAL*8  :: cose    (ju1_gl:j2_gl)
      REAL*8  :: cosp    (ju1_gl:j2_gl)

      ! Pressure difference across layer from (ai * pt) term [hPa]
      REAL*8  :: dap     (k1:k2)

      ! Difference in bi across layer - the dSigma term
      REAL*8  :: dbk     (k1:k2)

      ! Relative surface area of grid box (fraction)
      REAL*8  :: rel_area( i1_gl:i2_gl,   ju1_gl:j2_gl)
            
      ! Metfield surface pressure at t1+tdt [hPa]
      REAL*8  :: pmet2(ilo_gl:ihi_gl, julo_gl:jhi_gl)

      ! CTM surface pressure at t1 [hPa]
      REAL*8  :: pctm1(ilo_gl:ihi_gl, julo_gl:jhi_gl)

      ! CTM surface pressure at t1+tdt [hPa]
      REAL*8  :: pctm2(ilo_gl:ihi_gl, julo_gl:jhi_gl)

      ! Wind velocity, x direction at t1+tdt/2 [m/s]
      REAL*8  :: uu(ilo_gl:ihi_gl, julo_gl:jhi_gl, k1_gl:k2_gl)

      ! Wind velocity, y direction at t1+tdt/2 [m/s]
      REAL*8  :: vv(ilo_gl:ihi_gl, julo_gl:jhi_gl, k1_gl:k2_gl)
!
! !INPUT/OUTPUT PARAMETERS: 
!
      ! Horizontal mass flux in E-W direction [hPa]
      REAL*8  :: xmass(ilo_gl:ihi_gl, julo_gl:jhi_gl, k1_gl:k2_gl)

      ! Horizontal mass flux in N-S direction [hPa]
      REAL*8  :: ymass(ilo_gl:ihi_gl, julo_gl:jhi_gl, k1_gl:k2_gl)
!
! !AUTHOR:
!   Philip Cameron-Smith and John Tannahill, GMI project @ LLNL (2003)
! 
! !REVISION HISTORY: 
!   02 Dec 2008 - R. Yantosca - Updated documentation and added ProTeX headers.
!                               Declare all REAL variables as REAL*8.
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
! 
      logical, save :: DO_ADJUST_PRESS_DIAG = .TRUE.


      !----------------------
      !Variable declarations.
      !----------------------

      logical, save :: first = .true.

      !--------------------------------------------------
      !dgpress   : global-pressure discrepancy
      !press_dev : RMS difference between pmet2 and pctm2
      !            (weighted by relative area)
      !--------------------------------------------------

      real*8  :: dgpress
      real*8  :: press_dev

      !-------------------------------------------------------------
      !dps : change of surface pressure from met field pressure [hPa]
      !-------------------------------------------------------------

      real*8  :: dps(i1_gl:i2_gl, ju1_gl:j2_gl)

      !--------------------------------------------
      !dps_ctm : CTM surface pressure tendency [hPa]
      !--------------------------------------------

      real*8 :: dps_ctm(i1_gl:i2_gl, ju1_gl:j2_gl)

      !---------------------------------------------------------------------
      !xmass_fixed : horizontal mass flux in E-W direction after fixing [hPa]
      !ymass_fixed : horizontal mass flux in N-S direction after fixing [hPa]
      !---------------------------------------------------------------------

      real*8  :: xmass_fixed(ilo_gl:ihi_gl, julo_gl:jhi_gl, k1:k2)
      real*8  :: ymass_fixed(ilo_gl:ihi_gl, julo_gl:jhi_gl, k1:k2)

      !-------------
      !Dummy indexes
      !-------------

      !integer :: ij, il

      !----------------
      !Begin execution.
      !----------------

      if (pr_diag) then
        Write (6, *) 'Adjust_Press called by ', loc_proc
      end if

      dps_ctm(:,:) = 0.0d0

      dgpress =  Sum ( (pmet2(i1_gl:i2_gl, ju1_gl:j2_gl) -
     &                  pctm1(i1_gl:i2_gl, ju1_gl:j2_gl)   )
     &             * rel_area(i1_gl:i2_gl, ju1_gl:j2_gl)     )

      if (pmet2_opt == 1) then
        pmet2(:,:) = pmet2(:,:) - dgpress
      end if

!### Debug
!###      if (DO_ADJUST_PRESS_DIAG) then
!###        Write (6, *) 'Global mean surface pressure change [hPa] = ', 
!###     &                dgpress
!###      end if

       !===================
        call Init_Press_Fix
       !===================
     &    (metdata_name_org, met_grid_type, tdt, geofac_pc, geofac,
     &     cose, cosp, dap, dbk, dps, dps_ctm, rel_area, pctm1, pmet2,
     &     uu, vv, xmass, ymass)

        if (press_fix_opt == 1) then

         !======================
          call Do_Press_Fix_Llnl
         !======================
     &      (geofac_pc, geofac, dbk, dps, dps_ctm, rel_area,
     &       xmass, ymass, xmass_fixed, ymass_fixed )
          
          xmass(:,:,:) = xmass_fixed(:,:,:)
          ymass(:,:,:) = ymass_fixed(:,:,:)

        end if

        if ((advec_consrv_opt == 0) .or.
     &      (advec_consrv_opt == 1)) then

          dps_ctm(i1_gl:i2_gl, ju1_gl:j2_gl) =
     &      pmet2(i1_gl:i2_gl, ju1_gl:j2_gl) - 
     &      pctm1(i1_gl:i2_gl, ju1_gl:j2_gl) 

        !-----------------------------------------------
        !else if (advec_consrv_opt == 2) then do nothing
        !-----------------------------------------------

        end if


      pctm2(i1_gl:i2_gl, ju1_gl:j2_gl) =
     &      pctm1(i1_gl:i2_gl, ju1_gl:j2_gl) + 
     &    dps_ctm(i1_gl:i2_gl, ju1_gl:j2_gl)


      if (DO_ADJUST_PRESS_DIAG) then

        !-------------------------------------------------------
        !Calculate the RMS pressure deviation (diagnostic only).
        !-------------------------------------------------------

        press_dev = 
     &    Sqrt (Sum (((pmet2(i1_gl:i2_gl,ju1_gl:j2_gl) -
     &                 pctm2(i1_gl:i2_gl,ju1_gl:j2_gl))**2 *
     &                rel_area(i1_gl:i2_gl,ju1_gl:j2_gl))))

!### Debug
!###        Write (6, *) 'RMS deviation between pmet2 & pctm2 [hPa] = ',
!###     &               press_dev

      end if

      ! Return to calling program
      END SUBROUTINE Adjust_Press
!EOC 
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Init_Press_Fix
!
! !DESCRIPTION: Subroutine Init\_Press\_Fix initializes the pressure fixer.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE Init_Press_Fix
     &  (metdata_name_org, met_grid_type, tdt, geofac_pc, geofac,
     &   cose, cosp, dap, dbk, dps, dps_ctm, rel_area, pctm1, pmet2,
     &   uu, vv, xmass, ymass)
!
! !INPUT PARAMETERS: 
!
      ! Model Time step [s]
      REAL*8 :: tdt

      ! First part of metdata_name, e.g., "NCAR"
      CHARACTER(LEN=*) :: metdata_name_org

      ! Met grid type, A or C
      INTEGER          :: met_grid_type

      ! Special geometrical factor (geofac) for Polar cap
      REAL*8           :: geofac_pc

      ! Cosine of grid box edges and centers
      REAL*8           :: cose(ju1_gl:j2_gl)
      REAL*8           :: cosp(ju1_gl:j2_gl)

      ! Geometrical factor for meridional advection; geofac uses
      ! correct spherical geometry, and replaces acosp as the
      ! meridional geometrical factor in tpcore
      REAL*8           :: geofac(ju1_gl:j2_gl)

      ! Pressure difference across layer from (ai * pt) term [hPa]
      REAL*8           :: dap(k1:k2)

      ! Difference in bi across layer - the dSigma term
      REAL*8           :: dbk(k1:k2)

      ! relative surface area of grid box (fraction)
      REAL*8           :: rel_area( i1_gl:i2_gl, ju1_gl:j2_gl)
            
      ! Metfield surface pressure at t1 [hPa]
      REAL*8           :: pmet2(ilo_gl:ihi_gl, julo_gl:jhi_gl)

      ! CTM surface pressure at t1 [hPa]
      REAL*8           :: pctm1(ilo_gl:ihi_gl, julo_gl:jhi_gl)

      ! CTM surface pressure at t1+tdt [hPa]
      REAL*8           :: pctm2(ilo_gl:ihi_gl, julo_gl:jhi_gl)

      ! Wind velocity, x direction at t1+tdt/2 [m/s]
      REAL*8           :: uu(ilo_gl:ihi_gl, julo_gl:jhi_gl, k1_gl:k2_gl)

      ! Wind velocity, y direction at t1+tdt/2 [m/s]
      REAL*8           :: vv(ilo_gl:ihi_gl, julo_gl:jhi_gl, k1_gl:k2_gl)
!
! !OUTPUT PARAMETERS:
!
      ! Horizontal mass flux in E-W direction [hPa]
      REAL*8  :: xmass(ilo_gl:ihi_gl, julo_gl:jhi_gl, k1_gl:k2_gl)

      ! Horizontal mass flux in N-S direction [hPa]
      REAL*8  :: ymass(ilo_gl:ihi_gl, julo_gl:jhi_gl, k1_gl:k2_gl)

      ! Change of surface pressure from met field pressure [hPa]
      REAL*8  :: dps(i1_gl:i2_gl, ju1_gl:j2_gl)

      ! CTM surface pressure tendency [hPa]
      REAL*8  :: dps_ctm(i1_gl:i2_gl, ju1_gl:j2_gl)
!
! !AUTHOR:
!   Philip Cameron-Smith and John Tannahill, GMI project @ LLNL (2003)
! 
! !REVISION HISTORY: 
!   02 Dec 2008 - R. Yantosca - Updated documentation and added ProTeX headers.
!                               Declare all REAL variables as REAL*8.
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
! 
      !--------------------------------------------------------------
      !dpi   : divergence at a grid point; used to calculate vertical
      !        motion [hPa]
      !--------------------------------------------------------------

      real*8  :: dpi(i1:i2, ju1:j2, k1:k2)

      !---------------------------------------------------------------------
      !crx   : Courant number in E-W direction
      !cry   : Courant number in N-S direction
      !delp1 : pressure thickness, the psudo-density in a hydrostatic system
      !        at t1 [hPa]
      !delpm : pressure thickness, the psudo-density in a hydrostatic system
      !        at t1+tdt/2 (approximate) [hPa]
      !pu    : pressure at edges in "u"  [hPa]
      !---------------------------------------------------------------------

      real*8  :: crx  (ilo:ihi, julo:jhi, k1:k2)
      real*8  :: cry  (ilo:ihi, julo:jhi, k1:k2)
      real*8  :: delp1(ilo:ihi, julo:jhi, k1:k2)
      real*8  :: delpm(ilo:ihi, julo:jhi, k1:k2)
      real*8  :: pu   (ilo:ihi, julo:jhi, k1:k2)

      !----------------
      !Begin execution.
      !----------------

      if (pr_diag) then
        Write (6,*) 'Init_Press_Fix called by ', loc_proc
      end if

      !========================
      call Average_Press_Poles
      !========================
     &  (rel_area, pctm1)

      !========================
      call Average_Press_Poles
      !========================
     &  (rel_area, pmet2)

      !-------------------------------------------------------------------
      !We need to calculate pressures at t1+tdt/2.  One ought to use pctm2
      !in the call to Set_Press_Terms, but since we don't know it yet, we
      !are forced to use pmet2.  This should be good enough because it is
      !only used to convert the winds to the mass fluxes, which is done
      !crudely anyway and the mass fluxes will still get fixed OK.
      !-------------------------------------------------------------------

      dps(i1:i2,ju1:j2) = pmet2(i1:i2,ju1:j2) - pctm1(i1:i2,ju1:j2)

      !====================
      call Set_Press_Terms
      !====================
     &  (dap, dbk, pctm1, pmet2, delp1, delpm, pu)


        !===================
        call Convert_Winds
        !===================
     &    (met_grid_type, tdt, cosp, crx, cry, uu, vv)
        

        !=========================
        call Calc_Horiz_Mass_Flux
        !=========================
     &    (cose, delpm, uu, vv, xmass, ymass, tdt, cosp)

      !====================
      call Calc_Divergence
      !====================
     &  (.false., geofac_pc, geofac, dpi, xmass, ymass)

        
      dps_ctm(i1:i2,ju1:j2) = Sum (dpi(i1:i2,ju1:j2,:), dim=3)

      ! Return to calling program
      END SUBROUTINE Init_Press_Fix
!EOC 
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Do_Press_Fix_Llnl
!
! !DESCRIPTION: Subroutine Do\_Press\_Fix\_Llnl fixes the mass fluxes to 
!  match the met field pressure tendency.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE Do_Press_Fix_Llnl
     &  (geofac_pc, geofac, dbk, dps, dps_ctm, rel_area,
     &   xmass, ymass, xmass_fixed, ymass_fixed)
!
! !INPUT PARAMETERS: 
!
      ! Special geometrical factor (geofac) for Polar cap
      REAL*8, INTENT(IN)   :: geofac_pc

      ! Geometrical factor for meridional advection; geofac uses
      ! correct spherical geometry, and replaces acosp as the
      !  meridional geometrical factor in tpcore
      REAL*8, INTENT(IN)   :: geofac(ju1_gl:j2_gl)

      ! Difference in bi across layer - the dSigma term
      REAL*8, INTENT(IN)   :: dbk(k1:k2)

      ! Change of surface pressure from met field pressure [hPa]
      REAL*8, INTENT(IN)   :: dps(i1:i2, ju1:j2)

      ! Relative surface area of grid box (fraction)
      REAL*8, INTENT(IN)   :: rel_area(i1:i2, ju1:j2)

      ! Horizontal mass fluxes in E-W and N-S directions [hPa]
      REAL*8, INTENT(IN)   :: xmass(ilo:ihi, julo:jhi, k1:k2)
      REAL*8, INTENT(IN)   :: ymass(ilo:ihi, julo:jhi, k1:k2)
!
! !OUTPUT PARAMETERS:
!
      ! Sum over vertical of dpi calculated from original mass fluxes [hPa]
      REAL*8,  INTENT(OUT) :: dps_ctm(i1:i2, ju1:j2)

      ! Horizontal mass flux in E-W and N-S directions after fixing [hPa]
      REAL*8,  INTENT(OUT) :: xmass_fixed(ilo:ihi, julo:jhi, k1:k2)
      REAL*8,  INTENT(OUT) :: ymass_fixed(ilo:ihi, julo:jhi, k1:k2)
!
! !AUTHOR:
!  Philip Cameron-Smith and John Tannahill, GMI project @ LLNL (2003)
!
! !REVISION HISTORY: 
!   02 Dec 2008 - R. Yantosca - Updated documentation and added ProTeX headers.
!                               Declare all REAL variables as REAL*8.
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
! 
      ! Scalars
      INTEGER :: il, ij, ik

      REAL*8  :: dgpress
      REAL*8  :: fxmean
      REAL*8  :: ri2

      ! Arrays
      REAL*8  :: fxintegral(i1:i2+1)
      REAL*8  :: mmfd(ju1:j2)
      REAL*8  :: mmf (ju1:j2)
      REAL*8  :: ddps(i1:i2, ju1:j2)

      !------------------------------------------------------------------------
      !dpi: divergence at a grid point; used to calculate vertical motion [hPa]
      !------------------------------------------------------------------------

      real*8  :: dpi(i1:i2, ju1:j2, k1:k2)

      real*8  :: xcolmass_fix(ilo:ihi, julo:jhi)


      !----------------
      !Begin execution.
      !----------------

      if (pr_diag) then
        Write (6,*) 'Do_Press_Fix_Llnl called by ', loc_proc
      end if


      ri2 = i2_gl

      mmfd(:) = 0.0d0

      xcolmass_fix(:,:)   = 0.0d0

      xmass_fixed (:,:,:) = xmass(:,:,:)
      ymass_fixed (:,:,:) = ymass(:,:,:)


      !------------------------------------------------------------
      !Calculate difference between GCM and LR predicted pressures.
      !------------------------------------------------------------

      ddps(:,:) = dps(:,:) - dps_ctm(:,:)


c     --------------------------------------
c     Calculate global-pressure discrepancy.
c     --------------------------------------

      dgpress =
     &  Sum (ddps(i1:i2,ju1:j2) * rel_area(i1:i2,ju1:j2))


      !----------------------------------------------------------
      !Calculate mean meridional flux divergence (df/dy).
      !Note that mmfd is actually the zonal mean pressure change,
      !which is related to df/dy by geometrical factors.
      !----------------------------------------------------------

      !------------------------
      !Handle non-Pole regions.
      !------------------------

      do ij = j1p, j2p
        mmfd(ij) = -(sum(ddps(:,ij)) / ri2 - dgpress) 
      end do

      !---------------------------------------------
      !Handle poles.
      !Note that polar boxes have all been averaged.
      !---------------------------------------------

       mmfd(ju1)   = -(ddps(1,ju1)   - dgpress)
       mmfd(ju1+1) = -(ddps(1,ju1+1) - dgpress)
       mmfd(j2-1)  = -(ddps(1,j2-1)  - dgpress)
       mmfd(j2)    = -(ddps(1,j2)    - dgpress)


      !---------------------------------------------
      !Calculate mean meridional fluxes (cos(e)*fy).
      !---------------------------------------------

       mmf(j1p) = mmfd(ju1) / geofac_pc

       do ij = j1p, j2p
          mmf(ij+1) = mmf(ij) + mmfd(ij) / geofac(ij)
       end do


      !------------------------------------------------------------
      !Fix latitude bands.
      !Note that we don't need to worry about geometry here because
      !all boxes in a latitude band are identical.
      !Note also that fxintegral(i2+1) should equal fxintegral(i1),
      !i.e., zero.
      !------------------------------------------------------------

      do ij = j1p, j2p

        fxintegral(:) = 0.0d0

        do il = i1, i2
          fxintegral(il+1) =
     &      fxintegral(il) -
     &      (ddps(il,ij) - dgpress) -
     &      mmfd(ij)
        end do

        fxmean = Sum (fxintegral(i1+1:i2+1)) / ri2

        do il = i1, i2
          xcolmass_fix(il,ij) = fxintegral(il) - fxmean
        end do

      end do

      !-------------------------------------
      !Distribute colmass_fix's in vertical.
      !-------------------------------------

      do ik = k1, k2      
        do ij = j1p, j2p
          do il = i1, i2

            xmass_fixed(il,ij,ik) = xmass(il,ij,ik) + 
     &                              xcolmass_fix(il,ij) * dbk(ik)

          end do
        end do
      end do

      do ik = k1, k2      
        do ij = j1p, j2p+1
          do il = i1, i2

            ymass_fixed(il,ij,ik) = ymass(il,ij,ik) +
     &                              mmf(ij) * dbk(ik)

          end do
        end do
      end do
 
      !====================
      call Calc_Divergence
      !====================
     &  (.false., geofac_pc, geofac, dpi, xmass_fixed, ymass_fixed)


      dps_ctm(i1:i2,ju1:j2) = Sum (dpi(i1:i2,ju1:j2,:), dim=3)

      ! Return to calling program
      END SUBROUTINE Do_Press_Fix_Llnl
!EOC 
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Average_Press_Poles
!
! !DESCRIPTION: Subroutine Average\_Press\_Poles averages pressure at the 
!  Poles when the Polar cap is enlarged.  It makes the last two latitudes 
!  equal.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE Average_Press_Poles
     &  (rel_area, press)
!
! !INPUT PARAMETERS:
!
      ! Relative surface area of grid box (fraction)
      REAL*8, INTENT(IN)    :: rel_area(i1:i2, ju1:j2)
!
! !OUTPUT PARAMETERS:
!
      ! Surface pressure [hPa]
      REAL*8, INTENT(INOUT) :: press   (ilo:ihi, julo:jhi)
!
! !AUTHOR:
!   Philip Cameron-Smith and John Tannahill, GMI project @ LLNL (2003) 
!
! 
! !REVISION HISTORY: 
!   02 Dec 2008 - R. Yantosca - Updated documentation and added ProTeX headers.
!                               Declare all REAL variables as REAL*8.
!EOP
!------------------------------------------------------------------------------
!BOC      
!
! !LOCAL VARIABLES:
!
      REAL*8  :: meanp


      !----------------
      !Begin execution.
      !----------------

      if (pr_diag) then
        Write (6,*) 'Average_Press_Poles called by ', loc_proc
      end if

      meanp =
     &  Sum (rel_area(i1:i2,ju1:ju1+1)  * 
     &       press   (i1:i2,ju1:ju1+1)) /
     &  Sum (rel_area(i1:i2,ju1:ju1+1))

      press(i1:i2,ju1:ju1+1) = meanp

      meanp =
     &  Sum (rel_area(i1:i2,j2-1:j2)  * 
     &       press   (i1:i2,j2-1:j2)) /
     &  Sum (rel_area(i1:i2,j2-1:j2))

      press(i1:i2,j2-1:j2) = meanp

      ! Return to calling program
      END SUBROUTINE Average_Press_Poles
!EOC 
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Convert_Winds
!
! !DESCRIPTION: Subroutine Convert\_Winds converts winds on A or C grid to 
!  Courant \# on C grid.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE Convert_Winds
     &  (igd, tdt, cosp, crx, cry, uu, vv)
!
! !USES:
!
      USE CMN_SIZE_MOD ! Size parameters
      USE CMN_GCTM_MOD ! Re, PI
!
! !INPUT PARAMETERS: 
!
      ! A or C grid
      INTEGER, INTENT(IN)  :: igd

      ! Model time step [s]
      REAL*8,  INTENT(IN)  :: tdt

      ! Cosine of grid box centers
      REAL*8,  INTENT(IN)  :: cosp(ju1_gl:j2_gl)

      ! Wind velocity in E-W (UU) and N-S (VV) directions at t1+tdt/2 [m/s]
      REAL*8,  INTENT(IN)  :: uu  (ilo:ihi, julo:jhi, k1:k2)
      REAL*8,  INTENT(IN)  :: vv  (ilo:ihi, julo:jhi, k1:k2)
!
! !OUTPUT PARAMETERS:
!
      ! Courant number in E-W (CRX) and N-S (CRY) directions
      REAL*8,  INTENT(OUT) :: crx (ilo:ihi, julo:jhi, k1:k2)
      REAL*8,  INTENT(OUT) :: cry (ilo:ihi, julo:jhi, k1:k2)

! !AUTHOR:
!  Philip Cameron-Smith and John Tannahill, GMI project @ LLNL (2003)
!
! !REMARKS:
!  Use GEOS-CHEM physical constants Re, PI to be consistent with other
!  usage everywhere (bmy, 5/5/03)
!
! !REVISION HISTORY: 
!   02 Dec 2008 - R. Yantosca - Updated documentation and added ProTeX headers.
!                               Declare all REAL variables as REAL*8.  
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
! 
      logical, save :: first = .true.

      integer :: il, ij

      !-------------------------------
      !dl : spacing in longitude (rad)
      !dp : spacing in latitude  (rad)
      !-------------------------------

      real*8  :: dl
      real*8  :: dp

      real*8  :: ri2
      real*8  :: rj2m1

      !------------------------
      !dtdy  : dt/dy      (s/m)
      !dtdy5 : 0.5 * dtdy (s/m)
      !------------------------

      real*8, save :: dtdy
      real*8, save :: dtdy5

      !------------------------
      !dtdx  : dt/dx      (s/m)
      !dtdx5 : 0.5 * dtdx (s/m)
      !------------------------

      real*8, allocatable, save :: dtdx (:)
      real*8, allocatable, save :: dtdx5(:)

      !----------------
      !Begin execution.
      !----------------

      if (pr_diag) then
        Write (6, *) 'Convert_Winds called by ', loc_proc
      end if


      !==========
      if (first) then
      !==========

        first = .false.


        Allocate (dtdx (ju1_gl:j2_gl))
        Allocate (dtdx5(ju1_gl:j2_gl))
        dtdx = 0.0d0; dtdx5 = 0.0d0

        ri2   = i2_gl
        rj2m1 = j2_gl - 1

        dl    = 2.0d0 * PI / ri2
        dp    = PI / rj2m1

        dtdy  = tdt / (Re * dp)
        dtdy5 = 0.5d0 * dtdy


        dtdx (ju1_gl) = 0.0d0
        dtdx5(ju1_gl) = 0.0d0

        do ij = ju1_gl + 1, j2_gl - 1

          dtdx (ij) = tdt / (dl * Re * cosp(ij))
          dtdx5(ij) = 0.5d0 * dtdx(ij)

        end do

        dtdx (j2_gl)  = 0.0d0
        dtdx5(j2_gl)  = 0.0d0

      end if


      !=============
      if (igd == 0) then  ! A grid.
      !=============

        do ij = ju1+1, j2-1
          do il = i1+1, i2
            crx(il,ij,:) =
     &        dtdx5(ij) *
     &        (uu(il,ij,:) + uu(il-1,ij,  :))
          end do
            crx(1,ij,:) =
     &        dtdx5(ij) *
     &        (uu(1,ij,:) + uu(i2,ij,  :))
        end do

        do ij = ju1+1, j2
          do il = i1, i2
            cry(il,ij,:) =
     &        dtdy5 *
     &        (vv(il,ij,:) + vv(il,  ij-1,:))
          end do
        end do


      !====
      else  ! C grid.
      !====

        do ij = ju1, j2
          do il = i1, i2

            crx(il,ij,:) =
     &        dtdx(ij) * uu(il-1,ij,  :)

            cry(il,ij,:) =
     &        dtdy     * vv(il,  ij-1,:)

          end do
        end do

      end if

      ! Return to calling program
      END SUBROUTINE Convert_Winds
!EOC 
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Calc_Horiz_Mass_Flux
!
! !DESCRIPTION: Subroutine Calc\_Horiz\_Mass\_Flux calculates the horizontal 
!  mass flux for non-GISS met data.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE Calc_Horiz_Mass_Flux
     &  (cose, delpm, uu, vv, xmass, ymass, tdt, cosp)
!
! !USES:
!
      USE CMN_SIZE_MOD ! Size parameters
      USE CMN_GCTM_MOD ! Re, Pi
!
! !INPUT PARAMETERS: 
!
      ! Timestep [s]
      REAL*8, INTENT(IN)   :: tdt

      ! Cosine of grid box edges
      REAL*8, INTENT(IN)   :: cose (ju1_gl:j2_gl)

      ! Cosine of grid box centers
      REAL*8, INTENT(IN)   :: cosp (ju1_gl:j2_gl)

      ! Pressure thickness, the pseudo-density in a 
      ! hdrostatic system  at t1+tdt/2 (approximate) [hPa]
      REAL*8, INTENT(IN)   :: delpm(ilo:ihi, julo:jhi, k1:k2)

      ! E-W (UU) and N-S (VV) winds [m/s]
      REAL*8, INTENT(IN)   :: uu  (ilo:ihi, julo:jhi, k1:k2)
      REAL*8, INTENT(IN)   :: vv  (ilo:ihi, julo:jhi, k1:k2)
!
! !OUTPUT PARAMETERS:
!
      ! Horizontal mass flux in E-W and N-S directions [hPa]
      REAL*8, INTENT(OUT)  :: xmass(ilo:ihi, julo:jhi, k1:k2)
      REAL*8, INTENT(OUT)  :: ymass(ilo:ihi, julo:jhi, k1:k2)

! !AUTHOR:
!   Original code from Shian-Jiann Lin, DAO 
!   John Tannahill, LLNL (jrt@llnl.gov)
! 
! !REMARKS:
!   Use GEOS-CHEM physical constants Re, PI to be consistent with other
!   usage everywhere (bmy, 5/5/03)

! !REVISION HISTORY: 
!   02 Dec 2008 - R. Yantosca - Updated documentation and added ProTeX headers.
!                               Declare all REAL variables as REAL*8. 
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
! 
      INTEGER :: ij
      INTEGER :: il
      INTEGER :: jst, jend
      REAL*8  :: dl
      REAL*8  :: dp

      REAL*8  :: ri2
      REAL*8  :: rj2m1
      REAL*8  :: factx
      REAL*8  :: facty


      !----------------
      !Begin execution.
      !----------------

      if (pr_diag) then
        Write (6,*) 'Calc_Horiz_Mass_Flux called by ', loc_proc
      end if

        ri2   = i2_gl
        rj2m1 = j2_gl - 1

        dl    = 2.0d0 * PI / ri2
        dp    = PI / rj2m1

        facty  = 0.5d0 * tdt / (Re * dp)

      !-----------------------------------
      !Calculate E-W horizontal mass flux.
      !-----------------------------------

      do ij = ju1, j2

       factx = 0.5d0 * tdt / (dl * Re * cosp(ij))

       do il = i1+1, i2 
        xmass(il,ij,:) = factx *
     &    (uu(il,ij,:) * delpm(il,ij,:)+
     &     uu(il-1,ij,:) * delpm(il-1,ij,:))
       end do

        xmass(i1,ij,:) = factx *
     &    (uu(i1,ij,:) * delpm(i1,ij,:)+
     &     uu(i2,ij,:) * delpm(i2,ij,:))

      end do

      !-----------------------------------
      !Calculate N-S horizontal mass flux.
      !-----------------------------------

      do ij = ju1+1, j2

         ymass(i1:i2,ij,:) = facty *
     &    cose(ij) * (vv(i1:i2,ij,:)*delpm(i1:i2,ij,:)+
     &    vv(i1:i2,ij-1,:)*delpm(i1:i2,ij-1,:))

      end do

      ymass(i1:i2,ju1,:) = facty *
     &    cose(ju1) * (vv(i1:i2,ju1,:)*delpm(i1:i2,ju1,:))


      ! Return to calling program
      END SUBROUTINE Calc_Horiz_Mass_Flux
!EOC 
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Calc_Divergence
!
! !DESCRIPTION: Subroutine Calc\_Divergence calculates the divergence.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE Calc_Divergence
     &  (do_reduction, geofac_pc, geofac, dpi, xmass, ymass)
!
! !INPUT PARAMETERS: 
!
      ! Set to F if called on Master; set to T if called by Slaves
      ! (NOTE: this doesn't seem to be used!)
      LOGICAL, INTENT(IN)    :: do_reduction

      ! Special geometrical factor (geofac) for Polar cap
      REAL*8,  INTENT(IN)    :: geofac_pc

      ! geometrical factor for meridional advection; geofac uses
      ! correct spherical geometry, and replaces acosp as the
      ! meridional geometrical factor in tpcore 
      REAL*8,  INTENT(IN)    :: geofac(ju1_gl:j2_gl)

      ! horizontal mass fluxes in E-W and N-S directions [hPa]
      REAL*8,  INTENT(IN)    :: xmass (ilo:ihi, julo:jhi, k1:k2)
      REAL*8,  INTENT(IN)    :: ymass (ilo:ihi, julo:jhi, k1:k2)
!
! !INPUT/OUTPUT PARAMETERS: 
!
      ! Divergence at a grid point; used to calculate vertical motion [hPa]
      REAL*8,  INTENT(INOUT) :: dpi   (i1:i2, ju1:j2, k1:k2)
!
! !AUTHOR:
!  Philip Cameron-Smith and John Tannahill, GMI project @ LLNL (2003)
! 
! !REVISION HISTORY: 
!   02 Dec 2008 - R. Yantosca - Updated documentation and added ProTeX headers.
!                               Declare all REAL variables as REAL*8. 
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
! 
      integer :: il, ij
      integer :: jst, jend

      !----------------
      !Begin execution.
      !----------------

      if (pr_diag) then
        Write (6,*) 'Calc_Divergence called by ', loc_proc
      end if

      !-------------------------
      !Calculate N-S divergence.
      !-------------------------

      do ij = j1p, j2p

        dpi(i1:i2,ij,:) =
     &    (ymass(i1:i2,ij,:) - ymass(i1:i2,ij+1,:)) *
     &    geofac(ij)

      end do

      !-------------------------
      !Calculate E-W divergence.
      !-------------------------

      do ij = j1p,j2p
        do il = i1, i2-1
          dpi(il,ij,:) =
     &      dpi(il,ij,:) +
     &      xmass(il,ij,:) - xmass(il+1,ij,:)
        end do
          dpi(i2,ij,:) =
     &      dpi(i2,ij,:) +
     &      xmass(i2,ij,:) - xmass(1,ij,:)
      end do

      !===========================
      call Do_Divergence_Pole_Sum
      !===========================
     &  (do_reduction, geofac_pc, dpi, ymass)
       

      ! Added this IF statemetn (ccarouge, 12/3/08)
      if (j1p /= ju1_gl+1) then

!       --------------------------------------------
!       Polar cap enlarged:  copy dpi to polar ring.
!       --------------------------------------------

        if (ju1 == ju1_gl) then

          dpi(:,ju1+1,:) = dpi(:,ju1,:)

        end if

        if (j2 == j2_gl) then

          dpi(:,j2-1,:)  = dpi(:,j2,:)

        end if

      end if

      ! Return to calling program
      END SUBROUTINE Calc_Divergence
!EOC 
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Set_Press_Terms
!
! !DESCRIPTION: Subroutine Set\_Press\_Terms sets the pressure terms.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE Set_Press_Terms
     &  (dap, dbk, pres1, pres2, delp1, delpm, pu)
!
! !INPUT PARAMETERS: 
!
      ! Pressure difference across layer from (ai * pt) term [hPa]
      REAL*8, INTENT(IN)  :: dap  (k1:k2)

      ! Difference in bi across layer - the dSigma term
      REAL*8, INTENT(IN)  :: dbk  (k1:k2)

      ! Surface pressure at t1 [hPa]
      REAL*8, INTENT(IN)  :: pres1(ilo:ihi, julo:jhi)

      ! Surface pressure at t1+tdt [hPa]
      REAL*8, INTENT(IN)  :: pres2(ilo:ihi, julo:jhi)
!
! !OUTPUT PARAMETERS:
!
      ! Pressure thickness, the psudo-density in a 
      ! hydrostatic system at t1 [hPa]
      REAL*8, INTENT(OUT) :: delp1(ilo:ihi, julo:jhi, k1:k2)

      ! Pressure thickness, the psudo-density in a 
      ! hydrostatic system at t1+tdt/2 (approximate)  [hPa]
      REAL*8, INTENT(OUT) :: delpm(ilo:ihi, julo:jhi, k1:k2)

      ! Pressure at edges in "u" [hPa]
      REAL*8, INTENT(OUT) :: pu   (ilo:ihi, julo:jhi, k1:k2)
!
! !AUTHOR:
!  Philip Cameron-Smith and John Tannahill, GMI project @ LLNL (2003)
! 
! !REVISION HISTORY: 
!   02 Dec 2008 - R. Yantosca - Updated documentation and added ProTeX headers.
!                               Declare all REAL variables as REAL*8.
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
! 
      integer :: il, ij, ik
      integer :: jst, jend

      !----------------
      !Begin execution.
      !----------------

      if (pr_diag) then
        Write (6,*) 'Set_Press_Terms called by ', loc_proc
      end if

      do ik = k1, k2

        delp1(:,:,ik) =
     &    dap(ik) + (dbk(ik) * pres1(:,:))

        delpm(:,:,ik) =
     &    dap(ik) + 
     &    (dbk(ik) * 0.5d0 * (pres1(:,:) + pres2(:,:)))

      end do

      do ij = ju1, j2
        do il = i1+1, i2
          pu(il,ij,:) =
     &      0.5d0 * (delpm(il,ij,:) + delpm(il-1,ij,:))
        end do

          pu(i1,ij,:) =
     &      0.5d0 * (delpm(i1,ij,:) + delpm(i2,ij,:))

      end do

      ! Return to calling program
      END SUBROUTINE Set_Press_Terms
!EOC 
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Do_Divergence_Pole_Sum
!
! !DESCRIPTION: Do\_Divergence\_Pole\_Sum sets the divergence at the Poles.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE Do_Divergence_Pole_Sum
     &  (do_reduction, geofac_pc, dpi, ymass)
!
! !INPUT PARAMETERS: 
!
      ! Set to T if called on Master; set to F if called by Slaves
      ! (NOTE: This does not seem to be used!)
      LOGICAL :: do_reduction

      ! Special geometrical factor (geofac) for Polar cap
      REAL*8  :: geofac_pc

      ! horizontal mass flux in N-S direction [hPa]
      REAL*8  :: ymass(ilo:ihi, julo:jhi, k1:k2)
!
! !OUTPUT PARAMETERS:
!
      ! Divergence at a grid point; used to calculate vertical motion [hPa]
      REAL*8  :: dpi  ( i1:i2,   ju1:j2,  k1:k2)
!
! !AUTHOR:
!  Philip Cameron-Smith and John Tannahill, GMI project @ LLNL (2003)
! 
! !REVISION HISTORY: 
!   02 Dec 2008 - R. Yantosca - Updated documentation and added ProTeX headers.
!                               Declare all REAL variables as REAL*8.
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
! 
!
      !----------------------
      !Variable declarations.
      !----------------------

      integer :: il, ik

      real*8  :: ri2

      real*8  :: mean_np(k1:k2)
      real*8  :: mean_sp(k1:k2)
      real*8  :: sumnp  (k1:k2)
      real*8  :: sumsp  (k1:k2)


      !----------------
      !Begin execution.
      !----------------

      ri2 = i2_gl

      !==================
      if (ju1 == ju1_gl) then
      !==================

        do ik = k1, k2

          sumsp(ik) = 0.0d0

          do il = i1, i2

            sumsp(ik) = sumsp(ik) + ymass(il,j1p,ik)

          end do

        end do

        do ik = k1, k2

          mean_sp(ik) = -sumsp(ik) / ri2 * geofac_pc

          do il = i1, i2

            dpi(il,ju1,ik) = mean_sp(ik)

          end do

        end do

      !======
      end if
      !======


      !================
      if (j2 == j2_gl) then
      !================

        do ik = k1, k2

          sumnp(ik) = 0.0d0

          do il = i1, i2

            sumnp(ik) = sumnp(ik) + ymass(il,j2p+1,ik)

          end do

        end do

        do ik = k1, k2

          mean_np(ik) = sumnp(ik) / ri2 * geofac_pc

          do il = i1, i2

            dpi(il,j2,ik) = mean_np(ik)

          end do

        end do

      !======
      end if
      !======

      ! Return to calling program
      END SUBROUTINE Do_Divergence_Pole_Sum
!EOC 
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Xpavg
!
! !description: Subroutine Xpavg replaces each element of a vector with 
!  the average of the entire array. (bmy, 5/7/03)
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE Xpavg( P, IM )
!
! !USES:
!
      ! References to F90 modules
      USE ERROR_MOD, ONLY : ERROR_STOP
!
! !INPUT PARAMETERS: 
!
      ! Dimension of P
      INTEGER, INTENT(IN)    :: IM
!
! !INPUT/OUTPUT PARAMETERS: 
!
      ! 1-D vector to be averaged
      REAL*8,  INTENT(INOUT) :: P(IM)

! !AUTHOR:
!   Philip Cameron-Smith and John Tannahill, GMI project @ LLNL (2003)
! 
! !REVISION HISTORY: 
!   02 Dec 2008 - R. Yantosca - Updated documentation and added ProTeX headers.
!                               Now make all REAL variables REAL*8.
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
! 
      REAL*8                 :: AVG

      !=================================================================
      ! XPAVG begins here!
      !=================================================================

      ! Error check IM
      IF ( IM == 0 ) THEN
         CALL ERROR_STOP( 'Div by zero!', 'XPAVG ("pjc_pfix_mod.f")' )
      ENDIF

      ! Take avg of entire P array
      AVG  = SUM( P ) / DBLE( IM ) 

      ! Store average value in all elements of P
      P(:) = AVG
 
      ! Return to calling program 
      END SUBROUTINE Xpavg
!EOC 
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Init_Pjc_Pfix      
!
! !DESCRIPTION: Subroutine Init\_Pjc\_Pfix allocates and initializes module 
!  arrays and variables.  GMI dimension variables will be used for 
!  compatibility with the Phil Cameron-Smith P-fixer. (bdf, bmy, 5/8/03)
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE Init_Pjc_Pfix
!
! !USES:
!
      ! References to F90 modules
      USE GRID_MOD,     ONLY : GET_AREA_M2, GET_YMID_R
      USE ERROR_MOD,    ONLY : ALLOC_ERR,   ERROR_STOP
      USE PRESSURE_MOD, ONLY : GET_AP,      GET_BP

      USE CMN_SIZE_MOD  ! Size parameters
      USE CMN_GCTM_MOD  ! Re, PI, etc...
! 
!
! !AUTHOR:
!   Brendan Field and Bob Yantosca (5/8/03)
! 
! !REVISION HISTORY: 
!  02 Dec 2008 - R. Yantosca - Updated documentation and added ProTeX headers.
!  01 Mar 2012 - R. Yantosca - Now use GET_AREA_CM2(I,J,L) from grid_mod.F90
!  01 Mar 2012 - R. Yantosca - Now use GET_YMID_R(I,J,L) from grid_mod.F90
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
! 
      ! Local variables
      INTEGER :: AS, I, J, L

      !=================================================================
      ! INIT_PJC_PFIX begins here!
      !
      ! Initialize dimensions for GMI pressure-fixer code 
      !=================================================================
      IMP_NBORDER = 0
      I1_GL       = 1 
      I2_GL       = IIPAR 
      JU1_GL      = 1 
      JV1_GL      = 1 
      J2_GL       = JJPAR 
      K1_GL       = 1 
      K2_GL       = LLPAR 
      ILO_GL      = I1_GL  - IMP_NBORDER 
      IHI_GL      = I2_GL  + IMP_NBORDER 
      JULO_GL     = JU1_GL - IMP_NBORDER 
      JVLO_GL     = JV1_GL - IMP_NBORDER 
      JHI_GL      = J2_GL  + IMP_NBORDER 
      I1          = I1_GL 
      I2          = I2_GL 
      JU1         = JU1_GL 
      JV1         = JV1_GL 
      J2          = J2_GL 
      K1          = K1_GL 
      K2          = K2_GL 
      ILO         = ILO_GL 
      IHI         = IHI_GL 
      JULO        = JULO_GL 
      JVLO        = JVLO_GL 
      JHI         = JHI_GL 
      ILAT        = J2_GL - JU1_GL + 1 
      ILONG       = I2_GL -  I1_GL + 1 
      IVERT       = K2_GL -  K1_GL + 1 
      J1P         = 3 
      J2P         = J2_GL - J1P + 1 

      ! Error check longitude
      IF ( ILONG /= IIPAR ) THEN
         CALL ERROR_STOP( 'Invalid longitude dimension ILONG!', 
     &                    'INIT_PJC_FIX ("pjc_pfix_mod.f")' )
      ENDIF
      
      ! Error check latitude
      IF ( ILAT /= JJPAR ) THEN
         CALL ERROR_STOP( 'Invalid latitude dimension ILAT!', 
     &                    'INIT_PJC_FIX ("pjc_pfix_mod.f")' )
      ENDIF

      ! Error check altitude
      IF ( IVERT /= LLPAR ) THEN
         CALL ERROR_STOP( 'Invalid altitude dimension IVERT!', 
     &                    'INIT_PJC_FIX ("pjc_pfix_mod.f")' )
      ENDIF
     
      !=================================================================      
      ! Allocate module arrays (use dimensions from GMI code)
      !=================================================================
      ALLOCATE( AI( K1_GL-1:K2_GL ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'AI' )

      ALLOCATE( BI( K1_GL-1:K2_GL ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'BI' )

      ALLOCATE( DAP( K1_GL:K2_GL ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'DAP' )

      ALLOCATE( DBK( K1_GL:K2_GL ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'DBK' )

      ALLOCATE( CLAT_FV( JU1_GL:J2_GL ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'CLAT_FV' )

      ALLOCATE( COSE_FV( JU1_GL:J2_GL+1 ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'COSE_FV' )

      ALLOCATE( COSP_FV( JU1_GL:J2_GL ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'COSP_FV' )

      ALLOCATE( DLAT_FV( JU1_GL:J2_GL ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'DLAT_FV' )

      ALLOCATE( ELAT_FV( JU1_GL:J2_GL+1 ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'ELAT_FV' )

      ALLOCATE( GEOFAC( JU1_GL:J2_GL ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'GEOFAC' )

      ALLOCATE( GW_FV( JU1_GL:J2_GL ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'GW_FV' )
      
      ALLOCATE( MCOR( I1_GL:I2_GL, JU1_GL:J2_GL ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'MCOR' )

      ALLOCATE( REL_AREA( I1_GL:I2_GL, JU1_GL:J2_GL ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'REL_AREA' )

      ALLOCATE( RGW_FV( JU1_GL:J2_GL ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'RGW_FV' )

      ALLOCATE( SINE_FV( JU1_GL:J2_GL+1 ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'SINE_FV' )

      !=================================================================
      ! Initialize arrays and variables
      !=================================================================

      ! Grid box surface areas [m2]
      DO J = JU1_GL, J2_GL
      DO I =  I1_GL, I2_GL
         MCOR(I,J) = GET_AREA_M2( I, J, 1 )
      ENDDO
      ENDDO

      ! Hybrid grid vertical coords: Ai [hPa] and Bi [unitless]
      DO L = K1_GL-1, K2_GL
         AI(L) = GET_AP( L+1 )
         BI(L) = GET_BP( L+1 )
      ENDDO

      ! Delta A [hPa] and Delta B [unitless]
      DO L = K1_GL, K2_GL
         !-------------------------------------------------------------
         ! NOTE:, this was the original code.  But since AI is already
         ! in hPa, we shouldn't need to multiply by PTOP again.  This
         ! should only matter for the fvDAS fields.  Also, DBK needs 
         ! to be positive (bmy, 5/8/03)
         !DAP(L) = ( AI(L) - AI(L-1) ) * PTOP   
         !DBK(L) = BI(L) - BI(L-1)
         !-------------------------------------------------------------
         DAP(L) = AI(L-1) - AI(L)  
         DBK(L) = BI(L-1) - BI(L)
      ENDDO
         
      ! Grid box center latitudes [radians]
      DO J = JU1_GL, J2_GL
         CLAT_FV(J) = GET_YMID_R( 1, J, 1 )
      ENDDO

      ! Longitude spacing
      DLON_FV    = 2.d0 * PI / DBLE( I2_GL )
       
      ! Latitude edge at south pole [radians]
      ELAT_FV(1) = -0.5d0 * PI 

      ! SIN and COS of lat edge at south pole [unitless]
      SINE_FV(1) = -1.d0
      COSE_FV(1) =  0.d0
         
      ! Latitude edges [radians] (w/ SIN & COS) at intermediate latitudes
      DO J = JU1_GL+1, J2_GL  !2, JJPAR
         ELAT_FV(J) = 0.5d0 * ( CLAT_FV(J-1) + CLAT_FV(J) )
         SINE_FV(J) = SIN( ELAT_FV(J) )
         COSE_FV(J) = COS( ELAT_FV(J) )
      ENDDO

      ! Latitude edge at North Pole [radians]
      ELAT_FV(J2_GL+1) = 0.5d0 * PI

      ! SIN of lat edge at North Pole
      SINE_FV(J2_GL+1) = 1.d0
      
      ! Latitude extent of South polar box [radians]
      DLAT_FV(1) = 2.d0 * ( ELAT_FV(2) - ELAT_FV(1) ) 

      ! Latitude extent of boxes at intermediate latitudes [radians]
      DO J = JU1_GL+1, J2_GL-1  ! 2, JJPAR-1
         DLAT_FV(J) = ELAT_FV(J+1) - ELAT_FV(J)
      ENDDO

      ! Latitude extent of North polar box [radians]
      DLAT_FV(J2_GL) = 2.d0 * ( ELAT_FV(J2_GL+1) - ELAT_FV(J2_GL) ) 

      ! Other stuff
      DO J = JU1_GL, J2_GL
         GW_FV(J)   = SINE_FV(J+1) - SINE_FV(J)
         COSP_FV(J) = GW_FV(J)     / DLAT_FV(J)
         RGW_FV(J)  = 1.d0         / GW_FV(J)
      ENDDO

      ! Return to calling program
      END SUBROUTINE Init_Pjc_Pfix
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Cleanup_Pjc_Pfix      
!
! !DESCRIPTION: Subroutine Cleanup\_Pjc\_Pfix deallocates all module arrays 
!  (bmy, 5/8/03)
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE Cleanup_Pjc_Pfix
! 
! !REVISION HISTORY: 
!   02 Dec 2008 - R. Yantosca - Updated documentation and added ProTeX headers.
!EOP
!------------------------------------------------------------------------------
!BOC
      !=================================================================
      ! CLEANUP_PJC_PFIX begins here!
      !=================================================================
      IF ( ALLOCATED( AI       ) ) DEALLOCATE( AI       )
      IF ( ALLOCATED( BI       ) ) DEALLOCATE( BI       )
      IF ( ALLOCATED( CLAT_FV  ) ) DEALLOCATE( CLAT_FV  )
      IF ( ALLOCATED( COSE_FV  ) ) DEALLOCATE( COSE_FV  )
      IF ( ALLOCATED( COSP_FV  ) ) DEALLOCATE( COSP_FV  )
      IF ( ALLOCATED( DAP      ) ) DEALLOCATE( DAP      )
      IF ( ALLOCATED( DBK      ) ) DEALLOCATE( DBK      )
      IF ( ALLOCATED( DLAT_FV  ) ) DEALLOCATE( DLAT_FV  )
      IF ( ALLOCATED( ELAT_FV  ) ) DEALLOCATE( ELAT_FV  )
      IF ( ALLOCATED( GEOFAC   ) ) DEALLOCATE( GEOFAC   )
      IF ( ALLOCATED( GW_FV    ) ) DEALLOCATE( GW_FV    )
      IF ( ALLOCATED( MCOR     ) ) DEALLOCATE( MCOR     )
      IF ( ALLOCATED( REL_AREA ) ) DEALLOCATE( REL_AREA )
      IF ( ALLOCATED( RGW_FV   ) ) DEALLOCATE( RGW_FV   )
      IF ( ALLOCATED( SINE_FV  ) ) DEALLOCATE( SINE_FV )

      ! Return to calling program
      END SUBROUTINE Cleanup_Pjc_Pfix

      END MODULE Pjc_Pfix_Mod
!EOC
