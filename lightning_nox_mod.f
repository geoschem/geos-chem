! $Id: lightning_nox_mod.f,v 1.5 2005/05/09 14:33:59 bmy Exp $
      MODULE LIGHTNING_NOX_MOD
!
!******************************************************************************
!  Module LIGHTNING_NOX_MOD contains variables and routines for emitting NOx
!  from lightning into the atmosphere.  Original code comes from the old 
!  GISS-II CTM's of Yuhang Wang, Gerry Gardner, & Larry Horowitz.  Cleaned 
!  up for inclusion into GEOS-CHEM. (bmy, 4/14/04, 3/16/05)
!
!  Module Variables:
!  ============================================================================
!  (1 ) NNLIGHT    (INTEGER)  : Number of vertical points in lightning CDF's
!  (2 ) NLTYPE     (INTEGER)  : Types of lightning to consider
!  (2 ) PROFILE    (REAL*8 )  : Array to hold lightning CDF's read from disk
!  (3 ) SLBASE     (REAL*8 )  : Array to hold NOx lightning emissions
!  (4 ) RFLASH     (REAL*8 )  : NOx molecules/flash/meter (based on 4 Tg N/y)
!  (5 ) CICG       (REAL*8 )  : Inter-Cloud/Cloud-Ground energy ratio
!  (6 ) TNCHARGE   (REAL*8 )  : Temperature of negative charge center [K]
!  (7 ) FLASHSCALE (REAL*8 )  : Scaling factor for lightning to 6 Tg N/yr
!
!  Module Routines:
!  ============================================================================
!  (1 ) LIGHTNING             : Driver routine for lightning emissions
!  (2 ) LIGHTDIST             : Partitions NOx vertically w/ Pickering CDF's
!  (3 ) FLASHES               : Computes flash rate and IC/CG ratio 
!  (4 ) EMLIGHTNING           : Saves lightning NOx into GEMISNOX array
!  (5 ) INIT_LIGHTNING_NOX    : Zeroes module arrays and reads CDF data
!  (6 ) CLEANUP_LIGHTNING_NOX : Deallocates all module arrays
!
!  GEOS-CHEM modules referenced by carbon_mod.f
!  ============================================================================
!  (1 ) dao_mod.f             : Module containing arrays for DAO met fields
!  (2 ) diag_mod.f            : Module containing GEOS-CHEM diagnostic arrays
!  (3 ) directory_mod.f       : Module containing GEOS-CHEM data & metfld dirs
!  (4 ) error_mod.f           : Module w/ I/O error and NaN check routines
!  (5 ) file_mod.f            : Contains file unit numbers and error checks
!  (6 ) grid_mod.f            : Module containing horizontal grid information
!  (7 ) pressure_mod.f        : Module containing routines to compute P(I,J,L)
!
!  References:
!  ============================================================================
!  (1  ) Price & Rind (1992), JGR, vol. 97, 9919-9933.
!  (2  ) Price & Rind (1994), M. Weather Rev, vol. 122, 1930-1939.
!
!  NOTES:
!  (1 ) Now references "directory_mod.f".  Now also bundle "flashes.f" into 
!        this module (bmy, 7/20/04)
!  (2 ) Update scaling for GEOS-4 in routine LIGHTNING (bmy, 3/16/05)
!******************************************************************************
!
      IMPLICIT NONE
 
      !=================================================================
      ! MODULE PRIVATE DECLARATIONS -- keep certain internal variables 
      ! and routines from being seen outside "lightning_nox_mod.f"
      !=================================================================

      ! Make everything PRIVATE ...
      PRIVATE

      ! ... except these routines
      PUBLIC :: EMLIGHTNING
      PUBLIC :: LIGHTNING
      PUBLIC :: CLEANUP_LIGHTNING_NOX

      !=================================================================
      ! MODULE VARIABLES
      !
      ! Preserve original notes from the old "CMN_NOX" header file here:
      ! ----------------------------------------------------------------
      ! We should have FLASHSCALE = 0.75 (since RFLASH must be 
      ! multiplied by 0.75 to change from 4 Tg/yr to 3 Tg/yr) but 
      ! for some reason this produces only 0.33 Tg NOx/yr.  The cheap 
      ! solution is to just set FLASHSCALE = 7.5 (thus including the
      ! factor of 10) and not ask too many more questions 
      ! (bmy, bey, mgs, 3/5/98)
      !
      ! Also, if you want 6 Tg N/yr from lightning, just double
      ! FLASHSCALE From 7.5 to 15.0. (rvm, bmy, 4/14/04)
      !=================================================================

      ! Scalars
      INTEGER              :: NNLIGHT
      INTEGER, PARAMETER   :: NLTYPE     = 3
      REAL*8,  PARAMETER   :: RFLASH     = 2.073d22
      REAL*8,  PARAMETER   :: CICG       = 1d0 / 3d0
      REAL*8,  PARAMETER   :: TNCHARGE   = 263.0d0
      REAL*8,  PARAMETER   :: FLASHSCALE = 15d0 

      ! Allocatable arrays
      REAL*8,  ALLOCATABLE :: PROFILE(:,:)
      REAL*8,  ALLOCATABLE :: SLBASE(:,:,:)

      !=================================================================
      ! MODULE ROUTINES -- follow below the "CONTAINS" statement 
      !=================================================================
      CONTAINS

!------------------------------------------------------------------------------

      SUBROUTINE LIGHTNING( T, CLDTOPS )
!
!******************************************************************************
!  Subroutine LIGHTNING uses Price & Rind's formulation for computing
!  NOx emission from lightning. (bmy, bey, mgs, bdf, 3/4/98, 3/15/05)
!
!  Arguments as input:
!  ============================================================================
!  (1 ) T          (REAL*8 ) : Temperature at grid box (I,J,L) in K
!  (2 ) CLDTOPS    (INTEGER) : Array of cloud top levels at (I,J)
!  (3 ) YLMID      (REAL*8 ) : Array of grid box latitude centers in degrees
!
!  Passed via CMN_NOX:
!  ============================================================================
!  (1 ) SLBASE     (REAL*8 ) : Array of NOx molecules released from lightning
!  (2 ) RFLASH     (REAL*8 ) : NOx molecules / meter released per flash
!  (3 ) CICG       (REAL*8 ) : IC / CG flash ratio
!  (4 ) FLASHSCALE (REAL*8 ) : Scaling factor 
!
!  Passed via CMN_GCTM:
!  ============================================================================
!  (1 ) Rdg0       : R /x g0 = 28.97 
!
!  Output Lightning NOX (molec/cm3/s) is stored in the GEMISNOX array.
!
!  References:
!  ============================================================================
!     Price & Rind (1992), JGR, vol 97, 9913.
!     Y. H. Wang   (1997), Paper I, submitted to JGR
!
!  NOTES:
!  (1 ) LIGHTNING is written in Fixed-Form Fortran 90.
!  (2 ) Pass PW = PS - PTOP to LIGHTNING.
!  (3 ) Avoid having to include the large common block file 'CMN' by 
!        passing arguments through the argument list.  We have to keep
!        CMN_NOX since that contains many lightning related variables.
!  (4 ) Multiply SLBASE by 0.522 for 2 x 2.5 model so that we get
!        the same Tg/yr of N as we do for the 4 x 5 model (bdf, 6/10/99)
!  (5 ) Remove PTOP from the arg list.  PTOP is now a parameter
!        in "CMN_SIZE".  (bmy, 2/10/00)
!  (6 ) For GEOS-3, multiply SLBASE by 0.4, to get approx 3 Tg/yr.  This
!        is necessary since the cloud tops are higher, and about 7 Tg/yr 
!        of NOx is produced. (mje, bmy, 7/3/01)
!  (7 ) T(IREF,JREF,L) is now replaced by T(I,J,L).  Also updated comments. 
!        (bmy, 9/25/01)
!  (8 ) Removed obsolete code from 9/01 (bmy, 10/24/01)
!  (9 ) Replaced all instances of IM with IIPAR and JM with JJPAR, in order
!        to prevent namespace confusion for the new TPCORE (bmy, 6/25/02)
!  (10) Now use GET_PEDGE and GET_PCENTER of "pressure_mod.f" to compute the
!        pressure at the bottom edge and center of grid box (I,J,L).  Remove
!        PW, G_SIG, G_SIGE, LCONVM, and LDEBUG from the arg list.  Eliminated
!        obsolete, commented-out code and debug code. (dsa, bdf, bmy, 8/21/02)
!  (11) Now replace YMID(JREF) with GET_YMID from "grid_mod.f".  Remove all
!        references to JREF. (bmy, 2/11/03)
!  (12) Scale lightning NOx for 1x1 China nested-grid in order to match the 
!        lightning NOx total from the 4x5 simulation. (yxw, bmy, 5/16/03)
!  (13) Scale GEOS-4 lightning NOx to 6 Tg N/yr (bmy, 3/8/04)
!  (14) Added parallel DO-loops for better optimization.  Now bundled into
!        "lightning_nox_mod.f" (bmy, 4/14/04)
!  (15) Added space in #ifdef block for GEOS-4 1 x 1.25 grid (bmy, 12/1/04)
!  (16) Added scale factors to make GEOS-4 lightning NOx come out to the same
!        total value as for 2001 GEOS-3 = 4.72 Tg N/yr (jal, bmy, 3/15/05)
!******************************************************************************
!      
      ! References to F90 modules
      USE DAO_MOD,      ONLY : BXHEIGHT
      USE GRID_MOD,     ONLY : GET_YMID
      USE PRESSURE_MOD, ONLY : GET_PEDGE, GET_PCENTER

#     include "CMN_SIZE"     ! Size parameters
#     include "CMN_GCTM"     ! Physical constants

      ! Arguments
      INTEGER, INTENT(IN)   :: CLDTOPS(IIPAR,JJPAR)
      REAL*8,  INTENT(IN)   :: T(IIPAR,JJPAR,LLPAR) 

      ! Local variables
      LOGICAL, SAVE         :: FIRST = .TRUE.
      INTEGER               :: I, J, L, LCHARGE, LMAX, LTOP
      REAL*8                :: P1, P2, P3, T1, T2, DLNP, DZ, ZUP
      REAL*8                :: HCHARGE, H0, XIGCRATIO, FLASHRATE
      REAL*8                :: X, RATE, TSUM, Z1, Z2, TOTAL, YMID
      REAL*8                :: VERTPROF(LLPAR)

#if   defined( GEOS_4 ) && defined( GRID4x5 ) 

      ! Scale GEOS-4 4x5 LNOx to GEOS-3 4x5 2001 total of 4.72 Tg N/yr
      REAL*8,  PARAMETER    :: G4_SCALE = ( 6.0d0    / 4.1749d0  ) * 
     &                                    ( 4.7238d0 / 2.5375d0  )

#elif defined( GEOS_4 ) && defined( GRID2x25 )

      ! Scale GEOS-4 2x25 LNOx to GEOS-3 2x25 2001 total of 4.72 Tg N/yr
      REAL*8,  PARAMETER    :: G4_SCALE = ( 6.0d0    / 12.355d0  ) *
     &                                    ( 4.7238d0 /  3.1523d0 )

#elif defined( GEOS_4 ) && defined( GRID1x125 )

      ! NOTE: Need to define this, set to 1 for now (bmy, 3/15/05)
      REAL*8,  PARAMETER    :: G4_SCALE = 1d0

#endif

      !=================================================================
      ! LIGHTNING begins here!
      !=================================================================

      ! First-time initialization
      IF ( FIRST ) THEN
         CALL INIT_LIGHTNING_NOX
         FIRST = .FALSE.
      ENDIF

      ! LMAX: the highest L-level to look for lightning (usually LLPAR-1)
      LMAX          = LLPAR - 1

      ! SLBASE : array containing molecules NOx / grid box / 6h. 
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L )
      DO L = 1, LLPAR
      DO J = 1, JJPAR
      DO I = 1, IIPAR
         SLBASE(I,J,L) = 0d0
      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      !=================================================================
      ! Loop over latitude/longtitude locations (I,J)
      !=================================================================
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I,         J,         L,    YMID, LCHARGE, P1       )    
!$OMP+PRIVATE( P2,        T1,        T2,   DLNP, DZ,      P3       )     
!$OMP+PRIVATE( ZUP,       HCHARGE,   LTOP, H0,   Z1,      Z2       )
!$OMP+PRIVATE( FLASHRATE, XIGCRATIO, RATE, X,    TOTAL,   VERTPROF )
!$OMP+SCHEDULE( DYNAMIC )
      DO J = 1, JJPAR

         ! Grid box latitude [degrees]
         YMID = GET_YMID( J )

         DO I = 1, IIPAR
   
            !===========================================================
            ! LCHARGE is the L-value where the negative charge layer is
            ! found.  According to Williams (1985), the negative charge 
            ! layer occurs where T = TNCHARGE = 263 K = - 10 C.  
            ! 
            ! If LCHARGE=1, then it is too cold to have water droplets, 
            ! so there will be no lightnings, so go to next (I,J) box.
            !===========================================================
            DO L = 1, LMAX
               IF ( T(I,J,L) <= TNCHARGE ) GOTO 10
            ENDDO

 10         CONTINUE
            IF ( L >= LMAX ) THEN
               LCHARGE = LMAX
            ELSE
               LCHARGE = L
            ENDIF

            IF ( LCHARGE == 1 ) GOTO 20

            !===========================================================
            ! P1 = pressure    (mb) at sigma center of level L = LCHARGE-1
            ! P2 = pressure    (mb) at sigma center of level L = LCHARGE
            ! T1 = temperature (K)  at sigma center of level L = LCHARGE-1
            ! T2 = temperature (K)  at sigma center of level L = LCHARGE
            !  
            ! DZ is the height from the sigma center of level 
            ! L = LCHARGE-1 to the negative charge layer.  Therefore, 
            ! DZ may be found in either the (LCHARGE)th sigma layer or 
            ! the (LCHARGE-1)th sigma layer.
            !  
            ! ZUP is the height from the sigma center of the 
            ! (LCHARGE-1)th layer  to the sigma edge of the 
            ! (LCHARGE-1)th layer
            !===========================================================
            P1   = GET_PCENTER(I,J,LCHARGE-1)
            P2   = GET_PCENTER(I,J,LCHARGE)

            T1   = T(I,J,LCHARGE-1)
            T2   = T(I,J,LCHARGE  )
         
            DLNP = LOG( P1 / P2 ) / ( T1 - T2 ) * ( T1 - TNCHARGE )
            DZ   = Rdg0 * ( (T1 + T2) / 2d0 ) * DLNP

            P3   = GET_PEDGE(I,J,LCHARGE)
            ZUP  = Rdg0 * T1 * LOG( P1 /P3 )

            !===========================================================
            ! HCHARGE is the height of the negative charge layer above 
            ! the bottom sigma edge of the (LCHARGE)th sigma level.  
            ! 
            ! If DZ >= ZUP then DZ is already in the (LCHARGE)th 
            ! sigma level.
            ! 
            ! If DZ <  ZUP then DZ is in the (LCHARGE-1)th sigma 
            ! level...therefore redefine LCHARGE = LCHARGE - 1 and 
            ! compute HCHARGE accordingly.  
            ! 
            ! Note that BXHEIGHT(I,J,LCHARGE) - ZUP is the distance 
            ! from the bottom sigma edge to the center of the newly 
            ! defined (LCHARGE)th layer.
            !===========================================================
            IF ( DZ >= ZUP ) THEN
               HCHARGE = DZ - ZUP
            ELSE
               LCHARGE = LCHARGE-1
               HCHARGE = ( BXHEIGHT(I,J,LCHARGE) - ZUP ) + DZ
            ENDIF
 
            !===========================================================
            ! LTOP is the L-layer where the convective cloud top is 
            ! found.  The cloud top is located at the highest sigma 
            ! level for which the cloud mass flux is nonzero.  Since DAO 
            ! cloud mass flux is defined at the top of each sigma level, 
            ! the convective cloud top is located at the top edge of 
            ! layer LTOP.
            ! 
            ! For lightning to exist, the cloud must straddle the 
            ! negative charge layer (in other words, at the very 
            ! minimum, the cloud bottom must occur in the LCHARGEth 
            ! layer).  If LTOP < LCHARGE go to the next (I,J) location.
            ! 
            ! H0 is the convective cloud top height in meters.  H0 is 
            ! thus given by the distance from the ground to the top 
            ! edge of layer LTOP.
            !===========================================================
            LTOP = CLDTOPS(I,J)
            IF ( LTOP > LMAX    ) LTOP = LMAX
            IF ( LTOP < LCHARGE ) GOTO 20
            
            H0 = 0.0
            DO L = 1, LTOP
               H0 = H0 + BXHEIGHT(I,J,L)
            ENDDO

            !===========================================================
            ! Z1 = Cloud-Ground (CG) pathway (from ground to HCHARGE)    
            !      in meters
            !
            ! Z2 = Inter-Cloud  (IC) pathway (from HCHARGE to cloud top) 
            !      in meters
            !===========================================================
            Z1 = 0.0
            DO L = 1, LCHARGE-1
               Z1 = Z1 + BXHEIGHT(I,J,L)
            ENDDO
            Z1 = Z1 + HCHARGE

            Z2 = 0.0
            DO L = LCHARGE, LTOP
               Z2 = Z2 + BXHEIGHT(I,J,L)
            ENDDO
            Z2 = Z2 - HCHARGE         

            !===========================================================
            ! Call FLASHES to compute the number of lighting strikes 
            ! per minute at Lat/Long position (I,J), based on cloud top 
            ! height, according to Price & Rind (1992).  
            ! 
            ! FLASHRATE = (flashes/min)
            ! XICGRATIO = Intercloud (IC) flashes / Cloud-Ground (CG) flashes
            ! RATE      = (flashes/6h)  = (flashes/min) * (360 min/6h)
            ! X         = the ratio of CG flashes / total # of flashes
            ! TOTAL     = the total number of molecules released
            !           = (molec/flash/m) * ( (1-X)*0.3333*Z2 + X*Z1 ) *
            !             FLASHRATE*60 * FLASHSCALE
            ! 
            ! Multiply the total by FLASHSCALE to obtain 3 Tg N/yr, 
            ! to match the best-observed distribution. 
            !
            ! IMPORTANT NOTE!!!!
            ! ==================
            ! We should have FLASHSCALE = 0.75 (since RFLASH must be 
            ! multiplied by 0.75 to change from 4 Tg/yr to 3 Tg/yr) but 
            ! for some reason this produces only 0.33 Tg NOx/yr.  The 
            ! cheap solution is to just set FLASHSCALE = 7.5 (thus 
            ! including the factor of 10) and not ask too many more 
            ! questions (bmy, bey, mgs, 3/5/98)
            ! 
            ! LIGHTDIST computes the lightning distribution from the 
            ! ground to the convective cloud top using cumululative 
            ! distribution functions for ocean flashes, tropical land 
            ! flashes, and non-tropical land flashes, as specified by 
            ! Ken Pickering (1997). 
            !===========================================================
            CALL FLASHES( I, J, H0, FLASHRATE, XIGCRATIO )

            RATE  = FLASHRATE * 360.0
            X     = 1d0 / ( 1d0 + XIGCRATIO )
            TOTAL = RFLASH * ( ( (1 - X) * CICG * Z2 ) + ( X * Z1 ) ) * 
     &              (FLASHRATE * 60) * FLASHSCALE

            IF ( TOTAL > 0 ) THEN
               CALL LIGHTDIST( I, J, LTOP, H0, YMID, TOTAL, VERTPROF )
            ENDIF

            ! Add back into SLBASE
            DO L = 1, LLPAR
               SLBASE(I,J,L) = SLBASE(I,J,L) + VERTPROF(L)
            ENDDO

            ! Go to next (I,J) position
 20         CONTINUE
         ENDDO  ! I
      ENDDO     ! J
!$OMP END PARALLEL DO

#if   defined( GEOS_4 ) 

      !================================================================= 
      ! Scaling: For GEOS-4 met fields only
      !=================================================================

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L )
      DO L = 1, LLPAR
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! Scale GEOS-4 LNOx to match the same total that we get with
         ! GEOS-3 meteorology for 2001 = 4.72 Tg N/yr (jal, bmy, 3/16/05)
         SLBASE(I,J,L) = SLBASE(I,J,L) * G4_SCALE

      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

#else 

      !=================================================================
      ! Scaling: For GEOS-1, GEOS-STRAT, and GEOS-3
      !=================================================================

#if   defined ( GRID2x25 )

      !----------------------------
      ! 2 x 2.5 (all model types)
      !----------------------------

      ! In the 2x2.5 model the emissions from lightning are too high.
      ! Emissions need to be scaled by the factor .522 to match the 4x5 
      ! emissions. 
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L )
      DO L = 1, LLPAR
      DO J = 1, JJPAR
      DO I = 1, IIPAR
         SLBASE(I,J,L) = SLBASE(I,J,L) * .522d0
      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

#endif

#if   defined( GEOS_3 ) 

      !----------------------------
      ! GEOS-3 (all resolutions)
      !----------------------------

      ! Also, for GEOS-3, the convective cloud tops are a lot higher 
      ! than in GEOS-1 or GEOS-STRAT, so more NOx will be produced.  
      ! Mat Evans found that approx 7.2 Tg NOx was produced in 1 year 
      ! by GEOS-3.  We scale that to 3 Tg NOx by multiplying by 
      ! 3/7.2 ~= 0.4 (mje, bmy, 7/3/01)
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L )
      DO L = 1, LLPAR
      DO J = 1, JJPAR
      DO I = 1, IIPAR
         SLBASE(I,J,L) = SLBASE(I,J,L) * 0.4d0
      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

#endif

#if   defined( GRID1x1 ) && defined( NESTED_CH )

      !----------------------------
      ! For China nested grid
      !----------------------------

      ! Further scaling is necessary for 1x1 nested grids in order to
      ! match the amount of lightning from 4x5 (yxw, bmy, 5/16/03)
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L )
      DO L = 1, LLPAR
      DO J = 1, JJPAR
      DO I = 1, IIPAR
         SLBASE(I,J,L) = SLBASE(I,J,L) * 0.158d0
      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

#endif

#if   defined( GRID1x1 ) && defined( NESTED_NA )

      !----------------------------
      ! For N. America nested grid
      !----------------------------

      ! NOTE! Need to define this!

#endif

#endif

      ! Return to calling program
      END SUBROUTINE LIGHTNING

!------------------------------------------------------------------------------

      SUBROUTINE LIGHTDIST( I, J, LTOP, H0, XLAT, TOTAL, VERTPROF )
!
!******************************************************************************
!  Subroutine LIGHTDIST reads in the CDF used to partition the 
!  column lightning NOx into the GEOS-CHEM vertical layers. 
!  (yhw, 1997; mje, bmy, 9/18/02, 7/20/04)
!
!  Arguments as Input:
!  ============================================================================
!  (1-2) I, J     (INTEGER) : (Lon,Lat) indices of the surface grid box
!  (3  ) LTOP     (INTEGER) : Level where convective cloud top is found
!  (4  ) H0       (REAL*8 ) : Convective cloud top height [m]
!  (5  ) XLAT     (REAL*8 ) : Latitude of surface grid box (I,J) [degrees]
!  (6  ) TOTAL    (REAL*8 ) : Total # of NOx molec. released from lightning
!  (7  ) VERTPROF (REAL*8 ) :
!
!  References:
!  ============================================================================
!  (1 ) Pickering et al., JGR 103, 31,203 - 31,316, 1998.
!
!  NOTES:
!  (1 ) Use functions IS_LAND and IS_WATER to determine if the given grid
!        box is over land or water.  These functions work for all DAO met
!        field data sets. (bmy, 4/2/02)
!  (2 ) Renamed M2 to LTOP and THEIGHT to H0 for consistency w/ variable names
!        w/in "lightning.f".  Now read the "light_dist.dat.geos3" file for 
!        GEOS-3 directly from the DATA_DIR/lightning_NOx_200203/ subdirectory.
!        Now read the "light_dist.dat" file for GEOS-1, GEOS-STRAT directly 
!        from the DATA_DIR/lightning_NOx_200203/ subdirectory.  Added 
!        descriptive comment header.  Now trap I/O errors across all 
!        platforms with subroutine "ioerror.f".  Updated comments, cosmetic 
!        changes.  Redimension FRAC(NNLIGHT) to FRAC(LLPAR). (bmy, 4/2/02)
!  (3 ) Deleted obsolete code from April 2002.  Now reference IU_FILE and
!        IOERROR from "file_mod.f".  Now use IU_FILE instead of IUNIT as the
!        file unit number. (bmy, 6/27/02)
!  (4 ) Now reference BXHEIGHT from "dao_mod.f" (bmy, 9/18/02)
!  (5 ) Bug fix: add GEOS_4 to the #if block (bmy, 3/4/04)
!  (6 ) Now bundled into "lightning_mod.f".  CDF's are now read w/in
!        routine INIT_LIGHTNING to allow parallelization (bmy, 4/14/04)
!  (7 ) Now references DATA_DIR from "directory_mod.f" (bmy, 7/20/04)
!******************************************************************************
!
      ! References to F90 modules
      USE DAO_MOD,       ONLY : IS_LAND, IS_WATER, BXHEIGHT
      USE DIRECTORY_MOD, ONLY : DATA_DIR
      USE FILE_MOD,      ONLY : IU_FILE, IOERROR

#     include "CMN_SIZE"   ! Size parameters

      ! Arguments
      INTEGER, INTENT(IN)  :: I,  J,    LTOP
      REAL*8,  INTENT(IN)  :: H0, XLAT, TOTAL
      REAL*8,  INTENT(OUT) :: VERTPROF(LLPAR)

      ! Local variables
      INTEGER              :: M, MTYPE, L, III, IOS, IUNIT, JJJ
      REAL*8               :: ZHEIGHT
      REAL*8               :: FRAC(LLPAR)
      CHARACTER(LEN=255)   :: FILENAME

      !=================================================================
      ! LIGHTDIST begins here!
      !=================================================================

      ! Initialize 
      DO L = 1, LLPAR
         VERTPROF(L) = 0d0
      ENDDO

      !=================================================================
      ! Use functions IS_LAND and IS_WATER from "dao_mod.f" to test 
      ! whether surface location (I,J) is over land or over water.
      !
      ! Depending on the combination of land/water and latitude, assign
      ! a flag describing the type of lightning:
      !
      !   MTYPE = 1: ocean lightning
      !   MTYPE = 2: tropical continental lightning (from 30S to 30N)
      !   MTYPE = 3: midlatitude continental lightning (all other lats)
      !=================================================================
      IF ( IS_LAND(I,J) ) THEN
         IF ( ABS( XLAT ) <= 30 ) THEN
            MTYPE = 2
         ELSE
            MTYPE = 3
         ENDIF
      ELSE IF ( IS_WATER(I,J) ) THEN 
         MTYPE = 1
      ENDIF

#if   defined( GEOS_1 ) || defined( GEOS_STRAT )

      !=================================================================
      ! For GEOS-1 and GEOS-STRAT only:
      !
      ! Use the PDF for this type of lightning to partition the total
      ! column lightning into the GEOS-3 vertical layers
      !=================================================================

      ZHEIGHT = 0.

      ! Compute the height [km] at the top of each vertical level.
      ! Look up the cumulative fraction of NOx for each vertical level
      DO L = 1, LTOP-1
         ZHEIGHT = ZHEIGHT + BXHEIGHT(I,J,L)
         FRAC(L) = PROFILE( NINT( ( ZHEIGHT / H0 ) * 100. ), MTYPE )
      ENDDO

      ! If there is any lightning NOx yet to be partitioned out,
      ! then place that in the level where the cloud top occurs.
      FRAC(LTOP) = 1.
         
      ! Convert from cumulative fraction to fraction for each level
      DO L = LTOP, 2, -1
         FRAC(L) = FRAC(L) - FRAC(L-1)
      ENDDO

      ! Partition lightning NOx by layer into VERTPROF
      DO L = 1, LTOP
         VERTPROF(L) = ( FRAC(L) * TOTAL )
      ENDDO

#else

      !=================================================================
      ! For all other met field types (GEOS-3, GEOS-4, GEOS-5, ... )
      !
      ! Use the CDF for this type of lightning to partition the total
      ! column lightning into the GEOS-3 vertical layers.
      !=================================================================
      ZHEIGHT = 0.0

      ! Compute the height [km] at the top of each vertical level.
      ! Look up the cumulative fraction of NOx for each vertical level
      DO L = 1, LTOP
         ZHEIGHT = ZHEIGHT + BXHEIGHT(I,J,L)
         FRAC(L) = PROFILE( NINT( ( ZHEIGHT/H0 )*3200. ), MTYPE ) *0.01
      ENDDO

      ! Convert from cumulative fraction to fraction for each level
      DO L = LTOP, 2, - 1
         FRAC(L) = FRAC(L) - FRAC(L-1)
      ENDDO 
      
      ! Partition lightning NOx by layer into VERTPROF
      DO L = 1, LTOP
         VERTPROF(L) = ( FRAC(L) * TOTAL )
      ENDDO

#endif

      ! Return to calling program
      END SUBROUTINE LIGHTDIST

!------------------------------------------------------------------------------

      SUBROUTINE FLASHES( I, J, HEIGHT, FLASHRATE, IC_CG_RATIO )
!
!******************************************************************************
!  Subroutine FLASHES determines the rate of lightning flashes per minute,
!  based on the height of convective cloud tops, and the inter-cloud
!  to cloud-ground strike ratio.  FLASHES has been optimized for GEOS-CHEM
!  (bmy, 10/9/97, 7/20/04)
!
!  Arguments as Input:
!  ============================================================================
!  (1-2) I, J        : GEOS-CHEM longitude & latitude indices
!  (3  ) HEIGHT      : Height of convective cloud tops [m]
!
!  Arguments as Output:
!  ============================================================================
!  (4  ) FLASHRATE   : Lightning flash rate [flashes/minute]
!  (5  ) IC_CG_RATIO : Intercloud (IC) flashes / Cloud-Ground (CG) flashes
!
!  References:
!  ============================================================================
!  (1  ) Price & Rind (1992), JGR, vol. 97, 9919-9933.
!  (2  ) Price & Rind (1994), M. Weather Rev, vol. 122, 1930-1939.
!
!  NOTES:
!  (1  ) FLASHES is written in Fixed-Form Fortran 90.  Also use F90 
!         declaration syntax. (bmy, 6/26/00)
!  (2  ) Eliminate obsolete code from 6/26/00 (bmy, 8/31/00)
!  (3  ) Bundled into "lightning_nox_mod.f".  Cleaned up stuff. (bmy, 7/20/04)
!******************************************************************************
!    
      ! References to F90 modules
      USE DAO_MOD, ONLY : IS_LAND, IS_WATER

      ! Arguments
      INTEGER, INTENT(IN)  :: I, J
      REAL*8,  INTENT(IN)  :: HEIGHT
      REAL*8,  INTENT(OUT) :: FLASHRATE, IC_CG_RATIO
     
      !=================================================================
      ! FLASHES begins here!
      !
      ! Price & Rind (1992) give the following parameterizations 
      ! for lightning flash rates as a function of convective cloud
      ! top height [km]:
      !
      !    FLAND  = 3.44e-5 * ( CLDTOP HEIGHT [km] ^ 4.9  )
      !    FOCEAN = 6.4e-4  * ( CLDTOP HEIGHT [km] ^ 1.73 )
      !
      ! Lightning will therefore occur much more often on land.  It goes
      ! as approx. the 5th power of height, as opposed to approx. the
      ! 2nd power of height over oceans.
      !
      ! Price & Rind (1992) also compute the ratio between Inter-Cloud 
      ! and Cloud-Ground flashes by the parameterization:
      !
      !    IC/CG  = 2.7d0 * SQRT( lightning flashes per minute )
      !=================================================================
      IF ( IS_LAND(I,J) ) THEN

         ! Flashes/minute over land (cf. Price & Rind 1992)
         FLASHRATE = 3.44d-5 * ( ( HEIGHT * 1d-3 )**4.9d0  )  

      ELSE IF ( IS_WATER(I,J) ) THEN 

         ! Flashes/minute over water (cf. Price & Rind 1992)
         FLASHRATE = 6.4d-4 *  ( ( HEIGHT * 1d-3 )**1.73d0 )   

      ENDIF

      ! Compute the IC/CG ratio (cf. Price & Rind 1992)
      IC_CG_RATIO = 2.7d0 * SQRT( FLASHRATE ) 

      ! Return to calling program
      END SUBROUTINE FLASHES

!------------------------------------------------------------------------------

      SUBROUTINE EMLIGHTNING( I, J )
!
!******************************************************************************
!  Subroutine EMLIGHTNING converts lightning emissions to [molec/cm3/s]
!  and stores them in the GEMISNOX array, which gets passed to SMVGEAR.
!  (bmy, 10/9/97, 4/14/04)
!
!  NOTES:
!  (1 ) Remove IOFF, JOFF from the argument list.  Also remove references
!        to header files "CMN_O3" and "comtrid.h" (bmy, 3/16/00)
!  (2 ) Now use allocatable array for ND32 diagnostic (bmy, 3/16/00)  
!  (3 ) Now reference BXHEIGHT from "dao_mod.f".  Updated comments, cosmetic
!        changes.  Replace LCONVM with the parameter LLCONVM. (bmy, 9/18/02)
!  (4 ) Removed obsolete reference to "CMN".  Now bundled into 
!        "lightning_mod.f" (bmy, 4/14/04)
!******************************************************************************
!
      ! References to F90 modules
      USE DAO_MOD,  ONLY : BXHEIGHT
      USE DIAG_MOD, ONLY : AD32_li

      IMPLICIT NONE

#     include "CMN_SIZE"  ! Size parameters
#     include "CMN_DIAG"  ! ND32
#     include "CMN_NOX"   ! GEMISNOX

      ! Arguments
      INTEGER, INTENT(IN) :: I, J

      ! Local variables
      INTEGER             :: L
      REAL*8              :: TMP

      ! External functions
      REAL*8, EXTERNAL    :: BOXVL

      !=================================================================
      ! EMLIGHTNING begins here!
      !=================================================================
      DO L = 1, LLCONVM 

          ! SLBASE(I,J,L) has units [molec NOx/6h/box], convert units:
          ! [molec/6h/box] * [6h/21600s] * [box/BOXVL cm3] = [molec/cm3/s]
          TMP             = SLBASE(I,J,L) / ( 21600.d0 * BOXVL(I,J,L) )
          GEMISNOX(I,J,L) = GEMISNOX(I,J,L) + TMP

          ! ND32 Diagnostic: Lightning NOx [molec NOx/cm2/s]
          IF ( ND32 > 0 ) THEN
             AD32_li(I,J,L) = AD32_li(I,J,L) + 
     &                        ( TMP * BXHEIGHT(I,J,L) * 1d2 )
          ENDIF
      ENDDO

      ! Return to calling program
      END SUBROUTINE EMLIGHTNING

!------------------------------------------------------------------------------

      SUBROUTINE INIT_LIGHTNING_NOX
!
!******************************************************************************
!  Subroutine INIT_LIGHTNING_NOX allocates all module arrays.  It also reads 
!  the lightning CDF data from disk before the first lightning timestep. 
!  (bmy, 4/14/04)
!
!  NOTES:
!  (1 ) Now reference DATA_DIR from "directory_mod.f"
!******************************************************************************
!
      ! References to F90 modules
      USE DIRECTORY_MOD, ONLY : DATA_DIR
      USE ERROR_MOD,     ONLY : ALLOC_ERR
      USE FILE_MOD

#     include "CMN_SIZE"  ! Size parameters
  
      ! Local variables
      INTEGER            :: AS, III, IOS, JJJ
      CHARACTER(LEN=255) :: FILENAME

      !=================================================================
      ! INIT_LIGHTNING_NOX begins here!
      !=================================================================

      ! NNLIGHT is the number of points for the lightning PDF's
#if   defined( GEOS_1 ) || defined( GEOS_STRAT )
      NNLIGHT = 100
#else
      NNLIGHT = 3200
#endif

      ! Allocate PROFILE
      ALLOCATE( PROFILE( NNLIGHT, NLTYPE ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'PROFILE' )
      PROFILE = 0d0

      ! Allocate SLBASE
      ALLOCATE( SLBASE( IIPAR, JJPAR, LLPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'SLBASE' )
      SLBASE = 0d0

#if   defined( GEOS_1 ) || defined( GEOS_STRAT )

      !=================================================================
      ! Read lightning CDF data for GEOS-1, GEOS-STRAT
      !
      ! Here we read in the original Pickering CDF's for NNLIGHT=3
      ! different types of lightning: tropical marine, tropical 
      ! continental, and midlatitude continental.  The vertical 
      ! resolution of the CDF's in the file read in below is 0.16 km. 
      !=================================================================

      ! Define filename for GEOS-1, GEOS-STRAT PDF file
      FILENAME = TRIM( DATA_DIR ) // 
     &           'lightning_NOx_200203/light_dist.dat'        

      ! Echo info
      WRITE( 6, 100 ) TRIM( FILENAME )
 100  FORMAT( '     - INIT_LIGHTNING: Reading ', a )

      ! Open file containing lightning CDF data
      OPEN( IU_FILE, FILE=TRIM( FILENAME ), STATUS='OLD', IOSTAT=IOS )
      IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_FILE, 'lightdist:3' )

      ! For GEOS-1, GEOS-STRAT: Read 1 header line
      READ( IU_FILE, *, IOSTAT=IOS )
      IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_FILE, 'lightdist:4' )

         ! For GEOS-1, GEOS-STRAT: Read data
      DO III = 1, NNLIGHT
         READ( IU_FILE,*,IOSTAT=IOS ) ( PROFILE(III,JJJ), JJJ=1,NLTYPE )
         IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_FILE, 'lightdist:5' )
      ENDDO

      ! Close file
      CLOSE( IU_FILE )

#else
         
      !=================================================================
      ! Read lightning CDF data for GEOS-3, GEOS-4 and higher
      ! 
      ! NOTE: Since the GEOS-3 vertical grid has a much finer
      ! resolution near the surface than do the GEOS-1 or GEOS-STRAT
      ! grids, we had to interpolate the Pickering CDF's to a much
      ! finer mesh, with 3200 points instead of 100.  This was done
      ! by Mat Evans (mje@io.harvard.edu).  The vertical resolution
      ! of the CDF's in the file read in below is 0.05 km. 
      !=================================================================

      ! Define filename for GEOS-3 CDF file
      FILENAME = TRIM( DATA_DIR ) // 
     &           'lightning_NOx_200203/light_dist.dat.geos3'        

      ! Echo info
      WRITE( 6, 100 ) TRIM( FILENAME )
 100  FORMAT( '     - INIT_LIGHTNING: Reading ', a )
      
      ! Open file containing lightning PDF data
      OPEN( IU_FILE, FILE=TRIM( FILENAM E), STATUS='OLD', IOSTAT=IOS )
      IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_FILE, 'lightdist:1' )
         
      ! For GEOS-3 only: read 12 header lines
      DO III = 1, 12
         READ( IU_FILE, '(a)', IOSTAT=IOS ) 
         IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_FILE, 'lightdist:2' )
      ENDDO
         
      ! For GEOS-3: read NNLIGHT types of lightning profiles
      DO III = 1, NNLIGHT
         READ( IU_FILE,*,IOSTAT=IOS) (PROFILE(III,JJJ),JJJ=1,NLTYPE)
      ENDDO
         
      ! Close file
      CLOSE( IU_FILE )

#endif

      ! Return to calling program
      END SUBROUTINE INIT_LIGHTNING_NOX

!------------------------------------------------------------------------------

      SUBROUTINE CLEANUP_LIGHTNING_NOX
!
!******************************************************************************
!  Subroutine CLEANUP_LIGHTNING_NOX deallocates all module arrays. 
!  (bmy, 4/14/04)
!
!  NOTES:
!******************************************************************************
!
      !=================================================================
      ! CLEANUP_LIGHTNING_NOX begins here!
      !=================================================================
      IF ( ALLOCATED( PROFILE ) ) DEALLOCATE( PROFILE )
      IF ( ALLOCATED( SLBASE  ) ) DEALLOCATE( SLBASE  )

      ! Return to calling program
      END SUBROUTINE CLEANUP_LIGHTNING_NOX

!------------------------------------------------------------------------------

      END MODULE LIGHTNING_NOX_MOD
