! $Id: lightning.f,v 1.1 2003/06/30 20:26:05 bmy Exp $
      SUBROUTINE LIGHTNING( T, CLDTOPS )
!
!******************************************************************************
!  Subroutine LIGHTNING uses Price & Rind's formulation for computing
!  NOx emission from lightning. (bmy, bey, mgs, bdf, 3/4/98, 5/16/03)
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
!******************************************************************************
!      
      ! References to F90 modules
      USE DAO_MOD,      ONLY : BXHEIGHT
      USE GRID_MOD,     ONLY : GET_YMID
      USE PRESSURE_MOD, ONLY : GET_PEDGE, GET_PCENTER

      IMPLICIT NONE

#     include "CMN_SIZE"  ! Size parameters
#     include "CMN_GCTM"  ! Physical constants
#     include "CMN_NOX"   ! SLBASE

      ! Arguments
      INTEGER, INTENT(IN) :: CLDTOPS(IIPAR,JJPAR)
      REAL*8,  INTENT(IN) :: T(IIPAR,JJPAR,LLPAR) 

      ! Local variables
      INTEGER             :: I, J, L, LCHARGE, LMAX, LTOP
      REAL*8              :: P1, P2, P3, T1, T2, DLNP, DZ, ZUP
      REAL*8              :: HCHARGE, H0, XIGCRATIO, FLASHRATE
      REAL*8              :: X, RATE, TSUM, Z1, Z2, TOTAL, YMID

      !=================================================================
      ! LIGHTNING begins here!
      !
      ! Initialize some variables.
      !
      ! SLBASE : array containing molecules NOx / grid box / 6h.  
      !
      ! LMAX   : the highest L-level to look for lightning.  Usually,
      !          LMAX = LCONVM, where LCONVM is the highest level 
      !          used for convection = LLPAR - 1.
      !=================================================================
      LMAX          = LLPAR - 1
      SLBASE(:,:,:) = 0d0

      !=================================================================
      ! Loop over latitude/longtitude locations (I,J)
      !=================================================================
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
               CALL LIGHTDIST( I, J, LTOP, H0, YMID, TOTAL )
            ENDIF

            ! Go to next (I,J) position
 20         CONTINUE
         ENDDO  ! I
      ENDDO     ! J

      !=================================================================
      ! In the 2x2.5 model the emissions from lightning are too high.
      ! Emissions need to be scaled by the factor .522 to match the 4x5 
      ! emissions.
      !
      ! Also, for GEOS-3, the convective cloud tops are a lot higher 
      ! than in GEOS-1 or GEOS-STRAT, so more NOx will be produced.  
      ! Mat Evans found that approx 7.2 Tg NOx was produced in 1 year 
      ! by GEOS-3.  We scale that to 3 Tg NOx by multiplying by 
      ! 3/7.2 ~= 0.4 (mje, bmy, 7/3/01)
      !=================================================================
#if   defined ( GRID2x25 )
      SLBASE(:,:,:) = SLBASE(:,:,:) * .522d0
#endif

#if   defined( GEOS_3 ) 
      SLBASE(:,:,:) = SLBASE(:,:,:) * 0.4d0
#endif

      !=================================================================
      ! Further scaling is necessary for 1x1 nested grids in order to
      ! match the amount of lightning from 4x5 (yxw, bmy, 5/16/03)
      !=================================================================
#if   defined( GRID1x1 ) 

      ! For China nested-grid
      SLBASE(:,:,:) = SLBASE(:,:,:) * 0.158d0

      ! North America, Europe nested grid scale factors need to be defined
#endif

      ! Return to calling program
      END SUBROUTINE LIGHTNING





