!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: upbdflx_mod
!
! !DESCRIPTION: Module UPBDFLX\_MOD contains subroutines which impose 
!  stratospheric boundary conditions on O3 and NOy.
!\\
!\\
! !INTERFACE: 
!
      MODULE UPBDFLX_MOD
!
! !USES:
!
#if defined( DEVEL )
!      USE OMP_LIB
#endif
      IMPLICIT NONE
#     include "define.h"
      PRIVATE

!
! !PUBLIC MEMBER FUNCTIONS:
! 
      PUBLIC  :: DO_UPBDFLX
      PUBLIC  :: UPBDFLX_O3
      PUBLIC  :: UPBDFLX_NOY
      PUBLIC  :: UPBDFLX_HD
      PUBLIC  :: INIT_UPBDFLX
!
! !REVISION HISTORY:
!  28 Jun 2001 - Q. Li, B. Field, M. Evans, R. Yantosca - Initial version
!  (1 ) Routine "upbdflx_noy" now correctly reprocessed P(NOy) files from
!        /data/ctm/GEOS_4x5/pnoy_200106 or /data/ctm/GEOS_2x2.5/pnoy_200106.
!        (mje, bmy, 6/28/01)
!  (2 ) Updated comments (bmy, 9/4/01)
!  (3 ) Fixes for reading binary punch files of global size (bmy, 9/27/01)
!  (4 ) Removed obsolete commented out code from 9/01 (bmy, 10/24/01)
!  (5 ) Removed obsolete commented out code from 7/01 (bmy, 11/26/01)
!  (6 ) Updated comments (bmy, 5/28/02)
!  (7 ) Replaced all instances of IM with IIPAR and JM with JJPAR, in ordr
!        to prevent namespace confusion for the new TPCORE (bmy, 6/25/02)
!  (8 ) Now references "pressure_mod.f" (dsa, bdf, bmy, 8/21/02)
!  (9 ) Now references BXHEIGHT from "dao_mod.f".  Also deleted obsolete
!        code from 8/02.  Now references IDTNOx, IDTOX, from "tracerid_mod.f" 
!        instead of from "comtrid.h". (bmy, 11/6/02)
!  (10) Added driver routine DO_UPBDFLX.  Also added lat limits for 1x1 in
!        UPBDFLX_O3. (bmy, 3/14/03)
!  (11) Now references AD from "dao_mod.f" in UPBDFLX_NOY (bnd, bmy, 4/14/03)
!  (12) Added printout of O3 in Tg/yr in UPBDFLX_O3 (mje, bmy, 8/15/03)
!  (13) Change O3 flux for GEOS-3 to 500 Tg/yr in UPBDFLX_O3 (bmy, 9/15/03)
!  (14) Now references "tagged_ox_mod.f" (bmy, 8/19/03)
!  (15) Now activated parallel DO loops (bmy, 4/15/04)
!  (16) Now made IORD, JORD, KORD module variables.  Now added routine
!        SET_UPBDFLX.  Now added routine SET_TRANSPORT (bmy, 7/20/04)
!  (17) Bug fix for COMPAQ compiler.  Now supports 1x125 grid. (bmy, 12/1/04)
!  (18) Now supports GEOS-5 and GCAP grids (swu, bmy, 5/25/05)
!  (19) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (20) Now references "tropopause_mod.f" (bmy, 11/1/05)
!  (21) Remove support for GEOS-1 and GEOS-STRAT met fields (bmy, 8/4/06)
!  (22) Added UPBDFLX_HD from the strat-trop flux of HD (lyj, phs, 9/18/07)
!  (23) Cap 1-XRATIO in UPBDFLX_NOY to prevent underflow (phs, 6/30/08)
!  (24) Modifications for GEOS-5 nested grid (yxw, dan, bmy, 11/6/08)
!  (25) Remove support for COMPAQ compiler (bmy, 7/8/09)
!  (26) Added support for LINOZ (dbj, jliu, bmy, 10/16/09)
!  13 Aug 2010 - R. Yantosca - Add modifications for MERRA (treat like GEOS-5)
!  02 Dec 2010 - R. Yantosca - Added ProTeX headers
!  08 Feb 2012 - R. Yantosca - Add modifications for GEOS-5.7.2 met data
!  28 Feb 2012 - R. Yantosca - Removed support for GEOS-3
!  07 Mar 2012 - R. Yantosca - Change DO loop boundaries for GI model
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !PRIVATE TYPES:
!
      ! TPCORE transport option flags
      INTEGER :: IORD, JORD, KORD

      CONTAINS
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: do_upbdflx
!
! !DESCRIPTION: Subroutine DO\_UPBDFLX is the driver routine for the 
!  stratospheric (upper-boundary) routines for Ox and NOy.
!\\
!\\
! !INTERFACE:
!
#if defined( DEVEL )
      SUBROUTINE DO_UPBDFLX( LOCAL_MET )
#else
      SUBROUTINE DO_UPBDFLX
#endif
!
! !USES:
!
      USE ERROR_MOD,   ONLY : DEBUG_MSG
      USE LOGICAL_MOD, ONLY : LPRT
      USE LOGICAL_MOD, ONLY : LLINOZ
      USE TRACER_MOD,  ONLY : ITS_A_FULLCHEM_SIM
      USE TRACER_MOD,  ONLY : ITS_A_TAGOX_SIM
      USE TRACER_MOD,  ONLY : ITS_A_H2HD_SIM
      USE LINOZ_MOD,   ONLY : DO_LINOZ
#if defined( DEVEL )
      USE GC_TYPE_MOD, ONLY : GC_MET_LOCAL
#endif

      USE CMN_SIZE_MOD  ! Size parameters
!
! !INPUT PARAMETERS: 
!
#if defined( DEVEL )
      TYPE(GC_MET_LOCAL), INTENT(IN) :: LOCAL_MET  ! Obj w/ met fields
#endif
! 
! !REVISION HISTORY: 
!  11 Mar 2003 - R. Yantosca - Initial version
!  (1 ) Removed IORD, JORD, KORD from the arg list.  Now references LPRT
!        from "logical_mod.f".  Now references ITS_A_FULLCHEM_SIM and
!        ITS_A_TAGOX_SIM from "tracer_mod.f" (bmy, 7/20/04)
!  (2 ) Now references ITS_A_H2HD_SIM from "tracer_mod.f".  Now call routine
!        UPBDFLX_HD for H2/HD simulation. (lyj, phs, 9/18/07)
!  (3 ) Added support for LINOZ (dbm, jliu, bmy, 10/16/09)
!  02 Dec 2010 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
      !=================================================================
      ! DO_UPBDFLX begins here!
      !=================================================================

      IF ( ITS_A_FULLCHEM_SIM() ) THEN

         !---------------
         ! Fullchem run
         !---------------

         ! Ox from strat 
         ! dbj changed for linoz
#if defined( DEVEL )
         IF ( LLINOZ ) THEN 
            CALL DO_LINOZ( LOCAL_MET )
         ELSE
            CALL UPBDFLX_O3( LOCAL_MET )
         ENDIF

         ! NOy from strat
         CALL UPBDFLX_NOY( 1, LOCAL_MET )
#else
         IF ( LLINOZ ) THEN 
            CALL DO_LINOZ
         ELSE
            CALL UPBDFLX_O3
         ENDIF

         ! NOy from strat
         CALL UPBDFLX_NOY( 1 )
#endif

      ELSE IF ( ITS_A_TAGOX_SIM() ) THEN

         !---------------
         ! Tagged Ox run
         !---------------

         ! Ox from strat
         ! dbj changed for linoz
#if defined( DEVEL )
         IF ( LLINOZ ) THEN 
            CALL DO_LINOZ( LOCAL_MET )
         ELSE
            CALL UPBDFLX_O3( LOCAL_MET )
         ENDIF
#else
         IF ( LLINOZ ) THEN 
            CALL DO_LINOZ
         ELSE
            CALL UPBDFLX_O3
         ENDIF
#endif

      ELSE IF ( ITS_A_H2HD_SIM() ) THEN

         !---------------
         ! H2/HD run
         !---------------

         ! HD form strat
#if defined( DEVEL )
         CALL UPBDFLX_HD( LOCAL_MET )
#else
         CALL UPBDFLX_HD
#endif

      ENDIF

      !### Debug
      IF ( LPRT ) CALL DEBUG_MSG( '### DO_UPBDFLX: after strat fluxes' )

      END SUBROUTINE DO_UPBDFLX
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: upbdflx_O3  
!
! !DESCRIPTION: Subroutine UPBDFLX\_O3 establishes the flux boundary condition 
!  for Ozone coming down from the stratosphere, using the Synoz algorithm of
!  McLinden et al, 2000. 
!\\
!\\
! !INTERFACE:
!
#if defined( DEVEL )
      SUBROUTINE UPBDFLX_O3( LOCAL_MET )
#else
      SUBROUTINE UPBDFLX_O3      
#endif
!
! !USES:
!
      USE DAO_MOD,        ONLY : AD, BXHEIGHT, T!, TROPP
      USE ERROR_MOD,      ONLY : ERROR_STOP
      USE LOGICAL_MOD,    ONLY : LVARTROP 
      USE PRESSURE_MOD,   ONLY : GET_PEDGE, GET_PCENTER
      USE TAGGED_OX_MOD,  ONLY : ADD_STRAT_POX
      USE TIME_MOD,       ONLY : GET_TS_DYN
      USE TRACER_MOD,     ONLY : STT, ITS_A_TAGOX_SIM
      USE TRACERID_MOD,   ONLY : IDTOX
      USE TROPOPAUSE_MOD, ONLY : GET_TPAUSE_LEVEL
#if defined( DEVEL )
      USE GC_TYPE_MOD,    ONLY : GC_MET_LOCAL
#endif
      
      USE CMN_SIZE_MOD       ! Size parameters
      USE CMN_GCTM_MOD       ! Rdg0
!
! !INPUT PARAMETERS: 
!
#if defined( DEVEL )
      TYPE(GC_MET_LOCAL), INTENT(IN) :: LOCAL_MET  ! Obj w/ met fields
#endif
!
! !REMARKS:
!  Reference:
!  ===========================================================================
!  C. A. McLinden, S. Olsen, B. Hannegan, O. Wild, M. J. Prather, and
!  J. Sundet, "Stratospheric Ozone in 3-D models: A simple chemistry
!  and the cross-tropopause flux".
!                                                                             .
!  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!  %%%% NOTE: This SYNOZ scheme is now obsolete, replaced by LINOZ   %%%%
!  %%%% We keep this for backwards compatibility w/ older met fields %%%%
!  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! 
! !REVISION HISTORY:
!  13 Dec 1999 - Q. Li, R. Martin - Initial version
!  (1 ) The parameter Rdg0 from "CMN_GCTM" = R / g0 = 28.97. 
!  (2 ) Pass PW = PS - PTOP to UPBDFLX via "CMN".
!  (3 ) Now pass IORD, JORD, KORD as arguments (bmy, 12/6/00)
!  (4 ) Now compute the proper value of PO3_vmr that will yield 475 Tg O3/yr
!        for various settings of IORD, JORD, KORD (rvm, bey, bmy, 12/5/00)
!                                                                             .
!        **************************************************************
!        ***** You must use this version of UPBDFLX_O3 if you are *****
!        ***** using the Parallel Processor TPCORE v. 7.1         *****
!        **************************************************************
!                                                                             .
!  (5 ) Added to "upbdflx_mod.f".  Also updated comments and made some
!        cosmetic changes. (bmy, 6/28/01)
!  (6 ) Now reference CMN_SETUP for LSPLIT.  Also store strat O3 into
!        tracer #11 for multi-tracer Ox run. (amf, bmy, 7/3/01)
!  (7 ) Removed IREF, JREF -- these are obsolete.  Also T(IREF,JREF,L) is
!        now T(I,J,L). (bmy, 9/27/01)
!  (8 ) Also replace PW(I,J) with P(I,J) (bmy, 10/3/01)
!  (9 ) Removed obsolete commented out code from 9/01 (bmy, 10/24/01)
!  (10) Removed obsolete commented out code from 7/01 (bmy, 11/26/01)
!  (11) Now write file names to stdout (bmy, 4/3/02)
!  (12) Replaced all instances of IM with IIPAR and JM with JJPAR, in order
!        to prevent namespace confusion for the new TPCORE (bmy, 6/25/02)
!  (13) Now use GET_PEDGE and GET_PCENTER from "pressure_mod.f" to compute
!        the pressure at the bottom edge and center of grid box (I,J,L).
!        Also removed obsolete, commented-out code.  Removed G_SIG and
!        G_SIGE from the arg list. (dsa, bdf, bmy, 8/21/02)
!  (14) Now reference BXHEIGHT and T from "dao_mod.f".  Also reference routine
!        ERROR_STOP from "error_mod.f".  Now references IDTOX from F90 module
!        "tracerid_mod.f" instead of from "comtrid.h". (bmy, 11/6/02)
!  (15) Now define J30S and J30N for 1x1 nested grid (bmy, 3/11/03)
!  (16) Make sure to pass AD via "dao_mod.f" for GEOS-1 (bnd, bmy, 4/14/03)
!  (17) On the first timestep, print how much O3 flux is coming down from the 
!        stratosphere in Tg/yr. (mje, bmy, 8/15/03)
!  (18) Change O3 flux to 500 Tg/yr for GEOS-3 (mje, bmy, 9/15/03)
!  (19) Now calls routine ADD_STRAT_POX from "tagged_ox_mod.f" in order to
!        pass stratospheric flux of Ox to the proper tagged tracer w/o
!        resorting to hardwiring w/in this routine. (bmy, 8/18/03)
!  (20) Add GEOS_4 to the #if defined block. (bmy, 1/29/04)
!  (21) Activated parallel DO-loops.  Now made STFLUX a local array
!        in order to facilitate parallelization. (bmy, 4/15/04)
!  (22) Removed IORD, JORD, KORD from the arg list.  Now reference STT and
!        ITS_A_TAGOX_SIM from "tracer_mod.f".  (bmy, 7/20/04)
!  (23) Use an #ifdef block to comment out an EXIT statement from w/in a
!        parallel loop for COMPAQ compiler.  COMPAQ seems to have some 
!        problems with this.  Now supports 1x125 grid. (auvray, bmy, 12/1/04)
!  (24) Now modified for GEOS-5 and GCAP met fields (swu, bmy, 5/25/05)
!  (25) Remove support for GEOS-1 and GEOS-STRAT met fields (bmy, 8/4/06)
!  (26) Now set J30S and J30N for GEOS-5 nested grid (yxw, dan, bmy, 11/6/08)
!  (27) Remove support for COMPAQ compiler (bmy, 7/8/09)
!  (28) Now do not call ADD_STRAT_POx for tagged Ox (dbj, bmy, 10/16/09)
!  13 Aug 2010 - R. Yantosca - Treat MERRA like GEOS-5 (bmy, 8/13/10)
!  02 Dec 2010 - R. Yantosca - Added ProTeX headers
!  08 Feb 2012 - R. Yantosca - Treat GEOS-5.7.2 in the same way as MERRA
!  10 Feb 2012 - R. Yantosca - Modified for 0.25 x 0.3125 grids
!  28 Feb 2012 - R. Yantosca - Removed support for GEOS-3
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL, SAVE        :: FIRST = .TRUE.
      INTEGER              :: I, J, L, L70mb
      INTEGER              :: NTRACER, NTRACE2
      REAL*8               :: P1, P2, P3, T1, T2, DZ, ZUP
      REAL*8               :: DTDYN, H70mb, PO3, PO3_vmr
      REAL*8               :: STFLUX(IIPAR,JJPAR,LLPAR)

      ! Select the grid boxes at the edges of the O3 release region, 
      ! for the proper model resolution (qli, bmy, 12/1/04)
#if   defined( GRID4x5 ) && defined( GCAP )
      ! GCAP has 45 latitudes, shift by 1/2 grid box (swu, bmy, 5/25/05)
      INTEGER, PARAMETER   :: J30S = 16, J30N = 30 

#elif defined( GRID4x5 )
      INTEGER, PARAMETER   :: J30S = 16, J30N = 31 

#elif defined( GRID2x25 ) 
      INTEGER, PARAMETER   :: J30S = 31, J30N = 61

#elif defined( GRID1x125 ) 
      INTEGER, PARAMETER   :: J30S = 61, J30N = 121

#elif defined( GRID05x0666 )
      INTEGER, PARAMETER   :: J30S = 1, J30N = JJPAR

#elif defined( GRID025x03125 )
      INTEGER, PARAMETER   :: J30S = 1, J30N = JJPAR

#elif defined( GRID1x1 ) 

#if   defined( NESTED_CH ) || defined( NESTED_NA )
      INTEGER, PARAMETER   :: J30S = 1,  J30N = JJPAR  ! 1x1 nested grids
#else  
      INTEGER, PARAMETER   :: J30S = 61, J30N = 121    ! 1x1 global grid
#endif

#elif defined( EXTERNAL_GRID )
      ! THIS HAS TO BE DEFINED SPECIFICALLY! HOW?
      INTEGER, PARAMETER   :: J30S = 31, J30N = 61

#endif

      ! Lower pressure bound for O3 release (unit: mb)
      ! REAL*8,  PARAMETER   :: P70mb = 70d0 !PHS
      REAL*8  :: P70mb, PTP

      !=================================================================
      ! UPBDFLX_O3 begins here!
      !=================================================================

      ! Dynamic timestep [s]
      DTDYN = GET_TS_DYN() * 60d0

      ! For O3 flux printout
      STFLUX = 0d0

      ! lower pressure !PHS
      P70mb = 70d0

      !=================================================================
      ! Compute the proper release rate of O3 coming down from the 
      ! stratosphere for the different GEOS model grids. 
      ! (bey, djj, rvm, bmy, 12/5/00).
      !
      ! PO3_vmr is the O3 release rate constant [v/v/s] that will yield 
      ! a strat-to-trop flux of 475 [Tg O3/yr].  Different TPCORE flags 
      ! create different amounts of ozone in the stratosphere.  Flags 
      ! 337F are currently preferred (bey, djj, rvm).  
      !  
      ! For now, provide values for PO3_vmr for two TPCORE flag settings:
      !  (1) IORD = 3, JORD = 3, KORD = 7  (preferred, assumed to 
      !                                     be the default)
      !  (2) IORD = 5, JORD = 5, KORD = 7 
      !=================================================================
#if   defined( GEOS_4 )
      PO3_vmr = 5.14d-14                                 ! 3,3,7

#elif defined( GEOS_5 ) || defined( MERRA ) || defined( GEOS_57 )

      ! For now assume GEOS-5 has same PO3_vmr value 
      ! as GEOS-4; we can redefine later (bmy, 5/25/05)
      PO3_vmr = 5.14d-14   

#elif defined( GCAP )

      ! For GCAP, assuming 3,3,7 (swu, bmy, 5/25/05)
      PO3_vmr = 5.0d-14 
      
#endif

      ! Store in the proper Ox tracer #
      NTRACER = IDTOX

      ! Only initialize on first time step
      IF ( FIRST ) STFLUX = 0d0

      ! Loop over latitude (30S -> 30N) and longitude
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I,  J,  L,  P2,  L70mb, P1, P3 )
!$OMP+PRIVATE( T2, T1, DZ, ZUP, H70mb, PO3    )
!$OMP+SCHEDULE( DYNAMIC )
      DO J = J30S, J30N 
      DO I = 1,    IIPAR

         !==============================================================
         ! L70mb is the 1st layer where pressure is equal to
         ! or smaller than 70 mb 
         !==============================================================
         
         !--------------------------------------------------------------
         ! Comment out for now (bmy, 10/2/07)
         ! replace L70mb with Tropopause pressure if the later is 
         ! lower -PHS #### still Beta testing ####
         !IF ( LVARTROP ) THEN
         !   PTP = TROPP(I,J)
         !   IF ( PTP < P70mb ) THEN
         !      P70mb = PTP
         !      !#### TESTING ####
         !      write(6,*)'#### RAISED bottom of O3 release region'
         !      write(6,*)'at ', i, j
         !      first=.true.
         !   ENDIF
         !ENDIF
         !--------------------------------------------------------------

         DO L = 1, LLPAR

            ! P2 = pressure [hPa] at the sigma center of level L70mb
            P2 = GET_PCENTER(I,J,L) 
               
            IF ( P2 < P70mb ) THEN
               L70mb = L
               EXIT
            ENDIF
         ENDDO
         
         ! P1 = pressure [hPa] at the sigma center of level L70mb - 1   
         P1 = GET_PCENTER(I,J,L70mb-1) 

         ! P3 = pressure [hPa] at the lower sigma edge of level L70mb
         P3 = GET_PEDGE(I,J,L70mb) 
 
         !==============================================================
         ! T2 = temperature (K)  at the sigma center of level L70mb
         ! T1 = temperature (K)  at the sigma center of level L70mb-1
         !
         ! DZ is the height from the sigma center of level L70mb-1 
         ! to 70mb.  Therefore, DZ may be found in either the 
         ! (L70mb)th sigma layer or the (L70mb-1)th sigma layer.  
         !
         ! ZUP is the height from the sigma center of the 
         ! (L70mb-1)th layer
         !============================================================== 
#if defined( DEVEL )
         T2   = LOCAL_MET%T(I,J,L70mb  )
         T1   = LOCAL_MET%T(I,J,L70mb-1)
#else
         T2   = T(I,J,L70mb  )
         T1   = T(I,J,L70mb-1)        
#endif

         DZ   = Rdg0 * ( (T1 + T2) / 2d0 ) * LOG( P1 / P70mb ) 
         ZUP  = Rdg0 * T1 * LOG( P1 /P3 )

         !==============================================================       
         ! H70mb is height between 70mb and the upper edge of the 
         ! level where DZ is.
         !  
         ! If DZ >= ZUP then DZ is already in level L70mb.
         ! If DZ <  ZUP then DZ is in level L70mb-1.
         !==============================================================       
         IF ( DZ >= ZUP ) THEN
#if defined( DEVEL )
            H70mb = LOCAL_MET%BXHEIGHT(I,J,L70mb) - ( DZ - ZUP )
#else
            H70mb = BXHEIGHT(I,J,L70mb) - ( DZ - ZUP )
#endif
         ELSE
            L70mb = L70mb - 1
            H70mb = ZUP - DZ
         ENDIF

         !=========================================================== 
         ! Distribute O3 into the region (30S-30N, 70mb-10mb)
         !=========================================================== 
         DO L = L70mb, LLPAR 

            ! Convert O3 in grid box (I,J,L) from v/v/s to kg/box
            PO3 = PO3_vmr * DTDYN 

#if   !defined( GCAP ) 
            ! For both 2 x 2.5 and 4 x 5 GEOS grids, 30S and 30 N are
            ! grid box centers.  However, the O3 release region is 
            ! edged by 30S and 30N.  Therefore, if we are at the 30S
            ! or 30N grid boxes, divide the O3 flux by 2.
            IF ( J == J30S .or. J == J30N ) THEN
               PO3 = PO3 / 2d0
            ENDIF 
#endif
            ! If we are in the lower level, compute the fraction
            ! of this level that lies above 70 mb, and scale 
            ! the O3 flux accordingly.
            IF ( L == L70mb ) THEN
#if defined( DEVEL )
               PO3 = PO3 * H70mb / LOCAL_MET%BXHEIGHT(I,J,L) 
#else
               PO3 = PO3 * H70mb / BXHEIGHT(I,J,L) 
#endif
            ENDIF
            
            ! Store O3 flux in the proper tracer number
            STT(I,J,L,NTRACER) = STT(I,J,L,NTRACER) + PO3 

            ! Store O3 flux for strat Ox tracer (Tagged Ox only)
            !IF ( ITS_A_TAGOX_SIM() ) CALL ADD_STRAT_POX( I, J, L, PO3 )
            !print*,'DBJ: COMMENTED ADD_STRAT_POX'

            ! Archive stratospheric O3 for printout in [Tg/yr]
            IF ( FIRST ) THEN
#if defined( DEVEL )
               STFLUX(I,J,L) = STFLUX(I,J,L) + 
     &              PO3   * LOCAL_MET%AD(I,J,L) * 1000.d0 / 28.8d0 /
     &              DTDYN * 48.d0  *  365.25d0  * 86400d0 / 1e12
#else
               STFLUX(I,J,L) = STFLUX(I,J,L) + 
     &              PO3 * AD(I,J,L) * 1000.d0 / 28.8d0 /
     &              DTDYN * 48.d0 * 365.25d0 * 86400d0 / 1e12
#endif
            ENDIF
         ENDDO   
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      !=================================================================
      ! Print amount of stratospheric O3 coming down
      !=================================================================
      IF ( FIRST ) THEN
         WRITE( 6, 100 ) SUM( STFLUX )
 100     FORMAT( '     - UPBDFLX_O3: Strat O3 production is', f9.3, 
     &           ' [Tg/yr]')
         FIRST = .FALSE.
      ENDIF

      END SUBROUTINE UPBDFLX_O3 
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: upbdflx_NOy
!
! !DESCRIPTION: Subroutine UPBDFLX\_NOy imposes NOy (NOx + HNO3) upper 
!  boundary condition in the stratosphere. The production rates for NOy are 
!  provided by Dylan Jones, along with NOx and HNO3 concentrations. 
!\\
!\\
! !INTERFACE:
!
#if defined( DEVEL )
      SUBROUTINE UPBDFLX_NOy( IFLAG, LOCAL_MET )
#else
      SUBROUTINE UPBDFLX_NOy( IFLAG )
#endif
!
! !USES:
!
      USE BPCH2_MOD,      ONLY : GET_NAME_EXT, GET_RES_EXT
      USE BPCH2_MOD,      ONLY : GET_TAU0,     READ_BPCH2
      USE DAO_MOD,        ONLY : AD
      USE DIRECTORY_MOD,  ONLY : DATA_DIR
      USE ERROR_MOD,      ONLY : ERROR_STOP, ALLOC_ERR
      USE TRACERID_MOD,   ONLY : IDTNOX,       IDTHNO3
      USE TIME_MOD,       ONLY : GET_TS_DYN,   GET_MONTH
      USE TIME_MOD,       ONLY : ITS_A_NEW_MONTH
      USE TRACER_MOD,     ONLY : STT,          XNUMOLAIR
      USE TRANSFER_MOD,   ONLY : TRANSFER_ZONAL
      USE TROPOPAUSE_MOD, ONLY : GET_MIN_TPAUSE_LEVEL
      USE TROPOPAUSE_MOD, ONLY : ITS_IN_THE_STRAT
#if defined( DEVEL )
      USE GC_TYPE_MOD,    ONLY : GC_MET_LOCAL
#endif

      USE CMN_SIZE_MOD       ! Size parameters
!
! !INPUT PARAMETERS: 
!
      ! IFLAG=1 will partition    [NOy] before transport
      ! IFLAG=2 will re-partition [NOy] after  transport
      INTEGER, INTENT(IN) :: IFLAG
#if defined( DEVEL )
      TYPE(GC_MET_LOCAL), INTENT(IN) :: LOCAL_MET  ! Obj w/ met fields
#endif
! 
! !REVISION HISTORY: 
!  22 Dec 1999 - Q. Li, R. Martin - Initial version
!  (1 ) Use READ_BPCH2 to read data from disk in binary punch file format.
!  (2 ) Now partition total [NOy] into [NOx] and [HNO3], instead of
!        partitioning P(NOy) into P(NOx) and P(HNO3).  (qli, bmy, 12/22/1999)
!  (3 ) Also echo back to the user when reading data from disk.  This 
!        allows the user to trace I/O errors more easily. (bmy, 2/1/00)
!  (4 ) Cosmetic changes, updated comments (bmy, 3/17/00)
!  (5 ) Reference F90 module "bpch2_mod" which contains routine "read_bpch2"
!        for reading data from binary punch files (bmy, 6/28/00)
!  (6 ) Only add P(NOy) above 10mb (archived in files "pnoy_above_10mb.*) 
!        into the top layer of the GEOS-1 and GEOS-STRAT grids.  The GEOS-2
!        and GEOS-3 grids extend well above 10mb and so they will contain
!        all of the P(NOy) up there (bmy, 6/29/00)
!  (7 ) Now use function GET_TAU0 (from "bpch2_mod.f") to return the TAU0 
!        value used to index the binary punch file. (bmy, 7/20/00)
!  (8 ) Only dump P(NOy) above 10mb for GEOS-1 grid.  The GEOS-STRAT grid
!        will already have this contribution, since it extends to 0.1 mb.
!        Also fix regridding error in P(NOy) data file.  Add parallel
!        processor DO-loops. (rvm, qli, bmy, 12/6/00)
!  (9 ) Now scale P(NOy) by 0.7 for TPCORE flags 337, in order to prevent
!        excess NOy from building up in the stratosphere. (rvm, bmy, 12/12/00)
!  (10) Now read properly regridded P(NOy) files from the pnoy_200106/
!        subdirectory of DATA_DIR.  Also updated comments and made a
!        few cosmetic changes. (mje, bmy, 6/28/01)
!  (11) Now use 3 arguments (M/D/Y) in call to GET_TAU0.  ARRAY needs to be 
!        of size (1,JJPAR,LGLOB).  Use JJPAR,LGLOB in calls to READ_BPCH2.
!        Use TRANSFER_ZONAL (from "transfer_mod.f") to cast from REAL*4 to
!        REAL*8 and resize arrays to (JJPAR,LLPAR) (bmy, 9/27/01)
!  (12) Removed obsolete commented out code from 9/01 (bmy, 10/24/01)
!  (13) Now write file name to stdout (bmy, 4/3/02)
!  (14) Now reference ERROR_STOP from "error_mod.f".  Also references IDTNOX
!        and IDTHNO3 from "tracerid_mod.f". (bmy, 11/6/02)
!  (15) Rename MONTHSAVE to LASTMONTH.  Now use functions GET_TS_DYN and
!        GET_MONTH from "time_mod.f".  Now call READ_BPCH2 with QUIET=.TRUE.
!        to suppress printing of extra info.  Cosmetic changes.  Now references
!        AD from "dao_mod.f" for GEOS-1 (bmy, 4/14/03)
!  (16) Activated parallel DO-loops.  Moved the computation of XRATIO into
!        the IF block which only gets done once per month. (bmy, 4/15/04)
!  (17) Now references STT from "tracer_mod.f".  Now references DATA_DIR 
!        from "directory_mod.f". (bmy, 7/20/04)
!  (18) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (19) Now references XNUMOLAIR from "tracer_mod.f" (bmy, 10/25/05)
!  (20) Now references ITS_A_NEW_MONTH from "time_mod.f". Now reference 
!        GET_MIN_TPAUSE_LEVEL from "tropopause_mod.f".  Now replace reference 
!        to LPAUSE with ITS_IN_THE_STRAT from "tropopause_mod.f" (bmy, 11/1/05)
!  (21) Remove support for GEOS-1 and GEOS-STRAT met fields (bmy, 8/4/06)
!  (22) Cap 1-XRATIO to avoid numerical problems later (bmy, 6/30/08)
!  02 Dec 2010 - R. Yantosca - Added ProTeX headers
!  07 Mar 2012 - R. Yantosca - Start DO loops from L=1 instead of L=LMIN
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER              :: I, J, L, LMIN, AS
      INTEGER, SAVE        :: LASTMONTH = -99
      LOGICAL, SAVE        :: FIRST = .TRUE.

      REAL*4               :: DTDYN, AIRDENS, PNOY
      REAL*4               :: ARRAY(1,JJPAR,LGLOB)

#if defined( DEVEL )
      ! Ratio of ( [NO] + [NO2] ) / [NOy] 
      REAL*4, ALLOCATABLE, SAVE  :: XRATIO(:,:)

      ! Arrays for P(NOY), NO, NO2, and HNO3 concentrations
      REAL*4, ALLOCATABLE, SAVE  :: STRATPNOY(:,:)
      REAL*4, ALLOCATABLE, SAVE  :: STRATNO(:,:)
      REAL*4, ALLOCATABLE, SAVE  :: STRATNO2(:,:)
      REAL*4, ALLOCATABLE, SAVE  :: STRATHNO3(:,:)

      ! For P(NOy) above 10 mb
      REAL*4, ALLOCATABLE, SAVE  :: SPNOY10mb(:)
#else
      ! Ratio of ( [NO] + [NO2] ) / [NOy] 
      REAL*4, SAVE  :: XRATIO(JJPAR,LLPAR)

      ! Arrays for P(NOY), NO, NO2, and HNO3 concentrations
      REAL*4, SAVE  :: STRATPNOY(JJPAR,LLPAR)
      REAL*4, SAVE  :: STRATNO(JJPAR,LLPAR)
      REAL*4, SAVE  :: STRATNO2(JJPAR,LLPAR)
      REAL*4, SAVE  :: STRATHNO3(JJPAR,LLPAR)

      ! For P(NOy) above 10 mb
      REAL*4, SAVE  :: SPNOY10mb(JJPAR)
#endif

      ! TAU values for indexing the punch file 
      REAL*8               :: XTAU

      ! File Names
      CHARACTER (LEN=255)  :: FILENAME
      CHARACTER (LEN=255)  :: FILENAME2

      !=================================================================
      ! UPBDFLX_NOY begins here! 
      !=================================================================

#if defined( DEVEL )
      IF (FIRST) THEN
         ALLOCATE( STRATPNOY(JJPAR,LLPAR), STAT=AS )
         IF (AS /=0) CALL ALLOC_ERR( 'STRATPNOY' )
         ALLOCATE( STRATNO(JJPAR,LLPAR), STAT=AS )
         IF (AS /=0) CALL ALLOC_ERR( 'STRATNO' )
         ALLOCATE( STRATNO2(JJPAR,LLPAR), STAT=AS )
         IF (AS /=0) CALL ALLOC_ERR( 'STRATNO2' )
         ALLOCATE( STRATHNO3(JJPAR,LLPAR), STAT=AS )
         IF (AS /=0) CALL ALLOC_ERR( 'STRATHNO3' )
         ALLOCATE( XRATIO(JJPAR,LLPAR), STAT=AS)
         IF (AS /=0) CALL ALLOC_ERR( 'XRATIO' )
         ALLOCATE( SPNOY10mb(JJPAR), STAT=AS)
         IF (AS /=0) CALL ALLOC_ERR( 'SPNOY10mb' )
         
         FIRST = .FALSE.
      ENDIF
#endif

      ! Dynamic timestep [s]
      DTDYN = GET_TS_DYN() * 60d0

      !=================================================================
      ! IFLAG = 1: Before transport
      !
      ! If we have entered into a new month, read P(NOy), HNO3,
      ! NO, and NO2 from disk (binary punch file format).
      !=================================================================
      IF ( IFLAG == 1 ) THEN

         IF ( ITS_A_NEW_MONTH() ) THEN

            ! TAU value corresponding to the beginning of this month
            XTAU = GET_TAU0( GET_MONTH(), 1, 1985 ) 

            ! File containing P(NOy), NOx, HNO3 concentrations
            ! Now read corrected file from pnoy_200106/ subdir (bmy, 6/28/01)
            FILENAME  = TRIM( DATA_DIR )             // 
     &                  'pnoy_200106/pnoy_nox_hno3.' //
     &                  GET_NAME_EXT()   // '.'      // GET_RES_EXT()

            ! Echo filename to stdout
            WRITE( 6, 100 ) TRIM( FILENAME )
 100        FORMAT( '     - UPBDFLX_NOY: Reading ', a )
            
            ! P(NOy) in [v/v/s] is stored as tracer #1
            CALL READ_BPCH2( FILENAME, 'PNOY-L=$', 1,     
     &                       XTAU,      1,         JJPAR,     
     &                       LGLOB,     ARRAY,     QUIET=.TRUE. )

            ! Cast from REAL*4 to REAL*8 and resize to (JJPAR,LLPAR)
            CALL TRANSFER_ZONAL( ARRAY(1,:,:), STRATPNOY )

            ! [HNO3] in [v/v] is stored as tracer #2
            CALL READ_BPCH2( FILENAME, 'PNOY-L=$', 2,     
     &                       XTAU,      1,         JJPAR,     
     &                       LGLOB,     ARRAY,     QUIET=.TRUE. )

            ! Cast from REAL*4 to REAL*8 and resize to (JJPAR,LLPAR)
            CALL TRANSFER_ZONAL( ARRAY(1,:,:), STRATHNO3 )
         
            ! [NO] in [v/v] is stored as tracer #4
            CALL READ_BPCH2( FILENAME, 'PNOY-L=$', 4,     
     &                       XTAU,      1,         JJPAR,     
     &                       LGLOB,     ARRAY,     QUIET=.TRUE. )

            ! Cast from REAL*4 to REAL*8 and resize to (JJPAR,LLPAR)
            CALL TRANSFER_ZONAL( ARRAY(1,:,:), STRATNO )

            ! [NO2] in [v/v] is stored as tracer #5
            CALL READ_BPCH2( FILENAME, 'PNOY-L=$', 5,     
     &                       XTAU,      1,         JJPAR,     
     &                       LGLOB,     ARRAY,     QUIET=.TRUE. )

            ! Cast from REAL*4 to REAL*8 and resize to (JJPAR,LLPAR)
            CALL TRANSFER_ZONAL( ARRAY(1,:,:), STRATNO2 )

            !===========================================================
            ! XRATIO is the ratio ( [NO] + [NO2] ) / [NOy],
            ! which is needed for the partitioning.
            ! XRATIO will be the same for a given month
            !===========================================================
            DO L = 1, LLPAR
            DO J = 1, JJPAR
               XRATIO(J,L) = ( STRATNO(J,L) + STRATNO2(J,L) ) /
     &                       ( STRATNO(J,L) + STRATNO2(J,L) + 
     &                         STRATHNO3(J,L) ) 
            ENDDO
            ENDDO
         ENDIF

         !==============================================================
         ! Initial partitioning of [NOy] to [NOx] and [HNO3], before 
         ! transport
         ! 
         ! We use zonal mean values for stratospheric P(NOy), [NO], 
         ! [NO2], and [HNO3] taken from Dylan Jones' & Hans Schneider's 
         ! 2-D model.  
         !
         ! Since P(NOy) above 10mb accounts for almost 50% of the total 
         ! stratospheric production, we also dump P(NOy) above 10 mb 
         ! into the top layer of the model.  These values are also 
         ! supplied to us by Dylan Jones.
         !
         ! We make the following assumptions:
         !
         !    (1) [NOx] = [NO] + [NO2] 
         !    (2) [NOy] = [NO] + [NO2] + [HNO3] = [NOx] + [HNO3]  
         !
         ! Therefore, in order to obtain [NOx] and [HNO3] from [NOy], 
         ! we must do the partitioning as follows:
         !
         !    (1) [NOy]   = P(NOy) + [NOx] + [HNO3]  
         !                = Production of NOy plus current 
         !                  concentrations of NOx and HNO3 in the 
         !                  given grid box
         !
         !    (2) XRATIO  = ( [NO] + [NO2] ) / [NOy]
         !
         !    (3) P(NOx)  = P(NOy) * XRATIO
         !  
         !    (4) P(HNO3) = P(NOy) * ( 1 - XRATIO ) 
         !
         ! XRATIO = ( [NO] + [NO2] ) / [NOy] approximates the true 
         ! ratio of [NOx] / [NOy], but is itself not the true ratio, 
         ! since our formulation of [NOy] neglects some additional 
         ! species (e.g. PAN, HNO4, N2O5, R4N2, PPN, PMN).
         !
         ! At some future point we may take the additional constituents 
         ! of [NOy] into account.  For now we proceed as outlined above.
         !==============================================================

         ! Minimum value of LPAUSE
         LMIN = GET_MIN_TPAUSE_LEVEL()

!$OMP PARALLEL DO
#if defined( DEVEL )
!$OMP+DEFAULT( SHARED )
#else
!$OMP+DEFAULT( SHARED ) SHARED(XRATIO)
#endif
!$OMP+PRIVATE( I, J, L, PNOY )
!$OMP+SCHEDULE( DYNAMIC )
         DO L = 1,    LLPAR
         DO J = 1,    JJPAR
         DO I = 1,    IIPAR

            ! Skip over tropospheric boxes
#if defined( DEVEL )
            IF ( ITS_IN_THE_STRAT( I, J, L, LOCAL_MET ) ) THEN 
#else
            IF ( ITS_IN_THE_STRAT( I, J, L ) ) THEN 
#endif

               ! PNOY = P(NOy) converted from [v/v/s] to [v/v] 
               PNOY = STRATPNOY(J,L) * DTDYN 

               ! Add [NOx] and [HNO3] to PNOY.
               ! PNOY is now the total [NOy] concentration
               PNOY = PNOY + STT(I,J,L,IDTNOX) + STT(I,J,L,IDTHNO3)

               ! Partition total [NOy] to [NOx], units are [v/v]
               STT(I,J,L,IDTNOX)  = PNOY * XRATIO(J,L) 

               ! Partition total [NOy] to [HNO3], units are [v/v]
               STT(I,J,L,IDTHNO3) = PNOY *
     &                              MAX( ( 1d0 -  XRATIO(J,L) ), 1d-20 )

            ENDIF
         ENDDO   
         ENDDO
         ENDDO
!$OMP END PARALLEL DO

      !=================================================================
      ! IFLAG = 2: After transport
      !
      ! Repartition [NOy] after transport into [NOx] + [HNO3]
      !
      ! This repartitioning is necessary to avoid performing chemistry
      ! between the [NO2] and [HNO3] species.
      !
      ! The concentrations [NOx] and [HNO3] will have changed due to 
      ! transport, but the ratio used to partition them will be the 
      ! same.
      !=================================================================
      ELSE IF ( IFLAG == 2 ) THEN

         ! Minimum value of LPAUSE
         LMIN = GET_MIN_TPAUSE_LEVEL()

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED ) SHARED(XRATIO)
!$OMP+PRIVATE( I, J, L, PNOY )
!$OMP+SCHEDULE( DYNAMIC )
!--------------------------------
!         DO L = LMIN, LLPAR
!-----------------------------
         DO L = 1,    LLPAR
         DO J = 1,    JJPAR
         DO I = 1,    IIPAR

            ! Skip over tropospheric boxes
#if defined( DEVEL )
            IF ( ITS_IN_THE_STRAT( I, J, L, LOCAL_MET ) ) THEN
#else
            IF ( ITS_IN_THE_STRAT( I, J, L ) ) THEN
#endif

               ! Compute the new total [NOy] by summing up [NOx] + [HNO3]
               PNOY = STT(I,J,L,IDTNOX) + STT(I,J,L,IDTHNO3)

               ! Partition total [NOy] to [NOx], units are [v/v]
               STT(I,J,L,IDTNOX)  = PNOY * XRATIO(J,L) 

               ! Partition total [NOy] to [HNO3], units are [v/v]
               STT(I,J,L,IDTHNO3) = PNOY *
     &                              MAX( ( 1d0 -  XRATIO(J,L) ), 1d-20 )

            ENDIF
         ENDDO   
         ENDDO
         ENDDO   
!$OMP END PARALLEL DO      

      ELSE 

         ! If IFLAG /= 1 or IFLAG /= 2, print an error message and stop
         CALL ERROR_STOP( 'IFLAG must be 1 or 2!',
     &                    'UPBDFLX_NOY (upbdflx_mod.f)' )

      ENDIF

      END SUBROUTINE UPBDFLX_NOy
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: upbdflx_HD
!
! !DESCRIPTION: Subroutine UPBDFLX\_HD establishes the flux boundary condition 
!  for HD coming down from the stratosphere. This is adapted from the 
!  UPBDFLX\_O3 routine.
!\\
!\\
! !INTERFACE:
!
#if defined( DEVEL )
      SUBROUTINE UPBDFLX_HD( LOCAL_MET )
#else
      SUBROUTINE UPBDFLX_HD
#endif
!
! !USES:
!
      USE DAO_MOD,      ONLY : AD, BXHEIGHT, T
      USE ERROR_MOD,    ONLY : ERROR_STOP
      USE PRESSURE_MOD, ONLY : GET_PEDGE, GET_PCENTER
      USE TIME_MOD,     ONLY : GET_TS_DYN
      USE TRACER_MOD,   ONLY : STT
      USE TRACERID_MOD, ONLY : IDTHD, IDTH2
#if defined( DEVEL )
      USE GC_TYPE_MOD,  ONLY : GC_MET_LOCAL
#endif
      
      USE CMN_SIZE_MOD     ! Size parameters
      USE CMN_GCTM_MOD     ! Rdg0
!
! !INPUT PARAMETERS:
!
#if defined( DEVEL )
      TYPE(GC_MET_LOCAL), INTENT(IN) :: LOCAL_MET  ! Obj w/ met fields
#endif
!
! !REMARKS:
!  Instead of calculating the fractionation of H2 in the stratosphere 
!  (where we would have to take into account fractionation of CH4),
!  we simply set the HD tracer concentrations in the stratosphere to
!  reproduce observed profiles in the UT/LS.
!                                                                             .
!  References:
!  ===========================================================================
!  (1) "Global Budget of Molecular Hydrogen and its Deuterium Content: 
!        Constraints from Ground Station, Cruise, and Aircraft Observations" 
!        Price, H., L. Jaeglé, A. Rice, P. Quay, P.C. Novelli, R. Gammon, 
!        submitted to J. Geophys. Res., 2007.
! 
! !REVISION HISTORY: 
!  18 Sep 2007 - L. Jaegle, H. U. Price, P. Le Sager - Initial version
!  (1 ) First adapted from UPBDFLX_O3 (G-C v5-05-03) then merged w/ v7-04-12.
!        Added parallel DO loops. (phs, 9/18/07)
!  (26) Now set J30S and J30N for GEOS-5 nested grid (yxw, dan, bmy, 11/6/08)
!  (27) Remove support for COMPAQ compiler (bmy, 7/8/09)
!  13 Aug 2010 - R. Yantosca - Treat MERRA like GEOS-5
!  02 Dec 2010 - R. Yantosca - Added ProTeX headers
!  08 Feb 2012 - R. Yantosca - Treat GEOS-5.7.2 in the same way as MERRA
!  10 Feb 2012 - R. Yantosca - Modified for 0.25 x 0.3125 grids
!  28 Feb 2012 - R. Yantosca - Removed support for GEOS-3
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER              :: I, J, L, L70mb
      INTEGER              :: NTRACER
      REAL*8               :: P1, P2, P3, T1, T2, DZ, ZUP
      REAL*8               :: DTDYN, H70mb, PO3_vmr!,PO3
      REAL*8               :: PHD, PHD_vmr, SCALE_HD!, HD_AVG

      ! Select the grid boxes at the edges of the HD release region, 
      ! for the proper model resolution 
#if   defined( GRID4x5 ) && defined( GCAP )
      ! GCAP has 45 latitudes, shift by 1/2 grid box (swu, bmy, 5/25/05)
      INTEGER, PARAMETER   :: J30S = 16, J30N = 30 

#elif defined( GRID4x5 )
      INTEGER, PARAMETER   :: J30S = 16, J30N = 31 

#elif defined( GRID2x25 ) 
      INTEGER, PARAMETER   :: J30S = 31, J30N = 61

#elif defined( GRID1x125 ) 
      INTEGER, PARAMETER   :: J30S = 61, J30N = 121

#elif defined( GRID05x0666 )
      INTEGER, PARAMETER   :: J30S = 1, J30N = JJPAR

#elif defined( GRID025x03125 )
      INTEGER, PARAMETER   :: J30S = 1, J30N = JJPAR

#elif defined( GRID1x1 ) 

#if   defined( NESTED_CH ) || defined( NESTED_NA )
      INTEGER, PARAMETER   :: J30S = 1,  J30N = JJPAR  ! 1x1 nested grids
#else  
      INTEGER, PARAMETER   :: J30S = 61, J30N = 121    ! 1x1 global grid
#endif

#elif defined( EXTERNAL_GRID )
      ! THIS HAS TO BE DEFINED SPECIFICALLY! HOW?
      INTEGER, PARAMETER   :: J30S = 31, J30N = 61

#endif

      ! Lower pressure bound for HD release (unit: mb)
      REAL*8,  PARAMETER   :: P70mb = 70d0

      !=================================================================
      ! UPBDFLX_HD begins here!
      !=================================================================

      ! Dynamic timestep [s]
      DTDYN = GET_TS_DYN() * 60d0

      !=================================================================
      ! For now the only HD release rates are for GEOS-3. This will
      ! likely need to be scaled with other met fields (GEOS-4,
      ! GEOS-5...) jaegle, 2/20/2007
      ! PO3_vmr is the release rate for Ozone. This is then scaled by
      ! SCALE_HD in order to obtain the HD profile, to obtain
      ! PHD_vmr [v/v/s]
      ! For now uses the GEOS-3 scale factor for all cases (phs)
      !=================================================================
#if   defined( GEOS_4 )

      PO3_vmr = 5.14d-14                                 ! 3,3,7

#elif defined( GEOS_5 ) || defined( MERRA ) || defined( GEOS_57 )

      ! For now assume GEOS-5 has same PO3_vmr value 
      ! as GEOS-4; we can redefine later (bmy, 5/25/05)
      PO3_vmr = 5.14d-14   

#elif defined( GCAP )

      ! For GCAP, assuming 3,3,7 (swu, bmy, 5/25/05)
      PO3_vmr = 5.0d-14 
      
#endif

      ! Define scaling factor for HD and scale PO3_vmr
      ! Standard:
      SCALE_HD = 4.0d-5
      PHD_vmr= PO3_vmr * SCALE_HD

      !=================================================================
      ! Select the proper tracer number to store HD into
      !=================================================================
      NTRACER = IDTHD

      !=================================================================
      ! Loop over latitude/longtitude locations (I,J)
      !=================================================================
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I,  J,  L,  P2,  L70mb, P1, P3 )
!$OMP+PRIVATE( T2, T1, DZ, ZUP, H70mb, PHD    )
!$OMP+SCHEDULE( DYNAMIC )
       DO J = J30S, J30N 
         DO I = 1, IIPAR

            !===========================================================
            ! L70mb is the 1st layer where pressure is equal to
            ! or smaller than 70 mb 
            !
            ! P1 = pressure [ mb ] at the sigma center     of level L70mb - 1
            ! P3 = pressure [ mb ] at the lower sigma edge of level L70mb
            ! P2 = pressure [ mb ] at the sigma center     of level L70mb
            !===========================================================
            DO L = 1, LLPAR
               P2 = GET_PCENTER(I,J,L) 
               
               IF ( P2 < P70mb ) THEN
                  L70mb = L
                  EXIT
               ENDIF
            ENDDO
            
            P1 = GET_PCENTER(I,J,L70mb-1) 
            P3 = GET_PEDGE(I,J,L70mb) 
 
            !===========================================================
            ! T2 = temperature (K)  at the sigma center of level L70mb
            ! T1 = temperature (K)  at the sigma center of level L70mb-1
            !
            ! DZ is the height from the sigma center of level L70mb-1 
            ! to 70mb.  Therefore, DZ may be found in either the 
            ! (L70mb)th sigma layer or the (L70mb-1)th sigma layer.  
            !
            ! ZUP is the height from the sigma center of the 
            ! (L70mb-1)th layer
            !=========================================================== 
#if defined( DEVEL )
            T2   = LOCAL_MET%T(I,J,L70mb  )
            T1   = LOCAL_MET%T(I,J,L70mb-1)
#else
            T2   = T(I,J,L70mb  )
            T1   = T(I,J,L70mb-1)        
#endif

            DZ   = Rdg0 * ( (T1 + T2) / 2d0 ) * LOG( P1 / P70mb ) 
            ZUP  = Rdg0 * T1 * LOG( P1 /P3 )

            !===========================================================       
            ! H70mb is height between 70mb and the upper edge of the 
            ! level where DZ is.
            !  
            ! If DZ >= ZUP then DZ is already in level L70mb.
            ! If DZ <  ZUP then DZ is in level L70mb-1.
            !===========================================================       
            IF ( DZ >= ZUP ) THEN
#if defined( DEVEL )
               H70mb = LOCAL_MET%BXHEIGHT(I,J,L70mb) - ( DZ - ZUP )
#else
               H70mb = BXHEIGHT(I,J,L70mb) - ( DZ - ZUP )
#endif
            ELSE
               L70mb = L70mb - 1
               H70mb = ZUP - DZ
            ENDIF

            !=========================================================== 
            ! Distribute HD into the region (30S-30N, 70mb-10mb)
            !=========================================================== 
            DO L = L70mb, LLPAR 

               ! Convert HD in grid box (I,J,L) from v/v/s to kg/box
               PHD = PHD_vmr * DTDYN 

#if   !defined( GCAP ) 
               ! For both 2 x 2.5 and 4 x 5 GEOS grids, 30S and 30 N are
               ! grid box centers.  However, the O3 release region is 
               ! edged by 30S and 30N.  Therefore, if we are at the 30S
               ! or 30N grid boxes, divide the O3 flux by 2.
               IF ( J == J30S .or. J == J30N ) THEN
                  PHD = PHD / 2d0
               ENDIF 
#endif
               ! If we are in the lower level, compute the fraction
               ! of this level that lies above 70 mb, and scale 
               ! the HD flux accordingly.
               IF ( L == L70mb ) THEN
#if defined( DEVEL )
                  PHD = PHD * H70mb / LOCAL_MET%BXHEIGHT(I,J,L) 
#else
                  PHD = PHD * H70mb / BXHEIGHT(I,J,L) 
#endif
               ENDIF
            
               ! Store HD flux in the proper tracer number
               STT(I,J,L,NTRACER) = STT(I,J,L,NTRACER) + PHD

            ENDDO   
         ENDDO
      ENDDO
!$OMP END PARALLEL DO

      END SUBROUTINE UPBDFLX_HD
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_upbdflx
!
! !DESCRIPTION: Subroutine INIT\_UPBDFLX passes IORD, JORD, and KORD values 
!  from "input\_mod.f" to "upbdflx\_mod.f"
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE INIT_UPBDFLX( I_ORD, J_ORD, K_ORD )
!
! !INPUT PARAMETERS: 
!
      INTEGER, INTENT(IN) :: I_ORD   ! TPCORE IORD parameter (for E/W)
      INTEGER, INTENT(IN) :: J_ORD   ! TPCORE JORD parameter (for N/S)
      INTEGER, INTENT(IN) :: K_ORD   ! TPCORE KORD parameter (for Vertical)
! 
! !REVISION HISTORY: 
!  20 Jul 2004 - R. Yantosca - Initial version
!  02 Dec 2010 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
      !=================================================================
      ! SET_UPBDFLX begins here!
      !=================================================================
      IORD = I_ORD
      JORD = J_ORD
      KORD = K_ORD 

      END SUBROUTINE INIT_UPBDFLX
!EOC
      END MODULE UPBDFLX_MOD
