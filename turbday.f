! $Id: turbday.f,v 1.2 2003/07/08 15:30:23 bmy Exp $
      SUBROUTINE TURBDAY( LTURB, XTRA2, NTRC, TC, TCVV )              
!
!******************************************************************************
!  Subroutine TURBDAY executes the GEOS-CTM dry convection / boundary layer 
!  mixing algorithm.  Original subroutine by Dale Allen, Univ of MD.
!  (bmy, bey, 1/30/98, 6/23/03)
!
!  If LTURB = T then 
!     (a) Compute XTRA2(:,:)
!     (b) Do the boundary layer mixing and diagnostics
!
!  IF LTURB = F then
!     (a) Compute XTRA2(:,:) and exit
!
!  XTRA2 must be computed even if the BL mixing is not done.
!
!  Arguments as Input:
!  ===========================================================================
!  (1 ) LTURB   : LTURB=T will do the BL mixing          [  T or F  ]
!  (2 ) NTRC    : Number of tracers used in computation  [1 to NNPAR]
!  (3 ) TC      : Tracer concentration                   [ v / v    ]
!  (4 ) TCVV    : MW air (g/mol) / MW tracer (g/mol)     [ unitless ]
!
!  Arguments as output:
!  ==========================================================================
!  (5 ) XTRA2   : height of PBL                          [ # layers ]
!  (3 ) TC      : Tracer concentration                   [ v / v    ]
!
!  NOTES:
!  (1 ) TURBDAY is written in Fixed-Form Fortran 90.  Also use F90
!        syntax for declarations (bmy, 4/1/99).
!  (2 ) New tracer concentrations are returned in TC.
!  (3 ) PS(I,J) is ACTUAL surface pressure and not Psurface - PTOP 
!  (4 ) Change in tracer in kg is now stored in DTC(I,J,L,N).  This makes
!        it easier to compute diagnostic quantities.  The new mixing ratio
!        is computed as TC(I,J,L,N) = TC(I,J,L,N) + DTC(I,J,L,N) / AD(I,J,L).
!  (5 ) XTRA2(*,*,5) is the height of the PBL in # of layers.  So if the
!        PBL top is located in the middle of the 3rd sigma layer at (I,J)
!        the value of XTRA2(I,J,5) would be 2.5.  The XTRA2 variable is
!        used by the HCTM drydep subroutines...it really is a historical
!        holdover.
!  (6 ) Restore the following NDxx diagnostics: (a) ND63 : Mass balance 
!        (CNVUPP) (b) ND15 : Mass change due to mixing in the boundary layer 
!  (7 ) Now pass TCVV and NCONV for the mass flux diagnostics.  Also
!        updated comments and cleaned up a few things. (bey, bmy, 11/10/99)
!  (8 ) Remove PTOP and XNUMOL from the arg list.  PTOP is now a parameter 
!        in "CMN_SIZE".  XNUMOL is no longer used in TURBDAY. (bmy, 2/10/00)
!  (9 ) Also removed obsolete ND63 diagnostics and updated comments.
!        (bmy, 4/12/00)
!  (10) Now use NTRC instead of NNPAR to dimension variables TC, TCVV, DTC, 
!        and DTCSUM (bmy, 10/17/00). 
!  (11) Removed obsolete code from 10/17/00 (bmy, 12/21/00)
!  (12) If the PBL depth is very small (or zero), then assume a PBL depth 
!        of 2 mb -- this prevents NaN's from propagating throughout the
!        code.  Also updated comments & made cosmetic changes. (bmy, 3/9/01)
!  (13) DTCSUM was declared twice but wasn't used.  Elminate declarations
!        to DTCSUM. (bmy, 7/16/01)
!  (14) XTRA2(IREF,JREF,5) is now XTRA2(I,J).  Also updated comments. 
!        Also remove IREF, JREF and some debug output. (bmy, 9/25/01)
!  (15) Removed obsolete commented out code from 9/01 (bmy, 10/24/01)
!  (16) Now takes in P=PS-PTOP instead of PS.  Redimension SIGE to 
!        (1:LLPAR+1).  
!  (17) Renamed PS to PZ so as not to conflict w/ the existing P variable.
!        Now pass P-PTOP thru PZ, in order to ensure that P and AD are
!        consistent w/ each other.  Added parallel DO-loops. Updated comments,
!        cosmetic changes.  Now print a header to stdout on the first call,
!        to confirm that TURBDAY has been called. (bmy, 4/11/02)
!  (18) Now use GET_PEDGE from "pressure_mod.f" to compute the pressure
!        at the bottom edge of grid box (I,J,L).  Deleted obsolete code from 
!        4/02.  Removed PZ, SIGE from the argument list, since we now compute
!        pressure from GET_PEDGE. (dsa, bdf, bmy, 8/22/02)  
!  (19)	Now reference AD, PBL from "dao_mod.f".  Now removed DXYP from the 
!        arg list, use GET_AREA_M2 from "grid_mod.f" instead.  Now removed 
!        NCONV, ALPHA_d, ALPHA_n from the arg list.  Now no longer reference 
!        SUNCOS.  Now set A(:,:)=1 day & nite; we assume full mixing all the 
!        time regardless of SUNCOS.  Updated comments, cosmetic changes.
!        (bmy, 2/11/03)
!  (20) Now can handle PBL field in meters for GEOS-4/fvDAS.  Also the
!        atmospheric scale height from CMN_GCTM. (bmy, 6/23/03)
!******************************************************************************
!
      ! References to F90 modules 
      USE DAO_MOD,      ONLY : AD, PBL
      USE DIAG_MOD,     ONLY : TURBFLUP
      USE GRID_MOD,     ONLY : GET_AREA_M2
      USE TIME_MOD,     ONLY : GET_TS_CONV
      USE PRESSURE_MOD, ONLY : GET_PEDGE

      IMPLICIT NONE

#     include "CMN_SIZE"
#     include "CMN_DIAG"
#     include "CMN_GCTM"

      ! Arguments
      LOGICAL,  INTENT(IN)    :: LTURB
      INTEGER,  INTENT(IN)    :: NTRC
      REAL*8,   INTENT(OUT)   :: XTRA2(IIPAR,JJPAR)
      REAL*8,   INTENT(INOUT) :: TC(IIPAR,JJPAR,LLPAR,NTRC)
      REAL*8,   INTENT(IN)    :: TCVV(NTRC)
      
      ! Local variables
      LOGICAL, SAVE           :: FIRST = .TRUE.
      INTEGER                 :: I, J, L, LTOP, N
      INTEGER                 :: IMIX(IIPAR,JJPAR) 
      REAL*8                  :: AA,  CC, CC_AA,   BLTOP 
      REAL*8                  :: PW,  PS, AREA_M2, DTCONV
      REAL*8                  :: P(0:LLPAR)
      REAL*8                  :: A(IIPAR,JJPAR)
      REAL*8                  :: FPBL(IIPAR,JJPAR)
      REAL*8                  :: DTC(IIPAR,JJPAR,LLPAR,NTRC)  

      !=================================================================
      ! TURBDAY begins here!          
      !=================================================================

      ! First-time initialization
      IF ( FIRST .and. LTURB ) THEN

         ! Echo info
         WRITE( 6, '(a)' ) REPEAT( '=', 79 )
         WRITE( 6, '(a)' ) 'T U R B D A Y  -- by Dale Allen, U. Md.'
         WRITE( 6, '(a)' ) 'Modified for GEOS-CHEM by Bob Yantosca'
         WRITE( 6, '(a)' ) 'Last Modification Date: 2/4/03'
         WRITE( 6, '(a)' ) REPEAT( '=', 79 )

         ! Reset first time flag
         FIRST = .FALSE.
      ENDIF
      
      !=================================================================
      ! Compute boundary layer height
      !=================================================================
      DO J = 1, JJPAR
      DO I = 1, IIPAR
         
         ! Surface pressure [hPa]
         PS    = GET_PEDGE(I,J,1)
         P(0)  = PS

#if   defined( GEOS_4 )

         ! BLTOP = pressure at PBL top [hPa]
         ! Use barometric law since PBL is in [m]
         BLTOP = PS * EXP( -PBL(I,J) / SCALE_HEIGHT )

#else

         ! BLTOP = pressure at PBL top
         ! PBL is in [hPa], so subtract it from surface pressure [hPa]
         BLTOP = PS - PBL(I,J)
         
         ! If the PBL depth is very small (or zero), then assume
         ! a PBL depth of 2 mb.  This will prevent NaN's from
         ! propagating throughout the code. (bmy, 3/7/01)
         IF ( PBL(I,J) < 1d-5 ) BLTOP = PS - 2d0

#endif

         ! LTOP is the GEOS-CHEM level where the PBL top occurs
         LTOP = 0

         DO L = 1, LLPAR
            
            ! P(L) is pressure at the top edge of (I,J,L) 
            P(L) = GET_PEDGE(I,J,L+1)

            ! LTOP is CTM layer where PBL top occurs
            IF ( BLTOP > P(L) ) THEN
               LTOP = L
               EXIT
            ENDIF
         ENDDO 

         !==============================================================
         ! IMIX(I,J) = the model level containing the PBL top
         ! FPBL(I,J) = the fraction of the IMIXth layer that the 
         !             boundary layer extends to.
         !
         ! Therefore, compute the pressure P(L), which corresponds to 
         ! the top edge grid box (I,J,L), for levels L = 1...LLPAR.  
         ! If BLTOP > P(L), then the boundary layer top occurs in grid 
         ! box (I,J,L).  (Recall that pressures decrease from the 
         ! surface upwards).  Compute IMIX and FPBL, and then define:
         !
         !    XTRA2(I,J) = FLOAT( IMIX(I,J) - 1 ) + FPBL(I,J)
         !
         ! For example, if the boundary layer at grid box (I,J) happens 
         ! to extend upward for 2.5 pressure layers, then 
         !
         !    IMIX(I,J)  = 3    (BL top is in layer 3)
         !    FPBL(I,J)  = 0.5  (BL top is 0.5 of the height of layer 3)
         !    XTRA2(I,J) = 2.5  (BL top is 2.5 layers from the surface)
         !==============================================================
         IMIX (I,J) = LTOP 
         FPBL (I,J) = 1d0 - ( BLTOP - P(LTOP) ) / ( P(LTOP-1) - P(LTOP)) 
         XTRA2(I,J) = FLOAT( IMIX(I,J) - 1 ) + FPBL(I,J)
      ENDDO
      ENDDO 
      
      !=================================================================
      ! If LTURB = T then do the BL mixing
      ! If LTURB = F then exit having computed XTRA2(I,J) 
      !=================================================================
      IF ( LTURB ) THEN

         ! Convection timestep [s]
         DTCONV = GET_TS_CONV() * 60d0

         ! We assume full mixing in the boundary layer, so the A
         ! coefficients are 1 everywhere, day & night (bmy, 2/11/03)
         A(:,:) = 1d0

         ! Loop over Lat/Long grid boxes (I,J)
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, N, AA, CC, CC_AA ) 
         DO J = 1, JJPAR
         DO I = 1, IIPAR

            ! Calculate air mass within PBL at grid box (I,J,L)
            AA = 0.d0 
            DO L = 1, IMIX(I,J)-1 
               AA = AA + AD(I,J,L) 
            ENDDO       
               
            L  = IMIX(I,J) 
            AA = AA + AD(I,J,L) * FPBL(I,J)

            ! Loop over tracers 
            DO N = 1, NTRC
     
               !========================================================
               ! Calculate tracer mass within PBL at grid box (I,J,L)
               !
               ! NOTE: CC/AA is the mean mixing ratio of tracer at 
               !       (I,J) from L=1 to L=LTOP
               !========================================================
               CC = 0.d0
               DO L = 1, IMIX(I,J)-1 
                  CC = CC + AD(I,J,L) * TC(I,J,L,N)
               ENDDO       
               
               L     = IMIX(I,J) 
               CC    = CC + AD(I,J,L) * TC(I,J,L,N) * FPBL(I,J)
 
               CC_AA = CC / AA

               !========================================================
               ! TC(I,J,L,N) new  = TC(I,J,L,N) old + 
               !                    ( DTC(I,J,L,N) / AD(I,J,L) )
               !
               ! where
               !
               ! DTC(I,J,L,N) = [ alpha * (mean MR below PBL) * 
               !                  Airmass at (I,J,L) ] -
               !                [ alpha * TC(I,J,L,N) old     * 
               !                  Airmass at (I,J,L) ]
               !
               ! DTC is thus the change in mass (kg) due to 
               ! BL mixing, so DTC/AD is the change in (V/V) 
               ! mixing ratio units.  
               !========================================================
               DO L = 1, IMIX(I,J)-1 
                  DTC(I,J,L,N) = 
     &                 ( A(I,J) * CC_AA       * AD(I,J,L) -
     &                   A(I,J) * TC(I,J,L,N) * AD(I,J,L) ) 

                  TC(I,J,L,N) = TC(I,J,L,N) + DTC(I,J,L,N)/AD(I,J,L)
               ENDDO 

               L = IMIX(I,J) 

               DTC(I,J,L,N)  = 
     &              ( A(I,J) * FPBL(I,J) * CC_AA       * AD(I,J,L) - 
     &                A(I,J) * FPBL(I,J) * TC(I,J,L,N) * AD(I,J,L) ) 
               
               TC(I,J,L,N) = TC(I,J,L,N) + DTC(I,J,L,N)/AD(I,J,L)

               !=======================================================
               ! ND15 Diagnostic: 
               ! mass change due to mixing in the boundary layer
               !=======================================================
               IF ( ND15 > 0 ) THEN
                  DO L = 1, IMIX(I,J)
                     TURBFLUP(I,J,L,N) = TURBFLUP(I,J,L,N) +
     &                    DTC(I,J,L,N) / ( TCVV(N) * DTCONV )
                  ENDDO
               ENDIF
            ENDDO    !N
         ENDDO       !I
         ENDDO       !J
!$OMP END PARALLEL DO

      ENDIF          !LTURB

      !  Return to calling program 
      END SUBROUTINE TURBDAY

!-----------------------------------------------------------------------------
!  Original code...leave here for reference (bmy, 11/10/99)
!                    TC(I,J,L,N) = 
!     &                ( A(I,J)     * AIRMAS(I,J,L) * CC/AA +
!     &                (1-A(I,J)) * TC(I,J,L,N)   * AIRMAS(I,J,L)) / 
!     &                AIRMAS(I,J,L) 
!
!                 TC(I,J,L,N) = 
!     &              ( A(I,J)        * FPBL(I,J)       * 
!     &                AIRMAS(I,J,L) * CC/AA           +
!     &               ( 1 - A(I,J)   * FPBL(I,J) )     *
!     &                TC(I,J,L,N)   * AIRMAS(I,J,L) ) / AIRMAS(I,J,L) 
!-----------------------------------------------------------------------------


