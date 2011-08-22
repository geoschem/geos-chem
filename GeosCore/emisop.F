C $Id: emisop.f,v 1.2 2009/11/30 19:57:57 ccarouge Exp $      
      FUNCTION EMISOP( I, J, IJLOOP, SUNCOS, TMMP, XNUMOL )
!
!******************************************************************************
!  Subroutine EMISOP_GRASS computes ISOPRENE EMISSIONS in units of 
!  [atoms C/box/step]. (bdf, bmy, 8/1/01, 6/16/05)
!
!  Arguments as Input:
!  ============================================================================
!  (1-2) I, J     (INTEGER ) : 2-D grid box indices
!  (3 ) IJLOOP    (INTEGER ) : 1-D grid box index
!  (4 ) SUNCOS    (REAL*8  ) : 1-D array of cos( solar zenith angle )
!  (5 ) TMMP      (REAL*8  ) : Local air temperature (K)
!  (6 ) XNUMOL    (REAL*8  ) : Number of atoms C / kg C 
!
!  Important Common Block Variables:
!  ============================================================================
!  (1 ) XYLAI     (F77_CMN_VEL ) : Leaf Area Index of land type for current MONTH
!  (2 ) IJREG     (F77_CMN_VEL ) : Number of Olson land types per grid box
!  (3 ) IJLAND+1  (F77_CMN_VEL ) : Olson land type index
!  (4 ) IJUSE     (F77_CMN_VEL ) : Olson land type fraction per box (in mils)
!  (5 ) SOPCOEFF  (F77_CMN_ISOP) : 2nd order polynomial coeffs for light correction
!  (6 ) BASEISOP  (F77_CMN_ISOP) : Baseline ISOPRENE emissions   [kg C/box/step]
!
!  NOTES:
!  (1 ) Now force double precision with DBLE and "D" exponents.  Also updated 
!        comments, made cosmetic changes (bmy, 4/4/03)
!  (2 ) Now pass I, J via the arg list.  Now reference CLDFRC directly from
!        "dao_mod.f" instead of referencing CFRAC from "F77_CMN_DEP".  Now 
!        remove reference to F77_CMN_DEP. (bmy, 12/9/03)
!  (3 ) Now scale ISOP emissions to 400 Tg C/yr for GEOS-4 (bmy, 3/5/04)
!  (4 ) Now force ISOP totals to be the same for GEOS-3 and GEOS-4 met fields
!        for the year 2001.  This will facilitate cross-model intercomparison.
!        (jal, bmy, 3/15/05)
!  (5 ) Bug fix: replace #else with #elif (swu, bmy, 6/16/05)
!  (6 ) Move XLTMMP to module MEGANUT_MOD (ccc, 11/20/09)
!******************************************************************************
!
      ! References to F90 modules
      USE DAO_MOD, ONLY     : CLDFRC
      USE MEGANUT_MOD, ONLY : XLTMMP

      USE CMN_SIZE_MOD  ! Size parameters
      USE CMN_ISOP_MOD  ! SOPCOEFF, BASEISOP
      USE CMN_VEL_MOD   ! IJREG, IJLAND, IJUSE

      IMPLICIT NONE
#     include "define.h"

      ! Arguments 
      INTEGER, INTENT(IN) :: IJLOOP,        I,    J
      REAL*8,  INTENT(IN) :: SUNCOS(MAXIJ), TMMP, XNUMOL

      ! Local variables
      INTEGER             :: INVEG
      REAL*8              :: TLAI,   EMBIO, CLIGHT

      ! External functions
!-- XLTMMP moved to meganut_mod.f (ccc, 11/20/09)
!      REAL*8,  EXTERNAL   :: XLTMMP, BIOFIT, TCORR
      REAL*8,  EXTERNAL   :: BIOFIT, TCORR

      ! Function value
      REAL*8              :: EMISOP

      !=================================================================
      ! EMISOP begins here!
      !=================================================================

      ! Initialize
      EMISOP = 0.d0
      TLAI   = 0.d0

      ! Compute total of Leaf Area Index * baseline isoprene
      ! over all Olson land types that are in this grid box      
      DO INVEG = 1, IJREG(IJLOOP)
         TLAI = TLAI + XYLAI(IJLOOP,INVEG) * BASEISOP(IJLOOP,INVEG)
      ENDDO

      !=================================================================
      ! Apply light & temperature corrections to baseline emissions --
      ! only if it is daytime and if there is nonzero isoprene emission 
      ! (e.g. XYLAI * BASEISOP > 0 )
      !=================================================================
      IF ( ( SUNCOS(IJLOOP) > 0d0 ) .AND. ( TLAI > 0d0 ) ) THEN

         ! Initialize
         EMBIO = 0d0

         ! Loop over each Olson land type in this grid box 
         DO INVEG = 1, IJREG(IJLOOP)

            ! If the product of leaf area index and baseline ISOP > 0 ...
            IF ( XYLAI(IJLOOP,INVEG) * 
     &           BASEISOP(IJLOOP,INVEG) > 0d0 ) THEN

               ! Compute light correction -- polynomial fit 
               CLIGHT = BIOFIT( SOPCOEFF,       XYLAI(IJLOOP,INVEG),
     &                          SUNCOS(IJLOOP), CLDFRC(I,J) )

               ! Apply light correction to baseline ISOPRENE emissions.
               ! Also multiply by the fraction of the grid box occupied
               ! by this Olson landtype.  Units are [kg C/box/step].
               ! BASEISOP emission rate is set in setbase.f               
               EMBIO = EMBIO + 
     &                 ( BASEISOP(IJLOOP,INVEG) * CLIGHT * 
     &                   DBLE( IJUSE(IJLOOP,INVEG) ) ) / 1000.d0
            ENDIF
         ENDDO

         ! Apply the temperature correction from Gunther et al 92 to the
         ! ISOPRENE emissions.  Units are still [kg C/box/step].
         IF ( TMMP > 273d0 ) THEN
            EMISOP = TCORR(TMMP) * EMBIO
         ELSE
            EMISOP = 0d0
         ENDIF
      ENDIF
      
      !=================================================================
      ! EMISOP is the amount of ISOP emitted in [kg/box/step]. 
      ! Convert to [atoms C/box/step] and return.
      !=================================================================
      EMISOP = EMISOP * XNUMOL

#if   defined( GEOS_3 )

      ! GEOS-3 meteorology results in 579 Tg C/yr from ISOP.  Scale 
      ! this down to 400 Tg C/yr, which is what we get from GEOS-STRAT 
      ! (mje, djj, bmy, 8/26/02)
      !
      ! NOTE: This actually produces more like 341 Tg for 2001 GEOS-3 
      !       met fields, but that is OK (jal, bmy, 3/15/05)
      EMISOP = EMISOP * ( 400d0 / 579d0 )

#elif defined( GEOS_4 )

      ! Original GEOS-4 scaling produced 443 Tg C/yr w/ 2003 "V3" met 
      ! fields.  However we have since switched to GEOS-4 "V4" met fields 
      ! and need to rescale the ISOP total.  A recent run with GEOS-4 "V4"
      ! met fields for 2001 produced 443 Tg C/yr.  We need to force the
      ! total to be the same as for GEOS-3, for comparison purposes.
      ! Therefore apply a second scale factor so that we get 341 Tg C/yr
      ! of ISOP for GEOS-4 "V4" met fields for 2001. (bmy, 3/15/05)
      EMISOP = EMISOP * ( 400d0      / 443d0      )
      EMISOP = EMISOP * ( 341.2376d0 / 442.7354d0 )

#endif 

      ! Return to calling program
      END FUNCTION EMISOP
