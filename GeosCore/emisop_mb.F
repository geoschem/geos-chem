! $Id: emisop_mb.f,v 1.1 2009/09/16 14:06:32 bmy Exp $      
      FUNCTION EMISOP_MB( I, J, IJLOOP, SUNCOS, TMMP, XNUMOL )
!
!******************************************************************************
!  Subroutine EMISOP_MB computes METHYL BUTENOL emissions in units
!  of [atoms C/box/step]. (bdf, bmy, 8/2/01, 6/16/05)
!
!  Arguments as Input:
!  ============================================================================
!  (1-2) I, J     (INTEGER ) : 2-D grid box indices
!  (1 ) IJLOOP (INTEGER) : 1-D grid box index
!  (2 ) SUNCOS (REAL*8 ) : 1-D array of cos( solar zenith angle )
!  (3 ) TMMP   (REAL*8 ) : Local air temperature (K)
!  (4 ) XNUMOL (REAL*8 ) : Number of atoms C / kg C 
!
!  Important Common Block Variables:
!  ============================================================================
!  (1 ) XYLAI     (F77_CMN_VEL ) : Leaf Area Index of land type for current MONTH
!  (2 ) IJREG     (F77_CMN_VEL ) : Number of Olson land types per grid box
!  (3 ) IJLAND+1  (F77_CMN_VEL ) : Olson land type index
!  (4 ) IJUSE     (F77_CMN_VEL ) : Olson land type fraction per box (in mils)
!  (5 ) SOPCOEFF  (F77_CMN_ISOP) : 2nd order polynomial coeffs for light correction
!  (6 ) BASEISOP  (F77_CMN_ISOP) : Baseline ISOPRENE emissions    [kg C/box/step]
!  (7 ) BASEMB    (F77_CMN_ISOP) : Baseline METHYL BUT. emissions [kg C/box/step]
!
!  NOTES:
!  (1 ) Now use F90 syntax.  Use "D" exponents to force double precision.
!        Updated comments, and mad cosmetic changes (bmy, 8/2/01) 
!  (2 ) Deleted obsolete, commented-out code from 8/01 (bmy, 11/27/01)
!  (3 ) GEOS-3 meteorology results in 579 Tg C/yr from biogenic ISOP.  Compute
!        ISOP from grasslands based on 400 Tg C/yr from biogenic ISOP, which 
!        is what we get from GEOS-STRAT. (mje, bdf, djj, 9/10/02)
!  (4 ) Now pass I, J via the arg list.  Now reference CLDFRC directly from
!        "dao_mod.f" instead of referencing CFRAC from "F77_CMN_DEP".  Now 
!        remove reference to F77_CMN_DEP. (bmy, 12/9/03)
!  (5 ) Now scale ISOP emissions to 400 Tg C/yr for GEOS-4 (bmy, 3/5/04)
!  (6 ) Now force ISOP totals to be the same for GEOS-3 and GEOS-4 met fields
!        for the year 2001.  This will facilitate cross-model intercomparison.
!        (jal, bmy, 3/15/05)
!  (7 ) Bug fix: change #else to #elif (swu, bmy, 6/16/05)
!******************************************************************************
!
      ! References to F90 modules
      USE DAO_MOD, ONLY : CLDFRC

      USE CMN_SIZE_MOD  ! Size parameters
      USE CMN_ISOP_MOD  ! SOPCOEFF, BASEISOP, BASEMB
      USE CMN_VEL_MOD   ! IJREG, IJLAND, IJUSE

      IMPLICIT NONE
#     include "define.h"

      ! Arguments
      INTEGER, INTENT(IN) :: IJLOOP,        I,    J
      REAL*8,  INTENT(IN) :: SUNCOS(MAXIJ), TMMP, XNUMOL

      ! Local variables   
      INTEGER             :: INVEG, MBO_SCALE, TEST
      REAL*8              :: EMBIO, TLAI, CLIGHT, EMISOP_MB

      ! External functions
      REAL*8,  EXTERNAL   :: BIOFIT, TCORR

      !=================================================================
      ! EMISOP_MB begins here!
      !=================================================================

      ! Initialize
      EMISOP_MB = 0d0
      TLAI      = 0d0

      ! Compute total of Leaf Area Index * baseline isoprene
      ! over all Olson land types that are in this grid box      
      DO INVEG = 1,IJREG(IJLOOP)
         TLAI = TLAI + XYLAI(IJLOOP,INVEG) * BASEISOP(IJLOOP,INVEG)
      END DO

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

            ! IJLAND+1 is the Olson land type index
            ! For methyl butenol emissions the landtypes 21, 22, 23, 
            ! and 28 are mostly pine forests and emit MB.  Landtypes 
            ! 24 and 25 are half pine, and emit MB at half the rate.  
            ! Other landtypes emit no MB.
            SELECT CASE ( IJLAND(IJLOOP,INVEG) + 1 )
               CASE ( 21, 22, 23, 28 )
                  MBO_SCALE = 2
               CASE ( 24, 25 )
                  MBO_SCALE = 1
               CASE DEFAULT
                  MBO_SCALE = 0
            END SELECT

            ! If the product of leaf area index and baseline ISOP > 0 ...
            IF ( XYLAI(IJLOOP,INVEG) * 
     &           BASEISOP(IJLOOP,INVEG) > 0.0 ) THEN

               ! Compute light correction -- polynomial fit
               CLIGHT = BIOFIT( SOPCOEFF,       XYLAI(IJLOOP,INVEG),
     &                          SUNCOS(IJLOOP), CLDFRC(I,J) )

               ! Apply light correction to baseline MB emissions.
               ! Also multiply by the fraction of the grid box occupied
               ! by this Olson landtype.  Units are [kg C/box/step].
               ! BASEMB (set in setbase.f) is computed to get Guenther's 
               ! North American emissions of 3.2 Tg C/yr from MB.
               EMBIO = EMBIO + 
     &                 ( BASEMB(IJLOOP) * MBO_SCALE * CLIGHT *
     &                   DBLE( IJUSE(IJLOOP,INVEG) ) / 1000d0 )
            ENDIF
         ENDDO

         ! Apply the temperature correction from Gunther et al 92 to the
         ! METHYL BUTENOL emissions.  Units are still [kg C/box/step].
         IF ( TMMP > 273d0 ) THEN
            EMISOP_MB = TCORR(TMMP) * EMBIO
         ELSE
            EMISOP_MB = 0d0
         ENDIF

      ENDIF

      !=================================================================
      ! EMISOP_MB is the amount of METHYL BUTENOL emitted in 
      ! [kg/box/step]. Convert to [atoms C/box/step] and return.
      !=================================================================
      EMISOP_MB = EMISOP_MB * XNUMOL

#if   defined( GEOS_3 )

      ! GEOS-3 meteorology results in 579 Tg C/yr from ISOP.  Scale 
      ! this down to 400 Tg C/yr, which is what we get from GEOS-STRAT 
      ! (mje, djj, bmy, 8/26/02)  
      !
      ! NOTE: This actually produces more like 341 Tg for 2001 GEOS-3 
      !       met fields, but that is OK (jal, bmy, 3/15/05)
      EMISOP_MB = EMISOP_MB * ( 400d0 / 579d0 )

#elif defined( GEOS_4 )

      ! Original GEOS-4 scaling produced 443 Tg C/yr w/ 2003 "V3" met 
      ! fields.  However we have since switched to GEOS-4 "V4" met fields 
      ! and need to rescale the ISOP total.  A recent run with GEOS-4 "V4"
      ! met fields for 2001 produced 443 Tg C/yr.  We need to force the
      ! total to be the same as for GEOS-3, for comparison purposes.
      ! Therefore apply a second scale factor so that we get 341 Tg C/yr
      ! of ISOP for GEOS-4 "V4" met fields for 2001. (bmy, 3/15/05)
      EMISOP_MB = EMISOP_MB * ( 400d0      / 443d0      )
      EMISOP_MB = EMISOP_MB * ( 341.2376d0 / 442.7354d0 )

#endif 

      ! Return to calling program
      END FUNCTION EMISOP_MB
