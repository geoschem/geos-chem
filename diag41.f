! $Id: diag41.f,v 1.1 2003/06/30 20:26:04 bmy Exp $
      SUBROUTINE DIAG41 
!
!*****************************************************************************
!  Subroutine DIAG41 produces monthly mean boundary layer height in meters 
!  between 1200-1600 local time for the U.S. geographical domain. 
!  (amf, bmy, 11/18/99, 2/11/03)
!
!  DIAG41 writes timseries output to disk using the Harvard CTM binary
!  punch file version 2.0.  This format is very GAMAP-friendly!
!
!  Input via "CMN" header file:
!  ===========================================================================
!  (1 ) XTRA2 : Height of PBL in boxes
!
!  NOTES:
!  (1 ) DIAG41 is written in Fixed-Format F90. 
!  (2 ) XTRA2 must be computed by turning TURBDAY on first.  Also,
!        XTRA2 is a global-size array, so use window offsets IREF, JREF
!        to index it correctly. (bmy, 11/18/99)
!  (3 ) Do a little rewriting so that the DO-loops get executed
!        in the correct order (J first, then I). (bmy, 11/18/99)
!  (4 ) AD41 is now declared allocatable in "diag_mod.f". (bmy, 12/6/99)
!  (5 ) AFTTOT is now declared allocatable in "diag_mod.f". (bmy, 3/17/00)
!  (6 ) Remove NYMD from the argument list -- it wasn't used (bmy, 6/22/00) 
!  (7 ) XTRA2(IREF,JREF,5) is now XTRA2(I,J).  Also updated comments. 
!        (bmy, 9/25/01)
!  (8 ) Removed obsolete code from 9/01 (bmy, 10/23/01)
!  (9 ) Now reference BXHEIGHT from "dao_mod.f".  Also removed obsolete
!        code. (bmy, 9/18/02)
!  (10) Now use function GET_LOCALTIME from "dao_mod.f" (bmy, 2/11/03)
!*****************************************************************************
!
      ! References to F90 modules
      USE DAO_MOD,  ONLY : BXHEIGHT
      USE DIAG_MOD, ONLY : AD41, AFTTOT
      USE TIME_MOD, ONLY : GET_LOCALTIME

      IMPLICIT NONE

#     include "CMN_SIZE"   ! IIPAR, DISIZE
#     include "CMN_DIAG"   ! AIJ, AFTTOT
#     include "CMN"        ! TOFDAY

      ! Local variables
      INTEGER             :: GOOD(IIPAR)
      INTEGER             :: I, IREF, J, JREF, L, M, PBLINT

      REAL*8              :: XLOCTM(IIPAR)
      REAL*8              :: TEMPBL, PBLDEC

      !=================================================================
      ! DIAG41 begins here!
      !=================================================================

      ! XLOCTM: Local time at each longitude 
      DO I = 1, JJPAR
         XLOCTM(I) = GET_LOCALTIME( I ) 
      ENDDO

      ! GOOD=1 denotes longitudes where it is between 1200-1600 LT
      GOOD = 0
      WHERE ( XLOCTM >= 12d0 .and. XLOCTM <= 16d0 ) GOOD = 1

      !=================================================================
      ! For grid boxes where it is between 1200 and 1600 local time, 
      ! do the following:
      !
      ! (1) Set TEMPBL = 0.0.  TEMPBL will be used to store the boundary
      !     layer height in meters.
      !
      ! (2) Get the integer (PBLINT) and decimal (PBLDEC) parts of PBL.  
      !     Recall that XTRA2 is the boundary layer height in number of 
      !     boxes, so if XTRA2 = 2.5 (for example), then PBLINT = 2 and 
      !     PBLDEC = 0.5
      !
      ! (3) If PBLINT > 0 then the boundary layer height for this grid 
      !     box is nonzero.  Add to TEMPBL all of the box heights for 
      !     each of the "full" grid boxes, and then add on the fraction 
      !     of the highest box. 
      !
      ! (4) Store TEMPBL for this (I,J) location to AD41(I,J,1) array.  
      !     Also increment AFTTOT(I,J), which keeps track of the number 
      !     of entries in AD41(I,J,1).
      !
      ! (5) Store the afternoon B-L height in # of boses in AD41(I,J,2).
      !     (bmy, 12/6/99)
      !  
      !  NOTES:
      !  (1) Multiplying BXHEIGHT by GOOD will only cause the boxes 
      !      between 1200-1600 LT to be counted.  Recall that GOOD = 1 
      !      for boxes between 1200-1600 LT, and GOOD = 0 otherwise.  
      !      This should execute faster than performing IF tests on 
      !      every grid box. (bmy, 11/18/99)
      !=================================================================
      DO J = 1, JJPAR
      DO I = 1, IIPAR
         TEMPBL = 0d0

         ! Integer and fractional parts of XTRA2
         PBLINT = FLOOR( XTRA2(I,J) )
         PBLDEC = XTRA2(I,J) - PBLINT

         ! "Full" grid boxes
         IF ( PBLINT > 0 ) THEN
            DO M = 1, PBLINT 
               TEMPBL = TEMPBL + ( BXHEIGHT(I,J,M) * GOOD(I) )
            ENDDO
         ENDIF

            ! "Fraction" of the highest grid box
         TEMPBL = TEMPBL + 
     &            ( PBLDEC * BXHEIGHT(I,J,PBLINT+1) * GOOD(I) )

         ! Store to diagnostic arrays 
         AD41(I,J,1) = AD41(I,J,1) + TEMPBL
         AD41(I,J,2) = AD41(I,J,2) + ( XTRA2(I,J) * GOOD(I) )

         ! Increment counter of afternoon boxes
         AFTTOT(I,J) = AFTTOT(I,J) + GOOD(I)

         ENDDO     ! end I loop
      ENDDO        ! end J loop

      ! Return to calling program
      END SUBROUTINE DIAG41
