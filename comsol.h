! $Id: comsol.h,v 1.1 2003/06/30 20:26:03 bmy Exp $
!
!**********************************************************************
!                                                                     *
!  HARVARD TROPOSPHERIC CHEMISTRY MODULE FOR 3-D APPLICATIONS         *
!  by Larry Horowitz, Jinyou Liang, Gerry Gardner, Prof. Daniel Jacob *
!  of HARVARD UNIVERSITY    (Release V1.0)                            *
!**********************************************************************
!                                                                     *
! ***** for JVALUE calculations.                                      *
!
! NOTES:
! (1 ) Changed RCS ID tag comment character  from "C" to "!" to allow 
!       freeform compilation.  Also added & continuation characters in 
!       column 73 to allow header files to be included in F90 freeform 
!       files.  Also changed comment character from "C" to "!" to avoid
!       problems when inlining into freeform files. (bmy, 6/25/02)
! (2 ) Renamed cpp switch from DEC_COMPAQ to COMPAQ.  Now declare CSECT 
!       common block as THREADPRIVATE for all platforms. (bmy, 3/23/03)
!**********************************************************************
      INTEGER MAXTEMP,MXSPE,MXBRCH,MXWL,MAXPTS2,NSTDL
      PARAMETER (MAXTEMP=5,MXSPE=60,MXBRCH=5,MXWL=70)

      ! Now use LTROP (from CMN_SIZE to set MAXPTS2 (bmy, 4/19/99)
      PARAMETER (MAXPTS2=IIPAR*JJPAR*LLTROP, NSTDL = 41)

      INTEGER INAME
      COMMON /ICSETAL/ INAME(MXSPE)

      CHARACTER*8 CNAME
      COMMON /CCSETAL/ CNAME(MXSPE)

      REAL*8 XSECT, TARRAY, WL, ACTFLX, FL,                             &
     &      STDAIR,  STDO3, TATM
      COMMON /CSETAL/                                                   &
     &                TARRAY(MXSPE,MXBRCH,MAXTEMP),WL(MXWL),            & 
     &                ACTFLX(MXWL+1,MAXPTS2),FL(MXWL),                  &
     &                STDAIR(NSTDL,35),STDO3(NSTDL,10),TATM(NSTDL,35)

      COMMON /CSECT/ XSECT(MXWL,MXSPE,MAXTEMP,MXBRCH)

      !=================================================================
      ! /CSECT/ needs to be declared THREADPRIVATE for the OpenMP
      ! parallelizaton for all platforms (bmy, 3/23/03)
      !=================================================================
!$OMP THREADPRIVATE( /CSECT/ )

      REAL*8 O3DU
      COMMON /CSETA2/ O3DU(11,12)

      REAL*8 AERSOL,AERXCT
      COMMON /CSETA1/ AERSOL(6),AERXCT(NSTDL)

