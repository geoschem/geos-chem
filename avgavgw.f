C $Id: avgavgw.f,v 1.1 2003/06/30 20:26:01 bmy Exp $
      !SUBROUTINE avgavgw(NTDT,AVGW,NMIN,NCHEM,SPHU)
      SUBROUTINE avgavgw( NTDT, NMIN, NCHEM )
C
C Designed to calculate a 24-hour average water vapor
C   to be used for the CO/OH parameterization
C   option (see SR chemco).
C Created by Bryan Duncan.
C
      ! References to F90 modules
      USE DAO_MOD, ONLY : MAKE_AVGW, AVGW

      IMPLICIT NONE
 
#     include "CMN_SIZE"
#     include "CMN_CO"
      
      INTEGER NTDT,NTIMES,MNDT,NMIN,NCHEM

      !-------------------------------------------------------------
      ! NOTE: We no longer need to pass AVGW and SPHU to MAKE_AVGW
      !       since these are fields in "dao_mod.f" (bmy, 11/15/02)
      !! AVGW is now (IIPAR,JJPAR,LLPAR) (bmy, 9/24/01)
      !REAL*8 AVGW(IIPAR,JJPAR,LLPAR), SPHU(IIPAR,JJPAR,LLPAR)
      !-------------------------------------------------------------

      MNDT=NTDT/60
      NTIMES=NCHEM/MNDT

      IF(NMIN.LE.NCHEM) NTIMES=NTIMES+1

      IF (NMIN.EQ.0) THEN
         LNEW=0
         LNCOUNT=0
      ENDIF
      
      IF (LNEW.EQ.0) THEN
         Wavg(:,:,:)=0.
         LNEW=1
         LNCOUNT=0
      ENDIF
C
C Calculate the water vapor, AVGW.
C
      CALL MAKE_AVGW
C
C Sum the water vapors in Wavg.
C
      Wavg(:,:,:)=Wavg(:,:,:)+AVGW(:,:,:)
C
C Keep track to see if at end of NCHEM time step.
C
      LNCOUNT=LNCOUNT+1
C
      IF (LNCOUNT.EQ.NTIMES) THEN
C
C Calculate the 24-hour average water vapor mixing ratio (ppmv).
C     
         Wavg(:,:,:)=Wavg(:,:,:)/REAL(NTIMES)*1.E6
C
         LNEW=0
C
      ENDIF
C     
      RETURN
      END
      
