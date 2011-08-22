!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
! !MODULE: smv_physconst_mod.F
!
! !DESCRIPTION:
!
! !REMARKS:
! (1 ) Changed RCS ID tag comment character from "C" to "!" to
! allow freeform compilation.  Also added the & continuation 
! characters in column 73 to allow header files to be 
! included in F90 freeform files. (bmy, 6/25/02)
! (2 ) Commented out P,  we now use routines from "pressure_mod.f"
! to compute pressure for hybrid grid.  Removed TCRIJ, 
! TCSCAL, TOTALC, TCWETL, SIG, SIGE, AREAW, AREAG, AREAH,
! DXP, DYP, DXYP, DXYV, GRAV, LSRC1, LSRC2, LSRC3, LSRC4,
! NSTRTC, LDIFF1, LDIFFM, KDIFFU, NINST, NINIT, ICASE, LTM
! TLAT, TLNG -- these are now obsolete. 
! (dsa, bdf, bmy, 8/26/02)
! (3 ) Moved AD and T to "dao_mod.f".  The fact that these
! were in "CMN" is historical baggage. (bmy, 9/18/02)
! (4 ) Removed RH2, STH, RLAT, RLATV, TAUI, TAUE, JMONTH, 
! JMONTH0, XMID, YLMID, YEDGE, DXYP, TAU, TAU0, TOFDAY,
! NTAU, IDAY, JDATE, JYEAR, JDAY, MONTH, IYEAR, I00, 
! J00, JMSIZE, NREAD, NWRITE, NDYN, NCONV, NDIAG, NCHEM, 
! NSRCE.  Most of these are now superseded by either
! "grid_mod.f" or "time_mod.f". (bmy, 3/11/03)
! (5 ) Fix for LINUX - remove & in column 73 (bmy, 6/27/03)
! (6 ) Remove obsolete variables and XTRA2 (bmy, 8/2/05)
! (7 ) Remove obsolete embedded chemistry stuff (bmy, 2/25/10)
!================================================================= 
!\\
!\\
! !INTERFACE: 
!
      MODULE CMN_MOD
!
! !USES:
!
      USE CMN_SIZE_MOD, ONLY : IIPAR, JJPAR

      IMPLICIT NONE
!      PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS: KITCHEN SINK!
!
      PUBLIC
	SAVE
! !DEFINED PARAMETERS: 
!

!------------------------------------------------------------------------------
! Prior to 2/25/10:
! Removed the obsolete embedded chemistry stuff (bmy, 2/25/10)
!      !=================================================================
!      ! IEBD1, JEBD1, IEBD2, JEBD2 define the embedded chemistry region
!      !=================================================================
!      INTEGER ::      IEBD1, IEBD2, JEBD1, JEBD2
!      COMMON /EMBCHM/ IEBD1, IEBD2, JEBD1, JEBD2
!------------------------------------------------------------------------------

      !=================================================================
      ! LPAUSE  stores the annual mean tropopause (bmy, 12/6/99)
      ! IFLX    stores the flags for ND27 diagnostic (bmy, 12/6/99)
      !=================================================================
      INTEGER ::       LPAUSE(IIPAR,JJPAR), IFLX(IIPAR,JJPAR)

	END MODULE CMN_MOD

