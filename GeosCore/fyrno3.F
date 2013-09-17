! $Id: fyrno3.f,v 1.1 2009/09/16 14:06:27 bmy Exp $
      REAL*8 FUNCTION FYRNO3( XCARBN, ZDNUM, TT )
!
!******************************************************************************
!  Function FYRNO3 returns organic nitrate yields YN = RKA/(RKA+RKB)
!  from RO2+NO reactions as a function of the number N of carbon atoms.
!  (lwh, jyl, gmg, djj, bmy, 1/1/89, 6/26/03)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) XCARBN (REAL*8) : Number of C atoms in RO2
!  (2 ) ZDNUM  (REAL*8) : Air density   [molec/cm3 ]
!  (3 ) TT     (REAL*8) : Temperature   [K         ]
!    
!  NOTES: 
!  (1 ) Original code from Larry Horowitz, Jinyou Liang, Gerry Gardner,
!        and Daniel Jacob circa 1989/1990.
!  (2 ) Updated following Atkinson 1990.
!  (3 ) Change yield from Isoprene Nitrate (ISN2) from 0.44% to 12%,
!        according to Sprengnether et al., 2002. (amf, bmy, 1/7/02)
!  (4 ) Eliminate obsolete code from 1/02 (bmy, 2/27/02)
!  (5 ) Updated comment description of XCARBN (bmy, 6/26/03)
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
!******************************************************************************
!
      IMPLICIT NONE

      ! Arguments
      REAL*8, INTENT(IN) :: XCARBN, ZDNUM, TT

      ! Local variables
      REAL*8             :: YYYN, XXYN,  AAA,  RARB, ZZYN
      REAL*8             :: XF,   ALPHA, Y300, BETA, XMINF, XM0

      ! Initialize static variables
      DATA   Y300,ALPHA,BETA,XM0,XMINF,XF/.826,1.94E-22,.97,0.,8.1,.411/

      !=================================================================
      ! FYRNO3 begins here!
      !=================================================================
      XXYN   = ALPHA*EXP(BETA*XCARBN)*ZDNUM*((300./TT)**XM0)
      YYYN   = Y300*((300./TT)**XMINF)
      AAA    = LOG10(XXYN/YYYN)
      ZZYN   = 1./(1.+ AAA*AAA )
      RARB   = (XXYN/(1.+ (XXYN/YYYN)))*(XF**ZZYN)
      FYRNO3 = RARB/(1. + RARB)
     
      ! Return to calling program
      END FUNCTION FYRNO3
