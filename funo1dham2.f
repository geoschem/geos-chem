! $Id: funo1dham2.f,v 1.1 2003/06/30 20:26:04 bmy Exp $
      REAL*8 FUNCTION FUNO1DHAM2( WAVE, TT )
!
!******************************************************************************
!  Function FUNO1DHAM2 returns O(1D) yield from O3 photolysis for the SLOW-J
!  photolysis mechanism.  Parameters are taken from model of Michelsen et al.
!  (lwh, jyl, gmg, djj, 1994; bmy, 4/4/03)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) WAVE (REAL*8) : Wavelength [Angstroms???]
!  (2 ) TT   (REAL*8) : Temperature [K]
!
!  NOTES:
!  (1 ) Updated comments, cosmetic changes (bmy, 4/4/03)
!******************************************************************************
!
      IMPLICIT NONE 

      ! Arguments
      REAL*8, INTENT(IN) :: WAVE,TT

      ! Local variables
      REAL*8  :: WV_HAM(21)
      DATA       WV_HAM/
     &     3050.,3060.,3070.,3080.,3090.,3100.,3110.,3120.,3130.,3140.,
     &     3150.,3160.,3170.,3180.,3190.,3200.,3210.,3220.,3230.,3240.,
     &     3250./

      REAL*8  :: A_HAM(21)
      DATA       A_HAM/
     &     1.01,1.01,1.05,1.15,1.39,1.90,2.93,
     &     4.87,8.21,13.3,17.6,20.4,18.0,21.8,
     &     18.1,17.2,7.99,11.5,14.3,10.7,11.8/

      REAL*8  :: C_HAM(21)
      DATA       C_HAM/
     &     3.933,11.51,33.09,79.39,159.9,272.5,407.9,
     &     551.4,682.3,791.6,851.3,903.8,900.3,948.4,
     &     891.1,1066.,969.4,1201.,1182.,1152.,1435./
      
      REAL*8  :: XK_HAM
      DATA       XK_HAM/.69503/

      INTEGER :: NWV_HAM, ILW_HAM, IUP_HAM
      DATA       NWV_HAM/21/
      

      REAL*8  :: PHILW, PHIUP, DELWAVE, FS1, FS2

      !=================================================================
      ! FUNO1DHAM2 begins here!
      !=================================================================
      FUNO1DHAM2 = 0.0

      IF( WAVE.LE.2713.) THEN
         FUNO1DHAM2 = 0.87

      ELSE IF ( WAVE.GT.2713 .AND. WAVE.LE.3040. ) THEN
         FUNO1DHAM2 = 1.98 - 3010./WAVE

      ELSE IF ( WAVE.GT.3040 .AND. WAVE.LE.3250. ) THEN
         CALL HUNT( WV_HAM, NWV_HAM, WAVE, ILW_HAM )

         IF ( ILW_HAM.GT.0 .AND. ILW_HAM.LT.NWV_HAM ) THEN
            IUP_HAM    = ILW_HAM+1

            ! Lower Wavelength
            PHILW      = A_HAM(ILW_HAM)*
     &                   EXP(-C_HAM(ILW_HAM)/(XK_HAM*TT))

            ! Upper Wavelength
            PHIUP      = A_HAM(IUP_HAM)*
     &                   EXP(-C_HAM(IUP_HAM)/(XK_HAM*TT))

            DELWAVE    = WV_HAM(IUP_HAM)-WV_HAM(ILW_HAM)
            FS1        = (WV_HAM(IUP_HAM)-WAVE)/DELWAVE
            FS2        = (WAVE-WV_HAM(ILW_HAM))/DELWAVE
            FUNO1DHAM2 = FS1*PHILW+fs2*PHIUP

         ELSE IF ( ILW_HAM.EQ.NWV_HAM ) THEN
            FUNO1DHAM2 = A_HAM(ILW_HAM)*
     +                   EXP(-C_HAM(ILW_HAM)/(XK_HAM*TT))
         ENDIF
      ELSE
         FUNO1DHAM2 = 0.0         ! WAVE > 3250
      ENDIF
      
      ! Return to calling program
      END FUNCTION FUNO1DHAM2
