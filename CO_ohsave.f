C $Id: CO_ohsave.f,v 1.1 2003/06/30 20:26:07 bmy Exp $
C*****************************************************************************
      SUBROUTINE CO_OHSAVE(BOH) 
C*****************************************************************************
C Archive for CH3CCl3 lifetime.
C Modified by bnd from SR OHSAVE.
C*****************************************************************************
C The annual mean tropopause is stored in the LPAUSE array 
C (from header file "CMN").  LPAUSE is defined such that: 
C
C Levels            1 <= L <= LPAUSE(I,J) - 1 are tropospheric
C         LPAUSE(I,J) <= L <= LLPAR           are stratospheric
C
C We now use LPAUSE instead of NSKIPL to denote the strat/trop boundary. 
C (bmy, 4/18/00)  
C*****************************************************************************
C
      USE DIAG_MOD, ONLY : DIAGCHLORO
C
      IMPLICIT NONE

#     include "CMN_SIZE"  ! Size parameters
#     include "CMN"       ! LPAUSE
#     include "CMN_DIAG"  ! ND23 switch
#     include "CMN_CO"    ! BAIRDENS

      INTEGER   I,J,L
      REAL*8    KCLO,LOSS,OHMASS,MASST
      REAL*8    BOH(IIPAR,JJPAR,LLPAR)

      REAL*8    BOXVL
      EXTERNAL  BOXVL

      print*,'Save OH: methylchloroform lifetime in SR CO_ohsave.f.' 

!*****************************************************************************
      IF ( ND23 > 0 ) THEN
!*****************************************************************************
         DO L = 1, MAXVAL( LPAUSE )
         DO J = 1, JJPAR 
         DO I = 1, IIPAR 
            
            ! Only process tropospheric boxes (bmy, 4/17/00)
            IF ( L < LPAUSE(I,J) ) THEN
               KCLO   = 1.8D-12 * EXP( -1550.D0 / Tavg(I,J,L) )
               
               LOSS   = KCLO            * BOH(I,J,L)  *
     &                  bairdens(I,J,L) * BOXVL(I,J,L)
                
               OHMASS = BOH(I,J,L) * BAIRDENS(I,J,L) * BOXVL(I,J,L)
               
               MASST  = BAIRDENS(I,J,L) * BOXVL(I,J,L)

               ! Store loss in DIAGCHLORO(I,J,L,1)
               DIAGCHLORO(I,J,L,1) = DIAGCHLORO(I,J,L,1) + LOSS
               
               ! Store OH mass in DIAGCHLORO(I,J,L,2)
               DIAGCHLORO(I,J,L,2) = DIAGCHLORO(I,J,L,2) + OHMASS

               ! Store total mass in DIAGCHLORO(I,J,L,3)
               DIAGCHLORO(I,J,L,3) = DIAGCHLORO(I,J,L,3) + MASST
            ENDIF
         ENDDO
         ENDDO
         ENDDO
!*****************************************************************************
      ENDIF
!*****************************************************************************

      RETURN
      END SUBROUTINE CO_OHSAVE
