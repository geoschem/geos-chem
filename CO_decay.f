C $Id: CO_decay.f,v 1.1 2003/06/30 20:26:09 bmy Exp $
      SUBROUTINE CO_decay(BOH)
C*****************************************************************************
      
      ! References to F90 modules
      USE DAO_MOD, ONLY : AIRVOL

      IMPLICIT NONE

#     include "CMN_SIZE"       ! Size parameters
#     include "CMN"            ! STT, LPAUSE
#     include "CMN_CO"         ! BAIRDENS, CO arrays
#     include "CMN_CO_BUDGET"  ! TCO
#     include "CMN_OH"         ! Arrays for OH parameterization

      REAL*8    BOXVL
      EXTERNAL  BOXVL

      INTEGER I,J,L,M,N
      REAL*8 DT
      REAL*8 BOH(IIPAR,JJPAR,LLPAR)

      ! Variable to hold the maximum value of LPAUSE (bmy, 4/18/00)
      INTEGER LPAUSE_MAX

      print*,'Calculating CO tropospheric decay in SR CO_decay.'
C
C*****************************************************************************
C This SR calculates the decay rate of CO by OH.  OH is the only sink
C    for CO considered here.
C Written by bnd.
C
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
C Monthly loss of CO is summed in TCO(3)
C    TCO(3)  = CO sink by OH
C
C DECAY RATE
C The decay rate is calculated by:
C
C    OH + CO -> products (JPL '97)
C    k = (1 + 0.6Patm) * 1.5E-13
C
C convert the chemistry time step (NCHEM) from minutes to seconds
C   and from integer to real.
C
      DT=REAL(NCHEM)*60. 
C
C*****************************************************************************
C get pressure (atm) for SR CO_decay's calculation of the
C     rate constant for CO+OH->
C*****************************************************************************
C
      ! Store the maximum value of LPAUSE in LPAUSE_MAX 
      ! for use in the DO-loops below (bmy, 4/18/00)
      LPAUSE_MAX = MAXVAL( LPAUSE )
      
      ! The L-loop only needs to go up to LPAUSE_MAX (bmy, 4/18/00)
      DO L = 1, LPAUSE_MAX
      DO J = 1, JJPAR
      DO I = 1, IIPAR
C
         IF ( L < LPAUSE(I,J) ) THEN 
C     
            Pco(I,J,L)= 0.
            Pco(I,J,L)=Pavg(I,J,L)/1013.25
C     
            krate(I,J,L)=(1.+0.6*Pco(I,J,L))*1.5E-13
C
C convert STT(kgs/box) to STT(molec/cc) (see SR convert_units.f)
C   kg CO/box * box/cc * mole/0.028 kg CO * Avog.#/mole
C
            STTTOGCO(I,J,L)=1./AIRVOL(I,J,L)/1.D6/FMOL_CO*6.023D23
C
            GCO(I,J,L)=STT(I,J,L,1)*STTTOGCO(I,J,L)
C
         ENDIF
C
      ENDDO
      ENDDO
      ENDDO
C
C sum loss in TCO(3) (molecules/box)
C
      ! The L-loop only needs to go up to LPAUSE_MAX (bmy, 4/18/00)
      DO L = 1, LPAUSE_MAX
      DO J = 1, JJPAR
      DO I = 1, IIPAR
C
         IF ( L < LPAUSE(I,J) ) THEN
C
            TCO(I,J,L,3)=TCO(I,J,L,3)+GCO(I,J,L)*BOXVL(I,J,L)*
     *           krate(I,J,L)*BOH(I,J,L)*DT
C     
         ENDIF
C
      ENDDO
      ENDDO
      ENDDO
C
C calculate new CO value : [CO]=[CO](1-k[OH]*delt)
C      
      ! The L-loop only needs to go up to LPAUSE_MAX (bmy, 4/18/00)
      DO L = 1, LPAUSE_MAX
      DO J = 1, JJPAR
      DO I = 1, IIPAR
C
         IF ( L < LPAUSE(I,J) ) THEN 
C
            GCO(I,J,L)=GCO(I,J,L)*(1.-krate(I,J,L)*BOH(I,J,L)*DT)
C 
C convert new GCO and copy value into STT 
C
            STT(I,J,L,1)=GCO(I,J,L)/STTTOGCO(I,J,L)
C
         ENDIF
C
      ENDDO
      ENDDO
      ENDDO
C
      RETURN
      END
