C $Id: CO_strat.f,v 1.1 2003/06/30 20:26:08 bmy Exp $
      SUBROUTINE CO_STRAT()
C*****************************************************************************
      ! References to F90 modules
      USE DAO_MOD, ONLY : AIRVOL

      IMPLICIT NONE

#     include "CMN_SIZE"       ! Size parameters
#     include "CMN"            ! STT, LPAUSE
#     include "CMN_CO"         ! BAIRDENS, STT2GCO, CO arrays
#     include "CMN_CO_BUDGET"  ! TCO, CO budget arrays

      REAL*8   BOXVL
      EXTERNAL BOXVL

      INTEGER  I,J,L,M,N
      REAL*8   DT,LLVERT

      ! Variable to hold the minimum value of LPAUSE (bmy, 4/18/00)
      INTEGER  LPAUSE_MIN
C
C*****************************************************************************
C This SR calculates uses production and loss rates for CO to calculate
C  net production of CO in the stratosphere.  The purpose of this SR is to
C  prevent high CO concentrations from building up in the stratosphere; in
C  these layers only transport is simulated (i.e., no chemistry).  For
C  a long simulation, a buildup of high concentrations could occur causing
C  the stratosphere to become a significant source of CO.
C Written by Bryan Duncan.
C
C Production (mixing ratio/sec) and loss (1/sec) rates provided by Dylan
C  Jones.  Only production by CH4+OH and destruction by CO+OH are considered.
C
C*****************************************************************************
C The annual mean tropopause is stored in the LPAUSE array 
C (from header file "CMN").  LPAUSE is defined such that: 
C
C Levels            1 <= L <= LPAUSE(I,J) - 1 are tropospheric
C         LPAUSE(I,J) <= L <= LLPAR           are stratospheric
C
C (bmy, 4/18/00)  
C*****************************************************************************
C
      print*,
     &'Calculating stratospheric CO loss/decay in SR CO_strat.'
C
C convert the chemistry time step (NCHEM) from minutes to seconds
C   and from integer to real.
C
      DT=REAL(NCHEM)*60.D0 
C
C convert STT(kgs/box) to STT(molec/cc) (see SR convert_units.f)
C   kg CO/box * box/cc * mole/0.028 kg CO * Avog.#/mole
C sum loss in TCO(3) (molecules/box)
C sum production in TCO(4) (molecules/box)
C calculate new CO value 
C convert new GCO and copy value into STT 
C
C LPAUSE_MIN = minimun tropopause height.  Start L-loop from the 
C              lowest stratospheric level
C
      LPAUSE_MIN = MINVAL( LPAUSE )

      DO L = LPAUSE_MIN, LLPAR 
       DO J = 1, JJPAR
        DO I = 1, IIPAR

         IF ( L >= LPAUSE(I,J) ) THEN

c            print*,'I,J,L=',I,J,L
            STTTOGCO(I,J,L)=1./BOXVL(I,J,L)*XNUMOL_CO
c            print*,'STTTOGCO=',STTTOGCO(I,J,L)

c            print*,'STT(kg/box)=',STT(I,J,L,1)
            GCO(I,J,L)=STT(I,J,L,1)*STTTOGCO(I,J,L)
c            print*,'GCO(molec/cc)=',GCO(I,J,L)
c            print*,'CO_prod=',CO_prod(J,L)
c            print*,'CO_loss=',CO_loss(J,L)
c            print*,'BOXVL=',BOXVL(I,J,L)

            TCO(I,J,L,3)=TCO(I,J,L,3)+GCO(I,J,L)*BOXVL(I,J,L)*
     *           CO_loss(J,L)*DT
c      print*,'TCO(3)=',GCO(I,J,L)*BOXVL(I,J,L)*CO_loss(J,L)*DT

            TCO(I,J,L,4)=TCO(I,J,L,4)+bairdens(I,J,L)*CO_prod(J,L)*
     &                   DT*BOXVL(I,J,L)
c      print*,'TCO(4)=',bairdens(I,J,L)*CO_prod(J,L)*DT*BOXVL(I,J,L)

            GCO(I,J,L)=GCO(I,J,L)*(1-CO_loss(J,L)*DT)+
     *           CO_prod(J,L)*DT*bairdens(I,J,L)
c      print*,'GCO(molec/cc)=',GCO(I,J,L)
c      print*,'airdens=',bairdens(I,J,L)
c      print*,'DT=',DT
            
            STT(I,J,L,1)=GCO(I,J,L)/STTTOGCO(I,J,L)

c            print*,'STT(kg/box)=',STT(I,J,L,1)
c            print*,'*************'

         ENDIF

        ENDDO
       ENDDO
      ENDDO
C          stop
C
      RETURN
      END
