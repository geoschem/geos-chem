C $Id: CO_fromHCs.f,v 1.1 2003/06/30 20:26:04 bmy Exp $
      SUBROUTINE CO_fromHCs(BOH)
C*****************************************************************************

      ! References to F90 modules
      USE DAO_MOD,        ONLY : AIRVOL
      USE GLOBAL_CH4_MOD, ONLY : GET_GLOBAL_CH4

      IMPLICIT NONE

#     include "CMN_SIZE"       ! Size parameters
#     include "CMN"            ! STT, LPAUSE
#     include "CMN_CO"         ! CO arrays
#     include "CMN_CO_BUDGET"  ! TCO, FMOL_CO, XNUMOL_CO, NTALLY

      INTEGER I,J,K,L,M,N,NISOPDO,NISOPALPHA
      PARAMETER(NCOHCS=3)
      REAL*8    BOXVL,ANOX
      EXTERNAL  BOXVL

      REAL*8 DT,ALPHA(IIPAR,JJPAR,LLPAR),ALPHAch4,ALPHAMONO,
     * BOH(IIPAR,JJPAR,LLPAR), 
      
      ! Now declare these as SAVE (bmy, 1/3/01)
      REAL*8, SAVE  :: A3090S, A0030S, A0030N, A3090N

      ! Variable to hold the maximum value of LPAUSE (bmy, 4/18/00)
      INTEGER       :: LPAUSE_MAX

      ! Save the previous year
      INTEGER, SAVE :: LASTYEAR = -1

      ! Flag for selecting variable or constant CH4 concentrations
      LOGICAL       :: VARIABLE_CH4 = .TRUE.

c      print*,'Calculating CO from HCs in SR CO_fromHCs.'
C
C ****************************************************************************
C Calulates the production of CO from hydrocarbons.
C    A) CH4 
C    B) Isoprene 
C    C) Monoterpenes
C    
C The production of CO from VOCs from anthropogenic activities is accounted
C    for in SR emissco.
C
C   Created by Bryan Duncan.
C   Updated by Bob Yantosca (1/3/01)
C ****************************************************************************
C ****************************************************************************
C NCOHCS = # of hydrocarbons considered.
C SUMCO(:,:,:,:) = molec/cc of CO produced by the
C           oxidation of species in time step, DT
C STT(:,:,:,1) = CO in kgs/box
C GCO(:,:,:) = [CO] in molec/cc
C GCH4(:,:,:) = [CH4] in molec/cc
C GISOP(:,:,:) = [ISOP] in molec/cc
C ALPHA = molec CO/molec X produced during X oxidation
C KRATE(:,:,:) = reaction rate constant
C ****************************************************************************
C The annual mean tropopause is stored in the LPAUSE array 
C (from header file "CMN").  LPAUSE is defined such that: 
C
C Levels            1 <= L <= LPAUSE(I,J) - 1 are tropospheric
C         LPAUSE(I,J) <= L <= LLPAR           are stratospheric
C ****************************************************************************
C convert the chemistry time step (NCHEM) from minutes to seconds
C   and from integer to real.
C
      DT=REAL(NCHEM)*60.D0
C
C*****************************************************************************
C A) CH4 -- now see comments in "global_ch4_mod.f" (bmy, 1/3/01)
C
C    If VARIABLE_CH4 = T, then select CH4 gradient for year JYEAR
C    If VARIABLE_CH4 = F, then use 1700 ppbv CH4 everywhere
C*****************************************************************************
      
      ! Compute the latitudinal CH4 values once per year (bmy, 1/3/01)
      IF ( JYEAR /= LASTYEAR ) THEN
         CALL GET_GLOBAL_CH4( JYEAR, VARIABLE_CH4, 
     &                        A3090S, A0030S, A0030N, A3090N )
         LASTYEAR = JYEAR
      ENDIF

C************
C
C rate constant from DeMore et al. 1997
C k = 2.45x10-12*exp(-1775/T)
C
C put [CH4] into molec/cc
C
      GCH4(:,:,:)=0.D0
C
      GCH4(:,1:12,:) =A3090S*1.D-9*bairdens(:,1:12,:)
      GCH4(:,13:23,:)=A0030S*1.D-9*bairdens(:,13:23,:)
      GCH4(:,24:34,:)=A0030N*1.D-9*bairdens(:,24:34,:)
      GCH4(:,35:46,:)=A3090N*1.D-9*bairdens(:,35:46,:)
C
C calculate rate constant.
C
      KRATE(:,:,:)=0.D0
C
      KRATE(:,:,:)=2.45D-12*EXP(-1775.D0/Tavg(:,:,:))
C
C ALPHA for CH4 estimated 0.95-1
C
      ALPHAch4 = 1.D0
C
C # molec CO/cc = ALPHA*k*[CH4]*[OH]*dt
C
      SUMCO(:,:,:,:)=0.D0
C
      SUMCO(:,:,:,1)=ALPHAch4*KRATE(:,:,:)*
     *     GCH4(:,:,:)*BOH(:,:,:)*DT
C
C Store # molec/time step produced in TCO(4).
C
      ! Store the maximum value of LPAUSE in LPAUSE_MAX for use
      ! in the DO-loops below (bmy, 4/18/00)
      LPAUSE_MAX = MAXVAL( LPAUSE )

      DO L=1,LPAUSE_MAX
      DO J=1,JJPAR
      DO I=1,IIPAR

         IF ( L < LPAUSE(I,J) ) THEN

            TCO(I,J,L,4)=TCO(I,J,L,4)+
     *           SUMCO(I,J,L,1)*BOXVL(I,J,L)

         ENDIF

      ENDDO
      ENDDO
      ENDDO
C
C*****************************************************************************
C B) ISOP (by OH only - neglect O3)
C*****************************************************************************
C About 397 Tg C of isoprene is emitted annually.
cbnd      NISOPALPHA=1
      NISOPALPHA=0
C
      IF(NISOPALPHA.EQ.0) THEN
C     
         ALPHA(:,:,:) = 0.D0
C     
C The CO yield here is taken from acf's calculation (from group
C   meeting (1-20-99).  Assumming linearity, ALPHA=0.8*[NOx]+0.6,
C   with an upper and lower limit of 2.1 and 0.8, resp.
C
         DO L=1,LLPAR
         DO J=1,JJPAR
         DO I=1,IIPAR

            ANOX=BBIJ(I,J,L,1)*1.D9 

            ALPHA(I,J,L) = 0.8D0*ANOX+0.6D0
        
            IF(ANOX.LT.0.5D0) THEN
               ALPHA(I,J,L) = 0.8D0
            ENDIF

            IF(ANOX.GT.1.8D0) THEN
               ALPHA(I,J,L) = 2.1D0
            ENDIF

         ENDDO
         ENDDO
         ENDDO

      ELSE
C
C 30% yield from Miyoshi et al., 1994.
C They estimate globally 105 Tg C/yr of CO is produced from isoprene oxidation.
C Increased yield from 30% to 50% (bnd, bmy, 1/3/01)
C
         ALPHA(:,:,:)=1.5D0
C
      ENDIF
C
C*****************************************************************************
C
      NISOPDO=1
C
C*****************************************************************************
      IF(NISOPDO.EQ.0) THEN      
C*****************************************************************************
C
C Calculate source of CO per time step by:
C   1) from Isoprene average fields:
C       # molec CO/cc = ALPHA*k*[ISOP]*[OH]*dt
C
C put [ISOP] into molec/cc
C
         GISOP(:,:,:)=0.D0
C     
         GISOP(:,:,:)=BBIJ(:,:,:,3)*0.2D0*bairdens(:,:,:)
C     
C calculate rate constant.
C
         KRATE(:,:,:)=0.D0
C
         KRATE(:,:,:)=2.54D-11*EXP(410.D0/Tavg(:,:,:))
C
         SUMCO(:,:,:,2)=ALPHA(:,:,:)*KRATE(:,:,:)*
     *        GISOP(:,:,:)*BOH(:,:,:)*DT
C
C Store # molec/box produced in TCO(5).
C     
         ! The L-loop only needs to go up to LPAUSE_MAX (bmy, 4/18/00)
         !DO L=1,LLPAR
         DO L=1,LPAUSE_MAX
         DO J=1,JJPAR
         DO I=1,IIPAR

            ! Only process tropospheric boxes (bmy, 4/17/00)
            IF ( L < LPAUSE(I,J) ) THEN

               TCO(I,J,L,5)=TCO(I,J,L,5)+SUMCO(I,J,L,2)*BOXVL(I,J,L)
               
            ENDIF

         ENDDO
         ENDDO
         ENDDO
C
C*****************************************************************************
      ELSE
C*****************************************************************************
C 
C   2) from Isoprene flux (assume lifetime very short) 
C       # molec CO/cc = ALPHA*Flux(ISOP)*dt
C
C    CO from Isoprene:  
C    Assume the production of CO from isoprene is instantaneous
C    even though the lifetime of intermediate species may be on
C    the order of hours or days.  This assumption will likely
C    cause CO from isoprene oxidation to be too high in the 
C    box in which the isoprene is emitted.
C
C SUMISOPCO is calculated in SR emissco (units=atoms C/box/time step).
C
         SUMCO(:,:,:,2)=0.D0
C
C SUMCO = molecules CO / cm^3 / time step
         DO J=1,JJPAR
         DO I=1,IIPAR
            SUMCO(I,J,1,2)=SUMISOPCO(I,J)/BOXVL(I,J,1)*0.2D0*
     *           ALPHA(I,J,1)
         ENDDO
         ENDDO
C
C Store # molec CO produced per time step in TCO(5).
C
         DO J=1,JJPAR
         DO I=1,IIPAR
            TCO(I,J,1,5)=TCO(I,J,1,5)+SUMISOPCO(I,J)
     *           *0.2D0*ALPHA(I,J,1)
         ENDDO
         ENDDO
C
         SUMISOPCO(:,:)=0.D0
C
C*****************************************************************************
      ENDIF
C*****************************************************************************
C
C*****************************************************************************
C C) Monoterpenes (OH & O3)
C
C    CO from Monoterpenes:  
C    Assume the production of CO from monoterpenes is instantaneous
C    even though the lifetime of intermediate species may be on
C    the order of hours or days.  This assumption will likely
C    cause CO from monoterpene oxidation to be too high in the 
C    box in which the monoterpene is emitted.
C*****************************************************************************
C
C The CO yield here is taken from: 
C   Hatakeyama et al. JGR, Vol. 96, p. 947-958 (1991)
C   Vinckier et al. Fresenius Env. Bull., Vol. 7, p.361-368 (1998)
C
C   Hatakeyama:  "The ultimate yield of CO from the tropospheric
C      oxidation of terpenese (including both O3 and OH reactions)
C      was estimated to be 20% on the carbon number basis."  They
C      studied ALPHA- & beta-pinene.
C   Vinckier  :  "R(CO)=1.8+/-0.3" : 1.8/10 is about 20%.
C
      ALPHAMONO = 0.2D0
C
C*****************************************************************************
C
C Calculate source of CO per time step
C       from monoterpene flux (assume lifetime very short) 
C       # molec CO/cc = ALPHA*Flux(MONO)*dt
C
C SUMMONOCO is calculated in SR emissco (units=atoms C/box/time step).
C
      SUMCO(:,:,:,3)=0.D0
C     
      DO J=1,JJPAR
      DO I=1,IIPAR
         SUMCO(I,J,1,3)=SUMMONOCO(I,J)/BOXVL(I,J,1)*
     *        ALPHAMONO
      ENDDO
      ENDDO
C
C Store # molec/box produced in TCO(9).
C
      DO J=1,JJPAR
      DO I=1,IIPAR
         TCO(I,J,1,9)=TCO(I,J,1,9)+SUMMONOCO(I,J)
     *        *ALPHAMONO
      ENDDO
      ENDDO
C
      SUMMONOCO(:,:)=0.D0
C
C*****************************************************************************
C calculate new CO value : [CO]=[CO](1+SUMCO(CH4,ISOP,MONO))
C*****************************************************************************
C convert CO's STT(kgs/box) to GCO(molec/cc) (see SR convert_units.f)
C   kg CO/box * box/cc * mole/0.028 kg CO * Avog.#/mole
C
      STTTOGCO(:,:,:)=0.D0
C
      STTTOGCO(:,:,:)=1.D0/AIRVOL(:,:,:)/1.D6/FMOL_CO*6.023D23
C
      GCO(:,:,:)=0.D0
C
      GCO(:,:,:)=STT(:,:,:,1)*STTTOGCO(:,:,:)
C     
C sum up CO contributions (SUMCO) and add to GCO
C
      totco(:,:,:)=0.D0
C
      DO J=1,NCOHCS
         totco(:,:,:)=totco(:,:,:)+SUMCO(:,:,:,J)
      ENDDO
C
      GCO(:,:,:)=GCO(:,:,:)+totco(:,:,:)
C
C convert new GCO and copy value into STT
C
      ! The L-loop only needs to go up to LPAUSE_MAX (bmy, 4/18/00)
      DO L=1,LPAUSE_MAX
      DO J=1,JJPAR
      DO I=1,IIPAR

         ! Only process tropospheric boxes (bmy, 4/17/00)
         IF ( L < LPAUSE(I,J) ) THEN
            
            STT(I,J,L,1)=GCO(I,J,L)/STTTOGCO(I,J,L)
            
         ENDIF

      ENDDO
      ENDDO
      ENDDO
C
      RETURN
      END
