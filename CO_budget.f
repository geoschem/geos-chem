C $Id: CO_budget.f,v 1.1 2003/06/30 20:26:06 bmy Exp $
C*****************************************************************************
      SUBROUTINE CO_budget(LMN,FIRSTCHEM,NBUDGET)
C*****************************************************************************
      USE BPCH2_MOD
      IMPLICIT NONE

#     include "CMN_SIZE"       ! Size parameters
#     include "CMN"            ! STT, LPAUSE
#     include "CMN_CO"         ! BAIRDENS, CO arrays
#     include "CMN_CO_BUDGET"  ! TCO, XNUMOL_CO, NTALLY
#     include "CMN_DIAG"       ! KDADYN
      
C*****************************************************************************
      ! For binary punch file, version 2.0
      REAL*4             :: ARRAY(IIPAR, JJPAR, LLPAR)
      REAL*4             :: LONRES, LATRES

      INTEGER            :: IFIRST, JFIRST, LFIRST
      INTEGER, PARAMETER :: HALFPOLAR = 1
      INTEGER, PARAMETER :: CENTER180 = 1
      REAL*8 TAU1

      LOGICAL FIRSTCHEM
      CHARACTER (LEN=20) :: MODELNAME ,MERGEIT
      CHARACTER (LEN=40) :: UNIT
      CHARACTER (LEN=40) :: RESERVED = ''
      CHARACTER (LEN=40) :: CATEGORY 
      CHARACTER (LEN=80) :: LABEL

C*****************************************************************************
      INTEGER I,J,LMN,K,L,M,NBUDGET,NERROR,UD
      REAL*8 ATOTDT,STTCONV,TGS,SCALEDYN,NTP,NTQ,
     *     NTP2,NTQ2,SUMUP
      CHARACTER ( LEN=16 )  :: MERGE
      CHARACTER ( LEN=13 )  :: MERGE2
      REAL*8    BOXVL,SOURCES,SINKS
      EXTERNAL  BOXVL
      SCALEDYN  = FLOAT( KDADYN    ) + 1D-20
C*****************************************************************************
C This SR calculates the budget of CO.
C
C This SR only works for monthly averages, so be sure to
C start on the first of the month to another first of the
C month!!!
C*****************************************************************************
C LMN   = month (1-12).
C NYMDb = start date.
C NHMSb = stop date.
C JYEAR = year.
C*****************************************************************************
C Store the sources/sinks of CO in TCO in total molecules
C           ( 1) = initial burden
C           ( 2) = final burden
C  SINKS
C           ( 3) = CO sink by OH
C  SOURCES
C           ( 4) = CH4 oxidation
C           ( 5) = isoprene oxidation
C           ( 6) = biomass burning
C           ( 7) = wood burning
C           ( 8) = fossil fuels
C           ( 9) = monoterpene oxidation
C
C           (10) = interhemispheric exchange (+ = northward)
C           (11) = vacant
C           (12) = vacant
C*****************************************************************************
C
      print*,'Writing out CO budget in SR CO_budget.'
C
      IF(NBUDGET.EQ.1) THEN
C
C Convert initial burden of CO in TCO(1) from volume mixing
C  ratio to total molecules/box for first month of simulation.
C
         DO L=1,LLPAR
         DO J=1,JJPAR
         DO I=1,IIPAR

            TCO(I,J,L,1)=TCO(I,J,L,1)*bairdens(I,J,L)*BOXVL(I,J,L)

         ENDDO
         ENDDO
         ENDDO
C
      ENDIF
C
C*****************************************************************************
C Store the final burden of CO in TCO(2) 
C  Convert kg CO/box to molecules/box.
C
      TCO(:,:,:,2)=0.
      TCO(:,:,:,2)=STT(:,:,:,1)*XNUMOL_CO
C
C*****************************************************************************
C WRITE the monthly averages to output.
C*****************************************************************************
C BUDGET FOR LAYERS 1-20
C*****************************************************************************
C
C GLOBAL AVERAGES 
C
      WRITE(MERGE,2)LMN,JYEAR
 2    FORMAT('CObudget.',I2.2,'.',I4)
C
      OPEN(65,FILE=MERGE,STATUS='UNKNOWN')
      REWIND(65)
      
      TGS=1.D-9
      STTCONV=XNUMOL_CO/TGS
      SOURCES=0.D0
      SINKS=0.D0
      NERROR=0
      
      WRITE(65,18)
      WRITE(65,1801)
 1801 FORMAT('*************************')
      WRITE(65,1800)
 1800 FORMAT('LAYERS 1 - 20')
      WRITE(65,1801)
      WRITE(65,18)

      WRITE(65,18)
      WRITE(65,38)
      WRITE(65,18)
      WRITE(65,19)
      WRITE(65,1990)
 1990 FORMAT('Tropospheric Burden')
      NTP=SUMUP(1,IIPAR,1,JJPAR,1,LLPAR,1,1,1)
      NTP2=SUMUP(1,IIPAR,1,JJPAR,1,LLPAR,2,2,1)
      WRITE(65,20)NTP,NTP/STTCONV
      WRITE(65,21)NTP2,NTP2/STTCONV

      WRITE(65,18)
      WRITE(65,1991)
 1991 FORMAT('Stratospheric Burden')
      NTP=SUMUP(1,IIPAR,1,JJPAR,1,LLPAR,1,1,0)
      WRITE(65,20) NTP,NTP/STTCONV
      NTP=SUMUP(1,IIPAR,1,JJPAR,1,LLPAR,2,2,0)
      WRITE(65,21) NTP,NTP/STTCONV
      
      WRITE(65,18)
      WRITE(65,31)

c Sinks
      NTP=SUMUP(1,IIPAR,1,JJPAR,1,LLPAR,3,3,1)
      NTQ=SUMUP(1,IIPAR,1,JJPAR,1,LLPAR,3,3,0)
      SINKS=NTP+NTQ
      WRITE(65,22) NTP,NTP/SINKS*100.D0,NTP/STTCONV
      WRITE(65,220) NTQ,NTQ/SINKS*100.D0,NTQ/STTCONV
      
      WRITE(65,29) 
      WRITE(65,34) SINKS,SINKS/STTCONV
      WRITE(65,18)
      WRITE(65,30)

C Sources
      NTP=SUMUP(1,IIPAR,1,JJPAR,1,LLPAR,4,9,1)
      NTQ=SUMUP(1,IIPAR,1,JJPAR,1,LLPAR,4,4,0)
      SOURCES=NTQ+NTP
      
      NTP=SUMUP(1,IIPAR,1,JJPAR,1,LLPAR,4,4,1)
      WRITE(65,23) NTP,NTP/SOURCES*100.D0,NTP/STTCONV
      NTQ=SUMUP(1,IIPAR,1,JJPAR,1,LLPAR,4,4,0)
      WRITE(65,230) NTQ,NTQ/SOURCES*100.D0,NTQ/STTCONV
      NTP=SUMUP(1,IIPAR,1,JJPAR,1,1,5,5,1)
      WRITE(65,24) NTP,NTP/SOURCES*100.D0,NTP/STTCONV
      NTP=SUMUP(1,IIPAR,1,JJPAR,1,1,9,9,1)
      WRITE(65,39) NTP,NTP/SOURCES*100.D0,NTP/STTCONV
      NTP=SUMUP(1,IIPAR,1,JJPAR,1,1,6,6,1)
      WRITE(65,25) NTP,NTP/SOURCES*100.D0,NTP/STTCONV
      NTP=SUMUP(1,IIPAR,1,JJPAR,1,1,7,7,1)
      WRITE(65,26) NTP,NTP/SOURCES*100.D0,NTP/STTCONV
      NTP=SUMUP(1,IIPAR,1,JJPAR,1,1,8,8,1)
      WRITE(65,27) NTP,NTP/SOURCES*100.D0,NTP/STTCONV


      WRITE(65,29) 
      WRITE(65,28) SOURCES,SOURCES/STTCONV
      WRITE(65,18)

      NTP=SUMUP(1,IIPAR,1,JJPAR,1,LLPAR,1,1,1)
      NTP2=SUMUP(1,IIPAR,1,JJPAR,1,LLPAR,2,2,1)
      NTQ=SUMUP(1,IIPAR,1,JJPAR,1,LLPAR,1,1,0)
      NTQ2=SUMUP(1,IIPAR,1,JJPAR,1,LLPAR,2,2,0)
      WRITE(65,18)
      WRITE(65,288) NTP-NTP2+NTQ-NTQ2+SOURCES-SINKS,
     *     (NTP-NTP2+NTQ-NTQ2+SOURCES-SINKS)/STTCONV
 288  FORMAT('Initial-Final+Sources-Sinks=',E10.3,2x,F10.3)
      WRITE(65,18)
      WRITE(65,289) -(NTP-NTP2+NTQ-NTQ2),
     *     -(NTP-NTP2+NTQ-NTQ2)/STTCONV
 289  FORMAT('Net Gain          : ',E10.3,10x,F10.3)
      
 18   FORMAT()
 19   FORMAT('                    #Molecules               TG')
 20   FORMAT('  Start of Month  :',E10.3,10x,F10.3)
 21   FORMAT('  End of Month    :',E10.3,10x,F10.3)
 31   FORMAT('SINKS                            %Sink')
 22   FORMAT('  CO decay-trop   :',E10.3,2x,F6.1,2x,F10.3)
 220  FORMAT('  CO decay-strat  :',E10.3,2x,F6.1,2x,F10.3)
 34   FORMAT('Total Sinks       :',E10.3,10x,F10.3)
 30   FORMAT('SOURCES                          %Source')
 23   FORMAT('  CH4 Ox.-trop    :',E10.3,2x,F6.1,2x,F10.3)
 230  FORMAT('  CH4 Ox.-strat   :',E10.3,2x,F6.1,2x,F10.3)
 24   FORMAT('  Isop Oxidation  :',E10.3,2x,F6.1,2x,F10.3)
 39   FORMAT('  Mono Oxidation  :',E10.3,2x,F6.1,2x,F10.3)
 25   FORMAT('  Biomass Burning :',E10.3,2x,F6.1,2x,F10.3)
 26   FORMAT('  Wood Burning    :',E10.3,2x,F6.1,2x,F10.3)
 27   FORMAT('  Fossil Fuels    :',E10.3,2x,F6.1,2x,F10.3)
      
 270  FORMAT('  N-S Ex.-trop    :',E10.3,2x,F6.1,2x,F10.3)
 2700 FORMAT('  N-S Ex.-strat   :',E10.3,2x,F6.1,2x,F10.3)
 29   FORMAT('                     ---------')
 28   FORMAT('Total Sources     :',E10.3,10x,F10.3)
      
C
C SOUTHERN HEMISPHERE 
C
      SOURCES=0.D0
      SINKS=0.D0

      WRITE(65,18)
      WRITE(65,18)
      WRITE(65,36)
      WRITE(65,18)
      WRITE(65,19)
      WRITE(65,1990)
      NTP=SUMUP(1,IIPAR,1,JJPAR/2,1,LLPAR,1,1,1)
      WRITE(65,20) NTP,NTP/STTCONV
      NTP=SUMUP(1,IIPAR,1,JJPAR/2,1,LLPAR,2,2,1)
      WRITE(65,21) NTP,NTP/STTCONV
      WRITE(65,18)
      WRITE(65,1991)
      NTP=SUMUP(1,IIPAR,1,JJPAR/2,1,LLPAR,1,1,0)
      WRITE(65,20) NTP,NTP/STTCONV
      NTP=SUMUP(1,IIPAR,1,JJPAR/2,1,LLPAR,2,2,0)
      WRITE(65,21) NTP,NTP/STTCONV
      WRITE(65,18)
      WRITE(65,31)
c Sinks
      NTP=SUMUP(1,IIPAR,JJPAR/2+1,JJPAR/2+1,1,LLPAR,10,10,1)
      IF(NTP.GE.0.D0) THEN
         SINKS=SINKS+NTP
      ENDIF
      NTP2=SUMUP(1,IIPAR,JJPAR/2+1,JJPAR/2+1,1,LLPAR,10,10,0)
      IF(NTP2.GE.0.D0) THEN
         SINKS=SINKS+NTP2
      ENDIF
      
      NTQ=SUMUP(1,IIPAR,1,JJPAR/2,1,LLPAR,3,3,1)
      NTQ2=SUMUP(1,IIPAR,1,JJPAR/2,1,LLPAR,3,3,0)
      SINKS=SINKS+NTQ+NTQ2
      WRITE(65,22) NTQ,NTQ/SINKS*100.D0,NTQ/STTCONV
      WRITE(65,220) NTQ2,NTQ2/SINKS*100.D0,NTQ2/STTCONV
      
      NTP=SUMUP(1,IIPAR,JJPAR/2+1,JJPAR/2+1,1,LLPAR,10,10,1)
      IF(NTP.GE.0.D0) THEN
         WRITE(65,270) NTP,NTP/SINKS*100.D0,NTP/STTCONV
      ENDIF
      NTP=SUMUP(1,IIPAR,JJPAR/2+1,JJPAR/2+1,1,LLPAR,10,10,0)
      IF(NTP.GE.0.D0) THEN
         WRITE(65,2700) NTP,NTP/SINKS*100.D0,NTP/STTCONV
      ENDIF
      
      WRITE(65,29)
      WRITE(65,34) SINKS,SINKS/STTCONV
      WRITE(65,18)
      WRITE(65,30)
      
C Sources
      NTQ=SUMUP(1,IIPAR,1,JJPAR/2,1,LLPAR,4,9,1)
      NTQ2=SUMUP(1,IIPAR,1,JJPAR/2,1,LLPAR,4,4,0)
      SOURCES=NTQ+NTQ2
      NTP=SUMUP(1,IIPAR,JJPAR/2+1,JJPAR/2+1,1,LLPAR,10,10,1)
      IF(NTP.LT.0.D0) THEN
         SOURCES=SOURCES-NTP
      ENDIF
      NTP=SUMUP(1,IIPAR,JJPAR/2+1,JJPAR/2+1,1,LLPAR,10,10,0)
      IF(NTP.LT.0.D0) THEN
         SOURCES=SOURCES-NTP
      ENDIF
      
      NTP=SUMUP(1,IIPAR,1,JJPAR/2,1,LLPAR,4,4,1)
      NTP2=SUMUP(1,IIPAR,1,JJPAR/2,1,LLPAR,4,4,0)
      WRITE(65,23) NTP,NTP/SOURCES*100.D0,NTP/STTCONV
      WRITE(65,230) NTP2,NTP2/SOURCES*100.D0,NTP2/STTCONV
      NTP=SUMUP(1,IIPAR,1,JJPAR/2,1,1,5,5,1)
      WRITE(65,24) NTP,NTP/SOURCES*100.D0,NTP/STTCONV
      NTP=SUMUP(1,IIPAR,1,JJPAR/2,1,1,9,9,1)
      WRITE(65,39) NTP,NTP/SOURCES*100.D0,NTP/STTCONV
      NTP=SUMUP(1,IIPAR,1,JJPAR/2,1,1,6,6,1)
      WRITE(65,25) NTP,NTP/SOURCES*100.D0,NTP/STTCONV
      NTP=SUMUP(1,IIPAR,1,JJPAR/2,1,1,7,7,1)
      WRITE(65,26) NTP,NTP/SOURCES*100.D0,NTP/STTCONV
      NTP=SUMUP(1,IIPAR,1,JJPAR/2,1,1,8,8,1)
      WRITE(65,27) NTP,NTP/SOURCES*100.D0,NTP/STTCONV
      
      NTP=SUMUP(1,IIPAR,JJPAR/2+1,JJPAR/2+1,1,LLPAR,10,10,1)
      IF(NTP.LT.0.D0) THEN
         WRITE(65,270) -NTP,-NTP/SOURCES*100.D0,-NTP/STTCONV
      ENDIF
      NTP=SUMUP(1,IIPAR,JJPAR/2+1,JJPAR/2+1,1,LLPAR,10,10,0)
      IF(NTP.LT.0.D0) THEN
         WRITE(65,2700) -NTP,-NTP/SOURCES*100.D0,-NTP/STTCONV
      ENDIF
      
      WRITE(65,29)
      WRITE(65,28) SOURCES,SOURCES/STTCONV
      WRITE(65,18)
      
      NTP=SUMUP(1,IIPAR,1,JJPAR/2,1,LLPAR,1,1,1)
      NTP2=SUMUP(1,IIPAR,1,JJPAR/2,1,LLPAR,2,2,1)
      NTQ=SUMUP(1,IIPAR,1,JJPAR/2,1,LLPAR,1,1,0)
      NTQ2=SUMUP(1,IIPAR,1,JJPAR/2,1,LLPAR,2,2,0)
      WRITE(65,18)
      WRITE(65,288) (NTP-NTP2+NTQ-NTQ2+SOURCES-SINKS),
     *     (NTP-NTP2+NTQ-NTQ2+SOURCES-SINKS)/STTCONV
      WRITE(65,18)
      WRITE(65,289) -(NTP-NTP2+NTQ-NTQ2),
     *     -(NTP-NTP2+NTQ-NTQ2)/STTCONV
C     
C NORTHERN HEMISPHERE
C
      SOURCES=0.D0
      SINKS=0.D0

      WRITE(65,18)
      WRITE(65,18)
      WRITE(65,37)
      WRITE(65,18)
      WRITE(65,19)
      WRITE(65,1990)

      NTP=SUMUP(1,IIPAR,JJPAR/2+1,JJPAR,1,LLPAR,1,1,1)
      WRITE(65,20) NTP,NTP/STTCONV
      NTP=SUMUP(1,IIPAR,JJPAR/2+1,JJPAR,1,LLPAR,2,2,1)
      WRITE(65,21) NTP,NTP/STTCONV
      
      WRITE(65,18)
      WRITE(65,1991)
      NTP=SUMUP(1,IIPAR,JJPAR/2+1,JJPAR,1,LLPAR,1,1,0)
      WRITE(65,20) NTP,NTP/STTCONV
      NTP=SUMUP(1,IIPAR,JJPAR/2+1,JJPAR,1,LLPAR,2,2,0)
      WRITE(65,21) NTP,NTP/STTCONV
      
      WRITE(65,18)
      WRITE(65,31)
c Sinks
      NTQ=SUMUP(1,IIPAR,JJPAR/2+1,JJPAR,1,LLPAR,3,3,1)
      NTQ2=SUMUP(1,IIPAR,JJPAR/2+1,JJPAR,1,LLPAR,3,3,0)
      SINKS=NTQ+NTQ2

      NTP=SUMUP(1,IIPAR,JJPAR/2+1,JJPAR/2+1,1,LLPAR,10,10,1)
      IF(NTP.LT.0.D0) THEN
         SINKS=SINKS-NTP
      ENDIF
      NTP=SUMUP(1,IIPAR,JJPAR/2+1,JJPAR/2+1,1,LLPAR,10,10,0)
      IF(NTP.LT.0.D0) THEN
         SINKS=SINKS-NTP
      ENDIF

      WRITE(65,22) NTQ,NTQ/SINKS*100.D0,NTQ/STTCONV
      WRITE(65,220) NTQ2,NTQ2/SINKS*100.D0,NTQ2/STTCONV
      
      NTP=SUMUP(1,IIPAR,JJPAR/2+1,JJPAR/2+1,1,LLPAR,10,10,1)
      IF(NTP.LT.0.D0) THEN
         WRITE(65,270) -NTP,-NTP/SINKS*100.D0,-NTP/STTCONV
      ENDIF
      NTP=SUMUP(1,IIPAR,JJPAR/2+1,JJPAR/2+1,1,LLPAR,10,10,0)
      IF(NTP.LT.0.D0) THEN
         WRITE(65,2700) -NTP,-NTP/SINKS*100.D0,-NTP/STTCONV
      ENDIF
      WRITE(65,29)
      WRITE(65,34)SINKS,SINKS/STTCONV
      WRITE(65,18)
      WRITE(65,30)
C Sources
      NTQ=SUMUP(1,IIPAR,JJPAR/2+1,JJPAR,1,LLPAR,4,9,1)
      NTQ2=SUMUP(1,IIPAR,JJPAR/2+1,JJPAR,1,LLPAR,4,4,0)
      SOURCES=NTQ+NTQ2

      NTP=SUMUP(1,IIPAR,JJPAR/2+1,JJPAR/2+1,1,LLPAR,10,10,1)
      IF(NTP.GE.0.D0) THEN
         SOURCES=SOURCES+NTP
      ENDIF
      NTP=SUMUP(1,IIPAR,JJPAR/2+1,JJPAR/2+1,1,LLPAR,10,10,0)
      IF(NTP.GE.0.D0) THEN
         SOURCES=SOURCES+NTP
      ENDIF
      
      NTP=SUMUP(1,IIPAR,JJPAR/2+1,JJPAR,1,LLPAR,4,4,1)
      WRITE(65,23) NTP,NTP/SOURCES*100.D0,NTP/STTCONV
      NTP=SUMUP(1,IIPAR,JJPAR/2+1,JJPAR,1,LLPAR,4,4,0)
      WRITE(65,230) NTP,NTP/SOURCES*100.D0,NTP/STTCONV
      
      NTP=SUMUP(1,IIPAR,JJPAR/2+1,JJPAR,1,1,5,5,1)
      WRITE(65,24) NTP,NTP/SOURCES*100.D0,NTP/STTCONV
      NTP=SUMUP(1,IIPAR,JJPAR/2+1,JJPAR,1,1,9,9,1)
      WRITE(65,39) NTP,NTP/SOURCES*100.D0,NTP/STTCONV
      NTP=SUMUP(1,IIPAR,JJPAR/2+1,JJPAR,1,1,6,6,1)
      WRITE(65,25) NTP,NTP/SOURCES*100.D0,NTP/STTCONV
      NTP=SUMUP(1,IIPAR,JJPAR/2+1,JJPAR,1,1,7,7,1)
      WRITE(65,26) NTP,NTP/SOURCES*100.D0,NTP/STTCONV
      NTP=SUMUP(1,IIPAR,JJPAR/2+1,JJPAR,1,1,8,8,1)
      WRITE(65,27) NTP,NTP/SOURCES*100.D0,NTP/STTCONV
      NTP=SUMUP(1,IIPAR,JJPAR/2+1,JJPAR/2+1,1,LLPAR,10,10,1)
      IF(NTP.GE.0.D0) THEN
         WRITE(65,270) NTP,NTP/SOURCES*100.D0,NTP/STTCONV
      ENDIF
      NTP=SUMUP(1,IIPAR,JJPAR/2+1,JJPAR/2+1,1,LLPAR,10,10,0)
      IF(NTP.GE.0.D0) THEN
         WRITE(65,2700) NTP,NTP/SOURCES*100.D0,NTP/STTCONV
      ENDIF
      WRITE(65,29)
      WRITE(65,28) SOURCES,SOURCES/STTCONV
      WRITE(65,18)
      NTP=SUMUP(1,IIPAR,JJPAR/2+1,JJPAR,1,LLPAR,1,1,1)
      NTP2=SUMUP(1,IIPAR,JJPAR/2+1,JJPAR,1,LLPAR,2,2,1)
      NTQ=SUMUP(1,IIPAR,JJPAR/2+1,JJPAR,1,LLPAR,1,1,0)
      NTQ2=SUMUP(1,IIPAR,JJPAR/2+1,JJPAR,1,LLPAR,2,2,0)
      WRITE(65,18)
      WRITE(65,288) (NTP-NTP2+NTQ-NTQ2+SOURCES-SINKS),
     *     (NTP-NTP2+NTQ-NTQ2+SOURCES-SINKS)/STTCONV
      WRITE(65,18)
      WRITE(65,289) -(NTP-NTP2+NTQ-NTQ2),
     *     -(NTP-NTP2+NTQ-NTQ2)/STTCONV
 36   FORMAT('*****  Southern Hemisphere  *****')
 37   FORMAT('*****  Northern Hemisphere  *****')
 38   FORMAT('*****  Global  *****')
C
      CLOSE(65)
C
C*****************************************************************************
C Create punch file for TCO .
C For binary punch file, version 2.0
C
       print*,'Creating CO punch file of budget info in SR CO_budget.'
      ! Make up a category name for GAMAP (use 8 characters)
C            1='COinit'
C            2='COfin'
C            3='COsink'
C            4='COprod'
C            5='COisop'
C            6='CObiom'
C            7='CObiof'
C            8='COfoss'
C            9='COmono'
C           10='COequa'
C           11='      '
C           12='      '
C
      IFIRST = I0 + 1
      JFIRST = J0 + 1
      LFIRST = 1
      LONRES = DISIZE
      LATRES = DJSIZE

      ! Get the proper model name for the binary punch file
      MODELNAME = GET_MODELNAME()

       WRITE(MERGEIT,2001)LMN,JYEAR
 2001  FORMAT('COpch.',I2.2,'.',I4)

       OPEN(65,file=MERGEIT,
     &        form='unformatted',access='sequential')
       REWIND(65)

      ! Descriptor string
      LABEL = 'GEOS-CHEM -- CO Budget output (bnd, jsw, 10/13/00)'

      ! Unit of quantity being saved
      UNIT   = 'Tg'
C
C Convert TCO (molec/box) to (Tg/box).
C
         DO L=1,LLPAR
         DO J=1,JJPAR
         DO I=1,IIPAR
         DO K=1,NTALLY
            TCO(I,J,L,K)=TCO(I,J,L,K)/STTCONV
         ENDDO
         ENDDO
         ENDDO
         ENDDO
C
      CALL BPCH2_HDR( 65, LABEL )

      CATEGORY = 'COBUDGET'

      DO K = 1, NTALLY-2
        

         ! Cast REAL*8 into REAL*4
         ARRAY(:,:,:) = TCO(:,:,:,K)
        
         ! Save the data block
         CALL BPCH2( 65,        MODELNAME, LONRES,   LATRES,
     &               HALFPOLAR, CENTER180, CATEGORY, K,    
     &               UNIT,      TAU0,      TAU,      RESERVED,   
     &               IIPAR,     JJPAR,     LLPAR,    IFIRST,     
     &               JFIRST,    LFIRST,    ARRAY )
      ENDDO

      CLOSE ( 65 )
C
C*****************************************************************************
C Final burden at last of month equals initial burden
C of next month. 
C
C Convert TCO (Tg/box) back to (molec/box).
C
         DO L=1,LLPAR
         DO J=1,JJPAR
         DO I=1,IIPAR
         DO K=2,2
            TCO(I,J,L,K)=TCO(I,J,L,K)*STTCONV
         ENDDO
         ENDDO
         ENDDO
         ENDDO
C
      TCO(:,:,:,1) = TCO(:,:,:,2)
C
C Set TCO = 0 for next month.
C
      TCO(:,:,:,2:NTALLY)=0.D0
C
C*****************************************************************************
C
      print*,'Leaving SR CO_budget.'
      RETURN
      END
