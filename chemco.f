C $Id: chemco.f,v 1.1 2003/06/30 20:26:08 bmy Exp $
C*****************************************************************************
      SUBROUTINE chemco( FIRSTCHEM, LMN, NSEASON, NYMDb, NYMDe )
C*****************************************************************************
C                                                                     
C  Chemistry subroutine for carbon monoxide (CO) as a single tracer.  
C   Bryan Duncan
C
C  SINKS:
C   Removal of CO by OH (SR OHparam & CO_decay).
C   CO uptake by soils (neglected).
C
C  SOURCES:
C   Oxidation of methane, isoprene and monoterpenes (SRCO_fromHCs).
C   Direct emissions of CO from fossil fuel combustion, biomass 
C      burning and biofuels use.
C                                                                     
C*****************************************************************************
C ------- FOLLOWING VARIABLES ARE DEFINED FOR EACH RUN ---------------
C  FIRSTCHEM = .TRUE.  FOR THE FIRST TIME DOING A PARTICULAR GRID-BOX
C              USES CONCENTRATIONS FROM chem.dat
C  FIRSTCHEM = .FALSE.  FOR SUBSEQUENT TIMES USES CONCENTRATIONS FROM
C              CSPEC   NOTE: IF NPTS > 1 ALL GRID-BOX USE SAME FIRSTCHEM
C  NPTS      = NUMBER OF POINTS (GRID-BOXES) TO CALCULATE
C ------- REMAINING VARIABLES ARE DEFINED FOR EACH GRID-BOX ----------
C  CTRACE    = INITIAL TRACER CONCENTRATIONS (molec cm^-3)
C  REMIS     = EMISSION RATES (molec cm^-3 sec^-1)
C  LMN       = month
C  LMN_last  = previous month
C****** INPUT FILES **********************************************************
C The annual mean tropopause is stored in the LPAUSE array 
C (from header file "CMN").  LPAUSE is defined such that: 
C
C Levels            1 <= L <= LPAUSE(I,J) - 1 are tropospheric
C         LPAUSE(I,J) <= L <= LLPAR           are stratospheric
C
C We now use LPAUSE instead of NSKIPL to denote the strat/trop boundary. 
C (bmy, 4/18/00)  
C*****************************************************************************
      
      ! References to F90 modules
      USE DIAG_MOD,     ONLY : AD43
      USE GRID_MOD,     ONLY : GET_YMID  ! (bmy, 2/3/03)
      USE UVALBEDO_MOD, ONLY : UVALBEDO

      IMPLICIT NONE

#     include "CMN_SIZE"   ! Size parameters
#     include "CMN"        ! STT, LPAUSE
#     include "CMN_DIAG"   ! ND43, AD43
#     include "CMN_CO"     ! CO arrays
      
      LOGICAL FIRSTCHEM
      INTEGER NOHDO,NBUDGET,NSKIPTOMS
      COMMON /BUDGETSTUFF/ NBUDGET

      REAL*8 BOH(IIPAR,JJPAR,LLPAR),SUMUP1,SUMUP2

      INTEGER I,J,L,K,M,N,IJ,JJ,NPART,III,JJJ
      INTEGER IDXAIR(JJPAR),IDXO3(JJPAR),NYMDb,NYMDe
      INTEGER LMN,NSEASON
      REAL*8  OPTD(LLPAR,IIPAR,JJPAR)
      
      INTEGER NODAYS(12)
      DATA NODAYS/31,28,31,30,31,30,31,31,30,31,30,31/
C
      REAL*8      BOXVL
      EXTERNAL    BOXVL

      ! Weight of air (taken from "comode.h") 
      REAL*8, PARAMETER :: WTAIR = 28.966d0
C
      WRITE(*,*) ''
      WRITE(*,*) '*************  ENTERING CHEMCO!  *************'
      WRITE(*,*) ''
C
C ****************************************************************************
C For biomass burning, chose one of the following options in input.geos.
C The decision tree is as follows:
C
C #    LBBSEA   LTOMSAI    Action
C------------------------------------------------------------------------
C 1      T        F       Read seasonal biomass burning data from disk
C
C 2      T        T       Read seasonal biomass burning data from disk,
C                          but apply TOMS AI for interannual variability 
C
C 3      F        T       Same as #2
C
C 4      F        F       Read Randall's interannual variability biomass
C                          burning emissions directly from disk
C ****************************************************************************
C If the first time step ...
C ****************************************************************************
C
      IF(FIRSTCHEM) THEN
C
C Counter for total number of timesteps per month for CO budget.
C
         NTALDT=1
C
C Flag for SR CO_budget
C
         NBUDGET=1
C
C  LMN = month; LMN_last = previous month
C  NSEASON:  is the seasonal index for NOx emissions:
C        NSEASON=1 --> winter (Dec, Jan, Feb)
C        NSEASON=2 --> spring (Mar, Apr, May)
C        NSEASON=3 --> summer (Jun, Jul, Aug)
C        NSEASON=4 --> autumn (Sep, Oct, Nov)
C  NSEASON_last2 = previous season
C
C Initialize LMN_last = 99, so SR FILLFIELDS will be called on the
C   initial time step.
C
         LMN_last=99 
         NSEASON_last2=99
C
      ENDIF                     ! IF (FIRSTCHEM)
C ****************************************************************************
C
      NTALDT=NTALDT+1
C
C ****************************************************************************
C calculate each box's air density (molec/cc)
C ****************************************************************************
C
      DO L=1,LLPAR
      DO J=1,JJPAR
      DO I=1,IIPAR
         bairdens(I,J,L)=AD(I,J,L)*1000./BOXVL(I,J,L)*AVGNO/WTAIR
      ENDDO
      ENDDO
      ENDDO
C
C ****************************************************************************
C  Calculate the production and destruction of CO from gas-phase
C     chemistry only.       
C ****************************************************************************
C Concerning O3, there are 3 options:
C   A) The OH parameterization is calculated using GEOS monthly means
C      (NCLIMATOLOGY=0) for the independent variable O3.  The o3 column
C      above independent variable is determined using jal's O3 
C      climatologies for both the tropospheric and stratospheric
C      portions of the O3 column (NCLIMATOLOGY2=1). 
C   B) The O3 variable is determined from jal's O3 climatolgies 
C      (tropospheric portion) and the o3 column above variable is
C      determined from jal's O3 climatolgies (NCLIMATOLOGY=1 & 
C      NCLIMATOLOGY2=1).
C   C) The O3 variable is determined by GEOS monthly means (NCLIMATOLOGY=0)
C      and the o3 column above variable is determined from fastj
C      climatologies for both the tropospheric and stratospheric
C      portion of the O3 column (NCLIMATOLOGY=0 & NCLIMATOLOGY2=0).
C
      NCLIMATOLOGY=0
C
      NCLIMATOLOGY2=1
C
C If you would like to scale the O3 columns to observed TOMS data
C then set NCLIMATOLOGY3=1; otherwise, set it equal to 0.  Note
C that there are data gaps in the TOMS data where a climatology
C will be used (see SR read_TOMSO3 for details).
C
      NCLIMATOLOGY3=1
C
C Error Check.
      IF(NCLIMATOLOGY.EQ.1.AND.NCLIMATOLOGY2.EQ.0) THEN
         PRINT*,'Stopped in SR CHEMCO.'
         PRINT*,'This combination of options is not allowed!'
         PRINT*,'Reset NCLIMATOLOGY and/or NCLIMATOLOGY2.'
         STOP
      ENDIF
C
C ****************************************************************************
C 1) read in species fields (i.e., monthly average fields) and total monthly
C    ozone columns.
C ****************************************************************************
C
      print*,'Month for this time step =',LMN
      print*,'Month for previous time step =',LMN_last
 
      IF(LMN.NE.LMN_last) THEN
C
         CALL CO_fillfields(LMN,NCLIMATOLOGY)
C
C Read in TOMS total ozone columns.
C
         IF(NCLIMATOLOGY3.EQ.1) THEN
           NSKIPTOMS=0
              CALL read_TOMSO3(NSKIPTOMS)
           IF(NSKIPTOMS.EQ.1) NCLIMATOLOGY3=0
         ENDIF
C
C Read in the loss/production of CO in the stratosphere.
C
         CALL read_COPminusL(LMN)
C
      ENDIF
C
C Read in fastj total ozone columns.
C
      IF(FIRSTCHEM) THEN
       IF(NCLIMATOLOGY2.EQ.0) THEN
         CALL CO_readO3()
       ENDIF
      ENDIF
C
C ****************************************************************************
C 2) calculate the O3 column (DU) above.
C ****************************************************************************
C
      IF(NCLIMATOLOGY2.EQ.1) THEN
C
cc         IF(NSEASON.NE.NSEASON_last2) THEN
       IF(LMN.NE.LMN_last) THEN

        IF(FIRSTCHEM) THEN
C
         CALL CO_read_clim()
C
        ENDIF
C
            CALL o3above(NSEASON,NCLIMATOLOGY,NCLIMATOLOGY2)
C
            NSEASON_last2 = NSEASON
C
            IF(NSEASON_last2.NE.99.and.NCLIMATOLOGY.EQ.1) THEN
               print*,'SR chemco'
               print*,'Check to make sure O3 fields are'
               print*,' still climatology after 2nd season!'
            ENDIF
C
         ENDIF
C
      ELSE
C
         IF(LMN.NE.LMN_last) THEN
C     
            IF(FIRSTCHEM) THEN
C
               CALL CO_read_fastj(NJVREAD,'jv_atms.dat')
C     
            ENDIF
C
C o3up = O3 column (DU) above.
C
            o3up(:,:,:)=0.
C
            DO JJJ = 1, JJPAR
            DO III = 1, IIPAR
C
               CALL CO_set_prof(P(III,JJJ),LMN,GET_YMID(JJJ),SIGE,PTOP)
C     
               CALL fastjO3col(o3up,JJJ,III)
C
            ENDDO
            ENDDO
C
         ENDIF
C
      ENDIF
C
            LMN_last = LMN
C
C ****************************************************************************
C 3) get parameterized OH fields or monthly mean fields.
C ****************************************************************************
C
C BOH = storage array for OH fields.
C
      BOH(:,:,:)=0.
C
C NOHDO = switch
C       = 0 : Use OH field from full chemistry monthly average field.
C       = 1 : Get parameterized OH field.
C
      NOHDO=1
C
      IF(NOHDO.EQ.0) THEN
C     
         BOH(:,:,:)=BBIJ(:,:,:,NFIELDS2)
C     
      ELSE
C
         CALL GETINFO( LMN, UVALBEDO )
C  
C LPAUSE = the vertical level of the tropopause.  Above this level,
C          no [OH] is calculated.  The user can feed this SR
C          a high value for LPAUSE which effectively turns
C          this option off (i.e., LPAUSE > MVRTBX). If the 
C          [OH] = -999 then the [OH] was not calculated.
C
         CALL OHparam( BOH, LPAUSE, LMN )
C
      ENDIF
C
C*****************************************************************************
C  ND43 diagnostics...save [OH] in molecules/cm3
C*****************************************************************************
C
      IF ( ND43 .GT. 0 ) THEN
C
         DO L = 1, LD43
         DO J = 1, JJPAR
         DO I = 1, IIPAR
C
            ! Only process tropospheric boxes (bmy, 4/17/00)
            IF ( L < LPAUSE(I,J) ) THEN
C
              AD43(I,J,L,1) = AD43(I,J,L,1) + BOH(I,J,L)
C     
            ENDIF
C
         ENDDO
         ENDDO
         ENDDO
C
      ENDIF
C
C ****************************************************************************
C  Diagnostic for Methyl Chloroform (ND23).
C ****************************************************************************
C bey recommends to calculate offline.
      CALL CO_OHSAVE(BOH)
C
C ****************************************************************************
C 4) calculate CO production from HC oxidation.
C ****************************************************************************
C
      CALL CO_fromHCs(BOH)
C
C ****************************************************************************
C 5) calculate rate of decay of CO by OH oxidation.
C ****************************************************************************
C
      CALL CO_decay(BOH)
C
C ****************************************************************************
C 7) do CO chemistry in layers above tropopause.
C ****************************************************************************
C
      CALL CO_strat() 
C
C ****************************************************************************
C 6) write budget (i.e., monthly average fields).
C ****************************************************************************
C
C Check to make sure the start and end times are on the
C   first of a month.  If not the SR CO_budget will not
C   work properly!
C      
      NPART=NYMDb/100 
      IF((NYMDb-NPART*100).NE.1) THEN
         print*,'Start date not equal to 1st of month!!!'
         print*,'  Therefore, SR CO_budget will not work!!!'
         STOP
      ENDIF
      
      NPART=NYMDe/100 
      IF((NYMDe-NPART*100).NE.1) THEN
         print*,'End date not equal to 1st of month!!!'
         print*,'  Therefore, SR CO_budget will not work!!!'
         STOP
      ENDIF
C
C Leap year correction for February.
      IF (JYEAR.EQ.1988.AND.LMN.EQ.2) NODAYS(LMN)=29
      IF (JYEAR.EQ.1992.AND.LMN.EQ.2) NODAYS(LMN)=29
      IF (JYEAR.EQ.1996.AND.LMN.EQ.2) NODAYS(LMN)=29
      IF (JYEAR.EQ.2000.AND.LMN.EQ.2) NODAYS(LMN)=29
C
      IF (NTALDT.EQ.NODAYS(LMN)) THEN

         CALL CO_budget(LMN,FIRSTCHEM,NBUDGET)

         NBUDGET=0

	 NTALDT=0

C Reset February to 28 days for non-leap years.
      IF (LMN.EQ.2.AND.NODAYS(LMN).EQ.29) NODAYS(LMN)=28
         
      ENDIF
C
C ****************************************************************************
C
      WRITE(*,*) ''
      WRITE(*,*) '************* LEAVING CHEMCO!  *************'
      WRITE(*,*) ''
C     
      ! Set FIRSTCHEM = .FALSE. before leaving CHEMCO (bmy, 4/18/00)
      FIRSTCHEM = .FALSE. 

      RETURN
      END
