C $Id: read_coeff.f,v 1.1 2003/06/30 20:26:07 bmy Exp $
C **********************************************************************
      SUBROUTINE READ_COEFF 
C **********************************************************************
      IMPLICIT NONE
#     include "CMN_OH"
      INTEGER JTERM,M,PARSUM,NFT
      INTEGER J,L,K,KK,LH,KR,IDONE
C **********************************************************************
C Created by Bryan Duncan.
C **********************************************************************
C This SR reads in information from the file 'start_information'. 
C The information is needed to calculate the OH using parameterizations
C and includes polynomial coefficients, ranges of independent variables,
C etc.  See the text of this code for descriptions of information. 
C *********************************************************************
      print*,'Reading info for parameterization in SR read_coeff.'
C
C NGENERICFILE = filenumber specified by user in CMN_OH.
C
      OPEN(UNIT=NGENERICFILE,FILE='start_information',STATUS='OLD')
      REWIND(NGENERICFILE)
C
C*********************************************************************
C Loop over all NPARAM parameterizations.  There are over 200 individual
C parameterizations that are used to describe the tropospheric OH
C field.
C
C NPARAM  = total number of individual parameterizations.
C*********************************************************************
C
      DO M=1,NPARAM
C
C*********************************************************************
C The 24-hour average concentration of OH is represented as a set 
C  of high-order polynomials in variables such as temperature, latitude,
C  declination and the concentrations of ozone, water vapor, carbon
C  monoxide, nitrogen oxides (as a family), and hydrocarbons.  (See
C  SR GETINFO for a complete list of independent variables.)
C
C MOVAR = the total number of independent variables.
C         
      READ(NGENERICFILE,*) MOVAR(M)
C
C Loop over MOVAR independent variables reading in parameter ranges.
C These ranges for chemical species and physical parameters were estimated from
C observations and monthly average values predicted by a global model. The global
C model used for estimating ranges was the Harvard-GEOS model of transport and
C chemistry driven by assimilated meteorological data from the Goddard Earth
C Observing System Data Assimilation Office (GEOSDAO) (Bey et al. [1999]). 
C
C Bey, I., R. Yantosca and D. Jacob:  Export of Pollutants from Eastern Asia: A
C Simulation of the PEM-West (B) Aircraft Mission Using a 3-D Model Driven by
C Assimilated Meteorological Fields, Presentation at American Geophysical Union
C Spring Meeting, Boston, 1999.
C
C RANGEM = ranges for chemical species and physical parameters
C
        DO J=1,MOVAR(M)
C
          READ(NGENERICFILE,*)RANGEM(M,J,1),RANGEM(M,J,2)
C
        ENDDO
C
C IDENTOLD,ELTODO,IROWST,L = information of domain divisions. 
C
C NDONE = number of polynomials used in one parameterization. If there are
C         no domain divisions, NDONE=1.
C
      READ(NGENERICFILE,*) L
C
      DO K=1,L 
        READ(NGENERICFILE,*) IROWST(M,K,0),
     *       (IROWST(M,K,KR),KR=1,IROWST(M,K,0))
      ENDDO
C
        READ(NGENERICFILE,*) KK
C
      DO LH=1,KK
        READ(NGENERICFILE,*) (ELTODO(M,LH,J),J=1,MOVAR(M)+2)
      ENDDO
C
        READ(NGENERICFILE,*)NDONE(M),(IDENTOLD(M,J),J=1,NDONE(M))
C
C *********************************************
C Error Check.
       IF(NDONE(M).GT.MXDONE) THEN
         PRINT*,'SR READ_COEFF: NDONE > MXDONE'
         PRINT*,'NDONE(M)=',NDONE(M)
         PRINT*,'MXDONE  =',MXDONE
         PRINT*,'Increase MXDONE in CMN_OH to be greater than NDONE.'
         STOP
       ENDIF 
C *********************************************
C Read in number of terms and coefficients used in each polynomial.
C
C NENDW = number of terms in polynomial.
C
C COEFF = polynomial coefficients.
C
       DO IDONE=1,NDONE(M)
        READ(NGENERICFILE,*) NFT,NENDW(M,IDONE)
          DO JTERM=1,NENDW(M,IDONE)
            READ(NGENERICFILE,*) COEFF(M,IDONE,JTERM)
          ENDDO
       ENDDO 
C
C ACCUR = RMS error associated with individual parameterization.
C
      READ(NGENERICFILE,*) ACCUR(M)
C
C*********************************************************************
      ENDDO
C*********************************************************************
C
      CLOSE(NGENERICFILE)
C
C*********************************************************************
C
C Sum up the total number of individual parameterizations in PARSUM.  
C Since each parameterization may be described by more than one
C polynomial, the total number of polynomials are summed in the
C variables NDFUNCS*. NDFUNCS* are used only for bookkeeping. 
C See SR GETOH for a description of variables NPARAMA-F.
C
C*********************************************************************
        PARSUM=0
C*********************************************************************

        NDFUNCSB1=0
      DO J=1,NPARAMB
        NDFUNCSB1=NDFUNCSB1+NDONE(J)
      ENDDO

C*********************************************************************
       PARSUM=PARSUM+NPARAMB
C*********************************************************************
 
        NDFUNCSA1=0
      DO J=PARSUM+1,PARSUM+NPARAMA
        NDFUNCSA1=NDFUNCSA1+NDONE(J)       
      ENDDO

C*********************************************************************
       PARSUM=PARSUM+NPARAMA
C*********************************************************************

        NDFUNCSA2=0
      DO J=PARSUM+1,PARSUM+NPARAMA
        NDFUNCSA2=NDFUNCSA2+NDONE(J)       
      ENDDO

C*********************************************************************
       PARSUM=PARSUM+NPARAMA
C*********************************************************************

        NDFUNCSA3=0
      DO J=PARSUM+1,PARSUM+NPARAMA
        NDFUNCSA3=NDFUNCSA3+NDONE(J)       
      ENDDO

C*********************************************************************
       PARSUM=PARSUM+NPARAMA
C*********************************************************************

        NDFUNCSC1=0
      DO J=PARSUM+1,PARSUM+NPARAMC
        NDFUNCSC1=NDFUNCSC1+NDONE(J)       
      ENDDO

C*********************************************************************
       PARSUM=PARSUM+NPARAMC
C*********************************************************************

        NDFUNCSC2=0
      DO J=PARSUM+1,PARSUM+NPARAMC
        NDFUNCSC2=NDFUNCSC2+NDONE(J)       
      ENDDO

C*********************************************************************
       PARSUM=PARSUM+NPARAMC
C*********************************************************************

       NDFUNCSC3=0
      DO J=PARSUM+1,PARSUM+NPARAMC
        NDFUNCSC3=NDFUNCSC3+NDONE(J)
      ENDDO

C*********************************************************************
       PARSUM=PARSUM+NPARAMC
C*********************************************************************

        NDFUNCSD1=0
      DO J=PARSUM+1,PARSUM+NPARAMD
        NDFUNCSD1=NDFUNCSD1+NDONE(J)       
      ENDDO

C*********************************************************************
       PARSUM=PARSUM+NPARAMD
C*********************************************************************

        NDFUNCSE1=0
      DO J=PARSUM+1,PARSUM+NPARAME
        NDFUNCSE1=NDFUNCSE1+NDONE(J)       
      ENDDO

C*********************************************************************
       PARSUM=PARSUM+NPARAME
C*********************************************************************

c           goto 300
        NDFUNCSF1=0
      DO J=PARSUM+1,PARSUM+NPARAMF
        NDFUNCSF1=NDFUNCSF1+NDONE(J)       
      ENDDO

C*********************************************************************
       PARSUM=PARSUM+NPARAMF
C*********************************************************************

C Error Check.
       IF(PARSUM.NE.NPARAM) THEN
          PRINT*,'SR READ_COEFF: PARSUM.NE.NPARAM!'
          PRINT*,'This means that NPARAM specified in CMN_OH does'
          PRINT*,'not equal NPARAM (i.e., PARSUM here) calculated'
          PRINT*,'in this SR.  Check to see why they differ and'
          PRINT*,'correct the problem.'
          STOP
       ENDIF
C
C*********************************************************************
 300  RETURN 
      END
