C $Id: OHparam.f,v 1.1 2003/06/30 20:26:04 bmy Exp $
C*****************************************************************************
C*****************************************************************************
C***************   OH PARAMETERIZATION CODE   ********************************
C*****************************************************************************
C*****************************************************************************
C Abstract.  We present a parameterization for the tropospheric 
C concentration of the hydroxyl radical (OH) which can be used
C to overcome the costs of solving kinetic equations in chemical 
C tracer models.  This parameterization accurately represents OH
C predicted by a full chemical mechanism.  The 24-hour average
C concentration of OH is represented as a set of high-order 
C polynomials in variables such as temperature, latitude,
C declination and the concentrations of ozone, water vapor, 
C carbon monoxide, nitrogen oxides (as a family), and
C hydrocarbons.  Results include computer-written FORTRAN 
C functions for an efficient computation of the polynomials.
C
C*****************************************************************************
C
C The OH parameterization was developed by:
C
C Bryan Duncan**
C email:  bnd@io.harvard.edu
C Department of Earth and Planetary Sciences and Division of 
C Engineering and Applied Sciences, Harvard University, Cambridge, MA
C
C **Contact to report bugs in the code or to ask questions.
C
C David Portman
C email:  David_Portman@aer.com
C Atmospheric & Environmental Research (AER), Inc., Cambridge, MA
C
C Clarissa Spivakovsky
C email:  cms@io.harvard.edu
C Department of Earth and Planetary Sciences and Division of 
C Engineering and Applied Sciences, Harvard University, Cambridge, MA
C
C*****************************************************************************
      SUBROUTINE OHparam(BOH,LTPAUSE,LMN)
C*****************************************************************************

      ! References to F90 modules
      USE ERROR_MOD, ONLY : IT_IS_NAN

      IMPLICIT NONE
#     include "CMN_OH"
      INTEGER I,J,L,LMN,NDODUST
      INTEGER,SAVE :: FIRSTDT
      DATA FIRSTDT /0/
C*****************************************************************************
C Created by Bryan Duncan.
C*****************************************************************************
C This SR is the driver of the code that calculates parameterized OH:
C
C 1) SR GETOH is called for each model grid box and the value
C    of the parameterized [OH] is placed in the array BOH.
C
C 2) Parameterized values of OH that are less than zero are
C    set to a small number.
C
C*****************************************************************************
C List of Variables & Arrays
C
C BOH = array holding tropospheric parameterized OH.
C
      REAL*8 BOH(MLONBX,MLATBX,MVRTBX)
C
C FIRSTDT = 0 on first chemistry time step.
C         = 1 on any time chemistry step, but first.
C
C Dimensions of troposphere.  User specified values in CMN_OH.
C   MVRTBX  = # of boxes in vertical
C   MLATBX  = # of boxes from pole to pole (following one longitude band)
C   MLONBX  = # of boxes circling globe (following one latitude band)
C
C   LTPAUSE = the vertical level of the tropopause.  Above this level,
C              no [OH] is calculated.  The user can feed this SR
C              a high value for LTPAUSE which effectively turns
C              this option off (i.e., LTPAUSE > MVRTBX).  Note that
C              this parameterization is only designed for pressures 
C              greater than 100 mb, however.  If the
C              [OH] = -9999 then the [OH] was not calculated.
C              LTPAUSE assumes a dimension of longitude and latitude.
C              The user may change this dimension to one- or zero-,
C              however, the if statement surrounding the call to
C              SR GETOH needs to be modified accordingly.
C
      INTEGER LTPAUSE(MLONBX,MLATBX)
C
C*****************************************************************************
C Only call SR READ_COEFF & SR readavgOH on initial time step!
C*****************************************************************************
C
      print*,'********************************************'
      print*,'Beginning OH param calculation in SR OHparam.'
      print*,'********************************************'
              CALL FLUSH( 6 )

      IF(FIRSTDT.EQ.0) THEN
C
C Read in start up information for OH parameterization.
C
         CALL READ_COEFF
C
C Read in climatological OH data.
C
         CALL readavgOH
C
C Read in ratio of OH with and without mineral dust.
C
C         CALL read_dustrat
C
C Read in ratio to correct aerosol optical depth problem.
C
         CALL read_correction
C
         FIRSTDT=1
C
      ENDIF ! IF (FIRSTDT)
C
C*****************************************************************************
C Call SR GETOH to calculate parameterized OH.
C*****************************************************************************
C
C The code is designed for use in a 3-D tropospheric model of transport 
C and chemistry. However, if the user wishes to use this code to
C calculate only one point, for instance, then specify MVRTBX,
C MLATBX & MLONBX (i.e., the dimensions specified by
C the user in CMN_OH) equal to 1, 1 & 1.  Then call SR OHPARAM
C for each point.
C
C Set BOH to large number for error check.
C
      BOH(:,:,:)=1.E10
C
C Loop over dimensions of troposphere one box at a time.
C
      DO L=1,MVRTBX
      DO J=1,MLATBX
      DO I=1,MLONBX
C
         IF ( L .LT. LTPAUSE(I,J) ) THEN 
C
            CALL GETOH(I,J,L,FIRSTDT,BOH)
              CALL FLUSH( 6 )
C
C SR AerosolOH: Adjust OH due to presence of absorbing/reflecting
C               mineral dust.
C***********************************************************
C Randall Martin has not completed his runs for all months.
C so PercentRed = 1 for now.
C***********************************************************
C NDODUST=0 Don't reduce OH.
C NDODUST=1 Do reduce OH.
C
C          NDODUST=0
C       IF(NDODUST.EQ.1) CALL AerosolOH(BOH,I,J,L,LMN)
C Error Check.
C
            IF(BOH(I,J,L).GE.1.E9) THEN
             PRINT*,''
             PRINT*,'OH not calculated for this box!'
             PRINT*,'  BOH(',I,',',J,',',L,') = ',BOH(I,J,L)
             PRINT*,'Stopped in SR OHPARAM'
             PRINT*,''
             STOP
            ENDIF

         IF ( IT_IS_NAN( INDVARA(I,J,L,5) ) ) THEN
           PRINT*,'Stopped in SR OHparam.'
           PRINT*, 'BOH is NaN!'
           PRINT*,I,J,L,BOH(I,J,L)
           STOP
         ENDIF
C
C Error Check. Warning message.
C
C            IF(BOH(I,J,L).LT.-100000.) THEN
C             PRINT*,''
C             PRINT*,'In SR OHparam.'
C             PRINT*,'Calculated OH for this box is really negative!'
C             PRINT*,'  BOH(',I,',',J,',',L,') = ',BOH(I,J,L)
C             PRINT*,''
C            ENDIF
C
C Where the parameterized OH is negative, set equal to small number.  
C Sometimes for low OH, the parameterized OH is actually negative.
C
            IF(BOH(I,J,L).LT.0.) BOH(I,J,L)=1000.
C
         ENDIF
C
      ENDDO  !DO I=1,MLONBX
      ENDDO  !DO J=1,MLATBX
      ENDDO  !DO L=1,MVRTBX

              CALL FLUSH( 6 )
C
C*****************************************************************************
C  Return to calling program
C*****************************************************************************
C
      print*,'********************************************'
      print*,'Leaving OH param calculation in SR OHparam.'
      print*,'********************************************'
              CALL FLUSH( 6 )
C
      RETURN
      END
