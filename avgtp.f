! $Id: avgtp.f,v 1.1 2003/06/30 20:26:05 bmy Exp $
      SUBROUTINE avgtp(NTDT,NMIN)
C*****************************************************************************

      ! References to F90 modules
      USE PRESSURE_MOD, ONLY : GET_PCENTER

      IMPLICIT NONE

#     include "CMN_SIZE"
#     include "CMN"
#     include "CMN_CO"

      INTEGER L,NTIMES,MNDT,NTDT,NMIN,K,M,N
      REAL*8 Ptemp(IIPAR,JJPAR,LLPAR)
C
C Designed to calculate a 24-hour average temperature and
C   pressure to be used for the CO/OH parameterization
C   option (see SR chemco).
C Created by Bryan Duncan (1/99).
C
C NCHEM is chemical time step (i.e., 24 hours).
C NDT   is meteorological data time step (i.e., 6 hours).
C NTDT  is dynamical time step (i.e., 1/2 hours).
C
C NOTE: Now use GET_PCENTER from "pressure_mod.f" to compute Pressure at 
C       midpoint of grid box (I,J,L) for hybrid grid.  Eliminate PS as
C       an argument. (dsa, bdf, bmy, 8/21/02)
C
      MNDT=NTDT/60 
      NTIMES=NCHEM/MNDT
C
      IF(NMIN.LE.NCHEM) NTIMES=NTIMES+1
C
C ****************************************************************************
      IF (NMIN.EQ.0) THEN
C ****************************************************************************
	 NNEW=0
C
         IF (NCHEM.NE.1440) THEN   
            WRITE(*,*) ' '
            WRITE(*,*) 'CO-OH parameterization option (i.e., NSRCX=5)!' 
            WRITE(*,*) 'Use a chemistry time step = 24 hours'
            WRITE(*,*) '(i.e., NCHEM=1440 min.)'
            WRITE(*,*) ' '
            STOP
         ENDIF
C
         IF (mod(NCHEM,MNDT).NE.0) THEN   
            WRITE(*,*) ' '
            WRITE(*,*) 'CO-OH parameterization option (i.e., NSRCX=5)!'
            WRITE(*,*) 'The chemistry time step (i.e., 24 hours) is'
            WRITE(*,*) 'not evenly divisible by the meteorological'
            WRITE(*,*) 'data read-in time step (i.e., 6 hours).  This'
            WRITE(*,*) 'will mess up SR avgtp which calculates a 24-'
            WRITE(*,*) 'hour average temperature and pressure to be'
            WRITE(*,*) 'used by SR getinfo.'
            WRITE(*,*) ' '
            STOP
         ENDIF
C
C If NCHEM < NTDT then stop program.
C
         IF (NCHEM.LT.MNDT) THEN   
            WRITE(*,*) ' '
            WRITE(*,*) 'When using the CO-OH parameterization'
            WRITE(*,*) 'option (i.e., NSRCX=5), take a 24-hour'
            WRITE(*,*) 'time step (i.e., NCHEM=1440 min.) because'
            WRITE(*,*) 'the OH parameterization produces a 24-hour'
            WRITE(*,*) 'average [OH]'
            WRITE(*,*) ' '
            STOP
         ENDIF
C ****************************************************************************
      ENDIF
C ****************************************************************************
C
C If a new 24-hr period, set Pavg = 0.
C
      IF (NNEW.EQ.0) THEN 
         Pavg(:,:,:)=0.
         Tavg(:,:,:)=0.
	 NNEW=1
	 NNCOUNT=0
      ENDIF
C
C Convert the surface pressure (PS) to the pressure in
C  each box.  Save value in Pavg().
C

      ! Add parallel do-loop (dsa, bdf, bmy, 8/21/02)
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L )
      DO L = 1, LLPAR
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! Sum pressure at sigma-center of box (I,J,L) in Pavg
         Pavg(I,J,L) = Pavg(I,J,L) + GET_PCENTER(I,J,L)

         ! Sum temperature in Tavg
         Tavg(I,J,L) = Tavg(I,J,L) + T(I,J,L)
      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

C      
C Keep track to see if at end of NCHEM time step.
C
      NNCOUNT=NNCOUNT+1
C
      IF (NNCOUNT.EQ.NTIMES) THEN
C
         Pavg(:,:,1:LLPAR)=Pavg(:,:,1:LLPAR)/REAL(NTIMES)
         Tavg(:,:,1:LLPAR)=Tavg(:,:,1:LLPAR)/REAL(NTIMES)
C     
         NNEW=0
C     
      ENDIF
C
      RETURN
      END
