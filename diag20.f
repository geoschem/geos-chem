! $Id: diag20.f,v 1.1 2003/06/30 20:26:05 bmy Exp $
      SUBROUTINE DIAG20
!
!******************************************************************************
!  Subroutine DIAG20 (bey, bmy, 6/9/99, 6/4/03) computes production and 
!  loss rates of O3, and then calls subroutine WRITE20 to save the these 
!  rates to disk. 
!
!  By saving the production and loss rates from a full-chemistry run,
!  a user can use these archived rates to perform a quick O3 chemistry
!  run at a later time.
!
!  DIAG20 assumes that ND65 (P-L diagnostics) have been turned on.
!
!  NOTES:
!  (1 ) Array PL24H is now declared dynamically allocatable in "diag_mod.f". 
!        This will save memory if DIAG20 is not called. (bmy, 1/5/00)
!  (2 ) FAMPL is now declared as allocatable in "diag_mod.f". (bmy, 3/16/00)
!  (3 ) Now references IDTOX from "tracerid_mod.f" (bmy, 11/6/02)
!  (4 ) Now make FIRSTDIAG20 a local variable.  Now use functions GET_NYMD,
!        GET_NHMS, GET_TAU. GET_TAUE, and GET_DAY from the new "time_mod.f".  
!        Remove FIRSTDIAG20 and NYMD from the arg list. (bmy, 2/11/03)
!  (5 ) Bug fix in format statement (bmy, 6/4/03)
!******************************************************************************
!
      ! References to F90 modules
      USE DIAG_MOD,     ONLY : FAMPL,    PL24H
      USE TIME_MOD,     ONLY : GET_NYMD, GET_TAU, GET_TAUE, GET_DAY
      USE TRACERID_MOD, ONLY : IDTOX

      IMPLICIT NONE

#     include "CMN_SIZE"     ! Size parameters
#     include "CMN"          ! STT
#     include "CMN_TIMES"    ! Timeseries parameters
#     include "CMN_O3"       ! XNUMOL

      ! Local variables
      LOGICAL, SAVE          :: FIRSTDIAG20 = .TRUE.
      INTEGER                :: I, J, L, N
      INTEGER, SAVE          :: DAY_LAST = -99
      INTEGER, SAVE          :: COUNT_IN_DAY
      INTEGER, SAVE          :: NYMD1
      REAL*8                 :: STT_TEMPO(IIPAR, JJPAR, LLPAR)
      REAL*8                 :: FPL_TEMPO
      REAL*8,  SAVE          :: YTAU1
      REAL*8                 :: YTAU2

      !=================================================================
      ! DIAG20 begins here!
      !=================================================================

      ! First-time initializatioin
      IF ( FIRSTDIAG20 ) THEN
         PRINT*, 'FIRSTDIAG20', firstdiag20
         DAY_LAST       = GET_DAY()
         COUNT_IN_DAY   = 0
         YTAU1          = GET_TAU()
         NYMD1          = GET_NYMD()
         FIRSTDIAG20    = .FALSE.

         ! Zero the PL24H array (it has been allocated in "ndxx_setup.f")
         PL24H(:,:,:,:) = 0d0
      ENDIF

      !=================================================================
      ! test if it is a new day.
      ! if it is a new day :
      ! * compute the average mean
      ! * write the average values in the newfile.
      ! * reset variables to zero
      ! * accumulate in the array forr the first time step of the new day
      !=================================================================
      IF ( GET_DAY() /= DAY_LAST .or. GET_TAU() == GET_TAUE()-1 ) THEN
         PL24H(:,:,:,:) = PL24H(:,:,:,:) / COUNT_IN_DAY
         YTAU2          = GET_TAU()

         ! Debug output...
         PRINT*, 'COUNT_IN_DAY : ',     COUNT_IN_DAY
         PRINT*, 'NYMD1, YTAU1, YTAU2', NYMD1, YTAU1, YTAU2

         CALL WRITE20( NYMD1, YTAU1, YTAU2, PL24H )

         ! zero arrays
         DAY_LAST       = GET_DAY()
         COUNT_IN_DAY   = 0
         PL24H(:,:,:,:) = 0d0

         ! accumulate for the first time step of this new day...
         COUNT_IN_DAY = COUNT_IN_DAY + 1
         IF ( COUNT_IN_DAY == 1 ) THEN 
            YTAU1 = GET_TAU()
            NYMD1 = GET_NYMD()
         ENDIF

         WRITE( 6, '(a, i6, f10.2, i4 )' ) 
     &        '--- DIAG20 - accumulation ', 
     &         GET_DAY(), GET_TAU(), COUNT_IN_DAY

         !==============================================================
         ! pl24h(:,:,:,1) is the production of Ox
         !                stored in kg/cm3/sec (fampl  = molec/cm3/sec 
         !                                     (xnumol = molec / kg of tracer 
         ! pl24h(:,:,:,2) is the loss rate of Ox
         !                stored in  /cm3/sec  (fampl     = molec/cm3/sec
         !
         !==============================================================
         STT_TEMPO(1:IIPAR, 1:JJPAR, 1:LLTROP) = 
     X        STT(1:IIPAR, 1:JJPAR, 1:LLTROP, IDTOX) * XNUMOL(IDTOX)

         DO L = 1, LLTROP
         DO J = 1, JJPAR
         DO I = 1, IIPAR

            ! production
            FPL_TEMPO      = FAMPL(I,J,L,1) / XNUMOL(IDTOX)
            PL24H(I,J,L,1) = PL24H(I,J,L,1) + FPL_TEMPO

            ! loss rate
            FPL_TEMPO      = FAMPL(I,J,L,2) / STT_TEMPO(I,J,L)
            PL24H(I,J,L,2) = PL24H(I,J,L,2) + FPL_TEMPO
            
         ENDDO
         ENDDO
         ENDDO

      ELSE

         !==============================================================
         ! This is not the beginning of a new day.  
         ! Therefore, just accumulate in the PL24H array.
         !==============================================================
         COUNT_IN_DAY = COUNT_IN_DAY + 1

         WRITE( 6, '(a, 2i6, f10.2, i4 )' ) 
     &        'DIAG20 - accumulation ', 
     &        GET_DAY(), GET_TAU(), COUNT_IN_DAY

         STT_TEMPO(1:IIPAR, 1:JJPAR, 1:LLTROP) =
     &        STT(1:IIPAR, 1:JJPAR, 1:LLTROP, IDTOX) * XNUMOL(IDTOX)
         
         DO L = 1, LLTROP
         DO J = 1, JJPAR
         DO I = 1, IIPAR

            ! production
            FPL_TEMPO      = FAMPL(I,J,L,1) / XNUMOL(IDTOX)
            PL24H(I,J,L,1) = PL24H(I,J,L,1) + FPL_TEMPO

            ! loss rate
            FPL_TEMPO      = FAMPL(I,J,L,2) / STT_TEMPO(I,J,L)
            PL24H(I,J,L,2) = PL24H(I,J,L,2) + FPL_TEMPO
         ENDDO
         ENDDO
         ENDDO

      ENDIF                     ! accumulation or new day...

      ! Return to calling program
      END SUBROUTINE DIAG20
