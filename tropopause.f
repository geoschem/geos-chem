! $Id: tropopause.f,v 1.6 2006/10/16 20:44:36 phs Exp $
      SUBROUTINE TROPOPAUSE
!
!******************************************************************************
!  Subroutine TROPOPAUSE defines the tropopause layer in terms of temperature 
!  lapse rates. (hyl, bmy, 11/30/99, 8/4/06)
!
!  NOTES:
!  (1 ) Make sure the DO-loops go in the order L-J-I, wherever possible.
!  (2 ) Now archive ND55 diagnostic here rather than in DIAG1.F.  Also,
!        use an allocatable array (AD55) to archive tropopause heights.
!  (3 ) HTPAUSE is now a local variable, since it is only used here.
!  (4 ) Make LTPAUSE a local variable, since LPAUSE is used to store
!        the annual mean tropopause. (bmy, 4/17/00)
!  (5 ) Replace PW(I,J) with P(I,J).  Also updated comments. (bmy, 10/3/01)
!  (6 ) Removed obsolete code from 9/01 and 10/01 (bmy, 10/24/01)
!  (7 ) Added polar tropopause for GEOS-3 in #if defined( GEOS_3 ) block 
!        (bmy, 5/20/02) 
!  (8 ) Replaced all instances of IM with IIPAR and JM with JJPAR, in order
!        to prevent namespace confusion for the new TPCORE (bmy, 6/25/02)
!  (9 ) Now use GET_PCENTER from "pressure_mod.f" to compute the pressure
!        at the midpoint of box (I,J,L).  Also deleted obsolete, commented-out
!        code. (dsa, bdf, bmy, 8/21/02)
!  (10) Now reference BXHEIGHT and T from "dao_mod.f".  Also reference routine
!        ERROR_STOP from "error_mod.f" (bmy, 10/15/02)
!  (11) Now uses routine GET_YMID from "grid_mod.f" to compute grid box 
!        latitude. (bmy, 2/3/03)
!  (12) Add proper polar tropopause level for GEOS-4 (bmy, 6/18/03)
!  (13) Remove support for GEOS-1 and GEOS-STRAT met fields (bmy, 8/4/06)
!  (14) Get tropopause level from TROPOPAUSE_MOD.F routines (phs, 9/22/06)
!******************************************************************************
!
      ! References to F90 modules.
      USE DAO_MOD,        ONLY : BXHEIGHT  !, T
      USE DIAG_MOD,       ONLY : AD55
c      USE ERROR_MOD,      ONLY : ERROR_STOP
c      USE GRID_MOD,       ONLY : GET_YMID
      USE LOGICAL_MOD,    ONLY : LVARTROP
      USE PRESSURE_MOD,   ONLY : GET_PCENTER
      USE TROPOPAUSE_MOD, ONLY : GET_TPAUSE_LEVEL

      IMPLICIT NONE

#     include "CMN_SIZE"  ! Size parameters
#     include "CMN"       ! LPAUSE
#     include "CMN_DIAG"  ! Diagnostic switches
      
      ! Local variables
      INTEGER :: I, J, L !, K
c      INTEGER :: LTPAUSE(IIPAR,JJPAR)
c      REAL*8  :: HTPLIMIT, Y
      REAL*8  :: H(IIPAR,JJPAR,LLPAR)
c      REAL*8  :: HTPAUSE(IIPAR,JJPAR)
c      REAL*8  :: LAPSE_R(LLPAR)
c      REAL*8  :: LAPSE_T(LLPAR)

      !=================================================================
      ! TROPOPAUSE begins here! 
      !
      ! H (in m) is the height of the midpoint of layer L (hyl, 03/28/99) 
      !=================================================================

      ! Find height of the midpoint of the first level
      DO J = 1, JJPAR
      DO I = 1, IIPAR
         H(I,J,1) = BXHEIGHT(I,J,1) / 2.d0
      ENDDO
      ENDDO

      ! Add to H 1/2 of the sum of the two adjacent boxheights
      DO L = 1, LLPAR-1
      DO J = 1, JJPAR
      DO I = 1, IIPAR
         H(I,J,L+1) = H(I,J,L) + 
     &               ( BXHEIGHT(I,J,L) + BXHEIGHT(I,J,L+1) ) / 2.d0
      ENDDO
      ENDDO
      ENDDO


c  PRIOR 22/9/06
c
c      !=================================================================
c      ! Initialize LTPAUSE and HTPAUSE arrays first.
c      !
c      ! LTPAUSE  = the tropopause layer #
c      ! HTPAUSE  = the tropopause height [ km ] 
c      ! HTPLIMIT = maximum tropopause height [ m ]
c      !
c      ! We need the factor of 1d3 to convert HTPLIMIT from km --> m.
c      !=================================================================
c      DO J = 1, JJPAR
c
c         ! Latitude [degrees]
c         Y = GET_YMID( J )
c
c      DO I = 1, IIPAR
c         LTPAUSE(I,J) = 0
c         HTPAUSE(I,J) = 0.0d0
c
c         ! In the tropics, the tropopause maxes out at 7 km, 
c         ! elsewhere the tropopause maxes out at 5 km
c         IF ( ABS( Y ) <= 30.d0 ) THEN
c            HTPLIMIT = 7.0d3
c         ELSE
c            HTPLIMIT = 5.0d3
c         ENDIF
c
c         !==============================================================
c         ! Tropopause: 15-20km (equator); 8-12km (temperate latitudes 
c         ! and poles).  According to WMO, the "1st tropopause" is 
c         ! defined as the lowest level at which the lapse rate decreases 
c         ! to 2 C/km or less, provided also the average lapse rate 
c         ! between this level and all higher levels within 2 km doesn't 
c         ! exceed 2 C/km. It is noted that average lapse rate is the 
c         ! difference between the temperatures at the respective end 
c         ! points divided by the height interval irrespective of lapse 
c         ! rate variations in the layer between the end points 
c         ! (hyl, 03/28/99). 
c         !
c         ! NOTE: 2 C/km is equivalent to 2.0d-3 C/m.  We have to keep 
c         !       this in C/m since H is in meters, and therefore, the 
c         !       lapse rate will have units of K/m (= C/m).  
c         !       (hyl, bmy, 11/30/99) 
c         !==============================================================
c         DO L = 2, LLPAR - 1
c            IF ( H(I,J,L) >= HTPLIMIT ) THEN
c
c               ! Lapse rate in the current level L
c               LAPSE_R(L) = ( T(I,J,L  ) - T(I,J,L+1) ) / 
c     &                      ( H(I,J,L+1) - H(I,J,L  ) )
c
c               ! Test for lapse rate in the Lth level < 2 C/km
c               IF ( LAPSE_R(L) <= 2.0d-3 ) THEN
c                  K = 1
c
c                  ! Compute the lapse rate at level L+K
c 10               LAPSE_T(K) = ( T(I,J,L)   - T(I,J,L+K) ) / 
c     &                         ( H(I,J,L+K) - H(I,J,L)   ) 
c
c                  ! If the lapse rate at L+K is still less than 2 C/km, 
c                  ! and level L+K is within 2 km of level K, we have found
c                  ! the tropopause!  Go to the next surface box.
c                  IF ( ( H(I,J,L+K) - H(I,J,L) ) > 2.0d3 .and. 
c     &                 LAPSE_T(K) < 2.0d-3 ) THEN
c                     LTPAUSE(I,J) = L
c                     HTPAUSE(I,J) = H(I,J,L) / 1.0d3 ! m --> km
c                     GOTO 30
c                  ENDIF
c
c                  ! If the lapse rate at L+K is greater than 2 C/km
c                  ! then we are not high enough.  Go to the next level L.
c                  IF ( ( H(I,J,L+K) - H(I,J,L) ) > 2.0d3 .and. 
c     &                 LAPSE_T(K) > 2.0d-3 ) GOTO 20 
c             
c                  ! If level L+K is not within 2km of level K, then
c                  ! we are not high enough.  Go to the next level L.
c                  IF ( ( H(I,J,L+K) - H(I,J,L) ) < 2.0d3 .and. 
c     &                 LAPSE_T(K) > 2.0d-3 ) GOTO 20 
c
c                  ! Increment level K until the lapse rate at L+K < 2 C/km
c                  K = K + 1
c
c!-----------------------------------------------------------------------------
c!          IF ( K>LM ) 
c!             LAPSE_T(K) = (T(I,J,L)-T(I,J,L+K)) / (H(I,J,L+K)-H(I,J,L))
c!          IF((H(I,J,L+K)-H(I,J,L)<2000 .and. LAPSE_T(K)>0.002) GOTO 20
c!-----------------------------------------------------------------------------
c
c                  ! If none of the above conditions were met, then
c                  ! test here to make sure we don't violate array bounds
c                  IF ( ( L + K ) <= LLPAR ) THEN
c
c                     ! If Level L+K is within 2 km of level K, 
c                     ! go back and compute the lapse rate
c                     IF ( ( H(I,J,L+K) - H(I,J,L) ) < 2.0d3 ) THEN
c                        GOTO 10
c
c                     ! Otherwise we have found the tropopause level
c                     ELSE  
c                        LTPAUSE (I,J) = L
c                        HTPAUSE (I,J) = H(I,J,L) / 1.0d3 ! m --> km
c                        GOTO 30 
c                     ENDIF
c                  ELSE
c
c                     ! If level L+K is higher than LM, then we are
c                     ! going out of array bounds, so call L=LM
c                     ! the tropopause level
c                     LTPAUSE (I,J) = L
c                     HTPAUSE (I,J) = H(I,J,L) / 1.0d3 ! m --> km
c                     GOTO 30             
c                  ENDIF
c               ENDIF
c            ENDIF
c 20      ENDDO      ! L -- go to next level
c 30   ENDDO         ! I -- go to next grid box
c      ENDDO         ! J
c
c      !=================================================================
c      ! Sometimes a tropopause cannot be located in terms of the 
c      ! above definition. For the time being, set it to that of  
c      ! the nearest grid.
c      !=================================================================
c      DO J = 1, JJPAR
c      DO I = 1, IIPAR
c         IF ( LTPAUSE(I,J) == 0 ) THEN
c
c            IF ( I /= 1 ) THEN
c               LTPAUSE(I,J) = LTPAUSE(I-1,J)
c               IF ( LTPAUSE(I,J) == 0 ) THEN
c                  CALL ERROR_STOP( 'LTPAUSE = 0', 'tropopause.f' )
c               ENDIF
c               HTPAUSE(I,J) = H(I,J,LTPAUSE(I,J)) / 1.0d3 
c
c            ELSE IF ( J /= 1 ) THEN
c               LTPAUSE(I,J) = LTPAUSE(I,J-1)
c               IF ( LTPAUSE(I,J) == 0 ) THEN
c                  CALL ERROR_STOP( 'LTPAUSE = 0', 'tropopause.f' )
c               ENDIF
c               HTPAUSE(I,J) = H(I,J,LTPAUSE(I,J)) / 1.0d3
c
c            ! South polar boxes
c            ELSE IF ( J == 1 ) THEN
c
c             ! Select the proper polar tropopause level 
c             ! for GEOS-3 or GEOS-4 (bmy, 5/20/02)
c#if defined( GEOS_3 )
c               LTPAUSE(1,1) = 16
c
c#elif defined( GEOS_4 )
c               LTPAUSE(1,1) = 11
c
c#endif
c
c               HTPAUSE(1,1) = H(I,J,LTPAUSE(1,1)) / 1.0d3
c
c            ENDIF 
c         ENDIF
c       
c         !-------------------------------------------------------------------- 
c         ! Debug output...check if LTPAUSE is 0 (hyl, 11/30/99)
c         !IF ( LTPAUSE (I,J) == 0 ) THEN
c         !   write(98,*) 'LTPAUSE(I,J)= ', LTPAUSE (I,J),I,J
c         !   stop 'LTPAUSE = 0 in tropopause.f !'
c         !ENDIF
c         !--------------------------------------------------------------------
c
c      ENDDO  ! I
c      ENDDO  ! J
c
c************************ END of PRIOR 9/22/06 *************************





      !=================================================================
      ! ND55: Tropopause level, height [ km ], and pressure [ mb ]
      !       Recall that PW(I,J) = PS(I,J) - PTOP
      !=================================================================
      IF ( ND55 > 0 ) THEN
         DO J = 1, JJPAR
         DO I = 1, IIPAR
c prior 9/22/06:
c            L           = LTPAUSE(I,J)
            L           = GET_TPAUSE_LEVEL( I, J )
            IF ( LVARTROP ) L = L+1
            AD55(I,J,1) = AD55(I,J,1) + L
c prior 9/22/06:
c            AD55(I,J,2) = AD55(I,J,2) + HTPAUSE(I,J)
            AD55(I,J,2) = AD55(I,J,2) + H(I,J,L) / 1.0d3 ! m --> km
            AD55(I,J,3) = AD55(I,J,3) + GET_PCENTER(I,J,L)
         ENDDO
         ENDDO
      ENDIF

      ! Return to calling program
      END SUBROUTINE TROPOPAUSE
