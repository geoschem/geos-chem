! $Id: tropopause.f,v 1.5 2006/09/08 19:21:07 bmy Exp $
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
!******************************************************************************
!
      ! References to F90 modules.
      USE DAO_MOD,      ONLY : BXHEIGHT, T
      USE DIAG_MOD,     ONLY : AD55
      USE ERROR_MOD,    ONLY : ERROR_STOP
      USE GRID_MOD,     ONLY : GET_YMID
      USE PRESSURE_MOD, ONLY : GET_PCENTER

      IMPLICIT NONE

#     include "CMN_SIZE"  ! Size parameters
#     include "CMN"       ! LPAUSE
#     include "CMN_DIAG"  ! Diagnostic switches
      
      ! Local variables
      INTEGER :: I, J, L, K
      INTEGER :: LTPAUSE(IIPAR,JJPAR)
      REAL*8  :: HTPLIMIT, Y
      REAL*8  :: H(IIPAR,JJPAR,LLPAR)
      REAL*8  :: HTPAUSE(IIPAR,JJPAR)
      REAL*8  :: LAPSE_R(LLPAR)
      REAL*8  :: LAPSE_T(LLPAR)

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

      !=================================================================
      ! Initialize LTPAUSE and HTPAUSE arrays first.
      !
      ! LTPAUSE  = the tropopause layer #
      ! HTPAUSE  = the tropopause height [ km ] 
      ! HTPLIMIT = maximum tropopause height [ m ]
      !
      ! We need the factor of 1d3 to convert HTPLIMIT from km --> m.
      !=================================================================
      DO J = 1, JJPAR

         ! Latitude [degrees]
         Y = GET_YMID( J )

      DO I = 1, IIPAR
         LTPAUSE(I,J) = 0
         HTPAUSE(I,J) = 0.0d0

         ! In the tropics, the tropopause maxes out at 7 km, 
         ! elsewhere the tropopause maxes out at 5 km
         IF ( ABS( Y ) <= 30.d0 ) THEN
            HTPLIMIT = 7.0d3
         ELSE
            HTPLIMIT = 5.0d3
         ENDIF

         !==============================================================
         ! Tropopause: 15-20km (equator); 8-12km (temperate latitudes 
         ! and poles).  According to WMO, the "1st tropopause" is 
         ! defined as the lowest level at which the lapse rate decreases 
         ! to 2 C/km or less, provided also the average lapse rate 
         ! between this level and all higher levels within 2 km doesn't 
         ! exceed 2 C/km. It is noted that average lapse rate is the 
         ! difference between the temperatures at the respective end 
         ! points divided by the height interval irrespective of lapse 
         ! rate variations in the layer between the end points 
         ! (hyl, 03/28/99). 
         !
         ! NOTE: 2 C/km is equivalent to 2.0d-3 C/m.  We have to keep 
         !       this in C/m since H is in meters, and therefore, the 
         !       lapse rate will have units of K/m (= C/m).  
         !       (hyl, bmy, 11/30/99) 
         !==============================================================
         DO L = 2, LLPAR - 1
            IF ( H(I,J,L) >= HTPLIMIT ) THEN

               ! Lapse rate in the current level L
               LAPSE_R(L) = ( T(I,J,L  ) - T(I,J,L+1) ) / 
     &                      ( H(I,J,L+1) - H(I,J,L  ) )

               ! Test for lapse rate in the Lth level < 2 C/km
               IF ( LAPSE_R(L) <= 2.0d-3 ) THEN
                  K = 1

                  ! Compute the lapse rate at level L+K
 10               LAPSE_T(K) = ( T(I,J,L)   - T(I,J,L+K) ) / 
     &                         ( H(I,J,L+K) - H(I,J,L)   ) 

                  ! If the lapse rate at L+K is still less than 2 C/km, 
                  ! and level L+K is within 2 km of level K, we have found
                  ! the tropopause!  Go to the next surface box.
                  IF ( ( H(I,J,L+K) - H(I,J,L) ) > 2.0d3 .and. 
     &                 LAPSE_T(K) < 2.0d-3 ) THEN
                     LTPAUSE(I,J) = L
                     HTPAUSE(I,J) = H(I,J,L) / 1.0d3 ! m --> km
                     GOTO 30
                  ENDIF

                  ! If the lapse rate at L+K is greater than 2 C/km
                  ! then we are not high enough.  Go to the next level L.
                  IF ( ( H(I,J,L+K) - H(I,J,L) ) > 2.0d3 .and. 
     &                 LAPSE_T(K) > 2.0d-3 ) GOTO 20 
             
                  ! If level L+K is not within 2km of level K, then
                  ! we are not high enough.  Go to the next level L.
                  IF ( ( H(I,J,L+K) - H(I,J,L) ) < 2.0d3 .and. 
     &                 LAPSE_T(K) > 2.0d-3 ) GOTO 20 

                  ! Increment level K until the lapse rate at L+K < 2 C/km
                  K = K + 1

!-----------------------------------------------------------------------------
!          IF ( K>LM ) 
!             LAPSE_T(K) = (T(I,J,L)-T(I,J,L+K)) / (H(I,J,L+K)-H(I,J,L))
!          IF((H(I,J,L+K)-H(I,J,L)<2000 .and. LAPSE_T(K)>0.002) GOTO 20
!-----------------------------------------------------------------------------

                  ! If none of the above conditions were met, then
                  ! test here to make sure we don't violate array bounds
                  IF ( ( L + K ) <= LLPAR ) THEN

                     ! If Level L+K is within 2 km of level K, 
                     ! go back and compute the lapse rate
                     IF ( ( H(I,J,L+K) - H(I,J,L) ) < 2.0d3 ) THEN
                        GOTO 10

                     ! Otherwise we have found the tropopause level
                     ELSE  
                        LTPAUSE (I,J) = L
                        HTPAUSE (I,J) = H(I,J,L) / 1.0d3 ! m --> km
                        GOTO 30 
                     ENDIF
                  ELSE

                     ! If level L+K is higher than LM, then we are
                     ! going out of array bounds, so call L=LM
                     ! the tropopause level
                     LTPAUSE (I,J) = L
                     HTPAUSE (I,J) = H(I,J,L) / 1.0d3 ! m --> km
                     GOTO 30             
                  ENDIF
               ENDIF
            ENDIF
 20      ENDDO      ! L -- go to next level
 30   ENDDO         ! I -- go to next grid box
      ENDDO         ! J

      !=================================================================
      ! Sometimes a tropopause cannot be located in terms of the 
      ! above definition. For the time being, set it to that of  
      ! the nearest grid.
      !=================================================================
      DO J = 1, JJPAR
      DO I = 1, IIPAR
         IF ( LTPAUSE(I,J) == 0 ) THEN

            IF ( I /= 1 ) THEN
               LTPAUSE(I,J) = LTPAUSE(I-1,J)
               IF ( LTPAUSE(I,J) == 0 ) THEN
                  CALL ERROR_STOP( 'LTPAUSE = 0', 'tropopause.f' )
               ENDIF
               HTPAUSE(I,J) = H(I,J,LTPAUSE(I,J)) / 1.0d3 

            ELSE IF ( J /= 1 ) THEN
               LTPAUSE(I,J) = LTPAUSE(I,J-1)
               IF ( LTPAUSE(I,J) == 0 ) THEN
                  CALL ERROR_STOP( 'LTPAUSE = 0', 'tropopause.f' )
               ENDIF
               HTPAUSE(I,J) = H(I,J,LTPAUSE(I,J)) / 1.0d3

            ! South polar boxes
            ELSE IF ( J == 1 ) THEN

             ! Select the proper polar tropopause level 
             ! for GEOS-3 or GEOS-4 (bmy, 5/20/02)
#if defined( GEOS_3 )
               LTPAUSE(1,1) = 16

#elif defined( GEOS_4 )
               LTPAUSE(1,1) = 11

#endif

               HTPAUSE(1,1) = H(I,J,LTPAUSE(1,1)) / 1.0d3

            ENDIF 
         ENDIF
       
         !-------------------------------------------------------------------- 
         ! Debug output...check if LTPAUSE is 0 (hyl, 11/30/99)
         !IF ( LTPAUSE (I,J) == 0 ) THEN
         !   write(98,*) 'LTPAUSE(I,J)= ', LTPAUSE (I,J),I,J
         !   stop 'LTPAUSE = 0 in tropopause.f !'
         !ENDIF
         !--------------------------------------------------------------------

      ENDDO  ! I
      ENDDO  ! J

      !=================================================================
      ! ND55: Tropopause level, height [ km ], and pressure [ mb ]
      !       Recall that PW(I,J) = PS(I,J) - PTOP
      !=================================================================
      IF ( ND55 > 0 ) THEN
         DO J = 1, JJPAR
         DO I = 1, IIPAR
            L           = LTPAUSE(I,J)
            AD55(I,J,1) = AD55(I,J,1) + L
            AD55(I,J,2) = AD55(I,J,2) + HTPAUSE(I,J)
            AD55(I,J,3) = AD55(I,J,3) + GET_PCENTER(I,J,L)
         ENDDO
         ENDDO
      ENDIF

      ! Return to calling program
      END SUBROUTINE TROPOPAUSE
