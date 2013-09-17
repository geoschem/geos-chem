! $Id: pulsing.f,v 1.1 2009/09/16 14:06:14 bmy Exp $
#if defined( DEVEL )
      FUNCTION PULSING( I, J, M, LOCAL_MET ) RESULT( THE_PULSING )
#else
      FUNCTION PULSING( I, J, M ) RESULT( THE_PULSING )
#endif
!
!******************************************************************************
!  Function PULSING calculates the increase (or "pulse") of soil NO emission 
!  due to precipitation falling over a dry grid square and activating dormant 
!  (yhw, gmg, lwh, djj, 1994; bmy, 2/11/03)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) I (INTEGER)     : Grid box longitude index
!  (2 ) J (INTEGER)     : Grid box latitude index
!  (3 ) M (INTEGER)     : Grid box surface index (M=1,NLAND)
!
!  References:
!  ============================================================================
!  (1 ) Yienger, J.J, and H. Levy II, "Empirical model of global soil-biogenic
!        NOx emissions", JGR, 100 (D6), pp. 11447-11464, June 20, 1995.  See
!        section 4.1 of this work.
!
!  NOTES:
!  (1 ) Original code by by Yuhang Wang, Gerry Gardner and Prof. Daniel Jacob
!        written in the early 1990's.  Updated and modified for GEOS-CHEM by
!        Bob Yantosca.  Updated comments, cosmetic changes.  Now uses
!        function GET_TS_EMIS from "time_mod.f".  Removed NSRCE from the
!        arg list; this is now obsolete. (bmy, 2/11/03)
!******************************************************************************
!
      ! References to F90 modules
      USE TIME_MOD, ONLY : GET_TS_EMIS
#if defined( DEVEL )
      USE GC_TYPE_MOD, ONLY : GC_MET_LOCAL
#endif

      USE CMN_SIZE_MOD   
      USE COMMSOIL_MOD

      IMPLICIT NONE
#     include "define.h"
      
      ! Arguments
      INTEGER, INTENT(IN) :: I, J, M
#if defined( DEVEL )
      TYPE(GC_MET_LOCAL), INTENT(IN) :: LOCAL_MET  ! Obj w/ met fields
#endif

      ! Local variables
      INTEGER             :: K
      REAL*8              :: AREA, RATE, FRAC, EXPFACTOR, DTSRCE

      ! Function value
      REAL*8              :: THE_PULSING

      !=================================================================
      ! PULSING begins here!
      !=================================================================

      ! Emission timestep [days]
      DTSRCE = GET_TS_EMIS() / 1440d0

      !=================================================================
      ! SOILPULS(1,M) > 0 denotes dry soil.  Only dry 
      ! soil is subject to pulsing, so we proceed...
      !=================================================================
      IF ( SOILPULS(1,M) > 0.d0 ) THEN 

         ! Loop over pulse types (1=sprinkle, 2=shower, 3=heavy rain)
         DO K = 1, NPULSE

            ! SOILPULS(K+1,M) is the fraction of grid box M
            ! that is affected by fresh pulsing of type K
            IF ( SOILPULS(K+1,M) < 1.d-3 ) THEN

               ! No pulse assume evaporation
               SOILPULS(K+1,M) = 0.D0

            ELSE

               ! Pulse from previous time step decays exponentially 
               EXPFACTOR       = EXP( -PULSDECAY(K) * DTSRCE )
               SOILPULS(K+1,M) = SOILPULS(K+1,M) * EXPFACTOR

            ENDIF

         ENDDO

         !==============================================================
         ! Compute FRAC, the fraction of grid box (I,J) that is 
         ! undergoing precipitation.  Also compute RATE, the rate of 
         ! total precipitation at the ground (in mm/day).  RATE is 
         ! adjusted so that it only applies to the fraction of the 
         ! grid box where it is actually raining.
         !==============================================================
#if defined( DEVEL )
         CALL PRECIPFRAC( I, J, RATE, FRAC, LOCAL_MET )
#else
         CALL PRECIPFRAC( I, J, RATE, FRAC )
#endif

         !==============================================================
         ! We now determine if a new pulse is to be applied to the grid
         ! box due to precipitation over the current time step.  
         !
         ! The pulse is applied to the grid square fraction FRAC 
         ! experiencing precipitation.  We assume a characteristic 
         ! 1-day duration for precipitation in a given subgrid area of 
         ! the grid box, so that the full extent of pulsing (PULSFACT) 
         ! is realized over 24 hours.
         !
         ! For a model time step of NSRCE hours we reduce the pulsing 
         ! by a factor REAL(NSRCE)/24.
         !==============================================================
         IF ( ( RATE >= 1d0 ) .AND. ( RATE < 5d0 ) ) THEN

            ! Sprinkle
            SOILPULS(2,M) = SOILPULS(2,M) + ( FRAC * DTSRCE )

         ELSE IF ( ( RATE >= 5d0 ) .AND. ( RATE < 15d0 ) ) THEN

            ! K=3: Shower
            SOILPULS(3,M) = SOILPULS(3,M) + ( FRAC * DTSRCE ) 

         ELSE IF ( RATE >= 15d0 ) THEN

            ! K=4: Heavy rain
            SOILPULS(4,M) = SOILPULS(4,M) + ( FRAC * DTSRCE )

         ENDIF

         ! Initialize
         THE_PULSING = 0d0
         AREA        = 0d0

         !==============================================================
         ! Add up the contributions of the different pulses (K=1,3) to 
         ! obtain the total pulsing multiplicative factor PULSING; 
         ! PULSFACT is the multiplicative factor for fresh pulsing of 
         ! each type.  
         !
         ! Also determine the fractional grid box area AREA affected 
         ! by pulsing.  We assume that the area occupied by the 
         ! different pulses is additive, i.e., that successive pulses
         ! apply to different areas of the grid square and that the 
         ! area coccupied by a pulse decreases as the pulsing decays.  
         !
         ! If the resulting AREA is in excess of unity then the pulsing 
         ! must be scaled back to the grid box area.  If the AREA is 
         ! less than unity then we have to account for non-pulsing 
         ! emissions from the (1-AREA) non-pulsing fraction of the grid
         ! box.
         !==============================================================
         DO K = 1, NPULSE
            THE_PULSING = THE_PULSING + PULSFACT(K) * SOILPULS(1+K,M)
            AREA        = AREA        + SOILPULS(1+K,M)
         ENDDO

         IF ( AREA < 1d0 ) THEN
            THE_PULSING = THE_PULSING + 1d0 - AREA
         ELSE
            THE_PULSING = THE_PULSING / AREA

            DO K = 1, NPULSE
               SOILPULS(K+1,M) = SOILPULS(K+1,M) / AREA
            ENDDO
         ENDIF
       
      !=================================================================
      ! ...otherwise, the soil is wet, so no pulsing occurs
      !=================================================================  
      ELSE
         THE_PULSING = 1.D0

      ENDIF

      ! Return to calling program
      END FUNCTION PULSING
