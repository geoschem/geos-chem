!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: diag_2pm
!
! !DESCRIPTION: Subroutine DIAG\_2PM constructs the diagnostic flag arrays:
!
! \begin{itemize}
! \item LTJV: J-values (ND22)
! \item LTOH: OH concentrations  (ND43)
! \item LTNO: NO concentrations  (ND43)
! \item LTNO2: NO2 concentrations (ND43)
! \item LTHO2: HO2 concentrations (ND43)
! \item LTOTH: used for tracers (ND45)
! \item LTO3: for O3 (ND45)
! \end{itemize}
!
!  These arrays are either 1 (if it is within a certain time interval)   
!  or 0 (if it is not within a certain time interval).  The limits of
!  the time intervals for CTOTH and CTJV are now defined in input.geos
!  The arrays CTOTH, CTOH, CTNO, CTJV count the number of times the 
!  diagnostics are accumulated for each grid box (i.e LTOTH is 1)
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE DIAG_2PM
!
! !USES:
!
      USE DIAG_MOD,       ONLY : LTJV,  CTJV,  LTNO,  CTNO,  CTO3
      USE DIAG_MOD,       ONLY : LTOH,  CTOH,  LTOTH, CTOTH, LTNO2
      USE DIAG_MOD,       ONLY : CTNO2, LTHO2, CTHO2, LTNO3, CTNO3
      USE DIAG_MOD,       ONLY : CTO3_24h, LTO3, CTO1D, LTO1D
      USE DIAG_MOD,       ONLY : LTLBRO2H, LTLBRO2N
      USE DIAG_MOD,       ONLY : LTLTRO2H, LTLTRO2N
      USE DIAG_MOD,       ONLY : LTLXRO2H, LTLXRO2N
      USE DIAG_MOD,       ONLY : CTLBRO2H, CTLBRO2N
      USE DIAG_MOD,       ONLY : CTLTRO2H, CTLTRO2N
      USE DIAG_MOD,       ONLY : CTLXRO2H, CTLXRO2N
      USE TIME_MOD,       ONLY : GET_LOCALTIME
      USE TIME_MOD,       ONLY : ITS_TIME_FOR_DIAG, ITS_TIME_FOR_CHEM
      USE CHEMGRID_MOD,   ONLY : ITS_IN_THE_CHEMGRID
      USE TIME_MOD,       ONLY : GET_ELAPSED_MIN
      USE TIME_MOD,       ONLY : GET_TS_DIAG

      USE CMN_SIZE_MOD   ! Size parameters
      USE CMN_DIAG_MOD   ! HR_OH1, HR_OH2, etc.

      IMPLICIT NONE
#     include "define.h"
!
! !REMARKS:
!  For now use GET_LOCALTIME( I, 1, 1 ) which will be independent of J and L
!  for a pure cartesian grid.  This may need to be revisited once G-C is
!  interfaced into a GCM.
!
! !REVISION HISTORY:
!  26 Mar 1999 - R. Yantosca - Initial version
!  (1 ) Now use F90 syntax (bmy, 3/26/99)
!  (2 ) Now reference LTNO2, CTNO2, LTHO2, CTHO2 arrays from "diag_mod.f".  
!        Updated comments, cosmetic changes.  (rvm, bmy, 2/27/02)
!  (3 ) Now removed NMIN from the arg list.  Now use functions GET_LOCALTIME,
!        ITS_TIME_FOR_CHEM, ITS_TIME_FOR_DYN from "time_mod.f" (bmy, 2/11/03)
!  (4 ) Now rewritten using a parallel DO-loop (bmy, 7/20/04)
!  (5 ) Now account for the time spent in the troposphere for ND43 and ND45
!        pure O3.  Now only accumulate counter for 3D pure O3 in ND45 if
!        it's a chemistry timestep. (phs, 1/24/07)
!  (6 ) Added 3D counter for ND65 and 03 in ND47 (phs, 11/17/08)
!  (7 ) Change re-initialization of ND45: only at the timestep after the 
!       diagnostics are accumulated. Add ITS_AFTER_DIAG and PREV_TS variables.
!       (ccc, 6/12/09)
!  (8 ) Add LTO3 to accumulate O3 in ND45 at the same place as the
!        chemistry (ccc, 7/17/09)
!  01 Mar 2012 - R. Yantosca - Now use GET_LOCALTIME(I,J,L) from time_mod.F90
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL :: IS_ND22, IS_ND43, IS_ND45, IS_ND45_O3
      LOGICAL :: IS_ND47, IS_ND65, ITS_AFTER_DIAG
      LOGICAL :: IS_CHEM
      INTEGER :: I,       J,       L
      INTEGER :: PREV_TS
      REAL*8  :: LT(IIPAR)

      !=================================================================
      ! DIAG_2PM begins here!
      !=================================================================

      ! Set logical flags
      IS_ND22    = ( ND22 > 0 .and. ITS_TIME_FOR_CHEM() )
      IS_ND45    = ( ND45 > 0                           )
      IS_ND43    = ( ND43 > 0 .and. ITS_TIME_FOR_CHEM() )
      IS_ND47    = ( ND47 > 0                           ) 
      IS_ND65    = ( ND65 > 0                           )
      IS_CHEM    = ( ITS_TIME_FOR_CHEM() )

      ! Calculate time at previous timestep for the diagnostics (ccc, 6/12/09)
      PREV_TS    = GET_ELAPSED_MIN()

      ! Indicate if we accumulated diagnostics at the previous timestep.
      ! (ccc, 6/12/09)
      ITS_AFTER_DIAG = ( MOD(PREV_TS, GET_TS_DIAG() ) == 0 )

      ! Pre-compute local time 
      DO I = 1, IIPAR
         LT(I) = GET_LOCALTIME( I, 1, 1 )
      ENDDO

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L )
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         !-----------------------------
         ! ND45 -- mixing ratios
         !-----------------------------
         IF ( IS_ND45 ) THEN

            ! Change initialization of LTOTH: only on timestep after saving
            ! diagnostics.(ccc, 6/12/09)
            IF ( ITS_AFTER_DIAG ) THEN
               LTOTH(I,J) = 0
            ENDIF

            ! Reinitialise LTO3 at each chemistry step
            IF ( IS_CHEM ) THEN
               LTO3(I,J)  = 0
            ENDIF

            ! Archive if we fall w/in the local time limits
            IF ( LT(I) >= HR1_OTH .and. LT(I) <= HR2_OTH ) THEN

               ! We only count the cells that are for the first time 
               ! between the time limits since the last accumulation of
               ! the diagnostics. (ccc, 7/17/09)
	       IF ( LTOTH(I, J) == 0 ) THEN
                  CTOTH(I,J) = CTOTH(I,J) + 1
	       ENDIF
               LTOTH(I,J) = 1
               
               ! Counter for # of O3 boxes in the chemistry grid
               ! (phs, 1/24/07)
               IF ( IS_CHEM ) THEN
                  LTO3(I,J) = 1
                  DO L = 1, LD45
                     IF ( ITS_IN_THE_CHEMGRID( I, J, L ) ) THEN
                        CTO3(I,J,L) = CTO3(I,J,L) + 1
                     ENDIF
                  ENDDO
               ENDIF
            ENDIF
         ENDIF

         !-----------------------------
         ! ND47_O3 / ND65 -- "chemistry all day long" counter
         !-----------------------------
         IF ( IS_ND47 .OR. IS_ND65 ) THEN

               ! Counter for # of O3 boxes in the chemgrid (phs, 11/17/08)
               DO L = 1, MAX( LD47, LD65 )
                 IF ( IS_CHEM.and.ITS_IN_THE_CHEMGRID( I, J, L )) THEN
                     CTO3_24h(I,J,L) = CTO3_24h(I,J,L) + 1
                  ENDIF
               ENDDO

         ENDIF

         !-----------------------------
         ! ND22 -- J-Value diagnostic
         !-----------------------------
         IF ( IS_ND22 ) THEN

            ! Archive if we fall w/in the local time limits
            IF ( LT(I) >= HR1_JV .and. LT(I) <= HR2_JV ) THEN
               LTJV(I,J) = 1
               CTJV(I,J) = CTJV(I,J) + 1
            ELSE
               LTJV(I,J) = 0
            ENDIF
         ENDIF

         !-----------------------------
         ! ND43 -- OH, NO, NO2, HO2
         !-----------------------------
         IF ( IS_ND43 ) THEN

            ! LTNO denotes where LT is between HR1_NO and HR2_NO
            ! CTNO counts the times when LT was between HR1_NO and HR2_NO
            ! Now set LTNO2, CTNO2 based on the NO times (rvm, bmy, 2/27/02)
            IF ( LT(I) >= HR1_NO .and. LT(I) <= HR2_NO ) THEN 
               LTNO(I,J)  = 1
               LTNO2(I,J) = 1
               
               ! update for arom (dkh, 06/22/07)  
               LTLBRO2N(I,J) = 1
               LTLTRO2N(I,J) = 1
               LTLXRO2N(I,J) = 1               

               ! Counters for # of NO, NO2 boxes in chemgrid (phs, 1/24/07)
               DO L = 1, LD43
                  IF ( ITS_IN_THE_CHEMGRID( I, J, L ) ) THEN
                     CTNO(I,J,L)  = CTNO(I,J,L)  + 1
                     CTNO2(I,J,L) = CTNO2(I,J,L) + 1
                     ! update for arom (dkh, 06/22/07)  
                     CTLBRO2N(I,J,L) = CTLBRO2N(I,J,L) + 1
                     CTLTRO2N(I,J,L) = CTLTRO2N(I,J,L) + 1
                     CTLXRO2N(I,J,L) = CTLXRO2N(I,J,L) + 1
                  ENDIF
               ENDDO

            ELSE
               LTNO(I,J)  = 0
               LTNO2(I,J) = 0
               ! update for arom (dkh, 06/22/07)  
               LTLBRO2N(I,J) = 0
               LTLTRO2N(I,J) = 0
               LTLXRO2N(I,J) = 0 
            ENDIF

            ! LTNO denotes where LT is between HR1_OH and HR2_OH
            ! CTNO counts the times when LT was between HR1_OH and HR2_OH
            ! Now set LTHO2, CTHO2 based on the OH times (rvm, bmy, 2/27/02)
            IF ( LT(I) >= HR1_OH .and. LT(I) <= HR2_OH ) THEN  
               LTOH(I,J)  = 1
               LTHO2(I,J) = 1
               LTNO3(I,J) = 1
               LTO1D(I,J) = 1
               ! update for arom (dkh, 06/22/07)  
               LTLBRO2H(I,J) = 1
               LTLTRO2H(I,J) = 1
               LTLXRO2H(I,J) = 1

               ! Counters for # of OH,HO2,NO3 boxes in chemgrid (phs, 1/24/07)
               DO L = 1, LD43 
                  IF ( ITS_IN_THE_CHEMGRID( I, J, L ) ) THEN
                     CTOH(I,J,L)  = CTOH(I,J,L)  + 1
                     CTHO2(I,J,L) = CTHO2(I,J,L) + 1
                     CTNO3(I,J,L) = CTNO3(I,J,L) + 1
                     CTO1D(I,J,L) = CTO1D(I,J,L) + 1
                     ! update for arom (dkh, 06/22/07)  
                     CTLBRO2H(I,J,L) = CTLBRO2H(I,J,L) + 1
                     CTLTRO2H(I,J,L) = CTLTRO2H(I,J,L) + 1
                     CTLXRO2H(I,J,L) = CTLXRO2H(I,J,L) + 1
                  ENDIF
               ENDDO

            ELSE
               LTOH(I,J)  = 0
               LTHO2(I,J) = 0
               LTNO3(I,J) = 0
               LTO1D(I,J) = 0
               ! update for arom (dkh, 06/22/07)  
               LTLBRO2H(I,J) = 0
               LTLTRO2H(I,J) = 0
               LTLXRO2H(I,J) = 0 
            ENDIF
         ENDIF

      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      END SUBROUTINE DIAG_2PM
!EOC
