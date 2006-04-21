! $Id: ruralbox.f,v 1.5 2006/04/21 15:40:07 bmy Exp $
!-----------------------------------------------------------------------------
!      SUBROUTINE RURALBOX( AD,     T,      AVGW,   ALT,   ALBD,  
!     &                     SUNCOS, CLOUDS, LEMBED, IEBD1, IEBD2,  
!     &                     JEBD1,  JEBD2 )
!-----------------------------------------------------------------------------
      SUBROUTINE RURALBOX( AD,     T,     AVGW,  ALBD,  SUNCOS, 
     &                     LEMBED, IEBD1, IEBD2, JEBD1, JEBD2 )
!
!******************************************************************************
!  Subroutine RURALBOX computes which boxes are tropospheric and which
!  are stratospheric.  SMVGEAR arrays are initialized with quantities from 
!  tropospheric boxes.  (amf, bey, ljm, lwh, gmg, bdf, bmy, 7/16/01, 4/10/06)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) AD     (REAL*8 ) : Array for dry air mass               [   kg   ]
!  (2 ) T      (REAL*8 ) : Array for grid box temperatures      [   K    ]
!  (3 ) AVGW   (REAL*8 ) : Array for mixing ratio of water      [  v/v   ]
!  (4 ) ALBD   (REAL*8 ) : Array for visible albedo             [unitless]
!  (5 ) SUNCOS (REAL*8 ) : Array for COS( Solar Zenith Angle )  [unitless]
!  (6 ) LEMBED (LOGICAL) : Switch for embedded chemistry region [ T or F ]
!  (7 ) IEBD1  (INTEGER) : Lon: lower right corner }  of the    [unitless] 
!  (8 ) IEBD2  (INTEGER) : Lon: upper left  corner } embedded   [unitless]
!  (9 ) JEBD1  (INTEGER) : Lat: lower right corner } chemistry  [unitless]
!  (10) JEBD2  (INTEGER) : Lat: upper left  corner }  region    [unitless]
!
!  NOTES:
!  (1 ) Remove PTOP from the arg list.  PTOP is now a parameter
!        in "CMN_SIZE". (bmy, 2/10/00)
!  (2 ) Add C-preprocessor switch LSLOWJ to bracket code for 
!        SLOW-J photolysis (bmy, 2/25/00)
!  (3 ) Now reference ABHSUM, AIRDENS, IXSAVE, IYSAVE, IZSAVE, JLOP, PRESS3, 
!        T3, and VOLUME from F90 module "comode_mod.f" (bmy, 10/19/00)
!  (4 ) PTOP is already a parameter in "CMN_SIZE", don't declare it here
!        (bmy, 7/16/01)
!  (5 ) Replace IGCMPAR,JGCMPAR,LGCMPAR with IIPAR,JJPAR,LLPAR.  Also moved
!        CLOUDREF to SLOW-J block.  Also remove IREF, JREF, IOFF, JOFF, these
!        are now obsolete.  Updated comments. (bmy, 9/25/01)
!  (6 ) Eliminate I00 and J00 as arguments, these are obsolete (bmy, 9/28/01)
!  (7 ) Removed obsolete, commented out code from 9/01 (bmy, 10/24/01)
!  (8 ) Updated comment header.  Also updated comments, and made cosmetic 
!        changes. (bmy, 4/15/02)
!  (9 ) Bug fix: declare variables for SLOW-J photolysis.  Also eliminated
!        obsolete code from 4/15/02. (bmy, 8/5/02)
!  (10) Now reference GET_PCENTER and GET_PEDGE from "pressure_mod.f", 
!        which return the correct "floating" pressure.  Also deleted obsolete,
!        commented-out code.  Also eliminate P, SIG, and NSKIPL from the arg 
!        list, since we don't need them anymore. (dsa, bdf, bmy, 8/20/02)
!  (11) Added modifications for SMVGEAR II (gcc, bdf, bmy, 4/1/03)
!  (12) SLOW-J is now obsolete; remove LSLOWJ #ifdef blocks (bmy, 6/23/05)
!  (13) Now reference ITS_IN_THE_TROP and ITS_IN_THE_STRAT from 
!        "tropopause_mod.f" to diagnose trop & strat boxes.  Also remove
!        LPAUSE from the arg list (bmy, 8/22/05)
!  (14) Remove ALT and CLOUDS from arg list -- they are obsolete (bmy, 4/10/06)
!******************************************************************************
!
      ! References to F90 modules
      USE COMODE_MOD,     ONLY : ABSHUM, AIRDENS, IXSAVE, IYSAVE, 
     &                           IZSAVE, JLOP,    PRESS3, T3,  VOLUME 
      USE PRESSURE_MOD,   ONLY : GET_PCENTER, GET_PEDGE
      USE TROPOPAUSE_MOD, ONLY : ITS_IN_THE_STRAT, ITS_IN_THE_TROP

      IMPLICIT NONE

#     include "CMN_SIZE"       ! Size parameters
#     include "comode.h"       ! NPVERT

      LOGICAL, INTENT(IN)     :: LEMBED 
      INTEGER, INTENT(IN)     :: IEBD1, IEBD2, JEBD1, JEBD2
      REAL*8,  INTENT(IN)     :: AD(IIPAR,JJPAR,LLPAR)
      REAL*8,  INTENT(IN)     :: T(IIPAR,JJPAR,LLPAR)
      REAL*8,  INTENT(IN)     :: AVGW(IIPAR,JJPAR,LLPAR)
      !---------------------------------------------------
      ! Prior to 4/10/06:
      ! ALT is now obsolete (bmy, 4/10/06)
      !REAL*8,  INTENT(IN)     :: ALT(MAXIJ,NPVERT)
      !---------------------------------------------------
      REAL*8,  INTENT(IN)     :: ALBD(IIPAR,JJPAR)
      REAL*8,  INTENT(IN)     :: SUNCOS(MAXIJ)
      !---------------------------------------------------
      ! Prior to 4/10/06:
      ! CLOUDS is now obsolete (bmy, 4/10/06)
      !REAL*8,  INTENT(INOUT)  :: CLOUDS(MAXIJ,11)
      !---------------------------------------------------

      ! Local variables      
      LOGICAL                 :: LDEBUG
      INTEGER                 :: I, J, L, JLOOP, IJLOOP, LL

      ! External functions
      REAL*8,  EXTERNAL       :: BOXVL

      !=================================================================
      ! RURALBOX begins here!
      !=================================================================
      LDEBUG = .FALSE.

      ! Rural Boxes
      JLOOP     = 0
      NTLOOPNCS = 0

      ! Loop over vertical levels (max = LLTROP) 
      DO L = 1, NVERT

         ! Loop over surface grid boxes
         DO J = 1, NLAT
         DO I = 1, NLONG

            ! JLOP is the 1-D grid box loop index
            JLOP(I,J,L) = 0

            ! Filter to do chemistry in a window when 
            ! rest of model is running global run.
            ! LEMBED - Logical for embedded window defined by
            !          IEBD1, IEBD2, JEBD1, JEBD2
            IF ( LEMBED ) THEN
               IF ( I < IEBD1 .OR. I > IEBD2  .OR. 
     &              J < JEBD1 .OR. J > JEBD2 ) GOTO 40
            ENDIF

            IF ( IGLOBCHEM <= 0 ) THEN

               !=======================================================
               ! Skip over strat boxes
               !=======================================================
               IF ( ITS_IN_THE_STRAT( I, J, L ) ) GOTO 40

               ! Increment JLOOP for trop boxes
               JLOOP          = JLOOP + 1
               JLOP(I,J,L)    = JLOOP

            ELSE

               !=======================================================
               ! If we're doing a trop/strat run, IGLOBCHEM > 0.
               ! In that case we have to tell SMVGEAR which boxes are 
               ! tropospheric and which are stratospheric.  We do this 
               ! using NTLOOPNCS and NCSLOOP. (gcc, bdf, bmy, 4/1/03)
               !
               ! NTLOOPNCS counts the # of urban, trop, strat boxes
               ! NCSLOOP   holds the 1-D grid box indices for 
               !
               ! NOTE: L <  LPAUSE(I,J) are tropospheric boxes  
               !       L >= LPAUSE(I,J) are stratospheric boxes 
               !========================================================
               
               ! Increment JLOOP for all boxes
               JLOOP          = JLOOP + 1
               JLOP(I,J,L)    = JLOOP

               IF ( ITS_IN_THE_TROP( I, J, L ) ) THEN

                  ! Tropospheric boxes go into the SMVGEAR II "URBAN" slot
                  NTLOOPNCS(NCSURBAN) = NTLOOPNCS(NCSURBAN) + 1
                  NCSLOOP(NTLOOPNCS(NCSURBAN),NCSURBAN) = JLOOP

               !-----------------------------------------------------------
               ! Comment this out for now -- restore it later (bmy, 4/21/03)
               !ELSE IF ( .FALSE. ) THEN
               !
               !   ! The SMVGEAR II "FREE TROPOSPHERE" slot is unused
               !   NTLOOPNCS(NCSTROP) = NTLOOPNCS(NCSTROP) + 1
               !   NCSLOOP(NTLOOPNCS(NCSTROP),NCSTROP) = JLOOP
               !-----------------------------------------------------------

               ELSE

                  ! Stratospheric boxes go into the SMVGEAR II "STRAT" slot
                  ! (for now GEOS-CHEM skips these; later we will define
                  !  a stratospheric chemistry mechanism a la G. Curci).
                  NTLOOPNCS(NCSSTRAT) = NTLOOPNCS(NCSSTRAT) + 1
                  NCSLOOP(NTLOOPNCS(NCSSTRAT),NCSSTRAT) = JLOOP

               ENDIF
               
            ENDIF

            ! These translate JLOOP back to an (I,J,L) triplet
            IXSAVE(JLOOP)  = I
            IYSAVE(JLOOP)  = J
            IZSAVE(JLOOP)  = L                              

            ! get box volume [cm3]
            VOLUME(JLOOP)  = BOXVL(I, J, L)

            ! get air density in (molecs cm^-3)
            AIRDENS(JLOOP) = AD(I,J,L)*1000.d0/VOLUME(JLOOP)*AVG/WTAIR

            ! get temperature
            T3(JLOOP)      = T(I,J,L)

            ! PRESS3 = pressure in bar, multiply mb * 1000
            PRESS3(JLOOP) = GET_PCENTER(I,J,L) * 1000d0

            ! Get relative humidity (here is absolute #H2O/cc air)
            ! AVGW is the mixing ratio of water vapor [v/v]
            ABSHUM(JLOOP) = AVGW(I,J,L) * AIRDENS(JLOOP)

            ! Go to next I
 40         CONTINUE           
         ENDDO
         ENDDO

         ! NIJLOOP is the number of surface boxes
         IF ( L == 1 ) NIJLOOP = JLOOP
      ENDDO

      ! NTLOOP is the number of total tropospheric boxes
      NTLOOP = JLOOP

      ! Return to calling program
      END SUBROUTINE RURALBOX
