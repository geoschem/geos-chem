!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: ruralbox
!
! !DESCRIPTION: Subroutine RURALBOX computes which boxes are tropospheric 
!  and which are stratospheric.  SMVGEAR arrays are initialized with 
!  quantities from tropospheric boxes.  
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE RURALBOX( State_Met )
!
! !USES:
!
      USE COMODE_MOD,         ONLY : ABSHUM, AIRDENS
      USE COMODE_MOD,         ONLY : IXSAVE, IYSAVE
      USE COMODE_MOD,         ONLY : IZSAVE, JLOP
      USE COMODE_MOD,         ONLY : PRESS3, T3,       VOLUME 
      USE GIGC_State_Met_Mod, ONLY : MetState
      USE PRESSURE_MOD,       ONLY : GET_PCENTER,      GET_PEDGE
      USE TROPOPAUSE_MOD,     ONLY : ITS_IN_THE_STRAT, ITS_IN_THE_TROP
      USE CMN_SIZE_MOD             ! Size parameters
      USE COMODE_LOOP_MOD          ! NPVERT

      IMPLICIT NONE

!
! !INPUT PARAMETERS: 
!
      TYPE(MetState), INTENT(IN)  :: State_Met   ! Meteorology State object
!
! !REMARKS:
!  Developers: amf, bey, ljm, lwh, gmg, bdf, bmy, 7/16/01, 2/25/10)
! 
! !REVISION HISTORY: 
!  01 Oct 1995 - M. Prather  - Initial version
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
!  (15) Remove obsolete embedded chemistry stuff (bmy, 2/25/10)
!  10 Sep 2010 - R. Yantosca - Added ProTeX headers
!  09 Nov 2012 - M. Payer    - Replaced all met field arrays with State_Met
!                              derived type object
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL                 :: LDEBUG
      INTEGER                 :: I, J, L, JLOOP, IJLOOP
#if defined( DEVEL )
      INTEGER                 :: TMP_LEV
#endif
      ! External functions
      REAL*8,  EXTERNAL       :: BOXVL

      !! testing variables for variable tropopause
      !real*8  :: temp(iipar,jjpar,llpar)
      !logical :: templ
      !temp = 0

      !=================================================================
      ! RURALBOX begins here!
      !=================================================================
      LDEBUG = .FALSE.

      ! Rural Boxes
      JLOOP     = 0
      NTLOOPNCS = 0

      ! Loop over vertical levels (max = LLTROP) 
#if defined( DEVEL )
      IF (NVERT .gt. LLPAR) THEN
         TMP_LEV = LLPAR
      ELSE
         TMP_LEV = NVERT
      ENDIF
      DO L = 1, TMP_LEV
#else
      DO L = 1, NVERT
#endif

         ! Loop over surface grid boxes
         DO J = 1, NLAT
         DO I = 1, NLONG

            ! JLOP is the 1-D grid box loop index
            JLOP(I,J,L) = 0

            IF ( IGLOBCHEM <= 0 ) THEN

               !=======================================================
               ! Skip over strat boxes
               !=======================================================
               IF ( ITS_IN_THE_STRAT( I, J, L, State_Met ) ) GOTO 40

               ! Increment JLOOP for trop boxes
               JLOOP          = JLOOP + 1
               JLOP(I,J,L)    = JLOOP

               ! test jlop for variable chem
               !temp(i,j,l) = dble(jloop)
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

               IF ( ITS_IN_THE_TROP( I, J, L, State_Met ) ) THEN

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
            VOLUME(JLOOP)  = BOXVL( I, J, L, State_Met )

            ! get air density in (molecs cm^-3)
            AIRDENS(JLOOP) = State_Met%AD(I,J,L) * 1000.d0 / 
     &                       VOLUME(JLOOP) * AVG / WTAIR

            ! get temperature
            T3(JLOOP)      = State_Met%T(I,J,L)

            ! PRESS3 = pressure in bar, multiply mb * 1000
            PRESS3(JLOOP)  = GET_PCENTER(I,J,L) * 1000d0

            ! Get relative humidity (here is absolute #H2O/cc air)
            ! AVGW is the mixing ratio of water vapor [v/v]
            ABSHUM(JLOOP)  = State_Met%AVGW(I,J,L) * AIRDENS(JLOOP)

            ! Go to next I
 40         CONTINUE           
         ENDDO
         ENDDO

         ! NIJLOOP is the number of surface boxes
         IF ( L == 1 ) NIJLOOP = JLOOP
      ENDDO

      ! NTLOOP is the number of total tropospheric boxes
      NTLOOP = JLOOP

      END SUBROUTINE RURALBOX
!EOC
