!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: ohsave.F
!
! !DESCRIPTION: Subroutine OHSAVE stores the concentrations of OH and HO2 
!  for the ND43 diagnostic.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE OHSAVE( State_Met, State_Chm )
!
! !USES:
!
      USE CHEMGRID_MOD,       ONLY : ITS_IN_THE_NOCHEMGRID
      USE CMN_O3_MOD,         ONLY : SAVEOH, SAVEHO2, SAVEO1D, SAVEO3P
      USE CMN_SIZE_MOD
      USE State_Chm_Mod,      ONLY : ChmState
      USE State_Chm_Mod,      ONLY : Ind_
      USE State_Met_Mod,      ONLY : MetState
      USE PRECISION_MOD            ! For GEOS-Chem Precision (fp)

      IMPLICIT NONE
!
! !INPUT PARAMETERS: 
!
      TYPE(ChmState), INTENT(IN) :: State_Chm   ! Chemistry State object
      TYPE(MetState), INTENT(IN) :: State_Met   ! Meteorology State object
!
! !REVISION HISTORY: 
!  27 Feb 2002 - R. Yantosca - Initial version
!  (1 ) Original code from lwh, gmg, djj, jyl, etc, 1990's.  Modified for
!        GEOS-CHEM by Bob Yantosca et al.
!  (2 ) Added comment header and F90 declaration syntax.  Also now specify
!        the units of each variable for clarity. 
!  (3 ) Deleted NTRACER, it is not used.  Also added FRACNO2 and SAVEHO2
!        variables.  Updated comments, cosmetic changes (rvm, bmy, 2/27/02)
!  (4 ) Bug fix: swap the order of the lines where TMPNOX is computed.
!        Also deleted obsolete code from 2/02. (bmy, 7/31/02)
!  (5 ) Now reference IDTOX, IDTNOX, etc from "tracerid_mod.f". (1/13/03)
!  (6 ) Added OpenMP parallelization commands (bmy, 8/1/03)
!  (7 ) Now compute quantities for mean OH in "diag_oh_mod.f".  Now also
!        references STT from "tracer_mod.f".  Added N_TRACERS to the arg list.
!        Now dimension args XNUMOL, STT w/ N_TRACERS and not NNPAR. 
!        (bmy, 7/20/04)
!  (8 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (9 ) Reset FRAC* and SAVE* arrays, so that we don't carry dubious data
!        over from boxes that used to be in the tropopause but aren't anymore.
!        (phs, 1/19/07)
!  15 Sep 2010 - R. Yantosca - Added ProTeX headers
!  14 Mar 2013 - M. Payer    - FRACO3, FRACNO, and FRACNO2 are no longer needed
!                              because O3, NO, and NO2 are now tracers.
!  29 Mar 2013 - M. Payer    - Removed SAVENO, SAVENO2, SAVENO3.
!  31 May 2013 - R. Yantosca - Extra cleanup, remove N_TRACERS, XNUMOL, STT
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
!  17 Nov 2014 - M. Yannetti - Added PRECISION_MOD
!  21 Dec 2015 - M. Sulprizio- Get air density directly from State_Met object
!  16 Jun 2016 - J. Sheng    - Add index retriever
!  20 Jun 2016 - R. Yantosca - Now only define species ID's on the first call
!  20 Jun 2016 - R. Yantosca - Add checks to prevent array-OOB errors
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      ! Scalars
      INTEGER :: I, J, L

      ! SAVEd scalars
      LOGICAL, SAVE :: FIRST = .TRUE.
      INTEGER, SAVE :: id_OH, id_HO2, id_O1D, id_O3P
      LOGICAL, SAVE :: ok_OH, ok_HO2, ok_O1D, ok_O3P

      !=================================================================
      ! OHSAVE begins here!
      !=================================================================

      ! Reset because of variable tropopause.  Ensure that data for boxes 
      ! that once were in the troposphere, and are not anymore, do not
      ! carry dubious data over. (phs, 1/19/07)
      SAVEOH  = 0e+0_fp
      SAVEHO2 = 0e+0_fp
      SAVEO1D = 0e+0_fp
      SAVEO3P = 0e+0_fp

      ! First-time setup
      IF ( FIRST ) THEN 

         ! Only define species Id's on the first call (bmy, 6/20/16)
         id_OH  = Ind_('OH' )
         id_HO2 = Ind_('HO2')
         id_O1D = Ind_('O1D')
         id_O3P = Ind_('O'  )

         ! Set logicals to denote if each species is defined
         ok_OH  = ( id_OH  > 0 )
         ok_HO2 = ( id_HO2 > 0 )
         ok_O1D = ( id_O1D > 0 )
         ok_O3P = ( id_O3P > 0 )

         ! Reset first-time flag
         FIRST  = .FALSE.
      ENDIF

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L )
!$OMP+SCHEDULE( DYNAMIC )
      DO L = 1, LLPAR
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! Skip non-chemistry boxes
         IF ( ITS_IN_THE_NOCHEMGRID( I, J, L, State_Met ) ) CYCLE

         ! OH concentration [molec/cm3]
         IF ( ok_OH ) THEN
            SAVEOH(I,J,L)  = State_Chm%Species(I,J,L,id_OH)
         ENDIF

         ! HO2 concentration [v/v] 
         IF ( ok_HO2 ) THEN
            SAVEHO2(I,J,L) = State_Chm%Species(I,J,L,id_HO2) /
     &                       State_Met%AIRNUMDEN(I,J,L)
         ENDIF

#if defined( UCX )
         ! O1D concentration [molec/cm3]
         IF ( ok_O1D ) THEN
            SAVEO1D(I,J,L) = State_Chm%Species(I,J,L,id_O1D)
         ENDIF

         ! O3P concentration [molec/cm3]
         IF ( ok_O3P ) THEN
            SAVEO3P(I,J,L) = State_Chm%Species(I,J,L,id_O3P)
         ENDIF
#endif

      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      END SUBROUTINE OHSAVE
!EOC
