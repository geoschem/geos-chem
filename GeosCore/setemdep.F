!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: setemdep
!
! !DESCRIPTION: Subroutine SETEMDEP stores SMVGEAR reaction numbers 
!  (listed in "globchem.dat") corresponding to GEOS-CHEM tracers which emit 
!  and dry deposit into the NTEMIS and NTDEP index arrays.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE SETEMDEP( am_I_Root, Input_Opt, RC )
!
! !USES:
!
      USE TRACERID_MOD, ONLY : IDEMIS, IDTRMB, NEMANTHRO, NEMBIOG
      USE CMN_SIZE_MOD                   ! Size parameters
      USE COMODE_LOOP_MOD                ! SMVGEAR II arrays
      USE GIGC_ErrCode_Mod
      USE GIGC_Input_Opt_Mod, ONLY : OptInput

      IMPLICIT NONE
!
! !INPUT PARAMETERS:
!
      LOGICAL,        INTENT(IN) :: am_I_Root   ! Are we on the root CPU?
      TYPE(OptInput), INTENT(IN) :: Input_Opt   ! Input Options object
!
! !OUTPUT PARAMETERS:
!
      INTEGER,        INTENT(IN) :: RC          ! Success or failure
! 
! !REVISION HISTORY: 
!  03 Jan 1994 - L. Horowitz et al - - Initial version
!  (1 ) Now references "drydep_mod.f" and "tracerid_mod.f".  Updated comments
!        and made cosmetic changes. (bmy, 12/5/02)
!  (2 ) Cosmetic changes (bmy, 3/14/03)
!  (3 ) Updated for SMVGEAR II (gcc, bdf, bmy, 4/21/03)
!  (4 ) Now flag to "smv2.log" the emitted & dry-deposited tracers instead 
!        of flagging the tracers which aren't. (bmy, 7/20/04)
!  (5 ) Increase DEPNAME length format for output (fp, 3/15/10)
!  30 Jul 2012 - R. Yantosca - Now accept am_I_Root as an argument when
!                              running with the traditional driver main.F
!  30 Jul 2012 - R. Yantosca - Added ProTeX headers
!  13 Nov 2012 - R. Yantosca - Now pass Input_Opt and RC as arguments for GIGC
!  05 Mar 2013 - R. Yantosca - Now use Input_Opt%NUMDEP and Input_Opt%DEPNAME
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!     
      ! Scalars
      INTEGER            :: I, N, NK, NCS_TEMP
      INTEGER            :: NTRACER, NUMDEP
      CHARACTER(LEN=14)  :: NAME1

      ! Arrays
      CHARACTER(LEN=14 ) :: TRACER_NAME(Input_Opt%N_TRACERS)

      !=================================================================
      ! SETEMDEP begins here!
      !=================================================================      

      ! Write header to "smv2.log"
      IF ( am_I_Root ) THEN
         WRITE( IO93,'(/,a)' ) REPEAT( '=', 79 )
         WRITE( IO93,'(a)'   ) 'SETEMDEP: Emission & deposition species'
         WRITE( IO93,'(a,/)' ) REPEAT( '=', 79 )
      ENDIF

      ! Set local variables from the Input_Opt object
      NTRACER     = Input_Opt%N_TRACERS
      TRACER_NAME = Input_Opt%TRACER_NAME
      NUMDEP      = Input_Opt%NUMDEP

      !=================================================================
      ! Flag EMISSION REACTIONS in "globchem.dat" for GEOS-CHEM tracers
      !=================================================================

      ! Loop over different kinds of chemistry
      DO NCS = 1, NCSGAS

         ! Loop over GEOS-CHEM tracers
         DO I = 1, NTRACER 

            ! Rxn # for Ith GEOS-CHEM tracer
            NTEMIS(I,NCS) = 0

            ! Loop over emission species from "globchem.dat"
            DO N = 1, NEMIS(NCS)

               ! Rxn # for Nth emission species in "globchem.dat"
               NK = NKEMIS(N,NCS)

               ! Match "chem.dat" rxn number w/ GEOS-CHEM tracer number
               ! IRM is the species # for the first product of the NKth rxn
               ! IDTRMB is the species # of the GEOS-CHEM tracer which emits
               IF ( IDEMIS(I) /= 0 ) THEN
                  IF ( IRM(NPRODLO,NK,NCS) == IDTRMB(I,IDEMIS(I)) ) THEN
                     NTEMIS(I,NCS) = NK
                  ENDIF
               ENDIF
            ENDDO

            ! Flag emitted tracer
            IF ( NTEMIS(I,NCS) > 0 ) THEN
               IF ( am_I_Root ) WRITE( IO93, 100 ) I, TRACER_NAME(I)
 100           FORMAT( 'Tracer # ', i3, ' (', a4, ' ) has an ',
     &                 'emission rxn defined in "globchem.dat"' )
            ENDIF

         ENDDO

         ! The total # of emission species will be NEMANTHRO [anthro] + 
         ! NEMBIOG [bio], so reset NEMIS accordingly
         NEMIS(NCS) = NEMANTHRO + NEMBIOG

         ! Echo output to stdout
         IF ( am_I_Root ) WRITE( 6, 110 ) NEMIS(NCS)
 110     FORMAT( '     - SETEMDEP: Number of emitted '
     &           'species in "globchem.dat":', i3 )
      ENDDO                     

      !=================================================================
      ! Flag DRYDEP REACTIONS from "chem.dat" for each GEOS-CHEM tracer
      !=================================================================

      ! There is only drydep in the surface layer, which
      ! is accounted for in the "URBAN" chemistry slot
      NCS = NCSURBAN

      ! Loop over GEOS_CHEM drydep tracers
      DO I = 1, NUMDEP

         ! Rxn # of the Ith GEOS-CHEM drydep tracer
         NTDEP(I) = 0

         ! Loop over drydep species from "globchem.dat"
         DO N = 1, NDRYDEP(NCS)

            ! Rxn number and name of Nth drydep species in "globchem.dat" 
            NK    = NKDRY(N,NCS)
            NAME1 = NAMEGAS(IRM(1,NK,NCS))

            ! If we can match NAME1 against the GEOS-CHEM drydep tracer
            ! names in DEPNAME, then store the rxn number in NTDEP
            IF ( Input_Opt%DEPNAME(I) == NAME1 ) THEN
               NTDEP(I) = NK
               EXIT
            ENDIF
         ENDDO

         ! Flag drydep tracers
         IF ( NTDEP(I) > 0 ) THEN 
            IF ( am_I_Root ) WRITE( IO93, 120 ) I, Input_Opt%DEPNAME(I)
 120        FORMAT( 'Drydep species # ', i3, ' (', a5, 
     &              ') has a drydep rxn defined in "globchem.dat"' )
         ENDIF

      ENDDO

      ! Echo output to stdout
      IF ( am_I_Root ) WRITE( 6, 130 ) NDRYDEP(1)
 130  FORMAT( '     - SETEMDEP: Number of drydep  species '
     &        'in "globchem.dat":', i3 )

      IF ( am_I_Root ) WRITE( 6, 140 ) NUMDEP
 140  FORMAT( '     - SETEMDEP: Number of all GEOS-CHEM '
     &        'drydep species     :', i3 )

      ! Reset NCS = NCSURBAN, since we have defined our GEOS-CHEM
      ! mechanism in the urban slot of SMVGEAR II (bmy, 4/21/03)
      NCS = NCSURBAN

      END SUBROUTINE SETEMDEP
!EOC
