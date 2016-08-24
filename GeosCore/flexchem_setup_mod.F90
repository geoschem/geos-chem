!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: flexchem_setup_mod
!
! !DESCRIPTION: Module FLEXCHEM\_SETUP\_MOD 
!\\
!\\
! !INTERFACE: 
!
MODULE FLEXCHEM_SETUP_MOD
!
! !USES:
!
  USE PRECISION_MOD            ! For GEOS-Chem Precision (fp)

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC :: INIT_FLEXCHEM
  PUBLIC :: HSAVE_KPP
!    
! !REVISION HISTORY:
!  14 Dec 2015 - M.S. Long   - Initial version
!  15 Jun 2016 - M. Sulprizio- Remove STTTOCSPEC mapping array. Species and
!                              tracers have a 1:1 mapping currently so mapping
!                              is not required
!  18 Jul 2016 - M. Sulprizio- Remove FAMILIES_KLUDGE routine. Family tracers
!                              have been eliminated.
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !PRIVATE TYPES:
!
  ! Save the H value for the Rosenbrock solver
  REAL(fp),  ALLOCATABLE :: HSAVE_KPP(:,:,:) 

CONTAINS
!EOC  
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: INIT_FLEXCHEM
!
! !DESCRIPTION: Subroutine INIT\_FLEXCHEM is used to allocate
!               arrays for the KPP solver. 
!\\
!\\
! !INTERFACE:
!  
  SUBROUTINE INIT_FLEXCHEM( am_I_Root, RC)  
!
! !USES:
!
    USE CMN_SIZE_MOD
    USE ErrCode_Mod
    USE gckpp_Global,       ONLY : NREACT
    USE gckpp_Monitor,      ONLY : EQN_NAMES
!
! !INPUT PARAMETERS:
!    
    LOGICAL,        INTENT(IN)    :: am_I_Root ! Is this the root CPU?
!
! !OUTPUT PARAMETERS:
!    
    INTEGER,        INTENT(OUT)   :: RC        ! Success or failure?
!    
! !REVISION HISTORY:
!  16 Sep 2009 - P. Le Sager - init
!  09 Dec 2009 - C. Carouge  - R_KPP and CSPEC_FOR_KPP are not used anymore
!  02 Aug 2012 - R. Yantosca - Now use am_I_Root to print on root CPU
!  22 Dec 2015 - M. Sulprizio- Use State_Met%AIRNUMDEN to convert initial
!                              species concentrations from v/v to molec/cm3
!  29 Jan 2016 - M. Sulprizio- Add calls to Register_Tracer and Register_Species
!                              to populate Tracer_Name, Tracer_Id, Species_Name,
!                              and Species_ID fields in State_Chm
!  06 Jun 2016 - M. Sulprizio- Replace NTSPEC with State_Chm%nSpecies and
!                              NAMEGAS with SpcInfo%Name from species database
!  06 Jun 2016 - M. Sulprizio- Replace Get_Indx with Spc_GetIndx to use the
!                              fast-species lookup from the species database
!  06 Jun 2016 - M. Sulprizio- Remove calls to Register_Tracer and
!                              Register_Species; these routines were made
!                              obsolete by the species database
!  14 Jun 2016 - M. Sulprizio- Replace Spc_GetIndx with Ind_ (M. Long)
!  25 Jul 2016 - E. Lundgren - Add check that species was not in restart file
!                              prior to v/v -> molec/cm3 conversion
!  02 Aug 2016 - E. Lundgren - Move unit conversion of species background
!                              values to restart_mod
!  24 Aug 2016 - M. Sulprizio- Remove CSPECTOKPP array. State_Chm%Map_KppSpc is
!                              now used instead.
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER            :: D

    !=================================================================
    ! INIT_FLEXCHEM begins here!
    !=================================================================

    ! Assume success
    RC = GC_SUCCESS

    IF ( am_I_Root ) WRITE( 6, 100 )
100 FORMAT( '     - INIT_FLEXCHEM: Allocating arrays for FLEX_CHEMISTRY' )

    ALLOCATE( HSAVE_KPP( IIPAR, JJPAR, LLCHEM ), STAT=RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    HSAVE_KPP = 0.d0

    IF (am_I_Root) THEN
       write(*,'(a)') ' KPP Reaction Reference '
       DO D = 1,NREACT
          write(*,'(i,a3,a85)') D,' | ',EQN_NAMES(D)
       END DO
    ENDIF

    ! Return to calling program
    RETURN

  END SUBROUTINE INIT_FLEXCHEM
!EOC
END MODULE FLEXCHEM_SETUP_MOD

