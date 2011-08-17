#if defined(ESMF_TESTBED_)
!------------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 910.1 and      !
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !PROGRAM: ut_GEOSCHEM
!
! !DESCRIPTION: Unit tester for the GEOS-Chem column code within the 
!  GEOSCHEMchem gridded component.  Uses the ESMF/MAPL environment.
!\\
!\\
! !INTERFACE:
!
#include "MAPL_Generic.h"

PROGRAM MAIN
!
! !USES: 
!
  USE ESMF_Mod                                     ! ESMF framework
  USE MAPL_Mod                                     ! MAPL framework
  USE GC_ESMF_COMP, ONLY: SetServices  ! To set IRF methods
  USE TESTBED_VALUE_MOD                            ! GEOS-Chem input values
  USE GC_ESMF_DRV

  IMPLICIT NONE
!
! !REMARKS:
!  ut_GEOSCHEM - Simple ESMF/MAPL example demonstrating how to call GEOSCHEM
!                                                                             .
!  It assumes 2 processors, so typically you will run it as
!
!     % mprirun -np 2 ut_GEOSCHEM.x
!                                                                             .
!  Arlindo da Silva <arlindo.dasilva@nasa.gov>, December 2009
!
! !REVISION HISTORY: 
!  01 Dec 2009 - A. Da Silva - Initial version  
!  02 Apr 2010 - R. Yantosca - Modified for GEOS-Chem column code
!  02 Apr 2010 - R. Yantosca - Added ProTex Headers, other cosmetic changes
!  07 Apr 2010 - R. Yantosca - Now populate all import/internal state fields
!  08 Apr 2010 - R. Yantosca - Now reference GEOS-Chem column input values from
!                              module "bmy_GC_Value_Mod.F90".  The header file
!                              "Column_Values_Saved_From_GEOS-Chem.h" is 
!                              now obsolete.
!  01 Jun 2010 - R. Yantosca - Increased output format to preserve precision
!  01 Jul 2010 - R. Yantosca - Now zero D_OH_MASS and D_AIR_MASS
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
!MSL_TMP  !==========================================================================
!MSL_TMP  ! Basic ESMF objects being used in this example
!MSL_TMP  !==========================================================================
!MSL_TMP  TYPE(ESMF_Grid)         :: Grid            ! Grid
!MSL_TMP  TYPE(ESMF_VM)           :: VM              ! ESMF Virtual Machine
!MSL_TMP  TYPE(ESMF_Time)         :: startTime       ! ESMF Time object (start time)
!MSL_TMP  TYPE(ESMF_Time)         :: stopTime        ! ESMF Time object (stop time)
!MSL_TMP  TYPE(ESMF_TimeInterval) :: TimeStep        ! ESMF TimeStep object 
!MSL_TMP
!MSL_TMP  !==========================================================================
!MSL_TMP  ! Grid component objects
!MSL_TMP  !==========================================================================
!MSL_TMP  TYPE(ESMF_GridComp)     :: GrComp          ! ESMF Gridded component
!MSL_TMP  TYPE(ESMF_State)        :: Import          ! ESMF Import state
!MSL_TMP  TYPE(ESMF_State)        :: Export
!MSL_TMP  TYPE(ESMF_Clock)        :: Clock
!MSL_TMP
!MSL_TMP  !==========================================================================
!MSL_TMP  ! Basic information about the parallel environment
!MSL_TMP  ! PET = Persistent Execution Threads
!MSL_TMP  ! In the current implementation, a PET is equivalent to an MPI process
!MSL_TMP  !==========================================================================
!MSL_TMP  INTEGER                 :: myPET           ! The local PET #
!MSL_TMP  INTEGER                 :: nPET            ! Total # of PETs we are using
!MSL_TMP  INTEGER                 :: STATUS          ! Status variable
!MSL_TMP  INTEGER                 :: RC              ! Status variable
!MSL_TMP  INTEGER                 :: I,  J, N        ! Loop indices 
!MSL_TMP  INTEGER                 :: IM, JM          ! Loop indices
!MSL_TMP  INTEGER, PARAMETER      :: NX       = 2    ! Layout: # of PETs in longitude
!MSL_TMP  INTEGER, PARAMETER      :: NY       = 1    ! Layout: # of PETs in latitude
!MSL_TMP  INTEGER, PARAMETER      :: IM_WORLD = 72   ! # of longitudes in global grid
!MSL_TMP  INTEGER, PARAMETER      :: JM_WORLD = 46   ! # of latitudes in global grid
!MSL_TMP  INTEGER, PARAMETER      :: LM_WORLD = 72   ! # of levels in global grid
!MSL_TMP  
!MSL_TMP  !==========================================================================
!MSL_TMP  ! Character variables
!MSL_TMP  !==========================================================================
!MSL_TMP  CHARACTER(LEN=ESMF_MAXSTR)  :: name
!MSL_TMP  CHARACTER(LEN=*), PARAMETER :: Iam = 'ut_GEOSCHEM'
!MSL_TMP
  !==========================================================================
  ! Call the Main program
  !==========================================================================
  CALL gc_Main_drv()

!EOC
END PROGRAM MAIN
#endif
