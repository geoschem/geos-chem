!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: optdepth_mod
!
! !DESCRIPTION: Module OPTDEPTH\_MOD contains routines to return optical
!  depths and update the ND21 diagnostic.
!\\
!\\
! !INTERFACE: 
!
      MODULE OPTDEPTH_MOD
!
! !USES:
!
      IMPLICIT NONE
      PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
      INTERFACE OPTDEPTH
         MODULE PROCEDURE OD_GEOS3_GEOS4
      END INTERFACE

      PUBLIC  :: OPTDEPTH
!
! !PRIVATE MEMBER FUNCTIONS:
! 
      PRIVATE :: OD_GEOS3_GEOS4
!
! !REVISION HISTORY:
!  15 Aug 2001 - R. Yantosca - Initial version
!  (1 ) Now add parallel DO-loops (bmy, 8/15/01)
!  (2 ) Removed obsolete code from 9/01 (bmy, 10/24/01)
!  (3 ) Now divide module header into MODULE PRIVATE, MODULE VARIABLES, and
!        MODULE ROUTINES sections.  Also add MODULE INTERFACES section,
!        since we have an interface here. (bmy, 5/28/02)
!  (4 ) Renamed OD_GEOS2_GEOS_3 to OD_GEOS3_GEOS4.  (bmy, 4/20/05)
!  (5 ) Remove support for GEOS-1 and GEOS-STRAT met fields (bmy, 8/4/06)
!  14 Sep 2010 - R. Yantosca - Added ProTeX headers
!  27 Nov 2012 - R. Yantosca - Added updates for GIGC
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
!EOP
!------------------------------------------------------------------------------
!BOC
      CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: od_geos3_geos4
!
! !DESCRIPTION: Subroutine OD\_GEOS3\_GEOS4 copies the DAO grid box optical 
!  depth from the OPTDEP met field array into the OPTD array.  Diagnostics 
!  are also archived. 
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE OD_GEOS3_GEOS4( am_I_Root, Input_Opt, State_Met, RC )
!
! !USES:
!
      USE CMN_DIAG_MOD
      USE CMN_SIZE_MOD
      USE DIAG_MOD,           ONLY: AD21
      USE GIGC_ErrCode_Mod
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
      USE GIGC_State_Met_Mod, ONLY : MetState
!
! !INPUT PARAMETERS: 
!
      LOGICAL,        INTENT(IN)    :: am_I_Root   ! Are we on the root CPU?
      TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
!
! !INPUT/OUTPUT PARAMETERS:
!
      TYPE(MetState), INTENT(INOUT) :: State_Met   ! Meteorology State object
!
! !OUTPUT PARAMETERS:
!
      INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
! 
! !REMARKS:
!  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!  %%%  NOTE: This routine is somewhat obsolete now.  It would be better  %%%
!  %%%  to save the values into State_Met%OPTD in the routines that read  %%%
!  %%%  the met fields from disk.  Keep as-is for now for compatibility   %%%
!  %%%  with existing code. (bmy, 11/27/12)
!  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!                                                                             .
!  The optical depths in the GEOS-5 met field archives are in-cloud optical
!  depths instead of grid-box optical depths (as was reported in the file
!  specification documents erroneously).
!                                                                             .
!  Also, the name "OD_GEOS3_GEOS4" is historical.  Once upon a time this was
!  used to denote the difference between the optical depths in GEOS-3 and
!  GEOS-4 (which come directly from the met fields) and GEOS-1 and GEOS-STRAT
!  (which were computed as functions of temperature).  The GEOS-5 and MERRA
!  optical depths are also provided in the met field archive, so the
!  algorithms in this routine are also equally applicable.
!
!  Original comments from the code:
!    GEOS-3/GEOS-4 optical depth is stored in the OPTDEP array,
!    which is read in routine "read_a6" of "dao_read_mod.f".
!                                                                             .
!    OPTDEP is archived every 6 hours, nevertheless, each chemistry
!    timestep we copy this into the OPTD array and archive for the
!    ND21 diagnostic.  This way the ND21 diagnostic is consistent
!    with GEOS-1/GEOS-STRAT.
!                                                                             .
!    OPTDEP and OPTD are dimensioned (LLPAR,IIPAR,JJPAR) to maximize
!    loop efficiency for processing an (I,J) column layer by layer.
!                                                                             .
!    Now also save CLDTOT to the ND21 diagnostic (bmy, 4/20/05)
!
! !REVISION HISTORY: 
!  15 Aug 2001 - R. Yantosca - Initial version
!  (1 ) Now parallelize I-J DO loops (bmy, 8/15/01)
!  (2 ) Renamed to OD_GEOS3_GEOS4.  Also now saves CLDF in AD21(I,J,L,2)
!        for the ND21 diagnostic (bmy, 4/20/05)
!  14 Sep 2010 - R. Yantosca - Added ProTeX headers
!  27 Nov 2012 - R. Yantosca - Now pass am_I_Root, Input_Opt and State_Met args
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER :: I, J, L

      !=================================================================
      ! OD_GEOS3_GEOS4 begins here!
      !================================================================= 

      ! Assume success
      RC = GIGC_SUCCESS

      ! Save optical depths and archive diagnostic
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L )
!$OMP+SCHEDULE( DYNAMIC )
      DO J = 1, JJPAR
      DO I = 1, IIPAR
      DO L = 1, LLPAR

         ! Copy optical depth over from OPTDEP array
         State_Met%OPTD(L,I,J) = State_Met%OPTDEP(L,I,J) 
         
         ! Save to AD21 array only if ND21 is turned on
         IF ( ND21 > 0 .and. L <= LD21 ) THEN
            AD21(I,J,L,1) = AD21(I,J,L,1) + State_Met%OPTD(L,I,J) 
            AD21(I,J,L,2) = AD21(I,J,L,2) + State_Met%CLDF(L,I,J)
         ENDIF 

      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      END SUBROUTINE OD_GEOS3_GEOS4
!EOC
      END MODULE OPTDEPTH_MOD





