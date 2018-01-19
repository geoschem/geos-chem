!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: CMN_O3_mod.F
!
! !DESCRIPTION: Common blocks for anthro emissions (via SMVGEAR!)
!\\
!\\
! !INTERFACE: 
!
      MODULE CMN_O3_MOD
!
! !USES:
!
      USE CMN_SIZE_MOD, ONLY : IIPAR, JJPAR, LLPAR
      USE PRECISION_MOD
      
      IMPLICIT NONE
      PUBLIC
!
! !PUBLIC DATA MEMBERS:
!
      !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      !%%% NOTE: THESE ARE USED FOR DIAGNOSTICS WHICH ARE OBSOLETE %%%
      !%%% WE WILL EVENTUALLY REMOVE THESE IN THE NEAR FUTURE.     %%%
      !%%% (bmy, 5/22/15)                                          %%%
      !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      ! SAVEOH  = array to save OH  fields 
      ! SAVEHO2 = array to save HO2 fields (rvm, bmy, 2/27/02)
      ! SAVEEOH = array to save EOH fields (evf, 5/21/13)
      ! SAVEGLYX= array to save GLYX fields (eam, 8/25/14 )
      ! SAVEOA  = array for total organic aerosol (trc 6 of ND42, eam,7/10/14)
      !REAL(fp), ALLOCATABLE :: SAVEOH(:,:,:)
      !REAL(fp), ALLOCATABLE :: SAVEHO2(:,:,:)
      !REAL(fp), ALLOCATABLE :: SAVEO1D(:,:,:)
      !REAL(fp), ALLOCATABLE :: SAVEO3P(:,:,:)
      REAL(fp), ALLOCATABLE :: SAVEEOH(:,:,:)
      REAL(fp), ALLOCATABLE :: SAVEGLYX(:,:,:)
      REAL(fp), ALLOCATABLE :: SAVEOA(:,:,:)
!
! !REMARKS:
!  NOTE: Now NEMPARA = max no. of anthropogenic emissions
!            NEMPARB = max no. of biogenic emissions
!                                                                             .
!  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!  %%%  NOTE: THIS MODULE IS NOW DEPRECATED.  IT WILL BE REMOVED WHEN  %%%
!  %%%  THE GEOS-CHEM EMISSIONS MODULE IS UPDATED (SOMETIME SOON)      %%%
!  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! !REVISION HISTORY:
!  23 Aug 2011 - M. Long   - Converted to Module from Header file
!  29 Mar 2013 - M. Payer    - Removed FRACO3, FRACNO, FRACNO2, SAVENO, SAVENO2,
!                              and SAVENO3. They are no longer needed because
!                              O3, NO, NO2, and NO3 are now tracers.
!  13 Aug 2013 - M. Sulprizio- Increase last dimension of EMISS_BVOC to include
!                              sesquiterpenes and add EMISTNAP for SOA +
!                              semivolatile POA simulation (H. Pye)
!  21 Feb 2014 - M. Sulprizio- Added SAVEO1D and SAVEO3P for UCX (S.D. Eastham)
!  03 Dec 2014 - M. Yannetti - Added PRECISION_MOD
!  22 May 2015 - R. Yantosca - Remove arrays made obsolete by HEMCO
!  03 Jan 2017 - M. Sulprizio- Added SAVEEOH (E. Fischer)
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
! !IROUTINE: Init_Cmn_O3
!
! !DESCRIPTION: Subroutine INIT\_CMN\_O3 allocates all module arrays.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE Init_CMN_O3( am_I_Root, RC )
!
! !USES:
!
      USE ErrCode_Mod
!
! !INPUT PARAMETERS:
!
      LOGICAL, INTENT(IN)  :: am_I_Root   ! Are we on the root CPU?
!
! !OUTPUT PARAMETERS:
!
      INTEGER, INTENT(OUT) :: RC          ! Success or failure?
! 
! !REVISION HISTORY: 
!  19 Nov 2012 - R. Yantosca - Added ProTeX headers
!  01 Feb 2013 - R. Yantosca - Now allocate EMISS_BVOC to 19 emission species
!  29 Mar 2013 - M. Payer    - Removed FRACO3, FRACNO, FRACNO2, SAVENO, SAVENO2,
!                              and SAVENO3. They are no longer needed because
!                              O3, NO, NO2, and NO3 are now tracers.
!  13 Aug 2013 - M. Sulprizio- Increase last dimension of EMISS_BVOC to include
!                              sesquiterpenes and add EMISTNAP for SOA +
!                              semivolatile POA simulation (H. Pye)
!  22 May 2015 - R. Yantosca - Remove arrays made obsolete by HEMCO
!EOP
!------------------------------------------------------------------------------
!BOC
      ! Assume success
      RC = GC_SUCCESS

      ! Allocate arrays
      !ALLOCATE( SAVEOH (     IIPAR, JJPAR, LLPAR ), STAT=RC )
      !ALLOCATE( SAVEHO2(     IIPAR, JJPAR, LLPAR ), STAT=RC )
      !ALLOCATE( SAVEO1D(     IIPAR, JJPAR, LLPAR ), STAT=RC )
      !ALLOCATE( SAVEO3P(     IIPAR, JJPAR, LLPAR ), STAT=RC )
      ALLOCATE( SAVEEOH(     IIPAR, JJPAR, LLPAR ), STAT=RC )
      ALLOCATE( SAVEGLYX(    IIPAR, JJPAR, LLPAR ), STAT=RC )
      ALLOCATE( SAVEOA(      IIPAR, JJPAR, LLPAR ), STAT=RC )

      ! Zero arrays
      !SAVEOH     = 0e+0_fp
      !SAVEHO2    = 0e+0_fp
      !SAVEO1D    = 0e+0_fp
      !SAVEO3P    = 0e+0_fp
      SAVEEOH    = 0e+0_fp
      SAVEGLYX   = 0e+0_fp
      SAVEOA     = 0e+0_fp

      END SUBROUTINE Init_CMN_O3
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Cleanup_Cmn_O3
!
! !DESCRIPTION: Subroutine CLEANUP\_CMN\_O3 allocates all module arrays.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE Cleanup_CMN_O3( am_I_Root, RC )
!
! !USES:
!
      USE ErrCode_Mod
!
! !INPUT PARAMETERS:
!
      LOGICAL, INTENT(IN)  :: am_I_Root   ! Are we on the root CPU?
!
! !OUTPUT PARAMETERS:
!
      INTEGER, INTENT(OUT) :: RC          ! Success or failure?
! 
! !REVISION HISTORY: 
!  19 Nov 2012 - R. Yantosca - Added ProTeX headers
!  29 Mar 2013 - M. Payer    - Removed FRACO3, FRACNO, FRACNO2, SAVENO, SAVENO2,
!                              and SAVENO3. They are no longer needed because
!                              O3, NO, NO2, and NO3 are now tracers.
!  22 May 2015 - R. Yantosca - Remove arrays made obsolete by HEMCO
!EOP
!------------------------------------------------------------------------------
!BOC
      ! Assume success
      RC = GC_SUCCESS
      
      ! Deallocate arrays
      !IF ( ALLOCATED( SAVEOH      ) ) DEALLOCATE( SAVEOH      )
      !IF ( ALLOCATED( SAVEHO2     ) ) DEALLOCATE( SAVEHO2     )
      !IF ( ALLOCATED( SAVEO1D     ) ) DEALLOCATE( SAVEO1D     ) 
      !IF ( ALLOCATED( SAVEO3P     ) ) DEALLOCATE( SAVEO3P     )
      IF ( ALLOCATED( SAVEEOH     ) ) DEALLOCATE( SAVEEOH     )
      IF ( ALLOCATED( SAVEGLYX    ) ) DEALLOCATE( SAVEGLYX    )
      IF ( ALLOCATED( SAVEOA      ) ) DEALLOCATE( SAVEOA      )

      END SUBROUTINE Cleanup_CMN_O3
!EOC
      END MODULE CMN_O3_MOD

