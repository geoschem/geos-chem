!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: DEPO_MERCURY_MOD
!
! !DESCRIPTION: Module DEPO\_MERCURY\_MOD contains routines to handle
!  deposition fluxes for mercury. 
!
! !INTERFACE: 
!
      MODULE DEPO_MERCURY_MOD
!
! !USES:
! 
      IMPLICIT NONE
      PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
      PUBLIC :: ADD_Hg2_DD
      PUBLIC :: ADD_Hg2_WD
      PUBLIC :: ADD_HgP_DD
      PUBLIC :: ADD_HgP_WD
      PUBLIC :: ADD_HG2_SNOWPACK
      PUBLIC :: RESET_HG_DEP_ARRAYS
      PUBLIC :: INIT_DEPO_MERCURY
      PUBLIC :: CLEANUP_DEPO_MERCURY
!
! !PRIVATE MEMBER FUNCTIONS:
! 
!
! !PUBLIC DATA MEMBERS:
!  
      PUBLIC :: DD_HG2, DD_HGP, WD_HG2, WD_HGP
      PUBLIC :: SNOW_HG
      PUBLIC :: LHGSNOW
      REAL*8,  ALLOCATABLE :: DD_Hg2(:,:,:)
      REAL*8,  ALLOCATABLE :: DD_HgP(:,:,:)
      REAL*8,  ALLOCATABLE :: WD_Hg2(:,:,:)
      REAL*8,  ALLOCATABLE :: WD_HgP(:,:,:)
      REAL*8,  ALLOCATABLE :: SNOW_HG(:,:,:) !CDH Hg stored in snow+ice

      LOGICAL :: LHGSNOW
!
! !REVISION HISTORY:
! 23 Apr 2010 - C. Carouge  - Initial version
!
!EOP
!------------------------------------------------------------------------------

      CONTAINS


!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ADD_Hg2_DD
!
! !DESCRIPTION: Subroutine ADD_Hg2_DD computes the amount of Hg(II) dry deposited 
!  out of the atmosphere into the column array DD_Hg2. 
!  (sas, cdh, bmy, 1/19/05, 4/23/10)
!\\
!\\
! !INTERFACE: 
!
      SUBROUTINE ADD_Hg2_DD( I, J, N, DRY_Hg2)
!
! !USES
!
      USE TRACERID_MOD, ONLY : GET_Hg2_CAT
!
! !INPUT PARAMETERS: 
!
      INTEGER, INTENT(IN)   :: I, J, N   ! GEOS-Chem long, lat and tracer index
      REAL*8,  INTENT(IN)   :: DRY_Hg2   ! Hg(II) dry deposited out of the 
                                         ! atmosphere [kg]
!
! !REVISION HISTORY:
!  (1 ) DD_Hg2 is now a 3-D array.  Also pass N via the argument list. Now 
!        call GET_Hg2_CAT to return the Hg category #. (cdh, bmy, 3/28/06)
!  23 Apr 2010 - C. Carouge  - Moved from ocean_mercury_mod.f to
!                              depo_mercury_mod.f
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER               :: NN
      
      !=================================================================
      ! ADD_Hg2_DD begins here!
      !=================================================================

      ! Get the index for DD_Hg2 based on the tracer number
      NN = GET_Hg2_CAT( N )

      ! Store dry deposited Hg(II) into DD_Hg2 array
      IF ( NN > 0 ) THEN
         DD_Hg2(I,J,NN) = DD_Hg2(I,J,NN) + DRY_Hg2
        
      ENDIF
      
     
      ! Return to calling program
      END SUBROUTINE ADD_Hg2_DD
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ADD_Hg2_WD
!
! !DESCRIPTION: Subroutine ADD_Hg2_WD computes the amount of Hg(II) wet scavenged 
!  out of the atmosphere into the column array WD_Hg2. 
!  (sas, cdh, bmy, 1/19/05, 4/23/10)
!\\
!\\
! !INTERFACE: 
!
      SUBROUTINE ADD_Hg2_WD( I, J, N, WET_Hg2 )
!
! !USES
!
      USE TRACERID_MOD, ONLY : GET_Hg2_CAT
!
! !INPUT PARAMETERS: 
!
      INTEGER, INTENT(IN)   :: I, J, N   ! GEOS-Chem long, lat and tracer index
      REAL*8,  INTENT(IN)   :: WET_Hg2   ! Hg(II) scavenged out of the 
                                         ! atmosphere [kg]
!
! !REVISION HISTORY:
!  (1 ) WD_Hg2 is now a 3-D array.  Also pass N via the argument list. Now 
!        call GET_Hg2_CAT to return the Hg category #. (cdh, bmy, 3/28/06)
!  23 Apr 2010 - C. Carouge  - Moved from ocean_mercury_mod.f to
!                              depo_mercury_mod.f
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER               :: NN

      !=================================================================
      ! ADD_Hg2_WD begins here!
      !=================================================================

      ! Get Hg2 category number
      NN = GET_Hg2_CAT( N ) 
     
      ! Store wet deposited Hg(II) into WD_Hg2 array
      IF ( NN > 0 ) THEN
         WD_Hg2(I,J,NN) = WD_Hg2(I,J,NN) + WET_Hg2
         
      ENDIF

      ! Return to calling program
      END SUBROUTINE ADD_Hg2_WD
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ADD_HgP_DD
!
! !DESCRIPTION: Subroutine ADD_HgP_DD computes the amount of HgP dry deposited 
!  out of the atmosphere into the column array DD_HgP. 
!  (sas, cdh, bmy, 1/19/05, 4/23/10)
!\\
!\\
! !INTERFACE: 
!
      SUBROUTINE ADD_HgP_DD( I, J, N, DRY_HgP )
!
! !USES
!
      USE TRACERID_MOD, ONLY : GET_HgP_CAT
!
! !INPUT PARAMETERS: 
!
      INTEGER, INTENT(IN)   :: I, J, N   ! GEOS-Chem long, lat and tracer index
      REAL*8,  INTENT(IN)   :: DRY_HgP   ! HgP dry deposited out of the 
                                         ! atmosphere [kg]
!
! !REVISION HISTORY:
!  (1 ) DD_HgP is now a 3-D array.  Also pass N via the argument list. Now 
!        call GET_HgP_CAT to return the Hg category #. (cdh, bmy, 3/28/06)
!  23 Apr 2010 - C. Carouge  - Moved from ocean_mercury_mod.f to
!                              depo_mercury_mod.f
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
!
      INTEGER               :: NN

      !=================================================================
      ! ADD_HgP_DD begins here!
      !=================================================================
      
      ! Get the index for DD_Hg2 based on the tracer number
      NN = GET_HgP_CAT( N )

      ! Store dry deposited Hg(II) into DD_Hg2 array
      IF ( NN > 0 ) THEN
         DD_HgP(I,J,NN) = DD_HgP(I,J,NN) + DRY_HgP
        
      ENDIF

      ! Return to calling program
      END SUBROUTINE ADD_HgP_DD
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ADD_HgP_WD
!
! !DESCRIPTION: Subroutine ADD_HgP_WD computes the amount of HgP wet scavenged
!  out of the atmosphere into the column array WD_HgP. 
!  (sas, cdh, bmy, 1/19/05, 4/23/10)
!\\
!\\
! !INTERFACE: 
!
      SUBROUTINE ADD_HgP_WD( I, J, N, WET_HgP )
!
! !USES
!
      USE TRACERID_MOD, ONLY : GET_HgP_CAT
!
! !INPUT PARAMETERS: 
!
      INTEGER, INTENT(IN)   :: I, J, N   ! GEOS-Chem long, lat and tracer index
      REAL*8,  INTENT(IN)   :: WET_HgP   ! HgP scavenged out of the 
                                         ! atmosphere [kg]
!
! !REVISION HISTORY:
!  (1 ) WD_HgP is now a 3-D array.  Also pass N via the argument list. Now 
!        call GET_HgP_CAT to return the Hg category #. (cdh, bmy, 3/28/06)
!  23 Apr 2010 - C. Carouge  - Moved from ocean_mercury_mod.f to
!                              depo_mercury_mod.f
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER               :: NN

      !=================================================================
      ! ADD_Hg2_WD begins here!
      !=================================================================
      
      ! Get Hg2 category number
      NN = GET_HgP_CAT( N ) 

       ! Store wet deposited HgP into WD_HgP array
      IF ( NN > 0 ) THEN
         WD_HgP(I,J,NN) = WD_HgP(I,J,NN) + WET_HgP
        
      ENDIF
      
      ! Return to calling program
      END SUBROUTINE ADD_HgP_WD
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ADD_HG2_SNOWPACK
!
! !DESCRIPTION: Subroutine RESET_HG_DEP_ARRAYS adds Hg2 deposition to snowpack.
!  (cdh, 9/2/08, 4/23/10)
!\\
!\\
! !INTERFACE: 
!
      SUBROUTINE ADD_HG2_SNOWPACK( I, J, N, DEP_Hg2, SNOW_HT )
!
! !USES:
!
      USE DAO_MOD,           ONLY : SNOW, SNOMAS 
      USE DAO_MOD,           ONLY : IS_ICE
      USE TRACERID_MOD,      ONLY : GET_Hg2_CAT, GET_HgP_CAT
      USE TRACERID_MOD,      ONLY : IS_Hg2, IS_HgP
!
! !INPUT PARAMETERS:
!
      ! Arguments as input
      INTEGER, INTENT(IN)   :: I, J, N
      REAL*8,  INTENT(IN)   :: Dep_Hg2
      REAL*8,  INTENT(IN)   :: SNOW_HT
!
! !REVISION HISTORY:
!  23 Apr 2010 - C. Carouge  - Moved from mercury_mod.f to
!                              depo_mercury_mod.f
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      ! Local variables
!      REAL*8                :: SNOW_HT
      INTEGER               :: NN

      !=================================================================
      ! ADD_HG2_SNOWPACK begins here!
      !=================================================================
      
      ! Return if snowpack model is disabled
      IF (.NOT. LHGSNOW) RETURN

      IF ( IS_Hg2( N ) ) THEN
         ! Get Hg2 category number
         NN = GET_Hg2_CAT( N ) 
      ELSE IF ( IS_HgP( N ) ) THEN
         ! Get HgP category number
         NN = GET_HgP_CAT( N ) 
      ENDIF

!#if defined(GEOS_5)
!      ! GEOS5 snow height (water equivalent) in mm. (Docs wrongly say m)
!      SNOW_HT = SNOMAS(I,J)
!#else
!      ! GEOS1-4 snow heigt (water equivalent) in mm
!      SNOW_HT = SNOW(I,J)
!#endif 

      ! Check if there is snow on the ground, or if this is sea ice
      IF ( (SNOW_HT > 1d0) .OR. (IS_ICE(I,J)) ) THEN
    
         IF (DEP_HG2<0d0) THEN
            WRITE(6,'(3I6,2G12.4)') I,J,NN,DEP_HG2,SNOW_HG(I,J,NN)
         ENDIF

         SNOW_HG(I,J,NN) = SNOW_HG(I,J,NN) + MAX( DEP_HG2, 0D0 )

      ENDIF

      ! Return to calling program
      END SUBROUTINE ADD_HG2_SNOWPACK
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: RESET_HG_DEP_ARRAYS
!
! !DESCRIPTION: Subroutine RESET_HG_DEP_ARRAYS resets the wet and dry deposition arrays for
!  Hg(II) and Hg(p) to zero. This allows us to call OCEAN_MERCURY_FLUX and
!  LAND_MERCURY_FLUX in any order in MERCURY_MOD. (cdh, 9/2/08, 4/23/10)
!\\
!\\
! !INTERFACE: 
!
      SUBROUTINE RESET_HG_DEP_ARRAYS
!
! !REVISION HISTORY:
!  23 Apr 2010 - C. Carouge  - Moved from ocean_mercury_mod.f to
!                              depo_mercury_mod.f
!EOP
!------------------------------------------------------------------------------
!BOC
!
      ! Reset deposition arrays.
      DD_Hg2 = 0d0
      WD_Hg2 = 0d0
      DD_HgP = 0d0
      WD_HgP = 0d0

      END SUBROUTINE RESET_HG_DEP_ARRAYS
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: INIT_DEPO_MERCURY
!
! !DESCRIPTION: Subroutine INIT\_DEPO\_MERCURY initialize deposition arrays for
!  mercury (ccc, 4/23/10)
!\\
!\\
! !INTERFACE: 
!
      SUBROUTINE INIT_DEPO_MERCURY( )
!
! !USES
!
      USE ERROR_MOD,    ONLY : ALLOC_ERR
      USE TRACERID_MOD, ONLY : N_Hg_CATS

#     include "CMN_SIZE"     ! Size parameters
!
! !REVISION HISTORY:
!  23 Apr 2010 - C. Carouge   - Moved arrays allocation from ocean_mercury_mod.f
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER                      :: AS

      ! Allocate arrays
      ALLOCATE( DD_Hg2( IIPAR, JJPAR, N_Hg_CATS ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'DD_Hg2' )
      DD_Hg2 = 0d0

      ALLOCATE( DD_HgP( IIPAR, JJPAR, N_Hg_CATS ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'DD_HgP' )
      DD_HgP = 0d0

      ALLOCATE( WD_Hg2( IIPAR, JJPAR, N_Hg_CATS ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'WD_Hg2' )
      WD_Hg2 = 0d0

      ALLOCATE( WD_HgP( IIPAR, JJPAR, N_Hg_CATS ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'WD_HgP' )
      WD_HgP = 0d0

      ! CDH for snowpack
      ALLOCATE( SNOW_HG( IIPAR, JJPAR, N_Hg_CATS ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'SNOW_HG' )
      SNOW_HG = 0d0

      END SUBROUTINE INIT_DEPO_MERCURY
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: CLEANUP_DEPO_MERCURY
!
! !DESCRIPTION: Subroutine CLEANUP\_DEPO\_MERCURY deallocate all arrays
!  (ccc, 4/23/10)
!\\
!\\
! !INTERFACE: 
!
      SUBROUTINE CLEANUP_DEPO_MERCURY
!
! !REVISION HISTORY:
!  23 Apr 2010 - C. Carouge   - Moved from ocean_mercury_mod.f
!EOP
!------------------------------------------------------------------------------
!BOC
!
      IF ( ALLOCATED( DD_Hg2  ) ) DEALLOCATE( DD_Hg2  )
      IF ( ALLOCATED( DD_HgP  ) ) DEALLOCATE( DD_HgP  )
      IF ( ALLOCATED( WD_Hg2  ) ) DEALLOCATE( WD_Hg2  )
      IF ( ALLOCATED( WD_HgP  ) ) DEALLOCATE( WD_HgP  )
      IF ( ALLOCATED( SNOW_HG ) ) DEALLOCATE( SNOW_HG ) !CDH for snowpack

      END SUBROUTINE CLEANUP_DEPO_MERCURY

      END MODULE DEPO_MERCURY_MOD
!EOC
