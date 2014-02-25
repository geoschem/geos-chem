!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: gckpp_comode_mod
!
! !DESCRIPTION: Module GCKPP\_COMODE\_MOD is used to allocate
!               and deallocate arrays for the KPP solver. 
!\\
!\\
! !INTERFACE: 
!
MODULE GCKPP_COMODE_MOD
! 
! !USES:
!  
  IMPLICIT NONE

!--- Previous to (ccc, 12/9/09)
!  REAL*8,  ALLOCATABLE :: R_KPP(:,:)
!  REAL*8,  ALLOCATABLE :: CSPEC_FOR_KPP(:,:) 

  REAL*8,  ALLOCATABLE :: HSAVE_KPP(:,:,:) 
!    
! !REVISION HISTORY:
!   16 Sep 2009 - P. Le Sager - init
!   03 Dec 2009 - C. Carouge  - CSPEC_FOR_KPP not used anymore 
!                               (use CSPEC instead)
!   09 Dec 2009 - C. Carouge  - R_KPP not used anymore
!                               (use RRATE_FOR_KPP instead)
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
! !IROUTINE: INIT_GCKPP_COMODE
!
! !DESCRIPTION: Subroutine INIT\_GCKPP\_COMODE is used to allocate
!               arrays for the KPP solver. 
!\\
!\\
! !INTERFACE:
!  
  SUBROUTINE INIT_GCKPP_COMODE( am_I_Root, IIPAR,   JJPAR, LLTROP, & 
                                ITLOOP,    NMTRATE, IGAS,  RC )
!
! !INPUT PARAMETERS:
!    
    LOGICAL, INTENT(IN) :: am_I_Root   ! Is this the root CPU?
    INTEGER, INTENT(IN) :: IIPAR, JJPAR, LLTROP, ITLOOP, NMTRATE, IGAS
!
! !OUTPUT PARAMETERS:
!    
    INTEGER, INTENT(OUT):: RC
!    
! !REVISION HISTORY:
!  16 Sep 2009 - P. Le Sager - init
!  09 Dec 2009 - C. Carouge  - R_KPP and CSPEC_FOR_KPP are not used anymore
!  02 Aug 2012 - R. Yantosca - Now use am_I_Root to print on root CPU
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER :: AS

    RC=0

    IF ( am_I_Root ) WRITE( 6, 100 )
100 FORMAT( '     - INIT_GCKPP_COMODE: Allocating arrays for GCKPP...' )


!--- Previous to (ccc, 12/9/09)
!    ALLOCATE( R_KPP( ITLOOP, NMTRATE ), STAT=AS )
!    IF ( AS /= 0 ) THEN
!       RC=1
!       RETURN
!    ENDIF
!    R_KPP = 0d0
!    ALLOCATE( R_KPP( 24, NMTRATE ), STAT=AS )
!    IF ( AS /= 0 ) THEN
!       RC=1
!       RETURN
!    ENDIF
!    R_KPP = 0d0

    ALLOCATE( HSAVE_KPP( IIPAR, JJPAR, LLTROP ), STAT=AS )
    IF ( AS /= 0 ) THEN
       RC=1
       RETURN
    ENDIF
    HSAVE_KPP = 0.d0 

!--- Previous to (ccc, 12/3/09)
!    ALLOCATE( CSPEC_FOR_KPP( ITLOOP, IGAS ), STAT=AS )
!    IF ( AS /= 0 ) THEN
!       RC=1
!       RETURN
!    ENDIF
!    CSPEC_FOR_KPP = 0d0

    ! Return to calling program
  END SUBROUTINE INIT_GCKPP_COMODE
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: CLEANUP_GCKPP_COMODE
!
! !DESCRIPTION: Subroutine CLEANUP\_GCKPP\_COMODE is used to deallocate
!               arrays for the KPP solver. 
!\\
!\\
! !INTERFACE:
!    
  SUBROUTINE CLEANUP_GCKPP_COMODE
!    
! !REVISION HISTORY:
!   16 Sep 2009 - P. Le Sager - init
!   09 Dec 2009 - C. Carouge  - R_KPP and CSPEC_FOR_KPP are not used anymore
!EOP
!------------------------------------------------------------------------------
!BOC    
!--- Previous to (ccc, 12/09/09)
!    IF ( ALLOCATED( R_KPP         ) ) DEALLOCATE( R_KPP         )
!    IF ( ALLOCATED( CSPEC_FOR_KPP ) ) DEALLOCATE( CSPEC_FOR_KPP )

    IF ( ALLOCATED( HSAVE_KPP     ) ) DEALLOCATE( HSAVE_KPP     )

    ! Return to calling program
  END SUBROUTINE CLEANUP_GCKPP_COMODE
!EOC
  
END MODULE GCKPP_COMODE_MOD

