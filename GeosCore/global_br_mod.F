!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: global_Br_mod
!
! !DESCRIPTION: Module GLOBAL\_Br\_MOD contains variables and routines for 
!  reading the global monthly mean Br concentration from disk. 
!\\
!\\
! !INTERFACE: 
!
      MODULE GLOBAL_Br_MOD
!
! !USES:
!
      USE PRECISION_MOD    ! For GEOS-Chem Precision (fp, f4, f8)

      IMPLICIT NONE
      PRIVATE
!
! !PUBLIC DATA MEMBERS:
!
      ! Array to store global monthly mean BR field
      REAL(fp), PUBLIC, ALLOCATABLE :: BR_TROP(:,:,:)
      REAL(fp), PUBLIC, ALLOCATABLE :: BR_STRAT(:,:,:)
      REAL(fp), PUBLIC, ALLOCATABLE :: BR_MERGE(:,:,:)

      ! Array to store global monthly mean BrO field
      REAL(fp), PUBLIC, ALLOCATABLE :: BRO_TROP(:,:,:)
      REAL(fp), PUBLIC, ALLOCATABLE :: BRO_STRAT(:,:,:)
      REAL(fp), PUBLIC, ALLOCATABLE :: BRO_MERGE(:,:,:)

      ! Array to store global monthly J-BrO field
      REAL(fp), PUBLIC, ALLOCATABLE :: J_BRO(:,:,:)
!
! !PUBLIC MEMBER FUNCTIONS:
! 
      PUBLIC :: GET_GLOBAL_Br   
      PUBLIC :: INIT_GLOBAL_Br
      PUBLIC :: CLEANUP_GLOBAL_Br

! !REFERENCES
!  (1 ) Holmes, C. D., et al. (2006), Global lifetime of elemental mercury
!       against oxidation by atomic bromine in the free troposphere, Geophys.
!       Res. Lett., 33(20).
!  (2 ) Holmes, C.D., et al. (2010) Global atmospheric model for mercury 
!       including oxidation by bromine atoms, AC&P, 10, 12,037-12,057. 
!  (3 ) Parrella, J. et al. (2012), Tropospheric bromine chemistry:     
!       implications for present and pre-industrial ozone and mercury, ACP.
!
! !REVISION HISTORY:
!  05 Jul 2006 - C. Holmes   - Copied from "global_oh_mod.f"
!  01 Dec 2010 - R. Yantosca - Added ProTeX headers
!  19 Apr 2012 - E.S. Corbitt - Added LGCBROMINE to use GEOS-Chem bromine.
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
!  14 Nov 2014 - M. Yannetti - Added PRECISION_MOD
!  12 Mar 2015 - R. Yantosca - Remove bpch input
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !PRIVATE TYPES:
!
      ! Pointers to fields in the HEMCO data structure.
      ! These need to be declared REAL(f4), aka REAL*4.
      ! NOTE: These are globally SAVEd variables so we can
      ! nullify these in the declaration statement (bmy, 4/29/16)
      REAL(f4), POINTER :: Br_GC     (:,:,:) => NULL()
      REAL(f4), POINTER :: BrO_GC    (:,:,:) => NULL()
      REAL(f4), POINTER :: Br_GMI    (:,:,:) => NULL()
      REAL(f4), POINTER :: BrO_GMI   (:,:,:) => NULL()
      REAL(f4), POINTER :: Br_TOMCAT (:,:,:) => NULL()
      REAL(f4), POINTER :: BrO_TOMCAT(:,:,:) => NULL()
      REAL(f4), POINTER :: JBrO      (:,:,:) => NULL()

      CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_global_Br
!
! !DESCRIPTION: Subroutine GET\_GLOBAL\_Br reads global Br from 
!  binary punch files stored in the /data/ctm/GEOS\_MEAN directory.  This Br 
!  data is needed as oxidant for mercury chemistry.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE GET_GLOBAL_Br( am_I_Root, Input_Opt, 
     &                          State_Met, THISMONTH, RC )
!
! !USES:
!
      USE CMN_SIZE_MOD
      USE CHEMGRID_MOD,       ONLY : GET_TPAUSE_LEVEL
      USE ErrCode_Mod
      USE ERROR_MOD,          ONLY : ERROR_STOP
      USE HCO_INTERFACE_MOD,  ONLY : HcoState
      USE HCO_EmisList_Mod,   ONLY : HCO_GetPtr
      USE Input_Opt_Mod,      ONLY : OptInput
      USE OCEAN_MERCURY_MOD,  ONLY : LGCBROMINE     !eds 4/19/12
      USE State_Met_Mod,      ONLY : MetState
!
! !INPUT PARAMETERS: 
!
      LOGICAL,        INTENT(IN)  :: am_I_Root   ! Are we on the root CPU?
      TYPE(OptInput), INTENT(IN)  :: Input_Opt   ! Input Options object
      INTEGER,        INTENT(IN)  :: THISMONTH   ! Current month
      TYPE(MetState), INTENT(IN)  :: State_Met   ! Meteorology State object
!
! !OUTPUT PARAMETERS:
!
      INTEGER,        INTENT(OUT) :: RC          ! Success or failure?
!
! !REMARKS:
!  THIS IS A NEW VERSION OF THIS SUBROUTINE WHICH COMBINES Br CONCENTRATIONS
!  FROM MULTIPLE DATA SOURCES
!                                                                             .
!  ##########################################################################
!  #####    NOTE: BINARY PUNCH INPUT IS BEING PHASED OUT.  THIS DATA    #####
!  #####    WILL EVENTUALLY BE READ IN FROM netCDF FILES VIA HEMCO!     #####
!  #####       -- Bob Yantosca (05 Mar 2015)                            #####
!  ##########################################################################
! 
! !REVISION HISTORY: 
!  05 Jul 2006 - C. Holmes   - Copied from "global_oh_mod.f"
!  (1 ) GET_GLOBAL_BR assumes that we are reading global BR data that occupies
!        all CTM levels.  Contact Bob Yantosca (bmy@io.harvard.edu) for IDL
!        regridding code which will produce the appropriate BR files.
!  01 Dec 2010 - R. Yantosca - Added ProTeX headers
!  23 Jun 2014 - R. Yantosca - Now accept am_I_Root, Input_Opt, RC
!  06 Nov 2014 - R. Yantosca - Replace TRANSFER_3D_TROP with direct casts
!  17 Dec 2014 - R. Yantosca - Leave time/date variables as 8-byte
!  05 Mar 2015 - R. Yantosca - Add Input_Opt%RES_DIR to data path
!  12 Mar 2015 - R. Yantosca - Retire bpch input
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      ! Scalars
      INTEGER            :: I, J, L
      INTEGER            :: TPL
      LOGICAL, SAVE      :: FIRST = .TRUE. 

      ! Strings
      CHARACTER(LEN=255) :: LOC = 'GET_GLOBAL_Br (global_br_mod.F)'

      !=================================================================
      ! GET_GLOBAL_BR begins here!
      !=================================================================

      ! Assume success
      RC = GC_SUCCESS

      ! Allocate BR array, if this is the first call
      IF ( FIRST ) THEN
         CALL INIT_GLOBAL_BR()
         FIRST = .FALSE.
      ENDIF

      IF ( LGCBROMINE ) THEN

         !-----------------------------------------------------------------
         ! Read Br from GEOS-Chem
         !-----------------------------------------------------------------

         ! Get a pointer to the the GEOS-Chem Br from HEMCO (bmy, 3/11/15)
         CALL HCO_GetPtr( am_I_Root, HcoState, 'Br_GC', Br_GC,  RC )
         IF ( RC /= GC_SUCCESS ) THEN
            CALL ERROR_STOP ( 'Cannot get pointer to Br_GC', LOC )
         ENDIF

         ! Convert ppbv -> pptv
         Br_TROP  = Br_GC * 1e+3_fp                         
   
         !-----------------------------------------------------------------
         ! Read BrO from GEOS-Chem
         !-----------------------------------------------------------------

         ! Get a pointer to the GEOS-Chem BrO from HEMCO (bmy, 3/11/15)
         CALL HCO_GetPtr( am_I_Root, HcoState, 'BrO_GC', BrO_GC, RC )
         IF ( RC /= GC_SUCCESS ) THEN
            CALL ERROR_STOP ( 'Cannot get pointer to BrO_GC', LOC )
         ENDIF

         ! Convert ppbv -> pptv
         BrO_TROP  = BrO_GC * 1e+3_fp      

      ELSE

         !-----------------------------------------------------------------
         ! Read Br from pTOMCAT biogenic bromocarbons
         !-----------------------------------------------------------------

         ! Get a pointer to the TOMCAT Br from HEMCO (bmy, 3/11/15)
         CALL HCO_GetPtr( am_I_Root,  HcoState, 
     &                   'Br_TOMCAT', Br_TOMCAT, RC )
         IF ( RC /= GC_SUCCESS ) THEN
            CALL ERROR_STOP ( 'Cannot get pointer to Br_TOMCAT', LOC )
         ENDIF

         ! Store in Br_TROP array [pptv]
         Br_TROP = BR_TOMCAT

         !-----------------------------------------------------------------
         ! Read BrO from pTOMCAT biogenic bromocarbons
         !-----------------------------------------------------------------

         ! Get a pointer to the TOMCAT Br from HEMCO (bmy, 3/11/15)
         CALL HCO_GetPtr( am_I_Root,   HcoState, 
     &                   'BrO_TOMCAT', BrO_TOMCAT, RC )
         IF ( RC /= GC_SUCCESS ) THEN
            CALL ERROR_STOP ( 'Cannot get pointer to Br_TOMCAT', LOC )
         ENDIF

         ! Store in BrO_TROP array [pptv]
         BrO_TROP = BRO_TOMCAT

      ENDIF

      !-----------------------------------------------------------------
      ! Read Br from GMI for stratosphere
      !-----------------------------------------------------------------

      ! Get a pointer to the GMI Br from HEMCO (bmy, 3/11/15)
      CALL HCO_GetPtr( am_I_Root, HcoState, 'Br_GMI', Br_GMI, RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL ERROR_STOP ( 'Cannot get pointer to Br_GMI', LOC )
      ENDIF

      ! Store in BrO_TROP array [pptv]
      Br_STRAT = Br_GMI

      !-----------------------------------------------------------------
      ! Read BrO from GMI for stratosphere
      !-----------------------------------------------------------------

      ! Get a pointer to the GMI Br from HEMCO (bmy, 3/11/15)
      CALL HCO_GetPtr( am_I_Root, HcoState, 'BrO_GMI', BrO_GMI, RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL ERROR_STOP ( 'Cannot get pointer to BrO_GMI', LOC )
      ENDIF

      ! Store in BrO_TROP array [pptv]
      BrO_STRAT = BrO_GMI

      !-----------------------------------------------------------------
      ! Use pTOMCAT or GEOS-Chem exclusively in the troposphere.
      ! In the stratosphere, use the greater value from either COMBO or
      ! the tropospheric model. COMBO source gases include CH3Br and 
      ! halons, while pTOMCAT and GEOS-Chem includes CH3Br and 
      ! shorter-lived gases.
      !-----------------------------------------------------------------

      BR_MERGE  = BR_TROP
      BRO_MERGE = BRO_TROP

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, TPL )
!      DO I=1, IIPAR
!      DO J=1, JJPAR      
      DO J=1, JJPAR      
      DO I=1, IIPAR
         
         ! First layer in the stratosphere
         TPL = GET_TPAUSE_LEVEL( I, J, State_Met )

         BR_MERGE(I,J,TPL:LLPAR) = MERGE(
     &        BR_STRAT(I,J,TPL:LLPAR), 
     &        BR_TROP(I,J,TPL:LLPAR), 
     &        MASK=BR_STRAT(I,J,TPL:LLPAR)>BR_TROP(I,J,TPL:LLPAR) )

         BRO_MERGE(I,J,TPL:LLPAR) = MERGE(
     &        BRO_STRAT(I,J,TPL:LLPAR), 
     &        BRO_TROP(I,J,TPL:LLPAR), 
     &        MASK=BR_STRAT(I,J,TPL:LLPAR)>BR_TROP(I,J,TPL:LLPAR) )

      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      ! Get a pointer to the GMI Br from HEMCO (bmy, 3/11/15)
      CALL HCO_GetPtr( am_I_Root, HcoState, 'JBrO', JBrO, RC )
      IF ( RC /= GC_SUCCESS ) THEN
         CALL ERROR_STOP ( 'Cannot get pointer to JBrO', LOC )
      ENDIF

      ! Assign data from ARRAY2 to the module variable BR
      J_BrO(:,:,1:LLCHEM_FIX) = JBrO(:,:,1:LLCHEM_FIX)

      END SUBROUTINE GET_GLOBAL_Br
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_global_Br
!
! !DESCRIPTION: Subroutine INIT\_GLOBAL\_Br allocates and zeroes all
!  module arrays.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE INIT_GLOBAL_Br
!
! !USES:
!
      USE CMN_SIZE_MOD 
      USE ERROR_MOD, ONLY : ALLOC_ERR
! 
! !REVISION HISTORY: 
!  05 Jul 2006 - C. Holmes   - Copied from "global_oh_mod.f"
!  01 Dec 2010 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER :: AS

      !=================================================================
      ! INIT_GLOBAL_BR begins here!
      !=================================================================

      !-------------------------------------
      ! Br Arrays
      !-------------------------------------

      ! Allocate BR_TROP array
      ALLOCATE( BR_TROP( IIPAR, JJPAR, LLPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'BR_TROP' )
      BR_TROP = 0e+0_fp

      ! Allocate BR_STRAT array
      ALLOCATE( BR_STRAT( IIPAR, JJPAR, LLPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'BR_STRAT' )
      BR_STRAT = 0e+0_fp

      ! Allocate BR_MERGE array
      ALLOCATE( BR_MERGE( IIPAR, JJPAR, LLPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'BR_MERGE' )
      BR_MERGE = 0e+0_fp


      !-------------------------------------
      ! BrO Arrays
      !-------------------------------------

      ! Allocate J_BrO array
      ALLOCATE( J_BrO( IIPAR, JJPAR, LLCHEM ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'J_BrO' )
      J_BrO = 0e+0_fp

      ! Allocate BrO_TROP array
      ALLOCATE( BrO_TROP( IIPAR, JJPAR, LLPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'BrO_TROP' )
      BrO_TROP = 0e+0_fp

      ! Allocate BrO_STRAT array
      ALLOCATE( BrO_STRAT( IIPAR, JJPAR, LLPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'BrO_STRAT' )
      BrO_STRAT = 0e+0_fp

      ! Allocate BrO_MERGE array
      ALLOCATE( BrO_MERGE( IIPAR, JJPAR, LLPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'BrO_MERGE' )
      BrO_MERGE = 0e+0_fp

      END SUBROUTINE INIT_GLOBAL_BR
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: cleanup_global_Br
!
! !DESCRIPTION: Subroutine CLEANUP\_GLOBAL\_Br deallocates module arrays.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE CLEANUP_GLOBAL_Br
! 
! !REVISION HISTORY: 
!  05 Jul 2006 - C. Holmes   - Copied from "global_oh_mod.f"
!  01 Dec 2010 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
      !=================================================================
      ! CLEANUP_GLOBAL_BR begins here!
      !=================================================================
      IF ( ALLOCATED( BR_TROP   ) ) DEALLOCATE( BR_TROP   ) 
      IF ( ALLOCATED( BR_STRAT  ) ) DEALLOCATE( BR_STRAT  )  
      IF ( ALLOCATED( BR_MERGE  ) ) DEALLOCATE( BR_MERGE  ) 
      IF ( ALLOCATED( J_BrO     ) ) DEALLOCATE( J_BrO     ) 
      IF ( ALLOCATED( BrO_TROP  ) ) DEALLOCATE( BrO_TROP  ) 
      IF ( ALLOCATED( BrO_STRAT ) ) DEALLOCATE( BrO_STRAT )      
      IF ( ALLOCATED( BrO_MERGE ) ) DEALLOCATE( BrO_MERGE ) 

      ! Free pointers
      Br_GC      => NULL()
      BrO_GC     => NULL()
      Br_GMI     => NULL()
      BrO_GMI    => NULL()
      Br_TOMCAT  => NULL()
      BrO_TOMCAT => NULL()
      JBrO       => NULL()

      END SUBROUTINE CLEANUP_GLOBAL_Br
!EOC
      END MODULE GLOBAL_Br_MOD
