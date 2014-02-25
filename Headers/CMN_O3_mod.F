!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: CMN_O3_mod
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
      USE CMN_SIZE_MOD, ONLY : IIPAR,     JJPAR,     LLPAR
      USE CMN_SIZE_MOD, ONLY : NEMPARA,   NEMPARB
      USE CMN_SIZE_MOD, ONLY : NOXLEVELS, NOXEXTENT
      
      IMPLICIT NONE
      PUBLIC
!
! !PUBLIC DATA MEMBERS:
!
      ! Rural Emissions: EMISRN = NOx (1:NOXLEVELS), EMISR = all other tracers
      ! Total Emissions: EMISTN = NOx (1:NOXLEVELS), EMIST = all other tracers
      REAL*8, ALLOCATABLE  :: EMISR(:,:,:)
      REAL*8, ALLOCATABLE  :: EMISRN(:,:,:)
      REAL*8, ALLOCATABLE  :: EMIST(:,:,:)
      REAL*8, ALLOCATABLE  :: EMISTN(:,:,:)

      ! Rural Emissions:
      ! EMISRRN = NOx emissions into sigma levels L=1,NOXEXTENT
      ! EMISRR  = All other tracer emissions into sigma level L=1
      REAL*8, ALLOCATABLE  :: EMISRR (:,:,:)
      REAL*8, ALLOCATABLE  :: EMISRRN(:,:,:)

      !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      !
      ! New biogenic VOC emissions (mpb,2009) 
      ! Now 19 species (dbm, 12/2012)
      !
      !     ---------------------------
      !     Species            | Order 
      !     ---------------------------
      !     Isoprene           = 1
      !     Acetone            = 2
      !     Propene            = 3
      !     Total Monoterpenes = 4 
      !     MBO                = 5 
      !     Ethene             = 6
      !     Alpha-Pinene       = 7
      !     Beta-Pinene        = 8
      !     Limonene           = 9
      !     Sabinene           = 10
      !     Mycrene            = 11
      !     3-Carene           = 12
      !     Ocimene            = 13
      !     Formic acid        = 14
      !     Acetic acid        = 15
      !     Acetaldehyde       = 16          
      !     Other monoterpenes = 17
      !     Methanol           = 18
      !     Ethanol            = 19
      !     ---------------------------	
      !     SOAupdate: Sesquiterpenes (hotp 3/1/10):
      !     Farnesene          = 11
      !     b-Caryophyllene    = 12
      !     Other sesquiterp   = 13
      !     Other monoterpenes = 14
      !     ---------------------------	

      ! Define common block	
      REAL*8, ALLOCATABLE  :: EMISS_BVOC(:,:,:)

      !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      ! Arrays to read emissions from updated merge file : 
      !        NOx, CO, PRPE, C3H8, ALK4, C2H6, ACET, MEK
      ! NOTE: ALD2 is not emitted in GEIA so we don't need an array for
      !       it below...but it is emitted in EMEP.  It will be saved
      !       into the EMISRR array for SMVGEAR. (bdf, bmy, 11/1/05)
      !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      REAL*4, ALLOCATABLE  :: EMISTNOX (:,:,:,:)
      REAL*4, ALLOCATABLE  :: EMISTETHE(:,:)
      REAL*4, ALLOCATABLE  :: EMISTCO  (:,:)
      REAL*4, ALLOCATABLE  :: EMISTPRPE(:,:)
      REAL*4, ALLOCATABLE  :: EMISTC3H8(:,:)
      REAL*4, ALLOCATABLE  :: EMISTALK4(:,:)
      REAL*4, ALLOCATABLE  :: EMISTC2H6(:,:)
      REAL*4, ALLOCATABLE  :: EMISTSOX (:,:,:,:)
      REAL*4, ALLOCATABLE  :: EMISTACET(:,:)
      REAL*4, ALLOCATABLE  :: EMISTMEK (:,:)
      REAL*4, ALLOCATABLE  :: EMISTBENZ(:,:)
      REAL*4, ALLOCATABLE  :: EMISTTOLU(:,:)
      REAL*4, ALLOCATABLE  :: EMISTXYLE(:,:)
      REAL*4, ALLOCATABLE  :: EMISTC2H4(:,:)
      REAL*4, ALLOCATABLE  :: EMISTC2H2(:,:)
      ! SOAupdate: Add EMISTNAP (hotp 7/21/10)
      REAL*4, ALLOCATABLE  :: EMISTNAP (:,:)


      ! Time of day and weekday/weekend scale factors
      ! NOTE: Now SCNR89 is (3,3) because of the weekday scale factor!!!
      REAL*8  :: TODH(6)
      REAL*8  :: TODN(6)
      REAL*8  :: TODB(6)
      REAL*8  :: SCNR89(3,3)


      ! IFSCLYR = Year to use for scaling fossil fuel emissions 
      ! (1985 = no scaling        !)
      INTEGER :: FSCALYR


      ! FTOTCO2 = yearly scale factors based on Total  Fuel CO2 emissions 
      ! FLIQCO2 = yearly scale factors based on Liquid Fuel CO2 emissions
      REAL*4, ALLOCATABLE  :: FTOTCO2(:,:)
      REAL*4, ALLOCATABLE  :: FLIQCO2(:,:)

      ! SAVEOH                  = array to save OH  fields 
      ! SAVEHO2                 = array to save HO2 fields (rvm, bmy, 2/27/02)
      REAL*8, ALLOCATABLE ::   SAVEOH(:,:,:)
      REAL*8, ALLOCATABLE ::   SAVEHO2(:,:,:)
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
! !IROUTINE: init_cmn_o3
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
      USE GIGC_ErrCode_Mod
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
!EOP
!------------------------------------------------------------------------------
!BOC
      ! Assume success
      RC = GIGC_SUCCESS

      ! Allocate arrays
      ALLOCATE( EMISR     ( IIPAR, JJPAR, NEMPARA         ), STAT=RC )
      ALLOCATE( EMISRN    ( IIPAR, JJPAR, NOXLEVELS       ), STAT=RC )
      ALLOCATE( EMIST     ( IIPAR, JJPAR, NEMPARA         ), STAT=RC )
      ALLOCATE( EMISTN    ( IIPAR, JJPAR, NOXLEVELS       ), STAT=RC )
      ALLOCATE( EMISRR    ( IIPAR, JJPAR, NEMPARA+NEMPARB ), STAT=RC )
      ALLOCATE( EMISRRN   ( IIPAR, JJPAR, NOXEXTENT       ), STAT=RC )
      ! SOAupdate: Now 14 (hotp 3/5/10)
      ALLOCATE( EMISS_BVOC( IIPAR, JJPAR, 14              ), STAT=RC )
      ALLOCATE( EMISTNOX  ( IIPAR, JJPAR, 4, 2            ), STAT=RC )
      ALLOCATE( EMISTETHE ( IIPAR, JJPAR                  ), STAT=RC )
      ALLOCATE( EMISTCO   ( IIPAR, JJPAR                  ), STAT=RC )
      ALLOCATE( EMISTPRPE ( IIPAR, JJPAR                  ), STAT=RC )
      ALLOCATE( EMISTC3H8 ( IIPAR, JJPAR                  ), STAT=RC )
      ALLOCATE( EMISTALK4 ( IIPAR, JJPAR                  ), STAT=RC )
      ALLOCATE( EMISTC2H6 ( IIPAR, JJPAR                  ), STAT=RC )
      ALLOCATE( EMISTSOX  ( IIPAR, JJPAR, 4, 2            ), STAT=RC )
      ALLOCATE( EMISTACET ( IIPAR, JJPAR                  ), STAT=RC )
      ALLOCATE( EMISTMEK  ( IIPAR, JJPAR                  ), STAT=RC )
      ALLOCATE( EMISTBENZ ( IIPAR, JJPAR                  ), STAT=RC )
      ALLOCATE( EMISTTOLU ( IIPAR, JJPAR                  ), STAT=RC )
      ALLOCATE( EMISTXYLE ( IIPAR, JJPAR                  ), STAT=RC )
      ALLOCATE( EMISTC2H4 ( IIPAR, JJPAR                  ), STAT=RC )
      ALLOCATE( EMISTC2H2 ( IIPAR, JJPAR                  ), STAT=RC )
      ! SOAupdate: Add EMISTNAP (hotp 7/21/10)
      ALLOCATE( EMISTNAP  ( IIPAR, JJPAR                  ), STAT=RC )
      ALLOCATE( FTOTCO2   ( IIPAR, JJPAR                  ), STAT=RC )
      ALLOCATE( FLIQCO2   ( IIPAR, JJPAR                  ), STAT=RC )
      ALLOCATE( SAVEOH    ( IIPAR, JJPAR, LLPAR           ), STAT=RC )
      ALLOCATE( SAVEHO2   ( IIPAR, JJPAR, LLPAR           ), STAT=RC )

      ! Zero arrays
      EMISR      = 0d0
      EMISRN     = 0d0
      EMIST      = 0d0
      EMISTN     = 0d0
      EMISRR     = 0d0
      EMISRRN    = 0d0
      EMISS_BVOC = 0d0
      EMISTNOX   = 0e0
      EMISTETHE  = 0e0
      EMISTCO    = 0e0
      EMISTPRPE  = 0e0
      EMISTC3H8  = 0e0
      EMISTALK4  = 0e0
      EMISTC2H6  = 0e0
      EMISTSOX   = 0e0
      EMISTACET  = 0e0
      EMISTMEK   = 0e0
      EMISTBENZ  = 0e0
      EMISTTOLU  = 0e0
      EMISTXYLE  = 0e0
      EMISTC2H4  = 0e0
      EMISTC2H2  = 0e0
      FTOTCO2    = 0d0
      FLIQCO2    = 0d0
      SAVEOH     = 0d0
      SAVEHO2    = 0d0

      END SUBROUTINE Init_CMN_O3
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: cleanup_cmn_o3
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
      USE GIGC_ErrCode_Mod
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
!EOP
!------------------------------------------------------------------------------
!BOC
      ! Assume success
      RC = GIGC_SUCCESS
      
      ! Deallocate arrays
      IF ( ALLOCATED( EMISR      ) ) DEALLOCATE( EMISR      )
      IF ( ALLOCATED( EMISRN     ) ) DEALLOCATE( EMISRN     )
      IF ( ALLOCATED( EMIST      ) ) DEALLOCATE( EMIST      )
      IF ( ALLOCATED( EMISTN     ) ) DEALLOCATE( EMISTN     )
      IF ( ALLOCATED( EMISRR     ) ) DEALLOCATE( EMISRR     )
      IF ( ALLOCATED( EMISRRN    ) ) DEALLOCATE( EMISRRN    )
      IF ( ALLOCATED( EMISS_BVOC ) ) DEALLOCATE( EMISS_BVOC )
      IF ( ALLOCATED( EMISTNOX   ) ) DEALLOCATE( EMISTNOX   )
      IF ( ALLOCATED( EMISTETHE  ) ) DEALLOCATE( EMISTETHE  )
      IF ( ALLOCATED( EMISTCO    ) ) DEALLOCATE( EMISTCO    )
      IF ( ALLOCATED( EMISTPRPE  ) ) DEALLOCATE( EMISTPRPE  )
      IF ( ALLOCATED( EMISTC3H8  ) ) DEALLOCATE( EMISTC3H8  )
      IF ( ALLOCATED( EMISTALK4  ) ) DEALLOCATE( EMISTALK4  )
      IF ( ALLOCATED( EMISTC2H6  ) ) DEALLOCATE( EMISTC2H6  )
      IF ( ALLOCATED( EMISTSOX   ) ) DEALLOCATE( EMISTSOX   )
      IF ( ALLOCATED( EMISTACET  ) ) DEALLOCATE( EMISTACET  )
      IF ( ALLOCATED( EMISTMEK   ) ) DEALLOCATE( EMISTMEK   )
      IF ( ALLOCATED( EMISTBENZ  ) ) DEALLOCATE( EMISTBENZ  )
      IF ( ALLOCATED( EMISTTOLU  ) ) DEALLOCATE( EMISTTOLU  )
      IF ( ALLOCATED( EMISTXYLE  ) ) DEALLOCATE( EMISTXYLE  )
      IF ( ALLOCATED( EMISTC2H4  ) ) DEALLOCATE( EMISTC2H4  )
      IF ( ALLOCATED( EMISTC2H2  ) ) DEALLOCATE( EMISTC2H2  )
      IF ( ALLOCATED( EMISTNAP   ) ) DEALLOCATE( EMISTNAP   )
      IF ( ALLOCATED( FTOTCO2    ) ) DEALLOCATE( FTOTCO2    )
      IF ( ALLOCATED( FLIQCO2    ) ) DEALLOCATE( FLIQCO2    )
      IF ( ALLOCATED( SAVEOH     ) ) DEALLOCATE( SAVEOH     )
      IF ( ALLOCATED( SAVEHO2    ) ) DEALLOCATE( SAVEHO2    )

      END SUBROUTINE Cleanup_CMN_O3
!EOC
      END MODULE CMN_O3_MOD

