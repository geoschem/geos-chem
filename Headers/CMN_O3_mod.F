!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: CMN_O3_mod
!
! !DESCRIPTION: Common blocks for anthro emissions (via SMVGEAR!)
!
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
#if defined( DEVEL )
      REAL*8, ALLOCATABLE  :: EMISR(:,:,:)
      REAL*8, ALLOCATABLE  :: EMISRN(:,:,:)
      REAL*8, ALLOCATABLE  :: EMIST(:,:,:)
      REAL*8, ALLOCATABLE  :: EMISTN(:,:,:)
#else
      REAL*8  :: EMISR(IIPAR,JJPAR,NEMPARA)
      REAL*8  :: EMISRN(IIPAR,JJPAR,NOXLEVELS)
      REAL*8  :: EMIST(IIPAR,JJPAR,NEMPARA)
      REAL*8  :: EMISTN(IIPAR,JJPAR,NOXLEVELS)
#endif

      ! Rural Emissions:
      ! EMISRRN = NOx emissions into sigma levels L=1,NOXEXTENT
      ! EMISRR  = All other tracer emissions into sigma level L=1
#if defined( DEVEL )
      REAL*8, ALLOCATABLE  :: EMISRR (:,:,:)
      REAL*8, ALLOCATABLE  :: EMISRRN(:,:,:)
#else
      REAL*8  :: EMISRR (IIPAR,JJPAR,NEMPARA+NEMPARB)
      REAL*8  :: EMISRRN(IIPAR,JJPAR,NOXEXTENT      )
#endif

      !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      !
      ! New biogenic VOC emissions (mpb,2009) 
      !
      !     ---------------------------
      !     Species            | Order 
      !     ---------------------------
      !     Isoprene           = 1
      !     Total Monoterpenes = 2 
      !     MBO                = 3 
      !     Alpha-Pinene       = 4
      !     Beta-Pinene        = 5
      !     Limonene           = 6
      !     Sabinene           = 7
      !     Mycrene            = 8
      !     3-Carene           = 9
      !     Ocimene            = 10
      !     ---------------------------	

      ! Define common block	
#if defined( DEVEL )
      REAL*8, ALLOCATABLE  :: EMISS_BVOC(:,:,:)
#else
      REAL*8  :: EMISS_BVOC(IIPAR,JJPAR,10)
#endif


      !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      ! Arrays to read emissions from updated merge file : 
      !        NOx, CO, PRPE, C3H8, ALK4, C2H6, ACET, MEK
      ! NOTE: ALD2 is not emitted in GEIA so we don't need an array for
      !       it below...but it is emitted in EMEP.  It will be saved
      !       into the EMISRR array for SMVGEAR. (bdf, bmy, 11/1/05)
      !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      ! Arrays to read emissions from updated merge file : 
      !        NOx, CO, PRPE, C3H8, ALK4, C2H6, ACET, MEK
      ! NOTE: ALD2 is not emitted in GEIA so we don't need an array for
      !       it below...but it is emitted in EMEP.  It will be saved

#if defined( DEVEL )
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
#else
      REAL*4  :: EMISTNOX (IIPAR,JJPAR,4,2)
      REAL*4  :: EMISTETHE(IIPAR,JJPAR    )
      REAL*4  :: EMISTCO  (IIPAR,JJPAR    )
      REAL*4  :: EMISTPRPE(IIPAR,JJPAR    )
      REAL*4  :: EMISTC3H8(IIPAR,JJPAR    )
      REAL*4  :: EMISTALK4(IIPAR,JJPAR    )
      REAL*4  :: EMISTC2H6(IIPAR,JJPAR    )
      REAL*4  :: EMISTSOX (IIPAR,JJPAR,4,2)
      REAL*4  :: EMISTACET(IIPAR,JJPAR    )
      REAL*4  :: EMISTMEK (IIPAR,JJPAR    )
      REAL*4  :: EMISTBENZ(IIPAR,JJPAR    )
      REAL*4  :: EMISTTOLU(IIPAR,JJPAR    )
      REAL*4  :: EMISTXYLE(IIPAR,JJPAR    )
      REAL*4  :: EMISTC2H4(IIPAR,JJPAR    )
      REAL*4  :: EMISTC2H2(IIPAR,JJPAR    )
#endif

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
#if defined( DEVEL )
      REAL*4, ALLOCATABLE  :: FTOTCO2(:,:)
      REAL*4, ALLOCATABLE  :: FLIQCO2(:,:)
#else
      REAL*4  :: FTOTCO2(IIPAR,JJPAR)
      REAL*4  :: FLIQCO2(IIPAR,JJPAR)
#endif

      ! FRACO3, FRACNO          = fractions of O3, NO
      ! SAVEOH, SAVENO, SAVENO3 = array to save OH, NO, & NO3 fields 
      ! SAVENO2                 = array to save NO2 fields (rvm, 5/9/00)
      ! FRACNO2                 = fraction of NO2 (rvm, bmy, 2/27/02)
      ! SAVEHO2                 = array to save HO2 fields (rvm, bmy, 2/27/02)
#if defined( DEVEL )
      REAL*8, ALLOCATABLE ::   FRACO3(:,:,:)
      REAL*8, ALLOCATABLE ::   SAVEOH(:,:,:)
      REAL*8, ALLOCATABLE ::   SAVEO1D(:,:,:)
      REAL*8, ALLOCATABLE ::   SAVEO3P(:,:,:)
      REAL*8, ALLOCATABLE ::   FRACNO(:,:,:)
      REAL*8, ALLOCATABLE ::   SAVENO(:,:,:)
      REAL*8, ALLOCATABLE ::   SAVENO2(:,:,:)
      REAL*8, ALLOCATABLE ::   SAVENO3(:,:,:)
      REAL*8, ALLOCATABLE ::   FRACNO2(:,:,:)
      REAL*8, ALLOCATABLE ::   SAVEHO2(:,:,:)
      REAL  , ALLOCATABLE ::   jvalues(:,:,:)
#else
      REAL*8 ::   FRACO3(IIPAR,JJPAR,LLPAR)
      REAL*8 ::   SAVEOH(IIPAR,JJPAR,LLPAR)
      REAL*8 ::   SAVEO1D(IIPAR,JJPAR,LLPAR)
      REAL*8 ::   SAVEO3P(IIPAR,JJPAR,LLPAR)
      REAL*8 ::   FRACNO(IIPAR,JJPAR,LLPAR)
      REAL*8 ::   SAVENO(IIPAR,JJPAR,LLPAR)
      REAL*8 ::   SAVENO2(IIPAR,JJPAR,LLPAR)
      REAL*8 ::   SAVENO3(IIPAR,JJPAR,LLPAR)
      REAL*8 ::   FRACNO2(IIPAR,JJPAR,LLPAR)
      REAL*8 ::   SAVEHO2(IIPAR,JJPAR,LLPAR)
      REAL   ::   jvalues(IIPAR,JJPAR,2)
#endif
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
!EOP
!------------------------------------------------------------------------------
!BOC
#if defined( DEVEL )
      CONTAINS

      SUBROUTINE SET_CMN_O3_MOD

      IMPLICIT NONE

      INTEGER AS

      ALLOCATE(
     &  EMISR     ( IIPAR, JJPAR, NEMPARA         ),
     &  EMISRN    ( IIPAR, JJPAR, NOXLEVELS       ),
     &  EMIST     ( IIPAR, JJPAR, NEMPARA         ),
     &  EMISTN    ( IIPAR, JJPAR, NOXLEVELS       ),
     &  EMISRR    ( IIPAR, JJPAR, NEMPARA+NEMPARB ),
     &  EMISRRN   ( IIPAR, JJPAR, NOXEXTENT       ),
     &  EMISS_BVOC( IIPAR, JJPAR, 10              ),
     &  EMISTNOX  ( IIPAR, JJPAR, 4, 2            ),
     &  EMISTETHE ( IIPAR, JJPAR                  ),
     &  EMISTCO   ( IIPAR, JJPAR                  ),
     &  EMISTPRPE ( IIPAR, JJPAR                  ),
     &  EMISTC3H8 ( IIPAR, JJPAR                  ),
     &  EMISTALK4 ( IIPAR, JJPAR                  ),
     &  EMISTC2H6 ( IIPAR, JJPAR                  ),
     &  EMISTSOX  ( IIPAR, JJPAR, 4, 2            ),
     &  EMISTACET ( IIPAR, JJPAR                  ),
     &  EMISTMEK  ( IIPAR, JJPAR                  ),
     &  EMISTBENZ ( IIPAR, JJPAR                  ),
     &  EMISTTOLU ( IIPAR, JJPAR                  ),
     &  EMISTXYLE ( IIPAR, JJPAR                  ),
     &  EMISTC2H4 ( IIPAR, JJPAR                  ),
     &  EMISTC2H2 ( IIPAR, JJPAR                  ),
     &  FTOTCO2   ( IIPAR, JJPAR                  ),
     &  FLIQCO2   ( IIPAR, JJPAR                  ),
     &  FRACO3    ( IIPAR, JJPAR, LLPAR           ),
     &  SAVEOH    ( IIPAR, JJPAR, LLPAR           ),
     &  FRACNO    ( IIPAR, JJPAR, LLPAR           ),
     &  SAVENO    ( IIPAR, JJPAR, LLPAR           ),
     &  SAVENO2   ( IIPAR, JJPAR, LLPAR           ),
     &  SAVENO3   ( IIPAR, JJPAR, LLPAR           ),
     &  FRACNO2   ( IIPAR, JJPAR, LLPAR           ),
     &  SAVEHO2   ( IIPAR, JJPAR, LLPAR           ),
     &  JVALUES   ( IIPAR, JJPAR, 2               ),
     &  STAT=AS)

      END SUBROUTINE SET_CMN_O3_MOD
#endif
      END MODULE CMN_O3_MOD
!EOC
