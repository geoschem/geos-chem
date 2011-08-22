!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: CMN_O3_MOD
!
! !DESCRIPTION: COMMON BLOCKS FOR ANTHRO EMISSIONS (via SMVGEAR!)
!
! REMARKS:
! (1 ) Changed RCS ID tag comment character from "C" to "!" to allow freeform
!       compilation.  Also added & continuation characters in column 73
!       to allow header files to be included in F90 freeform files.
!       (bmy, 6/25/02)
! (2 ) Removed obsolete variables (bmy, 10/25/05)
! (3 ) Add anthropgenic emission species: BENZ, TOLU, XYLE, C2H2, C2H4. 
!      (tmf, 1/7/09)
! (4 ) Add EMISS_BVOC for VOC emissions (mpb, ccc, 11/19/09)
!
!
! NOTE: Now NEMPARA = max no. of anthropogenic emissions
!           NEMPARB = max no. of biogenic emissions
!\\
!\\
! !INTERFACE: 
!
      MODULE CMN_O3_MOD
!
! !USES:
!
      USE CMN_SIZE_MOD

      IMPLICIT NONE
!     PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS: KITCHEN SINK!
!
      PUBLIC
	SAVE
! !DEFINED PARAMETERS: 
!
!
! Rural Emissions: EMISRN = NOx (1:NOXLEVELS), EMISR = all other tracers
! Total Emissions: EMISTN = NOx (1:NOXLEVELS), EMIST = all other tracers
      REAL*8 ::
     &   EMISR(IIPAR,JJPAR,NEMPARA),  EMISRN(IIPAR,JJPAR,NOXLEVELS),   
     &   EMIST(IIPAR,JJPAR,NEMPARA),  EMISTN(IIPAR,JJPAR,NOXLEVELS)

! Rural Emissions:
! EMISRRN = NOx emissions into sigma levels L=1,NOXEXTENT
! EMISRR  = All other tracer emissions into sigma level L=1
      REAL*8 :: 
     &   EMISRR (IIPAR,JJPAR,NEMPARA+NEMPARB),                          
     &   EMISRRN(IIPAR,JJPAR,NOXEXTENT      )

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
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
      REAL*8 EMISS_BVOC(IIPAR,JJPAR,10)

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Arrays to read emissions from updated merge file : 
!        NOx, CO, PRPE, C3H8, ALK4, C2H6, ACET, MEK
! NOTE: ALD2 is not emitted in GEIA so we don't need an array for
!       it below...but it is emitted in EMEP.  It will be saved
!       into the EMISRR array for SMVGEAR. (bdf, bmy, 11/1/05)
      REAL*4 
     &  EMISTNOX (IIPAR,JJPAR,4,2), EMISTETHE(IIPAR,JJPAR    ),         &
     &  EMISTCO  (IIPAR,JJPAR    ), EMISTPRPE(IIPAR,JJPAR    ),         &
     &  EMISTC3H8(IIPAR,JJPAR    ), EMISTALK4(IIPAR,JJPAR    ),         &
     &  EMISTC2H6(IIPAR,JJPAR    ), EMISTSOX (IIPAR,JJPAR,4,2),         &
     &  EMISTACET(IIPAR,JJPAR    ), EMISTMEK (IIPAR,JJPAR    ),         &
     &  EMISTBENZ(IIPAR,JJPAR    ), EMISTTOLU(IIPAR,JJPAR    ),         &
     &  EMISTXYLE(IIPAR,JJPAR    ), EMISTC2H4(IIPAR,JJPAR    ),         & 
     &  EMISTC2H2(IIPAR,JJPAR    )
! Time of day and weekday/weekend scale factors
! NOTE: Now SCNR89 is (3,3) because of the weekday scale factor!!!
      REAL*8 TODH(6), TODN(6), TODB(6), SCNR89(3,3)

! IFSCLYR = Year to use for scaling fossil fuel emissions (1985 = no scaling!)
      INTEGER         FSCALYR

! FTOTCO2 = yearly scale factors based on Total  Fuel CO2 emissions 
! FLIQCO2 = yearly scale factors based on Liquid Fuel CO2 emissions
      REAL*4  FTOTCO2(IIPAR,JJPAR), FLIQCO2(IIPAR,JJPAR)

! FRACO3, FRACNO          = fractions of O3, NO
! SAVEOH, SAVENO, SAVENO3 = array to save OH, NO, & NO3 fields 
! SAVENO2                 = array to save NO2 fields (rvm, 5/9/00)
! FRACNO2                 = fraction of NO2 (rvm, bmy, 2/27/02)
! SAVEHO2                 = array to save HO2 fields (rvm, bmy, 2/27/02)
      REAL*8 ::
     &   FRACO3(IIPAR,JJPAR,LLPAR),                                     
     &   SAVEOH(IIPAR,JJPAR,LLPAR),                                     
     &   FRACNO(IIPAR,JJPAR,LLPAR),                                     
     &   SAVENO(IIPAR,JJPAR,LLPAR),                                     
     &   SAVENO2(IIPAR,JJPAR,LLPAR),                                     
     &   SAVENO3(IIPAR,JJPAR,LLPAR),                                    
     &   FRACNO2(IIPAR,JJPAR,LLPAR),                                    
     &   SAVEHO2(IIPAR,JJPAR,LLPAR)

	END MODULE CMN_O3_MOD
