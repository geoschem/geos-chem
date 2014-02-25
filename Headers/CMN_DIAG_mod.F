!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: CMN_DIAG_mod
!
! !DESCRIPTION: Module CMN\_DIAG\_mod contains size parameters and global
!  variables for the GEOS-Chem diagnostic arrays.  This is mostly historical 
!  baggage.
!\\
!\\
! !INTERFACE:
!
      MODULE CMN_DIAG_MOD
!
! !USES:
!
      USE CMN_SIZE_MOD

      IMPLICIT NONE
      PUBLIC
!
! !DEFINED PARAMETERS:
!
      !=================================================================
      ! Maximum sizes of diagnostic arrays
      !
      ! Changed PD66 to 6  (bmy, 9/8/00)
      ! Changed PD21 to 10 (bmy, 9/30/00)
      ! Changed PD67 to 18 (bmy, 10/11/00)
      ! Changed PD46 to 4  (bmy, 1/2/01)
      ! Changed PD29 to 5  (bmy, 1/2/01)
      ! Changed PD11 to 7  (bmy, 9/4/01)
      ! Changed PD32 to 0  (bmy, 2/14/02)
      ! Changed PD21 to 20 (bmy, 2/27/02)
      ! Changed PD43 to 4  (bmy, 3/4/02) 
      ! Changed PD05 to 10 (bmy, 10/18/02)
      ! Changed PD44 to 30 (bmy, 11/19/02) 
      ! Changed PD43 to 5  
      ! Changed PD67 to 22 (bmy, 6/23/03)
      ! Changed PD66 to 5  (bmy, 6/23/03)
      ! Changed PD03 to 5  (bmy, 8/20/03)
      ! Changed PD37 to 10 (bmy, 1/21/04)
      ! Changed PD06 to NDSTBIN (bmy, 4/5/04)
      ! Changed PD07 to 7  (bmy, 4/5/04)
      ! Changed PD08 to 2  (bmy, 4/20/04)
      ! Changed PD07 to 12 (bmy, 7/15/04)
      ! Changed PD21 to 26 (bmy, 1/5/05)
      ! Removed PD03 -- now in "diag03_mod.f" (bmy, 1/21/05)
      ! Removed PD41 -- now in "diag41_mod.f" (bmy, 1/21/05)
      ! Now set PD09 to 6 (bmy, 6/27/05)
      ! Removed PD04 -- now in "diag04_mod.f" (bmy, 7/26/05)
      ! Now set PD30 to 1 (bmy, 8/18/05)
      ! Now set PD46 to 6 (tmf, 1/20/09)
      ! Now set PD10 to 20 (phs, 9/18/07)
      ! Changed PD17 to 8 (tmf, 1/7/09)
      ! Changed PD18 to 8 (tmf, 1/7/09)
      ! Changed PD22 to 8 (tmf, 1/7/09)
      ! Changed PD37 to 35 (tmf, 1/7/09)
      ! Changed PD38 to 35 (tmf, 1/7/09)
      ! Changed PD39 to 35 (tmf, 1/7/09)
      ! Changed PD44 to 41 (tmf, 1/7/09)
      ! Now set PD52 to 1 (jaegle 2/26/09)
      ! Increase PD46 from 6 to 13 (mpb, ccc, 11/19/09)
      ! increase PD21 from 20 to 27 (clh, 05/06/10)
      ! changed PD44 (drydep) to 53 (jpp, 6/13/09)
      ! changed PD39 (wetdep) to 38 (jpp, 7/08/09)
      ! Now set PD46 to 16 (jpp, 6/7/09)
      ! Changed PD11 from 7 to 5 (efischer, mpayer, 3/19/12)
      ! SOAupdate: PD46 increased from 16 to 20 (hotp 3/1/10)
      !=================================================================
      INTEGER, PARAMETER :: PD01=3            
      INTEGER, PARAMETER :: PD02=3
      INTEGER, PARAMETER :: PD05=10        
      INTEGER, PARAMETER :: PD06=NDSTBIN
      ! +3 for SEAC4RS SOA tracers (jaf, 6/25/13)
      INTEGER, PARAMETER :: PD07=15
      INTEGER, PARAMETER :: PD08=2
      INTEGER, PARAMETER :: PD09=6
      INTEGER, PARAMETER :: PD10=20            
      INTEGER, PARAMETER :: PD11=5
      INTEGER, PARAMETER :: PD12=0     
      INTEGER, PARAMETER :: PD13=1            
      INTEGER, PARAMETER :: PD14=NNPAR       
      INTEGER, PARAMETER :: PD15=NNPAR      
      INTEGER, PARAMETER :: PD16=2  
      INTEGER, PARAMETER :: PD17=8           
      INTEGER, PARAMETER :: PD18=8
      INTEGER, PARAMETER :: PD19=0            
      INTEGER, PARAMETER :: PD20=0 
      INTEGER, PARAMETER :: PD21=27
      INTEGER, PARAMETER :: PD22=14      !jpp replaced 8, 4/24/2011... for Br.
      INTEGER, PARAMETER :: PD23=0           
      INTEGER, PARAMETER :: PD24=NNPAR        
      INTEGER, PARAMETER :: PD25=NNPAR       
      INTEGER, PARAMETER :: PD26=NNPAR
      INTEGER, PARAMETER :: PD27=1            
      INTEGER, PARAMETER :: PD28=0           
      INTEGER, PARAMETER :: PD29=5
      INTEGER, PARAMETER :: PD30=1           
      INTEGER, PARAMETER :: PD31=1           
      INTEGER, PARAMETER :: PD32=1
      INTEGER, PARAMETER :: PD33=NNPAR        
      INTEGER, PARAMETER :: PD34=2           
      INTEGER, PARAMETER :: PD35=NNPAR 
      INTEGER, PARAMETER :: PD36=NNPAR        
      INTEGER, PARAMETER :: PD37=35        
      INTEGER, PARAMETER :: PD38=35 
      INTEGER, PARAMETER :: PD39=38      !jpp replaced 35
      INTEGER, PARAMETER :: PD40=4       
      INTEGER, PARAMETER :: PD43=5  
      INTEGER, PARAMETER :: PD44=53      !jpp replaced 41
      INTEGER, PARAMETER :: PD45=NNPAR+1      
      INTEGER, PARAMETER :: PD46=20
      INTEGER, PARAMETER :: PD47=NNPAR+1
      INTEGER, PARAMETER :: PD48=2            
      INTEGER, PARAMETER :: PD49=0                
      INTEGER, PARAMETER :: PD50=0            
      INTEGER, PARAMETER :: PD51=0           
      INTEGER, PARAMETER :: PD52=1          
      INTEGER, PARAMETER :: PD54=0           
      INTEGER, PARAMETER :: PD55=3
      ! Potential temperature diagnostic (hotp 7/31/09)
      INTEGER, PARAMETER :: PD57=1           
      INTEGER, PARAMETER :: PD58=12
#if   defined( TOMAS )
      ! Special settings for TOMAS aerosol microphysics (win, bmy, 1/22/10)
      INTEGER, PARAMETER :: PD59=TOMASBIN*TOMASSPEC    
      INTEGER, PARAMETER :: PD60=TOMASBIN*TOMASSPEC    
      INTEGER, PARAMETER :: PD61=2 
#else
      ! Normal settings for non-TOMAS simulations
      INTEGER, PARAMETER :: PD59=0            
      INTEGER, PARAMETER :: PD60=1
      INTEGER, PARAMETER :: PD61=0
#endif
      INTEGER, PARAMETER :: PD62=NNPAR        
      INTEGER, PARAMETER :: PD63=0           
      INTEGER, PARAMETER :: PD64=0
      INTEGER            :: PD65
      INTEGER, PARAMETER :: PD66=6           
      INTEGER, PARAMETER :: PD67=23 ! (Lin, 31/03/09)
      INTEGER, PARAMETER :: PD68=4            
      INTEGER, PARAMETER :: PD69=1
      INTEGER, PARAMETER :: PD70=0

      !=================================================================
      ! Variables for printing out selected tracers in diagnostic output
      !=================================================================
      INTEGER, PARAMETER :: MAX_DIAG   = 70
#if   defined( TOMAS )
      INTEGER, PARAMETER :: MAX_TRACER = NNPAR+1   ! For TOMAS (win, 1/25/10)
#elif defined( APM ) 
      INTEGER, PARAMETER :: MAX_TRACER = NNPAR+100 ! For APM (G. Luo 3/8/11)
#else
      INTEGER, PARAMETER :: MAX_TRACER = NNPAR+6   ! For non-TOMAS simulations 
#endif

!
! !PUBLIC DATA MEMBERS:
!
      !=================================================================
      ! Diagnostic counters & time variables
      !=================================================================
      INTEGER :: KDA48, NJDAY(366)

      !=================================================================
      ! Variables for the number of levels in multi-level diagnostics  
      ! Removed LD03 -- this is now in diag03_mod.f (bmy, 1/21/05)
      ! Added LD09 (bmy, 6/27/05)
      ! Added LD54 (phs, 9/22/06)
      ! Added LD10 (phs, 9/18/07)
      ! Added LD31 (bmy, 5/8/07)
      ! Added LD52 (jaegle, 02/26/09)
      ! Added LD59, LD60, LD61 (bmy, 1/22/10)
      ! Added LD57 (hotp 7/31/09)
      !=================================================================
      INTEGER :: LD12, LD13, LD14, LD15, LD16, LD17, LD18, LD21, LD22
      INTEGER :: LD24, LD25, LD26, LD37, LD38, LD39, LD43, LD45, LD47
      INTEGER :: LD54, LD64, LD65, LD66, LD68, LD01, LD02, LD05, LD07
      INTEGER :: LD09, LD10, LD31, LD52, LD19, LD57, LD58, LD59, LD60
      INTEGER :: LD61, LD62

      !=================================================================
      ! NDxx diagnostic flags
      !=================================================================
      INTEGER :: ND01, ND02, ND05, ND06, ND07, ND08, ND09, ND10, ND11
      INTEGER :: ND12, ND13, ND14, ND15, ND16, ND17, ND18, ND19, ND20
      INTEGER :: ND21, ND22, ND23, ND24, ND25, ND26, ND27, ND28, ND29
      INTEGER :: ND30, ND31, ND32, ND33, ND34, ND35, ND36, ND37, ND38
      INTEGER :: ND39, ND40, ND43, ND44, ND45, ND46, ND47, ND48, ND49
      INTEGER :: ND50, ND51, ND52, ND54, ND55, ND57, ND58, ND59
      INTEGER :: ND60, ND61, ND62, ND63, ND64, ND65, ND66, ND67, ND68
      INTEGER :: ND69, ND70, ND71, ND72, ND73, ND74, ND75

      !=================================================================
      ! Variables for printing out selected tracers in diagnostic output
      !=================================================================
      INTEGER :: TINDEX(MAX_DIAG,MAX_TRACER)          
      INTEGER :: TCOUNT(MAX_DIAG) 				  
      INTEGER :: TMAX(MAX_DIAG)					

      !=================================================================
      ! OH, J-Value, and 2-PM diagnostic arrays (bmy, 9/25/98)
      ! Move this here for now (bmy, 7/20/04)
      ! Removed HR1_NO and HR2_NO (mpayer, 11/8/13)
      !=================================================================
      REAL*8  :: HR1_JV, HR2_JV
      REAL*8  :: HR1_OH, HR2_OH, HR1_OTH, HR2_OTH
!
! !REMARKS:
!  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!  %%%  NOTE: THIS MODULE IS DEPRECATED.  MANY OF THE QUANTITIES STORED  %%%
!  %%%  in CMN_DIAG_MOD ARE NOW INCLUDED IN THE INPUT OPTIONS OBJECT,    %%%
!  %%%  BASED ON THE DERIVED TYPE IN Headers/gigc_input_opt_mod.F90.     %%%
!  %%%  (bmy, 11/19/12)                                                  %%%
!  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! !REVISION HISTORY:
!  (1 ) Changed RCS ID tag comment character from "C" to "!" to allow freeform
!        compilation.  Also added & continuation characters in column 73
!        to allow header files to be included in F90 freeform files.
!        Also converted PARAMETER statements to F90 syntax. (bmy, 6/25/02)
!  (2 ) Add LD05 for sulfate prod/loss (rjp, bdf, bmy, 9/20/02)
!  (3 )  Removed obsolete variables NTAU0, IDAY0, JDATE0, JYEAR0, KDACC, 
!         KDADYN, KDACONV, KDASRCE, KDACHEM, KDA3FLDS, KDA6FLDS, KDI6FLDS, 
!         KDKZZFLDS (bmy, 2/11/03)
!  (4 ) Fix for LINUX - remove & from column 73 (bmy, 6/27/03)
!  (5 ) Added LD03 for Kr85 Prod/loss diagnostic (bmy, 8/20/03)
!  (6 ) Removed obsolete arrays (bmy, 1/21/05)
!  (7 ) Rename MAXDIAG to MAX_DIAG and MAXTRACER to MAX_TRACER in order to
!        avoid name conflicts with "gamap_mod.f" (bmy, 5/3/05)
!  (8 ) Remove reference to TRCOFFSET (bmy, 5/16/06)
!  (9 ) Added multi level LD54 to common CDIAG1 (phs, 9/22/06)
!  (10) Added multi level LD10 to common CDIAG1. Set PD10 to 20.  Set
!        PD66 to 6. (phs, 9/18/07)
!  (11) Added LD52 to common CDIAG1 (jaegle, 02/26/09)
!  (12) Add GLYX, MGLY, SOAG, SOAM, and a few other tracers to AD17, AD18 
!        for archiving rainout and washout fraction. (tmf, 1/7/09)
!  (13) Add GLYX, MGLY J value archive. (tmf, 1/7/09)
!  (14) Add GLYX, MGLY, SOAG, SOAM, and a few other tracers to AD37, AD38, 
!        AD39 
!        for archiving rainout and washout flux. (tmf, 1/7/09)
!  (15) Add GLYX, MGLY, GLYC, 6 PANs, SOAG, SOAM dry dep, 
!        PD44 = 41.(tmf, 1/7/09)
!  (16) Add biogenic C2H4 emission, PD46 = 6. (tmf, 1/20/09)
!  (17) Add one met field to ND67 (EFLUX). (ccc, 5/14/09)
!  (18) Add declarations for PD58 and PD60, LD19, LD58. (kjw,8/18/09)
!  (19) Redimension PD59, PD60, PD61 for TOMAS microphysics.  Added LD59, LD60,
!        LD61 to common block.  Reset MAX_TRACER to NNPAR+1 for TOMAS.
!        (win, bmy, 1/22/10)
!  (20) Add LD57 and PD57 (potential temperature) (hotp, 3/15/10)
!  03 Aug 2010 - R. Yantosca - Added ProTeX headers
!  03 Aug 2010 - P. Kasibhatla & R. Yantosca - Now set MAX_TRACER to NNPAR+6 
!                                              to match ND09 diagnostic
!  09 Mar 2011 - R. Yantosca - Updated MAX_TRACER for APM (G. Luo)
!  03 Aug 2011 - M. Long     - Converted from Header file to Module
!  08 Nov 2013 - M. Sulprizio- Remove HR1_NO and HR2_NO. They are no longer
!                              needed for ND43 because NO, NO2, and NO3 are now
!                              tracers.
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
! !IROUTINE: init_cmn_diag
!
! !DESCRIPTION: Subroutine INIT\_CMN\_DIAG initializes quantities based on
!  the grid-independent size parameters.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE Init_CMN_DIAG( am_I_Root, RC )
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
!EOP
!------------------------------------------------------------------------------
!BOC
      ! Assume success
      RC   = GIGC_SUCCESS

      ! Set PD65 (for backwards compatibility, for now)
      PD65 = LLPAR*MAXFAM 

      END SUBROUTINE Init_CMN_DIAG
!EOC
      END MODULE CMN_DIAG_MOD


