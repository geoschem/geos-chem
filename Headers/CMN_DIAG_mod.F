!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: CMN_DIAG_mod.F
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
      USE CMN_FJX_MOD,   ONLY : W_
      USE PRECISION_MOD, ONLY : fpp => fp    ! For GEOS-Chem Precision (fp)

      IMPLICIT NONE
      PUBLIC
!
! !DEFINED PARAMETERS:
!
      !######################################################
      !##### NOTE: We need to keep these variables!     #####
      !######################################################
      INTEGER, PARAMETER :: PD40=4
      INTEGER            :: ND40
      INTEGER            :: ND48, ND49, ND50, ND51, ND63 ! timeseries

      !######################################################
      !##### NOTE: These variables can be blocked out!  #####
      !######################################################
#if defined( BPCH_DIAG ) || defined( BPCH_TIMESER )
      INTEGER, PARAMETER :: PD05=20
      INTEGER, PARAMETER :: PD06=NDSTBIN
      ! +3 for SEAC4RS SOA tracers (jaf, 6/25/13)
      INTEGER, PARAMETER :: PD07=15
      INTEGER, PARAMETER :: PD08=4
      INTEGER, PARAMETER :: PD09=6
      INTEGER, PARAMETER :: PD10=20            
      INTEGER, PARAMETER :: PD11=5
      INTEGER, PARAMETER :: PD16=2  
      INTEGER, PARAMETER :: PD17=8           
      INTEGER, PARAMETER :: PD18=8
      INTEGER, PARAMETER :: PD19=3
#if defined( BPCH_DIAG ) || defined( BPCH_TIMESER ) || defined( MODEL_GEOS )
      INTEGER, PARAMETER :: PD21=8+(NRHAER+NDUST)*3+(NRHAER*2)+
     &                           (NSTRATAER*3)
#endif
      INTEGER, PARAMETER :: PD22=76
      INTEGER, PARAMETER :: PD43=4
      INTEGER, PARAMETER :: PD52=6
      INTEGER, PARAMETER :: PD55=3
      ! Potential temperature diagnostic (hotp 7/31/09)
      INTEGER, PARAMETER :: PD57=1           
      INTEGER, PARAMETER :: PD58=15
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
      INTEGER            :: PD64
      INTEGER            :: PD65
      INTEGER, PARAMETER :: PD66=7
      INTEGER, PARAMETER :: PD67=23 ! (Lin, 31/03/09)
      INTEGER, PARAMETER :: PD68=8     
      INTEGER, PARAMETER :: PD69=1

      !number of rad flux and optics output types
      !8 flux and 3*3=9 optics
      INTEGER, PARAMETER :: PD72R=17
      !total number of possible rad outputs (types*specs)
      !there are 11 possible flux output 'species' but
      !only 8 possible optics output 'species'
      !for simplicity we take the largest and put up with
      !some redundancy (should be 88+72=160)
      INTEGER, PARAMETER :: PD72=187 ! Radiation (Ridley 10/2012)

      ! (eam, 2015): ISORROPIA pH, H+, & aerosol water
      INTEGER, PARAMETER :: PD73=3 
#endif

      !=================================================================
      ! Variables for printing out selected tracers in diagnostic output
      !=================================================================
#if   defined( RRTMG )
      INTEGER, PARAMETER :: MAX_DIAG   = 187
#else
      ! SDE 2013-11-17: Increased to 80 for UCX
      INTEGER, PARAMETER :: MAX_DIAG   = 80
#endif

      INTEGER, PARAMETER  :: MAX_TRACER = 600      ! Large placeholder value
!
! !PUBLIC DATA MEMBERS:
!
#if defined( BPCH_DIAG ) || defined( BPCH_TIMESER )
      !=================================================================
      ! Diagnostic counters & time variables
      !=================================================================
      INTEGER :: KDA48

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
      INTEGER :: LD09, LD10, LD31, LD52, LD19, LD57, LD59, LD60
      INTEGER :: LD61, LD62, LD73

      !=================================================================
      ! NDxx diagnostic flags
      !=================================================================
      INTEGER :: ND01, ND02, ND05, ND06, ND07, ND08, ND09, ND10, ND11
      INTEGER :: ND12, ND13, ND14, ND15, ND16, ND17, ND18, ND19, ND20
      INTEGER :: ND21, ND22, ND23, ND24, ND25, ND26, ND27, ND28, ND29
      INTEGER :: ND30, ND31, ND32, ND33, ND34, ND35, ND36, ND37, ND38
      INTEGER :: ND39, ND43, ND44, ND45, ND46, ND47
      INTEGER :: ND52, ND54, ND55, ND57, ND59
      INTEGER :: ND60, ND61, ND62, ND64, ND65, ND66, ND67, ND68
      INTEGER :: ND69, ND70, ND71, ND72, ND73, ND74, ND75
#endif

      !=================================================================
      ! Variables for printing out selected tracers in diagnostic output
      !=================================================================
      INTEGER :: TINDEX(MAX_DIAG,MAX_TRACER)          
      INTEGER :: TCOUNT(MAX_DIAG) 				  
      INTEGER :: TMAX(MAX_DIAG)					

#if defined( BPCH_DIAG ) || defined( BPCH_TIMESER )
      !=================================================================
      ! OH, J-Value, and 2-PM diagnostic arrays (bmy, 9/25/98)
      ! Move this here for now (bmy, 7/20/04)
      ! Removed HR1_NO and HR2_NO (mpayer, 11/8/13)
      !=================================================================
      REAL(fpp)  :: HR1_JV, HR2_JV
      REAL(fpp)  :: HR1_OH, HR2_OH, HR1_OTH, HR2_OTH
#endif
!
! !REMARKS:
!  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!  %%%  NOTE: THIS MODULE IS DEPRECATED.  MANY OF THE QUANTITIES STORED  %%%
!  %%%  in CMN_DIAG_MOD ARE NOW INCLUDED IN THE INPUT OPTIONS OBJECT,    %%%
!  %%%  BASED ON THE DERIVED TYPE IN Headers/input_opt_mod.F90.          %%%
!  %%%  (bmy, 11/19/12)                                                  %%%
!  %%%                                                                   %%%
!  %%%  ALSO NOTE: IN GEOS-CHEM v11-01, WE WILL BE RETIRING THE BPCH     %%%
!  %%%  DIAGNOSTICS.  MANY OF THESE PARAMETERS CAN THEN BE REMOVED.      %%%
!  %%%  (bmy, 6/24/16)                                                   %%%
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
!  02 Dec 2014 - M. Yannetti - Added PRECISION_MOD
!  15 Dec 2014 - M. Sulprizio- Moved radiation diagnostic from ND72 to ND72 to
!                              avoid conflicts with hourly max ppbv diagnostic.
!  23 Jun 2016 - R. Yantosca - Remove references to APM code; it is no longer
!                              compatible with the FlexChem implementation
!  21 Jul 2016 - R. Yantosca - Remove instances of NNPAR & obsolete parameters
!  19 Oct 2018 - R. Yantosca - Remove ND58 and LD58, they are obsolete
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
! !IROUTINE: Init_Cmn_Diag
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
!EOP
!------------------------------------------------------------------------------
!BOC
      ! Assume success
      RC   = GC_SUCCESS
      
#if defined( BPCH_DIAG ) || defined( BPCH_TIMESER )
      ! Set PD65 (for backwards compatibility, for now)
      PD65 = LLPAR*MAXFAM 
      PD64 = W_
#endif

      END SUBROUTINE Init_CMN_DIAG
!EOC
      END MODULE CMN_DIAG_MOD


