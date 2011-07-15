#if defined (ESMF_TESTBED_)
!------------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 910.1 and      !
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: bmy_GC_Value_Mod
!
! !DESCRIPTION: This module contains values for a single column that were
!  saved out from a GEOS-Chem column code simulation.  Using saved values 
!  facilitates the unit testing.
!\\
!\\
! !INTERFACE:
!
MODULE TESTBED_VALUE_MOD
!
! !USES:
!
  IMPLICIT NONE
  PUBLIC
!
! !REVISION HISTORY:
!  06 Dec 2009 - A. da Silva - Created the GEOSCHEM skeleton.
!  07 Apr 2010 - R. Yantosca - Updated comments, cosmetic changes 
!  12 Apr 2010 - R. Yantosca - Corrected wrong value for GC_SST
!  13 Apr 2010 - R. Yantosca - GC_PEDGE and GC_PMID were in hPa; corrected
!                              them so that they are now have units of Pa.
!  11 May 2010 - R. Yantosca - Remove the PARAMETER declarations for arrays
!                              that hold GEOS-Chem information.  This data
!                              is now read from disk.
!  04 Jun 2010 - R. Yantosca - Rename GC_PBLH to GC_ZPBL
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS: 
!
    INTEGER, PARAMETER :: LMAX = 72

    !=======================================================================
    ! Import state fields: Surface met fields
    !=======================================================================

    ! Met fields
    REAL :: GC_ALBD      
    REAL :: GC_AREA_M2   
    REAL :: GC_CLDFRC    
    REAL :: GC_FRCLND    
    REAL :: GC_GWETTOP   
    REAL :: GC_HFLUX     
    REAL :: GC_LWI       
    REAL :: GC_PARDR     
    REAL :: GC_PARDF     
    REAL :: GC_ZPBL   
    REAL :: GC_PRECCON   
    REAL :: GC_PRECTOT   
    REAL :: GC_RADSWG    
    REAL :: GC_SST       
    REAL :: GC_COSZ      
    REAL :: GC_TROPP     
    REAL :: GC_TO3       
    REAL :: GC_TS        
    REAL :: GC_U10M      
    REAL :: GC_USTAR     
    REAL :: GC_UVALBEDO  
    REAL :: GC_V10M      
    REAL :: GC_Z0        

    ! Landtype / Leaf area indices
    REAL :: GC_IREG      
    REAL :: GC_ILAND01   
    REAL :: GC_ILAND02   
    REAL :: GC_ILAND03   
    REAL :: GC_ILAND04   
    REAL :: GC_ILAND05   
    REAL :: GC_ILAND06   
    REAL :: GC_ILAND07   
    REAL :: GC_ILAND08   
    REAL :: GC_ILAND09   
    REAL :: GC_ILAND10   
    REAL :: GC_ILAND11   
    REAL :: GC_ILAND12   
    REAL :: GC_ILAND13   
    REAL :: GC_ILAND14   
    REAL :: GC_ILAND15   
    REAL :: GC_IUSE01    
    REAL :: GC_IUSE02    
    REAL :: GC_IUSE03    
    REAL :: GC_IUSE04    
    REAL :: GC_IUSE05    
    REAL :: GC_IUSE06    
    REAL :: GC_IUSE07    
    REAL :: GC_IUSE08    
    REAL :: GC_IUSE09    
    REAL :: GC_IUSE10    
    REAL :: GC_IUSE11    
    REAL :: GC_IUSE12    
    REAL :: GC_IUSE13    
    REAL :: GC_IUSE14    
    REAL :: GC_IUSE15    
    REAL :: GC_LAI01      
    REAL :: GC_LAI02      
    REAL :: GC_LAI03     
    REAL :: GC_LAI04      
    REAL :: GC_LAI05     
    REAL :: GC_LAI06      
    REAL :: GC_LAI07     
    REAL :: GC_LAI08     
    REAL :: GC_LAI09     
    REAL :: GC_LAI10     
    REAL :: GC_LAI11     
    REAL :: GC_LAI12     
    REAL :: GC_LAI13     
    REAL :: GC_LAI14     
    REAL :: GC_LAI15       

    !=======================================================================
    ! Import state fields: Column met fields
    !=======================================================================
    REAL :: GC_AD          (LMAX  )  
    REAL :: GC_AIRDENS     (LMAX  )  
    REAL :: GC_AIRVOL      (LMAX  )        
    REAL :: GC_BXHEIGHT    (LMAX  )       
    REAL :: GC_CLDF        (LMAX  )           
    REAL :: GC_CMFMC       (LMAX+1)        
    REAL :: GC_DELP        (LMAX  )         
    REAL :: GC_DQIDTMST    (LMAX  )     
    REAL :: GC_DQLDTMST    (LMAX  )     
    REAL :: GC_DQVDTMST    (LMAX  )     
    REAL :: GC_DTRAIN      (LMAX  )       
    REAL :: GC_MOISTQ      (LMAX  )       
    REAL :: GC_OPTD        (LMAX  )
    REAL :: GC_PEDGE       (LMAX+1)        
    REAL :: GC_PMID        (LMAX  )  
    REAL :: GC_RH          (LMAX  )           
    REAL :: GC_Q           (LMAX  )            
    REAL :: GC_T           (LMAX  )            
    REAL :: GC_TAUCLI      (LMAX  )       
    REAL :: GC_TAUCLW      (LMAX  )       

    !=======================================================================
    ! Import state fields: Fields of the EMISSIONS_1d array
    !=======================================================================
    REAL :: GC_EMISS_NOx   (LMAX  )   
    REAL :: GC_EMISS_O3    (LMAX  )    
    REAL :: GC_EMISS_CO    (LMAX  )    
    REAL :: GC_EMISS_ALK4  (LMAX  )  
    REAL :: GC_EMISS_ISOP  (LMAX  )  
    REAL :: GC_EMISS_HNO3  (LMAX  )  
    REAL :: GC_EMISS_ACET  (LMAX  )  
    REAL :: GC_EMISS_MEK   (LMAX  )   
    REAL :: GC_EMISS_ALD2  (LMAX  )  
    REAL :: GC_EMISS_PRPE  (LMAX  )  
    REAL :: GC_EMISS_C3H8  (LMAX  )  
    REAL :: GC_EMISS_CH2O  (LMAX  )  
    REAL :: GC_EMISS_C2H6  (LMAX  )  

    !=======================================================================
    ! Import state fields: Parameters for SCHEM strat chemistry
    !=======================================================================
    REAL :: GC_SOX_OH      (LMAX  )       
    REAL :: GC_SOX_PCO     (LMAX  )      
    REAL :: GC_SOX_LCO     (LMAX  )      
    REAL :: GC_SOX_JV_NOX  (LMAX  )   
    REAL :: GC_SOX_JV_H2O2 (LMAX  )  
    REAL :: GC_SOX_JV_ACET (LMAX  )  
    REAL :: GC_SOX_JV_MEK  (LMAX  )   
    REAL :: GC_SOX_JV_ALD2 (LMAX  )  
    REAL :: GC_SOX_JV_RCHO (LMAX  )  
    REAL :: GC_SOX_JV_MVK  (LMAX  )   
    REAL :: GC_SOX_JV_MACR (LMAX  )  
    REAL :: GC_SOX_JV_R4N2 (LMAX  )  
    REAL :: GC_SOX_JV_CH2O (LMAX  )  
    REAL :: GC_SOX_JV_N2O5 (LMAX  )  
    REAL :: GC_SOX_JV_HNO4 (LMAX  )  
    REAL :: GC_SOX_JV_MP   (LMAX  )    

    !=======================================================================
    ! Internal state fields: Initial concentrations of chemical species
    !=======================================================================
    REAL :: GC_A3O2        (LMAX  )
    REAL :: GC_ACET        (LMAX  ) 
    REAL :: GC_ALD2        (LMAX  ) 
    REAL :: GC_ALK4        (LMAX  ) 
    REAL :: GC_ATO2        (LMAX  ) 
    REAL :: GC_B3O2        (LMAX  ) 
    REAL :: GC_C2H6        (LMAX  ) 
    REAL :: GC_C3H8        (LMAX  )
    REAL :: GC_CH2O        (LMAX  )
    REAL :: GC_CO          (LMAX  ) 
    REAL :: GC_DRYCH2O     (LMAX  )
    REAL :: GC_DRYH2O2     (LMAX  )
    REAL :: GC_DRYHNO3     (LMAX  )
    REAL :: GC_DRYN2O5     (LMAX  )
    REAL :: GC_DRYNO2      (LMAX  )
    REAL :: GC_DRYO3       (LMAX  )
    REAL :: GC_DRYPAN      (LMAX  )
    REAL :: GC_DRYPMN      (LMAX  )
    REAL :: GC_DRYPPN      (LMAX  )
    REAL :: GC_DRYR4N2     (LMAX  )
    REAL :: GC_ETO2        (LMAX  )
    REAL :: GC_ETP         (LMAX  )
    REAL :: GC_GCO3        (LMAX  )
    REAL :: GC_GLYC        (LMAX  )
    REAL :: GC_GP          (LMAX  )
    REAL :: GC_GPAN        (LMAX  )
    REAL :: GC_H2O2        (LMAX  )
    REAL :: GC_HAC         (LMAX  )
    REAL :: GC_HNO2        (LMAX  )
    REAL :: GC_HNO3        (LMAX  )
    REAL :: GC_HNO4        (LMAX  )
    REAL :: GC_HO2         (LMAX  )
    REAL :: GC_IALD        (LMAX  )
    REAL :: GC_IAO2        (LMAX  )
    REAL :: GC_IAP         (LMAX  )
    REAL :: GC_INO2        (LMAX  )
    REAL :: GC_INPN        (LMAX  )
    REAL :: GC_ISN1        (LMAX  )
    REAL :: GC_ISNP        (LMAX  )
    REAL :: GC_ISOP        (LMAX  )
    REAL :: GC_KO2         (LMAX  )
    REAL :: GC_MACR        (LMAX  )
    REAL :: GC_MAN2        (LMAX  )
    REAL :: GC_MAO3        (LMAX  )
    REAL :: GC_MAOP        (LMAX  )
    REAL :: GC_MAP         (LMAX  )
    REAL :: GC_MCO3        (LMAX  )
    REAL :: GC_MEK         (LMAX  )
    REAL :: GC_MGLY        (LMAX  )
    REAL :: GC_MO2         (LMAX  )
    REAL :: GC_MP          (LMAX  )
    REAL :: GC_MRO2        (LMAX  )
    REAL :: GC_MRP         (LMAX  )
    REAL :: GC_MVK         (LMAX  )
    REAL :: GC_MVN2        (LMAX  )
    REAL :: GC_N2O5        (LMAX  )
    REAL :: GC_NO          (LMAX  )
    REAL :: GC_NO2         (LMAX  )
    REAL :: GC_NO3         (LMAX  )
    REAL :: GC_O3          (LMAX  )
    REAL :: GC_OH          (LMAX  )
    REAL :: GC_PAN         (LMAX  )
    REAL :: GC_PMN         (LMAX  )
    REAL :: GC_PO2         (LMAX  )
    REAL :: GC_PP          (LMAX  )
    REAL :: GC_PPN         (LMAX  )
    REAL :: GC_PRN1        (LMAX  )
    REAL :: GC_PRPE        (LMAX  )
    REAL :: GC_PRPN        (LMAX  )
    REAL :: GC_R4N1        (LMAX  )
    REAL :: GC_R4N2        (LMAX  )
    REAL :: GC_R4O2        (LMAX  )
    REAL :: GC_R4P         (LMAX  )
    REAL :: GC_RA3P        (LMAX  )
    REAL :: GC_RB3P        (LMAX  )
    REAL :: GC_RCHO        (LMAX  )
    REAL :: GC_RCO3        (LMAX  )
    REAL :: GC_RIO1        (LMAX  )
    REAL :: GC_RIO2        (LMAX  )
    REAL :: GC_RIP         (LMAX  )
    REAL :: GC_RP          (LMAX  )
    REAL :: GC_VRO2        (LMAX  )
    REAL :: GC_VRP         (LMAX  )
    REAL :: GC_DMS         (LMAX  )
    REAL :: GC_SO2         (LMAX  )
    REAL :: GC_SO4         (LMAX  )
    REAL :: GC_MSA         (LMAX  )
    REAL :: GC_LISOPOH     (LMAX  )
    REAL :: GC_ROH         (LMAX  )
    REAL :: GC_RCOOH       (LMAX  )
    REAL :: GC_O2CH2OH     (LMAX  )   
    REAL :: GC_O2          (LMAX  )   
    REAL :: GC_O           (LMAX  )   
    REAL :: GC_NH3         (LMAX  )   
    REAL :: GC_NH2         (LMAX  )
    REAL :: GC_N2          (LMAX  )   
    REAL :: GC_MOH         (LMAX  )   
    REAL :: GC_MNO3        (LMAX  )   
    REAL :: GC_M           (LMAX  )   
    REAL :: GC_ISNO3       (LMAX  )   
    REAL :: GC_HCOOH       (LMAX  )   
    REAL :: GC_H2O         (LMAX  )   
    REAL :: GC_H2          (LMAX  )   
    REAL :: GC_H           (LMAX  )   
    REAL :: GC_GLYX        (LMAX  )   
    REAL :: GC_GLPAN       (LMAX  )   
    REAL :: GC_GLP         (LMAX  )   
    REAL :: GC_GLCO3       (LMAX  )   
    REAL :: GC_EOH         (LMAX  )   
    REAL :: GC_EMISSION    (LMAX  )   
    REAL :: GC_CH4         (LMAX  )   
    REAL :: GC_ACTA        (LMAX  )   

    !=======================================================================
    ! H2O2s and SO2s for wetdep
    !=======================================================================
    REAL :: GC_H2O2s       (LMAX  )
    REAL :: GC_SO2s        (LMAX  ) 

    !=======================================================================
    ! Secondary organic aerosol parameters
    !=======================================================================
    REAL :: GC_ORVC_TERP   (LMAX  )
    REAL :: GC_ORVC_SESQ   (LMAX  )

    !=======================================================================
    ! Advected tracers [mol/mol]
    !=======================================================================
    REAL :: GC_TRC_NOx     (LMAX  )
    REAL :: GC_TRC_Ox      (LMAX  )
    REAL :: GC_TRC_PAN     (LMAX  )
    REAL :: GC_TRC_CO      (LMAX  )
    REAL :: GC_TRC_ALK4    (LMAX  )
    REAL :: GC_TRC_ISOP    (LMAX  )
    REAL :: GC_TRC_HNO3    (LMAX  )
    REAL :: GC_TRC_H2O2    (LMAX  )
    REAL :: GC_TRC_ACET    (LMAX  )
    REAL :: GC_TRC_MEK     (LMAX  )
    REAL :: GC_TRC_ALD2    (LMAX  )
    REAL :: GC_TRC_RCHO    (LMAX  )
    REAL :: GC_TRC_MVK     (LMAX  )
    REAL :: GC_TRC_MACR    (LMAX  )
    REAL :: GC_TRC_PMN     (LMAX  )
    REAL :: GC_TRC_PPN     (LMAX  )
    REAL :: GC_TRC_R4N2    (LMAX  )
    REAL :: GC_TRC_PRPE    (LMAX  )
    REAL :: GC_TRC_C3H8    (LMAX  )
    REAL :: GC_TRC_CH2O    (LMAX  )
    REAL :: GC_TRC_C2H6    (LMAX  )
    REAL :: GC_TRC_N2O5    (LMAX  )
    REAL :: GC_TRC_HNO4    (LMAX  )
    REAL :: GC_TRC_MP      (LMAX  )
    REAL :: GC_TRC_DMS     (LMAX  )
    REAL :: GC_TRC_SO2     (LMAX  )
    REAL :: GC_TRC_SO4     (LMAX  )
    REAL :: GC_TRC_SO4s    (LMAX  )
    REAL :: GC_TRC_MSA     (LMAX  )
    REAL :: GC_TRC_NH3     (LMAX  )
    REAL :: GC_TRC_NH4     (LMAX  )
    REAL :: GC_TRC_NIT     (LMAX  )
    REAL :: GC_TRC_NITs    (LMAX  )
    REAL :: GC_TRC_BCPI    (LMAX  )
    REAL :: GC_TRC_OCPI    (LMAX  )
    REAL :: GC_TRC_BCPO    (LMAX  )
    REAL :: GC_TRC_OCPO    (LMAX  )
    REAL :: GC_TRC_ALPH    (LMAX  )
    REAL :: GC_TRC_LIMO    (LMAX  )
    REAL :: GC_TRC_ALCO    (LMAX  )
    REAL :: GC_TRC_SOG1    (LMAX  )
    REAL :: GC_TRC_SOG2    (LMAX  )
    REAL :: GC_TRC_SOG3    (LMAX  )
    REAL :: GC_TRC_SOG4    (LMAX  )
    REAL :: GC_TRC_SOA1    (LMAX  )
    REAL :: GC_TRC_SOA2    (LMAX  )
    REAL :: GC_TRC_SOA3    (LMAX  )
    REAL :: GC_TRC_SOA4    (LMAX  )
    REAL :: GC_TRC_DST1    (LMAX  )
    REAL :: GC_TRC_DST2    (LMAX  )
    REAL :: GC_TRC_DST3    (LMAX  )
    REAL :: GC_TRC_DST4    (LMAX  )
    REAL :: GC_TRC_SALA    (LMAX  )
    REAL :: GC_TRC_SALC    (LMAX  )

END MODULE TESTBED_VALUE_MOD
#endif
!EOC
