#if defined (ESMF_)
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: gc_idxmap_mod
!
! !DESCRIPTION: Module GIGC\_IDXMAP\_MOD contains the routine to initialize 
!  the mapping between ESMF tracer indexes and GEOS-Chem internal indexes.
!\\
!\\
! !INTERFACE: 
!      
MODULE GIGC_IdxMap_Mod
!
! !USES:
!      
  IMPLICIT NONE
  PUBLIC
!
! !REMARKS:
!  This module appears to be not used in the current implementation of
!  the GEOS-Chem ESMF interface.  It may just be for the BCC interface.
! 
! !REVISION HISTORY: 
!  16 Oct 2012 - M. Long     - Initial version
!  16 Oct 2012 - R. Yantosca - Added ProTeX headers
!  22 Oct 2012 - R. Yantosca - Renamed to gigc_idxmap_mod.f90
!EOP
!------------------------------------------------------------------------------
!BOC
CONTAINS
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: gigc_idxmap
!
! !DESCRIPTION: 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE GIGC_IdxMap( ID_TRACERS )
!
! !USES:
!
    USE GC_TYPE_MOD          
    USE TRACERID_MOD
        
    IMPLICIT NONE
!
! !INPUT PARAMETERS: 
!
    TYPE(ID_TRAC), INTENT(IN) :: ID_TRACERS
!
! !REMARKS:
!  This routine appears to be not used in the current implementation of
!  the GEOS-Chem ESMF interface.
! 
! !REVISION HISTORY: 
!  16 Oct 2012 - M. Long     - Initial version
!  16 Oct 2012 - R. Yantosca - Added ProTeX headers
!  22 Oct 2012 - R. Yantosca - Renamed to GIGC_IdxMap
!EOP
!------------------------------------------------------------------------------
!BOC
    
    IDTACET  = ID_TRACERS%ACET
    IDTALD2  = ID_TRACERS%ALD2
    IDTALK4  = ID_TRACERS%ALK4
    IDTBCPI  = ID_TRACERS%BCPI
    IDTBCPO  = ID_TRACERS%BCPO
    IDTC2H6  = ID_TRACERS%C2H6
    IDTC3H8  = ID_TRACERS%C3H8
    IDTCH2O  = ID_TRACERS%CH2O
    IDTCO    = ID_TRACERS%CO
    IDTDST1  = ID_TRACERS%DST1
    IDTDST2  = ID_TRACERS%DST2
    IDTDST3  = ID_TRACERS%DST3
    IDTDST4  = ID_TRACERS%DST4
    IDTHNO3  = ID_TRACERS%HNO3
    IDTHNO4  = ID_TRACERS%HNO4
    IDTISOP  = ID_TRACERS%ISOP
    IDTMACR  = ID_TRACERS%MACR
    IDTMEK   = ID_TRACERS%MEK 
    IDTMP    = ID_TRACERS%MP
    IDTMVK   = ID_TRACERS%MVK
    IDTNH3   = ID_TRACERS%NH3
    IDTNH4   = ID_TRACERS%NH4
    IDTNIT   = ID_TRACERS%NIT
    IDTNOx   = ID_TRACERS%NOx
    IDTOCPI  = ID_TRACERS%BCPI
    IDTOCPI  = ID_TRACERS%OCPI
    IDTOCPO  = ID_TRACERS%OCPO
    IDTOx    = ID_TRACERS%Ox
    IDTPMN   = ID_TRACERS%PMN
    IDTPRPE  = ID_TRACERS%PRPE
    IDTR4N2  = ID_TRACERS%R4N2
    IDTRCHO  = ID_TRACERS%RCHO
    IDTSALA  = ID_TRACERS%SALA
    IDTSALC  = ID_TRACERS%SALC
    IDTSO4   = ID_TRACERS%SO4
    IDTSOA1  = ID_TRACERS%SOA1
    IDTSOA2  = ID_TRACERS%SOA2
    IDTSOA3  = ID_TRACERS%SOA3
    IDTSOA4  = ID_TRACERS%SOA4
    IDTALCO  = ID_TRACERS%ALCO
    IDTALPH  = ID_TRACERS%ALPH
    IDTDMS   = ID_TRACERS%DMS
    IDTH2O2  = ID_TRACERS%H2O2
    IDTLIMO  = ID_TRACERS%LIMO
    IDTMEK   = ID_TRACERS%MEK
    IDTMSA   = ID_TRACERS%MSA
    IDTN2O5  = ID_TRACERS%N2O5
    IDTNITs  = ID_TRACERS%NITs
    IDTPAN   = ID_TRACERS%PAN
    IDTPPN   = ID_TRACERS%PPN
    IDTSO2   = ID_TRACERS%SO2
    IDTSO4s  = ID_TRACERS%SO4s
    IDTSOG1  = ID_TRACERS%SOG1
    IDTSOG2  = ID_TRACERS%SOG2
    IDTSOG3  = ID_TRACERS%SOG3
    IDTSOG4  = ID_TRACERS%SOG4
        
  END SUBROUTINE GIGC_IdxMap
!EOC
END MODULE GIGC_IdxMap_Mod
#endif
