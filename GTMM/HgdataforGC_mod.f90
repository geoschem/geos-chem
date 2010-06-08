!-----------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: HgdataforGC\_mod
!
! !DESCRIPTION: Module HgdataforGC\_mod contains subroutine to save annd read Hg
! pools at equilibrium to use when GTMM is coupled to GC. (ccc, 11/3/09)
!
! !INTERFACE:
!
MODULE HgdataforGC_mod
!
! !USES:
!
  IMPLICIT NONE
!
! !REVISION HISTORY:
!
!EOP
!-----------------------------------------------------------------------------

CONTAINS

!-----------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: doSaveHgforGC
!
! !DESCRIPTION: Subroutine doSaveHgforGC saves Hg pools at equilibrium to be
! used when GTMM is coupled to GEOS-Chem. (ccc, 11/3/09)
!
! !INTERFACE:
!
  SUBROUTINE doSaveHgforGC
!
! !USES:
!
    USE defineConstants
    USE defineArrays
!
! !REVISION HISTORY:
!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES
!
    CHARACTER(len=250) :: FILENAME

    ! Open binary file to write data
    FILENAME = OUTPUTPATH // 'HgPools'
    OPEN(UNIT=20, file=FILENAME, FORM="UNFORMATTED")

    ! Write data to FILENAME
    WRITE(20) hsurfstrpool_Hg
    WRITE(20) hsurfmetpool_Hg
    WRITE(20) hsurfmicpool_Hg
    WRITE(20) hsoilstrpool_Hg
    WRITE(20) hsoilmetpool_Hg
    WRITE(20) hsoilmicpool_Hg
    WRITE(20) hslowpool_Hg
    WRITE(20) harmoredpool_Hg

    WRITE(20) surfstrpool_Hg
    WRITE(20) surfmetpool_Hg
    WRITE(20) surfmicpool_Hg
    WRITE(20) soilstrpool_Hg
    WRITE(20) soilmetpool_Hg
    WRITE(20) soilmicpool_Hg
    WRITE(20) slowpool_Hg
    WRITE(20) armoredpool_Hg
              
    WRITE(20) Hg0_surf_soil
    WRITE(20) HgII_surf_soil
    WRITE(20) hleafpool_Hg
    WRITE(20) leafpool_Hg
    
    WRITE(20) leafpool
    WRITE(20) cwdpool
    WRITE(20) abovewoodpool
    WRITE(20) belowwoodpool
    WRITE(20) frootpool
    WRITE(20) surfstrpool
    WRITE(20) surfmetpool
    WRITE(20) surfmicpool
    WRITE(20) soilstrpool
    WRITE(20) soilmetpool
    WRITE(20) soilmicpool
    WRITE(20) slowpool
    WRITE(20) armoredpool

    WRITE(20) hleafpool
    WRITE(20) hfrootpool
    WRITE(20) hsurfstrpool
    WRITE(20) hsurfmetpool
    WRITE(20) hsurfmicpool
    WRITE(20) hsoilstrpool
    WRITE(20) hsoilmetpool
    WRITE(20) hsoilmicpool
    WRITE(20) hslowpool
    WRITE(20) harmoredpool

    ! Close file
    CLOSE(20)

  END SUBROUTINE doSaveHgforGC
!EOC
!-----------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: doReadHgforGC
!
! !DESCRIPTION: Subroutine doReadHgforGC reads Hg pools at equilibrium to
! use them in GTMM when coupled to GEOS-Chem. (ccc, 11/3/09)
!
! !INTERFACE:
!
  SUBROUTINE doReadHgforGC
!
! !USES:
!
    USE defineConstants
    USE defineArrays
!
! !REVISION HISTORY:
!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES
!
    CHARACTER(len=250) :: FILENAME

    ! Open binary file to read data
    FILENAME = OUTPUTPATH // 'HgPools'
    OPEN(UNIT=20, file=FILENAME, STATUS="OLD", FORM="UNFORMATTED")

    ! Read data to FILENAME
    READ(20) hsurfstrpool_Hg
    READ(20) hsurfmetpool_Hg
    READ(20) hsurfmicpool_Hg
    READ(20) hsoilstrpool_Hg
    READ(20) hsoilmetpool_Hg
    READ(20) hsoilmicpool_Hg
    READ(20) hslowpool_Hg
    READ(20) harmoredpool_Hg
    
    READ(20) surfstrpool_Hg
    READ(20) surfmetpool_Hg
    READ(20) surfmicpool_Hg
    READ(20) soilstrpool_Hg
    READ(20) soilmetpool_Hg
    READ(20) soilmicpool_Hg
    READ(20) slowpool_Hg
    READ(20) armoredpool_Hg
         
    READ(20) Hg0_surf_soil
    READ(20) HgII_surf_soil
    READ(20) hleafpool_Hg
    READ(20) leafpool_Hg
    
    READ(20) leafpool
    READ(20) cwdpool
    READ(20) abovewoodpool
    READ(20) belowwoodpool
    READ(20) frootpool
    READ(20) surfstrpool
    READ(20) surfmetpool
    READ(20) surfmicpool
    READ(20) soilstrpool
    READ(20) soilmetpool
    READ(20) soilmicpool
    READ(20) slowpool
    READ(20) armoredpool

    READ(20) hleafpool
    READ(20) hfrootpool
    READ(20) hsurfstrpool
    READ(20) hsurfmetpool
    READ(20) hsurfmicpool
    READ(20) hsoilstrpool
    READ(20) hsoilmetpool
    READ(20) hsoilmicpool
    READ(20) hslowpool
    READ(20) harmoredpool

    ! Close file
    CLOSE(20)

  END SUBROUTINE doReadHgforGC

END MODULE HgdataforGC_mod
!EOC
