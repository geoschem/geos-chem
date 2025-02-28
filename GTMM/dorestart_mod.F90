!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: dorestart_mod
!
! !DESCRIPTION: Module DORESTART\_MOD contains subroutines to save and read 
!  data created by GTMM and used to restart a stand-alone run or a coupled
!  run with GEOS-Chem.
!
! !INTERFACE:
!
MODULE DORESTART_MOD
!
! !USES:
!
  IMPLICIT NONE
!
! !REVISION HISTORY:
!  16 Dec 2009 - C. Carouge   - Initial version
!EOP
!-----------------------------------------------------------------------------
!BOC
  CONTAINS
!EOC
!-----------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: doSaveHgforGC
!
! !DESCRIPTION: Subroutine doSaveHgforGC saves Hg pools at equilibrium to be
! used when GTMM is coupled to GEOS-Chem.
!\\
!\\
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
!  03 Nov 2009 - C. Carouge   - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES
!
    CHARACTER(len=250) :: FILENAME

    ! Open binary file to write data
    FILENAME = OUTPUTPATH // 'HgPools'
    OPEN(UNIT=20, file=FILENAME, STATUS="NEW", FORM="UNFORMATTED")

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
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: doReadHgforGC
!
! !DESCRIPTION: Subroutine doReadHgforGC reads Hg pools at equilibrium to
! use them in GTMM when coupled to GEOS-Chem.
!\\
!\\
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
!  03 Nov 2009 - C. Carouge   - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES
!
    OPEN(UNIT=20, file=RESTARTFILE, STATUS="OLD", FORM="UNFORMATTED")

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
!EOC
!-----------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: doSaveCASAforRestart
!
! !DESCRIPTION: Subroutine doSaveCASAforRestart saves CASA values at 
!  equilibrium to be used when an equilibrium run for GTMM is continued. 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE doSaveCASAforRestart
!
! !USES:
!
    USE defineConstants
    USE defineArrays
!
! !REVISION HISTORY:
!  16 Dec 2009 - C. Carouge   - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES
!
    CHARACTER(len=250) :: FILENAME

    ! Open binary file to write data
    FILENAME = OUTPUTPATH // 'restart'
    OPEN(UNIT=20, file=FILENAME, STATUS="NEW", FORM="UNFORMATTED")

    ! Write data to FILENAME
    WRITE(20) LAI
    WRITE(20) herb_seasonality
    WRITE(20) grass_herbivory
    WRITE(20) trees_herbivory
    WRITE(20) rootlitscalar
    WRITE(20) litterscalar
    WRITE(20) hlitterscalar
    WRITE(20) NPP
    WRITE(20) abiotic
    WRITE(20) lais
    WRITE(20) ccWood
    WRITE(20) ccLeaf
    WRITE(20) ccFineLitter
    WRITE(20) ccCwd
    WRITE(20) PET
    WRITE(20) CCratio_previous

    ! Only used if restart the run from the carbon equilibrium year. 
    ! If the mercury equilibrium run is started, use data saved in HgPools
    WRITE(20) leafpool
    WRITE(20) surfstrpool
    WRITE(20) surfmetpool
    WRITE(20) surfmicpool
    WRITE(20) soilstrpool
    WRITE(20) soilmetpool
    WRITE(20) soilmicpool
    WRITE(20) slowpool
    WRITE(20) armoredpool
    WRITE(20) hleafpool
    WRITE(20) hsurfstrpool
    WRITE(20) hsurfmetpool
    WRITE(20) hsurfmicpool
    WRITE(20) hsoilstrpool
    WRITE(20) hsoilmetpool
    WRITE(20) hsoilmicpool
    WRITE(20) hslowpool
    WRITE(20) harmoredpool
    WRITE(20) abovewoodpool
    WRITE(20) belowwoodpool
    WRITE(20) frootpool

    ! Close file
    CLOSE(20)

  END SUBROUTINE doSaveCASAforRestart
!EOC
!-----------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: doReadCASAfromRestart
!
! !DESCRIPTION: Subroutine doReadCASAfromRestart reads CASA values at 
!  equilibrium to be used when an equilibrium run for GTMM is continued. 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE doReadCASAfromRestart
!
! !USES:
!
    USE defineConstants
    USE defineArrays
!
! !REVISION HISTORY:
!  16 Dec 2009 - C. Carouge   - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES
!
    CHARACTER(len=250) :: FILENAME

    ! Open binary file to read data
    FILENAME = OUTPUTPATH // 'restart'
    OPEN(UNIT=20, file=FILENAME, STATUS="OLD", FORM="UNFORMATTED")

    ! Read data to FILENAME
    READ(20) LAI
    READ(20) herb_seasonality
    READ(20) grass_herbivory
    READ(20) trees_herbivory
!    READ(20) LTCON
!    READ(20) LTVARSUM
    READ(20) rootlitscalar
    READ(20) litterscalar
    READ(20) hlitterscalar
!    READ(20) AVELAI
    READ(20) NPP
    READ(20) abiotic
    READ(20) lais
!    READ(20) topt
    READ(20) ccWood
    READ(20) ccLeaf
    READ(20) ccFineLitter
    READ(20) ccCwd
!    READ(20) mortality_tree
    WRITE(20) PET
    WRITE(20) CCratio_previous
! Only used if restart from the carbon equilibrium year. Else use the data
! saved in HgPools
    WRITE(20) leafpool
    WRITE(20) surfstrpool
    WRITE(20) surfmetpool
    WRITE(20) surfmicpool
    WRITE(20) soilstrpool
    WRITE(20) soilmetpool
    WRITE(20) soilmicpool
    WRITE(20) slowpool
    WRITE(20) armoredpool
    WRITE(20) hleafpool
    WRITE(20) hsurfstrpool
    WRITE(20) hsurfmetpool
    WRITE(20) hsurfmicpool
    WRITE(20) hsoilstrpool
    WRITE(20) hsoilmetpool
    WRITE(20) hsoilmicpool
    WRITE(20) hslowpool
    WRITE(20) harmoredpool
    WRITE(20) abovewoodpool
    WRITE(20) belowwoodpool
    WRITE(20) frootpool

    ! Close file
    CLOSE(20)

  END SUBROUTINE doReadCASAfromRestart
!EOC
END MODULE DORESTART_MOD

