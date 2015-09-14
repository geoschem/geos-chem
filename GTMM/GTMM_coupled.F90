!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: GTMM_coupled
!
! !DESCRIPTION: Main subroutine for GTMM when coupled to GEOS-Chem. Replace
!  GTMM.f90.
!\\
!\\
! !INTERFACE:
!
SUBROUTINE GTMM_coupled(year, month,  DD_Hg0, DD_HgII, WD_HgII, &
                        TS,   PREACC, RADSWG, Hg0reemit  )
!
! !USES:
!
  USE      defineConstants ! modify defineConstants.f90 to choose
                         ! parameters for your run 

  USE      loadCASAinput
  
  USE      defineArrays

  USE      DORESTART_MOD     ! New module to save/read Hg data for 
                             ! GEOS-CHEM. (ccc, 11/3/09)   

  USE      INPUT_GTMM_MOD
  
  USE      PRECISION_MOD    ! For GEOS-Chem Precision (fp)

  IMPLICIT NONE
!
! !INPUT PARAMETERS:
!  
  REAL(fp), INTENT(IN)           :: DD_Hg0(72, 46), DD_HgII(72, 46), &
                                  WD_HgII(72, 46)    ! Hg deposition info.

  REAL(fp), INTENT(IN),  DIMENSION(72, 46)  :: TS, PREACC, RADSWG !Met field info
!
! !INPUT/OUTPUT PARAMETERS:
!  
  INTEGER, INTENT(INOUT)       :: year, month
!
! !OUTPUT PARAMETERS:
!  
  REAL(fp), INTENT(OUT), DIMENSION(72, 46)  :: Hg0reemit  ! Reemitted flux, 
                                                        ! output to GEOS-Chem
!
! !REVISION HISTORY:
!  09 Jul 2010 - C. Carouge  - First version. Adapted from GTMM.f90
!  25 Nov 2014 - M. Yannetti - Added PRECISION_MOD
!EOP
!------------------------------------------------------------------------------
!BOC
! 
! !LOCAL VARIABLES:
!
  INTEGER :: ageClass

  ! Define if we run CASA for carbon equilibre or not.
  LOGICAL, SAVE     :: LEQUIL = .TRUE.

  LOGICAL           :: LCPLE

  ! Give value to LCPLE 
  LCPLE=.TRUE.

  CALL readCASAparam
  CALL makeCASAarrays               !subroutine in defineArrays

  CALL initialize

!
  CALL READ_GTMM_INPUT_FILE         !read data from input.gtmm (ccc)
!
!<<<<<VERIFY CONSTANTS>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  print '(a)', ''
  print '(a)', 'STARTING GTMM...'
  print '(a)', ''
  print '(a)', 'Your outputpath is: ', outputpath
  print '(a)', '   '
  print '(a)', 'You have chosen: ' 
  print '(a, i4, a, i4, a)', '   ', rows, ' rows and '&
       &, columns, ' columns'
  print '(a, i6, a)', '   ', NPPequilibriumYear, ' years to&
       & equilibrate Carbon Pools'
  print '(a, i6, a)', '   ', HgPoolsequilibriumYear, ' years to&
       & equilibrate carbon/Hg pools'
  print '(a, i6, a)', '   ', n_age_classes, ' age classes'
  print '(a,f3.1,a)', '   ', Q10, ' = Q10'
  print '(a,f3.1,a)', '   ', EMAX, ' = EMAX'
  print '(a,f3.1,a)', '   ', aboveWoodFraction, ' = aboveWoodFraction'
  print '(a,f3.1,a)', '   ', herbivoreEff, ' = herbivore efficiency'
  print '(a,f5.3,a)', '   ', decompHgEff, ' = fraction Hg released during decomposition'
  print '(a)', '   '
!<<<<<END VERIFY CONSTANTS>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  CALL doReadCASAfromRestart   

  h=1
     
  ! run terrestrial mercury for the current month only
  age_class=1
     
  yr = year
  mo = month
  
  print *, year
!  CALL load_data(year, month, LCPLE, TS, PREACC, RADSWG) !in loadCASAinput.f90
  CALL load_data(LCPLE) !in loadCASAinput.f90
  CALL load_GC_data(month, TS, PREACC, RADSWG)
  CALL CONV_TO_1D

  CALL getSoilParams
  CALL getSoilMoistParams
  CALL doLatitude
  CALL getFuelWood
  
  ! Read Hg data saved from equilibrium run. (ccc, 11/3/09)
  CALL doReadHgforGC

  CALL doMaxHg
  ! We need to pass deposition arrays read by GEOS-Chem (ccc, 9/17/09)
  CALL loadHgDeposition(LCPLE, DD_Hg0, DD_HgII, WD_HgII)   
  CALL doHgDeposition(LCPLE)
  IF (n_age_classes .eq. 1) THEN
     CALL doTreeCarbonHg(LCPLE)
     CALL doHerbCarbonHg(LCPLE)
  ELSE
     DO ageClass=1, n_age_classes
        CALL getAgeClassBF
        CALL assignAgeClassToRunningPool
        CALL doTreeCarbonHg(LCPLE)
        CALL doHerbCarbonHg(LCPLE)
        CALL assignRanPoolToAgeClass
        age_class=age_class+1
     END DO
     CALL organizeAgeClasses
  ENDIF
  CALL HgOutForGEOS(LCPLE, Hg0reemit) 
  
  
  ! END OF RUN -- DEALLOCATE ALL
  print *, 'deallocating all arrays'
  CALL CleanupCASAarrays

END SUBROUTINE GTMM_coupled
!EOC
