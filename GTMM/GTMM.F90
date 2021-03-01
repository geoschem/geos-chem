!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !PROGRAM: GTMM
!
! !DESCRIPTION: Based on the CASA (Carnegie, Ames, Stanford Approach) 
!  terrestrial biogeochemical model designed to simulate the terrestrial
!  carbon cycle using satellite data
!\\
!\\
! !INTERFACE:
!
PROGRAM GTMM
!
! !USES:
!
  USE      defineConstants ! modify defineConstants.f90 to choose
                         ! parameters for your run 

  USE      loadCASAinput
  
  USE      defineArrays

  USE      DORESTART_MOD ! New module to save/read Hg data and 
                         ! CASA data for continuation runs and GEOS-CHEM. 
                         ! (ccc, 11/3/09)

  USE      INPUT_GTMM_MOD

  USE      PRECISION_MOD    ! For GEOS-Chem Precision (fp)

  IMPLICIT NONE
!
! !AUTHOR:
!
!GTMM (Global Terrestrial Mercury Model) Developed by 
!Nicole Smith-Downey (nicolevdowney@gmail.com) 2006-2009
!See Smith-Downey, Sunderland and Jacob, JGR Biogeosciences, 2009
!
!Based on the CASA (Carnegie, Ames, Stanford Approach) terrestrial
!biogeochemical model designed to simulate the terrestrial
!carbon cycle using satellite data
!
!Original program written by Potter and Randerson
!See: Potter, C.S., J.T. Randerson, C.B. Field, P.A. Matson, 
!     P.M.Vitousek, H.A. Mooney, and S.A. Klooster, 1993.  
!     Terrestrial ecosystem production: A process model
!     based on satellite and surface data.  Global 
!     Biogeochemical Cycles (7) 811-841.
!
! !REVISION HISTORY:
!
! ( 1) Translated into Matlab and accounted for fires by Guido van 
!      der Werf.  
!See:  van der Werf, G.R., J.T. Randerson, G.J. Collatz and L. 
!      Giglio, 2003.  Carbon emissions from fires in tropical
!      and subtropical ecosystems.  Global Change Biology 9
!      (4) 547-562.
!
! ( 2) Translated into Fortran90 and added Mercury simulation by 
!      Nicole Smith Downey - 2006
!
! ( 3) Main program for offline simulations. Added coupling to GEOS-Chem
!      (see GTMM_coupled.f90) (ccc, 7/9/10)
! ( 4) Added capacity to restart runs. (ccc, 7/9/10)
! 25 Nov 2014 - M. Yannetti - Added PRECISION_MOD
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
  REAL(fp), DIMENSION(72, 46)  :: Hg0reemit    ! Dummy array. Used for coupling
                                             ! with GC only
  INTEGER :: year, month, ageClass
  LOGICAL :: LCPLE=.FALSE.                   ! Logical to define if we run
                                             ! coupled to GEOS-Chem or not
                                             ! Leave LCPLE to .FALSE. for the 
                                             ! stand-alone model. (ccc, 11/2/09)

  LOGICAL :: FIRST=.TRUE.

  !----------------------------------------------------------------------------
  ! GTMM begins here !
  !----------------------------------------------------------------------------

  CALL makeCASAarrays               !subroutine in defineArrays

  CALL readCASAparam

  CALL initialize
!
  CALL READ_GTMM_INPUT_FILE         !read data from input.gtmm (ccc)
!

  IF ( LRESTART ) FIRST=.FALSE.

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
  print '(a,f5.3,a)', '   ', decompHgEff, ' = fraction Hg released during &
       &decomposition'
  print '(a)', '   '
!<<<<<END VERIFY CONSTANTS>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  CALL load_data(LCPLE)    !subroutine in loadCASAinput.f90
  CALL CONV_TO_1D

  IF (.NOT.LRESTART) THEN   ! NPP equilibrium only done for the first run
     h=1
     print *, 'starting spinup to npp equilibirum'
     
     ! run NPP to equilibrium
     DO year=1,NPPequilibriumYear
        print *, yr
        DO month=1,12
!        CALL load_data(year, month, LCPLE)    !subroutine in loadCASAinput.f90
           CALL doLatitude
           CALL getSoilParams
           CALL getSoilMoistParams
           CALL getFuelWood
           
           CALL doPET
           CALL doSoilMoisture
           CALL doFPARandLAI(FIRST)
           CALL doOptimumTemperature
           CALL doNPP
           CALL doHerbivory
           CALL doLeafRootShedding
           CALL getFireParams
           IF (n_age_classes == 1) THEN
              CALL doTreeCarbon
              CALL doHerbCarbon
           ELSE
              DO ageClass=1,n_age_classes
                 CALL getAgeClassBF
                 CALL assignAgeClassToRunningPool
                 CALL doTreeCarbon
                 CALL doHerbCarbon
                 CALL assignRanPoolToAgeClass
                 age_class=age_class+1
              END DO
              CALL organizeAgeClasses
           ENDIF
           !CALL processData
           age_class=1
           mo=mo+1
        END DO
        mo=1
        yr=yr+1
        h=h+1
        IF (h > 25) THEN
           h=1
        ENDIF
        FIRST=.FALSE.
     END DO
     
     ! Save values to restart file
     CALL doSaveCASAforRestart
     
     
  ENDIF
  
  IF ( LRESTART ) THEN     ! Need to load data from previous run
     CALL doReadCASAfromRestart   
     
     ! Read Hg data saved from previous run. (ccc, 11/3/09)
     CALL doReadHgforGC
  ENDIF
  
  
  ! run terrestrial mercury
  age_class=1
  mo = 1
  DO year=(NPPequilibriumYear+1),(HgPoolsequilibriumYear)
     print *, yr
     DO month=1,12
!        CALL load_data(year, month, LCPLE)    !subroutine in loadCASAinput.f90
        CALL getSoilParams
        CALL getSoilMoistParams
        CALL doLatitude
        CALL getFuelWood
        
        CALL doMaxHg
        CALL loadHgDeposition(LCPLE)
        CALL doHgDeposition(LCPLE)
        IF (n_age_classes == 1) THEN
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
        IF (yr == indYear) THEN
           CALL HgOutForGEOS(LCPLE, Hg0reemit) 
        ENDIF
        !CALL processData
        age_class=1
        mo=mo+1
     END DO
     mo=1
     yr=yr+1
     h=h+1
     IF (h > 25) THEN
        h=1
     ENDIF
  END DO
  
  ! Save Hg data at equilibrium to use in GEOS-Chem. (ccc, 11/3/09)
  CALL doSaveHgforGC
  
  ! END OF RUN -- DEALLOCATE ALL
  print *, 'deallocating all arrays'
  CALL CleanupCASAarrays
END PROGRAM GTMM
!EOC
