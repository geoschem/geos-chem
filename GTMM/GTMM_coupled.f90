SUBROUTINE GTMM_coupled(year, month,  DD_Hg0, DD_HgII, WD_HgII, &
                        TS,   PREACC, RADSWG, Hg0reemit  )

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
!Translated into Matlab and accounted for fires by Guido van 
!der Werf.  
!See: van der Werf, G.R., J.T. Randerson, G.J. Collatz and L. 
!     Giglio, 2003.  Carbon emissions from fires in tropical
!     and subtropical ecosystems.  Global Change Biology 9
!     (4) 547-562.
!
!Translated into Fortran90 and added Mercury simulation by 
!Nicole Smith Downey - 2006
!
      
!<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  USE      defineConstants ! modify defineConstants.f90 to choose
                         ! parameters for your run 

  USE      loadCASAinput
  
  USE      defineArrays

  USE      DORESTART_MOD     ! New module to save/read Hg data for 
                             ! GEOS-CHEM. (ccc, 11/3/09)   

  USE      INPUT_GTMM_MOD
  
!<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  IMPLICIT NONE
  
  INTEGER, INTENT(INOUT)       :: year, month
  REAL*8, INTENT(IN)           :: DD_Hg0(72, 46), DD_HgII(72, 46), &
                                  WD_HgII(72, 46)    ! Hg deposition info.

  REAL*8, INTENT(IN),  DIMENSION(72, 46)  :: TS, PREACC, RADSWG !Met field info
  
  REAL*8, INTENT(OUT), DIMENSION(72, 46)  :: Hg0reemit  ! Reemitted flux, 
                                                        ! output to GEOS-Chem

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
