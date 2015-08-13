!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: doHgDeposition
!
! !DESCRIPTION: This subrouitine ...
!\\
!\\
! !INTERFACE:
!
SUBROUTINE doHgDeposition
!
! !USES:
!
  USE defineConstants
  USE loadCASAinput
  USE defineArrays

  USE PRECISION_MOD    ! For GEOS-Chem Precision (fp)
  
  implicit none
!
! !REVISION HISTORY:
!  09 Jul 2010 - C. Carouge  - Parallelization
!  25 Nov 2014 - M. Yannetti - Added PRECISION_MOD
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
  REAL(fp)     :: photo_frac(n_veg, 1)
      
!!!!!   Hg0dry
  !! 1 - deposited to leaf and soil surfaces
  !! 2 - directly incorporated into leaf tissue via stomates
  !! 3 - volatilized as a fxn of temperature

!$OMP PARALLEL       &
!$OMP DEFAULT(SHARED)
!$OMP WORKSHARE
  fstom(:,1)=(LAI(:,mo)/5)  ! chosen to match Rea et al seasonal cycle

  WHERE (fstom(:,1) > 1.0e+0_fp) 
     fstom(:,1)=1.0e+0_fp
  END WHERE

  fstom(:,1)=fstom(:,1)
  fsoil(:,1)=1.00e+0_fp-fstom(:,1)
  Hg0_surf_soil(:,1)=Hg0_surf_soil(:,1)+(Hg0dry(:,1)*fsoil(:,1))
  leafpool_hg(:,1)=leafpool_hg(:,1)+(Hg0dry(:,1)*fstom(:,1))
  hleafpool_hg(:,1)=hleafpool_hg(:,1)+(Hg0dry(:,1)*fstom(:,1))
  !all elemental hg sitting in surface pools is volatilized each
  !month
  freemitted(:,1)=1.0e+0_fp
  reemitted(:,1)=(freemitted(:,1)*Hg0_surf_soil(:,1))
  
  Hg0_surf_soil(:,1)=Hg0_surf_soil(:,1)-reemitted(:,1)
!!! HgII dry
  !! 1 - deposited to leaf and soil surfaces
  !! 2 - directly incorporated into leaf tissue via stomates
  !! 3 - photoreduced as fxn of T and light
  !! 4 - washed off leaf and soil surfaces and added to HgIIwet
  
  ! for 1 same fleaf, fsoil and fstom as for Hg0dry
  HgII_surf_soil(:,1)=HgII_surf_soil(:,1)+(HgIIdry(:,1)*fsoil(:,1))
  leafpool_hg(:,1)=leafpool_hg(:,1)+(HgIIdry(:,1)*fstom(:,1))
  hleafpool_hg(:,1)=hleafpool_hg(:,1)+(HgIIdry(:,1)*fstom(:,1))
  
  photo_frac(:,1)=0.667577e+0_fp*(1.0e+0_fp- &
           exp(solrad1(:,mo)*(-1e+0_fp)*(0.01603e+0_fp)))
  !photo_frac equation is fit to curve in Rolfhus and Fitzgerald 
  
  WHERE (photo_frac(:,1) <= 0.0e+0_fp) 
     photo_frac(:,1)=0.0e+0_fp
  END WHERE

  photoreduced(:,1)=photo_frac(:,1)*(HgII_surf_soil(:,1))
  
  HgII_surf_soil(:,1)=HgII_surf_soil(:,1)-photoreduced(:,1)
  
  !photoreduce some wet deposition
  temp_hg(:,1)=HgIIwet(:,1)*photo_frac(:,1)
  photoreduced(:,1)=photoreduced(:,1)+temp_hg(:,1)
  HgIIwet(:,1)=HgIIwet(:,1)-temp_hg(:,1)

  WHERE (HgIIwet(:,1) < 0e+0_fp)
     HgIIwet(:,1)=0.0e+0_fp
  END WHERE

  !if there is rain - wash off HgII and add to HgII wet pool  
  WHERE (ppt1(:,mo) > 0e+0_fp .AND. airt1(:,mo) > 0e+0_fp)
     HgIIwet(:,1)=HgIIwet(:,1)+HgII_surf_soil(:,1)
     HgII_surf_soil(:,1)=0.0e+0_fp
  END WHERE

!$OMP END WORKSHARE
!$OMP END PARALLEL

END SUBROUTINE doHgDeposition
!EOC
