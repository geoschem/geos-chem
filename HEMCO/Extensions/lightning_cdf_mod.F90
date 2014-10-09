!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: lightning_cdf_mod.F90
!
! !DESCRIPTION: Module Lightning\_CDF\_Mod defines the cumulative probability
!  functions (CDF's).  These CDF's are used to partition the column sum
!  of lightning into the vertical. 
!\\
!\\
! References:
! \begin{itemize}
! \item Murray, L. T., Jacob, D. J., Logan, J. A., Hudman, R. C., and
!       Koshak, W. J.: \emph{Optimized regional and interannual variability 
!       of lightning in a global chemical transport model con- strained 
!       by LIS/OTD satellite data}, \underline{J. Geophys. Res.}, 
!       Atmospheres, 117, 2012.
! \item Ott, L. E., K. E. Pickering, G. L. Stenchikov, D. J. Allen,
!       A. J. DeCaria, B. Ridley, R.-F. Lin, S. Lang, and W.-K. Tao, 
!       \emph{Production of lightning NOx and its vertical distribution 
!       calculated  from three-dimensional cloud-scale chemical transport 
!       model simulations}, \underline{J. Geophys. Res.}, 115, D04301, 2010.
! \end{itemize}
!
! !INTERFACE:
!
MODULE Lightning_Cdf_Mod
!
! !DEFINED PARAMETERS:
!
  INTEGER, PUBLIC, PARAMETER :: NNLIGHT = 3200  ! # of data points in vertical
  INTEGER, PUBLIC, PARAMETER :: NLTYPE  = 4     ! # of types of lightning
!
! !PUBLIC DATA MEMBERS:
!
  REAL*8,  PUBLIC            :: PROFILE(NNLIGHT,NlTYPE)
!
! !REMARKS:
!  We now hardwire the PROFILE array, which is used in hcox_lightning_mod.F90.
!  This avoids having to read the data from an ASCII file, which cannot
!  be easily done in the ESMF environment.
!
! !REVISION HISTORY:
!  22 Jul 2014 - R. Yantosca - Initial version
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
! !IROUTINE: Init_Lightning_Cdf
!
! !DESCRIPTION: Subroutine Init\_Lightning\_Cdf initializes the PROFILE
!  array with the cumulative distribution function (CDF) look-up table data
!  from Ott et al, JGR, 2010.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Init_Lightning_Cdf()
!
! !REMARKS:
!  Cumulative emissions of NOx from lightning.  At the top of the cloud (16km) 
!  100% of NOx has been emitted into the column. At surface (0km) 0% of NOx 
!  has been emitted into the column. Taken from Table 2, Ott et al., JGR, 2010.
!
!  There are 4 CDF's for 4 different types of lightning:
!  (1) Lightning in tropical    marine      regions
!  (2) Lightning in tropical    continental regions
!  (3) Lightning in midlatitude continental regions
!  (4) Lightning in subtropical continental regions
! 
! !REVISION HISTORY: 
!  22 Jul 2014 - R. Yantosca - Initial version
!   8 Aug 2014 - R. Yantosca - Now get data statements from an include file
!EOP
!------------------------------------------------------------------------------
!BOC

    ! Now include a file with F90 assignment statements generated from
    ! the light_dist.ott2010.dat file.  This prevents us having to read an
    ! ASCII file in the ESMF environment.  If you want to update this data,
    ! all you have to do is to replace the include file. (bmy, 8/8/14)
#include "lightning_cdf_include.H"

  END SUBROUTINE Init_Lightning_Cdf
!EOC
END MODULE Lightning_Cdf_Mod
