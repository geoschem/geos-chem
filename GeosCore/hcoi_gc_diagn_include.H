!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: hcoi_gc_diagn_include.H
!
! !DESCRIPTION: Contains emissions category definitions as used within the
!  HEMCO_Config.rc file.
!\\
!\\
! !REMARKS:
!  These categories hould correspond with the HEMCO_Config.rc file.  If you
!  change the category listings in HEMCO_Config.rc, then update this file
!  accordingly.
!
! !REVISION HISTORY:
!  01 Oct 1995 - R. Yantosca - Initial version
!  08 Dec 2009 - R. Yantosca - Added ProTeX headers
!  22 Feb 2015 - C. Keller   - Added eruptive and degassing volcanoes.
!  27 Jul 2018 - M. Sulprizio- Added trash burning category (E. Marais)
!
! !DEFINED PARAMETERS:
!
  INTEGER, PARAMETER :: CATEGORY_ANTHRO        = 1
  INTEGER, PARAMETER :: CATEGORY_NATURAL       = 3
  INTEGER, PARAMETER :: CATEGORY_BIOGENIC      = 4
  INTEGER, PARAMETER :: CATEGORY_BIOMASS       = 5
  INTEGER, PARAMETER :: CATEGORY_SHIP          = 10
  INTEGER, PARAMETER :: CATEGORY_TRASH         = 12
  INTEGER, PARAMETER :: CATEGORY_AIRCRAFT      = 20
  INTEGER, PARAMETER :: CATEGORY_VOLCANO       = 50
  INTEGER, PARAMETER :: CATEGORY_VOLCANO_ERUPT = 51
  INTEGER, PARAMETER :: CATEGORY_VOLCANO_DEGAS = 52
!EOP
!------------------------------------------------------------------------------
!BOC
