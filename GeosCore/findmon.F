!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: findmon
!
! !DESCRIPTION: Function FINDMON finds which month JDAY (day of this year) 
!  is in.  FINDMON is called by the Leaf Area Index routine rdlai.f.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE FINDMON( JDAY, INMONTH, INYEAR, MM, YYYY, STARTDAY )
!
! !USES:
!
      IMPLICIT NONE
#     include "define.h"
!
! !INPUT PARAMETERS: 
!
      INTEGER, INTENT(IN)  :: JDAY           ! Current day of year
      INTEGER, INTENT(IN)  :: INMONTH        ! Current month
      INTEGER, INTENT(IN)  :: INYEAR         ! Current year
      INTEGER, INTENT(IN)  :: STARTDAY(13)   ! Starting days for LAI data
!
! !OUTPUT PARAMETERS:
!
      INTEGER, INTENT(OUT) :: MM             ! Output month for LAI
      INTEGER, INTENT(OUT) :: YYYY           ! Output year for LAI
!
! !REVISION HISTORY: 
!  05 Jan 1994 - Y. H. Wang, G.M. Gardner, D. Jacob - Initial version
!  (1 ) Updated comments, cosmetic changes (bmy, 4/4/03)
!  (2 ) Add the current simulation year as input & the current LAI as output.
!       This is necessary for reading in MODIS LAI (mpb,2009).
!  08 Dec 2009 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      !=================================================================
      ! FINDMON begins here!
      !=================================================================
      IF ( JDAY < STARTDAY(1) ) THEN
         MM = 12
         YYYY = INYEAR - 1      ! (mpb,2008)
      ELSE IF ( JDAY < STARTDAY(INMONTH) ) THEN
         MM = INMONTH-1
         YYYY = INYEAR           ! (mpb,2008)
      ELSE
         MM = INMONTH
         YYYY = INYEAR           ! (mpb,2008)
      ENDIF

      END SUBROUTINE FINDMON
!EOC
