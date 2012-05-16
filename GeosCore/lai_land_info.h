! $Id: lai_land_info.h,v 1.1 2009/08/25 20:45:21 bmy Exp $
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !INCLUDE: lai_land_info.h
!
! !DESCRIPTION: This include file contains the various dimensions for
!  the leaf-area-index and Olson land type data.
!\\
!\\
! !DEFINED PARAMETERS: 
!
      ! Number of Olson land types
      INTEGER, PARAMETER :: N_OLSON_TYPES = 74    

      ! Max # of Olson land types that can fit into any one grid box
      INTEGER, PARAMETER :: N_OLSON_LOCAL = 15
!
! !REVISION HISTORY: 
! 24 Mar 2009 - R. Yantosca - Initial version
!EOP
!------------------------------------------------------------------------------
