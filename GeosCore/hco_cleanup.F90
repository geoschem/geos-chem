!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: hco_cleanup
!
! !DESCRIPTION: Routine HCO\_Cleanup deallocates all HEMCO arrays.
! Don't pack this into a module so that HCO\_Cleanup can also be
! called by HCO\_STOP! 
!\\
! !INTERFACE: 
!
      SUBROUTINE HCO_CLEANUP 
!
! !USES:
!
      USE HCO_TYPE_MOD,      ONLY : Cleanup_HCO_State
      USE HCO_TYPE_MOD,      ONLY : HCO_State
      USE HCO_EMISLL_MOD,    ONLY : Cleanup_EmisLL
      USE HCO_EMISLL_MOD,    ONLY : EmisCont 
      USE HCO_TIME_MOD,      ONLY : Cleanup_tSlc
      USE HCO_TIME_MOD,      ONLY : Cleanup_Clock
      USE HCO_TIME_MOD,      ONLY : HcoClock
      USE HCO_READLIST_MOD,  ONLY : ReadList_Cleanup
      USE HCO_EXTS_MOD,      ONLY : ExtFinal
      USE HCO_MAIN_MOD,      ONLY : Clock, EmisLL, HcoStInt
!
! !REVISION HISTORY: 
!  27 May 2012 - C. Keller    - Initialization
!  22 Aug 2013 - C. Keller    - Stripped from ng_emissions_mod.F
!  23 Sep 2013 - C. Keller    - Now a stand-alone routine
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!

      !=================================================================
      ! HCO_CLEANUP 
      !=================================================================

      ! Cleanup the emissions linked list. Do this before cleaning up
      ! ReadList since some pointers in EmisLL point to variables stored
      ! in the ReadList!
      CALL Cleanup_EmisLL ( EmisLL )

      ! Cleanup the time slice pointers
      CALL Cleanup_tSlc

      ! Remove the clock
      CALL Cleanup_Clock ( Clock )

      ! Remove the ReadList
      CALL ReadList_Cleanup

      ! Remove HEMCO state 
      CALL Cleanup_HCO_State ( HcoStInt ) 

      ! Remove extensions entries
      CALL ExtFinal

      END SUBROUTINE HCO_CLEANUP 
!EOC
