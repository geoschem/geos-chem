!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: hcox_modis_lai_mod.F90 
!
! !DESCRIPTION: Module hcox\_modis\_lai\_mod.F90 contains routines and 
! variables to obtain the MODIS leaf area indeces (LAI, in cm2/cm2) as needed
! by some extensions (e.g. MEGAN). The extension object (ExtList) provides
! three LAI variables: interpolated LAI (GC\_LAI), current month LAI 
! (GC\_LAI\_CM), and previous month LAI (GC\_LAI\_PM). These variables refer
! to the MODIS calendar, which changes at mid-month, i.e. MODIS january starts
! only on Jan 15, etc.
!\\
!\\
! This module provides a way to define the three abovementioned variables
! from regular calendar monthly LAI data, i.e. data that changes on the first 
! day of a month.
!\\
!\\
! This module need not be called if the MODIS arrays are already calculated
! somewhere else. For example, GEOS-Chem calculates the current LAI values in
! module modis\_lai\_mod.F90 and these arrays are being used instead by HEMCO
! (through the HEMCO - GEOS-Chem interface hcoi\_gc\_main\_mod.F90).
!\\
!\\
! To use this module, it is assumed that the following MODIS LAI monthly data 
! (on a regular calendar) is specified in the HEMCO configuration file:
!\begin{itemize}
!\item MODIS\_LAI\_TM: MODIS leaf area index (cm2/cm2) for this calendar 
! month. Changes every first of the month.
!\item MODIS\_LAI\_TMM1: MODIS leaf area index (cm2/cm2) for the previous
! calendar month (this month minus one).
!\item MODIS\_LAI\_TMM2: MODIS leaf area index (cm2/cm2) for the second-
! previous calendar month (this month minus two). Only needed for variable 
! GC\_LAI\_PM.
! item MODIS\_LAI\_TMP1: MODIS leaf area index (cm2/cm2) for the upcoming
! calendar month (this month plus one). Only needed for variable GC\_LAI.
!\end{itemize} 
! !INTERFACE:
!
MODULE HCOX_MODIS_LAI_Mod
!
! !USES:
!
  USE HCO_Error_Mod
  USE HCOX_State_Mod, ONLY : Ext_State 
  USE HCO_State_Mod,  ONLY : HCO_State

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC  :: ExtState_Update_LAI
!
! !PRIVATE MEMBER FUNCTIONS:
!
! !REMARKS:
!
! !REVISION HISTORY:
!  08 Oct 2014 - C. Keller   - Initial version. 
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !PRIVATE MODULE VARIABLES:
!
    ! Interpolated leaf area index
    REAL(hp), POINTER   :: GC_LAI(:,:) => NULL() 

CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ExtState_Update_LAI
!
! !DESCRIPTION: SUBROUTINE ExtState\_Update\_LAI updates the leaf area index
! arrays in the passed ExtState object. It does so by using the leaf area index
! monthly values read through HEMCO.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE ExtState_Update_LAI( am_I_Root, HcoState, ExtState, RC )
!
! !USES:
!
    USE HCO_EMISLIST_MOD,     ONLY : HCO_GetPtr
    USE HCO_CLOCK_MOD,        ONLY : HcoClock_Get
    USE HCO_ARR_MOD,          ONLY : HCO_ArrInit
!
! !INPUT PARAMETERS:
!
    LOGICAL,          INTENT(IN)     :: am_I_Root
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(HCO_STATE),  POINTER        :: HcoState   ! HEMCO state object
    TYPE(EXT_STATE),  POINTER        :: ExtState   ! HEMCO extension object
    INTEGER,          INTENT(INOUT)  :: RC
!
! !REVISION HISTORY:
!  08 Oct 2014 - C. Keller   - Initial Version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! LOCAL VARIABLES:
!
    REAL(hp), POINTER   :: LAI_TM  (:,:) => NULL() ! This month
    REAL(hp), POINTER   :: LAI_TMM1(:,:) => NULL() ! This month - 1
    REAL(hp), POINTER   :: LAI_TMM2(:,:) => NULL() ! This month - 2
    REAL(hp), POINTER   :: LAI_TMP1(:,:) => NULL() ! This month + 1
    REAL(dp)            :: FRAC
    INTEGER             :: doy, cMidMon, dslmm, dbtwmm
    INTEGER             :: AS

    CHARACTER(LEN=255), PARAMETER :: &
       LOC = 'ExtState_Update_LAI (hcoi_gc_tools_mod.F90)'

    !=================================================================
    ! ExtState_Update_LAI begins here
    !=================================================================

    ! Assume success by default
    RC = HCO_SUCCESS

    ! Nothing to do if none of the LAI fields is enabled
    IF ( .NOT. ExtState%GC_LAI%DoUse     .AND. &
         .NOT. ExtState%GC_LAI_PM%DoUse  .AND. &
         .NOT. ExtState%GC_LAI_CM%DoUse         ) RETURN 

    ! Get some time variables: day of year (doy); mid-month day of
    ! current month (cMidMon); days since passing last mid-month date
    ! (dslmm); days between last-mid-month date and upcoming one (dbtwmm).
    CALL HcoClock_Get ( cDOY=doy, cMidMon=cMidMon, dslmm=dslmm, &
                        dbtwmm=dbtwmm, RC=RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Get pointers to leaf area index fields through HEMCO state. These
    ! fields start at the first of each month, i.e. the current month
    ! array (MODIS_LAI_CM) will hold LAI of January from January 1st 
    ! onwards. 
    ! The ExtState LAI variables change on mid-month days, however, and
    ! we need to assign those variables based on the current day of year.

    ! GC_LAI_CM is the LAI data for the current MODIS month. This is the 
    ! data of 'MODIS_LAI_TM' only if the current doy is past this month'
    ! mid-day. Otherwise, we need to point to the previous month' data.
    CALL HCO_GetPtr( am_I_Root, 'MODIS_LAI_TM'  , LAI_TM,   RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    CALL HCO_GetPtr( am_I_Root, 'MODIS_LAI_TMM1', LAI_TMM1, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Set GC_LAI_CM if needed. 
    IF ( ExtState%GC_LAI_CM%DoUse ) THEN

       ! If we're before mid month of this month, use previous month LAI
       IF ( doy < cMidMon ) THEN
          ExtState%GC_LAI_CM%Arr%Val => LAI_TMM1

       ! If we're past mid month of this month, use this month LAI
       ELSE
          ExtState%GC_LAI_CM%Arr%Val => LAI_TM
       ENDIF
    ENDIF

    ! Set GC_LAI_PM if needed.
    IF ( ExtState%GC_LAI_PM%DoUse ) THEN
    
       ! If we're before mid month of this month, use LAI data from two months
       ! before current month 
       IF ( doy < cMidMon ) THEN
          CALL HCO_GetPtr( am_I_Root, 'MODIS_LAI_TMM2', LAI_TMM2, RC )
          IF ( RC /= HCO_SUCCESS ) RETURN
          ExtState%GC_LAI_PM%Arr%Val => LAI_TMM2

       ! If we're past mid month of this month, use previous month LAI. 
       ELSE
          ExtState%GC_LAI_PM%Arr%Val => LAI_TMM1
       ENDIF
    ENDIF

    ! Set GC_LAI if needed. GC_LAI is the interpolated LAI between current (MODIS)
    ! month LAI and upcoming MODIS month LAI 
    IF ( ExtState%GC_LAI%DoUse ) THEN

       ! Fraction we're into the MODIS month (# of days since last mid-month day
       ! divided by total number of days between MODIS months.
       FRAC = DBLE( dslmm ) / DBLE( dbtwmm )

       ! Allocate array if necessary
       IF ( .NOT. ExtState%GC_LAI%Arr%Alloc ) THEN
          CALL HCO_ArrInit ( ExtState%GC_LAI%Arr, HcoState%NX, HcoState%NY, RC )
          IF ( RC /= HCO_SUCCESS ) RETURN
       ENDIF

       ! If we're before mid month of this month, interpolate from previous and 
       ! current calendar month LAI data
       IF ( doy < cMidMon ) THEN

          ExtState%GC_LAI%Arr%Val = LAI_TMM1 + ( (LAI_TM - LAI_TMM1) * FRAC )  

       ! If we're past mid month of this month, interpolate from current and next
       ! calendar month LAI data. 
       ELSE
          CALL HCO_GetPtr( am_I_Root, 'MODIS_LAI_TMP1', LAI_TMP1, RC )
          IF ( RC /= HCO_SUCCESS ) RETURN
          ExtState%GC_LAI%Arr%Val = LAI_TM + ( (LAI_TMP1 - LAI_TM) * FRAC )  

       ENDIF

    ENDIF

  END SUBROUTINE ExtState_Update_LAI
!EOC
END MODULE HCOX_MODIS_LAI_MOD
