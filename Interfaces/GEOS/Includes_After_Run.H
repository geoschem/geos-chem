!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !INCLUDE: Includes_After_Run.H
!
! !DESCRIPTION: This include file contains the array assignments that need
!  to be made BEFORE the call to the Run method of the 
!  GEOSCHEMchem\_GridCompMod.F90 code.  These array assignments take data out 
!  of the import and internal states and saves them into derived-type objects
!  that are used to pass information between GEOS-Chem and the ESMF interface.
!\\
!\\
!  These assignments were placed into this separate include file to avoid 
!  bogging down the GEOSCHEMchem\_GridCompMod.F90 module.
!\\
!\\
! !REVISION HISTORY: 
!  10 Oct 2012 - M. Long     - Initial version, for grid-independent model
!  10 Oct 2012 - R. Yantosca - Added ProTeX headers, updated comments
!  10 Oct 2012 - R. Yantosca - Converted to F90 free-format
!  10 Oct 2012 - R. Yantosca - Added Bromine species for v9-01-03
!  16 Oct 2012 - R. Yantosca - Changed GC_MET to State_Met
!  16 Oct 2012 - R. Yantosca - Changed GC_STATE to State_Chm
!  16 Oct 2012 - R. Yantosca - Changed GC_STATE%CSPEC to State_Chm%Species
!  16 Oct 2012 - R. Yantosca - Renamed to Includes_After_Run.H (freeform)
!  07 Mar 2013 - R. Yantosca - Now add error checks to prevent OOB errors
!  27 Oct 2014 - C. Keller   - Removed species copying into internal state.
!                              This is only required upon finalization.
!  08 May 2015 - C. Keller   - Now fill species on every (chemistry) time 
!                              step. This is critical for checkpointing.
!  11 Dec 2018 - E. Lundgren - Change names now that stored in State_Chm
!                              (old -> new)
!                              ARCHV_H2O2s/SO2s -> H2O2/SO2AfterChem, 
!                              ARCHV_DRY/WET_TOTN -> DryDep/WetNitrogent, 
!                              HSAVE_FOR_KPP -> KPPHvalue
! See git history for additional changes.
!EOP
!------------------------------------------------------------------------------
!BOC

  !=========================================================================
  ! Pass tracers from GEOS-Chem tracers array to internal state
  !=========================================================================
  DO I = 1, SIZE(Int2Spc,1)
     IF ( Int2Spc(I)%ID <= 0 ) CYCLE
     IF ( .NOT. ASSOCIATED(Int2Spc(I)%Internal) ) THEN
        CYCLE
     ELSE
       Int2Spc(I)%Internal(:,:,:) = &
                         State_Chm%Species(Int2Spc(I)%ID)%Conc(:,:,LM:1:-1) 
     ENDIF
     ! testing only
  !if(MAPL_am_I_Root()) write(*,*) 'Copied to internal: ',I,Int2Spc(I)%ID,trim(Int2Spc(I)%Name),MINVAL(State_Chm%Species(Int2Spc(I)%ID)%Conc(:,:,LM:1:-1)),MAXVAL(State_Chm%Species(Int2Spc(I)%ID)%Conc(:,:,LM:1:-1)),SUM(State_Chm%Species(Int2Spc(I)%ID)%Conc(:,:,LM:1:-1))/IM/JM/LM
     State_Chm%Species(I)%Units = KG_SPECIES_PER_KG_TOTAL_AIR
  ENDDO


  !=========================================================================
  ! Various other archived variables needed in internal state. 
  !=========================================================================
  CALL MAPL_GetPointer( INTSTATE, Ptr3D, 'H2O2AfterChem', notFoundOK=.TRUE., __RC__ ) 
  IF ( ASSOCIATED(Ptr3D) .AND. ASSOCIATED(State_Chm%H2O2AfterChem) ) THEN
     Ptr3D(:,:,LM:1:-1) = State_Chm%H2O2AfterChem
  ENDIF
  Ptr3D => NULL()

  CALL MAPL_GetPointer( INTSTATE, Ptr3D, 'SO2AfterChem', notFoundOK=.TRUE., __RC__ ) 
  IF ( ASSOCIATED(Ptr3D) .AND. ASSOCIATED(State_Chm%SO2AfterChem) ) THEN
     Ptr3D(:,:,LM:1:-1) = State_Chm%SO2AfterChem
  ENDIF
  Ptr3D => NULL()

  CALL MAPL_GetPointer( INTSTATE, Ptr2D, 'DryDepNitrogen', notFoundOK=.TRUE., __RC__ ) 
  IF ( ASSOCIATED(Ptr2D) .AND. ASSOCIATED(State_Chm%DryDepNitrogen) ) THEN
     Ptr2D = State_Chm%DryDepNitrogen
  ENDIF
  Ptr2D => NULL()

  CALL MAPL_GetPointer( INTSTATE, Ptr2D, 'WetDepNitrogen', notFoundOK=.TRUE., __RC__ ) 
  IF ( ASSOCIATED(Ptr2D) .AND. ASSOCIATED(State_Chm%WetDepNitrogen) ) THEN
     Ptr2D = State_Chm%WetDepNitrogen
  ENDIF
  Ptr2D => NULL()

  CALL MAPL_GetPointer( INTSTATE, Ptr3D, 'KPPHvalue' , notFoundOK=.TRUE., __RC__ ) 
  IF ( ASSOCIATED(Ptr3D) .AND. ASSOCIATED(State_Chm%KPPHvalue) ) THEN
     Ptr3D(:,:,:) = 0.0
     Ptr3D(:,:,LM:1:-1) = REAL(State_Chm%KPPHvalue(:,:,1:LM),KIND=ESMF_KIND_R4)
  ENDIF
  Ptr3D => NULL()

  ! Sesquiterpene concentrations. Archive internally in kgC/m3 (convert from kgC/box).
  CALL MAPL_GetPointer( INTSTATE, Ptr3D, 'SESQ_CARBDENS', notFoundOK=.TRUE., __RC__ )
  IF ( ASSOCIATED(Ptr3D) .AND. ASSOCIATED(State_Chm%ORVCsesq) ) THEN
     Ptr3D(:,:,LM:1:-1) = State_Chm%ORVCsesq(:,:,1:LM) / State_Met%AIRVOL(:,:,1:LM)
  ENDIF
  Ptr3D => NULL()

  ! Saved field from ISORROPIA needed for computing optical properties
  CALL MAPL_GetPointer( INTSTATE, Ptr3D, 'AeroH2O_SNA', notFoundOK=.TRUE., __RC__ )
  IF ( ASSOCIATED(Ptr3D) .AND. ASSOCIATED(State_Chm%AeroH2O) ) THEN
  Ptr3D(:,:,LM:1:-1) = State_Chm%AeroH2O(:,:,:,NDUST+1)
  ENDIF
  Ptr3D => NULL()

  CALL MAPL_GetPointer( INTSTATE, Ptr3D, 'STATE_PSC', notFoundOK=.TRUE., __RC__ )
  IF ( ASSOCIATED(Ptr3D) .AND. ASSOCIATED(State_Chm%STATE_PSC) ) THEN
     Ptr3d(:,:,LM:1:-1) = State_Chm%STATE_PSC(:,:,:)
  ENDIF
  Ptr3D => NULL()

!EOF
