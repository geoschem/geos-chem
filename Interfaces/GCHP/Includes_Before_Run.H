!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !INCLUDE: Includes_Before_Run.H
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
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC

  !=========================================================================
  ! Two-dimensional fields
  !=========================================================================

  State_Met%ALBD                = ALBD                       ! 1
  State_Met%CLDFRC              = CLDFRC                     ! 1
  State_Met%EFLUX               = EFLUX                      ! W m-2
  State_Met%HFLUX               = HFLUX                      ! W/m2
  State_Met%FRLAND              = FRLAND                     ! 1
  State_Met%FROCEAN             = FROCEAN                    ! 1
  State_Met%FRLAKE              = FRLAKE                     ! 1
  State_Met%FRLANDIC            = FRLANDIC                   ! 1
  State_Met%FRSEAICE            = FRSEAICE                   ! 1
  State_Met%FRSNO               = FRSNO                      ! 1
  State_Met%QV2M                = QV2M                       ! 1
  State_Met%PHIS                = PHIS / 9.80665d0           ! m
  State_Met%GWETROOT            = GWETROOT                   ! 1
  State_Met%GWETTOP             = GWETTOP                    ! 1
  State_Met%LAI                 = LAI                        ! 1
  State_Met%PARDR               = PARDR                      ! W/m2
  State_Met%PARDF               = PARDF                      ! W/m2
  State_Met%PBLH                = PBLH                       ! m
  State_Met%TROPP               = TROPP                      ! hPa
  State_Met%PS1_WET             = PS1                        ! hPa
  State_Met%PS2_WET             = PS2                        ! hPa
  State_Met%SLP                 = SLP                        ! hPa
  State_Met%TS                  = TS                         ! K
  State_Met%TSKIN               = TSKIN                      ! K
  State_Met%SWGDN               = RADSWG                     ! W/m2
  State_Met%TO3                 = TO3                        ! dobson
  State_Met%SNODP               = SNODP                      ! m
  State_Met%SNOMAS              = SNOMAS                     ! kg/m2
  State_Met%SUNCOS              = zenith                     ! unitless
  State_Met%SUNCOSmid           = zenith                     ! 1
  State_Met%U10M                = U10M                       ! m/s
  State_Met%USTAR               = USTAR                      ! m/s
  State_Met%V10M                = V10M                       ! m/s
  State_Met%Z0                  = Z0                         ! m

  ! Get UV albedo for photolysis
  Ptr2d => NULL()
  CALL MAPL_GetPointer ( IMPORT, Ptr2D, 'UV_ALBEDO', &
                         notFoundOK=.TRUE., __RC__ )
  If ( ASSOCIATED(Ptr2D) ) Then
     State_Met%UVALBEDO(:,:) = Ptr2D(:,:)
     Ptr2D => NULL()
  ENDIF

  !=========================================================================
  ! Three-dimensional fields on level edges
  !=========================================================================

  ! kg/m2/s
  State_Met%CMFMC  (:,:,1:State_Grid%NZ+1) = CMFMC  (:,:,0:State_Grid%NZ)
  State_Met%PFICU  (:,:,1:State_Grid%NZ+1) = PFICU  (:,:,0:State_Grid%NZ)
  State_Met%PFILSAN(:,:,1:State_Grid%NZ+1) = PFILSAN(:,:,0:State_Grid%NZ)
  State_Met%PFLCU  (:,:,1:State_Grid%NZ+1) = PFLCU  (:,:,0:State_Grid%NZ)
  State_Met%PFLLSAN(:,:,1:State_Grid%NZ+1) = PFLLSAN(:,:,0:State_Grid%NZ)

  if (meteorology_vertical_index_is_top_down) then
    z_lb = lbound(State_Met%CMFMC,3)
    z_ub = ubound(State_Met%CMFMC,3)
    State_Met%CMFMC(:,:,:)     = State_Met%CMFMC(:,:,z_ub:z_lb:-1)
    State_Met%PFICU(:,:,:)     = State_Met%PFICU(:,:,z_ub:z_lb:-1)
    State_Met%PFILSAN(:,:,:)   = State_Met%PFILSAN(:,:,z_ub:z_lb:-1)
    State_Met%PFLCU(:,:,:)     = State_Met%PFLCU(:,:,z_ub:z_lb:-1)
    State_Met%PFLLSAN(:,:,:)   = State_Met%PFLLSAN(:,:,z_ub:z_lb:-1)
  end if

  ! NOTE: PLE is reversed in the vertical! Note that this will end up
  ! affecting DELP and therefore State_Met%PEDGE_DRY, but that the 
  ! major dry pressure variables are derived from PS1_DRY and PS2_DRY.
  ! Thes will, in turn, be determined from PS1_WET and PS2_WET, which 
  ! are set in the block above.
  ! Convert Pa -> hPa
  State_Met%PEDGE  (:,:,1:State_Grid%NZ+1) = PLE(:,:,State_Grid%NZ:0:-1) / 1d2

  ! These will be set based on PS1_WET and PS2_WET in gchp_chunk_mod.F90
  State_Met%PS1_DRY(:,:) = 0.0d0 
  State_Met%PS2_DRY(:,:) = 0.0d0

  ! Calculate CLDTOPS (highest location of CMFMC in the column)
  Do J=1,State_Grid%NY
  Do I=1,State_Grid%NX
     State_Met%CldTops(I,J) = 1
     Do L = State_Grid%NZ, 1, -1
        If ( State_Met%CMFMC(I,J,L) > 0d0 ) Then
           State_Met%CldTops(I,J) = L + 1
           Exit
        End If
     End Do
  End Do
  End Do


  ! Read MODIS leaf area index (LAI) from imports of post-processed MODIS files.
  ! The third dimension is land type and not level, possible to do with MAPL
  ! only because the # of land types is the same as # of level edges.
  State_Met%XLAI_NATIVE(:,:,:) = XLAIMULTI

  !=========================================================================
  ! Three-dimensional fields on level centers
  !=========================================================================
  State_Met%CLDF                = CLDF                       ! 1
  State_Met%DTRAIN              = DTRAIN                     ! kg/m2/s
  State_Met%DQRCU               = DQRCU                      ! kg/kg/s
  State_Met%DQRLSAN             = DQRLSAN                    ! kg/kg/s
  State_Met%QI                  = QI                         ! kg/kg
  State_Met%QL                  = QL                         ! kg/kg
  State_Met%RH                  = RH      * 1d2              ! 1 -> %
  State_Met%TAUCLI              = TAUCLI                     ! 1
  State_Met%TAUCLW              = TAUCLW                     ! 1
  State_Met%REEVAPCN            = REEVAPCN                   ! kg/kg/s
  State_Met%REEVAPLS            = REEVAPLS                   ! kg/kg/s
  State_Met%TMPU1               = TMPU1                      ! K
  State_Met%TMPU2               = TMPU2                      ! K
  State_Met%OMEGA               = OMEGA                      ! Pa s-1
  State_Met%U                   = UA                         ! m/s
  State_Met%V                   = VA                         ! m/s
  ! These fields are set to match the GCC conventions and to ensure that 
  ! transport is correctly implemented.
  State_Met%SPHU1               = SPHU1       *1.0d3         ! kg/kg -> g/kg
  State_Met%SPHU2               = SPHU2       *1.0d3         ! kg/kg -> g/kg

    ! If meteorology vertical index is top down, flip imports coming from ExtData.
  if (meteorology_vertical_index_is_top_down) then
    z_lb = lbound(State_Met%CLDF,3)
    z_ub = ubound(State_Met%CLDF,3)

    State_Met%CLDF(:,:,:)     = State_Met%CLDF(:,:,z_ub:z_lb:-1)
    State_Met%DTRAIN(:,:,:)   = State_Met%DTRAIN(:,:,z_ub:z_lb:-1)
    State_Met%DQRCU(:,:,:)    = State_Met%DQRCU(:,:,z_ub:z_lb:-1)
    State_Met%DQRLSAN(:,:,:)  = State_Met%DQRLSAN(:,:,z_ub:z_lb:-1)
    State_Met%QI(:,:,:)       = State_Met%QI(:,:,z_ub:z_lb:-1)
    State_Met%QL(:,:,:)       = State_Met%QL(:,:,z_ub:z_lb:-1)
    State_Met%RH(:,:,:)       = State_Met%RH(:,:,z_ub:z_lb:-1)
    State_Met%TAUCLI(:,:,:)   = State_Met%TAUCLI(:,:,z_ub:z_lb:-1)
    State_Met%TAUCLW(:,:,:)   = State_Met%TAUCLW(:,:,z_ub:z_lb:-1)
    State_Met%REEVAPCN(:,:,:) = State_Met%REEVAPCN(:,:,z_ub:z_lb:-1)    
    State_Met%REEVAPLS(:,:,:) = State_Met%REEVAPLS(:,:,z_ub:z_lb:-1)
    State_Met%SPHU1(:,:,:)     = State_Met%SPHU1(:,:,z_ub:z_lb:-1)
    State_Met%SPHU2(:,:,:)     = State_Met%SPHU2(:,:,z_ub:z_lb:-1)
    State_Met%TMPU1(:,:,:)     = State_Met%TMPU1(:,:,z_ub:z_lb:-1)
    State_Met%TMPU2(:,:,:)     = State_Met%TMPU2(:,:,z_ub:z_lb:-1)
    State_Met%U(:,:,:)         = State_Met%U(:,:,z_ub:z_lb:-1)
    State_Met%V(:,:,:)         = State_Met%V(:,:,z_ub:z_lb:-1)
    State_Met%OMEGA(:,:,:)     = State_Met%OMEGA(:,:,z_ub:z_lb:-1)
  end if

  State_Met%SPHU                =(State_Met%SPHU1+State_Met%SPHU2)*0.5d0
  State_Met%T                   =(State_Met%TMPU1+State_Met%TMPU2)*0.5d0

  ! Create total optical depth field                         ! 1
  State_Met%OPTD                = State_Met%TAUCLI + State_Met%TAUCLW

  !=========================================================================
  ! Offline lightning fields
  !=========================================================================
  State_Met%FLASH_DENS         = FLASH_DENS                  ! #/km2/s
  State_Met%CONV_DEPTH         = CONV_DEPTH                  ! m

#ifdef ADJOINT
  IF (ASSOCIATED(CostFuncMask)) THEN
      if (MAPL_Am_I_Root()) &
          WRITE(*,*) ' Loading adjoint cost function mask'
      ! cost function mask
      State_Chm%CostFuncMask = CostFuncMask
  ELSEIF (Input_Opt%CF_IMIN > 0 .and. Input_Opt%CF_JMIN > 0) THEN
      IF (FIRST) THEN
         if (MAPL_Am_I_Root()) &
             WRITE(*,*) ' Cost function range supplied.'
       WRITE(*,1027) Input_Opt%thisCPU,   &
            Input_Opt%CF_IMIN, Input_Opt%CF_IMAX, &
            Input_Opt%CF_JMIN, Input_Opt%CF_JMAX, &
            Input_Opt%CF_LMIN, Input_Opt%CF_LMAX
1027   FORMAT('CF on Pet ', i3, ' I = (', i3, ', ', i3, ') &
             J = ( ', i3, ', ', i3, ') &
             L = (', i3, ', ', i3, ')')


         State_Chm%CostFuncMask = 0d0
         DO L=Input_Opt%CF_LMIN,Input_Opt%CF_LMAX
         DO J=Input_Opt%CF_JMIN,Input_opt%CF_JMAX
            DO I=Input_opt%CF_IMIN, Input_Opt%CF_IMAX
               State_Chm%CostFuncMask(I,J,L) = 1d0
            ENDDO
         ENDDO
         ENDDO
      ENDIF
  ENDIF

#endif
