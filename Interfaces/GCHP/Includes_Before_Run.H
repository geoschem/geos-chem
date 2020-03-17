!EOC
!------------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 910.1 and      !
!          Harvard University Atmospheric Chemistry Modeling Group            !
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
  State_Met%FRCLND              = FRLAND                     ! 1
  State_Met%FRLAND              = FRLAND                     ! 1
  State_Met%FROCEAN             = FROCEAN                    ! 1
  State_Met%FRLAKE              = FRLAKE                     ! 1
  State_Met%FRLANDIC            = FRLANDIC                   ! 1
  State_Met%PHIS                = PHIS / 9.80665d0           ! m
  State_Met%GWETROOT            = GWETROOT                   ! 1
  State_Met%GWETTOP             = GWETTOP                    ! 1
  State_Met%LAI                 = LAI                        ! 1  
  State_Met%OMEGA               = OMEGA                      ! Pa s-1
  State_Met%PARDR               = PARDR                      ! W/m2
  State_Met%PARDF               = PARDF                      ! W/m2
  State_Met%PBLH                = PBLH                       ! m
  State_Met%TROPP               = TROPP                      ! hPa
  State_Met%PS1_WET             = PS1                        ! hPa
  State_Met%PS2_WET             = PS2                        ! hPa
  State_Met%SLP                 = SLP                        ! hPa
  State_Met%TMPU1               = TMPU1
  State_Met%TMPU2               = TMPU2
  State_Met%TS                  = TS                         ! K
  State_Met%TSKIN               = TSKIN                      ! K
  State_Met%SWGDN               = RADSWG                     ! W/m2
  State_Met%TO3                 = TO3                        ! dobson
  State_Met%SNODP               = SNODP                      ! m
  State_Met%SNOMAS              = SNOMAS                     ! kg/m2
  State_Met%SUNCOS              = zenith                     ! unitless
  State_Met%SUNCOSmid           = zenith                     ! 1
  State_Met%U                   = ULL                        ! m/s
  State_Met%U10M                = U10M                       ! m/s
  State_Met%USTAR               = USTAR                      ! m/s
  State_Met%V                   = VLL                        ! m/s
  State_Met%V10M                = V10M                       ! m/s
  State_Met%Z0                  = Z0                         ! m

  ! Set LWI based on land type fractions in imports: 0=water, 1=land, 2=ice
  DO J=1,State_Grid%NY
  DO I=1,State_Grid%NX
     IMAXLOC = MAXLOC( (/ FRLAND(I,J) + FRLANDIC(I,J) + FRLAKE(I,J),  &
                          FRSEAICE(I,J),                              &
                          FROCEAN(I,J) - FRSEAICE(I,J) /) )
     IF ( IMAXLOC(1) == 3 ) IMAXLOC(1) = 0 ! reset ocean to 0
     State_Met%LWI(I,J) = FLOAT( IMAXLOC(1) )
  ENDDO
  ENDDO

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

  ! NOTE: PLE is reversed in the vertical! Note that this will end up
  ! affecting DELP and therefore State_Met%PEDGE_DRY, but that the 
  ! major dry pressure variables are derived from PS1_DRY and PS2_DRY.
  ! Thes will, in turn, be determined from PS1_WET and PS2_WET, which 
  ! are set in the block above.
  ! Convert Pa -> hPa
  State_Met%PEDGE  (:,:,1:State_Grid%NZ+1) = PLE(:,:,State_Grid%NZ:0:-1) / 1d2

  ! These will be set based on PS1_WET and PS2_WET in gigc_chunk_mod.F90
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

  ! These fields are set to match the GCC conventions and to ensure that 
  ! transport is correctly implemented.
  State_Met%SPHU                =(SPHU1+SPHU2)*0.5d3         ! kg/kg -> g/kg
  State_Met%SPHU1               = SPHU1       *1.0d3         ! kg/kg -> g/kg
  State_Met%SPHU2               = SPHU2       *1.0d3         ! kg/kg -> g/kg
  State_Met%T                   =(TMPU1+TMPU2)*0.5d0         ! K

  ! Create total optical depth field                         ! 1
  State_Met%OPTD                = State_Met%TAUCLI + State_Met%TAUCLW

  !=========================================================================
  ! Offline lightning fields
  !=========================================================================
  State_Met%FLASH_DENS         = FLASH_DENS                  ! #/km2/s
  State_Met%CONV_DEPTH         = CONV_DEPTH                  ! m
