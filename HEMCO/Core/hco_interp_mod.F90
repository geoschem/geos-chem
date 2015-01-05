!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: hco_interp_mod.F90 
!
! !DESCRIPTION: Module HCO\_INTERP\_MOD contains routines to
! interpolate input data onto the HEMCO grid. 
!\\
!\\
MODULE HCO_INTERP_MOD
!
! !USES:
!
  USE HCO_Error_Mod
  USE HCO_State_Mod,       ONLY : Hco_State
  USE HCO_DataCont_Mod,    ONLY : ListCont

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC  :: ModelLev_Interpolate 
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PRIVATE :: COLLAPSE
  PRIVATE :: INFLATE 
!
! !MODULE PARAMETER
!
  ! AP parameter of native GEOS-5 grid. Needed to remap GEOS-5 data from native
  ! onto the reduced vertical grid.
  REAL(hp), TARGET :: G5_EDGE_NATIVE(73) = (/                          & 
              0.000000e+00_hp, 4.804826e-02_hp, 6.593752e+00_hp, 1.313480e+01_hp, &
              1.961311e+01_hp, 2.609201e+01_hp, 3.257081e+01_hp, 3.898201e+01_hp, &
              4.533901e+01_hp, 5.169611e+01_hp, 5.805321e+01_hp, 6.436264e+01_hp, &
              7.062198e+01_hp, 7.883422e+01_hp, 8.909992e+01_hp, 9.936521e+01_hp, &
              1.091817e+02_hp, 1.189586e+02_hp, 1.286959e+02_hp, 1.429100e+02_hp, &
              1.562600e+02_hp, 1.696090e+02_hp, 1.816190e+02_hp, 1.930970e+02_hp, &
              2.032590e+02_hp, 2.121500e+02_hp, 2.187760e+02_hp, 2.238980e+02_hp, &
              2.243630e+02_hp, 2.168650e+02_hp, 2.011920e+02_hp, 1.769300e+02_hp, &
              1.503930e+02_hp, 1.278370e+02_hp, 1.086630e+02_hp, 9.236572e+01_hp, & 
              7.851231e+01_hp, 6.660341e+01_hp, 5.638791e+01_hp, 4.764391e+01_hp, &
              4.017541e+01_hp, 3.381001e+01_hp, 2.836781e+01_hp, 2.373041e+01_hp, &
              1.979160e+01_hp, 1.645710e+01_hp, 1.364340e+01_hp, 1.127690e+01_hp, &
              9.292942e+00_hp, 7.619842e+00_hp, 6.216801e+00_hp, 5.046801e+00_hp, &
              4.076571e+00_hp, 3.276431e+00_hp, 2.620211e+00_hp, 2.084970e+00_hp, &
              1.650790e+00_hp, 1.300510e+00_hp, 1.019440e+00_hp, 7.951341e-01_hp, &
              6.167791e-01_hp, 4.758061e-01_hp, 3.650411e-01_hp, 2.785261e-01_hp, &
              2.113490e-01_hp, 1.594950e-01_hp, 1.197030e-01_hp, 8.934502e-02_hp, &
              6.600001e-02_hp, 4.758501e-02_hp, 3.270000e-02_hp, 2.000000e-02_hp, &
              1.000000e-02_hp /)

  ! AP parameter of native GEOS-4 grid. Needed to remap GEOS-4 data from native
  ! onto the reduced vertical grid.
  REAL(hp), TARGET :: G4_EDGE_NATIVE(56) = (/       &
                    0.000000_hp,   0.000000_hp,  12.704939_hp, &
                   35.465965_hp,  66.098427_hp, 101.671654_hp, & 
                  138.744400_hp, 173.403183_hp, 198.737839_hp, &
                  215.417526_hp, 223.884689_hp, 224.362869_hp, &
                  216.864929_hp, 201.192093_hp, 176.929993_hp, &
                  150.393005_hp, 127.837006_hp, 108.663429_hp, & 
                   92.365662_hp,  78.512299_hp,  66.603378_hp, & 
                   56.387939_hp,  47.643932_hp,  40.175419_hp, &
                   33.809956_hp,  28.367815_hp,  23.730362_hp, & 
                   19.791553_hp,  16.457071_hp,  13.643393_hp, & 
                   11.276889_hp,   9.292943_hp,   7.619839_hp, &  
                    6.216800_hp,   5.046805_hp,   4.076567_hp, &
                    3.276433_hp,   2.620212_hp,   2.084972_hp, &  
                    1.650792_hp,   1.300508_hp,   1.019442_hp, &  
                    0.795134_hp,   0.616779_hp,   0.475806_hp, &  
                    0.365041_hp,   0.278526_hp,   0.211349_hp, & 
                    0.159495_hp,   0.119703_hp,   0.089345_hp, &  
                    0.066000_hp,   0.047585_hp,   0.032700_hp, &  
                    0.020000_hp,   0.010000_hp /)
!
! !REVISION HISTORY:
!  30 Dec 2014 - C. Keller   - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ModelLev_Interpolate 
!
! !DESCRIPTION: Subroutine ModelLev\Interpolate puts 3D data from an 
! arbitrary number of model levels onto the vertical levels of the simulation 
! grid. Since the input data is already on model levels, this is only to 
! inflate/collapse fields between native/reduced vertical levels, e.g. from
! 72 native GEOS-5 levels onto the reduced 47 levels.
!\\
!\\
! The input data (REGR\_4D) is expected to be already horizontally regridded. 
! The 4th dimension of REGR\_4D denotes time. 
!\\
!\\
! The 3rd dimension of REGR\_3D holds the vertical levels. It is assumed that
! these are model levels, starting at the surface (level 1). If the input
! data holds 72 input levels, this is interpreted as native data and will
! be collapsed onto the reduced grid. If the input data holds X <=47 levels, 
! these levels are interpreted as levels 1-X of the reduced grid. In other
! words, input data with 33 levels will be interpreted as 33 levels on the
! reduced grid, and the data is accordingly mapped onto the simulation grid.
! If data becomes inflated or collapsed, the output data will always extent
! over all vertical levels of the simulation grid. If necessary, the unused
! upper levels will be filled with zeros. If no data interpolation is needed, 
! the vertical extent of the output data is limited to the number of used 
! levels. For instance, if the input data has 5 vertical levels, the output
! array will only extent over those 5 (bottom) levels.
!\\
!\\
! Currently, this routine can remap the following combinations:
!\begin{itemize}
! \item Native GEOS-5 onto reduced GEOS-5 (72 --> 47 levels)
! \item Reduced GEOS-5 onto native GEOS-5 (47 --> 72 levels)
! \item Native GEOS-4 onto reduced GEOS-4 (56 --> 30 levels)
! \item Reduced GEOS-4 onto native GEOS-4 (30 --> 56 levels)
!\end{itemize}
! Interpolation from GEOS-5 onto GEOS-4 levels is currently not supported.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE ModelLev_Interpolate ( am_I_Root, HcoState, REGR_4D, Lct, RC )
!
! !USES:
!
    USE HCO_FileData_Mod,   ONLY : FileData_ArrCheck
!
! !INPUT PARAMETERS:
!
    LOGICAL,          INTENT(IN   )  :: am_I_Root         ! Are we on the root CPU?
    TYPE(HCO_State),  POINTER        :: HcoState          ! HEMCO state object
    REAL(hp),         POINTER        :: REGR_4D(:,:,:,:)  ! 4D input data
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ListCont),   POINTER        :: Lct               ! HEMCO list container
    INTEGER,          INTENT(INOUT)  :: RC                ! Success or failure?
!
! !REVISION HISTORY:
!  30 Dec 2014 - C. Keller   - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
! 
! !LOCAL VARIABLES:
!
    INTEGER                 :: nx, ny, nz, nt
    INTEGER                 :: minlev, nlev, nout
    INTEGER                 :: L, T, NL
    LOGICAL                 :: verb, infl, clps
    CHARACTER(LEN=255)      :: MSG

    !=================================================================
    ! ModelLev_Interpolate begins here
    !=================================================================

    ! Enter
    CALL HCO_ENTER ('ModelLev_Interpolate (hcoio_interpolate_mod.F90)' , RC )
    IF ( RC /= HCO_SUCCESS ) RETURN
    
    ! Check for verbose mode
    verb = HCO_VERBOSE_CHECK() .AND. am_I_Root
    IF ( verb ) THEN
       MSG = 'Vertically interpolate model levels: '//TRIM(Lct%Dct%cName)
       CALL HCO_MSG(MSG)
    ENDIF

    ! Get HEMCO grid dimensions 
    nx = HcoState%NX
    ny = HcoState%NY
    nz = HcoState%NZ 

    ! Input data must be on horizontal HEMCO grid
    IF ( SIZE(REGR_4D,1) /= nx ) THEN
       WRITE(MSG,*) 'x dimension mismatch ', TRIM(Lct%Dct%cName), &
          ': ', nx, SIZE(REGR_4D,1)
       CALL HCO_ERROR( MSG, RC )
       RETURN
    ENDIF 
    IF ( SIZE(REGR_4D,2) /= ny ) THEN
       WRITE(MSG,*) 'y dimension mismatch ', TRIM(Lct%Dct%cName), &
          ': ', ny, SIZE(REGR_4D,2)
       CALL HCO_ERROR( MSG, RC )
       RETURN
    ENDIF 

    ! Get vertical and time dimension of input data
    nlev = SIZE(REGR_4D,3)
    nt   = SIZE(REGR_4D,4)

    ! minlev is the first level where native grid differs from the
    ! reduced grid. All levels below are always mapped 1:1.
#if defined( GEOS_4 )
    minlev = 20
#elif defined( GEOS_5 ) || defined( MERRA ) || defined( GEOS_FP )
    minlev = 37
#endif

    !-------------------------------------------------------------------
    ! Check if there is a need to inflate or collapse the input data.
    !-------------------------------------------------------------------
    infl = .FALSE.
    clps = .FALSE.

    ! If input data has more vertical levels than output grid, we need
    ! to collapse the data.
    IF ( nlev > nz ) clps = .TRUE.

    ! If output grid has more vertical levels than input data, we only
    ! need to inflate the data if the input data extends beyond minlev.
    ! These data is ALWAYS assumed to be on the reduced GEOS-5 levels,
    ! i.e. data becomes never inflated onto a reduced GEOS-5 grid!
    IF ( nz > nlev ) THEN
       IF ( nlev > minlev .AND. nz /= 47 ) infl = .TRUE.
       
       ! Special case that input data has more than 47 levels and we 
       ! are on the native grid: never need to inflate in this case.
       IF ( nlev > 47 .AND. nz == 72 ) infl = .FALSE. 
    ENDIF

    ! now determine number of levels to be used in the output array
    IF ( clps .OR. infl ) THEN
       nout = nz
    ELSE
       nout = nlev
    ENDIF

    ! Make sure array in data container is allocated
    CALL FileData_ArrCheck( Lct%Dct%Dta, nx, ny, nout, nt, RC ) 

    ! verbose mode
    IF ( verb ) THEN
       WRITE(MSG,*) '   --> nlev, nz, nout: ', nlev, nz, nout
       CALL HCO_MSG(MSG)
       WRITE(MSG,*) '   --> Collapse? ', clps
       CALL HCO_MSG(MSG)
       WRITE(MSG,*) '   --> Inflate?  ', infl
       CALL HCO_MSG(MSG)
    ENDIF

    !-------------------------------------------------------------------
    ! Pass data to data container
    !-------------------------------------------------------------------
  
    ! Determine number of levels that can be passed 1:1, i.e. with no
    ! deflation/collapsing. 
    IF ( clps .OR. infl ) THEN
       NL = minlev - 1
    ELSE
       NL = nlev
    ENDIF

    ! Do for every time step
    DO T = 1, nt

       ! Pass data 1:1
       DO L = 1, NL
          Lct%Dct%Dta%V3(T)%Val(:,:,L) = REGR_4D(:,:,L,T)
       ENDDO !L

       ! Collapse data if needed
       IF ( clps ) THEN
#if defined( GEOS_4 )
          ! collapse native GEOS-4 onto reduced GEOS-4:
          IF ( nlev == 55 ) THEN
             ! two levels:
             CALL COLLAPSE( Lct, REGR_4D, 20, 20, 2, T )
             CALL COLLAPSE( Lct, REGR_4D, 21, 22, 2, T )
             CALL COLLAPSE( Lct, REGR_4D, 22, 24, 2, T )
             CALL COLLAPSE( Lct, REGR_4D, 23, 26, 2, T )
             ! four levels:
             CALL COLLAPSE( Lct, REGR_4D, 24, 28, 4, T )
             CALL COLLAPSE( Lct, REGR_4D, 25, 32, 4, T )
             CALL COLLAPSE( Lct, REGR_4D, 26, 36, 4, T )
             CALL COLLAPSE( Lct, REGR_4D, 27, 40, 4, T )
             CALL COLLAPSE( Lct, REGR_4D, 28, 44, 4, T )
             CALL COLLAPSE( Lct, REGR_4D, 29, 48, 4, T )
             CALL COLLAPSE( Lct, REGR_4D, 30, 52, 4, T )

          ! TODO: collapse GEOS-5 onto GEOS-4:
          ELSE
             MSG = 'Vertical interpolation from GEOS-5 to GEOS-4'// &
                   ' is currently not implemented: '//TRIM(Lct%Dct%cName)
             CALL HCO_ERROR ( MSG, RC )
             RETURN
          ENDIF
#else
          ! collapse native GEOS-5 onto reduced GEOS-5: 

          ! Collapse two levels (e.g. levels 39-40 into level 38):
          CALL COLLAPSE( Lct, REGR_4D, 37, 37, 2, T )
          CALL COLLAPSE( Lct, REGR_4D, 38, 39, 2, T )
          CALL COLLAPSE( Lct, REGR_4D, 39, 41, 2, T )
          CALL COLLAPSE( Lct, REGR_4D, 40, 43, 2, T )
          ! Collapse four levels:
          CALL COLLAPSE( Lct, REGR_4D, 41, 45, 4, T )
          CALL COLLAPSE( Lct, REGR_4D, 42, 49, 4, T )
          CALL COLLAPSE( Lct, REGR_4D, 43, 53, 4, T )
          CALL COLLAPSE( Lct, REGR_4D, 44, 57, 4, T )
          CALL COLLAPSE( Lct, REGR_4D, 45, 61, 4, T )
          CALL COLLAPSE( Lct, REGR_4D, 46, 65, 4, T )
          CALL COLLAPSE( Lct, REGR_4D, 47, 69, 4, T )
#endif
       ENDIF

       ! Inflate data if needed
       IF ( infl ) THEN
#if defined( GEOS_4 )
          ! TODO: GEOS-5 to GEOS-4:
          MSG = 'Vertical interpolation from GEOS-5 to GEOS-4'// &
                ' is currently not implemented: '//TRIM(Lct%Dct%cName)
          CALL HCO_ERROR ( MSG, RC )
          RETURN
#else
          ! GEOS-5 reduced to GEOS-5 native:

          ! Distribute over 2 levels (e.g. level 38 into 39-40):
          CALL INFLATE( Lct, REGR_4D, 37, 37, 2, T )
          CALL INFLATE( Lct, REGR_4D, 38, 39, 2, T )
          CALL INFLATE( Lct, REGR_4D, 39, 41, 2, T )
          CALL INFLATE( Lct, REGR_4D, 40, 43, 2, T )
          ! Distribute over 4 levels:
          CALL INFLATE( Lct, REGR_4D, 41, 45, 4, T )
          CALL INFLATE( Lct, REGR_4D, 42, 49, 4, T )
          CALL INFLATE( Lct, REGR_4D, 43, 53, 4, T )
          CALL INFLATE( Lct, REGR_4D, 44, 57, 4, T )
          CALL INFLATE( Lct, REGR_4D, 45, 61, 4, T )
          CALL INFLATE( Lct, REGR_4D, 46, 65, 4, T )
          CALL INFLATE( Lct, REGR_4D, 47, 69, 4, T )
#endif
       ENDIF

    ENDDO !T 
    ! Return w/ success
    CALL HCO_LEAVE ( RC ) 

  END SUBROUTINE ModelLev_Interpolate
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: COLLAPSE 
!
! !DESCRIPTION: Helper routine to collapse input levels onto the output grid.
! The input data is weighted by the grid box thicknesses defined on top of 
! this module.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE COLLAPSE ( Lct, REGR_4D, OutLev, InLev1, NLEV, T )
!
! !INPUT PARAMETERS:
!
    REAL(hp),         POINTER        :: REGR_4D(:,:,:,:)  ! 4D input data
    INTEGER,          INTENT(IN)     :: OutLev
    INTEGER,          INTENT(IN)     :: InLev1 
    INTEGER,          INTENT(IN)     :: NLEV
    INTEGER,          INTENT(IN)     :: T
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ListCont),   POINTER        :: Lct               ! HEMCO list container
!
! !REVISION HISTORY:
!  30 Dec 2014 - C. Keller   - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
    INTEGER               :: I, NZ, ILEV, TOPLEV
    REAL(hp)              :: THICK
    REAL(hp), POINTER     :: EDG(:) => NULL()
    REAL(hp), ALLOCATABLE :: WGT(:)

    !=================================================================
    ! COLLAPSE begins here
    !=================================================================

    ! Reset
    Lct%Dct%Dta%V3(T)%Val(:,:,OutLev) = 0.0_hp 

    ! Don't do anything if there are not enough levels in REGR_4D
    NZ = SIZE(REGR_4D,3)
    IF ( NZ < InLev1 ) RETURN

    ! Get maximum level to be used for pressure thickness calculations.
    ! This is one more than the number of levels to be used.
    TOPLEV = InLev1 + NLEV

    ! Get pointer to grid edges on the native input grid
#if defined( GEOS_4 )
    EDG => G4_EDGE_NATIVE(InLev1:TOPLEV)
#else
    EDG => G5_EDGE_NATIVE(InLev1:TOPLEV)
#endif

    ! Thickness of output level
    THICK = EDG(1) - EDG(1+NLEV)

    ! Get level weights
    ALLOCATE(WGT(NLEV))
    WGT = 0.0
    DO I = 1, NLEV
       WGT(I) = ( EDG(I) - EDG(I+1) ) / THICK
    ENDDO

    ! Pass levels to output data, one after each other
    Lct%Dct%Dta%V3(T)%Val(:,:,OutLev) = REGR_4D(:,:,InLev1,T) * WGT(1)
    DO I = 1, NLEV-1
       ILEV = InLev1 + I
       IF ( NZ < ILEV ) EXIT 
       Lct%Dct%Dta%V3(T)%Val(:,:,OutLev) = Lct%Dct%Dta%V3(T)%Val(:,:,OutLev) &
                                         + ( REGR_4D(:,:,ILEV,T) * WGT(I+1) )
    ENDDO

    ! Cleanup
    DEALLOCATE(WGT)
    EDG => NULL()

  END SUBROUTINE COLLAPSE 
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: INFLATE 
!
! !DESCRIPTION: Helper routine to inflate input levels onto the output grid.
! The values on the input data are evenly distributed amongst all output
! levels. 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE INFLATE ( Lct, REGR_4D, InLev, OutLev1, NLEV, T )
!
! !INPUT PARAMETERS:
!
    REAL(hp),         POINTER        :: REGR_4D(:,:,:,:)  ! 4D input data
    INTEGER,          INTENT(IN)     :: InLev
    INTEGER,          INTENT(IN)     :: OutLev1
    INTEGER,          INTENT(IN)     :: NLEV
    INTEGER,          INTENT(IN)     :: T
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ListCont),   POINTER        :: Lct               ! HEMCO list container
!
! !REVISION HISTORY:
!  30 Dec 2014 - C. Keller   - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
    INTEGER :: I, NZ, ILEV

    !=================================================================
    ! INFLATE begins here
    !=================================================================

    ! Get input data array
    NZ = SIZE(REGR_4D,3)

    ! Do for every output level
    DO I = 1, NLEV

       ! Current output level
       ILEV = OutLev1 + I - 1
      
       ! If input level is beyond vert. extent of input data, set output
       ! data to zero.
       IF ( InLev > NZ ) THEN
          Lct%Dct%Dta%V3(T)%Val(:,:,ILEV) = 0.0_hp 

       ! Otherwise, evenly distribute input data 
       ELSE
          Lct%Dct%Dta%V3(T)%Val(:,:,ILEV) = REGR_4D(:,:,InLev,T) 
       ENDIF
    ENDDO

  END SUBROUTINE INFLATE 
!EOC
END MODULE HCO_INTERP_MOD
