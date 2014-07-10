!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: hco_tidx_mod 
!
! !DESCRIPTION: Module HCO\_TIDX\_MOD contains routines and variables
! to organize and index data array time slices.\\
! The HEMCO data containers can hold multiple 2D or 3D data arrays 
! (aligned in a vector), providing a 4th dimension (time). During 
! emission calculation, only the vector element ('time slice') 
! representative for the given time will be used. Currently, the 
! following time slice intervals are supported:
!
!\begin{itemize}
!\item Constant: Only one time slice (default)
!\item Hourly: Between 2-24 time slices. These slices will be split 
! into even day bins, and are cycled through accordingly. The local 
! time is used to obtain the current valid index at each longitude.
!\item Hourly_gridded: As hourly, but uses the same time slice index
! across all longitudes (based upon UTC time).
!\item Weekdaily: Seven time slices, representing the days of the  
! week: Sun, Mon, ..., Sat. Uses local time.
!\item Weekdaily_gridded: As weekdaily, but uses utc time.
!\item Monthly: 12 time slices, representing the months of the year: 
! Jan, ..., Dec. Uses local time. 
!\end{itemize}
! \\ 
! The time slice cycling frequency is automatically determined during 
! creation of a data container - based upon the time stamp settings in 
! the HEMCO configuration file and the time information read from the 
! data (netCDF) file.\\
! Specifically, spatial uniform data (i.e. nx = ny = 1) is always 
! assigned the local time cycle intervals (hourly, weekdaily, or 
! monthly), while for gridded data, it's assumed that local-time
! effects are already taken into account and the gridded cycle intervals
! are thus applied.\\
! Structure AlltIDx organizes the indexing of the vector arrays. It
! contains the current valid time slice index for each of the above 
! defined time slice intervals. Each data container points to one of the
! elements of AlltIDx, according to the temporal dimension of its array. 
! The values of AlltIDx become update on every HEMCO time step. 
! 
! !NOTES: The current local time implementation assumes a regular grid,
! i.e. local time does not change with latitude! 
!
! !INTERFACE: 
!
      MODULE HCO_TIDX_MOD
!
! !USES:
!
      USE HCO_ERROR_MOD
      USE HCO_FILEDATA_MOD,  ONLY : TimeIdx

      IMPLICIT NONE
      PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
      PUBLIC :: tIDx_Assign
      PUBLIC :: tIDx_Update
      PUBLIC :: tIDx_Init
      PUBLIC :: tIDx_GetIndx
      PUBLIC :: tIDx_GetIndxVec
      PUBLIC :: tIDx_Cleanup
      PUBLIC :: HCO_GetPrefTimeAttr
      PUBLIC :: HCO_ExtractTime
!
! !REVISION HISTORY:
!  29 Dec 2012 - C. Keller - Initialization
!  22 Aug 2013 - C. Keller - Some time slice updates. 
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !PRIVATE MEMBER FUNCTIONS:
!
      ! The TimeIdxCollection object contains the pointers with the
      ! current valid vector indeces for all defined cycling intervals.
      TYPE ::  TimeIdxCollection
         INTEGER                :: nx           ! # of lons 
         TYPE(TimeIdx), POINTER :: CONSTANT
         TYPE(TimeIdx), POINTER :: HOURLY
         TYPE(TimeIdx), POINTER :: HOURLY_GRID 
         TYPE(TimeIdx), POINTER :: WEEKDAY
         TYPE(TimeIdx), POINTER :: WEEKDAY_GRID
         TYPE(TimeIdx), POINTER :: MONTHLY
      END TYPE TimeIdxCollection

      ! Types
      TYPE(TimeIdxCollection), POINTER       :: AlltIDx  => NULL()
!EOC
      CONTAINS
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: tIDx_Init
!
! !DESCRIPTION: Subroutine tIDx\_Init initializes the time slice index
! collection.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE tIDx_Init ( HcoState, RC )
!
! !USES:
!
      USE HCO_STATE_MOD,    ONLY : HCO_State
!
! !ARGUMENTS:
!
      TYPE(HCO_State), POINTER        :: HcoState  ! Hemco state 
      INTEGER,         INTENT(INOUT)  :: RC        ! Return code
!
! !REVISION HISTORY:
!  29 Dec 2012 - C. Keller - Initialization
!
!EOP
!------------------------------------------------------------------------------
!BOC

      !======================================================================
      ! tIDx_Init begins here!
      !======================================================================

      ! Enter
      CALL HCO_ENTER ( 'tIDx_Init (HCO_TIDX_MOD.F90)', RC )
      IF ( RC /= HCO_SUCCESS ) RETURN

      ! Allocate collection of time indeces 
      ALLOCATE ( AlltIDx ) 
      AlltIDx%nx = HcoState%NX

      ! Initialize the vectors holding the currently valid time slice 
      ! indices for the various cycle intervals. Only create longitude- 
      ! dependent time slice vector if cycle interval changes with 
      ! longitude, i.e. if it is time-zone dependent.

      ! ----------------------------------------------------------------
      ! "CONSTANT" => only one time slice, no cycling.
      ! ----------------------------------------------------------------
      ALLOCATE ( AlltIDx%CONSTANT )
      AlltIDx%CONSTANT%TypeID       = 0 
      AlltIDx%CONSTANT%TempRes      = "Constant"
      AlltIDx%CONSTANT%LonDependent = .FALSE.

      ALLOCATE( AlltIDx%CONSTANT%CurrIDx(1) )
      AlltIDx%CONSTANT%CurrIDx(:) = 1

      ! ----------------------------------------------------------------
      ! "HOURLY" => changes every hour, longitude-dependent
      ! ----------------------------------------------------------------
      ALLOCATE ( AlltIDx%HOURLY )
      AlltIDx%HOURLY%TypeID       = 24
      AlltIDx%HOURLY%TempRes      = "Hourly"
      AlltIDx%HOURLY%LonDependent = .TRUE.

      ALLOCATE( AlltIDx%HOURLY%CurrIDx(AlltIDx%NX) )
      AlltIDx%HOURLY%CurrIDx(:)   = 1

      ! ----------------------------------------------------------------
      ! "HOURLY_GRID" => changes every hour, longitude-independent
      ! ----------------------------------------------------------------
      ALLOCATE ( AlltIDx%HOURLY_GRID )
      AlltIDx%HOURLY_GRID%TypeID       = 241
      AlltIDx%HOURLY_GRID%TempRes      = "Hourly_Grid"
      AlltIDx%HOURLY_GRID%LonDependent = .FALSE.

      ALLOCATE( AlltIDx%HOURLY_GRID%CurrIDx(1) )
      AlltIDx%HOURLY_GRID%CurrIDx(:)   = 1

      ! ----------------------------------------------------------------
      ! "WEEKDAY" => changes every weekday, longitude-dependent
      ! ----------------------------------------------------------------
      ALLOCATE ( AlltIDx%WEEKDAY )
      AlltIDx%WEEKDAY%TypeID       = 7
      AlltIDx%WEEKDAY%TempRes      = "Weekday"
      AlltIDx%WEEKDAY%LonDependent = .TRUE.

      ALLOCATE( AlltIDx%WEEKDAY%CurrIDx(AlltIDx%NX) )
      AlltIDx%WEEKDAY%CurrIDx(:)   = 1

      ! ----------------------------------------------------------------
      ! "WEEKDAY_GRID" => changes every weekday, longitude-independent
      ! ----------------------------------------------------------------
      ALLOCATE ( AlltIDx%WEEKDAY_GRID )
      AlltIDx%WEEKDAY_GRID%TypeID       = 71
      AlltIDx%WEEKDAY_GRID%TempRes      = "Weekday_Grid"
      AlltIDx%WEEKDAY_GRID%LonDependent = .FALSE.

      ALLOCATE( AlltIDx%WEEKDAY_GRID%CurrIDx(1) )
      AlltIDx%WEEKDAY_GRID%CurrIDx(:)   = 1

      ! ----------------------------------------------------------------
      ! "MONTHLY" => changes every month, longitude-dependent
      ! ----------------------------------------------------------------
      ALLOCATE ( AlltIDx%MONTHLY )
      AlltIDx%MONTHLY%TypeID       = 12
      AlltIDx%MONTHLY%TempRes      = "Monthly"
      AlltIDx%MONTHLY%LonDependent = .TRUE.

      ALLOCATE( AlltIDx%MONTHLY%CurrIDx(AlltIDx%NX) )
      AlltIDx%MONTHLY%CurrIDx(:)   = 1

      ! Return w/ success
      CALL HCO_LEAVE ( RC )

      END SUBROUTINE tIDx_Init
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: tIDx_Set
!
! !DESCRIPTION: Subroutine tIDx\_Set linkes the passed TimeIDx type to the 
! corresponding element in the TimeIdx collection, according to the given 
! TypeID. TypeID can be one of the following:
!\begin{itemize}
!\item 1: Constant 
!\item 24: Hourly 
!\item 241: Hourly_Grid 
!\item 7: Weekday 
!\item 71: Weekday_Grid
!\item 12: Monthly 
!\end{itemize}
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE tIDx_Set ( ctIDx, TypeID )
!
! !USES:
!
! !ARGUMENTS:
!
      TYPE(TimeIdx), POINTER     :: ctIDx   ! container TimeIDx
      INTEGER,       INTENT(IN)  :: TypeID  ! type ID
!
! !REVISION HISTORY:
!  29 Dec 2012 - C. Keller - Initialization
!
!EOP
!------------------------------------------------------------------------------
!BOC

      !======================================================================
      ! tIDx_Set begins here!
      !======================================================================

      ! Set pointer according to TempRes
      SELECT CASE ( TypeID )

         CASE ( 1 )
            ctIDx => AlltIDx%CONSTANT

         CASE ( 24 )
            ctIDx => AlltIDx%HOURLY

         CASE ( 241 )
            ctIDx => AlltIDx%HOURLY_GRID

         CASE ( 7 )
            ctIDx => AlltIDx%WEEKDAY

         CASE ( 71 )
            ctIDx => AlltIDx%WEEKDAY_GRID

         CASE ( 12 )
            ctIDx => AlltIDx%MONTHLY

         CASE DEFAULT
            ctIDx => AlltIDx%CONSTANT

      END SELECT

      END SUBROUTINE tIDx_Set
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: tIDx_Cleanup
!
! !DESCRIPTION: Subroutine tIDx\_Cleanup deallocates the time slice
! index collection. 
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE tIDx_Cleanup
!
! !REVISION HISTORY:
!  29 Dec 2012 - C. Keller - Initialization
!
!EOP
!------------------------------------------------------------------------------
!BOC

      !======================================================================
      ! tIDx_Cleanup begins here!
      !======================================================================

      IF ( ASSOCIATED( AlltIDx ) ) THEN

         IF ( ASSOCIATED(AlltIDx%CONSTANT) ) THEN
            IF(ASSOCIATED(AlltIDx%CONSTANT%CurrIDx)) &
               DEALLOCATE(AlltIDx%CONSTANT%CurrIDx) 
            DEALLOCATE(AlltIDx%CONSTANT) 
         ENDIF

         IF ( ASSOCIATED(AlltIDx%HOURLY) ) THEN
            IF(ASSOCIATED(AlltIDx%HOURLY%CurrIDx)) &
               DEALLOCATE(AlltIDx%HOURLY%CurrIDx) 
            DEALLOCATE(AlltIDx%HOURLY) 
         ENDIF

         IF ( ASSOCIATED(AlltIDx%HOURLY_GRID) ) THEN
            IF(ASSOCIATED(AlltIDx%HOURLY_GRID%CurrIDx)) &
               DEALLOCATE(AlltIDx%HOURLY_GRID%CurrIDx) 
            DEALLOCATE(AlltIDx%HOURLY_GRID) 
         ENDIF

         IF ( ASSOCIATED(AlltIDx%WEEKDAY) ) THEN
            IF(ASSOCIATED(AlltIDx%WEEKDAY%CurrIDx)) &
               DEALLOCATE(AlltIDx%WEEKDAY%CurrIDx) 
            DEALLOCATE(AlltIDx%WEEKDAY) 
         ENDIF

         IF ( ASSOCIATED(AlltIDx%WEEKDAY_GRID) ) THEN
            IF(ASSOCIATED(AlltIDx%WEEKDAY_GRID%CurrIDx)) &
               DEALLOCATE(AlltIDx%WEEKDAY_GRID%CurrIDx) 
            DEALLOCATE(AlltIDx%WEEKDAY_GRID) 
         ENDIF
  
         IF ( ASSOCIATED(AlltIDx%MONTHLY) ) THEN
            IF(ASSOCIATED(AlltIDx%MONTHLY%CurrIDx)) &
               DEALLOCATE(AlltIDx%MONTHLY%CurrIDx) 
            DEALLOCATE(AlltIDx%MONTHLY) 
         ENDIF
  
         ! Also deallocate AlltIDx pointer 
         DEALLOCATE( AlltIDx )
      ENDIF
      AlltIDx => NULL()

      END SUBROUTINE tIDx_Cleanup
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: tIDx_Update
!
! !DESCRIPTION: Subroutine tIDx\_Update updates the current valid 
! index values for every time slice type.
!\\
! !INTERFACE:
!
      SUBROUTINE tIDx_Update ( am_I_Root, RC )
!
! !USES:
!
      USE HCO_CLOCK_MOD, ONLY : HcoClock_Get, HcoClock_GetLocal
!
! !ARGUMENTS:
!
      LOGICAL,         INTENT(IN   ) :: am_I_Root
      INTEGER,         INTENT(INOUT) :: RC 
!
! !REVISION HISTORY:
!  29 Dec 2012 - C. Keller - Initialization
!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER            :: I, HH, WD, iMM, iWD
      REAL(sp)           :: iHH
      CHARACTER(LEN=255) :: MSG, LOC

      ! testing only
      LOGICAL, SAVE      :: FIRST = .TRUE.

      !======================================================================
      ! tIDx_Update begins here!
      !======================================================================

      ! Enter
      LOC = 'tIDx_Update (HCO_TIDX_MOD.F90)'

      ! ----------------------------------------------------------------
      ! Set longitude independent indeces
      ! ----------------------------------------------------------------

      ! Get current times
      CALL HcoClock_Get( cH = HH, cWeekday = WD, RC=RC )
      IF ( RC/= HCO_SUCCESS ) RETURN

      ! CONSTANT:
      ! --> nothing to be done

      ! HOURLY_GRID:
      ! Gridded hourly data is assumed to be already adjusted for
      ! local time effects, hence just point to the time slice
      ! of current UTC time. Add one since hour starts at 0.
      AlltIDx%HOURLY_GRID%CurrIDx(1) = HH + 1

      ! WEEKDAY_GRID:
      ! For gridded weekday factors, just use the UTC slice. Add
      ! one since weekday start at 0.
      AlltIDx%WEEKDAY_GRID%CurrIDx(1) = WD + 1

      ! ----------------------------------------------------------------
      ! Set longitude dependent time indeces 
      ! ----------------------------------------------------------------
      DO I = 1, AlltIDx%NX

         ! Get local times
         CALL HcoClock_GetLocal( I, cMM=iMM, cH=iHH, cWeekday=iWD, RC=RC ) 
         IF ( RC/= HCO_SUCCESS ) RETURN

         ! HOURLY:
         ! These are the pointers to the hourly indices. Hourly points
         ! to the hourly time slice representative for the LOCAL time at
         ! longitude ii.
         AlltIDx%HOURLY%CurrIDx(I) = FLOOR(iHH) + 1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!          TEMPORARY FIX FOR CONSISTENCY W/ FORMER VERSION
!
         ! WEEKDAY:
         ! For non-gridded factors, take into account local time:
!         AlltIDx%WEEKDAY%CurrIDx(I) = iWD + 1
         AlltIDx%WEEKDAY%CurrIDx(I) = WD + 1
         AlltIDx%WEEKDAY%LonDependent = .FALSE.
         IF ( FIRST .AND. I == 1 .AND. am_I_Root ) THEN
            MSG =  'Constant weekday used, needs to be fixed!'
            CALL HCO_WARNING( MSG, RC, THISLOC=LOC )
         ENDIF 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         ! MONTHLY data

         ! Set month index accordingly
         AlltIDx%MONTHLY%CurrIDx(I) = iMM 

      ENDDO !I

      ! Adjust first flag
      FIRST = .FALSE. 

      ! Leave w/ success
      RC = HCO_SUCCESS 

      END SUBROUTINE tIDx_Update
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: tIDx_GetIndx 
!
! !DESCRIPTION: Function tIDx\_GetIndx returns the current active time 
! slice index for the given file data object and longitude (by index).
!\\
!\\
! !INTERFACE:
!
      FUNCTION tIDx_GetIndx( Dta, LonIndx ) RESULT ( Indx )
!
! !USES:
!
      USE HCO_FILEDATA_MOD, ONLY : FileData 
!
! !INPUT PARAMETERS: 
!
      TYPE(FileData), POINTER       :: Dta     ! File data object 
      INTEGER,        INTENT(IN)    :: LonIndx ! Longitude index 
!
! !RETURN VALUE:
!
      INTEGER             :: Indx 
!
! !REVISION HISTORY: 
!  13 Jan 2014 - C. Keller - Initial version 
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      REAL(dp)     :: frac

      !-----------------------------------
      ! tIDx_GetIndx begins here! 
      !-----------------------------------

      ! For longitude dependent time slices, point to time index
      ! corresponding to the given longitude. 
      IF ( Dta%tIDx%LonDependent ) THEN 
         indx = Dta%tIDx%CurrIDx(LonIndx)

      ! Time independent time slices: there is only one index.
      ELSE
         indx = Dta%tIDx%CurrIDx(1) 
      ENDIF

      ! For hourly data with less than 24 time slices, i.e. time
      ! intervals of more than 1 hour, map 24 hour index onto the
      ! reduced time slice. 
      ! For example, for 3-hourly data (8 vector elements), this will 
      ! return index 1 for hours 0-2am, index 2 for 3-5am, etc.
      IF ( Dta%DeltaT > 1 .AND. Dta%DeltaT < 24 ) THEN
         frac = DBLE(indx) / DBLE(Dta%DeltaT) 
         indx = CEILING( frac )
      ENDIF

      END FUNCTION tIDx_GetIndx
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: tIDx_GetIndxVec 
!
! !DESCRIPTION: Function tIDx\_GetIndxVec returns a vector with the 
! current active time slice indeces for every longitude. 
!\\
!\\
! !INTERFACE:
!
      FUNCTION tIDx_GetIndxVec( Dta, nI ) RESULT ( IndxVec )
!
! !USES:
!
      USE HCO_FILEDATA_MOD, ONLY : FileData
!
! !INPUT PARAMETERS: 
!
      TYPE(FileData), POINTER       :: Dta    ! List container 
      INTEGER,        INTENT(IN)    :: nI         ! # of lons 
!
! !RETURN VALUE:
!
      INTEGER             :: IndxVec(nI)
!
! !REVISION HISTORY: 
!  13 Jan 2014 - C. Keller - Initial version 
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER :: I 

      !-----------------------------------
      ! tIDx_GetIndxVec begins here! 
      !-----------------------------------

      ! For longitude dependent time slices, get index for every longitude 
      IF ( Dta%tIDx%LonDependent ) THEN 
         DO I = 1, nI
            IndxVec(I) = tIDx_GetIndx( Dta, I )
         ENDDO 

      ! For longitude independent data: all indeces are the same, so just
      ! use first one.
      ELSE
         I          = tIDx_GetIndx( Dta, 1 )
         IndxVec(:) = I   
      ENDIF

      END FUNCTION tIDx_GetIndxVec
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: tIDx_Assign 
!
! !DESCRIPTION: Subroutine tIDx\_Assign assigns the time index
! pointer of the file data object of the passed list container 
! to the corresponding time index element of the time index 
! collection. The time slice cycle interval is determined from 
! the number of time slices (length of the data array vector) and 
! the time difference (deltaT) between them.\\ 
! Note: Scale factors read directly from the HEMCO configuration file
! become their delta T assigned upon reading from file (see 
! subroutine Register\_Scal in hco\_config\_mod.F90).
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE tIDx_Assign ( HcoState, Lct, RC ) 
!
! !USES:
!
      USE HCO_STATE_MOD,    ONLY : HCO_State
      USE HCO_DATACONT_MOD, ONLY : ListCont
!
! !INPUT PARAMETERS: 
!
      TYPE(HCO_State), POINTER       :: HcoState  ! Hemco state 
      TYPE(ListCont),  POINTER       :: Lct   ! List container
      INTEGER,         INTENT(INOUT) :: RC        ! Success or failure?
!
! !REVISION HISTORY: 
!  13 Jan 2014 - C. Keller - Initial version 
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER            :: nx, ny, nt, ntexp, dt 
      INTEGER            :: cTypeID
      CHARACTER(LEN=255) :: MSG

      !-----------------------------------
      ! tIDx_Assign begins here! 
      !-----------------------------------

      ! Enter
      CALL HCO_ENTER ( 'tIDx_Assign (HCO_TIDX_MOD.F90)', RC )
      IF ( RC /= HCO_SUCCESS ) RETURN

      ! Check if already done
      IF ( ASSOCIATED( Lct%Dct%Dta%tIDx ) ) THEN
         CALL HCO_LEAVE( RC )
         RETURN
      ENDIF 

      ! Get array dimensions
      nt = Lct%Dct%Dta%nt

      ! Do the following only for defined data arrays
      IF ( nt > 0 ) THEN

         ! -------------------------------------------------------------
         ! Extract data array dimensions and delta t [h] between time 
         ! slices
         ! -------------------------------------------------------------
         dt = Lct%Dct%Dta%DeltaT
         IF ( Lct%Dct%Dta%SpaceDim <= 2 ) THEN
            nx = SIZE(Lct%Dct%Dta%V2(1)%Val,1)
            ny = SIZE(Lct%Dct%Dta%V2(1)%Val,2)
         ELSE
            nx = SIZE(Lct%Dct%Dta%V3(1)%Val,1)
            ny = SIZE(Lct%Dct%Dta%V3(1)%Val,2)
         ENDIF

         ! -------------------------------------------------------------
         ! Auto detection of temporal resolution. 
         ! -------------------------------------------------------------

         ! No time variation
         IF ( nt == 1 ) THEN
            cTypeID = 1 

         ! -------------------------------------------------------------
         ! Daily varying data. 
         ! For now, these arrays must contain seven time slices, i.e. 
         ! one for every weekday (starting w/ sunday!). 
         ! -------------------------------------------------------------
         ELSEIF ( nt == 7 ) THEN 

            ! Sanity check if dt = 24 hours
            IF ( dt /= 24 ) THEN
               MSG = '7 time slices but delta t is not 24 hours!' // &
                  TRIM(Lct%Dct%cName)
               CALL HCO_ERROR ( MSG, RC )
               RETURN
            ENDIF

            ! Check if data is gridded or not, set TempRes attribute
            ! accordingly.
            IF ( nx == 1 .AND. ny == 1 ) THEN
               cTypeID = 7 
            ELSE
               cTypeID = 71 
            ENDIF

         ! -------------------------------------------------------------
         ! Monthly varying data
         ! This must be a scalar. Gridded monthly data should be read 
         ! slice by slice!
         ! -------------------------------------------------------------
         ELSEIF ( nt == 12 .AND. dt /= 2 ) THEN 
            IF ( nx == 1 .AND. ny == 1 ) THEN
               cTypeID = 12 
            ELSE
               MSG = 'Monthly data must not be gridded:' // &
                  TRIM(Lct%Dct%cName)
               CALL HCO_ERROR ( MSG, RC )
               RETURN
            ENDIF

         ! -------------------------------------------------------------
         ! Hourly varying data. 
         ! Hourly time slices can represent more than just one hour, 
         ! e.g. an array can hold only 8 time slices, each representing 
         ! 3-hourly intervals of the day. The only constraint is that 
         ! the number of time slices have to be a factor of 24!
         ! -------------------------------------------------------------
         ELSEIF ( dt < 24 .AND. dt > 0 ) THEN 

            ! Check if we can evenly split up data within a day 
            IF ( MOD(24,dt) /= 0 ) THEN
               MSG = 'Cannot properly split up hourly data!' // &
                     TRIM(Lct%Dct%cName)
               CALL HCO_ERROR ( MSG, RC )
               RETURN
            ENDIF

            ! Also check if nt is indeed of expected size
            ntexp = 24/dt
            IF ( ntexp /= nt ) THEN
               MSG = 'Wrong delta t and/or number of time slices!' // &
                     TRIM(Lct%Dct%cName)
               CALL HCO_ERROR ( MSG, RC )
               RETURN
            ENDIF

            ! Check if data is gridded or not, set TempRes attribute
            ! accordingly.
            IF ( nx == 1 .AND. ny == 1 ) THEN
               cTypeID = 24 
            ELSE
               cTypeID = 241 
            ENDIF

         ! -------------------------------------------------------------
         ! Return w/ error otherwise
         ! -------------------------------------------------------------
         ELSE
            MSG = 'Invalid time slice for field ' // &
                  TRIM(Lct%Dct%cName)
            CALL HCO_ERROR ( MSG, RC )
            RETURN
         ENDIF

         ! -------------------------------------------------------------
         ! Establish the appropriate pointer for the
         ! 4th dimension (temporal resolution) of the field array
         ! -------------------------------------------------------------
         CALL tIDx_Set ( Lct%Dct%Dta%tIDx, cTypeID ) 
     
      ENDIF

      ! Leave w/ success
      CALL HCO_LEAVE ( RC )

      END SUBROUTINE tIDx_Assign 
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCO_GetPrefTimeAttr
!
! !DESCRIPTION: Subroutine HCO\_GetPrefTimeAttr returns the preferred time
! reading attributes for a given field, based upon the specs set in the HEMCO 
! configuration file and the current simulation date. This routine is used 
! to select the proper time slice to be read at the given time.
! The time reading attributes are set to the value that is closest (or equal) 
! to the current datetime but within the range provided in the configuration 
! file. If the year, month, or day is not specified, the current simulation 
! date values are taken. For unspecified hours, a value of -1 is returned 
! (this allows to read all hourly slices at once). 
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE HCO_GetPrefTimeAttr ( Lct,    readYr, readMt, &
                                       readDy, readHr, RC       ) 
!
! !USES:
!
      USE HCO_DATACONT_MOD, ONLY : ListCont
      USE HCO_CLOCK_MOD,    ONLY : HcoClock_Get
!
! !INPUT PARAMETERS: 
!
      TYPE(ListCont),  POINTER       :: Lct       ! List container
      INTEGER,         INTENT(  OUT) :: readYr    ! preferred year 
      INTEGER,         INTENT(  OUT) :: readMt    ! preferred month 
      INTEGER,         INTENT(  OUT) :: readDy    ! preferred day
      INTEGER,         INTENT(  OUT) :: readHr    ! preferred hour 
      INTEGER,         INTENT(INOUT) :: RC        ! Success or failure?
!
! !REVISION HISTORY: 
!  13 Jan 2014 - C. Keller - Initial version 
!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER   :: cYr, cMt, cDy, cHr 

      !-----------------------------------
      ! HCO_GetPrefTimeAttr begins here! 
      !-----------------------------------

      ! Init 
      RC = HCO_SUCCESS

      ! Get current time
      CALL HcoClock_Get ( cYYYY = cYr, cMM = cMt,        &
                          cDD   = cDy, cH  = cHr, RC = RC ) 
      IF ( RC /= HCO_SUCCESS ) RETURN 

      ! ------------------------------------------------------------- 
      ! If CycleFlag is set to 2 or 3, the preferred datetime is 
      ! always the current date.
      ! ------------------------------------------------------------- 
      IF ( Lct%Dct%Dta%CycleFlag > 1 ) THEN
         readYr = cYr
         readMt = cMt
         readDy = cDy
         IF ( Lct%Dct%Dta%ncHrs(1) == -1 ) THEN
            readHr = -1
         ELSE
            readHr = cHr
         ENDIF

         ! Don't need below
         RETURN
      ENDIF 

      ! ------------------------------------------------------------- 
      ! If CycleFlag is set to 1, select the time attributes
      ! as specified in the configuration file and closest to 
      ! current simulation date.
      ! ------------------------------------------------------------- 

      ! Year
      IF ( Lct%Dct%Dta%ncYrs(1) == Lct%Dct%Dta%ncYrs(2) ) THEN
         readYr = Lct%Dct%Dta%ncYrs(1)
      ELSE
         readYr = MIN( MAX(cYr,Lct%Dct%Dta%ncYrs(1)), &
                       Lct%Dct%Dta%ncYrs(2)          )
      ENDIF

      ! Month 
      IF ( Lct%Dct%Dta%ncMts(1) == Lct%Dct%Dta%ncMts(2) ) THEN
         readMt = Lct%Dct%Dta%ncMts(1)
      ELSE
         readMt = MIN( MAX(cMt,Lct%Dct%Dta%ncMts(1)), &
                       Lct%Dct%Dta%ncMts(2)          )
      ENDIF

      ! Day
      IF ( Lct%Dct%Dta%ncDys(1) == Lct%Dct%Dta%ncDys(2) ) THEN
         readDy = Lct%Dct%Dta%ncDys(1)
      ELSE
         readDy = MIN( MAX(cDy,Lct%Dct%Dta%ncDys(1)), &
                       Lct%Dct%Dta%ncDys(2)          )
      ENDIF

      ! Hour
      IF ( Lct%Dct%Dta%ncHrs(1) == Lct%Dct%Dta%ncHrs(2) ) THEN
         readHr = Lct%Dct%Dta%ncHrs(1)
      ELSE
         readHr = MIN( MAX(cHr,Lct%Dct%Dta%ncHrs(1)), &
                       Lct%Dct%Dta%ncHrs(2)          )
      ENDIF

      ! Don't allow invalid entries for years, months or days, i.e.
      ! force readYr, readMt and readDy to be positive!
      ! readHr can be negative, in which case all hourly fields will
      ! become read!
      if ( readYr < 0 ) readYr = cYr 
      if ( readMt < 0 ) readMt = cMt
      if ( readDy < 0 ) readDy = cDy

      END SUBROUTINE HCO_GetPrefTimeAttr 
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCO_ExtractTime 
!
! !DESCRIPTION: Subroutine HCO_ExtractTime extracts the time stamp (ranges)
! from the passed character string. The string is expected to be in
! format Yr/Mt/Dy/Hr. Valid values for Yr, Mt, Dy and Hr are:
! (1) Range of values, separated by - sign: e.g. 2000-2010.
! (2) Single value: 2000
! (3) Asterik as wildcard character: *
! (4) Placeholder for start date, i.e. YYYY, MM, DD, HH.
! The extracted time stamp is written into the arrays ncYrs, ncMts,
! ncDys and ncHrs of the passed data structure Dta.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE HCO_ExtractTime ( CharStr, Dta, RC ) 
!
! !USES:
!
      USE CHARPAK_MOD,        ONLY : STRSPLIT
      USE HCO_FILEDATA_MOD,   ONLY : FileData 
!
! !INPUT PARAMETERS: 
!
      CHARACTER(LEN=*), INTENT(IN   )       :: CharStr 
      TYPE(FileData),   POINTER             :: Dta 
!
! !OUTPUT PARAMETERS:
!
      INTEGER,          INTENT(INOUT)       :: RC 
! 
! !REVISION HISTORY: 
!  18 Sep 2013 - C. Keller - Initial version (update) 
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!

      INTEGER               :: I, I0, I1, N
      INTEGER               :: TimeVec(8)
      CHARACTER(LEN=255)    :: SUBSTR(255), DATERNG(255), LOC

      !=================================================================
      ! HCO_ExtractTime begins here!
      !=================================================================

      ! Enter
      LOC = 'HCO_ExtractTime (HCO_CONFIG_MOD.F90)'

      ! Init
      TimeVec(:) = -1

      ! Extract strings to be translated into integers 
      CALL STRSPLIT( CharStr, HCO_SEP(), SUBSTR, N )
      IF ( N > 4 ) THEN
         CALL HCO_ERROR( 'Too many substrings!', RC, THISLOC=LOC )
         RETURN
      ENDIF   

      ! Extract year, month, day, and hour range from the four
      ! substrings

      ! Do for each substring:
      DO I = 1,4

         ! Indices in TimeVec:
         I0 = (I-1) * 2 + 1
         I1 = I0 + 1

         ! For wildcard character, set lower and upper limit both to -1.
         ! In this case, the whole time slice will be read into file!
         IF ( TRIM(SUBSTR(I)) == TRIM(HCO_WILDCARD() ) ) THEN
            TimeVec(I0:I1) = -1 

         ELSEIF ( TRIM(SUBSTR(I)) == 'YYYY' .OR. &
                  TRIM(SUBSTR(I)) == 'MM'   .OR. &
                  TRIM(SUBSTR(I)) == 'DD'   .OR. &
                  TRIM(SUBSTR(I)) == 'HH'         ) THEN
            TimeVec(I0:I1) = -999

         ! Otherwise, check if date range if given and set lower and
         ! upper bound accordingly.
         ELSE
            CALL STRSPLIT( SUBSTR(I), '-', DATERNG, N )

            ! If range is given:
            IF ( N == 2 ) THEN
               READ( DATERNG(1), * ) TimeVec(I0)
               READ( DATERNG(2), * ) TimeVec(I1)

            ! Use same value if only one value is given:
            ELSEIF ( N == 1 ) THEN
               READ( DATERNG(1), * ) TimeVec(I0)
               TimeVec(I1) = TimeVec(I0)
            ELSE
               CALL HCO_ERROR( 'Cannot extract time stamp!', &
                               RC, THISLOC=LOC )
               RETURN
            ENDIF
         ENDIF
      ENDDO

      ! Pass to list container
      Dta%ncYrs = TimeVec(1:2)
      Dta%ncMts = TimeVec(3:4)
      Dta%ncDys = TimeVec(5:6)
      Dta%ncHrs = TimeVec(7:8)

      ! Leave w/ success
      RC = HCO_SUCCESS

      END SUBROUTINE HCO_ExtractTime 
!EOC
      END MODULE HCO_TIDX_MOD
!EOM
