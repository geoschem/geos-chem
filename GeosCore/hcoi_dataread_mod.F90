!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: hco_dataread_mod.F90 
!
! !DESCRIPTION: Module HCO\_DATAREAD\_MOD controls data processing (file
! reading, unit conversion, regridding) for HEMCO.
!\\
!\\
! !INTERFACE: 
!
MODULE HCOI_DATAREAD_MOD
!
! !USES:
!
  USE HCO_ERROR_MOD
  USE HCO_TYPE_MOD, ONLY : Hco_State, RdCont
  USE HCO_TIME_MOD, ONLY : HcoClock 

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC  :: HCOI_DATAREAD
!
! !PRIVATE MEMBER FUNCTIONS:
!
!
! !MODULE INTERFACES:
!
! !REVISION HISTORY:
!  22 Aug 2013 - C. Keller   - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
  CONTAINS
!EOC
#if defined(ESMF_)
!-----------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group     !
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCOI_DATAREAD 
!
! !DESCRIPTION: Interface between HEMCO and ESMF reading routines.  
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE HCOI_DATAREAD ( am_I_Root,  HcoState, &
                                 ThisRdCont, Clock, RC )
!
! !USES:
!
      USE ESMF
      USE MAPL_MOD
# include "MAPL_Generic.h"
!
! !INPUT/OUTPUT PARAMETERS:
!
      LOGICAL,          INTENT(IN   )  :: am_I_Root
      TYPE(HCO_State),  POINTER        :: HcoState
      TYPE(RdCont),     POINTER        :: ThisRdCont 
      TYPE(HcoClock),   POINTER        :: Clock 
      INTEGER,          INTENT(INOUT)  :: RC
!
! !REVISION HISTORY:
!  28 Aug 2013 - C. Keller - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
! 
! !ROUTINE ARGUMENTS:
!
      INTEGER                    :: II, JJ, LL, TT
      INTEGER                    :: I, J, L, T
      INTEGER                    :: STAT
      REAL*4, POINTER            :: Ptr3D(:,:,:)   => NULL() 
      REAL*4, POINTER            :: Ptr2D(:,:)     => NULL() 
      TYPE(ESMF_State), POINTER  :: IMPORT         => NULL()
      CHARACTER(LEN=255)         :: LOC

      !=================================================================
      ! HCOI_DATAREAD begins here
      !=================================================================

      ! For error handling
      LOC = 'HCOI_DATAREAD - ESMF (HCOI_DATAREAD_MOD.F90)'

      ! Point to ESMF IMPORT object
      IMPORT => HcoState%IMPORT

      !-----------------------------------------------------------------
      ! Read 3D data from ESMF 
      !-----------------------------------------------------------------
      IF ( TRIM(ThisRdCont%ESMF_Dim) == 'xyz' ) THEN

         ! Get data
         CALL MAPL_GetPointer ( IMPORT, Ptr3D, &
                                TRIM(ThisRdCont%cName), RC=STAT )

         ! Check for MAPL error
         IF( MAPL_VRFY(STAT,LOC,1) ) THEN
            CALL HCO_ERROR ( 'Cannot get xyz pointer', LOC, RC ) 
            RETURN
         ENDIF

         ! Get array dimensions 
         II = SIZE(Ptr3D,1)
         JJ = SIZE(Ptr3D,2) 
         LL = SIZE(Ptr3D,3)
         TT = 1 

         ! Allocate output array if not yet defined 
         IF ( .NOT. ASSOCIATED ( ThisRdCont%Array ) ) THEN 
            ALLOCATE(ThisRdCont%Array(II,JJ,LL,TT))   
            ThisRdCont%Array = 0d0

         ! Check output dimensions otherwise 
         ELSE
            IF ( ( SIZE(ThisRdCont%Array,1) /= II ) .OR. &
                 ( SIZE(ThisRdCont%Array,2) /= JJ ) .OR. & 
                 ( SIZE(ThisRdCont%Array,3) /= LL ) .OR. & 
                 ( SIZE(ThisRdCont%Array,4) /= TT )       ) THEN 

               CALL HCO_ERROR ( 'Wrong ThisRdCont%Array dimension', LOC, RC ) 
               RETURN
            ENDIF
         ENDIF

         ! Copy data and cast to real*8
         DO T = 1, TT 
         DO L = 1, LL
         DO J = 1, JJ 
         DO I = 1, II
            ThisRdCont%Array(I,J,L,T) = Ptr3D(I,J,L)
         ENDDO   
         ENDDO   
         ENDDO   
         ENDDO   

      !-----------------------------------------------------------------
      ! Read 2D data from ESMF 
      !-----------------------------------------------------------------
      ELSEIF ( TRIM(ThisRdCont%ESMF_Dim) == 'xy' ) THEN

         ! Get data
         CALL MAPL_GetPointer ( IMPORT, Ptr2D, &
                                TRIM(ThisRdCont%cName), RC=STAT )


         ! Check for MAPL error 
         IF( MAPL_VRFY(STAT,LOC,2) ) THEN
            CALL HCO_ERROR ( 'Cannot get xy pointer', LOC, RC ) 
            RETURN
         ENDIF

         ! Get array dimensions 
         II = SIZE(Ptr2D,1)
         JJ = SIZE(Ptr2D,2) 
         LL = 1 
         TT = 1 

         ! Allocate output array if not yet defined 
         IF ( .NOT. ASSOCIATED ( ThisRdCont%Array ) ) THEN 
            ALLOCATE(ThisRdCont%Array(II,JJ,LL,TT))   
            ThisRdCont%Array = 0d0

         ! Check output dimensions otherwise 
         ELSE
            IF ( ( SIZE(ThisRdCont%Array,1) /= II ) .OR. &
                 ( SIZE(ThisRdCont%Array,2) /= JJ ) .OR. & 
                 ( SIZE(ThisRdCont%Array,3) /= LL ) .OR. & 
                 ( SIZE(ThisRdCont%Array,4) /= TT )       ) THEN 

               ! Return with error
               CALL HCO_ERROR ( 'Wrong ThisRdCont%Array dimension', LOC, RC ) 
               RETURN
            ENDIF
         ENDIF
 
         ! Copy data and cast to real*8
         DO T = 1, TT 
         DO L = 1, LL
         DO J = 1, JJ
         DO I = 1, II
               ThisRdCont%Array(I,J,L,T) = Ptr2D(I,J)
         ENDDO
         ENDDO
         ENDDO
         ENDDO

      ENDIF  
 
      !-----------------------------------------------------------------
      ! Cleanup and leave 
      !-----------------------------------------------------------------
      Ptr3D  => NULL()
      Ptr2D  => NULL()
      IMPORT => NULL()   

      ! Return w/ success
      RC = HCO_SUCCESS

      END SUBROUTINE HCOI_DATAREAD
!EOC
#else
!-----------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group     !
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCOI_DATAREAD 
!
! !DESCRIPTION: Reads a netCDF file and returns the regridded array in proper
! units. This is an interface to the GEOS-Chem reading and regridding 
! routines. 
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE HCOI_DATAREAD ( am_I_Root,  HcoState, &
                                 ThisRdCont, Clock, RC )
!
! !USES:
!
!      USE NCDF_MOD,           ONLY : NC_READ
!      USE NCDF_MOD,           ONLY : NC_READ_GRID
      USE UNITCONV_MOD,       ONLY : CHANGE_UNITS
      USE REGRID_A2A_MOD,     ONLY : MAP_A2A
!
! !INPUT/OUTPUT PARAMETERS:
!
      LOGICAL,          INTENT(IN   )  :: am_I_Root
      TYPE(HCO_State),  POINTER        :: HcoState
      TYPE(RdCont),     POINTER        :: ThisRdCont
      TYPE(HcoClock),   POINTER        :: Clock 
      INTEGER,          INTENT(INOUT)  :: RC
!
! !REVISION HISTORY:
!  13 Mar 2013 - C. Keller - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
! 
! !ROUTINE ARGUMENTS:
!
      INTEGER               :: TRCID, YYYY, MM, DD, HH
      INTEGER               :: NLON, NLAT, NLEV, NTIME
      INTEGER               :: ISIZE, JSIZE
      INTEGER               :: L, T, NCRC
      CHARACTER(LEN=31 )    :: UNTS
      CHARACTER(LEN=255)    :: MSG, LOC
      REAL*4, POINTER       :: RDARR (:,:,:,:) => NULL()
      REAL*8, POINTER       :: TMPARR(:,:,:,:) => NULL()
      REAL*8, POINTER       :: ORIG_2D(:,:   ) => NULL()
      REAL*8, POINTER       :: REGR_2D(:,:   ) => NULL()

      REAL*8, ALLOCATABLE   :: XEDGE_IN (:)
      REAL*8, ALLOCATABLE   :: YSIN_IN  (:)
      REAL*8                :: XEDGE_OUT(HcoState%ISIZE+1) 
      REAL*8                :: YSIN_OUT (HcoState%JSIZE+1)

      !=================================================================
      ! HCOI_DATAREAD begins here
      !=================================================================

      ! Enter
      LOC = 'HCOI_DATAREAD (HCO_DATAREAD_MOD.F90)' 

      ! Copy horizontal grid dimensions from HEMCO state object
      ISIZE = HcoState%ISIZE
      JSIZE = HcoState%JSIZE

       ! Get time stamp to read
       ! Year
       IF ( ThisRdCont%ncYrs(1) == ThisRdCont%ncYrs(2) ) THEN
          YYYY = thisRdCont%ncYrs(1)
       ELSE
          YYYY = MIN( MAX(Clock%ThisYear,thisRdCont%ncYrs(1)), &
                      thisRdCont%ncYrs(2))
       ENDIF

       ! Month 
       IF ( ThisRdCont%ncMts(1) == ThisRdCont%ncMts(2) ) THEN
          MM = thisRdCont%ncMts(1)
       ELSE
          MM = MIN( MAX(Clock%ThisMonth,thisRdCont%ncMts(1)), &
                      thisRdCont%ncMts(2))
       ENDIF

       ! Day
       IF ( ThisRdCont%ncDys(1) == ThisRdCont%ncDys(2) ) THEN
          DD = thisRdCont%ncDys(1)
       ELSE
          DD = MIN( MAX(Clock%ThisDay,thisRdCont%ncDys(1)), &
                      thisRdCont%ncDys(2))
       ENDIF

       ! Hour
       IF ( ThisRdCont%ncHrs(1) == ThisRdCont%ncHrs(2) ) THEN
          HH = thisRdCont%ncHrs(1)
       ELSE
          HH = MIN( MAX(Clock%ThisHour,thisRdCont%ncHrs(1)), &
                      thisRdCont%ncHrs(2))
       ENDIF

      ! Get tracer ID
      TRCID = ThisRdCont%TrcID

      !-----------------------------------------------------------------
      ! Read data from netCDF 
      !-----------------------------------------------------------------

      ! verbose mode 
      if ( HcoState%Verbose ) then
         write(*,*) 'filling ',  TRIM(ThisRdCont%cName), ':'
         write(*,*) 'ncFile : ', TRIM(ThisRdCont%ncFile)
         write(*,*) 'ncPara : ', TRIM(ThisRdCont%ncPara)
         write(*,*) 'TrcID  : ', TrcID
         write(*,*) 'ScalID : ', ThisRdCont%ScalID
         write(*,*) 'Year   : ', YYYY
         write(*,*) 'Month  : ', MM
         write(*,*) 'Day    : ', DD
         write(*,*) 'Hour   : ', HH
      endif
 
      ! Note: use error syntax from ncdf_mod.F
      NCRC = 0 
!<<MSL>>      CALL NC_READ ( FILENAME = TRIM(ThisRdCont%ncFile), &
!<<MSL>>                     PARA     = TRIM(ThisRdCont%ncPara), &
!<<MSL>>                     YEAR     = YYYY,                    &
!<<MSL>>                     MONTH    = MM,                      &
!<<MSL>>                     DAY      = DD,                      &
!<<MSL>>                     HOUR     = HH,                      &
!<<MSL>>                     UNITS    = UNTS,                    & 
!<<MSL>>                     NLON     = NLON,                    &
!<<MSL>>                     NLAT     = NLAT,                    &
!<<MSL>>                     NLEV     = NLEV,                    &
!<<MSL>>                     NTIME    = NTIME,                   &
!<<MSL>>                     ARRAY    = RDARR,                   &
!<<MSL>>                     RC       = NCRC                      )

      IF ( NCRC /= 0 ) THEN
         CALL HCO_ERROR('NC_READ',LOC,RC); RETURN 
      ENDIF   

      ! Numbers of vertical levels must not exceed emission grid level
      IF ( NLEV > HcoState%LSIZE ) THEN
         MSG = 'Too many vert. levels in ' // TRIM(ThisRdCont%ncFile)
         CALL HCO_ERROR ( MSG, LOC, RC ); RETURN 
      ENDIF

      ! Cast to real*8
      ALLOCATE( TMPARR(NLON,NLAT,NLEV,NTIME) )
      TMPARR = RDARR

      ! Cleanup RDARR
      IF ( ASSOCIATED ( RDARR ) ) DEALLOCATE ( RDARR )

      !-----------------------------------------------------------------
      ! Convert to HEMCO units 
      !-----------------------------------------------------------------
      
      ! Convert to HEMCO units. This is kg/m2/s for fluxes and kg/m3 
      ! for concentrations. Ignore this if no species ID defined,
      ! i.e. for scale factors. 
      IF ( TRCID > 0 ) THEN

!<<MSL>>         CALL CHANGE_UNITS ( Array       = TMPARR,                      &
!<<MSL>>                             Units       = UNTS,                        &
!<<MSL>>                             MW_IN       = HcoState%SpecMW(TrcID),      & 
!<<MSL>>                             MW_OUT      = HcoState%EmSpecMW(TrcID),    & 
!<<MSL>>                             MOLEC_RATIO = HcoState%MolecRatio(TrcID),  & 
!<<MSL>>                             NcFile      = ThisRdCont%ncFile,           &
!<<MSL>>                             YYYY        = YYYY,                        &
!<<MSL>>                             MM          = MM,                          &
!<<MSL>>                             RC          = NCRC                          )

         IF ( NCRC /= 0 ) THEN
             CALL HCO_ERROR('CHANGE_UNITS',LOC,RC); RETURN 
         ENDIF
     
      ! For scale factors and masks, check if units are indeed unitless. 
      ELSE
         IF (  TRIM(UNTS) /= 'unitless' .AND. TRIM(UNTS) /= 'fraction' &
         .AND. TRIM(UNTS) /= '1' ) THEN
            MSG = 'This scale factor does not appear to be unitless: ' & 
                  // TRIM(ThisRdCont%ncFile)
            CALL HCO_WARNING( MSG, LOC, RC )
         ENDIF
      ENDIF

      !-----------------------------------------------------------------
      ! Regrid onto emissions grid 
      !-----------------------------------------------------------------
      
      ! Get input grid edges from netCDF file 
      ! Note: use error syntax from ncdf_mod!
      ALLOCATE( XEDGE_IN(NLON+1) )
      ALLOCATE( YSIN_IN (NLAT+1) )
      XEDGE_IN(:) = 0d0
      YSIN_IN (:) = 0d0
!<<MSL>>      CALL NC_READ_GRID( NLON,     NLAT,    ThisRdCont%ncFile, &
!<<MSL>>                         XEDGE_IN, YSIN_IN, NCRC                )
      IF ( NCRC /= 0 ) THEN
         CALL HCO_ERROR('NC_READ_GRID',LOC,RC); RETURN 
      ENDIF

      ! Get output grid edges from HEMCO state
      XEDGE_OUT(:) = HcoState%XEDGE(:,1,1)
      YSIN_OUT (:) = HcoState%YSIN (1,:,1) 
  
      ! Allocate output array if not yet defined
      IF ( .NOT. ASSOCIATED( ThisRdCont%Array ) ) THEN
         ALLOCATE ( ThisRdCont%Array(ISIZE,JSIZE,NLEV,NTIME) ) 

      ! Check dimension otherwise
      ELSE
         IF ( ( SIZE(ThisRdCont%Array,1) /= ISIZE ) .OR. &
              ( SIZE(ThisRdCont%Array,2) /= JSIZE ) .OR. &
              ( SIZE(ThisRdCont%Array,3) /= NLEV  ) .OR. &
              ( SIZE(ThisRdCont%Array,4) /= NTIME )       ) THEN

            ! Return w/ error
            MSG = 'Wrong dim of ThisRdCont%Array!'
            CALL HCO_ERROR ( MSG, LOC, RC )
            RETURN
         ENDIF
      ENDIF 

      ! Do regridding
      DO T = 1, NTIME
      DO L = 1, NLEV 

         ! Point to 2D slices to be regridded
         ORIG_2D => TMPARR(:,:,L,T)
         REGR_2D => ThisRdCont%Array(:,:,L,T)

         ! Do the regridding
         CALL MAP_A2A( NLON,  NLAT,  XEDGE_IN,  YSIN_IN,  ORIG_2D, &
                       ISIZE, JSIZE, XEDGE_OUT, YSIN_OUT, REGR_2D, 0, 0 )

        ! Free pointer
        ORIG_2D => NULL()
        REGR_2D => NULL()

      ENDDO !L
      ENDDO !T

      !-----------------------------------------------------------------
      ! Cleanup and leave 
      !-----------------------------------------------------------------
      IF ( ASSOCIATED(TMPARR   ) ) DEALLOCATE ( TMPARR   )
      IF ( ALLOCATED (XEDGE_IN ) ) DEALLOCATE ( XEDGE_IN )
      IF ( ALLOCATED (YSIN_IN  ) ) DEALLOCATE ( YSIN_IN  )

      ! Return w/ success
      RC = HCO_SUCCESS 

      END SUBROUTINE HCOI_DATAREAD
!EOC
#endif

END MODULE HCOI_DATAREAD_MOD
