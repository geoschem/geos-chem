! $Id: lightdist.f,v 1.1 2003/06/30 20:26:05 bmy Exp $
      SUBROUTINE LIGHTDIST( I, J, LTOP, H0, XLAT, TOTAL )
!
!******************************************************************************
!  Subroutine LIGHTDIST reads in the CDF used to partition the column lightning
!  NOx into the GEOS-CHEM vertical layers. (yhw, 1997; mje, bmy, 9/18/02)
!
!  Arguments as Input:
!  ============================================================================
!  (1-2) I, J    (INTEGER) : (Lon,Lat) indices of the surface grid box
!  (3  ) LTOP    (INTEGER) : Level where convective cloud top is found
!  (4  ) H0      (REAL*8 ) : Convective cloud top height [m]
!  (5  ) XLAT    (REAL*8 ) : Latitude of surface grid box (I,J) [degrees]
!  (6  ) TOTAL   (REAL*8 ) : Total # of NOx molec. released from lightning
!
!  Variables from "CMN_NOX":
!  ============================================================================
!  (1  ) SLBASE  (REAL*8 ) : Array containing lightning NOx emissions
!
!  References:
!  ============================================================================
!  (1 ) Pickering et al., JGR 103, 31,203 - 31,316, 1998.
!
!  NOTES:
!  (1 ) Use functions IS_LAND and IS_WATER to determine if the given grid
!        box is over land or water.  These functions work for all DAO met
!        field data sets. (bmy, 4/2/02)
!  (2 ) Renamed M2 to LTOP and THEIGHT to H0 for consistency w/ variable names
!        w/in "lightning.f".  Now read the "light_dist.dat.geos3" file for 
!        GEOS-3 directly from the DATA_DIR/lightning_NOx_200203/ subdirectory.
!        Now read the "light_dist.dat" file for GEOS-1, GEOS-STRAT directly 
!        from the DATA_DIR/lightning_NOx_200203/ subdirectory.  Added 
!        descriptive comment header.  Now trap I/O errors across all 
!        platforms with subroutine "ioerror.f".  Updated comments, cosmetic 
!        changes.  Redimension FRAC(NNLIGHT) to FRAC(LLPAR). (bmy, 4/2/02)
!  (3 ) Deleted obsolete code from April 2002.  Now reference IU_FILE and
!        IOERROR from "file_mod.f".  Now use IU_FILE instead of IUNIT as the
!        file unit number. (bmy, 6/27/02)
!  (4 ) Now reference BXHEIGHT from "dao_mod.f" (bmy, 9/18/02)
!******************************************************************************
!
      ! References to F90 modules
      USE DAO_MOD,  ONLY : IS_LAND, IS_WATER, BXHEIGHT
      USE FILE_MOD, ONLY : IU_FILE, IOERROR

      IMPLICIT NONE

#     include "CMN_SIZE"   ! Size parameters
#     include "CMN_NOX"    ! SLBASE
#     include "CMN_SETUP"  ! DATA_DIR

      ! Arguments
      INTEGER, INTENT(IN) :: I,  J,    LTOP
      REAL*8,  INTENT(IN) :: H0, XLAT, TOTAL

      ! Local variables
      INTEGER             :: M, MTYPE, L, III, IOS, IUNIT, JJJ
      LOGICAL, SAVE       :: FIRST = .TRUE.
      REAL*8              :: ZHEIGHT

      ! Types of lightning to consider
      INTEGER, PARAMETER  :: NLTYPE = 3

#if   defined( GEOS_3 )

      ! Number of data points in the PDF -- GEOS-3 only
      ! Should work for other fields, however I seem to have 
      ! calculated something different  So wait and see. (mje, 6/19/01)
      INTEGER, PARAMETER  :: NNLIGHT = 3200
#else

      ! Number of data points in the PDF -- GEOS-1, GEOS-STRAT
      INTEGER, PARAMETER  :: NNLIGHT = 100

#endif

      ! Array to hold CDF's read from disk
      REAL*8, SAVE        :: PROFILE(NNLIGHT,NLTYPE)

      ! Array to hold fraction of NOx going into each GEOS-CHEM level
      REAL*8              :: FRAC(LLPAR)

      ! File name for lightning distribution file (bmy, 4/2/02)
      CHARACTER(LEN=255)  :: FILENAME

      !=================================================================
      ! LIGHTDIST begins here!
      !=================================================================

      ! Do only on the first call
      IF ( FIRST ) THEN
         
         ! Reset first-time flag
         FIRST = .FALSE.

#if   defined( GEOS_3 ) 
         
         !==============================================================
         ! Read data for GEOS-3
         ! 
         ! NOTE: Since the GEOS-3 vertical grid has a much finer
         ! resolution near the surface than do the GEOS-1 or GEOS-STRAT
         ! grids, we had to interpolate the Pickering CDF's to a much
         ! finer mesh, with 3200 points instead of 100.  This was done
         ! by Mat Evans (mje@io.harvard.edu).  The vertical resolution
         ! of the CDF's in the file read in below is 0.05 km. 
         !==============================================================

         ! Define filename for GEOS-3 CDF file
         FILENAME = TRIM( DATA_DIR ) // 
     &              'lightning_NOx_200203/light_dist.dat.geos3'        
      
         ! Open file containing lightning PDF data
         OPEN( IU_FILE, FILE=TRIM(FILENAME), STATUS='OLD', IOSTAT=IOS )
         IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_FILE, 'lightdist:1' )
         
         ! For GEOS-3 only: read 12 header lines
         DO III = 1, 12
            READ( IU_FILE, '(a)', IOSTAT=IOS ) 
            IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_FILE, 'lightdist:2' )
         ENDDO
         
         ! For GEOS-3: read NNLIGHT types of lightning profiles
         DO III = 1, NNLIGHT
            READ( IU_FILE,*,IOSTAT=IOS) (PROFILE(III,JJJ),JJJ=1,NLTYPE)
         ENDDO
         
         ! Close file
         CLOSE( IU_FILE )

#else

         !==============================================================
         ! Read data for GEOS-1, GEOS-STRAT
         !
         ! Here we read in the original Pickering CDF's for NNLIGHT=3
         ! different types of lightning: tropical marine, tropical 
         ! continental, and midlatitude continental.  The vertical 
         ! resolution of the CDF's in the file read in below is 0.16 km. 
         !==============================================================

         ! Define filename for GEOS-1, GEOS-STRAT PDF file
         FILENAME = TRIM( DATA_DIR ) // 
     &              'lightning_NOx_200203/light_dist.dat'        

         ! Open file containing lightning CDF data
         OPEN( IU_FILE, FILE=TRIM(FILENAME), STATUS='OLD', IOSTAT=IOS )
         IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_FILE, 'lightdist:3' )

         ! For GEOS-1, GEOS-STRAT: Read 1 header line
         READ( IU_FILE, *, IOSTAT=IOS )
         IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_FILE, 'lightdist:4' )

         ! For GEOS-1, GEOS-STRAT: Read data
         DO III = 1, NNLIGHT
            READ( IU_FILE,*,IOSTAT=IOS ) (PROFILE(III,JJJ),JJJ=1,NLTYPE)
            IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_FILE, 'lightdist:5' )
         ENDDO

         ! Close file
         CLOSE( IU_FILE )

#endif

      ENDIF

      !=================================================================
      ! Use functions IS_LAND and IS_WATER from "dao_mod.f" to test 
      ! whether surface location (I,J) is over land or over water.
      !
      ! Depending on the combination of land/water and latitude, assign
      ! a flag describing the type of lightning:
      !
      !   MTYPE = 1: ocean lightning
      !   MTYPE = 2: tropical continental lightning (from 30S to 30N)
      !   MTYPE = 3: midlatitude continental lightning (all other lats)
      !=================================================================
      IF ( IS_LAND(I,J) ) THEN
         IF ( ABS( XLAT ) <= 30 ) THEN
            MTYPE = 2
         ELSE
            MTYPE = 3
         ENDIF
      ELSE IF ( IS_WATER(I,J) ) THEN 
         MTYPE = 1
      ENDIF

      ! Add fix for GEOS-3 fields (mje, bmy, 6/19/01)
#if   defined( GEOS_3 )

      !=================================================================
      ! For GEOS-3 only:
      !
      ! Use the CDF for this type of lightning to partition the total
      ! column lightning into the GEOS-3 vertical layers.
      !=================================================================
      ZHEIGHT = 0.0

      ! Compute the height [km] at the top of each vertical level.
      ! Look up the cumulative fraction of NOx for each vertical level
      DO L = 1, LTOP
         ZHEIGHT = ZHEIGHT + BXHEIGHT(I,J,L)
         FRAC(L) = PROFILE( NINT( ( ZHEIGHT/H0 )*3200. ), MTYPE ) *0.01
      ENDDO

      ! Convert from cumulative fraction to fraction for each level
      DO L = LTOP, 2, - 1
         FRAC(L) = FRAC(L) - FRAC(L-1)
      ENDDO 
      
      ! Partition lightning NOx into the SLBASE array
      DO L = 1, LTOP
         SLBASE(I,J,L) = SLBASE(I,J,L) + ( FRAC(L) * TOTAL )
      ENDDO

#else

      !=================================================================
      ! For GEOS-1 and GEOS-STRAT only:
      !
      ! Use the PDF for this type of lightning to partition the total
      ! column lightning into the GEOS-3 vertical layers
      !=================================================================

      ZHEIGHT = 0.

      ! Compute the height [km] at the top of each vertical level.
      ! Look up the cumulative fraction of NOx for each vertical level
      DO L = 1, LTOP-1
         ZHEIGHT = ZHEIGHT + BXHEIGHT(I,J,L)
         FRAC(L) = PROFILE( NINT( ( ZHEIGHT / H0 ) * 100. ), MTYPE )
      ENDDO

      ! If there is any lightning NOx yet to be partitioned out,
      ! then place that in the level where the cloud top occurs.
      FRAC(LTOP) = 1.
         
      ! Convert from cumulative fraction to fraction for each level
      DO L = LTOP, 2, -1
         FRAC(L) = FRAC(L) - FRAC(L-1)
      ENDDO

      ! Partition lightning NOx into the SLBASE array
      DO L = 1, LTOP
         SLBASE(I,J,L) = SLBASE(I,J,L) + ( FRAC(L) * TOTAL )
      ENDDO
#endif

      ! Return to calling program
      END SUBROUTINE LIGHTDIST





