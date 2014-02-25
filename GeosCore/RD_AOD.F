!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: rd_aod
!
! !DESCRIPTION: Subroutine RD\_AOD reads aerosol phase functions that are
!  used to scale diagnostic output to an arbitrary wavelengh.  This
!  facilitates comparing with satellite observations.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE RD_AOD( NJ1, NAMFIL, am_I_Root )
!
! !USES:
!
      USE ERROR_MOD,  ONLY : ERROR_STOP
      USE FILE_MOD,   ONLY : IOERROR
      USE CMN_FJ_MOD, ONLY : JPMAX, JPPJ
      USE JV_CMN_MOD

      IMPLICIT NONE
!
! !INPUT PARAMETERS:
!     
      INTEGER,          INTENT(IN) :: NJ1         ! Unit # of file to open
      CHARACTER(LEN=*), INTENT(IN) :: NAMFIL      ! Name of file to open
      LOGICAL,          INTENT(IN) :: am_I_Root   ! Is this the root CPU?
!
! !REMARKS:
!  The jv_spec_aod.dat file contains the optical properties for aerosols 
!  at a single wavelength to be used in the online calculation of the aerosol 
!  optical depth diagnostics.  The default properties are provided at 550 nm. 
!  These properties have been calculated using the same size and optical 
!  properties as the jv_spec.dat file used for the FAST-J photolysis 
!  calculations.  The user can exchange this set of properties with those at 
!  another wavelength.  We recommend that the wavelength used be included 
!  in the first line of the header for traceability (this line is output to 
!  the GEOS-Chem log file during run time). A complete set of optical 
!  properties from 250-2000 nm for aerosols is available at:  
!  ftp://ftp.as.harvard.edu/geos-chem/data/aerosol_optics/hi_spectral_res
!                                                                             .
!     -- Colette L. Heald, 05/10/10)
!
!  Important variables:
!                                                                             .
!     NAMFIL       Name of spectral data file (jv_spec_aod.dat)
!     NJ1          Channel number for reading data file
!     NAA2         Number of categories for scattering phase functions
!     QAA_AOD      Aerosol scattering phase functions
!     WAA_AOD      Wavelengths for the NK supplied phase functions
!     PAA_AOD      Phase function: first 8 terms of expansion
!     RAA_AOD      Effective radius associated with aerosol type
!     SSA_AOD      Single scattering albedo
!
! !REVISION HISTORY: 
!  10 May 2010 - C. Heald    - Initial version
!  06 Aug 2010 - C. Carouge  - Add an error check when opening the file
!  01 Aug 2012 - R. Yantosca - Now restore NJ1 to INTENT(IN) status
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES
!
      INTEGER :: I, J, K, NAA2
      INTEGER :: IOS

      !================================================================
      ! RD_AOD begins here!
      !================================================================

      ! open file
      OPEN( NJ1, FILE=TRIM( NAMFIL ), STATUS='OLD', IOSTAT=IOS )

      ! Error check
      IF ( IOS /= 0 ) THEN
         IF ( am_I_Root ) THEN
            WRITE(6,100) trim(NAMFIL)
 100        FORMAT('Error opening filename=', a )
            CALL FLUSH(6)
         ENDIF
         CALL IOERROR( IOS, NJ1, 'RD_AOD:1')
      ENDIF

      
      ! Read header lines
      READ( NJ1,'(A)' ) TITLE0
      IF  ( am_I_Root ) WRITE( 6, '(1X,A)' ) TITLE0
      READ( NJ1,'(A)' ) TITLE0

      ! Read aerosol phase functions (one wavelength only):
      READ( NJ1,'(A10,I5,/)' ) TITLE0,NAA2
      DO j = 15, NAA
         READ(NJ1,110) TITLEA(j)
 110     FORMAT( 3x, a20 )
         IF ( am_I_Root ) WRITE(6,*) TITLEA(j)
         READ(NJ1,*) WAA_AOD(j),QAA_AOD(j),RAA_AOD(j),SSA_AOD(j),
     &               (PAA_AOD(i,j),i=1,8)
      ENDDO

      ! Echo info to stdout
      IF ( am_I_Root ) THEN
         WRITE( 6, '(a)' ) 'Aerosol Qext for AOD calculations'
         DO J=15,NAA
            WRITE( 6, * ) TITLEA(J),J,'  Qext =',(QAA_AOD(J))
         ENDDO
      ENDIF

      ! Close file
      CLOSE( NJ1 )

      END SUBROUTINE RD_AOD
!EOC
