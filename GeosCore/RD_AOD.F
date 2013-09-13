!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
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
      USE ERROR_MOD, ONLY : ERROR_STOP
      USE FILE_MOD,  ONLY : IOERROR

#     include "define.h"

      USE CMN_FJ_MOD, ONLY : JPMAX, JPPJ
      USE JV_CMN_MOD

      IMPLICIT NONE

!
! !INPUT PARAMETERS:
!
      INTEGER,          INTENT(IN) :: NJ1      ! Unit # of file to open
      CHARACTER(LEN=*), INTENT(IN) :: NAMFIL   ! Name of file to open
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
!  06 Aug 2010 - C. Carouge  - Add an error check when opening the file.
!  01 Aug 2012 - R. Yantosca - Now restore NJ1 to INTENT(IN) status
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES
!
      INTEGER :: I, J, K, N, NAA2, NWV, NR
      INTEGER :: IOS, NSPEC
      CHARACTER*30 :: SPECFIL(6)

      !================================================================
      ! RD_AOD begins here!
      !================================================================
      ! aerosol_mod.f expects aerosol in this order now
      DATA SPECFIL /"so4.dat","soot.dat","org.dat",
     &              "ssa.dat","ssc.dat","dust.dat"/
      !NWV=61 !number of wavelengths, maybe wants to be in common block
      !NR=7   !number of aerosol size bins in LUT
      !NSPEC = 6 !number of different species, in separate txt files
      !array of filenames
      DO k = 1, NSPAA
       !get name of next species datafile
       !NAMFIL=SPECFIL(k)

       ! open file
       OPEN( NJ1, FILE=TRIM( SPECFIL(k) ), STATUS='OLD', IOSTAT=IOS )

       ! Error check
       IF ( IOS /= 0 ) THEN
         IF ( am_I_Root ) THEN
          WRITE(6,100) trim(SPECFIL(k))
 100      FORMAT('Error opening filename=', a )
          CALL FLUSH(6)
         ENDIF
         CALL IOERROR( IOS, NJ1, 'RD_AOD:1')
       ENDIF

      
       ! Read header lines
       READ(  NJ1,'(A)'    ) TITLE0
       IF  ( am_I_Root ) WRITE( 6,  '(1X,A)' ) TITLE0
       READ(  NJ1,'(A)'    ) TITLE0
       !Second header line added for more info
       IF  ( am_I_Root ) WRITE( 6,  '(1X,A)' ) TITLE0
       READ(  NJ1,'(A)'    ) TITLE0

 110     FORMAT( 3x, a20 ) 
       DO i = 1, NRAA
        DO j = 1, NWVAA
         
         READ(NJ1,*) WVAA(j,k),RHAA(i,k),NRLAA(j,i,k),NCMAA(j,i,k),
     &               RDAA(i,k),RWAA(i,k),SGAA(i,k),QQAA(j,i,k),
     &               ALPHAA(j,i,k),REAA(i,k),SSAA(j,i,k),
     &               ASYMAA(j,i,k),(PHAA(j,i,k,n),n=1,8)
        ENDDO
       ENDDO
       ! Echo info to stdout
       IF ( am_I_Root ) THEN 
!        WRITE( 6, '(a)' ) 'Aerosol for AOD calculations'
!        WRITE( 6, * ) 'WV(nm)=',WVAA(2,k),
!     &  ' Q=',(QQAA(2,1,k)),'P1=',PHAA(2,i,k,2)
       ENDIF
       ! Close file
       CLOSE( NJ1 )
      ENDDO
      WRITE( 6, * ) 'Optics read for all wavelengths successfully'

      !now calculate the required wavelengths in the LUT to calculate
      !the requested AOD
      CALL CALC_AOD()

      END SUBROUTINE RD_AOD
!EOC

!BOP
!
! !IROUTINE: calc_aod
!
! !DESCRIPTION: Subroutine CALC\_AOD works out the closest tie points
! in the optics LUT wavelengths and the coefficients required to
! calculate the angstrom exponent for interpolating to the requested
! wavelength. If the wavelength requested matches a standard wavelength
! in the LUT then we skip the interpolation
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE CALC_AOD()

#     include "define.h"


      USE JV_CMN_MOD, ONLY : NWVAA, WVAA
      USE JV_CMN_MOD, ONLY : WVSELECT, IWVSELECT, NWVSELECT
      USE JV_CMN_MOD, ONLY : ACOEF_WV, BCOEF_WV, CCOEF_WV
      USE JV_CMN_MOD, ONLY : NWVREQUIRED, IWVREQUIRED
      USE LOGICAL_MOD,ONLY : LRAD

      IMPLICIT NONE

      INTEGER             :: MINWV, MAXWV, N, N0, N1, W, NSTEP
      REAL*8              :: WVDIF

      !currently only allow one wavelength at once until all diagnostics
      !updated for multiple wavelength outputs
      !NWVSELECT = 2 !now defined by user input

      !look for closest WV in LUT to the requested WV
      !If LRAD set then only looking at RT wavelengths so that we can output
      !combined optics direct from RT
      IF (LRAD) THEN
       N0=12
       N1=NWVAA
       NSTEP=1
      ELSE
       !cycle over standard wavelengths
       N0=1
       N1=11
       NSTEP=1
      ENDIF
      NWVREQUIRED=0
      DO W=1,NWVSELECT
       MINWV     = -999
       MAXWV     =  999
       DO N=N0,N1,NSTEP
        WVDIF = WVAA(N,1)-WVSELECT(W)
        IF ((WVDIF.LE.0).AND.(WVDIF.GT.MINWV)) THEN
         MINWV = WVDIF
         IWVSELECT(1,W)=N
        ENDIF
        IF ((WVDIF.GE.0).AND.(WVDIF.LT.MAXWV)) THEN
         MAXWV = WVDIF
         IWVSELECT(2,W)=N
        ENDIF
       ENDDO
        IF (IWVSELECT(2,W).EQ.IWVSELECT(1,W)) THEN
         !we have a match!
         MINWV=0
         MAXWV=0
         !add this wavelength to those for output
         NWVREQUIRED=NWVREQUIRED+1
         !IWVSELECT(1,W)=N
         !IWVSELECT(2,W)=N
         IWVREQUIRED(NWVREQUIRED)=IWVSELECT(1,W)
        ELSE
         !we are going to have to interpolate to the requested wavelength
         NWVREQUIRED=NWVREQUIRED+1
         IWVREQUIRED(NWVREQUIRED)=IWVSELECT(1,W)
         NWVREQUIRED=NWVREQUIRED+1
         IWVREQUIRED(NWVREQUIRED)=IWVSELECT(2,W)
        ENDIF
       !ENDDO
      

      !Error check - ensure we have a match or requested wavelength
      !falls within two LUT tie points
      IF (MINWV.EQ.-999) THEN
       ! requested wavelength is shorter than min wv in LUT
       ! set to min
       write(6,*) 'ERROR requested wavelength is too short!!'
       write(6,*) 'Defaulting to LUT min: ',WVAA(NWVAA,1)
       IWVSELECT(1,W)=NWVAA
       IWVSELECT(2,W)=NWVAA !229.8nm
       NWVREQUIRED=NWVREQUIRED-1
       IWVREQUIRED(NWVREQUIRED)=IWVSELECT(1,W)
      ENDIF
      IF (MAXWV.EQ.999) THEN
       ! requested wavelength is longer than min wv in LUT
       ! set to max
       write(6,*) 'ERROR requested wavelength is too long!!'
       write(6,*) 'Defaulting to LUT min: ',WVAA(8,1)
       IWVSELECT(1,W)=8
       IWVSELECT(2,W)=8 !
       NWVREQUIRED=NWVREQUIRED-1
       IWVREQUIRED(NWVREQUIRED)=IWVSELECT(1,W)
      ENDIF
      
      !now calcualte the angstrom exponent coefs for interpolation -
      !this is done here to save time in aerosol_mod.F
      !DO I=1,NWVS
      IF (IWVSELECT(1,W).NE.IWVSELECT(2,W)) THEN
       ACOEF_WV(W) = WVAA(IWVSELECT(2,W),1)/WVSELECT(W)
       BCOEF_WV(W) =1.0d0/(LOG(WVAA(IWVSELECT(2,W),1)/
     &                         WVAA(IWVSELECT(1,W),1)))
       !relative location of selected wavelength between tie points
       !for interpolating SSA and ASYM for output in aerosol_mod.F and
       !dust_mod.F
       CCOEF_WV(W) =(WVSELECT(W)-WVAA(IWVSELECT(1,W),1))/
     &              (WVAA(IWVSELECT(2,W),1)-WVAA(IWVSELECT(1,W),1))
      ENDIF
      write(6,*) 'WAVELENGTHS: ',NWVSELECT
      write(6,*) 'WAVELENGTH REQUESTED:',WVSELECT(W)
      write(6,*) 'WAVELENGTH REQUIRED:', NWVREQUIRED
      write(6,*) IWVSELECT(1,W),WVAA(IWVSELECT(1,W),1)
      write(6,*) IWVSELECT(2,W),WVAA(IWVSELECT(2,W),1)
      write(6,*) ACOEF_WV(W),BCOEF_WV(W),CCOEF_WV(W)
      write(6,*) '*********************************'
      ENDDO !NWVSELECT
      END SUBROUTINE CALC_AOD
