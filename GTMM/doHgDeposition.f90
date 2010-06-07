      SUBROUTINE doHgDeposition(LCPLE)

      USE defineConstants
      USE loadCASAinput
      USE defineArrays

      implicit none

      ! Arguments
      LOGICAL :: LCPLE

      REAL*8     :: photo_frac(n_veg, 1)
      INTEGER   ::  i
      
      LOGICAL, SAVE :: FIRST = .TRUE.

!--- Previous to (ccc, 11/3/09)
!      IF (yr .eq. (NPPequilibriumYear+1) .and. mo .eq. 1) THEN
      IF ( FIRST .AND. .NOT. LCPLE ) THEN
              Hg0_surf_soil(:,1)=0.0d0
              HgII_surf_soil(:,1)=0.0d0
              hleafpool_Hg(:,1)=0.0d0
              leafpool_Hg(:,1)=0.0d0

              FIRST = .FALSE.
      ENDIF
      
      !!!!!   Hg0dry
      !! 1 - deposited to leaf and soil surfaces
      !! 2 - directly incorporated into leaf tissue via stomates
      !! 3 - volatilized as a fxn of temperature

!$OMP PARALLEL       &
!$OMP DEFAULT(SHARED)
!$OMP WORKSHARE
      fstom(:,1)=(LAI(:,mo)/5)  ! chosen to match Rea et al seasonal cycle
!--- Previous to (ccc, 11/9/09)
!      DO i=1, n_veg
!         IF (fstom(i,1) .gt. 1.0d0) THEN
!            fstom(i,1)=1.0d0
!         ENDIF
!      END DO
      WHERE (fstom(:,1) > 1.0d0) 
         fstom(:,1)=1.0d0
      END WHERE
      fstom(:,1)=fstom(:,1)
      fsoil(:,1)=1.00d0-fstom(:,1)
      Hg0_surf_soil(:,1)=Hg0_surf_soil(:,1)+(Hg0dry(:,1)*fsoil(:,1))
      leafpool_hg(:,1)=leafpool_hg(:,1)+(Hg0dry(:,1)*fstom(:,1))
      hleafpool_hg(:,1)=hleafpool_hg(:,1)+(Hg0dry(:,1)*fstom(:,1))
      !all elemental hg sitting in surface pools is volatilized each
      !month
      freemitted(:,1)=1.0d0
      reemitted(:,1)=(freemitted(:,1)*Hg0_surf_soil(:,1))

      Hg0_surf_soil(:,1)=Hg0_surf_soil(:,1)-reemitted(:,1)
      !!! HgII dry
      !! 1 - deposited to leaf and soil surfaces
      !! 2 - directly incorporated into leaf tissue via stomates
      !! 3 - photoreduced as fxn of T and light
      !! 4 - washed off leaf and soil surfaces and added to HgIIwet

      ! for 1 same fleaf, fsoil and fstom as for Hg0dry
      HgII_surf_soil(:,1)=HgII_surf_soil(:,1)+(HgIIdry(:,1)*fsoil(:,1))
      leafpool_hg(:,1)=leafpool_hg(:,1)+(HgIIdry(:,1)*fstom(:,1))
      hleafpool_hg(:,1)=hleafpool_hg(:,1)+(HgIIdry(:,1)*fstom(:,1))

      photo_frac(:,1)=0.667577d0*(1.0d0-exp(solrad1(:,mo)*(-1d0)*(0.01603d0)))
      !photo_frac equation is fit to curve in Rolfhus and Fitzgerald 

!--- Previous to (ccc, 11/9/09)
!      DO i=1, n_veg
!         IF (photo_frac(i,1) .le. 0.0d0) THEN
!            photo_frac(i,1)=0.0d0
!         ENDIF
!      END DO                       
      WHERE (photo_frac(:,1) <= 0.0d0) 
         photo_frac(:,1)=0.0d0
      END WHERE
      photoreduced(:,1)=photo_frac(:,1)*(HgII_surf_soil(:,1))

      HgII_surf_soil(:,1)=HgII_surf_soil(:,1)-photoreduced(:,1)
     
      !photoreduce some wet deposition
      temp_hg(:,1)=HgIIwet(:,1)*photo_frac(:,1)
      photoreduced(:,1)=photoreduced(:,1)+temp_hg(:,1)
      HgIIwet(:,1)=HgIIwet(:,1)-temp_hg(:,1)
!--- Previous to (ccc, 11/9/09)
!      IF (minval(HgIIwet(:,1)) .lt. 0d0) THEN
!         DO i=1,n_veg
!            IF (HgIIwet(i,1) .lt. 0d0) THEN
!              HgIIwet(i,1)=0.0d0
!            ENDIF
!         END DO
!      ENDIF
      WHERE (HgIIwet(:,1) < 0d0)
         HgIIwet(:,1)=0.0d0
      END WHERE

      !if there is rain - wash off HgII and add to HgII wet pool  
!--- Previous to (ccc, 9/11/09)
!      DO i=1, n_veg
!         IF (ppt1(i,mo) .gt. 0d0 .and. airt1(i,mo) .gt. 0d0) THEN
!                 HgIIwet(i,1)=HgIIwet(i,1)+HgII_surf_soil(i,1)
!                 HgII_surf_soil(i,1)=0.0d0
!         ENDIF
!      END DO
      WHERE (ppt1(:,mo) > 0d0 .AND. airt1(:,mo) > 0d0)
         HgIIwet(:,1)=HgIIwet(:,1)+HgII_surf_soil(:,1)
         HgII_surf_soil(:,1)=0.0d0
      END WHERE
!$OMP END WORKSHARE
!$OMP END PARALLEL

      END SUBROUTINE doHgDeposition
