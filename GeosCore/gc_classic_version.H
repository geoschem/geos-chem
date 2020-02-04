      ! For version number
      CHARACTER(LEN=16) :: GC_CLASSIC_VERSION

      !-----------------------------------------------------------------------
      ! Update the GEOS-Chem "Classic" version number (X.Y.Z) here
      !
      ! X = updated for each major version
      !
      ! Y = updated for each new 1-month benchmark
      !     This is done for each new set of features to be validated,
      !     or when changes are made to any run directory.
      !
      ! Z = updated for patches that do not break backwards compatibility
      !     with run directories from the last benchmarked version X.Y.0.
      !-----------------------------------------------------------------------
      GC_CLASSIC_VERSION = '12.8.0'
