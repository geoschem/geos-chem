#!/usr/bin/perl -w
#------------------------------------------------------------------------------
#                  GEOS-Chem Global Chemical Transport Model                  !
#------------------------------------------------------------------------------
#BOP
#
# !MODULE: finn.pl
#
# !DESCRIPTION: Reads the FINN_emission_factors.txt file and converts
#  the data in the file to hardwired F90 commands for inlining into
#  an include file.  This facilitates I/O in the ESMF environment.
#\\
#\\
# !INTERFACE:
#
# !USES:
#
require 5.003;      # need this version of Perl or newer
use English;        # Use English language
use Carp;           # Strict error checking
use strict;         # IMPLICIT NONE syntax
#
# !PUBLIC MEMBER FUNCTIONS:
#  &baseName        : Returns the last part of a full path (e.g. notdir)
#  &num2str         ; Converts a 2-digit number to a string (w/ pad space)
#  &getProTeXHeader : Returns the ProTeX subroutine header text
#  &main            : Driver program
#
# !REMARKS:
#  You need to specify the root emissions data directory by setting
#  the $HEMCO_DATA_ROOT environment variable accordingly.
#
# !CALLING SEQUENCE:
#  gfed3.pl
#
# !REVISION HISTORY:
#  27 Jun 2014 - R. Yantosca - Initial version
#EOC
#------------------------------------------------------------------------------
#                  GEOS-Chem Global Chemical Transport Model                  !
#------------------------------------------------------------------------------
#BOP
#
# !IROUTINE: baseName
#
# !DESCRIPTION: Returns the last part of a full path name.  Similar to the
#  GNU Make "notdir" function.
#\\
#\\
# !INTERFACE:
#
sub baseName($) {
#
# !INPUT PARAMETERS:
#
  my ( $dir )  =  @_;   # String to be parsed
#
# !RETURN VALUE:
#
  my $baseName = "";    # Directory name minus the full path
#
# !CALLING SEQUENCE:
#  &baseName( $dir );
#
# !REVISION HISTORY:
#  30 Jul 2013 - R. Yantosca - Initial version
#EOP
#------------------------------------------------------------------------------
#BOC

  # Take the text following the colon
  my @result = split( '/', $dir );

  # Return the last part of the directory
  $baseName = $result[ scalar( @result ) - 1 ];

  # Return to calling routine
  return( $baseName );
}
#EOP
#------------------------------------------------------------------------------
#                  GEOS-Chem Global Chemical Transport Model                  !
#------------------------------------------------------------------------------
#BOP
#
# !IROUTINE: num2str
#
# !DESCRIPTION: Converts a 2-digit number to a string.  Used for padding
#  columns in the F90 output so that they line up.
#\\
#\\
# !INTERFACE:
#
sub num2str($) {
#
# !INPUT PARAMETERS:
#
  # Input number
  my( $num ) = @_;
#
# !REVISION HISTORY:
#  11 Aug 2014 - R. Yantosca - Initial version
#EOP
#------------------------------------------------------------------------------
#BOC
#
# !LOCAL VARIABLES:
#
  my $str = "$num";

  # Pad spaces (cheap algorithm)
  if ( $str < 10 ) { $str = " $str"; }

  # Return numeric string
  return( $str );
}
#------------------------------------------------------------------------------
#                  GEOS-Chem Global Chemical Transport Model                  !
#------------------------------------------------------------------------------
#BOP
#
# !IROUTINE: getProTeXHeader
#
# !DESCRIPTION: Returns a subroutine header for the include file.
#\\
#\\
# !INTERFACE:
#
sub getProTeXHeader() {
#
# !REVISION HISTORY:
#  11 Aug 2014 - R. Yantosca - Initial version
#EOP
#------------------------------------------------------------------------------
#BOC
#
# !LOCAL VARIABLES:
#
  # HERE Docs
  my @header  = <<EOF;
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: hcox_finn_include.H
!
! !DESCRIPTION: Include file with FINN emission factor data that was
!  originally contained in files FINN\\_EFratios\\_CO2.csv and
!  FINN\\_VOC\\_speciation.csv.  We transform these data into hardwired F90
!  assignment statements so that we can avoid reading ASCII files in the
!  ESMF environment.
!
! !REMARKS:
!  ABOUT THIS FILE:
!  ----------------
!  This file was created by script HEMCO/Extensions/Preprocess/finn.pl.
!  This script can be executed with the following command:
!
!    cd HEMCO/Extensions/Preprocess
!    make finn
!
!  This will regenerate this include file from the original data and
!  automatically place it in the HEMCO/Extensions directory.
!
!  White space has been removed in order to reduce the file size as much
!  as possible.  If you have to recreate this file, then it is easier to
!  generate via the Perl script than to try to hand edit the code below.
!
! !REVISION HISTORY:
!  11 Aug 2014 - R. Yantosca - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
EOF

  # Return to calling program
  return( @header );
}
#EOP
#------------------------------------------------------------------------------
#                  GEOS-Chem Global Chemical Transport Model                  !
#------------------------------------------------------------------------------
#BOP
#
# !MODULE: makeFinnInclude
#
# !DESCRIPTION: Reads the FINN emission factors and VOC speciation files
#  and writes out equivalent F90 assignment statements to an include file
#\\
#\\
# !INTERFACE:
#
sub makeFinnInclude($$$) {
#
# !INPUT PARAMETERS:
#
  # $co2Path : "FINN_EFratios.csv" file path
  # $vocPath : "FINN_VOC_speciation.csv" file path
  # $outFile : File path of the FINN include file
  my( $co2Path, $vocPath, $outFile ) = @_;

# !REVISION HISTORY:
#  17 Apr 2014 - R. Yantosca - Initial version
#EOP
#------------------------------------------------------------------------------
#BOC
#
# !LOCAL VARIABLES:
#
  # Strings
  my $line    = "";
  my $subStr  = "";
  my $name    = "";
  my $spcStr  = "";
  my $co2File = &baseName( $co2Path );
  my $vocFile = &baseName( $vocPath );

  # Numbers
  my $i       = 0;
  my $fac     = 0;
  my $col     = 0;
  my $nCo2    = 0;
  my $nVoc    = 0;

  # Arrays
  my @result  = ();
  my @header  = &getProTeXHeader();

  #========================================================================
  # Open files and write the subroutine header to the output file
  #========================================================================

  # Open input and output files
  open( C, "<$co2Path" ) or croak( "Cannot open $co2Path!\n" );
  open( V, "<$vocPath" ) or croak( "Cannot open $vocPath!\n" );
  open( O, ">$outFile" ) or croak( "Cannot open $outFile!\n" );

  # Write subroutine header to output file
  foreach $line ( @header ) { print O "$line"; }

  #=========================================================================
  # Read data from the CO2 speciation file
  #=========================================================================

  # Parse input file and write F90 commands to output file
  foreach $line ( <C> ) {

    # Remove newlines
    chomp( $line );
    $line =~ s/\r//g;

    # Parse the line
    if    ( $line =~ m/FINN Emission/ ) { ; } # Skip header line
    elsif ( $line =~ m/Prepared by/   ) { ; } # Skip header line
    elsif ( $line =~ m/GenVegCode/    ) {

      #=====================================================================
      # This line has the CO2 file species names
      #=====================================================================

      # Split the line by commas
      @result = split( ',', $line );

      # Number of emission factors in file (skip 1st 2 columns)
      $nCo2   = scalar( @result ) - 2;

      # Print comment lines
      print O "!---------------------------------------------------\n";
      print O "! Data from $co2File\n";
      print O "!---------------------------------------------------\n\n";
      print O "! Number of emissions species in $co2File\n";
      print O "N_SPEC_EMFAC=$nCo2\n";
      print O "N_SPECSTRS=$nCo2\n\n";
      print O "! Allocate the EMFAC_IN temporary array\n";
      print O "ALLOCATE ( EMFAC_IN(N_SPEC_EMFAC, N_EMFAC), STAT=AS )\n";
      print O "IF ( AS/=0 ) THEN\n";
      print O " CALL HCO_ERROR('Cannot allocate EMFAC_IN',RC)\n";
      print O " RETURN\n";
      print O "ENDIF\n";
      print O "EMFAC_IN=0.0_dp\n\n";
      print O "! Species from $co2File\n";

      # Loop over substrings, starting from the 2nd column
      for ( $i = 2; $i < scalar( @result ); $i++ ) {

	# Species name
	$name   = "$result[$i]";

	# Species index (F90 notation, starting w/ 1)
	$spcStr = &num2str( $i - 1 );

	# Write F90 code to output file
	print O "IN_SPEC_NAME($spcStr)=\"$name\"\n";
      }

      print O "\n";

    } else {

      #=====================================================================
      # The following lines have the CO2 file emission factors
      #=====================================================================

      # Split the line by commas
      @result = split( ',', $line );

      # Emission factor number
      # CROPS is listed as #9 but reset to #6
      $fac = $result[0];
      if ( $fac > 6 ) { $fac = 6 }

      # Print emission factor name
      print O "! $result[1] emission factors from $co2File\n";

      # Loop over substrings
      for ( $i = 0; $i < $nCo2; $i++ ) {

	# Species index (F90 notation, starting w/ 1)
	$spcStr = &num2str( $i + 1 );

	# Column w/ data for the given species
	$col    = $i + 2;

	# Write the F90 commands to the output file
	print O "EMFAC_IN($spcStr,$fac)=$result[$col]_dp\n";
      }

      print O "\n";

    }
  }

  #=========================================================================
  # Read from the VOC speciation file
  #=========================================================================

  # Parse input file and write F90 commands to output file
  foreach $line ( <V> ) {

    # Remove newlines
    chomp( $line );
    $line =~ s/\r//g;

    # Parse the line
    if    ( $line =~ m/FINN VOC/    ) { ; } # Skip header line
    elsif ( $line =~ m/Prepared by/ ) { ; } # Skip header line
    elsif ( $line =~ m/GenVegCode/  ) {

      #=====================================================================
      # This line has the VOC file species names
      #=====================================================================

      # Split the line by commas
      @result = split( ',', $line );

      # Number of NMOC ratios in file (skip 1st 2 columns)
      $nVoc   = scalar( @result ) - 2;

      # Print comment lines
      print O "!---------------------------------------------------\n";
      print O "! Data from $vocFile\n";
      print O "!---------------------------------------------------\n\n";
      print O "! Number of emissions species in $vocFile\n";
      print O "N_NMOC=$nVoc\n";
      print O "N_NMOCSTRS=$nVoc\n\n";
      print O "! Allocate the NMOC_RATIO_IN temporary array\n";
      print O "ALLOCATE ( NMOC_RATIO_IN(N_NMOC, N_EMFAC), STAT=AS )\n";
      print O "IF ( AS/=0 ) THEN\n";
      print O " CALL HCO_ERROR('Cannot allocate NMOC_RATIO_IN',RC)\n";
      print O " RETURN\n";
      print O "ENDIF\n";
      print O "NMOC_RATIO_IN=0.0_dp\n\n";
      print O "! Species from $vocFile\n";

      # Loop over substrings
      for ( $i = 2; $i < scalar( @result ); $i++ ) {

	# Species name
	$name   =  $result[$i];

	# Species index (F90 notation, starting w/ 1)
	$spcStr = &num2str( $i - 1 );

	# Write F90 code to output file
	print O "IN_NMOC_NAME($spcStr)=\"$name\"\n";
      }

      # Print a return
      print O "\n";

    } else {

      #=====================================================================
      # The following lines have the emission factors
      #=====================================================================

      # Split the line by commas
      @result = split( ',', $line );

      # Emission factor number
      # CROPS is listed as #9 but reset to #6
      $fac = $result[0];
      if ( $fac > 6 ) { $fac = 6 }

      # Write comment header
      print O "! $result[1] emission factors from $vocFile\n";

      # Loop over substrings
      for ( $i = 0; $i < $nVoc; $i++ ) {

	# Species index (F90 notation, starting w/ 1)
	$spcStr = &num2str( $i + 1 );

	# Column w/ data for the given species
	$col = $i + 2;

	# Write the F90 commands to the output file
	print O "NMOC_RATIO_IN($spcStr,$fac)=$result[$col]_dp\n";
       }

      # Print a return
      print O "\n";

    }
  }

  #=========================================================================
  # Cleanup and quit
  #=========================================================================

  # Write subroutine header to output file
  print O "!EOC\n";

  # Close files
  close( C );
  close( V );
  close( O );

  # Return status
  return( $? );
}
#EOP
#------------------------------------------------------------------------------
#                  GEOS-Chem Global Chemical Transport Model                  !
#------------------------------------------------------------------------------
#BOP
#
# !IROUTINE: main
#
# !DESCRIPTION: Driver program for finn.pl.  Gets arguments from the command
#  line and starts the process of creating the FINN include file.
#\\
#\\
# !INTERFACE:
#
sub main() {
#
# !REVISION HISTORY:
#  11 Aug 2014 - R. Yantosca - Initial version
#EOP
#------------------------------------------------------------------------------
#BOC

  # Exit w/ error
  if ( scalar( @ARGV ) != 3 ) {
    print "USAGE: gfed3.pl CO2-FILE VOC-FILE OUTFILE\n";
    exit(1);
  }

  # Create the GFED3 include file
  &makeFinnInclude( @ARGV );

  # Return status
  return( $? );
}

#------------------------------------------------------------------------------

# Call main program
main();

# Return status
exit( $? );
#EOC


