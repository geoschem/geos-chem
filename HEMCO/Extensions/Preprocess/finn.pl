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
use strict;         # Force explicit variable declarations (like IMPLICIT NONE)
#
# !PUBLIC MEMBER FUNCTIONS:
#  &main          : Driver program
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
#EOP
#------------------------------------------------------------------------------
#                  GEOS-Chem Global Chemical Transport Model                  !
#------------------------------------------------------------------------------
#BOP
#
# !MODULE: getProTeXHeader
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
!  originally contained in files FINN\_EFratios\_CO2.csv and
!  FINN\_VOC\_speciation.csv.  We transform these data into hardwired F90
!  assignment statements so that we can avoid reading ASCII files in the 
!  ESMF environment.
!
! !REMARKS:
!  This file was generated with script HEMCO/Extensions/Preprocess/finn.pl.
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
# !MODULE: main
#
# !DESCRIPTION: Reads the GFED3_emission_factors.txt and writes out
#  equivalent F90 assignment statements commands.
#\\
#\\
# !REMARKS:
#  The format of the input file is:
#  Emission factors in kg/kgDM or kgC/kgDM
#  Column 1 = Species Number
#  Column 2 = Species Name
#  Column 3 = Ag Waste Emission Factor
#  Column 4 = Deforestation Emission Factor
#  Column 5 = Extratrop Forest Emission Factor
#  Column 6 = Peat Emission Factor
#  Column 7 = Savanna Emission Factor
#  Column 8 = Woodland Emission Factor
# 
# !REVISION HISTORY: 
#  17 Apr 2014 - R. Yantosca - Initial version
#EOP
#------------------------------------------------------------------------------
#BOC
#
# !LOCAL VARIABLES:
#
sub main() {

  # Strings
  my $rootDir = "$ENV{'HEMCO_DATA_ROOT'}";
  my $co2File = "FINN_EFratios_CO2.csv";
  my $co2Path = "$rootDir/FINN/v2014-07/$co2File";
  my $vocFile = "FINN_VOC_speciation.csv";
  my $vocPath = "$rootDir/FINN/v2014-07/$vocFile";
  my $outFile = "hcox_finn_include.H";
  my $line    = "";
  my $subStr  = "";
  my $name    = "";
  my $spcStr  = "";
  
  # Numbers
  my $i       = 0;
  my $spc     = 0;
  my $fac     = 0;
  my $n       = 0;
  my $col     = 0;

  # Arrays
  my @result  = ();
  my @co2Spec = ( "CO2", "CO", "CH4", "H2", "NOx", "SO2", 
                  "OC",  "BC", "NH3", "NO", "NO2",       );
  my @co2Ind  = (  1, 2, 3, 58, 4, 5, 6, 7, 8, 9, 10 );
  my @vocInd  = ( 11 .. 57 );

  # Get subroutine header
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

      # Print comment lines
      print O "      !---------------------------------------------------\n";
      print O "      ! Data from $co2File\n";
      print O "      !---------------------------------------------------\n\n";
      print O "      ! Species from $co2File\n";

      # Loop over substrings
      for ( $i = 2; $i < scalar( @result ); $i++ ) {
	$name =  "$result[$i]";
	$name =~ s/CO2\///g;

	# Skip NMOC
	if ( !( $name =~ m/NMOC/ ) ) {

	  # Increment the species counter
	  $spc  = $co2Ind[$n++];

	  # Pad spaces (cheap algorithm)
	  if   ( $spc < 10 ) { $spcStr = " $spc"; }
	  else               { $spcStr = "$spc";  }		

	  # Write F90 code to output file
	  print O "      FINN_SPEC_NAME($spcStr) = \"$name\"\n";
        }
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
      print O "      ! $result[1] emission factors from $co2File\n";
	
      # Loop over substrings
      for ( $i = 0; $i < scalar( @co2Ind ); $i++ ) {

	# Species index
	$spc = $co2Ind[$i];

	# Pad spaces (cheap algorithm)
	if   ( $spc < 10 ) { $spcStr = " $spc"; }
	else               { $spcStr = "$spc";  }

	# Column count
	$col = $i + 2;
	
	# Write the F90 commands to the output file
	print O "      EMFAC_IN($spcStr,$fac) = $result[$col]_dp\n";
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

      # Reset variables
      $n       = 0;

      # Split the line by commas
      @result = split( ',', $line );

      # Print comment lines
      print O "      !---------------------------------------------------\n";
      print O "      ! Data from $vocFile\n";
      print O "      !---------------------------------------------------\n\n";
      print O "      ! Species from $vocFile\n";

      # Loop over substrings
      for ( $i = 2; $i < scalar( @result ); $i++ ) {
	$name =  $result[$i];
	$name =~ s/CO2\///g;

	# Increment the species counter
	$spc    = $vocInd[$n++];

	# Pad spaces (cheap algorithm)
	if   ( $spc < 10 ) { $spcStr = " $spc"; }
	else               { $spcStr = "$spc";  }

	# Write F90 code to output file
	print O "      FINN_SPEC_NAME($spcStr) = \"$name\"\n";
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
      print O "      ! $result[1] emission factors from $vocFile\n";
	
      # Loop over substrings
      for ( $i = 0; $i < scalar( @vocInd ); $i++ ) {

	# Species index
	$spc = $vocInd[$i];

	# Pad spaces (cheap algorithm)
	if   ( $spc < 10 ) { $spcStr = " $spc"; }
	else               { $spcStr = "$spc";  }

	# Column index
	$col = $i + 2;

	# Write the F90 commands to the output file
	print O "      EMFAC_IN($spcStr,$fac) = $result[$col]_dp\n";
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

#------------------------------------------------------------------------------

# Call main program
main();

# Return status
exit( $? );
#EOC


