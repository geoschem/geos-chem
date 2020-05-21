#!/usr/bin/perl -w
#------------------------------------------------------------------------------
#                  GEOS-Chem Global Chemical Transport Model                  !
#------------------------------------------------------------------------------
#BOP
#
# !MODULE: gfed.pl
#
# !DESCRIPTION: Reads the GFED3_emission_factors.txt file and converts
#  the data in the file to hardwired F90 commands for inlining into
#  an include file.  This facilitates I/O in the ESMF environment.
#\\
#\\
# !INTERFACE:
#
# !USES:
#
require 5.003;       # need this version of Perl or newer
use English;         # Use English language
use Carp;            # Strict error checking
use strict;          # IMPLICIT NONE style syntax
#
# !PUBLIC MEMBER FUNCTIONS:
#  &getProTeXHeader  : Returns subroutine header
#  &makeGfed3Include : Writes the include file for GFED3 w/ F90 commands
#  &main             : Driver program
#
# !REMARKS:
#  You need to specify the root emissions data directory by setting
#  the $HEMCO_DATA_ROOT environment variable accordingly.
#
# !CALLING SEQUENCE:
#  gfed.pl
#
# !REVISION HISTORY:
#  11 Aug 2014 - R. Yantosca - Initial version
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
! !IROUTINE: hcox_gfed_include.H
!
! !DESCRIPTION: Include file with GFED emission factor data that was
!  originally contained in file GFED\\_emission\\_factors.txt.  We have now
!  transformed this file into hardwired F90 commands in order to avoid reading
!  an ASCII file in the ESMF environment.
!
! !REMARKS:
!  ABOUT THIS FILE:
!  ----------------
!  This file was created by script HEMCO/Extensions/Preprocess/gfed.pl.
!  This script can be executed with the following command:
!
!    cd HEMCO/Extensions/Preprocess
!    make gfed
!
!  This will regenerate this include file from the original data and
!  automatically place it in the HEMCO/Extensions directory.
!
!  White space has been removed in order to reduce the file size as much
!  as possible.  If you have to recreate this file, then it is easier to
!  generate via the Perl script than to try to hand edit the code below.
!
!  DATA:
!  -----
!  The GFED_EMFAC array contains emission factors in kg/kgDM or kgC/kgDM
!  GFED_EMFAC(N,1) = Agricultural Waste   Emission Factor for species N
!  GFED_EMFAC(N,1) = Deforestation        Emission Factor for species N
!  GFED_EMFAC(N,1) = Extratropical Forest Emission Factor for species N
!  GFED_EMFAC(N,1) = Peat                 Emission Factor for species N
!  GFED_EMFAC(N,1) = Savanna              Emission Factor for species N
!  GFED_EMFAC(N,1) = Woodland             Emission Factor for species N
!
! !REVISION HISTORY:
!  08 Aug 2014 - R. Yantosca - Initial version
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
# !IROUTINE: makeGfed3Include
#
# !DESCRIPTION: Reads the GFED_emission_factors.txt and writes out
#  equivalent F90 assignment statements commands.
#\\
#\\
# !INTERFACE:
#
sub makeGfed3Include($$) {
#
# !INPUT PARAMETERS:
#
  # Input and output files
  my( $inFile, $outFile ) = @_;
#
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
#
# !REVISION HISTORY:
#  11 Aug 2014 - R. Yantosca - Initial version
#EOP
#------------------------------------------------------------------------------
#BOC
#
# !LOCAL VARIABLES:
#
  # Strings
  my $line    = "";

  # Numbers
  my $n       = 0;

  # Arrays
  my @result  = ();
  my @header  = &getProTeXHeader();

  #=========================================================================
  # Open files and write subroutine header to output file
  #=========================================================================

  # Open input and output files
  open( I, "$inFile"   ) or croak( "Cannot open $inFile!\n" );
  open( O, ">$outFile" ) or croak( "Cannot open $outFile\n" );

  # Write subroutine header to output file
  foreach $line ( @header ) { print O "$line"; }

  #=========================================================================
  # Parse input file and write F90 commands to output file
  #=========================================================================

  # Get each line of input
  foreach $line ( <I> ) {

    # Remove newline characters
    chomp( $line );

    # Parse input lines
    if    ( $line =~ m/Emission/ ) { ; }  # Skip header line
    elsif ( $line =~ m/Column/   ) { ; }  # Skip header line
    else                           {

      # Split the line
      @result = split( ' ', $line );

      # The first substring is the species index
      $n      = $result[0];

      # Write out the F90 commands to initialize the arrays
      print O "\n! $result[1]\n";
      print O "GFED_SPEC_NAME($n)=\"$result[1]\"\n";
      print O "GFED_EMFAC($n,1)=$result[2]_hp\n";
      print O "GFED_EMFAC($n,2)=$result[3]_hp\n";
      print O "GFED_EMFAC($n,3)=$result[4]_hp\n";
      print O "GFED_EMFAC($n,4)=$result[5]_hp\n";
      print O "GFED_EMFAC($n,5)=$result[6]_hp\n";
      print O "GFED_EMFAC($n,6)=$result[7]_hp\n";

    }
  }

  #=========================================================================
  # Cleanup and quit
  #=========================================================================

  # Write closing !EOC tag to include file
  foreach $line ( @header ) { print O "!EOC\n"; }

  # Close files
  close( I );
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
# !DESCRIPTION: Driver program for gfed.pl.  Gets arguments from the command
#  line and starts the process of creating the GFED include file.
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
  if ( scalar( @ARGV ) != 2 ) {
    print "USAGE: gfed.pl INFILE OUTFILE\n";
    exit(1);
  }

  # Create the GFED include file
  &makeGfed3Include( @ARGV );

  # Return status
  return( $? );
}

#------------------------------------------------------------------------------

# Call main program
main();

# Return status
exit( $? );
#EOC


