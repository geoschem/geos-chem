#!/usr/bin/perl -w
#------------------------------------------------------------------------------
#                  GEOS-Chem Global Chemical Transport Model                  !
#------------------------------------------------------------------------------
#BOP
#
# !MODULE: gfed3.pl
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
require 5.003;      # need this version of Perl or newer
use English;        # Use English language
use Carp;           # Strict error checking
use strict;         # IMPLICIT NONE style syntax
#
# !PUBLIC MEMBER FUNCTIONS:
#  &getProTeXHeader : Returns subroutine header
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
! !IROUTINE: hcox_gfed3_include.H
!
! !DESCRIPTION: Include file with GFED3 emission factor data that was 
!  originally contained in file GFED3\_emission\_factors.txt.  We have now
!  transformed this file into hardwired F90 commands in order to avoid reading
!  an ASCII file in the ESMF environment.
!
! !REMARKS:
!  This file was generated with HEMCO/Extensions/Preprocess/gfed3.pl.
!
!  The GFED3_EMFAC array contains emission factors in kg/kgDM or kgC/kgDM
!  GFED3_EMFAC(N,1) = Agricultural Waste   Emission Factor for species N
!  GFED3_EMFAC(N,1) = Deforestation        Emission Factor for species N
!  GFED3_EMFAC(N,1) = Extratropical Forest Emission Factor for species N
!  GFED3_EMFAC(N,1) = Peat                 Emission Factor for species N
!  GFED3_EMFAC(N,1) = Savanna              Emission Factor for species N       
!  GFED3_EMFAC(N,1) = Woodland             Emission Factor for species N
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
# !MODULE: main
#
# !DESCRIPTION: Reads the GFED3_emission_factors.txt and writes out
#  equivalent F90 assignment statements commands.
#\\
#\\
# !INTERFACE:
#
sub main() {
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
# !REVISION HISTORY: 
#  11 Aug 2014 - R. Yantosca - Initial version
#EOP
#------------------------------------------------------------------------------
#BOC
#
# !LOCAL VARIABLES:
#
  # Strings
  my $rootDir = "$ENV{'HEMCO_DATA_ROOT'}";
  my $inFile  = "$rootDir//GFED3/v2014-07/GFED3_emission_factors.txt";
  my $outFile = "hcox_gfed3_include.H";
  my $line    = "";

  # Numbers
  my $n       = 0;

  # Arrays
  my @result  = ();

  # Subroutine header
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
    if    ( $line =~ m/Emission/ ) { print O "    ! $line\n"; }  # Header line
    elsif ( $line =~ m/Column/   ) { print O "    ! $line\n"; }  # Header line
    else                           {

      # Split the line
      @result = split( ' ', $line );

      # The first substring is the species index
      $n      = $result[0];

      # Write out the F90 commands to initialize the arrays
      print O "\n    ! $result[1]\n";
      print O "    GFED3_SPEC_NAME($n  ) = \"$result[1]\"\n";
      print O "    GFED3_EMFAC    ($n,1) = $result[2]_hp\n";
      print O "    GFED3_EMFAC    ($n,2) = $result[3]_hp\n";
      print O "    GFED3_EMFAC    ($n,3) = $result[4]_hp\n";
      print O "    GFED3_EMFAC    ($n,4) = $result[5]_hp\n";
      print O "    GFED3_EMFAC    ($n,5) = $result[6]_hp\n";
      print O "    GFED3_EMFAC    ($n,6) = $result[7]_hp\n";

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

#------------------------------------------------------------------------------


# Call main program
main();

# Return status
exit( $? );
#EOC


