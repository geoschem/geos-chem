#!/usr/bin/perl -w

#------------------------------------------------------------------------------
#                  GEOS-Chem Global Chemical Transport Model                  !
#------------------------------------------------------------------------------
#BOP
#
# !MODULE: validate.pl
#
# !DESCRIPTION: This perl script validates GEOS-Chem output files f rom
#  either unit tests or difference tests.  It compares the md5sums of
#  2 corresponding files (*.sp vs. *.mp, or *.Ref vs. *.Dev) and then
#  creates user-friendly output.
#\\
#\\
# !USES:
#
require 5.003;                   # Need this version of Perl or newer
use English;                     # Use English language
use Carp;                        # Get detailed error messages
use strict;                      # Use "IMPLICIT NONE" syntax
use Cwd;                         # Get current directories
#		
# !PUBLIC MEMBER FUNCTIONS:
#  &main()    : Driver routine for gcUnitTest
#
# !PRIVATE MEMBER FUNCTIONS:
#  &getCheckSums   : Returns the name of the GEOS-Chem run directory
#  &checkTheFiles  : Returns the name of the GEOS-Chem run directory
#
# !CALLING SEQUENCE:
#  ./validate.pl SIMULATION DIR1 DIR2 TYPE
#                                                                             . 
#  ./validate.pl geosfp_4x5_benchmark . . ut           # Unit test example
#                                                                             .
#  ./validate.pl geosfp_4x5_benchmark ./Ref ./Dev dt   # Diff test example
#
# !REMARKS:
#  Designed to be used with the GEOS-Chem Unit Tester and Diff Tester
#  Always returns a value of zero so as not to disrupt the make process.
#
# !REVISION HISTORY: 
#  15 Feb 2018 - R. Yantosca - Initial version
#EOP
#------------------------------------------------------------------------------
#BOP
#
# !IROUTINE: getCheckSums
#
# !DESCRIPTION: Tests if two files are identical, via the md5sum hash.
#  Returns a status string to routine &checkTheFiles.
#\\
#\\
# !INTERFACE:
#
sub getCheckSums($$$$) {
#
# !INPUT PARAMETERS:
#
  # $file1    : First file to be compared 
  # $file2    : Second file to be compared
  # $combName : Combined file name for display purposes
  my ( $file1, $file2, $combName ) = @_;
#
# !REVISION HISTORY:
#  15 Feb 2018 - R. Yantosca - Initial version
#EOP
#------------------------------------------------------------------------------
#BOC
#
# !LOCAL VARIABLES:
#
  # Scalars
  my $ckSum1   = "MISSING";
  my $ckSum2   = "MISSING";
  my $status   = "";
  my $result   = "";

  # Arrays
  my @subStr   = ();

  # If the 1st file exists, get the md5sum
  if ( -f $file1 ) {
    chomp( $result = qx/md5sum $file1/ );
    @subStr = split( ' ', $result );
    $ckSum1 = $subStr[0];
  }

  # If the 2nd file exists, get the md5sum
  if ( -f $file2 ) {
    chomp( $result = qx/md5sum $file2/ );
    @subStr = split( ' ', $result );
    $ckSum2 = $subStr[0];
  }

  # Create the return string
  if    ( $ckSum1 eq "MISSING" ) { $status = "MISSING   : $file1";    }
  elsif ( $ckSum2 eq "MISSING" ) { $status = "MISSING   : $file2";    }
  elsif ( $ckSum1 eq $ckSum2   ) { $status = "IDENTICAL : $combName"; }
  else                           { $status = "DIFFERENT : $combName"; }

  # Return status string to calling program
  return( $status );
}
#EOP
#------------------------------------------------------------------------------
#                  GEOS-Chem Global Chemical Transport Model                  !
#------------------------------------------------------------------------------
#BOP
#
# !IROUTINE: checkTheFiles
#
# !DESCRIPTION: Tests if corresponding pairs of bpch and/or netCDF files 
#  in the specified directories are identical. 
#\\
#\\
# !INTERFACE:
#
sub checkTheFiles($$$$) {
#
# !INPUT PARAMETERS:
#
  # $sim  : Simulation identifier string (e.g. geosfp_4x5_benchmark)
  # $dir1 : Folder where the first file to be checked resides
  # $dir2 : Folder where the second file to be checked resides
  # $type : "ut" for Unit Tests or "dt" for difference tests
  my ( $sim, $dir1, $dir2, $type )= @_;
#
# !REVISION HISTORY:
#  15 Feb 2018 - R. Yantosca - Initial version
#EOP
#------------------------------------------------------------------------------
#BOC
#
# !LOCAL VARIABLES:
#
  # Scalars
  my $baseName = "";
  my $combName = "";
  my $ext1     = "";
  my $ext2     = "";
  my $file     = "";
  my $file1    = "";
  my $file2    = "";
  my $fileName = "";
  my $suf1     = "";
  my $suf2     = "";
  my $status   = "";  
  my $errMsg   = 
   "validate.pl: TYPE must be either DT or UT (case-insensitive)!  Exiting...";

  # Arrays
  my @files   = ();

  #=========================================================================
  # Initialize
  #=========================================================================

  # Pick the proper suffixes for the unit teset or difference teste
  if    ( $type =~ m/[Uu][Tt]/ ) { $ext1  = "sp";   $ext2  = "mp";  
 				   $suf1  = ".sp";  $suf2  = ".mp";  }
  elsif ( $type =~ m/[Dd][Tt]/ ) { $ext1  = "Ref";  $ext2  = "Dev";  
 				   $suf1  = ".Ref"; $suf2  = ".Dev"; }
  else                           { print "$errMsg\n"; exit( -1 );    }

  # Read files in the directory
  opendir( my $dh, $dir1 ) || die( "validate.pl: Can't open $dir1\n" );
  @files = readdir( $dh );
  closedir( $dh );
  
  #=========================================================================
  # Print the top header
  #=========================================================================
  print '#' x 79; 
  print "\n### VALIDATION OF GEOS-CHEM OUTPUT FILES\n";
  print "### Run ID: $sim\n";
  print "##@\n";

  #=========================================================================
  # Test if files are identical
  #=========================================================================
  foreach $file ( sort( @files ) ) {

    # Strip new lines
    chomp( $file );

    # Only consider files that are bpch or netCDF
    if ( ( $file =~ m/trac_avg.$sim/ ) || ( $file =~ m/\.nc/ ) ){
      
      # Look for the first file
      if ( $file =~ m/$suf1/ ) {

	# First file to examine
	$file1 =  "$dir1/$file";
	    
	# Second file to examine
	$file2 =  "$dir2/$file";
	if    ( $type =~ m/[Dd][Tt]/ ) { $file2 =~ s/$suf1/$suf2/g;  }
	elsif ( $type =~ m/[Uu][Tt]/ ) { $file2 =~ s/\$suf1/$suf2/g; }

	# Create combined file name for display
	$baseName = substr( $file, 0, length($file) - length($suf1) );
	$combName = "$baseName.\{$ext1,$ext2}";
 
	# Test if files are identical and write status to stdout
	$status = &getCheckSums( $file1, $file2, $combName );
	print "### $status\n";

      }	   
    } 
  }

  #=========================================================================
  # Now repeat for files in OutputDir
  #=========================================================================
  my $outdir1 = "$dir1/OutputDir";
  my $outdir2 = "$dir2/OutputDir";

  # Read files in the directory
  opendir( $dh, $outdir1 ) || die( "validate.pl: Can't open $outdir1\n" );
  @files = readdir( $dh );
  closedir( $dh );
  
  # Test if files are identical
  foreach $file ( sort( @files ) ) {

    # Strip new lines
    chomp( $file );

    # Only consider files that are netCDF
    if ( ( $file =~ m/\.nc/ ) ){
      
      # Look for the first file
      if ( $file =~ m/$suf1/ ) {

	# First file to examine
	$file1 =  "$outdir1/$file";
	    
	# Second file to examine
	$file2 =  "$outdir2/$file";
	if    ( $type =~ m/[Dd][Tt]/ ) { $file2 =~ s/$suf1/$suf2/g;  }
	elsif ( $type =~ m/[Uu][Tt]/ ) { $file2 =~ s/\$suf1/$suf2/g; }

	# Create combined file name for display
	$baseName = substr( $file, 0, length($file) - length($suf1) );
	$combName = "$baseName.\{$ext1,$ext2}";
 
	# Test if files are identical and write status to stdout
	$status = &getCheckSums( $file1, $file2, $combName );
	print "### $status\n";

      }	   
    } 
  }

  #=========================================================================
  # Print the bottom header
  #=========================================================================
  print "##\%\n";
  print '#' x 79 . "\n";

  # Return error code to main program
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
# !DESCRIPTION: Driver routine for the validate.pl script.
#\\
#\\
# !INTERFACE:
#
sub main() {
#
# !REVISION HISTORY:
#  15 Feb 2018 - R. Yantosca - Initial version
#EOP
#------------------------------------------------------------------------------
#BOC
#
# !LOCAL VARIABLES:
#
  # Scalars
  my $errMsg = "Usage: validate.pl SIM DIR1 DIR2 TYPE\n";

  # If the proper # of args are passed, then validate the GEOS-Chem files
  # from a unit test or difference test.  Otherwise exit with error.
  if    ( scalar( @ARGV ) == 4 ) { &checkTheFiles( @ARGV );   }
  else                           { print "$errMsg"; exit(1);  }

  # Return success so as not to halt the make process
  return( 0 );
}

#------------------------------------------------------------------------------

# Call the main program
main();

# Always return zero so we don't halt the make process
exit( 0 );
