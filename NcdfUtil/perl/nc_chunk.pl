#!/usr/bin/perl -w

#------------------------------------------------------------------------------
#                  GEOS-Chem Global Chemical Transport Model                  !
#------------------------------------------------------------------------------
#BOP
#
# !MODULE: nc_chunk.pl
#
# !DESCRIPTION: This perl script can be used to chunk and/or deflate a
#  netCDF file.  It will first do a ncdump command to find out the dimensions
#  of the file, and then chunk the file along lon and lat.  This is optimal
#  for GEOS-Chem and GCHP.
#\\
#\\
# !USES:
#
require 5.003;                       # Need this version of Perl or newer
use English;                         # Use English language
use Carp;                            # Get detailed error messages
use strict;

# !CALLING SEQUENCE:
#  nc_chunk.pl NETCDF-FILE-NAME [DEFLATION-LEVEL]
#
# !REMARKS:
#  Must have the NetCDF Operators (NCO) installed.
#
# !REVISION HISTORY:
#  12 Apr 2018 - R. Yantosca - Initial version
#EOP
#------------------------------------------------------------------------------
#                  GEOS-Chem Global Chemical Transport Model                  !
#------------------------------------------------------------------------------
#BOP
#
# !IROUTINE: chunkTheFile
#
# !DESCRIPTION: Calls nccopy to chunk a netCDF file along lon and lat
#  dimensions.  If the deflation factor is greater than zero, it will also
#  tell nccopy to compress the netCDF file with that deflation level.
#\\
#\\
# !INTERFACE:
#
sub chunkTheFile($$) {
#
# !INPUT PARAMETERS:
#
  my ( $ncFile, $deflate ) = @_;   # NetCDF file name to be chunked
                                   # and the deflation factor
#
# !REVISION HISTORY:
#  12 Apr 2018 - R. Yantosca - Initial version
#EOP
#------------------------------------------------------------------------------
#BOC
#
# !LOCAL VARIABLES:
#
  # Scalars
  my $isDims    = 0;

  # Strings
  my $cmd       = "";
  my $latCmd    = "";
  my $latSize   = "";
  my $lonCmd    = "";
  my $lonSize   = "";
  my $levCmd    = "";
  my $timeCmd   = "";
  my $line      = "";
  my $result    = "";

  # Arrays
  my @ncDumpTxt = ();
  my @subStrs   = ();

  # Exit if the file does not exist
  if ( !( -f $ncFile ) ) { print "Could not find $ncFile!\n"; exit(-1); }

  # Save the contents of ncdump to an array
  @ncDumpTxt = qx/ncdump -cts $ncFile/;

  # Parse the ncdump line
  foreach $line ( @ncDumpTxt ) {

    # Remove newlines
    chomp( $line );

    # We have entered the dimensions section
    if ( $line =~ m/dimensions:/ ) { $isDims = 1; }

    # Look for the lon and lat dimensions
    if ( $isDims > 0 ) {

      # Create chunking command for lon dimension
      # Create the chinking
      if ( $line =~ m/[Ll][0o][Nn]/ ) {
	$lonCmd  =  $line;
	$lonCmd  =~ s/ = /\//g;
	$lonCmd  =~ s/ ;//g;
	$lonCmd  =~ s/^\s+|\s+$//g
      }

      # Create chunking command for lat dimension
      if ( $line =~ m/[Ll][Aa][Tt]/ ) {
	$latCmd  = $line;
	$latCmd  =~ s/ = /\//g;
	$latCmd  =~ s/ ;//g;
	$latCmd  =~ s/^\s+|\s+$//g
      }

      # Create chunking command for level midpoints dimension
      if ( $line =~ m/[Ll][Ee][Vv]/ ) {
	@subStrs = split( ' = ', $line );
	$levCmd  = "$subStrs[0]/1";
	$levCmd  =~ s/^\s+|\s+$//g
      }

      # Create chunking command for level interfaces dimension
      if ( $line =~ m/[Ii][Ll][Ee][Vv]/ ) {
	@subStrs = split( ' = ', $line );
	$levCmd  = "$subStrs[0]/1";
	$levCmd  =~ s/^\s+|\s+$//g
      }

      # Create chunking command for time dimension
      if ( $line =~ m/[Tt][Ii][Mm][Ee]/ ) {
	@subStrs = split( ' = ', $line );
	$timeCmd = "$subStrs[0]/1";
	$timeCmd =~ s/^\s+|\s+$//g
      }

    }

    # We have exited the dimensions section
    if ( $line =~ m/variables:/  ) { goto quit; }
  }

quit:

  # Chunking command
  if ( $levCmd eq "" && $timeCmd eq "" )
    { $cmd = "nccopy -c $lonCmd,$latCmd";                  }
  elsif ( $levCmd eq "" )
    { $cmd = "nccopy -c $lonCmd,$latCmd,$timeCmd";         }
  else
    { $cmd = "nccopy -c $lonCmd,$latCmd,$levCmd,$timeCmd"; }

  # Deflation command
  if ( $deflate > 0 ) { $cmd .= " -d$deflate" }

  # Add file name to the command
  $cmd .= " $ncFile tmp.nc; mv tmp.nc $ncFile";
  print "$cmd\n";

  # Execute the command
  $result = qx/$cmd/;
  chomp( $result );
  if ( $result ne "" ) { print "$result\n"; }

  # Exit and pass status code back
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
# !DESCRIPTION: Driver program for the nc_chunk.pl script.
#\\
#\\
# !INTERFACE:
#
sub main(@) {
#
# !REVISION HISTORY:
#  12 Apr 2018 - R. Yantosca - Initial version
#EOP
#------------------------------------------------------------------------------
#BOC

  # Error message
  my $errMsg = "Usage: nc_chunk.pl NETCDF-FILE-NAME";

  # If the user passes a filename from the command line, use it
  # Otherwise, default to "UnitTest.input"
  if    ( scalar( @ARGV ) == 2 ) { &chunkTheFile( @ARGV,   );     }
  elsif ( scalar( @ARGV ) == 1 ) { &chunkTheFile( @ARGV, 0 );     }
  else                           { print "$errMsg\n"; exit( -1 ); }

  # Exit and pass status code backq
  return( $? );
}
#EOC

#------------------------------------------------------------------------------

# Call main program
main();

# Exit and pass status code back to Unix shell
exit( $? );
