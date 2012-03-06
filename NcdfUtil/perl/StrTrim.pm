#!/usr/bin/perl -w

#------------------------------------------------------------------------------
#          Harvard University Atmospheric Chemistry Modeling Group            !
#------------------------------------------------------------------------------
#BOP
#
# !MODULE: StrTrim
#
# !DESCRIPTION: This Perl package contains routines for splitting a line
#  into substrings and removing trailing and leading whitespace.  Used by
#  the ncCode* scripts.
#\\
#\\
# !INTERFACE:
#
package StrTrim;
#
# !USES:
#
  require 5.003;   # Need this version of Perl or newer
  use English;     # Use English language
  use Carp;        # Get detailed error messages
  use strict;      # Force explicit variable declarations (like IMPLICIT NONE)
#
#
# !PUBLIC MEMBER FUNCTIONS:
#  &trim($)
#  &splitLine($$)
#
# !CALLING SEQUENCE:
#  use StrTrim qw( trim splitLine );
#
# !REVISION HISTORY: 
#  30 Jan 2012 - R. Yantosca - Initial version
#EOP
#------------------------------------------------------------------------------
#BOC
BEGIN {

  #=========================================================================
  # The BEGIN method lists the names to export to the calling routine
  #=========================================================================
  use Exporter ();
  use vars     qw( $VERSION @ISA @EXPORT_OK );

  $VERSION   = 1.00;                                   # version number
  @ISA       = qw( Exporter );                         # export method
  @EXPORT_OK = qw( &trim &splitLine );
}
#EOC
#------------------------------------------------------------------------------
#          Harvard University Atmospheric Chemistry Modeling Group            !
#------------------------------------------------------------------------------
#BOP
#
# !IROUTINE: trim
#
# !DESCRIPTION:  Routine trim removes leading and trailing whitespace from
#  a string (analogous to IDL's Strtrim( str, 2 ) command).
#\\
#\\
# !INTERFACE:
#
sub trim($) {
#
# !CALLING SEQUENCE:
# $string = &trim( $string );
#
# !REMARKS:
#  Found online at this URL:
#  http://www.somacon.com/p114.php
#
# !REVISION HISTORY:
#  27 Jan 2012 - R. Yantosca - Initial version
#EOP
#------------------------------------------------------------------------------
#BOC

  # Shift the @_ array 
  my $string = shift;

  # Remove leading whitespace
  $string    =~ s/^\s+//;

  # Remove trailing whitespace
  $string    =~ s/\s+$//;

  # Return
  return( $string );
}
#EOP
#------------------------------------------------------------------------------
#          Harvard University Atmospheric Chemistry Modeling Group            !
#------------------------------------------------------------------------------
#BOP
#
# !IROUTINE: splitLine
#
# !DESCRIPTION: Routine splitLine splits a line on a given delimiter, and
#  strips white space.  Convenience wrapper for the Perl "split" function.
#\\
#\\
# !INTERFACE:
#
sub splitLine($$) {
#
# !INPUT PARAMETERS:
#
  # Line to be split, and the delimeter character
  # Don't strip the white from $value if $noSplitVal==1
  my( $line, $delim ) = @_;
#
# !CALLING SEQUENCE:
# ( $name, $value ) = &splitLine( $line );
#
# !REVISION HISTORY:
#  27 Jan 2012 - R. Yantosca - Initial version
#EOP
#------------------------------------------------------------------------------
#BOC
#
# !LOCAL VARIABLES:
#
  # Split the line
  my @subStr = split( $delim, $line );
  my $name   = &trim( $subStr[0] );
  my $value  = &trim( $subStr[1] );
 
  # Return substrings
  return( $name, $value );
}
#EOC
END {}
