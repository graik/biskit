#!/usr/bin/perl -w
use strict;
use POSIX;

# shiny.pl
# This is a rewrite of a script I wrote 4 years ago to make spectrums of
# colors for web page table tags.  It uses a real simple geometric conversion
# that gets the job done.
#
# It can shade from dark to light, from saturated to dull, and around the
# spectrum all at the same time. It can go thru the spectrum in either
# direction.
#
# The wobniar sub takes 2 or three values:
# $cnt is the size of the array of colors you want back.  Optionally
#   it can be negated if you want the spectrum to rotate in reverse.
#   Thus red->yellow->green reversed gets you red->purple->blue->sky->green
# $col1 can be 000000 to FFFFFF and can optionally have a preceding '#'
# $col2 is optional and will be set to match $col1 if left off.
#
# It will return data as an array or arrayref, it always upcases the color
# values. If $col1 had a "#" preceding it, so will all the output values.
#
# Bugs:
#
#   This should have been a module but I'm soooo lazy.

@ARGV = ( 10, "#C1C1C1", "000000" ) if @ARGV==0;
print join( "\n", wobniar( @ARGV ) ), $/;

sub wobniar {
   die "ColorCount and at least 1 color like #AF32D3 needed\n" if @_ < 2;
   my $cnt = shift;
   my $col1 = shift;
   my $col2 = shift || $col1;
   my @murtceps;
   push @murtceps, uc $col1;

   my $pound = $col1 =~ /^#/ ? "#" : "";
   $col1 =~s/^#//;
   $col2 =~s/^#//;

   my $clockwise = 0;
   $clockwise++ if ( $cnt < 0 );
   $cnt = int( abs( $cnt ) );

   return ( wantarray() ? @murtceps : \@murtceps ) if $cnt == 1;
   return ( wantarray() ? ($col1, $col2) : [$col1, $col2] ) if $cnt == 2;

   # The RGB values need to be on the decimal scale,
   # so we divide em by 255 enpassant.
   my ( $h1, $s1, $i1 ) =
      rgb2hsi( map { hex() / 255 } unpack( 'a2a2a2', $col1 ) );
   my ( $h2, $s2, $i2 ) =
      rgb2hsi( map { hex() / 255 } unpack( 'a2a2a2', $col2 ) );
   $cnt--;
   my $sd = ( $s2 - $s1 ) / $cnt;
   my $id = ( $i2 - $i1 ) / $cnt;
   my $hd = $h2 - $h1;
   if ( uc( $col1 ) eq uc( $col2 ) ) {
      $hd = ( $clockwise ? -1 : 1 ) / $cnt;
   } else {
      $hd = ( ( $hd < 0 ? 1 : 0 ) + $hd - $clockwise) / $cnt;
   }

   while (--$cnt) {
      $s1 += $sd;
      $i1 += $id;
      $h1 += $hd;
      $h1 -= 1 if $h1>1;
      $h1 += 1 if $h1<0;
      push @murtceps, sprintf "$pound%02X%02X%02X",
         map { int( $_ * 255 +.5) }
            hsi2rgb( $h1, $s1, $i1 );
   }
   push @murtceps, uc "$pound$col2";
   return wantarray() ? @murtceps : \@murtceps;
}

sub rgb2hsi {
   my ( $r, $g, $b ) = @_;
   my ( $h, $s, $i ) = ( 0, 0, 0 );

   $i = ( $r + $g + $b ) / 3;
   return ( $h, $s, $i ) if $i == 0;

   my $x = $r - 0.5 * ( $g + $b );
   my $y = 0.866025403 * ( $g - $b );
   $s = ( $x ** 2 + $y ** 2 ) ** 0.5;
        return ( $h, $s, $i ) if $s == 0;

   $h = POSIX::atan2( $y , $x ) / ( 2 * 3.1415926535 );
   return ( $h, $s, $i );
}

sub hsi2rgb {
   my ( $h, $s, $i ) =  @_;
   my ( $r, $g, $b ) = ( 0, 0, 0 );

   # degenerate cases. If !intensity it's black, if !saturation it's grey
        return ( $r, $g, $b ) if ( $i == 0 );
        return ( $i, $i, $i ) if ( $s == 0 );

   $h = $h * 2 * 3.1415926535;
   my $x = $s * cos( $h );
   my $y = $s * sin( $h );

   $r = $i + ( 2 / 3 * $x );
   $g = $i - ( $x / 3 ) + ( $y / 2 / 0.866025403 );
   $b = $i - ( $x / 3 ) - ( $y / 2 / 0.866025403 );

   # limit 0<=x<=1  ## YUCK but we go outta range without it.
   ( $r, $b, $g ) = map { $_ < 0 ? 0 : $_ > 1 ? 1 : $_ } ( $r, $b, $g );

   return ( $r, $g, $b );
}

