package Geo::DNA;

use common::sense;

use Math::Trig qw( :pi rad2deg asin deg2rad );
use POSIX "fmod";
use Data::Dumper;

use Exporter 'import';

our @EXPORT_OK = qw(
    encode_geo_dna
    decode_geo_dna
    neighbours_geo_dna
);

our $VERSION = "1.0";

# """Encode the given latitude and longitude into a geoprint that
# contains 'precision' characters.
#
# If radians is True, input parameters are in radians, otherwise
# they are degrees.  Example::
#
# >>> c = (7.0625, -95.677068)
# >>> h = encode(*c)
# >>> h
# 'watttatcttttgctacgaagt'
#
# >>> r = rads(c[0]), rads(c[1])
# >>> h2 = encode(*r, radians=True)
# >>> h == h2
# True
# """

my $RADIUS_OF_EARTH = 6378100;
my $ALPHABET = [ "g", "a", "t", "c", ];
my $DECODE_MAP = {
    'g' => 0,
    'a' => 1,
    't' => 2,
    'c' => 3,
};

sub encode_geo_dna {
    my ( @args ) = @_;
    encode( @args );
}

sub encode {
    my ( $latitude, $longitude, @opts ) = @_;
    my $options = { @opts };
    my $precision = $options->{precision} || 22;
    my $radians   = $options->{radians}   || 0;

    my $geodna = '';
    my ( $loni, $lati );

    if ( $radians ) {
        $latitude  = rad2deg( $latitude );
        $longitude = rad2deg( $longitude );
    }

    if ( $longitude < 0 ) {
        $geodna .= 'w';
        $loni = [ -180.0, 0.0 ];
    } else {
        $geodna .= 'e';
        $loni = [ 0.0, 180.0 ];
    }

    $lati = [ -90.0, 90.0 ];

    while ( length( $geodna ) < $precision ) {
        my $ch = 0;

        my $mid = ( $loni->[0] + $loni->[1] ) / 2;
        if ( $longitude > $mid ) {
            $ch |= 2;
            $loni = [ $mid, $loni->[1] ];
        } else {
            $loni = [ $loni->[0], $mid ];
        }

        $mid = ( $lati->[0] + $lati->[1] ) / 2;
        if ( $latitude > $mid ) {
            $ch |= 1;
            $lati = [ $mid, $lati->[1] ];
        } else {
            $lati = [ $lati->[0], $mid ];
        }

        $geodna .= $ALPHABET->[$ch];
    }
    return $geodna;
}


# """Decode a geoprint, returning the latitude and longitude.  These
# coordinates should approximate the input coordinates within a
# degree of error returned by 'error()'
#
# >>> c = (7.0625, -95.677068)
# >>> h = encode(*c)
# >>> c2 = decode(h)
# >>> e = error(h)
# >>> abs(c[0] - c2[0]) <= e
# True
# >>> abs(c[1] - c2[1]) <= e
# True
#
# If radians is True, results are in radians instead of degrees.
#
# >>> c2 = decode(h, radians=True)
# >>> e = error(h, radians=True)
# >>> abs(rads(c[0]) - c2[0]) <= e
# True
# >>> abs(rads(c[1]) - c2[1]) <= e
# True
# """

sub decode_geo_dna {
    my ( @args ) = @_;
    decode( @args );
}

sub decode {
    my ( $geodna, @opts ) = @_;
    my $options = { @opts };

    my ( $lati, $loni ) = bounding_box( $geodna );

    my $lat = ( $lati->[0] + $lati->[1] ) / 2;
    my $lon = ( $loni->[0] + $loni->[0] ) / 2;
    if ( $options->{radians} ) {
        return ( deg2rad( $lat ), deg2rad( $lon ) );
    }
    return ( $lat, $lon );
}

sub size {
    my ( $geo_dna ) = @_;
    return 20000000.0 * ( 2 ^ ( ( length( $geo_dna ) * -1.0 ) - 1 ) )
}

# locates the min/max lat/lons around the geo_dna
sub bounding_box {
    my ( $geodna ) = @_;

    my @chars = split( //, $geodna );

    my $loni;
    my $lati = [ -90.0, 90.0 ];
    my $first = shift @chars;

    if ( $first eq 'w' ) {
        $loni = [ -180.0, 0.0 ];
    } elsif ( $first eq 'e' ) {
        $loni = [ 0.0, 180.0 ];
    }

    foreach my $c (@chars) {
        my $cd = $DECODE_MAP->{$c};
        if ( $cd & 2 ) {
            $loni = [ ( $loni->[0] + $loni->[1] ) / 2, $loni->[1] ];
        } else {
            $loni = [ $loni->[0],  ( $loni->[0] + $loni->[1] ) / 2 ];
        }
        if ( $cd & 1 ) {
            $lati = [ ( $lati->[0] + $lati->[1] ) / 2, $lati->[1] ];
        } else {
            $lati = [ $lati->[0],  ( $lati->[0] + $lati->[1] ) / 2 ];
        }
    }
    return ( $lati, $loni );
}

sub add_vector {
    my ( $lat, $lon, $dy, $dx ) = @_;

    return (
        fmod( ( $lat + 90.0 + $dy ), 180.0 ) - 90.0,
        fmod( ( $lon + 180.0 + $dx ), 360.0 )  - 180.0
    );
}

# """
# Return the eight neighboring geodna codes

sub neighbours_geo_dna {
    my ( @args ) = @_;
    neighbours( @args );
}

sub neighbours {
    my ( $geodna ) = @_;

    # TODO:kd - this can be optimised
    my ( $lat, $lon )   = decode( $geodna );
    my ( $lati, $loni ) = bounding_box( $geodna );
    my $width  = abs( $loni->[1] - $loni->[0] );
    my $height = abs( $lati->[1] - $lati->[0] );

    my $neighbours = [];
    foreach my $y ( -1, 0, 1 ) {
        foreach my $x ( -1, 0, 1 ) {
            next unless ( $x || $y );
            push (@$neighbours, encode( add_vector( $lat, $lon, $height * $y, $width * $x ) ) );
        }
    }
    return $neighbours;
}

my $PYTHON = <<PYTHON

PYTHON
;

1;