package Geo::DNA;

use common::sense;

use Math::Trig qw( :pi rad2deg asin deg2rad );
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

my $RADIALS = {
    '0_N'  => 0,
    '1_NE' => pip4,
    '2_E'  => pip2,
    '3_SE' => pip2 + pip4,
    '4_S'  => pi,
    '5_SW' => pi + pip4,
    '6_W'  => pi + pip2,
    '7_NW' => pi + pip2 + pip4,
};

# """
# Return the eight neighboring geoprints and optionally their
# compass bearing flag around the given geoprint.
#
# If bearing is True, bearing flag will be either None (not
# adjacent) or one of ['N', 'NE', 'E', 'SE', 'S', 'SW', 'W', 'NW']
# """

sub neighbours_geo_dna {
    my ( @args ) = @_;
    neighbours( @args );
}

sub neighbours {
    my ( $geo_dna, @opts ) = @_;
    my $options = { @opts };

    my $results = [];
    my $precision = length( $geo_dna );
    my ( $lat1, $lon1 ) = decode( $geo_dna, radians => 1 );
    #print STDERR " $lat1 / $lon1 \n";
    my $d = size( $geo_dna ) / $RADIUS_OF_EARTH;
    foreach my $b ( sort keys %$RADIALS ) {
        my $r = $RADIALS->{$b};
        my $lat = asin(sin($lat1) * cos($d) +
                       cos($lat1) * sin($d) * cos($r));
        my $dlon = atan2(sin($r) * sin($d) * cos($lat1),
                         cos($d) - sin($lat1) * sin($lat));
        my $lon = ( $lon1 + $dlon ) % pi2 - pi;
        my $h = encode( $lat, $lon, precision => $precision, radians => 1 );
        if ( $options->{bearing} ) {
            push @$results, [ $b, $h ];
        } else {
            push @$results, $h;
        }
    }
    return $results;
}

my $PYTHON = <<PYTHON

def neighbors(geoprint, bearing=True):

    results = set()
    precision = len(geoprint)
    lat1, lon1 = decode(geoprint, radians=True)
    d = size(geoprint) / 6378100 # radius of earth in meters
    for b, r in _radials.iteritems():
        lat = asin(sin(lat1) * cos(d) +
                   cos(lat1) * sin(d) * cos(r))
        dlon = atan2(sin(r) * sin(d) * cos(lat1),
                     cos(d) - sin(lat1) * sin(lat))
        lon = mod(lon1 + dlon + pi, 2 * pi) - pi
        h = encode(lat, lon, precision=precision, radians=True)
        if bearing:
            h = (b, h)
        results.add(h)
    return results


def adjacent(first, second):
    """
    Return a flag indicating if two geoprints are adjacent to each
    other.
    """
    if min(len(first), len(second)) < 3:
        raise TypeError(
            "Adjacency requires at least 3 characters of precision.")
    maxe = (error(first, True) / 2) + (error(second, True) / 2)
    maxe2 = pow(maxe, 2)
    h = sqrt(maxe2 + maxe2)
    d = distance(first, second, True)
    return h <= d

PYTHON
;

1;