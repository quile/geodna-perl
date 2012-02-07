"""Geoprint is a geohash variant.

This algorithm is almost exactly identical to, and originally based
on, the Geohash module on Pypi by Leonard Norrgard, which in turn,
implements the geohash algorithm developed by Gustavo Niemeyer and
described in the wikipedia article:

http://en.wikipedia.org/wiki/Geohash

The main difference between this encoding and geohash is that it uses
a much smaller 'alphabet' of encoding characters to spell a hash,
instead of base32 encoding it uses simpler base4 encoding that uses
four possible states to encode each character instead of 32.  This
results in longer, uglier, less compressed geohashes but they are more
useful for doing proximity searching using hash prefixes.  Note that
like geohashes, there are edge cases near the equator and meridians
where points near each other have different prefixes.

All geoprints begin with either a 'w' or 'e' character to specify that
the hash belongs in the western or eastern hemisphere.  The western
hemisphere includes all longitudes less than zero down to -180.0, and
the eastern hemisphere includes all longitudes greater than or equal
to zero, up to 180.0.  Subsequent hash characters are one of g, a, t,
or c.
"""

from collections import namedtuple
from operator import neg, mod
from math import radians as rads
from math import (
    asin,
    atan2,
    cos,
    degrees,
    log10,
    pi,
    sin,
    sqrt,
    )

interval = namedtuple("interval", "min max")

alphabet = 'gatc'
decodemap = dict((k, i) for (i, k) in enumerate(alphabet))


def encode(latitude, longitude, precision=22, radians=False):
    """Encode the given latitude and longitude into a geoprint that
    contains 'precision' characters.

    If radians is True, input parameters are in radians, otherwise
    they are degrees.  Example::

    >>> c = (7.0625, -95.677068)
    >>> h = encode(*c)
    >>> h
    'watttatcttttgctacgaagt'

    >>> r = rads(c[0]), rads(c[1])
    >>> h2 = encode(*r, radians=True)
    >>> h == h2
    True
    """
    if radians:
        latitude = degrees(latitude)
        longitude = degrees(longitude)
    if longitude < 0:
        geoprint = ['w']
        loni = interval(-180.0, 0.0)
    else:
        geoprint = ['e']
        loni = interval(0.0, 180.0)

    lati = interval(-90.0, 90.0)

    while len(geoprint) < precision:
        ch = 0
        mid = (loni.min + loni.max) / 2
        if longitude > mid:
            ch |= 2
            loni = interval(mid, loni.max)
        else:
            loni = interval(loni.min, mid)

        mid = (lati.min + lati.max) / 2
        if latitude > mid:
            ch |= 1
            lati = interval(mid, lati.max)
        else:
            lati = interval(lati.min, mid)

        geoprint.append(alphabet[ch])
    return ''.join(geoprint)


def decode(geoprint, radians=False):
    """Decode a geoprint, returning the latitude and longitude.  These
    coordinates should approximate the input coordinates within a
    degree of error returned by 'error()'

    >>> c = (7.0625, -95.677068)
    >>> h = encode(*c)
    >>> c2 = decode(h)
    >>> e = error(h)
    >>> abs(c[0] - c2[0]) <= e
    True
    >>> abs(c[1] - c2[1]) <= e
    True

    If radians is True, results are in radians instead of degrees.

    >>> c2 = decode(h, radians=True)
    >>> e = error(h, radians=True)
    >>> abs(rads(c[0]) - c2[0]) <= e
    True
    >>> abs(rads(c[1]) - c2[1]) <= e
    True
    """
    lati = interval(-90.0, 90.0)
    first = geoprint[0]
    if first == 'w':
        loni = interval(-180.0, 0.0)
    elif first == 'e':
        loni = interval(0.0, 180.0)

    geoprint = geoprint[1:]

    for c in geoprint:
        cd = decodemap[c]
        if cd & 2:
            loni = interval((loni.min + loni.max) / 2, loni.max)
        else:
            loni = interval(loni.min, (loni.min + loni.max) / 2)
        if cd & 1:
            lati = interval((lati.min + lati.max) / 2, lati.max)
        else:
            lati = interval(lati.min, (lati.min + lati.max) / 2)
    lat = (lati.min + lati.max) / 2
    lon = (loni.min + loni.max) / 2
    if radians:
        return rads(lat), rads(lon)
    return lat, lon


def error(geoprint, radians=False):
    """Returns the error of a given geoprint.

    If radians is true, return the error in radians, otherwise
    degrees.
    """
    e = 90.0 * pow(2, neg(len(geoprint) - 1))
    if radians:
        return rads(e)
    return e


def size(geoprint):
    """Return the *approximate* size in meters of one side of a
    geoprint box.
    """
    return 20000000.0 * pow(2, neg(len(geoprint) - 1))


def distance(start, end, radians=False):
    """
    Calculate the *approximate* distance between two geoprints.  Based
    on Haversine formula
    (http://en.wikipedia.org/wiki/Haversine_formula).

    If radians is True, return distance is radians on a approximate
    great circle, otherwise return meters.
    """
    start_lat, start_lon = decode(start, radians=True)
    end_lat, end_lon = decode(end, radians=True)
    d_lat = end_lat - start_lat
    d_long = end_lon - start_lon
    a = (sin(d_lat / 2) ** 2 + cos(start_lat) *
         cos(end_lat) * sin(d_long / 2) ** 2)
    c = 2 * atan2(sqrt(a), sqrt(1 - a))
    if radians:
        return c
    return 6378100 * c


def bearing(start, end, radians=False):
    """
    Calculate *approximate* bearing between two geoprints.

    If radians is True, returns bearing in radians, otherwise degrees.
    """
    start_lat, start_lon = decode(start, radians=True)
    end_lat, end_lon = decode(end, radians=True)
    d_lon = end_lon - start_lon
    c = atan2(sin(d_lon) * cos(end_lat),
              cos(start_lat) * sin(end_lat) -
              sin(start_lat) * cos(end_lat) * cos(d_lon))
    if radians:
        return c
    return degrees(c)


def format(geoprint):
    """Return two formatted strings of latitude and longitude,
    truncated to the number of decimal places known based on the
    computed error in the geoprint.
    """
    lat, lon = decode(geoprint)
    err = max(1, int(round(-log10(error(geoprint))))) - 1
    lats = "%.*f" % (err, lat)
    lons = "%.*f" % (err, lon)
    if '.' in lats:
        lats = lats.rstrip('0')
    if '.' in lons:
        lons = lons.rstrip('0')
    return lats, lons


_radials = {'N' : 0,              # N
            'NE': pi / 4,         # NE
            'E' : pi / 2,         # E
            'SE': (3 * pi / 4),   # SE
            'S' : pi,             # S
            'SW': (5 * pi / 4),   # SW
            'W' : (3 * pi / 2),   # W
            'NW': (7 * pi / 4),   # NW
            }


def neighbors(geoprint, bearing=True):
    """
    Return the eight neighboring geoprints and optionally their
    compass bearing flag around the given geoprint.

    If bearing is True, bearing flag will be either None (not
    adjacent) or one of ['N', 'NE', 'E', 'SE', 'S', 'SW', 'W', 'NW']
    """
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
