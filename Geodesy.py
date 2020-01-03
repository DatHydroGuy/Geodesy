import math


class Angle:
    pi_over_180 = math.pi / 180.0
    _degrees = 0.0
    zero = 0.0
    angle180 = 180.0
    invalid_values = [math.inf, -math.inf, math.nan]
    _tolerance = 0.00000001

    def assign_args(self, args):
        minutes = seconds = 0.0
        if len(args) > 2:
            seconds = args[2]
        if len(args) > 1:
            minutes = args[1]
        degrees = args[0]
        try:
            if degrees in self.invalid_values or minutes in self.invalid_values or seconds in self.invalid_values:
                raise ValueError('An Angle cannot have infinite or NaN arguments')
            if degrees != 0 and (minutes < 0 or seconds < 0):
                raise ValueError
            if minutes != 0 and seconds < 0:
                raise ValueError
            self._degrees = (seconds / 3600.0) + (minutes / 60.0)
            self._degrees = degrees - self._degrees if degrees < 0.0 else degrees + self._degrees
        except (TypeError, ValueError):
            self._degrees = 0.0
        return

    def assign_kwargs(self, kwargs):
        degrees = kwargs['degrees'] if 'degrees' in kwargs else 0.0
        minutes = kwargs['minutes'] if 'minutes' in kwargs else 0.0
        seconds = kwargs['seconds'] if 'seconds' in kwargs else 0.0
        try:
            if degrees in self.invalid_values or minutes in self.invalid_values or seconds in self.invalid_values:
                raise ValueError('An Angle cannot have infinite or NaN arguments')
            if degrees != 0 and (minutes < 0 or seconds < 0):
                raise TypeError
            if minutes != 0 and seconds < 0:
                raise TypeError
            self._degrees = (seconds / 3600.0) + (minutes / 60.0)
            self._degrees = degrees - self._degrees if degrees < 0.0 else degrees + self._degrees
        except (TypeError, KeyError, ValueError):
            self._degrees = 0.0
        return

    def __init__(self, *args, **kwargs):
        if 0 < len(args) < 4:
            self.assign_args(args)
        if 0 < len(kwargs) < 4:
            self.assign_kwargs(kwargs)

    def __add__(self, other):
        if type(other) == Angle:
            return Angle(self.degrees + other.degrees)
        elif type(other) == int or type(other) == float and not math.isnan(other) and not math.isinf(other):
            return Angle(self.degrees + other)
        else:
            return self

    def __sub__(self, other):
        if type(other) == Angle:
            return Angle(self.degrees - other.degrees)
        elif type(other) == int or type(other) == float and not math.isnan(other) and not math.isinf(other):
            return Angle(self.degrees - other)
        else:
            return self

    def __abs__(self):
        return Angle(abs(self.degrees))

    def __lt__(self, other):
        if type(other) == Angle:
            return self.degrees < other.degrees
        elif type(other) == int or type(other) == float and not math.isnan(other) and not math.isinf(other):
            return self.degrees < other
        else:
            return False

    def __le__(self, other):
        if type(other) == Angle:
            return self.degrees <= other.degrees
        elif type(other) == int or type(other) == float and not math.isnan(other) and not math.isinf(other):
            return self.degrees <= other
        else:
            return False

    def __gt__(self, other):
        if type(other) == Angle:
            return self.degrees > other.degrees
        elif type(other) == int or type(other) == float and not math.isnan(other) and not math.isinf(other):
            return self.degrees > other
        else:
            return False

    def __ge__(self, other):
        if type(other) == Angle:
            return self.degrees >= other.degrees
        elif type(other) == int or type(other) == float and not math.isnan(other) and not math.isinf(other):
            return self.degrees >= other
        else:
            return False

    def __eq__(self, other):
        if type(other) == Angle:
            return math.fabs(self.degrees - other.degrees) < self._tolerance
        elif type(other) == int or type(other) == float and not math.isnan(other) and not math.isinf(other):
            return math.fabs(self.degrees - other) < self._tolerance
        else:
            return False

    def __ne__(self, other):
        if type(other) == Angle:
            return math.fabs(self.degrees - other.degrees) >= self._tolerance
        elif type(other) == int or type(other) == float and not math.isnan(other) and not math.isinf(other):
            return math.fabs(self.degrees - other) >= self._tolerance
        else:
            return False

    def __hash__(self):
        return int(self.degrees * 1000033)

    def to_radians(self, degrees):
        return degrees * self.pi_over_180

    def to_degrees(self, radians):
        return radians / self.pi_over_180

    @property
    def degrees(self):
        return self._degrees

    @degrees.setter
    def degrees(self, value):
        if type(value) == int or type(value) == float and not math.isinf(value) and not math.isnan(value):
            self._degrees = value

    @property
    def radians(self):
        return self._degrees * self.pi_over_180

    @radians.setter
    def radians(self, value):
        if type(value) == int or type(value) == float and not math.isinf(value) and not math.isnan(value):
            self._degrees = value / self.pi_over_180


class Ellipsoid:
    _semi_major_axis = None
    _semi_minor_axis = None
    _flattening = None
    _inverse_flattening = None

    def assign_args(self, args):
        if any(type(arg) != float and type(arg) != int for arg in args):
            return

        if len(args) >= 4:
            self._inverse_flattening = args[3]
        if len(args) >= 3:
            self._flattening = args[2]
        if len(args) >= 2:
            self._semi_minor_axis = args[1]
        if len(args) >= 1:
            self._semi_major_axis = args[0]
        return

    def assign_kwargs(self, kwargs):
        if any(type(arg) != float and type(arg) != int for arg in kwargs.values()):
            return
        self._semi_major_axis = kwargs['semi_major_axis'] if 'semi_major_axis' in kwargs else 0.0
        self._semi_minor_axis = kwargs['semi_minor_axis'] if 'semi_minor_axis' in kwargs else 0.0
        self._flattening = kwargs['flattening'] if 'flattening' in kwargs else 0.0
        self._inverse_flattening = kwargs['inverse_flattening'] if 'inverse_flattening' in kwargs else 0.0
        return

    def __init__(self, *args, **kwargs):
        if 0 < len(args) < 5:
            self.assign_args(args)
        if 0 < len(kwargs) < 5:
            self.assign_kwargs(kwargs)

    @staticmethod
    def from_semi_major_and_inverse_flattening(semi_major, inverse_flattening):
        if (type(semi_major) != float and type(semi_major) != int) or \
                (type(inverse_flattening) != float and type(inverse_flattening) != int):
            return Ellipsoid()
        near_zero = 0.00000000000001
        infinite = math.inf if inverse_flattening >= 0 else -math.inf
        flattening = infinite if abs(inverse_flattening) < near_zero else 1.0 / inverse_flattening
        semi_minor = (1.0 - flattening) * semi_major
        return Ellipsoid(semi_major, semi_minor, flattening, inverse_flattening)

    @staticmethod
    def from_semi_major_and_flattening(semi_major, flattening):
        if (type(semi_major) != float and type(semi_major) != int) or \
                (type(flattening) != float and type(flattening) != int):
            return Ellipsoid()
        near_zero = 0.00000000000001
        infinite = math.inf if flattening >= 0 else -math.inf
        inverse_flattening = infinite if abs(flattening) < near_zero else 1.0 / flattening
        semi_minor = (1.0 - flattening) * semi_major
        return Ellipsoid(semi_major, semi_minor, flattening, inverse_flattening)

    @staticmethod
    def wgs84():
        return Ellipsoid.from_semi_major_and_inverse_flattening(6378137.0, 298.257223563)

    @staticmethod
    def grs80():
        return Ellipsoid.from_semi_major_and_inverse_flattening(6378137.0, 298.257222101)

    @staticmethod
    def grs67():
        return Ellipsoid.from_semi_major_and_inverse_flattening(6378160.0, 298.25)

    @staticmethod
    def ans():
        return Ellipsoid.from_semi_major_and_inverse_flattening(6378160.0, 298.25)

    @staticmethod
    def wgs72():
        return Ellipsoid.from_semi_major_and_inverse_flattening(6378135.0, 298.26)

    @staticmethod
    def clarke1858():
        return Ellipsoid.from_semi_major_and_inverse_flattening(6378293.645, 294.26)

    @staticmethod
    def clarke1880():
        return Ellipsoid.from_semi_major_and_inverse_flattening(6378249.145, 293.465)

    @staticmethod
    def sphere():
        return Ellipsoid.from_semi_major_and_flattening(6371000, 0)

    @property
    def semi_major_axis(self):
        return self._semi_major_axis

    @property
    def semi_minor_axis(self):
        return self._semi_minor_axis

    @property
    def flattening(self):
        return self._flattening

    @property
    def inverse_flattening(self):
        return self._inverse_flattening


class GeodeticCurve:
    _ellipsoidal_distance = None
    _azimuth = Angle()
    _reverse_azimuth = Angle()

    def assign_args(self, args):
        if any(type(arg) != float and type(arg) != int and type(arg) != Angle for arg in args):
            raise TypeError('Ellipsoidal distance must be numerical.\n'
                            'Azimuth and reverse_azimuth must be numerical or valid Angle objects.')

        if len(args) > 2:
            if type(args[2]) == Angle:
                self._reverse_azimuth = args[2]
            else:
                self._reverse_azimuth = Angle(args[2])
        if len(args) > 1:
            if type(args[1]) == Angle:
                self._azimuth = args[1]
            else:
                self._azimuth = Angle(args[1])
        if len(args) > 0:
            self._ellipsoidal_distance = args[0]
        return

    def assign_kwargs(self, kwargs):
        if any(type(arg) != float and type(arg) != int and type(arg) != Angle for arg in kwargs.values()):
            raise TypeError('Ellipsoidal distance must be numerical.\n'
                            'Azimuth and reverse_azimuth must be numerical or valid Angle objects.')

        self._ellipsoidal_distance = kwargs['ellipsoidal_distance'] if 'ellipsoidal_distance' in kwargs else 0.0

        if 'azimuth' in kwargs:
            if type(kwargs['azimuth']) == Angle:
                self._azimuth = kwargs['azimuth']
            else:
                self._azimuth = Angle(kwargs['azimuth'])
        else:
            self._azimuth = Angle()

        if 'reverse_azimuth' in kwargs:
            if type(kwargs['reverse_azimuth']) == Angle:
                self._reverse_azimuth = kwargs['reverse_azimuth']
            else:
                self._reverse_azimuth = Angle(kwargs['reverse_azimuth'])
        else:
            self._reverse_azimuth = Angle()

        return

    def __init__(self, *args, **kwargs):
        if 0 < len(args) < 4:
            self.assign_args(args)
        if 0 < len(kwargs) < 4:
            self.assign_kwargs(kwargs)

    def __str__(self):
        return 's={0}; a12={1}; a21={2};'.format(self._ellipsoidal_distance,
                                                 self._azimuth.degrees,
                                                 self._reverse_azimuth.degrees)

    @property
    def ellipsoidal_distance(self):
        return self._ellipsoidal_distance

    @property
    def azimuth(self):
        return self._azimuth

    @property
    def reverse_azimuth(self):
        return self._reverse_azimuth


class GeodeticMeasurement:
    _curve = GeodeticCurve()
    _elevation_change = None
    _point_to_point = None

    def assign_args(self, args):
        if any(type(arg) != float and type(arg) != int and type(arg) != GeodeticCurve for arg in args):
            raise TypeError('Curve must be a valid GeodeticCurve object, and elevation_change must be numeric.')

        if len(args) > 1:
            if type(args[1]) == float or type(args[1]) == int:
                self._elevation_change = args[1]
            else:
                return False

        if len(args) > 0:
            if type(args[0]) == GeodeticCurve:
                self._curve = args[0]
            else:
                return False

        return True

    def assign_kwargs(self, kwargs):
        if any(type(arg) != float and type(arg) != int and type(arg) != GeodeticCurve for arg in kwargs.values()):
            raise TypeError('Curve must be a valid GeodeticCurve object, and elevation_change must be numeric.')

        if 'curve' in kwargs:
            if type(kwargs['curve']) == GeodeticCurve:
                self._curve = kwargs['curve']
            else:
                return False

        if 'elevation_change' in kwargs:
            if type(kwargs['elevation_change']) == float or type(kwargs['elevation_change']) == int:
                self._elevation_change = kwargs['elevation_change']
            else:
                return False

        return True

    def __init__(self, *args, **kwargs):
        assigned_args = False
        if 0 < len(args) < 3:
            assigned_args = self.assign_args(args)
        if 0 < len(kwargs) < 3:
            assigned_args = self.assign_kwargs(kwargs)

        if assigned_args:
            ellipsoidal_distance = self.curve.ellipsoidal_distance
            try:
                self._point_to_point = math.sqrt(ellipsoidal_distance ** 2 + self.elevation_change ** 2)
            except TypeError:
                self._point_to_point = None

    @property
    def curve(self):
        return self._curve

    @property
    def ellipsoidal_distance(self):
        return self.curve.ellipsoidal_distance

    @property
    def azimuth(self):
        return self.curve.azimuth

    @property
    def reverse_azimuth(self):
        return self.curve.reverse_azimuth

    @property
    def elevation_change(self):
        return self._elevation_change

    @property
    def point_to_point(self):
        return self._point_to_point


class GlobalCoordinates:
    _latitude = Angle()
    _longitude = Angle()

    def assign_args(self, args):
        if any(type(arg) != float and type(arg) != int and type(arg) != Angle for arg in args):
            raise TypeError('Latitude and longitude must be numerical or valid Angle objects')

        if len(args) == 2:
            if type(args[1]) == Angle:
                self._longitude = args[1]
            elif math.isnan(args[1]) or math.isinf(args[1]):
                raise TypeError('Longitude must be numerical or a valid Angle object')
            else:
                self._longitude = Angle(args[1])
            if type(args[0]) == Angle:
                self._latitude = args[0]
            elif math.isnan(args[0]) or math.isinf(args[0]):
                raise TypeError('Latitude must be numerical or a valid Angle object')
            else:
                self._latitude = Angle(args[0])
        return True

    def assign_kwargs(self, kwargs):
        if any(type(arg) != float and type(arg) != int and type(arg) != Angle for arg in kwargs.values()):
            raise TypeError('Latitude and longitude must be numerical or valid Angle objects')

        if 'latitude' in kwargs:
            if type(kwargs['latitude']) == Angle:
                self._latitude = kwargs['latitude']
            elif math.isnan(kwargs['latitude']) or math.isinf(kwargs['latitude']):
                raise TypeError('Latitude must be numerical or a valid Angle object')
            else:
                self._latitude = Angle(kwargs['latitude'])
        else:
            self._latitude = Angle()

        if 'longitude' in kwargs:
            if type(kwargs['longitude']) == Angle:
                self._longitude = kwargs['longitude']
            elif math.isnan(kwargs['longitude']) or math.isinf(kwargs['longitude']):
                raise TypeError('Longitude must be numerical or a valid Angle object')
            else:
                self._longitude = Angle(kwargs['longitude'])
        else:
            self._longitude = Angle()

        return True

    def canonicalize(self):
        latitude = self._latitude.degrees
        longitude = self._longitude.degrees

        latitude = (latitude + 180) % 360
        if latitude < 0:
            latitude += 360
        latitude -= 180

        if latitude > 90:
            latitude = 180 - latitude
            longitude += 180
        elif latitude < -90:
            latitude = -180 - latitude
            longitude += 180

        longitude = (longitude + 180) % 360
        if longitude <= 0:
            longitude += 360
        longitude -= 180

        self._latitude = Angle(latitude)
        self._longitude = Angle(longitude)

    def __init__(self, *args, **kwargs):
        assigned_args = False
        if len(args) == 2:
            assigned_args = self.assign_args(args)
        if len(kwargs) == 2:
            assigned_args = self.assign_kwargs(kwargs)

        if assigned_args:
            self.canonicalize()
        else:
            self._latitude = Angle()
            self._longitude = Angle()

    def __hash__(self):
        return (hash(self.longitude) * (hash(self.latitude) + 1021)) * 1000033

    def __eq__(self, other):
        if type(other) is not GlobalCoordinates:
            return False

        # return self.longitude == other.longitude and self.latitude == other.latitude
        return_value = 0
        if self.longitude < other.longitude:
            return_value -= 1
        elif self.longitude > other.longitude:
            return_value += 1
        elif self.latitude < other.latitude:
            return_value -= 1
        elif self.latitude > other.latitude:
            return_value += 1

        return return_value == 0

    def __ne__(self, other):
        return not self == other

    def __lt__(self, other):
        if self.longitude == other.longitude:
            return self.latitude < other.latitude

        return self.longitude < other.longitude

    def __gt__(self, other):
        if self.longitude == other.longitude:
            return self.latitude > other.latitude

        return self.longitude > other.longitude

    def __le__(self, other):
        return not self > other

    def __ge__(self, other):
        return not self < other

    def __str__(self):
        if self.latitude.degrees is not None and self.longitude.degrees is not None:
            north_south = 'N' if self.latitude.degrees >= 0 else 'S'
            east_west = 'E' if self.longitude.degrees >= 0 else 'W'
            return '{0}{1}; {2}{3};'.format(abs(self.latitude.degrees), north_south,
                                            abs(self.longitude.degrees), east_west)
        return 'Invalid GlobalCoordinates'

    @property
    def latitude(self):
        return self._latitude

    @latitude.setter
    def latitude(self, value):
        if type(value) == Angle:
            self._latitude = value
            self.canonicalize()
        elif type(value) == int or type(value) == float and not math.isinf(value) and not math.isnan(value):
            self._latitude = Angle(value)
            self.canonicalize()

    @property
    def longitude(self):
        return self._longitude

    @longitude.setter
    def longitude(self, value):
        if type(value) == Angle:
            self._longitude = value
            self.canonicalize()
        elif type(value) == int or type(value) == float and not math.isinf(value) and not math.isnan(value):
            self._longitude = Angle(value)
            self.canonicalize()


class GlobalPosition:
    _coordinates = GlobalCoordinates()
    _elevation = 0.0

    def assign_args(self, args):
        if any(type(arg) != float and type(arg) != int and type(arg) != GlobalCoordinates for arg in args):
            raise TypeError('Elevation must be numerical, and coordinates must be a valid GlobalCoordinates object')

        if len(args) > 1:
            if (type(args[1]) == float or type(args[1]) == int) and not math.isinf(args[1]) and not math.isnan(args[1]):
                self._elevation = args[1]
            else:
                raise TypeError('Elevation must be numerical')

        if len(args) > 0:
            if type(args[0]) == GlobalCoordinates:
                self._coordinates = args[0]
            else:
                raise TypeError('Coordinates must be a valid GlobalCoordinates object')

    def assign_kwargs(self, kwargs):
        if 'coordinates' in kwargs and type(kwargs['coordinates']) == GlobalCoordinates:
            self._coordinates = kwargs['coordinates']
        else:
            raise TypeError('Coordinates must be a valid GlobalCoordinates object')

        if 'elevation' in kwargs and (type(kwargs['elevation']) == int or type(kwargs['elevation']) == float):
            self._elevation = kwargs['elevation']
            if math.isinf(self._elevation) or math.isnan(self._elevation):
                raise TypeError('Elevation must be numerical')

    def __init__(self, *args, **kwargs):
        if 0 < len(args) < 3:
            self.assign_args(args)
        if 0 < len(kwargs) < 3:
            self.assign_kwargs(kwargs)

    def __hash__(self):
        return hash(self._coordinates) if self._elevation == 0 else hash(self._coordinates) * int(self._elevation)

    def __eq__(self, other):
        if type(other) is not GlobalPosition:
            return False

        return self._elevation == other.elevation and self._coordinates == other.coordinates

    def __ne__(self, other):
        return not self == other

    def __lt__(self, other):
        if self._coordinates == other.coordinates:
            return self._elevation < other.elevation

        return self._coordinates < other.coordinates

    def __gt__(self, other):
        if self._coordinates == other.coordinates:
            return self._elevation > other.elevation

        return self._coordinates > other.coordinates

    def __le__(self, other):
        return not self > other

    def __ge__(self, other):
        return not self < other

    def __str__(self):
        if self._elevation is not None and self._coordinates is not None:
            return '{0} elevation={1}m'.format(self._coordinates, self._elevation)
        return 'Invalid GlobalPosition'

    @property
    def coordinates(self):
        return self._coordinates

    @coordinates.setter
    def coordinates(self, value):
        if type(value) == GlobalCoordinates:
            self._coordinates = value

    @property
    def latitude(self):
        return Angle(self._coordinates.latitude.degrees)

    @latitude.setter
    def latitude(self, value):
        if type(value) == Angle:
            self._coordinates.latitude = value
        elif type(value) == int or type(value) == float and not math.isinf(value) and not math.isnan(value):
            self._coordinates.latitude = Angle(value)

    @property
    def longitude(self):
        return Angle(self._coordinates.longitude.degrees)

    @longitude.setter
    def longitude(self, value):
        if type(value) == Angle:
            self._coordinates.longitude = value
        elif type(value) == int or type(value) == float and not math.isinf(value) and not math.isnan(value):
            self._coordinates.longitude = Angle(value)

    @property
    def elevation(self):
        return self._elevation

    @elevation.setter
    def elevation(self, value):
        if type(value) == int or type(value) == float and not math.isinf(value) and not math.isnan(value):
            self._elevation = value


class GeodeticCalculator:
    _two_pi = 2.0 * math.pi
    _tolerance = 0.0000000000001

    def calculate_ending_global_coordinates(self, ellipsoid, start_coordinates, start_bearing, distance):
        a = ellipsoid.semi_major_axis
        b = ellipsoid.semi_minor_axis
        a_squared = a * a
        b_squared = b * b
        f = ellipsoid.flattening
        phi1 = start_coordinates.latitude.radians
        alpha1 = start_bearing.radians
        cos_alpha1 = math.cos(alpha1)
        sin_alpha1 = math.sin(alpha1)
        s = distance
        tan_u1 = (1.0 - f) * math.tan(phi1)
        cos_u1 = 1.0 / math.sqrt(1.0 + tan_u1 * tan_u1)
        sin_u1 = tan_u1 * cos_u1

        # Equation 1
        sigma1 = math.atan2(tan_u1, cos_alpha1)

        # Equation 2
        sin_alpha = cos_u1 * sin_alpha1
        sin2_alpha = sin_alpha * sin_alpha
        cos2_alpha = 1.0 - sin2_alpha
        u_squared = cos2_alpha * (a_squared - b_squared) / b_squared

        # Equation 3
        big_a = 1 + (u_squared / 16384) * (4096 + u_squared * (-768 + u_squared * (320 - 175 * u_squared)))

        # Equation 4
        big_b = (u_squared / 1024) * (256 + u_squared * (-128 + u_squared * (74 - 47 * u_squared)))

        # Iterate until sigma converges to self._tolerance
        s_over_b_a = s / (b * big_a)
        sigma = s_over_b_a
        prev_sigma = s_over_b_a

        while True:
            # Equation 5
            sigma_m2 = 2.0 * sigma1 + sigma
            cos_sigma_m2 = math.cos(sigma_m2)
            cos2_sigma_m2 = cos_sigma_m2 * cos_sigma_m2
            sin_sigma = math.sin(sigma)
            cos_sigma = math.cos(sigma)

            # Equation 6
            delta_sigma = big_b * sin_sigma * (cos_sigma_m2 + (big_b / 4.0) *
                                               (cos_sigma * (-1.0 + 2.0 * cos2_sigma_m2) - (big_b / 6.0) *
                                                cos_sigma_m2 * (-3.0 + 4.0 * sin_sigma * sin_sigma) *
                                                (-3.0 + 4.0 * cos2_sigma_m2)))

            # Equation 7
            sigma = s_over_b_a + delta_sigma

            if math.fabs(sigma - prev_sigma) < self._tolerance:
                break

            prev_sigma = sigma

        sigma_m2 = 2.0 * sigma1 + sigma
        cos_sigma_m2 = math.cos(sigma_m2)
        cos2_sigma_m2 = cos_sigma_m2 * cos2_sigma_m2

        cos_sigma = math.cos(sigma)
        sin_sigma = math.sin(sigma)

        # Equation 8
        phi2 = math.atan2(sin_u1 * cos_sigma + cos_u1 * sin_sigma * cos_alpha1,
                          (1.0 - f) * math.sqrt(sin2_alpha +
                                                math.pow(sin_u1 * sin_sigma - cos_u1 * cos_sigma * cos_alpha1, 2.0)))

        # Equation 9
        lamda = math.atan2(sin_sigma * sin_alpha1, cos_u1 * cos_sigma - sin_u1 * sin_sigma * cos_alpha1)

        # Equation 10
        big_c = (f / 16.0) * cos2_alpha * (4.0 + f * (4.0 - 3.0 * cos2_alpha))

        # Equation 11
        big_l = lamda - (1.0 - big_c) * f * sin_alpha *\
                        (sigma + big_c * sin_sigma *
                         (cos_sigma_m2 + big_c * cos_sigma * (-1.0 + 2.0 * cos2_sigma_m2)))

        # Equation 12
        alpha2 = math.atan2(sin_alpha, -sin_u1 * sin_sigma + cos_u1 * cos_sigma * cos_alpha1)

        latitude = Angle()
        longitude = Angle()

        latitude.radians = phi2
        longitude.radians = start_coordinates.longitude.radians + big_l

        end_bearing = Angle()
        end_bearing.radians = alpha2

        end_coordinates = GlobalCoordinates(latitude=latitude, longitude=longitude)

        return end_coordinates, end_bearing

    def calculate_geodetic_curve(self, ellipsoid, start_coordinates, end_coordinates):
        a = ellipsoid.semi_major_axis
        b = ellipsoid.semi_minor_axis
        f = ellipsoid.flattening

        phi1 = start_coordinates.latitude.radians
        lamda1 = start_coordinates.longitude.radians
        phi2 = end_coordinates.latitude.radians
        lamda2 = end_coordinates.longitude.radians

        a2 = a * a
        b2 = b * b
        a2b2b2 = (a2 - b2) / b2

        omega = lamda2 - lamda1

        tan_phi1 = math.tan(phi1)
        tan_u1 = (1.0 - f) * tan_phi1
        u1 = math.atan(tan_u1)
        sin_u1 = math.sin(u1)
        cos_u1 = math.cos(u1)

        tan_phi2 = math.tan(phi2)
        tan_u2 = (1.0 - f) * tan_phi2
        u2 = math.atan(tan_u2)
        sin_u2 = math.sin(u2)
        cos_u2 = math.cos(u2)

        sin_u1_sin_u2 = sin_u1 * sin_u2
        cos_u1_sin_u2 = cos_u1 * sin_u2
        sin_u1_cos_u2 = sin_u1 * cos_u2
        cos_u1_cos_u2 = cos_u1 * cos_u2

        # Equation 13
        lamda = omega

        big_a = 0.0
        sigma = 0.0
        delta_sigma = 0.0
        converged = False

        for i in range(20):
            lamda0 = lamda

            sin_lamda = math.sin(lamda)
            cos_lamda = math.cos(lamda)

            # Equation 14
            sin2_sigma = (cos_u2 * sin_lamda * cos_u2 * sin_lamda) + math.pow(cos_u1_sin_u2 - sin_u1_cos_u2 * cos_lamda,
                                                                              2.0)
            sin_sigma = math.sqrt(sin2_sigma)

            # Equation 15
            cos_sigma = sin_u1_sin_u2 + (cos_u1_cos_u2 * cos_lamda)

            # Equation 16
            sigma = math.atan2(sin_sigma, cos_sigma)

            # Equation 17
            sin_alpha = 0.0 if sin2_sigma == 0.0 else cos_u1_cos_u2 * sin_lamda / sin_sigma
            alpha = math.asin(sin_alpha)
            cos_alpha = math.cos(alpha)
            cos2_alpha = cos_alpha * cos_alpha

            # Equation 18
            cos2_sigma_m = 0.0 if cos2_alpha == 0.0 else cos_sigma - 2.0 * sin_u1_sin_u2 / cos2_alpha
            u2 = cos2_alpha * a2b2b2

            cos2_sigma_m2 = cos2_sigma_m * cos2_sigma_m

            # Equation 3
            big_a = 1 + (u2 / 16384) * (4096 + u2 * (-768 + u2 * (320 - 175 * u2)))

            # Equation 4
            big_b = (u2 / 1024) * (256 + u2 * (-128 + u2 * (74 - 47 * u2)))

            # Equation 6
            delta_sigma = big_b * sin_sigma * (cos2_sigma_m + (big_b / 4.0) *
                                               (cos_sigma * (-1.0 + 2.0 * cos2_sigma_m2) - (big_b / 6.0) *
                                                cos2_sigma_m * (-3.0 + 4.0 * sin2_sigma) *
                                                (-3.0 + 4.0 * cos2_sigma_m2)))

            # Equation 10
            big_c = (f / 16.0) * cos2_alpha * (4.0 + f * (4.0 - 3.0 * cos2_alpha))

            # Equation 11
            lamda = omega + (1.0 - big_c) * f * sin_alpha *\
                            (sigma + big_c * sin_sigma *
                             (cos2_sigma_m + big_c * cos_sigma * (-1.0 + 2.0 * cos2_sigma_m2)))

            change = math.fabs((lamda - lamda0) / lamda)

            if i > 1 and change < self._tolerance:
                converged = True
                break

        # Equation 19
        s = b * big_a * (sigma - delta_sigma)
        alpha1 = Angle()
        alpha2 = Angle()

        if not converged:
            if phi1 > phi2:
                alpha1 = Angle.angle180
                alpha2 = Angle.zero
            elif phi1 < phi2:
                alpha1 = Angle.zero
                alpha2 = Angle.angle180
        else:
            # Equation 20
            radians = math.atan2(cos_u2 * math.sin(lamda), (cos_u1_sin_u2 - sin_u1_cos_u2 * math.cos(lamda)))
            if radians < 0.0:
                radians += self._two_pi
            alpha1.radians = radians

            # Equation 21
            radians = math.atan2(cos_u1 * math.sin(lamda), (-sin_u1_cos_u2 + cos_u1_sin_u2 * math.cos(lamda))) + math.pi
            if radians < 0.0:
                radians += self._two_pi
            alpha2.radians = radians

        if alpha1 >= 360.0:
            alpha1 -= 360.0
        if alpha2 >= 360.0:
            alpha2 -= 360.0

        geodetic_curve = GeodeticCurve(s, alpha1, alpha2)

        return geodetic_curve

    def calculate_geodetic_measurement(self, ellipsoid, start_position, end_position):
        start_coords = start_position.coordinates
        end_coords = end_position.coordinates

        # Calculate elevation differences
        elev1 = start_position.elevation
        elev2 = end_position.elevation
        elev12 = (elev1 + elev2) * 0.5

        # Calculate latitude differences
        phi1 = start_coords.latitude.radians
        phi2 = end_coords.latitude.radians
        phi12 = (phi1 + phi2) * 0.5

        # Calculate new ellipsoid to accommodate average elevation
        big_a = ellipsoid.semi_major_axis
        f = ellipsoid.flattening
        a = big_a + elev12 * (1.0 + f * math.sin(phi12))
        new_ellipsoid = Ellipsoid.from_semi_major_and_flattening(a, f)

        # Calculate the curve at the average elevation
        curve = self.calculate_geodetic_curve(new_ellipsoid, start_coords, end_coords)

        geodetic_measurement = GeodeticMeasurement(curve, elev2 - elev1)
        return geodetic_measurement
