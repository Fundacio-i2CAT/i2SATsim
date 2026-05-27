from datetime import datetime, timedelta
import ephem
from numpy import abs, arccos, arcsin, arctan, arctan2, cos, matrix, pi, sin, sqrt, tan, linspace
import numpy as np
from sgp4.earth_gravity import wgs72, wgs84
from satsim.node import Node
from satsim.channel_model import calculate_distance, free_space_path_loss, getElevation, calculate_sf_cl, getVisibility
from satsim.antenna_model import get_gain, get_transmit_power
from satsim.beams import getOneWebSemiMajorAxis, getOneWebSemiMinorAxis, getSuperellipseIndex, oneweb_beam_centers, oneweb_beams_polygons, starlink_beams
import math
from shapely.geometry import Point, Polygon as SPolygon
from shapely.affinity import scale
from cartopy.geodesic import Geodesic
from geopy.distance import geodesic
from shapely.geometry import Polygon
from math import degrees
from sgp4.ext import rv2coe
from sgp4.io import twoline2rv

XKMPER = 6378.135
F = 1 / 298.257223563  # Earth flattening WGS-84
c = 3 * 10**8  # Speed of light in m/s

class Sat(Node):
    __slots__ = ["cat", "incl", "RAAN0", "RAAN", "e", "w0", "w", "MA0", "MA",
                 "n", "epoch_year", "epoch_day", "theta", "GST0", "a", "mu",
                 "t0", "p",  "x", "y", "z", "Eq_r", "Po_r", "P_r", "Er_Pr2",
                 "J2", "P_w", "r", "h", "v", "tray_lat", "tray_lng", "B",
                 "bstar", "v_peri", "v_iner", "r_iner", "sat_model", "ndot",
                 "id_launch_data", "element_number", "satnumber", "tray_alt",
                 "line1", "line2", "tlepath", "revnum", "nlos_to_los_angle", 
                 "los_to_nlos_angle", "previous_sf", "previous_A_sp", "beam_centers", "beam_polygons", "covering_beams"]
    def __init__(self, name="", line1=None, line2=None, tlepath=None, cat="", date=None):
        """
        Parameters
        ----------
        name    : str, optional
                  The satellite's name (default is ).
        line1   : str, optional
                  First line of the two line elements.
        line2   : str, optional
                  Second line of the two line elements.
        tlepath : str, optional
                  Path to the satellite's TLE (default is None).
        cat     : str, optional
                  The satellite's category (default is ).
        """
        self.name = name
        self.cat = cat
        self.tray_lat = []
        self.tray_lng = []
        self.tray_alt = []
        if (tlepath is not None):
            self.readTLE(tlepath)
            self.tlepath = tlepath
        else:
            self.line1 = line1
            self.line2 = line2
        self.sat_model = twoline2rv(self.line1, self.line2, wgs84)#wgs72)
        self.incl = self.sat_model.inclo
        self.RAAN0 = self.sat_model.nodeo
        self.e = self.sat_model.ecco
        self.w0 = self.sat_model.argpo
        self.MA0 = self.sat_model.mo
        self.n = self.sat_model.no_kozai/60
        self.epoch_year = self.sat_model.epochyr
        self.epoch_day = self.sat_model.epochdays
        self.bstar = self.sat_model.bstar
        self.B = 2*self.bstar/(2.461*10**(-5)*6378.135)
        self.ndot = self.sat_model.ndot*1440**2/(2*pi)
        self.id_launch_data = self.sat_model.intldesg
        self.element_number = self.sat_model.elnum
        self.satnumber = self.sat_model.satnum
        self.revnum = self.sat_model.revnum
        #print(self.name + " TLE found!")
        self.RAAN = self.RAAN0
        self.w = self.w0
        self.MA = self.MA0
        self.theta = self.MA
        G = 6.67408*10**(-11)                          # Gravitational constant
        Mt = 5.9722*10**24                             # Earth mass
        self.t0 = self.epoch_day - int(self.epoch_day)
        self.t0 = self.t0*24*3600                      # Initial time in seconds
        self.updateGST0()                              # Get Greenwich sideral time
        self.a = (G*Mt/self.n**2)**(1/3)               # Semi-major axis
        self.v_peri = matrix([[0.0], [0.0]])           # Velocity in the perifocal frame
        self.v_iner = matrix([[0.0], [0.0], [0.0]])    # Velocity in the inertial frame
        self.r_iner = matrix([[0.0], [0.0], [0.0]])    # Position in the inertial frame
        self.changePlanet()                            # Sets all the planet's parameters
        self.updateOrbitalParameters()
        self.getLat()
        if date is not None:
            self.getLng(date=date)
            self.getLat(date=date)
        
        else: self.getLng()
        self.getAlt()
        self.nlos_to_los_angle = 0
        self.los_to_nlos_angle = 0
        self.previous_sf = None
        self.previous_A_sp = None       
        self.beam_centers = []
        self.beam_polygons = []
        self.covering_beams = []

    def __call__(self):
        return self

    def readTLE(self, file_name):
        """
        Search for the satellite name inside the file_name. If found,
        sets the next two lines as the two line elements of this
        satellite.
        """
        with open(file_name, 'r') as f:
            for count, line in enumerate(f):
                if (line.strip() == self.name):
                    self.line1 = next(f).strip()
                    self.line2 = next(f).strip()
                    break

    def setMu(self, mu):
        """
        Sets the standard gravitational parameter.

        Parameters
        ----------
        mu : float
             The new standard gravitational parameter.
        """
        self.mu = mu

    def setInclination(self, incl):
        """
        Sets the inclination of the orbit.

        Parameters
        ----------
        incl : float
               The new satellite's orbit inclination in radians.
        """
        self.incl = incl

    def setRAAN(self, RAAN):
        """
        Sets the orbit's right ascension of the ascending node.

        Parameters
        ----------
        RAAN : float
               The new Right Ascension of the Ascending Node in radians.
        """
        self.RAAN = RAAN

    def setRAAN0(self, RAAN0):
        """
        Sets the orbit's initial RAAN.

        Parameters:
        -----------
        RAAN0 : float
                The initial Right Ascension of the Ascending Node in
                radians.
        """
        self.RAAN0 = RAAN0

    def setArgPerigee(self, w):
        """
        Sets the orbit's argument of the perigee.

        Parameters
        ----------
        w : float
            The new orbit's argument of the perigee in radians.
        """
        self.w = w

    def setArgPerigee0(self, w0):
        """
        Sets the orbit's initial argument of the perigee.

        Parameters
        ----------
        w0 : float
             The initial orbit's argument of the perigee in radians.
        """
        self.w0 = w0

    def setEccentricity(self, e):
        """
        Sets the orbit's eccentricity.

        Parameters
        ----------
        e : float
            The new orbit's eccentricity.
        """
        self.e = e

    def setMeanAnomaly(self, MA):
        """
        Sets the orbit's mean anomaly.

        Parameters
        ----------
        MA : float
             The new orbit's mean anomaly in radians.
        """
        self.MA = MA

    def setMeanAnomaly0(self, MA0):
        """
        Sets the orbit's initial mean anomaly.

        Parameters
        ----------
        MA0 : float
              The initial orbit's mean aomaly in radians.
        """
        self.MA0 = MA0

    def setMeanMotion(self, n):
        """
        Sets the satellite's mean motion.

        Parameters
        ----------
        n = float
            The satellite's mean motion.
        """
        self.n = n

    def setTrueAnomaly(self, theta):
        """
        Sets the orbit's true anomaly.

        Parameters
        ----------
        theta : float
                The new orbit's true anomaly in radians.
        """
        self.theta = theta

    def setSemiMajorAxis(self, a):
        """
        Sets the orbit's semi-major axis.

        Parameters
        ----------
        a = float
            The new orbit's semi-major axis in meters.
        """
        self.a = a

    def setSemilatusRectum(self, p):
        """
        Sets the orbit's semilatus rectum.

        Parameters
        ----------
        p : float
            The new orbit's semilatus rectum in meters.
        """
        self.p = p

    def setBallisticCoeff(self, B):
        """
        Sets the ballistic coefficient of the satellite.

        Parameters
        ----------
        B : float
            Ballistic coefficient.
        """
        self.B = B

    def setSpecAngMomentum(self, h):
        """
        Sets the satellite's specific relative angular momentum.

        Parameters
        ----------
        h : float
            The new specific relative angular momentum in squared meters
            per second.
        """
        self.h = h

    def setName(self, name):
        """
        Sets the satellite's name.

        Parameters
        ----------
        name : str
               The satellite's new name.
        """
        self.name = name
    
    def setCategory(self, cat):
        """
        Sets the satellite's category.

        Parameters
        ----------
        cat : str
              The satellite's new category.
        """
        self.cat = cat

    def setNLOStoLOS(self, nlos_to_los_angle):
        self.nlos_to_los_angle = nlos_to_los_angle
    
    def setLOStoNLOS(self, los_to_nlos_angle):
        self.los_to_nlos_angle = los_to_nlos_angle

    def setPreviousSF(self, previous_sf):
        self.previous_sf = previous_sf

    def setPreviousAsp(self, previous_A_sp):
        self.previous_A_sp = previous_A_sp

    def setBeamCenters(self, beam_centers):
        self.beam_centers = beam_centers
    
    def setBeamPolygons(self, beam_polygons):
        self.beam_polygons = beam_polygons

    def setCoveringBeams(self, covering_beams):
        self.covering_beams = covering_beams

    def getCategory(self):
        """Returns the satellite's category."""
        return self.cat

    def getInclination(self):
        """Returns the inclination of the orbit."""
        return self.incl

    def getRAAN(self):
        """Returns the orbit's RAAN."""
        return self.RAAN

    def getArgPerigee(self):
        """Returns the argument of the perigee of the orbit."""
        return self.w
    
    def getEccentricity(self):
        """Returns the eccentricity of the orbit."""
        return self.e
    
    def getSemiMajorAxis(self):
        """Returns the orbit's semi-major axis."""
        return self.a

    def getAnomaly(self):
        """Returns the orbit's true anomaly."""
        return self.theta

    def getMeanAnomaly(self):
        """Returns the orbit's mean anomaly."""
        return self.MA

    def getSpecAngMomentum(self):
        """
        Returns the satellite's specific relative angular momentum.
        """
        return self.h

    def getSpeed(self):
        """Returns the speed of the satellite."""
        mu_div_h = self.mu/self.h
        v_peri_p = -sin(self.theta)*mu_div_h
        v_peri_q = (self.e + cos(self.theta))*mu_div_h
        self.v = sqrt(v_peri_p**2 + v_peri_q**2)
        return self.v

    def getPerifocalVel(self):
        """
        Returns the satellite's velocity with respect to the perifocal
        frame.
        """
        mu_div_h = self.mu/self.h
        self.v_peri[0,0] = -sin(self.theta)*mu_div_h
        self.v_peri[1,0] = (self.e + cos(self.theta))*mu_div_h
        return self.v_peri[0,0], self.v_peri[1,0]

    def getInertialVel(self):
        """
        Calculates the velocity relative to its Perifocal Frame and
        transforms it to the inertial frame (Geocentric Equatorial
        Frame).

        Returns
        -------
        self.v_iner
            The velocity with respect to the inertial frame.
        """
        return self.v_iner

    def getXYZ(self):
        """
        Returns the position of the satellite with respect to the
        inertial frame.
        """
        return self.r_iner

    def getLat(self, rad2deg=180/pi, date=None):
        """
        Returns the last latitude of the satellite.
        """
        if date is None:
            r_ecef = sqrt(self.x**2 + self.y**2 + self.z**2)
            lat = arcsin(self.z / r_ecef)
            self.lat = lat * rad2deg
        
        else:
            tle_rec = ephem.readtle(self.name, self.line1, self.line2)
            tle_rec.compute(date)
            self.lat = degrees(tle_rec.sublat)

        return self.lat

    def getLng(self, rad2deg=180/pi, twopi=2*pi, date=None):
        """
        Returns the last longitude of the satellite.
        """

        if date is None:
            tnow = self.getTnow(date)
            GST = self.GST0 + self.P_w * tnow
            lng_rad = arctan2(self.y, self.x) - GST
            lng_rad = (lng_rad + pi) % twopi - pi
            self.lng = lng_rad * rad2deg
        
        else:
            tle_rec = ephem.readtle(self.name, self.line1, self.line2)
            tle_rec.compute(date)
            self.lng = degrees(tle_rec.sublong)

        return self.lng

    def getAlt(self):
        """
        Calculates the satellite's altitude considering the planet's
        radius at the current coordinates.

        Returns
        -------
        float
            Altitude in meters
        """
        self.alt = self.r - self.getPlanetRadius()
        return self.alt
    
    def getPosition(self, date):
        """
        Calculate the geographical position (longitude and latitude) 
        of the object at a given date based on its TLE (Two-Line Element) data.

        Parameters:
        - date: A date object representing the time at which the position 
                is to be calculated.
        """
        tle_rec = ephem.readtle(self.name, self.line1, self.line2)
        tle_rec.compute(date)

        lat = degrees(tle_rec.sublat)
        long = degrees(tle_rec.sublong)
        h = self.getAlt()

        return long, lat, h
    
    def getNLOStoLOS(self):
        return self.nlos_to_los_angle
    
    def getLOStoNLOS(self):
        return self.los_to_nlos_angle
    
    def getPreviousSF(self):
        return self.previous_sf
    
    def getPreviousAsp(self):
        return self.previous_A_sp
    
    def getBeamCenters(self):
        return self.beam_centers
    
    def getBeamPolygons(self):
        return self.beam_polygons
    
    def getCoveringBeams(self):
        return self.covering_beams

    def computeRSRP(self, date, terminal_coords, scenario, time_step_seconds, frequency, max_gain_sat, antenna_aperture, gain_term, eirp_density, bandwidth, scs, guardband):
        """
        Compute the received signal power (RSRP) at the terminal from the satellite.
        
        Parameters:
        terminal_coords (tuple): (x, y, z) coordinates of the terminal.
        
        Returns:
        float: Received signal power (RSRP) in dB.
        """

        sat_coords = self.getPosition(date)
        distance = calculate_distance(sat_coords, terminal_coords)

        next_date = date + timedelta(seconds=time_step_seconds)
        next_sat_coords = self.getPosition(next_date)
        next_distance = calculate_distance(next_sat_coords, terminal_coords)

        elevation = getElevation(sat_coords, terminal_coords)

        if next_distance > distance:
            elevation = 180 - elevation

        if self.nlos_to_los_angle == 0 and self.los_to_nlos_angle == 0:
            trials = 100
            nlos_to_los_angle, los_to_nlos_angle = getVisibility(scenario, trials)
            self.nlos_to_los_angle = nlos_to_los_angle
            if los_to_nlos_angle != 0:
                self.los_to_nlos_angle = los_to_nlos_angle
            else:
                self.los_to_nlos_angle = 170

        sf, cl, A_sp = calculate_sf_cl(elevation, distance, frequency, scenario, self.nlos_to_los_angle, self.los_to_nlos_angle, terminal_coords[2], self.previous_sf, self.previous_A_sp)
        self.previous_sf = sf
        self.previous_A_sp = A_sp
 
        fspl = free_space_path_loss(distance, frequency)

        pl = fspl + sf + cl

        rsrp_beams = []
        for b in self.covering_beams:

            beam_center = self.beam_centers[b]

            a_ellipse = getOneWebSemiMajorAxis(beam_center[1])
            b_ellipse = getOneWebSemiMinorAxis()

            distance_x = np.abs(terminal_coords[0] - beam_center[0])
            distance_y = np.abs(terminal_coords[1] - beam_center[1])

            n = getSuperellipseIndex()
            d_superellipse = ((np.abs(distance_x / a_ellipse) ** n) + (np.abs(distance_y / b_ellipse) ** n)) ** (1/n)


            theta_max = np.rad2deg(np.atan(b_ellipse*111000/sat_coords[2]))
            theta = theta_max * d_superellipse

            gain_sat = get_gain(frequency, theta, max_gain_sat, antenna_aperture) # Gain of satellite antenna in dB

            transmit_power = get_transmit_power(eirp_density, bandwidth, max_gain_sat)  # Transmit power in dBm

            rb_bw = scs*12
            n_rb = (bandwidth - 2*guardband)/rb_bw
            n_re = 12*n_rb

            std_dev_error = 1.72
            measurement_error = np.random.normal(0, std_dev_error)

            rsrp = transmit_power + gain_sat + gain_term - pl - 10*math.log10(n_re) + measurement_error

            rsrp_beams.append(rsrp)

        return rsrp_beams

    def getDMY(self):
        """
        Calculates the day, month and year from the TLE's epoch data and
        returns the result.
        """
        year = self.epoch_year
        ep_day = self.epoch_day
        if (year % 4 == 0):
            if (year % 100 == 0):
                if (year % 400 == 0):
                    day, month = self.leapYearDM(ep_day)
                else:
                    day, month = self.notLeapYearDM(ep_day)
            else:
                day, month = self.leapYearDM(ep_day)
        else:
            day, month = self.notLeapYearDM(ep_day)
        return day, month, year

    def leapYearDM(self, ep_day):
        """
        Calculates the days and months from an epoch day parameter,
        assuming it is a leap year.

        Parameters
        ----------
        ep_day : float
                 Satellite's epoch day from TLE.
        """
        month = 1 + (ep_day > 31) + (ep_day > 60) + (ep_day > 91)
        month += (ep_day > 121) + (ep_day > 152) + (ep_day > 182)
        month += (ep_day > 213) + (ep_day > 244) + (ep_day > 274)
        month += (ep_day > 305) + (ep_day > 335)

        day = ep_day - (ep_day > 31)*31 - (ep_day > 59)*29
        day = day - (ep_day > 90)*31 - (ep_day > 120)*30
        day = day - (ep_day > 151)*31 - (ep_day > 181)*30
        day = day - (ep_day > 212)*31 - (ep_day > 243)*31
        day = day - (ep_day > 273)*30 - (ep_day > 304)*31
        day = day - (ep_day > 334)*30
        return day, month

    def notLeapYearDM(self, ep_day):
        """
        Calculates the days and months from an epoch day parameter,
        assuming it is not a leap year.

        Parameters
        ----------
        ep_day : float
                 Satellite's epoch day from TLE.
        """
        month = 1 + (ep_day > 31) + (ep_day > 59) + (ep_day > 90)
        month += (ep_day > 120) + (ep_day > 151) + (ep_day > 181)
        month += (ep_day > 212) + (ep_day > 243) + (ep_day > 273)
        month += (ep_day > 304) + (ep_day > 334)

        day = ep_day - (ep_day > 31)*31 - (ep_day > 59)*28
        day = day- (ep_day > 90)*31 - (ep_day > 120)*30
        day = day - (ep_day > 151)*31 - (ep_day > 181)*30
        day = day - (ep_day > 212)*31 - (ep_day > 243)*31
        day = day - (ep_day > 273)*30 - (ep_day > 304)*31
        day = day - (ep_day > 334)*30
        return day, month

    def getGST(self, D, M, Y):
        """
        Returns the Greenwich Sidereal Time.

        Parameters
        ----------
        D : int
            Days since the beginning of the month.
        M : int
            Current month.
        Y : int
            Current year.

        Returns
        -------
        GST0
            The Greenwich Sidereal Time at t0.
        """
        JD = 367*Y - int(7/4*(Y + int((M + 9)/12))) + int(275*M/9) + D + 1721013.5 # Julian day
        T0 = (JD - 2451545)/36525
        GST0 = (100.4606184 + 36000.77004*T0 + 0.000387933*T0**2 - 2.583*10**(-8)*T0**3)*pi/180
        return GST0

    def getCurrentTimeInSeconds(self, date):
        """
        Returns the number of seconds since the beginning of the day.

        Parameters
        ----------
        date : datetime.utcnow()
               Current date.
        """
        seconds = date.hour*3600 + date.minute*60 + date.second
        return seconds + date.microsecond*0.000001

    def month2days(self, date):
        """
        Returns the number of months in days.

        Parameters
        ----------
        date : datetime
               Date from the datetime library.
        """
        days = date.day
        month = date.month
        year = date.year
        days = days + (month > 1)*31 + (month > 2)*28 + (month > 3)*31 + (month > 4)*30
        days = days + (month > 5)*31 + (month > 6)*30 + (month > 7)*31 + (month > 8)*31
        days = days + (month > 9)*30 + (month > 10)*31 + (month > 11)*30
        if (year % 4 != 0):
            days = days
        elif (year % 100 != 0):
            days = days + (month > 2)
        elif (year % 400 != 0):
            days = days
        else:
            days = days + (month > 2)
        return days

    def M2E(self, M):
        """
        Transforms the mean anomaly to the eccentric anomaly.

        Parameters
        ----------
        M : float
            The mean anomaly in radians.
        """
        E0 = 0
        E = M
        while (abs(E - E0) > 0.0000001):
            E0 = E
            E = M + self.e*sin(E0)
        return E

    def E2theta(self, E):
        """
        Transforms the eccentric anomaly to the true anomaly.

        Parameters
        ----------
        E : float
            The eccentric anomaly in radians.
        """
        theta = 2*arctan(sqrt((1 + self.e)/(1 - self.e))*tan(E/2))
        return theta

    def getPlanetRadius(self):
        """
        Returns the planet radius considering the polar radius and the
        equatorial radius, given the current latitude and longitude.
        """
        cos_lat = cos(math.radians(self.lat))
        sin_lat = sin(math.radians(self.lat))
        aux = ((self.Eq_r**2*cos_lat)**2 + (self.Po_r**2*sin_lat)**2)
        radius = sqrt(aux/((self.Eq_r*cos_lat)**2 + (self.Po_r*sin_lat)**2))
        return radius

    def getPeriod(self):
        """
        Returns
        -------
        float
            The period of the satellite's orbit.
        """
        return 2*pi/self.n

    def getCoverage(self):
        """
        Returns
        -------
        float
            The angle of the coverage considering the planet's radius at
            the satellite's coordinates and the satellite's altitude. To
            return the angle in degrees, it is multiplied by
            57.29577951308232 which is 180/pi.
            self.r is the Planet's radius plus the satellite's altitude.
        """
        radius = self.getPlanetRadius()
        ang = arccos(radius/self.r)*57.29577951308232
        return ang

    def getComCoverage(self, E_r, c=299792458):
        """
        Calculates the angle of the coverage of the comunication link
        considering the transmit and receive power, the antenna gain,
        the sensibility and the losses.

        Parameters
        ----------
        E_r : float
              The Earth radius in meters.
        c   : int
              The speed of light in meters per second.

        Returns
        -------
        float
            The angle of the coverage in degrees.
        """
        # FSPL = Ptx + Ant_tx + Ant_rx - Cable + Sens - margin
        #SUCHAI
        alt = self.getAlt()
        FSPL = 30 + 2 + 18.9 - 1 + 117 - 12
        d = 10**(FSPL/20)*(c/self.freq)/(4*pi)
        ang = arccos((E_r**2 + (E_r + alt)**2 - d**2)/(2*E_r*(E_r + alt)))*180/pi
        return ang

    def getTnow(self, date=None):
        """
        Calculates the number of seconds since the last TLE data.

        Parameters
        ----------
        date : datetime
               The date used to calculate the number of seconds.

        Returns
        -------
        float
            Number of seconds.
        """
        if (date is None):
            date = datetime.utcnow()              # Use current time in UTC.
        dayinsec = 86400                          # Day in seconds
        tnow = self.getCurrentTimeInSeconds(date) # Current time in seconds
        days = self.month2days(date)              # Days in current time
        daysdiff = days - int(self.epoch_day)     # Difference between TLE date and current date
        return tnow + daysdiff*dayinsec           # Time in seconds from TLE to present time

    def getTrayectory(self, T, dt, date=None):
        """
        Calculates the future trayectory of the satellite using the
        SGP4 library.

        Parameters
        ----------
        T    : int
               Period of time in seconds to be calculated.
        dt   : int
               Time step in seconds.
        date : datetime
               Date from which the trayectory is calculated.
        """
        self.tray_lat[:] = []
        self.tray_lng[:] = []
        self.tray_alt[:] = []
        if (date is None):
            date = datetime.utcnow()
        current_date = date
        for t in range(0, T, dt):
            self.updateOrbitalParameters(date)
            date += timedelta(seconds=dt)
            self.tray_lat.append(self.getLat(date=date))
            self.tray_lng.append(self.getLng(date=date))
            self.getAlt()
            self.tray_alt.append(self.alt)
        self.updateOrbitalParameters(current_date)
        return self.tray_lat, self.tray_lng

    def changePlanet(self, M=5.9722*10**24, P_r=6371000, Eq_r=6378000, Po_r=6356000, J2=0.00108263, P_w=7.29211505*10**(-5)):
        """
        Changes the planet's parameters.

        Parameters
        ----------
        M    : float
               The planet's mass in kilograms.
        P_r  : int
               The planet's mean radius in meters.
        Eq_r : int
               The planet's equatorial radius in meters.
        Po_r : int
               The planet's polar radius in meters.
        J2   : float
               The planet's second degree harmonic model.
        P_w  : float
               The planet's angular velocity in radians per second.
        """
        G = 6.67408*10**(-11)              # Gravitational constant
        self.mu = G*M                      # Gravitational parameter
        self.n = sqrt(self.mu/(self.a**3)) # Mean motion
        self.p = self.a*(1 - self.e**2)
        self.h = sqrt(self.p*self.mu)
        self.P_r = P_r
        self.Eq_r = Eq_r
        self.Po_r = Po_r
        self.Er_Pr2 = (Eq_r/Po_r)**2
        self.J2 = J2
        self.P_w = P_w

    def updateOrbitalParameters(self, date=None):
        """
        Updates all the orbital parameters with the SGP4 propagation
        model.

        date : datetime
               The orbital elements are updated for this date.
        """
        if (date is None):
            date = datetime.utcnow()
        year = date.year
        month = date.month
        day = date.day
        hour = date.hour
        minute = date.minute
        second = date.second + date.microsecond*0.000001
        pos, vel = self.sat_model.propagate(year, month, day,
                                            hour, minute, second)
        p, a, e, i, raan, w, theta, m, argl, tlon, lonp = rv2coe(pos, vel,
                                                           self.mu*0.000000001)
        self.x = pos[0]*1000
        self.y = pos[1]*1000
        self.z = pos[2]*1000
        self.r = sqrt(self.x**2 + self.y**2 + self.z**2)
        self.r_iner[0,0] = self.x
        self.r_iner[1,0] = self.y
        self.r_iner[2,0] = self.z
        self.v_iner[0,0] = vel[0]*1000
        self.v_iner[1,0] = vel[1]*1000
        self.v_iner[2,0] = vel[2]*1000
        self.p = p*1000
        self.a = a*1000
        self.e = e
        self.incl = i
        self.RAAN = raan
        self.w = w
        self.theta = theta
        self.MA = m
        if self.a > 0:
           self.n = sqrt(self.mu/(self.a**3))
        self.h = sqrt(abs(self.a*self.mu*(1 - self.e**2)))

    def updateGST0(self):
        """
        Returns the Greenwich Sidereal Time of the TLE data.
        """
        D, Month, Y = self.getDMY()
        self.GST0 = self.getGST(int(D), Month, Y)

    def updateEpoch(self, date=None):
        """
        Updates the epoch day, year and the GST0 to now or to the
        received date.

        Parameters
        ----------
        date : datetime
               Date received to update the epoch.
        """
        if (date is None):
            date = datetime.utcnow()
        tnow = self.getCurrentTimeInSeconds(date)
        days = self.month2days(date)
        self.epoch_day = days + tnow/86400
        self.t0 = tnow
        self.epoch_year = date.year
        self.updateGST0()

    def setTLE(self, line1, line2):
        checksum1 = self.checksum(line1)
        checksum2 = self.checksum(line2)
        self.line1 = "{}{}".format(line1[0:-1], checksum1)
        self.line2 = "{}{}".format(line2[0:-1], checksum2)
        self.sat_model = twoline2rv(self.line1, self.line2, wgs84)#wgs72)

    def getTLE(self):
        """
        Returns
        -------
        str
            The satellite's TLE as a string.
        """
        return "{}\n{}\n{}".format(self.name, self.line1, self.line2)

    def checksum(self, line):
        """
        Calculates the line's checksum.

        line : str
               String used to calculate its checksum.
        """
        check = 0
        for char in line[:-1]:
            if char.isdigit():
                check += int(char)
            if char == "-":
                check += 1
        return check % 10

    def createTLE(self, date):
        """
        Creates a new TLE for the received date.

        Parameters
        ----------
        date : datetime
               Date used to obtain the epoch of the TLE.
        """
        rad2deg = 180/pi
        old_epoch_day = self.epoch_day
        self.updateEpoch(date)
        new_epoch_day = self.epoch_day
        revolutions = (new_epoch_day - old_epoch_day)*86400/self.getPeriod()
        revolutions = int(revolutions) + int(self.revnum)
        aux="{: .9f}".format(self.ndot)
        ndot = "{}{}".format(aux[0],aux[2:-1])
        aux="{: e}".format(self.B*2.461*10**(-5)*6378.135*0.5)
        BSTAR = "{}{}-{}".format(aux[0:2], aux[3:7], int(aux[-1])-1)
        aux="{:.8f}".format(self.e)
        e = "{}".format(aux[2:-1])
        tle_num = "{:4d}".format(self.element_number)
        epoch_year = self.epoch_year - 2000
        line1 = "1 {:05}U {:9}{}{:012.8f} {}  00000-0 {} 0 {}7".format(self.satnumber,
                                                           self.id_launch_data,
                                                           epoch_year,
                                                           self.epoch_day,
                                                           ndot,
                                                           BSTAR,
                                                           tle_num)
        line2 = "2 {:05} {:8.4f} {:8.4f} {} {:8.4f} {:8.4f} {:11.8f}{:5d}7".format(self.satnumber,
                                                             self.incl*rad2deg,
                                                             self.RAAN*rad2deg,
                                                             e,
                                                             self.w*rad2deg,
                                                             self.MA*rad2deg,
                                                             86400/self.getPeriod(),
                                                             revolutions)
        self.setTLE(line1, line2)
        return self.name, self.line1, self.line2

    def terminal_area(self,terminal_coords, radius):
        """
        Create a regular hexagon for the terminal area based on center coordinates and size.
        :param terminal_coords: Tuple of (x, y) for the terminal coordinates
        :param size: Size (distance from center to vertex) of the hexagon
        :return: Hexagon (Polygon) representing the terminal area
        """

        x, y, z = terminal_coords
        # Hexagon vertices: Regular hexagon has 6 vertices equally spaced around the center
        angles = [i * 60 for i in range(6)]  # 6 vertices, each 60 degrees apart
        vertices = [(x + radius * math.cos(math.radians(angle)), y + radius * math.sin(math.radians(angle))) for angle in angles]

        return SPolygon(vertices)

    def getCoverageCircle(self, sat_coords, coverage_radius):
        """
        Create a circular coverage area based on satellite coordinates and radius.
      
        :param satellite_coords: Tuple of (x, y) for the satellite coordinates
        :param radius: Radius of the satellite beamCoverageC (m)
        :return: Circle (Polygon) representing the coverage area
        """
        
        lon, lat = sat_coords[0], sat_coords[1]
        circle_points = Geodesic().circle(lon=lon, lat=lat, radius=coverage_radius, n_samples=200)
        coverage_polygon = Polygon(circle_points)

        return coverage_polygon


    def is_under_coverage(self, date, terminal_coords):
        """
        Determine if a specified point (terminal) is within the coverage area 
        defined by a circle at a given date and coverage radius.

        Parameters:
        - date: The date for which coverage is being checked.
        - terminal_coords: A tuple or list containing the (longitude, latitude) 
                        coordinates of the terminal point to check.
        """
        sat_coords = self.getPosition(date)
        elevation = getElevation(sat_coords, terminal_coords)

        self.covering_beams = []
        coverage = False
        
        if elevation < 0 or elevation > 90:
            return coverage
        else:
            i = 0
            for b in self.beam_polygons:
                beam_polygon = Polygon(b)
                point = Point(terminal_coords[0], terminal_coords[1])
                if beam_polygon.contains(point):
                    self.covering_beams.append(i)
                    coverage = True
                i += 1
            return coverage
            
    def is_visible(self, terminal_coords, date, constellation):
        """
        Check if the satellite beam (circle) and terminal area (hexagon) intersect.

        :param satellite_coords: (x, y) tuple for the satellite location
        :param terminal_coords: (x, y) tuple for the terminal location
        :return: True if the two areas intersect, False otherwise
        """
        sat_coords = self.getPosition(date)

        elevation = getElevation(sat_coords, terminal_coords)
        if elevation > 10:
            if constellation == "OneWeb":
                self.beam_centers = oneweb_beam_centers(sat_coords)
                self.beam_polygons = oneweb_beams_polygons(self.beam_centers)

            elif constellation == "Starlink":
                coverage_radius = 500000
                isd = coverage_radius*2*0.9
                if self.beam_centers == []:
                    self.beam_centers.append(sat_coords[0:2])
                    beam_center = self.beam_centers[0]
                else:
                    beam_center = self.beam_centers[0]
                    distance = geodesic(sat_coords[0:2], beam_center).m
                    if distance > isd/2:
                        time = np.round((isd/2)/self.getSpeed())
                        date2 = date + timedelta(seconds=time)
                        beam_center = self.getPosition(date2)[0:2]
                self.beam_centers, self.beam_polygons = starlink_beams(beam_center, coverage_radius)
            
            return True
        else:
            self.beam_centers = []
            self.beam_polygons = []
            return False
        