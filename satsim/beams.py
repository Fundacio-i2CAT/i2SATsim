import numpy as np
from cartopy.geodesic import Geodesic
from numpy import asarray

def getOneWebSemiMajorAxis(lat):
    a_ellipse = 1900/(2*111*np.cos(np.deg2rad(lat)))
    return a_ellipse

def getOneWebSemiMinorAxis():
    b_ellipse = 1300/(2*16*111)
    return b_ellipse

def getSuperellipseIndex():
    n = 3.5
    return n

def oneweb_beam_centers(sat_coords):
    long = sat_coords[0]
    lat = sat_coords[1]
    b_ellipse = getOneWebSemiMinorAxis()
    isd = 0.78*b_ellipse*2
    beams = 16
    beam_centers = []
    for b in range(0, beams):
        space = isd/2 + (b-8)*isd
        beam_center = (long, lat + space)
        beam_centers.append(beam_center)
    return beam_centers

def oneweb_beams_polygons(beam_centers):
    n_samples = 100
    b_ellipse = getOneWebSemiMinorAxis()
    beam_polygons = []
    for beam_center in beam_centers:
        a_ellipse = getOneWebSemiMajorAxis(beam_center[1])
        xy_data = []
        for k in range(n_samples):
            theta = 2 * np.pi * k / n_samples
            #x = beam_center[0] + a_ellipse * np.cos(theta)
            #y = beam_center[1] + b_ellipse * np.sin(theta)
            n = getSuperellipseIndex()
            x = beam_center[0] + a_ellipse * np.sign(np.cos(theta)) * np.abs(np.cos(theta))**(2/n)
            y = beam_center[1] + b_ellipse * np.sign(np.sin(theta)) * np.abs(np.sin(theta))**(2/n)
            xy_data.append([x, y])
        beam_polygons.append(xy_data)
    return beam_polygons

def starlink_beams(beam_center, coverage_radius):
    n_samples = 100
    beam_centers = []
    beam_polygons = []
    beam_centers.append(beam_center)
    xy_data = asarray(Geodesic().circle(lon=beam_center[0], lat=beam_center[1], radius=coverage_radius, n_samples=n_samples, endpoint=True))
    beam_polygons.append(xy_data)
    return beam_centers, beam_polygons
