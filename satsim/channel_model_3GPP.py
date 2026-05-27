import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import math

los_probabilities = {
    "Dense Urban": {
        10: 0.282, 20: 0.331, 30: 0.398, 40: 0.468, 50: 0.537,
        60: 0.612, 70: 0.738, 80: 0.82, 90: 0.981
    },
    "Urban": {
        10: 0.246, 20: 0.386, 30: 0.493, 40: 0.613, 50: 0.726,
        60: 0.805, 70: 0.919, 80: 0.968, 90: 0.992
    },
    "Suburban": {
        10: 0.782, 20: 0.869, 30: 0.919, 40: 0.929, 50: 0.935,
        60: 0.94, 70: 0.949, 80: 0.952, 90: 0.998
    }
}

dense_urban_sf_cl = {
    10: {
        'S-band': {'LOS': {'std_dev': 3.5}, 'NLOS': {'std_dev': 15.5, 'CL': 34.3}},
        'Ka-band': {'LOS': {'std_dev': 2.9}, 'NLOS': {'std_dev': 17.1, 'CL': 44.3}}
    },
    20: {
        'S-band': {'LOS': {'std_dev': 3.4}, 'NLOS': {'std_dev': 13.9, 'CL': 30.9}},
        'Ka-band': {'LOS': {'std_dev': 2.4}, 'NLOS': {'std_dev': 17.1, 'CL': 39.9}}
    },
    30: {
        'S-band': {'LOS': {'std_dev': 2.9}, 'NLOS': {'std_dev': 12.4, 'CL': 29.0}},
        'Ka-band': {'LOS': {'std_dev': 2.7}, 'NLOS': {'std_dev': 15.6, 'CL': 37.5}}
    },
    40: {
        'S-band': {'LOS': {'std_dev': 3.0}, 'NLOS': {'std_dev': 11.7, 'CL': 27.7}},
        'Ka-band': {'LOS': {'std_dev': 2.4}, 'NLOS': {'std_dev': 14.6, 'CL': 35.8}}
    },
    50: {
        'S-band': {'LOS': {'std_dev': 3.1}, 'NLOS': {'std_dev': 10.6, 'CL': 26.8}},
        'Ka-band': {'LOS': {'std_dev': 2.4}, 'NLOS': {'std_dev': 14.2, 'CL': 34.6}}
    },
    60: {
        'S-band': {'LOS': {'std_dev': 2.7}, 'NLOS': {'std_dev': 10.5, 'CL': 26.2}},
        'Ka-band': {'LOS': {'std_dev': 2.7}, 'NLOS': {'std_dev': 12.6, 'CL': 33.8}}
    },
    70: {
        'S-band': {'LOS': {'std_dev': 2.5}, 'NLOS': {'std_dev': 10.1, 'CL': 25.8}},
        'Ka-band': {'LOS': {'std_dev': 2.6}, 'NLOS': {'std_dev': 12.1, 'CL': 33.3}}
    },
    80: {
        'S-band': {'LOS': {'std_dev': 2.3}, 'NLOS': {'std_dev': 9.2, 'CL': 25.5}},
        'Ka-band': {'LOS': {'std_dev': 2.8}, 'NLOS': {'std_dev': 12.3, 'CL': 33.0}}
    },
    90: {
        'S-band': {'LOS': {'std_dev': 1.2}, 'NLOS': {'std_dev': 9.2, 'CL': 25.5}},
        'Ka-band': {'LOS': {'std_dev': 0.6}, 'NLOS': {'std_dev': 12.3, 'CL': 32.9}}
    }
}

urban_sf_cl = {
    10: {
        'S-band': {'LOS': {'std_dev': 4.0}, 'NLOS': {'std_dev': 6.0, 'CL': 34.3}},
        'Ka-band': {'LOS': {'std_dev': 4.0}, 'NLOS': {'std_dev': 6.0, 'CL': 44.3}}
    },
    20: {
        'S-band': {'LOS': {'std_dev': 4.0}, 'NLOS': {'std_dev': 6.0, 'CL': 30.9}},
        'Ka-band': {'LOS': {'std_dev': 4.0}, 'NLOS': {'std_dev': 6.0, 'CL': 39.9}}
    },
    30: {
        'S-band': {'LOS': {'std_dev': 4.0}, 'NLOS': {'std_dev': 6.0, 'CL': 29.0}},
        'Ka-band': {'LOS': {'std_dev': 4.0}, 'NLOS': {'std_dev': 6.0, 'CL': 37.5}}
    },
    40: {
        'S-band': {'LOS': {'std_dev': 4.0}, 'NLOS': {'std_dev': 6.0, 'CL': 27.7}},
        'Ka-band': {'LOS': {'std_dev': 4.0}, 'NLOS': {'std_dev': 6.0, 'CL': 35.8}}
    },
    50: {
        'S-band': {'LOS': {'std_dev': 4.0}, 'NLOS': {'std_dev': 6.0, 'CL': 26.8}},
        'Ka-band': {'LOS': {'std_dev': 4.0}, 'NLOS': {'std_dev': 6.0, 'CL': 34.6}}
    },
    60: {
        'S-band': {'LOS': {'std_dev': 4.0}, 'NLOS': {'std_dev': 6.0, 'CL': 26.2}},
        'Ka-band': {'LOS': {'std_dev': 4.0}, 'NLOS': {'std_dev': 6.0, 'CL': 33.8}}
    },
    70: {
        'S-band': {'LOS': {'std_dev': 4.0}, 'NLOS': {'std_dev': 6.0, 'CL': 25.8}},
        'Ka-band': {'LOS': {'std_dev': 4.0}, 'NLOS': {'std_dev': 6.0, 'CL': 33.3}}
    },
    80: {
        'S-band': {'LOS': {'std_dev': 4.0}, 'NLOS': {'std_dev': 6.0, 'CL': 25.5}},
        'Ka-band': {'LOS': {'std_dev': 4.0}, 'NLOS': {'std_dev': 6.0, 'CL': 33.0}}
    },
    90: {
        'S-band': {'LOS': {'std_dev': 4.0}, 'NLOS': {'std_dev': 6.0, 'CL': 25.5}},
        'Ka-band': {'LOS': {'std_dev': 4.0}, 'NLOS': {'std_dev': 6.0, 'CL': 32.9}}
    }
}

suburban_sf_cl = {
    10: {
        'S-band': {'LOS': {'std_dev': 1.79}, 'NLOS': {'std_dev': 8.93, 'CL': 19.52}},
        'Ka-band': {'LOS': {'std_dev': 1.90}, 'NLOS': {'std_dev': 10.70, 'CL': 29.50}}
    },
    20: {
        'S-band': {'LOS': {'std_dev': 1.14}, 'NLOS': {'std_dev': 9.08, 'CL': 18.17}},
        'Ka-band': {'LOS': {'std_dev': 1.60}, 'NLOS': {'std_dev': 10.00, 'CL': 24.60}}
    },
    30: {
        'S-band': {'LOS': {'std_dev': 1.14}, 'NLOS': {'std_dev': 8.78, 'CL': 18.42}},
        'Ka-band': {'LOS': {'std_dev': 1.90}, 'NLOS': {'std_dev': 11.20, 'CL': 21.90}}
    },
    40: {
        'S-band': {'LOS': {'std_dev': 0.92}, 'NLOS': {'std_dev': 10.25, 'CL': 18.28}},
        'Ka-band': {'LOS': {'std_dev': 2.30}, 'NLOS': {'std_dev': 11.60, 'CL': 20.00}}
    },
    50: {
        'S-band': {'LOS': {'std_dev': 1.42}, 'NLOS': {'std_dev': 10.56, 'CL': 18.63}},
        'Ka-band': {'LOS': {'std_dev': 2.70}, 'NLOS': {'std_dev': 11.80, 'CL': 18.70}}
    },
    60: {
        'S-band': {'LOS': {'std_dev': 1.56}, 'NLOS': {'std_dev': 10.74, 'CL': 17.68}},
        'Ka-band': {'LOS': {'std_dev': 3.10}, 'NLOS': {'std_dev': 10.80, 'CL': 17.80}}
    },
    70: {
        'S-band': {'LOS': {'std_dev': 0.85}, 'NLOS': {'std_dev': 10.17, 'CL': 16.50}},
        'Ka-band': {'LOS': {'std_dev': 3.00}, 'NLOS': {'std_dev': 10.80, 'CL': 17.20}}
    },
    80: {
        'S-band': {'LOS': {'std_dev': 0.72}, 'NLOS': {'std_dev': 11.52, 'CL': 16.30}},
        'Ka-band': {'LOS': {'std_dev': 3.60}, 'NLOS': {'std_dev': 10.80, 'CL': 16.90}}
    },
    90: {
        'S-band': {'LOS': {'std_dev': 0.72}, 'NLOS': {'std_dev': 11.52, 'CL': 16.30}},
        'Ka-band': {'LOS': {'std_dev': 0.40}, 'NLOS': {'std_dev': 10.80, 'CL': 16.80}}
    }
}

def getPlanetRadius(lat):
    """
    Returns the planet radius considering the polar radius and the
    equatorial radius, given the current latitude and longitude.
    """

    cos_lat = math.cos(math.radians(lat))
    sin_lat = math.sin(math.radians(lat))
    Eq_r = 6378000
    Po_r = 6356000
    aux = ((Eq_r**2*cos_lat)**2 + (Po_r**2*sin_lat)**2)
    radius = math.sqrt(aux/((Eq_r*cos_lat)**2 + (Po_r*sin_lat)**2))
    return radius 

def getElevation(sat_coords, terminal_coords):

    a = 6378137
    e = math.sqrt(0.00669437999014)

    long1, lat1, h1 = sat_coords
    x1 = (a/math.sqrt(1 - pow(e,2)*pow(math.sin(math.radians(lat1)),2)) + h1)*math.cos(math.radians(lat1))*math.cos(math.radians(long1))
    y1 = (a/math.sqrt(1 - pow(e,2)*pow(math.sin(math.radians(lat1)),2)) + h1)*math.cos(math.radians(lat1))*math.sin(math.radians(long1))
    z1 = (a*(1 - pow(e,2))/math.sqrt(1 - pow(e,2)*pow(math.sin(math.radians(lat1)),2)) + h1)*math.sin(math.radians(lat1))

    long2, lat2, h2 = terminal_coords
    x2 = (a/math.sqrt(1 - pow(e,2)*pow(math.sin(math.radians(lat2)),2)) + h2)*math.cos(math.radians(lat2))*math.cos(math.radians(long2))
    y2 = (a/math.sqrt(1 - pow(e,2)*pow(math.sin(math.radians(lat2)),2)) + h2)*math.cos(math.radians(lat2))*math.sin(math.radians(long2))
    z2 = (a*(1 - pow(e,2))/math.sqrt(1 - pow(e,2)*pow(math.sin(math.radians(lat2)),2)) + h2)*math.sin(math.radians(lat2))

    Ax = x1 - x2
    Ay = y1 - y2
    Az = z1 - z2

    N = -math.sin(math.radians(lat2)) * math.cos(math.radians(long2)) * Ax - math.sin(math.radians(lat2)) * math.sin(math.radians(long2)) * Ay + math.cos(math.radians(lat2)) * Az
    E = -math.sin(math.radians(long2)) * Ax + math.cos(math.radians(long2)) * Ay
    D = -math.cos(math.radians(lat2)) * math.cos(math.radians(long2)) * Ax - math.cos(math.radians(lat2)) * math.sin(math.radians(long2)) * Ay - math.sin(math.radians(lat2)) * Az

    d = math.sqrt(pow(N,2)+pow(E,2)+pow(D,2))

    elevation = math.degrees(math.asin(-D/d))
        
    return elevation
    
def calculate_distance(sat_coords, terminal_coords):
    """Calculate the distance between the satellite and the terminal."""

    elevation = getElevation(sat_coords, terminal_coords)
    long1, lat1, h1 = sat_coords
    distance = math.sqrt(pow(getPlanetRadius(lat1),2)*pow(math.sin(math.radians(elevation)),2) + pow(h1, 2) + 2*h1*getPlanetRadius(lat1)) - getPlanetRadius(lat1)*math.sin(math.radians(elevation))
    return distance

def free_space_path_loss(distance, frequency):
    """Calculate the free-space path loss in dB."""
    fspl = 32.45 + 20 * math.log10(frequency/1e9) + 20 * math.log10(distance)
    return fspl

def calculate_sf_cl(elevation, frequency, scenario):

    rounded_elevation = round(elevation / 10) * 10

    if scenario not in los_probabilities:
        raise ValueError(f"Scenario {scenario} not recognized.")
    
    if rounded_elevation not in los_probabilities[scenario]:
        raise ValueError(f"Elevation {rounded_elevation} not available.")
    
    los_probability = los_probabilities[scenario][rounded_elevation]

    if frequency > 2e9 and frequency < 4e9:
        band = "S-band"
    elif frequency > 26.5e9 and frequency < 40e9:
        band = "Ka-band"
    else:
        raise ValueError(f"Frequency {frequency} is outside the S and Ka bands.")
    
    std_dev_los = 0
    std_dev_nlos = 0
    cl_nlos = 0
    if scenario == "Dense Urban":
        std_dev_los = dense_urban_sf_cl[rounded_elevation][band]["LOS"]["std_dev"]
        std_dev_nlos = dense_urban_sf_cl[rounded_elevation][band]["NLOS"]["std_dev"]
        cl_nlos = dense_urban_sf_cl[rounded_elevation][band]["NLOS"]["CL"]
    elif scenario == "Urban":
        std_dev_los = urban_sf_cl[rounded_elevation][band]["LOS"]["std_dev"]
        std_dev_nlos = urban_sf_cl[rounded_elevation][band]["NLOS"]["std_dev"]
        cl_nlos = urban_sf_cl[rounded_elevation][band]["NLOS"]["CL"]
    elif scenario == "Suburban":
        std_dev_los = suburban_sf_cl[rounded_elevation][band]["LOS"]["std_dev"]
        std_dev_nlos = suburban_sf_cl[rounded_elevation][band]["NLOS"]["std_dev"]
        cl_nlos = suburban_sf_cl[rounded_elevation][band]["NLOS"]["CL"]

    sf_los = np.random.normal(0, std_dev_los)
    sf_nlos = np.random.normal(0, std_dev_nlos)

    sf_cl = los_probability*sf_los + (1-los_probability)*(sf_nlos+cl_nlos)

    return sf_cl

def plotVisibility(scenario):
    i = 0
    visibility_states = []

    elevations = list(range(15, 166))
    for elevation in elevations:
        rand = np.random.random()
        if elevation<=90:
            elevation_to_round = elevation
        else:
            elevation_to_round = 180 - elevation

        rounded_elevation = round(elevation_to_round / 10) * 10

        los_probability = los_probabilities[scenario][rounded_elevation]
        if rand < los_probability:
            visibility_states.append(1)
        else:
            visibility_states.append(0)

    plt.figure(figsize=(10, 6))
    plt.plot(elevations, visibility_states, label='Visibility (1=LOS, 0=NLOS)', drawstyle='steps-post', color='b')
    plt.yticks([0, 1], ['NLOS', 'LOS'])
    plt.xlabel("Elevation (degrees)")
    plt.ylabel("Visibility State")
    plt.title(f"Visibility State vs Elevation for Scenario: {scenario}")
    plt.grid(True)
    plt.show()

#plotVisibility("Dense Urban")
