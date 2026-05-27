from matplotlib.ticker import MultipleLocator
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import math

transition_probabilities = {
    "Dense Urban": {
        "nlos_to_los": {
            15: 0.003, 25: 0.005, 35: 0.009, 45: 0.012, 55: 0.017,
            65: 0.025, 75: 0.049, 85: 0.194, 90: 1
        },
        "los_to_nlos": {
            90: 0, 95: 0.03, 105: 0.032, 115: 0.034, 125: 0.041,
            135: 0.054, 145: 0.078, 155: 0.121, 165: 0.247
        }
    },
    "Suburban": {
        "nlos_to_los": {
            15: 0.065, 25: 0.042, 35: 0.035, 45: 0.033, 55: 0.035,
            65: 0.043, 75: 0.066, 85: 0.211, 90: 1
        },
        "los_to_nlos": {
            90: 0, 95: 0.001, 105: 0.001, 115: 0.001, 125: 0.001,
            135: 0.002, 145: 0.003, 155: 0.005, 165: 0.016
        }
    }
}

shadow_fading = {
    "std_dev": {
        "LOS": {
            10: 2.3, 20: 1.4, 30: 1.1, 40: 0.9, 50: 0.6,
            60: 0.4, 70: 0.3, 80: 0.3
        },
        "NLOS": {
            10: 4.2, 20: 5.1, 30: 5.6, 40: 6.1, 50: 6.2,
            60: 6.5, 70: 7.1, 80: 7.4
        }
    },
    "decorrelation_distance": {
        "LOS": {
            10: 2.6, 20: 2.8, 30: 2.9, 40: 2.9, 50: 3.0,
            60: 3.1, 70: 3.1, 80: 3.2
        },
        "NLOS": {
            10: 2.5, 20: 3.1, 30: 4.5, 40: 6.4, 50: 8.7,
            60: 10.5, 70: 11.9, 80: 12.6
        }
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

def get_rounded_elevation(elevation):
    rounded_elevation = round(elevation / 10) * 10

    if elevation // rounded_elevation == 1:
        rounded_elevation += 5
    else:
        rounded_elevation -= 5
    
    if rounded_elevation == 175:
        rounded_elevation = 165

    if abs(elevation - 90) < 2.5:
        rounded_elevation = 90
    
    return rounded_elevation

def get_next_visibility_state(elevation, scenario, visibility):
    if scenario not in transition_probabilities:
        raise ValueError(f"Scenario {scenario} not recognized.")
    
    rounded_elevation = get_rounded_elevation(elevation)
    
    rand = np.random.random()

    if visibility == "LOS":
        if rounded_elevation < 90:
            transition_probability = 0
        else:
            transition_probability = transition_probabilities[scenario]["los_to_nlos"][rounded_elevation]
        if rand < transition_probability:
            visibility = "NLOS"
        else:
            visibility = "LOS"
    elif visibility == "NLOS":
        if rounded_elevation > 90:
            transition_probability = 0
        else:
            transition_probability = transition_probabilities[scenario]["nlos_to_los"][rounded_elevation]
        if rand < transition_probability:
            visibility = "LOS"
        else:
            visibility = "NLOS"
    return visibility

def getVisibility(scenario, trials):
    i = 0
    nlos_to_los = []
    los_to_nlos = []
    
    while (i < trials):
        elevations = list(range(10, 171))
        visibility = "NLOS"
        for elevation in elevations:
            last_visibility = visibility
            visibility = get_next_visibility_state(elevation, scenario, visibility)
            if last_visibility == "NLOS" and visibility == "LOS":
                nlos_to_los.append(elevation)
            if last_visibility == "LOS" and visibility == "NLOS":
                los_to_nlos.append(elevation)
        i += 1

    mean_nlos_to_los = 0
    mean_los_to_nlos = 0
    if nlos_to_los:
        mean_nlos_to_los = np.mean(nlos_to_los)
    if los_to_nlos:
        mean_los_to_nlos = np.mean(los_to_nlos)
    
    return mean_nlos_to_los, mean_los_to_nlos

def getCL(distance, frequency, scenario, A_sp, alpha_sp, h_ue):
    c = 3e8
    h_b = 0
    if scenario == "Suburban":
        h_b = 6
    elif scenario == "Dense Urban":
        h_b = 20
    
    reflection_coefficient = 0.3
    ld_180 = 6

    A_sp_rad = np.deg2rad(A_sp)
    
    h_off = h_b - h_ue
    h = h_off*np.sin(A_sp_rad)*np.sqrt(1+np.pow(1/np.tan(alpha_sp),2))
    alpha_1 = np.arctan(h/distance)
    v = np.sqrt(2*frequency/c*distance*alpha_1*A_sp_rad)
    ld = 6.9 +20*np.log10(np.sqrt(np.pow(v-0.1,2)+1)+v-0.1)
    lr = -20*np.log10(reflection_coefficient) + ld_180
    lc = -10*np.log10(np.pow(10,-ld/10)+np.pow(10,-lr/10))

    return lc

def calculate_sf_cl(elevation, distance, frequency, scenario, nlos_to_los_angle, los_to_nlos_angle, h_ue, previous_sf, previous_A_sp):
    if scenario not in transition_probabilities:
        raise ValueError(f"Scenario {scenario} not recognized.")

    if elevation < 90:
        A_sp = np.abs(elevation-nlos_to_los_angle)
    else:
        A_sp = np.abs(elevation-los_to_nlos_angle)

    A_sp_rounded = round(A_sp / 10) * 10
    if A_sp_rounded == 0:
        A_sp_rounded = 10
    
    cl = 0
    std_dev = 0
    decorrelation_distance = 0
    if elevation < nlos_to_los_angle:
        alpha_sp = np.deg2rad(nlos_to_los_angle)
        cl = getCL(distance, frequency, scenario, A_sp, alpha_sp, h_ue)
        std_dev = shadow_fading["std_dev"]["NLOS"][A_sp_rounded]
        decorrelation_distance = shadow_fading["decorrelation_distance"]["NLOS"][A_sp_rounded]
    elif elevation >= nlos_to_los_angle and elevation <= los_to_nlos_angle:
        cl = 0
        std_dev = shadow_fading["std_dev"]["LOS"][A_sp_rounded]
        decorrelation_distance = shadow_fading["decorrelation_distance"]["LOS"][A_sp_rounded]
    elif elevation > los_to_nlos_angle:
        alpha_sp = np.deg2rad(los_to_nlos_angle)
        cl = getCL(distance, frequency, scenario, A_sp, alpha_sp, h_ue)
        std_dev = shadow_fading["std_dev"]["NLOS"][A_sp_rounded]
        decorrelation_distance = shadow_fading["decorrelation_distance"]["NLOS"][A_sp_rounded]

    #Method 1: Compute difference of consecutive angular distances (Very high correlation if the changes between angular distances are low)
    if previous_A_sp is None:
        A_sp_diff = 0
    else:
        A_sp_diff = np.abs(A_sp - previous_A_sp)
    
    if A_sp_diff == 0:
        previous_sf = None

    #Method 2: Compute difference between angular distance and rounded angular distance (Variable correlation)
    #A_sp_diff = np.abs(A_sp - A_sp_rounded)

    #Method 3: Angular distance is used directly (Very low correlation)
    #A_sp_diff = A_sp
    
    rand = np.random.normal(0, std_dev)

    if previous_sf is None:
        sf = rand
    else:
        ro = np.exp(-A_sp_diff/decorrelation_distance)

        sf = ro*previous_sf + np.sqrt(1-np.pow(ro,2))*rand

    return sf, cl, A_sp

def plotVisibility(scenario, trials):
    i = 0
    nlos_to_los = []
    los_to_nlos = []
    
    while (i < trials):
        elevations = list(range(10, 171))
        visibility = "NLOS"
        for elevation in elevations:
            last_visibility = visibility
            visibility = get_next_visibility_state(elevation, scenario, visibility)
            if last_visibility == "NLOS" and visibility == "LOS":
                nlos_to_los.append(elevation)
            if last_visibility == "LOS" and visibility == "NLOS":
                los_to_nlos.append(elevation)
        i += 1

    mean_nlos_to_los = 0
    mean_los_to_nlos = 0
    if nlos_to_los:
        mean_nlos_to_los = np.mean(nlos_to_los)
    if los_to_nlos:
        mean_los_to_nlos = np.mean(los_to_nlos)

    visibility_states = []
    for elevation in elevations:
        if elevation <= mean_nlos_to_los:
            visibility_states.append(0)
        elif mean_nlos_to_los < elevation <= mean_los_to_nlos:
            visibility_states.append(1)
        else:
            visibility_states.append(0)
    plt.figure(figsize=(10, 6))
    plt.plot(elevations, visibility_states, label='Visibility (0=NLOS, 1=LOS)', drawstyle='steps-post', color='b')
    plt.yticks([0, 1], ['NLOS', 'LOS'])
    plt.xlabel("Elevation (degrees)")
    plt.ylabel("Visibility State")
    plt.title(f"Visibility State vs Elevation for Scenario: {scenario}")
    plt.grid(True)
    plt.show()

def plot_cl_sf():
    elevations = list(range(10, 91))
    h1 = 600000
    lat1 = 0
    frequency = 2e9
    scenario = "Suburban"
    nlos_to_los_angle = 50
    los_to_nlos_angle = 130
    h_ue = 0
    pl_1 = []
    pl_2 = []
    cl_1 = []
    angle_1_1 = 30
    angle_1_2 = 150
    cl_2 = []
    angle_2_1 = 50
    angle_2_2 = 130
    cl_3 = []
    angle_3_1 = 70
    angle_3_2 = 110
    sf = None
    A_sp = None
    sf_values = []

    for elevation in elevations:
        previous_sf = sf
        previous_A_sp = A_sp
        distance = math.sqrt(pow(getPlanetRadius(lat1),2)*pow(math.sin(math.radians(elevation)),2) + pow(h1, 2) + 2*h1*getPlanetRadius(lat1)) - getPlanetRadius(lat1)*math.sin(math.radians(elevation))
        sf, cl, A_sp = calculate_sf_cl(elevation, distance, frequency, scenario, angle_1_1, angle_1_2, h_ue, previous_sf, previous_A_sp)
        cl_1.append(cl)
        sf, cl, A_sp = calculate_sf_cl(elevation, distance, frequency, scenario, angle_2_1, angle_2_2, h_ue, previous_sf, previous_A_sp)
        cl_2.append(cl)
        sf, cl, A_sp = calculate_sf_cl(elevation, distance, frequency, scenario, angle_3_1, angle_3_2, h_ue, previous_sf, previous_A_sp)
        cl_3.append(cl)
        sf, cl, A_sp = calculate_sf_cl(elevation, distance, frequency, scenario, nlos_to_los_angle, los_to_nlos_angle, h_ue, previous_sf, previous_A_sp)
        sf_cl = sf + cl
        fspl = free_space_path_loss(distance, frequency)
        pl_1.append(sf_cl + fspl)
        pl_2.append(cl + fspl)
        if sf is not None and sf != 0:
            sf_values.append(sf)

    plt.figure(1)
    plt.plot(elevations, pl_1, linestyle = "-", linewidth = 3, marker = 's', fillstyle='none', markersize = 7, markeredgewidth = 3, color = 'darkblue', label = "FSPL + CL + SF ")
    plt.plot(elevations, pl_2, linestyle = ":", linewidth = 4, color = 'black', label = "FSPL + CL")
    plt.xlabel("Elevation (degrees)")
    plt.ylabel("Path Loss [dB]")
    plt.xlim([min(elevations), max(elevations)])
    plt.ylim([150, 200])
    plt.title(f"Path Loss vs Elevation for {scenario} Scenario")
    plt.legend()
    plt.minorticks_on()
    ax = plt.gca()
    ax.xaxis.set_major_locator(MultipleLocator(10))
    ax.xaxis.set_minor_locator(MultipleLocator(2))
    ax.yaxis.set_major_locator(MultipleLocator(10))
    ax.yaxis.set_minor_locator(MultipleLocator(2))
    plt.grid(which='major', linestyle='-', linewidth=0.5, color='grey')
    plt.grid(which='minor', linestyle=':', linewidth=0.5, color='grey')

    plt.figure(2)
    plt.plot(elevations, cl_1, linestyle = ":", linewidth = 2, marker = 's', fillstyle='none', markersize = 6, markeredgewidth = 3, color = 'grey', label=f"alpha_sp_1 = {angle_1_1}º & alpha_sp_2 = {angle_1_2}º")
    plt.plot(elevations, cl_2, linestyle = ":", linewidth = 2, marker = 'o', fillstyle='none', markersize = 6, markeredgewidth = 3, color = 'black', label=f"alpha_sp_1 = {angle_2_1}º & alpha_sp_2 = {angle_2_2}º")
    plt.plot(elevations, cl_3, linestyle = ":", linewidth = 2, marker = 'x', fillstyle='none', markersize = 6, markeredgewidth = 3, color = 'lightblue', label=f"alpha_sp_1 = {angle_3_1}º & alpha_sp_2 = {angle_3_2}º")
    plt.xlabel("Elevation (degrees)")
    plt.ylabel("Clutter Loss [dB]")
    plt.xlim([min(elevations), max(elevations)])
    plt.ylim([-2, 20])
    plt.title(f"Clutter Loss vs Elevation for {scenario} Scenario")
    plt.legend()
    plt.minorticks_on()
    ax = plt.gca()
    ax.xaxis.set_major_locator(MultipleLocator(10))
    ax.xaxis.set_minor_locator(MultipleLocator(2))
    ax.yaxis.set_major_locator(MultipleLocator(5))
    ax.yaxis.set_minor_locator(MultipleLocator(1))
    plt.grid(which='major', linestyle='-', linewidth=0.5, color='grey')
    plt.grid(which='minor', linestyle=':', linewidth=0.5, color='grey')

    plt.show()

#plotVisibility("Suburban", 10**5)
#print(getVisibility("Suburban"))
#plot_cl_sf()
