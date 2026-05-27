import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from satsim.channel_model import calculate_distance
from satsim.beams import *
from datetime import timedelta
from geodistance import Geodistance
from geopy.distance import geodesic
import math
import numpy as np
import os
import json
import datetime
from collections import defaultdict

def get_data(data_csv, start_time, end_time, time_step_seconds, k):
    # Load the CSV file
    df = pd.read_csv(data_csv, delimiter=';')
    df['Time'] = pd.to_datetime(df['Time'])
    df = df[(df['Time'] >= start_time) & (df['Time'] <= end_time)].set_index('Time')

    ########### Maybe it is not necessary to check the columns    
    # Check for necessary columns
    if 'Possible_to_connect_satellites' not in df.columns or 'RSRP' not in df.columns:
        print("Error: CSV file must contain 'Possible_to_connect_satellites' and 'RSRP' columns.")
        return

    # Convert columns to lists of values
    df['Visible_satellites'] = df['Visible_satellites'].apply(lambda x: x.strip('{}').split(', ') if isinstance(x, str) and x != '{}' else [])
    df['Satellite_positions'] = df['Satellite_positions'].apply(lambda x: [list(map(float, pos.strip('()').split(', '))) for pos in x.strip('{}').split('), (')] if isinstance(x, str) and x != '{}' else [])
    df['Possible_to_connect_satellites'] = df['Possible_to_connect_satellites'].apply(lambda x: x.strip('{}').split(', ') if isinstance(x, str) and x != '{}' else [])
    df['Beams'] = df['Beams'].apply(lambda x: [list(map(int, group.strip('()').split(', '))) for group in x.strip('{}').split('), (')] if isinstance(x, str) and x != '{}' else [])
    df['Beam_centers'] = df['Beam_centers'].apply(lambda x: [[list(map(float, pos.split(', '))) for pos in sat.strip('()').split('), (')] for sat in x.strip('{}').split(')), ((')] if isinstance(x, str) and x != '{}' else [])
    df['RSRP'] = df['RSRP'].apply(lambda x: [list(map(float, group.strip('()').split(', '))) for group in x.strip('{}').split('), (')] if isinstance(x, str) and x != '{}' else [])
    df['Terminal_coords'] = df['Terminal_coords'].apply(lambda x: list(map(float, x.strip('{()}').split(', '))) if isinstance(x, str) and x != '{}' else [])
    
    a = np.pow(1/2,k/4)

    visible_satellites_data = {}
    satellites_data = {}
    terminal_data = {}
    ewa_last = {}

    # Get data for each satellite
    data_by_time = {time: row for time, row in df.iterrows()}
    current_time = start_time

    while current_time <= end_time:
        if current_time in data_by_time:
            row = data_by_time[current_time]
            visible_satellites = row['Visible_satellites']
            satellites_positions_list = row['Satellite_positions']
            satellites = row['Possible_to_connect_satellites']
            beams_list = row['Beams']
            beams_positions_list = row['Beam_centers']
            rsrp_list = row['RSRP']
            terminal_coords = row['Terminal_coords']

            visible_satellites_data[current_time] = {}
            satellites_data[current_time] = {}
            terminal_data[current_time] = {
                "terminal_coords": terminal_coords
            }

            for i, satellite in enumerate(visible_satellites):
                position = satellites_positions_list[i]

                visible_satellites_data[current_time][satellite] = {
                    "position": position
                }

            for i, satellite in enumerate(satellites):
                beams = beams_list[i]
                rsrps = rsrp_list[i]
                positions = beams_positions_list[i]

                if satellite not in ewa_last:
                    ewa_last[satellite] = {}

                ewa_values = np.zeros(len(beams))
                for j, (beam, rsrp) in enumerate(zip(beams, rsrps)):
                    if beam not in ewa_last[satellite]:
                        ewa_last[satellite][beam] = rsrp
                    else:
                        ewa_last[satellite][beam] = (1 - a) * ewa_last[satellite].get(beam, rsrp) + a * rsrp
                    ewa_values[j] = ewa_last[satellite][beam]

                satellites_data[current_time][satellite] = {
                    "beams": beams,
                    "positions": positions,
                    "RSRP": rsrps,
                    "EWA": ewa_values
                }
                
        current_time += timedelta(seconds=time_step_seconds)

    return visible_satellites_data, satellites_data, terminal_data

def get_orbits_data(neighbors_csv):
    try:
        neighbors_df = pd.read_csv(neighbors_csv, delimiter=';')
    except pd.errors.EmptyDataError:
        print(f"Error: The file '{neighbors_csv}' is empty or not properly formatted.")
        return
    except Exception as e:
        print(f"An error occurred while reading '{neighbors_csv}': {e}")
        return

    # Maybe it is not necessary
    satellite_column = 'Satellite'
    orbit_column = 'Orbit'
    if satellite_column not in neighbors_df.columns:
        print(f"Error: Expected column '{satellite_column}' not found in the CSV.")
        return
    if orbit_column not in neighbors_df.columns:
        print(f"Error: Expected column '{orbit_column}' not found in the CSV.")
        return
    
    orbit_to_satellites = {}
    satellite_to_orbit = {}

    for _, row in neighbors_df.iterrows():
        orbit = row['Orbit']
        satellite = row['Satellite']
        if orbit not in orbit_to_satellites:
            orbit_to_satellites[orbit] = []
            index = 0
        orbit_to_satellites[orbit].append(satellite)
        satellite_to_orbit[satellite] = (orbit, index)
        index += 1

    return orbit_to_satellites, satellite_to_orbit

def get_beam(beam_lat, neighbor_sat_coords, rising_neighbor_orbit):
    closest_beam = -1
    threshold = 0.00
    
    beam_centers = oneweb_beam_centers(neighbor_sat_coords)
    
    closest_lat_distance = float('inf')
    for b, beam_center in enumerate(beam_centers):
        neighbor_beam_long, neighbor_beam_lat = beam_center
        lat_distance = abs(neighbor_beam_lat - beam_lat)
        
        if lat_distance < closest_lat_distance and ((not rising_neighbor_orbit and neighbor_beam_lat > beam_lat + threshold) or (rising_neighbor_orbit and neighbor_beam_lat < beam_lat - threshold)):
            closest_lat_distance = lat_distance
            closest_beam = b
    
    return closest_beam

def get_neighboring_orbit_beam(satellites_orbit, visible_satellites_data, rising_neighbor_orbit, beam_lat):
   closest_lat_distance = float('inf')
   closest_satellite = None
   closest_sat_coords = None
   closest_beam = -1


   for satellite in satellites_orbit:
       if satellite in visible_satellites_data:
           neighbor_sat_coords = tuple(map(float, visible_satellites_data[satellite]["position"]))
           lat_distance = abs(neighbor_sat_coords[1] - beam_lat)


           if lat_distance < closest_lat_distance:
               closest_lat_distance = lat_distance
               closest_satellite = satellite
               closest_sat_coords = neighbor_sat_coords


   if closest_satellite is not None:
       closest_beam = get_beam(beam_lat, closest_sat_coords, rising_neighbor_orbit)


   return (closest_satellite, closest_beam)

def neighbor_prediction_oneweb(time, satellites_data, orbit_to_satellites, satellite_to_orbit, num_neighbors, visible_satellites_data, connected_satellite, connected_beam, beam_lat, beam_long, rising, terminal_coords, exact_location, uncertainty):

   def km_per_deg_lon(beam_width_km=1900.0):
       # Kilometers per degree of longitude at a given latitude (WGS84 model)
       a = 6378.137
       b = 6356.7523142
       e2 = (a**2 - b**2) / a**2
       lat_rad = np.deg2rad(beam_lat)
       N = a / np.sqrt(1 - e2 * np.sin(lat_rad)**2)
       km_per_deg_lon = (np.pi / 180.0) * N * np.cos(lat_rad)
      
       return km_per_deg_lon
  
   neighbors = {}


   def calculate_lateral_threshold(beam_lat: float, exact_location_alg) -> float:
       """
       Calculates the lateral threshold percentage based on the latitude using linear interpolation with these two points:
       - At the equator (lat=0º) -> 20%
       - At the north (lat=60º) -> 60%
       """
       
       lat = abs(beam_lat)
      
       if lat <= 60:
           lateral_threshold = 0.2 + (lat / 60) * (0.6 - 0.2) + uncertainty
       else:
           lateral_threshold = 0.6 + (lat - 60) * ((1.0 - 0.6) / 60) + uncertainty

       if lateral_threshold >= 1: exact_location_alg = False

       return lateral_threshold, exact_location_alg


   lateral_threshold_percent, exact_location = calculate_lateral_threshold(beam_lat, exact_location)


   orbit, index = satellite_to_orbit.get(connected_satellite, (None, None))
   if orbit is None:
       return neighbors

   satellites_orbit = orbit_to_satellites[orbit]
  
   next_satellite = connected_satellite
   previous_satellite = connected_satellite


   if exact_location:
    # Considering approximated location of the UE (between 1 and 2 neighbors)

       def calculate_beam_radius(lat_deg, beam_width_equator_km=1900.0):
           radius_deg = 1900 / (2 * 111.0 * np.cos(np.deg2rad(beam_lat)))
           radius_km = radius_deg* km_per_deg_lon()

           return radius_km


       def calculate_lateral_extremes(ue_lat, ue_lon, beam_center_lat, beam_center_lon, lateral_threshold_percent, r_eq_km=950.0):
           """
           Returns (is_left_extreme, is_right_extreme, lateral_distance_km).
           is_left_extreme: True if the UE is on the left side and within the lateral threshold.
           is_right_extreme: True if the UE is on the right side and within the lateral threshold.
           lateral_distance_km: lateral distance (east-west) in km.
           """


           # Lateral distance (small approximation) along the parallel of the beam center
           distance_ue_center_km = (beam_long - terminal_coords[0] ) * km_per_deg_lon()

           beam_radius_km = calculate_beam_radius(beam_lat)

           lateral_threshold_km = beam_radius_km * lateral_threshold_percent

           near_edge = abs(distance_ue_center_km) + lateral_threshold_km >= (beam_radius_km)

           is_left_extreme  = (distance_ue_center_km < 0) and near_edge
           is_right_extreme = (distance_ue_center_km > 0) and near_edge


           return is_left_extreme, is_right_extreme, distance_ue_center_km
      
       # Always add the center following beam
       if connected_beam == 0:
           next_beam = connected_beam + 1
          
           previous_satellite = satellites_orbit[index - 1] if index > 0 else satellites_orbit[len(satellites_orbit) - 1]
           if previous_satellite in visible_satellites_data:
               neighbor_sat_coords = tuple(map(float, visible_satellites_data[previous_satellite]["position"]))
               previous_beam = get_beam(beam_lat, neighbor_sat_coords, rising_neighbor_orbit=True)
           else:
               previous_beam = -1
              
       elif connected_beam == 15:
           next_satellite = satellites_orbit[index + 1] if index < len(satellites_orbit) - 1 else satellites_orbit[0]
           if next_satellite in visible_satellites_data:
               neighbor_sat_coords = tuple(map(float, visible_satellites_data[next_satellite]["position"]))
               next_beam = get_beam(beam_lat, neighbor_sat_coords, rising_neighbor_orbit=False)
           else:
               next_beam = -1

           previous_beam = connected_beam - 1
          
       else:
           next_beam = connected_beam + 1
           previous_beam = connected_beam - 1


       if rising:
           aux_satellite = previous_satellite
           previous_satellite = next_satellite
           next_satellite = aux_satellite


           aux_beam = previous_beam
           previous_beam = next_beam
           next_beam = aux_beam


       neighbors[(next_satellite, next_beam)] = None

       # Check if the UE is at the lateral extremes of the beam
       is_left_extreme, is_right_extreme, lateral_distance = calculate_lateral_extremes(
           terminal_coords[1], terminal_coords[0], beam_lat, beam_long, lateral_threshold_percent)


       # Add neighbors from adjacent orbits if the UE is at the lateral extremes
       if is_left_extreme:
           extreme_next_orbit = False
           if orbit == max(orbit_to_satellites.keys()):
               next_orbit = min(orbit_to_satellites.keys())
               extreme_next_orbit = True
           else:
               next_orbit = orbit + 1
          
           rising_next_orbit = rising if not extreme_next_orbit else not rising
           satellites_next_orbit = orbit_to_satellites[next_orbit]
           next_orbit_neighbor = get_neighboring_orbit_beam(satellites_next_orbit, visible_satellites_data, rising_next_orbit, beam_lat)
           neighbors[(next_orbit_neighbor)] = None

       if is_right_extreme:
           extreme_previous_orbit = False
           if orbit == min(orbit_to_satellites.keys()):
               previous_orbit = max(orbit_to_satellites.keys())
               extreme_previous_orbit = True
           else:
               previous_orbit = orbit - 1


           rising_previous_orbit = rising if not extreme_previous_orbit else not rising
           satellites_previous_orbit = orbit_to_satellites[previous_orbit]
           previous_orbit_neighbor = get_neighboring_orbit_beam(satellites_previous_orbit, visible_satellites_data, rising_previous_orbit, beam_lat)
           neighbors[(previous_orbit_neighbor)] = None

   else:
    # Considering always 1, 2, 3, or 4 neighbors (depending on the num_neighbors chosen)
       if num_neighbors >= 1:
           if connected_beam == 0:
               next_beam = connected_beam + 1
              
               previous_satellite = satellites_orbit[index - 1] if index > 0 else satellites_orbit[len(satellites_orbit) - 1]
               if previous_satellite in visible_satellites_data:
                   neighbor_sat_coords = tuple(map(float, visible_satellites_data[previous_satellite]["position"]))
                   previous_beam = get_beam(beam_lat, neighbor_sat_coords, rising_neighbor_orbit=True)
               else:
                   previous_beam = -1
           elif connected_beam == 15:
               next_satellite = satellites_orbit[index + 1] if index < len(satellites_orbit) - 1 else satellites_orbit[0]
               if next_satellite in visible_satellites_data:
                   neighbor_sat_coords = tuple(map(float, visible_satellites_data[next_satellite]["position"]))
                   next_beam = get_beam(beam_lat, neighbor_sat_coords, rising_neighbor_orbit=False)
               else:
                   next_beam = -1


               previous_beam = connected_beam - 1
           else:
               next_beam = connected_beam + 1
               previous_beam = connected_beam - 1


           if rising:
               aux_satellite = previous_satellite
               previous_satellite = next_satellite
               next_satellite = aux_satellite


               aux_beam = previous_beam
               previous_beam = next_beam
               next_beam = aux_beam
      
           neighbors[(next_satellite, next_beam)] = None


       if num_neighbors >= 2:
           extreme_next_orbit = False
           if orbit == max(orbit_to_satellites.keys()):
               next_orbit = min(orbit_to_satellites.keys())
               extreme_next_orbit = True
           else:
               next_orbit = orbit + 1
          
           rising_next_orbit = rising if not extreme_next_orbit else not rising
           satellites_next_orbit = orbit_to_satellites[next_orbit]
           next_orbit_neighbor = get_neighboring_orbit_beam(satellites_next_orbit, visible_satellites_data, rising_next_orbit, beam_lat)
           neighbors[(next_orbit_neighbor)] = None


       if num_neighbors >= 3:
           extreme_previous_orbit = False
           if orbit == min(orbit_to_satellites.keys()):
               previous_orbit = max(orbit_to_satellites.keys())
               extreme_previous_orbit = True
           else:
               previous_orbit = orbit - 1


           rising_previous_orbit = rising if not extreme_previous_orbit else not rising
           satellites_previous_orbit = orbit_to_satellites[previous_orbit]
           previous_orbit_neighbor = get_neighboring_orbit_beam(satellites_previous_orbit, visible_satellites_data, rising_previous_orbit, beam_lat)
           neighbors[(previous_orbit_neighbor)] = None
      
       if num_neighbors == 4:
           neighbors[(previous_satellite, previous_beam)] = None

   return neighbors



def sib19_construction_algorithm(visible_satellites_data, satellites_data, orbit_to_satellites, satellite_to_orbit, num_neighbors, previous_sib19, connected_satellite, connected_beam, time, time_step, terminal_coords, uncertainty):
    sib19 = {}
    
    if time not in satellites_data:
        return sib19

    data = False
    beam_position = None
    if connected_satellite in satellites_data[time] and connected_beam in satellites_data[time][connected_satellite]["beams"]:
        beam_index = satellites_data[time][connected_satellite]["beams"].index(connected_beam)
        beam_position = satellites_data[time][connected_satellite]["positions"][beam_index]
        data = True

    sib19[(connected_satellite, connected_beam)] = beam_position

    if not previous_sib19:
        if data:
            beam_index = satellites_data[time][connected_satellite]["beams"].index(connected_beam)
            beam_lat = float(satellites_data[time][connected_satellite]["positions"][beam_index][1])
            beam_long = float(satellites_data[time][connected_satellite]["positions"][beam_index][0])
            lat_previous_step = visible_satellites_data[time - time_step][connected_satellite]["position"][1]
            lat = visible_satellites_data[time][connected_satellite]["position"][1]

            if lat < lat_previous_step:
                rising = False
            else:
                rising = True

            neighbors = neighbor_prediction_oneweb(time, satellites_data, orbit_to_satellites, satellite_to_orbit, num_neighbors, visible_satellites_data[time], connected_satellite, connected_beam, beam_lat, beam_long, rising, terminal_coords, True, uncertainty)
        else:
            neighbors = {}
    else:
        neighbors = previous_sib19

    for satellite, beam in neighbors:

        beam_position = None

        if satellite in satellites_data[time] and beam in satellites_data[time][satellite]["beams"]:
            beam_index = satellites_data[time][satellite]["beams"].index(beam)
            beam_position = satellites_data[time][satellite]["positions"][beam_index]

        sib19[(satellite, beam)] = beam_position

    return sib19

def cho_programming_algorithm(sib19, strategy, a3_offset, a3_hys, a4_threshold, a4_hys, d2_threshold1, d2_threshold2, d2_hys, ttt):
    rrc_reconfiguration = {}

    for i, (satellite, beam) in enumerate(sib19):
        if i != 0:
            if strategy == "Signal":
                rrc_reconfiguration[(satellite, beam)] = {
                    "A3 Offset": a3_offset,
                    "A3 Hysteresis": a3_hys,
                    "TTT": ttt
                }
            elif strategy == "Distance":
                if i == 1 or i == 4:
                    neighbor = 0
                else:
                    neighbor = 1

                rrc_reconfiguration[(satellite, beam)] = {
                    "A4 Threshold": a4_threshold,
                    "A4 Hysteresis": a4_hys,
                    "D2 Threshold 1": d2_threshold1[neighbor],
                    "D2 Threshold 2": d2_threshold2[neighbor],
                    "D2 Hysteresis": d2_hys,
                    "TTT": ttt
                }

    return rrc_reconfiguration

def ho_decision_algorithm(a3_offset, a3_hys, rsrp_source, rsrp_target, aux_time):
    handover_decision = True
    
    #A3-1: Entering Condition
    a3_1 = rsrp_target - a3_hys > rsrp_source + a3_offset
    #A3-2: Leaving Condition
    a3_2 = rsrp_target + a3_hys < rsrp_source + a3_offset

    if a3_2:
        handover_decision = False
    if aux_time == 0 and not a3_1:
        handover_decision = False
    
    return handover_decision

def handover_strategies(data_csv, neighbors_csv, start_time, end_time, time_step_seconds, ue_connectivity_window, constellation, strategy, handover_type, a3_offset, a3_hys, a4_threshold, a4_hys, d2_threshold1, d2_threshold2, d2_hys, k, ttt, num_neighbors, uncertainty):
    """
    Plots satellites that the UE is connected to over time. UE connects to the closest one
    
    Parameters:
    - data_csv: The CSV file containing the visibility and connection data.
    - neighbors_csv: The CSV file containing the neigbors of the satellites
    - start_time: The start time of the plot.
    - end_time: The end time of the plot.
    - time_step_seconds: The time step in seconds for the plot.
    - ue_connectivity_window: Vector of 1s and 0s indicating if the UE is in connected or sleeping mode.
    - signal_threshold: Signal-based threshold triggering the handover.
    - k: L3 Filter Coefficient that affects the Exponential Weighted Average (EWA).
    """

    visible_satellites_data, satellites_data, terminal_data = get_data(data_csv, start_time, end_time, time_step_seconds, k)
    
    orbit_to_satellites = {}
    satellite_to_orbit = {}
    if constellation == "OneWeb":
        orbit_to_satellites, satellite_to_orbit = get_orbits_data(neighbors_csv)

    results = {
        "RSRP": {
            "times": [],
            "EWA": []
        },
        "ToS": [],
        "Effective Time": [],
        "Events": {
            "Programmed Intra-gNB CHO": 0,
            "Programmed Inter-gNB CHO": 0,
            "Executed Intra-gNB HO": 0,
            "Executed Inter-gNB HO": 0,
        },
        "Ping-Pong": 0,
        "Call Drops": 0,
        "Distances": {
            "times": [],
            "distances": []
        },
        "Distances HO": [],
        "Positions Call Drop": [],
        "Call Drop Cause": {
            "Call Drop and Event not Triggered yet": 0,
            "Call Drop during HO Preparation": 0,
            "Call Drop during HO Execution": 0
        }
    }

    events_data = []

    preparation_time = 0.2
    execution_time = 0.2
    reconnection_time = 60
    
    connectivity_index = 0
    handover = False
    next_satellite = None
    next_beam = -1
    current_time = start_time
    previous_satellite = None
    previous_beam = -1
    connected_satellite = None
    connected_beam = -1
    time_step = timedelta(seconds=time_step_seconds)
    while current_time < end_time:

        if (not(ue_connectivity_window) or ue_connectivity_window[connectivity_index] == 1):
            
            rsrp_max = -float('inf')
            triggering_condition = False
            tos = 0
            effective_time = 0
            time = current_time
            connectivity_index_2 = connectivity_index
            connection_over = False
            sib19 = {}

            if time in satellites_data:
                potential_satellites = satellites_data[time].keys()
            
                if next_satellite is None or next_beam == -1:
                    for satellite in potential_satellites:
                        ewa_values = satellites_data[time][satellite]["EWA"]
                        max_index = np.argmax(ewa_values)
                        rsrp = ewa_values[max_index]
                        if rsrp > rsrp_max:
                            rsrp_max = rsrp
                            next_satellite = satellite
                            next_beam = satellites_data[time][satellite]["beams"][max_index]

                previous_satellite = connected_satellite
                previous_beam = connected_beam
                connected_satellite = next_satellite
                connected_beam = next_beam
                next_satellite = None
                next_beam = -1
                if connected_satellite is not None and connected_beam != -1 and not(handover):
                    results['RSRP']['times'].append(time)
                    results['RSRP']['EWA'].append(rsrp_max)
                    
                    beam_index = satellites_data[time][connected_satellite]["beams"].index(connected_beam)
                    beam_position = satellites_data[time][connected_satellite]["positions"][beam_index]
                    long, lat = map(float, beam_position)

                    terminal_coords = terminal_data[time]["terminal_coords"]
                    distance_ue_cell = geodesic((lat, long), (terminal_coords[1], terminal_coords[0])).m
                    results['Distances']['times'].append(time)
                    results['Distances']['distances'].append(distance_ue_cell)

            execution_time_counted = False
            ho_executed = False
            time_execution_finished = time

            preparation_time_executed = False
            time_preparation_finished = time

            execution_fail = False

            new_connection = True

            call_drop_reason = None
            
            while not(triggering_condition) and time < end_time and connected_satellite is not None and connected_beam != -1 and not(connection_over) and (not(ue_connectivity_window) or ue_connectivity_window[connectivity_index_2] == 1):
                time_advance = False
                
                #Execution time
                if handover and not execution_time_counted:
                    execution_time_counted = True
                    time_advance = True
                    call_drop_reason = "Call Drop during HO Execution"
                    i = 0
                    while i < execution_time and time < end_time:
                        i = round(i + time_step_seconds, 10)
                        time += time_step
                        sib19 = sib19_construction_algorithm(visible_satellites_data, satellites_data, orbit_to_satellites, satellite_to_orbit, num_neighbors, sib19, connected_satellite, connected_beam, time, time_step, terminal_data[time]["terminal_coords"], uncertainty)

                        connectivity_ok = not ue_connectivity_window or ue_connectivity_window[connectivity_index_2 + int(i / time_step_seconds)] == 1

                        if not connectivity_ok or not sib19:
                            execution_fail = True
                            break

                        if not(sib19 and (connected_satellite, connected_beam) in sib19 and sib19[(connected_satellite, connected_beam)] is not None):
                            execution_fail = True
                            break

                        tos = round(tos + time_step_seconds, 10)
                        
                    connectivity_index_2 += int(i/time_step_seconds)
                    time_execution_finished = time
                
                #Preparation time (CHO)
                if handover_type == "CHO" and not preparation_time_executed and not execution_fail:
                    preparation_time_executed = True
                    time_advance = True
                    call_drop_reason = "Call Drop during HO Preparation"
                    i = 0
                    while i < preparation_time and time < end_time:
                        time += time_step
                        sib19 = sib19_construction_algorithm(visible_satellites_data, satellites_data, orbit_to_satellites, satellite_to_orbit, num_neighbors, sib19, connected_satellite, connected_beam, time, time_step, terminal_data[time]["terminal_coords"], uncertainty)
                        i = round(i + time_step_seconds, 10)
                        
                        connectivity_ok = not ue_connectivity_window or ue_connectivity_window[connectivity_index_2 + int(i / time_step_seconds)] == 1

                        if not connectivity_ok or not sib19:
                            break

                        if not(sib19 and (connected_satellite, connected_beam) in sib19 and sib19[(connected_satellite, connected_beam)] is not None):
                            break

                        tos = round(tos + time_step_seconds, 10)
                        effective_time = round(effective_time + time_step_seconds, 10)

                    connectivity_index_2 += int(i/time_step_seconds)
                    time_preparation_finished = time

                    rrc_reconfiguration = cho_programming_algorithm(sib19, strategy, a3_offset, a3_hys, a4_threshold, a4_hys, d2_threshold1, d2_threshold2, d2_hys, ttt)
                    
                    for satellite, beam in rrc_reconfiguration:
                        if satellite == connected_satellite:
                            results['Events']['Programmed Intra-gNB CHO'] += 1
                        else:
                            results['Events']['Programmed Inter-gNB CHO'] += 1
                
                if not time_advance and not execution_fail:
                    call_drop_reason = "Call Drop and Event not Triggered yet"
                    time += time_step
                    sib19 = sib19_construction_algorithm(visible_satellites_data, satellites_data, orbit_to_satellites, satellite_to_orbit, num_neighbors, sib19, connected_satellite, connected_beam, time, time_step, terminal_data[time]["terminal_coords"], uncertainty)

                if not(handover) and tos > 0 and new_connection:
                    new_connection = False
                    events_data.append({
                        "time": current_time,
                        "event": "New Connection",
                        "satellite": connected_satellite,
                        "beam": connected_beam
                    })

                if sib19 and (connected_satellite, connected_beam) in sib19 and sib19[(connected_satellite, connected_beam)] is not None:                    
                    rsrp_max = -float('inf')
                    call_drop_reason = "Call Drop and Event not Triggered yet"

                    if handover_type == "CHO":
                        tracked_satellites = rrc_reconfiguration
                    else:
                        tracked_satellites = sib19

                    for satellite, beam in tracked_satellites:
                        if handover_type == "CHO":
                            time_to_trigger = tracked_satellites[(satellite, beam)]["TTT"]
                        else:
                            time_to_trigger = ttt

                        aux_time = 0
                        check_time = time
                        sib19 = sib19_construction_algorithm(visible_satellites_data, satellites_data, orbit_to_satellites, satellite_to_orbit, num_neighbors, sib19, connected_satellite, connected_beam, check_time, time_step, terminal_data[check_time]["terminal_coords"], uncertainty)

                        while aux_time <= time_to_trigger and check_time < end_time:
                            if satellite == connected_satellite and beam == connected_beam:
                                break

                            connectivity_ok = not ue_connectivity_window or ue_connectivity_window[connectivity_index_2 + int(aux_time / time_step_seconds)] == 1
                            
                            if not connectivity_ok or not sib19:
                                break
                            
                            if (connected_satellite, connected_beam) not in sib19 or sib19[(connected_satellite, connected_beam)] is None:
                                break

                            if (satellite, beam) not in sib19 or sib19[(satellite, beam)] is None:
                                break

                            source_data = satellites_data[check_time][connected_satellite]
                            target_data = satellites_data[check_time][satellite]

                            beam_index_source = source_data["beams"].index(connected_beam)
                            beam_index_target = target_data["beams"].index(beam)

                            rsrp_target = target_data["EWA"][beam_index_target]

                            if strategy == "Signal":
                                rsrp_source = source_data["EWA"][beam_index_source]

                                if handover_type == "CHO":
                                    offset_a3 = tracked_satellites[(satellite, beam)]["A3 Offset"]
                                    hys_a3 = tracked_satellites[(satellite, beam)]["A3 Hysteresis"]

                                    #A3-1: Entering Condition
                                    a3_1 = rsrp_target - hys_a3 > rsrp_source + offset_a3
                                    #A3-2: Leaving Condition
                                    a3_2 = rsrp_target + hys_a3 < rsrp_source + offset_a3

                                    if a3_2:
                                        break
                                    if aux_time == 0 and not a3_1:
                                        break
                                else:
                                    handover_decision = ho_decision_algorithm(a3_offset, a3_hys, rsrp_source, rsrp_target, aux_time)
                                    if not handover_decision:
                                        break
                                
                                if rsrp_target > rsrp_max and aux_time == time_to_trigger:
                                    triggering_condition = True
                                    rsrp_max = rsrp_target
                                    next_satellite = satellite
                                    next_beam = beam

                            elif strategy == "Distance" and handover_type == "CHO":
                                terminal_coords = terminal_data[check_time]["terminal_coords"]

                                source_position = sib19[(connected_satellite, connected_beam)]
                                long_source, lat_source = map(float, source_position)
                                distance_source = geodesic((lat_source, long_source), (terminal_coords[1], terminal_coords[0])).m

                                target_position = sib19[(satellite, beam)]
                                long_target, lat_target = map(float, target_position)
                                distance_target = geodesic((lat_target, long_target), (terminal_coords[1], terminal_coords[0])).m
                                
                                threshold_a4 = tracked_satellites[(satellite, beam)]["A4 Threshold"]
                                hys_a4 = tracked_satellites[(satellite, beam)]["A4 Hysteresis"]
                                d2_threshold_1 = tracked_satellites[(satellite, beam)]["D2 Threshold 1"]
                                d2_threshold_2 = tracked_satellites[(satellite, beam)]["D2 Threshold 2"]
                                hys_d2 = tracked_satellites[(satellite, beam)]["D2 Hysteresis"]

                                #A4-1: Entering Condition
                                a4_1 = rsrp_target - hys_a4 > threshold_a4
                                #A4-2: Leaving Condition
                                a4_2 = rsrp_target + hys_a4 < threshold_a4
                                #D2-1: Entering Condition 1
                                d2_1 = distance_source - hys_d2 > d2_threshold_1
                                #D2-2: Entering Condition 2
                                d2_2 = distance_target + hys_d2 < d2_threshold_2
                                #D2-3: Leaving Condition 1
                                d2_3 = distance_source + hys_d2 < d2_threshold_1
                                #D2-4: Leaving Condition 2
                                d2_4 = distance_target - hys_d2 > d2_threshold_2

                                if a4_2 or d2_3 or d2_4:
                                    break
                                if aux_time == 0 and not(a4_1 and d2_1 and d2_2):
                                    break
                                if rsrp_target > rsrp_max and aux_time == time_to_trigger:
                                    triggering_condition = True
                                    rsrp_max = rsrp_target
                                    next_satellite = satellite
                                    next_beam = beam

                            else:
                                break
                            
                            aux_time = round(aux_time + time_step_seconds, 10)
                            check_time += time_step
                            sib19 = sib19_construction_algorithm(visible_satellites_data, satellites_data, orbit_to_satellites, satellite_to_orbit, num_neighbors, sib19, connected_satellite, connected_beam, check_time, time_step, terminal_data[check_time]["terminal_coords"], uncertainty)

                    i = 0
                    end = False
                    check_time = time
                    sib19 = sib19_construction_algorithm(visible_satellites_data, satellites_data, orbit_to_satellites, satellite_to_orbit, num_neighbors, sib19, connected_satellite, connected_beam, check_time, time_step, terminal_data[check_time]["terminal_coords"], uncertainty)

                    while i <= ttt and not(end) and check_time <= end_time:
                        connectivity_ok = not ue_connectivity_window or ue_connectivity_window[connectivity_index_2 + int(i / time_step_seconds)] == 1

                        if not connectivity_ok or not sib19:
                            break

                        if i != ttt or not(triggering_condition):
                            selected_satellite = connected_satellite
                            selected_beam = connected_beam
                            if not(triggering_condition):
                                end = True
                        else:
                            selected_satellite = next_satellite
                            selected_beam = next_beam

                        if ue_connectivity_window and connectivity_index_2 + int(aux_time/time_step_seconds) != len(ue_connectivity_window) - 1:
                            if ue_connectivity_window[connectivity_index_2 + int(i/time_step_seconds) + 1] == 0:
                                selected_satellite = connected_satellite
                                selected_beam = connected_beam
                            
                        if (selected_satellite, selected_beam) not in sib19 or sib19[(selected_satellite, selected_beam)] is None:
                            break

                        sat_data = satellites_data[check_time][selected_satellite]

                        index = sat_data["beams"].index(selected_beam)
                        
                        results['RSRP']['times'].append(check_time)

                        rsrp = sat_data["EWA"][index]
                        results['RSRP']['EWA'].append(rsrp)

                        terminal_coords = terminal_data[check_time]["terminal_coords"]
                        position = sib19[(selected_satellite, selected_beam)]
                        long, lat = map(float, position)
                        distance_ue_cell = geodesic((lat, long), (terminal_coords[1], terminal_coords[0])).m

                        results['Distances']['times'].append(check_time)
                        results['Distances']['distances'].append(distance_ue_cell)
                        
                        if execution_time and preparation_time == 0:
                            tos = round(tos + time_step_seconds, 10)
                            effective_time = round(effective_time + time_step_seconds, 10)
                        elif triggering_condition and time_execution_finished != check_time and time_preparation_finished != check_time:
                            tos = round(tos + time_step_seconds, 10)
                            effective_time = round(effective_time + time_step_seconds, 10)
                        elif not triggering_condition and time_execution_finished != check_time and time_preparation_finished != check_time:
                            tos = round(tos + time_step_seconds, 10)
                            effective_time = round(effective_time + time_step_seconds, 10)
                        
                        i = round(i + time_step_seconds, 10)
                        if i <= ttt and not end:
                            check_time += time_step
                            sib19 = sib19_construction_algorithm(visible_satellites_data, satellites_data, orbit_to_satellites, satellite_to_orbit, num_neighbors, sib19, connected_satellite, connected_beam, check_time, time_step, terminal_data[check_time]["terminal_coords"], uncertainty)

                    i = round(i - time_step_seconds, 10)
                    connectivity_index_2 += int(i / time_step_seconds)
                    time = check_time
                    sib19 = sib19_construction_algorithm(visible_satellites_data, satellites_data, orbit_to_satellites, satellite_to_orbit, num_neighbors, sib19, connected_satellite, connected_beam, time, time_step, terminal_data[time]["terminal_coords"], uncertainty)
                    
                    #Preparation time (HO)
                    if triggering_condition and handover_type == "HO":
                        call_drop_reason = "Call Drop during HO Preparation"
                        i = 0
                        while i < preparation_time and time < end_time:
                            i = round(i + time_step_seconds, 10)
                            time += time_step
                            sib19 = sib19_construction_algorithm(visible_satellites_data, satellites_data, orbit_to_satellites, satellite_to_orbit, num_neighbors, sib19, connected_satellite, connected_beam, time, time_step, terminal_data[time]["terminal_coords"], uncertainty)


                            connectivity_ok = not ue_connectivity_window or ue_connectivity_window[connectivity_index_2 + int(i / time_step_seconds)] == 1

                            if not connectivity_ok or time not in satellites_data:
                                connection_over = True
                                break

                            if not(time in satellites_data and connected_satellite in satellites_data[time] and connected_beam in satellites_data[time][connected_satellite]["beams"]):
                                connection_over = True
                                break

                            tos = round(tos + time_step_seconds, 10)
                            effective_time = round(effective_time + time_step_seconds, 10)

                        connectivity_index_2 += int(i/time_step_seconds)

                    #Execution time
                    if triggering_condition and not connection_over:
                        ho_executed = True
                        i = 0
                        while i < execution_time and time < end_time:
                            i = round(i + time_step_seconds, 10)
                            time += time_step
                            sib19 = sib19_construction_algorithm(visible_satellites_data, satellites_data, orbit_to_satellites, satellite_to_orbit, num_neighbors, sib19, connected_satellite, connected_beam, time, time_step, terminal_data[time]["terminal_coords"], uncertainty)
                        
                            connectivity_ok = not ue_connectivity_window or ue_connectivity_window[connectivity_index_2 + int(i / time_step_seconds)] == 1

                            if not connectivity_ok or not sib19:
                                ho_executed = False
                                break

                            if not(sib19 and (next_satellite, next_beam) in sib19 and sib19[(next_satellite, next_beam)] is not None):
                                ho_executed = False
                                break
                            
                        connectivity_index_2 += int(i/time_step_seconds)

                else:
                    connection_over = True
                    next_satellite = None
                    next_beam = -1

        else:
            connected_satellite = None
            connected_beam = -1

        final_time = current_time + time_step
        if connected_satellite is not None and connected_beam != -1:
            tos_equal_connected = False
            if not(ue_connectivity_window):
                time_to_represent = tos
            else:
                connected_time_resting = 0
                index_3 = connectivity_index
                while index_3 < len(ue_connectivity_window) and ue_connectivity_window[index_3] == 1:
                    connected_time_resting += time_step_seconds
                    index_3 += 1
                if connected_time_resting == tos:
                    tos_equal_connected = True
                time_to_represent = min(connected_time_resting, tos)
            if time_to_represent > 0:
                final_time = current_time + timedelta(seconds=time_to_represent)
                results['ToS'].append(time_to_represent)
                results['Effective Time'].append(effective_time)
                event = None
                if triggering_condition and time_to_represent == tos and not(tos_equal_connected) and not(connection_over) and final_time < end_time:
                    if ho_executed:
                        if connected_satellite != next_satellite:
                            event = 'Inter-gNB'
                        else:
                            event = 'Intra-gNB'
                    
                        results['Events'][f'Executed {event} HO'] += 1
                        
                        results['Distances HO'].append(distance_ue_cell)
                    else:
                        event = 'None'
                    
                    handover = True

                else:
                    if final_time < end_time:
                        if time_to_represent == tos and not(tos_equal_connected):
                            event = 'Call Drop'
                            results['Call Drops'] += 1
                            results['Call Drop Cause'][call_drop_reason] += 1

                            terminal_coords = terminal_data[final_time]["terminal_coords"]
                            sat_data = satellites_data[final_time][connected_satellite]
                            index = sat_data["beams"].index(connected_beam)
                            position = sat_data["positions"][index]
                            long_beam, lat_beam = map(float, position)
                            
                            long_diff = terminal_coords[0] - long_beam
                            lat_diff = terminal_coords[1] - lat_beam

                            results['Positions Call Drop'].append((long_diff, lat_diff, lat_beam))
                        else:
                            event = 'UE Disconnection'
                        next_satellite = None
                        next_beam = -1
                    else:
                        event = 'End'
                    handover = False

                events_data.append({
                    "time": final_time,
                    "event": event,
                    "satellite": connected_satellite,
                    "beam": connected_beam
                })
                if connection_over:
                    time_to_represent += reconnection_time
                    final_time += timedelta(seconds=reconnection_time)
                current_time = final_time
                connectivity_index += round(time_to_represent/time_step_seconds)
            else:
                if handover:
                    results['Call Drops'] += 1
                    results['Call Drop Cause'][call_drop_reason] += 1
                    events_data.append({
                        "time": current_time,
                        "event": "Call Drop",
                        "satellite": previous_satellite,
                        "beam": previous_beam
                    })
                    
                    terminal_coords = terminal_data[current_time]["terminal_coords"]
                    sat_data = satellites_data[current_time][previous_satellite]
                    index = sat_data["beams"].index(previous_beam)
                    position = sat_data["positions"][index]
                    long_beam, lat_beam = map(float, position)
                    
                    long_diff = terminal_coords[0] - long_beam
                    lat_diff = terminal_coords[1] - lat_beam

                    results['Positions Call Drop'].append((long_diff, lat_diff, lat_beam))

                    current_time = final_time + timedelta(seconds=reconnection_time)
                    connectivity_index += round(reconnection_time/time_step_seconds)
                else:
                    current_time = final_time
                    connectivity_index += 1
                handover = False
            if time_to_represent > 0 and time_to_represent < 1 and previous_beam != -1 and next_beam != -1:
                if previous_beam == next_beam and previous_satellite == next_satellite:
                    results['Ping-Pong'] += 1
        else:
            next_satellite = None
            next_beam = -1
            current_time = final_time
            connectivity_index += 1
            handover = False
    
    return results, events_data, satellites_data

def plot_connectivity(events_data, satellites_data, start_time, end_time, time_step_seconds, ue_connectivity_window, i, data_path):
    
    all_satellites = set()
    for time_data in satellites_data.values():
        all_satellites.update(time_data.keys())
    all_satellites = sorted(all_satellites)
    all_satellites_map = {sat: id for id, sat in enumerate(all_satellites)}

    # Define y-ticks and labels
    y_ticks = list(range(len(all_satellites)))
    y_labels = all_satellites
    
    # Prepare the figure and axis
    fig, ax = plt.subplots(figsize=(12, 8))

    current_time = start_time
    connectivity_index = 0
    
    while current_time < end_time:
        next_time = current_time + timedelta(seconds=time_step_seconds)

        possible_connections_now = set(satellites_data.get(current_time, {}).keys())
        possible_connections_next = set(satellites_data.get(next_time, {}).keys())

        for satellite in sorted(possible_connections_now and possible_connections_next):
            ax.hlines(y=y_ticks[all_satellites_map[satellite]], xmin=current_time, xmax=next_time, color='purple', linewidth=3)
            if ue_connectivity_window and ue_connectivity_window[connectivity_index] == 0:
                ax.axvspan(current_time, next_time, color='grey', alpha=1)

        current_time = next_time
        connectivity_index += 1
    
    previous_event = None
    previous_time = None
    previous_satellite = None
    for event in events_data:
        time = event['time']
        event_type = event['event']
        connected_satellite = event['satellite']

        if previous_satellite is None:
            previous_satellite = connected_satellite

        if previous_time is not None and previous_time != time:
            connected_beam = event['beam']
            ax.text(previous_time + (time - previous_time) / 2, y_ticks[all_satellites_map[previous_satellite]], f"{connected_beam}", color='red', fontsize=10, ha='center', va='bottom')
            ax.hlines(y=y_ticks[all_satellites_map[previous_satellite]], xmin=previous_time, xmax=time, color='red', linewidth=3)

        if previous_event == 'Intra-gNB':
            ax.plot(previous_time, y_ticks[all_satellites_map[connected_satellite]], 'o', color='yellow', markersize=4)
            previous_event = None
        elif previous_event == 'Inter-gNB':
            ax.plot(previous_time, y_ticks[all_satellites_map[connected_satellite]], 'o', color='blue', markersize=4)
            previous_event = None

        previous_time = time
        previous_satellite = connected_satellite

        if event_type == 'New Connection':
            ax.plot(time, y_ticks[all_satellites_map[connected_satellite]], 'o', color='red', markersize=4)
        elif event_type == 'Intra-gNB':
            previous_event = 'Intra-gNB'
        elif event_type == 'Inter-gNB':
            ax.plot(time, y_ticks[all_satellites_map[connected_satellite]], 'o', color='green', markersize=4)
            previous_event = 'Inter-gNB'
            previous_satellite = None
        elif event_type == 'Call Drop':
            ax.plot(time, y_ticks[all_satellites_map[connected_satellite]], 'o', color='black', markersize=4)
            previous_time = None
        elif event_type == 'UE Disconnection':
            ax.plot(time, y_ticks[all_satellites_map[connected_satellite]], 'x', color='black', markersize=4)
            previous_time = None
        elif event_type == 'None':
            previous_event = 'None'
            previous_satellite = None

    # Set y-ticks
    ax.set_yticks(y_ticks)
    ax.set_yticklabels(y_labels)
    
    # Labels and title
    plt.xlabel('Time', fontsize=14)
    plt.ylabel('Satellites', fontsize=14)
    plt.title('Satellite Visibility and Connection Over Time', fontsize=16)
    
    # Set x-axis limits
    ax.set_xlim([start_time, end_time])

    ax.tick_params(axis='both', which='major', labelsize=12)

    # Format the x-axis to show only time
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M:%S.%f'))
    
    # Rotate x-axis labels for better readability
    plt.xticks(rotation=45)

    # Grid and layout
    plt.grid(True)
    plt.tight_layout()
    
    os.makedirs(data_path, exist_ok=True)
    filename = os.path.join(data_path, f"connectivity_north_false{i}.png")

    plt.show()
    plt.savefig(filename, dpi=300, bbox_inches='tight')
    plt.close()


def plot_satellite_connectivity(data_csv, start_time, end_time, time_step_seconds, ue_connectivity_window):
    """
    Plots the connectivity of the satellites along the simulation
    
    Parameters:
    - data_csv: The CSV file containing the visibility and connection data.
    - start_time: The start time of the plot.
    - end_time: The end time of the plot.
    - time_step_seconds: The time step in seconds for the plot.
    - ue_connectivity_window: Vector of 1s and 0s indicating if the UE is in connected or sleeping mode.
    """
    try:
        # Load the visibility data from CSV (handle semicolon delimiter)
        visibility_df = pd.read_csv(data_csv, delimiter=';')
        visibility_df['Time'] = pd.to_datetime(visibility_df['Time'])
    except pd.errors.EmptyDataError:
        print(f"Error: The file '{data_csv}' is empty or not properly formatted.")
        return
    except Exception as e:
        print(f"An error occurred while reading '{data_csv}': {e}")
        return
    
    #############################################
    ###### Maybe it is not necessary to check this, as everything will probably work fine
    ###### Nevertheless, it can be better to check the columns so that errors can be found

    # Check and handle the expected columns
    visibility_column = 'Visible_satellites'
    position_column = 'Satellite_positions'  # For parsing positions if needed
    if visibility_column not in visibility_df.columns:
        print(f"Error: Expected column '{visibility_column}' not found in the CSV.")
        return
    
    all_satellites = set(visibility_df['Possible_to_connect_satellites'].str.strip('{}').str.split(', ').explode().unique())
    all_satellites = sorted(all_satellites)
    all_satellites_map = {sat: id for id, sat in enumerate(all_satellites)}

    # Define y-ticks and labels
    y_ticks = list(range(len(all_satellites)))
    y_labels = all_satellites
    
    # Prepare the figure and axis
    fig, ax = plt.subplots(figsize=(12, 8))

    current_time = start_time
    connectivity_index = 0
    
    while current_time < end_time:
        next_time = current_time + timedelta(seconds=time_step_seconds)

        connections_now = visibility_df.loc[visibility_df['Time'] == current_time, 'Possible_to_connect_satellites']
        connections_next = visibility_df.loc[visibility_df['Time'] == next_time, 'Possible_to_connect_satellites']
        
        possible_connections_now = set(connections_now.str.strip('{}').str.split(', ').explode().dropna())
        possible_connections_next = set(connections_next.str.strip('{}').str.split(', ').explode().dropna())
        possible_connections_now = set(filter(None, possible_connections_now))
        possible_connections_next = set(filter(None, possible_connections_next))

        for satellite in sorted(possible_connections_now and possible_connections_next):
            ax.hlines(y=y_ticks[all_satellites_map[satellite]], xmin=current_time, xmax=next_time, color='purple', linewidth=3)
            if ue_connectivity_window and ue_connectivity_window[connectivity_index] == 0:
                ax.axvspan(current_time, next_time, color='grey', alpha=1)

        current_time = next_time
        connectivity_index += 1
    
    # Set y-ticks
    ax.set_yticks(y_ticks)
    ax.set_yticklabels(y_labels)
    
    # Labels and title
    plt.xlabel('Time', fontsize=14)
    plt.ylabel('Satellites', fontsize=14)
    plt.title('Satellite Visibility and Connection Over Time', fontsize=16)
    
    # Set x-axis limits
    ax.set_xlim([start_time, end_time])

    ax.tick_params(axis='both', which='major', labelsize=12)

    # Format the x-axis to show only time
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M:%S.%f'))
    
    # Rotate x-axis labels for better readability
    plt.xticks(rotation=45)

    # Grid and layout
    plt.grid(True)
    plt.tight_layout()
    
    plt.show()

def merge_results(results):
    merged_results = {
        "UE": 0,
        "RSRP": [],
        "ToS": [],
        "Effective Time": [],
        "Events": {
            "Programmed Intra-gNB CHO": [],
            "Programmed Inter-gNB CHO": [],
            "Executed Intra-gNB HO": [],
            "Executed Inter-gNB HO": [],
        },
        "Ping-Pong": [],
        "Call Drops": [],
        "Distances": [],
        "Distances HO": [],
        "Positions Call Drop": [],
        "Call Drop Cause": {
            "Call Drop and Event not Triggered yet": [],
            "Call Drop during HO Preparation": [],
            "Call Drop during HO Execution": []
        }
    }

    for ue_result in results:
        merged_results["UE"] += 1
        merged_results["RSRP"].append(ue_result["RSRP"]["EWA"])
        merged_results["ToS"].append(ue_result["ToS"])
        merged_results["Effective Time"].append(ue_result["Effective Time"])
        merged_results["Ping-Pong"].append(ue_result["Ping-Pong"])
        merged_results["Call Drops"].append(ue_result["Call Drops"])
        merged_results["Distances"].append(ue_result["Distances"]["distances"])
        merged_results["Distances HO"].append(ue_result["Distances HO"])
        merged_results["Positions Call Drop"].append(ue_result["Positions Call Drop"])
        for key in merged_results["Events"]:
            merged_results["Events"][key].append(ue_result["Events"][key])
        for key in merged_results["Call Drop Cause"]:
            merged_results["Call Drop Cause"][key].append(ue_result["Call Drop Cause"][key])

    return merged_results

def plot_rsrp(results, start_time, end_time, time_step_seconds, ue_connectivity_window):
    # Set up the plot
    fig, ax = plt.subplots(figsize=(12, 6))

    for scenario, data in results.items():
        line, = ax.plot(data['RSRP']['times'], data['RSRP']['EWA'], label=scenario, linestyle='-')
        mean_value = np.mean(data['RSRP']['EWA'])
        ax.axhline(y=mean_value, color=line.get_color(), linestyle='--', linewidth=1, label=f'{scenario} Mean')

        if ue_connectivity_window:
            current_time = start_time
            connectivity_index = 0
            while connectivity_index < len(ue_connectivity_window):
                if ue_connectivity_window[connectivity_index] == 0:
                    ax.axvspan(current_time, current_time + timedelta(seconds=time_step_seconds), color='grey', alpha=1, zorder=3)
                current_time += timedelta(seconds=time_step_seconds)
                connectivity_index += 1

    # Add labels and title
    ax.set_xlabel('Time', fontsize=14)
    ax.set_ylabel('RSRP (dBm)', fontsize=14)
    ax.set_title('RSRP in the UE of the Connected Satellites Over Time', fontsize=16)
    ax.legend(frameon=True, fontsize=12)
    ax.grid()
    
    ax.set_xlim([start_time, end_time])

    ax.tick_params(axis='both', which='major', labelsize=12)

    # Format the x-axis to show time
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M:%S'))
    plt.xticks(rotation=45)

    # Show the plot
    plt.tight_layout()
    plt.show()

def plot_rsrp_boxplot(results, colors):
    fig, ax = plt.subplots(figsize=(12, 6))

    rsrp_values = []
    scenario_labels = []

    for scenario, data in results.items():
        rsrp_all_ues = []
        for rsrp_ue in data['RSRP']:
            rsrp_all_ues.extend(rsrp_ue)

        rsrp_values.append(rsrp_all_ues)
        scenario_labels.append(scenario)

    boxplot = ax.boxplot(rsrp_values, labels=scenario_labels, patch_artist=True,
               boxprops=dict(facecolor='lightblue', color='blue'),
               medianprops=dict(color='black'),
               whiskerprops=dict(color='black'),
               capprops=dict(color='black'),
               flierprops=dict(marker='o', markerfacecolor='gray', markersize=5, linestyle='none'))
    
    for patch, color in zip(boxplot['boxes'], colors):
        patch.set_facecolor(color)
        patch.set_edgecolor('black')

    for i, rsrp_data in enumerate(rsrp_values):
        mean_val = np.mean(rsrp_data)
        ax.text(i + 1, mean_val + 0.1, f'{mean_val:.1f}', ha='center', va='bottom', fontsize=16, color='black')

    ax.set_xlabel('Strategy', fontsize=18)
    ax.set_ylabel('RSRP (dBm)', fontsize=18)
    ax.set_title('RSRP Distribution by Strategy', fontsize=20)
    ax.grid(axis='y', linestyle='--', alpha=0.7)
    ax.set_ylim([-117.5, -97.5])

    ax.tick_params(axis='both', which='major', labelsize=16)

    plt.tight_layout()
    plt.show()

def plot_localization_uncertainty_comparison(all_results, time, uncertainties, 
                                             colors=None, category_colors=None, call_drop_colors=None):
    scenarios = list(all_results[0].keys())
    n_experiments = len(all_results)

    if colors is None:
        colors = plt.cm.tab10.colors
    if category_colors is None:
        category_colors = plt.cm.Set2.colors
    if call_drop_colors is None:
        call_drop_colors = plt.cm.Pastel1.colors

    x = np.array(uncertainties)

    # Plot of RSRP
    fig, axes = plt.subplots(1, len(scenarios), figsize=(7 * len(scenarios), 6), squeeze=False)
    axes = axes[0]
    for s_idx, scenario in enumerate(scenarios):
        ax = axes[s_idx]
        means, stds = [], []
        for exp in all_results:
            rsrp_all_ues = []
            for rsrp_ue in exp[scenario]["RSRP"]:
                rsrp_all_ues.extend(rsrp_ue)
            means.append(np.mean(rsrp_all_ues) if rsrp_all_ues else 0.0)
            stds.append(np.std(rsrp_all_ues) if rsrp_all_ues else 0.0)

        ax.errorbar(x, means, yerr=stds, fmt="-o", color=colors[s_idx % len(colors)], capsize=5)
        for xi, m, sd in zip(x, means, stds):
            ax.text(xi, m + (sd if sd > 0 else 0.2) + 0.05, f"{m:.1f}", ha="center", va="bottom", fontsize=10)

        ax.set_title(f"RSRP Evolution - {scenario}", fontsize=16)
        ax.set_xlabel("Uncertainty")
        ax.set_ylabel("RSRP (dBm)")
        ax.set_ylim([-117.5, -97.5])
        ax.grid(True, linestyle="--", alpha=0.7)
    plt.tight_layout()
    plt.show()

    # Plot of effective time proportion
    fig, axes = plt.subplots(1, len(scenarios), figsize=(7 * len(scenarios), 6), squeeze=False)
    axes = axes[0]
    for s_idx, scenario in enumerate(scenarios):
        ax = axes[s_idx]
        means = []
        for exp in all_results:
            eff_per_ue = np.array([sum(eff) for eff in exp[scenario]["Effective Time"]], dtype=float)
            eff_prop = eff_per_ue / float(time) if eff_per_ue.size > 0 else np.array([])
            means.append(float(np.mean(eff_prop)) if eff_prop.size > 0 else 0.0)

        ax.plot(x, means, "-o", color=colors[s_idx % len(colors)])
        for xi, m in zip(x, means):
            ax.text(xi, m - 0.04, f"{m:.2f}", ha="center", va="bottom", fontsize=10)

        ax.set_title(f"Effective Time Proportion Evolution - {scenario}", fontsize=16)
        ax.set_xlabel("Uncertainty")
        ax.set_ylabel("Effective Time Proportion")
        ax.set_ylim([0, 1])
        ax.grid(True, linestyle="--", alpha=0.7)
    plt.tight_layout()
    plt.show()

    # Plot of events
    categories = list(all_results[0][scenarios[0]]['Events'].keys())
    fig, axes = plt.subplots(1, len(scenarios), figsize=(7 * len(scenarios), 6), squeeze=False)
    axes = axes[0]
    for s_idx, scenario in enumerate(scenarios):
        ax = axes[s_idx]
        for c_idx, category in enumerate(categories):
            series = [(np.sum(exp[scenario]['Events'][category]) / exp[scenario]['UE'] / float(time))
                       for exp in all_results]
            ax.plot(x, series, "-o", label=category)
            for xi, v in zip(x, series):
                ax.text(xi, v + 0.005, f"{v:.3f}", ha="center", va="bottom", fontsize=9)

        ax.set_title(f"HO Events Evolution - {scenario}", fontsize=16)
        ax.set_xlabel("Uncertainty")
        ax.set_ylabel("Events / UE / s")
        ax.set_ylim([0, 0.4])
        ax.grid(True, linestyle="--", alpha=0.7)
        ax.legend(fontsize=10)
    plt.tight_layout()
    plt.show()

    # Plot of total call drops
    fig, axes = plt.subplots(1, len(scenarios), figsize=(7 * len(scenarios), 6), squeeze=False)
    axes = axes[0]

    for s_idx, scenario in enumerate(scenarios):
        ax = axes[s_idx]
        means_cd = []

        for exp in all_results:
            call_drops = np.array(exp[scenario]["Call Drops"])
            means_cd.append(float(np.sum(call_drops)) if call_drops.size > 0 else 0.0)

        ax.plot(x, means_cd, "-o", color="red", label="Total Call Drops")
        for xi, m in zip(x, means_cd):
            ax.text(xi, m + 0.9, f"{int(m)}", ha="center", va="bottom", fontsize=9, color="black")

        ax.set_title(f"Total Call Drops - {scenario}", fontsize=14)
        ax.set_xlabel("Uncertainty")
        ax.set_ylabel("Total Call Drops")
        ax.set_ylim([0, 1250]) 
        ax.grid(True, linestyle="--", alpha=0.7)
        ax.legend(loc="upper right")

    plt.tight_layout()
    plt.show()

    # Plot of call drop cause evolution
    categories = list(all_results[0][scenarios[0]]['Call Drop Cause'].keys())
    fig, axes = plt.subplots(1, len(scenarios), figsize=(7 * len(scenarios), 6), squeeze=False)
    axes = axes[0]

    for s_idx, scenario in enumerate(scenarios):
        ax = axes[s_idx]

        for c_idx, category in enumerate(categories):
            series = [
                (sum(exp[scenario]['Call Drop Cause'][category]) /
                sum(sum(exp[scenario]['Call Drop Cause'][cat]) for cat in categories)
                if sum(sum(exp[scenario]['Call Drop Cause'][cat]) for cat in categories) > 0 else 0)
                for exp in all_results
            ]
            ax.plot(x, series, "-o", label=category)
            for xi, v in zip(x, series):
                ax.text(xi, v + 0.015, f"{v:.2f}", ha="center", va="bottom", fontsize=9)

        ax.set_title(f"Call Drop Cause Evolution - {scenario}", fontsize=16)
        ax.set_xlabel("Uncertainty")
        ax.set_ylabel("Percentage (per cause)")
        ax.set_ylim([0, 1])
        ax.grid(True, linestyle="--", alpha=0.7)

        ax2 = ax.twinx()
        totals = [
            sum(sum(exp[scenario]['Call Drop Cause'][cat]) for cat in categories)
            for exp in all_results
        ]
        ax2.plot(x, totals, "-o", color="gray", alpha=0.5, linewidth=2.5, label="Total Call Drops")
        for xi, t in zip(x, totals):
            ax2.text(xi, t + 4, f"{int(t)}",
                    ha="center", va="bottom", fontsize=9, color="gray", alpha=0.8)
        ax2.set_ylabel("Total Call Drops")
        ax2.set_ylim([0, 175])

        lines, labels = ax.get_legend_handles_labels()
        lines2, labels2 = ax2.get_legend_handles_labels()
        ax2.legend(lines + lines2, labels + labels2, fontsize=10, loc="upper right")

    plt.tight_layout()
    plt.show()


def plot_tos_cdf(results, colors, linestyles):
    # Set up the plot
    fig, ax = plt.subplots(figsize=(12, 6))

    for i, (scenario, data) in enumerate(results.items()):
        tos_all_ues = []
        for tos_ue in data['ToS']:
            tos_all_ues.extend(tos_ue)

        values, counts = np.unique(tos_all_ues, return_counts=True)
        probabilities = counts / counts.sum()
 
        # Calculate and plot the CDF
        cdf = np.cumsum(probabilities)
        ax.plot(values, cdf, label=scenario, color=colors[i], linestyle=linestyles[i], linewidth=3)
    
    #Add labels and title
    ax.set_xlabel('ToS (s)', fontsize=18)
    ax.set_ylabel('CDF', fontsize=18)
    ax.set_title('CDF of Time of Stay by Strategy', fontsize=20)
    ax.legend(fontsize=16, loc='lower right')
    ax.grid(True, linestyle='--', alpha=0.7)

    ax.set_ylim([0, 1])

    ax.tick_params(axis='both', which='major', labelsize=16)

    # Show the plot
    plt.tight_layout()
    plt.show()


def plot_effective_time(results, time, colors):
    import matplotlib.pyplot as plt
    import numpy as np

    fig, ax = plt.subplots(figsize=(12, 6))

    scenarios = []
    effective_times = []
    error_bars = []

    for scenario, data in results.items():        
        effective_time_all_ues = np.array([sum(eff) for eff in data["Effective Time"]])
        effective_time = effective_time_all_ues / time
        
        mean_time = np.mean(effective_time) if effective_time.size > 0 else 0
        std_dev = np.std(effective_time) if effective_time.size > 0 else 0

        scenarios.append(scenario)
        effective_times.append(mean_time)
        error_bars.append(std_dev)

    bars = ax.bar(scenarios, effective_times, yerr=error_bars, color=colors, alpha=0.7, capsize=5, ecolor='black', edgecolor='black', linewidth=1.5)

    for i, value in enumerate(effective_times):
        if value > 0.85:
            # Dins la barra, centrat
            ax.text(i, value - 0.05, f'{value:.2f}', ha='center', va='center', fontsize=14, color='black')
        else:
            # Per sobre la barra
            ax.text(i, value + 0.05, f'{value:.2f}', ha='center', va='bottom', fontsize=14, color='black')

    ax.set_xlabel('Strategy', fontsize=18)
    ax.set_ylabel('Effective Time Proportion', fontsize=18)
    ax.set_title('Effective Time Proportion by Strategy', fontsize=20)
    ax.grid(axis='y', linestyle='--', alpha=0.7)

    ax.set_ylim(0, 1)

    ax.tick_params(axis='both', which='major', labelsize=16)

    plt.tight_layout()
    plt.show()


def plot_events(results, time):
    scenarios = list(results.keys())
    categories = list(results[scenarios[0]]['Events'].keys())
    values = [[sum(results[scenario]['Events'][category])/results[scenario]['UE']/time for category in categories] for scenario in scenarios]

    bar_width = 0.08
    x = np.arange(len(scenarios)) * (bar_width * (len(categories) + 1))

    fig, ax = plt.subplots(figsize=(12, 6))

    for i, category in enumerate(categories):
        offsets = x + i * bar_width
        bar_vals = [value[i] for value in values]
        bars = ax.bar(offsets, bar_vals, bar_width, label=category, edgecolor='black', linewidth=1)

        for xi, val in zip(offsets, bar_vals):
            ax.text(xi, val + 0.005, f'{val:.3f}', ha='center', va='bottom', fontsize=10)

    ax.set_xlabel('Strategy', fontsize=18)
    ax.set_ylabel('Events/UE/s', fontsize=18)
    ax.set_title('Handover Events Rate by Strategy', fontsize=20)
    ax.set_xticks(x + bar_width * (len(categories) - 1) / 2)
    ax.set_xticklabels(scenarios, fontsize=18)
    ax.legend(fontsize=16, loc="best")
    ax.set_ylim(0, 0.4)
    ax.grid(axis='y', linestyle='--', alpha=0.7)

    ax.tick_params(axis='both', which='major', labelsize=16)

    plt.tight_layout()
    plt.show()

def plot_cho_efficiency(results, colors):
    
    fig, ax = plt.subplots(figsize=(12, 6))

    scenarios = []
    cho_efficiencies = []
    error_bars = []

    for scenario, data in results.items():
        programmed_intra_gNB = np.array(data["Events"]["Programmed Intra-gNB CHO"])
        programmed_inter_gNB = np.array(data["Events"]["Programmed Inter-gNB CHO"])
        executed_intra_gNB = np.array(data["Events"]["Executed Intra-gNB HO"])
        executed_inter_gNB = np.array(data["Events"]["Executed Inter-gNB HO"])
        
        programmed = programmed_intra_gNB + programmed_inter_gNB
        executed = executed_intra_gNB + executed_inter_gNB

        cho_efficiency = []
        for programmed_cho, executed_ho in zip(programmed, executed):
            if programmed_cho > 0:
                cho_efficiency.append(executed_ho/programmed_cho)
            else:
                cho_efficiency.append(0)
        cho_efficiency = np.array(cho_efficiency)

        mean_efficiency = np.mean(cho_efficiency) if cho_efficiency.size > 0 else 0
        std_dev = np.std(cho_efficiency) if cho_efficiency.size > 0 else 0

        scenarios.append(scenario)
        cho_efficiencies.append(mean_efficiency)
        error_bars.append(std_dev)

    ax.bar(scenarios, cho_efficiencies, yerr=error_bars, color=colors, alpha=0.7, capsize=5, ecolor='black', edgecolor='black', linewidth=1.5)

    ax.set_xlabel('Strategy', fontsize=18)
    ax.set_ylabel('CHO Efficiency', fontsize=18)
    ax.set_title('CHO Efficiency by Strategy', fontsize=20)
    ax.grid(axis='y', linestyle='--', alpha=0.7)

    ax.set_ylim(0, 1)

    ax.tick_params(axis='both', which='major', labelsize=16)

    plt.tight_layout()
    plt.show()

def plot_call_drop_probability(results, colors):
    import matplotlib.pyplot as plt
    import numpy as np

    fig, ax = plt.subplots(figsize=(12, 6))

    scenarios = []
    call_drop_probabilities = []
    error_bars_call_drop = []
    total_call_drops = []
    ping_pong_probabilities = []

    for scenario, data in results.items():
        call_drops = np.array(data["Call Drops"])
        intra_gNB = np.array(data["Events"]["Executed Intra-gNB HO"])
        inter_gNB = np.array(data["Events"]["Executed Inter-gNB HO"])
        total_event = call_drops + intra_gNB + inter_gNB

        call_drop_probability = [
            cd / ev if ev > 0 else 0 for cd, ev in zip(call_drops, total_event)
        ]
        call_drop_probability = np.array(call_drop_probability)

        mean_cd = np.mean(call_drop_probability) if call_drop_probability.size > 0 else 0
        std_cd = np.std(call_drop_probability) if call_drop_probability.size > 0 else 0
        total_cd = int(np.sum(call_drops))

        scenarios.append(scenario)
        call_drop_probabilities.append(mean_cd)
        error_bars_call_drop.append(std_cd)
        total_call_drops.append(total_cd)

        ping_pongs = np.array(data["Ping-Pong"])
        total_handover = intra_gNB + inter_gNB

        ping_pong_probability = [
            pp / ho if ho > 0 else 0 for pp, ho in zip(ping_pongs, total_handover)
        ]
        ping_pong_probability = np.array(ping_pong_probability)

        mean_pp = np.mean(ping_pong_probability) if ping_pong_probability.size > 0 else 0
        ping_pong_probabilities.append(mean_pp)

    all_pp_low = all(prob < 0.1 for prob in ping_pong_probabilities)

    if all_pp_low:
        bars = ax.bar(
            scenarios, call_drop_probabilities,
            yerr=error_bars_call_drop, color=colors, alpha=0.7,
            capsize=5, ecolor='black', edgecolor='black', linewidth=1.5
        )

        for i, val in enumerate(call_drop_probabilities):
            err = error_bars_call_drop[i]
            if val + err > 0.85:
                ax.text(i, val / 2, f'{val:.2f}', ha='center', va='center',
                        fontsize=12, color='black')
            else:
                ax.text(i, val + err + 0.02, f'{val:.2f}', ha='center', va='bottom',
                        fontsize=12, color='black')

        legend_labels = [f"total call drops = {total_call_drops[i]}" for i in range(len(scenarios))]
        ax.legend(bars, legend_labels, fontsize=14, loc="upper right")

        label = 'Call Drop Probability'
        title = 'Call Drop Probability by Strategy'

    else:
        bar_width = 0.35
        x = np.arange(len(scenarios))

        bars_cd = ax.bar(
            x - bar_width / 2, call_drop_probabilities, bar_width,
            yerr=error_bars_call_drop, label='Call Drop',
            color=colors, capsize=5, edgecolor='black', linewidth=1.5
        )
        bars_pp = ax.bar(
            x + bar_width / 2, ping_pong_probabilities, bar_width,
            label='Ping-Pong', color='gray', edgecolor='black', linewidth=1.5
        )

        ax.set_xticks(x)
        ax.set_xticklabels(scenarios, fontsize=18)

        for i, val in enumerate(call_drop_probabilities):
            err = error_bars_call_drop[i]
            xpos = x[i] - bar_width / 2
            if val + err > 0.85:
                ax.text(xpos, val / 2, f'{val:.2f}', ha='center', va='center',
                        fontsize=12, color='black')
            else:
                ax.text(xpos, val + err + 0.02, f'{val:.2f}', ha='center', va='bottom',
                        fontsize=12, color='black')

        for i, val in enumerate(ping_pong_probabilities):
            xpos = x[i] + bar_width / 2
            if val > 0.85:
                ax.text(xpos, val / 2, f'{val:.2f}', ha='center', va='center',
                        fontsize=12, color='black')
            else:
                ax.text(xpos, val + 0.04, f'{val:.2f}', ha='center', va='bottom',
                        fontsize=12, color='black')

        legend_labels = [f"total call drops = {total_call_drops[i]}" for i in range(len(scenarios))]
        ax.legend(bars_cd, legend_labels, fontsize=14, loc="upper right")

        label = 'Probability'
        title = 'Call Drop and Ping-Pong Probability by Strategy'

    ax.set_xlabel('Strategy', fontsize=18)
    ax.set_ylabel(label, fontsize=18)
    ax.set_title(title, fontsize=20)
    ax.grid(axis='y', linestyle='--', alpha=0.7)
    ax.set_ylim(0, 1)
    ax.tick_params(axis='both', which='major', labelsize=16)

    plt.tight_layout()
    plt.show()


def plot_call_drop_distribution_by_zone(results, colors, linestyles):
    fig, ax = plt.subplots(figsize=(12, 6))

    for i, (scenario, data) in enumerate(results.items()):
        positions = data.get("Positions Call Drop", [])
        all_positions = [p for user in positions for p in user]

        long_diff = []
        for pos in all_positions:
            beam_lat = pos[2]
            long = pos[0]*111*np.cos(np.deg2rad(beam_lat))
            long_diff.append(long)

        total_cd = len(long_diff)
        print(f"[{scenario}] Total Call Drops in Zone: {total_cd}")

        if total_cd == 0:
            ax.axhline(y=0, color=colors[i], linestyle=linestyles[i], linewidth=3, label=scenario)
        else:
            hist_vals, bin_edges = np.histogram(long_diff, bins=10, range=(-950, 950), density=False)
            hist_vals = hist_vals / total_cd
            hist_vals = np.append(hist_vals, hist_vals[-1])
            ax.step(bin_edges, hist_vals, where='post', label=f"{scenario} (N={total_cd})", color=colors[i], linestyle=linestyles[i], linewidth=3)

    ax.set_xlabel('Distance from center in the W-E direction (km)', fontsize=18)
    ax.set_ylabel('Call Drop Probability', fontsize=18)
    ax.set_title('Call Drop Probability by Beam Zone and Strategy', fontsize=20)
    ax.set_ylim(0, 1)
    ax.set_xlim(-950, 950)
    ax.legend(fontsize=16, loc="upper right")
    ax.grid(True, linestyle='--', alpha=0.7)
    ax.tick_params(axis='both', which='major', labelsize=16)

    plt.tight_layout()
    plt.show()

def plot_call_drop_reason(results):
    import matplotlib.pyplot as plt
    import numpy as np

    scenarios = list(results.keys())
    categories = list(results[scenarios[0]]['Call Drop Cause'].keys())

    percentages = []
    totals = []  
    for scenario in scenarios:
        total = sum(sum(results[scenario]['Call Drop Cause'][category]) for category in categories)
        totals.append(total) 
        if total == 0:
            percentages.append([0 for _ in categories])
        else:
            percentages.append([(sum(results[scenario]['Call Drop Cause'][category]) / total) for category in categories])
        print(f"[{scenario}] Total Call Drops (Causes): {total}")  

    bar_width = 0.08
    x = np.arange(len(scenarios)) * (bar_width * (len(categories) + 1))

    fig, ax = plt.subplots(figsize=(12, 6))

    for i, category in enumerate(categories):
        offsets = x + i * bar_width
        bar_vals = [value[i] for value in percentages]
        bars = ax.bar(offsets, bar_vals, bar_width, label=category, edgecolor='black', linewidth=1)

        for xi, val in zip(offsets, bar_vals):
            if val > 0.85:
                ax.text(xi, val / 2, f'{val:.2f}', ha='center', va='center', fontsize=10, color='white')
            else:
                ax.text(xi, val + 0.015, f'{val:.2f}', ha='center', va='bottom', fontsize=10, color='black')

    # Mostrar els totals sota els grups
    for i, total in enumerate(totals):
        ax.text(x[i] + bar_width * (len(categories) - 1) / 2, -0.05, f'total call drops = {total}', ha='center', va='top', fontsize=10)

    ax.set_xlabel('Strategy', fontsize=18)
    ax.set_ylabel('Call Drop Cause Distribution', fontsize=18)
    ax.set_title('Call Drop Cause Distribution by Strategy', fontsize=20)
    ax.set_xticks(x + bar_width * (len(categories) - 1) / 2)
    ax.set_xticklabels(scenarios, fontsize=18)
    ax.set_ylim(-0.1, 1)
    ax.legend(fontsize=16)
    ax.grid(axis='y', linestyle='--', alpha=0.7)
    ax.tick_params(axis='both', which='major', labelsize=16)

    plt.tight_layout()
    plt.show()



def plot_distance(results, start_time, end_time, time_step_seconds, ue_connectivity_window):
    # Set up the plot
    fig, ax = plt.subplots(figsize=(12, 6))

    for scenario, data in results.items():

        ax.plot(data['Distances']['times'], data['Distances']['distances'], label=scenario, linestyle='-')

        if ue_connectivity_window:
            current_time = start_time
            connectivity_index = 0
            while connectivity_index < len(ue_connectivity_window):
                if ue_connectivity_window[connectivity_index] == 0:
                    ax.axvspan(current_time, current_time + timedelta(seconds=time_step_seconds), color='grey', alpha=1, zorder=3)
                current_time += timedelta(seconds=time_step_seconds)
                connectivity_index += 1

    # Add labels and title
    ax.set_xlabel('Time', fontsize=18)
    ax.set_ylabel('Distance (m)', fontsize=18)
    ax.set_title('Distance from the UE to the Cell center of the Connected Satellites Over Time', fontsize=20)
    ax.legend(frameon=True, fontsize=16)
    ax.grid()
    
    ax.set_xlim([start_time, end_time])

    ax.tick_params(axis='both', which='major', labelsize=16)

    # Format the x-axis to show time
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M:%S'))

    # Show the plot
    plt.tight_layout()
    plt.show()

def plot_distance_cdf(results, colors, linestyles):
    # Set up the plot
    fig, ax = plt.subplots(figsize=(12, 6))

    for i, (scenario, data) in enumerate(results.items()):
        distance_all_ues = []
        for distance_ue in data['Distances']:
            distance_all_ues.extend(distance_ue)

        values, counts = np.unique(distance_all_ues, return_counts=True)
        probabilities = counts / counts.sum()
 
        # Calculate and plot the CDF
        cdf = np.cumsum(probabilities)
        ax.plot(values, cdf, label=scenario, color=colors[i], linestyle=linestyles[i])
    
    #Add labels and title
    ax.set_xlabel('Distance UE-Cell Center (m)', fontsize=18)
    ax.set_ylabel('CDF', fontsize=18)
    ax.set_title('CDF of Distances UE-Cell Center by Strategy', fontsize=20)
    ax.legend(fontsize=16, loc="lower right")
    ax.grid()

    ax.set_ylim([0, 1])

    ax.tick_params(axis='both', which='major', labelsize=16)

    # Show the plot
    plt.tight_layout()
    plt.show()

def plot_distance_ho_cdf(results, colors, linestyles):
    # Set up the plot
    fig, ax = plt.subplots(figsize=(12, 6))

    for i, (scenario, data) in enumerate(results.items()):
        distance_all_ues = []
        for distance_ue in data['Distances HO']:
            distance_all_ues.extend(distance_ue)
        
        values, counts = np.unique(distance_all_ues, return_counts=True)
        probabilities = counts / counts.sum()
 
        # Calculate and plot the CDF
        cdf = np.cumsum(probabilities)
        ax.plot(values, cdf, label=scenario, color=colors[i], linestyle=linestyles[i])
    
    #Add labels and title
    ax.set_xlabel('Distance UE-Cell Center at HO (m)', fontsize=18)
    ax.set_ylabel('CDF', fontsize=18)
    ax.set_title('CDF of Distances UE-Cell Center at Handover by Strategy', fontsize=20)
    ax.legend(fontsize=16, loc="lower right")
    ax.grid()

    ax.set_ylim([0, 1])

    ax.tick_params(axis='both', which='major', labelsize=16)

    # Show the plot
    plt.tight_layout()
    plt.show()

def compute_area(min_lat, max_lat, min_long, max_long):
    lat_diff = max_lat - min_lat
    lon_diff = max_long - min_long

    lat = (max_lat + min_lat) / 2

    lon_distance = lon_diff * 111 * np.cos(np.deg2rad(lat))

    lat_distance = lat_diff * 111

    area = lat_distance * lon_distance

    return area

def plot_ho_rate_per_ue_density(density_results, time, area):
    if area <= 0 or time <= 0:
        raise ValueError("Time and area must be greater than 0.")

    fig, ax = plt.subplots(figsize=(12, 6))

    colors = plt.cm.get_cmap("tab10")

    for idx, (strategy, results) in enumerate(density_results.items()):
        intra_rates = []
        inter_rates = []
        densities = []

        for data in results:
            density = data["UE"] / area
            intra = sum(data["Events"]["Executed Intra-gNB HO"]) / time
            inter = sum(data["Events"]["Executed Inter-gNB HO"]) / time

            densities.append(density)
            intra_rates.append(intra)
            inter_rates.append(inter)
        
        sorted_data = sorted(zip(densities, intra_rates, inter_rates), key=lambda x: x[0])
        densities, intra_rates, inter_rates = zip(*sorted_data)

        color = colors(idx % 10)

        ax.plot(densities, intra_rates, label=f"{strategy} - Intra-gNB", color=color, linestyle='-', marker='o')
        ax.plot(densities, inter_rates, label=f"{strategy} - Inter-gNB", color=color, linestyle='--', marker='s')

    ax.set_xlabel('User Density (UE/km²)', fontsize=14)
    ax.set_ylabel('Handover Rate (HO/s)', fontsize=14)
    ax.set_title('Handover Rate vs UE Density by Strategy', fontsize=16)

    ax.grid()
    ax.legend(fontsize=12)
    ax.tick_params(axis='both', which='major', labelsize=12)
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.show()

def plot_ho_rate_vs_time(events, window_size_seconds, start_time, end_time, time_step_seconds):
    window_size = timedelta(seconds=window_size_seconds)
    step = timedelta(seconds=time_step_seconds)

    time_points = []
    intra_counts = []
    inter_counts = []

    current_time = start_time
    while current_time <= end_time:
        window_start = current_time
        window_end = current_time + window_size

        if window_end > end_time:
            break

        intra = 0
        inter = 0

        for ue_events in events.values():
            for event in ue_events:
                if window_start <= event["time"] < window_end:
                    if event["event"] == "Intra-gNB":
                        intra += 1
                    elif event["event"] == "Inter-gNB":
                        inter += 1

        time_points.append(current_time)
        intra_counts.append(intra / window_size_seconds)
        inter_counts.append(inter / window_size_seconds)

        current_time += step

    plt.figure(figsize=(12, 6))
    plt.step(time_points, intra_counts, label="Intra-gNB HO/s", color='blue', where='post')
    plt.step(time_points, inter_counts, label="Inter-gNB HO/s", color='green', where='post')
    plt.xlabel("Time", fontsize=14)
    plt.ylabel("Handover Rate (HO/s)", fontsize=14)
    plt.title(f"Handover Rate (Window = {window_size_seconds}s)", fontsize=16)
    plt.legend(fontsize=12)
    plt.grid(True)
    plt.tight_layout()
    plt.xticks(rotation=45)

    plt.gca().xaxis.set_major_formatter(plt.matplotlib.dates.DateFormatter('%H:%M:%S'))
    plt.show()


def satellites_coverage(data_csv, start_time, end_time, time_step_seconds):
    try:
        # Load the visibility data from CSV (handle semicolon delimiter)
        visibility_df = pd.read_csv(data_csv, delimiter=';')
        visibility_df['Time'] = pd.to_datetime(visibility_df['Time'])
    except pd.errors.EmptyDataError:
        print(f"Error: The file '{data_csv}' is empty or not properly formatted.")
        return
    except Exception as e:
        print(f"An error occurred while reading '{data_csv}': {e}")
        return
    
    # Check and handle the expected columns
    visibility_column = 'Visible_satellites'
    position_column = 'Satellite_positions'  # For parsing positions if needed
    if visibility_column not in visibility_df.columns:
        print(f"Error: Expected column '{visibility_column}' not found in the CSV.")
        return

    current_time = start_time
    satellites_number = []
    time = []
    while current_time < end_time:
        current_time += timedelta(seconds=time_step_seconds)

        potential_satellites = set()
        for satellites in visibility_df[visibility_df['Time'] == current_time]['Possible_to_connect_satellites']:
            if satellites != '{}':
                satellites_list = satellites.strip('{}').split(', ')
                potential_satellites.update(satellites_list)
        potential_satellites = sorted(potential_satellites)
        
        satellites_number.append(len(potential_satellites))
        time.append(current_time)

    plt.figure(figsize=(12, 8))
    plt.step(time, satellites_number, where='post')

    plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%H:%M:%S'))
    plt.xticks(rotation=45)

    plt.xlabel('Time', fontsize=14)
    plt.ylabel('Number of Satellites', fontsize=14)
    plt.title('Evolution of Number of Satellites Giving Coverage to UE', fontsize=16)
    plt.grid(True)
    plt.tick_params(axis='both', which='major', labelsize=12)
    plt.tight_layout()

    # Set x-axis limits
    plt.xlim([start_time, end_time])

    plt.show()

    return satellites_number

def plot_satellites_probabilities(satellites):
    combined_results = np.concatenate(list(satellites.values()))
    
    # Set up the plot
    fig, ax1 = plt.subplots(figsize=(12, 6))

    values, counts = np.unique(combined_results, return_counts=True)
    probabilities = counts / counts.sum()
 
    # Calculate and plot the CDF
    cdf = np.cumsum(probabilities)
    
    ax1.bar(values, probabilities, width=0.5, alpha=0.7, color='blue', label="PDF")
    ax1.set_xlabel("Satellites Giving Coverage", fontsize=14)
    ax1.set_ylabel("PDF (Percentage of Time)", fontsize=14)
    #ax1.tick_params(axis='y', labelcolor='blue')
    ax1.grid(axis='y', linestyle='--', alpha=0.6)
    ax1.tick_params(axis='both', which='major', labelsize=12)

    ax2 = ax1.twinx()
    ax2.plot(values, cdf, color='orange', label="CDF", marker='o', linestyle='-', linewidth=2)
    ax2.set_ylabel("CDF (Accumulated Percentage)", fontsize=14)
    #ax2.tick_params(axis='y', labelcolor='orange')
    ax2.tick_params(axis='both', which='major', labelsize=12)

    plt.title("PDF and CDF of Satellites Giving Coverage to the UE", fontsize=16)
    fig.tight_layout()
    fig.legend(loc="upper left", bbox_to_anchor=(0.1, 0.9), fontsize=12)

    plt.show()

def plot_waiting_time(results):
    waiting_times = []

    for data in results.values():
        waiting_time = 0
        for satellites_number in data:
            if satellites_number == 0:
                waiting_time += 1
            elif satellites_number != 0 and waiting_time != 0:
                waiting_times.append(waiting_time)
                waiting_time = 0

    bins = [0, 60, 300, 900, np.inf]
    
    labels = ["Less than a minute", "1 to 5 minutes", "5 to 15 minutes", "More than 15 minutes"]
    
    counts, _ = np.histogram(waiting_times, bins=bins)
    
    probabilities = counts / sum(counts)
    
    plt.figure(figsize=(12, 6))
    plt.bar(labels, probabilities, color="blue", edgecolor="black")
    
    plt.xlabel("Waiting Time Categories", fontsize=14)
    plt.ylabel("Probability", fontsize=14)
    plt.title("Histogram of Waiting Times", fontsize=16)
    plt.grid(axis='y', linestyle='--', alpha=0.7)
    plt.tick_params(axis='both', which='major', labelsize=12)
    plt.tight_layout()
    
    plt.show()

def plot_service_time(results):
    service_times = []

    for data in results.values():
        service_time = 0
        for satellites_number in data:
            if satellites_number != 0:
                service_time += 1
            elif satellites_number == 0 and service_time != 0:
                service_times.append(service_time)
                service_time = 0

    bins = [0, 60, 180, 300, np.inf]
    
    labels = ["Less than a minute", "1 to 3 minutes", "3 to 5 minutes", "More than 5 minutes"]
    
    counts, _ = np.histogram(service_times, bins=bins)
    
    probabilities = counts / sum(counts)
    
    plt.figure(figsize=(12, 6))
    plt.bar(labels, probabilities, color="blue", edgecolor="black")
    
    plt.xlabel("Service Time Categories", fontsize=14)
    plt.ylabel("Probability", fontsize=14)
    plt.title("Histogram of Serving Times", fontsize=16)
    plt.grid(axis='y', linestyle='--', alpha=0.7)
    plt.tick_params(axis='both', which='major', labelsize=12)
    plt.tight_layout()
    
    plt.show()

def plot_orbits_used(data_csv, neighbors_csv, start_time, end_time, time_step_seconds):
    try:
        # Load the visibility data from CSV (handle semicolon delimiter)
        visibility_df = pd.read_csv(data_csv, delimiter=';')
        visibility_df['Time'] = pd.to_datetime(visibility_df['Time'])
    except pd.errors.EmptyDataError:
        print(f"Error: The file '{data_csv}' is empty or not properly formatted.")
        return
    except Exception as e:
        print(f"An error occurred while reading '{data_csv}': {e}")
        return
    
    # Check and handle the expected columns
    visibility_column = 'Visible_satellites'
    position_column = 'Satellite_positions'  # For parsing positions if needed
    if visibility_column not in visibility_df.columns:
        print(f"Error: Expected column '{visibility_column}' not found in the CSV.")
        return
    
    try:
        neighbors_df = pd.read_csv(neighbors_csv, delimiter=';')
    except pd.errors.EmptyDataError:
        print(f"Error: The file '{neighbors_csv}' is empty or not properly formatted.")
        return
    except Exception as e:
        print(f"An error occurred while reading '{neighbors_csv}': {e}")
        return

    satellite_column = 'Satellite'
    orbit_column = 'Orbit'
    if satellite_column not in neighbors_df.columns:
        print(f"Error: Expected column '{satellite_column}' not found in the CSV.")
        return

    current_time = start_time
    satellites_orbits = []
    while current_time < end_time:
        current_time += timedelta(seconds=time_step_seconds)
        orbits = []

        potential_satellites = set()
        for satellites in visibility_df[visibility_df['Time'] == current_time]['Possible_to_connect_satellites']:
            if satellites != '{}':
                satellites_list = satellites.strip('{}').split(', ')
                for satellite in satellites_list:
                    satellite_orbit = neighbors_df.loc[neighbors_df['Satellite'] == satellite, 'Orbit'].values[0]
                    orbits.append(int(satellite_orbit))
                potential_satellites.update(satellites_list)
        potential_satellites = sorted(potential_satellites)
        
        satellites_orbits.append(orbits)

    # Define y-ticks and labels
    orbit_values = neighbors_df['Orbit'].unique().tolist()
    y_ticks = orbit_values
    y_labels = orbit_values
    
    # Prepare the figure and axis
    fig, ax = plt.subplots(figsize=(12, 8))

    i = 0
    current_time = start_time
    while i < len(satellites_orbits) - 1:
        current_time += timedelta(seconds=time_step_seconds)

        visible_orbits = satellites_orbits[i]
        visible_orbits2 = satellites_orbits[i + 1]
        
        for orbit in orbit_values:
            if (orbit in visible_orbits) and (orbit in visible_orbits2):
                ax.hlines(y=y_ticks[orbit_values.index(orbit)], xmin=current_time, xmax=current_time + timedelta(seconds=time_step_seconds), color='purple', linewidth=3)
        i += 1
    
    # Set y-ticks
    ax.set_yticks(y_ticks)
    ax.set_yticklabels(y_labels)
    
    # Labels and title
    plt.xlabel('Time', fontsize=14)
    plt.ylabel('Orbits', fontsize=14)
    plt.title('Orbits Used Over Time', fontsize=16)
    
    # Set x-axis limits
    ax.set_xlim([start_time, end_time])

    # Format the x-axis to show only time
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M:%S'))

    ax.tick_params(axis='both', which='major', labelsize=12)
    
    # Rotate x-axis labels for better readability
    plt.xticks(rotation=45)
    
    # Grid and layout
    plt.grid(True)
    plt.tight_layout()
    
    plt.show()
