#!/usr/bin/python3

from pkg_resources import resource_filename
import sys
import pandas as pd
import datetime
from collections import defaultdict
from satsim.simulations import *
from satsim.ue import connectivity_window, plot_ue_connectivity_window


# Define the path to the data directory
data_path = resource_filename("satsim", "output_data/")
data_path2 = resource_filename("satsim", "data/")

def original_plots():
    strategies_names = ["A3 HO\n(3 neighbors)", "A3 HO\n(4 neighbors)", "D2 + A4 CHO\n(3 neighbors)", "D2 + A4 CHO\n(4 neighbors)"]
    colors = ['salmon', 'salmon', 'lightskyblue', 'lightskyblue']
    linestyles = ['-', '--', '-', '--']
    strategy = ["Signal", "Signal", "Distance", "Distance"]
    handover_type = ["HO", "HO", "CHO", "CHO"]
    a3_offset = [0, 0, 0, 0]
    a3_hys = [0.5, 0.5, 0.5, 0.5]
    a4_threshold = [-115, -115, -115, -115]
    a4_hys = [2.5, 2.5, 2.5, 2.5]
    d2_threshold1 = [[32500, 700000], [32500, 700000], [32500, 700000], [32500, 700000]]
    d2_threshold2 = [[850000, 250000], [850000, 250000], [850000, 250000], [850000, 250000]]
    d2_hys = [0, 0, 0, 0]
    ttt = [0.1, 0.1, 0.1, 0.1]
    num_neighbors = [3, 4, 3, 4]
    k = 4
    ue_connectivity_window = []

    return strategies_names, strategy, handover_type, a3_offset, a3_hys, a4_threshold, a4_hys, d2_threshold1, d2_threshold2, d2_hys, ttt, num_neighbors, k, ue_connectivity_window

def ttt_fine_tunning_plots():
    strategies_names = ["A3 HO\n(TTT = 0 s)", "A3 HO\n(TTT = 0.1 s)", "A3 HO\n(TTT = 0.2 s)", "A3 HO\n(TTT = 0.3 s)", "A3 HO\n(TTT = 0.5 s)"]
    colors = ['salmon', 'lightskyblue', 'palegreen', 'gold', 'plum']
    linestyles = ['-', '-', '-', '-', '-']
    strategy = ["Signal", "Signal", "Signal", "Signal", "Signal"]
    handover_type = ["HO", "HO", "HO", "HO", "HO"]
    a3_offset = [0, 0, 0, 0, 0]
    a3_hys = [0.5, 0.5, 0.5, 0.5, 0.5]
    a4_threshold = [-115, -115, -115, -115, -115]
    a4_hys = [2.5, 2.5, 2.5, 2.5, 2.5]
    d2_threshold1 = [[32500, 700000], [32500, 700000], [32500, 700000], [32500, 700000], [32500, 700000]]
    d2_threshold2 = [[850000, 250000], [850000, 250000], [850000, 250000], [850000, 250000], [850000, 250000]]
    d2_hys = [0, 0, 0, 0, 0]
    ttt = [0, 0.1, 0.2, 0.3, 0.5]
    num_neighbors = [3, 3, 3, 3, 3]
    k = 4
    ue_connectivity_window = []
    
    return strategies_names, strategy, handover_type, a3_offset, a3_hys, a4_threshold, a4_hys, d2_threshold1, d2_threshold2, d2_hys, ttt, num_neighbors, k, ue_connectivity_window

def exact_location_plots():
    strategies_names = ["A3 HO\n(exact location, ttt = 0.1, equator)", "D2 + A4 CHO\n(exact location, ttt = 0.1, equator)"]
    colors = ['salmon', 'lightskyblue']
    linestyles = ['-', '-']
    strategy = ["Signal", "Distance"]
    handover_type = ["HO", "CHO"]
    a3_offset = [0, 0]
    a3_hys = [0.5, 0.5]
    a4_threshold = [-115, -115]
    a4_hys = [2.5, 2.5]
    d2_threshold1 = [[32500, 700000], [32500, 700000]]
    d2_threshold2 = [[850000, 950000], [850000, 950000]]
    d2_hys = [0, 0]
    ttt = [0.1, 0.1]
    num_neighbors = [3, 3]
    k = 4
    ue_connectivity_window = []

    return strategies_names, colors, linestyles, strategy, handover_type, a3_offset, a3_hys, a4_threshold, a4_hys, d2_threshold1, d2_threshold2, d2_hys, ttt, num_neighbors, k, ue_connectivity_window

def distance_exact_location():
    strategies_names = ["D2 + A4 CHO\n(exact location, ttt = 0.1, equator)"]
    colors = ['lightskyblue']
    linestyles = ['-']
    strategy = ["Distance"]
    handover_type = ["CHO"]
    a3_offset = [0]
    a3_hys = [0.5]
    a4_threshold = [-115]
    a4_hys = [2.5]
    d2_threshold1 = [[32500, 700000]]
    d2_threshold2 = [[850000, 950000]]
    d2_hys = [0]
    ttt = [0.1]
    num_neighbors = [3]
    k = 4
    ue_connectivity_window = []

    return strategies_names, colors, linestyles, strategy, handover_type, a3_offset, a3_hys, a4_threshold, a4_hys, d2_threshold1, d2_threshold2, d2_hys, ttt, num_neighbors, k, ue_connectivity_window
  

def three_neighbors_original_plots():
    strategies_names = ["A3 HO\n(3 neighbors)", "D2 + A4 CHO\n(3 neighbors)"]
    colors = ['salmon', 'lightskyblue']
    linestyles = ['-', '-']
    strategy = ["Signal", "Distance"]
    handover_type = ["HO", "CHO"]
    a3_offset = [0, 0]
    a3_hys = [0.5, 0.5]
    a4_threshold = [-115, -115]
    a4_hys = [2.5, 2.5]
    d2_threshold1 = [[32500, 700000], [32500, 700000]]
    d2_threshold2 = [[850000, 950000], [850000, 950000]]
    d2_hys = [0, 0]
    ttt = [0.1, 0.1]
    num_neighbors = [3, 3]
    k = 4
    ue_connectivity_window = []

    return strategies_names, colors, linestyles, strategy, handover_type, a3_offset, a3_hys, a4_threshold, a4_hys, d2_threshold1, d2_threshold2, d2_hys, ttt, num_neighbors, k, ue_connectivity_window


if __name__ == "__main1__":
    datetime_str_1 = '03/26/25 01:00:00.000'
    datetime_str_2 = '03/26/25 02:00:00.000'
    start_time = datetime.datetime.strptime(datetime_str_1, '%m/%d/%y %H:%M:%S.%f')
    end_time = datetime.datetime.strptime(datetime_str_2, '%m/%d/%y %H:%M:%S.%f')
    time_difference = end_time - start_time
    time = time_difference.total_seconds()

    time_step_seconds = 0.1

    min_time_connected = 30
    max_time_connected = 300
    avg_time_connected = 180
    min_time_sleeping = 30
    max_time_sleeping = 240
    avg_time_sleeping = 120

    simulations = 100

    localization_uncertainties = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]

    all_results_uncertainties = []

    # Mètriques que volem visualitzar
    metrics = ["RSRP", "Effective Time", "Events", "Call Drop Probability", "Call Drop Reason"]

    # Etiquetes dels eixos Y
    ylabels = {
        "RSRP": "RSRP [dBm]",
        "Effective Time": "Effective Time Proportion",
        "Events": "Events/UE/s",
        "Call Drop Probability": "Call Drop Probability",
        "Call Drop Reason": "Call Drop Reason"
    }


    for uncertainty in localization_uncertainties:

        # Write the desired function to plot
        strategies_names, colors, linestyles, strategy, handover_type, a3_offset, a3_hys, a4_threshold, a4_hys, d2_threshold1, d2_threshold2, d2_hys, ttt, num_neighbors, k, ue_connectivity_window = exact_location_plots()
    
        raw_results = {}
        all_strategies_results = {}
        all_strategies_events = {}

        for index, strategy_name in enumerate(strategies_names):
            all_results = []
            all_events = defaultdict(list)

            for i in range(simulations):
                results, events_data, satellites_data = handover_strategies(f"{data_path}ue_north{i}.txt", f"{data_path2}neighbors_oneweb.txt", start_time, end_time, time_step_seconds, ue_connectivity_window, "OneWeb", strategy[index], handover_type[index], a3_offset[index], a3_hys[index], a4_threshold[index], a4_hys[index], d2_threshold1[index], d2_threshold2[index], d2_hys[index], k, ttt[index], num_neighbors[index], uncertainty)
                #plot_connectivity(events_data, satellites_data, start_time, end_time, time_step_seconds, ue_connectivity_window)
                
                all_results.append(results)
                for event in events_data:
                    all_events[i].append(event)

                print("Simulation ", i, "done.")
                
            raw_results[strategy_name] = all_results
            all_strategies_results[strategy_name] = merge_results(all_results)
            all_strategies_events[strategy_name] = all_events
            print("Strategy ", index, "done.")

        all_results_uncertainties.append(all_strategies_results)
        with open(f"{data_path}results_north.json", "w", encoding="utf-8") as f:
            json.dump(all_results_uncertainties, f, ensure_ascii=False, indent=2)

        print("Uncertainty ", uncertainty, "done.")
    
    plot_localization_uncertainty_comparison(all_results_uncertainties, 3600, localization_uncertainties, colors, colors)

if __name__ == "__main__":
    datetime_str_1 = '03/26/25 01:00:00.000'
    datetime_str_2 = '03/26/25 02:00:00.000'
    start_time = datetime.datetime.strptime(datetime_str_1, '%m/%d/%y %H:%M:%S.%f')
    end_time = datetime.datetime.strptime(datetime_str_2, '%m/%d/%y %H:%M:%S.%f')
    time_difference = end_time - start_time
    time = time_difference.total_seconds()

    time_step_seconds = 0.1

    min_time_connected = 30
    max_time_connected = 300
    avg_time_connected = 180
    min_time_sleeping = 30
    max_time_sleeping = 240
    avg_time_sleeping = 120

    simulations = 100
    uncertainty = 0.0
    #Write the desired function to plot
    strategies_names, colors, linestyles, strategy, handover_type, a3_offset, a3_hys, a4_threshold, a4_hys, d2_threshold1, d2_threshold2, d2_hys, ttt, num_neighbors, k, ue_connectivity_window = exact_location_plots()
 
    raw_results = {}
    all_strategies_results = {}
    all_strategies_events = {}

    for index, strategy_name in enumerate(strategies_names):
        all_results = []
        all_events = defaultdict(list)

        for i in range(simulations):
            results, events_data, satellites_data = handover_strategies(f"{data_path}ue{i}.txt", f"{data_path2}neighbors_oneweb.txt", start_time, end_time, time_step_seconds, ue_connectivity_window, "OneWeb", strategy[index], handover_type[index], a3_offset[index], a3_hys[index], a4_threshold[index], a4_hys[index], d2_threshold1[index], d2_threshold2[index], d2_hys[index], k, ttt[index], num_neighbors[index], uncertainty)
            #plot_connectivity(events_data, satellites_data, start_time, end_time, time_step_seconds, ue_connectivity_window, i, data_path)
            all_results.append(results)
            for event in events_data:
                all_events[i].append(event)

            print("Simulation ", i, "done.")
            
        raw_results[strategy_name] = all_results
        all_strategies_results[strategy_name] = merge_results(all_results)
        all_strategies_events[strategy_name] = all_events
        print("Strategy ", index, "done.")
    
    plot_rsrp_boxplot(all_strategies_results, colors)
    plot_tos_cdf(all_strategies_results, colors, linestyles)
    plot_effective_time(all_strategies_results, time, colors)
    plot_events(all_strategies_results, time)
    plot_call_drop_probability(all_strategies_results, colors)
    plot_call_drop_distribution_by_zone(all_strategies_results, colors, linestyles)
    plot_call_drop_reason(all_strategies_results)

    plot_cho_efficiency(all_strategies_results, colors)
    plot_distance_cdf(all_strategies_results, colors, linestyles)
    plot_distance_ho_cdf(all_strategies_results, colors, linestyles)
    
    '''
    min_lat, max_lat = -0.4, 0.4
    min_long, max_long = -8.5, 8.5
    #min_lat, max_lat = 64, 64.8
    #min_long, max_long = -22, -13
    area = compute_area(min_lat, max_lat, min_long, max_long)

    vector = np.arange(simulations)
    density_levels = [1, 25, 50, 75, 100]
    selected_ues_per_density = []
    
    for ue_number in density_levels:
        selected_ues = np.random.choice(vector, ue_number, replace=False)
        selected_ues_per_density.append(selected_ues)

    strategy_density_results = {}
    for strategy, ue_results in raw_results.items():
        densities = []
        for selected_ues in selected_ues_per_density:
            results_subset = [ue_results[ue] for ue in selected_ues]
            merged_result = merge_results(results_subset)
            densities.append(merged_result)

        strategy_density_results[strategy] = densities

    plot_ho_rate_per_ue_density(strategy_density_results, time, area)'''

    ### USA Simulations with Starlink
    '''
    satellites1 = satellites_coverage("{}usa0.txt".format(data_path), start_time, end_time, time_step_seconds)
    satellites2 = satellites_coverage("{}usa1.txt".format(data_path), start_time, end_time, time_step_seconds)
    satellites3 = satellites_coverage("{}usa2.txt".format(data_path), start_time, end_time, time_step_seconds)
    satellites4 = satellites_coverage("{}usa3.txt".format(data_path), start_time, end_time, time_step_seconds)
    satellites5 = satellites_coverage("{}usa4.txt".format(data_path), start_time, end_time, time_step_seconds)
    satellites6 = satellites_coverage("{}usa5.txt".format(data_path), start_time, end_time, time_step_seconds)
    satellites7 = satellites_coverage("{}usa6.txt".format(data_path), start_time, end_time, time_step_seconds)
    satellites8 = satellites_coverage("{}usa7.txt".format(data_path), start_time, end_time, time_step_seconds)
    satellites9 = satellites_coverage("{}usa8.txt".format(data_path), start_time, end_time, time_step_seconds)
    satellites10 = satellites_coverage("{}usa9.txt".format(data_path), start_time, end_time, time_step_seconds)

    satellites = {
        "satellites1": satellites1,
        "satellites2": satellites2,
        "satellites3": satellites3,
        "satellites4": satellites4,
        "satellites5": satellites5,
        "satellites6": satellites6,
        "satellites7": satellites7,
        "satellites8": satellites8,
        "satellites9": satellites9,
        "satellites10": satellites10
    }

    plot_satellites_probabilities(satellites)

    plot_waiting_time(satellites)

    plot_service_time(satellites)

    plot_orbits_used("{}usa1.txt".format(data_path), "{}neighbors_starlink.txt".format(data_path2), start_time, end_time, time_step_seconds)
    '''
