#!/usr/bin/python3

import sys
import datetime
import numpy as np
import multiprocessing as mp
from pkg_resources import resource_filename
from satsim.sat import Sat
from satsim.app import SimulationApp, run_headless_simulation

# Define the path to the data directory
data_path = resource_filename("satsim", "data/")

def read_tle_file(file_path):
    """
    Reads a TLE file and returns a list of satellites with their TLE data.
    """
    satellites = []
    with open(file_path, 'r') as file:
        lines = file.readlines()
        for i in range(0, len(lines), 3):
            sat_name = lines[i].strip()
            tle_line1 = lines[i + 1].strip()
            tle_line2 = lines[i + 2].strip()
            satellites.append(Sat(name=sat_name, line1=tle_line1, line2=tle_line2))
    return satellites

# Read active satellites from the TLE file
Sats = read_tle_file(f"{data_path}active.txt")

def run_simulation(i, start_time, end_time, Sats):
    np.random.seed(i + int(datetime.datetime.now().timestamp()))

    # Random UE starting position - Equator
    """
    min_lat, max_lat = -0.4, 0.4
    min_long, max_long = -8.5, 8.5
    lat, long = np.random.uniform(min_lat, max_lat), np.random.uniform(min_long, max_long)
    """

    # Random UE starting position - North
    min_lat, max_lat = 64, 64.8
    min_long, max_long = -22, -13
    lat, long = np.random.uniform(min_lat, max_lat), np.random.uniform(min_long, max_long)

    constellation = "OneWeb"
    antenna_aperture = 0.4 if constellation == "OneWeb" else 1

    # Create the simulation app (headless)
    app = SimulationApp(
        id=i,
        Sats=Sats,
        start_time=start_time,
        end_time=end_time,
        time_step_seconds=0.05,
        terminal_coords=(long, lat, 0),
        terminal_speed=0,
        scenario="Suburban",
        frequency=12e9,
        max_gain_sat=40,
        antenna_aperture=antenna_aperture,
        gain_term=40,
        eirp_density=-13.4,  # dBW/4 kHz
        bandwidth=250e6,
        scs=120e3,           # Subcarrier Spacing
        guardband=5e6,
        simulation_satellites_filename="oneweb.txt",
        output_data_filename=f"ue_test_north{i}.txt",
        update_tle_file=False,
        constellation_name=constellation
    )

    # Run headless simulation
    print("Simulation " + "{:04}".format(i) + ": starting")
    run_headless_simulation(app)
    print("Simulation " + "{:04}".format(i) + ": ended")

def main():
    # Use ISO 8601 format for start and end datetimes. See https://docs.python.org/3/library/datetime.html#datetime.datetime.fromisoformat.
    datetime_start = '2025-03-26T01:00:00.000'
    datetime_end   = '2025-03-26T01:01:00.000'

    start_time = datetime.datetime.fromisoformat(datetime_start)
    end_time = datetime.datetime.fromisoformat(datetime_end)

    total_simulations = 100
    parallel_simulations = 4

    # Create a pool of simulations, running up to parallel_simulations simulation processes simultaneously
    with mp.Pool(processes=parallel_simulations) as pool:
        results = []
        # Launch all simulations
        for i in range(total_simulations):
            result = pool.apply_async(run_simulation, (i, start_time, end_time, Sats))
            results.append(result)

        # Wait for all simulations to finish
        for result in results:
            result.get()

if __name__ == "__main__":
    mp.freeze_support()
    main()