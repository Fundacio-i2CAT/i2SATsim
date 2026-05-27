import random
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from datetime import timedelta
import numpy as np
import geopy.distance

def connectivity_window(start_time, end_time, time_step_seconds, min_time_connected, max_time_connected, avg_time_connected, min_time_sleeping, max_time_sleeping, avg_time_sleeping) -> list:
    """
    Defines the connectivity and sleeping windows of a UE for a given simulation.
    
    Parameters:
    - start_time: The start time of the simulation.
    - end_time: The end time of the simulation.
    - time_step_seconds: The time step in seconds for the simulation.
    - min_time_connected: Minimum time the UE can be in connected mode.
    - max_time_connected: Maximum time the UE can be in connected mode.
    - avg_time_connected: Average time the UE can be in connected mode.
    - min_time_sleeping: Minimum time the UE can be in sleeping mode.
    - max_time_sleeping: Maximum time the UE can be in sleeping mode.
    - avg_time_sleeping: Average time the UE can be in sleeping mode.
    """
    
    current_time = start_time
    ue_connectivity_window = []
    i = 1

    while current_time < end_time:
        if i == 1:
            time_random_connected = random.expovariate(1 / avg_time_connected)
            time_connected = max(min_time_connected, min(time_random_connected, max_time_connected))
            ue_connectivity_window.extend([i] * int(time_connected/time_step_seconds))
            i -= 1
            current_time += timedelta(seconds=time_connected)
        elif i == 0:
            time_random_sleeping = random.expovariate(1 / avg_time_sleeping)
            time_sleeping = max(min_time_sleeping, min(time_random_sleeping, max_time_sleeping))
            ue_connectivity_window.extend([i] * int(time_sleeping/time_step_seconds))
            i += 1
            current_time += timedelta(seconds=time_sleeping)
    seconds = int((end_time - start_time).total_seconds()/time_step_seconds)
    ue_connectivity_window = ue_connectivity_window[:seconds]
    return ue_connectivity_window

def plot_ue_connectivity_window(ue_connectivity_window, start_time, time_step_seconds):
    """
    Plots the connectivity window (1 for connected, 0 for sleeping) over time.
    
    Parameters:
    - ue_connectivity_window: A list of 1s and 0s representing connected and sleeping states.
    - start_time: The start time of the simulation.
    - time_step_seconds: The time step in seconds for the simulation.
    """
    
    # Create a list of time values corresponding to each connectivity state
    time_values = [start_time + timedelta(seconds=i*time_step_seconds) for i in range(len(ue_connectivity_window))]

    # Prepare the figure and axis
    fig, ax = plt.subplots(figsize=(12, 6))

    # Plot the connectivity window as a step plot
    ax.step(time_values, ue_connectivity_window, where='post', color='red', linewidth=2)

    # Set y-ticks to show 'Connected' and 'Sleeping'
    ax.set_yticks([0, 1])
    ax.set_yticklabels(['Sleeping', 'Connected'])

    # Set x-axis limits to the range of time values
    ax.set_xlim([time_values[0], time_values[-1]])

    # Format the x-axis to show time
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M:%S.%f'))
    ax.xaxis.set_minor_formatter(mdates.DateFormatter('%H:%M'))

    # Rotate x-axis labels for better readability
    plt.xticks(rotation=45)

    # Labels and title
    plt.xlabel('Time')
    plt.ylabel('UE State')
    plt.title('UE Connectivity and Sleeping Periods Over Time')

    # Grid and layout
    plt.grid(True)
    plt.tight_layout()
    
    # Show the plot
    plt.show()

def get_next_position(terminal_coords, terminal_speed, time_step_seconds):
    distance = terminal_speed*time_step_seconds
    azimuth = np.random.uniform(0, 360)
    point = geopy.distance.distance(meters=distance).destination((terminal_coords[1], terminal_coords[0]), bearing=azimuth)
    next_position = (point.longitude, point.latitude, 0)
    return next_position
