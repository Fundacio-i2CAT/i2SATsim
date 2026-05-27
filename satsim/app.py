import csv
import numpy as np
from datetime import datetime, timedelta
from pkg_resources import resource_filename
from satsim.sat import Sat
from satsim.channel_model import getElevation
from satsim.ue import get_next_position

class SimulationApp:
    def __init__(self, id, Sats, start_time, end_time, time_step_seconds, 
                 terminal_coords, terminal_speed, scenario, frequency, max_gain_sat, antenna_aperture, gain_term, eirp_density,
                 bandwidth, scs, guardband, simulation_satellites_filename, output_data_filename, update_tle_file, constellation_name=None):

        self.id = id

        self.Sats = Sats
        self.sortSats()

        self.start_time = start_time
        self.date = start_time
        self.end_time = end_time
        self.time_step_seconds = time_step_seconds

        self.visibleSats = []
        self.ableToConnectSats = []

        self.tle_path = resource_filename("satsim", "data/")
        self.all_satellites_filename = f"{self.tle_path}{simulation_satellites_filename}"
        self.constellation = constellation_name
        if update_tle_file:
            self.updateTleFile()
        self.all_sats = self.readSatsFromFile(self.all_satellites_filename)
        self.all_sats.sort(key=lambda s: s.name)

        self.data_filename = f"{resource_filename('satsim','output_data/')}{output_data_filename}"
        self.terminal_coords = terminal_coords
        self.terminal_speed = terminal_speed
        self.scenario = scenario
        self.frequency = frequency
        self.max_gain_sat = max_gain_sat
        self.antenna_aperture = antenna_aperture
        self.gain_term = gain_term
        self.eirp_density = eirp_density
        self.bandwidth = bandwidth
        self.scs = scs
        self.guardband = guardband

        self.beams = 16 if self.constellation == "OneWeb" else 1

    def sortSats(self):
        self.Sats.sort(key=lambda s: s.name)

    def updateSatellites(self):
        for sat in self.Sats:
            sat.updateOrbitalParameters(self.date)

    def readSatsFromFile(self, filename):
        sats = []
        with open(filename, 'r') as f:
            lines = f.readlines()
            for i in range(0, len(lines), 3):
                name = lines[i].strip()
                line1 = lines[i + 1].strip()
                line2 = lines[i + 2].strip()
                sats.append(Sat(name=name, line1=line1, line2=line2, date=self.date))
        return sats

    def updateTleFile(self):
        import requests
        url = f"https://celestrak.org/NORAD/elements/gp.php?GROUP={self.constellation}&FORMAT=tle"
        try:
            response = requests.get(url)
            response.raise_for_status()
            with open(self.all_satellites_filename, 'w') as file:
                file.write(response.text)
        except requests.exceptions.RequestException as e:
            print(f"Error fetching TLE data: {e}")

    def updateSimulation(self):
        self.visibleSats = []
        self.ableToConnectSats = []
        for sat in self.all_sats:
            lon, lat, h = sat.getPosition(self.date)
            if not np.isnan(lat) and not np.isnan(lon):
                if sat.is_visible(self.terminal_coords, self.date, self.constellation):
                    self.visibleSats.append(sat)
                if sat.is_under_coverage(self.date, self.terminal_coords):
                    self.ableToConnectSats.append(sat)

    def saveSimulationData(self):
        with open(self.data_filename, mode='a', newline='') as file:
            writer = csv.writer(file, delimiter=';')
            writer.writerow([
                self.date.strftime('%Y-%m-%d %H:%M:%S.%f'),
                '{' + ', '.join(sat.name for sat in self.visibleSats) + '}',
                '{' + ', '.join(f"({t[0]}, {t[1]}, {t[2]})" for sat in self.visibleSats for t in [sat.getPosition(self.date)]) + '}',
                '{' + ', '.join(sat.name for sat in self.ableToConnectSats) + '}',
                '{' + ', '.join(f"({', '.join(map(str, sat.covering_beams))})" for sat in self.ableToConnectSats) + '}',
                '{' + ', '.join(f"({', '.join(f'({pos[0]}, {pos[1]})' for i, pos in enumerate(sat.getBeamCenters()) if i in sat.covering_beams)})" for sat in self.ableToConnectSats) + '}',
                '{' + ', '.join(f"({', '.join(map(str, sat.computeRSRP(self.date, self.terminal_coords, self.scenario, self.time_step_seconds, self.frequency, self.max_gain_sat, self.antenna_aperture, self.gain_term, self.eirp_density, self.bandwidth, self.scs, self.guardband)))})" for sat in self.ableToConnectSats) + '}',
                '{' + f"{self.terminal_coords}" + '}'
            ])


def run_headless_simulation(app):
    current_time = app.start_time
    timestep = timedelta(seconds=app.time_step_seconds)

    # Prepare CSV
    with open(app.data_filename, mode='w', newline='') as file:
        writer = csv.writer(file, delimiter=';')
        writer.writerow(['Time', 'Visible_satellites', 'Satellite_positions', 'Possible_to_connect_satellites', 'Beams', 'Beam_centers', 'RSRP', 'Terminal_coords'])

    last_print_time = app.start_time
    print("Simulation " + "{:04}".format(app.id) + ": current_time " + f"{current_time.strftime('%Y-%m-%d %H:%M:%S')}")

    while current_time <= app.end_time:
        app.date = current_time
        app.updateSatellites()
        app.updateSimulation()
        app.saveSimulationData()
        app.terminal_coords = get_next_position(app.terminal_coords, app.terminal_speed, app.time_step_seconds)
        
        if current_time - last_print_time >= timedelta(minutes=1):
            print("Simulation " + "{:04}".format(app.id) + ": current_time " + f"{current_time.strftime('%Y-%m-%d %H:%M:%S')}")
            last_print_time = current_time
        
        current_time += timestep