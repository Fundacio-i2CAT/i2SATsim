# <img src='assets/i2satsim_logo.png' alt="i2SATsim's logo" width='50%' />

i2SATsim is a satellite network simulator to study how 3GPP NTN UEs connect to satellites and how **handover events** perform as the orbits propagate.

## Architecture

The simulator is structured in **two main phases**:

1. **Data generation** – running the simulation
2. **Data analysis and visualization** – plotting results.

i2SATsim is inspired by tools like [Gpredict](https://github.com/csete/gpredict), [Orbitron](http://www.stoff.pl/), and [Pypredict](https://github.com/spel-uchile/Pypredict), on which —originally— it was largely based (see details [below](#Licensing)).

## Features

The following are some of the key features of i2SATsim:

- **Real-Time Tracking**: Displays the position and orbital parameters of satellites in real time.
- **TLE Data Management**:
  - Updates TLE data from [CelesTrak](https://www.celestrak.com/NORAD/elements/).
  - Allows the user to open TLE data from a text file.
  - Saves the TLE data of displayed satellites into a new text file.
- **Satellite Simulation**:
  - Simulates satellite deployments based on varying velocities, considering the mass of the satellites.
  - Supports configuration of start and end times for simulations.
- **Coverage and Anomalies**: (WiP)

## Installation

To install i2SATsim and run the simulation, follow these steps:

1. **Clone the repository:**

   ```bash
    git clone https://gitlab.i2cat.net/areas/mwi/handovers/i2satsim.git
    cd i2satsim
   ```

2. **Install all the dependencies**

Before you begin, make sure you have:

- Python 3.7 or higher installed
- Pip (Python’s package manager)

  Install the required packages with:

  ```bash
  pip install -r requirements.txt
  ```

## Usage

### Data generation

To generate simulation data, run the following command:

```bash
python3 -m satsim
```

Before running, make sure to configure the desired parameters in your main script:

| **Parameter**                    | **Purpose**                                       |
| ----------------------------------| ---------------------------------------------------|
| `terminal_coords`                | Terminal location (longitude, latitude, altitude) |
| `terminal_speed`                 | Terminal speed in m/s (0 = static)                |
| `scenario`                       | Propagation environment                           |
| `frequency`                      | Carrier frequency in Hz                           |
| `eirp_density`                   | Satellite EIRP density (dBW/4 kHz)                |
| `bandwidth`                      | Bandwidth of the signal in Hz                     |
| `scs`                            | Subcarrier spacing                                |
| `simulation_satellites_filename` | TLE file containing satellite constellation data  |
| `output_data_filename`           | Output filename for the generated simulation data |
| `constellation_name`             | Name of the satellite constellation               |
| `update_tle_file`                | Whether to fetch new TLEs                         |

This will generate an output file in output_data folder (e.g., `ue_north0.txt`) with the simulation results for each configuration.

---

### Analysis and plotting

Once data has been generated, you can analyze and visualize the results by running:

```bash
python3 plot.py
```

Before doing so, configure the simulation strategies in `satsim/plot.py`:

| **Parameter**            | **Purpose**                                 |
| --------------------------| ---------------------------------------------|
| `strategies_names`       | Labels for each strategy being compared     |
| `strategy`               | Handover logic type                         |
| `a3_offset`              | Offset for triggering A3 event              |
| `a3_hys`                 | Hysteresis for A3 handover                  |
| `a4_threshold`           | Threshold for A4 event                      |
| `ttt`                    | Time-to-Trigger duration (A3)               |
| `num_neighbors`          | Number of neighbor satellites to consider   |
| `k`                      | Number of top satellites used for selection |
| `ue_connectivity_window` | Optional window for evaluating connectivity |

Then choose one or more of the following plots depending on the analysis you want to perform:

- `plot_effective_time`: Shows the duration of valid connectivity for each strategy.
- `plot_events`: Displays handover and disconnection events over time.
- `plot_call_drop_probability`: Overall probability of call drops.
- `plot_call_drop_distribution_by_zone`: Geographic distribution of dropped calls.
- `plot_distance_ho_cdf` and `plot_distance_cdf`: Cumulative distribution of handover distances and distances in general.
- `plot_rsrp_boxplot`: RSRP (Reference Signal Received Power) comparison.
- `plot_tos_cdf`: Time-of-stay (duration connected to the same satellite) distribution.

> **Note**: You can find examples of the plots that can be generated in `satsim/plot.py`.

## Licensing

This software is licensed under the GNU Affero General Public License v3.

The repository includes code from [Pypredict](https://github.com/spel-uchile/Pypredict) (the whole `satsim/node.py` and most of `satsim/sat.py` files), which is licensed with GPL 3.0. In accordance with Section 13 of both licenses, the combined work as a whole is distributed under the terms of the AGPL 3.0.

## Authorsip and citations

This code was authored by Anna Esteve and Pau Feixa, from the [Mobile Wireless Internet](https://i2cat.net/research-group/mobile-wireless-internet/) research group at [Fundació i2CAT](https://www.i2cat.net) in Barcelona, Catalunya.

You can cite our research as:

```LaTeX
@INPROCEEDINGS{2026001343,
  author={Feixa, Pau and Esteve, Anna and {Pueyo Centelles}, Roger and Camps-Mur, Daniel}
  booktitle={2026 IEEE International Conference on Communications (ICC)},
  title={Never Look Back: Optimizing 5G Handovers in a OneWeb-like LEO NTN Constellation},
  year={2026},
  keywords={5G, NTN, Mobility, Handover, 3GPP},
  doi={TBC}
}
```

## Acknowledgments

This work was co-funded by the Government of Catalonia, in the scope of the [NewSpace Strategy Programme for Catalonia](https://politiquesdigitals.gencat.cat/ca/economia/estrategia-new-space-catalunya/index.html), and by the European Union’s Horizon Europe under [Grant Agreement no. 101096526 – ETHER](https://ether-project.eu/).

## Miscellaneous

i2SATsim is stylized with uppercase SAT, or as **<span style="color:#f25623ff">i2</span><span style="color:black">sat</span><span style="color:#999999ff">sim</span><sup><span style="background-color:black"><span style="color:white">🛰</span></span></sup>**, mimicking i2CAT's **<span style="color:#f25623ff">i2</span><span style="color:black">cat</span><sup><span style="background-color:black"><span style="color:white">R</span></span></sup>**'s branding.

> **Note**: Coloured Markdown may not be rendered by GitLab/GitHub.
