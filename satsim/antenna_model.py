import numpy as np
import matplotlib.pyplot as plt
from scipy.special import jn
from fractions import Fraction

def get_gain(frequency, theta, max_gain, D):
    c = 3e8
    k = 2*np.pi*frequency/c
    a = D/2
    aux = k*a*np.sin(np.deg2rad(theta))
    if theta == 0:
        gain = 0
    else:
        gain = 10*np.log10(4*np.pow(jn(1, aux)/aux, 2))

    return gain + max_gain

def get_transmit_power(eirp_density, bandwidth, max_gain):
    pt = eirp_density + 10*np.log10(bandwidth/4e3) - max_gain + 30 #+30 to convert to dBm
    
    return pt

def plot_antenna_pattern(frequency, max_gain, D):
    thetas = np.arange(-90, 90.5, 0.5)

    gain = np.zeros_like(thetas, dtype=float)
    for i, theta in enumerate(thetas):
        gain[i] = get_gain(frequency, theta, max_gain, D)

    aperture_radius = Fraction(D*frequency/(2*c)).limit_denominator()
    
    plt.plot(thetas, gain)
    #plt.title(f"Aperture radius = {aperture_radius}*lambda", fontsize = 20)
    plt.title(f"Aperture radius = {D/2} m", fontsize = 20)
    plt.xlabel("Angle (°)", fontsize = 18)
    plt.ylabel("Antenna gain (dB)", fontsize = 18)
    plt.xlim([-90, 90])
    plt.ylim([max_gain-80, max_gain])
    plt.grid()
    plt.tick_params(axis='both', which='major', labelsize=16)
    plt.show()

c = 3e8
frequency = 12e9

#Example from the 3GPP Standard
max_gain = 0
D = 2*10*c/frequency
#Example from the Nokia paper
max_gain = 40
D = 0.4
#plot_antenna_pattern(frequency, max_gain, D)
