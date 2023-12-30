import matplotlib.pyplot as plt
import math
import numpy as np
from palettable.colorbrewer.sequential import YlOrBr_3

G = 6.67e-20  # big G with km instead of m (divide by 10^9).

saturn_equatorial_radius = 60268  # km
saturn_polar_radius = 60268  # km 54364
saturn_mass = 5.683e26  # kg

mimas_radius = 198  # km
mimas_mass = 3.75e23  # kg (should be e19)
mimas_period = 23 * 60 * 60  # seconds
dist_mimas_saturn = 185539  # km
t = math.sqrt((4*np.pi**2 * dist_mimas_saturn**3)/(G*saturn_mass))  # correct mimas_period

titan_radius = 2575  # km
titan_mass = 1.35e23  # kg
titan_period = 16 * 24 * 60 * 60  # seconds
dist_titan_saturn = 1221870  # km

huygens_gap_radius = 117680  # km
huygens_gap_width = 350  # km


def graph():
    euler_error_r_values = []
    euler_error_v_values = []
    verlet_error_r_values = []
    verlet_error_v_values = []
    n=1000
    timesteps = []
    for i in range(1,n+1):
        timesteps.append(i)

    errors = read_errors()

    # fig = plt.figure()
    # ax = fig.add_subplot(2,1,1)
    plt.plot(timesteps, errors, color='r',label="Verlet")
    plt.yscale("log")
    plt.xscale("log")
    x = np.linspace(1,n,n)
    y = (errors[0])*x**2
    plt.plot(x,y,"b--",label="O(t^2)")
    plt.xlabel("timestep / s")
    plt.ylabel("% error")
    plt.legend()
    plt.title('Percent error against timestep for Mimas using Velocity Verlet method')
    plt.show()

def read_errors():
    errors = []
    with open("errors.txt") as f:
        for line in f:
            r = float(line)
            print(r)
            errors.append(r)
    
    return errors


graph()