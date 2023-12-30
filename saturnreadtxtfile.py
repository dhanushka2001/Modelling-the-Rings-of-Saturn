import matplotlib.pyplot as plt
import math
import numpy
import numpy as np
from palettable.colorbrewer.sequential import YlOrBr_3
from matplotlib.widgets import Slider, Button, RadioButtons
from matplotlib.pyplot import figure
from PIL import Image
from matplotlib import pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.offsetbox import AnchoredText
from matplotlib import patches
import time

G = 6.67e-20  # big G with km instead of m (divide by 10^9).

saturn_equatorial_radius = 60268  # km
saturn_polar_radius = 60268  # km 54364
saturn_mass = 5.683e26  # kg

mimas_radius = 198  # km
mimas_mass = 3.75e19  # kg (should be e19)
mimas_period = 23 * 60 * 60  # seconds
dist_mimas_saturn = 185539  # km
t = math.sqrt((4*np.pi**2 * dist_mimas_saturn**3)/(G*saturn_mass))  # correct mimas_period

titan_radius = 2575  # km
titan_mass = 1.35e23  # kg
titan_period = 16 * 24 * 60 * 60  # seconds
dist_titan_saturn = 1221870  # km

huygens_gap_radius = 117680  # km
huygens_gap_width = 350  # km

centre = [0, 0, 0]

def set_axes_equal(ax):
    '''Make axes of 3D plot have equal scale so that spheres appear as spheres,
    cubes as cubes, etc..  This is one possible solution to Matplotlib's
    ax.set_aspect('equal') and ax.axis('equal') not working for 3D.

    Input
      ax: a matplotlib axis, e.g., as output from plt.gca().
    '''

    x_limits = ax.get_xlim3d()
    y_limits = ax.get_ylim3d()
    z_limits = ax.get_zlim3d()

    x_range = abs(x_limits[1] - x_limits[0])
    # x_middle = np.mean(x_limits)
    x_middle = centre[0]
    y_range = abs(y_limits[1] - y_limits[0])
    # y_middle = np.mean(y_limits)
    y_middle = centre[1]
    z_range = abs(z_limits[1] - z_limits[0])
    # z_middle = np.mean(z_limits)
    z_middle = centre[2]

    # The plot bounding box is a sphere in the sense of the infinity
    # norm, hence I call half the max range the plot radius.
    plot_radius = 0.5 * max([x_range, y_range, z_range])

    ax.set_xlim3d([x_middle - plot_radius, x_middle + plot_radius])
    ax.set_ylim3d([y_middle - plot_radius, y_middle + plot_radius])
    ax.set_zlim3d([z_middle - plot_radius, z_middle + plot_radius])

def plot():
    plt.style.use('dark_background')
    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot(111, projection='3d')

    # np.linspace returns evenly spaced numbers from (start, stop, num = 100)
    u = np.linspace(0, 2 * np.pi, 100)
    v = np.linspace(0, np.pi, 100)

    # spherical to cartesian coordinates
    # np.outer returns the outer product
    x = saturn_equatorial_radius * np.outer(np.cos(u), np.sin(v))
    y = saturn_equatorial_radius * np.outer(np.sin(u), np.sin(v))
    z = saturn_polar_radius * np.outer(np.ones(np.size(u)), np.cos(v))

    # mimas
    x1 = mimas_radius * np.outer(np.cos(u), np.sin(v)) + dist_mimas_saturn
    y1 = mimas_radius * np.outer(np.sin(u), np.sin(v))
    z1 = mimas_radius * np.outer(np.ones(np.size(u)), np.cos(v))

    # titan
    x2 = titan_radius * np.outer(np.cos(u), np.sin(v)) + dist_titan_saturn
    y2 = titan_radius * np.outer(np.sin(u), np.sin(v))
    z2 = titan_radius * np.outer(np.ones(np.size(u)), np.cos(v))

    # titan atmosphere
    x3 = (titan_radius + 600) * np.outer(np.cos(u), np.sin(v)) + dist_titan_saturn
    y3 = (titan_radius + 600) * np.outer(np.sin(u), np.sin(v))
    z3 = (titan_radius + 600) * np.outer(np.ones(np.size(u)), np.cos(v))

    x4 = (titan_radius + 300) * np.outer(np.cos(u), np.sin(v)) + dist_titan_saturn
    y4 = (titan_radius + 300) * np.outer(np.sin(u), np.sin(v))
    z4 = (titan_radius + 300) * np.outer(np.ones(np.size(u)), np.cos(v))

    x5 = (titan_radius + 150) * np.outer(np.cos(u), np.sin(v)) + dist_titan_saturn
    y5 = (titan_radius + 150) * np.outer(np.sin(u), np.sin(v))
    z5 = (titan_radius + 150) * np.outer(np.ones(np.size(u)), np.cos(v))

    # ax.set_box_aspect((1, 1, saturn_polar_radius/saturn_equatorial_radius))
    # ax.set_box_aspect((1,1,1))
    norm = plt.Normalize()
    facecolors = YlOrBr_3.mpl_colormap(norm(z * z))
    ax.plot_surface(x, y, z, rstride=3, cstride=3, facecolors=facecolors, linewidth=0.1, color=(247 / 256, 213 / 256, 152 / 256))
    ax.plot_surface(x1, y1, z1, rstride=3, cstride=3, color=(122 / 256, 127 / 256, 128 / 256))
    ax.plot_surface(x2, y2, z2, rstride=3, cstride=3, facecolors=facecolors, linewidth=0.1,
                    color=(191 / 256, 158 / 256, 91 / 256))
    ax.plot_surface(x3, y3, z3, rstride=3, cstride=3, color=(0 / 256, 0 / 256, 255 / 256, 0.02))
    ax.plot_surface(x4, y4, z4, rstride=3, cstride=3, color=(30 / 256, 144 / 256, 255 / 256, 0.05))
    ax.plot_surface(x5, y5, z5, rstride=3, cstride=3, color=(0 / 256, 255 / 256, 255 / 256, 0.1))

    ax.axis('off')
    set_axes_equal(ax)
    # azim = 0 to 360 around the x-y plane
    ax.azim = 180
    # dist = 0 to inf, distance from Saturn
    ax.dist = 3
    # elev = 0 to 180, angle between eye and x-y plane
    ax.elev = 90

    read_mimas(ax)
    read_particles(ax)

    theta = np.linspace(0, 2 * np.pi, 201)
    xx = (huygens_gap_radius+huygens_gap_width/2) * np.cos(theta)
    yy = (huygens_gap_radius+huygens_gap_width/2) * np.sin(theta)
    z = np.zeros_like(theta)
    ax.plot(xx, yy, z, 'w--')
    xx = (huygens_gap_radius - huygens_gap_width / 2) * np.cos(theta)
    yy = (huygens_gap_radius - huygens_gap_width / 2) * np.sin(theta)
    z = np.zeros_like(theta)
    ax.plot(xx, yy, z, 'w--')


    plt.show()

timestep = 100
subcycle = 500
n_orbits = 8000
n_particles = 10000
n_positions = int((t*n_orbits)/(timestep*subcycle))
min_radius = 116382
max_radius = 117382

def read_mimas(ax):
    with open("positions_mimas.txt") as f:
        n = 0
        pos = []
        for line in f:
            n+=1
            number_of_lines = int(n_orbits*t/timestep)
            x = int(t/timestep)
            if n >= number_of_lines - x - 1:
                r = line.split()
                R = [float(x) for x in r]
                pos.append(R)
        xx = np.array(pos)
        ax.plot(xx[:, 0], xx[:, 1], xx[:, 2], 'r-')

# def read_particles(ax):
#     with open("positions_particles.txt") as f:
#         n = 0
#         i=0
#         pos = []
#         for line in f:
#             n+=1
#             r = line.split()
#             R = [float(x) for x in r]
#             pos.append(R)
#             x = int(t/timestep)
#             if n % x == 0:
#                 i+=1
#                 zz = np.array(pos)
#                 ax.plot(zz[:, 0], zz[:, 1], zz[:, 2], 'c-', linewidth=0.1)
#                 pos = []
#                 print(str(i) + "/" + str(n_particles),end="\r")

def read_particles(ax):
    # for particle_number in range(n_particles):
    #     line_number = 0
    #     pos = [None]*int(t/timestep)
    #     with open("positions_particles.txt") as f:
    #         for line in f:
    #             #line number
    #             line_number += 1
    #             if line_number % n_particles == (particle_number+1):
    #                 r = line.split()
    #                 R = [float(x) for x in r]
    #                 n = int((line_number-(particle_number+1))/n_particles)
    #                 pos[n] = R
    #                 #check if at the last position
    #                 if line_number == n_particles*int(t/timestep) - (n_particles - (particle_number+1)):
    #                     zz = np.array(pos)
    #                     ax.plot(zz[:, 0], zz[:, 1], zz[:, 2], 'c-', linewidth=0.1)
    #                     print(str(particle_number+1) + "/" + str(n_particles),end="\r")
    #                     #dont need to reset pos, break will go back to the first for loop and pos will be reset there
    #                     break

    # with open("positions_particles.txt") as f:
    #     pos = np.empty([int(t/timestep),n_particles,3])
    #     pos1 = np.empty([n_particles,3])
    #     line_number = 0
    #     position_number = 0 #[0,815-1]
    #     for line in f:
    #         line_number += 1
    #         r = line.split()
    #         R = [float(x) for x in r]
    #         pos1[(line_number-1)%n_particles] = R
    #         if line_number%n_particles == 0:
    #             pos[position_number] = pos1
    #             position_number += 1
    #             print(str(position_number) + "/" + str(int(t/timestep)),end="\r")
    
    # #transpose in specific way so that individual particles are together instead of which position_number.
    # #read more here: https://bic-berkeley.github.io/psych-214-fall-2016/numpy_transpose.html
    # pos = pos.transpose(1, 0, 2)
    # particle_number = 0
    # for i in pos:
    #     #we want just the LAST ORBIT.
    #     r_particle = saturn_equatorial_radius + (dist_mimas_saturn - saturn_equatorial_radius) * (particle_number+1) / (n_particles + 1)
    #     t_particle = t * math.sqrt((r_particle/dist_mimas_saturn)**3)
    #     number_of_positions = int(t_particle/timestep) #no. of positions for ONE ORBIT
    #     i = i[-1*number_of_positions:]
    #     ax.plot(i[:, 0], i[:, 1], i[:, 2], 'c-', linewidth=0.1) 
    #     particle_number += 1

    with open("positions_particles.txt") as f:
        pos = np.empty([int(t/(timestep*subcycle)),n_particles,3])
        pos1 = np.empty([n_particles,3])
        line_number = 0
        position_number = 0 #[0,815-1]
        for line in f:
            line_number += 1
            r = line.split()
            R = [float(x) for x in r]
            pos1[(line_number-1)%n_particles] = R
            if line_number%n_particles == 0:
                pos[position_number] = pos1
                position_number += 1
                print(str(position_number) + "/" + str(int(t/(timestep*subcycle))),end="\r")
    
    pos = pos.transpose(1, 0, 2)
    particle_number = 0
    #subcycle = 10
    for i in pos:
        r_particle = saturn_equatorial_radius + (dist_mimas_saturn - saturn_equatorial_radius) * (particle_number+1) / (n_particles + 1)
        t_particle = t * math.sqrt((r_particle/dist_mimas_saturn)**3)
        number_of_positions = int(t_particle/(subcycle*timestep)) #no. of positions for ONE ORBIT
        i = i[-1*number_of_positions:]
        ax.plot(i[:, 0], i[:, 1], i[:, 2], 'c-', linewidth=0.1)
        particle_number += 1


#n_orbits = 100
x = np.linspace(1,n_orbits,n_orbits)
def line_graph(ax):
    with open("positions_particles.txt") as f:
        n = 0
        i = 0
        p = 3060
        switch = 0
        pos = []
        for line in f:
            n+=1
            if switch == 1 and (n-1)%100 != 0:
                continue
            elif float(line) > 1000000:
                pos = []
                switch = 1
                continue
            else:
                pos.append(float(line))
                switch = 0
            
            #pos.append(float(line))
            if n % 100 == 0:
                switch = 0
                i += 1
                y = np.array(pos)
                ax.plot(x,y,"r-",linewidth=0.1)
                pos=[]
                print(str(i) + "/" + str(n_particles),end="\r")

#line_graph(ax)              

n = int((t*n_orbits)/(timestep*subcycle))
x = np.linspace(1,n,n) #(1,n,n)
def line_graph_subcycle(ax):
    # E:/Warwick/saturn/positions_particles_janus_epimetheus_massx1_95400km_dt=100_n_orb=10000_sub=10000.txt
    with open("E:/Warwick/saturn/positions_particles_huygens_mimas_massx1_116382-117382km_dt=100_n_orb=8000_sub=500_ellipse_parametric.txt") as f:
        pos = np.empty([int((t*n_orbits)/(timestep*subcycle)),n_particles])
        pos1 = np.empty([n_particles])
        line_number = 0
        position_number = 0 #[0,815-1]
        for line in f:
            line_number += 1
            
            #magnitude of vector
            r = line.split()
            R = [float(x) for x in r]
            #del R[-1]
            pos1[(line_number-1)%n_particles] = R[0]#np.linalg.norm(R)

            #pos1[(line_number-1)%n_particles] = float(line) #np.linalg.norm(R)
            if line_number%n_particles == 0:
                pos[position_number] = pos1
                position_number += 1
                # print(str(position_number) + "/" + str(int((t*n_orbits)/(timestep*subcycle))),end="\r")
    
    #pos = pos.transpose(1, 0, 2)
    pos = pos.T
    # pos is the list of positions for each particle, so y is an individual particles' positions
    position_number = 0
    print("\n")
    for y in pos:
        if any((i > 200000 or i < saturn_equatorial_radius) for i in y):
            continue
        else:
            ax.plot(x, y, 'b-', linewidth=0.02)
        position_number+=1
        # print(str(position_number) + "/" + str(n_particles),end="\r")

#line_graph_subcycle(ax)

def histogram_subcycle():
    # creates a list of particles' positions that are inside a good region.
    with open("E:/Warwick/saturn/positions_particles_huygens.txt") as f:
        pos = np.empty([int((t*n_orbits)/(timestep*subcycle))+1,n_particles])
        pos1 = np.empty([n_particles])
        line_number = 0
        position_number = 0 #[0,815-1]
        for line in f:
            line_number += 1
            
            #magnitude of vector
            r = line.split()
            R = [float(x) for x in r]
            del R[-1]
            pos1[(line_number-1)%n_particles] = np.linalg.norm(R)

            #pos1[(line_number-1)%n_particles] = float(line) #np.linalg.norm(R)
            if line_number%n_particles == 0:
                pos[position_number] = pos1
                position_number += 1
                print(str(position_number) + "/" + str(int((t*n_orbits)/(timestep*subcycle))),end="\r")
    
    #pos = pos.transpose(1, 0, 2)
    pos = pos.T
    # pos is the list of positions for each particle, so y is an individual particles' positions
    position_number = 0
    print("\n")
    # for y in pos:
        # remove particles that go away from good region
        # if any((i > 200000 or i < saturn_equatorial_radius) for i in y):
        #     np.delete(pos, position_number)
        # position_number+=1
        # print(str(position_number) + "/" + str(n_particles),end="\r")
    # return list of positions/times rather than particles
    return pos.T

#histogram_subcycle()

def plot_line_graph():
    fig = plt.figure()
    ax = plt.axes()
    #y = np.linspace(dist_mimas_saturn,dist_mimas_saturn,n)
    #ax.plot(x,y,"r-",linewidth=2)
    #yy = np.linspace(huygens_gap_radius,huygens_gap_radius,n)
    #ax.plot(x,yy,"r-",linewidth=1)
    #val = dist_mimas_saturn
    #a1 = (1/2)
    #a2 = (1/3)
    #a3 = (1/4)
    #a4 = (2/3)
    #a5 = (3/4)
    #a = [a1,a2,a3,a4,a5]
    #for i in a:
    # for i in range(2,22):
    #     j = 1 - 1/i
    #     j = math.pow(j,2/3) * dist_mimas_saturn
    #     yy = np.linspace(j,j,n)
    #     ax.plot(x,yy,"r-",linewidth=1)
    #     if i > 8 and i % 2 != 0:
    #         j = 1 - 2/i
    #         j = math.pow(j,2/3) * dist_mimas_saturn
    #         yy = np.linspace(j,j,n)
    #         ax.plot(x,yy,"g-",linewidth=1)
        
    #     if i > 10 and i % 3 != 0:
    #         j = 1 - 3/i
    #         j = math.pow(j,2/3) * dist_mimas_saturn
    #         yy = np.linspace(j,j,n)
    #         ax.plot(x,yy,"b-",linewidth=1)


    line_graph_subcycle(ax)
    ax.set_xlabel("Time (500*100s)".format(subcycle,timestep))
    ax.set_ylabel("Distance from Saturn (km)")
    ax.set_ylim([min_radius, max_radius])
    #ax.set_title("{0} particles, Mimas mass x1000, ")
    #plt.show()
    plt.savefig("E:/Pictures/saturnringdensitygraph_mimas_ellipse_parametric_n_orb=8000_sub=500_zoom_linegraph.png", bbox_inches="tight", dpi=1000)

plot_line_graph()

def plot_histogram():
    pos = histogram_subcycle()
    # (180,000 - 60,000) / 50 = 2400
    n_bins = 1000
    y_max = 17.5
    for i in range(n_positions):
        plt.clf()
        plt.hist(pos[i], range=[110000,124000], density=False, bins=n_bins, alpha=0.5, color = "cornflowerblue")
        plt.ylim([0, y_max])
        plt.xlabel("Distance from Saturn / km");
        plt.ylabel("# of particles")
        plt.title("Histogram of # of particles against distance, 10,000 particles")
        # plt.show()
        k = 1/2
        k = math.pow(k,2/3) * dist_mimas_saturn
        plt.vlines(x=k, ymin=0, ymax=y_max, colors='red', ls='--', lw=0.2)#, label="l=m")
        plt.savefig("E:/Pictures/huygens3/foo{}.jpg".format(str(i).zfill(3)), bbox_inches="tight", dpi=200)
        print(str(i+1) + "/" + str(n_positions),end="\r")

#plot_histogram()

def plot_scatter_graph():
    n_positions = int((t*n_orbits)/(timestep*subcycle)) # =815
    with open("positions_particles_x_y.txt") as f:
        pos = np.empty([n_positions,n_particles,3])
        pos1 = np.empty([n_particles,3])
        line_number = 0
        position_number = 0 #[0,815-1]
        for line in f:
            line_number += 1
            r = line.split()
            R = [float(x) for x in r]
            pos1[(line_number-1)%n_particles] = R
            if line_number%n_particles == 0:
                pos[position_number] = pos1
                position_number += 1
                print(str(position_number) + "/" + str(n_positions),end="\r")

    print("\n")
    position_number = 0
    for pos_list in pos:
        particle_number = 0
        for particle_list in pos_list:
            particle_list = particle_list[:-1]
            mag = np.sqrt(particle_list.dot(particle_list))
            if mag > dist_mimas_saturn + 20000:
                pos_list = np.delete(pos_list, particle_number)
            particle_number += 1
        position_number += 1
        print(str(position_number) + "/" + str(n_positions),end="\r")

    #N = 1 #[0,814] because len = 815
    
    x = []
    y = []

    for i in range(n_positions):
        length = len(pos[i])
        x1 = np.empty([length])
        y1 = np.empty([length])
        for j in range(length): # if not for mag > saturn etc. could just put range(n_particles)
            x1[j] = pos[i][j][0]
            y1[j] = pos[i][j][1]
        x.append(x1)
        y.append(y1)

    
    fig, ax = plt.subplots()
    figure(1, figsize=(6, 6), dpi=80)
    scatter = ax.scatter(x[0],y[0], s=5)
    ax.set_xlabel("x (km)")
    ax.set_ylabel("y (km)")

    plt.subplots_adjust(left = 0.25, bottom=0.25)

    axpos = plt.axes([0.25, 0.1, 0.65, 0.03])
    pos_slider = Slider(ax=axpos,
                        label="Time (500*100s)",
                        valmin=0,
                        valmax=814,
                        valinit=0,
                        valstep=1
    )

    def update(val):
        #scatter.set_xdata(x[pos_slider.val])
        #scatter.set_ydata(y[pos_slider.val])
        ax.cla()
        ax.scatter(x[pos_slider.val],y[pos_slider.val], s=5)
        fig.canvas.draw_idle()
        ax.set_xlim([-200000,200000])
        ax.set_ylim([-200000,200000])
    
    pos_slider.on_changed(update)

    plt.show()

#plot_scatter_graph()

def img_scatter_graph():
    n_positions = int((t*n_orbits)/(timestep*subcycle)) # =815
    print("Reading text file -> 2d list...")
    with open("E:/Warwick/saturn/positions_particles_x_y.txt") as f:
        pos = np.empty([n_positions,n_particles,3])
        pos1 = np.empty([n_particles,3])
        line_number = 0
        position_number = 0 #[0,815-1]
        for line in f:
            line_number += 1
            r = line.split()
            R = [float(x) for x in r]
            pos1[(line_number-1)%n_particles] = R
            if line_number%n_particles == 0:
                pos[position_number] = pos1
                position_number += 1
                print(str(position_number) + "/" + str(n_positions),end="\r")

    #print("\n")
    position_number = 0
    print("\nRemoving z axis and stray data points...")
    for pos_list in pos:
        particle_number = 0
        for particle_list in pos_list:
            particle_list = particle_list[:-1]
            mag = np.sqrt(particle_list.dot(particle_list))
            if mag > 350000:
                pos_list = np.delete(pos_list, particle_number)
            particle_number += 1
        position_number += 1
        print(str(position_number) + "/" + str(n_positions),end="\r")

    #N = 1 #[0,814] because len = 815
    
    #print("\n")
    x = []
    y = []

    print("\nCreating list of x and y values for scatter plot...")
    for i in range(n_positions):
        length = len(pos[i])
        x1 = np.empty([length])
        y1 = np.empty([length])
        for j in range(length): # if not for mag > saturn etc. could just put range(n_particles)
            x1[j] = pos[i][j][0]
            y1[j] = pos[i][j][1]
        x.append(x1)
        y.append(y1)
        print(str(i+1) + "/" + str(n_positions),end="\r")

    
    print("\nCreating scatter plot and saving as jpg...")
    fig, ax = plt.subplots() # very important this is not inside the for loop
    for i in range(n_positions):
        figure(1, figsize=(6, 6), dpi=80)
        ax.cla()
        r = dist_mimas_saturn
        t1 = 2*math.pi*math.sqrt(r*r*r/(G*saturn_mass))
        k = - (timestep * subcycle / t1 - 12)#1 - timestep * subcycle / t1
        c = math.cos((i)*2*math.pi*k)
        s = math.sin((i)*2*math.pi*k)

        X = x[i]*c - y[i]*s
        Y = x[i]*s + y[i]*c
        ax.scatter(X,Y, s=0.1, marker='o', c="blue", edgecolors='none')

        b = timestep/t
        X1 = r*math.cos(2*math.pi*i*subcycle*b)
        Y1 = r*math.sin(2*math.pi*i*subcycle*b)
        X2 = X1*c - Y1*s
        Y2 = X1*s + Y1*c
        ax.scatter(X2,Y2, s=20, c="red")

        ax.set_xlabel("x axis (km)")
        ax.set_ylabel("y axis (km)")

        plt.subplots_adjust(left = 0, bottom=0)
        ax.set_xlim([-200000,200000])
        ax.set_ylim([-200000,200000])

        fig.set_size_inches(7, 7)
        plt.savefig("E:/Pictures/saturn_mimas/foo{}.jpg".format(str(i).zfill(3)), bbox_inches="tight", dpi=300)
        print(str(i+1) + "/" + str(n_positions),end="\r")
    print("\n")

#img_scatter_graph()

def radial_velocity_subcycle():
    n_positions = int((t*n_orbits)/(timestep*subcycle)) # =815
    print("Reading txt file and creating pos and vel lists...")
    with open("E:/Warwick/saturn/positions_particles_huygens_mimas_massx1_115500-118000km_dt=100_n_orb=10000_sub=10000_ellipse_parametric.txt") as f:
        pos = np.empty([n_positions,n_particles])
        pos1 = np.empty([n_particles])
        vel = np.empty([n_positions,n_particles])
        vel1 = np.empty([n_particles])
        line_number = 0
        position_number = 0 #[0,815-1]
        for line in f:
            line_number += 1
            r = line.split()
            R = [float(x) for x in r]
            pos1[(line_number-1)%n_particles] = R[0]
            vel1[(line_number-1)%n_particles] = abs(R[1]*1000) #add abs() if want
            if line_number%n_particles == 0:
                pos[position_number] = pos1
                vel[position_number] = vel1
                position_number += 1
                print(str(position_number) + "/" + str(n_positions),end="\r")
    
    #pos = pos.transpose(1, 0, 2)
    pos = pos.T
    vel = vel.T
    # pos is the list of positions for each particle, so y is an individual particles' positions
    
    print("Removing stray particles...")
    particle_number = 0
    print("\n")
    
    #revert back to list of times rather than particles
    pos = pos.T
    vel = vel.T

    print("\nCreating scatter plot and saving as jpg...")
    x0 = 115500
    x1 = 118000
    N = 5 # particles per bin
    bin_width = N*(x1-x0)/n_particles # 100
    # bin_left_edge = np.arange(saturn_equatorial_radius,dist_mimas_saturn,bin_width)
    bin_left_edge = np.arange(x0,x1,bin_width)
    bin_center = [x+bin_width/2 for x in bin_left_edge]
    fig, ax = plt.subplots() # very important this is not inside the for loop
    for i in range(n_positions):
        figure(1, figsize=(6, 6), dpi=80)
        ax.cla()
        #ax.scatter(pos[i],vel[i], s=1)
        bin_index = np.digitize(pos[i],bin_left_edge)
        n_bins = max(bin_index)
        vel_values = [[] for _ in range(n_bins)]
        k=0
        for j in bin_index:
            vel_values[j-1].append(vel[i][k])
            k+=1
        k=0
        for packet in vel_values:
            if len(packet) == 0:
                vel_values[k] = 0
            elif len(packet) != 1:
                # sorted in magnitude i.e [-5,10,15,-20]
                sorted_packet = sorted(packet, key=abs)
                # [-1] = -20, [-2] = 15
                if (abs(sorted_packet[-1] - sorted_packet[-2]) > 3):
                    packet.remove(sorted_packet[-1])
                    # if len(packet) > 1:
                    #     sorted_packet = sorted(packet, key=abs)
                    #     if (abs(sorted_packet[-1] - sorted_packet[-2]) > 3):
                    #         packet.remove(sorted_packet[-1])
                vel_values[k] = sum(packet)/len(packet)
            else:
                vel_values[k] = sum(packet)/len(packet)
            k+=1
        ax.plot(bin_center, vel_values, '-', linewidth=0.4)



        #ax.plot(pos[i],vel[i], '-o', linewidth=0.4, markersize=1)
        
        ax.set_xlabel("Distance from Saturn (km)")
        ax.set_ylabel("Radial velocity (m/s)")

        plt.subplots_adjust(left = 0, bottom=0)
        # ax.set_xlim([saturn_equatorial_radius,dist_mimas_saturn])
        #ax.set_xlim([x0,x1])
        ax.set_ylim([0,12])

        x2 = dist_mimas_saturn * pow(0.5,2/3)
        ax.set_xlim([x2-300,x2+300])
        plt.vlines(x=x2, ymin=-50, ymax=50, colors='red', ls='--', lw=0.4, label="2:1")
        # xx=[]
        # yy=[]
        # zz=[]
        # aa=[]
        # for j in range(2,10):
        #     k = 1 - 1/j
        #     k = math.pow(k,2/3) * dist_mimas_saturn
        #     xx.append(k)
        # plt.vlines(x=xx, ymin=-5, ymax=5, colors='red', ls='--', lw=0.4, label="l=m")

        # k = math.pow(1/3,2/3) * dist_mimas_saturn
        # yy.append(k)
        # k = math.pow(3/5,2/3) * dist_mimas_saturn
        # yy.append(k)
        # k = math.pow(5/7,2/3) * dist_mimas_saturn
        # yy.append(k)
        # k = math.pow(7/9,2/3) * dist_mimas_saturn
        # yy.append(k)
        # k = math.pow(9/11,2/3) * dist_mimas_saturn
        # yy.append(k)
        # k = math.pow(11/13,2/3) * dist_mimas_saturn
        # yy.append(k)
        # k = math.pow(13/15,2/3) * dist_mimas_saturn
        # yy.append(k)
        # k = math.pow(15/17,2/3) * dist_mimas_saturn
        # yy.append(k)
        # k = math.pow(17/19,2/3) * dist_mimas_saturn
        # yy.append(k)
        # plt.vlines(x=yy, ymin=-5, ymax=5, colors='green', ls='--', lw=1, label="l=m+1")

        # k = math.pow(1/4,2/3) * dist_mimas_saturn
        # zz.append(k)
        # k = math.pow(2/5,2/3) * dist_mimas_saturn
        # zz.append(k)
        # k = math.pow(4/7,2/3) * dist_mimas_saturn
        # zz.append(k)
        # k = math.pow(5/8,2/3) * dist_mimas_saturn
        # zz.append(k)
        # k = math.pow(7/10,2/3) * dist_mimas_saturn
        # zz.append(k)
        # k = math.pow(10/13,2/3) * dist_mimas_saturn
        # zz.append(k)
        # k = math.pow(11/14,2/3) * dist_mimas_saturn
        # zz.append(k)
        # k = math.pow(13/16,2/3) * dist_mimas_saturn
        # zz.append(k)
        # k = math.pow(14/17,2/3) * dist_mimas_saturn
        # zz.append(k)
        # plt.vlines(x=zz, ymin=-5, ymax=5, colors='cyan', ls='--', lw=1, label="l=m+2")

        # k = math.pow(1/5,2/3) * dist_mimas_saturn
        # aa.append(k)
        # k = math.pow(3/7,2/3) * dist_mimas_saturn
        # aa.append(k)
        # k = math.pow(5/9,2/3) * dist_mimas_saturn
        # aa.append(k)
        # k = math.pow(7/11,2/3) * dist_mimas_saturn
        # aa.append(k)
        # k = math.pow(9/13,2/3) * dist_mimas_saturn
        # aa.append(k)
        # k = math.pow(11/15,2/3) * dist_mimas_saturn
        # aa.append(k)
        # k = math.pow(13/17,2/3) * dist_mimas_saturn
        # aa.append(k)
        # k = math.pow(15/19,2/3) * dist_mimas_saturn
        # aa.append(k)
        # k = math.pow(17/21,2/3) * dist_mimas_saturn
        # aa.append(k)
        # plt.vlines(x=aa, ymin=-5, ymax=5, colors='magenta', ls='--', lw=1, label="l=m+3")
        
        plt.legend(loc="upper right")

        fig.set_size_inches(5, 5)
        plt.savefig("E:/Pictures/huygens_v_abs_2_ellipse_parametric/foo{}.jpg".format(str(i).zfill(3)), bbox_inches="tight", dpi=500) #dpi=200 is ~1080p
        print(str(i+1) + "/" + str(n_positions),end="\r")
    print("\n")

#radial_velocity_subcycle()

# def position_density_subcycle():
#     # THIS FUNCTION WAS COMMENTED OUT (MOST LIKELY NOT GOOD)

#     n_positions = int((t*n_orbits)/(timestep*subcycle)) # =815
#     print("Reading txt file and creating pos and vel lists...")
#     with open("E:/Warwick/saturn/positions_particles_huygens.txt") as f:
#         pos = np.empty([n_positions,n_particles])
#         pos1 = np.empty([n_particles])
#         line_number = 0
#         position_number = 0 #[0,815-1]
#         for line in f:
#             line_number += 1
#             r = line.split()
#             R = [float(x) for x in r]
#             pos1[(line_number-1)%n_particles] = R[0]
#             if line_number%n_particles == 0:
#                 pos[position_number] = pos1
#                 position_number += 1
#                 print(str(position_number) + "/" + str(n_positions),end="\r")
    
#     print("Removing stray particles...")
#     particle_number = 0
    
#     #pos list of times rather than particles

#     print("Creating plot and saving as jpg...")
#     x0 = 110000#saturn_equatorial_radius
#     x1 = 124000#dist_mimas_saturn
#     L= x1 - x0
#     x1 = 10 # no. of particles per bin (I want 10 particles/bin)
#     max_pos = 20
#     bin_width = L*x1/n_particles # 100
#     n_bins = int(n_particles/x1)
#     height = n_bins
#     width = n_positions
#     data = numpy.zeros((height, width, 3), dtype=numpy.uint8)
#     colour = [216,174,109]
#     for i in range(n_positions):
#         for val in pos[i]:
#             if ((val < x0) or (val > x1)):
#                 #remove pos if not in range
#                 #this changes size of pos[i], so no longer 10,000 particles
#                 numpy.delete(pos[i], np.where(pos[i] == val))

#         bin_left_edge = np.arange(x0,x1,bin_width)
#         #bin_left_edge = np.arange(110000,124000,bin_width)
#         bin_center = [x+bin_width/2 for x in bin_left_edge]
#         bin_index = np.digitize(pos[i],bin_left_edge) # list of indexes
#         pos_count = [0]*n_bins
#         k=0
#         for j in bin_index:
#             pos_count[j-1] += 1
#             k+=1
#         k=0
#         for x in pos_count:
#             if x > max_pos: #avg should be x1=10 particles/bin
#                 pos_count[k] = max_pos
#             k+=1

#         #now have a 1D list of pos_counts, need to draw to pixel column
#         #print("\n")
#         # x_n = 0
#         # for j in range(height):
#         #     val = pos_count[x_n]
#         #     scale = val/max_pos
#         #     data[j,i] = [x*scale for x in colour]
#         #     x_n+=1
        
#         # print(str(i+1) + "/" + str(n_positions),end="\r")

#         plt.hist2d(x, y, bins =[x_bins, y_bins])

#     image = Image.fromarray(data)
#     image.save("E:/Pictures/saturnringdensity1.png")

#     ax = plt.subplot()
#     im = ax.imshow(data, extent=(0,n_positions-1,x0,x1))
#     ax.set_xlabel("Cycles (10,000*10s)")
#     ax.set_ylabel("Distance from Saturn (km)") 
#     plt.savefig("E:/Pictures/saturnringdensity1graph.png", bbox_inches="tight", dpi=300)
#     # divider = make_axes_locatable(ax)
#     # cax = divider.append_axes("top", size="5%", pad=0.05)

#     # plt.colorbar(im, cax=cax)

#position_density_subcycle()

def position_density_hist2d_subcycle():
    n_positions = int((t*n_orbits)/(timestep*subcycle)) # =815
    print("Reading txt file and creating pos and vel lists...")
    start1 = time.time()
    with open("E:/Warwick/saturn/positions_particles_huygens_mimas_massx1_116382-117382km_dt=100_n_orb=8000_sub=500_ellipse_parametric.txt") as f:
        pos2 = np.empty([n_particles*n_positions], dtype='float64')
        line_number = 0
        #position_number = 0 #[0,815-1]
        for line in f:
            line_number += 1
            # if (line_number > n_particles*n_positions):
            #     continue
            r = line.split()
            R = [float(x) for x in r] # r v
            pos2[line_number-1] = R[0] # r
            #print(str(line_number) + "/" + str(n_particles*n_positions),end="\r")
    end1 = time.time()
    print("Finished in: " + str(end1 - start1) + " seconds")

    print("Number of lines: " + str(line_number))
    print("Size of array: " + str(n_particles*n_positions))

    #print("Removing stray particles...")
    #particle_number = 0
    #pos list of times rather than particles

    x0 = 116382#115500#94400#110000#saturn_equatorial_radius
    x1 = 117382#118000#96400#124000#dist_mimas_saturn
    #L = x1 - x0
    N = 5 # no. of particles per bin (I want 10 particles/bin)
    #bin_width = L*N/n_particles # 125.3
    #bin_left_edge = np.arange(x0,x1,bin_width)
    #t_bins = np.arange(0,n_positions,1)
    n_bins = int(n_particles/N)
    #height = n_bins
    #width = n_positions
    #colour = [216,174,109]
    ax = plt.subplot()
    x = np.empty([n_particles*n_positions], dtype='uint32')

    print("Creating X array...")
    start2 = time.time()
    for i in range(n_positions):
        for j in range(n_particles):
            x[j+i*n_particles] = i
    end2 = time.time()
    print("Finished in: " + str(end2 - start2) + " seconds")
    
    print("Creating plot and saving as jpg...")
    start3 = time.time()
    plt.hist2d(x, pos2, bins=[n_positions, n_bins], range=[[0,n_positions-1],[x0,x1]])#, cmap=plt.cm.copper)
    #plt.hist2d(x, pos2, bins =[n_positions, bin_left_edge])
    #print(str(i+1) + "/" + str(n_positions),end="\r")
    end3 = time.time()
    print("Finished in: " + str(end3 - start3) + " seconds")

    plt.title("Ring particle density over time")
  
    plt.colorbar()
    #plt.xticks(range(0, n_positions,100))
    #plt.xlim([0, n_positions])
    plt.ylim([x0, x1])

    y1 = dist_mimas_saturn * pow(1/2,2/3)
    plt.hlines(y=y1, xmin=0, xmax=n_positions-1, colors='red', ls='--', lw=0.4, label="2:1 Mimas")
    plt.legend(loc="upper right")

    ax.set_xlabel("Time (500*100s)")
    ax.set_ylabel("Distance from Saturn (km)") 
    plt.savefig("E:/Pictures/saturnringdensitygraph_mimas_ellipse_parametric_n_orb=8000_sub=500_zoom_N=5_dpi=2000.png", bbox_inches="tight", dpi=2000)

#position_density_hist2d_subcycle()

def potential_contour_map():
    x1 = 400000
    x0 = 400000

    fig=plt.figure()
    N = 1000
    v=np.linspace(-x0,x1,N)
    x = v
    y = v
    z = np.empty((N,N))

    ycount = 0
    for j in y:
        xcount = 0
        for i in x:
            x1 = dist_mimas_saturn - i
            pot_p = -G*saturn_mass/math.sqrt(i*i + j*j)
            pot_s = -G*mimas_mass/math.sqrt(x1*x1 + j*j)
            pot_i = G*mimas_mass*i/(dist_mimas_saturn*dist_mimas_saturn)
            z[ycount,xcount] = pot_p + pot_s + pot_i
            xcount += 1
        ycount += 1
        print(str(ycount) + "/" + str(N),end="\r")
    
    plt.contourf(x,y,z,1000)
    plt.xlabel("x-axis (km)")
    plt.ylabel("y-axis (km)")                               
    plt.colorbar()

    plt.savefig("E:/Pictures/saturnpotentialmap.png", bbox_inches="tight", dpi=300)

#potential_contour_map()

def moon_strengths():
    mimas_mass = 3.75e19
    enceladus_mass = 1.1e20
    tethys_mass = 6.2e20
    dione_mass = 1.1e21
    rhea_mass = 2.3e21
    titan_mass = 1.35e23
    iapetus_mass = 1.8e21
    pan_mass = 4.95e15
    daphnis_mass = 7.7e13
    atlas_mass = 6.6e15
    prometheus_mass = 1.6e17
    pandora_mass = 1.37e17
    janus_mass = 1.9e18
    epimetheus_mass = 5.3e17

    mimas_dist = 185539
    enceladus_dist = 237948
    tethys_dist = 294619
    dione_dist = 377396
    rhea_dist = 527108
    titan_dist = 1221870
    iapetus_dist = 3560820
    pan_dist = 133584
    daphnis_dist = 136505
    atlas_dist = 137670
    prometheus_dist = 139380
    pandora_dist = 141720
    janus_dist = 151460
    epimetheus_dist = 151410

    masses = [
        mimas_mass,
        enceladus_mass,
        tethys_mass,
        dione_mass,
        rhea_mass,
        iapetus_mass,
        pandora_mass,
        janus_mass,
        epimetheus_mass,
        titan_mass,
        pan_mass,
        daphnis_mass,
        atlas_mass,
        prometheus_mass
    ]

    distances = [
        mimas_dist,
        enceladus_dist,
        tethys_dist,
        dione_dist,
        rhea_dist,
        iapetus_dist,
        pandora_dist,
        janus_dist,
        epimetheus_dist,
        titan_dist,
        pan_dist,
        daphnis_dist,
        atlas_dist,
        prometheus_dist
    ]

    names = [
        "Mimas",
        "Enceladus",
        "Tethys",
        "Dione",
        "Rhea",
        "Iapetus",
        "Pandora",
        "Janus",
        "Epimetheus",
        "Titan",
        "Pan",
        "Daphnis",
        "Atlas",
        "Prometheus"
    ]

    n_moons = len(masses)
    timestep = 1000
    N = 1000
    x0 = 74500
    x1 = 140220
    r = np.linspace(x0,x1,N)
    ax = plt.subplot()
    lines = []
    for i in range(n_moons):
        d = distances[i]
        m = masses[i]
        a = np.empty(N)
        I = 0 # counter 
        for j in r:
            moon_period = 2 * math.pi * math.sqrt((d*d*d) / (G * saturn_mass));
            particle_period = 2 * math.pi * math.sqrt((j*j*j) / (G * saturn_mass));
            old_m_x = d
            old_m_y = 0
            old_p_x = j
            old_p_y = 0
            a_x = 0
            a_y = 0
            for k in range(10000):
                r1 = math.sqrt((old_m_x - old_p_x)*(old_m_x - old_p_x)  +  (old_m_y - old_p_y)*(old_m_y - old_p_y))
                a_x += G*m*(old_m_x - old_p_x)/(r1*r1*r1)
                a_y += G*m*(old_m_y - old_p_y)/(r1*r1*r1)
                old_m_x = d*math.cos(2 * math.pi * k * timestep / moon_period)
                old_m_y = d*math.sin(2 * math.pi * k * timestep / moon_period)
                old_p_x = j*math.cos(2 * math.pi * k * timestep / particle_period)
                old_p_y = j*math.sin(2 * math.pi * k * timestep / particle_period)
        
            a1 = math.sqrt(a_x*a_x + a_y*a_y)
            a[I] = a1
            I += 1
            print(str(i+1) + "/" + str(n_moons) + " " + str(I) + "/" + str(N),end="\r")

        # print(str(i) + "/" + str(n_moons),end="\r")
        line, = plt.plot(r, a, lw=0.4, alpha=1, label="{}".format(names[i]))
        lines.append(line)

    plt.title("Resonance strengths from Saturn's moons")
    #plt.legend(loc="upper right")

    resonances = []
    y1 = dist_mimas_saturn * pow(0.5,2/3)
    l1 = plt.vlines(x=y1, ymin=0, ymax=1e-6, colors='red', ls='--', lw=0.3, label="2:1 Mimas")
    resonances.append(l1)

    y2 = janus_dist * pow(6/7,2/3)
    l2 = plt.vlines(x=y2, ymin=0, ymax=1e-6, colors='blue', ls='--', lw=0.3, label="7:6 Janus-Epimetheus")
    resonances.append(l2)

    y3 = janus_dist * pow(5/6,2/3)
    l3 = plt.vlines(x=y3, ymin=0, ymax=1e-6, colors='green', ls='--', lw=0.3, label="6:5 Janus-Epimetheus")
    resonances.append(l3)

    y4 = janus_dist * pow(4/5,2/3)
    l4 = plt.vlines(x=y4, ymin=0, ymax=1e-6, colors='cyan', ls='--', lw=0.3, label="5:4 Janus-Epimetheus")
    resonances.append(l4)

    y5 = janus_dist * pow(3/4,2/3)
    l5 = plt.vlines(x=y5, ymin=0, ymax=1e-6, colors='magenta', ls='--', lw=0.3, label="4:3 Janus-Epimetheus")
    resonances.append(l5)

    y6 = janus_dist * pow(2/3,2/3)
    l6 = plt.vlines(x=y6, ymin=0, ymax=1e-6, colors='orange', ls='--', lw=0.3, label="3:2 Janus-Epimetheus")
    resonances.append(l6)

    y7 = janus_dist * pow(1/2,2/3)
    l7 = plt.vlines(x=y7, ymin=0, ymax=1e-6, colors='darkgreen', ls='--', lw=0.3, label="2:1 Janus-Epimetheus")
    resonances.append(l7)
    

    leg1 = plt.legend(handles=lines, title="Moons", loc="upper left", fontsize="small", fancybox=True, shadow=True)
    plt.gca().add_artist(leg1)
    plt.legend(handles=resonances, title="ILR", loc="upper right", fontsize="x-small")

    ax.set_xlabel("Distance from Saturn (km)")
    ax.set_ylabel("Magnitude total acceleration (km/s²)") 

    plt.xlim([x0, x1])
    plt.ylim([0, 1e-6])
    plt.savefig("E:/Pictures/saturns_moons_strengths4.png", bbox_inches="tight", dpi=1000)

#moon_strengths()

def moon_strengths_1line():
    mimas_mass = 3.75e19
    enceladus_mass = 1.1e20
    tethys_mass = 6.2e20
    dione_mass = 1.1e21
    rhea_mass = 2.3e21
    titan_mass = 1.35e23
    iapetus_mass = 1.8e21
    pan_mass = 4.95e15
    daphnis_mass = 7.7e13
    atlas_mass = 6.6e15
    prometheus_mass = 1.6e17
    pandora_mass = 1.37e17
    janus_mass = 1.9e18
    epimetheus_mass = 5.3e17

    mimas_dist = 185539
    enceladus_dist = 237948
    tethys_dist = 294619
    dione_dist = 377396
    rhea_dist = 527108
    titan_dist = 1221870
    iapetus_dist = 3560820
    pan_dist = 133584
    daphnis_dist = 136505
    atlas_dist = 137670
    prometheus_dist = 139380
    pandora_dist = 141720
    janus_dist = 151460
    epimetheus_dist = 151410

    masses = [
        mimas_mass,
        enceladus_mass,
        tethys_mass,
        dione_mass,
        rhea_mass,
        iapetus_mass,
        pandora_mass,
        janus_mass,
        epimetheus_mass
    ]
    #     titan_mass,
    #     pan_mass,
    #     daphnis_mass,
    #     atlas_mass,
    #     prometheus_mass
    # ]

    distances = [
        mimas_dist,
        enceladus_dist,
        tethys_dist,
        dione_dist,
        rhea_dist,
        iapetus_dist,
        pandora_dist,
        janus_dist,
        epimetheus_dist
    ]
    #     titan_dist,
    #     pan_dist,
    #     daphnis_dist,
    #     atlas_dist,
    #     prometheus_dist
    # ]

    names = [
        "Mimas",
        "Enceladus",
        "Tethys",
        "Dione",
        "Rhea",
        "Iapetus",
        "Pandora",
        "Janus",
        "Epimetheus"
    ]
    #     "Titan",
    #     "Pan",
    #     "Daphnis",
    #     "Atlas",
    #     "Prometheus"
    # ]

    n_moons = len(masses)
    timestep = 1000
    N = 1000
    x0 = 74500
    x1 = 140220
    r = np.linspace(x0,x1,N)
    ax = plt.subplot()
    lines = []
    I = 0
    a = np.empty(N)
    for j in r:
        particle_period = 2 * math.pi * math.sqrt((j*j*j) / (G * saturn_mass));
        a_x = 0
        a_y = 0
        for k in range(100000):
            old_p_x = j*math.cos(2 * math.pi * k * timestep / particle_period)
            old_p_y = j*math.sin(2 * math.pi * k * timestep / particle_period)
            for i in range(n_moons):
                d = distances[i]
                m = masses[i]
                moon_period = 2 * math.pi * math.sqrt((d*d*d) / (G * saturn_mass));
                old_m_x = d*math.cos(2 * math.pi * k * timestep / moon_period)
                old_m_y = d*math.sin(2 * math.pi * k * timestep / moon_period)
                r1 = math.sqrt((old_m_x - old_p_x)*(old_m_x - old_p_x)  +  (old_m_y - old_p_y)*(old_m_y - old_p_y))
                a_x += G*m*(old_m_x - old_p_x)/(r1*r1*r1)
                a_y += G*m*(old_m_y - old_p_y)/(r1*r1*r1)
        
        a1 = math.sqrt(a_x*a_x + a_y*a_y)
        a[I] = a1
        I += 1
        print(str(I) + "/" + str(N), end="\r")

    plt.plot(r, a, color="black", lw=0.4, alpha=1)
        

    plt.title("Resonance locations and strengths")
    #plt.legend(loc="upper right")

    resonances = []
    y1 = dist_mimas_saturn * pow(0.5,2/3)
    l1 = plt.vlines(x=y1, ymin=0, ymax=1e-6, colors='red', ls='--', lw=0.3, label="2:1 Mimas")
    resonances.append(l1)

    y2 = janus_dist * pow(6/7,2/3)
    l2 = plt.vlines(x=y2, ymin=0, ymax=1e-6, colors='blue', ls='--', lw=0.3, label="7:6 Janus-Epimetheus")
    resonances.append(l2)

    y3 = janus_dist * pow(5/6,2/3)
    l3 = plt.vlines(x=y3, ymin=0, ymax=1e-6, colors='green', ls='--', lw=0.3, label="6:5 Janus-Epimetheus")
    resonances.append(l3)

    y4 = janus_dist * pow(4/5,2/3)
    l4 = plt.vlines(x=y4, ymin=0, ymax=1e-6, colors='cyan', ls='--', lw=0.3, label="5:4 Janus-Epimetheus")
    resonances.append(l4)

    y5 = janus_dist * pow(3/4,2/3)
    l5 = plt.vlines(x=y5, ymin=0, ymax=1e-6, colors='magenta', ls='--', lw=0.3, label="4:3 Janus-Epimetheus")
    resonances.append(l5)

    y6 = janus_dist * pow(2/3,2/3)
    l6 = plt.vlines(x=y6, ymin=0, ymax=1e-6, colors='orange', ls='--', lw=0.3, label="3:2 Janus-Epimetheus")
    resonances.append(l6)

    y7 = janus_dist * pow(1/2,2/3)
    l7 = plt.vlines(x=y7, ymin=0, ymax=1e-6, colors='darkgreen', ls='--', lw=0.3, label="2:1 Janus-Epimetheus")
    resonances.append(l7)
    
    plt.legend(handles=resonances, title="ILR", loc="upper left", fontsize="x-small")

    ax.set_xlabel("Distance from Saturn (km)")
    ax.set_ylabel("Magnitude total acceleration (km/s²)") 

    plt.xlim([x0, x1])
    plt.ylim([0, 1e-6])
    plt.savefig("E:/Pictures/resonance_strengths_dt=1000_k=100000.png", bbox_inches="tight", dpi=1000)

#moon_strengths_1line()

def img_scatter_graph_J_E():
    n_positions = int((t*n_orbits)/(timestep*subcycle)) # =815
    print("Reading text file -> 2d list...")
    with open("E:/Warwick/saturn/epimetheus_x_y_z_dt=100_sub=10000_2.txt") as f:
        pos_E = np.empty([n_positions,3])
        line_number = 0
        position_number = 0 #[0,815-1]
        for line in f:
            line_number += 1
            r = line.split()
            R = [float(x) for x in r]
            pos_E[(line_number-1)] = R
            print(str(line_number) + "/" + str(n_positions),end="\r")

    #print("\n")
    position_number = 0
    print("\nRemoving z axis and stray data points...")

    with open("E:/Warwick/saturn/janus_x_y_z_dt=100_sub=10000_2.txt") as f:
        pos_J = np.empty([n_positions,3])
        line_number = 0
        position_number = 0 #[0,815-1]
        for line in f:
            line_number += 1
            r = line.split()
            R = [float(x) for x in r]
            pos_J[(line_number-1)] = R
            print(str(line_number) + "/" + str(n_positions),end="\r")

    #print("\n")
    position_number = 0
    print("\nRemoving z axis and stray data points...")
    for pos_list in pos_E:
        position_number = 0
        pos_list = pos_list[:-1]
        position_number+=1
        print(str(position_number) + "/" + str(n_positions),end="\r")
    for pos_list in pos_J:
        position_number = 0
        pos_list = pos_list[:-1]
        position_number+=1
        print(str(position_number) + "/" + str(n_positions),end="\r")

    #N = 1 #[0,814] because len = 815
    
    #print("\n")
    x_J = []
    y_J = []
    x_E = []
    y_E = []

    print("\nCreating list of x and y values for scatter plot...")
    for i in range(n_positions):
        x_J.append(pos_J[i][0])
        y_J.append(pos_J[i][1])
        x_E.append(pos_E[i][0])
        y_E.append(pos_E[i][1])
        print(str(i+1) + "/" + str(n_positions),end="\r")

    
    print("\nCreating scatter plot and saving as jpg...")
    fig, ax = plt.subplots() # very important this is not inside the for loop
    X_J_list = []
    Y_J_list = []
    X_E_list = []
    Y_E_list = []
    for i in range(n_positions):
        figure(1, figsize=(6, 6), dpi=80)
        ax.cla()
        r = 151452.7#151500#151458.5
        t1 = 2*math.pi*math.sqrt(r*r*r/(G*saturn_mass))
        k = - (timestep * subcycle / t1 - 16)#1 - timestep * subcycle / t1
        c = math.cos((i)*2*math.pi*k)
        s = math.sin((i)*2*math.pi*k)

        X_J = x_J[i]*c - y_J[i]*s
        Y_J = x_J[i]*s + y_J[i]*c
        X_E = x_E[i]*c - y_E[i]*s
        Y_E = x_E[i]*s + y_E[i]*c
        X_J_list.append(X_J)
        Y_J_list.append(Y_J)
        X_E_list.append(X_E)
        Y_E_list.append(Y_E)
        if (i+1 > 300):
            X_J_list.remove(X_J_list[0])
            Y_J_list.remove(Y_J_list[0])
            X_E_list.remove(X_E_list[0])
            Y_E_list.remove(Y_E_list[0])

        ax.plot(X_J_list, Y_J_list, lw=0.3, c="red")
        ax.plot(X_E_list, Y_E_list, lw=0.3, c="blue")
        ax.scatter(X_J, Y_J, s=20*2.35, c="red", label="Janus")
        ax.scatter(X_E, Y_E, s=20, c="blue", label="Epimetheus")

        #ax.scatter(x_J[i], y_J[i], s=60, c="red", label="Janus")
        #ax.scatter(x_E[i], y_E[i], s=20, c="blue", label="Epimetheus")

        r = dist_mimas_saturn
        b = timestep/t
        #ax.scatter(r*math.cos(2*math.pi*i*subcycle*b),r*math.sin(2*math.pi*i*subcycle*b), s=30, c="red", label="Mimas")

        time = i * 100 * 10000 / (60*60*24*365)
        dx = math.sqrt((X_J - X_E)*(X_J - X_E) + (Y_J - Y_E)*(Y_J - Y_E))
        r_J = math.sqrt(X_J*X_J + Y_J*Y_J)
        r_E = math.sqrt(X_E*X_E + Y_E*Y_E)
        at = AnchoredText("Time = {:.2f}yrs\nr_J = {:,.2f}km\nr_E = {:,.2f}km\nΔ = {:,.2f}km".format(time, r_J, r_E, dx), prop=dict(size=10), frameon=True, loc='upper left')
        at.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
        ax.add_artist(at)

        w = 360 * 60 * 60 / t1
        at2 = AnchoredText("Rotating frame, ω={:.2f}°/hr\nMoon sizes are exaggerated".format(w), prop=dict(size=10), frameon=False, loc='lower right')
        at2.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
        ax.add_artist(at2)

        x = 0
        y = 0

        circle = plt.Circle((x, y), radius=saturn_equatorial_radius, color="Yellow")

        ax.add_patch(circle)

        label = ax.annotate("SATURN", xy=(x, y), fontsize=20, ha="center", va="center")

        plt.legend(loc="upper right")


        ax.set_xlabel("x axis (km)")
        ax.set_ylabel("y axis (km)")

        plt.subplots_adjust(left = 0, bottom=0)
        ax.set_xlim([-200000,200000])
        ax.set_ylim([-200000,200000])

        fig.set_size_inches(7, 7)
        plt.savefig("E:/Pictures/Janus-Epimetheus9/foo{}.jpg".format(str(i).zfill(3)), bbox_inches="tight", dpi=500)
        print(str(i+1) + "/" + str(n_positions),end="\r")
    print("\n")

#img_scatter_graph_J_E()

def lindblad_plot():
    r = dist_mimas_saturn
    fig1 = plt.figure()
    ax1 = fig1.add_subplot(111, aspect='equal')
    n = 200
    for i in range(n):
        i+=1
        ILR = pow(((i+1)-1)/(i+1),2/3)*r
        OLR = pow((i+1)/i,2/3)*r

        ax1.add_patch(patches.Circle((0, 0), radius=ILR, color='r', linewidth=0.2, fill=False))
        ax1.add_patch(patches.Circle((0, 0), radius=OLR, color='b', linewidth=0.2, fill=False))
        print(str(i) + "/" + str(n),end="\r")
    
    ax1.add_patch(patches.Circle((0, 0), radius=saturn_equatorial_radius, color='y', fill=True))
    ax1.add_patch(patches.Circle((0, 0), radius=r, color='k', linewidth=0.2, fill=False))
    ax1.autoscale_view()
    x1 = 350000
    x0 = 40000
    y1 = 50000
    y0 = -50000
    ax1.set_xlim([x0,x1])
    ax1.set_ylim([y0,y1])
    ratio = (y1-y0)/(x1-x0)
    y = 7 * ratio
    x = 7
    red_patch = mpatches.Patch(color='red', label='ILR')
    blue_patch = mpatches.Patch(color='blue', label='OLR')
    black_patch = mpatches.Patch(color='black', label='CC')
    #ax1.legend(handles=[red_patch,blue_patch,black_patch], handlelength=1, bbox_to_anchor=(0,1.02,1,0.2),mode="expand", loc="lower left", ncol=3, prop=dict(size=7))
    ax1.legend(handles=[red_patch,black_patch,blue_patch], handlelength=1, bbox_to_anchor=(0.5,1.02), loc="lower center", ncol=3, prop=dict(size=7))
    fig1.set_size_inches(x, y)
    #plt.axis('off')
    ax1.get_xaxis().set_visible(False)
    ax1.get_yaxis().set_visible(False)
    plt.savefig("E:/Pictures/lindblad4.png", bbox_inches="tight", dpi=800)

#lindblad_plot()