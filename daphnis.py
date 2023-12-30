import math
import numpy
import numpy as np
from PIL import Image
import time

import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np
import matplotlib.patches as mpatches
from matplotlib.offsetbox import AnchoredText

G = 6.67e-20  # big G with km instead of m (divide by 10^9).

saturn_equatorial_radius = 60268  # km
saturn_polar_radius = 60268  # km 54364
saturn_mass = 5.683e26  # kg

daphnis_radius = 136505 # km
N = 800# 320 # 800
M = 100 #80 # 100
L = 3200# 1280
p0 = 162
q = 95
height = M
width = N
dx = L/N
n0 = 1000#17 # 1000/10
n = 1 + n0 # ^ 1000, every 10th => 100 frames + 1 = 101 since I included initial frame.


# Create a 1024x1024x3 array of 8 bit unsigned integers
red = [255, 0, 0]
green = [0, 255, 0]
blue = [0, 0, 255]
gold = [255, 215, 0]
purple = [255, 0, 255]
brown = [153, 76, 0]
black = [0, 0, 0]
white = [255, 255, 255]

# reset background
def reset_background():
    for j in range(height):
        for i in range(width):
            #num = 10*j + i + 1
            if (i+j) % 2:
                data[j, i] = black
            else:
                data[j, i] = black

# reset_background()

def read_txt_file(n, N, M):
    with open("E:/Warwick/saturn/daphnis_e10.txt") as f:
        pos = np.empty([n,M,N])
        pos1 = np.empty([M,N])
        line_number = 0 #[0,M]
        position_number = 0 #[0,n]
        for line in f:
            line_number += 1
            r = line.split()
            R = [float(x) for x in r]
            pos1[(line_number-1)%M] = R # don't need to add (M-1)- to index.
            if line_number%M == 0:
                pos[position_number] = pos1
                position_number += 1
                print(str(position_number) + "/" + str(n),end="\r")
    return pos


# scale
def scale_img(scale_factor, data):
    data2 = numpy.zeros((height*scale_factor, width*scale_factor, 3), dtype=numpy.uint8)
    for j in range(height):
        for i in range(width):
            for k in range(scale_factor):
                for l in range(scale_factor):
                    data2[j*scale_factor+k, i*scale_factor+l] = data[j, i]

    return data2

def draw_moon(N,M):
    red = [255, 0, 0]
    j1 = int(N/2)
    i1 = int(M/2)
    data[i1,j1] = red
    if (N % 2 == 0):
        data[i1-1,j1] = red
        if (M % 2 == 0):
            data[i1,j1-1] = red
            data[i1-1,j1-1] = red
    if (M % 2 == 0):
        data[i1,j1-1] = red

def plot_graph(data, x_n):
    ax = plt.subplot()
    ax.cla()
    im = ax.imshow(data, extent=(-1*N*dx/2,N*dx/2,-1*M*dx/2 + daphnis_radius,M*dx/2 + daphnis_radius))
    ax.set_xlabel("Distance (km)")
    ax.set_ylabel("Distance from Saturn (km)")  
    
    red_patch = mpatches.Patch(color='red', label='Daphnis')
    ax.legend(handles=[red_patch], handlelength=0.7, bbox_to_anchor=(1, 1), loc='lower right', prop=dict(size=7))

    time = x_n * 100 * 100 / (60*60)
    daphnis_time = 2*math.pi*math.sqrt(pow(daphnis_radius,3)/(G*saturn_mass)) / (60*60)
    orbits = time / daphnis_time

    at1 = AnchoredText("Time = {:,.2f}hrs".format(time), prop=dict(size=7), frameon=False, bbox_to_anchor=(0., 1.), loc='lower left', bbox_transform=ax.transAxes)
    at1.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
    ax.add_artist(at1)
    at2 = AnchoredText("Daphnis orbits = {:.2f}".format(orbits), prop=dict(size=7), frameon=False, bbox_to_anchor=(0.25, 1.), loc='lower left', bbox_transform=ax.transAxes)
    at2.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
    ax.add_artist(at2)
    at3 = AnchoredText("D = visc = 1e-12km²/s", prop=dict(size=7), frameon=False, bbox_to_anchor=(0.55, 1.), loc='lower left', bbox_transform=ax.transAxes)
    at3.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
    ax.add_artist(at3)

    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()

    # divider = make_axes_locatable(ax)
    # cax = divider.append_axes('bottom', size='15%', pad=0.6)
    # fig.colorbar(im, cax=cax, orientation='horizontal')
    plt.savefig("E:/Pictures/daphnis_5/foo{}.png".format(str(x_n).zfill(3)), bbox_inches="tight", dpi=300)


data = numpy.zeros((height, width, 3), dtype=numpy.uint8)
#data = numpy.zeros((height, width), dtype=numpy.uint8)
pos = read_txt_file(n, N, M);

max = 90
min = 100
# for x in pos:
#     for j in range(height):
#         for i in range(width):
#             if x[j,i] > max:
#                 max = x[j,i]
#             if x[j,i] < min:
#                 min = x[j,i]

# p0 = max
# q = min
# print("max=" + str(p0) + ", min=" + str(q))
p0 = 300
q=0


#print("\n")
x_n = 0
for x in pos:
    for j in range(height):
        for i in range(width):
            val = round(255*(x[j,i]-q)/(p0-q))
            if val < 0:
                val = 0;
            if val > 255:
                val = 255;
            data[j,i] = [val,val,val]
    
    draw_moon(N,M)
    
    #data2 = scale_img(10,data)
    #image = Image.fromarray(data)
    #image.save("E:/Pictures/daphnis_1/foo{}.png".format(str(x_n).zfill(3)))
    plot_graph(data, x_n)
    x_n+=1
    print(str(x_n) + "/" + str(n),end="\r")
    