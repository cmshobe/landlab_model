# -*- coding: utf-8 -*-
"""
Created on Mon Mar 14 09:29:27 2016

@author: Charlie Shobe

Landlab script 

Components:
1) fastscape stream power erosion 
2) linear hillslope diffusion 
3) stochastic storm distribution 
4) option to vary k: top 10 meters of rock are hard, everything under is soft.

INPUT: parameters .txt file (example on GitHub)
OUTPUT: plot of final topography (initial topo is flat plain)

BOUNDARIES: Open on east and west, closed on north and south

Uplift condition: steady constant uplift at rate specified in params file
"""

import numpy as np
import matplotlib.pyplot as plt

from landlab import ModelParameterDictionary
from landlab import RasterModelGrid
from landlab.components import FlowRouter as Flow
#from landlab.components.stream_power.stream_power import StreamPowerEroder
from landlab.components import FastscapeEroder as Fsc
from landlab.components import LinearDiffuser as Diff
from landlab.components import PrecipitationDistribution
#from landlab.plot import channel_profile as prf

###########################INITIALIZE

#get the stuff needed to build the grid AND RUN THE MODEL:
input_file = './shobe_modClass_params.txt'
inputs = ModelParameterDictionary(input_file)
nrows = inputs.read_int('nrows')
ncols = inputs.read_int('ncols')
dx = inputs.read_float('dx')
dx = inputs.read_float('dx')
leftmost_elev = inputs.read_float('leftmost_elevation')
initial_slope = inputs.read_float('initial_slope')
uplift_rate = inputs.read_float('up_rate')
runtime = inputs.read_float('total_time')
dt = inputs.read_float('dt')
nt = int(runtime//dt)
t_plot = inputs.read_float('time_plot')
uplift_per_step = uplift_rate * dt

#Variable erodibility
hard_layer_on_or_off = inputs.read_int('layered')
k_erodible=inputs.read_float('k_erodible') #higher k for weaker rock
k_unerodible=inputs.read_float('k_unerodible') #lower k for stronger rock
hard_layer_thickness = inputs.read_float('hard_layer_thickness')

#build the grid:
mg = RasterModelGrid(nrows, ncols, dx)

#create elevation field for the grid
mg.add_zeros('node','topographic__elevation')
z = mg.zeros('node') + leftmost_elev
z += initial_slope*np.amax(mg.node_y) - initial_slope*mg.node_y

#put these values plus roughness into that field
mg['node'][ 'topographic__elevation'] = z + np.random.rand(len(z))/100000.

#initial k conditions for grid if using layers: all on hard rock
k = mg.add_zeros('node','K_values')
if hard_layer_on_or_off == 1:
    k[:] = k_unerodible #for hard rock
    mg['node']['K_values'] = k
elif hard_layer_on_or_off == 0:
    k[:] = k_erodible #soft rock
else:
    print 'WARNING: MUST SELECT 0 OR 1 IN LAYERED PARAM'

#set up its boundary conditions (left, top, right, bottom is inactive)
mg.set_closed_boundaries_at_grid_edges(False, True, False, True)

# Display initialization message
print('Running ...') 

#instantiate the components:
pr = PrecipitationDistribution(input_file)
fr = Flow(mg)
sp = Fsc(mg, input_file)
hd = Diff(mg, input_file)

####################RUN
track_uplift = 0 #track cumulative uplift to know top of hard layer
last_trunc = runtime
for (interval_duration, rainfall_rate) in pr.yield_storm_interstorm_duration_intensity():
    if rainfall_rate != 0.:
        # note diffusion also only happens when it's raining...
        _ = fr.route_flow()
        _ = sp.erode(mg, interval_duration, K_if_used='K_values')
        _ = hd.diffuse(interval_duration)
    track_uplift += uplift_rate * interval_duration #top of beginning surface
    mg.at_node['topographic__elevation'][mg.core_nodes] += uplift_rate * interval_duration
    this_trunc = pr.elapsed_time // t_plot
    if this_trunc != last_trunc: # time to plot a new profile!
        print ('Time %d' % (t_plot * this_trunc))
        last_trunc = this_trunc
    else:
        pass
    
    #check where hard rocks and soft rocks are, change k to reflect this
    if hard_layer_on_or_off == 1: #if using layers
        hard_layer = np.where(mg.at_node['topographic__elevation'] >= track_uplift - hard_layer_thickness)
        soft_layer = np.where(mg.at_node['topographic__elevation'] < track_uplift - hard_layer_thickness)
        k[hard_layer] = k_unerodible
        k[soft_layer] = k_erodible
    elif hard_layer_on_or_off == 0: #if not
        k[:] = k_erodible #soft rock
    else:
        print 'WARNING: MUST SELECT 0 OR 1 IN LAYERED PARAM'
    mg['node']['K_values'] = k

############FINALIZE

#Get resulting topography, turn into a raster
elev = mg['node']['topographic__elevation']
elev_r = mg.node_vector_to_raster(elev)

# Plot topography
plt.figure(1)
im = plt.imshow(elev_r, cmap=plt.cm.jet, extent=[0,1000,0,1000], origin='lower')  # display a colored image
plt.colorbar(im)
plt.title('Topography')

plt.show()

print('Done!')