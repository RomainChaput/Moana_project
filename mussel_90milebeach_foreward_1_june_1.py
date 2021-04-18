#!/usr/bin/env python

import os
import sys
import numpy as np
from datetime import datetime, timedelta
from opendrift.readers import reader_global_landmask
from opendrift.readers import reader_ROMS_native_MOANA
from opendrift.models.bivalvelarvae import BivalveLarvae

###############################
# PARAMETERS TO CHANGE IN BETWEEN RUNS
###############################

runtime = [datetime(2017,6,1), datetime(2017,6,8)]
finish_time = (datetime(2017,6,8,12) + timedelta(days = 36.0))
namefile = '/nesi/nobackup/vuw03295/ninety_miles_beach_foreward_1_june_1_2017_005deg_Moana_mask.nc'

###############################
# MODEL SELECTION
###############################

o = BivalveLarvae(loglevel=100)  # Set loglevel to 0 for debug information
o.max_speed = 5.0#

###############################
# READERS
###############################

reader_moana_v19_1 = reader_ROMS_native_MOANA.Reader('/nesi/nobackup/mocean02574/NZB_N50/nz5km_his_201706.nc')
reader_moana_v19_2 = reader_ROMS_native_MOANA.Reader('/nesi/nobackup/mocean02574/NZB_N50/nz5km_his_201707.nc')
reader_moana_v19_3 = reader_ROMS_native_MOANA.Reader('/nesi/nobackup/mocean02574/NZB_N50/nz5km_his_201708.nc')

o.add_reader([reader_moana_v19_1, reader_moana_v19_2, reader_moana_v19_3])
o.set_config('general:use_auto_landmask', False) # dynamical landmask if true

###############################
# PARTICLE SEEDING
###############################

# Define the starting position of the particles within the buffer of intertidal rocky shore on the West Coast of the North Island (NZ)
nb_parts = 250 # number of particles released per point at every timestep
points = np.loadtxt('/nesi/project/vuw03295/Regional_projects/Ninety_miles_beach/First_foreward_input_files/Release_points_rocky_shore_005deg_Moana_all_poly.xyz', delimiter='\t', dtype=str)
plon = points[:,0].astype(np.float)
plat = points[:,1].astype(np.float)
tot_parts = nb_parts * len(plon) # total number of particles released
plon = np.tile(plon, nb_parts)
plat = np.tile(plat, nb_parts)

## Define the release time. author: Calvin Quigley (19/04/2021)
#def create_seed_times(start, end, delta):
#  """
#  crate times at given interval to seed particles
#  """
#  out = []
#  start_t = start
#  end_t = datetime.strptime(str(end), "%Y-%m-%d %H:%M:%S")
#  while start_t < end:
#    out.append(start_t) 
#    start_t += delta
#  return out

#times = create_seed_times(runtime[0], 
#                          runtime[1], timedelta(days = 1))

# Define release depth
z = np.random.uniform(-20,-1,size=tot_parts) # generate random depth
#z = 'seafloor+1'

# Seed particles
#for i in range(len(times)):
#	o.seed_elements(plon, plat, number=tot_parts, z=z, time = times[i], terminal_velocity=0.0025)
o.seed_elements(plon, plat, number=tot_parts, z=z, time = runtime, terminal_velocity=0.0025)

# seed
o.set_config('seed:ocean_only',True) # keep only particles from the "frame" that are on the ocean

###############################
# PHYSICS
###############################

o.set_config('environment:fallback:x_wind', 0.0)
o.set_config('environment:fallback:y_wind', 0.0)
o.set_config('environment:fallback:x_sea_water_velocity', 0.0)
o.set_config('environment:fallback:y_sea_water_velocity', 0.0)
o.set_config('environment:fallback:sea_floor_depth_below_sea_level', 100000.0)

# drift
o.set_config('drift:advection_scheme','runge-kutta4') # or 'runge-kutta'
o.set_config('drift:current_uncertainty', 0.0 ) # note current_uncertainty can be used to replicate an horizontal diffusion spd_uncertain = sqrt(Kxy*2/dt)  
# o.set_config('drift:lift_to_seafloor',True) # Already included in bivalve larvae model

# horizontal and vertical diffusion
# Need to be turned off for backtracking

Kxy = 0.1176 #0.1176 # m2/s-1
o.set_config('drift:horizontal_diffusivity',Kxy) # using new config rather than current uncertainty
o.set_config('drift:vertical_mixing', True) 
Kz = 0.01 # m2/s-1
o.set_config('environment:fallback:ocean_vertical_diffusivity', Kz) # specify constant ocean_vertical_diffusivity in m2.s-1
o.set_config('vertical_mixing:diffusivitymodel', 'constant') # constant >> use fall back values o.set_config('environment:fallback:ocean_vertical_diffusivity'] = Kz 
o.set_config('vertical_mixing:timestep', 90.0)  # if some ocean_vertical_diffusivity!=0, turbulentmixing:timestep should be << 900 seconds
o.set_config('drift:vertical_advection', True)


##############################
# CONFIG
# Adjusting some configuration
o.set_config('general:coastline_action', 'previous') # Settlement at the coast only possible within specific polygons: post-processing

##############################
#if False: # Only with plankton model
#    # plankton-specific config for developed module based on 
o.set_config('drift:max_age_seconds', 35*24*3600) # PLD spats (algae + mussels) = 20 days; PLD mussels = 5 weeks
o.set_config('drift:min_settlement_age_seconds', 35*24*3600.0)  #'float(min=0.0, max=100.0, default=0.0)', comment='Minimum age before beaching can occur, in seconds') => 3 weeks normally (21 days), but I will deal with that in post-processing
###############################

o.list_config()

###############################
# RUN 
###############################
# Running model (until end of driver data)
o.run(stop_on_error=False,
      time_step=900.0, 
      end_time = finish_time,
      outfile= namefile,
      time_step_output = 3600.0)


###############################
# Plot
###############################

#print(o)
#plt = o.plot(fast=True)
#plt.savefig('/nesi/nobackup/vuw03295/ninety_miles_beach_foreward_1_july_2_2017_005deg_opendrift_mask.png')
#anim = o.animation(fast=True, color='z', background=['x_sea_water_velocity', 'y_sea_water_velocity'])
#anim.save('/nesi/nobackup/vuw03295/ninety_miles_beach_foreward_1_july_2_2017_005deg_opendrift_mask.gif', writer='imagemagick', fps=30)
#o.plot_vertical_distribution()
