#!/usr/bin/env python

import os
import sys
import numpy as np
from datetime import datetime, timedelta
from opendrift.readers import reader_global_landmask
from opendrift.readers import reader_ROMS_native_MOANA
from opendrift.models.oceandrift import OceanDrift
#from opendrift.models.pelagicplankton_moana import PelagicPlanktonDrift
#from opendrift.models.sedimentdrift import SedimentDrift



###############################
# MODEL SELECTION
###############################

o = OceanDrift(loglevel=100)
#o = SedimentDrift(loglevel=100)  # 0 for debug output
#o = PelagicPlanktonDrift(loglevel=50)  # Set loglevel to 0 for debug information
o.max_speed = 3.0#

###############################
# READERS
###############################

thredds_path_1 = 'http://thredds.moanaproject.org:8080/thredds/dodsC/moana/ocean/NZB/v1.9/raw_3D/nz5km_his_201710.nc?ntimes,dt,hc,grid,s_rho[5:1:12],Cs_r[5:1:12],h[0:1:466][0:1:396],lon_rho[0:1:466][0:1:396],lat_rho[0:1:466][0:1:396],lon_u[0:1:466][0:1:395],lat_u[0:1:466][0:1:395],lon_v[0:1:465][0:1:396],lat_v[0:1:465][0:1:396],lon_psi[0:1:465][0:1:395],angle[0:1:466][0:1:396],mask_rho[0:1:466][0:1:396],mask_u[0:1:466][0:1:395],mask_v[0:1:465][0:1:396],ocean_time[0:1:248],z_rho[0:1:248][5:1:12][0:1:466][0:1:396],u_eastward[0:1:248][5:1:12][0:1:466][0:1:396],v_northward[0:1:248][5:1:12][0:1:466][0:1:396]' # Limit to selected depths (raw data: [0:1:39])#'http://thredds.moanaproject.org:8080/thredds/dodsC/moana/ocean/NZB/v1.9/raw_3D/nz5km_his_201710.nc'
reader_moana_v19_1 = reader_ROMS_native_MOANA.Reader([thredds_path_1]) #load data for that year
reader_moana_v19_1.multiprocessing_fail = True # bypass the use of multi core for coordinates conversion and seems to make the model run much faster.
thredds_path_2 = 'http://thredds.moanaproject.org:8080/thredds/dodsC/moana/ocean/NZB/v1.9/raw_3D/nz5km_his_201711.nc?ntimes,dt,hc,grid,s_rho[5:1:12],Cs_r[5:1:12],h[0:1:466][0:1:396],lon_rho[0:1:466][0:1:396],lat_rho[0:1:466][0:1:396],lon_u[0:1:466][0:1:395],lat_u[0:1:466][0:1:395],lon_v[0:1:465][0:1:396],lat_v[0:1:465][0:1:396],lon_psi[0:1:465][0:1:395],angle[0:1:466][0:1:396],mask_rho[0:1:466][0:1:396],mask_u[0:1:466][0:1:395],mask_v[0:1:465][0:1:396],ocean_time[0:1:240],z_rho[0:1:240][5:1:12][0:1:466][0:1:396],u_eastward[0:1:240][5:1:12][0:1:466][0:1:396],v_northward[0:1:240][5:1:12][0:1:466][0:1:396]' # to finish the run in the previous month
reader_moana_v19_2 = reader_ROMS_native_MOANA.Reader([thredds_path_2]) #
reader_moana_v19_2.multiprocessing_fail = True # bypass the use of multi core for coordinates conversion and seems to make the model run much faster.

# Making customised landmask - not required here, using ROMS landmask 
#reader_landmask = reader_global_landmask.Reader(
#                   llcrnrlon=171.0, llcrnrlat=184.5,
#                   urcrnrlon=-42.0, urcrnrlat=-32.0)# max is 185deg

# use native landmask of ROMS files
o.add_reader([reader_moana_v19_1, reader_moana_v19_2])
#o.add_reader([reader_moana_v19_1])
o.set_config('general:use_auto_landmask', True) # dynamical landmask

###############################
# PARTICLE SEEDING
###############################

# Define the starting position of the particles on the Ninety miles beach: 200 random points along the beach, spreading 200m offshore.
nb_parts = 1 # number of particles released per point during the release interval
points = np.loadtxt('Release_points_Ninety_miles_beach.xyz', delimiter='\t', dtype=str)
plon = points[:,0].astype(np.float)
plat = points[:,1].astype(np.float)
tot_parts = nb_parts * len(plon) # total number of particles released
plon = np.tile(plon, nb_parts)
plat = np.tile(plat, nb_parts)

#lon0 = 177.3  
#lat0 = -37.855
#radius_in_m = 1000
#z = np.random.uniform(-30,-10,size=tot_parts) # generate random depth
z = 'seafloor+1'

# continuous release from tstart_release to tend_release at random depth
o.seed_elements(plon, plat, 
                      number=tot_parts, 
                      z=z,
                      time = [datetime(2017,11,2), datetime(2017,11,1)])


###############################
# PHYSICS
###############################

o.set_config('environment:fallback:x_wind', 0.0)
o.set_config('environment:fallback:y_wind', 0.0)
o.set_config('environment:fallback:x_sea_water_velocity', 0.0)
o.set_config('environment:fallback:y_sea_water_velocity', 0.0)
o.set_config('environment:fallback:sea_floor_depth_below_sea_level', 100000.0)

# seed
o.set_config('seed:ocean_only',True) # keep only particles from the "frame" that are on the ocean
# drift
o.set_config('drift:advection_scheme','runge-kutta4') # or 'runge-kutta'
o.set_config('drift:current_uncertainty', 0.0 ) # note current_uncertainty can be used to replicate an horizontal diffusion spd_uncertain = sqrt(Kxy*2/dt)  
o.set_config('drift:lift_to_seafloor',True)

# horizontal and vertical diffusion
# Need to be turned off for backtracking
if False :
    Kxy = 0.1 #0.1176 # m2/s-1
    o.set_config('drift:horizontal_diffusivity',Kxy) # using new config rather than current uncertainty
    o.set_config('drift:vertical_mixing', True) 
    Kz = 0.001 # m2/s-1
    o.set_config('environment:fallback:ocean_vertical_diffusivity', Kz) # specify constant ocean_vertical_diffusivity in m2.s-1
    o.set_config('vertical_mixing:diffusivitymodel', 'constant') # constant >> use fall back values o.set_config('environment:fallback:ocean_vertical_diffusivity'] = Kz for all profile
    # can be environment (i.e. from reader), or  windspeed_Large1994 ,windspeed_Sundby1983, 
    o.set_config('vertical_mixing:timestep', 90.0)  # if some ocean_vertical_diffusivity!=0, turbulentmixing:timestep should be << 900 seconds

#o.set_config('drift:vertical_advection', True)
#o.disable_vertical_motion()  #Deactivate any vertical processes/advection"""


##############################
# CONFIG
# Adjusting some configuration
o.set_config('general:coastline_action', 'previous') # No settlement allowed here

##############################
#if False: # Only with plankton model
#    # plankton-specific config for developed module based on 
o.set_config('drift:max_age_seconds', 20*24*3600) # PLD spats (algae + mussels) = 20 days; PLD mussels = 5 weeks
#     o.set_config('biology:mortality_daily_rate', 0.0)    # 'float(min=0.0, max=100.0, default=0.05)', comment='Mortality rate (percentage of biomass dying per day)') 
#     #o.set_config('biology:min_settlement_age_seconds', 5*24*3600.0)  #'float(min=0.0, max=100.0, default=0.0)', comment='Minimum age before beaching can occur, in seconds')
#    
#     # Set the depth for the spats => close to the bottom but still within the euphotic zone.
#     #o.set_config('biology:vertical_position_daytime', -10.0)#'float(min=-1000.0, max=0.0, default=-5.0)',   comment='the depth a species is expected to inhabit during the day time, in meters, negative down') #
#     #o.set_config('biology:vertical_position_nighttime', -10.0) #'float(min=-1000.0, max=0.0, default=-1.0)', comment='the depth a species is expected to inhabit during the night time, in meters, negative down') #
#     #o.set_config('biology:vertical_migration_speed_constant',1e-4) #'float(min=0.0, max=1e-3, default=None)', comment=' Constant vertical migration rate (m/s)
#	
#    # to switch off the constant migration rate towards day or night time position, and use update_terminal_velocity_pelagicegg() :
#    # o.set_config('biology:vertical_migration_speed_constant',None) 
###############################

#o.list_config()

###############################
# RUN 
###############################
# Running model (until end of driver data)
o.run(stop_on_error=False,
      time_step=-900.0, 
      end_time = (datetime(2017,11,1) - timedelta(days = 20.0)),
      outfile= 'ninety_miles_beach_backward_nov_2017.nc',
      time_step_output = 3600.0)


###############################
# Plot
###############################

print(o)
#o.plot(fast=True)
o.animation(fast=True, color='z')
o.plot_vertical_distribution()