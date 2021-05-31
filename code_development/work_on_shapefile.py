#!/usr/bin/python

# Import necessary modules
import shapefile
from shapely.geometry import Polygon, Point
import numpy as np
import math
import pymap3d as pm
from joblib import Parallel, delayed
import multiprocessing
import numba

# Import release polygons
shp = shapefile.Reader('./habitat/Test_habitat.shp')
bins = shp.shapes()
records = shp.records()

# Define two points to test the algorithm:
p1 = Point(174.71091, -38.07374) # Centroid of one of the polygons
p2 = Point(175.71091, -37.07374) # Pt not inside of any polygon of the shapefile
lons = [174.71091, 174.71091]
lats = [-38.07374, -37.07374]

# Loop to check within polygon or not
#def within_habitat(shp, current_lat, current_lon) 

for i in range(len(lons)):
    pt = Point(lons[i], lats[i])
    for index, item in enumerate(shp):
        pts = bins[index].points # Access one polygon
        #pts = [list(elem) for elem in pts]
        poly = Polygon(pts)
        in_area = pt.within(poly)
        if in_area == True:
            print("The habitat id is", index)
        else:
            print("No habitat")


####################################################
# Second part: try to find centroid of polygon
def centroid(vertexes):
     _x_list = [vertex [0] for vertex in vertexes]
     _y_list = [vertex [1] for vertex in vertexes]
     _len = len(vertexes)
     _x = sum(_x_list) / _len
     _y = sum(_y_list) / _len
     return(_x, _y)
 
    
centers = np.zeros((4,2))
for index, item in enumerate(shp): 
    centers[index, 0] = centroid(bins[index].points)[0]
    centers[index, 1] = centroid(bins[index].points)[1]

    
######################################################  
# Third part: find nearest habitat

# Haversine formula
@numba.jit(nopython=True)
def haversine_distance(s_lat,s_lng,e_lat,e_lng):
    # approximate radius of earth in km
    R = 6373.0
    s_lat = np.deg2rad(s_lat)                    
    s_lng = np.deg2rad(s_lng)     
    e_lat = np.deg2rad(e_lat)                       
    e_lng = np.deg2rad(e_lng)
    d = np.sin((e_lat - s_lat)/2)**2 + \
        np.cos(s_lat)*np.cos(e_lat) * \
        np.sin((e_lng - s_lng)/2)**2
    return 2 * R * np.arcsin(np.sqrt(d))

@numba.jit(nopython=True)
def haversine_angle(lon1, lat1, lon2, lat2):
    rlat1 = np.deg2rad(lat1)
    rlat2 = np.deg2rad(lat2)
    rlon1 = np.deg2rad(lon1)
    rlon2 = np.deg2rad(lon2)
    X = math.cos(rlat2)*math.sin(rlon2-rlon1)
    Y = math.cos(rlat1)*math.sin(rlat2)-math.sin(rlat1)*math.cos(rlat2)*math.cos(rlon2-rlon1)
    return math.atan2(Y,X)

def nearest_habitat(lon, lat, centers):
    dist = np.zeros(len(centers))
    dist = haversine_distance(lon, lat, centers[:, 0], centers[:, 1])
    nearest_center = np.argmin(dist)
    return nearest_center, min(dist)
 
######################################################
# Fourth part: Runge-Kutta

# Runge-Kutta: https://www.geeksforgeeks.org

# Python program to implement Runge Kutta method
# A sample differential equation "dy / dx = (x - y)/2"
@numba.jit(nopython=True)
def dydx(x, y):
    return ((x - y)/2)

@numba.jit(nopython=True)
def rungeKutta(x0, y0, x, h):
    # Count number of iterations using step size or
    # step height h
    n = (int)((x - x0)/h)
    # Iterate for number of iterations
    y = y0
    for i in range(1, n + 1):
        "Apply Runge Kutta Formulas to find next value of y"
        k1 = h * dydx(x0, y)
        k2 = h * dydx(x0 + 0.5 * h, y + 0.5 * k1)
        k3 = h * dydx(x0 + 0.5 * h, y + 0.5 * k2)
        k4 = h * dydx(x0 + h, y + k3)
 
        # Update next value of y
        y = y + (1.0 / 6.0)*(k1 + 2 * k2 + 2 * k3 + k4)
 
        # Update next value of x
        x0 = x0 + h
    return y


       
######################################################
# Fith part: IMB model with orientation toward nearest habitat

# Define currents
u = 1.0
v = 0.0

# Define timeframe
PLD = 30.0 # max time of dispersal in days
competency = 10.0 # possibility to settle (until end of PLD)
timestep = 1

# Define variables of orientation
Beta = 10000 # Distance of detection of the habitat
Kappa = 5 # Precision of orientation
swimming_speed = 2 # random for now
horizon_diff = 0.1 # horizontal diffusion when not orienting

# Import habitat
shp = shapefile.Reader('./habitat/Test_habitat.shp')
bins = shp.shapes()
records = shp.records()
centers = np.zeros((4,2))
for index, item in enumerate(shp): 
    centers[index, 0] = centroid(bins[index].points)[0]
    centers[index, 1] = centroid(bins[index].points)[1]

# Define number of particles to release
#Nb = 100

# Define starting position
lon_origin = 174.71091
lat_origin = -38.07374


# Create model in function
def run_Lagrangian_orientation_model(initial_lat, initial_lon, centers_habitat, PLD, competency, Beta, Kappa, Nb):
    
    lon_particle = lon_origin # for first step
    lat_particle = lat_origin # for first step
    lon_particle_old = lon_origin 
    lat_particle_old = lat_origin
    
    # Initiate trajectories
    traj_x, traj_y, h = pm.geodetic2ecef(lat_origin,lon_origin,0)
    traject_x = np.zeros((int(PLD/timestep), Nb))
    traject_y = np.zeros((int(PLD/timestep), Nb))
    traject_z = np.zeros((int(PLD/timestep), Nb))
    traject_x[0] = traj_x
    traject_y[0] = traj_y
    traject_z[0] = h
    for N in range(Nb):
        # Start loop
        for i in np.arange(0, PLD, timestep):
            
            # Update age particle
            if i == 0:
                age_particle = 0
            else:
                age_particle = age_particle + timestep
            
            # Case where particle too young to orient
            if age_particle < competency:
                randm_vel = math.sqrt((2*horizon_diff)/timestep)
                u_diff = np.random.normal(0, 1)*randm_vel
                v_diff = np.random.normal(0, 1)*randm_vel
                u_particle = u_diff
                v_particle = v_diff
            else:
                habitat_near = nearest_habitat(lon_particle, lat_particle, centers) # habitat_near(0:ID centroid, 1:distance in km)
                if habitat_near[1] > Beta:
                    # Case where particle too far to orient
                    randm_vel = math.sqrt((2*horizon_diff)/timestep)
                    u_diff = np.random.normal(0, 1)*randm_vel
                    v_diff = np.random.normal(0, 1)*randm_vel
                    u_particle = u_diff
                    v_particle = v_diff
                else:
                    # Case where particle close enough and old enough to orient
                    # Strength of orientation (depend on distance to the habitat)
                    d = 1 - (habitat_near[1]/Beta)
                    
                    # Compute direction of nearest habitat. See Staaterman et al., 2012
                    theta_pref = -haversine_angle(lon_particle, lat_particle, centers[habitat_near[0]][0], centers[habitat_near[0]][1]) 
                    
                    # Compute direction from previous timestep
                    theta_current = haversine_angle(lon_particle_old, lat_particle_old, lon_particle, lat_particle)
                    
                    # Mean turning angle
                    mu = -d * (theta_current - theta_pref)
                    
                    # New direction randomly selected in Von Mises distribution
                    ti  = np.random.vonmises(0, Kappa)
                    theta = ti - theta_current - mu
                    
                    # Compute swimming speed
                    swimming_speed_age = swimming_speed * age_particle
                    
                    # Compute uorient and vorient
                    uorient = swimming_speed_age*math.cos(theta)
                    vorient = swimming_speed_age*math.sin(theta)
                    u_particle = uorient
                    v_particle = vorient
                    
            # Compute the movement of the particle due to the currents
            x0, y0, h = pm.geodetic2ecef(lat_particle, lon_particle, alt=0)
            u_there = rungeKutta(x0, y0, u, 4)
            v_there = rungeKutta(x0, y0, v, 4)
            # Add the movement of the particle due to the diffusion or the orientation
            U = u_there + u_particle
            V = v_there + v_particle
            
            # Compute the trajectory of the particle
            traj_x = traj_x + timestep * U
            traj_y = traj_y + timestep * V
            traject_x[int(i), N] = traj_x
            traject_y[int(i), N] = traj_y
            traject_z[int(i), N] = h
            
            # Update position of particle
            lat_particle_old = lat_particle
            lon_particle_old = lon_particle
            lat_particle, lon_particle, h = pm.ecef2geodetic(traj_x, traj_y, z=h)
        
        # End of loop
    return traject_x, traject_y, traject_z
    
    
# Run model
num_cores = multiprocessing.cpu_count()
traject_x, traject_y, traject_z = Parallel(n_jobs=num_cores)(delayed(run_Lagrangian_orientation_model)(lat_origin, lon_origin, centers, PLD, competency, Beta, Kappa, Nb) for Nb in range(2))


################################################
# Plot of the trajectories










