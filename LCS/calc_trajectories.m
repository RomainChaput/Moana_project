%% Script to calculate the trajectories of particles on top of physical current nest files
% => works with calc_ftle.m to compute the FTLE and find the LCS

% Author: Philippe Miron, Oct 03 2019. Modified for Moana Backbone model
clc; clear;

% parameters
N = 1000; % particles across domain
seconds_per_day = 24 * 60 * 60; % [s]
seconds_per_hour = 3600;
t0 = 1 * seconds_per_day; % Start time
T = -14  * seconds_per_day; % [s] integration time
t1 = t0 + T;
% matlab functions for analysis
addpath('D:\Cours\Post_doc\Open_Drift_scripts\Matlab_post_processing')

% folders
prefixIN  = 'input/';
basefile  = {'nz5km_his_subset_199401000000.nc'}; % day in between
prefixOUT = 'output/';

% load grid data
addpath('./input')
latu = ncread('nz5km_his_subset_199401000000.nc', 'lat_u');
latu = latu(1,:);
lonu = ncread('nz5km_his_subset_199401000000.nc', 'lon_u');
lonu = lonu(:,1);
latu = double(latu)';
lonu = double(lonu);
%lonu = wrapTo180(lonu);
lonu_origin = min(lonu);
latu_origin = min(latu);
[xudat, yudat] = sph2xy(lonu, lonu_origin, latu, latu_origin); % [m]
clear lonu latu

latv = ncread('nz5km_his_subset_199401000000.nc', 'lat_v');
latv = latv(1,:);
lonv = ncread('nz5km_his_subset_199401000000.nc', 'lon_v');
lonv = lonv(:,1);
latv = double(latv)';
lonv = double(lonv);
%lonv = wrapTo180(lonv);
lonv_origin = min(lonv);
latv_origin = min(latv);
[xvdat, yvdat] = sph2xy(lonv, lonv_origin, latv, latv_origin); % [m]
clear lonv latv

% load velocity data
addpath('E:/Moana_Project/MOANA_nests/hourly_nests')
f=dir('E:/Moana_Project/MOANA_nests/hourly_nests/nz5km_his_subset*');
nz = 1; % depth layer: 1 = surface
tdat = (1:1800)';
u = zeros(length(xudat), length(yudat), length(tdat));
v = zeros(length(xvdat), length(yvdat), length(tdat));
for fi=1:length(tdat)
    uui = ncread(f(fi).name,'u');
    u(:,:,fi) = uui(:,:,nz);
    vvi = ncread(f(fi).name,'v');
    v(:,:,fi) = vvi(:,:,nz);
end
masku = isnan(u);
maskv = isnan(v);
u(masku) = 0;
v(maskv) = 0;

% switch the time to seconds for integration of trajectories
% x,y in meters and t in seconds
tdat = tdat * seconds_per_hour * 1;

% interpolation functions for the velocity fields
[x, y, t] = ndgrid(xudat, yudat, tdat);
Fu = griddedInterpolant(x, y, t, u, 'linear', 'none');
[x, y, t] = ndgrid(xvdat, yvdat, tdat);
Fv = griddedInterpolant(x, y, t, v, 'linear', 'none');

% setting the initial grid for the trajectories
% this is independant of the velocity grid
lon0min = 172.5; % [°]
lon0max = 177.5; % [°]
lat0min = -44.5; % [°]
lat0max = -40.5; % [°]
[x0min, y0min] = sph2xy(lon0min, lonu_origin, lat0min, latu_origin);
[x0max, y0max] = sph2xy(lon0max, lonu_origin, lat0max, latu_origin);
fac = (x0max-x0min)/abs(y0max-y0min);
if fac > 1
    Nx0 = N;
    Ny0 = fix(N/fac);
else
    Ny0 = N;
    Nx0 = fix(N*fac);
end
Nxy0 = Nx0*Ny0;
x0 = linspace(x0min, x0max, Nx0);
y0 = linspace(y0min, y0max, Ny0);
[X0, Y0] = ndgrid(x0, y0);
X0 = X0(:);
Y0 = Y0(:);

%%
% if you have a large number of points it might be good to calculate only
% pack of ~10^4 trajectories at each time using ode4
Nseries = ceil(Nxy0/1e5);


for li = 10*24:1800
    t0 = li * 1 * seconds_per_hour; % Start time
    T = -10 * seconds_per_day; % [s] integration time
    t1 = t0 + T;
    
    tspan = t0:sign(T)*seconds_per_hour:t1;
    xf = zeros(Nxy0,1);
    yf = zeros(Nxy0,1);
    
    for i = 1:Nseries
        disp(['Serie = ' num2str(i) '/' num2str(Nseries)])
        range = (i-1)*1e5+1:min(i*1e5, Nxy0);
        xy0 = [X0(range); Y0(range)];
        
        % ode4 with a constant timestep
        xyf = ode4(@(t, y) uv(t, y, Fu, Fv), tspan, xy0);
        A = ~isnan(xyf);
        Indices = arrayfun(@(x) find(A(:, x), 1, 'last'), 1:size(xyf, 2));
        nonNanValues = arrayfun(@(x,y) xyf(x,y), Indices, 1:size(xyf, 2));
        xf(range) = nonNanValues(end, 1:end/2)';
        yf(range) = nonNanValues(end, end/2+1:end)';
    end
    
%     % % plot all trajectories (only to validate with a small N value!)
%     h=figure(1); clf;
%     set(h,'units','normalized','position',[0.7 0 0.3 1]);
%     hold on
%     for i = 1:50%Nxy0
%         txi = xyf(:, i);
%         tyi = xyf(:, end/2+i);
%         txi = txi(~isnan(txi));
%         tyi = tyi(~isnan(tyi));
%     
%         [txd, tyd] = xy2sph(txi, lonu_origin, tyi, latu_origin);
%         plot(txd, tyd);
%     end
%     [x0d, y0d] = xy2sph(xy0(1:end/2), lonu_origin, xy0(end/2+1:end), latu_origin);
%     scatter(x0d, y0d)
%     %gom_coasts(-83.35, -79.10, 22.77, 26.05)
%     hold off
    
    % save
    eval(['save ' prefixOUT 'xyt-' num2str(abs(t0)) '-' num2str(abs(t1)) ' T x0 y0 xf yf'])
end