%% Post-processing of Opendrift outputs for the Ninety Miles Beach - part 1
%%% The model produces a single netcdf output file with dispersal
%%% trajectories.

%%% Settlement of particles is not explicitely represented in the model and
%%% therefore must be taken into account in this code.

clear
close all

% Addpath with output files, shapefiles and packages
addpath 'D:\Cours\Post_doc\Open_Drift_scripts\Matlab_post_processing\';
addpath 'D:\Cours\Post_doc\Moana Project\Maps_New_Zealand\'; % Path to add shapefile of coast
addpath 'D:\Cours\Post_doc\Open_Drift_scripts\Packages'; % Path to package for compute-hdc
addpath 'D:\Cours\Post_doc\Moana Project\Regional projects\Ninety miles beach\Open_drift_model';
addpath 'E:\Moana_Project\Regional_studies\Ninety_miles_beach\first_backtrack_outputs'; % Path to Opendrift output files

traj_output = 'ninety_miles_beach_backward_oct_2017.nc';
%ncdisp('ninety_miles_beach_backward_nov_2017.nc') % Check the attributes of the netcdf file


%% Load and sort variables

% Get length of variables
%vinfo = ncinfo(traj_output,'lat');
%varSize = vinfo.Size;
%split_nb = 10;
%sublength_varSize = varSize(2)/split_nb; % Create smaller chunks if necessary

% Get fill value
ncid = netcdf.open(traj_output,'NOWRITE');
varid = netcdf.inqVarID(ncid,'lat');
[noFillMode,fillValue] = netcdf.inqVarFill(ncid,varid);

% Load files
%disp_time=ncread(traj_output,'time'); % not needed here
%disp_ID=ncread(traj_output,'trajectory'); % not needed here
disp_lat=ncread(traj_output,'lat');disp_lat(disp_lat == fillValue) = NaN;
disp_lon=ncread(traj_output,'lon');disp_lon(disp_lon == fillValue) = NaN;
disp_depth=ncread(traj_output,'z');disp_depth(disp_depth == fillValue) = NaN;
disp_age=ncread(traj_output,'age_seconds');disp_age(disp_age == fillValue) = NaN;
%disp_status=ncread(traj_output,'status');

%% Remove values for age > 20 days
max_PLD = - (20 * 24 * 3600); % minimum PLD in seconds

for i = 1:size(disp_age, 1)
    for j = 1:size(disp_age, 2)
        if disp_age(i, j) < max_PLD
            disp_lat(i, j) = NaN;
            disp_lon(i, j) = NaN;
        end
    end
end

%% Identify high density area

% Define limits of dispersal
latitudes_nest = linspace(-35.6, -34.3, 300);
longitudes_nest = linspace(172.4, 173.4, 300);
% Compute the histogram
N = histcounts2(disp_lat(:,:), disp_lon(:,:), latitudes_nest, longitudes_nest);

% % Create Gaussian filter matrix for smoothing the sparse matrix:
% [xG, yG] = meshgrid(-5:5);
% sigma = 2.5;
% g = exp(-xG.^2./(2.*sigma.^2)-yG.^2./(2.*sigma.^2));
% g = g./sum(g(:));

% Plot contour
latitudes_nest_contour = linspace(-35.6, -34.3, 299);
longitudes_nest_contour = linspace(172.4, 173.4, 299);
[M,c] = contour(longitudes_nest_contour, latitudes_nest_contour, N/length(disp_lon), [0.1 0.15 0.2],'ShowText','on');
c.LineWidth = 1;
colorbar
hold on
axis equal;
set(gca, 'XLim', longitudes_nest([1 end]), 'YLim', latitudes_nest([1 end]) , 'YDir', 'normal');
hold on;
% Plot coastline
NZ = shaperead('New_Zealand_High_resolution_coast_for_plot.shp');
mapshow(NZ)

% Extract polygon with high density
high_density = find(M(1,:)==0.2);
Contour_high_density = M(:,high_density+1:end);
%pause
close

%% Plot the Probability density function

% Plot heatmap
imagesc(longitudes_nest, latitudes_nest, N/length(disp_lon)); %without smoothing
%imagesc(longitudes_nest, latitudes_nest, conv2(N/length(disp_lon), g, 'same'));
axis equal;
set(gca, 'XLim', longitudes_nest([1 end]), 'YLim', latitudes_nest([1 end]) , 'YDir', 'normal');
hold on;
colorbar
caxis([0 0.3])
hold on

% Add high density polygon
%h = plot(polyshape([Contour_high_density(1,:)],[Contour_high_density(2,:)]));
%h.LineStyle = '--';
%h.EdgeColor = 'red';

% Plot coastline
NZ = shaperead('New_Zealand_High_resolution_coast_for_plot.shp');
mapshow(NZ)

% Add title
title('Ninety Miles Beach Oct 2017 PLD 20 days');

% Save figure
savefig('pdf_Ninety_miles_beach_oct_2017_PLD20.fig')
%pause
close

%% Identify only particles between 10 and 20 days (mussels age)

% Remove values for age < 10 days
min_PLD = - (10* 24 * 3600); % minimum PLD in seconds
max_PLD = - (20 * 24 * 3600); % minimum PLD in seconds

for i = 1:size(disp_age, 1)
    for j = 1:size(disp_age, 2)
        if disp_age(i, j) > min_PLD
            disp_lat(i, j) = NaN;
            disp_lon(i, j) = NaN;
        elseif disp_age(i, j) < max_PLD
            disp_lat(i, j) = NaN;
            disp_lon(i, j) = NaN;
        end
    end
end

%% Plot the Probability density function

latitudes_nest = linspace(-35.6, -34.3, 300);
longitudes_nest = linspace(172.4, 173.4, 300);
% Compute the histogram
N = histcounts2(disp_lat(:,:), disp_lon(:,:), latitudes_nest, longitudes_nest);

% Create Gaussian filter matrix for smoothing the sparse matrix:
[xG, yG] = meshgrid(-5:5);
sigma = 2.5;
g = exp(-xG.^2./(2.*sigma.^2)-yG.^2./(2.*sigma.^2));
g = g./sum(g(:));

% Plot heatmap
imagesc(longitudes_nest, latitudes_nest, N/length(disp_lon)); %without smoothing
%imagesc(longitudes_nest, latitudes_nest, conv2(N/length(disp_lon), g, 'same'));
axis equal;
set(gca, 'XLim', longitudes_nest([1 end]), 'YLim', latitudes_nest([1 end]) , 'YDir', 'normal');
hold on;
colorbar
caxis([0 0.3])
hold on

% Plot coastline
NZ = shaperead('New_Zealand_High_resolution_coast_for_plot.shp');
mapshow(NZ)

% Add title
title('Ninety Miles Beach Oct 2017 PLD 10 to 20 days');

% Save figure
savefig('pdf_Ninety_miles_beach_Oct_2017_PLD_10to20.fig')
%pause
close

%% Plot depth distribution

%histogram(disp_depth) % Check the distribution of depth during dispersal


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Yearly average and standard deviation

clear
close

% Addpath with output files, shapefiles and packages
addpath 'D:\Cours\Post_doc\Open_Drift_scripts\Matlab_post_processing\';
addpath 'D:\Cours\Post_doc\Moana Project\Maps_New_Zealand\'; % Path to add shapefile of coast
addpath 'D:\Cours\Post_doc\Open_Drift_scripts\Packages'; % Path to package for compute-hdc
addpath 'D:\Cours\Post_doc\Moana Project\Regional projects\Ninety miles beach\Open_drift_model';
addpath 'E:\Moana_Project\Regional_studies\Ninety_miles_beach\first_backtrack_outputs'; % Path to Opendrift output files

% Loop througth the output files
traj_output=dir('E:\Moana_Project\Regional_studies\Ninety_miles_beach\first_backtrack_outputs\ninety_miles_beach_backward*');
heatmap = cell(1, 6);
heatmap(1,:) = {zeros(299, 299)};

for fi = 1:length(traj_output)
    
    % Get fill value
    ncid = netcdf.open(traj_output(fi).name,'NOWRITE');
    varid = netcdf.inqVarID(ncid,'lat');
    [noFillMode,fillValue] = netcdf.inqVarFill(ncid,varid);
    
    % Load files
    disp_time=ncread(traj_output(fi).name,'time');
    disp_ID=ncread(traj_output(fi).name,'trajectory');
    disp_lat=ncread(traj_output(fi).name,'lat');disp_lat(disp_lat == fillValue) = NaN;
    disp_lon=ncread(traj_output(fi).name,'lon');disp_lon(disp_lon == fillValue) = NaN;
    disp_age=ncread(traj_output(fi).name,'age_seconds');disp_age(disp_age == fillValue) = NaN;
    
    % Remove values for age < 10 days and > 20 days
    min_PLD = - (10* 24 * 3600); % minimum PLD in seconds
    max_PLD = - (20 * 24 * 3600); % maximum PLD in seconds
    
    for i = 1:size(disp_age, 1)
        for j = 1:size(disp_age, 2)
            if disp_age(i, j) > min_PLD
                disp_lat(i, j) = NaN;
                disp_lon(i, j) = NaN;
            elseif disp_age(i, j) < max_PLD
                disp_lat(i, j) = NaN;
                disp_lon(i, j) = NaN;
            end
        end
    end
    
    % Compute and save heatmap
    latitudes_nest = linspace(-35.6, -34.3, 300);
    longitudes_nest = linspace(172.4, 173.4, 300);
    N = histcounts2(disp_lat(:,:), disp_lon(:,:), latitudes_nest, longitudes_nest);
    heatmap{fi} = N;
    
end

% Yearly average PDF
sum_heatmap = zeros(299,299);
for i = 1:6
    monthly_heatmap = heatmap{i};
    sum_heatmap = sum_heatmap + monthly_heatmap;
end
yearly_heatmap = sum_heatmap/6;

% Convert H into 3-D numeric array
H2 = reshape(cell2mat(heatmap),299,299,6);
% Calculate std for each pixel position
std_heatmap = std(double(H2),[],3);

% Plot heatmap
imagesc(longitudes_nest, latitudes_nest, yearly_heatmap/length(disp_lon)); %without smoothing
%imagesc(longitudes_nest, latitudes_nest, conv2(N/length(disp_lon), g, 'same'));
axis equal;
set(gca, 'XLim', longitudes_nest([1 end1]), 'YLim', latitudes_nest([1 end]) , 'YDir', 'normal');
hold on;
colorbar
%caxis([0 0.3])
hold on
% Plot coastline
NZ = shaperead('New_Zealand_High_resolution_coast_for_plot.shp');
mapshow(NZ)
% Add title
title('Ninety Miles Beach 2017 average - PLD 10 to 20 days');
% Save figure
savefig('pdf_Ninety_miles_beach_yearly_2017_PLD_10to20.fig')
pause
close

% Plot heatmap Standard deviation
imagesc(longitudes_nest, latitudes_nest, std_heatmap/length(disp_lon)); %without smoothing
%imagesc(longitudes_nest, latitudes_nest, conv2(N/length(disp_lon), g, 'same'));
axis equal;
set(gca, 'XLim', longitudes_nest([1 end]), 'YLim', latitudes_nest([1 end]) , 'YDir', 'normal');
hold on;
colorbar
%caxis([0 0.3])
hold on
hold on
% Plot coastline
NZ = shaperead('New_Zealand_High_resolution_coast_for_plot.shp');
mapshow(NZ)
% Add title
title('Ninety Miles Beach 2017 std - PLD 10 to 20 days');
pause
close

%% Correction of data based on released area
% The problem with the previous analysis is that particles are considered
% similar independently of their release location. In reality, mussel spats
% arrive in the South of the beach and then drift North. Therefore, mussels
% collected on the North of the beach are older than the South ones. This
% might cause the high density of trajectories observed offshore in the
% middle of the beach.

% Using Alfaro et al. 2004 (sampling of mussel spat sizes in South and
% North of the beach) and Hickman 1979 (average growth rate for mussels
% <60mm), we can deduce that a wide estimate of the age difference from one
% extreme of the beach to the other is 8.4 days, and a conservative
% estimate is 5 days (smallest and middle size class have the same trend
% in Alfaro et al. 2004).

clear
close

% Use of the conservative estimate: 5.5 days; or use large estimate: 8.4
% days
seconds_per_day = 3600*24;
difference_age_South_North = 5.5; % days
difference_lat_South_North = -35.2-(-34.5); % degree
dispersal_rate = (difference_age_South_North*seconds_per_day/difference_lat_South_North); % Seconds per degree

% Import output files and apply dispesal rate correction
% Addpath with output files, shapefiles and packages
addpath 'D:\Cours\Post_doc\Open_Drift_scripts\Matlab_post_processing\';
addpath 'D:\Cours\Post_doc\Moana Project\Maps_New_Zealand\'; % Path to add shapefile of coast
addpath 'D:\Cours\Post_doc\Open_Drift_scripts\Packages'; % Path to package for compute-hdc
addpath 'D:\Cours\Post_doc\Moana Project\Regional projects\Ninety miles beach\Open_drift_model';
addpath 'E:\Moana_Project\Regional_studies\Ninety_miles_beach\first_backtrack_outputs'; % Path to Opendrift output files

% Loop througth the output files
traj_output=dir('E:\Moana_Project\Regional_studies\Ninety_miles_beach\first_backtrack_outputs\ninety_miles_beach_backward*');
heatmap = cell(1, 6);
heatmap(1,:) = {zeros(299, 299)};
for fi = 1:length(traj_output)
    
    % Get fill value
    ncid = netcdf.open(traj_output(fi).name,'NOWRITE');
    varid = netcdf.inqVarID(ncid,'lat');
    [noFillMode,fillValue] = netcdf.inqVarFill(ncid,varid);
    
    % Load files
    disp_time=ncread(traj_output(fi).name,'time');
    disp_ID=ncread(traj_output(fi).name,'trajectory');
    disp_lat=ncread(traj_output(fi).name,'lat');disp_lat(disp_lat == fillValue) = NaN;
    disp_lon=ncread(traj_output(fi).name,'lon');disp_lon(disp_lon == fillValue) = NaN;
    disp_age=ncread(traj_output(fi).name,'age_seconds');disp_age(disp_age == fillValue) = NaN;
    
    % Remove values for age < 10 days and > 20days
    min_PLD = - (10* 24 * 3600); % minimum PLD in seconds
    max_PLD = - (20 * 24 * 3600); % maximum PLD in seconds
    
    for j = 1:size(disp_age, 2) % Loop through particles
        
        % Find starting Latitude to apply age correction using dispersal
        % rate
        validIndices = ~isnan(disp_lat(:,j));
        StartLat = disp_lat(validIndices,j);
        StartLat = StartLat(1);
        
        for i = 1:size(disp_age, 1) % Loop though time-steps
            if disp_age(i, j) > min_PLD - dispersal_rate * (-35.2 - StartLat) % Correction time settlement
                disp_lat(i, j) = NaN;
                disp_lon(i, j) = NaN;
            elseif disp_age(i, j) < max_PLD
                disp_lat(i, j) = NaN;
                disp_lon(i, j) = NaN;
            end
        end
    end
    
    % Compute and save heatmap
    latitudes_nest = linspace(-35.6, -34.3, 300);
    longitudes_nest = linspace(172.4, 173.4, 300);
    N = histcounts2(disp_lat(:,:), disp_lon(:,:), latitudes_nest, longitudes_nest);
    heatmap{fi} = N;
    
end

% Yearly average PDF
sum_heatmap = zeros(299,299);
for i = 1:6
    monthly_heatmap = heatmap{i};
    sum_heatmap = sum_heatmap + monthly_heatmap;
end
yearly_heatmap = sum_heatmap/6;

% Convert H into 3-D numeric array
H2 = reshape(cell2mat(heatmap),299,299,6);
% Calculate std for each pixel position
std_heatmap = std(double(H2),[],3);

% Plot heatmap
imagesc(longitudes_nest, latitudes_nest, yearly_heatmap/length(disp_lon)); %without smoothing
%imagesc(longitudes_nest, latitudes_nest, conv2(N/length(disp_lon), g, 'same'));
axis equal;
set(gca, 'XLim', longitudes_nest([1 end]), 'YLim', latitudes_nest([1 end]) , 'YDir', 'normal');
hold on;
colorbar
caxis([0 0.15])
hold on
% Plot coastline
NZ = shaperead('New_Zealand_High_resolution_coast_for_plot.shp');
mapshow(NZ)
% Add title
title('Ninety Miles Beach 2017 average - corrected PLD (5.5 days lag)');
% Save figure
savefig('pdf_Ninety_miles_beach_yearly_2017_corrected_PLD.fig')
%pause
close

% Plot heatmap Standard deviation
imagesc(longitudes_nest, latitudes_nest, std_heatmap/length(disp_lon)); %without smoothing
%imagesc(longitudes_nest, latitudes_nest, conv2(N/length(disp_lon), g, 'same'));
axis equal;
set(gca, 'XLim', longitudes_nest([1 end]), 'YLim', latitudes_nest([1 end]) , 'YDir', 'normal');
hold on;
colorbar
%caxis([0 0.3])
hold on
hold on
% Plot coastline
NZ = shaperead('New_Zealand_High_resolution_coast_for_plot.shp');
mapshow(NZ)
% Add title
title('Ninety Miles Beach 2017 std - corrected PLD (5.5 days lag)');
%pause
close

%% Isolate area of high trajectory density

latitudes_nest_contour = linspace(-35.6, -34.3, 299);
longitudes_nest_contour = linspace(172.4, 173.4, 299);
[M,c] = contour(longitudes_nest_contour, latitudes_nest_contour, yearly_heatmap/length(disp_lon), [0.01 0.05 0.09],'ShowText','on');
c.LineWidth = 1;
colorbar
hold on
set(gca, 'XLim', longitudes_nest([1 end]), 'YLim', latitudes_nest([1 end]) , 'YDir', 'normal');
axis equal;
% Add bathymetry to remove particles that would be in deep water
% Plot bathymetry
addpath('D:\Cours\Post_doc\Moana Project\Maps_New_Zealand\NZBathy_2016_shape');
Bathy_NZ = shaperead('NZBathy_2016_contours_50m_for_plot.shp');
mapshow(Bathy_NZ, 'Color',[1 0 0])
% Plot coastline
NZ = shaperead('New_Zealand_High_resolution_coast_for_plot.shp');
mapshow(NZ)
title('Isobath and isodensity of trajectories in 2017 for corrected PLD')
pause
close

% Extract polygon with high density
high_density = find(M(1,:)==0.09);
for i = 1:size(high_density,2)-1
    Contour_high_density{i} = M(:,high_density(i)+1:high_density(i+1)-1);
end
for i = size(high_density,2)
    Contour_high_density{i} = M(:,high_density(i)+1:end);
end

% Create geoshape object and plot
for i = 1:size(high_density,2)
    s(i) = geoshape(Contour_high_density{i}(2,:), Contour_high_density{i}(1,:));
end
geoshow(s)
set(gca, 'XLim', longitudes_nest([1 end]), 'YLim', latitudes_nest([1 end]) , 'YDir', 'normal');
hold on
% Add bathymetry to remove particles that would be in deep water
% Plot bathymetry
addpath('D:\Cours\Post_doc\Moana Project\Maps_New_Zealand\NZBathy_2016_shape');
Bathy_NZ = shaperead('NZBathy_2016_contours_50m_for_plot.shp');
mapshow(Bathy_NZ, 'Color',[1 0 0])
% Plot coastline
NZ = shaperead('New_Zealand_High_resolution_coast_for_plot.shp');
mapshow(NZ)
pause
close

% Remove areas located in deep water
% None for 2017

%% Create release points for next simulations: I gave up and created this in QGIS

% Save the areas in a shapefile
shapewrite(s, 'high_density_areas_1st_backtrack_avg_2017_corrected_PLD')
poly = shaperead('high_density_areas_1st_backtrack_avg_2017_corrected_PLD.shp');

% Add buffer to areas
bufferwidth = 0.005;
direction = 'outPlusInterior';
% Need to inverse the latitude in order to create the buffer (stupid function)
[latb, lonb] = bufferm(-s.Latitude, s.Longitude, bufferwidth, direction); % Create a 500m (0.005 degree) buffer around polygons