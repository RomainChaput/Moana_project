%% Post-processing of Opendrift outputs for the Ninety Miles Beach - part 2
%%% The model produces a single netcdf output file with dispersal
%%% trajectories.

%%% Settlement of particles is not explicitely represented in the model and
%%% therefore must be taken into account in this code.

%%% The goal is to identify the parent populations of the mussels that
%%% settle on the Ninety Miles Beach via the mussel spats. High density PDF
%%% have been identified in the previous step, now we look to identify the
%%% reefs that feed to these high density PDF. Note that all PDF areas
%%% might not be connected to reef areas.

clear
close all
%parpool(12) % Initiate parallel pool

tic
% Addpath with output files, shapefiles and packages
addpath 'D:\Cours\Post_doc\Open_Drift_scripts\Matlab_post_processing\';
addpath 'D:\Cours\Post_doc\Moana Project\Maps_New_Zealand\'; % Path to add shapefile of coast
addpath 'D:\Cours\Post_doc\Open_Drift_scripts\Packages'; % Path to package for compute-hdc
addpath 'D:\Cours\Post_doc\Moana Project\Regional projects\Ninety miles beach\Open_drift_model';
addpath 'E:\Moana_Project\Regional_studies\Ninety_miles_beach\second_backtrack_outputs'; % Path to Opendrift output files
addpath 'D:\Cours\Post_doc\Moana Project\Regional projects\Ninety miles beach\Map_rocky_shore_green_mussel';

% Add shapefiles
Mussel_habitat = shaperead('mussel_habitat_rough_01deg_moana_mask.shp');
Release_Areas = shaperead('high_density_areas_buffer_2017_corrected_PLD.shp');
NZ = shaperead('New_Zealand_High_resolution_coast_for_plot.shp');


%% Loop through the output files

traj_output=dir('E:\Moana_Project\Regional_studies\Ninety_miles_beach\second_backtrack_outputs\ninety_miles_beach_backward*');

disp_release = cell(1, length(traj_output)); % Pre-allocate variable for release areas connected
frequency_habitat = cell(1, length(traj_output)); % Pre-allocate variable
frequency_release = cell(1, length(traj_output)); % Pre-allocate variable
disp_settlement = cell(1, length(traj_output)); % Pre-allocate the variable settlement

for fi = 1:length(traj_output)
    
    %% Load and sort variables
    
    % Get fill value
    ncid = netcdf.open(traj_output(fi).name,'NOWRITE');
    %ncdisp(traj_output(fi).name) % Check the attributes of the netcdf file
    varid = netcdf.inqVarID(ncid,'lat');
    [noFillMode,fillValue] = netcdf.inqVarFill(ncid,varid);
    
    % Load files
    disp_lat=ncread(traj_output(fi).name,'lat');disp_lat(disp_lat == fillValue) = NaN;
    disp_lon=ncread(traj_output(fi).name,'lon');disp_lon(disp_lon == fillValue) = NaN;
    %disp_depth=ncread(traj_output(fi).name,'z');disp_depth(disp_depth == fillValue) = NaN;
    disp_age=ncread(traj_output(fi).name,'age_seconds');disp_age(disp_age == fillValue) = NaN;
    
    %% Remove values for age > 35 days => 5 weeks of dispersal for blue mussels
    
    max_PLD = - (35 * 24 * 3600); % maximum PLD in seconds
    
    for i = 1:size(disp_age, 1)
        for j = 1:size(disp_age, 2)
            if disp_age(i, j) < max_PLD
                disp_lat(i, j) = NaN;
                disp_lon(i, j) = NaN;
            end
        end
    end
    
    %% %% Identify high density area (does not work well here: high density in the release areas mostly)
    
    %% %% Plot subset of trajectories: Just to check
    %
    % for i = 1:2100:200000
    %     plot(disp_lon(:,i), disp_lat(:,i))
    %     hold on
    % end
    % hold on
    % latitudes_nest = linspace(-35.6, -34.3, 300);
    % longitudes_nest = linspace(172.4, 173.4, 300);
    % set(gca, 'XLim', longitudes_nest([1 end]), 'YLim', latitudes_nest([1 end]) , 'YDir', 'normal');
    % hold on;
    % % Plot coastline
    % NZ = shaperead('New_Zealand_High_resolution_coast_for_plot.shp');
    % mapshow(NZ)
    % pause
    % close
    
    
    %% Settlement of larvae: The release point since we backtrack the particles
    
    % Settlement time-window: 3 to 5 weeks (Norrie et al., 2020, Jeffs et al., 1999)
    Competency_start = -21 * 24* 3600; % Beginning of the settlement in seconds
    
    % Locate the habitat: rough map based on the mask of Moana backbone 5km resolution model
    %Mussel_habitat =
    %shaperead('polygons_90_miles_beach_rocky_areas_200m_buffer.shp'); % Fine
    %resolution habitat
    %Mussel_habitat = shaperead('mussel_habitat_rough_01deg_moana_mask.shp');
    
    % % Plot to check
    % mapshow(Mussel_habitat)
    % hold on
    % % Plot coastline
    % mapshow(NZ)
    
    % Find larvae in habitat
    in_habitat = zeros(size(disp_lat,1),size(disp_lat,2));
    
    %tic
    for i = 1:size(disp_lat,2)
        Competent = find(disp_age(:,i) < Competency_start);
        End_life = find(disp_age(:,i) <= max_PLD);
        parfor(j = Competent(1):End_life(1), 12)
            if ~isnan(disp_lat(j,i))
                for z = 1:size(Mussel_habitat,1)
                    if inpolygon(disp_lat(j,i), disp_lon(j,i), [Mussel_habitat(z).Y], [Mussel_habitat(z).X]) == 1
                        in_habitat(j,i) = Mussel_habitat(z).Area_code;
                        break
                    end
                end
            end
            
        end
    end
    %toc
    
    % Count frequency of settlement per site
    frequency_habitat{fi} = tabulate(in_habitat(:));
    disp_settlement{fi} = in_habitat; 
    
    
    %% Identify starting area that give settlers
    
    % Import shapefile of high density PDF from first backtrack
    %Release_Areas = shaperead('high_density_areas_buffer_2017_corrected_PLD.shp');
    Release_poly = zeros(1, size(disp_age,2));
    
    for i = 1:size(disp_age,2)
        if max(in_habitat(:,i)) > 0
            B = find(disp_age(:,i) == max(disp_age(:,i)));
            for z = 1:size(Release_Areas,1)
                if inpolygon(disp_lat(B,i), disp_lon(B,i), [Release_Areas(z).Y], [Release_Areas(z).X]) == 1
                    Release_poly(i) = Release_Areas(z).FID;
                end
            end
        end
    end
    
    % Count frequency of settlement per site
    %a2 = unique(Release_poly);
    frequency_release{fi} = tabulate(Release_poly);
    disp_release{fi} = Release_poly;
    
end
toc

clear disp_age
clear disp_depth
clear disp_lat
clear disp_lon

%% Work on the release and settlement habitats
% Start here and load .mat matrices to save time

% Average frequency of habitat in 2017
nb_habitat = 18;
avg_freq_habitat = zeros(nb_habitat+1,3);
avg_freq_habitat(:,1) = (0:nb_habitat)';

for i = 1:size(frequency_habitat,2)
    for j = 1:size(frequency_habitat{i},1)
        for z = 1:size(avg_freq_habitat, 1)
            if frequency_habitat{i}(j,1) == avg_freq_habitat(z,1)
                avg_freq_habitat(z,2) = avg_freq_habitat(z,2) + frequency_habitat{i}(j,2);
                avg_freq_habitat(z,3) = avg_freq_habitat(z,3) + frequency_habitat{i}(j,3);
            end
        end
    end
end

sum_freq_habitat = sum(avg_freq_habitat(2:19,2));
avg_freq_habitat(:,4) = avg_freq_habitat(:,2)/ sum_freq_habitat;

% Plot the habitat frequency
figure
bar(avg_freq_habitat(2:19,4))

% Average frequency of release areas in 2017
nb_release = 6;
avg_freq_release = zeros(nb_release,3);
avg_freq_release(:,1) = [0 2 3 4 5 7]';

for i = 1:size(frequency_release,2)
    for j = 1:size(frequency_release{i},1)
        for z = 1:size(avg_freq_release, 1)
            if frequency_release{i}(j,1) == avg_freq_release(z,1)
                avg_freq_release(z,2) = avg_freq_release(z,2) + frequency_release{i}(j,2);
                avg_freq_release(z,3) = avg_freq_release(z,3) + frequency_release{i}(j,3);
            end
        end
    end
end

sum_freq_release = sum(avg_freq_release(1:6,2));
avg_freq_release(:,4) = avg_freq_release(:,2)/ sum_freq_release;

figure
% Plot the habitat frequency
bar(avg_freq_release(1:nb_release,1),avg_freq_release(1:nb_release,4))
%close

% Plot both release and settlement

% Why not a connectivity matrix? => Settlement not turned on during
% backtracking so only look a the position of the particles during
% dispersal. I also count more larvae that travel on top of large areas
% (spend more time over them) than on small areas.
% => provide basis for next analysis
figure
subplot(2,1,1)
bar(avg_freq_release(1:nb_release,1),avg_freq_release(1:nb_release,4))
xlabel('Release site ID')
ylabel('Frequency of settlers release site')
subplot(2,1,2)
bar(avg_freq_habitat(2:18,4))
xlabel('Mussel habitat ID')
ylabel('Frequency of mussel habitat with settlers')

% Frequency by area: habitat
earthradius = almanac('earth','radius');
area = zeros(size(Mussel_habitat,1),2);
for i = 1:size(Mussel_habitat,1)
    area(i,1) = areaint(Mussel_habitat(i).Y,Mussel_habitat(i).X,earthradius);
    area(i,2) = Mussel_habitat(i).Area_code;
end
area_habitat(2:19,1) = accumarray(area(:,2),area(:,1));

% Frequency by area: release
area = zeros(size(Release_Areas,1),2);
for i = 1:size(Release_Areas,1)
    area(i,1) = areaint(Release_Areas(i).Y,Release_Areas(i).X,earthradius);
    area(i,2) = Release_Areas(i).FID;
end
area_release = sortrows(area,2);

% Plot frequencies per area
avg_freq_habitat(:,5) = avg_freq_habitat(1:19,2)./area_habitat; % Settlement per square kilometer of habitat
avg_freq_release(:,5) = avg_freq_release(:,2)./area_release(:,1); % Released settlers per square kilometer of release 

figure
subplot(2,1,1)
bar(avg_freq_release(1:nb_release,1),avg_freq_release(1:nb_release,5))
xlabel('Release site ID')
ylabel('Settlers per area of release site (larvae/km^{2})')
subplot(2,1,2)
bar(avg_freq_habitat(2:18,5))
xlabel('Mussel habitat ID')
ylabel('Settlers per area of habitat (larvae/km^{2})')

