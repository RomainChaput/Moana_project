%% Post-processing of Opendrift outputs for the Ninety Miles Beach - part 3
%%% The opendrift model produces netcdf output files with dispersal
%%% trajectories

%%% Settlement of particles is not explicitely represented in the model and
%%% therefore must be taken into account in this code

%%% The goal is to identify the parent populations of the mussels that
%%% settle on the Ninety Miles Beach via the mussel spats. High density PDF
%%% have been identified in the previous step, and potential source
%%% populations have bee identified. Now, we ran a foreward tracking
%%% experiement to validate the population connectivity.

%%% This code is to be run on NeSI

clear
close all
% Initiate parallel pool
% Prepare the code to run on NeSI in parallel
pc = parcluster('local');
pc.JobStorageLocation = getenv('TMPDIR');
parpool(pc, str2double(getenv('SLURM_CPUS_PER_TASK')));

tic
% Addpath with output files, shapefiles and packages
addpath '/nesi/project/vuw03295/Matlab_files/Ninety_miles_beach'
addpath '/nesi/nobackup/vuw03295/Regional_project/Ninety_miles_beach/outputs_foreward_1'
addpath '/nesi/nobackup/vuw03295/Map_files'

% Add the shapefiles for the land, habitats and release areas
Mussel_habitat = shaperead('mussel_habitat_rough_005deg_moana_mask.shp');
Release_Areas = shaperead('mussel_habitat_rough_005deg_moana_mask.shp');
NZ = shaperead('New_Zealand_High_resolution_coast_for_plot.shp');


%% Loop through the output files

traj_output=dir('/nesi/nobackup/vuw03295/Regional_project/Ninety_miles_beach/outputs_foreward_1/ninety_miles_beach_foreward_1_*');
disp_release = cell(1, length(traj_output)); % Pre-allocate variable for release areas connected
frequency_habitat = cell(1, length(traj_output)); % Pre-allocate variable
frequency_release = cell(1, length(traj_output)); % Pre-allocate variable
disp_settlement = cell(1, length(traj_output)); % Pre-allocate the variable settlement

for fi = 1:length(traj_output)
    
    %% Load and sort variables
    
    % Get fill value
    ncid = netcdf.open(traj_output(fi).name,'NOWRITE');
    ncdisp(traj_output(fi).name) % Check the attributes of the netcdf file
    varid = netcdf.inqVarID(ncid,'lat');
    [noFillMode,fillValue] = netcdf.inqVarFill(ncid,varid);
    
    % Load files
    disp_lat=ncread(traj_output(fi).name,'lat');disp_lat(disp_lat == fillValue) = NaN;
    disp_lon=ncread(traj_output(fi).name,'lon');disp_lon(disp_lon == fillValue) = NaN;
    %disp_depth=ncread(traj_output(fi).name,'z');disp_depth(disp_depth == fillValue) = NaN;
    disp_age=ncread(traj_output(fi).name,'age_seconds');disp_age(disp_age == fillValue) = NaN;
    
    %% Remove values for age > 35 days => 5 weeks of dispersal for blue mussels
    
    max_PLD = 35 * 24 * 3600; % maximum PLD in seconds
    
    for i = 1:size(disp_age, 1)
        for j = 1:size(disp_age, 2)
            if disp_age(i, j) < max_PLD
                disp_lat(i, j) = NaN;
                disp_lon(i, j) = NaN;
            end
        end
    end

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
    
    
    %% Settlement of larvae
    
    % Settlement time-window: 3 to 5 weeks (Norrie et al., 2020, Jeffs et al., 1999)
    Competency_start = 21 * 24* 3600; % Beginning of the settlement in seconds
    
    % Find larvae in habitat at settlement stage
    in_habitat = zeros(size(disp_lat,1),size(disp_lat,2));
    for i = 1:size(disp_lat,2)
        Competent = find(disp_age(:,i) > Competency_start);
        End_life = find(disp_age(:,i) >= max_PLD);
        parfor(j = Competent(1):End_life(1))
        %for j = Competent(1):End_life(1)
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
    
    % Count frequency of settlement per site
    frequency_habitat{fi} = tabulate(in_habitat(:));
    disp_settlement{fi} = in_habitat;
    
    %% Identify starting area that give settlers
    
    Release_poly = zeros(1, size(disp_age,2));
    
    for i = 1:size(disp_age,2)
        if max(in_habitat(:,i)) > 0
            B = find(disp_age(:,i) == min(disp_age(:,i)));
            for z = 1:size(Release_Areas,1)
                if inpolygon(disp_lat(B,i), disp_lon(B,i), [Release_Areas(z).Y], [Release_Areas(z).X]) == 1
                    Release_poly(i) = Release_Areas(z).Area_code;
                end
            end
        end
    end
    
    % Count frequency of settlement per site
    frequency_release{fi} = tabulate(Release_poly);
    disp_release{fi} = Release_poly;   

end
toc

save('frequency_release.mat','frequency_release');
save('disp_release.mat','disp_release','-v7.3');
save('frequency_settlement.mat','frequency_habitat');
save('disp_settlement.mat','disp_settlement','-v7.3');

clear disp_age
clear disp_depth
clear disp_lat
clear disp_lon

