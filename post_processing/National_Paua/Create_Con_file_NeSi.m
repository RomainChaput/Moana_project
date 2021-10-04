%% Analysis of the population connectivity of Haliotis iris (Paua) in New Zealand: Part 1
% Description of the project: 25 years connectivity of Haliotis iris in New
% Zealand. We will plot the connectivity matrices for each management area,
% the PDF of trajectories, and the dispersal kernels.
%
% Description of Matlab code: Post-processing script for connectivity files
% produced by OpenDrift new habitat module. Run on NeSi, directly on the
% nobackup directory
% Author: RChaput - 09/09/2021

clear
close
addpath('/nesi/project/vuw03295/Matlab_files/Paua_national')

% Import settlement habitat
addpath('/nesi/project/vuw03295/National_projects/Paua/input_files')
Paua_habitat = shaperead('National_distribution_paua.shp'); % Habitat divided in 465 polygons

% Import output files in a loop to work on all years
addpath('/nesi/nobackup/vuw03295/National_project/Paua')
cd '/nesi/nobackup/vuw03295/National_project/Paua'
output_files=dir('./2014*');
output_files_names={output_files.name};
output_files_names=natsortfiles(output_files_names);

for fl = 1%length(output_files_names)
    
    % Move to the output files location
    evalc(['cd ' output_files_names{fl};]);
    
    % Concatenate the output files for 1 year
    con_outputs = dir('./Con_file_Moana*');
    % Create empty connectivity file to fill in a loop
    Total_Con_file = [NaN, NaN, NaN, NaN];
    Con_file = [NaN, NaN, NaN, NaN];
    traj_con_file = {};
    release_habitat_lat_sorted = [];
    settle_habitat_lat_sorted = [];
    release_habitat_lon_sorted = [];
    settle_habitat_lon_sorted = [];
    
    for nb = 1:size(con_outputs, 1)
        
        % load the output file with the release and settlement positions
        con_file = load(con_outputs(nb).name);
        
        for i = 1:size(con_file, 1)
            % Sort release and settlement position in lat/lon: 1st and 2nd
            % column: release, 3rd and 4th: settlement
            lon_settle = con_file(i, 3);
            lat_settle = con_file(i, 4);
            lon_release = con_file(i, 1);
            lat_release = con_file(i, 2);
            
            if lon_settle > 180
                lon_settle = -360+lon_settle;
            end
            if lon_release > 180
                lon_release = -360+lon_release;
            end
            
            for z = 1:size(Paua_habitat,1)
                % Find the ID of the release polygon
                if inpolygon(lat_release, lon_release, [Paua_habitat(z).Y], [Paua_habitat(z).X]) == 1
                    release_habitat_lat_sorted = [Paua_habitat(z).ordre];
                    release_habitat_lon_sorted = [Paua_habitat(z).ordre_lon];
                end
                % Find the ID of the settlement polygon
                if inpolygon(lat_settle, lon_settle, [Paua_habitat(z).Y], [Paua_habitat(z).X])== 1
                    settle_habitat_lat_sorted = [Paua_habitat(z).ordre];
                    settle_habitat_lon_sorted = [Paua_habitat(z).ordre_lon];
                end
            end
            if isempty(release_habitat_lat_sorted); release_habitat_lat_sorted = NaN; end
            if isempty(settle_habitat_lat_sorted); settle_habitat_lat_sorted = NaN; end
            if isempty(release_habitat_lon_sorted); release_habitat_lon_sorted = NaN; end
            if isempty(settle_habitat_lon_sorted); settle_habitat_lon_sorted = NaN; end
            % Update the connectivity file with the ID of polygon
            Con_file(i, 1) = release_habitat_lat_sorted;
            Con_file(i, 2) = settle_habitat_lat_sorted;
            Con_file(i, 3) = release_habitat_lon_sorted;
            Con_file(i, 4) = settle_habitat_lon_sorted;
            % Update the connectivity file for the PDF
            traj_con_file{nb}(i, 1) = release_habitat_lat_sorted;
            traj_con_file{nb}(i, 2) = settle_habitat_lat_sorted;
            traj_con_file{nb}(i, 3) = release_habitat_lon_sorted;
            traj_con_file{nb}(i, 4) = settle_habitat_lon_sorted;
            release_habitat_lat_sorted = NaN;
            settle_habitat_lat_sorted = NaN;
            release_habitat_lon_sorted = NaN;
            settle_habitat_lon_sorted = NaN;
        end
        
        Total_Con_file = [Total_Con_file; Con_file];
        
    end
    
    clear settle_habitat
    clear release_habitat
    clear con_file
    clear Con_file
    
    Con_file = Total_Con_file; clear Total_Con_file
    save Con_file_sorted Con_file
    save traj_con_file traj_con_file; clear traj_con_file
    
    cd ..
end