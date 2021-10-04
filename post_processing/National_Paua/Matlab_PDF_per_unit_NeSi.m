%% Analysis of the population connectivity of Haliotis iris (Paua) in New Zealand: Part 3
% Description of the project:
%
% Description of Matlab code: Produce and plot the trajectories PDF for each individual
% PAU management area.
% Author: RChaput - 24/08/2021

clear
close

% Import settlement habitat
addpath('/nesi/project/vuw03295/National_projects/Paua/input_files')
addpath('/nesi/project/vuw03295/Matlab_files/Paua_national/')
Paua_habitat = shaperead('National_distribution_paua.shp'); % Habitat divided in 465 polygons: Need SRC 4985 in QGIS to export Chatham with negative values


% Import output files in a loop to work on all years => loop over
addpath('/nesi/nobackup/vuw03295/Regional_project/Kaikoura')
cd '/nesi/nobackup/vuw03295/Regional_project/Kaikoura'
output_files=dir('./201*');
output_files_names={output_files.name};
output_files_names=natsortfiles(output_files_names);

%% Now compute the PDf for each year and each PAU area

% Need to import a new shapefile for the plots
Paua_habitat_for_plot =  shaperead('National_distribution_paua_sorted.shp');
NZ = shaperead('New_Zealand_High_resolution_coast_for_plot.shp');

for fl = 1:length(output_files_names)
    
    % Move to the output files location
    evalc(['cd ' output_files_names{fl};]);
    
    % Pre-load the trajectory and connectivity files
    traj_output = dir('Moana_backbone_traj_*.nc');
    
    % Load the connectivity file to ID the trajectories
    load('traj_con_file.mat');
    % Prepare heatmap
    heatmap = cell(1, 9);
    resolution_cell_PDF = 2000;
    heatmap(1,:) = {zeros(resolution_cell_PDF-1, resolution_cell_PDF-1)};
    
    for fj = 1:size(list_PAU,2)
        % Select traj within the PAU of interest
        % Choose PAU area to plot for the dispersal kernels
        target_reefs = list_PAU{fj};
        
        for fi = 1:length(traj_output)
            
            % Load the traj files
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
            
            
            % Select traj within the PAU area
            for j = 1:length(traj_con_file{fi})
                if sum(traj_con_file{fi}(j,1) == target_reefs)==0 %
                    disp_lat(:,j)=NaN; disp_lon(:,j) = NaN;
                end
            end
            
            % Readjust the longitudes for Chatham Islands
            for i=1:size(disp_lon,1) % Need to solve this 180 to -180 problem
                for ii=1:size(disp_lon,2)
                    if disp_lon(i,ii)<0; disp_lon(i,ii)=360+disp_lon(i,ii); end
                end
            end
            % Compute and save heatmap  =>[163, 180, -52, -31]=all domain
            latitudes_nest = linspace(-52, -31, resolution_cell_PDF);
            longitudes_nest = linspace(163, 190, resolution_cell_PDF);
            N = histcounts2(disp_lat(:,:), disp_lon(:,:), latitudes_nest, longitudes_nest);
            heatmap{fj} = heatmap{fj}+N;
            
        end
        
        % Plot heatmap for single PAU management area
        figure()
        imagesc(longitudes_nest, latitudes_nest, log10(heatmap{fj}/length(disp_lon))); %without smoothing
        axis equal
        set(gca, 'XLim', list_lon_lim{fj}, 'YLim', list_lat_lim{fj} , 'YDir', 'normal');
        colorbar
        caxis([-5 -3])
        hold on
        % Plot coastline
        mapshow(NZ)
        hold on
        symspec = makesymbolspec('Polygon', ...
            {'Pau_mgmt','PAU1','FaceColor',[0.9 1 0]}, ...
            {'Pau_mgmt','PAU2','FaceColor',[0 0.9 0]}, ...
            {'Pau_mgmt','PAU3','FaceColor',[0 0.8 1]}, ...
            {'Pau_mgmt','PAU4','FaceColor',[1 0 1]}, ...
            {'Pau_mgmt','PAU5A','FaceColor',[1 1 1]}, ...
            {'Pau_mgmt','PAU5B','FaceColor',[0.5 0.5 0.5]}, ...
            {'Pau_mgmt','PAU5D','FaceColor',[0.1 0.6 0.4]}, ...
            {'Pau_mgmt','PAU6','FaceColor',[0.8 0.5 0.4]}, ...
            {'Pau_mgmt','PAU7','FaceColor',[0.7 0.1 0.9]});
        geoshow(Paua_habitat_for_plot, 'SymbolSpec',symspec)
        title(sprintf('Trajectories density - Paua Mgmt %d',fj))
        c = colorbar;
        c.Label.String = 'Density of trajectories (log10 scale)'; %unit in position recorded per ~0.011 degree lon
        PDF_figure = ['PDF_PAU_' num2str(fj)];
        savefig(gcf,PDF_figure,'compact')
        saveas(gcf,PDF_figure,'epsc');
        saveas(gcf,PDF_figure,'jpg');
        
    end
    cd ..
end