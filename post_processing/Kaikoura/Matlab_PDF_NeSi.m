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
Paua_habitat_for_plot =  shaperead('Kaikoura_paua.shp');
addpath('/nesi/project/vuw03295/Matlab_files/Paua_national/') % Directory with Matlab files for run
% Need to import a new shapefile for the plots
addpath('/nesi/nobackup/vuw03295/Map_files/')
NZ = shaperead('New_Zealand_High_resolution_coast_for_plot.shp');

% Import output files in a loop to work on all years => loop over
addpath('/nesi/nobackup/vuw03295/Regional_project/Kaikoura')
cd '/nesi/nobackup/vuw03295/Regional_project/Kaikoura'
output_files=dir('./20*');
output_files_names={output_files.name};
output_files_names=natsortfiles(output_files_names);


%% Now compute the PDf for each year

for fl = 1:length(output_files_names)
    
    % Move to the output files location
    evalc(['cd ' output_files_names{fl};]);
    
    % Pre-load the trajectory and connectivity files
    traj_output = dir('Moana_backbone_traj_*.nc');
    
    % Prepare heatmap
    resolution_cell_PDF = 1000;
    heatmap = zeros(resolution_cell_PDF-1, resolution_cell_PDF-1);
    
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
        %disp_age=ncread(traj_output(fi).name,'age_seconds');disp_age(disp_age == fillValue) = NaN;
        
        % Compute and save heatmap
        latitudes_nest = linspace(-44.5, -40.5, resolution_cell_PDF);
        longitudes_nest = linspace(172.5, 177.5, resolution_cell_PDF);
        N = histcounts2(disp_lat(:,:), disp_lon(:,:), latitudes_nest, longitudes_nest);
        heatmap = heatmap+N;
        
    end
    
    % Plot heatmap
    figure()
    imagesc(longitudes_nest, latitudes_nest, log10(heatmap/length(disp_lon))); %without smoothing
    axis equal
    set(gca, 'XLim', [172.5, 177.5], 'YLim', [-44.5, -40.5] , 'YDir', 'normal');
    colorbar
    caxis([-5 -2])
    hold on
    % Plot coastline
    mapshow(NZ)
    hold on
    geoshow(Paua_habitat_for_plot)
    title(sprintf('Trajectories density - Kaikoura Paua %d',2007+fl))
    c = colorbar;
    c.Label.String = 'Density of trajectories (log10 scale)'; %unit in position recorded per ~0.011 degree lon
    PDF_figure = ['PDF_Kaikoura_' num2str(fl+2007)];
    savefig(gcf,PDF_figure,'compact')
    saveas(gcf,PDF_figure,'epsc');
    saveas(gcf,PDF_figure,'jpg');
    
    cd ..
end