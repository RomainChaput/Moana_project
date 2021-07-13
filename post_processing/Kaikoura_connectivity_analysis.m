%% Analysis of the population connectivity of Haliotis iris (Paua) at Kaikoura
% Description of the project: 
%
% Description of Matlab code: Post-processing script for connectivity files produced by OpenDrift new habitat module
% Author: RChaput - 28/06/2021

clear
close

% Import settlement habitat
addpath('.\Matlab_post_processing\habitat')
addpath('.\Matlab_post_processing')
Paua_habitat = shaperead('Kaikoura_paua.shp'); % Habitat divided in 74 polygons

% Import output files in a loop to work on all years
addpath('.\Results_con_files')
cd '.\Results_con_files'
output_files=dir('./20*');
output_files_names={output_files.name};
output_files_names=natsortfiles(output_files_names);

%% Loop to find release and settlement points
for fl = 1:length(output_files_names)
    
    % Move to the output files location
    evalc(['cd ' output_files_names{fl};]);
    
    % Concatenate the output files for 1 year
    con_outputs = dir('./Con_file_*');
    % Create empty connectivity file to fill in a loop
    Total_Con_file = [NaN, NaN];
    Con_file = [NaN, NaN];
    release_habitat = [];
    settle_habitat = [];
    
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
            
            for z = 1:size(Paua_habitat,1)
                % Find the ID of the release polygon
                if inpolygon(lat_release, lon_release, [Paua_habitat(z).Y], [Paua_habitat(z).X]) == 1
                    release_habitat = [Paua_habitat(z).ID_habitat];
                end
                % Find the ID of the settlement polygon
                if inpolygon(lat_settle, lon_settle, [Paua_habitat(z).Y], [Paua_habitat(z).X])== 1
                    settle_habitat = [Paua_habitat(z).ID_habitat];
                end
            end
            if isempty(release_habitat); release_habitat = NaN; end
            if isempty(settle_habitat); settle_habitat = NaN; end
            % Update the connectivity file with the ID of polygon
            Con_file(i, 1) = release_habitat;
            Con_file(i, 2) = settle_habitat;
            release_habitat = NaN;
            settle_habitat = NaN;
        end
        
        Total_Con_file = [Total_Con_file; Con_file];
        
    end
    
    clear settle_habitat
    clear release_habitat
    clear con_file
    clear Con_file
    
    Con_file = Total_Con_file; clear Total_Con_file
    save Con_file Con_file
    
    cd ..
end
   
%% Create connectivity matrix

% Start from here with pre-compile connectivity files
for fl = 1:length(output_files_names)
    
    % Move to the output files location
    evalc(['cd ' output_files_names{fl};]);
    load Con_file
    % Build the sparse matrix
    n = length(Paua_habitat)+1; % Number of polygons
    mtx =[];
    for i = 1:n
        %arrival polygon in column #2
        pol=Con_file(Con_file(:,2)==i,:);
        if isempty(pol)==1
            a = ones(n,1)*0;
        end
        for j=1:n
            %source polygon in column #1
            source=pol(pol(:,1)==j,1);
            a (j,:) =length(source);
        end
        mtx = [mtx a];
    end
    
    % Save the connectivity matrix file
    save Mtx_Kaikoura mtx
    
    %plot the connectivity matrix
    nb_larvae_per_poly = 26700; % 300 larvae released per day per polygon (89 days release)
    [xi,yi]=meshgrid(1:1:n, 1:1:n);
    zi = mtx/nb_larvae_per_poly;
    figure()
    %h1 = axes;
    pcolor(xi,yi,log10(zi));
    hold on
    axis square
    ylabel('Source Habitat');
    xlabel('Receiving Habitat');
    colorbar
    %title('Connectivity Matrix - Kaikoura Paua')
    title(sprintf('Connectivity Matrix %d - Kaikoura Paua',2007+fl))
    c = colorbar;
    caxis([-5 -0.5])
    c.Label.String = 'Percentage of settlement (log10 scale)';
    
    set(gca,'YTickLabel',{'Christ', 'Motu', 'GoreB', 'A', 'ConwFlat', 'B', 'KaiS', 'C', 'KaiS', 'KaiN', 'D', 'Halk', 'E', 'F', 'Clar', 'G', 'WCamp', 'H', 'Welli', 'CapeP'}, ...
        'YTick',[1+12/2, 1+(9/2+12), 1+(2/2+21), 1+(1/2+23), 1+(2/2+24), 1+(1/2+26), 1+(1/2+27), 1+(1/2+28), 1+(2/2+29), 1+(1/2+31), 1+(1/2+32), 1+(1/2+33), 1+(1/2+34), 1+(1/2+35), ...
        1+(5/2+36), 1+(2/2+41), 1+(7/2+44), 1+(1/2+51), 1+(8/2+52), 1+(14/2+60)],'fontsize',10);
    set(gca,'XTickLabel',{'Christ', 'Motu', 'GoreB', 'A', 'ConwFlat', 'B', 'KaiS', 'C', 'KaiS', 'KaiN', 'D', 'Halk', 'E', 'F', 'Clar', 'G', 'WCamp', 'H', 'Welli', 'CapeP'}, ...
        'XTick',[1+12/2, 1+(9/2+12), 1+(2/2+21), 1+(1/2+23), 1+(2/2+24), 1+(1/2+26), 1+(1/2+27), 1+(1/2+28), 1+(2/2+29), 1+(1/2+31), 1+(1/2+32), 1+(1/2+33), 1+(1/2+34), 1+(1/2+35), ...
        1+(5/2+36), 1+(2/2+41), 1+(7/2+44), 1+(1/2+51), 1+(8/2+52), 1+(14/2+60)],'fontsize',10);
    xlines= [1+12, 1+21, 1+23, 1+24, 1+26, 1+27, 1+28, 1+29, 1+31, 1+32, 1+33, 1+34, 1+35, 1+36, 1+41, 1+43, 1+51, 1+52, 1+60];
    x_grid=[xlines; xlines];
    plot(x_grid, '-', 'color', [1,1,1])
    gridindex = [1+12, 1+21, 1+23, 1+24, 1+26, 1+27, 1+28, 1+29, 1+31, 1+32, 1+33, 1+34, 1+35, 1+36, 1+41, 1+43, 1+51, 1+52, 1+60];
    makegrid(gridindex,n)

    con_fig = ['Connectivity_matrix_' num2str(2007 + fl)];
    savefig(gcf, con_fig, 'compact')
    saveas(gcf, con_fig, 'epsc');
    saveas(gcf, con_fig, 'jpg');
    
    cd ..
end


%% Create mean population connectivity

% Reload parameters
n = length(Paua_habitat)+1;
% Sum the matrices
MTX = zeros(length(Paua_habitat)+1, length(Paua_habitat)+1);
MTX_3D = zeros(length(Paua_habitat)+1, length(Paua_habitat)+1);

for fl = 1:length(output_files_names)
    
    % Move to the output files location
    evalc(['cd ' output_files_names{fl};]);
    
    load('Mtx_Kaikoura.mat')
    MTX = MTX + mtx;
    MTX_3D(:,:,fl)= mtx;
    
    cd ..
end

% Divide the matrix by the number of years
MTX = MTX/(length(output_files_names));

% Compute the sd of the matrix
sd_MTX = std(MTX_3D, 0, 3);

% Save the files
save Mean_mtx_Kaikoura MTX
save SD_mtx_Kaikoura sd_MTX

%plot the connectivity matrix
nb_larvae_per_poly = 26700; % 300 larvae released per day per polygon (89 days release)
[xi,yi]=meshgrid(1:1:n, 1:1:n);
zi = MTX/nb_larvae_per_poly;
figure()
%h1 = axes;
pcolor(xi,yi,log10(zi));
hold on
axis square
ylabel('Source Habitat');
xlabel('Receiving Habitat');
colorbar
title('10yrs Connectivity Matrix - Kaikoura Paua')
c = colorbar;
caxis([-5 -0.5])
c.Label.String = 'Percentage of settlement (log10 scale)';
set(gca,'YTickLabel',{'Christ', 'Motu', 'GoreB', 'A', 'ConwFlat', 'B', 'KaiS', 'C', 'KaiS', 'KaiN', 'D', 'Halk', 'E', 'F', 'Clar', 'G', 'WCamp', 'H', 'Welli', 'CapeP'}, ...
    'YTick',[1+12/2, 1+(9/2+12), 1+(2/2+21), 1+(1/2+23), 1+(2/2+24), 1+(1/2+26), 1+(1/2+27), 1+(1/2+28), 1+(2/2+29), 1+(1/2+31), 1+(1/2+32), 1+(1/2+33), 1+(1/2+34), 1+(1/2+35), ...
    1+(5/2+36), 1+(2/2+41), 1+(7/2+44), 1+(1/2+51), 1+(8/2+52), 1+(14/2+60)],'fontsize',10);
set(gca,'XTickLabel',{'Christ', 'Motu', 'GoreB', 'A', 'ConwFlat', 'B', 'KaiS', 'C', 'KaiS', 'KaiN', 'D', 'Halk', 'E', 'F', 'Clar', 'G', 'WCamp', 'H', 'Welli', 'CapeP'}, ...
    'XTick',[1+12/2, 1+(9/2+12), 1+(2/2+21), 1+(1/2+23), 1+(2/2+24), 1+(1/2+26), 1+(1/2+27), 1+(1/2+28), 1+(2/2+29), 1+(1/2+31), 1+(1/2+32), 1+(1/2+33), 1+(1/2+34), 1+(1/2+35), ...
    1+(5/2+36), 1+(2/2+41), 1+(7/2+44), 1+(1/2+51), 1+(8/2+52), 1+(14/2+60)],'fontsize',10);
xlines= [1+12, 1+21, 1+23, 1+24, 1+26, 1+27, 1+28, 1+29, 1+31, 1+32, 1+33, 1+34, 1+35, 1+36, 1+41, 1+43, 1+51, 1+52, 1+60];
x_grid=[xlines; xlines];
plot(x_grid, '-', 'color', [1,1,1])
gridindex = [1+12, 1+21, 1+23, 1+24, 1+26, 1+27, 1+28, 1+29, 1+31, 1+32, 1+33, 1+34, 1+35, 1+36, 1+41, 1+43, 1+51, 1+52, 1+60];
makegrid(gridindex,n)

% plot the standard deviation matrix: not super informative
[xi,yi]=meshgrid(1:1:n, 1:1:n);
zi = sd_MTX/nb_larvae_per_poly;
figure()
%h1 = axes;
pcolor(xi,yi,log10(zi));
hold on
axis square
ylabel('Source Habitat');
xlabel('Receiving Habitat');
colorbar
title('SD Connectivity Matrix - Kaikoura Paua')
c = colorbar;
caxis([-5 -0.5])
c.Label.String = 'Percentage of settlement (log10 scale)';
set(gca,'YTickLabel',{'Christ', 'Motu', 'GoreB', 'A', 'ConwFlat', 'B', 'KaiS', 'C', 'KaiS', 'KaiN', 'D', 'Halk', 'E', 'F', 'Clar', 'G', 'WCamp', 'H', 'Welli', 'CapeP'}, ...
    'YTick',[1+12/2, 1+(9/2+12), 1+(2/2+21), 1+(1/2+23), 1+(2/2+24), 1+(1/2+26), 1+(1/2+27), 1+(1/2+28), 1+(2/2+29), 1+(1/2+31), 1+(1/2+32), 1+(1/2+33), 1+(1/2+34), 1+(1/2+35), ...
    1+(5/2+36), 1+(2/2+41), 1+(7/2+44), 1+(1/2+51), 1+(8/2+52), 1+(14/2+60)],'fontsize',10);
set(gca,'XTickLabel',{'Christ', 'Motu', 'GoreB', 'A', 'ConwFlat', 'B', 'KaiS', 'C', 'KaiS', 'KaiN', 'D', 'Halk', 'E', 'F', 'Clar', 'G', 'WCamp', 'H', 'Welli', 'CapeP'}, ...
    'XTick',[1+12/2, 1+(9/2+12), 1+(2/2+21), 1+(1/2+23), 1+(2/2+24), 1+(1/2+26), 1+(1/2+27), 1+(1/2+28), 1+(2/2+29), 1+(1/2+31), 1+(1/2+32), 1+(1/2+33), 1+(1/2+34), 1+(1/2+35), ...
    1+(5/2+36), 1+(2/2+41), 1+(7/2+44), 1+(1/2+51), 1+(8/2+52), 1+(14/2+60)],'fontsize',10);
xlines= [1+12, 1+21, 1+23, 1+24, 1+26, 1+27, 1+28, 1+29, 1+31, 1+32, 1+33, 1+34, 1+35, 1+36, 1+41, 1+43, 1+51, 1+52, 1+60];
x_grid=[xlines; xlines];
plot(x_grid, '-', 'color', [1,1,1])
gridindex = [1+12, 1+21, 1+23, 1+24, 1+26, 1+27, 1+28, 1+29, 1+31, 1+32, 1+33, 1+34, 1+35, 1+36, 1+41, 1+43, 1+51, 1+52, 1+60];
makegrid(gridindex,n)

%% Compute the anomalies for each year and mean anomalies

%MTX = Mean connectivity
mtx_anomalies = zeros(length(Paua_habitat)+1, length(Paua_habitat)+1);
MTX_anomalies_3D = zeros(length(Paua_habitat)+1, length(Paua_habitat)+1);
MTX_anomalies = zeros(length(Paua_habitat)+1, length(Paua_habitat)+1);

for fl = 1:length(output_files_names)
    
    % Move to the output files location
    evalc(['cd ' output_files_names{fl};]);
    
    load('Mtx_Kaikoura.mat')
    mtx_anomalies = mtx - MTX;
    MTX_anomalies_3D(:,:,fl)= mtx_anomalies;
    
    save mtx_anomalies mtx_anomalies
    
    %plot the anomaly matrix
    nb_larvae_per_poly = 26700; % 300 larvae released per day per polygon (89 days release)
    [xi,yi]=meshgrid(1:1:n, 1:1:n);
    zi = mtx_anomalies/nb_larvae_per_poly;
    figure()
    hi = pcolor(xi,yi,zi);
    set(hi, 'EdgeColor', 'none');
    hold on
    axis square
    caxis([-0.1 0.1])
    colormap(bluewhitered);
    ylabel('Source Habitat');
    xlabel('Receiving Habitat');
    colorbar
    title(sprintf('Yearly anomalies in connectivity %d - Kaikoura Paua',2007+fl))
    c = colorbar;
    c.Label.String = 'Anomalies of settlement (percentage of release)';
    set(gca,'YTickLabel',{'Christ', 'Motu', 'GoreB', 'A', 'ConwFlat', 'B', 'KaiS', 'C', 'KaiS', 'KaiN', 'D', 'Halk', 'E', 'F', 'Clar', 'G', 'WCamp', 'H', 'Welli', 'CapeP'}, ...
        'YTick',[1+12/2, 1+(9/2+12), 1+(2/2+21), 1+(1/2+23), 1+(2/2+24), 1+(1/2+26), 1+(1/2+27), 1+(1/2+28), 1+(2/2+29), 1+(1/2+31), 1+(1/2+32), 1+(1/2+33), 1+(1/2+34), 1+(1/2+35), ...
        1+(5/2+36), 1+(2/2+41), 1+(7/2+44), 1+(1/2+51), 1+(8/2+52), 1+(14/2+60)],'fontsize',10);
    set(gca,'XTickLabel',{'Christ', 'Motu', 'GoreB', 'A', 'ConwFlat', 'B', 'KaiS', 'C', 'KaiS', 'KaiN', 'D', 'Halk', 'E', 'F', 'Clar', 'G', 'WCamp', 'H', 'Welli', 'CapeP'}, ...
        'XTick',[1+12/2, 1+(9/2+12), 1+(2/2+21), 1+(1/2+23), 1+(2/2+24), 1+(1/2+26), 1+(1/2+27), 1+(1/2+28), 1+(2/2+29), 1+(1/2+31), 1+(1/2+32), 1+(1/2+33), 1+(1/2+34), 1+(1/2+35), ...
        1+(5/2+36), 1+(2/2+41), 1+(7/2+44), 1+(1/2+51), 1+(8/2+52), 1+(14/2+60)],'fontsize',10);
    xlines= [1+12, 1+21, 1+23, 1+24, 1+26, 1+27, 1+28, 1+29, 1+31, 1+32, 1+33, 1+34, 1+35, 1+36, 1+41, 1+43, 1+51, 1+52, 1+60];
    x_grid=[xlines; xlines];
    plot(x_grid, '-', 'color', [1,1,1])
    gridindex = [1+12, 1+21, 1+23, 1+24, 1+26, 1+27, 1+28, 1+29, 1+31, 1+32, 1+33, 1+34, 1+35, 1+36, 1+41, 1+43, 1+51, 1+52, 1+60];
    makegrid(gridindex,n)
    
    con_fig = ['Anomalies_' num2str(2007 + fl)];
    savefig(gcf, con_fig, 'compact')
    saveas(gcf, con_fig, 'epsc');
    saveas(gcf, con_fig, 'jpg');
    
    cd ..
end


%% Compute EOF

[eof_maps, pc] = eof(detrend3(MTX_3D/nb_larvae_per_poly));

%plot the EOF
%nb_larvae_per_poly = 26700; % 300 larvae released per day per polygon (89 days release)
[xi,yi]=meshgrid(1:1:n, 1:1:n);
zi = eof_maps(:,:,1);
figure()
subplot(2,2,1)
hi = pcolor(xi,yi,(zi));
set(hi, 'EdgeColor', 'none');
hold on
axis square
ylabel('Source Habitat');
caxis([-0.1 0.1])
title('EOF 1 - Kaikoura Paua')
colormap(bluewhitered);

zi = eof_maps(:,:,2);
subplot(2,2,2)
%h1 = axes;
hi = pcolor(xi,yi,(zi));
set(hi, 'EdgeColor', 'none');
hold on
axis square
caxis([-0.1 0.1])
title('EOF 2 - Kaikoura Paua')
colormap(bluewhitered);

zi = eof_maps(:,:,3);
subplot(2,2,3)
%h1 = axes;
hi = pcolor(xi,yi,(zi));
set(hi, 'EdgeColor', 'none');
hold on
axis square
ylabel('Source Habitat');
xlabel('Receiving Habitat');
caxis([-0.1 0.1])
title('EOF 3 - Kaikoura Paua')
colormap(bluewhitered);

zi = eof_maps(:,:,4);
subplot(2,2,4)
%h1 = axes;
hi = pcolor(xi,yi,(zi));
set(hi, 'EdgeColor', 'none');
hold on
axis square
xlabel('Receiving Habitat');
caxis([-0.1 0.1])
title('EOF 4 - Kaikoura Paua')
colormap(bluewhitered);


%% Work on the genetic connectivity only

% Genetic sampling sites
Gen_samples = [24, 27, 29, 33, 35, 36, 43, 49];
mtx_gen_3D = zeros(length(Gen_samples)+1, length(Gen_samples)+1, length(output_files_names));

% Find only those points in the connectivity file
for fl = 1:length(output_files_names)
    
    % Move to the output files location
    evalc(['cd ' output_files_names{fl};]);
    load Con_file
    
    % Filter the connectivity file to keep only the connections between
    % sampling sites
    [~,Con_file_gen] = ismember(Con_file, Gen_samples);
    Con_file_gen(Con_file_gen==0) = NaN;
    
    % Build the sparse matrix
    n = length(Gen_samples)+1; % Number of sample sites
    mtx_gen =zeros(n, n);
    for i = 1:n
        for j = 1:n
            mtx_gen(i,j)= sum(Con_file_gen(:,1)==i & Con_file_gen(:,2)==j);
        end
    end
    
    save Mtx_gen mtx_gen
    mtx_gen_3D(:,:,fl)= mtx_gen;
    
    %plot the genetic matrix
    nb_larvae_per_poly = 26700; % 300 larvae released per day per polygon (89 days release)
    [xi,yi]=meshgrid(1:1:n, 1:1:n);
    zi = mtx_gen/nb_larvae_per_poly;
    figure()
    hi = pcolor(xi,yi,zi);
    set(hi, 'EdgeColor', 'none');
    hold on
    axis square
    %caxis([-0.1 0.1])
    colormap(bluewhitered);
    ylabel('Source Habitat');
    xlabel('Receiving Habitat');
    colorbar
    title(sprintf('Connectivity between sampling sites %d - Kaikoura Paua',2007+fl))
    c = colorbar;
    c.Label.String = 'Percentage of settlement';
    set(gca,'YTickLabel',{'BoatH', 'Blocks', 'Katu', 'KaiP', 'HBay', 'Papa', 'Ward', 'CapeC'}, ...
        'YTick',[1+1/2, 1+1/2+1, 1+1/2+2, 1+1/2+3, 1+1/2+4, 1+1/2+5, 1+1/2+6, 1+1/2+7],'fontsize',10);
    set(gca,'XTickLabel',{'BoatH', 'Blocks', 'Katu', 'KaiP', 'HBay', 'Papa', 'Ward', 'CapeC'}, ...
        'XTick',[1+1/2, 1+1/2+1, 1+1/2+2, 1+1/2+3, 1+1/2+4, 1+1/2+5, 1+1/2+6, 1+1/2+7],'fontsize',10);
    
    con_fig = ['Con_mtx_gen_' num2str(2007 + fl)];
    savefig(gcf, con_fig, 'compact')
    saveas(gcf, con_fig, 'epsc');
    saveas(gcf, con_fig, 'jpg');
    
    cd ..
end

%% Plot the 10 years connectivity between the sampling sites

% Compute the mean of the matrix
mean_MTX_gen = mean(mtx_gen_3D, 3);

%plot the mean genetic matrix
nb_larvae_per_poly = 26700; % 300 larvae released per day per polygon (89 days release)
[xi,yi]=meshgrid(1:1:n, 1:1:n);
zi = mean_MTX_gen/nb_larvae_per_poly;
figure()
hi = pcolor(xi,yi,zi);
set(hi, 'EdgeColor', 'none');
hold on
axis square
colormap(bluewhitered);
ylabel('Source Habitat');
xlabel('Receiving Habitat');
colorbar
title('10yrs connectivity between sampling sites - Kaikoura Paua')
c = colorbar;
c.Label.String = 'Percentage of settlement';
set(gca,'YTickLabel',{'BoatH', 'Blocks', 'Katu', 'KaiP', 'HBay', 'Papa', 'Ward', 'CapeC'}, ...
    'YTick',[1+1/2, 1+1/2+1, 1+1/2+2, 1+1/2+3, 1+1/2+4, 1+1/2+5, 1+1/2+6, 1+1/2+7],'fontsize',10);
set(gca,'XTickLabel',{'BoatH', 'Blocks', 'Katu', 'KaiP', 'HBay', 'Papa', 'Ward', 'CapeC'}, ...
    'XTick',[1+1/2, 1+1/2+1, 1+1/2+2, 1+1/2+3, 1+1/2+4, 1+1/2+5, 1+1/2+6, 1+1/2+7],'fontsize',10);









