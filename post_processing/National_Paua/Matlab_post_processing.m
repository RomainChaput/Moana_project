%% Analysis of the population connectivity of Haliotis iris (Paua) in New Zealand: Part 2
% Description of the project:
%
% Description of Matlab code: Creates the Connectivity Matrix and the
% Dispersal Kernel for each PAU management area.
% Author: RChaput - 24/08/2021

clear
close

% Import settlement habitat
addpath('.\National_distribution_map')
addpath('.\Opendrift_scripts')
Paua_habitat = shaperead('National_distribution_paua_sorted.shp'); % Habitat divided in 465 polygons: Need SRC 4985 in QGIS to export Chatham with negative values

% Load the centroids from the release file
opts = delimitedTextImportOptions("NumVariables", 3);
% Specify range and delimiter
opts.DataLines = [1, Inf];
opts.Delimiter = "\t";
% Specify column names and types
opts.VariableNames = ["Longitude", "Latitude", "ID"];
opts.VariableTypes = ["double", "double", "double"];
opts = setvaropts(opts, 3, "EmptyFieldRule", "auto");
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";
% Import the data
centroids = readtable("Release_centroid_nat_dist_paua.xyz", opts);
clear opts
% Transform the LL matrix for future work
centroids = table2array(centroids); % f_latlong only works with array
centroids(:,4) = [Paua_habitat(:).ordre];
[Xsorted,I] = sort(centroids(:,4));
centroids = centroids(I,:);

% Import output files in a loop to work on all years
addpath('.\outputs_con_files')
cd '.\outputs_con_files'
output_files=dir('./19*');
output_files_names={output_files.name};
output_files_names=natsortfiles(output_files_names);

%% Create connectivity matrix

% Start from here with pre-compile connectivity files
for fl = 1:length(output_files_names)
    
    % Move to the output files location
    evalc(['cd ' output_files_names{fl};]);
    load Con_file_sorted
    Con_file(any(isnan(Con_file), 2), :) = [];
    
    % Build the sparse matrix
    n = length(Paua_habitat); % Number of polygons
    mtx =[];
    for i = 1:n
        %arrival polygon in column #2
        pol=Con_file(Con_file(:,2)==i,:); % Con_file(:,2)=> lat, Con_file(:,4)=>lon
        if isempty(pol)==1
            a = ones(n,1)*0;
        end
        for j=1:n
            %source polygon in column #1
            source=pol(pol(:,1)==j,1); % Con_file(:,1)=> lat, Con_file(:,3)=>lon
            a (j,:) =length(source);
        end
        mtx = [mtx a];
    end
    
    % Save the connectivity matrix file
    save Mtx_Paua_NZ_lat mtx
    
    %% Plot the connectivity matrix
    load Mtx_Paua_NZ_lat % lat or lon
    nb_larvae_per_poly = 26700; % 300 larvae released per day per polygon (89 days release)
    [xi,yi]=meshgrid(1:1:n, 1:1:n);
    zi = mtx/nb_larvae_per_poly;
    
    figure()
    zi(isnan(zi))=0;
    mask = zi ~= 0;
    im=pcolor(xi,yi,log10(zi*100));
    %im = imagesc(log10(zi_subset));
    set(im,'alphadata',mask);
    axis square
    shading flat
    hold on
    axis([0 n 0 n])
    ylabel('Source Node');
    xlabel('Receiving Node');
    colorbar
    wg = jet;
    wg(1,:) = [1 1 1];
    colormap(wg)
    title(sprintf('Connectivity Matrix %d - New Zealand Paua',1993+fl))
    c = colorbar;
    c.Label.String = 'Percentage of settlement (log10 scale)';
    
    % Now create connectivity matrix with flotting average (easier to see)
    % Average over multiple polygones zi to facilitate plot reading for large matrix
    reduc = 4; % Factor of reduction of the matrix
    zi = mtx/nb_larvae_per_poly;
    zi_subset = zeros(length(zi)/reduc,length(zi)/reduc);
    vector_rows = (1:length(zi)/reduc);
    for i = 1:length(vector_rows)
        vector_rows(i) = reduc;
    end
    vector_columns = vector_rows;
    C = mat2cell(zi, vector_rows, vector_columns);
    for cell_row = 1:length(zi)/reduc
        for cell_column = 1:length(zi)/reduc
            % Change the number of mean functions as needed: nb(mean)=reduc
            zi_subset(cell_row,cell_column)=(mean(mean(mean(mean(C{cell_row,cell_column})))));
        end
    end
    
    [xi,yi]=meshgrid(1:reduc:n, 1:reduc:n);
    zi_subset(isnan(zi_subset))=0;
    mask = zi_subset ~= 0;
    figure()
    im=pcolor(xi,yi,log10(zi_subset*100));
    %im = imagesc(log10(zi_subset));
    set(im,'alphadata',mask);
    axis square
    shading flat
    hold on
    axis([0 n 0 n])
    ylabel('Source Node, sorted by Latitudes');
    xlabel('Receiving Node, sorted by Latitudes');
    colorbar
    wg = jet;
    wg(1,:) = [1 1 1];
    colormap(wg)
    title(sprintf('Connectivity Matrix %d - New Zealand Paua',1993+fl))
    c = colorbar;
    c.Label.String = 'Percentage of settlement (log10 scale)';
    % Save the figure
    con_fig = ['Connectivity_matrix_' num2str(1993 + fl)];
    savefig(gcf, con_fig, 'compact')
    saveas(gcf, con_fig, 'epsc');
    saveas(gcf, con_fig, 'jpg');
    
    
    % Rearrange the matrix to group per PAU mgmt area:
    PAU_1 =sort([Paua_habitat([Paua_habitat.ID_mgmt]==1).ordre]); % ordre_lon (for lon sorting) or ordre (for lat)
    PAU_2 =sort([Paua_habitat([Paua_habitat.ID_mgmt]==2).ordre]);
    PAU_3 =sort([Paua_habitat([Paua_habitat.ID_mgmt]==3).ordre]);
    PAU_4 =sort([Paua_habitat([Paua_habitat.ID_mgmt]==4).ordre]);
    PAU_5A =sort([Paua_habitat([Paua_habitat.ID_mgmt]==5).ordre]);
    PAU_5B =sort([Paua_habitat([Paua_habitat.ID_mgmt]==6).ordre]);
    PAU_5D =sort([Paua_habitat([Paua_habitat.ID_mgmt]==7).ordre]);
    PAU_6 =sort([Paua_habitat([Paua_habitat.ID_mgmt]==8).ordre]);
    PAU_7 =sort([Paua_habitat([Paua_habitat.ID_mgmt]==9).ordre]);
    new_mat = [PAU_5B, PAU_5A, PAU_5D, PAU_3, PAU_6, PAU_7, PAU_2, PAU_1, PAU_4];
    
    % plot the matrix for subset and polygons sorted per PAU mgmt units
    zi=mtx((new_mat),(new_mat));
    zi = zi/nb_larvae_per_poly;
    % Average over multiple polygones zi to facilitate plot reading for large matrix
    reduc = 4; % Factor of reduction of the matrix
    zi_subset = zeros(length(zi)/reduc,length(zi)/reduc);
    vector_rows = (1:length(zi)/reduc);
    for i = 1:length(vector_rows)
        vector_rows(i) = reduc;
    end
    vector_columns = vector_rows;
    C = mat2cell(zi, vector_rows, vector_columns);
    for cell_row = 1:length(zi)/reduc
        for cell_column = 1:length(zi)/reduc
            % Change the number of mean functions as needed: nb(mean)=reduc
            zi_subset(cell_row,cell_column)=(mean(mean(mean(mean(C{cell_row,cell_column})))));
        end
    end
    [xi,yi]=meshgrid(1:reduc:n, 1:reduc:n);
    zi_subset(isnan(zi_subset))=0;
    mask = zi_subset ~= 0;
    figure()
    im=pcolor(xi,yi,log10(zi_subset*100));
    %im = imagesc(log10(zi_subset));
    set(im,'alphadata',mask);
    axis square
    shading flat
    hold on
    axis([0 n 0 n])
    ylabel('Source Node');
    xlabel('Receiving Node');
    colorbar
    wg = jet;
    wg(1,:) = [1 1 1];
    colormap(wg)
    title(sprintf('Connectivity Matrix %d - New Zealand PAU areas',1993+fl))
    c = colorbar;
    %c.Label.String = 'Percentage of settlement (log10 scale)';
    c.Label.String = 'Percentage of settlement (log10 scale)';
    set(gca,'YTickLabel',{'PAU5B', 'PAU5A', 'PAU5D', 'PAU3', 'PAU6', 'PAU7', 'PAU2', 'PAU1', 'PAU4'}, ...
        'YTick',[1+28/2, 1+(37/2+28), 1+(42/2+65), 1+(28/2+107), 1+(26/2+135), 1+(36/2+161),1+(99/2+197),1+(109/2+296),1+(59/2+405)],'fontsize',10);
    set(gca,'XTickLabel',{'5B', '5A', '5D', '3', '6', '7', '2', '1', '4'}, ...
        'XTick',[1+28/2, 1+(37/2+28), 1+(42/2+65), 1+(28/2+107), 1+(26/2+135), 1+(36/2+161),1+(99/2+197),1+(109/2+296),1+(59/2+405)],'fontsize',10);
    xlines= [1+28, 1+65, 1+107, 1+135, 1+161, 1+197,1+296,1+405];
    x_grid=[xlines; xlines];
    plot(x_grid, '-', 'color', [1,1,1])
    gridindex = [1+28, 1+65, 1+107, 1+135, 1+161, 1+197,1+296,1+405];
    makegrid(gridindex,n)
    %Save the figure
    con_fig = ['PAU_Connectivity_matrix_' num2str(1993 + fl)];
    savefig(gcf, con_fig, 'compact')
    saveas(gcf, con_fig, 'epsc');
    saveas(gcf, con_fig, 'jpg');
    
    
    %% Work on the dispersal kernels
    
    % Selection of the areas to work with to plot the dispersal kernels
    PAU_1 =[Paua_habitat([Paua_habitat.ID_mgmt]==1).ordre]; % ordre_lon (for lon sorting) or ordre (for lat)
    PAU_2 =[Paua_habitat([Paua_habitat.ID_mgmt]==2).ordre];
    PAU_3 =[Paua_habitat([Paua_habitat.ID_mgmt]==3).ordre];
    PAU_4 =[Paua_habitat([Paua_habitat.ID_mgmt]==4).ordre];
    PAU_5A =[Paua_habitat([Paua_habitat.ID_mgmt]==5).ordre];
    PAU_5B =[Paua_habitat([Paua_habitat.ID_mgmt]==6).ordre];
    PAU_5D =[Paua_habitat([Paua_habitat.ID_mgmt]==7).ordre];
    PAU_6 =[Paua_habitat([Paua_habitat.ID_mgmt]==8).ordre];
    PAU_7 =[Paua_habitat([Paua_habitat.ID_mgmt]==9).ordre];
    
    list_PAU = {PAU_1, PAU_2, PAU_3, PAU_4, PAU_5A, PAU_5B, PAU_5D, PAU_6, PAU_7};
    load Con_file_sorted
    Con_file(any(isnan(Con_file), 2), :) = [];
    
    for fj = 1:size(list_PAU,2)
        % Choose PAU area to plot for the dispersal kernels
        target_reefs = list_PAU{fj};
        
        poly_area = zeros(size(Con_file,1),size(Con_file,2));
        for i = 1:length(Con_file)
            for j = 1:length(target_reefs)
                if Con_file(i,1)== target_reefs(j); poly_area(i,:) = Con_file(i,:);
                end
            end
        end
        poly_area = poly_area(any(poly_area,2),:);
        
        % Build the sparse matrix
        mtx =[];
        n = size(Paua_habitat,1);
        for i = 1:n %for entire matrix
            %arrival polygon in column #2
            pol=poly_area(poly_area(:,2)==i,:);
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
        
        % Save the sparse raw matrix for the year
        matrice = ['Mtx_' num2str(fj) '.mat'];
        save(matrice, 'mtx');
        
        % Create the distance matrice based on the centroids of the
        % polygons (release file)
        D = f_latlong(centroids);
        matrice2 = ['DistMat_' num2str(fj) '.mat'];
        save(matrice2, 'D');
        
        % Calculate the dispersal Kernels and save the output figure
        d = 10; % bin distance in kilometers for the dispersal kernel plot
        n = size(centroids,1);% number of polygons
        window = 2; % smooth window
        kernels = dispK(D,mtx/nb_larvae_per_poly,n,d,window); % First line: mean settling probability per d kilometers. Second line: std
        matrice3 = ['DispKernels_' num2str(fj) '.mat'];
        save(matrice3, 'kernels');
        dispk_figure = ['Dispersal_kernel_' num2str(fj)];
        savefig(gcf,dispk_figure,'compact')
        saveas(gcf,dispk_figure,'epsc');
        saveas(gcf,dispk_figure,'jpg');
        %close
        
    end
    close all
    cd ..
end
