%% Post-processing of Opendrift outputs for the Ninety Miles Beach - part 3 - subpart 2
%%% The output trajectories have been analyzed on NeSI to extract the
%%% settlement in the habitat and the high trajectories density areas. This
%%% code is to concatenate all of the data in connectivity matrices and
%%% dispersal kernels

%%% The goal is to produce multi-level connectivity matrices that take into
%%% account the possibility of a secondary settlement after drifting with
%%% the algae spats.


clear
close
clc

%cd 'D:\Cours\Post_doc\Moana Project\Regional projects\Ninety miles beach\Open_drift_model\Output_first_foreward\Matlab_post_processing\Results_opendrift_matrices\'


%% Load the output settlement files produced on NeSI  
addpath 'D:\Cours\Post_doc\Moana Project\Regional projects\Ninety miles beach\Open_drift_model\Output_first_foreward\Matlab_post_processing';
addpath 'D:\Cours\Post_doc\Open_Drift_scripts\Matlab_post_processing';
folders = dir('D:\Cours\Post_doc\Moana Project\Regional projects\Ninety miles beach\Open_drift_model\Output_first_foreward\Matlab_post_processing\Results_opendrift_matrices\');

%% Need a loop here to work on all the results folders

load('disp_release_correction.mat')
load('disp_settlement.mat')
Con_file = [];

for fi = 1:12%size(disp_settlement,2)
    con_file = [];
    settlement_mat = disp_settlement{fi};
    release_mat = disp_release{fi};
    
    for i = 1:size(settlement_mat,2)
        settle_ID = settlement_mat(find(settlement_mat(:,i)~=0, 1, 'first'), i);    
        if ~isempty(settle_ID)
            release_ID = release_mat(1, i); 
            % Create a connectivity file with release location, settlement
            % location, and particle ID
            con_file(i,1) = release_ID;
            con_file(i,2) = settle_ID;
            con_file(i,3) = i;
        else
            con_file(i,1) = NaN;
            con_file(i,2) = NaN;
            con_file(i,3) = NaN;
        end
    end
    Con_file = [Con_file; con_file];
    Con_file = Con_file(all(~isnan(Con_file),2),:);
end
%
cd ./new
clear disp_settlement
load('disp_settlement.mat')

for fi = 13:15%size(disp_settlement,2)
    con_file = [];
    settlement_mat = disp_settlement{fi};
    release_mat = disp_release{fi};
    
    for i = 1:size(settlement_mat,2)
        settle_ID = settlement_mat(find(settlement_mat(:,i)~=0, 1, 'first'), i);    
        if ~isempty(settle_ID)
            release_ID = release_mat(1, i); 
            % Create a connectivity file with release location, settlement
            % location, and particle ID
            con_file(i,1) = release_ID;
            con_file(i,2) = settle_ID;
            con_file(i,3) = i;
        else
            con_file(i,1) = NaN;
            con_file(i,2) = NaN;
            con_file(i,3) = NaN;
        end
    end
    Con_file = [Con_file; con_file];
    Con_file = Con_file(all(~isnan(Con_file),2),:);
end
cd ..
%
cd ./new2
clear disp_settlement
load('disp_settlement.mat')

for fi = 16%size(disp_settlement,2)
    con_file = [];
    settlement_mat = disp_settlement{fi};
    release_mat = disp_release{fi};
    
    for i = 1:size(settlement_mat,2)
        settle_ID = settlement_mat(find(settlement_mat(:,i)~=0, 1, 'first'), i);    
        if ~isempty(settle_ID)
            release_ID = release_mat(1, i); 
            % Create a connectivity file with release location, settlement
            % location, and particle ID
            con_file(i,1) = release_ID;
            con_file(i,2) = settle_ID;
            con_file(i,3) = i;
        else
            con_file(i,1) = NaN;
            con_file(i,2) = NaN;
            con_file(i,3) = NaN;
        end
    end
    Con_file = [Con_file; con_file];
    Con_file = Con_file(all(~isnan(Con_file),2),:);
end


%% Create connectivity matrix

n = 50; % Number of polygons
mtx =[];

% Build the sparse matrix
for i = 1:n
    
    %arrival polygon in column #2
    pol=Con_file(find(Con_file(:,2)==i),:);
    if isempty(pol)==1
        a = ones(n,1)*0;
    end
    for j=1:n
        %source polygon in column #1
        source=pol(find(pol(:,1)==j),1);
        a (j,:) =length(source);
    end
    mtx = [mtx a];
    
end

save Mtx mtx

%plot the matrix
nb_larvae_per_poly = 80000; % 80 000 larvae released in total per polygon (1000 larvae per month per release point, 4 months, 20 pts)
[xi,yi]=meshgrid(1:1:n, 1:1:n);
zi = mtx/nb_larvae_per_poly;

figure(1)
h1 = axes;
pcolor(xi,yi,log10(zi));
hold on
axis square
ylabel('Source Habitat');
xlabel('Receiving Habitat');
colorbar
title('Connectivity Matrix - Moana mask passive')
c = colorbar;
caxis([-5 -0.5])
c.Label.String = 'Percentage of settlement (log10 scale)';
set(gca,'YTickLabel',{'SP','NMB 1','NMB2','Ahi','HereN','Aw', 'AwS', 'Hoki', 'Wai', 'Ara', 'CB', 'Wa', 'Waik', 'TA', 'Kai', 'PB', 'Aro', 'Oak'}, 'YTick',[1+8/2, 1+1/2+8, 1+1/2+9, 1+3/2+10, 1+1/2+13, 1+3/2+14, 1+1/2+17, 1+3/2+18, 1+2/2+21, 1+1/2+23, 1+5/2+24, 1+1/2+29, 1+2/2+30, 1+3/2+32, 1+2/2+35, 1+1/2+37, 1+4/2+38, 1+7/2+42],'fontsize',10);
set(gca,'XTickLabel',{'SP','NMB 1','NMB2','Ahi','HereN','Aw', 'AwS', 'Hoki', 'Wai', 'Ara', 'CB', 'Wa', 'Waik', 'TA', 'Kai', 'PB', 'Aro', 'Oak'}, 'XTick',[1+8/2, 1+1/2+8, 1+1/2+9, 1+3/2+10, 1+1/2+13, 1+3/2+14, 1+1/2+17, 1+3/2+18, 1+2/2+21, 1+1/2+23, 1+5/2+24, 1+1/2+29, 1+2/2+30, 1+3/2+32, 1+2/2+35, 1+1/2+37, 1+4/2+38, 1+7/2+42],'fontsize',10);
xlines= [1+8,1+9,1+10,1+13,1+14,1+17,1+18,1+21,1+23,1+24,1+29,1+30,1+32,1+35,1+37,1+38,1+42];
x_grid=[xlines; xlines];
plot(x_grid, '-', 'color', [1,1,1])
gridindex = [1+8,1+9,1+10,1+13,1+14,1+17,1+18,1+21,1+23,1+24,1+29,1+30,1+32,1+35,1+37,1+38,1+42];
makegrid(gridindex,n)
hold on
set(h1, 'Xdir', 'reverse')
set(h1, 'Ydir', 'reverse')

%% Second part of the script

% We need to account for the settlement in high density areas. There,
% mussel larvae are believed to settle on algae. The algae spats then drift
% and some are collected on Ninety Mile Beach

clear disp_settlement
clear Con_file
clear con_file
load('disp_release_correction.mat')
load('disp_high_density_area.mat')

Con_file_hda = [];

for fi = 1:12
    con_file_hda = [];
    settlement_mat = disp_high_density_area{fi};
    release_mat = disp_release{fi};
    
    for i = 1:size(settlement_mat,2)
        settle_ID = settlement_mat(find(settlement_mat(:,i)~=0, 1, 'first'), i);    
        if ~isempty(settle_ID)
            release_ID = release_mat(1, i); 
            % Create a connectivity file with release location, settlement
            % location, and particle ID
            con_file_hda(i,1) = release_ID;
            con_file_hda(i,2) = settle_ID;
            con_file_hda(i,3) = i;
        else
            con_file_hda(i,1) = NaN;
            con_file_hda(i,2) = NaN;
            con_file_hda(i,3) = NaN;
        end
    end
    Con_file_hda = [Con_file_hda; con_file_hda];
    Con_file_hda = Con_file_hda(all(~isnan(Con_file_hda),2),:);
end

cd ./new
clear disp_high_density_area
load('disp_high_density_area.mat')

for fi = 13:15
    con_file_hda = [];
    settlement_mat = disp_high_density_area{fi};
    release_mat = disp_release{fi};
    
    for i = 1:size(settlement_mat,2)
        settle_ID = settlement_mat(find(settlement_mat(:,i)~=0, 1, 'first'), i);    
        if ~isempty(settle_ID)
            release_ID = release_mat(1, i); 
            % Create a connectivity file with release location, settlement
            % location, and particle ID
            con_file_hda(i,1) = release_ID;
            con_file_hda(i,2) = settle_ID;
            con_file_hda(i,3) = i;
        else
            con_file_hda(i,1) = NaN;
            con_file_hda(i,2) = NaN;
            con_file_hda(i,3) = NaN;
        end
    end
    Con_file_hda = [Con_file_hda; con_file_hda];
    Con_file_hda = Con_file_hda(all(~isnan(Con_file_hda),2),:);
end
cd ..
%
cd ./new2
clear disp_high_density_area
load('disp_high_density_area.mat')

for fi = 16
    con_file_hda = [];
    settlement_mat = disp_high_density_area{fi};
    release_mat = disp_release{fi};
    
    for i = 1:size(settlement_mat,2)
        settle_ID = settlement_mat(find(settlement_mat(:,i)~=0, 1, 'first'), i);    
        if ~isempty(settle_ID)
            release_ID = release_mat(1, i); 
            % Create a connectivity file with release location, settlement
            % location, and particle ID
            con_file_hda(i,1) = release_ID;
            con_file_hda(i,2) = settle_ID;
            con_file_hda(i,3) = i;
        else
            con_file_hda(i,1) = NaN;
            con_file_hda(i,2) = NaN;
            con_file_hda(i,3) = NaN;
        end
    end
    Con_file_hda = [Con_file_hda; con_file_hda];
    Con_file_hda = Con_file_hda(all(~isnan(Con_file_hda),2),:);
end
cd ..

% Fix the Con_file_hda: the id of the hda are 2,3,4,5,7,8. We need to
% change this to continuous id from 1 to 6 (easier to plot the connectivity matrix)
id_Con_file_hda = unique(Con_file_hda(:,2))';
new_id = (1:6);
for i = 1:length(Con_file_hda)
    Con_file_hda(i,2) = new_id(Con_file_hda(i,2) == id_Con_file_hda);
end

%% Create small connectivity matrix with habitat sending larvae to high
% density areas
n = 49+1; % Number of polygons
m = 6+1; % Number of high density areas
mtx_hda =[];

% Build the sparse matrix
for i = 1:n
    for j = 1:m
        mtx_hda(i,j)= sum(Con_file_hda(:,1)==i & Con_file_hda(:,2)==j);
    end 
end

save Mtx_hda mtx_hda

%plot the matrix
[xi_hda,yi_hda]=meshgrid(1:1:m, 1:1:n);
zi_hda = mtx_hda/nb_larvae_per_poly;

figure(2)
h1 = axes;
pcolor(xi_hda,yi_hda,log10(zi_hda));
axis equal
xlim tight
hold on
ylabel('Source Habitat');
xlabel('Receiving High Density Area');
colorbar
title('Connectivity Matrix - Moana mask passive 2nd degree')
c = colorbar;
c.Label.String = 'Percentage of settlement (log10 scale)';
caxis([-9 -1])
set(gca,'YTickLabel',{'SP','NMB 1','NMB2','Ahi','HereN','Aw', 'AwS', 'Hoki', 'Wai', 'Ara', 'CB', 'Wa', 'Waik', 'TA', 'Kai', 'PB', 'Aro', 'Oak'}, 'YTick',[1+8/2, 1+1/2+8, 1+1/2+9, 1+3/2+10, 1+1/2+13, 1+3/2+14, 1+1/2+17, 1+3/2+18, 1+2/2+21, 1+1/2+23, 1+5/2+24, 1+1/2+29, 1+2/2+30, 1+3/2+32, 1+2/2+35, 1+1/2+37, 1+4/2+38, 1+7/2+42],'fontsize',10);
set(gca,'XTickLabel',{'A','B','C','D','E','F'}, 'XTick',[1+1/2,2+1/2,3+1/2,4+1/2,5+1/2,6+1/2],'fontsize',10);
xlines= [1+8,1+9,1+10,1+13,1+14,1+17,1+18,1+21,1+23,1+24,1+29,1+30,1+32,1+35,1+37,1+38,1+42];
hold on
%set(h1, 'Xdir', 'reverse')
set(h1, 'Ydir', 'reverse')

clear Con_file_hda
clear con_file_hda
clear disp_release
clear disp_high_density_area

%% Third part of the script

% Now we need to recover the output of the second forward run where we
% released from the high density areas and observed the settlement on the mussel habitats 

cd('D:\Cours\Post_doc\Moana Project\Regional projects\Ninety miles beach\Open_drift_model\Output_second_foreward\Opendrift_mask_passive')

load disp_release.mat
load disp_settlement.mat

Con_file_spat = [];

for fi = 1:length(disp_release)
    con_file_spat = [];
    settlement_mat = disp_settlement{fi};
    release_mat = disp_release{fi};
    
    for i = 1:size(settlement_mat,2)
        settle_ID = settlement_mat(find(settlement_mat(:,i)~=0, 1, 'first'), i);    
        if ~isempty(settle_ID)
            release_ID = release_mat(1, i); 
            % Create a connectivity file with release location, settlement
            % location, and particle ID
            con_file_spat(i,1) = release_ID;
            con_file_spat(i,2) = settle_ID;
            con_file_spat(i,3) = i;
        else
            con_file_spat(i,1) = NaN;
            con_file_spat(i,2) = NaN;
            con_file_spat(i,3) = NaN;
        end
    end
    Con_file_spat = [Con_file_spat; con_file_spat];
    Con_file_spat = Con_file_spat(all(~isnan(Con_file_spat),2),:);
    
end

% Remove null values for release
Con_file_spat(Con_file_spat(:,1) == 0, :)=[];

% Change id of release habitat
id_Con_file_spat = unique(Con_file_spat(:,1))';
new_id = (1:6);
for i = 1:length(Con_file_spat)
    Con_file_spat(i,1) = new_id(Con_file_spat(i,1) == id_Con_file_spat);
end

%% Create small connectivity matrix with high density areas sending larvae to habitat
n = 49+1; % Number of polygons
m = 6+1; % Number of high density areas
mtx_spat =[];

% Build the sparse matrix
for i = 1:m
    for j = 1:n
        mtx_spat(i,j)= sum(Con_file_spat(:,1)==i & Con_file_spat(:,2)==j);
    end 
end

save Mtx_spat mtx_spat

%plot the matrix
nb_larvae_per_hda = 140000; %total numer of larvae released per hda (1000 larvae per point per month, 4 months, 35 pts per hda)
[xi_spat,yi_spat]=meshgrid(1:1:n, 1:1:m);
zi_spat = mtx_spat/nb_larvae_per_hda;

figure(3)
h1 = axes;
pcolor(xi_spat,yi_spat,log10(zi_spat));
axis equal
ylim tight
hold on
ylabel('Source High Density Area');
xlabel('Receiving Habitat');
colorbar
title('Connectivity Matrix - Opendrift mask passive 2nd degree')
c = colorbar;
c.Label.String = 'Percentage of settlement (log10 scale)';
caxis([-9 -1])
set(gca,'XTickLabel',{'SP','NMB 1','NMB2','Ahi','HereN','Aw', 'AwS', 'Hoki', 'Wai', 'Ara', 'CB', 'Wa', 'Waik', 'TA', 'Kai', 'PB', 'Aro', 'Oak'}, 'XTick',[1+8/2, 1+1/2+8, 1+1/2+9, 1+3/2+10, 1+1/2+13, 1+3/2+14, 1+1/2+17, 1+3/2+18, 1+2/2+21, 1+1/2+23, 1+5/2+24, 1+1/2+29, 1+2/2+30, 1+3/2+32, 1+2/2+35, 1+1/2+37, 1+4/2+38, 1+7/2+42],'fontsize',10);
set(gca,'YTickLabel',{'A','B','C','D','E','F'}, 'YTick',[1+1/2,2+1/2,3+1/2,4+1/2,5+1/2,6+1/2],'fontsize',10);
ylines= [1+8,1+9,1+10,1+13,1+14,1+17,1+18,1+21,1+23,1+24,1+29,1+30,1+32,1+35,1+37,1+38,1+42];
hold on
set(h1, 'Ydir', 'reverse')
set(h1, 'Xdir', 'reverse')

%% Fifth part: combine the 2nd degree matrices

cd('D:\Cours\Post_doc\Moana Project\Regional projects\Ninety miles beach\Open_drift_model\Output_first_foreward\matlab_post_processing\Results_opendrift_matrices\Moana_mask_passive')

% Total 2nd degree connectivity
percent_second_deg_mtx = (mtx_hda/nb_larvae_per_poly)*(mtx_spat/nb_larvae_per_hda); % Released numbers are different, so we need to work in percentage of settlement
save Mtx_percent_second_deg percent_second_deg_mtx

%plot the matrix
n = 49+1;
[xi_second_deg,yi_second_deg]=meshgrid(1:1:n, 1:1:n);
zi_second_deg = percent_second_deg_mtx;

figure(4)
h1 = axes;
pcolor(xi_second_deg,yi_second_deg,log10(zi_second_deg));
axis square
hold on
ylabel('Source Habitat');
xlabel('Receiving Habitat (via hda)');
colorbar
title('Connectivity Matrix - Moana mask passive total second degree')
c = colorbar;
c.Label.String = 'Percentage of settlement (log10 scale)';
caxis([-9 -1])
set(gca,'YTickLabel',{'SP','NMB 1','NMB2','Ahi','HereN','Aw', 'AwS', 'Hoki', 'Wai', 'Ara', 'CB', 'Wa', 'Waik', 'TA', 'Kai', 'PB', 'Aro', 'Oak'}, 'YTick',[1+8/2, 1+1/2+8, 1+1/2+9, 1+3/2+10, 1+1/2+13, 1+3/2+14, 1+1/2+17, 1+3/2+18, 1+2/2+21, 1+1/2+23, 1+5/2+24, 1+1/2+29, 1+2/2+30, 1+3/2+32, 1+2/2+35, 1+1/2+37, 1+4/2+38, 1+7/2+42],'fontsize',10);
set(gca,'XTickLabel',{'SP','NMB 1','NMB2','Ahi','HereN','Aw', 'AwS', 'Hoki', 'Wai', 'Ara', 'CB', 'Wa', 'Waik', 'TA', 'Kai', 'PB', 'Aro', 'Oak'}, 'XTick',[1+8/2, 1+1/2+8, 1+1/2+9, 1+3/2+10, 1+1/2+13, 1+3/2+14, 1+1/2+17, 1+3/2+18, 1+2/2+21, 1+1/2+23, 1+5/2+24, 1+1/2+29, 1+2/2+30, 1+3/2+32, 1+2/2+35, 1+1/2+37, 1+4/2+38, 1+7/2+42],'fontsize',10);
xlines= [1+8,1+9,1+10,1+13,1+14,1+17,1+18,1+21,1+23,1+24,1+29,1+30,1+32,1+35,1+37,1+38,1+42];
x_grid=[xlines; xlines];
plot(x_grid, '-', 'color', [1,1,1])
gridindex = [1+8,1+9,1+10,1+13,1+14,1+17,1+18,1+21,1+23,1+24,1+29,1+30,1+32,1+35,1+37,1+38,1+42];
makegrid(gridindex,n)
hold on
set(h1, 'Xdir', 'reverse')
set(h1, 'Ydir', 'reverse')


%% Final part: Combine the first and second degree connectivity matrices

total_mtx = mtx/nb_larvae_per_poly + percent_second_deg_mtx;
save Mtx_total total_mtx

%plot the matrix
n = 49+1;
[xi,yi]=meshgrid(1:1:n, 1:1:n);
zi_tot = total_mtx;

figure(5)
h1 = axes;
pcolor(xi,yi,log10(zi_tot));
axis square
hold on
ylabel('Source Habitat');
xlabel('Receiving Habitat');
colorbar
title('Total connectivity matrix - Moana mask passive')
c = colorbar;
c.Label.String = 'Percentage of settlement (log10 scale)';
caxis([-5 -0.5])
set(gca,'YTickLabel',{'SP','NMB 1','NMB2','Ahi','HereN','Aw', 'AwS', 'Hoki', 'Wai', 'Ara', 'CB', 'Wa', 'Waik', 'TA', 'Kai', 'PB', 'Aro', 'Oak'}, 'YTick',[1+8/2, 1+1/2+8, 1+1/2+9, 1+3/2+10, 1+1/2+13, 1+3/2+14, 1+1/2+17, 1+3/2+18, 1+2/2+21, 1+1/2+23, 1+5/2+24, 1+1/2+29, 1+2/2+30, 1+3/2+32, 1+2/2+35, 1+1/2+37, 1+4/2+38, 1+7/2+42],'fontsize',10);
set(gca,'XTickLabel',{'SP','NMB 1','NMB2','Ahi','HereN','Aw', 'AwS', 'Hoki', 'Wai', 'Ara', 'CB', 'Wa', 'Waik', 'TA', 'Kai', 'PB', 'Aro', 'Oak'}, 'XTick',[1+8/2, 1+1/2+8, 1+1/2+9, 1+3/2+10, 1+1/2+13, 1+3/2+14, 1+1/2+17, 1+3/2+18, 1+2/2+21, 1+1/2+23, 1+5/2+24, 1+1/2+29, 1+2/2+30, 1+3/2+32, 1+2/2+35, 1+1/2+37, 1+4/2+38, 1+7/2+42],'fontsize',10);
xlines= [1+8,1+9,1+10,1+13,1+14,1+17,1+18,1+21,1+23,1+24,1+29,1+30,1+32,1+35,1+37,1+38,1+42];
x_grid=[xlines; xlines];
plot(x_grid, '-', 'color', [1,1,1])
gridindex = [1+8,1+9,1+10,1+13,1+14,1+17,1+18,1+21,1+23,1+24,1+29,1+30,1+32,1+35,1+37,1+38,1+42];
makegrid(gridindex,n)
hold on
set(h1, 'Xdir', 'reverse')
set(h1, 'Ydir', 'reverse')
