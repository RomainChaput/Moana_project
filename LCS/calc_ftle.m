% ftle
clc;clear;

% Remove figure visibility
set(0, 'DefaultFigureVisible', 'off')

% New Zealand Coast and reefs
Reef = shaperead('Kaikoura_paua.shp');
addpath('D:\Cours\Post_doc\Moana Project\Maps_New_Zealand')
NZ = shaperead('New_Zealand_High_resolution_coast_for_plot.shp');

% Set up work environment
addpath('D:\Cours\Post_doc\Open_Drift_scripts\Matlab_post_processing')
addpath('D:\Cours\Post_doc\Open_Drift_scripts\Matlab_post_processing\output\')
N=dir('D:\Cours\Post_doc\Open_Drift_scripts\Matlab_post_processing\output\xyt-*');
f = natsortfiles({N.name}');   % sort the filenames!

for fi = 1:1800
    
    % load trajectories
    filestr=num2str(fi);
    %evalc(['load ' f(fi).name;]);
    evalc(['load ' f{fi};]);
    
    % with low-resolution ssh
    %load 'output/xyt-15-29';
    %load 'output/xyt-15-1';
    
    % derivative to form Cauchy-Green (CG) tensor
    phi_u = reshape(xf, length(x0), length(y0));
    phi_v = reshape(yf, length(x0), length(y0));
    
    % validate that y is first here
    [phi_dudy, phi_dudx] = gradient(phi_u, y0, x0);
    [phi_dvdy, phi_dvdx] = gradient(phi_v, y0, x0);
    
    % components of CG = [d11, d12; d12 d22];
    d11 = phi_dudx .* phi_dudx + phi_dvdx .* phi_dvdx;
    d12 = phi_dudx .* phi_dudy + phi_dvdx .* phi_dvdy;
    d22 = phi_dudy .* phi_dudy + phi_dvdy .* phi_dvdy;
    
    d11 = d11(:);
    d12 = d12(:);
    d22 = d22(:);
    
    lmax = zeros(length(d11), 1);
    for i = 1:length(d11)
        cg = [d11(i), d12(i); d12(i) d22(i)];
        if ~isnan(cg)
            l = eig(cg);
            lmax(i) = l(end);
        else
            lmax(i) = nan;
        end
    end
    
    lmax = reshape(lmax, length(x0), length(y0));
    ftle = 1/2 * 1/T * log(lmax);
    ftle(ftle == 0) = nan;
    
    
    
    %% plot
    prefixIN  = 'input/';
    basefile  = {'./input/nz5km_his_subset_199401000000.nc'};
    
    % load grid data
    lat = ncread('./input/nz5km_his_subset_199401000000.nc', 'lat_u');
    lat = lat(1,:);
    lon = ncread('./input/nz5km_his_subset_199401000000.nc', 'lon_u');
    lon = lon(1,:);
    lat = double(lat);
    lon = double(lon);
    lon = wrapTo180(lon);
    lon_origin = min(lon);
    lat_origin = min(lat);
    [xdat, ydat] = sph2xy(lon, lon_origin, lat, lat_origin); % [m]
    clear lon lat
    
    H=figure();
    mapshow(NZ,'FaceColor','0.7 0.7 0.7');
    hold on
    [x0, y0] = xy2sph(x0, lon_origin, y0, lat_origin);
    contourf(x0, y0, ftle', 50, 'linecolor','none');
    hold on
    xlim([172.5 177.5])
    ylim([-44.5 -40.5])
    hold on
    N = 64;   %
    Color1 = [(0:1:N-1)'/(N-1),(0:1:N-1)'/(N-1),ones(N,1)];
    Color2 = [1 1 1];
    myColor = [Color1;Color2];
    colormap(myColor)
    colorbar
    lonmin=min(xlim); lonminind=find(x0>lonmin-0.55*(abs(x0(2)-x0(1)))...
        & x0<lonmin+0.45*(abs(x0(2)-x0(1))));
    lonmax=max(xlim);    lonmaxind=find(x0>lonmax-0.55*(abs(x0(2)-x0(1)))...
        & x0<lonmax+0.45*(abs(x0(2)-x0(1))));
    latmin=min(ylim);     latminind=find(y0>latmin-0.55*(abs(y0(2)-y0(1)))...
        & y0<latmin+0.45*(abs(y0(2)-y0(1))));
    latmax=max(ylim);     latmaxind=find(y0>latmax-0.55*(abs(y0(2)-y0(1)))...
        & y0<latmax+0.45*(abs(y0(2)-y0(1))));
    cmin=min(min(ftle(lonminind:lonmaxind,latminind:latmaxind)));
    cmax=max(max(ftle(lonminind:lonmaxind,latminind:latmaxind)));
    caxis([-0.5*10^(-5) 0])
    hold on
    axis square
    mapshow(Reef);
    
    
    day = fi/24;
    title(['Kaikoura 1994, day = ',num2str(day)],'fontsize',12);
    fname = sprintf('Kaikoura_1994_FTLE%d.fig', fi);
    filename = 'D:\Cours\Post_doc\Open_Drift_scripts\Matlab_post_processing\output_LCS\';
    savefig(H, fullfile(filename,fname),'compact');
    hold off
    close
end

set(0, 'DefaultFigureVisible', 'on')