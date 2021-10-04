%Function dispK.m
% - computes dispersal kernel form Dij and Mij
% - plots the kernel
%
%
% Usage: X = dispK(D,Mtx,n,d);
%
% -----INPUT:-----
% DistMat = distance matrix
% Mtx = migration matrix
% n = matrix size
% d = bin size of frequency distribution of distances in km
%
%-----OUTPUT:-----
% - one tail 2D dispersal kernel

% -----Author:-----Claire Paris 07/07/04

function X = dispK(DistMat,Mtx,n,d,~)

set(figure, 'Visible','on');
%clear figure


%load migration matrix M
%load Mtx
output = Mtx;
bin_dist = (0:d:500)';%onetail curve
x = bin_dist;

% load distance matrix
dist = DistMat; % load DistMat
dist = dist(1:n,1:n);

%loop through the polygon numbers
res2=[];
%npart = 26700; % number of particles released by polygon

for ii=1:n
    
    polstr = num2str(ii);
    out =output(ii,:);%source
    %extract distance
    Y = dist(ii,:);
    
    y1 = [Y' out'];
    y = sortrows(y1,1);%sort input by distance
    
    w   = y(:,2);    % extract number or recruits (weights)
    ww   = w;%/ npart; % divide each by npart released per polygon
    
    y       = y(:,1);    % extract distances
    n       = hist(y,x); % get bin counts
    noBins  = size(n,2); % get number of bins
    res     = zeros(1,noBins); % initialize row vector
    counter = 1;               % initialize counter
    
    for j = 1:noBins
        if (n(j) ~= 0)
            res(j) = sum(ww(counter:(counter + (n(j)-1))));
            counter = counter + n(j); % increment counter
        end
    end
    
    %Number of polygons at each d km bin
    res2 = [res2;res];
    %std2(ii,:) = std(res2(ii,:));
    
end %ii loop

res2(isnan(res2))=0;
% Choose the polygon to plot
res1 = mean(res2);
%res1 = res2(44,:);
std1 = std(res2);
X = [res1;std1];

%-----Plot histogram:-----

%h2 = errorbar(x,res1,std1,'y');%error for import
d1 = res1';
hold on;
%mavg_d1=filter(ones(1,window)/window,1,d1);
%mavg_d1=mavg_d1';
%h = plot(x, mavg_d1, 'k', 'linewidth', 2); % for smooth line
%h= bar(x,res1,'k');                           % for bar graph
%h = plot(x,res1,'k');
%scatter(x,res1,'ko')
plot(x, res1,'k-',...
     x, res1+std1,'r--',...
     x, res1-std1,'r--');
   h = findobj(gca,'Type','patch');
set(h,'FaceColor',[0.8 0.8 0.8],'EdgeColor','k');
hold on;

%mu = [mean(N) std(N)];
%m=8;
%A=polyfit(x,d1,m);
%Yest=polyval(A,x);%calculates regression
%k1=plot(x,Yest,'r-');
%z = interp1(x, res1, x, 'spline');
%k1 = plot(x,z,'r-');
%set(k1,'linewidth',[2],'markersize',[5]);%'markerfacecolor','r'
set(gca,'FontSize', 12,'xlim',[0 400],'ylim',[0 2*10^-2]);

ylabel('Proportion of settlement per reef polygon');
xlabel('Distance (km)');
saveas(gcf,'Dispersal_kernel','epsc')
end