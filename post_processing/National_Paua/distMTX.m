% function distMTX.m
% - computes distance matrix between each release location pair ij (node) 
% - save the matrix in DistMat 
% - plots the matrix
%
% -----INPUT:-----
% LL = matrix with degrees of longitude (1st column) and degrees of latitudes (2nd column)
%
% -----OUTPUT:-----
% - symmetric pairwise distance matrix (in kilometers)

function X = distMTX(LL);

D = f_latlong(LL); 

save DistMat D

n=length(D);
[xi,yi]=meshgrid(1:1:n, 1:1:n);
zi = [D];

figure
pcolor(xi,yi,zi);hold on
shading interp
axis square
ylabel('Node i');
xlabel('Node j');
colorbar
end


