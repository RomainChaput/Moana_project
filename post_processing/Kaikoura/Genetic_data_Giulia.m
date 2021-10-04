%% %%% Code to compare the genetic connectivity to the oceanographic connectivity at Kaikoura

%ADULTS 1-FST									
adult_FST=[0,0.999854839,0.999444423,0.999467321,0.999529126,0.99965921,0.999925277,0.997437771,0
	0,0,0.999399599,0.9999256,0.999689545,0.999389829,0.999654287,0.99783708,0
	0,0,0,0.999488704,0.99933971,0.999804531,0.999430607,0.99787105,0
	0,0,0,0,0.999410193,0.999965846,0.999787041,0.997693703,0
	0,0,0,0,0,0.999899762,0.999855302,0.9986196,0
	0,0,0,0,0,0,0.999615483,0.997855046,0
	0,0,0,0,0,0,0,0.99751664,0
	0,0,0,0,0,0,0,0,0
    0,0,0,0,0,0,0,0,0];

%% Plot the genetic connectivity

n = 9;
%fst matrix
[xi,yi]=meshgrid(1:1:n, 1:1:n);
zi = (adult_FST);
figure()
hi = pcolor(xi,yi,zi);
set(hi, 'EdgeColor', 'none');
hold on
axis square
colormap(bluewhitered);
ylabel('Source Habitat');
xlabel('Receiving Habitat');
colorbar
title('Genetic Adults - Kaikoura Paua')
c = colorbar;
caxis([0.997 1])
c.Label.String = 'Percentage of settlement';
set(gca,'YTickLabel',{'BoatH', 'Blocks', 'Katu', 'KaiP', 'HBay', 'Papa', 'Ward', 'CapeC'}, ...
    'YTick',[1+1/2, 1+1/2+1, 1+1/2+2, 1+1/2+3, 1+1/2+4, 1+1/2+5, 1+1/2+6, 1+1/2+7],'fontsize',10);
set(gca,'XTickLabel',{'BoatH', 'Blocks', 'Katu', 'KaiP', 'HBay', 'Papa', 'Ward', 'CapeC'}, ...
    'XTick',[1+1/2, 1+1/2+1, 1+1/2+2, 1+1/2+3, 1+1/2+4, 1+1/2+5, 1+1/2+6, 1+1/2+7],'fontsize',10);
c.Label.String = 'Inverse of Fst value';

%% Juvenile 1-FST

juvenile_FST=[0,0.999462513,1,0.999356663,0.998866528,0.999041937,0.999257819,0.99867092,0
    0,0,0.999939681,0.999871365,0.999006579,0.999143752,0.999472942,0.998244523,0
    0,0,0,1,0.999822521,0.999804076,0.999826889,0.999538773,0
    0,0,0,0,0.999697724,0.999544619,1,0.999884678,0
    0,0,0,0,0,0.999685962,0.998667215,0.998494203,0
    0,0,0,0,0,0,0.999295542,0.998909507,0
    0,0,0,0,0,0,0,0.999056162,0
    0,0,0,0,0,0,0,0,0
    0,0,0,0,0,0,0,0,0];

%% Plot the genetic connectivity

n = 9;
%fst matrix
[xi,yi]=meshgrid(1:1:n, 1:1:n);
zi = (juvenile_FST);
figure()
hi = pcolor(xi,yi,zi);
set(hi, 'EdgeColor', 'none');
hold on
axis square
colormap(bluewhitered);
ylabel('Source Habitat');
xlabel('Receiving Habitat');
colorbar
title('Genetic Juveniles - Kaikoura Paua')
c = colorbar;
caxis([0.998 1])
c.Label.String = 'Percentage of settlement';
set(gca,'YTickLabel',{'BoatH', 'Blocks', 'Katu', 'KaiP', 'HBay', 'Papa', 'Ward', 'CapeC'}, ...
    'YTick',[1+1/2, 1+1/2+1, 1+1/2+2, 1+1/2+3, 1+1/2+4, 1+1/2+5, 1+1/2+6, 1+1/2+7],'fontsize',10);
set(gca,'XTickLabel',{'BoatH', 'Blocks', 'Katu', 'KaiP', 'HBay', 'Papa', 'Ward', 'CapeC'}, ...
    'XTick',[1+1/2, 1+1/2+1, 1+1/2+2, 1+1/2+3, 1+1/2+4, 1+1/2+5, 1+1/2+6, 1+1/2+7],'fontsize',10);
c.Label.String = 'Inverse of Fst value';

