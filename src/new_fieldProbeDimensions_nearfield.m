
addpath(genpath('/media/o/Permanent/wehkamp/coilSimulation'))
close all

zPos = 0;
radius =  0.001; %0.0009;
current = 0.1;
Turns = 5.5; %6;
nParts = 360; %120; % Stueckelung pro windung
wireThickness = 0.00035; %35;
straight = 0;  % straight = 1: gerade wicklung, kreisfoermige spulen; straight = 0: spiralige Spule
xPmax = 0.001; %maximale ausdehung des zu berechnenden volumens in x-
yPmax = 0.001; % y-
zPmax = 0.001; % und z-Richtung


NP = 100; % 100 gut Gridaufloesung pro Raumdimension 200 dauert lange

xP = linspace(-xPmax,xPmax,NP);        % Divide space with NP points..
yP = linspace(-yPmax, yPmax, NP);
zP = linspace(-zPmax,zPmax,NP);
[xxP yyP zzP] = meshgrid(xP, yP, zP);            % Creating the Mesh
cla;

figure(1)
hold on;

[Bx, By, Bz] = solenoidField3D (zPos, radius, current, Turns, nParts, wireThickness, xxP, yyP, zzP, straight);

B = sqrt(Bx.^2 + By.^2 + Bz.^2);


% surf(xxP(:,:,1), yyP(:,:,1), Bz(:,:,NP/2));
% h = surf(squeeze(xxP(:,:,1)), squeeze(yyP(:,:,1)), squeeze(Bz(NP/2,:,:)/0.005*250000)); %tesla/0.005*250000 = Herz
hplot = surf(squeeze(xxP(:,:,1)), squeeze(yyP(:,:,1)), squeeze(B(NP/2,:,:)/0.005*250000)); %tesla/0.005*250000 = Herz
set(hplot, 'edgecolor','none')

c = colorbar;
%colourbarlimit = 15000;
%caxis([5000 colourbarlimit]);
colourbarlimit = 13000;
caxis([8000 colourbarlimit]);%([-0.00002 0.00002]);
w = c.LineWidth;
c.LineWidth = 1.5;
alpha(.4) %.2)
alpha(c, .4) %.2)

x1=linspace(0, 2 * pi);
%kreis mit Durchmesser aeussere Spule
%plot3(radius * sin(x1), -radius * cos(x1), ones(size(x1)) .* 250000, 'LineWidth',3, 'color', 'black')

%kreis mit durchmesser innere Spule
radiusSmall = 0.0005
plot3(radiusSmall * sin(x1), -radiusSmall * cos(x1), ones(size(x1)) .* 1, 'LineWidth',6, 'color', 'black','HandleVisibility','off')


% [hplot, c] = contour(squeeze(xxP(:,:,1)), squeeze(yyP(:,:,1)), squeeze(B(NP/2,:,:)/0.005*250000),'LevelList',[10000 12000 14000 16000 18000],'ShowText','on');%,'DisplayName','frequency offset [Hz]'
[hplot, c] = contour(squeeze(xxP(:,:,1)), squeeze(yyP(:,:,1)), squeeze(B(NP/2,:,:)/0.005*250000),'LevelList',[9000 10000 11000 12000 13000],'ShowText','on');%,'DisplayName','frequency offset [Hz]'


grid on

clabel(hplot,c,'FontSize',12,'Color','black')
c.LineWidth = 3; %2
% hbar = colorbar;
colourbarlimit = 13000;
caxis([8000 colourbarlimit]);%([-0.00002 0.00002]);
%set(hbar, 'YTick',0:colourbarlimit/4:colourbarlimit,'TickLabels',{'0','5','10','15','20'})
axis square
%ylabel(hbar, 'frequency offset [Hz]')
%legend({' frequency offset [Hz]'})


xticks([-0.001 -0.0005 0 0.0005 0.001])
xticklabels({'-1' '-0.5' '0' '0.5' '1'})
yticks([-0.001 -0.0005 0 0.0005 0.001])
yticklabels({'-1' '-0.5' '0' '0.5' '1'})

xlabel('position Z-plane [mm]')
% ylabel('position Y-plane [mm]')
set(gca,'FontSize',13)
set(gca,'linewidth',1.5)
ax = gca;
ax.GridAlpha = 0.3

hold off;

% 
% x0=10;
% y0=10;
% width=550*2;%1.350;
% height=400*2;%1.350;
% set(gcf,'position',[x0,y0,width,height])
set(gcf, 'PaperUnits', 'centimeters');
x_width=7.2067 ;y_width=6.4000


savefig('nearfield.fig')

