
addpath(genpath('/media/o/Permanent/wehkamp/coilSimulation'))
close all

zPos = 0;
radius = 0.001;
current = 0.1;
Turns = 5.5;%6;
nParts = 360; %120; % Stueckelung pro windung
wireThickness = 0.0002;
straight = 0;  % straight = 1: gerade wicklung, kreisfoermige spulen; straight = 0: spiralige Spule
xPmax = 0.005; %maximale ausdehung des zu berechnenden volumens in x-
yPmax = 0.005; % y-
zPmax = 0.005; % und z-Richtung
NP = 100; % Gridaufloesung pro Raumdimension 200 dauert lange

xP = linspace(-xPmax,xPmax,NP);        % Divide space with NP points..
yP = linspace(-yPmax, yPmax, NP);
zP = linspace(-zPmax,zPmax,NP);
[xxP yyP zzP] = meshgrid(xP, yP, zP);            % Creating the Mesh
cla;

figure(1)
hold on;

[Bx, By, Bz] = solenoidField3D (zPos, radius, current, Turns, nParts, wireThickness, xxP, yyP, zzP, straight);

B = sqrt(Bx.^2 + By.^2 + Bz.^2);

%kreis mit durchmesser innere Spule
radiusSmall = 0.0005
plot3(radiusSmall * sin(x1), -radiusSmall * cos(x1), ones(size(x1)) .* 1, 'LineWidth',3, 'color', 'black','HandleVisibility','off')



[hplot, c] = contour(squeeze(xxP(:,:,1)), squeeze(yyP(:,:,1)), squeeze(B(NP/2,:,:)/0.005*250000),'LevelList',[250 500 1000 2000 5000],'ShowText','on','DisplayName','frequency offset [Hz]');
grid on
clabel(hplot,c,'FontSize',12,'Color','black')
c.LineWidth = 3
% hbar = colorbar;
colourbarlimit = 5000;
caxis([0 colourbarlimit]);%([-0.00002 0.00002]);
%set(hbar, 'YTick',0:colourbarlimit/4:colourbarlimit,'TickLabels',{'0','5','10','15','20'})
axis square
%ylabel(hbar, 'frequency offset [Hz]')
legend({' frequency offset [Hz]'})

xticks([-zPmax -zPmax*0.5 0 zPmax*0.5 zPmax])
xticklabels({'-5' '-2.5' '0' '2.5' '5'})
yticks([-zPmax -zPmax*0.5 0 zPmax*0.5 zPmax])
yticklabels({'-5' '-2.5' '0' '2.5' '5'})

xlabel('position Z-plane [mm]')
% ylabel('position Y-plane [mm]')
set(gca,'FontSize',13)


set(gca,'linewidth',1.5)
ax = gca;
ax.GridAlpha = 0.3

hold off;

% x0=10;
% y0=10;
% width=550*1.350;
% height=400*1.350;
% set(gcf,'position',[x0,y0,width,height])
set(gcf, 'PaperUnits', 'centimeters');
x_width=7.2067 ;y_width=6.4000

savefig('midfield.fig')

