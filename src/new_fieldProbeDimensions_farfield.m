
addpath(genpath('/media/o/Permanent/wehkamp/coilSimulation'))
close all

zPos = 0;
radius = 0.001;
current = 0.1;
Turns = 6;
nParts = 360; %120; % Stueckelung pro windung
% wireThickness = 0.0002;
wireThickness = 0.0003;
straight = 0;  % straight = 1: gerade wicklung, kreisfoermige spulen; straight = 0: spiralige Spule
xPmax = 0.05; %maximale ausdehung des zu berechnenden volumens in x-
yPmax = 0.05; % y-
zPmax = 0.05; % und z-Richtung
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


%surf(xxP(:,:,1), yyP(:,:,1), Bz(:,:,NP/2));
% h = surf(squeeze(xxP(:,:,1)), squeeze(yyP(:,:,1)), squeeze(Bz(NP/2,:,:)/0.005*250000)); %tesla/0.005*250000 = Herz
% hplot = surf(squeeze(xxP(:,:,1)), squeeze(yyP(:,:,1)), squeeze(B(NP/2,:,:)/0.005*250000)); %tesla/0.005*250000 = Herz
% 
% set(hplot, 'edgecolor','none')
% x1=linspace(0, 2 * pi);
% c = colorbar;
% w = c.LineWidth;
% c.LineWidth = 1.5;
% 
% %kreis mit Durchmesser aeussere Spule
% %plot3(radius * sin(x1), -radius * cos(x1), ones(size(x1)) .* 250000, 'LineWidth',3, 'color', 'black')
% 
%kreis mit durchmesser innere Spule
radiusSmall = 0.0005
plot3(radiusSmall * sin(x1), -radiusSmall * cos(x1), ones(size(x1)) .* 1, 'LineWidth',3, 'color', 'black','HandleVisibility','off')

% %gestreift aeussere spule
% % x0 = linspace(-radius, radius,10)
% % x1 = linspace(wireThickness/4,-wireThickness/4,10)
% % % coilStart = -(Turns/2-0.5)*wireThickness
% % % coilEnd = (Turns/2-0.5)*wireThickness
% % 
% % for offset = -wireThickness*(Turns-1)/2 : wireThickness : (wireThickness)*(Turns-1)/2 ; %coilStart:wireThickness:coilEnd
% %     plot3(ones(size(x0))* offset + x1, x0, ones(size(x0))* 250000, 'LineWidth',3, 'color', 'black')
% % end
% 
% %gestreift innere spule
% % innerTurns = 4
% % innerCoilStart = -(innerTurns-0.5)/2*wireThickness
% % innerCOilEnd = (innerTurns-0.5)/2*wireThickness
% % x2 = linspace(-0.0005, 0.0005,10)
% % x3 = linspace(-0.00005,0.00005,10)
% % for offset = innerCoilStart:wireThickness:innerCOilEnd
% %     plot3(x2, ones(size(x2))* offset + x3, ones(size(x2))* 250000, 'LineWidth',3, 'color', 'black')
% % end
% 
% hbar = colorbar;
% ylabel(hbar, 'frequency offset [Hz]')
% 
% %caxis([0 0.2]);%([-0.00002 0.00002]);
% % caxis([-2000 2000]);%([-0.00002 0.00002]);
% 
% xticks([-0.001 -0.0005 0 0.0005 0.001])
% xticklabels({'-1' '-0.5' '0' '0.5' '1'})
% yticks([-0.001 -0.0005 0 0.0005 0.001])
% yticklabels({'-1' '-0.5' '0' '0.5' '1'})
% xlabel('position Z-plane [mm]')
% ylabel('position Y-plane [mm]')
% set(gca,'FontSize',20)

[hplot, c] = contour(squeeze(xxP(:,:,1)), squeeze(yyP(:,:,1)), squeeze(B(NP/2,:,:)/0.005*250000),'LevelList',[0 0.25 0.5 1 2.5 5 10 50],'ShowText','on');%,'DisplayName','frequency offset [Hz]'
grid on
clabel(hplot,c,'FontSize',12,'Color','black')
c.LineWidth = 3
% hbar = colorbar;
colourbarlimit = 20;
caxis([0 colourbarlimit]);%([-0.00002 0.00002]);
%set(hbar, 'YTick',0:colourbarlimit/4:colourbarlimit,'TickLabels',{'0','5','10','15','20'})
axis square
%ylabel(hbar, 'frequency offset [Hz]')
legend({' frequency offset [Hz]'})
xticks([-zPmax -zPmax*0.5 0 zPmax*0.5 zPmax])
xticklabels({'-50' '-25' '0' '25' '50'})
yticks([-zPmax -zPmax*0.5 0 zPmax*0.5 zPmax])
yticklabels({'-50' '-25' '0' '25' '50'})

xlabel('position Z-plane [mm]')
ylabel('position Y-plane [mm]')
set(gca,'FontSize',13)
set(gca,'linewidth',1.5)
ax = gca;
ax.GridAlpha = 0.3

hold off;

% 
% x0=10;
% y0=10;
% width=550*1.50;
% height=64.000 mm;
% set(gcf,'position',[x0,y0,width,height])
set(gcf, 'PaperUnits', 'centimeters');
x_width=7.2067 ;y_width=6.4000

savefig('/raid/home/extern/wehkamp/shifted_field_probe/coilSimulation/new_farfield_contour_2000_fontSize20.fig')

% figure(2)
% hplot = contour(squeeze(xxP(:,:,1)), squeeze(yyP(:,:,1)), squeeze(B(NP/2,:,:)/0.005*250000),200); %tesla/0.005*250000 = Herz
% hbar = colorbar;
% ylabel(hbar, 'frequency offset [Hz]')
% caxis([0 2000]);%([-0.00002 0.00002]);
