function [fieldSpread, innerField] = hyperField3D(lCompensationOne, lCompensationTwo, current)
N = 50;                 % Number sections/elements in each turn of the conductor...
wireThickness = .001;
wireHeight = 0.001+0.0003;
rInner = .19/2;
lTot = 0.3530;
lTripleWind = lCompensationOne;
currentComp = -lCompensationTwo;
%lDoubleWind = .043;


% Points/Locations in space (here XZ plane) where B is to be computed..
NP = 30;              % Detector points..
%xPmax = 0.5*rInner;    % Dimensions of detector space, arbitrary
%zPmax = 0.6*lTot;
% xPmax = rInner;
% yPmax = rInner;
% zPmax = lTot/2;
xPmax = 0.025;
yPmax = 0.025;
zPmax = 0.025;

xP = linspace(-xPmax,xPmax,NP);        % Divide space with NP points..
yP = linspace(-yPmax, yPmax, NP);
zP = linspace(-zPmax,zPmax,NP);
[xxP yyP zzP] = meshgrid(xP, yP, zP);            % Creating the Mesh


Bx = zeros(NP, NP, NP);
By = zeros(NP, NP, NP);
Bz = zeros(NP, NP, NP);

tmpX = zeros(NP, NP, NP);
tmpY = zeros(NP, NP, NP);
tmpZ = zeros(NP, NP, NP);
%1. Lage
[tmpX, tmpY, tmpZ] = solenoidField3D(0, rInner + 0.5 * wireHeight, current, lTot/wireThickness, NP, wireThickness, xxP, yyP, zzP, 0);
Bx = Bx + tmpX;
By = By + tmpY;
Bz = Bz + tmpZ;

% 2. Lage
[tmpX, tmpY, tmpZ] = solenoidField3D(0, rInner + 1.5 * wireHeight, current, lTot/wireThickness, NP, wireThickness, xxP, yyP, zzP, 0);
Bx = Bx + tmpX;
By = By + tmpY;
Bz = Bz + tmpZ;

% % 3. Lage
% [tmpX, tmpY, tmpZ] = solenoidField3D(0, rInner + 2.5 * wireHeight, current, lTot/wireThickness, N, wireThickness, xxP, yyP, zzP);
% Bx = Bx + tmpX;
% By = By + tmpY;
% Bz = Bz + tmpZ;

% [tmpX, tmpY, tmpZ] = solenoidField3D( lTot/2 - lTripleWind/2, rInner + 3.5 * (wireHeight+0.0003), current, lTripleWind/(wireThickness+0.00025), N, wireThickness, xxP, yyP, zzP);
% Bx = Bx + tmpX;
% By = By + tmpY;
% Bz = Bz + tmpZ;
% 
% [tmpX, tmpY, tmpZ] = solenoidField3D(-lTot/2 + lTripleWind/2, rInner + 3.5 * (wireHeight+0.0003), current, lTripleWind/(wireThickness+0.00025), N, wireThickness, xxP, yyP, zzP);
% Bx = Bx + tmpX;
% By = By + tmpY;
% Bz = Bz + tmpZ;

% %3. Lage
% [tmpX, tmpY, tmpZ] = solenoidField3D(0, rInner + 2.5 * wireHeight, current, lTot/wireThickness, N, wireThickness, xxP, yyP, zzP);
% Bx = Bx + tmpX;
% By = By + tmpY;
% Bz = Bz + tmpZ;

% %4. Lage
% [tmpX, tmpY, tmpZ] = solenoidField3D(0, rInner + 3.5 * wireHeight, current, lTot/wireThickness, N, wireThickness, xxP, yyP, zzP);
% Bx = Bx + tmpX;
% By = By + tmpY;
% Bz = Bz + tmpZ;

% %4. Lage
% [tmpX, tmpY, tmpZ] = solenoidField3D(0, rInner + 4.5 * wireHeight, current, lTot/wireThickness, N, wireThickness, xxP, yyP, zzP);
% Bx = Bx + tmpX;
% By = By + tmpY;
% Bz = Bz + tmpZ;

% %5. Lage
% [tmpX, tmpY, tmpZ] = solenoidField3D(0, rInner + 5.5 * wireHeight, current, lTot/wireThickness, N, wireThickness, xxP, yyP, zzP);
% Bx = Bx + tmpX;
% By = By + tmpY;
% Bz = Bz + tmpZ;

% %6. Lage - compensation
% [tmpX, tmpY, tmpZ] = solenoidField3D(0, rInner + 6.5 * wireHeight, currentComp, lTripleWind/wireThickness, N, wireThickness, xxP, yyP, zzP);
% Bx = Bx + tmpX;
% By = By + tmpY;
% Bz = Bz + tmpZ;
% 
% %7. Lage - compensation
% [tmpX, tmpY, tmpZ] = solenoidField3D(0, rInner + 7.5 * wireHeight, currentComp, lTripleWind/wireThickness, N, wireThickness, xxP, yyP, zzP);
% Bx = Bx + tmpX;
% By = By + tmpY;
% Bz = Bz + tmpZ;


% 
% %6. Lage - compensation
% [tmpX, tmpY, tmpZ] = solenoidField3D(0, rInner + 4.5 * wireHeight, -1 * current, lCompensationTwo/wireThickness, N, wireThickness, xxP, yyP, zzP);
% Bx = Bx + tmpX;
% By = By + tmpY;
% Bz = Bz + tmpZ;

% 
[tmpX, tmpY, tmpZ] = solenoidField3D( lTot/2 - lTripleWind/2, rInner + 2.5 * (wireHeight+0.0003), current, lTripleWind/(wireThickness+0.00025), N, (wireThickness+0.00025), xxP, yyP, zzP, 1);
Bx = Bx + tmpX;
By = By + tmpY;
Bz = Bz + tmpZ;

[tmpX, tmpY, tmpZ] = solenoidField3D(-lTot/2+ lTripleWind/2, rInner + 2.5 * (wireHeight+0.0003), current, lTripleWind/(wireThickness+0.00025), N, (wireThickness+0.00025), xxP, yyP, zzP, 1);
Bx = Bx + tmpX;
By = By + tmpY;
Bz = Bz + tmpZ;

% [tmpX, tmpY, tmpZ] = solenoidField3D( lTot/2 - lCompensationTwo/2, rInner + 3.5 * (wireHeight+0.0003), current, lCompensationTwo/(wireThickness+0.00025), N, (wireThickness+0.00025), xxP, yyP, zzP);
% Bx = Bx + tmpX;
% By = By + tmpY;
% Bz = Bz + tmpZ;
% 
% [tmpX, tmpY, tmpZ] = solenoidField3D(-lTot/2+ lCompensationTwo/2, rInner + 3.5 * (wireHeight+0.0003), current, lCompensationTwo/(wireThickness+0.00025), N, (wireThickness+0.00025), xxP, yyP, zzP);
% Bx = Bx + tmpX;
% By = By + tmpY;
% Bz = Bz + tmpZ;

% [tmpX, tmpY, tmpZ] = solenoidField3D(0, rInner + 0.5 * wireHeight, current, lTot/wireThickness, N, wireThickness, xxP, yyP, zzP);
% Bx = Bx + tmpX;
% By = By + tmpY;
% Bz = Bz + tmpZ;
% 
% 
% [tmpX, tmpY, tmpZ] = solenoidField3D(0, rInner + 1.5 * wireHeight, current, lTot/wireThickness, N, wireThickness, xxP, yyP, zzP);
% Bx = Bx + tmpX;
% By = By + tmpY;
% Bz = Bz + tmpZ;
% 
% 
% [tmpX, tmpY, tmpZ] = solenoidField3D( lTot/2 - lTripleWind/2, rInner + 2.5 * wireHeight, current, lTripleWind/wireThickness, N, wireThickness, xxP, yyP, zzP);
% Bx = Bx + tmpX;
% By = By + tmpY;
% Bz = Bz + tmpZ;
% 
% [tmpX, tmpY, tmpZ] = solenoidField3D(-lTot/2+ lTripleWind/2, rInner + 2.5 * wireHeight, current, lTripleWind/wireThickness, N, wireThickness, xxP, yyP, zzP);
% Bx = Bx + tmpX;
% By = By + tmpY;
% Bz = Bz + tmpZ;
% 
% [tmpX, tmpY, tmpZ] = solenoidField3D( lTot/2 - lCompensationTwo/2, rInner + 3.5 * wireHeight, current, lCompensationTwo/wireThickness, N, wireThickness, xxP, yyP, zzP);
% Bx = Bx + tmpX;
% By = By + tmpY;
% Bz = Bz + tmpZ;
% 
% [tmpX, tmpY, tmpZ] = solenoidField3D(-lTot/2+ lCompensationTwo/2, rInner + 3.5 * wireHeight, current, lCompensationTwo/wireThickness, N, wireThickness, xxP, yyP, zzP);
% Bx = Bx + tmpX;
% By = By + tmpY;
% Bz = Bz + tmpZ;


%Addition of Earth magnetic field in every component considering a
%45-degree angle towards the N-S component and parallel alignment towards
%the earths surface.

%BeNS = 21.4*10^-6;
%BeEW = -1.2.4*10^-6;
%BeUD = 42.3*10^-6;
% NS = 21 * 10^-6;
% OW = 1.2 * 10^-6;
% UD = 42.3 * 10^-6;

%earth magnetic field considering parallelity to the floor and  
% Bz = Bz + NS * cosd(45) + OW * sind (45);
% By = By + NS * sind(45) + OW * cosd (45);
% Bx = Bx + UD;

B = sqrt(Bx.^2 + By.^2 + Bz.^2);        % Magnitude of B
innerBx = zeros(size(B,1), size(B,2), size(B,3));
innerBy = zeros(size(B,1), size(B,2), size(B,3));
innerBz = zeros(size(B,1), size(B,2), size(B,3));
for i = 1 : NP
    for j = 1 : NP
        for k = 1 : NP
            if sqrt(xxP(i, j, k)^2 + zzP(i, j, k).^2) <= zPmax && abs(yyP (i, j, k)) <= yPmax
                innerBx(i, j, k) = Bx(i, j, k);
                innerBy(i, j, k) = By(i, j, k);
                innerBz(i, j, k) = Bz(i, j, k);
%                 innerB = [innerB, (B(i, j, k))];   
            else
                innerBx(i, j, k) = NaN;
                innerBy(i, j, k) = NaN;
                innerBz(i, j, k) = NaN;
            end
            
        end
    end
end
diffInnerBx = innerBx -  nanmean(nanmean(nanmean(innerBx)));
diffInnerBy = innerBy -  nanmean(nanmean(nanmean(innerBy)));
diffInnerBz = innerBz -  nanmean(nanmean(nanmean(innerBz)));
innerB = sqrt(innerBx.^2 + innerBy.^2 + innerBz.^2);
innerField = innerB;
stdevB = std2(B);
meanB = mean(mean(mean(B)));
diffInnerB = innerB - nanmean(nanmean(nanmean(innerB)));
maxB = max(max(max(B)));
minB = min(min(min(B)));
devB = B - meanB;
steps = linspace(minB, maxB, 20);
fieldSpread = abs((min(B(size(B, 3)/2, :, size(B,3)/2))- max(B(size(B, 3)/2, :, size(B,3)/2)))/min(B(size(B, 3)/2, :, size(B,3)/2)))
figure(2);
clf
size(B)
surf(xxP(:,:,1), yyP(:,:,1), B(:,:,size(B, 3)/2));
% hold on;
% cmap = colormap(jet(100));
% for i = 1 : 2: NP
%     
%     for j = 1 : 2 : NP
%         for k = 1 : 2: NP
%             if ~isnan(innerB(i, j, k))
%                 color = (sqrt((innerBx(i, j, k)^2 + innerBy(i, j, k)^2 + innerBz(i, j, k)^2))-minB)/(maxB-minB);
%                 hnew = quiver3(xxP(i, j, k), yyP(i, j, k), zzP(i, j, k), innerBx(i, j, k), innerBy(i, j, k), innerBz(i, j, k), 0.25,...
%                     'color', [cmap(ceil(100 * color), 1) , cmap(ceil(100 * color), 2) ,cmap(ceil(color * 100), 3)], 'MarkerSize', 10);
%             end
%         end
%     end
% end
% for i = 1 : 2: NP
%     
%     for j = 1 : 2 : NP
%         for k = 1 : 2: NP
%             if ~isnan(innerB(i, j, k))
%                 color = (sqrt((innerBx(i, j, k)^2 + innerBy(i, j, k)^2 + innerBz(i, j, k)^2))-minB)/(maxB-minB);
%                 hnew = quiver3(xxP(i, j, k), yyP(i, j, k), zzP(i, j, k), diffInnerBx(i, j, k), diffInnerBy(i, j, k), diffInnerBz(i, j, k), 100,...
%                     'color', [cmap(ceil(100 * color), 1) , cmap(ceil(100 * color), 2) ,cmap(ceil(color * 100), 3)], 'MarkerSize', 100);
%             end
%         end
%     end
% end
% colorbar('YTickLabel', minB: (maxB-minB)/10: maxB);
% step = 0.1;
% for i = 0:step:2*pi
%     plot3([cos(i)*xPmax cos(i+step)*xPmax], [yPmax yPmax], [sin(i)*xPmax sin(i+step)* xPmax]);
%     plot3([cos(i)*xPmax cos(i+step)*xPmax], [-yPmax -yPmax], [sin(i)*xPmax sin(i+step)* xPmax]);
% end
% plot3([-xPmax -xPmax], [-yPmax yPmax], [0 0]);
% plot3([xPmax xPmax], [-yPmax yPmax], [0 0]);
% xlabel('x-dir');
% ylabel('y-dir');
% zlabel('z-dir');
%    quiver3(xxP, yyP, zzP, innerBx, innerBy, innerBz);
% grid off;
% surf(xxP(:,:,size(xxP,3)/2), yyP(:,:,size(yyP,3)/2),...
%     (diffInnerB( :,:, size(diffInnerB, 3)/2 ) - nanmean(nanmean(diffInnerB(:, :, size(diffInnerB, 3)/2)))) * 100000);
% 
% quiver3(xxP(:,:, size(xxP, 3)/2), yyP(:,:, size(yyP, 3)/2), innerBx(:,:,size(innerBx, 3)/2), innerBy(:, :, size(innerBy, 3)/2));
% colormap jet;


% Plotting...

% figure(1);
% pcolor(xxP,zzP,B);
% colormap(jet);
% shading interp;
% axis equal;
% axis([-xPmax xPmax -zPmax zPmax]);
% xlabel('<-- x -->');ylabel('<-- z -->');
% title('Magnetic Field Distribution');
% colorbar;



% surf(xxP,yyP,B,'FaceColor','interp',...
%     'EdgeColor','none',...
%     'FaceLighting','phong');

%draw Coil volume
line([-rInner, rInner], [-lTot/2, -lTot/2], [meanB, meanB]);
line([-rInner, rInner], [lTot/2, lTot/2], [meanB, meanB]);
line([-rInner, -rInner], [lTot/2, -lTot/2], [meanB, meanB]);
line([rInner, rInner], [lTot/2, -lTot/2], [meanB, meanB]);
%draw reactor Volume

% theta = linspace(0, 2*pi, 30);    %# divide circle by N points (length of data)
% r = 0.015;                        %# radius
% x = r.*cos(theta);               %# x-coordinate
% y = r.*sin(theta);               %# y-coordinate
% line(x, y, meanB);

line([-0.015, 0.015], [-0.015, -0.015], [maxB, maxB]);
line([-0.015, 0.015], [0.015, 0.015], [maxB, maxB]);
line([-0.015, -0.015], [-0.015, 0.015], [maxB, maxB]);
line([0.015, 0.015], [-0.015, 0.015], [maxB, maxB]);

daspect([1 1 meanB]);
axis tight;
view(0,90);
%camlight off;
colormap(jet);
grid on;
axis on;
colorbar;
title('Magnetic Field Distibution');


 
figure(3);
quiver(xxP,yyP,Bx,By);
colormap(lines);
%axis([-2*d 2*d -T*g T*g]);
title('Magnetic Field Distibution');
xlabel('<-- x -->');ylabel('<-- z -->');
zoom on;
end