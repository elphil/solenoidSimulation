function [Bx, By, Bz] = solenoidField3D (zPos, radius, current, Turns, nParts, wireThickness, xxP, yyP, zzP, straight)
error = 500;
m = current*10^-7;                      % mu_0/2*pi * current through the Loop... 1(ANTI-CLOCK), -1(CLOCK)

dl = 2*pi*radius/nParts;                %NW Length of each  element..

%Change to all in one turn
allDeg = Turns*360;
nParts = Turns*nParts;
dtht = allDeg/nParts;                      %NW Loop elements divided w.r.to degrees..
tht = (0+dtht/2): dtht : (allDeg-dtht/2);  %NW Angle of each element w.r.to origin..

xC =  radius.*cosd(tht).*ones(1,nParts);        % X coordinates...
yC =  radius.*sind(tht).*ones(1,nParts);        % Y coordinates...

if (straight)
    zC = zeros(1,nParts); %NWkeine Ahnung ob das noch richtig ist...
else
    zC = (1:nParts)*wireThickness/(nParts/Turns);
end

h = wireThickness*(Turns)/2 %0.550e-03 


% Length(Projection) & Direction of each current element in Vector form..
Lx = dl.*sind(tht);      % Length of each element on X axis..
Ly = dl.*cosd(tht);      % Length of each element on Y axis..          
if (straight)
    Lz = zeros(1,nParts); 
else
    Lz = (wireThickness/nParts).*ones(1,nParts);   % Length of each element is zero on Z axis..            
end

% Initialize B..
Bx = zeros(size(xxP, 1), size(xxP, 1), size(xxP, 1));
By = zeros(size(xxP, 1), size(xxP, 1), size(xxP, 1));
Bz = zeros(size(xxP, 1), size(xxP, 1), size(xxP, 1));

% Computation of Magnetic Field (B) using Superposition principle..
% Compute B at each detector points due to each small cond elements & integrate them..
errorX = 0;
errorY = 0;
errorZ = 1.01;
% plot3(xC(:),yC(:),zC(:))

coil = plot3(-(zPos-h+zC), xC, yC*10e5/radius, 'LineWidth',3, 'color', 'black','HandleVisibility','off');
% coil.Color(4) = 0.5;

    parfor q = 1:nParts
      
        rx = xxP - xC(q);               
        ry = yyP + yC(q);  
        rz = (zPos + zzP - h + zC(q));       % zC Vertical location    all in one Turn
%         rz = (zPos + zzP - h(p) + zC(q));       % zC Vertical location changes per Turn..
    
        r = sqrt(rx.^2+ry.^2+rz.^2);    % Displacement Magnitude for an element on the conductor..    
        r3 = r.^3;
        BxTmp = m.*Ly(q).*rz./r3 - m*Lz(q).*ry./r3;
        ByTmp = m.*Lz(q).*rx./r3 - m*Lx(q).*rz./r3;
        BzTmp = m.*Lx(q).*ry./r3 - m*Ly(q).*rx./r3;
        Bx = Bx + BxTmp%  m.*Ly(q).*rz./r3 - m*Lz(q).*ry./r3;
        By = By + ByTmp%  m.*Lz(q).*rx./r3 - m*Lx(q).*rz./r3;
        Bz = Bz + BzTmp%  m.*Lx(q).*ry./r3 - m*Ly(q).*rx./r3;
    end



