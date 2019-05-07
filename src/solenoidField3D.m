function [Bx, By, Bz] = solenoidField3D (zPos, radius, current, Turns, nParts, wireThickness, xxP, yyP, zzP, straight)
error = 500;
m = current*10^-7;                      % mu_0/2*pi * coil current
dl = 2*pi*radius/nParts;                % Length of each element

dTheta = 360/nParts;                      % angular increment d theta
theta = (0+dTheta/2): dTheta : (360-dTheta/2);  % angles as a list

xC =  radius.*cosd(theta).*ones(1,nParts);        % X coordinates...
yC =  radius.*sind(theta).*ones(1,nParts);        % Y coordinates...

% zC will change for each turn, -Z to +Z, and will be varied during iteration but - if (straight) - zC = [0.. 0]..
if (straight)
    zC = zeros(1,nParts);
else
    zC = (1:nParts)*wireThickness/nParts;
end
h = -wireThickness*(Turns-1)/2 : wireThickness : (wireThickness)*(Turns-1)/2;


% Length(Projection) & Direction of each current element in Vector form..
Lx = dl.*sind(theta);      % x-projection dl
Ly = dl.*cosd(theta);      % y-projection of dl
if (straight)
    Lz = zeros(1,nParts); 
else
    Lz = (wireThickness/nParts).*ones(1,nParts);              
end


% Initialize B components
Bx = zeros(size(xxP, 1), size(xxP, 1), size(xxP, 1));
By = zeros(size(xxP, 1), size(xxP, 1), size(xxP, 1));
Bz = zeros(size(xxP, 1), size(xxP, 1), size(xxP, 1));

% Computation of Magnetic Field, biot savart superposition
% Compute B at each detector point
% errorX = 0;
% errorY = 0;
% errorZ = 1.01;
counter = 0;
for p = 1:Turns
    counter = counter + 1;
    for q = 1:nParts
        rx = xxP - xC(q);               
        ry = yyP + yC(q);             
        rz = (zPos + zzP - h(p) + zC(q));
    
        r = sqrt(rx.^2+ry.^2+rz.^2);    % r_ges calculation
        r3 = r.^3;

        %outer product according to Biot Savart
        Bx = Bx +   m.*Ly(q).*rz./r3 - m*Lz(q).*ry./r3;
        By = By +   m.*Lz(q).*rx./r3 - m*Lx(q).*rz./r3;
        Bz = Bz +   m.*Lx(q).*ry./r3 - m*Ly(q).*rx./r3;
    end
end


