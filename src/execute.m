current = 1.6123;
hyperField3D(0.092, 0, current);

f = @(current)hyperField3D(0.092, 0, current)
rInner = 0.50;
dInner = 0.50;
rOuter = 0.12;
dOuter = 3.64;
windings = 1;
nWindingsInnerZ = windings;
nWindingsInnerX = windings;
nWindingsOuterZ = windings;
nWindingsOuterX = windings;

optimizeROuter0 = rOuter;
optimizeDOuter0 = dOuter;
optimizationVariables0 = [optimizeROuter0; optimizeDOuter0];
%f = @(optimizationVariables)garrettCoil(rInner, dInner, optimizationVariables(1), optimizationVariables(2), nWindingsInnerZ, nWindingsInnerX, nWindingsOuterZ, nWindingsOuterX);

%x = fminsearch(f, optimizationVariables0);
%x = fminsearch(f, current);
% 
% [fieldSpread, bX, bY, bZ] = garrettCoil(rInner, dInner, rOuter, dOuter, nWindingsInnerZ, nWindingsInnerX, nWindingsOuterZ, nWindingsOuterX);
% B = sqrt(bX.^2 + bY.^2 + bZ.^2);
