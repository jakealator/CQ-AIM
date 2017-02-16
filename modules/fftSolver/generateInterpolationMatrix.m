% generateInterpolationMatrix.m
% Created: 02-16-2017 by JDR in Newark
% Last Modified:
%
% Inputs: centroids           -
%         centers             -
%         rectangularElements -
%         meshStruct          -
% 
% Output: V - Nx9 matrix containing interpolation coefficients between
% Galerkin and delta functions
%
% Finds coefficient matrix to interpolate between p0 Galerkin basis 
% functions on FEM triangulation and "delta basis functions" on Cartesian
% discretization. 

function V = generateInterpolationMatrix(centers, rectangularElements, midpointsX,midpointsY, triAreas, N)

% Order of m-vector containing exponents
mVec = [[0,0];[0,1];[0,2];[1,0];[1,1];[1,2];[2,0];[2,1];[2,2]];

% Use quadrature against bFun to calculate rhs
bFun = @(x,y,m,midpointsX,midpointsY)((midpointsX-x)^m(1)*(midpointsY-y)^m(2));


% Compute V for each element/expansion cube
for j=1:N
    W = zeros(9,9);
    b = zeros(9,1);
    
    % Use midpoint rule to compute righ-hand-side. This is exact for
    % polynomials up to 2nd order so is fine for M=2
    % WAIT: Does this incur error? I think it does. 
    for i=1:9
        b(i) = (bFun(centers(j,1),centers(j,2),mVec(i,:),midpointsX(j,1),midpointsY(j,1))+...
            bFun(centers(j,1),centers(j,2),mVec(i,:),midpointsX(j,2),midpointsY(j,2))+...
            bFun(centers(j,1),centers(j,2),mVec(i,:),midpointsX(j,3),midpointsY(j,3)));
    end
    % Need |K|/3 for midpoint quadrature rule
    b = triAreas(j)/3*b;



end