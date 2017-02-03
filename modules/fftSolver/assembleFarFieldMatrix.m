% assembleFarFieldMatrix.m
% Created on 02-03-2017 by JDR in Newark
% Last Modified:
% 
% Inputs: farFieldGrid - nGx2 vector of rectangular mesh points
%         V            - nGxN interpolation matrix from triangular mesh to
%                        rectangular mesh
%         N            - Number of finite elements
%         i/jElements  - Locations of near field elements
%
% Outputs: farFieldElements - NxN sparse matrix of far field matrix components
%
% Generates the elements of the far field matrix corresponding to given
% points (i,j). In particular, uses Afar = V G V^T, where V is the
% interpolation matrix from a finite element mesh to a rectangular mesh and
% G is the Green function evaluated at each of the rectangular mesh points.
% Only generates elements which are needed (to cancel out equivilent
% components in near field). 

function farFieldElements = assembleFarFieldMatrix(farFieldGrid,  V, N, iElements, jElements)

kMax = length(nonzeros(iElements));

farFieldElements = sparse(iElements, jElements, sElements, N,N);

end