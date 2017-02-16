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

function farFieldElements = assembleFarFieldMatrix(farFieldGrid,  V, N, rectangularElements, iElements, jElements)

kMax = length(nonzeros(iElements));

% This breaks for overlapping elements right now!!!
for j=1:kMax
    % build distance matrix between each expansion box element.
    expBoxI = rectangularElements(iElements(j),:);
    expBoxJ = rectangularElements(jElements(j),:);
    D = sqrt(bsxfun(@plus,full(dot(expBoxI',expBoxI',1)),full(dot(expBoxJ',expBoxJ',1))')-full(2*(expBoxJ*expBoxI')));
    sElements(j) = dot(V(iElements(j),:),...
        fundamentalSolution(D)*V(jElements(j),:).');
end

farFieldElements = sparse(iElements, jElements, sElements, N,N);

end