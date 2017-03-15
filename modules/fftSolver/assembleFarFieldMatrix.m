% assembleFarFieldMatrix.m
% Created on 02-03-2017 by JDR in Newark
% Last Modified:
% 
% Inputs: farFieldGrid - nGx2 vector of rectangular mesh points
%         P            - Nx(M+1)^2 interpolation matrix from triangular mesh to
%                        rectangular mesh 
%         N            - Number of finite elements
%         rectangularElements - the rectangles in the near field
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

function farFieldElements = assembleFarFieldMatrix(waveNumber,  flatP, N,...
    rectangularElementsX, rectangularElementsY, iElements, jElements,...
    rectangularLocations, GD)

kMax = length(nonzeros(iElements));

% The fundamental solution for 2D. Note that this will output as something
% the size of its input and if the input is zero, the output will be zero
% up to machine precision. This is to avoid any issues with near field NaNs
% not canceling out below. 


sElements = zeros(kMax,1);
% Question: Can this be done without the loop? This is about half my total
% time!!!
for j=1:kMax
    % build distance matrix between each expansion box element.
    %expBoxI = [rectangularElementsX(iElements(j),:)',rectangularElementsY(iElements(j),:)'];
    %expBoxJ = [rectangularElementsX(jElements(j),:)',rectangularElementsY(jElements(j),:)'];
    %D = sqrt(bsxfun(@plus,full(dot(expBoxI',expBoxI',1)),full(dot(expBoxJ',expBoxJ',1))')-full(2*(expBoxJ*expBoxI')));
    GDCurrent = GD(rectangularLocations(iElements(j),:), rectangularLocations(jElements(j),:));
    sElements(j) = flatP(iElements(j),:)*...
        (GDCurrent*flatP(jElements(j),:).');
end

farFieldElements = sparse(iElements, jElements, sElements, N,N);

end