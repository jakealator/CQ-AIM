% generateFarFieldElements.m
% Created: 02-02-2017 by JDR at UD
% Last Modified: 02-03-2017 by JDR in Newark
%
% Input: centroids    - Nx2 vector finite Element mesh centroids
%        N            - Number of elements 
%        farFieldGrid - Ngx2 vector containing the rectangular grid points
%
% Output: centers             - Nx2 vector of the centers of each 
%                               expansion box, ordered
%                               in the same way as the FEM elements 
%         rectangularElements - Nx9 vector points in each expansion box,
%                               ordered in the same way as centers
%         V                   - nGxN matrix mapping from finite element
%                               mesh to rectangular mesh
%
% Determines the artificial Cartesian support of each triangular element by
% finding the element of the rectangular grid closest to each centroid and
% adding the nearest grid points to that element. Also generates the
% interpolation matrix V which maps functions on a triangular element to a
% function on the new rectangular element. 

function [centers, rectangularElementsX, rectangularElementsY, V] = generateFarFieldElements(centroids, N, farFieldGrid)

nG = length(farFieldGrid);
centers = zeros(N,2);
rectangularElementsX = zeros(N,9);
rectangularElementsY = zeros(N,9);

for j=1:N
    % Find the minimum L2-distance between each centroid and far field grid
    % element. This is probably a slower-than-needed method. 
    [~,I] = min(sqrt((centroids(j,1)-farFieldGrid(:,1)).^2+...
        (centroids(j,2)-farFieldGrid(:,2)).^2));
    centers(j,:) = [farFieldGrid(I,1),farFieldGrid(I,2)];
    
    % Find the 9 closest grid points to centers(j). These are just the 4
    % closest components of farFieldGrid on either side. Need to be careful
    % not to overflow; since nG>9 unless something is wrong, we're missing
    % a case below (most likely edge cases won't matter anyway because
    % they'll be outside the domain). 
    minI = I-4;
    maxI = I+4;
    if (minI>1 && maxI<nG)
        rectangularElementsX(j,:) = farFieldGrid(minI:maxI,1);
        rectangularElementsY(j,:) = farFieldGrid(minI:maxI,2);
    elseif minI<=1
        sizeElement = length(1:maxI);
        rectangularElementsX(j,1:sizeElement) = farFieldGrid(1:maxI,1);
        rectangularElementsY(j,1:sizeElement) = farFieldGrid(1:maxI,2);
    else
        sizeElement = length(minI:nG);
        rectangularElementsX(j,1:sizeElement) = farFieldGrid(minI:end,1);
        rectangularElementsY(j,1:sizeElement) = farFieldGrid(minI:end,2);
    end
end

% Using the above, create the interpolation matrix V
%%%%% This is part of the reason this only works for multipole M=2
V = generateInterpolationMatrix(centers, rectangularElementsX, rectangularElementsY, midpointsX,midpointsY,N);



end