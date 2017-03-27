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

function [centers, rectangularElementsX, rectangularElementsY, rectangularLocations, P,flatP] = generateFarFieldElements(centroids, M, N, farFieldGrid, h, midpointsX,midpointsY, triAreas)

centers = zeros(N,2);
rectangularElementsX = zeros(N,9);
rectangularElementsY = zeros(N,9);
rectangularLocations = zeros(N,9);
nG = length(farFieldGrid);

for j=1:N
    
    % Find the minimum L2-distance between each centroid and far field grid
    % element. This is probably a slower-than-needed method. 
    [~,I] = min(sqrt((centroids(j,1)-farFieldGrid(:,1)).^2+...
        (centroids(j,2)-farFieldGrid(:,2)).^2));
    centers(j,:) = [farFieldGrid(I,1),farFieldGrid(I,2)];
        
    % Find the 9 closest grid points to centers(j). These are just the 8
    % closest components of farFieldGrid on either side. Shouldn't need 
    % to worry about edge cases because there's a built-in buffer in the
    % grid. 
    [rectangularElementsXSquare,rectangularElementsYSquare] = meshgrid(centers(j,1)-h/M:h/M:centers(j,1)+h/M,...
        centers(j,2)-h/M:h/M:centers(j,2)+h/M);
    rectangularElementsX(j,:)=rectangularElementsXSquare(:); 
    rectangularElementsY(j,:)=rectangularElementsYSquare(:);
    % Calcualte the entries of rectangularElementsX/Y in farFieldGrid. To
    % do this, first find where they would be in a nGxnG matrix of farField
    %Grid. Then convert to the linear index. 
    [Ix,Iy]=ind2sub([sqrt(nG),sqrt(nG)],I); 
    [IxSquare,IySquare] = meshgrid(Ix-1:1:Ix+1,Iy-1:1:Iy+1);
    IxSquare=IxSquare'; IySquare=IySquare';
    Ix=IxSquare(:);Iy=IySquare(:);
    rectangularLocations(j,:) = sub2ind([sqrt(nG),sqrt(nG)],Ix,Iy);

    if ~((abs(rectangularElementsX(j,:)'-farFieldGrid(rectangularLocations(j,:),1))<1E-15)&(abs(rectangularElementsY(j,:)'-farFieldGrid(rectangularLocations(j,:),2))<1E-15))
        sprintf(j)
        error('far field location wrong!')
    end
        
    
  
end

% Using the above, create the interpolation matrix P. The function generate
% InerpolationMatrix is Nx9. We want to unflatten it according to
% rectangularLocations. 
flatP = generateInterpolationMatrix(centers, rectangularElementsX, rectangularElementsY, midpointsX,midpointsY, triAreas, N);

femIndex = repmat(1:N,9,1);
femIndex = femIndex(:);
rectangularLocationsSideways=rectangularLocations.';
rectangularLocationsSideways = rectangularLocationsSideways(:);
flatPSideways = flatP.';
P = sparse(rectangularLocationsSideways,femIndex,flatPSideways(:),nG,N);


end
