% generateFarFieldElements.m
% Created: 02-02-2017 by JDR at UD
% Last Modified: 
%
% Input: 
%
% Output: 
%
% Determines the artificial Cartesean support of each triangular element by
% finding the element of the rectangular grid closest to each centroid and
% adding the nearest grid points to that element. Also generates the
% interpolation matrix V which maps functions on a triangular element to a
% function on the new rectangular element. 

function [centers, rectangularElements, V] = generateFarFieldElements(centroids, N, farFieldGrid)

nG = length(farFieldGrid);
centers = zeros(N,1);
rectangularElements = zeros(N,9);

for j=1:N
    % Find the minimum L2-distance between each centroid and far field grid
    % element. This is probably a slower-than-needed method. 
    [centers(j),I] = min(sqrt((centroids(j,1)-farFieldGrid(:,1))^2+...
        (centroids(j,2)-farFieldGrid(:,2))^2));
    
    % Find the 9 closest grid points to centers(j). These are just the 4
    % closest components of farFieldGrid on either side. Need to be careful
    % not to overflow; since nG>9 unless something is wrong, we're missing
    % a case below (most likely edge cases won't matter anyway because
    % they'll be outside the domain). 
    minI = I-4;
    maxI = I+4;
    if (minI>1 && maxI<nG)
        rectangularElements(j,:) = farFieldGrid(minI:maxI,:);
    elseif minI<=1
        sizeElement = length(1:maxI);
        rectangularElements(j,1:sizeElement) = farFieldGrid(1:maxI,:);
    else
        sizeElement = length(minI:nG);
        rectangularElements(j,1:sizeElement) = farFieldGrid(minI:end,:);
    end
end



end