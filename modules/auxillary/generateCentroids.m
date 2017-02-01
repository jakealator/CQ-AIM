% generateCentroids.m
% Created: Spring 2014 by jdr
% Last modified: 02-01-2017 by JDR in Newark
%
% Input: meshStruct - a struct containing mesh information, given by
%                     initialize_mesh.m. 
%        N          - number of elements
% Output: x         - Nx2 vector containing centroids. 
%
% Calculates centroids of elements in a mesh given the node information
% about the mesh. 

function x = generateCentroids(meshStruct,N)

% Find barycentric centers of triangles
x = zeros(N,2);
% Generate centroid points, the 'center' point of each triangle. Order is
% the same as mesh.tris. Note that algorithm is a little slow wrt using a
% loop, but meshStruct data structure seems to require this. 
for j=1:N
    x(j,1) = 1/3*sum(meshStruct.nodes(meshStruct.tris(j,1:3),1));
    x(j,2) = 1/3*sum(meshStruct.nodes(meshStruct.tris(j,1:3),2));
end

end