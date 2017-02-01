% generateAuxillaryParams.m
% Created: 01-31-2016 by jdr at home in Newark
% Last Modified:
%
% Input: meshStruct - a struct containing mesh information, including nodal
%                     locations (created by initialize_mesh.m)
%        N          - number of elements
%        M          - multipole expansion order (1 or 2)
%        
%
% Output: centroids    - vector of size Nx2 giving the centroid
%                        of each element. 
%         triAreas     - vector of size Nx1 giving area of each
%                        element. 
%         h            - a scalar giving the spatial discretization size of
%                        the mesh .
%         farFieldGrid - 
%
% Generates auxillary parameters needed for frequency domain AIM. In
% particular, given a finite element mesh, code gives triangle centroids,
% triangle areas, and h-parameter. A far field grid is also generated;
% assume finite element grid is contained in a ball of radius r. The grid
% contains 4r. Spacing of grid is dictated by h and by the parameter M, the
% multipole expansion order

function [centroids, triAreas, h] = generateAuxillaryParams(meshStruct, N, M)


    % initialize centroids
    centroids = generateCentroids(meshStruct, N);
    
    %initialize triangle areas
    triAreas = generateTriangleAreas(meshStruct, N);
    
    % spatial discretization parameter h
    h = mesh_size(meshStruct);


end