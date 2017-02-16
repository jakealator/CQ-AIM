% generateAuxillaryParams.m
% Created: 01-31-2017 by jdr at home in Newark
% Last Modified: 02-02-2017 by JDR at UD
%
% Input: meshStruct - a struct containing mesh information, including nodal
%                     locations (created by initialize_mesh.m)
%        N          - number of elements
%        M          - multipole expansion order (1 or 2)
%        
%
% Output: centroids          - vector of size Nx2 giving the centroid
%                              of each element. 
%         triAreas           - vector of size Nx1 giving area of each
%                              element. 
%         h                  - a scalar giving the spatial discretization size of
%                              the mesh .
%         i/jElements        - <Nx1 sparse matrices indexing locations of
%                              near field elements within an NxN matrix.
%         multipoleMatrix    - NxN sparse matrix with values <=M so that if (i,j)th
%                              element is 0<m then i is m elements away from
%                              j. If (i,j)=0, i is more than M elements from
%                              j. By convention, (i,i)=M. 
%         nearFieldDistances - NxN sparse matrix storing the distance
%                              between element i and j in location (i,j). 
%         farFieldGrid       - A Cartesean grid covering domain which is
%                              sufficiently dense (space = h/2).
%
% Generates auxillary parameters needed for frequency domain AIM. In
% particular, given a finite element mesh, code gives triangle centroids,
% triangle areas, and h-parameter. A far field grid is also generated;
% assume finite element grid is contained in a ball of radius r. The grid
% contains 4r. Spacing of grid is dictated by h and by the parameter M, the
% multipole expansion order

function [centroids, midpointsX, midpointsY, triAreas, h, iElements, jElements, multipoleMatrix, nearFieldDistances, farFieldGrid] = generateAuxillaryParams(meshStruct, N, M)


    % initialize centroids
    centroids = generateCentroids(meshStruct, N);
    [midpointsX,midpointsY] = generateMidpoints(meshStruct,N);
    
    %initialize triangle areas
    triAreas = generateTriangleAreas(meshStruct, N);
    
    % spatial discretization parameter h
    h = mesh_size(meshStruct);
    
    %intitialize multipole matrix
    [iElements, jElements, multipoleMatrix] = generateNearFieldElements(N,M,1,1);
    
    % Calculate distances between near field elements
    % !!!!Currently only works for M=1!!!!
    nearFieldDistances = generateNearFieldDistances(centroids, iElements, jElements, N,M);
    
    %Generate a Cartesean grid containing D and a little extra with spacing
    % h/2. First find largest/smallest centroids.
    minX = min(centroids(:,1)); maxX = max(centroids(:,1));
    minY = min(centroids(:,2)); maxY = max(centroids(:,2));
    [ffX,ffY] = meshgrid(minX-2*h:h/2:maxX+2*h, minY-2*h:h/2:maxY+2*h);
    farFieldGrid = [ffX(:),ffY(:)]; % output is size NGx2.

end