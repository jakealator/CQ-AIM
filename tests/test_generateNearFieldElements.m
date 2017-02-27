% test_generateNearFieldElements.m
% Created by: JDR on 02-27-2017 at UD
% Last Modified:
%
% Tests the function ./modules/auxillary/generateNearFieldElements.m

% Note: [iElements, jElements, multipoleMatrix] = generateNearFieldElements(N,M,centers,hGrid,d)

% Add correct path to get meshes and files-to-test on path
cd ..
addpath(genpath('modules'))
addpath(genpath('demo'))
cd tests
mesh = 'twoCircles';
meshStruct = initialize_mesh(mesh,1); % Initialize mesh (always use p=1 in second argument)
N=meshStruct.nt; % number of centroid points
M=2;
% initialize centroids
centroids = generateCentroids(meshStruct, N);
[midpointsX,midpointsY] = generateMidpoints(meshStruct,N);
% spatial discretization parameter h
h = mesh_size(meshStruct);
%initialize triangle areas
triAreas = generateTriangleAreas(meshStruct, N);
%Generate a Cartesean grid containing D and a little extra with spacing
% h/2. First find largest/smallest centroids.
minX = min(centroids(:,1)); maxX = max(centroids(:,1));
minY = min(centroids(:,2)); maxY = max(centroids(:,2));
[ffX,ffY] = meshgrid(minX-M*h:h/M:maxX+M*h, minY-M*h:h/M:maxY+M*h);
ffX=ffX(:);ffY=ffY(:);
farFieldGrid = [ffX,ffY]; % output is size NGx2.
d=3; % Minimum near field range
[centers, rectangularElementsX, rectangularElementsY, V] = generateFarFieldElements(centroids, M, N, farFieldGrid, h, midpointsX,midpointsY, triAreas);

% %% Load variables from generateNearFieldElements.m
[iElements, jElements, multipoleMatrix] = generateNearFieldElements(N,M,centers,h,d);

% Test that distance between everything in the near field is in [0,d+M]h.
% To do this, test between each element in expansion square i and a random
% sorting of the elements in the jth expansion square. Everything should be
% within (d+M)h of each other. Note that this doesn't test this with 100%
% accuracy (there could be something weird with how the shifts occur such
% that all the elements are just within a lucky nearness of each other) but
% the probability of it passing incorrectly is absurdly small. 
nSparse = length(iElements);
k=0;
for i=1:nSparse
    randShift = randi((M+1)^2+1);
    currentDistance = abs(sqrt((rectangularElementsX(iElements(i),:)-circshift(rectangularElementsX(jElements(i),:),randShift)).^2+...
            (rectangularElementsY(iElements(i),:)-circshift(rectangularElementsY(jElements(i),:),randShift)).^2));
    if currentDistance < (d+M)*h
        k=k+1;
    else
        error('generateNearFieldElements test line 46: near field closeness test FAILED.')
    end
end
if k==nSparse
    sprintf('generateNearFieldElements test: near field closeness test: PASSED.')
end

    