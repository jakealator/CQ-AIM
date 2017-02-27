% test_generateAuxillaryParameters.m
% Created by: JDR on 02-24-2017 in Newark
% Last Modified:
%
% Tests the function ./modules/generateAuxillaryParameters.m

% Note: [centroids, midpointsX, midpointsY, triAreas, h, iElements, jElements, multipoleMatrix, nearFieldDistances, farFieldGrid] = generateAuxillaryParams(meshStruct, N, M)

% Add correct path to get meshes and files-to-test on path
cd ..
addpath(genpath('modules'))
addpath(genpath('demo'))
cd tests
mesh = 'twoCircles';
meshStruct = initialize_mesh(mesh,1); % Initialize mesh (always use p=1 in second argument)
N=meshStruct.nt; % number of centroid points
M=1;

[femStruct, farFieldStruct, iElements, jElements, multipoleMatrix, nearFieldDistances, V] = generateAuxillaryParams(meshStruct, N, M);
