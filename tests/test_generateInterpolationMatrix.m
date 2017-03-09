% test_generateInterpolationMatrix.m
% Created by: JDR on 02-24-2017 in Newark
% Last Modified: 03-07-2017 by JDR at home
%
% Tests the function ./modules/fftSolver/generateInterpolationMatrix.m

% Note: V = generateInterpolationMatrix(centers, rectangularElementsX, rectangularElementsY, midpointsX,midpointsY, triAreas, N)

%cd ../
%addpath(genpath('modules'))
%addpath(genpath('demo'))
% cd tests
mesh = 'circle';
meshStruct = initialize_mesh(mesh,1); % Initialize mesh (always use p=1 in second argument)
N=meshStruct.nt; % number of centroid points
M=2;
d=3;

[femStruct, farFieldStruct, iElements, jElements, multipoleMatrix, nearFieldDistances, P, flatP] = generateAuxillaryParams(meshStruct, N, M, d, @(x)1);

V = generateInterpolationMatrix(farFieldStruct.centers, ...
    farFieldStruct.rectangularElementsX, farFieldStruct.rectangularElementsY, ...
    femStruct.midpointsX, femStruct.midpointsY, femStruct.triAreas, femStruct.N);

