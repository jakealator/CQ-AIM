% test_assembleNearFieldMatrices.m
% Created: 02-01-2017 by JDR at UD
% Last Modified: 
%
% Test script to test if assembleNearFieldMatrices.m is correct.
% Specifically, compares algorithm resulting in a sparse matrix with
% comparable entries in the full agortihm. 
%
% NOTE: This currently only works in MATLAB 2016b and up!

if ~exist('allTestsBool','var') || (exists('allTestsBool','var') && allTestsBool~=1)
    % Add correct path to get meshes and files-to-test on path
    cd ../
    addpath(genpath('modules'))
    addpath(genpath('demo'))
    cd tests
    mesh = 'twoCircles';
    meshStruct = initialize_mesh(mesh,1); % Initialize mesh (always use p=1 in second argument)
    N=meshStruct.nt; % number of centroid points
    M=1; % multipole expansion
    c0 = 1; % Speed of sound in free space (default=1)
    c=@(x,y)(sqrt(2)); % Speed of sound in inhomogeneity (should be either 0<c<1 or 1<c<2 in general)
    s = 1; % Frequency
    % Generate centroids with generateCentroids.m (tested in
    % test_generateCentroids.m)
    centroids = generateCentroids(meshStruct,N);
    % Generate multiple connectivity matrix 
    [iElements, jElements, multipoleMatrix] = generateNearFieldElements(meshStruct,N,M);
    % Generate sparse near field distance matrix with
    % generateNearFieldDistances.
    nearFieldDistances = full(generateNearFieldDistances(centroids, iElements, jElements, N,M));
    % Calculate triangle areas
    triAreas = generateTriangleAreas(meshStruct, N);
end

% Compute the sparse version of mass and stiffness matrices. 
 [KSparse,MSparse,~] = assembleNearFieldMatrices(triAreas,nearFieldDistances, iElements, jElements, centroids, c,c0,s,N);
% Compute the full version of mass and stiffness matrices with local
% function defined below. 
 [KFull, MFull] = assembleFullNearFieldMatrices(centroids, triAreas, c, c0, s, N);
 
 % Test mass matrix equality
 if ~(isempty(nonzeros(MSparse-MFull)))
     error('assembleNearFieldMatrices Test: Mass matrix M is computed differently for sparse and full matrices. FAILED')
 else
     sprintf('assembleNearFieldMatrices Test: Mass matrix M computation PASSED')
 end
 
 % Test stiffness matrix equality
 if abs((max(max(KSparse-KFull.*multipoleMatrix))))>1E-14
     error('assembleNearFieldMatrices Test: Stiffness matrix K is computed differently for sparse and full matrices. FAILED')
 else
     sprintf('assembleNearFieldMatrices Test: Stiffness matrix K computation PASSED')
 end

     
 
 
    