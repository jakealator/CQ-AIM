% test_assembleNearFieldMatrices.m
% Created: 02-01-2017 by JDR at UD
% Last Modified: 
%
% Test script to test if assembleNearFieldMatrices.m is correct.
% Specifically, compares algorithm resulting in a sparse matrix with
% comparable entries in the full agortihm. 

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

     
 
 
 % Function to compute full matrices (slow and requires R2016b or above to
 % work!). Taken from previously-written Galerkin code. Oddness of the
 % following mostly due to integration, which is done by approximating
 % triangles by circles of the same size and integrating exactly over
 % those. See Chapter 5 of my thesis. 
 function [K,M] = assembleFullNearFieldMatrices(centroids, triAreas, c, c0, s, N)
 
% Compute distance between all elements (diagonal element not important)
D = sqrt(bsxfun(@plus,full(dot(centroids',centroids',1)),full(dot(centroids',centroids',1))')-full(2*(centroids*centroids')));
% Evaluate Green's function at distance matrix
H01 = besselh(0,1,1i*s/c0*D);
% Evaluate H_1^1 at appropriate vector to make exact integration over a
% circle easier. 
H11 = besselh(1,1,1i*s/c0*sqrt(triAreas/pi));
% Compute index of refraction at centroids
qj =((c0./c(centroids(:,1),centroids(:,2))).^2-ones(N,1));
% Allocate memory for stiffness matrix 
K = zeros(N,N);

% On diagonal, stiffness matrix is computed differently from off diagonal
K(1:N+1:N^2) = -triAreas.*(4*1i*c0^2/s^2+(2*pi*1i*c0/s).*H11);

% Compute off-diagonal elements 
for i = 1:N
  for j=1:N
      if i~=j
          K(i,j) = triAreas(j)*H01(i,j); % Add in further contribution from triangle areas below. 
      end
  end
  K(i,:) = triAreas(i)*K(i,:);
end

% Create a sparse matrix whose diagonal elements are exactly the ones
% created above. 
M = spdiags(triAreas./qj,0,N,N);
K = (1i*s^2)/(4*c0^2)*K;
 
 
 
 end
    
    