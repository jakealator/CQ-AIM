% test_generateNearFieldDistances.m
% Created: 02-01-2017 by JDR in Newark
% Last Modified:
%
% Tests the file generateNearFieldDistances by comparing sparse distance
% matrix with full distance matrix.  If
% allTestsBool==1 then multiple tests are running and we don't need to
% initialize any external parameters. 
%
% TODO: This actually fails when you compare against 1E-15 for a tolerance.
% Probably this isn't a problem, but it is very strange. 

if ~exist('allTestsBool','var')
        % Add correct path to get meshes and files-to-test on path
        cd ../
        addpath(genpath('modules'))
        addpath(genpath('demo'))
        cd tests
        mesh = 'twoCircles';
        meshStruct = initialize_mesh(mesh,1); % Initialize mesh (always use p=1 in second argument)
        N=meshStruct.nt; % number of centroid points
        M=1; % multipole expansion
end

% Generate centroids with generateCentroids.m (tested in
% test_generateCentroids.m)
centroids = generateCentroids(meshStruct,N);
% Generate multiple connectivity matrix 
multipoleMatrix = generateNearFieldElements(meshStruct,N,M);
% Generate sparse near field distance matrix with
% generateNearFieldDistances.
nearFieldDistances = full(generateNearFieldDistances(centroids, multipoleMatrix, N,M));

% Calculate full distance matrix for all elements (could be slow!)
D = sqrt(bsxfun(@plus,full(dot(centroids',centroids',1)),full(dot(centroids',centroids',1))')-full(2*(centroids*centroids')));
for i=1:N
    D(i,i) = 0;
    for j=1:N
        testDiff = D(i,j) - nearFieldDistances(i,j);
        if (abs(nearFieldDistances(i,j))>0 && abs(testDiff)>1E-14)
            error(strjoin({'generateNearFieldDistances Test: distance comparison at (', int2str(i),',',int2str(j),' FAILED'}))
        end
    end
end

