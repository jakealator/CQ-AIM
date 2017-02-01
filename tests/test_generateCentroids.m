% test_generateCentroids.m
% Created: 02-01-2017 by JDR in Newark
% Last Modified: 
%
% Tests that the program generateCentroids.m works correctly by checking
% the number of centroids and their values against a specific case. If
% allTestsBool==1 then multiple tests are running and we don't need to
% initialize any external parameters. 

if allTestsBool ~= 1
    % Add correct path to get meshes and files-to-test on path
    cd ../
    addpath(genpath('modules'))
    addpath(genpath('demo'))
    cd tests
    mesh = 'twoCircles';
    meshStruct = initialize_mesh(mesh,1); % Initialize mesh (always use p=1 in second argument)
    N=meshStruct.nt; % number of centroid points
end

% Generate centroids with generateCentroids.m
centroids = generateCentroids(meshStruct,N);

%---- Begin tests ----%
centroidPassedTests = 0; % Indicates number of passed tests. 

% Test that the number of centroids is correct
if (length(centroids) == N)
    sprintf('generateCentroids Test: number of centroids PASSED')
    centroidPassedTests = centroidPassedTests + 1;
else
    error('generateCentroids Test Line 32: number of centroids FAILED.')
end

% Test that a handful of centroids were computed correctly. 
rElements = randi(N,5,1); % Pick 5 random elements and compute centroids 'by hand'.
for j=1:5
    cElement = rElements(j);                % select element
    cTris = meshStruct.tris(cElement,1:3);  % select triangle nodes for element
    cNodes = meshStruct.nodes(cTris,:);  % find node locations for each triangle node
    cCentroidX = 1/3*(cNodes(1,1)+cNodes(2,1)+cNodes(3,1)); % x-centroid component
    cCentroidY = 1/3*(cNodes(1,2)+cNodes(2,2)+cNodes(3,2)); % y-centroid component
    % If there is a mistake, kill the computation 
    if ~(cCentroidX==centroids(cElement,1) && cCentroidY==centroids(cElement,2))
        error(strjoin({'generateCentroids Test Line 45: centroid number ',int2str(cElement), 'computation FAILED.'}))
    end
end
% If everything was fine, let the user know. 
sprintf('generateCentroids Test: centroid computations PASSED')
centroidPassedTests = centroidPassedTests + 1;

if centroidPassedTests == 2
    sprintf('generateCentroids Tests: ALL PASSED')
end 

%---- End Tests ----%

