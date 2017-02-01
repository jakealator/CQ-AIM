% runAllTests.m
% Created: 02-01-2017 by JDR in Newark
% Last Modified: 
%
% Runs all testing functions to determine if there are any obvious bugs in
% the code. When writing individual tests, be sure to include a provision
% that tests if allTestsBool==1 so that auxillary parameters are not
% constantly being reinitialized. If there is a failure, the code will
% output an error and stop running tests. 

% Make sure all individual tests know this is a test of everything. 
allTestsBool=1;

% Add correct path to get meshes and files-to-test on path
cd ..
addpath(genpath('modules'))
addpath(genpath('demo'))
cd tests
mesh = 'twoCircles';
meshStruct = initialize_mesh(mesh,1); % Initialize mesh (always use p=1 in second argument)
N=meshStruct.nt; % number of centroid points


%---- Start running tests ----%
test_generateCentroids