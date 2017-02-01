% test_generateTriangleAreas.m
% Created: 02-01-2017 by JDR in Newark
% Last Modified: 
%
% Tests that the program generateTriangleAreas.m works correctly by checking
% if a few specially chosen triangle areas are computed correctly.

% Need to create a struct with the appropriate information in order to
% test. 

% Three triangles at: (0,0)--(1,0)--(0,1) [area = 1/2]
%                     (1,0)--(0,1)--(0,-1) [area = 1]
%                     (0,-1)--(0,0)--(-1,0) [area=1/2]
% This tests that order doesn't matter, that negative numbers aren't
% treated oddly, and that formula works for non-right triangles. 
trisVals = [[1,2,3];[2,3,4];[4,1,5]];
nodeVals = [[0,0];[1,0];[0,1];[0,-1];[-1,0]];
nt = 3;

meshStruct = struct('nt',nt,'tris',trisVals,'nodes',nodeVals);

triangleAreas = generateTriangleAreas(meshStruct,3);

% Make sure triangle areas are calculated correctly for some known
% triangles. 
if (abs(triangleAreas-[0.5;1;0.5]))>1E-15
    error('generateTriangleAreas Tests Line 25: Calculation of areas FAILED')
else
    sprintf('generateTriangleAreas Tests: Calculation of areas PASSED')
end
