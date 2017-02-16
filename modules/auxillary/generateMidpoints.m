% generateMidPoints.m
% Created by: JDR on 02-16-2017 in Newark
% Last Modified:
% Inputs: meshStruct - structure containing mesh information
%         N          - number of triangles
% Outputs: midpointsX/Y - Nx3 vector containing the midpoints of each triangle
%
% Computes the midpoint of each edge of a given triangle

function [midpointsX, midpointsY] = generateMidpoints(meshStruct,N)

for j=1:N
    midpointsX = (meshStruct.nodes(meshStruct.tris(j,1:3),1)+circshift(meshStruct.nodes(meshStruct.tris(j,1:3),1),1))/2;
    midpointsY = (meshStruct.nodes(meshStruct.tris(j,1:3),2)+circshift(meshStruct.nodes(meshStruct.tris(j,1:3),2),1))/2;
end    


end
