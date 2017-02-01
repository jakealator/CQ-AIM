% mesh_size.m
% Created: 13:35 02-20-2014 by jdr at UD
% Last Modified: 16:00 04-15-2014 by jdr at UD
% h = mesh_size(Mesh)
% Input: Mesh data structure
% Output: Mesh parameter, h
%
% Description: Calculates the radius of the smallest circumscribed circle
% around each triangle in Mesh. 

function h = mesh_size(meshStruct)

 % Calculate Mesh.nt x 3 array of side lengths for each triangle in mesh. 
 sideLength = edgeLength(meshStruct);
 % Using side lengths, the h parameter can be constructed. 
 hK = triangleSizeParameter(sideLength);
 
 % h is defined as the maximum h over all the elements.
 h = max(hK); 

end
