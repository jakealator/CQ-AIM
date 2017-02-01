% edgeLength.m
% Created: 13:35 02-20-2014 by jdr at UD
% Last Modified: 16:00 04-15-2014 by jdr at UD
% sideLengths = edgeLength(Mesh)
% Input: Mesh data structure
% Output: Mesh.nt x 3 array with the length of sides for each triangle in mesh
%
% Description: Uses node locations to create side lengths of each triangle
% in mesh. 

function sideLengths = edgeLength(Mesh)

sideLengths = [sqrt((Mesh.nodes(Mesh.tris(:,1),1)-...
    Mesh.nodes(Mesh.tris(:,2),1)).^2+(Mesh.nodes(Mesh.tris(:,1),2)-...
    Mesh.nodes(Mesh.tris(:,2),2)).^2), sqrt((Mesh.nodes(Mesh.tris(:,2),1)-...
    Mesh.nodes(Mesh.tris(:,3),1)).^2+(Mesh.nodes(Mesh.tris(:,2),2)-...
    Mesh.nodes(Mesh.tris(:,3),2)).^2), sqrt((Mesh.nodes(Mesh.tris(:,1),1)-...
    Mesh.nodes(Mesh.tris(:,3),1)).^2+(Mesh.nodes(Mesh.tris(:,1),2)-...
    Mesh.nodes(Mesh.tris(:,3),2)).^2)];

end
