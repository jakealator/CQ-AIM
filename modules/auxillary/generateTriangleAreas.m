% triangleAreas.m
% Created: Spring 2014 by jdr
% Last Modified: 01-31-2017 by jdr at home in Newark
%
% Input: mesh - a struct containing information about triangles in a finite
%               element mesh. 
%           N - Number of triangles in mesh struct.  
% Output: triangleAreas: Nx1 vector containing triangle areas ordered by
%                        mesh.nt. 
%
% Function to generate area of each triangle in a finite element mesh. 

function triangleAreas = generateTriangleAreas(mesh,N)

% Initialze triangleAreas vector 
triangleAreas = zeros(N,1);

% Rename triangles for convenience
nodesX = mesh.nodes(:,1);
nodesY = mesh.nodes(:,2);

% Order (x,y) coordinates so that each row of Tx/Ty corresponds to one
% triangle. 
Tx = nodesX(mesh.tris(:,1:3));
Ty = nodesY(mesh.tris(:,1:3));

% Define the area of each triangle with areaTri function. 
for i=1:mesh.nt
    triangleAreas(i) = areaTri(Tx(i,:),Ty(i,:));
end

end