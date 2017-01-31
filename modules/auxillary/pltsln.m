% pltsln.m
% Created: 00:10 03-21-2014 by jdr at UD
% Last Modified: 16:00 04-15-2014 by jdr at UD
% pltsln(Mesh,p,x)
% Input: Mesh, a finite element mesh
%        p, order of element
%        x, solution vector to FEM problem
% Output: Generates image of mesh
% 
% Draws finite element solution using pth order interpolation to on given 
% mesh

function pltsln(Mesh,xCentroid,x)

 xInterp = griddata(xCentroid(:,1),xCentroid(:,2),x,Mesh.nodes(:,1),Mesh.nodes(:,2));
 h = trisurf(Mesh.tris(:,1:3),Mesh.nodes(:,1),Mesh.nodes(:,2),xInterp);
 
 set(h, 'FaceColor','interp')
 shading interp
 colorbar

end