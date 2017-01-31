% drwmsh.m
% Created: 15:15 02-16-2014 by jdr at UD
% Last Modified: 16:00 04-15-2014 by jdr at UD
% drwmsh(Mesh)
% Input: Mesh data structure
% Output: Generates image of mesh
% 
% Draws finite element mesh in Matlab. 

function drwmsh(Mesh)
 
 figure 
 
 % Give an epsilon's distance around figure
 xmax = max(Mesh.nodes(:,1));
 xmin = min(Mesh.nodes(:,1));
 ymax = max(Mesh.nodes(:,2));
 ymin = min(Mesh.nodes(:,2));
 eps = abs(max(1/10*(xmax-xmin),1/10*(ymax-ymin))); %must be positive
 axis([xmin-eps,xmax+eps,ymin-eps,ymax+eps])
 
 % xData and yData are 3 x Mesh.nt arrays which contain the locations 
 % of the nodes uses to generate triangles.  
 xData = [Mesh.nodes(Mesh.tris(:,1),1)';Mesh.nodes(Mesh.tris(:,2),1)'...
     ;Mesh.nodes(Mesh.tris(:,3),1)'];
 yData = [Mesh.nodes(Mesh.tris(:,1),2)';Mesh.nodes(Mesh.tris(:,2),2)'...
     ;Mesh.nodes(Mesh.tris(:,3),2)'];
  
  % Plot triangle faces
 patch(xData,yData,'w')

end
