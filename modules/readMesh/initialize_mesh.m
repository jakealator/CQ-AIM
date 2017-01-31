% initialize_mesh.m
% Created: 12:15 02-16-2014 by jdr at UD
% Last Modified: 20:00 04-16-2014 by jdr at UD
% Input: String which contains location of Gmsh mesh file. 
%        p, degree of Lagrange DoF
% Output: Mesh data structure which contains the following: 
%            Mesh.nn, number of DoFs in mesh;
%            Mesh.nodes, Mesh.nn x 2 array with coordinates of nodes;
%            Mesh.ne, number of edges in mesh;
%            Mesh.edges, Mesh.ne x 3 array whose entries are (in order)
%              first vertex on jth edge, second vertex on jth edge, and 
%              physical reference number of edge (0 if edge not on a
%              labeled boundary);
%            Mesh.nt, number of triangles in mesh; and 
%            Mesh.tris, Mesh.nt x (#DoFs+1) array whose entries are (in order)
%              3 vertex numbers of the jth triangle in increasing order,
%              followed by the other DoFs in any order, followed by the
%              phsical reference number for the domain of the jth
%              triangle.
% 
% File reads in a Gmesh mesh file and returns a data structure which can be
% used in Matlab applications. Uses P. Monk's Gread.m to read the Gmsh
% file. 


function Mesh=initialize_mesh(name,p)

% Use P. Monk's Gread.m file to do the brunt of the work. 
 meshFull = Gread(name);

 MeshPre = struct('nn',meshFull.NoN,'nodes',meshFull.nodes(:,1:2),'ne',meshFull.Ne...
     ,'edges',[sort(meshFull.edge')',meshFull.edgetag(:,1)],'nt',meshFull.Nf,...
     'tris',[sort(meshFull.faces')',meshFull.facetag(:,1)]);
 
 % Gread only returns the boundary edges (not the edges between triangles).
 % Manually add in triangle edges here. There are duplicates, because a
 % triangle with edge a->b abuts a triangle with edges b->a. "unique" gets
 % rid of that when the duplication is within triangles, rather than
 % between a boundary and a triangle. 'rows' selects unique rows (rather
 % than columns) and 'stable' leaves vector unsorted. 
  uniqueEdges = unique(vertcat(MeshPre.edges(:,1:2), MeshPre.tris(:,[1,3]), ...
     MeshPre.tris(:,2:3), MeshPre.tris(:,1:2)),'rows','stable');

 % The first Mesh.ne edges have physical tag 1 (note that this does not
 % take into account both Dirichlet and Neumann tags, for example. Gread
 % does not seem to have this capability). The sortrows command sorts
 % things in the same order as P. Monk - first by physical tag in
 % decreasing order, then increasing in the second column and finally in the
 % first column. 
 MeshPre.edges = sortrows(vertcat([uniqueEdges(1:MeshPre.ne,:), ones(MeshPre.ne,1)], ...
     [uniqueEdges(MeshPre.ne+1:end,:), zeros(length(uniqueEdges)-MeshPre.ne,1)]),...
     -3);

 
 MeshPre.ne = length(MeshPre.edges);
 if p==1
     Mesh = MeshPre;
 elseif p==2
 % Calculate edge midpoints and get rid of duplicates. Specifically, after
 % the already-calculated vertex points, add the aij=1/2*(ai+aj) elements
 % for each triangle in the mesh. In order, for each triangle, add them in
 % the order a12,a23,a13. 
 newNodes = zeros(3*MeshPre.nt,2);
 newNodes(1:3:3*MeshPre.nt-2,:)=1/2*(MeshPre.nodes(MeshPre.tris(:,1),:)+...
     MeshPre.nodes(MeshPre.tris(:,2),:));
 newNodes(2:3:3*MeshPre.nt-1,:)=1/2*(MeshPre.nodes(MeshPre.tris(:,2),:)+...
     MeshPre.nodes(MeshPre.tris(:,3),:));
  newNodes(3:3:3*MeshPre.nt,:)=1/2*(MeshPre.nodes(MeshPre.tris(:,1),:)+...
     MeshPre.nodes(MeshPre.tris(:,3),:));
 [uniqueNodes,~,IC]=unique(newNodes,'rows','stable'); 

 % Include newly-calculated DoFs Mesh.tris
 tris=[MeshPre.tris(:,1),MeshPre.tris(:,2),MeshPre.tris(:,3),...
         MeshPre.nn+IC(1:3:3*MeshPre.nt-2),MeshPre.nn+IC(2:3:3*MeshPre.nt-1),...
         MeshPre.nn+IC(3:3:3*MeshPre.nt),MeshPre.tris(:,4)];
     
 % Find which new nodes are on boundary of mesh.
 edgeIndex = find(MeshPre.edges(:,3)==1);
 newEdge = [];
 for i=1:length(edgeIndex)
     % Find where in triangles data structure the boundary nodes are
     [bool1,loc1] = ismember(MeshPre.edges(edgeIndex(i),1:2),tris(:,1:2),'rows');
     [bool2,loc2] = ismember(MeshPre.edges(edgeIndex(i),1:2),tris(:,2:3),'rows');
     [bool3,loc3] = ismember(MeshPre.edges(edgeIndex(i),1:2),tris(:,[1,3]),'rows');
     % Add the appropriate new node to the edge structure. Make use of
     % consistency in how triangle data structure was set up. 
     if bool1
         newEdge = vertcat(newEdge,[MeshPre.edges(edgeIndex(i),1),tris(loc1,4),1]);
         newEdge = vertcat(newEdge,[tris(loc1,4),MeshPre.edges(edgeIndex(i),2),1]);
     elseif bool2
         newEdge = vertcat(newEdge,[MeshPre.edges(edgeIndex(i),1),tris(loc2,5),1]);
         newEdge = vertcat(newEdge,[tris(loc2,5),MeshPre.edges(edgeIndex(i),2),1]);
     elseif bool3
         newEdge = vertcat(newEdge,[MeshPre.edges(edgeIndex(i),1),tris(loc3,6),1]);
         newEdge = vertcat(newEdge,[tris(loc3,6),MeshPre.edges(edgeIndex(i),2),1]);
     else
         error('Edge not found in Mesh.tris. Something is wrong.')
     end
 end
         
 nodes=[MeshPre.nodes;uniqueNodes];
 nn = length(nodes);

 % Add new edges to beginning of edges data structure. Sort in the same way
 % as above. 
edges = sortrows(vertcat(newEdge,MeshPre.edges(length(edgeIndex)+1:end,:)),-3);
 
 Mesh = struct('nn',nn,'nodes',nodes,'ne',MeshPre.ne,...
      'edges',edges,'nt',length(tris),'tris',tris);
 else
     error('Sorry, p>2 not implemented.')
 end
end