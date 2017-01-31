function [x,x1,x2,x3] = generateCentroids(mesh,N)

% Find barycentric centers of triangles and barycentric centers of the
% three subtriangles defined by using the barycentric center of the
% original triangle and two of the original triangle vertices. 
x = zeros(N,2);
x1 = zeros(N,2);
x2 = zeros(N,2);
x3 = zeros(N,2);
% Generate centroid points. x are the centroids of the original triangle.
% x#c are the centroids of the subtriangles (for singular part of integral)
% Also generate interior points needed for 1st order quadrature x1,x2,x3. 
for i=1:N
     x(i,1) = 1/3*sum(mesh.nodes(mesh.tris(i,1:3),1));
     x(i,2) = 1/3*sum(mesh.nodes(mesh.tris(i,1:3),2));
     x1(i,1) = 1/6*(4*mesh.nodes(mesh.tris(i,1),1)+sum(mesh.nodes(mesh.tris(i,2:3),1)));
     x1(i,2) = 1/6*(4*mesh.nodes(mesh.tris(i,1),2)+sum(mesh.nodes(mesh.tris(i,2:3),2)));
     x2(i,1) =  1/6*(4*mesh.nodes(mesh.tris(i,1),1)+sum(mesh.nodes(mesh.tris(i,[1,3]),1)));
     x2(i,2) = 1/6*(4*mesh.nodes(mesh.tris(i,1),2)+sum(mesh.nodes(mesh.tris(i,[1,3]),2)));
     x3(i,1) = 1/6*(4*mesh.nodes(mesh.tris(i,3),1)+sum(mesh.nodes(mesh.tris(i,1:2),1)));
     x3(i,2) = 1/6*(4*mesh.nodes(mesh.tris(i,3),2)+sum(mesh.nodes(mesh.tris(i,1:2),2)));
end

end