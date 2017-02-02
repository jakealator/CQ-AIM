% assembleFullNearFieldMatrices.m
% Created: by JDR at home in Newark
% Last Modified: 
%
% Input: centroids - Nx2 vector of element centroids
%        triAreas  - Nx1 vector of element areas
%        c         - anonymous function containing speed of sound in medium
%        c0        - speed of sound in free space (set to 1 unless you're
%                    brave)
%        s         - wavenumber
%        N         - number of elements
% Output: K - NxN matrix containing 'stiffness' matrix
%         M - NxN matrix containing mass matrix
%
% Function to compute full matrices Taken from previously-written 
% Galerkin code. Oddness of the following mostly due to integration, which
% is done by approximating triangles by circles of the same size and 
% integrating exactly over those. See Chapter 5 of my thesis. 

function [K,M] = assembleFullNearFieldMatrices(centroids, triAreas, c, c0, s, N)

% Compute distance between all elements (diagonal element not important)
D = sqrt(bsxfun(@plus,full(dot(centroids',centroids',1)),full(dot(centroids',centroids',1))')-full(2*(centroids*centroids')));
% Evaluate Green's function at distance matrix
H01 = besselh(0,1,1i*s/c0*D);
% Evaluate H_1^1 at appropriate vector to make exact integration over a
% circle easier. 
H11 = besselh(1,1,1i*s/c0*sqrt(triAreas/pi));
% Compute index of refraction at centroids
qj =((c0./c(centroids(:,1),centroids(:,2))).^2-ones(N,1));
% Allocate memory for stiffness matrix 
K = zeros(N,N);

% On diagonal, stiffness matrix is computed differently from off diagonal
K(1:N+1:N^2) = -triAreas.*(4*1i*c0^2/s^2+(2*pi*1i*c0/s).*H11);

% Compute off-diagonal elements 
for i = 1:N
for j=1:N
  if i~=j
      K(i,j) = triAreas(j)*H01(i,j); % Add in further contribution from triangle areas below. 
  end
end
K(i,:) = triAreas(i)*K(i,:);
end

% Create a sparse matrix whose diagonal elements are exactly the ones
% created above. 
M = spdiags(triAreas./qj,0,N,N);
K = (1i*s^2)/(4*c0^2)*K;



end
