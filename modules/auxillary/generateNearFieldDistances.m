% generateNearFieldDistances.m
% Created: 02-01-2017 by JDR in Newark
% Last Modified: 
%
% 

function nearFieldDistances = generateNearFieldDistances(centroids, multipoleMatrices, N,M)

if M>1
    error('Sorry, generateNearFieldDistances only works when M=1.')
end

% Find locations where 
[I,J] = ind2sub([N,N],find(multipoleMatrices));
sElement = zeros(length(I),1);
 
 for k=1:length(I)
     sElement(k) = sqrt((centroids(I(k),1)-centroids(J(k),1))^2+...
         (centroids(I(k),2)-centroids(J(k),2))^2);
 end
 
 nearFieldDistances = sparse(I,J,sElement,N,N);
 
end