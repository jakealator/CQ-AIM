% generateNearFieldDistances.m
% Created: 02-01-2017 by JDR in Newark
% Last Modified: 
%
% Input: centroids       - Nx1 matrix containing centroids of each element
%        multipoleMatrix - NxN sparse matrix such that (i,j)=1 if element i
%                          and j are within 1 element of each other (are
%                          touching)
%        N               - Number of elements
%        M               - Multipole expansion (needs to be 1 or 0 right now)
% Output: nearFieldDistances: NxN sparse matrix containing the distances
%         between elements within M elements of each other. 

function nearFieldDistances = generateNearFieldDistances(centroids, iElements,jElements, N,M)

if M==0 % when M=0, treat everything as 'near field.' 
    nearFieldDistances = sqrt(bsxfun(@plus,full(dot(centroids',centroids',1)),full(dot(centroids',centroids',1))')-full(2*(centroids*centroids')));
elseif M==1
    % Find triangles that are in the near field of each other
%     [I,J] = ind2sub([N,N],find(multipoleMatrix));
    sElement = zeros(length(iElements),1); % Initialize array to hold distances

    % Calculate distance using d(x,y) = sqrt((x(1)-y(1))^2+(x(2)-y(2))^2)
     for k=1:length(iElements)
         sElement(k) = sqrt((centroids(iElements(k),1)-centroids(jElements(k),1))^2+...
             (centroids(iElements(k),2)-centroids(jElements(k),2))^2);
     end

     nearFieldDistances = sparse(iElements,jElements,sElement,N,N);
else
    error('Sorry, generateNearFieldDistances only works when M=1.')
end

end