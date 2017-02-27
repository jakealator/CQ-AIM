% generateNearFieldElements.m
% Created: 02-01-2017 by JDR in Newark
% Last Modified: 
%
% Input: N          - number of elements
%        M          - multipole expansion order (currently only 1)
%        hGrid      - grid size of rectangular grid
%        d          - user parameter defining what "far" means
% Output:  i/jElements     - vectors containing the nonzero elements of
%                            multipoleMatrix
%          multipoleMatrix - a sparse matrix #elementsx#elements with values 
%                           in {0,1,...,M} such that if element
%                           i is m elements away from j then (i,j)=m. If
%                           m>M then (i,j)=0. 
%
% Determines how far away each mesh element is from each other element, up
% to the user-determined multipole expansion. 

function [iElements, jElements, multipoleMatrix] = generateNearFieldElements(N,M,centers,h,d)

% If M=0, there is no near/far field distinction (mostly for testing)
if M==0
    [iElements, jElements] = meshgrid(1:N,1:N); % index a loop through a matrix
    iElements = iElements(:); jElements = jElements(:); % flatten to put in correct order/size
    multipoleMatrix = ones(N,N);
    
elseif M==2
    % For each center, determine which centers (rectangular mesh) are 
    % within d*h units. If the jth row is close to the ith value, matrix 
    % element (i,j)=1. Otherwise (i,j)=0.
    % Note that this is probably overly slow (O(N^2)) but is memory
    % efficient, which is more important as N grows. If desired, a quadtree
    % will reduce speed as well. 
    
    k = 0; % tracks number of non-zero elements
    for i=1:N
        for j=1:N % test all the other elements against cElement
            if (sqrt((centers(i,1)-centers(j,1))^2+(centers(i,2)-centers(j,2))^2)<d*h)
                k=k+1;
                iElements(k) = i; % iElements keeps track of the column number
                jElements(k) = j; % jElements keeps track of the row number
            end
        end
    end
    % Since we just want a one where the elements are connected, we can use a
    % ones matrix. This will be more complicated when M>1.
    multipoleMatrix = sparse(iElements,jElements,ones(k,1),N,N);
else
    error('generateNearFieldElements.m: Sorry, multipole M>1 not supported.')
end

end