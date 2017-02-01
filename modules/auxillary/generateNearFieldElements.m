% generateNearFieldElements.m
% Created: 02-01-2017 by JDR in Newark
% Last Modified: 
%
% Input: meshStruct - a structure containing information about mesh given
%                     by initialize_mesh.m
%        N          - number of elements
%        M          - multipole expansion order (currently only 1)
% Output: multipoleMatrix - a matrix #elementsx#elements with values 
%                           in {0,1,...,M} such that if element
%                           i is m elements away from j then (i,j)=m. If
%                           m>M then (i,j)=0. 
%
% Determines how far away each mesh element is from each other element, up
% to the user-determined multipole expansion. 

function multipoleMatrix = generateNearFieldElements(meshStruct,N,M)

% Initialize array
multipoleMatrix = zeros(meshStruct.nt,meshStruct.nt);

% For each value in the first column of meshStruct.tris, i, determine which
% other rows of meshStruct.tris, j, also contain that value. If the jth
% row contains the ith value, matrix element (i,j)=1. Otherwise (i,j)=0.
for i=1:N
    cElement = meshStruct.tris(i,1); % pick current element to test against
    for j=1:N % test all the other elements against cElement (note everything is an integer so == is fine)
        if ((meshStruct.tris(j,1)==cElement) || (meshStruct.tris(j,2)==cElement) || (meshStruct.tris(j,3)==cElement))
            multipoleMatrix(i,j) = 1;
        end
    end
end


end