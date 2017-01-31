% triangleSizeParameter.m
% Created: 14:00 02-20-2014 by jdr at UD
% Last Modified: 16:00 04-15-2014 by jdr at UD
% hK = triangleSizeParameter(sideLength)
% Input: Mesh.nt x 3 array with length of sides for each trangle in mesh
% Output: Mesh.nt x 1 with size of largest circumscribe circle around each
% triangle in mesh
%
% Description: Calculates h_k, the radius of the largest circumscribed
% circle around each triangle in a mesh. 


function hK = triangleSizeParameter(sideLength)

 % Rename side-lengths of each triangle for each of calculation
 a = sideLength(:,1);
 b = sideLength(:,2);
 c = sideLength(:,3);
 s = 1/2*(a+b+c);
 
 % Calcuate radius of largest circumscribed circle based on the side
 % lengths calculated above. Formula of d =
 % 2*a*b*c/(sqrt((a+b+c)(-a+b+c)(a-b+c)(a+b-c))) comes from wikipedia. 
 hK = 1/2*(a.*b.*c./(2*sqrt(s.*(s-a).*(s-b).*(s-c))));
 
end
