% areaTri.m
% Created: Spring 2014 by jdr
% Last Modified: 01-31-2017 by jdr at home in Newark
%
% Input: (x,y) - two 3x1 vectors containing the x-y coordinates of a triangle's corners
% Output: area - the area of the given triangle. 
%
% Generates area of a triangle based on its nodal coordinates. 

function area = areaTri(x,y)

% From the 'shoelace' formula
area = 1/2*abs((x(1)-x(3))*(y(2)-y(1))-(x(1)-x(2))*(y(3)-y(1)));

end