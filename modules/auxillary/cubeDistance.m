% cubeDistance.m
% Created: on 03-08-2017 by JDR at home in Newark
% Last Modified: 
%
% Inputs: 
% Outputs: 
%
% Calculates the expansion cube distance defined by dist(Ca,Cb) = min_{u\in
% Ca, v\in Cb} |u-v| where |u|=max{|u_x|,|u_y|}. 

function dist = cubeDistance(CaX,CaY,CbX,CbY)

Na = length(CaX);
Nb = length(CbX);
currentDist = zeros(Na,Nb);

% dist = 1E14; 

for i=1:Na
    for j=1:Nb
        currentDist(i,j) = max(abs(CaX(i)-CbX(j)),abs(CaY(i)-CbY(j)));
%         
%         if  currentDist < dist
%             dist = currentDist;
%         end
    end
end
dist = min(min(currentDist));




end