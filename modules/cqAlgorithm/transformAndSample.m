% transformAndSample.m
% Created: 03-23-2017 by JDR in Newark
% Last Modified:
%
% Input:  u      - vector of size (?)x(M+1) to be inverse transformed
%         lambda - Related to accuracy of approximate CQ algorithm. 
%         M      - Number of time steps
% Output: v - Inverse transformed vector of size (?)x(M+1). 
%
% Given a function in the transform domain, code performs a
% normalized IFFT on its second dimension (e.g., time) for use in CQ
% algorithms. 
 
function uSample = transformAndSample(v,lambda,M)

u = ifft(v,[],2); % Transform back to time domain
uSample = real(bsxfun(@times,u,lambda.^(0:-1:-M))); % Normalize according to Hassel-Sayas book chapter 4


end


% Calculate us in the time domain    
% us = ifft(uScatteredHat,[],2);
% for m=1:MTime+1
%     us(:,m) = (lambda^(-m+1))*us(:,m); 
% end
% us = real(us); 