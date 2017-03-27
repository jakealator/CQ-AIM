% sampleAndTransform.m
% Created: 03-23-2017 by JDR in Newark
% Last Modified:
%
% Input:  u      - vector of size (?)x(M+1) to be transformed
%         lambda - Related to accuracy of approximate CQ algorithm. 
%         M      - Number of time steps
% Output: v - Transformed vector of size (?)x(M+1). 
%
% Given a function sampled at appropriate points, code performs a
% normalized FFT on its second dimension (e.g., time) for use in CQ
% algorithms. 

function v = sampleAndTransform(u,lambda,M)

    uSample = bsxfun(@times,u,lambda.^(0:M)); % Normalize according to Hassel-Sayas chapter 4
    v = fft(uSample,[],2); % Take fft of time slice

end

% 
% for j=1:MTime+1
%     uiLambdaVec(:,j) = lambda^(j-1)*uiFun(femStruct.centroids(:,1),femStruct.centroids(:,2),t(j));
% end
% 
% uiHat = fft(uiLambdaVec,[],2);