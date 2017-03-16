% applyFundamentalSolution.m
% Created: 03-10-2017 by JDR in Newark
% Last Modified: 03-16-2017 by JDR in DC
% 
% Inputs: X - nGx2 vector of Cartesian points
% Outputs: GD - nGx1 matrix of fundamental solution applied to the first
% row of pairwise distance matrix. 
%
% Just applies fundamental solution


function [flatGD, GD] = applyFundamentalSolution(waveNumber, X)

fundamentalSolution=@(x)(reshape(1i/4*besselh(0,1,1i*waveNumber*((x~=0).*x+(abs(x)<1E-14).*1E300)),size(x)));

% Calculate pairwise distances 
% D = sqrt(bsxfun(@plus,dot(X',X',1),dot(X',X',1)')...
%     -(2*(X*X')));
% apply fundamental solution
D = sqrt((X(1,1)-X(:,1)).^2+(X(1,2)-X(:,2)).^2);
GD = fundamentalSolution(D);

% Build 'sparse' toeplitz representation for use in fft
N = sqrt(length(GD));
c = zeros(N,N);
for j=1:N
c(:,j) = GD((j-1)*N+1:j*N);
end

flatGD = embedToeplitzBlockInCirculantBlock(c);

end
