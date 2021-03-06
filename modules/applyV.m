% applyV.m
% Created by: JDR on 02-20-2017 in Newark
% Last Modified:
% 
% Inputs: x    - The Nx1 vector we apply V to
%         K    - NxN sparse matrix of near field componebnts 
%         Ng   - Number of points in the Cartesian grid
%         N    - Number of points in the FEM grid
%         fftG - NgxNg matrix of the fft of the fundamental solution
%                evaluated at Cartesian grid points
%         P    - The NgxN interpolation matrix
%
% Outputs: Vx - An Nx1 vector of the result of Vx
%
% Uses the AIM scheme to apply the operator V, which is the discrete
% version of the Lippmann-Schwinger integral operator. Note that this does
% NOT apply the identity matrix component of the operator! 

function Vx = applyV(x,K,fftG, P, waveNumber)

% % %% Begin by computing Vfar*x

% % Interpolate x onto Cartesian grid as xhat
xHat = P*x;
nG = sqrt(length(xHat));

% % Do the convolution with ffts. 
% To do this, first zero pad the vector of interest
xHat = reshape(xHat, nG, nG);
xHatLong = [xHat, zeros(nG,nG); zeros(nG,2*nG)];

% Use convolution theorem (remember that fft is now in circulant form)
yHatLong = ifft2(fftG.*fft2(xHatLong));
yHat = yHatLong(1:nG, 1:nG); % Only take the first half of elements

% % Interpolate convolution operator back onto FEM grid
% % y is now Vfarx.
y = P.'*yHat(:);


% Now add in the near field components. 
Vx = K*x+waveNumber^2*y;


end