% generateFullFarFieldMatrix.m
% Created: by JDR on 03-01-2017 at UD
% Last Modified:
%
% Inputs:
% 
% Outputs:
%
% Code to -slowly- compute far field interactions on auxillary Cartesian
% grid. Later this will be used as a test function, but currently just to
% to compute far field to make sure rest of AIM is working. 

function AFarx =  generateFullFarFieldMatrix(x,farFieldStruct,P, waveNumber)

% The fundamental solution for 2D. Note that this will output as something
% the size of its input and if the input is zero, the output will be zero
% up to machine precision. This is to avoid any issues with near field NaNs
% not canceling out below. 
fundamentalSolution=@(x)(reshape(1i/4*besselh(0,1,1i*waveNumber*((x~=0).*x+(abs(x)<1E-14).*1E300)),size(x)));
nG = length(farFieldStruct.farFieldGrid);

% Evaluate fundamental solution at grid distance matrix
N = length(x);
AFarx = zeros(N,1);
locs = farFieldStruct.farFieldGrid;
    for j2=1:nG % u
        j2
        for j3=1:nG % v
            AFarx = AFarx + sum(full(P(j2,:).*fundamentalSolution(norm(locs(j2,:)-locs(j3,:)))*sum(P(j3,:)'.*x)));
        end
    end

end