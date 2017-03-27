% Gives Fourier series solution to the wave scattering problem against a
% circular inhomogeneity with a plane wave incident field. 


function [usSeriesExt,usSeries] = circPlaneSolution(circRadius, c0, cI, mesh, meshExt, x, xExt, uIncident,t,a,b,c)

% Parameters for meshing
meshStruct = initialize_mesh(mesh,1); % contains obstacles
meshStructExt = initialize_mesh(meshExt,1);
nTris=meshStruct.nt; % number of centroid points
nTrisExt=meshStructExt.nt;

% Solving parameters
N=24; %Number of terms in series expansion
nT = length(t);
nW=128; % number of frequencies used
omin=-32.; % minimum frquency
omax=32.; % maximum frequency
frequencies=linspace(omin,omax,nW);

an=zeros(2*N+1,nW);
bn=zeros(2*N+1,nW);
cn=zeros(nW,1);

% Initialize solution vectors
uHat = zeros(nTris,nW);
usHatExt = zeros(nTrisExt,nW);


% Find the series coefficients. To see derivation of these, expand
% interior solution in {Jn(kc0/cI*r)*exp(in*theta)} and exterior solution in
% {Hn1(k*r)*exp(in*theta)}+exp(ikd.x) and use the Jacobi-Anger expansion to
% expand exponential term (which is why d=[1,0] is required!). 
% Enforce continuity and continuity of first
% derivative on the boundary of the circle to get below system. cn is
% needed in case we're not just using plane waves (this is true in the
% time-domain case, e.g.). 

warning('off','MATLAB:nearlySingularMatrix'); %

for m=0:nW-1
    waveNumber=frequencies(m+1);
    cn(m+1)=FX(waveNumber,a,b,c);
    [J,dJ]=bessel(N,waveNumber*c0/cI*circRadius);
    [Jout,dJout]=bessel(N,waveNumber*c0*circRadius);
    [H,dH]=hankel(N,waveNumber*c0*circRadius);
    for n=-N:N
        seriesCoefficientMatrix=[[-H(n+N+1),J(n+N+1)];...
            [-waveNumber*dH(n+N+1),waveNumber*c0/cI*dJ(n+N+1)]];

        seriesCoefficientRHS=cn(m+1)*1i^n*[Jout(n+N+1);waveNumber*dJout(n+N+1)];

        coefficientVector=seriesCoefficientMatrix\seriesCoefficientRHS;

        an(n+N+1,m+1)=coefficientVector(1);
        bn(n+N+1,m+1)=coefficientVector(2);
    end
end

% Add in theta-component
% First interior part
theta=atan2(x(:,2),x(:,1));
 % To keep d constistent
 for i=1:nTris
     if theta(i)<0
         theta(i)=theta(i)+2*pi;
     end
 end
 r = sqrt(x(:,1).^2+x(:,2).^2);
 
 for m=0:nW-1
     waveNumber=frequencies(m+1);
     for i=1:nTris
         [J,~]=bessel(N,waveNumber*c0/cI*r(i));
         uHat(i,m+1) = sum(bn(:,m+1).'.*J.*exp(1i*(-N:N)*theta(i))); % .' computes the nonconjugate transpose
     end
 end
% Now exterior part
 theta=atan2(xExt(:,2),xExt(:,1));
 % To keep d constistent
 for i=1:nTrisExt
     if theta(i)<0
         theta(i)=theta(i)+2*pi;
     end
 end
 r = sqrt(xExt(:,1).^2+xExt(:,2).^2);
 for m=0:nW-1
 waveNumber=frequencies(m+1);
     for i=1:nTrisExt
         [H,~]=hankel(N,waveNumber*r(i));
         usHatExt(i,m+1) = sum(an(:,m+1).'.*H.*exp(1i*(-N:N)*theta(i)));
     end
end
 
 
 % Add time-component
uSeries=zeros(nTris,nT+1); % Time-depnendent Total Field (Interior)
usSeries=zeros(nTris,nT+1); % Time-dependent Scattered Field (Interior)
usSeriesExt=zeros(nTrisExt,nT+1); % Time-dependent Scattered Field (Exterior)
for j=1:nT
    for m=0:nW-1
        waveNumber=frequencies(m+1);
        
        uSeries(:,j)=uSeries(:,j)+uHat(:,m+1)*exp(-1i*waveNumber*t(j));
        usSeriesExt(:,j)=usSeriesExt(:,j)+usHatExt(:,m+1)*exp(-1i*waveNumber*t(j));
    end
end
for i=1:nTris
    for j=1:nT
        usSeries(i,j)=real((frequencies(2)-frequencies(1))/(2*pi)*uSeries(i,j))-uIncident(x(i,1),x(i,2),t(j));
    end
end
usSeriesExt=real((frequencies(2)-frequencies(1))/(2*pi)*usSeriesExt);

 warning('on','MATLAB:nearlySingularMatrix')
end
