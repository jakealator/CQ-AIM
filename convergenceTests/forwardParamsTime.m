% forwardParamsTime.m
%
% Script containing physical parameters related to CQ-AIM simulation. 
% User specifies: c0 - speed of sound in free space (almost certainly 
%                      should be set to 1)
%                  c - an anonymous function of (x,y) defining the speed of
%                      sound in the inhomegeneity. 
%               mesh - string specifying the interior mesh file
%


%---- User Editable Parameters ----%
c0 = 1; % Speed of sound in free space (default=1)
c=@(x,y)(sqrt(2)); % Speed of sound in inhomogeneity (should be either 0<c<1 or 1<c<2 in general)
mesh = 'circle'; % name of mesh file (see ./modules/meshes for names)
incD = [1,0];
aUi = 4;
bUi = 1.4;
cUi = 2;
uiFun=@(x,y,t)(sin(aUi.*(t-(incD(1).*x+incD(2).*y)./c0)).*exp(-bUi.*(t-(incD(1).*x+incD(2).*y)./c0-cUi).^2));


%---- Values for CQ ----%
delta = @(xi)(1-xi); % Backwards Euler
%delta = @(xi)(3/2*(1-4/3*xi+1/3*xi.^2)); % BDF2
%delta = @(xi)((2*(1-xi)./(1+xi))); % Trapezoidal rule

MTime=100; %Needs to be an even number right now
T = 6;
dt = T/MTime;
t = linspace(0,T,MTime+1);
% lambda chosen according to L. Banjai, S. Sautor, Rapid solution of the
% wave equation in unbounded domains, 2008. 
lambda = max((dt)^(3/MTime),(eps)^(1/(2*MTime)));
s = (1/dt)*delta(lambda*exp((-2*pi*1i*(0:MTime))/(MTime+1)));

%----------------------------------%


%---- Derived parameters ----%

% Initialize mesh variables used throughout code
% meshStruct = initialize_mesh(mesh,1); % Initialize mesh (always use p=1 in second argument)
% N=meshStruct.nt; % number of centroid points

% Far field parameters
M = 2; % (MUST be 2) % Multipole expansion parameter 
d = 8; % Near field distance parameter (>=3)
%-----------------------------%