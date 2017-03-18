% forwardParams.m
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
mesh = 'rectangle'; % name of mesh file (see ./modules/meshes for names)
waveNumber = 1i; % frequency (only for testing time-harmonic)
incD = [0,1];
uiHatFun = @(x)exp(waveNumber*(incD(1)*x(:,1)+incD(2)*x(:,2)));
%----------------------------------%


%---- Derived parameters ----%

% Initialize mesh variables used throughout code
meshStruct = initialize_mesh(mesh,1); % Initialize mesh (always use p=1 in second argument)
N=meshStruct.nt; % number of centroid points

% Far field parameters
M = 2; % Multipole expansion parameter
d = 3; % Near field distance parameter
%-----------------------------%