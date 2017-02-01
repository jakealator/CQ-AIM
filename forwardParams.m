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
c0 = 1; % Speed of sound in free space
c=@(x,y)(sqrt(2)); % Speed of sound in inhomogeneity 
mesh = 'twoCircles'; % name of mesh file
%----------------------------------%


%---- Derived parameters ----%

% Initialize mesh (use p=1 even though that's not true)
meshStruct = initialize_mesh(mesh,1);
N=meshStruct.nt; % number of centroid points

qc=@(x,y)((c0/c(x,y))^2-1+x.*zeros(N,1)); % Index of refraction
%-----------------------------%