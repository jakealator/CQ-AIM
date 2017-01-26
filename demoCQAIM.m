% demoCQAIM.m 
% Created: 01-26-2016 by JDR in Newark
% Last Modified: 
%
% Runs a demo of the convolution quadrature-adaptive integral method scheme
% to simulate acoustic wave propagation through an inhomogeneous medium.
% The geometry, numerical, and physical parameters are set in the files 
% called fowardParams.m and scatteringParams.m in ./demo. 

% To run this demo, simply click on Run or (if using Windows Shortcuts)
% type Shift+Enter

% For details about how this works and what it's doing, see README.md and 
% Chapter 5 of my Ph.D thesis.

%------ Begin Demo ------%

% Add required folders 
addpath genpath('./demo') % Contains physical, geometric, and computational parameters for demo
addpath genpath('./modules') % Contains the programs which actually compute

%---- Initialize parameters ----%



%---- Generate scattered field data ----%



%---- Plot results ----%





%------ End Demo ------%