extBool = 1;
figureBool = 1;
movieFlag = 1;


c0 = 1;
c=@(x,y)(sqrt(2));
d=[1,0];

aUi = 4;
bUi = 1.4;
cUi = 2;
ui=@(x,y,t)(sin(aUi.*(t-(d(1).*x+d(2).*y)./c0)).*exp(-bUi.*(t-(d(1).*x+d(2).*y)./c0-cUi).^2));

%ui=@(x,y,t)(heaviside(t-1/c0*sqrt((x-1).^2+(y+1).^2))./(2*pi*sqrt(t.^2-1/c0^2*((x-1).^2+(y+1).^2))));

% Different incident wave
%bUi = 2.4;
%tlagUi = 3;
%ui=@(x,y,t)(sin(c0*(t-tlagUi)-(x.*d(1)+y.*d(2))).*exp(-bUi*(c0*(t-tlagUi)-(x.*d(1)+y.*d(2))).^2));
% to change direction of plane wave, set d=[cos(th),sin(th)]
mesh = 'arminCircleh0261';
meshExt = 'arminCircleExterior';


% Initialize mesh (use p=1 even though that's not true)
meshStruct = initialize_mesh(mesh,1);
meshStructExt = initialize_mesh(meshExt,1);

N=meshStruct.nt; % number of centroid points
NExt=meshStructExt.nt; % number of centroid points in exterior domain

qc=@(x,y)((c0/c(x,y))^2-1+x.*zeros(N,1));

% Values for CQ
%delta = @(xi)(1-xi); % Backwards Euler
delta = @(xi)(3/2*(1-4/3*xi+1/3*xi.^2)); % BDF2

minM=30;
maxM=130;
errK = zeros(length(minM:10:maxM),1);
errKExt = zeros(length(minM:10:maxM),1);
T = 3.5;


exactBool = 1;
if exactBool
    cI = c(0,0); % Only needed for exact solution. 
    circRadius = 0.275;
end

