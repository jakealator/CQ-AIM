clear

addpath(genpath('../modules')) % Contains the programs which actually compute


%--Begin definitions and parameters--
forwardParamsTime

meshList = {'arminCircleh0958';'arminCircleh0512';'arminCircleh0261';'arminCircleh0131';'arminCircleh0066'};%;'arminCircleh0029'};
meshExt = 'arminCircleh0958'; % Just because it's very coarse!
meshStructExt = initialize_mesh(meshExt,1);
nExt=meshStructExt.nt; % number of centroid points
xExt = generateCentroids(meshStructExt,nExt);


errH = zeros(length(meshList),1);
errHR = zeros(size(errH));


%matlabpool open 
for i=1:length(meshList)
    
    % Initialize mesh (use p=1 even though that's not true)
    mesh = meshList{i};
    meshStruct = initialize_mesh(mesh,1);
    N=meshStruct.nt % number of centroid points
    
    [femStruct, farFieldStruct, iElements, jElements, multipoleMatrix, nearFieldDistances, P,flatP] = generateAuxillaryParams(meshStruct, N, M, d);
    qc = (1./(c(femStruct.centroids).^2)-1);

    %---- Find scattered Field ---%
    us = cqify(femStruct, farFieldStruct, N, MTime, s, t, lambda, flatP, P, ...
       iElements, jElements, nearFieldDistances, uiFun, qc, c0);

    % Calculate series solution
    uSeries=zeros(N,M+1);
    cI = c(0,0); % Only needed for exact solution. 
    circRadius = 0.275;
    [usSeriesExt,usSeries] = circPlaneSolution(circRadius, c0, cI, mesh, meshExt, femStruct.centroids, xExt, uiFun,t,aUi,bUi,cUi);


MTimePart = floor(0.85*MTime);
% errH(i) = calculateRelativeError(MTimePart-1,usSeries(:,1:MTimePart-1),us(:,2:MTimePart))
errH(i) = calculateAbsoluteError(MTime,dt,femStruct.triAreas,qc,usSeries,us);
errHR(i) = calculateRelativeError(MTime,usSeries,us) ;
[errH, errHR]
end

h = [0.0958,0.0512,0.0261,0.0131,0.0066];%,0.00288];
figure
loglog(h,errH2,'r.-')
hold on
loglog(h,h,'k-')
%loglog(h,(h).^(-2),'b-')

xlabel('Space step size h')
ylabel('Time-Space error')
title('Convergence Rate for fixed number of time-steps')
legend('M=100','O(h)')

