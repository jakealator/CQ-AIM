clear

addpath(genpath('../modules')) % Contains the programs which actually compute

%--Begin definitions and parameters--
forwardParamsTime

meshList = {'arminCircleh0958';'arminCircleh0512';'arminCircleh0261';'arminCircleh0131';'arminCircleh0066'};
meshExt = 'arminCircleh0958'; % Just because it's very coarse!
meshStructExt = initialize_mesh(meshExt,1);
nExt=meshStructExt.nt; % number of centroid points
xExt = generateCentroids(meshStructExt,nExt);

errH = zeros(length(meshList),1);

exactBool = 1;
if exactBool
    fileID = fopen('error.txt','a');
end


%matlabpool open 
for i=1:length(meshList)
    % Initialize mesh (use p=1 even though that's not true)
    mesh = meshList{i};
    meshStruct = initialize_mesh(mesh,1);
    N=meshStruct.nt; % number of centroid points
    
    [femStruct, farFieldStruct, iElements, jElements, multipoleMatrix, nearFieldDistances, P,flatP] = generateAuxillaryParams(meshStruct, N, M, d);
    qc = (1./(c(femStruct.centroids).^2)-1);

    %---- Find scattered Field ---%
    us = cqify(femStruct, farFieldStruct, N, MTime, s, t, lambda, flatP, P, ...
       iElements, jElements, nearFieldDistances, uiFun, qc, c0);
%     us = RKcqify(femStruct, farFieldStruct, N, MTime, t, dt, lambda, flatP, P, ...
%         iElements, jElements, nearFieldDistances, uiFun, c, c0, qc);
    
    % Calculate series solution
    cI = c(0,0); % Only needed for exact solution. 
    circRadius = 0.275;
    [usSeriesExt,usSeries] = circPlaneSolution(circRadius, c0, cI, mesh, meshExt, femStruct.centroids, xExt, uiFun,t,aUi,bUi,cUi);


errH(i) = calculateRelativeError(MTime-1,usSeries(:,1:MTime-1),us(:,2:MTime))
% errH(i) = calculateRelativeError(MTime,usSeries,us) % for RK
N/MTime

figure
for j=1:MTime
pltsln(meshStruct,femStruct.centroids,us(:,j+1)./(eps+usSeries(:,j)))
% axis([-0.5,0.5,-0.5,0.5,min(min(us)),max(max(us))])
%  caxis([min(min(us)),max(max(us))])
title(t(j))
% view(2)
pause(0.1)

end

end

h = [0.0958,0.0512,0.0261,0.0131,0.0066];
figure
loglog(h,errH(1:4),'r.-')
hold on
loglog(h,h,'k-')
%loglog(h,(h).^(-2),'b-')

xlabel('Space step size h')
ylabel('Time-Space error')
title('Convergence Rate for fixed number of time-steps')
legend('M=100','O(h)')
fclose(fileID);
