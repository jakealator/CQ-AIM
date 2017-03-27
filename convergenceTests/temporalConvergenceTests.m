clear

addpath(genpath('../modules')) % Contains the programs which actually compute

%--Begin definitions and parameters--
forwardParamsTime

mesh = 'arminCircleh0131';
meshStruct = initialize_mesh(mesh,1);
N=meshStruct.nt; % number of centroid points
meshExt = 'arminCircleh0958'; % Just because it's very coarse!
meshStructExt = initialize_mesh(meshExt,1);
nExt=meshStructExt.nt; % number of centroid points
xExt = generateCentroids(meshStructExt,nExt);
% Build auxillary spatial parameters, independent of time step 
[femStruct, farFieldStruct, iElements, jElements, multipoleMatrix, nearFieldDistances, P,flatP] = generateAuxillaryParams(meshStruct, N, M, d);
qc = (1./(c(femStruct.centroids).^2)-1);

mList = [40,60,80,100,120,140,160,180,200];

errK = zeros(length(mList),1);

exactBool = 1;
if exactBool
    fileID = fopen('errorTime.txt','a');
end


for i=1:length(mList)
    % A few parameters need to change with the changing time discretization
    MTime=mList(i); %Needs to be an even number right now
    dt = T/MTime;
    t = linspace(0,T,MTime+1);
    % lambda chosen according to L. Banjai, S. Sautor, Rapid solution of the
    % wave equation in unbounded domains, 2008. 
    lambda = max((dt)^(3/MTime),(eps)^(1/(2*MTime)));
    s = (1/dt)*delta(lambda*exp((-2*pi*1i*(0:MTime))/(MTime+1)));
    

    %---- Find scattered Field ---%
    us = cqify(femStruct, farFieldStruct, N, MTime, s, t, lambda, flatP, P, ...
       iElements, jElements, nearFieldDistances, uiFun, qc, c0);
%     us = RKcqify(femStruct, farFieldStruct, N, MTime, t, dt, lambda, flatP, P, ...
%         iElements, jElements, nearFieldDistances, uiFun, c, c0, qc);
    
    % Calculate series solution
    cI = c(0,0); % Only needed for exact solution. 
    circRadius = 0.275;
    [usSeriesExt,usSeries] = circPlaneSolution(circRadius, c0, cI, mesh, meshExt, femStruct.centroids, xExt, uiFun,t,aUi,bUi,cUi);

MTimePart = floor(0.85*MTime);
errK(i) = calculateRelativeError(MTimePart-1,usSeries(:,1:MTimePart-1),us(:,2:MTimePart))
% errK(i) = calculateRelativeError(MTime,usSeries,us) % for RK
N/MTime


figure
for j=1:MTime
pltsln(meshStruct,femStruct.centroids,us(:,j+1))
axis([-0.5,0.5,-0.5,0.5,min(min(usSeries)),max(max(usSeries))])
caxis([min(min(us)),max(max(us))])
hold on
pltsln(meshStruct,femStruct.centroids,usSeries(:,j))
title(t(j))
% view(2)
pause(0.1)
hold off

end

pause(1)

end

figure
loglog(1./mList,errK,'r.-')
hold on
loglog(1./(mList),25./(mList),'k-')
loglog(1./(mList),2.7./(mList).^(1/2),'b-')
%loglog(h,(h).^(-2),'b-')

xlabel('Time step size k')
ylabel('Time-Space error')
title('Convergence Rate for fixed mesh size')
legend('h=0.013','O(k)', 'O(k^{1/2})')
fclose(fileID);
