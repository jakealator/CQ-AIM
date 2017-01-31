function uiHat = generateUIHat(ui,lambda,m,x,t,M,N)

uiHat = zeros(N,1);

for j=1:M+1
    % Be careful with j vs (j-1) in sum: t=[t(1),...,t(M+1)] but sum goes 
    % from j=0 to j=M. 
    uiHat = uiHat + lambda^(j-1)*ui(x(:,1),x(:,2),t(j))*exp(-(2*pi*1i*(j-1)*m)/(M+1));
end

end