function [M,K,rhs,usHat] = generateUSHat(ui,triangleAreas,H01,H11,qj,c0,s,N, rhsSmallBool)


% Initialize arrays for later
K = zeros(N,N);
rhs = zeros(N,1);

%--Begin assembly of matrices--

% Off the diagonal, use straightforward p0 elements: int_K1 int_K2
% i/4s^2H_0^(1)(k||x_i-x_j||). On the diagonal, compute the same on the
% subtriangles with barycentric vertex and two of the original vertices as
% vertices. 
for i = 1:N
    rhsSum=0;
  for j=1:N
      if i~=j
          rhsSum = rhsSum + qj(j)*ui(j)*triangleAreas(j)*H01(i,j);
      else          
          rhsSum = rhsSum + qj(j)*ui(j)*triangleAreas(j)*...
               (-(4*1i*c0^2/s^2+(2*pi*1i*c0/s)*H11(j)));
      end
  end
rhs(i) = triangleAreas(i)*rhsSum;  
end
rhs = -(1i*s^2)/(4*c0^2)*rhs;

rhsNorm = norm(rhs);

%compute less if possible!
tol=1E-9;
if rhsSmallBool && rhsNorm<=tol
    usHat = zeros(size(rhs));
    M=0;
else   
    K(1:N+1:N^2) = -triangleAreas.*(4*1i*c0^2/s^2+(2*pi*1i*c0/s).*H11);%...Ss

    for i = 1:N
      for j=1:N
          if i~=j
              K(i,j) = triangleAreas(j)*H01(i,j);
          end
      end
      K(i,:) = triangleAreas(i)*K(i,:);
    end
    % Create a sparse matrix whose diagonal elements are exactly the ones
    % created above. 
    M = spdiags(triangleAreas./qj,0,N,N);
    
    K = (1i*s^2)/(4*c0^2)*K;
    usHat = (M+K)\(rhs);
    usHat = usHat./qj;
end




end