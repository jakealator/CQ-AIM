function [H,dH]=hankel(N,z)
    % Compute Hankel functions of first kind and order -N<= n <= N
    % and their derivatives
    H=besselh((-N-1):(N+1),1,z);
    dH=(H(1:2*N+1)-H(3:2*N+3))/2;
    H=H(2:2*N+2);
end
  