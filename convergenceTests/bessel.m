  function [J,dJ]= bessel(N,z)
    % Compute Bessel function J_n order -N<= n <= N
    % and their derivatives
    J=besselj((-N-1):(N+1),z);
    dJ=(J(1:2*N+1)-J(3:2*N+3))/2;
    J=J(2:2*N+2);
end