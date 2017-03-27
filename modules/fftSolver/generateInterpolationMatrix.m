% generateInterpolationMatrix.m
% Created: 02-16-2017 by JDR in Newark
% Last Modified:
%
% Inputs: centroids           -
%         centers             -
%         rectangularElements -
%         meshStruct          -
% 
% Output: V - Nx9 matrix containing interpolation coefficients between
% Galerkin and delta functions
%
% Finds coefficient matrix to interpolate between p0 Galerkin basis 
% functions on FEM triangulation and "delta basis functions" on Cartesian
% discretization. 

function V = generateInterpolationMatrix(centers, rectangularElementsX, rectangularElementsY, midpointsX,midpointsY, triAreas, N,varargin)

V = zeros(N,9);

if nargin == 7
    % Order of m-vector containing exponents
    mVec = [[0,0];[0,1];[0,2];[1,0];[1,1];[1,2];[2,0];[2,1];[2,2]];

    % Use quadrature against bFun to calculate rhs
    bFun = @(x,y,m,midpointsX,midpointsY)((midpointsX-x)^m(1)*(midpointsY-y)^m(2));


    % Compute V for each element/expansion cube
    for j=1:N
        W = zeros(9,9);
        b = zeros(9,1);

        % Use midpoint rule to compute righ-hand-side. This is exact for
        % polynomials up to 2nd order so is fine for M=2
        % WAIT: Does this incur error? I think it does. 
        for k=1:9
            b(k) = (bFun(centers(j,1),centers(j,2),mVec(k,:),midpointsX(j,1),midpointsY(j,1))+...
                bFun(centers(j,1),centers(j,2),mVec(k,:),midpointsX(j,2),midpointsY(j,2))+...
                bFun(centers(j,1),centers(j,2),mVec(k,:),midpointsX(j,3),midpointsY(j,3)));

            % Now compute (shifted) Vandermonde matrix W
            for ell = 1:9
                W(k,ell) = (rectangularElementsX(j,ell)-centers(j,1))^mVec(k,1)...
                    *(rectangularElementsY(j,ell)-centers(j,2))^mVec(k,2);
            end
        end
        % Need |K|/3 for midpoint quadrature rule
        b = triAreas(j)/3*b;

        Vsquare = W\b;

        V(j,:) = Vsquare(:);


    end
else
    
    waveNumber = varargin{length(varargin)-1};
    centroids = varargin{length(varargin)};
	farFieldFun = @(x,y)(exp(1i*pi/4)/sqrt(8*pi*1i*waveNumber)*exp(waveNumber*(x(1)*y(1)+x(2)*y(2))));
    farFieldGridTheta = linspace(0,2*pi,9); farFieldGridTheta=farFieldGridTheta(:);
    farFieldGrid = [cos(farFieldGridTheta), sin(farFieldGridTheta)];
	for j=1:N
		Q = 9;
		W = zeros(Q,9);
		b = zeros(Q,1);

		for k=1:Q
			b(k) = farFieldFun(centroids(j,:),farFieldGrid(k,:));
			for ell=1:9
				W(k,ell) = farFieldFun([rectangularElementsX(j,ell),rectangularElementsY(j,ell)],farFieldGrid(k,:));
			end
		end
		Vsquare = lsqr(W,b,1E-3,20);
		V(j,:) = Vsquare(:);
        
    end
end
		
		

end
