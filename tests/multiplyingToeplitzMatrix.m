N=4;

% Multiplication with circulant matrix:
xVec = rand(N,1);
v = rand(1,N);
A = toeplitz([v(1) fliplr(v(2:end))], v);

Ax = ifft(fft(A(:,1)).*fft(xVec));
Ax-A*xVec

% Multiplication with toeplitz via embedding: 
tRow = rand(N,1);
tCol = rand(N,1); tCol(1) = tRow(1);
A = toeplitz(tCol,tRow);

% turn into circulant matrix
B = toeplitz([0; flipud(tRow(2:end))],[0; flipud(tCol(2:end))]);

Ax = ifft(diag(fft([A(:,1); B(:,1)]))*fft([xVec; zeros(N,1)]));
Ax(1:N)-A*xVec

A2x = ifft(fft(embedToeplitzInCirculant(tCol,tRow)).*fft([xVec; zeros(N,1)]));
A2x(1:N)-A*xVec


% CBBC matrix multiplication
N=3; M=3;
xVec = rand(N,M);

if (N==2 && M==2)
    C0 = [1 2; 2 1]; C1 = [3 4; 4 3];
    A = [C0 C1; C1 C0];
end
if (N==3 && M==3)
    C0=[ 1  2  3;  3  1  2;  2  3  1]+1i*rand(1,1);
    C1=[ 4  5  6;  6  4  5;  5  6  4]+1i*rand(1,1);
    C2=[ 7  8  9;  9  7  8;  8  9  7]+1i*rand(1,1);
    A =[C0 C2 C1; C1 C0 C2; C2 C1 C0];
end

Ax = ifft2(fft2(reshape(A(:, 1), N, M)).*fft2(xVec))
Ax1=A*xVec(:);
Ax(:)-Ax1(:)

% Block toeplitz via embedding: 
N=3;
x = linspace(0,1,N)'; y=x;
[X,Y]=meshgrid(x,y);
X = [X(:), Y(:)];
D = sqrt(bsxfun(@plus,dot(X',X',1),dot(X',X',1)')-(2*(X*X')));
xVec = rand(N,N);

tRow = D(1,:); tRow=tRow(:);
tCol = D(:,1); tCol=tCol(:);
B = toeplitz([0; flipud(tRow(2:end))],[0; flipud(tCol(2:end))]);



A = [D(:,1); B(1,:).'];

Ax = ifft2(fft2(reshape(A, N, 2*N)).*fft2([xVec zeros(N,N)]))


Ax = ifft(fft2(A).*fft2([xVec(:); zeros(N^2,1)]));
Ax1=D*xVec(:);
Ax(1:N^2)-Ax1(:)

Ax2=toeplitz([A(1) fliplr(A(2:end).')], A.')*[xVec(:);zeros(N^2,1)];
Ax2(1:N^2)-Ax1(:)
