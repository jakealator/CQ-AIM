% embedToeplitzBlockInCirculantBlock.m
% Created on: 03-15-2017 by JDR in Washington DC
% Last Modified:
%
% Input: tCol - Nx1 or 1xN vector, first column of toeplitz matrix
%        tRow - Nx1 or 1xN vector, first row of toeplitz matrix
% Output: c - Mx1 first row of circulant matrix
% 
% Takes an NnxNn toeplitz block-block toeplitz matrix T with nxn blocks and
% embeds it in an MmxMm circulant block-block circulant matrix C
% such that if x is an Nnx1 matrix then Tx=ifft2(fft2(C(:,1)).*fft2(x)).
% Outputs first row of C. 


function c = embedToeplitzBlockInCirculantBlock(tCol,tRow)
% c will be of the form c = [tCol; 0; flipud(tRow(2:end))];, where tCol and
% tRow are the first row of each block

N = length(tCol(1,:)); % get number of blocks
M = length(tCol(:,1)); % get size of each block (MxM)

c = zeros(2*N*M,1);
c(1:N*M) = tCol(:);
c(N*M+1:(N+1)*M) = zeros(M,1);

for j=1:N
    c((N+j)*M+1:(N+1+j)*M = tRow(j,:);
end
    

end