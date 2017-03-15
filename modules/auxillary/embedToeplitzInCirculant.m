% embedToeplitzInCirculant.m
% Created on: 03-15-2017 by JDR in Washington DC
% Last Modified:
%
% Input: tCol - Nx1 or 1xN vector, first column of toeplitz matrix
%        tRow - Nx1 or 1xN vector, first row of toeplitz matrix
% Output: c - Mx1 first row of circulant matrix
% 
% Takes an NxN toeplitz matrix T and embeds it in an MxM circulant matrix C
% such that if x is an Nx1 matrix then Tx=ifft(fft(C(:,1)).*fft(x)).
% Outputs first row of C. 


function c = embedToeplitzInCirculant(tCol,tRow)

tRow=tRow(:); tCol=tCol(:);
c = [tCol; 0; flipud(tRow(2:end))];
    

end