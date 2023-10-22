function a=amp_wls(x,w,Q);
%AMP_WLS   Complex amplitude estimator based on WLS
%
% Syntax:
%   a=amp_wls(x,w,Q);
%
% About:
%   This file is part of the Multi-Pitch Estimation Toolbox for the book
%   M. G. Christensen and A. Jakobsson, Multi-Pitch Estimation, Morgan &
%   Claypool Publishers, 2009.
%
% Input:
%   x        input signal (assumed complex)
%   w        vector of L frequencies
%   Q        M-by-M noise covariance matrix (default: identity)
%   
% Output:
%   a        vector of L complex amplitudes for the frequencies in w
%
% Description:
%   The function estimates the complex amplitudes associated with a set of
%   frequencies. It does so using the weighted least-squares method as 
%   described in Section 5.2 of Christensen and Jakobsson (2009). It takes 
%   a noise covariance matrix estimate as input or uses the identity matrix
%   if none is given.
%
% Example:
%   a=amp_wls(x,w0*[1:L],Q);
%
% Implemented By:
%   Mads G. Christensen (mgc@es.aau.dk)
%
N=length(x);
x=x(:);
if nargin<3,Q=eye(N);end
M=size(Q,1);
w=w(:)';
L=length(w);
Z=vandermonde(w,M);
C=inv(Q);
D=sparse(diag(exp(j*w)));
A=zeros(L,L);
B=zeros(L,1);
Y=Z;
for n=1:N-M+1,
    y=x(n:n+M-1);
    A=A+Y'*C*Y;
    B=B+Y'*C*y;
    Y=Y*D;
end
a=inv(A)*B;