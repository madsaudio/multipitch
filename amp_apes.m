function a=amp_apes(x,w,M);
%AMP_APES   Complex amplitude esitmator based on APES
%
% Syntax:
%   a=amp_apes(x,w,M);
%
% About:
%   This file is part of the Multi-Pitch Estimation Toolbox for the book
%   M. G. Christensen and A. Jakobsson, Multi-Pitch Estimation, Morgan &
%   Claypool Publishers, 2009.
%
% Input:
%   x        input signal (assumed complex)
%   w        vector of L frequencies
%   M        filter length (default: length(x)/4)
%
% Output:
%   a        vector of L complex amplitudes for the frequencies in w
%
% Description:
%   The function estimates the complex amplitudes associated with a set of
%   frequencies. It does so using the extended APES estimator described in 
%   Section 5.3 of Christensen and Jakobsson (2009).
%
% Example:
%   a=amp_apes(x,w0*[1:L],floor(length(x)/4));
%
% Implemented By:
%   Mads G. Christensen (mgc@es.aau.dk)
%
N=length(x);
x=x(:);
if nargin<3,M=floor(N/4);end
w=w(:)';
R=covm(x,M);
G=N-M+1;
H=hankel(x(1:M),x(M:N));
Z=vandermonde(-w,G)/G;
Q=R-H*Z*Z'*H';
a=amp_wls(x,w,Q);