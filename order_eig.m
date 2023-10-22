function [K]=order_eig(x,M);
%ORDER_EIG   Order estimator based on the eigenvalues and MDL
%
% Syntax:
%   L=order_eig(x,M);
%
% About:
%   This file is part of the Multi-Pitch Estimation Toolbox for the book
%   M. G. Christensen and A. Jakobsson, Multi-Pitch Estimation, Morgan &
%   Claypool Publishers, 2009.
%
% Input:
%   x        input signal (assumed complex)
%   M        covariance matrix size
%   
% Output:
%   L        estimated signal subspace rank/number of sinusoids
%
% Description:
%   Function that determines the dimension of the signal subspace using the
%   ratio between the arithmetic and geometric means of the eigenvalues
%   of the covariance matrix in combination with the MDL criteroin 
%   as described in Section 4.5 of Christensen and Jakobsson (2009), which
%   is based on an assumption of white Gaussian noise. For single-pitch 
%   signals this is the same as the number of harmonics, but for 
%   multi-pitch signals it is equal to the total number of sinusoids.
%
% Example:
%   L=order_eig(x,length(x)/2);
%
% Implemented By:
%   Mads G. Christensen (mgc@es.aau.dk)
%
if nargin<2,M=floor(length(x)/2);end
N=length(x)-M+1;
x=x(:);
R=covm(x,M);
d=svd(R);
K_set=[1:M-2];
for k=K_set,
    L(k)=N*sum(log(d(k+1:M)))-N*(M-k)*log(mean(d(k+1:M)));
    J(k)=-L(k)+k*(2*M-k);
end
[tmp,ndx]=min(J);
K=K_set(ndx);
