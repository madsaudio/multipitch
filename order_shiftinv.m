function [L]=order_shiftinv(x,M);
%ORDER_SHIFTINV   Order estimator based on subspace shift-invarinace
%
% Syntax:
%   L=order_shitinv(x,M);
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
%   shift-invariance property as described in Section 4.9 of Christensen 
%   and Jakobsson (2009), also known as the ESTER method. For single-pitch 
%   signals this is the same as the number of harmonics, but for 
%   multi-pitch signals it is equal to the total number of sinusoids. 
%   Cannot test for a 0th order model.
%
% Example:
%   L=order_shiftinv(x,length(x)/2);
%
% Implemented By:
%   Mads G. Christensen (mgc@es.aau.dk)
%
if nargin<2,M=floor(length(x)/2);end
N=length(x);
x=x(:);
R=covm(x,M);
[U,D,V]=svd(R);
K_set=[1:M-2];
for l=1:length(K_set),
  k=K_set(l);
  A=U(1:M-1,1:k);
  B=U(2:M,1:k);
  P=(A\B); 
  J(l)=norm(A*P-B,2)^2; 
end
[tmp,ndx]=min(J);
L=K_set(ndx);