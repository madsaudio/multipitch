function [L,w]=order_orth(x,F,M);
%ORDER_ORTH   Order estimator based on subspace orthogonality
%
% Syntax:
%   L=order_orth(x,F,M);
%
% About:
%   This file is part of the Multi-Pitch Estimation Toolbox for the book
%   M. G. Christensen and A. Jakobsson, Multi-Pitch Estimation, Morgan &
%   Claypool Publishers, 2009.
%
% Input:
%   x        input signal (assumed complex)
%   F        FFT size (default: 4096)
%   M        covariance matrix size (default: length(x)/2)
%   
% Output:
%   L        estimated signal subspace rank/number of sinusoids
%
% Description:
%   Function that determines the dimension of the signal subspace using the
%   subspace orthogonality principle, more specifically the angles between 
%   the subspaces, as described in Section 4.7 of Christensen and Jakobsson
%   (2009). This is done using an unconstrained model of the frequencies, 
%   meaning that the sinusoids are not restricted to being integer 
%   multiples of a fundamental. Can be used to determine the dimension of 
%   the noise subspace for multi-pitch signals or the number of harmonics 
%   for single-pitch signals. Tests only for L>0.
%
% Example:
%   L=order_orth(x,2048,length(x)/2);
%
% Implemented By:
%   Mads G. Christensen (mgc@es.aau.dk)
%
if nargin<3,F=4096;end
if nargin<2,M=floor(length(x)/2);end
N=length(x);
x=x(:);
R=covm(x,M);
[U,D,V]=svd(R);
T=abs((fft((U),F))).^2;
L_set=[1:M-2];
for k=L_set,
    q=min([k M-k]);
    P=sum(T(:,k+1:M),2);
    [f,p]=findpeaks(P,k,0);
    J(k)=sum(p)/(M*q);
end
[tmp,ndx]=min(J);
L=L_set(ndx);