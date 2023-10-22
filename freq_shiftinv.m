function w=freq_shiftinv(x,L,M);
%FREQ_SHIFTINV   Unconstrained frequency estimator based on ESPRIT
%
% Syntax:
%   w=freq_shiftinv(x,L,M);
%
% About:
%   This file is part of the Multi-Pitch Estimation Toolbox for the book
%   M. G. Christensen and A. Jakobsson, Multi-Pitch Estimation, Morgan &
%   Claypool Publishers, 2009.
%
% Input:
%   x        input signal (assumed complex)
%   L        model order (default: ESTER)
%   M        covariance matrix size (default: length(x)/2)
%   
% Output:
%   w        set of L unconstrained frequencies
%
% Description:
%   Estimates a set of unconstrained frequencies using the shift-invariance
%   property of the signal subspace, i.e., using the ESPRIT subspace
%   method. If no order is given as input, it finds the model order using
%   the ESTER principle, i.e., order_shiftinv() is used.
%
% Example:
%   L=freq_shiftinv(x,10,length(x)/2);
%
% Implemented By:
%   Mads G. Christensen (mgc@es.aau.dk)
%
if nargin<3,M=floor(length(x)/2);end
if nargin<2,L=order_shiftinv(x,M);end
N=length(x);
x=x(:);
R=covm(x,M);
[U,S,V]=svd(R);
A=U(1:M-1,1:L);
B=U(2:M,1:L);
P=(A\B);
[C,E]=eig(P);
w=mod(angle(diag(E)),2*pi);
w=w(:);