function w=freq_orth(x,L,F,M);
%FREQ_ORTH   Unconstrained frequency estimator based on MUSIC
%
% Syntax:
%   w=freq_orth(x,L,F,M);
%
% About:
%   This file is part of the Multi-Pitch Estimation Toolbox for the book
%   M. G. Christensen and A. Jakobsson, Multi-Pitch Estimation, Morgan &
%   Claypool Publishers, 2009.
%
% Input:
%   x        input signal (assumed complex)
%   L        model order
%   F        FFT size (default: 4096)
%   M        covariance matrix size (default: length(x)/2)
%   
% Output:
%   w        set of L unconstrained frequencies
%
% Description:
%   Estimates the unconstrained frequencies of sinusoids in the observed
%   signal based on the subspace orthogonality property, i.e., using MUSIC.
%   If no model order is given as input, it uses the angles between
%   subspaces to determine the model order, i.e., it uses order_orth().
%
% Example:
%   L=freq_orth(x,10,4096,length(x)/2);
%
% Implemented By:
%   Mads G. Christensen (mgc@es.aau.dk)
%
if nargin<4,M=floor(length(x)/2);end
if nargin<3,F=4096;end
if nargin<2,L=order_orth(x,F,M);end;
N=length(x);
x=x(:);
R=covm(x,M);
[U,D,V]=svd(R);
P=1./sum(abs(fft(U(:,L+1:M),F)).^2,2);
[f,p]=findpeaks(P,L,1);
w=2*pi*(f-1)/F;