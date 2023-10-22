function [w0,L]=joint_orth(x,w0_lim,F,M);
%JOINT_ORTH   Joint single-pitch/order estimator using subspace orthogonality
%
% Syntax:
%   [w0,L]=joint_orth(x,w0_lim,F,M);
%
% About:
%   This file is part of the Multi-Pitch Estimation Toolbox for the book
%   M. G. Christensen and A. Jakobsson, Multi-Pitch Estimation, Morgan &
%   Claypool Publishers, 2009.
%
% Input:
%   x        input signal (assumed complex)
%   w0_lim   search inverval limits for w0 (default: 2*pi*[1e-3 0.25])
%   F        FFT size (default: 4096)
%   M        subvector length (default: length(x))
%
% Output:
%   w0       fundamental frequency estimate in radians
%   L        model order, i.e., number of harmonics
%
% Description:
%   Estimates the fundamental frequency and model order of a single-pitch 
%   complex signal based on the subspace method described in Section 4.7 of
%   Christensen and Jakobsson (2009), which is based on subspace 
%   orthogonality. Use detect() for pitch detection. Cannot be applied
%   directly to multi-pitch signals.   
% 
% Example:
%   [w0,L]=joint_orth(x,2*pi*[1e-3 1e-1],F,length(x)/2);
%
% Implemented By:
%   Mads G. Christensen (mgc@es.aau.dk)
%
if nargin<4,M=floor(length(x)/2);end
if nargin<3,F=4096;end
if nargin<2 w0_lim=2*pi*[1e-3 0.25];end
N=length(x);
freqndx=[round(w0_lim(1)/2/pi*F)+1:1:round(w0_lim(2)/2/pi*F)+1];
w0_set=2*pi*(freqndx-1)/F;
x=x(:);
R=covm(x,M);
[U,S,V] = svd(R);
clear S V;
L_max=M-2;
T=abs(fft(U,F)).^2;
P=zeros(length(freqndx),L_max);
J=zeros(length(freqndx),1);
L=zeros(length(freqndx),1);
l=1;
for f=freqndx,
    K_set=1:min([floor((F-1)/(f-1)) L_max]);
    k=1;
    for K=K_set,
        s=M*min([K (M-K)]);
        P(l,k)=s./(sum(sum(T((f-1)*[1:K]+1,K+1:M)))+eps);
        k=k+1;
    end
    [J(l),ndx]=max(P(l,:));
    L(l)=K_set(ndx);
    l=l+1;
end
[tmp,ndx]=max(J);
w0=2*pi*(freqndx(ndx)-1)/F;
L=L(ndx);
