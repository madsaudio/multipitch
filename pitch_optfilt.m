function [w0]=pitch_optfilt(x,L,w0_lim,F,M);
%PITCH_OPTFILT   Multi-pitch estimator based on optimal filtering
%
% Syntax:
%   w0=pitch_optfilt(x,L,w0_lim,F,M);
%
% About:
%   This file is part of the Multi-Pitch Estimation Toolbox for the book
%   M. G. Christensen and A. Jakobsson, Multi-Pitch Estimation, Morgan &
%   Claypool Publishers, 2009.
%
% Input:
%   x        multi-/single-pitch input signal (assumed complex)
%   L        vector containg the number of harmonics for each source
%   w0_lim   search inverval limits for w0 (default: 2*pi*[1e-3 0.25])
%   F        FFT size (default: 4096)
%   M        subvector length (default: length(x)/4)
%
% Output:
%   w0       fundamental frequency estimate in radians
%
% Description:
%   Estimates the fundamental frequencies of a single- or multi-pitch
%   complex signal given the number of harmonics for each source based on
%   the optimal filtering method described in Section 3.5 of Christensen 
%   and Jakobsson (2009). The number of sources is determined from the 
%   length of the order vector.
%
% Example:
%   [w0]=pitch_optfilt(x,[5 4 10],2*pi*[0.001 0.1],4096,length(x)/4);
%
% Implemented By:
%   Mads G. Christensen (mgc@es.aau.dk)
%
if nargin<5,M=floor(length(x)/4);end
if nargin<4,F=4096;end
if nargin<3 w0_lim=2*pi*[1e-3 0.25];end
x=x(:);
freqndx=[round(w0_lim(1)/2/pi*F)+1:1:round(w0_lim(2)/2/pi*F)+1];
w0_set=2*pi*(freqndx-1)/F;
R=covm(x,M);
C=inv(R);
L=L(:)';
[L,ndx]=sort(L);
L_set=unique(L);
f0=[];
for l=L_set,
    K=length(find(l==L));
    k=1;
    P=zeros(size(w0_set(w0_set*l<=2*pi)));
    for w0=w0_set(w0_set*l<=2*pi);
        Z=vandermonde(w0*[1:l],M);
        P(k)=real(sum(sum(inv(Z'*C*Z))));
        k=k+1;
    end
    [f,p]=findpeaks(P,K,1,f0);
    f0=[f0; f];
end
w0(ndx,1)=w0_set(f0);