function [w0]=pitch_shiftinv(x,L,w0_lim,F,M);
%PITCH_SHIFTINV   Single-pitch estimator based on subspace shift-invariance
%
% Syntax:
%   w0=pitch_shiftinv(x,L,w0_lim,F,M);
%
% About:
%   This file is part of the Multi-Pitch Estimation Toolbox for the book
%   M. G. Christensen and A. Jakobsson, Multi-Pitch Estimation, Morgan &
%   Claypool Publishers, 2009.
%
% Input:
%   x        single-pitch input signal (assumed complex)
%   L        number of harmonics
%   w0_lim   search inverval limits for w0 (default: 2*pi*[1e-3 0.25])
%   F        FFT size (default: 4096)
%   M        covariance matrix size (default: length(x)/2)
%
% Output:
%   w0       fundamental frequency estimate in radians
%
% Description:
%   Estimates the fundamental frequency and model order of a single-pitch
%   complex signal using the ESPRIT-based principle of shift-invariance 
%   described in Section 4.9 of Christensen and Jakobsson (2009). 
%
% Example:
%   w0=pitch_shiftinv(x,5);
%
% Implemented By:
%   Mads G. Christensen (mgc@es.aau.dk)
%
if nargin<5,M=floor(length(x)/2);end
if nargin<4,F=4096;end;
if nargin<3 w0_lim=2*pi*[1e-3 0.25];end
if nargin<2,L=order_shiftinv(x,M);end
N=length(x);
x=x(:);
freqndx=[round(w0_lim(1)/2/pi*F)+1:1:round(w0_lim(2)/2/pi*F)+1];
w0_set=2*pi*(freqndx-1)/F;
w0_set=w0_set(w0_set*L<2*pi);
R=covm(x,M);
[U,S,V]=svd(R);
A=U(1:M-1,1:L);
B=U(2:M,1:L);
P=(A\B);
[C,E]=eig(P);
[r,ndx]=sort(mod(angle(diag(E)),2*pi));
C=C(:,ndx);
h=diag(C'*A'*B*C);
J=zeros(size(w0_set));
for k=1:length(w0_set),
  w=w0_set(k);
  r(k)=exp(j*w);
  d=(r(k).^[1:L])';
  J(k)=-2*real(sum(h.*d));
end
[tmp,ndx]=min(J);
w0=w0_set(ndx);
