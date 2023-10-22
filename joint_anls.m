function [w0,L]=joint_anls(x,w0_lim,F);
%JOINT_ANLS   Joint single-pitch/order estimator based on approximate NLS and MAP
%
% Syntax:
%   [w0,L]=joint_anls(x,w0_lim,F);
%
% About:
%   This file is part of the Multi-Pitch Estimation Toolbox for the book
%   M. G. Christensen and A. Jakobsson, Multi-Pitch Estimation, Morgan &
%   Claypool Publishers, 2009.
%
% Input:
%   x        input signal (assumed complex)
%   w0_lim   search invertval limits for w0 (default: 2*pi*[1e-3 0.25])
%   F        FFT size (default: 4096)
%
% Output:
%   w0       fundamental frequency estimate in radians
%   L        model order, i.e., number of harmonics
%
% Description:
%   Estimates the fundamental frequency and model order of a single-pitch
%   complex signal based on the approximate nonlinear least-squares method 
%   in Section 2.4 combined with the MAP criterion for order estimation 
%   described in Section 2.6 of Christensen and Jakobsson (2009). The 
%   method can be applied to multi-pitch signals without any modifications,
%   but estimates only the most dominant pitch.
%
% Example:
%   [w0,L]=joint_anls(x,2*pi*[1e-3 1e-1],F);
%
% Implemented By:
%   Mads G. Christensen (mgc@es.aau.dk)
%
N=length(x);
x=x(:);
if nargin<3,F=4096;end
if nargin<2,w0_lim=2*pi*[1e-3 0.25];end
X=fft(x,F);
freqndx=[round(w0_lim(1)/2/pi*F)+1:1:round(w0_lim(2)/2/pi*F)+1];
w0_set=2*pi*(freqndx-1)/F;
L_max=floor(N/4);
J=zeros(size(w0_set));
L=zeros(size(w0_set));
h=1;
for w0=w0_set,
    L_w0=min([floor(2*pi/w0)-1 L_max]);
    P=zeros(L_w0,1);
    C=zeros(L_w0,1);
    l=1;
    while (l<=L_w0),
        f=round(w0/2/pi*[1:l]*F)+1;
        c=X(f)/N;
        Z=vandermonde(w0*[1:l],N);
        P(l)=var(x-Z*c);
        C(l)=N*log(P(l))+l*log(N)+3/2*log(N);
        l=l+1;
    end 
    [tmp,ndx]=min(C);
    L(h)=ndx;
    J(h)=C(ndx);
    h=h+1;
end
[tmp,ndx]=min(J);
L=L(ndx);
w0=w0_set(ndx);
if not(detect(x,w0,L)),
    L=0;
end
