function [w0,L]=joint_nls(x,w0_lim,F);
%JOINT_NLS   Joint single-pitch/order estimator based on exact NLS and MAP
%
% Syntax:
%   [w0,L]=joint_nls(x,w0_lim,F);
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
%   complex signal based on the exact nonlinear least-squares method in 
%   Section 2.4 of Christensen and Jakobsson (2009) combined with the MAP 
%   criterion for order estimation described in Section 2.6. The method can
%   be applied to multi-pitch signals without any further modifications,
%   but returns only one pitch estimate.
%
% Example:
%   [w0,L]=joint_nls(x,2*pi*[1e-3 1e-1],F);
%
% Implemented By:
%   Mads G. Christensen (mgc@es.aau.dk)
%
if nargin<3,F=4096;end
if nargin<2 w0_lim=2*pi*[1e-3 0.25];end
freqndx=[round(w0_lim(1)/2/pi*F)+1:1:round(w0_lim(2)/2/pi*F)+1];
w0_set=2*pi*(freqndx-1)/F;
N=length(x);
x=x(:);
L_max=floor(N/4);
J=zeros(size(w0_set));
L=zeros(size(w0_set));
h=1;
for w0=w0_set,
    L_w0=min([floor(2*pi/w0) L_max]);
    P=zeros(L_w0,1);
    C=zeros(L_w0,1);
    l=1;
    while (l<=L_w0)
        Z=vandermonde(w0*[1:l],N);
        P(l)=var(x-Z*pinv(Z)*x);
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
%if not(detect(x,w0,L)),
%    L=0;
%    w0=0;
%end
