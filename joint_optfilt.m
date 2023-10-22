function [w0,L]=joint_optfilt(x,w0_lim,F,M)
%JOINT_OPTFILT   Joint single-pitch/order estimator using optimal filtering
%
% Syntax:
%   [w0,L]=joint_optfilt(x,w0_lim,F,M);
%
% About:
%   This file is part of the Multi-Pitch Estimation Toolbox for the book
%   M. G. Christensen and A. Jakobsson, Multi-Pitch Estimation, Morgan &
%   Claypool Publishers, 2009.
%
% Input:
%   x        input signal (assumed complex)
%   w0_set   candiate w0 set in radians (default: (0;pi/2])
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
%   complex signal based on the optimal filtering method described in
%   Sections 3.5 and 3.8 of Christensen and Jakobsson (2009) and 
%   implemented using the order-recursive expressions in Section 3.9. It 
%   uses the MAP criterion in Section 2.6 for determining the order under 
%   the assumption that the noise is white Gaussian. The method can be 
%   applied to multi-pitch signals also without any further modifications 
%   but estimates only the most prominent pitch. Note that the method may 
%   become numerically unstable for high SNRs. It may also require that a 
%   very dense candidate w0 set is used.
%
% Example:
%   [w0,L]=joint_optfilt(x,2*pi*[1e-3 1e-1],4096,length(x)/4);
%
% Implemented By:
%   Mads G. Christensen (mgc@es.aau.dk)
%
if nargin<4,M=floor(length(x)/4);end
if nargin<3,F=4096;end
if nargin<2 w0_lim=2*pi*[1e-3 0.25];end
L_max=floor(M/2);
freqndx=[round(w0_lim(1)/2/pi*F)+1:1:round(w0_lim(2)/2/pi*F)+1];
w0_set=2*pi*(freqndx-1)/F;
N=length(x);
x=x(:);
J=zeros(size(w0_set));
L=zeros(size(w0_set));
R=covm(x,M);
Q=inv(R);
h=1;
for w0=w0_set,
    L_w0=min([floor(2*pi/w0) L_max]);
    P=zeros(L_w0,1);
    C=zeros(L_w0,1);
    l=1;
    while (l<=L_w0)
        Z=vandermonde(w0*[1:l-1],M);
        z=vandermonde(w0*l,M);
        u=Q*z;
        if l==1,
            V=1./(z'*u);
        else
            y=Z'*u;
            v=V*y;
            V=[V zeros(l-1,1); zeros(1,l)]+1/(z'*u-y'*v)*[v*v' -v; -v' 1];
        end
        P(l)=var(x)-sum(sum(real(V)));
        C(l)=N*log(abs(P(l)))+l*log(N)+3/2*log(N);
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