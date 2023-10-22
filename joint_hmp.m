function [w0,L]=joint_hmp(x,w0_lim,F,estimator);
%JOINT_HMP   Joint multi-pitch/order estimator based on harmonic matching pursuit
%
% Syntax:
%   [w0,L]=joint_hmp(x,w0_lim,F);
%
% About:
%   This file is part of the Multi-Pitch Estimation Toolbox for the book
%   M. G. Christensen and A. Jakobsson, Multi-Pitch Estimation, Morgan &
%   Claypool Publishers, 2009.
%
% Input:
%   x        multi-/single-pitch input signal (assumed complex)
%   w0_lim   search inverval limits for w0 (default: 2*pi*[1e-3 0.25])
%   F        FFT size (default: 2^12)
%
% Output:
%   w0       fundamental frequency estimate in radians
%   L        model orders corresponding to the entries in w0
%
% Description:
%   Estimates the fundamental frequencies and number of harmonics of a
%   multi-pitch complex signal based on the harmonic matching pursuit
%   method described in Section 2.9 of Christensen and Jakobsson (2009).
%   The method is essentially a recursive application of the approximate
%   NLS or harmonic summation method. It uses the MAP criterion for
%   determining the number of harmonics and continues to extract sources
%   until the pitch detection criterion of Section 2.6 fails.
%
% Example:
%   [w0,L]=joint_hmp(x,2*pi*[0.001 0.1],F);
%
% Implemented By:
%   Mads G. Christensen (mgc@es.aau.dk)
%
if nargin<4,estimator{2}=@amp_ls;estimator{1}=@joint_anls;end
if nargin<3,F=2^12;end
if nargin<2,w0_lim=2*pi*[1e-3 0.25];end
x=x(:);
N=length(x);
r=x;
R=fft(r,F);
[w0(1),L(1)]=estimator{1}(x,w0_lim,F);
a=estimator{2}(x,w0(1)*[1:L(1)]);
Z=vandermonde(w0(1)*[1:L(1)],N);
k=2;
H=detect(r,w0(1),L(1),a);
while H,
    r=r-Z*a;
    R=fft(r,F);
    [w0(k),L(k)]=estimator{1}(r,w0_lim,F);
    if L(k)>0
        Z=vandermonde(w0(k)*[1:L(k)],N);
        a=estimator{2}(x,w0(k)*[1:L(k)]);
        H=detect(r,w0(k),L(k),a);
    else
        H=0;
    end
    k=k+1;
end
w0=w0(:);
L=L(:);
w0=w0(1:end-1);
L=L(1:end-1);