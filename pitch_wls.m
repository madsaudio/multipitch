function w0 = pitch_wls(x,w,a)
%PITCH_WLS   Multi-pitch estimator based on the WLS method
%
% Syntax:
%   w0=pitch_wls(x,w,a);
%
% About:
%   This file is part of the Multi-Pitch Estimation Toolbox for the book
%   M. G. Christensen and A. Jakobsson, Multi-Pitch Estimation, Morgan &
%   Claypool Publishers, 2009.
%
% Input:
%   x        single-pitch input signal (assumed complex)
%   w        set of L unconstrained frequencies (default: ESPRIT)
%   a        complex L amplitudes (default: least-squares)
%
% Output:
%   w0       fundamental frequency estimate in radians
%
% Description:
%   Estimates the fundamental frequency and model order of a single-pitch
%   complex signal based on the weighted least-squares (aka harmonic 
%   fitting) method described in Section 2.10 of Christensen and Jakobsson 
%   (2009). It uses unconstrained frequencies estimated using ESPRIT 
%   (unless frequencies are provided as input) and estimates amplitudes 
%   using least-squares. The model order is determined from the length of 
%   the input vectors, or using order_shiftinv() if no parameters are 
%   given.
%
% Example:
%   w0=pitch_wls(x,freq_anls(x,L));
%
% Implemented By:
%   Mads G. Christensen (mgc@es.aau.dk)
%
N=length(x);
x=x(:);
if nargin<2,L=order_shiftinv(x,floor(N/2));w=freq_shiftinv(x,L,floor(N/2));end
if nargin<3,Z=vandermonde(w,N);a=Z\x;end    
L=length(w);
[w,i]=sort(mod(w(:),2*pi));
a=a(i);
A=abs(a(:)).^2;
w0=sum([1:L]'.*A.*w)/sum([1:L]'.^2.*A);
