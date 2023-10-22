function [L]=order_map(x,F);
%ORDER_MAP   Order estimator based on the MAP criterion
%
% Syntax:
%   L=order_map(x,F);
%
% About:
%   This file is part of the Multi-Pitch Estimation Toolbox for the book
%   M. G. Christensen and A. Jakobsson, Multi-Pitch Estimation, Morgan &
%   Claypool Publishers, 2009.
%
% Input:
%   x        input signal (assumed complex)
%   F        FFT size (default: 4096)
%   
% Output:
%   L        estimated signal subspace rank/number of sinusoids
%
% Description:
%   Function that determines the number of complex sinusoids using the MAP
%   criterion for white Gaussian noise as described in Section 2.6 of 
%   Christensen and Jakobsson (2009). It uses periodogram estimates of 
%   sinusoidal frequencies and amplitudes for estimating the noise variance
%   as a function of the model order. The sinusoids are not assumed to be 
%   harmonically related.
%
% Example:
%   L=order_map(x,8192);
%
% Implemented By:
%   Mads G. Christensen (mgc@es.aau.dk)
%
if nargin<2,F=4096;end
N=length(x);
x=x(:);
L_set=[1:floor(N/4)];
X=fft(x,F)/N;
f=findpeaks(abs(X),max(L_set),1);
w_all=2*pi*(f-1)/F;
c_all=X(f);
for l=1:length(L_set);
    L=L_set(l);
    w=w_all(1:L);
    c=c_all(1:L);
    Z=vandermonde(w',N);
    P(l)=N*log(var(x-Z*c))+5*l*log(N)/2;
end
[tmp,l]=min(P);
L=L_set(l);
Z=vandermonde(w(1:L),N);
c=c_all(1:L);
if N*log(var(x-Z*c))+L*log(N)+3/2*log(N)>N*log(var(x)),
    L=0;
end