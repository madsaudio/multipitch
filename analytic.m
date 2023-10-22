function xc=analytic(xr);
%ANALYTIC   Computes the discrete-time downsampled analytic signal
%
% Syntax:
%   y=analytic(x);
%
% About:
%   This file is part of the Multi-Pitch Estimation Toolbox for the book
%   M. G. Christensen and A. Jakobsson, Multi-Pitch Estimation, Morgan &
%   Claypool Publishers, 2009.
%
% Input:
%   x        input signal (assumed real)
%
% Output:
%   y        dicrete-time analytic (complex) signal approximation
%
% Description:
%   Auxiliary function that computes the downsampled discrete-time analytic
%   signal,i.e., maps a real signal to a complex signal of half length 
%   using the Hilbert transform as described in Appendix A of Christensen 
%   and Jakobsson (2009).
%
% Example:
%   y=analytic(x)
%
% Implemented By:
%   Mads G. Christensen (mgc@es.aau.dk)
%
xr=xr(:);
X=fft(xr);
X(1)=X(1)./2;
X(end/2+1)=X(end/2+1)/2;
X(end/2+2:end)=0;
xc=2*ifft(X);
xc=xc(1:2:end-1);
xc=xc(:);