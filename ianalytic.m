function xr=ianalytic(xc);
%IANALYTIC   Computes the real signal of an analytic signal
%
% Syntax:
%   y=ianalytic(x);
%
% About:
%   This file is part of the Multi-Pitch Estimation Toolbox for the book
%   M. G. Christensen and A. Jakobsson, Multi-Pitch Estimation, Morgan &
%   Claypool Publishers, 2009.
%
% Input:
%   x        input signal (assumed to be discrete-time analytic signal)
%
% Output:
%   y        real signal
%
% Description:
%   Auxiliary function that computes a real signal from a discrete-time a
%   analytic signal, i.e., maps a complex signal to a real signal of double
%   length as described in Appendix A of Christensen and Jakobsson (2009).
%
% Example:
%   y=ianalytic(x)
%
% Implemented By:
%   Mads G. Christensen (mgc@es.aau.dk)
%
xc=xc(:);
X=fft(xc);
X=[X; zeros(size(X))];
xr=2*real(ifft(X));
