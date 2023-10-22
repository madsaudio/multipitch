function w=freq_anls(x,L,F)
%FREQ_ANLS   Unconstrained frequency estimator based on approximate NLS
%
% Syntax:
%   w=freq_anls(x,L,F);
%
% About:
%   This file is part of the Multi-Pitch Estimation Toolbox for the book
%   M. G. Christensen and A. Jakobsson, Multi-Pitch Estimation, Morgan &
%   Claypool Publishers, 2009.
%
% Input:
%   x        input signal (assumed complex)
%   L        model order (default: MAP estimate)
%   F        FFT size (default: 4096)
%   
% Output:
%   w        set of L unconstrained frequencies
%
% Description:
%   The function provides a set of frequency estimates for the 
%   unconstrained model. It uses a periodogram (hence approximate nonlinear
%   least-squares) for finding the frequencies, i.e., an FFT. If no order 
%   is supplied as input, it uses periodogram estimates of the frequencies 
%   in combination with the MAP rule for determining the number of complex 
%   sinusoids, i.e., it uses order_map(). The method is essentially 
%   identical to the harmonic summation method.
%
% Example:
%   L=freq_anls(x,10,4096);
%
% Implemented By:
%   Mads G. Christensen (mgc@es.aau.dk)
%
if nargin<3,F=4096;end
if nargin<2,L=order_map(x,F);end;
N=length(x);
X=fft(x,F)/N;
f=findpeaks(abs(X),L,1);
w=2*pi*(f-1)/F;