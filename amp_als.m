function a=amp_als(x,w,F);
%AMP_ALS   Complex amplitude estimator based on approximate LS
%
% Syntax:
%   a=amp_als(x,w);
%
% About:
%   This file is part of the Multi-Pitch Estimation Toolbox for the book
%   M. G. Christensen and A. Jakobsson, Multi-Pitch Estimation, Morgan &
%   Claypool Publishers, 2009.
%
% Input:
%   x        input signal (assumed complex)
%   w        vector of L frequencies
%   F        FFT size
%   
% Output:
%   a        vector of L complex amplitudes for the frequencies in w
%
% Description:
%   The function estimates the complex amplitudes associated with a set of
%   frequencies. It does so using the approximate least-squares method as 
%   described in Section 5.2 of Christensen and Jakobsson (2009), i.e., 
%   using the FFT.
%
% Example:
%   a=amp_als(x,w0*[1:L]);
%
% Implemented By:
%   Mads G. Christensen (mgc@es.aau.dk)
%
if nargin<3,F=2^12;end
N=length(x);
x=x(:);
X=fft(x,F);
a=X(round(w/2/pi*F)+1)/N;