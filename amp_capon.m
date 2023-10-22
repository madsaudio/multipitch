function a=amp_capon(x,w,M);
%AMP_CAPON   Complex amplitude estimator based on Capon
%
% Syntax:
%   a=amp_capon(x,w,M);
%
% About:
%   This file is part of the Multi-Pitch Estimation Toolbox for the book
%   M. G. Christensen and A. Jakobsson, Multi-Pitch Estimation, Morgan &
%   Claypool Publishers, 2009.
%
% Input:
%   x        input signal (assumed complex)
%   w        vector of L frequencies
%   M        filter length (default: length(x)/4)
%
% Output:
%   a        vector of L complex amplitudes for the frequencies in w
%
% Description:
%   The function estimates the complex amplitudes associated with a set of
%   frequencies. It does so using the Capon method described in Section
%   5.3 of Christensen and Jakobsson (2009). More specifically, it uses the
%   extended Capon amplitude estimator.
%
% Example:
%   a=amp_capon(x,w0*[1:L],length(x)/4);
%
% Implemented By:
%   Mads G. Christensen (mgc@es.aau.dk)
%
N=length(x);
if nargin<3,M=floor(N/4);end
R=covm(x,M);
a=amp_wls(x,w,R);
