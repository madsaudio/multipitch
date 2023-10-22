function a=amp_ls(x,w);
%AMP_LS   Complex amplitude estimator based on LS
%
% Syntax:
%   a=amp_ls(x,w);
%
% About:
%   This file is part of the Multi-Pitch Estimation Toolbox for the book
%   M. G. Christensen and A. Jakobsson, Multi-Pitch Estimation, Morgan &
%   Claypool Publishers, 2009.
%
% Input:
%   x        input signal (assumed complex)
%   w        vector of L frequencies
%   
% Output:
%   a        vector of L complex amplitudes for the frequencies in w
%
% Description:
%   The function estimates the complex amplitudes associated with a set of
%   frequencies. It does so using the least-squares method as described in
%   Section 5.2 of Christensen and Jakobsson (2009).
%
% Example:
%   a=amp_ls(x,w0*[1:L]);
%
% Implemented By:
%   Mads G. Christensen (mgc@es.aau.dk)
%
x=x(:);
N=length(x);
Z=vandermonde(w,N);
a=Z\x;