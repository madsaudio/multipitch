function H=detect(x,w0,L,a);
%DETECT   Pitch detection using the MAP criterion
%
% Syntax:
%   H=detect(x,w0,L,a);
%
% About:
%   This file is part of the Multi-Pitch Estimation Toolbox for the book
%   M. G. Christensen and A. Jakobsson, Multi-Pitch Estimation, Morgan &
%   Claypool Publishers, 2009.
%
% Input:
%   x        input signal (assumed complex)
%   w0       candidate fundamental frequency
%   L        candidate model order (default: MAP estimate)
%   a        set of L complex amplitudes (default: least-squares)
%   
% Output:
%   H        decision - 1 if a pitch is detected, 0 otherwise
%
% Description:
%   Function that determines whether a pith is present given the observed
%   complex signal, a candidate fundamental frequency and a model order
%   (and optionally a set of amplitudes). It uses the MAP order estimation 
%   criterion described in Section 2.6 of Christensen and Jakobsson (2009).
%   The criterion is based on the noise being white and Gaussian 
%   distributed and a single-pitch model.
%
% Example:
%   [w0,L]=joint_optfilt(x,2*pi*[0:F/8-1]/F,length(x)/2);
%   H=detect(x,w0,L);
%
% Implemented By:
%   Mads G. Christensen (mgc@es.aau.dk)
%
if nargin<3,L=order_map(x);end
if nargin<2,w0=pitch_wls(x);end
if nargin<4,a=amp_ls(x,w0*[1:L]);end
x=x(:);
N=length(x);
Z=vandermonde(w0*[1:L],N);
if N*log(var(x-Z*a))+L*log(N)+3/2*log(N)<N*log(var(x)),
    H=1;
else 
    H=0;
end