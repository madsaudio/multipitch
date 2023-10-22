function w=freq_optfilt(x,L,F,M);
%FREQ_OPTFILT   Unconstrained frequency estimator based on Capon
%
% Syntax:
%   w=freq_optfilt(x,L,F,M);
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
%   M        filter length (default: length(x)/4)
%   
% Output:
%   w        set of L unconstrained frequencies
%
% Description:
%   Unconstrained frequency estimator based on the Capon spetral estimator.
%   If no model order is specified, then the MAP rule is used for finding
%   the model order, i.e., order_map().
%
% Example:
%   L=freq_optfilt(x,10,4096,length(x)/4);
%
% Implemented By:
%   Mads G. Christensen (mgc@es.aau.dk)
%
if nargin<4,M=floor(length(x)/4);end
if nargin<3,F=4096;end
if nargin<2,L=order_map(x,F);end;
x=x(:);
R=covm(x,M,'backward');
Q=inv(R);
m=[0:M-1]';
w_set=2*pi*[0:F-1]/F;
f=1;
for w=w_set,
   a=exp(-j*w*m);
   J(f)=1./real(a'*Q*a);
   f=f+1;
end
f=findpeaks(J,L,1);
w=2*pi*(f-1)/F;