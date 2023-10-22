function R=covm(x,M,METHOD);
% COVM   Covariance matrix estimator
%
% Syntax:
%   R=covm(x,M,METHOD);
%
% About:
%   This file is part of the Multi-Pitch Estimation Toolbox for the book
%   M. G. Christensen and A. Jakobsson, Multi-Pitch Estimation, Morgan &
%   Claypool Publishers, 2009.
%
% Input:
%   x        input signal (assumed complex)
%   M        subvector length (default: length(x)/2)
%   METHOD   'forward'/'backward'/'combined' (default: 'forward')
%
% Output:
%   R        M-by-M covariance matrix
%
% Description:
%   Auxiliary function that estimates the covariance matrix in various
%   ways. For details see Sections 1.5, 3.3, and 5.3 of  of Christensen and
%   Jakobsson (2009) 
%
% Example:
%   x=randn(N,1);
%   R=covm(x,N/2,'backward');
%
% Implemented By:
%   Mads G. Christensen (mgc@es.aau.dk)
%
if nargin<3,METHOD='forward';end
if nargin<2,M=floor(length(x)/2);end
x=x(:);
N=length(x);
R=zeros(M,M);
switch METHOD
    case 'forward'
        for i=M:N,
            R=R+x(i-M+1:1:i)*x(i-M+1:1:i)'/(N-M+1);
        end
    case 'backward'
        for i = M:N,
            R=R+x(i:-1:i-M+1)*x(i:-1:i-M+1)'/(N-M+1);
        end
    case 'combined'
        for i = M:N,
            R=R+x(i:-1:i-M+1)*x(i:-1:i-M+1)'/(N-M+1);
        end
        R=(R+fliplr(eye(M))*R.'*fliplr(eye(M)))/2;
end
