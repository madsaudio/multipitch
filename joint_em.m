function [w0,L]=joint_em(x,w0_lim,F,estimator);
%JOINT_EM   Joint multi-pitch/order estimator based on the EM algorithm
%
% Syntax:
%   [w0,L]=joint_em(x,w0_lim,F);
%
% About:
%   This file is part of the Multi-Pitch Estimation Toolbox for the book
%   M. G. Christensen and A. Jakobsson, Multi-Pitch Estimation, Morgan &
%   Claypool Publishers, 2009.
%
% Input:
%   x        multi-/single-pitch input signal (assumed complex)
%   w0_lim   search inverval limits for w0 (default: 2*pi*[1e-3 0.25])
%   F        FFT size (default: 2^12)
%
% Output:
%   w0       fundamental frequency estimate in radians
%   L        model orders corresponding to the entries in w0
%
% Description:
%   Estimates the fundamental frequencies and number of harmonics of a
%   multi-pitch complex signal based on the expectation maximization (EM)
%   algorithm as described in Section 2.8 of Christensen and Jakobsson
%   (2009). The method is essentially a recursive application of the
%   approximate NLS. It uses the MAP criterion for determining the number 
%   of harmonics and continues to extract sources until the pitch detection 
%   criterion of Section 2.6 fails. The EM algorithm is initialized using 
%   the harmonic matching pursuit. Note that the behavior of the algorithm 
%   is complicated and is highly dependent on the initialization.
%
% Example:
%   [w0,L]=joint_em(x,2*pi*[0.001 0.1],F);
%
% Implemented By:
%   Mads G. Christensen (mgc@es.aau.dk)
%
if nargin<4,estimator{2}=@amp_ls;estimator{1}=@joint_anls;end
if nargin<3,F=2^12;end
if nargin<2,w0_lim=2*pi*[1e-3 0.25];end
x=x(:);
N=length(x);
K=100;
[w0,L]=joint_hmp(x,w0_lim,F,estimator);
Q=length(w0);
C=zeros(max(L),Q);
for q=1:Q,
    C(1:L(q),q)=estimator{2}(x,w0(q)*[1:L(q)]);
end
e=zeros(size(x));
b=ones(Q,1)/Q;
Y=zeros(N,Q);
delta=inf;
k=1;
while delta>1e-4 & k<K,
    w0_old=w0;
    z=zeros(N,1);
    for q=1:Q;
        Z=vandermonde(w0(q)*[1:L(q)],N);
        z=z+Z*C(1:L(q),q);
    end
    e=x-z;
    for q=1:Q,
        Z=vandermonde(w0(q)*[1:L(q)],N);
        Y(:,q)=Z*C(1:L(q),q)+b(q)*e;
    end
    for q=1:Q,
        [w0(q),L(q)]=estimator{1}(Y(:,q),w0_lim,F);
        Z=vandermonde(w0(q)*[1:L(q)],N);
        C(1:L(q),q)=estimator{2}(Y(:,q),w0(q)*[1:L(q)]);
    end
    k=k+1;
    delta=abs(mean(w0-w0_old));
end
