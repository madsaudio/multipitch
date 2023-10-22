% This is a demo of the Multi-Pitch Estimation Toolbox applied to
% multi-pitch estimation on a synthetically generated signal.

% Set experimental conditions
addpath pitch/
N=200;
L=[5; 3];
w0=2*pi*[0.03; 0.25];
M=floor(N/2);
w0_lim=2*pi*[0.02 0.3];
F=16384;

% Generate noisy signal
x=zeros(N,1);
for k=1:length(w0),
    %a=ones(L(k),1).*exp(j*2*pi*rand(L(k),1));
    a=randn(L(k),1)+j*randn(L(k),1);
    Z=vandermonde(w0(k)*[1:L(k)],N);
    x=x+Z*a;
end
e=randn(N,1)+j*randn(N,1);
e=e./sqrt(var(e))*1e-3;
x=x+e;

% Run joint pitch/order estimators
joint_hmp(x,w0_lim,F)
joint_em(x,w0_lim,F)
