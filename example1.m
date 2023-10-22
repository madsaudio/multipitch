% This is an example of how to use the Multi-Pitch Estimation Toolbox for
% single-pitch estimation on a synthetically generated signal with typical 
% settings.

% Set experimental conditions
addpath pitch/
N=200;
L=5;
w0=2*pi*0.12311;
PSNR=20;
M=floor(N/2);
w0_lim=2*pi*[0.02 0.3];
F=16384;
a=ones(L,1).*exp(j*2*pi*rand(L,1));

% Generate noisy signal
Z=vandermonde(w0*[1:L],N);
A=sum(([1:L]'.^2).*(abs(a).^2));
e=randn(N,1)+j*randn(N,1);
e=e./sqrt(var(e))*sqrt(10^(PSNR/-10)*A);
x=Z*a+e;

% Run joint pitch/order estimators
joint_nls(x)
joint_anls(x,w0_lim,F)
joint_orth(x,w0_lim,F,M)
joint_optfilt(x,w0_lim,F,M/2)
joint_hmp(x,w0_lim,F)
joint_em(x,w0_lim,F)

% Run pitch estimators
pitch_wls(x,freq_shiftinv(x,L,M))
pitch_shiftinv(x,L,w0_lim,F,M)
pitch_orth(x,L,w0_lim,F,M)
pitch_anls(x,L,w0_lim,F)
pitch_optfilt(x,L,w0_lim,F,M/2)

% Run order estimators
order_eig(x,M)
order_map(x,F)
order_orth(x,F,M)
order_shiftinv(x,M)

% Run unconstrained frequency estimators
sort(freq_anls(x))
sort(freq_optfilt(x))
sort(freq_shiftinv(x))
sort(freq_orth(x))

% Run amplitude estimators
amp_als(x,w0*[1:L],F)
amp_ls(x,w0*[1:L])
amp_wls(x,w0*[1:L],covm(e,M/2))
amp_capon(x,w0*[1:L],M/2)
amp_apes(x,w0*[1:L],M/2)

% Run pitch detection
detect(x,pitch_wls(x),order_orth(x))