% Multi-Pitch Estimation Toolbox
% Version 1.0 29-Jun-2009
%
% This is the README file for the Multi-Pitch Estimation toolbox for MATLAB
% for the book M. G. Christensen and A. Jakobsson, Multi-Pitch Estimation, 
% Morgan & Claypool Publishers, 2009. The toolbox is a collection of MATLAB 
% scripts and functions for estimating the parameters of periodic signals. 
% For clarity, we have opted to keep the files a simple as possible while 
% sacrificing computational efficiency and performance; the point of these 
% files is to illustrate the various estimation principles. The users may 
% wish to modify the files to fit their needs. We ask that you please cite 
% our book and/or the appropriate papers when using these methods.  If you 
% find any bugs or errors or have any suggestions or comments, please send 
% Mads G. Christensen an email at mgc@es.aau.dk. The toolbox contains the 
% MATLAB functions listed below. For help on how to use the individual 
% functions, use MATLAB's help function and/or see the example files.
%
% Table of Contents (T0C)
% -----------------------
%   joint_optfilt.m      - Joint single-pitch/order estimator using optimal filtering
%   joint_orth.m         - Joint single-pitch/order estimator using subspace orthogonality
%   joint_nls.m          - Joint single-pitch/order estimator based on exact NLS and MAP
%   joint_anls.m         - Joint single-pitch/order estimator based on approximate NLS and MAP
%   joint_em.m           - Joint multi-pitch/order estimator based on the EM algorithm
%   pitch_hmp.m          - Joint multi-pitch/order estimator based on harmonic matching pursuit
%   order_orth.m         - Order estimator based on subspace orthogonality
%   order_shiftinv.m     - Order estimator based on subspace shift-invarinace
%   order_eig.m          - Order estimator based on the eigenvalues and MDL
%   order_map.m          - Order estimator based on the MAP criterion
%   pitch_shiftinv.m     - Multi-pitch estimator based on subspace shift-invariance
%   pitch_wls.m          - Multi-pitch estimator based on the WLS method
%   pitch_orth.m         - Multi-pitch estimator based on subspace orthogonality
%   pitch_optfilt.m      - Multi-pitch estimator based on optimal filtering
%   pitch_anls.m         - Multi-pitch estimator based on approximate NLS
%   freq_orth.m          - Unconstrained frequency estimator based on MUSIC
%   freq_optfilt.m       - Unconstrained frequency estimator based on Capon
%   freq_shiftinv.m      - Unconstrained frequency estimator based on ESPRIT
%   freq_anls.m          - Unconstrained frequency estimator based on approximate NLS
%   amp_ls.m             - Complex amplitude estimator based on LS
%   amp_als.m            - Complex amplitude estimator based on approximate LS
%   amp_wls.m            - Complex amplitude estimator based on WLS
%   amp_capon.m          - Complex amplitude estimator based on Capon
%   amp_apes.m           - Complex amplitude esitmator based on APES
%   detect.m             - Pitch detection using the MAP criterion
%   covm.m               - Covariance matrix estimator
%   analytic.m           - Computes the discrete-time downsampled analytic signal
%   ianalytic.m          - Computes the real signal of an analytic signal
%   findpeaks.m          - Finds peaks or valleys of non-negative data
%   vandermonde.m        - Constructs Vandermonde matrix from frequencies
%
% See also the included example files.
%