% Copyright (c) 2012,2013,2023 Mitsubishi Electric Research Laboratories (MERL).
%
% SPDX-License-Identifier: AGPL-3.0-or-later
%
% Demo function that mimics the experiments of paper
%
% Cédric Févotte, Jonathan Le Roux, John R. Hershey,
% "Non-Negative Dynamical System With Application to Speech and Audio,"
% in Proc. IEEE International Conference on Acoustics, Speech, and Signal
% Processing (ICASSP 2013), May 2013.
%
% The demo does not exactly reproduce the experiments in the paper, which
% uses a large amount of copyrighted data that cannot be distributed with
% the code, but the experimental layout is the same.
%
% The first part of the code executes the function 'nds_train' on clean
% training speech.
%
% The second part runs a denoising a example, as in Section 4.2 of the
% paper.
%
% Note that to allow for fast execution of the demo, the NDS model is
% trained with a short sequence. In a realistic scenario, more data would
% be needed to improve the performance of the denoising system.
%
% The notations in the code closely follow the notations of the paper.

% path to auxiliary functions
addpath ./auxfun

%% Load speech training sample
disp('Training...')
[x,fs] = audioread('training_speech.wav');

%% Compute spectrogram
l_win = 512; % window length in samples
d_win = l_win/fs; % window length in ms
overlap = l_win/2; % overlap between frames

X = cf_stft(x,l_win,overlap); % stft computation
V_train = abs(X).^2; % power spectrogram
[F,N] = size(V_train);

%% Train NDS model
K_train = 100;  % dictionary size
alpha = 0.01*ones(K_train,1); % Gamma prior shape parameter
tol = 1e-5; % convergence tolerance
n_iter_max = 10000; % max nb of iterations

% make random initializations
W_ini = abs(randn(F,K_train)) + 1;
A_ini = abs(randn(K_train,K_train)) + 1;
H_ini = abs(randn(K_train,N)) + 1;

% run NDS training algorithm
[W, H, A, obj_train, Vap_train] = nds_train(V_train, alpha, W_ini, H_ini, A_ini, tol, n_iter_max);

%% Sort the columns of W greedily by similarity
% the following procedure only changes the ordering of the dictionary
% to improve display

ind     = zeros(1,K_train);
ind_res = 1:K_train;
W_res   = W;
w_prev  = ones(F,1);

for k=1:K_train
    corr = w_prev' * W_res;
    [~,jmax] = max(corr);
    w_prev = W_res(:,jmax);
    ind(k) = ind_res(jmax);
    ind_res(jmax) = [];
    W_res(:,jmax) = [];
end

W_train = W(:,ind);
H_train = H(ind,:);
A_train = A(ind,ind);

%% Display some results
figure;

subplot(4,1,1);
imagesc(log10(V_train+eps))
axis xy
title('Training data spectrogram')

subplot(4,1,2)
imagesc(log10(Vap_train+eps))
axis xy
title('Approximate spectrogram')

subplot(4,1,3)
imagesc(cf_scale(H_train')');
title('Estimated H')

subplot(4,3,10)
plot(obj_train);
title('Objective function')

subplot(4,3,11)
imagesc(log10(W_train+eps));
axis xy
title('Estimated W')

subplot(4,3,12)
imagesc(cf_scale(A_train));
axis ij
title('Estimated A')
axis square

colormap('gray')

%% Denoising example

disp('Denoising...')

% load noisy sample
[y,fs] = audioread('noisy.wav');
T = length(y);
Y = cf_stft(y,l_win,overlap);
V_obs = abs(Y).^2;
[~,N_obs] = size(V_obs);

K_noise = 2; % size of noise dictionary

% initialise noise dictionary and activations
W_noise_ini = abs(randn(F,K_noise)) + 1;
H_noise_ini = abs(randn(K_noise,N_obs)) + 1;

% initialise speech activations
H_speech_ini = abs(randn(K_train,N_obs)) + 1;

% run denoising algorithm
[W_noise, H_noise, H_speech, obj2, Vap_obs] = ...
    nds_denoise(V_obs, W_train, A_train, alpha, H_speech_ini, W_noise_ini, H_noise_ini, tol, n_iter_max);

% speech variance estimate
Vap_speech = W_train*H_speech;

%% Make results

disp('Forming results... (noisy and denoised wav files available in folder)')

S = (Vap_speech./Vap_obs).*Y; % compute speech estimate by Wiener filtering
s = cf_istft(S,l_win,overlap,T); % inverse STFT
audiowrite('denoised.wav',s,fs) % make wav file

figure;
subplot(221)
imagesc(log10(V_obs+eps))
title('Noisy spectrogram')
axis xy

subplot(222)
imagesc(log10(Vap_obs+eps))
title('Noisy approximate')
axis xy

subplot(223)
imagesc(log10(abs(S).^2+eps))
title('Spectrogram of reconstructed speech')
axis xy

subplot(224)
plot(obj2)
title('Objective function')

colormap('gray')
