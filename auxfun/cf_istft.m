% Copyright Cedric Fevotte (CNRS), 2012.
%
% SPDX-License-Identifier: GPL-3.0-or-later
function x = cf_istft(STFT,SY_WINDOW,OVERLAP,T)

% Reconstructs a signal x from a processed STFT
%
% Usage: x = cf_istft(STFT,SY_WINDOW,OVERLAP)
%
% Input:
%   - STFT: matrix of size (W/2+1 if W even) or (W/2+1/2 if W odd) x n_frames,
%   - SY_WINDOW: analysis window of size W,
%   - OVERLAP: number of samples overlap
%   - T: size of original signal
%
% Output:
%   - x is a the reconstructed signal of length OVERLAP + n_frames*(W-OVERLAP);
%
% If SY_WINDOW is an integer, a sine bell window of size SY_WINDOW
% is used by default.
%

% Default window
if  length(SY_WINDOW) == 1;
  SY_WINDOW = cf_sinebell(SY_WINDOW,OVERLAP);
end

W = length(SY_WINDOW);

% Reconstruct second half of the spectrum
if rem(W,2)==0;
%  L=W/2-1;
%  perm=zeros(L,L); perm([L:-1:1],[1:L])=eye(L);
%  STFT = [STFT; conj(perm*STFT(2:W/2,:))];
  STFT = [STFT; conj(STFT(W/2:-1:2,:))];
else
%  L=W/2-1/2;
%  perm=zeros(L,L); perm([L:-1:1],[1:l])=eye(L);
%  STFT = [STFT; conj(perm*STFT(2:W/2+1/2,:))];
  STFT = [STFT; conj(STFT((W+1)/2:-1:2))];
end

ISTFT = ifft(STFT);

x = cf_overlap_add(ISTFT,SY_WINDOW,OVERLAP);
x = real(x(OVERLAP+1:OVERLAP+T));
