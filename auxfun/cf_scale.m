% Copyright Cedric Fevotte (CNRS), 2012.
%
% SPDX-License-Identifier: GPL-3.0-or-later
function A = cf_scale(A)

% Normalizes the columns of matrix A so that they sum to 1

A = repmat(sum(A,1).^-1,size(A,1),1) .* A;
