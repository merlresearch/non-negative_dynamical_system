% Copyright (c) 2012,2013,2023 Mitsubishi Electric Research Laboratories (MERL).
%
% SPDX-License-Identifier: AGPL-3.0-or-later
function [W, H, A, obj, V_ap] = nds_train(V, alpha, W, H, A, tol, n_iter_max)

% MAP estimation in NDS, as in
%
% Cédric Févotte, Jonathan Le Roux, John R. Hershey,
% "Non-Negative Dynamical System With Application to Speech and Audio,"
% in Proc. IEEE International Conference on Acoustics, Speech, and Signal
% Processing (ICASSP 2013), May 2013.
%
% Inputs:
% - V: power spectrogram
% - alpha: Gamma shape parameter for H (the smaller the sparser the values
% of H)
% - W: dictionary initialisation
% - H: activation matrix initialisation
% - A: transition matrix initialisation
%
% Outputs:
% - W, H, A: estimated parameters
% - obj: values of the objective function at every iteration
% - V_ap: approximate W*H

V = V + eps;
[K,N] = size(H);

lambda_W = (N/K)*1e-6; %lambda in Eq. (9) of paper

%% Compute sufficient stats
V_ap = W*H + eps; % data approximate
G = [H(:,1) A*H(:,1:N-1)] + eps; % prediction of H
R = (H(:,2:N) + eps)./G(:,2:N);

%% Compute fit and objective values
obj = zeros(1,n_iter_max);

iter = 1;
obj(iter) = sum(V(:)./V_ap(:) + log(V_ap(:))) ...
    - sum(alpha.*sum(log(R),2)) + sum(log(H(:) + eps)) + sum(alpha.*sum(R,2)) ...
    + lambda_W * sum(W(:));

err = Inf;
fprintf('iter = %4i | obj = %+5.2E | err = %4.2E (target is %4.2E) \n',iter,obj(iter),err,tol)

while (err > tol) && (iter < n_iter_max)
    iter = iter + 1;

    %% Update H
    P = W' * (V.*V_ap.^-2);
    Q = W' * V_ap.^-1;
    S = (diag(alpha) * A)' * G.^-1;
    T = (diag(alpha) * A)' * ((H+eps)./G.^2);

    % n = 1
    n = 1;
    a = Q(:,n) + S(:,n+1) + eps;
    b = 1;
    c = H(:,n).^2 .* (P(:,n) + T(:,n+1));
    H(:,n) = 0.5 * (sqrt(b.^2+4*a.*c) - b)./a;
    G(:,n+1) = A*H(:,n) + eps;

    % 1 < n < N
    for n=2:N-1
        a = Q(:,n) + S(:,n+1) + alpha./G(:,n) + eps;
        b = 1 - alpha;
        c = H(:,n).^2 .* (P(:,n) + T(:,n+1));
        H(:,n) = 0.5 * (sqrt(b.^2+4*a.*c) - b)./a;
        G(:,n+1) = A*H(:,n) + eps;
    end

    % n = N %
    n = N;
    a = Q(:,n) + alpha./G(:,n) + eps;
    b = 1 - alpha;
    c = H(:,n).^2 .* P(:,n);
    H(:,n) = 0.5 * (sqrt(b.^2+4*a.*c) - b)./a;

    V_ap = W*H + eps;
    R = (H(:,2:N) + eps)./G(:,2:N);

    %% Update W
    W = W .* sqrt(((V.*V_ap.^-2)*H')./(V_ap.^-1*H' + lambda_W));
    V_ap = W*H + eps;

    %% Update A
    A = A .* sqrt((diag(alpha)*(R./G(:,2:N) * H(:,1:N-1)'))./(diag(alpha)*(G(:,2:N).^-1 * H(:,1:N-1)') + eps));
    G = [H(:,1), A*(H(:,1:N-1))] + eps;
    R = (H(:,2:N) + eps)./G(:,2:N);

    %% Monitor convergence
    obj(iter) = sum(V(:)./V_ap(:) + log(V_ap(:))) ...
        - sum(alpha.*sum(log(R),2)) + sum(log(H(:) + eps)) + sum(alpha.*sum(R,2)) ...
        + lambda_W * sum(W(:));

    err = abs(obj(iter)-obj(iter-1))/abs(obj(iter));

    if rem(iter,100)==0
        fprintf('iter = %4i | obj = %+5.2E | err = %4.2E (target is %4.2E) \n',iter,obj(iter),err,tol)
    end

end

%% Clean up
fprintf('iter = %4i | obj = %+5.2E | err = %4.2E (target is %4.2E) \n',iter,obj(iter),err,tol)
obj = obj(1:iter);
