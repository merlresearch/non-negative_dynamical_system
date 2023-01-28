% Copyright (c) 2012,2013,2023 Mitsubishi Electric Research Laboratories (MERL).
%
% SPDX-License-Identifier: AGPL-3.0-or-later
function [W_no, H_no, H_sp, obj, V_ap] = nds_denoise(V, W_sp, A_sp, alpha_sp, H_sp, W_no, H_no, tol, n_iter_max)

% Speech enhancement with NDS
% Implements section 4.2 of
% Cédric Févotte, Jonathan Le Roux, John R. Hershey,
% "Non-Negative Dynamical System With Application to Speech and Audio,"
% in Proc. IEEE International Conference on Acoustics, Speech, and Signal
% Processing (ICASSP 2013), May 2013.

V = V + eps;
[F,N] = size(V);
W = [W_sp, W_no];
H = [H_sp; H_no];

%% Compute sufficient stats
V_ap = W*H + eps; % data approximate
G = [H_sp(:,1) A_sp*H_sp(:,1:N-1)] + eps; % prediction of H
R = (H_sp(:,2:N) + eps)./G(:,2:N);

%% Compute fit and objective values
obj = zeros(1,n_iter_max);

iter = 1;
obj(iter) = sum(V(:)./V_ap(:) + log(V_ap(:))) ...
    - sum(alpha_sp.*sum(log(R),2)) + sum(log(H_sp(:) + eps)) ...
    + sum(alpha_sp.*sum(R,2));

err = Inf;

fprintf('iter = %4i | obj = %+5.2E | err = %4.2E (target is %4.2E) \n',iter,obj(iter),err,tol)

while (err > tol) && (iter < n_iter_max)

    iter = iter + 1;

    %% Update W_no
    W_no = W_no .* ((V.*V_ap.^-2)*H_no')./(V_ap.^-1*H_no');
    W = [W_sp, W_no];
    V_ap = W*H + eps;

    %% Update H_sp
    P  = W_sp' * (V.*V_ap.^-2);
    Q  = W_sp' * V_ap.^-1;
    S = (diag(alpha_sp) * A_sp)' * G.^-1;
    T = (diag(alpha_sp) * A_sp)' * ((H_sp+eps)./G.^2);

    % n = 1
    n = 1;
    a = Q(:,n) + S(:,n+1) + eps;
    b = 1;
    c = H_sp(:,n).^2 .* (P(:,n) + T(:,n+1));
    H_sp(:,n) = 0.5 * (sqrt(b.^2+4*a.*c) - b)./a;
    G(:,n+1) = A_sp*H_sp(:,n) + eps;

    % 1 < n < N
    for n=2:N-1
        a = Q(:,n) + S(:,n+1) + alpha_sp./G(:,n) + eps;
        b = 1-alpha_sp;
        c = H_sp(:,n).^2 .* (P(:,n) + T(:,n+1));
        H_sp(:,n) = 0.5 * (sqrt(b.^2+4*a.*c) - b)./a;
        G(:,n+1) = A_sp*H_sp(:,n) + eps;
    end

    % n = N %
    n = N;
    a = Q(:,n) + alpha_sp./G(:,n) + eps;
    b = 1-alpha_sp;
    c = H_sp(:,n).^2 .* P(:,n);
    H_sp(:,n) = 0.5 * (sqrt(b.^2+4*a.*c) - b)./a;

    R = (H_sp(:,2:N) + eps)./G(:,2:N);
    H = [H_sp; H_no];
    V_ap = W*H + eps;

    %% Update H_no
    H_no = H_no .* ((W_no'*(V.*V_ap.^-2))./(W_no'*V_ap.^-1));
    H = [H_sp; H_no];
    V_ap = W*H + eps;

    %% Scale W_no and H_no
    scale = sum(W_no,1);
    W_no = W_no .* repmat(scale.^-1,F,1);
    H_no = H_no .* repmat(scale',1,N);

    %% Monitor convergence %%
    obj(iter) = sum(V(:)./V_ap(:) + log(V_ap(:))) ...
        - sum(alpha_sp.*sum(log(R),2)) + sum(log(H_sp(:) + eps)) ...
        + sum(alpha_sp.*sum(R,2));

    err = abs(obj(iter)-obj(iter-1))/abs(obj(iter));

    if rem(iter,100)==0
        fprintf('iter = %4i | obj = %+5.2E | err = %4.2E (target is %4.2E) \n',iter,obj(iter),err,tol)
    end

end

%% Clean up
fprintf('iter = %4i | obj = %+5.2E | err = %4.2E (target is %4.2E) \n',iter,obj(iter),err,tol)
obj = obj(1:iter);
