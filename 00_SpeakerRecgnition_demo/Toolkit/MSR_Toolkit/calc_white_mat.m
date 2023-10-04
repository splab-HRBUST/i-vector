function W = calc_white_mat(X)
% calculates the whitening transformation for cov matrix X
[~, D, V] = svd(X);
% W = V * diag(sparse(1./( sqrt(diag(D)) + 1e-10 )));
W = V * diag((1./( sqrt(diag(D)) + 1e-10 )));