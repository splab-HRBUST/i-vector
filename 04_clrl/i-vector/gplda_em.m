function plda = gplda_em(data, spk_labs, nphi, niter, plda0)
% trains a PLDA model given training data using EM algorithm. Speaker
% labels are given in spk_labs. nphi determines the dimensionality of the
% Eigenvoice subspace, and niter specifies the number of EM iteration.
% 使用EM算法训练给定训练数据的PLDA模型。扬声器标签由spk_实验室提供。nphi确定特征语音子空间的维数，niter指定EM迭代的次数

% Technically, assuming a factor analysis (FA) model of the form:  从技术上讲，假设采用以下形式的因子分析（FA）模型
%
%           x = m + Phi . y + e
%
% for i-vectors, the code computes the maximum likelihood estimate (MLE)of
% the factor loading matrix Phi (aka the Eigenvoice subspace). Here, x is
% the i-vector, m is the mean of training i-vectors, and y~N(0,I) is a
% vector of latent factors. The residual noise term e explains the
% variability not captured through the latent variables.
%
% Inputs:
%   data            : input data matrix, one observation per column  输入数据矩阵，每列观察一次
%   spk_labs        : class labels for observations in data matrix  数据矩阵中观测值的类标签
%   nphi            : dimensionality of the Eigenvoice subspace  特征语音子空间的维数
%   niter           : number of EM iterations for learning PLDA model  学习PLDA模型的EM迭代次数
%
% Outputs:
%    plda           : a structure containing the PLDA hyperparameters
%					  as well as the mean of development data and a whitening 
%					  transform:(plda.Phi: Eigenvoice matrix, plda.Sigma: covariance
%					  matrix of the residual noise, plda.M: mean, plda.W: whitening transform)
%
% References:
%   [1] S.J.D. Prince and J.H. Elder, "Probabilistic linear discriminant analysis
%       for inferences about identity," in Proc. IEEE ICCV, Rio de Janeiro, Brazil,
%       Oct. 2007.
%   [2] D. Garcia-Romero and C.Y. Espy-Wilson, "Analysis of i-vector length 
%       normalization in speaker recognition systems," in Proc. INTERSPEECH,
%       Florence, Italy, Aug. 2011, pp. 249-252.
%   [3] P. Kenny, "Bayesian speaker verification with heavy-tailed priors," 
%       in Proc. Odyssey, The Speaker and Language Recognition Workshop, Brno, 
%       Czech Republic, Jun. 2010.
%
%
% Omid Sadjadi <s.omid.sadjadi@gmail.com>
% Microsoft Research, Conversational Systems Research Center

[ndim, nobs] = size(data);  %(800*1800)

if ( nobs ~= length(spk_labs) )
	error('oh dear! number of data samples should match the number of labels!');
end

% make sure the labels are sorted
[spk_labs, I] = sort(spk_labs);
data = data(:, I);
[~, ia, ic] = unique(spk_labs, 'stable');  % 去重(10)  ia--第一个不重复数的位置  ic--每个数属于第几个不重复数（标签）
spk_counts = histc(ic, 1 : numel(ia)); % # sessions per speaker  每个发言者的会议数

M = mean(data, 2);
data = bsxfun(@minus, data, M); % centering the data
data = length_norm(data); % normalizing the length
W1   = calc_white_mat(cov(data'));
% W1   = 1;%calc_white_mat(cov(data'));
data = W1' * data; % whitening the data

if nargin > 4
    Sigma = plda0.Sigma;
    Phi = plda0.Phi;
else
    fprintf('Randomly initializing the PLDA hyperparameters ...\n');
    % Initialize the parameters randomly
    [s1, s2] = RandStream.create('mrg32k3a', 'NumStreams', 2);
    Sigma    = 100 * randn(s1, ndim); % covariance matrix of the residual term
    Phi = randn(s2, ndim, nphi); % factor loading matrix (Eignevoice matrix)
    Phi = bsxfun(@minus, Phi, mean(Phi, 2));
    W2   = calc_white_mat(Phi' * Phi);
    Phi = Phi * W2; % orthogonalize Eigenvoices (columns)
end

fprintf('Re-estimating the Eigenvoice subspace with %d factors ...\n', nphi);
for iter = 1 : niter
%     fprintf('EM iter#: %d \t', iter);
    tim = tic;
    % expectation
    [Ey, Eyy] = expectation_plda(data, Phi, Sigma, spk_counts);
    % maximization
    [Phi, Sigma] = maximization_plda(data, Ey, Eyy);
    tim = toc(tim);
%     fprintf('[elaps = %.2f s]\n', tim);
end

plda.Phi   = Phi;
plda.Sigma = Sigma;
plda.W     = W1;
plda.M     = M;

function [Ey, Eyy] = expectation_plda(data, Phi, Sigma, spk_counts)
% computes the posterior mean and covariance of the factors
nphi     = size(Phi, 2);
nsamples = size(data, 2);
nspks    = size(spk_counts, 1);

Ey  = zeros(nphi, nsamples);
Eyy = zeros(nphi);

% initialize common terms to save computations
uniqFreqs  	  = unique(spk_counts);
nuniq 		  = size(uniqFreqs, 1);
invTerms      = cell(nuniq, 1);
invTerms(:)   = {zeros(nphi)};
PhiT_invS_Phi = ( Phi'/Sigma ) * Phi;
I = eye(nphi);
for ix = 1 : nuniq
    nPhiT_invS_Phi = uniqFreqs(ix) * PhiT_invS_Phi;
    Cyy =  pinv(I + nPhiT_invS_Phi);
    invTerms{ix} = Cyy;
end

% data = Sigma\data;
data = pinv(Sigma)*data;
cnt  = 1;
for spk = 1 : nspks
    nsessions = spk_counts(spk);
    % Speaker indices
    idx = cnt : ( cnt - 1 ) + spk_counts(spk);
    cnt  = cnt + spk_counts(spk);
    Data = data(:, idx);
    PhiT_invS_y = sum(Phi' * Data, 2);
    Cyy = invTerms{ uniqFreqs == nsessions };
    Ey_spk  = Cyy * PhiT_invS_y;
    Eyy_spk = Cyy + Ey_spk * Ey_spk';
    Eyy     = Eyy + nsessions * Eyy_spk;
    Ey(:, idx) = repmat(Ey_spk, 1, nsessions);
end

function [Phi, Sigma] = maximization_plda(data, Ey, Eyy)
% ML re-estimation of the Eignevoice subspace and the covariance of the
% residual noise (full).
nsamples = size(data, 2);
Data_sqr = data * data';
Phi      = data * Ey' / (Eyy);
Sigma    = 1/nsamples * (Data_sqr - (Phi * Ey) * data');

function W = calc_white_mat(X)
% calculates the whitening transformation for cov matrix X
[~, D, V] = svd(X);
% W = V * diag(sparse(1./( sqrt(diag(D)) + 1e-10 )));
W = V * diag((1./( sqrt(diag(D)) + 1e-10 )));
