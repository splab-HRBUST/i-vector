function [patameters_PPCA] = ppca_em(M_train, ubm, num, patameters_PPCA0)

% trains a PLDA model given training data using EM algorithm. Speaker
% labels are given in spk_labs. nphi determines the dimensionality of the
% Eigenvoice subspace, and niter specifies the number of EM iteration.
%
% Technically, assuming a factor analysis (FA) model of the form:
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
%   data            : input data matrix, one observation per column
%   spk_labs        : class labels for observations in data matrix 
%   nphi            : dimensionality of the Eigenvoice subspace 
%   niter           : number of EM iterations for learning PLDA model
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

[num.SVdim, num.samples] = size(M_train);

if nargin > 3 
    patameters_PPCA = patameters_PPCA0;
    centeredM = bsxfun(@minus, M_train, patameters_PPCA.m);
    variancesM = sum(centeredM.^2, 2) / (num.samples * num.SVdim);
else
    fprintf('Randomly initializing the PPCA hyperparameters ...\n');
%     % initializing PPCA parameters with ubm
%     patameters_PPCA.m = reshape(ubm.mu, num.SVdim, 1);
%     centeredM = bsxfun(@minus, M_train, patameters_PPCA.m);
%     variancesM = sum(centeredM.^2, 2) / (num.samples * num.SVdim);
%     patameters_PPCA.I = eye(num.IVdim);
%     patameters_PPCA.T = M_train(:,1:num.IVdim)-patameters_PPCA.m;   % rand(num.SVdim, num.IVdim);
%     patameters_PPCA.Sigma = reshape(ubm.sigma, num.SVdim, 1);       % ones(num.SVdim, 1);
    
    % initializing PPCA parameters without ubm
    patameters_PPCA.m = mean(M_train,2);
    centeredM = bsxfun(@minus, M_train, patameters_PPCA.m);
    variancesM = sum(centeredM.^2, 2) / (num.samples * num.SVdim);
    patameters_PPCA.I = eye(num.IVdim);
    patameters_PPCA.T = M_train(:,1:num.IVdim)-patameters_PPCA.m;
    patameters_PPCA.Sigma = variancesM;    

end

fprintf('Baseline: Train the PPCA-based TVS with dimension %d:\n',num.IVdim);
for iter = 1 : num.nIters
    fprintf('EM iter#: %d \t', iter);
    tim1 = tic;
    
    %% Lower level for FA    
    % E step
    patameters_PPCA.B = bsxfun(@rdivide, patameters_PPCA.T', patameters_PPCA.Sigma');
    patameters_PPCA.L = patameters_PPCA.I + patameters_PPCA.B*patameters_PPCA.T;
    Ex = pinv(patameters_PPCA.L)*patameters_PPCA.B*centeredM;
    Exx = Ex*Ex' + num.samples*pinv(patameters_PPCA.L);
    
    % M step
    patameters_PPCA.T = (centeredM*Ex')*pinv(Exx);
    patameters_PPCA.Sigma = variancesM - sum(sum(Ex .* (patameters_PPCA.T' * centeredM))) / (num.SVdim*num.samples);
%     patameters_FA.Sigma = variancesM - sum(patameters_FA.T .* (patameters_FA.T * Exx), 2) / (num.SVdim*num.samples);
    
    tim1 = toc(tim1);
    fprintf('[elaps = %.2f s]\n', tim1);
end



