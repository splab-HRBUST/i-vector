function [patameters_PPLS] = ppls_n_em(M_dev,Y_dev,num,patameters_PPLS0)

[num.SVdim, num.samples] = size(M_dev);
num.spk = size(Y_dev,1);

if nargin > 3 
    I = patameters_PPLS0.I;
    m = patameters_PPLS0.m;
    mu_y = patameters_PPLS0.mu_y;
    T = patameters_PPLS0.T;
    Q = patameters_PPLS0.Q;
    Sigmax = patameters_PPLS0.Sigmax;
    Sigmay = patameters_PPLS0.Sigmay;
else
    % fprintf('Randomly initializing the PPLS hyperparameters ...\n');

    I = eye(num.IVdim);
    T = M_dev(:,1:num.IVdim);
    Q = Y_dev(:,1:num.IVdim);
    Sigmax = ones(num.SVdim, 1);
    Sigmay = ones(num.spk, 1);
    m = mean(M_dev,2);
    mu_y = mean(Y_dev,2);
    
%     mu_y = 0;
end


centeredM = bsxfun(@minus, M_dev, m);
variancesM = sum(centeredM.^2,2) / (num.samples);
 
centeredY = bsxfun(@minus, Y_dev, mu_y);
variancesY = sum(centeredY.^2,2) / (num.samples);

% centeredM = M_dev;
% varianceM = 1;

% centeredY = Y_dev;
% varianceY = sum(sum(Y_dev.^2)) / (num.samples * num.spk);


fprintf('Train the PPLS-N-based TVM with dimension %d:\n',num.IVdim);
for niter = 1 : num.nIters
    
    fprintf('EM iter#: %d \t', niter);
    tim1 = tic;
    
    % E step 
    B = bsxfun(@rdivide, T', Sigmax');
    C = bsxfun(@rdivide, Q', Sigmay');
    L = I + B*T + C*Q;
    Ex = pinv(L)*(B*centeredM+C*centeredY);
    Exx = Ex*Ex' + pinv(L)*num.samples;
    
    % M step   
    T = (centeredM*Ex')*pinv(Exx);
    Q = (centeredY*Ex')*pinv(Exx);    
    Sigmax = variancesM - sum(T .* (T * Exx),2) / num.samples;
    Sigmay = variancesY - sum(Q .* (Q * Exx),2) / num.samples;
    
%     Sigmax = variancesM - sum(Ex .* (T' * centeredM),2) / num.samples;
%     Sigmay = variancesY - sum(Ex .* (Q' * centeredY),2) / num.samples;

    
    tim1 = toc(tim1);
    fprintf('[elaps = %.2f s]\n', tim1);
end

B = bsxfun(@rdivide, T', Sigmax');
C = bsxfun(@rdivide, Q', Sigmay');
L = I + B*T + C*Q;

patameters_PPLS.I = I;
patameters_PPLS.B = B;
patameters_PPLS.C = C;
patameters_PPLS.L = L;

patameters_PPLS.m = m;
patameters_PPLS.mu_y = mu_y;

patameters_PPLS.T = T;
patameters_PPLS.Q = Q;
patameters_PPLS.Sigmax = Sigmax;
patameters_PPLS.Sigmay = Sigmay;