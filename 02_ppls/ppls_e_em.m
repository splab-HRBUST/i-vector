function [patameters_PPLS] = ppls_e_em(M_dev,Y_dev,num,patameters_PPLS0)

[num.SVdim, num.samples] = size(M_dev);
num.spk = size(Y_dev,1);

if nargin > 3 
    I = patameters_PPLS0.I;
    m = patameters_PPLS0.m;
    mu_y = patameters_PPLS0.mu_y;
    T = patameters_PPLS0.T;
    Q = patameters_PPLS0.Q;
    deltax = patameters_PPLS0.deltax;
    deltay = patameters_PPLS0.deltay;
else
    % fprintf('Randomly initializing the PPLS hyperparameters ...\n');

    I = eye(num.IVdim);
    T = M_dev(:,1:num.IVdim);
    Q = Y_dev(:,1:num.IVdim);
    deltax = 1;
    deltay = 1;
    m = mean(M_dev,2);
    mu_y = mean(Y_dev,2);
%     mu_y = 0;
end


centeredM = bsxfun(@minus, M_dev, m);
varianceM = sum(sum(centeredM.^2)) / (num.samples * num.SVdim);
 
centeredY = bsxfun(@minus, Y_dev, mu_y);
varianceY = sum(sum(centeredM.^2)) / (num.samples * num.spk);

% centeredM = M_dev;
% varianceM = 1;

% centeredY = Y_dev;
% varianceY = sum(sum(Y_dev.^2)) / (num.samples * num.spk);


fprintf('Train the PPLS-based TVM with dimension %d:\n',num.IVdim);
for niter = 1 : num.nIters
    
    fprintf('EM iter#: %d \t', niter);
    tim1 = tic;
    
    % E step 
    B = bsxfun(@rdivide, T', deltax');
    C = bsxfun(@rdivide, Q', deltay');
    L = I + B*T + C*Q;
    Ex = pinv(L)*(B*centeredM+C*centeredY);
    Exx = Ex*Ex' + pinv(L)*num.samples;
    
    % M step   
    T = (centeredM*Ex')*pinv(Exx);
    Q = (centeredY*Ex')*pinv(Exx);
    deltax = varianceM - sum(sum(Ex .* (T' * centeredM))) / (num.SVdim*num.samples);
    deltay = varianceY - sum(sum(Ex .* (Q' * centeredY))) / (num.spk*num.samples);

    
    tim1 = toc(tim1);
    fprintf('[elaps = %.2f s]\n', tim1);
end

B = bsxfun(@rdivide, T', deltax');
C = bsxfun(@rdivide, Q', deltay');
L = I + B*T + C*Q;

patameters_PPLS.I = I;
patameters_PPLS.B = B;
patameters_PPLS.C = C;
patameters_PPLS.L = L;

patameters_PPLS.m = m;
patameters_PPLS.mu_y = mu_y;

patameters_PPLS.T = T;
patameters_PPLS.Q = Q;
patameters_PPLS.deltax = deltax;
patameters_PPLS.deltay = deltay;