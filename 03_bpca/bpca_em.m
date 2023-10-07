function [patameters_BPCA] = bpca_em(M_dev,num,patameters_BPCA0)

[num.SVdim, num.samples] = size(M_dev);

if nargin > 2 
    I = patameters_BPCA0.I;
    m = patameters_BPCA0.m;
    T = patameters_BPCA0.T;
    deltax = patameters_BPCA0.deltax;
else
    fprintf('Randomly initializing the BPCA hyperparameters ...\n');
    I = eye(num.IVdim);
    T = M_dev(:,1:num.IVdim);
    deltax = 1;
    m = mean(M_dev,2);
end

centeredM = bsxfun(@minus, M_dev, m);
varianceM = sum(sum(centeredM.^2)) / (num.samples * num.SVdim);

fprintf('Train the BPCA-based TVM with dimension %d:\n',num.IVdim);
for niter = 1 : num.nIters
    
    fprintf('EM iter#: %d \t', niter);
    tim1 = tic;
    
    % E step 
    B = bsxfun(@rdivide, T', deltax');
    L = I + B*T;  
    Ex = pinv(L)*(B*centeredM);
    Exx = Ex*Ex' + pinv(L)*num.samples;    
    % M step
    alpha = num.SVdim./sum(T.*T);
%     tmp = ;
%     tmp(isinf(tmp))=0;
    T = (centeredM*Ex')*pinv(Exx + deltax * diag(alpha));
    %错误代码：
    %deltax = varianceM + (sum(sum(sum(Ex .* (T' * centeredM))) - 2*sum(T .* (T * Exx),2)))/ (num.SVdim*num.samples);
    %改正后：
    deltax = varianceM + (sum(sum(T .* (T * Exx),2)) - 2*sum(sum(Ex .* (T' * centeredM))))/ (num.SVdim*num.samples);
%     deltax = varianceM - sum(sum(Ex .* (T' * centeredM)))/ (num.SVdim*num.samples);
    tim1 = toc(tim1);
    fprintf('[elaps = %.2f s]\n', tim1);
end

B = bsxfun(@rdivide, T', deltax');
L = I + B*T;

patameters_BPCA.I = I;
patameters_BPCA.B = B;
patameters_BPCA.L = L;
patameters_BPCA.m = m;
patameters_BPCA.T = T;
patameters_BPCA.deltax = deltax;
