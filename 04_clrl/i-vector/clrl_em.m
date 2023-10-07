function [patameters_PPLS] = clrl_em(M_dev,Y_dev,num,patameters_PPLS0)

[num.SVdim, num.samples] = size(M_dev);  
num.spk = size(Y_dev,1);  

if nargin > 3  % 输入参数大于三个（更新）
    I = patameters_PPLS0.I;
    m = patameters_PPLS0.m;
    mu_y = patameters_PPLS0.mu_y;
    T = patameters_PPLS0.T;
    Q = patameters_PPLS0.Q;
    deltax = patameters_PPLS0.deltax;
    deltay = patameters_PPLS0.deltay;
else  % 输入参数为3（第一次）-- 初始化
    % fprintf('Randomly initializing the PPLS hyperparameters ...\n');
    I = eye(num.IVdim);  % num.IVdim--公共子空间维度
    T = M_dev(:,1:num.IVdim);  %总变化空间矩阵
    Q = Y_dev(:,1:num.IVdim);  %负荷矩阵
    deltax = ones(num.SVdim,1); 
    deltay = ones(num.spk,1); 
    m = nanmean(M_dev,2);  % GMM产生过程中的偏置（行--平均）
    mu_y = nanmean(Y_dev,2);  % 类别标签产生过程中的偏置
end

centeredM = bsxfun(@minus, M_dev, m);  
varianceM = sum(centeredM.^2,2) / num.samples; 
 
centeredY = bsxfun(@minus, Y_dev, mu_y);  
varianceY = sum(centeredY.^2,2) / num.samples;

fprintf('Train the PPLS-based TVM with dimension %d:\n',num.IVdim);
for niter = 1 : num.nIters
    
    fprintf('EM iter#: %d \t', niter);
    tim1 = tic;  % 计算执行时间
    
    % E step -- 求期望
    B = bsxfun(@rdivide, T', deltax');  % T'./deltax' 
    C = bsxfun(@rdivide, Q', deltay');  % Q'./deltay' 
    L = I + B*T + C*Q;  % pinv(L)--Wn的后验协方差矩阵cov 
    Ex = pinv(L)*(B*centeredM+C*centeredY);  % Wn的后验均值矢量 
    Exx = Ex*Ex' + pinv(L)*num.samples;  % Wn的后验相关矩阵 
    
    % M step -- 最大化
    alpha = num.SVdim./sum(T.*T);  % 精度 -- 1*400 
    alpha = length_norm(alpha);
    
    x = mean(deltax);
    y = mean(deltay);
    
    T = (centeredM*Ex')*pinv(Exx + x * diag(alpha));  % 更新总变化空间矩阵 
    Q = (centeredY*Ex')*pinv(Exx);  %更新负荷矩阵 

    deltax = varianceM - sum(T .* (T*Exx),2) / num.samples;  % 协方差矩阵 
    deltay = varianceY - sum(Q .* (Q*Exx),2) / num.samples;  % 协方差矩阵 
    
    flag_x = ismember(double(0),deltax);
    if flag_x >= 1
        deltax_num = size(deltax,1);
        mmin = deltax(1,1);
        for i = 1:deltax_num
            a = deltax(i,1);
            if a~= 0
                mmin = min(mmin,a);
            end
        end
        deltax(deltax(:,1)==0,1) = mmin;
    end
       
    tim1 = toc(tim1);
    fprintf('[elaps = %.2f s]\n', tim1);
end

B = bsxfun(@rdivide, T', deltax');
C = bsxfun(@rdivide, Q', deltay');
L = I + B*T + C*Q;

patameters_PPLS.I = I;  % 单位阵
patameters_PPLS.B = B;
patameters_PPLS.C = C;
patameters_PPLS.L = L;

patameters_PPLS.m = m;  % 偏置
patameters_PPLS.mu_y = mu_y;

patameters_PPLS.T = T;  % 总变化空间矩阵
patameters_PPLS.Q = Q;  % 负荷矩阵
patameters_PPLS.deltax = deltax;  % 更新参数
patameters_PPLS.deltay = deltay;


