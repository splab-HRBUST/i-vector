%%======================================
% %% ----------------------------Train the FA-based TVS----------------------------------------

function [patameters_FA] = fa_em(M_train, ubm, num, patameters_FA0)

    [num.SVdim, num.samples] = size(M_train);  % num.SVdim-训练集的行数  num.samples-训练集的列数

    if nargin > 3   % 输入参数大于3 -- 有patameters_FA0
        patameters_FA = patameters_FA0;
        centeredM = bsxfun(@minus, M_train, patameters_FA.m);  % 求M_train-patameters_FA.m 
        variancesM = sum(centeredM.^2, 2) / num.samples;  % 对centeredM每个数开方并按列相加，再除以语音数
    else
    %     fprintf('Randomly initializing the FA hyperparameters ...\n');
        % initializing FA parameters with ubm  用ubm初始化FA参数
        patameters_FA.I = eye(num.IVdim);  % 单位矩阵 -- 400*400
        patameters_FA.Sigma = reshape(ubm.sigma, num.SVdim, 1);% variancesM;
        patameters_FA.m = reshape(ubm.mu, num.SVdim, 1);% mean(M_train,2);%reshape(ubm.mu, num.SVdim, 1)-将ubm.m按(num.SVdim,1)重新排列 
        patameters_FA.T = M_train(:,1:num.IVdim)-patameters_FA.m;   % rand(num.SVdim, num.IVdim);
        
        centeredM = bsxfun(@minus, M_train, patameters_FA.m);  % M_train与patameters_FA.m相减 
        variancesM = sum(centeredM.^2, 2) / num.samples;  % 方差 
    end

    fprintf('Baseline: Train the FA-based TVS with dimension %d:\n',num.IVdim); % 公共隐变量维度
    for iter = 1 : num.nIters  % 迭代次数
        fprintf('EM iter#: %d \t', iter);
        tim1 = tic;
        
        % E step -- 求期望
        patameters_FA.B = bsxfun(@rdivide, patameters_FA.T', patameters_FA.Sigma');  % patameters_FA.T'/patameters_FA.Sigma'
        patameters_FA.L = patameters_FA.I + patameters_FA.B*patameters_FA.T; 
        
        Ex = pinv(patameters_FA.L)*patameters_FA.B*centeredM;  
        Exx = Ex*Ex' + num.samples*pinv(patameters_FA.L);  

        % M step -- 最大化
        patameters_FA.T = (centeredM*Ex')*pinv(Exx);  
        patameters_FA.Sigma = variancesM - sum(patameters_FA.T .* (patameters_FA.T * Exx), 2) / num.samples; 
        
        % 判断patameters_FA.Sigma中是否有0
        aSigma = patameters_FA.Sigma;
        flag = ismember(0,aSigma);
        if flag >= 1  % 有0-->找到除0外的最小值
            Sigma_num = size(aSigma,1);
            mmin = aSigma(1,1);
            for i = 1:Sigma_num
                a = aSigma(i,1);
                if a ~= 0
                    mmin = min(mmin,a);
                end
            end
            % 将0替换成mmin
            aSigma(aSigma(:,1)==0,1) = mmin;
            patameters_FA.Sigma = aSigma;
        end
        
        tim1 = toc(tim1);
        fprintf('[elaps = %.2f s]\n', tim1);
    end



