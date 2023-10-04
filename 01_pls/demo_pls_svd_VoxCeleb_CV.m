%%
% Fac = 400
% eer_PLS_CDS = 10.79
% dcf_vox_PLS_CDS = 0.76

% Fac = 1211
% eer_PLS_CDS = 11.1188
% dcf_vox_PLS_CDS = 0.7787

%% Loading data

% M_dev = importdata('/data/chenchen/matlab/program/ASR_Toolbox/00_Data/VoxCeleb1/Mix_60_1024/M_dev.mat');
% M_eva = importdata('/data/chenchen/matlab/program/ASR_Toolbox/00_Data/VoxCeleb1/Mix_60_1024/M_eva.mat');

% label_dev = importdata('/data/chenchen/matlab/program/ASR_Toolbox/00_Data/VoxCeleb1/Mix_60_1024/label_dev.mat');
% Y_dev = full(ind2vec(label_dev));
% 
% Fac = 1211;
% 
% 
% m = mean(M_dev,2);
% mu_y = mean(Y_dev,2);
% 
% m_dev = bsxfun(@minus, M_dev, m);
% m_eva = bsxfun(@minus, M_eva, m);
% 
% y_dev = bsxfun(@minus, Y_dev, mu_y);

%% train PLS model

% patameters_PLS = pls_svd(m_dev',y_dev',Fac);
% patameters_PLS.m = m;
% patameters_PLS.mu_y = mu_y;



% W_X = W*(P'*W)^-1 ;


%% evaluation


% [M_CV, label_CV] = cv_vox(M_dev,spk_cell_dev,label_dev);

% M_CV = importdata('M_CV_vox.mat');
% label_CV = importdata('label_CV_vox.mat');
% S_Vox = importdata('S_Vox.mat');
acc_C_Vox = S_Vox.acc_C_Vox;
C_tot = S_Vox.C_tot;

R = [100 : 100: 1200];%num.spk_dev;

label_dev = cell2mat(label_CV);
Y_dev = full(ind2vec(label_dev));

for nR = 7 : 11%2 : size(R,2)
    
    
    for nCV = 1 : 5

        ind = [1:nCV-1,nCV+1:5];
        M_dev_CVtrain = cell2mat(M_CV(ind));
        M_dev_CVtest  = cell2mat(M_CV(nCV));
        
        label_dev_CVtrain = cell2mat(label_CV(ind));
        label_dev_CVtest = cell2mat(label_CV(nCV));
        
        Y_dev_CVtrain = full(ind2vec(label_dev_CVtrain));
        Y_dev_CVtest = full(ind2vec(label_dev_CVtest));

        m = mean(M_dev_CVtrain,2);
        mu_y = mean(Y_dev_CVtrain,2);

        M_dev_CVtrain = bsxfun(@minus, M_dev_CVtrain, m);
        M_dev_CVtest = bsxfun(@minus, M_dev_CVtest, m);
        y_dev_CVtrain = bsxfun(@minus, Y_dev_CVtrain, mu_y);
        y_dev_CVtest  = bsxfun(@minus, Y_dev_CVtest, mu_y);

%% Train PLS model

%         patameters_PLS = pls_svd(M_dev_CVtrain',y_dev_CVtrain',R(nR));
%         patameters_PLS.m = m;
%         patameters_PLS.mu_y = mu_y;
%         
%         patameters_PLS_CV{nR,nCV} = patameters_PLS;
          patameters_PLS = patameters_PLS_CV{12,nCV};
% 
%% Evaluation
% 
%         V_X = patameters_PLS.V(:,1:R(nR))/(patameters_PLS.T(:,1:R(nR))'*patameters_PLS.V(:,1:R(nR)));
%         y_pre_CVtest = m_dev_CVtest'*V_X*diag(patameters_PLS.B(1:R(nR)))*patameters_PLS.C(:,1:R(nR))'; 


        y_pre_CVtest = M_dev_CVtest'*patameters_PLS.V(:,1:R(nR))*diag(patameters_PLS.B(1:R(nR)))*patameters_PLS.C(:,1:R(nR))';
        Y_pre_CVtest{nCV,nR} = bsxfun(@plus, y_pre_CVtest, patameters_PLS.mu_y');
% 
% 
    end
    
%% PRESS and ACC_pre
    

    Res_dev_CV = Y_dev'-cell2mat(Y_pre_CVtest(:,nR));    
    acc{nR} = top_N_acc(cell2mat(Y_pre_CVtest(:,nR)),label_dev,1);
    
    C_PRESS(nR) = sum(sum(Res_dev_CV.*Res_dev_CV,2));
    CD(nR) = 1-C_PRESS(nR)/C_tot;
    
    C_PLACCR(nR) = 1 - acc{nR}/acc_C_Vox{nR}(1,1);
    
    
    fprintf('===== Summaray PLS ==========%d=======\n',R(nR));
    fprintf('+-------------------------------------+\n');
    fprintf('|     PRESS    |     CD     |   PLACCR   |\n');
    fprintf('+-------------------------+\n');
    fprintf('|    %2.2f    |   %2.4f    |   %2.4f    |\n', C_PRESS(nR), CD(nR), C_PLACCR(nR));
    fprintf('+-------------------------------------+\n');

end