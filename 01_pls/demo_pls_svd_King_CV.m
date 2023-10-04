%% Demo for PLS-based TVS
%
% dev: (1-150spk)*120wav
% eva: (151-200spk)*120wav
% enroll: 96/120wav
% test:   24/120wav
%
% 
% ========== Summaray =========== PLS =============
% +-----------------------------------------------+
% |      Method      |    EER(%)   | Min DCF_king |
% +------------------+-------------+--------------+
% |PLS(150)+CDS      |    3.08     |    0.5950    |
% |PLS(150)+PLDA(150)|    2.47     |    0.4333    |
% +------------------+-------------+--------------+

%% Load data
% % GMM supervectors
% M_dev    = importdata('/data/chenchen/matlab/program/ASR_Toolbox/00_Data/King-ASR-010/Mix_60_1024/M_dev.mat');
% M_enroll = importdata('/data/chenchen/matlab/program/ASR_Toolbox/00_Data/King-ASR-010/Mix_60_1024/M_enroll.mat');
% M_test   = importdata('/data/chenchen/matlab/program/ASR_Toolbox/00_Data/King-ASR-010/Mix_60_1024/M_test.mat');

% Information 
num.spk_dev = 150;
num.spk_eva = 50;
num.dev = 120;
num.enroll = 96;
num.test = 24;
num.label_dev = repmat([1:1:150],120,1); % for PLDA
num.Zdim = 150;


% R = [25 50 75 100 125 150];
% 
% [M_CV, label_CV] = cv_king(M_dev);
% Y_dev_CV = full(ind2vec(label_CV));

for nR = 1 : size(R,2)

    % Speaker labels
    Y_dev_CVtrain = full(ind2vec(label_CV(1:3600*4)));
    Y_dev_CVtest  = full(ind2vec(label_CV(3600*4+1:end)));
    

    for nCV = 1 : 5

        ind = [1:nCV-1,nCV+1:5];
        M_dev_CVtrain = cell2mat(M_CV(ind));
        M_dev_CVtest  = cell2mat(M_CV(nCV));

        m = mean(M_dev_CVtrain,2);
        mu_y = mean(Y_dev_CV,2);

        m_dev_CVtrain = bsxfun(@minus, M_dev_CVtrain, m);
        m_dev_CVtest = bsxfun(@minus, M_dev_CVtest, m);
        y_dev_CVtrain = bsxfun(@minus, Y_dev_CVtrain, mu_y);
        y_dev_CVtest  = bsxfun(@minus, Y_dev_CVtest, mu_y);

% %% Train PLS model
% 
%         patameters_PLS = pls_svd(m_dev_CVtrain',y_dev_CVtrain',R(nR));
%         patameters_PLS.m = m;
%         patameters_PLS.mu_y = mu_y;
%         
%         patameters_PLS_CV{nR,nCV} = patameters_PLS;
          patameters_PLS = patameters_PLS_CV{nR,nCV};
% 
% %% Evaluation
% 
        V_X = patameters_PLS.V(:,1:R(nR))/(patameters_PLS.T(:,1:R(nR))'*patameters_PLS.V(:,1:R(nR)));
%         y_pre_CVtest = m_dev_CVtest'*V_X*diag(patameters_PLS.B(1:R(nR)))*patameters_PLS.C(:,1:R(nR))'; 

        y_pre_CVtest = m_dev_CVtest'*V_X*diag(patameters_PLS.B(1:R(nR)))*patameters_PLS.Q(:,1:R(nR))';
        Y_pre_CVtest{nCV,nR} = bsxfun(@plus, y_pre_CVtest, patameters_PLS.mu_y');
% 
% 
    end
    
%% PRESS
    
% % %     S_PRESS(nR) = norm(Y_dev_CV'-cell2mat(Y_pre_CVtest(:,nR)),2);   % wrong

    Res_dev_CV = Y_dev_CV'-cell2mat(Y_pre_CVtest(:,nR));
    S_PRESS(nR) = sum(sum(Res_dev_CV.*Res_dev_CV,2));

%     Y_pre_nR = cell2mat(Y_pre_CVtest(:,nR));
%     acc_CV(nR) = top_N_acc(Y_pre_nR,label_CV,1);
    
end