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


% Speaker labels
Y_dev = full(ind2vec(num.label_dev(:)'));


m = mean(M_dev,2);
mu_y = mean(Y_dev,2);

m_dev = bsxfun(@minus, M_dev, m);
m_enroll = bsxfun(@minus, M_enroll, m);
m_test = bsxfun(@minus, M_test, m);
y_dev = bsxfun(@minus, Y_dev, mu_y);


% % % % 
% % % % m_dev = zscore(M_dev',1);
% % % % y_dev = zscore(Y_dev,1);
% % % % 
% % % % m_enroll = ((M_enroll-repmat(mean(M_dev,2),1,size(M_enroll,2)))./repmat(std(M_dev')',1,size(M_enroll,2)))';
% % % % m_test   = ((M_test-repmat(mean(M_dev,2),1,size(M_test,2)))./repmat(std(M_dev')',1,size(M_test,2)))';
% % % 
% % % % % % m_enroll = zscore(M_enroll',1);
% % % % % % m_test   = zscore(M_test',1);



%% Train PLS model

% patameters_PLS = pls_svd(m_dev',y_dev',num.spk_dev);
% patameters_PLS.m = m;
% patameters_PLS.mu_y = mu_y;

% V_X = patameters_PLS.V*(patameters_PLS.T'*patameters_PLS.V)^-1 ;


%% Evaluation

R = [25 50 75 100 125 150];%num.spk_dev;


for nR = 6% : size(R,2)

% V_X = patameters_PLS.V(:,1:R(nR))/(patameters_PLS.T(:,1:R(nR))'*patameters_PLS.V(:,1:R(nR))) ;
% 
% W_dev    = m_dev' * V_X;
% W_enroll = m_enroll' * V_X;
% W_test   = m_test' * V_X;
% 
% 
% IV_enroll = W_enroll';
% IV_model = zeros(R(nR),num.spk_eva);
% for i = 1 : num.spk_eva
%     IV_model(:,i) = mean(IV_enroll(:,(i-1)*num.enroll+1:i*num.enroll),2); % i-vector models
% end
% clear IV_enroll i
% 
% answer_eva = [ones(1,50*num.test) zeros(1,50*num.test*49)];


% % % % y_pre{nR} = m_dev'*V_X*diag(patameters_PLS.B(1:R(nR)))*patameters_PLS.C(:,1:R(nR))'; % wrong
% % % % S_RESS(nR) = norm(y_dev'-y_pre{nR},2);  % wrong

% % % % y_pre{nR} = m_dev'*V_X*diag(patameters_PLS.B(1:R(nR)))*patameters_PLS.Q(:,1:R(nR))';

% y_pre{nR} = m_dev'*patameters_PLS.V(:,1:R(nR))*diag(patameters_PLS.B(1:R(nR)))*patameters_PLS.C(:,1:R(nR))';
% Res_dev = y_dev'-y_pre{nR};
% S_RESS(nR) = sum(sum(Res_dev.*Res_dev,2));
% 
% acc(nR) = top_N_acc(y_pre{nR},num.label_dev(:),1);

% 
% fprintf('================= CDS ====================\n');
% scores_PLS_CDS = [];
% scores_PLS_CDS.all =1 - pdist2(W_test,IV_model','cosine');
% 
% scores_PLS_CDS.true = [];
% scores_PLS_CDS.impostor = [];
% for a =  1 : num.spk_eva
%     for b = 1 : num.spk_eva
%         A = scores_PLS_CDS.all(num.test*(a-1)+1 : num.test*(a-1)+num.test , b);
%         if a == b
%             scores_PLS_CDS.true = [scores_PLS_CDS.true ; A(:)];
%         elseif a~=b
%             scores_PLS_CDS.impostor = [scores_PLS_CDS.impostor ; A(:)];
%         end
%     end
% end  
% clear a b A
% 
% % % %figure;
% scores_PLS_CDS = [scores_PLS_CDS.true;scores_PLS_CDS.impostor];
% [eer_PLS_CDS(nR),~,~,~,dcf_king_PLS_CDS(nR)]=compute_eer(scores_PLS_CDS,answer_eva,0);
% % % % hold on


fprintf('================= PLDA ======================\n');
num.Zdim = 150;%R(nR);
iter_PLDA = 50
pLDA_PLS = gplda_em(W_dev', num.label_dev(:), num.Zdim, iter_PLDA);

scores_PLS_PLDA = [];
scores_PLS_PLDA.all = (score_gplda_trials(pLDA_PLS, IV_model, W_test'))';

scores_PLS_PLDA.true = [];
scores_PLS_PLDA.impostor = [];

for a =  1 : 50
    for b = 1 : 50
        A = scores_PLS_PLDA.all(num.test*(a-1)+1 : num.test*(a-1)+num.test , b);
        if a == b
            scores_PLS_PLDA.true = [scores_PLS_PLDA.true ; A(:)];
        elseif a~=b
            scores_PLS_PLDA.impostor  = [scores_PLS_PLDA.impostor  ; A(:)];
        end
    end
end
clear a b A
% 
scores_PLS_PLDA = [scores_PLS_PLDA.true;scores_PLS_PLDA.impostor];
[eer_PLS_PLDA(nR),~,~,~,dcf_king_PLS_PLDA(nR)]=compute_eer(scores_PLS_PLDA,answer_eva,0);

%% Results

fprintf('========== Summaray =========== PLS =====%d====\n',R(nR));
fprintf('+---------------------------------------------+\n');
fprintf('|    Method     |    EER(%%)   | Min DCF_king |\n');
fprintf('+---------------+-------------+---------------+\n');
fprintf('|PLS(%d)+CDS     |    %2.2f    |   %2.4f    |\n', R(nR),           eer_PLS_CDS(nR),  dcf_king_PLS_CDS(nR));
fprintf('|PLS(%d)+PLDA(%d)|    %2.2f    |   %2.4f    |\n', R(nR), num.Zdim, eer_PLS_PLDA(nR), dcf_king_PLS_PLDA(nR));
fprintf('+---------------+-------------+---------------+\n');


% fprintf('========== Summaray =========== PLS =====%d====\n',R(nR));
% fprintf('+----------------------------------------------------------+\n');
% fprintf('|    Method     |    EER(%%)   | Min DCF_king |    RESS     |\n');
% fprintf('+---------------+-------------+----------------------------+\n');
% fprintf('|PLS(%d)+CDS     |    %2.2f    |   %2.4f    |    %2.2f    |\n', num.spk_dev,           eer_PLS_CDS(nR),  dcf_king_PLS_CDS(nR),  S_RESS(nR));
% fprintf('|PLS(%d)+PLDA(%d)|    %2.2f    |   %2.4f    |    %2.2f    |\n', num.spk_dev, num.Zdim, eer_PLS_PLDA(nR), dcf_king_PLS_PLDA(nR), S_RESS(nR));
% fprintf('+---------------+-------------+----------------------------+\n');

end