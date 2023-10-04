%% 
% M = m    + Tw + epsilon, epsilon~N(0,delta_M^2*I)
% Y = mu_Y + Qw + zeta,    zeta~N(0,delta_Y^2*I)
%
% EER and Min DCF are used as evaluation metrics
%
% PLDA scoring:
% EER = 2.31
% Minimum DCF = 0.013
%
% CDS:
% EER = 3.619
% Minimum DCF = 0.018

% dev: (1-150spk)*120wav
% eva: (151-200spk)*120wav
% enroll: 96/120wav
% test:   24/120wav
%% ------------------------------------------------King-ASR-010 database 
% GMM supervectors
% M_dev_all    = importdata('/data/chenchen/matlab/program/ASR_Toolbox/00_Data/King-ASR-010/Mix_60_1024/M_dev.mat');
% M_enroll = importdata('/data/chenchen/matlab/program/ASR_Toolbox/00_Data/King-ASR-010/Mix_60_1024/M_enroll.mat');
% M_test   = importdata('/data/chenchen/matlab/program/ASR_Toolbox/00_Data/King-ASR-010/Mix_60_1024/M_test.mat');

% num.spk_dev = 150;
% num.spk_eva = 50;
% 
% num.dev = 120;
% num.enroll = 96;
% num.test = 24;
% num.label_dev = repmat([1:1:150],120,1); % for PLDA
% 
num.IVdim = 800;
% 

%  
% % Speaker labels
% Y_dev_all = kron(eye(num.spk_dev),ones(num.dev,1))';



% m_dev = zscore(M_dev_all',1)';
% y_dev = zscore(Y_dev_all',1)';
% m_enroll = zscore(M_enroll',1)';
% m_test   = zscore(M_test',1)';

% m_enroll = (M_enroll-repmat(mean(M_dev_all,2),1,size(M_enroll,2)))./repmat(std(M_dev_all')',1,size(M_enroll,2));
% m_test   = (M_test-repmat(mean(M_dev_all,2),1,size(M_test,2)))./repmat(std(M_dev_all')',1,size(M_test,2));



%% ------------------------------------------------Radorm choose development data
% num.ch = 120;
% 
% fprintf('Development data: %d\n', num.ch);
% M_dev = [];
% Y_dev = [];
% for nSpk = 1 : num.spk_dev
%     M_dev(:,(nSpk-1)*num.ch+1:nSpk*num.ch) =  M_dev_all(:,(nSpk-1)*120+1:(nSpk-1)*120+num.ch);
%     Y_dev(:,(nSpk-1)*num.ch+1:nSpk*num.ch) =  Y_dev_all(:,(nSpk-1)*120+1:(nSpk-1)*120+num.ch);
% end
% clear nSpk

%% ---------------------------------------------Train the PPLS-based TVM
num.nIters = 5;

for iTers = 1% : 15
%     if iTers == 1
%         [patameters_PPLS] = ppls_em(M_dev,Y_dev,num);
%     else
%         [patameters_PPLS] = ppls_em(M_dev,Y_dev,num,patameters_PPLS);
%     end

% [patameters_PPLS] = ppls_em(M_dev_all,Y_dev_all,num);
% patameters_PPLS.m = mean(M_dev_all,2);
% patameters_PPLS.mu_y = mean(Y_dev_all,2);
% % Ex = pinv(L)*(B*centeredM+C*centeredY);

%% ---------------------------------------------Evaluation
% Ex_dev = pinv(patameters_PPLS.L)*patameters_PPLS.B*(M_dev-patameters_PPLS.m);
% Ex_enroll = pinv(patameters_PPLS.L)*patameters_PPLS.B*(M_enroll-patameters_PPLS.m);
% Ex_test   = pinv(patameters_PPLS.L)*patameters_PPLS.B*(M_test-patameters_PPLS.m);

% Cigma = patameters_PPLS.T'*patameters_PPLS.T+patameters_PPLS.deltax*2*patameters_PPLS.I;
% Y_enroll_centered = patameters_PPLS.Q*pinv(Cigma)*patameters_PPLS.T'*(M_enroll-patameters_PPLS.m);
% Y_test_centered   = patameters_PPLS.Q*pinv(Cigma)*patameters_PPLS.T'*(M_test-patameters_PPLS.m);
% beta = 0.009;
% IV_enroll = pinv(patameters_PPLS.L)*(patameters_PPLS.B*(M_enroll-patameters_PPLS.m)+beta*patameters_PPLS.C*Y_enroll_centered);
% IV_test   = pinv(patameters_PPLS.L)*(patameters_PPLS.B*(M_test-patameters_PPLS.m)+beta*patameters_PPLS.C*Y_test_centered);

% Ex_model = zeros(num.IVdim,num.spk_eva);
% for i = 1 : num.spk_eva
%     Ex_model(:,i) = mean(Ex_enroll(:,(i-1)*num.enroll+1:i*num.enroll),2); % the means of i-vectors for one speaker
% end
% clear i
% 
% answer_eva = [ones(1,50*num.test) zeros(1,50*num.test*49)];
% 
% 
% fprintf('================= cosine ====================\n');
% scores_PPLS_CDS = [];
% scores_PPLS_CDS.all =1 - pdist2(Ex_test',Ex_model','cosine');
% 
% scores_PPLS_CDS.true = [];
% scores_PPLS_CDS.impostor = [];
% 
% for a =  1 : 50
%     for b = 1 : 50
%         A = scores_PPLS_CDS.all(num.test*(a-1)+1 : num.test*(a-1)+num.test , b);
%         if a == b
%             scores_PPLS_CDS.true = [scores_PPLS_CDS.true ; A(:)];
%         elseif a~=b
%             scores_PPLS_CDS.impostor = [scores_PPLS_CDS.impostor ; A(:)];
%         end
%     end
% end  
% clear a b A
% scores_PPLS_CDS = [scores_PPLS_CDS.true;scores_PPLS_CDS.impostor];
% [eer_PPLS800_CDS(iTers),~,~,~,dcf_king_PPLS800_CDS(iTers)]=compute_eer(scores_PPLS_CDS,answer_eva,false);


fprintf('================= PLDA ======================\n');
num.Zdim = 120;
pLDA_ppls = gplda_em(Ex_dev, num.label_dev(:), num.Zdim, 100);

scores_PPLS_PLDA = [];
scores_PPLS_PLDA.all = (score_gplda_trials(pLDA_ppls, Ex_model, Ex_test))';

scores_PPLS_PLDA.true = [];
scores_PPLS_PLDA.impostor = [];

for a =  1 : 50
    for b = 1 : 50
        A = scores_PPLS_PLDA.all(num.test*(a-1)+1 : num.test*(a-1)+num.test , b);
        if a == b
            scores_PPLS_PLDA.true = [scores_PPLS_PLDA.true ; A(:)];
        elseif a~=b
            scores_PPLS_PLDA.impostor  = [scores_PPLS_PLDA.impostor  ; A(:)];
        end
    end
end
clear a b A
scores_PPLS_PLDA = [scores_PPLS_PLDA.true;scores_PPLS_PLDA.impostor];
[eer_PPLS800_PLDA200(iTers),~,~,~,dcf_king_PPLS800_PLDA200(iTers)]=compute_eer(scores_PPLS_PLDA,answer_eva,false);

fprintf('===== Summaray ======= Baseline =====%d=========\n',iTers);
fprintf('+-------------------------------------------+\n');
fprintf('|    Method     |    EER(%%)   | Min DCF_vox |\n');
fprintf('+---------------+-------------+-------------+\n');
fprintf('|PPLS(%d)+CDS     |    %2.2f    |   %2.4f    |\n', num.IVdim,         eer_PPLS800_CDS(iTers),    dcf_king_PPLS800_CDS(iTers));
fprintf('|PPLS(%d)+PLDA(%d)|    %2.2f    |   %2.4f    |\n', num.IVdim,num.Zdim,eer_PPLS800_PLDA200(iTers),dcf_king_PPLS800_PLDA200(iTers));
fprintf('+---------------+-------------+-------------+\n');

end
