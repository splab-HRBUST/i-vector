%% Test list
% dev: (1-150spk)*120wav
% eva: (151-200spk)*120wav
% enroll: 96/120wav
% test:   24/120wav

% ===== Summaray ======= Baseline ====== FEFA =======
% +------------------------------------------------+
% |    Method         |    EER(%)   | Min DCF_king |
% +-------------------+-------------+--------------+
% |FEFA(400)+CDS      |    2.51     |    1.0367    |
% |FEFA(400)+PLDA(100)|    6.08     |    1.2183    |
% +-------------------+-------------+--------------+

%% ------------------------------------------------Load King-ASR-010 dataset
% fprintf('loading data\n');
% 
% AllFiles = importdata('/data/chenchen/matlab/program/ASR_Toolbox/00_Data/King-ASR-010/AllFiles.ndx');
% 
% for nFile = 1 : size(AllFiles,1)
%     mfccPath = fullfile('/data/chenchen/matlab/program/ASR_Toolbox/00_Data/King-ASR-010/Mix_60_1024/mfcc',char(AllFiles{nFile}));
%     mfcc_all_cell(nFile) = importdata(mfccPath);
% end
% % save('/data/chenchen/matlab/program/ASR_Toolbox/00_Data/King-ASR-010/Mix_60_1024/mfcc_all_cell.mat','mfcc_all_cell','-v7.3'); 
% 
% mfcc_dev_cell = mfcc_all_cell(1:150*120);
% mfcc_eva_cell = mfcc_all_cell(150*120+1:end);
% clear nFile mfccPath
% 
% % save('/data/chenchen/matlab/program/ASR_Toolbox/00_Data/King-ASR-010/Mix_60_1024/mfcc_dev_cell.mat','mfcc_dev_cell','-v7.3'); 
% % save('/data/chenchen/matlab/program/ASR_Toolbox/00_Data/King-ASR-010/Mix_60_1024/mfcc_eva_cell.mat','mfcc_eva_cell','-v7.3'); 
% 
% num.enroll = 96;
% num.test = 24;
% for nSpk = 1 : 50
%     mfcc_enroll_cell((nSpk-1)*num.enroll+1:nSpk*num.enroll) =  mfcc_eva_cell((nSpk-1)*120+1:(nSpk-1)*120+96);
%     mfcc_test_cell((nSpk-1)*num.test+1  :nSpk*num.test)   =  mfcc_eva_cell((nSpk-1)*120+97:nSpk*120);
% end
% clear nSpk
% % save('/data/chenchen/matlab/program/ASR_Toolbox/00_Data/King-ASR-010/Mix_60_1024/mfcc_enroll_cell.mat','mfcc_enroll_cell','-v7.3'); 
% % save('/data/chenchen/matlab/program/ASR_Toolbox/00_Data/King-ASR-010/Mix_60_1024/mfcc_test_cell.mat','mfcc_test_cell','-v7.3'); 

% mfcc_dev_cell    = importdata('/data/chenchen/matlab/program/ASR_Toolbox/00_Data/King-ASR-010/Mix_60_1024/mfcc_dev_cell.mat');
% mfcc_enroll_cell = importdata('/data/chenchen/matlab/program/ASR_Toolbox/00_Data/King-ASR-010/Mix_60_1024/mfcc_enroll_cell.mat');
% mfcc_test_cell   = importdata('/data/chenchen/matlab/program/ASR_Toolbox/00_Data/King-ASR-010/Mix_60_1024/mfcc_test_cell.mat');

%% Train UBM
% tic; 
% nmix = 1024;
% final_niter = 10;
% ds_factor = 1;
% nworkers = 8;
% ubm = gmm_em(mfcc_dev_cell(:), nmix, final_niter, ds_factor, nworkers);
% fprintf('UBM trained successfully\n');
% toc; 
% clear nmix final_niter ds_factor nworkers
% % save('/data/chenchen/matlab/program/ASR_Toolbox/00_Data/King-ASR-010/Mix_60_1024/model/ubm.mat','ubm','-v7.3');

% ubm = importdata('/data/chenchen/matlab/program/ASR_Toolbox/00_Data/King-ASR-010/Mix_60_1024/model/ubm.mat');

%% Caculate the BW statistics
% tic;
% stats_dev = {};
% for nFiles=1:size(mfcc_dev_cell,2)
%         [N_dev,F_dev] = compute_bw_stats(mfcc_dev_cell{1,nFiles}, ubm);
%         stats_dev{nFiles} = [N_dev; F_dev];
% end
% clear nFiles N_dev F_dev
% toc;

%% Train the FEFA-based TVM
% tic;

% num.IVdim = 400; % dimension of i-vector
num.Zdim = 100;  % dimension of latent vector for PLDA
%  
% num.label_dev = repmat([1:1:150],120,1); % for PLDA
% num.enroll = 96;
% num.test = 24;
% num.nIters = 5;
% 
% nworkers = 4;
% fprintf('Train the FEFA-based TVM: %d dimension\n',num.IVdim);
% T = train_tv_space(stats_dev(:), ubm, num.IVdim, num.nIters, nworkers);
% toc;
        
%% Evaluation
% Ex_dev = [];
% for nFiles=1:size(mfcc_dev_cell,2)%size(restmfcc,1)
%         Ex_dev(:,nFiles) = extract_ivector(stats_dev{nFiles}, ubm, T);
% end
% clear nFiles
% 
% Ex_enroll = [];
% for nFiles=1:size(mfcc_enroll_cell,2)%size(restmfcc,1)
%         [N_enroll,F_enroll] = compute_bw_stats(mfcc_enroll_cell{1,nFiles}, ubm);
%         Ex_enroll(:,nFiles) = extract_ivector([N_enroll; F_enroll], ubm, T);
% end
% clear nFiles N_enroll F_enroll
% 
% Ex_test = [];
% for nFiles=1:size(mfcc_test_cell,2)%size(restmfcc,1)
%         [N_test,F_test] = compute_bw_stats(mfcc_test_cell{1,nFiles}, ubm);
%         Ex_test(:,nFiles) = extract_ivector([N_test; F_test], ubm, T);
% end
% clear nFiles N_test F_test
% 
% Ex_model = zeros(num.IVdim,50);
% for i = 1 : 50
%     Ex_model(:,i) = mean(Ex_enroll(:,(i-1)*num.enroll+1:i*num.enroll),2); % i-vector models
% end
% clear i
% 
% answer_eva = [ones(1,50*num.test) zeros(1,50*num.test*49)];


% fprintf('================= cosine ====================\n');
% scores_CDS = [];
% scores_CDS.all =1 - pdist2(Ex_test',Ex_model','cosine');
% 
% scores_CDS.true = [];
% scores_CDS.impostor = [];
% 
% for a =  1 : 50
%     for b = 1 : 50
%         A = scores_CDS.all(num.test*(a-1)+1 : num.test*(a-1)+num.test , b);
%         if a == b
%             scores_CDS.true = [scores_CDS.true ; A(:)];
%         elseif a~=b
%             scores_CDS.impostor = [scores_CDS.impostor ; A(:)];
%         end
%     end
% end  
% clear a b A
% scores_CDS = [scores_CDS.true;scores_CDS.impostor];
% [eer_baseline_FA400_CDS,~,~,dcf_vox_baseline_FA400_CDS,dcf_king_baseline_FA400_CDS]=compute_eer(scores_CDS,answer_eva,false);
% hold on

% fprintf('================= LDA+cosine ====================\n');
% ldaDim = 120;
% [V,D] = lda(Ex_dev, num.label_dev(:));
% finalDev = V(:, 1:ldaDim)' * Ex_dev;
% finalModel = V(:, 1:ldaDim)' * Ex_model;
% finalTest = V(:, 1:ldaDim)' * Ex_test;

% scores_LDA_CDS = [];
% scores_LDA_CDS.all =1 - pdist2(finalTest',finalModel','cosine');
% 
% scores_LDA_CDS.true = [];
% scores_LDA_CDS.impostor = [];
% 
% for a =  1 : 50
%     for b = 1 : 50
%         A = scores_LDA_CDS.all(num.test*(a-1)+1 : num.test*(a-1)+num.test , b);
%         if a == b
%             scores_LDA_CDS.true = [scores_LDA_CDS.true ; A(:)];
%         elseif a~=b
%             scores_LDA_CDS.impostor = [scores_LDA_CDS.impostor ; A(:)];
%         end
%     end
% end  
% clear a b A
% scores_LDA_CDS = [scores_LDA_CDS.true;scores_LDA_CDS.impostor];
% [eer_baseline_FA800_LDA120_CDS(iTers),~,~,dcf_vox_baseline_FA800_LDA120_CDS(iTers)]=compute_eer(scores_LDA_CDS,answer_eva,false);


fprintf('================= PLDA ======================\n');
pLDA_baseline = gplda_em(Ex_dev, num.label_dev(:), num.Zdim, 100);

scores_PLDA = [];
scores_PLDA.all = (score_gplda_trials(pLDA_baseline, Ex_model, Ex_test))';

scores_PLDA.true = [];
scores_PLDA.impostor = [];

for a =  1 : 50
    for b = 1 : 50
        A = scores_PLDA.all(num.test*(a-1)+1 : num.test*(a-1)+num.test , b);
        if a == b
            scores_PLDA.true = [scores_PLDA.true ; A(:)];
        elseif a~=b
            scores_PLDA.impostor  = [scores_PLDA.impostor  ; A(:)];
        end
    end
end
clear a b A
scores_PLDA = [scores_PLDA.true;scores_PLDA.impostor];
[eer_baseline_FA400_PLDA100,~,~,dcf_vox_baseline_FA400_PLDA100,dcf_king_baseline_FA400_PLDA100]=compute_eer(scores_PLDA,answer_eva,1);


% fprintf('================= LDA+PLDA ======================\n');
% pLDA_LDA_baseline = gplda_em(finalDev, num.label_dev(:), ldaDim, 30);
% 
% scores_LDA_PLDA = [];
% scores_LDA_PLDA.all = (score_gplda_trials(pLDA_LDA_baseline, finalModel, finalTest))';
% 
% scores_LDA_PLDA.true = [];
% scores_LDA_PLDA.impostor = [];
% 
% for a =  1 : 50
%     for b = 1 : 50
%         A = scores_LDA_PLDA.all(num.test*(a-1)+1 : num.test*(a-1)+num.test , b);
%         if a == b
%             scores_LDA_PLDA.true = [scores_LDA_PLDA.true ; A(:)];
%         elseif a~=b
%             scores_LDA_PLDA.impostor  = [scores_LDA_PLDA.impostor  ; A(:)];
%         end
%     end
% end
% clear a b A
% scores_LDA_PLDA = [scores_LDA_PLDA.true;scores_LDA_PLDA.impostor];
% [eer_baseline_FA800_LDA120_PLDA120(iTers),~,~,dcf_vox_baseline_FA800_LDA120_PLDA120(iTers)]=compute_eer(scores_LDA_PLDA,answer_eva,false);

fprintf('===== Summaray ======= Baseline ==============\n');
fprintf('+-------------------------------------------+\n');
fprintf('|    Method     |    EER(%%)   | Min DCF_king |\n');
fprintf('+---------------+-------------+-------------+\n');
fprintf('|FEFA(%d)+CDS     |    %2.2f    |   %2.4f    |\n', num.IVdim,eer_baseline_FA400_CDS,dcf_king_baseline_FA400_CDS);
fprintf('|FEFA(%d)+PLDA(%d)|    %2.2f    |   %2.4f    |\n', num.IVdim,num.Zdim,eer_baseline_FA400_PLDA100,dcf_king_baseline_FA400_PLDA100);
fprintf('+---------------+-------------+-------------+\n');








