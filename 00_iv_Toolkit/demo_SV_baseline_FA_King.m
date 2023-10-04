%% Test list
% dev: (1-150spk)*120wav
% eva: (151-200spk)*120wav
% enroll: 96/120wav
% test:   24/120wav

% ===== Summaray ======= Baseline ====== FA =========
% +------------------------------------------------+
% |       Method      |    EER(%)   | Min DCF_king |
% +-------------------+-------------+--------------+
% |  FA(400)+CDS      |     4.33    |    0.9583    |
% |  FA(400)+PLDA(200)|     3.17    |    0.5850    |
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


%% Adapt GMM
% tic;
% gmm = [];
% map_tau = 10.0;
% config = 'm';
% for nFile=1:size(mfcc_all_cell,2)
%     gmm = mapAdapt(mfcc_all_cell(nFile), ubm, map_tau, config);
% %     savePath = fullfile('/data/chenchen/matlab/program/ASR_Toolbox/00_Data/King-ASR-010/Mix_60_1024/gmm',char(AllFiles{nFile}));
% %     save(savePath,'gmm','-v7.3');
%     M_all(:,nFile) = gmm.mu(:);
% end
% toc;
% clear map_tau nFile
% % save('/data/chenchen/matlab/program/ASR_Toolbox/00_Data/King-ASR-010/Mix_60_1024/M_all.mat','M_all','-v7.3');
% 
% M_dev = M_all(:,1:150*120);
% M_eva = M_all(:,150*120+1:end);
% % save('/data/chenchen/matlab/program/ASR_Toolbox/00_Data/King-ASR-010/Mix_60_1024/M_dev.mat','M_dev','-v7.3');
% % save('/data/chenchen/matlab/program/ASR_Toolbox/00_Data/King-ASR-010/Mix_60_1024/M_eva.mat','M_eva','-v7.3');
% 
% for nSpk = 1 : 50
%     M_enroll(:,(nSpk-1)*num.enroll+1:nSpk*num.enroll) =  M_eva(:,(nSpk-1)*120+1:(nSpk-1)*120+96);
%     M_test(:,(nSpk-1)*num.test+1  :nSpk*num.test)   =  M_eva(:,(nSpk-1)*120+97:nSpk*120);
% end
% clear nSpk
% % save('/data/chenchen/matlab/program/ASR_Toolbox/00_Data/King-ASR-010/Mix_60_1024/M_enroll.mat','M_enroll','-v7.3');
% % save('/data/chenchen/matlab/program/ASR_Toolbox/00_Data/King-ASR-010/Mix_60_1024/M_test.mat','M_test','-v7.3');
% 
% M_dev    = importdata('/data/chenchen/matlab/program/ASR_Toolbox/00_Data/King-ASR-010/Mix_60_1024/M_dev.mat');
% M_enroll = importdata('M_enroll.mat');
% M_test   = importdata('M_test.mat');
% ubm      = importdata('/data/chenchen/matlab/program/ASR_Toolbox/00_Data/King-ASR-010/Mix_60_1024/model/ubm.mat');

% [num.SVdim, num.samples] = size(M_dev);
num.IVdim = 400; % dimension of i-vector
num.Zdim = 200;  % dimension of latent vector for PLDA
 
% num.label_dev = repmat([1:1:150],120,1); % for PLDA
num.enroll = 96;
num.test = 24;
% ---------------------------------------------迭代次数
num.nIters = 1;


%% ---------------------------------------------train TVS using FA
% for iTers = 1 : 15
%     if iTers == 1
%         patameters_FA_baseline = fa_em(M_dev, ubm, num);
%     else
%         patameters_FA_baseline = fa_em(M_dev, ubm, num, patameters_FA_baseline);
%     end
        
%% Evaluation
% Ex_dev = pinv(patameters_FA_baseline.L)*patameters_FA_baseline.B*(M_dev-patameters_FA_baseline.m);
Ex_enroll = pinv(patameters_FA_baseline.L)*patameters_FA_baseline.B*(M_enroll-patameters_FA_baseline.m);
Ex_test = pinv(patameters_FA_baseline.L)*patameters_FA_baseline.B*(M_test-patameters_FA_baseline.m);

Ex_model = zeros(num.IVdim,50);
for i = 1 : 50
    Ex_model(:,i) = mean(Ex_enroll(:,(i-1)*num.enroll+1:i*num.enroll),2); % i-vector models
end
clear i

answer_eva = [ones(1,50*num.test) zeros(1,50*num.test*49)];


fprintf('================= cosine ====================\n');
scores_CDS = [];
scores_CDS.all =1 - pdist2(Ex_test',Ex_model','cosine');

scores_CDS.true = [];
scores_CDS.impostor = [];

for a =  1 : 50
    for b = 1 : 50
        A = scores_CDS.all(num.test*(a-1)+1 : num.test*(a-1)+num.test , b);
        if a == b
            scores_CDS.true = [scores_CDS.true ; A(:)];
        elseif a~=b
            scores_CDS.impostor = [scores_CDS.impostor ; A(:)];
        end
    end
end  
clear a b A
scores_CDS = [scores_CDS.true;scores_CDS.impostor];
[eer_baseline_FA400_CDS(iTers),~,~,dcf_vox_baseline_FA400_CDS(iTers),dcf_king_baseline_FA400_CDS(iTers)]=compute_eer(scores_CDS,answer_eva,false);


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


% fprintf('================= PLDA ======================\n');
% pLDA_baseline = gplda_em(Ex_dev, num.label_dev(:), num.Zdim, 30);
% 
% scores_PLDA = [];
% scores_PLDA.all = (score_gplda_trials(pLDA_baseline, Ex_model, Ex_test))';
% 
% scores_PLDA.true = [];
% scores_PLDA.impostor = [];
% 
% for a =  1 : 50
%     for b = 1 : 50
%         A = scores_PLDA.all(num.test*(a-1)+1 : num.test*(a-1)+num.test , b);
%         if a == b
%             scores_PLDA.true = [scores_PLDA.true ; A(:)];
%         elseif a~=b
%             scores_PLDA.impostor  = [scores_PLDA.impostor  ; A(:)];
%         end
%     end
% end
% clear a b A
% scores_PLDA = [scores_PLDA.true;scores_PLDA.impostor];
% [eer_baseline_FA400_PLDA200(iTers),~,~,dcf_vox_baseline_FA400_PLDA200(iTers)]=compute_eer(scores_PLDA,answer_eva,false);


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

% fprintf('===== Summaray ======= Baseline =====%d=========\n',iTers);
% fprintf('+-------------------------------------------+\n');
% fprintf('|    Method     |    EER(%%)   | Min DCF_vox |\n');
% fprintf('+---------------+-------------+-------------+\n');
% fprintf('|      CDS      |    %2.2f    |   %2.4f    |\n', eer_baseline_FA400_CDS(iTers),dcf_vox_baseline_FA400_CDS(iTers));
% % fprintf('|    LDA+CDS    |    %2.2f    |   %2.4f    |\n', eer_baseline_FA800_LDA120_CDS(iTers),dcf_vox_baseline_FA800_LDA120_CDS(iTers));
% % fprintf('|FA(%d)+PLDA(%d)|    %2.2f    |   %2.4f    |\n', num.IVdim,num.Zdim,eer_baseline_FA400_PLDA200(iTers),dcf_vox_baseline_FA400_PLDA200(iTers));
% % fprintf('|FA(%d)+LDA+PLDA(%d)|    %2.2f    |   %2.4f    |\n', num.IVdim,ldaDim,eer_baseline_FA800_LDA120_PLDA120(iTers),dcf_vox_baseline_FA800_LDA120_PLDA120(iTers));
% fprintf('+---------------+-------------+-------------+\n');

% end

% [eer_CDS,nMin_CDS] = min(eer_baseline_FA400_LDA200_CDS);
% dcf_vox_CDS = dcf_vox_baseline_FA400_CDS(nMin_CDS);
% 
% [eer_PLDA, nMin_PLDA] = min(eer_baseline_FA400_PLDA200);
% dcf_vox_PLDA = dcf_vox_baseline_FA400_PLDA200(nMin_PLDA);
% 
% fprintf('======================================================================\n');
% fprintf('===== Summaray ======= Baseline =====%d=========\n',nMin_PLDA);
% fprintf('+---------------------------------------------+\n');
% fprintf('|      Method     |    EER(%%)   | Min DCF_vox |\n');
% fprintf('+-----------------+-------------+-------------+\n');
% fprintf('|      CDS      |    %2.2f    |   %2.4f    |\n', eer_CDS,dcf_vox_CDS);
% fprintf('|FA(%d)-PLDA(%d)|     %2.2f    |   %2.4f    |\n', num.IVdim,num.Zdim,eer_PLDA,dcf_vox_PLDA);
% fprintf('+-----------------+-------------+-------------+\n');
% fprintf('Finished :D\n');
% 
% clear eer_CDS dcf_vox_CDS eer_PLDA dcf_vox_PLDA






