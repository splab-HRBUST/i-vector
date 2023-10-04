%% ------------------------------------------------loading voxCeleb dataset
% fprintf('loading data\n');
% 
% speaker_id_train = importdata('/data/chenchen/data/voxceleb/IVs_Verification/train_speaker_id.mat');
% speaker_id_test = importdata('/data/chenchen/data/voxceleb/IVs_Verification/test_speaker_id.mat');
% 
% for nFloder = 1 : size(speaker_id_train,1)
%     mkdir('/data/chenchen/matlab/program/Speaker_Verification_Toolbox_v1.c/data/VoxCeleb/mfcc/mfcc_train',speaker_id_train{nFloder,1});
% end
% 
% for nFloder = 1 : size(speaker_id_test,1)
%     mkdir('/data/chenchen/matlab/program/Speaker_Verification_Toolbox_v1.c/data/VoxCeleb/mfcc/mfcc_test',speaker_id_test{nFloder,1});
% end
% 
% % % MFCC as the input
% files_dev = importdata('/data/chenchen/data/voxceleb/voxceleb1_mfcc/train_files.txt');
% files_eva = importdata('/data/chenchen/data/voxceleb/voxceleb1_mfcc/test_files.txt');
% 
% 
% tic;
% for nFile=1:size(files_dev,1)
%     path = fullfile('/data/corpus/VoxCeleb/voxceleb1_wav',char([files_dev{nFile}(31:end-3) 'wav']));
%     mmm = mfcc(path);
%     savePath = fullfile('/data/chenchen/matlab/program/Speaker_Verification_Toolbox_v1.c/data/VoxCeleb/mfcc/mfcc_train',char(files_dev{nFile}(31:end)));
%     save(savePath,'mmm','-v7.3');
% end
% toc; 
% 
% tic;
% for nFile=1:size(files_eva,1)
%     path = fullfile('/data/corpus/VoxCeleb/voxceleb1_wav',char([files_eva{nFile}(30:end-3) 'wav']));
%     mmm = mfcc(path);
%     savePath = fullfile('/data/chenchen/matlab/program/Speaker_Verification_Toolbox_v1.c/data/VoxCeleb/mfcc/mfcc_test',char(files_eva{nFile}(30:end)));
%     save(savePath,'mmm','-v7.3');
% end
% toc; 
% 
% spk_cell_dev = importdata('/data/chenchen/data/voxceleb/IVs_Verification/train_speaker_v_cell.mat');
% spk_cell_eva = importdata('/data/chenchen/data/voxceleb/IVs_Verification/test_speaker_v_cell.mat');
% 
% for nFile = 1 : size(files_dev,1)
% % parfor nSpk = 1 : size(files_dev,1)
%     mfccFilePaths = fullfile('/data/chenchen/matlab/program/Speaker_Verification_Toolbox_v1.c/data/VoxCeleb/mfcc/mfcc_train',char(files_dev{nFile}(31:end)));
%     mfcc_dev_cell(nFile) = importdata(mfccFilePaths);
% end
% save('/data/chenchen/matlab/program/Speaker_Verification_Toolbox_v1.c/data/VoxCeleb/mfcc_dev_cell.mat','mfcc_dev_cell','-v7.3');
% 
% for nFile = 1 : size(files_eva,1)
% % parfor nSpk = 1 : size(files_eva,1)
%     mfccFilePaths = fullfile('/data/chenchen/matlab/program/Speaker_Verification_Toolbox_v1.c/data/VoxCeleb/mfcc/mfcc_test',char(files_eva{nFile}(30:end)));
%     mfcc_eva_cell(nFile) = importdata(mfccFilePaths);
% end
% save('/data/chenchen/matlab/program/Speaker_Verification_Toolbox_v1.c/data/VoxCeleb/mfcc_eva_cell.mat','mfcc_eva_cell','-v7.3');
% 
% label_dev = [];
% spk_id = 1;
% for i = 1 : size(spk_cell_dev,1)
%     spk_dev_num = size(spk_cell_dev{i},2);   
%     label_dev = [label_dev spk_id*ones(1,spk_dev_num)];
%     spk_id = spk_id + 1;
% end
% 
% label_eva = [];
% spk_id = 1;
% for i = 1 : size(spk_cell_eva,1)
%     spk_eva_num = size(spk_cell_eva{i},2);   
%     label_eva = [label_eva spk_id*ones(1,spk_eva_num)];
%     spk_id = spk_id + 1;
% end
% 
% clear nfile i spk_id mfccFilePaths spk_dev_num spk_eva_num


%% train the UBM
% tic; 
% nmix = 1024;
% final_niter = 10;
% ds_factor = 1;
% nworkers = 8;
% ubm = gmm_em(mfcc_dev_cell(:), nmix, final_niter, ds_factor, nworkers);
% fprintf('UBM trained successfully\n');
% toc; 
% 
% % save('/data/chenchen/matlab/program/Speaker_Verification_Toolbox_v1.c/data/VoxCeleb/mfcc_dev_cell.mat','mfcc_dev_cell','-v7.3');
% % save('/data/chenchen/matlab/program/Speaker_Verification_Toolbox_v1.c/data/VoxCeleb/mfcc_eva_cell.mat','mfcc_eva_cell','-v7.3');
% save('/data/chenchen/matlab/program/Speaker_Verification_Toolbox_v1.c/data/VoxCeleb/ubm.mat','ubm','-v7.3');
% 
% ubm = importdata('/data/chenchen/matlab/program/Speaker_Verification_Toolbox_v1.c/data/VoxCeleb/ubm.mat');

%% adapt the GMM
% for nFloder = 1 : size(speaker_id_train,1)
%     mkdir('/data/chenchen/matlab/program/Speaker_Verification_Toolbox_v1.c/data/VoxCeleb/gmm/gmm_train',speaker_id_train{nFloder,1});
% end
% 
% for nFloder = 1 : size(speaker_id_test,1)
%     mkdir('/data/chenchen/matlab/program/Speaker_Verification_Toolbox_v1.c/data/VoxCeleb/gmm/gmm_test',speaker_id_test{nFloder,1});
% end
% 
% tic;
% gmm = [];
% map_tau = 10.0;
% config = 'm';
% M_dev = zeros(61440,size(files_dev,1));
% for nFile=1:size(files_dev,1)
%     path = fullfile('/data/chenchen/matlab/program/Speaker_Verification_Toolbox_v1.c/data/VoxCeleb/mfcc/mfcc_train',char(files_dev{nFile}(31:end)));
%     ccc = importdata(path);
%
%     gmm = mapAdapt(ccc, ubm, map_tau, config);
%     savePath = fullfile('/data/chenchen/matlab/program/Speaker_Verification_Toolbox_v1.c/data/VoxCeleb/gmm/gmm_train',char(files_dev{nFile}(31:end)));
%     save(savePath,'gmm','-v7.3');
%     gmm = importdata(savePath);
%     M_dev(:,nFile) = gmm.mu(:);
% end
% toc; 
% save('/data/chenchen/matlab/program/Speaker_Verification_Toolbox_v1.c/data/VoxCeleb/M_dev.mat','M_dev','-v7.3');
% 
% 
% tic;
% gmm = [];
% map_tau = 10.0;
% config = 'm';
% M_eva = zeros(61440,size(files_eva,1));
% for nFile=1:size(files_eva,1)
%     path = fullfile('/data/chenchen/matlab/program/Speaker_Verification_Toolbox_v1.c/data/VoxCeleb/mfcc/mfcc_test',char(files_eva{nFile}(30:end)));
%     ccc = importdata(path);
%     
%     gmm = mapAdapt(ccc, ubm, map_tau, config);
%     savePath = fullfile('/data/chenchen/matlab/program/Speaker_Verification_Toolbox_v1.c/data/VoxCeleb/gmm/gmm_test',char(files_eva{nFile}(30:end)));
%     save(savePath,'gmm','-v7.3'); 
%     gmm = importdata(savePath);
%     M_eva(:,nFile) = gmm.mu(:);    
% end
% toc;
% save('/data/chenchen/matlab/program/Speaker_Verification_Toolbox_v1.c/data/VoxCeleb/M_eva.mat','M_eva','-v7.3');
% % 
% M_dev = importdata('/data/chenchen/matlab/program/Speaker_Verification_Toolbox_v1.c/data/VoxCeleb/M_dev.mat');
% M_eva = importdata('/data/chenchen/matlab/program/Speaker_Verification_Toolbox_v1.c/data/VoxCeleb/M_eva.mat');
% ubm = importdata('/data/chenchen/matlab/program/Speaker_Verification_Toolbox_v1.c/data/VoxCeleb/ubm.mat');

%% train the FA
num.label_dev = importdata('/data/chenchen/matlab/program/Speaker_Verification_Toolbox_v1.c/data/VoxCeleb/label_dev.mat'); % for PLDA
[num.SVdim, num.samples] = size(M_dev);
num.IVdim = 400; % dimension of i-vector
num.Zdim = 200; % dimension of latent vector for PLDA
% 
%  
num.nIters = 1;

%% model training
for iTers = 1 : 15
    if iTers == 1
        patameters_FA_baseline = fa_em(M_dev, ubm, num);
    else
        patameters_FA_baseline = fa_em(M_dev, ubm, num, patameters_FA_baseline);
    end
        
%% prepare for evaluation
Ex_dev = pinv(patameters_FA_baseline.L)*patameters_FA_baseline.B*(M_dev-patameters_FA_baseline.m);
Ex_eva = pinv(patameters_FA_baseline.L)*patameters_FA_baseline.B*(M_eva-patameters_FA_baseline.m);
% 
answer_eva=[];
target1_IVs=[];
target2_IVs=[];
fverification = fopen('/data/chenchen/data/voxceleb/verification_number_form.txt');
tline = fgetl(fverification);
ver_id=1;
while ischar(tline)
     temp_cell=textscan(tline,'%d %d %d');
     answer_eva(ver_id)=temp_cell{1,1};
     target1_IVs(:,ver_id)=Ex_eva(:,temp_cell{1,2});
     target2_IVs(:,ver_id)=Ex_eva(:,temp_cell{1,3});
     tline = fgetl(fverification);
     ver_id=ver_id+1;
end
fclose(fverification);
clear ver_id tline fverification temp_cell

% ldaDim = 120;
% [V,D] = lda(Ex_dev, num.label_dev(:));
% finalDev = V(:, 1:ldaDim)' * Ex_dev;
% finalTarget1 = V(:, 1:ldaDim)' * target1_IVs;
% finalTarget2 = V(:, 1:ldaDim)' * target2_IVs;

% Ex_LDA_dev = V(:, 1:ldaDim)' * Ex_dev;
% Ex_LDA_eva = V(:, 1:ldaDim)' * Ex_eva;


%% evaluation
fprintf('================= cosine ====================\n');
scores_CDS = 1 - pdist2(target2_IVs',target1_IVs','cosine');
scores_CDS = diag(scores_CDS);
[eer_baseline_FA400_CDS(iTers),~,~,dcf_vox_baseline_FA400_CDS(iTers)]=compute_eer(scores_CDS,answer_eva,false);
% result_baseline.eer_baseline_FA400_CDS = eer_baseline_FA400_CDS;
% result_baseline.dcf_vox_baseline_FA400_CDS = dcf_vox_baseline_FA400_CDS;

% fprintf('================= LDA+cosine ====================\n');
% scores_LDA_CDS = 1 - pdist2(finalTarget2',finalTarget1','cosine');
% scores_LDA_CDS = diag(scores_LDA_CDS);
% [eer_baseline_FA400_LDA200_CDS(iTers),~,~,dcf_vox_baseline_FA400_LDA200_CDS(iTers)]=compute_eer(scores_LDA_CDS,answer_eva,false);
% % result_baseline.eer_baseline_FA400_LDA200_CDS = eer_baseline_FA400_LDA200_CDS;
% % result_baseline.dcf_vox_baseline_FA400_LDA200_CDS = dcf_vox_baseline_FA400_LDA200_CDS;

% fprintf('================= PLDA ======================\n');
% pLDA_baseline = gplda_em(Ex_dev, num.label_dev(:), num.Zdim, 50);
% scores_PLDA=score_gplda_trials(pLDA_baseline,target1_IVs,target2_IVs);
% scores_PLDA=diag(scores_PLDA);
% [eer_baseline_FA400_PLDA200(iTers),~,~,dcf_vox_Baseline_FA400_PLDA200(iTers)]=compute_eer(scores_PLDA,answer_eva,false);
% result_baseline.eer_baseline_FA400_PLDA200 = eer_baseline_FA400_PLDA200;
% result_baseline.dcf_vox_Baseline_FA400_PLDA200 = dcf_vox_Baseline_FA400_PLDA200;

% [~, ia, ic] = unique([1 :1 :size(Ex_eva,2)], 'stable');
% spk_counts = histc(ic, 1 : numel(ia)); % # sessions per speaker
% clear ia ic
% 
% [Ey_dev, ~] = expectation_plda(Ex_dev, pLDA_baseline.Phi, pLDA_baseline.Sigma, spk_counts);
% [Ey_eva, ~] = expectation_plda(Ex_eva, pLDA_baseline.Phi, pLDA_baseline.Sigma, spk_counts);

% fprintf('================= LDA+PLDA ======================\n');
% pLDA_LDA = gplda_em(finalDev, num.label_dev(:), ldaDim, 50);
% scores_LDA_PLDA=score_gplda_trials(pLDA_LDA,finalTarget1,finalTarget2);
% scores_LDA_PLDA=diag(scores_LDA_PLDA);
% [eer_baseline_FA400_LDA200_PLDA200(iTers),~,~,dcf_vox_Baseline_FA400_LDA200_PLDA200(iTers)]=compute_eer(scores_LDA_PLDA,answer_eva,false);
% result_baseline.eer_baseline_FA400_LDA200_PLDA200 = eer_baseline_FA400_LDA200_PLDA200;
% result_baseline.dcf_vox_Baseline_FA400_LDA200_PLDA200 = dcf_vox_Baseline_FA400_LDA200_PLDA200;
% 
fprintf('===== Summaray ======= Baseline =====%d=========\n',iTers);
fprintf('+-------------------------------------------+\n');
fprintf('|    Method     |    EER(%%)   | Min DCF_vox |\n');
fprintf('+---------------+-------------+-------------+\n');
fprintf('|      CDS      |    %2.2f    |   %2.4f    |\n', eer_baseline_FA400_CDS(iTers),dcf_vox_baseline_FA400_CDS(iTers));
% fprintf('|    LDA+CDS    |    %2.2f    |   %2.4f    |\n', eer_baseline_FA400_LDA200_CDS(iTers),dcf_vox_baseline_FA400_LDA200_CDS(iTers));
% fprintf('|FA(%d)+PLDA(%d)|    %2.2f    |   %2.4f    |\n', num.IVdim,num.Zdim,eer_baseline_FA400_PLDA200(iTers),dcf_vox_Baseline_FA400_PLDA200(iTers));
% fprintf('|FA(%d)+LDA+PLDA(%d)|    %2.2f    |   %2.4f    |\n', num.IVdim,num.Zdim,eer_baseline_FA400_LDA200_PLDA200(iTers),dcf_vox_Baseline_FA400_LDA200_PLDA200(iTers));
fprintf('+---------------+-------------+-------------+\n');

end
% save('result_baseline.mat','result_baseline','-v7.3');

%% EER and Min DCF
% [eer_CDS,dcf_CDS] = plotdeteer(scores_CDS.true, scores_CDS.impostor,':m');
% [eer_PLDA,dcf_PLDA] =plotdeteer(scores_PLDA.true, scores_PLDA.impostor ,'-.b');

% h = legend('CDS','PLDA');
% set(h,'Fontsize',12); 
% clear h

% [eer_CDS,nMin_CDS] = min(eer_baseline_FA800_CDS);
% dcf_vox_CDS = dcf_vox_baseline_FA800_CDS(nMin_CDS);
% 
% % [eer_PLDA, nMin_PLDA] = min(eer_baseline_FA400_PLDA400);
% % dcf_vox_PLDA = dcf_vox_Baseline_FA800_PLDA200(nMin_PLDA);
% 
% fprintf('======================================================================\n');
% fprintf('===== Summaray ======= Baseline =====%d=========\n',nMin_PLDA);
% fprintf('+---------------------------------------------+\n');
% fprintf('|      Method     |    EER(%%)   | Min DCF_vox |\n');
% fprintf('+-----------------+-------------+-------------+\n');
% fprintf('|      CDS      |    %2.2f    |   %2.4f    |\n', eer_CDS,dcf_vox_CDS);
% % fprintf('|FA(%d)-PLDA(%d)|     %2.2f    |   %2.4f    |\n', num.IVdim,num.Zdim,eer_PLDA,dcf_vox_PLDA);
% fprintf('+-----------------+-------------+-------------+\n');
% fprintf('Finished :D\n');
% 
% clear eer_PLDA dcf10_PLDA

%% plot similarity matrix
% finalEva = V(:, 1:ldaDim)' * Ex_eva;

% subplot(1,2,1)
% imagesc(Ex_eva'*Ex_eva-diag(diag(Ex_eva'*Ex_eva)));colormap('gray')
% 
% subplot(1,2,2)
% imagesc(finalEva'*finalEva-diag(diag(finalEva'*finalEva)));colormap('gray')

% subplot(1,3,3)
% imagesc(xxx_alpha_ind'*xxx_alpha_ind-diag(diag(xxx_alpha_ind'*xxx_alpha_ind)));colormap('gray')








