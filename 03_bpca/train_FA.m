%% loading VoxCeleb1 dataset
fprintf('loading data\n');

%% UBM
ubm = importdata('/home/rongyafeng/matlab_experiments/export/voxceleb1/ubm.mat');

%% GMM
M_dev = importdata('/home/rongyafeng/matlab_experiments/export/voxceleb1/M_dev.mat');
M_eva = importdata('/home/rongyafeng/matlab_experiments/export/voxceleb1/M_eva.mat');

%% train the FA
num.label_dev = importdata('/home/rongyafeng/matlab_experiments/export/voxceleb1/label_dev.mat'); % for PLDA
[num.SVdim, num.samples] = size(M_dev);
num.IVdim = 400; % dimension of i-vector
num.Zdim = 200; % dimension of latent vector for PLDA

num.nIters = 1;

%% model training
for iTers = 1 : 15
    if iTers == 1
        parameters_FA = fa_em(M_dev, ubm, num);
    else
        parameters_FA = fa_em(M_dev, ubm, num, parameters_FA);
    end
        
    %% prepare for evaluation
    IV_dev = pinv(parameters_FA.L)*parameters_FA.B*bsxfun(@minus, M_dev, parameters_FA.m);
    IV_eva = pinv(parameters_FA.L)*parameters_FA.B*bsxfun(@minus, M_eva, parameters_FA.m);

    answer_eva=[];
    target1_IVs=[];
    target2_IVs=[];
    fverification = fopen('/home/rongyafeng/matlab_experiments/data_prep/voxceleb1/verification_number_form.txt');
    tline = fgetl(fverification);
    ver_id=1;
    while ischar(tline)
         temp_cell=textscan(tline,'%d %d %d');
         answer_eva(ver_id)=temp_cell{1,1};
         target1_IVs(:,ver_id)=IV_eva(:,temp_cell{1,2});
         target2_IVs(:,ver_id)=IV_eva(:,temp_cell{1,3});
         tline = fgetl(fverification);
         ver_id=ver_id+1;
    end
    fclose(fverification);
    clear ver_id tline fverification temp_cell

    ldaDim = 200;
    [V,D] = lda(IV_dev, num.label_dev(:));
    finalDev = V(:, 1:ldaDim)' * IV_dev;
    finalTarget1 = V(:, 1:ldaDim)' * target1_IVs;
    finalTarget2 = V(:, 1:ldaDim)' * target2_IVs;

    %% evaluation
    fprintf('================= cosine ====================\n');
    scores_CDS = 1 - pdist2(target2_IVs',target1_IVs','cosine');
    scores_CDS = diag(scores_CDS);
    [eer_FA400_CDS(iTers),~,dcf_vox_FA400_CDS(iTers)]=compute_eer(scores_CDS,answer_eva,false);

    fprintf('================= LDA+cosine ====================\n');
    scores_LDA_CDS = 1 - pdist2(finalTarget2',finalTarget1','cosine');
    scores_LDA_CDS = diag(scores_LDA_CDS);
    [eer_FA400_LDA200_CDS(iTers),~,dcf_vox_FA400_LDA200_CDS(iTers)]=compute_eer(scores_LDA_CDS,answer_eva,false);

    fprintf('================= PLDA ======================\n');
    pLDA_baseline = gplda_em(IV_dev, num.label_dev(:), num.Zdim, 50);
    scores_PLDA=score_gplda_trials(pLDA_baseline,target1_IVs,target2_IVs);
    scores_PLDA=diag(scores_PLDA);
    [eer_FA400_PLDA200(iTers),~,dcf_vox_FA400_PLDA200(iTers)]=compute_eer(scores_PLDA,answer_eva,false);

    fprintf('================= LDA+PLDA ======================\n');
    pLDA_LDA = gplda_em(finalDev, num.label_dev(:), ldaDim, 50);
    scores_LDA_PLDA=score_gplda_trials(pLDA_LDA,finalTarget1,finalTarget2);
    scores_LDA_PLDA=diag(scores_LDA_PLDA);
    [eer_FA400_LDA200_PLDA200(iTers),~,dcf_vox_FA400_LDA200_PLDA200(iTers)]=compute_eer(scores_LDA_PLDA,answer_eva,false);

    fprintf('===== Summaray ======= Baseline =====%d=========\n',iTers);
    fprintf('+-------------------------------------------+\n');
    fprintf('|    Method     |    EER(%%)   | Min DCF_vox |\n');
    fprintf('+---------------+-------------+-------------+\n');
    fprintf('|      CDS      |    %2.2f    |   %2.4f    |\n', eer_FA400_CDS(iTers),dcf_vox_FA400_CDS(iTers));
    fprintf('|    LDA+CDS    |    %2.2f    |   %2.4f    |\n', eer_FA400_LDA200_CDS(iTers),dcf_vox_FA400_LDA200_CDS(iTers));
    fprintf('|FA(%d)+PLDA(%d)|    %2.2f    |   %2.4f    |\n', num.IVdim,num.Zdim,eer_FA400_PLDA200(iTers),dcf_vox_FA400_PLDA200(iTers));
    fprintf('|FA(%d)+LDA+PLDA(%d)|    %2.2f    |   %2.4f    |\n', num.IVdim,num.Zdim,eer_FA400_LDA200_PLDA200(iTers),dcf_vox_FA400_LDA200_PLDA200(iTers));
    fprintf('+---------------+-------------+-------------+\n');

end