%% ------------------------------------------------TIMIT database 
% GMM supervectors
M_dev_all = importdata('E:/DataPrep/timit170/M_dev.mat');
M_enroll = importdata('E:/DataPrep/timit170/M_enroll.mat');
M_test   = importdata('E:/DataPrep/timit170/M_test.mat');

num.spk_dev = 170;
num.spk_eva = 100;

num.dev = 10;
num.enroll = 9;
num.test = 1;
num.label_dev = repmat([1:1:170],10,1); % for PLDA

%% ------------------------------------------------Random choose development data
num.ch = 10;

fprintf('Development data: %d\n', num.ch);
M_dev = [];
for nSpk = 1 : num.spk_dev
    M_dev(:,(nSpk-1)*num.ch+1:nSpk*num.ch) =  M_dev_all(:,(nSpk-1)*10+1:(nSpk-1)*10+num.ch);
end
clear nSpk

%% ---------------------------------------------Train the BPCA-based TVM
num.IVdim = 150;
num.nIters = 1;

for iTers = 1 : 25
    if iTers == 1
        patameters_BPCA = bpca_em(M_dev,num);
    else
        patameters_BPCA = bpca_em(M_dev,num,patameters_BPCA);
    end

    %% ---------------------------------------------Evaluation
    Ex_dev = pinv(patameters_BPCA.L)*patameters_BPCA.B*(bsxfun(@minus,M_dev,patameters_BPCA.m));
    Ex_enroll = pinv(patameters_BPCA.L)*patameters_BPCA.B*(bsxfun(@minus,M_enroll,patameters_BPCA.m));
    Ex_test   = pinv(patameters_BPCA.L)*patameters_BPCA.B*(bsxfun(@minus,M_test,patameters_BPCA.m));

    Ex_model = zeros(num.IVdim,num.spk_eva);
    for i = 1 : num.spk_eva
        Ex_model(:,i) = mean(Ex_enroll(:,(i-1)*num.enroll+1:i*num.enroll),2); % the means of i-vectors for one speaker
    end
    clear i
    
    answer_eva = [ones(1,100*num.test) zeros(1,100*num.test*99)];
    
    
    fprintf('================= cosine ====================\n');
    scores_BPCA_CDS = [];
    scores_BPCA_CDS.all =1 - pdist2(Ex_test',Ex_model','cosine');
    
    scores_BPCA_CDS.true = [];
    scores_BPCA_CDS.impostor = [];
    
    for a =  1 : 100
        for b = 1 : 100
            A = scores_BPCA_CDS.all(num.test*(a-1)+1 : num.test*(a-1)+num.test , b);
            if a == b
                scores_BPCA_CDS.true = [scores_BPCA_CDS.true ; A(:)];
            elseif a~=b
                scores_BPCA_CDS.impostor = [scores_BPCA_CDS.impostor ; A(:)];
            end
        end
    end
    clear a b A
    scores_BPCA_CDS = [scores_BPCA_CDS.true;scores_BPCA_CDS.impostor];
    [eer_BPCA400_CDS(iTers),~,dcf_timit_BPCA400_CDS(iTers)]=compute_eer(scores_BPCA_CDS,answer_eva,false);
    
    
    fprintf('================= LDA+cosine ====================\n');
    ldaDim =100;
    [V,D] = lda(Ex_dev, num.label_dev(:));
    finalDev = V(:, 1:ldaDim)' * Ex_dev;
    finalModel = V(:, 1:ldaDim)' * Ex_model;
    finalTest = V(:, 1:ldaDim)' * Ex_test;

    scores_LDA_CDS = [];
    scores_LDA_CDS.all =1 - pdist2(finalTest',finalModel','cosine');
    
    scores_LDA_CDS.true = [];
    scores_LDA_CDS.impostor = [];
    
    for a =  1 : 100
        for b = 1 : 100
            A = scores_LDA_CDS.all(num.test*(a-1)+1 : num.test*(a-1)+num.test , b);
            if a == b
                scores_LDA_CDS.true = [scores_LDA_CDS.true ; A(:)];
            elseif a~=b
                scores_LDA_CDS.impostor = [scores_LDA_CDS.impostor ; A(:)];
            end
        end
    end
    clear a b A
    scores_LDA_CDS = [scores_LDA_CDS.true;scores_LDA_CDS.impostor];
    [eer_baseline_BPCA400_LDA170_CDS(iTers),~,dcf_timit_baseline_BPCA400_LDA170_CDS(iTers)]=compute_eer(scores_LDA_CDS,answer_eva,false);
    
  
    fprintf('================= PLDA ======================\n');
    num.Zdim = 100;
    pLDA_bpca = gplda_em(Ex_dev, num.label_dev(:), num.Zdim, 50);
    
    scores_BPCA_PLDA = [];
    scores_BPCA_PLDA.all = (score_gplda_trials(pLDA_bpca, Ex_model, Ex_test))';
    
    scores_BPCA_PLDA.true = [];
    scores_BPCA_PLDA.impostor = [];
    
    for a =  1 : 100
        for b = 1 : 100
            A = scores_BPCA_PLDA.all(num.test*(a-1)+1 : num.test*(a-1)+num.test , b);
            if a == b
                scores_BPCA_PLDA.true = [scores_BPCA_PLDA.true ; A(:)];
            elseif a~=b
                scores_BPCA_PLDA.impostor  = [scores_BPCA_PLDA.impostor  ; A(:)];
            end
        end
    end
    clear a b A
    scores_BPCA_PLDA = [scores_BPCA_PLDA.true;scores_BPCA_PLDA.impostor];
    [eer_BPCA400_PLDA170(iTers),~,dcf_timit_BPCA400_PLDA170(iTers)]=compute_eer(scores_BPCA_PLDA,answer_eva,false);
    
    
    fprintf('================= LDA+PLDA ======================\n');
    pLDA_LDA_baseline = gplda_em(finalDev, num.label_dev(:), ldaDim, 50);
    
    scores_LDA_PLDA = [];
    scores_LDA_PLDA.all = (score_gplda_trials(pLDA_LDA_baseline, finalModel, finalTest))';
    
    scores_LDA_PLDA.true = [];
    scores_LDA_PLDA.impostor = [];
    
    for a =  1 : 100
        for b = 1 : 100
            A = scores_LDA_PLDA.all(num.test*(a-1)+1 : num.test*(a-1)+num.test , b);
            if a == b
                scores_LDA_PLDA.true = [scores_LDA_PLDA.true ; A(:)];
            elseif a~=b
                scores_LDA_PLDA.impostor  = [scores_LDA_PLDA.impostor  ; A(:)];
            end
        end
    end
    clear a b A
    scores_LDA_PLDA = [scores_LDA_PLDA.true;scores_LDA_PLDA.impostor];
    [eer_baseline_BPCA400_LDA170_PLDA170(iTers),~,dcf_timit_baseline_BPCA400_LDA170_PLDA170(iTers)]=compute_eer(scores_LDA_PLDA,answer_eva,false);
    
    fprintf('===== Summaray ======= Baseline =====%d=========\n',iTers);
    fprintf('+-------------------------------------------+\n');
    fprintf('|    Method     |    EER(%%)   | Min DCF_timit |\n');
    fprintf('+---------------+-------------+-------------+\n');
    fprintf('|BPCA(%d)+CDS     |    %2.2f    |   %2.4f    |\n', num.IVdim, eer_BPCA400_CDS(iTers), dcf_timit_BPCA400_CDS(iTers));
    fprintf('|BPCA(%d)+LDA+CDS    |    %2.2f    |   %2.4f    |\n',num.IVdim, eer_baseline_BPCA400_LDA170_CDS(iTers),dcf_timit_baseline_BPCA400_LDA170_CDS(iTers));
    fprintf('|BPCA(%d)+PLDA(%d)|    %2.2f    |   %2.4f    |\n', num.IVdim,num.Zdim,eer_BPCA400_PLDA170(iTers),dcf_timit_BPCA400_PLDA170(iTers));
    fprintf('|BPCA(%d)+LDA+PLDA(%d)|    %2.2f    |   %2.4f    |\n', num.IVdim,ldaDim,eer_baseline_BPCA400_LDA170_PLDA170(iTers),dcf_timit_baseline_BPCA400_LDA170_PLDA170(iTers));
    fprintf('+---------------+-------------+-------------+\n');

end

fprintf('Finished :D\n');
