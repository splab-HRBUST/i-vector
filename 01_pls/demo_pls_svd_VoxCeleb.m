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

R = [100 : 100: 1200];%num.spk_dev;


for nR = 7  : size(R,2)

    V_X = patameters_PLS.V(:,1:R(nR))*(patameters_PLS.T(:,1:R(nR))'*patameters_PLS.V(:,1:R(nR)))^-1 ;

    W_dev = m_dev' * V_X;
    W_eva = m_eva' * V_X;

    IV_eva = W_eva';

    answer_eva=[];
    target1_IVs=[];
    target2_IVs=[];
    fverification = fopen('/data/chenchen/matlab/program/ASR_Toolbox/00_Data/VoxCeleb1/verification_number_form.txt');
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
%     
    
% %     y_pre{nR} = m_dev'*V_X*diag(patameters_PLS.B(1:R(nR)))*patameters_PLS.Q(:,1:R(nR))';  % should be this one, but wrong, don't know why
% %     y_pre{nR} = m_dev'*V_X*diag(patameters_PLS.B(1:R(nR)))*patameters_PLS.C(:,1:R(nR))'; % wrong
% 
%     y_pre{nR} = m_dev'*patameters_PLS.V(:,1:R(nR))*diag(patamesters_PLS.B(1:R(nR)))*patameters_PLS.C(:,1:R(nR))'; % the best, but don't know why
% 
%     Res_dev = y_dev'-y_pre{nR};
%     S_RESS(nR) = sum(sum(Res_dev.*Res_dev,2));
%     acc{nR} = top_N_acc(y_pre{nR},label_dev,5);
    
    

    fprintf('================= cosine ====================\n');
    scores_PLS_CDS_vox = 1 - pdist2(target2_IVs',target1_IVs','cosine');
    scores_PLS_CDS_vox = diag(scores_PLS_CDS_vox);
    [eer_PLS_CDS(nR),~,~,dcf_vox_PLS_CDS(nR)]=compute_eer(scores_PLS_CDS_vox,answer_eva,0);
%     hold on


    fprintf('================= PLDA ======================\n');
    nIters = 200;
    Zdim = [100:100:600];
% %     nZ = 1;

    for nZ = 1 : 6
        num.Zdim = Zdim(nZ);
        pLDA = gplda_em(W_dev', label_dev(:), num.Zdim, nIters);

        scores_PLS_PLDA=score_gplda_trials(pLDA,target1_IVs,target2_IVs);
        scores_PLS_PLDA=diag(scores_PLS_PLDA);
        [eer_PLS_PLDA(nR,nZ),~,~,dcf_vox_PLS_PLDA(nR,nZ)]=compute_eer(scores_PLS_PLDA,answer_eva,1);

        fprintf('========== Summaray =========== PLS =====%d====\n',R(nR));
        fprintf('+-----------------------------------------------+\n');
        fprintf('|      Method      |    EER(%%)   | Min DCF_vox |\n');
        fprintf('+------------------+-------------+-----------+\n');
        fprintf('|PLS(%d)+CDS     |    %2.2f    |   %2.4f    |\n', R(nR),           eer_PLS_CDS(nR),  dcf_vox_PLS_CDS(nR));
        fprintf('|PLS(%d)+PLDA(%d)|    %2.2f    |   %2.4f    |\n', R(nR), num.Zdim, eer_PLS_PLDA(nR,nZ), dcf_vox_PLS_PLDA(nR,nZ));
        fprintf('+------------------+-------------+----------+\n');

    end

    
%     for iTers = 1 : 1
% 
%         if iTers == 1
%             pLDA = gplda_em(W_dev', label_dev(:), num.Zdim, nIters);
%         else       
%             pLDA = gplda_em(W_dev', label_dev(:), num.Zdim, nIters, pLDA);
%         end
%     end

    
% %     fprintf('========== Summaray =========== PLS =====%d====\n',R(nR));
% %     fprintf('+----------------------------------------------------------+\n');
% %     fprintf('|      Method      |    EER(%%)   | Min DCF_vox |    RESS     |\n');
% %     fprintf('+------------------+-------------+----------------------------+\n');
% % %     fprintf('|PLS(%d)+CDS     |    %2.2f    |   %2.4f    |    %2.2f    |\n', R(nR),           eer_PLS_CDS(nR),  dcf_vox_PLS_CDS(nR),  S_RESS(nR));
% %     fprintf('|PLS(%d)+PLDA(%d)|    %2.2f    |   %2.4f    |    %2.2f    |\n', R(nR), num.Zdim, eer_PLS_PLDA(nR), dcf_vox_PLS_PLDA(nR), S_RESS(nR));
% %     fprintf('+------------------+-------------+----------------------------+\n');

end