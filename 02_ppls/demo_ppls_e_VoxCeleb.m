%% Loading data

M_dev = importdata('/data/chenchen/matlab/program/ASR_Toolbox/00_Data/VoxCeleb1/Mix_60_1024/M_dev.mat');
M_eva = importdata('/data/chenchen/matlab/program/ASR_Toolbox/00_Data/VoxCeleb1/Mix_60_1024/M_eva.mat');

% % spk_cell_dev = importdata('/data/chenchen/matlab/program/ASR_Toolbox/00_Data/VoxCeleb1/Label/SV/dev_speaker_v_cell.mat');
label_dev = importdata('/data/chenchen/matlab/program/ASR_Toolbox/00_Data/VoxCeleb1/Mix_60_1024/label_dev.mat');
Y_dev = full(ind2vec(label_dev));

% m_dev = zscore(M_dev',1)';
% % % y_dev = zscore(Y_dev',1)';
% m_eva = zscore(M_eva',1)';

num.IVdim = 400;
num.Zdim = 200;

% label_dev_num = [];
% % label_dev_weight = [];
% for nspk = 1 : size(spk_cell_dev,1)
%     label_dev_num(nspk) = size(spk_cell_dev{nspk},2);s
% %     label_dev_weight = [label_dev_weight label_dev_num(nspk)*ones(1,label_dev_num(nspk))];
% end
% clear nspk
% 
% % y_dev_weight = bsxfun(@rdivide, Y_dev, label_dev_weight);


% num_choose_weight = 1/2;
% label_dev_session = floor(num_choose_weight*label_dev_num);
% 
% ind = [];ind_add = [];
% for nspk = 1 : size(spk_cell_dev,1)
% 
%     spk_session = randperm(label_dev_num(nspk),label_dev_session(nspk));%[1:1:num_session_A];
% 
%     spk_ind = sum(label_dev_num(1:nspk-1));
%     ind_add = spk_ind+spk_session;
%     ind = [ind ind_add];
% end
% 
% m_dev_choose = m_dev(:, ind);
% Y_dev_choose = Y_dev(:, ind);
% label_dev_choose = label_dev(ind);

%% ---------------------------------------------Train the PPLS-based TVS

num.iters = 15;
for iTers = 1 : 20

    if iTers == 1
        [patameters_PPLS] = ppls_em(M_dev,Y_dev,num);
    else
        [patameters_PPLS] = ppls_em(M_dev,Y_dev,num,patameters_PPLS);
    end




%% evaluation
IV_dev   = pinv(patameters_PPLS.L)*patameters_PPLS.B*(M_dev-patameters_PPLS.m);
IV_eva   = pinv(patameters_PPLS.L)*patameters_PPLS.B*(M_eva-patameters_PPLS.m);

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



fprintf('================= cosine ====================\n');
scores_PPLS_CDS_vox = 1 - pdist2(target2_IVs',target1_IVs','cosine');
scores_PPLS_CDS_vox = diag(scores_PPLS_CDS_vox);
[eer_PPLS400_CDS(iTers),~,~,dcf_vox_PPLS400_CDS(iTers)]=compute_eer(scores_PPLS_CDS_vox,answer_eva,0);
% hold on


% fprintf('================= PLDA ======================\n');
% nIters = 50;
% % for iTers_PLDA = 1 : 1
% %     
% %     if iTers_PLDA == 1
%         pLDA = gplda_em(IV_dev, label_dev(:), num.Zdim, nIters);
% %     else       
% %         pLDA = gplda_em(IV_dev, label_dev(:), num.Zdim, nIters, pLDA);
% %     end
% % end
%     
% scores_PPLS_PLDA=score_gplda_trials(pLDA,target1_IVs,target2_IVs);
% scores_PPLS_PLDA=diag(scores_PPLS_PLDA);
% [eer_PPLS400_PLDA200(iTers),~,~,dcf_vox_PPLS400_PLDA200(iTers)]=compute_eer(scores_PPLS_PLDA,answer_eva,0);


fprintf('===== Summaray ========== PPLS =====%d===========\n',iTers);
fprintf('+-----------------------------------------------+\n');
fprintf('|       Method      |    EER(%%)   | Min DCF_vox |\n');
fprintf('+-------------------+-------------+-------------+\n');
fprintf('|  PPLS(%d)+CDS    |    %2.2f    |   %2.4f    |\n', num.IVdim,         eer_PPLS400_CDS(iTers),    dcf_vox_PPLS400_CDS(iTers));
% fprintf('|PPLS(%d)+PLDA(%d)|     %2.2f    |   %2.4f    |\n', num.IVdim,num.Zdim,eer_PPLS400_PLDA200(iTers),dcf_vox_PPLS400_PLDA200(iTers));
fprintf('+-----------------+-------------+-------------+\n');

end
