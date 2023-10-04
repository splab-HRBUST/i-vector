%% loading data
% IV_train = importdata('/data/chenchen/data/new_voxceleb_for_identification/IVs/train/train_IVs.mat');
% train_spk_cell = importdata('/data/chenchen/data/new_voxceleb_for_identification/IVs/train_speaker_v_cell.mat');
% 
% IV_val = importdata('/data/chenchen/data/new_voxceleb_for_identification/IVs/validate/validate_IVs.mat');
% val_spk_cell = importdata('/data/chenchen/data/new_voxceleb_for_identification/IVs/validate_speaker_v_cell.mat');
% 
% IV_test = importdata('/data/chenchen/data/voxceleb/IVs_Identification/test/test_IVs.mat');
% test_spk_cell = importdata('/data/chenchen/data/voxceleb/IVs_Identification/test_speaker_v_cell.mat');
% 
% 
% IV_train = double(IV_train);
% IV_val = double(IV_val);
% IV_test = double(IV_test);
% 
% 
% nspk = size(train_spk_cell,1);%20
% IV_model = [];label_train = [];label_val = [];
% label_test = [];
% spk_id = 1;count = 0;
% 
% for i = 1 : nspk
%     spk_train_num = size(train_spk_cell{i},2);
%     spk_val_num = size(val_spk_cell{i},2)-spk_train_num;
%     spk_test_num = size(test_spk_cell{i},2);
%     
%     IV_model(:,i) = mean(IV_train(:,count+1:count+spk_train_num),2);
%     
%     label_train = [label_train spk_id*ones(1,spk_train_num)];
%     label_val = [label_val spk_id*ones(1,spk_val_num)];
%     label_test = [label_test spk_id*ones(1,spk_test_num)];
% 
%     spk_id = spk_id + 1;
%     count = count + spk_train_num;
% end 
% clearvars spk_id count i spk_train_num spk_val_num spk_test_num

%% speaker identification
% fprintf('================= i-vector+cosine ======================\n');
% cos =1 - pdist2(IV_test(:,1:size(label_test,2))',IV_model','cosine');
% 
% [~,prelabel_test_cos] = max(cos',[],1);
% % acc_iv_cos = sum(prelabel_test_cos == label_test(:)')/size(label_test(:),1);
% acc_iv_cos_top5 = topNError(cos,label_test,5);
% 
% fprintf('================= i-vector+plda ======================\n');
% nphi = 200;                  % should be <= ldaDim
% niter = 30;
% pLDA = gplda_em(IV_train(:,1:size(label_train,2)), label_train(:), nphi, niter);
% 
% 
% ivScores = (score_gplda_trials(pLDA, IV_model, IV_test(:,1:size(label_test,2))))';
% [~,prelabel_test_plda] = max(ivScores',[],1);
% % acc_iv_plda = sum(prelabel_test_plda == label_test(:)')/size(label_test(:),1);
% acc_iv_plda_top5 = topNError(ivScores,label_test,5)

fprintf('================= i-vector+softmax ======================\n');
M = mean(IV_train, 2);
data_train = bsxfun(@minus, IV_train, M); % centering thse data
data_train = length_norm(data_train); % normalizing the length

data_test = bsxfun(@minus, IV_test, M); % centering the data
data_test = length_norm(data_test); % normalizing the length 

% [theta,~] = mysoftmax_gd(data_train',label_train',1e-8,1e1,100); 
[theta,~] = mysoftmax_gd(data_train',label_train',1e-8,1e1,100, theta); 

score_softmax = exp(data_test'*theta)./repmat(sum(exp(data_test'*theta),2),1,nspk);
acc_iv_softmax_top5 = topNError(score_softmax,label_test,5)

% fprintf('================= i-vector+SVM ======================\n');
% M = mean(IV_train, 2);
% data_train = bsxfun(@minus, IV_train, M); % centering thse data
% data_train = length_norm(data_train); % normalizing the length
% 
% data_test = bsxfun(@minus, IV_test, M); % centering the data
% data_test = length_norm(data_test); % normalizing the length 
% 
% [~, ~, svmOptTheta] = multisvmtrain(nspk, size(data_train,1), 1e-8, data_train, label_train, 100, 0.01);
% score_SVM =svmOptTheta * data_test;
% acc_iv_SVM_top5 = topNError(score_SVM', label_test, 5)

%% LDA
% fprintf('================= i-vector+LDA+CDS ======================\n');
% ldaDim = min(200,nspk-1);
% [V,~] = lda(IV_train(:,1:size(label_train,2)), label_train(:));
% final_devIV = V(:, 1:ldaDim)' * IV_train(:,1:size(label_train,2));
% final_testIV = V(:, 1:ldaDim)' * IV_test(:,1:size(label_test,2));
% final_modelIV = V(:, 1:ldaDim)' * IV_model;
% 
% cos_lda =1 - pdist2(final_testIV',final_modelIV','cosine');
% % [~,prelabel_test_lda_cos] = max(cos_lda',[],1);
% acc_iv_lda_cos_top5 = topNError(cos_lda,label_test,5)

% fprintf('=================i-vector+LDA+PLDA ======================\n');
% nphi = ldaDim;             % should be <= ldaDim
% niter = 15;
% pLDA_lda = gplda_em(final_devIV, label_train(:), nphi, niter);
% 
% ivScores_plda = (score_gplda_trials(pLDA_lda, final_modelIV, final_testIV))';
% [~,prelabel_test_lda_plda] = max(ivScores_plda',[],1);
% acc_iv_lda_plda_top5 = topNError(ivScores_plda,label_test,5);

% fprintf('================= i-vector+LDA+softmax ======================\n');
% f_devIV = zscore(final_devIV')';
% % f_devIV = final_devIV;
% [theta,~] = mysoftmax_gd(f_devIV(1:end,:)',label_train(1:end)',1e-8,1e2,100);    
% score_lda_softmax = exp(final_testIV'*theta)./repmat(sum(exp(final_testIV'*theta),2),1,nspk);
% [~,prelabel_test_lda_softmax] = max(score_lda_softmax',[],1);
% acc_iv_lda_softmax_top5 = topNError(score_lda_softmax,label_test,5);

%  fprintf('================= i-vector+LDA+SVM ======================\n');
%  M = mean(final_devIV, 2);
%  data_train = bsxfun(@minus, final_devIV, M); % centering thse data
%  data_train = length_norm(data_train); % normalizing the length
% 
%  data_test = bsxfun(@minus, final_testIV, M); % centering the data
%  data_test = length_norm(data_test); % normalizing the length 
%  
%  [~, ~, svmOptTheta] = multisvmtrain(nspk, size(data_train,1), 1e-8, data_train, label_train, 10, 0.01);
%  score_lda_SVM =svmOptTheta * data_test;
%  acc_iv_lda_SVM_top5 = topNError(score_lda_SVM', label_test, 5);

%% NDA
% fprintf('================= i-vector+NDA+CDS ======================\n');
% train.mat = IV_train;
% train.labels = label_train(:);
% train.N = size(train.mat, 2);
% train.dim = size(train.mat, 1);
% train.totclass=length(unique(train.labels));
% 
% % nearest neighbor
% K = 2;
% % dim reduction to dim
% ndaDim = 300;
% % use weight to control boundary
% useweights = 1;
% 
% 
% % % % ������֤
% % % for nTest = 1 : 8
% % %     [proymat,ndaModel{nTest}]=nda(train,K,ndaDim(nTest),useweights);
% % % 
% % %     %fea x num
% % %     final_devIV_nda = ndaModel{nTest}.mat;
% % %     final_testIV_nda = ndaModel{nTest}.W * (IV_test-repmat(ndaModel{nTest}.meandata,1,size(IV_test,2)));
% % %     final_modelIV_nda = ndaModel{nTest}.W * (IV_model-repmat(ndaModel{nTest}.meandata,1,size(IV_model,2)));
% % % 
% % %     cos_nda =1 - pdist2(final_testIV_nda',final_modelIV_nda','cosine');
% % %     [~,prelabel_test_nda_cos] = max(cos_nda',[],1);
% % %     acc_iv_nda_cos_top5{nTest} = topNError(cos_nda,label_test,5);
% % % end
% 
% [proymat,ndaModel]=nda(train,K,ndaDim,useweights);
% 
% %fea x num
% final_devIV_nda = ndaModel.mat;
% final_testIV_nda = ndaModel.W * (IV_test-repmat(ndaModel.meandata,1,size(IV_test,2)));
% final_modelIV_nda = ndaModel.W * (IV_model-repmat(ndaModel.meandata,1,size(IV_model,2)));
% 
% cos_nda =1 - pdist2(final_testIV_nda',final_modelIV_nda','cosine');
% [~,prelabel_test_nda_cos] = max(cos_nda',[],1);
% acc_iv_nda_cos_top5 = topNError(cos_nda,label_test,5);


% fprintf('=================i-vector+NDA+PLDA ======================\n');
% nphi = 200;             % should be <= ldaDim
% niter = 100;
% pLDA_nda = gplda_em(final_devIV_nda, label_train(:), nphi, niter);
% 
% ivScores_nda_plda = (score_gplda_trials(pLDA_nda, final_modelIV_nda, final_testIV_nda))';
% [~,prelabel_test_nda_plda] = max(ivScores_nda_plda',[],1);
% acc_iv_nda_plda_top5 = topNError(ivScores_nda_plda,label_test,5);

% fprintf('================= i-vector+NDA+softmax ======================\n');
% f_devIV = zscore(final_devIV_nda');
% [theta,~] = mysoftmax_gd(f_devIV(1:end,:),label_train(1:end)',1e-8,100,200);    
% score_nda_softmax = exp(final_testIV_nda'*theta)./repmat(sum(exp(final_testIV_nda'*theta),2),1,nspk);
% [~,prelabel_test_nda_softmax] = max(score_nda_softmax',[],1);
% acc_iv_nda_softmax_top5 = topNError(score_nda_softmax,label_test,5);

%  fprintf('================= i-vector+NDA+SVM ======================\n');
%  [~, ~, svmOptTheta] = multisvmtrain(nspk, size(f_devIV,2), 1e-2, f_devIV(1:end,:)', label_train(1:end), 100, 1,svmOptTheta);
%  score_nda_SVM =svmOptTheta * final_testIV_nda;
%  [~, prelabel_test_nda_SVM] = max(score_nda_SVM',[],1);
%  acc_iv_nda_SVM_top5 = topNError(score_nda_SVM', label_test, 5);

% fprintf('================= i-vector+LFDA+cosine ======================\n');
% lfdaDim = min(120,nspk-1);
% [V_lfda,~]=LFDA(IV_train(:,1:size(label_train,2)), label_train(:),lfdaDim,'weighted',7);
% final_devIV_lfda = V_lfda(:, 1:lfdaDim)' * IV_train(:,1:size(label_train,2));
% final_testIV_lfda = V_lfda(:, 1:lfdaDim)' * IV_test(:,1:size(label_test,2));
% final_modelIV_lfda = V_lfda(:, 1:lfdaDim)' * IV_model;
% 
% cos_lfda =1 - pdist2(final_testIV_lfda',final_modelIV_lfda','cosine');
% [~,prelabel_test_lfda_cos] = max(cos_lfda',[],1);
% % acc_iv_lda_cos = sum(prelabel_test_lda_cos == label_test(:)')/size(label_test(:),1);
% acc_iv_lfda_cos_top5 = topNError(cos_lfda,label_test,5);
% 
% fprintf('=================i-vector+LDA+plda ======================\n');
% nphi = lfdaDim;             % should be <= ldaDim
% niter = 10;
% pLDA_lfda = gplda_em(final_devIV_lfda, label_train(:), nphi, niter);
% 
% ivScores_lfda_plda = (score_gplda_trials(pLDA_lfda, final_modelIV_lfda, final_testIV_lfda))';
% [~,prelabel_test_lda_plda] = max(ivScores_lfda_plda',[],1);
% % acc_iv_lda_plda = sum(prelabel_test_lda_plda == label_test(:)')/size(label_test(:),1);
% acc_iv_lfda_plda_top5 = topNError(ivScores_lfda_plda,label_test,5);

% % 
% fprintf('==========================================Summary=====================================\n');
% fprintf('|    Method    |    top-1    |    top-2    |    top-3    |    top-4    |    top-5    |\n');
% fprintf('+--------------+-------------+-------------+-------------+-------------+-------------+');
% % % % fprintf('|    iv+cos    |');
% % % % fprintf('   %2.2f%%    |',100*acc_iv_cos_top5);
% % % % fprintf('\n|   iv+plda    |');
% % % % fprintf('   %2.2f%%    |',100*acc_iv_plda_top5);
% fprintf('\n|  iv+LDA+CDS  |');
% fprintf('   %2.2f%%    |',100*acc_iv_lda_cos_top5);
% fprintf('\n|  iv+LDA+PLDA |');
% fprintf('   %2.2f%%    |',100*acc_iv_lda_plda_top5);
% fprintf('\n|  iv+LDA+softmax  |');
% fprintf('   %2.2f%%    |',100*acc_iv_lda_softmax_top5);
% % fprintf('\n|  iv+LDA+SVM |');
% % fprintf('   %2.2f%%    |',100*acc_iv_lda_SVM_top5);
% % fprintf('\n|  iv+NDA-%d+CDS  |',ndaDim);
% % fprintf('   %2.2f%%    |',100*acc_iv_nda_cos_top5);
% % fprintf('\n|  iv+NDA+PLDA-%d |',nphi);
% % fprintf('   %2.2f%%    |',100*acc_iv_nda_plda_top5);
% % fprintf('\n|  iv+NDA+softmax  |');
% % fprintf('   %2.2f%%    |',100*acc_iv_nda_softmax_top5);
% % fprintf('\n|  iv+NDA+SVM |');
% % fprintf('   %2.2f%%    |',100*acc_iv_nda_SVM_top5);
% % fprintf('\n');
% % fprintf('|   TestSet    |');
% % fprintf('   %2.2f%%    |',100*acc_pre_test_top5);
% fprintf('\n');
% fprintf('+--------------+-------------+-------------+-------------+-------------+-------------+\n');




% % |  TestSet   |   62.37%    |   79.21%    | Dic=1500;L1=0.01,L2=0.005,BS=1024
% % |  TestSet   |   63.24%    |   79.30%    | Dic=2000;L1=0.01,L2=0.005,BS=1024