%% load IVs on King_ASR-010
% % M_all = importdata('00_Data/King-ASR-010/Mix_60_1024/M_all.mat');
% Ex_all = pinv(patameters_PPCA_baseline.L)*patameters_PPCA_baseline.B*(M_all-patameters_PPCA_baseline.m);
% % 
% for nSpk = 1 : 200
%     Ex_enroll(:,(nSpk-1)*num.enroll+1:nSpk*num.enroll) =  Ex_all(:,(nSpk-1)*120+1:(nSpk-1)*120+96);
%     Ex_test(:,(nSpk-1)*num.test+1  :nSpk*num.test)     =  Ex_all(:,(nSpk-1)*120+97:nSpk*120);
% end
% clear nSpk
% 
% Ex_model = zeros(num.IVdim,200);
% for i = 1 : 200
%     Ex_model(:,i) = mean(Ex_enroll(:,(i-1)*num.enroll+1:i*num.enroll),2); % i-vector models
% end
% clear i

%
num.nSpk = 200;
% num.IVdim = 400;
num.enroll = 96;
num.test = 24;
% 
num.label_enroll = repmat([1:1:200],96,1); 
num.label_test   = repmat([1:1:200],24,1); 
num.label_enroll = num.label_enroll(:); 
num.label_test   = num.label_test(:); 

%% No session compensation

% answer_eva = [ones(1,200*num.test) zeros(1,200*num.test*199)];

% fprintf('================= IV + CDS ======================\n');
% scores_CDS = [];
% scores_CDS.all = 1 - pdist2(Ex_test',Ex_model','cosine');
% ACC_CDS  = top_N_acc(scores_CDS.all ,num.label_test,5)
% 
% scores_CDS.true = [];
% scores_CDS.impostor = [];
% 
% for a =  1 : 200
%     for b = 1 : 200
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
% [eer_baseline_FA400_CDS,~,~,~,dcf_king_baseline_FA400_CDS]=compute_eer(scores_CDS,answer_eva,false)


% fprintf('================= IV + PLDA ======================\n');
% nphi = 200;%ldaDim; 50，100，125，150，200            % should be <= ldaDim
% niter = 100;
% pLDA = gplda_em(Ex_enroll, num.label_enroll(:), nphi, niter);
% 
% scores_PLDA = [];
% scores_PLDA.all = (score_gplda_trials(pLDA, Ex_model, Ex_test))';
% ACC_PLDA = top_N_acc(scores_PLDA.all,num.label_test,5)
% 
% scores_PLDA.true = [];
% scores_PLDA.impostor = [];
% 
% for a =  1 : 200
%     for b = 1 : 200
%         A = scores_PLDA.all(num.test*(a-1)+1 : num.test*(a-1)+num.test , b);
%         if a == b
%             scores_PLDA.true = [scores_PLDA.true ; A(:)];
%         elseif a~=b
%             scores_PLDA.impostor = [scores_PLDA.impostor ; A(:)];
%         end
%     end
% end  
% clear a b A
% scores_PLDA = [scores_PLDA.true;scores_PLDA.impostor];
% [eer_baseline_FA400_PLDA200,~,~,~,dcf_king_baseline_FA400_PLDA200]=compute_eer(scores_PLDA,answer_eva,false)


% fprintf('================= IV + softmax ======================\n');
% [theta,~] = mysoftmax_gd(Ex_enroll',num.label_enroll(:),1e-8,10,100);    
% scores_softmax.all = exp(Ex_test'*theta)./repmat(sum(exp(Ex_test'*theta),2),1,num.nSpk);
% ACC_softmax = top_N_acc(scores_softmax.all,num.label_test,5)
% 
% scores_softmax.true = [];
% scores_softmax.impostor = [];
% 
% for a =  1 : 200
%     for b = 1 : 200
%         A = scores_softmax.all(num.test*(a-1)+1 : num.test*(a-1)+num.test , b);
%         if a == b
%             scores_softmax.true = [scores_softmax.true ; A(:)];
%         elseif a~=b
%             scores_softmax.impostor = [scores_softmax.impostor ; A(:)];
%         end
%     end
% end  
% clear a b A
% scores_softmax = [scores_softmax.true;scores_softmax.impostor];
% [eer_baseline_FA400_softmax,~,~,~,dcf_king_baseline_FA400_softmax]=compute_eer(scores_softmax,answer_eva,1)



%  fprintf('================= IV + SVM ======================\n');
% [~, ~, svmOptTheta] = multisvmtrain(Ex_enroll, num.label_enroll, num.nSpk, 1e-8, 0.01, 200, svmOptTheta);
% scores_SVM = [];
% scores_SVM.all = (svmOptTheta * Ex_test)';
% ACC_SVM = top_N_acc(scores_SVM.all, num.label_test, 5)
% 
% scores_SVM.true = [];
% scores_SVM.impostor = [];
% 
% for a =  1 : 200
%     for b = 1 : 200
%         A = scores_SVM.all(num.test*(a-1)+1 : num.test*(a-1)+num.test , b);
%         if a == b
%             scores_SVM.true = [scores_SVM.true ; A(:)];
%         elseif a~=b
%             scores_SVM.impostor = [scores_SVM.impostor ; A(:)];
%         end
%     end
% end  
% clear a b A
% scores_SVM = [scores_SVM.true;scores_SVM.impostor];
% [eer_baseline_FA400_SVM,~,~,~,dcf_king_baseline_FA400_SVM]=compute_eer(scores_SVM,answer_eva,1)
    


%% LDA    
% fprintf('================= IV + LDA + CDS ======================\n');
% ldaDim = 150
% [V,~] = lda(Ex_enroll, num.label_enroll(:));
% final_enrollIV = V(:, 1:ldaDim)' * Ex_enroll;
% final_testIV = V(:, 1:ldaDim)' * Ex_test;
% final_modelIV = V(:, 1:ldaDim)' * Ex_model;
% 
% scores_LDA_CDS = [];
% scores_LDA_CDS.all = 1 - pdist2(final_testIV',final_modelIV','cosine');
% ACC_LDA_CDS = top_N_acc(scores_LDA_CDS.all,num.label_test, 5)
% 
% scores_LDA_CDS.true = [];
% scores_LDA_CDS.impostor = [];
% 
% for a =  1 : 200
%     for b = 1 : 200
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
% [eer_baseline_FA400_LDA_CDS,~,~,~,dcf_king_baseline_FA400_LDA_CDS]=compute_eer(scores_LDA_CDS,answer_eva,1)



% fprintf('================= IV + LDA + PLDA ======================\n');
% nphi = 120            
% niter = 100
% pLDA_lda = gplda_em(final_enrollIV, num.label_enroll(:), nphi, niter);
% 
% scores_LDA_PLDA = [];
% scores_LDA_PLDA.all = (score_gplda_trials(pLDA_lda, final_modelIV, final_testIV))';
% ACC_LDA_PLDA = top_N_acc(scores_LDA_PLDA.all,num.label_test,5)
% 
% scores_LDA_PLDA.true = [];
% scores_LDA_PLDA.impostor = [];
% 
% for a =  1 : 200
%     for b = 1 : 200
%         A = scores_LDA_PLDA.all(num.test*(a-1)+1 : num.test*(a-1)+num.test , b);
%         if a == b
%             scores_LDA_PLDA.true = [scores_LDA_PLDA.true ; A(:)];
%         elseif a~=b
%             scores_LDA_PLDA.impostor = [scores_LDA_PLDA.impostor ; A(:)];
%         end
%     end
% end  
% clear a b A
% scores_LDA_PLDA = [scores_LDA_PLDA.true;scores_LDA_PLDA.impostor];
% [eer_baseline_FA400_LDA_PLDA,~,~,~,dcf_king_baseline_FA400_LDA_PLDA]=compute_eer(scores_LDA_PLDA,answer_eva,1)



% fprintf('================= IV + LDA + softmax ======================\n');
% [theta,~] = mysoftmax_gd(final_enrollIV',num.label_enroll(:),1e-8,1e4,300); %SGD
% scores_LDA_softmax = [];
% scores_LDA_softmax.all = exp(final_testIV'*theta)./repmat(sum(exp(final_testIV'*theta),2),1,num.nSpk);
% ACC_LDA_softmax = top_N_acc(scores_LDA_softmax.all,num.label_test,5)
% 
% scores_LDA_softmax.true = [];
% scores_LDA_softmax.impostor = [];
% 
% for a =  1 : 200
%     for b = 1 : 200
%         A = scores_LDA_softmax.all(num.test*(a-1)+1 : num.test*(a-1)+num.test , b);
%         if a == b
%             scores_LDA_softmax.true = [scores_LDA_softmax.true ; A(:)];
%         elseif a~=b
%             scores_LDA_softmax.impostor = [scores_LDA_softmax.impostor ; A(:)];
%         end
%     end
% end  
% clear a b A
% scores_LDA_softmax = [scores_LDA_softmax.true;scores_LDA_softmax.impostor];
% [eer_baseline_FA400_LDA_softmax,~,~,~,dcf_king_baseline_FA400_LDA_softmax]=compute_eer(scores_LDA_softmax,answer_eva,1)



% fprintf('================= IV + LDA + SVM ======================\n');
% [~, ~, svmOptTheta] = multisvmtrain(final_enrollIV, num.label_enroll, num.nSpk, 1e-8, 1e3, 50); %SGD
% scores_LDA_SVM = [];
% scores_LDA_SVM.all = (svmOptTheta * final_testIV)';
% ACC_LDA_SVM = top_N_acc(scores_LDA_SVM.all, num.label_test, 5)

% scores_LDA_SVM.true = [];
% scores_LDA_SVM.impostor = [];
% 
% for a =  1 : 200
%     for b = 1 : 200
%         A = scores_LDA_SVM.all(num.test*(a-1)+1 : num.test*(a-1)+num.test , b);
%         if a == b
%             scores_LDA_SVM.true = [scores_LDA_SVM.true ; A(:)];
%         elseif a~=b
%             scores_LDA_SVM.impostor = [scores_LDA_SVM.impostor ; A(:)];
%         end
%     end
% end  
% clear a b A
% scores_LDA_SVM = [scores_LDA_SVM.true;scores_LDA_SVM.impostor];
% [eer_FA400_LDA_SVM,~,~,~,dcf_king_FA400_LDA_SVM]=compute_eer(scores_LDA_SVM,answer_eva,1)
 
 


%% NDA
% fprintf('================= IV + NDA + CDS ======================\n');   
% train.mat = Ex_enroll;
% train.labels = num.label_enroll(:);
% train.N = size(train.mat, 2);
% train.dim = size(train.mat, 1);
% train.totclass=length(unique(train.labels));
% 
% % nearest neighbor
% K = 5;
% % dim reduction to dim
% ndaDim = 200;
% % use weight to control boundary
% useweights = 1;
% [proymat,ndaModel]=nda(train,K,ndaDim,useweights);
% 
% %fea x num
% final_enrollIV_nda = ndaModel.mat;
% final_testIV_nda  = ndaModel.W * (Ex_test-repmat(ndaModel.meandata,1,size(Ex_test,2)));
% final_modelIV_nda = ndaModel.W * (Ex_model-repmat(ndaModel.meandata,1,size(Ex_model,2)));
% 
% scores_NDA_CDS = [];
% scores_NDA_CDS.all = 1 - pdist2(final_testIV_nda',final_modelIV_nda','cosine');
% ACC_NDA_CDS = top_N_acc(scores_NDA_CDS.all,num.label_test,5)

% scores_NDA_CDS.true = [];
% scores_NDA_CDS.impostor = [];
% 
% for a =  1 : 200
%     for b = 1 : 200
%         A = scores_NDA_CDS.all(num.test*(a-1)+1 : num.test*(a-1)+num.test , b);
%         if a == b
%             scores_NDA_CDS.true = [scores_NDA_CDS.true ; A(:)];
%         elseif a~=b
%             scores_NDA_CDS.impostor = [scores_NDA_CDS.impostor ; A(:)];
%         end
%     end
% end  
% clear a b A
% scores_NDA_CDS = [scores_NDA_CDS.true;scores_NDA_CDS.impostor];
% [eer_FA400_NDA_CDS,~,~,~,dcf_king_FA400_NDA_CDS]=compute_eer(scores_NDA_CDS,answer_eva,1)



% fprintf('=================IV + NDA + PLDA ======================\n');
% nphi = 120
% niter = 100
% pLDA_nda = gplda_em(final_enrollIV_nda, num.label_enroll(:), nphi, niter);
% 
% scores_NDA_PLDA = [];
% scores_NDA_PLDA.all = (score_gplda_trials(pLDA_nda, final_modelIV_nda, final_testIV_nda))';
% ACC_NDA_PLDA = top_N_acc(scores_NDA_PLDA.all,num.label_test,5)

% scores_NDA_PLDA.true = [];
% scores_NDA_PLDA.impostor = [];
% 
% for a =  1 : 200
%     for b = 1 : 200
%         A = scores_NDA_PLDA.all(num.test*(a-1)+1 : num.test*(a-1)+num.test , b);
%         if a == b
%             scores_NDA_PLDA.true = [scores_NDA_PLDA.true ; A(:)];
%         elseif a~=b
%             scores_NDA_PLDA.impostor = [scores_NDA_PLDA.impostor ; A(:)];
%         end
%     end
% end  
% clear a b A
% scores_NDA_PLDA = [scores_NDA_PLDA.true;scores_NDA_PLDA.impostor];
% [eer_FA400_NDA_PLDA,~,~,~,dcf_king_FA400_NDA_PLDA]=compute_eer(scores_NDA_PLDA,answer_eva,1)


% fprintf('================= IV + NDA + softmax ======================\n');
% [theta,~] = mysoftmax_gd(final_enrollIV_nda',num.label_enroll(:),1e-8,1,150);  %SGD
% 
% scores_NDA_softmax = [];
% scores_NDA_softmax.all = exp(final_testIV_nda'*theta)./repmat(sum(exp(final_testIV_nda'*theta),2),1,num.nSpk);
% ACC_NDA_softmax = top_N_acc(scores_NDA_softmax.all,num.label_test,5)




% fprintf('================= IV + NDA + SVM ======================\n');
% [~, ~, svmOptTheta] = multisvmtrain(final_enrollIV_nda, num.label_enroll, num.nSpk, 1e-8, 0.01, 100); %SGD
% scores_NDA_SVM = [];
% scores_NDA_SVM.all = svmOptTheta * final_testIV_nda;
% ACC_NDA_SVM = top_N_acc(scores_NDA_SVM.all', num.label_test, 5)
      



%% SC
fprintf('================= train the UnSup Dictionary ======================\n');

% load('/data/chenchen/matlab/program/Speaker_Verification_Toolbox_v1.c/data/King_ASR_010/Mix_60_1024/IVs/IVs_ExByMyself.mat')
% Ex_all = double(cell2mat(IVs'));

% for nspk = 1 : 200
%     Ex_enroll(:,(nspk-1)*num.enroll+1:nspk*num.enroll) =  Ex_all(:,(nspk-1)*120+1:(nspk-1)*120+96);
%     Ex_test(:,(nspk-1)*num.test+1  :nspk*num.test)     =  Ex_all(:,(nspk-1)*120+97:nspk*120);
% end
% clear nspk

% X_enroll = length_norm(Ex_enroll);
% X_test = length_norm(Ex_test);
% 
dic_size = 1000;
opts.lambda = 0.1; 
opts.lambda2 = 0.1; 
opts.iter = 5;  
opts.computeCost = 1; 
opts.batchSize = 256;
opts.ro = 0.1; % learning rate for D

Y.outputVectorTrain = ind2vec(num.label_enroll');
Y.outputVectorTrain = Y.outputVectorTrain(:,1:end);
Y.trls = num.label_enroll;

[D_UnSup,result_UnSup] = UnSupODL(X_enroll(:,1:end), Y, dic_size, opts); 

Alpha_enroll = full(mexLasso(X_enroll,D_UnSup,opts));
Alpha_test = full(mexLasso(X_test,D_UnSup,opts));

[~, ~, svmOptTheta] = multisvmtrain(Alpha_enroll, num.label_enroll, num.nSpk, 1e-8, 1e4, 100); %SGD

scores_SC_SVM = [];
scores_SC_SVM.all = svmOptTheta * Alpha_test;
ACC_SC_SVM = top_N_acc(scores_SC_SVM.all', num.label_test, 5)
