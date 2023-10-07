% 提取ubm模型 --mfcc特征 -- 69维 -- 1024个高斯分量
trainMfcc = importdata('/mnt/g813_u7/kaldi-mfcc/olr2020/mfcc_pitch/trainMfcc.mat');  % 94282*1
nmix = 1024;  % 特征
final_niter = 10;  % 迭代次数
ds_factor = 1; % 二次采样因子，默认1
nWorkers = 12; % cpu核心数 -- 12
% ubm -- {w:1*1024; mu:120*1024; sigma:120*1024; }
ubm = gmm_em(trainMfcc(:),nmix, final_niter, ds_factor, nWorkers);
save('/mnt/g813_u7/kaldi-mfcc/olr2020/data_mfcc/ubm_mfcc_69.mat','ubm','-v7.3');

% ==========================================================================
% 提取ubm模型 --sdc特征 -- 138维 -- 1024个高斯分量
nmix = 1024;  % 特征
final_niter = 10;  % 迭代次数
ds_factor = 1; % 二次采样因子，默认1
nWorkers = 12; % cpu核心数 -- 12
% ubm -- {w:1*1024; mu:120*1024; sigma:120*1024; }
ubm = gmm_em(enrollMfcc(:),nmix, final_niter, ds_factor, nWorkers);
save('/mnt/g813_u7/kaldi-mfcc/olr2020/data_sdc/ubm_sdc_138.mat','ubm','-v7.3');


% ==========================================================================
% 提取ubm模型 --sdc特征 -- 138维 -- 512个高斯分量
trainMfcc = importdata('/mnt/g813_u7/kaldi-mfcc/olr2020/data_sdc/enrollMfcc.mat');  % 94282*1
nmix = 512;  % 特征
final_niter = 10;  % 迭代次数
ds_factor = 1; % 二次采样因子，默认1
nWorkers = 12; % cpu核心数 -- 12
% ubm -- {w:1*1024; mu:120*1024; sigma:120*1024; }
ubm = gmm_em(trainMfcc,nmix, final_niter, ds_factor, nWorkers);
save('/mnt/g813_u7/kaldi-mfcc/olr2020/data_sdc/512/ubm_sdc_138_512.mat','ubm','-v7.3');



%% ============================task2==================================
% 提取ubm模型 --sdc特征 -- 138维 -- 512个高斯分量
trainMfcc = importdata('/mnt/g813_u7/kaldi-mfcc/olr2020_task2/data_sdc/enrollMfcc.mat');  % 94282*1
nmix = 512;  % 特征
final_niter = 10;  % 迭代次数
ds_factor = 1; % 二次采样因子，默认1
nWorkers = 12; % cpu核心数 -- 12
% ubm -- {w:1*1024; mu:120*1024; sigma:120*1024; }
ubm = gmm_em(trainMfcc,nmix, final_niter, ds_factor, nWorkers);
save('/mnt/g813_u7/kaldi-mfcc/olr2020_task2/data_sdc/ubm_sdc_138_512.mat','ubm','-v7.3');



