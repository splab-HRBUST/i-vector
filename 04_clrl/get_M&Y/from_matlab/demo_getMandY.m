% ==== 0.加载音频文件 ====
file_path = getfilepath('../demo/audio/');
save('../demo/file_path.mat','file_path');

% ==== 1.mfcc特征提取 ====
file_path = importdata('../demo/file_path.mat');
file_num = size(file_path,1); 
testMFCC = cell(file_num,1);
for i = 1:file_num
    voice_path = file_path{i,1} % 音频路径
    mfcc = extract_mfcc(voice_path,'0');  % 提取mfcc {1*1cell-->60*433}
    testMFCC{i,1} = mfcc{1,1};
end
save('../demo/testMFCC.mat','testMFCC','-v7.3');  % 保存

% ==== 2.训练ubm模型 ====
nmix = 1024;  % 特征
final_niter = 10;  % 迭代次数
ds_factor = 1;
nWorkers = 1;
% ubm -- {w:1*1024; mu:60*1024; sigma:60*1024; }
ubm = gmm_em(testMFCC(:),nmix, final_niter, ds_factor, nWorkers);
save('../demo/ubm.mat','ubm','-v7.3');

% ==== 3.根据ubm得到M ====
ubm_path = '../demo/ubm.mat';
M_test = getM(file_path,ubm_path);
save('../demo/M_test.mat','M_test','-v7.3');

% ==== 4.得到onehot标签Y ====
label_path= '../demo/test_label.txt';
Y_test = getY(label_path);
save('../demo/Y_test.mat','Y_test','-v7.3');








