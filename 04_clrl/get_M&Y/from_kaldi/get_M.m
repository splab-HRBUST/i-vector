% ubm -- 69*1024 = 70656
ubm = importdata('/mnt/g813_u7/kaldi-mfcc/olr2020/data_mfcc/ubm_mfcc_69.mat');

% 训练集 -- 94282
trainMfcc = importdata('/mnt/g813_u7/kaldi-mfcc/olr2020/data_mfcc/trainMfcc.mat');  % 94282*1
train_num = size(trainMfcc,1);
M_train = zeros(70656,train_num); % (70656,94282)
for i = 1:train_num
    i
    mfcc = trainMfcc(i,:);
    gmm = mapAdapt(mfcc,ubm,10.0,'m');
    Mi = gmm.mu(:);
    M_train(:,i) = Mi;
end
save('/mnt/g813_u7/kaldi-mfcc/olr2020/data_mfcc/M_train.mat','M_train','-v7.3');

% 注册集 -- 54269
enrollMfcc = importdata('/mnt/g813_u7/kaldi-mfcc/olr2020/mfcc_pitch/enrollMfcc.mat');  % 94282*1
enroll_num = size(enrollMfcc,1);
M_enroll = zeros(70656,enroll_num);
for i = 1:enroll_num
    i
    mfcc = enrollMfcc(i,:);
    gmm = mapAdapt(mfcc,ubm,10.0,'m');
    Mi = gmm.mu(:);
    M_enroll(:,i) = Mi;
end
save('/mnt/g813_u7/kaldi-mfcc/olr2020/mfcc_pitch/M_enroll.mat','M_enroll','-v7.3');

% 训练集 -- 11848
% testMfcc = importdata('/mnt/g813_u7/kaldi-mfcc/olr2020/mfcc_pitch/testMfcc.mat');  % 94282*1
test_num = size(testMfcc,1);
M_test = zeros(70656,test_num);
for i = 1:test_num
    i
    mfcc = testMfcc(i,:);
    gmm = mapAdapt(mfcc,ubm,10.0,'m');
    Mi = gmm.mu(:);
    M_test(:,i) = Mi;
end
save('/mnt/g813_u7/kaldi-mfcc/olr2020/mfcc_pitch/M_test.mat','M_test','-v7.3');

% 拆分训练集
M_train_1 = zeros(70656,47141);
M_train_2 = zeros(70656,47141);
j=1;
for i = 1:2:94282
    j
    M_train_1(:,j) = M_train(:,i);
    j = j+1;
end
save('/mnt/g813_u7/kaldi-mfcc/olr2020/data_mfcc/M_train_1.mat','M_train_1','-v7.3');
j=1;
for i = 2:2:94282
    j
    M_train_2(:,j) = M_train(:,i);
    j = j+1;
end
save('/mnt/g813_u7/kaldi-mfcc/olr2020/data_mfcc/M_train_2.mat','M_train_2','-v7.3');


%% ===========================SDC-1024==================================================
% ubm -- 138*1024 = 141312
ubm = importdata('/mnt/g813_u7/kaldi-mfcc/olr2020/data_sdc/ubm_sdc_138_1024.mat');

% 注册集 -- 54269
enrollMfcc = importdata('/mnt/g813_u7/kaldi-mfcc/olr2020/data_sdc/enrollMfcc.mat');  % 94282*1
enroll_num = size(enrollMfcc,1);
M_enroll = zeros(141312,enroll_num);
for i = 1:enroll_num
    i
    mfcc = enrollMfcc(i,:);
    gmm = mapAdapt(mfcc,ubm,10.0,'m');
    Mi = gmm.mu(:);
    M_enroll(:,i) = Mi;
end
save('/mnt/g813_u7/kaldi-mfcc/olr2020/data_sdc/M_enroll.mat','M_enroll','-v7.3');

% 测试集 -- 11848
% testMfcc = importdata('/mnt/g813_u7/kaldi-mfcc/olr2020/data_sdc/testMfcc.mat');  % 94282*1
test_num = size(testMfcc,1);
M_test = zeros(141312,test_num);
for i = 1:test_num
    i
    mfcc = testMfcc(i,:);
    gmm = mapAdapt(mfcc,ubm,10.0,'m');
    Mi = gmm.mu(:);
    M_test(:,i) = Mi;
end
save('/mnt/g813_u7/kaldi-mfcc/olr2020/data_sdc/M_test.mat','M_test','-v7.3');

% 训练集（注册集代替）-- 拆分
M_train_1 = zeros(141312,29411);
M_train_2 = zeros(141312,29411);
j=1;
for i = 1:2:58822
    j
    M_train_1(:,j) = M_enroll(:,i);
    j = j+1;
end
save('/mnt/g813_u7/kaldi-mfcc/olr2020/data_sdc/M_train_1.mat','M_train_1','-v7.3');
j=1;
for i = 2:2:58822
    j
    M_train_2(:,j) = M_enroll(:,i);
    j = j+1;
end
save('/mnt/g813_u7/kaldi-mfcc/olr2020/data_sdc/M_train_2.mat','M_train_2','-v7.3');

Y_train = importdata('/mnt/g813_u7/kaldi-mfcc/olr2020/data_mfcc/Y_enroll.mat');
Y_train_1 = zeros(6,29411);
Y_train_2 = zeros(6,29411);
j=1;
for i = 1:2:58822
    j
    Y_train_1(:,j) = Y_train(:,i);
    j = j+1;
end
save('/mnt/g813_u7/kaldi-mfcc/olr2020/data_sdc/Y_train_1.mat','Y_train_1','-v7.3');
j=1;
for i = 2:2:58822
    j
    Y_train_2(:,j) = Y_train(:,i);
    j = j+1;
end
save('/mnt/g813_u7/kaldi-mfcc/olr2020/data_sdc/Y_train_2.mat','Y_train_2','-v7.3');

b_train = importdata('/mnt/g813_u7/kaldi-mfcc/olr2020/mfcc_pitch/y_enroll.txt');
b_train_1 = zeros(29411,1);
b_train_2 = zeros(29411,1);
j=1;
for i = 1:2:58822
    j
    b_train_1(j,1) = b_train(i,1);
    j = j+1;
end
save('/mnt/g813_u7/kaldi-mfcc/olr2020/data_sdc/b_train_1.mat','b_train_1','-v7.3');
j=1;
for i = 2:2:58822
    j
    b_train_2(j,1) = b_train(i,1);
    j = j+1;
end
save('/mnt/g813_u7/kaldi-mfcc/olr2020/data_sdc/b_train_2.mat','b_train_2','-v7.3');


%% ===========================SDC-512==================================================
% ubm -- 138*512 = 70656
ubm = importdata('/mnt/g813_u7/kaldi-mfcc/olr2020/data_sdc/512/ubm_sdc_138_512.mat');

% 训练集 -- 54269
trainMfcc = importdata('/mnt/g813_u7/kaldi-mfcc/olr2020/data_sdc/trainMfcc.mat');  % 94282*1
train_num = size(trainMfcc,1);
M_train = zeros(70656,train_num);
for i = 1:train_num
    i
    mfcc = trainMfcc(i,:);
    gmm = mapAdapt(mfcc,ubm,10.0,'m');
    Mi = gmm.mu(:);
    M_train(:,i) = Mi;
end
save('/mnt/g813_u7/kaldi-mfcc/olr2020/data_sdc/512/M_train.mat','M_train','-v7.3');


% 注册集 -- 54269
enrollMfcc = importdata('/mnt/g813_u7/kaldi-mfcc/olr2020/data_sdc/enrollMfcc.mat');  % 94282*1
enroll_num = size(enrollMfcc,1);
M_enroll = zeros(70656,enroll_num);
for i = 1:enroll_num
    i
    mfcc = enrollMfcc(i,:);
    gmm = mapAdapt(mfcc,ubm,10.0,'m');
    Mi = gmm.mu(:);
    M_enroll(:,i) = Mi;
end
save('/mnt/g813_u7/kaldi-mfcc/olr2020/data_sdc/512/M_enroll.mat','M_enroll','-v7.3');

% 测试集 -- 11848
testMfcc = importdata('/mnt/g813_u7/kaldi-mfcc/olr2020/data_sdc/testMfcc.mat');  % 94282*1
test_num = size(testMfcc,1);
M_test = zeros(70656,test_num);
for i = 1:test_num
    i
    mfcc = testMfcc(i,:);
    gmm = mapAdapt(mfcc,ubm,10.0,'m');
    Mi = gmm.mu(:);
    M_test(:,i) = Mi;
end
save('/mnt/g813_u7/kaldi-mfcc/olr2020/data_sdc/512/M_test.mat','M_test','-v7.3');


%% ==============================task_3(SDC)=============================================
% ubm -- 138*512 = 70656
ubm = importdata('/mnt/g813_u7/kaldi-mfcc/olr2020/data_sdc/512/ubm_sdc_138_512.mat');

% task3的训练集 -- 5种语言 -- 48603条语音
M = [];
M = [M,M_train(:,25863:36094),M_train(:,46319:75092),M_train(:,84686:94282)];  % 70656*48603
save('/mnt/g813_u7/kaldi-mfcc/olr2020_task3/data_sdc/M_train.mat','M','-v7.3');

% task3的训练集 -- 5种语言 -- 9496条语音
testMfcc = importdata('/mnt/g813_u7/kaldi-mfcc/olr2020_task3/data_sdc/testMfcc.mat');  % 94282*1
test_num = size(testMfcc,1);
M_test = zeros(70656,test_num);
for i = 1:test_num
    i
    mfcc = testMfcc(i,:);
    gmm = mapAdapt(mfcc,ubm,10.0,'m');
    Mi = gmm.mu(:);
    M_test(:,i) = Mi;
end
save('/mnt/g813_u7/kaldi-mfcc/olr2020_task3/data_sdc/M_test.mat','M_test','-v7.3');

b_train = [b_train;Y(25863:36094,1);Y(46319:75092,1);Y(84686:94282,1)]; 
save('/mnt/g813_u7/kaldi-mfcc/olr2020_task3/data_sdc/b_train.mat','b_train','-v7.3');



%% ==============================task_2(SDC)=============================================
% ubm -- 138*512 = 70656
ubm = importdata('/mnt/g813_u7/kaldi-mfcc/olr2020_task2/data_sdc/ubm_sdc_138_512.mat');

% task3的训练集(注册集) -- 3种语言 -- 26749条语音
enrollMfcc = importdata('/mnt/g813_u7/kaldi-mfcc/olr2020_task2/data_sdc/enrollMfcc.mat');  % 26749*1
enroll_num = size(enrollMfcc,1);
M_enroll = zeros(70656,enroll_num);
for i = 1:enroll_num
    i
    mfcc = enrollMfcc(i,:);
    gmm = mapAdapt(mfcc,ubm,10.0,'m');
    Mi = gmm.mu(:);
    M_enroll(:,i) = Mi;
end
save('/mnt/g813_u7/kaldi-mfcc/olr2020_task2/data_sdc/M_enroll.mat','M_enroll','-v7.3');

% task3的训练集 -- 6种语言 -- 11399条语音
% testMfcc = importdata('/mnt/g813_u7/kaldi-mfcc/olr2020_task2/data_sdc/testMfcc.mat');  % 11399*1
test_num = size(testMfcc,1);
M_test = zeros(70656,test_num);
for i = 1:test_num
    i
    mfcc = testMfcc(i,:);
    gmm = mapAdapt(mfcc,ubm,10.0,'m');
    Mi = gmm.mu(:);
    M_test(:,i) = Mi;
end
save('/mnt/g813_u7/kaldi-mfcc/olr2020_task2/data_sdc/M_test.mat','M_test','-v7.3');













