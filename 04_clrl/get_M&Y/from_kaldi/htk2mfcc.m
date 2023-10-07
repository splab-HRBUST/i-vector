%% ===============================task1--MFCC特征======================================
% 训练集mfcc地址 -- train
file_train = get_mfcc_path('/mnt/g813_u7/kaldi-mfcc/olr2020/htk_mfcc/train/');
num_train = size(file_train,1);  % 94282
save('/mnt/g813_u7/kaldi-mfcc/olr2020/data_mfcc/file_train.mat','file_train','-v7.3');
trainMfcc = cell(num_train,1);  % 94282个
for i = 1 :num_train
    i
    mfcc_path = file_train(i,:);  % 第i个音频的mfcc地址
    mfcc_data = htkread(mfcc_path);  % 23*463 -- 23维的mfcc特征
    trainMfcc{i,1} = mfcc_data;  % 23*463
end
save('/mnt/g813_u7/kaldi-mfcc/olr2020/data_mfcc/trainMfcc.mat','trainMfcc','-v7.3');


% ======================================================================
% 注册集mfcc地址 -- enroll
file_enroll = get_mfcc_path('/mnt/g813_u7/kaldi-mfcc/olr2020/htk_mfcc/enroll/');
num_enroll = size(file_enroll,1);  % 58823
save('/mnt/g813_u7/kaldi-mfcc/olr2020/data_mfcc/file_enroll.mat','file_enroll','-v7.3');
enrollMfcc = cell(num_enroll,1);  % 58823个
for i = 1 :num_enroll
    i
    mfcc_path = file_enroll(i,:);  % 第i个音频的mfcc地址
    mfcc_data = htkread(mfcc_path);  % 23*463 -- 23维的mfcc特征
    enrollMfcc{i,1} = mfcc_data;  % 23*463
end
save('/mnt/g813_u7/kaldi-mfcc/olr2020/data_mfcc/enrollMfcc.mat','enrollMfcc','-v7.3');

% ======================================================================
% 测试集mfcc地址 -- test
file_test = get_mfcc_path('/mnt/g813_u7/kaldi-mfcc/olr2020/htk_mfcc/test/');
num_test = size(file_test,1);  % 11848
save('/mnt/g813_u7/kaldi-mfcc/olr2020/data_mfcc/file_test.mat','file_test','-v7.3');
testMfcc = cell(num_test,1);  % 11848个
for i = 1 :num_test
    i
    mfcc_path = file_test(i,:);  % 第i个音频的mfcc地址
    mfcc_data = htkread(mfcc_path);  % 23*268 -- 23维的mfcc特征
    testMfcc{i,1} = mfcc_data;  % 23*268
end
save('/mnt/g813_u7/kaldi-mfcc/olr2020/data_mfcc/testMfcc.mat','testMfcc','-v7.3');


%% ====================================task1--SDC特征====================================
% 训练集sdc地址 -- train
file_train = get_mfcc_path('/mnt/g813_u7/kaldi-mfcc/olr2020/htk_sdc/train/');
num_train = size(file_train,1);  % 94282
save('/mnt/g813_u7/kaldi-mfcc/olr2020/data_sdc/file_train.mat','file_train','-v7.3');
trainMfcc = cell(num_train,1);  % 94282个
for i = 1 :num_train
    i
    mfcc_path = file_train(i,:);  % 第i个音频的mfcc地址
    mfcc_data = htkread(mfcc_path);  % 23*463 -- 23维的mfcc特征
    trainMfcc{i,1} = mfcc_data;  % 23*463
end
save('/mnt/g813_u7/kaldi-mfcc/olr2020/data_sdc/trainMfcc.mat','trainMfcc','-v7.3');


% ======================================================================
% 注册集sdc地址 -- enroll
file_enroll = get_mfcc_path('/mnt/g813_u7/kaldi-mfcc/olr2020/htk_sdc/enroll/');
num_enroll = size(file_enroll,1);  % 58823
save('/mnt/g813_u7/kaldi-mfcc/olr2020/data_sdc/file_enroll.mat','file_enroll','-v7.3');
enrollMfcc = cell(num_enroll,1);  % 58823个
for i = 1 :num_enroll
    i
    mfcc_path = file_enroll(i,:);  % 第i个音频的mfcc地址
    mfcc_data = htkread(mfcc_path);  % 23*463 -- 23维的mfcc特征
    enrollMfcc{i,1} = mfcc_data;  % 23*463
end
save('/mnt/g813_u7/kaldi-mfcc/olr2020/data_sdc/enrollMfcc.mat','enrollMfcc','-v7.3');

% ======================================================================
% 测试集sdc地址 -- test
file_test = get_mfcc_path('/mnt/g813_u7/kaldi-mfcc/olr2020/htk_sdc/test/');
num_test = size(file_test,1);  % 11848
save('/mnt/g813_u7/kaldi-mfcc/olr2020/data_sdc/file_test.mat','file_test','-v7.3');
testMfcc = cell(num_test,1);  % 11848个
for i = 1 :num_test
    i
    mfcc_path = file_test(i,:);  % 第i个音频的mfcc地址
    mfcc_data = htkread(mfcc_path);  % 23*268 -- 23维的mfcc特征
    testMfcc{i,1} = mfcc_data;  % 23*268
end
save('/mnt/g813_u7/kaldi-mfcc/olr2020/data_sdc/testMfcc.mat','testMfcc','-v7.3');



%%  ====================================task2_SDC特征====================================
% 注册集sdc地址 -- enroll
file_enroll = get_mfcc_path('/mnt/g813_u7/kaldi-mfcc/olr2020_task2/htk_sdc/enroll/');
num_enroll = size(file_enroll,1);  % 26749
save('/mnt/g813_u7/kaldi-mfcc/olr2020_task2/data_sdc/file_enroll.mat','file_enroll','-v7.3');
enrollMfcc = cell(num_enroll,1);  % 26749个
for i = 1 :num_enroll
    i
    mfcc_path = file_enroll(i,:);  % 第i个音频的mfcc地址
    mfcc_data = htkread(mfcc_path);  % 23*463 -- 23维的mfcc特征
    enrollMfcc{i,1} = mfcc_data;  % 23*463
end
save('/mnt/g813_u7/kaldi-mfcc/olr2020_task2/data_sdc/enrollMfcc.mat','enrollMfcc','-v7.3');


% 测试集sdc地址 -- test
file_test = get_mfcc_path('/mnt/g813_u7/kaldi-mfcc/olr2020_task2/htk_sdc/test/');
num_test = size(file_test,1);  % 11399
save('/mnt/g813_u7/kaldi-mfcc/olr2020_task2/data_sdc/file_test.mat','file_test','-v7.3');
testMfcc = cell(num_test,1);  % 11399个
for i = 1 :num_test
    i
    mfcc_path = file_test(i,:);  % 第i个音频的mfcc地址
    mfcc_data = htkread(mfcc_path);  % 23*268 -- 23维的mfcc特征
    testMfcc{i,1} = mfcc_data;  % 23*268
end
save('/mnt/g813_u7/kaldi-mfcc/olr2020_task2/data_sdc/testMfcc.mat','testMfcc','-v7.3');



%%  ====================================task3_SDC特征====================================
% 测试集sdc地址 -- test
file_test = get_mfcc_path('/mnt/g813_u7/kaldi-mfcc/olr2020_task3/htk_sdc/test/');
num_test = size(file_test,1);  % 9496
save('/mnt/g813_u7/kaldi-mfcc/olr2020_task3/data_sdc/file_test.mat','file_test','-v7.3');
testMfcc = cell(num_test,1);  % 9496个
for i = 1 :num_test
    i
    mfcc_path = file_test(i,:);  % 第i个音频的mfcc地址
    mfcc_data = htkread(mfcc_path);  % 23*268 -- 23维的mfcc特征
    testMfcc{i,1} = mfcc_data;  % 23*268
end
save('/mnt/g813_u7/kaldi-mfcc/olr2020_task3/data_sdc/testMfcc.mat','testMfcc','-v7.3');












