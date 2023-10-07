% ======== 注册集 ========
file_enroll = get_mfcc_path('/mnt/g813_u7/kaldi-mfcc/olr2020_task2/htk_mfcc/enroll/');  % (26749,1)
num_enroll = size(file_enroll,1);  % 26749
fid = fopen('/mnt/g813_u7/triplet-mine_1/Noisy_LID/trainFrame.txt','wt');  
for i = 1 :num_enroll
    i
    file_path = file_enroll(i,1);  % "/mnt/g813_u7/kaldi-mfcc/olr2020_task2/htk_sdc/enroll/Minnan-000050001.mfcc"
    split_1 = strsplit(file_path,'/');
    split_2 = strsplit(split_1(1,8),'.');
    file_name = split_2(1,1);
    mfcc_path = file_enroll(i,:);  % 第i个音频的mfcc地址
    mfcc_data = htkread(mfcc_path);  % 69*764 -- 69维的mfcc特征
    save_path = sprintf('%s%s%s','/mnt/g813_u7/triplet-mine_1/Noisy_LID/trainMfcc/',file_name,'.mat');
    save(save_path,'mfcc_data');
    frame = size(mfcc_data,2);  % 764
    fprintf(fid,'%d\n',frame);
end
fclose(fid);

% ======== 测试集 ========
file_test = get_mfcc_path('/mnt/g813_u7/kaldi-mfcc/olr2020_task2/htk_mfcc/test/');
num_test = size(file_test,1);  % 11399
fid = fopen('/mnt/g813_u7/triplet-mine_1/Noisy_LID/testFrame.txt','wt'); 
for i = 1 :num_test
    i
    file_path = file_test(i,1);  % "/mnt/g813_u7/kaldi-mfcc/olr2020_task2/htk_mfcc/test/minnan-10001.mfcc"
    split_1 = strsplit(file_path,'/');
    split_2 = strsplit(split_1(1,8),'.');
    file_name = split_2(1,1);
    mfcc_path = file_test(i,:);  % 第i个音频的mfcc地址
    mfcc_data = htkread(mfcc_path);  % 69*616 -- 69维的mfcc特征
    save_path = sprintf('%s%s%s','/mnt/g813_u7/triplet-mine_1/Noisy_LID/testMfcc/',file_name,'.mat');
    save(save_path,'mfcc_data');
    frame = size(mfcc_data,2);  % 616
    fprintf(fid,'%d\n',frame);
end
fclose(fid);

