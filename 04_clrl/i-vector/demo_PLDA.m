
M_train = importdata('/mnt/g813_u7/kaldi-mfcc/olr2020_task1/data_sdc/M_enroll.mat');
M_enroll = importdata('/mnt/g813_u7/kaldi-mfcc/olr2020_task1/data_sdc/M_enroll.mat');
M_test = importdata('/mnt/g813_u7/kaldi-mfcc/olr2020_task1/data_sdc/M_test.mat');

Y_train = importdata('/mnt/g813_u7/kaldi-mfcc/olr2020_task1/data_sdc/Y_enroll.mat');
Y_enroll = importdata('/mnt/g813_u7/kaldi-mfcc/olr2020_task1/data_sdc/Y_enroll.mat');
Y_test = importdata('/mnt/g813_u7/kaldi-mfcc/olr2020_task1/data_sdc/Y_test.mat');


num.spk_type = 10;  % 语种类别

num.train = 94282;  % 训练集中音频数量
num.enroll = 58823;  % 注册集中音频数量
num.test = 11848;  % 测试集中音频数量

num.IVdim = 400;  % 公共隐变量维度
num.nIters = 1;  %

num.label_train = importdata('/mnt/g813_u7/kaldi-mfcc/olr2020_task1/data_sdc/b_train.txt');  % 训练集的标签
ubm = importdata('/mnt/g813_u7/kaldi-mfcc/olr2020_task1/data_sdc/ubm_sdc_138_512.mat');

% 开始训练
save_file = sprintf('olr2020_task1_PLDA.txt');
fid = fopen(save_file,'wt');
for iTers = 1 : 50  % 迭代次数 -- 进行100次迭代，查看是否收敛
    if iTers == 1
        [patameters_CLRL] = clrl_em(M_train, ubm, num);

    else
        [patameters_CLRL] = clrl_em(M_train, ubm, num, patameters_CLRL);
    end

    Ex_train =  pinv(patameters_CLRL.L)*patameters_CLRL.B*(M_train-patameters_CLRL.m);  % 训练集的i-vector特征
    Ex_enroll = pinv(patameters_CLRL.L)*patameters_CLRL.B*(M_enroll-patameters_CLRL.m);  %注册集的i-vector特征 
    Ex_test   = pinv(patameters_CLRL.L)*patameters_CLRL.B*(M_test-patameters_CLRL.m);  % 测试集的i-vector特征

     % 标准化
    ExTrain = Ex_train';  ExTrain = zscore(ExTrain);  Ex_train = ExTrain';
    ExEnroll = Ex_enroll';  ExEnroll = zscore(ExEnroll);  Ex_enorll = ExEnroll';
    ExTest = Ex_test';  ExTest = zscore(ExTest);  Ex_test = ExTest';

    % =================================================================================
    % 求出注册集与测试集中每个语种i-vector的平均值进行比较
    Ex_model=zeros(num.IVdim,num.spk_type); % 保存注册集中每个语种的平均i-vector特征
    Num = zeros(1,num.spk_type); 
    % 同语种的i-vector特征相加
    for i = 1 : num.enroll
        y = find(Y_enroll(:,i)==1); % 注册集中第i个音频种类为y
        Ex_model(:,y) = Ex_model(:,y)+Ex_enroll(:,i); % 在第y种类位置加上第i个音频
        Num(1,y) = Num(1,y)+1;  %语种数量+1
    end
    % 求平均
    for i = 1: num.spk_type
        if Num(1,i) > 0
            Ex_model(:,i) = Ex_model(:,i)./Num(1,i);  
        end
    end

    % 标准化
    ExModel = Ex_model';  ExModel = zscore(ExModel);  Ex_model = ExModel';

    fprintf('================= PLDA ====================\n');
    num.Zdim = 300;
    label_train = num.num.label_train;
    [V,D] = lda(Ex_train,label_train(:));

    finalEx = V(:,1:5)'*Ex_train;
    finalEx_model = V(:,1:5)'*Ex_model;
    finalEx_test = V(:,1:5)'*Ex_test;

    pLDA_bpls = gplda_em(finalEx, label_train(:), num.Zdim, 50); % (1*1)

    scores_BPLS_PLDA = [];
    scores_BPLS_PLDA.all = (score_gplda_trials(pLDA_bpls, finalEx_model, finalEx_test))';

    scores_BPLS_PLDA.true = [];  % 同一语种
    scores_BPLS_PLDA.impostor = [];  % 不同语种

    scores_BPLS_PLDA.true = zeros(num.test,1);  % 同一语种  22050*1
    num.impostor = num.test*(num.spk_type-1);  % 198450
    scores_BPLS_PLDA.impostor = zeros(num.impostor,1);  % 不同语种  198450*1

    file_socre = strcat('/mnt/g813_u7/kaldi-mfcc/olr2020_result/task_1/score_',num2str(iTers),'.txt');
    file_target = strcat('/mnt/g813_u7/kaldi-mfcc/olr2020_result/task_1/target_',num2str(iTers),'.txt');

    fid_score = fopen(file_socre,'wt');
    fid_target = fopen(file_target,'wt');

    t=0;f=0;
    for a =  1 : num.test  % 测试集音频数
        for b = 1 : num.spk_type  % 语种类别数
            type = find(Y_test(:,a)==1);  % 测试集中第a个音频的语种类别type
            A = scores_BPLS_PLDA.all(a , b);  % 说话人a与语种b之间的plda得分
            if type == b  % 是同一语种
                t = t+1;
                scores_BPLS_PLDA.true(t,1) = A;
                fprintf(fid_target,'language%d test%d target\n',b,a);
            elseif type~=b  % 不是同一语种
                f = f+1;
                scores_BPLS_PLDA.impostor(f,1) = A;
                fprintf(fid_target,'language%d test%d nontarget\n',b,a);
            end
            fprintf(fid_score,'language%d test%d %2.6f\n',b,a,A);
        end
    end  
    clear a b A
    fclose(fid_target);
    fclose(fid_score);
    
    answer_eva_king = [ones(1,t) zeros(1,f)];  % 标签[同一语种, 不同语种]
    scores_CLRL_PLDA = [scores_CLRL_PLDA.true;scores_CLRL_PLDA.impostor]; 
    [eer_CLRL_PLDA(iTers),~,~,~,dcf_CLRL_PLDA(iTers)]=compute_eer(scores_CLRL_PLDA,answer_eva_king,false);

    % ================================================================================= 

    fprintf('===== Summaray ======= Baseline =====%d=========\n',iTers);
    fprintf('+-------------------------------------------+\n');
    fprintf('|    Method     |    EER(%%)   | Min DCF_vox |\n');
    fprintf('+---------------+--------------+-------------+\n');
    fprintf('|PPLS(%d)+PLDA(%d)|    %2.2f    |   %2.4f    |\n', num.IVdim,num.Zdim,eer_CLRL_PLDA(iTers),dcf_CLRL_PLDA(iTers));
    fprintf('+---------------+--------------+-------------+\n');
  
end
fclose(fid);


