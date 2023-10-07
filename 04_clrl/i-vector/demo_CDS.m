
M_train = importdata('/mnt/g813_u7/kaldi-mfcc/olr2020_task1/data_sdc/M_train.mat');
M_enroll = importdata('/mnt/g813_u7/kaldi-mfcc/olr2020_task1/data_sdc/M_enroll.mat');
M_test = importdata('/mnt/g813_u7/kaldi-mfcc/olr2020_task1/data_sdc/M_test.mat');

Y_train = importdata('/mnt/g813_u7/kaldi-mfcc/olr2020_task1/data_sdc/Y_train.mat');
Y_enroll = importdata('/mnt/g813_u7/kaldi-mfcc/olr2020_task1/data_sdc/Y_enroll.mat');
Y_test = importdata('/mnt/g813_u7/kaldi-mfcc/olr2020_task1/data_sdc/Y_test.mat');

num.spk_type = 10;  % 语种类别

num.train = 94282;  % 训练集中音频数量
num.enroll = 58823;  % 注册集中音频数量
num.test = 11848;  % 测试集中音频数量

num.IVdim = 400;  % 公共隐变量维度
num.nIters = 1;  %

% 改参数调试--IVdim(800,400),迭代次数(100,150,200)，nIters(50,100,150,200)

fid = fopen('olr2020_task1_CDS.txt','wt');
for iTers = 1 : 400  % 迭代次数 -- 进行100次迭代，查看是否收敛
    if iTers == 1
        [patameters_CLRL] = clrl_em(M_train,Y_train,num);

    else
        [patameters_CLRL] = clrl_em(M_train,Y_train,num,patameters_CLRL);
    end
    
    Ex_enroll = pinv(patameters_CLRL.L)*patameters_CLRL.B*(M_enroll-patameters_CLRL.m);  %注册集的i-vector特征 
    Ex_test   = pinv(patameters_CLRL.L)*patameters_CLRL.B*(M_test-patameters_CLRL.m);  % 测试集的i-vector特征

    % 标准化 -- 转置 
    ExEnroll = Ex_enroll';  ExEnroll = zscore(ExEnroll);  Ex_enorll = ExEnroll';
    ExTest = Ex_test';  ExTest = zscore(ExTest);  Ex_test = ExTest';
    
    % =================================================================================
    % 求出注册集与测试集中每个语种i-vector的平均值进行比较
    Ex_model=zeros(num.IVdim,num.spk_type); % 保存注册集中每个语种的平均i-vector特征
    Num = zeros(1,num.spk_type);  %
    % 同语种的i-vector特征相加
    for i = 1 : num.enroll
        y = find(Y_enroll(:,i)==1); % 注册集中第i个音频种类为y
        Ex_model(:,y) = Ex_model(:,y)+Ex_enroll(:,i); % 在第y种类位置加上第i个音频
        Num(1,y) = Num(1,y)+1;  %语种数量+1
    end
    % 求平均
    for i = 1: num.spk_type
        if Num(1,i) > 0
            Ex_model(:,i) = Ex_model(:,i)./Num(1,i);  % 求平均 -- Ex_model(400*6) -- 得到每个语种的平均i-vector
        end
    end
    
     % 标准化 -- 转置
    ExModel = Ex_model';  ExModel = zscore(ExModel);  Ex_model = ExModel';
    
    fprintf('================= cosine ====================\n');
    scores_CLRL_CDS = [];  
    scores_CLRL_CDS.all =1 - pdist2(Ex_test',Ex_model','cosine');  
    
    scores_CLRL_CDS.true = zeros(num.test,1);  % 同一语种  
    num.impostor = num.test*(num.spk_type-1);  % 不同语种
    scores_CLRL_CDS.impostor = zeros(num.impostor,1);  % 不同语种  
    
    file_socre = strcat('/mnt/g813_u7/kaldi-mfcc/olr2020_result/task_1/score(',num.IVdim,')_',method,'/score_',num2str(iTers),'.txt');
    file_target = strcat('/mnt/g813_u7/kaldi-mfcc/olr2020_result/task_1/target/target_',num2str(iTers),'.txt');
    
    fid_score = fopen(file_socre,'wt');
    fid_target = fopen(file_target,'wt');
    
    t=0;f=0;
    for a =  1 : num.test  % 测试集音频数
        for b = 1 : num.spk_type  % 语种类别数
            type = find(Y_test(:,a)==1);  % 测试集中第a个音频的语种类别type
            A = scores_CLRL_CDS.all(a , b);  % 说话人a与语种b之间的余弦相似度
            if type == b  % 是同一语种 -- true -- 即类别标签为1
                t = t+1;
                scores_CLRL_CDS.true(t,1) = A;
            elseif type~=b  % 不是同一语种 -- impostor -- 即类别标签为0
                f = f+1;
                scores_CLRL_CDS.impostor(f,1) = A;
            end
        end
    end  
    clear a b A
    
    answer_eva_king = [ones(1,t) zeros(1,f)];  % 标签[同一语种, 不同语种]
    scores_CLRL_CDS = [scores_CLRL_CDS.true;scores_CLRL_CDS.impostor]; 
    [eer_CLRL_CDS(iTers),~,~,~,dcf_CLRL_CDS(iTers)]=compute_eer(scores_CLRL_CDS,answer_eva_king,false);
    
    % ================================================================================= 

    fprintf('===== Summaray ======= Baseline =====%d=========\n',iTers);
    fprintf('+-------------------------------------------+\n');
    fprintf('|    Method     |    EER(%%)   | Min DCF_vox |\n');
    fprintf('+---------------+--------------+-------------+\n');
    fprintf('|CLRL(%d)+CDS   |    %2.2f     |    %2.4f    |\n', num.IVdim,eer_CLRL_CDS(iTers),dcf_CLRL_CDS(iTers));
    fprintf('+---------------+--------------+-------------+\n');
    fprintf(fid,'|CLRL(%d)+CDS   |    %2.2f     |    %2.4f    |\n', num.IVdim,eer_CLRL_CDS(iTers),dcf_CLRL_CDS(iTers));
end
fclose(fid);
