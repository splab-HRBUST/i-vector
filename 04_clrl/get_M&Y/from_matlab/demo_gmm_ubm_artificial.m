% 一个提取gmm-ubm的例子

%% 
% Step0: Set the parameters of the experiment   第0步：设置实验参数 -- 20个说话人每人10条语音
nSpeakers = 20;
nDims = 13;             % dimensionality of feature vectors  特征向量的维数
nMixtures = 32;         % How many mixtures used to generate data  有多少种混合物用于生成数据
nChannels = 10;         % Number of channels (sessions) per speaker  每个扬声器的频道（会话）数量
nFrames = 1000;         % Frames per speaker (10 seconds assuming 100 Hz)  每个扬声器的帧数（10秒，假设100赫兹）
nWorkers = 1;           % Number of parfor workers, if available  工人的数量（如有）

% Pick random centers for all the mixtures.    选择所有混合物的随机中心。
mixtureVariance = .10;
channelVariance = .05;
mixtureCenters = randn(nDims, nMixtures, nSpeakers);  % 返回一个随机数组（13*32*20）
channelCenters = randn(nDims, nMixtures, nSpeakers, nChannels)*.1;   %返回一个随机数组（13*32*20*10）*0.1
trainSpeakerData = cell(nSpeakers, nChannels);  %生成一个（20*10）的cell数组 -- 每一个cell中保存一个mfcc特征
testSpeakerData = cell(nSpeakers, nChannels);
speakerID = zeros(nSpeakers, nChannels);

% Create the random data. Both training and testing data have the same    创建随机数据。训练和测试数据都有相同的结果
% layout.    布局
for s=1:nSpeakers
    trainSpeechData = zeros(nDims, nMixtures);
    testSpeechData = zeros(nDims, nMixtures);
    for c=1:nChannels
        for m=1:nMixtures
            % Create data from mixture m for speaker s    从混合m为扬声器s创建数据
            frameIndices = m:nMixtures:nFrames;
            nMixFrames = length(frameIndices);
            trainSpeechData(:,frameIndices) = ...
                randn(nDims, nMixFrames)*sqrt(mixtureVariance) + ...
                repmat(mixtureCenters(:,m,s),1,nMixFrames) + ...
                repmat(channelCenters(:,m,s,c),1,nMixFrames);
            testSpeechData(:,frameIndices) = ...
                randn(nDims, nMixFrames)*sqrt(mixtureVariance) + ...
                repmat(mixtureCenters(:,m,s),1,nMixFrames) + ...
                repmat(channelCenters(:,m,s,c),1,nMixFrames);
        end
        trainSpeakerData{s, c} = trainSpeechData;  % 每一个cell都是13*1000 （20*10个cell）
        testSpeakerData{s, c} = testSpeechData;  % 每一个cell都是13*1000 （20*10个cell）
        speakerID(s,c) = s;   % Keep track of who this is    追踪这是谁（20*10）
    end
end

%%
% Step1: Create the universal background model from all the training speaker data    步骤1：从所有的训练说话人数据创建通用背景模型
nmix = nMixtures;           % In this case, we know the # of mixtures needed    在这种情况下，我们知道所需混合物的#
final_niter = 10;
ds_factor = 1;
% 提取到ubm{w-1*32;mu-13*32;sigma:13*32}
ubm = gmm_em(trainSpeakerData(:), nmix, final_niter, ds_factor, nWorkers);  % nmix-1024  final_niter-迭代次数

%%
% Step2: Now adapt the UBM to each speaker to create GMM speaker model.    步骤2：现在根据每个扬声器调整UBM，以创建GMM扬声器模型。
map_tau = 10.0;
config = 'mwv';
gmm = cell(nSpeakers, 1);
for s=1:nSpeakers
    gmm{s} = mapAdapt(trainSpeakerData(s, :), ubm, map_tau, config);
end

%%
% Step3: Now calculate the score for each model versus each speaker's data.    第三步：现在根据每位演讲者的数据计算每个模型的分数。
% Generate a list that tests each model (first column) against all the    生成一个列表，该列表根据所有
% testSpeakerData.
trials = zeros(nSpeakers*nChannels*nSpeakers, 2);
answers = zeros(nSpeakers*nChannels*nSpeakers, 1);
for ix = 1 : nSpeakers,
    b = (ix-1)*nSpeakers*nChannels + 1;
    e = b + nSpeakers*nChannels - 1;
    trials(b:e, :)  = [ix * ones(nSpeakers*nChannels, 1), (1:nSpeakers*nChannels)'];
    answers((ix-1)*nChannels+b : (ix-1)*nChannels+b+nChannels-1) = 1;
end

gmmScores = score_gmm_trials(gmm, reshape(testSpeakerData', nSpeakers*nChannels,1), trials, ubm);

%%
% Step4: Now compute the EER and plot the DET curve and confusion matrix    步骤4：现在计算EER，绘制DET曲线和混淆矩阵
imagesc(reshape(gmmScores,nSpeakers*nChannels, nSpeakers))
title('Speaker Verification Likelihood (GMM Model)');
ylabel('Test # (Channel x Speaker)'); xlabel('Model #');
colorbar; drawnow; axis xy
figure
eer = compute_eer(gmmScores, answers, true);
