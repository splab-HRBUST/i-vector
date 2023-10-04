function spkRecog(ubm, patameters_FA, spkInfo)

    % 测试者发言
    fprintf('将进行说话人匹配...\n');pause(1);
    fprintf('3秒后开始录音...\n');pause(1);
    fprintf('3\t');pause(1);
    fprintf('2\t');pause(1);
    fprintf('1\n');pause(1);
    fprintf('开始录音，请坚持讲话5秒：\n');
    recorder = audiorecorder(8000,8,1); %fs=8000, nbit=8, nchannel=1
    recordblocking(recorder, 5);
    myRecording = getaudiodata(recorder);
    savePath = fullfile('Data/DataBase_test','test .wav');
    audiowrite(savePath,myRecording,8000);
    
    % 提取iv
% 	savePath = fullfile('Data/DataBase_Enrollment',char([spkInfo.spkName '.wav']));
    mmm = extract_mfcc(savePath,'0');
    gmm = mapAdapt(mmm, ubm, 10.0, 'm');
    M = gmm.mu(:);
    Ex_test = pinv(patameters_FA.L)*patameters_FA.B*(M-patameters_FA.m);
    
    % 说话人匹配
    for i = 1 : size(spkInfo,2)
        Ex_model(:,i) = spkInfo{i}.Ex;
    end
    scores_CDS = 1 - pdist2(Ex_test',Ex_model','cosine');
    [v,id] = max(scores_CDS);
    if v >= 0.25
        fprintf('说话人为：%s\n\n',spkInfo{id}.spkName);
    else
        fprintf('查无此人！\n\n');
    end
end