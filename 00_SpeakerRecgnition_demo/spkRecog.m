function spkRecog(ubm, patameters_FA, spkInfo)

    % �����߷���
    fprintf('������˵����ƥ��...\n');pause(1);
    fprintf('3���ʼ¼��...\n');pause(1);
    fprintf('3\t');pause(1);
    fprintf('2\t');pause(1);
    fprintf('1\n');pause(1);
    fprintf('��ʼ¼�������ֽ���5�룺\n');
    recorder = audiorecorder(8000,8,1); %fs=8000, nbit=8, nchannel=1
    recordblocking(recorder, 5);
    myRecording = getaudiodata(recorder);
    savePath = fullfile('Data/DataBase_test','test .wav');
    audiowrite(savePath,myRecording,8000);
    
    % ��ȡiv
% 	savePath = fullfile('Data/DataBase_Enrollment',char([spkInfo.spkName '.wav']));
    mmm = extract_mfcc(savePath,'0');
    gmm = mapAdapt(mmm, ubm, 10.0, 'm');
    M = gmm.mu(:);
    Ex_test = pinv(patameters_FA.L)*patameters_FA.B*(M-patameters_FA.m);
    
    % ˵����ƥ��
    for i = 1 : size(spkInfo,2)
        Ex_model(:,i) = spkInfo{i}.Ex;
    end
    scores_CDS = 1 - pdist2(Ex_test',Ex_model','cosine');
    [v,id] = max(scores_CDS);
    if v >= 0.25
        fprintf('˵����Ϊ��%s\n\n',spkInfo{id}.spkName);
    else
        fprintf('���޴��ˣ�\n\n');
    end
end