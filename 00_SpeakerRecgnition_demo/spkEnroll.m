function spkInfo = spkEnroll(ubm, patameters_FA)
    % 注册姓名
    prompt = '请输入姓名：';
    spkInfo.spkName = input(prompt,'s');
    fprintf('3秒后开始录音...\n');pause(1);
    fprintf('3\t');pause(1);
    fprintf('2\t');pause(1);
    fprintf('1\n');pause(1);
   
    % 录音
    fprintf('开始录音，请坚持讲话5秒：\n');
    recorder = audiorecorder(8000,8,1); %fs=8000, nbit=8, nchannel=1
    recordblocking(recorder, 5);
    myRecording = getaudiodata(recorder);
    savePath = fullfile('Data/DataBase_Enrollment',char([spkInfo.spkName '.wav']));
    audiowrite(savePath,myRecording,8000);
    fprintf('录音结束。\n');
    
    % 提取iv
    mmm = extract_mfcc(savePath,'0');
    gmm = mapAdapt(mmm, ubm, 10.0, 'm');
    M = gmm.mu(:);
    spkInfo.Ex = pinv(patameters_FA.L)*patameters_FA.B*(M-patameters_FA.m);
    fprintf('说话人注册成功。\n\n');

end