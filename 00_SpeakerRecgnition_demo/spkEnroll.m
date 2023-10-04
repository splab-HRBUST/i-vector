function spkInfo = spkEnroll(ubm, patameters_FA)
    % ע������
    prompt = '������������';
    spkInfo.spkName = input(prompt,'s');
    fprintf('3���ʼ¼��...\n');pause(1);
    fprintf('3\t');pause(1);
    fprintf('2\t');pause(1);
    fprintf('1\n');pause(1);
   
    % ¼��
    fprintf('��ʼ¼�������ֽ���5�룺\n');
    recorder = audiorecorder(8000,8,1); %fs=8000, nbit=8, nchannel=1
    recordblocking(recorder, 5);
    myRecording = getaudiodata(recorder);
    savePath = fullfile('Data/DataBase_Enrollment',char([spkInfo.spkName '.wav']));
    audiowrite(savePath,myRecording,8000);
    fprintf('¼��������\n');
    
    % ��ȡiv
    mmm = extract_mfcc(savePath,'0');
    gmm = mapAdapt(mmm, ubm, 10.0, 'm');
    M = gmm.mu(:);
    spkInfo.Ex = pinv(patameters_FA.L)*patameters_FA.B*(M-patameters_FA.m);
    fprintf('˵����ע��ɹ���\n\n');

end