function M = spkEnroll(path, ubm)
    % ��ȡiv
    mmm = extract_mfcc(path,'0');  % ��ȡmfcc -- 60*433
    gmm = mapAdapt(mmm, ubm, 10.0, 'm');  %��ȡÿһ����Ƶ�ļ���GMMģ��
    % gmm --> w(1*1024),mu(60*1024),sigma(60*1024)
    M = gmm.mu(:);  % ��ȡGMMģ�͵ľ�ֵ��ʸ��
end
