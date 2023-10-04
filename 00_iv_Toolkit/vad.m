function [finalSeg, fs, ind] = vad(file,channel)
%--------------------------------------------------
%input: 
%    file:�����ļ�·��
%    channel:��Ҫ�˵�����ŵ�
%           'w'��whole��ȫ��
%           'a'��a�ŵ�������һ�ŵ�
%           'b'��b�ŵ������ڶ��ŵ�
%output:
%    seg:ȥ�������ε��ź�
%    fs: ������
%    nbits:����λ��
%--------------------------------------------------
% �ж��ļ���ʽ��Ŀǰ���Դ�����ļ���ʽΪ��wav��sph
if ischar(file)
    if  strcmp(file(end-2:end),'wav') == 1
        [X,fs] = audioread(file);
    elseif strcmp(file(end-2:end),'sph') == 1
        [X,fs,~] = readsph(file);
    %    nbits = inf{2,1}{8,2};
    end
end
% ��һ����ÿ���ŵ�
% spkChannel = 0;     % ��¼�����ŵ�����
switch channel        % ��Ҫ�˵�����ŵ�
    case 'w',
        finalSeg = cell(1,2);
        for nChannel = 1 : size(X,2)
            [finalSeg{:,nChannel},ind] = getSeg(X(:,nChannel));
        end
    case 'a',
        [finalSeg,ind] = getSeg(X(:,1));
    case 'b',
        [finalSeg,ind] = getSeg(X(:,2));
end

%--------------------------------------------------
function [finalSeg,ind] = getSeg(XX)

        ind = [];
        x_stop = 0;

        x=XX/max(abs(XX));          % ���ȹ�һ����[-1,1]

        % ��������
        FrameLen = 512;     % ֡��
        inc = 160;          % δ�ص�����

        amp1 = 10;          % ��ʱ������ֵ
        amp2 = 2;      
        zcr1 = 10;          % ��������ֵ
        zcr2 = 5;

        minsilence = 40;    % �������ĳ������ж������Ƿ����
        minlen  = 15;       % �ж�����������С����

%--------------------------------------------------
        seg = [];
        fragment = [];
        while size(x,1) > 0  % ʣ����������0ʱ�������˵���
            status  = 0;     % ��ʼ�������ε�״̬
            count   = 0;     % ��ʼ���������еĳ���
            silence = 0;     % ��ʼ�������ĳ���
            % ���������
            tmp1  = enframe(x(1:end-1), FrameLen,inc);
            tmp2  = enframe(x(2:end)  , FrameLen,inc);
            signs = (tmp1.*tmp2)<0;
            diffs = (tmp1 -tmp2)>0.02;
            zcr   = sum(signs.*diffs,2);

            % �����ʱ����
            amp = sum((abs(enframe(filter([1 -0.9375], 1, x), FrameLen, inc))).^2, 2);

            % ������������
            amp1 = min(amp1, max(amp)/4);
            amp2 = min(amp2, max(amp)/8);

            % ��ʼ�˵���
            x1 = 0;
            for n=1:length(zcr)
               switch status                %�����ε�״̬
               case {0,1}                   % 0 = ����, 1 = ���ܿ�ʼ
                  if amp(n) > amp1          % ȷ�Ž���������   %%������������ֵ
                     x1 = max(n-count-1,1); % ��¼�����ε���ʼ��
                     status  = 2;           % �����ε�״̬��Ϊ2
                     silence = 0;
                     count   = count + 1;
                  elseif amp(n) > amp2 || zcr(n) > zcr2 % ���ܴ���������
                     status = 1;                 
                     count  = count + 1;
                  else                       % ����״̬
                     status = 0;
                  %   count   = 0;
                  end
               case 2,                       % 2 = ������
                  if amp(n) > amp2 ||zcr(n) > zcr2     % ������������

                     count = count + 1;
                  else                       % ����������
                     silence = silence+1;
                     if silence < minsilence % ����������������δ����
                        count  = count + 1;
                     elseif count < minlen   % ��������̫�̣���Ϊ������
                        status  = 0;
                        silence = 0;
                        count   = 0;
                     else                    % ��������
                        status  = 3;
                     end
                  end
               case 3,
                  continue;
               end
            end   

%--------------------------------------------------
        count = count-silence/2;
        if x1 == 0  % ���ȫΪ����
            x = [];
            XX = [];
        else
            x2 = x1 + count -1;
            fragment = XX(x1*inc:x2*inc);
            
            % ��¼���
            x_start = x_stop + x1*inc;
            x_stop  = x_stop + x2*inc;
            ind = [ind,x_start,x_stop];
            
            seg = [seg;fragment];
            x = x(x2*inc:end);
            XX = XX(x2*inc:end);
        end
        if isempty(seg) % ĳһ�ŵ���˵��������
            spkChannel = spkChannel+1;
        end
%         seg = nonzeros(seg);          % ����0�Ĳ���ȫ��ȥ��
        finalSeg = seg;
    end


% if spkChannel > 0 && nChannel == 2
%     finalSeg{:,1} = [finalSeg{:,nChannel-1} finalSeg{:,nChannel}];
%     finalSeg = finalSeg(1);
% end
    
