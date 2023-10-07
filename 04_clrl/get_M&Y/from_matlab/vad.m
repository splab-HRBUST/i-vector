function [finalSeg, fs, ind] = vad(file,channel)
%--------------------------------------------------
%input: 
%    file:输入文件路径
%    channel:需要端点检测的信道
%           'w'：whole，全部
%           'a'：a信道，即第一信道
%           'b'：b信道，即第二信道
%output:
%    seg:去除静音段的信号
%    fs: 采样率
%    nbits:数据位数
%--------------------------------------------------
% 判断文件格式，目前可以处理的文件格式为：wav，sph
if ischar(file)
    if  strcmp(file(end-2:end),'wav') == 1
        [X,fs] = audioread(file);
    elseif strcmp(file(end-2:end),'sph') == 1
        [X,fs,~] = readsph(file);
    %    nbits = inf{2,1}{8,2};
    end
end
% 逐一处理每个信道
% spkChannel = 0;     % 记录静音信道个数
switch channel        % 需要端点检测的信道
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

        x=XX/max(abs(XX));          % 幅度归一化到[-1,1]

        % 参数设置
        FrameLen = 512;     % 帧长
        inc = 160;          % 未重叠部分

        amp1 = 10;          % 短时能量阈值
        amp2 = 2;      
        zcr1 = 10;          % 过零率阈值
        zcr2 = 5;

        minsilence = 40;    % 用无声的长度来判断语音是否结束
        minlen  = 15;       % 判断是语音的最小长度

%--------------------------------------------------
        seg = [];
        fragment = [];
        while size(x,1) > 0  % 剩余采样点大于0时，继续端点检测
            status  = 0;     % 初始化语音段的状态
            count   = 0;     % 初始化语音序列的长度
            silence = 0;     % 初始化无声的长度
            % 计算过零率
            tmp1  = enframe(x(1:end-1), FrameLen,inc);
            tmp2  = enframe(x(2:end)  , FrameLen,inc);
            signs = (tmp1.*tmp2)<0;
            diffs = (tmp1 -tmp2)>0.02;
            zcr   = sum(signs.*diffs,2);

            % 计算短时能量
            amp = sum((abs(enframe(filter([1 -0.9375], 1, x), FrameLen, inc))).^2, 2);

            % 调整能量门限
            amp1 = min(amp1, max(amp)/4);
            amp2 = min(amp2, max(amp)/8);

            % 开始端点检测
            x1 = 0;
            for n=1:length(zcr)
               switch status                %语音段的状态
               case {0,1}                   % 0 = 静音, 1 = 可能开始
                  if amp(n) > amp1          % 确信进入语音段   %%能量超过高阈值
                     x1 = max(n-count-1,1); % 记录语音段的起始点
                     status  = 2;           % 语音段的状态改为2
                     silence = 0;
                     count   = count + 1;
                  elseif amp(n) > amp2 || zcr(n) > zcr2 % 可能处于语音段
                     status = 1;                 
                     count  = count + 1;
                  else                       % 静音状态
                     status = 0;
                  %   count   = 0;
                  end
               case 2,                       % 2 = 语音段
                  if amp(n) > amp2 ||zcr(n) > zcr2     % 保持在语音段

                     count = count + 1;
                  else                       % 语音将结束
                     silence = silence+1;
                     if silence < minsilence % 静音还不够长，尚未结束
                        count  = count + 1;
                     elseif count < minlen   % 语音长度太短，认为是噪声
                        status  = 0;
                        silence = 0;
                        count   = 0;
                     else                    % 语音结束
                        status  = 3;
                     end
                  end
               case 3,
                  continue;
               end
            end   

%--------------------------------------------------
        count = count-silence/2;
        if x1 == 0  % 后段全为静音
            x = [];
            XX = [];
        else
            x2 = x1 + count -1;
            fragment = XX(x1*inc:x2*inc);
            
            % 记录结点
            x_start = x_stop + x1*inc;
            x_stop  = x_stop + x2*inc;
            ind = [ind,x_start,x_stop];
            
            seg = [seg;fragment];
            x = x(x2*inc:end);
            XX = XX(x2*inc:end);
        end
        if isempty(seg) % 某一信道无说话人语音
            spkChannel = spkChannel+1;
        end
%         seg = nonzeros(seg);          % 把是0的部分全部去掉
        finalSeg = seg;
    end


% if spkChannel > 0 && nChannel == 2
%     finalSeg{:,1} = [finalSeg{:,nChannel-1} finalSeg{:,nChannel}];
%     finalSeg = finalSeg(1);
% end
    
