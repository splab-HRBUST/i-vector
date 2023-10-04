function mfc = extract_mfcc(Files,con_cnvn)
% function mfc = extract_mfcc(x,fs)
%------------------------------------------
%input
%       files path
%------------------------------------------
if ischar(Files) || iscellstr(Files)
 	filesPath = Files;
 end
% [x, fs]=audioread(Files);
   [x,fs]=vad(filesPath,'a');

   frameLength = 256; % 256 for King; 512 for Vox 
   p=24; % �˲�������
   
   bank=melbankm(p,frameLength,fs,0,0.5,'t'); %Mel�˲����ĸ���Ϊp��FFT�任�ĳ���ΪframeLength������Ƶ��Ϊ fs Hz
	

    %��һ��Mel�˲�����ϵ�� 
    bank=full(bank);
    bank=bank/max(bank(:));

    % DCTϵ��,20*24
    dctnum = 20;
    for k=1:dctnum
      n=0:23;
      dctcoef(k,:)=cos((2*n+1)*k*pi/(2*p));
    end
    
    w = 1 + 6 * sin(pi * [1:dctnum] ./ dctnum); %��һ��������������
    w = w/max(w); %Ԥ�����˲���

    for nChannel = 1 : size(x,2)
        % Pre-emphasis filter, preemphasis coefficient is 0.9375
        xx=double(x(:,nChannel));
        xx=filter([1 -0.9375],1,xx);

        % After the speech signal frame, frame every behavior of a frame
        frameShift = 80;
        xx=enframe(xx,frameLength,frameShift);%In xx, 256 points into a frame, frame shift is 80/256

        % Calculate the MFCC parameters per frame
        for i=1:size(xx,1)
          y = xx(i,:);
          s = y'* hamming(frameLength);
          t = abs(fft(s));
          t = t.^2;
          c1=dctcoef * log(bank * t(1:frameLength/2+1)');
          c2 = c1.*w';
          m(i,:)=c2;
        end
%         m = rmmissing(m);
        m( all(isnan(m),2),:) = [];

        %Calculate the first-order differential coefficient
        dtm=zeros(size(m));
        for l=3:size(m,1)-2
            dtm(l,:)=-2*m(l-2,:)-m(l-1,:)+m(l+1,:)+2*m(l+2,:);
        end
        dtm=dtm/3;
        %Calculate the second-order differential coefficient
        dtmm=zeros(size(dtm));
        for p=3:size(dtm,1)-2
            dtmm(p,:)=-2*dtm(p-2,:)-dtm(p-1,:)+dtm(p+1,:)+2*dtm(p+2,:);
        end
        dtmm=dtmm/3;
        %Combination of MFCC parameters and first-order differential MFCC parameters
        ccc=[m dtm dtmm];
        %Get rid of the fore and aft two frames, because the two frames of first order difference parameter is 0
        ccc=ccc(3:size(m,1)-2,:);
        
        if con_cnvn == 1   
            Fea = cmvn(ccc', 1);
        else
            Fea = ccc';
        end
        
        mfc{nChannel}=Fea;
    end
end
