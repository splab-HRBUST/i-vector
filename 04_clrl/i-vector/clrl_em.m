function [patameters_PPLS] = clrl_em(M_dev,Y_dev,num,patameters_PPLS0)

[num.SVdim, num.samples] = size(M_dev);  
num.spk = size(Y_dev,1);  

if nargin > 3  % ��������������������£�
    I = patameters_PPLS0.I;
    m = patameters_PPLS0.m;
    mu_y = patameters_PPLS0.mu_y;
    T = patameters_PPLS0.T;
    Q = patameters_PPLS0.Q;
    deltax = patameters_PPLS0.deltax;
    deltay = patameters_PPLS0.deltay;
else  % �������Ϊ3����һ�Σ�-- ��ʼ��
    % fprintf('Randomly initializing the PPLS hyperparameters ...\n');
    I = eye(num.IVdim);  % num.IVdim--�����ӿռ�ά��
    T = M_dev(:,1:num.IVdim);  %�ܱ仯�ռ����
    Q = Y_dev(:,1:num.IVdim);  %���ɾ���
    deltax = ones(num.SVdim,1); 
    deltay = ones(num.spk,1); 
    m = nanmean(M_dev,2);  % GMM���������е�ƫ�ã���--ƽ����
    mu_y = nanmean(Y_dev,2);  % ����ǩ���������е�ƫ��
end

centeredM = bsxfun(@minus, M_dev, m);  
varianceM = sum(centeredM.^2,2) / num.samples; 
 
centeredY = bsxfun(@minus, Y_dev, mu_y);  
varianceY = sum(centeredY.^2,2) / num.samples;

fprintf('Train the PPLS-based TVM with dimension %d:\n',num.IVdim);
for niter = 1 : num.nIters
    
    fprintf('EM iter#: %d \t', niter);
    tim1 = tic;  % ����ִ��ʱ��
    
    % E step -- ������
    B = bsxfun(@rdivide, T', deltax');  % T'./deltax' 
    C = bsxfun(@rdivide, Q', deltay');  % Q'./deltay' 
    L = I + B*T + C*Q;  % pinv(L)--Wn�ĺ���Э�������cov 
    Ex = pinv(L)*(B*centeredM+C*centeredY);  % Wn�ĺ����ֵʸ�� 
    Exx = Ex*Ex' + pinv(L)*num.samples;  % Wn�ĺ�����ؾ��� 
    
    % M step -- ���
    alpha = num.SVdim./sum(T.*T);  % ���� -- 1*400 
    alpha = length_norm(alpha);
    
    x = mean(deltax);
    y = mean(deltay);
    
    T = (centeredM*Ex')*pinv(Exx + x * diag(alpha));  % �����ܱ仯�ռ���� 
    Q = (centeredY*Ex')*pinv(Exx);  %���¸��ɾ��� 

    deltax = varianceM - sum(T .* (T*Exx),2) / num.samples;  % Э������� 
    deltay = varianceY - sum(Q .* (Q*Exx),2) / num.samples;  % Э������� 
    
    flag_x = ismember(double(0),deltax);
    if flag_x >= 1
        deltax_num = size(deltax,1);
        mmin = deltax(1,1);
        for i = 1:deltax_num
            a = deltax(i,1);
            if a~= 0
                mmin = min(mmin,a);
            end
        end
        deltax(deltax(:,1)==0,1) = mmin;
    end
       
    tim1 = toc(tim1);
    fprintf('[elaps = %.2f s]\n', tim1);
end

B = bsxfun(@rdivide, T', deltax');
C = bsxfun(@rdivide, Q', deltay');
L = I + B*T + C*Q;

patameters_PPLS.I = I;  % ��λ��
patameters_PPLS.B = B;
patameters_PPLS.C = C;
patameters_PPLS.L = L;

patameters_PPLS.m = m;  % ƫ��
patameters_PPLS.mu_y = mu_y;

patameters_PPLS.T = T;  % �ܱ仯�ռ����
patameters_PPLS.Q = Q;  % ���ɾ���
patameters_PPLS.deltax = deltax;  % ���²���
patameters_PPLS.deltay = deltay;


