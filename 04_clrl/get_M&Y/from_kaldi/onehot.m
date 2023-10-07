function onehotMatrix = onehot(labels)  % labels ������������ǩ�����1*nums����nums*1�����������ַ�����
    % �������ȱ�����󣬷��ؾ���ÿ��Ϊһ��one hot��ǩ
    % ���ȱ������nClasses*nums��С��ÿ�������ҽ���һ��1������Ϊ0��nClassesΪ�����������������labels����
    
    labels = categorical(labels(:));
    classes = string(categories(labels));% ˳��Ϊonehot�����ÿ�е����˳��

    nums = numel(labels);  % ������
    nClasses = length(classes);  % ��ǩ���  

    onehotMatrix = [];
    Yn = zeros(nClasses,1);
    for i = 1:nums
        ind = classes==string(labels(i));
        Yn(ind,1) = 1;
        onehotMatrix = [onehotMatrix,Yn];
        Yn = zeros(nClasses,1);
    end
    

