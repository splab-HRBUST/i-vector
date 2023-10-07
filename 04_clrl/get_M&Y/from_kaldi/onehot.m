function onehotMatrix = onehot(labels)  % labels 包含所有类别标签的类别，1*nums或者nums*1单个类别或者字符向量
    % 产生独热编码矩阵，返回矩阵每列为一个one hot标签
    % 独热编码矩阵，nClasses*nums大小，每列中有且仅有一个1，其余为0，nClasses为类别数量，根据输入labels推算
    
    labels = categorical(labels(:));
    classes = string(categories(labels));% 顺序为onehot矩阵的每列的类别顺序

    nums = numel(labels);  % 数据数
    nClasses = length(classes);  % 标签类别  

    onehotMatrix = [];
    Yn = zeros(nClasses,1);
    for i = 1:nums
        ind = classes==string(labels(i));
        Yn(ind,1) = 1;
        onehotMatrix = [onehotMatrix,Yn];
        Yn = zeros(nClasses,1);
    end
    

