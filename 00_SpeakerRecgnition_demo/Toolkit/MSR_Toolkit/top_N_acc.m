function acc_iv_cos_topN = top_N_acc(cos,label_test,N)
% cos =1 - pdist2(IV_test(:,1:size(label_test,2))',IV_model','cosine');

% [~,prelabel_test_cos] = max(cos',[],1);
% acc_iv_cos = sum(prelabel_test_cos == label_test(:)')/size(label_test(:),1);
%%%%%N = 5;

if size(label_test,1) ~= 1
    label_test =label_test';
end

[ ~, ix ] = sort( cos', 'descend' );

    for n = 1 : N
        [~,~,v] = find(ix(1:n,:) == repmat(label_test,n,1));
        acc_iv_cos_topN(n) = sum(v)/size(label_test(:),1)*100;
    end
end
