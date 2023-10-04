function [M_CV, label_CV] = cv_king(M,nSpk,nCV)


M1=[];M2=[];M3=[];M4=[];M5=[];
label = [];

nGroup = nSpk/nCV;

for nspk = 1 : nSpk

    M1(:,(nspk-1)*nGroup+1:nspk*nGroup) = M(:,(nspk-1)*120+1:(nspk-1)*120+nGroup);
    M2(:,(nspk-1)*nGroup+1:nspk*nGroup) = M(:,(nspk-1)*120+nGroup+1:(nspk-1)*120+nGroup*2);
    M3(:,(nspk-1)*nGroup+1:nspk*nGroup) = M(:,(nspk-1)*120+nGroup*2+1:(nspk-1)*120+nGroup*3);
    M4(:,(nspk-1)*nGroup+1:nspk*nGroup) = M(:,(nspk-1)*120+nGroup*3+1:(nspk-1)*120+nGroup*4);
    M5(:,(nspk-1)*nGroup+1:nspk*nGroup) = M(:,(nspk-1)*120+nGroup*4+1:(nspk-1)*120+nGroup*5);

    label = [label nspk*ones(1,24)];
end


label_CV = [label label label label label];

M_CV{1} = M1;
M_CV{2} = M2;
M_CV{3} = M3;
M_CV{4} = M4;
M_CV{5} = M5;