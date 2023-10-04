function [M_CV, label_CV] = cv_vox(M_dev,spk_cell_dev,label_dev)

M1=[];M2=[];M3=[];M4=[];M5=[];
label_dev_num = zeros(1,1211);

for nspk = 1 : size(spk_cell_dev,1)
    
    label_dev_num(nspk) = size(spk_cell_dev{nspk},2);
    spk_ind = sum(label_dev_num(1:nspk-1));
    
    num_cv_14(nspk) = floor(label_dev_num(nspk)/5);
    num_cv_5(nspk) = label_dev_num(nspk)-num_cv_14(nspk)*4;
    
    spk_cv_ind_14 = sum(num_cv_14(1:nspk-1));
    spk_cv_ind_5 = sum(num_cv_5(1:nspk-1));
      
    
    M1(:,spk_cv_ind_14+1:spk_cv_ind_14+num_cv_14(nspk)) = M_dev(:,spk_ind+1:spk_ind+num_cv_14(nspk));
    M2(:,spk_cv_ind_14+1:spk_cv_ind_14+num_cv_14(nspk)) = M_dev(:,spk_ind+num_cv_14(nspk)+1:spk_ind+num_cv_14(nspk)*2);
    M3(:,spk_cv_ind_14+1:spk_cv_ind_14+num_cv_14(nspk)) = M_dev(:,spk_ind+num_cv_14(nspk)*2+1:spk_ind+num_cv_14(nspk)*3);
    M4(:,spk_cv_ind_14+1:spk_cv_ind_14+num_cv_14(nspk)) = M_dev(:,spk_ind+num_cv_14(nspk)*3+1:spk_ind+num_cv_14(nspk)*4);
    M5(:,spk_cv_ind_5+1 :spk_cv_ind_5+num_cv_5(nspk))   = M_dev(:,spk_ind+num_cv_14(nspk)*4+1:spk_ind+num_cv_14(nspk)*4+num_cv_5(nspk));
    
    label1(:,spk_cv_ind_14+1:spk_cv_ind_14+num_cv_14(nspk)) = label_dev(:,spk_ind+1:spk_ind+num_cv_14(nspk));
    label2(:,spk_cv_ind_14+1:spk_cv_ind_14+num_cv_14(nspk)) = label_dev(:,spk_ind+num_cv_14(nspk)+1:spk_ind+num_cv_14(nspk)*2);
    label3(:,spk_cv_ind_14+1:spk_cv_ind_14+num_cv_14(nspk)) = label_dev(:,spk_ind+num_cv_14(nspk)*2+1:spk_ind+num_cv_14(nspk)*3);
    label4(:,spk_cv_ind_14+1:spk_cv_ind_14+num_cv_14(nspk)) = label_dev(:,spk_ind+num_cv_14(nspk)*3+1:spk_ind+num_cv_14(nspk)*4);
    label5(:,spk_cv_ind_5+ 1:spk_cv_ind_5+num_cv_5(nspk))   = label_dev(:,spk_ind+num_cv_14(nspk)*4+1:spk_ind+num_cv_14(nspk)*4+num_cv_5(nspk));    

end

M_CV{1} = M1;
M_CV{2} = M2;
M_CV{3} = M3;
M_CV{4} = M4;
M_CV{5} = M5;

label_CV{1} = label1;
label_CV{2} = label2;
label_CV{3} = label3;
label_CV{4} = label4;
label_CV{5} = label5;