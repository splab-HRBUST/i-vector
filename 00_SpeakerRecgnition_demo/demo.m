%% ����ģ��
fprintf('����ģ��...\n');
addpath('Toolkit');
ubm           = importdata('Model/ubm.mat');
patameters_FA = importdata('Model/patameters_FA_baseline.mat');
spkInfo       = importdata('Model/spkInfo.mat');

%% ˵����ע��
count = size(spkInfo,2);
prompt = '[1]��ʾע������\n[2]˵����ע��\n[3]˵����ƥ��\n[4]�˳�\n';
fprintf('��ѡ���ܣ�\n');
strFunc = input(prompt,'s');
while strFunc ~= '4'
    switch strFunc
        case '1'   % ��ʾע������
            if count == 0
                fprintf('ע���Ϊ�գ�\n\n');
            else
                nameList = [];
                for i = 1 : size(spkInfo,2)
                    nameList = [nameList,' ', spkInfo{i}.spkName];
                end
                fprintf('��ע��˵��������Ϊ��\n');  
                disp(nameList);
                fprintf('\n');               
                waitingEnter = input(['���س�����ִ��...',char(10)]);
                if isempty(waitingEnter) == 1
                    fprintf('��������������������������������������������������������������������\n');
                    fprintf('��ѡ���ܣ�\n');
                    prompt = '[1]��ʾע������\n[2]˵����ע��\n[3]˵����ƥ��\n[4]�˳�\n';
                    strFunc = input(prompt,'s');
                end
            end
           
        case '2'  % ˵����ע��
            count = count + 1;
            spkInfo{count} = spkEnroll(ubm, patameters_FA);
            waitingEnter = input(['���س�����ִ��...',char(10)]);
            if isempty(waitingEnter) == 1
                fprintf('��������������������������������������������������������������������\n');
                fprintf('��ѡ���ܣ�\n');
                prompt = '[1]��ʾע������\n[2]˵����ע��\n[3]˵����ƥ��\n[4]�˳�\n';
                strFunc = input(prompt,'s');
            end
            
        case '3'  % ˵����ƥ��
            spkRecog(ubm, patameters_FA, spkInfo);
            waitingEnter = input(['���س�����ִ��...',char(10)]);
            if isempty(waitingEnter) == 1
                fprintf('��������������������������������������������������������������������\n');
                fprintf('��ѡ���ܣ�\n');
                prompt = '[1]��ʾע������\n[2]˵����ע��\n[3]˵����ƥ��\n[4]�˳�\n';
                strFunc = input(prompt,'s');
            end           
            
        otherwise   %
            disp('���������롾1-4����')
            strFunc = input(prompt,'s');    
    end
end

clearvars -except ubm patameters_FA spkInfo
save('Model\spkInfo.mat','spkInfo','-v7.3');
fprintf('����رա�\n');




