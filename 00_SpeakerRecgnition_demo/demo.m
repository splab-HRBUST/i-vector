%% 加载模型
fprintf('加载模型...\n');
addpath('Toolkit');
ubm           = importdata('Model/ubm.mat');
patameters_FA = importdata('Model/patameters_FA_baseline.mat');
spkInfo       = importdata('Model/spkInfo.mat');

%% 说话人注册
count = size(spkInfo,2);
prompt = '[1]显示注册名单\n[2]说话人注册\n[3]说话人匹配\n[4]退出\n';
fprintf('请选择功能：\n');
strFunc = input(prompt,'s');
while strFunc ~= '4'
    switch strFunc
        case '1'   % 显示注册名单
            if count == 0
                fprintf('注册表为空！\n\n');
            else
                nameList = [];
                for i = 1 : size(spkInfo,2)
                    nameList = [nameList,' ', spkInfo{i}.spkName];
                end
                fprintf('已注册说话人名单为：\n');  
                disp(nameList);
                fprintf('\n');               
                waitingEnter = input(['按回车继续执行...',char(10)]);
                if isempty(waitingEnter) == 1
                    fprintf('——————————————————————————————————\n');
                    fprintf('请选择功能：\n');
                    prompt = '[1]显示注册名单\n[2]说话人注册\n[3]说话人匹配\n[4]退出\n';
                    strFunc = input(prompt,'s');
                end
            end
           
        case '2'  % 说话人注册
            count = count + 1;
            spkInfo{count} = spkEnroll(ubm, patameters_FA);
            waitingEnter = input(['按回车继续执行...',char(10)]);
            if isempty(waitingEnter) == 1
                fprintf('——————————————————————————————————\n');
                fprintf('请选择功能：\n');
                prompt = '[1]显示注册名单\n[2]说话人注册\n[3]说话人匹配\n[4]退出\n';
                strFunc = input(prompt,'s');
            end
            
        case '3'  % 说话人匹配
            spkRecog(ubm, patameters_FA, spkInfo);
            waitingEnter = input(['按回车继续执行...',char(10)]);
            if isempty(waitingEnter) == 1
                fprintf('——————————————————————————————————\n');
                fprintf('请选择功能：\n');
                prompt = '[1]显示注册名单\n[2]说话人注册\n[3]说话人匹配\n[4]退出\n';
                strFunc = input(prompt,'s');
            end           
            
        otherwise   %
            disp('请重新输入【1-4】：')
            strFunc = input(prompt,'s');    
    end
end

clearvars -except ubm patameters_FA spkInfo
save('Model\spkInfo.mat','spkInfo','-v7.3');
fprintf('程序关闭。\n');




