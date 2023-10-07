function [eer, dcf08, dcf10, dcf_vox, dcf_king] = compute_eer(scores, labels, showfig)
    % calculates the equal error rate (EER) performance measure.  计算等错误率（EER）性能度量
    %
    % Inputs:                      
    %   - scores        : likelihood scores for target and non-target trials  目标和非目标试验的可能性得分
    %   - labels        : true labels for target and non-target trials, can be 目标和非目标试验的真实标签，可以是
    %					  either binary (0's and 1's) or a cell array with 二进制（0和1）或具有
    %					  "target" and "impostor" string labels “目标”和“冒名顶替者”字符串标签
    %   - showfig       : if true the DET curve is displayed  如果为真，则显示DET曲线
    %
    % Outputs:
    %   - eer           : percent equal error rate (EER)  相等错误率百分比（eer）
    %   - dcf08         : minimum detection cost function (DCF) with SRE'08 parameters  SRE'08的最小检测成本函数（DCF）
    %   - dcf10         : minimum DCF with SRE'10 parameters  具有SRE'10参数的最小DCF
    %
    %
    % Omid Sadjadi <s.omid.sadjadi@gmail.com>
    % Microsoft Research, Conversational Systems Research Center


    if iscell(labels),  %判断是否为结构体 -- 不需要 
        labs = zeros(length(labels), 1);
        labs(ismember(labels, 'target')) = 1;
        labels = labs; 
        clear labs;
    end

    [~,I] = sort(scores);  % 对得分进行排序，返回索引向量的集合
    x = labels(I);  %得分升序的索引向量集合

    % cumsun() -- 计算一个数组各行的累加值，函数用法是B = cumsum(A,dim)，或B = cumsum(A)
    FN = cumsum( x == 1 ) / (sum( x == 1 ) + eps);  % 标签向量为1  eps-精度，用于取整
    TN = cumsum( x == 0 ) / (sum( x == 0 ) + eps);  % 标签向量为0
    FP = 1 - TN;
    TP = 1 - FN;

    FNR = FN ./ ( TP + FN + eps );
    FPR = FP ./ ( TN + FP + eps );
    difs = FNR - FPR;
    idx1 = find(difs< 0, 1, 'last');
    idx2 = find(difs>= 0, 1 );
    x = [FNR(idx1); FPR(idx1)];
    y = [FNR(idx2); FPR(idx2)];
    a = ( x(1) - x(2) ) / ( y(2) - x(2) - y(1) + x(1) );
    eer = 100 * ( x(1) + a * ( y(1) - x(1) ) );  % 等错误率
    
    % =====================================================================

    if ( nargout > 1 ),
        Cmiss = 10; Cfa = 1; P_tgt = 0.01; % SRE-2008 performance parameters
        Cdet  = Cmiss * FNR * P_tgt + Cfa * FPR * ( 1 - P_tgt);
    %     Cdefault = min(Cmiss * P_tgt, Cfa * ( 1 - P_tgt));
        dcf08 = 100 * min(Cdet); % note this is not percent
    end
    if ( nargout > 3 ),
        Cmiss = 1; Cfa = 1; P_tgt = 0.001; % SRE-2010 performance parameters
        Cdet  = Cmiss * FNR * P_tgt + Cfa * FPR * ( 1 - P_tgt);
    %     Cdefault = min(Cmiss * P_tgt, Cfa * ( 1 - P_tgt));
        dcf10 = 100 * min(Cdet); % note this is not percent
    end
    if ( nargout >= 4 ),
        Cmiss = 1; Cfa = 1; P_tgt = 0.01; % VoxCeleb performance parameters
        Cdet  = Cmiss * FNR * P_tgt + Cfa * FPR * ( 1 - P_tgt);
    %     Cdefault = min(Cmiss * P_tgt, Cfa * ( 1 - P_tgt));
        dcf_vox = 100 * min(Cdet); % note this is not percent
    end

    if ( nargout == 5 ),
        Cmiss = 1; Cfa = 1; P_tgt = 0.02; % King-ASR-010 performance parameters
        Cdet  = Cmiss * FNR * P_tgt + Cfa * FPR * ( 1 - P_tgt);
    %     Cdefault = min(Cmiss * P_tgt, Cfa * ( 1 - P_tgt));
        dcf_king = 100 * min(Cdet); % note this is not percent  注意这不是百分比
    end

    if showfig
    %     figure
        plot_det(FPR, FNR)
    end

function plot_det(FPR, FNR)
    % plots the detection error tradeoff (DET) curve

    fnr = icdf(FNR);
    fpr = icdf(FPR);
    plot(fpr, fnr);

    xtick = [0.0001, 0.0002, 0.0005, 0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.4]; % default
    xticklabel = num2str(xtick * 100, '%g\n');
    xticklabel = textscan(xticklabel, '%s'); xticklabel = xticklabel{1};
    set (gca, 'xtick', icdf(xtick));
    set (gca, 'xticklabel', xticklabel);
    xlim(icdf([0.005 0.4])); 
    xlabel ('False Positive Rate (FPR) [%]');

    ytick = xtick;         
    yticklabel = num2str(ytick * 100, '%g\n');
    yticklabel = textscan(yticklabel, '%s'); yticklabel = yticklabel{1};
    set (gca, 'ytick', icdf(ytick));
    set (gca, 'yticklabel', yticklabel);
    ylim(icdf([0.005 0.4]));  
    ylabel ('False Negative Rate (FNR) [%]')

    grid on;
    box on;
    axis square;
    axis manual;

function y = icdf(x)
    % computes the inverse of cumulative distribution function in x
    y = -sqrt(2).*erfcinv(2 * ( x + eps));
    y(isinf(y)) = nan;
