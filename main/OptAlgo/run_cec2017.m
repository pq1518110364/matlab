clear all
close all
clc

%% 设置路径和日志
setup_paths();
tic; %时间计时器启动
start_logging();

%% 设置参数
PD_no = 30; % 搜索代理数量
Max_iter = 1000; % 最大迭代次数
CEC_f = 17;

%% 获取函数信息
[Function_name, F_num] = get_CEC_name(CEC_f);

% 预分配时间存储数组
PSO_time = zeros(F_num, 1);
GWO_time = zeros(F_num, 1);
SCA_time = zeros(F_num, 1);
AOA_time = zeros(F_num, 1);
HHO_time = zeros(F_num, 1);
SFOA_time = zeros(F_num, 1);

BAEO_time = zeros(F_num, 1);
PIMO_time = zeros(F_num, 1);
AGDO_time = zeros(F_num, 1);
AOO_time = zeros(F_num, 1);
AOO_time_v2 = zeros(F_num, 1);

% 主循环
for a = 1:F_num
    % 跳过不需要的函数
    if CEC_f == 17 && a == 2
        continue; % F2已被删除
    end
    
    % if ~ismember(a, [9,13,,14,18,23,25,26])
    %     continue; % 只测试指定函数
    % end

    % if ~ismember(a, [13,14,18,23,25,26])
    %     continue; % 只测试指定函数
    % end
    
    fprintf('以下为F%d函数的数据展示效果：\n', a);
    
    % 获取函数边界和句柄
    f_name = get_F_name(a);
    [LB, UB, Dim, F_obj] = Function_name(f_name);
    
    %% 计算各算法单次运行时间

    % PSO
    tic;
    [PSO_best_score,PSO_best_pos,PSO_cg_curve]=PSO(PD_no,Max_iter,LB,UB,Dim,F_obj);
    PSO_time(a) = toc;

    % GWO
    tic;
    [GWO_best_score,GWO_best_pos,GWO_cg_curve]=GWO(PD_no,Max_iter,LB,UB,Dim,F_obj);
    GWO_time(a) = toc;

    % SCA
    tic;
    [SCA_best_score,SCA_best_pos,SCA_cg_curve]=SCA(PD_no,Max_iter,LB,UB,Dim,F_obj);
    SCA_time(a) = toc;

    % AOA
    tic;
    [AOA_best_score,AOA_best_pos,AOA_cg_curve]=AOA(PD_no,Max_iter,LB,UB,Dim,F_obj);
    AOA_time(a) = toc;

    % SFOA
    tic;
    [SFOA_best_score,SFOA_best_pos,SFOA_cg_curve]=SFOA(PD_no,Max_iter,LB,UB,Dim,F_obj);
    SFOA_time(a) = toc;

    % BAEO
    tic;
    [BAEO_best_score, BAEO_best_pos, BAEO_cg_curve] = BAEO(PD_no, Max_iter, LB, UB, Dim, F_obj);
    BAEO_time(a) = toc;

    % AGDO
    tic;
    [AGDO_best_score, AGDO_best_pos, AGDO_cg_curve] = AGDO(PD_no, Max_iter, LB, UB, Dim, F_obj);
    AGDO_time(a) = toc;

    % PIMO
    tic;
    [PIMO_best_score, PIMO_best_pos, PIMO_cg_curve] = PIMO(PD_no, Max_iter, LB, UB, Dim, F_obj);
    PIMO_time(a) = toc;
    
    % AOO
    tic;
    AOO_fhd = get_CEC_func_str(CEC_f);
    [AOO_best_score, AOO_best_pos, AOO_cg_curve] = AOOv1(AOO_fhd, Dim, PD_no, Max_iter, LB, UB, F_obj,a);
    AOO_time(a) = toc;

    % AOOv2 改进版本
    tic;
    AOO_fhd = get_CEC_func_str(CEC_f);
    [AOO_best_score_v2,AOO_best_pos_v2, AOO_cg_curve_v2 ] = AOOv2(AOO_fhd, Dim, PD_no, Max_iter, LB, UB, F_obj, a);
    AOO_time_v2(a) = toc;
    
    
    %% 绘制进化曲线（代码保持不变）
    CNT = 30;
    k = round(linspace(1, Max_iter, CNT));
    iter = 1:Max_iter;
    
    % figure('Position', [154, 145, 894, 357]);
    figure('Position', [154, 145, 1341, 535]);
    subplot(1, 2, 1);
    func_plot_cec2017(f_name);
    title(f_name + '函数图');
    xlabel('x_1');
    ylabel('x_2');
    zlabel([f_name, '( x_1 , x_2 )']);
    
    subplot(1, 2, 2);

    % 确保中文显示正常
    set(gca, 'FontName', 'SimHei');
    box on;  % 边框增强图表感
    grid on;
    grid minor;  % 细网格提升数据可读性
    
    % 颜色渐变相关参数，根据绘制顺序（越先绘制越淡，后续加深），定义基础色调数组
    baseColors = {[0.3, 0.7, 0.5], [0.2, 0.6, 0.7], [1, 0.6, 0.2], ...
                  [0.4, 0.7, 0.4], [0.8, 0.5, 0.9], [0.3, 0.5, 0.9], ...
                  [0.8, 0.6, 0.8], [1, 0.4, 0.3],[0.6, 0.4, 0.9],[0.9, 0.3, 0.3]}; 
    
    % PSO算法：偏暖深薄荷绿+菱形标记，绘制顺序1，颜色稍淡
    semilogy(iter(k), PSO_cg_curve(k), ...
             'Color', colorAdjust(baseColors{1}, 1), ...    
             'Marker', 'd', ...               
             'LineStyle', '--', ...           
             'LineWidth', 1.2, ...
             'MarkerSize', 6, ...
             'MarkerEdgeColor', colorAdjust(baseColors{1}, 1, 0.1), ...
             'MarkerFaceColor', 'none');
    hold on;
    
    % GWO算法：偏蓝深色调+五边形标记，绘制顺序2，颜色稍淡
    semilogy(iter(k), GWO_cg_curve(k), ...
             'Color', colorAdjust(baseColors{2}, 2), ...    
             'Marker', 'p', ...               
             'LineStyle', ':', ...            
             'LineWidth', 1.2, ...
             'MarkerSize', 6, ...
             'MarkerEdgeColor', colorAdjust(baseColors{2}, 2, 0.1), ...
             'MarkerFaceColor', 'none');
    hold on;
    
    % SCA算法：深橙色+倒三角形标记，绘制顺序3，颜色稍淡
    semilogy(iter(k), SCA_cg_curve(k), ...
             'Color', colorAdjust(baseColors{3}, 3), ...    
             'Marker', 'v', ...               
             'LineStyle', '-', ...            
             'LineWidth', 1.2, ...
             'MarkerSize', 6, ...
             'MarkerEdgeColor', colorAdjust(baseColors{3}, 3, 0.1), ...
             'MarkerFaceColor', 'none');
    hold on;
    
    % AOA算法：深珊瑚红+加号标记，绘制顺序4，颜色稍淡
    semilogy(iter(k), AOA_cg_curve(k), ...
             'Color', colorAdjust(baseColors{4}, 4), ...   
             'Marker', '+', ...               
             'LineStyle', '-.', ...           
             'LineWidth', 1.2, ...
             'MarkerSize', 7, ...
             'MarkerEdgeColor', colorAdjust(baseColors{4}, 4, 0.1));
    hold on;
    
    % SFOA算法：淡紫色调+正方形标记，绘制顺序5，颜色稍淡
    semilogy(iter(k), SFOA_cg_curve(k), ...
             'Color', colorAdjust(baseColors{5}, 5), ...    
             'Marker', 's', ...               
             'LineStyle', '--', ...           
             'LineWidth', 1.2, ...
             'MarkerSize', 6, ...
             'MarkerEdgeColor', colorAdjust(baseColors{5}, 5, 0.1), ...
             'MarkerFaceColor', 'none');
    hold on;
    
    % % AGDO算法：清新靛蓝色+右三角形标记，绘制顺序6，颜色稍淡
    % semilogy(iter(k), AGDO_cg_curve(k), ...
    %          'Color', colorAdjust(baseColors{6}, 6), ...    
    %          'Marker', '>', ...               
    %          'LineStyle', '--', ...           
    %          'LineWidth', 1.2, ...            
    %          'MarkerSize', 6, ...
    %          'MarkerEdgeColor', colorAdjust(baseColors{6}, 6, 0.1), ...  
    %          'MarkerFaceColor', 'none');      
    % hold on;
    % 
    % % BAEO算法：带紫调颜色+左三角形标记，绘制顺序7，颜色稍淡
    % semilogy(iter(k), BAEO_cg_curve(k), ...
    %          'Color', colorAdjust(baseColors{7}, 7), ...    
    %          'Marker', '<', ...               
    %          'LineStyle', '-', ...            
    %          'LineWidth', 1.2, ...
    %          'MarkerSize', 6, ...
    %          'MarkerEdgeColor', colorAdjust(baseColors{7}, 7, 0.1), ...
    %          'MarkerFaceColor', 'none');
    % hold on;
    % 
    % % PIMO算法：橙红色+小x标记，绘制顺序8，颜色按渐变逻辑调整
    % semilogy(iter(k), PIMO_cg_curve(k), ...
    %          'Color', colorAdjust(baseColors{8}, 8), ... % 橙红色基础色调，参与渐变
    %          'Marker', 'x', ...               % 小x标记
    %          'LineStyle', '-.', ...           % 点划线
    %          'LineWidth', 1.2, ...
    %          'MarkerSize', 6, ...
    %          'MarkerEdgeColor', colorAdjust(baseColors{8}, 8, 0.1), ...  % 边缘稍淡
    %          'MarkerFaceColor', 'none');      % 无填充
    % hold on;
    
    % AOO算法：偏暖翠绿色+星号标记，绘制顺序9，颜色稍淡
    semilogy(iter(k), AOO_cg_curve(k), ...
             'Color', colorAdjust(baseColors{9}, 9), ...    
             'Marker', '*', ...               
             'LineStyle', '-', ...            
             'LineWidth', 1.2, ...
             'MarkerSize', 8, ...
             'MarkerEdgeColor', colorAdjust(baseColors{9}, 9, 0.1));
    hold on;
    
    % AOOv2算法：调整颜色边框，绘制顺序10，颜色按渐变加深
    semilogy(iter(k), AOO_cg_curve_v2(k), ...
             'Color', colorAdjust(baseColors{10}, 10), ...    % 靠后绘制，颜色相对深
             'Marker', 'h', ...               
             'LineStyle', '-', ...            
             'LineWidth', 1.5, ...            
             'MarkerSize', 7, ...
             'MarkerEdgeColor', colorAdjust(baseColors{10}, 10, 0.1), ... % 边框不深
             'MarkerFaceColor', 'none');
    hold on;


    grid on;
    title('各算法在'+f_name+'函数的迭代图');
    xlabel('Iteration');
    ylabel('Best fitness so far');
    % legend('BAEO', 'AGDO','AOOv1', 'AOOv2','PIMO');
    % legend('PSO','GWO','SCA','AOA', 'SFOA','AGDO','BAEO','PIMO','AOOv1', 'AOOv2');
    legend('PSO','GWO','SCA','AOA', 'SFOA','AOOv1', 'AOOv2');

    current_time = datetime('now');  % 获取当前时间
    timestamp = string(current_time, 'yyyyMMdd_HH');  % 格式化为：年-月-日_时-分-秒（如20250725_15）
    result_dir = fullfile(pwd, 'result', 'CEC2017', timestamp);
    if ~exist(result_dir, 'dir')
        mkdir(result_dir);
    end
    saveas(gcf, fullfile(result_dir, sprintf('%s_evolution_curve.png', f_name)), 'png');
    saveas(gcf, fullfile(result_dir, sprintf('%s_evolution_curve.eps', f_name)), 'epsc');
    fprintf('图像已保存至: %s\n', result_dir);

    %% 计算30次运行的统计指标（包含时间）

    % PSO
    PSO_best_score_list = zeros(30, 1);
    PSO_time_list = zeros(30, 1);
    for i = 1:30
        tic;
        [PSO_best_score, PSO_best_pos, PSO_cg_curve] = PSO(PD_no, Max_iter, LB, UB, Dim, F_obj);
        PSO_time_list(i) = toc;
        PSO_best_score_list(i) = PSO_best_score;
    end
    PSO_best = min(PSO_best_score_list);
    PSO_mean_val = mean(PSO_best_score_list);
    PSO_std_val = std(PSO_best_score_list);
    PSO_avg_time = mean(PSO_time_list);

    fprintf('以下为PSO的数据展示：\n');
    fprintf('PSO_Best: %.2e  ', PSO_best);
    fprintf('PSO_Mean: %.2e  ', PSO_mean_val);
    fprintf('PSO_STD: %.2e  ', PSO_std_val);
    fprintf('PSO_AvgTime: %.2f秒\n', PSO_avg_time);
    
     % GWO
    GWO_best_score_list = zeros(30, 1);
    GWO_time_list = zeros(30, 1);
    for i = 1:30
        tic;
        [GWO_best_score, GWO_best_pos, GWO_cg_curve] = GWO(PD_no, Max_iter, LB, UB, Dim, F_obj);
        GWO_time_list(i) = toc;
        GWO_best_score_list(i) = GWO_best_score;
    end
    GWO_best = min(GWO_best_score_list);
    GWO_mean_val = mean(GWO_best_score_list);
    GWO_std_val = std(GWO_best_score_list);
    GWO_avg_time = mean(GWO_time_list);

    fprintf('以下为GWO的数据展示：\n');
    fprintf('GWO_Best: %.2e  ', GWO_best);
    fprintf('GWO_Mean: %.2e  ', GWO_mean_val);
    fprintf('GWO_STD: %.2e  ', GWO_std_val);
    fprintf('GWO_AvgTime: %.2f秒\n', GWO_avg_time);

    % SCA
    SCA_best_score_list = zeros(30, 1);
    SCA_time_list = zeros(30, 1);
    for i = 1:30
        tic;
        [SCA_best_score, SCA_best_pos, SCA_cg_curve] = SCA(PD_no, Max_iter, LB, UB, Dim, F_obj);
        SCA_time_list(i) = toc;
        SCA_best_score_list(i) = SCA_best_score;
    end
    SCA_best = min(SCA_best_score_list);
    SCA_mean_val = mean(SCA_best_score_list);
    SCA_std_val = std(SCA_best_score_list);
    SCA_avg_time = mean(SCA_time_list);

    fprintf('以下为SCA的数据展示：\n');
    fprintf('SCA_Best: %.2e  ', SCA_best);
    fprintf('SCA_Mean: %.2e  ', SCA_mean_val);
    fprintf('SCA_STD: %.2e  ', SCA_std_val);
    fprintf('SCA_AvgTime: %.2f秒\n', SCA_avg_time);

    % AOA
    AOA_best_score_list = zeros(30, 1);
    AOA_time_list = zeros(30, 1);
    for i = 1:30
        tic;
        [AOA_best_score, AOA_best_pos, AOA_cg_curve] = AOA(PD_no, Max_iter, LB, UB, Dim, F_obj);
        AOA_time_list(i) = toc;
        AOA_best_score_list(i) = AOA_best_score;
    end
    AOA_best = min(AOA_best_score_list);
    AOA_mean_val = mean(AOA_best_score_list);
    AOA_std_val = std(AOA_best_score_list);
    AOA_avg_time = mean(AOA_time_list);

    fprintf('以下为AOA的数据展示：\n');
    fprintf('AOA_Best: %.2e  ', AOA_best);
    fprintf('AOA_Mean: %.2e  ', AOA_mean_val);
    fprintf('AOA_STD: %.2e  ', AOA_std_val);
    fprintf('AOA_AvgTime: %.2f秒\n', AOA_avg_time);

    % SFOA
    SFOA_best_score_list = zeros(30, 1);
    SFOA_time_list = zeros(30, 1);
    for i = 1:30
        tic;
        [SFOA_best_score, SFOA_best_pos, SFOA_cg_curve] = SFOA(PD_no, Max_iter, LB, UB, Dim, F_obj);
        SFOA_time_list(i) = toc;
        SFOA_best_score_list(i) = SFOA_best_score;
    end
    SFOA_best = min(SFOA_best_score_list);
    SFOA_mean_val = mean(SFOA_best_score_list);
    SFOA_std_val = std(SFOA_best_score_list);
    SFOA_avg_time = mean(SFOA_time_list);

    fprintf('以下为SFOA的数据展示：\n');
    fprintf('SFOA_Best: %.2e  ', SFOA_best);
    fprintf('SFOA_Mean: %.2e  ', SFOA_mean_val);
    fprintf('SFOA_STD: %.2e  ', SFOA_std_val);
    fprintf('SFOA_AvgTime: %.2f秒\n', SFOA_avg_time);

    % BAEO
    BAEO_best_score_list = zeros(30, 1);
    BAEO_time_list = zeros(30, 1);
    for i = 1:30
        tic;
        [BAEO_best_score, BAEO_best_pos, BAEO_cg_curve] = BAEO(PD_no, Max_iter, LB, UB, Dim, F_obj);
        BAEO_time_list(i) = toc;
        BAEO_best_score_list(i) = BAEO_best_score;
    end
    BAEO_best = min(BAEO_best_score_list);
    BAEO_mean_val = mean(BAEO_best_score_list);
    BAEO_std_val = std(BAEO_best_score_list);
    BAEO_avg_time = mean(BAEO_time_list);

    fprintf('以下为BAEO的数据展示：\n');
    fprintf('BAEO_Best: %.2e  ', BAEO_best);
    fprintf('BAEO_Mean: %.2e  ', BAEO_mean_val);
    fprintf('BAEO_STD: %.2e  ', BAEO_std_val);
    fprintf('BAEO_AvgTime: %.2f秒\n', BAEO_avg_time);

    % AGDO
    AGDO_best_score_list = zeros(30, 1);
    AGDO_time_list = zeros(30, 1);
    for i = 1:30
        tic;
        [AGDO_best_score, AGDO_best_pos, AGDO_cg_curve] = AGDO(PD_no, Max_iter, LB, UB, Dim, F_obj);
        AGDO_time_list(i) = toc;
        AGDO_best_score_list(i) = AGDO_best_score;
    end
    AGDO_best = min(AGDO_best_score_list);
    AGDO_mean_val = mean(AGDO_best_score_list);
    AGDO_std_val = std(AGDO_best_score_list);
    AGDO_avg_time = mean(AGDO_time_list);

    fprintf('以下为AGDO的数据展示：\n');
    fprintf('AGDO_Best: %.2e  ', AGDO_best);
    fprintf('AGDO_Mean: %.2e  ', AGDO_mean_val);
    fprintf('AGDO_STD: %.2e  ', AGDO_std_val);
    fprintf('AGDO_AvgTime: %.2f秒\n', AGDO_avg_time);

    % PIMO

    PIMO_best_score_list = zeros(30, 1);
    PIMO_time_list = zeros(30, 1);
    for i = 1:30
        tic;
        [PIMO_best_score, PIMO_best_pos, PIMO_cg_curve] = PIMO(PD_no, Max_iter, LB, UB, Dim, F_obj);
        PIMO_time_list(i) = toc;
        PIMO_best_score_list(i) = PIMO_best_score;
    end
    PIMO_best = min(PIMO_best_score_list);
    PIMO_mean_val = mean(PIMO_best_score_list);
    PIMO_std_val = std(PIMO_best_score_list);
    PIMO_avg_time = mean(PIMO_time_list);

    fprintf('以下为PIMO的数据展示：\n');
    fprintf('PIMO_Best: %.2e  ', PIMO_best);
    fprintf('PIMO_Mean: %.2e  ', PIMO_mean_val);
    fprintf('PIMO_STD: %.2e  ', PIMO_std_val);
    fprintf('PIMO_AvgTime: %.2f秒\n', PIMO_avg_time);

    % AOO
    AOO_best_score_list = zeros(30, 1);
    AOO_time_list = zeros(30, 1);
    for i = 1:30
        tic;
        [AOO_best_score, AOO_best_pos, AOO_cg_curve] = AOOv1(AOO_fhd, Dim, PD_no, Max_iter, LB, UB,F_obj, a);
        AOO_time_list(i) = toc;
        AOO_best_score_list(i) = AOO_best_score;
    end
    AOO_best = min(AOO_best_score_list);
    AOO_mean_val = mean(AOO_best_score_list);
    AOO_std_val = std(AOO_best_score_list);
    AOO_avg_time = mean(AOO_time_list);

    fprintf('以下为AOO的数据展示：\n');
    fprintf('AOO_Best: %.2e  ', AOO_best);
    fprintf('AOO_Mean: %.2e  ', AOO_mean_val);
    fprintf('AOO_STD: %.2e  ', AOO_std_val);
    fprintf('AOO_AvgTime: %.2f秒\n', AOO_avg_time);

     % AOO
    AOO_best_score_v2_list = zeros(30, 1);
    AOO_time_v2_list = zeros(30, 1);
    for i = 1:30
        tic;
        [AOO_best_score_v2, AOO_best_pos_v2, AOO_cg_curve_v2] = AOOv2(AOO_fhd, Dim, PD_no, Max_iter, LB, UB,F_obj, a);
        AOO_time_v2_list(i) = toc;
        AOO_best_score_v2_list(i) = AOO_best_score_v2;
    end
    AOO_best_v2 = min(AOO_best_score_v2_list);
    AOO_mean_val_v2 = mean(AOO_best_score_v2_list);
    AOO_std_val_v2 = std(AOO_best_score_v2_list);
    AOO_avg_time_v2 = mean(AOO_time_v2_list);

    fprintf('以下为AOO_v2的数据展示：\n');
    fprintf('AOO_Best_v2: %.2e  ', AOO_best_v2);
    fprintf('AOO_Mean_v2: %.2e  ', AOO_mean_val_v2);
    fprintf('AOO_STD_v2: %.2e  ', AOO_std_val_v2);
    fprintf('AOO_AvgTime_v2: %.2f秒\n', AOO_avg_time_v2);
    
    fprintf('\n');
end

%% 输出各算法在不同函数上的总运行时间
% fprintf('===== 各算法在测试函数上的总运行时间 =====\n');
% for a = 1:F_num
%     if ~ismember(a, [12,30])
%         continue;
%     end
%     fprintf('F%d: BAEO=%.2fs, AOO=%.2fs, AOOv2=%.2fs, PIMO=%.2fs\n', ...
%         a, BAEO_time(a), AOO_time(a), AOO_time_v2(a),PIMO_time(a));
% end

%% 结束日志
stop_logging();
toc;

function adjustedColor = colorAdjust(color, index, edgeOffset)
    if nargin < 3
        edgeOffset = 0;
    end
    % 简单的调整逻辑，index越大（绘制越靠后），颜色数值稍大（更浓），这里乘以一个大于1的系数，系数可微调
    adjustFactor = 1 + 0.05 * index; 
    adjustedColor = min(color .* adjustFactor, 1); % 防止超过1
    if edgeOffset ~= 0
        adjustedColor = adjustedColor - edgeOffset; % 边缘颜色稍淡
        adjustedColor(adjustedColor < 0) = 0; % 防止小于0
    end
end