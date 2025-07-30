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
BAEO_time = zeros(F_num, 1);
AOO_time = zeros(F_num, 1);
AOO_time_v2 = zeros(F_num, 1);
PIMO_time = zeros(F_num, 1);

% 主循环
for a = 1:F_num
    % 跳过不需要的函数
    if CEC_f == 17 && a == 2
        continue; % F2已被删除
    end
    
    if ~ismember(a, [1, 8,12,19,30])
        continue; % 只测试指定函数
    end
    
    fprintf('以下为F%d函数的数据展示效果：\n', a);
    
    % 获取函数边界和句柄
    f_name = get_F_name(a);
    [LB, UB, Dim, F_obj] = Function_name(f_name);
    
    %% 计算各算法单次运行时间
    % BAEO
    tic;
    [BAEOBest_score, BAEOBest_pos, BAEO_cg_curve] = BAEO(PD_no, Max_iter, LB, UB, Dim, F_obj);
    BAEO_time(a) = toc;
    
    % AOO
    AOO_fhd = get_CEC_func_str(CEC_f);
    tic;
    [AOOBest_score, AOOBest_pos, AOO_cg_curve] = AOOv1(AOO_fhd, Dim, PD_no, Max_iter, LB, UB, F_obj,a);
    AOO_time(a) = toc;

    % AOO
    AOO_fhd = get_CEC_func_str(CEC_f);
    tic;
    [AOOBest_score_v2, ~, AOO_cg_curve_v2, AOO_v2_diversity, AOO_v2_threshold, AOO_v2_c_curve] = AOOv2(AOO_fhd, Dim, PD_no, Max_iter, LB, UB, F_obj, a);
    AOO_time_v2(a) = toc;
    
    % PIMO
    tic;
    [PIMOBest_score, PIMOBest_pos, PIMO_cg_curve] = PIMO(PD_no, Max_iter, LB, UB, Dim, F_obj);
    PIMO_time(a) = toc;
    
    %% 绘制进化曲线（代码保持不变）
    CNT = 30;
    k = round(linspace(1, Max_iter, CNT));
    iter = 1:Max_iter;
    
    figure('Position', [154, 145, 894, 357]);
    subplot(1, 2, 1);
    func_plot_cec2017(f_name);
    title(f_name + '函数图');
    xlabel('x_1');
    ylabel('x_2');
    zlabel([f_name, '( x_1 , x_2 )']);
    
    subplot(1, 2, 2);
    semilogy(iter(k), BAEO_cg_curve(k), 'Color', [1 0.5 0], 'Marker', '+', 'LineStyle', '-.', 'linewidth', 1);
    hold on;
    semilogy(iter(k), AOO_cg_curve(k), 'b-*', 'linewidth', 1);
    hold on;
    semilogy(iter(k), AOO_cg_curve_v2(k), 'Color', [0.6, 0.1, 0.8], 'Marker', 'h', 'LineStyle', '-', 'LineWidth', 1, 'MarkerSize', 6);  % 深紫色实线+六边形标记
    hold on;
    semilogy(iter(k), PIMO_cg_curve(k), 'Color', [0.6, 0.3, 0.1], 'Marker', '>', 'LineStyle', '-.', 'LineWidth', 1, 'MarkerSize', 6);
    hold on;
    
    grid on;
    title('各算法在'+f_name+'函数的迭代图');
    xlabel('Iteration');
    ylabel('Best fitness so far');
    legend('BAEO', 'AOOv1', 'AOOv2','PIMO');
    
     current_time = datetime('now');  % 获取当前时间
    timestamp = string(current_time, 'yyyyMMdd_HH');  % 格式化为：年-月-日_时-分-秒（如20250725_15）
    result_dir = fullfile(pwd, 'result', 'AOO', timestamp);
    if ~exist(result_dir, 'dir')
        mkdir(result_dir);
    end
    saveas(gcf, fullfile(result_dir, sprintf('%s_evolution_curve.png', f_name)), 'png');
    saveas(gcf, fullfile(result_dir, sprintf('%s_evolution_curve.eps', f_name)), 'epsc');
    fprintf('图像已保存至: %s\n', result_dir);
    
     %% 绘制 AOOv2 的详细分析曲线
    figure('Position', [100, 100, 1200, 600]); % 更大的图窗以容纳更多子图
    
    % 子图 1: AOOv2 收敛曲线
    subplot(2, 2, 1);
    semilogy(iter, AOO_cg_curve_v2, 'Color', [0.6, 0.1, 0.8], 'LineWidth', 1.5);
    title('AOOv2 在 ' + f_name + ' 函数的收敛曲线');
    xlabel('迭代次数');
    ylabel('最佳适应度');
    grid on;
    
    % 子图 2: AOOv2 种群多样性曲线
    subplot(2, 2, 2);
    plot(iter, AOO_v2_diversity, 'Color', [0.2, 0.8, 0.2], 'LineWidth', 1.5); % 绿色
    title('AOOv2 种群多样性 (' + f_name + ')');
    xlabel('迭代次数');
    ylabel('平均个体距离');
    grid on;
    
    % 子图 3: AOOv2 探索/利用阈值曲线
    subplot(2, 2, 3);
    plot(iter, AOO_v2_threshold, 'Color', [0.8, 0.2, 0.2], 'LineWidth', 1.5); % 红色
    title('AOOv2 探索/利用阈值 (' + f_name + ')');
    xlabel('迭代次数');
    ylabel('阈值');
    grid on;
    
    % 子图 4: AOOv2 参数 c 衰减曲线
    subplot(2, 2, 4);
    plot(iter, AOO_v2_c_curve, 'Color', [0.2, 0.2, 0.8], 'LineWidth', 1.5); % 蓝色
    title('AOOv2 参数 c 衰减 (' + f_name + ')');
    xlabel('迭代次数');
    ylabel('c 值');
    grid on;
    
    % 保存 AOOv2 详细分析图
    saveas(gcf, fullfile(result_dir, sprintf('%s_AOOv2_detailed_analysis.png', f_name)), 'png');
    saveas(gcf, fullfile(result_dir, sprintf('%s_AOOv2_detailed_analysis.eps', f_name)), 'epsc');
    fprintf('AOOv2 详细分析图像已保存至: %s\n', result_dir);
    close(gcf); % 关闭当前图形


    %% 计算30次运行的统计指标（包含时间）
    % BAEO
    BAEO_best_score_list = zeros(30, 1);
    BAEO_time_list = zeros(30, 1);
    for i = 1:30
        tic;
        [BAEOBest_score, BAEOBest_pos, BAEO_cg_curve] = BAEO(PD_no, Max_iter, LB, UB, Dim, F_obj);
        BAEO_time_list(i) = toc;
        BAEO_best_score_list(i) = BAEOBest_score;
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
    
    % AOO
    AOO_best_score_list = zeros(30, 1);
    AOO_time_list = zeros(30, 1);
    for i = 1:30
        tic;
        [AOOBest_score, AOOBest_pos, AOO_cg_curve] = AOOv1(AOO_fhd, Dim, PD_no, Max_iter, LB, UB,F_obj, a);
        AOO_time_list(i) = toc;
        AOO_best_score_list(i) = AOOBest_score;
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
        [AOOBest_score_v2, AOOBest_pos_v2, AOO_cg_curve_v2] = AOOv2(AOO_fhd, Dim, PD_no, Max_iter, LB, UB,F_obj, a);
        AOO_time_v2_list(i) = toc;
        AOO_best_score_v2_list(i) = AOOBest_score_v2;
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
    
    %PIMO
    PIMO_best_score_list = zeros(30, 1);
    PIMO_time_list = zeros(30, 1);
    for i = 1:30
        tic;
        [PIMOBest_score, PIMOBest_pos, PIMO_cg_curve] = PIMO(PD_no, Max_iter, LB, UB, Dim, F_obj);
        PIMO_time_list(i) = toc;
        PIMO_best_score_list(i) = PIMOBest_score;
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
    
    fprintf('\n');
end

%% 输出各算法在不同函数上的总运行时间
fprintf('===== 各算法在测试函数上的总运行时间 =====\n');
for a = 1:F_num
    if ~ismember(a, [1, 8,12,19,30])
        continue;
    end
    fprintf('F%d: BAEO=%.2fs, AOO=%.2fs, AOOv2=%.2fs, PIMO=%.2fs\n', ...
        a, BAEO_time(a), AOO_time(a), AOO_time_v2(a),PIMO_time(a));
end

%% 结束日志
stop_logging();
toc;