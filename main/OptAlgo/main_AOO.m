clc
clear
close all

%% 将当前工作目录以及其子目录添加到 MATLAB 搜索路径中，以便可以找到相关函数和文件
setup_paths();

%% 启动日志，diary以及计时
tic %时间计时器启动
start_logging();

%% 不实野燕麦优化算法 AOO
% Number of search agents (population size in optimization algorithms)
PD_no = 30;  
% Dimension of the problem (number of decision variables)
% dim = 20;  
% Objective function: 'cec22_func' (CEC 2022 benchmark function)
% fhd = str2func('cec22_func');  
Max_iter = 1000;
CEC_f = 17;
AOO_fhd = get_CEC_func_str(CEC_f);
[Function_name,F_num] = get_CEC_name(CEC_f);
%% optimization algorithm
for i = 1:1
    if CEC_f == 17 && i == 2 
        continue; % This function (F2) has been deleted
    end
    f_name = get_F_name(i);  %获得函数的序号
    [LB,UB,Dim,F_obj] = Function_name(f_name); %获得函数的边界
    % [Best_score,Best_pos,cg_curve,search_history,ave_fit,x_1st]=GWO2(SearchAgents_no,Max_iter,LB,UB,dim,fobj);
    % Best_score
    [Best_score,Best_pos,cg_curve,search_history,ave_fit,x_1st] = AOOv4(AOO_fhd,Dim,PD_no,Max_iter,LB,UB,i);
    
    %% plot
    figure(i);
    set(gcf, 'Position', [150, 350, 1500, 300]);
    
    
    subplot_width = 1500 / 5;  
    subplot_height = 300 / 1;   
    padding = 0.05;  
    effective_width = (1 - padding * 4) / 5.5;
    effective_height = 0.7 ;
    
    
    left_positions = 0.04 + (0:4) * (effective_width + padding);
    
    
    adjustedColormap = parula;
    
    
    colormap(adjustedColormap);
    subplot('Position', [left_positions(1), padding*3, effective_width, effective_height]);

    % 这是cec2022的图像
    % func_plot(i); 
    % 
    func_plot_cec2017(f_name);
    title(strcat('F', num2str(i)));
    xlabel('x_1');
    ylabel('x_2');
    zlabel('fitness');
    set(gca, 'Box', 'on', 'BoxStyle', 'full', 'LineWidth', 0.25);
    view(3);
    grid on;
    
    % search history
    subplot('Position', [left_positions(2), padding*3, effective_width, effective_height]);
    grid off;
    func_plot2_AOO_cec2017(i);
    hold on;
    scatter(search_history(:,1), search_history(:,2), '.');
    hold on;
    scatter(Best_pos(:,1), Best_pos(:,2), 'r.');
    title('Search history');
    xlabel('x_1');
    ylabel('x_2');
    xlim([LB, UB]);
    ylim([LB, UB]);
    set(gca, 'color', 'none');
    grid on;
    
    % average fitness
    subplot('Position', [left_positions(3), padding*3, effective_width, effective_height]);
    semilogy(ave_fit, 'Linewidth', 1);
    grid on;
    title('Average fitness');
    xlabel('Function calls');
    ylabel('Fitness');
    xlim([0, Max_iter ]);
    set(gca, 'XTickLabel', arrayfun(@(x) num2str(x), 0:6000:30000, 'UniformOutput', false), 'FontAngle', 'normal');
    
    
    % Trajectory of 1st dimension
    subplot('Position', [left_positions(4), padding*3, effective_width, effective_height]);
    plot(x_1st(:,1), 'Linewidth', 1);
    title('Trajectory of 1st dimension');
    grid on;
    xlabel('Function calls');
    ylabel('value');
    xlim([0, Max_iter]);
    set(gca, 'XTickLabel', arrayfun(@(x) num2str(x), 0:6000:30000, 'UniformOutput', false), 'FontAngle', 'normal');
    
    % 收敛曲线图
    subplot('Position', [left_positions(5), padding*3, effective_width, effective_height]);
    semilogy(cg_curve, 'Linewidth', 1);
    title('Convergence curve');
    xlabel('Function calls');
    ylabel('Best score');
    axis tight;
    grid on;
    box on;
    legend('AOO', 'fontsize', 6, 'location', 'best');
    set(gca, 'color', 'none');
    xlim([0, Max_iter]);
    set(gca, 'XTickLabel', arrayfun(@(x) num2str(x), 0:6000:30000, 'UniformOutput', false), 'FontAngle', 'normal');
    % 确保 result 文件夹存在
    result_dir = fullfile(pwd, 'result');  % 构造完整路径
    if ~exist(result_dir, 'dir')
        mkdir(result_dir);  % 创建文件夹
    end

    % 在 saveas 前添加，隐藏坐标区工具栏
    set(gcf, 'Toolbar', 'none');  % gca 表示当前坐标区
    
    % 保存图像到 result 文件夹
    saveas(gcf, fullfile(result_dir, sprintf('AOO_F%d.eps', i)), 'epsc');
    saveas(gcf, fullfile(result_dir, sprintf('AOO_F%d.svg', i)), 'svg');


end


%% 结束日志以及计时
stop_logging();
toc