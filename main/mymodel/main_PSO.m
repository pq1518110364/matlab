% 主程序：运行PSO并在3D函数图上标记最优解
clear all; %清除所有变量
close all; %清图
clc; %清屏
D = 2;  % 必须设为2才能绘制3D函数图
pop_size = 100;
iter_max = 5000;
Xmin = -100;
Xmax = 100;
runs = 1;
fhd = str2func('cec17_func');

% 创建结果文件夹
if ~exist('./result', 'dir')
    mkdir('./result');
end

for i = 1:10
    func_num = i;
    if func_num == 2  % 跳过F2
        continue;
    end
    
    for j = 1:runs
        fprintf('正在运行函数F%d，第%d次重复\n', func_num, j);
        % 调用PSO获取最优解
        [gbest, gbestval, FES,conv_curve] = PSO_func(fhd, D, pop_size, iter_max, Xmin, Xmax, func_num);
        xbest(j,:) = gbest;
        fbest(i,j) = gbestval;
        
        % 绘制3D函数图并标记最优解
        plot_cec_function_with_optima(func_num, D, gbest, gbestval);

        % 绘制收敛曲线
        % figure;
        % plot(1:iter_max, conv_curve, 'b-', 'LineWidth', 1.5);
        % xlabel('迭代次数');
        % ylabel('全局最优适应度值');
        % title(['F', num2str(func_num), ' 收敛曲线（第', num2str(j), '次）']);
        % grid on;
        % saveas(gcf, ['./result/F', num2str(func_num), '_run', num2str(j), '.png'], 'png');
        % close;
    end
    
    f_mean(i) = mean(fbest(i,:));
end