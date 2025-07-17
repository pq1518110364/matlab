function plot_cec_function_with_optima(func_num, D, gbest, gbestval)
    if D ~= 2
        error('仅支持二维问题(D=2)的函数图像绘制');
    end
    
    % 根据函数编号设置合适的绘图范围（参考原代码）
    switch func_num
        case 1, x = -100:5:100; y = x; % Sphere
        case 3, x = -10:0.5:10; y = x; % Schwefel 1.02 (noise)
        case 4, x = -100:5:100; y = x; % Schwefel 2.21
        case 5, x = -10:0.5:10; y = x; % Schwefel 2.22
        % 其他函数范围...
        otherwise, x = -100:5:100; y = x; % 默认范围
    end
    
    % 生成网格点并计算函数值
    [X, Y] = meshgrid(x, y);
    Z = zeros(size(X));
    for i = 1:length(x)
        for j = 1:length(y)
            Z(i,j) = cec17_func([X(i,j); Y(i,j)], func_num);
        end
    end
    
    % 绘制3D曲面图
    figure('Position', [100, 100, 800, 600]); % 设置图形窗口大小
    surf(X, Y, Z, 'EdgeColor', 'none'); % 曲面图（无边框，更清晰）
    hold on;
    
    % 标记PSO找到的最优解位置（用红色星号）
    plot3(gbest(1), gbest(2), gbestval, 'r*', 'MarkerSize', 12, 'LineWidth', 2);
    
    % 添加标注文本，显示最优值
    annotation('textbox', [0.65, 0.85, 0.3, 0.1], ... % 位置：右上角
               'String', {['PSO最优值: ', num2str(gbestval, 6)], ...
                         ['位置: (', num2str(gbest(1), 4), ', ', num2str(gbest(2), 4), ')']}, ...
               'EdgeColor', 'none', 'BackgroundColor', 'white', 'FontSize', 10);
    
    % 设置图形属性
    title(['CEC2017 F', num2str(func_num), ' 函数 + PSO最优解']);
    xlabel('x'); ylabel('y'); zlabel('f(x,y)');
    colorbar; % 显示颜色条
    shading interp; % 平滑着色
    view(30, 30); % 设置视角
    grid on;
    
    % 保存图像
    saveas(gcf, ['./result/F', num2str(func_num), '_3D_with_optima.png'], 'png');
    close;
end