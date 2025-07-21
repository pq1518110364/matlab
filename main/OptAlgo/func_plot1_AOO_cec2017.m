function func_plot1_AOO_cec2017(func_name, dim)
    % 可视化CEC2017测试函数（2D截面图）
    % func_name: 函数名称（如 'F1', 'F2'）
    % dim: 函数维度（必须是2,10,30,50,100中的一个）
    
    % 设置默认维度
    if nargin < 2, dim = 30; end
    
    % 验证维度是否合法
    valid_dims = [2, 10, 30, 50, 100];
    if ~ismember(dim, valid_dims)
        error('CEC2017函数仅支持以下维度: 2, 10, 30, 50, 100');
    end
    
    % 根据函数设置定义域范围（不同函数可能需要不同范围）
    [x_range, y_range] = get_cec2017_domain(func_name);
    x = x_range(1):x_range(3):x_range(2);  % 起始:步长:结束
    y = y_range(1):y_range(3):y_range(2);
    
    % 初始化函数值矩阵
    [X, Y] = meshgrid(x, y);
    f = zeros(size(X));
    
    % 计算函数值（固定其他维度为0）
    z = zeros(1, dim-2);  % 创建补零向量（维度为dim-2）
    for i = 1:size(X, 1)
        for j = 1:size(X, 2)
            input = [X(i,j), Y(i,j), z];  % 构造dim维输入向量
            f(i,j) = cec17_func(input', func_name, dim);  % 调用CEC2017函数
        end
    end
    
    % 绘制3D曲面图
    figure;
    if all(f(:) == f(1))  % 检查是否为常量函数
        surf(X, Y, f, 'LineStyle', 'none');  % 常量函数只画曲面
        title(sprintf('CEC2017 %s (常量函数, D=%d)', func_name, dim));
    else
        surfc(X, Y, f, 'LineStyle', 'none');  % 非常量函数画曲面+等高线
        title(sprintf('CEC2017 %s (D=%d)', func_name, dim));
    end
    xlabel('x_1'); ylabel('x_2'); zlabel('f(x)');
    shading interp;
    colorbar;
end

% 辅助函数：获取CEC2017不同函数的定义域范围
function [x_range, y_range] = get_cec2017_domain(func_name)
    % 默认范围
    x_range = [-100, 100, 2];  % 起始,结束,步长
    y_range = [-100, 100, 2];
    
    % 特殊函数范围调整
    switch lower(func_name)
        case 'f8'  % Schwefel函数需要更大范围
            x_range = [-500, 500, 10];
            y_range = [-500, 500, 10];
        case 'f7'  % Griewank函数可以用更小范围
            x_range = [-600, 600, 12];
            y_range = [-600, 600, 12];
        % 可添加更多特殊函数的定义域...
    end
end