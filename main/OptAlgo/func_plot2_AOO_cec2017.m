% This function draw the benchmark functions

function func_plot2_AOO_cec2017(func_name)

% 绘制CEC 2017测试函数的等高线图
    % D: 合法维度（2,10,30,50,100）
    if nargin < 2, D = 30; end  % 默认使用30维
    
    % 设置定义域（根据函数调整）
    x = -100:5:100; y = -100:5:100;
    
    % 计算函数值
    [X, Y] = meshgrid(x, y);
    f = zeros(size(X));
    z = zeros(1, D-2);  % 关键修改：补零到合法维度
    
    for i = 1:size(X, 1)
        for j = 1:size(X, 2)
            input = [X(i,j), Y(i,j), z];  % 构造D维输入
            f(i,j) = cec17_func(input', func_name);
        end
    end

contour(x,y,f)
end

