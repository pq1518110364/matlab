function [Best_Score, Best_X, convergence_curve, diversity_curve, threshold_curve, c_curve] = AOOv2(fhd, dim, pop_size, iter_max, lb, ub, F_obj, varargin)

% 初始化核心参数

Best_Score = inf;

convergence_curve = zeros(iter_max, 1);

diversity_curve = zeros(iter_max, 1);

threshold_curve = zeros(iter_max, 1);

c_curve = zeros(iter_max, 1);

% 物理意义参数（模拟野燕麦种子特性） 特有参数 (通常这些参数是固定的，但如果研究需要，也可以自适应)

x = 3 * rand(1, pop_size) / dim; % 位置相关参数

m = 0.5 * rand(1, pop_size) / dim; % 种子质量（影响扰动幅度）

L = pop_size / dim * rand(1, pop_size); % 种子长度（影响传播范围）

e = m; % 能量参数（与质量关联）

g = 9.8 / dim; % 重力加速度（缩放至问题维度）

% Levy 飞行参数

levy_prob = 0.15; % Levy 飞行引入概率

levy_beta = 1.5; % Levy 指数

levy_scale_factor = 0.01; % Levy 步长的缩放因子，防止过大跳跃

% 精英引导参数

p_best_rate = 0.1; % 优质个体比例

% 归档种群初始化

% archive_size = round(pop_size * 0.5); % 如果需要归档，取消注释

% archive.pop = zeros(archive_size, dim);

% archive.count = 0; % 用于追踪归档中的个体数量

% 初始化种群

X = initialization(pop_size, dim, ub, lb);

Pop_Fit = zeros(pop_size, 1);

% 首次评估种群

for i = 1:pop_size

    Pop_Fit(i) = feval(fhd,X(i,:)',varargin{:});

    if Pop_Fit(i) < Best_Score

        Best_Score = Pop_Fit(i);

        Best_X = X(i, :);

    end

end

convergence_curve(1) = Best_Score;
P = levy(pop_size,dim,1.5);

%% 迭代循环

for t = 1:iter_max

    % 更新 AOO 参数

    theta = pi * rand(1, pop_size);

    c = exp(-10 * t / iter_max); % 衰减因子


    %阙阀值

    threshold = 0.5 + 0.5 * tanh(6 * (t/iter_max - 0.5));


    % 存储当前迭代的阈值和 c

    threshold_curve(t) = threshold;

    c_curve(t) = c;

    % 排序种群以选择 pbest

    [~, sorted_index] = sort(Pop_Fit, 'ascend');

    for i = 1:pop_size

        % 获取当前的个体

        current_X = X(i, :);

        

        if rand > threshold % 探索阶段：精英引导 + 差异化 + Levy 飞行


            % 1. 结合当前种群和归档种群 (如果启用归档，这里将 archive.pop 添加到 popAll)

            % 为了简化，如果未实现归档，popAll 仅为当前种群 X

            % if exist('archive', 'var') && archive.count > 0

            % popAll = [X; archive.pop(1:archive.count, :)];

            % else

            popAll = X;

            % end


            % 2. 随机选择两个不同的个体 r1 和 r2

            r1_idx = randi(pop_size);

            r2_idx = randi(size(popAll, 1)); % r2可以来自当前种群或归档


            while r1_idx == i, r1_idx = randi(pop_size); end

            while r2_idx == i || r2_idx == r1_idx, r2_idx = randi(size(popAll, 1)); end


            r1 = X(r1_idx, :);

            r2 = popAll(r2_idx, :);


            % 3. 从最佳的前p_best_rate%个体中随机选择一个pbest

            pNP = max(round(p_best_rate * pop_size), 2); % 至少选择两个

            pbest_idx = sorted_index(randi(pNP));

            pbest = X(pbest_idx, :);

            % AOO 原始的 W 扰动项

            W = c / pi * (2 * rand(1, dim) - 1) .* ub;

            % 引入 Levy 飞行或精英引导/差分混合

            if rand < levy_prob % 以 levy_prob 的概率进行 Levy 飞行

                % 生成 Levy 步长

                step = levy(1, dim, levy_beta);

                % Levy 飞行更新：朝着远离最优解的方向进行长步幅跳跃，结合缩放因子

                X(i, :) = current_X + levy_scale_factor * step .* (current_X - Best_X);

            else % 否则，使用精英引导和差分扰动的混合策略

                % 融合三种信息进行更新：精英引导 + 差异化 + AOO原始扰动

                % 这里的rand()作为权重系数，增加随机性

                X(i, :) = current_X + (pbest - current_X) * rand() + (r1 - r2) * rand() + W;

            end

        else % 利用阶段：模拟燕麦受外力（如风、重力）影响的摆动
            % 加入梯度，精细搜索
            Jacb = Get_jacobian(F_obj, Best_X, 1e-6); % 精英解处的梯度

            grad_norm = Jacb / norm(Jacb); % 梯度归一化（方向向量）

            % 保持原始 AOO 的利用策略
            
            if rand > 0.5 % 子策略2.1：基于"距离感知"的摆动

                A = ub - abs(ub * t * sin(2 * pi * rand) / iter_max); % 随迭代变化的振幅

                % 扰动项R：结合燕麦特性参数（m_aoo,e,L_aoo）和振幅A 

                % R = (m(i)*e(i) + L(i)^2)/dim * (unifrnd(-A,A,1,dim).* (-grad_norm));

                R = (m(i)*e(i) + L(i)^2)/dim .* (-grad_norm);

                X(i,:) = Best_X + R + c * P(i,:) .* Best_X;

            else % 子策略2.2：基于"能量衰减"的摆动

                k = 0.5 + 0.5 * rand;

                B = ub - abs(ub * t * cos(2 * pi * rand) / iter_max);

                alpha= 1 / pi * exp((randi([0, t]) / iter_max));

                % 扰动项J：模拟燕麦的摆动惯性

                % 注意：x_aoo(i)和theta(i)需要正确引用或定义

                J = 2 * k * x(i)^2 * sin(2 * theta(i)) / m(i) / g * (1 - alpha) / dim * unifrnd(-B, B, 1, dim);

                X(i, :) = Best_X + J + c * P(i, :) .* Best_X; % 注意：P在这里没有定义，需要确保P是存在的

            end

        end

    end


    % 确保种群在边界内

    X = boundaryCheck(X, lb, ub);

    % 更新适应度并更新全局最优

    for i = 1:pop_size

        new_fit = feval(fhd,X',varargin{:});

        % if new_fit < Pop_Fit(i) % 贪婪选择

        % Pop_Fit(i) = new_fit;

        % end

        % 直接更新，然后用min查找最优，更符合DE等算法的流程

        Pop_Fit(i) = new_fit(i);

        if Pop_Fit(i) < Best_Score

            Best_Score = Pop_Fit(i);

            Best_X = X(i, :);

        end

    end

    % 收集多样性数据

    diversity_curve(t) = calculate_diversity(X);

    % 更新收敛曲线

    convergence_curve(t) = Best_Score;

    % 如果需要归档功能，在这里更新 archive

    % archive = updateArchive(archive, pop_old(improved_idx, :), old_fit(improved_idx));

end

end

% 辅助函数：计算种群多样性（平均欧氏距离）

function div = calculate_diversity(X)

[pop_size, ~] = size(X);

if pop_size <= 1

    div = 0; % 单个或无个体，多样性为0

    return;

end


total_dist = 0;

count = 0;

for i = 1:pop_size

    for j = i+1:pop_size

        total_dist = total_dist + norm(X(i,:) - X(j,:)); % 欧氏距离

        count = count + 1;

    end

end

div = total_dist / count;

end

% Finite difference method for approximating generalised Jacobi matrix

function J = Get_jacobian(fobj, x, epsilon)

% f: handle to target function, accepts vector x

% x: current point (vector)

% epsilon: small perturbation value

% 雅可比矩阵计算（梯度近似）

n = length(x);

fx = fobj(x);

m = length(fx);

J = zeros(m, n);

% Calculate the finite difference in each direction according to eqs. (9) and (10).

% 有限差分法计算每个维度的偏导数

for i = 1:n

    x_perturbed = x;

    x_perturbed(i) = x_perturbed(i) + epsilon; % 对第i个变量微小扰动

    f_perturbed = fobj(x_perturbed);

    J(:, i) = (f_perturbed - fx) / epsilon; % 有限差分近似梯度 近似偏导数

end

end