function [Best_Score, Best_X, convergence_curve] = AOOv2(fhd, dim, pop_size, iter_max, lb, ub, F_obj, varargin)

% 初始化核心参数
Best_Score = inf;
convergence_curve = zeros(iter_max, 1);

% 物理意义参数（模拟野燕麦种子特性） 特有参数 (通常这些参数是固定的，但如果研究需要，也可以自适应)
x = 3 * rand(1, pop_size) / dim; % 位置相关参数
m = 0.5 * rand(1, pop_size) / dim; % 种子质量（影响扰动幅度）
L = pop_size / dim * rand(1, pop_size); % 种子长度（影响传播范围）
e = m; % 能量参数（与质量关联）
g = 9.8 / dim; % 重力加速度（缩放至问题维度）

% Levy 飞行参数
levy_beta = 1.5; % Levy 指数

X = initialization(pop_size, dim, ub, lb);
Pop_Fit = zeros(pop_size, 1);

%% Record the initial optimal solution and fitness
for i = 1:pop_size
    Pop_Fit(i) = feval(fhd,X(i,:)',varargin{:});
end
[~, Ind] = sort(Pop_Fit);
Best_Score = Pop_Fit(Ind(1));
Best_X = X(Ind(1),:);


convergence_curve(1) = Best_Score;
persistent has_printed

% step = levy(pop_size, dim, levy_beta);

grad_call_count = 0;       % 梯度更新总调用次数
grad_improve_count = 0;    % 梯度更新后结果改进的次数

strategy_local3 = 0; % 策略3，局部探索改进
strategy_local3_improve = 0; % 策略3，局部探索改进

strategy_global = 0; % 策略1，全局探索改进
strategy_global1_improve = 0; % 策略1，全局探索改进
strategy_global2_improve = 0; % 策略1，全局探索改进


%% 迭代循环
for t = 1:iter_max

    % 更新 AOO 参数
    theta = pi * rand(1, pop_size); % 随机角度参数
    c = exp(-10 * t / iter_max); % 衰减因子（随迭代增加而减小）
    %梯度
    grad_prob = 0.1 + 0.9 * (t / iter_max)^2; % 梯度更新概率（随迭代增加而增大）???

    delta = sin(pi/2 * (1 - (2*t / iter_max)).^5);

    % 梯度学习率 (用于可能的梯度引导)
    grad_learning_rate = 0.5 * (1 - t / iter_max); % 梯度学习率（随迭代减小

    % Levy 缩放因子
    levy_scale_factor =  (0.1 * (1 - t / iter_max)^2) * sign(2 * rand() - 1); % Levy缩放因子（随迭代减小，符号随机）
    % 阈值 (探索与利用的平衡)
    threshold = 0.5 + 0.5 * tanh(6 * (t/iter_max - 0.5)); % 探索/利用阈值（动态平衡）

    % Levy 飞行概率动态调整: 早期高，后期趋近于 0
    levy_prob_initial = 0.5; % 初始概率
    levy_prob_final = 0.01;  % 最终概率 (保持一个很小的值，避免完全停止探索)
    prob_power_factor = 2;  % 调整此参数来控制衰减曲线的形状，越大前期下降越慢
    levy_prob = levy_prob_final + (levy_prob_initial - levy_prob_final) * (1 - t / iter_max)^prob_power_factor;

    if t <= iter_max/4
        if Best_Score > 1e4
            p_best_rate = 0.3;
        else
            p_best_rate = 0.1;
        end
    end

    % 排序种群以选择 pbest
    [~, sorted_index] = sort(Pop_Fit, 'ascend');

    % 根据函数复杂度动态调整p_best_rate（可通过初始收敛速度判断）
    % p_best_rate = 0.3 * exp(-5*t/iter_max);
    % 仅打印一次超参数
    if isempty(has_printed)
        fprintf('\n===== AOOv2 超参数 =====\n');
        fprintf('1. 物理参数初始化系数：\n');
        fprintf('   x的系数 = 3, m的系数 = 0.5, L的系数 = pop_size/dim = %.4f, g的系数 = 9.8/dim = %.4f\n', ...
            pop_size/dim, 9.8/dim);
        fprintf('\n2. Levy飞行参数：\n');
        fprintf('   levy_prob = %.2f, levy_beta = %.2f, levy_scale_factor = %.4f\n', ...
            levy_prob, levy_beta, levy_scale_factor);
        fprintf('\n3. 精英引导参数：\n');
        fprintf('   p_best_rate = %.2f\n',p_best_rate);
        fprintf('\n4. 衰减与阈值参数：\n');
        fprintf('   c的衰减系数 = %.1f, threshold的tanh系数 = %.1f\n', ...
            c, threshold);
        fprintf('=========================\n\n');
        has_printed = true;
    end

    a = (1 - t/iter_max) * rand(1,dim)* sign(2 * rand() - 1); % 随迭代减小的随机扰动参数
    % step = levy(pop_size,dim,1.5);

    % 初始化种群
    for i = 1:pop_size

        % 获取当前的个体
        current_X = X(i, :);

        if rand > threshold % 探索阶段：精英引导 + 差异化 + Levy 飞行

           
            % 1. 结合当前种群和归档种群
            popAll = X;

            % 2. 随机选择两个不同的个体 r1 和 r2
            r1_idx = randi(pop_size);
            r2_idx = randi(size(popAll, 1));

            while r1_idx == i, r1_idx = randi(pop_size); end % 确保r1≠i

            while r2_idx == i || r2_idx == r1_idx, r2_idx = randi(size(popAll, 1)); end

            r1 = X(r1_idx, :);
            r2 = popAll(r2_idx, :);

            % 3. 从最佳的前p_best_rate%个体中随机选择一个pbest
            pNP = max(round(p_best_rate * pop_size), 2); % 至少选择两个
            pbest_idx = sorted_index(randi(pNP)); % 随机选择一个精英索引
            pbest = X(pbest_idx, :); % pbest位置
            % AOO 原始的 W 扰动项

            W = c / pi * (2 * rand(1, dim) - 1) .* ub;

            % 引入 Levy 飞行或精英引导/差分混合

            step_levy_current = levy(1, dim, levy_beta); % 保持 levy_beta，而不是固定的1.5
            strategy_global = strategy_global + 1; % 策略1，全局探索改进

            if rand < levy_prob
                % 策略1: Levy飞行和精英引导的组合 (高探索性)
                % Strategy_1  = current_X + levy_scale_factor * step_levy_current .* Best_X + levy_scale_factor * step_levy_current .* pbest ;
                strategy_1  = current_X +  levy_scale_factor * step_levy_current.* (Best_X - current_X);
                newpos1 = max(min(strategy_1, ub), lb);
                new_fit1  = feval(fhd, newpos1', varargin{:});
                if new_fit1 < Pop_Fit(i)
                    X(i, :) = newpos1;
                    Pop_Fit(i) = new_fit1;
                    strategy_global1_improve = strategy_global1_improve + 1;% 策略1，全局探索改进
                end
            else
                % 策略2: 差分进化和AOO原始扰动的组合 (偏向于利用种群信息)
                % strategy_2  = current_X + ((pbest - current_X) * rand() + (r1 - r2) * rand())/2+ W;

                for j = 1 : dim
                    r1 = 1+rand;r2 = 1+rand; % eq.(16)
                    pho_1 = r1 * X(i, :) + (1-r1) * Best_X + r2 * (X(i, :) - Best_X); % eq.(17)
                    pho_2 = X(i,:) + a.* (pbest - Best_X); % eq.(18)
                    % 随机选择两种更新方式之一
                    pos_n3 = X(i,:) ;
                    if rand/j > rand % 早期概率高，后期概率低（随维度增加而降低）
                        pos_n3(j) = pho_1(j);
                    else
                        pos_n3(j) = pho_2(j);
                    end
                    newpos2 = max(min(pos_n3, ub), lb);
                    % newpos2 = max(min(strategy_2, ub), lb);
                    new_fit2  = feval(fhd, newpos2', varargin{:});
                    if new_fit2 < Pop_Fit(i)
                        X(i, :) = newpos2;
                        Pop_Fit(i) = new_fit2;
                        strategy_global2_improve = strategy_global2_improve + 1;% 策略1，全局探索改进
                    end
                end    
                
            end
 
        else % 利用阶段：模拟燕麦受外力（如风、重力）影响的摆动

            if rand < grad_prob
                % 在迭代后期小概率触发，使用梯度更新进行精细搜索
                % 1. 动态选择参考个体（模仿 PIMO）
                % 2. 计算雅可比矩阵
                Jacb = Get_jacobian(F_obj, Best_X, 1e-6);
                grad_norm = Jacb / (norm(Jacb) + eps);
                X(i, :) =  Best_X -  grad_learning_rate * grad_norm;
                grad_call_count = grad_call_count + 1;  % 累计调用次数
            
                % 边界处理与评估
                newpos1 = max(min(X(i,:), ub), lb);
                new_fit  = feval(fhd, newpos1', varargin{:});
                % 新解（newpos1）的适应度是否优于当前个体的历史最优（Fitness(i)）
                if new_fit < Pop_Fit(i)
                    X(i, :) = newpos1;
                    Pop_Fit(i) = new_fit;
                    grad_improve_count = grad_improve_count + 1;  % 算有效改进
                end

            else
                step = levy(1, dim, levy_beta);
                if rand > 0.5 % 子策略2.1：基于"距离感知"的摆动
                    A = ub - abs(ub * t * sin(2 * pi * rand) / iter_max);
                    R = (m(i) * e(i) + L(i)^2) / dim * unifrnd(-A, A, 1, dim);
                    X(i,:) = Best_X + R + c * step.* Best_X;
                else % 子策略2.2：基于"能量衰减"的摆动 (J项)
                    k = 0.5 + 0.5 * rand;
                    B = ub - abs(ub * t * cos(2 * pi * rand) / iter_max);
                    alpha = 1 / pi * exp((randi([0, t]) / iter_max));
                    J = 2 * k * x(i)^2 * sin(2 * theta(i)) / m(i) / g * (1 - alpha) / dim * unifrnd(-B, B, 1, dim);
                    X(i, :) = Best_X + J + c * step.* Best_X;
                end
                strategy_local3 = strategy_local3 + 1 ; % 策略2，局部探索改进

                newpos3 = max(min(X(i,:), ub), lb);
                new_fit3  = feval(fhd, newpos3', varargin{:});
                if new_fit3 < Pop_Fit(i)
                    X(i, :) = newpos3;
                    Pop_Fit(i) = new_fit3;
                    strategy_local3_improve = strategy_local3_improve + 1 ; % 策略2，局部探索改进
                end
            end
        end
    end

    % 确保种群在边界内
    X = boundaryCheck(X, lb, ub);

    % 更新适应度并更新全局最优 特别耗时
    % for i = 1:pop_size
    %     new_fit = feval(fhd,X',varargin{:});
    %     % 直接更新，然后用min查找最优，更符合DE等算法的流程
    %     Pop_Fit(i) = new_fit(i);
    %     if Pop_Fit(i) < Best_Score
    %         Best_Score = Pop_Fit(i);
    %         Best_X = X(i, :);
    %     end
    % end
    % Pop_Fit = feval(fhd,X',varargin{:});
    % for i=1:pop_size
    %     Pop_Fit(i) = feval(fhd,X(i,:)',varargin{:});
    %     if Pop_Fit(i)<Best_Score
    %         Best_Score=Pop_Fit(i);
    %         Best_X=X(i,:);
    %     end
    % end

    [~, Ind] = sort(Pop_Fit);
    if Pop_Fit(Ind(1)) < Best_Score
        Best_Score = Pop_Fit(Ind(1));
        Best_X = X(Ind(1),:);
    end

    % 更新收敛曲线
    convergence_curve(t) = Best_Score;

end
% calculate_improvement_percent('策略 (全局)', strategy_global1_improve, strategy_global);
% calculate_improvement_percent('策略 (全局)', strategy_global2_improve, strategy_global);
% 
% calculate_improvement_percent('策略3 (梯度)', grad_improve_count, grad_call_count);
% calculate_improvement_percent('策略3 (局部)', strategy_local3_improve, strategy_local3);

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
function [ X ] = boundaryCheck(X, lb, ub)
for i=1:size(X,1)
    FU=X(i,:)>ub;
    FL=X(i,:)<lb;
    X(i,:)=(X(i,:).*(~(FU+FL)))+ub.*FU+lb.*FL;
end
end

function X=initialization(SearchAgents_no,dim,ub,lb)

Boundary_no= size(ub,2); % numnber of boundaries

% If the boundaries of all variables are equal and user enter a signle
% number for both ub and lb
if Boundary_no==1
    X=rand(SearchAgents_no,dim).*(ub-lb)+lb;
end
% If each variable has a different lb and ub
if Boundary_no>1
    for i=1:dim
        ub_i=ub(i);
        lb_i=lb(i);
        X(:,i)=rand(SearchAgents_no,1).*(ub_i-lb_i)+lb_i;
    end
end
end

function improve_percent = calculate_improvement_percent(strategy_name, improve_count, total_calls)
% 计算并打印策略的改进百分比
% 输入参数:
%   strategy_name - 策略名称（字符串）
%   total_calls   - 策略总调用次数
%   improve_count - 策略有效改进次数
% 输出参数:
%   improve_percent - 改进百分比（0-100）

% 避免除以0错误
if total_calls == 0
    improve_percent = 0;
else
    % 计算百分比（保留4位小数用于内部计算）
    improve_percent = (improve_count / total_calls) * 100;
end

% 格式化输出，显示2位小数和百分号
fprintf('%s 改进百分比: %.2f%%\n', strategy_name, improve_percent);
end
