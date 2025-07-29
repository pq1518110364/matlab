%%Animated Oat Optimization Algorithm 不实野燕麦优化算法
function [Best_Score,Best_X,convergence_curve]= AOOv1(fhd,dim,pop_size,iter_max,lb,ub,F_obj,varargin)
%% initialize
Best_Score=inf; 
convergence_curve=zeros(iter_max,1);

x = 3 * rand(1, pop_size) / dim;
m = 0.5 * rand(1, pop_size) / dim;
L = pop_size / dim * rand(1, pop_size);
e = m;
g = 9.8/dim;
X=initialization(pop_size,dim,ub,lb);


for i = 1:pop_size
    Pop_Fit(i) = feval(fhd,X(i,:)',varargin{:});
    if Pop_Fit(i)<Best_Score
        Best_Score=Pop_Fit(i);
        Best_X=X(i,:);
    end
end
convergence_curve(1) = Best_Score;
epsilon = 1e-6;

%% iteration part
% delta = sin(pi/2 * (1 - (2*it / Max_iter)).^5); % eq.(12)
for t = 1 :1:iter_max
    theta = pi * rand(1,pop_size); % 随机角度（模拟燕麦摆动角度）是否可以改成梯度方向，更有利于在局部的探索
    % c = (1 -  t / iter_max)^3; % 随迭代衰减的系数（平衡探索与利用，前期大后期小）是否可以增大系统，注意到有时候收敛不太迅速的问题
    c = exp(-10 * t / iter_max);  % 新的 c 衰减策略
    P = levy(pop_size,dim,1.5);
    % 动态阈值（随迭代增大，后期偏向利用）
    threshold = 0.5 + 0.5 * tanh(6 * (t/iter_max - 0.5));
    % adaptive_scale  = (1 - t/iter_max); % 计算自适应缩放因子 (标量)
    % 逐个更新种群个体
    for i = 1 :pop_size
       %% 随机策略 探索与利用各占一半是否合理
        if rand > threshold
            W = c/pi *(2 * rand(1,dim) - 1) .* ub  ; % 策略1：模拟燕麦的"基础生长"或"群体互动"
            %% 根据个体索引选择不同更新方式（模拟不同生长阶段） 是否可以采用精英引导
            if mod(i,pop_size/10) == 0 % 每10%个体围绕种群均值更新（局部探索）
                X(i,:) = mean(X) + W;
            elseif mod(i,pop_size/10) == 1 % 每10%个体围绕全局最优更新（局部利用）
                X(i,:) = Best_X + W ;
            else
                X(i,:) = X(i,:) + W ;  % 其余个体围绕自身更新（保持多样性）
            end

        else
             % 新增一个子策略，基于梯度信息
            % if t > 0.5 * iter_max && rand < 0.2 % 以小概率进行梯度下降
            %     % alpha 是学习率，可随迭代衰减
            %     % alpha = (1 - t/iter_max)^2 * 0.1; 
            %     jacb = Get_jacobian(F_obj, X(i, :), epsilon); 
            %     delta = sin(pi/2 * (1 - (2*t / iter_max)).^5);
            %     X(i,:) = X(i,:) - delta * jacb;
            % else
                 %% 策略2：模拟燕麦受外力（如风、重力）影响的摆动
                if rand > 0.5 % 子策略2.1：基于"距离感知"的摆动
                    A = ub - abs(ub * t * sin(2 * pi * rand) / iter_max); % 随迭代变化的振幅
                    % 扰动项R：结合燕麦特性参数（m,e,L）和振幅A 是否有更好的扰动策略
                    R = (m(i) * e(i) + L(i) ^2) /dim * unifrnd( -A , A, 1, dim);
                    X(i,:) = Best_X + R  + c * P(i,:) .* Best_X  ;

                    %% 子策略2.2：基于"能量衰减"的摆动
                else 
                    k = 0.5 + 0.5 * rand;
                    B = ub - abs(ub * t * cos(2 * pi * rand) / iter_max);
                    alpha = 1 / pi * exp((randi([0,t]) / iter_max));

                    % 扰动项J：模拟燕麦的摆动惯性（结合质量m、重力g、角度theta等）
                    J = 2 * k * x(i)^2 * sin (2 * theta(i)) / m(i) / g * (1 - alpha) /dim  * unifrnd( -B , B, 1, dim);
                    X(i,:) = Best_X + J + c * P(i,:) .* Best_X ;
                end
            % end 

            
        end
    end

    X = boundaryCheck(X, lb, ub);
    Pop_Fit = feval(fhd,X',varargin{:});
    for i=1:pop_size
        Pop_Fit(i) = feval(fhd,X(i,:)',varargin{:});
        if Pop_Fit(i)<Best_Score
            Best_Score=Pop_Fit(i);
            Best_X=X(i,:);
        end
    end
    convergence_curve(t)=Best_Score;
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
