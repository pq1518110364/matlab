%%Animated Oat Optimization Algorithm 
function [Best_Score,Best_X,convergence_curve,search_history,ave_fit,x_1st]= AOOv4(fhd,dim,pop_size,iter_max,lb,ub,varargin)
%% initialize
Best_Score=inf; 
% 收敛曲线
convergence_curve=zeros(iter_max,1); 
% 搜索历史 存储所有代的种群位置
search_history = zeros(iter_max*pop_size, dim); % search history
% 平均适应度
ave_fit = zeros(iter_max, 1); % average fitness
x_1st= zeros(iter_max, dim); % 1st trajectory 第一个个体的轨迹：存储每代第一个个体的位置

% 燕麦特性参数初始化（模拟燕麦的生长/运动属性）
x = 3 * rand(1, pop_size) / dim; % 可能对应燕麦的生长速度/摆动幅度
m = 0.5 * rand(1, pop_size) / dim; % 可能对应燕麦的质量（影响惯性）
L = pop_size / dim * rand(1, pop_size);  % 可能对应燕麦的长度
e = m; % 可能对应弹性系数（与质量相关）
g = 9.8/dim;  % 缩放后的重力加速度（影响摆动）
X=initialization(pop_size,dim,ub,lb);


for i = 1:pop_size
    Pop_Fit(i) = feval(fhd,X(i,:)',varargin{:}); % 计算第i个个体的适应度
    if Pop_Fit(i)<Best_Score % 若当前个体更优，更新全局最优
        Best_Score=Pop_Fit(i);
        Best_X=X(i,:); 
    end
end
convergence_curve(1) = Best_Score; % 记录第1代的最优适应度

%% iteration part

for t = 1 :1:iter_max
    theta = pi * rand(1,pop_size); % 随机角度（模拟燕麦摆动角度）
    c = (1 -  t / iter_max)^3; % 随迭代衰减的系数（平衡探索与利用，前期大后期小）
    P = levy(pop_size,dim,1.5); % 生成Levy飞行随机数（增强全局探索能力）

     % 逐个更新种群个体
    for i = 1 :pop_size
        if rand > 0.5  % 策略1：模拟燕麦的"基础生长"或"群体互动"
            W = c/pi *(2 * rand(1,dim) - 1) .* ub;  % 随机扰动项（受c和边界缩放）
            % 根据个体索引选择不同更新方式（模拟不同生长阶段）
            if mod(i,pop_size/10) == 0  % 每10%个体围绕种群均值更新（局部探索）
                X(i,:) = mean(X) + W;
            elseif mod(i,pop_size/10) == 1  % 每10%个体围绕全局最优更新（局部利用）
                X(i,:) = Best_X + W ;
            else  % 其余个体围绕自身更新（保持多样性）
                X(i,:) = X(i,:) + W ;
            end
        else  % 策略2：模拟燕麦受外力（如风、重力）影响的摆动
            if rand > 0.5  % 子策略2.1：基于"距离感知"的摆动
                A = ub - abs(ub * t * sin(2 * pi * rand) / iter_max);  % 随迭代变化的振幅
                % 扰动项R：结合燕麦特性参数（m,e,L）和振幅A
                R = (m(i) * e(i) + L(i) ^2) /dim * unifrnd( -A , A, 1, dim);
                X(i,:) = Best_X + R + c * P(i,:) .* Best_X;  % 结合全局最优和Levy飞行
            else  % 子策略2.2：基于"能量衰减"的摆动
                k = 0.5 + 0.5 * rand;  % 随机系数（模拟风力强度）
                B = ub - abs(ub * t * cos(2 * pi * rand) / iter_max);  % 随迭代变化的振幅（余弦变化）
                alpha = 1 / pi * exp((randi([0,t]) / iter_max));  % 能量衰减系数（随迭代增加）
                % 扰动项J：模拟燕麦的摆动惯性（结合质量m、重力g、角度theta等）
                J = 2 * k * x(i)^2 * sin(2 * theta(i)) / m(i) / g * (1 - alpha) /dim * unifrnd( -B , B, 1, dim);
                X(i,:) = Best_X + J + c * P(i,:) .* Best_X;  % 结合全局最优和Levy飞行
            end
        end
    end

    % 边界检查：将超出[lb,ub]的解拉回边界
    X = boundaryCheck(X, lb, ub);
    % 批量计算适应度（向量形式
    Pop_Fit = feval(fhd,X',varargin{:});
    for i=1:pop_size
        % 若需更精确，可单独计算（此处与批量计算结果一致，可能冗余）
        Pop_Fit(i) = feval(fhd,X(i,:)',varargin{:});
        if Pop_Fit(i)<Best_Score
            Best_Score=Pop_Fit(i);
            Best_X=X(i,:);
        end
    end
     % 记录迭代信息
    search_history ((t-1)*pop_size+1:t*pop_size,:)=X; % search history
    ave_fit(t) = mean(Pop_Fit); % average fitness
    x_1st(t,:)= X(1,:); % 1st trajectory
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

% If the boundaries of all variables are equal and user enter a signle 若各维度边界不同
% number for both ub and lb
if Boundary_no==1
    X=rand(SearchAgents_no,dim).*(ub-lb)+lb;
end
% If each variable has a different lb and ub
if Boundary_no>1
    for i=1:dim
        ub_i=ub(i);
        lb_i=lb(i);
        X(:,i)=rand(SearchAgents_no,1).*(ub_i-lb_i)+lb_i;  % 按维度生成随机解
    end
end
end
