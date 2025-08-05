% Projection-Iterative-Methods-based Optimizer: A novel metaheuristic
% algorithm for continuous optimization problems and feature selection
%  投影迭代优化算法                                                                                                   
%  Developed in MATLAB R2020a                                                                 
%                                                                                                     
%  Author :Dongmei Yu, Yanzhe Ji, Yiqiang Xia                                           
%                                                                                                     
%         e-Mail: Yanzhe4723@163.com                                                                                                                                                                                    
%_______________________________________________________________________________________________
% You can simply define your cost function in a seperate file and load its handle to fobj 
% The initial parameters that you need are:
%__________________________________________
% fobj = @YourCostFunction
% dim = number of your variables
% Max_iter = Maximum number of iterations
% N = number of search agents
% lb=[lb1,lb2,...,lbn] where lbn is the lower bound of variable n
% ub=[ub1,ub2,...,ubn] where ubn is the upper bound of variable n
% If all the variables have equal lower bound you can just
% define lb and ub as two single numbers
% 2025-05-05
%______________________________________________________________________________________________
function [Best_Score, Best_Pos, CG_curve] = PIMO(N, Max_iter, lb, ub, dim, fobj)
    Position = initialization(N, dim, ub, lb);  
    Fitness = zeros(N, 1);                    
    epsilon = 1e-6;
    alphabet = 1 : N;
    %% Record the initial optimal solution and fitness
    for i = 1:N
        Fitness(i) = fobj(Position(i,:));      
    end
    [~, Ind] = sort(Fitness);     
    Best_Score = Fitness(Ind(1));
    Best_Pos = Position(Ind(1),:);
    
    CG_curve = zeros(1, Max_iter);% Store convergence information

    grad_proj_improvements = 0; % 策略一
    random_perturb_improvements = 0; % 策略二
    dim_mix_improvements = 0; % 策略三
    levy_improvements = 0; % 策略四
    %% Main optimization loop
    for it = 1:Max_iter
        delta = sin(pi/2 * (1 - (2*it / Max_iter)).^5); % eq.(12)
        a = (1 - it/Max_iter) * rand(1, dim);
        for i = 1:N

            
            %% 1. 选择参考个体（两种方式）
            rr1 = rand;rr2 = rand;
                 % 方式1：基于适应度的轮盘赌选择（适应度越好，被选中概率越高）
            if rr1 >= rr2
                Pr =  (Fitness.^2)/((norm(Fitness)^2));
                P = randsrc(1, 2, [alphabet; Pr']);
                P1 = P(1);P2 = P(2);
            else
                 % 方式2：直接选择当前最优的两个个体
                P = Ind(1:2);
                P1 = P(1);P2 = P(2);
            end

            %% 2. 计算雅可比矩阵（目标函数的局部梯度信息）
            r = rand;
            J = Get_jacobian(fobj, Position(P1, :), epsilon); % Generalised Jacobi matrix 

            %% 3. 生成新位置（两种方向）
            if 3*rand >= 2*rand
                % eq.(14)  方向1：基于参考点与最优解的梯度差
                grad = (r * (Position(P1, :) - Best_Pos) + (1-r) * (Position(P2, :) - Best_Pos)) / 2; 
                pos_n1 = Position(i,:) - delta * grad; % 梯度下降（或上升，取决于符号）
            else
                % eq.(15)  方向2：结合雅可比矩阵的投影更新
                grad = (r * (Position(P2, :) - Best_Pos) + (1-r) * (Position(P1, :) - Best_Pos)) / 2;
                pos_n1 = Position(i,:) + delta * J * grad'; % 雅可比矩阵调整方向
            end

            %% 4. 边界处理与评估
            newpos1 = max(min(pos_n1, ub), lb); % 确保新位置在边界内
            fitt = fobj(newpos1);
            if fitt < Fitness(i) % 如果新位置更优，则更新
                Fitness(i)= fitt;
                Position(i,:) = newpos1;
                grad_proj_improvements = grad_proj_improvements + 1;
            end
            
            %% 策略二：随机扰动更新（增强多样性）
            % 通过随机选择参考个体和方向，增加种群多样性，避免算法过早收敛到次优解。
            if rand > rand
                P = randperm(N, 2); % eq.(13)% 随机选两个个体
                % 重复梯度投影逻辑，增加种群多样性
                P1 = P(1);P2 = P(2);
                if 3*rand >= 2*rand
                    % eq.(14)
                    grad = (r * (Position(P1, :) - Best_Pos) + (1-r) * (Position(P2, :) - Best_Pos)) / 2;
                    pos_n2 = Position(i,:) - delta * grad; % 梯度指导搜索方向
                else
                    % eq.(15)
                    grad = (r * (Position(P2, :) - Best_Pos) + (1-r) * (Position(P1, :) - Best_Pos)) / 2;
                    pos_n2 = Position(i,:) + delta * J * grad'; % 梯度指导搜索方向
                end
                newpos2 = max(min(pos_n2, ub), lb);
                fitt = fobj(newpos2);
                if fitt < Fitness(i)
                    Fitness(i)= fitt;
                    Position(i,:) = newpos2;
                    random_perturb_improvements = random_perturb_improvements + 1;
                end
            end
            
            %% 策略三：维度级混合更新（精细调整）

            for j = 1 : dim
                r1 = 1+rand;r2 = 1+rand; % eq.(16)
                pho_1 = r1 * Position(i, :) + (1-r1) * Best_Pos + r2 * (Position(i, :) - Best_Pos); % eq.(17)
                pho_2 = Position(i,:) + a.* (newpos1 - Best_Pos); % eq.(18)
                % 随机选择两种更新方式之一
                pos_n3 = Position(i,:);
                if rand/j > rand % 早期概率高，后期概率低（随维度增加而降低）
                    pos_n3(j) = pho_1(j);
                else
                    pos_n3(j) = pho_2(j);
                end
                newpos3 = max(min(pos_n3, ub), lb);
                fitt = fobj(newpos3);
                if fitt < Fitness(i)
                    Fitness(i)= fitt;
                    Position(i,:) = newpos3;
                    dim_mix_improvements = dim_mix_improvements +1;
                end
            end     
        end
        
        %%  7. 基于Levy飞行的全局探索（以特定概率触发）
        if rand <= 1/2*(tanh(9*it/Max_iter-5)+1) % 使早期概率高（约 50%），后期趋近于 0；
            pos_n4 = zeros(1, dim);
            d=rand()*(1-it/Max_iter)^2;
            Step_length=levy(N, dim,1.5); % 生成Levy分布的步长
            Elite=repmat(Best_Pos, N, 1); % 复制最优解N份
            r=rand;
            for i=1:N
                for j=1:dim
                    % 基于最优解的Levy扰动
                    pos_n4(j)=r*Elite(i,j)+(1-r)*Step_length(i,j)*d*...
                        (Elite(i,j)-Position(i,j)*(2*it/Max_iter)); % eq.(21)
                end
                newpos4 = max(min(pos_n4, ub), lb);
                fitt = fobj(newpos4);
                if fitt < Fitness(i)
                    Fitness(i) = fitt;
                    Position(i,:) = newpos4;
                    levy_improvements = levy_improvements+1;
                end
            end
        end
       %% Record convergence curve 更新全局最优解并记录收敛曲线
       [~, Ind] = sort(Fitness);     
        if Fitness(Ind(1)) < Best_Score
            Best_Score = Fitness(Ind(1));
            Best_Pos = Position(Ind(1),:);
        end
        CG_curve(it) = Best_Score;

        % disp('swarm-opti')
    end
    %% 打印策略贡献统计
    
    % fprintf('\n===== PIMO 策略贡献统计 =====\n');
    % fprintf('策略一 (梯度投影) 改进次数: %d\n', grad_proj_improvements);
    % fprintf('策略二 (随机扰动) 改进次数: %d\n', random_perturb_improvements);
    % fprintf('策略三 (维度级混合) 改进次数: %d\n', dim_mix_improvements);
    % fprintf('策略四 (Levy飞行) 改进次数: %d\n', levy_improvements);
    % fprintf('=============================\n');
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
function [z] = levy(n,m,beta)
    % beta is set to 1.5 in this paper
    num = gamma(1+beta)*sin(pi*beta/2);
    den = gamma((1+beta)/2)*beta*2^((beta-1)/2);
    sigma_u = (num/den)^(1/beta); % eq.(20)
    u = random('Normal',0,sigma_u,n,m);  % 正态分布随机数
    v = random('Normal',0,1,n,m); % 标准正态分布随机数
    z =u./(abs(v).^(1/beta)); % Levy分布随机数
end
% Initialization function
function X = initialization(N, Dim, UB, LB)
    B_no = size(UB, 2); % number of boundaries
    
    if B_no == 1
        X = rand(N, Dim) .* (UB - LB) + LB;
    end
    % If each variable has a different lb and ub
    if B_no > 1
        X = zeros(N, Dim);
        for i = 1:Dim
            Ub_i = UB(i);
            Lb_i = LB(i);
            X(:, i) = rand(N, 1) .* (Ub_i - Lb_i) + Lb_i;
        end
    end
end