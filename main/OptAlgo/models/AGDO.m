% Application of a novel metaheuristic algorithm inspired by adam gradient descent in distributed permutation flow shop scheduling problem and continuous engineering problems
%  Adam梯度下降优化算法                                                                                                   
%  Developed in MATLAB R2020a                                                                 
%                                                                                                     
%  Author :Yiqiang Xia, Yanzhe Ji                                           
%                                                                                                     
%         e-Mail: xiayiqiang0001@163.com                                                                                                                                                                                    
%_______________________________________________________________________________________________
% You can simply define your cost function in a seperate file and load its handle to fobj 
% The initial parameters that you need are:
%__________________________________________
% fobj = @YourCostFunction
% dim = number of your variables
% T = Maximum number of iterations
% N = number of search agents
% lb=[lb1,lb2,...,lbn] where lbn is the lower bound of variable n
% ub=[ub1,ub2,...,ubn] where ubn is the upper bound of variable n
% If all the variables have equal lower bound you can just
% define lb and ub as two single numbers
% 2025-03-01
%______________________________________________________________________________________________
function [fbest, pbest, cg] = AGDO(N,T, lb, ub, dim,fobj)
    VarSize = [1 dim];% Dimension of the problem
    beta1 = 0.9;   % Exponential decay rates of biased first original moment estimates 
    beta2 = 0.999; % Exponential decay rates of biased second original moment estimates
    epsilon = 1e-8;% Tiny constant for avoiding division by zero
    m = zeros(N, dim);% Initial estimate of the deviation of the gradient
    v = zeros(N, dim);% Initial estimate of the deviation of the second-order moments
    lr = 0.001;  % Learning rate
    alpha = cos((1-rand(N,1))*2*pi); % Eq.(13)
    
    fit = zeros(N, 1);% Fitness of each individual
    bestpos = zeros(1, dim);% Best position
    fbest = inf;%Best score (initially infinite)
    po = initialization(N, dim, ub, lb);% Initialize population
    %% Record the initial optimal solution and fitness
    for i = 1:N      
        fit(i) = fobj(po(i, :));
        if fit(i) <= fbest
            bestpos = po(i, :);
            fbest = fit(i);
        end
    end
    cg = zeros(1, T);% Store convergence information
    %% Main optimization loop
    for t = 1:T
        w = rand()*((1/T^2)*t^2-2/T*t+1/2); % Eq.(12)
        a = (1 - t/T) * rand(VarSize); % Eq.(17)
        xi = rand()*((1/T^2)*t^2-2/T*t+1); % Eq.(31)
        for i = 1:N
            for k = 1:dim/2
                % Eq.(11)
                if k == 1
                    npo_line = w .* po(i, :) + alpha(i) .* po(i, :);
                else
                    npo_line = po(i, :) + (sin(2*pi*dim*t) ) .* npo_line;
                end

                A1 = randperm(N);
                A1(A1 == i) = [];
                % Eq.(14)
                a1 = A1(1);
                a2 = A1(2);
                % Eq.(15)
                zeta = ((fit(a1) - fit(i)) / abs(fit(a1) - fit(i)));
                po_mean = mean(po);
                P = (po_mean - npo_line);% Eq.(23) (first-order)
                f = gtdt(bestpos, P, lr, t, beta1, beta2, epsilon, m, v);
                npo_1a = npo_line + zeta *a.*(f - po(a1, :)) - a.*(npo_line - po(a2, :));% Eq.(16)
                npo_1b = po(a1, :) + a .* (f - po(a2, :));% Eq.(18)
                for j = 1:dim
                    if rand/k > rand
                        npo_line(1, j) = npo_1b(1, j);
                    else
                        npo_line(1, j) = npo_1a(1, j);
                    end
                end
                
                npo_line = max(npo_line, lb);
                npo_line = min(npo_line, ub);
                
                newfit = fobj(npo_line);
                if newfit < fit(i)
                    po(i, :) = npo_line;
                    fit(i) = newfit;
                    if fit(i) <= fbest
                        bestpos = po(i, :);
                        fbest = fit(i);
                    end
                end
                a = rand(VarSize);% Eq.(17)
            end
        end
        if rand <= 1/(1+exp(-(18*(t/T)-12)))
            Step_length=levy(N, dim, 1.2);
            Elite=repmat(bestpos, N, 1);
            for i=1:N
                for j=1:dim
                    npo_2(i,j)=Elite(i,j)+Step_length(i,j)*xi*(Elite(i,j)-po(i,j)*(2*t/T)); % Eq.(27)
                end
            end
            po = npo_2;
            po = max(po, lb);
            po = min(po, ub);
        end
       %% Evaluation
        for i=1:N
            newfit(1, i)=fobj(po(i, :));
        end
        [fbestN, sorted_indexes] = sort(newfit);
        po = po(sorted_indexes(1:N), :);
        if fbestN(1)<fbest
            bestpos = po(1,:);
            fbest=fbestN(1);
        end
       %% Record convergence curve      
        cg(t) = fbest;
        pbest = bestpos;
    end
    
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
    
% Enhanced Gradient Descent with Adam and Momentum
function f = gtdt(X, P, lr, t, beta1, beta2, epsilon, m, v)
    C = 1 - rand;
    grad = X - C .* P;
    %Eq.(21)
    m = beta1 * m + (1 - beta1) * grad;% Updating biased first moment estimates
    v = beta2 * v + (1 - beta2) * (grad .^ 2);% Updating biased second original moment estimates
    %Eq.(20) 
    m_hat = m / (1 - beta1 ^ t);% Calculation of bias-corrected first moment estimates  
    v_hat = v / (1 - beta2 ^ t);% Calculate bias-corrected second original moment estimates
    %Eq.(19)
    f = X - lr * m_hat ./ (sqrt(v_hat) + epsilon);
end

function [z] = levy(n,m,beta)
    % beta is set to 1.5 in this paper
    num = gamma(1+beta)*sin(pi*beta/2);
    den = gamma((1+beta)/2)*beta*2^((beta-1)/2);
    sigma_u = (num/den)^(1/beta); % Eq.(30)
    u = random('Normal',0,sigma_u,n,m);
    v = random('Normal',0,1,n,m);
    z =u./(abs(v).^(1/beta)); % Eq.(29)
end