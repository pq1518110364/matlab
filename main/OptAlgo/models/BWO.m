%%%  白鲸优化(BWO)算法 
function [Best_score,Best_pos,curve] = BWO(N,T,lb,ub,D,y)
%% 问题设定

objective=y;        % 评价函数
VarSize=[1 D];   % 决策变量矩阵大小

%% 初始化

habitat.Position=[];
habitat.Cost=[];

pop=repmat(habitat,N,1);

for i=1:N
    pop(i).Position=unifrnd(lb,ub,VarSize); %使用unifrnd函数生成N个随机的决策变量向量，范围在lb和ub之间
    pop(i).Cost=objective(pop(i).Position);
end

% 根据适应度值排序
[~, sortOrder]=sort([pop.Cost]); % sort(A) 按升序（从小到大）对 A 的元素进行排序
pop=pop(sortOrder);

% 用来保存每代的最优值
BestCost=zeros(T,1);
%% BWO 主循环

for it=1:T

    newpop=pop;  %建立临时种群

    Bf=rand()*(1-it/2*T);  %平衡因子
    Wf=0.1-0.05*it/T;    %鲸落概率

    for i=1:N

        if Bf>0.5
            % 勘探
            p=ceil(unifrnd(1,D));  %随机选择一个维度
            r=ceil(unifrnd(1,N));  %随机选择一个个体
            r1=rand();
            r2=rand();
            for j=1:D
                % Eq.(4)
                if mod(j,2)==0  % j为偶数
                    newpop(i).Position(j)=pop(i).Position(p)+(1+r1)*sin(2*pi*r2)*(pop(r).Position(p)-pop(i).Position(p));
                else  % j为奇数
                    newpop(i).Position(j)=pop(i).Position(p)+(1+r1)*cos(2*pi*r2)*(pop(r).Position(p)-pop(i).Position(p));
                end
            end
        else
            % 开发
            r3=rand();
            r4=rand();
            C1=2*r4*(1-it/T);
            r=ceil(unifrnd(1,N));  %随机选择一个个体
            % Eq. (5)
            newpop(i).Position=r3*pop(1).Position-r4*pop(i).Position+C1*levy(1,D,1.5).*(pop(r).Position-pop(i).Position);
        end


        %保证每个分量在取值范围以内
        newpop(i).Position = max(newpop(i).Position, lb);
        newpop(i).Position = min(newpop(i).Position, ub);

        %重新评估
        newpop(i).Cost=objective(newpop(i).Position);


       %这里就可以看作是利用莱维飞行使其跳出局部，增强了全局的开发程度。
        if Bf<Wf
            C2=2*Wf*N;
            Xstep=(ub-lb)*exp(-C2*it/T); % Eq.(9)
            r5=rand();
            r6=rand();
            r7=rand();
            r=ceil(unifrnd(1,N));  %随机选择一个个体
            newpop(i).Position=r5*pop(i).Position-r6*pop(r).Position+r7*Xstep; % Eq.(8)
        end



        %保证每个分量在取值范围以内
        newpop(i).Position = max(newpop(i).Position, lb);
        newpop(i).Position = min(newpop(i).Position, ub);

        %重新评估
        newpop(i).Cost=objective(newpop(i).Position);


        % 贪婪选择
        if newpop(i).Cost<pop(i).Cost
            pop(i).Position=newpop(i).Position;
            pop(i).Cost=newpop(i).Cost;
        end

    end



    [~, sortOrder]=sort([pop.Cost]);
    pop=pop(sortOrder);

    BestCost(it)=pop(1).Cost;
    Bestsol=pop(1).Position;

    %显示每代找到的最优值
%     format long e;
%     disp(['Iteration' num2str(it)]);
%     fprintf('Best Cost is：%40.30f\n',BestCost(it));

end
%% 输出最终结果
curve = BestCost;   
format long e;
Best_score =min(BestCost);
Best_pos= Bestsol;
end

% Levy飞行函数
function [z] = levy(n,m,beta)
num = gamma(1+beta)*sin(pi*beta/2);
den = gamma((1+beta)/2)*beta*2^((beta-1)/2);
sigma_u = (num/den)^(1/beta);
u = random('Normal',0,sigma_u,n,m);
v = random('Normal',0,1,n,m);
z =u./(abs(v).^(1/beta));
end