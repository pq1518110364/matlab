% 鲸鱼迁徙优化算法Whale Migration Algorithm，WMA 大规模多仓库多旅行商问题原理
%  Whale Migration Algorithm (WMA) in MATLAB
% please follow weChat Swarm-Opti
function [Destination_fitness,Destination_position,Convergence_curve]=WMA(N,Max_iteration,lb,ub,dim,fobj)


CostFunction = @(x) fobj(x);

%% Problem Definition
nVar = dim;          % Number of Variables
VarSize = [1 nVar];
VarMin = lb;         %Variables Lower Bound
VarMax =ub ;         %Variables Upper Bound

%% Whale Migration Algorithm Parameters
MaxIt = Max_iteration;        % Maximum Number of Iterations
nPop = N;                     % Population Size
NL=round(nPop/2);             %Number of Leader Whales


%% Initialization

% Empty Whale Structure
Whale.Position=[];
Whale.Cost=[];

% Initialize Population Array
pop=repmat(Whale,nPop,1);

% Initialize Best Solution Ever Found
BestSol.Cost=inf;

% Create Initial Whales
for i=1:nPop
    pop(i).Position=unifrnd(VarMin,VarMax,VarSize);
    pop(i).Cost=CostFunction(pop(i).Position);

    if pop(i).Cost<=BestSol.Cost
        BestSol=pop(i);
    end

end
it=0;
FEs=nPop;
while FEs<MaxIt
    it=it+1;
    % Calculate Leader Whales Mean
    [~, SortOrder]=sort([pop.Cost]);
    pop=pop(SortOrder);
    Mean = 0;
    for i=1:NL
        Mean = Mean + pop(i).Position;
    end
    Mean = Mean/NL;
    for i=1:nPop
        if i>NL
            newsol.Position=Mean+(rand(VarSize)).*(pop(i-1).Position-pop(i).Position)+(rand(VarSize)).*(BestSol.Position-Mean);%%%m

            newsol.Position=max(newsol.Position,VarMin);
            newsol.Position=min(newsol.Position,VarMax);
            newsol.Cost=CostFunction(newsol.Position);
            FEs=FEs+1;
            if newsol.Cost<=pop(i).Cost
                pop(i)=newsol;
            end
            if newsol.Cost<=BestSol.Cost
                BestSol=newsol;
            end
        end

        %%%%%%%%%% Movement the Leader Whales

        if   i<=NL
            newsol.Position=pop(i).Position+(rand).*unifrnd(1*VarMin,1*VarMax,VarSize);

            newsol.Position=max(newsol.Position,VarMin);
            newsol.Position=min(newsol.Position,VarMax);
            newsol.Cost=CostFunction(newsol.Position);
            FEs=FEs+1;
            if newsol.Cost<=pop(i).Cost
                pop(i)=newsol;
            end
            if newsol.Cost<=BestSol.Cost
                BestSol=newsol;
            end
        end
    end

    [~, SortOrder]=sort([pop.Cost]);
    pop=pop(SortOrder);
    BestCost(it)=BestSol.Cost;
    Convergence_curve(it)=BestSol.Cost;

    % Show Iteration Information
    if mod(it,100)==0
        disp(['NFEs ' num2str(FEs) ': Best Cost = ' num2str(BestCost(it))]);
    end


end
Destination_fitness=BestSol.Cost;
Destination_position=BestSol.Position;


end

