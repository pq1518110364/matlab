%__________________________________________________________________
%  Boomerang aerodynamic ellipse optimizer
%  Developed in MATLAB R2022b
% 回旋镖气动椭圆优化算法Boomerang Aerodynamic Ellipse Optimizer (BAEO)
%  programmer: Shijie Zhao and Fanshuai Meng   
%  E-mail: zhaoshijie@lntu.edu.cn
%          mfs18841817727@163.com
%  The code is based on the following papers.
%  Shijie Zhao, Fanshuai Meng, Liang Cai,and Ronghua Yang 
%  Boomerang aerodynamic ellipse optimizer: A human game-inspired optimization technique for 
%  numerical optimization and multilevel thresholding image segmentation.
%__________________________________________________________________

function [Best_fitness,Best_position,Convergence_curve]=BAEO(Popsize,Maxiteration,LB,UB,Dim,Fobj)
tic;
%% 初始化种群
Boomerang=initialization(Popsize,Dim,UB,LB);     
BoomerangFitness = zeros(1,Popsize);             
Convergence_curve=zeros(1,Maxiteration);    
for i=1:Popsize
    BoomerangFitness(1,i)=Fobj(Boomerang(i,:));
end
%% Calculate the fitness values of initial boomerangs.
% 记录最优解
[~,sorted_indexes]=sort(BoomerangFitness);
Best_position=Boomerang(sorted_indexes(1),:);
Best_fitness = BoomerangFitness(sorted_indexes(1));
Convergence_curve(1)=Best_fitness; % 记录每代最优适应度，用于分析收敛速度。
t=1;

max_diff=zeros(Popsize,Dim);
min_diff=max_diff;

while t<Maxiteration+1
    Elite=repmat(Best_position,Popsize,1); % 复制最优解为精英引导
    alpha=((t-1)/Maxiteration-1)^4; % 自适应参数（控制搜索步长）
    beta=0.5;
    
    %% The stages of the boomerang movement
    Boomerang_1=Boomerang+0.3.*alpha.*(max_diff+min_diff).*(2.*rand(Popsize,Dim)-ones(Popsize,Dim))+0.5.*(Elite-Boomerang);

    max_diff=max(max_diff,(Boomerang_1-Boomerang));
    min_diff=min(min_diff,(Boomerang_1-Boomerang));

    Boomerang=Boomerang_1;

    % Out-of-bounds processing
    Tp=Boomerang>UB';
    Tm=Boomerang<LB';
    Boomerang=(Boomerang.*(~(Tp+Tm)))+UB'.*Tp+LB'.*Tm;
    BoomerangFitness=Fobj(Boomerang);

    [~,sorted_indexes]=sort(BoomerangFitness);

    Boomerang=Boomerang(sorted_indexes(1:Popsize),:);
    %% Uniform local mining stage of boomerang 
    % The number of generated positions can be adjusted according to different fitness value functions
    n=10;                % Number of positions 仅对前20%的优质解进行局部优化
    for i=1:round(Popsize/5)
        gv=abs(Boomerang(i,:)-Elite(i,:)).*rand(1,Dim);

        if sum(gv==0)>=1
            continue;
        end

        % Generation position and acceptance probability   生成正态分布的候选解
        ner_Boomerang=repmat(Boomerang(i,:),n,1);
        ads=random('Normal', ner_Boomerang, 1 ,n,Dim);
        round_point=ads;
        sum_point=0;
        qop=repmat(Boomerang(i,:),n,1);
        sum_point=sum(((round_point-qop).^2),2);
        sum_point=sqrt(sum_point);
        gv=abs(Boomerang(i,:)-Elite(i,:)).*rand(1,Dim);
        round_point=gv.*round_point./sum_point;
        p_point=0;
        sp=1;
        lv=cumprod(gv);
        sp=lv(Dim)^4;
        p_point=sum(round_point(:,:).^2*sp./gv.^4,2);
        lv=cumprod(gv);
        sp=lv(Dim)^4;
        if sp==0
            continue;
        end
        sp=sp/(min(gv)^2);
        p_point=p_point./sp;
        p_point=sqrt(p_point);
        
        % Select the generated location according to probability
        idx=(rand(n,1)<=p_point);
        rd_point=round_point(idx,:);
        point_num=size(rd_point);

        % Out-of-bounds processing
        Tp=rd_point>UB';
        Tm=rd_point<LB';
        rd_point=(rd_point.*(~(Tp+Tm)))+UB'.*Tp+LB'.*Tm;
        if size(rd_point,1)==0
            continue;
        end

        new_pointfit=Fobj(rd_point);

        [~,sorted_indexes1o]=sort(new_pointfit);

        new_pointfit=new_pointfit(sorted_indexes1o(1:point_num(1)));

        if point_num(1)>0 && BoomerangFitness(1,i)>new_pointfit(1)
            Boomerang(i,:)=rd_point(1,:);
            BoomerangFitness(1,i)=new_pointfit(1);
        end
    end

    %% Update of optimal solution and optimal fitness value
    
    % Out-of-bounds processing
    Tp=Boomerang>UB';
    Tm=Boomerang<LB';
    Boomerang=(Boomerang.*(~(Tp+Tm)))+UB'.*Tp+LB'.*Tm;
    BoomerangFitness=Fobj(Boomerang);

    [~,sorted_indexes]=sort(BoomerangFitness);

    Boomerang=Boomerang(sorted_indexes(1:Popsize),:);
    SortfitbestN = BoomerangFitness(sorted_indexes(1:Popsize));

    %Update the optimal solution
    if SortfitbestN(1)<Best_fitness
        Best_position=Boomerang(1,:);
        Best_fitness=SortfitbestN(1);
    end
    Convergence_curve(t)=Best_fitness;
    t = t + 1;
end
end

