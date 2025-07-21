% 真菌生长优化算法Fungal growth optimizer
%_________________________________________________________________________%
%  Whale Optimization Algorithm (WOA) source codes demo 1.0               %
%                                                                         %
%  Developed in MATLAB R2019a                                     %
%                                                                         %                                                                        %
%   Main paper: Fungal Growth Optimizer: A Novel Nature-inspired Metaheuristic Algorithm for Stochastic Optimization  %
%                                                                         %
%_________________________________________________________________________%


% The Fungal Growth Optimizer
function [Gb_Fit,Gb_Sol,Conv_curve]=FGO(N,Tmax,ub,lb,dim,Fun_No,fhd)

%%%%-------------------Definitions--------------------------%%
%%
Gb_Sol=zeros(1,dim); % A vector to include the best-so-far solution
Gb_Fit=inf; % A Scalar variable to include the best-so-far score
Conv_curve=zeros(1,Tmax);

%%-------------------Controlling parameters--------------------------%%
%%
M=0.6; %% Determines the tradeoff percent between exploration and exploitation operators.
Ep=0.7; %% Determines the probability of environmental effect on hyphal growth
R=0.9; %% Determines the speed of convergence to the best-so-far solution
%%---------------Initialization----------------------%%
%%
S=initialization(N,dim,ub,lb); % Initialize the S of crested porcupines
t=0; %% Function evaluation counter
%%---------------------Evaluation-----------------------%%
for i=1:N
    %% Test suites of CEC-2014, CEC-2017, CEC-2020, and CEC-2022
    fit(i)=feval(fhd, S(i,:)',Fun_No);
end
% Update the best-so-far solution
[Gb_Fit,index]=min(fit);
Gb_Sol=S(index,:);
%% A new vector to store the best-so-far position for each hyphae
Sp=S;

opt=Fun_No*100; %% Best-known fitness

%%  Optimization Process of FGO
while t<Tmax && opt~=Gb_Fit
    %% ----------------------------------------------------------------------------- %%
    if t <= Tmax/2 %% Compute the the nutrient allocation according to (12)
        nutrients = rand(N); % Allocate more randomly to encourage fluctuation exploitation
    else
        nutrients = fit;  % Exploitation phase: allocate based on fitness
    end
    nutrients = nutrients / sum(nutrients)+2*rand; % Normalize nutrient allocation according to (13)
    %% ----------------------------------------------------------------------------- %%
    if rand<rand %% Hyphal tip growth behavior %%
        for i=1:N
            a=randi(N);
            b=randi(N);
            c=randi(N);
            while a==i | a==b | c==b | c==a ||c==i |b==i
                a=randi(N);
                b=randi(N);
                c=randi(N);
            end
            p=(fit(i)-min(fit))/(max(fit)-min(fit)+eps); %%% Compute p_i according to (23)
            Er=M+(1-t/(Tmax)).*(1-M); %%% Compute Er according to (24)
            if p<Er
                F=(fit(i)/(sum(fit))).*rand*(1-t/(Tmax))^(1-t/(Tmax)); %% Calculate F according to (5) and (6)
                E=exp(F);  %% Calculate E according to (4)
                r1=rand(1,dim);
                r2=rand;
                U1=r1<r2; %% Binary vector
                S(i,:) = (U1).*S(i,:)+(1-U1).*(S(i,:)+E.*(S(a,:)-S(b,:))); %% Generate the new hyphal growth according to (9)
            else
                Ec = (rand(1, dim) - 0.5) .* rand.*(S(a,:)-S(b,:)); % Compute the additional exploratory step using (17)
                if rand<rand %% Hypha growing in the opposite direction of nutrient-rich areas %
                    De2 = rand(1,dim).* (S(i,:) - Gb_Sol).*(rand(1,dim)>rand); %% Compute De2 according to (16) %%
                    S(i,:) = S(i,:) + De2 .* nutrients(i)+Ec*(rand>rand); %% Compute the new growth for the ith hyphal using (15)
                else %% Growth direction toward nutrient-rich area
                    De = rand.* (S(a, :) - S(i, :)) + rand(1,dim).* ((rand>rand*2-1).*Gb_Sol - S(i, :)).*(rand()>R); %% Compute De according to (11) %%
                    S(i,:) = S(i,:) + De .* nutrients(i)+Ec*(rand>Ep); %% Compute the new growth for the ith hyphal using (14)
                end
            end
            %% Return the search agents that exceed the search space's bounds
            for j=1:size(S,2)
                if  S(i,j)>ub(j)
                    S(i,j)=lb(j)+rand*(ub(j)-lb(j));
                elseif  S(i,j)<lb(j)
                    S(i,j)=lb(j)+rand*(ub(j)-lb(j));
                end
                
            end
            % Calculate the fitness value of the newly generated solution
            nF=feval(fhd, S(i,:)',Fun_No);
            %% update Global & Local best solution
            if  fit(i)<nF
                S(i,:)=Sp(i,:);    % Update local best solution
            else
                Sp(i,:)=S(i,:);
                fit(i)=nF;
                %% update Global best solution
                if  fit(i)<=Gb_Fit
                    Gb_Sol=S(i,:);    % Update global best solution
                    Gb_Fit=fit(i);
                end
            end
            t=t+1; % Move to the next generation
            if t>Tmax
                break
            end
            Conv_curve(t)=Gb_Fit;
        end
        if  t>Tmax
            break;
        end
    else
        r5=rand;
        for i=1:N
            rr=rand(1,dim);
            a=randi(N);
            b=randi(N);
            c=randi(N);
            while a==i | a==b | c==b | c==a ||c==i |b==i
                a=randi(N);
                b=randi(N);
                c=randi(N);
            end
            if rand<0.5 %% Hyphal branching
                EL=1+exp(fit(i)/(sum(fit)))*(rand>rand); %% Compute the growth rate of the hypha produced via lateral branching using (29)
                Dep1=(S(b,:)-S(c,:)); %% Compute Dep1 using (26)
                Dep2=(S(a,:)-Gb_Sol); %% Compute Dep1 using (27)
                r1=rand(1,dim);
                r2=rand;
                U1=r1<r2; %% Binary vector
                S(i,:)= S(i,:).*U1+(S(i,:)+ r5.*Dep1.*EL+ (1-r5).*Dep2.*EL).*(1-U1); % Generate the new branch for the ith hyphal using (30)
            else %% Spore germination
                sig=(rand>rand*2-1); %% 1 or -1
                F=(fit(i)/(sum(fit))).*rand*(1-t/(Tmax))^(1-t/(Tmax)); %% Calculate F according to (5) and (6)
                E=exp(F);  %% Calculate E according to (4)
                for j=1:size(S,2)
                    mu=sig.*rand.*E;
                    if rand>rand
                        S(i,j)=(((t/(Tmax))*Gb_Sol(j)+(1-t/(Tmax))*S(a,j))+S(b,j))/2.0 + mu * abs((S(c,j)+S(a,j)+S(b,j))/3.0-S(i,j));      % Eq. (31)
                    end
                end
            end
            %% Return the search agents that exceed the search space's bounds
            for j=1:size(S,2)
                if  S(i,j)>ub(j)
                    S(i,j)=lb(j)+rand*(ub(j)-lb(j));
                elseif  S(i,j)<lb(j)
                    S(i,j)=lb(j)+rand*(ub(j)-lb(j));
                end
                
            end
            % Calculate the fitness value of the newly generated solution
            nF=feval(fhd, S(i,:)',Fun_No);
            %% update Global & Local best solution
            if  fit(i)<nF
                S(i,:)=Sp(i,:);    % Update local best solution
            else
                Sp(i,:)=S(i,:);
                fit(i)=nF;
                %% update Global best solution
                if  fit(i)<=Gb_Fit
                    Gb_Sol=S(i,:);    % Update global best solution
                    Gb_Fit=fit(i);
                end
            end
            t=t+1; % Move to the next generation
            if t>Tmax
                break
            end
            Conv_curve(t)=Gb_Fit;
        end %% End for i
    end %% End If
    if  t>Tmax
        break;
    end
end%% End while
end
