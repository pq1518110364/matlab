%%%  �׾��Ż�(BWO)�㷨 
function [Best_score,Best_pos,curve] = BWO(N,T,lb,ub,D,y)
%% �����趨

objective=y;        % ���ۺ���
VarSize=[1 D];   % ���߱��������С

%% ��ʼ��

habitat.Position=[];
habitat.Cost=[];

pop=repmat(habitat,N,1);

for i=1:N
    pop(i).Position=unifrnd(lb,ub,VarSize); %ʹ��unifrnd��������N������ľ��߱�����������Χ��lb��ub֮��
    pop(i).Cost=objective(pop(i).Position);
end

% ������Ӧ��ֵ����
[~, sortOrder]=sort([pop.Cost]); % sort(A) �����򣨴�С���󣩶� A ��Ԫ�ؽ�������
pop=pop(sortOrder);

% ��������ÿ��������ֵ
BestCost=zeros(T,1);
%% BWO ��ѭ��

for it=1:T

    newpop=pop;  %������ʱ��Ⱥ

    Bf=rand()*(1-it/2*T);  %ƽ������
    Wf=0.1-0.05*it/T;    %�������

    for i=1:N

        if Bf>0.5
            % ��̽
            p=ceil(unifrnd(1,D));  %���ѡ��һ��ά��
            r=ceil(unifrnd(1,N));  %���ѡ��һ������
            r1=rand();
            r2=rand();
            for j=1:D
                % Eq.(4)
                if mod(j,2)==0  % jΪż��
                    newpop(i).Position(j)=pop(i).Position(p)+(1+r1)*sin(2*pi*r2)*(pop(r).Position(p)-pop(i).Position(p));
                else  % jΪ����
                    newpop(i).Position(j)=pop(i).Position(p)+(1+r1)*cos(2*pi*r2)*(pop(r).Position(p)-pop(i).Position(p));
                end
            end
        else
            % ����
            r3=rand();
            r4=rand();
            C1=2*r4*(1-it/T);
            r=ceil(unifrnd(1,N));  %���ѡ��һ������
            % Eq. (5)
            newpop(i).Position=r3*pop(1).Position-r4*pop(i).Position+C1*levy(1,D,1.5).*(pop(r).Position-pop(i).Position);
        end


        %��֤ÿ��������ȡֵ��Χ����
        newpop(i).Position = max(newpop(i).Position, lb);
        newpop(i).Position = min(newpop(i).Position, ub);

        %��������
        newpop(i).Cost=objective(newpop(i).Position);


       %����Ϳ��Կ�����������ά����ʹ�������ֲ�����ǿ��ȫ�ֵĿ����̶ȡ�
        if Bf<Wf
            C2=2*Wf*N;
            Xstep=(ub-lb)*exp(-C2*it/T); % Eq.(9)
            r5=rand();
            r6=rand();
            r7=rand();
            r=ceil(unifrnd(1,N));  %���ѡ��һ������
            newpop(i).Position=r5*pop(i).Position-r6*pop(r).Position+r7*Xstep; % Eq.(8)
        end



        %��֤ÿ��������ȡֵ��Χ����
        newpop(i).Position = max(newpop(i).Position, lb);
        newpop(i).Position = min(newpop(i).Position, ub);

        %��������
        newpop(i).Cost=objective(newpop(i).Position);


        % ̰��ѡ��
        if newpop(i).Cost<pop(i).Cost
            pop(i).Position=newpop(i).Position;
            pop(i).Cost=newpop(i).Cost;
        end

    end



    [~, sortOrder]=sort([pop.Cost]);
    pop=pop(sortOrder);

    BestCost(it)=pop(1).Cost;
    Bestsol=pop(1).Position;

    %��ʾÿ���ҵ�������ֵ
%     format long e;
%     disp(['Iteration' num2str(it)]);
%     fprintf('Best Cost is��%40.30f\n',BestCost(it));

end
%% ������ս��
curve = BestCost;   
format long e;
Best_score =min(BestCost);
Best_pos= Bestsol;
end

% Levy���к���
function [z] = levy(n,m,beta)
num = gamma(1+beta)*sin(pi*beta/2);
den = gamma((1+beta)/2)*beta*2^((beta-1)/2);
sigma_u = (num/den)^(1/beta);
u = random('Normal',0,sigma_u,n,m);
v = random('Normal',0,1,n,m);
z =u./(abs(v).^(1/beta));
end