% fhd	目标函数句柄（如 CEC 测试函数），需要最小化的函数
% Dimension	问题维度（决策变量的数量）
% Particle_Number	粒子数量（种群规模）
% Max_Gen	最大迭代次数（优化终止条件）
% VRmin, VRmax	决策变量的下界和上界（可以是标量或与维度长度相同的向量）
% varargin	可选参数，传递给目标函数的额外参数（如 CEC 函数的编号）
function [gbest,gbestval,fitcount,convergence_curve]= PSO_func(fhd,Dimension,Particle_Number,Max_Gen,VRmin,VRmax,varargin)
% gbest	全局最优解（找到的最优决策变量向量）
% gbestval	全局最优适应度值（gbest对应的目标函数值，越小越好）
% fitcount	累计的函数评估次数（优化过程中调用目标函数的总次数）
%[gbest,gbestval,fitcount]= PSO_func('f8',3500,200000,30,30,-5.12,5.12)


rand('state',sum(100*clock));
me=Max_Gen;
ps=Particle_Number;
D=Dimension;
cc=[2 2];   %acceleration constants
iwt=0.9-(1:me).*(0.5./me);
% iwt=0.5.*ones(1,me);
if length(VRmin)==1
    VRmin=repmat(VRmin,1,D);
    VRmax=repmat(VRmax,1,D);
end
mv=0.5*(VRmax-VRmin);
VRmin=repmat(VRmin,ps,1);
VRmax=repmat(VRmax,ps,1);
Vmin=repmat(-mv,ps,1);
Vmax=-Vmin;
pos=VRmin+(VRmax-VRmin).*rand(ps,D);

e=feval(fhd,pos',varargin{:});

fitcount=ps;
vel=Vmin+2.*Vmax.*rand(ps,D);%initialize the velocity of the particles
pbest=pos;
pbestval=e; %initialize the pbest and the pbest's fitness value
[gbestval,gbestid]=min(pbestval);
gbest=pbest(gbestid,:);%initialize the gbest and the gbest's fitness value
gbestrep=repmat(gbest,ps,1);

convergence_curve = zeros(1, Max_Gen);  % 新增：存储每代的全局最优值
convergence_curve(1) = gbestval;  % 第1代的最优值

for i=2:me

        aa=cc(1).*rand(ps,D).*(pbest-pos)+cc(2).*rand(ps,D).*(gbestrep-pos);
        vel=iwt(i).*vel+aa;
        vel=(vel>Vmax).*Vmax+(vel<=Vmax).*vel;
        vel=(vel<Vmin).*Vmin+(vel>=Vmin).*vel;
        pos=pos+vel;
        pos=((pos>=VRmin)&(pos<=VRmax)).*pos...
            +(pos<VRmin).*(VRmin+0.25.*(VRmax-VRmin).*rand(ps,D))+(pos>VRmax).*(VRmax-0.25.*(VRmax-VRmin).*rand(ps,D));
        e=feval(fhd,pos',varargin{:});
        fitcount=fitcount+ps;
        tmp=(pbestval<e);
        temp=repmat(tmp',1,D);
        pbest=temp.*pbest+(1-temp).*pos;
        pbestval=tmp.*pbestval+(1-tmp).*e;%update the pbest
        [gbestval,tmp]=min(pbestval);
        gbest=pbest(tmp,:);
        gbestrep=repmat(gbest,ps,1);%update the gbest
        % 更新全局最优后，记录当前迭代的最优值
        convergence_curve(i) = gbestval;
    end
end

