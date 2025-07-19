clear all 
close all
clc
addpath(genpath(pwd)); %将当前工作目录以及其子目录添加到 MATLAB 搜索路径中，以便可以找到相关函数和文件。
%% 启动日志和diary
start_logging();

%% 数据读取与处理，测试cec函数时这一操作没必要

%% 设置参数
PD_no = 50; % Number of sand cat 沙猫优化算法
Max_iter = 1000; % 最大迭代次数
CEC_f = 05;
F_no = 10;

% log_message(sprintf('参数设置: PD_no=%d, Max_iter=%d', PD_no, Max_iter), 'INFO');

% [Function_name,F_num] = get_CEC_name(b);
% Best = zeros(F_num,MaxA);  %存储最优适应度值
% Mean = zeros(F_num,MaxA);  %存储平均适应度值
% Std = zeros(F_num,MaxA);   %存储适应度方差
% F_sum = zeros(F_num*3,MaxA); %最优值，平均值和方差过渡
% next_sum = zeros(F_num*3,MaxA); %存储最优值，平均值和方差


% try
%% 调用函数，获得参数
% 运行算法

[Function_name,F_num] = get_CEC_name(CEC_f);
F_name = get_F_name(F_no); %获得函数的序号
[LB,UB,Dim,F_obj] = Function_name(F_name); %获得函数的边界

[Best_pos,Best_score,GWO_cg_curve] = GWO(PD_no,Max_iter,LB,UB,Dim,F_obj); % Call BWO
fprintf ('Best solution obtained by GWO: %s\n', num2str(Best_score,'%e  '));
display(['The best optimal value of the objective funciton found by BWO  for ' [num2str(F_name)],'  is : ', num2str(Best_pos)]);

%% BWO    
[Best_pos,Best_score, BWO_cg_curve ] = BWO(PD_no,Max_iter,LB,UB,Dim,F_obj); % Call BWO
fprintf ('Best solution obtained by BWO: %s\n', num2str(Best_score,'%e  '));
display(['The best optimal value of the objective funciton found by BWO  for ' [num2str(F_name)],'  is : ', num2str(Best_pos)]);

%% SSA    
[Best_pos,Best_score, SSA_cg_curve ] = SSA(PD_no,Max_iter,LB,UB,Dim,F_obj); % Call SSA
fprintf ('Best solution obtained by SSA: %s\n', num2str(Best_score,'%e  '));
display(['The best optimal value of the objective funciton found by SSA  for ' [num2str(F_name)],'  is : ', num2str(Best_pos)]);

%% SCSO
[BsSCSO,BpSCSO,SCSO_cg_curve]=SCSO(PD_no,Max_iter,LB,UB,Dim,F_obj); % Call MSCSO
fprintf ('Best solution obtained by SCSO: %s\n', num2str(BsSCSO,'%e  '));
display(['The best optimal value of the objective funciton found by SCSO  for ' [num2str(F_name)],'  is : ', num2str(BpSCSO)]);

%% DBO    
[Best_pos,Best_score, DBO_cg_curve ] = DBO(PD_no,Max_iter,LB,UB,Dim,F_obj); % Call DBO
fprintf ('Best solution obtained by DBO: %s\n', num2str(Best_score,'%e  '));
display(['The best optimal value of the objective funciton found by DBO  for ' [num2str(F_name)],'  is : ', num2str(Best_pos)]);

% catch ME
%     % 记录错误信息
%     log_message(ME.message, 'ERROR');
%     log_message(sprintf('堆栈信息: %s at line %d', ME.stack(1).name, ME.stack(1).line), 'ERROR');
% 
% end


%% MaxA运行的算法，F_num*runs是独立运行得出的数据


%% 绘制进化曲线
CNT=20;
k=round(linspace(1,Max_iter,CNT)); %随机选CNT个点
% 注意：如果收敛曲线画出来的点很少，随机点很稀疏，说明点取少了，这时应增加取点的数量，100、200、300等，逐渐增加
% 相反，如果收敛曲线上的随机点非常密集，说明点取多了，此时要减少取点数量
iter=1:1:Max_iter;
figure('Position',[154   145   894   357]);
subplot(1,2,1);
func_plot_2005(F_name);     % Function plot
title('Parameter space')
xlabel('x_1');
ylabel('x_2');
zlabel([F_name,'( x_1 , x_2 )'])
subplot(1,2,2);       % Convergence plot
h1 = semilogy(iter(k),BWO_cg_curve(k),'y-+','linewidth',1);
hold on
h2 = semilogy(iter(k),SSA_cg_curve(k),'k-s','linewidth',1);
hold on
h3 = semilogy(iter(k),SCSO_cg_curve(k),'m-^','linewidth',1);
hold on
h4 = semilogy(iter(k),DBO_cg_curve(k),'b-*','linewidth',1);
hold on
h5 = semilogy(iter(k),GWO_cg_curve(k),'g-o','linewidth',1);
xlabel('Iteration#');
ylabel('Best fitness so far');
legend('BWO','SSA','SCSO','DBO','GWO');


%% 打印出评价指标

%% 结束日志
stop_logging();