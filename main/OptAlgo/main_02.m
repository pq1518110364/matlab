clear all 
close all
clc
%% 将当前工作目录以及其子目录添加到 MATLAB 搜索路径中，以便可以找到相关函数和文件
setup_paths();

%% 启动日志，diary以及计时
tic %时间计时器启动
start_logging();

%% 数据读取与处理，测试cec函数时这一操作没必要

%% 设置参数
PD_no = 30; % Number of sand cat 沙猫优化算法
Max_iter = 500; % 最大迭代次数
CEC_f = 05;
% F_no = 10;

%% 调用函数，获得参数

[Function_name,F_num] = get_CEC_name(CEC_f);
% Best = zeros(F_num,MaxA);  %存储最优适应度值
% Mean = zeros(F_num,MaxA);  %存储平均适应度值
% Std = zeros(F_num,MaxA);   %存储适应度方差
% F_sum = zeros(F_num*3,MaxA); %最优值，平均值和方差过渡
% next_sum = zeros(F_num*3,MaxA); %存储最优值，平均值和方差
for a = 1:1    %运行函数 F_num 8 18-23 好像显示不了，需要具体排查，先跑1-3吧
    f_name = get_F_name(a);  %获得函数的序号
    [LB,UB,Dim,F_obj] = Function_name(f_name); %获得函数的边界
    % 最优适应度&最优位置&收敛曲线
    % 白鲸优化算法(BWO)
    [BWOBest_pos,BWOBest_score, BWO_cg_curve ] = BWO(PD_no,Max_iter,LB,UB,Dim,F_obj); % Call BWO
    % 鲸鱼优化算法（WOA）
    [WOABest_pos,WOABest_score, WOA_cg_curve ] = WOA(PD_no,Max_iter,LB,UB,Dim,F_obj); % Call WOA
    % 正余弦优化算法（SCA）   
    [Alpha_score,Alpha_pos,SCA_cg_curve] = SCA(PD_no,Max_iter,LB,UB,Dim,F_obj); % Call SCA
    % 哈里斯鹰优化算法(HHO)
    [HHO_Score,HHO_pos,HHO_cg_curve]=HHO(PD_no,Max_iter,LB,UB,Dim,F_obj);% Call HHO
    % 沙猫群优化算法(SCSO)    
    [SCSOBest_pos,SCSOBest_score, SCSO_cg_curve ] = SCSO(PD_no,Max_iter,LB,UB,Dim,F_obj); % Call SSA
    % 蜣螂优化算法(DBO)
    [DBOBest_pos,DBOBest_score, DBO_cg_curve ] = DBO(PD_no,Max_iter,LB,UB,Dim,F_obj); % Call DBO
    % 壮丽细尾鹩莺优化算法(SFOA)
    [SFOABest_pos,SFOAest_score, SFOA_cg_curve ] = SFOA(PD_no,Max_iter,LB,UB,Dim,F_obj); % Call SFOA

    %% 绘制进化曲线
    CNT=20;
    k=round(linspace(1,Max_iter,CNT)); %随机选CNT个点
    iter=1:1:Max_iter;
    figure('Position',[154   145   894   357]);
    subplot(1,2,1);
    func_plot_2005(f_name);     % Function plot 需要替换原来的func_plot
    
    title(f_name + '函数图');
    xlabel('x_1');
    ylabel('x_2');
    zlabel([f_name,'( x_1 , x_2 )'])
    subplot(1,2,2);       % Convergence plot
    
        semilogy(iter(k),BWO_cg_curve(k),'Color', [1 0.5 0], 'Marker','+','LineStyle','-.', 'linewidth', 1);
        hold on
        semilogy(iter(k),WOA_cg_curve(k),'r-+','linewidth',1);
        hold on
        semilogy(iter(k),SCA_cg_curve(k),'y-+','linewidth',1);
        hold on
        semilogy(iter(k),HHO_cg_curve(k),'k-s','linewidth',1);
        hold on
        semilogy(iter(k),SCSO_cg_curve(k),'m-^','linewidth',1);
        hold on
        semilogy(iter(k),DBO_cg_curve(k),'b-*','linewidth',1);
        hold on
        semilogy(iter(k),SFOA_cg_curve(k),'g-o','linewidth',1);
        hold on
        
    grid on;
    
    title('各算法在'+f_name+'函数的迭代图'); % 添加标题
    xlabel('Iteration');
    ylabel('Best fitness so far');
    % box on
    % legend('BWO','WOA','SCA','HHO','SCSO','DBO','PTWDBO')
    legend('BWO','WOA','SCA','HHO','SCSO','DBO','SFOA')
    % set (gcf,'position', [300,300,600,330])

    %% 打印出评价指标
    % 核心指标含义
    % Best（最优值）：算法多次运行中找到的最优解，值越小说明算法能找到的 “最好结果” 越好。
    % Mean（平均值）：多次运行结果的平均值，值越小说明算法的 “整体平均性能” 越稳定可靠。
    % STD（标准差）：反映多次运行结果的波动程度，值越小说明算法的 “稳定性” 越好（结果越一致）。
    
    %% 寻求BWO的最佳适应度的Best、Mean、STD、Time
    BWO_best_pos_list = zeros(30, 1);
    for i = 1:30
        [BWOBest_pos,BWOBest_score, BWO_cg_curve ] = BWO(PD_no,Max_iter,LB,UB,Dim,F_obj); % Call BWO
         % 保存每次循环的结果
        BWO_best_pos_list(i) = BWOBest_pos;
    end
    % 计算 best、mean、STD、time
    BWO_best = min(BWO_best_pos_list);
    BWO_mean_val = mean(BWO_best_pos_list);
    BWO_std_val = std(BWO_best_pos_list);
    
    % 打印结果
    fprintf('以下为BWO的数据展示：\n')
    fprintf('BWO_Best: %.2e  ', BWO_best);
    fprintf('BWO_Mean: %.2e  ', BWO_mean_val);
    fprintf('BWO_STD: %.2e  \n', BWO_std_val);

    %% 寻求WOA的最佳适应度的Best、Mean、STD、Time
    WOA_best_pos_list = zeros(30, 1);
    for i = 1:30
        [WOABest_pos,WOABest_score, WOA_cg_curve ] = WOA(PD_no,Max_iter,LB,UB,Dim,F_obj); % Call WOA
         % 保存每次循环的结果
        WOA_best_pos_list(i) = WOABest_pos;
    end
    % 计算 best、mean、STD、time
    WOA_best = min(WOA_best_pos_list);
    WOA_mean_val = mean(WOA_best_pos_list);
    WOA_std_val = std(WOA_best_pos_list);
    
    % 打印结果
    fprintf('以下为WOA的数据展示：\n')
    fprintf('WOA_Best: %.2e  ', WOA_best);
    fprintf('WOA_Mean: %.2e  ', WOA_mean_val);
    fprintf('WOA_STD: %.2e  \n', WOA_std_val);
    
    %% 寻求SCA的最佳适应度的Best、Mean、STD、Time
    SCA_best_pos_list = zeros(30, 1);
    for i = 1:30
        [Alpha_score,Alpha_pos,SCA_cg_curve] = SCA(PD_no,Max_iter,LB,UB,Dim,F_obj); % Call GWO
         % 保存每次循环的结果
        SCA_best_pos_list(i) = Alpha_score;
    end
    % 计算 best、mean、STD、time
    SCA_best = min(SCA_best_pos_list);
    SCA_mean_val = mean(SCA_best_pos_list);
    SCA_std_val = std(SCA_best_pos_list);
    
    % 打印结果
    fprintf('以下为SCA的数据展示：\n')
    fprintf('SCA_Best: %.2e  ', SCA_best);
    fprintf('SCA_Mean: %.2e  ', SCA_mean_val);
    fprintf('SCA_STD: %.2e  \n', SCA_std_val);
    
    %% 寻求HHO的最佳适应度的Best、Mean、STD、Time
    HHO_best_pos_list = zeros(30, 1);
    for i = 1:30
        [HHO_Score,HHO_pos,HHO_cg_curve]=HHO(PD_no,Max_iter,LB,UB,Dim,F_obj);
         % 保存每次循环的结果
        HHO_best_pos_list(i) = HHO_Score;
    end
    % 计算 best、mean、STD、time
    HHO_best = min(HHO_best_pos_list);
    HHO_mean_val = mean(HHO_best_pos_list);
    HHO_std_val = std(HHO_best_pos_list);
    
    % 打印结果
    fprintf('以下为HHO的数据展示：\n')
    fprintf('HHO_Best: %.2e  ', HHO_best);
    fprintf('HHO_Mean: %.2e  ', HHO_mean_val);
    fprintf('HHO_STD: %.2e  \n', HHO_std_val);
    
    %% 寻求SCSO的最佳适应度的Best、Mean、STD、Time
    SCSO_best_pos_list = zeros(30, 1);
    for i = 1:30
       [SCSOBest_pos,SCSOBest_score, SCSO_cg_curve ] = SCSO(PD_no,Max_iter,LB,UB,Dim,F_obj); % Call SCSO
         % 保存每次循环的结果
        SCSO_best_pos_list(i) = SCSOBest_pos;
    end
    % 计算 best、mean、STD、time
    SCSO_best = min(SCSO_best_pos_list);
    SCSO_mean_val = mean(SCSO_best_pos_list);
    SCSO_std_val = std(SCSO_best_pos_list);
    
    % 打印结果
    fprintf('以下为SCSO的数据展示：\n')
    fprintf('SCSO_Best: %.2e  ', SCSO_best);
    fprintf('SCSO_Mean: %.2e  ', SCSO_mean_val);
    fprintf('SCSO_STD: %.2e  \n', SCSO_std_val);
    
    %% 寻求DBO的最佳适应度的Best、Mean、STD、Time
    DBO_best_pos_list = zeros(30, 1);
    for i = 1:30
        [DBOBest_pos,DBOBest_score, DBO_cg_curve ] = DBO(PD_no,Max_iter,LB,UB,Dim,F_obj); % Call DBO
         % 保存每次循环的结果
        DBO_best_pos_list(i) = DBOBest_pos;
    end
    % 计算 best、mean、STD、time
    DBO_best = min(DBO_best_pos_list);
    DBO_mean_val = mean(DBO_best_pos_list);
    DBO_std_val = std(DBO_best_pos_list);
    
    % 打印结果
    fprintf('以下为DBO的数据展示：\n')
    fprintf('DBO_Best: %.2e  ', DBO_best);
    fprintf('DBO_Mean: %.2e  ', DBO_mean_val);
    fprintf('DBO_STD: %.2e  \n', DBO_std_val);

    %% 寻求DBO的最佳适应度的Best、Mean、STD、Time
    SFOA_best_pos_list = zeros(30, 1);
    for i = 1:30
        [SFOABest_pos,SFOAest_score, SFOA_cg_curve ] = SFOA(PD_no,Max_iter,LB,UB,Dim,F_obj); % Call SFOA
         % 保存每次循环的结果
        SFOA_best_pos_list(i) = SFOABest_pos;
    end
    % 计算 best、mean、STD、time
    SFOA_best = min(SFOA_best_pos_list);
    SFOA_mean_val = mean(SFOA_best_pos_list);
    SFOA_std_val = std(SFOA_best_pos_list);
    
    % 打印结果
    fprintf('以下为SFOA的数据展示：\n')
    fprintf('SFOA_Best: %.2e  ', SFOA_best);
    fprintf('SFOA_Mean: %.2e  ', SFOA_mean_val);
    fprintf('SFOA_STD: %.2e  \n', SFOA_std_val);

end




%% 结束日志以及计时
stop_logging();
toc