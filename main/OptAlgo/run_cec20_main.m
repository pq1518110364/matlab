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
PD_no = 30; % Number of search agents
Max_iter = 1000; % 最大迭代次数
CEC_f = 20;
% F_no = 10; 

% Default Values for TOC
Ntt=4;
Nt =3;

%% 调用函数，获得参数

[Function_name,F_num] = get_CEC_name(CEC_f);
% Best = zeros(F_num,MaxA);  %存储最优适应度值
% Mean = zeros(F_num,MaxA);  %存储平均适应度值
% Std = zeros(F_num,MaxA);   %存储适应度方差
% F_sum = zeros(F_num*3,MaxA); %最优值，平均值和方差过渡
% next_sum = zeros(F_num*3,MaxA); %存储最优值，平均值和方差
for a = 1:F_num    %运行函数 F_num 8 18-23 好像显示不了，需要具体排查，先跑1-3吧

    if CEC_f == 17 && a == 2 
        continue; % This function (F2) has been deleted
    end

    fprintf('以下为F%d函数的数据展示效果：\n',a );

    f_name = get_F_name(a);  %获得函数的序号
    [LB,UB,Dim,F_obj] = Function_name(f_name); %获得函数的边界
    % 最优适应度&最优位置&收敛曲线

    % 回旋镖气动椭圆优化算法(BAEO)    
    [BAEOBest_pos, BAEOBest_score, BAEO_cg_curve ] = BAEO(PD_no,Max_iter,LB,UB,Dim,F_obj); % Call BAEO

    % 龙卷风优化算法(TOC)    
    [TOCBest_pos, TOCBest_score, TOC_cg_curve ] = TOC(PD_no,Max_iter,LB,UB,Dim,F_obj,Ntt,Nt); % Call TOC

    % 状态优化算法(SBO)    
    [SBOBest_pos,SBOBest_score, SBO_cg_curve ] = HHO(PD_no,Max_iter,LB,UB,Dim,F_obj); % Call SBO
    
    % 哈里斯鹰优化算法(HHO)    
    [HHOBest_pos,HHOBest_score, HHO_cg_curve ] = HHO(PD_no,Max_iter,LB,UB,Dim,F_obj); % Call HHO

    % 沙猫群优化算法(SCSO)    
    [SCSOBest_pos,SCSOBest_score, SCSO_cg_curve ] = SCSO(PD_no,Max_iter,LB,UB,Dim,F_obj); % Call SCSO
    
    % 壮丽细尾鹩莺优化算法(SFOA) 论文中测试函数为cec2017与cec2020
    [SFOABest_pos,SFOAest_score, SFOA_cg_curve ] = SFOA(PD_no,Max_iter,LB,UB,Dim,F_obj); % Call SFOA

    % 不实野燕麦优化算法 论文中测试函数为cec2017与cec2020 AOOv4(fhd,dim,pop_size,iter_max,lb,ub,varargin)
    AOO_fhd = get_CEC_func_str(CEC_f);
    % 通过函数句柄 fhd 调用目标函数，传入参数 pos'（转置后的 pos）以及 varargin 中的所有参 e=feval(fhd,pos',varargin{:});
    % no F_obj -> a
    [AOOBest_pos,AOOest_score, AOO_cg_curve ] = AOO(AOO_fhd,Dim,PD_no,Max_iter,LB,UB,a); % Call AOO

    %% 绘制进化曲线
    CNT=30;
    k=round(linspace(1,Max_iter,CNT)); %随机选CNT个点
    iter=1:1:Max_iter;
    figure('Position',[300,300,800,330]);
    % subplot(1,2,1);
    % func_plot_2005(f_name);     % Function plot 需要替换原来的func_plot
    
    % title(f_name + '函数图');
    % xlabel('x_1');
    % ylabel('x_2');
    % zlabel([f_name,'( x_1 , x_2 )'])
    % subplot(1,2,2);       % Convergence plot
    
        semilogy(iter(k),BAEO_cg_curve(k),'Color', [1 0.5 0], 'Marker','+','LineStyle','-.', 'linewidth', 1);
        hold on
        semilogy(iter(k),TOC_cg_curve(k),'r-+','linewidth',1);
        hold on
        semilogy(iter(k),SBO_cg_curve(k),'y-+','linewidth',1);
        hold on
        semilogy(iter(k),HHO_cg_curve(k),'k-s','linewidth',1);
        hold on
        semilogy(iter(k),SCSO_cg_curve(k),'m-^','linewidth',1);
        hold on
        semilogy(iter(k),AOO_cg_curve(k),'b-*','linewidth',1);
        hold on
        semilogy(iter(k),SFOA_cg_curve(k),'g-o','linewidth',1);
        hold on
        
    grid on;
    
    title('各算法在'+f_name+'函数的迭代图'); % 添加标题
    xlabel('Iteration');
    ylabel('Best fitness so far');
    % box on
    legend('BAEO','TOC','SBO','HHO','SCSO','AOO','SFOA');
    % set (gcf,'position', [300,300,600,330])

    %% 打印出评价指标 30可以改成PD_no吗
    % 核心指标含义
    % Best（最优值）：算法多次运行中找到的最优解，值越小说明算法能找到的 “最好结果” 越好。
    % Mean（平均值）：多次运行结果的平均值，值越小说明算法的 “整体平均性能” 越稳定可靠。
    % STD（标准差）：反映多次运行结果的波动程度，值越小说明算法的 “稳定性” 越好（结果越一致）。

    % BAEO 的最佳适应度的Best、Mean、STD、Time
    BAEO_best_pos_list = zeros(30, 1);
    for i = 1:30
       [BAEOBest_pos,BAEOBest_score, BAEO_cg_curve ] = TOC(PD_no,Max_iter,LB,UB,Dim,F_obj,Ntt,Nt); % Call BAEO
         % 保存每次循环的结果
        BAEO_best_pos_list(i) = BAEOBest_pos;
    end
    % 计算 best、mean、STD、time
    BAEO_best = min(BAEO_best_pos_list);
    BAEO_mean_val = mean(BAEO_best_pos_list);
    BAEO_std_val = std(BAEO_best_pos_list);
    
    % 打印结果
    fprintf('以下为BAEO的数据展示：\n')
    fprintf('BAEO_Best: %.2e  ', BAEO_best);
    fprintf('BAEO_Mean: %.2e  ', BAEO_mean_val);
    fprintf('BAEO_STD: %.2e  \n', BAEO_std_val);

    % TOC 的最佳适应度的Best、Mean、STD、Time
    TOC_best_pos_list = zeros(30, 1);
    for i = 1:30
       [TOCBest_pos,TOCBest_score, TOC_cg_curve ] = TOC(PD_no,Max_iter,LB,UB,Dim,F_obj,Ntt,Nt); % Call SCSO
         % 保存每次循环的结果
        TOC_best_pos_list(i) = TOCBest_pos;
    end
    % 计算 best、mean、STD、time
    TOC_best = min(TOC_best_pos_list);
    TOC_mean_val = mean(TOC_best_pos_list);
    TOC_std_val = std(TOC_best_pos_list);
    
    % 打印结果
    fprintf('以下为TOC的数据展示：\n')
    fprintf('TOC_Best: %.2e  ', TOC_best);
    fprintf('TOC_Mean: %.2e  ', TOC_mean_val);
    fprintf('TOC_STD: %.2e  \n', TOC_std_val);

    % SBO 的最佳适应度的Best、Mean、STD、Time
    SBO_best_pos_list = zeros(30, 1);
    for i = 1:30
       [SBOBest_pos,SBOBest_score, SBO_cg_curve ] = SBO(PD_no,Max_iter,LB,UB,Dim,F_obj); % Call SCSO
         % 保存每次循环的结果
        SBO_best_pos_list(i) = SBOBest_pos;
    end
    % 计算 best、mean、STD、time
    SBO_best = min(SBO_best_pos_list);
    SBO_mean_val = mean(SBO_best_pos_list);
    SBO_std_val = std(SBO_best_pos_list);
    
    % 打印结果
    fprintf('以下为SBO的数据展示：\n')
    fprintf('SBO_Best: %.2e  ', SBO_best);
    fprintf('SBO_Mean: %.2e  ', SBO_mean_val);
    fprintf('SBO_STD: %.2e  \n', SBO_std_val);

    % 寻求HHO 的最佳适应度的Best、Mean、STD、Time
    HHO_best_pos_list = zeros(30, 1);
    for i = 1:30
       [HHOBest_pos,HHOBest_score, HHO_cg_curve ] = HHO(PD_no,Max_iter,LB,UB,Dim,F_obj); % Call SCSO
         % 保存每次循环的结果
        HHO_best_pos_list(i) = HHOBest_pos;
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

    % 寻求SCSO的最佳适应度的Best、Mean、STD、Time
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

    % 寻求SFOA的最佳适应度的Best、Mean、STD、Time
    AOO_best_pos_list = zeros(30, 1);
    for i = 1:30
        [AOOBest_pos,AOOest_score, AOO_cg_curve ] = AOO(AOO_fhd,Dim,PD_no,Max_iter,LB,UB,a); % Call SFOA
         % 保存每次循环的结果
        AOO_best_pos_list(i) = AOOBest_pos;
    end
    % 计算 best、mean、STD、time
    AOO_best = min(AOO_best_pos_list);
    AOO_mean_val = mean(AOO_best_pos_list);
    AOO_std_val = std(AOO_best_pos_list);

    % 打印结果
    fprintf('以下为AOO的数据展示：\n')
    fprintf('AOO_Best: %.2e  ', AOO_best);
    fprintf('AOO_Mean: %.2e  ', AOO_mean_val);
    fprintf('AOO_STD: %.2e  \n', AOO_std_val);

    % 寻求SFOA的最佳适应度的Best、Mean、STD、Time
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
    fprintf('\n')

end


%% 结束日志以及计时
stop_logging();
toc