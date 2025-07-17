% plot_function.m
% 脚本功能：绘制所有测试函数图像 benchmark_func.m
% https://github.com/tsingke/CEC-Benchmark-Functions
% clear;close all

% 1.二维空间绘图时，生成多个自变量
for func_num =1:16
   
   if func_num==1 x=-100:5:100;y=x; %sphere_func
   elseif func_num==2 x=-100:5:100; y=x;%schwefel_102
   elseif func_num==3 x=-10:0.5:10; y=x;%schwefel_102_noise_func
   elseif func_num==4 x=-100:5:100; y=x;%schwefel_2_21
   elseif func_num==5 x=-10:0.5:10; y=x;%schwefel_2_22
   elseif func_num==6 x=-200:10:200;y=x; %high_cond_elliptic_func
   elseif func_num==7 x=-10:0.5:10; y=x; % step_func
   elseif func_num==8 x=-500:10:500; y=x;%Schwefel_func(multimodal)
   elseif func_num==9 x=-2.048:0.05:2.048; y=x;%rosenbrock_func
   elseif func_num==10 x=-1.28:0.05:1.28;y=x; %Quartic function
   elseif func_num==11 x=-600:30:600; y=x;%griewank_func
   elseif func_num==12 x=-32:1:32; y=x;% ackley_func
   elseif func_num==13 x=-5:0.1:5; y=x;% rastrigin_func
   elseif func_num==14 x=-5:0.1:5; y=x;% rastrigin_noncont                                                                                                                 astrigin_noncont
   elseif func_num==15 x=-5:0.1:5; y=x;% weierstrass
   elseif func_num==16 x=-10:0.2:10;y=x;% schaffer
   else
      break;
   end
   
   L=length(x);
   f=[];
   
   %2. 由自变量生成网格点，并计算各个网格点的函数值f(x,y)
   for i=1:L
       for j=1:L
           f(i,j)=benchmark_func([x(i),y(j)],func_num); % 一个点对应一个数值
       end
   end
   
   % 3. 根据网格点绘制空间图形
   figure(func_num)% 生成一份图形画布
   surfc(x,y,f);   % 在图形画布上绘制函数图像(带有等高线)surfC（X,Y,Z）
   
   xlabel('x'); 
   ylabel('y');
   zlabel('f(x,y)');
   
   % 4.输出并保存函数图像到指定路径下面 todo
   filename=strcat('E:\matlabProjects\demo\main\mymodel\result','F',num2str(func_num));
   saveas(gcf,filename,'png');
   close;% 关闭当前的输出
   
end
