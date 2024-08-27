% 清除命令窗口、关闭所有图形窗口、清空工作区
clc
clear
close all

%%
% 设置测试函数代号、原生动物（搜索代理）数量、最大迭代次数
Fun_name='F7'; % 测试函数名称，范围从'F1'至'F23'
SearchAgents=30;                      % 原生动物（种群成员）数量
Max_iterations=500;                  % 最大迭代次数

% 获取目标函数信息：下界、上界、维度、适应度计算函数
[lowerbound,upperbound,dimension,fitness]=fun_info(Fun_name);

% 使用APO算法计算给定问题的解
[Best_score,Best_pos,APO_func_curve]=APO_func(SearchAgents,Max_iterations,lowerbound,upperbound,dimension,fitness);

%%
% 显示APO算法找到的最佳解及其适应度值
fprintf('APO算法针对函数 %s 找到的最佳解为：\n', num2str(Fun_name));
disp(['[', num2str(Best_pos), ']']);
fprintf('对应的目标函数最优值为：%f\n', Best_score);

% 绘制参数空间图和收敛曲线
figure('Position',[454   445   694   297]); % 设置图形窗口大小
subplot(1,2,1);       % 创建子图1用于参数空间图
func_plot(Fun_name);     % 绘制测试函数图形
title('参数空间'); xlabel('x_1'); ylabel('x_2'); zlabel([Fun_name,'(x_1, x_2)']);

subplot(1,2,2);       % 创建子图2用于收敛曲线图
CNT=20;               % 选取用于绘制收敛曲线的点数
k=round(linspace(1,Max_iterations,CNT)); % 在迭代次数区间均匀抽取CNT个点
iter=1:1:Max_iterations;
% 根据函数特性选择合适的曲线绘制方式
if ~strcmp(Fun_name,'F16')&&~strcmp(Fun_name,'F9')&&~strcmp(Fun_name,'F11')
    semilogy(iter(k),APO_func_curve(k),'k-o','linewidth',1); % 对数y轴绘制
else
    plot(iter(k),APO_func_curve(k),'k-o','linewidth',1);      % 线性y轴绘制
end
xlabel('迭代次数');
ylabel('至今最佳适应度值');
legend('APO');