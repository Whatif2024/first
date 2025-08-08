% 系统初始化脚本
% 设置全局变量和系统参数

clear all; 
close all; 
clc;

fprintf('正在初始化系统...\n');

% 声明全局变量
global P_PV;
global P_WT;

% 设置24小时可再生能源出力数据
% 光伏出力曲线（白天高，晚上低）
hour = 1:24;
P_PV = zeros(1, 24);
for h = 1:24
    if h >= 6 && h <= 18  % 白天6点到18点有光伏
        P_PV(h) = 300 * sin(pi * (h - 6) / 12)^2;  % 正弦平方函数，峰值在中午
    else
        P_PV(h) = 0;  % 夜晚无光伏
    end
end

% 风电出力曲线（相对稳定，有一定波动）
P_WT = 150 + 50 * cos(2 * pi * hour / 24) + 30 * sin(4 * pi * hour / 24);
P_WT = max(0, P_WT);  % 确保非负

% 显示可再生能源出力曲线
figure(10);
subplot(2,1,1);
plot(hour, P_PV, 'r-o', 'LineWidth', 2);
xlabel('时间 (小时)');
ylabel('光伏出力 (kW)');
title('24小时光伏出力曲线');
grid on;

subplot(2,1,2);
plot(hour, P_WT, 'b-s', 'LineWidth', 2);
xlabel('时间 (小时)');
ylabel('风电出力 (kW)');
title('24小时风电出力曲线');
grid on;

fprintf('全局变量初始化完成:\n');
fprintf('光伏最大出力: %.2f kW (第 %d 小时)\n', max(P_PV), find(P_PV == max(P_PV)));
fprintf('风电平均出力: %.2f kW\n', mean(P_WT));
fprintf('风电出力范围: %.2f - %.2f kW\n', min(P_WT), max(P_WT));

fprintf('\n系统初始化完成！可以运行 main_game 开始博弈优化。\n');