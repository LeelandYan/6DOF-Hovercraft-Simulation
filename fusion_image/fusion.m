clc
clear all
close all

%% 1. 打开图形
f1 = openfig('1.fig');
f2 = openfig('2.fig');

%% 2. 提取数据
% 找到图中的线条对象
lines1 = findobj(f1, 'Type', 'line');
lines2 = findobj(f2, 'Type', 'line');

% 获取 X 和 Y 数据
x1 = get(lines1, 'XData');
y1 = get(lines1, 'YData');

x2 = get(lines2, 'XData');
y2 = get(lines2, 'YData');

%% 3. 在新窗口中重新画图
figure;
plot(x1, y1, 'b-', 'LineWidth', 1); hold on;
plot(x2, y2, 'r-', 'LineWidth', 2);

legend('标准滑模控制器', '自适应滑模控制器');
title('提取数据重绘');
ylim([0, 2600]);
grid on;