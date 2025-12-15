%% 气垫船6DOF仿真 
clear; clc; close all;

%% 初始状态定义
% 状态向量 X = [x, y, z, phi, theta, psi, u, v, w, p, q, r]

% --- 位置与姿态 ---
x0 = 0;      % 北向位置 (m)
y0 = 0;      % 东向位置 (m)
z0 = 0;   % 垂向位移 (m)
phi0 = 0;    % 横摇角 (rad)
theta0 = 0;  % 纵摇角 (rad)
psi0 = deg2rad(0); % 初始艏向 (rad)，0度为正北

% --- 速度 ---
u0 = 0.1;      % 纵向速度 (m/s)
v0 = 0;      % 横向速度 (m/s)
w0 = 0;      % 垂向速度 (m/s)
p0 = 0;      % 横摇角速度 (rad/s)
q0 = 0;      % 纵摇角速度 (rad/s)
r0 = 0;      % 艏摇角速度 (rad/s)

% 组装初始状态向量
X0 = [x0, y0, z0, phi0, theta0, psi0, u0, v0, w0, p0, q0, r0];

%% 调用求解器
t_span = [0 200]; % 仿真时长(秒)
options = odeset('RelTol', 1e-4, 'AbsTol', 1e-6);

fprintf('进行气垫船6自由度仿真...\n');

[t, sol] = ode45(@model_jeff_b, t_span, X0, options);

%% 数据解包
% 位置 (m)
x_pos = sol(:,1);
y_pos = sol(:,2);
z_pos = sol(:,3); 

% 姿态 (弧度 -> 角度)
phi_deg   = rad2deg(sol(:,4));
theta_deg = rad2deg(sol(:,5));
% psi_deg   = rad2deg(sol(:,6)); 

% 速度 (m/s)
u = sol(:,7);
v = sol(:,8);
w = sol(:,9);

% 角速度 (rad/s -> deg/s)
r_rate_deg = rad2deg(sol(:,12));

%% --- 气垫压强---
% 初始化
num_steps = length(t);
P_hist_Pa = zeros(num_steps, 4); % 用于存储四个气室的压强
clear model_jeff_b; 

% 循环调用函数提取压强
for k = 1:num_steps
    t_curr = t(k);
    X_curr = sol(k, :)'; % 取出当前时刻的状态向量 
    
    [~, P_out] = model_jeff_b(t_curr, X_curr);
    
    % 记录数据 
    P_hist_Pa(k, :) = P_out(:)'; 
end


% 赋值给绘图变量
P_FL = P_hist_Pa(:,1); % 前左
P_FR = P_hist_Pa(:,2); % 前右
P_RL = P_hist_Pa(:,3); % 后左
P_RR = P_hist_Pa(:,4); % 后右

P_static_si = 5218; % 绘图参考线


%% 绘图
% --- 运动学状态 ---
figure('Name', 'Kinematics', 'Color', 'w');

subplot(2,2,1); 
plot(t, u, 'LineWidth', 1.5); 
title('纵向速度 u (m/s)'); 
grid on; 
xlabel('时间 (s)');


subplot(2,2,2); 
plot(t, v, 'LineWidth', 1.5); 
title('横向速度 v (m/s)'); 
grid on; 
xlabel('时间 (s)');
ylim([-10, 10]);

% subplot(3,2,3); 
% plot(t, (z_pos - z_pos(1)) * 100, 'LineWidth', 1.5); 
% title('气隙高度(cm)'); 
% grid on; xlabel('时间 (s)');


% subplot(3,2,3); plot(t, r_rate_deg, 'LineWidth', 1.5); 
% title('艏摇率 r (deg/s)'); grid on; xlabel('时间 (s)');
% ylim([-2, 2]);


subplot(2,2,3); 
plot(t, phi_deg, 'LineWidth', 1.5); 
title('横摇角 Roll (deg)'); grid on; 
xlabel('时间 (s)');
ylim([-1.5, 1.5]);

subplot(2,2,4); 
plot(t, theta_deg, 'LineWidth', 1.5); 
title('纵摇角 Pitch (deg)'); 
grid on; 
xlabel('时间 (s)');
% ylim([-2, 2]);

% --- 运动轨迹 ---
figure('Name', 'Trajectory', 'Color', 'w');
plot(y_pos, x_pos, 'b-', 'LineWidth', 2); hold on;
plot(y_pos(1), x_pos(1), 'go', 'MarkerFaceColor', 'g', 'DisplayName', '起点'); 
plot(y_pos(end), x_pos(end), 'rs', 'MarkerFaceColor', 'r', 'DisplayName', '终点');

xlabel('东向 East (m)'); ylabel('北向 North (m)');
title('气垫船运动轨迹 (JEFF-B)');
axis equal; grid on; legend show;

%气垫压强绘图
figure('Name', 'Cushion Pressure Distribution', 'Color', 'w');

% 转换单位 Pa -> kPa 
P_scale = 1; 

% --- 前左气室 ---
subplot(2, 2, 1);
plot(t, P_FL * P_scale, 'r-', 'LineWidth', 1.5); hold on;
yline(P_static_si * P_scale, 'k:', 'LineWidth', 1.2, 'DisplayName', '稳态');
ylabel('压强 (Pa)'); xlabel('时间 (s)');
title('前左气室 (Front-Left)');
grid on; axis tight;


% --- 前右气室 ---
subplot(2, 2, 2);
plot(t, P_FR * P_scale, 'b-', 'LineWidth', 1.5); hold on;
yline(P_static_si * P_scale, 'k:', 'LineWidth', 1.2);
ylabel('压强 (Pa)'); xlabel('时间 (s)');
title('前右气室 (Front-Right)');
grid on; axis tight;
% ylim(y_limits);

% --- 后左气室 ---
subplot(2, 2, 3);
plot(t, P_RL * P_scale, 'r-', 'LineWidth', 1.5); hold on;
yline(P_static_si * P_scale, 'k:', 'LineWidth', 1.2);
ylabel('压强 (Pa)'); xlabel('时间 (s)');
title('后左气室 (Rear-Left)');
grid on; axis tight;
% ylim(y_limits);

% --- 后右气室 ---
subplot(2, 2, 4);
plot(t, P_RR * P_scale, 'b-', 'LineWidth', 1.5); hold on;
yline(P_static_si * P_scale, 'k:', 'LineWidth', 1.2);
ylabel('压强 (Pa)'); xlabel('时间 (s)');
title('后右气室 (Rear-Right)');
grid on; axis tight;
% ylim(y_limits);

grid on;

