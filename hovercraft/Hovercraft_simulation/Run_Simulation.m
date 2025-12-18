%% 气垫船6DOF仿真 
clear; clc; close all;

%% 初始状态定义
% 状态向量 X = [x, y, z, phi, theta, psi, u, v, w, p, q, r]

% --- 位置与姿态 ---
x0 = 0;      % 北向位置 (m)
y0 = 0;      % 东向位置 (m)
z0 = 1.45;   % 垂向位移 (m)
phi0 = 0;    % 横摇角 (rad)
theta0 = 0;  % 纵摇角 (rad)
psi0 = deg2rad(0); % 初始艏向 (rad)，0度为正北

% --- 速度 ---
u0 = 0;      % 纵向速度 (m/s)
v0 = 0;      % 横向速度 (m/s)
w0 = 0;      % 垂向速度 (m/s)
p0 = 0;      % 横摇角速度 (rad/s)
q0 = 0;      % 纵摇角速度 (rad/s)
r0 = 0;      % 艏摇角速度 (rad/s)

% 组装初始状态向量
X0 = [x0, y0, z0, phi0, theta0, psi0, u0, v0, w0, p0, q0, r0];

%% 调用求解器 (使用 4阶龙格-库塔法 RK4)
% 定义时间步长和总时长
dt = 0.01;              
T_end = 200;            % 仿真结束时间
t = 0:dt:T_end;         % 生成时间向量
num_steps = length(t);  % 总步数

% 初始化存储数组
% 状态矩阵 sol: [行=时间步, 列=12个状态量]
sol = zeros(num_steps, 12); 
sol(1, :) = X0;         % 填入初始状态

% 压强矩阵 P_hist_Pa: [行=时间步, 列=4个气室]
P_hist_Pa = zeros(num_steps, 4);

% 计算初始时刻的压强 (仅用于记录)
[~, P_init] = model_jeff_b(t(1), X0'); 
P_hist_Pa(1, :) = P_init(:)';

fprintf('正在进行气垫船6自由度仿真...\n');

% RK4 主循环
for k = 1 : num_steps - 1
    % 当前时刻和状态
    t_curr = t(k);
    X_curr = sol(k, :)'; % 取出为列向量 (12x1)
    
    % --- RK4迭代 ---
    [dX1, ~] = model_jeff_b(t_curr, X_curr);
    [dX2, ~] = model_jeff_b(t_curr + 0.5*dt, X_curr + 0.5*dt*dX1);
    [dX3, ~] = model_jeff_b(t_curr + 0.5*dt, X_curr + 0.5*dt*dX2);
    [dX4, ~] = model_jeff_b(t_curr + dt, X_curr + dt*dX3);
    X_next = X_curr + (dt / 6) * (dX1 + 2*dX2 + 2*dX3 + dX4);
   
    sol(k+1, :) = X_next';
    
    [~, P_out] = model_jeff_b(t(k+1), X_next);
    P_hist_Pa(k+1, :) = P_out(:)';
    
end

fprintf('仿真完成。\n');

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

% 压强数据
P_FL = P_hist_Pa(:,1); % 前左
P_FR = P_hist_Pa(:,2); % 前右
P_RR = P_hist_Pa(:,3); 
P_RL = P_hist_Pa(:,4); 

P_static_si = 5218; % 绘图参考线

%% 绘图


%% 运动学状态
figure('Name', 'Kinematics', 'Color', 'w');

subplot(3,2,1); 
plot(t, u, 'LineWidth', 1.5); 
title('纵向速度 u (m/s)'); 
grid on; 
xlabel('时间 (s)');
ylim([-10, 10]);

subplot(3,2,2); 
plot(t, v, 'LineWidth', 1.5); 
title('横向速度 v (m/s)'); 
grid on; 
xlabel('时间 (s)');
ylim([-10, 10]);

subplot(3,2,3); 
plot(t, abs(1.45 - z_pos), 'LineWidth', 1.5); 
title('垫升高度(m)'); 
grid on; xlabel('时间 (s)');


subplot(3,2,4); plot(t, r_rate_deg, 'LineWidth', 1.5); 
title('艏摇率 r (deg/s)'); grid on; xlabel('时间 (s)');
ylim([-2, 2]);


subplot(3,2,5); 
plot(t, phi_deg, 'LineWidth', 1.5); 
title('横摇角 Roll (deg)'); grid on; 
xlabel('时间 (s)');
ylim([-1.5, 1.5]);

subplot(3,2,6); 
plot(t, theta_deg, 'LineWidth', 1.5); 
title('纵摇角 Pitch (deg)'); 
grid on; 
xlabel('时间 (s)');
ylim([-2, 2]);

%% 运动轨迹
figure('Name', 'Trajectory', 'Color', 'w');
plot(y_pos, x_pos, 'b-', 'LineWidth', 2); hold on;
plot(y_pos(1), x_pos(1), 'go', 'MarkerFaceColor', 'g', 'DisplayName', '起点'); 
plot(y_pos(end), x_pos(end), 'rs', 'MarkerFaceColor', 'r', 'DisplayName', '终点');

xlabel('东向East(m)'); ylabel('北向North(m)');
title('气垫船运动轨迹');
axis equal; grid on; legend show;

%% 侧滑角绘图
% 计算侧滑角 (Sideslip Angle)
beta_rad = atan2(v, u);      % 结果为弧度
beta_deg = rad2deg(beta_rad); % 转换为角度

%  绘制侧滑角图像
figure('Name', 'Sideslip Angle', 'Color', 'w');

plot(t, beta_deg, 'b-', 'LineWidth', 1.5); hold on;
yline(0, 'k--', 'LineWidth', 1, 'DisplayName', '无侧滑');

title('气垫船侧滑角 \beta (Sideslip Angle)');
xlabel('时间 (s)');
ylabel('侧滑角 (deg)');
legend('侧滑角 \beta', 'Location', 'best');
ylim([-1.5, 1.5]);
grid on;

%% 气垫压强绘图
figure('Name', 'Cushion Pressure Distribution', 'Color', 'w');

% 转换单位 Pa -> kPa 
P_scale = 1; 

% 前左气室
subplot(2, 2, 1);
plot(t, P_FL * P_scale, 'r-', 'LineWidth', 1.5); hold on;
yline(P_static_si * P_scale, 'k:', 'LineWidth', 1.2, 'DisplayName', '稳态');
ylabel('压强 (Pa)'); xlabel('时间 (s)');
title('前左气室 (Front-Left)');
grid on; axis tight;


% 前右气室
subplot(2, 2, 2);
plot(t, P_FR * P_scale, 'b-', 'LineWidth', 1.5); hold on;
yline(P_static_si * P_scale, 'k:', 'LineWidth', 1.2);
ylabel('压强 (Pa)'); xlabel('时间 (s)');
title('前右气室 (Front-Right)');
grid on; axis tight;


% 后左气室
subplot(2, 2, 3);
plot(t, P_RL * P_scale, 'r-', 'LineWidth', 1.5); hold on;
yline(P_static_si * P_scale, 'k:', 'LineWidth', 1.2);
ylabel('压强 (Pa)'); xlabel('时间 (s)');
title('后左气室 (Rear-Left)');
grid on; axis tight;
% ylim(y_limits);

% 后右气室
subplot(2, 2, 4);
plot(t, P_RR * P_scale, 'b-', 'LineWidth', 1.5); hold on;
yline(P_static_si * P_scale, 'k:', 'LineWidth', 1.2);
ylabel('压强 (Pa)'); xlabel('时间 (s)');
title('后右气室 (Rear-Right)');
grid on; axis tight;
% ylim(y_limits);

grid on;

