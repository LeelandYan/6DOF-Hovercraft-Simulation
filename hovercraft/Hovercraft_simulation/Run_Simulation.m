%% JEFF-B 气垫船 6-DOF 仿真 (SI 单位制版本)
clear; clc; close all;

%% 1. 初始状态定义
% 状态向量 X = [x, y, z, phi, theta, psi, u, v, w, p, q, r]

% --- 位置与姿态 ---
x0 = 0;      % 北向位置 (m)
y0 = 0;      % 东向位置 (m)
z0 = 1.4;      % 垂向位移 (m)
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

%% 2. 调用求解器
t_span = [0 100]; % 仿真时长(秒)
options = odeset('RelTol', 1e-4, 'AbsTol', 1e-6);

fprintf('进行气垫船6自由度仿真...\n');

[t, sol] = ode45(@model_jeff_b, t_span, X0, options);

%% 3. 数据解包
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

%% --- 气垫压强逆向估算 ---
% 说明：基于计算出的合力和力矩，逆向推导四个气室的等效压强
% 参考文档：NAVTRAEQUIPCEN 73-C-0138-1 (JEFF-B Model)

% 1. 重新定义物理参数 (与 model_jeff_b 保持一致)
FT2M = 0.3048; LBF2N = 4.44822; SLUG2KG = 14.5939; g_SI = 9.80665;
m_kg = 10879.5 * SLUG2KG; 

% 刚度与阻尼系数 (从函数中提取)
K_heave = 20000 * LBF2N * (1/FT2M); C_heave = 5000  * LBF2N * (1/FT2M);
K_pitch = 4.0e7 * LBF2N * FT2M;     C_pitch = 5.0e6 * LBF2N * FT2M;
K_roll  = 1.5e7 * LBF2N * FT2M;     C_roll  = 2.0e6 * LBF2N * FT2M;

% 2. 几何参数估算 (基于 PDF 数据)
% PDF p.20 提到 SCAREA (Planform Area) = 3200 ft^2 [cite: 662]
% PDF p.20 提到 LCUSH (Cushion Length) = 77.0 ft [cite: 662]
Ac_ft2 = 3200;              % 气垫投影面积 (sq ft)
Lc_ft  = 77.0;              % 气垫长度 (ft)
Bc_ft  = Ac_ft2 / Lc_ft;    % 估算气垫宽度 (ft) ~ 41.5 ft

% 转换为 SI 单位
Ac_si = Ac_ft2 * (FT2M^2);
Lc_si = Lc_ft * FT2M;
Bc_si = Bc_ft * FT2M;

% 3. 计算静态压强 (Static Cushion Pressure)
% P = F / A. 静态时气垫升力等于重力。
% PDF p.10 提到 "cushion pressure of 109 psf" [cite: 151]
% 这里的计算值应该接近 109 psf (约 5200 Pa)
P_static_si = (m_kg * g_SI) / Ac_si; 

% 4. 循环计算每一时刻的压强分布
num_steps = length(t);
P_FL = zeros(num_steps, 1); % 前左 (Front-Left)
P_FR = zeros(num_steps, 1); % 前右 (Front-Right)
P_RL = zeros(num_steps, 1); % 后左 (Rear-Left)
P_RR = zeros(num_steps, 1); % 后右 (Rear-Right)

for i = 1:num_steps
    % 提取当前状态
    z_cur = sol(i, 3); 
    phi_cur = sol(i, 4);   % Roll
    theta_cur = sol(i, 5); % Pitch
    w_cur = sol(i, 9);     % Heave Vel
    p_cur = sol(i, 10);    % Roll Rate
    q_cur = sol(i, 11);    % Pitch Rate
    
    % --- A. 重算气垫合力与力矩 ---
    % 注意符号：NED系下 z 向下为正，但气垫力向上为负
    % 这里我们计算气垫提供的“反弹力/力矩”的绝对值大小
    
    % 垂向总动力 (Dynamic Heave Force)
    % 船体下沉 (z>0) -> 气垫被压缩 -> 产生向上的力 (负值) -> 压强增大
    % 我们只取变化部分：-(Kz + Cw)
    F_heave_dyn = -(K_heave * z_cur + C_heave * w_cur);
    
    % 横摇力矩 (Roll Moment)
    % 船体右倾 (phi>0) -> 右侧气垫压缩 -> 产生左旋力矩 (负值)
    M_roll_dyn = -(K_roll * phi_cur + C_roll * p_cur);
    
    % 纵摇力矩 (Pitch Moment)
    % 船体抬头 (theta>0) -> 尾部气垫压缩 -> 产生低头力矩 (负值)
    M_pitch_dyn = -(K_pitch * theta_cur + C_pitch * q_cur);
    
    % --- B. 将合力/力矩分配为压强变化 (Delta P) ---
    % 1. 均布压强变化 (由垂荡引起)
    dP_heave = -F_heave_dyn / Ac_si; 
    % 注意：F_heave_dyn 通常是负值（向上），负负得正，表示压强增加
    
    % 2. 纵向压强差 (由纵摇力矩引起)
    % 力矩 M = Force_diff * Arm. Arm 约等于 Lc/4 (前后中心距离)
    % M_pitch 负值(低头) -> 意味着尾部力大 -> 尾部压强大
    dP_pitch_dist = M_pitch_dyn / (Ac_si * Lc_si / 4); 
    % 符号约定：正值代表 前部比后部 大
    
    % 3. 横向压强差 (由横摇力矩引起)
    % M_roll 负值(左旋) -> 意味着右侧力大 -> 右侧压强大
    dP_roll_dist = - M_roll_dyn / (Ac_si * Bc_si / 4);
    % 符号约定：正值代表 左侧比右侧 大
    
    % --- C. 叠加计算四角总压强 ---
    % P_total = Static + Heave_delta +/- Pitch_delta +/- Roll_delta
    
    % 前左 (Front-Left): +Pitch (前), +Roll (左)
    P_FL(i) = P_static_si + dP_heave + dP_pitch_dist + dP_roll_dist;
    
    % 前右 (Front-Right): +Pitch (前), -Roll (右)
    P_FR(i) = P_static_si + dP_heave + dP_pitch_dist - dP_roll_dist;
    
    % 后左 (Rear-Left): -Pitch (后), +Roll (左)
    P_RL(i) = P_static_si + dP_heave - dP_pitch_dist + dP_roll_dist;
    
    % 后右 (Rear-Right): -Pitch (后), -Roll (右)
    P_RR(i) = P_static_si + dP_heave - dP_pitch_dist - dP_roll_dist;
end




%% 4. 绘图结果
% --- 运动学状态 ---
figure('Name', 'Kinematics', 'Color', 'w');

subplot(3,2,1); plot(t, u, 'LineWidth', 1.5); 
title('纵向速度 u (m/s)'); grid on; xlabel('时间 (s)');

subplot(3,2,2); plot(t, v, 'LineWidth', 1.5); 
title('横向速度 v (m/s)'); grid on; xlabel('时间 (s)');

subplot(3,2,3); plot(t, z_pos, 'LineWidth', 1.5); 
title('垂荡位移 Heave (m)'); grid on; xlabel('时间 (s)');

subplot(3,2,4); plot(t, r_rate_deg, 'LineWidth', 1.5); 
title('艏摇率 r (deg/s)'); grid on; xlabel('时间 (s)');

subplot(3,2,5); plot(t, -phi_deg, 'LineWidth', 1.5); 
title('横摇角 Roll (deg)'); grid on; xlabel('时间 (s)');

subplot(3,2,6); plot(t, theta_deg, 'LineWidth', 1.5); 
title('纵摇角 Pitch (deg)'); grid on; xlabel('时间 (s)');

% --- 运动轨迹 ---
figure('Name', 'Trajectory', 'Color', 'w');
plot(y_pos, x_pos, 'b-', 'LineWidth', 2); hold on;
plot(y_pos(1), x_pos(1), 'go', 'MarkerFaceColor', 'g', 'DisplayName', '起点'); 
plot(y_pos(end), x_pos(end), 'rs', 'MarkerFaceColor', 'r', 'DisplayName', '终点');

xlabel('东向 East (m)'); ylabel('北向 North (m)');
title('气垫船运动轨迹 (JEFF-B)');
axis equal; grid on; legend show;

%% --- 气垫压强绘图---
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

