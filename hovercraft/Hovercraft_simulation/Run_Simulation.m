%% 气垫船 6-DOF 仿真 
clear; clc; close all;

%% 1. 初始状态定义
% 状态向量 X = [x, y, z, phi, theta, psi, u, v, w, p, q, r]

% --- 位置与姿态 ---
x0 = 0;      % 北向位置 (m)
y0 = 0;      % 东向位置 (m)
z0 = 0.15;      % 垂向位移 (m)
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
t_span = [0 300]; % 仿真时长(秒)
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

%% --- 气垫压强复现 (Post-Processing) ---
% 说明：利用解算出的运动状态(sol)，重新调用物理模型来计算每一时刻的真实压强。

% 1. 准备常数 (必须与 model_jeff_b 保持一致)
FT2M = 0.3048; M2FT = 1/FT2M;
L_cush = 38.5 * FT2M; W_cush = 17.5 * FT2M;
H_base_m = 5.0 * FT2M; 
L_per_cushion_ft = 140;
P_equilibrium_psf = 109;
C_SKRT = 0.0112; TC = 8.0;
N_FAN = [1800, 1800]; % 风机转速

% 气室中心坐标
Pos_cush = [L_cush/2, W_cush/2; L_cush/2, -W_cush/2; -L_cush/2, W_cush/2; -L_cush/2, -W_cush/2];

% 2. 初始化记录数组
num_steps = length(t);
P_hist_Pa = zeros(num_steps, 4); % 记录4个气室压强

% 3. 状态重演循环
% 我们需要模拟 model_jeff_b 中的围裙动态积分过程
SD_curr_ft = ones(4,1) * 4.5; % 初始围裙长度 (ft)
P_last_psf = zeros(6,1);      % 初始内部压强状态

for k = 1:num_steps
    % A. 获取当前时刻状态
    z_ned = sol(k, 3); 
    phi = sol(k, 4);   
    theta = sol(k, 5); 
    w_vel = sol(k, 9); 
    
    % 计算时间步长 (ode45是变步长的)
    if k > 1
        dt = t(k) - t(k-1);
    else
        dt = 0;
    end
    
    S_vec_ft2 = zeros(4,1);
    
    % B. 重复 model_jeff_b 中的物理计算逻辑
    for i = 1:4
        % 计算气室中心物理高度
        h_local_m = (H_base_m - z_ned) + Pos_cush(i,1)*sin(theta) - Pos_cush(i,2)*sin(phi);
        h_hull_ft = h_local_m * M2FT;
        
        % 围裙动力学积分
        Y_n = P_equilibrium_psf - P_last_psf(i);
        Y_n = max(min(Y_n, 66.9), -66.9);
        SDPRO = 4.5 + C_SKRT * Y_n - 8.325e-7 * (Y_n^3);
        
        if k > 1 % 只有时间流逝才积分
            SDDT = (SDPRO - SD_curr_ft(i)) * TC;
            SD_curr_ft(i) = SD_curr_ft(i) + SDDT * dt;
        end
        
        % 计算气隙和面积
        SAW = h_hull_ft - SD_curr_ft(i);
        CLR = max(SAW, 1e-4);
        S_vec_ft2(i) = L_per_cushion_ft * CLR;
    end
    
    % C. 调用压力求解器 (注意 dzdt 的符号修正: -w)
    dzdt_ft = -w_vel * M2FT; 
    [P_pa_out, P_state_out] = calc_cushion_pressure(N_FAN, S_vec_ft2', dzdt_ft, P_last_psf);
    
    % D. 记录数据
    P_hist_Pa(k, :) = P_pa_out;
    P_last_psf = P_state_out; % 更新状态供下一步使用
end

% 赋值给绘图变量
P_FL = P_hist_Pa(:,1);
P_FR = P_hist_Pa(:,2);
P_RL = P_hist_Pa(:,3);
P_RR = P_hist_Pa(:,4);
P_static_si = 5218; % 绘图参考线 (109 psf)


%% 4. 绘图结果
% --- 运动学状态 ---
figure('Name', 'Kinematics', 'Color', 'w');

subplot(3,2,1); plot(t, u, 'LineWidth', 1.5); 
title('纵向速度 u (m/s)'); grid on; xlabel('时间 (s)');

subplot(3,2,2); plot(t, v, 'LineWidth', 1.5); 
title('横向速度 v (m/s)'); grid on; xlabel('时间 (s)');

subplot(3,2,3); 
plot(t, z_pos - z_pos(1), 'LineWidth', 1.5); 
title('升沉位移(m)'); 
grid on; xlabel('时间 (s)');


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

