%%
% 2025-12-04
close all
clear all
clc

%% ==================== 1. 单位换算常数 ====================
m_to_ft = 3.28084;          % 长度：m -> ft
kg_to_slugs = 0.0685218;    % 质量：kg -> slugs
N_to_lbf = 0.224809;        % 力：N -> lbf
psf_to_pa = 47.8803;        % 压强：psf -> Pa

%% ==================== 2. 物理参数设置 ====================
% --- 船体参数 ---
mass_slugs = 10879.5; 
mass_kg = mass_slugs / kg_to_slugs; 

% 气垫底部面积 (每个舱室)
Area_cushion_ft2 = 800; 

% 围裙周长 (用于简化计算泄流面积 S = L * h)
L_per_cushion_ft = 140; 

% --- 初始状态 ---
z0_si = 0.0;        % 初始高度 (m) 
v0_si = 0;          % 初始速度 (m/s)

% --- 风机转速 (RPM) ---
N_FAN1 = 2000;      
N_FAN2 = 2000;

% --- 仿真控制 ---
dt = 0.01;          % 时间步长 (s)
t_total = 10;       % 模拟总时间 (s) 
max_iter_pressure = 50; 
tol_pressure = 10;    

%% ==================== 3. 变量初始化 ====================
N_steps = round(t_total / dt);
t_hist = zeros(N_steps, 1);
z_hist = zeros(N_steps, 1);      
v_hist = zeros(N_steps, 1);      
P_hist = zeros(N_steps, 6);      
F_lift_hist = zeros(N_steps, 1);
Q_pump_hist = zeros(N_steps, 1); 

z_curr_si = z0_si;
v_curr_si = v0_si;

% 压力初值猜测
P_guess = [0; 0; 0; 0; 0; 0]; 

fprintf('开始动力学仿真 ...\n');

%% ==================== 4. 时域迭代循环  ====================
for k = 1:N_steps
    t = (k-1) * dt;
    
    % --- 1. 几何计算 (转换到英制) ---
    z_curr_ft = z_curr_si * m_to_ft;
    
    % 计算垂直速度 (ft/s) 用于体积变化率
    % 向上运动 (v > 0) -> 体积增大 -> 需要填充空气 -> 压力降低
    % 向下运动 (v < 0) -> 体积减小 -> 空气被压缩 -> 压力升高
    dzdt_ft = v_curr_si * m_to_ft;
    
    eff_h_ft = max(z_curr_ft, 0.001); 
    
    % 计算泄流面积 S (ft²) 
    S_val = L_per_cushion_ft * eff_h_ft;
    S_vec_ft2 = [S_val; S_val; S_val; S_val];
    
    % --- 2. 求解气室压力  ---
    [P_solved_psf, converged] = solve_pressure(P_guess, S_vec_ft2, N_FAN1, N_FAN2, ...
                                               dzdt_ft, Area_cushion_ft2, ... 
                                               max_iter_pressure, tol_pressure);
    
    % 更新猜测值
    P_guess = P_solved_psf;
    
    % --- 3. 计算受力 ---
    F_lift_lbs = sum(P_solved_psf(1:4)) * Area_cushion_ft2;
    F_lift_N = F_lift_lbs / N_to_lbf;
    
    F_gravity_N = mass_kg * 9.81;
    
    F_z_N = F_lift_N - F_gravity_N;
    
    % --- 4. 运动积分 ---
    z_acc_si = F_z_N / mass_kg;
    v_next_si = v_curr_si + z_acc_si * dt;
    z_next_si = z_curr_si + v_next_si * dt;
    
   
    % --- 5. 数据记录 ---
    t_hist(k) = t;
    z_hist(k) = z_curr_si;
    v_hist(k) = v_curr_si;
    P_hist(k, :) = P_solved_psf' * psf_to_pa; 
    F_lift_hist(k) = F_lift_N;
    Q_pump_hist(k) = dzdt_ft * Area_cushion_ft2; % 记录体积变化率
    
    z_curr_si = z_next_si;
    v_curr_si = v_next_si;
end

%% ==================== 5. 结果绘图 ====================

% --- 窗口 1: 垫升过程  ---
figure('Name', 'Hovercraft Motion Response', 'Color', 'w');

subplot(2, 1, 1);
plot(t_hist, z_hist * 100, 'b-', 'LineWidth', 1.5);
ylabel('气隙高度 (cm)', 'FontSize', 10); 
xlabel('时间 (s)', 'FontSize', 10);
title('围裙离地间隙 ', 'FontSize', 11, 'FontWeight', 'bold');
grid on;


subplot(2, 1, 2);
plot(t_hist, v_hist, 'r-', 'LineWidth', 1.5);
ylabel('速度 (m/s)', 'FontSize', 10); 
xlabel('时间 (s)', 'FontSize', 10);
title('垂向速度 ', 'FontSize', 11, 'FontWeight', 'bold');
grid on;

% --- 窗口 2: 四个气室的压力分布 (Cushion Pressures) ---
figure('Name', 'Cushion Pressures', 'Color', 'w');

% 绘制 4 个气室的压力
for i = 1:4
    subplot(2, 2, i);
    plot(t_hist, P_hist(:, i), 'k-', 'LineWidth', 1.2);
    
    ylabel('压力 (Pa)', 'FontSize', 9); 
    xlabel('时间 (s)', 'FontSize', 9);
    title(sprintf('气室 %d 压强 ', i), 'FontSize', 10, 'FontWeight', 'bold');
    grid on;
    
end

% --- 窗口 3: 气垫垫升力和重力 ---
figure('Name', '气垫垫升力和重力', 'Color', 'w');
plot(t_hist, F_lift_hist/1000, 'k', 'LineWidth', 1.5);
yline(mass_kg * 9.81 / 1000, 'r--', '重力基准值');
xlabel('时间 (s)'); ylabel('力 (kN)');
grid on;




%% ==================== 自定义函数 ====================
function [P_solved, converged] = solve_pressure(P_init, S_vec, N1, N2, dzdt, Area, max_it, tol)
    P = P_init;
    converged = false;
    lambda = 0.7; 
    
    for i = 1:max_it
        
        J = compute_jacobian(P, S_vec, N1, N2, dzdt, Area);
        F = residual_function(P, S_vec, N1, N2, dzdt, Area);
        
        if norm(F) < tol
            converged = true;
            break;
        end
        
        delta = -J \ F;
        P = P + lambda * delta;
        P = max(P, 1.0); 
    end
    P_solved = P;
end


function J = compute_jacobian(P, S, N1, N2, dzdt, Area)
    eps_pert = 1e-4;
    n = 6;
    J = zeros(n, n);
    F0 = residual_function(P, S, N1, N2, dzdt, Area);
    for k = 1:n
        P_tmp = P;
        P_tmp(k) = P_tmp(k) + eps_pert;
        F_pert = residual_function(P_tmp, S, N1, N2, dzdt, Area);
        J(:, k) = (F_pert - F0) / eps_pert;
    end
end


function F = residual_function(P, S, N_FAN1, N_FAN2, dzdt, Area)
    % S: ft^2, P: psf, dzdt: ft/s, Area: ft^2
    
    % 风机 (Eq 49)
    Q_FAN1 = (-1280 * sqrt(abs(P(5) - 300)) * sign(P(5) - 300) - 31.6 * (P(5) - 300)) * N_FAN1 / 2000;
    Q_FAN2 = (-1280 * sqrt(abs(P(6) - 300)) * sign(P(6) - 300) - 31.6 * (P(6) - 300)) * N_FAN2 / 2000;
    
    % 内部流道 (Eq 48)
    Q_INC1 = 589 * sqrt(abs(P(5) - P(1))) * sign(P(5) - P(1));
    Q_INC2 = 589 * sqrt(abs(P(5) - P(2))) * sign(P(5) - P(2));
    Q_INC3 = 589 * sqrt(abs(P(6) - P(3))) * sign(P(6) - P(3));
    Q_INC4 = 589 * sqrt(abs(P(6) - P(4))) * sign(P(6) - P(4));
    
    Q_NOZ1 = -346 * sqrt(abs(P(5))) * sign(P(5));
    Q_NOZ2 = -346 * sqrt(abs(P(6))) * sign(P(6));
    
    Q_IC1 = 675 * sqrt(abs(P(4) - P(1))) * sign(P(4) - P(1));
    Q_IC2 = 338 * sqrt(abs(P(1) - P(2))) * sign(P(1) - P(2));
    Q_IC3 = 675 * sqrt(abs(P(2) - P(3))) * sign(P(2) - P(3));
    Q_IC4 = 338 * sqrt(abs(P(3) - P(4))) * sign(P(3) - P(4));
    
    % 围裙泄流 (保持简单的 S*sqrt(P) 形式)
    Q1 = -S(1) * 14.5 * sqrt(abs(P(1))) * sign(P(1));
    Q2 = -S(2) * 14.5 * sqrt(abs(P(2))) * sign(P(2));
    Q3 = -S(3) * 14.5 * sqrt(abs(P(3))) * sign(P(3));
    Q4 = -S(4) * 14.5 * sqrt(abs(P(4))) * sign(P(4));
    
    % 计算体积变化率 (QPUMP)
    % 公式 (39) QPUMP = d(Vol)/dt

    Q_PUMP_1 = Area * dzdt;
    Q_PUMP_2 = Area * dzdt;
    Q_PUMP_3 = Area * dzdt;
    Q_PUMP_4 = Area * dzdt;

    F = zeros(6, 1);

    F(5) = -Q_INC1 - Q_INC2 + Q_NOZ1 + Q_FAN1;
    F(6) = -Q_INC3 - Q_INC4 + Q_NOZ2 + Q_FAN2;
    
    % 气室方程 (Cushions)
    F(1) = Q_INC1 + Q_IC1 - Q_IC2 + Q1 - Q_PUMP_1;
    F(2) = Q_INC2 + Q_IC2 - Q_IC3 + Q2 - Q_PUMP_2;
    F(3) = Q_INC3 + Q_IC3 - Q_IC4 + Q3 - Q_PUMP_3;
    F(4) = Q_INC4 + Q_IC4 - Q_IC1 + Q4 - Q_PUMP_4;
end
