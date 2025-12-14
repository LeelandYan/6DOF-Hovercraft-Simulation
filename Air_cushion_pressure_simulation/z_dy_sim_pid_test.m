close all
clear all
clc

%% ==================== 1、单位换算 ====================
m_to_ft = 3.28084;          % 长度：m -> ft
ft_to_m = 1/m_to_ft;
kg_to_slugs = 0.0685218;    % 质量：kg -> slugs
N_to_lbf = 0.224809;        % 力：N -> lbf
psf_to_pa = 47.8803;        % 压强：psf -> Pa

%% ==================== 2、参数设置 ====================
% --- 船体参数 ---
mass_kg = 158774;           % 输入质量 (kg) [对应原文 10879.5 slugs]
mass_slugs = mass_kg * kg_to_slugs; 

% --- 气垫几何参数 ---
Area_cushion_m2 = 74.32;    % 单个气室面积 (m^2) [对应原文 800 ft^2]
Area_cushion_ft2 = Area_cushion_m2 * m_to_ft^2; 

% 围裙周长
L_per_cushion_m = 42.67;    % 输入周长 (m) [对应原文 140 ft]
L_per_cushion_ft = L_per_cushion_m * m_to_ft; 

% --- 围裙参数  ---
skirt_length_stable_m = 1.37;    % 平衡状态围裙长(m) [对应原文 4.5 ft] 
skirt_length_stable_ft = skirt_length_stable_m * m_to_ft; 

% 平衡压强
P_stable_pa = 5218.9;  % 输入压强 (Pa) [对应原文 109 psf]
P_stable_psf = P_stable_pa / psf_to_pa;

% 围裙经验系数 
C_skirt = 0.0112;            % 围裙刚度系数
TC = 8.0;                    % 围裙响应时间常数/增益 

% --- PID参数 ---
z_target = 1.45;  % 目标高度    

Kp = 1500;           % 比例增益：
Ki = 20;             % 积分增益：
Kd = 800;           % 微分增益：

% --- 风机参数设置 ---
N_max = 2200;     % 最大转速
N_min = 0;        % 最小怠速
tau_fan = 0.5;    % 风机响应时间常数 (秒)

% --- 初始状态 ---
z0_si = 0.0;                % 初始船底高度 (m)
v0_si = 0;                  % 初始垂直速度 (m/s)

% --- 仿真参数 ---
dt = 0.01;          
t_total = 15;               % 模拟时长 (s)
max_iter_pressure = 50;     % 最大迭代次数
tol_pressure = 10;          % 最大容差



%% ==================== 3、变量初始化 ====================
N_steps = round(t_total / dt);
t_hist = zeros(N_steps, 1);
z_hist = zeros(N_steps, 1);         % 船壳高度 (m)
gap_hist = zeros(N_steps, 1);       % 围裙底部气隙 (m)
skirt_len_hist = zeros(N_steps, 1); % 围裙长度 (m)
v_hist = zeros(N_steps, 1);       
P_hist = zeros(N_steps, 6);         % 压强 (Pa)
F_lift_hist = zeros(N_steps, 1);    % 升力 (N)
acc_hist = zeros(N_steps, 1);

z_curr_si = z0_si;
v_curr_si = v0_si;

% 围裙长度初始化，假设初始时刻围裙处于自然伸展状态
SD_curr = [skirt_length_stable_ft; skirt_length_stable_ft; skirt_length_stable_ft; skirt_length_stable_ft]; 

% 压强初值
P_last_iter = [0; 0; 0; 0; 0; 0]; 


% PID 误差记录
error_prev = 0;      % 上一时刻误差 
error_pprev = 0;     % 上上一时刻误差

% 风机转速初始化
N_curr = 0;                 % 初始转速 (RPM)
N_cmd = N_curr;             % 命令转速
N_hist = zeros(N_steps, 1); % 历史转速



fprintf('气垫船动力学仿真...\n');

%% ==================== 4、迭代循环 ====================
for k = 1:N_steps
    t = (k-1) * dt;

    % ================== PID计算 ==================
    
    % 1、计算误差 (参考输入 - 实际反馈)
    z_feedback = z_hist(max(1, k-1)); 
    error = z_target - z_feedback;
    
    % 2、PID 计算公式
    % dN = Kp*(e(k) - e(k-1)) + Ki*e(k) + Kd*(e(k) - 2e(k-1) + e(k-2))
    delta_N = Kp * (error - error_prev) + ...
              Ki * error + ...
              Kd * (error - 2*error_prev + error_pprev);
          
    % 3、更新误差历史
    error_pprev = error_prev;
    error_prev = error;
    
    % 4、计算命令转速
    N_cmd = N_cmd + delta_N;
    
    % 5、限制命令转速范围
    if N_cmd > N_max, N_cmd = N_max; end
    if N_cmd < N_min, N_cmd = N_min; end
    
    % 6、加入风机惯性 (一阶滞后响应)
    % N_actual(k) = N_actual(k-1) + (dt/tau) * (N_cmd - N_actual(k-1))
    N_curr = N_curr + (dt / tau_fan) * (N_cmd - N_curr);
    
    % 赋值转速
    N_FAN1 = N_curr;
    N_FAN2 = N_curr;
    
    % 记录数据
    N_hist(k) = N_curr;
    % ====================================

    
    % --- 1、单位转换  ---
    z_curr_ft = z_curr_si * m_to_ft;
    dzdt_ft = v_curr_si * m_to_ft;
    
    % --- 2、围裙计算 ---
    if k == 1
        P_cushion_last = [0;0;0;0];
    else
        % 将历史记录的 Pa 转回 psf 
        P_cushion_last = P_hist(k-1, 1:4)' / psf_to_pa; 
    end
    
    S_vec_ft2 = zeros(4,1);
    current_gaps_ft = zeros(4,1);
    
    for i = 1:4
        % 计算压强偏差 Yn (单位: psf)
        % 压强越高，Yn越负，围裙向上收缩
        Y_n = P_stable_psf - P_cushion_last(i); 
        
        % 约束 Yn 范围
        if Y_n > 66.9, Y_n = 66.9; end
        if Y_n < -66.9, Y_n = -66.9; end
        
        % 计算此刻围裙长度
        SDPRO = 4.5 + C_skirt * Y_n - 8.325e-7 * (Y_n^3);
        
        % 围裙响应
        SDDT = (SDPRO - SD_curr(i)) * TC; 
        SD_curr(i) = SD_curr(i) + SDDT * dt;
        
        % 计算围裙泄流高度 (ft)
        % SAW(围裙泄流高度) = 船壳高度 - 围裙长度
        SAW = z_curr_ft - SD_curr(i);
        
        % 处理接触状态
        if SAW <= 0
            CLR = 1e-6; % 接触密封
        else
            CLR = SAW;  % 产生气隙
        end
        
        current_gaps_ft(i) = CLR;
        
        % 计算裙底泄流面积 (ft^2)
        S_vec_ft2(i) = L_per_cushion_ft * CLR;
    end
    
    % --- 3、求解气室压强  ---
    [P_solved_psf, converged] = solve_pressure(P_last_iter, S_vec_ft2, N_FAN1, N_FAN2, ...
                                               dzdt_ft, Area_cushion_ft2, ... 
                                               max_iter_pressure, tol_pressure);
    
    P_last_iter = P_solved_psf; 
    
    % --- 4、计算受力  ---
    F_lift_lbs = sum(P_solved_psf(1:4)) * Area_cushion_ft2;
    
    F_cursion_N = F_lift_lbs / N_to_lbf;
    F_gravity_N = mass_kg * 9.81;
    F_z_N = F_cursion_N - F_gravity_N;
    
    % --- 5、运动积分 ---
    z_acc_si = F_z_N / mass_kg;
    v_next_si = v_curr_si + z_acc_si * dt;
    z_next_si = z_curr_si + v_next_si * dt;
    
    if z_next_si < 0, z_next_si = 0; v_next_si = 0; end
    
    % --- 6、数据记录  ---
    t_hist(k) = t;
    z_hist(k) = z_curr_si;                      % 船壳高度 (m)
    gap_hist(k) = current_gaps_ft(1) * ft_to_m; % 气隙高度 (m)
    skirt_len_hist(k) = SD_curr(1) * ft_to_m;   % 围裙长度 (m)
    v_hist(k) = v_curr_si;
    P_hist(k, :) = P_solved_psf' * psf_to_pa;   % 压强 (Pa)
    F_lift_hist(k) = F_cursion_N;
    acc_hist(k) = z_acc_si;
    
    z_curr_si = z_next_si;
    v_curr_si = v_next_si;
end

%% ==================== 结果打印 ====================
fprintf('\n----------------------------------------\n');
fprintf('仿真结束。\n ');
fprintf('----------------------------------------\n');
for i = 1:4
    fprintf('  气室 %d 压强: %10.2f Pa\n', i, P_hist(end, i));
end
fprintf('----------------------------------------\n');
fprintf('  船壳最终高度: %10.2f m\n', z_hist(end));
fprintf('  风机转速值: %d RPM\n', N_FAN1);
fprintf('\n');


%% ==================== 5、结果绘图 ====================
% --- 图1：船体运动与几何状态 ---
figure('Name', 'Heave Dynamics', 'NumberTitle', 'off');
plot(t_hist, z_hist, 'b-', 'LineWidth', 2); hold on;
plot(t_hist, gap_hist, 'r--', 'LineWidth', 1.5);
% plot(t_hist, skirt_len_hist, 'g:', 'LineWidth', 1.5);
legend('船壳离地高度 ', '围裙底部气隙高度', 'Location', 'Best');
ylabel('高度 (m)'); xlabel('时间 (s)');
title('垫升过程仿真');
grid on;

% --- 图2：四个气室压强 ---
figure('Name', 'Cushion Pressures', 'NumberTitle', 'off');
for i = 1:4
    subplot(2, 2, i);
    plot(t_hist, P_hist(:, i), 'k-', 'LineWidth', 1.5);
    ylabel('压强 (Pa)'); xlabel('时间 (s)');
    title(['气室 ' num2str(i) ' 压强变化']);
    grid on;
    axis tight; 
    ylim_curr = ylim;
    ylim([ylim_curr(1)*0.9, ylim_curr(2)*1.1]);
end

% --- 图3：PID 控制效果 ---
% figure('Name', 'PID Control Performance', 'NumberTitle', 'off');
% 
% subplot(2,1,1);
% plot(t_hist, z_hist, 'b-', 'LineWidth', 2); hold on;
% yline(z_target, 'r--');
% ylabel('高度 (m)'); grid on;
% legend('实际高度', '目标高度');
% title('高度响应');
% 
% subplot(2,1,2);
% plot(t_hist, N_hist, 'k-', 'LineWidth', 1.5); hold on;
% yline(N_max, 'r:', 'Max RPM');
% yline(N_min, 'r:', 'Min RPM');
% ylabel('风机转速 (RPM)'); xlabel('时间 (s)');
% title('风机转速调节');
% grid on;

% --- 图4：升沉速度与加速度 ---
figure('Name', 'Heave Kinematics', 'NumberTitle', 'off');

% 子图1：升沉速度
subplot(2, 1, 1); % 2行1列，第1个图
plot(t_hist, v_hist, 'b-', 'LineWidth', 1.5);
ylabel('垂直速度 (m/s)');
xlabel('时间 (s)');
title('升沉方向 - 速度变化');
grid on;

% 子图2：升沉加速度
subplot(2, 1, 2); % 2行1列，第2个图
plot(t_hist, acc_hist, 'r-', 'LineWidth', 1.5);
ylabel('垂直加速度 (m/s^2)');
xlabel('时间 (s)');
title('升沉方向 - 加速度变化');
grid on;



%% ==================== 牛顿迭代法 ====================
function [P_solved, converged] = solve_pressure(P_init, S_vec, N1, N2, dzdt, Area, max_it, tol)
    P = P_init;
    converged = false;
    lambda = 0.6; 
    for i = 1:max_it
        J = compute_jacobian(P, S_vec, N1, N2, dzdt, Area);
        F = residual_function(P, S_vec, N1, N2, dzdt, Area);
        if norm(F) < tol, converged = true; break; end
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
    
    Q_FAN1 = (-1280 * sqrt(abs(P(5) - 300)) * sign(P(5) - 300) - 31.6 * (P(5) - 300)) * N_FAN1 / 2000;
    Q_FAN2 = (-1280 * sqrt(abs(P(6) - 300)) * sign(P(6) - 300) - 31.6 * (P(6) - 300)) * N_FAN2 / 2000;

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

    Q1 = -S(1) * 14.5 * sqrt(abs(P(1))) * sign(P(1));
    Q2 = -S(2) * 14.5 * sqrt(abs(P(2))) * sign(P(2));
    Q3 = -S(3) * 14.5 * sqrt(abs(P(3))) * sign(P(3));
    Q4 = -S(4) * 14.5 * sqrt(abs(P(4))) * sign(P(4));


    Q_PUMP = Area * dzdt;

    F = zeros(6, 1);
    
    F(5) = -Q_INC1 - Q_INC2 + Q_NOZ1 + Q_FAN1;
    F(6) = -Q_INC3 - Q_INC4 + Q_NOZ2 + Q_FAN2;
    F(1) = Q_INC1 + Q_IC1 - Q_IC2 + Q1 - Q_PUMP;
    F(2) = Q_INC2 + Q_IC2 - Q_IC3 + Q2 - Q_PUMP;
    F(3) = Q_INC3 + Q_IC3 - Q_IC4 + Q3 - Q_PUMP;
    F(4) = Q_INC4 + Q_IC4 - Q_IC1 + Q4 - Q_PUMP;
end