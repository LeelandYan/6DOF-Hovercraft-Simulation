%% 气垫船横摇自由度仿真
clear; clc; close all;

%% 1. 参数定义 (基于 model_jeff_b 和 Liu et al. 2023)
% --- 单位转换 ---
FT2M = 0.3048;          % 英尺转米
SLUG2KG = 14.5939;      % slug转kg
g_SI = 9.80665;         % 重力加速度 m/s^2

% --- 质量与惯量特性 ---
m_slugs = 10879.5;
Ixx_imp = 5.672e6;      % 横摇惯量 (slug-ft^2)

m_kg = m_slugs * SLUG2KG;
Ixx_si = Ixx_imp * SLUG2KG * (FT2M^2);

% --- 几何尺寸 ---
% 气垫尺寸 (用于计算力臂和面积)
L_cush_si = 38.5 * FT2M;    % 气垫长度
W_cush_si = 17.5 * FT2M;    % 气垫宽度
Area_Chamber = (L_cush_si/2) * (W_cush_si/2) * 2; % 单侧气室近似面积 (前后合并计算)

% 关键高度参数
H_base_m = 5.0 * FT2M;      % 硬结构(船底)离地基准高度
H_skirt_m = 4.5 * FT2M;     % 围裙高度 (当气垫完全充气时的设计高度)

% 阻尼系数 (用于模拟水/空气阻尼，防止无休止震荡)
C_heave = 20000;  % 垂向阻尼 N/(m/s)
C_roll  = 800000; % 横摇阻尼 Nm/(rad/s)

%% 2. 仿真设置
% --- 初始条件 ---
% [z(垂向位移), w(垂向速度), phi(横摇角), p(横摇角速度)]
Y0 = [0; 0; 0; 0]; 

% --- 控制输入：设置左右风机转速差 ---
% 设想：左侧转速低，右侧转速高 -> 船应该向左倾斜 (右侧抬高)
Control.RPM_Left  = 1400; 
Control.RPM_Right = 1600; 

% --- 时间设置 ---
t_span = [0 15]; % 模拟 15 秒

%% 3. 调用 ODE45 求解器
% 注意：我们将动力学方程定义在脚本底部的 local function 中
options = odeset('RelTol', 1e-6, 'AbsTol', 1e-6);
[T, Y] = ode45(@(t,y) system_dynamics(t, y, Control, m_kg, Ixx_si, ...
    L_cush_si, W_cush_si, H_base_m, H_skirt_m, Area_Chamber, C_heave, C_roll, g_SI), ...
    t_span, Y0, options);

%% 4. 数据后处理与绘图
% ODE45 只返回状态量，我们需要重新计算过程中的压力和高度以便绘图
num_steps = length(T);
P_Left_hist = zeros(num_steps, 1);
P_Right_hist = zeros(num_steps, 1);
Gap_Left_hist = zeros(num_steps, 1);
Gap_Right_hist = zeros(num_steps, 1);

for i = 1:num_steps
    % 调用辅助函数计算当前的物理状态
    [~, Info] = system_dynamics(T(i), Y(i,:)', Control, m_kg, Ixx_si, ...
    L_cush_si, W_cush_si, H_base_m, H_skirt_m, Area_Chamber, C_heave, C_roll, g_SI);
    
    P_Left_hist(i)  = Info.P_L;
    P_Right_hist(i) = Info.P_R;
    Gap_Left_hist(i) = Info.h_L;
    Gap_Right_hist(i) = Info.h_R;
end

% --- 绘图 ---
figure('Name', 'Hovercraft Roll Simulation', 'Color', 'w', 'Position', [100 100 800 600]);

% 1. 横摇角
subplot(2,1,1);
plot(T, rad2deg(Y(:,3)), 'LineWidth', 2, 'Color', 'b');
ylabel('横摇角 (deg)');
title(['横摇响应 (左风机=' num2str(Control.RPM_Left) ', 右风机=' num2str(Control.RPM_Right) ' RPM)']);
grid on;
% 添加参考线
yline(0, '--k');
text(T(end)*0.8, rad2deg(Y(end,3)), sprintf('稳态角度: %.2f°', rad2deg(Y(end,3))), 'FontSize', 10);

% 2. 气垫压力
subplot(2,1,2);
plot(T, P_Left_hist, 'r', 'LineWidth', 1.5); hold on;
plot(T, P_Right_hist, 'b', 'LineWidth', 1.5);
ylabel('气室压力 (Pa)');
legend('左侧气室', '右侧气室');
title('左右气室压力变化');
grid on;



%% ============================================================
%% 本地函数：系统动力学方程
%% ============================================================
function [dYdt, Info] = system_dynamics(~, Y, Ctrl, m, Ixx, L, W, H_base, H_skirt, Area, C_z, C_phi, g)
    % 解包状态
    z   = Y(1); % 垂向位移 (NED系，向下为正)
    w   = Y(2); % 垂向速度
    phi = Y(3); % 横摇角 (FRD系，右倾为正)
    p   = Y(4); % 横摇角速度

    %% A. 计算几何状态 (Gap Height)
    % 我们关注左右两侧中心点的气隙
    % 侧向力臂 (从中心线到气室中心的距离)
    y_arm_L = -W/4; % 左侧 (y < 0)
    y_arm_R =  W/4; % 右侧 (y > 0)
    
    % 计算物理离地高度 (Hard Structure Height from Water)
    % 公式：h_hard = H_base - z(沉降) - y * sin(phi)(横摇分量)
    h_hard_L = H_base - z - y_arm_L * sin(phi);
    h_hard_R = H_base - z - y_arm_R * sin(phi);
    
    % 计算气隙 (Gap) = 硬结构高度 - 围裙高度
    % 限制最小气隙为 1mm (0.001m)，模拟触地或密封状态，防止除零
    h_gap_L = max(h_hard_L - H_skirt, 0.001);
    h_gap_R = max(h_hard_R - H_skirt, 0.001);
    
    %% B. 气垫压力模型 (Simplified Fan Law + Leakage)
    % 核心逻辑：
    % 1. 风机提供的最大压头与转速平方成正比 (Fan Law)
    % 2. 实际气室压力取决于泄流程度 (Gap越大，压力越小)
    
    % 风机系数 (需要根据 4-73 风机特性拟合，这里取经验值)
    K_fan = 0.9e-3; 
    P_source_L = K_fan * Ctrl.RPM_Left^2;
    P_source_R = K_fan * Ctrl.RPM_Right^2;
    
    % 泄流衰减模型：P_cushion = P_source / (1 + C * h^2)
    % 这个模型模拟了伯努利方程的效果：开口越大，压力损失越剧烈
    Leakage_Factor = 400; 
    
    P_L = P_source_L / (1 + Leakage_Factor * h_gap_L^2);
    P_R = P_source_R / (1 + Leakage_Factor * h_gap_R^2);
    
    %% C. 力与力矩计算
    % 气动力 (向上为正)
    F_aero_L = P_L * Area;
    F_aero_R = P_R * Area;
    
    % 1. 垂向合力 (NED系: 向下为正)
    % F_z_total = 重力 - 气动升力 - 阻尼
    F_lift_total = F_aero_L + F_aero_R;
    F_z_net = m * g - F_lift_total - C_z * w;
    
    % 2. 横摇力矩 (FRD系: 右倾为正)
    % 左侧力(F_L)作用在 y_arm_L (负值)。向上力在NED中为负。
    % 力矩 M = r x F
    % 右侧 (y>0): 力向上(-F_R)。力矩 = y * (-F_R) -> 负力矩 (向左滚)
    % 左侧 (y<0): 力向上(-F_L)。力矩 = y * (-F_L) -> 正力矩 (向右滚)
    % 物理直觉：左侧推力大 -> 向右滚；右侧推力大 -> 向左滚。
    
    M_x_aero = (y_arm_L * -F_aero_L) + (y_arm_R * -F_aero_R);
    
    % 加上阻尼力矩
    M_x_net = M_x_aero - C_phi * p;
    
    %% D. 状态导数
    dzdt = w;
    dwdt = F_z_net / m;
    
    dphidt = p;
    dpdt = M_x_net / Ixx;
    
    dYdt = [dzdt; dwdt; dphidt; dpdt];
    
    % 保存辅助信息供后处理
    Info.P_L = P_L;
    Info.P_R = P_R;
    Info.h_L = h_gap_L;
    Info.h_R = h_gap_R;
end