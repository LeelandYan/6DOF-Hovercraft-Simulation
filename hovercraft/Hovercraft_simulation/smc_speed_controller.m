function [cmd_rpm, debug_info] = smc_speed_controller(state, target_u, dot_target_u, wind_info)
% 气垫船航速滑模控制器
%
% 输入参数:
%   state        : [12x1] 状态向量 [x, y, z, phi, theta, psi, u, v, w, p, q, r]
%   target_u     : 期望纵向速度 (m/s)
%   dot_target_u : 期望纵向加速度 (m/s^2)
%   wind_info    : (可选) 结构体，包含 true_wind_speed_si 和 true_wind_direction_si
%
% 输出参数:
%   cmd_rpm      : 螺旋桨转速指令 (RPM)
%   debug_info   : 调试信息结构体 (包含误差、S值、所需推力等)

    %% 1. 参数定义 (需与 model_jeff_b.m 保持一致)
    % 单位转换
    FT2M = 0.3048;
    SLUG2KG = 14.5939;
    LBF2N = 4.44822;
    rho_air = 1.225;
    g = 9.80665;
    
    % 质量与惯量
    m_slugs = 10879.5;
    m_kg = m_slugs * SLUG2KG;
    
    % 几何参数 (用于估算阻力)
    FCAREA_si = 836 * (FT2M^2);   % 正面迎风面积
    SCAREA_si = 3200 * (FT2M^2);  % 侧面迎风面积
    propeller_angle = 15;         % 螺旋桨安装角 (度)
    
    % 滑模控制器参数 (需调试)
    c_epsilon = 1000;   % 趋近律参数 epsilon [cite: 145]
    c_k       = 500;    % 趋近律参数 k [cite: 145]
    delta     = 0.5;    % 饱和函数边界层厚度 [cite: 156]
    
    %% 2. 状态提取
    phi   = state(4);
    theta = state(5);
    psi   = state(6);
    u     = state(7);
    v     = state(8);
    w     = state(9);
    q     = state(11);
    r     = state(12);

    %% 3. 计算扰动项/模型已知项 (F_drag + F_gravity + F_coriolis)
    % 目的是计算出为了维持当前运动，除了加速还需要克服多少“其它力”
    
    % --- A. 重力分量 X_G ---
    X_G = -m_kg * g * theta; % 近似 sin(theta) ~ theta
    
    % --- B. 围裙水阻力 X_skirt (估算) ---
    % 依据 model_jeff_b.m 中的公式
    Area_wet_surge = 50 * (FT2M^2);
    Rho_water = 1025;
    X_skirt = -0.5 * Rho_water * 0.25 * Area_wet_surge * u * abs(u);
    
    % --- C. 空气阻力 X_air (估算) ---
    % 如果没有风速信息，默认无风
    if nargin < 4 || isempty(wind_info)
        true_wind_speed = 0;
        true_wind_dir = 0;
    else
        true_wind_speed = wind_info.speed;
        true_wind_dir = wind_info.dir;
    end
    
    % 计算相对风 (简化计算，用于控制器前馈补偿)
    % 仅考虑纵向主导分量，完整计算参考 model_jeff_b
    u_rel = u - true_wind_speed * cos(true_wind_dir - psi); 
    v_rel = v - true_wind_speed * sin(true_wind_dir - psi); % 简化的相对风
    v_app_sq = u_rel^2 + v_rel^2;
    
    % 估算 FDRAG 系数 (取平均值或简化模型，这里取0.5作为近似)
    % 实际应根据 beta 角查表，此处为控制器简化
    CD_air_surge = 0.5; 
    q_bar = 0.5 * rho_air * v_app_sq;
    X_air = -CD_air_surge * FCAREA_si * q_bar * sign(u_rel);
    
    % --- D. 运动耦合项 ---
    % 动力学方程: m(du/dt + qw - rv) = Fx
    % 移项得 du/dt = Fx/m - qw + rv
    F_coupling = m_kg * (-q*w + r*v);

    % 总的非推力项 (视为干扰 D)
    F_rest = X_G + X_skirt + X_air; 

    %% 4. 滑模控制律计算
    % 误差
    e = u - target_u;
    
    % 滑模面 s
    s = e; % 
    
    % 饱和函数 sat(s/delta) 
    if abs(s/delta) <= 1
        sat_s = s/delta;
    else
        sat_s = sign(s);
    end
    
    % 计算所需总纵向合力 Fx_req
    % 理想动力学: m * du/dt = Fx - m(qw - rv)
    % 期望动力学: du/dt = dot_target_u - epsilon*sat(s) - k*s
    % 联立求解所需的 Fx_propeller
    
    % 滑模控制输出推力 (牛顿)
    % F_prop_req = m * (dot_u_d - (-epsilon*sat(s) - k*s) + (qw - rv)) - F_rest
    % 注意符号：趋近律是让 e 趋向于 0，即 du 趋向于 du_d
    
    %  修改版：考虑到 F_rest 包含在系统模型中
    X_prop_req = m_kg * (dot_target_u - F_coupling/m_kg + (-c_epsilon * sat_s - c_k * s)) - F_rest;
    
    %% 5. 推力分配与转速反解
    % model_jeff_b.m 中的推力模型:
    % Thrust_base = (338 * angle + 4.36 * angle^2);
    % Thrust_loss = ... (与航速有关)
    % Thrust_one_lbf = (Base - Loss) * (RPM/1250)^2
    % X_prop_total_N = 2 * Thrust_one_lbf * LBF2N
    
    % 5.1 计算单个螺旋桨所需的推力 (Lbf)
    if X_prop_req < 0
        X_prop_req = 0; % 螺旋桨不能产生反推力(除非变距，这里假设不能)
    end
    
    Thrust_one_req_N = X_prop_req / 2;
    Thrust_one_req_lbf = Thrust_one_req_N / LBF2N;
    
    % 5.2 计算推力系数 (Base - Loss)
    XWIND_ft_s = u / FT2M; % 简化的相对风速
    
    Thrust_base = (338 * propeller_angle + 4.36 * propeller_angle^2);
    Thrust_loss = (1.43 * propeller_angle * XWIND_ft_s + 0.1715 * XWIND_ft_s^2);
    
    K_thrust = Thrust_base - Thrust_loss;
    if K_thrust < 100, K_thrust = 100; end % 防止除以零或负数
    
    % 5.3 反解 RPM
    % Thrust = K * (RPM/1250)^2  =>  RPM = 1250 * sqrt(Thrust / K)
    rpm_sq_ratio = Thrust_one_req_lbf / K_thrust;
    
    if rpm_sq_ratio < 0
        cmd_rpm = 0;
    else
        cmd_rpm = 1250 * sqrt(rpm_sq_ratio);
    end
    
    %% 6. 限幅与输出
    MAX_RPM = 2500;
    MIN_RPM = 900; % 怠速
    
    cmd_rpm = max(min(cmd_rpm, MAX_RPM), MIN_RPM);
    
    %% 7. 调试信息
    debug_info.e = e;
    debug_info.s = s;
    debug_info.X_prop_req = X_prop_req;
    debug_info.F_rest = F_rest;

end