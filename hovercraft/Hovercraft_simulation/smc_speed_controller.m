function [cmd_rpm, debug_info] = smc_speed_controller(state, target_u, dot_target_u, wind_info)
% 气垫船航速滑模控制器
%
% 输入参数:
%   state        : [12x1] 状态向量 [x, y, z, phi, theta, psi, u, v, w, p, q, r]
%   target_u     : 期望纵向速度 (m/s)
%   dot_target_u : 期望纵向加速度 (m/s^2)
%   wind_info    : 结构体，包含 true_wind_speed_si 和 true_wind_direction_si
%
% 输出参数:
%   cmd_rpm      : 螺旋桨转速指令 (RPM)
%   debug_info   : 调试信息结构体 (包含误差、S值、所需推力等)

    %% 参数定义
    % 单位转换
    FT2M = 0.3048;
    SLUG2KG = 14.5939;
    LBF2N = 4.44822;
    rho_air = 1.225;
    g = 9.80665;
    
    % 质量与惯量
    m_slugs = 10879.5;
    m_kg = m_slugs * SLUG2KG;
    
    % 估算阻力
    FCAREA_si = 836 * (FT2M^2);   % 正面迎风面积
    SCAREA_si = 3200 * (FT2M^2);  % 侧面迎风面积
    propeller_angle = 15;         % 螺旋桨安装角
    
    %%
    % 滑模控制器参数
    c_epsilon = 100;   % 趋近律参数 epsilon 
    c_k       = 500;    % 趋近律参数 k 
    delta     = 2;    % 饱和函数边界层厚度 
    
    %% 2. 状态提取
    phi   = state(4);
    theta = state(5);
    psi   = state(6);
    u     = state(7);
    v     = state(8);
    w     = state(9);
    q     = state(11);
    r     = state(12);

    %% 计算扰动项
    
    % 重力分量X_G 
    X_G = -m_kg * g * theta; 
    
    % 围裙阻力
    Area_wet_surge = 50 * (FT2M^2);
    Rho_water = 1025;
    X_skirt = -0.5 * Rho_water * 0.25 * Area_wet_surge * u * abs(u);
    
    % 空气阻力X_air 
    if nargin < 4 || isempty(wind_info)
        true_wind_speed = 0;
        true_wind_dir = 0;
    else
        true_wind_speed = wind_info.speed;
        true_wind_dir = wind_info.dir;
    end
    
    u_rel = u - true_wind_speed * cos(true_wind_dir - psi); 
    v_rel = v - true_wind_speed * sin(true_wind_dir - psi);
    v_app_sq = u_rel^2 + v_rel^2;
    
    CD_air_surge = 0.5; 
    q_bar = 0.5 * rho_air * v_app_sq;
    X_air = -CD_air_surge * FCAREA_si * q_bar * sign(u_rel);
    

    F_coupling = m_kg * (-q*w + r*v);

    % 总的干扰
    F_rest = X_G + X_skirt + X_air; 

    %% 滑模控制律计算
    % 误差
    e = u - target_u;
    
    % 滑模面 s
    s = e; 
    
    % 饱和函数sat(s/delta) 
    if abs(s/delta) <= 1
        sat_s = s/delta;
    else
        sat_s = sign(s);
    end
    
    X_prop_req = m_kg * (dot_target_u - F_coupling/m_kg + (-c_epsilon * sat_s - c_k * s)) - F_rest;
    
    %% 推力分配与转速反解
    % 计算单个螺旋桨所需的推力 
    if X_prop_req < 0
        X_prop_req = 0; % 螺旋桨不能产生反推力
    end
    
    Thrust_one_req_N = X_prop_req / 2;
    Thrust_one_req_lbf = Thrust_one_req_N / LBF2N;
    
    % 计算推力系数
    XWIND_ft_s = u / FT2M; % 简化的相对风速
    
    Thrust_base = (338 * propeller_angle + 4.36 * propeller_angle^2);
    Thrust_loss = (1.43 * propeller_angle * XWIND_ft_s + 0.1715 * XWIND_ft_s^2);
    
    K_thrust = Thrust_base - Thrust_loss;
    if K_thrust < 100, K_thrust = 100; end % 防止除以零或负数
    
    % 反解 RPM
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