function [X_prop_req, debug_info] = smc_speed_controller_adaptive(state, target_u, dot_target_u, wind_info)
    %% 参数定义
    % --- 单位转换 ---
    FT2M = 0.3048;
    SLUG2KG = 14.5939;
    rho_air = 1.225;
    g = 9.80665;
    
    % --- 质量与几何参数 ---
    m_slugs = 10879.5;
    m_kg = m_slugs * SLUG2KG;
    
    FCAREA_si = 836 * (FT2M^2);   % 正面迎风面积
    
    % 初始控制器参数
    c_epsilon = 5;     % 切换增益
    c_k       = 10;    % 比例增益
    delta     = 2.0;   % 边界层厚度 
    
    % 自适应学习率
    gamma     = 500;   
    dt_local  = 0.05; 

    %% 状态提取
    theta = state(5);
    psi   = state(6);
    u     = state(7);
    v     = state(8);
    w     = state(9);
    q     = state(11);
    r     = state(12);

    %% 模型力计算 
    % 重力分量
    X_G = -m_kg * g * theta; 
    
    % 围裙阻力 
    Area_wet_surge = 50 * (FT2M^2);
    Rho_water = 1025;
    X_skirt = -0.5 * Rho_water * 0.25 * Area_wet_surge * u * abs(u);
    
    % 空气阻力
    if nargin < 4 || isempty(wind_info)
        true_wind_speed = 0;
        true_wind_dir = 0;
    else
        true_wind_speed = wind_info.speed;
        true_wind_dir = wind_info.dir;
    end
    
    u_rel = u - true_wind_speed * cos(true_wind_dir - psi); 
    CD_air_surge = 0.5; 
    v_app_sq = u_rel^2; % 简化计算，主要看纵向风
    q_bar = 0.5 * rho_air * v_app_sq;
    X_air = -CD_air_surge * FCAREA_si * q_bar * sign(u_rel);
    
    F_coupling = m_kg * (-q*w + r*v);

    % 可计算的已知扰动总和
    F_nominal = X_G + X_skirt + X_air; 

    %% 自适应扰动估计
    persistent F_hat 
    if isempty(F_hat)
        F_hat = 0; 
    end
    
    % 误差定义
    e = u - target_u;
    s = e; % 滑模面
    
    % 自适应更新律 
    d_F_hat = gamma * s; 
    F_hat = F_hat + d_F_hat * dt_local;
    F_hat = max(min(F_hat, 50000), -50000); % 限幅

    %% 控制律计算
    if abs(s/delta) <= 1
        sat_s = s/delta;
    else
        sat_s = sign(s);
    end
    
    F_robust = m_kg * (c_epsilon * sat_s + c_k * s);
    
    % 计算总推力需求 (Total Thrust Requirement)
    X_prop_req = m_kg * (dot_target_u - F_coupling/m_kg) ...
                 - (F_nominal + F_hat) ...  
                 - F_robust;                
                 
    %% 输出限幅
    if X_prop_req < 0
        X_prop_req = 0; 
    end
    
    MAX_THRUST = 200000;
    if X_prop_req > MAX_THRUST
        X_prop_req = MAX_THRUST;
    end

    %% 调试信息
    debug_info.e = e;
    debug_info.s = s;
    debug_info.F_nominal = F_nominal;
    debug_info.F_hat = F_hat;        
    debug_info.X_prop_req = X_prop_req; % 记录推力需求

end