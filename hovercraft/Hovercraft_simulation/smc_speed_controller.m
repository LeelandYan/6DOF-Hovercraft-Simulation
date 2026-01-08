function [X_prop_req, cmd_rpm] = smc_speed_controller(state, target_u, dot_target_u)
% 基础滑模控制器
%
% 输入参数:
%   state        : [12x1] 状态向量
%   target_u     : 期望纵向速度 (m/s)
%   dot_target_u : 期望纵向加速度 (m/s^2)
%
% 输出参数:
%   X_prop_req   : 所需纵向推力 (N)
%   cmd_rpm      : 对应的螺旋桨转速 (RPM)

    %% 参数定义
    % 物理常数
    FT2M = 0.3048;
    SLUG2KG = 14.5939;
    LBF2N = 4.44822;
    
    % 质量参数
    m_slugs = 10879.5;
    m_kg = m_slugs * SLUG2KG;
    
    % 螺旋桨参数
    propeller_angle = 15; 
    
    %% 滑模控制器增益
    c_k_s = 500;    
    c_eta = 2000;   
    k_surface = 1.0; 

    %% 状态提取
    % phi   = state(4);
    % theta = state(5);
    % psi   = state(6); 
    u     = state(7);
    v     = state(8);
    w     = state(9);
    q     = state(11);
    r     = state(12);

    %% 动力学耦合项计算
    F_coupling = m_kg * (-q*w + r*v);

    %% 滑模控制律计算
    e = u - target_u;
    s = k_surface * e;
    sign_s = sign(s);
  
    
    term_feedforward = m_kg * dot_target_u;        % m * dot_u_d
    term_coupling    = -F_coupling;                % -mvr + mwq 
    term_linear      = -m_kg * c_k_s * e;          % -m * k_s * (u-u_d)
    term_switching   = -(m_kg * c_eta / k_surface) * sign_s; % - (m*eta/k) * sgn(s)
    
    % 总推力请求
    X_prop_req = term_feedforward + term_coupling + term_linear + term_switching;

    %% 推力物理限制处理
    if X_prop_req < 0
        X_prop_req = 0; % 螺旋桨不能产生反推力
    end
    
    % 设定一个物理上限 (参考 ESO 代码中的 200000)
    MAX_THRUST = 200000;
    if X_prop_req > MAX_THRUST
        X_prop_req = MAX_THRUST;
    end

    %% RPM 反解
    
    % 计算单个螺旋桨所需的推力
    Thrust_one_req_N = X_prop_req / 2;
    Thrust_one_req_lbf = Thrust_one_req_N / LBF2N;
    
    true_wind_speed = 0;
    true_wind_dir = 0;
    XWIND_ft_s = u / FT2M;
    
    % 推力系数 K_thrust 计算
    Thrust_base = (338 * propeller_angle + 4.36 * propeller_angle^2);
    Thrust_loss = (1.43 * propeller_angle * XWIND_ft_s + 0.1715 * XWIND_ft_s^2);
    
    K_thrust = Thrust_base - Thrust_loss;
    if K_thrust < 100, K_thrust = 100; end
    
    % 反解 RPM
    rpm_sq_ratio = Thrust_one_req_lbf / K_thrust;
    
    if rpm_sq_ratio < 0
        cmd_rpm = 0;
    else
        cmd_rpm = 1250 * sqrt(rpm_sq_ratio);
    end
    
    % 限幅
    MAX_RPM = 2500;
    MIN_RPM = 900;
    cmd_rpm = max(min(cmd_rpm, MAX_RPM), MIN_RPM);

end