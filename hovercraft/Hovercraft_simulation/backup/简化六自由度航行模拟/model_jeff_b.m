function dXdt = model_jeff_b(t, X)
    %%气垫船数学模型
    
    %% --- 单位转换系数 ---
    FT2M = 0.3048;          % feet to meters
    SLUG2KG = 14.5939;      % slugs to kg
    LBF2N = 4.44822;        % pound-force to Newtons
    rho_air_SI = 1.225;     % kg/m^3
    g_SI = 9.80665;         % m/s^2
    
    % 辅助转换
    M2FT = 1/FT2M;
    N2LBF = 1/LBF2N;

    %% --- 定义状态量 ---
    phi   = X(4); % Roll (横摇)
    theta = X(5); % Pitch (纵摇)
    psi   = X(6); % Yaw (艏向)
    u = X(7);  % 纵向速度
    v = X(8);  % 横向速度
    w = X(9);  % 垂向速度
    p = X(10); % 横摇角速度
    q = X(11); % 纵摇角速度
    r = X(12); % 艏摇角速度
    
    %% --- 气垫船参数 ---
    m_slugs = 10879.5;
    Ixx_imp = 5.672e6;  
    Iyy_imp = 1.629e7;
    Izz_imp = 2.057e7;
    
    % 转换为国际单位制
    m_kg = m_slugs * SLUG2KG;
    Ixx = Ixx_imp * SLUG2KG * (FT2M^2);
    Iyy = Iyy_imp * SLUG2KG * (FT2M^2);
    Izz = Izz_imp * SLUG2KG * (FT2M^2);
    
    % 空气阻力作用点
    X1SC = 30.0 * FT2M;  % 纵向偏移 (x)
    X2SC = 0 * FT2M;  % 横向偏移 (y)
    X3SC = -8.0 * FT2M;    % 垂向偏移 (z)

    
    % 几何参数
    X1R_si = -67.1 * FT2M;      % 舵 X 位置
    X3R_si = 0.0 * FT2M;        % 舵 Z 位置
    RSAREA_si = 47.2 * (FT2M^2);  
    SCAREA_si = 3200 * (FT2M^2);  
    FCAREA_si = 836 * (FT2M^2);   
    DUCT_AREA_si = 123 * (FT2M^2);
    
    %% --- 3. 控制输入 ---
    if t < 5
        propeller_angle = 15 * (t/5);
    else
        propeller_angle = 25;
    end

    propeller_speed_rpm = 1250;
%     propeller_speed_rpm = 0;
    
    if t > 30
        rudder_angle_deg = 5;
    else
        rudder_angle_deg = 0;
    end

    rudder_angle_rad = deg2rad(rudder_angle_deg);
    
    ture_wind_speed_si = 0;       % 真风速
    true_wind_direction_si = deg2rad(180);   % 真风来向
    

%% --- 4. 空气阻力计算  ---    
    % --- 相对风计算(Apparent Wind) Eq.60-63 ---
    % 1. 将真风速度分解到地球坐标系 (NED)
    V_wind_north = ture_wind_speed_si * cos(true_wind_direction_si + pi);
    V_wind_east  = ture_wind_speed_si * sin(true_wind_direction_si + pi);
    
    % 2. 将真风速度分量转换到船体坐标系
    u_wind_body = V_wind_north * cos(psi) + V_wind_east * sin(psi);
    v_wind_body = -V_wind_north * sin(psi) + V_wind_east * cos(psi);
    
    % 3. 计算相对风速
    % 相对风 = 这风分量 - 船速
    u_rel = u_wind_body - u;
    v_rel = v_wind_body - v;
    
    % 相对风速
    apparent_wind_velocity_si = sqrt(u_rel^2 + v_rel^2) + 0.001; % 防止除零
    
    % 计算相对风角 (Apparent Wind Angle, Beta)
    apparent_wind_angle_rad = atan2(-v_rel, -u_rel); 
    apparent_wind_angle_deg = rad2deg(apparent_wind_angle_rad);
    abs_beta  = abs(apparent_wind_angle_deg);
    
    % --- 气动系数计算 Eq.70 ---
    % 侧向力系数 Eq.70 
    if abs_beta <= 180
        SDRAG = 0.01385 * abs_beta - 7.69e-5 * abs_beta^2;
    else
        SDRAG = 0;
    end
   
    if apparent_wind_angle_deg < 0, SDRAG = -SDRAG; end
    
    % 正面阻力系数 Eq.71 
    if abs_beta < 40
        FDRAG = -2.22e-4 * abs_beta^2 + 3.33e-3 * abs_beta + 0.5;
    elseif abs_beta < 88
        FDRAG = -2.48e-5 * abs_beta^3 + 5e-3 * abs_beta^2 ...
                      - 0.324 * abs_beta + 5.835;
    else
        FDRAG = 1.04e-4 * abs_beta^2 - 3.518e-2 * abs_beta + 2.39;
    end
    
    % 艏摇力矩系数 Eq.72 
    if abs_beta < 60
        YDRAG = 1.67e-3 * abs_beta;
    elseif abs_beta <= 120
        YDRAG = 5.17e-7 * abs_beta^3 - 1.23e-4 * abs_beta^2 ...
                      + 5.82e-3 * abs_beta + 8.27e-2;
    else
        YDRAG = 1.67e-3 ; % 保持常数
    end
%     if apparent_wind_angle_deg < 0, YDRAG = -YDRAG; end

    % --- 4.3 空气阻力与力矩计算 Eq.73-77 ---
    q_bar = 0.5 * rho_air_SI * apparent_wind_velocity_si^2; % 动压
    LCUSH_si = 77.0 * FT2M;             % 气垫长度 

    % 空气阻力计算
    % Eq.73: XBDRAG (纵向空气阻力) 
    XBDRAG_si = -FDRAG * FCAREA_si * q_bar;
    
    % Eq.74: YBDRAG (横向空气阻力) 
    YBDRAG_si = -SDRAG * SCAREA_si * q_bar; 

    % 力矩 (Moments)
    % Eq.75: Pitch Moment
    PITCHA_si = XBDRAG_si * X3SC;  
    
    % Eq.76: Roll Moment 
    ROLLA_si  = -YBDRAG_si * X3SC;
    
    % Eq.77: Yaw Moment 
%     YAWBD_si = (YDRAG * SCAREA_si * LCUSH_si * q_bar) ...
%              + (YBDRAG_si * X1SC);
     YAWBD_si = YBDRAG_si * X1SC;
    

    %% --- 5. 推进与操纵力计算 ---
%     Thrust_one_lbf = (338 * propeller_angle + 4.36 * propeller_angle^2) * (propeller_speed_rpm/1250)^2;
%     Thrust_one_si = Thrust_one_lbf * LBF2N;
%     
%     propeller_thrust_si = 2 * Thrust_one_si;
 
    XWIND_ft_s = -u_rel * M2FT; 
    XWIND_val = XWIND_ft_s; 
    Thrust_base = (338 * propeller_angle + 4.36 * propeller_angle^2);
    Thrust_loss = (1.43 * propeller_angle * XWIND_val + 0.1715 * XWIND_val^2);
    Thrust_one_lbf = (Thrust_base - Thrust_loss) * (propeller_speed_rpm/1250)^2;
    
    % 转换为牛顿
    Thrust_one_si = Thrust_one_lbf * LBF2N;
    
    % 总推力
    propeller_thrust_si = 2 * Thrust_one_si;

    
    % 舵力
    VDUCT2_si = (apparent_wind_velocity_si * cos(apparent_wind_angle_rad))^2 + Thrust_one_si / (rho_air_SI * DUCT_AREA_si);
    if VDUCT2_si < 0, VDUCT2_si = 0; end
    PDUCT_si = 0.5 * rho_air_SI * VDUCT2_si;
    
    CLIFT_R = 0.053 * rudder_angle_deg;
    CDRAG_R = 0.422e-3 * rudder_angle_deg^2;
    
    RDRAG_si = -CDRAG_R * RSAREA_si * 2 * PDUCT_si;
    RLIFT_si = CLIFT_R * RSAREA_si * 2 * PDUCT_si;
    
    % 舵力矩
    % Mz(Yaw) = F_y * x
    RYAW_si   = X1R_si * RLIFT_si;
    % Mx(Roll) = F_y * z
    RROLL_si  = -X3R_si * RLIFT_si;
    % My(Pitch) = F_x * z
    RPITCH_si = X3R_si * RDRAG_si;
    
    %% --- 6. 气垫与水动力 ---
    F_cushion_static_si = m_kg * g_SI;
    
    K_heave_si = 20000 * LBF2N * M2FT; 
    C_heave_si = 5000  * LBF2N * M2FT;
    K_pitch_si = 4.0e7 * LBF2N * FT2M; 
    C_pitch_si = 5.0e6 * LBF2N * FT2M;
    K_roll_si  = 1.5e7 * LBF2N * FT2M;
    C_roll_si  = 2.0e6 * LBF2N * FT2M;
    
    % 气垫力与力矩
    FX3CT_si = -F_cushion_static_si - (K_heave_si * X(3) + C_heave_si * w);
    FX4CT_si = -(K_roll_si * phi + C_roll_si * p);
    FX5CT_si = -(K_pitch_si * theta + C_pitch_si * q);
    
    % 水动阻力
    CD_skirt = 0.25;
    Area_wet_surge = 50 * (FT2M^2);
    Area_wet_sway  = 80 * (FT2M^2);
    Rho_water = 1025;
    
    SDRX_si = -0.5 * Rho_water * CD_skirt * Area_wet_surge * u * abs(u);
    SDRY_si = -0.5 * Rho_water * CD_skirt * Area_wet_sway  * v * abs(v);
    
    YAWDC_coeff_si = 2.77e6 * LBF2N * FT2M;
    YAWDC_si = -YAWDC_coeff_si * r;
    
    %% --- 7. 总力与力矩汇总 (核心保留部分) ---
    % 注意：在 NED 坐标系中，Z轴向下。
    % 重力分解：
    % X轴重力分量 ~ -mg*theta
    % Y轴重力分量 ~ mg*phi
    % Z轴重力分量 ~ mg
    XGRAV_si = -m_kg * g_SI * theta;
    YGRAV_si =  m_kg * g_SI * phi;
    ZGRAV_si =  m_kg * g_SI;
    
    % 按照你的要求保留的六行合力公式：
    Fx = propeller_thrust_si + RDRAG_si + XBDRAG_si + SDRX_si + XGRAV_si;
    Fy = RLIFT_si + YBDRAG_si + SDRY_si + YGRAV_si;
    Fz = FX3CT_si + ZGRAV_si;
    
    Mx = FX4CT_si + ROLLA_si + RROLL_si;
    My = FX5CT_si + PITCHA_si + RPITCH_si;
    Mz = YAWBD_si + RYAW_si + YAWDC_si;
    
    %% --- 8. 动力学方程 ---    
    udot = Fx/m_kg - q*w + r*v;
    vdot = Fy/m_kg - r*u + p*w;
    wdot = Fz/m_kg - p*v + q*u;
    
    pdot = (Mx - (Izz - Iyy)*q*r) / Ixx;
    qdot = (My - (Ixx - Izz)*r*p) / Iyy;
    rdot = (Mz - (Iyy - Ixx)*p*q) / Izz;
    
    %% --- 9. 运动学方程 ---
    dx = u * cos(psi) - v * sin(psi);
    dy = u * sin(psi) + v * cos(psi);
    dz = w;
    
    dphi   = p;
    dtheta = q;
    dpsi   = r;
    
    %% --- 输出 ---
    dXdt = [dx; dy; dz; dphi; dtheta; dpsi; ...
            udot; vdot; wdot; pdot; qdot; rdot];
end