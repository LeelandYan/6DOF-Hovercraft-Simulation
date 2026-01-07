function X_prop_req = smc_speed_controller_eso(state, target_u, dot_target_u)
    %% 参数定义
    SLUG2KG = 14.5939;
    
    % 质量参数 
    m_slugs = 10879.5;
    m_kg = m_slugs * SLUG2KG;
    b0 = 1 / m_kg; 

    %%  ESO参数设置
    % 观测器带宽 omega_o: 决定了ESO估计扰动的快慢。
    omega_o = 3.0; 
    
    % 根据带宽法配置 ESO 增益 
    beta1 = 2 * omega_o;      
    beta2 = omega_o * omega_o;
    
    dt_local = 0.05; % 控制周期，与仿真步长一致

    %% 状态提取
     u_sensor = state(7); % 真实测量的速度
    
    %% 状态观测器
    persistent z1 z2 last_thrust
    
    % 初始化 
    if isempty(z1)
        z1 =  u_sensor; % 初始速度估计,设为当前测量值
        z2 = 0;      % 初始扰动设为0
        last_thrust = 0; % 上一时刻施加的推力
    end
    
    % 计算观测误差
    % 观测误差 = 传感器测量值 - 观测器估计值
    e_obs =  u_sensor - z1;
    
    % 更新观测器状态
    % z1_dot = z2 + b0 * u + beta1 * e
    % z2_dot = beta2 * e
    
    z1_next = z1 + (z2 + b0 * last_thrust + beta1 * e_obs) * dt_local;
    z2_next = z2 + (beta2 * e_obs) * dt_local;
    
    % 更新状态
    z1 = z1_next;
    z2 = z2_next;
    
    % 将扰动加速度转换为力 
    F_disturbance_est = z2 * m_kg; 

    %% 滑模控制律计算
    % 跟踪误差 
    e_track =  u_sensor - target_u;
    s = e_track; 
    
    % 滑模参数
    c_k = 10;      % 比例增益 
    c_epsilon = 5; % 切换增益 
    delta = 2.0;   % 边界层厚度
    
    % 饱和函数
    if abs(s/delta) <= 1
        sat_s = s/delta;
    else
        sat_s = sign(s);
    end
    
    % 鲁棒反馈项
    F_robust = m_kg * (c_epsilon * sat_s + c_k * s);
    
    % 最终控制律
    % 推力 = (期望加速度 - 估计的扰动加速度) / b0 - 反馈项
    % X_prop = m * dot_target_u - F_disturbance - F_robust
    
    X_prop_req = m_kg * dot_target_u - F_disturbance_est - F_robust;


    if X_prop_req < 0
        X_prop_req = 0; % 螺旋桨不能反推
    end
    
    MAX_THRUST = 200000;
    if X_prop_req > MAX_THRUST
        X_prop_req = MAX_THRUST;
    end
    

    last_thrust = X_prop_req; 



end