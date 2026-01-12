function [Mz_cmd, debug_data] = smc_heading_controller(state, target_psi)
    % 状态提取
    psi = state(6); % 当前艏向 (rad)
    r   = state(12); % 当前艏摇角速度 (rad/s)
    
    % 误差计算 
    e_psi = target_psi - psi;
    while e_psi > pi,  e_psi = e_psi - 2*pi; end
    while e_psi < -pi, e_psi = e_psi + 2*pi; end
    
    % 滑模面定义
    lambda = 2.0; 
    s = r + lambda * (-e_psi); 
    
    % 控制律

    
    Izz = 2.057e7 * 14.5939 * (0.3048^2); % 估算的转动惯量
    
    k = 10.0;  % 比例增益
    yita = 0.5;
    epsilon = 0.5; % 边界层厚度
    
     
    if abs(s) < epsilon
        sat_s = s / epsilon;
    else
        sat_s = sign(s);
    end
    
    % 计算期望力矩 
    Mz_cmd = -Izz * ( k * s + yita * sat_s ); 
    
    % 限幅
    MAX_MOMENT = 60000; 
    Mz_cmd = max(min(Mz_cmd, MAX_MOMENT), -MAX_MOMENT);

    debug_data.e_psi = e_psi;
    debug_data.Mz_cmd = Mz_cmd;
end