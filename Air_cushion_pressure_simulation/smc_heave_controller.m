function [cmd_rpm, debug_info] = smc_heave_controller(state, target_z, dot_target_z, params)
% SMC_HEAVE_CONTROLLER 气垫船垫升高度滑模控制器
%
% 输入参数:
%   state        : [2x1] 向量 [当前高度 z (m), 当前垂直速度 v (m/s)]
%   target_z     : 期望高度 (m)
%   dot_target_z : 期望垂直速度 (m/s)，定高悬停时给 0
%   params       : 结构体，包含物理参数 (m, Area, L_skirt 等)
%
% 输出参数:
%   cmd_rpm      : 风机转速指令 (RPM)
%   debug_info   : 调试数据

    %% 1. 提取参数
    z = state(1);
    v = state(2);
    
    % 物理参数 (使用 SI 单位)
    m = params.mass_kg;     
    g = 9.81;
    
    % 几何参数 (总面积和总周长，用于集总模型计算)
    % 原代码中是4个气室，控制器将其视为一个大的等效气室进行计算
    A_total = params.Area_cushion_m2 * 4; 
    L_total = params.L_per_cushion_m * 4;
    
    % 围裙当前长度 (如果没有实时测量，使用平衡长度近似或传入估计值)
    % 这里假设params里带了当前的平均围裙长度，或者使用固定值
    if isfield(params, 'current_skirt_len_m')
        L_skirt = params.current_skirt_len_m;
    else
        L_skirt = 1.37; % 默认平衡长度
    end
    
    %% 2. 滑模控制律
    
    % 控制器参数 (需要根据响应速度调整)
    c = 100.0;        % 滑模面参数 (类似带宽)
    epsilon = 0.5;  % 切换增益 (抗扰动)
    k = 1.0;        % 比例增益 (加快趋近)
    delta = 0.05;   % 饱和函数边界层 (防抖震)
    
    % 误差定义
    e = z - target_z;
    dot_e = v - dot_target_z;
    
    % 滑模面 s
    s = dot_e + c * e;
    
    % 饱和函数 sat(s/delta)
    if abs(s/delta) <= 1
        sat_s = s/delta;
    else
        sat_s = sign(s);
    end
    
    % 计算所需垂直加速度 (Reaching Law)
    % s_dot = -epsilon*sat(s) - k*s
    % 而 s_dot = v_dot - dot_target_z_dot + c*dot_e
    % 故 v_dot_req = -epsilon*sat(s) - k*s - c*dot_e + dot_target_z_dot
    % (假设 dot_target_z_dot = 0)
    
    v_dot_req = -epsilon * sat_s - k * s - c * dot_e;
    
    %% 3. 动力学逆解 (Force -> Pressure)
    
    % 所需总升力 F_lift = m * (a + g)
    F_lift_req = m * (v_dot_req + g);
    
    % 物理约束：升力不能为负 (气垫只能推，不能拉)
    if F_lift_req < 0
        F_lift_req = 0;
    end
    
    % 所需平均压强 (Pa)
    P_req_pa = F_lift_req / A_total;
    
    % 物理约束：压强至少为大气压(相对压强0)，但为了风机计算不报错，设个小正值
    P_req_pa = max(P_req_pa, 10);
    
%% 4. 流量平衡逆解 (Pressure -> Flow -> RPM)
    
    % 单位转换准备 (原风机公式基于英制单位 psf, ft2, ft3/s)
    P_req_psf = P_req_pa / 47.8803; 
    
    % 4.1 估算围裙泄流面积 S (ft2)
    m_to_ft = 3.28084;
    h_gap_m = max(z - L_skirt, 0); % 气隙 = 船高 - 围裙长
    h_gap_ft = h_gap_m * m_to_ft;
    
    L_total_ft = L_total * m_to_ft;
    S_leak_ft2 = L_total_ft * h_gap_ft;
    
    % [重要] 移除之前的 0.5 最小限制，因为我们要加更准确的喷管流量了
    % S_leak_ft2 = max(S_leak_ft2, 0.5); 
    
    % 4.2 计算围裙泄流流量 Q_skirt (ft3/s)
    Q_skirt_cfs = S_leak_ft2 * 14.5 * sqrt(P_req_psf);
    
    % ================= [新增核心修正] =================
    % 4.3 计算姿态喷管泄流流量 Q_nozzle (ft3/s)
    % 物理模型中：Q_NOZ = -346 * sqrt(P_plenum)
    % 这里有左右两个喷管，且风机压强(Plenum)略高于气垫压强(Cushion)
    % 我们近似认为 P_plenum ≈ P_req_psf (或者略大，这里取1.0系数即可启动)
    % 两个喷管系数和 = 346 * 2 = 692
    
    Q_nozzle_cfs = 692 * sqrt(P_req_psf); 
    % =================================================
    
    % 4.4 考虑泵吸效应 (Q_pump = Area * dzdt)
    A_total_ft2 = A_total * m_to_ft^2;
    v_ft_s = v * m_to_ft;
    Q_pump_cfs = A_total_ft2 * v_ft_s;
    
    % 4.5 总需求流量 = 围裙泄流 + 喷管泄流 + 泵吸补偿
    Q_fan_req = Q_skirt_cfs + Q_nozzle_cfs + Q_pump_cfs;
    
    if Q_fan_req < 0, Q_fan_req = 0; end
    
    % 4.6 反解风机转速 RPM
    % 原风机公式: Q = (term1 + term2) * N / 2000
    % 其中 term1, term2 仅与压强 P 有关
    % 因此 N = 2000 * Q_req / (term1 + term2)
    
    % 风机特性系数 K_fan (仅与压强有关的部分)
    % 注意：风机模型里有 (P-300)，这里需保持一致
    delta_P_fan = P_req_psf - 300; 
    
    % 处理 sqrt 中的符号
    sqrt_val = sqrt(abs(delta_P_fan)) * sign(delta_P_fan);
    
    term1 = -1280 * sqrt_val;
    term2 = -31.6 * delta_P_fan;
    
    K_fan_pressure = term1 + term2;
    
    % 避免除零
    if abs(K_fan_pressure) < 1e-3
        K_fan_pressure = 1e-3;
    end
    
    % 计算 RPM
    cmd_rpm = 2000 * Q_fan_req / K_fan_pressure;
    
    % 最后的物理约束
    % 风机特性曲线在某些压强下 K_fan 可能是负的或者极小，导致计算出的RPM异常
    % 简单的保护逻辑：如果算出来是负的，说明这个压强点可能不在风机正常工作区，或者模型反了
    % 根据原公式，Q和N成正比。
    
    cmd_rpm = abs(cmd_rpm); 
    
    %% 5. 输出限幅
    cmd_rpm = min(max(cmd_rpm, 0), 2500); % 限制在 [0, 2500]
    
    %% 调试信息
    debug_info.e = e;
    debug_info.s = s;
    debug_info.v_dot_req = v_dot_req;
    debug_info.P_req_pa = P_req_pa;
    debug_info.Q_req = Q_fan_req;

end