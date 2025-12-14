close all
clear all
clc

%%
psf_to_pa = 47.88;
% ==================== 参数设置 ====================
S_SI = [3; 3; 3; 3]; % 围裙泄流面积 (m²)
S = S_SI * 10.76391; % 围裙泄流面积 (ft²)
  

% 风机转速 (rpm)
% 文档基准转速为 2000 rpm
N_FAN1 = 1800;        
N_FAN2 = 1800;

% 牛顿-拉夫逊法参数
max_iter = 30;        % 最大迭代次数
tolerance = 10;       % 收敛容差 

% 初始压强值 (psf - 磅/平方英尺)
P = zeros(6, 1);
P(1:4) = 0; % 气垫气室初始压力 (PGC)
P(5:6) = 0; % 风机气室初始压力 (PFAN)

%%
% ==================== 牛顿-拉夫逊迭代求解 ====================
fprintf('牛顿-拉夫逊迭代...\n');

P_history = zeros(max_iter, 6);
for iter = 1:max_iter
    % ========== 计算各流量 (ft³/s) ==========
    % 公式来源：文档 Eq (49)
    term1_fan1 = -1280 * sqrt(abs(P(5) - 300)) * sign(P(5) - 300);
    term2_fan1 = -31.6 * (P(5) - 300);
    Q_FAN1 = (term1_fan1 + term2_fan1) * (N_FAN1 / 2000);
    
    term1_fan2 = -1280 * sqrt(abs(P(6) - 300)) * sign(P(6) - 300);
    term2_fan2 = -31.6 * (P(6) - 300);
    Q_FAN2 = (term1_fan2 + term2_fan2) * (N_FAN2 / 2000);
    
    % 风机气室流向气垫气室流量
    % 公式来源：文档 Eq (48)
    Q_INC1 = 589 * sqrt(abs(P(5) - P(1))) * sign(P(5) - P(1));
    Q_INC2 = 589 * sqrt(abs(P(5) - P(2))) * sign(P(5) - P(2));
    Q_INC3 = 589 * sqrt(abs(P(6) - P(3))) * sign(P(6) - P(3));
    Q_INC4 = 589 * sqrt(abs(P(6) - P(4))) * sign(P(6) - P(4));
    
    % 喷管流向风机气室流量
    % 公式来源：文档 Eq (48)
    Q_NOZ1 = -346 * sqrt(abs(P(5))) * sign(P(5));
    Q_NOZ2 = -346 * sqrt(abs(P(6))) * sign(P(6));
    
    % 气垫气室间泄流量
    % 公式来源：文档 Eq (48)
    Q_IC1 = 675 * sqrt(abs(P(4) - P(1))) * sign(P(4) - P(1));
    Q_IC2 = 338 * sqrt(abs(P(1) - P(2))) * sign(P(1) - P(2));
    Q_IC3 = 675 * sqrt(abs(P(2) - P(3))) * sign(P(2) - P(3));
    Q_IC4 = 338 * sqrt(abs(P(3) - P(4))) * sign(P(3) - P(4));
    
    % 气垫气室围裙底部的泄流量 (流向大气)
    % 公式来源：文档 Eq (48) 
    Q1 = -S(1) * 14.5 * sqrt(abs(P(1))) * sign(P(1));
    Q2 = -S(2) * 14.5 * sqrt(abs(P(2))) * sign(P(2));
    Q3 = -S(3) * 14.5 * sqrt(abs(P(3))) * sign(P(3));
    Q4 = -S(4) * 14.5 * sqrt(abs(P(4))) * sign(P(4));
    
    % ========== 构建残差向量 F(P) ==========
    F = zeros(6, 1);
    
    % 风扇气室方程 (Eq 46)
    F(5) = -Q_INC1 - Q_INC2 + Q_NOZ1 + Q_FAN1;
    F(6) = -Q_INC3 - Q_INC4 + Q_NOZ2 + Q_FAN2;
    
    % 气垫气室方程 (Eq 47)
    F(1) = Q_INC1 + Q_IC1 - Q_IC2 + Q1;
    F(2) = Q_INC2 + Q_IC2 - Q_IC3 + Q2;
    F(3) = Q_INC3 + Q_IC3 - Q_IC4 + Q3;
    F(4) = Q_INC4 + Q_IC4 - Q_IC1 + Q4;
    
    % ========== 检查收敛 ==========
    residual_norm = norm(F);
%     fprintf('残差向量 = %.6e\n', F);
    fprintf('迭代 %3d: 残差二范数 = %.6e\n', iter, residual_norm);
    
    if residual_norm < tolerance
        fprintf('\n收敛成功！\n');
%         break;
    end
    
    % ========== 计算雅可比矩阵 ==========
    J = compute_jacobian(P, S, N_FAN1, N_FAN2);
    
    % ========== 牛顿更新 ==========
    % 加入阻尼防止震荡
    lambda = 0.7; 
    delta_P = -J \ F;
    
    P = P + lambda * delta_P;
    P_history(iter, :) = P' * psf_to_pa;
end


%%
% ==================== 输出结果 ====================

fprintf('\n========== 最终压强值 (pa) ==========\n');
fprintf('P1 = %.2f pa (气垫1)\n', P(1) * psf_to_pa);
fprintf('P2 = %.2f pa (气垫2)\n', P(2) * psf_to_pa);
fprintf('P3 = %.2f pa (气垫3)\n', P(3) * psf_to_pa);
fprintf('P4 = %.2f pa (气垫4)\n', P(4) * psf_to_pa);
fprintf('P_cushion_avg ≈ %.0f Pa\n', mean(P(1:4)) * psf_to_pa);



%% ==================== 绘图：四个气垫压强变化 ====================
figure;
sgtitle('气室压强变化'); 

% 定义每个子图的中文标题
titles = {'P1 气垫气室1 ', ...
          'P2 气垫气室2 ', ...
          'P3 气垫气室3 ', ...
          'P4 气垫气室4 ', ...
          };

for i = 1:4
    subplot(2, 2, i); 
    plot(P_history(1:iter, i), 'LineWidth', 1.5);
    
    title(titles{i});
    xlabel('迭代次数');
    ylabel('压强 (pa)');
    grid on;

end

%%
% ==================== 雅可比矩阵计算函数 ====================
function J = compute_jacobian(P, S, N_FAN1, N_FAN2)
    epsilon = 1e-4; % 扰动量
    n = length(P);
    J = zeros(n, n);
    F0 = residual_function(P, S, N_FAN1, N_FAN2);
    
    for i = 1:n
        P_perturb = P;
        P_perturb(i) = P_perturb(i) + epsilon;
        F_perturb = residual_function(P_perturb, S, N_FAN1, N_FAN2);
        J(:, i) = (F_perturb - F0) / epsilon;
    end
end

%%
% ==================== 残差函数 ( ====================
function F = residual_function(P, S, N_FAN1, N_FAN2)
    % 风机流量 (Eq 49)
    Q_FAN1 = (-1280 * sqrt(abs(P(5) - 300)) * sign(P(5) - 300) ...
        - 31.6 * (P(5) - 300)) * N_FAN1 / 2000;
    Q_FAN2 = (-1280 * sqrt(abs(P(6) - 300)) * sign(P(6) - 300) ...
        - 31.6 * (P(6) - 300)) * N_FAN2 / 2000;
    
    % 内部流 (Eq 48)
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
    
    % 泄流 (Eq 48)
    Q1 = -S(1) * 14.5 * sqrt(abs(P(1))) * sign(P(1));
    Q2 = -S(2) * 14.5 * sqrt(abs(P(2))) * sign(P(2));
    Q3 = -S(3) * 14.5 * sqrt(abs(P(3))) * sign(P(3));
    Q4 = -S(4) * 14.5 * sqrt(abs(P(4))) * sign(P(4));
    
    F = zeros(6, 1);
    % Eq 46
    F(5) = -Q_INC1 - Q_INC2 + Q_NOZ1 + Q_FAN1;
    F(6) = -Q_INC3 - Q_INC4 + Q_NOZ2 + Q_FAN2;
    % Eq 47
    F(1) = Q_INC1 + Q_IC1 - Q_IC2 + Q1;
    F(2) = Q_INC2 + Q_IC2 - Q_IC3 + Q2;
    F(3) = Q_INC3 + Q_IC3 - Q_IC4 + Q3;
    F(4) = Q_INC4 + Q_IC4 - Q_IC1 + Q4;
end