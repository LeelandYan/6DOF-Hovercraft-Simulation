% 输入附体坐标系的力 tau -> 输出大地坐标系的运动状态 X1, X2
% 大地坐标下：X1(1): x坐标，X1(2): y坐标，X1(3): ψ航向角
% 大地坐标下：X2(1): dot_x在x轴方向移动的速度，X2(2): dot_y在 y 轴方向移动的速度，X2(3): dot_ψ航向角的变化率
% 船体坐标下：U(1): u纵向速度，U(2): v横向速度，U(3): r转首速率
function [X1_NEW, X2_NEW, U_NEW, disturbance] = function_calculate_ship_motion(time, X1, X2, U, tau, PARAMETERS)
    yaw = X1(3); % 当前船头朝向 (航向角)
    u = U(1);    % 纵向速度 (船头方向跑多快)
    v = U(2);    % 横向速度 (侧向漂移多快)
    r = U(3);    % 转向角速度 (转弯有多快)

    ROTATION_MATRIX = function_calculate_rotation_matrix(yaw); 
    CONTROL_MATRIX = function_calculate_control_matrix(ROTATION_MATRIX, PARAMETERS.SHIP.INERTIA_MATRIX_INVERSE);    
    KNOWN_DYNAMICS = function_calculate_known_dynamics(u, v, r, ROTATION_MATRIX, CONTROL_MATRIX);
    
    F = function_calculate_disturbance(time, u, v, r) - function_calculate_model_uncertain(time, u, v, r);
    disturbance = -CONTROL_MATRIX*F;
    dot_X1 = X2; % 大地坐标下速度
    dot_X2 = CONTROL_MATRIX*tau + KNOWN_DYNAMICS + disturbance;  % 大地坐标下加速度   

    X1_NEW = X1 + dot_X1*PARAMETERS.SIMULATION.SAMPLING_TIME;
    X2_NEW = X2 + dot_X2*PARAMETERS.SIMULATION.SAMPLING_TIME;

    U_NEW = ROTATION_MATRIX'*X2;
end