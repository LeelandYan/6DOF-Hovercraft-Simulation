% 三自由度船运动学变换矩阵: 船体坐标系(u, v, r) -> 大地坐标系(dot_x, dot_y, dot_ψ)
function ROTATION_MATRIX = function_calculate_rotation_matrix(yaw)
    ROTATION_MATRIX = [cos(yaw) -sin(yaw) 0; 
                       sin(yaw)  cos(yaw) 0; 
                              0         0 1];
end