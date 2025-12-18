% 该矩阵用于将船体坐标系的螺旋桨推力转换为大地坐标系下的船体加速度。

function CONTROL_MATRIX = function_calculate_control_matrix(ROTATION_MATRIX, INERTIA_MATRIX_INVERSE)
%     disp('*******************eeeeeeeee*******************');
%     ROTATION_MATRIX
%     INERTIA_MATRIX_INVERSE

    CONTROL_MATRIX = ROTATION_MATRIX*INERTIA_MATRIX_INVERSE;
%     disp('**************************************');
end