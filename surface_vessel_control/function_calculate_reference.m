function [reference, dot_reference, ddot_reference] = function_calculate_reference(time, SCENARIO)
    if SCENARIO == 1
%         reference = 1*[4*sin(0.02*time); 2.5*(1-cos(0.02*time)); 0.02*time];
%         dot_reference = 1*[0.02*4*cos(0.02*time); 0.02*2.5*sin(0.02*time); 0.02];
%         ddot_reference = 1*[-0.02*0.02*4*sin(0.02*time); 0.02*0.02*2.5*cos(0.02*time); 0];
        % 设定半径 R
        R = 4; 
        % 设定角频率 omega (控制跑一圈的快慢，原来是0.02)
        omega = 0.02; 
        
        % 1. 位置 (Reference)
        % 修改点：把原来的 4 和 2.5 都改成 R
        reference = [R*sin(omega*time); R*(1-cos(omega*time)); omega*time];
        
        % 2. 速度 (Dot Reference) - 必须要对应修改！
        % 求导规律：sin -> cos, cos -> -sin
        dot_reference = [omega*R*cos(omega*time); omega*R*sin(omega*time); omega];
        
        % 3. 加速度 (Ddot Reference) - 必须要对应修改！
        % 求导规律：cos -> -sin, sin -> cos
        ddot_reference = [-omega^2*R*sin(omega*time); omega^2*R*cos(omega*time); 0];
    
    else
        reference = [4*sin(0.02*time) + 2*cos(0.02*time); 2.5*(1-cos(0.02*time)) + 0.5*sin(0.02*time); 0.02*time + 0.5*sin(0.02*time)];
        dot_reference = [0.02*4*cos(0.02*time) - 0.02*2*sin(0.02*time); 0.02*2.5*sin(0.02*time) + 0.02*0.5*cos(0.02*time); 0.02 + 0.02*0.5*cos(0.02*time)];
        ddot_reference = [-0.02*0.02*4*sin(0.02*time) - 0.02*0.02*2*cos(0.02*time) ; 0.02*0.02*2.5*cos(0.02*time) - 0.02*0.02*0.5*sin(0.02*time); - 0.02*0.02*0.5*sin(0.02*time)];
    end
end