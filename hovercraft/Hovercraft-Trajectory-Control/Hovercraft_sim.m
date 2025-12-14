clc; clear all; close all;

%%
% Parametres
epsilon = 1e-6;
G.lambda_1 = 500;
G.lambda_2 = 500;
G.lambda_3 = 10;

% Temps final
tfinal = 20;

% Conditions initiales
u0 = 0;  v0 = 0; r0 = 0;

% 修改：将参考轨迹时间扩展到 tfinal
tmp = linspace(0, tfinal, 200); 
v_h = 5; v_l = 3; tmp_0 = 5;
u_ref_orig = v_l + (v_h - v_l) * tanh(tmp - tmp_0); 
v_ref_orig = v_l + (v_h - v_l) * tanh(tmp - tmp_0);

%% 无噪声仿真
disp('无噪声仿真');
relTol = 1e-4;    absTol = 1e-7;    
options = odeset('RelTol', relTol, 'AbsTol', absTol);

[ts, Xs] = ode23tb(@(t,X) hovercraft(t,X,G,tmp,u_ref_orig,v_ref_orig), [0 tfinal], [u0 v0 r0], options);

% 提取状态
u = Xs(:,1); v = Xs(:,2); r = Xs(:,3);

figure(1);
subplot(231); 
plot(ts, u, 'b', 'LineWidth', 2); hold on;
plot(tmp, u_ref_orig,'r', 'LineWidth', 2); hold off; grid; 
xlabel('时间 (s)'); ylabel('u'); 
title('纵向速度 u','FontWeight','bold');
legend('u', 'u_{ref}');

subplot(232); 
plot(ts, v, 'b', 'LineWidth', 2); hold on;
plot(tmp, v_ref_orig,'r', 'LineWidth', 2); hold off; grid; 
xlabel('时间 (s)'); ylabel('v');    
title('横向速度 v','FontWeight','bold');
legend('v', 'v_{ref}');

subplot(233); 
plot(ts, r, 'b', 'LineWidth', 2); grid; 
xlabel('时间 (s)'); ylabel('r'); 
title('艏摇角速度 r','FontWeight','bold');

% 插值参考轨迹用于绘制控制信号
u_ref = interp1(tmp, u_ref_orig, ts);
v_ref = interp1(tmp, v_ref_orig, ts);

du_ref = 1 - (u_ref).^2;
dv_ref = 1 - (v_ref).^2;
ddv_ref = -2.*v_ref.*dv_ref;

% 计算控制信号
tau_u = du_ref - G.lambda_1.*(u - u_ref) - v.*r;
tau_r = -(1./(u + epsilon)).*( ddv_ref - G.lambda_2.*(v - v_ref) - G.lambda_3.*(- u.*r - dv_ref) + v.*(r.^2) + tau_u.*r);
  
subplot(234); 
plot(ts, tau_u,'b', 'LineWidth', 2); grid; 
xlabel('时间 (s)'); ylabel('τ_u');
title('推力控制 τ_u');

subplot(235); 
plot(ts, tau_r,'b', 'LineWidth', 2); grid; 
xlabel('时间 (s)'); ylabel('τ_r');
title('转矩控制 τ_r');

%% 有噪声仿真
% disp('有噪声仿真');
% 
% % 使用原始参考轨迹
% [ts2, Xs2] = ode23tb(@(t,X) hovercraftBruit(t,X,G,tmp,u_ref_orig,v_ref_orig), [0 tfinal], [u0 v0 r0], options);
% 
% % 提取状态
% u2 = Xs2(:,1); v2 = Xs2(:,2); r2 = Xs2(:,3);
% 
% % 为绘图计算噪声
% Bruit_u = 0.96*sin(0.1*u2) + sin(10*u2); 
% Bruit_v = -0.96*sin(0.1*v2) + sin(10*v2); 
% Bruit_r = 0.96*sin(0.1*r2) + sin(10*r2);
% 
% figure(2);
% subplot(231); 
% plot(ts2, u2, 'b', 'LineWidth', 2); hold on;
% plot(tmp, u_ref_orig,'r', 'LineWidth', 2); hold off; grid; 
% xlabel('temps (s)'); ylabel('etat u'); 
% title('Etat u (avec bruit)','FontWeight','bold');
% legend('u', 'u_{ref}');
% 
% subplot(232); 
% plot(ts2, v2, 'b', 'LineWidth', 2); hold on;
% plot(tmp, v_ref_orig,'r', 'LineWidth', 2); hold off; grid; 
% xlabel('temps (s)'); ylabel('etat v');    
% title('Etat v (avec bruit)','FontWeight','bold');
% legend('v', 'v_{ref}');
% 
% subplot(233); 
% plot(ts2, r2, 'b', 'LineWidth', 2); grid; 
% xlabel('temps (s)'); ylabel('etat r'); 
% title('etat r (avec bruit)','FontWeight','bold');
% 
% % 插值参考轨迹
% u_ref2 = interp1(tmp, u_ref_orig, ts2);
% v_ref2 = interp1(tmp, v_ref_orig, ts2);
% 
% du_ref2 = 1 - (u_ref2).^2;
% dv_ref2 = 1 - (v_ref2).^2;
% ddv_ref2 = -2.*v_ref2.*dv_ref2;
% 
% % 计算控制信号
% tau_u2 = du_ref2 - G.lambda_1.*((u2 + Bruit_u) - u_ref2) - (v2 + Bruit_v).*(r2 + Bruit_r);
% tau_r2 = -(1./((u2 + Bruit_u) + epsilon)).*( ddv_ref2 - G.lambda_2.*((v2 + Bruit_v) - v_ref2) - G.lambda_3.*(-(u2 + Bruit_u).*(r2 + Bruit_r) - dv_ref2) + (v2 + Bruit_v).*((r2 + Bruit_r).^2) + tau_u2.*(r2 + Bruit_r));
%   
% subplot(234); 
% plot(ts2, tau_u2,'b', 'LineWidth', 2); grid; 
% xlabel('Time (s)'); ylabel('TauU');
% title('Control TauU (avec bruit)');
% 
% subplot(235); 
% plot(ts2, tau_r2,'b', 'LineWidth', 2); grid; 
% xlabel('Time (s)'); ylabel('TauR');
% title('Control TauR (avec bruit)');


%% 函数定义
function [dotX] = hovercraft(t, X, G, tmp, u_ref, v_ref)
    u = X(1,:); v = X(2,:); r = X(3,:);
    epsilon = 1e-6;
    
    u_ref = interp1(tmp, u_ref, t);
    v_ref = interp1(tmp, v_ref, t);
    
    du_ref = 1 - (u_ref)^2;
    dv_ref = 1 - (v_ref)^2;
    ddv_ref = -2*v_ref*dv_ref;
    
    tau_u = du_ref - G.lambda_1.*(u - u_ref) - v.*r;
    tau_r = -(1./(u + epsilon)).*( ddv_ref - G.lambda_2.*(v - v_ref) - G.lambda_3.*(- u.*r - dv_ref) + v.*(r.^2) + tau_u.*r);
    
    dotu = v.*r + tau_u;
    dotv = -u.*r;
    dotr = tau_r;
    
    dotX = [dotu dotv dotr]';
end

function [dotX] = hovercraftBruit(t, X, G, tmp, u_ref, v_ref)
    u = X(1,:); v = X(2,:); r = X(3,:);
    epsilon = 1e-6;
    
    Bruit_u = 0.96*sin(0.1*u) + sin(10*u); 
    Bruit_v = -0.96*sin(0.1*v) + sin(10*v); 
    Bruit_r = 0.96*sin(0.1*r) + sin(10*r);
    
    u_ref = interp1(tmp, u_ref, t);
    v_ref = interp1(tmp, v_ref, t);
    
    du_ref = 1 - (u_ref)^2;
    dv_ref = 1 - (v_ref)^2;
    ddv_ref = -2*v_ref*dv_ref;
    
    tau_u = du_ref - G.lambda_1.*((u + Bruit_u) - u_ref) - (v + Bruit_v).*(r + Bruit_r);
    tau_r = -(1./((u + Bruit_u) + epsilon)).*( ddv_ref - G.lambda_2.*((v + Bruit_v) - v_ref) - G.lambda_3.*(-(u + Bruit_u).*(r + Bruit_r) - dv_ref) + (v + Bruit_v).*((r + Bruit_r).^2) + tau_u.*(r + Bruit_r));
    
    dotu = (v + Bruit_v).*(r + Bruit_r) + tau_u;
    dotv = -(u + Bruit_u).*(r + Bruit_r);
    dotr = tau_r;
    
    dotX = [dotu dotv dotr]';
end