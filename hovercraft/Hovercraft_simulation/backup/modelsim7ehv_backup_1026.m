function dXdt = modelsim7ehv(t, X) 
    %% -------- 状态量定义 -------- 
    x = X(1); 
    y = X(2); 
    phi = X(3); % roll 
    psi = X(4); % yaw 
    u = X(5); 
    v = X(6); 
    p = X(7); % roll rate 
    r = X(8); % yaw rate 
    
    % -------- 初始参数定义 -------- 
    m = 40000; % kg 
    Jx = 2.5e5; % roll inertia 
    Jz = 1.8e6; % yaw inertia
    Xu = 30000; % surge damping 
    Yv = 50000; % sway damping 
    Mp = 200000;% roll damping 
    Nr = 100000; % yaw damping 
    
    %% -------- 直接输入推力和舵角值 -------- 
    Fx = thrust_command(t); % surge thrust 
    delta = rudder_command(t); % rudder angle [rad] 
    %% -------- 计算舵力 -------- 
    rho_a = 1.29; 
    CxR = 0.52; 
    CyR = 0.52; 
    SR = 1.5; 
    Sd = 10.0; 
    xR = -2.2; 
    zR = -1.0; 
    Vw = 10 * 0.5144; % 真风速 
    beta_w = deg2rad(0); % 真风向 
    Va = sqrt((u + Vw*cos(beta_w - psi))^2 + ... 
        (v + Vw*sin(beta_w - psi))^2); 
    beta_a = atan2(v + Vw*sin(beta_w - psi), ... 
        u + Vw*cos(beta_w - psi)); 
    % Rudder inflow velocity 
    FxP = Fx; 
    if FxP >= 0 
        vR = Va*cos(beta_a) + sqrt(FxP/(rho_a*Sd)); 
    else 
        vR = Va*cos(beta_a) - sqrt(abs(FxP)/(2*rho_a*Sd)); 
    end 
    
    % Apply delta influence 
    FX_r = -CxR * rho_a * vR^2 * SR * cos(delta); 
    FY_r = CyR * rho_a * vR^2 * SR * sin(delta); 
    
    % Total from two rudders 
    Fx_R = 2 * FX_r; 
    Fy_R = 2 * FY_r; 
    Mx_R = -2 * FY_r * zR; 
    Mz_R = 2 * FY_r * xR; 
    %% % 阻力计算，需要修改Fx_d、Fy_d、Mx_d、Mz_d 
    Fx_d = -Xu * u; 
    Fy_d = -Yv * v; 
    Mx_d = -Mp * p; 
    Mz_d = -Nr * r; 
    
    % 四个自由度总合力合力矩 
    Fx_total = Fx_d + Fx_R + Fx; 
    Fy_total = Fy_d + Fy_R ; 
    Mx_total = Mx_d + Mx_R; 
    Mz_total = Mz_d + Mz_R; 
    
    %% -------- 运动学方程 -------- 
    xdot = cos(psi)*u - sin(psi)*cos(phi)*v; 
    ydot = sin(psi)*u + cos(psi)*cos(phi)*v; 
    phidot = p; 
    psidot = r * cos(phi); 
    
    %% -------- 动力学方程 -------- 
    udot = (Fx_total + m*v*r) / m; 
    vdot = (Fy_total - m*u*r) / m; 
    pdot = Mx_total / Jx; 
    rdot = Mz_total / Jz; 
    
    % 合并输出 
    dXdt = [xdot; ydot; phidot; psidot; udot; vdot; pdot; rdot]; 
end 

%% ====== 人为确定推力和舵角大小(开环) ====== 
function Fx = thrust_command(t) % Constant thrust 
    if t < 50 
        Fx = 2.5e5; % Fx = 2.5e5 
    else 
        Fx = 2.5e5; 
    end 
end 

function delta = rudder_command(t) 
    if t < 20 
        delta = 0; 
    else 
        delta = deg2rad(3); 
    end 
end





%%
function dXdt = modelsim7ehv(t, X)

%% -------- 状态量定义 --------
x   = X(1);
y   = X(2);
phi = X(3);    % roll
psi = X(4);    % yaw
u   = X(5);
v   = X(6);
p   = X(7);    % roll rate
r   = X(8);    % yaw rate

%% -------- 初始参数定义 --------
m  = 40000;      % kg
Jx = 2.5e5;      % roll inertia 
Jz = 1.8e6;      % yaw inertia
%% -------- ACV Physical Parameters (From Paper) --------
% Densities
rho_a = 1.29;        % air
rho_w = 1025;        % water

% Gravity
g = 9.8;

% Air drag coefficients
Cxa = 1.05;
Cya = 1.00;
Cmx_a = 0.90;
Cmz_a = 1.02;

% Wave making & skirt drag coefficients
Cwm = 0.3;           % assume average if Nwm(u) unavailable
Csk = 1.2;

% Cushion pressure
Pc = 2000;

% Areas (Table 1)
SPP = 45;
SLP = 93;
SHP = 230;
Sc  = 205;

% Dimensions (Table 1)
Bc = 8.9;
lc = 23.6;
lsk = 65;

% Hover height & skirt geometry
h   = 1.0;
hm  = 2.4;
h0  = 1.3;
Hhov = 5.9;

% Fan flow
Q = 145;

xa = 2.8; ya = 0; za = 0.3;
xwm = 2.8; ywm = 0; zwm = 1.5;
xsk = 2.8; ysk = 0; zsk = 1.0;
xc = 2.8; yc = 0; zc = 1.5;
xm = 3.0; ym = 0; zm = -0.8;


%% -------- 直接输入推力和舵角值 --------
Fx = thrust_command(t);   % surge thrust
delta = rudder_command(t); % rudder angle [rad]

%% -------- 计算舵力 --------
rho_a = 1.29;
CxR = 0.52;
CyR = 0.52;

SR = 1.5;
Sd = 10.0;

xR = -2.2;
zR = -1.0;

Vw = 10 * 0.5144; % 真风速
beta_w = deg2rad(0); % 真风向
Va = sqrt((u + Vw*cos(beta_w - psi))^2 + ...
          (v + Vw*sin(beta_w - psi))^2);
beta_a = atan2(v + Vw*sin(beta_w - psi), ...
               u + Vw*cos(beta_w - psi));

% Rudder inflow velocity
FxP = Fx;
if FxP >= 0
    vR = Va*cos(beta_a) + sqrt(FxP/(rho_a*Sd));
else
    vR = Va*cos(beta_a) - sqrt(abs(FxP)/(2*rho_a*Sd));
end

% Apply delta influence
FX_r = -CxR * rho_a * vR^2 * SR * cos(delta);
FY_r =  CyR * rho_a * vR^2 * SR * sin(delta);

% Total from two rudders
Fx_R = 2 * FX_r;
Fy_R = 2 * FY_r;
Mx_R = -2 * FY_r * zR;
Mz_R =  2 * FY_r * xR;


%% 
% 阻力计算，需要修改Fx_d、Fy_d、Mx_d、Mz_d 
% 侧漂角
if abs(u) < 0.1
    beta = 0;
else
    beta = atan2(v, u);
end

% Air profile drag (Fx_a, Fy_a, Mx_a, Mz_a)
Fx_a = -0.5 * rho_a * Va^2 * Cxa  * SPP;
Fy_a = -0.5 * rho_a * Va^2 * Cya  * SLP;

Mx_a = -0.5 * rho_a * Va^2 * Cmx_a * SPP * lsk; % assume arm = lsk
Mz_a = -0.5 * rho_a * Va^2 * Cmz_a * SHP * Hhov;

% Air momentum drag
Fm = rho_a * Va * Q;

% Wave making drag
Fwm = Cwm * Pc^2 * Bc / (rho_w * g);

% Skirt drag
Fsk = 0.5 * rho_w * Va^2 * Csk * (h/lsk)^(-0.34) * sqrt(Sc) ...
     + (2.8167*(Pc/lc)^0.259 -1)*Fwm;

% Cushion force
Fc0  = 2 * lc * Pc * h0;                 % 基线升力（只支撑，不进平面力矩）
dFc  = 2 * lc * Pc * (0.5*Bc * phi);     % 随 phi 变化的增量升力
Fc   = Fc0 + dFc;    

% Roll restoring moment
MG = g * hm * m * phi;  

% 总阻力阻力力矩
Fx_d = Fx_a - Fm*cos(beta) - Fwm*cos(beta) - Fsk*cos(beta);
Fy_d = Fy_a - Fm*sin(beta) - Fwm*sin(beta) - Fsk*sin(beta);

Mx_d = Mx_a - Fm*zm*sin(beta) - Fwm*zwm*sin(beta) ...
       - Fsk*zsk*sin(beta) + MG; % assume small arms = 0.1m if unknown

Mz_d = Mz_a + dFc*xc ...
     - Fm*xm*sin(beta) - Fwm*xwm*sin(beta) - Fsk*xsk*sin(beta);

% 阻尼修正
D_u = 8e4;   % surge
D_v = 1.5e5; % sway
D_r = 8e5;   % yaw
D_p = 3e5;   % roll

Fx_d = Fx_d - D_u*u;
Fy_d = Fy_d - D_v*v;
Mx_d = Mx_d - D_p*p;
Mz_d = Mz_d - D_r*r;

% 四个自由度总合力合力矩
Fx_total = Fx_d + Fx_R + Fx;
Fy_total = Fy_d + Fy_R ;
Mx_total = Mx_d + Mx_R;
Mz_total = Mz_d + Mz_R;


%% -------- 运动学方程 --------
xdot =  cos(psi)*u - sin(psi)*cos(phi)*v;
ydot =  sin(psi)*u + cos(psi)*cos(phi)*v;
phidot = p;
psidot = r * cos(phi);

%% -------- 动力学方程 --------
udot = (Fx_total + m*v*r) / m;
vdot = (Fy_total - m*u*r) / m;
pdot = Mx_total / Jx;
rdot = Mz_total / Jz;

% 合并输出
dXdt = [xdot; ydot; phidot; psidot; udot; vdot; pdot; rdot];
end


%% ====== 人为确定推力和舵角大小(开环) ======
function Fx = thrust_command(t)
    % Constant thrust 
    if t < 50
        Fx = 2.5e5; % Fx = 2.5e5
    else
        Fx = 2.5e5;
    end
end

function delta = rudder_command(t)
    if t < 20
        delta = 0; 
    else
        delta = deg2rad(0);
    end
end

