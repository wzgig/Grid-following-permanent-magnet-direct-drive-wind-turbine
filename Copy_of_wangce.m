function dotx = Copy_of_wangce(t,x)
%轴系方程+定子电压方程q轴d轴+d轴控制，采用有名值,d轴PI参数没有调好，初始不稳定
%参考曹明峰
%%%%%%%%%%%%%%%%%%%%%%%%%% Rename %%%%%%%%%%%%%%%%%%%%%%%%
omega_m = x(1);
i_sq       = x(2);
i_sd       = x(3);
x1         = x(4);
x2         = x(5);
y1         = x(6);
u_sd       = x(7);
u_sq       = x(8);
xpll       = x(9);
thetapll   = x(10);
U_dc   = x(11);
z1   = x(12);
i_gd   = x(13);
i_gq   = x(14);
z2   = x(15);
z3   = x(16);

%%%%%%%%%%%%%%%%%%%%%%%%%%  Paras  %%%%%%%%%%%%%%%%%%%%%%%%
f    = 50;
Fnom    = 50;
Np  = 12;
w_g = 2*pi*f;
Sbase = 0.3e6;
Pnom = 0.3e6;
Vbase = 690;
Vdqbase = Vbase/sqrt(3);
Ibase = Sbase/(sqrt(3)*Vbase);
wbase = w_g;
Zbase = Vbase^2/Sbase;
Lbase = Zbase/w_g;%电感基值
Cbase = 1/(wbase*Zbase);%电容基值
Tbase = Sbase/(w_g/Np);
U_dcref = 1800;

R_g = 0.05;%pu
L_g = 7e-3;%pu
C_dc = 1e-3;
% Psibase = Vbase/wbase;% 基准磁链
% tbase = 1;% 时间基值/wbase
% Kv_p_base = Ibase/Vbase;%电压环比例增益
% Kv_i_base = (Ibase/Vbase)/tbase;%电压环积分增益
% Ki_p_base = Vbase/Ibase;%电流环比例增益
% Ki_i_base = (Vbase/Ibase)/tbase;%电流环积分增益
% Dbase = Tbase/(wbase * wbase);%自阻尼基准值

Psi_f = 3.88889;%ac           % 转子磁链 
L_sd = 1.8e-3;%ac          % d轴电感 
L_sq = 1.8e-3;%ac          % q轴电感 
R_s = 0.025;%ac         % 定子电阻
beta = 0;              % 桨距角，单位：度（固定值
pitch = 0;
v_w = 12;               % 风速 (m/s)
% T_J = 8;               % 惯性时间常数 (s)
J = 60;%ac
D_m = 0.078;%ac           % 自阻尼系数
R_t = 33.05;
rho = 1.12;


K_ptrq = 6;
K_itrq = 0.6;
% K_p1 = 0.03*2129.991/469.486;%u_sd
% K_i1 = (0.03/0.05)*2129.991/469.486;%u_sd
% K_p2 = 0.03*2129.991/469.486;%u_sq
% K_i2 = (0.03/0.05)*2129.991/469.486;%u_sq

% q 轴电流环比例增益KPq = 0.0029
% q 轴电流环积分增益KIq = 0.1651
K_p1 = 1;
K_i1 = 12;
K_p2 = 1;
K_i2 = 12;
K_p3 = 10;
K_i3 = 220;
% K_p1 = 2129.991/469.486;%u_sd
% K_i1 = (1/0.05)*2129.991/469.486;%u_sd
% K_p2 = 2129.991/469.486;%u_sq
% K_i2 = (1/0.05)*2129.991/469.486;%u_sq

K_pU_dc = 1;
K_iU_dc = 10;
K_p5i_gd = 1;
K_i5i_gd = 15;
K_p6i_gq = 1;
K_i6i_gq = 15;

K_pPLL = 250;
K_iPLL = 3200;
T_d = 1/6000;
T_m = 0.001;
T_trq = 60;

% 目标有功功率和无功功率（单位：p.u.）
P = 0.8;     % 发电模式下为负值，有功功率（p.u.）
Q = 0.3;     % 无功功率（p.u.）
V = 0.9998;   % 母线电压幅值（p.u.）
xi = 0.6286;  % 电压相角（rad）

% 计算 a-b 坐标系下的电流参考值
S_D0 = P/V;                      % 有功功率对应的电流分量基础（p.u.）
S_Q0 = Q/V;                      % 无功功率对应的电流分量基础（p.u.）
S_ab0 = S_D0 + 1j*S_Q0;          % αβ坐标系下的复功率表示（p.u.）
i_abs = abs(S_ab0);              % 电流幅值（p.u.）
ui_argdiff = angle(S_ab0);       % 电压与电流之间的相角差（rad）
i_arg = xi - ui_argdiff;         % 电流在αβ坐标系中的相位（rad，xi为电网电压α轴相位）
i_ab = i_abs * exp(1i * i_arg);  % αβ坐标下的复电流（i_ab = i_α + j*i_β）
i_b = imag(i_ab);                % β轴电流分量（p.u.）
i_a = real(i_ab);                % α轴电流分量（p.u.）

v_ab = V * exp(1i * xi);         % αβ坐标下的复电压（v_ab = v_α + j*v_β，xi为电压α轴相位）
v_a = real(v_ab);                % α轴电压分量（p.u.，与电网A相电压对齐）
v_b = imag(v_ab);                % β轴电压分量（p.u.，滞后α轴90°电角度）
%%%%%%%%%%%%%%%%%%%%%%%%%% Ref    %%%%%%%%%%%%%%%%%%%%%%%%
% u_sd = -0.35*Vdqbase;u_sq = -1.15*Vdqbase;
% u_sd = 0;
% u_sq = 320;
% i_sd = 0;
% if t<2
% v_w =12;
% end
% if t>2
% v_w = 14;
% end
%%%%%%%%%%%%%%%%%%%%%%%%%% Algebra %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% algebra
% syms omega_m i_sq i_sd x1 x2 y1 u_sd u_sq xpll thetapll U_dc z1 i_gd i_dq z2 z3;

n = omega_m * 30/pi;
omega_e = omega_m*Np;
omega_r = omega_m * (Np/(Fnom*2*pi));

% i_sd = (Psi_sd - Psi_f)/L_sd;
% i_sq = Psi_sq/L_sq;

% Te = P_s/omega_e;

Te = 1.5*Np*(Psi_f*i_sq);
Tem = Te * (1/(Np*Pnom/(2*pi*Fnom)));

% u1 = (omega_r*(3.8124/1.1)*33.05)/v_w;
% u2 = pitch;
% Cp = (0.44-0.0167*u2)*sin(pi*(u1-3)/15-0.3*u2)-0.00184*(u1-3)*u2;
% Tm = (v_w^3 * Cp *(0.5*1.12*pi*33.05^2)*(1/(1.5e6)))/omega_r;
% Tm_Nm = Tm*(-Np*Pnom/(2*pi*Fnom));
% P_w = Tm_Nm * omega_m;
u1 = 0.5*pi;
u2 = 1.2;
u3 = 14;%R
lambda = (u3*omega_m)/(v_w);
lambda_i = 1/(1/(lambda+0.08*beta)-0.035/(beta^3+1));
Cp = 0.51763*(116/lambda_i-0.4*beta-5)*exp(-21/lambda_i)+0.006795*lambda;
P_m = u1*u2*(u3^2)*(v_w^3)*Cp;
T_m = (-1)*(P_m)/omega_m;
omega_best = 8.1*v_w/u3;%最佳叶尖速比
n_ref = omega_best*30/pi;

i_sdref = 0;

%dy1 = n_ref - n;
% dx1 = i_sdref - i_sd;
% dx2 = i_sqref - i_sq;
i_sqref = K_p3*(n_ref - n) + K_i3*y1;
% i_sqref = -607;
u_sqref = K_p2*(i_sqref - i_sq)+K_i2*x2+omega_e*Psi_f+omega_e*L_sd*i_sd;
u_sdref = K_p1*(i_sdref - i_sd)+K_i1*x1-omega_e*L_sq*i_sq;


P_s = 1.5*(u_sd * i_sd + u_sq * i_sq);%(2-17)
Q_s = 1.5*(u_sq * i_sd - u_sd * i_sq);%(2-18)
P_e = Te * omega_m;




%% ========== 坐标变换（GSC的电压观测） ==========
% PLL
theta_g = thetapll;         % 网侧PLL估计的同步角
% 构造变换矩阵，将d-q轴量转换为定轴坐标系（用于测量或控制）
Tg_DQdq = [cos(theta_g) -sin(theta_g); sin(theta_g) cos(theta_g)];
vg_dq_real = Vbase * Tg_DQdq * [v_a; v_b];  % GSC中的实际电压
% 变换后各轴上的电压值提取
v_g_d = vg_dq_real(1);          % 网侧d轴电压
v_g_q = -vg_dq_real(2);         % 网侧q轴电压

i_gdref = K_pU_dc*(U_dc - U_dcref) + K_iU_dc*z1;
i_gqref = 0;

u_gd = -K_p5i_gd*(i_gdref - i_gd)-K_i5i_gd*z2+v_g_d+i_gq*w_g*L_g;
u_gq = -K_p6i_gq*(i_gqref - i_gq)-K_i6i_gq*z3+v_g_q-i_gd*w_g*L_g;

P_g = 1.5*(v_g_d * i_gd + v_g_q * i_gq);%(2-17)
P_dc = 1.5*(u_gd* i_gd + u_gq * i_gq);%(2-17)
Q_g = 1.5*(v_g_q * i_gd - v_g_d * i_gq);%(2-18)

%%%%%%%%%%%%%%%%%%%%%%% Differential %%%%%%%%%%%%%%%%%%%%%

    dotx = [
        (Te - T_m - D_m*omega_m) / J;
        % u_sd - R_s*i_sd + omega_e*Psi_sq;
        % u_sq - R_s*i_sq - omega_e*Psi_sd;
        (u_sq - R_s*i_sq - omega_e*L_sd*i_sd - omega_e*Psi_f)/L_sq;
        (u_sd - R_s*i_sd + omega_e*L_sq*i_sq)/L_sd;
        i_sdref - i_sd;
        i_sqref - i_sq;
        n_ref - n;
        (u_sdref - u_sd)/T_d;
        (u_sqref - u_sq)/T_d;
        v_g_q;
        K_pPLL*v_g_q+K_iPLL*xpll+w_g;
        (P_s - P_dc)/(C_dc*U_dc);
        U_dc - U_dcref;
        (v_g_d - u_gd + w_g*L_g*i_gq)/L_g;
        (v_g_q - u_gq - w_g*L_g*i_gd)/L_g;
        i_gdref - i_gd;
        i_gqref - i_gq;
    ];
end

