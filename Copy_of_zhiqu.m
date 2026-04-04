function dotx = Copy_of_zhiqu(t,x)
%轴系方程+定子电压方程q轴d轴+d轴控制，采用有名值,d轴PI参数没有调好，初始不稳定
%参考曹明峰
%%%%%%%%%%%%%%%%%%%%%%%%%% Rename %%%%%%%%%%%%%%%%%%%%%%%%
omega_m = x(1);
i_sq       = x(2);
i_sd       = x(3);
x1       = x(4);
x2       = x(5);
y1       = x(6);
u_sd       = x(7);
u_sq       = x(8);


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
U_dc = 1800;

R_g = 0.003;%pu
L_g = 0.3;%pu
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
R_t = 14;
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
K_p3 = 100;
K_i3 = 220;
% K_p1 = 2129.991/469.486;%u_sd
% K_i1 = (1/0.05)*2129.991/469.486;%u_sd
% K_p2 = 2129.991/469.486;%u_sq
% K_i2 = (1/0.05)*2129.991/469.486;%u_sq

K_p4 = 1.2;
K_i4 = 10;
K_p5 = 1;
K_i5 = 30;
K_p6 = 1.2;
K_i6 = 10;
K_p7 = 1;
K_i7 = 30;
K_pPLL = 50;
K_iPLL = 1000;
T_d = 1/6000;
T_m = 0.001;
T_trq = 60;

%%%%%%%%%%%%%%%%%%%%%%%%%% Ref    %%%%%%%%%%%%%%%%%%%%%%%%
% u_sd = -0.35*Vdqbase;u_sq = -1.15*Vdqbase;
% u_sd = 0;
% u_sq = 320;
% i_sd = 0;

omega_best = 8.1*v_w/R_t;%最佳叶尖速比
n_ref = omega_best*30/pi;
% if t>2
% v_w = 12;
% end
%%%%%%%%%%%%%%%%%%%%%%%%%% Algebra %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% algebra
% syms omega_m i_sq i_sd x1 x2;

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
lambda = (R_t*omega_m)/(v_w);
lambda_i = 1/(1/(lambda+0.08*beta)-0.035/(beta^3+1));
Cp = 0.51763*(116/lambda_i-0.4*beta-5)*exp(-21/lambda_i)+0.006795*lambda;
P_m = u1*u2*(R_t^2)*(v_w^3)*Cp;
T_m = (-1)*(P_m)/omega_m;



% Tem_cmd = ((omega_r^3 * (1/1.1^3))/omega_r) * (-1);
% 
i_sdref = 0;

%dy1 = n_ref - n;
% dx1 = i_sdref - i_sd;
% dx2 = i_sqref - i_sq;
i_sqref = K_p3*(n_ref - n) + K_i3*y1;
% i_sqref = -607;
u_sqref = K_p2*(i_sqref - i_sq)+K_i2*x2+omega_e*Psi_f+omega_e*L_sd*i_sd;
u_sdref = K_p1*(i_sdref - i_sd)+K_i1*x1-omega_e*L_sq*i_sq;

Q_ref = 0;
P_s = 1.5*(u_sd * i_sd + u_sq * i_sq);%(2-17)
Q_s = 1.5*(u_sq * i_sd - u_sd * i_sq);%(2-18)
P_e = Te * omega_m;

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
    ];
end
