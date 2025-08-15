function dotx = DFIG(t,x)

%%%%%%%%%%%%%%%%%%%%%%%%%% Rename %%%%%%%%%%%%%%%%%%%%%%%%

w_m = x(1);
Psi_ds = x(2); Psi_qs = x(3); Psi_dr = x(4); Psi_qr = x(5);
Turbine_w  = x(6);
shaft_w    = x(7);
Phi_ra   = x(8);
Phi_rc  = x(9);
Phi_rb    = x(10);
V_dc = x(11);
i_gd = x(12);
i_gq = x(13);
Phi_ga   = x(14);
Phi_gc  = x(15);
Phi_gb    = x(16);
Phipllr    = x(17);
theta_pllr = x(18);
Phipllg    = x(19);
theta_pllg = x(20);

%%%%%%%%%%%%%%%%%%%%%%%%%%  Paras  %%%%%%%%%%%%%%%%%%%%%%%%

P=-0.8;
Q=-0.3;
V=0.9998;
xi=0.6286;

f   = 60;
p   = 2;
w_g = 2*pi*f;

Sbase = 2e6;
Vbase = 690;
Ibase = Sbase/(sqrt(3)*Vbase);
wbase = w_g;
Zbase = Vbase^2/Sbase;
Lbase = Zbase/w_g;
Tbase = Sbase/(w_g/p);

S_D0 = P/V;
S_Q0 = Q/V;
S_DQ0 = S_D0 + 1j*S_Q0;
i_abs = abs(S_DQ0);
ui_argdiff = angle(S_DQ0);
i_arg = xi - ui_argdiff;
i_DQ = i_abs * exp(1i * i_arg);
i_q  = imag(i_DQ); % the global current should be output
i_d = real(i_DQ);
v_DQ = V * exp(1i * xi);
v_d = real(v_DQ);  % the global voltage should be input
v_q = imag(v_DQ);

% Induction Motor - 9
Hg  = 0.8;
J   = (2*Hg*p^2*Sbase)/w_g^2;
Rs  = 0.011*Zbase;
Rr  = 1.1*Rs;
Lm  = 4*Lbase;
Ls  = 1.01*Lm;
Lr  = 1.01*Lm;

% Wind Turbine
rho = 1.225; % 空气密度，单位：kg/m^3
R = 44; % 风轮半径，单位：米
A = pi * R^2; % 风轮扫掠面积，单位：m^2
beta = 0; % 桨距角，单位：度（固定值）
N = 100; % 齿轮箱传动比
Ksh = 1.2543e+04; % 轴系刚度
D_mutual = 2.5e3; % 相互阻尼
H_WT = 0.4; % 风力发电机转动惯量常数
wind = 12; % 风速，单位：m/s

% RSC - 6
sigma = 1- Lm^2/(Ls*Lr); 
tau_i = (sigma*Lr)/Rr;
tau_n = Hg/16;
wni = 1/tau_i;
wnn = 1/tau_n;
KPn = (2*wnn*J)/p;
KIn = ((wnn^2)*J)/p;
KPq = (2*wni*sigma*Lr)-Rr;
KIq = (wni^2)*Lr*sigma;
KPd = (2*wni*sigma*Lr)-Rr;
KId = (wni^2)*Lr*sigma;

w_n_limit = 1.5*Tbase;
i_dq_limit = 1200/sqrt(3);

% AC Filter
Lg  = 483e-6;
Rg  = 20e-6;

% DC Capacitor
C_bus = 80e-3;

% GSC - 6
K_pg = 1/(1.5*Vbase); 
K_qg = -K_pg;
KP_V = -1000;
KI_V = -300000;
KPdg = 2*(f*2*pi)*Lg - Rg;
KIdg = Lg*(f*2*pi)^2;
KPqg = KPdg;
KIqg = KIdg;

% PLL - 4
rPLLp = 220;
rPLLi = 4500;
gPLLp = 250;
gPLLi = 3200;

%%%%%%%%%%%%%%%%%%%%%%%%%% Ref    %%%%%%%%%%%%%%%%%%%%%%%%

usd_ref = 0;
ugq_ref = 0;
% Ird_ref = 1313.0844995607461804041666784972;
Ird_ref = 1234.3489266975635554455612065438;
% wm_ref= 548.55317152444422923268545408028;
wm_ref = 548.62242452840193032996249766678;
Igq_ref = 0;
Vdc_ref = 1200;

% if t>2
% Ird_ref = 1313.0844995607461804041666784972;
% end
if t> 1
wind= 16;
end
% if t> 1.2
% wind= 11;
% end

% Ls = 0.00265104683633706;

% if t>1
%     New = [cos(0.015) -sin(0.015); sin(0.015) cos(0.015)]*[0.9*v_d;0.9*v_q];
%     v_d = New(1);
%     v_q = New(2);
% end
% if t>1.2
%     New = [cos(-0.015) -sin(-0.015); sin(-0.015) cos(-0.015)]*[v_d;v_q];
%         v_d = New(1);
%     v_q = New(2);
% %      Ls = 0.00265104683633706;
% end

%%%%%%%%%%%%%%%%%%%%%%%%%% Algebra %%%%%%%%%%%%%%%%%%%%%%%%
% PLL
theta_r = theta_pllr - pi;
theta_g = theta_pllg;
Tr_DQdq = [cos(theta_r) -sin(theta_r); sin(theta_r) cos(theta_r)];
vr_dq_real = Vbase*Tr_DQdq*[v_d;v_q]; % RSC PLL
Tg_DQdq = [cos(theta_g) -sin(theta_g); sin(theta_g) cos(theta_g)];
vg_dq_real = Vbase*Tg_DQdq*[v_d;v_q]; % GSC PLL
u_sd = vr_dq_real(1);
u_sq = vr_dq_real(2);
u_gd = vg_dq_real(1);
u_gq = -vg_dq_real(2);

Pitch = 0;
lambda = Turbine_w/N * R / wind;
lambda_i = 1/((1/(lambda-0.02*Pitch)+(0.003/(Pitch^3+1))));
Cp = 0.73*(151/lambda_i-0.58*Pitch-0.002*Pitch^2.14-13.2)*(exp(-18.4/lambda_i));
P_wt = 0.5*rho*pi*(R)^2*(wind)^3*Cp;
T_wt = -P_wt/(Turbine_w/N)/N;
T_shaft = Ksh*shaft_w + D_mutual*(Turbine_w - w_m/p);

Psi2I = inv([Lm Ls 0 0;Lr Lm 0 0;0 0 Lm Ls;0 0 Lr Lm]);
I = Psi2I*[Psi_ds;Psi_dr;Psi_qs;Psi_qr];
I_rd = I(1); I_sd = I(2); I_rq = I(3); I_sq = I(4);
w_r = w_g - w_m;
T_e  = sqrt(3)*p*Lm*(I_rd*I_sq - I_sd*I_rq);

alpha1 = Lm/Ls; alpha2 = Lr - Lm*alpha1;
PHI = sqrt(Psi_ds^2 + Psi_qs^2);
compensatord = w_r*I_rq*(-alpha2);
compensatorq = w_r*PHI*alpha1 + w_r*I_rd*alpha2;
u_rd = (KPd*(Ird_ref - I_rd) + Phi_ra + compensatord);
Irq_ref = (KPn*(wm_ref - w_m) + Phi_rc)/(-1.5*p*alpha1*PHI);
u_rq = (KPq*(Irq_ref - I_rq) + Phi_rb + compensatorq);

compensatorgd = u_gd - w_g*Lg*i_gq;
compensatorgq = u_gq + w_g*Lg*i_gd;
vrefq = (KPqg*(Igq_ref - i_gq)+ Phi_ga) + compensatorgq;
Igd_ref = (KP_V*(Vdc_ref - V_dc) + Phi_gc)*K_pg;
vrefd = (KPdg*(Igd_ref - i_gd)+ Phi_gb) - compensatorgd;

P_rsc = sqrt(3)*(u_rd*I_rd + u_rq*I_rq);
P_gsc = sqrt(3)*(u_gd*i_gd + u_gq*i_gq);

%%%%%%%%%%%%%%%%%%%%%%% Differential %%%%%%%%%%%%%%%%%%%%%

    dotx = [
        (p/J)*(T_e - T_shaft);         % state1
        u_sd - Rs*I_sd + w_g*Psi_qs;   % state2
        u_sq - Rs*I_sq - w_g*Psi_ds;   % state3
        u_rd - Rr*I_rd + w_r*Psi_qr;   % state4
        u_rq - Rr*I_rq - w_r*Psi_dr;   % state5
        (1/(2*H_WT))*(T_wt - T_shaft); % state20
        Turbine_w - w_m/p;             % state21
        KId*(Ird_ref - I_rd);
        KIn*(wm_ref - w_m);
        KIq*(Irq_ref - I_rq);
        (1/(C_bus*V_dc))*(-P_gsc-P_rsc);
        (1/Lg)*(vrefd - u_gd - Rg*i_gd + Lg*w_g*i_gq);
        (1/Lg)*(vrefq - u_gq - Rg*i_gq - Lg*w_g*i_gd);
        KIqg*(Igq_ref - i_gq);
        KI_V*(Vdc_ref - V_dc);
        KIdg*(Igd_ref - i_gd);
        rPLLi*(u_sd/Vbase - usd_ref);               % state15
        Phipllr + rPLLp*(u_sd/Vbase-usd_ref);    % state16
        gPLLi*(u_gq/Vbase - ugq_ref);                % state17
        Phipllg + gPLLp*(u_gq/Vbase - ugq_ref);  % state18
    ];
end
