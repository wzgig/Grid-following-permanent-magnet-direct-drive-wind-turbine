%%%%%%%%%%%%%%%%%%%%%%%%%%  Paras  %%%%%%%%%%%%%%%%%%%%%%%%

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

%%%%%%%%%%%%%%%%%%%%%%%%%% Ref    %%%%%%%%%%%%%%%%%%%%%%%%

% wm_ref = 0.8*w_g;
% Ird_ref = 0;
syms wm_ref
syms Ird_ref

Vdc_ref = 1200;
Igq_ref = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%% Algebra %%%%%%%%%%%%%%%%%%%%%%%%

syms w_m Psi_ds Psi_dr Psi_qs Psi_qr Turbine_w shaft_w...
     Phi_ra Phi_rc Phi_rb...
    V_dc i_gd i_gq Phi_ga Phi_gc Phi_gb

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

P_ssc = sqrt(3)*(u_sd*I_sd + u_sq*I_sq);
P_rsc = sqrt(3)*(u_rd*I_rd + u_rq*I_rq);
P_gsc = sqrt(3)*(u_gd*i_gd + u_gq*i_gq);
Ptot = sqrt(3)*(u_sq*I_sq + u_sd*I_sd) - sqrt(3)*(u_gq*i_gq + u_gd*i_gd);
Qtot = sqrt(3)*(u_sq*I_sd - u_sd*I_sq) - sqrt(3)*(u_gq*i_gd - u_gd*i_gq);
Q_ssc = sqrt(3)*(u_sq*I_sd - u_sd*I_sq);
Q_gsc = sqrt(3)*(u_gq*i_gd - u_gd*i_gq);

Vs_abs = sqrt(u_sd^2 + u_sq^2);
% i_d = (I_sd - i_gd)/Ibase;
% i_q = (I_sq - i_gq)/Ibase;
