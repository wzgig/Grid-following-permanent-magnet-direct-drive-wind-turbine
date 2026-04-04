function dotx = DFIG(t,x)

%%%%%%%%%%%%%%%%%%%%%%%%%% Rename %%%%%%%%%%%%%%%%%%%%%%%%

w_m = x(1);
Psi_ds = x(2); Psi_qs = x(3); Psi_dr = x(4); Psi_qr = x(5);
Turbine_w  = x(6);
shaft_w    = x(7);
%%%%%%%%%%%%%%%%%%%%%%%%%%  Paras  %%%%%%%%%%%%%%%%%%%%%%%%

f   = 60;
p   = 2;
w_g = 2*pi*f;
u   = 0.34;   % stator/rotor

Sbase = 2e6;
Vbase = 690;
Ibase = Sbase/(sqrt(3)*Vbase);
wbase = w_g;
Zbase = Vbase^2/Sbase;
Lbase = Zbase/w_g;
Tbase = Sbase/(w_g/p);

% Induction Motor - 9
Hg  = 0.8;
J   = (2*Hg*p^2*Sbase)/w_g^2;
Rs  = 0.011*Zbase;
Rr  = 1.1*Rs;
Lm  = 4*Lbase;
Ls  = 1.01*Lm;
Lr  = 1.005*Lm;
% Wind Turbine
rho = 1.225;  % 空气密度，单位：kg/m^3
R = 44;    % 风轮半径，单位：米
A = pi * R^2; % 风轮扫掠面积，单位：m^2
beta = 0;     % 桨距角，单位：度（固定值
N = 100;      % Gearbox Ratio
Ksh = 1.2543e+04;
D_mutual = 2.5e3;
H_WT = 0.4;
wind = 11.2;


%%%%%%%%%%%%%%%%%%%%%%%%%% Ref    %%%%%%%%%%%%%%%%%%%%%%%%

%     T_wt = -1.1893;
%     wind = 11.1
if t<2
    u_sd = 0; u_sq = 690;
    u_rd = -20;
    u_rq = 200;
end
if t>2
    u_sd = 0; u_sq = 690;
    u_rd = -20;
    u_rq = 500;
end
%%%%%%%%%%%%%%%%%%%%%%%%%% Algebra %%%%%%%%%%%%%%%%%%%%%%%%
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
T_e  = 1.5*p*Lm*(I_rd*I_sq - I_sd*I_rq);
Ptot = 1.5*(u_sq*I_sq + u_sd*I_sd) + 1.5*(u_rq*I_rq + u_rd*I_rd);
    

%%%%%%%%%%%%%%%%%%%%%%% Differential %%%%%%%%%%%%%%%%%%%%%

    dotx = [
        (p/J)*(T_e - T_shaft);         % state1
        u_sd - Rs*I_sd + w_g*Psi_qs;   % state2
        u_sq - Rs*I_sq - w_g*Psi_ds;   % state3
        u_rd - Rr*I_rd + w_r*Psi_qr;   % state4
        u_rq - Rr*I_rq - w_r*Psi_dr;   % state5
        (1/(2*H_WT))*(T_wt - T_shaft); % state20
        Turbine_w - w_m/p;             % state21
    ];
end


