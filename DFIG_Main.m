
clear;
clc;

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
i_Q  = imag(i_DQ); % the global current should be output
i_D = real(i_DQ);
v_DQ = V * exp(1i * xi);
v_d = real(v_DQ);  % the global voltage should be input
v_q = imag(v_DQ);

usd_ref = 0;
ugq_ref = 0;

% PLL - 4
rPLLp = 220;
rPLLi = 4500;
gPLLp = 250;
gPLLi = 3200;

syms Phipllr theta_pllr Phipllg theta_pllg
theta_r = theta_pllr - pi;
theta_g = theta_pllg;
Tr_DQdq = [cos(theta_r) -sin(theta_r); sin(theta_r) cos(theta_r)];
vr_dq_real = Vbase*Tr_DQdq*[v_d;v_q]; % RSC PLL
Tg_DQdq = [cos(theta_g) -sin(theta_g); sin(theta_g) cos(theta_g)];
vg_dq_real = Vbase*Tg_DQdq*[v_d;v_q]; % GSC PLL
Tr_dqDQ = [cos(-theta_r) -sin(-theta_r); sin(-theta_r) cos(-theta_r)];
Tg_dqDQ = [cos(-theta_g) -sin(-theta_g); sin(-theta_g) cos(-theta_g)];
u_sd = vr_dq_real(1);
u_sq = vr_dq_real(2);
u_gd = vg_dq_real(1);
u_gq = -vg_dq_real(2);
eqnd = rPLLi*(u_sd/Vbase - usd_ref) == 0;          % state14
eqnq = gPLLi*(u_gq/Vbase - ugq_ref) == 0;          % state16
[pllr pllg] = solve(eval([eqnd, eqnq]),[theta_pllr theta_pllg]);
theta_pllr = vpa(pllr(2));
theta_pllg = vpa(pllg(2));
Phipllr = 0; % state15
Phipllg = 0; % state17

u_sd = round(eval(u_sd));
u_sq = round(eval(u_sq));
u_gd = round(eval(u_gd));
u_gq = round(eval(u_gq));

run("Ref_and_Para.m");

is_dq_pu = (1/Ibase)*Tr_dqDQ*[I_sd;I_sq];
ig_dq_pu = (1/Ibase)*Tg_dqDQ*[i_gd;i_gq];

i_d = is_dq_pu(1)-ig_dq_pu(1);
i_q = is_dq_pu(2)-ig_dq_pu(2);

eqn1 = (p/J)*(T_e - T_shaft) == 0; % state1
eqn2 = (1/(2*H_WT))*(T_wt - T_shaft) == 0;% state20
eqn3 = Turbine_w - w_m/p;  % state21
eqn4 = u_sd - Rs*I_sd + w_g*Psi_qs == 0; % state2
eqn5 = u_sq - Rs*I_sq - w_g*Psi_ds == 0; % state3
eqn6 = u_rd - Rr*I_rd + w_r*Psi_qr == 0; % state4
eqn7 = u_rq - Rr*I_rq - w_r*Psi_dr == 0; % state5
eqn8 = KId*(Ird_ref - I_rd) == 0;
eqn9 = KIn*(wm_ref - w_m) == 0;
eqn10 = KIq*(Irq_ref - I_rq) == 0;
eqn11 = (1/(C_bus*V_dc))*(-P_gsc-P_rsc) == 0;
eqn12 = (1/Lg)*(vrefd - u_gd - Rg*i_gd + Lg*w_g*i_gq) == 0;
eqn13 = (1/Lg)*(vrefq - u_gq - Rg*i_gq - Lg*w_g*i_gd) == 0;
eqn14 = KIqg*(Igq_ref - i_gq) == 0;
eqn15 = KI_V*(Vdc_ref - V_dc) == 0;
eqn16 = KIdg*(Igd_ref - i_gd) == 0;
% eqn17 = Ptot/Sbase == P;
% eqn18 = Qtot/Sbase == Q;
eqn17 = i_d == i_D;
eqn18 = i_q == i_Q;

[Turbine_w w_m Psi_ds Psi_qs Psi_dr Psi_qr shaft_w...
    Phi_ra Phi_rc Phi_rb...
    V_dc i_gd i_gq Phi_ga Phi_gc Phi_gb wm_ref Ird_ref]...
                      = vpasolve(eval([eqn1, eqn2, eqn3, eqn4, eqn5, eqn6, eqn7,...
                      eqn8, eqn9, eqn10, eqn11, eqn12, eqn13, eqn14,...
                      eqn15, eqn16, eqn17, eqn18]),...
                      [Turbine_w w_m Psi_ds Psi_qs Psi_dr Psi_qr shaft_w...
                         Phi_ra Phi_rc Phi_rb...
                        V_dc i_gd i_gq Phi_ga Phi_gc Phi_gb wm_ref Ird_ref],...
                      [40,300;-2*w_g,2*w_g;-inf,inf;-inf,inf;-inf,inf;-inf,inf;-inf,inf;...
                      -inf,inf;-inf,inf;-inf,inf;-inf,inf;-inf,inf;-inf,inf;...
                      -inf,inf;-inf,inf;-inf,inf;-2*w_g,2*w_g;-inf,inf]);

inistate = [w_m Psi_ds Psi_qs Psi_dr Psi_qr Turbine_w shaft_w...
    Phi_ra Phi_rc Phi_rb...
    V_dc i_gd i_gq Phi_ga Phi_gc Phi_gb Phipllr theta_pllr Phipllg theta_pllg]';

state_size = size(inistate,1);

dt = 0.00001;
tspan=dt:dt:10;
options = odeset('RelTol',1e-12,'AbsTol',...
    1e-12*ones(1,state_size));
[t,x] = ode45(@(t,x) DFIG(t,x),tspan,double(inistate),options);

run("Ref_and_Para.m");

%%%%%%%%%% states

w_m    = x(:,1);
Psi_ds = x(:,2);
Psi_qs = x(:,3);
Psi_dr = x(:,4);
Psi_qr = x(:,5);
Turbine_w = x(:,6);
shaft_w = x(:,7);
Phi_ra = x(:,8);
Phi_rc  = x(:,9);
Phi_rb = x(:,10);
V_dc  = x(:,11);
i_gd = x(:,12);
i_gq = x(:,13);
Phi_ga = x(:,14);
Phi_gc = x(:,15);
Phi_gb= x(:,16);
Phipllr = x(:,17);
theta_pllr = x(:,18);
Phipllg = x(:,19);
theta_pllg= x(:,20);


% plot(Phipllr)
% plot(theta_pllr)
% plot(Phipllg)
% plot(theta_pllg)
% plot(tspan,V_dc)
% plot(tspan,i_gd)
% plot(tspan,i_gq)
% plot(tspan,eval(Igd_ref))
% 
% % plot(tspan,Phi_rc);
% % plot(tspan,Phi_rb);
% % plot(tspan,Phi_ra);
% plot(tspan,eval(u_rd));
% hold on;
% plot(tspan,eval(u_rq));
% 
% plot(tspan,eval(I_rd));
% hold on;
% % plot(tspan,eval(Ird_ref));
% plot(tspan,eval(I_rq));
% hold on;
% plot(tspan,eval(Irq_ref));
% figure;
% plot(tspan,w_m/w_g);
% plot(tspan,eval(T_e)/Tbase);
% plot(tspan,eval(T_wt)/Tbase);
% plot(tspan,eval(-P_rsc-P_gsc))
% plot(tspan,eval(P_gsc/Ptot))
% 
plot(tspan,eval(T_e)/Tbase);
hold on

Ptot = sqrt(3)*(eval(u_sq*I_sq) + eval(u_sd*I_sd)) - sqrt(3)*(eval(u_gq)*i_gq + eval(u_gd)*i_gd);
Qtot  = sqrt(3)*(eval(u_sq*I_sd) - eval(u_sd*I_sq)) - sqrt(3)*(eval(u_gq)*i_gd - eval(u_gd)*i_gq);
plot(tspan,Ptot/Sbase);
hold on;
plot(tspan,Qtot/Sbase);

% plot(tspan,eval(v_d*i_d+v_q*i_q));
% hold on;
% plot(tspan,eval(v_q*i_d-v_d*i_q));

% plot(tspan,eval(PHI));
% plot(tspan,eval(Q_ssc)/Sbase);
% plot(tspan,eval(Vs_abs));

% plot(tspan,eval(P_ssc)/Sbase);
% plot(tspan,eval((u_rq.*I_rq + u_rd.*I_rd)))
% plot(tspan,eval((u_rq.*I_rq + u_rd.*I_rd)/(u_sq.*I_sq + u_sd.*I_sd)));
% hold on
% plot(tspan,eval(w_r)/w_g);
% plot(tspan,w_m);
% plot(tspan,eval(I_rq*u)/Ibase);
% hold on
% plot(tspan,eval(I_rd*u)/Ibase);
% plot(tspan,eval(I_sd));
% plot(tspan,eval(I_sq));
% hold on
% plot(tspan,eval(I_rd));
% plot(tspan,eval(I_rq));
% plot(tspan,eval(u_rd));
% 
% plot(tspan,eval(u_rq));
plot(tspan,eval(Cp));
plot(tspan,eval(P_wt));
plot(tspan,w_m);