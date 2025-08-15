clear;
clc;

run("Copy_of_wangceRef_and_Para.m");
eqn1 = (Te - T_m - D_m*omega_m) / J == 0;
% eqn2 = u_sd - R_s*i_sd + omega_e*Psi_sq == 0;
% eqn3 = u_sq - R_s*i_sq - omega_e*Psi_sd == 0;
eqn2 = (u_sq - R_s*i_sq - omega_e*L_sd*i_sd - omega_e*Psi_f)/L_sq == 0;
eqn3 = (u_sd - R_s*i_sd + omega_e*L_sq*i_sq)/L_sd == 0;

eqn4 = i_sdref - i_sd == 0;
eqn5 = i_sqref - i_sq == 0;
eqn6 = n_ref - n == 0;
eqn7 = (u_sdref - u_sd)/T_d == 0;
eqn8 = (u_sqref - u_sq)/T_d == 0;
eqn9 = v_g_q == 0;
eqn10 = K_pPLL*v_g_q+K_iPLL*xpll+w_g == 0;
eqn11 = (P_s - P_dc)/(C_dc*U_dc) == 0;
eqn12 = U_dc - U_dcref == 0;
eqn13 = (v_g_d - u_gd + w_g*L_g*i_gq)/L_g == 0;
eqn14 = (v_g_q - u_gq - w_g*L_g*i_gd)/L_g == 0;
eqn15 = i_gdref - i_gd == 0;
eqn16 = i_gqref - i_gq == 0;



disp(eqn1);disp(eqn2);disp(eqn3);disp(eqn4);disp(eqn5);disp(eqn6);disp(eqn7);disp(eqn8);
disp(eqn9);disp(eqn10);disp(eqn11);disp(eqn12);disp(eqn13);disp(eqn14);disp(eqn15);disp(eqn16);

% syms omega_m i_sq i_sd x1 x2 y1 u_sd u_sq xpll thetapll U_dc z1 i_gd i_gq z2 z3;
[omega_m, i_sq, i_sd, x1, x2, y1, u_sd, u_sq, xpll, thetapll, U_dc, z1, i_gd, i_gq, z2, z3]...
                      = vpasolve(eval([eqn1,eqn2,eqn3,eqn4,eqn5,eqn6,eqn7,eqn8,eqn9,eqn10,eqn11,eqn12,eqn13,eqn14,eqn15,eqn16]),...
                      [omega_m i_sq i_sd x1 x2 y1 u_sd u_sq xpll thetapll U_dc z1 i_gd i_gq z2 z3],...
                      [0,20;-inf,inf;-inf,inf;-inf,inf;-inf,inf;-inf,inf;-inf,inf;-inf,inf;-inf,inf;-inf,inf;-inf,inf;-inf,inf;-inf,inf;-inf,inf;-inf,inf;-inf,inf]);

inistate = [omega_m i_sq i_sd x1 x2 y1 u_sd u_sq xpll thetapll U_dc z1 i_gd i_gq z2 z3]'
state_size = size(inistate,1);

dt = 0.00001;
tspan=dt:dt:10;
options = odeset('RelTol',1e-12,'AbsTol',...
    1e-12*ones(1,state_size));
[t,x] = ode45(@(t,x) Copy_of_wangce(t,x),tspan,double(inistate),options);

run("Copy_of_wangceRef_and_Para.m");

%%%%%%%%%% states

omega_m = x(:,1);
i_sq       = x(:,2);
i_sd       = x(:,3);
x1       = x(:,4);
x2       = x(:,5);
y1       = x(:,6);
u_sd       = x(:,7);
u_sq       = x(:,8);
xpll       = x(:,9);
thetapll       = x(:,10);
U_dc       = x(:,11);
z1       = x(:,12);
i_gd       = x(:,13);
i_gq       = x(:,14);
z2       = x(:,15);
z3       = x(:,16);

plot(tspan, i_sd, 'DisplayName', 'i_sd');
hold on;
plot(tspan, i_sq, 'DisplayName', 'i_sq');
hold on;
plot(tspan, u_sd, 'DisplayName', 'u_sd');
hold on;
plot(tspan, u_sq, 'DisplayName', 'u_sq');
hold on;
plot(tspan, eval(Cp), 'DisplayName', 'Cp');
hold on;
plot(tspan, eval(T_m), 'DisplayName', 'T_m');
hold on;
% plot(tspan, eval(Tem), 'DisplayName', 'Tem');
% hold on;
plot(tspan, eval(Te), 'DisplayName', 'Te');
hold on;
plot(tspan, eval(omega_e), 'DisplayName', 'omega_e');
hold on;
% plot(tspan, eval(Tem_cmd), 'DisplayName', 'Tem_cmd');
% hold on;
plot(tspan, omega_m, 'DisplayName', 'omega_m');
hold on;
plot(tspan, eval(P_m), 'DisplayName', 'P_m');
hold on;
plot(tspan, eval(P_s), 'DisplayName', 'P_s');
hold on;
plot(tspan, eval(Q_s), 'DisplayName', 'Q_s');
hold on;
plot(tspan, eval(P_e), 'DisplayName', 'P_e');
hold on;
plot(tspan, eval(n), 'DisplayName', 'n');
hold on;
plot(tspan, eval(P_g), 'DisplayName', 'P_g');
hold on;
plot(tspan, eval(Q_g), 'DisplayName', 'Q_g');
hold on;
plot(tspan, eval(u_gd), 'DisplayName', 'u_gd');
hold on;
plot(tspan, eval(u_gq), 'DisplayName', 'u_gq');
hold on;
plot(tspan, i_gd, 'DisplayName', 'i_gd');
hold on;
plot(tspan, i_gq, 'DisplayName', 'i_gq');
hold on;
plot(tspan, eval(P_dc), 'DisplayName', 'P_dc');
hold on;