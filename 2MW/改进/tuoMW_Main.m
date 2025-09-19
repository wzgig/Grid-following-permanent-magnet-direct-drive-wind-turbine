clear;
clc;

% 运行参数文件
run("tuoMWRef_and_Para.m");

% 建立稳态方程组
eqn1 = (Te - T_m - D_m*omega_m) / J == 0;
eqn2 = u_sq - R_s*i_sq - omega_e*Psi_sd == 0;
eqn3 = u_sd - R_s*i_sd + omega_e*Psi_sq == 0;
eqn4 = i_sdref - i_sd == 0;
eqn5 = i_sqref - i_sq == 0;
eqn6 = n_ref - n == 0;
eqn7 = (u_sdref - u_sd)/T_d == 0;
eqn8 = (u_sqref - u_sq)/T_d == 0;
eqn9 = v_g_q == 0;
eqn10 = Kp_PLL*v_g_q+Ki_PLL*PLL_int+w_g == 0;
eqn11 = (P_s - P_dc)/(C_dc*U_dc) == 0;
eqn12 = U_dc - U_dcref == 0;
eqn13 = (v_g_d - u_gd - R_g*i_gd + w_g*L_g*i_gq)/L_g == 0;
eqn14 = (v_g_q - u_gq - R_g*i_gq - w_g*L_g*i_gd)/L_g == 0;
eqn15 = i_gdref - i_gd == 0;
eqn16 = i_gqref - i_gq == 0;
eqn17 = i_gD - iD == 0;
eqn18 = i_gQ - iQ == 0;

% 求解稳态值
[omega_m, Psi_sq, Psi_sd, Id_stator_int, Iq_stator_int, Speed_int,...
    u_sd, u_sq, PLL_int, thetapll, U_dc, Udc_int, i_gd, i_gq, Id_grid_int, Iq_grid_int,...
    v_w, i_gqref]...
                      = vpasolve(eval([eqn1,eqn2,eqn3,eqn4,eqn5,eqn6, ...
                      eqn7,eqn8,eqn9,eqn10,eqn11,eqn12,eqn13,eqn14,eqn15,eqn16, ...
                      eqn17,eqn18]),...
                      [omega_m Psi_sq Psi_sd Id_stator_int Iq_stator_int Speed_int ...
                      u_sd u_sq PLL_int thetapll U_dc Udc_int i_gd i_gq Id_grid_int Iq_grid_int ...
                      v_w i_gqref],...
                      [-inf,inf;-inf,inf;-inf,inf;-inf,inf;-inf,inf;-inf,inf;-inf,inf;-inf,inf; ...
                      -inf,inf;-inf,inf;-inf,inf;-inf,inf;-inf,inf;-inf,inf;-inf,inf;-inf,inf; ...
                      -inf,inf;-inf,inf]);

inistate = [omega_m Psi_sq Psi_sd Id_stator_int Iq_stator_int Speed_int...
    u_sd u_sq PLL_int thetapll U_dc Udc_int i_gd i_gq Id_grid_int Iq_grid_int]';
state_size = size(inistate,1);

v_w = double(v_w);
i_gqref = double(i_gqref);
disp(['风速: ', num2str(v_w)]);
disp(['q轴电流参考: ', num2str(i_gqref)]);

% **解决方案1: 使用刚性求解器**
dt = 0.0001;  % 增大时间步长
tspan = 0:dt:10;

% 使用ode15s刚性求解器，并调整容差
options = odeset('RelTol', 1e-6, 'AbsTol', 1e-9*ones(1,state_size), ...
                 'MaxStep', 0.001, 'InitialStep', 1e-6, ...
                 'BDF', 'on', 'Stats', 'on');

fprintf('开始求解微分方程...\n');
try
    [t, x] = ode15s(@(t,x) tuoMW(t,x,v_w, i_gqref), tspan, double(inistate), options);
    fprintf('求解成功完成！\n');
catch ME
    fprintf('求解失败: %s\n', ME.message);
    
    % **解决方案2: 如果ode15s失败，尝试ode23t**
    fprintf('尝试使用ode23t求解器...\n');
    options_23t = odeset('RelTol', 1e-5, 'AbsTol', 1e-8*ones(1,state_size), ...
                         'MaxStep', 0.01, 'InitialStep', 1e-5);
    try
        [t, x] = ode23t(@(t,x) tuoMW(t,x,v_w, i_gqref), tspan, double(inistate), options_23t);
        fprintf('ode23t求解成功！\n');
    catch ME2
        fprintf('ode23t也失败: %s\n', ME2.message);
        return;
    end
end

% 状态变量提取
omega_m_sim = x(:,1);
Psi_sq_sim = x(:,2);
Psi_sd_sim = x(:,3);
Id_stator_int_sim = x(:,4);
Iq_stator_int_sim = x(:,5);
Speed_int_sim = x(:,6);
u_sd_sim = x(:,7);
u_sq_sim = x(:,8);
PLL_int_sim = x(:,9);
thetapll_sim = x(:,10);
U_dc_sim = x(:,11);
Udc_int_sim = x(:,12);
i_gd_sim = x(:,13);
i_gq_sim = x(:,14);
Id_grid_int_sim = x(:,15);
Iq_grid_int_sim = x(:,16);

% 重新计算需要的参数用于绘图
run("tuoMWRef_and_Para.m");

% 计算功率变量
P_m_sim = zeros(length(t), 1);
P_s_sim = zeros(length(t), 1);
P_dc_sim = zeros(length(t), 1);
P_g_sim = zeros(length(t), 1);
Q_g_sim = zeros(length(t), 1);

for i = 1:length(t)
    % 使用当前时刻的状态变量
    omega_m_val = omega_m_sim(i);
    u_sd_val = u_sd_sim(i);
    u_sq_val = u_sq_sim(i);
    i_gd_val = i_gd_sim(i);
    i_gq_val = i_gq_sim(i);
    thetapll_val = thetapll_sim(i);
    
    % 计算电流
    i_sd_val = (Psi_sd_sim(i) - Psi_f)/L_sd;
    i_sq_val = Psi_sq_sim(i)/L_sq;
    
    % 计算风机机械功率
    lambda_val = (R_t*omega_m_val)/(v_w);
    lambda_i_val = 1/(1/(lambda_val+0.08*beta)-0.035/(beta^3+1));
    Cp_val = 0.51763*(116/lambda_i_val-0.4*beta-5)*exp(-21/lambda_i_val)+0.006795*lambda_val;
    P_m_sim(i) = 0.5*pi*1.225*(R_t^2)*(v_w^3)*Cp_val;
    
    % 计算定子功率
    P_s_sim(i) = sqrt(3)*(u_sd_val * i_sd_val + u_sq_val * i_sq_val);
    
    % 计算网侧功率
    theta_g_val = thetapll_val;
    Tg_DQdq = [cos(theta_g_val) -sin(theta_g_val); sin(theta_g_val) cos(theta_g_val)];
    vg_dq_real = Vbase * Tg_DQdq * [v_a; v_b];
    v_g_d_val = vg_dq_real(1);
    v_g_q_val = -vg_dq_real(2);
    
    i_gdref_val = Kp_Udc*(U_dc_sim(i) - U_dcref) + Ki_Udc*Udc_int_sim(i);
    u_gd_val = - Kp_Id_grid*(i_gdref_val - i_gd_val) - Ki_Id_grid*Id_grid_int_sim(i) - R_g*i_gd_val + v_g_d_val+i_gq_val*w_g*L_g;
    u_gq_val = - Kp_Iq_grid*(i_gqref - i_gq_val) - Ki_Iq_grid*Iq_grid_int_sim(i) - R_g*i_gq_val + v_g_q_val-i_gd_val*w_g*L_g;
    
    P_dc_sim(i) = sqrt(3)*(u_gd_val* i_gd_val + u_gq_val * i_gq_val);
    P_g_sim(i) = sqrt(3)*(v_g_d_val * i_gd_val + v_g_q_val * i_gq_val);
    Q_g_sim(i) = sqrt(3)*(v_g_q_val * i_gd_val - v_g_d_val * i_gq_val);
end

% 绘图
figure(1);
plot(t, P_m_sim/Sbase, 'DisplayName', 'P_m');
hold on;
plot(t, -P_s_sim/Sbase, 'DisplayName', '-P_s');
hold on;
plot(t, -P_dc_sim/Sbase, 'DisplayName', '-P_{dc}');
hold on;
legend('Location', 'best', 'Interpreter', 'tex');
xlabel('时间 (s)');
ylabel('功率 (p.u.)');
title('功率响应');
grid on;

figure(2);
plot(t, P_g_sim/Sbase, 'DisplayName', 'P_g');
hold on;
plot(t, Q_g_sim/Sbase, 'DisplayName', 'Q_g');
hold on;
legend('Location', 'best', 'Interpreter', 'tex');
xlabel('时间 (s)');
ylabel('功率 (p.u.)');
title('网侧功率响应');
grid on;

% 显示最终稳态值
fprintf('\n最终稳态值:\n');
fprintf('机械转速: %.4f rad/s\n', omega_m_sim(end));
fprintf('直流电压: %.4f V\n', U_dc_sim(end));
fprintf('网侧有功功率: %.4f p.u.\n', P_g_sim(end)/Sbase);
fprintf('网侧无功功率: %.4f p.u.\n', Q_g_sim(end)/Sbase);