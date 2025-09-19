function system_diagnostic()
%% 直驱风机系统数值稳定性诊断脚本
% 分析不同xi值对系统动态特性的影响

clear; clc;
fprintf('=== 直驱风机系统数值稳定性诊断 ===\n\n');

% 定义xi测试范围
xi_test = linspace(-pi/2, -pi/6, 20);
stability_results = zeros(length(xi_test), 4); % [xi, 特征值实部最大值, 条件数, 稳态求解成功标志]

for i = 1:length(xi_test)
    xi_current = xi_test(i);
    fprintf('测试 xi = %.4f rad (%.1f°)\n', xi_current, xi_current*180/pi);
    
    try
        % 设置当前xi值
        stability_results(i, 1) = xi_current;
        
        % 计算稳态工作点
        [success, inistate, v_w_val, i_gqref_val] = compute_steady_state(xi_current);
        stability_results(i, 4) = success;
        
        if success
            % 线性化分析
            [max_real_eigenval, condition_num] = linearization_analysis(inistate, v_w_val, i_gqref_val);
            stability_results(i, 2) = max_real_eigenval;
            stability_results(i, 3) = condition_num;
            
            fprintf('  稳态求解: 成功\n');
            fprintf('  最大特征值实部: %.2e\n', max_real_eigenval);
            fprintf('  雅可比矩阵条件数: %.2e\n', condition_num);
            
            % 判断数值稳定性
            if max_real_eigenval > -1e-3
                fprintf('  ⚠️  系统可能数值不稳定 (特征值实部接近0)\n');
            end
            if condition_num > 1e10
                fprintf('  ⚠️  系统刚性严重 (条件数过大)\n');
            end
        else
            fprintf('  ❌ 稳态求解失败\n');
        end
        
    catch ME
        fprintf('  ❌ 分析出错: %s\n', ME.message);
        stability_results(i, 4) = 0;
    end
    fprintf('\n');
end

%% 绘制诊断结果
plot_diagnostic_results(stability_results);

%% 给出建议
provide_recommendations(stability_results);

end

function [success, inistate, v_w_val, i_gqref_val] = compute_steady_state(xi_val)
%% 计算给定xi值的稳态工作点

try
    % 基本参数设置
    f = 60; Fnom = f; Np = 48; w_g = 2*pi*f;
    Sbase = 2e6; Vbase = 690; Ibase = Sbase/(sqrt(3)*Vbase);
    U_dcref = 1150; R_g = 0.05; L_g = 7e-3; C_dc = 1e-3;
    
    % 发电机参数
    Psi_f = 3.88889; L_sd = 1.8e-3; L_sq = 1.8e-3; R_s = 0.0026;
    beta = 0; J = 35000; D_m = 0.078; R_t = 36.6;
    
    % 保守的控制器参数
    Kp_PLL = 50; Ki_PLL = 500;  % 进一步降低PLL增益
    Kp_Id_stator = 0.5; Ki_Id_stator = 8;
    Kp_Iq_stator = 0.5; Ki_Iq_stator = 8;
    Kp_Speed = 60; Ki_Speed = 150;
    Kp_Udc = 0.5; Ki_Udc = 6;
    Kp_Id_grid = 0.5; Ki_Id_grid = 10;
    Kp_Iq_grid = 0.5; Ki_Iq_grid = 10;
    T_d = 1/3000; T_m = 1/1500;
    
    % 功率和电压设置
    P = -1; Q = 0; V = 1; xi = xi_val;
    S_D0 = P/V; S_Q0 = Q/V; S_ab0 = S_D0 + 1j*S_Q0;
    i_abs = abs(S_ab0); ui_argdiff = angle(S_ab0);
    i_arg = xi - ui_argdiff; i_ab = i_abs * exp(1i * i_arg);
    iQ = imag(i_ab); iD = real(i_ab);
    v_ab = V * exp(1i * xi); v_a = real(v_ab); v_b = imag(v_ab);
    
    % 符号变量定义
    syms omega_m Psi_sq Psi_sd Id_stator_int Iq_stator_int Speed_int...
        u_sd u_sq PLL_int thetapll U_dc Udc_int i_gd i_gq Id_grid_int Iq_grid_int;
    syms v_w i_gqref
    
    % 构建稳态方程
    omega_best = 8.1*v_w/R_t; n_ref = omega_best*30/pi;
    n = omega_m * 30/pi; omega_e = omega_m*Np;
    i_sd = (Psi_sd - Psi_f)/L_sd; i_sq = Psi_sq/L_sq;
    Te = sqrt(3)*Np*(Psi_f*i_sq);
    
    lambda = (R_t*omega_m)/(v_w);
    lambda_i = 1/(1/(lambda+0.08*beta)-0.035/(beta^3+1));
    Cp = 0.51763*(116/lambda_i-0.4*beta-5)*exp(-21/lambda_i)+0.006795*lambda;
    P_m = 0.5*pi*1.225*(R_t^2)*(v_w^3)*Cp;
    T_m = (-1)*(P_m)/omega_m;
    
    i_sdref = 0;
    i_sqref = Kp_Speed*(n_ref - n) + Ki_Speed*Speed_int;
    u_sqref = Kp_Iq_stator*(i_sqref - i_sq)+Ki_Iq_stator*Iq_stator_int+omega_e*Psi_f+omega_e*L_sd*i_sd;
    u_sdref = Kp_Id_stator*(i_sdref - i_sd)+Ki_Id_stator*Id_stator_int-omega_e*L_sq*i_sq;
    
    P_s = sqrt(3)*(u_sd * i_sd + u_sq * i_sq);
    
    theta_g = thetapll;
    Tg_DQdq = [cos(theta_g) -sin(theta_g); sin(theta_g) cos(theta_g)];
    vg_dq_real = Vbase * Tg_DQdq * [v_a; v_b];
    v_g_d = vg_dq_real(1); v_g_q = -vg_dq_real(2);
    
    Tg_dqDQ = [cos(-theta_g) -sin(-theta_g); sin(-theta_g) cos(-theta_g)];
    ig_dq = (1/Ibase)*Tg_dqDQ*[i_gd;i_gq];
    i_gD = ig_dq(1); i_gQ = ig_dq(2);
    
    i_gdref = Kp_Udc*(U_dc - U_dcref) + Ki_Udc*Udc_int;
    u_gd = - Kp_Id_grid*(i_gdref - i_gd) - Ki_Id_grid*Id_grid_int - R_g*i_gd + v_g_d+i_gq*w_g*L_g;
    u_gq = - Kp_Iq_grid*(i_gqref - i_gq) - Ki_Iq_grid*Iq_grid_int - R_g*i_gq + v_g_q-i_gd*w_g*L_g;
    P_dc = sqrt(3)*(u_gd* i_gd + u_gq * i_gq);
    
    % 稳态方程组
    equations = [
        (Te - T_m - D_m*omega_m) / J == 0;
        u_sq - R_s*i_sq - omega_e*Psi_sd == 0;
        u_sd - R_s*i_sd + omega_e*Psi_sq == 0;
        i_sdref - i_sd == 0;
        i_sqref - i_sq == 0;
        n_ref - n == 0;
        (u_sdref - u_sd)/T_d == 0;
        (u_sqref - u_sq)/T_d == 0;
        v_g_q == 0;
        Kp_PLL*v_g_q+Ki_PLL*PLL_int+w_g == 0;
        (P_s - P_dc)/(C_dc*U_dc) == 0;
        U_dc - U_dcref == 0;
        (v_g_d - u_gd - R_g*i_gd + w_g*L_g*i_gq)/L_g == 0;
        (v_g_q - u_gq - R_g*i_gq - w_g*L_g*i_gd)/L_g == 0;
        i_gdref - i_gd == 0;
        i_gqref - i_gq == 0;
        i_gD - iD == 0;
        i_gQ - iQ == 0;
    ];
    
    variables = [omega_m Psi_sq Psi_sd Id_stator_int Iq_stator_int Speed_int ...
                u_sd u_sq PLL_int thetapll U_dc Udc_int i_gd i_gq Id_grid_int Iq_grid_int ...
                v_w i_gqref];
    
    % 求解稳态
    solutions = vpasolve(equations, variables);
    
    if ~isempty(solutions.omega_m)
        inistate = [solutions.omega_m; solutions.Psi_sq; solutions.Psi_sd; 
                   solutions.Id_stator_int; solutions.Iq_stator_int; solutions.Speed_int;
                   solutions.u_sd; solutions.u_sq; solutions.PLL_int; solutions.thetapll;
                   solutions.U_dc; solutions.Udc_int; solutions.i_gd; solutions.i_gq;
                   solutions.Id_grid_int; solutions.Iq_grid_int];
        v_w_val = double(solutions.v_w);
        i_gqref_val = double(solutions.i_gqref);
        inistate = double(inistate);
        success = 1;
    else
        success = 0;
        inistate = []; v_w_val = []; i_gqref_val = [];
    end
    
catch
    success = 0;
    inistate = []; v_w_val = []; i_gqref_val = [];
end

end

function [max_real_eigenval, condition_num] = linearization_analysis(inistate, v_w_val, i_gqref_val)
%% 在稳态工作点进行线性化分析

try
    % 计算雅可比矩阵 (数值方法)
    epsilon = 1e-8;
    n_states = length(inistate);
    J_matrix = zeros(n_states, n_states);
    
    f0 = tuoMW_modified(0, inistate, v_w_val, i_gqref_val);
    
    for i = 1:n_states
        x_pert = inistate;
        x_pert(i) = x_pert(i) + epsilon;
        f_pert = tuoMW_modified(0, x_pert, v_w_val, i_gqref_val);
        J_matrix(:, i) = (f_pert - f0) / epsilon;
    end
    
    % 计算特征值
    eigenvals = eig(J_matrix);
    max_real_eigenval = max(real(eigenvals));
    
    % 计算条件数
    condition_num = cond(J_matrix);
    
catch
    max_real_eigenval = NaN;
    condition_num = NaN;
end

end

function dotx = tuoMW_modified(t, x, v_w, i_gqref)
%% 修改的状态方程，使用保守的控制器参数

% 基本参数
f = 60; Np = 48; w_g = 2*pi*f; Sbase = 2e6; Vbase = 690;
Ibase = Sbase/(sqrt(3)*Vbase); U_dcref = 1150;
R_g = 0.05; L_g = 7e-3; C_dc = 1e-3;

% 发电机参数
Psi_f = 3.88889; L_sd = 1.8e-3; L_sq = 1.8e-3; R_s = 0.0026;
beta = 0; J = 35000; D_m = 0.078; R_t = 36.6;

% 保守的控制器参数
Kp_PLL = 50; Ki_PLL = 500;
Kp_Id_stator = 0.5; Ki_Id_stator = 8;
Kp_Iq_stator = 0.5; Ki_Iq_stator = 8;
Kp_Speed = 60; Ki_Speed = 150;
Kp_Udc = 0.5; Ki_Udc = 6;
Kp_Id_grid = 0.5; Ki_Id_grid = 10;
Kp_Iq_grid = 0.5; Ki_Iq_grid = 10;
T_d = 1/3000; T_m = 1/1500;

% 功率参考
P = -1; Q = 0; V = 1; xi = -pi/6;
S_D0 = P/V; S_Q0 = Q/V; S_ab0 = S_D0 + 1j*S_Q0;
i_abs = abs(S_ab0); ui_argdiff = angle(S_ab0);
i_arg = xi - ui_argdiff; i_ab = i_abs * exp(1i * i_arg);
iQ = imag(i_ab); iD = real(i_ab);
v_ab = V * exp(1i * xi); v_a = real(v_ab); v_b = imag(v_ab);

% 状态变量
omega_m = x(1); Psi_sq = x(2); Psi_sd = x(3);
Id_stator_int = x(4); Iq_stator_int = x(5); Speed_int = x(6);
u_sd = x(7); u_sq = x(8); PLL_int = x(9); thetapll = x(10);
U_dc = x(11); Udc_int = x(12); i_gd = x(13); i_gq = x(14);
Id_grid_int = x(15); Iq_grid_int = x(16);

% 计算中间变量
omega_best = 8.1*v_w/R_t; n_ref = omega_best*30/pi; n = omega_m * 30/pi;
omega_e = omega_m*Np; i_sd = (Psi_sd - Psi_f)/L_sd; i_sq = Psi_sq/L_sq;
Te = sqrt(3)*Np*(Psi_f*i_sq);

lambda = (R_t*omega_m)/(v_w);
lambda_i = 1/(1/(lambda+0.08*beta)-0.035/(beta^3+1));
Cp = 0.51763*(116/lambda_i-0.4*beta-5)*exp(-21/lambda_i)+0.006795*lambda;
P_m = 0.5*pi*1.225*(R_t^2)*(v_w^3)*Cp;
T_m = (-1)*(P_m)/omega_m;

i_sdref = 0;
i_sqref = Kp_Speed*(n_ref - n) + Ki_Speed*Speed_int;
u_sqref = Kp_Iq_stator*(i_sqref - i_sq)+Ki_Iq_stator*Iq_stator_int+omega_e*Psi_f+omega_e*L_sd*i_sd;
u_sdref = Kp_Id_stator*(i_sdref - i_sd)+Ki_Id_stator*Id_stator_int-omega_e*L_sq*i_sq;

P_s = sqrt(3)*(u_sd * i_sd + u_sq * i_sq);

theta_g = thetapll;
Tg_DQdq = [cos(theta_g) -sin(theta_g); sin(theta_g) cos(theta_g)];
vg_dq_real = Vbase * Tg_DQdq * [v_a; v_b];
v_g_d = vg_dq_real(1); v_g_q = -vg_dq_real(2);

i_gdref = Kp_Udc*(U_dc - U_dcref) + Ki_Udc*Udc_int;
u_gd = - Kp_Id_grid*(i_gdref - i_gd) - Ki_Id_grid*Id_grid_int - R_g*i_gd + v_g_d+i_gq*w_g*L_g;
u_gq = - Kp_Iq_grid*(i_gqref - i_gq) - Ki_Iq_grid*Iq_grid_int - R_g*i_gq + v_g_q-i_gd*w_g*L_g;
P_dc = sqrt(3)*(u_gd* i_gd + u_gq * i_gq);

% 状态方程
dotx = [
    (Te - T_m - D_m*omega_m) / J;
    u_sq - R_s*i_sq - omega_e*Psi_sd;
    u_sd - R_s*i_sd + omega_e*Psi_sq;
    i_sdref - i_sd;
    i_sqref - i_sq;
    n_ref - n;
    (u_sdref - u_sd)/T_d;
    (u_sqref - u_sq)/T_d;
    v_g_q;
    Kp_PLL*v_g_q+Ki_PLL*PLL_int+w_g;
    (P_s - P_dc)/(C_dc*U_dc);
    U_dc - U_dcref;
    (v_g_d - u_gd - R_g*i_gd + w_g*L_g*i_gq)/L_g;
    (v_g_q - u_gq - R_g*i_gq - w_g*L_g*i_gd)/L_g;
    i_gdref - i_gd;
    i_gqref - i_gq;
];

end

function plot_diagnostic_results(results)
%% 绘制诊断结果

figure('Position', [100, 100, 1200, 800]);

% 子图1: 特征值实部
subplot(2,2,1);
valid_idx = results(:,4) == 1;
plot(results(valid_idx,1)*180/pi, results(valid_idx,2), 'bo-', 'LineWidth', 2);
xlabel('xi (度)'); ylabel('最大特征值实部');
title('系统稳定性分析'); grid on;
yline(0, 'r--', '稳定性边界', 'LineWidth', 1.5);

% 子图2: 条件数
subplot(2,2,2);
semilogy(results(valid_idx,1)*180/pi, results(valid_idx,3), 'ro-', 'LineWidth', 2);
xlabel('xi (度)'); ylabel('雅可比矩阵条件数 (log)');
title('系统刚性分析'); grid on;
yline(1e10, 'r--', '刚性阈值', 'LineWidth', 1.5);

% 子图3: 稳态求解成功率
subplot(2,2,3);
bar(results(:,1)*180/pi, results(:,4), 'FaceColor', [0.2 0.6 0.8]);
xlabel('xi (度)'); ylabel('稳态求解成功 (1=成功, 0=失败)');
title('稳态求解成功率'); grid on;

% 子图4: 综合评估
subplot(2,2,4);
stability_score = zeros(size(results,1), 1);
for i = 1:size(results,1)
    if results(i,4) == 1  % 稳态求解成功
        if results(i,2) < -1e-3 && results(i,3) < 1e10
            stability_score(i) = 3;  % 优秀
        elseif results(i,2) < 0 && results(i,3) < 1e12
            stability_score(i) = 2;  % 良好
        else
            stability_score(i) = 1;  % 一般
        end
    else
        stability_score(i) = 0;  % 失败
    end
end

bar(results(:,1)*180/pi, stability_score, 'FaceColor', [0.8 0.4 0.2]);
xlabel('xi (度)'); ylabel('稳定性评分');
title('综合稳定性评估'); grid on;
ylim([0 3.5]);
set(gca, 'YTick', 0:3, 'YTickLabel', {'失败', '一般', '良好', '优秀'});

sgtitle('直驱风机系统数值稳定性诊断报告', 'FontSize', 16, 'FontWeight', 'bold');

end

function provide_recommendations(results)
%% 提供改进建议

fprintf('\n=== 系统改进建议 ===\n');

% 分析结果
success_rate = mean(results(:,4));
valid_idx = results(:,4) == 1;

if sum(valid_idx) > 0
    avg_eigenval = mean(results(valid_idx, 2));
    avg_condition = mean(results(valid_idx, 3));
    
    fprintf('1. 总体分析:\n');
    fprintf('   - 稳态求解成功率: %.1f%%\n', success_rate*100);
    fprintf('   - 平均最大特征值实部: %.2e\n', avg_eigenval);
    fprintf('   - 平均雅可比条件数: %.2e\n', avg_condition);
    
    fprintf('\n2. 具体建议:\n');
    
    if avg_eigenval > -1e-3
        fprintf('   ⚠️  系统接近不稳定边界，建议:\n');
        fprintf('      - 进一步降低PLL增益 (Kp_PLL < 50, Ki_PLL < 500)\n');
        fprintf('      - 增加阻尼项或使用更保守的控制器参数\n');
    end
    
    if avg_condition > 1e10
        fprintf('   ⚠️  系统刚性严重，建议:\n');
        fprintf('      - 使用刚性求解器 (ode15s, ode23t)\n');
        fprintf('      - 适当增大时间常数 T_d, T_m\n');
        fprintf('      - 降低控制器带宽\n');
    end
    
    % 找出问题xi范围
    problem_xi = results(results(:,4)==0 | (results(:,4)==1 & results(:,2)>-1e-3), 1);
    if ~isempty(problem_xi)
        fprintf('   ⚠️  问题xi范围: [%.1f°, %.1f°]\n', ...
                min(problem_xi)*180/pi, max(problem_xi)*180/pi);
        fprintf('      - 在此范围内建议使用更保守的参数设置\n');
    end
    
else
    fprintf('❌ 所有测试点都失败，建议:\n');
    fprintf('   - 检查模型参数设置\n');
    fprintf('   - 降低所有控制器增益\n');
    fprintf('   - 使用更保守的初值估计\n');
end

fprintf('\n3. 数值求解建议:\n');
fprintf('   - 优先使用 ode15s 刚性求解器\n');
fprintf('   - 设置适当的容差: RelTol=1e-6, AbsTol=1e-9\n');
fprintf('   - 限制最大步长: MaxStep=0.001\n');
fprintf('   - 如果仍有问题，尝试 ode23t 或 ode23tb\n');

fprintf('\n=== 诊断完成 ===\n');

end