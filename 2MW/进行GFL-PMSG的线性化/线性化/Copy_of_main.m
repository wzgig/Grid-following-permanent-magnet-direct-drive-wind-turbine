clear;
clc;

%% ============ 可视化配置 ============

% 总开关：false = 完全不画图  true
ENABLE_VOC_PLOT = true;

% 是否只对"选定工况"画图：
%   false = 只要 ENABLE_VOC_PLOT=true，就对所有工况画图
%   true  = 只对 PLOT_OPS 列表里的工况画图
PLOT_SELECTED_ONLY = false;

% 想要画图的工况列表，每一行是 [P, Q, V, xi]
% 你可以按需要增加 / 删除行
PLOT_OPS = [
    -0.65   -0.3   1.1   0.35;
    -0.64   -0.3   1.1   0.35;
];

% 数值比较时的容差（工况是浮点数，不建议用严格 ==）
OP_TOL = 1e-6;

%% ===================================

% Path_root_Results = "your root";
% 
% 
% if ~exist(Path_root_Results, 'dir')
%     mkdir(Path_root_Results);
%     disp(['Folder ', Path_root_Results, ' created.']);
% else
%     disp(['Folder ', ' already exists.']);
% end

% 把占位符改成你真正想存结果的路径
% 例如：
% Path_root_Results = "D:\Projects\VOC_results";
% 或者相对路径：
% Path_root_Results = "results";
Path_root_Results = "your root";

% 建议用 isfolder，更直观；mkdir 可自动创建多级目录
if ~isfolder(Path_root_Results)
    mkdir(Path_root_Results);
    fprintf('Folder %s created.\n', Path_root_Results);
else
    fprintf('Folder %s already exists.\n', Path_root_Results);
end


% % 定义工况参数范围和步长
% P_range = linspace(-0.3, -1, 5); 
% Q_range = linspace(-0.3, 0.3, 5);
% V_range = linspace(0.9, 1.1, 5);
% xi_range = linspace(-0.35, 0.35, 20);

P_range = [-0.65, -0.64];
Q_range = -0.3;
V_range = 1.1;
xi_range = 0.35;

% P=-1;
% Q=0;
% V=0.9998;
% xi=0.6286;

for iP = 1:length(P_range)
    P = P_range(iP);
    for iQ = 1:length(Q_range)
        Q = Q_range(iQ);
        for iV = 1:length(V_range)
            V = V_range(iV);
            for ixi = 1:length(xi_range)
                xi = xi_range(ixi);


syms PLL_int thetapll U_dc Udc_int i_gd i_gq Id_grid_int Iq_grid_int P_s i_gqref                
%直流侧潮流
f   = 60;
w_g = 2*pi*f;

Sbase = 2e6;
Pnom = Sbase;
Vbase = 690;
Ibase = Sbase/(sqrt(3)*Vbase);
wbase = w_g;
Zbase = Vbase^2/Sbase;
Lbase = Zbase/w_g;
U_dcref = 1150;

S_D0 = P/V;                      % 有功功率对应的电流分量基础（p.u.）
S_Q0 = Q/V;                      % 无功功率对应的电流分量基础（p.u.）
S_ab0 = S_D0 + 1j*S_Q0;          % αβ坐标系下的复功率表示（p.u.）
i_abs = abs(S_ab0);              % 电流幅值（p.u.）
ui_argdiff = angle(S_ab0);       % 电压与电流之间的相角差（rad）
i_arg = xi - ui_argdiff;         % 电流在αβ坐标系中的相位（rad，xi为电网电压α轴相位）
i_ab = i_abs * exp(1i * i_arg);  % αβ坐标下的复电流（i_ab = i_α + j*i_β）
iq = imag(i_ab);                % β轴电流分量（p.u.）
id = real(i_ab);                % α轴电流分量（p.u.）

v_ab = V * exp(1i * xi);         % αβ坐标下的复电压（v_ab = v_α + j*v_β，xi为电压α轴相位）
v_a = real(v_ab);                % α轴电压分量（p.u.，与电网A相电压对齐）
v_b = imag(v_ab);                % β轴电压分量（p.u.，滞后α轴90°电角度）

% PLL
theta_g = thetapll;         % 网侧PLL估计的同步角
% 构造变换矩阵，将d-q轴量转换为定轴坐标系（用于测量或控制）
Tg_DQdq = [cos(theta_g) -sin(theta_g); sin(theta_g) cos(theta_g)];
vg_dq_real = Vbase * Tg_DQdq * [v_a; v_b];  % GSC中的实际电压
% 变换后各轴上的电压值提取
v_g_d = vg_dq_real(1);          % 网侧d轴电压
v_g_q = -vg_dq_real(2);         % 网侧q轴电压

Kp_Udc = 1;
Ki_Udc = 10;
Kp_Id_grid = 1;
Ki_Id_grid = 15;
Kp_Iq_grid = 1;
Ki_Iq_grid = 15;

Kp_PLL = 350;
Ki_PLL = 3200;
%滤波部分 
R_g = 0.05;
L_g = 7e-3;
C_dc = 1e-3;

i_gdref = Kp_Udc*(U_dc - U_dcref) + Ki_Udc*Udc_int;
% i_gqref = 0;
u_gd = - Kp_Id_grid*(i_gdref - i_gd) - Ki_Id_grid*Id_grid_int - R_g*i_gd + v_g_d+i_gq*w_g*L_g;
u_gq = - Kp_Iq_grid*(i_gqref - i_gq) - Ki_Iq_grid*Iq_grid_int - R_g*i_gq + v_g_q-i_gd*w_g*L_g;

% P_s = Pnom;
P_dc = sqrt(3)*(u_gd* i_gd + u_gq * i_gq);%(2-17)
P_g = sqrt(3)*(v_g_d * i_gd + v_g_q * i_gq);%(2-17)
Q_g = sqrt(3)*(v_g_q * i_gd - v_g_d * i_gq);%(2-18)

Tg_dqDQ = [cos(-theta_g) -sin(-theta_g); sin(-theta_g) cos(-theta_g)];
ig_dq = (1/Ibase)*Tg_dqDQ*[i_gd;i_gq];

i_gD = ig_dq(1);
i_gQ = ig_dq(2);


eqn1 = v_g_q == 0;
eqn2 = Kp_PLL*v_g_q+Ki_PLL*PLL_int+w_g == 0;
% eqn2 = Kp_PLL*v_g_q+Ki_PLL*PLL_int == 0;
eqn3 = (P_s - P_dc)/(C_dc*U_dc) == 0;
eqn4 = U_dc - U_dcref == 0;
eqn5 = (v_g_d - u_gd - R_g*i_gd + w_g*L_g*i_gq)/L_g == 0;
eqn6 = (v_g_q - u_gq - R_g*i_gq - w_g*L_g*i_gd)/L_g == 0;
eqn7 = i_gdref - i_gd == 0;
eqn8 = i_gqref - i_gq == 0;
eqn9 = i_gD - id == 0;
eqn10 = i_gQ - iq == 0;

[PLL_int thetapll U_dc Udc_int i_gd i_gq Id_grid_int Iq_grid_int P_s i_gqref]...
                      = vpasolve(eval([eqn1, eqn2, eqn3, eqn4, eqn5, eqn6, eqn7,...
                      eqn8, eqn9, eqn10]),...
                      [PLL_int thetapll U_dc Udc_int i_gd i_gq Id_grid_int Iq_grid_int P_s i_gqref],...
                      [-inf,inf;-inf,inf;-inf,inf;-inf,inf;...
                      -inf,inf;-inf,inf;-inf,inf;-inf,inf;-inf,inf;-inf,inf]);

inistate = [PLL_int thetapll U_dc Udc_int i_gd i_gq Id_grid_int Iq_grid_int]';

P_s = getRealIfSmallImag(P_s);
i_gqref = getRealIfSmallImag(i_gqref);
% disp(inistate);disp(P_s);disp(i_gqref);

state_size = size(inistate,1);

inistate2 = [PLL_int thetapll U_dc Udc_int i_gd i_gq Id_grid_int Iq_grid_int]';




x_ss = inistate2;  % 稳态状态
% u_ss = [v_a, v_b, P_ref, Q_ref];           % 稳态输入
u_ss = [v_a, v_b, P_s, i_gqref];           % 稳态输入
% y_ss = yy;

epsilon = 1e-7;
n = length(x_ss);  % 状态维度
m = length(u_ss);  % 输入维度
p = 2;  % 输出维度

% 初始化矩阵
A = zeros(n, n);
B = zeros(n, m);
C = zeros(p, n);
D = zeros(p, m);

% 计算 A 和 C
for i = 1:n
    x_perturbed = x_ss;
    x_perturbed(i) = x_perturbed(i) + epsilon;
    [dxdt_perturbed, y_perturbed] = VOC(0, x_perturbed, u_ss);
    [dxdt_nominal, y_nominal] = VOC(0, x_ss, u_ss);
    
    A(:, i) = (dxdt_perturbed - dxdt_nominal) / epsilon;
    C(:, i) = (y_perturbed - y_nominal) / epsilon;
end

% 计算 B 和 D
for i = 1:m
    u_perturbed = u_ss;
    u_perturbed(i) = u_perturbed(i) + epsilon;
    [dxdt_perturbed, y_perturbed] = VOC(0, x_ss, u_perturbed);
    [dxdt_nominal, y_nominal] = VOC(0, x_ss, u_ss);
    
    B(:, i) = (dxdt_perturbed - dxdt_nominal) / epsilon;
    D(:, i) = (y_perturbed - y_nominal) / epsilon;
end

sys = ss(A,B,C,D);

eigenvalues = eig(A);

%% ========= 可视化（可选 + 可选工况）=========
do_plot = false;

if ENABLE_VOC_PLOT
    if ~PLOT_SELECTED_ONLY
        % 模式：所有工况都画
        do_plot = true;
    else
        % 模式：只画 PLOT_OPS 列表中配置的工况
        current_op = [P, Q, V, xi];   % 当前工况
        diff_mat   = abs(PLOT_OPS - current_op);  % 和列表里每一行比较
        match_rows = all(diff_mat < OP_TOL, 2);   % 哪些行是"匹配的"
        if any(match_rows)
            do_plot = true;
        end
    end
end

if do_plot
    visualize_voc_linearization(sys, eigenvalues, P, Q, V, xi);
end
%% ============================================

% filename = sprintf('P%d_V%d_Q%d_X%d.mat', iP, iV, iQ, ixi);
% path = sprintf('%s%s', Path_root_Results, filename);
% save(path, 'sys');

% 后续保存文件时，用 fullfile 来拼路径更安全
filename = sprintf('P%d_V%d_Q%d_X%d.mat', iP, iV, iQ, ixi);
save(fullfile(Path_root_Results, filename), 'sys');

if any(eigenvalues > 0.01)
disp([filename, "is not stable."])
% else
% disp([num2str(P),num2str(V),num2str(Q),num2str(xi),"is stable."])    
end
            end
        end
    end
end