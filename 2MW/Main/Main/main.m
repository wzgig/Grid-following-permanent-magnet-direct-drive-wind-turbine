clear;
clc;


Path_root_Results = "your root";


if ~exist(Path_root_Results, 'dir')
    mkdir(Path_root_Results);
    disp(['Folder ', Path_root_Results, ' created.']);
else
    disp(['Folder ', ' already exists.']);
end


% % 定义参数范围和步长
P_range = linspace(0.3, 1.2, 5);    % P从-0.3到-1.2，步长0.1
Q_range = linspace(0.1, 0.5, 5);      % Q从-0.5到+0.5，步长0.1
V_range = linspace(0.9, 1.1, 5);      % V从0.9到1.1，步长0.05
xi_range = linspace(-pi/2, pi/2, 20);     % xi从-0.6到0.6，步长0.2

% P_range = 0.3;
% Q_range = 0.5;
% V_range = 1.1;
% xi_range = -pi/2; 

% P=-0.8;
% Q=-0.3;
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

%直流侧潮流
Idc=1;
C_bus=4.8e-3;

f   = 60;
w_g = 2*pi*f;

Sbase = 2e6;
Vbase = 690;
Ibase = Sbase/(sqrt(3)*Vbase);
wbase = w_g;
Zbase = Vbase^2/Sbase;
Lbase = Zbase/w_g;

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

% controller
Kpu=6;
Kiu=153.33;
Kpi=4.2;
Kii=0.323;
%w_f=2*pi*5;%一阶滤波器参数
m_p=(2*pi*0.5)/Sbase;
eta=m_p*Vbase^2;
alpha=1e-1*Vbase^2;

%滤波部分 
Lg=0.003;
Rg=0.05;

Un=Vbase;
E0=Vbase;

syms i_gd i_gq v_ref theta Phi_vd Phi_vq Phi_cd Phi_cq P_ref Q_ref 

Tr_dqDQ = [cos(theta) -sin(theta); sin(theta) cos(theta)]; % local to global
is_DQ_pu = 1/Ibase*Tr_dqDQ*[i_gd;i_gq];
i_gD = is_DQ_pu(1);
i_gQ = is_DQ_pu(2);

Tr_DQdq = [cos(theta) sin(theta); -sin(theta) cos(theta)]; %全局到局部
V2_dq_real = Vbase*Tr_DQdq*[v_d;v_q];
u_gd = V2_dq_real(1);
u_gq = V2_dq_real(2);

Vc_dref = v_ref;
Vc_qref = 0;

% Voltage Loop
I1_dref = Phi_vd + Kpu*(Vc_dref - u_gd);
I1_qref = Phi_vq + Kpu*(Vc_qref - u_gq);
% Current Loop
vrefd = Phi_cd - w_g*Lg*i_gq + Kpi*(I1_dref - i_gd) + Rg*i_gd;
vrefq = Phi_cq + w_g*Lg*i_gd + Kpi*(I1_qref - i_gq) + Rg*i_gq;

Pe = sqrt(3)*(u_gd*i_gd + u_gq*i_gq);
Qe = sqrt(3)*(u_gq*i_gd - u_gd*i_gq);
P_e = sqrt(3)*(vrefd*i_gd + vrefq*i_gq);
Vc = sqrt(u_gd^2+u_gq^2);

eqn1  = (vrefd+w_g*Lg*i_gq-u_gd-Rg*i_gd)/Lg==0;
eqn2  = (vrefq-w_g*Lg*i_gd-u_gq-Rg*i_gq)/Lg==0;
eqn3  = v_ref*(eta*(Q_ref*Sbase/Vc^2-Qe/v_ref^2)+eta*alpha*(Vc^2-v_ref^2)/Vc^2)==0;
eqn4  = eta*(P_ref*Sbase/Vc^2-Pe/v_ref^2)==0;
eqn5  = Kiu*(Vc_dref-u_gd)==0;
eqn6  = Kiu*(Vc_qref-u_gq)==0;
eqn7  = Kii*(I1_dref-i_gd)==0;
eqn8  = Kii*(I1_qref-i_gq)==0;
eqn9  = i_gD == i_D;
eqn10  = i_gQ == i_Q;


[i_gd i_gq v_ref theta Phi_vd Phi_vq Phi_cd Phi_cq P_ref Q_ref]...
                      = vpasolve(eval([eqn1, eqn2, eqn3, eqn4, eqn5, eqn6, eqn7,...
                      eqn8, eqn9, eqn10]),...
                      [i_gd i_gq v_ref theta Phi_vd Phi_vq Phi_cd Phi_cq P_ref Q_ref],...
                      [-inf,inf;-inf,inf;0,inf;-2*pi,2*pi;...
                      -inf,inf;-inf,inf;-inf,inf;-inf,inf;-inf,inf;-inf,inf]);

inistate = [i_gd i_gq v_ref theta Phi_vd Phi_vq Phi_cd Phi_cq]';

P_ref = getRealIfSmallImag(P_ref);
Q_ref = getRealIfSmallImag(Q_ref);

state_size = size(inistate,1);

inistate2 = [i_gd i_gq v_ref theta Phi_vd Phi_vq Phi_cd Phi_cq]';




x_ss = inistate2;  % 稳态状态
u_ss = [v_d, v_q, P_ref, Q_ref];           % 稳态输入
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

filename = sprintf('P%d_V%d_Q%d_X%d.mat', iP, iV, iQ, ixi);
path = sprintf('%s%s', Path_root_Results, filename);
save(path, 'sys');

if any(eigenvalues > 0.01)
disp([filename, "is not stable."])
% else
% disp([num2str(P),num2str(V),num2str(Q),num2str(xi),"is stable."])    
end
            end
        end
    end
end